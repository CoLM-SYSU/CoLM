#include <define.h>

#if (defined MaterialTransport) && (defined BGC)
MODULE MOD_MT_Methane_GIEMS

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Loads GIEMS-MC v1.1 monthly inundated+saturated wetland fraction
!   (Prigent+2020) and exposes a per-patch true monthly time series to
!   the methane physics. Activated through `DEF_wetland_finundation_scheme = 5`.
!
! INPUT FILE:
!   `DEF_file_GIEMS` (namelist); required when scheme 5 is selected.
!   Expected variable: `inund_sat_wetland_frac(time, latitude, longitude)`
!   - time: months since 1992-01-01 (monthly)
!            348 months = 1992-01 to 2020-12
!   - latitude/longitude: degrees, regular grid (any resolution)
!   - Sentinel values: -999 ocean, -998 snow, -997 urban; treated as 0.
!
! INIT ALGORITHM (called once during methane init):
!   1. Master opens file, reads dims + lat/lon coords (small).
!   2. Broadcast dims + coords; allocate per-rank ts (small).
!   3. Each rank computes nearest-pixel (best_ix, best_iy) per patch and
!      deduplicates pixels shared by its local patches.
!   4. Streaming loop in chunks of at most 12 months:
!      a. Master reads one (nlon, nlat) slab (~4 MB) per month.
!      b. Master packs only unique values requested by each rank.
!      c. One MPI_Scatterv per chunk sends each rank its unique values.
!   5. Build 12-month climatology fallback from full ts.
!
! RUNTIME:
!   FUNCTION giems_finundated(year, day_of_year) returns:
!   - The actual (year, month) value when 1992 <= year <= 2020
!     (captures El Nino/La Nina anomalies, long-term trends)
!   - Else: 12-month climatology fallback
!
! Notes:
!   - Master holds one flattened pixel index per rank-local unique request and a
!     bounded packed-value chunk (target <=64 MiB when one month fits).
!   - The source series uses real(r4); climatology sums remain real(r8).
!   - GIEMS-MC v1.1 is a 0.25 x 0.25 degree product.  This interface
!     receives patch centres only, so it samples the nearest source cell.
!     A conservative area average cannot be reconstructed from a centre:
!     it needs each patch's underlying mesh-pixel bounds and area weights.
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only: PI
   USE, INTRINSIC :: IEEE_ARITHMETIC, only: ieee_is_finite, ieee_is_nan

   IMPLICIT NONE
   SAVE
   PRIVATE

   ! Public state
   ! True monthly time series, 1992-01-01 .. 2020-12-01 (348 months)
   real(r4), allocatable, public :: giems_ts_wetland_frac(:,:)   ! (ntime, numpatch) [0-1]
   ! 12-month climatology fallback for years outside 1992-2020
   real(r8), allocatable, public :: giems_clim_wetland_frac(:,:) ! (12, numpatch) [0-1]
   integer,  public :: giems_year_start = 1992    ! GIEMS-MC v1.1 first year
   integer,  public :: giems_year_end   = 2020    ! GIEMS-MC v1.1 last year
   integer,  public :: giems_ntime      = 0       ! actual N months read
   logical,  public :: giems_active = .false.

   ! Public API
   PUBLIC :: allocate_methane_giems
   PUBLIC :: deallocate_methane_giems
   PUBLIC :: read_methane_giems
   PUBLIC :: giems_finundated

CONTAINS

   SUBROUTINE allocate_methane_giems(numpatch)
      integer, intent(in) :: numpatch
      IF (allocated(giems_clim_wetland_frac)) THEN
         IF (size(giems_clim_wetland_frac,1) == 12 .and. size(giems_clim_wetland_frac,2) == numpatch) RETURN
         deallocate(giems_clim_wetland_frac)
      ENDIF
      allocate(giems_clim_wetland_frac(12, numpatch))
      giems_clim_wetland_frac(:,:) = 0._r8
   END SUBROUTINE allocate_methane_giems

   SUBROUTINE deallocate_methane_giems()
      IF (allocated(giems_ts_wetland_frac))   deallocate(giems_ts_wetland_frac)
      IF (allocated(giems_clim_wetland_frac)) deallocate(giems_clim_wetland_frac)
      giems_active = .false.
   END SUBROUTINE deallocate_methane_giems

   SUBROUTINE read_methane_giems(file_giems, patchlatr_in, patchlonr_in, numpatch)
      ! Stream-read full GIEMS-MC monthly time series + build fallback
      ! climatology.  Memory-efficient: master holds only one slab at
      ! a time; each rank stores its own per-patch time series.
#ifdef USEMPI
      USE MPI
#endif
      USE netcdf

      character(len=*), intent(in) :: file_giems
      real(r8), intent(in) :: patchlatr_in(:), patchlonr_in(:)   ! radians
      integer, intent(in)  :: numpatch

      character(len=*), parameter :: varname = 'inund_sat_wetland_frac'
      character(len=*), parameter :: latname = 'latitude'
      character(len=*), parameter :: lonname = 'longitude'
      integer, parameter :: giems_expected_months = 348
      integer, parameter :: giems_chunk_months = 12
      integer, parameter :: giems_max_packed_values = 16 * 1024 * 1024

      integer :: ncid, vid, ierr, ierr_bc
      integer :: ntime, nlat, nlon
      integer :: t, mo, ilat, ilon, ipatch, irequest, month_in_chunk
      integer :: ndims, vdims(3), dlen
      integer :: giems_metadata_error, giems_block_error
      integer :: giems_count_error, giems_count_error_glb
      integer :: giems_mapping_error, giems_mapping_error_glb
      integer :: giems_value_error, giems_value_error_glb
      integer :: total_requests, chunk_n, chunk_max, failed_t
      integer :: pixel, n_unique, comm_size
      integer :: metadata(4)
      real(r4), allocatable :: slab(:,:)
      real(r4), allocatable :: unique_values(:,:), requested_values(:,:)
      real(r8), allocatable :: lat_g(:), lon_g(:)
      integer,  allocatable :: best_ix(:), best_iy(:)
      integer,  allocatable :: patch_to_unique(:), pixel_index_unique(:)
      integer,  allocatable :: pixel_to_unique(:), all_pixel_index(:)
      integer,  allocatable :: request_counts(:), request_displs(:)
      integer,  allocatable :: chunk_counts(:), chunk_displs(:)
      integer,  allocatable :: ccnt(:,:)            ! valid-count per (month, patch) for climatology
      real(r8) :: lat_deg, lon_deg, dlatd, dlond, dmin, dlat_min, dlon_min
      logical  :: fexists, found_match
      character(len=256) :: dname
      real(r4) :: v

      ! IMPORTANT: do NOT RETURN here on numpatch<=0.  Master rank often
      ! has no patches but must still participate in the MPI_Bcast
      ! collectives below; early return on master causes worker MPI
      ! deadlock.
      giems_active = .false.
      giems_ntime = 0
      giems_metadata_error = 0
      ntime = 0
      nlat = 0
      nlon = 0
      metadata = 0

      ! ---- Step 1: master opens file, reads coords ----
      IF (p_is_master) THEN
         INQUIRE(file=trim(file_giems), exist=fexists)
         IF (.not. fexists) THEN
            write(*,'(A,A)') ' ERROR: GIEMS file not found: ', trim(file_giems)
            giems_metadata_error = 1
            ntime = 0; nlat = 0; nlon = 0
         ELSE
            ierr = nf90_open(trim(file_giems), NF90_NOWRITE, ncid)
            IF (ierr /= NF90_NOERR) THEN
               write(*,'(A,A)') ' ERROR: GIEMS open failed: ', trim(nf90_strerror(ierr))
               giems_metadata_error = 1
               ntime = 0; nlat = 0; nlon = 0
            ELSE
               CALL get_dim(ncid, 'time',  ntime)
               CALL get_dim(ncid, latname, nlat)
               CALL get_dim(ncid, lonname, nlon)
               IF (ntime <= 0 .or. nlat <= 0 .or. nlon <= 0) THEN
                  write(*,'(A)') ' ERROR: GIEMS is missing required time/latitude/longitude dimensions.'
                  giems_metadata_error = 1
                  ierr = nf90_close(ncid)
               ELSEIF (ntime /= giems_expected_months) THEN
                  write(*,'(A,I0,A,I0,A)') ' ERROR: GIEMS time dimension has ', ntime, &
                     ' months; expected exactly ', giems_expected_months, ' (1992-01 through 2020-12).'
                  giems_metadata_error = 1
                  ierr = nf90_close(ncid)
               ELSEIF (int(nlat, i8) * int(nlon, i8) > int(huge(0), i8)) THEN
                  write(*,'(A)') ' ERROR: GIEMS grid is too large for the default-integer flattened pixel index.'
                  giems_metadata_error = 1
                  ierr = nf90_close(ncid)
               ELSE
                  write(*,'(A,I0,A,I0,A,I0)') ' GIEMS dims: time=', ntime, ' lat=', nlat, ' lon=', nlon

                  CALL validate_giems_time_axis(ncid, ntime, giems_metadata_error)

                  allocate(lat_g(nlat), lon_g(nlon))
                  ierr = nf90_inq_varid(ncid, latname, vid)
                  IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, lat_g)
                  IF (ierr /= NF90_NOERR) THEN
                     write(*,'(A,A)') ' ERROR: GIEMS latitude read failed: ', trim(nf90_strerror(ierr))
                     giems_metadata_error = 1
                  ELSE
                     ierr = nf90_inq_varid(ncid, lonname, vid)
                     IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, lon_g)
                     IF (ierr /= NF90_NOERR) THEN
                        write(*,'(A,A)') ' ERROR: GIEMS longitude read failed: ', trim(nf90_strerror(ierr))
                        giems_metadata_error = 1
                     ELSEIF (any(.not. ieee_is_finite(lat_g)) .or. &
                             any(.not. ieee_is_finite(lon_g)) .or. &
                             any(abs(lat_g) > 90._r8)) THEN
                        write(*,'(A)') ' ERROR: GIEMS latitude/longitude coordinates are invalid.'
                        giems_metadata_error = 1
                     ELSE
                        ierr = nf90_inq_varid(ncid, varname, vid)
                        IF (ierr /= NF90_NOERR) THEN
                           write(*,'(A,A)') ' ERROR: GIEMS variable not found: ', trim(nf90_strerror(ierr))
                           giems_metadata_error = 1
                        ELSE
                           ierr = nf90_inquire_variable(ncid, vid, ndims=ndims)
                           IF (ierr /= NF90_NOERR) THEN
                              write(*,'(A,A)') ' ERROR: GIEMS variable metadata read failed: ', &
                                 trim(nf90_strerror(ierr))
                              giems_metadata_error = 1
                           ELSEIF (ndims /= 3) THEN
                              write(*,'(A,I0,A)') ' ERROR: GIEMS variable has ', ndims, &
                                 ' dimensions; expected longitude, latitude, time.'
                              giems_metadata_error = 1
                           ELSE
                              ierr = nf90_inquire_variable(ncid, vid, dimids=vdims)
                              IF (ierr /= NF90_NOERR) THEN
                                 write(*,'(A,A)') ' ERROR: GIEMS dimension-id read failed: ', &
                                    trim(nf90_strerror(ierr))
                                 giems_metadata_error = 1
                              ENDIF

                              IF (ierr == NF90_NOERR) THEN
                                 dname = ''
                                 dlen = -1
                                 ierr = nf90_inquire_dimension(ncid, vdims(1), name=dname, len=dlen)
                                 IF (ierr /= NF90_NOERR .or. &
                                     trim(dname) /= 'longitude' .or. dlen /= nlon) THEN
                                    write(*,'(A,A,A,I0,A,I0)') ' ERROR: GIEMS dimension 1 is "', &
                                       trim(dname), '" len=', dlen, '; expected longitude len=', nlon
                                    giems_metadata_error = 1
                                 ENDIF

                                 dname = ''
                                 dlen = -1
                                 ierr = nf90_inquire_dimension(ncid, vdims(2), name=dname, len=dlen)
                                 IF (ierr /= NF90_NOERR .or. &
                                     trim(dname) /= 'latitude' .or. dlen /= nlat) THEN
                                    write(*,'(A,A,A,I0,A,I0)') ' ERROR: GIEMS dimension 2 is "', &
                                       trim(dname), '" len=', dlen, '; expected latitude len=', nlat
                                    giems_metadata_error = 1
                                 ENDIF

                                 dname = ''
                                 dlen = -1
                                 ierr = nf90_inquire_dimension(ncid, vdims(3), name=dname, len=dlen)
                                 IF (ierr /= NF90_NOERR .or. &
                                     trim(dname) /= 'time' .or. dlen /= ntime) THEN
                                    write(*,'(A,A,A,I0,A,I0)') ' ERROR: GIEMS dimension 3 is "', &
                                       trim(dname), '" len=', dlen, '; expected time len=', ntime
                                    giems_metadata_error = 1
                                 ENDIF
                              ENDIF
                           ENDIF
                           giems_active = giems_metadata_error == 0
                        ENDIF
                     ENDIF
                  ENDIF
                  IF (.not. giems_active) THEN
                     IF (allocated(lat_g)) deallocate(lat_g)
                     IF (allocated(lon_g)) deallocate(lon_g)
                     ierr = nf90_close(ncid)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      IF (p_is_master) metadata = [giems_metadata_error, ntime, nlat, nlon]
#ifdef USEMPI
      CALL MPI_Bcast(metadata, 4, MPI_INTEGER, p_address_master, p_comm_glb, ierr)
      CALL check_giems_mpi(ierr, 'metadata broadcast')
      giems_metadata_error = metadata(1)
      ntime = metadata(2)
      nlat = metadata(3)
      nlon = metadata(4)
      giems_active = giems_metadata_error == 0
      IF (giems_active) THEN
         IF (.not. p_is_master) THEN
            allocate(lat_g(nlat), lon_g(nlon))
         ENDIF
         CALL MPI_Bcast(lat_g, nlat, MPI_DOUBLE_PRECISION, p_address_master, p_comm_glb, ierr)
         CALL check_giems_mpi(ierr, 'latitude broadcast')
         CALL MPI_Bcast(lon_g, nlon, MPI_DOUBLE_PRECISION, p_address_master, p_comm_glb, ierr)
         CALL check_giems_mpi(ierr, 'longitude broadcast')
      ENDIF
#endif

      IF (giems_metadata_error /= 0) THEN
         IF (p_is_master) THEN
            CALL CoLM_stop ('ERROR: invalid GIEMS metadata; expected longitude, latitude, time.')
         ELSE
            CALL CoLM_stop ()
         ENDIF
      ENDIF

      IF (.not. giems_active) RETURN
      giems_ntime = ntime
      giems_year_end = giems_year_start + (ntime - 1) / 12

      ! patchlatr_in/patchlonr_in are worker-local coordinate vectors.
      ! Non-worker ranks call this collective routine with numpatch=0; worker
      ! ranks must provide arrays long enough for their local patch count.
      IF (numpatch < 0) THEN
         write(*,'(A,I0)') ' ERROR: read_methane_giems received negative numpatch=', numpatch
         CALL CoLM_stop ()
      ENDIF
      IF (numpatch > 0) THEN
         IF (size(patchlatr_in) < numpatch .or. size(patchlonr_in) < numpatch) THEN
            write(*,'(A,I0,A,I0,A,I0)') &
               ' ERROR: read_methane_giems patch-coordinate size mismatch: numpatch=', numpatch, &
               ' size(patchlatr_in)=', size(patchlatr_in), &
               ' size(patchlonr_in)=', size(patchlonr_in)
            CALL CoLM_stop ()
         ENDIF
      ENDIF

      ! ---- Step 2: each rank computes nearest-pixel per patch (cached) ----
      ! Split the lookup into independent nearest-lat and nearest-lon scans:
      ! O(npatch*(nlat+nlon)) rather than the old nested
      ! O(npatch*nlat*nlon).  This removes a global-run init hot spot while
      ! preserving the same nearest-gridcell mapping on a regular lon/lat
      ! GIEMS grid.
      allocate(best_ix(numpatch))
      allocate(best_iy(numpatch))
      best_ix(:) = -1
      best_iy(:) = -1
      giems_mapping_error = 0
      DO ipatch = 1, numpatch
         IF (.not. ieee_is_finite(patchlatr_in(ipatch)) .or. &
             .not. ieee_is_finite(patchlonr_in(ipatch)) .or. &
             abs(patchlatr_in(ipatch)) > 0.5_r8 * PI .or. &
             abs(patchlonr_in(ipatch)) > 2._r8 * PI) THEN
            giems_mapping_error = 1
            CYCLE
         ENDIF
         lat_deg = patchlatr_in(ipatch) * 180._r8 / PI
         lon_deg = patchlonr_in(ipatch) * 180._r8 / PI
         found_match = .false.

         dmin = huge(0._r8)
         DO ilat = 1, nlat
            dlatd = abs(lat_deg - lat_g(ilat))
            IF (dlatd < dmin) THEN
               dmin = dlatd
               best_iy(ipatch) = ilat
            ENDIF
         ENDDO
         dlat_min = dmin
         IF (dlat_min > 5._r8) THEN
            best_iy(ipatch) = -1
            giems_mapping_error = 1
            CYCLE
         ENDIF

         dmin = huge(0._r8)
         DO ilon = 1, nlon
            dlond = modulo(abs(lon_deg - lon_g(ilon)), 360._r8)
            dlond = min(dlond, 360._r8 - dlond)
            IF (dlond < dmin) THEN
               dmin = dlond
               best_ix(ipatch) = ilon
            ENDIF
         ENDDO
         dlon_min = dmin
         IF (sqrt(dlat_min**2 + dlon_min**2) <= 5._r8) found_match = .true.
         IF (.not. found_match) THEN
            best_ix(ipatch) = -1
            best_iy(ipatch) = -1
            giems_mapping_error = 1
         ENDIF
      ENDDO

#ifdef USEMPI
      CALL MPI_Allreduce(giems_mapping_error, giems_mapping_error_glb, 1, MPI_INTEGER, &
         MPI_MAX, p_comm_glb, ierr)
      CALL check_giems_mpi(ierr, 'patch mapping validation')
#else
      giems_mapping_error_glb = giems_mapping_error
#endif
      IF (giems_mapping_error_glb /= 0) THEN
         IF (p_is_master) write(*,'(A)') &
            ' ERROR: at least one CoLM patch has invalid coordinates or no matching GIEMS cell.'
         CALL CoLM_stop ()
      ENDIF

      ! ---- Step 3: allocate per-rank time series + climatology accumulator ----
      IF (allocated(giems_ts_wetland_frac)) deallocate(giems_ts_wetland_frac)
      allocate(giems_ts_wetland_frac(ntime, numpatch))
      giems_ts_wetland_frac(:,:) = 0._r4

      IF (allocated(giems_clim_wetland_frac)) THEN
         IF (size(giems_clim_wetland_frac,1) /= 12 .or. size(giems_clim_wetland_frac,2) /= numpatch) THEN
            deallocate(giems_clim_wetland_frac)
         ENDIF
      ENDIF
      IF (.not. allocated(giems_clim_wetland_frac)) THEN
         allocate(giems_clim_wetland_frac(12, numpatch))
      ENDIF
      giems_clim_wetland_frac(:,:) = 0._r8
      allocate(ccnt(12, numpatch))
      ccnt(:,:) = 0

      ! Deduplicate rank-local requests before they cross the network.
      allocate(patch_to_unique(max(1, numpatch)))
      allocate(pixel_index_unique(max(1, numpatch)))
      allocate(pixel_to_unique(max(1, nlon*nlat)))
      patch_to_unique(:) = 0
      pixel_index_unique(:) = 0
      pixel_to_unique(:) = 0
      n_unique = 0
      DO ipatch = 1, numpatch
         IF (best_ix(ipatch) > 0 .and. best_iy(ipatch) > 0) THEN
            pixel = (best_iy(ipatch) - 1) * nlon + best_ix(ipatch)
            IF (pixel_to_unique(pixel) == 0) THEN
               n_unique = n_unique + 1
               pixel_index_unique(n_unique) = pixel
               pixel_to_unique(pixel) = n_unique
            ENDIF
            patch_to_unique(ipatch) = pixel_to_unique(pixel)
         ENDIF
      ENDDO
      deallocate(pixel_to_unique, best_ix, best_iy)
      deallocate(lat_g, lon_g)

      comm_size = 1
#ifdef USEMPI
      CALL MPI_Comm_size(p_comm_glb, comm_size, ierr)
      CALL check_giems_mpi(ierr, 'communicator size query')
      IF (comm_size /= p_np_glb) THEN
         write(*,'(A,I0,A,I0)') ' ERROR: GIEMS communicator/task-count mismatch: MPI size=', &
            comm_size, ' p_np_glb=', p_np_glb
         CALL CoLM_stop ()
      ENDIF
#endif
      allocate(request_counts(comm_size), request_displs(comm_size))
      allocate(chunk_counts(comm_size), chunk_displs(comm_size))
      request_counts(:) = 0
      request_displs(:) = 0
      chunk_counts(:) = 0
      chunk_displs(:) = 0
#ifdef USEMPI
      CALL MPI_Gather(n_unique, 1, MPI_INTEGER, request_counts, 1, MPI_INTEGER, &
         p_address_master, p_comm_glb, ierr)
      CALL check_giems_mpi(ierr, 'unique request count gather')
#else
      request_counts(1) = n_unique
#endif

      total_requests = 0
      giems_count_error = 0
      chunk_max = min(giems_chunk_months, ntime)
      IF (p_is_master) THEN
         DO irequest = 1, comm_size
            IF (request_counts(irequest) < 0) THEN
               giems_count_error = 1
               EXIT
            ELSEIF (total_requests > huge(total_requests) - request_counts(irequest)) THEN
               giems_count_error = 1
               EXIT
            ENDIF
            request_displs(irequest) = total_requests
            total_requests = total_requests + request_counts(irequest)
         ENDDO
         IF (giems_count_error == 0) THEN
            ! Keep the root's real(r4) packing buffer near 64 MiB without
            ! giving up the 12x collective reduction for ordinary layouts.
            chunk_max = min(chunk_max, max(1, giems_max_packed_values / max(1, total_requests)))
            IF (total_requests > huge(total_requests) / chunk_max) giems_count_error = 1
         ENDIF
      ENDIF
#ifdef USEMPI
      CALL MPI_Bcast(chunk_max, 1, MPI_INTEGER, p_address_master, p_comm_glb, ierr)
      CALL check_giems_mpi(ierr, 'chunk size broadcast')
#endif
      IF (n_unique > huge(n_unique) / chunk_max) giems_count_error = 1
#ifdef USEMPI
      CALL MPI_Allreduce(giems_count_error, giems_count_error_glb, 1, MPI_INTEGER, &
         MPI_MAX, p_comm_glb, ierr)
      CALL check_giems_mpi(ierr, 'request count validation')
#else
      giems_count_error_glb = giems_count_error
#endif
      IF (giems_count_error_glb /= 0) THEN
         IF (p_is_master) write(*,'(A,I0,A)') &
            ' ERROR: GIEMS patch distribution exceeds the MPI default-integer count limit for a ', &
            chunk_max, '-month chunk.'
         CALL CoLM_stop ()
      ENDIF

      IF (p_is_master) THEN
         allocate(all_pixel_index(max(1, total_requests)))
      ELSE
         ! MPI ignores root-only receive arguments away from the root, but
         ! allocated one-element buffers keep the legacy `use mpi` interface
         ! valid on implementations that still inspect the actual argument.
         allocate(all_pixel_index(1))
      ENDIF
      all_pixel_index(:) = 0

#ifdef USEMPI
      CALL MPI_Gatherv(pixel_index_unique, n_unique, MPI_INTEGER, all_pixel_index, &
         request_counts, request_displs, MPI_INTEGER, p_address_master, &
         p_comm_glb, ierr)
      CALL check_giems_mpi(ierr, 'unique pixel index gather')
#else
      IF (n_unique > 0) all_pixel_index(1:n_unique) = pixel_index_unique(1:n_unique)
#endif

      ! ---- Step 4: streaming monthly reads + chunked directed distribution ----
      ! Dimension-order note: NetCDF metadata advertises the variable as
      ! inund_sat_wetland_frac(time, latitude, longitude) in C order, but
      ! Fortran reverses dimension order for nf90_get_var so the buffer is
      ! declared slab(nlon, nlat) and start/count are [lon, lat, time].
      ! Do NOT "fix" this by swapping nlat/nlon -- silently transposes the
      ! lookup.  init_methane_giems verifies dim names against the file at
      ! Step 1 (nf90_inq_dim_ids/names), and a future maintainer should keep
      ! that check before touching this read.
      IF (p_is_master) allocate(slab(nlon, nlat))

      DO t = 1, ntime, chunk_max
         chunk_n = min(chunk_max, ntime - t + 1)
         allocate(unique_values(chunk_n, max(1, n_unique)))
         IF (p_is_master) THEN
            allocate(requested_values(chunk_n, max(1, total_requests)))
         ELSE
            allocate(requested_values(chunk_n, 1))
         ENDIF
         unique_values(:,:) = 0._r4
         requested_values(:,:) = 0._r4
         chunk_counts(:) = request_counts(:) * chunk_n
         chunk_displs(:) = request_displs(:) * chunk_n

         giems_block_error = NF90_NOERR
         failed_t = 0
         IF (p_is_master) THEN
            DO month_in_chunk = 1, chunk_n
               ierr = nf90_get_var(ncid, vid, slab, &
                  start=[1, 1, t + month_in_chunk - 1], count=[nlon, nlat, 1])
               IF (ierr /= NF90_NOERR) THEN
                  giems_block_error = ierr
                  failed_t = t + month_in_chunk - 1
                  EXIT
               ENDIF

               DO irequest = 1, total_requests
                  IF (all_pixel_index(irequest) <= 0) THEN
                     requested_values(month_in_chunk, irequest) = -999._r4
                  ELSE
                     ilon = mod(all_pixel_index(irequest) - 1, nlon) + 1
                     ilat = (all_pixel_index(irequest) - 1) / nlon + 1
                     requested_values(month_in_chunk, irequest) = slab(ilon, ilat)
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
#ifdef USEMPI
         CALL MPI_Bcast(giems_block_error, 1, MPI_INTEGER, p_address_master, p_comm_glb, ierr_bc)
         CALL check_giems_mpi(ierr_bc, 'monthly NetCDF status broadcast')
#endif
         IF (giems_block_error /= NF90_NOERR) THEN
            ! Silent zeroing of a failed monthly slab would look like "no
            ! inundation that month" — produces wrong CH4 budget.
            ! Fail loud instead (scheme 5 is observation-driven; partial
            ! failure is not acceptable).
            IF (p_is_master) THEN
               write(*,'(A,I0,A)') ' ERROR: GIEMS slab read failed at time index ', failed_t, &
                  '; aborting.'
               write(*,'(A,A)') '   nf90_strerror: ', trim(nf90_strerror(giems_block_error))
            ENDIF
            CALL CoLM_stop ()
         ENDIF

#ifdef USEMPI
         ! Keep the NetCDF slab ABI explicit and send only unique requests.
         CALL MPI_Scatterv(requested_values, chunk_counts, chunk_displs, MPI_REAL4, &
            unique_values, chunk_n * n_unique, MPI_REAL4, p_address_master, p_comm_glb, ierr)
         CALL check_giems_mpi(ierr, 'monthly unique value scatter')
#else
         IF (n_unique > 0) unique_values(:,1:n_unique) = requested_values(:,1:n_unique)
#endif

         giems_value_error = 0
         DO ipatch = 1, numpatch
            IF (patch_to_unique(ipatch) <= 0) CYCLE
            DO month_in_chunk = 1, chunk_n
               mo = mod(t + month_in_chunk - 2, 12) + 1
               v = unique_values(month_in_chunk, patch_to_unique(ipatch))
               ! NaN is a documented fill value: test it first.  The ordered
               ! comparisons below are signalling compares, so with the
               ! standard -ffpe-trap=invalid build a NaN reaching them aborts
               ! the run instead of being counted as a physical zero.
               IF (ieee_is_nan(v)) THEN
                  CONTINUE
               ELSEIF (v >= 0._r4 .and. v <= 1._r4) THEN
                  giems_ts_wetland_frac(t + month_in_chunk - 1, ipatch) = v
                  giems_clim_wetland_frac(mo, ipatch) = &
                     giems_clim_wetland_frac(mo, ipatch) + real(v, r8)
               ELSEIF (v /= -999._r4 .and. v /= -998._r4 .and. v /= -997._r4) THEN
                  giems_value_error = 1
               ENDIF
               ! Valid observations and documented ocean/snow/urban/NaN
               ! fill values all represent this month.  Count the latter as
               ! physical zero in the climatology instead of excluding them,
               ! which would bias seasonally snow-covered cells high.
               ccnt(mo, ipatch) = ccnt(mo, ipatch) + 1
            ENDDO
         ENDDO
#ifdef USEMPI
         CALL MPI_Allreduce(giems_value_error, giems_value_error_glb, 1, MPI_INTEGER, &
            MPI_MAX, p_comm_glb, ierr)
         CALL check_giems_mpi(ierr, 'wetland value validation')
#else
         giems_value_error_glb = giems_value_error
#endif
         IF (giems_value_error_glb /= 0) THEN
            IF (p_is_master) write(*,'(A)') &
               ' ERROR: GIEMS contains a selected wetland fraction outside [0,1] or documented fill flags.'
            CALL CoLM_stop ()
         ENDIF
         deallocate(unique_values, requested_values)
      ENDDO

      ! Finalize climatology (per-patch monthly mean)
      DO ipatch = 1, numpatch
         DO mo = 1, 12
            IF (ccnt(mo, ipatch) > 0) THEN
               giems_clim_wetland_frac(mo, ipatch) = &
                  giems_clim_wetland_frac(mo, ipatch) / real(ccnt(mo, ipatch), r8)
            ENDIF
         ENDDO
      ENDDO

      IF (p_is_master) ierr = nf90_close(ncid)

      IF (allocated(slab)) deallocate(slab)
      deallocate(ccnt, patch_to_unique, pixel_index_unique)
      deallocate(all_pixel_index, request_counts, request_displs, chunk_counts, chunk_displs)

      IF (p_is_master) write(*,'(A,I0,A,I0,A,I0,A)') &
         ' GIEMS monthly time series loaded (', ntime, ' months, ', &
         giems_year_start, '-', giems_year_end, &
         '); years outside this range use 12-month climatology fallback.'
   END SUBROUTINE read_methane_giems

   SUBROUTINE validate_giems_time_axis(ncid, ntime, metadata_error)
      ! The runtime lookup maps array index 1 to 1992-01.  Verify that the
      ! coordinate actually advertises that origin instead of accepting any
      ! same-shaped but shifted or unordered 348-sample file.
      USE netcdf
      integer, intent(in)    :: ncid, ntime
      integer, intent(inout) :: metadata_error

      integer, parameter :: month_days(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
      integer :: ierr, time_vid, i, c, year, month, days_this_month, expected_day
      real(r8), allocatable :: time_values(:)
      character(len=256) :: time_units, units_lower
      character(len=256) :: time_calendar, calendar_lower
      character(len=256) :: time_resolution, resolution_lower

      ierr = nf90_inq_varid(ncid, 'time', time_vid)
      IF (ierr /= NF90_NOERR) THEN
         write(*,'(A,A)') ' ERROR: GIEMS time coordinate not found: ', trim(nf90_strerror(ierr))
         metadata_error = 1
         RETURN
      ENDIF

      allocate(time_values(ntime))
      ierr = nf90_get_var(ncid, time_vid, time_values)
      IF (ierr /= NF90_NOERR) THEN
         write(*,'(A,A)') ' ERROR: GIEMS time coordinate read failed: ', trim(nf90_strerror(ierr))
         metadata_error = 1
         deallocate(time_values)
         RETURN
      ENDIF

      time_units = ''
      ierr = nf90_get_att(ncid, time_vid, 'units', time_units)
      IF (ierr /= NF90_NOERR) THEN
         write(*,'(A)') ' ERROR: GIEMS time coordinate must define its day-based 1992-01-01 origin.'
         metadata_error = 1
      ELSE
         units_lower = adjustl(time_units)
         DO i = 1, len_trim(units_lower)
            c = iachar(units_lower(i:i))
            IF (c >= iachar('A') .and. c <= iachar('Z')) &
               units_lower(i:i) = achar(c + iachar('a') - iachar('A'))
         ENDDO
         IF (index(trim(units_lower), 'days since 1992-01-01') /= 1) THEN
            write(*,'(A,A)') ' ERROR: GIEMS time units must be days since 1992-01-01: ', trim(time_units)
            metadata_error = 1
         ENDIF
      ENDIF

      time_calendar = ''
      ierr = nf90_get_att(ncid, time_vid, 'calendar', time_calendar)
      calendar_lower = adjustl(time_calendar)
      DO i = 1, len_trim(calendar_lower)
         c = iachar(calendar_lower(i:i))
         IF (c >= iachar('A') .and. c <= iachar('Z')) &
            calendar_lower(i:i) = achar(c + iachar('a') - iachar('A'))
      ENDDO
      IF (ierr /= NF90_NOERR .or. &
          (trim(calendar_lower) /= 'proleptic_gregorian' .and. &
           trim(calendar_lower) /= 'gregorian' .and. trim(calendar_lower) /= 'standard')) THEN
         write(*,'(A,A)') ' ERROR: GIEMS time calendar is missing or not Gregorian: ', trim(time_calendar)
         metadata_error = 1
      ENDIF

      time_resolution = ''
      ierr = nf90_get_att(ncid, time_vid, 'time_resolution', time_resolution)
      resolution_lower = adjustl(time_resolution)
      DO i = 1, len_trim(resolution_lower)
         c = iachar(resolution_lower(i:i))
         IF (c >= iachar('A') .and. c <= iachar('Z')) &
            resolution_lower(i:i) = achar(c + iachar('a') - iachar('A'))
      ENDDO
      IF (ierr /= NF90_NOERR .or. trim(resolution_lower) /= 'monthly_average') THEN
         write(*,'(A,A)') ' ERROR: GIEMS time_resolution must be monthly_average: ', trim(time_resolution)
         metadata_error = 1
      ENDIF

      expected_day = 0
      DO i = 1, ntime
         IF (.not. ieee_is_finite(time_values(i)) .or. &
             abs(time_values(i) - real(expected_day, r8)) > 0.25_r8) THEN
            write(*,'(A,I0,A,F20.6,A,I0)') ' ERROR: GIEMS time index ', i, &
               ' has day=', time_values(i), '; expected month-start day=', expected_day
            metadata_error = 1
            EXIT
         ENDIF
         year = 1992 + (i - 1) / 12
         month = mod(i - 1, 12) + 1
         days_this_month = month_days(month)
         IF (month == 2 .and. &
             (mod(year, 400) == 0 .or. (mod(year, 4) == 0 .and. mod(year, 100) /= 0))) &
            days_this_month = days_this_month + 1
         expected_day = expected_day + days_this_month
      ENDDO
      deallocate(time_values)
   END SUBROUTINE validate_giems_time_axis

   real(r8) FUNCTION giems_finundated(ipatch, year, day_of_year)
      ! Returns per-patch finundated for the (year, day_of_year):
      ! - True monthly value when 1992 <= year <= 2020 (captures El Nino,
      !   long-term trend, inter-annual variability).
      ! - 12-month climatology fallback when year outside that range.
      USE MOD_TimeManager, only: julian2monthday
      integer, intent(in) :: ipatch, year, day_of_year
      integer :: mon, mday, t

      giems_finundated = 0._r8
      IF (.not. giems_active) RETURN
      IF (ipatch < 1) RETURN

      CALL julian2monthday(year, day_of_year, mon, mday)
      IF (mon < 1 .or. mon > 12) RETURN

      ! Preferred: true monthly time series for years in data range
      IF (allocated(giems_ts_wetland_frac) .and. &
          ipatch <= size(giems_ts_wetland_frac, 2)) THEN
         IF (year >= giems_year_start .and. year <= giems_year_end) THEN
            t = (year - giems_year_start) * 12 + mon
            IF (t >= 1 .and. t <= size(giems_ts_wetland_frac, 1)) THEN
               giems_finundated = giems_ts_wetland_frac(t, ipatch)
               RETURN
            ENDIF
         ENDIF
      ENDIF

      ! Fallback: 12-month climatology (year < 1992 or > 2020)
      IF (allocated(giems_clim_wetland_frac) .and. &
          ipatch <= size(giems_clim_wetland_frac, 2)) THEN
         giems_finundated = giems_clim_wetland_frac(mon, ipatch)
      ENDIF
   END FUNCTION giems_finundated

#ifdef USEMPI
   SUBROUTINE check_giems_mpi(ierr, operation)
      USE MPI
      integer, intent(in) :: ierr
      character(len=*), intent(in) :: operation
      integer :: error_length, ierr_string
      character(len=MPI_MAX_ERROR_STRING) :: error_message

      IF (ierr == MPI_SUCCESS) RETURN
      error_message = 'unknown MPI error'
      error_length = len_trim(error_message)
      CALL MPI_Error_string(ierr, error_message, error_length, ierr_string)
      IF (ierr_string == MPI_SUCCESS) THEN
         write(*,'(A,A,A,A)') ' ERROR: GIEMS MPI failure during ', trim(operation), ': ', &
            trim(error_message(1:max(1, error_length)))
      ELSE
         write(*,'(A,A,A,I0,A,I0)') ' ERROR: GIEMS MPI failure during ', trim(operation), &
            ': code=', ierr, '; MPI_Error_string code=', ierr_string
      ENDIF
      CALL CoLM_stop ()
   END SUBROUTINE check_giems_mpi
#endif

   SUBROUTINE get_dim(ncid, name, n)
      USE netcdf
      integer, intent(in)  :: ncid
      character(len=*), intent(in) :: name
      integer, intent(out) :: n
      integer :: did, ierr
      n = 0
      ierr = nf90_inq_dimid(ncid, trim(name), did)
      IF (ierr /= NF90_NOERR) RETURN
      ierr = nf90_inquire_dimension(ncid, did, len=n)
      IF (ierr /= NF90_NOERR) n = 0
   END SUBROUTINE get_dim

END MODULE MOD_MT_Methane_GIEMS
#endif
