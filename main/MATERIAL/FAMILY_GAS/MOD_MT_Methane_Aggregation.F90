#include <define.h>

#if (defined MaterialTransport) && (defined BGC)
MODULE MOD_MT_Methane_Aggregation
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: methane_aggregate_surface_data, methane_clip_surface_data
CONTAINS

   SUBROUTINE methane_clip_surface_data(dir_landdata_in, dir_landdata_out)
      USE MOD_MT_Methane_Preprocessing, only: methane_preprocessing_requirements
      USE MOD_RegionClip, only: srfdata_region_clip
      character(len=*), intent(in) :: dir_landdata_in, dir_landdata_out
      logical :: requires_lake_soilc, requires_spatial_ph

      CALL methane_preprocessing_requirements(requires_lake_soilc, requires_spatial_ph)
      CALL srfdata_region_clip(dir_landdata_in, dir_landdata_out, requires_lake_soilc, requires_spatial_ph)
   END SUBROUTINE methane_clip_surface_data

   SUBROUTINE methane_aggregate_surface_data(gland, dir_rawdata, dir_landdata, lc_year)
      USE MOD_Grid, only: grid_type
      USE MOD_MT_Methane_Preprocessing, only: methane_preprocessing_requirements
      type(grid_type), intent(in) :: gland
      character(len=*), intent(in) :: dir_rawdata, dir_landdata
      integer, intent(in) :: lc_year
      logical :: requires_lake_soilc, requires_spatial_ph

      CALL methane_preprocessing_requirements(requires_lake_soilc, requires_spatial_ph)
      IF (requires_lake_soilc) CALL Aggregation_LakeSoilC(gland, dir_rawdata, dir_landdata, lc_year)
      IF (requires_spatial_ph) CALL Aggregation_MethanePH(dir_rawdata, dir_landdata, lc_year)
   END SUBROUTINE methane_aggregate_surface_data

SUBROUTINE Aggregation_LakeSoilC ( &
      gland, dir_rawdata, dir_model_landdata, lc_year)

!-----------------------------------------------------------------------
! !DESCRIPTION:
!  Aggregate lake-sediment organic carbon for Methane lake production.
!
!  Prefer lake_soilc(soil, lon, lat) [gC/m3] from either raw path below.
!  Only when both files are absent, estimate from aggregated soil organic
!  matter density: 580 [gC/kg OM] * OM_density [kg OM/m3].
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   USE MOD_Vars_Global, only: nl_soil, spval
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_LandPatch
   USE MOD_NetCDFVector
   USE MOD_MT_NcIO, only: ncio_read_vector_complete
   USE MOD_NetCDFBlock
#ifdef RangeCheck
   USE MOD_RangeCheck
#endif
   USE MOD_AggregationRequestData

   IMPLICIT NONE

   integer, intent(in) :: lc_year
   type(grid_type),  intent(in) :: gland
   character(len=*), intent(in) :: dir_rawdata
   character(len=*), intent(in) :: dir_model_landdata

   character(len=256) :: landdir, lndname, cyear, c
   integer :: ipatch, L, j, src_layer, n_source_layers
   integer :: lake_local, lake_global, invalid_local, invalid_global
   logical :: raw_exists
   real(r8) :: valid_area, om_value
   real(r8), parameter :: carbon_per_kg_om = 580._r8    ! gC/kg OM

   type(block_data_real8_3d) :: lake_soilc_grid
   real(r8), allocatable :: lake_soilc_patches(:,:), lake_soilc_one(:,:), area_one(:)
   real(r8), allocatable :: om_density_patches(:)
#ifdef SrfdataDiag
   ! No gridded diagnostic is written here because srfdata_map_and_write is
   ! scalar-patch oriented; history/restart expose the 3-D soil-by-patch field.
#endif

      write(cyear,'(i4.4)') lc_year
      landdir = trim(dir_model_landdata) // '/soil/' // trim(cyear)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,'(/, A)') 'Aggregate lake sediment organic carbon ...'
         CALL system('mkdir -p ' // trim(adjustl(landdir)))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      lake_local = 0
      IF (p_is_worker) THEN
         DO ipatch = 1, numpatch
            IF (landpatch%settyp(ipatch) == WATERBODY) lake_local = lake_local + 1
         ENDDO
         allocate(lake_soilc_patches(nl_soil,numpatch))
         lake_soilc_patches(:,:) = 0._r8
      ENDIF
      lake_global = lake_local
#ifdef USEMPI
      CALL mpi_allreduce (lake_local, lake_global, 1, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif

      IF (lake_global > 0) THEN
         lndname = trim(dir_rawdata)//'/lake_soilc.nc'
         raw_exists = .false.
         IF (p_is_master) inquire(file=trim(lndname), exist=raw_exists)
#ifdef USEMPI
         CALL mpi_bcast (raw_exists, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif
         IF (.not. raw_exists) THEN
            lndname = trim(dir_rawdata)//'/soil/lake_soilc.nc'
            IF (p_is_master) inquire(file=trim(lndname), exist=raw_exists)
#ifdef USEMPI
            CALL mpi_bcast (raw_exists, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif
         ENDIF

         IF (raw_exists) THEN
            IF (p_is_master) write(*,*) '  using raw lake_soilc dataset: ', trim(lndname)
            IF (p_is_io) THEN
               CALL allocate_block_data (gland, lake_soilc_grid, nl_soil)
               CALL ncio_read_block (lndname, 'lake_soilc', gland, nl_soil, lake_soilc_grid)
#ifdef USEMPI
               CALL aggregation_data_daemon (gland, data_r8_3d_in1 = lake_soilc_grid, n1_r8_3d_in1 = nl_soil)
#endif
            ENDIF

            IF (p_is_worker) THEN
               DO ipatch = 1, numpatch
                  L = landpatch%settyp(ipatch)
                  IF (L == WATERBODY) THEN
                     CALL aggregation_request_data (landpatch, ipatch, gland, &
                        zip = USE_zip_for_aggregation, area = area_one, &
                        data_r8_3d_in1 = lake_soilc_grid, data_r8_3d_out1 = lake_soilc_one, &
                        n1_r8_3d_in1 = nl_soil, lb1_r8_3d_in1 = 1)

                     DO j = 1, nl_soil
                        valid_area = sum(area_one, mask = area_one > 0._r8 .and. &
                           ieee_is_finite(lake_soilc_one(j,:)) .and. lake_soilc_one(j,:) >= 0._r8 .and. &
                           abs(lake_soilc_one(j,:)) < 0.5_r8*abs(spval))
                        IF (valid_area > tiny(1._r8)) THEN
                           lake_soilc_patches(j,ipatch) = sum(lake_soilc_one(j,:)*area_one, &
                              mask = area_one > 0._r8 .and. ieee_is_finite(lake_soilc_one(j,:)) .and. &
                              lake_soilc_one(j,:) >= 0._r8 .and. &
                              abs(lake_soilc_one(j,:)) < 0.5_r8*abs(spval)) / valid_area
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO

#ifdef USEMPI
               CALL aggregation_worker_done ()
#endif
            ENDIF
         ELSE
            IF (p_is_master) THEN
               write(*,'(A)') '  raw lake_soilc.nc absent; using organic-matter proxy, not measured lake carbon.'
               write(*,'(A)') '  proxy = 580 [gC/kg OM] * OM_density [kg OM/m3]; '// &
                  'lake OM_density may include existing soil-aggregation waterbody fill.'
            ENDIF
            ! SoilParametersReadin maps raw layer 1 to model layers 1-2,
            ! layers 2-8 to 3-9, and repeats raw layer 8 below layer 9.
            ! Read each unique aggregated OM layer only once.
            n_source_layers = min(8, max(1, nl_soil-1))
            DO src_layer = 1, n_source_layers
               write(c,'(i0)') src_layer
               lndname = trim(landdir)//'/OM_density_s_l'//trim(c)//'_patches.nc'
               CALL ncio_read_vector_complete (lndname, 'OM_density_s_l'//trim(c)//'_patches', &
                  landpatch, om_density_patches)

               invalid_local = 0
               IF (p_is_worker) THEN
                  DO ipatch = 1, numpatch
                     IF (landpatch%settyp(ipatch) /= WATERBODY) CYCLE
                     om_value = om_density_patches(ipatch)
                     IF (.not. ieee_is_finite(om_value)) THEN
                        invalid_local = invalid_local + 1
                     ELSEIF (om_value < 0._r8 .or. om_value >= 0.5_r8*abs(spval)) THEN
                        invalid_local = invalid_local + 1
                     ELSE
                        DO j = 1, nl_soil
                           IF (min(8, max(1, j-1)) == src_layer) &
                              lake_soilc_patches(j,ipatch) = carbon_per_kg_om * om_value
                        ENDDO
                     ENDIF
                  ENDDO
               ENDIF
               invalid_global = invalid_local
#ifdef USEMPI
               CALL mpi_allreduce (invalid_local, invalid_global, 1, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif
               IF (invalid_global > 0) THEN
                  IF (p_is_master) write(*,'(A,I0,A,I0)') &
                     ' ***** ERROR: invalid lake OM_density_s layer ', src_layer, &
                     ' patch values: ', invalid_global
                  CALL CoLM_stop ('lake CH4 organic-matter proxy requires finite non-negative OM_density for every lake patch')
               ENDIF
            ENDDO
            IF (allocated(om_density_patches)) deallocate(om_density_patches)
         ENDIF
      ELSEIF (p_is_master) THEN
         write(*,'(A)') '  no lake patches in domain; writing zero lake sediment carbon.'
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('lake_soilc_patches [gC/m3]', lake_soilc_patches)
#endif

      lndname = trim(landdir)//'/lake_soilc_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'soil', nl_soil)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'lake_soilc_patches', 'soil', nl_soil, 'patch', &
         landpatch, lake_soilc_patches, DEF_Srfdata_CompressLevel)

      IF (p_is_worker) THEN
         deallocate(lake_soilc_patches)
      ENDIF

END SUBROUTINE Aggregation_LakeSoilC

SUBROUTINE Aggregation_MethanePH (dir_rawdata, dir_model_landdata, lc_year)

!-----------------------------------------------------------------------
! !DESCRIPTION:
!  Build a methane patch-level soil pH surface file from the native
!  PHH2O1.nc raw dataset.  This is intentionally done during mksrfdata
!  rather than at runtime: PHH2O1.nc is a 30 arc-sec global file
!  (43200 x 16800 x 4, ~2.8 GB), and broadcasting a gridded pH field
!  during colm.x initialization can deadlock or exhaust memory.
!
!  Output:
!    <dir_model_landdata>/soil/<lc_year>/methane_ph_patches.nc
!       variable methane_ph_patches(patch), pH units
!
!  The raw dataset is required only when the target domain contains soil or
!  wetland patches; domains without those patch types receive the fallback
!  value solely to preserve the surface-file schema.
!
!  Method:
!    Build a methane-private exact patch/grid intersection map, read only the
!    native cells requested by workers, and form a missing-aware mean hydrogen
!    ion activity over depth and patch area.  The final pH is -log10(activity).
!    The source depth coordinate is validated against the PHH2O1 contract, then
!    converted from layer bottoms to thickness weights for the top four layers.
!    PHH2O is declared in CDL order (depth,lat,lon); the NetCDF Fortran API
!    exposes its dimension IDs and start/count vectors as (lon,lat,depth).
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Namelist, only: DEF_Srfdata_CompressLevel
#ifdef USEMPI
   USE MOD_SPMD_Task, only: p_comm_glb, p_err, p_is_master, p_is_io, p_is_worker, &
      p_address_master, p_np_worker, p_np_io, p_address_worker, p_address_io, &
      p_stat, mpi_tag_data, MPI_LOGICAL, MPI_INTEGER, MPI_REAL8, MPI_SUM, CoLM_Stop
#else
   USE MOD_SPMD_Task, only: p_is_master, p_is_io, p_is_worker, p_np_worker, p_np_io, CoLM_Stop
#endif
   USE MOD_Const_LC, only: patchtypes
   USE MOD_LandPatch
   USE MOD_Land2mWMO, only: wmo_patch, wmo_source
   USE MOD_Grid
   USE MOD_MT_Methane_PHMapping, only: methane_ph_mapping_type, &
      build_methane_ph_areal_mapping
   USE MOD_DataType, only: pointer_real8_1d
   USE MOD_NetCDFVector
   USE netcdf
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
#ifdef RangeCheck
   USE MOD_RangeCheck
#endif

   IMPLICIT NONE

   integer, intent(in) :: lc_year
   character(len=*), intent(in) :: dir_rawdata
   character(len=*), intent(in) :: dir_model_landdata

   real(r8), parameter :: methane_ph_fallback = 6.2_r8
   real(r8), parameter :: invalid_ph = -huge(1._r8)
   real(r8), parameter :: expected_depth_bottom_cm(4) = [4.5_r8, 9.1_r8, 16.6_r8, 28.9_r8]
   integer(1), parameter :: missing_byte = -100_1

   character(len=256) :: landdir, lndname, cyear, rawfile
   character(len=64) :: depth_units
   logical :: raw_exists, requires_spatial_ph
   integer :: ncid, vid, depth_vid, ierr
   integer :: nlat, nlon, ndepth
   integer :: ipatch, invalid_local, invalid_global, relevant_local, relevant_global
   real(r8), allocatable :: methane_ph_patches(:)
   logical, allocatable :: methane_ph_valid(:)
   real(r8), allocatable :: lat_g(:), lon_g(:), depth_g(:)
   real(r8) :: dlon, depth_scale_to_cm, south_edge, north_edge
   real(r8) :: ph_layer_thickness_cm(4)
   type(grid_type) :: grid_ph
   type(methane_ph_mapping_type) :: map_ph

      write(cyear,'(i4.4)') lc_year
      landdir = trim(dir_model_landdata) // '/soil/' // trim(cyear)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,'(/, A)') 'Aggregate methane soil pH patches ...'
         CALL system('mkdir -p ' // trim(adjustl(landdir)))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      relevant_local = 0
      IF (p_is_worker) THEN
         DO ipatch = 1, numpatch
            IF (landpatch%settyp(ipatch) >= 1 .and. &
                landpatch%settyp(ipatch) <= size(patchtypes)) THEN
               requires_spatial_ph = nint(patchtypes(landpatch%settyp(ipatch))) == 0 .or. &
                  nint(patchtypes(landpatch%settyp(ipatch))) == 2
               IF (requires_spatial_ph) relevant_local = relevant_local + 1
            ENDIF
         ENDDO
         allocate(methane_ph_patches(numpatch))
         allocate(methane_ph_valid(numpatch))
         methane_ph_patches(:) = methane_ph_fallback
         methane_ph_valid(:) = .false.
      ENDIF
      relevant_global = relevant_local
#ifdef USEMPI
      CALL mpi_allreduce (relevant_local, relevant_global, 1, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif

      IF (relevant_global > 0) THEN
      rawfile = trim(dir_rawdata)//'/soil/PHH2O1.nc'
      raw_exists = .false.
      IF (p_is_master) inquire(file=trim(rawfile), exist=raw_exists)
#ifdef USEMPI
      CALL mpi_bcast (raw_exists, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif

      IF (.not. raw_exists) THEN
         CALL CoLM_Stop (' ***** ERROR: required PHH2O raw file not found: '//trim(rawfile))
      ENDIF

      nlat = 0
      nlon = 0
      ndepth = 0
      IF (p_is_master) THEN
         ierr = nf90_open(trim(rawfile), NF90_NOWRITE, ncid)
         IF (ierr /= NF90_NOERR) THEN
            write(*,'(A,A)') ' ERROR: PHH2O open failed in Aggregation_MethanePH: ', &
               trim(nf90_strerror(ierr))
            CALL CoLM_Stop (' ***** ERROR: cannot open required PHH2O source')
         ENDIF

         CALL get_dim(ncid, 'lat',   nlat)
         CALL get_dim(ncid, 'lon',   nlon)
         CALL get_dim(ncid, 'depth', ndepth)
         IF (nlat < 2 .or. nlon < 2 .or. ndepth < size(expected_depth_bottom_cm)) THEN
            write(*,'(A)') ' ERROR: PHH2O missing or invalid lat/lon/depth dimensions.'
            CALL CoLM_Stop (' ***** ERROR: invalid PHH2O dimensions')
         ENDIF

         allocate(lat_g(nlat), lon_g(nlon), depth_g(ndepth))
         ierr = nf90_inq_varid(ncid, 'lat', vid)
         IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, lat_g)
         IF (ierr /= NF90_NOERR) THEN
            write(*,'(A,A)') ' ERROR: PHH2O latitude read failed: ', trim(nf90_strerror(ierr))
            CALL CoLM_Stop (' ***** ERROR: cannot read PHH2O latitude')
         ENDIF

         ierr = nf90_inq_varid(ncid, 'lon', vid)
         IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, lon_g)
         IF (ierr /= NF90_NOERR) THEN
            write(*,'(A,A)') ' ERROR: PHH2O longitude read failed: ', trim(nf90_strerror(ierr))
            CALL CoLM_Stop (' ***** ERROR: cannot read PHH2O longitude')
         ENDIF

         ierr = nf90_inq_varid(ncid, 'depth', depth_vid)
         IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, depth_vid, depth_g)
         IF (ierr /= NF90_NOERR) THEN
            write(*,'(A,A)') ' ERROR: PHH2O depth-coordinate read failed: ', trim(nf90_strerror(ierr))
            CALL CoLM_Stop (' ***** ERROR: cannot read PHH2O depth coordinate')
         ENDIF
         depth_units = ''
         ierr = nf90_get_att(ncid, depth_vid, 'units', depth_units)
         IF (ierr /= NF90_NOERR) CALL CoLM_Stop (' ***** ERROR: PHH2O depth units are required')
         SELECT CASE (trim(adjustl(depth_units)))
         CASE ('cm', 'CM', 'centimeter', 'centimeters', 'centimetre', 'centimetres')
            depth_scale_to_cm = 1._r8
         CASE ('m', 'M', 'meter', 'meters', 'metre', 'metres')
            depth_scale_to_cm = 100._r8
         CASE ('mm', 'MM', 'millimeter', 'millimeters', 'millimetre', 'millimetres')
            depth_scale_to_cm = 0.1_r8
         CASE DEFAULT
            CALL CoLM_Stop (' ***** ERROR: unsupported PHH2O depth units: '//trim(depth_units))
         END SELECT
         depth_g = depth_g * depth_scale_to_cm
         IF (any(.not. ieee_is_finite(depth_g(1:4))) .or. &
             depth_g(1) <= 0._r8 .or. any(depth_g(2:4) <= depth_g(1:3))) &
            CALL CoLM_Stop (' ***** ERROR: PHH2O depth bottoms must be finite, positive and increasing')
         IF (any(abs(depth_g(1:4) - expected_depth_bottom_cm) > 0.05_r8)) &
            CALL CoLM_Stop (' ***** ERROR: PHH2O top-four depth coordinate is incompatible')
         ph_layer_thickness_cm(1) = depth_g(1)
         ph_layer_thickness_cm(2:4) = depth_g(2:4) - depth_g(1:3)

         ierr = nf90_inq_varid(ncid, 'PHH2O', vid)
         IF (ierr /= NF90_NOERR) THEN
            write(*,'(A,A)') ' ERROR: PHH2O variable read failed: ', trim(nf90_strerror(ierr))
            CALL CoLM_Stop (' ***** ERROR: required PHH2O variable is missing')
         ENDIF
         CALL validate_ph_variable_metadata(ncid, vid, nlon, nlat, ndepth)

         IF (any(.not. ieee_is_finite(lat_g)) .or. any(.not. ieee_is_finite(lon_g))) &
            CALL CoLM_Stop (' ***** ERROR: non-finite PHH2O coordinates')
         IF (.not. (all(lat_g(2:nlat) > lat_g(1:nlat-1)) .or. &
                    all(lat_g(2:nlat) < lat_g(1:nlat-1)))) &
            CALL CoLM_Stop (' ***** ERROR: PHH2O latitude must be strictly monotonic')
         IF (.not. all(lon_g(2:nlon) > lon_g(1:nlon-1))) &
            CALL CoLM_Stop (' ***** ERROR: PHH2O longitude must be strictly increasing')
         ! Derive the mean spacing from the full span instead of one adjacent
         ! difference.  PHH2O stores lon as float32, whose resolution near 180
         ! degrees is about 1.5e-5, so a single difference of 0.0083 keeps only
         ! three significant digits; scaling that by nlon amplifies the
         ! quantisation noise roughly 4e4 times and rejects a perfectly valid
         ! cyclic grid.  Comparing the span itself keeps the same intent without
         ! that amplification.
         dlon = (lon_g(nlon) - lon_g(1)) / real(nlon - 1, r8)
         IF (any(abs((lon_g(2:nlon) - lon_g(1:nlon-1)) - dlon) > &
                 max(1.e-6_r8, 1.e-2_r8 * dlon))) &
            CALL CoLM_Stop (' ***** ERROR: PHH2O longitude spacing must be regular')
         IF (abs(lon_g(nlon) - lon_g(1) + dlon - 360._r8) > max(1.e-4_r8, dlon)) &
            CALL CoLM_Stop (' ***** ERROR: PHH2O longitude does not cover a cyclic global grid')
         deallocate(depth_g)
         ierr = nf90_close(ncid)
         IF (ierr /= NF90_NOERR) CALL CoLM_Stop (' ***** ERROR: cannot close PHH2O metadata source')
      ENDIF

#ifdef USEMPI
      CALL mpi_bcast (nlat,   1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (nlon,   1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (ndepth, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      IF (.not. p_is_master) allocate(lat_g(nlat), lon_g(nlon))
      CALL mpi_bcast (lat_g, nlat, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (lon_g, nlon, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (ph_layer_thickness_cm, 4, MPI_REAL8, p_address_master, p_comm_glb, p_err)
#endif

      IF (lat_g(1) < lat_g(nlat)) THEN
         south_edge = lat_g(1) - 0.5_r8 * (lat_g(2) - lat_g(1))
         north_edge = lat_g(nlat) + 0.5_r8 * (lat_g(nlat) - lat_g(nlat-1))
      ELSE
         south_edge = lat_g(nlat) - 0.5_r8 * (lat_g(nlat-1) - lat_g(nlat))
         north_edge = lat_g(1) + 0.5_r8 * (lat_g(1) - lat_g(2))
      ENDIF
      south_edge = max(-90._r8, south_edge)
      north_edge = min( 90._r8, north_edge)
      CALL grid_ph%define_by_center(lat_g, lon_g, south=south_edge, north=north_edge)
      CALL build_methane_ph_areal_mapping(map_ph, grid_ph, landpatch)

      CALL aggregate_sparse_ph(rawfile, nlat, nlon, size(ph_layer_thickness_cm), &
         ph_layer_thickness_cm, map_ph, &
         methane_ph_patches, methane_ph_valid)

      deallocate(lat_g, lon_g)
      ELSEIF (p_is_master) THEN
         write(*,'(A)') '  no soil or wetland patches in domain; writing fallback pH only.'
      ENDIF

      invalid_local = 0
      IF (p_is_worker) THEN
         DO ipatch = 1, numpatch
            requires_spatial_ph = .false.
            IF (landpatch%settyp(ipatch) >= 1 .and. &
                landpatch%settyp(ipatch) <= size(patchtypes)) THEN
               requires_spatial_ph = nint(patchtypes(landpatch%settyp(ipatch))) == 0 .or. &
                  nint(patchtypes(landpatch%settyp(ipatch))) == 2
            ENDIF
            IF (requires_spatial_ph .and. .not. methane_ph_valid(ipatch)) &
               invalid_local = invalid_local + 1
         ENDDO
      ENDIF
      invalid_global = invalid_local
#ifdef USEMPI
      CALL mpi_allreduce (invalid_local, invalid_global, 1, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif
      ! PHH2O only spans latitude -56..84, so polar islands and a few remote
      ! islets carry no soil pH at all -- the abort that used to stand here was
      ! an unsatisfiable precondition for any global run.  Those patches keep
      ! the neutral default assigned at allocation instead, because disabling
      ! spatial pH globally would bias every other patch far more: pH 4.5 peat
      ! produces ~4x less methane than the pH 6.2 optimum.  The count is
      ! reported so the fallback share stays traceable in the log.
      IF (invalid_global > 0 .and. p_is_master) write(*,'(A,I0,A,I0,A,F4.1,A)') &
         '  NOTE: PHH2O has no spatial pH for ', invalid_global, ' of ', relevant_global, &
         ' soil/wetland patches (outside the -56..84 data window or all-missing); ', &
         methane_ph_fallback, ' is used there.'

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('methane_ph_patches [pH]', methane_ph_patches)
#endif

      lndname = trim(landdir)//'/methane_ph_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'methane_ph_patches', 'patch', &
         landpatch, methane_ph_patches, DEF_Srfdata_CompressLevel)

      IF (p_is_worker) THEN
         deallocate(methane_ph_patches)
         deallocate(methane_ph_valid)
      ENDIF

CONTAINS

   SUBROUTINE validate_ph_variable_metadata(ncid, vid, nlon, nlat, ndepth)
      integer, intent(in) :: ncid, vid, nlon, nlat, ndepth

      integer :: ierr_loc, ndims_loc, xtype_loc, idim, attr_value
      integer :: dimids(3), dimlen(3)
      real(r8) :: scale_factor, add_offset
      logical :: has_missing_marker
      character(len=NF90_MAX_NAME) :: dimname(3)
      character(len=64) :: ph_units

      ierr_loc = nf90_inquire_variable(ncid, vid, ndims=ndims_loc)
      IF (ierr_loc /= NF90_NOERR .or. ndims_loc /= 3) &
         CALL CoLM_Stop (' ***** ERROR: PHH2O must have exactly three dimensions')
      ierr_loc = nf90_inquire_variable(ncid, vid, xtype=xtype_loc, dimids=dimids)
      IF (ierr_loc /= NF90_NOERR) CALL CoLM_Stop (' ***** ERROR: cannot inspect PHH2O metadata')
      IF (xtype_loc /= NF90_BYTE) CALL CoLM_Stop (' ***** ERROR: PHH2O must use signed byte encoding')

      DO idim = 1, 3
         dimname(idim) = ''
         dimlen(idim) = -1
         ierr_loc = nf90_inquire_dimension(ncid, dimids(idim), name=dimname(idim), len=dimlen(idim))
         IF (ierr_loc /= NF90_NOERR) CALL CoLM_Stop (' ***** ERROR: cannot inspect PHH2O dimensions')
      ENDDO
      IF ((trim(dimname(1)) /= 'lon' .and. trim(dimname(1)) /= 'longitude') .or. &
          dimlen(1) /= nlon .or. &
          (trim(dimname(2)) /= 'lat' .and. trim(dimname(2)) /= 'latitude') .or. &
          dimlen(2) /= nlat .or. trim(dimname(3)) /= 'depth' .or. dimlen(3) /= ndepth) &
         CALL CoLM_Stop (' ***** ERROR: PHH2O Fortran dimension order must be lon,lat,depth')

      ph_units = ''
      ierr_loc = nf90_get_att(ncid, vid, 'units', ph_units)
      IF (ierr_loc /= NF90_NOERR) CALL CoLM_Stop (' ***** ERROR: PHH2O units are required')
      SELECT CASE (trim(adjustl(ph_units)))
      CASE ('1/10', '0.1', 'pH/10', 'ph/10')
         CONTINUE
      CASE DEFAULT
         CALL CoLM_Stop (' ***** ERROR: unsupported PHH2O units: '//trim(ph_units))
      END SELECT

      ierr_loc = nf90_get_att(ncid, vid, 'scale_factor', scale_factor)
      IF (ierr_loc == NF90_NOERR) THEN
         IF (.not. ieee_is_finite(scale_factor) .or. abs(scale_factor - 0.1_r8) > 1.e-12_r8) &
            CALL CoLM_Stop (' ***** ERROR: PHH2O scale_factor must be 0.1')
      ELSEIF (ierr_loc /= NF90_ENOTATT) THEN
         CALL CoLM_Stop (' ***** ERROR: cannot inspect PHH2O scale_factor')
      ENDIF
      ierr_loc = nf90_get_att(ncid, vid, 'add_offset', add_offset)
      IF (ierr_loc == NF90_NOERR) THEN
         IF (.not. ieee_is_finite(add_offset) .or. abs(add_offset) > 1.e-12_r8) &
            CALL CoLM_Stop (' ***** ERROR: PHH2O add_offset must be zero')
      ELSEIF (ierr_loc /= NF90_ENOTATT) THEN
         CALL CoLM_Stop (' ***** ERROR: cannot inspect PHH2O add_offset')
      ENDIF

      has_missing_marker = .false.
      ierr_loc = nf90_get_att(ncid, vid, 'missing_value', attr_value)
      IF (ierr_loc == NF90_NOERR) THEN
         has_missing_marker = .true.
         IF (attr_value /= int(missing_byte)) &
            CALL CoLM_Stop (' ***** ERROR: PHH2O missing_value must be -100')
      ELSEIF (ierr_loc /= NF90_ENOTATT) THEN
         CALL CoLM_Stop (' ***** ERROR: cannot inspect PHH2O missing_value')
      ENDIF
      ierr_loc = nf90_get_att(ncid, vid, '_FillValue', attr_value)
      IF (ierr_loc == NF90_NOERR) THEN
         has_missing_marker = .true.
         IF (attr_value /= int(missing_byte)) &
            CALL CoLM_Stop (' ***** ERROR: PHH2O _FillValue must be -100')
      ELSEIF (ierr_loc /= NF90_ENOTATT) THEN
         CALL CoLM_Stop (' ***** ERROR: cannot inspect PHH2O _FillValue')
      ENDIF
      IF (.not. has_missing_marker) CALL CoLM_Stop (' ***** ERROR: PHH2O missing marker is required')
   END SUBROUTINE validate_ph_variable_metadata

   SUBROUTINE aggregate_sparse_ph(rawfile, nlat, nlon, ndepth, layer_weight, mapping, patch_ph, patch_valid)
      character(len=*), intent(in) :: rawfile
      integer, intent(in) :: nlat, nlon, ndepth
      real(r8), intent(in) :: layer_weight(ndepth)
      type(methane_ph_mapping_type), intent(in) :: mapping
      real(r8), allocatable, intent(inout) :: patch_ph(:)
      logical, allocatable, intent(inout) :: patch_valid(:)

      integer :: ncid_loc, vid_loc, ierr_loc
      integer :: iproc, ng, iset, ipart, iloc, ie, wmo_src
      real(r8) :: sum_h_activity_area, valid_area_depth, value, area, depth_weight
      real(r8), allocatable :: grid_values(:), grid_weights(:)
      type(pointer_real8_1d), allocatable :: recv_values(:), recv_weights(:)

      IF (p_is_io) THEN
         ierr_loc = nf90_open(trim(rawfile), NF90_NOWRITE, ncid_loc)
         IF (ierr_loc /= NF90_NOERR) THEN
            write(*,'(A,A)') ' ERROR: PHH2O I/O-rank open failed: ', trim(nf90_strerror(ierr_loc))
            CALL CoLM_Stop (' ***** ERROR: cannot open PHH2O on I/O rank')
         ENDIF
         ierr_loc = nf90_inq_varid(ncid_loc, 'PHH2O', vid_loc)
         IF (ierr_loc /= NF90_NOERR) CALL CoLM_Stop (' ***** ERROR: PHH2O variable is missing on I/O rank')

         DO iproc = 0, p_np_worker-1
            ng = mapping%glist(iproc)%ng
            IF (ng <= 0) CYCLE
            allocate(grid_values(ng), grid_weights(ng))
            CALL read_ph_runs(ncid_loc, vid_loc, nlat, nlon, ndepth, layer_weight, &
               mapping%glist(iproc)%ilon, mapping%glist(iproc)%ilat, grid_values, grid_weights)
#ifdef USEMPI
            CALL mpi_send (grid_values, ng, MPI_REAL8, p_address_worker(iproc), &
               mpi_tag_data, p_comm_glb, p_err)
            CALL mpi_send (grid_weights, ng, MPI_REAL8, p_address_worker(iproc), &
               mpi_tag_data, p_comm_glb, p_err)
            deallocate(grid_values, grid_weights)
#endif
         ENDDO

         ierr_loc = nf90_close(ncid_loc)
         IF (ierr_loc /= NF90_NOERR) CALL CoLM_Stop (' ***** ERROR: cannot close PHH2O on I/O rank')
      ENDIF

      IF (p_is_worker) THEN
         allocate(recv_values(0:p_np_io-1), recv_weights(0:p_np_io-1))
         DO iproc = 0, p_np_io-1
            ng = mapping%glist(iproc)%ng
            IF (ng <= 0) CYCLE
            allocate(recv_values(iproc)%val(ng), recv_weights(iproc)%val(ng))
#ifdef USEMPI
            CALL mpi_recv (recv_values(iproc)%val, ng, MPI_REAL8, p_address_io(iproc), &
               mpi_tag_data, p_comm_glb, p_stat, p_err)
            CALL mpi_recv (recv_weights(iproc)%val, ng, MPI_REAL8, p_address_io(iproc), &
               mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
            recv_values(0)%val = grid_values
            recv_weights(0)%val = grid_weights
            deallocate(grid_values, grid_weights)
#endif
         ENDDO

         DO iset = 1, numpatch
            sum_h_activity_area = 0._r8
            valid_area_depth = 0._r8
            DO ipart = 1, mapping%npart(iset)
               iproc = mapping%address(iset)%val(1,ipart)
               iloc = mapping%address(iset)%val(2,ipart)
               value = recv_values(iproc)%val(iloc)
               depth_weight = recv_weights(iproc)%val(iloc)
               area = mapping%areapart(iset)%val(ipart)
               IF (value < 2._r8 .or. value > 10._r8 .or. area <= 0._r8 .or. depth_weight <= 0._r8) CYCLE
               sum_h_activity_area = sum_h_activity_area + 10._r8 ** (-value) * area * depth_weight
               valid_area_depth = valid_area_depth + area * depth_weight
            ENDDO
            IF (valid_area_depth > 0._r8) THEN
               patch_ph(iset) = -log10(sum_h_activity_area / valid_area_depth)
               patch_valid(iset) = .true.
            ENDIF
         ENDDO

         IF (allocated(wmo_patch) .and. allocated(wmo_source)) THEN
            DO iset = 1, numpatch
               ie = landpatch%ielm(iset)
               IF (ie < 1 .or. ie > size(wmo_patch)) CYCLE
               IF (iset /= wmo_patch(ie)) CYCLE
               wmo_src = wmo_source(ie)
               IF (wmo_src < 1 .or. wmo_src > numpatch) CYCLE
               patch_ph(iset) = patch_ph(wmo_src)
               patch_valid(iset) = patch_valid(wmo_src)
            ENDDO
         ENDIF

         DO iproc = 0, p_np_io-1
            IF (allocated(recv_values(iproc)%val)) deallocate(recv_values(iproc)%val)
            IF (allocated(recv_weights(iproc)%val)) deallocate(recv_weights(iproc)%val)
         ENDDO
         deallocate(recv_values, recv_weights)
      ENDIF
   END SUBROUTINE aggregate_sparse_ph

   SUBROUTINE read_ph_runs(ncid, vid, nlat, nlon, ndepth, layer_weight, ilon, ilat, values, valid_weights)
      integer, intent(in) :: ncid, vid, nlat, nlon, ndepth
      real(r8), intent(in) :: layer_weight(ndepth)
      integer, intent(in) :: ilon(:), ilat(:)
      real(r8), intent(out) :: values(:)
      real(r8), intent(out) :: valid_weights(:)

      integer :: irun, iend, i, k, d, nrun, maxrun, ierr_loc, ival
      integer(1), allocatable :: slab(:,:,:)
      real(r8) :: sum_h_activity_weight, valid_weight

      values(:) = invalid_ph
      valid_weights(:) = 0._r8
      IF (size(values) == 0) RETURN
      IF (size(ilon) /= size(values) .or. size(ilat) /= size(values) .or. &
          size(valid_weights) /= size(values)) &
         CALL CoLM_Stop (' ***** ERROR: inconsistent PHH2O sparse request')
      IF (any(ilon < 1) .or. any(ilon > nlon) .or. any(ilat < 1) .or. any(ilat > nlat)) &
         CALL CoLM_Stop (' ***** ERROR: out-of-range PHH2O sparse request')

      maxrun = 1
      irun = 1
      DO WHILE (irun <= size(values))
         iend = irun
         DO WHILE (iend < size(values))
            IF (ilat(iend+1) /= ilat(irun)) EXIT
            IF (ilon(iend+1) /= ilon(iend) + 1) EXIT
            iend = iend + 1
         ENDDO
         maxrun = max(maxrun, iend - irun + 1)
         irun = iend + 1
      ENDDO

      allocate(slab(maxrun,1,ndepth))
      irun = 1
      DO WHILE (irun <= size(values))
         iend = irun
         DO WHILE (iend < size(values))
            IF (ilat(iend+1) /= ilat(irun)) EXIT
            IF (ilon(iend+1) /= ilon(iend) + 1) EXIT
            iend = iend + 1
         ENDDO
         nrun = iend - irun + 1

         ierr_loc = nf90_get_var(ncid, vid, slab(1:nrun,1:1,1:ndepth), &
            start=[ilon(irun), ilat(irun), 1], count=[nrun, 1, ndepth])
         IF (ierr_loc /= NF90_NOERR) THEN
            write(*,'(A,A)') ' ERROR: PHH2O sparse slab read failed: ', trim(nf90_strerror(ierr_loc))
            CALL CoLM_Stop (' ***** ERROR: cannot read PHH2O sparse slab')
         ENDIF

         DO i = irun, iend
            k = i - irun + 1
            sum_h_activity_weight = 0._r8
            valid_weight = 0._r8
            DO d = 1, ndepth
               IF (slab(k,1,d) == missing_byte) CYCLE
               IF (slab(k,1,d) >= 0_1) THEN
                  ival = int(slab(k,1,d))
               ELSE
                  ival = int(slab(k,1,d)) + 256
               ENDIF
               IF (ival < 20 .or. ival > 100) CYCLE
               sum_h_activity_weight = sum_h_activity_weight + &
                  10._r8 ** (-0.1_r8 * real(ival, r8)) * layer_weight(d)
               valid_weight = valid_weight + layer_weight(d)
            ENDDO
            IF (valid_weight > 0._r8) THEN
               values(i) = -log10(sum_h_activity_weight / valid_weight)
               valid_weights(i) = valid_weight
            ENDIF
         ENDDO

         irun = iend + 1
      ENDDO
      deallocate(slab)
   END SUBROUTINE read_ph_runs

   SUBROUTINE get_dim(ncid, name, n)
      integer, intent(in)  :: ncid
      character(len=*), intent(in) :: name
      integer, intent(out) :: n
      integer :: did, ierr_loc

      n = 0
      ierr_loc = nf90_inq_dimid(ncid, trim(name), did)
      IF (ierr_loc /= NF90_NOERR) RETURN
      ierr_loc = nf90_inquire_dimension(ncid, did, len=n)
      IF (ierr_loc /= NF90_NOERR) n = 0
   END SUBROUTINE get_dim

END SUBROUTINE Aggregation_MethanePH

END MODULE MOD_MT_Methane_Aggregation
#endif
