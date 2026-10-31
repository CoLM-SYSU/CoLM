#include <define.h>

#if (defined MaterialTransport) && (defined BGC)

MODULE MOD_MT_Methane_PHMapping

!-----------------------------------------------------------------------
! Methane-private sparse areal mapping used while aggregating PHH2O.
!
! The general spatial-mapping builder also constructs a grid-wide sum-area
! field.  PHH2O needs only the patch/grid intersections and sparse request
! lists, so this narrow builder deliberately owns that methane-specific path
! and never allocates a 43200 x 16800 real(r8) grid field.
!-----------------------------------------------------------------------

   USE MOD_Precision, only: r8
   USE MOD_Grid, only: grid_type, grid_list_type
   USE MOD_DataType, only: pointer_int32_1d, pointer_int32_2d, pointer_real8_1d
   IMPLICIT NONE
   PRIVATE

   type, PUBLIC :: methane_ph_mapping_type
      type(grid_list_type), allocatable :: glist(:)
      integer, allocatable :: npart(:)
      type(pointer_int32_2d), allocatable :: address(:)
      type(pointer_real8_1d), allocatable :: areapart(:)
   END type methane_ph_mapping_type

   PUBLIC :: build_methane_ph_areal_mapping

CONTAINS

   SUBROUTINE build_methane_ph_areal_mapping(this, fgrid, pixelset)
      USE MOD_Namelist, only: SITE_lat_location, SITE_lon_location
      USE MOD_Block, only: gblock
      USE MOD_Pixel, only: pixel
      USE MOD_Pixelset, only: pixelset_type
      USE MOD_Mesh, only: mesh
      USE MOD_Utils, only: areaquad, expand_list, find_in_sorted_list2, &
         find_nearest_east, find_nearest_north, find_nearest_south, find_nearest_west, &
         insert_into_sorted_list1, insert_into_sorted_list2, lon_between_ceil, lon_between_floor
#ifdef USEMPI
      USE MOD_SPMD_Task, only: p_address_io, p_comm_glb, p_err, p_iam_glb, p_is_io, &
         p_is_master, p_is_worker, p_itis_io, p_itis_worker, p_np_io, p_np_worker, &
         p_stat, mpi_tag_data, mpi_tag_mesg, MPI_ANY_SOURCE, MPI_INTEGER, CoLM_Stop
#else
      USE MOD_SPMD_Task, only: p_is_io, p_is_master, p_is_worker, p_np_io, p_np_worker, CoLM_Stop
#endif
      IMPLICIT NONE

      type(methane_ph_mapping_type), intent(inout) :: this
      type(grid_type), intent(in) :: fgrid
      type(pixelset_type), intent(in) :: pixelset

      type(pointer_real8_1d), allocatable :: afrac(:)
      type(grid_list_type), allocatable :: gfrom(:)
      type(pointer_int32_1d), allocatable :: list_lat(:)
      integer, allocatable :: ng_lat(:), ys(:), yn(:), xw(:), xe(:)
      integer, allocatable :: xlist(:), ylist(:)
#ifdef USEMPI
      integer, allocatable :: ipt(:)
      logical, allocatable :: msk(:)
#endif
      integer :: ie, iset, ng, ig, ng_all, iloc, npxl, capacity
      integer :: ipxl, ilat, ilon, iworker, iproc, idest, isrc, nrecv
      integer :: iy, ix, xblk, yblk, ipxstt, ipxend
#ifdef USEMPI
      integer :: rmesg(2), smesg(2)
#endif
      real(r8) :: lat_s, lat_n, lon_w, lon_e, area
      logical :: skip, is_new

      IF (fgrid%nlat <= 0 .or. fgrid%nlon <= 0) &
         CALL CoLM_Stop(' ***** ERROR: invalid PHH2O grid in methane pH mapping')
      IF (p_np_io <= 0 .or. p_np_worker <= 0) &
         CALL CoLM_Stop(' ***** ERROR: methane pH mapping requires I/O and worker tasks')
      IF (allocated(this%glist) .or. allocated(this%npart) .or. &
          allocated(this%address) .or. allocated(this%areapart)) &
         CALL CoLM_Stop(' ***** ERROR: methane pH mapping object cannot be rebuilt')

      IF (p_is_master) THEN
         write(*,"(A, I0, A, I0, A)") &
            'Making methane pH sparse mapping: ', fgrid%nlat, &
            ' grids in latitude ', fgrid%nlon, ' grids in longitude.'

#ifndef SinglePoint
         IF (.not. (lon_between_floor(pixel%edgew, fgrid%lon_w(1), fgrid%lon_e(fgrid%nlon)) &
            .and. lon_between_ceil(pixel%edgee, fgrid%lon_w(1), fgrid%lon_e(fgrid%nlon)))) THEN
            write(*,'(A)') 'Warning: PHH2O grid does not cover the modeling longitude range.'
         ENDIF

         IF (fgrid%yinc == 1) THEN
            IF (.not. (pixel%edges >= fgrid%lat_s(1) .and. &
                       pixel%edgen <= fgrid%lat_n(fgrid%nlat))) &
               write(*,'(A)') 'Warning: PHH2O grid does not cover the modeling latitude range.'
         ELSE
            IF (.not. (pixel%edges >= fgrid%lat_s(fgrid%nlat) .and. &
                       pixel%edgen <= fgrid%lat_n(1))) &
               write(*,'(A)') 'Warning: PHH2O grid does not cover the modeling latitude range.'
         ENDIF
#endif
      ENDIF

#ifdef SinglePoint
      allocate(this%glist(0:0))
      allocate(this%glist(0)%ilat(1), this%glist(0)%ilon(1))
      this%glist(0)%ng = 1
      this%glist(0)%ilat(1) = find_nearest_south(SITE_lat_location, fgrid%nlat, fgrid%lat_s)
      this%glist(0)%ilon(1) = find_nearest_west(SITE_lon_location, fgrid%nlon, fgrid%lon_w)

      allocate(this%npart(pixelset%nset))
      allocate(this%address(pixelset%nset), this%areapart(pixelset%nset))
      this%npart(:) = 1
      DO iset = 1, pixelset%nset
         allocate(this%address(iset)%val(2,1), this%areapart(iset)%val(1))
         this%address(iset)%val = reshape([0, 1], [2, 1])
         this%areapart(iset)%val = 1._r8
      ENDDO
      RETURN
#endif

      IF (p_is_worker) THEN
         allocate(afrac(pixelset%nset), gfrom(pixelset%nset))
         allocate(ys(pixel%nlat), yn(pixel%nlat), xw(pixel%nlon), xe(pixel%nlon))

         DO ilat = 1, pixel%nlat
            ys(ilat) = find_nearest_south(pixel%lat_s(ilat), fgrid%nlat, fgrid%lat_s)
            yn(ilat) = find_nearest_north(pixel%lat_n(ilat), fgrid%nlat, fgrid%lat_n)
         ENDDO
         DO ilon = 1, pixel%nlon
            xw(ilon) = find_nearest_west(pixel%lon_w(ilon), fgrid%nlon, fgrid%lon_w)
            xe(ilon) = find_nearest_east(pixel%lon_e(ilon), fgrid%nlon, fgrid%lon_e)
         ENDDO

         allocate(list_lat(fgrid%nlat), ng_lat(fgrid%nlat))
         ng_lat(:) = 0
         capacity = min(100, max(1, fgrid%nlon))
         DO iy = 1, fgrid%nlat
            allocate(list_lat(iy)%val(capacity))
         ENDDO

         DO iset = 1, pixelset%nset
            ie = pixelset%ielm(iset)
            IF (ie < 1 .or. ie > size(mesh)) &
               CALL CoLM_Stop(' ***** ERROR: invalid element in methane pH mapping')

            ipxstt = pixelset%ipxstt(iset)
            ipxend = pixelset%ipxend(iset)
            IF (ipxstt == -1 .and. ipxend == -1) THEN
               ipxstt = 1
               ipxend = mesh(ie)%npxl
            ENDIF
            IF (ipxstt < 1 .or. ipxend < ipxstt .or. ipxend > mesh(ie)%npxl) &
               CALL CoLM_Stop(' ***** ERROR: invalid patch pixels in methane pH mapping')

            npxl = ipxend - ipxstt + 1
            capacity = max(1, npxl)
            allocate(afrac(iset)%val(capacity))
            allocate(gfrom(iset)%ilat(capacity), gfrom(iset)%ilon(capacity))
            gfrom(iset)%ng = 0

            DO ipxl = ipxstt, ipxend
               ilat = mesh(ie)%ilat(ipxl)
               ilon = mesh(ie)%ilon(ipxl)
               IF (ilat < 1 .or. ilat > pixel%nlat .or. ilon < 1 .or. ilon > pixel%nlon) &
                  CALL CoLM_Stop(' ***** ERROR: invalid mesh pixel in methane pH mapping')

               DO iy = ys(ilat), yn(ilat), fgrid%yinc
                  lat_s = max(fgrid%lat_s(iy), pixel%lat_s(ilat))
                  lat_n = min(fgrid%lat_n(iy), pixel%lat_n(ilat))
                  IF (lat_n - lat_s < 1.e-6_r8) CYCLE

                  ix = xw(ilon)
                  DO
                     IF (ix == xw(ilon)) THEN
                        lon_w = pixel%lon_w(ilon)
                     ELSE
                        lon_w = fgrid%lon_w(ix)
                     ENDIF
                     IF (ix == xe(ilon)) THEN
                        lon_e = pixel%lon_e(ilon)
                     ELSE
                        lon_e = fgrid%lon_e(ix)
                     ENDIF

                     skip = .not. (lon_between_floor(lon_w, pixel%lon_w(ilon), lon_e) .and. &
                                    lon_between_ceil(lon_e, lon_w, pixel%lon_e(ilon)))
                     IF (.not. skip) THEN
                        IF (lon_e > lon_w) THEN
                           skip = lon_e - lon_w < 1.e-6_r8
                        ELSE
                           skip = lon_e + 360._r8 - lon_w < 1.e-6_r8
                        ENDIF
                     ENDIF

                     IF (.not. skip) THEN
                        area = areaquad(lat_s, lat_n, lon_w, lon_e)

                        IF (gfrom(iset)%ng == size(gfrom(iset)%ilat)) THEN
                           CALL expand_list(gfrom(iset)%ilat, 0.2_r8)
                           CALL expand_list(gfrom(iset)%ilon, 0.2_r8)
                           CALL expand_list(afrac(iset)%val, 0.2_r8)
                        ENDIF
                        CALL insert_into_sorted_list2(ix, iy, gfrom(iset)%ng, &
                           gfrom(iset)%ilon, gfrom(iset)%ilat, iloc, is_new)
                        IF (is_new) THEN
                           IF (iloc < gfrom(iset)%ng) afrac(iset)%val(iloc+1:gfrom(iset)%ng) = &
                              afrac(iset)%val(iloc:gfrom(iset)%ng-1)
                           afrac(iset)%val(iloc) = area
                        ELSE
                           afrac(iset)%val(iloc) = afrac(iset)%val(iloc) + area
                        ENDIF

                        IF (ng_lat(iy) == size(list_lat(iy)%val)) &
                           CALL expand_list(list_lat(iy)%val, 0.2_r8)
                        CALL insert_into_sorted_list1(ix, ng_lat(iy), list_lat(iy)%val, iloc)
                     ENDIF

                     IF (ix == xe(ilon)) EXIT
                     ix = mod(ix, fgrid%nlon) + 1
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         deallocate(ys, yn, xw, xe)
         ng_all = sum(ng_lat)
         allocate(xlist(ng_all), ylist(ng_all))
         ig = 0
         DO iy = 1, fgrid%nlat
            DO ix = 1, ng_lat(iy)
               ig = ig + 1
               xlist(ig) = list_lat(iy)%val(ix)
               ylist(ig) = iy
            ENDDO
         ENDDO
         deallocate(ng_lat, list_lat)

#ifdef USEMPI
         allocate(ipt(ng_all), msk(ng_all))
         DO ig = 1, ng_all
            xblk = fgrid%xblk(xlist(ig))
            yblk = fgrid%yblk(ylist(ig))
            ipt(ig) = gblock%pio(xblk,yblk)
         ENDDO
#endif

         allocate(this%glist(0:p_np_io-1))
         this%glist(:)%ng = 0
         DO iproc = 0, p_np_io-1
#ifdef USEMPI
            msk = ipt == p_address_io(iproc)
            ng = count(msk)
#else
            ng = ng_all
#endif
            this%glist(iproc)%ng = ng
            IF (ng > 0) THEN
               allocate(this%glist(iproc)%ilat(ng), this%glist(iproc)%ilon(ng))
#ifdef USEMPI
               this%glist(iproc)%ilon = pack(xlist, msk)
               this%glist(iproc)%ilat = pack(ylist, msk)
#else
               this%glist(iproc)%ilon = xlist
               this%glist(iproc)%ilat = ylist
#endif
            ENDIF
         ENDDO

#ifdef USEMPI
         deallocate(ipt, msk)
#endif
         allocate(this%npart(pixelset%nset))
         allocate(this%address(pixelset%nset), this%areapart(pixelset%nset))
         DO iset = 1, pixelset%nset
            ng = gfrom(iset)%ng
            this%npart(iset) = ng
            allocate(this%address(iset)%val(2,ng), this%areapart(iset)%val(ng))
            IF (ng > 0) this%areapart(iset)%val = afrac(iset)%val(1:ng)
            IF (pixelset%has_shared .and. ng > 0) &
               this%areapart(iset)%val = this%areapart(iset)%val * pixelset%pctshared(iset)

            DO ig = 1, ng
               ilon = gfrom(iset)%ilon(ig)
               ilat = gfrom(iset)%ilat(ig)
               xblk = fgrid%xblk(ilon)
               yblk = fgrid%yblk(ilat)
#ifdef USEMPI
               iproc = p_itis_io(gblock%pio(xblk,yblk))
#else
               iproc = 0
#endif
               IF (iproc < 0 .or. iproc >= p_np_io) &
                  CALL CoLM_Stop(' ***** ERROR: invalid I/O owner in methane pH mapping')
               iloc = find_in_sorted_list2(ilon, ilat, this%glist(iproc)%ng, &
                  this%glist(iproc)%ilon, this%glist(iproc)%ilat)
               IF (iloc < 1 .or. iloc > this%glist(iproc)%ng) &
                  CALL CoLM_Stop(' ***** ERROR: inconsistent sparse methane pH mapping')
               this%address(iset)%val(1,ig) = iproc
               this%address(iset)%val(2,ig) = iloc
            ENDDO
         ENDDO

         deallocate(xlist, ylist, afrac, gfrom)
      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN
         DO iproc = 0, p_np_io-1
            idest = p_address_io(iproc)
            smesg = [p_iam_glb, this%glist(iproc)%ng]
            CALL mpi_send(smesg, 2, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)
            IF (this%glist(iproc)%ng > 0) THEN
               CALL mpi_send(this%glist(iproc)%ilon, this%glist(iproc)%ng, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
               CALL mpi_send(this%glist(iproc)%ilat, this%glist(iproc)%ng, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
            ENDIF
         ENDDO
      ENDIF

      IF (p_is_io) THEN
         allocate(this%glist(0:p_np_worker-1))
         this%glist(:)%ng = 0
         DO iworker = 0, p_np_worker-1
            CALL mpi_recv(rmesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, mpi_tag_mesg, &
               p_comm_glb, p_stat, p_err)
            isrc = rmesg(1)
            nrecv = rmesg(2)
            iproc = p_itis_worker(isrc)
            IF (iproc < 0 .or. iproc >= p_np_worker .or. nrecv < 0) &
               CALL CoLM_Stop(' ***** ERROR: invalid worker request in methane pH mapping')
            this%glist(iproc)%ng = nrecv
            IF (nrecv > 0) THEN
               allocate(this%glist(iproc)%ilon(nrecv), this%glist(iproc)%ilat(nrecv))
               CALL mpi_recv(this%glist(iproc)%ilon, nrecv, MPI_INTEGER, isrc, mpi_tag_data, &
                  p_comm_glb, p_stat, p_err)
               CALL mpi_recv(this%glist(iproc)%ilat, nrecv, MPI_INTEGER, isrc, mpi_tag_data, &
                  p_comm_glb, p_stat, p_err)
            ENDIF
         ENDDO
      ENDIF
#endif

   END SUBROUTINE build_methane_ph_areal_mapping

END MODULE MOD_MT_Methane_PHMapping
#endif
