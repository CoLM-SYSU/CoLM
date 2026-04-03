#include <define.h>

MODULE MOD_SpatialMapping

!--------------------------------------------------------------------------------
! !DESCRIPTION:
!
!    Spatial Mapping module.
!
!  Created by Shupeng Zhang, May 2024
!--------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   ! ------
   type :: spatial_mapping_type

      type(grid_type) :: grid

      type(grid_list_type), allocatable :: glist (:)

      integer :: npset
      integer, allocatable :: npart(:)
      type(pointer_int32_2d), allocatable :: address (:)

      logical  :: has_missing_value = .false.
      real(r8) :: missing_value     = spval

      type(pointer_real8_1d), allocatable :: areapart(:) ! intersection area
      real(r8), allocatable               :: areapset(:)
      type(block_data_real8_2d)           :: areagrid

   CONTAINS

      procedure, PUBLIC :: build_arealweighted => spatial_mapping_build_arealweighted
      procedure, PUBLIC :: build_bilinear      => spatial_mapping_build_bilinear

      procedure, PUBLIC :: set_missing_value   => spatial_mapping_set_missing_value

      ! 1) from pixelset to grid
      procedure, PRIVATE :: pset2grid_2d => spatial_mapping_pset2grid_2d
      procedure, PRIVATE :: pset2grid_3d => spatial_mapping_pset2grid_3d
      procedure, PRIVATE :: pset2grid_4d => spatial_mapping_pset2grid_4d
      generic,   PUBLIC  :: pset2grid    => pset2grid_2d, pset2grid_3d, pset2grid_4d

      procedure, PUBLIC  :: pset2grid_max   => spatial_mapping_pset2grid_max
      procedure, PUBLIC  :: pset2grid_split => spatial_mapping_pset2grid_split

      procedure, PUBLIC  :: get_sumarea  => spatial_mapping_get_sumarea

      ! 2) from grid to pixelset
      procedure, PRIVATE :: grid2pset_2d => spatial_mapping_grid2pset_2d
      procedure, PRIVATE :: grid2pset_3d => spatial_mapping_grid2pset_3d
      generic,   PUBLIC  :: grid2pset    => grid2pset_2d, grid2pset_3d

      procedure, PUBLIC  :: grid2pset_dominant => spatial_mapping_dominant_2d
      procedure, PUBLIC  :: grid2pset_varvalue => spatial_mapping_varvalue_2d

      ! 3) between grid and intersections
      procedure, PUBLIC  :: grid2part => spatial_mapping_grid2part
      procedure, PUBLIC  :: part2grid => spatial_mapping_part2grid
      procedure, PUBLIC  :: normalize => spatial_mapping_normalize

      ! 4) intersections to pixelset
      procedure, PUBLIC  :: part2pset => spatial_mapping_part2pset

      procedure, PUBLIC  :: allocate_part => spatial_mapping_allocate_part
      procedure, PUBLIC  :: deallocate_part => spatial_mapping_deallocate_part
      procedure, PUBLIC  :: forc_free_mem => forc_free_mem_spatial_mapping

      final :: spatial_mapping_free_mem

   END type spatial_mapping_type

!-----------------------
CONTAINS

   !------------------------------------------
   SUBROUTINE spatial_mapping_build_arealweighted (this, fgrid, pixelset)

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_Block
   USE MOD_Pixel
   USE MOD_Grid
   USE MOD_Pixelset
   USE MOD_DataType
   USE MOD_Mesh
   USE MOD_Utils
   USE MOD_SPMD_Task
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   type(grid_type),     intent(in) :: fgrid
   type(pixelset_type), intent(in) :: pixelset

   ! Local variables
   type(pointer_real8_1d), allocatable :: afrac(:)
   type(grid_list_type),   allocatable :: gfrom(:)
   type(pointer_int32_1d), allocatable :: list_lat(:)
   integer,  allocatable :: ng_lat(:)
   integer,  allocatable :: ys(:), yn(:), xw(:), xe(:)
   integer,  allocatable :: xlist(:), ylist(:)
   integer,  allocatable :: ipt(:)
   logical,  allocatable :: msk(:)

   integer  :: ie, iset, iblkme
   integer  :: ng, ig, ng_all, iloc
   integer  :: npxl, ipxl, ilat, ilon
   integer  :: iworker, iproc, idest, isrc, nrecv
   integer  :: rmesg(2), smesg(2)
   integer  :: iy, ix, xblk, yblk, xloc, yloc
   integer  :: ipxstt, ipxend
   real(r8) :: lat_s, lat_n, lon_w, lon_e, area
   logical  :: skip, is_new


#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN

         write(*,"(A, I0, A, I0, A)") &
            'Making areal weighted mapping between pixel set and grid: ', &
            fgrid%nlat, ' grids in latitude ', fgrid%nlon, ' grids in longitude.'

#ifndef SinglePoint
         IF (.not. (lon_between_floor(pixel%edgew, fgrid%lon_w(1), fgrid%lon_e(fgrid%nlon)) &
            .and. lon_between_ceil(pixel%edgee, fgrid%lon_w(1), fgrid%lon_e(fgrid%nlon)))) THEN
            write(*,'(A)') 'Warning: Grid does not cover longitude range of modeling region.'
         ENDIF

         IF (fgrid%yinc == 1) THEN
            IF (.not. ((pixel%edges >= fgrid%lat_s(1)) &
               .and. (pixel%edgen <= fgrid%lat_n(fgrid%nlat)))) THEN
               write(*,'(A)') 'Warning: Grid does not cover latitude range of modeling region.'
            ENDIF
         ELSE
            IF (.not. ((pixel%edges >= fgrid%lat_s(fgrid%nlat)) &
               .and. (pixel%edgen <= fgrid%lat_n(1)))) THEN
               write(*,'(A)') 'Warning: Grid does not cover latitude range of modeling region.'
            ENDIF
         ENDIF
#endif

      ENDIF

      allocate (this%grid%xblk (size(fgrid%xblk)));  this%grid%xblk = fgrid%xblk
      allocate (this%grid%yblk (size(fgrid%yblk)));  this%grid%yblk = fgrid%yblk
      allocate (this%grid%xloc (size(fgrid%xloc)));  this%grid%xloc = fgrid%xloc
      allocate (this%grid%yloc (size(fgrid%yloc)));  this%grid%yloc = fgrid%yloc
      allocate (this%grid%xcnt (size(fgrid%xcnt)));  this%grid%xcnt = fgrid%xcnt
      allocate (this%grid%ycnt (size(fgrid%ycnt)));  this%grid%ycnt = fgrid%ycnt

#ifdef SinglePoint
      allocate (this%glist (0:0))
      allocate (this%glist(0)%ilat (1))
      allocate (this%glist(0)%ilon (1))

      allocate (this%npart   (pixelset%nset))
      allocate (this%address (pixelset%nset))
      allocate (this%areapset(pixelset%nset))
      allocate (this%areapart(pixelset%nset))
      DO iset = 1, pixelset%nset
         allocate (this%address(iset)%val (2,1))
         allocate (this%areapart(iset)%val  (1))
      ENDDO

      this%glist(0)%ng = 1
      this%glist(0)%ilat(1) = find_nearest_south (SITE_lat_location, fgrid%nlat, fgrid%lat_s)
      this%glist(0)%ilon(1) = find_nearest_west  (SITE_lon_location, fgrid%nlon, fgrid%lon_w)

      this%npset = pixelset%nset
      this%npart   (:) = 1
      this%areapset(:) = 1.

      DO iset = 1, pixelset%nset
         this%address(iset)%val  = reshape((/0,1/), (/2,1/))
         this%areapart(iset)%val = 1.
      ENDDO

      CALL allocate_block_data (fgrid, this%areagrid)
      DO iblkme = 1, gblock%nblkme
         xblk = gblock%xblkme(iblkme)
         yblk = gblock%yblkme(iblkme)
         this%areagrid%blk(xblk,yblk)%val = 1.
      ENDDO

      RETURN
#endif


      IF (p_is_worker) THEN

         this%npset = pixelset%nset

         allocate (afrac (pixelset%nset))
         allocate (gfrom (pixelset%nset))

         allocate (ys (pixel%nlat))
         allocate (yn (pixel%nlat))
         allocate (xw (pixel%nlon))
         allocate (xe (pixel%nlon))

         DO ilat = 1, pixel%nlat
            ys(ilat) = find_nearest_south (pixel%lat_s(ilat), fgrid%nlat, fgrid%lat_s)
            yn(ilat) = find_nearest_north (pixel%lat_n(ilat), fgrid%nlat, fgrid%lat_n)
         ENDDO

         DO ilon = 1, pixel%nlon
            xw(ilon) = find_nearest_west (pixel%lon_w(ilon), fgrid%nlon, fgrid%lon_w)
            xe(ilon) = find_nearest_east (pixel%lon_e(ilon), fgrid%nlon, fgrid%lon_e)
         ENDDO

         allocate (list_lat (fgrid%nlat))
         DO iy = 1, fgrid%nlat
            allocate (list_lat(iy)%val (100))
         ENDDO

         allocate (ng_lat (fgrid%nlat)); ng_lat(:) = 0

         DO iset = 1, pixelset%nset

            ie = pixelset%ielm(iset)
            npxl = pixelset%ipxend(iset) - pixelset%ipxstt(iset) + 1

            ipxstt = pixelset%ipxstt(iset)
            ipxend = pixelset%ipxend(iset)

            ! deal with 2m WMO patch
            IF (ipxstt==-1 .and. ipxend==-1) THEN
               ipxstt = 1
               ipxend = mesh(ie)%npxl
               npxl   = mesh(ie)%npxl
            ENDIF

            allocate (afrac(iset)%val (npxl))
            allocate (gfrom(iset)%ilat(npxl))
            allocate (gfrom(iset)%ilon(npxl))

            gfrom(iset)%ng = 0

            DO ipxl = ipxstt, ipxend

               ilat = mesh(ie)%ilat(ipxl)
               ilon = mesh(ie)%ilon(ipxl)

               DO iy = ys(ilat), yn(ilat), fgrid%yinc

                  lat_s = max(fgrid%lat_s(iy), pixel%lat_s(ilat))
                  lat_n = min(fgrid%lat_n(iy), pixel%lat_n(ilat))

                  IF ((lat_n-lat_s) < 1.0e-6_r8) THEN
                     CYCLE
                  ENDIF

                  ix = xw(ilon)
                  DO WHILE (.true.)

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

                     skip = .false.
                     IF (.not. (lon_between_floor (lon_w, pixel%lon_w(ilon), lon_e) &
                        .and. lon_between_ceil (lon_e, lon_w, pixel%lon_e(ilon)))) THEN
                        skip = .true.
                     ELSE
                        IF (lon_e > lon_w) THEN
                           IF ((lon_e-lon_w) < 1.0e-6_r8) THEN
                              skip = .true.
                           ENDIF
                        ELSE
                           IF ((lon_e+360.0_r8-lon_w) < 1.0e-6_r8) THEN
                              skip = .true.
                           ENDIF
                        ENDIF
                     ENDIF

                     IF (.not. skip) THEN

                        area = areaquad (lat_s, lat_n, lon_w, lon_e)

                        CALL insert_into_sorted_list2 ( ix, iy, &
                           gfrom(iset)%ng, gfrom(iset)%ilon, gfrom(iset)%ilat, &
                           iloc, is_new)

                        IF (is_new) THEN
                           IF (iloc < gfrom(iset)%ng) THEN
                              afrac(iset)%val(iloc+1:gfrom(iset)%ng) &
                                 = afrac(iset)%val(iloc:gfrom(iset)%ng-1)
                           ENDIF

                           afrac(iset)%val(iloc) = area
                        ELSE
                           afrac(iset)%val(iloc) = afrac(iset)%val(iloc) + area
                        ENDIF

                        IF (gfrom(iset)%ng == size(gfrom(iset)%ilat)) THEN
                           CALL expand_list (gfrom(iset)%ilat, 0.2_r8)
                           CALL expand_list (gfrom(iset)%ilon, 0.2_r8)
                           CALL expand_list (afrac(iset)%val,  0.2_r8)
                        ENDIF

                        CALL insert_into_sorted_list1 ( &
                           ix, ng_lat(iy), list_lat(iy)%val, iloc)

                        IF (ng_lat(iy) == size(list_lat(iy)%val)) THEN
                           CALL expand_list (list_lat(iy)%val, 0.2_r8)
                        ENDIF

                     ENDIF

                     IF (ix == xe(ilon))  EXIT
                     ix = mod(ix,fgrid%nlon) + 1
                  ENDDO
               ENDDO

            ENDDO
         ENDDO

         deallocate (ys)
         deallocate (yn)
         deallocate (xw)
         deallocate (xe)

         ng_all = sum(ng_lat)
         allocate (xlist(ng_all))
         allocate (ylist(ng_all))

         ig = 0
         DO iy = 1, fgrid%nlat
            DO ix = 1, ng_lat(iy)
               ig = ig + 1
               xlist(ig) = list_lat(iy)%val(ix)
               ylist(ig) = iy
            ENDDO
         ENDDO

         deallocate (ng_lat)
         DO iy = 1, fgrid%nlat
            deallocate (list_lat(iy)%val)
         ENDDO
         deallocate (list_lat)

#ifdef USEMPI
         allocate (ipt (ng_all))
         allocate (msk (ng_all))
         DO ig = 1, ng_all
            xblk = fgrid%xblk(xlist(ig))
            yblk = fgrid%yblk(ylist(ig))
            ipt(ig) = gblock%pio(xblk,yblk)
         ENDDO
#endif

         allocate (this%glist (0:p_np_io-1))
         DO iproc = 0, p_np_io-1
#ifdef USEMPI
            msk = (ipt == p_address_io(iproc))
            ng  = count(msk)
#else
            ng  = ng_all
#endif

            this%glist(iproc)%ng = ng

            IF (ng > 0) THEN
               allocate (this%glist(iproc)%ilat (ng))
               allocate (this%glist(iproc)%ilon (ng))

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
         deallocate (ipt)
         deallocate (msk)
#endif

         allocate (this%address  (pixelset%nset))
         allocate (this%areapart (pixelset%nset))

         allocate (this%npart (pixelset%nset))

         DO iset = 1, pixelset%nset

            ng = gfrom(iset)%ng

            this%npart(iset) = ng

            allocate (this%address(iset)%val (2,ng))
            allocate (this%areapart(iset)%val (ng))

            this%areapart(iset)%val = afrac(iset)%val(1:ng)

            IF (pixelset%has_shared) THEN
               this%areapart(iset)%val = this%areapart(iset)%val * pixelset%pctshared(iset)
            ENDIF

            DO ig = 1, gfrom(iset)%ng
               ilon = gfrom(iset)%ilon(ig)
               ilat = gfrom(iset)%ilat(ig)
               xblk = fgrid%xblk(ilon)
               yblk = fgrid%yblk(ilat)

#ifdef USEMPI
               iproc = p_itis_io(gblock%pio(xblk,yblk))
#else
               iproc = 0
#endif

               this%address(iset)%val(1,ig) = iproc
               this%address(iset)%val(2,ig) = find_in_sorted_list2 ( &
                  ilon, ilat, this%glist(iproc)%ng, this%glist(iproc)%ilon, this%glist(iproc)%ilat)
            ENDDO
         ENDDO

         deallocate (xlist)
         deallocate (ylist)

         DO iset = 1, pixelset%nset
            deallocate (afrac(iset)%val )
            deallocate (gfrom(iset)%ilon)
            deallocate (gfrom(iset)%ilat)
         ENDDO

         deallocate (afrac)
         deallocate (gfrom)

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN

         DO iproc = 0, p_np_io-1
            idest = p_address_io(iproc)
            smesg = (/p_iam_glb, this%glist(iproc)%ng/)

            CALL mpi_send (smesg, 2, MPI_INTEGER, &
               idest, mpi_tag_mesg, p_comm_glb, p_err)

            IF (this%glist(iproc)%ng > 0) THEN
               CALL mpi_send (this%glist(iproc)%ilon, this%glist(iproc)%ng, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
               CALL mpi_send (this%glist(iproc)%ilat, this%glist(iproc)%ng, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
            ENDIF
         ENDDO

      ENDIF

      IF (p_is_io) THEN

         allocate (this%glist (0:p_np_worker-1))

         DO iworker = 0, p_np_worker-1

            CALL mpi_recv (rmesg, 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = rmesg(1)
            nrecv = rmesg(2)
            iproc = p_itis_worker(isrc)

            this%glist(iproc)%ng = nrecv

            IF (nrecv > 0) THEN
               allocate (this%glist(iproc)%ilon (nrecv))
               allocate (this%glist(iproc)%ilat (nrecv))

               CALL mpi_recv (this%glist(iproc)%ilon, nrecv, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL mpi_recv (this%glist(iproc)%ilat, nrecv, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
            ENDIF
         ENDDO

      ENDIF
#endif

      IF (p_is_worker) THEN
         IF (this%npset > 0) THEN
            allocate (this%areapset (this%npset))
            this%areapset(:) = 0.
         ENDIF
         DO iset = 1, this%npset
            IF (this%npart(iset) > 0) THEN
               this%areapset(iset) = sum(this%areapart(iset)%val)
            ENDIF
         ENDDO
      ENDIF

      IF (p_is_io) CALL allocate_block_data (fgrid, this%areagrid)
      IF (p_is_worker) THEN
         IF (this%npset > 0) THEN
            allocate (msk (this%npset))
            msk = pixelset%ipxstt > 0 .and. pixelset%ipxend > 0
         ENDIF
      ENDIF
      CALL this%get_sumarea (this%areagrid, msk)

      IF (allocated(msk)) deallocate(msk)


#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE spatial_mapping_build_arealweighted

   !------------------------------------------
   SUBROUTINE spatial_mapping_build_bilinear (this, fgrid, pixelset)

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_Block
   USE MOD_Pixel
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_Mesh
   USE MOD_Pixelset
   USE MOD_Utils
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only: pi
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   type(grid_type),     intent(in) :: fgrid
   type(pixelset_type), intent(in) :: pixelset

   ! Local variables
   integer,  allocatable :: ys(:), yn(:), xw(:), xe(:)
   integer,  allocatable :: xlist(:), ylist(:), ipt(:)

   real(r8), allocatable :: rlon_pset(:), rlat_pset(:)
   real(r8), allocatable :: nwgt(:), swgt(:), wwgt(:), ewgt(:)

   logical,  allocatable :: msk(:)

   integer  :: iset, ilat, ilon, iwest, ieast, ie, ipxl
   integer  :: nglist, iloc, ng, ig
   integer  :: iworker, iproc, iio, idest, isrc, nrecv
   integer  :: rmesg(2), smesg(2)
   integer  :: iy, ix, xblk, yblk, xloc, yloc
   integer  :: ipxstt, ipxend

   real(r8) :: lon, lonw, lone, latn, lats
   real(r8) :: distn, dists, distw, diste, diffw, diffe, areathis

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN

         write(*,*)
         write(*,"(A, I0, A, I0, A)") &
            'Building bilinear interpolation from grid to pixel set: ', &
            fgrid%nlat, ' grids in latitude ', fgrid%nlon, ' grids in longitude.'
         write(*,*)

         IF (.not. (lon_between_floor(pixel%edgew, fgrid%lon_w(1), fgrid%lon_e(fgrid%nlon)) &
            .and. lon_between_ceil(pixel%edgee, fgrid%lon_w(1), fgrid%lon_e(fgrid%nlon)))) THEN
            write(*,'(A)') 'Warning: Grid does not cover longitude range of modeling region.'
         ENDIF

         IF (fgrid%yinc == 1) THEN
            IF (.not. ((pixel%edges >= fgrid%lat_s(1)) &
               .and. (pixel%edgen <= fgrid%lat_n(fgrid%nlat)))) THEN
               write(*,'(A)') 'Warning: Grid does not cover latitude range of modeling region.'
            ENDIF
         ELSE
            IF (.not. ((pixel%edges >= fgrid%lat_s(fgrid%nlat)) &
               .and. (pixel%edgen <= fgrid%lat_n(1)))) THEN
               write(*,'(A)') 'Warning: Grid does not cover latitude range of modeling region.'
            ENDIF
         ENDIF

      ENDIF

      this%grid%nlat = fgrid%nlat
      this%grid%nlon = fgrid%nlon

      allocate (this%grid%xblk (size(fgrid%xblk)));  this%grid%xblk = fgrid%xblk
      allocate (this%grid%yblk (size(fgrid%yblk)));  this%grid%yblk = fgrid%yblk
      allocate (this%grid%xloc (size(fgrid%xloc)));  this%grid%xloc = fgrid%xloc
      allocate (this%grid%yloc (size(fgrid%yloc)));  this%grid%yloc = fgrid%yloc
      allocate (this%grid%xcnt (size(fgrid%xcnt)));  this%grid%xcnt = fgrid%xcnt
      allocate (this%grid%ycnt (size(fgrid%ycnt)));  this%grid%ycnt = fgrid%ycnt

      IF (p_is_worker) THEN

         allocate (this%grid%lat_s(this%grid%nlat));  this%grid%lat_s = fgrid%lat_s
         allocate (this%grid%lat_n(this%grid%nlat));  this%grid%lat_n = fgrid%lat_n
         allocate (this%grid%lon_w(this%grid%nlon));  this%grid%lon_w = fgrid%lon_w
         allocate (this%grid%lon_e(this%grid%nlon));  this%grid%lon_e = fgrid%lon_e
         allocate (this%grid%rlon (this%grid%nlon));  CALL this%grid%set_rlon ()
         allocate (this%grid%rlat (this%grid%nlat));  CALL this%grid%set_rlat ()

         this%npset = pixelset%nset

         allocate (yn (this%npset))
         allocate (ys (this%npset))
         allocate (xw (this%npset))
         allocate (xe (this%npset))
         allocate (rlon_pset (this%npset))
         allocate (rlat_pset (this%npset))

         CALL pixelset%get_lonlat_radian (rlon_pset, rlat_pset)

         allocate (xlist(4*this%npset))
         allocate (ylist(4*this%npset))

         allocate (nwgt (this%npset))
         allocate (swgt (this%npset))
         allocate (wwgt (this%npset))
         allocate (ewgt (this%npset))

         nglist = 0

         DO iset = 1, this%npset

            IF (this%grid%rlat(1) > this%grid%rlat(this%grid%nlat)) THEN
               ! from north to south
               ilat = 1
               DO WHILE ((rlat_pset(iset) < this%grid%rlat(ilat)) .and. (ilat < this%grid%nlat))
                  ilat = ilat + 1
               ENDDO

               IF (rlat_pset(iset) >= this%grid%rlat(ilat)) THEN
                  yn(iset) = max(ilat-1,1)
                  ys(iset) = ilat
               ELSE
                  yn(iset) = this%grid%nlat
                  ys(iset) = this%grid%nlat
               ENDIF
            ELSE
               ! from south to north
               ilat = this%grid%nlat
               DO WHILE ((rlat_pset(iset) < this%grid%rlat(ilat)) .and. (ilat > 1))
                  ilat = ilat - 1
               ENDDO

               IF (rlat_pset(iset) >= this%grid%rlat(ilat)) THEN
                  yn(iset) = min(ilat+1,this%grid%nlat)
                  ys(iset) = ilat
               ELSE
                  yn(iset) = 1
                  ys(iset) = 1
               ENDIF
            ENDIF

            IF (yn(iset) /= ys(iset)) THEN
               latn = this%grid%rlat(yn(iset))
               lats = this%grid%rlat(ys(iset))
               distn = arclen(rlat_pset(iset), rlon_pset(iset), latn, rlon_pset(iset))
               dists = arclen(rlat_pset(iset), rlon_pset(iset), lats, rlon_pset(iset))
               nwgt(iset) = dists/(dists+distn)
               swgt(iset) = distn/(dists+distn)
            ELSE
               nwgt(iset) = 1.0
               swgt(iset) = 0.0
            ENDIF


            lon = rlon_pset(iset)*180.0/pi
            CALL normalize_longitude (lon)

            DO iwest = 1, this%grid%nlon
               lonw = this%grid%rlon(iwest) *180.0/pi
               CALL normalize_longitude (lonw)

               ieast = mod(iwest,this%grid%nlon) + 1
               lone  = this%grid%rlon(ieast)*180.0/pi
               CALL normalize_longitude (lone)

               IF (lon_between_floor(lon, lonw, lone)) EXIT
            ENDDO

            xw(iset) = iwest
            xe(iset) = ieast

            ! for the case grid does not cover [-180,180)
            IF ((iwest == this%grid%nlon) .and. (this%grid%nlon > 1)) THEN
               IF (lon_between_floor( &
                  this%grid%lon_e(this%grid%nlon), lonw, this%grid%lon_w(1))) THEN

                  diffw = lon - lonw;  IF (diffw < 0) diffw = diffw + 360.0
                  diffe = lone - lon;  IF (diffe < 0) diffe = diffe + 360.0

                  IF (diffw > diffe) THEN
                     xw(iset) = ieast
                     xe(iset) = ieast
                  ELSE
                     xw(iset) = iwest
                     xe(iset) = iwest
                  ENDIF

               ENDIF
            ENDIF

            IF (xw(iset) /= xe(iset)) THEN
               lonw = this%grid%rlon(xw(iset))
               lone = this%grid%rlon(xe(iset))
               distw = arclen(rlat_pset(iset), rlon_pset(iset), rlat_pset(iset), lonw)
               diste = arclen(rlat_pset(iset), rlon_pset(iset), rlat_pset(iset), lone)
               wwgt(iset) = diste/(distw+diste)
               ewgt(iset) = distw/(distw+diste)
            ELSE
               wwgt(iset) = 1.0
               ewgt(iset) = 0.0
            ENDIF

            CALL insert_into_sorted_list2 ( xw(iset), yn(iset), nglist, xlist, ylist, iloc)
            CALL insert_into_sorted_list2 ( xe(iset), yn(iset), nglist, xlist, ylist, iloc)
            CALL insert_into_sorted_list2 ( xw(iset), ys(iset), nglist, xlist, ylist, iloc)
            CALL insert_into_sorted_list2 ( xe(iset), ys(iset), nglist, xlist, ylist, iloc)

         ENDDO

#ifdef USEMPI
         allocate (ipt (nglist))
         allocate (msk (nglist))
         DO ig = 1, nglist
            xblk = this%grid%xblk(xlist(ig))
            yblk = this%grid%yblk(ylist(ig))
            ipt(ig) = gblock%pio(xblk,yblk)
         ENDDO
#endif

         allocate (this%glist (0:p_np_io-1))
         DO iproc = 0, p_np_io-1
#ifdef USEMPI
            msk = (ipt == p_address_io(iproc))
            ng  = count(msk)
#else
            ng  = nglist
#endif

            this%glist(iproc)%ng = ng

            IF (ng > 0) THEN
               allocate (this%glist(iproc)%ilat (ng))
               allocate (this%glist(iproc)%ilon (ng))

#ifdef USEMPI
               this%glist(iproc)%ilon = pack(xlist(1:nglist), msk)
               this%glist(iproc)%ilat = pack(ylist(1:nglist), msk)
#else
               this%glist(iproc)%ilon = xlist(1:nglist)
               this%glist(iproc)%ilat = ylist(1:nglist)
#endif
            ENDIF

         ENDDO

         deallocate (xlist)
         deallocate (ylist)

#ifdef USEMPI
         deallocate (ipt)
         deallocate (msk)
#endif
      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN
         DO iproc = 0, p_np_io-1
            idest = p_address_io(iproc)
            smesg = (/p_iam_glb, this%glist(iproc)%ng/)

            CALL mpi_send (smesg, 2, MPI_INTEGER, &
               idest, mpi_tag_mesg, p_comm_glb, p_err)

            IF (this%glist(iproc)%ng > 0) THEN
               CALL mpi_send (this%glist(iproc)%ilon, this%glist(iproc)%ng, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
               CALL mpi_send (this%glist(iproc)%ilat, this%glist(iproc)%ng, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
            ENDIF
         ENDDO
      ENDIF

      IF (p_is_io) THEN

         allocate (this%glist (0:p_np_worker-1))

         DO iworker = 0, p_np_worker-1

            CALL mpi_recv (rmesg, 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = rmesg(1)
            nrecv = rmesg(2)
            iproc = p_itis_worker(isrc)

            this%glist(iproc)%ng = nrecv

            IF (nrecv > 0) THEN
               allocate (this%glist(iproc)%ilon (nrecv))
               allocate (this%glist(iproc)%ilat (nrecv))

               CALL mpi_recv (this%glist(iproc)%ilon, nrecv, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL mpi_recv (this%glist(iproc)%ilat, nrecv, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
            ENDIF

         ENDDO

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_worker) THEN

         allocate (this%address (this%npset))
         allocate (this%npart   (this%npset))
         allocate (this%areapart(this%npset))

         DO iset = 1, pixelset%nset

            this%npart(iset) = 4

            allocate (this%address (iset)%val(2,4))
            allocate (this%areapart(iset)%val(4))

            areathis = 0.

            ie = pixelset%ielm(iset)

            ipxstt = pixelset%ipxstt(iset)
            ipxend = pixelset%ipxend(iset)

            ! deal with 2m WMO patch
            IF (ipxstt==-1 .and. ipxend==-1) THEN
               ipxstt = 1
               ipxend = mesh(ie)%npxl
            ENDIF

            DO ipxl = ipxstt, ipxend
               areathis = areathis + areaquad (&
                  pixel%lat_s(mesh(ie)%ilat(ipxl)), pixel%lat_n(mesh(ie)%ilat(ipxl)), &
                  pixel%lon_w(mesh(ie)%ilon(ipxl)), pixel%lon_e(mesh(ie)%ilon(ipxl)) )
            ENDDO

            IF (pixelset%has_shared) THEN
               areathis = areathis * pixelset%pctshared(iset)
            ENDIF

            ! northwest grid
            ix = xw(iset); iy = yn(iset);
#ifdef USEMPI
            xblk = this%grid%xblk(ix)
            yblk = this%grid%yblk(iy)
            iproc = p_itis_io(gblock%pio(xblk,yblk))
#else
            iproc = 0
#endif
            this%address(iset)%val(1,1) = iproc
            this%address(iset)%val(2,1) = find_in_sorted_list2 ( ix, iy, &
               this%glist(iproc)%ng, this%glist(iproc)%ilon, this%glist(iproc)%ilat)

            this%areapart(iset)%val(1) = areathis * nwgt(iset) * wwgt(iset)

            ! northeast grid
            ix = xe(iset); iy = yn(iset);
#ifdef USEMPI
            xblk = this%grid%xblk(ix)
            yblk = this%grid%yblk(iy)
            iproc = p_itis_io(gblock%pio(xblk,yblk))
#else
            iproc = 0
#endif
            this%address(iset)%val(1,2) = iproc
            this%address(iset)%val(2,2) = find_in_sorted_list2 ( ix, iy, &
               this%glist(iproc)%ng, this%glist(iproc)%ilon, this%glist(iproc)%ilat)

            this%areapart(iset)%val(2) = areathis * nwgt(iset) * ewgt(iset)

            ! southwest
            ix = xw(iset); iy = ys(iset);
#ifdef USEMPI
            xblk = this%grid%xblk(ix)
            yblk = this%grid%yblk(iy)
            iproc = p_itis_io(gblock%pio(xblk,yblk))
#else
            iproc = 0
#endif
            this%address(iset)%val(1,3) = iproc
            this%address(iset)%val(2,3) = find_in_sorted_list2 ( ix, iy, &
               this%glist(iproc)%ng, this%glist(iproc)%ilon, this%glist(iproc)%ilat)

            this%areapart(iset)%val(3) = areathis * swgt(iset) * wwgt(iset)

            ! southeast
            ix = xe(iset); iy = ys(iset);
#ifdef USEMPI
            xblk = this%grid%xblk(ix)
            yblk = this%grid%yblk(iy)
            iproc = p_itis_io(gblock%pio(xblk,yblk))
#else
            iproc = 0
#endif
            this%address(iset)%val(1,4) = iproc
            this%address(iset)%val(2,4) = find_in_sorted_list2 ( ix, iy, &
               this%glist(iproc)%ng, this%glist(iproc)%ilon, this%glist(iproc)%ilat)

            this%areapart(iset)%val(4) = areathis * swgt(iset) * ewgt(iset)

         ENDDO

      ENDIF

      IF (p_is_worker) THEN
         IF (this%npset > 0) THEN
            allocate (this%areapset (this%npset))
         ENDIF
         DO iset = 1, this%npset
            this%areapset(iset) = sum(this%areapart(iset)%val)
         ENDDO
      ENDIF

      IF (p_is_io)  CALL allocate_block_data (fgrid, this%areagrid)
      IF (p_is_worker) THEN
         IF (this%npset > 0) THEN
            allocate (msk (this%npset))
            msk = pixelset%ipxstt > 0 .and. pixelset%ipxend > 0
         ENDIF
      ENDIF
      CALL this%get_sumarea (this%areagrid, msk)


      IF (allocated(this%grid%lat_s)) deallocate(this%grid%lat_s)
      IF (allocated(this%grid%lat_n)) deallocate(this%grid%lat_n)
      IF (allocated(this%grid%lon_w)) deallocate(this%grid%lon_w)
      IF (allocated(this%grid%lon_e)) deallocate(this%grid%lon_e)
      IF (allocated(this%grid%rlon )) deallocate(this%grid%rlon )
      IF (allocated(this%grid%rlat )) deallocate(this%grid%rlat )

      IF (allocated(yn)) deallocate(yn)
      IF (allocated(ys)) deallocate(ys)
      IF (allocated(xw)) deallocate(xw)
      IF (allocated(xe)) deallocate(xe)

      IF (allocated(rlon_pset)) deallocate(rlon_pset)
      IF (allocated(rlat_pset)) deallocate(rlat_pset)

      IF (allocated(nwgt)) deallocate(nwgt)
      IF (allocated(swgt)) deallocate(swgt)
      IF (allocated(wwgt)) deallocate(wwgt)
      IF (allocated(ewgt)) deallocate(ewgt)

      IF (allocated(msk)) deallocate(msk)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE spatial_mapping_build_bilinear

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_set_missing_value (this, gdata, missing_value, pmask)

   USE MOD_Precision
   USE MOD_Block
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_SPMD_Task
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   type(block_data_real8_2d), intent(in) :: gdata
   real(r8), intent(in) :: missing_value

   logical,  intent(inout), optional :: pmask(:)

   ! Local variables
   integer :: iproc, idest, isrc
   integer :: ig, ilon, ilat, xblk, yblk, xloc, yloc, iloc, iset, ipart, iblkme

   real(r8), allocatable :: gbuff(:)
   type(pointer_real8_1d), allocatable :: pbuff(:)

      this%has_missing_value = .true.
      this%missing_value = missing_value

      IF (p_is_io) THEN

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (this%glist(iproc)%ng))

               DO ig = 1, this%glist(iproc)%ng
                  ilon = this%glist(iproc)%ilon(ig)
                  ilat = this%glist(iproc)%ilat(ig)
                  xblk = this%grid%xblk (ilon)
                  yblk = this%grid%yblk (ilat)
                  xloc = this%grid%xloc (ilon)
                  yloc = this%grid%yloc (ilat)

                  gbuff(ig) = gdata%blk(xblk,yblk)%val(xloc,yloc)

               ENDDO

#ifdef USEMPI
               idest = p_address_worker(iproc)
               CALL mpi_send (gbuff, this%glist(iproc)%ng, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err)

               deallocate (gbuff)
#endif
            ENDIF
         ENDDO

         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            WHERE (gdata%blk(xblk,yblk)%val == missing_value)
               this%areagrid%blk(xblk,yblk)%val = 0.
            ENDWHERE
         ENDDO

      ENDIF

      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (pbuff(iproc)%val (this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_io(iproc)
               CALL mpi_recv (pbuff(iproc)%val, this%glist(iproc)%ng, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               pbuff(0)%val = gbuff
               deallocate (gbuff)
#endif
            ENDIF
         ENDDO


         DO iset = 1, this%npset

            this%areapset(iset) = 0.

            DO ipart = 1, this%npart(iset)
               iproc = this%address(iset)%val(1,ipart)
               iloc  = this%address(iset)%val(2,ipart)

               IF (pbuff(iproc)%val(iloc) == missing_value) THEN
                  this%areapart(iset)%val(ipart) = 0.
               ELSE
                  this%areapset(iset) = this%areapset(iset) + this%areapart(iset)%val(ipart)
               ENDIF
            ENDDO

            IF (present(pmask)) THEN
               pmask(iset) = (this%areapset(iset) > 0.)
            ENDIF

         ENDDO

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO
         deallocate (pbuff)

      ENDIF

   END SUBROUTINE spatial_mapping_set_missing_value

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_pset2grid_2d (this, pdata, gdata, spv, msk, input_mode)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_SPMD_Task
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   real(r8), intent(in) :: pdata(:)
   type(block_data_real8_2d), intent(inout) :: gdata

   real(r8), intent(in), optional :: spv
   logical,  intent(in), optional :: msk(:)

   character(len=*), intent(in), optional :: input_mode

   ! Local variables
   integer :: iproc, idest, isrc
   integer :: ig, ilon, ilat, xblk, yblk, xloc, yloc, iloc, iset, ipart

   real(r8), allocatable :: gbuff(:)
   type(pointer_real8_1d), allocatable :: pbuff(:)
   character(len=256) :: inmode
   real(r8) :: sumwt

      IF (p_is_worker) THEN

         inmode = 'average'
         IF (present(input_mode)) inmode = trim(input_mode)

         allocate (pbuff (0:p_np_io-1))

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               allocate (pbuff(iproc)%val (this%glist(iproc)%ng))

               IF (present(spv)) THEN
                  pbuff(iproc)%val(:) = spv
               ELSE
                  pbuff(iproc)%val(:) = 0.0
               ENDIF
            ENDIF
         ENDDO

         DO iset = 1, this%npset

            IF (present(spv)) THEN
               IF (pdata(iset) == spv) CYCLE
            ENDIF

            IF (present(msk)) THEN
               IF (.not. msk(iset)) CYCLE
            ENDIF

            IF ((this%npart(iset) > 0) .and. (trim(inmode) == 'total')) THEN
               sumwt = sum(this%areapart(iset)%val)
            ELSE
               sumwt = 1.
            ENDIF

            DO ipart = 1, this%npart(iset)
               iproc = this%address(iset)%val(1,ipart)
               iloc  = this%address(iset)%val(2,ipart)

               IF (present(spv)) THEN
                  IF (pbuff(iproc)%val(iloc) /= spv) THEN
                     pbuff(iproc)%val(iloc) = pbuff(iproc)%val(iloc) &
                        + pdata(iset)/sumwt * this%areapart(iset)%val(ipart)
                  ELSE
                     pbuff(iproc)%val(iloc) = &
                        pdata(iset)/sumwt * this%areapart(iset)%val(ipart)
                  ENDIF
               ELSE
                  pbuff(iproc)%val(iloc) = pbuff(iproc)%val(iloc) &
                     + pdata(iset)/sumwt * this%areapart(iset)%val(ipart)
               ENDIF
            ENDDO
         ENDDO

#ifdef USEMPI
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               idest = p_address_io(iproc)
               CALL mpi_send (pbuff(iproc)%val, this%glist(iproc)%ng, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
            ENDIF
         ENDDO
#endif

      ENDIF

      IF (p_is_io) THEN

         IF (present(spv)) THEN
            CALL flush_block_data (gdata, spv)
         ELSE
            CALL flush_block_data (gdata, 0.0_r8)
         ENDIF

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_worker(iproc)
               CALL mpi_recv (gbuff, this%glist(iproc)%ng, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               gbuff = pbuff(0)%val
#endif

               DO ig = 1, this%glist(iproc)%ng
                  IF (present(spv)) THEN
                     IF (gbuff(ig) /= spv) THEN
                        ilon = this%glist(iproc)%ilon(ig)
                        ilat = this%glist(iproc)%ilat(ig)
                        xblk = this%grid%xblk (ilon)
                        yblk = this%grid%yblk (ilat)
                        xloc = this%grid%xloc (ilon)
                        yloc = this%grid%yloc (ilat)

                        IF (gdata%blk(xblk,yblk)%val(xloc,yloc) /= spv) THEN
                           gdata%blk(xblk,yblk)%val(xloc,yloc) = &
                              gdata%blk(xblk,yblk)%val(xloc,yloc) + gbuff(ig)
                        ELSE
                           gdata%blk(xblk,yblk)%val(xloc,yloc) = gbuff(ig)
                        ENDIF
                     ENDIF
                  ELSE
                     ilon = this%glist(iproc)%ilon(ig)
                     ilat = this%glist(iproc)%ilat(ig)
                     xblk = this%grid%xblk (ilon)
                     yblk = this%grid%yblk (ilat)
                     xloc = this%grid%xloc (ilon)
                     yloc = this%grid%yloc (ilat)

                     gdata%blk(xblk,yblk)%val(xloc,yloc) = &
                        gdata%blk(xblk,yblk)%val(xloc,yloc) + gbuff(ig)
                  ENDIF
               ENDDO

               deallocate (gbuff)
            ENDIF
         ENDDO

      ENDIF

      IF (p_is_worker) THEN
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO
         deallocate (pbuff)
      ENDIF


   END SUBROUTINE spatial_mapping_pset2grid_2d

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_pset2grid_3d (this, pdata, gdata, spv, msk)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_SPMD_Task
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   real(r8), intent(in) :: pdata(:,:)
   type(block_data_real8_3d), intent(inout) :: gdata

   real(r8), intent(in), optional :: spv
   logical,  intent(in), optional :: msk(:)

   ! Local variables
   integer :: iproc, idest, isrc
   integer :: ig, ilon, ilat, iloc, iset, ipart
   integer :: xblk, yblk, xloc, yloc
   integer :: lb1, ub1, i1

   real(r8), allocatable :: gbuff(:,:)
   type(pointer_real8_2d), allocatable :: pbuff(:)


      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))

         lb1 = lbound(pdata,1)
         ub1 = ubound(pdata,1)

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               allocate (pbuff(iproc)%val (lb1:ub1, this%glist(iproc)%ng))

               IF (present(spv)) THEN
                  pbuff(iproc)%val(:,:) = spv
               ELSE
                  pbuff(iproc)%val(:,:) = 0.0
               ENDIF
            ENDIF
         ENDDO

         DO iset = 1, this%npset
            IF (present(msk)) THEN
               IF (.not. msk(iset)) CYCLE
            ENDIF

            DO ipart = 1, this%npart(iset)
               iproc = this%address(iset)%val(1,ipart)
               iloc  = this%address(iset)%val(2,ipart)

               DO i1 = lb1, ub1
                  IF (present(spv)) THEN
                     IF (pdata(i1,iset) /= spv) THEN
                        IF (pbuff(iproc)%val(i1,iloc) /= spv) THEN
                           pbuff(iproc)%val(i1,iloc) = pbuff(iproc)%val(i1,iloc) &
                              + pdata(i1,iset) * this%areapart(iset)%val(ipart)
                        ELSE
                           pbuff(iproc)%val(i1,iloc) = &
                              pdata(i1,iset) * this%areapart(iset)%val(ipart)
                        ENDIF
                     ENDIF
                  ELSE
                     pbuff(iproc)%val(i1,iloc) = pbuff(iproc)%val(i1,iloc) &
                        + pdata(i1,iset) * this%areapart(iset)%val(ipart)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

#ifdef USEMPI
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               idest = p_address_io(iproc)
               CALL mpi_send (pbuff(iproc)%val, &
                  (ub1-lb1+1) * this%glist(iproc)%ng, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err)

            ENDIF
         ENDDO
#endif

      ENDIF

      IF (p_is_io) THEN

         lb1 = gdata%lb1
         ub1 = gdata%ub1

         IF (present(spv)) THEN
            CALL flush_block_data (gdata, spv)
         ELSE
            CALL flush_block_data (gdata, 0.0_r8)
         ENDIF

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (lb1:ub1, this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_worker(iproc)
               CALL mpi_recv (gbuff, &
                  (ub1-lb1+1) * this%glist(iproc)%ng, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               gbuff = pbuff(0)%val
#endif

               DO ig = 1, this%glist(iproc)%ng
                  ilon = this%glist(iproc)%ilon(ig)
                  ilat = this%glist(iproc)%ilat(ig)
                  xblk = this%grid%xblk (ilon)
                  yblk = this%grid%yblk (ilat)
                  xloc = this%grid%xloc (ilon)
                  yloc = this%grid%yloc (ilat)

                  DO i1 = lb1, ub1
                     IF (present(spv)) THEN
                        IF (gbuff(i1,ig) /= spv) THEN
                           IF (gdata%blk(xblk,yblk)%val(i1,xloc,yloc) /= spv) THEN
                              gdata%blk(xblk,yblk)%val(i1,xloc,yloc) = &
                                 gdata%blk(xblk,yblk)%val(i1,xloc,yloc) + gbuff(i1,ig)
                           ELSE
                              gdata%blk(xblk,yblk)%val(i1,xloc,yloc) = gbuff(i1,ig)
                           ENDIF
                        ENDIF
                     ELSE
                        gdata%blk(xblk,yblk)%val(i1,xloc,yloc) = &
                           gdata%blk(xblk,yblk)%val(i1,xloc,yloc) + gbuff(i1,ig)
                     ENDIF
                  ENDDO
               ENDDO

               deallocate (gbuff)
            ENDIF

         ENDDO

      ENDIF

      IF (p_is_worker) THEN
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO
         deallocate (pbuff)
      ENDIF

   END SUBROUTINE spatial_mapping_pset2grid_3d

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_pset2grid_4d (this, pdata, gdata, spv, msk)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_SPMD_Task
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   real(r8), intent(in) :: pdata(:,:,:)
   type(block_data_real8_4d), intent(inout) :: gdata

   real(r8), intent(in), optional :: spv
   logical,  intent(in), optional :: msk(:)

   ! Local variables
   integer :: iproc, idest, isrc
   integer :: ig, ilon, ilat, iloc, iset, ipart
   integer :: xblk, yblk, xloc, yloc
   integer :: lb1, ub1, i1, ndim1, lb2, ub2, i2, ndim2

   real(r8), allocatable :: gbuff(:,:,:)
   type(pointer_real8_3d), allocatable :: pbuff(:)

      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))

         lb1 = lbound(pdata,1)
         ub1 = ubound(pdata,1)
         ndim1 = ub1 - lb1 + 1

         lb2 = lbound(pdata,2)
         ub2 = ubound(pdata,2)
         ndim2 = ub2 - lb2 + 1

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               allocate (pbuff(iproc)%val (lb1:ub1, lb2:ub2, this%glist(iproc)%ng))

               IF (present(spv)) THEN
                  pbuff(iproc)%val(:,:,:) = spv
               ELSE
                  pbuff(iproc)%val(:,:,:) = 0.0
               ENDIF
            ENDIF
         ENDDO

         DO iset = 1, this%npset
            IF (present(msk)) THEN
               IF (.not. msk(iset)) CYCLE
            ENDIF

            DO ipart = 1, this%npart(iset)
               iproc = this%address(iset)%val(1,ipart)
               iloc  = this%address(iset)%val(2,ipart)

               DO i1 = lb1, ub1
                  DO i2 = lb2, ub2
                     IF (present(spv)) THEN
                        IF (pdata(i1,i2,iset) /= spv) THEN
                           IF (pbuff(iproc)%val(i1,i2,iloc) /= spv) THEN
                              pbuff(iproc)%val(i1,i2,iloc) = pbuff(iproc)%val(i1,i2,iloc) &
                                 + pdata(i1,i2,iset) * this%areapart(iset)%val(ipart)
                           ELSE
                              pbuff(iproc)%val(i1,i2,iloc) = &
                                 pdata(i1,i2,iset) * this%areapart(iset)%val(ipart)
                           ENDIF
                        ENDIF
                     ELSE
                        pbuff(iproc)%val(i1,i2,iloc) = pbuff(iproc)%val(i1,i2,iloc) &
                           + pdata(i1,i2,iset) * this%areapart(iset)%val(ipart)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

#ifdef USEMPI
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               idest = p_address_io(iproc)
               CALL mpi_send (pbuff(iproc)%val, ndim1 * ndim2 * this%glist(iproc)%ng, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
            ENDIF
         ENDDO
#endif

      ENDIF

      IF (p_is_io) THEN

         lb1 = gdata%lb1
         ub1 = gdata%ub1
         ndim1 = ub1 - lb1 + 1

         lb2 = gdata%lb2
         ub2 = gdata%ub2
         ndim2 = ub2 - lb2 + 1

         IF (present(spv)) THEN
            CALL flush_block_data (gdata, spv)
         ELSE
            CALL flush_block_data (gdata, 0.0_r8)
         ENDIF

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (lb1:ub1, lb2:ub2, this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_worker(iproc)
               CALL mpi_recv (gbuff, ndim1 * ndim2 * this%glist(iproc)%ng, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               gbuff = pbuff(0)%val
#endif

               DO ig = 1, this%glist(iproc)%ng
                  ilon = this%glist(iproc)%ilon(ig)
                  ilat = this%glist(iproc)%ilat(ig)
                  xblk = this%grid%xblk (ilon)
                  yblk = this%grid%yblk (ilat)
                  xloc = this%grid%xloc (ilon)
                  yloc = this%grid%yloc (ilat)

                  DO i1 = lb1, ub1
                     DO i2 = lb2, ub2
                        IF (present(spv)) THEN
                           IF (gbuff(i1,i2,ig) /= spv) THEN
                              IF (gdata%blk(xblk,yblk)%val(i1,i2,xloc,yloc) /= spv) THEN
                                 gdata%blk(xblk,yblk)%val(i1,i2,xloc,yloc) = &
                                    gdata%blk(xblk,yblk)%val(i1,i2,xloc,yloc) + gbuff(i1,i2,ig)
                              ELSE
                                 gdata%blk(xblk,yblk)%val(i1,i2,xloc,yloc) = gbuff(i1,i2,ig)
                              ENDIF
                           ENDIF
                        ELSE
                           gdata%blk(xblk,yblk)%val(i1,i2,xloc,yloc) = &
                              gdata%blk(xblk,yblk)%val(i1,i2,xloc,yloc) + gbuff(i1,i2,ig)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO

               deallocate (gbuff)
            ENDIF
         ENDDO
      ENDIF

      IF (p_is_worker) THEN
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO
         deallocate (pbuff)
      ENDIF

   END SUBROUTINE spatial_mapping_pset2grid_4d

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_pset2grid_max (this, pdata, gdata, spv, msk)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   real(r8), intent(in) :: pdata(:)
   type(block_data_real8_2d), intent(inout) :: gdata

   real(r8), intent(in), optional :: spv
   logical,  intent(in), optional :: msk(:)

   ! Local variables
   integer  :: iproc, idest, isrc
   integer  :: ig, ilon, ilat, xblk, yblk, xloc, yloc, iloc, iset, ipart

   real(r8), allocatable :: gbuff(:)
   type(pointer_real8_1d), allocatable :: pbuff(:)

      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               allocate (pbuff(iproc)%val (this%glist(iproc)%ng))
               pbuff(iproc)%val(:) = spval
            ENDIF
         ENDDO

         DO iset = 1, this%npset

            IF (present(spv)) THEN
               IF (pdata(iset) == spv) CYCLE
            ENDIF

            IF (present(msk)) THEN
               IF (.not. msk(iset)) CYCLE
            ENDIF

            DO ipart = 1, this%npart(iset)
               iproc = this%address(iset)%val(1,ipart)
               iloc  = this%address(iset)%val(2,ipart)

               IF (pbuff(iproc)%val(iloc) /= spval) THEN
                  pbuff(iproc)%val(iloc) = max(pdata(iset), pbuff(iproc)%val(iloc))
               ELSE
                  pbuff(iproc)%val(iloc) = pdata(iset)
               ENDIF
            ENDDO
         ENDDO

#ifdef USEMPI
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               idest = p_address_io(iproc)
               CALL mpi_send (pbuff(iproc)%val, this%glist(iproc)%ng, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
            ENDIF
         ENDDO
#endif

      ENDIF

      IF (p_is_io) THEN

         CALL flush_block_data (gdata, spval)

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_worker(iproc)
               CALL mpi_recv (gbuff, this%glist(iproc)%ng, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               gbuff = pbuff(0)%val
#endif

               DO ig = 1, this%glist(iproc)%ng
                  IF (gbuff(ig) /= spval) THEN
                     ilon = this%glist(iproc)%ilon(ig)
                     ilat = this%glist(iproc)%ilat(ig)
                     xblk = this%grid%xblk (ilon)
                     yblk = this%grid%yblk (ilat)
                     xloc = this%grid%xloc (ilon)
                     yloc = this%grid%yloc (ilat)

                     IF (gdata%blk(xblk,yblk)%val(xloc,yloc) /= spval) THEN
                        gdata%blk(xblk,yblk)%val(xloc,yloc) = &
                           max(gdata%blk(xblk,yblk)%val(xloc,yloc), gbuff(ig))
                     ELSE
                        gdata%blk(xblk,yblk)%val(xloc,yloc) = gbuff(ig)
                     ENDIF
                  ENDIF
               ENDDO

               deallocate (gbuff)
            ENDIF
         ENDDO

      ENDIF

      IF (p_is_worker) THEN
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO
         deallocate (pbuff)
      ENDIF

   END SUBROUTINE spatial_mapping_pset2grid_max

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_pset2grid_split (this, pdata, settyp, typidx, gdata, spv)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_SPMD_Task
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   real(r8), intent(in) :: pdata (:)
   integer , intent(in) :: settyp(:)
   integer , intent(in) :: typidx(:)
   type(block_data_real8_3d), intent(inout) :: gdata

   real(r8), intent(in) :: spv

   ! Local variables
   integer :: iproc, idest, isrc
   integer :: ig, ilon, ilat, iloc, iset, ipart, ityp, ntyps
   integer :: xblk, yblk, xloc, yloc

   real(r8), allocatable :: gbuff(:)
   type(pointer_real8_1d), allocatable :: pbuff (:)

      IF (p_is_worker) THEN
         allocate (pbuff (0:p_np_io-1))
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               allocate (pbuff(iproc)%val  (this%glist(iproc)%ng))
            ENDIF
         ENDDO
      ENDIF

      IF (p_is_io) THEN
         CALL flush_block_data (gdata, spv)
      ENDIF

      ntyps = size(typidx)

      DO ityp = 1, ntyps

         IF (p_is_worker) THEN

            DO iproc = 0, p_np_io-1
               IF (this%glist(iproc)%ng > 0) THEN
                  pbuff(iproc)%val(:) = spv
               ENDIF
            ENDDO

            DO iset = 1, this%npset
               IF ((settyp(iset) == typidx(ityp)) .and. (pdata(iset) /= spv)) THEN
                  DO ipart = 1, this%npart(iset)
                     iproc = this%address(iset)%val(1,ipart)
                     iloc  = this%address(iset)%val(2,ipart)

                     IF (pbuff(iproc)%val(iloc) /= spv) THEN
                        pbuff(iproc)%val(iloc) = pbuff(iproc)%val(iloc) &
                           + pdata(iset) * this%areapart(iset)%val(ipart)
                     ELSE
                        pbuff(iproc)%val(iloc) = &
                           pdata(iset) * this%areapart(iset)%val(ipart)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO

#ifdef USEMPI
            DO iproc = 0, p_np_io-1
               IF (this%glist(iproc)%ng > 0) THEN
                  idest = p_address_io(iproc)
                  CALL mpi_send (pbuff(iproc)%val, this%glist(iproc)%ng, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
               ENDIF
            ENDDO
#endif

         ENDIF

         IF (p_is_io) THEN

            DO iproc = 0, p_np_worker-1
               IF (this%glist(iproc)%ng > 0) THEN

                  allocate (gbuff (this%glist(iproc)%ng))

#ifdef USEMPI
                  isrc = p_address_worker(iproc)
                  CALL mpi_recv (gbuff, this%glist(iproc)%ng, MPI_REAL8, &
                     isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
                  gbuff = pbuff(0)%val
#endif

                  DO ig = 1, this%glist(iproc)%ng
                     IF (gbuff(ig) /= spv) THEN
                        ilon = this%glist(iproc)%ilon(ig)
                        ilat = this%glist(iproc)%ilat(ig)
                        xblk = this%grid%xblk (ilon)
                        yblk = this%grid%yblk (ilat)
                        xloc = this%grid%xloc (ilon)
                        yloc = this%grid%yloc (ilat)

                        IF (gdata%blk(xblk,yblk)%val(ityp,xloc,yloc) /= spv) THEN
                           gdata%blk(xblk,yblk)%val(ityp,xloc,yloc) = &
                              gdata%blk(xblk,yblk)%val(ityp,xloc,yloc) + gbuff(ig)
                        ELSE
                           gdata%blk(xblk,yblk)%val(ityp,xloc,yloc) = gbuff(ig)
                        ENDIF
                     ENDIF
                  ENDDO

                  deallocate (gbuff)
               ENDIF

            ENDDO

         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif
      ENDDO

      IF (p_is_worker) THEN
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO
         deallocate (pbuff)
      ENDIF

   END SUBROUTINE spatial_mapping_pset2grid_split

   ! ------------------------------
   SUBROUTINE spatial_mapping_get_sumarea (this, sumarea, filter)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_SPMD_Task
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   type(block_data_real8_2d), intent(inout) :: sumarea
   logical, intent(in), optional :: filter(:)

   ! Local variables
   integer :: iproc, idest, isrc
   integer :: ig, ilon, ilat, xblk, yblk, xloc, yloc, iloc, iset, ipart

   real(r8), allocatable :: gbuff(:)
   type(pointer_real8_1d), allocatable :: pbuff(:)

      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               allocate (pbuff(iproc)%val (this%glist(iproc)%ng))
               pbuff(iproc)%val(:) = 0.0
            ENDIF
         ENDDO

         DO iset = 1, this%npset

            IF (present(filter)) THEN
               IF (.not. filter(iset)) CYCLE
            ENDIF

            DO ipart = 1, this%npart(iset)
               iproc = this%address(iset)%val(1,ipart)
               iloc  = this%address(iset)%val(2,ipart)
               pbuff(iproc)%val(iloc) = pbuff(iproc)%val(iloc) + this%areapart(iset)%val(ipart)
            ENDDO
         ENDDO

#ifdef USEMPI
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               idest = p_address_io(iproc)
               CALL mpi_send (pbuff(iproc)%val, this%glist(iproc)%ng, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
            ENDIF
         ENDDO
#endif

      ENDIF

      IF (p_is_io) THEN

         CALL flush_block_data (sumarea, 0.0_r8)

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_worker(iproc)
               CALL mpi_recv (gbuff, this%glist(iproc)%ng, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               gbuff = pbuff(0)%val
#endif

               DO ig = 1, this%glist(iproc)%ng
                  ilon = this%glist(iproc)%ilon(ig)
                  ilat = this%glist(iproc)%ilat(ig)
                  xblk = this%grid%xblk (ilon)
                  yblk = this%grid%yblk (ilat)
                  xloc = this%grid%xloc (ilon)
                  yloc = this%grid%yloc (ilat)

                  sumarea%blk(xblk,yblk)%val(xloc,yloc) = &
                     sumarea%blk(xblk,yblk)%val(xloc,yloc) + gbuff(ig)
               ENDDO

               deallocate (gbuff)
            ENDIF
         ENDDO

      ENDIF

      IF (p_is_worker) THEN
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO
         deallocate (pbuff)
      ENDIF

   END SUBROUTINE spatial_mapping_get_sumarea

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_grid2pset_2d (this, gdata, pdata)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_Pixelset
   USE MOD_DataType
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   type(block_data_real8_2d), intent(in) :: gdata
   real(r8), intent(out) :: pdata(:)

   ! Local variables
   integer :: iproc, idest, isrc
   integer :: ig, ilon, ilat, xblk, yblk, xloc, yloc, iloc, iset, ipart

   real(r8), allocatable :: gbuff(:)
   type(pointer_real8_1d), allocatable :: pbuff(:)

      IF (p_is_io) THEN

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (this%glist(iproc)%ng))

               DO ig = 1, this%glist(iproc)%ng
                  ilon = this%glist(iproc)%ilon(ig)
                  ilat = this%glist(iproc)%ilat(ig)
                  xblk = this%grid%xblk (ilon)
                  yblk = this%grid%yblk (ilat)
                  xloc = this%grid%xloc (ilon)
                  yloc = this%grid%yloc (ilat)

                  gbuff(ig) = gdata%blk(xblk,yblk)%val(xloc,yloc)

               ENDDO

#ifdef USEMPI
               idest = p_address_worker(iproc)
               CALL mpi_send (gbuff, this%glist(iproc)%ng, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err)

               deallocate (gbuff)
#endif
            ENDIF
         ENDDO

      ENDIF

      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (pbuff(iproc)%val (this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_io(iproc)
               CALL mpi_recv (pbuff(iproc)%val, this%glist(iproc)%ng, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               pbuff(0)%val = gbuff
               deallocate (gbuff)
#endif
            ENDIF
         ENDDO

         DO iset = 1, this%npset

            IF (this%areapset(iset) > 0.) THEN

               pdata(iset) = 0.

               DO ipart = 1, this%npart(iset)
                  iproc = this%address(iset)%val(1,ipart)
                  iloc  = this%address(iset)%val(2,ipart)

                  IF (this%areapart(iset)%val(ipart) > 0) THEN
                     pdata(iset) = pdata(iset) &
                        + pbuff(iproc)%val(iloc) * this%areapart(iset)%val(ipart)
                  ENDIF
               ENDDO

               pdata(iset) = pdata(iset) / this%areapset(iset)

            ELSE
               pdata(iset) = spval
            ENDIF

         ENDDO

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO
         deallocate (pbuff)

      ENDIF

   END SUBROUTINE spatial_mapping_grid2pset_2d

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_grid2pset_3d (this, gdata, ndim1, pdata)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_Pixelset
   USE MOD_DataType
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   type(block_data_real8_3d), intent(in) :: gdata
   integer, intent(in) :: ndim1
   real(r8), intent(out) :: pdata(:,:)

   ! Local variables
   integer :: iproc, idest, isrc
   integer :: ig, ilon, ilat, xblk, yblk, xloc, yloc, iloc, iset, ipart, i

   real(r8), allocatable :: gbuff(:,:)
   type(pointer_real8_2d), allocatable :: pbuff(:)


      IF (p_is_io) THEN

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (ndim1, this%glist(iproc)%ng))

               DO ig = 1, this%glist(iproc)%ng
                  ilon = this%glist(iproc)%ilon(ig)
                  ilat = this%glist(iproc)%ilat(ig)
                  xblk = this%grid%xblk (ilon)
                  yblk = this%grid%yblk (ilat)
                  xloc = this%grid%xloc (ilon)
                  yloc = this%grid%yloc (ilat)

                  gbuff(:,ig) = gdata%blk(xblk,yblk)%val(:,xloc,yloc)
               ENDDO

#ifdef USEMPI
               idest = p_address_worker(iproc)
               CALL mpi_send (gbuff, ndim1 * this%glist(iproc)%ng, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err)

               deallocate (gbuff)
#endif
            ENDIF
         ENDDO

      ENDIF

      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (pbuff(iproc)%val (ndim1, this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_io(iproc)
               CALL mpi_recv (pbuff(iproc)%val, ndim1 * this%glist(iproc)%ng, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               pbuff(0)%val = gbuff
               deallocate (gbuff)
#endif
            ENDIF
         ENDDO


         DO iset = 1, this%npset

            IF (this%areapset(iset) > 0.) THEN

               pdata(:,iset) = 0.

               DO ipart = 1, this%npart(iset)
                  iproc = this%address(iset)%val(1,ipart)
                  iloc  = this%address(iset)%val(2,ipart)

                  IF (this%areapart(iset)%val(ipart) > 0) THEN
                     pdata(:,iset) = pdata(:,iset) &
                        + pbuff(iproc)%val(:,iloc) * this%areapart(iset)%val(ipart)
                  ENDIF
               ENDDO

               pdata(:,iset) = pdata(:,iset) / this%areapset(iset)

            ELSE
               pdata(:,iset) = spval
            ENDIF

         ENDDO

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO
         deallocate (pbuff)

      ENDIF

   END SUBROUTINE spatial_mapping_grid2pset_3d

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_dominant_2d (this, gdata, pdata)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_Pixelset
   USE MOD_DataType
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   type(block_data_int32_2d), intent(in) :: gdata
   integer, intent(out) :: pdata(:)

   ! Local variables
   integer :: iproc, idest, isrc
   integer :: ig, ilon, ilat, xblk, yblk, xloc, yloc, iloc, iset, ipart

   integer, allocatable :: gbuff(:)
   type(pointer_int32_1d), allocatable :: pbuff(:)

      IF (p_is_io) THEN

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (this%glist(iproc)%ng))

               DO ig = 1, this%glist(iproc)%ng
                  ilon = this%glist(iproc)%ilon(ig)
                  ilat = this%glist(iproc)%ilat(ig)
                  xblk = this%grid%xblk (ilon)
                  yblk = this%grid%yblk (ilat)
                  xloc = this%grid%xloc (ilon)
                  yloc = this%grid%yloc (ilat)

                  gbuff(ig) = gdata%blk(xblk,yblk)%val(xloc,yloc)

               ENDDO

#ifdef USEMPI
               idest = p_address_worker(iproc)
               CALL mpi_send (gbuff, this%glist(iproc)%ng, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err)

               deallocate (gbuff)
#endif
            ENDIF
         ENDDO

      ENDIF

      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (pbuff(iproc)%val (this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_io(iproc)
               CALL mpi_recv (pbuff(iproc)%val, this%glist(iproc)%ng, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               pbuff(0)%val = gbuff
               deallocate (gbuff)
#endif
            ENDIF
         ENDDO

         DO iset = 1, this%npset
            IF (this%areapset(iset) > 0.) THEN
               ipart = maxloc(this%areapart(iset)%val, dim=1)
               iproc = this%address(iset)%val(1,ipart)
               iloc  = this%address(iset)%val(2,ipart)
               pdata(iset) = pbuff(iproc)%val(iloc)
            ELSE
               pdata(iset) = -9999
            ENDIF
         ENDDO

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO

         deallocate (pbuff)

      ENDIF

   END SUBROUTINE spatial_mapping_dominant_2d

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_varvalue_2d (this, gdata, pdata)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_Pixelset
   USE MOD_DataType
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   type(block_data_real8_2d), intent(in) :: gdata
   real(r8), intent(inout) :: pdata(:)

   ! Local variables
   integer :: iproc, idest, isrc
   integer :: ig, ilon, ilat, xblk, yblk, xloc, yloc, iloc, iset, ipart

   real(r8), allocatable :: gbuff(:)
   type(pointer_real8_1d), allocatable :: pbuff(:)
   real(r8), allocatable :: pdata_tem(:)

      IF (p_is_io) THEN

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (this%glist(iproc)%ng))

               DO ig = 1, this%glist(iproc)%ng
                  ilon = this%glist(iproc)%ilon(ig)
                  ilat = this%glist(iproc)%ilat(ig)
                  xblk = this%grid%xblk (ilon)
                  yblk = this%grid%yblk (ilat)
                  xloc = this%grid%xloc (ilon)
                  yloc = this%grid%yloc (ilat)

                  gbuff(ig) = gdata%blk(xblk,yblk)%val(xloc,yloc)

               ENDDO

#ifdef USEMPI
               idest = p_address_worker(iproc)
               CALL mpi_send (gbuff, this%glist(iproc)%ng, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err)

               deallocate (gbuff)
#endif
            ENDIF
         ENDDO

      ENDIF

      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))
         allocate (pdata_tem (size(pdata)))

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (pbuff(iproc)%val (this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_io(iproc)
               CALL mpi_recv (pbuff(iproc)%val, this%glist(iproc)%ng, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               pbuff(0)%val = gbuff
               deallocate (gbuff)
#endif
            ENDIF
         ENDDO

         DO iset = 1, this%npset

            IF (this%areapset(iset) > 0.) THEN

               pdata_tem(iset) = 0._r8

               DO ipart = 1, this%npart(iset)
                  iproc = this%address(iset)%val(1,ipart)
                  iloc  = this%address(iset)%val(2,ipart)

                  pdata_tem(iset) = pdata_tem(iset) &
                     + pdata(iset) * pbuff(iproc)%val(iloc) * this%areapart(iset)%val(ipart)
               ENDDO

               pdata_tem(iset) = pdata_tem(iset) / this%areapset(iset)

            ELSE
               pdata_tem(iset) = 0._r8
            ENDIF

         ENDDO

         pdata = pdata_tem

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO
         deallocate (pbuff)
         deallocate (pdata_tem)

      ENDIF

   END SUBROUTINE spatial_mapping_varvalue_2d

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_grid2part (this, gdata, sdata)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_Pixelset
   USE MOD_DataType
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   type(block_data_real8_2d), intent(in)    :: gdata
   type(pointer_real8_1d),    intent(inout) :: sdata(:)

   ! Local variables
   integer :: iproc, idest, isrc
   integer :: ig, ilon, ilat, xblk, yblk, xloc, yloc, iloc, iset, ipart

   real(r8), allocatable :: gbuff(:)
   type(pointer_real8_1d), allocatable :: pbuff(:)

      IF (p_is_io) THEN

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (this%glist(iproc)%ng))

               DO ig = 1, this%glist(iproc)%ng
                  ilon = this%glist(iproc)%ilon(ig)
                  ilat = this%glist(iproc)%ilat(ig)
                  xblk = this%grid%xblk (ilon)
                  yblk = this%grid%yblk (ilat)
                  xloc = this%grid%xloc (ilon)
                  yloc = this%grid%yloc (ilat)

                  gbuff(ig) = gdata%blk(xblk,yblk)%val(xloc,yloc)

               ENDDO

#ifdef USEMPI
               idest = p_address_worker(iproc)
               CALL mpi_send (gbuff, this%glist(iproc)%ng, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err)

               deallocate (gbuff)
#endif
            ENDIF
         ENDDO

      ENDIF

      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (pbuff(iproc)%val (this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_io(iproc)
               CALL mpi_recv (pbuff(iproc)%val, this%glist(iproc)%ng, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               pbuff(0)%val = gbuff
               deallocate (gbuff)
#endif
            ENDIF
         ENDDO

         DO iset = 1, this%npset
            DO ipart = 1, this%npart(iset)
               iproc = this%address(iset)%val(1,ipart)
               iloc  = this%address(iset)%val(2,ipart)

               sdata(iset)%val(ipart) = pbuff(iproc)%val(iloc)
            ENDDO
         ENDDO

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO
         deallocate (pbuff)

      ENDIF

   END SUBROUTINE spatial_mapping_grid2part

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_part2grid (this, sdata, gdata)

   USE MOD_Precision
   USE MOD_Block
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_SPMD_Task
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   type(pointer_real8_1d),    intent(in)    :: sdata(:)
   type(block_data_real8_2d), intent(inout) :: gdata

   ! Local variables
   integer :: iproc, idest, isrc
   integer :: ig, ilon, ilat, xblk, yblk, xloc, yloc, iloc, iset, ipart
   integer :: iblkme

   real(r8), allocatable :: gbuff(:)
   type(pointer_real8_1d), allocatable :: pbuff(:)

      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               allocate (pbuff(iproc)%val (this%glist(iproc)%ng))
               pbuff(iproc)%val(:) = 0.0
            ENDIF
         ENDDO

         DO iset = 1, this%npset
            DO ipart = 1, this%npart(iset)
               iproc = this%address(iset)%val(1,ipart)
               iloc  = this%address(iset)%val(2,ipart)

               pbuff(iproc)%val(iloc) = pbuff(iproc)%val(iloc) &
                  + sdata(iset)%val(ipart) * this%areapart(iset)%val(ipart)
            ENDDO
         ENDDO

#ifdef USEMPI
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               idest = p_address_io(iproc)
               CALL mpi_send (pbuff(iproc)%val, this%glist(iproc)%ng, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
            ENDIF
         ENDDO
#endif

      ENDIF

      IF (p_is_io) THEN

         CALL flush_block_data (gdata, 0.0_r8)

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_worker(iproc)
               CALL mpi_recv (gbuff, this%glist(iproc)%ng, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               gbuff = pbuff(0)%val
#endif

               DO ig = 1, this%glist(iproc)%ng
                  ilon = this%glist(iproc)%ilon(ig)
                  ilat = this%glist(iproc)%ilat(ig)
                  xblk = this%grid%xblk (ilon)
                  yblk = this%grid%yblk (ilat)
                  xloc = this%grid%xloc (ilon)
                  yloc = this%grid%yloc (ilat)

                  gdata%blk(xblk,yblk)%val(xloc,yloc) = &
                     gdata%blk(xblk,yblk)%val(xloc,yloc) + gbuff(ig)
               ENDDO

               deallocate (gbuff)
            ENDIF
         ENDDO

         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            WHERE (this%areagrid%blk(xblk,yblk)%val > 0)
               gdata%blk(xblk,yblk)%val = &
                  gdata%blk(xblk,yblk)%val / this%areagrid%blk(xblk,yblk)%val
            ELSEWHERE
               gdata%blk(xblk,yblk)%val = this%missing_value
            ENDWHERE
         ENDDO

      ENDIF

      IF (p_is_worker) THEN
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO
         deallocate (pbuff)
      ENDIF

   END SUBROUTINE spatial_mapping_part2grid

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_normalize (this, gdata, sdata)

   USE MOD_Precision
   USE MOD_Block
   USE MOD_Grid
   USE MOD_Pixelset
   USE MOD_DataType
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   type(block_data_real8_2d), intent(in)    :: gdata
   type(pointer_real8_1d),    intent(inout) :: sdata(:)

   ! Local variables
   integer :: iblkme, xblk, yblk, iset, ipart

   type(block_data_real8_2d) :: sumdata
   type(pointer_real8_1d), allocatable :: scaldata(:)


      IF (p_is_io)     CALL allocate_block_data (this%grid, sumdata)
      IF (p_is_worker) CALL this%allocate_part  (scaldata)

      CALL this%part2grid (sdata, sumdata)

      IF (p_is_io) THEN

         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            WHERE (sumdata%blk(xblk,yblk)%val /= this%missing_value)
               sumdata%blk(xblk,yblk)%val = gdata%blk(xblk,yblk)%val / sumdata%blk(xblk,yblk)%val
            ENDWHERE
         ENDDO

      ENDIF

      CALL this%grid2part (sumdata, scaldata)

      IF (p_is_worker) THEN

         DO iset = 1, this%npset
            DO ipart = 1, this%npart(iset)
               IF (this%areapart(iset)%val(ipart) > 0.) THEN
                  sdata(iset)%val(ipart) = sdata(iset)%val(ipart) * scaldata(iset)%val(ipart)
               ELSE
                  sdata(iset)%val(ipart) = this%missing_value
               ENDIF
            ENDDO
         ENDDO

      ENDIF

      IF (p_is_worker) deallocate(scaldata)

   END SUBROUTINE spatial_mapping_normalize

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_part2pset (this, sdata, pdata)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   type(pointer_real8_1d), intent(in)  :: sdata(:)
   real(r8),               intent(out) :: pdata(:)

   ! Local variables
   integer :: iset

      IF (p_is_worker) THEN

         pdata(:) = spval

         DO iset = 1, this%npset
            IF (this%areapset(iset) > 0) THEN
               pdata(iset) = sum(sdata(iset)%val * this%areapart(iset)%val) / this%areapset(iset)
            ENDIF
         ENDDO

      ENDIF

   END SUBROUTINE spatial_mapping_part2pset

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_allocate_part (this, datapart)

   USE MOD_SPMD_Task
   USE MOD_DataType
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   type(pointer_real8_1d), allocatable :: datapart (:)

   ! Local variables
   integer :: iset

      IF (p_is_worker) THEN

         IF (this%npset > 0) THEN
            allocate (datapart (this%npset))
         ENDIF

         DO iset = 1, this%npset
            IF (this%npart(iset) > 0) THEN
               allocate (datapart(iset)%val (this%npart(iset)))
            ENDIF
         ENDDO

      ENDIF

   END SUBROUTINE spatial_mapping_allocate_part

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_deallocate_part (this, datapart)

   USE MOD_SPMD_Task
   USE MOD_DataType
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   type(pointer_real8_1d), allocatable :: datapart (:)

   ! Local variables
   integer :: iset

      IF (p_is_worker) THEN

         DO iset = 1, this%npset
            IF (this%npart(iset) > 0) THEN
               deallocate (datapart(iset)%val)
            ENDIF
         ENDDO

         IF (this%npset > 0) THEN
            deallocate (datapart)
         ENDIF

      ENDIF

   END SUBROUTINE spatial_mapping_deallocate_part

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_free_mem (this)

   USE MOD_SPMD_Task
   IMPLICIT NONE

   type (spatial_mapping_type) :: this

   ! Local variables
   integer :: iproc, iset

      IF (allocated (this%grid%xblk))   deallocate (this%grid%xblk)
      IF (allocated (this%grid%yblk))   deallocate (this%grid%yblk)

      IF (allocated (this%grid%xloc))   deallocate (this%grid%xloc)
      IF (allocated (this%grid%yloc))   deallocate (this%grid%yloc)

      IF (allocated (this%grid%xcnt))   deallocate (this%grid%xcnt)
      IF (allocated (this%grid%ycnt))   deallocate (this%grid%ycnt)

      IF (allocated(this%glist)) THEN
         DO iproc = lbound(this%glist,1), ubound(this%glist,1)
            IF (allocated(this%glist(iproc)%ilat)) deallocate (this%glist(iproc)%ilat)
            IF (allocated(this%glist(iproc)%ilon)) deallocate (this%glist(iproc)%ilon)
         ENDDO

         deallocate (this%glist)
      ENDIF

      IF (p_is_worker) THEN

         IF (allocated(this%npart)) deallocate(this%npart)

         IF (allocated(this%address)) THEN
            DO iset = lbound(this%address,1), ubound(this%address,1)
               IF (allocated(this%address(iset)%val)) THEN
                  deallocate (this%address(iset)%val)
               ENDIF
            ENDDO

            deallocate (this%address)
         ENDIF

         IF (allocated(this%areapart)) THEN
            DO iset = lbound(this%areapart,1), ubound(this%areapart,1)
               IF (allocated(this%areapart(iset)%val)) THEN
                  deallocate (this%areapart(iset)%val)
               ENDIF
            ENDDO

            deallocate (this%areapart)
         ENDIF

         IF (allocated(this%areapset)) deallocate(this%areapset)

      ENDIF

   END SUBROUTINE spatial_mapping_free_mem

   SUBROUTINE forc_free_mem_spatial_mapping(this)

   USE MOD_SPMD_Task
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   ! Local variables
   integer :: iproc, iset

      IF (allocated (this%grid%xblk))   deallocate (this%grid%xblk)
      IF (allocated (this%grid%yblk))   deallocate (this%grid%yblk)

      IF (allocated (this%grid%xloc))   deallocate (this%grid%xloc)
      IF (allocated (this%grid%yloc))   deallocate (this%grid%yloc)

      IF (allocated (this%grid%xcnt))   deallocate (this%grid%xcnt)
      IF (allocated (this%grid%ycnt))   deallocate (this%grid%ycnt)

      IF (allocated(this%glist)) THEN
         DO iproc = lbound(this%glist,1), ubound(this%glist,1)
            IF (allocated(this%glist(iproc)%ilat)) deallocate (this%glist(iproc)%ilat)
            IF (allocated(this%glist(iproc)%ilon)) deallocate (this%glist(iproc)%ilon)
         ENDDO

         deallocate (this%glist)
      ENDIF

      IF (p_is_worker) THEN

         IF (allocated(this%npart)) deallocate(this%npart)

         IF (allocated(this%address)) THEN
            DO iset = lbound(this%address,1), ubound(this%address,1)
               IF (allocated(this%address(iset)%val)) THEN
                  deallocate (this%address(iset)%val)
               ENDIF
            ENDDO

            deallocate (this%address)
         ENDIF

         IF (allocated(this%areapart)) THEN
            DO iset = lbound(this%areapart,1), ubound(this%areapart,1)
               IF (allocated(this%areapart(iset)%val)) THEN
                  deallocate (this%areapart(iset)%val)
               ENDIF
            ENDDO

            deallocate (this%areapart)
         ENDIF

         IF (allocated(this%areapset)) deallocate(this%areapset)

      ENDIF

   END SUBROUTINE forc_free_mem_spatial_mapping

END MODULE MOD_SpatialMapping
