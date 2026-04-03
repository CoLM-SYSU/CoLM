#include <define.h>

#ifdef SrfdataDiag
MODULE MOD_SrfdataDiag
!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!    This module includes subroutines for checking the results of making
!    surface data.
!
!    The surface data in vector form is mapped to gridded data with last three
!    dimensions of [type,longitude,latitude], which can be viewed by other
!    softwares.
!
!    In GRIDBASED, the grid of gridded data is just the grid of the mesh.  In
!    UNSTRUCTURED or CATCHMENT, the grid is user defined and the mapping uses
!    area weighted scheme.
!
!  Created by Shupeng Zhang, May 2023
!
! !REVISIONS:
! TODO
!-----------------------------------------------------------------------

   USE MOD_Grid
   USE MOD_SpatialMapping

   IMPLICIT NONE

   ! PUBLIC variables and subroutines
   type(grid_type) :: gdiag

   type(spatial_mapping_type) :: m_elm2diag

   type(spatial_mapping_type) :: m_patch2diag
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   type(spatial_mapping_type) :: m_pft2diag
#endif
#ifdef URBAN_MODEL
   type(spatial_mapping_type) :: m_urb2diag
#endif

   PUBLIC :: srfdata_diag_init
   PUBLIC :: srfdata_map_and_write

   ! PRIVATE
   type(grid_concat_type), PRIVATE :: srf_concat
   integer, PRIVATE :: srf_data_id

CONTAINS

   ! ------ SUBROUTINE ------
   SUBROUTINE srfdata_diag_init (dir_landdata, lc_year)

   USE MOD_SPMD_Task
   USE MOD_LandElm
   USE MOD_LandPatch
#ifdef CROP
   USE MOD_LandCrop
#endif
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_LandPFT
#endif
#ifdef URBAN_MODEL
   USE MOD_LandUrban
#endif

   IMPLICIT NONE

   character(len=256), intent(in) :: dir_landdata
   integer           , intent(in) :: lc_year

   ! Local Variables
   character(len=256) :: landdir, landname
   integer :: ityp
   integer :: typindex(N_land_classification+1)
   real(r8), allocatable :: elmid_r8(:)
   character(len=4) :: cyear

      landdir = trim(dir_landdata) // '/diag/'
      IF (p_is_master) THEN
         CALL system('mkdir -p ' // trim(adjustl(landdir)))
      ENDIF

      CALL srf_concat%set (gdiag)

      CALL m_elm2diag%build_arealweighted (gdiag, landelm)

      CALL m_patch2diag%build_arealweighted (gdiag, landpatch)

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      CALL m_pft2diag%build_arealweighted (gdiag, landpft)
#endif

#ifdef URBAN_MODEL
      CALL m_urb2diag%build_arealweighted (gdiag, landurban)
#endif

      srf_data_id = 666

      IF (p_is_worker) THEN
         allocate (elmid_r8 (landelm%nset)); elmid_r8 = real(landelm%eindex, r8)
      ENDIF

      write(cyear,'(i4.4)') lc_year
      landname = trim(dir_landdata)//'/diag/element_'//trim(cyear)//'.nc'
      CALL srfdata_map_and_write (elmid_r8, landelm%settyp, (/0,1/), m_elm2diag, &
         -1.0e36_r8, landname, 'element', compress = 1, write_mode = 'one', &
         defval=0._r8, create_mode=.true.)

      IF (p_is_worker) deallocate (elmid_r8)

      typindex = (/(ityp, ityp = 0, N_land_classification)/)
      landname = trim(dir_landdata)//'/diag/patchfrac_elm_'//trim(cyear)//'.nc'

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            IF (.not. allocated(landpatch%pctshared)) THEN
               allocate (landpatch%pctshared (numpatch))
               landpatch%pctshared = 1.
            ENDIF
         ENDIF
      ENDIF

      CALL srfdata_map_and_write (landpatch%pctshared, landpatch%settyp, typindex, m_patch2diag, &
         -1.0e36_r8, landname, 'patchfrac_elm', compress = 1, write_mode = 'one', defval=0._r8, &
         stat_mode = 'fraction', pctshared = landpatch%pctshared, create_mode=.true.)

   END SUBROUTINE srfdata_diag_init

   ! ------ SUBROUTINE ------
   SUBROUTINE srfdata_map_and_write ( &
         vsrfdata,  settyp,   typindex,   m_srf,       spv,          filename, &
         dataname,  compress, write_mode, lastdimname, lastdimvalue, defval,   &
         stat_mode, pctshared, create_mode)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Block
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_NetCDFSerial
   IMPLICIT NONE

   real(r8), intent(in) :: vsrfdata (:)
   integer , intent(in) :: settyp   (:)
   integer , intent(in) :: typindex (:)

   type(spatial_mapping_type), intent(in) :: m_srf

   real(r8), intent(in) :: spv

   character (len=*), intent(in) :: filename
   character (len=*), intent(in) :: dataname
   integer, intent(in) :: compress

   character (len=*), intent(in), optional :: write_mode

   character (len=*), intent(in), optional :: lastdimname
   integer,  intent(in), optional :: lastdimvalue

   real(r8), intent(in), optional :: defval

   character (len=*), intent(in), optional :: stat_mode
   real(r8), intent(in), optional :: pctshared (:)
   logical , intent(in), optional :: create_mode

   ! Local variables
   type(block_data_real8_3d) :: wdata, wtone
   type(block_data_real8_2d) :: wdsum, wtsum
   real(r8), allocatable :: vecone   (:)
   real(r8), allocatable :: areafrac (:)
   real(r8), allocatable :: vdata (:,:,:), vdsum (:,:)
   real(r8), allocatable :: rbuf  (:,:,:), sbuf  (:,:,:)

   character(len=10) :: wmode, smode
   integer :: iblkme, ib, jb, iblk, jblk, idata, ixseg, iyseg
   integer :: ntyps, ityp, xcnt, ycnt, xbdsp, ybdsp, xgdsp, ygdsp
   integer :: rmesg(3), smesg(3), isrc
   character(len=256) :: fileblock
   logical :: fexists, cmode
   integer :: ilastdim

      IF (present(write_mode)) THEN
         wmode = trim(write_mode)
      ELSE
         wmode = 'one'
      ENDIF

      IF (present(stat_mode)) THEN
         smode = trim(stat_mode)
      ELSE
         smode = 'mean'
      ENDIF

      IF (present(create_mode)) THEN
         cmode = create_mode
      ELSE
         cmode = .false.
      ENDIF

      ntyps = size(typindex)

      IF (p_is_io) THEN
         CALL allocate_block_data (gdiag, wtone, ntyps)
         CALL allocate_block_data (gdiag, wdata, ntyps)
         CALL allocate_block_data (gdiag, wdsum)
         CALL allocate_block_data (gdiag, wtsum)
      ENDIF

      IF (p_is_worker) THEN
         IF (m_srf%npset > 0) THEN
            allocate (vecone (m_srf%npset))
            vecone(:) = 1.0
         ENDIF
      ENDIF

      CALL m_srf%pset2grid_split (vecone, settyp, typindex, wtone, spv)

      IF (trim(smode) == 'mean') THEN
         CALL m_srf%pset2grid_split (vsrfdata, settyp, typindex, wdata, spv)
      ELSEIF (trim(smode) == 'fraction') THEN
         IF (p_is_worker) THEN
            IF (m_srf%npset > 0) THEN
               allocate (areafrac (m_srf%npset))
               areafrac = vsrfdata
               IF (present(pctshared)) THEN
                  WHERE (pctshared > 0.) areafrac = areafrac / pctshared
               ENDIF
            ENDIF
         ENDIF
         CALL m_srf%pset2grid_split (areafrac, settyp, typindex, wdata, spv)
      ENDIF

      IF (p_is_io) THEN
         DO iblkme = 1, gblock%nblkme
            ib = gblock%xblkme(iblkme)
            jb = gblock%yblkme(iblkme)

            wtsum%blk(ib,jb)%val = sum(wtone%blk(ib,jb)%val, &
               mask = (wtone%blk(ib,jb)%val > 0.) .and. (wdata%blk(ib,jb)%val /= spv), dim = 1)
            wdsum%blk(ib,jb)%val = sum(wdata%blk(ib,jb)%val, &
               mask = (wtone%blk(ib,jb)%val > 0.) .and. (wdata%blk(ib,jb)%val /= spv), dim = 1)

            IF (trim(smode) == 'mean') THEN

               WHERE ((wtone%blk(ib,jb)%val > 0.) .and. (wdata%blk(ib,jb)%val /= spv))
                  wdata%blk(ib,jb)%val = wdata%blk(ib,jb)%val / wtone%blk(ib,jb)%val
               ELSEWHERE
                  wdata%blk(ib,jb)%val = spv
               END WHERE

            ELSEIF (trim(smode) == 'fraction') THEN

               DO ityp = 1, ntyps
                  WHERE ((wtsum%blk(ib,jb)%val(:,:) > 0.) .and. (wtone%blk(ib,jb)%val(ityp,:,:) /= spv))
                     wdata%blk(ib,jb)%val(ityp,:,:) = wtone%blk(ib,jb)%val(ityp,:,:) / wtsum%blk(ib,jb)%val
                  ELSEWHERE
                     wdata%blk(ib,jb)%val(ityp,:,:) = spv
                  END WHERE
               ENDDO

            ENDIF

            WHERE (wtsum%blk(ib,jb)%val > 0.)
               wdsum%blk(ib,jb)%val = wdsum%blk(ib,jb)%val / wtsum%blk(ib,jb)%val
            ELSEWHERE
               wdsum%blk(ib,jb)%val = spv
            END WHERE

         ENDDO
      ENDIF

      IF (trim(wmode) == 'one') THEN

         IF (p_is_master) THEN

            allocate (vdata (srf_concat%ginfo%nlon, srf_concat%ginfo%nlat, ntyps))
            vdata(:,:,:) = spv

            allocate (vdsum (srf_concat%ginfo%nlon, srf_concat%ginfo%nlat))
            vdsum(:,:) = spv

#ifdef USEMPI
            DO idata = 1, srf_concat%ndatablk

               CALL mpi_recv (rmesg, 3, MPI_INTEGER, MPI_ANY_SOURCE, &
                  srf_data_id, p_comm_glb, p_stat, p_err)

               isrc  = rmesg(1)
               ixseg = rmesg(2)
               iyseg = rmesg(3)

               xgdsp = srf_concat%xsegs(ixseg)%gdsp
               ygdsp = srf_concat%ysegs(iyseg)%gdsp
               xcnt  = srf_concat%xsegs(ixseg)%cnt
               ycnt  = srf_concat%ysegs(iyseg)%cnt

               allocate (rbuf (ntyps,xcnt,ycnt))

               CALL mpi_recv (rbuf, ntyps * xcnt * ycnt, MPI_REAL8, &
                  isrc, srf_data_id, p_comm_glb, p_stat, p_err)

               DO ityp = 1, ntyps
                  vdata (xgdsp+1:xgdsp+xcnt,ygdsp+1:ygdsp+ycnt,ityp) = rbuf(ityp,:,:)
               ENDDO

               CALL mpi_recv (rbuf(1,:,:), xcnt * ycnt, MPI_REAL8, &
                  isrc, srf_data_id, p_comm_glb, p_stat, p_err)

               vdsum (xgdsp+1:xgdsp+xcnt,ygdsp+1:ygdsp+ycnt) = rbuf(1,:,:)

               deallocate (rbuf)
            ENDDO
#else
            DO iyseg = 1, srf_concat%nyseg
               DO ixseg = 1, srf_concat%nxseg
                  iblk  = srf_concat%xsegs(ixseg)%blk
                  jblk  = srf_concat%ysegs(iyseg)%blk
                  xbdsp = srf_concat%xsegs(ixseg)%bdsp
                  ybdsp = srf_concat%ysegs(iyseg)%bdsp
                  xgdsp = srf_concat%xsegs(ixseg)%gdsp
                  ygdsp = srf_concat%ysegs(iyseg)%gdsp
                  xcnt  = srf_concat%xsegs(ixseg)%cnt
                  ycnt  = srf_concat%ysegs(iyseg)%cnt

                  DO ityp = 1, ntyps
                     vdata (xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt,ityp) = &
                        wdata%blk(iblk,jblk)%val(ityp,xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)
                  ENDDO

                  vdsum (xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt) = &
                     wdsum%blk(iblk,jblk)%val(xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)
               ENDDO
            ENDDO
#endif

            IF (present(defval)) THEN
               WHERE (vdata == spv)  vdata = defval
               WHERE (vdsum == spv)  vdsum = defval
            ENDIF

            write(*,*) 'Please check gridded data < ', trim(dataname), ' > in ', trim(filename)

            inquire (file=trim(filename), exist=fexists)
            IF (.not. fexists .or. cmode) THEN

               CALL ncio_create_file (filename)

               IF (ntyps > 1) THEN
                  CALL ncio_define_dimension (filename, 'TypeIndex', ntyps)
                  CALL ncio_write_serial (filename, 'TypeIndex', typindex, 'TypeIndex')
               ENDIF

               CALL ncio_define_dimension (filename, 'lon' , srf_concat%ginfo%nlon)
               CALL ncio_define_dimension (filename, 'lat' , srf_concat%ginfo%nlat)

               CALL ncio_write_serial (filename, 'lat', srf_concat%ginfo%lat_c, 'lat')
               CALL ncio_put_attr (filename, 'lat', 'long_name', 'latitude')
               CALL ncio_put_attr (filename, 'lat', 'units', 'degrees_north')

               CALL ncio_write_serial (filename, 'lon', srf_concat%ginfo%lon_c, 'lon')
               CALL ncio_put_attr (filename, 'lon', 'long_name', 'longitude')
               CALL ncio_put_attr (filename, 'lon', 'units', 'degrees_east')

               CALL ncio_write_serial (filename, 'lat_s', srf_concat%ginfo%lat_s, 'lat')
               CALL ncio_put_attr (filename, 'lat_s', 'long_name', 'southern latitude boundary')
               CALL ncio_put_attr (filename, 'lat_s', 'units', 'degrees_north')

               CALL ncio_write_serial (filename, 'lat_n', srf_concat%ginfo%lat_n, 'lat')
               CALL ncio_put_attr (filename, 'lat_n', 'long_name', 'northern latitude boundary')
               CALL ncio_put_attr (filename, 'lat_n', 'units', 'degrees_north')

               CALL ncio_write_serial (filename, 'lon_w', srf_concat%ginfo%lon_w, 'lon')
               CALL ncio_put_attr (filename, 'lon_w', 'long_name', 'western longitude boundary')
               CALL ncio_put_attr (filename, 'lon_w', 'units', 'degrees_east')

               CALL ncio_write_serial (filename, 'lon_e', srf_concat%ginfo%lon_e, 'lon')
               CALL ncio_put_attr (filename, 'lon_e', 'long_name', 'eastern longitude boundary')
               CALL ncio_put_attr (filename, 'lon_e', 'units', 'degrees_east')

            ENDIF

            IF (present(lastdimname) .and. present(lastdimvalue)) THEN

               CALL ncio_write_lastdim (filename, lastdimname, lastdimvalue, ilastdim)

               IF (ntyps > 1) THEN
                  CALL ncio_write_serial_time (filename, dataname, ilastdim, vdata, &
                     'lon', 'lat', 'TypeIndex', trim(lastdimname), compress)
                  CALL ncio_write_serial_time (filename, trim(dataname)//'_grid', ilastdim, vdsum, &
                     'lon', 'lat', trim(lastdimname), compress)
               ELSE
                  CALL ncio_write_serial_time (filename, dataname, ilastdim, vdata(:,:,1), &
                     'lon', 'lat', trim(lastdimname), compress)
               ENDIF

            ELSE

               IF (ntyps > 1) THEN
                  CALL ncio_write_serial (filename, dataname, vdata, 'lon', 'lat', 'TypeIndex', compress)
                  CALL ncio_write_serial (filename, trim(dataname)//'_grid', vdsum, 'lon', 'lat', compress)
               ELSE
                  CALL ncio_write_serial (filename, dataname, vdata(:,:,1), 'lon', 'lat', compress)
               ENDIF

            ENDIF

            CALL ncio_put_attr (filename, dataname, 'missing_value', spv)

            IF (ntyps > 1) THEN
               CALL ncio_put_attr (filename, trim(dataname)//'_grid', 'missing_value', spv)
            ENDIF

            deallocate (vdata)
            deallocate (vdsum)

         ENDIF

#ifdef USEMPI
         IF (p_is_io) THEN

            DO iyseg = 1, srf_concat%nyseg
               DO ixseg = 1, srf_concat%nxseg

                  iblk = srf_concat%xsegs(ixseg)%blk
                  jblk = srf_concat%ysegs(iyseg)%blk

                  IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN

                     xbdsp = srf_concat%xsegs(ixseg)%bdsp
                     ybdsp = srf_concat%ysegs(iyseg)%bdsp
                     xcnt  = srf_concat%xsegs(ixseg)%cnt
                     ycnt  = srf_concat%ysegs(iyseg)%cnt

                     smesg = (/p_iam_glb, ixseg, iyseg/)
                     CALL mpi_send (smesg, 3, MPI_INTEGER, &
                        p_address_master, srf_data_id, p_comm_glb, p_err)

                     allocate (sbuf (ntyps,xcnt,ycnt))

                     sbuf = wdata%blk(iblk,jblk)%val(:,xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)
                     CALL mpi_send (sbuf, ntyps*xcnt*ycnt, MPI_REAL8, &
                        p_address_master, srf_data_id, p_comm_glb, p_err)

                     sbuf(1,:,:) = wdsum%blk(iblk,jblk)%val(xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)
                     CALL mpi_send (sbuf(1,:,:), xcnt*ycnt, MPI_REAL8, &
                        p_address_master, srf_data_id, p_comm_glb, p_err)

                     deallocate (sbuf)
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
#endif

         srf_data_id = srf_data_id + 1

      ELSEIF (trim(wmode) == 'block') THEN

         IF (p_is_io) THEN

            DO iblkme = 1, gblock%nblkme
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)

               IF ((gdiag%xcnt(iblk) == 0) .or. (gdiag%ycnt(jblk) == 0)) CYCLE

               CALL get_filename_block (filename, iblk, jblk, fileblock)

               inquire (file=trim(filename), exist=fexists)
               IF (.not. fexists) THEN
                  CALL ncio_create_file (fileblock)
                  CALL ncio_define_dimension (fileblock, 'TypeIndex', ntyps)
                  CALL srf_write_grid_info   (fileblock, gdiag, iblk, jblk)
               ENDIF

               IF (present(lastdimname) .and. present(lastdimvalue)) THEN
                  CALL ncio_write_lastdim (fileblock, lastdimname, lastdimvalue, ilastdim)
                  CALL ncio_write_serial_time (fileblock, dataname, ilastdim, wdata%blk(iblk,jblk)%val, &
                     'TypeIndex', 'lon', 'lat', trim(lastdimname), compress)
               ELSE
                  CALL ncio_write_serial (fileblock, dataname, &
                     wdata%blk(iblk,jblk)%val, 'TypeIndex', 'lon', 'lat', compress)
               ENDIF

               CALL ncio_put_attr (fileblock, dataname, 'missing_value', spv)

            ENDDO

         ENDIF
      ENDIF

      IF (allocated(vecone)) deallocate(vecone)

   END SUBROUTINE srfdata_map_and_write

   !------------------
   SUBROUTINE srf_write_grid_info (fileblock, grid, iblk, jblk)

   USE MOD_Block
   USE MOD_Grid
   USE MOD_NetCDFSerial
   IMPLICIT NONE

   character(len=*), intent(in) :: fileblock
   type (grid_type), intent(in) :: grid
   integer, intent(in) :: iblk, jblk

   ! Local variable
   integer :: yl, yu, xl, xu, nx
   real(r8), allocatable :: lat_s(:), lat_n(:), lon_w(:), lon_e(:)

      allocate (lon_w (grid%xcnt(iblk)))
      allocate (lon_e (grid%xcnt(iblk)))
      allocate (lat_s (grid%ycnt(jblk)))
      allocate (lat_n (grid%ycnt(jblk)))

      yl = grid%ydsp(jblk) + 1
      yu = grid%ydsp(jblk) + grid%ycnt(jblk)

      lat_s = grid%lat_s(yl:yu)
      lat_n = grid%lat_n(yl:yu)

      IF (grid%xdsp(iblk) + grid%xcnt(iblk) > grid%nlon) THEN
         xl = grid%xdsp(iblk) + 1
         xu = grid%nlon
         nx = grid%nlon - grid%xdsp(iblk)
         lon_w(1:nx) = grid%lon_w(xl:xu)
         lon_e(1:nx) = grid%lon_e(xl:xu)

         xl = 1
         xu = grid%xcnt(iblk) - nx
         lon_w(nx+1:grid%xcnt(iblk)) = grid%lon_w(xl:xu)
         lon_e(nx+1:grid%xcnt(iblk)) = grid%lon_e(xl:xu)
      ELSE
         xl = grid%xdsp(iblk) + 1
         xu = grid%xdsp(iblk) + grid%xcnt(iblk)
         lon_w = grid%lon_w(xl:xu)
         lon_e = grid%lon_e(xl:xu)
      ENDIF

      CALL ncio_define_dimension (fileblock, 'lat', grid%ycnt(jblk))
      CALL ncio_define_dimension (fileblock, 'lon', grid%xcnt(iblk))
      CALL ncio_write_serial (fileblock, 'lat_s', lat_s, 'lat')
      CALL ncio_write_serial (fileblock, 'lat_n', lat_n, 'lat')
      CALL ncio_write_serial (fileblock, 'lon_w', lon_w, 'lon')
      CALL ncio_write_serial (fileblock, 'lon_e', lon_e, 'lon')

   END SUBROUTINE srf_write_grid_info

END MODULE MOD_SrfdataDiag
#endif
