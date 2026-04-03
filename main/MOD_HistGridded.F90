#include <define.h>

MODULE MOD_HistGridded

!----------------------------------------------------------------------------
! !DESCRIPTION:
!
!     Write out gridded model results to history files.
!
!  Original version: Yongjiu Dai, September 15, 1999, 03/2014
!
! !REVISIONS:
!  Shupeng Zhang, 05/2023: 1) porting codes to MPI parallel version
!
!  TODO...(need complement)
!----------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_SpatialMapping
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_SPMD_Task
#ifdef USEMPI
   USE MOD_HistWriteBack
#endif

   type(grid_type), target :: ghist
   type(spatial_mapping_type) :: mp2g_hist
   type(spatial_mapping_type) :: mp2g_hist_urb

   type(block_data_real8_2d) :: landfraction

   type(grid_concat_type) :: hist_concat

   integer :: hist_data_id

!--------------------------------------------------------------------------
CONTAINS

   SUBROUTINE hist_gridded_init (dir_hist, lulcc_call)

   USE MOD_SPMD_Task
   USE MOD_Vars_Global
   USE MOD_Block
   USE MOD_LandPatch
#ifdef URBAN_MODEL
   USE MOD_LandUrban
#endif
   USE MOD_Vars_1DAccFluxes
   USE MOD_Forcing, only: gforc
#ifdef SinglePoint
   USE MOD_SingleSrfData
#endif
   USE MOD_Utils
   IMPLICIT NONE

   character(len=*) , intent(in) :: dir_hist
   logical, optional, intent(in) :: lulcc_call

   ! Local Variables
   type(block_data_real8_2d) :: gridarea
   integer :: iblkme, xblk, yblk, xloc, yloc, xglb, yglb

      IF (DEF_hist_grid_as_forcing) THEN
         CALL ghist%define_by_copy (gforc)
      ELSE
         CALL ghist%define_by_res (DEF_hist_lon_res, DEF_hist_lat_res)
      ENDIF

      IF (present(lulcc_call)) CALL mp2g_hist%forc_free_mem
      CALL mp2g_hist%build_arealweighted (ghist, landpatch)

#ifdef URBAN_MODEL
      IF (present(lulcc_call)) CALL mp2g_hist_urb%forc_free_mem
      CALL mp2g_hist_urb%build_arealweighted (ghist, landurban)
#endif

      IF (p_is_io) THEN
         CALL allocate_block_data (ghist, landfraction)
         CALL allocate_block_data (ghist, gridarea)

         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)
            DO yloc = 1, ghist%ycnt(yblk)
               DO xloc = 1, ghist%xcnt(xblk)
                  xglb = ghist%xdsp(xblk) + xloc
                  yglb = ghist%ydsp(yblk) + yloc
                  gridarea%blk(xblk,yblk)%val(xloc,yloc) = areaquad ( &
                     ghist%lat_s(yglb), ghist%lat_n(yglb), ghist%lon_w(xglb), ghist%lon_e(xglb))
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      CALL mp2g_hist%get_sumarea (landfraction)
      CALL block_data_division   (landfraction, gridarea)

      CALL hist_concat%set (ghist)
#ifdef SinglePoint
      hist_concat%ginfo%lat_c(:) = SITE_lat_location
      hist_concat%ginfo%lon_c(:) = SITE_lon_location
#endif

      IF (trim(DEF_HIST_mode) == 'one') THEN
         hist_data_id = 10001
      ENDIF

   END SUBROUTINE hist_gridded_init


   SUBROUTINE flux_map_and_write_2d ( &
         acc_vec, file_hist, varname, itime_in_file, sumarea, filter, &
         longname, units, input_mode)

   USE MOD_Block
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   real(r8), intent(inout) :: acc_vec(:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer,          intent(in) :: itime_in_file
   character(len=*), intent(in) :: longname
   character(len=*), intent(in) :: units

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)
   character(len=*), intent(in), optional :: input_mode

   ! Local variables
   type(block_data_real8_2d) :: flux_xy_2d
   integer :: iblkme, xblk, yblk, xloc, yloc
   integer :: compress

      IF (p_is_io)      CALL allocate_block_data (ghist, flux_xy_2d)

      IF (present(input_mode)) THEN
         CALL mp2g_hist%pset2grid (acc_vec, flux_xy_2d, spv = spval, msk = filter, input_mode = input_mode)
      ELSE
         CALL mp2g_hist%pset2grid (acc_vec, flux_xy_2d, spv = spval, msk = filter)
      ENDIF

      IF (p_is_io) THEN
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            DO yloc = 1, ghist%ycnt(yblk)
               DO xloc = 1, ghist%xcnt(xblk)

                  IF (sumarea%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) THEN
                     IF (flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) /= spval) THEN
                        flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) &
                           = flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) &
                           / sumarea%blk(xblk,yblk)%val(xloc,yloc)
                     ENDIF
                  ELSE
                     flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) = spval
                  ENDIF

               ENDDO
            ENDDO

         ENDDO
      ENDIF

      compress = DEF_HIST_CompressLevel
      CALL hist_write_var_real8_2d (file_hist, varname, ghist, itime_in_file, &
         flux_xy_2d, compress, longname, units)

   END SUBROUTINE flux_map_and_write_2d


   SUBROUTINE flux_map_and_write_urb_2d ( &
         acc_vec, file_hist, varname, itime_in_file, sumarea, filter, &
         longname, units)

   USE MOD_Block
   USE MOD_Vars_1DAccFluxes,  only: nac
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   real(r8), intent(inout) :: acc_vec(:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer,          intent(in) :: itime_in_file
   character(len=*), intent(in) :: longname
   character(len=*), intent(in) :: units

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)

   ! Local variables
   type(block_data_real8_2d) :: flux_xy_2d
   integer :: iblkme, xblk, yblk, xloc, yloc
   integer :: compress

      IF (p_is_worker)  WHERE (acc_vec /= spval)  acc_vec = acc_vec / nac
      IF (p_is_io)      CALL allocate_block_data (ghist, flux_xy_2d)

      CALL mp2g_hist_urb%pset2grid (acc_vec, flux_xy_2d, spv = spval, msk = filter)

      IF (p_is_io) THEN
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            DO yloc = 1, ghist%ycnt(yblk)
               DO xloc = 1, ghist%xcnt(xblk)

                  IF (sumarea%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) THEN
                     IF (flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) /= spval) THEN
                        flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) &
                           = flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) &
                           / sumarea%blk(xblk,yblk)%val(xloc,yloc)
                     ENDIF
                  ELSE
                     flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) = spval
                  ENDIF

               ENDDO
            ENDDO

         ENDDO
      ENDIF

      compress = DEF_HIST_CompressLevel
      CALL hist_write_var_real8_2d (file_hist, varname, ghist, itime_in_file, flux_xy_2d, &
         compress, longname, units)

   END SUBROUTINE flux_map_and_write_urb_2d


   SUBROUTINE flux_map_and_write_3d ( &
         acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, sumarea, filter, &
         longname, units)

   USE MOD_Block
   USE MOD_Vars_1DAccFluxes,  only: nac
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   real(r8), intent(inout) :: acc_vec(:,:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer, intent(in) :: itime_in_file
   character(len=*), intent(in) :: dim1name
   integer, intent(in) :: lb1, ndim1

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)
   character (len=*), intent(in) :: longname
   character (len=*), intent(in) :: units

   ! Local variables
   type(block_data_real8_3d) :: flux_xy_3d
   integer :: iblkme, xblk, yblk, xloc, yloc, i1
   integer :: compress

      IF (p_is_worker)  THEN
         WHERE (acc_vec /= spval)  acc_vec = acc_vec / nac
      ENDIF
      IF (p_is_io) THEN
         CALL allocate_block_data (ghist, flux_xy_3d, ndim1, lb1)
      ENDIF

      CALL mp2g_hist%pset2grid (acc_vec, flux_xy_3d, spv = spval, msk = filter)

      IF (p_is_io) THEN
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            DO yloc = 1, ghist%ycnt(yblk)
               DO xloc = 1, ghist%xcnt(xblk)

                  IF (sumarea%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) THEN
                     DO i1 = flux_xy_3d%lb1, flux_xy_3d%ub1
                        IF (flux_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc) /= spval) THEN
                           flux_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc) &
                              = flux_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc) &
                              / sumarea%blk(xblk,yblk)%val(xloc,yloc)
                        ENDIF
                     ENDDO
                  ELSE
                     flux_xy_3d%blk(xblk,yblk)%val(:,xloc,yloc) = spval
                  ENDIF

               ENDDO
            ENDDO

         ENDDO
      ENDIF

      compress = DEF_HIST_CompressLevel
      CALL hist_write_var_real8_3d (file_hist, varname, dim1name, ghist, &
         itime_in_file, flux_xy_3d, compress, longname, units)

   END SUBROUTINE flux_map_and_write_3d


   SUBROUTINE flux_map_and_write_4d ( &
         acc_vec, file_hist, varname, itime_in_file, &
         dim1name, lb1, ndim1, dim2name, lb2, ndim2, &
         sumarea, filter, longname, units)

   USE MOD_Block
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   real(r8), intent(inout) :: acc_vec(:,:,:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer, intent(in) :: itime_in_file
   character(len=*), intent(in) :: dim1name, dim2name
   integer, intent(in) :: lb1, ndim1, lb2, ndim2

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)
   character (len=*), intent(in) :: longname
   character (len=*), intent(in) :: units

   ! Local variables
   type(block_data_real8_4d) :: flux_xy_4d
   integer :: iblkme, xblk, yblk, xloc, yloc, i1, i2
   integer :: compress

      IF (p_is_io) THEN
         CALL allocate_block_data (ghist, flux_xy_4d, ndim1, ndim2, lb1 = lb1, lb2 = lb2)
      ENDIF

      CALL mp2g_hist%pset2grid (acc_vec, flux_xy_4d, spv = spval, msk = filter)

      IF (p_is_io) THEN
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            DO yloc = 1, ghist%ycnt(yblk)
               DO xloc = 1, ghist%xcnt(xblk)

                  IF (sumarea%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) THEN
                     DO i1 = flux_xy_4d%lb1, flux_xy_4d%ub1
                        DO i2 = flux_xy_4d%lb2, flux_xy_4d%ub2
                           IF (flux_xy_4d%blk(xblk,yblk)%val(i1,i2,xloc,yloc) /= spval) THEN
                              flux_xy_4d%blk(xblk,yblk)%val(i1,i2,xloc,yloc) &
                                 = flux_xy_4d%blk(xblk,yblk)%val(i1,i2,xloc,yloc) &
                                 / sumarea%blk(xblk,yblk)%val(xloc,yloc)
                           ENDIF
                        ENDDO
                     ENDDO
                  ELSE
                     flux_xy_4d%blk(xblk,yblk)%val(:,:,xloc,yloc) = spval
                  ENDIF

               ENDDO
            ENDDO

         ENDDO
      ENDIF

      compress = DEF_HIST_CompressLevel
      CALL hist_write_var_real8_4d (file_hist, varname, dim1name, dim2name, &
         ghist, itime_in_file, flux_xy_4d, compress, longname, units)

   END SUBROUTINE flux_map_and_write_4d


   SUBROUTINE hist_gridded_write_time ( &
         filename, filelast, dataname, time, itime)

   USE MOD_Block
   IMPLICIT NONE

   character (len=*), intent(in) :: filename
   character (len=*), intent(in) :: filelast
   character (len=*), intent(in) :: dataname

   integer, intent(in)  :: time(3)
   integer, intent(out) :: itime

   ! Local variables
   character(len=256) :: fileblock
   integer :: iblkme, iblk, jblk
   logical :: fexists

      IF (trim(DEF_HIST_mode) == 'one') THEN
         IF (p_is_master) THEN
#ifdef USEMPI
            IF (DEF_HIST_WriteBack) THEN
               CALL hist_writeback_latlon_time (filename, filelast, dataname, time, hist_concat)
               itime = 1
            ELSE
#endif
            inquire (file=filename, exist=fexists)
            IF ((.not. fexists) .or. (trim(filename) /= trim(filelast))) THEN

               CALL ncio_create_file (trim(filename))
               CALL ncio_define_dimension(filename, 'time', 0)
               CALL ncio_define_dimension(filename, 'lat' , hist_concat%ginfo%nlat)
               CALL ncio_define_dimension(filename, 'lon' , hist_concat%ginfo%nlon)

               CALL ncio_write_serial (filename, 'lat', hist_concat%ginfo%lat_c, 'lat')
               CALL ncio_put_attr (filename, 'lat', 'long_name', 'latitude')
               CALL ncio_put_attr (filename, 'lat', 'units', 'degrees_north')

               CALL ncio_write_serial (filename, 'lon', hist_concat%ginfo%lon_c, 'lon')
               CALL ncio_put_attr (filename, 'lon', 'long_name', 'longitude')
               CALL ncio_put_attr (filename, 'lon', 'units', 'degrees_east')

#ifndef SinglePoint
               CALL ncio_write_serial (filename, 'lat_s', hist_concat%ginfo%lat_s, 'lat')
               CALL ncio_write_serial (filename, 'lat_n', hist_concat%ginfo%lat_n, 'lat')
               CALL ncio_write_serial (filename, 'lon_w', hist_concat%ginfo%lon_w, 'lon')
               CALL ncio_write_serial (filename, 'lon_e', hist_concat%ginfo%lon_e, 'lon')
#endif

               CALL ncio_write_colm_dimension (filename)

            ENDIF

            CALL ncio_write_time (filename, dataname, time, itime, DEF_HIST_FREQ)

#ifdef USEMPI
            ENDIF
#endif
         ENDIF

#ifdef USEMPI
         CALL mpi_bcast (itime, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif

      ELSEIF (trim(DEF_HIST_mode) == 'block') THEN

         IF (p_is_io) THEN

            DO iblkme = 1, gblock%nblkme
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)
               IF (ghist%ycnt(jblk) <= 0) CYCLE
               IF (ghist%xcnt(iblk) <= 0) CYCLE

               CALL get_filename_block (filename, iblk, jblk, fileblock)

               inquire (file=fileblock, exist=fexists)
               IF ((.not. fexists) .or. (trim(filename) /= trim(filelast))) THEN
                  CALL ncio_create_file (trim(fileblock))
                  CALL ncio_define_dimension (fileblock, 'time', 0)
                  CALL hist_write_grid_info  (fileblock, ghist, iblk, jblk)
               ENDIF

               CALL ncio_write_time (fileblock, dataname, time, itime, DEF_HIST_FREQ)

            ENDDO

         ENDIF
#ifdef USEMPI
         IF (.not. p_is_master) CALL mpi_bcast (itime, 1, MPI_INTEGER, p_root, p_comm_group, p_err)
#endif

      ENDIF

   END SUBROUTINE hist_gridded_write_time


   SUBROUTINE hist_write_var_real8_2d ( &
         filename, dataname, grid, itime, wdata, compress, longname, units)

   USE MOD_Block
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   character (len=*), intent(in) :: filename
   character (len=*), intent(in) :: dataname
   type (grid_type),  intent(in) :: grid
   integer, intent(in) :: itime

   type (block_data_real8_2d), intent(in) :: wdata

   integer, intent(in) :: compress
   character(len=*), intent(in) :: longname
   character(len=*), intent(in) :: units

   ! Local variables
   integer :: iblkme, iblk, jblk, idata, ixseg, iyseg
   integer :: xcnt, ycnt, xbdsp, ybdsp, xgdsp, ygdsp
   integer :: rmesg(3), smesg(3), isrc
   character(len=256) :: fileblock
   real(r8), allocatable :: rbuf(:,:), sbuf(:,:), vdata(:,:)

      IF (trim(DEF_HIST_mode) == 'one') THEN

         IF (p_is_master) THEN

#ifdef USEMPI
            IF (.not. DEF_HIST_WriteBack) THEN

               allocate (vdata (hist_concat%ginfo%nlon, hist_concat%ginfo%nlat))
               vdata(:,:) = spval

               DO idata = 1, hist_concat%ndatablk
                  CALL mpi_recv (rmesg, 3, MPI_INTEGER, MPI_ANY_SOURCE, &
                     hist_data_id, p_comm_glb, p_stat, p_err)

                  isrc  = rmesg(1)
                  ixseg = rmesg(2)
                  iyseg = rmesg(3)

                  xgdsp = hist_concat%xsegs(ixseg)%gdsp
                  ygdsp = hist_concat%ysegs(iyseg)%gdsp
                  xcnt = hist_concat%xsegs(ixseg)%cnt
                  ycnt = hist_concat%ysegs(iyseg)%cnt

                  allocate (rbuf(xcnt,ycnt))

                  CALL mpi_recv (rbuf, xcnt*ycnt, MPI_REAL8, &
                     isrc, hist_data_id, p_comm_glb, p_stat, p_err)

                  vdata (xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt) = rbuf
                  deallocate (rbuf)

               ENDDO

            ELSE
               IF (itime == -1) THEN
                  CALL hist_writeback_var_header (hist_data_id, filename, dataname, &
                     2, 'lon', 'lat', '', '', '', compress, longname, units)
               ELSE
                  CALL hist_writeback_var_header (hist_data_id, filename, dataname, &
                     3, 'lon', 'lat', 'time', '', '', compress, longname, units)
               ENDIF
            ENDIF
#else
            allocate (vdata (hist_concat%ginfo%nlon, hist_concat%ginfo%nlat))
            vdata(:,:) = spval

            DO iyseg = 1, hist_concat%nyseg
               DO ixseg = 1, hist_concat%nxseg
                  iblk = hist_concat%xsegs(ixseg)%blk
                  jblk = hist_concat%ysegs(iyseg)%blk
                  IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                     xbdsp = hist_concat%xsegs(ixseg)%bdsp
                     ybdsp = hist_concat%ysegs(iyseg)%bdsp
                     xgdsp = hist_concat%xsegs(ixseg)%gdsp
                     ygdsp = hist_concat%ysegs(iyseg)%gdsp
                     xcnt = hist_concat%xsegs(ixseg)%cnt
                     ycnt = hist_concat%ysegs(iyseg)%cnt

                     vdata (xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt) = &
                        wdata%blk(iblk,jblk)%val(xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)
                  ENDIF
               ENDDO
            ENDDO
#endif

#ifdef USEMPI
            IF (.not. DEF_HIST_WriteBack) THEN
#endif
               IF (itime >= 1) THEN

                  CALL ncio_write_serial_time (filename, dataname, itime, vdata, &
                     'lon', 'lat', 'time', compress)

               ELSEIF (itime == -1) THEN
                  CALL ncio_write_serial (filename, dataname, vdata, 'lon', 'lat', compress)
               ENDIF

               IF (itime <= 1) THEN
                  CALL ncio_put_attr (filename, dataname, 'long_name', longname)
                  CALL ncio_put_attr (filename, dataname, 'units', units)
                  CALL ncio_put_attr (filename, dataname, 'missing_value', spval)
               ENDIF

               deallocate (vdata)
#ifdef USEMPI
            ENDIF
#endif

         ENDIF

#ifdef USEMPI
         IF (p_is_io) THEN
            DO iyseg = 1, hist_concat%nyseg
               DO ixseg = 1, hist_concat%nxseg

                  iblk = hist_concat%xsegs(ixseg)%blk
                  jblk = hist_concat%ysegs(iyseg)%blk

                  IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN

                     xbdsp = hist_concat%xsegs(ixseg)%bdsp
                     ybdsp = hist_concat%ysegs(iyseg)%bdsp
                     xcnt = hist_concat%xsegs(ixseg)%cnt
                     ycnt = hist_concat%ysegs(iyseg)%cnt

                     allocate (sbuf (xcnt,ycnt))
                     sbuf = wdata%blk(iblk,jblk)%val(xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)

                     IF (.not. DEF_HIST_WriteBack) THEN
                        smesg = (/p_iam_glb, ixseg, iyseg/)
                        CALL mpi_send (smesg, 3, MPI_INTEGER, &
                           p_address_master, hist_data_id, p_comm_glb, p_err)
                        CALL mpi_send (sbuf, xcnt*ycnt, MPI_REAL8, &
                           p_address_master, hist_data_id, p_comm_glb, p_err)
                     ELSE
                        CALL hist_writeback_var (hist_data_id, ixseg, iyseg, wdata2d = sbuf)
                     ENDIF

                     deallocate (sbuf)

                  ENDIF
               ENDDO
            ENDDO
         ENDIF
#endif

         hist_data_id = mod(hist_data_id-10000,10000) + 100001

      ELSEIF (trim(DEF_HIST_mode) == 'block') THEN

         IF (p_is_io) THEN

            DO iblkme = 1, gblock%nblkme
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)

               IF ((grid%xcnt(iblk) == 0) .or. (grid%ycnt(jblk) == 0)) CYCLE

               CALL get_filename_block (filename, iblk, jblk, fileblock)

               IF (itime >= 1) THEN

                  CALL ncio_write_serial_time (fileblock, dataname, itime, &
                     wdata%blk(iblk,jblk)%val, 'lon', 'lat', 'time', compress)

               ELSEIF (itime == -1) THEN
                  CALL ncio_write_serial (fileblock, dataname, &
                     wdata%blk(iblk,jblk)%val, 'lon', 'lat', compress)
               ENDIF

            ENDDO

         ENDIF
      ENDIF

   END SUBROUTINE hist_write_var_real8_2d


   SUBROUTINE hist_write_var_real8_3d ( &
         filename, dataname, dim1name, grid, itime, wdata, compress, longname, units)

   USE MOD_Block
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   character (len=*), intent(in) :: filename
   character (len=*), intent(in) :: dataname
   character (len=*), intent(in) :: dim1name
   type (grid_type),  intent(in) :: grid
   integer, intent(in) :: itime

   type (block_data_real8_3d), intent(in) :: wdata

   integer, intent(in) :: compress
   character(len=*), intent(in) :: longname
   character(len=*), intent(in) :: units

   ! Local variables
   integer :: iblkme, iblk, jblk, idata, ixseg, iyseg
   integer :: xcnt, ycnt, ndim1, xbdsp, ybdsp, xgdsp, ygdsp, idim1
   integer :: rmesg(4), smesg(4), isrc
   character(len=256) :: fileblock
   real(r8), allocatable :: rbuf(:,:,:), sbuf(:,:,:), vdata(:,:,:)

      IF (trim(DEF_HIST_mode) == 'one') THEN

         IF (p_is_master) THEN

#ifdef USEMPI
            IF (.not. DEF_HIST_WriteBack) THEN

               DO idata = 1, hist_concat%ndatablk

                  CALL mpi_recv (rmesg, 4, MPI_INTEGER, MPI_ANY_SOURCE, &
                     hist_data_id, p_comm_glb, p_stat, p_err)

                  isrc  = rmesg(1)
                  ixseg = rmesg(2)
                  iyseg = rmesg(3)
                  ndim1 = rmesg(4)

                  xgdsp = hist_concat%xsegs(ixseg)%gdsp
                  ygdsp = hist_concat%ysegs(iyseg)%gdsp
                  xcnt = hist_concat%xsegs(ixseg)%cnt
                  ycnt = hist_concat%ysegs(iyseg)%cnt

                  allocate (rbuf (ndim1,xcnt,ycnt))

                  CALL mpi_recv (rbuf, ndim1*xcnt*ycnt, MPI_REAL8, &
                     isrc, hist_data_id, p_comm_glb, p_stat, p_err)

                  IF (idata == 1) THEN
                     allocate (vdata (hist_concat%ginfo%nlon, hist_concat%ginfo%nlat, ndim1))
                     vdata(:,:,:) = spval
                  ENDIF

                  DO idim1 = 1, ndim1
                     vdata (xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt, idim1) = rbuf(idim1,:,:)
                  ENDDO

                  deallocate (rbuf)
               ENDDO

            ELSE
               CALL hist_writeback_var_header (hist_data_id, filename, dataname, &
                  4, 'lon', 'lat', dim1name, 'time', '', compress, longname, units)
            ENDIF
#else
            ndim1 = wdata%ub1 - wdata%lb1 + 1
            allocate (vdata (hist_concat%ginfo%nlon, hist_concat%ginfo%nlat, ndim1))
            vdata(:,:,:) = spval

            DO iyseg = 1, hist_concat%nyseg
               DO ixseg = 1, hist_concat%nxseg
                  iblk = hist_concat%xsegs(ixseg)%blk
                  jblk = hist_concat%ysegs(iyseg)%blk
                  IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                     xbdsp = hist_concat%xsegs(ixseg)%bdsp
                     ybdsp = hist_concat%ysegs(iyseg)%bdsp
                     xgdsp = hist_concat%xsegs(ixseg)%gdsp
                     ygdsp = hist_concat%ysegs(iyseg)%gdsp
                     xcnt = hist_concat%xsegs(ixseg)%cnt
                     ycnt = hist_concat%ysegs(iyseg)%cnt

                     DO idim1 = 1, ndim1
                        vdata (xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt,idim1) = &
                           wdata%blk(iblk,jblk)%val(wdata%lb1+idim1-1,xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
#endif


#ifdef USEMPI
            IF (.not. DEF_HIST_WriteBack) THEN
#endif
               CALL ncio_define_dimension (filename, dim1name, ndim1)

               CALL ncio_write_serial_time (filename, dataname, itime, &
                  vdata, 'lon', 'lat', dim1name, 'time', compress)

               IF (itime == 1) THEN
                  CALL ncio_put_attr (filename, dataname, 'long_name', longname)
                  CALL ncio_put_attr (filename, dataname, 'units', units)
                  CALL ncio_put_attr (filename, dataname, 'missing_value', spval)
               ENDIF

               deallocate (vdata)
#ifdef USEMPI
            ENDIF
#endif
         ENDIF

#ifdef USEMPI
         IF (p_is_io) THEN

            DO iyseg = 1, hist_concat%nyseg
               DO ixseg = 1, hist_concat%nxseg

                  iblk = hist_concat%xsegs(ixseg)%blk
                  jblk = hist_concat%ysegs(iyseg)%blk

                  IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN

                     xbdsp = hist_concat%xsegs(ixseg)%bdsp
                     ybdsp = hist_concat%ysegs(iyseg)%bdsp
                     xcnt = hist_concat%xsegs(ixseg)%cnt
                     ycnt = hist_concat%ysegs(iyseg)%cnt
                     ndim1 = size(wdata%blk(iblk,jblk)%val,1)

                     allocate (sbuf (ndim1,xcnt,ycnt))
                     sbuf = wdata%blk(iblk,jblk)%val(:, xbdsp+1:xbdsp+xcnt, ybdsp+1:ybdsp+ycnt)

                     IF (.not. DEF_HIST_WriteBack) THEN
                        smesg = (/p_iam_glb, ixseg, iyseg, ndim1/)
                        CALL mpi_send (smesg, 4, MPI_INTEGER, &
                           p_address_master, hist_data_id, p_comm_glb, p_err)
                        CALL mpi_send (sbuf, ndim1*xcnt*ycnt, MPI_REAL8, &
                           p_address_master, hist_data_id, p_comm_glb, p_err)
                     ELSE
                        CALL hist_writeback_var (hist_data_id, ixseg, iyseg, wdata3d = sbuf)
                     ENDIF

                     deallocate (sbuf)
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
#endif

         hist_data_id = mod(hist_data_id-10000,10000) + 100001

      ELSEIF (trim(DEF_HIST_mode) == 'block') THEN

         IF (p_is_io) THEN

            DO iblkme = 1, gblock%nblkme
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)

               IF ((grid%xcnt(iblk) == 0) .or. (grid%ycnt(jblk) == 0)) CYCLE

               CALL get_filename_block (filename, iblk, jblk, fileblock)

               CALL ncio_define_dimension (fileblock, dim1name, wdata%ub1-wdata%lb1+1)

               CALL ncio_write_serial_time (fileblock, dataname, itime, &
                  wdata%blk(iblk,jblk)%val, dim1name, 'lon', 'lat', 'time', compress)

            ENDDO

         ENDIF
      ENDIF

   END SUBROUTINE hist_write_var_real8_3d


   SUBROUTINE hist_write_var_real8_4d ( &
         filename, dataname, dim1name, dim2name, grid, itime, wdata, compress, longname, units)

   USE MOD_Block
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   character (len=*), intent(in) :: filename
   character (len=*), intent(in) :: dataname
   character (len=*), intent(in) :: dim1name, dim2name
   type (grid_type),  intent(in) :: grid
   integer, intent(in) :: itime

   type (block_data_real8_4d), intent(in) :: wdata

   integer, intent(in) :: compress
   character(len=*), intent(in) :: longname
   character(len=*), intent(in) :: units

   ! Local variables
   integer :: iblkme, iblk, jblk, idata, ixseg, iyseg
   integer :: xcnt, ycnt, ndim1, ndim2, xbdsp, ybdsp, xgdsp, ygdsp, idim1, idim2
   integer :: rmesg(5), smesg(5), isrc
   character(len=256) :: fileblock
   real(r8), allocatable :: rbuf(:,:,:,:), sbuf(:,:,:,:), vdata(:,:,:,:)

      IF (trim(DEF_HIST_mode) == 'one') THEN

         IF (p_is_master) THEN

#ifdef USEMPI
            IF (.not. DEF_HIST_WriteBack) THEN

               DO idata = 1, hist_concat%ndatablk

                  CALL mpi_recv (rmesg, 5, MPI_INTEGER, MPI_ANY_SOURCE, &
                     hist_data_id, p_comm_glb, p_stat, p_err)

                  isrc  = rmesg(1)
                  ixseg = rmesg(2)
                  iyseg = rmesg(3)
                  ndim1 = rmesg(4)
                  ndim2 = rmesg(5)

                  xgdsp = hist_concat%xsegs(ixseg)%gdsp
                  ygdsp = hist_concat%ysegs(iyseg)%gdsp
                  xcnt = hist_concat%xsegs(ixseg)%cnt
                  ycnt = hist_concat%ysegs(iyseg)%cnt

                  allocate (rbuf (ndim1,ndim2,xcnt,ycnt))

                  CALL mpi_recv (rbuf, ndim1*ndim2*xcnt*ycnt, MPI_REAL8, &
                     isrc, hist_data_id, p_comm_glb, p_stat, p_err)

                  IF (idata == 1) THEN
                     allocate (vdata (hist_concat%ginfo%nlon,hist_concat%ginfo%nlat,ndim1,ndim2))
                     vdata(:,:,:,:) = spval
                  ENDIF

                  DO idim1 = 1, ndim1
                     DO idim2 = 1, ndim2
                        vdata (xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt, idim1, idim2) = rbuf(idim1,idim2,:,:)
                     ENDDO
                  ENDDO

                  deallocate (rbuf)
               ENDDO

            ELSE
               CALL hist_writeback_var_header (hist_data_id, filename, dataname, &
                  5, 'lon', 'lat', dim1name, dim2name, 'time', compress, longname, units)
            ENDIF
#else
            ndim1 = wdata%ub1 - wdata%lb1 + 1
            ndim2 = wdata%ub2 - wdata%lb2 + 1
            allocate (vdata (hist_concat%ginfo%nlon,hist_concat%ginfo%nlat,ndim1,ndim2))
            vdata(:,:,:,:) = spval

            DO iyseg = 1, hist_concat%nyseg
               DO ixseg = 1, hist_concat%nxseg
                  iblk = hist_concat%xsegs(ixseg)%blk
                  jblk = hist_concat%ysegs(iyseg)%blk
                  IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                     xbdsp = hist_concat%xsegs(ixseg)%bdsp
                     ybdsp = hist_concat%ysegs(iyseg)%bdsp
                     xgdsp = hist_concat%xsegs(ixseg)%gdsp
                     ygdsp = hist_concat%ysegs(iyseg)%gdsp
                     xcnt = hist_concat%xsegs(ixseg)%cnt
                     ycnt = hist_concat%ysegs(iyseg)%cnt

                     DO idim1 = 1, ndim1
                        DO idim2 = 1, ndim2
                           vdata (xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt, idim1, idim2) = &
                              wdata%blk(iblk,jblk)%val(idim1,idim2,xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO

#endif

#ifdef USEMPI
            IF (.not. DEF_HIST_WriteBack) THEN
#endif
               CALL ncio_define_dimension (filename, dim1name, ndim1)
               CALL ncio_define_dimension (filename, dim2name, ndim2)

               CALL ncio_write_serial_time (filename, dataname, itime, vdata, &
                  'lon', 'lat', dim1name, dim2name, 'time', compress)

               IF (itime == 1) THEN
                  CALL ncio_put_attr (filename, dataname, 'long_name', longname)
                  CALL ncio_put_attr (filename, dataname, 'units', units)
                  CALL ncio_put_attr (filename, dataname, 'missing_value', spval)
               ENDIF

               deallocate (vdata)
#ifdef USEMPI
            ENDIF
#endif
         ENDIF

#ifdef USEMPI
         IF (p_is_io) THEN

            DO iyseg = 1, hist_concat%nyseg
               DO ixseg = 1, hist_concat%nxseg

                  iblk = hist_concat%xsegs(ixseg)%blk
                  jblk = hist_concat%ysegs(iyseg)%blk

                  IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN

                     xbdsp = hist_concat%xsegs(ixseg)%bdsp
                     ybdsp = hist_concat%ysegs(iyseg)%bdsp
                     xcnt = hist_concat%xsegs(ixseg)%cnt
                     ycnt = hist_concat%ysegs(iyseg)%cnt

                     ndim1 = size(wdata%blk(iblk,jblk)%val,1)
                     ndim2 = size(wdata%blk(iblk,jblk)%val,2)
                     allocate (sbuf (ndim1,ndim2,xcnt,ycnt))
                     sbuf = wdata%blk(iblk,jblk)%val(:,:,xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)

                     IF (.not. DEF_HIST_WriteBack) THEN
                        smesg = (/p_iam_glb, ixseg, iyseg, ndim1, ndim2/)
                        CALL mpi_send (smesg, 5, MPI_INTEGER, &
                           p_address_master, hist_data_id, p_comm_glb, p_err)
                        CALL mpi_send (sbuf, ndim1*ndim2*xcnt*ycnt, MPI_REAL8, &
                           p_address_master, hist_data_id, p_comm_glb, p_err)
                     ELSE
                        CALL hist_writeback_var (hist_data_id, ixseg, iyseg, wdata4d = sbuf)
                     ENDIF

                     deallocate (sbuf)
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
#endif

         hist_data_id = mod(hist_data_id-10000,10000) + 100001

      ELSEIF (trim(DEF_HIST_mode) == 'block') THEN
         IF (p_is_io) THEN

            DO iblkme = 1, gblock%nblkme
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)

               IF ((grid%xcnt(iblk) == 0) .or. (grid%ycnt(jblk) == 0)) CYCLE

               CALL get_filename_block (filename, iblk, jblk, fileblock)

               CALL ncio_define_dimension (fileblock, dim1name, wdata%ub1-wdata%lb1+1)
               CALL ncio_define_dimension (fileblock, dim2name, wdata%ub2-wdata%lb2+1)

               CALL ncio_write_serial_time (fileblock, dataname, itime, &
                  wdata%blk(iblk,jblk)%val, dim1name, dim2name, 'lon', 'lat', 'time', compress)

            ENDDO

         ENDIF
      ENDIF

   END SUBROUTINE hist_write_var_real8_4d


   SUBROUTINE hist_write_grid_info (fileblock, grid, iblk, jblk)

   USE MOD_Block
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

      CALL ncio_write_colm_dimension (fileblock)

   END SUBROUTINE hist_write_grid_info

END MODULE MOD_HistGridded
