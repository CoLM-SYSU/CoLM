#include <define.h>

MODULE MOD_CheckEquilibrium

!----------------------------------------------------------------------------
! !DESCRIPTION:
!
!   Check equilibrium state.
!
!  Created by Shupeng Zhang, 10/2024
!----------------------------------------------------------------------------

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Grid
   USE netcdf
   USE MOD_NetCDFSerial
   USE MOD_SpatialMapping
   USE MOD_Forcing,     only: forcmask_pch
   USE MOD_Vars_Global, only: spval

   ! ----- Variables -----
   integer :: nyearcheck

   integer :: timestrlen
   character(len=24) :: timeform

   real(r8), allocatable :: tws_last  (:)
   real(r8), allocatable :: tws_this  (:)

   real(r8), allocatable :: prcp_year (:)
   real(r8), allocatable :: et_year   (:)
   real(r8), allocatable :: rnof_year (:)
   real(r8), allocatable :: rchg_year (:)

   real(r8), allocatable :: patcharea (:)  ! m^2

   character(len=256) :: mesg_equilibrium

#ifndef SinglePoint
   type(grid_type)            :: gridcheck
   type(grid_concat_type)     :: gcheck_concat
   type(spatial_mapping_type) :: map_check

   integer :: check_data_id = 0
#endif

   PUBLIC :: CheckEqb_init
   PUBLIC :: CheckEquilibrium
   PUBLIC :: CheckEqb_final

CONTAINS

   !-----------------------------------------------------------------------
   SUBROUTINE CheckEqb_init (n_spinupcycle, lc_year)

   USE MOD_Utils
   USE MOD_Vars_Global,        only: nl_soil
   USE MOD_Forcing,            only: gforc
   USE MOD_LandPatch,          only: numpatch, landpatch
   USE MOD_Pixel,              only: pixel
   USE MOD_Mesh,               only: mesh
   USE MOD_Vars_TimeVariables, only: wdsrf, ldew, scv, wetwat, &
                                     wliq_soisno, wice_soisno, wa
   IMPLICIT NONE

   integer, intent(in) :: n_spinupcycle
   integer, intent(in) :: lc_year

   ! Local Variable
   integer :: ilev, ip, ie, ipxl
   character(len=256) :: filename, cyear
   real(r8) :: totaldtws, totalprcp, pct_dtws_prcp


      IF (.not. DEF_CheckEquilibrium) RETURN

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN

            allocate (tws_last  (numpatch));   tws_last (:) = spval;
            allocate (tws_this  (numpatch));   tws_this (:) = spval;
            allocate (prcp_year (numpatch));   prcp_year(:) = spval;
            allocate (et_year   (numpatch));   et_year  (:) = spval;
            allocate (rnof_year (numpatch));   rnof_year(:) = spval;
            allocate (rchg_year (numpatch));   rchg_year(:) = spval;

         ENDIF
      ENDIF

      IF (n_spinupcycle >= 10000) THEN
         timestrlen = 16
         timeform = "('spinup',I5.5,'-',I4.4)"
      ELSEIF (n_spinupcycle >= 1000) THEN
         timestrlen = 15
         timeform = "('spinup',I4.4,'-',I4.4)"
      ELSEIF (n_spinupcycle >= 100) THEN
         timestrlen = 14
         timeform = "('spinup',I3.3,'-',I4.4)"
      ELSEIF (n_spinupcycle >= 10) THEN
         timestrlen = 13
         timeform = "('spinup',I2.2,'-',I4.4)"
      ELSE
         timestrlen = 12
         timeform = "('spinup',I1.1,'-',I4.4)"
      ENDIF

      mesg_equilibrium = ''

#ifndef SinglePoint
      ! grid
#ifdef GRIDBASED
      write(cyear,'(i4.4)') lc_year
      filename = trim(DEF_dir_landdata) // '/mesh/' //trim(cyear) // '/mesh.nc'
      CALL gridcheck%define_from_file (filename)
#else
      CALL gridcheck%define_by_copy (gforc)
#endif
      ! grid info for output
      CALL gcheck_concat%set (gridcheck)
      ! mapping from patch to grid
      CALL map_check%build_arealweighted (gridcheck, landpatch)
#endif

      IF ((p_is_worker) .and. (numpatch > 0)) THEN
         tws_last = wdsrf                                 ! 1. surface water
         CALL add_spv (ldew, tws_last)                    ! 2. water on foliage
         CALL add_spv (scv , tws_last)                    ! 3. snow cover water equivalent
         IF (DEF_USE_VariablySaturatedFlow) THEN
            CALL add_spv (wetwat, tws_last)               ! 4. water in wetland
         ENDIF
         DO ilev = 1, nl_soil
            CALL add_spv (wliq_soisno(ilev,:), tws_last)  ! 5. liquid water in soil
            CALL add_spv (wice_soisno(ilev,:), tws_last)  ! 6. ice in soil
         ENDDO
         CALL add_spv (wa, tws_last)                      ! 7. water in aquifer
      ENDIF

      IF (p_is_worker) THEN

         IF (numpatch > 0) THEN
            allocate (patcharea (numpatch))
            patcharea(:) = 0.
            DO ip = 1, numpatch
               ie = landpatch%ielm(ip)
               DO ipxl = landpatch%ipxstt(ip), landpatch%ipxend(ip)
                  patcharea(ip) = patcharea(ip) + 1.0e6 * areaquad ( &
                     pixel%lat_s(mesh(ie)%ilat(ipxl)), pixel%lat_n(mesh(ie)%ilat(ipxl)), &
                     pixel%lon_w(mesh(ie)%ilon(ipxl)), pixel%lon_e(mesh(ie)%ilon(ipxl)) )
               ENDDO
               IF (landpatch%has_shared) THEN
                  patcharea(ip) = patcharea(ip) * landpatch%pctshared(ip)
               ENDIF
            ENDDO
         ENDIF

      ENDIF

      nyearcheck = 0

   END SUBROUTINE CheckEqb_init

   !-----------------------------------------------------------------------
   SUBROUTINE CheckEqb_final ()

   IMPLICIT NONE

      IF (.not. DEF_CheckEquilibrium) RETURN

      IF (allocated(tws_last )) deallocate(tws_last )
      IF (allocated(tws_this )) deallocate(tws_this )
      IF (allocated(prcp_year)) deallocate(prcp_year)
      IF (allocated(et_year  )) deallocate(et_year  )
      IF (allocated(rnof_year)) deallocate(rnof_year)
      IF (allocated(rchg_year)) deallocate(rchg_year)
      IF (allocated(patcharea)) deallocate(patcharea)

   END SUBROUTINE CheckEqb_final

   !-----------------------------------------------------------------------
   SUBROUTINE CheckEquilibrium ( &
         idate, deltim, i_spinupcycle, is_spinup, dir_out, casename)

   USE MOD_Precision
   USE MOD_TimeManager
   USE MOD_DataType
   USE MOD_LandPatch,           only: numpatch
   USE MOD_Vars_Global,         only: nl_soil
   USE MOD_Vars_1DForcing,      only: forc_prc, forc_prl
   USE MOD_Vars_1DFluxes,       only: fevpa, rnof, rsur
   USE MOD_Vars_TimeInvariants, only: patchtype, patchmask
   USE MOD_Vars_TimeVariables,  only: wdsrf, ldew, scv, wetwat, wliq_soisno, &
                                      wice_soisno, wa, zwt

   IMPLICIT NONE

   integer,  intent(in) :: idate(3)
   real(r8), intent(in) :: deltim
   integer,  intent(in) :: i_spinupcycle
   logical,  intent(in) :: is_spinup

   character(len=*), intent(in) :: dir_out
   character(len=*), intent(in) :: casename

   ! Local variables
   logical  :: docheck
   integer  :: ilev
   character(len=256) :: filename, timestr
   integer  :: ncid, time_id, str_id, varid
   real(r8) :: totaldtws, totalprcp, pct_dtws_prcp

   real(r8), allocatable     :: rchg   (:)
   real(r8), allocatable     :: dtws   (:)
   logical,  allocatable     :: filter (:)
   real(r8), allocatable     :: vecone (:)
   type(block_data_real8_2d) :: sumarea


      IF (.not. DEF_CheckEquilibrium) RETURN

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            CALL add_spv (forc_prc, prcp_year, deltim)
            CALL add_spv (forc_prl, prcp_year, deltim)
            CALL add_spv (fevpa,    et_year,   deltim)
            CALL add_spv (rnof,     rnof_year, deltim)

            allocate (rchg (numpatch))
            rchg = spval
            WHERE ((forc_prc /= spval) .and. (forc_prl /= spval) .and. (fevpa /= spval) .and. (rsur /= spval))
               rchg = forc_prc + forc_prl - fevpa - rsur
            END WHERE

            CALL add_spv (rchg, rchg_year, deltim)

            deallocate (rchg)
         ENDIF
      ENDIF

      docheck = isendofyear (idate, deltim)

      IF (docheck) THEN

         IF ((p_is_worker) .and. (numpatch > 0)) THEN
            tws_this = wdsrf                                 ! 1. surface water
            CALL add_spv (ldew, tws_this)                    ! 2. water on foliage
            CALL add_spv (scv , tws_this)                    ! 3. snow cover water equivalent
            IF (DEF_USE_VariablySaturatedFlow) THEN
               CALL add_spv (wetwat, tws_this)               ! 4. water in wetland
            ENDIF
            DO ilev = 1, nl_soil
               CALL add_spv (wliq_soisno(ilev,:), tws_this)  ! 5. liquid water in soil
               CALL add_spv (wice_soisno(ilev,:), tws_this)  ! 6. ice in soil
            ENDDO
            CALL add_spv (wa, tws_this)                      ! 7. water in aquifer
         ENDIF

         nyearcheck = nyearcheck + 1

         IF (nyearcheck >= 1) THEN

            IF (p_is_worker) THEN

               totaldtws = 0.
               totalprcp = 0.

               IF (numpatch > 0) THEN

                  allocate (filter (numpatch))
                  filter = patchmask
                  IF (DEF_forcing%has_missing_value) THEN
                     filter = filter .and. forcmask_pch
                  ENDIF
#ifdef CatchLateralFlow
                  filter = filter .and. (patchtype <= 4)
#else
                  filter = filter .and. (patchtype <= 2)
#endif

                  allocate (dtws (numpatch))
                  WHERE (filter)
                     dtws = tws_this - tws_last
                  ELSEWHERE
                     dtws = spval
                  END WHERE

                  IF (any(filter)) THEN
                     totaldtws = sum(dtws*patcharea, mask=filter)
                     totalprcp = sum(prcp_year*patcharea, mask=filter)
                  ENDIF

                  allocate (vecone (numpatch))
                  vecone(:) = 1.

               ENDIF

#ifdef USEMPI
               CALL mpi_allreduce (MPI_IN_PLACE, totaldtws, 1, MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
               CALL mpi_allreduce (MPI_IN_PLACE, totalprcp, 1, MPI_REAL8, MPI_SUM, p_comm_worker, p_err)

               IF (totalprcp > 0.) THEN
                  pct_dtws_prcp = totaldtws/totalprcp
               ELSE
                  pct_dtws_prcp = spval
               ENDIF

               IF (p_iam_worker == p_root) THEN
                  CALL mpi_send (pct_dtws_prcp, 1, MPI_REAL8, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
               ENDIF
#endif
            ENDIF

            IF (p_is_master) THEN

               IF (is_spinup) THEN
                  write(timestr,timeform) i_spinupcycle, idate(1)
               ELSE
                  write(timestr,'(I4.4)') idate(1)
               ENDIF

#ifdef USEMPI
               CALL mpi_recv (pct_dtws_prcp, 1, MPI_REAL8, p_address_worker(p_root), &
                  mpi_tag_mesg, p_comm_glb, p_stat, p_err)
#endif
               IF (pct_dtws_prcp /= spval) THEN
                  write(mesg_equilibrium,'(A,F0.2,3A)') 'Total delTWS/precipitation is ', &
                     pct_dtws_prcp*100., '% in ', trim(timestr), '. Check history for detail.'
               ELSE
                  write(mesg_equilibrium,'(3A)') 'Precipitation is 0 ', trim(timestr), '.'
               ENDIF

               filename = trim(dir_out) // '/' // trim(casename) //'_check_equilibrium.nc'

               IF (nyearcheck == 1) THEN

                  CALL ncio_create_file (trim(filename))

                  CALL ncio_define_dimension(filename, 'iyear',   0)
                  CALL ncio_define_dimension(filename, 'timestr', timestrlen)

                  CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
                  CALL nccheck( nf90_inq_dimid(ncid, 'iyear',   time_id) )
                  CALL nccheck( nf90_inq_dimid(ncid, 'timestr', str_id ) )
                  CALL nccheck( nf90_redef(ncid) )
                  CALL nccheck( nf90_def_var(ncid, 'iyear', NF90_CHAR, (/str_id,time_id/), varid) )
                  CALL nccheck( nf90_put_att(ncid, varid, 'long_name', 'iyear in all spinup cycles') )
                  CALL nccheck( nf90_enddef(ncid) )
                  CALL nccheck( nf90_close(ncid) )

#ifndef SinglePoint
                  CALL ncio_define_dimension(filename, 'lat' , gcheck_concat%ginfo%nlat)
                  CALL ncio_define_dimension(filename, 'lon' , gcheck_concat%ginfo%nlon)

                  CALL ncio_write_serial (filename, 'lat', gcheck_concat%ginfo%lat_c, 'lat')
                  CALL ncio_put_attr (filename, 'lat', 'long_name', 'latitude')
                  CALL ncio_put_attr (filename, 'lat', 'units', 'degrees_north')

                  CALL ncio_write_serial (filename, 'lon', gcheck_concat%ginfo%lon_c, 'lon')
                  CALL ncio_put_attr (filename, 'lon', 'long_name', 'longitude')
                  CALL ncio_put_attr (filename, 'lon', 'units', 'degrees_east')
#else
                  CALL ncio_define_dimension(filename, 'patch', numpatch)
#endif

               ENDIF

               CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
               CALL nccheck( nf90_inq_varid(ncid, 'iyear', varid) )
               CALL nccheck( nf90_put_var(ncid, varid, timestr(1:timestrlen), &
                  (/1,nyearcheck/), (/timestrlen,1/)) )
               CALL nccheck( nf90_close(ncid) )

            ENDIF

#ifndef SinglePoint
            IF (p_is_io) CALL allocate_block_data (gridcheck, sumarea)

            CALL map_check%get_sumarea (sumarea, filter)

            IF (nyearcheck == 1) THEN
               CALL map_and_write_check_var ( &
                  vecone, filename, 'landarea', -1, sumarea, filter, &
                  'area of land excluding water bodies and glaciers in grid', 'km^2', &
                  amount_in_grid = .true.)
            ENDIF

            CALL map_and_write_check_var ( &
               dtws, filename, 'tws_change', nyearcheck, sumarea, filter, &
               'Change in terrestrial water storage', 'mm')
            CALL map_and_write_check_var ( &
               prcp_year, filename, 'total_precipitation', nyearcheck, sumarea, filter, &
               'total precipitation in a year', 'mm')
            CALL map_and_write_check_var ( &
               et_year, filename, 'total_evapotranspiration', nyearcheck, sumarea, filter, &
               'total evapotranspiration in a year', 'mm')
            CALL map_and_write_check_var ( &
               rnof_year, filename, 'total_runoff', nyearcheck, sumarea, filter, &
               'total runoff in a year', 'mm')
            CALL map_and_write_check_var ( &
               rchg_year, filename, 'total_recharge', nyearcheck, sumarea, filter, &
               'total recharge to ground water in a year', 'mm')
            CALL map_and_write_check_var ( &
               zwt, filename, 'zwt', nyearcheck, sumarea, filter, &
               'depth to water table', 'm')
#else
            CALL ncio_write_serial_time (filename, 'tws_change', &
                                         nyearcheck, dtws, 'patch', 'iyear')
            IF (nyearcheck == 1) THEN
               CALL ncio_put_attr (filename, 'tws_change', 'long_name', &
                  'Change in terrestrial water storage')
               CALL ncio_put_attr (filename, 'tws_change', 'units', 'mm')
               CALL ncio_put_attr (filename, 'tws_change', 'missing_value', spval)
            ENDIF

            CALL ncio_write_serial_time (filename, 'total_precipitation', &
                                         nyearcheck, prcp_year, 'patch', 'iyear')
            IF (nyearcheck == 1) THEN
               CALL ncio_put_attr (filename, 'total_precipitation', 'long_name', &
                  'total precipitation in a year')
               CALL ncio_put_attr (filename, 'total_precipitation', 'units', 'mm')
               CALL ncio_put_attr (filename, 'total_precipitation', 'missing_value', spval)
            ENDIF

            CALL ncio_write_serial_time (filename, 'total_evapotranspiration', &
                                         nyearcheck, et_year, 'patch', 'iyear')
            IF (nyearcheck == 1) THEN
               CALL ncio_put_attr (filename, 'total_evapotranspiration', 'long_name', &
                  'total evapotranspiration in a year')
               CALL ncio_put_attr (filename, 'total_evapotranspiration', 'units', 'mm')
               CALL ncio_put_attr (filename, 'total_evapotranspiration', 'missing_value', spval)
            ENDIF

            CALL ncio_write_serial_time (filename, 'total_runoff', &
                                         nyearcheck, rnof_year, 'patch', 'iyear')
            IF (nyearcheck == 1) THEN
               CALL ncio_put_attr (filename, 'total_runoff', 'long_name', &
                  'total runoff in a year')
               CALL ncio_put_attr (filename, 'total_runoff', 'units', 'mm')
               CALL ncio_put_attr (filename, 'total_runoff', 'missing_value', spval)
            ENDIF

            CALL ncio_write_serial_time (filename, 'total_recharge', &
                                         nyearcheck, rchg_year, 'patch', 'iyear')
            IF (nyearcheck == 1) THEN
               CALL ncio_put_attr (filename, 'total_recharge', 'long_name', &
                  'total recharge to ground water in a year')
               CALL ncio_put_attr (filename, 'total_recharge', 'units', 'mm')
               CALL ncio_put_attr (filename, 'total_recharge', 'missing_value', spval)
            ENDIF

            CALL ncio_write_serial_time (filename, 'zwt', &
                                         nyearcheck, zwt, 'patch', 'iyear')
            IF (nyearcheck == 1) THEN
               CALL ncio_put_attr (filename, 'zwt', 'long_name', 'depth to water table')
               CALL ncio_put_attr (filename, 'zwt', 'units', 'm')
               CALL ncio_put_attr (filename, 'zwt', 'missing_value', spval)
            ENDIF
#endif
         ENDIF

         IF ((p_is_worker) .and. (numpatch > 0)) THEN
            prcp_year(:) = spval
            et_year  (:) = spval
            rnof_year(:) = spval
            rchg_year(:) = spval
            tws_last = tws_this
         ENDIF

         IF (allocated(dtws  )) deallocate(dtws  )
         IF (allocated(filter)) deallocate(filter)
         IF (allocated(vecone)) deallocate(vecone)

      ENDIF

   END SUBROUTINE CheckEquilibrium

   !-----------------------------------------------------------------------
#ifndef SinglePoint
   SUBROUTINE map_and_write_check_var ( &
         vector, filename, varname, itime_in_file, sumarea, filter, &
         longname, units, amount_in_grid)

   USE MOD_Block
   IMPLICIT NONE

   real(r8), intent(in) :: vector(:)

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: varname
   integer,          intent(in) :: itime_in_file
   character(len=*), intent(in) :: longname
   character(len=*), intent(in) :: units

   logical, intent(in), optional :: amount_in_grid

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)

   ! Local variables
   type(block_data_real8_2d) :: data_xy_2d
   integer :: xblk, yblk, xloc, yloc, xcnt, ycnt, xbdsp, ybdsp, xgdsp, ygdsp
   integer :: iblkme, iblk, jblk, idata, ixseg, iyseg
   integer :: rmesg(3), smesg(3), isrc
   real(r8), allocatable :: rbuf(:,:), sbuf(:,:), vdata(:,:)
   logical :: amount

      IF (p_is_io) CALL allocate_block_data (gridcheck, data_xy_2d)
      CALL map_check%pset2grid (vector, data_xy_2d, spv = spval, msk = filter)

      amount = .false.
      IF (present(amount_in_grid)) amount = amount_in_grid

      IF (.not. amount) THEN
         IF (p_is_io) THEN
            DO iblkme = 1, gblock%nblkme
               xblk = gblock%xblkme(iblkme)
               yblk = gblock%yblkme(iblkme)

               DO yloc = 1, gridcheck%ycnt(yblk)
                  DO xloc = 1, gridcheck%xcnt(xblk)

                     IF (sumarea%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) THEN
                        IF (data_xy_2d%blk(xblk,yblk)%val(xloc,yloc) /= spval) THEN
                           data_xy_2d%blk(xblk,yblk)%val(xloc,yloc) &
                              = data_xy_2d%blk(xblk,yblk)%val(xloc,yloc) &
                              / sumarea%blk(xblk,yblk)%val(xloc,yloc)
                        ENDIF
                     ELSE
                        data_xy_2d%blk(xblk,yblk)%val(xloc,yloc) = spval
                     ENDIF

                  ENDDO
               ENDDO

            ENDDO
         ENDIF
      ENDIF

      check_data_id = mod(check_data_id,100) + 1
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN

         allocate (vdata (gcheck_concat%ginfo%nlon, gcheck_concat%ginfo%nlat))
         vdata(:,:) = spval

#ifdef USEMPI
         DO idata = 1, gcheck_concat%ndatablk
            CALL mpi_recv (rmesg, 3, MPI_INTEGER, MPI_ANY_SOURCE, &
               check_data_id, p_comm_glb, p_stat, p_err)

            isrc  = rmesg(1)
            ixseg = rmesg(2)
            iyseg = rmesg(3)

            xgdsp = gcheck_concat%xsegs(ixseg)%gdsp
            ygdsp = gcheck_concat%ysegs(iyseg)%gdsp
            xcnt  = gcheck_concat%xsegs(ixseg)%cnt
            ycnt  = gcheck_concat%ysegs(iyseg)%cnt

            allocate (rbuf(xcnt,ycnt))

            CALL mpi_recv (rbuf, xcnt*ycnt, MPI_REAL8, &
               isrc, check_data_id, p_comm_glb, p_stat, p_err)

            vdata (xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt) = rbuf
            deallocate (rbuf)

         ENDDO
#else
         DO iyseg = 1, gcheck_concat%nyseg
            DO ixseg = 1, gcheck_concat%nxseg
               iblk = gcheck_concat%xsegs(ixseg)%blk
               jblk = gcheck_concat%ysegs(iyseg)%blk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  xbdsp = gcheck_concat%xsegs(ixseg)%bdsp
                  ybdsp = gcheck_concat%ysegs(iyseg)%bdsp
                  xgdsp = gcheck_concat%xsegs(ixseg)%gdsp
                  ygdsp = gcheck_concat%ysegs(iyseg)%gdsp
                  xcnt  = gcheck_concat%xsegs(ixseg)%cnt
                  ycnt  = gcheck_concat%ysegs(iyseg)%cnt

                  vdata (xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt) = &
                     data_xy_2d%blk(iblk,jblk)%val(xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)
               ENDIF
            ENDDO
         ENDDO
#endif

         IF (itime_in_file >= 1) THEN
            CALL ncio_write_serial_time (filename, varname, itime_in_file, vdata, &
               'lon', 'lat', 'iyear', compress = 1)
         ELSE
            CALL ncio_write_serial (filename, varname, vdata, 'lon', 'lat', compress = 1)
         ENDIF

         IF (itime_in_file <= 1) THEN
            CALL ncio_put_attr (filename, varname, 'long_name', longname)
            CALL ncio_put_attr (filename, varname, 'units', units)
            CALL ncio_put_attr (filename, varname, 'missing_value', spval)
         ENDIF

         deallocate (vdata)

      ENDIF

#ifdef USEMPI
      IF (p_is_io) THEN
         DO iyseg = 1, gcheck_concat%nyseg
            DO ixseg = 1, gcheck_concat%nxseg

               iblk = gcheck_concat%xsegs(ixseg)%blk
               jblk = gcheck_concat%ysegs(iyseg)%blk

               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN

                  xbdsp = gcheck_concat%xsegs(ixseg)%bdsp
                  ybdsp = gcheck_concat%ysegs(iyseg)%bdsp
                  xcnt  = gcheck_concat%xsegs(ixseg)%cnt
                  ycnt  = gcheck_concat%ysegs(iyseg)%cnt

                  allocate (sbuf (xcnt,ycnt))
                  sbuf = data_xy_2d%blk(iblk,jblk)%val(xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)

                  smesg = (/p_iam_glb, ixseg, iyseg/)
                  CALL mpi_send (smesg, 3, MPI_INTEGER, &
                     p_address_master, check_data_id, p_comm_glb, p_err)
                  CALL mpi_send (sbuf, xcnt*ycnt, MPI_REAL8, &
                     p_address_master, check_data_id, p_comm_glb, p_err)

                  deallocate (sbuf)

               ENDIF
            ENDDO
         ENDDO
      ENDIF
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE map_and_write_check_var
#endif

   !-----------------------------------------------------------------------
   SUBROUTINE add_spv (var, s, dt)

   USE MOD_Precision

   IMPLICIT NONE

   real(r8), intent(in)    :: var(:)
   real(r8), intent(inout) :: s  (:)
   real(r8), intent(in), optional :: dt
   ! Local variables
   integer :: i

      IF (present(dt)) THEN
         DO i = lbound(var,1), ubound(var,1)
            IF (var(i) /= spval) THEN
               IF (s(i) /= spval) THEN
                  s(i) = s(i) + var(i)*dt
               ELSE
                  s(i) = var(i)*dt
               ENDIF
            ENDIF
         ENDDO
      ELSE
         DO i = lbound(var,1), ubound(var,1)
            IF (var(i) /= spval) THEN
               IF (s(i) /= spval) THEN
                  s(i) = s(i) + var(i)
               ELSE
                  s(i) = var(i)
               ENDIF
            ENDIF
         ENDDO
      ENDIF

   END SUBROUTINE add_spv

END MODULE MOD_CheckEquilibrium
