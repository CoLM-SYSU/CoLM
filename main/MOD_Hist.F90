#include <define.h>

MODULE MOD_Hist

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

   USE MOD_Vars_1DAccFluxes
   USE MOD_Vars_Global, only: spval
   USE MOD_NetCDFSerial

   USE MOD_HistGridded
#if (defined UNSTRUCTURED || defined CATCHMENT)
   USE MOD_HistVector
#endif
#ifdef SinglePoint
   USE MOD_HistSingle
#endif
#ifdef CatchLateralFlow
   USE MOD_Catch_Hist
#endif
#ifdef GridRiverLakeFlow
   USE MOD_Grid_RiverLakeHist
#endif
#ifdef EXTERNAL_LAKE
   USE MOD_Lake_Hist
#endif

   PUBLIC :: hist_init
   PUBLIC :: hist_out
   PUBLIC :: hist_final

   character(len=10) :: HistForm ! 'Gridded', 'Vector', 'Single'

   character(len=256) :: file_last = 'null'

!--------------------------------------------------------------------------
CONTAINS

   SUBROUTINE hist_init (dir_hist, lulcc_call)

   IMPLICIT NONE

   character(len=*) , intent(in) :: dir_hist
   logical, optional, intent(in) :: lulcc_call

      CALL allocate_acc_fluxes ()
      CALL FLUSH_acc_fluxes ()

      HistForm = 'Gridded'
#if (defined UNSTRUCTURED || defined CATCHMENT)
      IF (DEF_HISTORY_IN_VECTOR) THEN
         HistForm = 'Vector'
      ENDIF
#endif
#ifdef SinglePoint
      HistForm = 'Single'
#endif

      IF (HistForm == 'Gridded') THEN
         IF (present(lulcc_call)) THEN
            CALL hist_gridded_init (dir_hist, lulcc_call)
         ELSE
            CALL hist_gridded_init (dir_hist)
         ENDIF
#ifdef SinglePoint
      ELSEIF (HistForm == 'Single') THEN
         CALL hist_single_init  ()
#endif
      ENDIF

#ifdef CatchLateralFlow
      CALL hist_basin_init ()
#endif

#ifdef GridRiverLakeFlow
      CALL hist_grid_riverlake_init (HistForm)
#endif

   END SUBROUTINE hist_init


   SUBROUTINE hist_final ()

   IMPLICIT NONE

      CALL deallocate_acc_fluxes ()

#ifdef SinglePoint
      CALL hist_single_final ()
#endif

#ifdef CatchLateralFlow
      CALL hist_basin_final ()
#endif

#ifdef GridRiverLakeFlow
      CALL hist_grid_riverlake_final ()
#endif

   END SUBROUTINE hist_final


   SUBROUTINE hist_out (idate, deltim, itstamp, etstamp, ptstamp, &
         dir_hist, casename)

!=======================================================================
!  Original version: Yongjiu Dai, September 15, 1999, 03/2014
!=======================================================================

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_TimeManager
   USE MOD_SPMD_Task
   USE MOD_Vars_1DAccFluxes, filter_dt_ => filter_dt
   USE MOD_Vars_1DFluxes, only: nsensor
   USE MOD_Vars_TimeVariables, only: wa, wat, wetwat, wdsrf
   USE MOD_Block
   USE MOD_DataType
   USE MOD_LandPatch
   USE MOD_SpatialMapping
   USE MOD_Vars_TimeInvariants, only: patchtype, patchclass, patchmask
#ifdef URBAN_MODEL
   USE MOD_LandUrban
#endif
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_Vars_PFTimeInvariants, only: pftclass
   USE MOD_LandPFT, only: patch_pft_s
#endif
#if (defined CaMa_Flood)
   USE MOD_CaMa_Vars !definition of CaMa variables
#endif
   USE MOD_Forcing, only: forcmask_pch
#ifdef DataAssimilation
   USE MOD_DA_TWS, only: fslp_k_mon
   USE MOD_Vars_Global
   USE MOD_DA_Vars_TimeVariables
   USE MOD_Const_Physical, only: denh2o
#endif

   IMPLICIT NONE

   integer,  intent(in) :: idate(3)
   real(r8), intent(in) :: deltim
   type(timestamp), intent(in) :: itstamp
   type(timestamp), intent(in) :: etstamp
   type(timestamp), intent(in) :: ptstamp

   character(len=*), intent(in) :: dir_hist
   character(len=*), intent(in) :: casename

   ! Local variables
   logical :: lwrite
   character(len=256) :: file_hist
   integer :: itime_in_file
#if (defined CaMa_Flood)
   character(len=256) :: file_hist_cama
   integer :: itime_in_file_cama
#endif
   integer :: month, day
   integer :: days_month(1:12)
   character(len=10) :: cdate

   type(block_data_real8_2d) :: sumarea
   type(block_data_real8_2d) :: sumarea_dt
   type(block_data_real8_2d) :: sumarea_urb
   type(block_data_real8_2d) :: sumarea_one
   real(r8), allocatable ::  vecacc     (:)
   real(r8), allocatable ::  nac_one    (:)
   logical,  allocatable ::  filter     (:)
   logical,  allocatable ::  filter_dt  (:)

#ifdef CROP
   type(block_data_real8_2d) :: sumarea_crop
   type(block_data_real8_2d) :: sumarea_irrig
   logical,  allocatable ::  filter_crop (:)
   logical,  allocatable ::  filter_irrig (:)
#endif

   integer i, u
#ifdef URBAN_MODEL
   logical,  allocatable ::  filter_urb (:)
#endif

#ifdef DataAssimilation
   integer :: np
   real(r8), allocatable ::  a_wliq_h2osoi_5cm (:)
   real(r8), allocatable ::  a_t_soisno_5cm (:)
   real(r8), allocatable ::  a_wliq_soisno_ens_mean (:,:)
   real(r8), allocatable ::  a_wliq_soisno_5cm_ens (:,:)
   real(r8), allocatable ::  a_wliq_h2osoi_5cm_a (:)
   real(r8), allocatable ::  a_t_soisno_ens_mean (:,:)
   real(r8), allocatable ::  a_t_soisno_5cm_ens (:,:)
   real(r8), allocatable ::  a_t_soisno_5cm_a (:)
   real(r8), allocatable ::  a_t_brt_smap_a (:,:)
   real(r8), allocatable ::  a_t_brt_fy3d_a (:,:)
   real(r8), allocatable ::  a_wliq_soisno_5cm_ens_std (:)
   real(r8), allocatable ::  a_t_soisno_5cm_ens_std (:)
   real(r8), allocatable ::  a_t_brt_smap_ens_std (:,:)
   real(r8), allocatable ::  a_t_brt_fy3d_ens_std (:,:)
#endif

      IF (itstamp <= ptstamp) THEN
         CALL FLUSH_acc_fluxes ()
         RETURN
      ELSE
         CALL accumulate_fluxes ()
      ENDIF

      select CASE (trim(adjustl(DEF_HIST_FREQ)))
      CASE ('TIMESTEP')
         lwrite = .true.
      CASE ('HOURLY')
         lwrite = isendofhour (idate, deltim) .or. (.not. (itstamp < etstamp))
      CASE ('DAILY')
         lwrite = isendofday  (idate, deltim) .or. (.not. (itstamp < etstamp))
      CASE ('MONTHLY')
         lwrite = isendofmonth(idate, deltim) .or. (.not. (itstamp < etstamp))
      CASE ('YEARLY')
         lwrite = isendofyear (idate, deltim) .or. (.not. (itstamp < etstamp))
      CASE default
         lwrite = .false.
         write(*,*) &
         'Warning : Please USE one of TIMESTEP/HOURLY/DAILY/MONTHLY/YEARLY for history frequency.'
         write(*,*) &
         '          Set to FALSE by default.                                                     '
      END select

      IF (lwrite) THEN

         CALL julian2monthday(idate(1), idate(2), month, day)

         days_month = (/31,28,31,30,31,30,31,31,30,31,30,31/)
         IF (isleapyear(idate(1))) days_month(2) = 29

         IF ( trim(DEF_HIST_groupby) == 'YEAR' ) THEN
            write(cdate,'(i4.4)') idate(1)
#ifdef SinglePoint
            IF (USE_SITE_HistWriteBack) THEN
               memory_to_disk = isendofyear(idate,deltim) .or. (.not. (itstamp < etstamp))
            ENDIF
#endif
         ELSEIF ( trim(DEF_HIST_groupby) == 'MONTH' ) THEN
            write(cdate,'(i4.4,"-",i2.2)') idate(1), month
#ifdef SinglePoint
            IF (USE_SITE_HistWriteBack) THEN
               memory_to_disk = isendofmonth(idate,deltim) .or. (.not. (itstamp < etstamp))
            ENDIF
#endif
         ELSEIF ( trim(DEF_HIST_groupby) == 'DAY' ) THEN
            write(cdate,'(i4.4,"-",i2.2,"-",i2.2)') idate(1), month, day
#ifdef SinglePoint
            IF (USE_SITE_HistWriteBack) THEN
               memory_to_disk = isendofday(idate,deltim) .or. (.not. (itstamp < etstamp))
            ENDIF
#endif
         ELSE
            write(*,*) 'Warning : Please USE one of DAY/MONTH/YEAR for history group.'
         ENDIF

#if (defined CaMa_Flood)
         ! add variables to write cama-flood output.
         ! file name of cama-flood output
         file_hist_cama = trim(dir_hist) // '/' // trim(casename) //'_hist_cama_'//trim(cdate)//'.nc'
         ! write CaMa-Flood output
         CALL hist_write_cama_time (file_hist_cama, 'time', idate, itime_in_file_cama)
#endif

         file_hist = trim(dir_hist) // '/' // trim(casename) //'_hist_'//trim(cdate)//'.nc'

         CALL hist_write_time (file_hist, file_last, 'time', idate, itime_in_file)

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               allocate (filter    (numpatch))
               allocate (filter_dt (numpatch))
               allocate (vecacc    (numpatch))
            ENDIF
#ifdef URBAN_MODEL
            IF (numurban > 0) THEN
               allocate (filter_urb (numurban))
            ENDIF
#endif
#ifdef CROP
            IF (numpatch > 0) THEN
               allocate (filter_crop (numpatch))
               allocate (filter_irrig (numpatch))
            ENDIF
#endif
         ENDIF

         IF (HistForm == 'Gridded') THEN
            IF (p_is_io) THEN
               CALL allocate_block_data (ghist, sumarea)
               CALL allocate_block_data (ghist, sumarea_dt)
#ifdef URBAN_MODEL
               CALL allocate_block_data (ghist, sumarea_urb)
#endif
#ifdef CROP
               CALL allocate_block_data (ghist, sumarea_crop)
               CALL allocate_block_data (ghist, sumarea_irrig)
#endif
            ENDIF
         ENDIF

         ! ---------------------------------------------------
         ! Meteorological forcing and patch mask filter applying.
         ! ---------------------------------------------------
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN

               filter(:) = patchtype < 99

               IF (DEF_forcing%has_missing_value) THEN
                  filter = filter .and. forcmask_pch
               ENDIF

               filter    = filter .and. patchmask
               filter_dt = filter .and. patchmask .and. (nac_dt > 0)

            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
            CALL mp2g_hist%get_sumarea (sumarea_dt, filter_dt)
         ENDIF

         IF (HistForm == 'Gridded') THEN
            IF (trim(file_hist) /= trim(file_last)) THEN
               CALL hist_write_var_real8_2d (file_hist, 'landarea', ghist, -1, sumarea, &
                  compress = 1, longname = 'land area', units = 'km2')
               CALL hist_write_var_real8_2d (file_hist, 'landfraction', ghist, -1, landfraction, &
                  compress = 1, longname = 'land fraction', units = '-')
            ENDIF
         ENDIF

#ifdef CROP
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               filter_crop(:) = patchtype < 99
               filter_crop(:) = patchclass == 12

               IF (DEF_forcing%has_missing_value) THEN
                  filter_crop = filter_crop .and. forcmask_pch
               ENDIF

               filter_crop = filter_crop .and. patchmask
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea_crop, filter_crop)
         ENDIF

         IF (HistForm == 'Gridded') THEN
            IF (trim(file_hist) /= trim(file_last)) THEN
               CALL hist_write_var_real8_2d (file_hist, 'croparea', ghist, 1, sumarea_crop, &
                  compress = 1, longname = 'crop area', units = 'km2')
            ENDIF
         ENDIF

         IF (DEF_USE_IRRIGATION) THEN
            IF (p_is_worker) THEN
               IF (numpatch > 0) THEN
                  DO i=1,numpatch
                     IF(patchclass(i) == 12)THEN
                        IF((pftclass(patch_pft_s(i)).GE.npcropmin).and.(MOD(pftclass(patch_pft_s(i)),2).EQ.0))THEN
                           filter_irrig(i) = .true.
                        ELSE
                           filter_irrig(i) = .false.
                        ENDIF
                     ELSE
                        filter_irrig(i) = .false.
                     ENDIF
                  ENDDO
                  IF (DEF_forcing%has_missing_value) THEN
                     filter_irrig = filter_irrig .and. forcmask_pch
                  ENDIF
                  filter_irrig = filter_irrig .and. patchmask
               ENDIF
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea_irrig, filter_irrig)
         ENDIF

         IF (HistForm == 'Gridded') THEN
            IF (trim(file_hist) /= trim(file_last)) THEN
               CALL hist_write_var_real8_2d (file_hist, 'irrigarea', ghist, 1, sumarea_irrig, &
                  compress = 1, longname = 'irrigation area', units = 'km2')
            ENDIF
         ENDIF
#endif

         ! wind in eastward direction [m/s]
         CALL write_history_variable_2d ( DEF_hist_vars%xy_us, &
            a_us, file_hist, 'f_xy_us', itime_in_file, sumarea, filter, &
            'wind in eastward direction', 'm/s')

         ! wind in northward direction [m/s]
         CALL write_history_variable_2d ( DEF_hist_vars%xy_vs, &
            a_vs, file_hist, 'f_xy_vs', itime_in_file, sumarea, filter, &
            'wind in northward direction','m/s')

         ! temperature at reference height [kelvin]
         CALL write_history_variable_2d ( DEF_hist_vars%xy_t, &
            a_t, file_hist, 'f_xy_t', itime_in_file, sumarea, filter, &
            'temperature at reference height','kelvin')

         ! specific humidity at reference height [kg/kg]
         CALL write_history_variable_2d ( DEF_hist_vars%xy_q, &
            a_q, file_hist, 'f_xy_q', itime_in_file, sumarea, filter, &
            'specific humidity at reference height','kg/kg')

         ! convective precipitation [mm/s]
         CALL write_history_variable_2d ( DEF_hist_vars%xy_prc, &
            a_prc, file_hist, 'f_xy_prc', itime_in_file, sumarea, filter, &
            'convective precipitation','mm/s')

         ! large scale precipitation [mm/s]
         CALL write_history_variable_2d ( DEF_hist_vars%xy_prl, &
            a_prl, file_hist, 'f_xy_prl', itime_in_file, sumarea, filter, &
            'large scale precipitation','mm/s')

         ! atmospheric pressure at the surface [pa]
         CALL write_history_variable_2d ( DEF_hist_vars%xy_pbot, &
            a_pbot, file_hist, 'f_xy_pbot', itime_in_file, sumarea, filter, &
            'atmospheric pressure at the surface','pa')

         ! atmospheric infrared (longwave) radiation [W/m2]
         CALL write_history_variable_2d ( DEF_hist_vars%xy_frl, &
            a_frl, file_hist, 'f_xy_frl', itime_in_file, sumarea, filter, &
            'atmospheric infrared (longwave) radiation','W/m2')

         ! downward solar radiation at surface [W/m2]
         CALL write_history_variable_2d ( DEF_hist_vars%xy_solarin, &
            a_solarin, file_hist, 'f_xy_solarin', itime_in_file, sumarea, filter, &
            'downward solar radiation at surface','W/m2')

         ! rain [mm/s]
         CALL write_history_variable_2d ( DEF_hist_vars%xy_rain, &
            a_rain, file_hist, 'f_xy_rain', itime_in_file, sumarea, filter, &
            'rain','mm/s')

         ! snow [mm/s]
         CALL write_history_variable_2d ( DEF_hist_vars%xy_snow, &
            a_snow, file_hist, 'f_xy_snow', itime_in_file, sumarea, filter, &
            'snow','mm/s')

         IF (DEF_USE_CBL_HEIGHT) THEN
         ! atmospheric boundary layer height [m]
           CALL write_history_variable_2d ( DEF_hist_vars%xy_hpbl, &
              a_hpbl, file_hist, 'f_xy_hpbl', itime_in_file, sumarea, filter, &
              'boundary layer height','m')
         ENDIF

         IF (DEF_USE_OZONESTRESS)THEN
            CALL write_history_variable_2d ( DEF_hist_vars%o3uptakesun, &
              a_o3uptakesun, file_hist, 'f_o3uptakesun', itime_in_file, sumarea, filter, &
              'Accumulated ozone uptake by sunlit leaf','mmol O3 m-2')

            CALL write_history_variable_2d ( DEF_hist_vars%o3uptakesha, &
              a_o3uptakesha, file_hist, 'f_o3uptakesha', itime_in_file, sumarea, filter, &
              'Accumulated ozone uptake by shaded leaf','mmol O3 m-2')
         ENDIF

         ! ------------------------------------------------------------------------------------------
         ! Mapping the fluxes and state variables at patch [numpatch] to grid
         ! ------------------------------------------------------------------
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN

               filter(:) = patchtype < 99

               IF (DEF_forcing%has_missing_value) THEN
                  filter = filter .and. forcmask_pch
               ENDIF

               filter = filter .and. patchmask
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! wind stress: E-W [kg/m/s2]
         CALL write_history_variable_2d ( DEF_hist_vars%taux, &
            a_taux, file_hist, 'f_taux', itime_in_file, sumarea, filter, &
            'wind stress: E-W','kg/m/s2')

         ! wind stress: N-S [kg/m/s2]
         CALL write_history_variable_2d ( DEF_hist_vars%tauy, &
            a_tauy, file_hist, 'f_tauy', itime_in_file, sumarea, filter, &
            'wind stress: N-S','kg/m/s2')

         ! sensible heat from canopy height to atmosphere [W/m2]
         CALL write_history_variable_2d ( DEF_hist_vars%fsena, &
            a_fsena, file_hist, 'f_fsena', itime_in_file, sumarea, filter, &
            'sensible heat from canopy height to atmosphere','W/m2')

         ! latent heat flux from canopy height to atmosphere [W/m2]
         CALL write_history_variable_2d ( DEF_hist_vars%lfevpa, &
            a_lfevpa, file_hist, 'f_lfevpa', itime_in_file, sumarea, filter, &
            'latent heat flux from canopy height to atmosphere','W/m2')

         ! evapotranspiration from canopy to atmosphere [mm/s]
         CALL write_history_variable_2d ( DEF_hist_vars%fevpa, &
            a_fevpa, file_hist, 'f_fevpa', itime_in_file, sumarea, filter, &
            'evapotranspiration from canopy height to atmosphere','mm/s')

         ! sensible heat from leaves [W/m2]
         CALL write_history_variable_2d ( DEF_hist_vars%fsenl, &
            a_fsenl, file_hist, 'f_fsenl', itime_in_file, sumarea, filter, &
            'sensible heat from leaves','W/m2')

         ! evaporation+transpiration from leaves [mm/s]
         CALL write_history_variable_2d ( DEF_hist_vars%fevpl, &
            a_fevpl, file_hist, 'f_fevpl', itime_in_file, sumarea, filter, &
            'evaporation+transpiration from leaves','mm/s')

         ! transpiration rate [mm/s]
         CALL write_history_variable_2d ( DEF_hist_vars%etr, &
            a_etr, file_hist, 'f_etr', itime_in_file, sumarea, filter, &
            'transpiration rate','mm/s')

         ! sensible heat flux from ground [W/m2]
         CALL write_history_variable_2d ( DEF_hist_vars%fseng, &
            a_fseng, file_hist, 'f_fseng', itime_in_file, sumarea, filter, &
            'sensible heat flux from ground','W/m2')

         ! evaporation heat flux from ground [mm/s]
         CALL write_history_variable_2d ( DEF_hist_vars%fevpg, &
            a_fevpg, file_hist, 'f_fevpg', itime_in_file, sumarea, filter, &
            'evaporation flux from ground','mm/s')

         ! ground heat flux [W/m2]
         CALL write_history_variable_2d ( DEF_hist_vars%fgrnd, &
            a_fgrnd, file_hist, 'f_fgrnd', itime_in_file, sumarea, filter, &
            'ground heat flux','W/m2')

         ! solar absorbed by sunlit canopy [W/m2]
         CALL write_history_variable_2d ( DEF_hist_vars%sabvsun, &
            a_sabvsun, file_hist, 'f_sabvsun', itime_in_file, sumarea, filter, &
            'solar absorbed by sunlit canopy','W/m2')

         ! solar absorbed by shaded [W/m2]
         CALL write_history_variable_2d ( DEF_hist_vars%sabvsha, &
            a_sabvsha, file_hist, 'f_sabvsha', itime_in_file, sumarea, filter, &
            'solar absorbed by shaded','W/m2')

         ! solar absorbed by ground  [W/m2]
         CALL write_history_variable_2d ( DEF_hist_vars%sabg, &
            a_sabg, file_hist, 'f_sabg', itime_in_file, sumarea, filter, &
            'solar absorbed by ground','W/m2')

         ! outgoing long-wave radiation from ground+canopy [W/m2]
         CALL write_history_variable_2d ( DEF_hist_vars%olrg, &
            a_olrg, file_hist, 'f_olrg', itime_in_file, sumarea, filter, &
            'outgoing long-wave radiation from ground+canopy','W/m2')

         ! net radiation [W/m2]
         CALL write_history_variable_2d ( DEF_hist_vars%rnet, &
            a_rnet, file_hist, 'f_rnet', itime_in_file, sumarea, filter, &
            'net radiation','W/m2')

         ! the error of water balance [mm/s]
         CALL write_history_variable_2d ( DEF_hist_vars%xerr, &
            a_xerr, file_hist, 'f_xerr', itime_in_file, sumarea, filter, &
            'the error of water banace','mm/s')

         ! the error of energy balance [W/m2]
         CALL write_history_variable_2d ( DEF_hist_vars%zerr, &
            a_zerr, file_hist, 'f_zerr', itime_in_file, sumarea, filter, &
            'the error of energy balance','W/m2')

         ! surface runoff [mm/s]
         CALL write_history_variable_2d ( DEF_hist_vars%rsur, &
            a_rsur, file_hist, 'f_rsur', itime_in_file, sumarea, filter, &
            'surface runoff','mm/s')

#ifndef CatchLateralFlow
         ! saturation excess surface runoff [mm/s]
         CALL write_history_variable_2d ( DEF_hist_vars%rsur_se, &
            a_rsur_se, file_hist, 'f_rsur_se', itime_in_file, sumarea, filter, &
            'saturation excess surface runoff','mm/s')

         ! infiltration excess surface runoff [mm/s]
         CALL write_history_variable_2d ( DEF_hist_vars%rsur_ie, &
            a_rsur_ie, file_hist, 'f_rsur_ie', itime_in_file, sumarea, filter, &
            'infiltration excess surface runoff','mm/s')
#endif

         ! subsurface runoff [mm/s]
         CALL write_history_variable_2d ( DEF_hist_vars%rsub, &
            a_rsub, file_hist, 'f_rsub', itime_in_file, sumarea, filter, &
            'subsurface runoff','mm/s')

         ! total runoff [mm/s]
         CALL write_history_variable_2d ( DEF_hist_vars%rnof, &
            a_rnof, file_hist, 'f_rnof', itime_in_file, sumarea, filter, &
            'total runoff','mm/s')

#ifdef DataAssimilation
         IF (DEF_DA_TWS_GRACE) THEN
            ! slope factors for runoff [-]
            IF (p_is_worker .and. (numpatch > 0)) THEN
               vecacc = fslp_k_mon(month, :)
               WHERE (vecacc /= spval) vecacc = vecacc*nac
            ENDIF
            CALL write_history_variable_2d(.true., &
               vecacc, file_hist, 'f_slope_factor_k', itime_in_file, sumarea, filter, &
               'slope factor [k] for runoff', '-')
         ENDIF
#endif

#ifdef CatchLateralFlow
         ! rate of surface water depth change [mm/s]
         CALL write_history_variable_2d ( DEF_hist_vars%xwsur, &
            a_xwsur, file_hist, 'f_xwsur', itime_in_file, sumarea, filter, &
            'rate of surface water depth change','mm/s')

         ! rate of ground water change [mm/s]
         CALL write_history_variable_2d ( DEF_hist_vars%xwsub, &
            a_xwsub, file_hist, 'f_xwsub', itime_in_file, sumarea, filter, &
            'rate of ground water change','mm/s')

         ! fraction of flooded area [-]
         CALL write_history_variable_2d ( DEF_hist_vars%fldarea, &
            a_fldarea, file_hist, 'f_fldarea', itime_in_file, sumarea, filter, &
            'fraction of flooded area','-')
#endif

         ! interception [mm/s]
         CALL write_history_variable_2d ( DEF_hist_vars%qintr, &
            a_qintr, file_hist, 'f_qintr', itime_in_file, sumarea, filter, &
            'interception','mm/s')

         ! infiltration [mm/s]
         CALL write_history_variable_2d ( DEF_hist_vars%qinfl, &
            a_qinfl, file_hist, 'f_qinfl', itime_in_file, sumarea, filter, &
            'f_qinfl','mm/s')

         ! throughfall [mm/s]
         CALL write_history_variable_2d ( DEF_hist_vars%qdrip, &
            a_qdrip, file_hist, 'f_qdrip', itime_in_file, sumarea, filter, &
            'total throughfall','mm/s')

         ! total water storage [mm]
         CALL write_history_variable_2d ( DEF_hist_vars%wat, &
            a_wat, file_hist, 'f_wat', itime_in_file, sumarea, filter, &
            'total water storage','mm')

         ! instantaneous total water storage [mm]
         IF (p_is_worker .and. (numpatch > 0)) THEN
            vecacc = wat
            WHERE(vecacc /= spval) vecacc = vecacc * nac
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%wat_inst, &
            vecacc, file_hist, 'f_wat_inst', itime_in_file, sumarea, filter, &
            'instantaneous total water storage','mm')

         ! canopy assimilation rate [mol m-2 s-1]
         CALL write_history_variable_2d ( DEF_hist_vars%assim, &
            a_assim, file_hist, 'f_assim', itime_in_file, sumarea, filter, &
            'canopy assimilation rate','mol m-2 s-1')

         ! respiration (plant+soil) [mol m-2 s-1]
         CALL write_history_variable_2d ( DEF_hist_vars%respc, &
            a_respc, file_hist, 'f_respc', itime_in_file, sumarea, filter, &
            'respiration (plant+soil)','mol m-2 s-1')

         ! groundwater recharge rate [mm/s]
         CALL write_history_variable_2d ( DEF_hist_vars%qcharge &
            .and. (.not.DEF_USE_VariablySaturatedFlow), &
            a_qcharge, file_hist, 'f_qcharge', itime_in_file, sumarea, filter, &
            'groundwater recharge rate','mm/s')

         ! ground surface temperature [K]
         CALL write_history_variable_2d ( DEF_hist_vars%t_grnd, &
            a_t_grnd, file_hist, 'f_t_grnd', itime_in_file, sumarea, filter, &
            'ground surface temperature','K')

         ! leaf temperature [K]
         CALL write_history_variable_2d ( DEF_hist_vars%tleaf, &
            a_tleaf, file_hist, 'f_tleaf', itime_in_file, sumarea, filter, &
            'leaf temperature','K')

         ! depth of water on foliage [mm]
         CALL write_history_variable_2d ( DEF_hist_vars%ldew, &
            a_ldew, file_hist, 'f_ldew', itime_in_file, sumarea, filter, &
            'depth of water on foliage','mm')

         ! snow cover, water equivalent [mm]
         CALL write_history_variable_2d ( DEF_hist_vars%scv, &
            a_scv, file_hist, 'f_scv', itime_in_file, sumarea, filter, &
            'snow cover, water equivalent','mm')

         ! snow depth [meter]
         CALL write_history_variable_2d ( DEF_hist_vars%snowdp, &
            a_snowdp, file_hist, 'f_snowdp', itime_in_file, sumarea, filter, &
            'snow depth','meter')

         ! fraction of snow cover on ground
         CALL write_history_variable_2d ( DEF_hist_vars%fsno, &
            a_fsno, file_hist, 'f_fsno', itime_in_file, sumarea, filter, &
            'fraction of snow cover on ground','-')

         ! fraction of saturation area on ground
         CALL write_history_variable_2d ( DEF_hist_vars%frcsat, &
            a_frcsat, file_hist, 'f_frcsat', itime_in_file, sumarea, filter, &
            'fraction of saturation area on ground','-')

         ! fraction of veg cover, excluding snow-covered veg [-]
         CALL write_history_variable_2d ( DEF_hist_vars%sigf, &
            a_sigf, file_hist, 'f_sigf', itime_in_file, sumarea, filter, &
            'fraction of veg cover, excluding snow-covered veg','-')

         ! leaf greenness
         CALL write_history_variable_2d ( DEF_hist_vars%green, &
            a_green, file_hist, 'f_green', itime_in_file, sumarea, filter, &
            'leaf greenness','-')

         ! leaf area index
         CALL write_history_variable_2d ( DEF_hist_vars%lai, &
            a_lai, file_hist, 'f_lai', itime_in_file, sumarea, filter, &
            'leaf area index','m2/m2')

         ! leaf area index
         CALL write_history_variable_2d ( DEF_hist_vars%laisun, &
            a_laisun, file_hist, 'f_laisun', itime_in_file, sumarea, filter, &
            'sunlit leaf area index','m2/m2')

         ! leaf area index
         CALL write_history_variable_2d ( DEF_hist_vars%laisha, &
            a_laisha, file_hist, 'f_laisha', itime_in_file, sumarea, filter, &
            'shaded leaf area index','m2/m2')

         ! stem area index
         CALL write_history_variable_2d ( DEF_hist_vars%sai, &
            a_sai, file_hist, 'f_sai', itime_in_file, sumarea, filter, &
            'stem area index','m2/m2')

         ! averaged albedo [visible, direct; direct, diffuse]
         CALL write_history_variable_4d ( DEF_hist_vars%alb, &
            a_alb, file_hist, 'f_alb', itime_in_file, &
            'band', 1, 2, 'rtyp', 1, 2, sumarea_dt, filter_dt, &
            'averaged albedo','-',nac_dt)

         ! averaged bulk surface emissivity
         CALL write_history_variable_2d ( DEF_hist_vars%emis, &
            a_emis, file_hist, 'f_emis', itime_in_file, sumarea, filter, &
            'averaged bulk surface emissivity','-')

         ! effective roughness [m]
         CALL write_history_variable_2d ( DEF_hist_vars%z0m, &
            a_z0m, file_hist, 'f_z0m', itime_in_file, sumarea, filter, &
            'effective roughness','m')

         ! radiative temperature of surface [K]
         CALL write_history_variable_2d ( DEF_hist_vars%trad, &
            a_trad, file_hist, 'f_trad', itime_in_file, sumarea, filter, &
            'radiative temperature of surface','kelvin')

         ! 2 m height air temperature [kelvin]
         CALL write_history_variable_2d ( DEF_hist_vars%tref, &
            a_tref, file_hist, 'f_tref', itime_in_file, sumarea, filter, &
            '2 m height air temperature','kelvin')

IF (DEF_Output_2mWMO) THEN
         ! 2 m WMO air temperature [kelvin]
         CALL write_history_variable_2d ( DEF_hist_vars%t2m_wmo, &
            a_t2m_wmo, file_hist, 'f_t2m_wmo', itime_in_file, sumarea, filter, &
            '2 m WMO air temperature','kelvin')
ENDIF

         ! 2 m height air specific humidity [kg/kg]
         CALL write_history_variable_2d ( DEF_hist_vars%qref, &
            a_qref, file_hist, 'f_qref', itime_in_file, sumarea, filter, &
            '2 m height air specific humidity','kg/kg')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN

               filter(:) = patchtype == 2

               IF (DEF_forcing%has_missing_value) THEN
                  filter = filter .and. forcmask_pch
               ENDIF

               filter = filter .and. patchmask
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         IF (HistForm == 'Gridded') THEN
            IF (trim(file_hist) /= trim(file_last)) THEN
               CALL hist_write_var_real8_2d (file_hist, 'area_wetland', ghist, -1, sumarea, &
                  compress = 1, longname = 'area of wetland', units = 'km2')
            ENDIF
         ENDIF

         ! wetland water storage [mm]
         IF (DEF_USE_Dynamic_Wetland) THEN
            IF (p_is_worker .and. (numpatch > 0)) THEN
               vecacc = a_wdsrf
            ENDIF
            CALL write_history_variable_2d ( DEF_hist_vars%wetwat, &
               vecacc, file_hist, 'f_wetwat', itime_in_file, sumarea, filter, &
               'wetland water storage','mm')
         ELSE
            CALL write_history_variable_2d ( DEF_hist_vars%wetwat, &
               a_wetwat, file_hist, 'f_wetwat', itime_in_file, sumarea, filter, &
               'wetland water storage','mm')
         ENDIF

         ! instantaneous wetland water storage [mm]
         IF (p_is_worker .and. (numpatch > 0)) THEN
            vecacc = wetwat
            WHERE(vecacc /= spval) vecacc = vecacc * nac
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%wetwat_inst, &
            vecacc, file_hist, 'f_wetwat_inst', itime_in_file, sumarea, filter, &
            'instantaneous wetland water storage','mm')

         IF (p_is_worker .and. (numpatch > 0)) THEN
            vecacc = a_zwt
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%wetzwt, &
            vecacc, file_hist, 'f_wetzwt', itime_in_file, sumarea, filter, &
            'the depth to water table in wetland','m')

         ! ------------------------------------------------------------------
         ! Mapping the urban variables at patch [numurban] to grid
         ! ------------------------------------------------------------------

#ifdef URBAN_MODEL
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i = 1, numpatch
                  IF (patchtype(i) == 1) THEN
                     u = patch2urban(i)

                     filter_urb(u) = .true.

                     IF (DEF_forcing%has_missing_value) THEN
                        filter_urb(u) = filter_urb(u) .and. forcmask_pch(i)
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist_urb%get_sumarea (sumarea_urb, filter_urb)
         ENDIF

         ! sensible heat from building roof [W/m2]
         CALL write_history_variable_urb_2d ( DEF_hist_vars%fsen_roof, &
            a_senroof, file_hist, 'f_fsenroof', itime_in_file, sumarea_urb, filter_urb, &
            'sensible heat from urban roof [W/m2]','W/m2')

         ! sensible heat from building sunlit wall [W/m2]
         CALL write_history_variable_urb_2d ( DEF_hist_vars%fsen_wsun, &
            a_senwsun, file_hist, 'f_fsenwsun', itime_in_file, sumarea_urb, filter_urb, &
            'sensible heat from urban sunlit wall [W/m2]','W/m2')

         ! sensible heat from building shaded wall [W/m2]
         CALL write_history_variable_urb_2d ( DEF_hist_vars%fsen_wsha, &
            a_senwsha, file_hist, 'f_fsenwsha', itime_in_file, sumarea_urb, filter_urb, &
            'sensible heat from urban shaded wall [W/m2]','W/m2')

         ! sensible heat from impervious ground [W/m2]
         CALL write_history_variable_urb_2d ( DEF_hist_vars%fsen_gimp, &
            a_sengimp, file_hist, 'f_fsengimp', itime_in_file, sumarea_urb, filter_urb, &
            'sensible heat from urban impervious ground [W/m2]','W/m2')

         ! sensible heat from pervious ground [W/m2]
         CALL write_history_variable_urb_2d ( DEF_hist_vars%fsen_gper, &
            a_sengper, file_hist, 'f_fsengper', itime_in_file, sumarea_urb, filter_urb, &
            'sensible heat from urban pervious ground [W/m2]','W/m2')

         ! sensible heat from urban tree [W/m2]
         CALL write_history_variable_urb_2d ( DEF_hist_vars%fsen_urbl, &
            a_senurbl, file_hist, 'f_fsenurbl', itime_in_file, sumarea_urb, filter_urb, &
            'sensible heat from urban tree [W/m2]','W/m2')

         ! latent heat flux from building roof [W/m2]
         CALL write_history_variable_urb_2d ( DEF_hist_vars%lfevp_roof, &
            a_lfevproof, file_hist, 'f_lfevproof', itime_in_file, sumarea_urb, filter_urb, &
            'latent heat from urban roof [W/m2]','W/m2')

         ! latent heat flux from impervious ground [W/m2]
         CALL write_history_variable_urb_2d ( DEF_hist_vars%lfevp_gimp, &
            a_lfevpgimp, file_hist, 'f_lfevpgimp', itime_in_file, sumarea_urb, filter_urb, &
            'latent heat from urban impervious ground [W/m2]','W/m2')

         ! latent heat flux from pervious ground [W/m2]
         CALL write_history_variable_urb_2d ( DEF_hist_vars%lfevp_gper, &
            a_lfevpgper, file_hist, 'f_lfevpgper', itime_in_file, sumarea_urb, filter_urb, &
            'latent heat from urban pervious ground [W/m2]','W/m2')

         ! latent heat flux from urban tree [W/m2]
         CALL write_history_variable_urb_2d ( DEF_hist_vars%lfevp_urbl, &
            a_lfevpurbl, file_hist, 'f_lfevpurbl', itime_in_file, sumarea_urb, filter_urb, &
            'latent heat from urban tree [W/m2]','W/m2')

         ! sensible flux from heat or cool AC [W/m2]
         CALL write_history_variable_urb_2d ( DEF_hist_vars%fhac, &
            a_fhac, file_hist, 'f_fhac', itime_in_file, sumarea_urb, filter_urb, &
            'sensible flux from heat or cool AC [W/m2]','W/m2')

         ! waste heat flux from heat or cool AC [W/m2]
         CALL write_history_variable_urb_2d ( DEF_hist_vars%fwst, &
            a_fwst, file_hist, 'f_fwst', itime_in_file, sumarea_urb, filter_urb, &
            'waste heat flux from heat or cool AC [W/m2]','W/m2')

         ! flux from inner and outer air exchange [W/m2]
         CALL write_history_variable_urb_2d ( DEF_hist_vars%fach, &
            a_fach, file_hist, 'f_fach', itime_in_file, sumarea_urb, filter_urb, &
            'flux from inner and outter air exchange [W/m2]','W/m2')

         ! flux from total heating/cooling [W/m2]
         CALL write_history_variable_urb_2d ( DEF_hist_vars%fhah, &
            a_fhah, file_hist, 'f_fhah', itime_in_file, sumarea_urb, filter_urb, &
            'flux from heating/cooling [W/m2]','W/m2')

         ! flux from metabolism [W/m2]
         CALL write_history_variable_urb_2d ( DEF_hist_vars%meta, &
            a_meta, file_hist, 'f_fmeta', itime_in_file, sumarea_urb, filter_urb, &
            'flux from human metabolism [W/m2]','W/m2')

         ! flux from vehicle [W/m2]
         CALL write_history_variable_urb_2d ( DEF_hist_vars%vehc, &
            a_vehc, file_hist, 'f_fvehc', itime_in_file, sumarea_urb, filter_urb, &
            'flux from traffic [W/m2]','W/m2')

         ! temperature of inner building [K]
         CALL write_history_variable_urb_2d ( DEF_hist_vars%t_room, &
            a_t_room, file_hist, 'f_t_room', itime_in_file, sumarea_urb, filter_urb, &
            'temperature of inner building [K]','kelvin')

         ! temperature of outer building [K]
         CALL write_history_variable_urb_2d ( DEF_hist_vars%tafu, &
            a_tafu, file_hist, 'f_tafu', itime_in_file, sumarea_urb, filter_urb, &
            'temperature of outer building [K]','kelvin')

         ! temperature of building roof [K]
         CALL write_history_variable_urb_2d ( DEF_hist_vars%t_roof, &
            a_troof, file_hist, 'f_t_roof', itime_in_file, sumarea_urb, filter_urb, &
            'temperature of urban roof [K]','kelvin')

         ! temperature of building wall [K]
         CALL write_history_variable_urb_2d ( DEF_hist_vars%t_wall, &
            a_twall, file_hist, 'f_t_wall', itime_in_file, sumarea_urb, filter_urb, &
            'temperature of urban wall [K]','kelvin')
#endif

         ! ------------------------------------------------------------------
         ! Mapping the fluxes and state variables at patch [numpatch] to grid
         ! ------------------------------------------------------------------
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               filter(:) = patchtype < 99
               IF (DEF_forcing%has_missing_value) THEN
                  filter = filter .and. forcmask_pch
               ENDIF
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! 1: assimsun enf temperate
         CALL write_history_variable_2d ( DEF_hist_vars%assimsun, &
             a_assimsun, file_hist, 'f_assimsun', itime_in_file, sumarea, filter, &
         'Photosynthetic assimilation rate of sunlit leaf for needleleaf evergreen temperate tree',&
         'mol m-2 s-1')

         ! 1: assimsha enf temperate
         CALL write_history_variable_2d ( DEF_hist_vars%assimsha, &
             a_assimsha, file_hist, 'f_assimsha', itime_in_file, sumarea, filter, &
         'Photosynthetic assimilation rate of shaded leaf for needleleaf evergreen temperate tree',&
         'mol m-2 s-1')

         ! 1: etrsun enf temperate
         CALL write_history_variable_2d ( DEF_hist_vars%etrsun, &
             a_etrsun, file_hist, 'f_etrsun', itime_in_file, sumarea, filter, &
             'Transpiration rate of sunlit leaf for needleleaf evergreen temperate tree','mm s-1')

         ! 1: etrsha enf temperate
         CALL write_history_variable_2d ( DEF_hist_vars%etrsha, &
             a_etrsha, file_hist, 'f_etrsha', itime_in_file, sumarea, filter, &
             'Transpiration rate of shaded leaf for needleleaf evergreen temperate tree','mm s-1')

         ! rstfacsun
         CALL write_history_variable_2d ( DEF_hist_vars%rstfacsun, &
             a_rstfacsun, file_hist, 'f_rstfacsun', itime_in_file, sumarea, filter, &
             'Ecosystem level Water stress factor on sunlit canopy','unitless')

         ! rstfacsha
         CALL write_history_variable_2d ( DEF_hist_vars%rstfacsha, &
             a_rstfacsha, file_hist, 'f_rstfacsha', itime_in_file, sumarea, filter, &
             'Ecosystem level Water stress factor on shaded canopy','unitless')

         ! gssun
         CALL write_history_variable_2d ( DEF_hist_vars%gssun, &
             a_gssun, file_hist, 'f_gssun', itime_in_file, sumarea, filter, &
             'Ecosystem level canopy conductance on sunlit canopy','mol m-2 s-1')

         ! gssha
         CALL write_history_variable_2d ( DEF_hist_vars%gssha, &
             a_gssha, file_hist, 'f_gssha', itime_in_file, sumarea, filter, &
             'Ecosystem level canopy conductance on shaded canopy','mol m-2 s-1')

         ! soil resistance [m/s]
         CALL write_history_variable_2d ( DEF_hist_vars%rss, &
             a_rss, file_hist, 'f_rss', itime_in_file, sumarea, filter, &
             'soil surface resistance','s/m')

#ifdef BGC
         ! leaf carbon display pool
         CALL write_history_variable_2d ( DEF_hist_vars%leafc, &
             a_leafc, file_hist, 'f_leafc', itime_in_file, sumarea, filter, &
             'leaf carbon display pool','gC/m2')

         ! leaf carbon storage pool
         CALL write_history_variable_2d ( DEF_hist_vars%leafc_storage, &
             a_leafc_storage, file_hist, 'f_leafc_storage', itime_in_file, sumarea, filter, &
             'leaf carbon storage pool','gC/m2')

         ! leaf carbon transfer pool
         CALL write_history_variable_2d ( DEF_hist_vars%leafc_xfer, &
             a_leafc_xfer, file_hist, 'f_leafc_xfer', itime_in_file, sumarea, filter, &
             'leaf carbon transfer pool','gC/m2')

         ! fine root carbon display pool
         CALL write_history_variable_2d ( DEF_hist_vars%frootc, &
             a_frootc, file_hist, 'f_frootc', itime_in_file, sumarea, filter, &
             'fine root carbon display pool','gC/m2')

         ! fine root carbon storage pool
         CALL write_history_variable_2d ( DEF_hist_vars%frootc_storage, &
             a_frootc_storage, file_hist, 'f_frootc_storage', itime_in_file, sumarea, filter, &
             'fine root carbon storage pool','gC/m2')

         ! fine root carbon transfer pool
         CALL write_history_variable_2d ( DEF_hist_vars%frootc_xfer, &
             a_frootc_xfer, file_hist, 'f_frootc_xfer', itime_in_file, sumarea, filter, &
             'fine root carbon transfer pool','gC/m2')

         ! live stem carbon display pool
         CALL write_history_variable_2d ( DEF_hist_vars%livestemc, &
             a_livestemc, file_hist, 'f_livestemc', itime_in_file, sumarea, filter, &
             'live stem carbon display pool','gC/m2')

         ! live stem carbon storage pool
         CALL write_history_variable_2d ( DEF_hist_vars%livestemc_storage, &
             a_livestemc_storage, file_hist, 'f_livestemc_storage', itime_in_file, sumarea, filter,&
             'live stem carbon storage pool','gC/m2')

         ! live stem carbon transfer pool
         CALL write_history_variable_2d ( DEF_hist_vars%livestemc_xfer, &
             a_livestemc_xfer, file_hist, 'f_livestemc_xfer', itime_in_file, sumarea, filter, &
             'live stem carbon transfer pool','gC/m2')

         ! dead stem carbon display pool
         CALL write_history_variable_2d ( DEF_hist_vars%deadstemc, &
             a_deadstemc, file_hist, 'f_deadstemc', itime_in_file, sumarea, filter, &
             'dead stem carbon display pool','gC/m2')

         ! dead stem carbon storage pool
         CALL write_history_variable_2d ( DEF_hist_vars%deadstemc_storage, &
             a_deadstemc_storage, file_hist, 'f_deadstemc_storage', itime_in_file, sumarea, filter,&
             'dead stem carbon storage pool','gC/m2')

         ! dead stem carbon transfer pool
         CALL write_history_variable_2d ( DEF_hist_vars%deadstemc_xfer, &
             a_deadstemc_xfer, file_hist, 'f_deadstemc_xfer', itime_in_file, sumarea, filter, &
             'dead stem carbon transfer pool','gC/m2')

         ! live coarse root carbon display pool
         CALL write_history_variable_2d ( DEF_hist_vars%livecrootc, &
             a_livecrootc, file_hist, 'f_livecrootc', itime_in_file, sumarea, filter, &
             'live coarse root carbon display pool','gC/m2')

         ! live coarse root carbon storage pool
         CALL write_history_variable_2d ( DEF_hist_vars%livecrootc_storage, &
             a_livecrootc_storage, file_hist, 'f_livecrootc_storage', &
             itime_in_file, sumarea, filter, &
             'live coarse root carbon storage pool','gC/m2')

         ! live coarse root carbon transfer pool
         CALL write_history_variable_2d ( DEF_hist_vars%livecrootc_xfer, &
             a_livecrootc_xfer, file_hist, 'f_livecrootc_xfer', itime_in_file, sumarea, filter, &
             'live coarse root carbon transfer pool','gC/m2')

         ! dead coarse root carbon display pool
         CALL write_history_variable_2d ( DEF_hist_vars%deadcrootc, &
             a_deadcrootc, file_hist, 'f_deadcrootc', itime_in_file, sumarea, filter, &
             'dead coarse root carbon display pool','gC/m2')

         ! dead coarse root carbon storage pool
         CALL write_history_variable_2d ( DEF_hist_vars%deadcrootc_storage, &
             a_deadcrootc_storage, file_hist, 'f_deadcrootc_storage', &
             itime_in_file, sumarea, filter, &
             'dead coarse root carbon storage pool','gC/m2')

         ! dead coarse root carbon transfer pool
         CALL write_history_variable_2d ( DEF_hist_vars%deadcrootc_xfer, &
             a_deadcrootc_xfer, file_hist, 'f_deadcrootc_xfer', itime_in_file, sumarea, filter, &
             'dead coarse root carbon transfer pool','gC/m2')

#ifdef CROP
         ! grain carbon display pool
         CALL write_history_variable_2d ( DEF_hist_vars%grainc, &
             a_grainc, file_hist, 'f_grainc', itime_in_file, sumarea, filter, &
             'grain carbon display pool','gC/m2')

         ! grain carbon storage pool
         CALL write_history_variable_2d ( DEF_hist_vars%grainc_storage, &
             a_grainc_storage, file_hist, 'f_grainc_storage', itime_in_file, sumarea, filter, &
             'grain carbon storage pool','gC/m2')

         ! grain carbon transfer pool
         CALL write_history_variable_2d ( DEF_hist_vars%grainc_xfer, &
             a_grainc_xfer, file_hist, 'f_grainc_xfer', itime_in_file, sumarea, filter, &
             'grain carbon transfer pool','gC/m2')
#endif

         ! leaf nitrogen display pool
         CALL write_history_variable_2d ( DEF_hist_vars%leafn, &
             a_leafn, file_hist, 'f_leafn', itime_in_file, sumarea, filter, &
             'leaf nitrogen display pool','gN/m2')

         ! leaf nitrogen storage pool
         CALL write_history_variable_2d ( DEF_hist_vars%leafn_storage, &
             a_leafn_storage, file_hist, 'f_leafn_storage', itime_in_file, sumarea, filter, &
             'leaf nitrogen storage pool','gN/m2')

         ! leaf nitrogen transfer pool
         CALL write_history_variable_2d ( DEF_hist_vars%leafn_xfer, &
             a_leafn_xfer, file_hist, 'f_leafn_xfer', itime_in_file, sumarea, filter, &
             'leaf nitrogen transfer pool','gN/m2')

         ! fine root nitrogen display pool
         CALL write_history_variable_2d ( DEF_hist_vars%frootn, &
             a_frootn, file_hist, 'f_frootn', itime_in_file, sumarea, filter, &
             'fine root nitrogen display pool','gN/m2')

         ! fine root nitrogen storage pool
         CALL write_history_variable_2d ( DEF_hist_vars%frootn_storage, &
             a_frootn_storage, file_hist, 'f_frootn_storage', itime_in_file, sumarea, filter, &
             'fine root nitrogen storage pool','gN/m2')

         ! fine root nitrogen transfer pool
         CALL write_history_variable_2d ( DEF_hist_vars%frootn_xfer, &
             a_frootn_xfer, file_hist, 'f_frootn_xfer', itime_in_file, sumarea, filter, &
             'fine root nitrogen transfer pool','gN/m2')

         ! live stem nitrogen display pool
         CALL write_history_variable_2d ( DEF_hist_vars%livestemn, &
             a_livestemn, file_hist, 'f_livestemn', itime_in_file, sumarea, filter, &
             'live stem nitrogen display pool','gN/m2')

         ! live stem nitrogen storage pool
         CALL write_history_variable_2d ( DEF_hist_vars%livestemn_storage, &
             a_livestemn_storage, file_hist, 'f_livestemn_storage', &
             itime_in_file, sumarea, filter, &
             'live stem nitrogen storage pool','gN/m2')

         ! live stem nitrogen transfer pool
         CALL write_history_variable_2d ( DEF_hist_vars%livestemn_xfer, &
             a_livestemn_xfer, file_hist, 'f_livestemn_xfer', itime_in_file, sumarea, filter, &
             'live stem nitrogen transfer pool','gN/m2')

         ! dead stem nitrogen display pool
         CALL write_history_variable_2d ( DEF_hist_vars%deadstemn, &
             a_deadstemn, file_hist, 'f_deadstemn', itime_in_file, sumarea, filter, &
             'dead stem nitrogen display pool','gN/m2')

         ! dead stem nitrogen storage pool
         CALL write_history_variable_2d ( DEF_hist_vars%deadstemn_storage, &
             a_deadstemn_storage, file_hist, 'f_deadstemn_storage', &
             itime_in_file, sumarea, filter, &
             'dead stem nitrogen storage pool','gN/m2')

         ! dead stem nitrogen transfer pool
         CALL write_history_variable_2d ( DEF_hist_vars%deadstemn_xfer, &
             a_deadstemn_xfer, file_hist, 'f_deadstemn_xfer', itime_in_file, sumarea, filter, &
             'dead stem nitrogen transfer pool','gN/m2')

         ! live coarse root nitrogen display pool
         CALL write_history_variable_2d ( DEF_hist_vars%livecrootn, &
             a_livecrootn, file_hist, 'f_livecrootn', itime_in_file, sumarea, filter, &
             'live coarse root nitrogen display pool','gN/m2')

         ! live coarse root nitrogen storage pool
         CALL write_history_variable_2d ( DEF_hist_vars%livecrootn_storage, &
             a_livecrootn_storage, file_hist, 'f_livecrootn_storage', &
             itime_in_file, sumarea, filter, &
             'live coarse root nitrogen storage pool','gN/m2')

         ! live coarse root nitrogen transfer pool
         CALL write_history_variable_2d ( DEF_hist_vars%livecrootn_xfer, &
             a_livecrootn_xfer, file_hist, 'f_livecrootn_xfer', itime_in_file, sumarea, filter, &
             'live coarse root nitrogen transfer pool','gN/m2')

         ! dead coarse root nitrogen display pool
         CALL write_history_variable_2d ( DEF_hist_vars%deadcrootn, &
             a_deadcrootn, file_hist, 'f_deadcrootn', itime_in_file, sumarea, filter, &
             'dead coarse root nitrogen display pool','gN/m2')

         ! dead coarse root nitrogen storage pool
         CALL write_history_variable_2d ( DEF_hist_vars%deadcrootn_storage, &
             a_deadcrootn_storage, file_hist, 'f_deadcrootn_storage', &
             itime_in_file, sumarea, filter, &
             'dead coarse root nitrogen storage pool','gN/m2')

         ! dead coarse root nitrogen transfer pool
         CALL write_history_variable_2d ( DEF_hist_vars%deadcrootn_xfer, &
             a_deadcrootn_xfer, file_hist, 'f_deadcrootn_xfer', itime_in_file, sumarea, filter, &
             'dead coarse root nitrogen transfer pool','gN/m2')

#ifdef CROP
         ! grain nitrogen display pool
         CALL write_history_variable_2d ( DEF_hist_vars%grainn, &
             a_grainn, file_hist, 'f_grainn', itime_in_file, sumarea, filter, &
             'grain nitrogen display pool','gN/m2')

         ! grain nitrogen storage pool
         CALL write_history_variable_2d ( DEF_hist_vars%grainn_storage, &
             a_grainn_storage, file_hist, 'f_grainn_storage', itime_in_file, sumarea, filter, &
             'grain nitrogen storage pool','gN/m2')

         ! grain nitrogen transfer pool
         CALL write_history_variable_2d ( DEF_hist_vars%grainn_xfer, &
             a_grainn_xfer, file_hist, 'f_grainn_xfer', itime_in_file, sumarea, filter, &
             'grain nitrogen transfer pool','gN/m2')
#endif

         ! retranslocation nitrogen pool
         CALL write_history_variable_2d ( DEF_hist_vars%retrasn, &
             a_retransn, file_hist, 'f_retrasn', itime_in_file, sumarea, filter, &
             'retranslocation nitrogen pool','gN/m2')

         ! gross primary productivity
         CALL write_history_variable_2d ( DEF_hist_vars%gpp, &
             a_gpp, file_hist, 'f_gpp', itime_in_file, sumarea, filter, &
             'gross primary productivity','gC/m2/s')

         ! gross primary productivity
         CALL write_history_variable_2d ( DEF_hist_vars%downreg, &
             a_downreg, file_hist, 'f_downreg', itime_in_file, sumarea, filter, &
             'gpp downregulation due to N limitation','unitless')

         CALL write_history_variable_2d ( DEF_hist_vars%fpg, &
             a_fpg, file_hist, 'f_fpg', itime_in_file, sumarea, filter, &
             'fraction of gpp potential','unitless')

         CALL write_history_variable_2d ( DEF_hist_vars%fpi, &
             a_fpi, file_hist, 'f_fpi', itime_in_file, sumarea, filter, &
             'fraction of immobalization','unitless')

         CALL write_history_variable_2d ( DEF_hist_vars%totvegc, &
             a_totvegc, file_hist, 'f_totvegc', itime_in_file, sumarea, filter, &
             'total vegetation carbon','gC m-2')

         CALL write_history_variable_2d ( DEF_hist_vars%totlitc, &
             a_totlitc, file_hist, 'f_totlitc', itime_in_file, sumarea, filter, &
             'total litter carbon','gC m-2')

         CALL write_history_variable_2d ( DEF_hist_vars%totsomc, &
             a_totsomc, file_hist, 'f_totsomc', itime_in_file, sumarea, filter, &
             'total soil organic carbon','gC m-2')

         CALL write_history_variable_2d ( DEF_hist_vars%totcwdc, &
             a_totcwdc, file_hist, 'f_totcwdc', itime_in_file, sumarea, filter, &
             'total coarse woody debris carbon','gC m-2')

         CALL write_history_variable_2d ( DEF_hist_vars%totcolc, &
             a_totcolc, file_hist, 'f_totcolc', itime_in_file, sumarea, filter, &
             'total ecosystem carbon','gC m-2')

         CALL write_history_variable_2d ( DEF_hist_vars%totvegn, &
             a_totvegn, file_hist, 'f_totvegn', itime_in_file, sumarea, filter, &
             'total vegetation nitrogen','gN m-2')

         CALL write_history_variable_2d ( DEF_hist_vars%totlitn, &
             a_totlitn, file_hist, 'f_totlitn', itime_in_file, sumarea, filter, &
             'total litter nitrogen','gN m-2')

         CALL write_history_variable_2d ( DEF_hist_vars%totsomn, &
             a_totsomn, file_hist, 'f_totsomn', itime_in_file, sumarea, filter, &
             'total soil organic nitrogen','gN m-2')

         CALL write_history_variable_2d ( DEF_hist_vars%totcwdn, &
             a_totcwdn, file_hist, 'f_totcwdn', itime_in_file, sumarea, filter, &
             'total coarse woody debris nitrogen','gN m-2')

         CALL write_history_variable_2d ( DEF_hist_vars%totcoln, &
             a_totcoln, file_hist, 'f_totcoln', itime_in_file, sumarea, filter, &
             'total ecosystem nitrogen','gN m-2')

         ! autotrophic respiration
         CALL write_history_variable_2d ( DEF_hist_vars%ar , &
             a_ar, file_hist, 'f_ar', itime_in_file, sumarea, filter, &
             'autotrophic respiration','gC/m2/s')

!         ! CWD production
!         CALL write_history_variable_2d ( DEF_hist_vars%cwdprod , &
!             a_cwdprod, file_hist, 'f_cwdprod', itime_in_file, sumarea, filter, &
!             'CWD production','gC/m2/s')
!
!         ! CWD decomposition
!         CALL write_history_variable_2d ( DEF_hist_vars%cwddecomp , &
!             a_cwddecomp, file_hist, 'f_cwddecomp', itime_in_file, sumarea, filter, &
!             'CWD decomposition','gC/m2/s')

         ! heterotrophic respiration
         CALL write_history_variable_2d ( DEF_hist_vars%hr , &
             a_hr, file_hist, 'f_hr', itime_in_file, sumarea, filter, &
             'heterotrophic respiration','gC/m2/s')

#ifdef CROP
         ! crop phase
         CALL write_history_variable_2d ( DEF_hist_vars%cphase, &
             a_cphase, file_hist, 'f_cphase', itime_in_file, sumarea, filter, &
             'crop phase','unitless')

        ! heat unit index
        IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_hui (:)
            ENDIF
         ENDIF

        CALL write_history_variable_2d ( DEF_hist_vars%hui, &
             vecacc, file_hist, 'f_hui', itime_in_file, sumarea, filter, &
             'heat unit index','unitless')

         ! gdd needed to harvest
        CALL write_history_variable_2d ( DEF_hist_vars%gddmaturity, &
             a_gddmaturity, file_hist, 'f_gddmaturity', itime_in_file, sumarea, filter, &
             'gdd needed to harvest','ddays')

        ! gdd past planting date for crop
        CALL write_history_variable_2d ( DEF_hist_vars%gddplant, &
             a_gddplant, file_hist, 'f_gddplant', itime_in_file, sumarea, filter, &
             'gdd past planting date for crop','ddays')

        ! vernalization response
        CALL write_history_variable_2d ( DEF_hist_vars%vf, &
             a_vf, file_hist, 'f_vf', itime_in_file, sumarea, filter, &
             'vernalization response', 'unitless')

         ! 1-yr crop production carbon
         CALL write_history_variable_2d ( DEF_hist_vars%cropprod1c, &
             a_cropprod1c, file_hist, 'f_cropprod1c', itime_in_file, sumarea, filter, &
             '1-yr crop production carbon','gC/m2')

         ! loss rate of 1-yr crop production carbon
         CALL write_history_variable_2d ( DEF_hist_vars%cropprod1c_loss, &
             a_cropprod1c_loss, file_hist, 'f_cropprod1c_loss', itime_in_file, sumarea, filter, &
             'loss rate of 1-yr crop production carbon','gC/m2/s')

         ! crop seed deficit
         CALL write_history_variable_2d ( DEF_hist_vars%cropseedc_deficit, &
             a_cropseedc_deficit, file_hist, 'f_cropseedc_deficit', &
             itime_in_file, sumarea, filter, &
             'crop seed deficit','gC/m2/s')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_grainc_to_cropprodc (:)
            ENDIF
         ENDIF
         ! grain to crop production carbon
         CALL write_history_variable_2d ( DEF_hist_vars%grainc_to_cropprodc, &
             vecacc, file_hist, 'f_grainc_to_cropprodc', itime_in_file, sumarea, filter, &
             'grain to crop production carbon','gC/m2/s')

         ! grain to crop seed carbon
         CALL write_history_variable_2d ( DEF_hist_vars%grainc_to_seed, &
             a_grainc_to_seed, file_hist, 'f_grainc_to_seed', itime_in_file, sumarea, filter, &
             'grain to crop seed carbon','gC/m2/s')
         ! grain to crop seed carbon
         CALL write_history_variable_2d ( DEF_hist_vars%fert_to_sminn, &
             a_fert_to_sminn, file_hist, 'f_fert_to_sminn', itime_in_file, sumarea, filter, &
             'fertilization','gN/m2/s')

         IF (DEF_USE_IRRIGATION) THEN
            IF (p_is_worker) THEN
               IF (numpatch > 0) THEN
                  DO i=1,numpatch
                     IF(patchclass(i) == 12)THEN
                        IF((pftclass(patch_pft_s(i)).GE.npcropmin).and.(MOD(pftclass(patch_pft_s(i)),2).EQ.0))THEN
                           filter_irrig(i) = .true.
                        ELSE
                           filter_irrig(i) = .false.
                        ENDIF
                     ELSE
                        filter_irrig(i) = .false.
                     ENDIF
                  ENDDO
                  IF (DEF_forcing%has_missing_value) THEN
                     filter_irrig = filter_irrig .and. forcmask_pch
                  ENDIF
                  filter_irrig = filter_irrig .and. patchmask
               ENDIF
            ENDIF

            IF (HistForm == 'Gridded') THEN
               CALL mp2g_hist%get_sumarea (sumarea_irrig, filter_irrig)
            ENDIF

            !  total irrigation amounts at growing season
            IF (p_is_worker) THEN
               IF (numpatch > 0) THEN
                  vecacc (:) = a_sum_irrig (:)
               ENDIF
            ENDIF
            CALL write_history_variable_2d ( DEF_hist_vars%sum_irrig, &
               vecacc, file_hist, 'f_sum_irrig', itime_in_file, sumarea_irrig, filter_irrig, &
               'total irrigation amounts at growing season','kg/m2')

            !  total irrigation amounts demand at growing season
            IF (p_is_worker) THEN
               IF (numpatch > 0) THEN
                  vecacc (:) = a_sum_deficit_irrig (:)
               ENDIF
            ENDIF
            CALL write_history_variable_2d ( DEF_hist_vars%sum_deficit_irrig, &
               vecacc, file_hist, 'f_sum_deficit_irrig', itime_in_file, sumarea_irrig, filter_irrig, &
               'total irrigation amounts demand at growing season','kg/m2')

            ! total irrigation times at growing season
            CALL write_history_variable_2d ( DEF_hist_vars%sum_irrig_count, &
               a_sum_irrig_count, file_hist, 'f_sum_irrig_count', itime_in_file, sumarea_irrig, filter_irrig, &
               'total irrigation times at growing season','-')

            ! irrigation waterstorage [kg/m2]
            IF (p_is_worker) THEN
               IF (numpatch > 0) THEN
                  vecacc (:) = a_waterstorage (:)
               ENDIF
            ENDIF
            CALL write_history_variable_2d ( DEF_hist_vars%waterstorage, &
               vecacc, file_hist, 'f_waterstorage', itime_in_file, sumarea_irrig, filter_irrig, &
               'irrigation waterstorage','kg/m2')

            ! irrigation demand for ground water [kg/m2]
            IF (p_is_worker) THEN
               IF (numpatch > 0) THEN
                  vecacc (:) = a_groundwater_demand (:)
               ENDIF
            ENDIF
            CALL write_history_variable_2d ( DEF_hist_vars%groundwater_demand, &
               vecacc, file_hist, 'f_groundwater_demand', itime_in_file, sumarea_irrig, filter_irrig, &
               'irrigation demand for ground water','kg/m2')

            ! irrigation supply from ground water [kg/m2]
            IF (p_is_worker) THEN
               IF (numpatch > 0) THEN
                  vecacc (:) = a_groundwater_supply (:)
               ENDIF
            ENDIF
            CALL write_history_variable_2d ( DEF_hist_vars%groundwater_supply, &
               vecacc, file_hist, 'f_groundwater_supply', itime_in_file, sumarea_irrig, filter_irrig, &
               'irrigation supply from ground water','kg/m2')

            ! irrigation demand for reservoir or river [kg/m2]
            IF (p_is_worker) THEN
               IF (numpatch > 0) THEN
                  vecacc (:) = a_reservoirriver_demand (:)
               ENDIF
            ENDIF
            CALL write_history_variable_2d ( DEF_hist_vars%reservoirriver_demand, &
               vecacc, file_hist, 'f_reservoirriver_demand', itime_in_file, sumarea_irrig, filter_irrig, &
               'irrigation demand for reservoir or river','kg/m2')

            ! irrigation supply from reservoir or river [kg/m2]
            IF (p_is_worker) THEN
               IF (numpatch > 0) THEN
                  vecacc (:) = a_reservoirriver_supply (:)
               ENDIF
            ENDIF
            CALL write_history_variable_2d ( DEF_hist_vars%reservoirriver_supply, &
               vecacc, file_hist, 'f_reservoirriver_supply', itime_in_file, sumarea_irrig, filter_irrig, &
               'irrigation supply from reservoir or river','kg/m2')

            ! irrigation supply from reservoir [kg/m2]
            IF (p_is_worker) THEN
               IF (numpatch > 0) THEN
                  vecacc (:) = a_reservoir_supply (:)
               ENDIF
            ENDIF
            CALL write_history_variable_2d ( DEF_hist_vars%reservoirriver_supply, &
               vecacc, file_hist, 'f_reservoir_supply', itime_in_file, sumarea_irrig, filter_irrig, &
               'irrigation supply from reservoir','kg/m2')

            ! irrigation supply from river [kg/m2]
            IF (p_is_worker) THEN
               IF (numpatch > 0) THEN
                  vecacc (:) = a_river_supply (:)
               ENDIF
            ENDIF
            CALL write_history_variable_2d ( DEF_hist_vars%reservoirriver_supply, &
               vecacc, file_hist, 'f_river_supply', itime_in_file, sumarea_irrig, filter_irrig, &
               'irrigation supply from river','kg/m2')

            ! irrigation supply from runoff [kg/m2]
            IF (p_is_worker) THEN
               IF (numpatch > 0) THEN
                  vecacc (:) = a_runoff_supply (:)
               ENDIF
            ENDIF
            CALL write_history_variable_2d ( DEF_hist_vars%reservoirriver_supply, &
               vecacc, file_hist, 'f_runoff_supply', itime_in_file, sumarea_irrig, filter_irrig, &
               'irrigation supply from runoff','kg/m2')
         ENDIF
#endif

         ! grain to crop seed carbon
         CALL write_history_variable_2d ( DEF_hist_vars%ndep_to_sminn, &
             a_ndep_to_sminn, file_hist, 'f_ndep_to_sminn', itime_in_file, sumarea, filter, &
             'nitrogen deposition','gN/m2/s')

         IF(DEF_USE_DiagMatrix)THEN
            ! leaf carbon display pool
            CALL write_history_variable_2d ( DEF_hist_vars%leafcCap, &
                a_leafcCap, file_hist, 'f_leafcCap', itime_in_file, sumarea, filter, &
                'leaf carbon display pool Capacity','gC/m2')

            ! leaf carbon storage pool
            CALL write_history_variable_2d ( DEF_hist_vars%leafc_storageCap, &
                a_leafc_storageCap, file_hist, 'f_leafc_storageCap', &
                itime_in_file, sumarea, filter, &
                'leaf carbon storage pool capacity','gC/m2')

            ! leaf carbon transfer pool
            CALL write_history_variable_2d ( DEF_hist_vars%leafc_xferCap, &
                a_leafc_xferCap, file_hist, 'f_leafc_xferCap', itime_in_file, sumarea, filter, &
                'leaf carbon transfer pool capacity','gC/m2')

            ! fine root carbon display pool
            CALL write_history_variable_2d ( DEF_hist_vars%frootcCap, &
                a_frootcCap, file_hist, 'f_frootcCap', itime_in_file, sumarea, filter, &
                'fine root carbon display pool Capacity','gC/m2')

            ! fine root carbon storage pool
            CALL write_history_variable_2d ( DEF_hist_vars%frootc_storageCap, &
                a_frootc_storageCap, file_hist, 'f_frootc_storageCap', &
                itime_in_file, sumarea, filter, &
                'fine root carbon storage pool capacity','gC/m2')

            ! fine root carbon transfer pool
            CALL write_history_variable_2d ( DEF_hist_vars%frootc_xferCap, &
                a_frootc_xferCap, file_hist, 'f_frootc_xferCap', itime_in_file, sumarea, filter, &
                'fine root carbon transfer pool capacity','gC/m2')

            ! live stem carbon display pool
            CALL write_history_variable_2d ( DEF_hist_vars%livestemcCap, &
                a_livestemcCap, file_hist, 'f_livestemcCap', itime_in_file, sumarea, filter, &
                'live stem carbon display pool Capacity','gC/m2')

            ! live stem carbon storage pool
            CALL write_history_variable_2d ( DEF_hist_vars%livestemc_storageCap, &
                a_livestemc_storageCap, file_hist, 'f_livestemc_storageCap', &
                itime_in_file, sumarea, filter, &
                'live stem carbon storage pool capacity','gC/m2')

            ! live stem carbon transfer pool
            CALL write_history_variable_2d ( DEF_hist_vars%livestemc_xferCap, &
                a_livestemc_xferCap, file_hist, 'f_livestemc_xferCap', &
                itime_in_file, sumarea, filter, &
                'live stem carbon transfer pool capacity','gC/m2')

            ! dead stem carbon display pool
            CALL write_history_variable_2d ( DEF_hist_vars%deadstemcCap, &
                a_deadstemcCap, file_hist, 'f_deadstemcCap', itime_in_file, sumarea, filter, &
                'dead stem carbon display pool Capacity','gC/m2')

            ! dead stem carbon storage pool
            CALL write_history_variable_2d ( DEF_hist_vars%deadstemc_storageCap, &
                a_deadstemc_storageCap, file_hist, 'f_deadstemc_storageCap', &
                itime_in_file, sumarea, filter, &
                'dead stem carbon storage pool capacity','gC/m2')

            ! dead stem carbon transfer pool
            CALL write_history_variable_2d ( DEF_hist_vars%deadstemc_xferCap, &
                a_deadstemc_xferCap, file_hist, 'f_deadstemc_xferCap', &
                itime_in_file, sumarea, filter, &
                'dead stem carbon transfer pool capacity','gC/m2')

            ! live coarse root carbon display pool
            CALL write_history_variable_2d ( DEF_hist_vars%livecrootcCap, &
                a_livecrootcCap, file_hist, 'f_livecrootcCap', itime_in_file, sumarea, filter, &
                'live coarse root carbon display pool Capacity','gC/m2')

            ! live coarse root carbon storage pool
            CALL write_history_variable_2d ( DEF_hist_vars%livecrootc_storageCap, &
                a_livecrootc_storageCap, file_hist, 'f_livecrootc_storageCap', &
                itime_in_file, sumarea, filter, &
                'live coarse root carbon storage pool capacity','gC/m2')

            ! live coarse root carbon transfer pool
            CALL write_history_variable_2d ( DEF_hist_vars%livecrootc_xferCap, &
                a_livecrootc_xferCap, file_hist, 'f_livecrootc_xferCap', &
                itime_in_file, sumarea, filter, &
                'live coarse root carbon transfer pool capacity','gC/m2')

            ! dead coarse root carbon display pool
            CALL write_history_variable_2d ( DEF_hist_vars%deadcrootcCap, &
                a_deadcrootcCap, file_hist, 'f_deadcrootcCap', itime_in_file, sumarea, filter, &
                'dead coarse root carbon display pool Capacity','gC/m2')

            ! dead coarse root carbon storage pool
            CALL write_history_variable_2d ( DEF_hist_vars%deadcrootc_storageCap, &
                a_deadcrootc_storageCap, file_hist, 'f_deadcrootc_storageCap', &
                itime_in_file, sumarea, filter, &
                'dead coarse root carbon storage pool capacity','gC/m2')

            ! dead coarse root carbon transfer pool
            CALL write_history_variable_2d ( DEF_hist_vars%deadcrootc_xferCap, &
                a_deadcrootc_xferCap, file_hist, 'f_deadcrootc_xferCap', &
                itime_in_file, sumarea, filter, &
                'dead coarse root carbon transfer pool capacity','gC/m2')

            ! leaf nitrogen display pool
            CALL write_history_variable_2d ( DEF_hist_vars%leafnCap, &
                a_leafnCap, file_hist, 'f_leafnCap', itime_in_file, sumarea, filter, &
                'leaf nitrogen display pool Capacity','gC/m2')

            ! leaf nitrogen storage pool
            CALL write_history_variable_2d ( DEF_hist_vars%leafn_storageCap, &
                a_leafn_storageCap, file_hist, 'f_leafn_storageCap', &
                itime_in_file, sumarea, filter, &
                'leaf nitrogen storage pool capacity','gC/m2')

            ! leaf nitrogen transfer pool
            CALL write_history_variable_2d ( DEF_hist_vars%leafn_xferCap, &
                a_leafn_xferCap, file_hist, 'f_leafn_xferCap', itime_in_file, sumarea, filter, &
                'leaf nitrogen transfer pool capacity','gC/m2')

            ! fine root nitrogen display pool
            CALL write_history_variable_2d ( DEF_hist_vars%frootnCap, &
                a_frootnCap, file_hist, 'f_frootnCap', itime_in_file, sumarea, filter, &
                'fine root nitrogen display pool Capacity','gC/m2')

            ! fine root nitrogen storage pool
            CALL write_history_variable_2d ( DEF_hist_vars%frootn_storageCap, &
                a_frootn_storageCap, file_hist, 'f_frootn_storageCap', &
                itime_in_file, sumarea, filter, &
                'fine root nitrogen storage pool capacity','gC/m2')

            ! fine root nitrogen transfer pool
            CALL write_history_variable_2d ( DEF_hist_vars%frootn_xferCap, &
                a_frootn_xferCap, file_hist, 'f_frootn_xferCap', itime_in_file, sumarea, filter, &
                'fine root nitrogen transfer pool capacity','gC/m2')

            ! live stem nitrogen display pool
            CALL write_history_variable_2d ( DEF_hist_vars%livestemnCap, &
                a_livestemnCap, file_hist, 'f_livestemnCap', itime_in_file, sumarea, filter, &
                'live stem nitrogen display pool Capacity','gC/m2')

            ! live stem nitrogen storage pool
            CALL write_history_variable_2d ( DEF_hist_vars%livestemn_storageCap, &
                a_livestemn_storageCap, file_hist, 'f_livestemn_storageCap', &
                itime_in_file, sumarea, filter, &
                'live stem nitrogen storage pool capacity','gC/m2')

            ! live stem nitrogen transfer pool
            CALL write_history_variable_2d ( DEF_hist_vars%livestemn_xferCap, &
                a_livestemn_xferCap, file_hist, 'f_livestemn_xferCap', &
                itime_in_file, sumarea, filter, &
                'live stem nitrogen transfer pool capacity','gC/m2')

            ! dead stem nitrogen display pool
            CALL write_history_variable_2d ( DEF_hist_vars%deadstemnCap, &
                a_deadstemnCap, file_hist, 'f_deadstemnCap', itime_in_file, sumarea, filter, &
                'dead stem nitrogen display pool Capacity','gC/m2')

            ! dead stem nitrogen storage pool
            CALL write_history_variable_2d ( DEF_hist_vars%deadstemn_storageCap, &
                a_deadstemn_storageCap, file_hist, 'f_deadstemn_storageCap', &
                itime_in_file, sumarea, filter, &
                'dead stem nitrogen storage pool capacity','gC/m2')

            ! dead stem nitrogen transfer pool
            CALL write_history_variable_2d ( DEF_hist_vars%deadstemn_xferCap, &
                a_deadstemn_xferCap, file_hist, 'f_deadstemn_xferCap', &
                itime_in_file, sumarea, filter, &
                'dead stem nitrogen transfer pool capacity','gC/m2')

            ! live coarse root nitrogen display pool
            CALL write_history_variable_2d ( DEF_hist_vars%livecrootnCap, &
                a_livecrootnCap, file_hist, 'f_livecrootnCap', itime_in_file, sumarea, filter, &
                'live coarse root nitrogen display pool Capacity','gC/m2')

            ! live coarse root nitrogen storage pool
            CALL write_history_variable_2d ( DEF_hist_vars%livecrootn_storageCap, &
                a_livecrootn_storageCap, file_hist, 'f_livecrootn_storageCap', &
                itime_in_file, sumarea, filter, &
                'live coarse root nitrogen storage pool capacity','gC/m2')

            ! live coarse root nitrogen transfer pool
            CALL write_history_variable_2d ( DEF_hist_vars%livecrootn_xferCap, &
                a_livecrootn_xferCap, file_hist, 'f_livecrootn_xferCap', &
                itime_in_file, sumarea, filter, &
                'live coarse root nitrogen transfer pool capacity','gC/m2')

            ! dead coarse root nitrogen display pool
            CALL write_history_variable_2d ( DEF_hist_vars%deadcrootnCap, &
                a_deadcrootnCap, file_hist, 'f_deadcrootnCap', itime_in_file, sumarea, filter, &
                'dead coarse root nitrogen display pool Capacity','gC/m2')

            ! dead coarse root nitrogen storage pool
            CALL write_history_variable_2d ( DEF_hist_vars%deadcrootn_storageCap, &
                a_deadcrootn_storageCap, file_hist, 'f_deadcrootn_storageCap', &
                itime_in_file, sumarea, filter, &
                'dead coarse root nitrogen storage pool capacity','gC/m2')

            ! dead coarse root nitrogen transfer pool
            CALL write_history_variable_2d ( DEF_hist_vars%deadcrootn_xferCap, &
                a_deadcrootn_xferCap, file_hist, 'f_deadcrootn_xferCap', &
                itime_in_file, sumarea, filter, &
                'dead coarse root nitrogen transfer pool capacity','gC/m2')

         ENDIF

         IF(DEF_USE_OZONESTRESS)THEN
         ! ozone concentration
            CALL write_history_variable_2d ( DEF_hist_vars%xy_ozone, &
               a_ozone, file_hist, 'f_xy_ozone', itime_in_file, sumarea, filter, &
               'Ozone concentration','mol/mol')
         ENDIF

         ! litter 1 carbon density in soil layers
         CALL write_history_variable_3d ( DEF_hist_vars%litr1c_vr, &
            a_litr1c_vr, file_hist, 'f_litr1c_vr', itime_in_file, 'soil', 1, nl_soil, &
            sumarea, filter,'litter 1 carbon density in soil layers','gC/m3')

         ! litter 2 carbon density in soil layers
         CALL write_history_variable_3d ( DEF_hist_vars%litr2c_vr, &
            a_litr2c_vr, file_hist, 'f_litr2c_vr', itime_in_file, 'soil', 1, nl_soil, &
            sumarea, filter,'litter 2 carbon density in soil layers','gC/m3')

         ! litter 3 carbon density in soil layers
         CALL write_history_variable_3d ( DEF_hist_vars%litr3c_vr, &
            a_litr3c_vr, file_hist, 'f_litr3c_vr', itime_in_file, 'soil', 1, nl_soil, &
            sumarea, filter,'litter 3 carbon density in soil layers','gC/m3')

         ! soil 1 carbon density in soil layers
         CALL write_history_variable_3d ( DEF_hist_vars%soil1c_vr, &
            a_soil1c_vr, file_hist, 'f_soil1c_vr', itime_in_file, 'soil', 1, nl_soil, &
            sumarea, filter,'soil 1 carbon density in soil layers','gC/m3')

         ! soil 2 carbon density in soil layers
         CALL write_history_variable_3d ( DEF_hist_vars%soil2c_vr, &
            a_soil2c_vr, file_hist, 'f_soil2c_vr', itime_in_file, 'soil', 1, nl_soil, &
            sumarea, filter,'soil 2 carbon density in soil layers','gC/m3')

         ! soil 3 carbon density in soil layers
         CALL write_history_variable_3d ( DEF_hist_vars%soil3c_vr, &
            a_soil3c_vr, file_hist, 'f_soil3c_vr', itime_in_file, 'soil', 1, nl_soil, &
            sumarea, filter,'soil 3 carbon density in soil layers','gC/m3')

         ! coarse woody debris carbon density in soil layers
         CALL write_history_variable_3d ( DEF_hist_vars%cwdc_vr, &
            a_cwdc_vr, file_hist, 'f_cwdc_vr', itime_in_file, 'soil', 1, nl_soil, &
            sumarea, filter,'coarse woody debris carbon density in soil layers','gC/m3')

         ! litter 1 nitrogen density in soil layers
         CALL write_history_variable_3d ( DEF_hist_vars%litr1n_vr, &
            a_litr1n_vr, file_hist, 'f_litr1n_vr', itime_in_file, 'soil', 1, nl_soil, &
            sumarea, filter,'litter 1 nitrogen density in soil layers','gN/m3')

         ! litter 2 nitrogen density in soil layers
         CALL write_history_variable_3d ( DEF_hist_vars%litr2n_vr, &
            a_litr2n_vr, file_hist, 'f_litr2n_vr', itime_in_file, 'soil', 1, nl_soil, &
            sumarea, filter,'litter 2 nitrogen density in soil layers','gN/m3')

         ! litter 3 nitrogen density in soil layers
         CALL write_history_variable_3d ( DEF_hist_vars%litr3n_vr, &
            a_litr3n_vr, file_hist, 'f_litr3n_vr', itime_in_file, 'soil', 1, nl_soil, &
            sumarea, filter,'litter 3 nitrogen density in soil layers','gN/m3')

         ! soil 1 nitrogen density in soil layers
         CALL write_history_variable_3d ( DEF_hist_vars%soil1n_vr, &
            a_soil1n_vr, file_hist, 'f_soil1n_vr', itime_in_file, 'soil', 1, nl_soil, &
            sumarea, filter,'soil 1 nitrogen density in soil layers','gN/m3')

         ! soil 2 nitrogen density in soil layers
         CALL write_history_variable_3d ( DEF_hist_vars%soil2n_vr, &
            a_soil2n_vr, file_hist, 'f_soil2n_vr', itime_in_file, 'soil', 1, nl_soil, &
            sumarea, filter,'soil 2 nitrogen density in soil layers','gN/m3')

         ! soil 3 nitrogen density in soil layers
         CALL write_history_variable_3d ( DEF_hist_vars%soil3n_vr, &
            a_soil3n_vr, file_hist, 'f_soil3n_vr', itime_in_file, 'soil', 1, nl_soil, &
            sumarea, filter,'soil 3 nitrogen density in soil layers','gN/m3')

         ! coarse woody debris nitrogen density in soil layers
         CALL write_history_variable_3d ( DEF_hist_vars%cwdn_vr, &
            a_cwdn_vr, file_hist, 'f_cwdn_vr', itime_in_file, 'soil', 1, nl_soil, &
            sumarea, filter,'coarse woody debris nitrogen density in soil layers','gN/m3')

         ! mineral nitrogen density in soil layers
         CALL write_history_variable_3d ( DEF_hist_vars%sminn_vr, &
            a_sminn_vr, file_hist, 'f_sminn_vr', itime_in_file, 'soil', 1, nl_soil, &
            sumarea, filter,'mineral nitrogen density in soil layers','gN/m3')

         ! total nitrogen percentage in soil layers
         CALL write_history_variable_3d ( DEF_hist_vars%totsoiln_vr, &
            a_totsoiln_vr, file_hist, 'f_totsoiln_vr', itime_in_file, 'soil', 1, nl_soil, &
            sumarea, filter,&
         'Total nitrogen in soil layers, percentage of total soil nitrogen in total soil mass','%')

         ! bulk density in soil layers
         CALL write_history_variable_3d ( DEF_hist_vars%BD_all, &
            a_BD_all, file_hist, 'f_BD_all', itime_in_file, 'soil', 1, nl_soil, &
            sumarea, filter,'bulk density in soil layers','kg/m3')

         ! field capacity in soil layers
         CALL write_history_variable_3d ( DEF_hist_vars%wfc, &
            a_wfc, file_hist, 'f_wfc', itime_in_file, 'soil', 1, nl_soil, &
            sumarea, filter,'field capacity in soil layers','m3/m3')

         ! organic matter density in soil layers
         CALL write_history_variable_3d ( DEF_hist_vars%OM_density, &
            a_OM_density, file_hist, 'f_OM_density', itime_in_file, 'soil', 1, nl_soil, &
            sumarea, filter,'organic matter density in soil layers','kg/m3')

         IF (DEF_USE_NITRIF) THEN
            ! O2 soil Concentration for non-inundated area
            CALL write_history_variable_3d ( DEF_hist_vars%CONC_O2_UNSAT, &
               a_conc_o2_unsat, file_hist, 'f_CONC_O2_UNSAT', itime_in_file, 'soil', 1, nl_soil, &
               sumarea, filter,'O2 soil Concentration for non-inundated area','mol/m3')

            ! O2 consumption from HR and AR for non-inundated area
            CALL write_history_variable_3d ( DEF_hist_vars%O2_DECOMP_DEPTH_UNSAT, &
               a_o2_decomp_depth_unsat, file_hist, 'f_O2_DECOMP_DEPTH_UNSAT', &
               itime_in_file, 'soil', 1, nl_soil, &
               sumarea, filter,'O2 consumption from HR and AR for non-inundated area','mol/m3/s')
         ENDIF

         IF (DEF_USE_FIRE) THEN
            CALL write_history_variable_2d ( DEF_hist_vars%abm, &
                 vecacc, file_hist, 'f_abm', itime_in_file, sumarea, filter, &
                 'peak crop fire month','unitless')

            CALL write_history_variable_2d ( DEF_hist_vars%gdp, &
                 vecacc, file_hist, 'f_gdp', itime_in_file, sumarea, filter, &
                 'gdp','unitless')

            CALL write_history_variable_2d ( DEF_hist_vars%peatf, &
                 vecacc, file_hist, 'f_peatf', itime_in_file, sumarea, filter, &
                 'peatf','unitless')

            CALL write_history_variable_2d ( DEF_hist_vars%hdm, &
                 vecacc, file_hist, 'f_hdm', itime_in_file, sumarea, filter, &
                 'hdm','unitless')

            CALL write_history_variable_2d ( DEF_hist_vars%lnfm, &
                 vecacc, file_hist, 'f_lnfm', itime_in_file, sumarea, filter, &
                 'lnfm','unitless')
         ENDIF

         IF(DEF_USE_DiagMatrix)THEN
            ! litter 1 carbon capacity density in soil layers
            CALL write_history_variable_3d ( DEF_hist_vars%litr1cCap_vr, &
               a_litr1cCap_vr, file_hist, 'f_litr1cCap_vr', itime_in_file, 'soil', 1, nl_soil, &
               sumarea, filter,'litter 1 carbon capacity density in soil layers','gC/m3')

            ! litter 2 carbon capacity density in soil layers
            CALL write_history_variable_3d ( DEF_hist_vars%litr2cCap_vr, &
               a_litr2cCap_vr, file_hist, 'f_litr2cCap_vr', itime_in_file, 'soil', 1, nl_soil, &
               sumarea, filter,'litter 2 carbon capacity density in soil layers','gC/m3')

            ! litter 3 carbon capacity density in soil layers
            CALL write_history_variable_3d ( DEF_hist_vars%litr3cCap_vr, &
               a_litr3cCap_vr, file_hist, 'f_litr3cCap_vr', itime_in_file, 'soil', 1, nl_soil, &
               sumarea, filter,'litter 3 carbon capacity density in soil layers','gC/m3')

            ! soil 1 carbon capacity density in soil layers
            CALL write_history_variable_3d ( DEF_hist_vars%soil1cCap_vr, &
               a_soil1cCap_vr, file_hist, 'f_soil1cCap_vr', itime_in_file, 'soil', 1, nl_soil, &
               sumarea, filter,'soil 1 carbon capacity density in soil layers','gC/m3')

            ! soil 2 carbon capacity density in soil layers
            CALL write_history_variable_3d ( DEF_hist_vars%soil2cCap_vr, &
               a_soil2cCap_vr, file_hist, 'f_soil2cCap_vr', itime_in_file, 'soil', 1, nl_soil, &
               sumarea, filter,'soil 2 carbon capacity density in soil layers','gC/m3')

            ! soil 3 carbon capacity density in soil layers
            CALL write_history_variable_3d ( DEF_hist_vars%soil3cCap_vr, &
               a_soil3cCap_vr, file_hist, 'f_soil3cCap_vr', itime_in_file, 'soil', 1, nl_soil, &
               sumarea, filter,'soil 3 carbon capacity density in soil layers','gC/m3')

            ! coarse woody debris carbon capacity density in soil layers
            CALL write_history_variable_3d ( DEF_hist_vars%cwdcCap_vr, &
               a_cwdcCap_vr, file_hist, 'f_cwdcCap_vr', itime_in_file, 'soil', 1, nl_soil, &
               sumarea, filter,'coarse woody debris carbon capacity density in soil layers','gC/m3')

            ! litter 1 nitrogen capacity density in soil layers
            CALL write_history_variable_3d ( DEF_hist_vars%litr1nCap_vr, &
               a_litr1nCap_vr, file_hist, 'f_litr1nCap_vr', itime_in_file, 'soil', 1, nl_soil, &
               sumarea, filter,'litter 1 nitrogen capacity density in soil layers','gN/m3')

            ! litter 2 nitrogen capacity density in soil layers
            CALL write_history_variable_3d ( DEF_hist_vars%litr2nCap_vr, &
               a_litr2nCap_vr, file_hist, 'f_litr2nCap_vr', itime_in_file, 'soil', 1, nl_soil, &
               sumarea, filter,'litter 2 nitrogen capacity density in soil layers','gN/m3')

            ! litter 3 nitrogen capacity density in soil layers
            CALL write_history_variable_3d ( DEF_hist_vars%litr3nCap_vr, &
               a_litr3nCap_vr, file_hist, 'f_litr3nCap_vr', itime_in_file, 'soil', 1, nl_soil, &
               sumarea, filter,'litter 3 nitrogen capacity density in soil layers','gN/m3')

            ! soil 1 nitrogen capacity density in soil layers
            CALL write_history_variable_3d ( DEF_hist_vars%soil1nCap_vr, &
               a_soil1nCap_vr, file_hist, 'f_soil1nCap_vr', itime_in_file, 'soil', 1, nl_soil, &
               sumarea, filter,'soil 1 nitrogen capacity density in soil layers','gN/m3')

            ! soil 2 nitrogen capacity density in soil layers
            CALL write_history_variable_3d ( DEF_hist_vars%soil2nCap_vr, &
               a_soil2nCap_vr, file_hist, 'f_soil2nCap_vr', itime_in_file, 'soil', 1, nl_soil, &
               sumarea, filter,'soil 2 nitrogen capacity density in soil layers','gN/m3')

            ! soil 3 nitrogen capacity density in soil layers
            CALL write_history_variable_3d ( DEF_hist_vars%soil3nCap_vr, &
               a_soil3nCap_vr, file_hist, 'f_soil3nCap_vr', itime_in_file, 'soil', 1, nl_soil, &
               sumarea, filter,'soil 3 nitrogen capacity density in soil layers','gN/m3')

            ! coarse woody debris nitrogen capacity density in soil layers
            CALL write_history_variable_3d ( DEF_hist_vars%cwdnCap_vr, &
               a_cwdnCap_vr, file_hist, 'f_cwdnCap_vr', itime_in_file, 'soil', 1, nl_soil, &
               sumarea, filter,&
               'coarse woody debris nitrogen capacity density in soil layers','gN/m3')

            ! Temperature environmental scalar
            CALL write_history_variable_3d ( DEF_hist_vars%t_scalar, &
                a_t_scalar, file_hist, 'f_t_scalar', itime_in_file, 'soil', 1, nl_soil, &
                sumarea, filter, 'Temperature environmental scalar','unitless')

            ! Water environmental scalar
            CALL write_history_variable_3d ( DEF_hist_vars%w_scalar, &
                a_w_scalar, file_hist, 'f_w_scalar', itime_in_file, 'soil', 1, nl_soil, &
                sumarea, filter, 'Water environmental scalar','unitless')

         ENDIF



         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) .ne. 12 .and. patchtype(i) .eq. 0)THEN
                     filter(i) = .true.
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! 1: gpp enf temperate
         CALL write_history_variable_2d ( DEF_hist_vars%gpp_enftemp, &
             a_gpp_enftemp, file_hist, 'f_gpp_enftemp', itime_in_file, sumarea, filter, &
             'gross primary productivity for needleleaf evergreen temperate tree','gC/m2/s')

         ! 1: leaf carbon display pool enf temperate
         CALL write_history_variable_2d ( DEF_hist_vars%leafc_enftemp, &
             a_leafc_enftemp, file_hist, 'f_leafc_enftemp', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for needleleaf evergreen temperate tree','gC/m2')

         ! 1: leaf area index enf temperate
         CALL write_history_variable_2d ( DEF_hist_vars%lai_enftemp, &
             a_lai_enftemp, file_hist, 'f_lai_enftemp', itime_in_file, sumarea, filter, &
             'leaf area index for needleleaf evergreen temperate tree','m2/m2')

         ! 1: npp enf temperate
         CALL write_history_variable_2d ( DEF_hist_vars%npp_enftemp, &
             a_npp_enftemp, file_hist, 'f_npp_enftemp', itime_in_file, sumarea, filter, &
             'npp for needleleaf evergreen temperate tree','m2/m2')

         ! 1: npp to leafc enf temperate
         CALL write_history_variable_2d ( DEF_hist_vars%npptoleafc_enftemp, &
             a_npptoleafc_enftemp, file_hist, 'f_npptoleafc_enftemp', itime_in_file, sumarea, filter, &
             'npp to leafc for needleleaf evergreen temperate tree','m2/m2')

         ! 2: gpp enf boreal
         CALL write_history_variable_2d ( DEF_hist_vars%gpp_enfboreal, &
             a_gpp_enfboreal, file_hist, 'f_gpp_enfboreal', itime_in_file, sumarea, filter, &
             'gross primary productivity for needleleaf evergreen boreal tree','gC/m2/s')

         ! 2: leaf carbon display pool enf boreal
         CALL write_history_variable_2d ( DEF_hist_vars%leafc_enfboreal, &
             a_leafc_enfboreal, file_hist, 'f_leafc_enfboreal', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for needleleaf evergreen boreal tree','gC/m2')

         ! 2: leaf area index enf boreal
         CALL write_history_variable_2d ( DEF_hist_vars%lai_enfboreal, &
             a_lai_enfboreal, file_hist, 'f_lai_enfboreal', itime_in_file, sumarea, filter, &
             'leaf area index for needleleaf evergreen boreal tree','m2/m2')

         ! 2: npp enf boreal
         CALL write_history_variable_2d ( DEF_hist_vars%npp_enfboreal, &
             a_npp_enfboreal, file_hist, 'f_npp_enfboreal', itime_in_file, sumarea, filter, &
             'npp for needleleaf evergreen boreal tree','m2/m2')

         ! 2: npp to leafc enf boreal
         CALL write_history_variable_2d ( DEF_hist_vars%npptoleafc_enfboreal, &
             a_npptoleafc_enfboreal, file_hist, 'f_npptoleafc_enfboreal', itime_in_file, sumarea, filter, &
             'npp to leafc for needleleaf evergreen boreal tree','m2/m2')

         ! 3: gpp dnf boreal
         CALL write_history_variable_2d ( DEF_hist_vars%gpp_dnfboreal, &
             a_gpp_dnfboreal, file_hist, 'f_gpp_dnfboreal', itime_in_file, sumarea, filter, &
             'gross primary productivity for needleleaf deciduous boreal tree','gC/m2/s')

         ! 3: leaf carbon display pool dnf boreal
         CALL write_history_variable_2d ( DEF_hist_vars%leafc_dnfboreal, &
             a_leafc_dnfboreal, file_hist, 'f_leafc_dnfboreal', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for needleleaf deciduous boreal tree','gC/m2')

         ! 3: leaf area index dnf boreal
         CALL write_history_variable_2d ( DEF_hist_vars%lai_dnfboreal, &
             a_lai_dnfboreal, file_hist, 'f_lai_dnfboreal', itime_in_file, sumarea, filter, &
             'leaf area index for needleleaf deciduous boreal tree','m2/m2')

         ! 3: npp dnf boreal
         CALL write_history_variable_2d ( DEF_hist_vars%npp_dnfboreal, &
             a_npp_dnfboreal, file_hist, 'f_npp_dnfboreal', itime_in_file, sumarea, filter, &
             'npp for needleleaf deciduous boreal tree','m2/m2')

         ! 3: npp to leafc dnf boreal
         CALL write_history_variable_2d ( DEF_hist_vars%npptoleafc_dnfboreal, &
             a_npptoleafc_dnfboreal, file_hist, 'f_npptoleafc_dnfboreal', itime_in_file, sumarea, filter, &
             'npp to leafc for needleleaf deciduous boreal tree','m2/m2')

         ! 4: gpp ebf trop
         CALL write_history_variable_2d ( DEF_hist_vars%gpp_ebftrop, &
             a_gpp_ebftrop, file_hist, 'f_gpp_ebftrop', itime_in_file, sumarea, filter, &
             'gross primary productivity for broadleaf evergreen tropical tree','gC/m2/s')

         ! 4: leaf carbon display pool ebf trop
         CALL write_history_variable_2d ( DEF_hist_vars%leafc_ebftrop, &
             a_leafc_ebftrop, file_hist, 'f_leafc_ebftrop', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for broadleaf evergreen tropical tree','gC/m2')

         ! 4: leaf area index ebf trop
         CALL write_history_variable_2d ( DEF_hist_vars%lai_ebftrop, &
             a_lai_ebftrop, file_hist, 'f_lai_ebftrop', itime_in_file, sumarea, filter, &
             'leaf area index for broadleaf evergreen tropical tree','m2/m2')

         ! 4: npp ebf trop
         CALL write_history_variable_2d ( DEF_hist_vars%npp_ebftrop, &
             a_npp_ebftrop, file_hist, 'f_npp_ebftrop', itime_in_file, sumarea, filter, &
             'npp for broadleaf evergreen tropical tree','m2/m2')

         ! 4: npp to leafc ebf trop
         CALL write_history_variable_2d ( DEF_hist_vars%npptoleafc_ebftrop, &
             a_npptoleafc_ebftrop, file_hist, 'f_npptoleafc_ebftrop', itime_in_file, sumarea, filter, &
             'npp to leafc for broadleaf evergreen tropical tree','m2/m2')

         ! 5: gpp ebf temp
         CALL write_history_variable_2d ( DEF_hist_vars%gpp_ebftemp, &
             a_gpp_ebftemp, file_hist, 'f_gpp_ebftemp', itime_in_file, sumarea, filter, &
             'gross primary productivity for broadleaf evergreen temperate tree','gC/m2/s')

         ! 5: leaf carbon display pool ebf temp
         CALL write_history_variable_2d ( DEF_hist_vars%leafc_ebftemp, &
             a_leafc_ebftemp, file_hist, 'f_leafc_ebftemp', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for broadleaf evergreen temperate tree','gC/m2')

         ! 5: leaf area index ebf temp
         CALL write_history_variable_2d ( DEF_hist_vars%lai_ebftemp, &
             a_lai_ebftemp, file_hist, 'f_lai_ebftemp', itime_in_file, sumarea, filter, &
             'leaf area index for broadleaf evergreen temperate tree','m2/m2')

         ! 5: npp ebf temp
         CALL write_history_variable_2d ( DEF_hist_vars%npp_ebftemp, &
             a_npp_ebftemp, file_hist, 'f_npp_ebftemp', itime_in_file, sumarea, filter, &
             'npp for broadleaf evergreen temperate tree','m2/m2')

         ! 5: npp to leafc ebf temp
         CALL write_history_variable_2d ( DEF_hist_vars%npptoleafc_ebftemp, &
             a_npptoleafc_ebftemp, file_hist, 'f_npptoleafc_ebftemp', itime_in_file, sumarea, filter, &
             'npp to leafc for broadleaf evergreen temperate tree','m2/m2')

         ! 6: gpp dbf trop
         CALL write_history_variable_2d ( DEF_hist_vars%gpp_dbftrop, &
             a_gpp_dbftrop, file_hist, 'f_gpp_dbftrop', itime_in_file, sumarea, filter, &
             'gross primary productivity for broadleaf deciduous tropical tree','gC/m2/s')

         ! 6: leaf carbon display pool dbf trop
         CALL write_history_variable_2d ( DEF_hist_vars%leafc_dbftrop, &
             a_leafc_dbftrop, file_hist, 'f_leafc_dbftrop', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for broadleaf deciduous tropical tree','gC/m2')

         ! 6: leaf area index dbf trop
         CALL write_history_variable_2d ( DEF_hist_vars%lai_dbftrop, &
             a_lai_dbftrop, file_hist, 'f_lai_dbftrop', itime_in_file, sumarea, filter, &
             'leaf area index for broadleaf evergreen temperate tree','m2/m2')

         ! 6: npp dbf trop
         CALL write_history_variable_2d ( DEF_hist_vars%npp_dbftrop, &
             a_npp_dbftrop, file_hist, 'f_npp_dbftrop', itime_in_file, sumarea, filter, &
             'npp for broadleaf evergreen temperate tree','m2/m2')

         ! 6: npp to leafc dbf trop
         CALL write_history_variable_2d ( DEF_hist_vars%npptoleafc_dbftrop, &
             a_npptoleafc_dbftrop, file_hist, 'f_npptoleafc_dbftrop', itime_in_file, sumarea, filter, &
             'npp to leafc for broadleaf evergreen temperate tree','m2/m2')

         ! 7: gpp dbf temp
         CALL write_history_variable_2d ( DEF_hist_vars%gpp_dbftemp, &
             a_gpp_dbftemp, file_hist, 'f_gpp_dbftemp', itime_in_file, sumarea, filter, &
             'gross primary productivity for broadleaf deciduous temperate tree','gC/m2/s')

         ! 7: leaf carbon display pool dbf temp
         CALL write_history_variable_2d ( DEF_hist_vars%leafc_dbftemp, &
             a_leafc_dbftemp, file_hist, 'f_leafc_dbftemp', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for broadleaf deciduous temperate tree','gC/m2')

         ! 7: leaf area index dbf temp
         CALL write_history_variable_2d ( DEF_hist_vars%lai_dbftemp, &
             a_lai_dbftemp, file_hist, 'f_lai_dbftemp', itime_in_file, sumarea, filter, &
             'leaf area index for broadleaf decidous temperate tree','m2/m2')

         ! 7: npp dbf temp
         CALL write_history_variable_2d ( DEF_hist_vars%npp_dbftemp, &
             a_npp_dbftemp, file_hist, 'f_npp_dbftemp', itime_in_file, sumarea, filter, &
             'npp for broadleaf decidous temperate tree','m2/m2')

         ! 7: npp to leafc dbf temp
         CALL write_history_variable_2d ( DEF_hist_vars%npptoleafc_dbftemp, &
             a_npptoleafc_dbftemp, file_hist, 'f_npptoleafc_dbftemp', itime_in_file, sumarea, filter, &
             'npp to leafc for broadleaf decidous temperate tree','m2/m2')

         ! 8: gpp dbf boreal
         CALL write_history_variable_2d ( DEF_hist_vars%gpp_dbfboreal, &
             a_gpp_dbfboreal, file_hist, 'f_gpp_dbfboreal', itime_in_file, sumarea, filter, &
             'gross primary productivity for broadleaf deciduous boreal tree','gC/m2/s')

         ! 8: leaf carbon display pool dbf boreal
         CALL write_history_variable_2d ( DEF_hist_vars%leafc_dbfboreal, &
             a_leafc_dbfboreal, file_hist, 'f_leafc_dbfboreal', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for broadleaf deciduous boreal tree','gC/m2')

         ! 8: leaf area index dbf boreal
         CALL write_history_variable_2d ( DEF_hist_vars%lai_dbfboreal, &
             a_lai_dbfboreal, file_hist, 'f_lai_dbfboreal', itime_in_file, sumarea, filter, &
             'leaf area index for broadleaf decidous boreal tree','m2/m2')

         ! 8: npp dbf boreal
         CALL write_history_variable_2d ( DEF_hist_vars%npp_dbfboreal, &
             a_npp_dbfboreal, file_hist, 'f_npp_dbfboreal', itime_in_file, sumarea, filter, &
             'npp for broadleaf decidous boreal tree','m2/m2')

         ! 8: npp to leafc dbf boreal
         CALL write_history_variable_2d ( DEF_hist_vars%npptoleafc_dbfboreal, &
             a_npptoleafc_dbfboreal, file_hist, 'f_npptoleafc_dbfboreal', itime_in_file, sumarea, filter, &
             'npp to leafc for broadleaf decidous boreal tree','m2/m2')

         ! 9: gpp ebs temp
         CALL write_history_variable_2d ( DEF_hist_vars%gpp_ebstemp, &
             a_gpp_ebstemp, file_hist, 'f_gpp_ebstemp', itime_in_file, sumarea, filter, &
             'gross primary productivity for broadleaf evergreen temperate shrub','gC/m2/s')

         ! 9: leaf carbon display pool ebs temp
         CALL write_history_variable_2d ( DEF_hist_vars%leafc_ebstemp, &
             a_leafc_ebstemp, file_hist, 'f_leafc_ebstemp', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for broadleaf evergreen temperate shrub','gC/m2')

         ! 9: leaf area index ebs temp
         CALL write_history_variable_2d ( DEF_hist_vars%lai_ebstemp, &
             a_lai_ebstemp, file_hist, 'f_lai_ebstemp', itime_in_file, sumarea, filter, &
             'leaf area index for broadleaf evergreen temperate shrub','m2/m2')

         ! 9: npp ebs temp
         CALL write_history_variable_2d ( DEF_hist_vars%npp_ebstemp, &
             a_npp_ebstemp, file_hist, 'f_npp_ebstemp', itime_in_file, sumarea, filter, &
             'npp for broadleaf evergreen temperate shrub','m2/m2')

         ! 9: npp to leafc ebs temp
         CALL write_history_variable_2d ( DEF_hist_vars%npptoleafc_ebstemp, &
             a_npptoleafc_ebstemp, file_hist, 'f_npptoleafc_ebstemp', itime_in_file, sumarea, filter, &
             'npp to leafc for broadleaf evergreen temperate shrub','m2/m2')

         ! 10: gpp dbs temp
         CALL write_history_variable_2d ( DEF_hist_vars%gpp_dbstemp, &
             a_gpp_dbstemp, file_hist, 'f_gpp_dbstemp', itime_in_file, sumarea, filter, &
             'gross primary productivity for broadleaf deciduous temperate shrub','gC/m2/s')

         ! 10: leaf carbon display pool dbs temp
         CALL write_history_variable_2d ( DEF_hist_vars%leafc_dbstemp, &
             a_leafc_dbstemp, file_hist, 'f_leafc_dbstemp', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for broadleaf deciduous temperate shrub','gC/m2')

         ! 10: leaf area index dbs temp
         CALL write_history_variable_2d ( DEF_hist_vars%lai_dbstemp, &
             a_lai_dbstemp, file_hist, 'f_lai_dbstemp', itime_in_file, sumarea, filter, &
             'leaf area index for broadleaf deciduous temperate shrub','m2/m2')

         ! 10: npp dbs temp
         CALL write_history_variable_2d ( DEF_hist_vars%npp_dbstemp, &
             a_npp_dbstemp, file_hist, 'f_npp_dbstemp', itime_in_file, sumarea, filter, &
             'npp for broadleaf deciduous temperate shrub','m2/m2')

         ! 10: npp to leafc dbs temp
         CALL write_history_variable_2d ( DEF_hist_vars%npptoleafc_dbstemp, &
             a_npptoleafc_dbstemp, file_hist, 'f_npptoleafc_dbstemp', itime_in_file, sumarea, filter, &
             'npp to leafc for broadleaf deciduous temperate shrub','m2/m2')

         ! 11: gpp dbs boreal
         CALL write_history_variable_2d ( DEF_hist_vars%gpp_dbsboreal, &
             a_gpp_dbsboreal, file_hist, 'f_gpp_dbsboreal', itime_in_file, sumarea, filter, &
             'gross primary productivity for broadleaf deciduous boreal shrub','gC/m2/s')

         ! 11: leaf carbon display pool dbs boreal
         CALL write_history_variable_2d ( DEF_hist_vars%leafc_dbsboreal, &
             a_leafc_dbsboreal, file_hist, 'f_leafc_dbsboreal', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for broadleaf deciduous boreal shrub','gC/m2')

         ! 11: leaf area index dbs boreal
         CALL write_history_variable_2d ( DEF_hist_vars%lai_dbsboreal, &
             a_lai_dbsboreal, file_hist, 'f_lai_dbsboreal', itime_in_file, sumarea, filter, &
             'leaf area index for broadleaf deciduous boreal shrub','m2/m2')

         ! 11: npp dbs boreal
         CALL write_history_variable_2d ( DEF_hist_vars%npp_dbsboreal, &
             a_npp_dbsboreal, file_hist, 'f_npp_dbsboreal', itime_in_file, sumarea, filter, &
             'npp for broadleaf deciduous boreal shrub','m2/m2')

         ! 11: npp to leafc dbs boreal
         CALL write_history_variable_2d ( DEF_hist_vars%npptoleafc_dbsboreal, &
             a_npptoleafc_dbsboreal, file_hist, 'f_npptoleafc_dbsboreal', itime_in_file, sumarea, filter, &
             'npp to leafc for broadleaf deciduous boreal shrub','m2/m2')

         ! 12: gpp arctic c3 grass
         CALL write_history_variable_2d ( DEF_hist_vars%gpp_c3arcgrass, &
             a_gpp_c3arcgrass, file_hist, 'f_gpp_c3arcgrass', itime_in_file, sumarea, filter, &
             'gross primary productivity for c3 arctic grass','gC/m2/s')

         ! 12: leaf carbon display pool c3 arctic grass
         CALL write_history_variable_2d ( DEF_hist_vars%leafc_c3arcgrass, &
             a_leafc_c3arcgrass, file_hist, 'f_leafc_c3arcgrass', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for c3 arctic grass','gC/m2')

         ! 12: leaf area index c3 arctic grass
         CALL write_history_variable_2d ( DEF_hist_vars%lai_c3arcgrass, &
             a_lai_c3arcgrass, file_hist, 'f_lai_c3arcgrass', itime_in_file, sumarea, filter, &
             'leaf area index for c3 arctic grass','gC/m2')

         ! 12: npp c3 arctic grass
         CALL write_history_variable_2d ( DEF_hist_vars%npp_c3arcgrass, &
             a_npp_c3arcgrass, file_hist, 'f_npp_c3arcgrass', itime_in_file, sumarea, filter, &
             'npp for c3 arctic grass','gC/m2')

         ! 12: npp to leafc c3 arctic grass
         CALL write_history_variable_2d ( DEF_hist_vars%npptoleafc_c3arcgrass, &
             a_npptoleafc_c3arcgrass, file_hist, 'f_npptoleafc_c3arcgrass', itime_in_file, sumarea, filter, &
             'npp to leafc for c3 arctic grass','gC/m2')

         ! 13: gpp c3 grass
         CALL write_history_variable_2d ( DEF_hist_vars%gpp_c3grass, &
             a_gpp_c3grass, file_hist, 'f_gpp_c3grass', itime_in_file, sumarea, filter, &
             'gross primary productivity for c3 grass','gC/m2/s')

         ! 13: leaf carbon display pool arctic c3 grass
         CALL write_history_variable_2d ( DEF_hist_vars%leafc_c3grass, &
             a_leafc_c3grass, file_hist, 'f_leafc_c3grass', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for c3 arctic grass','gC/m2')

         ! 13: leaf area index arctic c3 grass
         CALL write_history_variable_2d ( DEF_hist_vars%lai_c3grass, &
             a_lai_c3grass, file_hist, 'f_lai_c3grass', itime_in_file, sumarea, filter, &
             'leaf area index for c3 arctic grass','gC/m2')

         ! 13: npp arctic c3 grass
         CALL write_history_variable_2d ( DEF_hist_vars%npp_c3grass, &
             a_npp_c3grass, file_hist, 'f_npp_c3grass', itime_in_file, sumarea, filter, &
             'npp for c3 arctic grass','gC/m2')

         ! 13: npp to leafc arctic c3 grass
         CALL write_history_variable_2d ( DEF_hist_vars%npptoleafc_c3grass, &
             a_npptoleafc_c3grass, file_hist, 'f_npptoleafc_c3grass', itime_in_file, sumarea, filter, &
             'npp to leafc for c3 arctic grass','gC/m2')

         ! 14: gpp c4 grass
         CALL write_history_variable_2d ( DEF_hist_vars%gpp_c4grass, &
             a_gpp_c4grass, file_hist, 'f_gpp_c4grass', itime_in_file, sumarea, filter, &
             'gross primary productivity for c4 grass','gC/m2/s')

         ! 14: leaf carbon display pool arctic c4 grass
         CALL write_history_variable_2d ( DEF_hist_vars%leafc_c4grass, &
             a_leafc_c4grass, file_hist, 'f_leafc_c4grass', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for c4 arctic grass','gC/m2')

         ! 14: leaf area index arctic c4 grass
         CALL write_history_variable_2d ( DEF_hist_vars%lai_c4grass, &
             a_lai_c4grass, file_hist, 'f_lai_c4grass', itime_in_file, sumarea, filter, &
             'leaf area index for c4 arctic grass','gC/m2')

         ! 14: npp arctic c4 grass
         CALL write_history_variable_2d ( DEF_hist_vars%npp_c4grass, &
             a_npp_c4grass, file_hist, 'f_npp_c4grass', itime_in_file, sumarea, filter, &
             'npp for c4 arctic grass','gC/m2')

         ! 14: npp to leafc arctic c4 grass
         CALL write_history_variable_2d ( DEF_hist_vars%npptoleafc_c4grass, &
             a_npptoleafc_c4grass, file_hist, 'f_npptoleafc_c4grass', itime_in_file, sumarea, filter, &
             'npp to leafc for c4 arctic grass','gC/m2')

#ifdef CROP
!*****************************************
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     filter(i) = .true.
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_manunitro (:)
            ENDIF
         ENDIF

         CALL write_history_variable_2d ( DEF_hist_vars%manunitro, &
             vecacc, file_hist, 'f_manunitro', itime_in_file, sumarea, filter, &
             'nitrogen in manure','gN/m2/yr')


         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                      IF(pftclass(patch_pft_s(i)) .eq. 19 .or. pftclass(patch_pft_s(i)) .eq. 20)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_hui (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%huiswheat, &
             vecacc, file_hist, 'f_huiswheat', itime_in_file, sumarea, filter, &
             'heat unit index  (rainfed spring wheat)','unitless')

!************************************************************
!         IF (p_is_worker) THEN
!            IF (numpatch > 0) THEN
!               DO i=1,numpatch
!                  IF(patchclass(i) == 12)THEN
!                     IF(pftclass(patch_pft_s(i)) .eq. 17 .or. pftclass(patch_pft_s(i)) .eq. 18 &
!                   .or. pftclass(patch_pft_s(i)) .eq. 75 .or. pftclass(patch_pft_s(i)) .eq. 76)THEN
!                        filter(i) = .true.
!                     ELSE
!                        filter(i) = .false.
!                     ENDIF
!                  ELSE
!                     filter(i) = .false.
!                  ENDIF
!               ENDDO
!            ENDIF
!         ENDIF
!
!         IF (HistForm == 'Gridded') THEN
!            CALL mp2g_hist%get_sumarea (sumarea, filter)
!         ENDIF
!
!         CALL write_history_variable_2d ( DEF_hist_vars%pdcorn, &
!            a_pdcorn, file_hist, 'f_pdcorn', &
!            itime_in_file, sumarea, filter, 'planting date of corn', 'day')
!
!         IF (p_is_worker) THEN
!            IF (numpatch > 0) THEN
!               DO i=1,numpatch
!                  IF(patchclass(i) == 12)THEN
!                     IF(pftclass(patch_pft_s(i)) .eq. 19 .or. pftclass(patch_pft_s(i)) .eq. 20)THEN
!                        filter(i) = .true.
!                     ELSE
!                        filter(i) = .false.
!                     ENDIF
!                  ELSE
!                     filter(i) = .false.
!                  ENDIF
!               ENDDO
!            ENDIF
!         ENDIF
!
!         IF (HistForm == 'Gridded') THEN
!            CALL mp2g_hist%get_sumarea (sumarea, filter)
!         ENDIF
!
!         CALL write_history_variable_2d ( DEF_hist_vars%pdswheat, &
!            a_pdswheat, file_hist, 'f_pdswheat', &
!            itime_in_file, sumarea, filter,'planting date of spring wheat','day')
!
!         IF (p_is_worker) THEN
!            IF (numpatch > 0) THEN
!               DO i=1,numpatch
!                  IF(patchclass(i) == 12)THEN
!                     IF(pftclass(patch_pft_s(i)) .eq. 21 .or. pftclass(patch_pft_s(i)) .eq. 22)THEN
!                        filter(i) = .true.
!                     ELSE
!                        filter(i) = .false.
!                     ENDIF
!                  ELSE
!                     filter(i) = .false.
!                  ENDIF
!               ENDDO
!            ENDIF
!         ENDIF
!
!         IF (HistForm == 'Gridded') THEN
!            CALL mp2g_hist%get_sumarea (sumarea, filter)
!         ENDIF
!
!         CALL write_history_variable_2d ( DEF_hist_vars%pdwwheat, &
!            a_pdwwheat, file_hist, 'f_pdwwheat', &
!            itime_in_file, sumarea, filter,'planting date of winter wheat','day')
!
!         IF (p_is_worker) THEN
!            IF (numpatch > 0) THEN
!               DO i=1,numpatch
!                  IF(patchclass(i) == 12)THEN
!                     IF(pftclass(patch_pft_s(i)) .eq. 23 .or. pftclass(patch_pft_s(i)) .eq. 24 &
!                   .or. pftclass(patch_pft_s(i)) .eq. 77 .or. pftclass(patch_pft_s(i)) .eq. 78)THEN
!                        filter(i) = .true.
!                     ELSE
!                        filter(i) = .false.
!                     ENDIF
!                  ELSE
!                     filter(i) = .false.
!                  ENDIF
!               ENDDO
!            ENDIF
!         ENDIF
!
!         IF (HistForm == 'Gridded') THEN
!            CALL mp2g_hist%get_sumarea (sumarea, filter)
!         ENDIF
!
!         CALL write_history_variable_2d ( DEF_hist_vars%pdsoybean, &
!            a_pdsoybean, file_hist, 'f_pdsoybean', &
!            itime_in_file, sumarea, filter,'planting date of soybean','day')
!
!         IF (p_is_worker) THEN
!            IF (numpatch > 0) THEN
!               DO i=1,numpatch
!                  IF(patchclass(i) == 12)THEN
!                     IF(pftclass(patch_pft_s(i)) .eq. 41 .or. pftclass(patch_pft_s(i)) .eq. 42)THEN
!                        filter(i) = .true.
!                     ELSE
!                        filter(i) = .false.
!                     ENDIF
!                  ELSE
!                     filter(i) = .false.
!                  ENDIF
!               ENDDO
!            ENDIF
!         ENDIF
!
!         IF (HistForm == 'Gridded') THEN
!            CALL mp2g_hist%get_sumarea (sumarea, filter)
!         ENDIF
!
!         CALL write_history_variable_2d ( DEF_hist_vars%pdcotton, &
!            a_pdcotton, file_hist, 'f_pdcotton', &
!            itime_in_file, sumarea, filter,'planting date of cotton','day')
!
!         IF (p_is_worker) THEN
!            IF (numpatch > 0) THEN
!               DO i=1,numpatch
!                  IF(patchclass(i) == 12)THEN
!                     IF(pftclass(patch_pft_s(i)) .eq. 61 .or. pftclass(patch_pft_s(i)) .eq. 62)THEN
!                        filter(i) = .true.
!                     ELSE
!                        filter(i) = .false.
!                     ENDIF
!                  ELSE
!                     filter(i) = .false.
!                  ENDIF
!               ENDDO
!            ENDIF
!         ENDIF
!
!         IF (HistForm == 'Gridded') THEN
!            CALL mp2g_hist%get_sumarea (sumarea, filter)
!         ENDIF
!
!         CALL write_history_variable_2d ( DEF_hist_vars%pdrice1, &
!            a_pdrice1, file_hist, 'f_pdrice1', &
!            itime_in_file, sumarea, filter,'planting date of rice1','day')
!
!         IF (p_is_worker) THEN
!            IF (numpatch > 0) THEN
!               DO i=1,numpatch
!                  IF(patchclass(i) == 12)THEN
!                     IF(pftclass(patch_pft_s(i)) .eq. 67 .or. pftclass(patch_pft_s(i)) .eq. 68)THEN
!                        filter(i) = .true.
!                     ELSE
!                        filter(i) = .false.
!                     ENDIF
!                  ELSE
!                     filter(i) = .false.
!                  ENDIF
!               ENDDO
!            ENDIF
!         ENDIF
!
!         IF (HistForm == 'Gridded') THEN
!            CALL mp2g_hist%get_sumarea (sumarea, filter)
!         ENDIF
!
!         CALL write_history_variable_2d ( DEF_hist_vars%pdrice2, &
!            a_pdrice2, file_hist, 'f_pdrice2', &
!            itime_in_file, sumarea, filter,'planting date of rice2','day')
!
!         IF (p_is_worker) THEN
!            IF (numpatch > 0) THEN
!               DO i=1,numpatch
!                  IF(patchclass(i) == 12)THEN
!                     IF(pftclass(patch_pft_s(i)) .eq. 67 .or. pftclass(patch_pft_s(i)) .eq. 68)THEN
!                        filter(i) = .true.
!                     ELSE
!                        filter(i) = .false.
!                     ENDIF
!                  ELSE
!                     filter(i) = .false.
!                  ENDIF
!               ENDDO
!            ENDIF
!         ENDIF
!
!         IF (HistForm == 'Gridded') THEN
!            CALL mp2g_hist%get_sumarea (sumarea, filter)
!         ENDIF
!
!         CALL write_history_variable_2d ( DEF_hist_vars%pdsugarcane, &
!            a_pdsugarcane, file_hist, 'f_pdsugarcane', &
!            itime_in_file, sumarea, filter,'planting date of sugarcane','day')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 17 .or. pftclass(patch_pft_s(i)) .eq. 18 &
                   .or. pftclass(patch_pft_s(i)) .eq. 75 .or. pftclass(patch_pft_s(i)) .eq. 76)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         CALL write_history_variable_2d ( DEF_hist_vars%fertnitro_corn, &
            a_fertnitro_corn, file_hist, 'f_fertnitro_corn', &
            itime_in_file, sumarea, filter,'nitrogen fertilizer for corn','gN/m2/yr')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 19 .or. pftclass(patch_pft_s(i)) .eq. 20)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         CALL write_history_variable_2d ( DEF_hist_vars%fertnitro_swheat, &
            a_fertnitro_swheat, file_hist, 'f_fertnitro_swheat', &
            itime_in_file, sumarea, filter,'nitrogen fertilizer for spring wheat','gN/m2/yr')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 21 .or. pftclass(patch_pft_s(i)) .eq. 22)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         CALL write_history_variable_2d ( DEF_hist_vars%fertnitro_wwheat, &
            a_fertnitro_wwheat, file_hist, 'f_fertnitro_wwheat', &
            itime_in_file, sumarea, filter,'nitrogen fertilizer for winter wheat','gN/m2/yr')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 23 .or. pftclass(patch_pft_s(i)) .eq. 24 &
                   .or. pftclass(patch_pft_s(i)) .eq. 77 .or. pftclass(patch_pft_s(i)) .eq. 78)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         CALL write_history_variable_2d ( DEF_hist_vars%fertnitro_soybean, &
            a_fertnitro_soybean, file_hist, 'f_fertnitro_soybean', &
            itime_in_file, sumarea, filter,'nitrogen fertilizer for soybean','gN/m2/yr')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 41 .or. pftclass(patch_pft_s(i)) .eq. 42)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         CALL write_history_variable_2d ( DEF_hist_vars%fertnitro_cotton, &
            a_fertnitro_cotton, file_hist, 'f_fertnitro_cotton', &
            itime_in_file, sumarea, filter,'nitrogen fertilizer for cotton','gN/m2/yr')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 61 .or. pftclass(patch_pft_s(i)) .eq. 62)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         CALL write_history_variable_2d ( DEF_hist_vars%fertnitro_rice1, &
            a_fertnitro_rice1, file_hist, 'f_fertnitro_rice1', &
            itime_in_file, sumarea, filter,'nitrogen fertilizer for rice1','gN/m2/yr')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 61 .or. pftclass(patch_pft_s(i)) .eq. 62)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         CALL write_history_variable_2d ( DEF_hist_vars%fertnitro_rice2, &
            a_fertnitro_rice2, file_hist, 'f_fertnitro_rice2', &
            itime_in_file, sumarea, filter,'nitrogen fertilizer for rice2','gN/m2/yr')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 67 .or. pftclass(patch_pft_s(i)) .eq. 68)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         CALL write_history_variable_2d ( DEF_hist_vars%fertnitro_sugarcane, &
            a_fertnitro_sugarcane, file_hist, 'f_fertnitro_sugarcane', &
            itime_in_file, sumarea, filter,'nitrogen fertilizer for sugarcane','gN/m2/yr')

         IF(DEF_USE_IRRIGATION)THEN
            IF (p_is_worker) THEN
               IF (numpatch > 0) THEN
                  DO i=1,numpatch
                     IF(patchclass(i) == 12)THEN
                        IF(pftclass(patch_pft_s(i)) .eq. 17)THEN
                           filter(i) = .true.
                        ELSE
                           filter(i) = .false.
                        ENDIF
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF

            IF (HistForm == 'Gridded') THEN
               CALL mp2g_hist%get_sumarea (sumarea, filter)
            ENDIF

            CALL write_history_variable_2d ( DEF_hist_vars%irrig_method_corn, &
               a_irrig_method_corn, file_hist, 'f_irrig_method_corn', &
               itime_in_file, sumarea, filter,'irrigation method for corn','-')

            IF (p_is_worker) THEN
               IF (numpatch > 0) THEN
                  DO i=1,numpatch
                     IF(patchclass(i) == 12)THEN
                        IF(pftclass(patch_pft_s(i)).eq.19 .or. pftclass(patch_pft_s(i)).eq.20)THEN
                           filter(i) = .true.
                        ELSE
                           filter(i) = .false.
                        ENDIF
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF

            IF (HistForm == 'Gridded') THEN
               CALL mp2g_hist%get_sumarea (sumarea, filter)
            ENDIF

            CALL write_history_variable_2d ( DEF_hist_vars%irrig_method_swheat, &
               a_irrig_method_swheat, file_hist, 'f_irrig_method_swheat', &
               itime_in_file, sumarea, filter,'irrigation method for spring wheat','-')

            IF (p_is_worker) THEN
               IF (numpatch > 0) THEN
                  DO i=1,numpatch
                     IF(patchclass(i) == 12)THEN
                        IF(pftclass(patch_pft_s(i)).eq.21 .or. pftclass(patch_pft_s(i)).eq.22)THEN
                           filter(i) = .true.
                        ELSE
                           filter(i) = .false.
                        ENDIF
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF

            IF (HistForm == 'Gridded') THEN
               CALL mp2g_hist%get_sumarea (sumarea, filter)
            ENDIF

            CALL write_history_variable_2d ( DEF_hist_vars%irrig_method_wwheat, &
               a_irrig_method_wwheat, file_hist, 'f_irrig_method_wwheat', &
               itime_in_file, sumarea, filter,'irrigation method for winter wheat','-')

            IF (p_is_worker) THEN
               IF (numpatch > 0) THEN
                  DO i=1,numpatch
                     IF(patchclass(i) == 12)THEN
                        IF(pftclass(patch_pft_s(i)).eq.23 .or. pftclass(patch_pft_s(i)).eq.24 &
                      .or. pftclass(patch_pft_s(i)).eq.77 .or. pftclass(patch_pft_s(i)).eq.78)THEN
                           filter(i) = .true.
                        ELSE
                           filter(i) = .false.
                        ENDIF
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF

            IF (HistForm == 'Gridded') THEN
               CALL mp2g_hist%get_sumarea (sumarea, filter)
            ENDIF

            CALL write_history_variable_2d ( DEF_hist_vars%irrig_method_soybean, &
               a_irrig_method_soybean, file_hist, 'f_irrig_method_soybean', &
               itime_in_file, sumarea, filter,'irrigation method for soybean','-')

            IF (p_is_worker) THEN
               IF (numpatch > 0) THEN
                  DO i=1,numpatch
                     IF(patchclass(i) == 12)THEN
                        IF(pftclass(patch_pft_s(i)).eq.41 .or. pftclass(patch_pft_s(i)).eq.42)THEN
                           filter(i) = .true.
                        ELSE
                           filter(i) = .false.
                        ENDIF
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF

            IF (HistForm == 'Gridded') THEN
               CALL mp2g_hist%get_sumarea (sumarea, filter)
            ENDIF

            CALL write_history_variable_2d ( DEF_hist_vars%irrig_method_cotton, &
               a_irrig_method_cotton, file_hist, 'f_irrig_method_cotton', &
               itime_in_file, sumarea, filter,'irrigation method for cotton','-')

            IF (p_is_worker) THEN
               IF (numpatch > 0) THEN
                  DO i=1,numpatch
                     IF(patchclass(i) == 12)THEN
                        IF(pftclass(patch_pft_s(i)).eq.61 .or. pftclass(patch_pft_s(i)).eq.62)THEN
                           filter(i) = .true.
                        ELSE
                           filter(i) = .false.
                        ENDIF
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF

            IF (HistForm == 'Gridded') THEN
               CALL mp2g_hist%get_sumarea (sumarea, filter)
            ENDIF

            CALL write_history_variable_2d ( DEF_hist_vars%irrig_method_rice1, &
               a_irrig_method_rice1, file_hist, 'f_irrig_method_rice1', &
               itime_in_file, sumarea, filter,'irrigation method for rice1','-')

            IF (p_is_worker) THEN
               IF (numpatch > 0) THEN
                  DO i=1,numpatch
                     IF(patchclass(i) == 12)THEN
                        IF(pftclass(patch_pft_s(i)).eq.61 .or. pftclass(patch_pft_s(i)).eq.62)THEN
                           filter(i) = .true.
                        ELSE
                           filter(i) = .false.
                        ENDIF
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF

            IF (HistForm == 'Gridded') THEN
               CALL mp2g_hist%get_sumarea (sumarea, filter)
            ENDIF

            CALL write_history_variable_2d ( DEF_hist_vars%irrig_method_rice2, &
               a_irrig_method_rice2, file_hist, 'f_irrig_method_rice2', &
               itime_in_file, sumarea, filter,'irrigation method for rice2','-')

            IF (p_is_worker) THEN
               IF (numpatch > 0) THEN
                  DO i=1,numpatch
                     IF(patchclass(i) == 12)THEN
                        IF(pftclass(patch_pft_s(i)).eq.67 .or. pftclass(patch_pft_s(i)).eq.68)THEN
                           filter(i) = .true.
                        ELSE
                           filter(i) = .false.
                        ENDIF
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF

            IF (HistForm == 'Gridded') THEN
               CALL mp2g_hist%get_sumarea (sumarea, filter)
            ENDIF

            CALL write_history_variable_2d ( DEF_hist_vars%irrig_method_sugarcane, &
               a_irrig_method_sugarcane, file_hist, 'f_irrig_method_sugarcane', &
               itime_in_file, sumarea, filter,'irrigation method for sugarcane','-')

         ENDIF

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 17)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! planting date of rainfed temperate corn
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_plantdate (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%plantdate_rainfed_temp_corn, &
             vecacc, file_hist, 'f_plantdate_rainfed_temp_corn', itime_in_file, sumarea, filter, &
             'Crop planting date (rainfed temperate corn)','day')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 18)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! planting date of irrigated temperate corn
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_plantdate (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%plantdate_irrigated_temp_corn, &
             vecacc, file_hist, 'f_plantdate_irrigated_temp_corn', itime_in_file, sumarea, filter, &
             'Crop planting date (irrigated temperate corn)','day')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 19)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF


         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! planting date of rainfed spring wheat
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_plantdate (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%plantdate_rainfed_spwheat, &
             vecacc, file_hist, 'f_plantdate_rainfed_spwheat', itime_in_file, sumarea, filter, &
             'Crop planting date (rainfed spring wheat)','day')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 20)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF


         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! planting date of irrigated spring wheat
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_plantdate (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%plantdate_irrigated_spwheat, &
             vecacc, file_hist, 'f_plantdate_irrigated_spwheat', itime_in_file, sumarea, filter, &
             'Crop planting date (irrigated spring wheat)','day')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 21)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! planting date of rainfed winter wheat
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_plantdate (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%plantdate_rainfed_wtwheat, &
             vecacc, file_hist, 'f_plantdate_rainfed_wtwheat', itime_in_file, sumarea, filter, &
             'Crop planting date (rainfed winter wheat)','day')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 22)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! planting date of irrigated winter wheat
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_plantdate (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%plantdate_irrigated_wtwheat, &
             vecacc, file_hist, 'f_plantdate_irrigated_wtwheat', itime_in_file, sumarea, filter, &
             'Crop planting date (irrigated winter wheat)','day')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 23)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! planting date of rainfed temperate soybean
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_plantdate (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%plantdate_rainfed_temp_soybean, &
             vecacc, file_hist, 'f_plantdate_rainfed_temp_soybean', &
             itime_in_file, sumarea, filter, &
             'Crop planting date (rainfed temperate soybean)','day')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 24)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! planting date of irrigated temperate soybean
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_plantdate (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%plantdate_irrigated_temp_soybean, &
             vecacc, file_hist, 'f_plantdate_irrigated_temp_soybean', &
             itime_in_file, sumarea, filter, &
             'Crop planting date (irrigated temperate soybean)','day')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 41)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! planting date of rainfed cotton
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_plantdate (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%plantdate_rainfed_cotton, &
             vecacc, file_hist, 'f_plantdate_rainfed_cotton', itime_in_file, sumarea, filter, &
             'Crop planting date (rainfed cotton)','day')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 42)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! planting date of irrigated cotton
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_plantdate (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%plantdate_irrigated_cotton, &
             vecacc, file_hist, 'f_plantdate_irrigated_cotton', itime_in_file, sumarea, filter, &
             'Crop planting date (irrigated cotton)','day')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 61)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! planting date of rainfed rice
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_plantdate (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%plantdate_rainfed_rice, &
             vecacc, file_hist, 'f_plantdate_rainfed_rice', itime_in_file, sumarea, filter, &
             'Crop planting date (rainfed rice)','day')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 62)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! planting date of irrigated rice
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_plantdate (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%plantdate_irrigated_rice, &
             vecacc, file_hist, 'f_plantdate_irrigated_rice', itime_in_file, sumarea, filter, &
             'Crop planting date (irrigated rice)','day')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 67)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! planting date of rainfed sugarcane
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_plantdate (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%plantdate_rainfed_sugarcane, &
             vecacc, file_hist, 'f_plantdate_rainfed_sugarcane', itime_in_file, sumarea, filter, &
             'Crop planting date (rainfed sugarcane)','day')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 68)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! planting date of irrigated sugarcane
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_plantdate (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%plantdate_irrigated_sugarcane, &
             vecacc, file_hist, 'f_plantdate_irrigated_sugarcane', itime_in_file, sumarea, filter, &
             'Crop planting date (irrigated sugarcane)','day')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 75)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! planting date of rainfed trop corn
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_plantdate (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%plantdate_rainfed_trop_corn, &
             vecacc, file_hist, 'f_plantdate_rainfed_trop_corn', itime_in_file, sumarea, filter, &
             'Crop planting date (rainfed trop corn)','day')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 76)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! planting date of irrigated trop corn
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_plantdate (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%plantdate_irrigated_trop_corn, &
             vecacc, file_hist, 'f_plantdate_irrigated_trop_corn', itime_in_file, sumarea, filter, &
             'Crop planting date (irrigated trop corn)','day')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 77)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! planting date of rainfed trop soybean
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_plantdate (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%plantdate_rainfed_trop_soybean, &
             vecacc, file_hist, 'f_plantdate_rainfed_trop_soybean', &
             itime_in_file, sumarea, filter, &
             'Crop planting date (rainfed trop soybean)','day')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 78)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! planting date of irrigated trop soybean
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_plantdate (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%plantdate_irrigated_trop_soybean, &
             vecacc, file_hist, 'f_plantdate_irrigated_trop_soybean', &
             itime_in_file, sumarea, filter, &
             'Crop planting date (irrigated trop soybean)','day')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 15)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF


         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! planting date of unmanaged crop production
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_plantdate (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%plantdate_unmanagedcrop, &
             vecacc, file_hist, 'f_plantdate_unmanagedcrop', itime_in_file, sumarea, filter, &
             'Crop planting date (unmanaged crop production)','day')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 17)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! grain to corn production carbon
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_grainc_to_cropprodc (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%cropprodc_rainfed_temp_corn, &
             vecacc, file_hist, 'f_cropprodc_rainfed_temp_corn', itime_in_file, sumarea, filter, &
             'Crop production (rainfed temperate corn)','gC/m2/s')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 18)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! grain to corn production carbon
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_grainc_to_cropprodc (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%cropprodc_irrigated_temp_corn, &
             vecacc, file_hist, 'f_cropprodc_irrigated_temp_corn', itime_in_file, sumarea, filter, &
             'Crop production (irrigated temperate corn)','gC/m2/s')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 19)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! grain to spring wheat production carbon
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_grainc_to_cropprodc (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%cropprodc_rainfed_spwheat, &
             vecacc, file_hist, 'f_cropprodc_rainfed_spwheat', itime_in_file, sumarea, filter, &
             'Crop production (rainfed spring wheat)','gC/m2/s')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 20)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! grain to spring wheat production carbon
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_grainc_to_cropprodc (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%cropprodc_irrigated_spwheat, &
             vecacc, file_hist, 'f_cropprodc_irrigated_spwheat', itime_in_file, sumarea, filter, &
             'Crop production (irrigated spring wheat)','gC/m2/s')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 21)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! grain to winter wheat production carbon
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_grainc_to_cropprodc (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%cropprodc_rainfed_wtwheat, &
             vecacc, file_hist, 'f_cropprodc_rainfed_wtwheat', itime_in_file, sumarea, filter, &
             'Crop production (rainfed winter wheat)','gC/m2/s')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 22)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! grain to winter wheat production carbon
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_grainc_to_cropprodc (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%cropprodc_irrigated_wtwheat, &
             vecacc, file_hist, 'f_cropprodc_irrigated_wtwheat', itime_in_file, sumarea, filter, &
             'Crop production (irrigated winter wheat)','gC/m2/s')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 23)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! grain to soybean production carbon
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_grainc_to_cropprodc (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%cropprodc_rainfed_temp_soybean, &
             vecacc, file_hist, 'f_cropprodc_rainfed_temp_soybean', &
             itime_in_file, sumarea, filter, &
             'Crop production (rainfed temperate soybean)','gC/m2/s')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 24)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! grain to soybean production carbon
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_grainc_to_cropprodc (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%cropprodc_irrigated_temp_soybean, &
             vecacc, file_hist, 'f_cropprodc_irrigated_temp_soybean', &
             itime_in_file, sumarea, filter, &
             'Crop production (irrigated temperate soybean)','gC/m2/s')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 41)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! grain to cotton production carbon
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_grainc_to_cropprodc (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%cropprodc_rainfed_cotton, &
             vecacc, file_hist, 'f_cropprodc_rainfed_cotton', itime_in_file, sumarea, filter, &
             'Crop production (rainfed cotton)','gC/m2/s')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 42)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! grain to cotton production carbon
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_grainc_to_cropprodc (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%cropprodc_irrigated_cotton, &
             vecacc, file_hist, 'f_cropprodc_irrigated_cotton', itime_in_file, sumarea, filter, &
             'Crop production (irrigated cotton)','gC/m2/s')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 61)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! grain to rice production carbon
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_grainc_to_cropprodc (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%cropprodc_rainfed_rice, &
             vecacc, file_hist, 'f_cropprodc_rainfed_rice', itime_in_file, sumarea, filter, &
             'Crop production (rainfed rice)','gC/m2/s')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 62)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! grain to rice production carbon
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_grainc_to_cropprodc (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%cropprodc_irrigated_rice, &
             vecacc, file_hist, 'f_cropprodc_irrigated_rice', itime_in_file, sumarea, filter, &
             'Crop production (irrigated rice)','gC/m2/s')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 67)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! grain to sugarcane production carbon
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_grainc_to_cropprodc (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%cropprodc_rainfed_sugarcane, &
             vecacc, file_hist, 'f_cropprodc_rainfed_sugarcane', itime_in_file, sumarea, filter, &
             'Crop production (rainfed sugarcane)','gC/m2/s')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 68)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! grain to sugarcane production carbon
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_grainc_to_cropprodc (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%cropprodc_irrigated_sugarcane, &
             vecacc, file_hist, 'f_cropprodc_irrigated_sugarcane', itime_in_file, sumarea, filter, &
             'Crop production (irrigated sugarcane)','gC/m2/s')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 75)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! grain to sugarcane production carbon
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_grainc_to_cropprodc (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%cropprodc_rainfed_trop_corn, &
             vecacc, file_hist, 'f_cropprodc_rainfed_trop_corn', itime_in_file, sumarea, filter, &
             'Crop production (rainfed_trop_corn)','gC/m2/s')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 76)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! grain to sugarcane production carbon
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_grainc_to_cropprodc (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%cropprodc_irrigated_trop_corn, &
             vecacc, file_hist, 'f_cropprodc_irrigated_trop_corn', itime_in_file, sumarea, filter, &
             'Crop production (irrigated_trop_corn)','gC/m2/s')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 77)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! grain to sugarcane production carbon
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_grainc_to_cropprodc (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%cropprodc_rainfed_trop_soybean, &
             vecacc, file_hist, 'f_cropprodc_rainfed_trop_soybean', &
             itime_in_file, sumarea, filter, &
             'Crop production (rainfed trop soybean)','gC/m2/s')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 78)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! grain to sugarcane production carbon
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_grainc_to_cropprodc (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%cropprodc_irrigated_trop_soybean, &
             vecacc, file_hist, 'f_cropprodc_irrigated_trop_soybean', &
             itime_in_file, sumarea, filter, &
             'Crop production (irrigated trop soybean)','gC/m2/s')

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO i=1,numpatch
                  IF(patchclass(i) == 12)THEN
                     IF(pftclass(patch_pft_s(i)) .eq. 15)THEN
                        filter(i) = .true.
                     ELSE
                        filter(i) = .false.
                     ENDIF
                  ELSE
                     filter(i) = .false.
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! grain to unmanaged crop production carbon
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               vecacc (:) = a_grainc_to_cropprodc (:)
            ENDIF
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%cropprodc_unmanagedcrop, &
             vecacc, file_hist, 'f_cropprodc_unmanagedcrop', itime_in_file, sumarea, filter, &
             'Crop production (unmanaged crop production)','gC/m2/s')
#endif
#endif
         ! --------------------------------------------------------------------
         ! Temperature and water (excluding land water bodies and ocean patches)
         ! [soil => 0; urban and built-up => 1; wetland => 2; land ice => 3;
         !  land water bodies => 4; ocean => 99]
         ! --------------------------------------------------------------------

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN

               filter(:) = patchtype <= 3

               IF (DEF_forcing%has_missing_value) THEN
                  filter = filter .and. forcmask_pch
               ENDIF

               filter = filter .and. patchmask
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! soil temperature [K]
         CALL write_history_variable_3d ( DEF_hist_vars%t_soisno, &
            a_t_soisno, file_hist, 'f_t_soisno', &
            itime_in_file, 'soilsnow', maxsnl+1, nl_soil-maxsnl, &
            sumarea, filter, 'soil temperature','K')

         ! liquid water in soil layers [kg/m2]
         CALL write_history_variable_3d ( DEF_hist_vars%wliq_soisno, &
            a_wliq_soisno, file_hist, 'f_wliq_soisno', &
            itime_in_file, 'soilsnow', maxsnl+1, nl_soil-maxsnl, &
            sumarea, filter,'liquid water in soil layers','kg/m2')

         ! ice lens in soil layers [kg/m2]
         CALL write_history_variable_3d ( DEF_hist_vars%wice_soisno, &
            a_wice_soisno, file_hist, 'f_wice_soisno', &
            itime_in_file, 'soilsnow', maxsnl+1, nl_soil-maxsnl, &
            sumarea, filter, 'ice lens in soil layers', 'kg/m2')


#ifdef DataAssimilation
         IF (p_is_worker) THEN
            allocate (a_wliq_h2osoi_5cm     (numpatch                 )); a_wliq_h2osoi_5cm         = spval
            allocate (a_t_soisno_5cm        (numpatch                 )); a_t_soisno_5cm            = spval

            allocate (a_wliq_soisno_ens_mean(maxsnl+1:nl_soil,numpatch)); a_wliq_soisno_ens_mean    = spval
            allocate (a_wliq_soisno_5cm_ens (DEF_DA_ENS_NUM,numpatch  )); a_wliq_soisno_5cm_ens     = spval
            allocate (a_wliq_h2osoi_5cm_a   (numpatch                 )); a_wliq_h2osoi_5cm_a       = spval

            allocate (a_t_soisno_ens_mean   (maxsnl+1:nl_soil,numpatch)); a_t_soisno_ens_mean       = spval
            allocate (a_t_soisno_5cm_ens    (DEF_DA_ENS_NUM,numpatch  )); a_t_soisno_5cm_ens        = spval
            allocate (a_t_soisno_5cm_a      (numpatch                 )); a_t_soisno_5cm_a          = spval

            allocate (a_t_brt_smap_a        (2,numpatch               )); a_t_brt_smap_a            = spval
            allocate (a_t_brt_fy3d_a        (2,numpatch               )); a_t_brt_fy3d_a            = spval

            allocate (a_wliq_soisno_5cm_ens_std(numpatch              )); a_wliq_soisno_5cm_ens_std = spval
            allocate (a_t_soisno_5cm_ens_std   (numpatch              )); a_t_soisno_5cm_ens_std    = spval
            allocate (a_t_brt_smap_ens_std     (2,numpatch            )); a_t_brt_smap_ens_std      = spval
            allocate (a_t_brt_fy3d_ens_std     (2,numpatch            )); a_t_brt_fy3d_ens_std      = spval
         END IF

         IF (p_is_worker) THEN
!#############################################################################
! States before DA
!#############################################################################
            ! calculate surface liquid soil moisture (0-5cm) before DA
            a_wliq_h2osoi_5cm = (a_wliq_soisno(1,:) + a_wliq_soisno(2,:) + &
               a_wliq_soisno(3, :)*(0.05 - 0.0451)/(0.0906 - 0.0451))/(0.05*denh2o)

            ! calculate surface liquid soil moisture (0-5cm) before DA
            a_t_soisno_5cm = (a_t_soisno(1,:)*0.0175 + a_t_soisno(2,:)*(0.0451 - 0.0175))/(0.0451)

!#############################################################################
! States after DA
!#############################################################################
            ! calculate surface liquid soil moisture (0-5cm) after DA
            a_wliq_soisno_ens_mean = sum(a_wliq_soisno_ens, dim=2) / DEF_DA_ENS_NUM
            a_wliq_soisno_5cm_ens = (a_wliq_soisno_ens(1,:,:) + a_wliq_soisno_ens(2,:,:) + &
               a_wliq_soisno_ens(3,:,:)*(0.05-0.0451)/(0.0906-0.0451))/(0.05*denh2o)
            a_wliq_h2osoi_5cm_a = (a_wliq_soisno_ens_mean(1,:) + a_wliq_soisno_ens_mean(2,:) + &
               a_wliq_soisno_ens_mean(3,:)*(0.05 - 0.0451)/(0.0906 - 0.0451))/(0.05*denh2o)

            ! calculate surface soil temperature (0-5cm) before DA & after DA
            a_t_soisno_ens_mean = sum(a_t_soisno_ens, dim=2) / DEF_DA_ENS_NUM
            a_t_soisno_5cm_ens = (a_t_soisno_ens(1,:,:)*0.0175 + a_t_soisno_ens(2,:,:)*(0.0451 - 0.0175))/(0.0451)
            a_t_soisno_5cm_a = (a_t_soisno_ens_mean(1,:)*0.0175 + a_t_soisno_ens_mean(2,:)*(0.0451 - 0.0175))/(0.0451)

!#############################################################################
! brightness temperature after DA
!#############################################################################
            a_t_brt_smap_a = sum(a_t_brt_smap_ens, dim=2) / DEF_DA_ENS_NUM
            a_t_brt_fy3d_a = sum(a_t_brt_fy3d_ens, dim=2) / DEF_DA_ENS_NUM

!#############################################################################
! Standard deviation of states and brightness temperature
!#############################################################################
            ! calculate standard deviation of surface soil moisture, temperature and brightness temperature
            DO np = 1, numpatch
               a_wliq_soisno_5cm_ens_std(np) = &
                  sqrt(sum((a_wliq_soisno_5cm_ens(:,np)-a_wliq_h2osoi_5cm_a(np))**2)/real(DEF_DA_ENS_NUM-1))
               a_t_soisno_5cm_ens_std(np) = &
                  sqrt(sum((a_t_soisno_5cm_ens(:,np)-a_t_soisno_5cm_a(np))**2)/real(DEF_DA_ENS_NUM-1))
               IF (DEF_DA_SM_SMAP) THEN
                  IF (patchtype(np) >= 3) cycle
                  a_t_brt_smap_ens_std(1,np) = &
                     sqrt(sum((a_t_brt_smap_ens(1,:,np)-a_t_brt_smap_a(1,np))**2)/real(DEF_DA_ENS_NUM-1))
                  a_t_brt_smap_ens_std(2,np) = &
                     sqrt(sum((a_t_brt_smap_ens(2,:,np)-a_t_brt_smap_a(2,np))**2)/real(DEF_DA_ENS_NUM-1))
               ENDIF
               IF (DEF_DA_SM_FY) THEN
                  IF (patchtype(np) >= 3) cycle
                  a_t_brt_fy3d_ens_std(1,np) = &
                     sqrt(sum((a_t_brt_fy3d_ens(1,:,np)-a_t_brt_fy3d_a(1,np))**2)/real(DEF_DA_ENS_NUM-1))
                  a_t_brt_fy3d_ens_std(2,np) = &
                     sqrt(sum((a_t_brt_fy3d_ens(2,:,np)-a_t_brt_fy3d_a(2,np))**2)/real(DEF_DA_ENS_NUM-1))
               ENDIF
            ENDDO
         ENDIF

         ! surface soil moisture (0-5cm) before and after DA
         CALL write_history_variable_2d(DEF_hist_vars%DA_wliq_h2osoi_5cm, &
            a_wliq_h2osoi_5cm, file_hist, 'f_wliq_h2osoi_5cm', itime_in_file, &
            sumarea, filter, 'Volumetric liquid water content in 0-5cm', 'm3/m3')
         CALL write_history_variable_2d(DEF_hist_vars%DA_wliq_h2osoi_5cm_a, &
            a_wliq_h2osoi_5cm_a, file_hist, 'f_wliq_h2osoi_5cm_a', itime_in_file, &
            sumarea, filter, 'Analysis volumetric liquid water content in 0-5cm', 'm3/m3')

         ! surface soil temperature (0-5cm) before and after DA
         CALL write_history_variable_2d(DEF_hist_vars%DA_t_soisno_5cm, &
            a_t_soisno_5cm, file_hist, 'f_t_soisno_5cm', itime_in_file, &
            sumarea, filter, 'Soil temperature in 0-5cm', 'K')
         CALL write_history_variable_2d(DEF_hist_vars%DA_t_soisno_5cm_a, &
            a_t_soisno_5cm_a, file_hist, 'f_t_soisno_5cm_a', itime_in_file, &
            sumarea, filter, 'Analysis soil temperature in 0-5cm', 'K')

         ! ensemble soil moisture & temperature in soil layers [kg/m2]
         IF (DEF_DA_ENS_NUM > 1) THEN
            CALL write_history_variable_4d(DEF_hist_vars%DA_wliq_soisno_ens, &
               a_wliq_soisno_ens, file_hist, 'f_wliq_soisno_ens', itime_in_file, 'soilsnow', maxsnl + 1, nl_soil - maxsnl, &
               'ens', 1, DEF_DA_ENS_NUM, sumarea, filter, 'ensemble liquid water in soil layers', 'kg/m2')
            CALL write_history_variable_4d(DEF_hist_vars%DA_t_soisno_ens, &
               a_t_soisno_ens, file_hist, 'f_t_soisno_ens', itime_in_file, 'soilsnow', maxsnl + 1, nl_soil - maxsnl, &
               'ens', 1, DEF_DA_ENS_NUM, sumarea, filter, 'ensemble soil temperature', 'K')
         ENDIF

         ! standard deviation of ensemble surface soil moisture and temperature (0-5cm)
         CALL write_history_variable_2d(DEF_hist_vars%DA_wliq_soisno_5cm_ens_std, &
            a_wliq_soisno_5cm_ens_std, file_hist, 'f_wliq_soisno_ens_5cm_ens_std', itime_in_file, &
            sumarea, filter, 'Standard deviation of ensemble volumetric liquid water content in 0-5cm', 'm3/m3')
         CALL write_history_variable_2d(DEF_hist_vars%DA_t_soisno_5cm_ens_std, &
            a_t_soisno_5cm_ens_std, file_hist, 'f_t_soisno_ens_5cm_ens_std', itime_in_file, &
            sumarea, filter, 'Standard deviation of ensemble soil temperature in 0-5cm', 'K')

         ! --------------------------------------------------------------------
         ! brightness temperature (excluding land ice, land water bodies and ocean patches)
         ! [soil => 0; urban and built-up => 1; wetland => 2; land ice => 3;
         !  land water bodies => 4; ocean => 99]
         ! --------------------------------------------------------------------

         IF (HistForm == 'Gridded') THEN
            IF (p_is_io) THEN
               CALL allocate_block_data(ghist, sumarea)
            END IF
         END IF

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN

               filter(:) = patchtype <= 2

               IF (DEF_forcing%has_missing_value) THEN
                  filter = filter .and. forcmask_pch
               END IF

               filter = filter .and. patchmask
            END IF
         END IF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea(sumarea, filter)
         END IF

         ! brightness temperature for SMAP and FY satellites
         IF (DEF_DA_SM_SMAP) THEN
            CALL write_history_variable_3d(DEF_hist_vars%DA_t_brt_smap, &
               a_t_brt_smap, file_hist, 'f_t_brt_smap', itime_in_file, 'band', 1, 2, sumarea, filter, &
               'H- & V- polarized brightness temperature for SMAP satellite (L-band, 1.4GHz)', 'K')
            CALL write_history_variable_3d(DEF_hist_vars%DA_t_brt_smap_a, &
               a_t_brt_smap_a, file_hist, 'f_t_brt_smap_a', itime_in_file, 'band', 1, 2, sumarea, filter, &
               'Analysis H- & V- polarized brightness temperature for SMAP satellite (L-band,1.4GHz)', 'K')
            IF (DEF_DA_ENS_NUM > 1) THEN
               CALL write_history_variable_4d(DEF_hist_vars%DA_t_brt_smap_ens, &
                  a_t_brt_smap_ens, file_hist, 'f_t_brt_smap_ens', itime_in_file, 'band', 1, 2, 'ens', 1, DEF_DA_ENS_NUM, &
                  sumarea, filter, 'ensemble H- & V- polarized brightness temperature for SMAP satellite (L-band,1.4GHz)', 'K')
            END IF
            CALL write_history_variable_3d(DEF_hist_vars%DA_t_brt_smap_ens_std, &
               a_t_brt_smap_ens_std, file_hist, 'f_t_brt_smap_ens_std', itime_in_file, 'band', 1, 2, &
               sumarea, filter, 'Standard deviation of H- & V- polarized brightness temperature for SMAP satellite (L-band,1.4GHz)', 'K')
         ENDIF

         IF (DEF_DA_SM_FY) THEN
            CALL write_history_variable_3d(DEF_hist_vars%DA_t_brt_fy3d, &
               a_t_brt_fy3d, file_hist, 'f_t_brt_fy3d', itime_in_file, 'band', 1, 2, sumarea, filter, &
               'H- & V- polarized brightness temperature for FY satellite (X-band, 10.65GHz)', 'K')
            CALL write_history_variable_3d(DEF_hist_vars%DA_t_brt_fy3d_a, &
               a_t_brt_fy3d_a, file_hist, 'f_t_brt_fy3d_a', itime_in_file, 'band', 1, 2, sumarea, filter, &
               'Analysis H- & V- polarized brightness temperature for FY satellite (X-band, 10.65GHz)', 'K')
            IF (DEF_DA_ENS_NUM > 1) THEN
               CALL write_history_variable_4d(DEF_hist_vars%DA_t_brt_fy3d_ens, &
                  a_t_brt_fy3d_ens, file_hist, 'f_t_brt_fy3d_ens', itime_in_file, 'band', 1, 2, 'ens', 1, DEF_DA_ENS_NUM, &
                  sumarea, filter, 'ensemble H- & V- polarized brightness temperature for FY satellite (X-band, 10.65GHz)', 'K')
            END IF
            CALL write_history_variable_3d(DEF_hist_vars%DA_t_brt_fy3d_ens_std, &
               a_t_brt_fy3d_ens_std, file_hist, 'f_t_brt_fy3d_ens_std', itime_in_file, 'band', 1, 2, &
               sumarea, filter, 'Standard deviation of H- & V- polarized brightness temperature for FY satellite (X-band, 10.65GHz)', 'K')
         ENDIF

         IF (p_is_worker) THEN
            deallocate (a_wliq_h2osoi_5cm)
            deallocate (a_t_soisno_5cm)
            deallocate (a_wliq_soisno_ens_mean)
            deallocate (a_wliq_soisno_5cm_ens)
            deallocate (a_wliq_h2osoi_5cm_a)
            deallocate (a_t_soisno_ens_mean)
            deallocate (a_t_soisno_5cm_ens)
            deallocate (a_t_soisno_5cm_a)
            deallocate (a_t_brt_smap_a)
            deallocate (a_t_brt_fy3d_a)
            deallocate (a_wliq_soisno_5cm_ens_std)
            deallocate (a_t_soisno_5cm_ens_std)
            deallocate (a_t_brt_smap_ens_std)
            deallocate (a_t_brt_fy3d_ens_std)
         END IF

#endif

         ! --------------------------------------------------------------------
         ! additional diagnostic variables for output (vegetated land only <=2)
         ! [soil => 0; urban and built-up => 1; wetland => 2; land ice => 3;
         !  land water bodies => 4; ocean => 99]
         ! --------------------------------------------------------------------

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN

               filter(:) = patchtype <= 2

               IF (DEF_forcing%has_missing_value) THEN
                  filter = filter .and. forcmask_pch
               ENDIF

               filter = filter .and. patchmask
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! volumetric soil water in layers [m3/m3]
         CALL write_history_variable_3d ( DEF_hist_vars%h2osoi, &
            a_h2osoi, file_hist, 'f_h2osoi', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
            'volumetric water in soil layers','m3/m3')

         IF (DEF_USE_VariablySaturatedFlow) THEN
            ! water flux between water layers [mm h2o/s]
            CALL write_history_variable_3d ( DEF_hist_vars%qlayer, &
               a_qlayer, file_hist, 'f_qlayer', itime_in_file, 'soilinterface', 0, nl_soil+1, &
               sumarea, filter, 'water flux between soil layers','mm/s')
         ENDIF

         ! fraction of root water uptake from each soil layer, all layers add to 1,
         ! when PHS is not defined water exchange between soil layers and root.
         ! Positive: soil->root [mm h2o/s], when PHS is defined
         CALL write_history_variable_3d ( DEF_hist_vars%rootr, &
            a_rootr, file_hist, 'f_rootr', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
            'root water uptake', 'mm h2o/s')

         IF (DEF_USE_PLANTHYDRAULICS) THEN
         ! vegetation water potential [mm]
            CALL write_history_variable_3d ( DEF_hist_vars%vegwp, &
               a_vegwp, file_hist, 'f_vegwp', &
               itime_in_file, 'vegnodes', 1, nvegwcs, sumarea, filter, &
               'vegetation water potential', 'mm')
         ENDIF

         ! water table depth [m]
         CALL write_history_variable_2d ( DEF_hist_vars%zwt, &
            a_zwt, file_hist, 'f_zwt', itime_in_file, sumarea, filter, &
            'the depth to water table','m')

         ! --------------------------------------------------------------------
         ! depth of surface water (including land ice and ocean patches)
         ! --------------------------------------------------------------------
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN

               filter(:) = (patchtype <= 4)

               IF (DEF_forcing%has_missing_value) THEN
                  filter = filter .and. forcmask_pch
               ENDIF

               filter = filter .and. patchmask
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! water storage in aquifer [mm]
         CALL write_history_variable_2d ( DEF_hist_vars%wa, &
            a_wa, file_hist, 'f_wa', itime_in_file, sumarea, filter, &
            'water storage in aquifer','mm')

         ! instantaneous water storage in aquifer [mm]
         IF (p_is_worker) THEN
            vecacc = wa
            WHERE(vecacc /= spval) vecacc = vecacc * nac
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%wa_inst, &
            vecacc, file_hist, 'f_wa_inst', itime_in_file, sumarea, filter, &
            'instantaneous water storage in aquifer','mm')

         ! depth of surface water [mm]
         CALL write_history_variable_2d ( DEF_hist_vars%wdsrf, &
            a_wdsrf, file_hist, 'f_wdsrf', itime_in_file, sumarea, filter, &
            'depth of surface water','mm')

         ! instantaneous depth of surface water [mm]
         IF (p_is_worker) THEN
            vecacc = wdsrf
            WHERE(vecacc /= spval) vecacc = vecacc * nac
         ENDIF
         CALL write_history_variable_2d ( DEF_hist_vars%wdsrf_inst, &
            vecacc, file_hist, 'f_wdsrf_inst', itime_in_file, sumarea, filter, &
            'instantaneous depth of surface water','mm')

         ! -----------------------------------------------
         ! Land water bodies' ice fraction and temperature
         ! -----------------------------------------------

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               filter(:) = patchtype == 4
               IF (DEF_forcing%has_missing_value) THEN
                  filter = filter .and. forcmask_pch
               ENDIF
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         IF (HistForm == 'Gridded') THEN
            IF (trim(file_hist) /= trim(file_last)) THEN
               CALL hist_write_var_real8_2d (file_hist, 'area_lake', ghist, -1, sumarea, &
                  compress = 1, longname = 'area of lake', units = 'km2')
            ENDIF
         ENDIF

         ! lake layer depth [m]
         CALL write_history_variable_3d ( DEF_hist_vars%dz_lake .and. DEF_USE_Dynamic_Lake, &
            a_dz_lake, file_hist, 'f_dz_lake', itime_in_file, 'lake', 1, nl_lake, sumarea, filter, &
            'lake layer thickness','m')

         ! lake temperature [K]
         CALL write_history_variable_3d ( DEF_hist_vars%t_lake, &
            a_t_lake, file_hist, 'f_t_lake', itime_in_file, 'lake', 1, nl_lake, sumarea, filter, &
            'lake temperature','K')

         ! lake ice fraction cover [0-1]
         CALL write_history_variable_3d ( DEF_hist_vars%lake_icefrac, &
            a_lake_icefrac, file_hist, 'f_lake_icefrac', itime_in_file, 'lake', 1, nl_lake, &
            sumarea, filter, 'lake ice fraction cover','0-1')

         ! lake water deficit due to evaporation [mm/s]
         IF (.not. DEF_USE_Dynamic_Lake) THEN
            CALL write_history_variable_2d ( DEF_hist_vars%lake_deficit, &
               a_lake_deficit, file_hist, 'f_lake_deficit', itime_in_file, sumarea, filter, &
               'lake water deficit due to evaporation','mm/s')
         ENDIF

#ifdef EXTERNAL_LAKE
         CALL LakeVarsSaveHist (nl_lake, file_hist, HistForm, itime_in_file, sumarea, filter)
#endif

         ! --------------------------------
         ! Retrieve through averaged fluxes
         ! --------------------------------
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN

               filter(:) = patchtype < 99

               IF (DEF_forcing%has_missing_value) THEN
                  filter = filter .and. forcmask_pch
               ENDIF

               filter = filter .and. patchmask
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! u* in similarity theory [m/s]
         CALL write_history_variable_2d ( DEF_hist_vars%ustar, &
            a_ustar, file_hist, 'f_ustar', itime_in_file, sumarea, filter, &
            'u* in similarity theory based on patch','m/s')

         ! u* in similarity theory [m/s]
         CALL write_history_variable_2d ( DEF_hist_vars%ustar2, &
            a_ustar2, file_hist, 'f_ustar2', itime_in_file, sumarea, filter, &
            'u* in similarity theory based on grid','m/s')

         ! t* in similarity theory [K]
         CALL write_history_variable_2d ( DEF_hist_vars%tstar, &
            a_tstar, file_hist, 'f_tstar', itime_in_file, sumarea, filter, &
            't* in similarity theory','K')

         ! q* in similarity theory [kg/kg]
         CALL write_history_variable_2d ( DEF_hist_vars%qstar, &
            a_qstar, file_hist, 'f_qstar', itime_in_file, sumarea, filter, &
            'q* in similarity theory', 'kg/kg')

         ! dimensionless height (z/L) used in Monin-Obukhov theory
         CALL write_history_variable_2d ( DEF_hist_vars%zol, &
            a_zol, file_hist, 'f_zol', itime_in_file, sumarea, filter, &
            'dimensionless height (z/L) used in Monin-Obukhov theory','-')

         ! bulk Richardson number in surface layer
         CALL write_history_variable_2d ( DEF_hist_vars%rib, &
            a_rib, file_hist, 'f_rib', itime_in_file, sumarea, filter, &
            'bulk Richardson number in surface layer','-')

         ! integral of profile FUNCTION for momentum
         CALL write_history_variable_2d ( DEF_hist_vars%fm, &
            a_fm, file_hist, 'f_fm', itime_in_file, sumarea, filter, &
            'integral of profile FUNCTION for momentum','-')

         ! integral of profile FUNCTION for heat
         CALL write_history_variable_2d ( DEF_hist_vars%fh, &
            a_fh, file_hist, 'f_fh', itime_in_file, sumarea, filter, &
            'integral of profile FUNCTION for heat','-')

         ! integral of profile FUNCTION for moisture
         CALL write_history_variable_2d ( DEF_hist_vars%fq, &
            a_fq, file_hist, 'f_fq', itime_in_file, sumarea, filter, &
            'integral of profile FUNCTION for moisture','-')

         ! 10m u-velocity [m/s]
         CALL write_history_variable_2d ( DEF_hist_vars%us10m, &
            a_us10m, file_hist, 'f_us10m', itime_in_file, sumarea, filter, &
            '10m u-velocity','m/s')

         ! 10m v-velocity [m/s]
         CALL write_history_variable_2d ( DEF_hist_vars%vs10m, &
            a_vs10m, file_hist, 'f_vs10m', itime_in_file, sumarea, filter, &
            '10m v-velocity','m/s')

         ! integral of profile FUNCTION for momentum at 10m [-]
         CALL write_history_variable_2d ( DEF_hist_vars%fm10m, &
            a_fm10m, file_hist, 'f_fm10m', itime_in_file, sumarea, filter, &
            'integral of profile FUNCTION for momentum at 10m','-')

         ! total reflected solar radiation (W/m2)
         CALL write_history_variable_2d ( DEF_hist_vars%sr, &
            a_sr, file_hist, 'f_sr', itime_in_file, sumarea, filter, &
            'reflected solar radiation at surface [W/m2]','W/m2')

         ! incident direct beam vis solar radiation (W/m2)
         CALL write_history_variable_2d ( DEF_hist_vars%solvd, &
            a_solvd, file_hist, 'f_solvd', itime_in_file, sumarea, filter, &
            'incident direct beam vis solar radiation (W/m2)','W/m2')

         ! incident diffuse beam vis solar radiation (W/m2)
         CALL write_history_variable_2d ( DEF_hist_vars%solvi, &
            a_solvi, file_hist, 'f_solvi', itime_in_file, sumarea, filter, &
            'incident diffuse beam vis solar radiation (W/m2)','W/m2')

         ! incident direct beam nir solar radiation (W/m2)
         CALL write_history_variable_2d ( DEF_hist_vars%solnd, &
            a_solnd, file_hist, 'f_solnd', itime_in_file, sumarea, filter, &
            'incident direct beam nir solar radiation (W/m2)','W/m2')

         ! incident diffuse beam nir solar radiation (W/m2)
         CALL write_history_variable_2d ( DEF_hist_vars%solni, &
            a_solni, file_hist, 'f_solni', itime_in_file, sumarea, filter, &
            'incident diffuse beam nir solar radiation (W/m2)','W/m2')

         ! reflected direct beam vis solar radiation (W/m2)
         CALL write_history_variable_2d ( DEF_hist_vars%srvd, &
            a_srvd, file_hist, 'f_srvd', itime_in_file, sumarea, filter, &
            'reflected direct beam vis solar radiation (W/m2)','W/m2')

         ! reflected diffuse beam vis solar radiation (W/m2)
         CALL write_history_variable_2d ( DEF_hist_vars%srvi, &
            a_srvi, file_hist, 'f_srvi', itime_in_file, sumarea, filter, &
            'reflected diffuse beam vis solar radiation (W/m2)','W/m2')

         ! reflected direct beam nir solar radiation (W/m2)
         CALL write_history_variable_2d ( DEF_hist_vars%srnd, &
            a_srnd, file_hist, 'f_srnd', itime_in_file, sumarea, filter, &
            'reflected direct beam nir solar radiation (W/m2)','W/m2')

         ! reflected diffuse beam nir solar radiation (W/m2)
         CALL write_history_variable_2d ( DEF_hist_vars%srni, &
            a_srni, file_hist, 'f_srni', itime_in_file, sumarea, filter, &
            'reflected diffuse beam nir solar radiation (W/m2)','W/m2')

         ! local noon fluxes
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               filter(:) = nac_ln > 0
               IF (DEF_forcing%has_missing_value) THEN
                  filter = filter .and. forcmask_pch
               ENDIF
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! incident direct beam vis solar radiation at local noon (W/m2)
         CALL write_history_variable_2d ( DEF_hist_vars%solvdln, &
            a_solvdln, file_hist, 'f_solvdln', itime_in_file, sumarea, filter, &
            'incident direct beam vis solar radiation at local noon(W/m2)','W/m2',nac_ln)

         ! incident diffuse beam vis solar radiation at local noon (W/m2)
         CALL write_history_variable_2d ( DEF_hist_vars%solviln, &
            a_solviln, file_hist, 'f_solviln', itime_in_file, sumarea, filter, &
            'incident diffuse beam vis solar radiation at local noon(W/m2)','W/m2',nac_ln)

         ! incident direct beam nir solar radiation at local noon (W/m2)
         CALL write_history_variable_2d ( DEF_hist_vars%solndln, &
            a_solndln, file_hist, 'f_solndln', itime_in_file, sumarea, filter, &
            'incident direct beam nir solar radiation at local noon(W/m2)','W/m2',nac_ln)

         ! incident diffuse beam nir solar radiation at local noon (W/m2)
         CALL write_history_variable_2d ( DEF_hist_vars%solniln, &
            a_solniln, file_hist, 'f_solniln', itime_in_file, sumarea, filter, &
            'incident diffuse beam nir solar radiation at local noon(W/m2)','W/m2',nac_ln)

         ! reflected direct beam vis solar radiation at local noon (W/m2)
         CALL write_history_variable_2d ( DEF_hist_vars%srvdln, &
            a_srvdln, file_hist, 'f_srvdln', itime_in_file, sumarea, filter, &
            'reflected direct beam vis solar radiation at local noon(W/m2)','W/m2',nac_ln)

         ! reflected diffuse beam vis solar radiation at local noon (W/m2)
         CALL write_history_variable_2d ( DEF_hist_vars%srviln, &
            a_srviln, file_hist, 'f_srviln', itime_in_file, sumarea, filter, &
            'reflected diffuse beam vis solar radiation at local noon(W/m2)','W/m2',nac_ln)

         ! reflected direct beam nir solar radiation at local noon (W/m2)
         CALL write_history_variable_2d ( DEF_hist_vars%srndln, &
            a_srndln, file_hist, 'f_srndln', itime_in_file, sumarea, filter, &
            'reflected direct beam nir solar radiation at local noon(W/m2)','W/m2',nac_ln)

         ! reflected diffuse beam nir solar radiation at local noon(W/m2)
         CALL write_history_variable_2d ( DEF_hist_vars%srniln, &
            a_srniln, file_hist, 'f_srniln', itime_in_file, sumarea, filter, &
            'reflected diffuse beam nir solar radiation at local noon(W/m2)','W/m2',nac_ln)


         IF ((p_is_worker) .and. (numpatch > 0)) THEN
            filter = (patchtype == 0) .and. patchmask
            IF (DEF_forcing%has_missing_value) filter = filter .and. forcmask_pch
         ENDIF
         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         CALL write_history_variable_3d ( DEF_hist_vars%sensors, &
            a_sensors, file_hist, 'f_sensors', itime_in_file, 'sensor', 1, nsensor, &
            sumarea, filter, 'variable sensors','user defined')

#if (defined CaMa_Flood)
#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif
         IF (p_is_master) THEN
            CALL hist_out_cama (file_hist_cama, itime_in_file_cama)
         ENDIF
#endif

#ifdef CatchLateralFlow
         CALL hist_basin_out (file_hist, idate)
#endif

#ifdef GridRiverLakeFlow
         CALL hist_grid_riverlake_out (file_hist, HistForm, idate, &
            itime_in_file, trim(file_hist)/=trim(file_last))

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               allocate (nac_one (numpatch))
               nac_one = 1.
            ENDIF
         ENDIF

         IF (p_is_io) CALL allocate_block_data (ghist, sumarea_one)
         IF (p_is_io) CALL flush_block_data (sumarea_one, 1.)

         CALL write_history_variable_2d ( DEF_hist_vars%riv_height, a_wdsrf_ucat_pch,   &
            file_hist, 'f_wdpth_ucat_regrid', itime_in_file, sumarea_ucat, filter_ucat, &
            'regridded deepest water depth in river and flood plain', 'm', nac_one)

         CALL write_history_variable_2d ( DEF_hist_vars%riv_veloct, a_veloc_riv_pch,    &
            file_hist, 'f_veloc_riv_regrid', itime_in_file, sumarea_ucat, filter_ucat,  &
            'regridded water velocity in river', 'm/s', nac_one)

         CALL write_history_variable_2d ( DEF_hist_vars%discharge, a_discharge_pch,     &
            file_hist, 'f_discharge', itime_in_file, sumarea_one, filter_ucat,          &
            'regridded discharge in river and flood plain', 'm^3/s',                    &
            nac_one, input_mode = 'total')

         CALL write_history_variable_2d ( DEF_hist_vars%discharge, a_dis_rmth_pch,      &
            file_hist, 'f_discharge_rivermouth_regrid', itime_in_file, sumarea_one,     &
            filter_ucat, 'regridded river mouth discharge into ocean', 'm^3/s',         &
            nac_one, input_mode = 'total')

         CALL write_history_variable_2d ( DEF_hist_vars%floodfrc, a_floodfrc_pch,       &
            file_hist, 'f_floodfrc', itime_in_file, sumarea_inpm, filter_inpm,          &
            'flooded area fraction', '100%', nac_one)

         IF (trim(HistForm) == 'Gridded') THEN
            CALL write_history_variable_2d ( DEF_hist_vars%floodarea, a_floodfrc_pch,   &
               file_hist, 'f_floodarea', itime_in_file, sumarea_one, filter_inpm,       &
               'flooded area', 'km^2', nac_one)
         ENDIF

         IF (allocated(nac_one   )) deallocate (nac_one   )
#endif

         IF (allocated(filter    )) deallocate (filter    )
         IF (allocated(filter_dt )) deallocate (filter_dt )
#ifdef URBAN_MODEL
         IF (allocated(filter_urb)) deallocate (filter_urb)
#endif

         CALL FLUSH_acc_fluxes ()

#ifdef SinglePoint
         IF (USE_SITE_HistWriteBack .and. memory_to_disk) THEN
            itime_mem = 0
         ENDIF
#endif

         file_last = file_hist

      ENDIF

   END SUBROUTINE hist_out


   SUBROUTINE write_history_variable_2d ( is_hist, &
         acc_vec, file_hist, varname, itime_in_file, sumarea, filter, &
         longname, units, acc_num, input_mode)

   USE MOD_Vars_1DAccFluxes, only: nac

   IMPLICIT NONE

   logical, intent(in) :: is_hist

   real(r8),         intent(inout) :: acc_vec(:)
   character(len=*), intent(in)    :: file_hist
   character(len=*), intent(in)    :: varname
   integer,          intent(in)    :: itime_in_file
   character(len=*), intent(in)    :: longname
   character(len=*), intent(in)    :: units

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)
   real(r8), intent(in), optional  :: acc_num(:)

   character(len=*), intent(in), optional :: input_mode

      IF (.not. is_hist) RETURN

#ifndef SinglePoint
      IF ( .not. present(acc_num) ) THEN
         IF (p_is_worker) &
            WHERE (acc_vec /= spval) acc_vec = acc_vec / nac
      ELSE
         IF (p_is_worker)  &
            WHERE (acc_vec/=spval .and. acc_num>0) acc_vec = acc_vec / acc_num
      ENDIF
#else
      IF ( .not. present(acc_num) ) THEN
         WHERE (acc_vec /= spval)  acc_vec = acc_vec / nac
      ELSE
         WHERE (acc_vec/=spval .and. acc_num>0) acc_vec = acc_vec / acc_num
      ENDIF
#endif

      select CASE (HistForm)
      CASE ('Gridded')
         IF (present(input_mode)) THEN
            CALL flux_map_and_write_2d ( &
               acc_vec, file_hist, varname, itime_in_file, sumarea, filter, longname, units, input_mode)
         ELSE
            CALL flux_map_and_write_2d ( &
               acc_vec, file_hist, varname, itime_in_file, sumarea, filter, longname, units)
         ENDIF
#if (defined UNSTRUCTURED || defined CATCHMENT)
      CASE ('Vector')
         IF (present(input_mode)) THEN
            CALL aggregate_to_vector_and_write_2d ( &
               acc_vec, file_hist, varname, itime_in_file, filter, longname, units, input_mode)
         ELSE
            CALL aggregate_to_vector_and_write_2d ( &
               acc_vec, file_hist, varname, itime_in_file, filter, longname, units)
         ENDIF
#endif
#ifdef SinglePoint
      CASE ('Single')
         CALL single_write_2d ( &
            acc_vec, file_hist, varname, itime_in_file, longname, units)
#endif
      END select

   END SUBROUTINE write_history_variable_2d


#ifdef URBAN_MODEL
   SUBROUTINE write_history_variable_urb_2d ( is_hist, &
         acc_vec, file_hist, varname, itime_in_file, sumarea, filter, &
         longname, units)

   IMPLICIT NONE

   logical, intent(in) :: is_hist

   real(r8), intent(inout) :: acc_vec(:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer,          intent(in) :: itime_in_file
   character(len=*), intent(in) :: longname
   character(len=*), intent(in) :: units

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)

      IF (.not. is_hist) RETURN

      select CASE (HistForm)
      CASE ('Gridded')
         CALL flux_map_and_write_urb_2d ( &
            acc_vec, file_hist, varname, itime_in_file, sumarea, filter, longname, units)
#if (defined UNSTRUCTURED || defined CATCHMENT)
      CASE ('Vector')
         !TODO: currently, it is not applicable to urban variables
         CALL aggregate_to_vector_and_write_2d ( &
            acc_vec, file_hist, varname, itime_in_file, filter, longname, units)
#endif
#ifdef SinglePoint
      CASE ('Single')
         CALL single_write_2d ( &
            acc_vec, file_hist, varname, itime_in_file, longname, units)
#endif
      END select

   END SUBROUTINE write_history_variable_urb_2d
#endif


   SUBROUTINE write_history_variable_3d ( is_hist, &
         acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, &
         sumarea, filter, longname, units)

   IMPLICIT NONE

   logical, intent(in) :: is_hist

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
   integer :: iblkme, xblk, yblk, xloc, yloc, i1
   integer :: compress

      IF (.not. is_hist) RETURN

      select CASE (HistForm)
      CASE ('Gridded')
         CALL flux_map_and_write_3d ( &
            acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, &
            sumarea, filter, longname, units)
#if (defined UNSTRUCTURED || defined CATCHMENT)
      CASE ('Vector')
         CALL aggregate_to_vector_and_write_3d ( &
            acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, &
            filter, longname, units)
#endif
#ifdef SinglePoint
      CASE ('Single')
         CALL single_write_3d (acc_vec, file_hist, varname, itime_in_file, &
            dim1name, ndim1, longname, units)
#endif
      END select

   END SUBROUTINE write_history_variable_3d


   SUBROUTINE write_history_variable_4d ( is_hist,   &
         acc_vec, file_hist, varname, itime_in_file, &
         dim1name, lb1, ndim1, dim2name, lb2, ndim2, &
         sumarea, filter, longname, units, acc_num)

   IMPLICIT NONE

   logical, intent(in) :: is_hist

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
   real(r8), intent(in), optional:: acc_num(:)

   ! Local variables
   integer :: i1, i2

      IF (.not. is_hist) RETURN

#ifndef SinglePoint
      IF ( .not. present(acc_num) ) THEN
         IF (p_is_worker) &
            WHERE (acc_vec /= spval) acc_vec = acc_vec / nac
      ELSE
         IF (p_is_worker) THEN
            DO i1 = lbound(acc_vec,1), ubound(acc_vec,1)
               DO i2 = lbound(acc_vec,2), ubound(acc_vec,2)
                  WHERE (acc_vec(i1,i2,:)/=spval .and. acc_num>0) &
                        acc_vec(i1,i2,:) = acc_vec(i1,i2,:) / acc_num
               ENDDO
            ENDDO
         ENDIF
      ENDIF
#else
      IF ( .not. present(acc_num) ) THEN
         WHERE (acc_vec /= spval)  acc_vec = acc_vec / nac
      ELSE
         DO i1 = lbound(acc_vec,1), ubound(acc_vec,1)
            DO i2 = lbound(acc_vec,2), ubound(acc_vec,2)
               WHERE (acc_vec(i1,i2,:)/=spval .and. acc_num>0) &
                     acc_vec(i1,i2,:) = acc_vec(i1,i2,:) / acc_num
            ENDDO
         ENDDO
      ENDIF
#endif

      select CASE (HistForm)
      CASE ('Gridded')
         CALL flux_map_and_write_4d ( &
            acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, &
            dim2name, lb2, ndim2, sumarea, filter, longname, units)
#if (defined UNSTRUCTURED || defined CATCHMENT)
      CASE ('Vector')
         CALL aggregate_to_vector_and_write_4d ( &
            acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, &
            dim2name, lb2, ndim2, filter, longname, units)
#endif
#ifdef SinglePoint
      CASE ('Single')
         CALL single_write_4d (acc_vec, file_hist, varname, itime_in_file, &
            dim1name, ndim1, dim2name, ndim2, longname, units)
#endif
      END select

   END SUBROUTINE write_history_variable_4d


   SUBROUTINE hist_write_time (filename, filelast, dataname, time, itime)

   IMPLICIT NONE

   character (len=*), intent(in) :: filename
   character (len=*), intent(in) :: filelast
   character (len=*), intent(in) :: dataname
   integer, intent(in)  :: time(3)
   integer, intent(out) :: itime

      select CASE (HistForm)
      CASE ('Gridded')
         CALL hist_gridded_write_time (filename, filelast, dataname, time, itime)
#if (defined UNSTRUCTURED || defined CATCHMENT)
      CASE ('Vector')
         CALL hist_vector_write_time  (filename, filelast, dataname, time, itime)
#endif
#ifdef SinglePoint
      CASE ('Single')
         CALL hist_single_write_time  (filename, filelast, dataname, time, itime)
#endif
      END select

   END SUBROUTINE hist_write_time

END MODULE MOD_Hist
