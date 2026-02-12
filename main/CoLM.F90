#include <define.h>

PROGRAM CoLM
!-----------------------------------------------------------------------------
!  Description:
!    This is the main program for the Common Land Model (CoLM)
!
!    Copyright © Yongjiu Dai Land Modeling Group at the School of Atmospheric Sciences
!    of the Sun Yat-sen University, Guangdong, CHINA.
!    All rights reserved.
!
!  Initial : Yongjiu Dai, 1998-2014
!  Revised : Hua Yuan, Shupeng Zhang, Nan Wei, Xingjie Lu, Zhongwang Wei, Yongjiu Dai
!            2014-2024
!-----------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Vars_Global
   USE MOD_Const_LC
   USE MOD_Const_PFT
   USE MOD_Const_Physical
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_TimeVariables
   USE MOD_Vars_1DForcing
   USE MOD_Vars_2DForcing
   USE MOD_Vars_1DFluxes
   USE MOD_Vars_1DAccFluxes
   USE MOD_Forcing
   USE MOD_Hist
   USE MOD_CheckEquilibrium
   USE MOD_TimeManager
   USE MOD_RangeCheck

   USE MOD_Block
   USE MOD_Pixel
   USE MOD_Mesh
   USE MOD_LandElm
#ifdef CATCHMENT
   USE MOD_LandHRU
#endif
   USE MOD_LandPatch
#ifdef URBAN_MODEL
   USE MOD_LandUrban
   USE MOD_Urban_LAIReadin
#endif
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_LandPFT
#endif
#if (defined UNSTRUCTURED || defined CATCHMENT)
   USE MOD_ElmVector
#endif
#ifdef CATCHMENT
   USE MOD_HRUVector
#endif
#if (defined CaMa_Flood)
   USE MOD_CaMa_colmCaMa
#endif
#ifdef SinglePoint
   USE MOD_SingleSrfdata
#endif
#if (defined CatchLateralFlow)
   USE MOD_Catch_BasinNetwork
   USE MOD_Catch_LateralFlow
#endif

   USE MOD_Ozone, only: init_ozone_data, update_ozone_data

   USE MOD_SrfdataRestart
   USE MOD_LAIReadin

#ifdef BGC
   USE MOD_NitrifData
   USE MOD_NdepData
   USE MOD_FireData
   USE MOD_LightningData
#endif

#ifdef CROP
   USE MOD_CropReadin
#endif

#ifdef LULCC
   USE MOD_Lulcc_Driver
#endif

#ifdef CoLMDEBUG
   USE MOD_Hydro_SoilWater
#endif

   ! SNICAR model
   USE MOD_SnowSnicar, only: SnowAge_init, SnowOptics_init
   USE MOD_Aerosol, only: AerosolDepInit, AerosolDepReadin

#ifdef DataAssimilation
   USE MOD_DA_Main
   USE MOD_DA_Vars_TimeVariables
   USE MOD_DA_Vars_1DFluxes
#endif

#ifdef USEMPI
   USE MOD_HistWriteBack
#endif

#ifdef EXTERNAL_LAKE
   USE MOD_Lake_Namelist
#endif

   IMPLICIT NONE

   character(len=256) :: nlfile
   character(len=256) :: casename
   character(len=256) :: dir_landdata
   character(len=256) :: dir_forcing
   character(len=256) :: dir_hist
   character(len=256) :: dir_restart
   character(len=256) :: fsrfdata

   real(r8) :: deltim       ! time step (seconds)
   integer  :: sdate(3)     ! calendar (year, julian day, seconds)
   integer  :: idate(3)     ! calendar (year, julian day, seconds)
   integer  :: edate(3)     ! calendar (year, julian day, seconds)
   integer  :: pdate(3)     ! calendar (year, julian day, seconds)
   integer  :: jdate(3)     ! calendar (year, julian day, seconds), year beginning style
   logical  :: greenwich    ! greenwich time

   logical :: doalb         ! true => start up the surface albedo calculation
   logical :: dolai         ! true => start up the time-varying vegetation parameter
   logical :: dosst         ! true => update sst/ice/snow

   integer :: Julian_1day_p, Julian_1day
   integer :: Julian_8day_p, Julian_8day
   integer :: s_year, s_month, s_day, s_seconds, s_julian
   integer :: e_year, e_month, e_day, e_seconds, e_julian
   integer :: p_year, p_month, p_day, p_seconds, p_julian
   integer :: lc_year, lai_year
   integer :: month, mday, year_p, month_p, mday_p, month_prev, mday_prev
   integer :: spinup_repeat, istep

   type(timestamp) :: ststamp, itstamp, etstamp, ptstamp, time_prev

   integer*8 :: start_time, end_time, c_per_sec, time_used
!-----------------------------------------------------------------------

#ifdef USEMPI
#ifdef USESplitAI
      integer :: num_procs, my_rank, ierr, color, new_comm
      CALL MPI_Init(ierr) ! Initialize MPI
      CALL MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr) ! Get the total number of processes
      CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr) ! Get the rank of the current process
      color = 1 ! The pyroot process will be in its own communicator
      print*, 'before split I am process', my_rank, 'of', num_procs
      CALL MPI_Comm_split(MPI_COMM_WORLD, color, my_rank, new_comm, ierr) ! Split the communicator
      print*, 'after split I am process', my_rank, 'of', num_procs
      CALL MPI_Comm_size(new_comm, num_procs, ierr) ! Get the total number of processes
      CALL MPI_Comm_rank(new_comm, my_rank, ierr) ! Get the rank of the current process
      print*,num_procs,"for CoLM"
      CALL spmd_init (new_comm)
#else
      CALL spmd_init ()
#endif
#endif

      CALL getarg (1, nlfile)

      CALL read_namelist (nlfile)

#ifdef EXTERNAL_LAKE
      CALL read_lake_namelist (nlfile)
#endif

#ifdef USEMPI
      IF (DEF_HIST_WriteBack) THEN
         CALL spmd_assign_writeback ()
      ENDIF

      IF (p_is_writeback) THEN
         CALL hist_writeback_daemon ()
      ELSE
#endif

      IF (p_is_master) THEN
         CALL system_clock (start_time)
      ENDIF

      casename     = DEF_CASE_NAME
      dir_landdata = DEF_dir_landdata
      dir_forcing  = DEF_dir_forcing
      dir_hist     = DEF_dir_history
      dir_restart  = DEF_dir_restart

#ifdef SinglePoint
      fsrfdata = trim(dir_landdata) // '/srfdata.nc'
#ifndef URBAN_MODEL
      CALL read_surface_data_single (fsrfdata, mksrfdata=.false.)
#else
      CALL read_urban_surface_data_single (fsrfdata, mksrfdata=.false., mkrun=.true.)
#endif
#endif

      deltim    = DEF_simulation_time%timestep
      greenwich = DEF_simulation_time%greenwich
      s_year    = DEF_simulation_time%start_year
      s_month   = DEF_simulation_time%start_month
      s_day     = DEF_simulation_time%start_day
      s_seconds = DEF_simulation_time%start_sec
      e_year    = DEF_simulation_time%end_year
      e_month   = DEF_simulation_time%end_month
      e_day     = DEF_simulation_time%end_day
      e_seconds = DEF_simulation_time%end_sec
      p_year    = DEF_simulation_time%spinup_year
      p_month   = DEF_simulation_time%spinup_month
      p_day     = DEF_simulation_time%spinup_day
      p_seconds = DEF_simulation_time%spinup_sec

      spinup_repeat = DEF_simulation_time%spinup_repeat

      CALL initimetype(greenwich)
      CALL monthday2julian(s_year,s_month,s_day,s_julian)
      CALL monthday2julian(e_year,e_month,e_day,e_julian)
      CALL monthday2julian(p_year,p_month,p_day,p_julian)

      sdate(1) = s_year; sdate(2) = s_julian; sdate(3) = s_seconds
      edate(1) = e_year; edate(2) = e_julian; edate(3) = e_seconds
      pdate(1) = p_year; pdate(2) = p_julian; pdate(3) = p_seconds

      CALL Init_GlobalVars
      CALL Init_LC_Const
      CALL Init_PFT_Const

#ifdef LULCC
      lc_year = s_year
      DEF_LC_YEAR = lc_year
#else
      lc_year = DEF_LC_YEAR
#endif

#ifndef SinglePoint
      CALL pixel%load_from_file    (dir_landdata)
      CALL gblock%load_from_file   (dir_landdata)

      CALL mesh_load_from_file (dir_landdata, lc_year)

      CALL pixelset_load_from_file (dir_landdata, 'landelm'  , landelm  , numelm  , lc_year)

#ifdef CATCHMENT
      CALL pixelset_load_from_file (dir_landdata, 'landhru'  , landhru  , numhru  , lc_year)
#endif

      CALL pixelset_load_from_file (dir_landdata, 'landpatch', landpatch, numpatch, lc_year)

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      CALL pixelset_load_from_file (dir_landdata, 'landpft'  , landpft  , numpft  , lc_year)
      CALL map_patch_to_pft
#endif

#ifdef URBAN_MODEL
      CALL pixelset_load_from_file (dir_landdata, 'landurban', landurban, numurban, lc_year)
      CALL map_patch_to_urban
#endif

#if (defined UNSTRUCTURED || defined CATCHMENT)
      CALL elm_vector_init ()
#ifdef CATCHMENT
      CALL hru_vector_init ()
#endif
#endif

#ifdef CatchLateralFlow
      CALL build_basin_network ()
#endif
#endif

      CALL adj2end(sdate)
      CALL adj2end(edate)
      CALL adj2end(pdate)

      ststamp = sdate
      etstamp = edate
      ptstamp = pdate

      ! date in beginning style
      jdate = sdate
      CALL adj2begin(jdate)

      IF (ptstamp <= ststamp) THEN
         spinup_repeat = 0
      ELSE
         spinup_repeat = max(1, spinup_repeat)
      ENDIF

      ! ----------------------------------------------------------------------
      ! Read in the model time invariant constant data
      CALL allocate_TimeInvariants ()
      CALL READ_TimeInvariants (lc_year, casename, dir_restart)

      ! Read in the model time varying data (model state variables)
      CALL allocate_TimeVariables  ()
      CALL READ_TimeVariables (jdate, lc_year, casename, dir_restart)

      ! Read in the model time varying data (model state variables) for ensemble
#ifdef DataAssimilation
      CALL allocate_TimeVariables_ens()
      CALL READ_TimeVariables_ens(jdate, lc_year, casename, dir_restart)
#endif

      ! Read in SNICAR optical and aging parameters
      IF (DEF_USE_SNICAR) THEN
         CALL SnowOptics_init( DEF_file_snowoptics ) ! SNICAR optical parameters
         CALL SnowAge_init( DEF_file_snowaging )     ! SNICAR aging   parameters
      ENDIF

      ! ----------------------------------------------------------------------
      doalb = .true.
      dolai = .true.
      dosst = .false.

      ! Initialize meteorological forcing data module
      CALL allocate_1D_Forcing ()
      CALL forcing_init (dir_forcing, deltim, ststamp, lc_year, etstamp)
      CALL allocate_2D_Forcing (gforc)

      ! Initialize history data module
      CALL hist_init (dir_hist)
      CALL allocate_1D_Fluxes ()
#ifdef DataAssimilation
      CALL allocate_1D_Fluxes_ens ()
#endif

      CALL CheckEqb_init ()

#if (defined CaMa_Flood)
#ifdef USEMPI
            CALL mpi_barrier (p_comm_glb, p_err)
#endif
      CALL colm_CaMa_init !initialize CaMa-Flood
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
#endif

      IF(DEF_USE_OZONEDATA)THEN
         CALL init_Ozone_data (sdate)
      ENDIF

      ! Initialize aerosol deposition forcing data
      IF (DEF_Aerosol_Readin) THEN
         CALL AerosolDepInit ()
      ENDIF

#ifdef BGC
      IF (DEF_USE_NITRIF) THEN
         CALL init_nitrif_data (ststamp)
      ENDIF

      IF (DEF_NDEP_FREQUENCY==1)THEN ! Initial annual ndep data readin
         CALL init_ndep_data_annually (sdate(1))
      ELSEIF(DEF_NDEP_FREQUENCY==2)THEN ! Initial monthly ndep data readin
         CALL init_ndep_data_monthly (sdate(1),s_month)
      ELSE
         write(6,*) 'ERROR: DEF_NDEP_FREQUENCY should be only 1-2, Current is:', &
                     DEF_NDEP_FREQUENCY
         CALL CoLM_stop ()
      ENDIF

      IF (DEF_USE_FIRE) THEN
         CALL init_fire_data (sdate(1))
         CALL init_lightning_data (sdate)
      ENDIF
#endif

#ifdef CROP
   CALL CROP_readin ()
#endif

#if (defined CatchLateralFlow)
      CALL lateral_flow_init (lc_year)
#endif

#ifdef DataAssimilation
      ! initialize data assimilation
      CALL init_DA ()
#endif

      ! ======================================================================
      ! begin time stepping loop
      ! ======================================================================

      istep   = 1
      idate   = sdate
      itstamp = ststamp

      TIMELOOP : DO WHILE (itstamp < etstamp)

         CALL julian2monthday (jdate(1), jdate(2), month_p, mday_p)

         year_p = jdate(1)

         IF (p_is_master) THEN
            IF (itstamp < ptstamp) THEN
               write(*, 99) istep, jdate(1), month_p, mday_p, jdate(3), spinup_repeat
            ELSE
               write(*,100) istep, jdate(1), month_p, mday_p, jdate(3)
            ENDIF
         ENDIF


         Julian_1day_p = int(calendarday(jdate)-1)/1*1 + 1
         Julian_8day_p = int(calendarday(jdate)-1)/8*8 + 1

         ! Read in the meteorological forcing
         ! ----------------------------------------------------------------------
         CALL read_forcing (jdate, dir_forcing)

         IF(DEF_USE_OZONEDATA)THEN
            CALL update_Ozone_data(itstamp, deltim)
         ENDIF

#ifdef BGC
         IF(DEF_USE_NITRIF) THEN
            time_prev = itstamp + int(-deltim)
            CALL julian2monthday(time_prev%year,time_prev%day,month_prev,mday_prev)
            if(month_p /= month_prev)then
               CALL update_nitrif_data (month_p)
            end if
         ENDIF
         IF(DEF_USE_FIRE)THEN
            CALL update_lightning_data (itstamp, deltim)
         ENDIF
#endif

         ! Read in aerosol deposition forcing data
         IF (DEF_Aerosol_Readin) THEN
            CALL AerosolDepReadin (jdate)
         ENDIF

         ! Calendar for NEXT time step
         ! ----------------------------------------------------------------------
         CALL TICKTIME (deltim,idate)
         itstamp = itstamp + int(deltim)
         jdate = idate
         CALL adj2begin(jdate)

         CALL julian2monthday (jdate(1), jdate(2), month, mday)

#ifdef BGC

         IF (DEF_NDEP_FREQUENCY==1)THEN ! Read Annual Ndep data
            IF (jdate(1) /= year_p) THEN
               CALL update_ndep_data_annually (idate(1), iswrite = .true.)
            ENDIF
         ELSEIF(DEF_NDEP_FREQUENCY==2)THEN! Read Monthly Ndep data
            IF (jdate(1) /= year_p .or. month /= month_p) THEN
               CALL update_ndep_data_monthly (jdate(1), month, iswrite = .true.)
            ENDIF
         ELSE
            write(6,*) 'ERROR: DEF_NDEP_FREQUENCY should be only 1-2, Current is:',&
                        DEF_NDEP_FREQUENCY
            CALL CoLM_stop ()
         ENDIF

         IF(DEF_USE_FIRE)THEN
            IF (jdate(1) /= year_p) THEN
               CALL update_hdm_data (idate(1))
            ENDIF
         ENDIF
#endif


         ! Call CoLM driver
         ! ----------------------------------------------------------------------
         IF (p_is_worker) THEN
#ifdef DataAssimilation
            CALL DADRIVER (idate, deltim, dolai, doalb, dosst, oroflag)
#else
            CALL CoLMDRIVER (idate,deltim,dolai,doalb,dosst,oroflag)
#endif
         ENDIF


#if (defined CatchLateralFlow)
         CALL lateral_flow (idate(1), deltim)
#endif
#if (defined CaMa_Flood)
#ifdef USEMPI
            CALL mpi_barrier (p_comm_glb, p_err)
#endif
         CALL colm_CaMa_drv(idate(3)) ! run CaMa-Flood
#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif
#endif
#ifdef DataAssimilation
         CALL run_DA (idate, deltim)
#endif

         ! Write out the model histroy file
         ! ----------------------------------------------------------------------
         CALL hist_out (idate, deltim, itstamp, etstamp, ptstamp, dir_hist, casename)

         CALL CheckEquilibrium (idate, deltim, itstamp, dir_hist, casename)

         ! DO land use and land cover change simulation
         ! ----------------------------------------------------------------------
#ifdef LULCC
         IF ( isendofyear(idate, deltim) .and. &
            ( jdate(1)>=2000 .or. (jdate(1)>1985 .and. MOD(jdate(1),5)==0) ) ) THEN

            ! Deallocate all Forcing and Fluxes variable of last year
            CALL deallocate_1D_Forcing
            CALL deallocate_1D_Fluxes
#ifdef DataAssimilation
            CALL deallocate_1D_Fluxes_ens
#endif

            CALL forcing_final ()
            CALL hist_final    ()

            ! Call LULCC driver
            CALL LulccDriver (casename, dir_landdata, dir_restart, jdate, greenwich)

            ! Allocate Forcing and Fluxes variable of next year
            CALL allocate_1D_Forcing
            CALL forcing_init (dir_forcing, deltim, itstamp, jdate(1), lulcc_call=.true.)

            CALL hist_init (dir_hist, lulcc_call=.true.)
            CALL allocate_1D_Fluxes
         ENDIF
#endif

         ! Get leaf area index
         ! ----------------------------------------------------------------------
#if (defined DYN_PHENOLOGY)
         ! Update once a day
         dolai = .false.
         Julian_1day = int(calendarday(jdate)-1)/1*1 + 1
         IF(Julian_1day /= Julian_1day_p)THEN
            dolai = .true.
         ENDIF
#else
         ! READ in Leaf area index and stem area index
         ! ----------------------------------------------------------------------
         ! NOTES: Should be caution for setting DEF_LAI_CHANGE_YEARLY to true in non-LULCC
         ! case, that means the LAI changes without consideration of land cover change.

         IF (DEF_LAI_CHANGE_YEARLY) THEN
            lai_year = jdate(1)
         ELSE
            lai_year = DEF_LC_YEAR
         ENDIF

         IF (DEF_LAI_MONTHLY) THEN
            IF (month /= month_p) THEN
               CALL LAI_readin (lai_year, month, dir_landdata)
#ifdef URBAN_MODEL
               CALL UrbanLAI_readin(lai_year, month, dir_landdata)
#endif
            ENDIF
         ELSE
            ! Update every 8 days (time interval of the MODIS LAI data)
            Julian_8day = int(calendarday(jdate)-1)/8*8 + 1
            IF (Julian_8day /= Julian_8day_p) THEN
               CALL LAI_readin (jdate(1), Julian_8day, dir_landdata)
            ENDIF
         ENDIF
#endif

         ! Write out the model state variables for restart run
         ! ----------------------------------------------------------------------
         IF (save_to_restart (idate, deltim, itstamp, ptstamp, etstamp)) THEN
#ifdef LULCC
            IF (jdate(1) >= 2000) THEN
               CALL WRITE_TimeVariables (jdate, jdate(1), casename, dir_restart)
            ELSE
               CALL WRITE_TimeVariables (jdate, (jdate(1)/5)*5, casename, dir_restart)
            ENDIF
#else
            CALL WRITE_TimeVariables (jdate, lc_year,  casename, dir_restart)
#endif
#ifdef DataAssimilation
            CALL WRITE_TimeVariables_ens (jdate, lc_year, casename, dir_restart)
#endif
#if (defined CaMa_Flood)
#ifdef USEMPI
            CALL mpi_barrier (p_comm_glb, p_err)
#endif
            IF (p_is_master) THEN
               CALL colm_cama_write_restart (jdate, lc_year,  casename, dir_restart)
            ENDIF
#ifdef USEMPI
            CALL mpi_barrier (p_comm_glb, p_err)
#endif
#endif
         ENDIF

#ifdef RangeCheck
         CALL check_TimeVariables ()
#ifdef DataAssimilation
         CALL check_TimeVariables_ens()
#endif
#endif

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CoLMDEBUG
         IF (DEF_USE_VariablySaturatedFlow) THEN
            CALL print_VSF_iteration_stat_info ()
         ENDIF
#endif

         IF (p_is_master) THEN
            CALL system_clock (end_time, count_rate = c_per_sec)
            time_used = (end_time - start_time) / c_per_sec
            IF (time_used >= 3600) THEN
               write(*,101) time_used/3600, mod(time_used,3600)/60, mod(time_used,60)
            ELSEIF (time_used >= 60) THEN
               write(*,102) time_used/60, mod(time_used,60)
            ELSE
               write(*,103) time_used
            ENDIF
         ENDIF

         IF ((spinup_repeat > 1) .and. (ptstamp <= itstamp)) THEN
            spinup_repeat = spinup_repeat - 1
            idate   = sdate
            jdate   = sdate
            itstamp = ststamp
            CALL adj2begin(jdate)
            CALL forcing_reset ()
         ENDIF

         istep = istep + 1

      ENDDO TIMELOOP

      CALL deallocate_TimeInvariants ()
      CALL deallocate_TimeVariables  ()
#ifdef DataAssimilation
      CALL deallocate_TimeVariables_ens()
      CALL end_DA()
#endif
      CALL deallocate_1D_Forcing     ()
      CALL deallocate_1D_Fluxes      ()
      CALL mesh_free_mem             ()

#if (defined CatchLateralFlow)
      CALL lateral_flow_final ()
#endif

      CALL forcing_final ()
      CALL hist_final    ()
      CALL CheckEqb_final()

#ifdef SinglePoint
      CALL single_srfdata_final ()
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#if (defined CaMa_Flood)
      CALL colm_cama_exit ! finalize CaMa-Flood
#endif

      IF (p_is_master) THEN
         write(*,'(/,A25)') 'CoLM Execution Completed.'
      ENDIF

      99  format(/, 'TIMESTEP = ', I0, ' | DATE = ', I4.4, '-', I2.2, '-', I2.2, '-', I5.5, &
          ' Spinup (', I0, ' repeat left)')
      100 format(/, 'TIMESTEP = ', I0, ' | DATE = ', I4.4, '-', I2.2, '-', I2.2, '-', I5.5)
      101 format(/, 'Time elapsed : ', I4, ' hours', I3, ' minutes', I3, ' seconds.')
      102 format(/, 'Time elapsed : ', I3, ' minutes', I3, ' seconds.')
      103 format(/, 'Time elapsed : ', I3, ' seconds.')

#ifdef USEMPI
      ENDIF

      IF (DEF_HIST_WriteBack) THEN
         CALL hist_writeback_exit ()
      ENDIF

      CALL spmd_exit
#endif

END PROGRAM CoLM
! ---------- EOP ------------
