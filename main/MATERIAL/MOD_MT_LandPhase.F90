#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_LandPhase

   USE MOD_Precision
   USE MOD_MT_LandWork, only: mt_land_work
   USE MOD_MT_Defs
   USE MOD_MT_HydroDiagnostics, only: mt_water_diagnostics
   USE MOD_MT_Solute, only: mt_solute_read_parameters
   USE MOD_MT_Vars
   USE MOD_MT_Precip
   USE MOD_MT_Evapo
   USE MOD_MT_SoilWater
   USE MOD_MT_Snow
   USE MOD_MT_Conservation
   USE MOD_MT_Rest
   USE MOD_MT_SoilInit
#ifdef GridRiverLakeFlow
   USE MOD_MT_RiverLakeCoupling, only: mt_rlc_land_flood_input, mt_rlc_land_publish_input
#endif
   USE MOD_SPMD_Task
   USE MOD_Namelist, only: DEF_USE_IRRIGATION, DEF_SPLIT_SOILSNOW, DEF_USE_Dynamic_Wetland
   USE MOD_MT_Lifecycle, only: material_lifecycle_init, &
      material_lifecycle_land_init, material_lifecycle_land_final, &
      material_lifecycle_lake_step, &
      material_lifecycle_wetland_decomp, material_lifecycle_soil_step, &
      material_lifecycle_land_report, material_lifecycle_land_flush_acc_fluxes, &
      material_lifecycle_land_accumulate_fluxes, material_lifecycle_land_read_restart

   IMPLICIT NONE
   PRIVATE

   ! MT_RegisterProviders.F90
   INTERFACE
      SUBROUTINE mt_register_providers ()
      END SUBROUTINE mt_register_providers
   END INTERFACE

   ! Land-side multi-phase MATERIAL facade. Keeps the implementation split
   ! across physics modules but exposes one import surface per host, so
   ! drivers do not USE every low-level material module directly. NOTE: this
   ! is a facade, not the engine -- per-step driving lives in the called
   ! modules. River/lake ROUTING materials are driven separately in HYDRO
   ! (MOD_Grid_RiverLakeFlow -> MOD_MT_Lifecycle) and are NOT re-exported
   ! here. PUBLIC declarations below are grouped by consuming host phase.

   ! -- Lifecycle: one-shot init/final, called from CoLM.F90 --------------
   PUBLIC :: mt_land_init, mt_land_final

   PUBLIC :: mt_land_work, mt_land_interception, mt_land_water_fluxes
   PUBLIC :: mt_land_thermal_step, mt_land_newsnow, mt_land_after_bgc
#ifdef CaMa_Flood
   PUBLIC :: mt_land_check_flood_feedback
#endif

   ! -- Column phase: per-patch soil-water/precip/evap/snow transport,
   !    called inside the CoLMMAIN patch loop (per ipatch, snl, nl_soil) --
   PUBLIC :: nmts, mt_tiny, mt_uses_generic_water
   PUBLIC :: mt_precip, mt_evapo, mt_land_water_step
   PUBLIC :: mt_flood_evap_loss
   PUBLIC :: mt_newsnow, mt_snowmelt_carry, mt_save_storage, mt_balance_check
   PUBLIC :: mt_apply_decay
   PUBLIC :: mt_wliq_soisno, mt_wice_soisno, mt_solid_soisno, mt_scv
   PUBLIC :: mt_ldew_rain, mt_ldew_snow, mt_sm_carry

   ! -- Element phase: per-step orchestration, called from CoLMDRIVER;
   !    these delegate to the methane family ------------------
   PUBLIC :: mt_resolve_step, mt_lake_step, mt_wetland_decomp
   PUBLIC :: mt_soil_step, mt_report

   ! -- Acc-flux: called from MOD_Vars_1DAccFluxes -----------------------
   PUBLIC :: mt_flush_acc_fluxes, mt_accumulate_fluxes

CONTAINS

   SUBROUTINE mt_land_thermal_step(work, ipatch, deltim, snl, nl_soil, &
      ldew, ldew_rain, ldew_snow, tleaf, wliq_soisno, wice_soisno, t_soisno, forc_q, forc_psrf, scv)
      type(mt_land_work), intent(inout) :: work
      integer, intent(in) :: ipatch, snl, nl_soil
      real(r8), intent(in) :: deltim, ldew, ldew_rain, ldew_snow, tleaf, forc_q, forc_psrf, scv
      real(r8), intent(in) :: wliq_soisno(snl+1:nl_soil), wice_soisno(snl+1:nl_soil)
      real(r8), intent(in) :: t_soisno(snl+1:nl_soil)

      CALL mt_canopy_parts(ldew, ldew_rain, ldew_snow, tleaf, &
         work%canopy_rain_postthermal, work%canopy_snow_postthermal)
      CALL mt_evapo(ipatch, deltim, snl, nl_soil, &
         work%canopy_rain_postthermal, work%canopy_snow_postthermal, &
         work%canopy_rain_prethermal, work%canopy_snow_prethermal, &
         wliq_soisno, wice_soisno, work%wliq_before, work%wice_before, &
         canopy_smelt_mass_th=work%thermal%canopy%melt, &
         canopy_frzc_mass_th=work%thermal%canopy%freeze, &
         soil_thaw_mass_th=work%thermal%soil_thaw, soil_frzc_mass_th=work%thermal%soil_freeze, &
         tleaf_frac=tleaf, t_soisno_frac=t_soisno, forc_q_frac=forc_q, forc_psrf_frac=forc_psrf)
      CALL mt_snowmelt_carry(ipatch, snl, scv, work%snow_before)
      CALL work%save_water(snl+1, nl_soil, wliq_soisno, wice_soisno)
   END SUBROUTINE mt_land_thermal_step

   SUBROUTINE mt_land_newsnow(work, ipatch, patchtype, snl, snl_bef, nl_soil, pg_snow, &
      deltim, scv, wetwat, wliq_soisno, wice_soisno)
      type(mt_land_work), intent(inout) :: work
      integer, intent(in) :: ipatch, patchtype, snl, snl_bef, nl_soil
      real(r8), intent(in) :: pg_snow, deltim, scv, wetwat
      real(r8), intent(in) :: wliq_soisno(snl+1:nl_soil), wice_soisno(snl+1:nl_soil)

      work%snow_before = scv
      IF (snl < 0) THEN
         CALL mt_newsnow(ipatch, patchtype, snl, snl_bef, pg_snow, deltim, &
            scv, work%snow_before, wetwat, wliq_soisno(snl+1:0), wice_soisno(snl+1:0), &
            work%snow_ice_before(snl+1:0))
      ELSE
         CALL mt_newsnow(ipatch, patchtype, snl, snl_bef, pg_snow, deltim, &
            scv, work%snow_before, wetwat)
      ENDIF
   END SUBROUTINE mt_land_newsnow

   SUBROUTINE mt_land_after_bgc(istep_in, ipatch, idate, deltim, patchtype)
#ifdef BGC
      USE MOD_MT_NitrogenLeaching, only: mt_bgc_n_leaching
#endif
      integer, intent(in), optional :: istep_in
      integer, intent(in) :: ipatch, idate(3), patchtype
      integer :: istep_local
      real(r8), intent(in) :: deltim

      CALL mt_resolve_step(istep_in, istep_local)
      CALL mt_wetland_decomp(ipatch, deltim)
      CALL mt_soil_step(istep_local, ipatch, idate, deltim)
#ifdef BGC
      CALL mt_bgc_n_leaching(ipatch, patchtype, deltim)
#endif
   END SUBROUTINE mt_land_after_bgc


   SUBROUTINE mt_land_interception (ipatch, deltim, snl, nl_soil, work, &
      ldew, ldew_rain, ldew_snow, tleaf, forc_rain, forc_snow, &
      qintr, qintr_rain, qintr_snow, pg_rain, pg_snow, qflx_irrig_sprinkler, &
      canopy_drain_rain, canopy_drain_snow)

      integer, intent(in) :: ipatch, snl, nl_soil
      type(mt_land_work), intent(inout) :: work
      real(r8), intent(in) :: deltim, ldew, ldew_rain, ldew_snow, tleaf
      real(r8), intent(in) :: forc_rain, forc_snow, qintr, qintr_rain, qintr_snow
      real(r8), intent(in) :: pg_rain, pg_snow, qflx_irrig_sprinkler
      real(r8), intent(in) :: canopy_drain_rain, canopy_drain_snow

      ! Interception parts derived from its outputs (rates, mm/s):
      ! throughfall = input - qintr_*, so old-canopy release xsc = pg_* - throughfall.
      work%rain_release = max(pg_rain - (forc_rain + qflx_irrig_sprinkler - qintr_rain), 0._r8)
      work%snow_release = max(pg_snow - (forc_snow - qintr_snow), 0._r8)
      ! qintr_* are NET interception (input - direct throughfall - drainage);
      ! mt_precip needs the GROSS interception (input - direct throughfall)
      ! to split pure-gap throughfall from post-mix canopy drip.  The leaf
      ! interception call above exposes the real drainage rate for MT builds
      ! (mt_leaf_qdrp_*); using the net value here lost all old-canopy solute
      ! flushing (N02).
      ! Old-water unloading can make net interception negative.
      work%gross_rain = max(qintr_rain + max(canopy_drain_rain, 0._r8), 0._r8)
      work%gross_snow = max(qintr_snow + max(canopy_drain_snow, 0._r8), 0._r8)
      CALL mt_canopy_lump (ipatch, ldew, ldew_rain, ldew_snow, tleaf, &
         work%canopy_rain_intercepted, work%canopy_snow_intercepted)
#ifdef CROP
      IF (DEF_USE_IRRIGATION) THEN
         CALL mt_save_storage(ipatch, snl, nl_soil, work%irrigation_storage_before)
      ELSE
         CALL mt_save_storage(ipatch, snl, nl_soil)
      ENDIF
#else
      CALL mt_save_storage(ipatch, snl, nl_soil)
#endif
      ! the host interception has no canopy phase change (melt/freeze masses are zero)
      CALL mt_precip(ipatch, deltim, &
         forc_rain, forc_snow, qintr, qintr_rain, qintr_snow, &
         pg_rain, pg_snow, work%canopy_rain_intercepted, work%canopy_snow_intercepted, &
         work%canopy_rain_before, work%canopy_snow_before, &
         qflx_irrig_sprinkler, &
         work%gross_rain, work%gross_snow, &
         work%rain_release, work%snow_release, &
         0._r8, 0._r8, work%irrigation_storage_before)
   END SUBROUTINE mt_land_interception

   SUBROUTINE mt_land_water_fluxes (work, water_diag, wa, qinfl, deltim)
      type(mt_land_work), intent(inout) :: work
      type(mt_water_diagnostics), intent(in) :: water_diag
      real(r8), intent(in) :: wa, qinfl, deltim

      ! Water recharging the aquifer, consistent with the root and baseflow
      ! parts taken from it.
      work%recharge = (wa - work%wa_before + water_diag%etroot_aquifer + &
         max(water_diag%rsub_source_aquifer, 0._r8)) / max(deltim, mt_tiny)
      ! The host qinfl also counts the ponding water that the aquifer
      ! exchange took for roots and baseflow; those parts are booked by
      ! etroot_surface and rsub_source_surface, so leave them out here.
      work%infiltration = qinfl - (water_diag%etroot_surface + water_diag%rsub_source_surface) / deltim

   END SUBROUTINE mt_land_water_fluxes

#ifdef CaMa_Flood
   SUBROUTINE mt_land_check_flood_feedback (patchtype, enabled)
      integer, intent(in) :: patchtype
      logical, intent(in) :: enabled
      integer :: imt_loc

      ! CaMa-Flood exchanges no materials with the land; stop instead of
      ! booking flood water without composition.
      IF (enabled .and. patchtype == 0 .and. nmts > 0) THEN
         DO imt_loc = 1, nmts
            IF (.not. mt_uses_generic_water(imt_loc)) CYCLE
            CALL CoLM_stop('Material transport: CaMa flood feedback (LWINFILT) is not supported.')
         ENDDO
      ENDIF
   END SUBROUTINE mt_land_check_flood_feedback
#endif

   SUBROUTINE mt_land_water_step (ipatch, deltim, snl, nl_soil, patchtype, work, water_diag, &
      qlayer, rsur, rsub, qseva, qsdew, qsubl, qfros, &
      qseva_soil, qsdew_soil, qsubl_soil, qfros_soil, &
      qseva_snow, qsdew_snow, qsubl_snow, qfros_snow, etr, sm, fsno, &
      wliq_soisno, wice_soisno, wa, wdsrf, wetwat, pg_rain, pg_snow, &
      t_soisno, porsl, dz_soisno, forc_us, forc_vs, forc_q, forc_psrf, &
      tleaf, lai, rst, rss, qinfl, qflx_irrig_drip, qflx_irrig_flood, qflx_irrig_paddy &
#ifdef GridRiverLakeFlow
      ,flood_credit_patch, fevpg_fld, qinfl_fld &
#endif
      )

      integer, intent(in) :: ipatch, snl, nl_soil, patchtype
      type(mt_land_work), intent(inout) :: work
      type(mt_water_diagnostics), intent(in) :: water_diag
      real(r8), intent(in) :: deltim, qlayer(0:nl_soil), rsur, rsub
      real(r8), intent(in) :: qseva, qsdew, qsubl, qfros
      real(r8), intent(in) :: qseva_soil, qsdew_soil, qsubl_soil, qfros_soil
      real(r8), intent(in) :: qseva_snow, qsdew_snow, qsubl_snow, qfros_snow
      real(r8), intent(in) :: etr, sm, fsno, wa, wdsrf, wetwat, pg_rain, pg_snow
      real(r8), intent(in) :: wliq_soisno(snl+1:nl_soil), wice_soisno(snl+1:nl_soil)
      real(r8), intent(in) :: t_soisno(snl+1:nl_soil), porsl(nl_soil)
      real(r8), intent(in) :: dz_soisno(min(snl+1,0):nl_soil)
      real(r8), intent(in) :: forc_us, forc_vs, forc_q, forc_psrf, tleaf, lai, rst, rss
      real(r8), intent(in) :: qinfl, qflx_irrig_drip, qflx_irrig_flood, qflx_irrig_paddy
#ifdef GridRiverLakeFlow
      real(r8), allocatable, intent(in) :: flood_credit_patch(:)
      real(r8), intent(in) :: fevpg_fld, qinfl_fld
      real(r8) :: flood_input_material(nmts)
#endif
      integer :: lbsn

      work%irrigation = qflx_irrig_drip + qflx_irrig_flood + qflx_irrig_paddy
      CALL mt_land_water_fluxes(work, water_diag, wa, qinfl, deltim)
#ifdef GridRiverLakeFlow
      CALL mt_rlc_land_flood_input(ipatch, patchtype, flood_credit_patch, &
         fevpg_fld, qinfl_fld, deltim, flood_input_material)
#endif
      lbsn = min(snl+1,0)
      IF (patchtype == 2 .and. .not. DEF_USE_Dynamic_Wetland) THEN
         CALL mt_wetland(ipatch, deltim, snl, nl_soil, &
            rsur, &
            qseva, qsdew, qsubl, qfros, &
            qseva_soil, qsdew_soil, qsubl_soil, qfros_soil, &
            qseva_snow, qsdew_snow, qsubl_snow, qfros_snow, &
            etr, sm, fsno, DEF_SPLIT_SOILSNOW, &
            wliq_soisno(snl+1:nl_soil), wice_soisno(snl+1:nl_soil), &
            work%wliq_before, work%wice_before, &
            wa, work%wa_before, wdsrf, work%wdsrf_before, &
            wetwat, work%wetwat_before, pg_rain, pg_snow, &
            t_soisno(snl+1:nl_soil), porsl(1:nl_soil), &
            dz_soisno(snl+1:nl_soil), &
            work%irrigation, &
            forc_us, forc_vs, work%irrigation_storage_ground, &
            snow_qout_layer = water_diag%snow_qout_layer(lbsn:0), &
            forc_q_frac = forc_q, &
            forc_psrf_frac = forc_psrf, &
            tleaf_frac = tleaf, &
            lai_frac = lai, &
            rst_frac = rst, &
            ra_frac = work%thermal%canopy%raw, &
            dz_sno_frac = dz_soisno(lbsn:0))
      ELSE
         CALL mt_soil_water(ipatch, deltim, snl, nl_soil, &
            qlayer, work%infiltration, work%recharge, rsur, rsub, &
            qseva, qsdew, qsubl, qfros, &
            qseva_soil, qsdew_soil, qsubl_soil, qfros_soil, &
            qseva_snow, qsdew_snow, qsubl_snow, qfros_snow, &
            sm, fsno, DEF_SPLIT_SOILSNOW, &
            wliq_soisno(snl+1:nl_soil), wice_soisno(snl+1:nl_soil), &
            work%wliq_before, work%wice_before, &
            wa, work%wa_before, wdsrf, work%wdsrf_before, &
            wetwat, work%wetwat_before, pg_rain, pg_snow, &
            water_diag%wblc_ice_sink, &
            water_diag%etroot_actual, water_diag%etroot_aquifer, &
            work%irrigation, &
            work%irrigation_storage_ground, &
            etroot_surface = water_diag%etroot_surface, dew_overflow = water_diag%dew_overflow, &
            frost_displaced = water_diag%frost_displaced, late_surface_runoff = water_diag%late_runoff, &
            rsub_source_layer = water_diag%rsub_source_layer, &
            rsub_source_surface = water_diag%rsub_source_surface, &
            rsub_source_aquifer = water_diag%rsub_source_aquifer, &
            imperv_evap_wdsrf = water_diag%imperv_evap_wdsrf, &
            imperv_evap_soil = water_diag%imperv_evap_soil, &
            imperv_subl_soil = water_diag%imperv_subl_soil, &
            snow_qout_layer = water_diag%snow_qout_layer(lbsn:0), &
            qgtop_solver = water_diag%qgtop, &
            tleaf_frac = tleaf, &
            t_soisno_frac = t_soisno(snl+1:nl_soil), &
            forc_q_frac = forc_q, &
            forc_psrf_frac = forc_psrf, &
            lai_frac = lai, &
            rst_frac = rst, &
            ra_frac = work%thermal%canopy%raw, &
            rss_frac = rss, &
            dz_soi_frac = dz_soisno(1:nl_soil), &
            porsl_frac = porsl(1:nl_soil), &
            dz_sno_frac = dz_soisno(lbsn:0), &
            permeable_soil = water_diag%permeable_soil &
#ifdef GridRiverLakeFlow
            ,flood_mt_input = flood_input_material, flood_infil_water = qinfl_fld*deltim &
#endif
            )
      ENDIF
#ifdef GridRiverLakeFlow
      CALL mt_rlc_land_publish_input(ipatch, patchtype, DEF_USE_Dynamic_Wetland, flood_input_material)
#endif
   END SUBROUTINE mt_land_water_step

   SUBROUTINE mt_land_init (numpatch, maxsnl, nl_soil, init_month, lc_year, jdate, &
      casename, dir_restart, dir_landdata, ldew_rain, ldew_snow, wliq_soisno, &
      wice_soisno, wa, wdsrf, wetwat, scv, waterstorage, loaded_restart, restart_file)

      IMPLICIT NONE
      integer, intent(in) :: numpatch, maxsnl, nl_soil
      integer, intent(in) :: init_month
      integer, intent(in) :: lc_year
      integer, intent(in) :: jdate(3)
      character(len=*), intent(in) :: casename
      character(len=*), intent(in) :: dir_restart
      character(len=*), intent(in) :: dir_landdata
      real(r8), allocatable, intent(in), optional :: ldew_rain(:)
      real(r8), allocatable, intent(in), optional :: ldew_snow(:)
      real(r8), allocatable, intent(in), optional :: wliq_soisno(:,:)
      real(r8), allocatable, intent(in), optional :: wice_soisno(:,:)
      real(r8), allocatable, intent(in), optional :: wa(:)
      real(r8), allocatable, intent(in), optional :: wdsrf(:)
      real(r8), allocatable, intent(in), optional :: wetwat(:)
      real(r8), allocatable, intent(in), optional :: scv(:)
      real(r8), allocatable, intent(in), optional :: waterstorage(:)
      logical, intent(out), optional :: loaded_restart
      character(len=*), intent(out), optional :: restart_file

      character(len=256) :: file_restart_mt
      character(len=14)  :: cdate_restart
      character(len=256) :: cyear_restart
      real(r8), allocatable :: material_dummy_patch(:)
      real(r8), allocatable :: material_dummy_soisno(:,:)
      logical :: have_patch_state
      logical :: have_waterstorage

      write(cyear_restart,'(i4.4)') lc_year
      write(cdate_restart,'(i4.4,"-",i3.3,"-",i5.5)') jdate(1), jdate(2), jdate(3)
      file_restart_mt = trim(dir_restart)//'/'//trim(cdate_restart)//'/'//trim(casename)// &
                         '_restart_'//trim(cdate_restart)//'_lc'//trim(cyear_restart)//'.nc'
      IF (present(restart_file)) restart_file = file_restart_mt
      IF (present(loaded_restart)) loaded_restart = .false.

      have_patch_state = present(ldew_rain) .and. present(ldew_snow) .and. &
         present(wliq_soisno) .and. present(wice_soisno) .and. present(wa) .and. &
         present(wdsrf) .and. present(wetwat) .and. present(scv)
      IF (have_patch_state) THEN
         have_patch_state = allocated(ldew_rain) .and. allocated(ldew_snow) .and. &
            allocated(wliq_soisno) .and. allocated(wice_soisno) .and. allocated(wa) .and. &
            allocated(wdsrf) .and. allocated(wetwat) .and. allocated(scv)
      ENDIF

      IF (have_patch_state) THEN
         have_waterstorage = DEF_USE_IRRIGATION .and. present(waterstorage)
         IF (have_waterstorage) have_waterstorage = allocated(waterstorage)
         IF (have_waterstorage) THEN
            CALL material_init_from_arrays (numpatch, maxsnl, nl_soil, &
               ldew_rain, ldew_snow, wliq_soisno, wice_soisno, &
               wa, wdsrf, wetwat, scv, file_restart_mt, waterstorage, &
               init_month=init_month, lc_year=lc_year, jdate=jdate, &
               casename=casename, dir_restart=dir_restart, dir_landdata=dir_landdata, &
               loaded_restart=loaded_restart)
         ELSE
            CALL material_init_from_arrays (numpatch, maxsnl, nl_soil, &
               ldew_rain, ldew_snow, wliq_soisno, wice_soisno, &
               wa, wdsrf, wetwat, scv, file_restart_mt, init_month=init_month, &
               lc_year=lc_year, jdate=jdate, casename=casename, &
               dir_restart=dir_restart, dir_landdata=dir_landdata, loaded_restart=loaded_restart)
         ENDIF
      ELSE
         ! Vector restart I/O is collective over IO/worker groups. Non-worker
         ! ranks do not own patch water arrays, but they must still enter the
         ! MATERIAL init path so vector restart reads can complete.
         allocate(material_dummy_patch(0))
         allocate(material_dummy_soisno(maxsnl+1:nl_soil, 0))
         CALL material_init_from_arrays (0, maxsnl, nl_soil, &
            material_dummy_patch, material_dummy_patch, &
            material_dummy_soisno, material_dummy_soisno, &
            material_dummy_patch, material_dummy_patch, material_dummy_patch, &
            material_dummy_patch, file_restart_mt, init_month=init_month, &
            lc_year=lc_year, jdate=jdate, casename=casename, &
            dir_restart=dir_restart, dir_landdata=dir_landdata, loaded_restart=loaded_restart)
         deallocate(material_dummy_patch, material_dummy_soisno)
      ENDIF

   END SUBROUTINE mt_land_init

   SUBROUTINE material_init_from_arrays (numpatch, maxsnl, nl_soil, &
      ldew_rain, ldew_snow, wliq_soisno, wice_soisno, &
      wa, wdsrf, wetwat, scv, file_restart, waterstorage, init_month, &
      lc_year, jdate, casename, dir_restart, dir_landdata, loaded_restart)

      IMPLICIT NONE
      integer,  intent(in) :: numpatch, maxsnl, nl_soil
      real(r8), intent(in) :: ldew_rain(numpatch)
      real(r8), intent(in) :: ldew_snow(numpatch)
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wa(numpatch)
      real(r8), intent(in) :: wdsrf(numpatch)
      real(r8), intent(in) :: wetwat(numpatch)
      real(r8), intent(in) :: scv(numpatch)
      character(len=*), intent(in), optional :: file_restart
      ! Optional: irrigation reservoir from MOD_Vars_TimeVariables.
      ! Passed by CoLM.F90 only under CROP+DEF_USE_IRRIGATION. Used for
      ! cold-start initialisation of mt_waterstorage; the restart path
      ! already handles persistence via read_land_mt_restart.
      real(r8), intent(in), optional :: waterstorage(numpatch)
      integer,  intent(in), optional :: init_month
      integer, intent(in), optional :: lc_year
      integer, intent(in), optional :: jdate(3)
      character(len=*), intent(in), optional :: casename
      character(len=*), intent(in), optional :: dir_restart
      character(len=*), intent(in), optional :: dir_landdata
      logical, intent(out), optional :: loaded_restart
      logical :: found_restart, scv_missing, waterstorage_missing
      integer :: soil_init_month

      ! Defensive: if a previous mt_land_init left state behind (e.g. a
      ! LULCC transition that rebuilt numpatch without an explicit
      ! mt_land_final call), release it so the re-allocation below uses
      ! the new numpatch. Without this the allocate() call would trip
      ! Fortran's "already allocated" runtime check, OR the save-level
      ! snap_* arrays in MOD_MT_Conservation would keep the old
      ! numpatch and then misalign with the new material arrays.
      CALL material_lifecycle_land_final ()
      CALL deallocate_Material_Vars()
      CALL deallocate_mt_conservation()
      CALL material_defs_init()
      CALL mt_solute_read_parameters()
      CALL material_defs_validate()
      CALL material_lifecycle_init (mt_register_providers)
      IF (present(loaded_restart)) loaded_restart = .false.
      IF (present(lc_year) .and. present(jdate) .and. present(casename) .and. &
          present(dir_restart) .and. present(dir_landdata)) THEN
         CALL material_lifecycle_land_init (numpatch, lc_year, jdate, casename, dir_restart, dir_landdata, file_restart)
      ENDIF
      IF (nmts <= 0) RETURN
      CALL allocate_Material_Vars(numpatch, maxsnl, nl_soil)
      found_restart = .false.
      scv_missing   = .false.
      waterstorage_missing = .false.
      IF (present(file_restart)) THEN
         CALL read_land_mt_restart(file_restart, maxsnl, nl_soil, &
            found_restart, scv_missing, waterstorage_missing, wa)
      ENDIF
#ifdef USEMPI
      ! The control/master rank does not own vector restart blocks, so
      ! read_land_mt_restart can only make a local decision there.
      ! Make restart usability global before any cold-start path that
      ! contains MPI collectives, including material soil initialisation.
      CALL mpi_allreduce(MPI_IN_PLACE, found_restart, 1, MPI_LOGICAL, &
         MPI_LAND, p_comm_glb, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, scv_missing, 1, MPI_LOGICAL, &
         MPI_LOR, p_comm_glb, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, waterstorage_missing, 1, MPI_LOGICAL, &
         MPI_LOR, p_comm_glb, p_err)
      scv_missing = scv_missing .and. found_restart
      waterstorage_missing = waterstorage_missing .and. found_restart
#endif
      IF (present(loaded_restart)) loaded_restart = found_restart
      mt_land_restart_loaded = found_restart .and. nmts > 0
      IF (.not. found_restart) THEN
         soil_init_month = 1
         IF (present(init_month)) soil_init_month = init_month
         IF (present(waterstorage)) THEN
            CALL mt_init_from_water(numpatch, maxsnl, nl_soil, &
               ldew_rain, ldew_snow, wliq_soisno, wice_soisno, &
               wa, wdsrf, wetwat, scv, waterstorage)
         ELSE
            CALL mt_init_from_water(numpatch, maxsnl, nl_soil, &
               ldew_rain, ldew_snow, wliq_soisno, wice_soisno, &
               wa, wdsrf, wetwat, scv)
         ENDIF
         CALL material_soil_init_from_file(soil_init_month, numpatch, maxsnl, nl_soil, &
            wliq_soisno, wice_soisno)
      ELSE
         IF (waterstorage_missing .and. present(waterstorage)) THEN
            CALL material_init_waterstorage_from_ratio(numpatch, maxsnl, nl_soil, &
               ldew_rain, ldew_snow, wliq_soisno, wice_soisno, &
               wa, wdsrf, wetwat, scv, waterstorage)
         ENDIF
         IF (scv_missing) THEN
            ! Old-format restart: every other pool loaded from file, only
            ! mt_scv needs to be rebuilt from the hydrology scv.
            CALL material_init_scv_from_water(numpatch, maxsnl, nl_soil, &
               wliq_soisno, wice_soisno, scv)
         ENDIF
      ENDIF
      ! Provider-owned state has its own restart schema and must not be gated
      ! by the generic land-water descriptor.  Each provider probes its own
      ! mandatory fields and no-ops when the checkpoint has no provider state.
      IF (present(file_restart)) CALL material_lifecycle_land_read_restart(file_restart)
      ! A committed restart already contains every dissolved/solid pool. Do not
      ! re-equilibrate it: that also mixes the distinct wetland/surface pools.
      ! Legacy/incompatible checkpoints take the cold-start branch above.
      IF (.not. found_restart) THEN
         IF (present(waterstorage)) THEN
            CALL material_enforce_solubility_from_water(numpatch, maxsnl, nl_soil, &
               ldew_rain, wliq_soisno, wa, wdsrf, wetwat, waterstorage)
         ELSE
            CALL material_enforce_solubility_from_water(numpatch, maxsnl, nl_soil, &
               ldew_rain, wliq_soisno, wa, wdsrf, wetwat)
         ENDIF
      ENDIF
   END SUBROUTINE material_init_from_arrays

   SUBROUTINE mt_resolve_step (istep_in, istep_local)

      IMPLICIT NONE
      integer, intent(in),  optional :: istep_in
      integer, intent(out)           :: istep_local

      istep_local = 1
      IF (present(istep_in)) istep_local = istep_in

   END SUBROUTINE mt_resolve_step

   SUBROUTINE mt_lake_step (istep_in, ipatch, idate, deltim_phy, isub, nsub)

      IMPLICIT NONE
      integer, intent(in), optional :: istep_in
      integer :: istep_local
      integer,  intent(in) :: ipatch
      integer,  intent(in) :: idate(3)
      real(r8), intent(in) :: deltim_phy
      integer,  intent(in) :: isub
      integer,  intent(in) :: nsub

      CALL mt_resolve_step(istep_in, istep_local)
      CALL material_lifecycle_lake_step (istep_local, ipatch, idate, deltim_phy, isub, nsub)

   END SUBROUTINE mt_lake_step

   SUBROUTINE mt_wetland_decomp (ipatch, deltim)

      IMPLICIT NONE
      integer, intent(in) :: ipatch
      real(r8), intent(in) :: deltim

      CALL material_lifecycle_wetland_decomp (ipatch, deltim)

   END SUBROUTINE mt_wetland_decomp

   SUBROUTINE mt_soil_step (istep_local, ipatch, idate, deltim)

      IMPLICIT NONE
      integer,  intent(in) :: istep_local
      integer,  intent(in) :: ipatch
      integer,  intent(in) :: idate(3)
      real(r8), intent(in) :: deltim

      CALL material_lifecycle_soil_step (istep_local, ipatch, idate, deltim)

   END SUBROUTINE mt_soil_step

   SUBROUTINE mt_report ()

      IMPLICIT NONE

      CALL material_balance_report ()
      CALL material_lifecycle_land_report ()

   END SUBROUTINE mt_report

   SUBROUTINE mt_flush_acc_fluxes ()

      IMPLICIT NONE

      CALL material_lifecycle_land_flush_acc_fluxes ()

   END SUBROUTINE mt_flush_acc_fluxes

   SUBROUTINE mt_accumulate_fluxes ()

      IMPLICIT NONE

      CALL material_lifecycle_land_accumulate_fluxes ()

   END SUBROUTINE mt_accumulate_fluxes

   SUBROUTINE mt_land_final ()
      CALL material_lifecycle_land_final ()
      CALL deallocate_Material_Vars()
      CALL deallocate_mt_conservation()
   END SUBROUTINE mt_land_final

END MODULE MOD_MT_LandPhase
#endif
