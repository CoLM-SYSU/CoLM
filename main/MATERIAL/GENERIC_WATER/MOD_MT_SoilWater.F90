#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_SoilWater

   USE, INTRINSIC :: IEEE_ARITHMETIC, only: ieee_is_finite
   USE MOD_Precision
   USE MOD_SPMD_Task, only: CoLM_stop
   USE MOD_Vars_Global, only: spval
   USE MOD_MT_Defs, only: nmts, materials, mt_tiny, mt_water_min_for_ratio, &
      mt_delta_sanity_max, material_init_water_ratio, material_is_isotope, &
      material_aquifer_isotope_state_valid, material_aquifer_actual_water, &
      material_aquifer_actual_mass, material_aquifer_isotope_ratio, mt_uses_generic_water, &
      material_is_nonvolatile_solute
   USE MOD_MT_Solute, only: mt_has_dissolved_limit, mt_equilibrate_dissolved
   USE MOD_MT_Forcing, only: material_forcing_rain_value, material_forcing_vapor_value, &
      material_forcing_has_vapor
   USE MOD_MT_Frac, only: material_fractionation_active, material_surface_relhum, &
      material_transpiration_nss_ratio, material_alpha_kinetic_craig_gordon, &
      material_alpha_kinetic_open_water, material_diffusivity_ratio_air, &
      material_alpha_kinetic_soil, &
      material_soil_effective_diffusivity, material_soil_diffusive_transfer, &
      material_isotope_soil_vapor_equivalent_diffusivity, &
      material_isotope_snow_vapor_equivalent_diffusivity, material_alpha_ice_vap, &
      material_alpha_ice_liq, material_equilibration_exchange, &
      material_liquid_self_diffusivity, &
      material_craig_gordon_evap_ratio, material_craig_gordon_evap_affine, material_alpha_liq_vap, &
      material_rayleigh_freezing_loss, material_equilibrium_deposition_ratio
   USE MOD_MT_Namelist, only: DEF_ISO_SUBL_SKIN_MM, DEF_ISO_SOIL_KINETIC, DEF_ISO_SOIL_DIFFUSION, &
      DEF_ISO_SOIL_VAPOR_DIFFUSION, DEF_ISO_SNOWMELT_EQUILIBRATION
   USE MOD_MT_EvapLimit, only: material_atmospheric_mt_loss
   USE MOD_MT_Isotope_SoilFlow, only: isotope_sweep_deficit_from_host_water, isotope_wetland_vapor_uptake
#ifdef MATERIAL_SIGNED_SOIL_DEBT_EXPERIMENT
   USE MOD_MT_Defs, only: material_has_decay
#endif
      USE MOD_MT_Vars, only: mt_wliq_soisno, mt_wice_soisno, &
         mt_wa, mt_aquifer_ref_water, mt_aquifer_ref_mass, &
         mt_wdsrf, mt_wetwat, mt_surface_residue, mt_subsurface_residue, &
         mt_waterstorage, mt_solid_soisno, mt_surface_solid, &
         mt_subsurface_solid, mt_waterstorage_solid, &
            a_mt_precip, a_water_precip, a_mt_transp_src, a_mt_evap, a_mt_transp, &
            a_water_transp, a_water_evap_gross, material_book_evap_loss, &
            MT_EVAP_KIND_TRANSP, MT_EVAP_KIND_SOILEVAP, MT_EVAP_KIND_SUBL, &
            MT_EVAP_KIND_WETLAND, &
            a_mt_qinfl, a_mt_qcharge, a_mt_rsur, a_mt_rsub, a_mt_rnof, &
         mt_pg_rain_ground, mt_rnof_step, mt_sm_carry, &
         mt_leaf_delta_e, mt_leaf_delta_b, mt_leaf_peclet, mt_leaf_water_moles, &
         mt_leaf_iso_storage, mt_numerical_residual_step, mt_numerical_water_step

   IMPLICIT NONE

CONTAINS

   SUBROUTINE check_isotope_aquifer (imt, ipatch, water_mass, isotope_mass, stage)
      integer, intent(in) :: imt, ipatch
      real(r8), intent(in) :: water_mass, isotope_mass
      character(len=*), intent(in) :: stage

      IF (.not. material_is_isotope(imt)) RETURN
      IF (material_aquifer_isotope_state_valid(water_mass, isotope_mass, &
            materials(imt)%ref_ratio, mt_aquifer_ref_water(ipatch), &
            mt_aquifer_ref_mass(imt, ipatch))) RETURN
      WRITE(*,'(2A,2(A,I0),2(A,ES16.8))') &
         'Unresolved isotope aquifer debt at ', trim(stage), &
         ' material=', imt, ' patch=', ipatch, ' water=', water_mass, ' isotope=', isotope_mass
      CALL CoLM_stop('isotope aquifer debt has no physical water carrier')
   END SUBROUTINE check_isotope_aquifer

   ! Material leaving the aquifer with `water` mm, `flux` being water times its
   ! concentration as computed by the caller.  Only the part of the
   ! withdrawal that the aquifer actually holds (`actual_water`) is capped by
   ! the actual material mass; the rest is the host driving wa below that store
   ! (an aquifer debt) and carries the same concentration, as when the aquifer
   ! is empty.  Capping the whole withdrawal by the mass of a dust-level store
   ! (e.g. wa=0 plus a 1e-14 mm root return on a patch without mixing
   ! reference water) exported almost no material with a full TOPMODEL
   ! baseflow.  When actual_water >= water this is exactly the previous
   ! min(flux, mass).
   pure real(r8) FUNCTION aquifer_export_material (water, flux, actual_water, actual_mass)
      real(r8), intent(in) :: water, flux, actual_water, actual_mass
      real(r8) :: covered

      IF (actual_water <= mt_tiny) THEN
         aquifer_export_material = flux
      ELSEIF (actual_water >= water) THEN
         aquifer_export_material = min(flux, max(actual_mass, 0._r8))
      ELSE
         covered = actual_water / water
         aquifer_export_material = min(flux * covered, max(actual_mass, 0._r8)) &
            + flux * (1._r8 - covered)
      ENDIF
   END FUNCTION aquifer_export_material

   !---------------------------------------------------------------
   ! Flux-driven material update after WATER.
   !
   ! Surface water: mixed-pool approach (throughfall + old wdsrf + soil upflow)
   ! Soil layers: qlayer-driven (exact inter-layer flux tracking)
   ! Aquifer: qcharge-driven
   ! Subsurface runoff: pre-WATER ratio
   ! Freeze/thaw: delta-based (same-layer internal transfer)
   !
   ! WATER_VSF supplies resolved qlayer and ET/aquifer diagnostics.  WATER_2014
   ! supplies its Campbell/Richards interface fluxes from the tridiagonal
   ! soilwater solve; any residual storage mismatch is still reconciled below.
   !---------------------------------------------------------------
   ! REQUIRES: DEF_USE_VariablySaturatedFlow = .true.
   ! The soil-water material closure uses VSF-only diagnostics (aquifer exchange,
   ! actual root uptake, and wblc ice sink).  MATERIAL runs must be rejected at
   ! namelist validation time when DEF_USE_VariablySaturatedFlow is false.
   SUBROUTINE mt_soil_water (ipatch, deltim, snl, nl_soil, &
      qlayer, qinfl, qcharge, rsur, rsub, &
      qseva_in, qsdew_in, qsubl_in, qfros_in, &
      qseva_soil, qsdew_soil, qsubl_soil, qfros_soil, &
      qseva_snow, qsdew_snow, qsubl_snow, qfros_snow, &
      sm, fsno, split_soilsnow, &
      wliq_soisno, wice_soisno, &
      wliq_soisno_bef, wice_soisno_bef, &
      wa, wa_bef, wdsrf, wdsrf_bef, &
      wetwat, wetwat_bef, pg_rain, pg_snow, &
      wblc_ice_sink, &
      etroot_actual, etroot_aquifer, &
      qflx_irrig_ground, waterstorage_patch, &
      imperv_evap_wdsrf, imperv_evap_soil, imperv_subl_soil, &
      snow_qout_layer, qgtop_solver, tleaf_frac, t_soisno_frac, forc_q_frac, forc_psrf_frac, lai_frac, rst_frac, ra_frac, &
      rss_frac, dz_soi_frac, porsl_frac, dz_sno_frac, flood_mt_input, flood_infil_water, &
      etroot_surface, dew_overflow, frost_displaced, late_surface_runoff, &
      rsub_source_layer, rsub_source_surface, rsub_source_aquifer, permeable_soil)

      IMPLICIT NONE
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: deltim
      integer,  intent(in) :: snl, nl_soil
      real(r8), intent(in) :: qlayer(0:nl_soil)     ! inter-layer flux [mm/s]
      real(r8), intent(in) :: qinfl                 ! diagnostic infiltration flux [mm/s]
      real(r8), intent(in) :: qcharge               ! groundwater recharge [mm/s]
      real(r8), intent(in) :: rsur, rsub            ! runoff [mm/s]
      real(r8), intent(in) :: qseva_in              ! total evaporation (no _soil suffix) [mm/s]
      real(r8), intent(in) :: qsdew_in              ! total dew [mm/s]
      real(r8), intent(in) :: qsubl_in              ! total sublimation [mm/s]
      real(r8), intent(in) :: qfros_in              ! total frost [mm/s]
      real(r8), intent(in) :: qseva_soil            ! soil-only evaporation [mm/s]
      real(r8), intent(in) :: qsdew_soil            ! soil-only dew [mm/s]
      real(r8), intent(in) :: qsubl_soil            ! soil-only sublimation [mm/s]
      real(r8), intent(in) :: qfros_soil            ! soil-only frost [mm/s]
      real(r8), intent(in) :: qseva_snow            ! snow-only evaporation [mm/s]
      real(r8), intent(in) :: qsdew_snow            ! snow-only dew [mm/s]
      real(r8), intent(in) :: qsubl_snow            ! snow-only sublimation [mm/s]
      real(r8), intent(in) :: qfros_snow            ! snow-only frost [mm/s]
      real(r8), intent(in) :: sm                    ! snowmelt [mm/s]
      real(r8), intent(in) :: fsno                  ! snow fraction [-]
      logical,  intent(in) :: split_soilsnow        ! DEF_SPLIT_SOILSNOW
      real(r8), intent(in) :: wliq_soisno(snl+1:nl_soil)     ! post-WATER
      real(r8), intent(in) :: wice_soisno(snl+1:nl_soil)     ! post-WATER
      real(r8), intent(in) :: wliq_soisno_bef(snl+1:nl_soil) ! pre-WATER (post-THERMAL)
      real(r8), intent(in) :: wice_soisno_bef(snl+1:nl_soil) ! pre-WATER
      real(r8), intent(in) :: wa, wa_bef, wdsrf, wdsrf_bef
      real(r8), intent(in) :: wetwat, wetwat_bef
      real(r8), intent(in) :: pg_rain, pg_snow
      ! Per-layer ice mass (kg/m2) withdrawn by WATER_VSF's wblc
      ! reconciliation, which offsets the solver's excess liquid (no net
      ! water leaves the column). Zero for wetland branches that skip wblc.
      real(r8), intent(in) :: wblc_ice_sink(1:nl_soil)
      ! Per-layer actual ET water (mm, not mm/s) after the deficit
      ! cascade in MOD_Hydro_SoilWater.F90:soil_water_vertical_movement,
      ! and the aquifer share absorbing any remaining deficit (mm).
      ! Material path uses these instead of raw pre-cascade demand values,
      ! so that a shallow dry layer reports 0 (its demand cascaded to
      ! a deeper layer, which now reports that extra uptake) and the
      ! aquifer deduction is applied via mt_wa. Without this, dry
      ! layers would get spurious material removal clipped by min() and
      ! deep-cascade / aquifer ET would silently leave material behind.
      real(r8), intent(in) :: etroot_actual(1:nl_soil)
      real(r8), intent(in) :: etroot_aquifer

      ! Combined drip + flood + paddy irrigation rate (mm/s). WATER_VSF
      ! adds this to gwat at MOD_SoilSnowHydrology.F90:828-832 AFTER the
      ! qseva subtraction, so post-WATER wdsrf + rsur*dt + qlayer(0)*dt
      ! implicitly contains the irrigation water. The material path must
      ! inject the matching R*dt into the mixed surface pool or the
      ! downstream wdsrf/rsur/qlayer(0) distribution is under-supplied
      ! and the material budget shows a spurious step-level dilution.
      ! Sprinkler is handled separately in mt_precip (merged into
      ! the rain/canopy path by LEAF_INTERCEPTION). 0 when irrigation
      ! is disabled.
      real(r8), intent(in) :: qflx_irrig_ground
      real(r8), intent(in), optional :: waterstorage_patch
      ! WATER_VSF has a special impermeable-top-layer branch where qgtop<0
      ! atmospheric loss is removed directly from wdsrf and then soil layer 1
      ! before soil_water_vertical_movement. The soil-layer loss is phase-split
      ! by WATER into liquid evaporation and ice sublimation; neither appears as
      ! negative qinfl, so pass the exact water-side removals (mm/step) and
      ! mirror them explicitly in material mass.
      real(r8), intent(in), optional :: imperv_evap_wdsrf
      real(r8), intent(in), optional :: imperv_evap_soil
      real(r8), intent(in), optional :: imperv_subl_soil
      real(r8), intent(in), optional :: snow_qout_layer(snl+1:0)
      real(r8), intent(in), optional :: qgtop_solver
      real(r8), intent(in), optional :: tleaf_frac
      real(r8), intent(in), optional :: t_soisno_frac(snl+1:nl_soil)
      real(r8), intent(in), optional :: forc_q_frac
      real(r8), intent(in), optional :: forc_psrf_frac
      real(r8), intent(in), optional :: lai_frac
      real(r8), intent(in), optional :: rst_frac
      real(r8), intent(in), optional :: ra_frac
      ! Soil surface (pore) resistance [s/m].  With ra_frac it selects the
      ! resistance-weighted soil kinetic law, so the kinetic effect grows as
      ! the evaporation front retreats into the pores.
      real(r8), intent(in), optional :: rss_frac
      ! Soil layer thickness [m] and porosity [-], needed by the liquid-phase
      ! isotope diffusion term.  Optional so third-party callers that do not
      ! supply them simply run without diffusion.
      real(r8), intent(in), optional :: dz_soi_frac(1:nl_soil)
      real(r8), intent(in), optional :: porsl_frac(1:nl_soil)
      ! Snow layer thickness [m], needed by firn vapour diffusion.  Separate
      ! from dz_soi_frac because snow layer count varies with snl.
      real(r8), intent(in), optional :: dz_sno_frac(snl+1:0)
      ! External river-flood infiltration already included in qinfl [mm/s].
      ! This is its matching dissolved material amount [R*mm] for this step.
      real(r8), intent(in), optional :: flood_mt_input(:)
      real(r8), intent(in), optional :: flood_infil_water
      real(r8), intent(in), optional :: etroot_surface, dew_overflow, frost_displaced, late_surface_runoff ! [mm]
      real(r8), intent(in), optional :: rsub_source_layer(1:nl_soil), rsub_source_surface, rsub_source_aquifer
      ! Exact Richards connectivity, not inferred from zero reported qlayer.
      ! Absent for older callers: leave layer residuals on the explicit ledger.
      logical, intent(in), optional :: permeable_soil(1:nl_soil)

      integer  :: imt, j, lb, lb_snow
      real(r8) :: R_precip, R_atm
      real(r8) :: mt_flux, ratio, ratio_src, storage_ratio
      real(r8) :: qcharge_eff
      real(r8) :: d_wice, d_wetwat
      real(r8) :: mt_pool_total, water_pool_total
      real(r8) :: mt_surface_collapse_residual
         real(r8) :: ratio_layer(1:nl_soil)  ! pre-WATER material ratio per layer
         real(r8) :: water_shadow(1:nl_soil)
         real(r8) :: remap_face_water(1:nl_soil-1), remap_trial_water(1:nl_soil)
         real(r8) :: layer_transport_ratio(1:nl_soil)
         real(r8) :: water_resid, water_shadow_ratio
      real(r8) :: mt_soil_upflow         ! material from true soil-to-surface upflow
      real(r8) :: mt_soil_evap           ! material removed by top-boundary evaporation deficit
      real(r8) :: mt_upflow_deferred     ! upflow material layer 1 pays after the qlayer transfers
      real(r8) :: evap_pool_water, evap_deficit_mt   ! layer-1 evaporation beyond its liquid
      real(r8) :: mt_wliq_qlayer_bef(1:nl_soil), water_shadow_qlayer_bef(1:nl_soil)
      logical  :: sweep_ok(1:nl_soil)
      real(r8) :: qlayer_through(1:nl_soil)
      real(r8) :: evap_base_water         ! surface water excluding this step's exfiltration
      integer  :: qlayer_pass
      logical  :: qlayer_clamp
      real(r8) :: qgtop_est
      ! Snow top-layer effective flux selectors
      real(r8) :: eff_qseva_snow, eff_qsdew_snow, eff_qsubl_snow, eff_qfros_snow
      ! Cross-phase deficit variables
      real(r8) :: water_ice_pool, water_liq_pool  ! post-deposit pool sizes
      real(r8) :: subl_water, evap_water           ! water to remove
      real(r8) :: deficit_water                    ! excess beyond primary phase
      real(r8) :: water_ice_pool_prefrost          ! pre-subl ice pool for deficit check
      real(r8) :: gwat_evap              ! evaporation subtracted from gwat [mm]
      real(r8) :: mt_gwat_evap         ! corresponding material removed
      real(r8) :: eff_qseva             ! effective evap used in gwat
      real(r8) :: eff_qsdew_topliq      ! effective dew on soil layer 1 liquid
      real(r8) :: eff_qsubl_top         ! effective sublimation on soil layer 1 ice
      real(r8) :: eff_qfros_top         ! effective frost on soil layer 1 ice
      real(r8) :: pool_ratio            ! actual ratio of mixed pool
      real(r8) :: top_infil_water       ! surface-to-soil infiltration water [mm/step]
      real(r8) :: flood_water, surface_base_water, surface_base_balance
      real(r8) :: flood_ground_evap_water, flood_ground_evap_material
      real(r8) :: late_water, late_material, late_ratio, late_runoff_water, early_runoff_water, late_scale
      real(r8) :: early_runoff_mt       ! early runoff drawn from the ordinary surface pool [mm/step]
      real(r8) :: flood_overflow_water   ! pond overflow of the late flood pool [mm/step]
      real(r8) :: flood_destination_water
      real(r8) :: top_exfil_water       ! soil-to-surface exfiltration water [mm/step]
      real(r8) :: top_boundary_out_water ! total negative qinfl water [mm/step]
      real(r8) :: top_soil_evap_water   ! negative qinfl caused by qseva deficit [mm/step]
      real(r8) :: imperv_wdsrf_loss      ! qgtop<0 loss from surface storage [mm/step]
      real(r8) :: imperv_soil_loss       ! qgtop<0 loss from layer-1 liquid [mm/step]
      real(r8) :: imperv_subl_loss       ! qgtop<0 loss from layer-1 ice [mm/step]
      real(r8) :: mt_gwat_snow         ! material leaving snow bottom [kg/m2 per step]
      real(r8) :: snow_rain_input       ! rain entering snowwater [mm]
      real(r8) :: snow_dew_input        ! dew entering top snow liquid [mm]
      real(r8) :: snow_evap_output      ! evap removed from top snow liquid [mm]
      real(r8) :: snow_frost_input      ! frost added to top snow ice [mm]
      real(r8) :: snow_subl_output      ! sublimation removed from top snow ice [mm]
      real(r8) :: snow_wliq_before_flow ! top snow liquid before percolation [mm]
      real(r8) :: qin_snow, qout_snow   ! reconstructed snowwater qin/qout [mm]
      real(r8) :: mt_qin_snow, mt_qout_snow
      real(r8) :: water_before_flow, mt_before_flow
      real(r8) :: wliq_pre_phase, wice_pre_phase   ! phase-change basis
      ! Actual external wice change on top snow layer after deficit handling:
      ! nominal (qfros - qsubl)*dt understates removal when qsubl
      ! exceeds post-frost ice, because the deficit is taken from liquid
      ! instead. Section 5b uses this to isolate internal freeze/thaw.
      real(r8) :: d_wice_ext_snow
      ! Actual external wice change on soil layer 1 after the
      ! max(0., wice_bef + (qfros-qsubl)*dt) clamp in
      ! MOD_SoilSnowHydrology.F90:1101. Symmetric to d_wice_ext_snow for
      ! top snow layer: when qsubl*dt exceeds wice_bef+qfros*dt,
      ! the water side floors wice at 0, so the true external change is
      ! -wice_bef rather than the nominal (qfros-qsubl)*dt. Without this
      ! correction, Section 5b's d_wice residual becomes spuriously
      ! positive and triggers a phantom internal freeze, moving liquid
      ! material into the already-empty ice slot.
      real(r8) :: d_wice_ext_soil1
      real(r8) :: wice_soil1_after_imperv
      real(r8) :: transp_water_total, xylem_mt_total, xylem_ratio
      real(r8) :: root_return_water, root_return_material, root_return_mt_total
      real(r8) :: root_gross_water, root_gross_material, return_ratio
      real(r8) :: root_return_excess, excess_ratio
      real(r8) :: aquifer_ratio
      real(r8) :: aquifer_water_pre_qcharge, aquifer_orphan_mass
      real(r8) :: aquifer_ref_water, aquifer_ref_mass, aquifer_actual_mass
      real(r8) :: surface_et_water, surface_root_return, surface_et_ratio
      real(r8) :: dew_surface_water, frost_surface_water
      real(r8) :: late_surface_water, pending_surface_material, late_surface_ratio
      real(r8) :: rsub_donor_water, rsub_donor_ratio
      logical :: resolved_rsub
      real(r8) :: transp_source_mt_total, transp_output_material
      real(r8) :: transp_ratio, remove_ratio, R_vapor
      ! True only around the exposed-soil-surface evaporation calls, so
      ! evap_cg_inputs applies the resistance-weighted soil kinetic law there
      ! without touching snow-layer or surface-water evaporation.
      logical  :: kinetic_on_soil_surface
      real(r8) :: diff_liquid, d_eff_up, d_eff_dn, diff_transfer
      real(r8) :: diff_face(snl+1:nl_soil), diff_out(snl+1:nl_soil), diff_scale(snl+1:nl_soil)
      real(r8) :: ratio_up, ratio_dn
      real(r8) :: melt_exchange
      real(r8) :: source_fallback_ratio
      real(r8) :: relhum_leaf, leaf_area_use, rst_use, ra_use, tleaf_use
      real(r8) :: leaf_delta_e_new, leaf_delta_b_new, leaf_peclet_new, leaf_moles_new
      real(r8) :: soil_resid_mt, soil_resid_water
      logical  :: transp_frac_active
#ifdef MATERIAL_SIGNED_SOIL_DEBT_EXPERIMENT
      logical :: signed_debt_step
      real(r8) :: soil_debt_limit
#endif

      IF (nmts <= 0) RETURN
#ifdef MATERIAL_SIGNED_SOIL_DEBT_EXPERIMENT
      ! shortcut: direct closed-advection experiment only; full-model debt needs lifecycle support.
      soil_debt_limit = 1.e-3_r8/(2._r8*real(nl_soil,r8))
      IF (.not. all(ieee_is_finite(mt_wliq_soisno(:,1:nl_soil,ipatch)))) &
         CALL CoLM_stop('signed soil debt: rejected nonfinite inventory')
      IF (.not. all(ieee_is_finite(wliq_soisno)) .or. &
          .not. all(ieee_is_finite(wliq_soisno_bef))) &
         CALL CoLM_stop('signed soil debt: rejected nonfinite carrier')
      signed_debt_step = any(wliq_soisno_bef(1:nl_soil) < 0._r8) .or. &
         any(wliq_soisno(1:nl_soil) < 0._r8) .or. &
         any(mt_wliq_soisno(:,1:nl_soil,ipatch) < 0._r8)
      IF (signed_debt_step) THEN
         CALL check_signed_soil_experiment()
         CALL signed_soil_advection()
         RETURN
      ENDIF
#endif
      lb = snl + 1
      kinetic_on_soil_surface = .false.
      qcharge_eff = qcharge
      flood_water = 0._r8
      IF (present(flood_infil_water)) flood_water = max(flood_infil_water, 0._r8)
      surface_et_water = 0._r8
      surface_root_return = 0._r8
      dew_surface_water = 0._r8
      frost_surface_water = 0._r8
      IF (present(etroot_surface)) THEN
         surface_et_water = max(etroot_surface, 0._r8)
         surface_root_return = max(-etroot_surface, 0._r8)
      ENDIF
      IF (present(dew_overflow)) dew_surface_water = max(dew_overflow, 0._r8)
      IF (present(frost_displaced)) frost_surface_water = max(frost_displaced, 0._r8)
      late_surface_water = dew_surface_water + frost_surface_water
      late_runoff_water = 0._r8
      IF (present(late_surface_runoff)) late_runoff_water = max(late_surface_runoff, 0._r8)
      early_runoff_water = max(rsur,0._r8)*deltim - late_runoff_water
      IF (early_runoff_water < -1.e-9_r8) CALL CoLM_stop('late runoff exceeds total surface runoff')
      early_runoff_water = max(early_runoff_water,0._r8)
      flood_destination_water = max(wdsrf,0._r8)+late_runoff_water-late_surface_water &
         + max(qinfl,0._r8)*deltim
      flood_destination_water = max(flood_destination_water,0._r8)
      resolved_rsub = .false.
      IF (present(rsub_source_layer) .and. present(rsub_source_surface) .and. present(rsub_source_aquifer)) THEN
         resolved_rsub = sum(rsub_source_layer) + rsub_source_surface + rsub_source_aquifer > mt_tiny
      ENDIF
      IF (flood_water > 0._r8 .and. .not. present(flood_mt_input)) &
         CALL CoLM_stop('grid flood material: water input has no material composition')
      IF (.not. ieee_is_finite(qcharge_eff) .or. &
          abs(qcharge_eff) > 0.5_r8 * abs(spval)) qcharge_eff = 0._r8

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         aquifer_ref_water = 0._r8
         aquifer_ref_mass = 0._r8
         IF (material_is_isotope(imt)) THEN
            aquifer_ref_water = mt_aquifer_ref_water(ipatch)
            aquifer_ref_mass = mt_aquifer_ref_mass(imt, ipatch)
         ENDIF
         ! Liquid proxy for the irrigation storage ratio (see below).
         R_precip = material_forcing_rain_value(imt, ipatch)
         R_atm = material_forcing_vapor_value(imt, ipatch)
         soil_resid_mt = 0._r8
         soil_resid_water = 0._r8
         water_shadow(1:nl_soil) = wliq_soisno_bef(1:nl_soil)

         IF (material_is_nonvolatile_solute(imt) .and. &
             wa_bef > mt_water_min_for_ratio .and. &
             mt_subsurface_residue(imt, ipatch) > mt_tiny) THEN
            mt_wa(imt, ipatch) = mt_wa(imt, ipatch) + &
               mt_subsurface_residue(imt, ipatch)
            mt_subsurface_residue(imt, ipatch) = 0._r8
         ENDIF
         CALL mt_equilibrate_dissolved(imt, wa_bef, mt_wa(imt, ipatch), &
            mt_subsurface_solid(imt, ipatch))
         DO j = lb, nl_soil
            CALL mt_equilibrate_dissolved(imt, max(wliq_soisno_bef(j), 0._r8), &
               mt_wliq_soisno(imt, j, ipatch), mt_solid_soisno(imt, j, ipatch))
         ENDDO

         ! The NSS leaf pool stores only an isotope anomaly, not water.
         ! When leaves disappear, return that signed anomaly to the
         ! root-zone liquid pools before computing this step's ratios.
         IF (present(lai_frac)) THEN
            IF (allocated(mt_leaf_iso_storage)) THEN
               IF (lai_frac <= mt_tiny .or. &
                   mt_leaf_water_moles(imt, ipatch) <= mt_tiny) THEN
                  CALL release_leaf_iso_storage(imt, ipatch, nl_soil, &
                     wliq_soisno_bef(1:nl_soil), wa_bef)
               ENDIF
            ENDIF
         ENDIF
         CALL check_isotope_aquifer(imt, ipatch, wa_bef, mt_wa(imt, ipatch), 'soil start')

         ! ============================================================
         ! 0. Compute pre-WATER material ratios for all soil layers.
         ! Dry-source fallbacks must not use current precipitation forcing:
         ! IsoGSM normalized ratios are noisy when precipitation is tiny,
         ! and root/aquifer uptake should inherit existing source water.
         ! ============================================================
         source_fallback_ratio = material_init_water_ratio(imt)
         DO j = nl_soil, 1, -1
            IF (wliq_soisno_bef(j) > mt_water_min_for_ratio) THEN
               source_fallback_ratio = mt_wliq_soisno(imt, j, ipatch) / wliq_soisno_bef(j)
               EXIT
            ENDIF
         ENDDO
         DO j = 1, nl_soil
            IF (wliq_soisno_bef(j) > mt_water_min_for_ratio) THEN
               ratio_layer(j) = mt_wliq_soisno(imt, j, ipatch) / wliq_soisno_bef(j)
            ELSEIF (material_is_nonvolatile_solute(imt)) THEN
               ! Dry solute inventory is immobile until this layer has a
               ! physically resolved liquid carrier; never borrow another
               ! layer's concentration to export it.
               ratio_layer(j) = 0._r8
            ELSE
               ratio_layer(j) = source_fallback_ratio
            ENDIF
         ENDDO

            ! ============================================================
            ! 0a. Transpiration per soil layer.
         !
         !    WATER_VSF's soil_water_vertical_movement removes water for
         !    transpiration with a cascade: a dry shallow layer's unmet
         !    demand propagates to deeper layers, and any residual is
         !    absorbed by the aquifer via soilwater_aquifer_exchange
         !    (MOD_Hydro_SoilWater.F90:323-351). `etroot_actual(j)` is the
         !    mm of water actually drawn from layer j after that cascade;
         !    `etroot_aquifer` is the mm absorbed by the aquifer. Using
         !    these — instead of raw pre-cascade demand — matches the
         !    water side layer-by-layer and lets mt_wa absorb the deep
         !    ET material that would otherwise be silently orphaned.
         !    No-fractionation assumption: ratio_layer stays unchanged
         !    after this step, so section 2's qlayer transfer can still
         !    use the same ratio_layer.
            ! ============================================================
               aquifer_ratio = source_fallback_ratio
               IF (abs(material_aquifer_actual_water(wa_bef, aquifer_ref_water)) > mt_water_min_for_ratio) THEN
                  CALL check_isotope_aquifer(imt, ipatch, wa_bef, mt_wa(imt, ipatch), 'root uptake')
                  aquifer_ratio = material_aquifer_isotope_ratio(wa_bef, mt_wa(imt, ipatch), &
                     aquifer_ref_water, aquifer_ref_mass, source_fallback_ratio)
               ELSEIF (material_is_nonvolatile_solute(imt)) THEN
                  aquifer_ratio = 0._r8
               ENDIF
               transp_water_total = max(etroot_aquifer, 0._r8)
               xylem_mt_total = transp_water_total * aquifer_ratio
               surface_et_ratio = source_fallback_ratio
               IF (wdsrf_bef > mt_water_min_for_ratio) THEN
                  surface_et_ratio = mt_wdsrf(imt, ipatch) / wdsrf_bef
               ELSEIF (material_is_nonvolatile_solute(imt)) THEN
                  surface_et_ratio = 0._r8
               ENDIF
               transp_water_total = transp_water_total + surface_et_water
               xylem_mt_total = xylem_mt_total + surface_et_water * surface_et_ratio
            DO j = 1, nl_soil
               IF (etroot_actual(j) > mt_tiny) THEN
                  transp_water_total = transp_water_total + etroot_actual(j)
                  xylem_mt_total = xylem_mt_total + etroot_actual(j) * ratio_layer(j)
               ENDIF
            ENDDO
               xylem_ratio = source_fallback_ratio
               IF (transp_water_total > mt_tiny) xylem_ratio = xylem_mt_total / transp_water_total
               root_gross_water = transp_water_total
               root_return_water = -sum(min(etroot_actual, 0._r8)) &
                  + max(-etroot_aquifer, 0._r8) + surface_root_return
               root_return_excess = max(root_return_water - root_gross_water, 0._r8)
               transp_water_total = max(root_gross_water - root_return_water, 0._r8)
               root_gross_material = 0._r8
                  transp_ratio = xylem_ratio
               transp_source_mt_total = 0._r8
               transp_output_material = 0._r8
               transp_frac_active = transp_water_total > mt_tiny .and. material_fractionation_active(imt)

            DO j = 1, nl_soil
               IF (etroot_actual(j) > mt_tiny .and. wliq_soisno_bef(j) > mt_tiny) THEN
                  IF (.not. material_is_nonvolatile_solute(imt)) THEN
                     mt_flux = etroot_actual(j) * ratio_layer(j)
                     mt_flux = min(mt_flux, max(mt_wliq_soisno(imt, j, ipatch), 0._r8))
                     mt_wliq_soisno(imt, j, ipatch) = mt_wliq_soisno(imt, j, ipatch) - mt_flux
                     root_gross_material = root_gross_material + mt_flux
                     IF (transp_frac_active) THEN
                        transp_source_mt_total = transp_source_mt_total + mt_flux
                     ELSE
                        CALL material_book_evap_loss(imt, ipatch, mt_flux, &
                           etroot_actual(j), MT_EVAP_KIND_TRANSP)
#if 0
                           a_mt_evap(imt, ipatch) = a_mt_evap(imt, ipatch) + mt_flux
#endif
                        a_mt_transp_src(imt, ipatch) = a_mt_transp_src(imt, ipatch) + mt_flux
                     ENDIF
                  ENDIF
                  water_shadow(j) = water_shadow(j) - etroot_actual(j)
               ENDIF
            ENDDO

         ! Aquifer share of ET (deficit cascade residual). Use SIGNED
         ! wa_bef so a pre-existing aquifer debt (wa_bef<0) still yields
         ! a well-defined ratio; the min/max clamp prevents mt_wa from
         ! going more negative than its current signed value.
         IF (etroot_aquifer > mt_tiny) THEN
            IF (.not. material_is_nonvolatile_solute(imt)) THEN
               remove_ratio = aquifer_ratio
               aquifer_actual_mass = material_aquifer_actual_mass(mt_wa(imt, ipatch), aquifer_ref_mass)
               mt_flux = aquifer_export_material(etroot_aquifer, etroot_aquifer * remove_ratio, &
                  material_aquifer_actual_water(wa_bef, aquifer_ref_water), aquifer_actual_mass)
               mt_wa(imt, ipatch) = mt_wa(imt, ipatch) - mt_flux
               root_gross_material = root_gross_material + mt_flux
               IF (transp_frac_active) THEN
                  transp_source_mt_total = transp_source_mt_total + mt_flux
               ELSE
                  CALL material_book_evap_loss(imt, ipatch, mt_flux, &
                     etroot_aquifer, MT_EVAP_KIND_TRANSP)
                  a_mt_transp_src(imt, ipatch) = a_mt_transp_src(imt, ipatch) + mt_flux
               ENDIF
            ENDIF
         ENDIF

         IF (surface_et_water > mt_tiny .and. .not. material_is_nonvolatile_solute(imt)) THEN
            mt_flux = min(surface_et_water * surface_et_ratio, &
               max(mt_wdsrf(imt, ipatch), 0._r8))
            mt_wdsrf(imt, ipatch) = mt_wdsrf(imt, ipatch) - mt_flux
            root_gross_material = root_gross_material + mt_flux
            IF (transp_frac_active) THEN
               transp_source_mt_total = transp_source_mt_total + mt_flux
            ELSE
               CALL material_book_evap_loss(imt, ipatch, mt_flux, surface_et_water, MT_EVAP_KIND_TRANSP)
               a_mt_transp_src(imt, ipatch) = a_mt_transp_src(imt, ipatch) + mt_flux
            ENDIF
         ENDIF

         ! Use the material mass actually removed from positive root donors;
         ! a finite-pool cap must never make reverse flow mint isotope mass.
         return_ratio = xylem_ratio
         excess_ratio = 0._r8
         IF (root_return_water > mt_tiny .and. .not. material_is_nonvolatile_solute(imt)) THEN
            IF (present(tleaf_frac)) THEN
               excess_ratio = deposition_ratio_for(tleaf_frac, .false.)
            ELSE
               excess_ratio = deposition_ratio_for(layer_temp(1), .false.)
            ENDIF
            IF (root_gross_water > mt_tiny .and. root_return_excess <= 0._r8) THEN
               return_ratio = root_gross_material/root_gross_water
            ELSEIF (root_gross_water > mt_tiny) THEN
               return_ratio = (root_gross_material*(root_return_water-root_return_excess)/root_gross_water &
                  + root_return_excess*excess_ratio) / root_return_water
            ELSE
               return_ratio = excess_ratio
            ENDIF
            xylem_ratio = return_ratio
            transp_ratio = return_ratio
         ENDIF
               IF (transp_frac_active .and. &
                   present(tleaf_frac) .and. present(forc_q_frac) .and. present(forc_psrf_frac) .and. &
                   present(lai_frac) .and. present(rst_frac) .and. allocated(mt_leaf_delta_e)) THEN
               tleaf_use = tleaf_frac
               leaf_area_use = lai_frac
               rst_use = rst_frac
               ra_use = 0._r8
               IF (present(ra_frac)) ra_use = max(ra_frac, 0._r8)
               IF (material_forcing_has_vapor(imt, ipatch)) THEN
                  R_vapor = R_atm
               ELSE
                  R_vapor = xylem_ratio / max(material_alpha_liq_vap(imt, tleaf_use), mt_tiny)
               ENDIF
               relhum_leaf = material_surface_relhum(forc_q_frac, forc_psrf_frac, tleaf_use, .false.)
               CALL material_transpiration_nss_ratio(imt, xylem_ratio, R_vapor, &
                  tleaf_use, relhum_leaf, forc_psrf_frac, transp_water_total, deltim, leaf_area_use, ra_use, rst_use, &
                  mt_leaf_delta_e(imt, ipatch), mt_leaf_delta_b(imt, ipatch), &
                  mt_leaf_peclet(imt, ipatch), mt_leaf_water_moles(imt, ipatch), &
                  transp_ratio, leaf_delta_e_new, leaf_delta_b_new, leaf_peclet_new, leaf_moles_new)
               mt_leaf_delta_e(imt, ipatch) = leaf_delta_e_new
               mt_leaf_delta_b(imt, ipatch) = leaf_delta_b_new
               mt_leaf_peclet(imt, ipatch) = leaf_peclet_new
               mt_leaf_water_moles(imt, ipatch) = leaf_moles_new
            ENDIF

         ! Negative plant-hydraulic root flux is xylem water returned to a
         ! soil layer, not an atmospheric input. The gross positive uptake
         ! above supplies its mixed isotope composition; excluded solutes
         ! remain in their donor soil layers and return with zero solute.
         IF (root_return_water > mt_tiny) THEN
            root_return_mt_total = 0._r8
            DO j = 1, nl_soil
               IF (etroot_actual(j) >= -mt_tiny) CYCLE
               water_shadow(j) = water_shadow(j) - etroot_actual(j)
               IF (material_is_nonvolatile_solute(imt)) CYCLE
               root_return_material = -etroot_actual(j) * return_ratio
               mt_wliq_soisno(imt, j, ipatch) = mt_wliq_soisno(imt, j, ipatch) + root_return_material
               root_return_mt_total = root_return_mt_total + root_return_material
            ENDDO
            IF (etroot_aquifer < -mt_tiny .and. .not. material_is_nonvolatile_solute(imt)) THEN
               root_return_material = -etroot_aquifer*return_ratio
               mt_wa(imt,ipatch) = mt_wa(imt,ipatch) + root_return_material
               root_return_mt_total = root_return_mt_total + root_return_material
            ENDIF
            IF (surface_root_return > mt_tiny .and. .not. material_is_nonvolatile_solute(imt)) THEN
               root_return_material = surface_root_return*return_ratio
               mt_wdsrf(imt,ipatch) = mt_wdsrf(imt,ipatch) + root_return_material
               root_return_mt_total = root_return_mt_total + root_return_material
            ENDIF
            IF (.not. material_is_nonvolatile_solute(imt)) THEN
               IF (root_return_mt_total > root_gross_material + root_return_excess*excess_ratio + &
                   max(1.e-12_r8, 1.e-9_r8*abs(return_ratio))) &
                  CALL CoLM_stop('plant hydraulic isotope return exceeds actual donor isotope')
               IF (transp_frac_active) THEN
                  transp_source_mt_total = transp_source_mt_total - root_return_mt_total
               ELSE
                  a_mt_evap(imt,ipatch) = a_mt_evap(imt,ipatch) - root_return_mt_total
                  a_mt_transp(imt,ipatch) = a_mt_transp(imt,ipatch) - root_return_mt_total
                  a_mt_transp_src(imt,ipatch) = a_mt_transp_src(imt,ipatch) - root_return_mt_total
                  a_water_transp(imt,ipatch) = a_water_transp(imt,ipatch) - root_return_water
                  a_water_evap_gross(imt,ipatch) = a_water_evap_gross(imt,ipatch) - root_return_water
               ENDIF
            ENDIF
         ENDIF

         CALL check_isotope_aquifer(imt, ipatch, wa_bef-etroot_aquifer, &
            mt_wa(imt, ipatch), 'after aquifer root exchange')

         ! Baseflow leaves from the same pond/soil/aquifer donors used by
         ! the VSF exchange, not always from the deepest soil layer.
         IF (resolved_rsub) THEN
            DO j = 1, nl_soil
               rsub_donor_water = max(rsub_source_layer(j), 0._r8)
               IF (rsub_donor_water <= mt_tiny) CYCLE
               mt_flux = min(rsub_donor_water * current_liq_ratio(j), &
                  max(mt_wliq_soisno(imt,j,ipatch), 0._r8))
               mt_wliq_soisno(imt,j,ipatch) = mt_wliq_soisno(imt,j,ipatch) - mt_flux
               water_shadow(j) = water_shadow(j) - rsub_donor_water
               a_mt_rsub(imt,ipatch) = a_mt_rsub(imt,ipatch) + mt_flux
               a_mt_rnof(imt,ipatch) = a_mt_rnof(imt,ipatch) + mt_flux
               mt_rnof_step(imt,ipatch) = mt_rnof_step(imt,ipatch) + mt_flux
            ENDDO
            rsub_donor_water = max(rsub_source_surface, 0._r8)
            IF (rsub_donor_water > mt_tiny) THEN
               rsub_donor_ratio = surface_et_ratio
               IF (wdsrf_bef-surface_et_water+surface_root_return > mt_water_min_for_ratio) THEN
                  rsub_donor_ratio = mt_wdsrf(imt,ipatch)/ &
                     (wdsrf_bef-surface_et_water+surface_root_return)
               ENDIF
               mt_flux = min(rsub_donor_water * rsub_donor_ratio, &
                  max(mt_wdsrf(imt,ipatch), 0._r8))
               mt_wdsrf(imt,ipatch) = mt_wdsrf(imt,ipatch) - mt_flux
               a_mt_rsub(imt,ipatch) = a_mt_rsub(imt,ipatch) + mt_flux
               a_mt_rnof(imt,ipatch) = a_mt_rnof(imt,ipatch) + mt_flux
               mt_rnof_step(imt,ipatch) = mt_rnof_step(imt,ipatch) + mt_flux
            ENDIF
            rsub_donor_water = max(rsub_source_aquifer, 0._r8)
            IF (rsub_donor_water > mt_tiny) THEN
               rsub_donor_ratio = aquifer_ratio
               CALL check_isotope_aquifer(imt, ipatch, wa_bef-etroot_aquifer, &
                  mt_wa(imt,ipatch), 'before aquifer baseflow')
               IF (abs(material_aquifer_actual_water(wa_bef-etroot_aquifer, &
                                                   aquifer_ref_water)) > mt_water_min_for_ratio) &
                  rsub_donor_ratio = material_aquifer_isotope_ratio(wa_bef-etroot_aquifer, &
                     mt_wa(imt,ipatch), aquifer_ref_water, aquifer_ref_mass, aquifer_ratio)
               mt_flux = aquifer_export_material(rsub_donor_water, rsub_donor_water * rsub_donor_ratio, &
                  material_aquifer_actual_water(wa_bef-etroot_aquifer, aquifer_ref_water), &
                  material_aquifer_actual_mass(mt_wa(imt,ipatch), aquifer_ref_mass))
               mt_wa(imt,ipatch) = mt_wa(imt,ipatch) - mt_flux
               CALL check_isotope_aquifer(imt, ipatch, &
                  wa_bef-etroot_aquifer-rsub_donor_water, mt_wa(imt,ipatch), &
                  'after aquifer baseflow')
               a_mt_rsub(imt,ipatch) = a_mt_rsub(imt,ipatch) + mt_flux
               a_mt_rnof(imt,ipatch) = a_mt_rnof(imt,ipatch) + mt_flux
               mt_rnof_step(imt,ipatch) = mt_rnof_step(imt,ipatch) + mt_flux
            ENDIF
         ENDIF

                  IF (transp_frac_active) THEN
                     a_mt_transp_src(imt, ipatch) = a_mt_transp_src(imt, ipatch) &
                        + transp_source_mt_total
                     ! NSS changes only the transpired vapor signature. Root uptake
                  ! still removes source/xylem water; the source-vapor difference
                  ! is held as leaf isotope isostorage to close material mass.
                  transp_output_material = transp_water_total * transp_ratio
                     CALL material_book_evap_loss(imt, ipatch, transp_output_material, &
                        transp_water_total, MT_EVAP_KIND_TRANSP)
                  IF (allocated(mt_leaf_iso_storage)) THEN
                     mt_leaf_iso_storage(imt, ipatch) = mt_leaf_iso_storage(imt, ipatch) &
                        + transp_source_mt_total - transp_output_material
                  ENDIF
               ENDIF

            ! ============================================================
         ! 0b. Snow layer external fluxes (when snow exists)
         !     snowwater applies qseva_snow/qsdew_snow/qsubl_snow/qfros_snow
         !     to the top snow layer. Track these as system I/O.
         ! ============================================================
         IF (snl < 0) THEN
            lb_snow = snl + 1  ! top snow layer index

            ! Select flux variables based on split path
            IF (split_soilsnow) THEN
               eff_qseva_snow = qseva_snow
               eff_qsdew_snow = qsdew_snow
               eff_qsubl_snow = qsubl_snow
               eff_qfros_snow = qfros_snow
            ELSE
               eff_qseva_snow = qseva_in
               eff_qsdew_snow = qsdew_in
               eff_qsubl_snow = qsubl_in
               eff_qfros_snow = qfros_in
            ENDIF

            ! Step 1: Ice phase — frost deposits, sublimation removes
            ! snowwater does: wice(lb) += (qfros - qsubl)*dt
            ! If wice goes negative, deficit removes from wliq.
            mt_flux = max(eff_qfros_snow, 0._r8) * &
               deposition_ratio_for(layer_temp(lb_snow), .true.) * deltim
            mt_wice_soisno(imt, lb_snow, ipatch) = mt_wice_soisno(imt, lb_snow, ipatch) + mt_flux
            a_mt_precip(imt, ipatch) = a_mt_precip(imt, ipatch) + mt_flux
            a_water_precip(imt, ipatch) = a_water_precip(imt, ipatch) + max(eff_qfros_snow, 0._r8) * deltim

            ! Post-frost ice pool (water and material)
            water_ice_pool = wice_soisno_bef(lb_snow) + max(eff_qfros_snow, 0._r8) * deltim
            water_ice_pool_prefrost = water_ice_pool  ! save pre-subl value for Step 2 deficit check
            subl_water = max(eff_qsubl_snow, 0._r8) * deltim

            IF (subl_water > mt_tiny) THEN
               IF (water_ice_pool - subl_water > mt_water_min_for_ratio) THEN
                  ! Normal case: enough ice to cover sublimation
                     mt_flux = atmospheric_loss_material(mt_wice_soisno(imt, lb_snow, ipatch), &
                        water_ice_pool, subl_water, layer_temp(lb_snow), .true.)
                     mt_wice_soisno(imt, lb_snow, ipatch) = mt_wice_soisno(imt, lb_snow, ipatch) - mt_flux
                     CALL material_book_evap_loss(imt, ipatch, mt_flux, subl_water, &
                        MT_EVAP_KIND_SUBL)
               ELSE
                  ! Ice insufficient (or empty): drain all ice, remainder from liquid
                     ! Part 1: drain ice completely
                     CALL exhaust_surface_phase(imt, ipatch, &
                        mt_wice_soisno(imt, lb_snow, ipatch), &
                        min(subl_water, max(water_ice_pool, 0._r8)), MT_EVAP_KIND_SUBL)
                  ! Part 2: excess from liquid at LIQUID ratio
                  deficit_water = subl_water - max(water_ice_pool, 0._r8)
                  IF (deficit_water > mt_tiny .and. wliq_soisno_bef(lb_snow) > mt_tiny) THEN
                        mt_flux = atmospheric_loss_material(mt_wliq_soisno(imt, lb_snow, ipatch), &
                           wliq_soisno_bef(lb_snow), deficit_water, layer_temp(lb_snow), .false.)
                        mt_wliq_soisno(imt, lb_snow, ipatch) = mt_wliq_soisno(imt, lb_snow, ipatch) - mt_flux
                        CALL material_book_evap_loss(imt, ipatch, mt_flux, deficit_water, &
                           MT_EVAP_KIND_SUBL)
                  ENDIF
               ENDIF
            ENDIF
            ! Update water_ice_pool to post-subl value for use in Step 2's ice deficit
            water_ice_pool = max(water_ice_pool - subl_water, 0._r8)

            ! Record the ACTUAL external wice change for Section 5b.
            ! In no-deficit case this equals (qfros - qsubl)*dt as before;
            ! in deficit case the excess sublimation was drawn from liquid,
            ! so the net wice change is only -(wice_bef + qfros*dt).
            d_wice_ext_snow = max(eff_qfros_snow, 0._r8) * deltim &
                            - min(subl_water, water_ice_pool_prefrost)

            ! Step 2: Liquid phase — rain + dew deposits, evaporation removes
            ! snowwater does: wliq(lb) += (pg_rain + qsdew - qseva)*dt
            ! If wliq goes negative, deficit removes from wice.
            !
            ! Add RAIN material here (before evaporation), mirroring snowwater's
            ! order: wliq += (pg_rain + qsdew - qseva)*dt applied together.
            ! This ensures evaporation can draw from rain that arrived this step.
            ! The rain material is then NOT re-added in section 0c percolation.
            IF (split_soilsnow) THEN
               mt_flux = mt_pg_rain_ground(imt, ipatch) * fsno
            ELSE
               mt_flux = mt_pg_rain_ground(imt, ipatch)
            ENDIF
            mt_wliq_soisno(imt, lb_snow, ipatch) = mt_wliq_soisno(imt, lb_snow, ipatch) + mt_flux

            ! Add dew
            mt_flux = max(eff_qsdew_snow, 0._r8) * &
               deposition_ratio_for(layer_temp(lb_snow), .false.) * deltim
            mt_wliq_soisno(imt, lb_snow, ipatch) = mt_wliq_soisno(imt, lb_snow, ipatch) + mt_flux
            a_mt_precip(imt, ipatch) = a_mt_precip(imt, ipatch) + mt_flux
            a_water_precip(imt, ipatch) = a_water_precip(imt, ipatch) + max(eff_qsdew_snow, 0._r8) * deltim

            ! Evaporate from the full pool (wliq_bef + qsubl_deficit + rain + dew)
            IF (eff_qseva_snow > mt_tiny) THEN
               ! Reconstruct pre-evap water pool (matches snowwater):
               ! wliq_bef + (pg_rain + qsdew)*dt + any ice deficit transferred earlier
               IF (split_soilsnow) THEN
                  water_liq_pool = wliq_soisno_bef(lb_snow) + (pg_rain*fsno + max(eff_qsdew_snow,0._r8)) * deltim
               ELSE
                  water_liq_pool = wliq_soisno_bef(lb_snow) + (pg_rain + max(eff_qsdew_snow,0._r8)) * deltim
               ENDIF
               ! Account for qsubl deficit transferred from ice (Step 1)
               ! Use PRE-subl ice pool to check if deficit actually occurred
               IF (water_ice_pool_prefrost < subl_water) THEN
                  water_liq_pool = water_liq_pool - (subl_water - water_ice_pool_prefrost)
               ENDIF
               water_liq_pool = max(water_liq_pool, 0._r8)
               evap_water = eff_qseva_snow * deltim

               IF (water_liq_pool - evap_water > mt_water_min_for_ratio) THEN
                  ! Normal case: enough liquid
                     mt_flux = atmospheric_loss_material(mt_wliq_soisno(imt, lb_snow, ipatch), &
                        water_liq_pool, evap_water, layer_temp(lb_snow), .false.)
                     mt_wliq_soisno(imt, lb_snow, ipatch) = mt_wliq_soisno(imt, lb_snow, ipatch) - mt_flux
                     CALL material_book_evap_loss(imt, ipatch, mt_flux, evap_water, &
                        MT_EVAP_KIND_SOILEVAP)
               ELSE
                  ! Liquid insufficient (or empty): drain liquid, remainder from ice
                     ! Part 1: drain liquid completely
                     CALL exhaust_surface_phase(imt, ipatch, &
                        mt_wliq_soisno(imt, lb_snow, ipatch), &
                        min(evap_water, max(water_liq_pool, 0._r8)), MT_EVAP_KIND_SOILEVAP)
                     ! Part 2: excess from ice at post-frost-post-subl ICE ratio
                     deficit_water = evap_water - max(water_liq_pool, 0._r8)
                     IF (deficit_water > mt_tiny .and. water_ice_pool > mt_tiny) THEN
                        deficit_water = min(deficit_water, water_ice_pool)
                        d_wice_ext_snow = d_wice_ext_snow - deficit_water
                        ! water_ice_pool is already post-frost, post-subl (updated at end of Step 1)
                           mt_flux = atmospheric_loss_material(mt_wice_soisno(imt, lb_snow, ipatch), &
                              water_ice_pool, deficit_water, layer_temp(lb_snow), .true.)
                           mt_wice_soisno(imt, lb_snow, ipatch) = mt_wice_soisno(imt, lb_snow, ipatch) - mt_flux
                           CALL material_book_evap_loss(imt, ipatch, mt_flux, deficit_water, &
                              MT_EVAP_KIND_SOILEVAP)
                        water_ice_pool = water_ice_pool - deficit_water
                     ENDIF
                  ENDIF
               ENDIF

         ENDIF  ! snl < 0

         ! ============================================================
         ! 0c. Snow percolation material tracking
         !
         !     Reconstruct the explicit snowwater loop:
         !       1. Apply top-layer external liquid/ice fluxes
         !       2. March qin/qout from top snow layer to bottom
         !       3. Bottom outflow becomes gwat material into surface pool
         !
         !     Use the separated rain material from mt_precip so only
         !     liquid rain enters snowwater; snowfall material stays in snow.
         ! ============================================================
         mt_gwat_snow = 0._r8
         IF (snl < 0) THEN
            IF (split_soilsnow) THEN
               snow_rain_input  = max(pg_rain * fsno, 0._r8) * deltim
               snow_dew_input   = max(qsdew_snow, 0._r8) * deltim
               snow_evap_output = max(qseva_snow, 0._r8) * deltim
               snow_frost_input = max(qfros_snow, 0._r8) * deltim
               snow_subl_output = max(qsubl_snow, 0._r8) * deltim
            ELSE
               snow_rain_input  = max(pg_rain, 0._r8) * deltim
               snow_dew_input   = max(qsdew_in, 0._r8) * deltim
               snow_evap_output = max(qseva_in, 0._r8) * deltim
               snow_frost_input = max(qfros_in, 0._r8) * deltim
               snow_subl_output = max(qsubl_in, 0._r8) * deltim
            ENDIF
            ! Rain material was already added to mt_wliq in Step 2 (section 0b).
            ! Top snow layer starts percolation with mt_qin_snow = 0.
            mt_qin_snow = 0._r8

            snow_wliq_before_flow = wliq_soisno_bef(lb_snow)
            IF (wice_soisno_bef(lb_snow) + snow_frost_input - snow_subl_output < 0._r8) THEN
               snow_wliq_before_flow = snow_wliq_before_flow + &
                  (wice_soisno_bef(lb_snow) + snow_frost_input - snow_subl_output)
            ENDIF
            snow_wliq_before_flow = snow_wliq_before_flow + snow_rain_input + snow_dew_input - snow_evap_output
            IF (snow_wliq_before_flow < 0._r8) THEN
               snow_wliq_before_flow = 0._r8
            ENDIF

            qin_snow = 0._r8
            DO j = lb_snow, 0
               IF (j == lb_snow) THEN
                  water_before_flow = max(snow_wliq_before_flow, 0._r8)
                  mt_before_flow = max(mt_wliq_soisno(imt, j, ipatch), 0._r8) + mt_qin_snow
               ELSE
                  water_before_flow = max(wliq_soisno_bef(j), 0._r8) + qin_snow
                  mt_before_flow = max(mt_wliq_soisno(imt, j, ipatch), 0._r8) + mt_qin_snow
               ENDIF

               IF (present(snow_qout_layer)) THEN
                  qout_snow = max(snow_qout_layer(j), 0._r8)
               ELSE
                  qout_snow = max(water_before_flow - max(wliq_soisno(j), 0._r8), 0._r8)
               ENDIF
               qout_snow = min(qout_snow, water_before_flow)

               ! --- meltwater <-> layer ice isotopic exchange ---
               ! Percolating meltwater trades isotopes with the ice it passes
               ! through, relaxing toward ice-water equilibrium
               ! R_liq = R_ice / alpha_ice_liq.  Because alpha_ice_liq > 1 the
               ! water leaves LIGHTER than the ice, which is the observed
               ! early-melt depletion and matching residual-pack enrichment
               ! (Taylor et al. 2001).  Internal ice<->water exchange within
               ! the layer: equal and opposite, so no accumulator and no
               ! numerical residual.
               IF (DEF_ISO_SNOWMELT_EQUILIBRATION > 0._r8 .and. &
                   material_fractionation_active(imt) .and. &
                   .not. material_is_nonvolatile_solute(imt) .and. &
                   water_before_flow > mt_water_min_for_ratio .and. &
                   wice_soisno(j) > mt_water_min_for_ratio) THEN
                  melt_exchange = material_equilibration_exchange(mt_before_flow, &
                     water_before_flow, &
                     max(mt_wice_soisno(imt, j, ipatch), 0._r8) / wice_soisno(j), &
                     1._r8 / max(material_alpha_ice_liq(imt, layer_temp(j)), mt_tiny), &
                     DEF_ISO_SNOWMELT_EQUILIBRATION)
                  IF (melt_exchange > 0._r8) THEN
                     melt_exchange = min(melt_exchange, &
                        max(mt_wice_soisno(imt, j, ipatch), 0._r8))
                  ELSE
                     melt_exchange = -min(-melt_exchange, max(mt_before_flow, 0._r8))
                  ENDIF
                  mt_before_flow = max(mt_before_flow + melt_exchange, 0._r8)
                  mt_wice_soisno(imt, j, ipatch) = &
                     max(mt_wice_soisno(imt, j, ipatch) - melt_exchange, 0._r8)
               ENDIF

               IF (qout_snow > mt_tiny .and. water_before_flow > mt_water_min_for_ratio .and. &
                   mt_before_flow > mt_tiny) THEN
                  ratio_src = mt_before_flow / water_before_flow
                  mt_qout_snow = qout_snow * ratio_src
                  mt_qout_snow = min(mt_qout_snow, mt_before_flow)
               ELSE
                  mt_qout_snow = 0._r8
               ENDIF

               mt_wliq_soisno(imt, j, ipatch) = max(mt_before_flow - mt_qout_snow, 0._r8)

               qin_snow = qout_snow
               mt_qin_snow = mt_qout_snow
            ENDDO

            mt_gwat_snow = mt_qin_snow
         ENDIF

         ! ============================================================
         ! 1. Surface water: mixed-pool approach
         !
         !    Negative qinfl has two physically different meanings:
         !      (a) qseva exceeds available surface input, so WATER extracts
         !          the deficit directly from layer 1 as evaporation output.
         !      (b) soil water truly exfiltrates into the surface pool.
         !    Only (b) should enter the surface mixed pool. Case (a) must be
         !    booked to a_mt_evap, otherwise active fractionation shows a
         !    storage loss equal to the negative-qinfl material diagnostic.
         ! ============================================================

         mt_soil_upflow = 0._r8
         mt_upflow_deferred = 0._r8
         mt_soil_evap = 0._r8
         top_exfil_water = 0._r8
         top_soil_evap_water = 0._r8
         imperv_wdsrf_loss = 0._r8
         imperv_soil_loss = 0._r8
         imperv_subl_loss = 0._r8
         IF (present(imperv_evap_wdsrf)) imperv_wdsrf_loss = max(imperv_evap_wdsrf, 0._r8)
         IF (present(imperv_evap_soil))  imperv_soil_loss  = max(imperv_evap_soil,  0._r8)
         IF (present(imperv_subl_soil))  imperv_subl_loss  = max(imperv_subl_soil,  0._r8)

         ! ---- Determine effective fluxes applied to gwat / soil layer 1 ----
         !
         ! Three cases for how WATER handles evaporation:
         !
         ! Case A: no snow (lb>=1)
         !   WATER uses the unsuffixed fluxes without split, but *_soil with
         !   split; THERMAL leaves unsuffixed fluxes zero in the split branch.
         !
         ! Case B: snow + split_soilsnow
         !   snowwater handles snow portion with qseva_snow/qsdew_snow/etc.
         !   gwat from snow bottom → then gwat += pg_rain*(1-fsno) - qseva_soil
         !   ice fluxes on soil layer 1: qfros_soil/qsubl_soil
         !   → surface pool deducts qseva_soil (soil portion only)
         !
         ! Case C: snow + NOT split_soilsnow
         !   snowwater handles ALL qseva/qsdew/qsubl/qfros on snow layers
         !   gwat is snowwater's OUTPUT (already post-evaporation)
         !   → surface pool should NOT deduct qseva again (already in gwat)
         !   → ice fluxes on soil layer 1: none separate (snowwater handled them)
         !     but qfros/qsubl still may affect soil layer 1 via wice update
         !
         IF (snl < 0 .and. .not. split_soilsnow) THEN
            ! Case C: evap/dew handled inside snowwater, not in surface pool
            eff_qseva = 0._r8
            eff_qsdew_topliq = 0._r8
            eff_qsubl_top = 0._r8
            eff_qfros_top = 0._r8
         ELSEIF (split_soilsnow) THEN
            ! Case B, or no snow with split: soil fluxes affect surface pool
            eff_qseva = qseva_soil
            eff_qsdew_topliq = qsdew_soil
            eff_qsubl_top = qsubl_soil
            eff_qfros_top = qfros_soil
         ELSE
            ! Case A: no snow, all fluxes apply
            eff_qseva = qseva_in
            eff_qsdew_topliq = qsdew_in
            eff_qsubl_top = qsubl_in
            eff_qfros_top = qfros_in
            ENDIF
            eff_qsdew_topliq = max(0._r8, eff_qsdew_topliq - dew_surface_water/max(deltim,mt_tiny))
            top_infil_water = max(qinfl, 0._r8) * deltim
            surface_base_water = max(wdsrf, 0._r8) + max(rsur, 0._r8)*deltim + top_infil_water
            surface_base_balance = surface_base_water
            IF (flood_water > 0._r8) THEN
               surface_base_balance = max(wdsrf, 0._r8) + max(rsur, 0._r8)*deltim &
                  + qinfl*deltim - flood_water
               surface_base_water = max(0._r8, surface_base_balance)
            ENDIF
            ! Dew overflow and ice-displaced old liquid arrive after the
            ! Richards solve. Neither may enter this step's infiltration.
            IF (late_surface_water > mt_tiny) THEN
               surface_base_water = max(surface_base_water-late_surface_water, 0._r8)
               surface_base_balance = surface_base_balance-late_surface_water
            ENDIF
            ! The pond overflow after the Richards solve (wdsrf > pondmx) is
            ! part of rsur; WATER_VSF reports no late runoff. Without flood
            ! water the overflow leaves the single mixed pool that also holds
            ! the exfiltrated water. With flood water, exfiltration and flood
            ! form the separate late pool below, so runoff beyond the ordinary
            ! pool can only be overflow of that late pool. Drawing it from the
            ! ordinary pool emptied the pool and the max(.,0) clamp created
            ! material (NEC: 13 mm exfiltration, 0.8 mm ordinary pool, 4 mm rsur).
            early_runoff_mt = early_runoff_water
            flood_overflow_water = 0._r8
            IF (flood_water > 0._r8 .and. early_runoff_water > surface_base_water) THEN
               early_runoff_mt = surface_base_water
               flood_overflow_water = early_runoff_water - surface_base_water
            ENDIF

         top_boundary_out_water = max(-qinfl, 0._r8) * deltim
         IF (top_boundary_out_water > mt_tiny) THEN
            ! qinfl = qgtop - d(wdsrf)/dt in WATER. Reconstruct qgtop from
            ! the final ponding change; for the evaporation-deficit cases
            ! that dominate active ptype=0 residuals, rsur is zero and this
            ! identifies negative qgtop exactly.
            IF (present(qgtop_solver)) THEN
               qgtop_est = qgtop_solver - flood_water/max(deltim,mt_tiny)
            ELSE
               qgtop_est = qinfl + (wdsrf - wdsrf_bef) / max(deltim, mt_tiny) &
                  + (late_runoff_water-flood_water-late_surface_water)/max(deltim,mt_tiny)
            ENDIF

            IF (eff_qseva > mt_tiny .and. qgtop_est < -mt_tiny) THEN
               ! Only the qseva deficit leaves layer 1 as
               ! evaporation. Negative qinfl can be far larger than qseva*deltim
               ! (a saturated column exfiltrating to the surface); booking all of
               ! it as evaporation removed up to ~24x the host's evaporated water
               ! from the material pool. The rest is ordinary exfiltration.
               top_soil_evap_water = min(top_boundary_out_water, eff_qseva * deltim)
               top_exfil_water = top_boundary_out_water - top_soil_evap_water
                  ! The host takes qseva from layer 1 before root uptake, so
                  ! the evaporation can exceed the liquid left in layer 1 here;
                  ! the Richards solve supplies the rest from below. Evaporate
                  ! what layer 1 holds and charge the rest at the layer's ratio
                  ! to the deferred upflow, paid after the qlayer transfers.
                  evap_pool_water = min(max(water_shadow(1), 0._r8), top_soil_evap_water)
                  kinetic_on_soil_surface = .true.
                  mt_soil_evap = atmospheric_loss_material(mt_wliq_soisno(imt, 1, ipatch), &
                     max(water_shadow(1), 0._r8), evap_pool_water, layer_temp(1), .false.)
                  evap_deficit_mt = 0._r8
                  IF (top_soil_evap_water > evap_pool_water .and. &
                      .not. material_is_nonvolatile_solute(imt)) THEN
                     evap_deficit_mt = (top_soil_evap_water - evap_pool_water) &
                        * evap_ratio_for(ratio_layer(1), layer_temp(1), .false.)
                  ENDIF
                  kinetic_on_soil_surface = .false.
                  mt_wliq_soisno(imt, 1, ipatch) = mt_wliq_soisno(imt, 1, ipatch) - mt_soil_evap
                  mt_upflow_deferred = mt_upflow_deferred + evap_deficit_mt
                  mt_soil_evap = mt_soil_evap + evap_deficit_mt
                  CALL material_book_evap_loss(imt, ipatch, mt_soil_evap, top_soil_evap_water, &
                     MT_EVAP_KIND_SOILEVAP)
               water_shadow(1) = water_shadow(1) - top_soil_evap_water
            ELSE
               top_exfil_water = top_boundary_out_water
            ENDIF
            IF (top_exfil_water > 0._r8) THEN
               mt_soil_upflow = top_exfil_water * ratio_layer(1)
               ! The host refills layer 1 from below within
               ! the same Richards solve, so exfiltration can exceed what layer 1
               ! holds before the qlayer transfers. Clamping here left the
               ! missing upflow material in the soil. Send the full upflow to the
               ! surface and let layer 1 pay the shortfall after section 2.
               ! Nonvolatile solutes keep the original clamp: they carry real
               ! gradients, and an unpaid deferral would create solute mass.
               IF (material_is_nonvolatile_solute(imt)) THEN
                  mt_soil_upflow = min(mt_soil_upflow, max(mt_wliq_soisno(imt, 1, ipatch), 0._r8))
                  mt_wliq_soisno(imt, 1, ipatch) = mt_wliq_soisno(imt, 1, ipatch) - mt_soil_upflow
               ELSEIF (mt_soil_upflow > max(mt_wliq_soisno(imt, 1, ipatch), 0._r8)) THEN
                  mt_upflow_deferred = mt_upflow_deferred &
                     + mt_soil_upflow - max(mt_wliq_soisno(imt, 1, ipatch), 0._r8)
                  mt_wliq_soisno(imt, 1, ipatch) = mt_wliq_soisno(imt, 1, ipatch) &
                     - max(mt_wliq_soisno(imt, 1, ipatch), 0._r8)
               ELSE
                  mt_wliq_soisno(imt, 1, ipatch) = mt_wliq_soisno(imt, 1, ipatch) - mt_soil_upflow
               ENDIF
               water_shadow(1) = water_shadow(1) - top_exfil_water
            ENDIF
         ENDIF

            ! Build pre-evaporation surface material pool.
         IF (snl < 0) THEN
            mt_pool_total = mt_wdsrf(imt, ipatch) + mt_gwat_snow
            IF (split_soilsnow) THEN
               mt_pool_total = mt_pool_total + mt_pg_rain_ground(imt, ipatch) * (1._r8 - fsno)
            ENDIF
         ELSE
            ! Only RAIN throughfall enters the surface mixed pool; snowfall
            ! throughfall is captured in mt_scv by mt_newsnow (mirroring
            ! the water side where `gwat = pg_rain + sm - qseva` never
            ! contains pg_snow — see MOD_SoilSnowHydrology.F90:789,809).
            ! Using mt_pg_to_ground here would double-count pg_snow material
            ! (once in mt_scv, once in mt_pool_total), breaking
            ! conservation on snl==0 patches with simultaneous rain+snow.
            mt_pool_total = mt_wdsrf(imt, ipatch) + mt_pg_rain_ground(imt, ipatch)
         ENDIF

         ! WATER uses gwat = pg_rain + sm - qseva when snl>=0, so snowmelt
         ! participates in the same pre-evaporation surface pool as rain.
         IF (snl >= 0 .and. sm > mt_tiny) THEN
            mt_pool_total = mt_pool_total + mt_sm_carry(imt, ipatch)
         ENDIF

         ! ---- Remove evaporation material from mixed pool ----
         ! gwat -= eff_qseva, so mt_pool must also be reduced.
         ! Use the ACTUAL pool ratio: mt_pool / (mt_pool / pool_ratio)
         ! where pool_ratio = mt_pool / water_pool_before_evap.
         ! Simplification: evap material = evap_water * (mt_pool / water_pool_before_evap)
         ! water_pool_before_evap ≈ mt_pool / pool_ratio. We don't know the exact
         ! water pool size at this point, but we can use:
         !   pool_ratio = mt_pool / water_pool_before_evap
         !   mt_evap = surface_evap_water * pool_ratio
         !            = surface_evap_water * mt_pool / water_bef
         ! water_bef = water_pool_after_evap + surface_evap_water
         ! If qseva exceeds the surface input/ponding supply, WATER_VSF
         ! draws the deficit directly from layer 1 as negative qgtop/qinfl;
         ! that deficit is already booked above as top_soil_evap_water.
         ! Only the remaining qseva can carry the surface mixed-pool
         ! signature. Using the full qseva here double-counts water and
         ! dilutes the surface material flux delta.
         gwat_evap = max(max(eff_qseva, 0._r8) * deltim - top_soil_evap_water &
                       - imperv_wdsrf_loss - imperv_soil_loss - imperv_subl_loss, 0._r8)
         ! The host adds qinfl_fld to qgtop after ordinary qseva. When the
         ! ordinary post-evaporation surface balance is negative, some newly
         ! arriving flood water satisfies that qseva before net infiltration.
         ! Remove its isotope/solute from the flood boundary input, not from
         ! soil or the earlier rain pool. E.g. qseva=1, flood=2, qinfl=1.
         flood_ground_evap_water = 0._r8
         IF (flood_water > 0._r8) THEN
            flood_ground_evap_water = min(flood_water, max(0._r8,-surface_base_balance), gwat_evap)
            ! If the solver leaves no ponding, runoff, or downward infiltration,
            ! the whole flood input was consumed at the evaporating boundary.
            ! Use that resolved destination rather than the cancellation-prone
            ! qseva-top_soil difference for a near-zero flood input.
            IF (flood_destination_water <= 0._r8) flood_ground_evap_water = flood_water
         ENDIF
         gwat_evap = max(gwat_evap - flood_ground_evap_water, 0._r8)
         flood_ground_evap_material = 0._r8
         IF (flood_ground_evap_water > mt_tiny .and. present(flood_mt_input)) THEN
            flood_ground_evap_material = atmospheric_loss_material( &
               flood_mt_input(imt), flood_water, flood_ground_evap_water, &
               layer_temp(1), .false.)
            CALL material_book_evap_loss(imt, ipatch, flood_ground_evap_material, &
               flood_ground_evap_water, MT_EVAP_KIND_SOILEVAP)
         ENDIF
         IF (imperv_soil_loss > mt_tiny) THEN
               kinetic_on_soil_surface = .true.
               mt_soil_evap = atmospheric_loss_material(mt_wliq_soisno(imt, 1, ipatch), &
                  max(water_shadow(1), 0._r8), imperv_soil_loss, layer_temp(1), .false.)
               kinetic_on_soil_surface = .false.
               mt_wliq_soisno(imt, 1, ipatch) = mt_wliq_soisno(imt, 1, ipatch) - mt_soil_evap
               CALL material_book_evap_loss(imt, ipatch, mt_soil_evap, imperv_soil_loss, &
                  MT_EVAP_KIND_SOILEVAP)
            water_shadow(1) = water_shadow(1) - imperv_soil_loss
         ENDIF
         IF (imperv_subl_loss > mt_tiny) THEN
               mt_soil_evap = atmospheric_loss_material(mt_wice_soisno(imt, 1, ipatch), &
                  max(wice_soisno_bef(1), 0._r8), imperv_subl_loss, layer_temp(1), .true.)
               mt_wice_soisno(imt, 1, ipatch) = mt_wice_soisno(imt, 1, ipatch) - mt_soil_evap
               CALL material_book_evap_loss(imt, ipatch, mt_soil_evap, imperv_subl_loss, &
                  MT_EVAP_KIND_SUBL)
         ENDIF
         ! Exfiltration joins the pool only after this
         ! evaporation (see below), but surface_base_water is built from the
         ! final wdsrf/rsur and already contains it. Leave it out of the
         ! evaporating pool too, or the evaporated material is diluted.
         evap_base_water = surface_base_water
         IF (flood_water <= 0._r8 .and. top_exfil_water > 0._r8) &
            evap_base_water = max(surface_base_water - top_exfil_water, 0._r8)
         IF (gwat_evap > mt_tiny .and. mt_pool_total > mt_tiny) THEN
            ! Reconstruct pre-evap water pool
            water_pool_total = evap_base_water + gwat_evap
               mt_gwat_evap = atmospheric_loss_material(mt_pool_total, water_pool_total, &
                  gwat_evap, layer_temp(1), .false.)
               mt_pool_total = mt_pool_total - mt_gwat_evap
               CALL material_book_evap_loss(imt, ipatch, mt_gwat_evap, gwat_evap, &
                  MT_EVAP_KIND_SOILEVAP)
         ENDIF
         IF (imperv_wdsrf_loss > mt_tiny .and. mt_pool_total > mt_tiny) THEN
            water_pool_total = evap_base_water + imperv_wdsrf_loss
               mt_gwat_evap = atmospheric_loss_material(mt_pool_total, water_pool_total, &
                  imperv_wdsrf_loss, layer_temp(1), .false.)
               mt_pool_total = mt_pool_total - mt_gwat_evap
               CALL material_book_evap_loss(imt, ipatch, mt_gwat_evap, imperv_wdsrf_loss, &
                  MT_EVAP_KIND_SOILEVAP)
         ENDIF

         ! Soil exfiltration is produced by the VSF solve after qseva has
         ! already been applied to gwat, so it should feed runoff/wdsrf but
         ! should not be available to this step's ground evaporation.
         IF (flood_water <= 0._r8) mt_pool_total = mt_pool_total + mt_soil_upflow

         ! Pore-excess dew and frost-displaced liquid arrive only after the
         ! soil solve; Section 4a's late pool handles both below.

         ! ---- Add drip/flood/paddy irrigation material ----
         ! Water side (MOD_SoilSnowHydrology.F90:828-832) adds the
         ! irrigation mass to gwat AFTER qseva, so none of it is
         ! evaporated at this step. Irrigation is an INTERNAL transfer
         ! from the per-patch waterstorage reservoir to the surface
         ! pool — water-side totwb/endwb include waterstorage
         ! (CoLMMAIN:806), so the atmospheric input term does not count
         ! irrigation. The material path mirrors this: debit
         ! mt_waterstorage at its current ratio (Phase-1 invariant
         ! keeps it at R_init, so we use R_precip as a proxy) and do
         ! NOT add to a_mt_precip. Without this mirroring the material
         ! budget claims more atmospheric input than the water budget.
         IF (qflx_irrig_ground > mt_tiny) THEN
            storage_ratio = R_precip
            IF (allocated(mt_waterstorage) .and. present(waterstorage_patch)) THEN
               CALL mt_equilibrate_dissolved(imt, max(waterstorage_patch, 0._r8), &
                  mt_waterstorage(imt, ipatch), mt_waterstorage_solid(imt, ipatch))
               IF (waterstorage_patch > mt_water_min_for_ratio) THEN
                  storage_ratio = mt_waterstorage(imt, ipatch) / waterstorage_patch
               ENDIF
            ENDIF
            mt_flux = qflx_irrig_ground * deltim * storage_ratio
            IF (allocated(mt_waterstorage)) THEN
               mt_flux = min(max(mt_flux, 0._r8), max(mt_waterstorage(imt, ipatch), 0._r8))
            ENDIF
            mt_pool_total = mt_pool_total + mt_flux
            IF (allocated(mt_waterstorage)) THEN
               mt_waterstorage(imt, ipatch) = mt_waterstorage(imt, ipatch) - mt_flux
            ENDIF
         ENDIF

         ! Use OUTPUT-side water pool for ratio computation.
         ! This guarantees: distributed = water_pool * ratio = mt_pool (exact).
         water_pool_total = surface_base_water
         mt_surface_collapse_residual = 0._r8

         ! Rewetting dissolves the numerical surface residue before runoff or
         ! infiltration is partitioned. It remains immobile while no resolved
         ! carrier water exists.
         IF (material_is_nonvolatile_solute(imt) .and. &
             (water_pool_total > mt_water_min_for_ratio .or. &
              mt_has_dissolved_limit(imt)) .and. &
             mt_surface_residue(imt, ipatch) > mt_tiny) THEN
            mt_pool_total = mt_pool_total + mt_surface_residue(imt, ipatch)
            mt_surface_residue(imt, ipatch) = 0._r8
         ENDIF
         CALL mt_equilibrate_dissolved(imt, water_pool_total, mt_pool_total, &
            mt_surface_solid(imt, ipatch))

         ! Compute mixed-pool ratio
         IF (water_pool_total > mt_water_min_for_ratio .and. mt_pool_total > mt_tiny) THEN
            ! The surface mixed-pool concentration is bounded at source: the
            ! evaporative-loss routine (material_evaporative_mt_loss) caps the
            ! residual pool R at the physical ceiling, so mt_pool_total can no
            ! longer run away near drydown. No export-side clamp is needed here.
            ratio = mt_pool_total / water_pool_total
         ELSEIF (mt_pool_total <= mt_tiny) THEN
            ratio = 0._r8
         ELSE
            ! No physically resolved water destination remains for this
            ! mixed-pool material.  mt_tiny is only a numerical nonzero test;
            ! the ratio denominator must use mt_water_min_for_ratio so a
            ! near-dry water sliver cannot export a non-physical
            ! material/water concentration.
            !
            ! Do not overwrite the snow/soil/throughfall signature with
            ! precipitation. A water isotope has no immobile phase of its own,
            ! so it must not be parked in a dry residue: when a resolved
            ! surface carrier water exists (flood infiltration water and/or
            ! the post-solver dew/frost late pool), the ordinary mixed-pool
            ! material stays in the pool and is mixed into that carrier below -
            ! the surface water body is one well-mixed pond, so the material
            ! follows the water that is actually resolved. Only a genuinely
            ! carrierless collapse (no flood water and no late surface water)
            ! is booked as a numerical residual, which the balance report
            ! exposes; that case is a real inconsistency, not this fix's
            ! target.
            IF (material_is_nonvolatile_solute(imt)) THEN
               mt_surface_residue(imt, ipatch) = mt_surface_residue(imt, ipatch) + &
                  mt_pool_total
               mt_pool_total = 0._r8
            ELSEIF (flood_water > 0._r8 .or. late_surface_water > mt_tiny) THEN
               ! Keep mt_pool_total for the late-pool mixing step below.
               ! ratio stays zero, so no ordinary runoff/infiltration export
               ! claims material without its carrier water.
            ELSE
               mt_surface_collapse_residual = mt_pool_total
               IF (allocated(mt_numerical_residual_step)) THEN
                  mt_numerical_residual_step(imt, ipatch) = &
                     mt_numerical_residual_step(imt, ipatch) - mt_surface_collapse_residual
               ENDIF
               mt_pool_total = 0._r8
            ENDIF
            ratio = 0._r8
         ENDIF

         ! Distribute: wdsrf residual, surface runoff, infiltration
         IF (early_runoff_mt > mt_tiny) THEN
            mt_flux = early_runoff_mt * ratio
            a_mt_rsur(imt, ipatch) = a_mt_rsur(imt, ipatch) + mt_flux
            a_mt_rnof(imt, ipatch) = a_mt_rnof(imt, ipatch) + mt_flux
            mt_rnof_step(imt, ipatch) = mt_rnof_step(imt, ipatch) + mt_flux
            mt_pool_total = max(mt_pool_total-mt_flux,0._r8)
         ENDIF

         late_ratio = ratio
         IF (flood_water > 0._r8) THEN
            ! Flood water joins qgtop only after ordinary evaporation and
            ! runoff. VSF may pond part of it, so mix at that late boundary
            ! and partition by the solver's actual wdsrf/qinfl outputs.
            late_material = mt_pool_total + mt_soil_upflow
            late_scale = abs(mt_pool_total) + abs(mt_soil_upflow)
            IF (present(flood_mt_input)) THEN
               late_material = late_material + flood_mt_input(imt) - flood_ground_evap_material
               late_scale = late_scale + abs(flood_mt_input(imt)) + abs(flood_ground_evap_material)
            ENDIF
            ! Directly reconstruct the carrier after ordinary runoff and
            ! before post-solver frost/dew.  The pre-flood ordinary balance
            ! may be negative when flood water first pays soil evaporation.
            late_water = flood_destination_water + flood_overflow_water
            ! A positive flood-water microcarrier can still have a bounded
            ! isotope ratio below the generic dry-pool denominator floor.
            ! Partition it conservatively; keep genuinely dry or excessive
            ! isotope mass on the fatal path.
            IF (late_water > mt_water_min_for_ratio .or. &
                (material_is_isotope(imt) .and. late_water > 0._r8 .and. late_material > mt_tiny .and. &
                 late_material <= materials(imt)%ref_ratio * &
                    (1._r8 + mt_delta_sanity_max/1000._r8) * late_water)) THEN
               ! A dry surface residue may have survived the earlier
               ! ordinary-water pool. Flood-only wetting dissolves it here,
               ! after ordinary runoff, before pond/infiltration partition.
               IF (material_is_nonvolatile_solute(imt) .and. &
                   mt_surface_residue(imt,ipatch) > mt_tiny) THEN
                  late_material = late_material + mt_surface_residue(imt,ipatch)
                  mt_surface_residue(imt,ipatch) = 0._r8
               ENDIF
               ! N03: the flood/late water is the first mobile carrier of a
               ! possibly dry surface pool.  Re-run the same-position
               ! liquid/solid equilibrium BEFORE computing late_ratio and
               ! partitioning infiltration/ponding/overflow: a dry surface
               ! solid must dissolve up to its water-based solubility, and a
               ! saturated pool must leave the excess in mt_surface_solid
               ! instead of carrying an unphysical concentration into the
               ! soil.  Ordinary-pool equilibrium above cannot see this late
               ! water because it entered after ordinary evaporation/runoff.
               IF (material_is_nonvolatile_solute(imt) .and. mt_has_dissolved_limit(imt)) THEN
                  CALL mt_equilibrate_dissolved(imt, late_water, late_material, &
                     mt_surface_solid(imt,ipatch))
               ENDIF
               late_ratio = max(late_material,0._r8)/late_water
            ELSE
               IF (material_is_nonvolatile_solute(imt)) THEN
                  mt_surface_residue(imt,ipatch) = mt_surface_residue(imt,ipatch) + max(late_material,0._r8)
               ELSEIF (late_material > max(mt_tiny, 1.e-10_r8*late_scale, mt_water_min_for_ratio * &
                       materials(imt)%ref_ratio*(1._r8 + mt_delta_sanity_max/1000._r8))) THEN
                  CALL CoLM_stop('grid flood material: unresolved dry isotope surface pool')
               ELSEIF (late_material > 0._r8 .and. allocated(mt_numerical_residual_step)) THEN
                  mt_numerical_residual_step(imt, ipatch) = &
                     mt_numerical_residual_step(imt, ipatch) - late_material
               ENDIF
               late_ratio = 0._r8
               late_material = 0._r8
            ENDIF
            mt_pool_total = late_material
         ENDIF
         mt_wdsrf(imt, ipatch) = max(wdsrf, 0._r8) * late_ratio
         IF (flood_overflow_water > mt_tiny) THEN
            mt_flux = flood_overflow_water * late_ratio
            a_mt_rsur(imt, ipatch) = a_mt_rsur(imt, ipatch) + mt_flux
            a_mt_rnof(imt, ipatch) = a_mt_rnof(imt, ipatch) + mt_flux
            mt_rnof_step(imt, ipatch) = mt_rnof_step(imt, ipatch) + mt_flux
         ENDIF

         pending_surface_material = 0._r8
         IF (late_surface_water > mt_tiny) THEN
            pending_surface_material = max(mt_pool_total-max(qinfl,0._r8)*late_ratio*deltim &
               - flood_overflow_water*late_ratio, 0._r8)
         ENDIF

         IF (qinfl > mt_tiny) THEN
            a_mt_qinfl(imt, ipatch) = a_mt_qinfl(imt, ipatch) + qinfl * late_ratio * deltim
         ELSEIF (qinfl < -mt_tiny) THEN
            a_mt_qinfl(imt, ipatch) = a_mt_qinfl(imt, ipatch) + qinfl * ratio_layer(1) * deltim
         ENDIF
         IF (present(flood_mt_input)) THEN
            a_mt_precip(imt,ipatch) = a_mt_precip(imt,ipatch) + flood_mt_input(imt)
         ENDIF

         ! ============================================================
         ! 2. Soil layers: qlayer-driven flux tracking
         !
         !    qlayer(j) > 0: downward from j to j+1, source = layer j
         !    qlayer(j) < 0: upward from j+1 to j, source = layer j+1
         !
         !    Each flux removes from source, adds to destination.
         ! ============================================================

         ! --- surface → layer 1 infiltration (downward only) ---
         ! WATER_VSF reports this as qinfl = qgtop - d(wdsrf)/dt. qlayer(0)
         ! can be zero when the top-layer storage changes through the VSF
         ! wetting-front representation, so using qlayer(0) here misses the
         ! actual water entering layer 1.
            IF (qinfl > mt_tiny) THEN
               mt_flux = qinfl * late_ratio * deltim
               mt_wliq_soisno(imt, 1, ipatch) = mt_wliq_soisno(imt, 1, ipatch) + mt_flux
               water_shadow(1) = water_shadow(1) + qinfl * deltim
            ENDIF

         ! --- qlayer(j): layer j ↔ layer j+1 ---
         DO j = 1, nl_soil
            layer_transport_ratio(j) = current_liq_ratio(j)
         ENDDO
         ! The faces are applied top-down, but a layer
         ! that donates upward is refilled from below in the same solve.
         ! Clamping each flux to the donor's current stock left material behind
         ! whenever water moved through several layers in one step. First try
         ! the faces unclamped with the snapshot ratios (identical arithmetic
         ! whenever no clamp binds); fall back to the clamped sweep only if a
         ! layer would end negative by more than rounding of what passed it.
         ! A layer that passes on more water than it holds at the snapshot is
         ! refilled by its supplier within the step; its stock alone cannot
         ! give the outgoing ratio (with a near-empty stock it is just noise).
         ! Mix stock and inflow instead, following the flow direction. The old
         ! clamped sweep always truncated such a flux, so no other step changes.
         IF (.not. material_is_nonvolatile_solute(imt)) THEN
            DO j = nl_soil - 1, 2, -1
               IF (qlayer(j) < -mt_tiny .and. qlayer(j-1) < -mt_tiny) THEN
                  IF (abs(qlayer(j-1)) * deltim > water_shadow(j)) THEN
                     layer_transport_ratio(j) = (max(mt_wliq_soisno(imt, j, ipatch), 0._r8) &
                        + abs(qlayer(j)) * deltim * layer_transport_ratio(j+1)) &
                        / (max(water_shadow(j), 0._r8) + abs(qlayer(j)) * deltim)
                  ENDIF
               ENDIF
            ENDDO
            DO j = 2, nl_soil - 1
               IF (qlayer(j-1) > mt_tiny .and. qlayer(j) > mt_tiny) THEN
                  IF (qlayer(j) * deltim > water_shadow(j)) THEN
                     layer_transport_ratio(j) = (max(mt_wliq_soisno(imt, j, ipatch), 0._r8) &
                        + qlayer(j-1) * deltim * layer_transport_ratio(j-1)) &
                        / (max(water_shadow(j), 0._r8) + qlayer(j-1) * deltim)
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
         mt_wliq_qlayer_bef(1:nl_soil) = mt_wliq_soisno(imt, 1:nl_soil, ipatch)
         water_shadow_qlayer_bef(1:nl_soil) = water_shadow(1:nl_soil)
         ! Nonvolatile solutes go straight to the original clamped sweep.
         DO qlayer_pass = merge(2, 1, material_is_nonvolatile_solute(imt)), 2
            qlayer_clamp = qlayer_pass == 2
            IF (qlayer_clamp) THEN
               mt_wliq_soisno(imt, 1:nl_soil, ipatch) = mt_wliq_qlayer_bef(1:nl_soil)
               water_shadow(1:nl_soil) = water_shadow_qlayer_bef(1:nl_soil)
            ENDIF
            qlayer_through(1:nl_soil) = abs(mt_wliq_qlayer_bef(1:nl_soil))
         DO j = 1, nl_soil - 1
            IF (qlayer(j) > mt_tiny) THEN
               ! Downward: j → j+1
               ratio_src = layer_transport_ratio(j)
               mt_flux = qlayer(j) * ratio_src * deltim
                  IF (qlayer_clamp) mt_flux = min(mt_flux, max(mt_wliq_soisno(imt, j, ipatch), 0._r8))
                  qlayer_through(j) = qlayer_through(j) + abs(mt_flux)
                  qlayer_through(j+1) = qlayer_through(j+1) + abs(mt_flux)
                  mt_wliq_soisno(imt, j,   ipatch) = mt_wliq_soisno(imt, j,   ipatch) - mt_flux
                  mt_wliq_soisno(imt, j+1, ipatch) = mt_wliq_soisno(imt, j+1, ipatch) + mt_flux
                  water_shadow(j)   = water_shadow(j)   - qlayer(j) * deltim
                  water_shadow(j+1) = water_shadow(j+1) + qlayer(j) * deltim
               ELSEIF (qlayer(j) < -mt_tiny) THEN
                  ! Upward: j+1 → j
                  ratio_src = layer_transport_ratio(j+1)
                  mt_flux = abs(qlayer(j)) * ratio_src * deltim
                  IF (qlayer_clamp) mt_flux = min(mt_flux, max(mt_wliq_soisno(imt, j+1, ipatch), 0._r8))
                  qlayer_through(j) = qlayer_through(j) + abs(mt_flux)
                  qlayer_through(j+1) = qlayer_through(j+1) + abs(mt_flux)
                  mt_wliq_soisno(imt, j+1, ipatch) = mt_wliq_soisno(imt, j+1, ipatch) - mt_flux
                  mt_wliq_soisno(imt, j,   ipatch) = mt_wliq_soisno(imt, j,   ipatch) + mt_flux
                  water_shadow(j+1) = water_shadow(j+1) + qlayer(j) * deltim
                  water_shadow(j)   = water_shadow(j)   - qlayer(j) * deltim
               ENDIF
            ENDDO
            IF (qlayer_clamp) EXIT
            sweep_ok(1:nl_soil) = mt_wliq_soisno(imt, 1:nl_soil, ipatch) >= -1.e-12_r8 * qlayer_through(1:nl_soil) &
                    .or. mt_wliq_soisno(imt, 1:nl_soil, ipatch) >= mt_wliq_qlayer_bef(1:nl_soil)
            ! Isotopes: a deficit that only mirrors the host's own rounding of
            ! the layer water is rounding too (see MOD_MT_Isotope_SoilFlow).
            IF (material_is_isotope(imt)) THEN
               DO j = 1, nl_soil
                  IF (.not. sweep_ok(j)) sweep_ok(j) = isotope_sweep_deficit_from_host_water( &
                     mt_wliq_soisno(imt, j, ipatch), water_shadow(j), layer_transport_ratio(j))
               ENDDO
            ENDIF
            IF (all(sweep_ok(1:nl_soil))) THEN
               ! Rounding dust left by the unclamped sweep: zero it as an
               ! explicit numerical source so no negative stock survives.
               DO j = 1, nl_soil
                  IF (mt_wliq_soisno(imt, j, ipatch) < 0._r8 .and. &
                      mt_wliq_qlayer_bef(j) >= 0._r8) THEN
                     IF (allocated(mt_numerical_residual_step)) &
                        mt_numerical_residual_step(imt, ipatch) = &
                           mt_numerical_residual_step(imt, ipatch) - mt_wliq_soisno(imt, j, ipatch)
                     mt_wliq_soisno(imt, j, ipatch) = 0._r8
                  ENDIF
               ENDDO
               EXIT
            ENDIF
         ENDDO

         ! Layer 1 now holds the water that reached it
         ! from below; it pays the deferred exfiltration material. What it
         ! cannot cover yet stays deferred: in a frozen top layer the host
         ! takes the evaporation from the ice (wblc_ice_sink), and that
         ! water reaches the liquid only in the phase-change section below,
         ! where the rest is paid (see there).
         IF (mt_upflow_deferred > 0._r8) THEN
            mt_flux = min(mt_upflow_deferred, max(mt_wliq_soisno(imt, 1, ipatch), 0._r8))
            mt_wliq_soisno(imt, 1, ipatch) = mt_wliq_soisno(imt, 1, ipatch) - mt_flux
            mt_upflow_deferred = mt_upflow_deferred - mt_flux
         ENDIF

         ! --- liquid-phase molecular diffusion between soil layers ---
         ! Internal, equal-and-opposite exchange: it moves no water and no
         ! material across the patch boundary, so it needs no flux accumulator
         ! and the balance check closes unchanged.  Runs on the post-advection
         ! state (water_shadow) so advection and diffusion see the same water.
         ! Compute all faces from a timestep snapshot, then scale each source
         ! layer's combined outgoing flux before applying: a layer can donate
         ! to both neighbours in the same step without going negative.
         IF ((DEF_ISO_SOIL_DIFFUSION .or. DEF_ISO_SOIL_VAPOR_DIFFUSION) .and. present(dz_soi_frac) &
             .and. present(porsl_frac)) THEN
            diff_face = 0._r8
            diff_out = 0._r8
            diff_scale = 1._r8
            DO j = 1, nl_soil - 1
               IF (water_shadow(j) <= mt_water_min_for_ratio) CYCLE
               IF (water_shadow(j+1) <= mt_water_min_for_ratio) CYCLE
               diff_liquid = material_liquid_self_diffusivity(imt, layer_temp(j))
               IF (diff_liquid <= 0._r8 .and. .not. DEF_ISO_SOIL_VAPOR_DIFFUSION) EXIT
               ! Liquid film + pore vapour.  Their moisture dependences are
               ! opposite (vapour uses AIR-filled porosity), so the sum keeps a
               ! working pathway at every wetness: liquid dominates wet soil,
               ! vapour dominates once the evaporation front retreats into the
               ! pores.  Vapour is folded in as an equivalent liquid
               ! diffusivity, so the flux form and overshoot limiter below are
               ! unchanged.  Needs forc_psrf, so it is gated on it.
               d_eff_up = 0._r8
               d_eff_dn = 0._r8
               IF (DEF_ISO_SOIL_DIFFUSION .and. diff_liquid > 0._r8) THEN
                  d_eff_up = material_soil_effective_diffusivity(water_shadow(j), &
                     dz_soi_frac(j), porsl_frac(j), diff_liquid)
                  d_eff_dn = material_soil_effective_diffusivity(water_shadow(j+1), &
                     dz_soi_frac(j+1), porsl_frac(j+1), &
                     material_liquid_self_diffusivity(imt, layer_temp(j+1)))
               ENDIF
               IF (DEF_ISO_SOIL_VAPOR_DIFFUSION .and. material_is_isotope(imt) &
                   .and. present(forc_psrf_frac)) THEN
                  ! Vapour transport is an isotope mechanism (pore-vapour /
                  ! liquid equilibrium + light/heavy air diffusivity).  The
                  ! family-gated helper also rejects nonvolatile solutes, which
                  ! have no vapour carrier; they must not inherit this path
                  ! through alpha=1/diff_ratio=1.  Pass ice as well: it blocks
                  ! pore air just as liquid does.
                  d_eff_up = d_eff_up + material_isotope_soil_vapor_equivalent_diffusivity( &
                     imt, water_shadow(j), max(wice_soisno(j), 0._r8), &
                     dz_soi_frac(j), porsl_frac(j), layer_temp(j), forc_psrf_frac)
                  d_eff_dn = d_eff_dn + material_isotope_soil_vapor_equivalent_diffusivity( &
                     imt, water_shadow(j+1), max(wice_soisno(j+1), 0._r8), &
                     dz_soi_frac(j+1), porsl_frac(j+1), &
                     layer_temp(j+1), forc_psrf_frac)
               ENDIF
               IF (d_eff_up <= 0._r8 .or. d_eff_dn <= 0._r8) CYCLE
               ratio_up = mt_wliq_soisno(imt, j,   ipatch) / water_shadow(j)
               ratio_dn = mt_wliq_soisno(imt, j+1, ipatch) / water_shadow(j+1)
               diff_face(j) = material_soil_diffusive_transfer(ratio_up, ratio_dn, &
                  water_shadow(j), water_shadow(j+1), d_eff_up, d_eff_dn, &
                  dz_soi_frac(j), dz_soi_frac(j+1), deltim)
               IF (diff_face(j) > 0._r8) THEN
                  diff_out(j) = diff_out(j) + diff_face(j)
               ELSE
                  diff_out(j+1) = diff_out(j+1) - diff_face(j)
               ENDIF
            ENDDO
            DO j = 1, nl_soil
               IF (diff_out(j) > max(mt_wliq_soisno(imt, j, ipatch), 0._r8) .and. &
                   diff_out(j) > mt_tiny) THEN
                  diff_scale(j) = max(mt_wliq_soisno(imt, j, ipatch), 0._r8) / diff_out(j)
               ENDIF
            ENDDO
            DO j = 1, nl_soil - 1
               IF (diff_face(j) > 0._r8) THEN
                  diff_transfer = diff_face(j) * diff_scale(j)
               ELSE
                  diff_transfer = diff_face(j) * diff_scale(j+1)
               ENDIF
               mt_wliq_soisno(imt, j,   ipatch) = &
                  mt_wliq_soisno(imt, j,   ipatch) - diff_transfer
               mt_wliq_soisno(imt, j+1, ipatch) = &
                  mt_wliq_soisno(imt, j+1, ipatch) + diff_transfer
            ENDDO
         ENDIF

         ! --- firn vapour diffusion between snow layers ---
         ! The carrier here is ICE, not the liquid film: isotopes barely move
         ! through the lattice, so a pack smooths its own vertical signal
         ! through pore vapour.  Same internal, equal-and-opposite exchange as
         ! the soil sweep, so no accumulator and no residual.  Implied
         ! timescale is years at seasonal layer thickness, so this matters for
         ! multi-year packs far more than for snow that melts out annually.
         IF (DEF_ISO_SOIL_VAPOR_DIFFUSION .and. material_is_isotope(imt) .and. snl < 0 &
             .and. present(dz_sno_frac) .and. present(forc_psrf_frac)) THEN
            diff_face = 0._r8
            diff_out = 0._r8
            diff_scale = 1._r8
            DO j = lb, -1
               IF (dz_sno_frac(j) <= 0._r8 .or. dz_sno_frac(j+1) <= 0._r8) CYCLE
               IF (wice_soisno(j) <= mt_water_min_for_ratio) CYCLE
               IF (wice_soisno(j+1) <= mt_water_min_for_ratio) CYCLE
               d_eff_up = material_isotope_snow_vapor_equivalent_diffusivity( &
                  imt, max(wliq_soisno(j), 0._r8), max(wice_soisno(j), 0._r8), &
                  dz_sno_frac(j), layer_temp(j), forc_psrf_frac)
               d_eff_dn = material_isotope_snow_vapor_equivalent_diffusivity( &
                  imt, max(wliq_soisno(j+1), 0._r8), max(wice_soisno(j+1), 0._r8), &
                  dz_sno_frac(j+1), layer_temp(j+1), forc_psrf_frac)
               IF (d_eff_up <= 0._r8 .or. d_eff_dn <= 0._r8) CYCLE
               ratio_up = mt_wice_soisno(imt, j,   ipatch) / wice_soisno(j)
               ratio_dn = mt_wice_soisno(imt, j+1, ipatch) / wice_soisno(j+1)
               diff_face(j) = material_soil_diffusive_transfer(ratio_up, ratio_dn, &
                  wice_soisno(j), wice_soisno(j+1), d_eff_up, d_eff_dn, &
                  dz_sno_frac(j), dz_sno_frac(j+1), deltim)
               IF (diff_face(j) > 0._r8) THEN
                  diff_out(j) = diff_out(j) + diff_face(j)
               ELSE
                  diff_out(j+1) = diff_out(j+1) - diff_face(j)
               ENDIF
            ENDDO
            DO j = lb, 0
               IF (diff_out(j) > max(mt_wice_soisno(imt, j, ipatch), 0._r8) .and. &
                   diff_out(j) > mt_tiny) THEN
                  diff_scale(j) = max(mt_wice_soisno(imt, j, ipatch), 0._r8) / diff_out(j)
               ENDIF
            ENDDO
            DO j = lb, -1
               IF (diff_face(j) > 0._r8) THEN
                  diff_transfer = diff_face(j) * diff_scale(j)
               ELSE
                  diff_transfer = diff_face(j) * diff_scale(j+1)
               ENDIF
               mt_wice_soisno(imt, j,   ipatch) = &
                  mt_wice_soisno(imt, j,   ipatch) - diff_transfer
               mt_wice_soisno(imt, j+1, ipatch) = &
                  mt_wice_soisno(imt, j+1, ipatch) + diff_transfer
            ENDDO
         ENDIF

         ! ============================================================
         ! 3. Groundwater: qcharge-driven
         ! ============================================================
         ! Root exchange and the resolved aquifer share of baseflow have
         ! already changed this donor.  wa_bef is not its current carrier;
         ! using it after a hydraulic return can even reverse the isotope
         ! flux into the positive bottom-soil pool.
         aquifer_water_pre_qcharge = wa_bef - etroot_aquifer
         IF (resolved_rsub) aquifer_water_pre_qcharge = aquifer_water_pre_qcharge &
            - max(rsub_source_aquifer, 0._r8)
         IF (qcharge_eff > mt_tiny) THEN
            j = nl_soil
            ratio_src = layer_transport_ratio(j)
            mt_flux = qcharge_eff * ratio_src * deltim
               mt_flux = min(mt_flux, max(mt_wliq_soisno(imt, j, ipatch), 0._r8))
               mt_wliq_soisno(imt, j, ipatch) = mt_wliq_soisno(imt, j, ipatch) - mt_flux
               mt_wa(imt, ipatch) = mt_wa(imt, ipatch) + mt_flux
               a_mt_qcharge(imt, ipatch) = a_mt_qcharge(imt, ipatch) + mt_flux
               water_shadow(j) = water_shadow(j) - qcharge_eff * deltim
            ELSEIF (qcharge_eff < -mt_tiny) THEN
            CALL check_isotope_aquifer(imt, ipatch, aquifer_water_pre_qcharge, &
               mt_wa(imt, ipatch), 'qcharge export')
            IF (abs(material_aquifer_actual_water(aquifer_water_pre_qcharge, &
                                                aquifer_ref_water)) > mt_water_min_for_ratio) THEN
               ratio_src = material_aquifer_isotope_ratio(aquifer_water_pre_qcharge, &
                  mt_wa(imt, ipatch), aquifer_ref_water, aquifer_ref_mass, source_fallback_ratio)
            ELSEIF (material_is_nonvolatile_solute(imt)) THEN
               ratio_src = 0._r8
            ELSE
               ratio_src = source_fallback_ratio
            ENDIF
            mt_flux = aquifer_export_material(abs(qcharge_eff) * deltim, abs(qcharge_eff) * ratio_src * deltim, &
               material_aquifer_actual_water(aquifer_water_pre_qcharge, aquifer_ref_water), &
               material_aquifer_actual_mass(mt_wa(imt, ipatch), aquifer_ref_mass))
            mt_wa(imt, ipatch) = mt_wa(imt, ipatch) - mt_flux
               j = nl_soil
               mt_wliq_soisno(imt, j, ipatch) = mt_wliq_soisno(imt, j, ipatch) + mt_flux
               a_mt_qcharge(imt, ipatch) = a_mt_qcharge(imt, ipatch) - mt_flux
               water_shadow(j) = water_shadow(j) - qcharge_eff * deltim
            ENDIF

         ! A previously dry aquifer can be recharged during this step.  Make
         ! its conservative residue mobile immediately after recharge, before
         ! any subsequent aquifer-state diagnostics or exports, rather than
         ! delaying the first flush until the next timestep.
         IF (material_is_nonvolatile_solute(imt) .and. &
             wa > mt_water_min_for_ratio .and. &
             mt_subsurface_residue(imt, ipatch) > mt_tiny) THEN
            mt_wa(imt, ipatch) = mt_wa(imt, ipatch) + &
               mt_subsurface_residue(imt, ipatch)
            mt_subsurface_residue(imt, ipatch) = 0._r8
         ENDIF
         CALL mt_equilibrate_dissolved(imt, wa, mt_wa(imt, ipatch), &
            mt_subsurface_solid(imt, ipatch))
         IF (material_is_isotope(imt)) THEN
            IF (.not. material_aquifer_isotope_state_valid(wa, mt_wa(imt, ipatch), &
                  materials(imt)%ref_ratio, aquifer_ref_water, aquifer_ref_mass)) THEN
               aquifer_orphan_mass = material_aquifer_actual_mass(mt_wa(imt, ipatch), aquifer_ref_mass)
               IF (material_aquifer_actual_water(wa, aquifer_ref_water) < -mt_water_min_for_ratio) &
                  aquifer_orphan_mass = aquifer_orphan_mass - &
                     material_aquifer_actual_water(wa, aquifer_ref_water) * aquifer_ratio
               mt_wa(imt, ipatch) = mt_wa(imt, ipatch) - aquifer_orphan_mass
               IF (allocated(mt_numerical_residual_step)) &
                  mt_numerical_residual_step(imt, ipatch) = &
                     mt_numerical_residual_step(imt, ipatch) - aquifer_orphan_mass
            ENDIF
         ENDIF
         CALL check_isotope_aquifer(imt, ipatch, wa, mt_wa(imt, ipatch), 'after qcharge')

         ! ============================================================
         ! 4. Subsurface runoff: post qlayer/qcharge bottom-layer ratio
         !    rsub occurs after infiltration, inter-layer exchange, and
         !    groundwater exchange have updated the bottom layer's material.
         !    Use the current (post-update) bottom-layer concentration.
         ! ============================================================
         IF (rsub > mt_tiny .and. .not. resolved_rsub) THEN
            j = nl_soil
            ! Reconstruct pre-rsub water: WATER's final wliq already has rsub
            ! removed, so add it back to get the state rsub drew from.
            ratio_src = wliq_soisno(j) + rsub * deltim
            IF (ratio_src > mt_water_min_for_ratio) THEN
               ratio_src = mt_wliq_soisno(imt, j, ipatch) / ratio_src
            ELSEIF (material_is_nonvolatile_solute(imt)) THEN
               ratio_src = 0._r8
            ELSE
               ratio_src = ratio_layer(j)
            ENDIF
            mt_flux = rsub * ratio_src * deltim
            mt_flux = min(mt_flux, max(mt_wliq_soisno(imt, j, ipatch), 0._r8))
            mt_wliq_soisno(imt, j, ipatch) = mt_wliq_soisno(imt, j, ipatch) - mt_flux
            a_mt_rsub(imt, ipatch) = a_mt_rsub(imt, ipatch) + mt_flux
               a_mt_rnof(imt, ipatch) = a_mt_rnof(imt, ipatch) + mt_flux
            mt_rnof_step(imt, ipatch) = mt_rnof_step(imt, ipatch) + mt_flux
            water_shadow(j) = water_shadow(j) - rsub * deltim
         ENDIF

         IF (late_surface_water > mt_tiny) THEN
            ! At the post-solver instant, draw the displaced *old* liquid
            ! from the transported top-soil pool, not the pre-WATER ratio.
            IF (frost_surface_water > mt_tiny) THEN
               mt_flux = min(frost_surface_water*current_liq_ratio(1), &
                  max(mt_wliq_soisno(imt,1,ipatch), 0._r8))
               mt_wliq_soisno(imt,1,ipatch) = mt_wliq_soisno(imt,1,ipatch) - mt_flux
               water_shadow(1) = water_shadow(1) - frost_surface_water
               pending_surface_material = pending_surface_material + mt_flux
            ENDIF
            IF (dew_surface_water > mt_tiny) THEN
               mt_flux = dew_surface_water*deposition_ratio_for(layer_temp(1), .false.)
               pending_surface_material = pending_surface_material + mt_flux
               a_mt_precip(imt,ipatch) = a_mt_precip(imt,ipatch) + mt_flux
               a_water_precip(imt,ipatch) = a_water_precip(imt,ipatch) + dew_surface_water
            ENDIF
            late_water = max(wdsrf,0._r8) + late_runoff_water
            late_surface_ratio = 0._r8
            IF (late_water > mt_water_min_for_ratio) THEN
               IF (material_is_nonvolatile_solute(imt) .and. &
                   mt_surface_residue(imt,ipatch) > mt_tiny) THEN
                  pending_surface_material = pending_surface_material + mt_surface_residue(imt,ipatch)
                  mt_surface_residue(imt,ipatch) = 0._r8
               ENDIF
               CALL mt_equilibrate_dissolved(imt, late_water, pending_surface_material, &
                  mt_surface_solid(imt,ipatch))
               late_surface_ratio = pending_surface_material/late_water
            ELSEIF (pending_surface_material > mt_tiny) THEN
               IF (.not. material_is_nonvolatile_solute(imt)) &
                  CALL CoLM_stop('late surface isotope water has no resolved store')
               mt_surface_residue(imt,ipatch) = mt_surface_residue(imt,ipatch) + pending_surface_material
            ENDIF
            mt_wdsrf(imt,ipatch) = max(wdsrf,0._r8)*late_surface_ratio
            mt_flux = late_runoff_water*late_surface_ratio
            a_mt_rsur(imt,ipatch) = a_mt_rsur(imt,ipatch) + mt_flux
            a_mt_rnof(imt,ipatch) = a_mt_rnof(imt,ipatch) + mt_flux
            mt_rnof_step(imt,ipatch) = mt_rnof_step(imt,ipatch) + mt_flux
         ENDIF

         ! ============================================================
         ! 4b. WATER-phase liquid external flux on soil layer 1: dew
         !     qsdew / qsdew_soil are applied directly to wliq_soisno(1)
         !     by WATER, not to the surface mixed pool.
         ! ============================================================
         IF (eff_qsdew_topliq > mt_tiny) THEN
               mt_flux = eff_qsdew_topliq * deposition_ratio_for(layer_temp(1), .false.) * deltim
               mt_wliq_soisno(imt, 1, ipatch) = mt_wliq_soisno(imt, 1, ipatch) + mt_flux
               a_mt_precip(imt, ipatch) = a_mt_precip(imt, ipatch) + mt_flux
               a_water_precip(imt, ipatch) = a_water_precip(imt, ipatch) + eff_qsdew_topliq * deltim
               water_shadow(1) = water_shadow(1) + eff_qsdew_topliq * deltim
            ENDIF

         ! ============================================================
         ! 5. WATER-phase ice external fluxes: qfros_soil / qsubl_soil
         !    These are atmosphere↔soil ice exchanges applied by WATER
         !    directly to wice_soisno(1). They are NOT internal freeze/thaw.
         !
         !    Internal freeze/thaw during WATER is not yet tracked here
         !    (needs separating d_wice into external + internal components).
         !    THERMAL-phase freeze/thaw is handled in mt_evapo.
         ! ============================================================
         ! Record the actual external wice change on soil layer 1 BEFORE
         ! mutating mt_wice, so Section 5b can subtract exactly what the
         ! water side did (including the max(0,...) clamp). Formula:
         !   imperv_subl_loss is an earlier qgtop<0 ice loss from the same
         !   WATER step, so the qfros/qsubl basis is already reduced by it.
         !   qfros_applied = max(qfros, 0) * dt            (always adds)
         !   qsubl_actual  = min(max(qsubl,0)*dt, wice_bef + qfros_applied)
         !   d_wice_ext    = qfros_applied - qsubl_actual
         ! Normal case: qsubl*dt <= wice_bef + qfros*dt  →  d_wice_ext = (qfros-qsubl)*dt
         ! Clamped case: qsubl*dt  > wice_bef + qfros*dt  →  d_wice_ext = -wice_bef
         wice_soil1_after_imperv = max(wice_soisno_bef(1) - imperv_subl_loss, 0._r8)
         d_wice_ext_soil1 = -imperv_subl_loss + max(eff_qfros_top, 0._r8) * deltim &
            - min(max(eff_qsubl_top, 0._r8) * deltim, &
                  wice_soil1_after_imperv + max(eff_qfros_top, 0._r8) * deltim)

         ! Use effective fluxes (handles both split and non-split paths)
         IF (eff_qfros_top > mt_tiny) THEN
            mt_flux = eff_qfros_top * deposition_ratio_for(layer_temp(1), .true.) * deltim
            mt_wice_soisno(imt, 1, ipatch) = mt_wice_soisno(imt, 1, ipatch) + mt_flux
            a_mt_precip(imt, ipatch) = a_mt_precip(imt, ipatch) + mt_flux
            a_water_precip(imt, ipatch) = a_water_precip(imt, ipatch) + eff_qfros_top * deltim
         ENDIF
         IF (eff_qsubl_top > mt_tiny) THEN
            ! mt_wice was just updated by frost above, so the right
            ! denominator is the post-frost ice pool, not wice_soisno_bef(1).
            wice_pre_phase = max(wice_soil1_after_imperv + max(eff_qfros_top, 0._r8) * deltim, 0._r8)
            subl_water = eff_qsubl_top * deltim
            IF (wice_pre_phase - subl_water > mt_water_min_for_ratio) THEN
               ! Normal case: ice covers sublimation.
                  mt_flux = atmospheric_loss_material(mt_wice_soisno(imt, 1, ipatch), &
                     wice_pre_phase, subl_water, layer_temp(1), .true.)
                  mt_wice_soisno(imt, 1, ipatch) = mt_wice_soisno(imt, 1, ipatch) - mt_flux
                  CALL material_book_evap_loss(imt, ipatch, mt_flux, subl_water, &
                     MT_EVAP_KIND_SUBL)
            ELSEIF (subl_water > mt_tiny) THEN
                  CALL exhaust_surface_phase(imt, ipatch, &
                     mt_wice_soisno(imt, 1, ipatch), &
                     min(subl_water, max(wice_pre_phase, 0._r8)), MT_EVAP_KIND_SUBL)
               ENDIF
            ENDIF

         ! ============================================================
         ! 5b. Internal freeze/thaw during WATER
         !
         !     d_wice_total = wice_post - wice_pre (total ice change)
         !     d_wice_external = (qfros_soil - qsubl_soil) * dt  (layer 1 only)
         !     d_wice_internal = d_wice_total - d_wice_external  (pure phase change)
         !
         !     d_wice_internal > 0 → freeze (liquid → ice)
         !     d_wice_internal < 0 → thaw  (ice → liquid)
         !
         ! ACCURACY LIMIT (H8 reconstruction):
         !     The denominators below are reconstructed from POST-WATER state:
         !         pre-freeze wliq = wliq_soisno(j) + d_wice
         !         pre-thaw  wice = wice_soisno(j) + |d_wice|
         !     This recovers the exact pre-phase pool size ONLY if freeze/thaw
         !     was the LAST operation in WATER's substep order. If WATER does
         !     percolation AFTER phase change in some layer, the reconstruction
         !     under/over-states the actual pool at phase-change time.
         !
         !     Non-fractionating materials still use conservative ratio
         !     transfer. Fractionating isotope materials now use a Rayleigh
         !     finite-pool freezing update with this reconstructed pre-freeze
         !     liquid pool. That is exact only if the reconstruction matches
         !     WATER's phase-change intermediate state. If a future WATER_VSF
         !     moves phase change before/among flow steps, pass explicit
         !     phase-change masses or move this material logic into WATER_VSF.
         !     Thaw remains conservative melting of existing ice.
         ! ============================================================
         DO j = lb, nl_soil
            d_wice = wice_soisno(j) - wice_soisno_bef(j)

            ! Subtract external ice fluxes to isolate internal freeze/thaw
            IF (j == 1) THEN
               ! Soil layer 1: use the ACTUAL external change recorded above
               ! in Section 5. Using the nominal (qfros-qsubl)*dt would
               ! misattribute a `max(0, ...)` clamp residual to internal
               ! freezing under extreme sublimation (qsubl*dt > wice_bef+qfros*dt).
               d_wice = d_wice - d_wice_ext_soil1
            ENDIF
            IF (j < 1 .and. j == lb .and. snl < 0) THEN
               ! Top snow layer: external ice flux is applied in snowwater.
               ! Use actual external change recorded in Section 0b, which
               ! correctly accounts for qsubl deficit drawn from liquid.
               d_wice = d_wice - d_wice_ext_snow
            ENDIF

            ! wblc reconciliation in WATER_VSF: the Richards solver ends
            ! with wblc more liquid than its fluxes allow, and the host
            ! removes the same amount of ice. All root and evaporation
            ! demand is already taken by the layers and the aquifer
            ! exchange (and booked from the liquid above), so this is a
            ! liquid/ice closure with no net water, not an ET output. Move
            ! the material of the removed ice into the liquid of the same
            ! layer, where the solver left the extra water; a remainder in
            ! other layers is handled by the shadow residual below. Then
            ! add wblc_ice_sink back into d_wice so the freeze/thaw
            ! classification below sees only the internal phase change.
            IF (j >= 1) THEN
               IF (wblc_ice_sink(j) > mt_tiny) THEN
                  wice_pre_phase = max(wice_soisno(j) + wblc_ice_sink(j), 0._r8)
                  IF (wice_pre_phase > mt_tiny) THEN
                     ratio_src = mt_wice_soisno(imt, j, ipatch) / wice_pre_phase
                     mt_flux = wblc_ice_sink(j) * ratio_src
                     mt_flux = min(mt_flux, max(mt_wice_soisno(imt, j, ipatch), 0._r8))
                     mt_wice_soisno(imt, j, ipatch) = mt_wice_soisno(imt, j, ipatch) - mt_flux
                     mt_wliq_soisno(imt, j, ipatch) = mt_wliq_soisno(imt, j, ipatch) + mt_flux
                     water_shadow(j) = water_shadow(j) + wblc_ice_sink(j)
                  ENDIF
                  d_wice = d_wice + wblc_ice_sink(j)
               ENDIF
            ENDIF

            ! Use the pool size at freeze/thaw time, not pre-WATER. By the time
            ! we reach this block mt_wliq has been updated by percolation /
            ! qcharge / rsub / dew, so wliq_soisno_bef no longer matches the
            ! material pool that participates in the phase change. Reconstruct
            ! pre-freeze wliq = post-WATER wliq + d_wice (the part that froze
            ! out of liquid into ice); pre-thaw wice = post-WATER wice + |d_wice|.
            IF (d_wice > mt_tiny) THEN
               ! Internal freeze: liquid → ice
               wliq_pre_phase = max(wliq_soisno(j) + d_wice, 0._r8)
               IF (wliq_pre_phase > mt_tiny) THEN
                  mt_flux = material_rayleigh_freezing_loss(imt, &
                     mt_wliq_soisno(imt, j, ipatch), wliq_pre_phase, &
                     d_wice, layer_temp(j))
                  mt_flux = min(mt_flux, max(mt_wliq_soisno(imt, j, ipatch), 0._r8))
                  mt_wliq_soisno(imt, j, ipatch) = mt_wliq_soisno(imt, j, ipatch) - mt_flux
                  mt_wice_soisno(imt, j, ipatch) = mt_wice_soisno(imt, j, ipatch) + mt_flux
                  IF (j >= 1) water_shadow(j) = water_shadow(j) - d_wice
               ENDIF
            ELSEIF (d_wice < -mt_tiny) THEN
               ! Internal thaw: ice → liquid
               wice_pre_phase = max(wice_soisno(j) - d_wice, 0._r8)   ! d_wice < 0, so -d_wice > 0
               IF (wice_pre_phase > mt_tiny) THEN
                  ratio_src = mt_wice_soisno(imt, j, ipatch) / wice_pre_phase
                  mt_flux = abs(d_wice) * ratio_src
                     mt_flux = min(mt_flux, max(mt_wice_soisno(imt, j, ipatch), 0._r8))
                     mt_wice_soisno(imt, j, ipatch) = mt_wice_soisno(imt, j, ipatch) - mt_flux
                     mt_wliq_soisno(imt, j, ipatch) = mt_wliq_soisno(imt, j, ipatch) + mt_flux
                     IF (j >= 1) water_shadow(j) = water_shadow(j) - d_wice
                  ENDIF
               ENDIF
            ENDDO

            ! The deferred exfiltration material layer 1 could not pay after the
            ! qlayer transfers. With a frozen top layer that holds almost no
            ! liquid, the host takes the soil evaporation from the ice of
            ! layer 1 (wblc_ice_sink), while above it was booked from the
            ! liquid and the shortfall deferred. The ice-sink transfer just
            ! above has moved that water and its material into the liquid, so
            ! layer 1 can pay now. Paying it after the qlayer transfers only
            ! booked the shortfall as a numerical source and left the same
            ! amount of material in the liquid: a NEC December run then carried
            ! 6.2e-8 of excess H2_18O into the next step's soil evaporation
            ! and stopped on the fixed-signature check. Only a deferral that
            ! is still unpaid here is an explicit numerical source; a
            ! deferral paid in full above leaves nothing to do, so every
            ! other step is unchanged.
            IF (mt_upflow_deferred > 0._r8) THEN
               mt_flux = min(mt_upflow_deferred, max(mt_wliq_soisno(imt, 1, ipatch), 0._r8))
               mt_wliq_soisno(imt, 1, ipatch) = mt_wliq_soisno(imt, 1, ipatch) - mt_flux
               IF (mt_upflow_deferred - mt_flux > 0._r8) THEN
                  IF (allocated(mt_numerical_residual_step)) THEN
                     mt_numerical_residual_step(imt, ipatch) = &
                        mt_numerical_residual_step(imt, ipatch) + (mt_upflow_deferred - mt_flux)
                  ELSE
                     mt_wliq_soisno(imt, 1, ipatch) = mt_wliq_soisno(imt, 1, ipatch) &
                        - (mt_upflow_deferred - mt_flux)
                  ENDIF
               ENDIF
               mt_upflow_deferred = 0._r8
            ENDIF

            ! VSF can redistribute finite water between connected saturated layers
            ! without reporting that transfer in qlayer.  Move its dissolved
            ! material across the omitted faces before accounting for any true
            ! numerical water residual.  A nonclosing column is left to the
            ! explicit residual ledger below.
            CALL reconcile_internal_soil_flow()

            ! VSF can also move water between layers that are not adjacent
            ! (for example infiltration through frozen layers into the first
            ! layer below them) without reporting it in qlayer.  The column
            ! residuals then cancel: layers that lost water against the
            ! shadow and layers that gained it.  Booking both at each layer's
            ! own concentration creates or destroys material whenever the
            ! layers differ (a solute concentrated near the surface, an empty
            ! receiving layer at zero).  Move the part that cancels as a
            ! transfer from the losing to the gaining layers at the losing
            ! layers' concentration; only what does not cancel is booked as a
            ! numerical source/sink below.
            CALL move_matched_layer_residual()

            ! Carry the remaining numerical storage residual at the layer's
            ! current concentration so it cannot appear as fractionation.
            DO j = 1, nl_soil
               water_resid = wliq_soisno(j) - water_shadow(j)
               IF (abs(water_resid) > mt_tiny) THEN
                  water_shadow_ratio = current_liq_ratio(j)
                  soil_resid_water = soil_resid_water + water_resid
                  IF (water_resid >= 0._r8) THEN
                     mt_flux = water_resid * water_shadow_ratio
                     mt_wliq_soisno(imt, j, ipatch) = mt_wliq_soisno(imt, j, ipatch) + mt_flux
                     soil_resid_mt = soil_resid_mt + mt_flux
                  ELSE
                     mt_flux = min(-water_resid * water_shadow_ratio, max(mt_wliq_soisno(imt, j, ipatch), 0._r8))
                     mt_wliq_soisno(imt, j, ipatch) = mt_wliq_soisno(imt, j, ipatch) - mt_flux
                     soil_resid_mt = soil_resid_mt - mt_flux
                  ENDIF
               ENDIF
            ENDDO

               ! The VSF layer residual is a numerical closure term, not a
               ! physical aquifer / runoff / evaporation flux. Keep it out of
               ! mt_wa and expose it to the conservation checker as an
               ! explicit numerical source/sink.
                  IF (abs(soil_resid_mt) > mt_tiny) THEN
                     IF (allocated(mt_numerical_residual_step)) THEN
                        mt_numerical_residual_step(imt, ipatch) = &
                           mt_numerical_residual_step(imt, ipatch) + soil_resid_mt
                        mt_numerical_water_step(imt, ipatch) = &
                           mt_numerical_water_step(imt, ipatch) + soil_resid_water
                     ELSE
                        mt_wa(imt, ipatch) = mt_wa(imt, ipatch) - soil_resid_mt
                     ENDIF
                  ENDIF
                     IF (abs(material_aquifer_actual_water(wa, aquifer_ref_water)) <= mt_water_min_for_ratio) THEN
                        CALL check_isotope_aquifer(imt, ipatch, wa, mt_wa(imt, ipatch), 'soil end')
                        IF (material_is_isotope(imt) .and. aquifer_ref_water > 0._r8) THEN
                           ! The relative mass is not an orphan: -Mref is
                           ! the zero-actual-mass state when the finite
                           ! reference carrier is genuinely exhausted.
                           mt_wa(imt, ipatch) = -aquifer_ref_mass
                        ELSEIF (mt_has_dissolved_limit(imt) .and. &
                            mt_wa(imt, ipatch) > mt_tiny) THEN
                           mt_subsurface_solid(imt, ipatch) = &
                              mt_subsurface_solid(imt, ipatch) + mt_wa(imt, ipatch)
                        ELSEIF (material_is_nonvolatile_solute(imt) .and. &
                               mt_wa(imt, ipatch) > mt_tiny) THEN
                           mt_subsurface_residue(imt, ipatch) = &
                              mt_subsurface_residue(imt, ipatch) + mt_wa(imt, ipatch)
                        ELSEIF (abs(mt_wa(imt, ipatch)) > mt_tiny) THEN
                           ! With exactly zero aquifer water there is no signed
                           ! storage/debt state to mirror. Keep conservation
                           ! explicit, but do not persist an orphan mt_wa.
                           IF (allocated(mt_numerical_residual_step)) THEN
                              mt_numerical_residual_step(imt, ipatch) = &
                                 mt_numerical_residual_step(imt, ipatch) &
                                 - mt_wa(imt, ipatch)
                           ENDIF
                        ENDIF
                        IF (.not. material_is_isotope(imt) .or. aquifer_ref_water <= 0._r8) &
                           mt_wa(imt, ipatch) = 0._r8
                     ENDIF

               ! ============================================================
               ! 6. Wetland water: standard soil path should not mutate it.
            ! ============================================================
            d_wetwat = wetwat - wetwat_bef
#if (defined CoLMDEBUG)
            IF (abs(d_wetwat) > mt_tiny) THEN
               write(*,'(A,I8,A,I3,A,E12.5,A,E12.5,A,E12.5)') &
                  ' WARNING mt_soil_water: wetwat changed in standard soil path ipatch=', &
                  ipatch, ' imt=', imt, ' d_wetwat=', d_wetwat, &
                  ' wetwat_bef=', wetwat_bef, ' wetwat=', wetwat
            ENDIF
#endif

            ! N04: a nonvolatile solute is retained in mt_wice while ice
            ! sublimates (it must not escape with water vapour).  Once the ice
            ! carrier is gone, that inventory is invisible to the liquid/solid
            ! equilibrium and can never re-dissolve on rewetting.  Move it to
            ! the same layer's rewettable liquid/solid inventory, never twice
            ! (the ice pool is cleared), then equilibrate.
            DO j = lb, nl_soil
               IF (material_is_nonvolatile_solute(imt) .and. &
                   mt_wice_soisno(imt, j, ipatch) > mt_tiny) THEN
                  IF (max(wice_soisno(j), 0._r8) <= mt_water_min_for_ratio) THEN
                     IF (max(wliq_soisno(j), 0._r8) > mt_water_min_for_ratio) THEN
                        mt_wliq_soisno(imt, j, ipatch) = &
                           mt_wliq_soisno(imt, j, ipatch) + mt_wice_soisno(imt, j, ipatch)
                     ELSE
                        mt_solid_soisno(imt, j, ipatch) = &
                           mt_solid_soisno(imt, j, ipatch) + mt_wice_soisno(imt, j, ipatch)
                     ENDIF
                     mt_wice_soisno(imt, j, ipatch) = 0._r8
                  ENDIF
               ENDIF
               CALL mt_equilibrate_dissolved(imt, max(wliq_soisno(j), 0._r8), &
                  mt_wliq_soisno(imt, j, ipatch), mt_solid_soisno(imt, j, ipatch))
            ENDDO
            CALL mt_equilibrate_dissolved(imt, wa, mt_wa(imt, ipatch), &
               mt_subsurface_solid(imt, ipatch))
            CALL mt_equilibrate_dissolved(imt, max(wdsrf, 0._r8), &
               mt_wdsrf(imt, ipatch), mt_surface_solid(imt, ipatch))
            IF (allocated(mt_waterstorage) .and. present(waterstorage_patch)) THEN
               CALL mt_equilibrate_dissolved(imt, &
                  max(waterstorage_patch - max(qflx_irrig_ground, 0._r8) * deltim, 0._r8), &
                  mt_waterstorage(imt, ipatch), mt_waterstorage_solid(imt, ipatch))
            ENDIF

         ENDDO

   CONTAINS

#ifdef MATERIAL_SIGNED_SOIL_DEBT_EXPERIMENT
      SUBROUTINE check_signed_soil_experiment()
         integer :: species
         real(r8) :: closed_inputs(28)

         IF (snl /= 0 .or. DEF_ISO_SOIL_DIFFUSION .or. DEF_ISO_SOIL_VAPOR_DIFFUSION) &
            CALL CoLM_stop('signed soil debt: rejected snow or diffusion')
         DO species = 1, nmts
            IF (.not. mt_uses_generic_water(species) .or. material_has_decay(species) .or. &
                materials(species)%decay_rate /= 0._r8) &
               CALL CoLM_stop('signed soil debt: rejected non-passive species')
            ! Production keeps solutes on the clamped, unmixed sweep; signed transport would change that.
            IF (.not. material_is_isotope(species)) &
               CALL CoLM_stop('signed soil debt: rejected non-isotope species')
            IF (material_fractionation_active(species)) &
               CALL CoLM_stop('signed soil debt: rejected fractionation')
         ENDDO
         closed_inputs = [qinfl,qcharge,rsur,rsub,qseva_in,qsdew_in,qsubl_in,qfros_in, &
            qseva_soil,qsdew_soil,qsubl_soil,qfros_soil,qseva_snow,qsdew_snow, &
            qsubl_snow,qfros_snow,sm,fsno,wa,wa_bef,wdsrf,wdsrf_bef,wetwat,wetwat_bef, &
            pg_rain,pg_snow,etroot_aquifer,qflx_irrig_ground]
         IF (.not. all(ieee_is_finite(closed_inputs)) .or. &
             .not. all(ieee_is_finite(wice_soisno)) .or. &
             .not. all(ieee_is_finite(wice_soisno_bef)) .or. &
             .not. all(ieee_is_finite(wblc_ice_sink)) .or. &
             .not. all(ieee_is_finite(etroot_actual))) &
            CALL CoLM_stop('signed soil debt: rejected nonfinite coupling')
         IF (.not. all(ieee_is_finite(qlayer)) .or. .not. ieee_is_finite(deltim)) &
            CALL CoLM_stop('signed soil debt: rejected nonfinite flux')
         IF (deltim <= 0._r8 .or. any(abs(qlayer([0,nl_soil])) > mt_tiny) .or. &
             any(abs(closed_inputs) > mt_tiny) .or. &
             any(abs(wice_soisno) > mt_tiny) .or. any(abs(wice_soisno_bef) > mt_tiny) .or. &
             any(abs(wblc_ice_sink) > mt_tiny) .or. any(abs(etroot_actual) > mt_tiny)) &
            CALL CoLM_stop('signed soil debt: rejected nonclosed or nonliquid step')
         IF (deltim > 1._r8) THEN
            IF (any(abs(qlayer) > huge(1._r8)/deltim)) &
               CALL CoLM_stop('signed soil debt: rejected flux overflow')
         ENDIF
         IF (present(qgtop_solver)) THEN
            IF (.not. ieee_is_finite(qgtop_solver)) &
               CALL CoLM_stop('signed soil debt: rejected nonfinite coupling')
            IF (abs(qgtop_solver) > mt_tiny) CALL CoLM_stop('signed soil debt: rejected boundary flux')
         ENDIF
         IF (present(flood_mt_input) .or. present(flood_infil_water) .or. &
             present(waterstorage_patch) .or. present(etroot_surface) .or. present(dew_overflow) .or. &
             present(frost_displaced) .or. present(late_surface_runoff) .or. &
             present(rsub_source_layer) .or. present(rsub_source_surface) .or. present(rsub_source_aquifer) .or. &
             present(imperv_evap_wdsrf) .or. present(imperv_evap_soil) .or. present(imperv_subl_soil)) &
            CALL CoLM_stop('signed soil debt: rejected production coupling arguments')
         IF (any(wliq_soisno < -soil_debt_limit) .or. &
             any(wliq_soisno_bef < -soil_debt_limit)) &
            CALL CoLM_stop('signed soil debt: rejected oversized carrier debt')
      END SUBROUTINE check_signed_soil_experiment

      SUBROUTINE signed_soil_advection()
         integer :: species, level, donor, receiver
         real(r8) :: water(nl_soil), amount(nl_soil), signature(nl_soil)
         real(r8) :: water_bef(nl_soil), amount_bef(nl_soil)
         real(r8) :: inflow(nl_soil), inflow_mass(nl_soil), passed(nl_soil), passed_water(nl_soil)
         logical  :: unresolved(nl_soil), borrows(nl_soil)
         real(r8) :: lower_ratio, upper_ratio, fallback, volume, moved, mismatch
         real(r8) :: tolerance, low_mass, high_mass, settled_amount

         DO species = 1, nmts
            water_bef = wliq_soisno_bef(1:nl_soil)
            amount_bef = mt_wliq_soisno(species,1:nl_soil,ipatch)
            water = water_bef
            amount = amount_bef
            fallback = material_init_water_ratio(species)
            DO level = nl_soil, 1, -1
               IF (abs(water(level)) <= mt_water_min_for_ratio) CYCLE
               fallback = signed_soil_ratio(water(level), amount(level), fallback)
               EXIT
            ENDDO
            DO level = 1, nl_soil
               signature(level) = signed_soil_ratio(water(level), amount(level), fallback)
            ENDDO
            lower_ratio = minval(signature)
            upper_ratio = maxval(signature)
            ! Reject ratios whose product with any water or face volume below could overflow.
            IF (upper_ratio > (huge(1._r8)/32._r8)/max(1._r8, maxval(abs(water_bef)), &
                maxval(abs(qlayer(1:nl_soil-1))*deltim))) &
               CALL CoLM_stop('signed soil debt: rejected ratio overflow')
            ! An unresolved layer has no composition; what it passes on takes the receiver's.
            unresolved = abs(water_bef) <= mt_water_min_for_ratio
            borrows = unresolved
            ! Include the same-step supplier before a thin layer passes water onward.
            DO level = nl_soil-1, 2, -1
               IF (qlayer(level) >= -mt_tiny .or. qlayer(level-1) >= -mt_tiny) CYCLE
               IF (abs(qlayer(level-1))*deltim <= water(level)) CYCLE
               signature(level) = signed_through_ratio(water(level), amount(level), signature(level), &
                  abs(qlayer(level))*deltim, &
                  merge(signature(level), signature(level+1), borrows(level+1) .and. .not. borrows(level)))
               borrows(level) = .false.
            ENDDO
            DO level = 2, nl_soil-1
               IF (qlayer(level-1) <= mt_tiny .or. qlayer(level) <= mt_tiny) CYCLE
               IF (qlayer(level)*deltim <= water(level)) CYCLE
               signature(level) = signed_through_ratio(water(level), amount(level), signature(level), &
                  qlayer(level-1)*deltim, &
                  merge(signature(level), signature(level-1), borrows(level-1) .and. .not. borrows(level)))
               borrows(level) = .false.
            ENDDO
            inflow = 0._r8
            inflow_mass = 0._r8
            passed = abs(amount_bef)
            passed_water = abs(water_bef)
            DO level = 1, nl_soil-1
               volume = abs(qlayer(level))*deltim
               IF (qlayer(level) >= 0._r8) THEN
                  donor = level
                  receiver = level+1
               ELSE
                  donor = level+1
                  receiver = level
               ENDIF
               moved = volume*merge(signature(receiver), signature(donor), &
                  borrows(donor) .and. .not. borrows(receiver))
               amount(donor) = amount(donor)-moved
               amount(receiver) = amount(receiver)+moved
               water(donor) = water(donor)-volume
               water(receiver) = water(receiver)+volume
               inflow(receiver) = inflow(receiver)+volume
               inflow_mass(receiver) = inflow_mass(receiver)+moved
               passed(level:level+1) = passed(level:level+1)+moved
               passed_water(level:level+1) = passed_water(level:level+1)+volume
            ENDDO
            ! Inflow repays a debt first, at the debt's own composition. Booking the
            ! mismatch keeps the rest a nonnegative mixture instead of an extrapolation.
            DO level = 1, nl_soil
               IF (water_bef(level) >= 0._r8 .or. inflow(level) <= 0._r8) CYCLE
               IF (inflow(level) >= -water_bef(level)) THEN
                  mismatch = inflow_mass(level)*(-water_bef(level)/inflow(level)) + amount_bef(level)
               ELSE
                  mismatch = inflow_mass(level) + amount_bef(level)* &
                     (1._r8-(water_bef(level)+inflow(level))/water_bef(level))
               ENDIF
               amount(level) = amount(level)-mismatch
               passed(level) = passed(level)+abs(mismatch)
               IF (allocated(mt_numerical_residual_step)) &
                  mt_numerical_residual_step(species,ipatch) = mt_numerical_residual_step(species,ipatch)-mismatch
            ENDDO
            DO level = 1, nl_soil
               IF (abs(water(level)-wliq_soisno(level)) > &
                   64._r8*epsilon(1._r8)*max(1._r8,passed_water(level))) &
                  CALL CoLM_stop('signed soil debt: rejected unreported water transfer')
               IF (.not. ieee_is_finite(amount(level))) &
                  CALL CoLM_stop('signed soil debt: rejected nonfinite inventory')
               ! Bound the mass, not mass/water, so rounding over a tiny carrier cannot reject.
               low_mass = min(lower_ratio*water(level), upper_ratio*water(level))
               high_mass = max(lower_ratio*water(level), upper_ratio*water(level))
               tolerance = 64._r8*epsilon(1._r8)*(passed(level)+abs(water(level))*upper_ratio)
               IF (unresolved(level)) tolerance = tolerance+abs(water_bef(level))*upper_ratio
               IF (amount(level) < low_mass-tolerance .or. amount(level) > high_mass+tolerance) &
                  CALL CoLM_stop('signed soil debt: rejected concentration overshoot')
               ! Book that rounding so the next snapshot is sign-compatible and in range.
               settled_amount = merge(0._r8, min(max(amount(level),low_mass),high_mass), &
                  abs(water(level)) <= mt_water_min_for_ratio)
               IF (allocated(mt_numerical_residual_step)) &
                  mt_numerical_residual_step(species,ipatch) = mt_numerical_residual_step(species,ipatch) &
                     + (settled_amount-amount(level))
               amount(level) = settled_amount
            ENDDO
            mt_wliq_soisno(species,1:nl_soil,ipatch) = amount
         ENDDO
      END SUBROUTINE signed_soil_advection

      real(r8) FUNCTION signed_through_ratio(water, amount, own_ratio, volume, inflow_ratio)
         ! Inflow repays debt first, so only nonnegative stock mixes with what is left.
         real(r8), intent(in) :: water, amount, own_ratio, volume, inflow_ratio

         IF (abs(water) <= mt_water_min_for_ratio .or. (water < 0._r8 .and. volume >= -water)) THEN
            signed_through_ratio = inflow_ratio
         ELSEIF (water < 0._r8) THEN
            signed_through_ratio = own_ratio
         ELSE
            signed_through_ratio = signed_soil_ratio(water+volume, amount+volume*inflow_ratio, own_ratio)
         ENDIF
      END FUNCTION signed_through_ratio

      real(r8) FUNCTION signed_soil_ratio(water, amount, fallback)
         real(r8), intent(in) :: water, amount, fallback

         IF (.not. ieee_is_finite(water) .or. .not. ieee_is_finite(amount)) &
            CALL CoLM_stop('signed soil debt: rejected nonfinite inventory')
         IF (water < -soil_debt_limit .or. (water < 0._r8 .and. amount > 0._r8) .or. &
             (water >= 0._r8 .and. amount < 0._r8)) &
            CALL CoLM_stop('signed soil debt: rejected incompatible inventory signs')
         IF (abs(water) <= mt_water_min_for_ratio) THEN
            IF (abs(amount) > mt_tiny) CALL CoLM_stop('signed soil debt: rejected unresolved carrier')
            IF (.not. ieee_is_finite(fallback)) CALL CoLM_stop('signed soil debt: rejected invalid fallback')
            IF (fallback < 0._r8) CALL CoLM_stop('signed soil debt: rejected invalid fallback')
            signed_soil_ratio = fallback
            RETURN
         ENDIF
         IF (abs(water) < 1._r8) THEN
            IF (abs(amount) > huge(1._r8)*abs(water)) &
               CALL CoLM_stop('signed soil debt: rejected ratio overflow')
         ENDIF
         signed_soil_ratio = amount/water
         IF (.not. ieee_is_finite(signed_soil_ratio)) &
            CALL CoLM_stop('signed soil debt: rejected invalid concentration')
         IF (signed_soil_ratio < 0._r8) CALL CoLM_stop('signed soil debt: rejected invalid concentration')
      END FUNCTION signed_soil_ratio
#endif

      SUBROUTINE reconcile_internal_soil_flow()
         ! For each closed permeable block, the unreported face volume is fixed
         ! by the cumulative layer water mismatch.  Downward faces are
         ! processed top-down and upward faces bottom-up, so incoming water
         ! and dissolved material are available before a layer passes them on.
         integer :: iface, donor, receiver, block_begin, block_end
         real(r8) :: column_resid, closure_tol, cumulative, moved_water
         logical :: feasible

         IF (.not. present(permeable_soil)) RETURN
         block_begin = 1
         DO WHILE (block_begin <= nl_soil)
            IF (.not. permeable_soil(block_begin)) THEN
               block_begin = block_begin + 1
               CYCLE
            ENDIF
            block_end = block_begin
            DO WHILE (block_end < nl_soil)
               IF (.not. permeable_soil(block_end+1)) EXIT
               block_end = block_end + 1
            ENDDO
            IF (block_end == block_begin) THEN
               block_begin = block_end + 1
               CYCLE
            ENDIF
            IF (.not. all(ieee_is_finite(wliq_soisno(block_begin:block_end))) .or. &
                .not. all(ieee_is_finite(water_shadow(block_begin:block_end)))) THEN
               block_begin = block_end + 1
               CYCLE
            ENDIF
            IF (any(wliq_soisno(block_begin:block_end) < 0._r8) .or. &
                any(water_shadow(block_begin:block_end) < 0._r8)) THEN
               block_begin = block_end + 1
               CYCLE
            ENDIF
            closure_tol = 64._r8 * epsilon(1._r8) * &
               max(1._r8, sum(abs(wliq_soisno(block_begin:block_end))) + &
                            sum(abs(water_shadow(block_begin:block_end))))
            column_resid = sum(wliq_soisno(block_begin:block_end) - water_shadow(block_begin:block_end))
            IF (abs(column_resid) > closure_tol .or. &
                maxval(abs(wliq_soisno(block_begin:block_end) - water_shadow(block_begin:block_end))) &
                   <= closure_tol) THEN
               block_begin = block_end + 1
               CYCLE
            ENDIF

            cumulative = 0._r8
            DO iface = block_begin, block_end-1
               cumulative = cumulative + wliq_soisno(iface) - water_shadow(iface)
               remap_face_water(iface) = -cumulative
            ENDDO

            ! Reject an inconsistent closed-block reconstruction before any
            ! material mutation.  The original explicit residual ledger remains
            ! responsible; no face volume is silently capped.
            remap_trial_water(block_begin:block_end) = water_shadow(block_begin:block_end)
            feasible = .true.
            DO iface = block_begin, block_end-1
               IF (remap_face_water(iface) <= 0._r8) CYCLE
               IF (remap_face_water(iface) > remap_trial_water(iface)) THEN
                  feasible = .false.
                  EXIT
               ENDIF
               remap_trial_water(iface) = remap_trial_water(iface) - remap_face_water(iface)
               remap_trial_water(iface+1) = remap_trial_water(iface+1) + remap_face_water(iface)
            ENDDO
            IF (feasible) THEN
               DO iface = block_end-1, block_begin, -1
                  IF (remap_face_water(iface) >= 0._r8) CYCLE
                  moved_water = -remap_face_water(iface)
                  IF (moved_water > remap_trial_water(iface+1)) THEN
                     feasible = .false.
                     EXIT
                  ENDIF
                  remap_trial_water(iface+1) = remap_trial_water(iface+1) - moved_water
                  remap_trial_water(iface) = remap_trial_water(iface) + moved_water
               ENDDO
            ENDIF
            IF (feasible) THEN
               DO iface = block_begin, block_end-1
                  IF (remap_face_water(iface) <= 0._r8) CYCLE
                  donor = iface
                  receiver = iface+1
                  moved_water = remap_face_water(iface)
                  CALL move_dissolved_face(donor, receiver, moved_water)
               ENDDO
               DO iface = block_end-1, block_begin, -1
                  IF (remap_face_water(iface) >= 0._r8) CYCLE
                  donor = iface+1
                  receiver = iface
                  moved_water = -remap_face_water(iface)
                  CALL move_dissolved_face(donor, receiver, moved_water)
               ENDDO
            ENDIF
            block_begin = block_end + 1
         ENDDO
      END SUBROUTINE reconcile_internal_soil_flow

      SUBROUTINE move_matched_layer_residual()
         ! The losing layers give their share of the matched volume at their
         ! own concentration (limited by what they hold); the gaining layers
         ! share the pooled material by the volume they receive.  Material
         ! is conserved exactly and nothing is booked.
         real(r8) :: resid(1:nl_soil), lost, gained, matched, w, pooled, moved
         integer  :: jl

         DO jl = 1, nl_soil
            resid(jl) = wliq_soisno(jl) - water_shadow(jl)
         ENDDO
         lost   = -sum(min(resid, 0._r8), mask = resid < -mt_tiny)
         gained =  sum(max(resid, 0._r8), mask = resid >  mt_tiny)
         matched = min(lost, gained)
         IF (matched <= mt_tiny) RETURN

         pooled = 0._r8
         DO jl = 1, nl_soil
            IF (resid(jl) >= -mt_tiny) CYCLE
            w = -resid(jl) * (matched / lost)
            moved = min(w * max(current_liq_ratio(jl), 0._r8), &
               max(mt_wliq_soisno(imt, jl, ipatch), 0._r8))
            mt_wliq_soisno(imt, jl, ipatch) = mt_wliq_soisno(imt, jl, ipatch) - moved
            water_shadow(jl) = water_shadow(jl) - w
            pooled = pooled + moved
         ENDDO
         DO jl = 1, nl_soil
            IF (resid(jl) <= mt_tiny) CYCLE
            w = resid(jl) * (matched / gained)
            mt_wliq_soisno(imt, jl, ipatch) = mt_wliq_soisno(imt, jl, ipatch) + pooled * (w / matched)
            water_shadow(jl) = water_shadow(jl) + w
         ENDDO
      END SUBROUTINE move_matched_layer_residual

      SUBROUTINE move_dissolved_face(donor, receiver, moved_water)
         integer, intent(in) :: donor, receiver
         real(r8), intent(in) :: moved_water
         real(r8) :: moved_material

         CALL mt_equilibrate_dissolved(imt, water_shadow(donor), &
            mt_wliq_soisno(imt, donor, ipatch), mt_solid_soisno(imt, donor, ipatch))
         moved_material = min(moved_water / water_shadow(donor), 1._r8) * &
            max(mt_wliq_soisno(imt, donor, ipatch), 0._r8)
         mt_wliq_soisno(imt, donor, ipatch) = mt_wliq_soisno(imt, donor, ipatch) - moved_material
         mt_wliq_soisno(imt, receiver, ipatch) = mt_wliq_soisno(imt, receiver, ipatch) + moved_material
         water_shadow(donor) = water_shadow(donor) - moved_water
         water_shadow(receiver) = water_shadow(receiver) + moved_water
      END SUBROUTINE move_dissolved_face

      real(r8) FUNCTION layer_temp (jlay)
         integer, intent(in) :: jlay

         layer_temp = 273.15_r8
         IF (present(t_soisno_frac)) layer_temp = t_soisno_frac(jlay)
      END FUNCTION layer_temp

      real(r8) FUNCTION current_liq_ratio (jlay)
         integer, intent(in) :: jlay

         IF (water_shadow(jlay) > mt_water_min_for_ratio) THEN
            current_liq_ratio = mt_wliq_soisno(imt, jlay, ipatch) / water_shadow(jlay)
         ELSEIF (material_is_nonvolatile_solute(imt)) THEN
            current_liq_ratio = 0._r8
         ELSE
            current_liq_ratio = ratio_layer(jlay)
         ENDIF
      END FUNCTION current_liq_ratio

      real(r8) FUNCTION atmospheric_loss_material (pool_mt, pool_water, water_loss, temp_k, from_ice)
         real(r8), intent(in) :: pool_mt
         real(r8), intent(in) :: pool_water
         real(r8), intent(in) :: water_loss
         real(r8), intent(in) :: temp_k
         logical,  intent(in) :: from_ice

         ! from_ice routes through the surface-skin limiter: a snow layer is
         ! not a well-mixed reservoir for sublimation.
         atmospheric_loss_material = material_atmospheric_mt_loss(pool_mt, pool_water, &
            water_loss, temp_k, from_ice, evap_affine_for, mt_tiny, &
            merge(materials(imt)%ref_ratio * (1._r8 + mt_delta_sanity_max / 1000._r8), &
                  0._r8, material_fractionation_active(imt)), material_is_nonvolatile_solute(imt), &
            skin_mass = DEF_ISO_SUBL_SKIN_MM)
      END FUNCTION atmospheric_loss_material

      ! Humidity and kinetic factor of the Craig-Gordon evaporate, shared by
      ! the ratio form (evap_ratio_for) and the affine form (evap_affine_for).
      SUBROUTINE evap_cg_inputs (temp_k, from_ice, active, relhum, alpha_k)
         real(r8), intent(in)  :: temp_k
         logical,  intent(in)  :: from_ice
         logical,  intent(out) :: active
         real(r8), intent(out) :: relhum, alpha_k

         relhum = 0._r8
         alpha_k = 1._r8
         active = material_fractionation_active(imt) .and. &
            present(forc_q_frac) .and. present(forc_psrf_frac)
         IF (.not. active) RETURN

         relhum = material_surface_relhum(forc_q_frac, forc_psrf_frac, temp_k, from_ice)
         ! Bare-soil liquid evaporation: weight the turbulent and pore-diffusion
         ! exponents by the actual resistances instead of pinning n = 2/3.
         ! Snow-layer and surface-water evaporation keep the exponent form --
         ! rss is a soil pore resistance and has no meaning for them.
         IF (kinetic_on_soil_surface .and. .not. from_ice .and. &
             trim(DEF_ISO_SOIL_KINETIC) == 'RESISTANCE' .and. &
             present(ra_frac) .and. present(rss_frac)) THEN
            alpha_k = material_alpha_kinetic_soil(imt, ra_frac, rss_frac)
         ELSE
            alpha_k = material_alpha_kinetic_craig_gordon(imt, from_ice)
         ENDIF
      END SUBROUTINE evap_cg_inputs

      ! Evaporate ratio of water that does not come from a finite pool (the
      ! layer-1 evaporation deficit drawn from below).
      real(r8) FUNCTION evap_ratio_for (source_ratio, temp_k, from_ice)
         real(r8), intent(in) :: source_ratio
         real(r8), intent(in) :: temp_k
         logical,  intent(in) :: from_ice
         real(r8) :: relhum, alpha_k
         logical  :: active

         evap_ratio_for = source_ratio
         CALL evap_cg_inputs(temp_k, from_ice, active, relhum, alpha_k)
         IF (.not. active) RETURN
         evap_ratio_for = material_craig_gordon_evap_ratio(imt, source_ratio, R_atm, &
            temp_k, relhum, alpha_k, from_ice)
      END FUNCTION evap_ratio_for

      ! Affine evaporate ratio of a finite pool, for atmospheric_loss_material.
      SUBROUTINE evap_affine_for (temp_k, from_ice, slope, offset)
         real(r8), intent(in)  :: temp_k
         logical,  intent(in)  :: from_ice
         real(r8), intent(out) :: slope, offset
         real(r8) :: relhum, alpha_k
         logical  :: active

         slope  = 1._r8
         offset = 0._r8
         CALL evap_cg_inputs(temp_k, from_ice, active, relhum, alpha_k)
         IF (.not. active) RETURN
         CALL material_craig_gordon_evap_affine(imt, R_atm, temp_k, relhum, &
            alpha_k, from_ice, slope, offset)
      END SUBROUTINE evap_affine_for

         real(r8) FUNCTION deposition_ratio_for (temp_k, from_ice)
            real(r8), intent(in) :: temp_k
            logical,  intent(in) :: from_ice

            deposition_ratio_for = R_atm
            ! Nonvolatile solutes do not enter through dew/frost
            ! water-vapour deposition.  They require an explicit deposition
            ! forcing/source instead of borrowing the atmospheric water ratio.
            IF (material_is_nonvolatile_solute(imt)) THEN
               deposition_ratio_for = 0._r8
               RETURN
            ENDIF
            IF (.not. material_fractionation_active(imt)) RETURN
            deposition_ratio_for = material_equilibrium_deposition_ratio(imt, R_atm, temp_k, from_ice)
         END FUNCTION deposition_ratio_for

   END SUBROUTINE mt_soil_water

   SUBROUTINE release_leaf_iso_storage (imt, ipatch, nl_soil, wliq_soil, wa_liq)
      integer,  intent(in) :: imt, ipatch, nl_soil
      real(r8), intent(in) :: wliq_soil(1:nl_soil), wa_liq
      integer :: j
      real(r8) :: anomaly, pool_total, release_fraction
      real(r8) :: aquifer_ref_water, aquifer_ref_mass, aquifer_water, aquifer_mass

      aquifer_ref_water = 0._r8
      aquifer_ref_mass = 0._r8
      IF (material_is_isotope(imt)) THEN
         aquifer_ref_water = mt_aquifer_ref_water(ipatch)
         aquifer_ref_mass = mt_aquifer_ref_mass(imt, ipatch)
      ENDIF
      aquifer_water = material_aquifer_actual_water(wa_liq, aquifer_ref_water)
      aquifer_mass = material_aquifer_actual_mass(mt_wa(imt, ipatch), aquifer_ref_mass)

      anomaly = mt_leaf_iso_storage(imt, ipatch)

      ! Leaf-off invalidates the prior NSS state even when a negative
      ! anomaly cannot yet be fully absorbed by an empty root zone.
      mt_leaf_delta_e(imt, ipatch) = 0._r8
      mt_leaf_delta_b(imt, ipatch) = 0._r8
      mt_leaf_peclet(imt, ipatch) = 1._r8
      mt_leaf_water_moles(imt, ipatch) = 0._r8
      IF (abs(anomaly) <= mt_tiny) THEN
         mt_leaf_iso_storage(imt, ipatch) = 0._r8
         RETURN
      ENDIF

      IF (anomaly > 0._r8) THEN
         ! Mix an enriched anomaly across existing liquid water. This is
         ! isotope exchange only: the virtual NSS pool carries no water.
         pool_total = sum(max(wliq_soil, 0._r8)) + max(aquifer_water, 0._r8)
         IF (pool_total <= mt_tiny) RETURN
         DO j = 1, nl_soil
            IF (wliq_soil(j) > 0._r8) THEN
               mt_wliq_soisno(imt, j, ipatch) = mt_wliq_soisno(imt, j, ipatch) + &
                  anomaly * wliq_soil(j) / pool_total
            ENDIF
         ENDDO
         IF (aquifer_water > 0._r8) THEN
            mt_wa(imt, ipatch) = mt_wa(imt, ipatch) + anomaly * aquifer_water / pool_total
         ENDIF
         mt_leaf_iso_storage(imt, ipatch) = 0._r8
      ELSE
         ! A depleted anomaly must remove material from the receiving liquid.
         ! Scale all nonnegative pools together so none can cross zero. If
         ! the column lacks enough material, retain the unreleased remainder
         ! and retry on the next leaf-off step.
         pool_total = sum(max(mt_wliq_soisno(imt, 1:nl_soil, ipatch), 0._r8))
         IF (aquifer_water > 0._r8) pool_total = pool_total + max(aquifer_mass, 0._r8)
         IF (pool_total <= mt_tiny) RETURN
         release_fraction = min(-anomaly / pool_total, 1._r8)
         DO j = 1, nl_soil
            IF (mt_wliq_soisno(imt, j, ipatch) > 0._r8) THEN
               mt_wliq_soisno(imt, j, ipatch) = &
                  mt_wliq_soisno(imt, j, ipatch) * (1._r8 - release_fraction)
            ENDIF
         ENDDO
         IF (aquifer_water > 0._r8 .and. aquifer_mass > 0._r8) THEN
            mt_wa(imt, ipatch) = mt_wa(imt, ipatch) - aquifer_mass * release_fraction
         ENDIF
         mt_leaf_iso_storage(imt, ipatch) = anomaly + release_fraction * pool_total
         IF (abs(mt_leaf_iso_storage(imt, ipatch)) <= mt_tiny) &
            mt_leaf_iso_storage(imt, ipatch) = 0._r8
      ENDIF
   END SUBROUTINE release_leaf_iso_storage

   !---------------------------------------------------------------
   ! Wetland-specific material update (patchtype==2 .and.
   ! .not. DEF_USE_Dynamic_Wetland).
   !
   ! WATER_VSF's wetland branch (MOD_SoilSnowHydrology.F90 L1170+) merges
   ! wdsrf + wa + wetwat + (gwat-etr+dew+frost-subl)*dt + sum(wresi) into a
   ! single wetwat pool, then overflows back to wdsrf / deficits into wa.
   ! The normal soil-water material path (mt_soil_water) assumes an
   ! independent wdsrf/wa/wetwat accounting and silently loses the
   ! pool-to-wetwat mass transfer. This routine mirrors the water merge
   ! exactly:
   !   pool_water  = wdsrf_bef + wa_bef + wetwat_bef + wresi + inputs
   !   pool_material = mt_wdsrf + mt_wa + mt_wetwat + mt_wresi + mt_inputs
   !   loss        = (qseva + qsubl + etr)*dt  removed at pool ratio
   !   remaining   = redistributed to post-WATER wdsrf / wa / wetwat / rsur
   !
   ! Snow layers (lb<0) are handled up front with the same top-layer flux +
   ! percolation logic as mt_soil_water sections 0b/0c, producing
   ! mt_gwat_snow which joins the wetland pool instead of the surface pool.
   !---------------------------------------------------------------
   SUBROUTINE mt_wetland (ipatch, deltim, snl, nl_soil, &
      rsur, &
      qseva_in, qsdew_in, qsubl_in, qfros_in, &
      qseva_soil, qsdew_soil, qsubl_soil, qfros_soil, &
      qseva_snow, qsdew_snow, qsubl_snow, qfros_snow, &
      etr, sm, fsno, split_soilsnow, &
      wliq_soisno, wice_soisno, &
      wliq_soisno_bef, wice_soisno_bef, &
      wa, wa_bef, wdsrf, wdsrf_bef, &
      wetwat, wetwat_bef, pg_rain, pg_snow, &
      t_soisno, porsl, dz_soisno, qflx_irrig_ground, forc_us_frac, forc_vs_frac, &
      waterstorage_patch, snow_qout_layer, forc_q_frac, forc_psrf_frac, &
      tleaf_frac, lai_frac, rst_frac, ra_frac, dz_sno_frac)

      USE MOD_Const_Physical, only: tfrz
      IMPLICIT NONE
      integer,  intent(in) :: ipatch, snl, nl_soil
      real(r8), intent(in) :: deltim, rsur
      real(r8), intent(in) :: qseva_in, qsdew_in, qsubl_in, qfros_in
      real(r8), intent(in) :: qseva_soil, qsdew_soil, qsubl_soil, qfros_soil
      real(r8), intent(in) :: qseva_snow, qsdew_snow, qsubl_snow, qfros_snow
      real(r8), intent(in) :: etr, sm, fsno
      logical,  intent(in) :: split_soilsnow
      real(r8), intent(in) :: wliq_soisno(snl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno(snl+1:nl_soil)
      real(r8), intent(in) :: wliq_soisno_bef(snl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno_bef(snl+1:nl_soil)
      real(r8), intent(in) :: wa, wa_bef, wdsrf, wdsrf_bef
      real(r8), intent(in) :: wetwat, wetwat_bef, pg_rain, pg_snow
      real(r8), intent(in) :: t_soisno(snl+1:nl_soil)
      real(r8), intent(in) :: porsl(1:nl_soil), dz_soisno(snl+1:nl_soil)
      ! Combined drip + flood + paddy irrigation rate (mm/s). Added to
      ! gwat in WATER_VSF (MOD_SoilSnowHydrology.F90:828-832) ahead of
      ! the wetland-merge branch, so wetwat_post already reflects this
      ! water. Without mirroring it here, pool_material is short by
      ! `qflx_irrig_ground*dt*R_precip` while pool_water is short by
      ! the same volume — the ratio stays roughly right but the
      ! absolute material mass shrinks and xerr_material grows by the
      ! irrigation input per step.
      real(r8), intent(in) :: qflx_irrig_ground
      real(r8), intent(in) :: forc_us_frac, forc_vs_frac
      real(r8), intent(in), optional :: waterstorage_patch
      real(r8), intent(in), optional :: snow_qout_layer(snl+1:0)
      real(r8), intent(in), optional :: forc_q_frac
      real(r8), intent(in), optional :: forc_psrf_frac
      real(r8), intent(in), optional :: tleaf_frac
      real(r8), intent(in), optional :: lai_frac
      real(r8), intent(in), optional :: rst_frac
      real(r8), intent(in), optional :: ra_frac
      real(r8), intent(in), optional :: dz_sno_frac(snl+1:0)

      integer  :: imt, j, lb, lb_snow
      real(r8) :: R_atm, ratio_src, mt_flux, storage_ratio
      real(r8) :: wresi_j, wresi_sum, mt_wresi_j, mt_wresi_sum
      real(r8) :: q_rain_in, q_dew_in, q_frost_in
      real(r8) :: q_evap_out, q_subl_out, q_etr_out, q_sm_in, q_vapor_uptake
      real(r8) :: pool_water, pool_material, pool_ratio
      real(r8) :: mt_loss, loss_water, mt_rsur_local
      real(r8) :: aquifer_ref_water, aquifer_ref_mass
      real(r8) :: mt_dew_input, mt_frost_input, mt_evap_loss, mt_subl_loss, mt_etr_loss
      real(r8) :: redist_storage, redist_target, redist_resid, redist_fix
      real(r8) :: mt_gwat_snow_local
      ! Water mass (mm) leaving the snow bottom during percolation — the
      ! hydrology side feeds this into wetwat via `gwat*deltim`. Material
      ! pool_water must mirror it so the mixed-pool ratio stays consistent
      ! with the water balance (otherwise mt_gwat_snow_local inflates
      ! pool_material without a matching pool_water contribution).
      real(r8) :: gwat_snow_local
      real(r8) :: eff_qseva_snow, eff_qsdew_snow, eff_qsubl_snow, eff_qfros_snow
      real(r8) :: water_ice_pool, water_ice_pool_prefrost, water_liq_pool
      real(r8) :: subl_water, evap_water, deficit_water
      real(r8) :: snow_rain_input, snow_dew_input, snow_evap_output
      real(r8) :: snow_frost_input, snow_subl_output, snow_wliq_before_flow
      real(r8) :: qin_snow, qout_snow, water_before_flow, mt_before_flow
      real(r8) :: mt_qin_snow, mt_qout_snow, melt_exchange
      real(r8) :: loss_liq_avail, loss_ice_avail
      real(r8) :: pool_water_loss, pool_mt_loss, pool_ratio_loss
      real(r8) :: relhum_leaf, leaf_area_use, rst_use, ra_use, tleaf_use, R_vapor
      real(r8) :: transp_ratio, transp_source_material, transp_output_material
      real(r8) :: leaf_delta_e_new, leaf_delta_b_new, leaf_peclet_new, leaf_moles_new
      real(r8) :: d_eff_up, d_eff_dn, diff_transfer, ratio_up, ratio_dn
      real(r8) :: diff_face(snl+1:nl_soil), diff_out(snl+1:nl_soil), diff_scale(snl+1:nl_soil)
      logical  :: fractionate_pool_loss, transp_frac_active
      IF (nmts <= 0) RETURN
      lb = snl + 1

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         aquifer_ref_water = 0._r8
         aquifer_ref_mass = 0._r8
         IF (material_is_isotope(imt)) THEN
            aquifer_ref_water = mt_aquifer_ref_water(ipatch)
            aquifer_ref_mass = mt_aquifer_ref_mass(imt, ipatch)
         ENDIF
         R_atm = material_forcing_vapor_value(imt, ipatch)
         DO j = lb, nl_soil
            CALL mt_equilibrate_dissolved(imt, max(wliq_soisno_bef(j), 0._r8), &
               mt_wliq_soisno(imt, j, ipatch), mt_solid_soisno(imt, j, ipatch))
         ENDDO
         IF (present(lai_frac)) THEN
            IF (allocated(mt_leaf_iso_storage)) THEN
               IF (lai_frac <= mt_tiny .or. &
                   mt_leaf_water_moles(imt, ipatch) <= mt_tiny) THEN
                  CALL release_leaf_iso_storage(imt, ipatch, nl_soil, &
                     wliq_soisno_bef(1:nl_soil), wa_bef)
               ENDIF
            ENDIF
         ENDIF
         IF (mt_has_dissolved_limit(imt)) THEN
            mt_surface_solid(imt, ipatch) = mt_surface_solid(imt, ipatch) &
               + mt_subsurface_solid(imt, ipatch)
            mt_subsurface_solid(imt, ipatch) = 0._r8
         ENDIF
         CALL check_isotope_aquifer(imt, ipatch, wa_bef, mt_wa(imt, ipatch), 'wetland start')

         !--------------------------------------------------------
         ! 1) Snow-top external fluxes + percolation → mt_gwat_snow
         !    (mirror of mt_soil_water sections 0b / 0c)
         !--------------------------------------------------------
         mt_gwat_snow_local = 0._r8
         gwat_snow_local     = 0._r8
         IF (snl < 0) THEN
            lb_snow = lb

            IF (split_soilsnow) THEN
               eff_qseva_snow = qseva_snow
               eff_qsdew_snow = qsdew_snow
               eff_qsubl_snow = qsubl_snow
               eff_qfros_snow = qfros_snow
            ELSE
               eff_qseva_snow = qseva_in
               eff_qsdew_snow = qsdew_in
               eff_qsubl_snow = qsubl_in
               eff_qfros_snow = qfros_in
            ENDIF

            ! Ice phase: frost in, sublim out
            mt_flux = max(eff_qfros_snow, 0._r8) * &
               deposition_ratio_for(layer_temp(lb_snow), .true.) * deltim
            mt_wice_soisno(imt, lb_snow, ipatch) = &
               mt_wice_soisno(imt, lb_snow, ipatch) + mt_flux
            a_mt_precip(imt, ipatch) = a_mt_precip(imt, ipatch) + mt_flux
            a_water_precip(imt, ipatch) = a_water_precip(imt, ipatch) + max(eff_qfros_snow, 0._r8) * deltim

            water_ice_pool = wice_soisno_bef(lb_snow) + max(eff_qfros_snow, 0._r8) * deltim
            water_ice_pool_prefrost = water_ice_pool
            subl_water = max(eff_qsubl_snow, 0._r8) * deltim

            IF (subl_water > mt_tiny) THEN
               IF (water_ice_pool - subl_water > mt_water_min_for_ratio) THEN
                     mt_flux = atmospheric_loss_material(mt_wice_soisno(imt, lb_snow, ipatch), &
                        water_ice_pool, subl_water, layer_temp(lb_snow), .true.)
                     mt_wice_soisno(imt, lb_snow, ipatch) = &
                        mt_wice_soisno(imt, lb_snow, ipatch) - mt_flux
                     CALL material_book_evap_loss(imt, ipatch, mt_flux, subl_water, &
                        MT_EVAP_KIND_SUBL)
               ELSE
                     CALL exhaust_surface_phase(imt, ipatch, &
                        mt_wice_soisno(imt, lb_snow, ipatch), &
                        min(subl_water, max(water_ice_pool, 0._r8)), MT_EVAP_KIND_SUBL)
                  deficit_water = subl_water - max(water_ice_pool, 0._r8)
                  IF (deficit_water > mt_tiny .and. wliq_soisno_bef(lb_snow) > mt_tiny) THEN
                        mt_flux = atmospheric_loss_material(mt_wliq_soisno(imt, lb_snow, ipatch), &
                           wliq_soisno_bef(lb_snow), deficit_water, layer_temp(lb_snow), .false.)
                        mt_wliq_soisno(imt, lb_snow, ipatch) = &
                           mt_wliq_soisno(imt, lb_snow, ipatch) - mt_flux
                        CALL material_book_evap_loss(imt, ipatch, mt_flux, deficit_water, &
                           MT_EVAP_KIND_SUBL)
                  ENDIF
               ENDIF
            ENDIF
            water_ice_pool = max(water_ice_pool - subl_water, 0._r8)

            ! Liquid phase: rain + dew in, evap out
            IF (split_soilsnow) THEN
               mt_flux = mt_pg_rain_ground(imt, ipatch) * fsno
            ELSE
               mt_flux = mt_pg_rain_ground(imt, ipatch)
            ENDIF
            mt_wliq_soisno(imt, lb_snow, ipatch) = &
               mt_wliq_soisno(imt, lb_snow, ipatch) + mt_flux

            mt_flux = max(eff_qsdew_snow, 0._r8) * &
               deposition_ratio_for(layer_temp(lb_snow), .false.) * deltim
            mt_wliq_soisno(imt, lb_snow, ipatch) = &
               mt_wliq_soisno(imt, lb_snow, ipatch) + mt_flux
            a_mt_precip(imt, ipatch) = a_mt_precip(imt, ipatch) + mt_flux
            a_water_precip(imt, ipatch) = a_water_precip(imt, ipatch) + max(eff_qsdew_snow, 0._r8) * deltim

            IF (eff_qseva_snow > mt_tiny) THEN
               IF (split_soilsnow) THEN
                  water_liq_pool = wliq_soisno_bef(lb_snow) &
                     + (pg_rain*fsno + max(eff_qsdew_snow, 0._r8)) * deltim
               ELSE
                  water_liq_pool = wliq_soisno_bef(lb_snow) &
                     + (pg_rain + max(eff_qsdew_snow, 0._r8)) * deltim
               ENDIF
               IF (water_ice_pool_prefrost < subl_water) THEN
                  water_liq_pool = water_liq_pool - (subl_water - water_ice_pool_prefrost)
               ENDIF
               water_liq_pool = max(water_liq_pool, 0._r8)
               evap_water = eff_qseva_snow * deltim

               IF (water_liq_pool - evap_water > mt_water_min_for_ratio) THEN
                     mt_flux = atmospheric_loss_material(mt_wliq_soisno(imt, lb_snow, ipatch), &
                        water_liq_pool, evap_water, layer_temp(lb_snow), .false.)
                     mt_wliq_soisno(imt, lb_snow, ipatch) = &
                        mt_wliq_soisno(imt, lb_snow, ipatch) - mt_flux
                     CALL material_book_evap_loss(imt, ipatch, mt_flux, evap_water, &
                        MT_EVAP_KIND_SOILEVAP)
               ELSE
                     CALL exhaust_surface_phase(imt, ipatch, &
                        mt_wliq_soisno(imt, lb_snow, ipatch), &
                        min(evap_water, max(water_liq_pool, 0._r8)), MT_EVAP_KIND_SOILEVAP)
                  deficit_water = evap_water - max(water_liq_pool, 0._r8)
                  IF (deficit_water > mt_tiny .and. water_ice_pool > mt_tiny) THEN
                        deficit_water = min(deficit_water, water_ice_pool)
                        mt_flux = atmospheric_loss_material(mt_wice_soisno(imt, lb_snow, ipatch), &
                           water_ice_pool, deficit_water, layer_temp(lb_snow), .true.)
                        mt_wice_soisno(imt, lb_snow, ipatch) = &
                           mt_wice_soisno(imt, lb_snow, ipatch) - mt_flux
                        CALL material_book_evap_loss(imt, ipatch, mt_flux, deficit_water, &
                           MT_EVAP_KIND_SOILEVAP)
                  ENDIF
               ENDIF
            ENDIF

            ! Snow percolation → mt_gwat_snow_local
            IF (split_soilsnow) THEN
               snow_rain_input  = max(pg_rain * fsno, 0._r8) * deltim
               snow_dew_input   = max(qsdew_snow, 0._r8) * deltim
               snow_evap_output = max(qseva_snow, 0._r8) * deltim
               snow_frost_input = max(qfros_snow, 0._r8) * deltim
               snow_subl_output = max(qsubl_snow, 0._r8) * deltim
            ELSE
               snow_rain_input  = max(pg_rain, 0._r8) * deltim
               snow_dew_input   = max(qsdew_in, 0._r8) * deltim
               snow_evap_output = max(qseva_in, 0._r8) * deltim
               snow_frost_input = max(qfros_in, 0._r8) * deltim
               snow_subl_output = max(qsubl_in, 0._r8) * deltim
            ENDIF

            mt_qin_snow = 0._r8
            snow_wliq_before_flow = wliq_soisno_bef(lb_snow)
            IF (wice_soisno_bef(lb_snow) + snow_frost_input - snow_subl_output < 0._r8) THEN
               snow_wliq_before_flow = snow_wliq_before_flow &
                  + (wice_soisno_bef(lb_snow) + snow_frost_input - snow_subl_output)
            ENDIF
            snow_wliq_before_flow = snow_wliq_before_flow &
               + snow_rain_input + snow_dew_input - snow_evap_output
            snow_wliq_before_flow = max(snow_wliq_before_flow, 0._r8)

            qin_snow = 0._r8
            DO j = lb_snow, 0
               IF (j == lb_snow) THEN
                  water_before_flow = max(snow_wliq_before_flow, 0._r8)
                  mt_before_flow = max(mt_wliq_soisno(imt, j, ipatch), 0._r8) + mt_qin_snow
               ELSE
                  water_before_flow = max(wliq_soisno_bef(j), 0._r8) + qin_snow
                  mt_before_flow = max(mt_wliq_soisno(imt, j, ipatch), 0._r8) + mt_qin_snow
               ENDIF

               IF (present(snow_qout_layer)) THEN
                  qout_snow = max(snow_qout_layer(j), 0._r8)
               ELSE
                  qout_snow = max(water_before_flow - max(wliq_soisno(j), 0._r8), 0._r8)
               ENDIF
               qout_snow = min(qout_snow, water_before_flow)

               ! --- meltwater <-> layer ice isotopic exchange ---
               ! Keep the wetland snow-percolation mirror on the same
               ! fractionating internal exchange as the ordinary soil-water
               ! path.  The transfer is internal: equal and opposite between
               ! percolating liquid and layer ice, with target
               ! R_liq = R_ice / alpha_ice_liq.
               IF (DEF_ISO_SNOWMELT_EQUILIBRATION > 0._r8 .and. &
                   material_fractionation_active(imt) .and. &
                   .not. material_is_nonvolatile_solute(imt) .and. &
                   water_before_flow > mt_water_min_for_ratio .and. &
                   wice_soisno(j) > mt_water_min_for_ratio) THEN
                  melt_exchange = material_equilibration_exchange(mt_before_flow, &
                     water_before_flow, &
                     max(mt_wice_soisno(imt, j, ipatch), 0._r8) / wice_soisno(j), &
                     1._r8 / max(material_alpha_ice_liq(imt, layer_temp(j)), mt_tiny), &
                     DEF_ISO_SNOWMELT_EQUILIBRATION)
                  IF (melt_exchange > 0._r8) THEN
                     melt_exchange = min(melt_exchange, &
                        max(mt_wice_soisno(imt, j, ipatch), 0._r8))
                  ELSE
                     melt_exchange = -min(-melt_exchange, max(mt_before_flow, 0._r8))
                  ENDIF
                  mt_before_flow = max(mt_before_flow + melt_exchange, 0._r8)
                  mt_wice_soisno(imt, j, ipatch) = &
                     max(mt_wice_soisno(imt, j, ipatch) - melt_exchange, 0._r8)
               ENDIF

               IF (qout_snow > mt_tiny .and. water_before_flow > mt_water_min_for_ratio &
                   .and. mt_before_flow > mt_tiny) THEN
                  ratio_src = mt_before_flow / water_before_flow
                  mt_qout_snow = qout_snow * ratio_src
                  mt_qout_snow = min(mt_qout_snow, mt_before_flow)
               ELSE
                  mt_qout_snow = 0._r8
               ENDIF

               mt_wliq_soisno(imt, j, ipatch) = max(mt_before_flow - mt_qout_snow, 0._r8)

               qin_snow = qout_snow
               mt_qin_snow = mt_qout_snow
            ENDDO

            mt_gwat_snow_local = mt_qin_snow
            gwat_snow_local     = qin_snow

            ! Wetland snowpacks keep the same firn vapour diffusion as the
            ! ordinary soil/snow path; it is an internal ice exchange and must
            ! run before the wetland pool merge strips meltwater out.
            IF (DEF_ISO_SOIL_VAPOR_DIFFUSION .and. material_is_isotope(imt) &
                .and. present(dz_sno_frac) .and. present(forc_psrf_frac)) THEN
               diff_face = 0._r8
               diff_out = 0._r8
               diff_scale = 1._r8
               DO j = lb, -1
                  IF (dz_sno_frac(j) <= 0._r8 .or. dz_sno_frac(j+1) <= 0._r8) CYCLE
                  IF (wice_soisno(j) <= mt_water_min_for_ratio) CYCLE
                  IF (wice_soisno(j+1) <= mt_water_min_for_ratio) CYCLE
                  d_eff_up = material_isotope_snow_vapor_equivalent_diffusivity( &
                     imt, max(wliq_soisno(j), 0._r8), max(wice_soisno(j), 0._r8), &
                     dz_sno_frac(j), layer_temp(j), forc_psrf_frac)
                  d_eff_dn = material_isotope_snow_vapor_equivalent_diffusivity( &
                     imt, max(wliq_soisno(j+1), 0._r8), max(wice_soisno(j+1), 0._r8), &
                     dz_sno_frac(j+1), layer_temp(j+1), forc_psrf_frac)
                  IF (d_eff_up <= 0._r8 .or. d_eff_dn <= 0._r8) CYCLE
                  ratio_up = mt_wice_soisno(imt, j,   ipatch) / wice_soisno(j)
                  ratio_dn = mt_wice_soisno(imt, j+1, ipatch) / wice_soisno(j+1)
                  diff_face(j) = material_soil_diffusive_transfer(ratio_up, ratio_dn, &
                     wice_soisno(j), wice_soisno(j+1), d_eff_up, d_eff_dn, &
                     dz_sno_frac(j), dz_sno_frac(j+1), deltim)
                  IF (diff_face(j) > 0._r8) THEN
                     diff_out(j) = diff_out(j) + diff_face(j)
                  ELSE
                     diff_out(j+1) = diff_out(j+1) - diff_face(j)
                  ENDIF
               ENDDO
               DO j = lb, 0
                  IF (diff_out(j) > max(mt_wice_soisno(imt, j, ipatch), 0._r8) .and. &
                      diff_out(j) > mt_tiny) THEN
                     diff_scale(j) = max(mt_wice_soisno(imt, j, ipatch), 0._r8) / diff_out(j)
                  ENDIF
               ENDDO
               DO j = lb, -1
                  IF (diff_face(j) > 0._r8) THEN
                     diff_transfer = diff_face(j) * diff_scale(j)
                  ELSE
                     diff_transfer = diff_face(j) * diff_scale(j+1)
                  ENDIF
                  mt_wice_soisno(imt, j,   ipatch) = &
                     mt_wice_soisno(imt, j,   ipatch) - diff_transfer
                  mt_wice_soisno(imt, j+1, ipatch) = &
                     mt_wice_soisno(imt, j+1, ipatch) + diff_transfer
               ENDDO
            ENDIF
         ENDIF   ! snl < 0

         !--------------------------------------------------------
         ! 2) Strip wresi (soil layer excess over saturation) into pool
         !    Mirrors MOD_SoilSnowHydrology.F90 L1187-1192.
         !--------------------------------------------------------
         mt_wresi_sum = 0._r8
         wresi_sum     = 0._r8
         DO j = 1, nl_soil
            IF (t_soisno(j) > tfrz) THEN
               wresi_j = max(wliq_soisno_bef(j) - porsl(j)*dz_soisno(j)*1000._r8, 0._r8)
               IF (wresi_j > mt_tiny) THEN
                  wresi_sum = wresi_sum + wresi_j
                  IF (wliq_soisno_bef(j) > mt_water_min_for_ratio) THEN
                     ratio_src = mt_wliq_soisno(imt, j, ipatch) / wliq_soisno_bef(j)
                     mt_wresi_j = min(wresi_j * ratio_src, &
                        max(mt_wliq_soisno(imt, j, ipatch), 0._r8))
                     mt_wliq_soisno(imt, j, ipatch) = &
                        mt_wliq_soisno(imt, j, ipatch) - mt_wresi_j
                     mt_wresi_sum = mt_wresi_sum + mt_wresi_j
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

         !--------------------------------------------------------
         ! 3) Build mixed wetland pool (water + material) PRE loss
         !
         !    WATER_VSF wetwat formula variants (lb = snl+1):
         !      lb>=1, no split : gwat = pg_rain + sm - qseva
         !                        wetwat += (gwat - etr + qsdew + qfros - qsubl)*dt + wresi
         !      lb>=1, split    : gwat = pg_rain + sm - qseva_soil
         !                        wetwat += (gwat - etr + qsdew_s + qfros_s - qsubl_s)*dt + wresi
         !      lb<0 , no split : gwat = snowwater_output (rain/dew/frost/subl/evap all absorbed)
         !                        wetwat += (gwat - etr)*dt + wresi
         !      lb<0 , split    : gwat = snowwater_out + pg_rain*(1-fsno) - qseva_soil
         !                        wetwat += (gwat - etr + qsdew_s + qfros_s - qsubl_s)*dt + wresi
         !--------------------------------------------------------
         IF (snl < 0) THEN
            IF (split_soilsnow) THEN
               q_rain_in  = max(pg_rain * (1._r8 - fsno), 0._r8) * deltim
               q_dew_in   = max(qsdew_soil, 0._r8) * deltim
               q_frost_in = max(qfros_soil, 0._r8) * deltim
               q_evap_out = max(qseva_soil, 0._r8) * deltim
               q_subl_out = max(qsubl_soil, 0._r8) * deltim
            ELSE
               ! All atmospheric exchange absorbed by snowwater; no extra terms.
               q_rain_in  = 0._r8
               q_dew_in   = 0._r8
               q_frost_in = 0._r8
               q_evap_out = 0._r8
               q_subl_out = 0._r8
            ENDIF
            q_sm_in = 0._r8   ! snowmelt contributes via mt_gwat_snow_local
         ELSE
            IF (split_soilsnow) THEN
               q_rain_in  = max(pg_rain, 0._r8) * deltim
               q_dew_in   = max(qsdew_soil, 0._r8) * deltim
               q_frost_in = max(qfros_soil, 0._r8) * deltim
               q_evap_out = max(qseva_soil, 0._r8) * deltim
               q_subl_out = max(qsubl_soil, 0._r8) * deltim
            ELSE
               q_rain_in  = max(pg_rain, 0._r8) * deltim
               q_dew_in   = max(qsdew_in, 0._r8) * deltim
               q_frost_in = max(qfros_in, 0._r8) * deltim
               q_evap_out = max(qseva_in, 0._r8) * deltim
               q_subl_out = max(qsubl_in, 0._r8) * deltim
            ENDIF
            q_sm_in = max(sm, 0._r8) * deltim
         ENDIF
         q_etr_out = max(etr, 0._r8) * deltim

         ! Use SIGNED arithmetic for wdsrf_bef / wa_bef / wetwat_bef so the
         ! aquifer "debt" (wa<0 from a prior wetwat<0 deficit branch at
         ! MOD_SoilSnowHydrology.F90 L1200-1203) is mirrored exactly. Using
         ! max(.,0) here previously dropped wa_bef<0 from pool_water while
         ! keeping its positive inputs, which silently lost material when the
         ! debt was repaid (e.g. wa_bef=-10, input=+20 -> pool_water was 20
         ! instead of 10, so only half the 20R_atm input ended up in wetwat).
         ! wdsrf / wetwat are non-negative by WATER_VSF construction, but
         ! writing them signed keeps the formula uniform.
         pool_water = wdsrf_bef + wa_bef + wetwat_bef + aquifer_ref_water &
                    + wresi_sum + q_rain_in + q_sm_in + q_dew_in + q_frost_in

         pool_material = mt_wdsrf(imt, ipatch) + mt_wa(imt, ipatch) + aquifer_ref_mass &
                     + mt_wetwat(imt, ipatch) + mt_wresi_sum

         ! Wetland hydrology mixes the aquifer into the signed bulk pool.  If
         ! either the beginning or ending aquifer is physically resolved,
         ! return a dry-period conservative residue before this step's bulk
         ! losses and redistribution.
         IF (material_is_nonvolatile_solute(imt) .and. &
             (wa_bef > mt_water_min_for_ratio .or. &
              wa > mt_water_min_for_ratio) .and. &
             mt_subsurface_residue(imt, ipatch) > mt_tiny) THEN
            pool_material = pool_material + mt_subsurface_residue(imt, ipatch)
            mt_subsurface_residue(imt, ipatch) = 0._r8
         ENDIF

         ! Drip / flood / paddy irrigation (if enabled) joins gwat in
         ! WATER_VSF ahead of the wetland merge, so it lands inside
         ! wetwat_post. Water side treats this as an internal transfer
         ! from the per-patch waterstorage reservoir (CoLMMAIN:806
         ! includes waterstorage in totwb/endwb), so the atmospheric
         ! input term ignores it. The material path mirrors that — debit
         ! mt_waterstorage at R_atm (Phase-1 invariant), add the
         ! matching mass to the wetland pool, and do NOT increment
         ! a_mt_precip. See MOD_MT_SoilWater.F90 main
         ! mt_soil_water branch for the symmetric treatment.
         IF (qflx_irrig_ground > mt_tiny) THEN
            pool_water  = pool_water  + qflx_irrig_ground * deltim
            storage_ratio = R_atm
            IF (allocated(mt_waterstorage) .and. present(waterstorage_patch)) THEN
               CALL mt_equilibrate_dissolved(imt, max(waterstorage_patch, 0._r8), &
                  mt_waterstorage(imt, ipatch), mt_waterstorage_solid(imt, ipatch))
               IF (waterstorage_patch > mt_water_min_for_ratio) THEN
                  storage_ratio = mt_waterstorage(imt, ipatch) / waterstorage_patch
               ENDIF
            ENDIF
            mt_flux    = qflx_irrig_ground * deltim * storage_ratio
            IF (allocated(mt_waterstorage)) THEN
               mt_flux = min(max(mt_flux, 0._r8), max(mt_waterstorage(imt, ipatch), 0._r8))
            ENDIF
            pool_material = pool_material + mt_flux
            IF (allocated(mt_waterstorage)) THEN
               mt_waterstorage(imt, ipatch) = mt_waterstorage(imt, ipatch) - mt_flux
            ENDIF
         ENDIF

         IF (snl < 0) THEN
            ! Rain that reaches wetland pool directly is only (1-fsno) portion
            ! in split mode; in non-split mode all rain went through snow.
            IF (split_soilsnow) THEN
               pool_material = pool_material + mt_pg_rain_ground(imt, ipatch) * (1._r8 - fsno)
            ENDIF
            pool_material = pool_material + mt_gwat_snow_local
            ! Mirror mt_gwat_snow_local with its water mass so
            ! pool_ratio = pool_material / pool_water stays consistent.
            ! Water side: wetwat += gwat*dt with gwat including the
            ! snow-bottom outflow (snowwater's qout at layer 0). Without
            ! this companion term, pool_water is systematically low and
            ! pool_ratio is overstated, over-draining mt_loss for snowy
            ! wetland patches (non-split mode is the worst case — no
            ! atmospheric term enters pool_water otherwise).
            pool_water  = pool_water  + gwat_snow_local
         ELSE
            pool_material = pool_material + mt_pg_rain_ground(imt, ipatch) &
                                      + mt_sm_carry(imt, ipatch)
         ENDIF

         ! Atmospheric input material (dew + frost)
         IF (q_dew_in + q_frost_in > mt_tiny) THEN
            mt_dew_input = q_dew_in * deposition_ratio_for(layer_temp(1), .false.)
            mt_frost_input = q_frost_in * deposition_ratio_for(layer_temp(1), .true.)
            pool_material = pool_material + mt_dew_input + mt_frost_input
            a_mt_precip(imt, ipatch) = a_mt_precip(imt, ipatch) + mt_dew_input + mt_frost_input
            a_water_precip(imt, ipatch) = a_water_precip(imt, ipatch) + q_dew_in + q_frost_in
         ENDIF

         ! Isotopes: vapour gained through negative transpiration, which the
         ! host wetland water includes (see MOD_MT_Isotope_SoilFlow).
         IF (material_is_isotope(imt)) THEN
            q_vapor_uptake = isotope_wetland_vapor_uptake(etr, deltim)
            IF (q_vapor_uptake > mt_tiny) THEN
               mt_flux = q_vapor_uptake * deposition_ratio_for(layer_temp(1), .false.)
               pool_water  = pool_water  + q_vapor_uptake
               pool_material = pool_material + mt_flux
               a_mt_precip(imt, ipatch) = a_mt_precip(imt, ipatch) + mt_flux
               a_water_precip(imt, ipatch) = a_water_precip(imt, ipatch) + q_vapor_uptake
            ENDIF
         ENDIF

         IF (material_is_nonvolatile_solute(imt) .and. &
             (pool_water > mt_water_min_for_ratio .or. &
              mt_has_dissolved_limit(imt)) .and. &
             mt_surface_residue(imt, ipatch) > mt_tiny) THEN
            pool_material = pool_material + mt_surface_residue(imt, ipatch)
            mt_surface_residue(imt, ipatch) = 0._r8
         ENDIF
         CALL mt_equilibrate_dissolved(imt, pool_water, pool_material, &
            mt_surface_solid(imt, ipatch))

         ! The mixed pool already includes the fixed aquifer reference once.
         ! The guard takes relative (wa, mt_wa), so remove that reference
         ! before validating rather than counting it twice.
         CALL check_isotope_aquifer(imt, ipatch, pool_water-aquifer_ref_water, &
            pool_material-aquifer_ref_mass, 'wetland mixed pool')

         !--------------------------------------------------------
         ! 4) Remove evap/subl at atmospheric-exchange signatures and
         !    transpiration at mixed pool ratio.
         !    mt_loss is UNCAPPED — over-drawn losses push pool_material
         !    negative, which becomes material debt assigned to mt_wa at
         !    step 5 (matching the water-side debt wa<0). Capping here
         !    would break conservation across the deficit→recovery cycle.
         !--------------------------------------------------------
         loss_water = q_evap_out + q_subl_out + q_etr_out
         IF (abs(pool_water) > mt_water_min_for_ratio) THEN
            pool_ratio = pool_material / pool_water
         ELSE
            pool_ratio = R_atm
         ENDIF

         fractionate_pool_loss = material_fractionation_active(imt) .and. &
            pool_water > mt_water_min_for_ratio .and. pool_material > mt_tiny
         IF (fractionate_pool_loss) THEN
            ! Apply isotope fractionation to wetland open-water/ice
            ! atmospheric losses.  Do this sequentially against a
            ! temporary finite pool so large single-step losses integrate
            ! the enrichment of the residual water.  Transpiration remains
            ! at the mixed-pool source ratio.
            pool_water_loss  = pool_water
            pool_mt_loss = pool_material

            loss_liq_avail = min(q_evap_out, max(pool_water_loss, 0._r8))
            mt_evap_loss = atmospheric_loss_material(pool_mt_loss, pool_water_loss, &
               loss_liq_avail, layer_temp(1), .false.)
            IF (q_evap_out > loss_liq_avail) THEN
               mt_evap_loss = mt_evap_loss + (q_evap_out - loss_liq_avail) * pool_ratio
            ENDIF
            pool_mt_loss = pool_mt_loss - mt_evap_loss
            pool_water_loss  = pool_water_loss  - q_evap_out

            IF (pool_water_loss > mt_water_min_for_ratio) THEN
               pool_ratio_loss = pool_mt_loss / pool_water_loss
            ELSE
               pool_ratio_loss = pool_ratio
            ENDIF
            loss_ice_avail = min(q_subl_out, max(pool_water_loss, 0._r8))
            mt_subl_loss = atmospheric_loss_material(pool_mt_loss, pool_water_loss, &
               loss_ice_avail, layer_temp(1), .true.)
            IF (q_subl_out > loss_ice_avail) THEN
               mt_subl_loss = mt_subl_loss + (q_subl_out - loss_ice_avail) * pool_ratio_loss
            ENDIF
            pool_mt_loss = pool_mt_loss - mt_subl_loss
            pool_water_loss  = pool_water_loss  - q_subl_out

            IF (pool_water_loss > mt_water_min_for_ratio) THEN
               pool_ratio_loss = pool_mt_loss / pool_water_loss
            ELSE
               pool_ratio_loss = pool_ratio
            ENDIF
            transp_source_material = q_etr_out * pool_ratio_loss
            transp_output_material = transp_source_material
            transp_frac_active = q_etr_out > mt_tiny .and. &
               present(tleaf_frac) .and. present(forc_q_frac) .and. present(forc_psrf_frac) .and. &
               present(lai_frac) .and. present(rst_frac) .and. allocated(mt_leaf_delta_e)
            IF (transp_frac_active) THEN
               tleaf_use = tleaf_frac
               leaf_area_use = lai_frac
               rst_use = rst_frac
               ra_use = 0._r8
               IF (present(ra_frac)) ra_use = max(ra_frac, 0._r8)
               IF (material_forcing_has_vapor(imt, ipatch)) THEN
                  R_vapor = R_atm
               ELSE
                  R_vapor = pool_ratio_loss / max(material_alpha_liq_vap(imt, tleaf_use), mt_tiny)
               ENDIF
               relhum_leaf = material_surface_relhum(forc_q_frac, forc_psrf_frac, tleaf_use, .false.)
               CALL material_transpiration_nss_ratio(imt, pool_ratio_loss, R_vapor, &
                  tleaf_use, relhum_leaf, forc_psrf_frac, q_etr_out, deltim, leaf_area_use, ra_use, rst_use, &
                  mt_leaf_delta_e(imt, ipatch), mt_leaf_delta_b(imt, ipatch), &
                  mt_leaf_peclet(imt, ipatch), mt_leaf_water_moles(imt, ipatch), &
                  transp_ratio, leaf_delta_e_new, leaf_delta_b_new, leaf_peclet_new, leaf_moles_new)
               mt_leaf_delta_e(imt, ipatch) = leaf_delta_e_new
               mt_leaf_delta_b(imt, ipatch) = leaf_delta_b_new
               mt_leaf_peclet(imt, ipatch) = leaf_peclet_new
               mt_leaf_water_moles(imt, ipatch) = leaf_moles_new
               transp_output_material = q_etr_out * transp_ratio
               IF (allocated(mt_leaf_iso_storage)) THEN
                  mt_leaf_iso_storage(imt, ipatch) = mt_leaf_iso_storage(imt, ipatch) &
                     + transp_source_material - transp_output_material
               ENDIF
            ENDIF
            mt_etr_loss = transp_source_material
            mt_loss = mt_evap_loss + mt_subl_loss + mt_etr_loss

         ELSEIF (material_is_nonvolatile_solute(imt)) THEN
            mt_loss = 0._r8
            mt_evap_loss = 0._r8
            mt_subl_loss = 0._r8
            mt_etr_loss  = 0._r8
            transp_output_material = 0._r8
         ELSE
            mt_loss = loss_water * pool_ratio
            mt_evap_loss = q_evap_out * pool_ratio
            mt_subl_loss = q_subl_out * pool_ratio
            mt_etr_loss  = q_etr_out  * pool_ratio
            transp_output_material = mt_etr_loss
         ENDIF
         ! Book signed components independently. Atmospheric isotope uptake
         ! can cancel a positive transpiration source in mt_loss; gating both
         ! diagnostics on that net would hide real fluxes and break the NSS
         ! storage/output ledger.
         IF (q_evap_out + q_subl_out > mt_tiny .or. &
             abs(mt_evap_loss + mt_subl_loss) > mt_tiny) THEN
            CALL material_book_evap_loss(imt, ipatch, mt_evap_loss + mt_subl_loss, &
               q_evap_out + q_subl_out, MT_EVAP_KIND_WETLAND)
         ENDIF
         IF (q_etr_out > mt_tiny .or. abs(transp_output_material) > mt_tiny) THEN
            CALL material_book_evap_loss(imt, ipatch, transp_output_material, q_etr_out, &
               MT_EVAP_KIND_TRANSP)
            a_mt_transp_src(imt, ipatch) = a_mt_transp_src(imt, ipatch) + mt_etr_loss
         ENDIF
         pool_material = pool_material - mt_loss
         pool_water = pool_water - loss_water
         CALL mt_equilibrate_dissolved(imt, pool_water, pool_material, &
            mt_surface_solid(imt, ipatch))
         CALL check_isotope_aquifer(imt, ipatch, pool_water-aquifer_ref_water, &
            pool_material-aquifer_ref_mass, 'wetland after loss')

         !--------------------------------------------------------
         ! 5) Redistribute to post-WATER wdsrf / wa / wetwat / rsur.
         !    Water balance holds in signed form:
         !      pool_water_post_loss = wetwat + wdsrf + wa + rsur*dt
         !    with wa<0 in the deficit branch. wetwat and wdsrf are non-
         !    negative by WATER_VSF construction, so only wa carries the
         !    debt sign. pool_ratio is the mixed ratio of the post-loss
         !    pool and applies uniformly to all three post-WATER pools.
         !--------------------------------------------------------
         mt_rsur_local = 0._r8
         IF (pool_material < 0._r8) THEN
            ! Deficit branch: negative material debt belongs on the signed
            ! aquifer pool only. Do not spread debt into non-negative wetland,
            ! surface-water, or runoff pools.
            mt_wetwat(imt, ipatch) = 0._r8
            mt_wdsrf (imt, ipatch) = 0._r8
            mt_rsur_local = 0._r8
            mt_wa (imt, ipatch) = pool_material - mt_rsur_local
         ELSE
            IF (pool_water > mt_water_min_for_ratio) THEN
               ! Wetland pool R is bounded at source by the evaporative-loss
               ! cap (material_evaporative_mt_loss r_max), so pool_ratio can
               ! no longer run away near drydown.  Use the physical
               ! denominator floor here as a second guard before exporting to
               ! wetland runoff or redistributing storage.
               pool_ratio = pool_material / pool_water
            ELSEIF (pool_material <= mt_tiny) THEN
               pool_ratio = 0._r8
            ELSE
               ! Positive material with no physically resolved post-loss
               ! wetland water would reappear as an unbounded concentration
               ! once routed to rsur/wdsrf/wetwat.  Keep the accounting
               ! explicit instead: this residual is reported by the
               ! conservation diagnostics and is subject to the independent
               ! residual hard cap.
               IF (material_is_nonvolatile_solute(imt)) THEN
                  mt_surface_residue(imt, ipatch) = mt_surface_residue(imt, ipatch) + &
                     pool_material
               ELSE
                  IF (allocated(mt_numerical_residual_step)) THEN
                     mt_numerical_residual_step(imt, ipatch) = &
                        mt_numerical_residual_step(imt, ipatch) - pool_material
                  ENDIF
               ENDIF
               pool_material = 0._r8
               pool_ratio = 0._r8
            ENDIF

            mt_wetwat(imt, ipatch) = wetwat * pool_ratio
            mt_wdsrf (imt, ipatch) = wdsrf  * pool_ratio
            mt_wa    (imt, ipatch) = &
               material_aquifer_actual_water(wa, aquifer_ref_water) * pool_ratio - aquifer_ref_mass

            IF (rsur > mt_tiny) THEN
               mt_rsur_local = rsur * deltim * pool_ratio
               a_mt_rsur   (imt, ipatch) = a_mt_rsur   (imt, ipatch) + mt_rsur_local
               a_mt_rnof   (imt, ipatch) = a_mt_rnof   (imt, ipatch) + mt_rsur_local
               mt_rnof_step(imt, ipatch) = mt_rnof_step(imt, ipatch) + mt_rsur_local
            ENDIF
         ENDIF

         ! The wetland water solver can leave a tiny signed mismatch between
         ! pool_water and the post-WATER destinations. Keep that internal
         ! redistribution residual in storage, not in the material budget. This
         ! applies to all material modes; dynamic forcing invalidates the old
         ! no-fractionation shortcut that reset wa to the current vapor ratio.
         redist_storage = mt_wetwat(imt, ipatch) + mt_wdsrf(imt, ipatch) &
                        + mt_wa(imt, ipatch) + aquifer_ref_mass
         redist_target = pool_material - mt_rsur_local
         redist_resid = redist_storage - redist_target
         IF (abs(redist_resid) > mt_tiny) THEN
            IF (redist_resid > 0._r8) THEN
               redist_fix = min(redist_resid, mt_wetwat(imt, ipatch))
               mt_wetwat(imt, ipatch) = mt_wetwat(imt, ipatch) - redist_fix
               redist_resid = redist_resid - redist_fix

               redist_fix = min(redist_resid, mt_wdsrf(imt, ipatch))
               mt_wdsrf(imt, ipatch) = mt_wdsrf(imt, ipatch) - redist_fix
               redist_resid = redist_resid - redist_fix

               mt_wa(imt, ipatch) = mt_wa(imt, ipatch) - redist_resid
            ELSE
               IF (wetwat > mt_tiny) THEN
                  mt_wetwat(imt, ipatch) = mt_wetwat(imt, ipatch) - redist_resid
               ELSEIF (wdsrf > mt_tiny) THEN
                  mt_wdsrf(imt, ipatch) = mt_wdsrf(imt, ipatch) - redist_resid
               ELSE
                  mt_wa(imt, ipatch) = mt_wa(imt, ipatch) - redist_resid
               ENDIF
            ENDIF
         ENDIF
         CALL check_isotope_aquifer(imt, ipatch, wa, mt_wa(imt, ipatch), 'wetland end')
         ! Recover retained nonvolatile solute when its ice carrier is exhausted.
         DO j = lb, nl_soil
            IF (material_is_nonvolatile_solute(imt) .and. &
                mt_wice_soisno(imt, j, ipatch) > mt_tiny) THEN
               IF (max(wice_soisno(j), 0._r8) <= mt_water_min_for_ratio) THEN
                  IF (max(wliq_soisno(j), 0._r8) > mt_water_min_for_ratio) THEN
                     mt_wliq_soisno(imt, j, ipatch) = &
                        mt_wliq_soisno(imt, j, ipatch) + mt_wice_soisno(imt, j, ipatch)
                  ELSE
                     mt_solid_soisno(imt, j, ipatch) = &
                        mt_solid_soisno(imt, j, ipatch) + mt_wice_soisno(imt, j, ipatch)
                  ENDIF
                  mt_wice_soisno(imt, j, ipatch) = 0._r8
               ENDIF
            ENDIF
            CALL mt_equilibrate_dissolved(imt, max(wliq_soisno(j), 0._r8), &
               mt_wliq_soisno(imt, j, ipatch), mt_solid_soisno(imt, j, ipatch))
         ENDDO
         IF (allocated(mt_waterstorage) .and. present(waterstorage_patch)) THEN
            CALL mt_equilibrate_dissolved(imt, &
               max(waterstorage_patch - max(qflx_irrig_ground, 0._r8) * deltim, 0._r8), &
               mt_waterstorage(imt, ipatch), mt_waterstorage_solid(imt, ipatch))
         ENDIF
      ENDDO

   CONTAINS

      real(r8) FUNCTION layer_temp (jlay)
         integer, intent(in) :: jlay

         layer_temp = 273.15_r8
         IF (jlay >= lbound(t_soisno, 1) .and. jlay <= ubound(t_soisno, 1)) THEN
            layer_temp = t_soisno(jlay)
         ENDIF
      END FUNCTION layer_temp

      real(r8) FUNCTION atmospheric_loss_material (pool_mt, pool_water, water_loss, temp_k, from_ice)
         real(r8), intent(in) :: pool_mt
         real(r8), intent(in) :: pool_water
         real(r8), intent(in) :: water_loss
         real(r8), intent(in) :: temp_k
         logical,  intent(in) :: from_ice

         ! from_ice routes through the surface-skin limiter: a snow layer is
         ! not a well-mixed reservoir for sublimation.
         atmospheric_loss_material = material_atmospheric_mt_loss(pool_mt, pool_water, &
            water_loss, temp_k, from_ice, evap_affine_for, mt_tiny, &
            merge(materials(imt)%ref_ratio * (1._r8 + mt_delta_sanity_max / 1000._r8), &
                  0._r8, material_fractionation_active(imt)), material_is_nonvolatile_solute(imt), &
            skin_mass = DEF_ISO_SUBL_SKIN_MM)
      END FUNCTION atmospheric_loss_material

      SUBROUTINE evap_affine_for (temp_k, from_ice, slope, offset)
         real(r8), intent(in) :: temp_k
         logical,  intent(in) :: from_ice
         real(r8), intent(out) :: slope, offset
         real(r8) :: relhum, alpha_k

         slope  = 1._r8
         offset = 0._r8
         IF (.not. material_fractionation_active(imt)) RETURN
         IF (.not. present(forc_q_frac) .or. .not. present(forc_psrf_frac)) RETURN

         relhum = material_surface_relhum(forc_q_frac, forc_psrf_frac, temp_k, from_ice)
         IF (from_ice) THEN
            alpha_k = material_alpha_kinetic_craig_gordon(imt, .true.)
         ELSE
            alpha_k = material_alpha_kinetic_open_water(imt, &
               sqrt(max(forc_us_frac*forc_us_frac + forc_vs_frac*forc_vs_frac, 0._r8)))
         ENDIF
         CALL material_craig_gordon_evap_affine(imt, R_atm, temp_k, relhum, &
            alpha_k, from_ice, slope, offset)
      END SUBROUTINE evap_affine_for

         real(r8) FUNCTION deposition_ratio_for (temp_k, from_ice)
            real(r8), intent(in) :: temp_k
            logical,  intent(in) :: from_ice

            deposition_ratio_for = R_atm
            ! Nonvolatile solutes do not enter through wetland
            ! dew/frost water-vapour deposition.  Keep wet deposition as a
            ! separate explicit source rather than implicit vapour uptake.
            IF (material_is_nonvolatile_solute(imt)) THEN
               deposition_ratio_for = 0._r8
               RETURN
            ENDIF
            IF (.not. material_fractionation_active(imt)) RETURN
            deposition_ratio_for = material_equilibrium_deposition_ratio(imt, R_atm, temp_k, from_ice)
         END FUNCTION deposition_ratio_for

   END SUBROUTINE mt_wetland

   SUBROUTINE exhaust_surface_phase (imt, ipatch, phase_material, water_loss, evap_kind)
      integer,  intent(in)    :: imt, ipatch, evap_kind
      real(r8), intent(inout) :: phase_material
      real(r8), intent(in)    :: water_loss
      real(r8) :: material_loss

      material_loss = max(phase_material, 0._r8)
      IF (material_loss <= mt_tiny) RETURN

      IF (material_is_nonvolatile_solute(imt)) THEN
         ! Sublimation/evaporation removes carrier water, not a nonvolatile
         ! solute.  Retain its mass in this exact snow/ice/liquid layer so a
         ! later melt or rewet event releases it from the correct reservoir.
         RETURN
      ENDIF

      phase_material = 0._r8
      CALL material_book_evap_loss(imt, ipatch, material_loss, water_loss, evap_kind)
   END SUBROUTINE exhaust_surface_phase


END MODULE MOD_MT_SoilWater
#endif
