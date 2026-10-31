#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_Vars

   USE MOD_Precision
   USE MOD_MT_Namelist, only: DEF_MT_LULCC_ABORT_NBAD
   USE MOD_SPMD_Task, only: p_is_io, CoLM_stop
   USE MOD_MT_Defs, only: nmts, materials, material_build_descriptor_identity, &
      mt_uses_generic_water, material_is_isotope, material_aquifer_actual_water
   USE MOD_MT_Solute, only: mt_equilibrate_dissolved

   IMPLICIT NONE
   SAVE
   PRIVATE

   logical :: lulcc_area_fallback_warned = .false.
   real(r8), parameter :: lulcc_mass_abs_tol = 1.e-8_r8
   real(r8), parameter :: lulcc_mass_rel_tol = 1.e-10_r8

   real(r8), allocatable :: mt_ldew_rain  (:,:)
   real(r8), allocatable :: mt_ldew_snow  (:,:)
   real(r8), allocatable :: mt_wliq_soisno(:,:,:)
   real(r8), allocatable :: mt_wice_soisno(:,:,:)
   real(r8), allocatable :: mt_wa         (:,:)
   ! Fixed isotope-only reference carrier. mt_wa remains a relative mass
   ! so its constant reference cancels from existing conservation budgets.
   real(r8), allocatable :: mt_aquifer_ref_water(:)
   real(r8), allocatable :: mt_aquifer_ref_mass(:,:)
   real(r8), allocatable :: mt_wdsrf      (:,:)
   real(r8), allocatable :: mt_wetwat     (:,:)
   ! Waterless surface residue for nonvolatile solutes.
   ! This is a numerical immobile store, not a precipitated mineral phase.
   real(r8), allocatable :: mt_surface_residue(:,:)
   ! Positive solute orphaned when aquifer water vanishes. It redissolves only
   ! when the aquifer returns, preventing first-rain release at the surface.
   real(r8), allocatable :: mt_subsurface_residue(:,:)
   ! Immobile precipitated solute. Layer, surface, aquifer, and irrigation
   ! reservoir inventories stay separate so saturation cannot move mass between
   ! hydrologically disconnected pools.
   real(r8), allocatable :: mt_solid_soisno(:,:,:)
   real(r8), allocatable :: mt_canopy_solid(:,:)
   real(r8), allocatable :: mt_surface_solid(:,:)
   real(r8), allocatable :: mt_subsurface_solid(:,:)
   real(r8), allocatable :: mt_waterstorage_solid(:,:)

   ! Material mirror of MOD_Vars_TimeVariables::waterstorage. Irrigation
   ! (drip/flood/paddy/sprinkler) is an internal transfer from the
   ! per-patch waterstorage reservoir to soil / canopy; water-side
   ! balance (CoLMMAIN:806) treats it that way, so the atmospheric-input
   ! term (forc_prc+forc_prl) does NOT include irrigation. Without a
   ! matching material reservoir the material side would have to record
   ! irrigation as external input, creating a systematic offset between
   ! the water and material conservation accountings. mt_waterstorage
   ! carries the isotope composition of the reservoir and is debited
   ! by exactly the same volume the water side withdraws each step.
   real(r8), allocatable :: mt_waterstorage(:,:)     ! (nmts, numpatch)

   ! Snow water equivalent material for pre-layer accumulation [kg/m2]
   ! Tracks material in scv before a snow layer is created (snl=0).
   ! When the first snow layer forms, mt_scv is transferred to mt_wice_soisno.
   real(r8), allocatable :: mt_scv        (:,:)  ! (nmts, numpatch)

   ! Throughfall material reaching ground [kg/m2 per step]
   ! Rain and snow kept separate: rain enters the surface mixed pool via
   ! mt_soil_water; snow enters mt_scv via mt_newsnow (mirrors
   ! water side where gwat only carries pg_rain). A previous "mt_pg_to_ground
   ! = rain+snow" combined var double-counted snow material when used in the
   ! surface pool — removed to prevent reintroduction.
   real(r8), allocatable :: mt_pg_rain_ground(:,:)    ! (nmts, numpatch)
   real(r8), allocatable :: mt_pg_snow_ground(:,:)    ! (nmts, numpatch)
   real(r8), allocatable :: mt_rnof_step(:,:)         ! (nmts, numpatch), accumulated runoff mass of
                                                       ! THIS MAIN STEP for the river hand-off. Reset once per
                                                       ! main step (CoLM.F90, before CoLMDRIVER); every lake/
                                                       ! glacier/waterbody substep adds to it. Transient state,
                                                       ! not in restart.

   ! Material carried by snowmelt sm from thin-snow (snl==0) scv during THERMAL.
   ! Set in CoLMMAIN after THERMAL based on the actual mt_scv reduction;
   ! consumed by mt_soil_water when snl==0 .and. sm>0. Replaces the old
   ! "sm * R_precip" injection which used a fixed init delta instead of the
   ! real scv ratio.
   real(r8), allocatable :: mt_sm_carry(:,:)          ! (nmts, numpatch)

   ! Non-steady-state leaf-water isotope memory for transpiration.
   ! These are diagnostic state variables for the isotope fractionation
   ! calculation; hydrology does not carry a prognostic leaf-tissue water pool.
   real(r8), allocatable :: mt_leaf_delta_e(:,:)      ! evaporation-site leaf water delta [permil]
   real(r8), allocatable :: mt_leaf_delta_b(:,:)      ! bulk leaf water delta [permil]
   real(r8), allocatable :: mt_leaf_peclet (:,:)      ! previous Peclet mixing factor [-]
   real(r8), allocatable :: mt_leaf_water_moles(:,:)  ! previous leaf water content [mol/m2]
   real(r8), allocatable :: mt_leaf_iso_storage(:,:)  ! NSS isotope storage anomaly [R*mm]
   logical,  allocatable :: mt_runtime_forced(:)      ! material has runtime atmospheric forcing

   ! LULCC snapshots for prognostic land-water material pools. These are
   ! saved before LULCC rebuilds landpatch and remapped after the new patch
   ! layout exists, mirroring the methane LULCC path.
   logical :: land_mt_lulcc_snapshot_valid = .false.
   real(r8), allocatable :: lulcc_mt_ldew_rain_old(:,:)
   real(r8), allocatable :: lulcc_mt_ldew_snow_old(:,:)
   real(r8), allocatable :: lulcc_mt_wliq_soisno_old(:,:,:)
   real(r8), allocatable :: lulcc_mt_wice_soisno_old(:,:,:)
   real(r8), allocatable :: lulcc_mt_wa_old(:,:)
   real(r8), allocatable :: lulcc_mt_aquifer_ref_water_old(:,:)
   real(r8), allocatable :: lulcc_mt_aquifer_ref_mass_old(:,:)
   real(r8), allocatable :: lulcc_mt_wdsrf_old(:,:)
   real(r8), allocatable :: lulcc_mt_wetwat_old(:,:)
   real(r8), allocatable :: lulcc_mt_surface_residue_old(:,:)
   real(r8), allocatable :: lulcc_mt_subsurface_residue_old(:,:)
   real(r8), allocatable :: lulcc_mt_solid_soisno_old(:,:,:)
   real(r8), allocatable :: lulcc_mt_canopy_solid_old(:,:)
   real(r8), allocatable :: lulcc_mt_surface_solid_old(:,:)
   real(r8), allocatable :: lulcc_mt_subsurface_solid_old(:,:)
   real(r8), allocatable :: lulcc_mt_waterstorage_solid_old(:,:)
   real(r8), allocatable :: lulcc_mt_waterstorage_old(:,:)
   real(r8), allocatable :: lulcc_mt_scv_old(:,:)
   real(r8), allocatable :: lulcc_mt_leaf_delta_e_old(:,:)
   real(r8), allocatable :: lulcc_mt_leaf_delta_b_old(:,:)
   real(r8), allocatable :: lulcc_mt_leaf_peclet_old(:,:)
   real(r8), allocatable :: lulcc_mt_leaf_water_moles_old(:,:)
   real(r8), allocatable :: lulcc_mt_leaf_iso_storage_old(:,:)

   real(r8), allocatable :: a_mt_precip   (:,:)
   ! Signed material mass exchanged with atmospheric vapour with NO net water
   ! flux (two-way equilibrium exchange on wet surfaces).  It is a genuine
   ! boundary flux, so it must be tracked separately from a_mt_precip: the
   ! latter is paired with a water flux and is cross-checked against
   ! water_input * R_init on fixed-signature materials.
   real(r8), allocatable :: a_mt_vapor_exchange(:,:)
   real(r8), allocatable :: a_water_precip (:,:)
	   integer, parameter :: MT_EVAP_KIND_TOTAL       = 0
	   integer, parameter :: MT_EVAP_KIND_TRANSP      = 1
	   integer, parameter :: MT_EVAP_KIND_SOILEVAP    = 2
	   integer, parameter :: MT_EVAP_KIND_CANOPYEVAP  = 3
	   integer, parameter :: MT_EVAP_KIND_SUBL        = 4
	   integer, parameter :: MT_EVAP_KIND_WETLAND     = 5
	   real(r8), allocatable :: a_mt_evap     (:,:)
	   real(r8), allocatable :: a_water_evap_gross(:,:)
	   real(r8), allocatable :: a_mt_transp   (:,:)
	   real(r8), allocatable :: a_mt_transp_src(:,:)
	   real(r8), allocatable :: a_water_transp (:,:)
	   real(r8), allocatable :: a_mt_soilevap (:,:)
	   real(r8), allocatable :: a_water_soilevap(:,:)
	   real(r8), allocatable :: a_mt_canopyevap(:,:)
	   real(r8), allocatable :: a_water_canopyevap(:,:)
	   real(r8), allocatable :: a_mt_subl(:,:)
	   real(r8), allocatable :: a_water_subl(:,:)
	   real(r8), allocatable :: a_mt_wetland_evap(:,:)
	   real(r8), allocatable :: a_water_wetland_evap(:,:)
   real(r8), allocatable :: a_mt_rsur     (:,:)
   real(r8), allocatable :: a_mt_rsub     (:,:)
   ! External nitrogen input from BGC leaching, booked separately so the
   ! balance check can count it as an input and the fixed-signature check can
   ! exclude it (it is not precipitation).  kg N m-2 per step.
   real(r8), allocatable :: a_mt_bgc_leach(:,:)
   real(r8), allocatable :: a_mt_rnof     (:,:)
   real(r8), allocatable :: a_water_rnof   (:,:)
   real(r8), allocatable :: a_mt_qinfl    (:,:)
   real(r8), allocatable :: a_mt_qcharge  (:,:)

   real(r8), allocatable :: mt_storage_beg(:,:)
   real(r8), allocatable :: mt_balance_err(:,:)
	   ! Net decay source/sink applied during the current land step.
	   ! Positive adds material to storage; negative removes material.
	   real(r8), allocatable :: mt_decay_source_step(:,:)
	   ! Explicit numerical source/sink used only to reconcile solver residuals
	   ! that have no physical flux signature. Positive adds material to storage.
	   real(r8), allocatable :: mt_numerical_residual_step(:,:)
	   real(r8), allocatable :: mt_numerical_water_step(:,:)

   real(r8), allocatable :: a_mt_ldew_mass (:,:)
   real(r8), allocatable :: a_water_ldew    (:)
   real(r8), allocatable :: a_mt_soil_mass (:,:,:)  ! (nmts, nl_soil, numpatch)
   real(r8), allocatable :: a_water_soil    (:,:)     ! (nl_soil, numpatch)
   real(r8), allocatable :: a_mt_snow_mass (:,:,:)  ! (nmts, abs(maxsnl), numpatch) snow layers
   real(r8), allocatable :: a_water_snow    (:,:)     ! (abs(maxsnl), numpatch) snow layers

   ! Surface / aquifer / wetland / thin-snow diagnostic accumulators.
   ! mt_balance_check already includes these four pools in the conservation
   ! sum, but the history writer only had canopy + soil + snow, so users could
   ! not see the δ values of wa (aquifer), wdsrf (surface water), wetwat
   ! (wetland pool) or scv (pre-layer thin snow). These accumulators let
   ! MOD_Hist emit a per-pool material ratio diagnostic alongside the existing
   ! ldew / soisno outputs.
	   real(r8), allocatable :: a_mt_wa_mass    (:,:)   ! (nmts, numpatch)
	   real(r8), allocatable :: a_water_wa       (:)     ! (numpatch)
	   real(r8), allocatable :: a_mt_wa_debt_mass(:,:)  ! signed-debt magnitude diagnostics
	   real(r8), allocatable :: a_water_wa_debt   (:)    ! aquifer debt magnitude
   real(r8), allocatable :: a_water_aquifer_actual(:)
   real(r8), allocatable :: a_mt_aquifer_actual_mass(:,:)
   real(r8), allocatable :: a_mt_wdsrf_mass (:,:)
   real(r8), allocatable :: a_water_wdsrf    (:)
   real(r8), allocatable :: a_mt_wetwat_mass(:,:)
   real(r8), allocatable :: a_water_wetwat   (:)
   real(r8), allocatable :: a_mt_surface_residue_mass(:,:)
   real(r8), allocatable :: a_mt_subsurface_residue_mass(:,:)
   real(r8), allocatable :: a_mt_layer_dry_mass(:,:)
   real(r8), allocatable :: a_mt_solid_mass(:,:)
   real(r8), allocatable :: a_mt_scv_mass   (:,:)
   real(r8), allocatable :: a_water_scv      (:)

   ! External state contract.  LULCC snapshots and remapping work arrays stay
   ! private; only pools and diagnostics consumed by the land physics, restart,
   ! history, and host facades are exposed.
   PUBLIC :: allocate_Material_Vars, deallocate_Material_Vars, flush_Material_Acc
   PUBLIC :: save_land_mt_lulcc_state, remap_land_mt_lulcc_state
   PUBLIC :: mt_patch_land_water_mass
   PUBLIC :: zero_provider_owned_land_mt_state, sync_mt_patch_ratio
   PUBLIC :: material_book_evap_loss
   PUBLIC :: mt_add_rnof_step, mt_reset_rnof_step

   PUBLIC :: mt_ldew_rain, mt_ldew_snow
   PUBLIC :: mt_wliq_soisno, mt_wice_soisno
   PUBLIC :: mt_wa, mt_wdsrf, mt_wetwat
   PUBLIC :: mt_aquifer_ref_water, mt_aquifer_ref_mass
   PUBLIC :: mt_surface_residue, mt_subsurface_residue, mt_waterstorage
   PUBLIC :: mt_solid_soisno, mt_surface_solid, mt_subsurface_solid
   PUBLIC :: mt_canopy_solid, mt_waterstorage_solid
   PUBLIC :: mt_scv, mt_pg_rain_ground, mt_pg_snow_ground, mt_rnof_step
   PUBLIC :: mt_sm_carry, mt_runtime_forced
   PUBLIC :: mt_leaf_delta_e, mt_leaf_delta_b, mt_leaf_peclet
   PUBLIC :: mt_leaf_water_moles, mt_leaf_iso_storage
   PUBLIC :: mt_storage_beg, mt_balance_err
   PUBLIC :: mt_decay_source_step, mt_numerical_residual_step, mt_numerical_water_step

   PUBLIC :: MT_EVAP_KIND_TRANSP, MT_EVAP_KIND_SOILEVAP
   PUBLIC :: MT_EVAP_KIND_CANOPYEVAP, MT_EVAP_KIND_SUBL, MT_EVAP_KIND_WETLAND
   PUBLIC :: a_mt_precip, a_water_precip, a_mt_evap, a_water_evap_gross
   PUBLIC :: a_mt_vapor_exchange
   PUBLIC :: a_mt_transp, a_mt_transp_src, a_water_transp
   PUBLIC :: a_mt_soilevap, a_water_soilevap
   PUBLIC :: a_mt_canopyevap, a_water_canopyevap
   PUBLIC :: a_mt_subl, a_water_subl, a_mt_wetland_evap, a_water_wetland_evap
   PUBLIC :: a_mt_rsur, a_mt_rsub, a_mt_rnof, a_water_rnof, a_mt_qinfl, a_mt_qcharge
   PUBLIC :: a_mt_bgc_leach
   PUBLIC :: a_mt_ldew_mass, a_water_ldew
   PUBLIC :: a_mt_soil_mass, a_water_soil, a_mt_snow_mass, a_water_snow
   PUBLIC :: a_mt_wa_mass, a_water_wa, a_mt_wa_debt_mass, a_water_wa_debt
   PUBLIC :: a_water_aquifer_actual, a_mt_aquifer_actual_mass
   PUBLIC :: a_mt_wdsrf_mass, a_water_wdsrf
   PUBLIC :: a_mt_wetwat_mass, a_water_wetwat
   PUBLIC :: a_mt_surface_residue_mass, a_mt_subsurface_residue_mass
   PUBLIC :: a_mt_layer_dry_mass, a_mt_solid_mass, a_mt_scv_mass, a_water_scv

   ! Whether the land material state came from the restart file (not a
   ! cold start); only then do the history sums saved with it belong to the
   ! current species (MOD_HistAccRestart).
   logical, PUBLIC :: mt_land_restart_loaded = .false.

CONTAINS

   SUBROUTINE allocate_Material_Vars (numpatch, maxsnl, nl_soil)
      IMPLICIT NONE
      integer, intent(in) :: numpatch, maxsnl, nl_soil
      integer :: imt
      ! History writers slice material rows on every rank, including the master
      ! and empty workers. Keep their second dimension empty, not unallocated.
      IF (nmts <= 0) RETURN

      allocate(mt_ldew_rain   (nmts, numpatch));           mt_ldew_rain   = 0._r8
      allocate(mt_ldew_snow   (nmts, numpatch));           mt_ldew_snow   = 0._r8
      allocate(mt_wliq_soisno (nmts, maxsnl+1:nl_soil, numpatch)); mt_wliq_soisno = 0._r8
      allocate(mt_wice_soisno (nmts, maxsnl+1:nl_soil, numpatch)); mt_wice_soisno = 0._r8
      allocate(mt_wa          (nmts, numpatch));           mt_wa          = 0._r8
      allocate(mt_aquifer_ref_water(numpatch));                mt_aquifer_ref_water = 0._r8
      allocate(mt_aquifer_ref_mass(nmts, numpatch));       mt_aquifer_ref_mass = 0._r8
      allocate(mt_wdsrf       (nmts, numpatch));           mt_wdsrf       = 0._r8
      allocate(mt_wetwat      (nmts, numpatch));           mt_wetwat      = 0._r8
      allocate(mt_surface_residue(nmts, numpatch));        mt_surface_residue = 0._r8
      allocate(mt_subsurface_residue(nmts, numpatch));     mt_subsurface_residue = 0._r8
      allocate(mt_solid_soisno(nmts, maxsnl+1:nl_soil, numpatch)); mt_solid_soisno = 0._r8
      allocate(mt_canopy_solid(nmts, numpatch));          mt_canopy_solid = 0._r8
      allocate(mt_surface_solid(nmts, numpatch));          mt_surface_solid = 0._r8
      allocate(mt_subsurface_solid(nmts, numpatch));       mt_subsurface_solid = 0._r8
      allocate(mt_waterstorage_solid(nmts, numpatch));     mt_waterstorage_solid = 0._r8
      allocate(mt_scv         (nmts, numpatch));           mt_scv         = 0._r8
      allocate(mt_waterstorage(nmts, numpatch));           mt_waterstorage = 0._r8
      allocate(mt_pg_rain_ground(nmts, numpatch));      mt_pg_rain_ground = 0._r8
      allocate(mt_pg_snow_ground(nmts, numpatch));      mt_pg_snow_ground = 0._r8
      allocate(mt_rnof_step(nmts, numpatch));           mt_rnof_step = 0._r8
      allocate(mt_sm_carry (nmts, numpatch));           mt_sm_carry  = 0._r8
      allocate(mt_leaf_delta_e(nmts, numpatch));        mt_leaf_delta_e = 0._r8
      allocate(mt_leaf_delta_b(nmts, numpatch));        mt_leaf_delta_b = 0._r8
      allocate(mt_leaf_peclet (nmts, numpatch));        mt_leaf_peclet  = 1._r8
      allocate(mt_leaf_water_moles(nmts, numpatch));    mt_leaf_water_moles = 0._r8
      allocate(mt_leaf_iso_storage(nmts, numpatch));    mt_leaf_iso_storage = 0._r8
      allocate(mt_runtime_forced(nmts));                mt_runtime_forced = .false.
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         mt_leaf_delta_e(imt, :) = materials(imt)%init_delta
         mt_leaf_delta_b(imt, :) = materials(imt)%init_delta
      ENDDO

	      allocate(a_mt_precip    (nmts, numpatch));           a_mt_precip    = 0._r8
	      allocate(a_mt_vapor_exchange(nmts, numpatch));       a_mt_vapor_exchange = 0._r8
	      allocate(a_water_precip  (nmts, numpatch));           a_water_precip  = 0._r8
		      allocate(a_mt_evap      (nmts, numpatch));           a_mt_evap      = 0._r8
		      allocate(a_water_evap_gross(nmts, numpatch));         a_water_evap_gross = 0._r8
		      allocate(a_mt_transp    (nmts, numpatch));           a_mt_transp    = 0._r8
		      allocate(a_mt_transp_src(nmts, numpatch));           a_mt_transp_src = 0._r8
		      allocate(a_water_transp  (nmts, numpatch));           a_water_transp  = 0._r8
		      allocate(a_mt_soilevap  (nmts, numpatch));           a_mt_soilevap  = 0._r8
		      allocate(a_water_soilevap(nmts, numpatch));           a_water_soilevap = 0._r8
		      allocate(a_mt_canopyevap(nmts, numpatch));           a_mt_canopyevap = 0._r8
		      allocate(a_water_canopyevap(nmts, numpatch));         a_water_canopyevap = 0._r8
		      allocate(a_mt_subl      (nmts, numpatch));           a_mt_subl      = 0._r8
		      allocate(a_water_subl    (nmts, numpatch));           a_water_subl    = 0._r8
		      allocate(a_mt_wetland_evap(nmts, numpatch));         a_mt_wetland_evap = 0._r8
		      allocate(a_water_wetland_evap(nmts, numpatch));       a_water_wetland_evap = 0._r8
      allocate(a_mt_rsur      (nmts, numpatch));           a_mt_rsur      = 0._r8
      allocate(a_mt_rsub      (nmts, numpatch));           a_mt_rsub      = 0._r8
      allocate(a_mt_bgc_leach (nmts, numpatch));           a_mt_bgc_leach = 0._r8
      allocate(a_mt_rnof      (nmts, numpatch));           a_mt_rnof      = 0._r8
      allocate(a_water_rnof    (nmts, numpatch));           a_water_rnof    = 0._r8
      allocate(a_mt_qinfl     (nmts, numpatch));           a_mt_qinfl     = 0._r8
      allocate(a_mt_qcharge   (nmts, numpatch));           a_mt_qcharge   = 0._r8

	      allocate(mt_storage_beg (nmts, numpatch));           mt_storage_beg = 0._r8
	      allocate(mt_balance_err (nmts, numpatch));           mt_balance_err = 0._r8
	      allocate(mt_decay_source_step(nmts, numpatch));   mt_decay_source_step = 0._r8
	      allocate(mt_numerical_residual_step(nmts, numpatch)); mt_numerical_residual_step = 0._r8
	      allocate(mt_numerical_water_step(nmts, numpatch)); mt_numerical_water_step = 0._r8

	      allocate(a_mt_ldew_mass (nmts, numpatch));           a_mt_ldew_mass = 0._r8
      allocate(a_water_ldew    (numpatch));                     a_water_ldew    = 0._r8
      allocate(a_mt_soil_mass (nmts, nl_soil, numpatch));  a_mt_soil_mass = 0._r8
      allocate(a_water_soil    (nl_soil, numpatch));             a_water_soil    = 0._r8
      allocate(a_mt_snow_mass (nmts, abs(maxsnl), numpatch)); a_mt_snow_mass = 0._r8
      allocate(a_water_snow    (abs(maxsnl), numpatch));            a_water_snow    = 0._r8

	      allocate(a_mt_wa_mass    (nmts, numpatch));          a_mt_wa_mass     = 0._r8
	      allocate(a_water_wa       (numpatch));                    a_water_wa        = 0._r8
	      allocate(a_mt_wa_debt_mass(nmts, numpatch));         a_mt_wa_debt_mass = 0._r8
	      allocate(a_water_wa_debt   (numpatch));                   a_water_wa_debt    = 0._r8
      allocate(a_water_aquifer_actual(numpatch));              a_water_aquifer_actual = 0._r8
      allocate(a_mt_aquifer_actual_mass(nmts, numpatch)); a_mt_aquifer_actual_mass = 0._r8
      allocate(a_mt_wdsrf_mass (nmts, numpatch));          a_mt_wdsrf_mass  = 0._r8
      allocate(a_water_wdsrf    (numpatch));                    a_water_wdsrf     = 0._r8
      allocate(a_mt_wetwat_mass(nmts, numpatch));          a_mt_wetwat_mass = 0._r8
      allocate(a_water_wetwat   (numpatch));                    a_water_wetwat    = 0._r8
      allocate(a_mt_surface_residue_mass(nmts, numpatch)); a_mt_surface_residue_mass = 0._r8
      allocate(a_mt_subsurface_residue_mass(nmts, numpatch)); a_mt_subsurface_residue_mass = 0._r8
      allocate(a_mt_layer_dry_mass(nmts, numpatch));       a_mt_layer_dry_mass = 0._r8
      allocate(a_mt_solid_mass(nmts, numpatch));           a_mt_solid_mass = 0._r8
      allocate(a_mt_scv_mass   (nmts, numpatch));          a_mt_scv_mass    = 0._r8
      allocate(a_water_scv      (numpatch));                    a_water_scv       = 0._r8
   END SUBROUTINE allocate_Material_Vars

   SUBROUTINE deallocate_Material_Vars ()
      IMPLICIT NONE
      IF (allocated(mt_ldew_rain  )) deallocate(mt_ldew_rain  )
      IF (allocated(mt_ldew_snow  )) deallocate(mt_ldew_snow  )
      IF (allocated(mt_wliq_soisno)) deallocate(mt_wliq_soisno)
      IF (allocated(mt_wice_soisno)) deallocate(mt_wice_soisno)
      IF (allocated(mt_wa         )) deallocate(mt_wa         )
      IF (allocated(mt_aquifer_ref_water)) deallocate(mt_aquifer_ref_water)
      IF (allocated(mt_aquifer_ref_mass)) deallocate(mt_aquifer_ref_mass)
      IF (allocated(mt_wdsrf      )) deallocate(mt_wdsrf      )
      IF (allocated(mt_wetwat     )) deallocate(mt_wetwat     )
      IF (allocated(mt_surface_residue)) deallocate(mt_surface_residue)
      IF (allocated(mt_subsurface_residue)) deallocate(mt_subsurface_residue)
      IF (allocated(mt_solid_soisno)) deallocate(mt_solid_soisno)
      IF (allocated(mt_canopy_solid)) deallocate(mt_canopy_solid)
      IF (allocated(mt_surface_solid)) deallocate(mt_surface_solid)
      IF (allocated(mt_subsurface_solid)) deallocate(mt_subsurface_solid)
      IF (allocated(mt_waterstorage_solid)) deallocate(mt_waterstorage_solid)
      IF (allocated(mt_scv        )) deallocate(mt_scv        )
      IF (allocated(mt_waterstorage)) deallocate(mt_waterstorage)
      IF (allocated(mt_pg_rain_ground)) deallocate(mt_pg_rain_ground)
      IF (allocated(mt_pg_snow_ground)) deallocate(mt_pg_snow_ground)
      IF (allocated(mt_rnof_step    )) deallocate(mt_rnof_step    )
      IF (allocated(mt_sm_carry     )) deallocate(mt_sm_carry     )
      IF (allocated(mt_leaf_delta_e )) deallocate(mt_leaf_delta_e )
      IF (allocated(mt_leaf_delta_b )) deallocate(mt_leaf_delta_b )
      IF (allocated(mt_leaf_peclet  )) deallocate(mt_leaf_peclet  )
      IF (allocated(mt_leaf_water_moles)) deallocate(mt_leaf_water_moles)
      IF (allocated(mt_leaf_iso_storage)) deallocate(mt_leaf_iso_storage)
      IF (allocated(mt_runtime_forced)) deallocate(mt_runtime_forced)
		      IF (allocated(a_mt_precip   )) deallocate(a_mt_precip   )
		      IF (allocated(a_mt_vapor_exchange)) deallocate(a_mt_vapor_exchange)
		      IF (allocated(a_water_precip )) deallocate(a_water_precip )
		      IF (allocated(a_mt_evap     )) deallocate(a_mt_evap     )
		      IF (allocated(a_water_evap_gross)) deallocate(a_water_evap_gross)
		      IF (allocated(a_mt_transp   )) deallocate(a_mt_transp   )
		      IF (allocated(a_mt_transp_src)) deallocate(a_mt_transp_src)
		      IF (allocated(a_water_transp )) deallocate(a_water_transp )
		      IF (allocated(a_mt_soilevap )) deallocate(a_mt_soilevap )
		      IF (allocated(a_water_soilevap)) deallocate(a_water_soilevap)
		      IF (allocated(a_mt_canopyevap)) deallocate(a_mt_canopyevap)
		      IF (allocated(a_water_canopyevap)) deallocate(a_water_canopyevap)
		      IF (allocated(a_mt_subl     )) deallocate(a_mt_subl     )
		      IF (allocated(a_water_subl   )) deallocate(a_water_subl   )
		      IF (allocated(a_mt_wetland_evap)) deallocate(a_mt_wetland_evap)
		      IF (allocated(a_water_wetland_evap)) deallocate(a_water_wetland_evap)
      IF (allocated(a_mt_rsur     )) deallocate(a_mt_rsur     )
      IF (allocated(a_mt_rsub     )) deallocate(a_mt_rsub     )
      IF (allocated(a_mt_bgc_leach)) deallocate(a_mt_bgc_leach)
      IF (allocated(a_mt_rnof     )) deallocate(a_mt_rnof     )
      IF (allocated(a_water_rnof   )) deallocate(a_water_rnof   )
      IF (allocated(a_mt_qinfl    )) deallocate(a_mt_qinfl    )
      IF (allocated(a_mt_qcharge  )) deallocate(a_mt_qcharge  )
	      IF (allocated(mt_storage_beg)) deallocate(mt_storage_beg)
	      IF (allocated(mt_balance_err)) deallocate(mt_balance_err)
	      IF (allocated(mt_decay_source_step)) deallocate(mt_decay_source_step)
	      IF (allocated(mt_numerical_residual_step)) deallocate(mt_numerical_residual_step)
	      IF (allocated(mt_numerical_water_step)) deallocate(mt_numerical_water_step)
	      IF (allocated(a_mt_ldew_mass)) deallocate(a_mt_ldew_mass)
      IF (allocated(a_water_ldew   )) deallocate(a_water_ldew   )
      IF (allocated(a_mt_soil_mass)) deallocate(a_mt_soil_mass)
      IF (allocated(a_water_soil   )) deallocate(a_water_soil   )
      IF (allocated(a_mt_snow_mass)) deallocate(a_mt_snow_mass)
      IF (allocated(a_water_snow   )) deallocate(a_water_snow   )
	      IF (allocated(a_mt_wa_mass    )) deallocate(a_mt_wa_mass    )
	      IF (allocated(a_water_wa       )) deallocate(a_water_wa       )
	      IF (allocated(a_mt_wa_debt_mass)) deallocate(a_mt_wa_debt_mass)
	      IF (allocated(a_water_wa_debt   )) deallocate(a_water_wa_debt   )
      IF (allocated(a_water_aquifer_actual)) deallocate(a_water_aquifer_actual)
      IF (allocated(a_mt_aquifer_actual_mass)) deallocate(a_mt_aquifer_actual_mass)
      IF (allocated(a_mt_wdsrf_mass )) deallocate(a_mt_wdsrf_mass )
      IF (allocated(a_water_wdsrf    )) deallocate(a_water_wdsrf    )
      IF (allocated(a_mt_wetwat_mass)) deallocate(a_mt_wetwat_mass)
      IF (allocated(a_water_wetwat   )) deallocate(a_water_wetwat   )
      IF (allocated(a_mt_surface_residue_mass)) deallocate(a_mt_surface_residue_mass)
      IF (allocated(a_mt_subsurface_residue_mass)) deallocate(a_mt_subsurface_residue_mass)
      IF (allocated(a_mt_layer_dry_mass)) deallocate(a_mt_layer_dry_mass)
      IF (allocated(a_mt_solid_mass)) deallocate(a_mt_solid_mass)
      IF (allocated(a_mt_scv_mass   )) deallocate(a_mt_scv_mass   )
      IF (allocated(a_water_scv      )) deallocate(a_water_scv      )
   END SUBROUTINE deallocate_Material_Vars

   SUBROUTINE save_land_mt_lulcc_state ()
      IMPLICIT NONE

      IF (.not. allocated(mt_ldew_rain) .or. .not. allocated(mt_wliq_soisno)) THEN
         land_mt_lulcc_snapshot_valid = .false.
         RETURN
      ENDIF

      CALL clear_land_mt_lulcc_snapshot()

      allocate(lulcc_mt_ldew_rain_old(size(mt_ldew_rain,1), size(mt_ldew_rain,2)))
      lulcc_mt_ldew_rain_old = mt_ldew_rain
      allocate(lulcc_mt_ldew_snow_old(size(mt_ldew_snow,1), size(mt_ldew_snow,2)))
      lulcc_mt_ldew_snow_old = mt_ldew_snow
      allocate(lulcc_mt_wliq_soisno_old(size(mt_wliq_soisno,1), &
         lbound(mt_wliq_soisno,2):ubound(mt_wliq_soisno,2), size(mt_wliq_soisno,3)))
      lulcc_mt_wliq_soisno_old = mt_wliq_soisno
      allocate(lulcc_mt_wice_soisno_old(size(mt_wice_soisno,1), &
         lbound(mt_wice_soisno,2):ubound(mt_wice_soisno,2), size(mt_wice_soisno,3)))
      lulcc_mt_wice_soisno_old = mt_wice_soisno
      allocate(lulcc_mt_wa_old(size(mt_wa,1), size(mt_wa,2)))
      lulcc_mt_wa_old = mt_wa
      allocate(lulcc_mt_aquifer_ref_water_old(1, size(mt_aquifer_ref_water)))
      lulcc_mt_aquifer_ref_water_old(1, :) = mt_aquifer_ref_water
      allocate(lulcc_mt_aquifer_ref_mass_old(size(mt_aquifer_ref_mass,1), size(mt_aquifer_ref_mass,2)))
      lulcc_mt_aquifer_ref_mass_old = mt_aquifer_ref_mass
      allocate(lulcc_mt_wdsrf_old(size(mt_wdsrf,1), size(mt_wdsrf,2)))
      lulcc_mt_wdsrf_old = mt_wdsrf
      allocate(lulcc_mt_wetwat_old(size(mt_wetwat,1), size(mt_wetwat,2)))
      lulcc_mt_wetwat_old = mt_wetwat
      allocate(lulcc_mt_surface_residue_old(size(mt_surface_residue,1), size(mt_surface_residue,2)))
      lulcc_mt_surface_residue_old = mt_surface_residue
      allocate(lulcc_mt_subsurface_residue_old(size(mt_subsurface_residue,1), size(mt_subsurface_residue,2)))
      lulcc_mt_subsurface_residue_old = mt_subsurface_residue
      allocate(lulcc_mt_solid_soisno_old(size(mt_solid_soisno,1), &
         lbound(mt_solid_soisno,2):ubound(mt_solid_soisno,2), size(mt_solid_soisno,3)))
      lulcc_mt_solid_soisno_old = mt_solid_soisno
      allocate(lulcc_mt_canopy_solid_old(size(mt_canopy_solid,1), size(mt_canopy_solid,2)))
      lulcc_mt_canopy_solid_old = mt_canopy_solid
      allocate(lulcc_mt_surface_solid_old(size(mt_surface_solid,1), size(mt_surface_solid,2)))
      lulcc_mt_surface_solid_old = mt_surface_solid
      allocate(lulcc_mt_subsurface_solid_old(size(mt_subsurface_solid,1), size(mt_subsurface_solid,2)))
      lulcc_mt_subsurface_solid_old = mt_subsurface_solid
      allocate(lulcc_mt_waterstorage_solid_old(size(mt_waterstorage_solid,1), size(mt_waterstorage_solid,2)))
      lulcc_mt_waterstorage_solid_old = mt_waterstorage_solid
      allocate(lulcc_mt_waterstorage_old(size(mt_waterstorage,1), size(mt_waterstorage,2)))
      lulcc_mt_waterstorage_old = mt_waterstorage
      allocate(lulcc_mt_scv_old(size(mt_scv,1), size(mt_scv,2)))
      lulcc_mt_scv_old = mt_scv
      allocate(lulcc_mt_leaf_delta_e_old(size(mt_leaf_delta_e,1), size(mt_leaf_delta_e,2)))
      lulcc_mt_leaf_delta_e_old = mt_leaf_delta_e
      allocate(lulcc_mt_leaf_delta_b_old(size(mt_leaf_delta_b,1), size(mt_leaf_delta_b,2)))
      lulcc_mt_leaf_delta_b_old = mt_leaf_delta_b
      allocate(lulcc_mt_leaf_peclet_old(size(mt_leaf_peclet,1), size(mt_leaf_peclet,2)))
      lulcc_mt_leaf_peclet_old = mt_leaf_peclet
      allocate(lulcc_mt_leaf_water_moles_old(size(mt_leaf_water_moles,1), size(mt_leaf_water_moles,2)))
      lulcc_mt_leaf_water_moles_old = mt_leaf_water_moles
      allocate(lulcc_mt_leaf_iso_storage_old(size(mt_leaf_iso_storage,1), size(mt_leaf_iso_storage,2)))
      lulcc_mt_leaf_iso_storage_old = mt_leaf_iso_storage

      land_mt_lulcc_snapshot_valid = .true.
   END SUBROUTINE save_land_mt_lulcc_state

   SUBROUTINE remap_land_mt_lulcc_state (patchclass_new, eindex_new, patchclass_old, eindex_old, &
                                             wa_new, host_water_old_2d, host_water_old_3d, &
                                             host_sync, lccpct_patches, new_patch_area, old_patch_area)
#ifdef LULCC
      USE MOD_MT_Isotope_Lulcc, only: mt_isotope_lulcc_aquifer_reference
#endif
      IMPLICIT NONE
      integer,   intent(in) :: patchclass_new(:), patchclass_old(:)
      integer*8, intent(in) :: eindex_new(:), eindex_old(:)
      real(r8),  intent(in) :: wa_new(:)      ! host aquifer water on the new patches
      ! Host water of the old patches, one row per compartment (MOD_MT_Lulcc
      ! builds and reads the rows). It is remapped with the same weights as
      ! the material pools, so host_sync can see where the host did not carry
      ! its water over the same way.
      real(r8),  intent(in) :: host_water_old_2d(:,:)    ! (compartment, old patch)
      real(r8),  intent(in) :: host_water_old_3d(:,:,:)  ! (phase, layer, old patch)
      ! Brings the remapped pools in line with what the host LULCC did to its
      ! water (MOD_MT_Lulcc): cold start of patches without a source patch,
      ! wetland water reset, material following the host water wherever the
      ! host water differs from the remapped old host water. It adds the
      ! material mass it changes, per unit area times patch area, to change.
      INTERFACE
         SUBROUTINE host_sync (no_source, patch_area, change, src_water_2d, src_water_3d)
            USE MOD_Precision, only: r8
            logical,  intent(in)    :: no_source(:)
            real(r8), intent(in)    :: patch_area(:)
            real(r8), intent(inout) :: change(:)
            real(r8), intent(in)    :: src_water_2d(:,:)
            real(r8), intent(in)    :: src_water_3d(:,:,:)
         END SUBROUTINE host_sync
      END INTERFACE
      real(r8),  intent(in), optional :: lccpct_patches(:,0:)
      real(r8),  intent(in), optional :: new_patch_area(:), old_patch_area(:)

      integer :: nnew, maxsnl_saved, nl_soil_saved, np
      real(r8), allocatable :: mass_before(:), mass_after(:), carrier_change(:), host_change(:)
      real(r8), allocatable :: remapped_ref_water(:,:), sync_area(:)
      real(r8), allocatable :: src_water_2d(:,:), src_water_3d(:,:,:)
      logical,  allocatable :: no_source(:)
      logical :: check_lulcc_mass

      IF (.not. land_mt_lulcc_snapshot_valid) RETURN
      IF (.not. allocated(lulcc_mt_wliq_soisno_old)) THEN
         CALL clear_land_mt_lulcc_snapshot()
         RETURN
      ENDIF

      nnew = size(patchclass_new)
      maxsnl_saved = lbound(lulcc_mt_wliq_soisno_old, 2) - 1
      nl_soil_saved = ubound(lulcc_mt_wliq_soisno_old, 2)
      check_lulcc_mass = present(lccpct_patches) .and. present(new_patch_area) .and. &
         present(old_patch_area)
      IF (check_lulcc_mass) THEN
         allocate(mass_before(nmts), mass_after(nmts))
         CALL compute_lulcc_snapshot_land_water_mass(old_patch_area, mass_before)
      ENDIF

      IF (allocated(mt_ldew_rain)) CALL deallocate_Material_Vars()
      IF (nnew <= 0) THEN
         IF (allocated(mass_before)) deallocate(mass_before, mass_after)
         CALL clear_land_mt_lulcc_snapshot()
         RETURN
      ENDIF

      CALL allocate_Material_Vars(nnew, maxsnl_saved, nl_soil_saved)
      IF (.not. allocated(mt_ldew_rain)) THEN
         IF (allocated(mass_before)) deallocate(mass_before, mass_after)
         CALL clear_land_mt_lulcc_snapshot()
         RETURN
      ENDIF

      ! Extensive material pools are stored per unit patch area.  When both old
      ! and new patch areas are available, remap them by conserving area-mass:
      ! old_pool * old_area is partitioned into the target class fractions and
      ! divided by the new patch area.  Dimensionless/intensive diagnostics
      ! remain area-weighted means.
      CALL remap2d_mass(lulcc_mt_ldew_rain_old,         mt_ldew_rain)
      CALL remap2d_mass(lulcc_mt_ldew_snow_old,         mt_ldew_snow)
      CALL remap3d_mass(lulcc_mt_wliq_soisno_old,       mt_wliq_soisno)
      CALL remap3d_mass(lulcc_mt_wice_soisno_old,       mt_wice_soisno)
      CALL remap2d_mass(lulcc_mt_wa_old,                mt_wa)
      allocate(remapped_ref_water(1, nnew))
      remapped_ref_water = 0._r8
      CALL remap2d_mass(lulcc_mt_aquifer_ref_water_old, remapped_ref_water)
      mt_aquifer_ref_water = remapped_ref_water(1, :)
      deallocate(remapped_ref_water)
      CALL remap2d_mass(lulcc_mt_aquifer_ref_mass_old, mt_aquifer_ref_mass)
      CALL remap2d_mass(lulcc_mt_wdsrf_old,             mt_wdsrf)
      CALL remap2d_mass(lulcc_mt_wetwat_old,            mt_wetwat)
      CALL remap2d_mass(lulcc_mt_surface_residue_old,    mt_surface_residue)
      CALL remap2d_mass(lulcc_mt_subsurface_residue_old, mt_subsurface_residue)
      CALL remap3d_mass(lulcc_mt_solid_soisno_old,       mt_solid_soisno)
      CALL remap2d_mass(lulcc_mt_canopy_solid_old,       mt_canopy_solid)
      CALL remap2d_mass(lulcc_mt_surface_solid_old,      mt_surface_solid)
      CALL remap2d_mass(lulcc_mt_subsurface_solid_old,   mt_subsurface_solid)
      CALL remap2d_mass(lulcc_mt_waterstorage_solid_old, mt_waterstorage_solid)
      CALL remap2d_mass(lulcc_mt_waterstorage_old,      mt_waterstorage)
      CALL remap2d_mass(lulcc_mt_scv_old,               mt_scv)
      CALL remap2d_intensive(lulcc_mt_leaf_delta_e_old,      mt_leaf_delta_e)
      CALL remap2d_intensive(lulcc_mt_leaf_delta_b_old,      mt_leaf_delta_b)
      CALL remap2d_intensive(lulcc_mt_leaf_peclet_old,       mt_leaf_peclet)
      CALL remap2d_mass(lulcc_mt_leaf_water_moles_old,  mt_leaf_water_moles)
      CALL remap2d_mass(lulcc_mt_leaf_iso_storage_old,  mt_leaf_iso_storage)

      ! The host re-initialises every patch for the new year and copies back
      ! only some of its water (MOD_Lulcc_Initialize, REST_LulccTimeVariables,
      ! LulccMassEnergyConserve). The material pools follow it here; what
      ! that adds or removes is counted by the mass check below.
      allocate(no_source(nnew), host_change(nmts))
      DO np = 1, nnew
         no_source(np) = .not. patch_has_source(np)
      ENDDO
      ! The old host water on the new patches, moved exactly as the pools.
      allocate(src_water_2d(size(host_water_old_2d,1), nnew))
      allocate(src_water_3d(size(host_water_old_3d,1), size(host_water_old_3d,2), nnew))
      src_water_2d = 0._r8
      src_water_3d = 0._r8
      CALL remap2d_mass(host_water_old_2d, src_water_2d)
      CALL remap3d_mass(host_water_old_3d, src_water_3d)
      host_change = 0._r8
      IF (check_lulcc_mass) THEN
         CALL host_sync (no_source, new_patch_area, host_change, src_water_2d, src_water_3d)
      ELSE
         allocate(sync_area(nnew)); sync_area = 0._r8
         CALL host_sync (no_source, sync_area, host_change, src_water_2d, src_water_3d)
         deallocate(sync_area)
      ENDIF
      deallocate(src_water_2d, src_water_3d)

      ! A special patch that received reference water, or a soil/wetland patch
      ! that received none, gets its isotope reference carrier fixed in the
      ! isotope family; the carrier mass that leaves or enters is counted by
      ! the mass check below. It runs after the host sync, so a
      ! cold-started patch already has its reference.
      allocate(carrier_change(nmts))
      carrier_change = 0._r8
#ifdef LULCC
      IF (check_lulcc_mass) THEN
         CALL mt_isotope_lulcc_aquifer_reference (patchclass_new, wa_new, &
            mt_aquifer_ref_water, mt_aquifer_ref_mass, mt_wa, new_patch_area, carrier_change)
      ELSE
         CALL mt_isotope_lulcc_aquifer_reference (patchclass_new, wa_new, &
            mt_aquifer_ref_water, mt_aquifer_ref_mass, mt_wa)
      ENDIF
#endif
      IF (check_lulcc_mass) THEN
         CALL compute_current_lulcc_land_water_mass(new_patch_area, mass_after)
         CALL assert_lulcc_land_water_mass_conserved(mass_before + host_change + carrier_change, mass_after)
         deallocate(mass_before, mass_after)
      ENDIF

      deallocate(carrier_change, host_change, no_source)
      CALL clear_land_mt_lulcc_snapshot()

   CONTAINS
      SUBROUTINE compute_lulcc_snapshot_land_water_mass(area, total)
         real(r8), intent(in)  :: area(:)
         real(r8), intent(out) :: total(:)

         total = 0._r8
         IF (allocated(lulcc_mt_ldew_rain_old))    CALL accumulate_lulcc_mass_2d(lulcc_mt_ldew_rain_old, area, total)
         IF (allocated(lulcc_mt_ldew_snow_old))    CALL accumulate_lulcc_mass_2d(lulcc_mt_ldew_snow_old, area, total)
         IF (allocated(lulcc_mt_wliq_soisno_old))  CALL accumulate_lulcc_mass_3d(lulcc_mt_wliq_soisno_old, area, total)
         IF (allocated(lulcc_mt_wice_soisno_old))  CALL accumulate_lulcc_mass_3d(lulcc_mt_wice_soisno_old, area, total)
         IF (allocated(lulcc_mt_wa_old))           CALL accumulate_lulcc_mass_2d(lulcc_mt_wa_old, area, total)
         IF (allocated(lulcc_mt_aquifer_ref_mass_old)) &
            CALL accumulate_lulcc_mass_2d(lulcc_mt_aquifer_ref_mass_old, area, total)
         IF (allocated(lulcc_mt_wdsrf_old))        CALL accumulate_lulcc_mass_2d(lulcc_mt_wdsrf_old, area, total)
         IF (allocated(lulcc_mt_wetwat_old))       CALL accumulate_lulcc_mass_2d(lulcc_mt_wetwat_old, area, total)
         IF (allocated(lulcc_mt_surface_residue_old)) &
            CALL accumulate_lulcc_mass_2d(lulcc_mt_surface_residue_old, area, total)
         IF (allocated(lulcc_mt_subsurface_residue_old)) &
            CALL accumulate_lulcc_mass_2d(lulcc_mt_subsurface_residue_old, area, total)
         IF (allocated(lulcc_mt_solid_soisno_old)) &
            CALL accumulate_lulcc_mass_3d(lulcc_mt_solid_soisno_old, area, total)
         IF (allocated(lulcc_mt_canopy_solid_old)) &
            CALL accumulate_lulcc_mass_2d(lulcc_mt_canopy_solid_old, area, total)
         IF (allocated(lulcc_mt_surface_solid_old)) &
            CALL accumulate_lulcc_mass_2d(lulcc_mt_surface_solid_old, area, total)
         IF (allocated(lulcc_mt_subsurface_solid_old)) &
            CALL accumulate_lulcc_mass_2d(lulcc_mt_subsurface_solid_old, area, total)
         IF (allocated(lulcc_mt_waterstorage_solid_old)) &
            CALL accumulate_lulcc_mass_2d(lulcc_mt_waterstorage_solid_old, area, total)
         IF (allocated(lulcc_mt_waterstorage_old)) CALL accumulate_lulcc_mass_2d(lulcc_mt_waterstorage_old, area, total)
         IF (allocated(lulcc_mt_scv_old))          CALL accumulate_lulcc_mass_2d(lulcc_mt_scv_old, area, total)
         IF (allocated(lulcc_mt_leaf_iso_storage_old)) &
            CALL accumulate_lulcc_mass_2d(lulcc_mt_leaf_iso_storage_old, area, total)
      END SUBROUTINE compute_lulcc_snapshot_land_water_mass

      SUBROUTINE compute_current_lulcc_land_water_mass(area, total)
         real(r8), intent(in)  :: area(:)
         real(r8), intent(out) :: total(:)

         total = 0._r8
         IF (allocated(mt_ldew_rain))    CALL accumulate_lulcc_mass_2d(mt_ldew_rain, area, total)
         IF (allocated(mt_ldew_snow))    CALL accumulate_lulcc_mass_2d(mt_ldew_snow, area, total)
         IF (allocated(mt_wliq_soisno))  CALL accumulate_lulcc_mass_3d(mt_wliq_soisno, area, total)
         IF (allocated(mt_wice_soisno))  CALL accumulate_lulcc_mass_3d(mt_wice_soisno, area, total)
         IF (allocated(mt_wa))           CALL accumulate_lulcc_mass_2d(mt_wa, area, total)
         IF (allocated(mt_aquifer_ref_mass)) &
            CALL accumulate_lulcc_mass_2d(mt_aquifer_ref_mass, area, total)
         IF (allocated(mt_wdsrf))        CALL accumulate_lulcc_mass_2d(mt_wdsrf, area, total)
         IF (allocated(mt_wetwat))       CALL accumulate_lulcc_mass_2d(mt_wetwat, area, total)
         IF (allocated(mt_surface_residue)) CALL accumulate_lulcc_mass_2d(mt_surface_residue, area, total)
         IF (allocated(mt_subsurface_residue)) CALL accumulate_lulcc_mass_2d(mt_subsurface_residue, area, total)
         IF (allocated(mt_solid_soisno)) CALL accumulate_lulcc_mass_3d(mt_solid_soisno, area, total)
         IF (allocated(mt_canopy_solid)) CALL accumulate_lulcc_mass_2d(mt_canopy_solid, area, total)
         IF (allocated(mt_surface_solid)) CALL accumulate_lulcc_mass_2d(mt_surface_solid, area, total)
         IF (allocated(mt_subsurface_solid)) CALL accumulate_lulcc_mass_2d(mt_subsurface_solid, area, total)
         IF (allocated(mt_waterstorage_solid)) CALL accumulate_lulcc_mass_2d(mt_waterstorage_solid, area, total)
         IF (allocated(mt_waterstorage)) CALL accumulate_lulcc_mass_2d(mt_waterstorage, area, total)
         IF (allocated(mt_scv))          CALL accumulate_lulcc_mass_2d(mt_scv, area, total)
         IF (allocated(mt_leaf_iso_storage)) &
            CALL accumulate_lulcc_mass_2d(mt_leaf_iso_storage, area, total)
      END SUBROUTINE compute_current_lulcc_land_water_mass

      SUBROUTINE accumulate_lulcc_mass_2d(pool, area, total)
         real(r8), intent(in)    :: pool(:,:)
         real(r8), intent(in)    :: area(:)
         real(r8), intent(inout) :: total(:)
         integer :: ip, imt, n1
         real(r8) :: patch_area

         n1 = min(size(pool,1), size(total), nmts)
         DO ip = 1, min(size(pool,2), size(area))
            patch_area = max(0._r8, area(ip))
            IF (patch_area <= tiny(1._r8)) CYCLE
            DO imt = 1, n1
               IF (.not. mt_uses_generic_water(imt)) CYCLE
               total(imt) = total(imt) + pool(imt, ip) * patch_area
            ENDDO
         ENDDO
      END SUBROUTINE accumulate_lulcc_mass_2d

      SUBROUTINE accumulate_lulcc_mass_3d(pool, area, total)
         real(r8), intent(in)    :: pool(:,:,:)
         real(r8), intent(in)    :: area(:)
         real(r8), intent(inout) :: total(:)
         integer :: ip, imt, iz, n1
         real(r8) :: patch_area

         n1 = min(size(pool,1), size(total), nmts)
         DO ip = 1, min(size(pool,3), size(area))
            patch_area = max(0._r8, area(ip))
            IF (patch_area <= tiny(1._r8)) CYCLE
            DO iz = 1, size(pool,2)
               DO imt = 1, n1
                  IF (.not. mt_uses_generic_water(imt)) CYCLE
                  total(imt) = total(imt) + pool(imt, iz, ip) * patch_area
               ENDDO
            ENDDO
         ENDDO
      END SUBROUTINE accumulate_lulcc_mass_3d

      SUBROUTINE assert_lulcc_land_water_mass_conserved(before, after)
         real(r8), intent(in) :: before(:), after(:)
         integer :: imt, nbad, worst_imt
         real(r8) :: diff, scale, tol, abs_err, rel_err, worst_abs

         nbad = 0
         worst_imt = 0
         worst_abs = 0._r8
         DO imt = 1, min(size(before), size(after), nmts)
            IF (.not. mt_uses_generic_water(imt)) CYCLE
            diff = after(imt) - before(imt)
            abs_err = abs(diff)
            scale = max(abs(before(imt)), abs(after(imt)), 1._r8)
            tol = max(lulcc_mass_abs_tol, lulcc_mass_rel_tol * scale)
            IF (abs_err > worst_abs) THEN
               worst_abs = abs_err
               worst_imt = imt
               rel_err = abs_err / scale
            ENDIF
            IF (abs_err > tol) nbad = nbad + 1
         ENDDO

         IF (nbad > DEF_MT_LULCC_ABORT_NBAD) THEN
            IF (p_is_io) THEN
               WRITE(*,'(A,I0,A,E12.5,A,E12.5)') &
                  'MT_LULCC_BAL failed: nbad=', nbad, ' worst_abs=', worst_abs, &
                  ' worst_rel=', rel_err
               IF (worst_imt > 0) THEN
                  WRITE(*,'(A,I0,3A,2(A,E12.5))') '  material=', worst_imt, &
                     ' name=', trim(materials(worst_imt)%name), ' ', &
                     ' before=', before(worst_imt), ' after=', after(worst_imt)
               ENDIF
            ENDIF
            CALL CoLM_stop('MT LULCC remap land-water mass conservation failed')
         ENDIF
      END SUBROUTINE assert_lulcc_land_water_mass_conserved

      SUBROUTINE remap2d_mass(old, new)
         real(r8), intent(in)    :: old(:,:)
         real(r8), intent(inout) :: new(:,:)
         integer :: np, op, src, n1
         real(r8) :: w, wsum, denom
         logical :: conserve_area_mass
         real(r8) :: default_vals(size(new,1))

         n1 = min(size(new,1), size(old,1))
         DO np = 1, min(size(new,2), nnew)
            wsum = 0._r8
            default_vals(:) = new(:,np)
            IF (present(lccpct_patches)) THEN
               new(:,np) = 0._r8
               conserve_area_mass = area_mass_remap_available(np)
               DO op = 1, min(size(old,2), size(patchclass_old), size(eindex_old))
                  IF (eindex_old(op) /= eindex_new(np)) CYCLE
                  IF (patchclass_old(op) < lbound(lccpct_patches,2) .or. &
                      patchclass_old(op) > ubound(lccpct_patches,2)) CYCLE
                  IF (conserve_area_mass) THEN
                     w = lulcc_mass_transfer_area(np, op)
                  ELSE
                     w = lulcc_source_area_weight(np, op)
                  ENDIF
                  IF (w <= 0._r8) CYCLE
                  new(1:n1,np) = new(1:n1,np) + w * old(1:n1,op)
                  wsum = wsum + w
               ENDDO
            ENDIF
            IF (wsum > 0._r8) THEN
               denom = remap_denominator(np, wsum, conserve_area_mass)
               new(:,np) = new(:,np) / denom
            ELSE
               src = fallback_source(np, size(old,2))
               IF (src > 0) THEN
                  new(:,np) = 0._r8
                  new(1:n1,np) = old(1:n1,src)
               ELSE
                  new(:,np) = default_vals(:)
               ENDIF
            ENDIF
         ENDDO
      END SUBROUTINE remap2d_mass

      SUBROUTINE remap3d_mass(old, new)
         real(r8), intent(in)    :: old(:,:,:)
         real(r8), intent(inout) :: new(:,:,:)
         integer :: np, op, src, n1, n2
         real(r8) :: w, wsum, denom
         logical :: conserve_area_mass
         real(r8) :: default_vals(size(new,1), size(new,2))

         n1 = min(size(new,1), size(old,1))
         n2 = min(size(new,2), size(old,2))
         DO np = 1, min(size(new,3), nnew)
            wsum = 0._r8
            default_vals(:,:) = new(:,:,np)
            IF (present(lccpct_patches)) THEN
               new(:,:,np) = 0._r8
               conserve_area_mass = area_mass_remap_available(np)
               DO op = 1, min(size(old,3), size(patchclass_old), size(eindex_old))
                  IF (eindex_old(op) /= eindex_new(np)) CYCLE
                  IF (patchclass_old(op) < lbound(lccpct_patches,2) .or. &
                      patchclass_old(op) > ubound(lccpct_patches,2)) CYCLE
                  IF (conserve_area_mass) THEN
                     w = lulcc_mass_transfer_area(np, op)
                  ELSE
                     w = lulcc_source_area_weight(np, op)
                  ENDIF
                  IF (w <= 0._r8) CYCLE
                  new(1:n1,1:n2,np) = new(1:n1,1:n2,np) + w * old(1:n1,1:n2,op)
                  wsum = wsum + w
               ENDDO
            ENDIF
            IF (wsum > 0._r8) THEN
               denom = remap_denominator(np, wsum, conserve_area_mass)
               new(:,:,np) = new(:,:,np) / denom
            ELSE
               src = fallback_source(np, size(old,3))
               IF (src > 0) THEN
                  new(:,:,np) = 0._r8
                  new(1:n1,1:n2,np) = old(1:n1,1:n2,src)
               ELSE
                  new(:,:,np) = default_vals(:,:)
               ENDIF
            ENDIF
         ENDDO
      END SUBROUTINE remap3d_mass

      SUBROUTINE remap2d_intensive(old, new)
         real(r8), intent(in)    :: old(:,:)
         real(r8), intent(inout) :: new(:,:)
         integer :: np, op, src, n1
         real(r8) :: w, wsum
         real(r8) :: default_vals(size(new,1))

         n1 = min(size(new,1), size(old,1))
         DO np = 1, min(size(new,2), nnew)
            wsum = 0._r8
            default_vals(:) = new(:,np)
            IF (present(lccpct_patches)) THEN
               new(:,np) = 0._r8
               DO op = 1, min(size(old,2), size(patchclass_old), size(eindex_old))
                  IF (eindex_old(op) /= eindex_new(np)) CYCLE
                  IF (patchclass_old(op) < lbound(lccpct_patches,2) .or. &
                      patchclass_old(op) > ubound(lccpct_patches,2)) CYCLE
                  w = lulcc_source_area_weight(np, op)
                  IF (w <= 0._r8) CYCLE
                  new(1:n1,np) = new(1:n1,np) + w * old(1:n1,op)
                  wsum = wsum + w
               ENDDO
            ENDIF
            IF (wsum > 0._r8) THEN
               new(:,np) = new(:,np) / remap_denominator(np, wsum, .false.)
            ELSE
               src = fallback_source(np, size(old,2))
               IF (src > 0) THEN
                  new(:,np) = 0._r8
                  new(1:n1,np) = old(1:n1,src)
               ELSE
                  new(:,np) = default_vals(:)
               ENDIF
            ENDIF
         ENDDO
      END SUBROUTINE remap2d_intensive

      integer FUNCTION fallback_source(np, old_n) RESULT(src)
         integer, intent(in) :: np, old_n
         integer :: op
         src = 0
         DO op = 1, min(old_n, size(patchclass_old), size(eindex_old))
            IF (eindex_old(op) == eindex_new(np) .and. &
                patchclass_old(op) == patchclass_new(np)) THEN
               src = op
               RETURN
            ENDIF
         ENDDO
         ! Do not fall back across land-cover classes; a missing class match
         ! should cold-start from allocation defaults instead of contaminating
         ! a newly created patch with another surface type's material pools.
      END FUNCTION fallback_source

      ! Whether the remap gives new patch np any source: a positive transfer
      ! weight (MEC) or an old patch of the same class in the element. The
      ! same tests as remap2d_mass.
      logical FUNCTION patch_has_source(np) RESULT(has)
         integer, intent(in) :: np
         integer :: op
         real(r8) :: w
         logical :: conserve_area_mass

         has = .false.
         IF (present(lccpct_patches)) THEN
            conserve_area_mass = area_mass_remap_available(np)
            DO op = 1, min(size(patchclass_old), size(eindex_old))
               IF (eindex_old(op) /= eindex_new(np)) CYCLE
               IF (patchclass_old(op) < lbound(lccpct_patches,2) .or. &
                   patchclass_old(op) > ubound(lccpct_patches,2)) CYCLE
               IF (conserve_area_mass) THEN
                  w = lulcc_mass_transfer_area(np, op)
               ELSE
                  w = lulcc_source_area_weight(np, op)
               ENDIF
               IF (w > 0._r8) THEN
                  has = .true.
                  RETURN
               ENDIF
            ENDDO
         ENDIF
         has = fallback_source(np, size(patchclass_old)) > 0
      END FUNCTION patch_has_source

      REAL(r8) FUNCTION lulcc_source_area_weight(np, op) RESULT(w)
         integer, intent(in) :: np, op
         integer :: oq
         real(r8) :: class_area

         w = 0._r8
         IF (.not. present(lccpct_patches)) RETURN
         IF (np > size(lccpct_patches,1)) RETURN
         IF (patchclass_old(op) < lbound(lccpct_patches,2) .or. &
             patchclass_old(op) > ubound(lccpct_patches,2)) RETURN
         w = max(0._r8, lccpct_patches(np, patchclass_old(op)))
         IF (w <= 0._r8) RETURN
         IF (.not. present(old_patch_area)) RETURN
         IF (op > size(old_patch_area)) RETURN

         class_area = 0._r8
         DO oq = 1, min(size(patchclass_old), size(eindex_old), size(old_patch_area))
            IF (eindex_old(oq) == eindex_new(np) .and. &
                patchclass_old(oq) == patchclass_old(op)) THEN
               class_area = class_area + max(0._r8, old_patch_area(oq))
            ENDIF
         ENDDO
         IF (class_area > tiny(1._r8)) THEN
            w = w * max(0._r8, old_patch_area(op)) / class_area
         ENDIF
      END FUNCTION lulcc_source_area_weight

      LOGICAL FUNCTION area_mass_remap_available(np) RESULT(ok)
         integer, intent(in) :: np

         ok = present(lccpct_patches) .and. present(new_patch_area) .and. present(old_patch_area)
         IF (.not. ok) THEN
            CALL warn_lulcc_area_fallback('missing old/new patch area arrays')
            RETURN
         ENDIF
         ok = np <= size(new_patch_area)
         IF (.not. ok) THEN
            CALL warn_lulcc_area_fallback('new patch index exceeds new_patch_area size')
            RETURN
         ENDIF
         ok = new_patch_area(np) > tiny(1._r8)
         IF (.not. ok) CALL warn_lulcc_area_fallback('new patch area is zero')
      END FUNCTION area_mass_remap_available

      SUBROUTINE warn_lulcc_area_fallback(reason)
         character(len=*), intent(in) :: reason

         IF (lulcc_area_fallback_warned) RETURN
         lulcc_area_fallback_warned = .true.
         IF (p_is_io) THEN
            WRITE(*,'(A,A)') 'MT LULCC remap falling back to source-area weighting: ', &
               trim(reason)
            WRITE(*,'(A)') '  Land-water material area-mass conservation requires old and new patch areas.'
         ENDIF
      END SUBROUTINE warn_lulcc_area_fallback

      REAL(r8) FUNCTION lulcc_mass_transfer_area(np, op) RESULT(w)
         integer, intent(in) :: np, op
         integer :: c, nq
         real(r8) :: target_area, class_target_area

         w = 0._r8
         IF (.not. area_mass_remap_available(np)) RETURN
         IF (op > size(old_patch_area)) RETURN
         IF (op > size(patchclass_old) .or. op > size(eindex_old)) RETURN
         c = patchclass_old(op)
         IF (c < lbound(lccpct_patches,2) .or. c > ubound(lccpct_patches,2)) RETURN

         target_area = lulcc_target_class_area(np, c)
         IF (target_area <= tiny(1._r8)) RETURN

         class_target_area = 0._r8
         DO nq = 1, min(nnew, size(eindex_new), size(new_patch_area))
            IF (eindex_new(nq) == eindex_new(np)) THEN
               class_target_area = class_target_area + lulcc_target_class_area(nq, c)
            ENDIF
         ENDDO
         IF (class_target_area <= tiny(1._r8)) RETURN

         ! Area of old patch op that is transferred into target patch np.
         ! Summing this weight over all target patches of the same element and
         ! source class returns old_patch_area(op), so extensive material mass is
         ! conserved when target area is non-zero.
         w = max(0._r8, old_patch_area(op)) * target_area / class_target_area
      END FUNCTION lulcc_mass_transfer_area

      REAL(r8) FUNCTION lulcc_target_class_area(np, c) RESULT(area)
         integer, intent(in) :: np, c
         integer :: cc
         real(r8) :: class_sum

         area = 0._r8
         IF (.not. present(lccpct_patches)) RETURN
         IF (.not. present(new_patch_area)) RETURN
         IF (np > size(new_patch_area)) RETURN
         IF (np > size(lccpct_patches,1)) RETURN
         IF (c < lbound(lccpct_patches,2) .or. c > ubound(lccpct_patches,2)) RETURN

         class_sum = 0._r8
         DO cc = lbound(lccpct_patches,2), ubound(lccpct_patches,2)
            class_sum = class_sum + max(0._r8, lccpct_patches(np, cc))
         ENDDO
         IF (class_sum <= tiny(1._r8)) RETURN

         area = max(0._r8, new_patch_area(np)) * max(0._r8, lccpct_patches(np, c)) / class_sum
      END FUNCTION lulcc_target_class_area

      REAL(r8) FUNCTION remap_denominator(np, wsum, conserve_mass) RESULT(denom)
         integer, intent(in) :: np
         real(r8), intent(in) :: wsum
         logical, intent(in) :: conserve_mass

         denom = max(wsum, tiny(1._r8))
         IF (conserve_mass .and. present(new_patch_area)) THEN
            IF (np <= size(new_patch_area) .and. new_patch_area(np) > tiny(1._r8)) THEN
               denom = new_patch_area(np)
            ENDIF
         ENDIF
      END FUNCTION remap_denominator
   END SUBROUTINE remap_land_mt_lulcc_state

   ! Land-water material mass of one patch per unit area, over the same pools
   ! as the LULCC mass check (remap_land_mt_lulcc_state).
   SUBROUTINE mt_patch_land_water_mass (ip, total)
      IMPLICIT NONE
      integer,  intent(in)  :: ip
      real(r8), intent(out) :: total(:)
      integer :: imt

      total = 0._r8
      DO imt = 1, min(size(total), nmts)
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         IF (allocated(mt_ldew_rain))    total(imt) = total(imt) + mt_ldew_rain(imt,ip)
         IF (allocated(mt_ldew_snow))    total(imt) = total(imt) + mt_ldew_snow(imt,ip)
         IF (allocated(mt_wliq_soisno))  total(imt) = total(imt) + sum(mt_wliq_soisno(imt,:,ip))
         IF (allocated(mt_wice_soisno))  total(imt) = total(imt) + sum(mt_wice_soisno(imt,:,ip))
         IF (allocated(mt_wa))           total(imt) = total(imt) + mt_wa(imt,ip)
         IF (allocated(mt_aquifer_ref_mass)) total(imt) = total(imt) + mt_aquifer_ref_mass(imt,ip)
         IF (allocated(mt_wdsrf))        total(imt) = total(imt) + mt_wdsrf(imt,ip)
         IF (allocated(mt_wetwat))       total(imt) = total(imt) + mt_wetwat(imt,ip)
         IF (allocated(mt_surface_residue)) total(imt) = total(imt) + mt_surface_residue(imt,ip)
         IF (allocated(mt_subsurface_residue)) total(imt) = total(imt) + mt_subsurface_residue(imt,ip)
         IF (allocated(mt_solid_soisno)) total(imt) = total(imt) + sum(mt_solid_soisno(imt,:,ip))
         IF (allocated(mt_canopy_solid)) total(imt) = total(imt) + mt_canopy_solid(imt,ip)
         IF (allocated(mt_surface_solid)) total(imt) = total(imt) + mt_surface_solid(imt,ip)
         IF (allocated(mt_subsurface_solid)) total(imt) = total(imt) + mt_subsurface_solid(imt,ip)
         IF (allocated(mt_waterstorage_solid)) total(imt) = total(imt) + mt_waterstorage_solid(imt,ip)
         IF (allocated(mt_waterstorage)) total(imt) = total(imt) + mt_waterstorage(imt,ip)
         IF (allocated(mt_scv))          total(imt) = total(imt) + mt_scv(imt,ip)
      ENDDO
   END SUBROUTINE mt_patch_land_water_mass

   SUBROUTINE clear_land_mt_lulcc_snapshot ()
      IMPLICIT NONE

      IF (allocated(lulcc_mt_ldew_rain_old)) deallocate(lulcc_mt_ldew_rain_old)
      IF (allocated(lulcc_mt_ldew_snow_old)) deallocate(lulcc_mt_ldew_snow_old)
      IF (allocated(lulcc_mt_wliq_soisno_old)) deallocate(lulcc_mt_wliq_soisno_old)
      IF (allocated(lulcc_mt_wice_soisno_old)) deallocate(lulcc_mt_wice_soisno_old)
      IF (allocated(lulcc_mt_wa_old)) deallocate(lulcc_mt_wa_old)
      IF (allocated(lulcc_mt_aquifer_ref_water_old)) deallocate(lulcc_mt_aquifer_ref_water_old)
      IF (allocated(lulcc_mt_aquifer_ref_mass_old)) deallocate(lulcc_mt_aquifer_ref_mass_old)
      IF (allocated(lulcc_mt_wdsrf_old)) deallocate(lulcc_mt_wdsrf_old)
      IF (allocated(lulcc_mt_wetwat_old)) deallocate(lulcc_mt_wetwat_old)
      IF (allocated(lulcc_mt_surface_residue_old)) deallocate(lulcc_mt_surface_residue_old)
      IF (allocated(lulcc_mt_subsurface_residue_old)) deallocate(lulcc_mt_subsurface_residue_old)
      IF (allocated(lulcc_mt_solid_soisno_old)) deallocate(lulcc_mt_solid_soisno_old)
      IF (allocated(lulcc_mt_canopy_solid_old)) deallocate(lulcc_mt_canopy_solid_old)
      IF (allocated(lulcc_mt_surface_solid_old)) deallocate(lulcc_mt_surface_solid_old)
      IF (allocated(lulcc_mt_subsurface_solid_old)) deallocate(lulcc_mt_subsurface_solid_old)
      IF (allocated(lulcc_mt_waterstorage_solid_old)) deallocate(lulcc_mt_waterstorage_solid_old)
      IF (allocated(lulcc_mt_waterstorage_old)) deallocate(lulcc_mt_waterstorage_old)
      IF (allocated(lulcc_mt_scv_old)) deallocate(lulcc_mt_scv_old)
      IF (allocated(lulcc_mt_leaf_delta_e_old)) deallocate(lulcc_mt_leaf_delta_e_old)
      IF (allocated(lulcc_mt_leaf_delta_b_old)) deallocate(lulcc_mt_leaf_delta_b_old)
      IF (allocated(lulcc_mt_leaf_peclet_old)) deallocate(lulcc_mt_leaf_peclet_old)
      IF (allocated(lulcc_mt_leaf_water_moles_old)) deallocate(lulcc_mt_leaf_water_moles_old)
      IF (allocated(lulcc_mt_leaf_iso_storage_old)) deallocate(lulcc_mt_leaf_iso_storage_old)
      land_mt_lulcc_snapshot_valid = .false.
   END SUBROUTINE clear_land_mt_lulcc_snapshot

   SUBROUTINE flush_Material_Acc ()
      IMPLICIT NONE
	      IF (allocated(a_mt_precip )) a_mt_precip  = 0._r8
	      IF (allocated(a_mt_vapor_exchange)) a_mt_vapor_exchange = 0._r8
	      IF (allocated(a_water_precip)) a_water_precip = 0._r8
		      IF (allocated(a_mt_evap   )) a_mt_evap    = 0._r8
		      IF (allocated(a_water_evap_gross)) a_water_evap_gross = 0._r8
		      IF (allocated(a_mt_transp )) a_mt_transp  = 0._r8
		      IF (allocated(a_mt_transp_src)) a_mt_transp_src = 0._r8
		      IF (allocated(a_water_transp)) a_water_transp = 0._r8
		      IF (allocated(a_mt_soilevap)) a_mt_soilevap = 0._r8
		      IF (allocated(a_water_soilevap)) a_water_soilevap = 0._r8
		      IF (allocated(a_mt_canopyevap)) a_mt_canopyevap = 0._r8
		      IF (allocated(a_water_canopyevap)) a_water_canopyevap = 0._r8
		      IF (allocated(a_mt_subl    )) a_mt_subl = 0._r8
		      IF (allocated(a_water_subl  )) a_water_subl = 0._r8
		      IF (allocated(a_mt_wetland_evap)) a_mt_wetland_evap = 0._r8
		      IF (allocated(a_water_wetland_evap)) a_water_wetland_evap = 0._r8
      IF (allocated(a_mt_rsur   )) a_mt_rsur    = 0._r8
      IF (allocated(a_mt_rsub   )) a_mt_rsub    = 0._r8
      IF (allocated(a_mt_bgc_leach)) a_mt_bgc_leach = 0._r8
      IF (allocated(a_mt_rnof   )) a_mt_rnof    = 0._r8
      IF (allocated(a_water_rnof )) a_water_rnof  = 0._r8
      IF (allocated(a_mt_qinfl  )) a_mt_qinfl   = 0._r8
      IF (allocated(a_mt_qcharge)) a_mt_qcharge  = 0._r8
      IF (allocated(a_mt_ldew_mass)) a_mt_ldew_mass = 0._r8
      IF (allocated(a_water_ldew   )) a_water_ldew    = 0._r8
      IF (allocated(a_mt_soil_mass)) a_mt_soil_mass = 0._r8
      IF (allocated(a_water_soil   )) a_water_soil    = 0._r8
      IF (allocated(a_mt_snow_mass)) a_mt_snow_mass = 0._r8
      IF (allocated(a_water_snow   )) a_water_snow    = 0._r8
	      IF (allocated(a_mt_wa_mass    )) a_mt_wa_mass     = 0._r8
	      IF (allocated(a_water_wa       )) a_water_wa        = 0._r8
	      IF (allocated(a_mt_wa_debt_mass)) a_mt_wa_debt_mass = 0._r8
	      IF (allocated(a_water_wa_debt   )) a_water_wa_debt    = 0._r8
      IF (allocated(a_water_aquifer_actual)) a_water_aquifer_actual = 0._r8
      IF (allocated(a_mt_aquifer_actual_mass)) a_mt_aquifer_actual_mass = 0._r8
      IF (allocated(a_mt_wdsrf_mass )) a_mt_wdsrf_mass  = 0._r8
      IF (allocated(a_water_wdsrf    )) a_water_wdsrf     = 0._r8
      IF (allocated(a_mt_wetwat_mass)) a_mt_wetwat_mass = 0._r8
      IF (allocated(a_mt_surface_residue_mass)) a_mt_surface_residue_mass = 0._r8
      IF (allocated(a_mt_subsurface_residue_mass)) a_mt_subsurface_residue_mass = 0._r8
      IF (allocated(a_mt_layer_dry_mass)) a_mt_layer_dry_mass = 0._r8
      IF (allocated(a_mt_solid_mass)) a_mt_solid_mass = 0._r8
      IF (allocated(a_water_wetwat   )) a_water_wetwat    = 0._r8
      IF (allocated(a_mt_scv_mass   )) a_mt_scv_mass    = 0._r8
      IF (allocated(a_water_scv      )) a_water_scv       = 0._r8
      IF (allocated(mt_pg_rain_ground)) mt_pg_rain_ground = 0._r8
      IF (allocated(mt_pg_snow_ground)) mt_pg_snow_ground = 0._r8
      IF (allocated(mt_rnof_step    )) mt_rnof_step     = 0._r8
	      IF (allocated(mt_sm_carry     )) mt_sm_carry      = 0._r8
	      IF (allocated(mt_decay_source_step)) mt_decay_source_step = 0._r8
	      IF (allocated(mt_numerical_residual_step)) mt_numerical_residual_step = 0._r8
	      IF (allocated(mt_numerical_water_step)) mt_numerical_water_step = 0._r8
	   END SUBROUTINE flush_Material_Acc

   SUBROUTINE zero_provider_owned_land_mt_state ()
      IMPLICIT NONE
      integer :: imt

      IF (nmts <= 0) RETURN

      DO imt = 1, nmts
         IF (mt_uses_generic_water(imt)) CYCLE

         IF (allocated(mt_ldew_rain  )) mt_ldew_rain  (imt, :)    = 0._r8
         IF (allocated(mt_ldew_snow  )) mt_ldew_snow  (imt, :)    = 0._r8
         IF (allocated(mt_wliq_soisno)) mt_wliq_soisno(imt, :, :) = 0._r8
         IF (allocated(mt_wice_soisno)) mt_wice_soisno(imt, :, :) = 0._r8
         IF (allocated(mt_wa         )) mt_wa         (imt, :)    = 0._r8
         IF (allocated(mt_aquifer_ref_mass)) mt_aquifer_ref_mass(imt, :) = 0._r8
         IF (allocated(mt_wdsrf      )) mt_wdsrf      (imt, :)    = 0._r8
         IF (allocated(mt_wetwat     )) mt_wetwat     (imt, :)    = 0._r8
         IF (allocated(mt_surface_residue)) mt_surface_residue(imt, :) = 0._r8
         IF (allocated(mt_subsurface_residue)) mt_subsurface_residue(imt, :) = 0._r8
         IF (allocated(mt_solid_soisno)) mt_solid_soisno(imt, :, :) = 0._r8
         IF (allocated(mt_canopy_solid)) mt_canopy_solid(imt, :) = 0._r8
         IF (allocated(mt_surface_solid)) mt_surface_solid(imt, :) = 0._r8
         IF (allocated(mt_subsurface_solid)) mt_subsurface_solid(imt, :) = 0._r8
         IF (allocated(mt_waterstorage_solid)) mt_waterstorage_solid(imt, :) = 0._r8
         IF (allocated(mt_scv        )) mt_scv        (imt, :)    = 0._r8
         IF (allocated(mt_waterstorage)) mt_waterstorage(imt, :)  = 0._r8
         IF (allocated(mt_pg_rain_ground)) mt_pg_rain_ground(imt, :) = 0._r8
         IF (allocated(mt_pg_snow_ground)) mt_pg_snow_ground(imt, :) = 0._r8
         IF (allocated(mt_rnof_step    )) mt_rnof_step    (imt, :) = 0._r8
         IF (allocated(mt_sm_carry     )) mt_sm_carry     (imt, :) = 0._r8
         IF (allocated(mt_leaf_delta_e )) mt_leaf_delta_e (imt, :) = 0._r8
         IF (allocated(mt_leaf_delta_b )) mt_leaf_delta_b (imt, :) = 0._r8
         IF (allocated(mt_leaf_peclet  )) mt_leaf_peclet  (imt, :) = 1._r8
         IF (allocated(mt_leaf_water_moles)) mt_leaf_water_moles(imt, :) = 0._r8
         IF (allocated(mt_leaf_iso_storage)) mt_leaf_iso_storage(imt, :) = 0._r8

         IF (allocated(a_mt_precip   )) a_mt_precip   (imt, :) = 0._r8
         IF (allocated(a_mt_vapor_exchange)) a_mt_vapor_exchange(imt, :) = 0._r8
         IF (allocated(a_water_precip )) a_water_precip (imt, :) = 0._r8
         IF (allocated(a_mt_evap     )) a_mt_evap     (imt, :) = 0._r8
         IF (allocated(a_water_evap_gross)) a_water_evap_gross(imt, :) = 0._r8
         IF (allocated(a_mt_transp   )) a_mt_transp   (imt, :) = 0._r8
         IF (allocated(a_mt_transp_src)) a_mt_transp_src(imt, :) = 0._r8
         IF (allocated(a_water_transp )) a_water_transp (imt, :) = 0._r8
         IF (allocated(a_mt_soilevap )) a_mt_soilevap (imt, :) = 0._r8
         IF (allocated(a_water_soilevap)) a_water_soilevap(imt, :) = 0._r8
         IF (allocated(a_mt_canopyevap)) a_mt_canopyevap(imt, :) = 0._r8
         IF (allocated(a_water_canopyevap)) a_water_canopyevap(imt, :) = 0._r8
         IF (allocated(a_mt_subl     )) a_mt_subl     (imt, :) = 0._r8
         IF (allocated(a_water_subl   )) a_water_subl   (imt, :) = 0._r8
         IF (allocated(a_mt_wetland_evap)) a_mt_wetland_evap(imt, :) = 0._r8
         IF (allocated(a_water_wetland_evap)) a_water_wetland_evap(imt, :) = 0._r8
         IF (allocated(a_mt_rsur     )) a_mt_rsur     (imt, :) = 0._r8
         IF (allocated(a_mt_rsub     )) a_mt_rsub     (imt, :) = 0._r8
         IF (allocated(a_mt_bgc_leach)) a_mt_bgc_leach(imt, :) = 0._r8
         IF (allocated(a_mt_rnof     )) a_mt_rnof     (imt, :) = 0._r8
         IF (allocated(a_water_rnof   )) a_water_rnof   (imt, :) = 0._r8
         IF (allocated(a_mt_qinfl    )) a_mt_qinfl    (imt, :) = 0._r8
         IF (allocated(a_mt_qcharge  )) a_mt_qcharge  (imt, :) = 0._r8
         IF (allocated(mt_storage_beg)) mt_storage_beg(imt, :) = 0._r8
         IF (allocated(mt_balance_err)) mt_balance_err(imt, :) = 0._r8
         IF (allocated(mt_decay_source_step)) mt_decay_source_step(imt, :) = 0._r8
         IF (allocated(mt_numerical_residual_step)) mt_numerical_residual_step(imt, :) = 0._r8
         IF (allocated(mt_numerical_water_step)) mt_numerical_water_step(imt, :) = 0._r8
         IF (allocated(a_mt_ldew_mass)) a_mt_ldew_mass(imt, :) = 0._r8
         IF (allocated(a_mt_soil_mass)) a_mt_soil_mass(imt, :, :) = 0._r8
         IF (allocated(a_mt_snow_mass)) a_mt_snow_mass(imt, :, :) = 0._r8
         IF (allocated(a_mt_wa_mass    )) a_mt_wa_mass    (imt, :) = 0._r8
         IF (allocated(a_mt_aquifer_actual_mass)) a_mt_aquifer_actual_mass(imt, :) = 0._r8
         IF (allocated(a_mt_wa_debt_mass)) a_mt_wa_debt_mass(imt, :) = 0._r8
         IF (allocated(a_mt_wdsrf_mass )) a_mt_wdsrf_mass (imt, :) = 0._r8
         IF (allocated(a_mt_wetwat_mass)) a_mt_wetwat_mass(imt, :) = 0._r8
         IF (allocated(a_mt_surface_residue_mass)) a_mt_surface_residue_mass(imt, :) = 0._r8
         IF (allocated(a_mt_subsurface_residue_mass)) a_mt_subsurface_residue_mass(imt, :) = 0._r8
         IF (allocated(a_mt_layer_dry_mass)) a_mt_layer_dry_mass(imt, :) = 0._r8
         IF (allocated(a_mt_solid_mass)) a_mt_solid_mass(imt, :) = 0._r8
         IF (allocated(a_mt_scv_mass   )) a_mt_scv_mass   (imt, :) = 0._r8
      ENDDO
   END SUBROUTINE zero_provider_owned_land_mt_state

	   SUBROUTINE material_book_evap_loss (imt, ipatch, material_mass, water_mass, evap_kind)
	      IMPLICIT NONE
	      integer,  intent(in) :: imt, ipatch
	      real(r8), intent(in) :: material_mass, water_mass
	      integer,  intent(in), optional :: evap_kind
	      integer :: kind

	      IF (nmts <= 0) RETURN
	      IF (.not. allocated(a_mt_evap)) RETURN
	      kind = MT_EVAP_KIND_TOTAL
	      IF (present(evap_kind)) kind = evap_kind

	      a_mt_evap(imt, ipatch) = a_mt_evap(imt, ipatch) + material_mass
	      IF (allocated(a_water_evap_gross)) THEN
	         a_water_evap_gross(imt, ipatch) = a_water_evap_gross(imt, ipatch) &
	            + max(water_mass, 0._r8)
	      ENDIF

	      SELECT CASE (kind)
	      CASE (MT_EVAP_KIND_TRANSP)
	         IF (allocated(a_mt_transp)) a_mt_transp(imt, ipatch) = &
	            a_mt_transp(imt, ipatch) + material_mass
	         IF (allocated(a_water_transp)) a_water_transp(imt, ipatch) = &
	            a_water_transp(imt, ipatch) + max(water_mass, 0._r8)
	      CASE (MT_EVAP_KIND_SOILEVAP)
	         IF (allocated(a_mt_soilevap)) a_mt_soilevap(imt, ipatch) = &
	            a_mt_soilevap(imt, ipatch) + material_mass
	         IF (allocated(a_water_soilevap)) a_water_soilevap(imt, ipatch) = &
	            a_water_soilevap(imt, ipatch) + max(water_mass, 0._r8)
	      CASE (MT_EVAP_KIND_CANOPYEVAP)
	         IF (allocated(a_mt_canopyevap)) a_mt_canopyevap(imt, ipatch) = &
	            a_mt_canopyevap(imt, ipatch) + material_mass
	         IF (allocated(a_water_canopyevap)) a_water_canopyevap(imt, ipatch) = &
	            a_water_canopyevap(imt, ipatch) + max(water_mass, 0._r8)
	      CASE (MT_EVAP_KIND_SUBL)
	         IF (allocated(a_mt_subl)) a_mt_subl(imt, ipatch) = &
	            a_mt_subl(imt, ipatch) + material_mass
	         IF (allocated(a_water_subl)) a_water_subl(imt, ipatch) = &
	            a_water_subl(imt, ipatch) + max(water_mass, 0._r8)
	      CASE (MT_EVAP_KIND_WETLAND)
	         IF (allocated(a_mt_wetland_evap)) a_mt_wetland_evap(imt, ipatch) = &
	            a_mt_wetland_evap(imt, ipatch) + material_mass
	         IF (allocated(a_water_wetland_evap)) a_water_wetland_evap(imt, ipatch) = &
	            a_water_wetland_evap(imt, ipatch) + max(water_mass, 0._r8)
	      END SELECT
	   END SUBROUTINE material_book_evap_loss

   !---------------------------------------------------------------
   ! Main-step runoff hand-off to the river routing.
   !
   ! Lake/waterbody patches are integrated in sub-steps of 1800 s inside
   ! one land step; each substep computes its own material runoff mass
   ! ``mt_rnof``.  The river coupling reads ``mt_rnof_step`` once per land
   ! step, so every substep must ADD its mass here.  The counter is cleared
   ! exactly once per main step (CoLM.F90 before CoLMDRIVER); it must not be
   ! cleared by the per-substep conservation snapshot or by the special-patch
   ! save/reconstruct routines.
   !---------------------------------------------------------------
   SUBROUTINE mt_add_rnof_step (imt, ipatch, mt_rnof)
      IMPLICIT NONE
      integer,  intent(in) :: imt, ipatch
      real(r8), intent(in) :: mt_rnof

      IF (nmts <= 0) RETURN
      IF (.not. allocated(mt_rnof_step)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN
      IF (.not. mt_uses_generic_water(imt)) RETURN
      IF (ipatch < 1 .or. ipatch > size(mt_rnof_step, 2)) RETURN
      mt_rnof_step(imt, ipatch) = mt_rnof_step(imt, ipatch) + max(mt_rnof, 0._r8)
   END SUBROUTINE mt_add_rnof_step

   SUBROUTINE mt_reset_rnof_step ()
      IMPLICIT NONE
      integer :: imt

      IF (nmts <= 0) RETURN
      IF (.not. allocated(mt_rnof_step)) RETURN
      DO imt = 1, nmts
         IF (mt_uses_generic_water(imt)) mt_rnof_step(imt, :) = 0._r8
      ENDDO
   END SUBROUTINE mt_reset_rnof_step

   !---------------------------------------------------------------
   ! Rebuild one material's prognostic storage from a caller-owned
   ! mixed-box ratio. This is used by patch branches that do not run
   ! the full material hydrology pipeline (glacier / waterbody) but can
   ! no longer assume a fixed R_init signature when runtime forcing or
   ! fractionation is enabled.
   !
   ! `snl` is the accounting lower bound used by mt_save_storage:
   ! layers below snl+1 are zeroed so hidden snow-layer material cannot
   ! leak outside the caller's conservation account.
   !---------------------------------------------------------------
   SUBROUTINE sync_mt_patch_ratio (imt, ipatch, snl, maxsnl, nl_soil, &
      wliq_soisno, wice_soisno, wa, wdsrf, scv, R_mix)
      IMPLICIT NONE
      integer,  intent(in) :: imt, ipatch, snl, maxsnl, nl_soil
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil)
      real(r8), intent(in) :: wa, wdsrf, scv, R_mix

      integer :: j

      IF (nmts <= 0) RETURN
      IF (.not. allocated(mt_wliq_soisno)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN

      IF (allocated(mt_ldew_rain)) mt_ldew_rain(imt, ipatch) = 0._r8
      IF (allocated(mt_ldew_snow)) mt_ldew_snow(imt, ipatch) = 0._r8

      DO j = maxsnl + 1, nl_soil
         IF (j >= snl + 1) THEN
            mt_wliq_soisno(imt, j, ipatch) = max(wliq_soisno(j), 0._r8) * R_mix
            mt_wice_soisno(imt, j, ipatch) = max(wice_soisno(j), 0._r8) * R_mix
            CALL mt_equilibrate_dissolved(imt, max(wliq_soisno(j), 0._r8), &
               mt_wliq_soisno(imt, j, ipatch), mt_solid_soisno(imt, j, ipatch))
         ELSE
            mt_wliq_soisno(imt, j, ipatch) = 0._r8
            mt_wice_soisno(imt, j, ipatch) = 0._r8
         ENDIF
      ENDDO

      IF (material_is_isotope(imt)) THEN
         mt_wa(imt, ipatch) = material_aquifer_actual_water( &
            wa, mt_aquifer_ref_water(ipatch)) * R_mix - mt_aquifer_ref_mass(imt, ipatch)
      ELSE
         mt_wa(imt, ipatch) = wa * R_mix
      ENDIF
      mt_wdsrf(imt, ipatch) = max(wdsrf, 0._r8) * R_mix
      CALL mt_equilibrate_dissolved(imt, wa, mt_wa(imt, ipatch), &
         mt_subsurface_solid(imt, ipatch))
      CALL mt_equilibrate_dissolved(imt, max(wdsrf, 0._r8), &
         mt_wdsrf(imt, ipatch), mt_surface_solid(imt, ipatch))

      IF (snl < 0) THEN
         mt_scv(imt, ipatch) = 0._r8
      ELSE
         mt_scv(imt, ipatch) = max(scv, 0._r8) * R_mix
      ENDIF
   END SUBROUTINE sync_mt_patch_ratio


END MODULE MOD_MT_Vars
#endif
