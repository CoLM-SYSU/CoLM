#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_Conservation

   USE MOD_Precision
   USE MOD_MT_Namelist, only: DEF_MT_BALANCE_ABORT_NBAD, DEF_MT_RESID_ABORT_NBAD
   USE MOD_MT_Defs, only: nmts, mt_tiny, material_init_water_ratio, &
      material_can_use_fixed_signature, mt_uses_generic_water
   USE MOD_MT_Solute, only: material_decay_fraction, mt_has_dissolved_limit
   USE MOD_MT_Frac, only: material_fractionation_active
   USE MOD_MT_Vars
#ifdef MATERIAL_SIGNED_SOIL_DEBT_EXPERIMENT
   USE MOD_SPMD_Task, only: CoLM_stop
#endif

   IMPLICIT NONE

   ! Hard balance tolerance for the water-corrected per-patch material residual.
   ! Keep ordinary land patches strict without using an unrealistically tiny
   ! absolute floor: double-precision budget sums over many storage/flux
   ! components can leave bounded roundoff residuals even when accounting is
   ! closed.  This floor is intentionally below the old 5e-5 value so stale-slot
   ! or missing-flux errors are not hidden by a broad absolute tolerance.
   ! Urban simple patches use the same threshold as ordinary land patches;
   ! their impermeable qgtop<0 evaporation branch is mirrored explicitly in
   ! MOD_MT_SoilWater instead of hidden behind a local tolerance waiver.
   real(r8), parameter :: mt_balance_abs_tol = 1.0e-12_r8
   real(r8), parameter :: mt_balance_rel_tol = 1.0e-10_r8
   integer, parameter :: n_storage_diag = 12
   integer, parameter :: n_flux_diag = 7

   ! Per-step accumulator snapshots (saved at start of each timestep)
   real(r8), allocatable, save :: snap_precip(:,:)  ! (nmts, numpatch)
   real(r8), allocatable, save :: snap_vapor_exchange(:,:)
   ! BGC nitrogen leaching input, kept out of the precipitation-signature check.
   real(r8), allocatable, save :: snap_bgc_leach(:,:)
   real(r8), allocatable, save :: snap_evap  (:,:)
   real(r8), allocatable, save :: snap_rnof  (:,:)
   real(r8), allocatable, save :: snap_rsur  (:,:)
   real(r8), allocatable, save :: snap_rsub  (:,:)
   real(r8), allocatable, save :: snap_qinfl (:,:)
   real(r8), allocatable, save :: snap_qcharge(:,:)
   real(r8), allocatable, save :: snap_storage_comp(:,:,:) ! (component, material, patch)

   ! Worker-local worst-patch tracker for mt_balance_check. Previously
   ! the routine only printed when `ipatch <= 1 .and. imt == 1`, so every
   ! other failing patch was invisible. These counters accumulate across
   ! patches within a step; caller drains them via material_balance_report
   ! at step end.
   real(r8), save :: balance_worst_err   = 0._r8
   real(r8), save :: balance_worst_diag(19) = 0._r8
   real(r8), save :: balance_worst_sbeg(n_storage_diag) = 0._r8
   real(r8), save :: balance_worst_send(n_storage_diag) = 0._r8
   real(r8), save :: balance_worst_sds (n_storage_diag) = 0._r8
   real(r8), save :: balance_worst_fcomp(n_flux_diag) = 0._r8
   integer,  save :: balance_worst_ipatch = 0
   integer,  save :: balance_worst_imt   = 0
   integer,  save :: balance_worst_ptype  = -1
   integer,  save :: balance_nbad         = 0
   ! Fixed-signature materials have a second, independent contract: each water-
   ! carrying boundary flux must equal the matching water flux times R_init.
   ! Keep this separate from mass conservation so neither check can hide the
   ! other by substituting theoretical fluxes for the amounts actually booked.
   real(r8), save :: signature_worst_abs      = 0._r8
   real(r8), save :: signature_worst_tol      = 0._r8
   real(r8), save :: signature_worst_terms(3) = 0._r8 ! input, evap, runoff
   integer,  save :: signature_worst_ipatch = 0
   integer,  save :: signature_worst_imt   = 0
   integer,  save :: signature_worst_ptype  = -1
   integer,  save :: signature_nbad         = 0
   ! P1 diagnostic: the per-layer mt_wliq<->wliq_soisno reconciliation
   ! (MOD_MT_SoilWater) injects/removes material booked as numerical_source_sink.
   ! This is an explicit solver-closure term, not an unbooked transport leak, so
   ! the hard transport-balance check subtracts it like decay_source_sink.
   ! It still has its own independent hard cap below; otherwise a large
   ! near-dry collapse could be "balanced" by the residual account forever.
   real(r8), parameter :: mt_resid_warn_frac = 1.0e-6_r8
   real(r8), parameter :: mt_resid_abs_tol   = 5.0e-5_r8
   real(r8), parameter :: mt_resid_rel_tol   = 1.0e-6_r8
   real(r8), save :: resid_worst_abs      = 0._r8
   real(r8), save :: resid_hard_worst_abs = 0._r8
   real(r8), save :: resid_hard_worst_tol = 0._r8
   integer,  save :: resid_hard_worst_ipatch = 0
   integer,  save :: resid_hard_worst_imt   = 0
   integer,  save :: resid_hard_worst_ptype  = -1
   integer,  save :: resid_nbad           = 0
   integer,  save :: resid_hard_nbad      = 0

   PUBLIC :: deallocate_mt_conservation
   PUBLIC :: mt_apply_decay
   PUBLIC :: material_balance_report

CONTAINS

   !---------------------------------------------------------------
   ! Release the save-level snap_* snapshots. Must be called from
   ! mt_land_final so that a subsequent mt_land_init on the same
   ! process (regression harness, embedded driver) does not re-use
   ! stale arrays sized to the previous numpatch — the lazy
   ! `IF (.not. allocated) allocate` guard in mt_save_storage
   ! only fires on the first call ever, so without this deallocate
   ! a changed numpatch would either truncate or overrun snap_*.
   !---------------------------------------------------------------
   SUBROUTINE deallocate_mt_conservation ()
      IMPLICIT NONE
      IF (allocated(snap_precip)) deallocate(snap_precip)
      IF (allocated(snap_vapor_exchange)) deallocate(snap_vapor_exchange)
      IF (allocated(snap_bgc_leach)) deallocate(snap_bgc_leach)
      IF (allocated(snap_evap  )) deallocate(snap_evap  )
      IF (allocated(snap_rnof  )) deallocate(snap_rnof  )
      IF (allocated(snap_rsur  )) deallocate(snap_rsur  )
      IF (allocated(snap_rsub  )) deallocate(snap_rsub  )
      IF (allocated(snap_qinfl )) deallocate(snap_qinfl )
      IF (allocated(snap_qcharge)) deallocate(snap_qcharge)
      IF (allocated(snap_storage_comp)) deallocate(snap_storage_comp)
   END SUBROUTINE deallocate_mt_conservation

   SUBROUTINE mt_save_storage (ipatch, snl, nl_soil, waterstorage)
      IMPLICIT NONE
      integer, intent(in) :: ipatch, snl, nl_soil
      ! Optional: patch-level irrigation reservoir (mm of water equivalent).
      ! Present only when DEF_USE_IRRIGATION is on under CROP. When passed,
      ! mt_waterstorage is re-synced to waterstorage*R_init (Phase-1
      ! invariant) and included in the storage balance; this mirrors the
      ! water-side convention at CoLMMAIN:806 that treats waterstorage as
      ! part of the column inventory (so irrigation becomes an internal
      ! transfer instead of an external atmospheric input).
      real(r8), intent(in), optional :: waterstorage
         integer :: imt, j, lb_store
      real(r8) :: R_init
      real(r8) :: storage_comp(n_storage_diag)
      logical  :: fixed_signature_storage

         IF (nmts <= 0) RETURN
#ifdef MATERIAL_SIGNED_SOIL_DEBT_EXPERIMENT
         IF (any(mt_wliq_soisno(:,1:nl_soil,ipatch) < 0._r8)) &
            CALL CoLM_stop('signed soil debt: rejected full timestep; direct-kernel experiment only')
#endif
         lb_store = max(lbound(mt_wliq_soisno, 2), snl + 1)

      ! Allocate snapshots on first call
      IF (.not. allocated(snap_precip)) THEN
         allocate(snap_precip(nmts, size(mt_storage_beg,2))); snap_precip = 0._r8
         allocate(snap_vapor_exchange(nmts, size(mt_storage_beg,2)))
         snap_vapor_exchange = 0._r8
         allocate(snap_bgc_leach(nmts, size(mt_storage_beg,2)))
         snap_bgc_leach = 0._r8
         allocate(snap_evap  (nmts, size(mt_storage_beg,2))); snap_evap   = 0._r8
         allocate(snap_rnof  (nmts, size(mt_storage_beg,2))); snap_rnof   = 0._r8
         allocate(snap_rsur  (nmts, size(mt_storage_beg,2))); snap_rsur   = 0._r8
         allocate(snap_rsub  (nmts, size(mt_storage_beg,2))); snap_rsub   = 0._r8
         allocate(snap_qinfl (nmts, size(mt_storage_beg,2))); snap_qinfl  = 0._r8
         allocate(snap_qcharge(nmts, size(mt_storage_beg,2))); snap_qcharge = 0._r8
         allocate(snap_storage_comp(n_storage_diag, nmts, size(mt_storage_beg,2)))
         snap_storage_comp = 0._r8
      ENDIF

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
            ! mt_rnof_step is the MAIN-STEP river hand-off accumulator: it
            ! is reset once before CoLMDRIVER, not by the per-substep
            ! conservation snapshot.
            IF (allocated(mt_decay_source_step)) mt_decay_source_step(imt, ipatch) = 0._r8
            IF (allocated(mt_numerical_residual_step)) mt_numerical_residual_step(imt, ipatch) = 0._r8
            IF (allocated(mt_numerical_water_step)) mt_numerical_water_step(imt, ipatch) = 0._r8

         ! Phase-1 re-sync of the irrigation reservoir material to current
         ! waterstorage. Under no-fractionation tests all refills arrive at
         ! R_init so the ratio is constant. When fractionation is active, do
         ! not overwrite the reservoir isotope state with the initial delta.
         fixed_signature_storage = material_can_use_fixed_signature(imt) .and. &
            .not. material_fractionation_active(imt) .and. &
            .not. mt_has_dissolved_limit(imt)
         IF (allocated(mt_runtime_forced)) THEN
            fixed_signature_storage = fixed_signature_storage .and. .not. mt_runtime_forced(imt)
         ENDIF
         IF (present(waterstorage) .and. allocated(mt_waterstorage) .and. &
             fixed_signature_storage) THEN
            R_init = material_init_water_ratio(imt)
            mt_waterstorage(imt, ipatch) = max(waterstorage, 0._r8) * R_init
         ENDIF

         ! Save current storage split by component:
         ! 1 ldew, 2 active snow/soil liquid, 3 active snow/soil ice,
         ! 4 wa, 5 wdsrf, 6 wetwat, 7 scv, 8 waterstorage,
         ! 9 leaf NSS isostorage.  Use the active lower bound (snl+1)
         ! so stale inactive snow slots left by snow layer topology changes
         ! are not double-counted after they have been merged into soil.
         storage_comp = 0._r8
         storage_comp(1) = mt_ldew_rain(imt, ipatch) + mt_ldew_snow(imt, ipatch)
            DO j = lb_store, nl_soil
               storage_comp(2) = storage_comp(2) + mt_wliq_soisno(imt, j, ipatch)
               storage_comp(3) = storage_comp(3) + mt_wice_soisno(imt, j, ipatch)
            ENDDO
         ! The aquifer reference isotope mass is fixed for this patch. Keep
         ! the relative component here: Mref cancels from the step's delta,
         ! matching hydrology's wa-only budget (the fixed Vref also cancels).
         storage_comp(4) = mt_wa(imt, ipatch)
         storage_comp(5) = mt_wdsrf(imt, ipatch)
         storage_comp(6) = mt_wetwat(imt, ipatch)
         storage_comp(7) = mt_scv(imt, ipatch)
         IF (allocated(mt_waterstorage)) THEN
            storage_comp(8) = mt_waterstorage(imt, ipatch)
         ENDIF
         IF (allocated(mt_leaf_iso_storage)) THEN
            storage_comp(9) = mt_leaf_iso_storage(imt, ipatch)
         ENDIF
         IF (allocated(mt_surface_residue)) THEN
            storage_comp(10) = mt_surface_residue(imt, ipatch)
         ENDIF
         IF (allocated(mt_subsurface_residue)) THEN
            storage_comp(11) = mt_subsurface_residue(imt, ipatch)
         ENDIF
         IF (allocated(mt_solid_soisno)) THEN
            DO j = lb_store, nl_soil
               storage_comp(12) = storage_comp(12) + mt_solid_soisno(imt, j, ipatch)
            ENDDO
         ENDIF
         IF (allocated(mt_canopy_solid)) &
            storage_comp(12) = storage_comp(12) + mt_canopy_solid(imt, ipatch)
         IF (allocated(mt_surface_solid)) &
            storage_comp(12) = storage_comp(12) + mt_surface_solid(imt, ipatch)
         IF (allocated(mt_subsurface_solid)) &
            storage_comp(12) = storage_comp(12) + mt_subsurface_solid(imt, ipatch)
         IF (allocated(mt_waterstorage_solid)) &
            storage_comp(12) = storage_comp(12) + mt_waterstorage_solid(imt, ipatch)
         snap_storage_comp(:, imt, ipatch) = storage_comp
         mt_storage_beg(imt, ipatch) = sum(storage_comp)

         ! Snapshot accumulators at start of this step
         snap_precip(imt, ipatch) = a_mt_precip(imt, ipatch)
         IF (allocated(a_mt_vapor_exchange)) &
            snap_vapor_exchange(imt, ipatch) = a_mt_vapor_exchange(imt, ipatch)
         IF (allocated(a_mt_bgc_leach)) &
            snap_bgc_leach(imt, ipatch) = a_mt_bgc_leach(imt, ipatch)
         snap_evap  (imt, ipatch) = a_mt_evap  (imt, ipatch)
         snap_rnof  (imt, ipatch) = a_mt_rnof  (imt, ipatch)
         snap_rsur  (imt, ipatch) = a_mt_rsur  (imt, ipatch)
         snap_rsub  (imt, ipatch) = a_mt_rsub  (imt, ipatch)
         snap_qinfl (imt, ipatch) = a_mt_qinfl (imt, ipatch)
         snap_qcharge(imt, ipatch) = a_mt_qcharge(imt, ipatch)
      ENDDO
   END SUBROUTINE mt_save_storage

   ! Contract: this routine applies only first-order decay of the generic land-water pools
   ! that live in MOD_MT_Vars and share the ordinary water-transport
   ! storage/flux budget.  species-owned reacting state (for example CH4's
   ! conc_methane* prognostic pools) must keep its own source/sink accounting
   ! and restart callbacks; it must not be double-booked here through
   ! mt_decay_source_step.
   SUBROUTINE mt_apply_decay (ipatch, snl, nl_soil, deltim)
      IMPLICIT NONE
      integer,  intent(in) :: ipatch, snl, nl_soil
      real(r8), intent(in) :: deltim
         integer :: imt, j, lb_store
      real(r8) :: decay_fraction, source_sink

         IF (nmts <= 0) RETURN
         IF (.not. allocated(mt_decay_source_step)) RETURN
         lb_store = max(lbound(mt_wliq_soisno, 2), snl + 1)

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         decay_fraction = material_decay_fraction(imt, deltim)
         IF (decay_fraction <= 0._r8) CYCLE

         source_sink = 0._r8
         CALL decay_pool(mt_ldew_rain(imt, ipatch), decay_fraction, source_sink)
         CALL decay_pool(mt_ldew_snow(imt, ipatch), decay_fraction, source_sink)
            DO j = lb_store, nl_soil
               CALL decay_pool(mt_wliq_soisno(imt, j, ipatch), decay_fraction, source_sink)
               CALL decay_pool(mt_wice_soisno(imt, j, ipatch), decay_fraction, source_sink)
            ENDDO
         CALL decay_pool(mt_wa(imt, ipatch), decay_fraction, source_sink)
         CALL decay_pool(mt_wdsrf(imt, ipatch), decay_fraction, source_sink)
         CALL decay_pool(mt_wetwat(imt, ipatch), decay_fraction, source_sink)
         CALL decay_pool(mt_scv(imt, ipatch), decay_fraction, source_sink)
         IF (allocated(mt_waterstorage)) THEN
            CALL decay_pool(mt_waterstorage(imt, ipatch), decay_fraction, source_sink)
         ENDIF
         IF (allocated(mt_leaf_iso_storage)) THEN
            CALL decay_pool(mt_leaf_iso_storage(imt, ipatch), decay_fraction, source_sink)
         ENDIF
         IF (allocated(mt_surface_residue)) THEN
            CALL decay_pool(mt_surface_residue(imt, ipatch), decay_fraction, source_sink)
         ENDIF
         IF (allocated(mt_subsurface_residue)) THEN
            CALL decay_pool(mt_subsurface_residue(imt, ipatch), decay_fraction, source_sink)
         ENDIF
         IF (allocated(mt_solid_soisno)) THEN
            DO j = lb_store, nl_soil
               CALL decay_pool(mt_solid_soisno(imt, j, ipatch), decay_fraction, source_sink)
            ENDDO
         ENDIF
         IF (allocated(mt_canopy_solid)) &
            CALL decay_pool(mt_canopy_solid(imt, ipatch), decay_fraction, source_sink)
         IF (allocated(mt_surface_solid)) &
            CALL decay_pool(mt_surface_solid(imt, ipatch), decay_fraction, source_sink)
         IF (allocated(mt_subsurface_solid)) &
            CALL decay_pool(mt_subsurface_solid(imt, ipatch), decay_fraction, source_sink)
         IF (allocated(mt_waterstorage_solid)) &
            CALL decay_pool(mt_waterstorage_solid(imt, ipatch), decay_fraction, source_sink)

         mt_decay_source_step(imt, ipatch) = mt_decay_source_step(imt, ipatch) &
            + source_sink
      ENDDO

   CONTAINS

      SUBROUTINE decay_pool (pool, fraction, source_sink)
         real(r8), intent(inout) :: pool
         real(r8), intent(in)    :: fraction
         real(r8), intent(inout) :: source_sink
         real(r8) :: before

         ! Negative material pools can occur only in signed water-debt states
         ! such as wa. Do not turn those numerical debts into chemical source.
         IF (pool <= mt_tiny) RETURN
         before = pool
         pool = pool * (1._r8 - fraction)
         source_sink = source_sink + pool - before
      END SUBROUTINE decay_pool

   END SUBROUTINE mt_apply_decay

   SUBROUTINE mt_balance_check (ipatch, snl, nl_soil, deltim, xerr_material, &
                                    patchtype_in, water_err_in, water_dS_in, &
                                    water_input_in, water_output_in, &
                                    water_evap_in, water_rnof_in, flood_heterogeneous_in)
      IMPLICIT NONE
      integer,  intent(in)  :: ipatch, snl, nl_soil
      real(r8), intent(in)  :: deltim
      real(r8), intent(out) :: xerr_material
      ! Optional patchtype passthrough so MT_BAL report can identify the
      ! offending branch (0=soil, 1=urban, 2=wetland, 3=glacier, 4=waterbody).
      ! Kept optional so old call sites (if any) still compile.
      integer,  intent(in), optional :: patchtype_in
      real(r8), intent(in), optional :: water_err_in
      real(r8), intent(in), optional :: water_dS_in
      real(r8), intent(in), optional :: water_input_in
      real(r8), intent(in), optional :: water_output_in
      real(r8), intent(in), optional :: water_evap_in
      real(r8), intent(in), optional :: water_rnof_in
      logical, intent(in), optional :: flood_heterogeneous_in

      integer  :: imt, j, lb_store
      real(r8) :: storage_end, step_input, step_evap, step_rnof, step_output, err
      real(r8) :: step_input_check, step_output_check
      real(r8) :: step_rsur, step_rsub, step_qinfl, step_qcharge
      real(r8) :: step_vapor_exchange
      real(r8) :: step_bgc_leach
      real(r8) :: R_init, water_err, water_err_R, err_minus_water, check_err, booked_host_water
      real(r8) :: decay_source_sink, numerical_source_sink
      real(r8) :: water_dS, water_input, water_output, water_evap, water_rnof
      real(r8) :: dS_minus_water_R, in_minus_water_R, out_minus_water_R
      real(r8) :: evap_minus_water_R, rnof_minus_water_R
      real(r8) :: balance_scale, balance_tol, resid_scale, resid_tol
      real(r8) :: signature_error, signature_scale, signature_tol
      real(r8) :: storage_comp_end(n_storage_diag)
      real(r8) :: storage_comp_beg(n_storage_diag)
      real(r8) :: storage_comp_dS (n_storage_diag)
      logical  :: fixed_signature_step, water_corrected_check

      xerr_material = 0._r8
      IF (nmts <= 0) RETURN
      lb_store = max(lbound(mt_wliq_soisno, 2), snl + 1)

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         ! Current storage split by component; keep this in the same
         ! component order as mt_save_storage.
         storage_comp_end = 0._r8
         storage_comp_end(1) = mt_ldew_rain(imt, ipatch) + mt_ldew_snow(imt, ipatch)
         DO j = lb_store, nl_soil
            storage_comp_end(2) = storage_comp_end(2) + mt_wliq_soisno(imt, j, ipatch)
            storage_comp_end(3) = storage_comp_end(3) + mt_wice_soisno(imt, j, ipatch)
         ENDDO
         storage_comp_end(4) = mt_wa(imt, ipatch)
         storage_comp_end(5) = mt_wdsrf(imt, ipatch)
         storage_comp_end(6) = mt_wetwat(imt, ipatch)
         storage_comp_end(7) = mt_scv(imt, ipatch)
         IF (allocated(mt_waterstorage)) THEN
            storage_comp_end(8) = mt_waterstorage(imt, ipatch)
         ENDIF
         IF (allocated(mt_leaf_iso_storage)) THEN
            storage_comp_end(9) = mt_leaf_iso_storage(imt, ipatch)
         ENDIF
         IF (allocated(mt_surface_residue)) THEN
            storage_comp_end(10) = mt_surface_residue(imt, ipatch)
         ENDIF
         IF (allocated(mt_subsurface_residue)) THEN
            storage_comp_end(11) = mt_subsurface_residue(imt, ipatch)
         ENDIF
         IF (allocated(mt_solid_soisno)) THEN
            DO j = lb_store, nl_soil
               storage_comp_end(12) = storage_comp_end(12) + mt_solid_soisno(imt, j, ipatch)
            ENDDO
         ENDIF
         IF (allocated(mt_canopy_solid)) &
            storage_comp_end(12) = storage_comp_end(12) + mt_canopy_solid(imt, ipatch)
         IF (allocated(mt_surface_solid)) &
            storage_comp_end(12) = storage_comp_end(12) + mt_surface_solid(imt, ipatch)
         IF (allocated(mt_subsurface_solid)) &
            storage_comp_end(12) = storage_comp_end(12) + mt_subsurface_solid(imt, ipatch)
         IF (allocated(mt_waterstorage_solid)) &
            storage_comp_end(12) = storage_comp_end(12) + mt_waterstorage_solid(imt, ipatch)
         storage_end = sum(storage_comp_end)
         storage_comp_beg = snap_storage_comp(:, imt, ipatch)
         storage_comp_dS = storage_comp_end - storage_comp_beg

         R_init = material_init_water_ratio(imt)
         water_corrected_check = material_can_use_fixed_signature(imt)
         fixed_signature_step = water_corrected_check .and. .not. material_fractionation_active(imt)
         IF (allocated(mt_runtime_forced)) THEN
            fixed_signature_step = fixed_signature_step .and. .not. mt_runtime_forced(imt)
         ENDIF
         ! Floodwater can carry a river signature different from R_init even
         ! when precipitation is fixed-signature. That signature persists in
         ! soil after a wet step and across restart. With feedback enabled,
         ! retain raw conservation but do not use R_init as a soil-patch
         ! flux or host-closure correction on later dry steps either.
         IF (present(flood_heterogeneous_in)) THEN
            IF (flood_heterogeneous_in) THEN
               fixed_signature_step = .false.
               ! The raw material budget already subtracts its explicitly
               ! booked numerical source/sink. Applying R_init times the
               ! host water closure as well would double-correct a flooded
               ! column whose incoming isotope ratio is not R_init.
               water_corrected_check = .false.
            ENDIF
         ENDIF

         ! Per-step fluxes = current accumulator - snapshot at step start.
         ! Conservation always uses the amounts actually booked. Fixed-
         ! signature agreement with water_flux*R_init is checked independently
         ! below; substituting theoretical values here would hide accounting
         ! errors in precipitation or evaporation.
         step_input  = a_mt_precip(imt, ipatch) - snap_precip(imt, ipatch)
         step_evap   = a_mt_evap(imt, ipatch) - snap_evap(imt, ipatch)
         step_rsur   = a_mt_rsur(imt, ipatch) - snap_rsur(imt, ipatch)
         step_rsub   = a_mt_rsub(imt, ipatch) - snap_rsub(imt, ipatch)
         step_qinfl  = a_mt_qinfl(imt, ipatch) - snap_qinfl(imt, ipatch)
         step_qcharge = a_mt_qcharge(imt, ipatch) - snap_qcharge(imt, ipatch)
         step_rnof   = a_mt_rnof(imt, ipatch) - snap_rnof(imt, ipatch)
         step_output = step_evap + step_rnof
         step_input_check = step_input
         step_output_check = step_output

         ! Two-way equilibrium exchange with atmospheric vapour moves material
         ! across the patch boundary with NO accompanying water flux, so it
         ! is included in mass conservation but excluded from the separate
         ! fixed-signature check. Signed: positive = net uptake.
         step_vapor_exchange = 0._r8
         IF (allocated(a_mt_vapor_exchange)) THEN
            step_vapor_exchange = a_mt_vapor_exchange(imt, ipatch) &
                                - snap_vapor_exchange(imt, ipatch)
         ENDIF
         step_input_check = step_input_check + step_vapor_exchange

         ! Nitrogen leached in BGC enters the material budget as an external
         ! input (it left the BGC mineral-N pool, not the atmosphere).  The
         ! injection happens after this step's check (CoLMDRIVER, after BGC),
         ! so the next step's snapshot already holds it together with the
         ! matching a_mt_rnof output; the two cancel and the balance stays
         ! exact.  It is counted in mass conservation but kept out of the
         ! precipitation signature check for the same reason as the vapour
         ! exchange above.
         step_bgc_leach = 0._r8
         IF (allocated(a_mt_bgc_leach)) &
            step_bgc_leach = a_mt_bgc_leach(imt, ipatch) - snap_bgc_leach(imt, ipatch)
         step_input_check = step_input_check + step_bgc_leach

         ! Conservation: Δstorage = input - output + source_sink.
         ! Decay and numerical source/sink terms are explicit state mutations
         ! booked by their owning processes.  Subtract both from the hard check
         ! so solver-closure residuals are reported separately instead of being
         ! misclassified as unbooked transport leaks.
         decay_source_sink = 0._r8
         IF (allocated(mt_decay_source_step)) THEN
            decay_source_sink = mt_decay_source_step(imt, ipatch)
         ENDIF
         numerical_source_sink = 0._r8
         IF (allocated(mt_numerical_residual_step)) THEN
            numerical_source_sink = mt_numerical_residual_step(imt, ipatch)
         ENDIF
         err = storage_end - mt_storage_beg(imt, ipatch) - step_input_check + step_output_check
         mt_balance_err(imt, ipatch) = err - decay_source_sink - numerical_source_sink
         IF (present(water_err_in)) THEN
            water_err = water_err_in
         ELSE
            water_err = 0._r8
         ENDIF
         booked_host_water = 0._r8
         IF (allocated(mt_numerical_water_step)) THEN
            IF (water_err * mt_numerical_water_step(imt, ipatch) > 0._r8) &
               booked_host_water = sign(min(abs(water_err), &
                  abs(mt_numerical_water_step(imt, ipatch))), water_err)
         ENDIF
         water_err_R = (water_err - booked_host_water) * R_init
         err_minus_water = err - water_err_R
         ! The hard material check should ignore host water-budget non-closure
         ! for isotope materials. Fractionation and runtime atmospheric forcing
         ! keep their process-owned flux signatures above, but they should not
         ! make the hard material check count host water non-closure as material
         ! non-conservation. A real material accounting bug remains visible in
         ! err_minus_water.
         IF (present(water_err_in) .and. water_corrected_check) THEN
            check_err = err_minus_water
         ELSE
            check_err = err
         ENDIF
         check_err = check_err - decay_source_sink - numerical_source_sink
         IF (present(water_dS_in)) THEN
            water_dS = water_dS_in
         ELSE
            water_dS = 0._r8
         ENDIF
         IF (present(water_input_in)) THEN
            water_input = water_input_in
         ELSE
            water_input = 0._r8
         ENDIF
         IF (present(water_output_in)) THEN
            water_output = water_output_in
         ELSE
            water_output = 0._r8
         ENDIF
         IF (present(water_evap_in)) THEN
            water_evap = water_evap_in
         ELSE
            water_evap = 0._r8
         ENDIF
         IF (present(water_rnof_in)) THEN
            water_rnof = water_rnof_in
         ELSE
            water_rnof = 0._r8
         ENDIF
         a_water_precip(imt, ipatch) = a_water_precip(imt, ipatch) + max(water_input, 0._r8)
         a_water_rnof(imt, ipatch) = a_water_rnof(imt, ipatch) + max(water_rnof, 0._r8)
         dS_minus_water_R  = (storage_end - mt_storage_beg(imt, ipatch)) &
                           - water_dS * R_init
         in_minus_water_R  = step_input  - water_input  * R_init
         out_minus_water_R = step_output - water_output * R_init
         evap_minus_water_R = step_evap - water_evap * R_init
         rnof_minus_water_R = step_rnof - water_rnof * R_init

         balance_scale = max(abs(storage_end), abs(mt_storage_beg(imt, ipatch)), &
            abs(step_input_check), abs(step_output_check), &
            abs(decay_source_sink), abs(numerical_source_sink))
         balance_tol = mt_balance_abs_tol + mt_balance_rel_tol * balance_scale
         resid_scale = max(1._r8, abs(storage_end), abs(mt_storage_beg(imt, ipatch)), &
            abs(step_input_check), abs(step_output_check), abs(decay_source_sink))
         resid_tol = max(mt_resid_abs_tol, mt_resid_rel_tol * resid_scale)

         signature_error = 0._r8
         signature_scale = 0._r8
         IF (fixed_signature_step) THEN
            IF (present(water_input_in) .or. present(water_evap_in)) THEN
               signature_error = max(signature_error, abs(in_minus_water_R - evap_minus_water_R))
               signature_scale = max(signature_scale, abs(step_input), abs(water_input * R_init), &
                  abs(step_evap), abs(water_evap * R_init))
            ENDIF
            IF (present(water_rnof_in)) THEN
               signature_error = max(signature_error, abs(rnof_minus_water_R))
               signature_scale = max(signature_scale, abs(step_rnof), abs(water_rnof * R_init))
            ENDIF
         ENDIF
         signature_tol = mt_balance_abs_tol + mt_balance_rel_tol * signature_scale
         IF (fixed_signature_step .and. signature_error > signature_tol) THEN
            signature_nbad = signature_nbad + 1
            IF (signature_error > signature_worst_abs) THEN
               signature_worst_abs = signature_error
               signature_worst_tol = signature_tol
               signature_worst_terms = (/ in_minus_water_R, evap_minus_water_R, rnof_minus_water_R /)
               signature_worst_ipatch = ipatch
               signature_worst_imt = imt
               IF (present(patchtype_in)) THEN
                  signature_worst_ptype = patchtype_in
               ELSE
                  signature_worst_ptype = -1
               ENDIF
            ENDIF
         ENDIF

         IF (abs(check_err) > balance_tol) THEN
            ! Track worst-ever failure within this step (worker-local).
            ! material_balance_report surfaces the accumulated bulk at end of step.
            IF (abs(check_err) > abs(balance_worst_err)) THEN
               balance_worst_err    = check_err
               balance_worst_diag(1) = err
               balance_worst_diag(2) = storage_end - mt_storage_beg(imt, ipatch)
               balance_worst_diag(3) = step_input
               balance_worst_diag(4) = step_output
               balance_worst_diag(5) = water_err
               balance_worst_diag(6) = water_err_R
               balance_worst_diag(7) = err_minus_water
               balance_worst_diag(8) = water_dS
               balance_worst_diag(9) = water_input
               balance_worst_diag(10) = water_output
               balance_worst_diag(11) = dS_minus_water_R
               balance_worst_diag(12) = in_minus_water_R
               balance_worst_diag(13) = out_minus_water_R
               balance_worst_diag(14) = step_evap
               balance_worst_diag(15) = step_rnof
               balance_worst_diag(16) = water_evap
               balance_worst_diag(17) = water_rnof
               balance_worst_diag(18) = evap_minus_water_R
               balance_worst_diag(19) = rnof_minus_water_R
               balance_worst_sbeg  = storage_comp_beg
               balance_worst_send  = storage_comp_end
               balance_worst_sds   = storage_comp_dS
               balance_worst_fcomp = (/ step_input, step_evap, step_rsur, step_rsub, &
                  step_rnof, step_qinfl, step_qcharge /)
               balance_worst_ipatch = ipatch
               balance_worst_imt   = imt
               IF (present(patchtype_in)) THEN
                  balance_worst_ptype = patchtype_in
               ELSE
                  balance_worst_ptype = -1
               ENDIF
            ENDIF
            balance_nbad = balance_nbad + 1
         ENDIF

         ! P1: surface a non-trivial numerical residual. This is context only;
         ! the residual is explicitly booked and excluded from the hard tol.
         IF (abs(numerical_source_sink) > mt_resid_warn_frac * &
                max(abs(storage_end), abs(mt_storage_beg(imt, ipatch)))) THEN
            resid_nbad = resid_nbad + 1
            resid_worst_abs = max(resid_worst_abs, abs(numerical_source_sink))
         ENDIF
         IF (abs(numerical_source_sink) > resid_tol) THEN
            resid_hard_nbad = resid_hard_nbad + 1
            IF (abs(numerical_source_sink) > resid_hard_worst_abs) THEN
               resid_hard_worst_abs = abs(numerical_source_sink)
               resid_hard_worst_tol = resid_tol
               resid_hard_worst_ipatch = ipatch
               resid_hard_worst_imt   = imt
               IF (present(patchtype_in)) THEN
                  resid_hard_worst_ptype = patchtype_in
               ELSE
                  resid_hard_worst_ptype = -1
               ENDIF
            ENDIF
         ENDIF

         ! deltim is intent(in) and never checked. Every Makeoptions builds
         ! with -ffpe-trap=zero, so a caller passing 0 would not produce a
         ! harmless Inf in a diagnostic -- it would abort the model inside the
         ! conservation check, pointing the backtrace at the wrong place
         ! entirely. Report the raw residual instead; a rate is meaningless
         ! without a timestep anyway.
         IF (deltim > 0._r8) THEN
            xerr_material = max(xerr_material, abs(check_err) / deltim)
         ELSE
            xerr_material = max(xerr_material, abs(check_err))
         ENDIF
      ENDDO
   END SUBROUTINE mt_balance_check

   !---------------------------------------------------------------
   ! Drain worker-local balance counters. Call once per CoLM step
   ! after all patches have run mt_balance_check. Reduces the
   ! worst-err across workers so a single-patch failure elsewhere
   ! still lands in the per-rank log. Zeros the state so next step
   ! starts fresh.
   !
   ! This runs from inside CoLMDRIVER (worker-only gated at
   ! CoLM.F90:528). Reducing over p_comm_glb here would hang — master
   ! and IO never enter. Use p_comm_worker so the collective stays on
   ! the same ranks that produced the counts. Worker rank 0 prints the
   ! summary (not p_is_master, which lives on a different comm).
   !---------------------------------------------------------------
   SUBROUTINE material_balance_report ()
      USE MOD_SPMD_Task, only: CoLM_stop
#ifdef USEMPI
      USE MOD_SPMD_Task, only: p_comm_worker, p_iam_worker, p_err, p_is_worker
#endif
      IMPLICIT NONE
      real(r8) :: worst_abs, reduced_abs
      integer  :: nbad_total
      real(r8) :: resid_abs_total
      real(r8) :: resid_hard_abs_total, resid_hard_tol_total
      integer  :: resid_nbad_total
      integer  :: resid_hard_nbad_total
      real(r8) :: signature_abs_total, signature_tol_total
      real(r8) :: signature_terms_total(3)
      integer  :: signature_nbad_total
      integer  :: winner_rank
      integer  :: resid_winner_rank
      integer  :: signature_winner_rank
      integer  :: resid_hard_ipatch_total, resid_hard_imt_total, resid_hard_ptype_total
      integer  :: signature_ipatch_total, signature_imt_total, signature_ptype_total
      logical  :: print_me
      real(r8) :: maxima_local(4), maxima_total(4)
      integer  :: counts_local(4), counts_total(4)
      ! MAXLOC reduction: pair |err| with rank so we can ship ipatch/imt
      ! from the rank that owns the worst patch. Standard MPI_MAXLOC on
      ! real+int works with MPI_2DOUBLE_PRECISION wrapping (rank as real).
      real(r8) :: in_pair(2), out_pair(2)
#ifdef USEMPI
      include 'mpif.h'
#endif

      IF (nmts <= 0) RETURN

      worst_abs = abs(balance_worst_err)
      nbad_total  = balance_nbad
      resid_abs_total  = resid_worst_abs
      resid_nbad_total = resid_nbad
      resid_hard_abs_total  = resid_hard_worst_abs
      resid_hard_tol_total  = resid_hard_worst_tol
      resid_hard_nbad_total = resid_hard_nbad
      resid_hard_ipatch_total = resid_hard_worst_ipatch
      resid_hard_imt_total   = resid_hard_worst_imt
      resid_hard_ptype_total  = resid_hard_worst_ptype
      signature_abs_total     = signature_worst_abs
      signature_tol_total     = signature_worst_tol
      signature_terms_total   = signature_worst_terms
      signature_nbad_total    = signature_nbad
      signature_ipatch_total  = signature_worst_ipatch
      signature_imt_total    = signature_worst_imt
      signature_ptype_total   = signature_worst_ptype
      winner_rank = 0
      resid_winner_rank = 0
      signature_winner_rank = 0

#ifdef USEMPI
      reduced_abs = worst_abs
      IF (p_is_worker) THEN
         maxima_local = (/ worst_abs, resid_worst_abs, resid_hard_worst_abs, signature_worst_abs /)
         maxima_total = maxima_local
         counts_local = (/ balance_nbad, resid_nbad, resid_hard_nbad, signature_nbad /)
         counts_total = counts_local
         CALL mpi_reduce(maxima_local, maxima_total, 4, MPI_REAL8, MPI_MAX, 0, &
                         p_comm_worker, p_err)
         CALL mpi_allreduce(counts_local, counts_total, 4, MPI_INTEGER, MPI_SUM, &
                            p_comm_worker, p_err)

         reduced_abs          = maxima_total(1)
         resid_abs_total      = maxima_total(2)
         resid_hard_abs_total = maxima_total(3)
         signature_abs_total  = maxima_total(4)
         nbad_total            = counts_total(1)
         resid_nbad_total      = counts_total(2)
         resid_hard_nbad_total = counts_total(3)
         signature_nbad_total  = counts_total(4)

         IF (resid_hard_nbad_total > 0) THEN
            in_pair = (/ resid_hard_worst_abs, real(p_iam_worker, r8) /)
            CALL mpi_allreduce(in_pair, out_pair, 1, MPI_2DOUBLE_PRECISION, &
                               MPI_MAXLOC, p_comm_worker, p_err)
            resid_winner_rank = nint(out_pair(2))
            IF (resid_winner_rank /= 0) THEN
               IF (p_iam_worker == resid_winner_rank) THEN
                  CALL mpi_send(resid_hard_worst_ipatch, 1, MPI_INTEGER, 0, 101, p_comm_worker, p_err)
                  CALL mpi_send(resid_hard_worst_imt,   1, MPI_INTEGER, 0, 102, p_comm_worker, p_err)
                  CALL mpi_send(resid_hard_worst_ptype,  1, MPI_INTEGER, 0, 103, p_comm_worker, p_err)
                  CALL mpi_send(resid_hard_worst_tol,    1, MPI_REAL8,   0, 104, p_comm_worker, p_err)
               ELSEIF (p_iam_worker == 0) THEN
                  CALL mpi_recv(resid_hard_ipatch_total, 1, MPI_INTEGER, resid_winner_rank, 101, &
                                p_comm_worker, MPI_STATUS_IGNORE, p_err)
                  CALL mpi_recv(resid_hard_imt_total,   1, MPI_INTEGER, resid_winner_rank, 102, &
                                p_comm_worker, MPI_STATUS_IGNORE, p_err)
                  CALL mpi_recv(resid_hard_ptype_total,  1, MPI_INTEGER, resid_winner_rank, 103, &
                                p_comm_worker, MPI_STATUS_IGNORE, p_err)
                  CALL mpi_recv(resid_hard_tol_total,    1, MPI_REAL8,   resid_winner_rank, 104, &
                                p_comm_worker, MPI_STATUS_IGNORE, p_err)
               ENDIF
            ENDIF
         ENDIF

         IF (signature_nbad_total > 0) THEN
            in_pair = (/ signature_worst_abs, real(p_iam_worker, r8) /)
            CALL mpi_allreduce(in_pair, out_pair, 1, MPI_2DOUBLE_PRECISION, &
                               MPI_MAXLOC, p_comm_worker, p_err)
            signature_winner_rank = nint(out_pair(2))
            IF (signature_winner_rank /= 0) THEN
               IF (p_iam_worker == signature_winner_rank) THEN
                  CALL mpi_send(signature_worst_ipatch, 1, MPI_INTEGER, 0, 105, p_comm_worker, p_err)
                  CALL mpi_send(signature_worst_imt,   1, MPI_INTEGER, 0, 106, p_comm_worker, p_err)
                  CALL mpi_send(signature_worst_ptype,  1, MPI_INTEGER, 0, 107, p_comm_worker, p_err)
                  CALL mpi_send(signature_worst_tol,    1, MPI_REAL8,   0, 108, p_comm_worker, p_err)
                  CALL mpi_send(signature_worst_terms,  3, MPI_REAL8,   0, 109, p_comm_worker, p_err)
               ELSEIF (p_iam_worker == 0) THEN
                  CALL mpi_recv(signature_ipatch_total, 1, MPI_INTEGER, signature_winner_rank, 105, &
                                p_comm_worker, MPI_STATUS_IGNORE, p_err)
                  CALL mpi_recv(signature_imt_total,   1, MPI_INTEGER, signature_winner_rank, 106, &
                                p_comm_worker, MPI_STATUS_IGNORE, p_err)
                  CALL mpi_recv(signature_ptype_total,  1, MPI_INTEGER, signature_winner_rank, 107, &
                                p_comm_worker, MPI_STATUS_IGNORE, p_err)
                  CALL mpi_recv(signature_tol_total,    1, MPI_REAL8, signature_winner_rank, 108, &
                                p_comm_worker, MPI_STATUS_IGNORE, p_err)
                  CALL mpi_recv(signature_terms_total,  3, MPI_REAL8, signature_winner_rank, 109, &
                                p_comm_worker, MPI_STATUS_IGNORE, p_err)
               ENDIF
            ENDIF
         ENDIF

         IF (nbad_total > 0) THEN
            in_pair = (/ worst_abs, real(p_iam_worker, r8) /)
            CALL mpi_allreduce(in_pair, out_pair, 1, MPI_2DOUBLE_PRECISION, &
                               MPI_MAXLOC, p_comm_worker, p_err)
            winner_rank = nint(out_pair(2))
            IF (winner_rank /= 0) THEN
               IF (p_iam_worker == winner_rank) THEN
                  CALL mpi_send(balance_worst_ipatch, 1, MPI_INTEGER, 0, 91, p_comm_worker, p_err)
                  CALL mpi_send(balance_worst_imt,   1, MPI_INTEGER, 0, 92, p_comm_worker, p_err)
                  CALL mpi_send(balance_worst_ptype,  1, MPI_INTEGER, 0, 93, p_comm_worker, p_err)
                  CALL mpi_send(balance_worst_diag,  19, MPI_REAL8,   0, 94, p_comm_worker, p_err)
                  CALL mpi_send(balance_worst_sds, n_storage_diag, MPI_REAL8, 0, 95, p_comm_worker, p_err)
                  CALL mpi_send(balance_worst_sbeg, n_storage_diag, MPI_REAL8, 0, 96, p_comm_worker, p_err)
                  CALL mpi_send(balance_worst_send, n_storage_diag, MPI_REAL8, 0, 97, p_comm_worker, p_err)
                  CALL mpi_send(balance_worst_fcomp, n_flux_diag, MPI_REAL8, 0, 98, p_comm_worker, p_err)
               ELSEIF (p_iam_worker == 0) THEN
                  CALL mpi_recv(balance_worst_ipatch, 1, MPI_INTEGER, winner_rank, 91, &
                                p_comm_worker, MPI_STATUS_IGNORE, p_err)
                  CALL mpi_recv(balance_worst_imt,   1, MPI_INTEGER, winner_rank, 92, &
                                p_comm_worker, MPI_STATUS_IGNORE, p_err)
                  CALL mpi_recv(balance_worst_ptype,  1, MPI_INTEGER, winner_rank, 93, &
                                p_comm_worker, MPI_STATUS_IGNORE, p_err)
                  CALL mpi_recv(balance_worst_diag,  19, MPI_REAL8,   winner_rank, 94, &
                                p_comm_worker, MPI_STATUS_IGNORE, p_err)
                  CALL mpi_recv(balance_worst_sds, n_storage_diag, MPI_REAL8, winner_rank, 95, &
                                p_comm_worker, MPI_STATUS_IGNORE, p_err)
                  CALL mpi_recv(balance_worst_sbeg, n_storage_diag, MPI_REAL8, winner_rank, 96, &
                                p_comm_worker, MPI_STATUS_IGNORE, p_err)
                  CALL mpi_recv(balance_worst_send, n_storage_diag, MPI_REAL8, winner_rank, 97, &
                                p_comm_worker, MPI_STATUS_IGNORE, p_err)
                  CALL mpi_recv(balance_worst_fcomp, n_flux_diag, MPI_REAL8, winner_rank, 98, &
                                p_comm_worker, MPI_STATUS_IGNORE, p_err)
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      print_me = (p_is_worker .and. p_iam_worker == 0)
#else
      reduced_abs = worst_abs
      print_me    = .true.
#endif

      IF (print_me .and. nbad_total > 0) THEN
         ! "nbad_entries" = (patch, material) pairs exceeding tol.
         ! Upper bound on unique patches is nbad_entries / nmts
         ! (each patch contributes at most nmts entries).
         ! ptype: 0=soil 1=urban 2=wetland 3=glacier 4=waterbody, -1=unknown.
         WRITE(*,'(A,I8,A,E12.5,A,I8,A,I4,A,I3,A,I4,A,E12.5,A,E12.5,A,E12.5,&
                  &A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5,&
                  &A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5,&
                  &A,E12.5,A,E12.5,A,E12.5)') &
            'MT_BAL step report: nbad_entries=', nbad_total, &
            ' worst_abs_err=', reduced_abs, &
            ' @ipatch=', balance_worst_ipatch, &
            ' imt=', balance_worst_imt, &
            ' ptype=', balance_worst_ptype, &
            ' owner=', winner_rank, &
            ' err=', balance_worst_diag(1), &
            ' dS=', balance_worst_diag(2), &
            ' in=', balance_worst_diag(3), &
            ' out=', balance_worst_diag(4), &
            ' water_err=', balance_worst_diag(5), &
            ' water_err_R=', balance_worst_diag(6), &
            ' err_minus_water=', balance_worst_diag(7), &
            ' water_dS=', balance_worst_diag(8), &
            ' water_input=', balance_worst_diag(9), &
            ' water_output=', balance_worst_diag(10), &
            ' dS_minus_water_R=', balance_worst_diag(11), &
            ' in_minus_water_R=', balance_worst_diag(12), &
            ' out_minus_water_R=', balance_worst_diag(13), &
            ' step_evap=', balance_worst_diag(14), &
            ' step_rnof=', balance_worst_diag(15), &
            ' water_evap=', balance_worst_diag(16), &
            ' water_rnof=', balance_worst_diag(17), &
            ' evap_minus_water_R=', balance_worst_diag(18), &
            ' rnof_minus_water_R=', balance_worst_diag(19)
         WRITE(*,'(A,I8,A,I4,A,I3,A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5,&
                  &A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5)') &
            'MT_BAL_DCOMP @ipatch=', balance_worst_ipatch, &
            ' imt=', balance_worst_imt, &
            ' ptype=', balance_worst_ptype, &
            ' d_ldew=', balance_worst_sds(1), &
            ' d_soil_liq=', balance_worst_sds(2), &
            ' d_soil_ice=', balance_worst_sds(3), &
            ' d_wa=', balance_worst_sds(4), &
            ' d_wdsrf=', balance_worst_sds(5), &
            ' d_wetwat=', balance_worst_sds(6), &
            ' d_scv=', balance_worst_sds(7), &
            ' d_waterstorage=', balance_worst_sds(8), &
            ' d_leaf_iso=', balance_worst_sds(9), &
            ' d_surface_residue=', balance_worst_sds(10), &
            ' d_subsurface_residue=', balance_worst_sds(11), &
            ' d_solid=', balance_worst_sds(12)
         WRITE(*,'(A,I8,A,I4,A,I3,A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5,&
                  &A,E12.5,A,E12.5)') &
            'MT_BAL_FCOMP @ipatch=', balance_worst_ipatch, &
            ' imt=', balance_worst_imt, &
            ' ptype=', balance_worst_ptype, &
            ' precip=', balance_worst_fcomp(1), &
            ' evap=', balance_worst_fcomp(2), &
            ' rsur=', balance_worst_fcomp(3), &
            ' rsub=', balance_worst_fcomp(4), &
            ' rnof=', balance_worst_fcomp(5), &
            ' qinfl=', balance_worst_fcomp(6), &
            ' qcharge=', balance_worst_fcomp(7)
         WRITE(*,'(A,I8,A,I4,A,I3,A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5,&
                  &A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5)') &
            'MT_BAL_SEND @ipatch=', balance_worst_ipatch, &
            ' imt=', balance_worst_imt, &
            ' ptype=', balance_worst_ptype, &
            ' ldew=', balance_worst_send(1), &
            ' soil_liq=', balance_worst_send(2), &
            ' soil_ice=', balance_worst_send(3), &
            ' wa=', balance_worst_send(4), &
            ' wdsrf=', balance_worst_send(5), &
            ' wetwat=', balance_worst_send(6), &
            ' scv=', balance_worst_send(7), &
            ' waterstorage=', balance_worst_send(8), &
            ' leaf_iso=', balance_worst_send(9), &
            ' surface_residue=', balance_worst_send(10), &
            ' subsurface_residue=', balance_worst_send(11), &
            ' solid=', balance_worst_send(12)
      ENDIF

      IF (print_me .and. resid_nbad_total > 0) THEN
         ! Fires independently of nbad_total to identify entries where the
         ! numerical residual was explicitly excluded from hard balance accounting.
         WRITE(*,'(A,I8,A,E12.5)') &
            'MT_BAL residual note (excluded from tol check): n_entries=', &
            resid_nbad_total, ' worst_abs=', resid_abs_total
      ENDIF

      IF (print_me .and. resid_hard_nbad_total > DEF_MT_RESID_ABORT_NBAD) THEN
         WRITE(*,'(A,I8,A,E12.5,A,E12.5,A,I8,A,I4,A,I3,A,I4)') &
            'MT_BAL residual hard failure: n_entries=', resid_hard_nbad_total, &
            ' worst_abs=', resid_hard_abs_total, &
            ' tol=', resid_hard_tol_total, &
            ' @ipatch=', resid_hard_ipatch_total, &
            ' imt=', resid_hard_imt_total, &
            ' ptype=', resid_hard_ptype_total, &
            ' owner=', resid_winner_rank
         CALL CoLM_stop ('MT_BAL hard failure: numerical residual exceeds independent threshold')
      ENDIF

      IF (print_me .and. signature_nbad_total > 0) THEN
         WRITE(*,'(A,I8,A,E12.5,A,E12.5,A,I8,A,I4,A,I3,A,I4,A,E12.5,A,E12.5,A,E12.5)') &
            'MT_SIG step report: nbad_entries=', signature_nbad_total, &
            ' worst_abs=', signature_abs_total, &
            ' tol=', signature_tol_total, &
            ' @ipatch=', signature_ipatch_total, &
            ' imt=', signature_imt_total, &
            ' ptype=', signature_ptype_total, &
            ' owner=', signature_winner_rank, &
            ' input_minus_water_R=', signature_terms_total(1), &
            ' evap_minus_water_R=', signature_terms_total(2), &
            ' rnof_minus_water_R=', signature_terms_total(3)
      ENDIF

      IF (print_me .and. signature_nbad_total > DEF_MT_BALANCE_ABORT_NBAD) THEN
         CALL CoLM_stop ('MT_SIG hard failure: fixed-signature flux differs from water flux times R_init')
      ENDIF

      IF (print_me .and. nbad_total > DEF_MT_BALANCE_ABORT_NBAD) THEN
         CALL CoLM_stop ('MT_BAL hard failure: material balance error exceeds aggregate threshold')
      ENDIF

      balance_worst_err    = 0._r8
      balance_worst_diag   = 0._r8
      balance_worst_sbeg   = 0._r8
      balance_worst_send   = 0._r8
      balance_worst_sds    = 0._r8
      balance_worst_fcomp  = 0._r8
      balance_worst_ipatch = 0
      balance_worst_imt   = 0
      balance_worst_ptype  = -1
      balance_nbad         = 0
      signature_worst_abs      = 0._r8
      signature_worst_tol      = 0._r8
      signature_worst_terms    = 0._r8
      signature_worst_ipatch   = 0
      signature_worst_imt     = 0
      signature_worst_ptype    = -1
      signature_nbad           = 0
      resid_worst_abs      = 0._r8
      resid_nbad           = 0
      resid_hard_worst_abs = 0._r8
      resid_hard_worst_tol = 0._r8
      resid_hard_worst_ipatch = 0
      resid_hard_worst_imt   = 0
      resid_hard_worst_ptype  = -1
      resid_hard_nbad      = 0
   END SUBROUTINE material_balance_report

END MODULE MOD_MT_Conservation
#endif
