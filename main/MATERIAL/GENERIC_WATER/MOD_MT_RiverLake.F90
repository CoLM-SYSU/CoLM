#include <define.h>

#ifdef GridRiverLakeFlow
#ifdef MaterialTransport
MODULE MOD_MT_RiverLake
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   Passive material transport module for GridRiverLakeFlow.
!   Supports arbitrary number of materials (e.g., delta18O, deltaD).
!
!   Physics:
!     - Upwind advection: material flux = concentration × water flux
!     - Flux limiter for mass conservation (prevents negative storage)
!     - Fully-mixed assumption within each unit catchment
!
!   Reference: CaMa-Flood material module (cmf_ctrl_mt_mod.F90)
!
! Created by CoLM, Apr 2026
!-------------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   ! Share the single authoritative `nmts` from MOD_MT_Defs to
   ! avoid the land/river material modules holding two independent copies
   ! that could silently diverge on re-init.
   USE MOD_MT_Defs, only: nmts, mt_uses_generic_water, material_concentration_units, &
      material_build_descriptor_identity, MATERIAL_DESCRIPTOR_IDENTITY_WIDTH
   USE MOD_MT_Solute, only: mt_has_dissolved_limit, mt_equilibrate_dissolved
   IMPLICIT NONE

   !-------------------------------------------------------------------------------------
   ! Module variables
   !-------------------------------------------------------------------------------------
   character(len=32), allocatable :: material_names(:)  ! Material names

   ! Schema 1 is the first river/lake material restart contract that carries
   ! the complete physical descriptor, rather than only count + name hash.
   integer, parameter :: RIVER_MT_RESTART_SCHEMA_VERSION = 2
   integer, save :: river_restart_schema_loaded = 0
   real(r8), parameter :: MT_RESTART_NEGATIVE_DUST = 1.0e-12_r8
   ! Relative roundoff admitted after a donor-limited update.  A cell drained
   ! exactly (rate = (mass + inflow)/outflow) lands within a few ulp of the
   ! gross terms of that update, which for a river cell (mass ~ volume*conc)
   ! is far above the absolute dust above: a spatial run stopped at step 2 on
   ! a cell of 3.5e4 that ended at -3.5e-12 (one ulp).
   real(r8), parameter :: MT_UPDATE_ROUNDOFF = 1.0e-12_r8
   ! The coupled donor limiter iterates dimensionless rates in [0,1] through a
   ! monotone non-decreasing map (rate_new = max(rate_old, feasible)): every
   ! iterate is feasible, so a cell can never export more material than it holds.
   ! The loop ends when the largest per-iteration rate increment is at or below
   ! MT_LIMITER_RATE_TOL, and is otherwise capped at the 2*N+1 topology bound
   ! (exhausting the cap is fatal).  The tolerance is ABSOLUTE: rates are
   ! clamped to [0,1], so a relative scale would always be 1.  Loosening it
   ! (e.g. to 1e-9) lets a slowly propagating chain stop early and leave material
   ! for the dry-cell sink to remove instead of delivering it downstream.
   real(r8), parameter :: MT_LIMITER_RATE_TOL = 1.0e-12_r8
   ! Fast-path DIAGNOSTIC threshold, not a truncation: the loop still runs to
   ! the 2*N+1 topology bound, so the fixed point is always reached.  A long
   ! river chain legitimately needs one round per edge, and truncating would
   ! under-export material that the dry-cell sink can later remove instead of
   ! delivering it downstream (mt_dry_drain is a budgeted sink, not a delay),
   ! so crossing this threshold is only counted and reported.
   integer,  parameter :: MT_LIMITER_SOFT_ITER = 32
   real(r8), parameter :: MT_LIMITER_OUT_TINY = 1.0e-30_r8

   ! Diagnostics: iterations the donor limiter actually needed (per call).
   integer, save :: mt_limiter_calls     = 0
   integer, save :: mt_limiter_iter_sum  = 0
   integer, save :: mt_limiter_iter_peak = 0
   integer, save :: mt_limiter_over_soft = 0

   ! State variables (prognostic)
   ! mt_mass is built from R_default * water_volume, so its units are
   ! [R*m3] (isotope ratio times cell water volume), NOT kg of heavy water.
   ! mt_conc = mt_mass / volwater is therefore a dimensionless ratio [R].
   real(r8), allocatable :: mt_mass  (:,:)   ! Material amount [R*m3] (nmts, numucat)
   real(r8), allocatable :: mt_conc  (:,:)   ! Material ratio  [R]    (nmts, numucat)

   ! Protected-side (behind-levee) material pool paired with
   ! `levsto` in MOD_Grid_RiverLakeLevee. Levee repartition moves water
   ! between visible-side (volwater_ucat, driving mt_mass) and
   ! protected-side storage; without its own material pool the mass would
   ! stay pinned to the visible side while water moves, systematically
   ! inflating visible-side concentration (and orphaning mass when
   ! visible-side goes dry).
   real(r8), allocatable :: mt_levsto(:,:)   ! Protected-side material amount [R*m3] (nmts, numucat)
   ! Only allocated when at least one generic material has finite solubility.
   real(r8), allocatable :: mt_solid(:,:), mt_levsto_solid(:,:)

   ! Flux variables (diagnostic, per routing period)
   real(r8), allocatable :: mt_flux_out  (:,:) ! Material outflux [R*m3/s] (nmts, numucat)
   real(r8), allocatable :: mt_dry_drain(:,:) ! Terminal (mouth) sink from forced drain [R*m3]
   ! Running totals of the dry-cell bookkeeping, per species:
   !   1 = drained mass handed to the split delivery (internal cells)
   !   2 = mass actually delivered to receivers
   !   3 = current pending inventory (computed when reporting, not cumulative)
   !   4 = mass booked as the terminal sink
   ! mt_dry_drain itself is cleared at the end of every routing call, so the
   ! totals have to be accumulated here to be reportable at the end.
   real(r8), allocatable :: mt_dry_stat(:,:)  ! [R*m3]
   real(r8), allocatable :: mt_dry_fwd(:,:)   ! Internal drained mass forwarded downstream [R*m3]
   real(r8), allocatable :: mt_dry_rev(:,:)   ! (legacy) unique-reverse buffer
   real(r8), allocatable :: mt_dry_split(:,:) ! Internal drained mass split by actual outflow flux [R*m3]
   real(r8), allocatable :: mt_dry_pending(:,:) ! Drained mass awaiting a known outflow direction [R*m3]
   real(r8), allocatable :: mt_decay_source(:,:) ! Net reactive source/sink over routing period
   ! Nitrogen converted to N2O at the land-to-river entry, positive mass
   ! [kg N] per routing period.  Diagnostic only: the mass is removed from
   ! the river input before it enters acc_mt_inp, so the river balance stays
   ! closed without an extra source term.
   real(r8), allocatable :: mt_n2o_source(:,:)

   ! Input: material amount flux from runoff
   real(r8), allocatable :: acc_mt_inp  (:,:) ! Accumulated material input [R*m3] (nmts, numucat)
   real(r8), allocatable :: mt_inp_buf  (:,:) ! Buffered runoff material awaiting release [R*m3]
   real(r8), allocatable :: acc_rnof_ref (:)   ! Accumulated runoff water volume paired with acc_mt_inp [m3]

   ! Bifurcation net flux (saved from last mt_riverlake_substep for diagnostics)
   real(r8), allocatable :: mt_bif_net_saved (:,:) ! Per-cell net bif flux [mass/s] (nmts, numucat)

   ! Near-dry volume floor used only for transport concentration
   ! stabilization; prognostic material mass always remains single-pool.
   real(r8), parameter :: mt_v_dry_off = 1.e-6_r8
   ! Display/history diagnostics should not let bathtub-scale residual
   ! volumes dominate isotope delta min/max. This does not alter state.
   real(r8), parameter :: mt_delta_diag_vmin = 1._r8
   ! History accumulators
   real(r8), allocatable :: a_mt_storage_mass(:,:) ! Accumulated visible storage material [R*m3*s]
   real(r8), allocatable :: a_water_storage(:)       ! Accumulated visible water storage [m3*s]
   real(r8), allocatable :: a_mt_levsto_mass(:,:)   ! Accumulated protected storage material [R*m3*s]
   real(r8), allocatable :: a_levsto_water(:)        ! Accumulated protected water storage [m3*s]
   real(r8), allocatable :: a_mt_out    (:,:) ! Accumulated material outflux [mass/s * s] (nmts, numucat)
   real(r8), allocatable :: a_mt_bifout (:,:) ! Accumulated material bif net flux [mass/s * s] (nmts, numucat)
   real(r8), allocatable :: a_mt_acctime(:)   ! Material history elapsed time [s]

   ! N05 terminal-window snapshot.  When a run stops inside a history window,
   ! hist_out writes the partial window and flush_acc then zeroes the live
   ! accumulators before the physical restart is written.  The snapshot keeps
   ! the pre-flush rows so the restart continues the same window (matching the
   ! host sidecar) instead of losing the first part of the window.
   real(r8), allocatable, save :: histsnap_stor(:,:), histsnap_levsto(:,:)
   real(r8), allocatable, save :: histsnap_out(:,:), histsnap_bifout(:,:)
   ! N2O converted at the river entry (kg N) must travel with the other
   ! history accumulators across a restart; otherwise a split run divides a
   ! window whose N2O numerator restarts from zero.
   real(r8), allocatable, save :: histsnap_n2o(:,:)
   real(r8), allocatable, save :: histsnap_wstor(:), histsnap_wlevsto(:), histsnap_acctime(:)
   logical, save :: histsnap_valid = .false.

   ! Routing material substep workspace. These buffers are reused across the
   ! many adaptive routing substeps to avoid repeated heap allocate/free churn.
   real(r8), allocatable, save :: conc_next(:)
   real(r8), allocatable, save :: flux_ups(:)
   real(r8), allocatable, save :: mt_flux(:)
   real(r8), allocatable, save :: mt_conc_flux(:)
   real(r8), allocatable, save :: mt_prot_conc_flux(:)
   real(r8), allocatable, save :: conc_dn_pth(:)
   real(r8), allocatable, save :: mt_prot_conc_dn_pth(:)
   real(r8), allocatable, save :: mt_pth_single(:)
   real(r8), allocatable, save :: mt_pth_lev(:)
   real(r8), allocatable, save :: mt_pth_levee(:,:)
   real(r8), allocatable, save :: bif_recv(:)
   real(r8), allocatable, save :: bif_lev_recv(:)
   real(r8), allocatable, save :: bif_net(:)
   real(r8), allocatable, save :: mt_bif_lev_net(:)
   real(r8), allocatable, save :: has_levee_r8(:)
   real(r8), allocatable, save :: has_levee_dn_pth(:)
   real(r8), allocatable, save :: dt_ucat(:)
   real(r8), allocatable, save :: dt_dn_pth(:)
   real(r8), allocatable, save :: mt_in_mass(:)
   real(r8), allocatable, save :: mt_in_mass_lev(:)
   real(r8), allocatable, save :: mt_out_mass(:)
   real(r8), allocatable, save :: mt_out_mass_lev(:)
   real(r8), allocatable, save :: rate_cell(:)
   real(r8), allocatable, save :: rate_cell_lev(:)
   real(r8), allocatable, save :: rate_next(:)
   real(r8), allocatable, save :: rate_dn_pth(:)
   real(r8), allocatable, save :: rate_dn_pth_lev(:)
   real(r8), allocatable, save :: mt_inp_step(:)
   real(r8), allocatable, save :: mt_dn_out_vis_pth(:)
   real(r8), allocatable, save :: mt_dn_out_lev_pth(:)
   real(r8), allocatable, save :: mt_dn_out_vis_recv(:)
   real(r8), allocatable, save :: mt_dn_out_lev_recv(:)

   PRIVATE :: conc_next, flux_ups, mt_flux, mt_conc_flux, mt_prot_conc_flux
   PRIVATE :: conc_dn_pth, mt_prot_conc_dn_pth, mt_pth_single, mt_pth_lev
   PRIVATE :: mt_pth_levee, bif_recv, bif_lev_recv, bif_net, mt_bif_lev_net
   PRIVATE :: has_levee_r8, has_levee_dn_pth, dt_ucat, dt_dn_pth
   PRIVATE :: mt_in_mass, mt_in_mass_lev, mt_out_mass, mt_out_mass_lev
   PRIVATE :: rate_cell, rate_cell_lev, rate_next, rate_dn_pth, rate_dn_pth_lev
   PRIVATE :: mt_inp_step, mt_dn_out_vis_pth, mt_dn_out_lev_pth
   PRIVATE :: mt_dn_out_vis_recv, mt_dn_out_lev_recv
   PRIVATE :: ensure_mt_substep_workspace, release_mt_substep_workspace
   PRIVATE :: ensure_real1_workspace, ensure_real2_workspace
   PRIVATE :: material_bif_path_levee_sides
   PRIVATE :: riverlake_restart_descriptor_compatible
   PRIVATE :: probe_riverlake_restart_vector
   PRIVATE :: validate_riverlake_restart_state

   PUBLIC :: mt_riverlake_init
   PUBLIC :: mt_riverlake_init_from_water
   PUBLIC :: mt_riverlake_input_from_runoff
   PUBLIC :: mt_riverlake_flush_acc
   PUBLIC :: mt_riverlake_read_restart
   PUBLIC :: mt_riverlake_drop_hist_window
   PUBLIC :: mt_riverlake_write_restart
   PUBLIC :: mt_riverlake_write_hist_acc_restart, mt_riverlake_read_hist_acc_restart
   PUBLIC :: mt_riverlake_snapshot_hist_acc
   PUBLIC :: mt_riverlake_write_history
   PUBLIC :: mt_riverlake_final
   PUBLIC :: check_mt_state
   PUBLIC :: mt_riverlake_substep
   PUBLIC :: mt_riverlake_limiter_stats
   PUBLIC :: mt_riverlake_diag_accumulate
   ! Exposed so the inland-depression overflow fix in
   ! MOD_Grid_RiverLakeFlow can use the same reservoir/levee/floodplain
   ! volume selection as the material transport instead of hard-wiring
   ! topo_rivstomax (which disagrees on leveed cells).
   PUBLIC :: get_cell_volume
   PUBLIC :: mt_inp_buf
   PUBLIC :: acc_rnof_ref
   PUBLIC :: mt_levsto
   PUBLIC :: mt_solid, mt_levsto_solid, mt_riverlake_equilibrate_cell
   PUBLIC :: mt_dry_drain
   PRIVATE :: mt_bif_dry_prepare, mt_bif_dry_deliver
   PUBLIC :: mt_dry_stat
   PUBLIC :: mt_riverlake_release_pending_drain
   PUBLIC :: mt_decay_source
   PUBLIC :: mt_riverlake_levee_repartition

CONTAINS

   !-------------------------------------------------------------------------------------
   ! Initialize material module
   !-------------------------------------------------------------------------------------
   SUBROUTINE mt_riverlake_init ()

   USE MOD_Grid_RiverLakeNetwork, only: numucat
   USE MOD_MT_Defs, only: material_defs_init, materials
   IMPLICIT NONE

      integer :: i
      logical :: has_finite_solute

      CALL mt_riverlake_final()
      CALL material_defs_init()
      IF (nmts <= 0) RETURN

      ! MOD_MT_Defs owns parsing, sanitisation, defaults, and de-duplication.
      ! River/lake keeps a local copy only because restart/history names use it.
      allocate (material_names(nmts))
      DO i = 1, nmts
         material_names(i) = materials(i)%name
      ENDDO
      has_finite_solute = .false.
      DO i = 1, nmts
         has_finite_solute = has_finite_solute .or. mt_has_dissolved_limit(i)
      ENDDO

      IF (p_is_worker) THEN
         ! Allocate on ALL workers (zero-length if numucat=0) for MPI safety.
         allocate (mt_mass     (nmts, numucat))
         allocate (mt_conc     (nmts, numucat))
         allocate (mt_flux_out (nmts, numucat))
         allocate (mt_dry_drain(nmts, numucat))
         allocate (mt_dry_stat(nmts, 4))
         allocate (mt_dry_fwd  (nmts, numucat))
         allocate (mt_dry_rev  (nmts, numucat))
         allocate (mt_dry_split(nmts, numucat))
         allocate (mt_dry_pending(nmts, numucat))
         allocate (mt_decay_source(nmts, numucat))
         allocate (mt_n2o_source  (nmts, numucat))
         allocate (acc_mt_inp        (nmts, numucat))
         allocate (mt_inp_buf        (nmts, numucat))
         allocate (acc_rnof_ref       (numucat))
         allocate (mt_bif_net_saved  (nmts, numucat))
         allocate (a_mt_storage_mass (nmts, numucat))
         allocate (a_water_storage    (numucat))
         allocate (a_mt_levsto_mass (nmts, numucat))
         allocate (a_levsto_water    (numucat))
         allocate (a_mt_out          (nmts, numucat))
         allocate (a_mt_bifout       (nmts, numucat))
         allocate (a_mt_acctime      (numucat))
         allocate (mt_levsto         (nmts, numucat))
         IF (has_finite_solute) THEN
            allocate (mt_solid(nmts, numucat), mt_levsto_solid(nmts, numucat))
            mt_solid = 0._r8
            mt_levsto_solid = 0._r8
         ENDIF

         mt_mass     = 0._r8
         mt_conc     = 0._r8
         mt_flux_out  = 0._r8
         mt_dry_drain = 0._r8
         mt_dry_stat  = 0._r8
         mt_dry_fwd   = 0._r8
         mt_dry_rev   = 0._r8
         mt_dry_split = 0._r8
         mt_dry_pending = 0._r8
         mt_decay_source = 0._r8
         mt_n2o_source   = 0._r8
         acc_mt_inp  = 0._r8
         mt_inp_buf  = 0._r8
         acc_rnof_ref = 0._r8
         mt_bif_net_saved = 0._r8
         a_mt_storage_mass = 0._r8
         a_water_storage = 0._r8
         a_mt_levsto_mass = 0._r8
         a_levsto_water = 0._r8
         a_mt_out    = 0._r8
         a_mt_bifout = 0._r8
         a_mt_acctime = 0._r8
         mt_levsto   = 0._r8
      ENDIF

      IF (p_is_master) THEN
         write(*,'(A,I4,A)') ' Material module initialised with ', nmts, ' materials:'
         write(*,'(A,*(A,:,", "))') '   Names: ', (trim(material_names(i)), i=1, nmts)
      ENDIF

   END SUBROUTINE mt_riverlake_init


   !-------------------------------------------------------------------------------------
   ! Cold-start river material mass from the current water volume.
   ! Mirrors the land-side mt_riverlake_init_from_water logic: when no material
   ! restart variable is present in the gridriver restart file but the
   ! hydrology restart provides a non-zero wdsrf_ucat/volresv, seed
   ! mt_mass = volwater * R_init. Without this, any hot-start of water
   ! starts with material=0 and bleeds spurious dilution into downstream
   ! concentrations until runoff "re-colours" the reservoir.
   !-------------------------------------------------------------------------------------
   SUBROUTINE mt_riverlake_init_from_water (wdsrf, volresv_in, ucat2resv_in, missing_mask, is_built_resv)

   USE MOD_Grid_RiverLakeNetwork, only: numucat, lake_type
   USE MOD_Grid_RiverLakeLevee,   only: has_levee, levsto
   USE MOD_MT_Defs,           only: material_init_water_ratio
   IMPLICIT NONE

   real(r8), intent(in) :: wdsrf(:)
   real(r8), intent(in) :: volresv_in(:)
   integer,  intent(in) :: ucat2resv_in(:)
   ! Optional per-material selector. When present, only materials with
   ! missing_mask(imt) == .true. are re-seeded from the water volume;
   ! this avoids overwriting materials that loaded successfully from a
   ! partially-complete restart file (e.g. user renamed one material or
   ! added a new material to DEF_MT_NAMES without re-writing the file).
   ! Without the mask every material is cold-started (equivalent to
   ! passing an all-.true. mask).
   logical,  optional, intent(in) :: missing_mask(:)
   logical,  optional, intent(in) :: is_built_resv(:)

   integer  :: i, imt
   real(r8) :: volwater, R_init
   logical  :: do_init
   integer  :: n_init

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN
      IF (.not. allocated(mt_mass)) RETURN

      n_init = 0
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         do_init = .true.
         IF (present(missing_mask)) THEN
            IF (imt <= size(missing_mask)) do_init = missing_mask(imt)
         ENDIF
         IF (.not. do_init) CYCLE
         R_init = material_init_water_ratio(imt)
         DO i = 1, numucat
            IF (present(is_built_resv)) THEN
               CALL get_cell_volume(i, wdsrf(i), volresv_in, ucat2resv_in, volwater, is_built_resv(i))
            ELSE
               CALL get_cell_volume(i, wdsrf(i), volresv_in, ucat2resv_in, volwater)
            ENDIF
            mt_mass(imt, i) = max(volwater, 0._r8) * R_init
            ! Seed protected-side pool from levsto so a cold
            ! start with water already behind the levee doesn't produce
            ! a phantom conservation jump on the first repartition.
            IF (DEF_USE_LEVEE .and. allocated(mt_levsto) .and. allocated(has_levee)) THEN
               IF (has_levee(i) .and. lake_type(i) /= 2) THEN
                  mt_levsto(imt, i) = max(levsto(i), 0._r8) * R_init
               ELSE
                  mt_levsto(imt, i) = 0._r8
               ENDIF
            ENDIF
            IF (allocated(mt_solid) .and. mt_has_dissolved_limit(imt)) THEN
               mt_solid(imt, i) = 0._r8
               mt_levsto_solid(imt, i) = 0._r8
               CALL mt_equilibrate_dissolved(imt, volwater, mt_mass(imt, i), mt_solid(imt, i))
               IF (allocated(levsto)) CALL mt_equilibrate_dissolved(imt, max(levsto(i), 0._r8), &
                  mt_levsto(imt, i), mt_levsto_solid(imt, i))
            ENDIF
         ENDDO
         ! Recompute concentration so diagnostics consistent before first step.
         DO i = 1, numucat
            IF (present(is_built_resv)) THEN
               CALL get_cell_volume(i, wdsrf(i), volresv_in, ucat2resv_in, volwater, is_built_resv(i))
            ELSE
               CALL get_cell_volume(i, wdsrf(i), volresv_in, ucat2resv_in, volwater)
            ENDIF
            CALL update_mt_concentration(imt, i, volwater)
         ENDDO
         n_init = n_init + 1
      ENDDO

      IF (p_is_master .and. n_init > 0) THEN
         write(*,'(A,I0,A)') '  River material cold-started from water volume for ', n_init, ' material(s).'
      ENDIF

   END SUBROUTINE mt_riverlake_init_from_water


   !-------------------------------------------------------------------------------------
   ! Accumulate heavy-water mass input from runoff.
   ! Called each land-model timestep (before routing accumulation threshold).
   !
   ! rnof_uc_depth(i):           catchment runoff volume [m3] after area-weighted mapping.
   ! mt_rnof_ext(imt, i):      (optional) material amount from land material
   !                             in matching [R*m3] units.
   !
   ! If mt_rnof_ext is present, use it directly (coupled to land material system).
   ! Otherwise, compute default heavy-water mass from ref_ratio and init_delta.
   !-------------------------------------------------------------------------------------
   SUBROUTINE mt_riverlake_input_from_runoff (rnof_uc_depth, numucat_in, mt_rnof_ext)

      ! Re-uses the per-material ref_ratio / init_delta cached by
      ! material_defs_init (MOD_MT_Defs.F90:82-84). Previously this
      ! routine re-parsed the CSV namelist strings via a local helper on
      ! every land timestep — wasteful and a divergence hazard (two
      ! parse code paths for the same inputs).
      USE MOD_MT_Defs, only: material_init_water_ratio
      USE MOD_MT_NitrogenN2O, only: mt_n2o_ef5_for
      IMPLICIT NONE
      real(r8), intent(in) :: rnof_uc_depth(:)
      integer,  intent(in) :: numucat_in
      real(r8), intent(in), optional :: mt_rnof_ext(:,:)

      integer :: i, imt
      real(r8), allocatable :: R_default(:)
      real(r8) :: mass_in, ef5, n2o_in

      ! MATERIAL may be compiled while the runtime registry intentionally
      ! contains zero materials (DEF_MT_NUM=0).  In that mode mt_riverlake_init()
      ! returns without allocating the routing buffers, so the input path must
      ! be a no-op instead of touching acc_rnof_ref / acc_mt_inp.
      IF (nmts <= 0) RETURN
      IF (.not. allocated(acc_rnof_ref)) RETURN
      IF (.not. allocated(acc_mt_inp)) RETURN
      IF (numucat_in <= 0) RETURN

      ! Cache per-material R_default once; absent mt_rnof_ext it is used
      ! for every cell and was previously recomputed inside the inner
      ! loop (delta_to_R is pure but costs trig/div per call).
      IF (.not. present(mt_rnof_ext)) THEN
         allocate(R_default(nmts))
         DO imt = 1, nmts
            IF (.not. mt_uses_generic_water(imt)) CYCLE
            R_default(imt) = material_init_water_ratio(imt)
         ENDDO
      ENDIF

      ! Mirror the water-side `acc_rnof_uc` accumulation
      ! at MOD_Grid_RiverLakeFlow.F90:329 exactly. That line has no
      ! sign filter — both positive and negative rnof contributions are
      ! summed. The previous `IF (rnof_uc_depth > 0)` guard here let
      ! FP-noise negatives drift the water/material accumulators apart,
      ! which would make restart recovery and levee pending-pool
      ! repartition use inconsistent water and material windows.
      DO i = 1, numucat_in
         acc_rnof_ref(i) = acc_rnof_ref(i) + rnof_uc_depth(i)
         DO imt = 1, nmts
            IF (.not. mt_uses_generic_water(imt)) CYCLE
            IF (present(mt_rnof_ext)) THEN
               mass_in = mt_rnof_ext(imt, i)
            ELSE
               mass_in = R_default(imt) * rnof_uc_depth(i)
            ENDIF
            ! Indirect N2O (design A): at the land-to-river entry, a
            ! configured fraction of the dissolved nitrogen is converted to
            ! N2O and leaves to the atmosphere.  Water is unchanged and only
            ! the dissolved mass is reduced, so the river mass balance
            ! closes without a separate source term.  The per-species policy
            ! (which solutes, and the factor: parameter file or
            ! &nl_colm_material override) lives in FAMILY_SOLUTE's
            ! MOD_MT_NitrogenN2O; this generic river code only applies the
            ! returned fraction.  0 disables the conversion for a species,
            ! which keeps those configurations bit-identical.
            ! This is the inventory arithmetic (N_leached x EF5), not a
            ! mechanistic water-body model: it does not represent where the
            ! emission occurs or the water residence time, and it must not be
            ! applied again downstream.  Reference values and sources are in
            ! MOD_MT_NitrogenN2O and the standard solute parameter files.
            ef5 = mt_n2o_ef5_for(imt)
            IF (ef5 > 0._r8 .and. mass_in > 0._r8) THEN
               n2o_in = mass_in * ef5
               mass_in = mass_in - n2o_in
               IF (allocated(mt_n2o_source)) &
                  mt_n2o_source(imt, i) = mt_n2o_source(imt, i) + n2o_in
            ENDIF
            acc_mt_inp(imt, i) = acc_mt_inp(imt, i) + mass_in
         ENDDO
      ENDDO

      IF (allocated(R_default)) deallocate(R_default)

   END SUBROUTINE mt_riverlake_input_from_runoff


   !-------------------------------------------------------------------------------------
   ! Helper: get water volume for a unit catchment cell, respecting reservoir
   ! state and the prognostic routing volume. For levee cells volwater_ucat
   ! is the visible-side volume; for ordinary cells it is the total volume.
   !-------------------------------------------------------------------------------------
   SUBROUTINE get_cell_volume (icell, wdsrf_cell, volresv_in, ucat2resv_in, volwater, is_built_resv_cell)

   USE MOD_Grid_RiverLakeNetwork, only: floodplain_curve, lake_type
   USE MOD_Grid_RiverLakeLevee,   only: has_levee, levsto, levee_visible_volume_from_stage
   USE MOD_Grid_RiverLakeTimeVars, only: volwater_ucat, volwater_ucat_valid
   USE MOD_Vars_Global,           only: spval
   IMPLICIT NONE

   integer,  intent(in)  :: icell
   real(r8), intent(in)  :: wdsrf_cell
   real(r8), intent(in)  :: volresv_in(:)
   integer,  intent(in)  :: ucat2resv_in(:)
   real(r8), intent(out) :: volwater
   logical, intent(in), optional :: is_built_resv_cell
   logical               :: has_levee_cell
   logical               :: use_reservoir
   real(r8), parameter   :: stage_restart_tol = 1.e-5_r8

      has_levee_cell = .false.
      IF (DEF_USE_LEVEE .and. allocated(has_levee)) THEN
         IF (icell >= 1 .and. icell <= size(has_levee)) has_levee_cell = has_levee(icell)
      ENDIF

      use_reservoir = lake_type(icell) == 2
      IF (present(is_built_resv_cell)) use_reservoir = is_built_resv_cell
      IF (use_reservoir .and. size(volresv_in) > 0 .and. size(ucat2resv_in) > 0) THEN
         IF (ucat2resv_in(icell) > 0 .and. ucat2resv_in(icell) <= size(volresv_in) &
            .and. volresv_in(ucat2resv_in(icell)) /= spval) THEN
            volwater = volresv_in(ucat2resv_in(icell))
         ELSE
            volwater = floodplain_curve(icell)%volume(wdsrf_cell)
         ENDIF
      ELSEIF (has_levee_cell) THEN
         IF (volwater_ucat_valid .and. allocated(volwater_ucat)) THEN
            IF (icell <= size(volwater_ucat)) THEN
               ! Backward compatibility: old restarts may carry
               ! volwater_ucat as a zero placeholder even for wet cells.
               ! Reconstruct those from stage instead of treating levee
               ! visible storage as dry.
               IF (volwater_ucat(icell) > 0._r8 .or. wdsrf_cell <= stage_restart_tol) THEN
                  volwater = volwater_ucat(icell)
                  volwater = max(volwater, 0._r8)
                  RETURN
               ENDIF
            ENDIF
         ENDIF

         IF (allocated(levsto)) THEN
            volwater = levee_visible_volume_from_stage(icell, wdsrf_cell, levsto(icell))
         ELSE
            volwater = floodplain_curve(icell)%volume(wdsrf_cell)
         ENDIF
         volwater = max(volwater, 0._r8)
      ELSE
         IF (volwater_ucat_valid .and. allocated(volwater_ucat)) THEN
            IF (icell <= size(volwater_ucat)) THEN
               IF (volwater_ucat(icell) > 0._r8 .or. wdsrf_cell <= stage_restart_tol) THEN
                  volwater = volwater_ucat(icell)
                  volwater = max(volwater, 0._r8)
                  RETURN
               ENDIF
            ENDIF
         ENDIF

         volwater = floodplain_curve(icell)%volume(wdsrf_cell)
         volwater = max(volwater, 0._r8)
      ENDIF

   END SUBROUTINE get_cell_volume

   ! Keep only dissolved material mobile; the solid inventories stay in their
   ! respective compartments through drydown and can dissolve on rewetting.
   SUBROUTINE mt_riverlake_equilibrate_cell (icell, visible_water, protected_water)
      USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
      integer, intent(in) :: icell
      real(r8), intent(in) :: visible_water, protected_water
      integer :: imt

      IF (.not. allocated(mt_solid)) RETURN
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         IF (.not. mt_has_dissolved_limit(imt)) CYCLE
         IF (.not. ieee_is_finite(mt_mass(imt, icell)) .or. &
             .not. ieee_is_finite(mt_levsto(imt, icell)) .or. &
             .not. ieee_is_finite(mt_solid(imt, icell)) .or. &
             .not. ieee_is_finite(mt_levsto_solid(imt, icell))) &
            CALL CoLM_stop('non-finite river solute pool before dissolution equilibrium')
         IF (min(mt_mass(imt, icell), mt_levsto(imt, icell), &
                 mt_solid(imt, icell), mt_levsto_solid(imt, icell)) < -MT_RESTART_NEGATIVE_DUST) &
            CALL CoLM_stop('negative river solute pool before dissolution equilibrium')
         mt_mass(imt, icell) = max(0._r8, mt_mass(imt, icell))
         mt_levsto(imt, icell) = max(0._r8, mt_levsto(imt, icell))
         mt_solid(imt, icell) = max(0._r8, mt_solid(imt, icell))
         mt_levsto_solid(imt, icell) = max(0._r8, mt_levsto_solid(imt, icell))
         CALL mt_equilibrate_dissolved(imt, max(visible_water, 0._r8), &
            mt_mass(imt, icell), mt_solid(imt, icell))
         CALL mt_equilibrate_dissolved(imt, max(protected_water, 0._r8), &
            mt_levsto(imt, icell), mt_levsto_solid(imt, icell))
      ENDDO
   END SUBROUTINE mt_riverlake_equilibrate_cell


   !-------------------------------------------------------------------------------------
   ! Repartition material mass between visible-side (`mt_mass`)
   ! and protected-side (`mt_levsto`) pools in lockstep with the water
   ! redistribution that `levee_fldstg` just performed.
   !
   ! Inputs are the visible-side and protected-side water volumes BEFORE
   ! and AFTER the levee flood-staging call. Two conventions:
   !   d_lev > 0 : water moved visible -> protected. Material transferred
   !               at the visible-side's current ratio.
   !   d_lev < 0 : water moved protected -> visible. Material transferred
   !               at the protected-side's current ratio.
   ! This preserves asymmetric concentrations across the levee (unlike
   ! a uniform pool-ratio split which homogenises the two compartments).
   !-------------------------------------------------------------------------------------
   SUBROUTINE mt_riverlake_levee_repartition (icell, vis_vol_bef, levsto_bef, vis_vol_aft, levsto_aft, &
                                        pending_mt_pool)

   USE MOD_MT_Defs, only: mt_tiny, materials
   IMPLICIT NONE
   integer,  intent(in) :: icell
   real(r8), intent(in) :: vis_vol_bef, levsto_bef
   real(r8), intent(in) :: vis_vol_aft, levsto_aft
   ! Optional pool of period-accumulated runoff material that is conceptually on
   ! the visible side (its water has already been folded into vis_vol_bef) but
   ! has not yet been merged into mt_mass. Caller passes mt_inp_buf(:, icell)
   ! at the pre-material-substep call site so the visible-side ratio reflects
   ! the post-input state. Without this, the ratio is diluted
   ! (denominator includes runoff water but numerator does not include the
   ! matching runoff material), leaking the Phase-1 R_init invariant.
   real(r8), intent(inout), optional :: pending_mt_pool(:)

   integer  :: imt
   real(r8) :: d_lev, ratio, mt_move
   real(r8) :: vis_pending, vis_total
   real(r8) :: debit_mass, debit_pending

      IF (.not. allocated(mt_mass) .or. .not. allocated(mt_levsto)) RETURN
      IF (icell < 1 .or. icell > size(mt_mass, 2)) RETURN

      d_lev = levsto_aft - levsto_bef

      CALL mt_riverlake_equilibrate_cell(icell, vis_vol_bef, levsto_bef)
      IF (abs(d_lev) < mt_tiny) THEN
         CALL mt_riverlake_equilibrate_cell(icell, vis_vol_aft, levsto_aft)
         RETURN
      ENDIF

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         vis_pending = 0._r8
         IF (present(pending_mt_pool)) THEN
            IF (imt <= size(pending_mt_pool)) &
               vis_pending = max(pending_mt_pool(imt), 0._r8)
         ENDIF
         IF (d_lev > 0._r8) THEN
            ! Visible -> protected. Use post-input ratio (mt_mass + pending) /
            ! vis_vol_bef so the runoff material that already shares vis_vol_bef
            ! contributes to the concentration moved across the levee. Debit
            ! mt_mass first; if it would underrun zero, spill the residual
            ! into pending so the visible-side pool stays non-negative without
            ! double-counting.
            vis_total = mt_mass(imt, icell) + vis_pending
            IF (vis_vol_bef > mt_tiny) THEN
               ratio = vis_total / vis_vol_bef
               IF (mt_has_dissolved_limit(imt)) ratio = min(ratio, materials(imt)%max_dissolved_conc)
            ELSE
               ratio = 0._r8
            ENDIF
            mt_move = d_lev * ratio
            ! Cap by total visible mass (signed) so we never transfer more
            ! than the donor pool actually contains.
            IF (mt_move > 0._r8) THEN
               mt_move = min(mt_move, max(vis_total, 0._r8))
            ELSE
               mt_move = max(mt_move, min(vis_total, 0._r8))
            ENDIF
            ! Apportion the debit: take from mt_mass up to its current value,
            ! spill any residual into pending (only if pending is provided).
            IF (mt_move > 0._r8) THEN
               debit_mass    = min(mt_move, max(mt_mass(imt, icell), 0._r8))
               debit_pending = mt_move - debit_mass
            ELSE
               debit_mass    = max(mt_move, min(mt_mass(imt, icell), 0._r8))
               debit_pending = mt_move - debit_mass
            ENDIF
            mt_mass  (imt, icell) = mt_mass  (imt, icell) - debit_mass
            mt_levsto(imt, icell) = mt_levsto(imt, icell) + mt_move
            IF (present(pending_mt_pool) .and. abs(debit_pending) > 0._r8) THEN
               IF (imt <= size(pending_mt_pool)) THEN
                  pending_mt_pool(imt) = pending_mt_pool(imt) - debit_pending
               ENDIF
            ENDIF
         ELSE
            ! Protected -> visible. Protected side has no pending runoff input,
            ! so the original ratio = mt_levsto / levsto_bef is correct.
            IF (levsto_bef > mt_tiny) THEN
               ratio = mt_levsto(imt, icell) / levsto_bef
            ELSE
               ratio = 0._r8
            ENDIF
            mt_move = abs(d_lev) * ratio
            IF (mt_move > 0._r8) THEN
               mt_move = min(mt_move, max(mt_levsto(imt, icell), 0._r8))
            ELSE
               mt_move = max(mt_move, min(mt_levsto(imt, icell), 0._r8))
            ENDIF
            mt_levsto(imt, icell) = mt_levsto(imt, icell) - mt_move
            mt_mass  (imt, icell) = mt_mass  (imt, icell) + mt_move
         ENDIF
      ENDDO

      CALL mt_riverlake_equilibrate_cell(icell, vis_vol_aft, levsto_aft)

   END SUBROUTINE mt_riverlake_levee_repartition


   !-------------------------------------------------------------------------------------
   ! Refresh concentration from single-pool mass and current water volume.
   !-------------------------------------------------------------------------------------
   SUBROUTINE update_mt_concentration (imt, icell, volwater)

   IMPLICIT NONE

   integer,  intent(in) :: imt, icell
   real(r8), intent(in) :: volwater

      ! Intentionally does NOT clamp mt_mass here, even for tiny
      ! negatives (e.g. -1e-15 R*m³). Clamping state silently would
      ! swallow any future transport / flux-limiter bug that drifts
      ! mt_mass negative progressively — each substep would re-zero
      ! the residual below an absolute threshold and hide the drift.
      ! Instead, display-side reports (check_mt_state) dust-clamp
      ! for log cleanliness, and report a WARNING when magnitudes
      ! exceed the FP-noise floor. See material_mass_fp_dust in
      ! check_mt_state below for the cosmetic threshold.
      IF (volwater <= mt_v_dry_off) THEN
         mt_conc(imt, icell) = 0._r8
      ELSE
         mt_conc(imt, icell) = mt_mass(imt, icell) / volwater
      ENDIF

   END SUBROUTINE update_mt_concentration



   !-------------------------------------------------------------------------------------
   ! Accumulate material history diagnostics for one routing substep using the final
   ! post-update water state and any post-transport flux corrections already applied
   ! by the flow solver (for example inland-depression overflow adjustments).
   !-------------------------------------------------------------------------------------
   SUBROUTINE mt_riverlake_diag_accumulate (dt_all, irivsys, ucatfilter, wdsrf, volresv_in, &
      ucat2resv_in, is_built_resv, hflux_fc_in, allow_reverse_flow, bif_hflux_lev_in)

   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_next, push_ups2ucat, push_next2ucat
   USE MOD_Grid_RiverLakeLevee,   only: has_levee, levsto
   USE MOD_WorkerPushData
   USE MOD_MT_Defs, only: mt_tiny
   IMPLICIT NONE

   real(r8), intent(in) :: dt_all(:)
   integer,  intent(in) :: irivsys(:)
   logical,  intent(in) :: ucatfilter(:)
   real(r8), intent(in) :: wdsrf(:)
   real(r8), intent(in) :: volresv_in(:)
   integer,  intent(in) :: ucat2resv_in(:)
   logical,  intent(in), optional :: is_built_resv(:)
   real(r8), intent(in), optional :: hflux_fc_in(:)  ! downstream face flux [m3/s]
   logical,  intent(in), optional :: allow_reverse_flow  ! actual main-channel reverse edges

   real(r8), intent(in), optional :: bif_hflux_lev_in(:,:)

   integer :: i, imt
   real(r8) :: dt_i, volwater, dry_drain
   real(r8), allocatable :: fwd_send(:), fwd_recv(:), rev_send(:), rev_recv(:)
   real(r8), allocatable :: out_total(:), up_rev(:), donor_out(:), rev_sum(:), bif_out(:)
   logical :: reverse_allowed
   real(r8) :: vol_j, dt_j, rev_share, dt_don
   real(r8), allocatable :: bif_dn_levee(:), dry_coeff(:), bif_add(:), bif_add_lev(:), bif_hist(:)
   real(r8), allocatable :: received_visible(:)
   real(r8) :: pending_receipt
   logical :: bif_active
   real(r8), allocatable :: active_flux(:)

      IF (nmts <= 0) RETURN
      IF (.not. p_is_worker) RETURN
      bif_active = present(bif_hflux_lev_in)
      IF (numucat <= 0 .and. .not. bif_active) RETURN

      IF (allocated(mt_dry_fwd)) mt_dry_fwd = 0._r8
      IF (allocated(mt_dry_rev)) mt_dry_rev = 0._r8
      reverse_allowed = .false.
      IF (present(allow_reverse_flow)) reverse_allowed = allow_reverse_flow
      IF (.not. present(hflux_fc_in)) reverse_allowed = .false.

      ! A cell is only visited by the substeps of its own routing group.  A
      ! group that has already finished this timestep keeps the face fluxes of
      ! its last substep, so those cells must neither be dry-cleaned again nor
      ! contribute those stale faces to any split below: active_flux is the
      ! per-call copy that is zero where the group is not active.
      allocate(active_flux(numucat))
      active_flux = 0._r8

      DO i = 1, numucat
         IF (.not. ucatfilter(i)) CYCLE
         IF (irivsys(i) <= 0 .or. irivsys(i) > size(dt_all)) CYCLE
         dt_i = dt_all(irivsys(i))
         IF (dt_i <= 0._r8) CYCLE
         IF (present(hflux_fc_in)) active_flux(i) = hflux_fc_in(i)
         IF (present(is_built_resv)) THEN
            CALL get_cell_volume(i, wdsrf(i), volresv_in, ucat2resv_in, volwater, is_built_resv(i))
         ELSE
            CALL get_cell_volume(i, wdsrf(i), volresv_in, ucat2resv_in, volwater)
         ENDIF
         IF (allocated(mt_solid)) THEN
            IF (allocated(levsto) .and. allocated(has_levee)) THEN
               IF (has_levee(i)) THEN
                  CALL mt_riverlake_equilibrate_cell(i, volwater, levsto(i))
               ELSE
                  CALL mt_riverlake_equilibrate_cell(i, volwater, 0._r8)
               ENDIF
            ELSE
               CALL mt_riverlake_equilibrate_cell(i, volwater, 0._r8)
            ENDIF
         ENDIF
         IF (volwater <= mt_v_dry_off) THEN
            DO imt = 1, nmts
               IF (.not. mt_uses_generic_water(imt)) CYCLE
               ! Dry cleanup is a budgeted sink for non-negative physical mass,
               ! not a way to erase an invalid state before the watchdog sees it.
               IF (.not. ieee_is_finite(mt_mass(imt, i))) THEN
                  CALL CoLM_stop('non-finite river material mass reached dry-cell cleanup')
               ELSEIF (mt_mass(imt, i) < -MT_RESTART_NEGATIVE_DUST) THEN
                  CALL CoLM_stop('negative river material mass reached dry-cell cleanup')
               ELSEIF (mt_mass(imt, i) < 0._r8) THEN
                  mt_mass(imt, i) = 0._r8
               ENDIF
               IF (.not. ieee_is_finite(mt_inp_buf(imt, i))) &
                  CALL CoLM_stop('non-finite signed river material buffer reached dry-cell cleanup')
               IF (allocated(mt_solid)) THEN
                  IF (mt_has_dissolved_limit(imt)) THEN
                     mt_solid(imt, i) = mt_solid(imt, i) + &
                        max(mt_mass(imt, i), 0._r8) + max(mt_inp_buf(imt, i), 0._r8)
                     mt_mass(imt, i) = 0._r8
                     mt_inp_buf(imt, i) = min(mt_inp_buf(imt, i), 0._r8)
                     CYCLE
                  ENDIF
               ENDIF
               ! NEG_RUNOFF_DEBT: a negative pending runoff material is a
               ! signed correction tied to future same-cell runoff input, not
               ! a physical negative river outflow.  Dry-cell cleanup drains
               ! only positive orphan mass/pending input and preserves a
               ! negative mt_inp_buf debt so the next positive runoff can
               ! cancel it with the correct source signature.
               dry_drain = max(mt_mass(imt, i), 0._r8) + max(mt_inp_buf(imt, i), 0._r8)
               IF (dry_drain > mt_tiny) THEN
                  ! A terminal can also drain toward upstream cells through
                  ! reverse faces. Split every residual over actual outgoing
                  ! faces; only its terminal forward share leaves the domain.
                  ! With no outgoing face, the split keeps pending inventory.
                  IF (present(hflux_fc_in)) THEN
                     ! N07: hand the drained mass to the split delivery below.
                     ! It is distributed over EVERY actual outflow edge during
                     ! this substep, weighted by that edge's water flux:
                     !   forward share = max(hflux(i),0)/out_total(i)
                     !   reverse share_j = |hflux(j)|/out_total(i),  next(j)=i
                     ! with out_total(i) = max(hflux(i),0) + sum_upstream(-hflux).
                     ! This conserves mass and matches the host flux split also
                     ! for confluences (several upstream branches).
                     mt_dry_split(imt, i) = mt_dry_split(imt, i) + dry_drain
                     IF (allocated(mt_dry_stat)) mt_dry_stat(imt, 1) = mt_dry_stat(imt, 1) + dry_drain
                  ELSE
                     IF (allocated(mt_dry_drain)) mt_dry_drain(imt, i) = &
                        mt_dry_drain(imt, i) + dry_drain
                     IF (allocated(mt_dry_stat)) mt_dry_stat(imt, 4) = mt_dry_stat(imt, 4) + dry_drain
                     ! Keep advective outflux and the terminal dry-cell sink
                     ! separate. CoLMDEBUG books mt_flux_out at river mouths
                     ! and adds mt_dry_drain globally; folding the sink into
                     ! both arrays counted a dry mouth twice.
                     IF (dt_i > 0._r8 .and. ucatfilter(i)) &
                        a_mt_out(imt, i) = a_mt_out(imt, i) + dry_drain
                  ENDIF
               ENDIF
               mt_mass(imt, i)    = 0._r8
               mt_inp_buf(imt, i) = min(mt_inp_buf(imt, i), 0._r8)
            ENDDO
         ENDIF
         DO imt = 1, nmts
            IF (.not. mt_uses_generic_water(imt)) CYCLE
            CALL update_mt_concentration(imt, i, volwater)
         ENDDO
         a_mt_acctime(i) = a_mt_acctime(i) + dt_i
         IF (volwater > mt_v_dry_off) THEN
            IF (allocated(a_water_storage)) a_water_storage(i) = a_water_storage(i) + volwater * dt_i
         ENDIF
         IF (allocated(a_levsto_water) .and. allocated(levsto) .and. allocated(has_levee)) THEN
            IF (i <= size(levsto) .and. i <= size(has_levee)) THEN
               IF (has_levee(i) .and. levsto(i) > mt_v_dry_off) THEN
                  a_levsto_water(i) = a_levsto_water(i) + levsto(i) * dt_i
               ENDIF
            ENDIF
         ENDIF

         DO imt = 1, nmts
            IF (.not. mt_uses_generic_water(imt)) CYCLE
            IF (volwater > mt_v_dry_off .and. allocated(a_mt_storage_mass)) THEN
               a_mt_storage_mass(imt, i) = a_mt_storage_mass(imt, i) &
                  + mt_mass(imt, i) * dt_i
            ENDIF
            IF (allocated(a_mt_levsto_mass) .and. allocated(mt_levsto) .and. &
                allocated(levsto) .and. allocated(has_levee)) THEN
               IF (i <= size(mt_levsto, 2) .and. i <= size(levsto) .and. i <= size(has_levee)) THEN
                  IF (has_levee(i) .and. levsto(i) > mt_v_dry_off) THEN
                     a_mt_levsto_mass(imt, i) = a_mt_levsto_mass(imt, i) &
                        + mt_levsto(imt, i) * dt_i
                  ENDIF
               ENDIF
            ENDIF
            a_mt_out   (imt, i) = a_mt_out   (imt, i) + mt_flux_out(imt, i) * dt_i
            a_mt_bifout(imt, i) = a_mt_bifout(imt, i) + mt_bif_net_saved(imt, i) * dt_i
         ENDDO
      ENDDO

      IF (allocated(mt_dry_split) .and. present(hflux_fc_in) .and. (numucat > 0 .or. bif_active)) THEN
         allocate(fwd_send(numucat), fwd_recv(numucat), rev_send(numucat), rev_recv(numucat), received_visible(numucat))
         allocate(out_total(numucat), up_rev(numucat), donor_out(numucat), rev_sum(numucat), &
            bif_out(numucat))
         ! reverse-outflow demand of every donor: sum of |negative upstream face
         ! fluxes|.  `push_ups2ucat` aggregates upstream-cell values at the
         ! downstream cell, so this is available locally after one push.
         IF (reverse_allowed) THEN
            up_rev(:) = max(-active_flux(:), 0._r8)
         ELSE
            ! The caller explicitly disabled main-channel reverse delivery.
            up_rev(:) = 0._r8
         ENDIF
         rev_sum(:) = 0._r8
         CALL worker_push_data (push_ups2ucat, up_rev, rev_sum, fillvalue = 0._r8, mode = 'sum')
         bif_out = 0._r8
         IF (bif_active) THEN
            CALL mt_bif_dry_prepare(bif_hflux_lev_in, bif_out, bif_dn_levee)
            allocate(dry_coeff(numucat), bif_add(numucat), bif_add_lev(numucat), bif_hist(numucat))
         ENDIF
         DO i = 1, numucat
            out_total(i) = max(active_flux(i), 0._r8) + rev_sum(i) + bif_out(i)
         ENDDO

         ! A donor with no actual main-channel outflow has nowhere physical to
         ! deliver: keep the conservative pending state instead of losing it.
         DO i = 1, numucat
            IF (out_total(i) <= mt_tiny) THEN
               DO imt = 1, nmts
                  IF (mt_dry_split(imt, i) <= mt_tiny) CYCLE
                  IF (allocated(mt_dry_pending)) THEN
                     mt_dry_pending(imt, i) = mt_dry_pending(imt, i) + mt_dry_split(imt, i)
                  ENDIF
                  mt_dry_split(imt, i) = 0._r8
               ENDDO
            ENDIF
         ENDDO
         ! upstream receivers get their donor's total outflow capacity
         donor_out(:) = 0._r8
         CALL worker_push_data (push_next2ucat, out_total, donor_out, fillvalue = 0._r8)
         DO imt = 1, nmts
            IF (.not. mt_uses_generic_water(imt)) CYCLE
            ! N08: reset BOTH send buffers for every species.  A position
            ! where this species has no drained mass is skipped below, so a
            ! stale value from the previous species would otherwise be sent
            ! as this species' mass (cross-species contamination).
            received_visible = 0._r8
            fwd_send(:) = 0._r8
            rev_send(:) = 0._r8
            fwd_recv(:) = 0._r8
            rev_recv(:) = 0._r8
            DO i = 1, numucat
               IF (mt_dry_split(imt, i) <= mt_tiny) CYCLE
               IF (out_total(i) > mt_tiny) THEN
                  fwd_send(i) = mt_dry_split(imt, i) * max(active_flux(i), 0._r8) / out_total(i)
                  IF (ucat_next(i) <= 0) THEN
                     mt_dry_drain(imt,i) = mt_dry_drain(imt,i) + fwd_send(i)
                     IF (allocated(mt_dry_stat)) mt_dry_stat(imt,4) = mt_dry_stat(imt,4) + fwd_send(i)
                  ENDIF
                  ! The delivered mass leaves this cell through the drained-mass
                  ! path, so it belongs in the cell's outflow history as well.
                  ! Both parts are known here: the forward share is fwd_send,
                  ! and the reverse receivers take rev_sum/out_total of the
                  ! split.  Recorded as mass (the history accumulator is a
                  ! mass integral), so no division by dt is needed.
                  ! Only the FORWARD share belongs to this cell's face.  The
                  ! reverse share is delivered to upstream receivers and is
                  ! recorded there (see the receive loop), with the sign the
                  ! advective ledger uses for a reverse face flux (negative:
                  ! mass entering this cell from its downstream neighbour).
                  IF (p_is_worker .and. allocated(a_mt_out) .and. ucatfilter(i)) THEN
                     dt_don = 0._r8
                     IF (irivsys(i) > 0 .and. irivsys(i) <= size(dt_all)) dt_don = dt_all(irivsys(i))
                     IF (dt_don > 0._r8) a_mt_out(imt, i) = a_mt_out(imt, i) + fwd_send(i)
                  ENDIF
               ENDIF
               ! rev_send carries the WHOLE split, so it must stay zero when
               ! reverse flow is disabled: the forward share is already the
               ! entire drained mass (out_total = max(hflux,0)), and a single
               ! negative hflux downstream would deliver it a second time.
               IF (reverse_allowed) rev_send(i) = mt_dry_split(imt, i)
            ENDDO
            CALL worker_push_data (push_ups2ucat, fwd_send, fwd_recv, fillvalue = 0._r8, mode = 'sum')
            CALL worker_push_data (push_next2ucat, rev_send, rev_recv, fillvalue = 0._r8)
            DO i = 1, numucat
               IF (fwd_recv(i) <= mt_tiny .and. rev_recv(i) <= mt_tiny) CYCLE
               IF (allocated(mt_dry_stat)) mt_dry_stat(imt, 2) = &
                  mt_dry_stat(imt, 2) + max(fwd_recv(i), 0._r8)
               IF (present(is_built_resv)) THEN
                  CALL get_cell_volume(i, wdsrf(i), volresv_in, ucat2resv_in, vol_j, is_built_resv(i))
               ELSE
                  CALL get_cell_volume(i, wdsrf(i), volresv_in, ucat2resv_in, vol_j)
               ENDIF
               dt_j = 0._r8
               IF (irivsys(i) > 0 .and. irivsys(i) <= size(dt_all)) dt_j = dt_all(irivsys(i))
               ! forward mass arrives as an upstream sum
               IF (fwd_recv(i) > mt_tiny) THEN
                  mt_mass(imt, i) = mt_mass(imt, i) + fwd_recv(i)
                  received_visible(i) = received_visible(i) + fwd_recv(i)
                  CALL update_mt_concentration(imt, i, vol_j)
                  IF (vol_j > mt_v_dry_off .and. allocated(a_mt_storage_mass) .and. dt_j > 0._r8) &
                     a_mt_storage_mass(imt, i) = a_mt_storage_mass(imt, i) + fwd_recv(i) * dt_j
               ENDIF
               ! reverse share: receiver pulls its donor's drained mass and
               ! weights it by its own reverse edge flux (known locally)
               IF (reverse_allowed .and. rev_recv(i) > mt_tiny .and. &
                   active_flux(i) < 0._r8 .and. donor_out(i) > mt_tiny) THEN
                  rev_share = rev_recv(i) * abs(active_flux(i)) / donor_out(i)
                  IF (allocated(mt_dry_stat)) mt_dry_stat(imt, 2) = mt_dry_stat(imt, 2) + rev_share
                  IF (allocated(a_mt_out) .and. ucatfilter(i)) THEN
                     dt_j = 0._r8
                     IF (irivsys(i) > 0 .and. irivsys(i) <= size(dt_all)) dt_j = dt_all(irivsys(i))
                     IF (dt_j > 0._r8) a_mt_out(imt, i) = a_mt_out(imt, i) - rev_share
                  ENDIF
                  mt_mass(imt, i) = mt_mass(imt, i) + rev_share
                  received_visible(i) = received_visible(i) + rev_share
                  CALL update_mt_concentration(imt, i, vol_j)
                  IF (vol_j > mt_v_dry_off .and. allocated(a_mt_storage_mass) .and. dt_j > 0._r8) &
                     a_mt_storage_mass(imt, i) = a_mt_storage_mass(imt, i) + rev_share * dt_j
               ENDIF
            ENDDO
            IF (bif_active) THEN
               dry_coeff = 0._r8
               WHERE (out_total > mt_tiny) dry_coeff = mt_dry_split(imt,:) / max(out_total, tiny(1._r8))
               CALL mt_bif_dry_deliver(bif_hflux_lev_in, bif_dn_levee, dry_coeff, &
                  bif_add, bif_add_lev, bif_hist)
               received_visible = received_visible + bif_add
               mt_mass(imt,:) = mt_mass(imt,:) + bif_add
               mt_levsto(imt,:) = mt_levsto(imt,:) + bif_add_lev
               IF (allocated(mt_dry_stat)) mt_dry_stat(imt,2) = mt_dry_stat(imt,2) + sum(bif_add+bif_add_lev)
               DO i = 1, numucat
                  dt_j = 0._r8
                  IF (irivsys(i) > 0 .and. irivsys(i) <= size(dt_all)) dt_j = dt_all(irivsys(i))
                  IF (.not. ucatfilter(i) .or. dt_j <= 0._r8) CYCLE
                  a_mt_bifout(imt,i) = a_mt_bifout(imt,i) + bif_hist(i)
                  IF (present(is_built_resv)) THEN
                     CALL get_cell_volume(i, wdsrf(i), volresv_in, ucat2resv_in, vol_j, is_built_resv(i))
                  ELSE
                     CALL get_cell_volume(i, wdsrf(i), volresv_in, ucat2resv_in, vol_j)
                  ENDIF
                  CALL update_mt_concentration(imt, i, vol_j)
                  IF (vol_j > mt_v_dry_off) &
                     a_mt_storage_mass(imt,i) = a_mt_storage_mass(imt,i) + bif_add(i)*dt_j
                  IF (allocated(levsto)) THEN
                     IF (levsto(i) > mt_v_dry_off) &
                        a_mt_levsto_mass(imt,i) = a_mt_levsto_mass(imt,i) + bif_add_lev(i)*dt_j
                  ENDIF
               ENDDO
            ENDIF
            ! A dry receiver has no carrier concentration. Keep this already
            ! booked transfer in explicit inventory until a later active step;
            ! do not chain it through another edge with the same water flux.
            DO i = 1, numucat
               IF (received_visible(i) <= mt_tiny) CYCLE
               IF (present(is_built_resv)) THEN
                  CALL get_cell_volume(i, wdsrf(i), volresv_in, ucat2resv_in, vol_j, is_built_resv(i))
               ELSE
                  CALL get_cell_volume(i, wdsrf(i), volresv_in, ucat2resv_in, vol_j)
               ENDIF
               IF (vol_j > mt_v_dry_off) CYCLE
               pending_receipt = min(received_visible(i), mt_mass(imt,i))
               mt_dry_pending(imt,i) = mt_dry_pending(imt,i) + pending_receipt
               mt_mass(imt,i) = mt_mass(imt,i) - pending_receipt
               CALL update_mt_concentration(imt, i, vol_j)
            ENDDO
            mt_dry_split(imt, :) = 0._r8
         ENDDO
         deallocate(fwd_send, fwd_recv, rev_send, rev_recv)
         deallocate(out_total, up_rev, donor_out, rev_sum, bif_out)
      ENDIF

   END SUBROUTINE mt_riverlake_diag_accumulate


   ! Only visible dry mass is extracted by diag_accumulate. Protected donors
   ! therefore contribute no capacity; receivers may be in either water pool.
   SUBROUTINE mt_bif_dry_prepare(q, gross, down_levee)
      USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
      USE MOD_Grid_RiverLakeNetwork, only: numucat, pth_upst_local, push_bif_dn2pth, push_bif_influx
      USE MOD_Grid_RiverLakeLevee, only: has_levee
      USE MOD_WorkerPushData
      real(r8), intent(in) :: q(:,:)
      real(r8), intent(out) :: gross(:)
      real(r8), allocatable, intent(out) :: down_levee(:)
      real(r8) :: flags(numucat), reverse_path(size(q,2)), reverse_recv(numucat)
      integer :: p, k, u
      logical :: up_protected, dn_protected

      IF (any(.not. ieee_is_finite(q))) CALL CoLM_stop('nonfinite dry BIF water flux')
      IF (size(q,2) /= size(pth_upst_local)) CALL CoLM_stop('dry BIF pathway shape mismatch')
      allocate(down_levee(size(q,2)))
      flags = 0._r8
      IF (DEF_USE_LEVEE .and. allocated(has_levee)) WHERE (has_levee) flags = 1._r8
      CALL worker_push_data(push_bif_dn2pth, flags, down_levee, fillvalue=0._r8)
      gross = 0._r8
      reverse_path = 0._r8
      DO p = 1, size(q,2)
         u = pth_upst_local(p)
         IF (u < 1 .or. u > numucat) CALL CoLM_stop('dry BIF invalid upstream endpoint')
         DO k = 1, size(q,1)
            up_protected = k > 1 .and. flags(u) > 0.5_r8
            dn_protected = k > 1 .and. down_levee(p) > 0.5_r8
            IF (.not. up_protected) gross(u) = gross(u) + max(q(k,p),0._r8)
            IF (.not. dn_protected) reverse_path(p) = reverse_path(p) + max(-q(k,p),0._r8)
         ENDDO
      ENDDO
      CALL worker_push_data(push_bif_influx, reverse_path, reverse_recv, fillvalue=0._r8, mode='sum')
      gross = gross + reverse_recv
   END SUBROUTINE mt_bif_dry_prepare

   SUBROUTINE mt_bif_dry_deliver(q, down_levee, coefficient, visible_add, protected_add, history_net)
      USE MOD_Grid_RiverLakeNetwork, only: numucat, pth_upst_local, push_bif_dn2pth, push_bif_influx
      USE MOD_Grid_RiverLakeLevee, only: has_levee
      USE MOD_WorkerPushData
      real(r8), intent(in) :: q(:,:), down_levee(:), coefficient(:)
      real(r8), intent(out) :: visible_add(:), protected_add(:), history_net(:)
      real(r8) :: down_coeff(size(q,2)), vis_path(size(q,2)), prot_path(size(q,2))
      real(r8) :: vis_recv(numucat), prot_recv(numucat), hist_path(size(q,2)), hist_recv(numucat)
      real(r8) :: mass
      integer :: p, k, u
      logical :: up_levee, up_protected, dn_protected

      CALL worker_push_data(push_bif_dn2pth, coefficient, down_coeff, fillvalue=0._r8)
      visible_add = 0._r8
      protected_add = 0._r8
      history_net = 0._r8
      vis_path = 0._r8
      prot_path = 0._r8
      hist_path = 0._r8
      DO p = 1, size(q,2)
         u = pth_upst_local(p)
         up_levee = .false.
         IF (DEF_USE_LEVEE .and. allocated(has_levee)) up_levee = has_levee(u)
         DO k = 1, size(q,1)
            up_protected = k > 1 .and. up_levee
            dn_protected = k > 1 .and. down_levee(p) > 0.5_r8
            mass = 0._r8
            IF (q(k,p) >= 0._r8) THEN
               IF (.not. up_protected) mass = q(k,p)*coefficient(u)
               IF (dn_protected) THEN
                  prot_path(p) = prot_path(p) + mass
               ELSE
                  vis_path(p) = vis_path(p) + mass
               ENDIF
            ELSE
               IF (.not. dn_protected) mass = q(k,p)*down_coeff(p)
               IF (up_protected) THEN
                  protected_add(u) = protected_add(u) - mass
               ELSE
                  visible_add(u) = visible_add(u) - mass
               ENDIF
            ENDIF
            history_net(u) = history_net(u) + mass
            hist_path(p) = hist_path(p) + mass
         ENDDO
      ENDDO
      CALL worker_push_data(push_bif_influx, vis_path, vis_recv, fillvalue=0._r8, mode='sum')
      CALL worker_push_data(push_bif_influx, prot_path, prot_recv, fillvalue=0._r8, mode='sum')
      CALL worker_push_data(push_bif_influx, hist_path, hist_recv, fillvalue=0._r8, mode='sum')
      visible_add = visible_add + vis_recv
      protected_add = protected_add + prot_recv
      history_net = history_net - hist_recv
   END SUBROUTINE mt_bif_dry_deliver

   SUBROUTINE ensure_mt_substep_workspace (nucat, npth, nlev_bif, bif_active)
      IMPLICIT NONE
      integer, intent(in) :: nucat, npth, nlev_bif
      logical, intent(in) :: bif_active

      CALL ensure_real1_workspace(conc_next, nucat)
      CALL ensure_real1_workspace(flux_ups, nucat)
      CALL ensure_real1_workspace(mt_flux, nucat)
      CALL ensure_real1_workspace(mt_conc_flux, nucat)
      CALL ensure_real1_workspace(mt_prot_conc_flux, nucat)
      CALL ensure_real1_workspace(bif_net, nucat)
      CALL ensure_real1_workspace(mt_bif_lev_net, nucat)
      CALL ensure_real1_workspace(mt_in_mass, nucat)
      CALL ensure_real1_workspace(mt_in_mass_lev, nucat)
      CALL ensure_real1_workspace(mt_out_mass, nucat)
      CALL ensure_real1_workspace(mt_out_mass_lev, nucat)
      CALL ensure_real1_workspace(rate_cell, nucat)
      CALL ensure_real1_workspace(rate_cell_lev, nucat)
      CALL ensure_real1_workspace(rate_next, nucat)
      CALL ensure_real1_workspace(mt_inp_step, nucat)

      IF (bif_active) THEN
         CALL ensure_real1_workspace(conc_dn_pth, npth)
         CALL ensure_real1_workspace(mt_prot_conc_dn_pth, npth)
         CALL ensure_real1_workspace(mt_pth_single, npth)
         CALL ensure_real1_workspace(mt_pth_lev, npth)
         CALL ensure_real2_workspace(mt_pth_levee, nlev_bif, npth)
         CALL ensure_real1_workspace(bif_recv, nucat)
         CALL ensure_real1_workspace(bif_lev_recv, nucat)
         CALL ensure_real1_workspace(rate_dn_pth, npth)
         CALL ensure_real1_workspace(rate_dn_pth_lev, npth)
         CALL ensure_real1_workspace(has_levee_r8, nucat)
         CALL ensure_real1_workspace(has_levee_dn_pth, npth)
         CALL ensure_real1_workspace(dt_ucat, nucat)
         CALL ensure_real1_workspace(dt_dn_pth, npth)
         CALL ensure_real1_workspace(mt_dn_out_vis_pth, npth)
         CALL ensure_real1_workspace(mt_dn_out_lev_pth, npth)
         CALL ensure_real1_workspace(mt_dn_out_vis_recv, nucat)
         CALL ensure_real1_workspace(mt_dn_out_lev_recv, nucat)
      ENDIF

   END SUBROUTINE ensure_mt_substep_workspace


   SUBROUTINE ensure_real1_workspace (buf, n)
      IMPLICIT NONE
      real(r8), allocatable, intent(inout) :: buf(:)
      integer, intent(in) :: n

      IF (allocated(buf)) THEN
         IF (size(buf) /= n) deallocate(buf)
      ENDIF
      IF (.not. allocated(buf)) allocate(buf(n))

   END SUBROUTINE ensure_real1_workspace


   SUBROUTINE ensure_real2_workspace (buf, n1, n2)
      IMPLICIT NONE
      real(r8), allocatable, intent(inout) :: buf(:,:)
      integer, intent(in) :: n1, n2

      IF (allocated(buf)) THEN
         IF (size(buf, 1) /= n1 .or. size(buf, 2) /= n2) deallocate(buf)
      ENDIF
      IF (.not. allocated(buf)) allocate(buf(n1, n2))

   END SUBROUTINE ensure_real2_workspace


   SUBROUTINE material_bif_path_levee_sides (can_use_levee_material, i_up, i_dn, ipth, nucat, &
      upstream_has_levee, downstream_has_levee)
      USE MOD_Grid_RiverLakeLevee, only: has_levee
      IMPLICIT NONE
      logical, intent(in) :: can_use_levee_material
      integer, intent(in) :: i_up, i_dn, ipth, nucat
      logical, intent(out) :: upstream_has_levee, downstream_has_levee

      upstream_has_levee = .false.
      downstream_has_levee = .false.
      IF (.not. can_use_levee_material .or. .not. allocated(has_levee)) RETURN

      IF (i_up > 0 .and. i_up <= size(has_levee)) upstream_has_levee = has_levee(i_up)
      IF (i_dn > 0 .and. i_dn <= nucat) THEN
         IF (i_dn <= size(has_levee)) downstream_has_levee = has_levee(i_dn)
      ELSEIF (allocated(has_levee_dn_pth)) THEN
         IF (ipth > 0 .and. ipth <= size(has_levee_dn_pth)) &
            downstream_has_levee = has_levee_dn_pth(ipth) > 0.5_r8
      ENDIF
   END SUBROUTINE material_bif_path_levee_sides


   SUBROUTINE release_mt_substep_workspace ()
      IMPLICIT NONE

      IF (allocated(conc_next             )) deallocate(conc_next             )
      IF (allocated(flux_ups              )) deallocate(flux_ups              )
      IF (allocated(mt_flux              )) deallocate(mt_flux              )
      IF (allocated(mt_conc_flux         )) deallocate(mt_conc_flux         )
      IF (allocated(mt_prot_conc_flux    )) deallocate(mt_prot_conc_flux    )
      IF (allocated(conc_dn_pth           )) deallocate(conc_dn_pth           )
      IF (allocated(mt_prot_conc_dn_pth  )) deallocate(mt_prot_conc_dn_pth  )
      IF (allocated(mt_pth_single          )) deallocate(mt_pth_single          )
      IF (allocated(mt_pth_lev           )) deallocate(mt_pth_lev           )
      IF (allocated(mt_pth_levee        )) deallocate(mt_pth_levee        )
      IF (allocated(bif_recv              )) deallocate(bif_recv              )
      IF (allocated(bif_lev_recv          )) deallocate(bif_lev_recv          )
      IF (allocated(bif_net               )) deallocate(bif_net               )
      IF (allocated(mt_bif_lev_net       )) deallocate(mt_bif_lev_net       )
      IF (allocated(has_levee_r8          )) deallocate(has_levee_r8          )
      IF (allocated(has_levee_dn_pth      )) deallocate(has_levee_dn_pth      )
      IF (allocated(dt_ucat               )) deallocate(dt_ucat               )
      IF (allocated(dt_dn_pth             )) deallocate(dt_dn_pth             )
      IF (allocated(mt_in_mass           )) deallocate(mt_in_mass           )
      IF (allocated(mt_in_mass_lev       )) deallocate(mt_in_mass_lev       )
      IF (allocated(mt_out_mass          )) deallocate(mt_out_mass          )
      IF (allocated(mt_out_mass_lev      )) deallocate(mt_out_mass_lev      )
      IF (allocated(rate_cell             )) deallocate(rate_cell             )
      IF (allocated(rate_cell_lev         )) deallocate(rate_cell_lev         )
      IF (allocated(rate_next             )) deallocate(rate_next             )
      IF (allocated(rate_dn_pth           )) deallocate(rate_dn_pth           )
      IF (allocated(rate_dn_pth_lev       )) deallocate(rate_dn_pth_lev       )
      IF (allocated(mt_inp_step          )) deallocate(mt_inp_step          )
      IF (allocated(mt_dn_out_vis_pth    )) deallocate(mt_dn_out_vis_pth    )
      IF (allocated(mt_dn_out_lev_pth    )) deallocate(mt_dn_out_lev_pth    )
      IF (allocated(mt_dn_out_vis_recv   )) deallocate(mt_dn_out_vis_recv   )
      IF (allocated(mt_dn_out_lev_recv   )) deallocate(mt_dn_out_lev_recv   )

   END SUBROUTINE release_mt_substep_workspace


   !-------------------------------------------------------------------------------------
   ! Per-sub-step material transport (called inside the DO WHILE routing loop).
   ! Uses instantaneous water fluxes, not time-averaged, so material and water
   ! advance in lockstep.  This avoids the "water left but material stayed"
   ! artefact of the old single-step-per-period approach.
   !-------------------------------------------------------------------------------------
   SUBROUTINE mt_riverlake_substep (dt_ref, dt_all, irivsys, hflux_fc, sum_hflux_riv, wdsrf, ucatfilter, &
      volresv, ucat2resv, is_built_resv, &
      do_bif, bif_hflux_lev_in, npthout_local_in)

   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_next, &
      floodplain_curve, lake_type, push_ups2ucat, push_next2ucat, &
      npthlev_bif, pth_upst_local, pth_down_local, &
      push_bif_influx, push_bif_dn2pth, rivsys_by_multiple_procs
#ifdef USEMPI
   USE MOD_Grid_RiverLakeNetwork, only: p_comm_rivsys
#endif
   USE MOD_Grid_RiverLakeLevee, only: has_levee, levsto
   USE MOD_Grid_RiverLakeTimeVars, only: volwater_ucat
   USE MOD_WorkerPushData
   USE MOD_MT_Defs, only: mt_tiny, materials, material_init_water_ratio
   USE MOD_MT_Solute, only: material_decay_fraction
   IMPLICIT NONE

   real(r8), intent(in) :: dt_ref
   real(r8), intent(in) :: dt_all(:)
   integer,  intent(in) :: irivsys(:)
   real(r8), intent(in) :: hflux_fc(:)
   real(r8), intent(in) :: sum_hflux_riv(:)
   real(r8), intent(in) :: wdsrf(:)
   logical,  intent(in) :: ucatfilter(:)
   real(r8), intent(in) :: volresv(:)
   integer,  intent(in) :: ucat2resv(:)
   logical,  intent(in) :: is_built_resv(:)
   logical,  intent(in), optional :: do_bif
   real(r8), intent(in), optional :: bif_hflux_lev_in(:,:)
   integer,  intent(in), optional :: npthout_local_in

   integer  :: i, imt, ipth, i_up, i_dn, ilev, limiter_iter, limiter_max_iter
   integer  :: limiter_iters_used
   real(r8) :: volwater, volflux, dt_i, dt_donor, mt_pth_fl
   real(r8) :: layer_wflux, mt_rate, limiter_rate_new
   real(r8) :: limiter_delta, limiter_delta_global
   logical  :: upstream_has_levee, downstream_has_levee, can_use_levee_material
   logical  :: limiter_converged
   real(r8) :: release, R_fill
   real(r8) :: mt_mass_new, levsto_old
   real(r8) :: decay_fraction, decay_src
   logical  :: bif_workspace_active
   integer  :: npth_bif, nlev_bif

      ! Initialisation intentionally leaves all river material arrays unallocated
      ! when DEF_MT_NUM=0. Skip transport on every rank before workspace or
      ! per-cell state can be touched (including the BIF collective path).
      IF (nmts <= 0) RETURN

      npth_bif = 0
      nlev_bif = 0
      bif_workspace_active = .false.
      IF (present(do_bif) .and. present(bif_hflux_lev_in) .and. present(npthout_local_in)) THEN
         npth_bif = npthout_local_in
         nlev_bif = size(bif_hflux_lev_in, 1)
         bif_workspace_active = do_bif
      ENDIF

      ! A Jacobi fixed point can move donor-rate information only one graph
      ! edge per iteration. The water-side BIF aggregate limiter reserves
      ! ordinary gross outflow before authorising BIF outflow, separately for
      ! visible/protected pools. With the same well-mixed concentration used
      ! here, every nonzero BIF sender therefore starts at material rate one;
      ! sub-unit-rate dependencies propagate only along the ordinary river
      ! forest, so 2*N+1 is a conservative topology-derived cap on the number of
      ! rounds.  The loop exits on convergence
      ! (MT_LIMITER_RATE_TOL) and is fatal if the cap is exhausted; it is never
      ! truncated early, because the dry-cell sink (mt_dry_drain) would then
      ! remove the under-exported material instead of delivering it downstream.
      ! MT_LIMITER_SOFT_ITER is only a fast-path threshold used to report how
      ! often long chains occur.
      IF (totalnumucat < 0 .or. totalnumucat > (huge(limiter_max_iter) - 1) / 2) &
         CALL CoLM_stop('invalid river topology size for material donor limiter')
      limiter_max_iter = 2 * totalnumucat + 1

      IF (numucat <= 0 .and. npth_bif <= 0 .and. .not. bif_workspace_active) RETURN

      CALL ensure_mt_substep_workspace(numucat, npth_bif, nlev_bif, bif_workspace_active)

      can_use_levee_material = DEF_USE_LEVEE .and. allocated(has_levee) .and. &
         allocated(mt_levsto) .and. allocated(levsto)

      IF (bif_workspace_active) THEN
         has_levee_r8(:) = 0._r8
         IF (can_use_levee_material) THEN
            DO i = 1, numucat
               IF (i <= size(has_levee)) THEN
                  IF (has_levee(i)) has_levee_r8(i) = 1._r8
               ENDIF
            ENDDO
         ENDIF
         CALL worker_push_data (push_bif_dn2pth, has_levee_r8, has_levee_dn_pth, fillvalue = 0._r8)
      ENDIF
      IF (bif_workspace_active) THEN
         dt_ucat(:) = 0._r8
         DO i = 1, numucat
            IF (irivsys(i) > 0 .and. irivsys(i) <= size(dt_all)) dt_ucat(i) = dt_all(irivsys(i))
         ENDDO
         CALL worker_push_data (push_bif_dn2pth, dt_ucat, dt_dn_pth, fillvalue = 0._r8)
      ENDIF

         DO imt = 1, nmts
            IF (.not. mt_uses_generic_water(imt)) CYCLE
            R_fill = material_init_water_ratio(imt)
            IF (mt_has_dissolved_limit(imt)) R_fill = min(R_fill, materials(imt)%max_dissolved_conc)

         ! --- 1. Concentration from pre-update single-pool state ---
         DO i = 1, numucat
            mt_conc_flux(i) = 0._r8
            mt_prot_conc_flux(i) = 0._r8
            IF (.not. ucatfilter(i)) CYCLE
            dt_i = 0._r8
            IF (irivsys(i) > 0 .and. irivsys(i) <= size(dt_all)) dt_i = dt_all(irivsys(i))
            IF (dt_i <= 0._r8) CYCLE
            CALL get_cell_volume(i, wdsrf(i), volresv, ucat2resv, volwater, is_built_resv(i))
            ! The restart/state contract admits only sub-tolerance negative
            ! roundoff. Canonicalize it before forming a concentration: a
            ! negative flux would violate the limiter's nonnegative monotone
            ! map even if the underlying mass error is numerically tiny.
            IF (.not. ieee_is_finite(mt_mass(imt, i))) THEN
               CALL CoLM_stop('non-finite river material mass before transport')
            ELSEIF (mt_mass(imt, i) < -MT_RESTART_NEGATIVE_DUST) THEN
               CALL CoLM_stop('negative river material mass before transport')
            ELSEIF (mt_mass(imt, i) < 0._r8) THEN
               mt_mass(imt, i) = 0._r8
            ENDIF
            IF (allocated(mt_levsto)) THEN
               IF (i <= size(mt_levsto, 2)) THEN
                  IF (.not. ieee_is_finite(mt_levsto(imt, i))) THEN
                     CALL CoLM_stop('non-finite protected material mass before transport')
                  ELSEIF (mt_levsto(imt, i) < -MT_RESTART_NEGATIVE_DUST) THEN
                     CALL CoLM_stop('negative protected material mass before transport')
                  ELSEIF (mt_levsto(imt, i) < 0._r8) THEN
                     mt_levsto(imt, i) = 0._r8
                  ENDIF
               ENDIF
            ENDIF
            ! Release every positive pending input in the first substep,
            ! even when the pre-update cell is dry. Water can arrive during
            ! this same substep; if it remains dry, the post-update dry-drain
            ! path removes the orphan mass. Signed debt stays in the buffer
            ! to cancel later same-cell runoff rather than borrowing from
            ! the current river inventory.
            release = max(mt_inp_buf(imt, i), 0._r8)
            mt_inp_buf(imt, i) = mt_inp_buf(imt, i) - release
            mt_mass(imt, i) = mt_mass(imt, i) + release
            ! Pending mass from a reverse-flow dry cleanup: put it back in
            ! the cell inventory so it is transported once the host water
            ! returns; never deleted, never sent against the flow.
            IF (allocated(mt_dry_pending)) THEN
               IF (mt_dry_pending(imt, i) > mt_tiny) THEN
                  mt_mass(imt, i) = mt_mass(imt, i) + mt_dry_pending(imt, i)
                  mt_dry_pending(imt, i) = 0._r8
               ENDIF
            ENDIF
            IF (allocated(mt_solid)) THEN
               IF (mt_has_dissolved_limit(imt)) THEN
                  CALL mt_equilibrate_dissolved(imt, volwater, mt_mass(imt, i), mt_solid(imt, i))
                  IF (can_use_levee_material) THEN
                     IF (has_levee(i)) CALL mt_equilibrate_dissolved(imt, levsto(i), &
                        mt_levsto(imt, i), mt_levsto_solid(imt, i))
                  ENDIF
               ENDIF
            ENDIF
            CALL update_mt_concentration(imt, i, volwater)
            IF (volwater > mt_v_dry_off) THEN
               ! Use the actual well-mixed pool concentration. Inflating the
               ! denominator to an outgoing edge volume breaks the constant-
               ! state invariant whenever same-step inflow balances gross
               ! outflow; the coupled donor limiter now handles that demand.
               mt_conc_flux(i) = mt_mass(imt, i) / volwater
            ELSE
               ! A dry cell has no physical concentration. Retain the bounded
               ! edge-volume fallback only to transport explicitly queued mass
               ! without division by zero; the limiter remains authoritative.
               volflux = max(abs(hflux_fc(i)) * dt_i, mt_v_dry_off)
               mt_conc_flux(i) = mt_mass(imt, i) / volflux
               IF (mt_has_dissolved_limit(imt)) mt_conc_flux(i) = 0._r8
            ENDIF
            mt_prot_conc_flux(i) = mt_conc_flux(i)
            IF (can_use_levee_material) THEN
               IF (i <= size(has_levee) .and. i <= size(levsto) .and. i <= size(mt_levsto, 2)) THEN
                  IF (has_levee(i)) THEN
                     IF (.not. ieee_is_finite(levsto(i)) .or. levsto(i) < 0._r8) &
                        CALL CoLM_stop('invalid protected water storage before material transport')
                     ! A real protected compartment keeps its own signature even
                     ! below the dry threshold. Falling back to visible C here
                     ! can make protected material demand exceed protected mass
                     ! and creates a spurious sub-unit BIF donor dependency.
                     IF (levsto(i) > 0._r8) THEN
                        ! Use the true pool volume even below the diagnostic dry
                        ! threshold. The water limiter guarantees Q*dt<=V, so a
                        ! complete tiny-pool drain exports exactly its material mass
                        ! instead of leaving an orphan inventory behind.
                        mt_prot_conc_flux(i) = mt_levsto(imt, i) / levsto(i)
                        IF (.not. ieee_is_finite(mt_prot_conc_flux(i))) &
                           CALL CoLM_stop('non-finite protected material concentration')
                     ELSE
                        mt_prot_conc_flux(i) = 0._r8
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

         ! --- 3. Get downstream concentration (main channel) ---
         CALL worker_push_data (push_next2ucat, mt_conc_flux, conc_next, fillvalue = R_fill)

         ! --- 4. Upwind main-channel material flux ---
         DO i = 1, numucat
            IF (.not. ucatfilter(i)) THEN
               mt_flux(i) = 0._r8
               CYCLE
            ENDIF
            IF (hflux_fc(i) >= 0._r8) THEN
               mt_flux(i) = mt_conc_flux(i) * hflux_fc(i)
            ELSE
               mt_flux(i) = conc_next(i) * hflux_fc(i)
            ENDIF
         ENDDO

         ! --- 5. Aggregate upstream material flux ---
         CALL worker_push_data (push_ups2ucat, mt_flux, flux_ups, fillvalue = 0._r8, mode = 'sum')

         ! --- 6. Bifurcation pathway material flux ---
         bif_net(:) = 0._r8
         mt_bif_lev_net(:) = 0._r8
         IF (bif_workspace_active) THEN
            CALL worker_push_data (push_bif_dn2pth, mt_conc_flux, conc_dn_pth, fillvalue = 0._r8)
            CALL worker_push_data (push_bif_dn2pth, mt_prot_conc_flux, mt_prot_conc_dn_pth, fillvalue = 0._r8)

            mt_pth_single(:) = 0._r8
            mt_pth_lev(:) = 0._r8
            mt_pth_levee(:,:) = 0._r8
            DO ipth = 1, npth_bif
               i_up = pth_upst_local(ipth)
               IF (i_up < 1 .or. i_up > numucat) CYCLE
               IF (.not. ucatfilter(i_up)) CYCLE

               i_dn = pth_down_local(ipth)
               CALL material_bif_path_levee_sides(can_use_levee_material, i_up, i_dn, ipth, numucat, &
                  upstream_has_levee, downstream_has_levee)

               DO ilev = 1, nlev_bif
                  layer_wflux = bif_hflux_lev_in(ilev, ipth)
                  IF (abs(layer_wflux) <= mt_tiny) CYCLE

                  IF (layer_wflux >= 0._r8) THEN
                     IF (ilev > 1 .and. upstream_has_levee) THEN
                        mt_pth_fl = mt_prot_conc_flux(i_up) * layer_wflux
                     ELSE
                        mt_pth_fl = mt_conc_flux(i_up) * layer_wflux
                     ENDIF
                  ELSE
                     IF (ilev > 1 .and. downstream_has_levee) THEN
                        mt_pth_fl = mt_prot_conc_dn_pth(ipth) * layer_wflux
                     ELSE
                        mt_pth_fl = conc_dn_pth(ipth) * layer_wflux
                     ENDIF
                  ENDIF
                  mt_pth_levee(ilev, ipth) = mt_pth_fl

                  IF (ilev > 1 .and. upstream_has_levee) THEN
                     mt_bif_lev_net(i_up) = mt_bif_lev_net(i_up) + mt_pth_fl
                  ELSE
                     bif_net(i_up) = bif_net(i_up) + mt_pth_fl
                  ENDIF

                  ! Every downstream cell, local or remote, receives its share
                  ! through push_bif_influx below, summed in global pathway
                  ! order as the water flux is (MOD_Grid_RiverLakeBifurcation
                  ! step 6). Adding local downstream cells here directly made
                  ! the result depend on the number of processes.
                  IF (ilev > 1 .and. downstream_has_levee) THEN
                     mt_pth_lev(ipth) = mt_pth_lev(ipth) + mt_pth_fl
                  ELSE
                     mt_pth_single(ipth) = mt_pth_single(ipth) + mt_pth_fl
                  ENDIF
               ENDDO
            ENDDO

            ! Bif scatter to all downstream cells.  mt_pth_* keep the
            ! upstream-sign pathway flux; the downstream cell's net
            ! contribution is therefore subtracted after push.
            CALL worker_push_data (push_bif_influx, mt_pth_single, bif_recv, &
               fillvalue = 0._r8, mode = 'sum')
            CALL worker_push_data (push_bif_influx, mt_pth_lev, bif_lev_recv, &
               fillvalue = 0._r8, mode = 'sum')
            DO i = 1, numucat
               bif_net(i) = bif_net(i) - bif_recv(i)
               mt_bif_lev_net(i) = mt_bif_lev_net(i) - bif_lev_recv(i)
            ENDDO
         ENDIF

         ! --- 7. Coupled mass-based flux limiter ---
         ! Water CFL limits NET storage change, while a material donor can have
         ! several simultaneous outgoing faces.  Build the raw gross donor
         ! demand first; the coupled iteration below authorises it only from
         ! mass that upstream donor rates actually deliver in this substep.

         ! 7a: P2STOOUT — total outgoing |material| per cell (water-direction)
         ! Guard dt_all(irivsys(i)) the same way as section 0/8 so a
         ! corrupt/uninitialised irivsys doesn't segfault here even if the
         ! water-flow main loop happens to skip it via ucatfilter. See
         ! MOD_Grid_RiverLakeFlow.F90:585 for the primary defense.
         DO i = 1, numucat
            IF (hflux_fc(i) >= 0._r8 .and. &
                irivsys(i) > 0 .and. irivsys(i) <= size(dt_all)) THEN
               mt_out_mass(i) = abs(mt_flux(i)) * dt_all(irivsys(i))
            ELSE
               mt_out_mass(i) = 0._r8
            ENDIF
         ENDDO
         ! Reverse flow: downstream cell is donor → push to its P2STOOUT.
         ! A normal river edge cannot cross river systems: Network assigns an
         ! upstream cell the same rivermouth label as ucat_next, and Flow uses
         ! one dt_all entry per label (reduced on p_comm_rivsys when split).
         ! Therefore this edge/receiver dt is also the downstream donor dt.
         DO i = 1, numucat
            IF (hflux_fc(i) < 0._r8 .and. &
                irivsys(i) > 0 .and. irivsys(i) <= size(dt_all)) THEN
               rate_cell(i) = abs(mt_flux(i)) * dt_all(irivsys(i))
            ELSE
               rate_cell(i) = 0._r8
            ENDIF
         ENDDO
         CALL worker_push_data (push_ups2ucat, rate_cell, rate_next, fillvalue = 0._r8, mode = 'sum')
         DO i = 1, numucat
            mt_out_mass(i) = mt_out_mass(i) + rate_next(i)
            mt_out_mass_lev(i) = 0._r8
         ENDDO
         ! Bif pathways: add sender-side contribution
         IF (bif_workspace_active) THEN
            mt_dn_out_vis_pth(:) = 0._r8
            mt_dn_out_lev_pth(:) = 0._r8
            DO ipth = 1, npth_bif
               i_up = pth_upst_local(ipth)
               IF (i_up < 1 .or. i_up > numucat) CYCLE
               IF (.not. (irivsys(i_up) > 0 .and. irivsys(i_up) <= size(dt_all))) CYCLE
               dt_i = dt_all(irivsys(i_up))
               i_dn = pth_down_local(ipth)
               CALL material_bif_path_levee_sides(can_use_levee_material, i_up, i_dn, ipth, numucat, &
                  upstream_has_levee, downstream_has_levee)

               DO ilev = 1, nlev_bif
                  layer_wflux = bif_hflux_lev_in(ilev, ipth)
                  IF (abs(layer_wflux) <= mt_tiny) CYCLE
                  mt_pth_fl = mt_pth_levee(ilev, ipth)
                  IF (layer_wflux >= 0._r8) THEN
                     IF (ilev > 1 .and. upstream_has_levee) THEN
                        mt_out_mass_lev(i_up) = mt_out_mass_lev(i_up) + abs(mt_pth_fl) * dt_i
                     ELSE
                        mt_out_mass(i_up) = mt_out_mass(i_up) + abs(mt_pth_fl) * dt_i
                     ENDIF
                  ELSE
                     dt_donor = 0._r8
                     IF (i_dn > 0 .and. i_dn <= numucat) THEN
                        IF (irivsys(i_dn) > 0 .and. irivsys(i_dn) <= size(dt_all)) &
                           dt_donor = dt_all(irivsys(i_dn))
                     ELSEIF (bif_workspace_active) THEN
                        dt_donor = dt_dn_pth(ipth)
                     ENDIF
                     IF (dt_donor <= 0._r8) CYCLE
                     IF (ilev > 1 .and. downstream_has_levee) THEN
                        mt_dn_out_lev_pth(ipth) = mt_dn_out_lev_pth(ipth) + abs(mt_pth_fl) * dt_donor
                     ELSE
                        mt_dn_out_vis_pth(ipth) = mt_dn_out_vis_pth(ipth) + abs(mt_pth_fl) * dt_donor
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
            CALL worker_push_data (push_bif_influx, mt_dn_out_vis_pth, mt_dn_out_vis_recv, &
               fillvalue = 0._r8, mode = 'sum')
            CALL worker_push_data (push_bif_influx, mt_dn_out_lev_pth, mt_dn_out_lev_recv, &
               fillvalue = 0._r8, mode = 'sum')
            DO i = 1, numucat
               mt_out_mass(i) = mt_out_mass(i) + mt_dn_out_vis_recv(i)
               mt_out_mass_lev(i) = mt_out_mass_lev(i) + mt_dn_out_lev_recv(i)
            ENDDO
         ENDIF

         ! 7b: D2RATE fixed point per cell
         ! Start at T(0): only pre-substep donor mass is authorised.  This is
         ! the first iterate of a monotone-from-zero fixed point, not a guess
         ! based on unscaled same-step inflow.  Given a conservative old rate,
         ! T(old) is also conservative because its larger rate can only increase
         ! every receiver's incoming mass.  Therefore every iterate is safe;
         ! convergence recovers same-step pass-through without ever permitting
         ! a downstream cell to spend inflow that its upstream donor later cuts.
         DO i = 1, numucat
            IF (.not. ieee_is_finite(mt_mass(imt, i)) .or. &
                .not. ieee_is_finite(mt_out_mass(i)) .or. &
                .not. ieee_is_finite(mt_out_mass_lev(i))) &
               CALL CoLM_stop('non-finite river material donor limiter state')
            IF (mt_mass(imt, i) < -MT_RESTART_NEGATIVE_DUST) &
               CALL CoLM_stop('negative river material mass entered donor limiter')
            IF (mt_out_mass(i) < 0._r8 .or. mt_out_mass_lev(i) < 0._r8) &
               CALL CoLM_stop('negative river material gross outflow demand')

            IF (mt_out_mass(i) > MT_LIMITER_OUT_TINY) THEN
               rate_cell(i) = min(max(mt_mass(imt, i), 0._r8) / mt_out_mass(i), 1._r8)
            ELSE
               rate_cell(i) = 1._r8
            ENDIF
            IF (mt_out_mass_lev(i) > MT_LIMITER_OUT_TINY) THEN
               IF (.not. allocated(mt_levsto)) THEN
                  CALL CoLM_stop('protected material outflow has no protected donor pool')
               ELSEIF (i > size(mt_levsto, 2)) THEN
                  CALL CoLM_stop('protected material outflow has no protected donor pool')
               ELSE
                  IF (.not. ieee_is_finite(mt_levsto(imt, i))) &
                     CALL CoLM_stop('non-finite protected material donor limiter state')
                  IF (mt_levsto(imt, i) < -MT_RESTART_NEGATIVE_DUST) &
                     CALL CoLM_stop('negative protected material mass entered donor limiter')
                  rate_cell_lev(i) = min(max(mt_levsto(imt, i), 0._r8) / mt_out_mass_lev(i), 1._r8)
               ENDIF
            ELSE
               rate_cell_lev(i) = 1._r8
            ENDIF
         ENDDO

         ! Skip communication entirely when T(0)==1 everywhere.  BIF mode is
         ! globally lock-stepped and includes empty workers.  Without BIF,
         ! independent river systems can take different numbers of substeps;
         ! only a river system actually split across ranks may reduce here.
         limiter_delta = 0._r8
         DO i = 1, numucat
            limiter_delta = max(limiter_delta, 1._r8 - rate_cell(i), &
               1._r8 - rate_cell_lev(i))
         ENDDO
         limiter_iters_used = 0
         limiter_delta_global = limiter_delta
#ifdef USEMPI
         IF (bif_workspace_active) THEN
            CALL mpi_allreduce(limiter_delta, limiter_delta_global, 1, MPI_REAL8, MPI_MAX, &
               p_comm_worker, p_err)
         ELSEIF (rivsys_by_multiple_procs) THEN
            CALL mpi_allreduce(limiter_delta, limiter_delta_global, 1, MPI_REAL8, MPI_MAX, &
               p_comm_rivsys, p_err)
         ENDIF
#endif
         IF (.not. ieee_is_finite(limiter_delta_global)) &
            CALL CoLM_stop('non-finite river material donor limiter residual')
         limiter_converged = limiter_delta_global <= MT_LIMITER_RATE_TOL

         DO limiter_iter = 1, limiter_max_iter
            IF (limiter_converged) EXIT
            limiter_iters_used = limiter_iter
            ! Fast-path diagnostic only: this call needs more than the soft
            ! threshold of rounds (long chain); iteration continues regardless.
            IF (limiter_iter == MT_LIMITER_SOFT_ITER) &
               mt_limiter_over_soft = mt_limiter_over_soft + 1

            ! Main-channel actual incoming mass under the current sender rates.
            CALL worker_push_data(push_next2ucat, rate_cell, rate_next, fillvalue = 1._r8)
            mt_in_mass(:) = 0._r8
            mt_in_mass_lev(:) = 0._r8
            mt_inp_step(:) = 0._r8
            DO i = 1, numucat
               dt_i = 0._r8
               IF (irivsys(i) > 0 .and. irivsys(i) <= size(dt_all)) dt_i = dt_all(irivsys(i))
               IF (dt_i <= 0._r8) CYCLE
               IF (hflux_fc(i) >= 0._r8) THEN
                  mt_inp_step(i) = max(mt_flux(i) * rate_cell(i), 0._r8)
               ELSE
                  mt_in_mass(i) = mt_in_mass(i) + &
                     max(-mt_flux(i) * rate_next(i), 0._r8) * dt_i
               ENDIF
            ENDDO
            CALL worker_push_data(push_ups2ucat, mt_inp_step, flux_ups, &
               fillvalue = 0._r8, mode = 'sum')
            DO i = 1, numucat
               dt_i = 0._r8
               IF (irivsys(i) > 0 .and. irivsys(i) <= size(dt_all)) dt_i = dt_all(irivsys(i))
               IF (dt_i > 0._r8) mt_in_mass(i) = mt_in_mass(i) + max(flux_ups(i), 0._r8) * dt_i
            ENDDO

            ! Bifurcation actual incoming mass.  Rebuild gross receiver-side
            ! amounts path by path: signed net BIF flux can hide simultaneous
            ! incoming and outgoing transfers at the same cell.
            IF (bif_workspace_active) THEN
               CALL worker_push_data(push_bif_dn2pth, rate_cell, rate_dn_pth, fillvalue = 1._r8)
               CALL worker_push_data(push_bif_dn2pth, rate_cell_lev, rate_dn_pth_lev, fillvalue = 1._r8)
               mt_pth_single(:) = 0._r8
               mt_pth_lev(:) = 0._r8
               DO ipth = 1, npth_bif
                  i_up = pth_upst_local(ipth)
                  IF (i_up < 1 .or. i_up > numucat) CYCLE
                  IF (.not. (irivsys(i_up) > 0 .and. irivsys(i_up) <= size(dt_all))) CYCLE
                  dt_i = dt_all(irivsys(i_up))
                  IF (dt_i <= 0._r8) CYCLE

                  i_dn = pth_down_local(ipth)
                  CALL material_bif_path_levee_sides(can_use_levee_material, i_up, i_dn, ipth, numucat, &
                     upstream_has_levee, downstream_has_levee)

                  DO ilev = 1, nlev_bif
                     layer_wflux = bif_hflux_lev_in(ilev, ipth)
                     mt_pth_fl = mt_pth_levee(ilev, ipth)
                     IF (abs(layer_wflux) <= mt_tiny .or. abs(mt_pth_fl) <= mt_tiny) CYCLE

                     IF (layer_wflux >= 0._r8) THEN
                        IF (ilev > 1 .and. upstream_has_levee) THEN
                           mt_rate = rate_cell_lev(i_up)
                        ELSE
                           mt_rate = rate_cell(i_up)
                        ENDIF
                        mt_pth_fl = mt_pth_fl * mt_rate
                        ! upstream -> downstream: receiver is i_dn, local or
                        ! remote; both go through push_bif_influx (global
                        ! pathway order, see step 6). dt_dn_pth is dt_all of
                        ! the receiver, pushed from dt_ucat.
                        dt_donor = dt_dn_pth(ipth)
                        IF (dt_donor <= 0._r8) CYCLE
                        IF (ilev > 1 .and. downstream_has_levee) THEN
                           mt_pth_lev(ipth) = mt_pth_lev(ipth) + max(mt_pth_fl, 0._r8) * dt_donor
                        ELSE
                           mt_pth_single(ipth) = mt_pth_single(ipth) + max(mt_pth_fl, 0._r8) * dt_donor
                        ENDIF
                     ELSE
                        IF (ilev > 1 .and. downstream_has_levee) THEN
                           mt_rate = rate_dn_pth_lev(ipth)
                        ELSE
                           mt_rate = rate_dn_pth(ipth)
                        ENDIF
                        mt_pth_fl = mt_pth_fl * mt_rate
                        ! downstream -> upstream: receiver is local i_up.
                        IF (ilev > 1 .and. upstream_has_levee) THEN
                           mt_in_mass_lev(i_up) = mt_in_mass_lev(i_up) + max(-mt_pth_fl, 0._r8) * dt_i
                        ELSE
                           mt_in_mass(i_up) = mt_in_mass(i_up) + max(-mt_pth_fl, 0._r8) * dt_i
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
               CALL worker_push_data(push_bif_influx, mt_pth_single, bif_recv, &
                  fillvalue = 0._r8, mode = 'sum')
               CALL worker_push_data(push_bif_influx, mt_pth_lev, bif_lev_recv, &
                  fillvalue = 0._r8, mode = 'sum')
               DO i = 1, numucat
                  mt_in_mass(i) = mt_in_mass(i) + max(bif_recv(i), 0._r8)
                  mt_in_mass_lev(i) = mt_in_mass_lev(i) + max(bif_lev_recv(i), 0._r8)
               ENDDO
            ENDIF

            limiter_delta = 0._r8
            DO i = 1, numucat
               IF (.not. ieee_is_finite(mt_in_mass(i)) .or. &
                   .not. ieee_is_finite(mt_in_mass_lev(i)) .or. &
                   mt_in_mass(i) < 0._r8 .or. mt_in_mass_lev(i) < 0._r8) &
                  CALL CoLM_stop('invalid actual incoming mass in river material donor limiter')

               IF (mt_out_mass(i) > MT_LIMITER_OUT_TINY) THEN
                  limiter_rate_new = min((max(mt_mass(imt, i), 0._r8) + mt_in_mass(i)) / &
                     mt_out_mass(i), 1._r8)
                  IF (.not. ieee_is_finite(limiter_rate_new)) &
                     CALL CoLM_stop('non-finite visible river material donor rate')
                  IF (limiter_rate_new + MT_LIMITER_RATE_TOL < rate_cell(i)) &
                     CALL CoLM_stop('river material donor limiter lost monotonicity')
                  limiter_rate_new = max(rate_cell(i), limiter_rate_new)
                  limiter_delta = max(limiter_delta, limiter_rate_new - rate_cell(i))
                  rate_cell(i) = limiter_rate_new
               ENDIF

               IF (mt_out_mass_lev(i) > MT_LIMITER_OUT_TINY) THEN
                  limiter_rate_new = min((max(mt_levsto(imt, i), 0._r8) + mt_in_mass_lev(i)) / &
                     mt_out_mass_lev(i), 1._r8)
                  IF (.not. ieee_is_finite(limiter_rate_new)) &
                     CALL CoLM_stop('non-finite protected river material donor rate')
                  IF (limiter_rate_new + MT_LIMITER_RATE_TOL < rate_cell_lev(i)) &
                     CALL CoLM_stop('protected material donor limiter lost monotonicity')
                  limiter_rate_new = max(rate_cell_lev(i), limiter_rate_new)
                  limiter_delta = max(limiter_delta, limiter_rate_new - rate_cell_lev(i))
                  rate_cell_lev(i) = limiter_rate_new
               ENDIF
            ENDDO

            limiter_delta_global = limiter_delta
#ifdef USEMPI
            IF (bif_workspace_active) THEN
               CALL mpi_allreduce(limiter_delta, limiter_delta_global, 1, MPI_REAL8, MPI_MAX, &
                  p_comm_worker, p_err)
            ELSEIF (rivsys_by_multiple_procs) THEN
               CALL mpi_allreduce(limiter_delta, limiter_delta_global, 1, MPI_REAL8, MPI_MAX, &
                  p_comm_rivsys, p_err)
            ENDIF
#endif
            IF (.not. ieee_is_finite(limiter_delta_global)) &
               CALL CoLM_stop('non-finite river material donor limiter residual')
            limiter_converged = limiter_delta_global <= MT_LIMITER_RATE_TOL
         ENDDO

         mt_limiter_calls = mt_limiter_calls + 1
         mt_limiter_iter_sum = mt_limiter_iter_sum + limiter_iters_used
         IF (limiter_iters_used > mt_limiter_iter_peak) &
            mt_limiter_iter_peak = limiter_iters_used

         IF (.not. limiter_converged) THEN
            ! The 2*N+1 topology bound guarantees the monotone fixed point is
            ! reachable, so exhausting it means the limiter itself is broken
            ! (a truncated exit would silently under-export into the dry-cell
            ! sink).  Keep this fatal.
            WRITE(*,'(A,I0,A,I0,A,ES12.4)') 'ERROR material donor limiter: material=', imt, &
               ' iterations=', limiter_iters_used, ' residual=', limiter_delta_global
            CALL CoLM_stop('river material donor limiter did not converge')
         ENDIF

         ! 7c: Get downstream cell rate (for reverse main flow)
         CALL worker_push_data (push_next2ucat, rate_cell, rate_next, fillvalue = 1._r8)

         ! 7d: Scale main-channel flux by sender rate (water direction)
         DO i = 1, numucat
            IF (hflux_fc(i) >= 0._r8) THEN
               mt_flux(i) = mt_flux(i) * rate_cell(i)
            ELSE
               mt_flux(i) = mt_flux(i) * rate_next(i)
            ENDIF
         ENDDO

         ! 7e: Scale bif pathway flux by sender rate (water direction)
         IF (bif_workspace_active) THEN
            CALL worker_push_data (push_bif_dn2pth, rate_cell, rate_dn_pth, fillvalue = 1._r8)
            CALL worker_push_data (push_bif_dn2pth, rate_cell_lev, rate_dn_pth_lev, fillvalue = 1._r8)
            DO ipth = 1, npth_bif
               i_up = pth_upst_local(ipth)
               IF (i_up < 1 .or. i_up > numucat) CYCLE
               i_dn = pth_down_local(ipth)
               CALL material_bif_path_levee_sides(can_use_levee_material, i_up, i_dn, ipth, numucat, &
                  upstream_has_levee, downstream_has_levee)
               DO ilev = 1, nlev_bif
                  layer_wflux = bif_hflux_lev_in(ilev, ipth)
                  IF (abs(layer_wflux) <= mt_tiny) CYCLE
                  IF (layer_wflux >= 0._r8) THEN
                     IF (ilev > 1 .and. upstream_has_levee) THEN
                        mt_rate = rate_cell_lev(i_up)
                     ELSE
                        mt_rate = rate_cell(i_up)
                     ENDIF
                  ELSE
                     IF (ilev > 1 .and. downstream_has_levee) THEN
                        mt_rate = rate_dn_pth_lev(ipth)
                     ELSE
                        mt_rate = rate_dn_pth(ipth)
                     ENDIF
                  ENDIF
                  mt_pth_levee(ilev, ipth) = mt_pth_levee(ilev, ipth) * mt_rate
               ENDDO
            ENDDO
            ! Rebuild visible/protected bif nets from scaled layer fluxes.
            bif_net(:) = 0._r8
            mt_bif_lev_net(:) = 0._r8
            mt_pth_single(:) = 0._r8
            mt_pth_lev(:) = 0._r8
            DO ipth = 1, npth_bif
               i_up = pth_upst_local(ipth)
               IF (i_up < 1 .or. i_up > numucat) CYCLE
               i_dn = pth_down_local(ipth)
               CALL material_bif_path_levee_sides(can_use_levee_material, i_up, i_dn, ipth, numucat, &
                  upstream_has_levee, downstream_has_levee)

               DO ilev = 1, nlev_bif
                  mt_pth_fl = mt_pth_levee(ilev, ipth)
                  IF (abs(mt_pth_fl) <= mt_tiny) CYCLE
                  IF (ilev > 1 .and. upstream_has_levee) THEN
                     mt_bif_lev_net(i_up) = mt_bif_lev_net(i_up) + mt_pth_fl
                  ELSE
                     bif_net(i_up) = bif_net(i_up) + mt_pth_fl
                  ENDIF

                  ! all downstream cells through the push (see step 6)
                  IF (ilev > 1 .and. downstream_has_levee) THEN
                     mt_pth_lev(ipth) = mt_pth_lev(ipth) + mt_pth_fl
                  ELSE
                     mt_pth_single(ipth) = mt_pth_single(ipth) + mt_pth_fl
                  ENDIF
               ENDDO
            ENDDO
            CALL worker_push_data (push_bif_influx, mt_pth_single, bif_recv, &
               fillvalue = 0._r8, mode = 'sum')
            CALL worker_push_data (push_bif_influx, mt_pth_lev, bif_lev_recv, &
               fillvalue = 0._r8, mode = 'sum')
            DO i = 1, numucat
               bif_net(i) = bif_net(i) - bif_recv(i)
               mt_bif_lev_net(i) = mt_bif_lev_net(i) - bif_lev_recv(i)
            ENDDO
         ENDIF

         ! Re-aggregate upstream material flux after limiting
         CALL worker_push_data (push_ups2ucat, mt_flux, flux_ups, fillvalue = 0._r8, mode = 'sum')

         ! --- 8. Update material mass ---
         DO i = 1, numucat
            IF (.not. ucatfilter(i)) CYCLE
            dt_i = 0._r8
            IF (irivsys(i) > 0 .and. irivsys(i) <= size(dt_all)) dt_i = dt_all(irivsys(i))
            IF (dt_i <= 0._r8) CYCLE

            mt_mass_new = mt_mass(imt, i) &
               + (- mt_flux(i) + flux_ups(i) - bif_net(i)) * dt_i
            IF (.not. ieee_is_finite(mt_mass_new)) &
               CALL CoLM_stop('non-finite river material mass after coupled donor limiter')
            IF (mt_mass_new < -max(MT_RESTART_NEGATIVE_DUST, MT_UPDATE_ROUNDOFF * &
                (abs(mt_mass(imt, i)) + (abs(mt_flux(i)) + abs(flux_ups(i)) + abs(bif_net(i))) * dt_i))) &
               CALL CoLM_stop('negative river material mass after coupled donor limiter')

            ! Transport owns only conservative edge transfers.  Never snap the
            ! result back to R_fill: doing so is an unbooked source/sink and can
            ! erase a legitimate restart or upstream isotope gradient.
            mt_mass(imt, i) = max(mt_mass_new, 0._r8)
            IF (allocated(mt_levsto)) THEN
               IF (i <= size(mt_levsto, 2)) THEN
                  levsto_old = mt_levsto(imt, i)
                  mt_levsto(imt, i) = mt_levsto(imt, i) - mt_bif_lev_net(i) * dt_i
                  IF (.not. ieee_is_finite(mt_levsto(imt, i))) &
                     CALL CoLM_stop('non-finite protected material mass after coupled donor limiter')
                  IF (mt_levsto(imt, i) < -max(MT_RESTART_NEGATIVE_DUST, MT_UPDATE_ROUNDOFF * &
                      (abs(levsto_old) + abs(mt_bif_lev_net(i)) * dt_i))) &
                     CALL CoLM_stop('negative protected material mass after coupled donor limiter')
                  mt_levsto(imt, i) = max(mt_levsto(imt, i), 0._r8)
               ENDIF
            ENDIF

            decay_fraction = material_decay_fraction(imt, dt_i)
            IF (decay_fraction > 0._r8) THEN
               decay_src = 0._r8
               CALL decay_river_pool(mt_mass(imt, i), decay_fraction, decay_src)
               CALL decay_river_pool(mt_inp_buf(imt, i), decay_fraction, decay_src)
               IF (allocated(mt_solid)) THEN
                  IF (mt_has_dissolved_limit(imt)) THEN
                     CALL decay_river_pool(mt_solid(imt, i), decay_fraction, decay_src)
                     CALL decay_river_pool(mt_levsto_solid(imt, i), decay_fraction, decay_src)
                  ENDIF
               ENDIF
               IF (allocated(mt_levsto)) THEN
                  IF (i <= size(mt_levsto, 2)) THEN
                     CALL decay_river_pool(mt_levsto(imt, i), decay_fraction, decay_src)
                  ENDIF
               ENDIF
#ifdef CoLMDEBUG
               ! Only the CoLMDEBUG mass-balance check in MOD_Grid_RiverLakeFlow
               ! resets and reads this accumulator.
               IF (allocated(mt_decay_source)) THEN
                  mt_decay_source(imt, i) = mt_decay_source(imt, i) + decay_src
               ENDIF
#endif
            ENDIF
         ENDDO

         ! --- 9. Save flux for diagnostics (only for active cells) ---
         DO i = 1, numucat
            IF (ucatfilter(i)) THEN
               mt_flux_out(imt, i) = mt_flux(i)
               mt_bif_net_saved(imt, i) = bif_net(i) + mt_bif_lev_net(i)
            ENDIF
         ENDDO

      ENDDO  ! imt

   END SUBROUTINE mt_riverlake_substep

   ! ------------------------------------------------------------------
   ! Donor-limiter work statistics (per rank).  calls = limiter invocations
   ! (one per water-transported material per substep), iter_sum/iter_peak = Jacobi
   ! iterations actually spent, over_soft = invocations that needed more than
   ! MT_LIMITER_SOFT_ITER rounds (long-chain indicator; no truncation happens).
   ! Pass reset=.true. to read and clear in one call; callers that print a
   ! global summary must reduce (SUM for calls/iter_sum/over_soft, MAX for peak).
   SUBROUTINE mt_riverlake_limiter_stats (calls, iter_sum, iter_peak, over_soft, reset)
   IMPLICIT NONE
   integer, intent(out) :: calls, iter_sum, iter_peak, over_soft
   logical, intent(in), optional :: reset

      calls     = mt_limiter_calls
      iter_sum  = mt_limiter_iter_sum
      iter_peak = mt_limiter_iter_peak
      over_soft = mt_limiter_over_soft
      IF (present(reset)) THEN
         IF (reset) THEN
            mt_limiter_calls     = 0
            mt_limiter_iter_sum  = 0
            mt_limiter_iter_peak = 0
            mt_limiter_over_soft = 0
         ENDIF
      ENDIF
   END SUBROUTINE mt_riverlake_limiter_stats

   ! Move any direction-pending dry-cleanup mass back into the cell inventory
   ! before a restart/window boundary, so a pending buffer can never be lost
   ! across restart.  Released mass is available to the next transport step.
   SUBROUTINE mt_riverlake_release_pending_drain ()
      USE MOD_Grid_RiverLakeNetwork, only: numucat
      USE MOD_MT_Defs, only: mt_tiny
      integer :: i, imt
      IF (.not. allocated(mt_dry_pending)) RETURN
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         DO i = 1, numucat
            IF (mt_dry_pending(imt, i) > mt_tiny) THEN
               mt_mass(imt, i) = mt_mass(imt, i) + mt_dry_pending(imt, i)
               mt_dry_pending(imt, i) = 0._r8
            ENDIF
         ENDDO
      ENDDO
   END SUBROUTINE mt_riverlake_release_pending_drain

   SUBROUTINE decay_river_pool (pool, fraction, source_sink)
      USE MOD_MT_Defs, only: mt_tiny
      IMPLICIT NONE
      real(r8), intent(inout) :: pool
      real(r8), intent(in)    :: fraction
      real(r8), intent(inout) :: source_sink
      real(r8) :: before

      IF (pool <= mt_tiny) RETURN
      before = pool
      pool = pool * (1._r8 - fraction)
      source_sink = source_sink + pool - before
   END SUBROUTINE decay_river_pool

   logical FUNCTION riverlake_restart_descriptor_compatible (file_restart)
      USE MOD_NetCDFSerial, only: ncio_var_exist, ncio_inquire_varsize, ncio_read_serial
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart

      integer :: flags(4), status, schema, descriptor_count, identity_count, complete
      integer, allocatable :: expected(:,:), ondisk(:,:), varsize(:)

      ! status: 1=current and matching, 0=legacy/incompatible (cold start),
      !        -1=malformed/incomplete metadata (fatal). A complete descriptor
      ! from another schema or material configuration is safe but incompatible;
      ! inconsistent ranks, fixed width, or internal counts indicate damage.
      CALL material_build_descriptor_identity(expected, transport_only=.true.)
      river_restart_schema_loaded = 0
      schema = 0
      flags = 0
      status = 0
      IF (p_is_master) THEN
         IF (ncio_var_exist(file_restart, 'mt_river_restart_schema', readflag=.false.)) flags(1) = 1
         IF (ncio_var_exist(file_restart, 'mt_river_descriptor_count', readflag=.false.)) flags(2) = 1
         IF (ncio_var_exist(file_restart, 'mt_river_descriptor_identity', readflag=.false.)) flags(3) = 1
         IF (ncio_var_exist(file_restart, 'mt_river_restart_complete', readflag=.false.)) flags(4) = 1
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast(flags, 4, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif

      ! No commit marker means a pre-transaction/legacy schema. Even if an
      ! earlier build wrote descriptor-like fields, their atomicity cannot be
      ! established, so never scatter their state.
      IF (flags(4) == 0) THEN
         IF (p_is_master) WRITE(*,'(A)') &
            '  NOTE mt_riverlake_read_restart: restart has no material commit marker; cold-starting river materials.'
         riverlake_restart_descriptor_compatible = .false.
         deallocate(expected)
         RETURN
      ENDIF

      IF (flags(1) /= 1 .or. flags(2) /= 1) THEN
         status = -1
      ELSEIF (p_is_master) THEN
         status = 1
         identity_count = -1

         CALL ncio_inquire_varsize(file_restart, 'mt_river_restart_schema', varsize)
         IF (.not. allocated(varsize)) THEN
            status = -1
         ELSE
            IF (size(varsize) /= 0) status = -1
            deallocate(varsize)
         ENDIF
         CALL ncio_inquire_varsize(file_restart, 'mt_river_descriptor_count', varsize)
         IF (.not. allocated(varsize)) THEN
            status = -1
         ELSE
            IF (size(varsize) /= 0) status = -1
            deallocate(varsize)
         ENDIF
         CALL ncio_inquire_varsize(file_restart, 'mt_river_restart_complete', varsize)
         IF (.not. allocated(varsize)) THEN
            status = -1
         ELSE
            IF (size(varsize) /= 0) status = -1
            deallocate(varsize)
         ENDIF

         IF (status == 1) THEN
            CALL ncio_read_serial(file_restart, 'mt_river_restart_complete', complete)
            IF (complete /= 1) status = -1
         ENDIF
         IF (status == 1) THEN
            CALL ncio_read_serial(file_restart, 'mt_river_restart_schema', schema)
            CALL ncio_read_serial(file_restart, 'mt_river_descriptor_count', descriptor_count)
            IF (descriptor_count < 0 .or. descriptor_count > 1000) status = -1
         ENDIF
         IF (status == 1) THEN
            IF (schema /= 1 .and. schema /= RIVER_MT_RESTART_SCHEMA_VERSION) THEN
               status = 0
            ELSEIF (descriptor_count == 0) THEN
               ! Coherent empty transaction: provider-only configurations have
               ! no generic river state and therefore no zero-length identity
               ! variable. A stale identity variable in a reused NetCDF target
               ! is ignored because committed count=0 is authoritative.
               status = merge(1, 0, size(expected, 2) == 0)
            ELSEIF (flags(3) /= 1) THEN
               status = -1
            ELSE
               CALL ncio_inquire_varsize(file_restart, 'mt_river_descriptor_identity', varsize)
               IF (.not. allocated(varsize)) THEN
                  status = -1
               ELSE
                  IF (size(varsize) /= 2) THEN
                     status = -1
                  ELSEIF (varsize(1) /= MATERIAL_DESCRIPTOR_IDENTITY_WIDTH) THEN
                     status = -1
                  ELSE
                     identity_count = varsize(2)
                     IF (identity_count /= descriptor_count) status = -1
                  ENDIF
                  deallocate(varsize)
               ENDIF
               IF (status == 1) THEN
                  CALL ncio_read_serial(file_restart, 'mt_river_descriptor_identity', ondisk)
                  IF (.not. allocated(ondisk)) THEN
                     status = -1
                  ELSE
                     IF (size(ondisk, 1) /= MATERIAL_DESCRIPTOR_IDENTITY_WIDTH .or. &
                         size(ondisk, 2) /= descriptor_count) THEN
                        status = -1
                     ELSEIF (descriptor_count /= size(expected, 2)) THEN
                        status = 0
                     ELSEIF (.not. all(ondisk == expected)) THEN
                        status = 0
                     ENDIF
                     deallocate(ondisk)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast(status, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(schema, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif

      IF (status < 0) THEN
         CALL CoLM_stop('malformed/incomplete river material restart descriptor metadata')
      ELSEIF (status == 0) THEN
         IF (p_is_master) WRITE(*,'(A)') &
            '  NOTE mt_riverlake_read_restart: material descriptor changed; cold-starting river materials.'
      ENDIF
      riverlake_restart_descriptor_compatible = status == 1
      IF (status == 1) river_restart_schema_loaded = schema
      deallocate(expected)
   END FUNCTION riverlake_restart_descriptor_compatible

   SUBROUTINE probe_riverlake_restart_vector (file_restart, varname, required, present_on_disk)
      USE MOD_NetCDFSerial, only: ncio_var_exist, ncio_inquire_varsize
      USE MOD_Grid_RiverLakeNetwork, only: totalnumucat
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart, varname
      logical, intent(in) :: required
      logical, intent(out) :: present_on_disk

      integer :: status
      integer, allocatable :: varsize(:)

      ! status: 1=present with the exact ucatch-vector shape, 0=absent,
      !        -1=present but malformed. Probe every required state field before
      ! the first scatter so a complete=1 file is loaded atomically or not at all.
      status = 0
      IF (p_is_master) THEN
         IF (ncio_var_exist(file_restart, trim(varname), readflag=.false.)) THEN
            status = 1
            CALL ncio_inquire_varsize(file_restart, trim(varname), varsize)
            IF (.not. allocated(varsize)) THEN
               status = -1
            ELSE
               IF (size(varsize) /= 1) THEN
                  status = -1
               ELSEIF (varsize(1) /= totalnumucat) THEN
                  status = -1
               ENDIF
               deallocate(varsize)
            ENDIF
         ENDIF
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast(status, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif

      IF (status < 0 .or. (required .and. status == 0)) THEN
         IF (p_is_master) WRITE(*,'(3A)') &
            'ERROR mt_riverlake_read_restart: required/current river material vector "', &
            trim(varname), '" is absent or has the wrong rank/ucatch extent.'
         CALL CoLM_stop('incomplete or malformed committed river material restart')
      ENDIF
      present_on_disk = status == 1
   END SUBROUTINE probe_riverlake_restart_vector

   SUBROUTINE validate_riverlake_restart_state ()
      USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
      USE MOD_Grid_RiverLakeNetwork, only: numucat
      IMPLICIT NONE

      integer :: invalid_counts(11)
      integer :: i, imt

      ! [nonfinite mass, negative mass, nonfinite signed buffer,
      !  nonfinite signed material accumulator, nonfinite runoff reference,
      !  nonfinite levee mass, negative levee mass]
      invalid_counts = 0
      IF (p_is_worker .and. numucat > 0) THEN
         DO i = 1, numucat
            IF (.not. ieee_is_finite(acc_rnof_ref(i))) invalid_counts(5) = invalid_counts(5) + 1
            DO imt = 1, nmts
               IF (.not. mt_uses_generic_water(imt)) CYCLE
               IF (.not. ieee_is_finite(mt_mass(imt, i))) THEN
                  invalid_counts(1) = invalid_counts(1) + 1
               ELSEIF (mt_mass(imt, i) < -MT_RESTART_NEGATIVE_DUST) THEN
                  invalid_counts(2) = invalid_counts(2) + 1
               ENDIF
               IF (.not. ieee_is_finite(mt_inp_buf(imt, i))) &
                  invalid_counts(3) = invalid_counts(3) + 1
               IF (.not. ieee_is_finite(acc_mt_inp(imt, i))) &
                  invalid_counts(4) = invalid_counts(4) + 1
               IF (.not. ieee_is_finite(mt_levsto(imt, i))) THEN
                  invalid_counts(6) = invalid_counts(6) + 1
               ELSEIF (mt_levsto(imt, i) < -MT_RESTART_NEGATIVE_DUST) THEN
                  invalid_counts(7) = invalid_counts(7) + 1
               ENDIF
               IF (allocated(mt_solid)) THEN
                  IF (mt_has_dissolved_limit(imt)) THEN
                     IF (.not. ieee_is_finite(mt_solid(imt, i))) THEN
                        invalid_counts(8) = invalid_counts(8) + 1
                     ELSEIF (mt_solid(imt, i) < -MT_RESTART_NEGATIVE_DUST) THEN
                        invalid_counts(9) = invalid_counts(9) + 1
                     ENDIF
                     IF (.not. ieee_is_finite(mt_levsto_solid(imt, i))) THEN
                        invalid_counts(10) = invalid_counts(10) + 1
                     ELSEIF (mt_levsto_solid(imt, i) < -MT_RESTART_NEGATIVE_DUST) THEN
                        invalid_counts(11) = invalid_counts(11) + 1
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDIF
#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, invalid_counts, size(invalid_counts), &
         MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif

      IF (any(invalid_counts > 0)) THEN
         IF (p_is_master) THEN
            WRITE(*,'(A,11(I0,1X))') 'ERROR mt_riverlake_read_restart invalid state counts '// &
               '[mass_nan mass_neg buffer_nan accinp_nan runoff_nan levee_nan levee_neg '// &
               'solid_nan solid_neg levee_solid_nan levee_solid_neg]: ', invalid_counts
         ENDIF
         CALL CoLM_stop('invalid river/lake material restart state')
      ENDIF

      ! Canonicalize only sub-picomass negative roundoff admitted by the
      ! domain check above. Signed runoff buffers/accumulators are prognostic
      ! correction state and must never be clipped.
      IF (p_is_worker .and. numucat > 0) THEN
         DO imt = 1, nmts
            IF (.not. mt_uses_generic_water(imt)) CYCLE
            WHERE (mt_mass(imt, :) < 0._r8) mt_mass(imt, :) = 0._r8
            WHERE (mt_levsto(imt, :) < 0._r8) mt_levsto(imt, :) = 0._r8
            IF (allocated(mt_solid)) THEN
               IF (mt_has_dissolved_limit(imt)) THEN
                  WHERE (mt_solid(imt, :) < 0._r8) mt_solid(imt, :) = 0._r8
                  WHERE (mt_levsto_solid(imt, :) < 0._r8) mt_levsto_solid(imt, :) = 0._r8
               ENDIF
            ENDIF
         ENDDO
      ENDIF
   END SUBROUTINE validate_riverlake_restart_state


   !-------------------------------------------------------------------------------------
   ! Flush accumulated material diagnostics
   !-------------------------------------------------------------------------------------
   SUBROUTINE mt_riverlake_flush_acc ()

   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

      ! Reset the *history* accumulators only. acc_mt_inp / acc_rnof_ref are
      ! routing-period state, not history state, and are reset at the end of
      ! each routing period in MOD_Grid_RiverLakeFlow.F90:1263-1264. Resetting
      ! acc_rnof_ref here too would leave acc_mt_inp / acc_rnof_ref with
      ! mismatched windows when history flushes happen between routing periods
      ! (acctime_rnof_max > deltim), corrupting restart recovery.
      IF (numucat > 0) THEN
         IF (allocated(a_mt_storage_mass)) a_mt_storage_mass = 0._r8
         IF (allocated(a_water_storage)) a_water_storage = 0._r8
         IF (allocated(a_mt_levsto_mass)) a_mt_levsto_mass = 0._r8
         IF (allocated(a_levsto_water)) a_levsto_water = 0._r8
         IF (allocated(a_mt_out   )) a_mt_out    = 0._r8
         IF (allocated(a_mt_bifout)) a_mt_bifout = 0._r8
         IF (allocated(a_mt_acctime)) a_mt_acctime = 0._r8
         IF (allocated(mt_n2o_source)) mt_n2o_source = 0._r8
      ENDIF

   END SUBROUTINE mt_riverlake_flush_acc

   !-------------------------------------------------------------------------------------
   ! Read material restart
   !-------------------------------------------------------------------------------------
   !-------------------------------------------------------------------------------------
   ! Drop the restored material history window and start a new one at the
   ! restart.  Called when MOD_HistAccRestart decided that the restart closes a
   ! complete history window: the rows a restart carries at a window boundary
   ! cover the window that just ended, and keeping them would make the first
   ! window after the restart span two windows.  Sharing the host's decision
   ! keeps both histories on the same window without duplicating the calendar
   ! logic here.
   !-------------------------------------------------------------------------------------
   SUBROUTINE mt_riverlake_drop_hist_window ()

   IMPLICIT NONE

      IF (.not. p_is_worker) RETURN
      IF (allocated(a_mt_storage_mass)) a_mt_storage_mass = 0._r8
      IF (allocated(a_water_storage)) a_water_storage = 0._r8
      IF (allocated(a_mt_levsto_mass)) a_mt_levsto_mass = 0._r8
      IF (allocated(a_levsto_water)) a_levsto_water = 0._r8
      IF (allocated(a_mt_out)) a_mt_out = 0._r8
      IF (allocated(a_mt_bifout)) a_mt_bifout = 0._r8
      IF (allocated(a_mt_acctime)) a_mt_acctime = 0._r8
      IF (allocated(mt_n2o_source)) mt_n2o_source = 0._r8

   END SUBROUTINE mt_riverlake_drop_hist_window


   SUBROUTINE mt_riverlake_read_restart (file_restart, found_restart, missing_mask)

   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_data_address, &
      ucat_gdid, ucat_next
   USE MOD_Grid_RiverLakeLevee, only: has_levee_bf => has_levee
   USE MOD_Grid_RiverLakeHistState, only: acctime_ucat
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart
   ! .true. when no generic transport row exists, or when the complete current
   ! transaction was loaded. Legacy/descriptor-incompatible files cold-start
   ! every generic row and report .false. only when such rows are configured.
   logical, optional, intent(out) :: found_restart
   ! Retained for caller ABI compatibility. A cold start marks every generic
   ! transport row; committed transactions never default individual rows.
   logical, optional, intent(out) :: missing_mask(:)

   integer :: imt
   integer :: ii_bf, imt_bf
   logical :: network_meta_matches, field_present, hist_complete, hist_time_present
   logical :: has_transport_material, descriptor_compatible
   character(len=64) :: varname
   real(r8), allocatable :: tmpvec(:)

      IF (present(missing_mask)) missing_mask = .false.
      IF (nmts <= 0) THEN
         IF (present(found_restart)) found_restart = .true.
         RETURN
      ENDIF
      ! Provider-owned materials have no generic river vectors, but an existing
      ! generic-river transaction must still be structurally valid.  Validate
      ! its commit metadata before the provider-only early return so complete=0
      ! or partial schema/count metadata cannot be silently accepted.  A
      ! coherent empty, legacy, or descriptor-incompatible transaction remains
      ! harmless here because no generic state will be probed or consumed.
      descriptor_compatible = riverlake_restart_descriptor_compatible(file_restart)
      has_transport_material = .false.
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         has_transport_material = .true.
      ENDDO
      IF (.not. has_transport_material) THEN
         IF (present(found_restart)) found_restart = .true.
         RETURN
      ENDIF

      ! Count/name-only metadata cannot prove that the loaded numbers retain
      ! the same units and physics. Older schemas are therefore an explicit
      ! whole-domain cold start; a complete current descriptor is required
      ! before any river material state is scattered.
      IF (.not. descriptor_compatible) THEN
         IF (present(missing_mask)) THEN
            DO imt = 1, nmts
               IF (.not. mt_uses_generic_water(imt)) CYCLE
               missing_mask(imt) = .true.
            ENDDO
         ENDIF
         IF (p_is_worker) THEN
            IF (allocated(mt_mass)) mt_mass = 0._r8
            IF (allocated(mt_conc)) mt_conc = 0._r8
            IF (allocated(mt_inp_buf)) mt_inp_buf = 0._r8
            IF (allocated(acc_mt_inp)) acc_mt_inp = 0._r8
            IF (allocated(acc_rnof_ref)) acc_rnof_ref = 0._r8
            IF (allocated(mt_levsto)) mt_levsto = 0._r8
            IF (allocated(mt_solid)) THEN
               mt_solid = 0._r8
               mt_levsto_solid = 0._r8
            ENDIF
            IF (allocated(a_mt_storage_mass)) a_mt_storage_mass = 0._r8
            IF (allocated(a_water_storage)) a_water_storage = 0._r8
            IF (allocated(a_mt_levsto_mass)) a_mt_levsto_mass = 0._r8
            IF (allocated(a_levsto_water)) a_levsto_water = 0._r8
            IF (allocated(a_mt_out)) a_mt_out = 0._r8
            IF (allocated(a_mt_bifout)) a_mt_bifout = 0._r8
            IF (allocated(a_mt_acctime)) a_mt_acctime = 0._r8
         ENDIF
         IF (present(found_restart)) found_restart = .false.
         RETURN
      ENDIF

      ! A committed current transaction must contain every required vector with
      ! the exact ucatch extent. Complete this whole preflight before the first
      ! scatter so no rank can observe a partially loaded material state.
      CALL probe_riverlake_restart_vector(file_restart, 'mt_numucat_meta', .true., field_present)
      CALL probe_riverlake_restart_vector(file_restart, 'mt_ucat_gdid_meta', .true., field_present)
      CALL probe_riverlake_restart_vector(file_restart, 'mt_ucat_next_meta', .true., field_present)
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         write(varname, '(A,A)') 'mt_mass_', trim(material_names(imt))
         CALL probe_riverlake_restart_vector(file_restart, trim(varname), .true., field_present)
         write(varname, '(A,A)') 'mt_inpbuf_', trim(material_names(imt))
         CALL probe_riverlake_restart_vector(file_restart, trim(varname), .true., field_present)
         write(varname, '(A,A)') 'mt_accinp_', trim(material_names(imt))
         CALL probe_riverlake_restart_vector(file_restart, trim(varname), .true., field_present)
         write(varname, '(A,A)') 'mt_levsto_', trim(material_names(imt))
         CALL probe_riverlake_restart_vector(file_restart, trim(varname), .true., field_present)
         IF (river_restart_schema_loaded >= 2 .and. mt_has_dissolved_limit(imt)) THEN
            write(varname, '(A,A)') 'mt_solid_', trim(material_names(imt))
            CALL probe_riverlake_restart_vector(file_restart, trim(varname), .true., field_present)
            write(varname, '(A,A)') 'mt_levsto_solid_', trim(material_names(imt))
            CALL probe_riverlake_restart_vector(file_restart, trim(varname), .true., field_present)
         ENDIF
      ENDDO
      CALL probe_riverlake_restart_vector(file_restart, 'acc_rnof_ref', .true., field_present)

      ! Schema 1 committed all mobile rows but had no finite-solute solids.
      ! Clear stale in-memory solids even when this reader is called again.
      IF (river_restart_schema_loaded == 1 .and. p_is_worker .and. allocated(mt_solid)) THEN
         mt_solid = 0._r8
         mt_levsto_solid = 0._r8
      ENDIF

      ! vector_read_and_scatter contains mpi_barrier(p_comm_glb), so
      ! ALL ranks (master, workers, IO) must enter this routine.
      IF (p_is_worker .and. numucat > 0) THEN
         allocate (tmpvec(numucat))
      ELSE
         allocate (tmpvec(0))
      ENDIF

      ! Stable grid-cell identity plus downstream global UCID catch both
      ! reordered cells and changed routing topology at equal network size.
      network_meta_matches = .true.
      CALL vector_read_and_scatter(file_restart, tmpvec, numucat, &
         'mt_ucat_gdid_meta', ucat_data_address)
      IF (p_is_worker .and. numucat > 0) &
         network_meta_matches = all(nint(tmpvec(:)) == ucat_gdid(:))
      CALL vector_read_and_scatter(file_restart, tmpvec, numucat, &
         'mt_ucat_next_meta', ucat_data_address)
      IF (p_is_worker .and. numucat > 0) &
         network_meta_matches = network_meta_matches .and. all(nint(tmpvec(:)) == ucat_next(:))
#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, network_meta_matches, 1, MPI_LOGICAL, MPI_LAND, p_comm_glb, p_err)
#endif
      IF (.not. network_meta_matches) THEN
         IF (p_is_master) WRITE(*,'(A)') &
            'ERROR: river/lake material restart belongs to a different catchment network; aborting.'
         CALL CoLM_stop()
      ENDIF

      CALL vector_read_and_scatter(file_restart, tmpvec, numucat, &
         'mt_numucat_meta', ucat_data_address)
      field_present = .true.
      IF (p_is_worker .and. numucat > 0) &
         field_present = all(nint(tmpvec(:)) == totalnumucat)
#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, field_present, 1, MPI_LOGICAL, MPI_LAND, p_comm_glb, p_err)
#endif
      IF (.not. field_present) THEN
         IF (p_is_master) WRITE(*,'(A)') &
            'ERROR: river/lake material restart catchment-count metadata is inconsistent.'
         CALL CoLM_stop()
      ENDIF

      ! The full descriptor is the single material-identity source of truth.
      ! Network size/identity are independently guarded above; duplicating
      ! count/name hashes here creates a second, weaker schema that can drift.

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         write(varname, '(A,A)') 'mt_mass_', trim(material_names(imt))
         CALL vector_read_and_scatter(file_restart, tmpvec, numucat, trim(varname), ucat_data_address)
         IF (p_is_worker .and. numucat > 0) mt_mass(imt, :) = tmpvec(:)

         write(varname, '(A,A)') 'mt_inpbuf_', trim(material_names(imt))
         CALL vector_read_and_scatter(file_restart, tmpvec, numucat, trim(varname), ucat_data_address)
         IF (p_is_worker .and. numucat > 0) mt_inp_buf(imt, :) = tmpvec(:)

         write(varname, '(A,A)') 'mt_levsto_', trim(material_names(imt))
         CALL vector_read_and_scatter(file_restart, tmpvec, numucat, trim(varname), ucat_data_address)
         IF (p_is_worker .and. numucat > 0) mt_levsto(imt, :) = tmpvec(:)
         IF (river_restart_schema_loaded >= 2 .and. mt_has_dissolved_limit(imt)) THEN
            write(varname, '(A,A)') 'mt_solid_', trim(material_names(imt))
            CALL vector_read_and_scatter(file_restart, tmpvec, numucat, trim(varname), ucat_data_address)
            IF (p_is_worker .and. numucat > 0) mt_solid(imt, :) = tmpvec(:)
            write(varname, '(A,A)') 'mt_levsto_solid_', trim(material_names(imt))
            CALL vector_read_and_scatter(file_restart, tmpvec, numucat, trim(varname), ucat_data_address)
            IF (p_is_worker .and. numucat > 0) mt_levsto_solid(imt, :) = tmpvec(:)
         ENDIF

         write(varname, '(A,A)') 'mt_accinp_', trim(material_names(imt))
         CALL vector_read_and_scatter(file_restart, tmpvec, numucat, trim(varname), ucat_data_address)
         IF (p_is_worker .and. numucat > 0) acc_mt_inp(imt, :) = tmpvec(:)

      ENDDO

      CALL vector_read_and_scatter(file_restart, tmpvec, numucat, 'acc_rnof_ref', ucat_data_address)
      IF (p_is_worker .and. numucat > 0) acc_rnof_ref(:) = tmpvec(:)

      ! History window.  The rows persisted by mt_riverlake_write_restart are
      ! restored when the restart sits in the MIDDLE of a window (monthly
      ! history with daily restarts, say): the first record after the restart
      ! then covers the same period a continuous run would report.
      !
      ! A restart on a window boundary carries a COMPLETE window instead, and
      ! restoring that would make the first window after the restart cover two
      ! windows.  Measured: before this rule, a daily restart reported
      ! 0.5876 kg N/s where the continuous run reported 0.7247, i.e. exactly
      ! (0.4505 + 0.7247)/2 = day1 and day2 averaged, and the restart file
      ! itself carried mt_hist_acctime = 86400 (172800 on the second day) --
      ! the window had not been flushed at the moment the restart was written,
      ! although hist_out runs before WRITE_TimeVariables in the main loop.
      !
      ! Which case applies is decided once, by MOD_HistAccRestart, from the
      ! restart date (a restart on mday == 1 closes a month, and so on); the
      ! material rows simply follow that verdict through
      ! mt_riverlake_drop_hist_window, so the two histories cannot disagree.
      ! History is one unit: partial optional legacy rows must not be mixed
      ! with a full-window denominator or with another row's history window.
      hist_complete = .true.
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         write(varname, '(A,A)') 'mt_hist_stor_', trim(material_names(imt))
         CALL probe_riverlake_restart_vector(file_restart, trim(varname), .false., field_present)
         hist_complete = hist_complete .and. field_present
         write(varname, '(A,A)') 'mt_hist_levsto_', trim(material_names(imt))
         CALL probe_riverlake_restart_vector(file_restart, trim(varname), .false., field_present)
         hist_complete = hist_complete .and. field_present
         write(varname, '(A,A)') 'mt_hist_out_', trim(material_names(imt))
         CALL probe_riverlake_restart_vector(file_restart, trim(varname), .false., field_present)
         hist_complete = hist_complete .and. field_present
         write(varname, '(A,A)') 'mt_hist_bifout_', trim(material_names(imt))
         CALL probe_riverlake_restart_vector(file_restart, trim(varname), .false., field_present)
         hist_complete = hist_complete .and. field_present
      ENDDO
      CALL probe_riverlake_restart_vector(file_restart, 'mt_hist_water_storage', .false., field_present)
      hist_complete = hist_complete .and. field_present
      CALL probe_riverlake_restart_vector(file_restart, 'mt_hist_levsto_water', .false., field_present)
      hist_complete = hist_complete .and. field_present
      IF (p_is_worker) THEN
         a_mt_storage_mass = 0._r8
         a_water_storage = 0._r8
         a_mt_levsto_mass = 0._r8
         a_levsto_water = 0._r8
         a_mt_out = 0._r8
         a_mt_bifout = 0._r8
         a_mt_acctime = 0._r8
      ENDIF
      IF (.not. hist_complete) THEN
         IF (p_is_master) WRITE(*,'(A)') &
            '  NOTE mt_riverlake_read_restart: legacy material history starts a post-restart window.'
      ELSE
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         write(varname, '(A,A)') 'mt_hist_stor_', trim(material_names(imt))
         CALL read_mt_hist_acc (file_restart, trim(varname), acc2=a_mt_storage_mass, imt=imt)
         write(varname, '(A,A)') 'mt_hist_levsto_', trim(material_names(imt))
         CALL read_mt_hist_acc (file_restart, trim(varname), acc2=a_mt_levsto_mass, imt=imt)
         write(varname, '(A,A)') 'mt_hist_out_', trim(material_names(imt))
         CALL read_mt_hist_acc (file_restart, trim(varname), acc2=a_mt_out, imt=imt)
         write(varname, '(A,A)') 'mt_hist_bifout_', trim(material_names(imt))
         CALL read_mt_hist_acc (file_restart, trim(varname), acc2=a_mt_bifout, imt=imt)
         write(varname, '(A,A)') 'mt_hist_n2o_', trim(material_names(imt))
         CALL probe_riverlake_restart_vector(file_restart, trim(varname), .false., hist_time_present)
         IF (hist_time_present) THEN
            CALL read_mt_hist_acc (file_restart, trim(varname), acc2=mt_n2o_source, imt=imt)
         ELSEIF (allocated(mt_n2o_source)) THEN
            ! Restarts written before the N2O diagnostic carry no numerator;
            ! start its window at zero instead of failing the reader.
            mt_n2o_source(imt, :) = 0._r8
         ENDIF
      ENDDO
      CALL read_mt_hist_acc (file_restart, 'mt_hist_water_storage', acc1=a_water_storage)
      CALL read_mt_hist_acc (file_restart, 'mt_hist_levsto_water', acc1=a_levsto_water)
      CALL probe_riverlake_restart_vector(file_restart, 'mt_hist_acctime', .false., hist_time_present)
      IF (hist_time_present) THEN
         CALL read_mt_hist_acc (file_restart, 'mt_hist_acctime', acc1=a_mt_acctime)
      ELSEIF (p_is_worker .and. allocated(acctime_ucat)) THEN
         ! Transitional files with all six numerators but no local clock used
         ! the water history clock; preserve that complete window on upgrade.
         a_mt_acctime = acctime_ucat
      ENDIF
      ENDIF

      ! Reject corrupt loaded rows before the levee-configuration migration can
      ! fold protected mass into the visible pool and thereby hide its origin.
      ! This same validator is run again after levee migration below.
      CALL validate_riverlake_restart_state()

      ! The water restart reader folds protected storage into visible storage
      ! whenever the current configuration has no levee for a cell.  Mirror
      ! that transfer for material mass; simply zeroing mt_levsto would either
      ! strand mass (global LEVEE off) or lose it (per-cell mask change).
      IF (p_is_worker .and. numucat > 0 .and. allocated(mt_mass) &
          .and. allocated(mt_levsto) .and. allocated(has_levee_bf)) THEN
         DO ii_bf = 1, numucat
            IF (ii_bf <= size(has_levee_bf)) THEN
               IF (.not. has_levee_bf(ii_bf)) THEN
                  DO imt_bf = 1, nmts
                     IF (.not. mt_uses_generic_water(imt_bf)) CYCLE
                     mt_mass(imt_bf, ii_bf) = mt_mass(imt_bf, ii_bf) &
                        + mt_levsto(imt_bf, ii_bf)
                     mt_levsto(imt_bf, ii_bf) = 0._r8
                     IF (allocated(mt_solid)) THEN
                        IF (mt_has_dissolved_limit(imt_bf)) THEN
                           mt_solid(imt_bf, ii_bf) = mt_solid(imt_bf, ii_bf) &
                              + mt_levsto_solid(imt_bf, ii_bf)
                           mt_levsto_solid(imt_bf, ii_bf) = 0._r8
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
      ENDIF

      ! Levee migration may have added two large finite pools. Recheck for
      ! overflow after the transformation before transport can consume state.
      CALL validate_riverlake_restart_state()

      deallocate (tmpvec)

      IF (present(found_restart)) found_restart = .true.

   END SUBROUTINE mt_riverlake_read_restart


   !-------------------------------------------------------------------------------------
   ! One history accumulator row <-> one ucatch restart vector.  Every rank must
   ! enter (vector_gather_and_write / vector_read_and_scatter are collective).
   ! An accumulator that is not allocated is written as zero.
   !-------------------------------------------------------------------------------------
   SUBROUTINE write_mt_hist_acc (file_restart, varname, acc2, acc1, imt)

   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_data_address
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart, varname
   real(r8), allocatable, intent(in), optional :: acc2(:,:)
   real(r8), allocatable, intent(in), optional :: acc1(:)
   integer, intent(in), optional :: imt

   real(r8), allocatable :: v(:)

      IF (p_is_worker .and. numucat > 0) THEN
         allocate (v(numucat)); v = 0._r8
         IF (present(acc2) .and. present(imt)) THEN
            IF (allocated(acc2)) v(:) = acc2(imt, :)
         ELSEIF (present(acc1)) THEN
            IF (allocated(acc1)) v(:) = acc1(:)
         ENDIF
      ELSE
         allocate (v(0))
      ENDIF

      CALL vector_gather_and_write (v, numucat, totalnumucat, ucat_data_address, &
         file_restart, varname, 'ucatch')
      deallocate (v)

   END SUBROUTINE write_mt_hist_acc

   SUBROUTINE read_mt_hist_acc (file_restart, varname, acc2, acc1, imt)

   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_data_address
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart, varname
   real(r8), allocatable, intent(inout), optional :: acc2(:,:)
   real(r8), allocatable, intent(inout), optional :: acc1(:)
   integer, intent(in), optional :: imt

   logical :: on_disk
   real(r8), allocatable :: v(:)

      ! Older restarts do not carry these rows; the accumulators then simply
      ! start from zero as before.
      CALL probe_riverlake_restart_vector (file_restart, varname, .false., on_disk)
      IF (.not. on_disk) RETURN

      IF (p_is_worker .and. numucat > 0) THEN
         allocate (v(numucat))
      ELSE
         allocate (v(0))
      ENDIF

      CALL vector_read_and_scatter (file_restart, v, numucat, varname, ucat_data_address)

      IF (p_is_worker .and. numucat > 0) THEN
         ! Diagnostics only: a non-finite value must not poison later records.
         WHERE (.not. ieee_is_finite(v)) v = 0._r8
         IF (present(acc2) .and. present(imt)) THEN
            IF (allocated(acc2)) acc2(imt, :) = v(:)
         ELSEIF (present(acc1)) THEN
            IF (allocated(acc1)) acc1(:) = v(:)
         ENDIF
      ENDIF
      deallocate (v)

   END SUBROUTINE read_mt_hist_acc

   SUBROUTINE mt_riverlake_snapshot_hist_acc ()
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE
   integer :: nmt, nuc

      IF (.not. p_is_worker) RETURN
      IF (nmts <= 0) RETURN
      nmt = nmts
      nuc = 0
      IF (numucat > 0) nuc = numucat

      IF (allocated(histsnap_stor)) deallocate(histsnap_stor, histsnap_levsto, histsnap_out, &
         histsnap_bifout, histsnap_n2o, histsnap_wstor, histsnap_wlevsto, histsnap_acctime)
      allocate(histsnap_stor(nmt,nuc), histsnap_levsto(nmt,nuc), histsnap_out(nmt,nuc), &
         histsnap_bifout(nmt,nuc), histsnap_n2o(nmt,nuc), histsnap_wstor(nuc), &
         histsnap_wlevsto(nuc), histsnap_acctime(nuc))
      histsnap_stor = 0._r8; histsnap_levsto = 0._r8; histsnap_out = 0._r8
      histsnap_bifout = 0._r8; histsnap_n2o = 0._r8
      histsnap_wstor = 0._r8; histsnap_wlevsto = 0._r8
      histsnap_acctime = 0._r8
      IF (nuc > 0) THEN
         IF (allocated(a_mt_storage_mass)) histsnap_stor = a_mt_storage_mass
         IF (allocated(a_mt_levsto_mass)) histsnap_levsto = a_mt_levsto_mass
         IF (allocated(a_mt_out)) histsnap_out = a_mt_out
         IF (allocated(a_mt_bifout)) histsnap_bifout = a_mt_bifout
         IF (allocated(mt_n2o_source)) histsnap_n2o = mt_n2o_source
         IF (allocated(a_water_storage)) histsnap_wstor = a_water_storage
         IF (allocated(a_levsto_water)) histsnap_wlevsto = a_levsto_water
         IF (allocated(a_mt_acctime)) histsnap_acctime = a_mt_acctime
      ENDIF
      histsnap_valid = .true.

   END SUBROUTINE mt_riverlake_snapshot_hist_acc

   SUBROUTINE mt_riverlake_write_hist_acc_restart(file)
   IMPLICIT NONE
   character(len=*), intent(in) :: file
   character(len=128) :: name
   integer :: imt

      IF (nmts <= 0) RETURN
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         write(name, '(A,A)') 'mt_hist_stor_', trim(material_names(imt))
         CALL write_mt_hist_acc(file, trim(name), acc2=a_mt_storage_mass, imt=imt)
         write(name, '(A,A)') 'mt_hist_levsto_', trim(material_names(imt))
         CALL write_mt_hist_acc(file, trim(name), acc2=a_mt_levsto_mass, imt=imt)
         write(name, '(A,A)') 'mt_hist_out_', trim(material_names(imt))
         CALL write_mt_hist_acc(file, trim(name), acc2=a_mt_out, imt=imt)
         write(name, '(A,A)') 'mt_hist_bifout_', trim(material_names(imt))
         CALL write_mt_hist_acc(file, trim(name), acc2=a_mt_bifout, imt=imt)
         write(name, '(A,A)') 'mt_hist_n2o_', trim(material_names(imt))
         CALL write_mt_hist_acc(file, trim(name), acc2=mt_n2o_source, imt=imt)
      ENDDO
      CALL write_mt_hist_acc(file, 'mt_hist_water_storage', acc1=a_water_storage)
      CALL write_mt_hist_acc(file, 'mt_hist_levsto_water', acc1=a_levsto_water)
      CALL write_mt_hist_acc(file, 'mt_hist_acctime', acc1=a_mt_acctime)
   END SUBROUTINE mt_riverlake_write_hist_acc_restart

   SUBROUTINE mt_riverlake_read_hist_acc_restart(file)
   IMPLICIT NONE
   character(len=*), intent(in) :: file
   character(len=128) :: name
   logical :: present_on_disk
   integer :: imt

      IF (nmts <= 0) RETURN
      ! A committed sidecar is indivisible: fail before scattering any row.
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         write(name, '(A,A)') 'mt_hist_stor_', trim(material_names(imt))
         CALL probe_riverlake_restart_vector(file, trim(name), .true., present_on_disk)
         write(name, '(A,A)') 'mt_hist_levsto_', trim(material_names(imt))
         CALL probe_riverlake_restart_vector(file, trim(name), .true., present_on_disk)
         write(name, '(A,A)') 'mt_hist_out_', trim(material_names(imt))
         CALL probe_riverlake_restart_vector(file, trim(name), .true., present_on_disk)
         write(name, '(A,A)') 'mt_hist_bifout_', trim(material_names(imt))
         CALL probe_riverlake_restart_vector(file, trim(name), .true., present_on_disk)
      ENDDO
      CALL probe_riverlake_restart_vector(file, 'mt_hist_water_storage', .true., present_on_disk)
      CALL probe_riverlake_restart_vector(file, 'mt_hist_levsto_water', .true., present_on_disk)
      CALL probe_riverlake_restart_vector(file, 'mt_hist_acctime', .true., present_on_disk)
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         write(name, '(A,A)') 'mt_hist_stor_', trim(material_names(imt))
         CALL read_mt_hist_acc(file, trim(name), acc2=a_mt_storage_mass, imt=imt)
         write(name, '(A,A)') 'mt_hist_levsto_', trim(material_names(imt))
         CALL read_mt_hist_acc(file, trim(name), acc2=a_mt_levsto_mass, imt=imt)
         write(name, '(A,A)') 'mt_hist_out_', trim(material_names(imt))
         CALL read_mt_hist_acc(file, trim(name), acc2=a_mt_out, imt=imt)
         write(name, '(A,A)') 'mt_hist_bifout_', trim(material_names(imt))
         CALL read_mt_hist_acc(file, trim(name), acc2=a_mt_bifout, imt=imt)
         write(name, '(A,A)') 'mt_hist_n2o_', trim(material_names(imt))
         CALL probe_riverlake_restart_vector(file, trim(name), .false., present_on_disk)
         IF (present_on_disk) THEN
            CALL read_mt_hist_acc(file, trim(name), acc2=mt_n2o_source, imt=imt)
         ELSEIF (allocated(mt_n2o_source)) THEN
            ! Older sidecars predate the N2O diagnostic; start its window
            ! at zero instead of failing the whole history transaction.
            mt_n2o_source(imt, :) = 0._r8
         ENDIF
      ENDDO
      CALL read_mt_hist_acc(file, 'mt_hist_water_storage', acc1=a_water_storage)
      CALL read_mt_hist_acc(file, 'mt_hist_levsto_water', acc1=a_levsto_water)
      CALL read_mt_hist_acc(file, 'mt_hist_acctime', acc1=a_mt_acctime)
   END SUBROUTINE mt_riverlake_read_hist_acc_restart

   !-------------------------------------------------------------------------------------
   ! Write material restart
   !-------------------------------------------------------------------------------------
   SUBROUTINE mt_riverlake_write_restart (file_restart)

   USE MOD_Vector_ReadWrite
   USE MOD_NetCDFSerial, only: ncio_define_dimension, ncio_write_serial
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_data_address, &
      ucat_gdid, ucat_next
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

   integer :: imt
   character(len=64) :: varname
   real(r8), allocatable :: tmpvec(:)
   real(r8), allocatable :: t2a(:,:), t2b(:,:), t2c(:,:), t2d(:,:), t2e(:,:)
   real(r8), allocatable :: t1a(:), t1b(:), t1c(:)
   integer, allocatable :: descriptor_identity(:,:)

      ! Guard: material module may not be initialised (e.g. mkinidata).
      ! material_names is allocated on ALL ranks by mt_riverlake_init (before the
      ! p_is_worker guard), so this check is safe for master/IO.
      ! Do NOT check allocated(mt_mass) here — mt_mass is worker-only,
      ! but vector_gather_and_write has mpi_barrier(p_comm_glb) that
      ! requires ALL ranks to participate.
      IF (.not. allocated(material_names)) RETURN
      ! Never return on a worker-only allocation predicate: every rank must
      ! enter the gather collectives below in the same order. Initialization
      ! allocates zero-length routing state on empty workers.

      ! A count/name hash cannot distinguish materials whose units, ownership, or
      ! physical parameters changed. Build the exact shared descriptor once;
      ! it is committed only after every distributed state field is durable.
      CALL material_build_descriptor_identity(descriptor_identity, transport_only=.true.)
      IF (size(descriptor_identity, 2) <= 0) THEN
         ! Provider-only configurations carry no generic state. Commit an
         ! explicit empty transaction instead of leaving complete=0: a later
         ! run that adds a transport material must cold-start cleanly rather than
         ! misclassifying this intentional absence as an interrupted write.
         IF (p_is_master) THEN
            CALL ncio_write_serial(file_restart, 'mt_river_restart_complete', 0)
            CALL ncio_write_serial(file_restart, 'mt_river_descriptor_count', 0)
            CALL ncio_write_serial(file_restart, 'mt_river_restart_schema', &
               RIVER_MT_RESTART_SCHEMA_VERSION)
            CALL ncio_write_serial(file_restart, 'mt_river_restart_complete', 1)
         ENDIF
         deallocate(descriptor_identity)
         RETURN
      ENDIF

      ! Never persist a poisoned checkpoint. This is the same collective state
      ! contract enforced after reads, so all ranks must enter before master I/O.
      CALL validate_riverlake_restart_state()
      IF (p_is_master) THEN
         CALL ncio_write_serial(file_restart, 'mt_river_restart_complete', 0)
      ENDIF

      IF (p_is_worker .and. numucat > 0) THEN
         allocate (tmpvec(numucat))
      ELSE
         allocate (tmpvec(0))
      ENDIF

      ! Restart metadata guards. Store as ucatch vectors so the existing
      ! vector I/O path can scatter and verify them collectively on read.
      IF (p_is_worker .and. numucat > 0) tmpvec(:) = real(totalnumucat, r8)
      CALL vector_gather_and_write ( &
         tmpvec, numucat, totalnumucat, ucat_data_address, file_restart, 'mt_numucat_meta', 'ucatch')

      IF (p_is_worker .and. numucat > 0) tmpvec(:) = real(ucat_gdid(:), r8)
      CALL vector_gather_and_write ( &
         tmpvec, numucat, totalnumucat, ucat_data_address, file_restart, 'mt_ucat_gdid_meta', 'ucatch')

      IF (p_is_worker .and. numucat > 0) tmpvec(:) = real(ucat_next(:), r8)
      CALL vector_gather_and_write ( &
         tmpvec, numucat, totalnumucat, ucat_data_address, file_restart, 'mt_ucat_next_meta', 'ucatch')

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         IF (p_is_worker .and. numucat > 0) THEN
            tmpvec(:) = mt_mass(imt, :)
         ENDIF

         write(varname, '(A,A)') 'mt_mass_', trim(material_names(imt))
         CALL vector_gather_and_write ( &
            tmpvec, numucat, totalnumucat, ucat_data_address, file_restart, trim(varname), 'ucatch')

         IF (p_is_worker .and. numucat > 0) THEN
            tmpvec(:) = mt_inp_buf(imt, :)
         ENDIF
         write(varname, '(A,A)') 'mt_inpbuf_', trim(material_names(imt))
         CALL vector_gather_and_write ( &
            tmpvec, numucat, totalnumucat, ucat_data_address, file_restart, trim(varname), 'ucatch')

         ! Per-material routing-period accumulator. Paired with acc_rnof_ref
         ! below; both must persist so a restart written mid-period
         ! recovers the queued land-material input instead of dropping it.
         IF (p_is_worker .and. numucat > 0) THEN
            tmpvec(:) = acc_mt_inp(imt, :)
         ENDIF
         write(varname, '(A,A)') 'mt_accinp_', trim(material_names(imt))
         CALL vector_gather_and_write ( &
            tmpvec, numucat, totalnumucat, ucat_data_address, file_restart, trim(varname), 'ucatch')

         ! Always persist the protected-side row, even when levees are disabled
         ! (then it is zero). This makes schema completeness independent of a
         ! runtime levee toggle and preserves mass when a leveed restart is read
         ! by a non-leveed configuration that folds the pool into visible water.
         IF (p_is_worker .and. numucat > 0) THEN
            tmpvec(:) = mt_levsto(imt, :)
         ENDIF
         write(varname, '(A,A)') 'mt_levsto_', trim(material_names(imt))
         CALL vector_gather_and_write ( &
            tmpvec, numucat, totalnumucat, ucat_data_address, file_restart, trim(varname), 'ucatch')
         IF (mt_has_dissolved_limit(imt)) THEN
            IF (p_is_worker .and. numucat > 0) tmpvec(:) = mt_solid(imt, :)
            write(varname, '(A,A)') 'mt_solid_', trim(material_names(imt))
            CALL vector_gather_and_write ( &
               tmpvec, numucat, totalnumucat, ucat_data_address, file_restart, trim(varname), 'ucatch')
            IF (p_is_worker .and. numucat > 0) tmpvec(:) = mt_levsto_solid(imt, :)
            write(varname, '(A,A)') 'mt_levsto_solid_', trim(material_names(imt))
            CALL vector_gather_and_write ( &
               tmpvec, numucat, totalnumucat, ucat_data_address, file_restart, trim(varname), 'ucatch')
         ENDIF
      ENDDO

      ! Shared runoff reference paired with acc_mt_inp for the active
      ! routing period. Written once (not per material).
      IF (p_is_worker .and. numucat > 0) THEN
         tmpvec(:) = acc_rnof_ref(:)
      ENDIF
      CALL vector_gather_and_write ( &
         tmpvec, numucat, totalnumucat, ucat_data_address, file_restart, 'acc_rnof_ref', 'ucatch')

      ! Persist the numerator rows and their material-local elapsed time as one
      ! history window. Older files without rows start a post-restart window.
      ! At a terminal (partial-window) stop the live rows were cleared by the
      ! history flush; swap the pre-flush snapshot in for these rows only.
      IF (histsnap_valid .and. p_is_worker .and. numucat > 0 .and. &
          allocated(histsnap_stor) .and. allocated(a_mt_storage_mass)) THEN
         allocate(t2a(nmts,numucat), t2b(nmts,numucat), t2c(nmts,numucat), &
            t2d(nmts,numucat), t2e(nmts,numucat), t1a(numucat), t1b(numucat), t1c(numucat))
         t2a = a_mt_storage_mass; t2b = a_mt_levsto_mass; t2c = a_mt_out; t2d = a_mt_bifout
         IF (allocated(mt_n2o_source)) THEN
            t2e = mt_n2o_source
            IF (allocated(histsnap_n2o)) mt_n2o_source = histsnap_n2o
         ENDIF
         t1a = a_water_storage; t1b = a_levsto_water; t1c = a_mt_acctime
         a_mt_storage_mass = histsnap_stor; a_mt_levsto_mass = histsnap_levsto
         a_mt_out = histsnap_out; a_mt_bifout = histsnap_bifout
         a_water_storage = histsnap_wstor; a_levsto_water = histsnap_wlevsto
         a_mt_acctime = histsnap_acctime
      ENDIF
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         write(varname, '(A,A)') 'mt_hist_stor_', trim(material_names(imt))
         CALL write_mt_hist_acc (file_restart, trim(varname), acc2=a_mt_storage_mass, imt=imt)
         write(varname, '(A,A)') 'mt_hist_levsto_', trim(material_names(imt))
         CALL write_mt_hist_acc (file_restart, trim(varname), acc2=a_mt_levsto_mass, imt=imt)
         write(varname, '(A,A)') 'mt_hist_out_', trim(material_names(imt))
         CALL write_mt_hist_acc (file_restart, trim(varname), acc2=a_mt_out, imt=imt)
         write(varname, '(A,A)') 'mt_hist_bifout_', trim(material_names(imt))
         CALL write_mt_hist_acc (file_restart, trim(varname), acc2=a_mt_bifout, imt=imt)
         write(varname, '(A,A)') 'mt_hist_n2o_', trim(material_names(imt))
         CALL write_mt_hist_acc (file_restart, trim(varname), acc2=mt_n2o_source, imt=imt)
      ENDDO
      CALL write_mt_hist_acc (file_restart, 'mt_hist_water_storage', acc1=a_water_storage)
      CALL write_mt_hist_acc (file_restart, 'mt_hist_levsto_water', acc1=a_levsto_water)
      CALL write_mt_hist_acc (file_restart, 'mt_hist_acctime', acc1=a_mt_acctime)
      IF (allocated(t2a)) THEN
         a_mt_storage_mass = t2a; a_mt_levsto_mass = t2b; a_mt_out = t2c; a_mt_bifout = t2d
         IF (allocated(mt_n2o_source)) mt_n2o_source = t2e
         a_water_storage = t1a; a_levsto_water = t1b; a_mt_acctime = t1c
         deallocate(t2a, t2b, t2c, t2d, t2e, t1a, t1b, t1c)
      ENDIF
      IF (histsnap_valid) THEN
         histsnap_valid = .false.
         IF (allocated(histsnap_stor)) deallocate(histsnap_stor, histsnap_levsto, histsnap_out, &
            histsnap_bifout, histsnap_n2o, histsnap_wstor, histsnap_wlevsto, histsnap_acctime)
      ENDIF

      ! Commit metadata last. If any state/descriptor write above is interrupted,
      ! complete remains zero and the reader refuses the mixed-generation file.
      IF (p_is_master) THEN
         CALL ncio_define_dimension(file_restart, 'mt_river_descriptor_field', &
            MATERIAL_DESCRIPTOR_IDENTITY_WIDTH)
         CALL ncio_define_dimension(file_restart, 'mt_river_transport_material', &
            size(descriptor_identity, 2))
         CALL ncio_write_serial(file_restart, 'mt_river_descriptor_count', &
            size(descriptor_identity, 2))
         CALL ncio_write_serial(file_restart, 'mt_river_descriptor_identity', &
            descriptor_identity, 'mt_river_descriptor_field', 'mt_river_transport_material')
         CALL ncio_write_serial(file_restart, 'mt_river_restart_schema', &
            RIVER_MT_RESTART_SCHEMA_VERSION)
         CALL ncio_write_serial(file_restart, 'mt_river_restart_complete', 1)
      ENDIF

      deallocate (tmpvec)
      deallocate (descriptor_identity)

   END SUBROUTINE mt_riverlake_write_restart


   !-------------------------------------------------------------------------------------
   ! Write material history output
   !-------------------------------------------------------------------------------------
   SUBROUTINE mt_riverlake_write_history (file_hist_ucat, itime_in_file_ucat, acctime_ucat_hist)

   USE MOD_Grid_RiverLakeNetwork, only: numucat, allups_mask_ucat
   ! Route history is dispatched, not written directly, so 'one' and
   ! 'block' cannot drift apart variable by variable.
   USE MOD_Grid_RiverLakeHistRoute, only: route_hist_write_ucat
   USE MOD_MT_Defs, only: material_uses_delta_diagnostics, materials, mt_tiny, &
      mt_delta_sanity_max
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   character(len=*), intent(in) :: file_hist_ucat
   integer,  intent(in) :: itime_in_file_ucat
   real(r8), intent(in) :: acctime_ucat_hist(:)  ! Water history time; retained for caller compatibility

   integer :: imt, i
   character(len=64) :: varname
   character(len=128) :: longname
   character(len=32) :: conc_word, conc_units, mass_units, flux_units
   real(r8), allocatable :: tmpvec(:)
   real(r8) :: ratio_loc, delta_loc
   real(r8), parameter :: mt_hist_fp_dust = 1.0e-12_r8

      IF (p_is_worker .and. numucat > 0) THEN
         allocate (tmpvec(numucat))
      ELSE
         allocate (tmpvec(0))
      ENDIF
      ! Output-side FP-dust clamp threshold. Same rationale as
      ! check_mt_state: mt_mass can carry sub-picomass alternating
      ! residue from upwind flux add/subtract, which propagates to
      ! mt_conc / a_mt_out / a_mt_bifout. The underlying
      ! state is preserved intact (so transport / limiter regressions
      ! remain visible via the WARNING watchdog); only this display
      ! slice hides the noise floor.
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         IF (material_uses_delta_diagnostics(imt)) THEN
            conc_word  = 'ratio'
            conc_units = 'R'
            mass_units = 'R*m3'
            flux_units = 'R*m3/s'
         ELSE
            conc_word  = 'concentration'
            conc_units = material_concentration_units(imt)
            mass_units = 'material'
            flux_units = 'material/s'
         ENDIF

         ! --- Material concentration ---
         ! Water-weighted mean: sum(mass*dt) / sum(volume*dt), the same quotient
         ! the delta diagnostic below is built from.  A plain time average of the
         ! per-step ratio would disagree with that delta for the very same field
         ! (for isotopes conc_units is 'R', i.e. they report the same quantity),
         ! and would fold dry sub-steps -- where mt_conc is forced to zero and
         ! concentration is undefined, not zero -- into the mean.  Both
         ! accumulators below are gated on volwater > mt_v_dry_off, so dry time
         ! leaves numerator and denominator alike.
         IF (p_is_worker .and. numucat > 0) THEN
            tmpvec(:) = spval
            IF (allocated(a_mt_storage_mass) .and. allocated(a_water_storage)) THEN
               DO i = 1, numucat
                  IF (i > size(a_mt_acctime)) CYCLE
                  IF (a_mt_acctime(i) <= 0._r8) CYCLE
                  IF (a_water_storage(i) <= mt_delta_diag_vmin * a_mt_acctime(i)) CYCLE
                  tmpvec(i) = a_mt_storage_mass(imt, i) / a_water_storage(i)
               ENDDO
            ENDIF
            WHERE (abs(tmpvec) < mt_hist_fp_dust) tmpvec = 0._r8
         ENDIF

         write(varname, '(A,A)') 'f_mt_conc_', trim(material_names(imt))
         write(longname, '(5A)') 'material ', trim(conc_word), &
            ' (', trim(material_names(imt)), ')'

         CALL route_hist_write_ucat (tmpvec, trim(varname), &
            longname = trim(longname), &
            units = trim(conc_units))

         IF (material_uses_delta_diagnostics(imt)) THEN
            IF (p_is_worker .and. numucat > 0) THEN
               tmpvec(:) = spval
               IF (materials(imt)%ref_ratio > mt_tiny) THEN
                  DO i = 1, numucat
                     IF (allocated(allups_mask_ucat)) THEN
                        IF (i > size(allups_mask_ucat)) CYCLE
                        IF (allups_mask_ucat(i) < 0.5_r8) CYCLE
                     ENDIF
                     IF (i > size(a_mt_acctime)) CYCLE
                     IF (a_mt_acctime(i) <= 0._r8) CYCLE
                     IF (.not. allocated(a_water_storage)) CYCLE
                     IF (.not. allocated(a_mt_storage_mass)) CYCLE
                     IF (a_water_storage(i) <= mt_delta_diag_vmin * a_mt_acctime(i)) CYCLE
                     ratio_loc = a_mt_storage_mass(imt, i) / a_water_storage(i)
                     IF (ratio_loc <= mt_tiny) CYCLE
                     delta_loc = (ratio_loc / materials(imt)%ref_ratio - 1.0_r8) * 1000.0_r8
                     IF (abs(delta_loc) <= mt_delta_sanity_max) tmpvec(i) = delta_loc
                  ENDDO
               ENDIF
            ENDIF

            write(varname, '(A,A)') 'f_mt_delta_', trim(material_names(imt))
            write(longname, '(A,A,A)') 'material delta (', trim(material_names(imt)), ')'
            CALL route_hist_write_ucat (tmpvec, trim(varname), &
               longname = trim(longname), &
               units = 'permil')
         ENDIF

         ! --- Material outflux ---
         IF (p_is_worker .and. numucat > 0) THEN
            tmpvec(:) = spval
            DO i = 1, numucat
               IF (i <= size(a_mt_acctime)) THEN
                  IF (a_mt_acctime(i) > 0._r8) tmpvec(i) = a_mt_out(imt, i) / a_mt_acctime(i)
               ENDIF
            ENDDO
            WHERE (abs(tmpvec) < mt_hist_fp_dust) tmpvec = 0._r8
         ENDIF

         write(varname, '(A,A)') 'f_mt_flux_', trim(material_names(imt))
         write(longname, '(A,A,A)') 'material outflux (', trim(material_names(imt)), ')'

         CALL route_hist_write_ucat (tmpvec, trim(varname), &
            longname = trim(longname), &
            units = trim(flux_units))

         ! The row is always written (interior values may be zero) so a
         ! configured EF5 is visible in the history rather than silently
         ! absent.  route_hist_write_ucat is collective: master/IO ranks must
         ! enter it too, so it is NOT guarded by allocated(mt_n2o_source)
         ! (that array exists only on workers and only for numucat > 0).
         IF (p_is_worker .and. numucat > 0) THEN
            tmpvec(:) = spval
            DO i = 1, numucat
               IF (i <= size(a_mt_acctime)) THEN
                  IF (a_mt_acctime(i) > 0._r8) &
                     tmpvec(i) = mt_n2o_source(imt, i) / a_mt_acctime(i)
               ENDIF
            ENDDO
            WHERE (abs(tmpvec) < mt_hist_fp_dust) tmpvec = 0._r8
         ENDIF

         write(varname, '(A,A)') 'f_mt_n2o_flux_', trim(material_names(imt))
         write(longname, '(A,A,A)') 'N2O-N emitted at river entry (', &
            trim(material_names(imt)), ')'
         CALL route_hist_write_ucat (tmpvec, trim(varname), &
            longname = trim(longname), &
            units = 'kg N/s')

         IF (DEF_USE_LEVEE) THEN
            ! --- Protected-side material storage behind levees ---
            IF (p_is_worker .and. numucat > 0) THEN
               tmpvec(:) = 0._r8
               IF (allocated(a_mt_levsto_mass)) THEN
                  tmpvec(:) = spval
                  DO i = 1, numucat
                     IF (i <= size(a_mt_acctime)) THEN
                        IF (a_mt_acctime(i) > 0._r8) &
                           tmpvec(i) = a_mt_levsto_mass(imt, i) / a_mt_acctime(i)
                     ENDIF
                  ENDDO
               ENDIF
               WHERE (abs(tmpvec) < mt_hist_fp_dust) tmpvec = 0._r8
            ENDIF

            write(varname, '(A,A)') 'f_mt_levsto_', trim(material_names(imt))
            write(longname, '(A,A,A)') 'protected-side levee material storage (', trim(material_names(imt)), ')'

            CALL route_hist_write_ucat (tmpvec, trim(varname), &
               longname = trim(longname), &
               units = trim(mass_units))

            IF (material_uses_delta_diagnostics(imt)) THEN
               IF (p_is_worker .and. numucat > 0) THEN
                  tmpvec(:) = spval
                  IF (materials(imt)%ref_ratio > mt_tiny .and. &
                      allocated(a_mt_levsto_mass) .and. allocated(a_levsto_water)) THEN
                     DO i = 1, numucat
                        IF (allocated(allups_mask_ucat)) THEN
                           IF (i > size(allups_mask_ucat)) CYCLE
                           IF (allups_mask_ucat(i) < 0.5_r8) CYCLE
                        ENDIF
                        IF (i > size(a_mt_acctime)) CYCLE
                        IF (a_mt_acctime(i) <= 0._r8) CYCLE
                        IF (a_levsto_water(i) <= mt_delta_diag_vmin * a_mt_acctime(i)) CYCLE
                        ratio_loc = a_mt_levsto_mass(imt, i) / a_levsto_water(i)
                        IF (ratio_loc <= mt_tiny) CYCLE
                        delta_loc = (ratio_loc / materials(imt)%ref_ratio - 1.0_r8) * 1000.0_r8
                        IF (abs(delta_loc) <= mt_delta_sanity_max) tmpvec(i) = delta_loc
                     ENDDO
                  ENDIF
               ENDIF

               write(varname, '(A,A)') 'f_mt_levdelta_', trim(material_names(imt))
               write(longname, '(A,A,A)') 'protected-side levee material delta (', trim(material_names(imt)), ')'
               CALL route_hist_write_ucat (tmpvec, trim(varname), &
                  longname = trim(longname), &
                  units = 'permil')
            ENDIF
         ENDIF

         IF (DEF_USE_BIFURCATION) THEN
            ! --- Material bifurcation net flux ---
            IF (p_is_worker .and. numucat > 0) THEN
               tmpvec(:) = spval
               DO i = 1, numucat
                  IF (i <= size(a_mt_acctime)) THEN
                     IF (a_mt_acctime(i) > 0._r8) tmpvec(i) = a_mt_bifout(imt, i) / a_mt_acctime(i)
                  ENDIF
               ENDDO
               WHERE (abs(tmpvec) < mt_hist_fp_dust) tmpvec = 0._r8
            ENDIF

            write(varname, '(A,A)') 'f_mt_bifout_', trim(material_names(imt))
            write(longname, '(A,A,A)') 'material net bifurcation outflux (', trim(material_names(imt)), ')'

            CALL route_hist_write_ucat (tmpvec, trim(varname), &
               longname = trim(longname), &
               units = trim(flux_units))
         ENDIF

      ENDDO

      deallocate (tmpvec)

   END SUBROUTINE mt_riverlake_write_history


   !-------------------------------------------------------------------------------------
   ! Diagnostic check: print min/max of material state variables
   !-------------------------------------------------------------------------------------
#ifdef RangeCheck
   SUBROUTINE check_mt_state ()

   USE MOD_SPMD_Task
   USE MOD_RangeCheck
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   USE MOD_Grid_RiverLakeTimeVars, only: wdsrf_ucat, volresv_state => volresv
   USE MOD_Grid_Reservoir, only: ucat2resv_state => ucat2resv
   USE MOD_MT_Defs, only: materials, mt_tiny, mt_delta_sanity_max, &
      material_uses_delta_diagnostics
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   integer :: imt
   integer :: icell_worst, nbad, i
   character(len=64) :: label
   character(len=16) :: mt_mass_units, mt_conc_units, mt_flux_units
   real(r8), allocatable :: tmp(:)
   real(r8), allocatable :: volresv_check(:)
   integer, allocatable :: ucat2resv_check(:)
   real(r8) :: worst_neg_mass
   real(r8) :: volwater, delta_loc
   real(r8), parameter :: mt_mass_fp_dust  = 1.0e-12_r8   ! [R*m3]
   ! Same magnitude, different dimension. Declared separately so the flux
   ! path cannot be read as comparing a rate against an absolute mass, and
   ! so the two can be tuned independently if the residue ever differs.
   real(r8), parameter :: mt_flux_fp_dust  = 1.0e-12_r8   ! [R*m3/s]
   real(r8), parameter :: mt_mass_neg_warn = -1.0e-6_r8

      ! Workers that never ran mt_riverlake_init have no data to check.
      ! Master and IO ranks don't allocate mt_mass but MUST still enter
      ! this routine: master participates in check_vector_data's recv/print
      ! side, and all ranks share the same call count to avoid orphaning
      ! worker-to-master mpi_send messages.
      IF (p_is_worker .and. (.not. allocated(mt_mass) .or. .not. allocated(mt_conc))) RETURN

      IF (p_is_master) THEN
         write(*,'(/,A)') 'Checking Material Variables ...'
      ENDIF

      IF (p_is_worker .and. numucat > 0) THEN
         allocate (tmp(numucat))
      ELSE
         allocate (tmp(0))
      ENDIF

      ! Display-only thresholds (state itself is never modified here):
      !   mt_mass_fp_dust: absolute mt_mass below this is FP residue
      !       from alternating substep add/subtract (order 1e-15 R*m³);
      !       dust-clamp for log cleanliness only. Well below any
      !       physical input (rnof*R_init is typically >= 1e-9).
      !   mt_mass_neg_warn: mt_mass below this signed threshold is
      !       NOT noise — indicates transport/limiter drift. Emit a
      !       diagnostic once so the root cause can be traced instead
      !       of silently swallowed by a state-side clamp.
      IF (p_is_worker) THEN
         IF (allocated(volresv_state)) THEN
            allocate (volresv_check(size(volresv_state)))
            volresv_check(:) = volresv_state(:)
         ELSE
            allocate (volresv_check(0))
         ENDIF
         IF (allocated(ucat2resv_state)) THEN
            allocate (ucat2resv_check(size(ucat2resv_state)))
            ucat2resv_check(:) = ucat2resv_state(:)
         ELSE
            allocate (ucat2resv_check(0))
         ENDIF
      ENDIF

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         IF (material_uses_delta_diagnostics(imt)) THEN
            mt_mass_units = 'R*m3'
            mt_conc_units = 'R'
            mt_flux_units = 'R*m3/s'
         ELSE
            mt_mass_units = 'material'
            mt_conc_units = material_concentration_units(imt)
            mt_flux_units = 'material/s'
         ENDIF

         IF (p_is_worker .and. numucat > 0) THEN
            tmp = mt_mass(imt,:)
            ! Significant-negative watchdog (state unchanged)
            worst_neg_mass = 0._r8
            icell_worst    = 0
            nbad           = 0
            IF (numucat > 0) THEN
               icell_worst = minloc(tmp, dim=1)
               IF (icell_worst >= 1 .and. icell_worst <= numucat) THEN
                  IF (tmp(icell_worst) < mt_mass_neg_warn) THEN
                     worst_neg_mass = tmp(icell_worst)
                     nbad = count(tmp < mt_mass_neg_warn)
                  ENDIF
               ENDIF
            ENDIF
            IF (worst_neg_mass < mt_mass_neg_warn) THEN
               write(*,'(A,A,A,I6,A,I8,A,E12.5)') &
                  ' WARNING check_mt_state: ', trim(material_names(imt)), &
                  ' mt_mass < -1e-6 in ', nbad, ' cell(s); worst @ucat=', &
                  icell_worst, ' mass=', worst_neg_mass
            ENDIF
            ! Cosmetic dust clamp for the min/max line only
            WHERE (abs(tmp) < mt_mass_fp_dust) tmp = 0._r8
         ENDIF
         write(label,'(5A)') 'mt_mass_', trim(material_names(imt)), ' [', trim(mt_mass_units), ']'
         CALL check_vector_data (label, tmp)

         IF (p_is_worker .and. numucat > 0) THEN
            tmp = mt_conc(imt,:)
            ! Cosmetic dust clamp on the conc display path, keyed on the MASS
            ! that produced each concentration rather than on the
            ! concentration itself. What this suppresses is a mt_mass of
            ! FP-residue size propagating into mt_conc, with the sign flipped
            ! by division -- so the test belongs on mt_mass.
            !
            ! It used to compare mt_conc against mt_mass_fp_dust directly.
            ! That threshold is an absolute mass in R*m3; mt_conc is a ratio
            ! or a concentration per unit water, so the comparison had no
            ! defined meaning and the cells it caught depended on cell volume.
            ! Display only either way -- mt_conc itself is never modified.
            WHERE (abs(mt_mass(imt,:)) < mt_mass_fp_dust) tmp = 0._r8
         ENDIF
         write(label,'(5A)') 'mt_conc_', trim(material_names(imt)), ' [', trim(mt_conc_units), ']'
         CALL check_vector_data (label, tmp)

         IF (material_uses_delta_diagnostics(imt)) THEN
            IF (p_is_worker .and. numucat > 0) THEN
               tmp = spval

               IF (materials(imt)%ref_ratio > mt_tiny) THEN
                  ! tmp is rebuilt element by element in the loop below; the
                  ! whole-array delta that used to be computed here was
                  ! overwritten on the next line without ever being read.
                  tmp = spval
                  DO i = 1, numucat
                     volwater = 0._r8
                     IF (allocated(wdsrf_ucat)) THEN
                        IF (i <= size(wdsrf_ucat)) THEN
                           CALL get_cell_volume(i, wdsrf_ucat(i), volresv_check, ucat2resv_check, volwater)
                        ENDIF
                     ENDIF

                     IF (mt_conc(imt, i) > mt_tiny) THEN
                        delta_loc = (mt_conc(imt, i) / materials(imt)%ref_ratio - 1.0_r8) * 1000.0_r8

                        IF (volwater > mt_delta_diag_vmin) THEN
                           IF (abs(delta_loc) <= mt_delta_sanity_max) tmp(i) = delta_loc
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF
            write(label,'(A,A,A)') 'mt_delta_', trim(material_names(imt)), ' [permil]'
            CALL check_vector_data (label, tmp)
         ENDIF

         IF (p_is_worker .and. numucat > 0) THEN
            tmp = mt_flux_out(imt,:)
            ! FP-dust clamp on the flux display path (same rationale
            ! as mt_mass / mt_conc above). Real flux magnitudes are
            ! many orders above the residue so this only hides
            ! the alternating +/- 1e-15 numerical noise.
            WHERE (abs(tmp) < mt_flux_fp_dust) tmp = 0._r8
         ENDIF
         write(label,'(5A)') 'mt_outflux_', trim(material_names(imt)), ' [', trim(mt_flux_units), ']'
         CALL check_vector_data (label, tmp)

         IF (p_is_worker .and. numucat > 0) tmp = acc_mt_inp(imt,:)
         write(label,'(5A)') 'mt_inp_', trim(material_names(imt)), ' [', trim(mt_mass_units), ']'
         CALL check_vector_data (label, tmp)

         IF (p_is_worker .and. numucat > 0) tmp = mt_inp_buf(imt,:)
         write(label,'(5A)') 'mt_inpbuf_', trim(material_names(imt)), ' [', trim(mt_mass_units), ']'
         CALL check_vector_data (label, tmp)

         ! Surface protected-side pool so the levee-on case has a
         ! min/max line alongside mt_mass. Same FP-dust clamp rationale.
         IF (p_is_worker .and. numucat > 0) THEN
            IF (allocated(mt_levsto)) THEN
               tmp = mt_levsto(imt,:)
               WHERE (abs(tmp) < mt_mass_fp_dust) tmp = 0._r8
            ELSE
               tmp = 0._r8
            ENDIF
         ENDIF
         write(label,'(5A)') 'mt_levsto_', trim(material_names(imt)), ' [', trim(mt_mass_units), ']'
         CALL check_vector_data (label, tmp)
      ENDDO
      IF (allocated(volresv_check)) deallocate (volresv_check)
      IF (allocated(ucat2resv_check)) deallocate (ucat2resv_check)

      deallocate (tmp)

   END SUBROUTINE check_mt_state
#else
   SUBROUTINE check_mt_state ()
   IMPLICIT NONE
   END SUBROUTINE check_mt_state
#endif


   !-------------------------------------------------------------------------------------
   ! Deallocate material module
   !-------------------------------------------------------------------------------------
   SUBROUTINE mt_riverlake_final ()

   IMPLICIT NONE

   integer :: k
   real(r8) :: stats(4)

      ! Processed/delivered are cumulative handling throughput, so a mass
      ! parcel handled on multiple edges/steps contributes on each occasion.
      ! Pending is current inventory, not the sum of repeated parking events.
      !
      ! Every rank must enter the allreduce: mt_dry_stat exists only on
      ! workers, so guarding the collective with allocated() would hang the
      ! run at shutdown with the results already written.  Non-workers
      ! contribute zeros.
      IF (nmts > 0 .and. allocated(material_names)) THEN
         DO k = 1, nmts
            IF (.not. mt_uses_generic_water(k)) CYCLE
            stats(:) = 0._r8
            IF (p_is_worker .and. allocated(mt_dry_stat)) stats(1:4) = mt_dry_stat(k, 1:4)
            stats(3) = 0._r8
            IF (p_is_worker .and. allocated(mt_dry_pending)) stats(3) = sum(mt_dry_pending(k,:))
#ifdef USEMPI
            CALL mpi_allreduce (MPI_IN_PLACE, stats, 4, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
#endif
            IF (p_is_master) WRITE(*,'(A,A,A,4ES14.6)') ' Dry-cell accounting ', &
               trim(material_names(k)), ' processed/delivered/pending_inventory/sink:', &
               stats(1), stats(2), stats(3), stats(4)
         ENDDO
      ENDIF

      IF (allocated(material_names )) deallocate (material_names )
      IF (allocated(mt_mass     )) deallocate (mt_mass     )
      IF (allocated(mt_conc     )) deallocate (mt_conc     )
      IF (allocated(mt_flux_out )) deallocate (mt_flux_out )
      IF (allocated(acc_mt_inp       )) deallocate (acc_mt_inp       )
      IF (allocated(mt_inp_buf       )) deallocate (mt_inp_buf       )
      IF (allocated(acc_rnof_ref      )) deallocate (acc_rnof_ref      )
      IF (allocated(mt_bif_net_saved)) deallocate (mt_bif_net_saved)
      IF (allocated(a_mt_storage_mass)) deallocate (a_mt_storage_mass)
      IF (allocated(a_water_storage   )) deallocate (a_water_storage   )
      IF (allocated(a_mt_levsto_mass )) deallocate (a_mt_levsto_mass )
      IF (allocated(a_levsto_water    )) deallocate (a_levsto_water    )
      IF (allocated(a_mt_out        )) deallocate (a_mt_out        )
      IF (allocated(a_mt_bifout     )) deallocate (a_mt_bifout     )
      IF (allocated(a_mt_acctime    )) deallocate (a_mt_acctime    )
      IF (allocated(mt_dry_drain    )) deallocate (mt_dry_drain    )
      IF (allocated(mt_dry_stat     )) deallocate (mt_dry_stat     )
      IF (allocated(mt_dry_fwd      )) deallocate (mt_dry_fwd      )
      IF (allocated(mt_dry_rev      )) deallocate (mt_dry_rev      )
      IF (allocated(mt_dry_split    )) deallocate (mt_dry_split    )
      IF (allocated(mt_dry_pending  )) deallocate (mt_dry_pending  )
      IF (allocated(mt_decay_source)) deallocate (mt_decay_source)
      IF (allocated(mt_n2o_source  )) deallocate (mt_n2o_source  )
      IF (allocated(mt_levsto       )) deallocate (mt_levsto       )
      IF (allocated(mt_solid)) deallocate (mt_solid, mt_levsto_solid)
      CALL release_mt_substep_workspace()

   END SUBROUTINE mt_riverlake_final

END MODULE MOD_MT_RiverLake
#endif
#endif
