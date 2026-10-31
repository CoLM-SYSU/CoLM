#include <define.h>

#ifdef GridRiverLakeFlow
MODULE MOD_Grid_RiverLakeBifurcation
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   Bifurcation (multi-channel flow) module for grid-based river-lake routing.
!   Computes water exchange through bifurcation pathways with a CaMa-style
!   local-inertial shallow-water update and storage limiters.
!
! Created by CoLM team, April 2026
!-------------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Namelist, only: DEF_USE_LEVEE
   USE MOD_SPMD_Task
   USE MOD_Grid_Reservoir, only: ucat2resv
   USE MOD_Grid_RiverLakeLevee, only: has_levee, levdph, levsto, levee_bashgt, &
      levee_visible_volume_from_stage
   USE MOD_Grid_RiverLakeNetwork, only: floodplain_curve
   IMPLICIT NONE
   PRIVATE

   real(r8), parameter :: RIVERLAKE_DRY_DEPTH = 1.e-5_r8
   real(r8), parameter :: BIFMIN = RIVERLAKE_DRY_DEPTH
   real(r8), parameter :: BIF_MISSING_VALUE = -9999._r8
   integer, parameter :: BIF_RESTART_SIGNATURE_VERSION = 1

   ! ----- State variables -----
   real(r8), allocatable :: pth_veloc     (:,:)  ! velocity (npthlev, npthout_local) [m/s]
   real(r8), allocatable :: pth_momen     (:,:)  ! momentum (npthlev, npthout_local) [m^2/s]
   real(r8), allocatable :: bif_hflux_lev (:,:)  ! effective volume flux per pathway layer [m^3/s]
   real(r8), allocatable :: bif_hflux_sum (:)    ! net volume flux per ucat [m^3/s]
   logical, allocatable :: bif_path_active (:)    ! true when pathway was eligible this call
   ! Per-ucat layer-2+ ("overland" / above-bankfull) bif water flux;
   ! the water-side accumulator `sum_hflux_riv` lumps every layer onto the
   ! visible side, so the caller needs this to find the protected-side part.
   real(r8), allocatable :: bif_lev_hflux_sum(:) ! net layer-2+ vol flux per ucat [m^3/s]
   logical, save :: dbg_bif_restart_checked = .false.
   integer, save :: n_dt_mismatch_warn_count = 0
   ! Static downstream-pushed path fields (riverbed elevation, levee mask,
   ! reservoir mask) do not change within a routing call.  bifurcation_calc
   ! pushes them only on the first sub-step, selected by its explicit
   ! first_substep argument (tied to the globally synchronized substep
   ! counter), so the push sequence is collective-matched by construction and
   ! cannot depend on stale module state.

   ! ----- Persistent scratch buffers (allocated once in bifurcation_init,
   !       reused every bifurcation_calc call to avoid alloc/dealloc churn).
   !       TARGET attribute lets bifurcation_calc alias them via pointers.  -----
   real(r8), allocatable, target, save :: wdsrf_dn_pth_buf     (:)
   real(r8), allocatable, target, save :: wdsrf_prev_dn_pth_buf (:)
   real(r8), allocatable, target, save :: rivelv_dn_pth_buf    (:)
   real(r8), allocatable, target, save :: protected_wdsrf_ucat_buf  (:)
   real(r8), allocatable, target, save :: protected_wdsrf_dn_pth_buf (:)
   real(r8), allocatable, target, save :: has_levee_dn_pth_buf (:)
   real(r8), allocatable, target, save :: storage_dn_pth_buf   (:)
   real(r8), allocatable, target, save :: pth_hflux_total_buf  (:)
   real(r8), allocatable, target, save :: pth_limiter_rate_buf (:)
   real(r8), allocatable, target, save :: bif_influx_buf       (:)
      real(r8), allocatable, target, save :: has_levee_r8_buf     (:)
      real(r8), allocatable, target, save :: storage_ucat_buf     (:)
      real(r8), allocatable, target, save :: visible_storage_ucat_buf   (:)
      real(r8), allocatable, target, save :: protected_storage_ucat_buf (:)
      real(r8), allocatable, target, save :: visible_storage_dn_pth_buf   (:)
      real(r8), allocatable, target, save :: protected_storage_dn_pth_buf (:)
      real(r8), allocatable, target, save :: limiter_outgoing_buf (:)
      real(r8), allocatable, target, save :: limiter_out_rate_buf (:)
      real(r8), allocatable, target, save :: protected_outgoing_buf (:)
      real(r8), allocatable, target, save :: protected_out_rate_buf (:)
      real(r8), allocatable, target, save :: layer_limiter_rate_buf(:,:)

   ! ----- cross-rank limiter buffers -----
   real(r8), allocatable, target, save :: neg_pth_buf              (:)
   real(r8), allocatable, target, save :: bif_neg_recv_buf         (:)
   real(r8), allocatable, target, save :: limiter_out_rate_pth_buf (:)
   real(r8), allocatable, target, save :: protected_neg_pth_buf     (:)
   real(r8), allocatable, target, save :: protected_out_recv_buf    (:)
   real(r8), allocatable, target, save :: protected_out_rate_pth_buf(:)

   ! ----- per-pathway layer-2+ flux + remote influx buffers -----
   real(r8), allocatable, target, save :: pth_lev_hflux_total_buf  (:)
   real(r8), allocatable, target, save :: bif_lev_influx_buf       (:)

   ! ----- cross-river-system sync buffers -----
   ! ucatfilter materialised as real8 so we can push it to pathway space
   ! and let source-rank pathways see whether their downstream cell is
   ! active in the current adaptive sub-step. Real-valued mask follows the
   ! same pattern as has_levee_r8.
   real(r8), allocatable, target, save :: ucatfilter_r8_buf        (:)
   real(r8), allocatable, target, save :: ucatfilter_dn_pth_buf    (:)
   real(r8), allocatable, target, save :: is_resv_r8_buf           (:)
   real(r8), allocatable, target, save :: is_resv_dn_pth_buf       (:)

   ! Local unit catchments that take part in bifurcation on this worker, ascending:
   ! endpoints of local pathways, sources of push_bif_dn2pth and receivers of
   ! push_bif_influx.  Built once in bifurcation_init (the network is fixed).
   ! Reason: bifurcation_calc used to reset and fill about a dozen numucat-sized
   ! buffers and evaluate the storage of every unit catchment on every routing
   ! substep, although only these cells are ever read (17242 pathways on the
   ! 15min network against 252383 unit catchments); this was ~15% of routing CPU
   ! with bifurcation on.  It now does that work for bif_ucat only.
   ! Results are unchanged, bit for bit: the per-cell buffers are only read at
   ! pathway endpoints (i_up, i_dn), at the cells whose values are pushed, and
   ! in the per-cell loops, which run over the same cells in the same order;
   ! values at other cells were written but never used.  bif_hflux_sum and
   ! bif_lev_hflux_sum, which are read by the routing for every cell, are still
   ! reset in full, and the remote inflow they lose is 0 outside bif_ucat.
   integer, allocatable, save :: bif_ucat (:)

   PUBLIC :: bifurcation_init
   PUBLIC :: bifurcation_calc
   PUBLIC :: read_bifurcation_restart
   PUBLIC :: write_bifurcation_restart
   PUBLIC :: bifurcation_final
   PUBLIC :: bif_hflux_sum
   PUBLIC :: bif_hflux_lev
   PUBLIC :: bif_lev_hflux_sum
   PUBLIC :: bif_path_active

CONTAINS

   ! =========================================================================
   SUBROUTINE bifurcation_init ()
   ! =========================================================================
   !
   ! Allocate and initialize bifurcation state arrays.
   !
   ! =========================================================================

   USE MOD_Grid_RiverLakeNetwork, only: numucat, npthout_local, npthlev_bif, &
      pth_upst_local, pth_down_local, push_bif_dn2pth, push_bif_influx, cache_bif_influx

   IMPLICIT NONE

   logical, allocatable :: is_bif_ucat (:)
   integer :: ipth, i, j
#ifdef USEMPI
   integer :: iworker
#endif

      ! bifurcation_init is contractually one-shot per run, but the
      ! allocates below will blow up with a runtime error if the module is
      ! ever re-initialised (e.g. during a future restart-in-place refactor).
      ! A single allocated() check on the first state array is enough to
      ! make re-entry a no-op instead of a hard stop, since all buffers
      ! below are allocated as a group.
      IF (allocated(pth_veloc)) RETURN

      IF (.not. p_is_worker) THEN
         CALL allocate_bifurcation_arrays (0, 0, 0)
         allocate (bif_ucat (0))
         RETURN
      ENDIF

      CALL allocate_bifurcation_arrays (numucat, npthlev_bif, npthout_local)

      ! cells taking part in bifurcation on this worker (see bif_ucat)
      allocate (is_bif_ucat (numucat))
      is_bif_ucat(:) = .false.
      DO ipth = 1, npthout_local
         i = pth_upst_local(ipth)
         IF (i >= 1 .and. i <= numucat) is_bif_ucat(i) = .true.
         i = pth_down_local(ipth)
         IF (i >= 1 .and. i <= numucat) is_bif_ucat(i) = .true.
      ENDDO
      IF (push_bif_dn2pth%nself > 0) is_bif_ucat(push_bif_dn2pth%self_from) = .true.
#ifdef USEMPI
      IF (allocated(push_bif_dn2pth%n_to_other)) THEN
         DO iworker = 0, p_np_worker-1
            IF (push_bif_dn2pth%n_to_other(iworker) > 0) &
               is_bif_ucat(push_bif_dn2pth%to_other(iworker)%val) = .true.
         ENDDO
      ENDIF
#endif
      IF (allocated(cache_bif_influx%ptr)) THEN
         DO j = 1, min(numucat, size(cache_bif_influx%ptr)-1)
            IF (cache_bif_influx%ptr(j+1) > cache_bif_influx%ptr(j)) is_bif_ucat(j) = .true.
         ENDDO
      ENDIF
      bif_ucat = pack((/(i, i = 1, numucat)/), is_bif_ucat)
      deallocate (is_bif_ucat)

   END SUBROUTINE bifurcation_init



   SUBROUTINE bif_path_levee_sides (i_up, i_dn, ipth, nucat, upstream_has_levee, downstream_has_levee)
      IMPLICIT NONE
      integer, intent(in) :: i_up, i_dn, ipth, nucat
      logical, intent(out) :: upstream_has_levee, downstream_has_levee

      upstream_has_levee = .false.
      downstream_has_levee = .false.
      IF (.not. DEF_USE_LEVEE) RETURN

      IF (allocated(has_levee)) THEN
         IF (i_up > 0 .and. i_up <= size(has_levee)) upstream_has_levee = has_levee(i_up)
      ENDIF
      IF (i_dn > 0 .and. i_dn <= nucat) THEN
         IF (allocated(has_levee)) THEN
            IF (i_dn <= size(has_levee)) downstream_has_levee = has_levee(i_dn)
         ENDIF
         RETURN
      ENDIF

      IF (allocated(has_levee_dn_pth_buf)) THEN
         IF (ipth > 0 .and. ipth <= size(has_levee_dn_pth_buf)) &
            downstream_has_levee = has_levee_dn_pth_buf(ipth) > 0.5_r8
      ENDIF
   END SUBROUTINE bif_path_levee_sides


   SUBROUTINE allocate_bifurcation_arrays (nucat, nlev, npth)
   ! =========================================================================
   !
   ! Allocate all persistent bifurcation buffers. Zero-length dimensions are
   ! intentional on non-worker ranks so module state is always allocated safely.
   !
   ! =========================================================================

   IMPLICIT NONE

   integer, intent(in) :: nucat
   integer, intent(in) :: nlev
   integer, intent(in) :: npth

      allocate (pth_veloc     (nlev, npth))
      allocate (pth_momen     (nlev, npth))
      allocate (bif_hflux_lev (nlev, npth))
      allocate (bif_hflux_sum (nucat))
      allocate (bif_path_active (npth))
      allocate (bif_lev_hflux_sum (nucat))

      pth_veloc     (:,:) = 0._r8
      pth_momen     (:,:) = 0._r8
      bif_hflux_lev (:,:) = 0._r8
      bif_hflux_sum (:)   = 0._r8
      bif_path_active(:)  = .false.
      bif_lev_hflux_sum(:) = 0._r8

      ! Persistent scratch buffers used by bifurcation_calc. Sized once here
      ! using the network-fixed dimensions; zero-length is valid.
      allocate (wdsrf_dn_pth_buf     (npth))
      allocate (wdsrf_prev_dn_pth_buf (npth))
      allocate (rivelv_dn_pth_buf    (npth))
      allocate (protected_wdsrf_ucat_buf   (nucat))
      allocate (protected_wdsrf_dn_pth_buf (npth))
      allocate (has_levee_dn_pth_buf (npth))
      allocate (storage_dn_pth_buf   (npth))
      allocate (pth_hflux_total_buf  (npth))
      allocate (pth_limiter_rate_buf (npth))
         allocate (bif_influx_buf       (nucat))
         allocate (has_levee_r8_buf     (nucat))
         allocate (storage_ucat_buf     (nucat))
         allocate (visible_storage_ucat_buf   (nucat))
         allocate (protected_storage_ucat_buf (nucat))
         allocate (visible_storage_dn_pth_buf   (npth))
         allocate (protected_storage_dn_pth_buf (npth))
         allocate (limiter_outgoing_buf (nucat))
         allocate (limiter_out_rate_buf (nucat))
         allocate (protected_outgoing_buf (nucat))
         allocate (protected_out_rate_buf (nucat))
         allocate (layer_limiter_rate_buf(nlev, npth))

      allocate (neg_pth_buf              (npth))
      allocate (bif_neg_recv_buf         (nucat))
      allocate (limiter_out_rate_pth_buf (npth))
      allocate (protected_neg_pth_buf     (npth))
      allocate (protected_out_recv_buf    (nucat))
      allocate (protected_out_rate_pth_buf(npth))

      allocate (ucatfilter_r8_buf        (nucat))
      allocate (ucatfilter_dn_pth_buf    (npth))
      allocate (is_resv_r8_buf           (nucat))
      allocate (is_resv_dn_pth_buf       (npth))

      allocate (pth_lev_hflux_total_buf  (npth))
      allocate (bif_lev_influx_buf       (nucat))

   END SUBROUTINE allocate_bifurcation_arrays


   ! =========================================================================
   SUBROUTINE bifurcation_calc (wdsrf_ucat, wdsrf_ucat_prev, volwater_ucat_in, volwater_ucat_valid_in, &
      volresv, is_built_resv, dt_all, irivsys, ucatfilter, normal_outgoing_rate, first_substep)
   ! =========================================================================
   !
   ! Compute bifurcation fluxes for one sub-timestep. The net volume flux per
   ! unit catchment is accumulated into bif_hflux_sum(:).
   !
   ! =========================================================================

   USE MOD_Const_Physical, only: grav
   USE MOD_WorkerPushData
   USE MOD_Grid_RiverLakePush, only: riverlake_push_field_type, riverlake_push_batch
   USE MOD_Grid_RiverLakeNetwork, only: &
      numucat, npthout_local, npthlev_bif, &
      pth_upst_local, pth_down_local, &
      pth_dst, pth_elv, pth_wth, pth_man, &
      push_bif_dn2pth, push_bif_influx, cache_bif_dn2pth, cache_bif_influx, &
      topo_rivelv

   IMPLICIT NONE

   real(r8), intent(in), target :: wdsrf_ucat (:)   ! water depth above riverbed [m]
   real(r8), intent(in), target :: wdsrf_ucat_prev (:) ! depth used by previous BIF substep [m]
   real(r8), intent(in) :: volwater_ucat_in(:) ! TimeVars-owned routing volume [m^3]
   logical,  intent(in) :: volwater_ucat_valid_in
   real(r8), intent(in) :: volresv    (:)   ! reservoir water volume [m^3]
   logical,  intent(in) :: is_built_resv (:) ! reservoir-built mask
   real(r8), intent(in) :: dt_all     (:)   ! timestep per ucat [s]
   integer,  intent(in) :: irivsys    (:)   ! river system id per ucat
   logical,  intent(in) :: ucatfilter (:)   ! active ucat filter
   real(r8), intent(in) :: normal_outgoing_rate(:) ! ordinary routing gross donor outflow [m3/s]
   logical,  intent(in) :: first_substep          ! first sub-step of this routing call

   ! Pointer aliases to persistent module-level buffers (no per-call alloc).
   real(r8), pointer :: wdsrf_dn_pth     (:) => null()
   real(r8), pointer :: wdsrf_prev_dn_pth (:) => null()
   real(r8), pointer :: rivelv_dn_pth    (:) => null()
   real(r8), pointer :: protected_wdsrf_ucat  (:) => null()
   real(r8), pointer :: protected_wdsrf_dn_pth (:) => null()
   real(r8), pointer :: has_levee_dn_pth (:) => null()
   real(r8), pointer :: storage_dn_pth   (:) => null()
   real(r8), pointer :: storage_ucat     (:) => null()
   real(r8), pointer :: pth_hflux_total  (:) => null()
      real(r8), pointer :: bif_influx       (:) => null()
      real(r8), pointer :: has_levee_r8     (:) => null()
      real(r8), pointer :: limiter_outgoing (:) => null()
      real(r8), pointer :: limiter_out_rate (:) => null()
      real(r8), pointer :: protected_outgoing (:) => null()
      real(r8), pointer :: protected_out_rate (:) => null()
      real(r8), pointer :: visible_storage_ucat   (:) => null()
      real(r8), pointer :: protected_storage_ucat (:) => null()
      real(r8), pointer :: visible_storage_dn_pth   (:) => null()
      real(r8), pointer :: protected_storage_dn_pth (:) => null()
      real(r8), pointer :: layer_limiter_rate(:,:) => null()
   real(r8), pointer :: pth_limiter_rate (:) => null()
   real(r8), pointer :: neg_pth              (:) => null()
   real(r8), pointer :: bif_neg_recv         (:) => null()
   real(r8), pointer :: limiter_out_rate_pth (:) => null()
   real(r8), pointer :: protected_neg_pth     (:) => null()
   real(r8), pointer :: protected_out_recv    (:) => null()
   real(r8), pointer :: protected_out_rate_pth(:) => null()
   real(r8), pointer :: ucatfilter_r8        (:) => null()
   real(r8), pointer :: ucatfilter_dn_pth    (:) => null()
   real(r8), pointer :: is_resv_r8           (:) => null()
   real(r8), pointer :: is_resv_dn_pth       (:) => null()
   real(r8), pointer :: pth_lev_hflux_total  (:) => null()
   real(r8), pointer :: bif_lev_influx       (:) => null()

   real(r8) :: rivelv_up, rivelv_dn, wdsrf_up, wdsrf_dn, wdsrf_up_eff, wdsrf_dn_eff
   real(r8) :: wdsrf_up_prev, wdsrf_dn_prev
   real(r8) :: zsurf_up, zsurf_dn, slope_lev
   real(r8) :: height_up, height_dn, height_up_prev, height_dn_prev
   real(r8) :: h_face, h_face_current, h_face_previous
   real(r8) :: hflux_lev, mflux_lev, zgrad_lev  ! flux/pressure terms for this layer
   real(r8) :: width_pth, pth_are
   real(r8) :: friction, momen_trial, veloc_trial
      real(r8) :: dt, dt_dn, dt_cell, storage_up, storage_dn, storage_ref, rate
      real(r8) :: donor_storage, layer_transfer
      real(r8) :: normal_outflow, bif_outflow, remaining_capacity
      logical  :: upstream_has_levee, downstream_has_levee
      logical  :: use_protected_fields
      type(riverlake_push_field_type) :: dynamic_state_fields(7)
      type(riverlake_push_field_type) :: donor_reverse_fields(2)
      type(riverlake_push_field_type) :: donor_rate_fields(2)
      type(riverlake_push_field_type) :: final_influx_fields(2)
      integer :: n_dynamic_state_fields, n_coupled_fields, n_final_influx_fields
   integer  :: ipth, ilev, i_up, i_dn, i_ucat, iu
   integer  :: n_dt_mismatch_skip, n_dt_mismatch_skip_glb

      IF (.not. p_is_worker) RETURN

      ! Contract: bifurcation_init must have been called before bifurcation_calc.
      ! Guarding so a missing init shows up as an explicit message instead of a
      ! silent segfault on the zero-assign below.
      IF (.not. allocated(bif_hflux_sum) .or. .not. allocated(bif_hflux_lev)) THEN
         write(*,'(A,I0)') 'ERROR: bifurcation_calc called before bifurcation_init on glb=', p_iam_glb
         call flush(6)
         CALL CoLM_stop ()
      ENDIF

      n_dt_mismatch_skip = 0
      n_dt_mismatch_skip_glb = 0
      use_protected_fields = DEF_USE_LEVEE .and. npthlev_bif >= 2
      n_final_influx_fields = merge(2, 1, npthlev_bif >= 2)

      ! Reset net flux accumulator
      bif_hflux_sum(:) = 0._r8
      bif_hflux_lev(:,:) = 0._r8
      IF (allocated(bif_path_active)) bif_path_active(:) = .false.
      IF (allocated(bif_lev_hflux_sum)) bif_lev_hflux_sum(:) = 0._r8

      ! Alias persistent module-level buffers (allocated once in
      ! bifurcation_init). Zero-length buffers are still valid targets.
      wdsrf_dn_pth     => wdsrf_dn_pth_buf
      wdsrf_prev_dn_pth => wdsrf_prev_dn_pth_buf
      rivelv_dn_pth    => rivelv_dn_pth_buf
      protected_wdsrf_ucat   => protected_wdsrf_ucat_buf
      protected_wdsrf_dn_pth => protected_wdsrf_dn_pth_buf
      has_levee_dn_pth => has_levee_dn_pth_buf
      storage_dn_pth   => storage_dn_pth_buf
      pth_hflux_total  => pth_hflux_total_buf
      pth_limiter_rate => pth_limiter_rate_buf
         bif_influx       => bif_influx_buf
         has_levee_r8     => has_levee_r8_buf
         storage_ucat     => storage_ucat_buf
         visible_storage_ucat   => visible_storage_ucat_buf
         protected_storage_ucat => protected_storage_ucat_buf
         visible_storage_dn_pth   => visible_storage_dn_pth_buf
         protected_storage_dn_pth => protected_storage_dn_pth_buf
         limiter_outgoing => limiter_outgoing_buf
         limiter_out_rate => limiter_out_rate_buf
         protected_outgoing => protected_outgoing_buf
         protected_out_rate => protected_out_rate_buf
         layer_limiter_rate => layer_limiter_rate_buf
      neg_pth              => neg_pth_buf
      bif_neg_recv         => bif_neg_recv_buf
      limiter_out_rate_pth => limiter_out_rate_pth_buf
      protected_neg_pth    => protected_neg_pth_buf
      protected_out_recv   => protected_out_recv_buf
      protected_out_rate_pth => protected_out_rate_pth_buf
      ucatfilter_r8        => ucatfilter_r8_buf
      ucatfilter_dn_pth    => ucatfilter_dn_pth_buf
      is_resv_r8           => is_resv_r8_buf
      is_resv_dn_pth       => is_resv_dn_pth_buf
      pth_lev_hflux_total  => pth_lev_hflux_total_buf
      bif_lev_influx       => bif_lev_influx_buf

      wdsrf_dn_pth     (:) = 0._r8
      wdsrf_prev_dn_pth(:) = BIF_MISSING_VALUE
      IF (first_substep) rivelv_dn_pth(:) = 0._r8
      protected_wdsrf_dn_pth(:) = BIF_MISSING_VALUE
      IF (first_substep) has_levee_dn_pth(:) = 0._r8
         storage_dn_pth   (:) = 0._r8
         visible_storage_dn_pth(:) = 0._r8
         protected_storage_dn_pth(:) = 0._r8
         pth_hflux_total(:) = 0._r8
         pth_limiter_rate(:) = 1._r8
         layer_limiter_rate(:,:) = 1._r8
      pth_lev_hflux_total(:) = 0._r8
      neg_pth             (:) = 0._r8
      limiter_out_rate_pth(:) = 1._r8
      protected_neg_pth(:) = 0._r8
      protected_out_rate_pth(:) = 1._r8
      ucatfilter_dn_pth(:) = 0._r8
      IF (first_substep) is_resv_dn_pth(:) = 0._r8

      ! Per-cell buffers are only read at bifurcation cells, so they are reset
      ! and filled for bif_ucat only (see its declaration; values identical).
      ! The receive buffers bif_neg_recv / protected_out_recv / bif_lev_influx
      ! are overwritten in full by their pushes.
      ! storage_ucat is fully overwritten for these cells below.
      ! has_levee_r8 is referenced only under DEF_USE_LEVEE guards, but zero it
      ! unconditionally so a future refactor that decouples the guards can't
      ! read uninitialised values.  ucatfilter / is_built_resv are materialised
      ! as r8 so the push machinery can ship them to each pathway's
      ! destination cell (the dn_pth views are filled by the Step 2 push).
      DO iu = 1, size(bif_ucat)
         i_ucat = bif_ucat(iu)
         protected_wdsrf_ucat(i_ucat) = 0._r8
         visible_storage_ucat(i_ucat) = 0._r8
         protected_storage_ucat(i_ucat) = 0._r8
         limiter_outgoing(i_ucat) = 0._r8
         limiter_out_rate(i_ucat) = 1._r8
         protected_outgoing(i_ucat) = 0._r8
         protected_out_rate(i_ucat) = 1._r8
         bif_lev_influx(i_ucat) = 0._r8
         has_levee_r8(i_ucat) = 0._r8
         bif_neg_recv(i_ucat) = 0._r8
         protected_out_recv(i_ucat) = 0._r8
         ucatfilter_r8(i_ucat) = merge(1._r8, 0._r8, ucatfilter(i_ucat))
         is_resv_r8   (i_ucat) = merge(1._r8, 0._r8, is_built_resv(i_ucat))
      ENDDO

      IF (DEF_USE_LEVEE) THEN
         IF (allocated(has_levee)) THEN
            DO iu = 1, size(bif_ucat)
               i_ucat = bif_ucat(iu)
               IF (has_levee(i_ucat)) THEN
                  has_levee_r8(i_ucat) = 1._r8
                  protected_wdsrf_ucat(i_ucat) = floodplain_curve(i_ucat)%rivhgt + levdph(i_ucat)
                  ! A dry protected pool has zero depth but its bed is above
                  ! the floodplain datum. Wet levdph already includes this bed.
                  IF (levsto(i_ucat) <= 0._r8) protected_wdsrf_ucat(i_ucat) = &
                     floodplain_curve(i_ucat)%rivhgt + levee_bashgt(i_ucat)
               ELSE
                  protected_wdsrf_ucat(i_ucat) = wdsrf_ucat(i_ucat)
               ENDIF
            ENDDO
         ENDIF
      ENDIF

      DO iu = 1, size(bif_ucat)
         i_ucat = bif_ucat(iu)
         storage_ucat(i_ucat) = available_storage_ucat(i_ucat, wdsrf_ucat(i_ucat), &
            volwater_ucat_in, volwater_ucat_valid_in, volresv, is_built_resv)
         IF (DEF_USE_LEVEE) THEN
            visible_storage_ucat(i_ucat) = available_visible_storage_ucat(i_ucat, wdsrf_ucat(i_ucat), &
               volwater_ucat_in, volwater_ucat_valid_in, volresv, is_built_resv)
            protected_storage_ucat(i_ucat) = available_protected_storage_ucat(i_ucat, is_built_resv)
         ELSE
            visible_storage_ucat(i_ucat) = storage_ucat(i_ucat)
            protected_storage_ucat(i_ucat) = 0._r8
         ENDIF
      ENDDO

      ! ----- Step 2: Get downstream cell state via push objects -----
      ! Dynamic fields share one fixed-order communication phase.  All ranks use
      ! the same field count because DEF_USE_LEVEE and npthlev_bif are global.
      ! Single-layer levee runs still need visible storage for the channel donor
      ! limiter, but omit protected stage/storage fields that no layer can use.
      dynamic_state_fields(1)%send => wdsrf_ucat
      dynamic_state_fields(1)%recv => wdsrf_dn_pth
      dynamic_state_fields(1)%fillvalue = BIF_MISSING_VALUE
      dynamic_state_fields(2)%send => wdsrf_ucat_prev
      dynamic_state_fields(2)%recv => wdsrf_prev_dn_pth
      dynamic_state_fields(2)%fillvalue = BIF_MISSING_VALUE
      dynamic_state_fields(3)%send => storage_ucat
      dynamic_state_fields(3)%recv => storage_dn_pth
      dynamic_state_fields(3)%fillvalue = BIF_MISSING_VALUE
      dynamic_state_fields(4)%send => ucatfilter_r8
      dynamic_state_fields(4)%recv => ucatfilter_dn_pth
      dynamic_state_fields(4)%fillvalue = 0._r8
      n_dynamic_state_fields = 4

      IF (DEF_USE_LEVEE) THEN
         dynamic_state_fields(5)%send => visible_storage_ucat
         dynamic_state_fields(5)%recv => visible_storage_dn_pth
         dynamic_state_fields(5)%fillvalue = 0._r8
         n_dynamic_state_fields = 5
      ENDIF
      IF (use_protected_fields) THEN
         dynamic_state_fields(6)%send => protected_wdsrf_ucat
         dynamic_state_fields(6)%recv => protected_wdsrf_dn_pth
         dynamic_state_fields(6)%fillvalue = BIF_MISSING_VALUE
         dynamic_state_fields(7)%send => protected_storage_ucat
         dynamic_state_fields(7)%recv => protected_storage_dn_pth
         dynamic_state_fields(7)%fillvalue = 0._r8
         n_dynamic_state_fields = 7
      ENDIF

      CALL riverlake_push_batch (push_bif_dn2pth, cache_bif_dn2pth, &
         dynamic_state_fields(1:n_dynamic_state_fields))

      IF (.not. DEF_USE_LEVEE) THEN
         visible_storage_dn_pth(:) = storage_dn_pth(:)
         protected_storage_dn_pth(:) = 0._r8
      ELSEIF (.not. use_protected_fields) THEN
         protected_storage_dn_pth(:) = 0._r8
      ENDIF

      ! Static fields are pushed on the first sub-step only (explicit argument
      ! instead of a module flag, so every worker makes this decision from the
      ! same globally synchronized substep counter).
      IF (first_substep) &
         CALL worker_push_data (push_bif_dn2pth, topo_rivelv, rivelv_dn_pth, fillvalue = BIF_MISSING_VALUE)
      IF (use_protected_fields .and. first_substep) &
         CALL worker_push_data (push_bif_dn2pth, has_levee_r8, has_levee_dn_pth, fillvalue = 0._r8)
      IF (first_substep) &
         CALL worker_push_data (push_bif_dn2pth, is_resv_r8, is_resv_dn_pth, fillvalue = 0._r8)

      ! Static downstream path fields (rivelv/has_levee/is_resv) are now cached
      ! in their persistent path buffers for the remaining sub-steps of this
      ! routing call.
      ! ----- Step 3: local-inertial solver for each pathway / layer -----
      DO ipth = 1, npthout_local

         i_up = pth_upst_local(ipth)

         ! Skip if upstream ucat is not active
         IF (i_up < 1 .or. i_up > numucat) CYCLE
         IF (.not. ucatfilter(i_up)) CYCLE
         IF (i_up <= size(is_built_resv) .and. is_built_resv(i_up)) CYCLE

         ! Skip if downstream cell is outside domain (CaMa: CYCLE when JSEQP<=0)
         IF (wdsrf_dn_pth(ipth) < -9000._r8 .or. rivelv_dn_pth(ipth) < -9000._r8) CYCLE

         rivelv_up = topo_rivelv(i_up)
         wdsrf_up  = wdsrf_ucat(i_up)
         wdsrf_up_prev = wdsrf_ucat_prev(i_up)
         rivelv_dn = rivelv_dn_pth(ipth)
         wdsrf_dn  = wdsrf_dn_pth(ipth)
         wdsrf_dn_prev = wdsrf_prev_dn_pth(ipth)
         dt        = dt_all(irivsys(i_up))

         pth_hflux_total(ipth) = 0._r8

         i_dn = pth_down_local(ipth)

         ! skip pathway when destination cell is in a river system
         ! that has exhausted its adaptive sub-step (dt = 0). Applying the
         ! flux would update bif_hflux_sum(i_up) on the source side, but
         ! the destination contribution is dropped by WHERE(ucatfilter) in
         ! MOD_Grid_RiverLakeFlow, producing a cross-system conservation
         ! violation. Deferring the transfer to the next top-level routing
         ! step is safe: pathway momentum state is carried forward
         ! unchanged, and the two systems resynchronise at the next outer
         ! routing step.
         IF (i_dn > 0 .and. i_dn <= numucat) THEN
            IF (.not. ucatfilter(i_dn)) CYCLE
            IF (i_dn <= size(is_built_resv) .and. is_built_resv(i_dn)) CYCLE
            dt_dn = dt_all(irivsys(i_dn))
         ELSE
            IF (ucatfilter_dn_pth(ipth) < 0.5_r8) CYCLE
            IF (is_resv_dn_pth(ipth) > 0.5_r8) CYCLE
            ! Flow synchronizes the dt of each group of river systems connected
            ! by bifurcation pathways after all local stability constraints, and
            ! both ends of a pathway are in the same group.  The remote active
            ! mask is therefore sufficient; sending the identical destination dt
            ! was a redundant MPI phase.
            dt_dn = dt
         ENDIF

         IF (dt <= 0._r8 .or. dt_dn <= 0._r8) CYCLE
         IF (abs(dt - dt_dn) > max(1.e-9_r8, 1.e-12_r8 * max(abs(dt), abs(dt_dn)))) THEN
            n_dt_mismatch_skip = n_dt_mismatch_skip + 1
            CYCLE
         ENDIF

         ! Zero-length pathways have no physical length scale for the
         ! momentum update; clear all layer states and skip this pathway.
         IF (pth_dst(ipth) <= 0._r8) THEN
            pth_momen(:, ipth) = 0._r8
            pth_veloc(:, ipth) = 0._r8
            CYCLE
         ENDIF

         IF (allocated(bif_path_active)) bif_path_active(ipth) = .true.
         CALL bif_path_levee_sides(i_up, i_dn, ipth, numucat, upstream_has_levee, downstream_has_levee)

         DO ilev = 1, npthlev_bif

            width_pth = pth_wth(ilev, ipth)
            IF (width_pth <= 0._r8) THEN
               pth_momen(ilev, ipth) = 0._r8
               pth_veloc(ilev, ipth) = 0._r8
               CYCLE
            ENDIF

            IF (ilev == 1) THEN
               wdsrf_up_eff = wdsrf_up
            ELSE
               IF (DEF_USE_LEVEE .and. upstream_has_levee) THEN
                  wdsrf_up_eff = protected_wdsrf_ucat(i_up)
               ELSE
                  wdsrf_up_eff = wdsrf_up
               ENDIF
            ENDIF

            IF (ilev == 1) THEN
               wdsrf_dn_eff = wdsrf_dn
            ELSE
               IF (i_dn > 0 .and. i_dn <= numucat) THEN
                  IF (DEF_USE_LEVEE .and. downstream_has_levee) THEN
                     wdsrf_dn_eff = protected_wdsrf_ucat(i_dn)
                  ELSE
                     wdsrf_dn_eff = wdsrf_dn
                  ENDIF
               ELSEIF (DEF_USE_LEVEE .and. downstream_has_levee) THEN
                  IF (protected_wdsrf_dn_pth(ipth) < -9000._r8) CYCLE
                  wdsrf_dn_eff = protected_wdsrf_dn_pth(ipth)
               ELSE
                  wdsrf_dn_eff = wdsrf_dn
               ENDIF
            ENDIF

            ! CaMa PTH_ELV semantics: path layer elevation is the hydraulic
            ! sill; h_face below takes max(surface_up, surface_down) - PTH_ELV.
            height_up = max(0._r8, wdsrf_up_eff + rivelv_up - pth_elv(ilev, ipth))
            height_dn = max(0._r8, wdsrf_dn_eff + rivelv_dn - pth_elv(ilev, ipth))

            ! Skip if both ends are dry
            IF (height_up < BIFMIN .and. height_dn < BIFMIN) THEN
               hflux_lev = 0._r8
               pth_momen(ilev, ipth) = 0._r8
               pth_veloc(ilev, ipth) = 0._r8
               CYCLE
            ENDIF

            ! Local-inertial pathway solver.  The old Godunov momentum
            ! update used the absolute HLL pressure flux (0.5*g*h**2)
            ! directly, which accelerates a dry-still link even when the two
            ! water surfaces are level.  Drive pathway momentum by water-surface
            ! slope instead; equal H_up/H_dn then leaves only frictional decay.
            ! CaMa-Flood levee alignment: channel flow uses the river-side
            ! surface, while overland flow uses the protected-side surface
            ! wherever a levee exists.  Use the same effective surface for
            ! both flow depth and hydraulic slope.
            zsurf_up = wdsrf_up_eff + rivelv_up
            zsurf_dn = wdsrf_dn_eff + rivelv_dn
            slope_lev = (zsurf_dn - zsurf_up) / pth_dst(ipth)
            slope_lev = max(-0.005_r8, min(0.005_r8, slope_lev))
            h_face_current = max(height_up, height_dn)
            IF (DEF_USE_LEVEE) THEN
               IF (ilev == 1) THEN
                  ! CaMa-Flood levee scheme: stabilize the in-channel local-
                  ! inertial update with the geometric mean of current and
                  ! previous-substep flow depth. If either side was dry, CaMa
                  ! falls back to the current depth.
                  height_up_prev = max(0._r8, wdsrf_up_prev + rivelv_up - pth_elv(ilev, ipth))
                  height_dn_prev = max(0._r8, wdsrf_dn_prev + rivelv_dn - pth_elv(ilev, ipth))
                  h_face_previous = max(height_up_prev, height_dn_prev)
                  h_face = sqrt(h_face_current * h_face_previous)
                  IF (h_face <= 0._r8) h_face = h_face_current
               ELSE
                  ! Above-bankfull layers use the protected-side surface and
                  ! stay explicit because the simplified levee scheme has no
                  ! matching previous protected-depth field.
                  h_face = h_face_current
               ENDIF
            ELSE
               ! Standard CaMa BIF applies the semi-implicit depth to every
               ! layer.  Its 0.01 m previous-depth floor preserves finite
               ! conductance when a previously dry pathway becomes wet.
               height_up_prev = max(0._r8, wdsrf_up_prev + rivelv_up - pth_elv(ilev, ipth))
               height_dn_prev = max(0._r8, wdsrf_dn_prev + rivelv_dn - pth_elv(ilev, ipth))
               h_face_previous = max(height_up_prev, height_dn_prev)
               h_face = sqrt(h_face_current * h_face_previous)
               h_face = max(h_face, sqrt(h_face_current * 0.01_r8))
            ENDIF

            IF (h_face <= BIFMIN) THEN
               hflux_lev = 0._r8
               pth_momen(ilev, ipth) = 0._r8
               pth_veloc(ilev, ipth) = 0._r8
               CYCLE
            ENDIF

            pth_are = width_pth * pth_dst(ipth)
            friction = grav * pth_man(ilev)**2 / h_face**(7._r8/3._r8) * abs(pth_momen(ilev, ipth))

            ! zgrad_lev is the single-link hydrostatic balance term: this is
            ! algebraically equivalent to using (mflux_lev - zgrad_lev) in the
            ! momentum equation, but avoids carrying the HLL advective pressure
            ! flux into a link-centered momentum state.  The expression below is
            ! g*h*dH/dx in single-width momentum units [m2/s2].
            mflux_lev = width_pth * grav * h_face * slope_lev * pth_dst(ipth)
            zgrad_lev = 0._r8
            momen_trial = (pth_momen(ilev, ipth) - (mflux_lev - zgrad_lev) / pth_are * dt) &
               / (1._r8 + friction * dt)
            veloc_trial = momen_trial / h_face

            ! Clamp velocity and keep the trial momentum in sync.
            veloc_trial = max(-20._r8, min(20._r8, veloc_trial))
            momen_trial = veloc_trial * h_face

            pth_momen(ilev, ipth) = momen_trial
            pth_veloc(ilev, ipth) = veloc_trial

            hflux_lev = width_pth * veloc_trial * h_face

            ! Accumulate total flux for this pathway [m^3/s]
            pth_hflux_total(ipth) = pth_hflux_total(ipth) + hflux_lev
            bif_hflux_lev(ilev, ipth) = hflux_lev

         ENDDO  ! ilev

         ! ----- Step 4: 5% storage limiter -----
            IF (abs(pth_hflux_total(ipth)) > 0._r8 .and. dt > 0._r8) THEN
               storage_up = storage_ucat(i_up)
               storage_dn = storage_dn_pth(ipth)
               ! CaMa path limiter: one path may move at most 5% of the
               ! smaller endpoint storage in this substep.
               storage_ref = max(min(storage_up, storage_dn), 0._r8)
               pth_limiter_rate(ipth) = bif_limiter_fraction(0.05_r8 * storage_ref, &
                  abs(pth_hflux_total(ipth)) * dt)
            ENDIF

            ! Split-pool no-overdraft limiter for CoLM's levee compartment.
            ! CaMa's original 5% limiter is total storage based; CoLM
            ! additionally has a protected levee pool (`levsto`).  Limit
            ! each layer by the donor compartment only so a dry protected side
            ! cannot be silently clipped by the caller.  CaMa's path-level
            ! endpoint-storage limiter above remains the stability constraint.
            IF (dt > 0._r8) THEN
               DO ilev = 1, npthlev_bif
                  layer_transfer = abs(bif_hflux_lev(ilev, ipth))
                  IF (layer_transfer <= 0._r8) CYCLE

                  IF (bif_hflux_lev(ilev, ipth) >= 0._r8) THEN
                     ! upstream -> downstream: local upstream cell is donor
                     IF (ilev > 1 .and. upstream_has_levee) THEN
                        donor_storage = protected_storage_ucat(i_up)
                        protected_outgoing(i_up) = protected_outgoing(i_up) + layer_transfer
                     ELSE
                        donor_storage = visible_storage_ucat(i_up)
                        limiter_outgoing(i_up) = limiter_outgoing(i_up) + layer_transfer
                     ENDIF
                  ELSE
                     ! downstream -> upstream: destination cell is donor.
                     ! Its outgoing total is gathered through push_bif_influx
                     ! also when it is local (neg_pth / protected_neg_pth), so
                     ! it is summed in global pathway order whatever the
                     ! decomposition (adding local terms here directly made
                     ! the results depend on the number of processes).
                     IF (i_dn > 0 .and. i_dn <= numucat) THEN
                        IF (ilev > 1 .and. downstream_has_levee) THEN
                           donor_storage = protected_storage_ucat(i_dn)
                           protected_neg_pth(ipth) = protected_neg_pth(ipth) + layer_transfer
                        ELSE
                           donor_storage = visible_storage_ucat(i_dn)
                           neg_pth(ipth) = neg_pth(ipth) + layer_transfer
                        ENDIF
                     ELSE
                        IF (ilev > 1 .and. downstream_has_levee) THEN
                           donor_storage = protected_storage_dn_pth(ipth)
                           protected_neg_pth(ipth) = protected_neg_pth(ipth) + layer_transfer
                        ELSE
                           donor_storage = visible_storage_dn_pth(ipth)
                           neg_pth(ipth) = neg_pth(ipth) + layer_transfer
                        ENDIF
                     ENDIF
                  ENDIF

                  layer_limiter_rate(ilev, ipth) = min(layer_limiter_rate(ilev, ipth), &
                     bif_limiter_fraction(max(donor_storage, 0._r8), layer_transfer * dt))
               ENDDO
            ENDIF

      ENDDO  ! ipth

      ! Reverse-flow donor totals share one multi-mapping phase.  Runs without
      ! protected overland flow send only the visible field.
      n_coupled_fields = merge(2, 1, use_protected_fields)
      donor_reverse_fields(1)%send => neg_pth
      donor_reverse_fields(1)%recv => bif_neg_recv
      donor_reverse_fields(1)%fillvalue = 0._r8
      IF (use_protected_fields) THEN
         donor_reverse_fields(2)%send => protected_neg_pth
         donor_reverse_fields(2)%recv => protected_out_recv
         donor_reverse_fields(2)%fillvalue = 0._r8
      ENDIF
      CALL riverlake_push_batch (push_bif_influx, cache_bif_influx, &
         donor_reverse_fields(1:n_coupled_fields), mode = 'sum')

      ! Per-cell loops below run over bif_ucat only: elsewhere the received donor
      ! totals and the remote inflow are 0 and the per-cell rates are never read
      ! (see bif_ucat), so the results are identical.
      IF (use_protected_fields) THEN
         DO iu = 1, size(bif_ucat)
            i_ucat = bif_ucat(iu)
            protected_outgoing(i_ucat) = protected_outgoing(i_ucat) + protected_out_recv(i_ucat)
            IF (.not. ucatfilter(i_ucat)) CYCLE
            dt_cell = dt_all(irivsys(i_ucat))
            IF (protected_outgoing(i_ucat) > 0._r8 .and. dt_cell > 0._r8) THEN
               protected_out_rate(i_ucat) = bif_limiter_fraction( &
                  max(protected_storage_ucat(i_ucat), 0._r8), protected_outgoing(i_ucat) * dt_cell)
            ENDIF
         ENDDO
      ENDIF

#ifdef CoLMDEBUG
      ! Diagnostic only: the global skip count just feeds the warning below.
      ! Keep this allreduce inside the CoLM-DEBUG guard so production builds do
      ! not pay a global collective on every routing substep for an unused count.
      n_dt_mismatch_skip_glb = n_dt_mismatch_skip
#ifdef USEMPI
      CALL mpi_allreduce (n_dt_mismatch_skip, n_dt_mismatch_skip_glb, 1, MPI_INTEGER, MPI_SUM, p_comm_worker, p_err)
#endif
      IF (p_iam_worker == p_root .and. n_dt_mismatch_skip_glb > 0 .and. n_dt_mismatch_warn_count < 5) THEN
         write(*,'(A,I0,A)') 'WARNING bifurcation_calc: skipped ', n_dt_mismatch_skip_glb, &
            ' pathway(s) because source/destination routing dt differ; cross-river-system bifurcation may be inactive.'
         n_dt_mismatch_warn_count = n_dt_mismatch_warn_count + 1
      ENDIF
#endif
      ! CaMa-style aggregate donor limiter.  CaMa applies the cell rate to
      ! ordinary routing and BIF together; CoLM's ordinary routing flux has
      ! already fixed dt and is not rescaled here.  For split-pool levee cells,
      ! this limiter covers the visible pool; protected pool has its own
      ! aggregate limiter above.
      DO iu = 1, size(bif_ucat)
         i_ucat = bif_ucat(iu)
         limiter_outgoing(i_ucat) = limiter_outgoing(i_ucat) + bif_neg_recv(i_ucat)
      ENDDO

      DO iu = 1, size(bif_ucat)
         i_ucat = bif_ucat(iu)
         IF (.not. ucatfilter(i_ucat)) CYCLE
         dt_cell = dt_all(irivsys(i_ucat))
         storage_ref = max(visible_storage_ucat(i_ucat), 0._r8)
         normal_outflow = 0._r8
         IF (i_ucat <= size(normal_outgoing_rate)) &
            normal_outflow = max(normal_outgoing_rate(i_ucat), 0._r8)
         bif_outflow = limiter_outgoing(i_ucat)
         IF (bif_outflow > 0._r8 .and. dt_cell > 0._r8) THEN
            remaining_capacity = max(storage_ref / dt_cell - normal_outflow, 0._r8)
            IF (remaining_capacity < bif_outflow) THEN
               limiter_out_rate(i_ucat) = min(1._r8, remaining_capacity / bif_outflow)
            ENDIF
         ENDIF
      ENDDO

      ! Push visible/protected downstream donor rates in one single-mapping
      ! phase. Fillvalue=1 keeps boundary pathways unconstrained.
      donor_rate_fields(1)%send => limiter_out_rate
      donor_rate_fields(1)%recv => limiter_out_rate_pth
      donor_rate_fields(1)%fillvalue = 1._r8
      IF (use_protected_fields) THEN
         donor_rate_fields(2)%send => protected_out_rate
         donor_rate_fields(2)%recv => protected_out_rate_pth
         donor_rate_fields(2)%fillvalue = 1._r8
      ENDIF
      CALL riverlake_push_batch (push_bif_dn2pth, cache_bif_dn2pth, &
         donor_rate_fields(1:n_coupled_fields))

      DO ipth = 1, npthout_local

         IF (allocated(bif_path_active)) THEN
            IF (ipth <= size(bif_path_active)) THEN
               IF (.not. bif_path_active(ipth)) CYCLE
            ENDIF
         ENDIF

         i_up = pth_upst_local(ipth)
         IF (i_up < 1 .or. i_up > numucat) CYCLE
         IF (.not. ucatfilter(i_up)) CYCLE
         IF (i_up <= size(is_built_resv) .and. is_built_resv(i_up)) CYCLE
         IF (wdsrf_dn_pth(ipth) < -9000._r8 .or. rivelv_dn_pth(ipth) < -9000._r8) CYCLE
         ! Zero-length pathways are skipped up in the first DO ipth loop, so
         ! pth_hflux_total(ipth) stays 0 here and this loop's scaling is a
         ! no-op. No need to re-check pth_dst.

         i_dn = pth_down_local(ipth)
         ! mirror the first DO loop's cross-system sync guard here
         ! so we never scale pth_momen / pth_veloc by a rate that the
         ! upstream system's other pathways produced for a pathway whose
         ! destination is not participating in this sub-step.
         IF (i_dn > 0 .and. i_dn <= numucat) THEN
            IF (.not. ucatfilter(i_dn)) CYCLE
            IF (i_dn <= size(is_built_resv) .and. is_built_resv(i_dn)) CYCLE
            dt_dn = dt_all(irivsys(i_dn))
         ELSE
            IF (ucatfilter_dn_pth(ipth) < 0.5_r8) CYCLE
            IF (is_resv_dn_pth(ipth) > 0.5_r8) CYCLE
            dt_dn = dt_all(irivsys(i_up))
         ENDIF
         dt = dt_all(irivsys(i_up))
         IF (dt <= 0._r8 .or. dt_dn <= 0._r8) CYCLE
         IF (abs(dt - dt_dn) > max(1.e-9_r8, 1.e-12_r8 * max(abs(dt), abs(dt_dn)))) CYCLE
            CALL bif_path_levee_sides(i_up, i_dn, ipth, numucat, upstream_has_levee, downstream_has_levee)

            pth_hflux_total(ipth) = 0._r8
            DO ilev = 1, npthlev_bif
               rate = min(pth_limiter_rate(ipth), layer_limiter_rate(ilev, ipth))
               IF (DEF_USE_LEVEE .and. ilev > 1) THEN
                  IF (bif_hflux_lev(ilev, ipth) > 0._r8 .and. upstream_has_levee) THEN
                     rate = min(rate, protected_out_rate(i_up))
                  ELSEIF (bif_hflux_lev(ilev, ipth) < 0._r8 .and. downstream_has_levee) THEN
                     IF (i_dn > 0 .and. i_dn <= numucat) THEN
                        rate = min(rate, protected_out_rate(i_dn))
                     ELSE
                        rate = min(rate, protected_out_rate_pth(ipth))
                     ENDIF
                  ELSEIF (bif_hflux_lev(ilev, ipth) > 0._r8) THEN
                     rate = min(rate, limiter_out_rate(i_up))
                  ELSEIF (bif_hflux_lev(ilev, ipth) < 0._r8) THEN
                     IF (i_dn > 0 .and. i_dn <= numucat) THEN
                        rate = min(rate, limiter_out_rate(i_dn))
                     ELSE
                        rate = min(rate, limiter_out_rate_pth(ipth))
                     ENDIF
                  ENDIF
               ELSEIF (bif_hflux_lev(ilev, ipth) > 0._r8) THEN
                  rate = min(rate, limiter_out_rate(i_up))
               ELSEIF (bif_hflux_lev(ilev, ipth) < 0._r8) THEN
                  IF (i_dn > 0 .and. i_dn <= numucat) THEN
                     rate = min(rate, limiter_out_rate(i_dn))
                  ELSE
                     rate = min(rate, limiter_out_rate_pth(ipth))
                  ENDIF
               ENDIF
               IF (rate < 1._r8) THEN
                  bif_hflux_lev(ilev, ipth) = bif_hflux_lev(ilev, ipth) * rate
                  pth_momen(ilev, ipth) = pth_momen(ilev, ipth) * rate
                  pth_veloc(ilev, ipth) = pth_veloc(ilev, ipth) * rate
               ENDIF
               pth_hflux_total(ipth) = pth_hflux_total(ipth) + bif_hflux_lev(ilev, ipth)
            ENDDO

            ! Donor limits can remove cancellation between opposed layers.
            ! Reapply the path cap to the final net flux, scaling state together.
            IF (abs(pth_hflux_total(ipth)) > 0._r8) THEN
               storage_ref = max(min(storage_ucat(i_up), storage_dn_pth(ipth)), 0._r8)
               rate = bif_limiter_fraction(0.05_r8 * storage_ref, abs(pth_hflux_total(ipth)) * dt)
               IF (rate < 1._r8) THEN
                  bif_hflux_lev(:, ipth) = bif_hflux_lev(:, ipth) * rate
                  pth_momen(:, ipth) = pth_momen(:, ipth) * rate
                  pth_veloc(:, ipth) = pth_veloc(:, ipth) * rate
                  pth_hflux_total(ipth) = sum(bif_hflux_lev(:, ipth))
               ENDIF
            ENDIF

            ! ----- Step 5: Accumulate to upstream ucat (local) -----
         bif_hflux_sum(i_up) = bif_hflux_sum(i_up) + pth_hflux_total(ipth)

         IF (npthlev_bif >= 2) THEN
            pth_lev_hflux_total(ipth) = sum(bif_hflux_lev(2:npthlev_bif, ipth))
            bif_lev_hflux_sum(i_up) = bif_lev_hflux_sum(i_up) + pth_lev_hflux_total(ipth)
         ENDIF

      ENDDO  ! ipth

      ! ----- Step 6: Scatter flux to the downstream ucats -----
      ! Every downstream ucat, local or remote, receives its inflow through this
      ! push, summed in global pathway order.  Local downstream ucats used to be
      ! updated directly in step 5 and their pathways subtracted again from the
      ! push result ((a+b)-a /= b in floating point), which made the results
      ! depend on the number of processes.
      IF (numucat > 0) THEN
         bif_influx(:) = 0._r8
         bif_lev_influx(:) = 0._r8
      ENDIF
      final_influx_fields(1)%send => pth_hflux_total
      final_influx_fields(1)%recv => bif_influx
      final_influx_fields(1)%fillvalue = 0._r8
      IF (npthlev_bif >= 2) THEN
         final_influx_fields(2)%send => pth_lev_hflux_total
         final_influx_fields(2)%recv => bif_lev_influx
         final_influx_fields(2)%fillvalue = 0._r8
      ENDIF
      CALL riverlake_push_batch (push_bif_influx, cache_bif_influx, &
         final_influx_fields(1:n_final_influx_fields), mode = 'sum')

      ! ----- Step 7: Apply the inflow to bif_hflux_sum (it is 0 outside bif_ucat)
      DO iu = 1, size(bif_ucat)
         i_ucat = bif_ucat(iu)
         bif_hflux_sum(i_ucat) = bif_hflux_sum(i_ucat) - bif_influx(i_ucat)
      ENDDO
      IF (npthlev_bif >= 2) THEN
         DO iu = 1, size(bif_ucat)
            i_ucat = bif_ucat(iu)
            bif_lev_hflux_sum(i_ucat) = bif_lev_hflux_sum(i_ucat) - bif_lev_influx(i_ucat)
         ENDDO
      ENDIF

      ! Disassociate pointer aliases (targets remain allocated, owned by module)
      nullify (wdsrf_dn_pth, wdsrf_prev_dn_pth, rivelv_dn_pth, &
                  protected_wdsrf_ucat, protected_wdsrf_dn_pth, has_levee_dn_pth, &
                  storage_dn_pth, visible_storage_ucat, protected_storage_ucat, &
                  visible_storage_dn_pth, protected_storage_dn_pth, &
                  pth_hflux_total, pth_limiter_rate, bif_influx, &
                  has_levee_r8, storage_ucat, limiter_outgoing, limiter_out_rate, &
                  protected_outgoing, protected_out_rate, &
                  layer_limiter_rate, &
               neg_pth, bif_neg_recv, &
               limiter_out_rate_pth, protected_neg_pth, &
               protected_out_recv, protected_out_rate_pth, &
               ucatfilter_r8, ucatfilter_dn_pth, is_resv_r8, is_resv_dn_pth, &
               pth_lev_hflux_total, bif_lev_influx)

   END SUBROUTINE bifurcation_calc


   ! =========================================================================
   PURE FUNCTION bif_limiter_fraction (available, transfer) RESULT(rate)
   ! Bound the transfer without evaluating a potentially overflowing ratio
   ! when the available water already exceeds the requested transfer.
   ! An empty donor must still suppress flux, even if flux*dt underflows.
      real(r8), intent(in) :: available, transfer
      real(r8) :: rate

      rate = 1._r8
      IF (available <= 0._r8) THEN
         rate = 0._r8
      ELSEIF (transfer > available) THEN
         rate = available / transfer
      ENDIF
   END FUNCTION bif_limiter_fraction


   ! =========================================================================
   FUNCTION available_storage_ucat (i, wdsrf, volwater_ucat_in, volwater_ucat_valid_in, &
      volresv_in, is_built_resv_in) RESULT(storage)
   ! =========================================================================
   !
   ! Current host-model available storage semantics for limiter reference.
   !
   ! =========================================================================

   IMPLICIT NONE

   integer,  intent(in) :: i
   real(r8), intent(in) :: wdsrf
   real(r8), intent(in) :: volwater_ucat_in(:)
   logical,  intent(in) :: volwater_ucat_valid_in
      real(r8), intent(in) :: volresv_in(:)
      logical,  intent(in) :: is_built_resv_in(:)
      real(r8)             :: storage
      logical              :: has_levee_cell
      real(r8), parameter  :: stage_restart_tol = 1.e-5_r8

      storage = 0._r8

      IF (i < 1) RETURN

      ! Size-check is_built_resv_in before dereferencing it; the dummy is
      ! caller-provided and we cannot rely on a fixed length contract.
      IF (allocated(ucat2resv) .and. i <= size(is_built_resv_in)) THEN
         IF (is_built_resv_in(i)) THEN
            IF (i <= size(ucat2resv)) THEN
               IF (ucat2resv(i) > 0 .and. ucat2resv(i) <= size(volresv_in)) THEN
                  storage = volresv_in(ucat2resv(i))
                  RETURN
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      has_levee_cell = .false.
      IF (DEF_USE_LEVEE .and. allocated(has_levee)) THEN
         IF (i <= size(has_levee)) has_levee_cell = has_levee(i)
      ENDIF

         IF (volwater_ucat_valid_in) THEN
            IF (i <= size(volwater_ucat_in)) THEN
               IF (volwater_ucat_in(i) > 0._r8 .or. wdsrf <= stage_restart_tol) THEN
                  IF (has_levee_cell) THEN
                     storage = volwater_ucat_in(i) + MAX(levsto(i), 0._r8)
                  ELSE
                     storage = volwater_ucat_in(i) + 0._r8
                  ENDIF
               ELSEIF (has_levee_cell) THEN
                  ! Old restart compatibility: some restarts store
                  ! volwater_ucat as zero placeholders for wet levee cells.
                  storage = levee_visible_volume_from_stage(i, wdsrf, levsto(i)) + MAX(levsto(i), 0._r8)
               ELSE
                  storage = storage_from_stage(i, wdsrf)
               ENDIF
            ELSE
               storage = storage_from_stage(i, wdsrf)
            ENDIF
         ELSEIF (has_levee_cell .and. allocated(levsto)) THEN
            storage = levee_visible_volume_from_stage(i, wdsrf, levsto(i)) + MAX(levsto(i), 0._r8)
         ELSE
            storage = storage_from_stage(i, wdsrf)
         ENDIF

      END FUNCTION available_storage_ucat


      ! =========================================================================
      FUNCTION available_visible_storage_ucat (i, wdsrf, volwater_ucat_in, volwater_ucat_valid_in, &
      volresv_in, is_built_resv_in) RESULT(storage)
      ! =========================================================================
      !
      ! Visible-side donor storage for layer-1 and non-levee layer-2+ bif
      ! fluxes.  Reservoirs remain a single visible compartment.
      !
      ! =========================================================================

      IMPLICIT NONE

      integer,  intent(in) :: i
      real(r8), intent(in) :: wdsrf
   real(r8), intent(in) :: volwater_ucat_in(:)
   logical,  intent(in) :: volwater_ucat_valid_in
      real(r8), intent(in) :: volresv_in(:)
      logical,  intent(in) :: is_built_resv_in(:)
      real(r8)             :: storage
      logical              :: has_levee_cell
      real(r8), parameter  :: stage_restart_tol = 1.e-5_r8

         storage = 0._r8
         IF (i < 1) RETURN

         IF (allocated(ucat2resv) .and. i <= size(is_built_resv_in)) THEN
            IF (is_built_resv_in(i)) THEN
               IF (i <= size(ucat2resv)) THEN
                  IF (ucat2resv(i) > 0 .and. ucat2resv(i) <= size(volresv_in)) THEN
                     storage = max(volresv_in(ucat2resv(i)), 0._r8)
                     RETURN
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         has_levee_cell = .false.
         IF (DEF_USE_LEVEE .and. allocated(has_levee)) THEN
            IF (i <= size(has_levee)) has_levee_cell = has_levee(i)
         ENDIF

         IF (volwater_ucat_valid_in .and. i <= size(volwater_ucat_in)) THEN
            IF (volwater_ucat_in(i) > 0._r8 .or. wdsrf <= stage_restart_tol) THEN
               storage = volwater_ucat_in(i)
            ELSEIF (has_levee_cell) THEN
               storage = levee_visible_volume_from_stage(i, wdsrf, levsto(i))
            ELSE
               storage = storage_from_stage(i, wdsrf)
            ENDIF
         ELSEIF (has_levee_cell .and. allocated(levsto)) THEN
            storage = levee_visible_volume_from_stage(i, wdsrf, levsto(i))
         ELSE
            storage = storage_from_stage(i, wdsrf)
         ENDIF
         storage = max(storage, 0._r8)

      END FUNCTION available_visible_storage_ucat


      ! =========================================================================
      FUNCTION available_protected_storage_ucat (i, is_built_resv_in) RESULT(storage)
      ! =========================================================================
      !
      ! Protected-side donor storage for leveed layer-2+ bif fluxes.
      !
      ! =========================================================================

      IMPLICIT NONE

      integer, intent(in) :: i
      logical, intent(in) :: is_built_resv_in(:)
      real(r8)            :: storage

         storage = 0._r8
         IF (i < 1) RETURN
         IF (i <= size(is_built_resv_in)) THEN
            IF (is_built_resv_in(i)) RETURN
         ENDIF
         IF (DEF_USE_LEVEE .and. allocated(has_levee) .and. allocated(levsto)) THEN
            IF (i <= size(has_levee) .and. i <= size(levsto)) THEN
               IF (has_levee(i)) storage = max(levsto(i), 0._r8)
            ENDIF
         ENDIF

      END FUNCTION available_protected_storage_ucat


      ! =========================================================================
      FUNCTION storage_from_stage (i, wdsrf) RESULT(storage)
   ! =========================================================================
   !
   ! Ordinary cells must use the same stage-volume curve as the host.
   ! Leveed cells use levee_visible_volume_from_stage at their call sites.
   !
   ! =========================================================================

   IMPLICIT NONE

   integer,  intent(in) :: i
   real(r8), intent(in) :: wdsrf
   real(r8)             :: storage

      storage = floodplain_curve(i)%volume(wdsrf)
      storage = MAX(storage, 0._r8)

   END FUNCTION storage_from_stage


   ! =========================================================================
   SUBROUTINE build_bifurcation_path_signature (signature)
   ! =========================================================================
   !
   ! Build lossless restart metadata for each pathway.  The rows are:
   !   schema version, upstream global UCID, downstream global UCID, distance,
   !   elevation(:), width(:), Manning(:).
   ! Values are copied directly from the network input, so exact comparison
   ! after a NetCDF real8 round trip is intentional.
   !
   ! =========================================================================

   USE MOD_Grid_RiverLakeNetwork, only: npthlev_bif, npthout_local, &
      pth_upst_local, pth_down_ucid, pth_dst, pth_elv, pth_wth, pth_man, ucat_ucid
   IMPLICIT NONE

   real(r8), allocatable, intent(out) :: signature(:,:)
   integer :: ipth, iup, npath

      npath = 0
      IF (p_is_worker) npath = npthout_local
      allocate (signature(4 + 3*npthlev_bif, npath))
      signature(:,:) = 0._r8

      IF (.not. p_is_worker) RETURN

      DO ipth = 1, npthout_local
         iup = pth_upst_local(ipth)
         IF (iup < 1 .or. iup > size(ucat_ucid)) THEN
            CALL CoLM_stop ('build_bifurcation_path_signature: invalid upstream local index')
         ENDIF

         signature(1, ipth) = real(BIF_RESTART_SIGNATURE_VERSION, r8)
         signature(2, ipth) = real(ucat_ucid(iup), r8)
         signature(3, ipth) = real(pth_down_ucid(ipth), r8)
         signature(4, ipth) = pth_dst(ipth)
         signature(5:4+npthlev_bif, ipth) = pth_elv(:, ipth)
         signature(5+npthlev_bif:4+2*npthlev_bif, ipth) = pth_wth(:, ipth)
         signature(5+2*npthlev_bif:4+3*npthlev_bif, ipth) = pth_man(:)
      ENDDO

   END SUBROUTINE build_bifurcation_path_signature


   ! =========================================================================
   SUBROUTINE read_bifurcation_restart (file_restart, previous_depth_restart_found, restart_loaded, &
      restart_transaction_validated_in, restart_feature_manifest_present_in, &
      restart_bifurcation_enabled_in, restart_levee_enabled_in)
   ! =========================================================================
   !
   ! Read bifurcation pathway state from restart in global pathway order.
   !
   ! =========================================================================

   USE MOD_NetCDFSerial, only: ncio_var_exist
   USE MOD_Grid_RiverLakeVectorIO, only: vector_read_matrix_and_scatter
   USE MOD_Grid_RiverLakeNetwork, only: npthlev_bif, npthout_local, totalnpthout, pth_global_id
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart
   logical, intent(in) :: previous_depth_restart_found
   logical, intent(out) :: restart_loaded
   logical, intent(in) :: restart_transaction_validated_in
   logical, intent(in) :: restart_feature_manifest_present_in
   logical, intent(in) :: restart_bifurcation_enabled_in
   logical, intent(in) :: restart_levee_enabled_in
   logical :: has_pth_veloc, has_pth_momen, has_path_signature
   logical :: restart_feature_present, strict_bif_restart, state_allocated
   integer, allocatable :: global_id_read(:)
   real(r8), allocatable :: matrix_dummy(:,:), path_signature(:,:), current_signature(:,:)
   integer :: has_flags(3)
   integer :: ipth, ilev, mismatch_count, first_mismatch_gid, invalid_state_count

      restart_loaded = .false.

      IF (.not. restart_transaction_validated_in) THEN
         CALL CoLM_stop ('read_bifurcation_restart requires validated transaction context')
      ENDIF

      ! vector_read_matrix_and_scatter contains mpi_barrier(p_comm_glb)
      ! and mpi_recv/mpi_send on p_comm_glb, so every rank in that
      ! communicator (master, workers AND IO ranks) must enter this
      ! routine or the barrier hangs forever. The IO-only ranks take the
      ! "ELSE" branch below along with the master-when-non-worker case,
      ! passing zero-length dummy arrays and participating only in the
      ! collectives. Do NOT short-circuit with `p_is_master .or.
      ! p_is_worker` here.
      IF (npthlev_bif <= 0 .or. totalnpthout <= 0) RETURN
      state_allocated = allocated(pth_veloc) .and. allocated(pth_momen)

      ! Restart state is opt-in only after its pathway identity has been
      ! validated. This also makes repeated reads of a legacy/incomplete file
      ! a deterministic cold start instead of retaining stale in-memory state.
      IF (state_allocated .and. p_is_worker) THEN
         pth_veloc(:,:) = 0._r8
         pth_momen(:,:) = 0._r8
      ENDIF

      ! Master probes restart variable existence then broadcasts the
      ! result. Having all ranks call ncio_var_exist directly would mean
      ! p_np_glb concurrent nf90_open calls on the same file; on an
      ! HDF5-backed NetCDF4 build the default file locking serialises or
      ! deadlocks those opens. One master-only probe + broadcast avoids
      ! the storm entirely and keeps all ranks on the same boolean.
      has_flags(:) = 0
      IF (p_is_master) THEN
         has_pth_veloc = ncio_var_exist(file_restart, 'pth_veloc', readflag = .false.)
         has_pth_momen = ncio_var_exist(file_restart, 'pth_momen', readflag = .false.)
         has_path_signature = ncio_var_exist(file_restart, 'bif_path_signature', readflag = .false.)
         has_flags(1) = merge(1, 0, has_pth_veloc)
         has_flags(2) = merge(1, 0, has_pth_momen)
         has_flags(3) = merge(1, 0, has_path_signature)
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast (has_flags, 3, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif
      has_pth_veloc = (has_flags(1) /= 0)
      has_pth_momen = (has_flags(2) /= 0)
      has_path_signature = (has_flags(3) /= 0)
      restart_feature_present = restart_feature_manifest_present_in
      strict_bif_restart = restart_feature_present .and. restart_bifurcation_enabled_in

      IF (strict_bif_restart .and. &
          (.not. state_allocated .or. .not. has_pth_veloc .or. .not. has_pth_momen .or. &
           .not. has_path_signature .or. .not. previous_depth_restart_found)) THEN
         CALL CoLM_stop (&
            'GridRiverLake restart declares bifurcation enabled but pathway state is incomplete')
      ENDIF
      IF (.not. state_allocated) RETURN

      IF ((.not. restart_feature_present .or. strict_bif_restart) .and. &
          has_pth_veloc .and. has_pth_momen .and. has_path_signature .and. &
          previous_depth_restart_found) THEN
         ! Validate identity before reading ordinally indexed momentum.
         CALL build_bifurcation_path_signature (current_signature)

         IF (p_is_worker) THEN
            CALL vector_read_matrix_and_scatter ( &
               file_restart, path_signature, 4 + 3*npthlev_bif, npthout_local, &
               'bif_path_signature', pth_global_id, totalnpthout)
         ELSE
            allocate (global_id_read(0))
            CALL vector_read_matrix_and_scatter ( &
               file_restart, matrix_dummy, 4 + 3*npthlev_bif, 0, &
               'bif_path_signature', global_id_read, totalnpthout)
            deallocate (global_id_read)
         ENDIF

         mismatch_count = 0
         first_mismatch_gid = huge(first_mismatch_gid)
         IF (p_is_worker) THEN
            DO ipth = 1, npthout_local
               IF (any(.not. ieee_is_finite(path_signature(:, ipth)))) THEN
                  mismatch_count = mismatch_count + 1
                  first_mismatch_gid = min(first_mismatch_gid, pth_global_id(ipth))
               ELSEIF (any(path_signature(:, ipth) /= current_signature(:, ipth))) THEN
                  mismatch_count = mismatch_count + 1
                  first_mismatch_gid = min(first_mismatch_gid, pth_global_id(ipth))
               ENDIF
            ENDDO
         ENDIF
#ifdef USEMPI
         CALL mpi_allreduce (MPI_IN_PLACE, mismatch_count, 1, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
         CALL mpi_allreduce (MPI_IN_PLACE, first_mismatch_gid, 1, MPI_INTEGER, MPI_MIN, p_comm_glb, p_err)
#endif
         deallocate (current_signature)
         IF (allocated(path_signature)) deallocate (path_signature)

         IF (mismatch_count > 0) THEN
            IF (p_is_master) THEN
               write(*,'(A,I0,A,I0)') 'ERROR: bifurcation restart pathway identity mismatch; count=', &
                  mismatch_count, ', first global pathway=', first_mismatch_gid
               call flush(6)
            ENDIF
            CALL CoLM_stop ('Refusing to load bifurcation momentum for a different pathway network')
         ENDIF

         IF (p_is_worker) THEN
            CALL vector_read_matrix_and_scatter ( &
               file_restart, pth_veloc, npthlev_bif, npthout_local, 'pth_veloc', pth_global_id, totalnpthout)
            CALL vector_read_matrix_and_scatter ( &
               file_restart, pth_momen, npthlev_bif, npthout_local, 'pth_momen', pth_global_id, totalnpthout)
         ELSE
            allocate (global_id_read(0))
            CALL vector_read_matrix_and_scatter ( &
               file_restart, matrix_dummy, npthlev_bif, 0, 'pth_veloc', global_id_read, totalnpthout)
            CALL vector_read_matrix_and_scatter ( &
               file_restart, matrix_dummy, npthlev_bif, 0, 'pth_momen', global_id_read, totalnpthout)
            deallocate (global_id_read)
         ENDIF
         restart_loaded = .true.
      ELSE
         IF (restart_feature_present .and. .not. strict_bif_restart) THEN
            IF (p_is_master) THEN
               write(*,'(A)') 'GridRiverLake restart declares bifurcation disabled; cold-starting pathway state.'
               call flush(6)
            ENDIF
         ELSE
            ! Old files have no pathway identity metadata. Loading their
            ! ordinally indexed state could bind momentum to a different path
            ! after a same-size network reorder, so cold-start BIF safely.
            IF (p_is_master) THEN
               write(*,'(A,L1,A,L1,A,L1,A,L1,A)') &
                  'WARNING: incomplete/legacy bifurcation restart (pth_veloc=', has_pth_veloc, &
                  ', pth_momen=', has_pth_momen, ', bif_path_signature=', has_path_signature, &
                  ', wdsrf_ucat_prev=', previous_depth_restart_found, &
                  '); cold-starting pathway velocity and momentum at zero.'
               call flush(6)
            ENDIF
         ENDIF
      ENDIF

      ! Validate without compound logical expressions: Fortran need not
      ! short-circuit `.or.`, so comparing abs(NaN) under -ffpe-trap=invalid
      ! can abort before the restart guard runs. Velocity is clamped to
      ! +/-20 m/s during routing; 50 m/s and 1e4 m2/s are generous restart
      ! corruption bounds. Any bad pathway cold-starts the whole paired
      ! (previous depth, pathway momentum) state unit via restart_loaded.
      invalid_state_count = 0
      IF (restart_loaded .and. p_is_worker .and. npthout_local > 0) THEN
         DO ipth = 1, npthout_local
            DO ilev = 1, npthlev_bif
               IF (.not. ieee_is_finite(pth_veloc(ilev, ipth))) THEN
                  invalid_state_count = invalid_state_count + 1
               ELSEIF (abs(pth_veloc(ilev, ipth)) > 50._r8) THEN
                  invalid_state_count = invalid_state_count + 1
               ELSEIF (.not. ieee_is_finite(pth_momen(ilev, ipth))) THEN
                  invalid_state_count = invalid_state_count + 1
               ELSEIF (abs(pth_momen(ilev, ipth)) > 1.e4_r8) THEN
                  invalid_state_count = invalid_state_count + 1
               ENDIF
            ENDDO
         ENDDO
      ENDIF
#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, invalid_state_count, 1, MPI_INTEGER, &
         MPI_SUM, p_comm_glb, p_err)
#endif
      IF (invalid_state_count > 0) THEN
         IF (p_is_worker) THEN
            pth_veloc = 0._r8
            pth_momen = 0._r8
         ENDIF
         restart_loaded = .false.
         IF (strict_bif_restart) THEN
            CALL CoLM_stop (&
               'GridRiverLake restart declares bifurcation enabled but pathway state is invalid')
         ENDIF
         IF (p_is_master) THEN
            write(*,'(A,I0,A)') 'WARNING: invalid bifurcation restart state (count=', &
               invalid_state_count, '); cold-starting paired pathway state.'
            call flush(6)
         ENDIF
      ENDIF

      ! A levee-mode change changes pathway water surfaces and depth rules.
      ! Keep all identity/corruption checks above, then cold-start the paired
      ! momentum/previous-depth state without discarding stored water or history.
      IF (restart_loaded .and. restart_feature_present .and. &
          (restart_levee_enabled_in .neqv. DEF_USE_LEVEE)) THEN
         IF (p_is_worker) THEN
            pth_veloc = 0._r8
            pth_momen = 0._r8
         ENDIF
         restart_loaded = .false.
         IF (p_is_master) THEN
            write(*,'(A)') 'WARNING: levee mode changed; cold-starting paired bifurcation state.'
            call flush(6)
         ENDIF
      ENDIF

      CALL debug_check_bifurcation_restart_state (restart_loaded, restart_loaded)

   END SUBROUTINE read_bifurcation_restart

   SUBROUTINE debug_check_bifurcation_restart_state (has_pth_veloc, has_pth_momen)

   USE MOD_Grid_RiverLakeNetwork, only: npthlev_bif, npthout_local, totalnpthout, pth_global_id
   IMPLICIT NONE

   logical, intent(in) :: has_pth_veloc
   logical, intent(in) :: has_pth_momen

   integer :: nbad_gid
   integer :: nnan_veloc
   integer :: nnan_momen
   real(r8) :: maxabs_veloc
   real(r8) :: maxabs_momen

      IF (dbg_bif_restart_checked) RETURN
      dbg_bif_restart_checked = .true.

      IF (.not. p_is_worker) RETURN

      nbad_gid = 0
      nnan_veloc = 0
      nnan_momen = 0
      maxabs_veloc = 0._r8
      maxabs_momen = 0._r8

      IF (npthout_local > 0) THEN
         nbad_gid = count(pth_global_id < 1 .or. pth_global_id > totalnpthout)
         nnan_veloc = count(pth_veloc /= pth_veloc)
         nnan_momen = count(pth_momen /= pth_momen)
         maxabs_veloc = maxval(abs(pth_veloc))
         maxabs_momen = maxval(abs(pth_momen))
      ENDIF

      IF (nbad_gid > 0 .or. nnan_veloc > 0 .or. nnan_momen > 0) THEN
         write(*,'(A,I0,A,L1,A,L1,A,I0,A,I0,A,I0,A,ES12.4,A,ES12.4)') &
            'WARNING bifurcation restart state glb=', p_iam_glb, &
            ' has_veloc=', has_pth_veloc, &
            ' has_momen=', has_pth_momen, &
            ' bad_gid=', nbad_gid, &
            ' nan_vel=', nnan_veloc, &
            ' nan_mom=', nnan_momen, &
            ' maxabs_vel=', maxabs_veloc, &
            ' maxabs_mom=', maxabs_momen
         call flush(6)
      ENDIF

   END SUBROUTINE debug_check_bifurcation_restart_state


   ! =========================================================================
   SUBROUTINE write_bifurcation_restart (file_restart)
   ! =========================================================================
   !
   ! Write bifurcation pathway state in global pathway order.
   !
   ! =========================================================================

   USE MOD_Namelist, only: DEF_REST_CompressLevel
   USE MOD_NetCDFSerial, only: ncio_define_dimension, ncio_write_serial
   USE MOD_Grid_RiverLakeVectorIO, only: vector_gather_matrix_to_master
   USE MOD_Grid_RiverLakeNetwork, only: npthlev_bif, npthout_local, totalnpthout, pth_global_id
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

   integer :: ncol_local_write
   integer, allocatable :: pth_global_id_write(:)
   real(r8), allocatable :: wdata(:,:), path_signature_write(:,:)
   real(r8), allocatable :: pth_veloc_write(:,:), pth_momen_write(:,:)

      IF (totalnpthout <= 0 .or. npthlev_bif <= 0) RETURN

      ! Guard: bifurcation state may not be initialised (e.g. mkinidata).
      ! Write zero state so the restart file has valid variables.
      IF (p_is_worker) THEN
         IF (allocated(pth_veloc) .and. allocated(pth_momen)) THEN
            ncol_local_write = npthout_local
            pth_veloc_write = pth_veloc
            pth_momen_write = pth_momen
            allocate (pth_global_id_write (ncol_local_write))
            pth_global_id_write(:) = pth_global_id(:)
         ELSE
            ! Not initialised: contribute zero-state columns
            ncol_local_write = npthout_local
            allocate (pth_veloc_write     (npthlev_bif, ncol_local_write))
            allocate (pth_momen_write     (npthlev_bif, ncol_local_write))
            allocate (pth_global_id_write (ncol_local_write))
            pth_veloc_write = 0._r8
            pth_momen_write = 0._r8
            pth_global_id_write(:) = pth_global_id(:)
         ENDIF
      ELSE
         ncol_local_write = 0
         allocate (pth_veloc_write     (npthlev_bif, 0))
         allocate (pth_momen_write     (npthlev_bif, 0))
         allocate (pth_global_id_write (0))
      ENDIF

      IF (p_is_master) THEN
         CALL ncio_define_dimension(file_restart, 'bifurcation_signature_field', 4 + 3*npthlev_bif)
         CALL ncio_define_dimension(file_restart, 'bifurcation_level',   npthlev_bif)
         CALL ncio_define_dimension(file_restart, 'bifurcation_pathway', totalnpthout)
      ENDIF

      ! Persist the identity before the ordinally indexed state. A reader must
      ! validate this lossless metadata before loading velocity or momentum.
      CALL build_bifurcation_path_signature (path_signature_write)
      CALL vector_gather_matrix_to_master ( &
         path_signature_write, 4 + 3*npthlev_bif, ncol_local_write, &
         totalnpthout, pth_global_id_write, wdata)

      IF (p_is_master) THEN
         CALL ncio_write_serial (file_restart, 'bif_path_signature', wdata, &
            'bifurcation_signature_field', 'bifurcation_pathway', DEF_REST_CompressLevel)
         deallocate (wdata)
      ENDIF
      deallocate (path_signature_write)

      ! Gather local (npthlev_bif, npthout_local) state into global
      ! (npthlev_bif, totalnpthout) pathway order using pth_global_id.
      CALL vector_gather_matrix_to_master ( &
         pth_veloc_write, npthlev_bif, ncol_local_write, totalnpthout, pth_global_id_write, wdata)

      IF (p_is_master) THEN
         CALL ncio_write_serial (file_restart, 'pth_veloc', wdata, &
            'bifurcation_level', 'bifurcation_pathway', DEF_REST_CompressLevel)
         deallocate (wdata)
      ENDIF

      CALL vector_gather_matrix_to_master ( &
         pth_momen_write, npthlev_bif, ncol_local_write, totalnpthout, pth_global_id_write, wdata)

      IF (p_is_master) THEN
         CALL ncio_write_serial (file_restart, 'pth_momen', wdata, &
            'bifurcation_level', 'bifurcation_pathway', DEF_REST_CompressLevel)
         deallocate (wdata)
      ENDIF

      deallocate (pth_veloc_write)
      deallocate (pth_momen_write)
      deallocate (pth_global_id_write)

   END SUBROUTINE write_bifurcation_restart


   ! =========================================================================
   SUBROUTINE bifurcation_final ()
   ! =========================================================================
   !
   ! Deallocate all bifurcation arrays.
   !
   ! =========================================================================

   IMPLICIT NONE

      IF (allocated(pth_veloc))     deallocate (pth_veloc)
      IF (allocated(bif_ucat))      deallocate (bif_ucat)
      IF (allocated(pth_momen))     deallocate (pth_momen)
      IF (allocated(bif_hflux_lev)) deallocate (bif_hflux_lev)
      IF (allocated(bif_hflux_sum)) deallocate (bif_hflux_sum)
      IF (allocated(bif_path_active)) deallocate (bif_path_active)
      IF (allocated(bif_lev_hflux_sum)) deallocate (bif_lev_hflux_sum)

      IF (allocated(wdsrf_dn_pth_buf    )) deallocate (wdsrf_dn_pth_buf    )
      IF (allocated(wdsrf_prev_dn_pth_buf)) deallocate (wdsrf_prev_dn_pth_buf)
      IF (allocated(rivelv_dn_pth_buf   )) deallocate (rivelv_dn_pth_buf   )
      IF (allocated(protected_wdsrf_ucat_buf  )) deallocate (protected_wdsrf_ucat_buf  )
      IF (allocated(protected_wdsrf_dn_pth_buf)) deallocate (protected_wdsrf_dn_pth_buf)
      IF (allocated(has_levee_dn_pth_buf)) deallocate (has_levee_dn_pth_buf)
      IF (allocated(storage_dn_pth_buf  )) deallocate (storage_dn_pth_buf  )
      IF (allocated(pth_hflux_total_buf )) deallocate (pth_hflux_total_buf )
      IF (allocated(pth_limiter_rate_buf)) deallocate (pth_limiter_rate_buf)
         IF (allocated(bif_influx_buf      )) deallocate (bif_influx_buf      )
         IF (allocated(has_levee_r8_buf    )) deallocate (has_levee_r8_buf    )
         IF (allocated(storage_ucat_buf    )) deallocate (storage_ucat_buf    )
         IF (allocated(visible_storage_ucat_buf  )) deallocate (visible_storage_ucat_buf  )
         IF (allocated(protected_storage_ucat_buf)) deallocate (protected_storage_ucat_buf)
         IF (allocated(visible_storage_dn_pth_buf  )) deallocate (visible_storage_dn_pth_buf  )
         IF (allocated(protected_storage_dn_pth_buf)) deallocate (protected_storage_dn_pth_buf)
         IF (allocated(limiter_outgoing_buf)) deallocate (limiter_outgoing_buf)
         IF (allocated(limiter_out_rate_buf)) deallocate (limiter_out_rate_buf)
         IF (allocated(protected_outgoing_buf)) deallocate (protected_outgoing_buf)
         IF (allocated(protected_out_rate_buf)) deallocate (protected_out_rate_buf)
         IF (allocated(layer_limiter_rate_buf)) deallocate (layer_limiter_rate_buf)

      IF (allocated(neg_pth_buf             )) deallocate (neg_pth_buf             )
      IF (allocated(bif_neg_recv_buf        )) deallocate (bif_neg_recv_buf        )
      IF (allocated(limiter_out_rate_pth_buf)) deallocate (limiter_out_rate_pth_buf)
      IF (allocated(protected_neg_pth_buf   )) deallocate (protected_neg_pth_buf   )
      IF (allocated(protected_out_recv_buf  )) deallocate (protected_out_recv_buf  )
      IF (allocated(protected_out_rate_pth_buf)) deallocate (protected_out_rate_pth_buf)

      IF (allocated(ucatfilter_r8_buf       )) deallocate (ucatfilter_r8_buf       )
      IF (allocated(ucatfilter_dn_pth_buf   )) deallocate (ucatfilter_dn_pth_buf   )
      IF (allocated(is_resv_r8_buf          )) deallocate (is_resv_r8_buf          )
      IF (allocated(is_resv_dn_pth_buf      )) deallocate (is_resv_dn_pth_buf      )

      IF (allocated(pth_lev_hflux_total_buf )) deallocate (pth_lev_hflux_total_buf )
      IF (allocated(bif_lev_influx_buf      )) deallocate (bif_lev_influx_buf      )

   END SUBROUTINE bifurcation_final

END MODULE MOD_Grid_RiverLakeBifurcation
#endif
