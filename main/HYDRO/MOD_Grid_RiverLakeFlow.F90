#include <define.h>

#ifdef GridRiverLakeFlow
MODULE MOD_Grid_RiverLakeFlow
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   River Lake flow.
!
!-------------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Grid_RiverLakeNetwork
   USE MOD_Grid_RiverLakeTimeVars
   USE MOD_Grid_Reservoir
   ! Outflow of the cells directly upstream of a built reservoir (CaMa-Flood's
   ! kinematic treatment against the reservoir backwater; see that module).
   USE MOD_Grid_ReservoirUpstream, only: reservoir_upstream_init, reservoir_upstream_mark, &
      reservoir_upstream_outflow, reservoir_upstream_final
   USE MOD_Grid_RiverLakeHist
   USE MOD_Grid_RiverLakeLevee, only: has_levee, levsto, levdph, levee_init, read_levee_restart, &
      levee_visible_volume_from_stage, levee_repartition_storage, levee_apply_protected_flux, levee_final, &
      levee_fldstg
   USE MOD_Grid_RiverLakeBifurcation, only: bifurcation_init, read_bifurcation_restart, &
      bifurcation_final, bifurcation_calc, &
      bif_hflux_sum, bif_hflux_lev, bif_lev_hflux_sum, bif_path_active
#ifdef MaterialTransport
   ! Material transport: the mt_rlc_* calls below (one per marked block) change
   ! only material pools, never the river water state, so the host solution
   ! is the same as without them.
   USE MOD_MT_RiverLakeCoupling, only: mt_rlc_init, mt_rlc_runoff_input, mt_rlc_begin_routing, &
      mt_rlc_levee_repartition, mt_rlc_equilibrate_cell, mt_rlc_substep, mt_rlc_inland_depression, &
      mt_rlc_diag_accumulate, mt_rlc_debug_discharge, mt_rlc_end_routing, mt_rlc_flood_init, &
      mt_rlc_flood_withhold, mt_rlc_flood_publish_reset, mt_rlc_flood_publish_cell, &
      mt_rlc_flood_publish_push, mt_rlc_flood_debit, mt_rlc_final
#endif
   IMPLICIT NONE

   real(r8), parameter :: RIVERMIN  = 1.e-5_r8
   real(r8), parameter :: RIVERLAKE_FLOOD_MISSING_VALUE = -1.e30_r8

   ! --- Adaptive routing substep protection ---------------------------------
   ! A near-dry unit catchment with a positive net outflow would otherwise set
   ! dt = volwater/sum_hflux_riv -> 0.  With DEF_USE_BIFURCATION the smallest
   ! dt of any cell on any worker is synchronized to every worker, so the
   ! substep loop would need acctime_rnof/dt iterations and print nothing
   ! until it finishes -- the job then looks hung for hours or days.  Under
   ! bifurcation keep the storage constraint out of the near-empty regime (the
   ! gross outgoing-face limiter caps those cells; without bifurcation there is
   ! no such limiter, so the constraint stays in force), and in every mode fail
   ! loudly if the substep count explodes.  The substep limits are given as a
   ! mean substep dt, so they scale with the routing interval (a healthy
   ! 0.1-degree US run averages ~8.6 s, i.e. ~420 substeps per 3600 s).
   real(r8), parameter :: ROUTING_STORAGE_DT_EPS   = 1.e-6_r8  ! [m3]
   real(r8), parameter :: ROUTING_SUBSTEP_WARN_DT  = 2._r8     ! [s] mean dt that warns
   real(r8), parameter :: ROUTING_SUBSTEP_ABORT_DT = 0.2_r8    ! [s] mean dt that aborts
   ! Global adaptive dt of the last completed substep (diagnostics only).
   real(r8), save :: routing_last_global_dt = 0._r8

   real(r8), save :: acctime_rnof_max
   integer, save :: routing_zero_dt_warn_count = 0

   logical,  allocatable :: filter_rnof (:)
   real(r8), allocatable :: flood_depth_patch(:), flood_fraction_patch(:)
   real(r8), allocatable :: flood_evap_patch(:), flood_infil_patch(:)
   real(r8), allocatable :: flood_visible_uc(:), flood_protected_uc(:)
   logical,  allocatable :: flood_reservoir_uc(:)
   real(r8), allocatable :: flood_credit_patch(:)
   real(r8), allocatable :: flood_grid_area(:)
   real(r8), allocatable :: flood_evap_acc(:), flood_infil_acc(:)
   real(r8), save :: flood_evap_period = 0._r8, flood_infil_period = 0._r8
#ifdef CoLMDEBUG
   ! land-side flood sink [m3] over the same period, to check it against the river-side debit
   real(r8), save :: flood_land_period = 0._r8
#endif

CONTAINS

   ! ---------
   SUBROUTINE grid_riverlake_flow_init (start_year, spinup)

   USE MOD_LandPatch,           only: numpatch
   USE MOD_Forcing,             only: forcmask_pch
   USE MOD_Vars_TimeInvariants, only: patchtype, patchmask
   USE MOD_Vars_Global,         only: spval
   IMPLICIT NONE

      integer, intent(in) :: start_year
      logical, intent(in) :: spinup

      integer :: i, irsv
      logical :: bif_restart_loaded
      real(r8) :: migrated_stage
      logical, allocatable :: levee_storage_folded(:)

      acctime_rnof_max = DEF_GRIDBASED_ROUTING_MAX_DT

      IF (p_is_master .and. (.not. DEF_GRIDBASED_ROUTING_STORAGE_DT_LIMIT) &
          .and. (.not. DEF_USE_BIFURCATION)) THEN
         write(*,'(A)') 'Warning: DEF_GRIDBASED_ROUTING_STORAGE_DT_LIMIT = .false. needs ' // &
            'DEF_USE_BIFURCATION (its outgoing-flux limiter); the storage dt bound stays on.'
      ENDIF

      ! excluding (patchtype >= 99), virtual patches and those forcing missed
      IF (p_is_worker) THEN
         allocate (filter_rnof (numpatch))
         IF (numpatch > 0) THEN
            filter_rnof = patchtype < 99
            filter_rnof = filter_rnof .and. patchmask
            IF (DEF_forcing%has_missing_value) THEN
               filter_rnof = filter_rnof .and. forcmask_pch
            ENDIF
         ENDIF
      ENDIF


      ! static slope of the cells upstream of reservoirs (reads topo_distance)
      IF (DEF_Reservoir_Method > 0) CALL reservoir_upstream_init ()

      allocate (levee_storage_folded(numucat))
      levee_storage_folded = .false.
      CALL levee_init()
      IF (len_trim(gridriver_restart_file) > 0) THEN
         CALL read_levee_restart(gridriver_restart_file, &
            restart_transaction_validated, restart_feature_manifest_present, &
            restart_levee_enabled, fold_protected_to_visible = .not. DEF_USE_LEVEE, &
            volwater_ucat_io = volwater_ucat, volwater_ucat_valid_io = volwater_ucat_valid, &
            wdsrf_ucat_in = wdsrf_ucat, storage_folded_out = levee_storage_folded)
      ENDIF

      IF (DEF_USE_BIFURCATION) THEN
         CALL bifurcation_init()
         IF (len_trim(gridriver_restart_file) > 0) THEN
            CALL read_bifurcation_restart(gridriver_restart_file, &
               wdsrf_ucat_prev_restart_found, bif_restart_loaded, &
               restart_transaction_validated, restart_feature_manifest_present, &
               restart_bifurcation_enabled, restart_levee_enabled)
            IF (.not. bif_restart_loaded) THEN
               wdsrf_ucat_prev = wdsrf_ucat
               wdsrf_ucat_prev_valid = .true.
            ENDIF
         ENDIF
      ENDIF

      ! Folding changes storage; disabling levees also changes the stage-volume
      ! curve when the old protected pool was empty. Preserve extrapolation
      ! history for ordinary cells whose stored stage already matches the curve.
      IF (volwater_ucat_valid .and. p_is_worker) THEN
         DO i = 1, numucat
            IF (lake_type(i) == 2) CYCLE
            IF (.not. levee_storage_folded(i)) THEN
               IF (.not. restart_levee_enabled) CYCLE
               IF (DEF_USE_LEVEE .and. has_levee(i)) CYCLE
               ! Cold-start zero storage is rebuilt from its positive stage below.
               IF (volwater_ucat(i) <= 0._r8 .and. wdsrf_ucat(i) > RIVERMIN) CYCLE
               migrated_stage = floodplain_curve(i)%depth(volwater_ucat(i))
               IF (abs(migrated_stage - wdsrf_ucat(i)) <= &
                   64._r8 * epsilon(1._r8) * max(1._r8, abs(migrated_stage), abs(wdsrf_ucat(i)))) CYCLE
            ELSE
               migrated_stage = floodplain_curve(i)%depth(volwater_ucat(i))
            ENDIF
            wdsrf_ucat(i) = migrated_stage
            wdsrf_ucat_prev(i) = wdsrf_ucat(i)
         ENDDO
      ENDIF
      deallocate (levee_storage_folded)

      CALL rebuild_volwater_ucat ()

#ifdef MaterialTransport
      ! material pools from the river restart, else from the water state
      CALL mt_rlc_init (start_year, gridriver_restart_file)
#endif

      gridriver_restart_file = ''

      CALL allocate_flood_patch ()
      IF (DEF_GridRiverLake_FloodFeedback) THEN
         CALL init_flood_feedback ()
         flood_evap_period = 0._r8
         flood_infil_period = 0._r8
         CALL publish_flood_feedback(start_year)
         IF (spinup) CALL withhold_flood_feedback ()
      ENDIF

   END SUBROUTINE grid_riverlake_flow_init

   ! ---------
   SUBROUTINE allocate_flood_patch ()

   USE MOD_LandPatch, only: numpatch
   IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (allocated(flood_depth_patch)) deallocate(flood_depth_patch, flood_fraction_patch, &
            flood_evap_patch, flood_infil_patch)
         allocate(flood_depth_patch(numpatch), flood_fraction_patch(numpatch), &
            flood_evap_patch(numpatch), flood_infil_patch(numpatch))
         flood_depth_patch = 0._r8
         flood_fraction_patch = 0._r8
         flood_evap_patch = 0._r8
         flood_infil_patch = 0._r8
      ENDIF

   END SUBROUTINE allocate_flood_patch

   ! ---------
   SUBROUTINE rebuild_volwater_ucat ()

   IMPLICIT NONE
   integer :: i

      IF (p_is_worker) THEN
         DO i = 1, numucat
            IF (.not. volwater_ucat_valid .or. &
                (volwater_ucat(i) <= 0._r8 .and. wdsrf_ucat(i) > RIVERMIN)) THEN
               IF (DEF_USE_LEVEE .and. has_levee(i)) THEN
                  volwater_ucat(i) = levee_visible_volume_from_stage(i, wdsrf_ucat(i), levsto(i))
               ELSE
                  volwater_ucat(i) = floodplain_curve(i)%volume(wdsrf_ucat(i))
               ENDIF
            ENDIF
         ENDDO
      ENDIF
      volwater_ucat_valid = .true.

   END SUBROUTINE rebuild_volwater_ucat

   ! ---------
   SUBROUTINE init_flood_feedback ()

   ! Patch and grid arrays of grid flood feedback, built on the current land
   ! patches and patch-to-grid mapping: at start-up and again after LULCC, which
   ! rebuilds both (see MOD_Initialize).

   USE MOD_LandPatch, only: numpatch
   IMPLICIT NONE

   real(r8), allocatable :: grid_area_local(:)

      IF (.not. p_is_worker) RETURN

      IF (allocated(flood_visible_uc)) deallocate(flood_visible_uc, flood_protected_uc, &
         flood_credit_patch, flood_evap_acc, flood_infil_acc, flood_reservoir_uc, flood_grid_area)

      allocate(flood_visible_uc(numucat), flood_protected_uc(numucat), &
         flood_credit_patch(numpatch), flood_evap_acc(numpatch), flood_infil_acc(numpatch), &
         flood_reservoir_uc(numucat), flood_grid_area(numinpm))
      flood_evap_acc = 0._r8
      flood_infil_acc = 0._r8

      ! land area of each runoff-input grid cell
      allocate(grid_area_local(numinpm))
      grid_area_local = 0._r8
      flood_credit_patch = 1._r8
      IF (numpatch > 0) CALL worker_remap_data_pset2grid(remap_patch2inpm, &
         flood_credit_patch, grid_area_local, fillvalue=0._r8, filter=flood_credit_patch>0._r8)
      CALL worker_push_data(allreduce_inpm, grid_area_local, flood_grid_area, fillvalue=0._r8)
      ! sum_area is allocated only on workers with requests (numinpm > 0)
      IF (numinpm > 0) flood_grid_area = max(flood_grid_area, push_ucat2inpm%sum_area)
      deallocate(grid_area_local)

#ifdef MaterialTransport
      ! material of the flood credit, on the same patches
      CALL mt_rlc_flood_init ()
#endif

   END SUBROUTINE init_flood_feedback

   ! ---------
   SUBROUTINE grid_riverlake_flow_lulcc (year, spinup)

   ! After LULCC.  MOD_Initialize has rebuilt the network on the new land patches
   ! (patch-to-grid mapping, pushes) and the river state has been carried over
   ! (hold/restore_GridRiverLakeTimeVars_lulcc); rebuild what this module keeps on
   ! the land patches.

   USE MOD_LandPatch, only: numpatch
   IMPLICIT NONE

   integer, intent(in) :: year   ! the new year
   logical, intent(in) :: spinup ! the routing is not called during spin-up

      IF (p_is_worker) THEN
         ! the runoff filter is reset every routing call; size it for the new patches
         IF (allocated(filter_rnof)) deallocate(filter_rnof)
         allocate (filter_rnof (numpatch))
         filter_rnof = .false.
      ENDIF

      ! the network was rebuilt: rebuild the slope of the cells upstream of reservoirs
      IF (DEF_Reservoir_Method > 0) THEN
         CALL reservoir_upstream_final ()
         CALL reservoir_upstream_init ()
      ENDIF

      CALL allocate_flood_patch ()
      CALL rebuild_volwater_ucat ()
      IF (DEF_USE_BIFURCATION) THEN
         wdsrf_ucat_prev = wdsrf_ucat
         wdsrf_ucat_prev_valid = .true.
      ENDIF

      IF (DEF_GridRiverLake_FloodFeedback) THEN
         ! The land sinks of the old patches were debited at the last routing call
         ! of the year (flood_evap_acc/flood_infil_acc are settled every call), so
         ! nothing is pending: rebuild the arrays on the new patches and publish the
         ! flood state for them.
         CALL init_flood_feedback ()
         CALL publish_flood_feedback (year)
         ! As at initialisation: during spin-up grid_riverlake_flow is not called,
         ! so what the land took from a published flood would never be debited
         ! from the river storage.
         IF (spinup) CALL withhold_flood_feedback ()
      ENDIF

   END SUBROUTINE grid_riverlake_flow_lulcc

   ! ---------
   SUBROUTINE withhold_flood_feedback ()

   ! Publish no flood to the land (used while the routing is not running).

      IF (p_is_worker) THEN
         flood_credit_patch = 0._r8
         flood_depth_patch = 0._r8
         flood_fraction_patch = 0._r8
      ENDIF
#ifdef MaterialTransport
      CALL mt_rlc_flood_withhold ()
#endif

   END SUBROUTINE withhold_flood_feedback

   ! ---------
   SUBROUTINE grid_riverlake_flow (year, deltime)

   USE MOD_Utils
   USE MOD_Namelist,       only: DEF_Reservoir_Method, &
      DEF_GRIDBASED_ROUTING_MOMENTUM_DT_LIMIT, DEF_GRIDBASED_ROUTING_STORAGE_DT_LIMIT
   USE MOD_Vars_1DFluxes,  only: rnof
   USE MOD_Forcing,        only: forcmask_pch
   USE MOD_Vars_TimeInvariants, only: patchtype, patchmask
   USE MOD_Mesh,           only: numelm
   USE MOD_LandPatch,      only: elm_patch, numpatch
   USE MOD_Const_Physical, only: grav
   USE MOD_Vars_Global,    only: spval
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   IMPLICIT NONE

   integer,  intent(in) :: year
   real(r8), intent(in) :: deltime

   ! Local Variables
   integer  :: i, j, irsv, ntimestep, ipth, i_up
   integer  :: i_local_min, nsub_warn, nsub_abort
   real(r8) :: dt_this
   real(r8) :: dt_local_min, sub_span, dt_report

   real(r8), allocatable :: rnof_gd(:)
   real(r8), allocatable :: rnof_uc(:)


   logical,  allocatable :: is_built_resv(:)
   ! downstream cell is a built reservoir: kinematic outflow (MOD_Grid_ReservoirUpstream)
   logical,  allocatable :: is_resv_upstream(:)
   real(r8), allocatable :: veloc_kine(:)

   real(r8), allocatable :: wdsrf_next(:)
   real(r8), allocatable :: veloc_next(:)

   real(r8), allocatable :: hflux_fc(:)
   real(r8), allocatable :: mflux_fc(:)
   real(r8), allocatable :: zgrad_dn(:)

   real(r8), allocatable :: hflux_resv(:)
   real(r8), allocatable :: mflux_resv(:)

   real(r8), allocatable :: hflux_sumups(:)
   real(r8), allocatable :: mflux_sumups(:)
   real(r8), allocatable :: zgrad_sumups(:)

   real(r8), allocatable :: sum_hflux_riv(:)
   real(r8), allocatable :: sum_hflux_base(:)
   real(r8), allocatable :: sum_mflux_riv(:)
   real(r8), allocatable :: sum_zgrad_riv(:)
   real(r8), allocatable :: normal_outgoing_rate(:), ordinary_scale(:), ordinary_scale_next(:)
   real(r8), allocatable :: volresv_safe(:)
   real(r8), allocatable :: levee_floodarea(:)
   ! volwater_curve(i) = floodplain_curve(i)%volume(wdsrf_ucat(i)), stored right after
   ! every assignment of wdsrf_ucat(i) in this routine.
   ! Reason: without bifurcation and levee the water volume is not carried but
   ! recomputed from the depth, and that was done 2-3 times per unit catchment and
   ! substep (CFL/negative-water constraint, storage update, accumulation) for the
   ! same depth.  The table lookup was the largest single cost of routing on the
   ! global network (~30% of routing CPU), so it is now done once per depth change.
   ! Results are unchanged, bit for bit: wdsrf_ucat(i) of a river catchment is only
   ! assigned at the two places that refresh volwater_curve(i) (before the substep
   ! loop and in the storage update), so every read gets the value the same
   ! function would return for the same depth.
   ! With bifurcation or levee, volwater_ucat is used as before and this is unused.
   real(r8), allocatable :: volwater_curve(:)

   real(r8) :: veloct_fc, height_fc, momen_fc, zsurf_fc
   real(r8) :: bedelv_fc, height_up, height_dn
   real(r8) :: vwave_up, vwave_dn, hflux_up, hflux_dn, mflux_up, mflux_dn
   real(r8) :: volwater, friction, floodarea, rivsto_hist, depression_capacity
   real(r8) :: fldfrc_levee, visible_hflux, protected_hflux, protected_clip
   real(r8) :: vis_vol_bef_lv, levsto_bef_lv, vis_vol_bef_lv2, levsto_bef_lv2
   real(r8) :: totalflood_evap, totalflood_infil, flood_vec(3)
   real(r8),  allocatable :: dt_res(:), dt_all(:)
   logical,   allocatable :: ucatfilter(:)
   logical :: loop_active, next_loop_active
#ifdef CoLMDEBUG
   real(r8) :: totalvol_bef, totalvol_aft, totalrnof, totaldis
#endif


      IF (p_is_worker) THEN

         allocate (rnof_gd (numinpm))
         allocate (rnof_uc (numucat))
         rnof_gd = 0.

         IF (numpatch > 0) THEN
            filter_rnof = patchtype < 99 .and. patchmask
            IF (DEF_forcing%has_missing_value) filter_rnof = filter_rnof .and. forcmask_pch
            CALL worker_remap_data_pset2grid (remap_patch2inpm, rnof, rnof_gd, &
               fillvalue = 0., filter = filter_rnof)
         ENDIF

         IF (numinpm > 0) THEN
            WHERE (push_ucat2inpm%sum_area > 0)
               rnof_gd = rnof_gd / push_ucat2inpm%sum_area
            END WHERE
         ENDIF

         CALL worker_push_data (push_inpm2ucat, rnof_gd, rnof_uc, &
            fillvalue = 0., mode = 'sum')


         IF (numucat > 0) THEN
            acc_rnof_uc = acc_rnof_uc + rnof_uc*1.e-3*deltime
         ENDIF
#ifdef MaterialTransport
         ! material of the runoff, mapped as the runoff water
         CALL mt_rlc_runoff_input (filter_rnof, rnof_uc, deltime)
#endif

         deallocate(rnof_gd)
         deallocate(rnof_uc)

      ENDIF

      IF (DEF_GridRiverLake_FloodFeedback) THEN
         IF (p_is_worker) THEN
            flood_evap_acc = flood_evap_acc + flood_evap_patch*deltime
            flood_infil_acc = flood_infil_acc + flood_infil_patch*deltime
            flood_evap_patch = 0._r8
            flood_infil_patch = 0._r8
         ENDIF
         CALL debit_flood_feedback(totalflood_evap, totalflood_infil)
         flood_evap_period = flood_evap_period + totalflood_evap
         flood_infil_period = flood_infil_period + totalflood_infil
         CALL publish_flood_feedback(year)
      ENDIF

      acctime_rnof = acctime_rnof + deltime

      IF (acctime_rnof+0.01 < acctime_rnof_max) THEN
         RETURN
      ENDIF


      IF (p_is_worker) THEN

         allocate (is_built_resv (numucat))
         allocate (is_resv_upstream (numucat))
         allocate (veloc_kine    (numucat))
         is_resv_upstream = .false.
         veloc_kine = 0._r8
         allocate (wdsrf_next    (numucat))
         allocate (veloc_next    (numucat))
         allocate (hflux_fc      (numucat))
         allocate (mflux_fc      (numucat))
         allocate (zgrad_dn      (numucat))
         allocate (sum_hflux_riv (numucat))
         IF (DEF_USE_BIFURCATION) THEN
            allocate (sum_hflux_base(numucat), normal_outgoing_rate(numucat))
            allocate (ordinary_scale(numucat), ordinary_scale_next(numucat))
         ENDIF
         allocate (sum_mflux_riv (numucat))
         allocate (sum_zgrad_riv (numucat))
         allocate (ucatfilter    (numucat))
         allocate (levee_floodarea(numucat))
         allocate (volwater_curve (numucat))
         levee_floodarea = 0._r8

         allocate (hflux_sumups  (numucat))
         allocate (mflux_sumups  (numucat))
         allocate (zgrad_sumups  (numucat))

         IF (DEF_Reservoir_Method > 0) THEN
            allocate (hflux_resv (numucat))
            allocate (mflux_resv (numucat))
         ENDIF

         allocate (dt_res (numrivsys))
         allocate (dt_all (numrivsys))

         IF (allocated(volresv)) THEN
            allocate (volresv_safe(size(volresv)))
            volresv_safe = volresv
         ELSE
            allocate (volresv_safe(0))
         ENDIF

#ifdef MaterialTransport
         ! the accumulated runoff material joins the input, as acc_rnof_uc below
         CALL mt_rlc_begin_routing ()
#endif

#ifdef CoLMDEBUG
         totalrnof = sum(acc_rnof_uc)
         totalvol_bef = 0.
#endif

         DO i = 1, numucat

            is_built_resv(i) = .false.
            IF (lake_type(i) == 2) THEN
               irsv = ucat2resv(i)
               IF (year >= dam_build_year(irsv)) THEN
                  is_built_resv(i) = .true.
                  IF (volresv(irsv) == spval) THEN
                     volresv(irsv) = floodplain_curve(i)%volume (wdsrf_ucat(i))
                  ELSE
                     wdsrf_ucat(i) = floodplain_curve(i)%depth (volresv(irsv))
                  ENDIF
               ENDIF
            ENDIF

            IF (.not. is_built_resv(i)) THEN
               momen_riv(i) = wdsrf_ucat(i) * veloc_riv(i)
               IF (DEF_USE_BIFURCATION .or. DEF_USE_LEVEE) THEN
                  volwater = volwater_ucat(i)
               ELSE
                  volwater = floodplain_curve(i)%volume (wdsrf_ucat(i))
               ENDIF
            ELSE
               ! water in reservoirs is assumued to be stationary.
               momen_riv(i) = 0
               veloc_riv(i) = 0
               volwater = volresv(ucat2resv(i))
            ENDIF

#ifdef CoLMDEBUG
            totalvol_bef = totalvol_bef + volwater
            IF (DEF_USE_LEVEE .and. has_levee(i) .and. (.not. is_built_resv(i))) &
               totalvol_bef = totalvol_bef + levsto(i)
#endif

            volwater = volwater + acc_rnof_uc(i)

            IF (.not. is_built_resv(i)) THEN
               IF (DEF_USE_LEVEE .and. has_levee(i)) THEN
                  CALL levee_repartition_storage(i, volwater, wdsrf_ucat(i), fldfrc_levee, &
                     vis_vol_bef_lv, levsto_bef_lv)
                  levee_floodarea(i) = fldfrc_levee * topo_area(i)
#ifdef MaterialTransport
                  CALL mt_rlc_levee_repartition (i, vis_vol_bef_lv, levsto_bef_lv, &
                     volwater, levsto(i), with_input = .true.)
#endif
               ELSE
                  wdsrf_ucat(i) = floodplain_curve(i)%depth (volwater)
                  ! wdsrf_ucat(i) changed: refresh the stored volume read below
                  ! (see volwater_curve; identical to computing it at each read)
                  IF (.not. (DEF_USE_BIFURCATION .or. DEF_USE_LEVEE)) &
                     volwater_curve(i) = floodplain_curve(i)%volume (wdsrf_ucat(i))
               ENDIF
               volwater_ucat(i) = volwater
               IF (wdsrf_ucat(i) > RIVERMIN) THEN
                   veloc_riv(i) = momen_riv(i) / wdsrf_ucat(i)
                ELSE
                   veloc_riv(i) = 0.
               ENDIF
            ELSE
               volresv(ucat2resv(i)) = volwater
               wdsrf_ucat(i) = floodplain_curve(i)%depth (volwater)
            ENDIF

         ENDDO


         ! cells whose downstream cell is a built reservoir, for this call
         IF (DEF_Reservoir_Method > 0) CALL reservoir_upstream_mark (is_built_resv, is_resv_upstream)

         ntimestep = 0
#ifdef CoLMDEBUG
         totaldis  = 0.
#endif

         dt_res(:) = acctime_rnof
         ! Cell that binds the adaptive dt (diagnostics for the substep-count
         ! warning below); reset again at the start of every substep.
         dt_local_min = huge(1._r8)
         i_local_min  = 0
         ! Substep-count limits of this call.  acctime_rnof may exceed
         ! acctime_rnof_max when deltime is longer than the routing interval.
         sub_span = acctime_rnof_max
         IF (ieee_is_finite(acctime_rnof)) sub_span = max(sub_span, acctime_rnof)
         nsub_warn  = max(2, int(min(sub_span / ROUTING_SUBSTEP_WARN_DT,  1.e8_r8)) + 1)
         nsub_abort = max(nsub_warn + 1, int(min(sub_span / ROUTING_SUBSTEP_ABORT_DT, 1.e9_r8)) + 1)

         IF (DEF_USE_BIFURCATION) THEN
            IF (.not. wdsrf_ucat_prev_valid) THEN
               wdsrf_ucat_prev = wdsrf_ucat
               wdsrf_ucat_prev_valid = .true.
            ENDIF
            IF (ieee_is_finite(acctime_rnof)) THEN
               loop_active = acctime_rnof > 0._r8
            ELSE
               loop_active = .true.
            ENDIF
         ELSE
            loop_active = any(dt_res > 0._r8)
         ENDIF

         DO WHILE (loop_active)

            ntimestep = ntimestep + 1

            ! Fail loudly instead of appearing to hang if the adaptive dt
            ! collapses: warn once, then abort, when the substep count implies
            ! a mean dt below ROUTING_SUBSTEP_WARN_DT / ROUTING_SUBSTEP_ABORT_DT
            ! (otherwise the process sits in an effectively unbounded loop
            ! with no output).  The binding cell and its state are those of the
            ! previous, completed substep.  Under bifurcation all workers run
            ! the substeps in lockstep (the dt is per group of connected river
            ! systems, dt_global below is the smallest one), so all workers
            ! reach the limits together and the owner(s) of the binding cell
            ! report it; otherwise each worker counts its own river systems and
            ! reports its own binding cell.
            IF (DEF_USE_BIFURCATION) THEN
               dt_report = routing_last_global_dt
            ELSE
               dt_report = dt_local_min
            ENDIF
            IF (ntimestep == nsub_warn) THEN
               IF (DEF_USE_BIFURCATION .and. p_iam_worker == p_root) THEN
                  write(*,'(A,I0,A,ES12.4,A)') 'WARNING grid_riverlake_flow: substep count ', &
                     ntimestep, ', dt_global=', dt_report, ' s; adaptive routing dt is very small'
               ENDIF
               IF (i_local_min > 0) THEN
                  IF ((.not. DEF_USE_BIFURCATION) .or. &
                      dt_local_min <= routing_last_global_dt * (1._r8 + 1.e-10_r8)) THEN
                     write(*,'(A,I0,A,ES12.4,A,I0,A,I0,A,I0,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4,A,L1)') &
                        'WARNING grid_riverlake_flow: substep ', ntimestep, &
                        ' binding dt=', dt_local_min, ' s; cell rank=', p_iam_glb, &
                        ' ucat=', i_local_min, ' gid=', ucat_gdid(i_local_min), &
                        ' rivlen=', topo_rivlen(i_local_min), ' rivare=', topo_rivare(i_local_min), &
                        ' wdsrf=', wdsrf_ucat(i_local_min), &
                        ' netflux=', sum_hflux_riv(i_local_min), &
                        ' built_dam=', is_built_resv(i_local_min)
                  ENDIF
               ENDIF
               ! Worker stdout is buffered when redirected; a hung job is
               ! usually killed, which would lose this line.
               flush(6)
            ELSEIF (ntimestep > nsub_abort) THEN
               IF (p_iam_worker == p_root .or. .not. DEF_USE_BIFURCATION) &
                  write(*,'(A,I0,A,ES12.4,A)') 'ERROR grid_riverlake_flow: substep count exceeded ', &
                     nsub_abort, ', binding dt=', dt_report, ' s'
               CALL CoLM_stop('grid_riverlake_flow: adaptive routing dt collapsed')
            ENDIF

            CALL worker_push_data (push_next2ucat, wdsrf_ucat, wdsrf_next, fillvalue = spval)
            ! velocity in ocean or inland depression is assumed to be 0.
            CALL worker_push_data (push_next2ucat, veloc_riv,  veloc_next, fillvalue = 0.)

            dt_all(:) = min(dt_res(:), 60.)

            DO i = 1, numucat

               ucatfilter(i) = dt_all(irivsys(i)) > 0

               IF (.not. ucatfilter(i)) CYCLE

               sum_hflux_riv(i) = 0.
               sum_mflux_riv(i) = 0.
               sum_zgrad_riv(i) = 0.

               ! reservoir
               IF (is_built_resv(i)) THEN
                  hflux_fc(i) = 0.
                  mflux_fc(i) = 0.
                  zgrad_dn(i) = 0.
                  CYCLE
               ENDIF

               IF (is_resv_upstream(i)) THEN
                  ! Outflow into a built reservoir: kinematic wave of this cell, not
                  ! an HLL face with the reservoir level, which is the level of the
                  ! whole reservoir volume in one cell and backs water up into the
                  ! cells upstream (see MOD_Grid_ReservoirUpstream).  The reservoir
                  ! does not use momentum and this cell's momentum is replaced by
                  ! the kinematic velocity below, so no pressure term is added.
                  IF (DEF_USE_BIFURCATION .or. DEF_USE_LEVEE) THEN
                     volwater = volwater_ucat(i)
                  ELSE
                     volwater = volwater_curve(i)
                  ENDIF
                  CALL reservoir_upstream_outflow (i, wdsrf_ucat(i), volwater, hflux_fc(i), veloc_kine(i))
                  mflux_fc(i) = hflux_fc(i) * veloc_kine(i)
                  zgrad_dn(i) = 0.

               ELSEIF ((ucat_next(i) > 0) .or. (ucat_next(i) == -9)) THEN

                  IF (ucat_next(i) > 0) THEN
                     ! both rivers are dry.
                     IF ((wdsrf_ucat(i) < RIVERMIN) .and. (wdsrf_next(i) < RIVERMIN)) THEN
                        hflux_fc(i) = 0
                        mflux_fc(i) = 0
                        zgrad_dn(i) = 0
                        CYCLE
                     ENDIF
                  ENDIF

                  ! reconstruction of height of water near interface
                  IF (ucat_next(i) > 0) THEN
                     bedelv_fc = max(topo_rivelv(i), bedelv_next(i))
                     height_up = max(0., wdsrf_ucat(i)+topo_rivelv(i)-bedelv_fc)
                     height_dn = max(0., wdsrf_next(i)+bedelv_next(i)-bedelv_fc)
                  ELSEIF (ucat_next(i) == -9) THEN ! for river mouth
                     bedelv_fc = topo_rivelv(i)
                     height_up = wdsrf_ucat (i)
                     ! sea level is assumed to be 0. and sea bed is assumed to be negative infinity.
                     height_dn = max(0., - bedelv_fc)
                  ENDIF

                  ! velocity at river downstream face (middle region in Riemann problem)
                  veloct_fc = 0.5 * (veloc_riv(i) + veloc_next(i)) &
                     + sqrt(grav * height_up) - sqrt(grav * height_dn)

                  ! height of water at downstream face (middle region in Riemann problem)
                  height_fc = 1/grav * (0.5*(sqrt(grav*height_up) + sqrt(grav*height_dn)) &
                     + 0.25 * (veloc_riv(i) - veloc_next(i))) ** 2

                  IF (height_up > 0) THEN
                     vwave_up = min(veloc_riv(i)-sqrt(grav*height_up), veloct_fc-sqrt(grav*height_fc))
                  ELSE
                     vwave_up = veloc_next(i) - 2.0 * sqrt(grav*height_dn)
                  ENDIF

                  IF (height_dn > 0) THEN
                     vwave_dn = max(veloc_next(i)+sqrt(grav*height_dn), veloct_fc+sqrt(grav*height_fc))
                  ELSE
                     vwave_dn = veloc_riv(i) + 2.0 * sqrt(grav*height_up)
                  ENDIF

                  hflux_up = veloc_riv(i)  * height_up
                  hflux_dn = veloc_next(i) * height_dn
                  mflux_up = veloc_riv(i)**2  * height_up + 0.5*grav * height_up**2
                  mflux_dn = veloc_next(i)**2 * height_dn + 0.5*grav * height_dn**2

                  IF (vwave_up >= 0.) THEN
                     hflux_fc(i) = outletwth(i) * hflux_up
                     mflux_fc(i) = outletwth(i) * mflux_up
                  ELSEIF (vwave_dn <= 0.) THEN
                     hflux_fc(i) = outletwth(i) * hflux_dn
                     mflux_fc(i) = outletwth(i) * mflux_dn
                  ELSE
                     hflux_fc(i) = outletwth(i) * (vwave_dn*hflux_up - vwave_up*hflux_dn &
                        + vwave_up*vwave_dn*(height_dn-height_up)) / (vwave_dn-vwave_up)
                     mflux_fc(i) = outletwth(i) * (vwave_dn*mflux_up - vwave_up*mflux_dn &
                        + vwave_up*vwave_dn*(hflux_dn-hflux_up)) / (vwave_dn-vwave_up)
                  ENDIF

                  sum_zgrad_riv(i) = sum_zgrad_riv(i) + outletwth(i) * 0.5*grav * height_up**2

                  zgrad_dn(i) = outletwth(i) * 0.5*grav * height_dn**2

               ELSEIF (ucat_next(i) == -99) THEN
                  ! downstream is not in model region.
                  ! assume: 1. downstream river bed is equal to this river bed.
                  !         2. downstream water surface is equal to this river depth.
                  !         3. downstream water velocity is equal to this velocity.
                  ! No downstream face exists, so its pressure (zgrad) contribution
                  ! is zero; leaving it unset fed an uninitialised value into
                  ! push_ups2ucat_sum -> sum_zgrad_riv.
                  zgrad_dn(i) = 0.

                  veloc_riv(i) = max(veloc_riv(i), 0.)

                  IF (wdsrf_ucat(i) > topo_rivhgt(i)) THEN

                     ! reconstruction of height of water near interface
                     height_up = wdsrf_ucat (i)
                     height_dn = topo_rivhgt(i)

                     veloct_fc = veloc_riv(i) + sqrt(grav * height_up) - sqrt(grav * height_dn)
                     height_fc = 1/grav * (0.5*(sqrt(grav*height_up) + sqrt(grav*height_dn))) ** 2

                     vwave_up = min(veloc_riv(i)-sqrt(grav*height_up), veloct_fc-sqrt(grav*height_fc))
                     vwave_dn = max(veloc_riv(i)+sqrt(grav*height_dn), veloct_fc+sqrt(grav*height_fc))

                     hflux_up = veloc_riv(i) * height_up
                     hflux_dn = veloc_riv(i) * height_dn
                     mflux_up = veloc_riv(i)**2 * height_up + 0.5*grav * height_up**2
                     mflux_dn = veloc_riv(i)**2 * height_dn + 0.5*grav * height_dn**2

                     IF (vwave_up >= 0.) THEN
                        hflux_fc(i) = outletwth(i) * hflux_up
                        mflux_fc(i) = outletwth(i) * mflux_up
                     ELSEIF (vwave_dn <= 0.) THEN
                        hflux_fc(i) = outletwth(i) * hflux_dn
                        mflux_fc(i) = outletwth(i) * mflux_dn
                     ELSE
                        hflux_fc(i) = outletwth(i) * (vwave_dn*hflux_up - vwave_up*hflux_dn &
                           + vwave_up*vwave_dn*(height_dn-height_up)) / (vwave_dn-vwave_up)
                        mflux_fc(i) = outletwth(i) * (vwave_dn*mflux_up - vwave_up*mflux_dn &
                           + vwave_up*vwave_dn*(hflux_dn-hflux_up)) / (vwave_dn-vwave_up)
                     ENDIF

                     sum_zgrad_riv(i) = sum_zgrad_riv(i) + outletwth(i) * 0.5*grav * height_up**2

                  ELSE
                     hflux_fc(i) = 0
                     mflux_fc(i) = 0
                  ENDIF

               ELSEIF (ucat_next(i) == -10) THEN ! inland depression
                  hflux_fc(i) = 0
                  mflux_fc(i) = 0
                  ! No outflow face; keep the upstream pressure sum deterministic.
                  zgrad_dn(i) = 0.
               ENDIF

               sum_hflux_riv(i) = sum_hflux_riv(i) + hflux_fc(i)
               sum_mflux_riv(i) = sum_mflux_riv(i) + mflux_fc(i)

            ENDDO

            ! the three upstream sums in one exchange (see push_ups2ucat_sum)
            CALL push_ups2ucat_sum (hflux_fc, hflux_sumups, mflux_fc, mflux_sumups, zgrad_dn, zgrad_sumups)

            IF (numucat > 0) THEN
               WHERE (ucatfilter)
                  sum_hflux_riv = sum_hflux_riv - hflux_sumups
                  sum_mflux_riv = sum_mflux_riv - mflux_sumups
                  sum_zgrad_riv = sum_zgrad_riv - zgrad_sumups
               END WHERE
            ENDIF

            ! reservoir operation.
            ! [GridRiverLakeFlow] Bug fix for cascades: the release of a reservoir is
            ! pushed downstream only after all reservoirs are operated, and the outlet
            ! flux of every reservoir was zeroed before the upstream sum, so a
            ! reservoir directly below another never saw that release in the inflow
            ! given to reservoir_operation (the water still arrived; only the
            ! operation was wrong, e.g. a flood release of the upper dam stored by the
            ! lower one).  The release that reached the reservoir in the previous
            ! sub-step (qresv_ups_prev, kept in the restart) is added to that inflow.
            ! A same-step value would need one push per cascade level; sub-steps are
            ! at most 60 s, so the lag is small.  With bifurcation the value is the
            ! release before the outgoing-flux limiter.  Results change only where a
            ! built reservoir is directly downstream of another.
            IF (DEF_Reservoir_Method > 0) THEN

               hflux_resv = 0.
               mflux_resv = 0.

               DO i = 1, numucat

                  IF (.not. ucatfilter(i)) CYCLE

                  IF (ucat_next(i) == -10) THEN
                     IF (is_built_resv(i)) THEN
                        qresv_in (ucat2resv(i)) = - sum_hflux_riv(i) + qresv_ups_prev(ucat2resv(i))
                        qresv_out(ucat2resv(i)) = 0.
                     ENDIF
                     CYCLE
                  ENDIF

                  IF (is_built_resv(i)) THEN

                     irsv = ucat2resv(i)
                     qresv_in(irsv) = - sum_hflux_riv(i) + qresv_ups_prev(irsv)

                     IF (volresv(irsv) > 1.e-4 * volresv_total(irsv) .or. &
                         (DEF_Reservoir_Method == 2 .and. volresv(irsv) > 0._r8)) THEN
                        CALL reservoir_operation (DEF_Reservoir_Method, &
                           irsv, qresv_in(irsv), volresv(irsv), qresv_out(irsv))
                     ELSE
                        qresv_out (irsv) = 0.
                     ENDIF

                     hflux_fc(i) = qresv_out(irsv)
                     mflux_fc(i) = qresv_out(irsv) * sqrt(2*grav*wdsrf_ucat(i))

                     sum_hflux_riv(i) = sum_hflux_riv(i) + hflux_fc(i)
                     sum_mflux_riv(i) = sum_mflux_riv(i) + mflux_fc(i)

                     hflux_resv(i) = hflux_fc(i)
                     mflux_resv(i) = mflux_fc(i)
                  ENDIF

               ENDDO

               CALL push_ups2ucat_sum (hflux_resv, hflux_sumups, mflux_resv, mflux_sumups)

               IF (numucat > 0) THEN
                  WHERE (ucatfilter)
                     sum_hflux_riv = sum_hflux_riv - hflux_sumups
                     sum_mflux_riv = sum_mflux_riv - mflux_sumups
                  END WHERE
               ENDIF

               ! hflux_sumups now holds the release of the reservoirs directly upstream
               DO i = 1, numucat
                  IF (ucatfilter(i) .and. is_built_resv(i)) THEN
                     qresv_ups_prev(ucat2resv(i)) = hflux_sumups(i)
                     ! Dispatch above intentionally uses the lagged upstream
                     ! release; history records the current water balance.
                     qresv_in(ucat2resv(i)) = hflux_fc(i) - sum_hflux_riv(i)
                  ENDIF
               ENDDO

            ENDIF

            dt_local_min = huge(1._r8)
            i_local_min  = 0

            DO i = 1, numucat

               IF (.not. ucatfilter(i)) CYCLE

               dt_this = dt_all(irivsys(i))

               ! constraint 1: CFL condition (only for rivers)
               IF (.not. is_built_resv(i)) THEN
                  IF ((veloc_riv(i) /= 0.) .or. (wdsrf_ucat(i) > 0.)) THEN
                     dt_this = min(dt_this, topo_rivlen(i)/(abs(veloc_riv(i))+sqrt(grav*wdsrf_ucat(i)))*0.8)
                  ENDIF
               ENDIF

               ! constraint 2: Avoid negative values of water
               IF (sum_hflux_riv(i) > 0) THEN
                  IF (.not. is_built_resv(i)) THEN
                     ! for river or lake catchment
                     IF (DEF_USE_BIFURCATION .or. DEF_USE_LEVEE) THEN
                        volwater = volwater_ucat(i)
                     ELSE
                        ! was floodplain_curve(i)%volume(wdsrf_ucat(i)), same value
                        volwater = volwater_curve(i)
                     ENDIF
                  ELSE
                     ! for reservoir
                     volwater = volresv(ucat2resv(i))
                  ENDIF

                  ! With bifurcation, tiny storage is protected by the gross
                  ! outgoing-face limiter below rather than by forcing a
                  ! pathological dt: a nearly dry cell with a small positive
                  ! net outflow would drive the globally synchronized dt
                  ! towards 0 and the substep loop would effectively never
                  ! terminate.  This restores the guard the fork had before
                  ! the SYSU integration merge.  That limiter only runs with
                  ! bifurcation, so otherwise the bound always applies.
                  ! With DEF_GRIDBASED_ROUTING_STORAGE_DT_LIMIT = .false. the limiter
                  ! protects every cell and the bound is skipped: within a CFL step
                  ! the outflow Q*dt = v*A*dt stays below the channel volume A*L, so
                  ! the limiter rarely has to cut a flux.
                  IF ((.not. DEF_USE_BIFURCATION) .or. &
                      (DEF_GRIDBASED_ROUTING_STORAGE_DT_LIMIT .and. volwater > ROUTING_STORAGE_DT_EPS)) &
                     dt_this = min(dt_this, volwater / sum_hflux_riv(i))

               ENDIF

               ! constraint 3: Avoid change of flow direction (only for rivers)
               IF (DEF_GRIDBASED_ROUTING_MOMENTUM_DT_LIMIT .and. (.not. is_built_resv(i)) &
                   .and. (.not. is_resv_upstream(i))) THEN
                  IF (abs(veloc_riv(i)) > 0.1_r8 .and. &
                      veloc_riv(i) * (sum_mflux_riv(i)-sum_zgrad_riv(i)) > 0._r8) THEN
                     dt_this = min(dt_this, &
                        abs(momen_riv(i) * topo_rivare(i) / (sum_mflux_riv(i)-sum_zgrad_riv(i))))
                  ENDIF
               ENDIF

               dt_all(irivsys(i)) = min(dt_this, dt_all(irivsys(i)))

               ! Record the binding cell of this sub-step (diagnostics only).
               IF (dt_this < dt_local_min) THEN
                  dt_local_min = dt_this
                  i_local_min  = i
               ENDIF

            ENDDO

            IF (DEF_USE_BIFURCATION) THEN
               ! Bug fix: synchronize dt BEFORE the outgoing-face limiter below.
               ! The limiter scales the fluxes by volwater / (outflow * dt), and the
               ! dt was the bound of this worker's cells only; the synchronization
               ! (formerly after the limiter) can only shorten it, and the fluxes
               ! already cut with the longer dt stayed cut.  So where the limiter
               ! acted (ordinary_scale < 1) the result depended on how the river
               ! systems were split over the workers.  Results change only in such
               ! cells, and there the flux is cut less; the synchronized dt, and so
               ! the number of substeps, is the same as before.  Without the
               ! momentum bound below this is the same single synchronization,
               ! moved; with it, a second one follows the bound.
               CALL sync_global_routing_dt(dt_res, dt_all, next_loop_active)
               normal_outgoing_rate = 0._r8
               hflux_sumups = 0._r8
               DO i = 1, numucat
                  IF (.not. ucatfilter(i)) CYCLE
                  IF (hflux_fc(i) >= 0._r8) THEN
                     normal_outgoing_rate(i) = hflux_fc(i)
                  ELSE
                     hflux_sumups(i) = -hflux_fc(i)
                  ENDIF
               ENDDO
               CALL push_ups2ucat_sum (hflux_sumups, mflux_sumups)
               WHERE (ucatfilter)
                  normal_outgoing_rate = normal_outgoing_rate + mflux_sumups
               END WHERE
               ordinary_scale = 1._r8
               DO i = 1, numucat
                  IF (.not. ucatfilter(i)) CYCLE
                  IF (normal_outgoing_rate(i) <= 0._r8) CYCLE
                  IF (is_built_resv(i)) THEN
                     volwater = volresv(ucat2resv(i))
                  ELSE
                     volwater = volwater_ucat(i)
                  ENDIF
                  IF (.not. ieee_is_finite(dt_all(irivsys(i)))) THEN
                     ordinary_scale(i) = 0._r8
                  ELSEIF (dt_all(irivsys(i)) <= 0._r8) THEN
                     ordinary_scale(i) = 0._r8
                  ELSE
                     ordinary_scale(i) = min(1._r8, max(volwater, 0._r8) / &
                        (normal_outgoing_rate(i) * dt_all(irivsys(i))))
                  ENDIF
                  normal_outgoing_rate(i) = normal_outgoing_rate(i) * ordinary_scale(i)
               ENDDO
               CALL worker_push_data(push_next2ucat, ordinary_scale, ordinary_scale_next, fillvalue=1._r8)
               DO i = 1, numucat
                  IF (.not. ucatfilter(i)) CYCLE
                  IF (hflux_fc(i) >= 0._r8) THEN
                     dt_this = ordinary_scale(i)
                  ELSE
                     dt_this = ordinary_scale_next(i)
                  ENDIF
                  hflux_fc(i) = hflux_fc(i) * dt_this
                  mflux_fc(i) = mflux_fc(i) * dt_this
                  IF (is_built_resv(i)) qresv_out(ucat2resv(i)) = hflux_fc(i)
               ENDDO
               IF (DEF_Reservoir_Method > 0) THEN
                  ! Cache the actual limited reservoir release for the next dispatch.
                  hflux_resv = 0._r8
                  WHERE (ucatfilter .and. is_built_resv) hflux_resv = hflux_fc
                  CALL push_ups2ucat_sum (hflux_fc, hflux_sumups, mflux_fc, mflux_sumups, &
                     hflux_resv, mflux_resv)
               ELSE
                  CALL push_ups2ucat_sum (hflux_fc, hflux_sumups, mflux_fc, mflux_sumups)
               ENDIF
               WHERE (ucatfilter)
                  sum_hflux_riv = hflux_fc - hflux_sumups
                  sum_mflux_riv = mflux_fc - mflux_sumups
               END WHERE
               DO i = 1, numucat
                  IF (.not. ucatfilter(i)) CYCLE
                  IF (is_built_resv(i)) THEN
                     qresv_in(ucat2resv(i)) = hflux_sumups(i)
                     IF (DEF_Reservoir_Method > 0) qresv_ups_prev(ucat2resv(i)) = mflux_resv(i)
                  ENDIF
               ENDDO
               IF (DEF_GRIDBASED_ROUTING_MOMENTUM_DT_LIMIT) THEN
                  DO i = 1, numucat
                     IF (.not. ucatfilter(i) .or. is_built_resv(i) .or. is_resv_upstream(i)) CYCLE
                     IF (abs(veloc_riv(i)) > 0.1_r8 .and. &
                         veloc_riv(i) * (sum_mflux_riv(i)-sum_zgrad_riv(i)) > 0._r8) THEN
                        dt_this = abs(momen_riv(i) * topo_rivare(i) / &
                           (sum_mflux_riv(i)-sum_zgrad_riv(i)))
                        IF (.not. ieee_is_finite(dt_this) .or. dt_this <= 0._r8) &
                           CALL CoLM_stop('grid_riverlake_flow: invalid final momentum dt')
                        dt_all(irivsys(i)) = min(dt_all(irivsys(i)), dt_this)
                        ! The final momentum bound can bind the step too.
                        IF (dt_this < dt_local_min) THEN
                           dt_local_min = dt_this
                           i_local_min  = i
                        ENDIF
                     ENDIF
                  ENDDO
                  ! the bound may have shortened dt on this worker only
                  CALL sync_global_routing_dt(dt_res, dt_all, next_loop_active)
               ENDIF
               sum_hflux_base = sum_hflux_riv
               IF (allocated(volresv)) volresv_safe = volresv
               CALL bifurcation_calc(wdsrf_ucat, wdsrf_ucat_prev, volwater_ucat, &
                  volwater_ucat_valid, volresv_safe, is_built_resv, dt_all, &
                  irivsys, ucatfilter, normal_outgoing_rate, first_substep = (ntimestep == 1))
               wdsrf_ucat_prev = wdsrf_ucat
               wdsrf_ucat_prev_valid = .true.
               IF (allocated(a_bifflw_lev) .and. allocated(a_bifflw_acctime)) THEN
                  DO ipth = 1, npthout_local
                     i_up = pth_upst_local(ipth)
                     IF (i_up < 1 .or. i_up > numucat) CYCLE
                     IF (.not. ucatfilter(i_up)) CYCLE
                     IF (allocated(bif_path_active)) THEN
                        IF (.not. bif_path_active(ipth)) CYCLE
                     ENDIF
                     a_bifflw_lev(:, ipth) = a_bifflw_lev(:, ipth) &
                        + bif_hflux_lev(:, ipth) * dt_all(irivsys(i_up))
                     a_bifflw_acctime(ipth) = a_bifflw_acctime(ipth) + dt_all(irivsys(i_up))
                  ENDDO
               ENDIF
               WHERE (ucatfilter)
                  sum_hflux_riv = sum_hflux_base + bif_hflux_sum
               END WHERE
            ELSE
#ifdef USEMPI
               IF (rivsys_by_multiple_procs) THEN
                  CALL mpi_allreduce (MPI_IN_PLACE, dt_all, 1, MPI_REAL8, MPI_MIN, p_comm_rivsys, p_err)
               ENDIF
#endif
            ENDIF

#ifdef MaterialTransport
            ! material transport with this substep's water fluxes
            CALL mt_rlc_substep (acctime_rnof, dt_all, irivsys, hflux_fc, sum_hflux_riv, &
               ucatfilter, is_built_resv)
#endif

            DO i = 1, numucat

               IF (.not. ucatfilter(i)) CYCLE

               IF (.not. is_built_resv(i)) THEN
                  IF (DEF_USE_BIFURCATION .or. DEF_USE_LEVEE) THEN
                     volwater = volwater_ucat(i)
                  ELSE
                     ! was floodplain_curve(i)%volume(wdsrf_ucat(i)), same value
                     volwater = volwater_curve(i)
                  ENDIF
               ELSE
                  volwater = volresv(ucat2resv(i))
               ENDIF


               visible_hflux = sum_hflux_riv(i)
               protected_hflux = 0._r8
               IF (DEF_USE_LEVEE .and. DEF_USE_BIFURCATION) THEN
                  IF (has_levee(i) .and. (.not. is_built_resv(i)) .and. allocated(bif_lev_hflux_sum)) THEN
                     protected_hflux = bif_lev_hflux_sum(i)
                     visible_hflux = visible_hflux - protected_hflux
                  ENDIF
               ENDIF
               volwater = volwater - visible_hflux * dt_all(irivsys(i))
               IF (DEF_USE_LEVEE .and. has_levee(i) .and. (.not. is_built_resv(i))) THEN
                  CALL levee_apply_protected_flux(i, protected_hflux, dt_all(irivsys(i)), protected_clip)
                  IF (protected_clip > 0._r8) CALL CoLM_stop('BIF protected-side limiter failed')
               ENDIF
               volwater = max(volwater, 0.)

               ! for inland depression, remove excess water (to be optimized)
               IF (ucat_next(i) == -10) THEN
                  depression_capacity = topo_rivstomax(i)
                  IF (is_built_resv(i)) depression_capacity = volresv_total(ucat2resv(i))
                  IF (volwater > depression_capacity) THEN
                     hflux_fc(i) = (volwater - depression_capacity) / dt_all(irivsys(i))
                     IF (is_built_resv(i)) qresv_out(ucat2resv(i)) = hflux_fc(i)
#ifdef MaterialTransport
                     CALL mt_rlc_inland_depression (i, volwater, depression_capacity, dt_all(irivsys(i)))
#endif
                     volwater = depression_capacity
                  ENDIF
               ENDIF

               IF (DEF_USE_LEVEE .and. has_levee(i) .and. (.not. is_built_resv(i))) THEN
                  CALL levee_repartition_storage(i, volwater, wdsrf_ucat(i), fldfrc_levee, &
                     vis_vol_bef_lv2, levsto_bef_lv2)
                  levee_floodarea(i) = fldfrc_levee * topo_area(i)
#ifdef MaterialTransport
                  CALL mt_rlc_levee_repartition (i, vis_vol_bef_lv2, levsto_bef_lv2, volwater, levsto(i))
#endif
               ELSE
                  wdsrf_ucat(i) = floodplain_curve(i)%depth (volwater)
                  ! wdsrf_ucat(i) changed: refresh the stored volume read below
                  ! (see volwater_curve; identical to computing it at each read)
                  IF (.not. (DEF_USE_BIFURCATION .or. DEF_USE_LEVEE)) &
                     volwater_curve(i) = floodplain_curve(i)%volume (wdsrf_ucat(i))
               ENDIF

               IF (is_built_resv(i)) THEN
                  volresv(ucat2resv(i)) = volwater
               ELSE
                  volwater_ucat(i) = volwater
               ENDIF

               IF ((.not. is_built_resv(i)) .and. (wdsrf_ucat(i) >= RIVERMIN)) THEN
                  friction = grav * topo_rivman(i)**2 / wdsrf_ucat(i)**(7.0/3.0) * abs(momen_riv(i))
                  momen_riv(i) = (momen_riv(i) &
                     - (sum_mflux_riv(i) - sum_zgrad_riv(i)) / topo_rivare(i) * dt_all(irivsys(i))) &
                     / (1 + friction * dt_all(irivsys(i)))
                  veloc_riv(i) = momen_riv(i) / wdsrf_ucat(i)
               ELSE
                  momen_riv(i) = 0
                  veloc_riv(i) = 0
               ENDIF

               ! upstream of a built reservoir: the kinematic channel velocity of the
               ! outflow (D2RIVVEL in CaMa-Flood), as the momentum is not used there
               IF (is_resv_upstream(i) .and. (wdsrf_ucat(i) >= RIVERMIN)) THEN
                  veloc_riv(i) = veloc_kine(i)
                  momen_riv(i) = veloc_riv(i) * wdsrf_ucat(i)
               ENDIF

               ! inland depression river
               IF ((.not. is_built_resv(i)) .and. (ucat_next(i) == -10)) THEN
                  momen_riv(i) = min(0., momen_riv(i))
                  veloc_riv(i) = min(0., veloc_riv(i))
               ENDIF

               ! momen_riv is defined as wdsrf_ucat * veloc_riv.  The clamp
               ! must update it as well: otherwise the next substep reads a
               ! velocity of 20 m/s next to a momentum of 30 (depth 1 m), so
               ! the face flux, the CFL estimate and the friction term use
               ! inconsistent state.  Cells that are not clamped keep exactly
               ! the previous behaviour.
               IF (veloc_riv(i) > 20._r8) THEN
                  veloc_riv(i) = 20._r8
                  momen_riv(i) = wdsrf_ucat(i) * veloc_riv(i)
               ELSEIF (veloc_riv(i) < -20._r8) THEN
                  veloc_riv(i) = -20._r8
                  momen_riv(i) = wdsrf_ucat(i) * veloc_riv(i)
               ENDIF

            ENDDO

#ifdef MaterialTransport
            CALL mt_rlc_diag_accumulate (dt_all, irivsys, ucatfilter, is_built_resv, &
               hflux_fc, levee_floodarea)
#endif

            DO i = 1, numucat
               IF (ucatfilter(i)) THEN

#ifdef CoLMDEBUG
                  IF (ucat_next(i) <= 0) THEN
                     totaldis = totaldis + hflux_fc(i)*dt_all(irivsys(i))
#ifdef MaterialTransport
                     CALL mt_rlc_debug_discharge (i, dt_all(irivsys(i)))
#endif
                  ENDIF
#endif

                  acctime_ucat(i) = acctime_ucat(i) + dt_all(irivsys(i))

                  a_wdsrf_ucat(i) = a_wdsrf_ucat(i) + wdsrf_ucat(i) * dt_all(irivsys(i))
                  a_veloc_riv (i) = a_veloc_riv (i) + veloc_riv (i) * dt_all(irivsys(i))
                  a_discharge (i) = a_discharge (i) + hflux_fc  (i) * dt_all(irivsys(i))

                  ! Bug fix: a leveed cell takes the levee result even when it is
                  ! zero.  levee_fldstg returns fldfrc = levee_frc between bank top
                  ! and levee top, and levee_frc = 0 (all floodplain protected) is
                  ! valid, so zero means "not flooded"; falling back to the
                  ! floodplain curve then put the protected floodplain under water
                  ! in the history (the flood fraction published to the land kept
                  ! the zero).  Only the history floodarea changes, and only for
                  ! leveed cells; reservoirs have no levee (has_levee = .false.)
                  ! and levee_floodarea was 0 for every other cell already.
                  IF (DEF_USE_LEVEE .and. has_levee(i) .and. (.not. is_built_resv(i))) THEN
                     floodarea = levee_floodarea(i)
                  ELSE
                     floodarea = floodplain_curve(i)%floodarea (wdsrf_ucat(i))
                  ENDIF
                  a_floodarea (i) = a_floodarea (i) + floodarea * dt_all(irivsys(i))

                  IF (is_built_resv(i)) THEN
                     volwater = volresv(ucat2resv(i))
                  ELSE
                     IF (DEF_USE_BIFURCATION .or. DEF_USE_LEVEE) THEN
                        volwater = volwater_ucat(i)
                     ELSE
                        ! was floodplain_curve(i)%volume(wdsrf_ucat(i)), same value
                        volwater = volwater_curve(i)
                     ENDIF
                  ENDIF
                  rivsto_hist = min(volwater, floodplain_curve(i)%rivstomax)
                  a_rivsto(i) = a_rivsto(i) + rivsto_hist * dt_all(irivsys(i))
                  a_fldsto(i) = a_fldsto(i) + (volwater - rivsto_hist) * dt_all(irivsys(i))
                  a_flddph(i) = a_flddph(i) + &
                     max(wdsrf_ucat(i) - floodplain_curve(i)%rivhgt, 0._r8) * dt_all(irivsys(i))
                  IF (DEF_USE_LEVEE .and. has_levee(i) .and. (.not. is_built_resv(i))) THEN
                     a_storge(i) = a_storge(i) + (volwater + levsto(i)) * dt_all(irivsys(i))
                     a_levsto(i) = a_levsto(i) + levsto(i) * dt_all(irivsys(i))
                     a_levdph(i) = a_levdph(i) + levdph(i) * dt_all(irivsys(i))
                  ELSE
                     a_storge(i) = a_storge(i) + volwater * dt_all(irivsys(i))
                  ENDIF
                  a_sfcelv(i) = a_sfcelv(i) + &
                     (topo_rivelv(i) + wdsrf_ucat(i)) * dt_all(irivsys(i))

                  IF (DEF_USE_BIFURCATION) &
                     a_bifout(i) = a_bifout(i) + bif_hflux_sum(i) * dt_all(irivsys(i))

                  IF (is_built_resv(i)) THEN
                     irsv = ucat2resv(i)
                     acctime_resv(irsv) = acctime_resv(irsv) + dt_all(irivsys(i))
                     a_volresv   (irsv) = a_volresv  (irsv) + volresv  (irsv) * dt_all(irivsys(i))
                     a_qresv_in  (irsv) = a_qresv_in (irsv) + qresv_in (irsv) * dt_all(irivsys(i))
                     a_qresv_out (irsv) = a_qresv_out(irsv) + qresv_out(irsv) * dt_all(irivsys(i))
                  ENDIF

               ENDIF
            ENDDO

            dt_res = dt_res - dt_all

            IF (DEF_USE_BIFURCATION) THEN
               loop_active = next_loop_active
            ELSE
               loop_active = any(dt_res > 0._r8)
            ENDIF


         ENDDO

#ifdef CoLMDEBUG
         totalvol_aft = 0.
         DO i = 1, numucat
            IF (.not. is_built_resv(i)) THEN
               IF (DEF_USE_BIFURCATION .or. DEF_USE_LEVEE) THEN
                  totalvol_aft = totalvol_aft + volwater_ucat(i)
               ELSE
                  totalvol_aft = totalvol_aft + volwater_curve(i)
               ENDIF
               IF (DEF_USE_LEVEE .and. has_levee(i)) totalvol_aft = totalvol_aft + levsto(i)
            ELSE
               totalvol_aft = totalvol_aft + volresv(ucat2resv(i))
            ENDIF
         ENDDO
#endif
      ENDIF

#ifdef CoLMDEBUG
#ifdef USEMPI
      IF (.not. p_is_worker) ntimestep = 0
      CALL mpi_allreduce (MPI_IN_PLACE, ntimestep, 1, MPI_INTEGER, MPI_MAX, p_comm_glb, p_err)

      IF (.not. p_is_worker) totalvol_bef = 0.
      IF (.not. p_is_worker) totalvol_aft = 0.
      IF (.not. p_is_worker) totalrnof    = 0.
      IF (.not. p_is_worker) totaldis     = 0.

      CALL mpi_allreduce (MPI_IN_PLACE, totalvol_bef, 1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, totalvol_aft, 1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, totalrnof,    1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, totaldis,     1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,'(/,A)') 'Checking River Routing Flow ...'
         write(*,'(A,F12.5,A)') 'River Lake Flow minimum average timestep: ', acctime_rnof/ntimestep, ' seconds'
         write(*,'(A,ES8.1,A)') 'Total water before :  ', totalvol_bef,  ' m^3'
         write(*,'(A,ES8.1,A)') 'Total runoff :        ', totalrnof, ' m^3'
         write(*,'(A,ES8.1,A)') 'Total discharge :     ', totaldis,  ' m^3'
         write(*,'(A,ES8.1,A)') 'Total water change :  ', totalvol_aft-totalvol_bef,  ' m^3'
         write(*,'(A,ES8.1,A)') 'Total water balance : ', totalvol_aft-totalvol_bef-totalrnof+totaldis,  ' m^3'
      ENDIF
#endif

#ifdef CoLMDEBUG
      IF (DEF_GridRiverLake_FloodFeedback) THEN
         flood_vec = (/flood_evap_period, flood_infil_period, flood_land_period/)
#ifdef USEMPI
         CALL mpi_allreduce (MPI_IN_PLACE, flood_vec, 3, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
#endif
         IF (p_is_master .and. any(flood_vec > 0._r8)) THEN
            write(*,'(A,2(1X,ES24.16))') 'Grid flood feedback evap/infil [m3]:', flood_vec(1:2)
            ! the river must lose what land took: relative difference of the two totals
            write(*,'(A,2(1X,ES24.16),A,ES10.3)') 'Grid flood feedback land/river [m3]:', &
               flood_vec(3), flood_vec(1)+flood_vec(2), ' reldiff', &
               abs(flood_vec(3)-flood_vec(1)-flood_vec(2)) / max(flood_vec(3), tiny(1._r8))
         ENDIF
      ENDIF
      flood_land_period = 0._r8
#endif
      flood_evap_period = 0._r8
      flood_infil_period = 0._r8

#ifdef MaterialTransport
      ! material balance (CoLMDEBUG), limiter statistics, reset of the accumulations
      CALL mt_rlc_end_routing ()
#endif

      acctime_rnof = 0.

      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            acc_rnof_uc = 0.
         ENDIF
      ENDIF

      IF (DEF_GridRiverLake_FloodFeedback) CALL publish_flood_feedback(year)

      IF (allocated(is_built_resv)) deallocate(is_built_resv)
      IF (allocated(is_resv_upstream)) deallocate(is_resv_upstream)
      IF (allocated(veloc_kine)) deallocate(veloc_kine)
      IF (allocated(wdsrf_next   )) deallocate(wdsrf_next   )
      IF (allocated(veloc_next   )) deallocate(veloc_next   )
      IF (allocated(hflux_fc     )) deallocate(hflux_fc     )
      IF (allocated(mflux_fc     )) deallocate(mflux_fc     )
      IF (allocated(zgrad_dn     )) deallocate(zgrad_dn     )
      IF (allocated(hflux_resv   )) deallocate(hflux_resv   )
      IF (allocated(mflux_resv   )) deallocate(mflux_resv   )
      IF (allocated(hflux_sumups )) deallocate(hflux_sumups )
      IF (allocated(mflux_sumups )) deallocate(mflux_sumups )
      IF (allocated(zgrad_sumups )) deallocate(zgrad_sumups )
      IF (allocated(sum_hflux_riv)) deallocate(sum_hflux_riv)
      IF (allocated(sum_hflux_base)) deallocate(sum_hflux_base)
      IF (allocated(sum_mflux_riv)) deallocate(sum_mflux_riv)
      IF (allocated(sum_zgrad_riv)) deallocate(sum_zgrad_riv)
      IF (allocated(normal_outgoing_rate)) deallocate(normal_outgoing_rate)
      IF (allocated(ordinary_scale)) deallocate(ordinary_scale)
      IF (allocated(ordinary_scale_next)) deallocate(ordinary_scale_next)
      IF (allocated(volresv_safe)) deallocate(volresv_safe)
      IF (allocated(levee_floodarea)) deallocate(levee_floodarea)
      IF (allocated(volwater_curve )) deallocate(volwater_curve )
      IF (allocated(ucatfilter   )) deallocate(ucatfilter   )
      IF (allocated(dt_res       )) deallocate(dt_res       )
      IF (allocated(dt_all       )) deallocate(dt_all       )

   END SUBROUTINE grid_riverlake_flow

   ! ---------
   SUBROUTINE push_ups2ucat_sum (send1, recv1, send2, recv2, send3, recv3)

   ! Sum over the upstream unit catchments of one to three fields.
   ! Reason: these sums were pushed one field at a time with worker_push_data,
   ! which costs one exchange (one message per neighbouring worker) and one
   ! pass over a dense (largest upstream count) x (unit catchments) table per
   ! field, several times per routing substep.  They are now pushed together in
   ! one batched exchange over compressed rows (MOD_Grid_RiverLakePush).
   ! Results are unchanged, bit for bit: the batch is used only when
   ! cache_ups2ucat%batch_exact says it reproduces worker_push_data (all
   ! upstream weights are 1 in the routing network); each field is then summed
   ! with the same entries in the same order.  Otherwise the fields are pushed
   ! one by one exactly as before.

   USE MOD_Grid_RiverLakePush, only: riverlake_push_field_type, riverlake_push_batch
   IMPLICIT NONE

   real(r8), target :: send1(:), recv1(:)
   real(r8), target, optional :: send2(:), recv2(:), send3(:), recv3(:)

   type(riverlake_push_field_type) :: fields(3)
   integer :: nfield

      IF (cache_ups2ucat%batch_exact) THEN
         fields(1)%send => send1;  fields(1)%recv => recv1;  fields(1)%fillvalue = 0._r8
         nfield = 1
         IF (present(send2)) THEN
            fields(2)%send => send2;  fields(2)%recv => recv2;  fields(2)%fillvalue = 0._r8
            nfield = 2
         ENDIF
         IF (present(send3)) THEN
            fields(3)%send => send3;  fields(3)%recv => recv3;  fields(3)%fillvalue = 0._r8
            nfield = 3
         ENDIF
         CALL riverlake_push_batch (push_ups2ucat, cache_ups2ucat, fields(1:nfield), mode = 'sum')
      ELSE
         CALL worker_push_data (push_ups2ucat, send1, recv1, fillvalue = 0., mode = 'sum')
         IF (present(send2)) &
            CALL worker_push_data (push_ups2ucat, send2, recv2, fillvalue = 0., mode = 'sum')
         IF (present(send3)) &
            CALL worker_push_data (push_ups2ucat, send3, recv3, fillvalue = 0., mode = 'sum')
      ENDIF

   END SUBROUTINE push_ups2ucat_sum


   ! Synchronise the adaptive dt under bifurcation.
   ! dt_all(k) holds the dt that the local stability constraints allow for time
   ! step group k (irivsys).  In bifurcation mode a group is a cluster of river
   ! systems connected by bifurcation pathways (see irivsys_shared in
   ! MOD_Grid_RiverLakeNetwork).  Each group keeps its own dt; a group that lies
   ! on several workers takes the minimum over them, so both ends of every
   ! pathway always use the same dt and the pathway flux stays conservative.
   ! Formerly this routine replaced every dt by the global minimum over all
   ! groups of all workers, which made the whole domain advance with the dt of
   ! its single most restrictive cell.  Changing that changes the time
   ! discretisation only (results differ at the level of the time-stepping
   ! error).  All workers still run the sub-step loop in lockstep, until no
   ! group anywhere has time left: the bifurcation and push routines inside the
   ! loop are collective, so a worker must not leave the loop on its own.
   SUBROUTINE sync_global_routing_dt(dt_res, dt_all, next_loop_active)

   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   IMPLICIT NONE

   real(r8), intent(in) :: dt_res(:)
   real(r8), intent(inout) :: dt_all(:)
   logical, intent(out) :: next_loop_active
   ! dt_reduce(1): smallest dt of the active groups (-huge flags a non-finite
   ! residual); dt_reduce(1+k): dt of shared group k
   real(r8) :: dt_reduce(1+numrivsys_shared)
   logical :: pathological_mask(size(dt_res))
   integer :: i, k

      pathological_mask = .false.
      DO i = 1, size(dt_res)
         IF (.not. ieee_is_finite(dt_res(i))) CYCLE
         IF (dt_res(i) > 0._r8) THEN
            IF (.not. ieee_is_finite(dt_all(i))) THEN
               pathological_mask(i) = .true.
            ELSEIF (dt_all(i) <= 0._r8) THEN
               pathological_mask(i) = .true.
            ENDIF
         ENDIF
      ENDDO

      IF (any(pathological_mask)) THEN
         WHERE (pathological_mask)
            dt_all = min(10._r8, dt_res)
         END WHERE
         IF (p_is_worker .and. p_iam_worker == p_root .and. routing_zero_dt_warn_count < 5) THEN
            routing_zero_dt_warn_count = routing_zero_dt_warn_count + 1
            write(*,'(A)') 'WARNING grid_riverlake_flow: invalid adaptive dt; using bounded fallback.'
         ENDIF
      ENDIF

      dt_reduce(:) = huge(1._r8)
      IF (any(.not. ieee_is_finite(dt_res))) THEN
         dt_reduce(1) = -huge(1._r8)
      ELSE
         DO i = 1, size(dt_res)
            IF (dt_res(i) > 0._r8) THEN
               dt_reduce(1) = min(dt_reduce(1), dt_all(i))
               k = irivsys_shared(i)
               IF (k > 0) dt_reduce(1+k) = min(dt_reduce(1+k), dt_all(i))
            ENDIF
         ENDDO
      ENDIF

#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, dt_reduce, size(dt_reduce), MPI_REAL8, MPI_MIN, &
         p_comm_worker, p_err)
#endif

      IF (.not. ieee_is_finite(dt_reduce(1))) &
         CALL CoLM_stop('grid_riverlake_flow: invalid synchronized routing dt')
      IF (dt_reduce(1) <= -0.5_r8 * huge(1._r8)) &
         CALL CoLM_stop('grid_riverlake_flow: non-finite routing residual')

      routing_last_global_dt = dt_reduce(1)
      next_loop_active = .false.

      DO i = 1, size(dt_res)
         IF (dt_res(i) > 0._r8) THEN
            k = irivsys_shared(i)
            IF (k > 0) dt_all(i) = dt_reduce(1+k)
            IF (dt_all(i) <= 0._r8) &
               CALL CoLM_stop('grid_riverlake_flow: non-positive synchronized routing dt')
            IF (dt_res(i) - dt_all(i) >= dt_res(i)) &
               CALL CoLM_stop('grid_riverlake_flow: routing dt makes no numerical progress')
            IF (dt_res(i) - dt_all(i) > 0._r8) next_loop_active = .true.
         ELSE
            dt_all(i) = 0._r8
         ENDIF
      ENDDO

#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, next_loop_active, 1, MPI_LOGICAL, MPI_LOR, &
         p_comm_worker, p_err)
#endif

   END SUBROUTINE sync_global_routing_dt

   SUBROUTINE publish_flood_feedback(year)

   USE MOD_LandPatch, only: numpatch
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE
   integer, intent(in) :: year
   integer :: i
   real(r8) :: visible, protected, stage, protected_stage, protected_new, fraction
   real(r8), allocatable :: density(:), grid_visible(:), grid_protected(:)
   real(r8), allocatable :: frac_uc(:), frac_grid(:)

      IF (.not. p_is_worker) RETURN
      allocate(density(numucat), grid_visible(numinpm), grid_protected(numinpm), &
         frac_uc(numucat), frac_grid(numinpm))
      flood_visible_uc = 0._r8
      flood_protected_uc = 0._r8
      flood_reservoir_uc = .false.
#ifdef MaterialTransport
      CALL mt_rlc_flood_publish_reset ()
#endif
      frac_uc = 0._r8
      DO i = 1, numucat
         IF (lake_type(i) == 2 .and. allocated(ucat2resv)) THEN
            IF (ucat2resv(i) > 0 .and. allocated(dam_build_year)) THEN
               flood_reservoir_uc(i) = year >= dam_build_year(ucat2resv(i))
            ENDIF
         ENDIF
         visible = max(0._r8, volwater_ucat(i))
         IF (.not. volwater_ucat_valid) visible = max(0._r8, floodplain_curve(i)%volume(wdsrf_ucat(i)))
         IF (flood_reservoir_uc(i)) THEN
            IF (volresv(ucat2resv(i)) == spval) &
               volresv(ucat2resv(i)) = floodplain_curve(i)%volume(wdsrf_ucat(i))
            visible = max(0._r8, volresv(ucat2resv(i)))
         ENDIF
         protected = 0._r8
         IF (DEF_USE_LEVEE .and. has_levee(i)) protected = max(0._r8, levsto(i))
#ifdef MaterialTransport
         CALL mt_rlc_equilibrate_cell (i, visible, DEF_USE_LEVEE .and. has_levee(i))
#endif
         IF (DEF_USE_LEVEE .and. has_levee(i)) THEN
            CALL levee_fldstg(i, visible+protected, stage, protected_new, protected_stage, fraction)
         ELSE
            stage = floodplain_curve(i)%depth(visible)
            fraction = floodplain_curve(i)%floodarea(stage) / max(topo_area(i), 1._r8)
         ENDIF
         IF (fraction <= 0._r8) CYCLE
         flood_visible_uc(i) = max(0._r8, visible - topo_rivstomax(i))
         flood_protected_uc(i) = protected
         frac_uc(i) = min(1._r8, fraction)
#ifdef MaterialTransport
         CALL mt_rlc_flood_publish_cell (i, visible, protected, flood_visible_uc(i))
#endif
      ENDDO

      ! sum_area is allocated only on workers with requests (numucat > 0)
      density = 0._r8
      IF (numucat > 0) THEN
         WHERE (push_inpm2ucat%sum_area > 0._r8) &
            density = flood_visible_uc / max(push_inpm2ucat%sum_area, tiny(1._r8))
      ENDIF
      CALL worker_push_data(push_ucat2inpm, density, grid_visible, fillvalue=0._r8, mode='sum')
      density = 0._r8
      IF (numucat > 0) THEN
         WHERE (push_inpm2ucat%sum_area > 0._r8) &
            density = flood_protected_uc / max(push_inpm2ucat%sum_area, tiny(1._r8))
      ENDIF
      CALL worker_push_data(push_ucat2inpm, density, grid_protected, fillvalue=0._r8, mode='sum')
      ! A zero fraction is a valid dry state and must keep its area in the
      ! average: worker_push_data skips values equal to fillvalue together with
      ! their area, so a missing-value marker distinct from 0 is used here
      ! (with fillvalue = 0, cells half covered by dry unit catchments got the
      ! fraction of the wet part only, i.e. too large a fraction, too small a depth).
      CALL worker_push_data(push_ucat2inpm, frac_uc, frac_grid, &
         fillvalue=RIVERLAKE_FLOOD_MISSING_VALUE, mode='average')
      WHERE (frac_grid == RIVERLAKE_FLOOD_MISSING_VALUE) frac_grid = 0._r8
      WHERE (flood_grid_area > 0._r8)
         frac_grid = min(1._r8, max(0._r8, frac_grid))
         grid_visible = grid_visible / max(flood_grid_area, tiny(1._r8))
         grid_protected = grid_protected / max(flood_grid_area, tiny(1._r8))
      ELSEWHERE
         frac_grid = 0._r8
         grid_visible = 0._r8
         grid_protected = 0._r8
      END WHERE
      IF (numpatch > 0) THEN
         CALL worker_remap_data_grid2pset(remap_patch2inpm, grid_visible, &
            flood_credit_patch, fillvalue=RIVERLAKE_FLOOD_MISSING_VALUE, mode='average')
         WHERE (flood_credit_patch == RIVERLAKE_FLOOD_MISSING_VALUE) flood_credit_patch = 0._r8
         CALL worker_remap_data_grid2pset(remap_patch2inpm, grid_protected, &
            flood_depth_patch, fillvalue=RIVERLAKE_FLOOD_MISSING_VALUE, mode='average')
         WHERE (flood_depth_patch == RIVERLAKE_FLOOD_MISSING_VALUE) flood_depth_patch = 0._r8
         flood_credit_patch = max(0._r8, flood_credit_patch) + max(0._r8, flood_depth_patch)
         CALL worker_remap_data_grid2pset(remap_patch2inpm, frac_grid, &
            flood_fraction_patch, fillvalue=RIVERLAKE_FLOOD_MISSING_VALUE, mode='average')
         WHERE (flood_fraction_patch == RIVERLAKE_FLOOD_MISSING_VALUE) flood_fraction_patch = 0._r8
         flood_fraction_patch = min(1._r8, max(0._r8, flood_fraction_patch))
         flood_depth_patch = 0._r8
         WHERE (flood_fraction_patch > epsilon(1._r8)) &
            flood_depth_patch = 1000._r8 * flood_credit_patch / max(flood_fraction_patch, tiny(1._r8))
      ENDIF
#ifdef MaterialTransport
      ! material of the flood credit, mapped to the patches as the water
      CALL mt_rlc_flood_publish_push (flood_grid_area)
#endif
   END SUBROUTINE publish_flood_feedback


   SUBROUTINE debit_flood_feedback(evap_volume, infil_volume)

   USE MOD_LandPatch, only: numpatch
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   IMPLICIT NONE
   real(r8), intent(out) :: evap_volume, infil_volume
   real(r8), allocatable :: ratio_patch(:), ratio_grid(:), ratio_uc(:), fraction_sum(:)
   real(r8) :: debit_fraction, sink, credit, tol, fraction
   real(r8) :: visible_before, protected_before
   integer :: i, j
#ifdef CoLMDEBUG
   real(r8), allocatable :: sink_patch(:), sink_grid(:)
#endif

      evap_volume = 0._r8
      infil_volume = 0._r8
      IF (.not. p_is_worker) RETURN
      allocate(ratio_patch(numpatch), ratio_grid(numinpm), ratio_uc(numucat), fraction_sum(numucat))
      fraction_sum = 0._r8
      DO i = 1, numpatch
         IF (.not. ieee_is_finite(flood_evap_acc(i)) .or. &
             .not. ieee_is_finite(flood_infil_acc(i))) &
            CALL CoLM_stop('grid flood feedback: nonfinite land sink')
         sink = (flood_evap_acc(i) + flood_infil_acc(i)) * 1.e-3_r8
         credit = flood_credit_patch(i)
         tol = 1.e-10_r8 * max(credit, 1.e-6_r8)
         IF (sink < -tol .or. sink > credit + tol) &
            CALL CoLM_stop('grid flood feedback: land exceeded published patch credit')
      ENDDO
#ifdef CoLMDEBUG
      ! land-side sink volume [m3] of this step, summed over this worker's patches
      allocate(sink_patch(numpatch), sink_grid(numinpm))
      sink_patch = (max(0._r8, flood_evap_acc) + max(0._r8, flood_infil_acc)) * 1.e-3_r8
      sink_grid = 0._r8
      IF (numpatch > 0) CALL worker_remap_data_pset2grid(remap_patch2inpm, sink_patch, sink_grid, &
         fillvalue=0._r8, filter=flood_credit_patch>0._r8)
      flood_land_period = flood_land_period + sum(sink_grid)
      deallocate(sink_patch, sink_grid)
#endif
      DO i = 1, 2
         ratio_patch = 0._r8
         IF (i == 1) THEN
            WHERE (flood_credit_patch > 0._r8) &
               ratio_patch = max(0._r8, flood_evap_acc*1.e-3_r8) / max(flood_credit_patch, tiny(1._r8))
         ELSE
            WHERE (flood_credit_patch > 0._r8) &
               ratio_patch = max(0._r8, flood_infil_acc*1.e-3_r8) / max(flood_credit_patch, tiny(1._r8))
         ENDIF
         CALL worker_remap_data_pset2grid(remap_patch2inpm, ratio_patch, ratio_grid, &
            fillvalue=0._r8, filter=flood_credit_patch>0._r8)
         WHERE (flood_grid_area > 0._r8)
            ratio_grid = ratio_grid / max(flood_grid_area, tiny(1._r8))
         ELSEWHERE
            ratio_grid = 0._r8
         END WHERE
         IF (any(ratio_grid < -1.e-10_r8) .or. any(ratio_grid > 1._r8+1.e-10_r8)) &
            CALL CoLM_stop('grid flood feedback: grid sink exceeded donor credit')
         CALL worker_push_data(push_inpm2ucat, ratio_grid, ratio_uc, fillvalue=0._r8, mode='sum')
         DO j = 1, numucat
            IF (push_inpm2ucat%sum_area(j) > 0._r8) THEN
               debit_fraction = min(1._r8, max(0._r8, ratio_uc(j) / push_inpm2ucat%sum_area(j)))
            ELSE
               debit_fraction = 0._r8
            ENDIF
            IF (i == 1) THEN
               evap_volume = evap_volume + debit_fraction * &
                  (flood_visible_uc(j) + flood_protected_uc(j))
            ELSE
               infil_volume = infil_volume + debit_fraction * &
                  (flood_visible_uc(j) + flood_protected_uc(j))
            ENDIF
            fraction_sum(j) = fraction_sum(j) + debit_fraction
         ENDDO
      ENDDO
      IF (any(fraction_sum > 1._r8+1.e-10_r8)) &
         CALL CoLM_stop('grid flood feedback: donor overdraft')
#ifdef MaterialTransport
      ! take from the river material pools what the land took; check the ledger
      CALL mt_rlc_flood_debit (flood_credit_patch, flood_evap_acc, flood_infil_acc, &
         flood_grid_area, flood_visible_uc, flood_protected_uc)
#endif
      DO j = 1, numucat
         IF (flood_reservoir_uc(j)) THEN
            volresv(ucat2resv(j)) = max(0._r8, volresv(ucat2resv(j)) - fraction_sum(j)*flood_visible_uc(j))
            volwater_ucat(j) = volresv(ucat2resv(j))
         ELSE
            volwater_ucat(j) = max(0._r8, volwater_ucat(j) - fraction_sum(j)*flood_visible_uc(j))
         ENDIF
         IF (DEF_USE_LEVEE .and. has_levee(j)) &
            levsto(j) = max(0._r8, levsto(j) - fraction_sum(j)*flood_protected_uc(j))
#ifdef MaterialTransport
         CALL mt_rlc_equilibrate_cell (j, volwater_ucat(j), DEF_USE_LEVEE .and. has_levee(j))
#endif
         IF (fraction_sum(j) > 0._r8) THEN
            IF (DEF_USE_LEVEE .and. has_levee(j) .and. .not. flood_reservoir_uc(j)) THEN
               CALL levee_repartition_storage(j, volwater_ucat(j), wdsrf_ucat(j), &
                  fraction, visible_before, protected_before)
#ifdef MaterialTransport
               CALL mt_rlc_levee_repartition (j, visible_before, protected_before, &
                  volwater_ucat(j), levsto(j))
#endif
            ELSEIF (flood_reservoir_uc(j)) THEN
               wdsrf_ucat(j) = floodplain_curve(j)%depth(volresv(ucat2resv(j)))
            ELSE
               wdsrf_ucat(j) = floodplain_curve(j)%depth(volwater_ucat(j))
            ENDIF
         ENDIF
      ENDDO
      flood_evap_acc = 0._r8
      flood_infil_acc = 0._r8
   END SUBROUTINE debit_flood_feedback

   ! ---------
   SUBROUTINE grid_riverlake_flow_final ()

      CALL riverlake_network_final ()

      IF (DEF_Reservoir_Method > 0) THEN
         CALL reservoir_final ()
         CALL reservoir_upstream_final ()
      ENDIF


      CALL levee_final()
      CALL bifurcation_final()
#ifdef MaterialTransport
      CALL mt_rlc_final ()
#endif

      IF (allocated(filter_rnof)) deallocate(filter_rnof)
      IF (allocated(flood_depth_patch)) deallocate(flood_depth_patch, flood_fraction_patch, &
         flood_evap_patch, flood_infil_patch)
      IF (allocated(flood_visible_uc)) deallocate(flood_visible_uc, flood_protected_uc, &
         flood_credit_patch, flood_evap_acc, flood_infil_acc, flood_reservoir_uc, flood_grid_area)

   END SUBROUTINE grid_riverlake_flow_final

END MODULE MOD_Grid_RiverLakeFlow
#endif
