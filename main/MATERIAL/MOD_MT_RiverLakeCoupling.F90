#include <define.h>

#if (defined GridRiverLakeFlow) && (defined MaterialTransport)
MODULE MOD_MT_RiverLakeCoupling
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   Materials in the grid river-lake flow: the calls made from
!   MOD_Grid_RiverLakeFlow (each one marked "#ifdef MaterialTransport" there).
!
!   The routines read the host water state (passed in, or from the river
!   modules) and change only the material pools of MOD_MT_RiverLake and the
!   flood-feedback material arrays below.  No host variable is changed, so
!   the host solution is the same as without them.
!
!   Units: river pools are R*m^3 (ratio times water volume); land-side flood
!   arrays are R*mm over the patch, as the land material pools.
!
!   Materials with their own river state (sediment) are reached through the
!   route hooks of MOD_MT_Lifecycle, called from the same routines.
!-------------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_MT_Defs, only: nmts, materials, mt_uses_generic_water, material_is_nonvolatile_solute
   USE MOD_MT_Solute, only: mt_has_dissolved_limit, mt_equilibrate_dissolved
   USE MOD_MT_RiverLake
   USE MOD_MT_Lifecycle, only: mt_lifecycle_route_has_active, mt_lifecycle_route_init, &
      mt_lifecycle_route_read_restart, mt_lifecycle_route_forcing_put, &
      mt_lifecycle_route_diag_accumulate, mt_lifecycle_route_calc, &
      mt_lifecycle_route_sediment_bif_accumulate, mt_lifecycle_route_sediment_levee_repartition, &
      mt_lifecycle_route_final, material_lifecycle_has_levee_flood_publisher, &
      material_lifecycle_has_flood_publisher, mt_lifecycle_publish_levee_flood_patch, &
      mt_lifecycle_publish_flood_patch
   IMPLICIT NONE

   ! ----- grid flood feedback -----
   ! Land side (nmts, numpatch), read and written in CoLMMAIN:
   !   flood_mt_credit_patch : material of the flood water published to the patch
   !   flood_mt_evap_patch   : material lost (>0) or gained (<0) by flood evaporation
   !   flood_mt_land_patch   : material that entered the soil with flood infiltration
   real(r8), allocatable :: flood_mt_credit_patch(:,:)
   real(r8), allocatable :: flood_mt_evap_patch  (:,:)
   real(r8), allocatable :: flood_mt_land_patch  (:,:)
   ! River side (nmts, numucat): material of the published visible and
   ! protected (levee) flood water of each unit catchment
   real(r8), allocatable :: flood_visible_mt_uc  (:,:)
   real(r8), allocatable :: flood_protected_mt_uc(:,:)
   integer, save :: flood_mt_ledger_reports = 0

   ! reservoir state passed to the transport (empty when there are no reservoirs)
   real(r8), allocatable :: volresv_mt  (:)
   integer,  allocatable :: ucat2resv_mt(:)

   integer, save :: lim_diag_printed = 0

   ! water carrier of the provider materials in one substep (numucat):
   ! storage at the start (before the host update) and at the end, the levee
   ! protected part of both, the flooded and the protected wet area
   real(r8), allocatable :: prov_storage_start  (:)
   real(r8), allocatable :: prov_protected_start(:)
   real(r8), allocatable :: prov_storage_end    (:)
   real(r8), allocatable :: prov_protected_end  (:)
   real(r8), allocatable :: prov_floodarea      (:)
   real(r8), allocatable :: prov_protected_area (:)

   ! flood of the last substep of a routing call, for land providers that
   ! read the river flood (methane inundation scheme 7), numucat:
   ! flooded area as in the river history, its water depth, and the flooded
   ! area of leveed cells
   real(r8), allocatable :: pub_floodarea  (:)
   real(r8), allocatable :: pub_flooddepth (:)
   real(r8), allocatable :: pub_leveearea  (:)

#ifdef CoLMDEBUG
   ! material balance of one routing call (cf. the water balance printed there)
   real(r8), allocatable :: mt_mass_bef(:), mt_mass_inp(:), mt_mass_dis(:)
#endif

   PUBLIC :: flood_mt_credit_patch, flood_mt_evap_patch, flood_mt_land_patch

CONTAINS

   SUBROUTINE mt_rlc_write_history (file_hist_ucat, itime_in_file_ucat, acctime_ucat)
      USE MOD_MT_Lifecycle, only: mt_lifecycle_route_write_history
      character(len=*), intent(in) :: file_hist_ucat
      integer, intent(in) :: itime_in_file_ucat
      real(r8), intent(in) :: acctime_ucat(:)

      CALL mt_riverlake_write_history (file_hist_ucat, itime_in_file_ucat, acctime_ucat)
      CALL mt_lifecycle_route_write_history (file_hist_ucat, itime_in_file_ucat)
   END SUBROUTINE mt_rlc_write_history

   SUBROUTINE mt_rlc_flush_history ()
      USE MOD_MT_Lifecycle, only: mt_lifecycle_route_flush_history

      CALL mt_riverlake_flush_acc ()
      CALL mt_lifecycle_route_flush_history ()
   END SUBROUTINE mt_rlc_flush_history

   SUBROUTINE mt_rlc_land_flood_evap (ipatch, patchtype, water_credit, &
         evaporation, deltim, temperature, forc_q, forc_psrf, forc_us, forc_vs)

      USE MOD_Namelist, only: LWINFILT => DEF_GridRiverLake_FloodFeedback
      USE MOD_MT_Evapo, only: mt_flood_evap_loss
      integer, intent(in) :: ipatch, patchtype
      real(r8), allocatable, intent(in) :: water_credit(:) ! published patch water [m]
      real(r8), intent(in) :: evaporation ! water rate [mm/s]
      real(r8), intent(in) :: deltim, temperature, forc_q, forc_psrf, forc_us, forc_vs

      ! Vapour exchange is debited from the river pools at the next routing call.
      IF (allocated(flood_mt_evap_patch)) THEN
         flood_mt_evap_patch(:,ipatch) = 0._r8
         IF (LWINFILT .and. patchtype == 0 .and. evaporation > 0._r8) &
            CALL mt_flood_evap_loss(ipatch, water_credit(ipatch)*1000._r8, &
               flood_mt_credit_patch(:,ipatch), evaporation*deltim, temperature, &
               forc_q, forc_psrf, forc_us, forc_vs, flood_mt_evap_patch(:,ipatch))
      ENDIF
   END SUBROUTINE mt_rlc_land_flood_evap

   SUBROUTINE mt_rlc_land_flood_input (ipatch, patchtype, water_credit, &
         evaporation, infiltration, deltim, material_input)

      USE MOD_Namelist, only: LWINFILT => DEF_GridRiverLake_FloodFeedback
      integer, intent(in) :: ipatch, patchtype
      real(r8), allocatable, intent(in) :: water_credit(:) ! published patch water [m]
      real(r8), intent(in) :: evaporation, infiltration ! water rates [mm/s]
      real(r8), intent(in) :: deltim
      real(r8), intent(out) :: material_input(nmts)
      integer :: imt

      ! Infiltration carries its share of the credit remaining after evaporation.
      material_input = 0._r8
      IF (LWINFILT .and. patchtype == 0 .and. infiltration > 0._r8 .and. &
          allocated(flood_mt_credit_patch)) THEN
         IF (water_credit(ipatch)*1000._r8 <= evaporation*deltim) THEN
            IF (infiltration*deltim > 1.e-12_r8*max(water_credit(ipatch)*1000._r8, 1._r8)) &
               CALL CoLM_stop('grid flood feedback: infiltration without published material credit')
         ELSE
            material_input = (flood_mt_credit_patch(:,ipatch) - &
               flood_mt_evap_patch(:,ipatch)) * &
               (infiltration*deltim/(water_credit(ipatch)*1000._r8-evaporation*deltim))
            DO imt = 1, nmts
               IF (.not. mt_has_dissolved_limit(imt)) CYCLE
               material_input(imt) = min(max(material_input(imt),0._r8), &
                  materials(imt)%max_dissolved_conc*infiltration*deltim)
            ENDDO
         ENDIF
      ENDIF
   END SUBROUTINE mt_rlc_land_flood_input

   SUBROUTINE mt_rlc_land_publish_input (ipatch, patchtype, dynamic_wetland, material_input)
      integer, intent(in) :: ipatch, patchtype
      logical, intent(in) :: dynamic_wetland
      real(r8), intent(in) :: material_input(:)

      IF (.not. (patchtype == 2 .and. .not. dynamic_wetland)) THEN
         IF (allocated(flood_mt_land_patch)) flood_mt_land_patch(:,ipatch) = material_input
      ENDIF
   END SUBROUTINE mt_rlc_land_publish_input

   ! ---------
   SUBROUTINE mt_rlc_init (start_year, restart_file)

   ! River material pools: from the river restart file, else (or for the
   ! materials missing in it) from the water state at the default ratio.
   ! Called after the river water state is complete (rebuild_volwater_ucat).

   USE MOD_Grid_RiverLakeNetwork,  only: lake_type
   USE MOD_Grid_RiverLakeTimeVars, only: wdsrf_ucat, volresv
   USE MOD_Grid_Reservoir,         only: ucat2resv, dam_build_year
   IMPLICIT NONE

   integer,          intent(in) :: start_year
   character(len=*), intent(in) :: restart_file

   integer :: i, irsv
   logical :: mt_restart_found
   logical,  allocatable :: mt_missing(:), is_built_resv_init(:)
   real(r8), allocatable :: wdsrf_safe(:), volresv_safe(:)
   integer,  allocatable :: ucat2resv_safe(:)

      CALL mt_riverlake_init ()
      IF (nmts <= 0) RETURN

      CALL mt_lifecycle_route_init ()
      IF (len_trim(restart_file) > 0) CALL mt_lifecycle_route_read_restart (restart_file)

      allocate (mt_missing(nmts))
      mt_missing = .true.
      mt_restart_found = .false.
      IF (len_trim(restart_file) > 0) THEN
         CALL mt_riverlake_read_restart (restart_file, mt_restart_found, mt_missing)
      ENDIF

      IF (any(mt_missing)) THEN
         IF (allocated(wdsrf_ucat)) THEN
            allocate (wdsrf_safe(size(wdsrf_ucat)));    wdsrf_safe = wdsrf_ucat
         ELSE
            allocate (wdsrf_safe(0))
         ENDIF
         IF (allocated(volresv)) THEN
            allocate (volresv_safe(size(volresv)));     volresv_safe = volresv
         ELSE
            allocate (volresv_safe(0))
         ENDIF
         IF (allocated(ucat2resv)) THEN
            allocate (ucat2resv_safe(size(ucat2resv))); ucat2resv_safe = ucat2resv
         ELSE
            allocate (ucat2resv_safe(0))
         ENDIF
         ! reservoirs already built in the start year hold their own volume
         allocate (is_built_resv_init(size(wdsrf_safe)))
         is_built_resv_init = .false.
         IF (allocated(lake_type) .and. allocated(dam_build_year)) THEN
            DO i = 1, min(size(is_built_resv_init), size(lake_type), size(ucat2resv_safe))
               IF (lake_type(i) /= 2) CYCLE
               irsv = ucat2resv_safe(i)
               IF (irsv < 1 .or. irsv > size(dam_build_year)) CYCLE
               is_built_resv_init(i) = start_year >= dam_build_year(irsv)
            ENDDO
         ENDIF
         CALL mt_riverlake_init_from_water (wdsrf_safe, volresv_safe, ucat2resv_safe, &
            mt_missing, is_built_resv_init)
         deallocate (wdsrf_safe, volresv_safe, ucat2resv_safe, is_built_resv_init)
      ENDIF

      deallocate (mt_missing)

   END SUBROUTINE mt_rlc_init

   ! ---------
   SUBROUTINE mt_rlc_runoff_input (filter_rnof, rnof_uc, deltime)

   ! Material of the land runoff of this land step, mapped to the unit
   ! catchments as the runoff water is (patch -> input grid -> unit catchment),
   ! into the accumulation that is merged at the next routing call; and the
   ! precipitation forcing of the provider materials (sediment yield).
   ! Called on all workers (the pushes are collective).

   USE MOD_LandPatch,             only: numpatch
   USE MOD_Grid_RiverLakeNetwork, only: numinpm, numucat, remap_patch2inpm, &
      push_inpm2ucat, push_ucat2inpm
   USE MOD_WorkerPushData
   USE MOD_MT_Vars,               only: mt_rnof_step
   IMPLICIT NONE

   logical,  intent(in) :: filter_rnof(:)   ! patches whose runoff enters the rivers
   real(r8), intent(in) :: rnof_uc(:)       ! runoff water [mm/s] per unit catchment
   real(r8), intent(in) :: deltime          ! land time step [s]

   integer :: imt
   real(r8), allocatable :: mt_rnof_gd(:,:), mt_rnof_uc(:,:)

      IF (.not. p_is_worker) RETURN
      IF (nmts <= 0) RETURN

      allocate (mt_rnof_gd (nmts, numinpm))
      allocate (mt_rnof_uc (nmts, numucat))
      mt_rnof_gd = 0._r8
      mt_rnof_uc = 0._r8

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         IF (numpatch > 0) &
            CALL worker_remap_data_pset2grid (remap_patch2inpm, mt_rnof_step(imt,:), &
               mt_rnof_gd(imt,:), fillvalue = 0._r8, filter = filter_rnof)
         IF (numinpm > 0) THEN
            WHERE (push_ucat2inpm%sum_area > 0._r8)
               mt_rnof_gd(imt,:) = mt_rnof_gd(imt,:) / push_ucat2inpm%sum_area
            END WHERE
         ENDIF
         CALL worker_push_data (push_inpm2ucat, mt_rnof_gd(imt,:), mt_rnof_uc(imt,:), &
            fillvalue = 0._r8, mode = 'sum')
      ENDDO

      ! water reference [m^3] as acc_rnof_uc; material [R*mm*m^2 -> R*m^3]
      IF (numucat > 0) &
         CALL mt_riverlake_input_from_runoff (rnof_uc*1.e-3*deltime, numucat, mt_rnof_uc*1.e-3)

      deallocate (mt_rnof_gd, mt_rnof_uc)

      IF (mt_lifecycle_route_has_active()) CALL provider_precip_input (filter_rnof, deltime)

   END SUBROUTINE mt_rlc_runoff_input

   ! ---------
   SUBROUTINE provider_precip_input (filter_rnof, deltime)

   ! Precipitation of this land step per unit catchment [mm/s], averaged over
   ! the patches with valid forcing, and the fraction of the catchment area
   ! those patches cover; mapped as the runoff (patch -> input grid -> unit
   ! catchment).  Workers.

   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   USE MOD_Vars_Global,           only: spval
   USE MOD_LandPatch,             only: numpatch
   USE MOD_Vars_1DForcing,        only: forc_prc, forc_prl
   USE MOD_Grid_RiverLakeNetwork, only: numinpm, numucat, topo_area, remap_patch2inpm, &
      push_inpm2ucat, push_ucat2inpm
   USE MOD_WorkerPushData
   IMPLICIT NONE

   logical,  intent(in) :: filter_rnof(:)   ! patches whose runoff enters the rivers
   real(r8), intent(in) :: deltime          ! land time step [s]

   integer :: i
   logical,  allocatable :: filter_prcp(:)
   real(r8), allocatable :: prcp_pch(:), prcp_gd(:), prcp_area_gd(:), prcp_uc(:), prcp_area_uc(:)

      allocate (prcp_pch(numpatch), filter_prcp(numpatch))
      allocate (prcp_gd(numinpm), prcp_area_gd(numinpm))
      allocate (prcp_uc(numucat), prcp_area_uc(numucat))
      prcp_pch     = 0._r8
      prcp_gd      = 0._r8
      prcp_area_gd = 0._r8

      IF (numpatch > 0) THEN
         filter_prcp = filter_rnof
         DO i = 1, numpatch
            IF (.not. filter_prcp(i)) CYCLE
            filter_prcp(i) = ieee_is_finite(forc_prc(i)) .and. ieee_is_finite(forc_prl(i))
            IF (filter_prcp(i)) filter_prcp(i) = forc_prc(i) /= spval .and. forc_prl(i) /= spval
            IF (filter_prcp(i)) filter_prcp(i) = forc_prc(i) >= 0._r8 .and. forc_prl(i) >= 0._r8
         ENDDO
         WHERE (filter_prcp) prcp_pch = forc_prc + forc_prl
         CALL worker_remap_data_pset2grid (remap_patch2inpm, prcp_pch, prcp_gd, &
            fillvalue = 0._r8, filter = filter_prcp)
         prcp_pch = 1._r8
         CALL worker_remap_data_pset2grid (remap_patch2inpm, prcp_pch, prcp_area_gd, &
            fillvalue = 0._r8, filter = filter_prcp)
      ENDIF

      IF (numinpm > 0) THEN
         WHERE (push_ucat2inpm%sum_area > 0._r8)
            prcp_gd      = prcp_gd      / push_ucat2inpm%sum_area
            prcp_area_gd = prcp_area_gd / push_ucat2inpm%sum_area
         ELSEWHERE
            prcp_gd      = 0._r8
            prcp_area_gd = 0._r8
         END WHERE
      ENDIF

      CALL worker_push_data (push_inpm2ucat, prcp_gd, prcp_uc, fillvalue = 0._r8, mode = 'sum')
      CALL worker_push_data (push_inpm2ucat, prcp_area_gd, prcp_area_uc, fillvalue = 0._r8, mode = 'sum')

      IF (numucat > 0) THEN
         WHERE (prcp_area_uc > 0._r8)
            prcp_uc = prcp_uc / prcp_area_uc
         ELSEWHERE
            prcp_uc = 0._r8
         END WHERE
         WHERE (topo_area > 0._r8)
            prcp_area_uc = prcp_area_uc / topo_area
         ELSEWHERE
            prcp_area_uc = 0._r8
         END WHERE
      ENDIF

      CALL mt_lifecycle_route_forcing_put (prcp_uc, deltime, prcp_area_uc)

      deallocate (prcp_pch, filter_prcp, prcp_gd, prcp_area_gd, prcp_uc, prcp_area_uc)

   END SUBROUTINE provider_precip_input

   ! ---------
   SUBROUTINE mt_rlc_begin_routing ()

   ! Start of a routing call (workers): the runoff material accumulated over
   ! the land steps joins the input buffer, as acc_rnof_uc joins the water.

   USE MOD_Grid_RiverLakeNetwork,  only: numucat
   USE MOD_Grid_RiverLakeTimeVars, only: volresv
   USE MOD_Grid_Reservoir,         only: ucat2resv
   IMPLICIT NONE

   integer :: imt

      IF (.not. p_is_worker) RETURN
      IF (nmts <= 0) RETURN

      IF (allocated(volresv_mt)) deallocate (volresv_mt)
      IF (allocated(ucat2resv_mt)) deallocate (ucat2resv_mt)
      IF (allocated(volresv)) THEN
         allocate (volresv_mt(size(volresv)))
         volresv_mt = volresv
      ELSE
         allocate (volresv_mt(0))
      ENDIF
      IF (allocated(ucat2resv)) THEN
         allocate (ucat2resv_mt(size(ucat2resv)))
         ucat2resv_mt = ucat2resv
      ELSE
         allocate (ucat2resv_mt(0))
      ENDIF

      IF (allocated(mt_decay_source)) mt_decay_source = 0._r8

#ifdef CoLMDEBUG
      CALL ensure_debug_arrays ()
      IF (numucat > 0) THEN
         DO imt = 1, nmts
            IF (.not. mt_uses_generic_water(imt)) CYCLE
            mt_mass_bef(imt) = pool_total(imt)
            mt_mass_inp(imt) = sum(acc_mt_inp(imt,:))
         ENDDO
      ENDIF
#endif

      IF (numucat > 0) THEN
         DO imt = 1, nmts
            IF (.not. mt_uses_generic_water(imt)) CYCLE
            mt_inp_buf(imt,:) = mt_inp_buf(imt,:) + acc_mt_inp(imt,:)
            acc_mt_inp(imt,:) = 0._r8
         ENDDO
      ENDIF

   END SUBROUTINE mt_rlc_begin_routing

   ! ---------
   SUBROUTINE mt_rlc_levee_repartition (i, vis_bef, lev_bef, vis_aft, lev_aft, with_input)

   ! The host has moved water between the visible and the protected (levee)
   ! side of unit catchment i; move the material with it.  with_input: the
   ! water of this call's runoff is already in vis_bef, its material still
   ! in the input buffer (before the first substep).

   IMPLICIT NONE
   integer,  intent(in) :: i
   real(r8), intent(in) :: vis_bef, lev_bef, vis_aft, lev_aft   ! [m^3]
   logical,  intent(in), optional :: with_input

      IF (nmts <= 0) RETURN
      CALL mt_lifecycle_route_sediment_levee_repartition (i, vis_bef, lev_bef, vis_aft, lev_aft)
      IF (present(with_input)) THEN
         IF (with_input) THEN
            CALL mt_riverlake_levee_repartition (i, vis_bef, lev_bef, vis_aft, lev_aft, &
               pending_mt_pool = mt_inp_buf(:,i))
            RETURN
         ENDIF
      ENDIF
      CALL mt_riverlake_levee_repartition (i, vis_bef, lev_bef, vis_aft, lev_aft)

   END SUBROUTINE mt_rlc_levee_repartition

   ! ---------
   SUBROUTINE mt_rlc_equilibrate_cell (i, visible, with_levee)

   ! Dissolution equilibrium of unit catchment i for the given visible water and
   ! the protected water of its levee (if any), before the host takes water out.

   USE MOD_Grid_RiverLakeLevee, only: levsto
   IMPLICIT NONE
   integer,  intent(in) :: i
   real(r8), intent(in) :: visible      ! [m^3]
   logical,  intent(in) :: with_levee

   real(r8) :: protected

      IF (nmts <= 0) RETURN
      protected = 0._r8
      IF (with_levee) protected = levsto(i)
      CALL mt_riverlake_equilibrate_cell (i, visible, protected)

   END SUBROUTINE mt_rlc_equilibrate_cell

   ! ---------
   SUBROUTINE mt_rlc_substep (dt_ref, dt_all, irivsys, hflux_fc, sum_hflux_riv, ucatfilter, is_built_resv)

   ! Material transport of one routing substep, with the water fluxes of the
   ! substep, before the host updates the storages with them.  For the provider
   ! materials: the bifurcation fluxes, and the water storage at the start.

   USE MOD_Grid_RiverLakeNetwork,     only: numucat, npthout_local, floodplain_curve
   USE MOD_Grid_RiverLakeTimeVars,    only: wdsrf_ucat, volresv, volwater_ucat
   USE MOD_Grid_RiverLakeBifurcation, only: bif_hflux_lev
   USE MOD_Grid_RiverLakeLevee,       only: has_levee, levsto
   IMPLICIT NONE

   real(r8), intent(in) :: dt_ref
   real(r8), intent(in) :: dt_all(:)
   integer,  intent(in) :: irivsys(:)
   real(r8), intent(in) :: hflux_fc(:)
   real(r8), intent(in) :: sum_hflux_riv(:)
   logical,  intent(in) :: ucatfilter(:)
   logical,  intent(in) :: is_built_resv(:)

   integer  :: i
   real(r8) :: volwater

      IF (.not. p_is_worker) RETURN
      IF (nmts <= 0) RETURN

      IF (DEF_USE_BIFURCATION) &
         CALL mt_lifecycle_route_sediment_bif_accumulate (dt_all, irivsys, ucatfilter, bif_hflux_lev)

      IF (mt_lifecycle_route_has_active()) THEN
         CALL ensure_provider_arrays (numucat)
         prov_storage_start   = 0._r8
         prov_protected_start = 0._r8
         DO i = 1, numucat
            IF (.not. ucatfilter(i)) CYCLE
            ! the water the host starts the substep from (cf. its update loop)
            IF (is_built_resv(i)) THEN
               volwater = volresv(ucat2resv_mt(i))
            ELSEIF (DEF_USE_BIFURCATION .or. DEF_USE_LEVEE) THEN
               volwater = volwater_ucat(i)
            ELSE
               volwater = floodplain_curve(i)%volume (wdsrf_ucat(i))
            ENDIF
            prov_storage_start(i) = max(volwater, 0._r8)
            IF (DEF_USE_LEVEE .and. has_levee(i) .and. (.not. is_built_resv(i))) THEN
               prov_storage_start(i)   = prov_storage_start(i) + max(levsto(i), 0._r8)
               prov_protected_start(i) = max(levsto(i), 0._r8)
            ENDIF
         ENDDO
      ENDIF

      IF (allocated(volresv)) volresv_mt = volresv
      IF (DEF_USE_BIFURCATION) THEN
         CALL mt_riverlake_substep (dt_ref, dt_all, irivsys, hflux_fc, sum_hflux_riv, &
            wdsrf_ucat, ucatfilter, volresv_mt, ucat2resv_mt, is_built_resv, &
            do_bif = .true., bif_hflux_lev_in = bif_hflux_lev, npthout_local_in = npthout_local)
      ELSE
         CALL mt_riverlake_substep (dt_ref, dt_all, irivsys, hflux_fc, sum_hflux_riv, &
            wdsrf_ucat, ucatfilter, volresv_mt, ucat2resv_mt, is_built_resv)
      ENDIF

   END SUBROUTINE mt_rlc_substep

   ! ---------
   SUBROUTINE mt_rlc_inland_depression (i, volwater, volmax, dt)

   ! The host removes the water above volmax from an inland depression; the
   ! material leaves with it at the cell's dissolved ratio.

   IMPLICIT NONE

   integer,  intent(in) :: i
   real(r8), intent(in) :: volwater   ! water before the removal [m^3]
   real(r8), intent(in) :: volmax     ! water kept [m^3]
   real(r8), intent(in) :: dt         ! substep [s]

   integer  :: imt
   real(r8) :: frac_remove, mt_removed

      IF (nmts <= 0) RETURN
      IF (volwater <= 1.e-6_r8) RETURN

      frac_remove = (volwater - volmax) / volwater
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         IF (allocated(mt_solid)) THEN
            IF (mt_has_dissolved_limit(imt)) &
               CALL mt_equilibrate_dissolved (imt, volwater, mt_mass(imt,i), mt_solid(imt,i))
         ENDIF
         mt_removed = mt_mass(imt,i) * frac_remove
         mt_mass(imt,i) = mt_mass(imt,i) - mt_removed
         mt_flux_out(imt,i) = mt_removed / dt
      ENDDO

   END SUBROUTINE mt_rlc_inland_depression

   ! ---------
   SUBROUTINE mt_rlc_diag_accumulate (dt_all, irivsys, ucatfilter, is_built_resv, &
         hflux_fc, levee_floodarea)

   ! History accumulation of the substep, after the host storage update; and
   ! the substep's water state for the provider materials.

   USE MOD_Grid_RiverLakeNetwork,  only: numucat, topo_area, floodplain_curve, levee_frc_data
   USE MOD_Grid_RiverLakeTimeVars, only: wdsrf_ucat, volresv, volwater_ucat, veloc_riv
   USE MOD_Grid_RiverLakeBifurcation, only: bif_hflux_lev
   USE MOD_Grid_RiverLakeLevee,    only: has_levee, levsto, levdph
   IMPLICIT NONE

   real(r8), intent(in) :: dt_all(:)
   integer,  intent(in) :: irivsys(:)
   logical,  intent(in) :: ucatfilter(:)
   logical,  intent(in) :: is_built_resv(:)
   real(r8), intent(in) :: hflux_fc(:)          ! downstream face flux [m^3/s]
   real(r8), intent(in) :: levee_floodarea(:)   ! flooded area of leveed cells [m^2]

   integer  :: i
   logical  :: leveed

      IF (.not. p_is_worker) RETURN
      IF (nmts <= 0) RETURN
      IF (numucat <= 0 .and. .not. DEF_USE_BIFURCATION) RETURN

      IF (allocated(volresv)) volresv_mt = volresv
      IF (DEF_USE_BIFURCATION) THEN
         CALL mt_riverlake_diag_accumulate (dt_all, irivsys, ucatfilter, wdsrf_ucat, &
         volresv_mt, ucat2resv_mt, is_built_resv, hflux_fc_in = hflux_fc, &
         allow_reverse_flow = .true., bif_hflux_lev_in = bif_hflux_lev)
      ELSE
         CALL mt_riverlake_diag_accumulate (dt_all, irivsys, ucatfilter, wdsrf_ucat, &
         volresv_mt, ucat2resv_mt, is_built_resv, hflux_fc_in = hflux_fc, &
         allow_reverse_flow = .true.)
      ENDIF
      IF (numucat <= 0) RETURN

      IF (material_lifecycle_has_flood_publisher() .or. material_lifecycle_has_levee_flood_publisher()) THEN
         ! CoLM202X total_floodarea / total_flooddepth / levee_floodarea of
         ! the last substep.  The area is the one of the river history
         ! (a leveed cell takes the levee result even when it is zero).
         CALL ensure_publish_arrays (numucat)
         DO i = 1, numucat
            IF (.not. ucatfilter(i)) CYCLE
            leveed = DEF_USE_LEVEE .and. has_levee(i) .and. (.not. is_built_resv(i))
            IF (leveed) THEN
               pub_floodarea (i) = levee_floodarea(i)
               pub_flooddepth(i) = max(levdph(i), max(wdsrf_ucat(i) - floodplain_curve(i)%rivhgt, 0._r8))
            ELSE
               pub_floodarea (i) = floodplain_curve(i)%floodarea (wdsrf_ucat(i))
               pub_flooddepth(i) = max(wdsrf_ucat(i) - floodplain_curve(i)%rivhgt, 0._r8)
            ENDIF
            pub_leveearea(i) = levee_floodarea(i)
         ENDDO
      ENDIF

      IF (.not. mt_lifecycle_route_has_active()) RETURN

      CALL ensure_provider_arrays (numucat)
      DO i = 1, numucat
         prov_storage_end   (i) = 0._r8
         prov_protected_end (i) = 0._r8
         prov_floodarea     (i) = 0._r8
         prov_protected_area(i) = 0._r8
         IF (.not. ucatfilter(i)) CYCLE

         leveed = DEF_USE_LEVEE .and. has_levee(i) .and. (.not. is_built_resv(i))
         ! flooded area as in the river history (a_floodarea)
         IF (leveed) THEN
            prov_floodarea(i) = levee_floodarea(i)
            prov_protected_area(i) = min(max(levee_floodarea(i) - levee_frc_data(i)*topo_area(i), 0._r8), &
               (1._r8 - levee_frc_data(i)) * topo_area(i))
         ELSE
            prov_floodarea(i) = floodplain_curve(i)%floodarea (wdsrf_ucat(i))
         ENDIF

         IF (is_built_resv(i)) THEN
            prov_storage_end(i) = max(volresv(ucat2resv_mt(i)), 0._r8)
         ELSE
            prov_storage_end(i) = max(volwater_ucat(i), 0._r8)
            IF (leveed) THEN
               prov_storage_end(i)   = prov_storage_end(i) + max(levsto(i), 0._r8)
               prov_protected_end(i) = max(levsto(i), 0._r8)
            ENDIF
         ENDIF
      ENDDO

      CALL mt_lifecycle_route_diag_accumulate (dt_all, irivsys, ucatfilter, &
         veloc_riv, wdsrf_ucat, prov_storage_start, prov_storage_end, hflux_fc, prov_floodarea, &
         prov_protected_start, prov_protected_end, prov_protected_area)

   END SUBROUTINE mt_rlc_diag_accumulate

   ! ---------
   SUBROUTINE mt_rlc_debug_discharge (i, dt)

   ! CoLMDEBUG: material leaving the network at outlet i in this substep.

   IMPLICIT NONE
   integer,  intent(in) :: i
   real(r8), intent(in) :: dt
#ifdef CoLMDEBUG
   integer :: imt

      IF (nmts <= 0) RETURN
      CALL ensure_debug_arrays ()
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         mt_mass_dis(imt) = mt_mass_dis(imt) + mt_flux_out(imt,i) * dt
      ENDDO
#endif

   END SUBROUTINE mt_rlc_debug_discharge

   ! ---------
   SUBROUTINE mt_rlc_end_routing ()

   ! End of a routing call, on all processes: CoLMDEBUG material balance,
   ! donor-limiter statistics, reset of the per-call accumulations, and the
   ! provider materials of the call (acctime_rnof is still its length).

   USE MOD_Grid_RiverLakeNetwork,  only: numucat
   USE MOD_Grid_RiverLakeTimeVars, only: acctime_rnof
   IMPLICIT NONE

   integer :: lim_calls, lim_iter_sum, lim_iter_peak, lim_over_soft
#ifdef CoLMDEBUG
   integer :: imt
   real(r8), allocatable :: mt_mass_aft(:), mt_mass_decay(:)
#endif

      CALL mt_riverlake_release_pending_drain ()
      IF (nmts <= 0) RETURN

#ifdef CoLMDEBUG
      CALL ensure_debug_arrays ()
      allocate (mt_mass_aft(nmts), mt_mass_decay(nmts))
      mt_mass_aft = 0._r8
      mt_mass_decay = 0._r8
      IF (p_is_worker .and. numucat > 0) THEN
         DO imt = 1, nmts
            IF (.not. mt_uses_generic_water(imt)) CYCLE
            mt_mass_aft(imt) = pool_total(imt)
            mt_mass_dis(imt) = mt_mass_dis(imt) + sum(mt_dry_drain(imt,:))
            IF (allocated(mt_decay_source)) mt_mass_decay(imt) = sum(mt_decay_source(imt,:))
         ENDDO
      ELSE
         mt_mass_bef = 0._r8
         mt_mass_inp = 0._r8
         mt_mass_dis = 0._r8
      ENDIF
#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, mt_mass_bef,   nmts, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, mt_mass_aft,   nmts, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, mt_mass_inp,   nmts, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, mt_mass_dis,   nmts, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, mt_mass_decay, nmts, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         DO imt = 1, nmts
            IF (.not. mt_uses_generic_water(imt)) CYCLE
            write(*,'(A,I0,A,ES12.4,A)') 'Material(', imt, ') mass before  : ', mt_mass_bef(imt), ' R*m3'
            write(*,'(A,I0,A,ES12.4,A)') 'Material(', imt, ') mass input   : ', mt_mass_inp(imt), ' R*m3'
            write(*,'(A,I0,A,ES12.4,A)') 'Material(', imt, ') mass discharge: ', mt_mass_dis(imt), ' R*m3'
            write(*,'(A,I0,A,ES12.4,A)') 'Material(', imt, ') mass decay   : ', mt_mass_decay(imt), ' R*m3'
            write(*,'(A,I0,A,ES12.4,A)') 'Material(', imt, ') mass after   : ', mt_mass_aft(imt), ' R*m3'
            write(*,'(A,I0,A,ES12.4,A)') 'Material(', imt, ') mass balance : ', &
               mt_mass_aft(imt) - mt_mass_bef(imt) - mt_mass_inp(imt) + mt_mass_dis(imt) &
               - mt_mass_decay(imt), ' R*m3'
         ENDDO
      ENDIF
      mt_mass_bef = 0._r8
      mt_mass_inp = 0._r8
      mt_mass_dis = 0._r8
      deallocate (mt_mass_aft, mt_mass_decay)
#endif

      IF (p_is_worker) THEN
         CALL mt_riverlake_limiter_stats (lim_calls, lim_iter_sum, lim_iter_peak, lim_over_soft, &
            reset = .true.)
#ifdef USEMPI
         CALL mpi_allreduce (MPI_IN_PLACE, lim_calls,     1, MPI_INTEGER, MPI_SUM, p_comm_worker, p_err)
         CALL mpi_allreduce (MPI_IN_PLACE, lim_iter_sum,  1, MPI_INTEGER, MPI_SUM, p_comm_worker, p_err)
         CALL mpi_allreduce (MPI_IN_PLACE, lim_over_soft, 1, MPI_INTEGER, MPI_SUM, p_comm_worker, p_err)
         CALL mpi_allreduce (MPI_IN_PLACE, lim_iter_peak, 1, MPI_INTEGER, MPI_MAX, p_comm_worker, p_err)
#endif
         IF (p_iam_worker == p_root) THEN
            IF (lim_calls > 0 .and. (lim_diag_printed < 5 .or. lim_over_soft > 0)) THEN
               IF (lim_diag_printed < 5) lim_diag_printed = lim_diag_printed + 1
               write(*,'(A,I0,A,F9.3,A,I0,A,I0)') 'River material donor limiter: calls=', lim_calls, &
                  ' mean iter=', real(lim_iter_sum, r8) / real(lim_calls, r8), &
                  ' peak=', lim_iter_peak, ' over_soft=', lim_over_soft
            ENDIF
         ENDIF

         IF (numucat > 0) THEN
            IF (allocated(acc_mt_inp))   acc_mt_inp   = 0._r8
            IF (allocated(acc_rnof_ref))  acc_rnof_ref  = 0._r8
            IF (allocated(mt_dry_drain)) mt_dry_drain = 0._r8
         ENDIF

         ! provider materials (sediment) over the whole routing call
         IF (mt_lifecycle_route_has_active()) CALL mt_lifecycle_route_calc (acctime_rnof)

         CALL publish_flood_to_patches ()
      ENDIF

   END SUBROUTINE mt_rlc_end_routing


   ! ---------
   SUBROUTINE publish_flood_to_patches ()

   ! Flood fraction (and depth) of the last substep mapped to the patches, as
   ! CoLM202X publish_levee_fldfrc_to_patches / publish_fldfrc_to_patches did,
   ! for land providers that read the river flood.  Workers only.

   USE MOD_Grid_RiverLakeNetwork, only: numucat, numinpm, topo_area, push_ucat2inpm, remap_patch2inpm
   USE MOD_LandPatch, only: numpatch
   USE MOD_WorkerPushData
   IMPLICIT NONE

   real(r8), parameter :: missing = -1.e30_r8
   real(r8), allocatable :: frc_uc(:), wat_uc(:), frc_gd(:), wat_gd(:)
   real(r8), allocatable :: frc_patch(:), wat_patch(:), dph_patch(:)
   integer :: i
   logical :: levee_pub, flood_pub

      levee_pub = material_lifecycle_has_levee_flood_publisher()
      flood_pub = material_lifecycle_has_flood_publisher()
      IF (.not. (levee_pub .or. flood_pub)) RETURN
      CALL ensure_publish_arrays (numucat)

      allocate (frc_uc(numucat), wat_uc(numucat), frc_gd(numinpm), wat_gd(numinpm))
      allocate (frc_patch(numpatch), wat_patch(numpatch), dph_patch(numpatch))

      IF (levee_pub) THEN
         DO i = 1, numucat
            frc_uc(i) = 0._r8
            IF (topo_area(i) > 0._r8) frc_uc(i) = min(1._r8, max(0._r8, pub_leveearea(i) / topo_area(i)))
         ENDDO
         CALL worker_push_data (push_ucat2inpm, frc_uc, frc_gd, fillvalue = missing, mode = 'average')
         CALL worker_remap_data_grid2pset (remap_patch2inpm, frc_gd, frc_patch, &
            fillvalue = missing, mode = 'average')
         WHERE (frc_patch == missing) frc_patch = 0._r8
         CALL mt_lifecycle_publish_levee_flood_patch (frc_patch)
      ENDIF

      IF (flood_pub) THEN
         DO i = 1, numucat
            frc_uc(i) = 0._r8
            wat_uc(i) = 0._r8
            IF (topo_area(i) > 0._r8) THEN
               frc_uc(i) = min(1._r8, max(0._r8, pub_floodarea(i) / topo_area(i)))
               wat_uc(i) = frc_uc(i) * max(0._r8, pub_flooddepth(i))
            ENDIF
         ENDDO
         CALL worker_push_data (push_ucat2inpm, frc_uc, frc_gd, fillvalue = missing, mode = 'average')
         CALL worker_remap_data_grid2pset (remap_patch2inpm, frc_gd, frc_patch, &
            fillvalue = missing, mode = 'average')
         WHERE (frc_patch == missing) frc_patch = 0._r8
         CALL worker_push_data (push_ucat2inpm, wat_uc, wat_gd, fillvalue = missing, mode = 'average')
         CALL worker_remap_data_grid2pset (remap_patch2inpm, wat_gd, wat_patch, &
            fillvalue = missing, mode = 'average')
         WHERE (wat_patch == missing) wat_patch = 0._r8
         dph_patch = 0._r8
         WHERE (frc_patch > 0._r8) dph_patch = max(0._r8, wat_patch) / frc_patch
         CALL mt_lifecycle_publish_flood_patch (frc_patch, dph_patch)
      ENDIF

      deallocate (frc_uc, wat_uc, frc_gd, wat_gd, frc_patch, wat_patch, dph_patch)

      ! a new routing call starts from no flood, as CoLM202X allocated these per call
      pub_floodarea  = 0._r8
      pub_flooddepth = 0._r8
      pub_leveearea  = 0._r8

   END SUBROUTINE publish_flood_to_patches


   ! ---------
   SUBROUTINE ensure_publish_arrays (n)

   integer, intent(in) :: n

      IF (allocated(pub_floodarea)) THEN
         IF (size(pub_floodarea) == n) RETURN
         deallocate (pub_floodarea, pub_flooddepth, pub_leveearea)
      ENDIF
      allocate (pub_floodarea(n), pub_flooddepth(n), pub_leveearea(n))
      pub_floodarea  = 0._r8
      pub_flooddepth = 0._r8
      pub_leveearea  = 0._r8

   END SUBROUTINE ensure_publish_arrays

   ! ---------
   SUBROUTINE mt_rlc_flood_init ()

   ! Flood-feedback material arrays on the current land patches: at start-up
   ! and after LULCC (called from init_flood_feedback, workers).

   USE MOD_LandPatch,             only: numpatch
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

   integer :: imt
   logical :: has_flood_material

      IF (.not. p_is_worker) RETURN

      IF (allocated(flood_mt_credit_patch)) deallocate (flood_mt_credit_patch, &
         flood_mt_evap_patch, flood_mt_land_patch, &
         flood_visible_mt_uc, flood_protected_mt_uc)

      has_flood_material = .false.
      DO imt = 1, nmts
         IF (mt_uses_generic_water(imt)) has_flood_material = .true.
      ENDDO
      IF (.not. has_flood_material) RETURN

      allocate (flood_mt_credit_patch(nmts,numpatch), flood_mt_evap_patch(nmts,numpatch), &
         flood_mt_land_patch(nmts,numpatch), &
         flood_visible_mt_uc(nmts,numucat), flood_protected_mt_uc(nmts,numucat))
      flood_mt_credit_patch = 0._r8
      flood_mt_evap_patch   = 0._r8
      flood_mt_land_patch   = 0._r8
      flood_visible_mt_uc   = 0._r8
      flood_protected_mt_uc = 0._r8

   END SUBROUTINE mt_rlc_flood_init

   ! ---------
   SUBROUTINE mt_rlc_flood_withhold ()

   ! No flood is published to the land: no material either.

      IF (.not. p_is_worker) RETURN
      IF (allocated(flood_mt_credit_patch)) flood_mt_credit_patch = 0._r8

   END SUBROUTINE mt_rlc_flood_withhold

   ! ---------
   SUBROUTINE mt_rlc_flood_publish_reset ()

   ! Start of publish_flood_feedback: a new flood state replaces the old one.

      IF (.not. p_is_worker) RETURN
      IF (.not. allocated(flood_visible_mt_uc)) RETURN
      flood_visible_mt_uc   = 0._r8
      flood_protected_mt_uc = 0._r8
      flood_mt_credit_patch = 0._r8
      flood_mt_evap_patch   = 0._r8
      flood_mt_land_patch   = 0._r8

   END SUBROUTINE mt_rlc_flood_publish_reset

   ! ---------
   SUBROUTINE mt_rlc_flood_publish_cell (i, visible, protected, visible_flood)

   ! Material of the flood water published from unit catchment i:
   ! the visible part above the river channel at the visible dissolved ratio,
   ! and all protected (levee) water.

   IMPLICIT NONE
   integer,  intent(in) :: i
   real(r8), intent(in) :: visible        ! visible water [m^3]
   real(r8), intent(in) :: protected      ! protected water [m^3]
   real(r8), intent(in) :: visible_flood  ! published visible flood water [m^3]

      IF (.not. allocated(flood_visible_mt_uc)) RETURN
      IF (visible > 0._r8) flood_visible_mt_uc(:,i) = &
         max(mt_mass(:,i), 0._r8) * (visible_flood / visible)
      IF (protected > 0._r8) flood_protected_mt_uc(:,i) = max(mt_levsto(:,i), 0._r8)

   END SUBROUTINE mt_rlc_flood_publish_cell

   ! ---------
   SUBROUTINE mt_rlc_flood_publish_push (flood_grid_area)

   ! Map the published flood material to the patches, as the flood water is
   ! (unit catchment -> input grid -> patch), into flood_mt_credit_patch.
   ! Called on all workers (the pushes are collective).

   USE MOD_LandPatch,             only: numpatch
   USE MOD_Grid_RiverLakeNetwork, only: numinpm, numucat, remap_patch2inpm, &
      push_inpm2ucat, push_ucat2inpm
   USE MOD_WorkerPushData
   IMPLICIT NONE

   real(r8), intent(in) :: flood_grid_area(:)   ! land area of each input grid cell

   real(r8), parameter :: MISSING = -1.e30_r8
   integer :: imt, iside
   real(r8), allocatable :: density(:), grid_material(:), patch_material(:)

      IF (.not. p_is_worker) RETURN
      IF (.not. allocated(flood_visible_mt_uc)) RETURN

      allocate (density(numucat), grid_material(numinpm), patch_material(numpatch))
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         DO iside = 1, 2
            ! sum_area is allocated only on workers with requests (numucat > 0)
            density = 0._r8
            IF (numucat > 0) THEN
               IF (iside == 1) THEN
                  WHERE (push_inpm2ucat%sum_area > 0._r8) &
                     density = flood_visible_mt_uc(imt,:) / max(push_inpm2ucat%sum_area, tiny(1._r8))
               ELSE
                  WHERE (push_inpm2ucat%sum_area > 0._r8) &
                     density = flood_protected_mt_uc(imt,:) / max(push_inpm2ucat%sum_area, tiny(1._r8))
               ENDIF
            ENDIF
            CALL worker_push_data (push_ucat2inpm, density, grid_material, fillvalue=0._r8, mode='sum')
            WHERE (flood_grid_area > 0._r8)
               grid_material = grid_material / max(flood_grid_area, tiny(1._r8))
            ELSEWHERE
               grid_material = 0._r8
            END WHERE
            IF (numpatch > 0) THEN
               CALL worker_remap_data_grid2pset (remap_patch2inpm, grid_material, &
                  patch_material, fillvalue=MISSING, mode='average')
               WHERE (patch_material == MISSING) patch_material = 0._r8
               IF (iside == 1) THEN
                  flood_mt_credit_patch(imt,:) = patch_material
               ELSE
                  ! [m] -> [mm], as flood_credit_patch*1000 is the water credit [mm]
                  flood_mt_credit_patch(imt,:) = 1000._r8 * &
                     (flood_mt_credit_patch(imt,:) + patch_material)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      deallocate (density, grid_material, patch_material)

   END SUBROUTINE mt_rlc_flood_publish_push

   ! ---------
   SUBROUTINE mt_rlc_flood_debit (flood_credit_patch, flood_evap_acc, flood_infil_acc, &
         flood_grid_area, flood_visible_uc, flood_protected_uc)

   ! Material side of debit_flood_feedback (workers, collective): take from
   ! the river pools what the land took with flood evaporation (vapour
   ! exchange, flood_mt_evap_patch) and infiltration (flood_mt_land_patch),
   ! and stop if the land-side and river-side totals do not close.

   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   USE MOD_LandPatch,             only: numpatch
   USE MOD_Grid_RiverLakeNetwork, only: numinpm, numucat, remap_patch2inpm, push_inpm2ucat
   USE MOD_WorkerPushData
   IMPLICIT NONE

   real(r8), intent(in) :: flood_credit_patch(:)   ! published water credit [m]
   real(r8), intent(in) :: flood_evap_acc(:)       ! land flood evaporation [mm]
   real(r8), intent(in) :: flood_infil_acc(:)      ! land flood infiltration [mm]
   real(r8), intent(in) :: flood_grid_area(:)
   real(r8), intent(in) :: flood_visible_uc(:)     ! published visible flood water [m^3]
   real(r8), intent(in) :: flood_protected_uc(:)   ! published protected flood water [m^3]

   integer  :: i, j, imt
   real(r8) :: patch_area, material_cell_debit, ledger_tol
   real(r8) :: water_credit, material_credit, vapor_loss, positive_loss, vapor_gain
   real(r8) :: infiltrated_fraction, evaporated_fraction, material_after
   real(r8), allocatable :: ratio_patch(:), ratio_grid(:), ratio_uc(:)
   real(r8), allocatable :: gain_patch(:), coefficient_uc(:), gain_uc(:)
   real(r8), allocatable :: material_exchange_ledger(:,:)   ! (land, vapour, river donor)

      IF (.not. p_is_worker) RETURN
      IF (.not. allocated(flood_visible_mt_uc)) RETURN

      allocate (ratio_patch(numpatch), ratio_grid(numinpm), ratio_uc(numucat))
      allocate (gain_patch(numpatch), coefficient_uc(numucat), gain_uc(numucat))
      allocate (material_exchange_ledger(3,nmts))
      material_exchange_ledger = 0._r8

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         ratio_patch = 0._r8
         gain_patch = 0._r8
         DO i = 1, numpatch
            water_credit = flood_credit_patch(i)*1000._r8
            IF (water_credit <= 0._r8) CYCLE
            material_credit = flood_mt_credit_patch(imt,i)
            vapor_loss = flood_mt_evap_patch(imt,i)
            IF (.not. ieee_is_finite(vapor_loss)) &
               CALL CoLM_stop('grid flood feedback: nonfinite isotope vapor exchange')
            positive_loss = max(vapor_loss,0._r8)
            vapor_gain = max(-vapor_loss,0._r8)
            ledger_tol = max(1.e-12_r8, 1.e-10_r8*abs(material_credit))
            IF (positive_loss > material_credit+ledger_tol) &
               CALL CoLM_stop('grid flood feedback: isotope evaporation overdrew material')
            evaporated_fraction = 0._r8
            IF (material_credit > 0._r8) evaporated_fraction = min(1._r8,positive_loss/material_credit)
            infiltrated_fraction = 0._r8
            IF (flood_infil_acc(i) > 0._r8) THEN
               IF (water_credit <= flood_evap_acc(i)) THEN
                  IF (flood_infil_acc(i) > 1.e-12_r8*max(water_credit, 1._r8)) &
                     CALL CoLM_stop('grid flood feedback: infiltration after full evaporation')
               ELSE
                  infiltrated_fraction = flood_infil_acc(i)/(water_credit-flood_evap_acc(i))
               ENDIF
            ENDIF
            IF (infiltrated_fraction > 1._r8+1.e-10_r8) &
               CALL CoLM_stop('grid flood feedback: material infiltration overdraft')
            infiltrated_fraction = min(1._r8,max(0._r8,infiltrated_fraction))
            ratio_patch(i) = evaporated_fraction + (1._r8-evaporated_fraction)*infiltrated_fraction
            IF (mt_has_dissolved_limit(imt) .and. material_credit > 0._r8) &
               ratio_patch(i) = min(1._r8,max(0._r8,flood_mt_land_patch(imt,i)/material_credit))
            gain_patch(i) = vapor_gain*(1._r8-infiltrated_fraction)/water_credit
            patch_area = sum(remap_patch2inpm%areapart(i)%val)*1.e-3_r8
            IF (mt_has_dissolved_limit(imt)) THEN
               IF (flood_mt_land_patch(imt,i) < -ledger_tol .or. &
                   flood_mt_land_patch(imt,i) > &
                   (material_credit-vapor_loss)*infiltrated_fraction + ledger_tol) &
                  CALL CoLM_stop('grid flood feedback: finite solute input exceeds donor exchange')
            ELSE
               IF (abs(flood_mt_land_patch(imt,i) - &
                   (material_credit-vapor_loss)*infiltrated_fraction) > &
                   max(1.e-12_r8,1.e-10_r8*max(abs(material_credit),abs(flood_mt_land_patch(imt,i))))) &
                  CALL CoLM_stop('grid flood feedback: land material input differs from donor exchange')
            ENDIF
            material_exchange_ledger(1,imt) = material_exchange_ledger(1,imt) + &
               flood_mt_land_patch(imt,i)*patch_area
            material_exchange_ledger(2,imt) = material_exchange_ledger(2,imt) + vapor_loss*patch_area
         ENDDO

         ! fraction of the published material taken, per unit catchment
         ratio_grid = 0._r8
         IF (numpatch > 0) CALL worker_remap_data_pset2grid (remap_patch2inpm, ratio_patch, ratio_grid, &
            fillvalue=0._r8, filter=flood_credit_patch>0._r8)
         WHERE (flood_grid_area > 0._r8)
            ratio_grid = ratio_grid/max(flood_grid_area,tiny(1._r8))
         ELSEWHERE
            ratio_grid = 0._r8
         END WHERE
         CALL worker_push_data (push_inpm2ucat, ratio_grid, ratio_uc, fillvalue=0._r8, mode='sum')
         coefficient_uc = 0._r8
         IF (numucat > 0) THEN
            WHERE (push_inpm2ucat%sum_area > 0._r8) &
               coefficient_uc = min(1._r8,max(0._r8,ratio_uc/push_inpm2ucat%sum_area))
         ENDIF

         ! vapour gained by the flood water still in the river (isotopes only)
         gain_uc = 0._r8
         IF (.not. material_is_nonvolatile_solute(imt)) THEN
            ratio_grid = 0._r8
            IF (numpatch > 0) CALL worker_remap_data_pset2grid (remap_patch2inpm, gain_patch, ratio_grid, &
               fillvalue=0._r8, filter=flood_credit_patch>0._r8)
            WHERE (flood_grid_area > 0._r8)
               ratio_grid = ratio_grid/max(flood_grid_area,tiny(1._r8))
            ELSEWHERE
               ratio_grid = 0._r8
            END WHERE
            CALL worker_push_data (push_inpm2ucat, ratio_grid, ratio_uc, fillvalue=0._r8, mode='sum')
            IF (numucat > 0) THEN
               WHERE (push_inpm2ucat%sum_area > 0._r8) gain_uc = ratio_uc/push_inpm2ucat%sum_area
            ENDIF
         ENDIF

         DO j = 1, numucat
            material_after = mt_mass(imt,j) - coefficient_uc(j)*flood_visible_mt_uc(imt,j) &
               + gain_uc(j)*flood_visible_uc(j)
            IF (material_after < -1.e-10_r8*max(1._r8,mt_mass(imt,j))) &
               CALL CoLM_stop('grid flood feedback: negative visible material')
            material_cell_debit = coefficient_uc(j)*flood_visible_mt_uc(imt,j) &
               - gain_uc(j)*flood_visible_uc(j) + min(material_after,0._r8)
            mt_mass(imt,j) = max(0._r8,material_after)
            material_after = mt_levsto(imt,j) - coefficient_uc(j)*flood_protected_mt_uc(imt,j) &
               + gain_uc(j)*flood_protected_uc(j)
            IF (material_after < -1.e-10_r8*max(1._r8,mt_levsto(imt,j))) &
               CALL CoLM_stop('grid flood feedback: negative protected material')
            material_cell_debit = material_cell_debit + coefficient_uc(j)*flood_protected_mt_uc(imt,j) &
               - gain_uc(j)*flood_protected_uc(j) + min(material_after,0._r8)
            mt_levsto(imt,j) = max(0._r8,material_after)
            material_exchange_ledger(3,imt) = material_exchange_ledger(3,imt) + material_cell_debit
         ENDDO
      ENDDO

#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, material_exchange_ledger, 3*nmts, MPI_REAL8, &
         MPI_SUM, p_comm_worker, p_err)
#endif
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         ledger_tol = max(1.e-8_r8, 1.e-10_r8*maxval(abs(material_exchange_ledger(:,imt))))
         IF (abs(material_exchange_ledger(1,imt)+material_exchange_ledger(2,imt) - &
                material_exchange_ledger(3,imt)) > ledger_tol) &
            CALL CoLM_stop('grid flood feedback: land/river material exchange ledger mismatch')
         IF (p_iam_worker == p_root .and. flood_mt_ledger_reports < 5 .and. &
             material_exchange_ledger(1,imt) > 0._r8) THEN
            write(*,'(A,I0,4(1X,ES16.8))') 'Grid flood material land/vapor/donor/residual:', imt, &
               material_exchange_ledger(1,imt), material_exchange_ledger(2,imt), &
               material_exchange_ledger(3,imt), &
               material_exchange_ledger(1,imt)+material_exchange_ledger(2,imt)-material_exchange_ledger(3,imt)
            flood_mt_ledger_reports = flood_mt_ledger_reports + 1
         ENDIF
      ENDDO

      deallocate (ratio_patch, ratio_grid, ratio_uc, gain_patch, coefficient_uc, gain_uc, &
         material_exchange_ledger)

   END SUBROUTINE mt_rlc_flood_debit

   ! ---------
   SUBROUTINE mt_rlc_final ()

      CALL mt_lifecycle_route_final ()
      CALL mt_riverlake_final ()

      IF (allocated(flood_mt_credit_patch)) deallocate (flood_mt_credit_patch, &
         flood_mt_evap_patch, flood_mt_land_patch, &
         flood_visible_mt_uc, flood_protected_mt_uc)
      IF (allocated(volresv_mt))   deallocate (volresv_mt)
      IF (allocated(ucat2resv_mt)) deallocate (ucat2resv_mt)
      IF (allocated(prov_storage_start)) deallocate (prov_storage_start, prov_protected_start, &
         prov_storage_end, prov_protected_end, prov_floodarea, prov_protected_area)
      IF (allocated(pub_floodarea)) deallocate (pub_floodarea, pub_flooddepth, pub_leveearea)
#ifdef CoLMDEBUG
      IF (allocated(mt_mass_bef)) deallocate (mt_mass_bef, mt_mass_inp, mt_mass_dis)
#endif

   END SUBROUTINE mt_rlc_final

   ! ---------
   SUBROUTINE ensure_provider_arrays (n)

   integer, intent(in) :: n

      IF (allocated(prov_storage_start)) THEN
         IF (size(prov_storage_start) == n) RETURN
         deallocate (prov_storage_start, prov_protected_start, prov_storage_end, &
            prov_protected_end, prov_floodarea, prov_protected_area)
      ENDIF
      allocate (prov_storage_start(n), prov_protected_start(n), prov_storage_end(n), &
         prov_protected_end(n), prov_floodarea(n), prov_protected_area(n))
      prov_storage_start   = 0._r8
      prov_protected_start = 0._r8
      prov_storage_end     = 0._r8
      prov_protected_end   = 0._r8
      prov_floodarea       = 0._r8
      prov_protected_area  = 0._r8

   END SUBROUTINE ensure_provider_arrays

#ifdef CoLMDEBUG
   ! ---------
   SUBROUTINE ensure_debug_arrays ()

      IF (allocated(mt_mass_bef)) RETURN
      allocate (mt_mass_bef(nmts), mt_mass_inp(nmts), mt_mass_dis(nmts))
      mt_mass_bef = 0._r8
      mt_mass_inp = 0._r8
      mt_mass_dis = 0._r8

   END SUBROUTINE ensure_debug_arrays

   ! ---------
   real(r8) FUNCTION pool_total (imt)

   ! all river pools of material imt on this worker [R*m^3]

   integer, intent(in) :: imt

      pool_total = sum(mt_mass(imt,:)) + sum(mt_inp_buf(imt,:))
      IF (allocated(mt_levsto)) pool_total = pool_total + sum(mt_levsto(imt,:))
      IF (allocated(mt_solid))  pool_total = pool_total + sum(mt_solid(imt,:)) &
         + sum(mt_levsto_solid(imt,:))

   END FUNCTION pool_total
#endif

END MODULE MOD_MT_RiverLakeCoupling
#endif
