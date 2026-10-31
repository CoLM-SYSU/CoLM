#include <define.h>

#if (defined MaterialTransport) && (defined BGC)
MODULE MOD_MT_Methane_BgcLink
!=======================================================================
! BGC coupling bridge for the methane reactive material.
!
! This is the single place where Methane reaches into BGC internals.
! CoLM202X stores the required carbon fluxes in pool/PFT-level arrays, so
! this bridge reconstructs patch-level inputs for Methane without moving
! material-specific fields back into MOD_BGC_* modules.
!=======================================================================

   USE MOD_Precision
   USE MOD_SPMD_Task, only: CoLM_stop
   USE MOD_TimeManager, only: isleapyear
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   USE MOD_Vars_Global, only: nl_soil, dz_soi, spval
   USE MOD_MT_Methane_Const, only: DEF_METHANE, catomw, gc_per_kg_om, &
      METHANE_COMP_SOIL, METHANE_COMP_RICE, N_METHANE_COMP
   USE MOD_MT_Methane_pH, only: get_ph_for_patch
   USE MOD_LandPFT, only: patch_pft_s, patch_pft_e
   USE MOD_Vars_PFTimeInvariants,  only: pftfrac, pftclass
   USE MOD_Vars_PFTimeVariables,   only: lai_p, irrig_method_p
   USE MOD_BGC_Vars_TimeInvariants, only: organic_max, &
      i_met_lit, i_cel_lit, i_lig_lit, i_cwd, i_soil1, i_soil2, i_soil3, &
      donor_pool, is_litter, is_soil, is_cwd
   USE MOD_BGC_Vars_1DFluxes,       only: decomp_hr_vr, pot_f_nit_vr
   USE MOD_BGC_Vars_TimeVariables,  only: decomp_cpools_vr, o_scalar
   ! Wetland patches: pools and respiration come from the methane shadow
   ! pools (MOD_MT_Methane_WetlandBGC); the host BGC arrays are not changed.
   USE MOD_MT_Methane_WetlandBGC, only: methane_wetland_bgc_view_active, methane_wetland_bgc_apply, &
      wl_cpools, wv_hr_vr, wv_o_scalar, wv_pot_f_nit_vr
   USE MOD_BGC_Vars_PFTimeVariables, only: annsum_npp_p, cinput_rootfr_p
   USE MOD_MT_Methane_VegOverride, only: wetland_aere_poros, wetland_aere_radius, &
                                              wetland_aere_tillerC, wetland_aere_scale, &
                                              wetland_aere_active
   USE MOD_BGC_Vars_1DPFTFluxes,    only: froot_mr_p, &
      cpool_to_leafc_p, cpool_to_leafc_storage_p, &
      cpool_to_livestemc_p, cpool_to_livestemc_storage_p, &
      cpool_to_deadstemc_p, cpool_to_deadstemc_storage_p, &
      cpool_to_frootc_p, cpool_to_frootc_storage_p, &
      cpool_to_livecrootc_p, cpool_to_livecrootc_storage_p, &
      cpool_to_deadcrootc_p, cpool_to_deadcrootc_storage_p, &
      cpool_froot_gr_p, cpool_froot_storage_gr_p, &
      cpool_livecroot_gr_p, cpool_livecroot_storage_gr_p, &
      cpool_deadcroot_gr_p, cpool_deadcroot_storage_gr_p, &
      transfer_froot_gr_p, transfer_livecroot_gr_p, transfer_deadcroot_gr_p

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: mt_ch4_bgc_patch_inputs
   PUBLIC :: mt_ch4_bgc_component_veg_inputs
   PUBLIC :: mt_ch4_bgc_finalize_step
   PUBLIC :: get_wetland_veg_proxy
   PUBLIC :: get_rice_veg_proxy
   PUBLIC :: paddy_rice_fraction
   PUBLIC :: has_paddy_rice_tile
   PUBLIC :: is_paddy_rice_live
   PUBLIC :: rice_days_past_planting
   PUBLIC :: rice_days_since_harvest
   PUBLIC :: methane_patch_active_mask
   PUBLIC :: organic_max
   PUBLIC :: get_biome_f_methane
   PUBLIC :: get_biome_redoxlag

   ! Driver gate threshold for "this patch hosts a paddy rice tile".
   ! Exposed so CoLMDRIVER / AccFlux / Hist filter all use the SAME
   ! cutoff via paddy_rice_fraction / has_paddy_rice_tile and never
   ! drift out of sync.
   real(r8), public, parameter :: PADDY_RICE_FRAC_MIN = 0.01_r8
   real(r8), parameter :: O_SCALAR_CONTRACT_TOL = 1.e-8_r8

CONTAINS

   SUBROUTINE mt_ch4_bgc_finalize_step(ipatch, patchtype, deltim, net_methane, co2_hr_out)
   ! CO2-C part of the heterotrophic respiration, total respiration plus the
   ! net methane flux.  CoLM202X wrote it into the host decomp_hr / er (and
   ! set ar = 0 on wetlands) and updated the wetland pools in the host BGC
   ! arrays; here the host keeps its upstream values (decomp_hr is the total
   ! decomposition respiration) and the CO2 part is methane output
   ! (co2_hr_out, history f_co2_hr_methane).  Wetland pools are the methane
   ! shadow pools, updated here.

      integer,  intent(in)  :: ipatch, patchtype
      real(r8), intent(in)  :: deltim, net_methane
      real(r8), intent(out) :: co2_hr_out
      real(r8) :: co2_hr, total_hr

      co2_hr_out = 0._r8
      IF (patchtype /= 0 .and. patchtype /= 2) RETURN

      IF (patchtype == 2) THEN
         IF (.not. methane_wetland_bgc_view_active(ipatch)) &
            CALL CoLM_stop(' ***** ERROR: wetland CH4 step without its decomposition cascade')
         total_hr = sum(sum(wv_hr_vr(1:nl_soil,:), dim=2) * dz_soi(1:nl_soil))
         CALL methane_wetland_bgc_apply (ipatch, deltim)
      ELSE
         total_hr = sum(sum(decomp_hr_vr(1:nl_soil,:,ipatch), dim=2) * dz_soi(1:nl_soil))
      ENDIF
      IF (.not. ieee_is_finite(total_hr) .or. total_hr < -1.e-12_r8 .or. &
          .not. ieee_is_finite(net_methane)) THEN
         CALL CoLM_stop(' ***** ERROR: CH4/BGC carbon partition received invalid respiration')
      ENDIF
      co2_hr = total_hr + catomw * net_methane
      IF (.not. ieee_is_finite(co2_hr) .or. co2_hr < -1.e-12_r8) THEN
         CALL CoLM_stop(' ***** ERROR: CH4/BGC carbon partition produced negative CO2 respiration')
      ENDIF
      co2_hr_out = max(co2_hr, 0._r8)

   END SUBROUTINE mt_ch4_bgc_finalize_step

   SUBROUTINE mt_ch4_bgc_patch_inputs (ipatch, rootfr, crootfr, pH, cellorg, &
      somhr, lithr, hr_vr, rr, agnpp, bgnpp, annsum_npp, fphr, &
      o_scalar_out, pot_f_nit_vr_out, bgc_inputs_ready)
      integer,  intent(in)  :: ipatch
      real(r8), intent(in)  :: rootfr(1:nl_soil)
      real(r8), intent(out) :: crootfr(1:nl_soil)
      real(r8), intent(out) :: pH
      real(r8), intent(out) :: cellorg(1:nl_soil)
      real(r8), intent(out) :: somhr
      real(r8), intent(out) :: lithr
      real(r8), intent(out) :: hr_vr(1:nl_soil)
      real(r8), intent(out) :: rr
      real(r8), intent(out) :: agnpp
      real(r8), intent(out) :: bgnpp
      real(r8), intent(out) :: annsum_npp
      real(r8), intent(out) :: fphr(1:nl_soil)
      real(r8), intent(out) :: o_scalar_out(1:nl_soil)
      real(r8), intent(out) :: pot_f_nit_vr_out(1:nl_soil)
      logical,  intent(out), optional :: bgc_inputs_ready

      integer :: ps, pe, j, ipool, m, k, donor
      integer :: decomp_pool_ids(7)
      real(r8) :: frac, profile_sum, unclassified_hr
      logical :: donor_soil, donor_litter_like
      logical :: have_decomp_state, have_pft_state
      logical :: wetland_view
      real(r8), allocatable :: cpools(:,:), hrv(:,:), oscal(:), pfnit(:)

      ps = -1
      pe = -1
      IF (allocated(patch_pft_s) .and. allocated(patch_pft_e)) THEN
         IF (ipatch >= 1 .and. ipatch <= size(patch_pft_s) .and. &
             ipatch <= size(patch_pft_e)) THEN
            ps = patch_pft_s(ipatch)
            pe = patch_pft_e(ipatch)
         ENDIF
      ENDIF

	      crootfr(:) = 0._r8
		      pH = get_ph_for_patch(ipatch)
      cellorg(:) = 0._r8
      somhr = 0._r8
      lithr = 0._r8
      hr_vr(:) = 0._r8
      rr = 0._r8
      agnpp = 0._r8
      bgnpp = 0._r8
      annsum_npp = 0._r8
      fphr(:) = 1.0_r8
      o_scalar_out(:) = 1.0_r8
      pot_f_nit_vr_out(:) = 0._r8
      unclassified_hr = 0._r8

      decomp_pool_ids = (/ i_met_lit, i_cel_lit, i_lig_lit, i_cwd, &
         i_soil1, i_soil2, i_soil3 /)

      ! decomp_cpools_vr is [gC/m3].  Methane physics expects organic
      ! matter density [kg OM/m3].  Keep this conversion consistent with
      ! the wetland/lake surface-data fallback convention:
      !   580 gC per kg OM.
      ! Pools, respiration, o_scalar and pot_f_nit_vr of this patch: the
      ! methane shadow pools for a wetland patch, the host BGC otherwise.
      wetland_view = methane_wetland_bgc_view_active(ipatch)
      have_decomp_state = .false.
      IF (wetland_view) THEN
         cpools = wl_cpools(:,:,ipatch)
         hrv    = wv_hr_vr
         oscal  = wv_o_scalar
         pfnit  = wv_pot_f_nit_vr
         have_decomp_state = .true.
      ELSEIF (allocated(decomp_cpools_vr) .and. allocated(decomp_hr_vr)) THEN
         have_decomp_state = ipatch >= 1 .and. &
                             ipatch <= size(decomp_cpools_vr,3) .and. &
                             ipatch <= size(decomp_hr_vr,3)
         IF (have_decomp_state) THEN
            cpools = decomp_cpools_vr(:,:,ipatch)
            hrv    = decomp_hr_vr(:,:,ipatch)
         ENDIF
      ENDIF
      IF (have_decomp_state) THEN
         have_decomp_state = &
            any(valid_bgc_value(cpools) .and. cpools > 0._r8) .and. &
            any(valid_bgc_value(hrv) .and. hrv > 0._r8)
      ENDIF
      IF (.not. wetland_view) THEN
         IF (allocated(o_scalar) .and. ipatch >= 1 .and. ipatch <= size(o_scalar,2)) &
            oscal = o_scalar(1:nl_soil, ipatch)
         IF (allocated(pot_f_nit_vr) .and. ipatch >= 1 .and. ipatch <= size(pot_f_nit_vr,2)) &
            pfnit = pot_f_nit_vr(1:nl_soil, ipatch)
      ENDIF
	      have_pft_state = pft_arrays_ready(ps, pe)
	      IF (present(bgc_inputs_ready)) bgc_inputs_ready = have_decomp_state .and. have_pft_state

	      IF (have_pft_state) THEN
	         DO m = ps, pe
	            frac = safe_nonnegative(pftfrac(m))
	            DO j = 1, nl_soil
	               IF (valid_bgc_value(cinput_rootfr_p(j,m))) THEN
	                  crootfr(j) = crootfr(j) + frac * max(0._r8, cinput_rootfr_p(j,m)) * dz_soi(j)
	               ENDIF
	            ENDDO
	         ENDDO
	      ENDIF
	      IF (sum(crootfr) <= 0._r8) THEN
	         DO j = 1, nl_soil
	            IF (valid_bgc_value(rootfr(j))) crootfr(j) = max(0._r8, rootfr(j))
	         ENDDO
	      ENDIF
	      profile_sum = sum(crootfr)
	      IF (profile_sum > 0._r8) THEN
	         crootfr(:) = crootfr(:) / profile_sum
	      ELSE
	         crootfr(1) = 1._r8
	      ENDIF

      IF (allocated(cpools)) THEN
         DO j = 1, nl_soil
            DO m = 1, size(decomp_pool_ids)
               ipool = decomp_pool_ids(m)
               IF (ipool < 1 .or. ipool > size(cpools,2)) CYCLE
               IF (valid_bgc_value(cpools(j, ipool))) THEN
                  cellorg(j) = cellorg(j) + max(0._r8, cpools(j, ipool))
               ENDIF
            END DO
         END DO
         cellorg(:) = cellorg(:) / gc_per_kg_om
      ENDIF

      IF (allocated(hrv)) THEN
         DO j = 1, nl_soil
            DO k = 1, size(hrv,2)
               IF (valid_bgc_value(hrv(j, k))) THEN
                  hr_vr(j) = hr_vr(j) + max(0._r8, hrv(j, k))

                  ! decomp_hr_vr's second dimension is decomposition
                  ! transition, not pool.  Classify the flux by the donor
                  ! pool of that transition before building CH4's litter/SOM
                  ! column respiration inputs.
                  donor_soil = .false.
                  donor_litter_like = .false.
                  IF (allocated(donor_pool) .and. k >= lbound(donor_pool,1) .and. &
                      k <= ubound(donor_pool,1)) THEN
                     donor = donor_pool(k)

                     IF (allocated(is_soil)) THEN
                        IF (donor >= lbound(is_soil,1) .and. donor <= ubound(is_soil,1)) &
                           donor_soil = is_soil(donor)
                     ENDIF
                     IF (allocated(is_litter)) THEN
                        IF (donor >= lbound(is_litter,1) .and. donor <= ubound(is_litter,1)) &
                           donor_litter_like = donor_litter_like .or. is_litter(donor)
                     ENDIF
                     IF (allocated(is_cwd)) THEN
                        IF (donor >= lbound(is_cwd,1) .and. donor <= ubound(is_cwd,1)) &
                           donor_litter_like = donor_litter_like .or. is_cwd(donor)
                     ENDIF
                  ENDIF

                  IF (donor_soil) THEN
                     somhr = somhr + max(0._r8, hrv(j, k)) * dz_soi(j)
                  ELSEIF (donor_litter_like) THEN
                     lithr = lithr + max(0._r8, hrv(j, k)) * dz_soi(j)
                  ELSE
                     unclassified_hr = unclassified_hr + max(0._r8, hrv(j, k)) * dz_soi(j)
                  ENDIF
               ENDIF
            END DO
         END DO
      ENDIF
      CALL methane_check_hr_classification_contract (unclassified_hr)

      IF (allocated(oscal)) THEN
         DO j = 1, nl_soil
            IF (valid_bgc_value(oscal(j))) THEN
               o_scalar_out(j) = min(1._r8, max(0._r8, oscal(j)))
            ENDIF
         END DO
      ENDIF
      CALL methane_check_anoxia_contract (o_scalar_out)

      IF (allocated(pfnit)) THEN
         DO j = 1, nl_soil
            IF (valid_bgc_value(pfnit(j))) THEN
               pot_f_nit_vr_out(j) = max(0._r8, pfnit(j))
            ENDIF
         END DO
      ENDIF

      IF (have_pft_state) THEN
         DO m = ps, pe
            frac = safe_nonnegative(pftfrac(m))

            agnpp = agnpp + frac * ( &
               safe_nonnegative(cpool_to_leafc_p(m)) &
             + safe_nonnegative(cpool_to_leafc_storage_p(m)) &
             + safe_nonnegative(cpool_to_livestemc_p(m)) &
             + safe_nonnegative(cpool_to_livestemc_storage_p(m)) &
             + safe_nonnegative(cpool_to_deadstemc_p(m)) &
             + safe_nonnegative(cpool_to_deadstemc_storage_p(m)) )

            bgnpp = bgnpp + frac * ( &
               safe_nonnegative(cpool_to_frootc_p(m)) &
             + safe_nonnegative(cpool_to_frootc_storage_p(m)) &
             + safe_nonnegative(cpool_to_livecrootc_p(m)) &
             + safe_nonnegative(cpool_to_livecrootc_storage_p(m)) &
             + safe_nonnegative(cpool_to_deadcrootc_p(m)) &
             + safe_nonnegative(cpool_to_deadcrootc_storage_p(m)) )

            rr = rr + frac * ( &
               safe_nonnegative(froot_mr_p(m)) &
             + safe_nonnegative(cpool_froot_gr_p(m)) &
             + safe_nonnegative(cpool_froot_storage_gr_p(m)) &
             + safe_nonnegative(cpool_livecroot_gr_p(m)) &
             + safe_nonnegative(cpool_livecroot_storage_gr_p(m)) &
             + safe_nonnegative(cpool_deadcroot_gr_p(m)) &
             + safe_nonnegative(cpool_deadcroot_storage_gr_p(m)) &
             + safe_nonnegative(transfer_froot_gr_p(m)) &
             + safe_nonnegative(transfer_livecroot_gr_p(m)) &
             + safe_nonnegative(transfer_deadcroot_gr_p(m)) )

            annsum_npp = annsum_npp + frac * safe_nonnegative(annsum_npp_p(m))
         END DO
      ENDIF

      ! Source CH4 treats fphr as a per-layer scaling factor. CoLM202X does
      ! not expose an equivalent field, so use full HR contribution while the
      ! rest of BGC coupling remains field-for-field reconstructed above.
   END SUBROUTINE mt_ch4_bgc_patch_inputs


   SUBROUTINE mt_ch4_bgc_component_veg_inputs(ipatch, rootfr_fallback, &
      component_fraction, lai_component, crootfr_component, rr_component, &
      agnpp_component, bgnpp_component, annsum_npp_component, ready)
      ! Reconstruct intensive PFT-driven vegetation forcing for the independent
      ! non-paddy and paddy methane columns.  Patch-soil quantities stay in the
      ! patch bridge above; only fields that actually retain a PFT dimension are
      ! split here.  One PFT scan computes both columns and introduces no MPI.
      integer, intent(in) :: ipatch
      real(r8), intent(in) :: rootfr_fallback(1:nl_soil)
      real(r8), intent(out) :: component_fraction(N_METHANE_COMP)
      real(r8), intent(out) :: lai_component(N_METHANE_COMP)
      real(r8), intent(out) :: crootfr_component(1:nl_soil,N_METHANE_COMP)
      real(r8), intent(out) :: rr_component(N_METHANE_COMP)
      real(r8), intent(out) :: agnpp_component(N_METHANE_COMP)
      real(r8), intent(out) :: bgnpp_component(N_METHANE_COMP)
      real(r8), intent(out) :: annsum_npp_component(N_METHANE_COMP)
      logical, intent(out) :: ready(N_METHANE_COMP)

      integer :: ps, pe, m, j, component
      real(r8) :: frac, rice_area, total_pft_area, profile_sum
      real(r8) :: agnpp_pft, bgnpp_pft, rr_pft
      logical :: is_rice

      component_fraction = 0._r8
      component_fraction(METHANE_COMP_SOIL) = 1._r8
      lai_component = 0._r8
      crootfr_component = 0._r8
      rr_component = 0._r8
      agnpp_component = 0._r8
      bgnpp_component = 0._r8
      annsum_npp_component = 0._r8
      ready = .false.

      ps = -1
      pe = -1
      IF (allocated(patch_pft_s) .and. allocated(patch_pft_e)) THEN
         IF (ipatch >= 1 .and. ipatch <= size(patch_pft_s) .and. &
             ipatch <= size(patch_pft_e)) THEN
            ps = patch_pft_s(ipatch)
            pe = patch_pft_e(ipatch)
         ENDIF
      ENDIF
      IF (.not. pft_arrays_ready(ps, pe)) RETURN
      IF (.not. allocated(pftclass) .or. .not. allocated(lai_p) .or. &
          .not. allocated(irrig_method_p)) RETURN
      IF (pe > size(pftclass) .or. pe > size(lai_p) .or. pe > size(irrig_method_p)) RETURN

      rice_area = 0._r8
      total_pft_area = 0._r8
      DO m = ps, pe
         ! A partially initialized PFT record must not be renormalized away;
         ! fall back to the already-sanitized patch inputs instead.
         IF (.not. valid_bgc_value(pftfrac(m)) .or. pftfrac(m) < 0._r8 .or. &
             .not. valid_bgc_value(lai_p(m)) .or. lai_p(m) < 0._r8 .or. &
             .not. valid_bgc_value(annsum_npp_p(m)) .or. &
             any(.not. valid_bgc_value(cinput_rootfr_p(1:nl_soil,m))) .or. &
             .not. valid_component_pft_fluxes(m)) RETURN
         frac = pftfrac(m)
         total_pft_area = total_pft_area + frac
         is_rice = ch4_rice_pft_is_paddy(pftclass(m), irrig_method_p(m))
         IF (is_rice) rice_area = rice_area + frac
         component = merge(METHANE_COMP_RICE, METHANE_COMP_SOIL, is_rice)

         agnpp_pft = &
            safe_nonnegative(cpool_to_leafc_p(m)) + safe_nonnegative(cpool_to_leafc_storage_p(m)) + &
            safe_nonnegative(cpool_to_livestemc_p(m)) + safe_nonnegative(cpool_to_livestemc_storage_p(m)) + &
            safe_nonnegative(cpool_to_deadstemc_p(m)) + safe_nonnegative(cpool_to_deadstemc_storage_p(m))
         bgnpp_pft = &
            safe_nonnegative(cpool_to_frootc_p(m)) + safe_nonnegative(cpool_to_frootc_storage_p(m)) + &
            safe_nonnegative(cpool_to_livecrootc_p(m)) + safe_nonnegative(cpool_to_livecrootc_storage_p(m)) + &
            safe_nonnegative(cpool_to_deadcrootc_p(m)) + safe_nonnegative(cpool_to_deadcrootc_storage_p(m))
         rr_pft = safe_nonnegative(froot_mr_p(m)) + safe_nonnegative(cpool_froot_gr_p(m)) + &
            safe_nonnegative(cpool_froot_storage_gr_p(m)) + safe_nonnegative(cpool_livecroot_gr_p(m)) + &
            safe_nonnegative(cpool_livecroot_storage_gr_p(m)) + safe_nonnegative(cpool_deadcroot_gr_p(m)) + &
            safe_nonnegative(cpool_deadcroot_storage_gr_p(m)) + safe_nonnegative(transfer_froot_gr_p(m)) + &
            safe_nonnegative(transfer_livecroot_gr_p(m)) + safe_nonnegative(transfer_deadcroot_gr_p(m))

         lai_component(component) = lai_component(component) + frac * lai_p(m)
         agnpp_component(component) = agnpp_component(component) + frac * agnpp_pft
         bgnpp_component(component) = bgnpp_component(component) + frac * bgnpp_pft
         rr_component(component) = rr_component(component) + frac * rr_pft
         annsum_npp_component(component) = annsum_npp_component(component) + &
            frac * safe_nonnegative(annsum_npp_p(m))
         DO j = 1, nl_soil
            crootfr_component(j,component) = crootfr_component(j,component) + &
               frac * max(cinput_rootfr_p(j,m), 0._r8) * dz_soi(j)
         ENDDO
      ENDDO
      IF (total_pft_area > 1._r8 + 1.e-8_r8) RETURN

      rice_area = min(max(rice_area, 0._r8), 1._r8)
      component_fraction(METHANE_COMP_RICE) = rice_area
      component_fraction(METHANE_COMP_SOIL) = 1._r8 - rice_area

      DO component = 1, N_METHANE_COMP
         IF (component_fraction(component) <= 1.e-14_r8) CYCLE
         lai_component(component) = lai_component(component) / component_fraction(component)
         agnpp_component(component) = agnpp_component(component) / component_fraction(component)
         bgnpp_component(component) = bgnpp_component(component) / component_fraction(component)
         rr_component(component) = rr_component(component) / component_fraction(component)
         annsum_npp_component(component) = annsum_npp_component(component) / component_fraction(component)

         profile_sum = sum(crootfr_component(:,component))
         IF (profile_sum <= 0._r8) THEN
            DO j = 1, nl_soil
               IF (valid_bgc_value(rootfr_fallback(j))) &
                  crootfr_component(j,component) = max(rootfr_fallback(j), 0._r8)
            ENDDO
            profile_sum = sum(crootfr_component(:,component))
         ENDIF
         IF (profile_sum > 0._r8) THEN
            crootfr_component(:,component) = crootfr_component(:,component) / profile_sum
         ELSE
            crootfr_component(1,component) = 1._r8
         ENDIF
         ready(component) = .true.
      ENDDO
   END SUBROUTINE mt_ch4_bgc_component_veg_inputs


   LOGICAL FUNCTION valid_component_pft_fluxes(m)
      integer, intent(in) :: m
      real(r8) :: values(22)

      values = [froot_mr_p(m), &
         cpool_to_leafc_p(m), cpool_to_leafc_storage_p(m), &
         cpool_to_livestemc_p(m), cpool_to_livestemc_storage_p(m), &
         cpool_to_deadstemc_p(m), cpool_to_deadstemc_storage_p(m), &
         cpool_to_frootc_p(m), cpool_to_frootc_storage_p(m), &
         cpool_to_livecrootc_p(m), cpool_to_livecrootc_storage_p(m), &
         cpool_to_deadcrootc_p(m), cpool_to_deadcrootc_storage_p(m), &
         cpool_froot_gr_p(m), cpool_froot_storage_gr_p(m), &
         cpool_livecroot_gr_p(m), cpool_livecroot_storage_gr_p(m), &
         cpool_deadcroot_gr_p(m), cpool_deadcroot_storage_gr_p(m), &
         transfer_froot_gr_p(m), transfer_livecroot_gr_p(m), transfer_deadcroot_gr_p(m)]
      valid_component_pft_fluxes = all(valid_bgc_value(values))
   END FUNCTION valid_component_pft_fluxes


   SUBROUTINE methane_check_anoxia_contract (o_scalar_patch)
      real(r8), intent(in) :: o_scalar_patch(1:nl_soil)

      ! Contract with the BGC decomposition side:
      ! * current CoLM202X BGC exposes a stub o_scalar==1, so CH4 applies its
      !   own seasonal inundation factor (use_ch4_sif=.true.).
      ! * if BGC later provides a real anoxia limiter (o_scalar<1), BGC must
      !   own that decomp limitation and CH4 must not apply the seasonal
      !   inundation factor a second time.
      IF (minval(o_scalar_patch) < 1._r8 - O_SCALAR_CONTRACT_TOL) THEN
         IF ((.not. DEF_METHANE%bgc_anoxia_limits_decomp) .or. DEF_METHANE%use_ch4_sif) THEN
            CALL CoLM_stop (' ***** ERROR: CH4/BGC anoxia contract violation: set ' // &
               'bgc_anoxia_limits_decomp=.true. and use_ch4_sif=.false. when BGC o_scalar<1')
         ENDIF
      ENDIF

   END SUBROUTINE methane_check_anoxia_contract


   SUBROUTINE methane_check_hr_classification_contract (unclassified_hr)
      real(r8), intent(in) :: unclassified_hr

      ! CH4 production partitions base_decomp=(somhr+lithr) by hr_vr. If any
      ! positive decomp_hr_vr donor is neither litter/CWD nor soil, then hr_vr
      ! contains carbon that is absent from somhr+lithr, making partition_z
      ! over-normalize the classified substrate and inflate CH4 production.
      IF (unclassified_hr > 1.e-20_r8) THEN
         CALL CoLM_stop (' ***** ERROR: CH4/BGC HR classification contract violation: '// &
            'positive decomp_hr_vr has donor outside soil/litter/CWD; somhr+lithr partition is unsafe')
      ENDIF

   END SUBROUTINE methane_check_hr_classification_contract


   REAL(r8) FUNCTION safe_nonnegative (x)
      real(r8), intent(in) :: x

      IF (valid_bgc_value(x)) THEN
         safe_nonnegative = max(0._r8, x)
      ELSE
         safe_nonnegative = 0._r8
      ENDIF
   END FUNCTION safe_nonnegative


   LOGICAL FUNCTION pft_arrays_ready (ps, pe)
      integer, intent(in) :: ps, pe

      pft_arrays_ready = .false.
      IF (ps <= 0 .or. pe < ps) RETURN
	      IF (.not. allocated(pftfrac)) RETURN
	      IF (.not. allocated(annsum_npp_p)) RETURN
	      IF (.not. allocated(cinput_rootfr_p)) RETURN
      IF (.not. allocated(froot_mr_p)) RETURN
      IF (.not. allocated(cpool_to_leafc_p)) RETURN
      IF (.not. allocated(cpool_to_leafc_storage_p)) RETURN
      IF (.not. allocated(cpool_to_livestemc_p)) RETURN
      IF (.not. allocated(cpool_to_livestemc_storage_p)) RETURN
      IF (.not. allocated(cpool_to_deadstemc_p)) RETURN
      IF (.not. allocated(cpool_to_deadstemc_storage_p)) RETURN
      IF (.not. allocated(cpool_to_frootc_p)) RETURN
      IF (.not. allocated(cpool_to_frootc_storage_p)) RETURN
      IF (.not. allocated(cpool_to_livecrootc_p)) RETURN
      IF (.not. allocated(cpool_to_livecrootc_storage_p)) RETURN
      IF (.not. allocated(cpool_to_deadcrootc_p)) RETURN
      IF (.not. allocated(cpool_to_deadcrootc_storage_p)) RETURN
      IF (.not. allocated(cpool_froot_gr_p)) RETURN
      IF (.not. allocated(cpool_froot_storage_gr_p)) RETURN
      IF (.not. allocated(cpool_livecroot_gr_p)) RETURN
      IF (.not. allocated(cpool_livecroot_storage_gr_p)) RETURN
      IF (.not. allocated(cpool_deadcroot_gr_p)) RETURN
      IF (.not. allocated(cpool_deadcroot_storage_gr_p)) RETURN
      IF (.not. allocated(transfer_froot_gr_p)) RETURN
      IF (.not. allocated(transfer_livecroot_gr_p)) RETURN
      IF (.not. allocated(transfer_deadcroot_gr_p)) RETURN

	      IF (pe > size(pftfrac)) RETURN
	      IF (pe > size(annsum_npp_p)) RETURN
	      IF (nl_soil > size(cinput_rootfr_p,1)) RETURN
	      IF (pe > size(cinput_rootfr_p,2)) RETURN
      IF (pe > size(froot_mr_p)) RETURN
      IF (pe > size(cpool_to_leafc_p)) RETURN
      IF (pe > size(cpool_to_leafc_storage_p)) RETURN
      IF (pe > size(cpool_to_livestemc_p)) RETURN
      IF (pe > size(cpool_to_livestemc_storage_p)) RETURN
      IF (pe > size(cpool_to_deadstemc_p)) RETURN
      IF (pe > size(cpool_to_deadstemc_storage_p)) RETURN
      IF (pe > size(cpool_to_frootc_p)) RETURN
      IF (pe > size(cpool_to_frootc_storage_p)) RETURN
      IF (pe > size(cpool_to_livecrootc_p)) RETURN
      IF (pe > size(cpool_to_livecrootc_storage_p)) RETURN
      IF (pe > size(cpool_to_deadcrootc_p)) RETURN
      IF (pe > size(cpool_to_deadcrootc_storage_p)) RETURN
      IF (pe > size(cpool_froot_gr_p)) RETURN
      IF (pe > size(cpool_froot_storage_gr_p)) RETURN
      IF (pe > size(cpool_livecroot_gr_p)) RETURN
      IF (pe > size(cpool_livecroot_storage_gr_p)) RETURN
      IF (pe > size(cpool_deadcroot_gr_p)) RETURN
      IF (pe > size(cpool_deadcroot_storage_gr_p)) RETURN
      IF (pe > size(transfer_froot_gr_p)) RETURN
      IF (pe > size(transfer_livecroot_gr_p)) RETURN
      IF (pe > size(transfer_deadcroot_gr_p)) RETURN

      pft_arrays_ready = .true.
   END FUNCTION pft_arrays_ready


	   ELEMENTAL LOGICAL FUNCTION valid_bgc_value (x)
	      real(r8), intent(in) :: x

	      valid_bgc_value = (.not. ieee_is_nan(x)) .and. &
         abs(x) < 0.5_r8 * abs(spval)
	   END FUNCTION valid_bgc_value

   !---------------------------------------------------------------------------
   SUBROUTINE get_wetland_veg_proxy (dlat, cellorg_top, lai_in, ipatch, &
      lai_out, annsum_npp_out, agnpp_out, bgnpp_out, rootfr_out)
   !
   ! Climate-zone wetland vegetation proxy for patchtype==2 patches.
   !
   ! IGBP class 11 ("permanent wetland") has no PFT attached in CoLM, so
   ! NPP / rootfr arrive as 0 and methane_aere (Wania 2010 aerenchyma
   ! transport) cannot fire fully -- losing the plant-mediated CH4
   ! pathway that observations show carries 50-90% of total wetland
   ! CH4 efflux (Bridgham et al. 2013 GCB).
   !
   ! LAI handling: CoLM mksrfdata does aggregate the Yuan+2011 (MODIS)
   ! LAI dataset onto wetland patches (fveg0_igbp(11)=1), so the input
   ! lai_in has real data + monthly seasonality.  Trust it when valid
   ! (0.1 < lai_in < 20).  Fall back to a climate-zone peak only when
   ! lai_in is missing (0 or sentinel), e.g. in regions where MODIS LAI
   ! has no signal over a wetland pixel.
   !
   ! NPP / rootfr always injected (no data source for wetland patches).
   !
   ! Climate zones:
   !   - |lat| <= 23.5 deg           -> tropical emergent grass/sedge
   !                                    (varzea, papyrus, Phragmites)
   !   - lat  >  50 deg AND high SOC -> boreal/arctic Sphagnum bog
   !                                    (LAI fallback=0 disables aerenchyma)
   !   - |lat| >  50 deg             -> boreal/subarctic sedge fen (Carex)
   !   - otherwise                   -> temperate marsh
   !
   ! Refs:
   !   Bridgham et al. 2013 GCB (wetland CH4 review)
   !   Wania et al. 2010 GMD   (Sphagnum/sedge parameterisation)
   !   Yuan et al. 2011 RSE    (CoLM LAI dataset)
   !   Bona et al. 2020 GCB    (boreal peatland PFT)
   !   Jackson et al. 1996 Oecologia (global root distributions)
   !---------------------------------------------------------------------------
   real(r8), intent(in)  :: dlat            ! patch latitude [deg, signed]
   real(r8), intent(in)  :: cellorg_top     ! top-soil organic matter density [kg OM/m3]
   real(r8), intent(in)  :: lai_in          ! input LAI from CoLM data [m2/m2] (may be 0/spval)
   integer,  intent(in)  :: ipatch          ! patch index, used to write override arrays
   real(r8), intent(out) :: lai_out         ! effective LAI for aerenchyma [m2/m2]
   real(r8), intent(out) :: annsum_npp_out  ! annual NPP [gC/m2/yr]
   real(r8), intent(out) :: agnpp_out       ! aboveground NPP rate [gC/m2/s]
   real(r8), intent(out) :: bgnpp_out       ! belowground NPP rate [gC/m2/s]
   real(r8), intent(out) :: rootfr_out(1:nl_soil)

   real(r8) :: lai_fallback, anpp_tot, bg_ratio
   real(r8) :: aere_poros_z, aere_radius_z, aere_tillerC_z, aere_scale_z
   real(r8), parameter :: secsperyear = 365._r8 * 86400._r8
   ! Histosol-style peat threshold: cellorg ~ 150 kg OM/m3 corresponds to
   ! ~75 kgC/m3, well above typical mineral soil (< 50 kg OM/m3).
   real(r8), parameter :: peat_om_threshold = 150._r8
   ! Tropical wetland discriminator: high SOC indicates peat-forming
   ! emergent vasculars (papyrus, Phragmites); low SOC indicates mineral
   ! soil under varzea grass / swamp forest.
   real(r8), parameter :: tropical_peat_threshold = 80._r8
   ! Bounds for trusting lai_in: 0.1 filters out empty cells; 20 filters
   ! spval (~1e36) and pathological aggregation outputs.
   real(r8), parameter :: lai_in_min = 0.1_r8, lai_in_max = 20._r8

   ! ---- 5-zone climate-based wetland classification --------------------
   ! Zone parameter values from:
   !   Brix 2001 Aquat Bot (Phragmites)
   !   Wania 2010 GMD     (Carex sedge defaults)
   !   Saunders 2007      (Papyrus)
   !   Bridgham 2013 GCB  (review of plant-mediated CH4)
   !   Pangala 2017 Nature(tropical swamp forest stem conduits)
   IF (abs(dlat) <= 23.5_r8 .and. cellorg_top >= tropical_peat_threshold) THEN
      ! Zone 1: Tropical reed/papyrus (peat-rich, large emergent vascular)
      lai_fallback   = 4.0_r8
      anpp_tot       = 600._r8
      bg_ratio       = 0.5_r8
      aere_poros_z   = 0.50_r8       ! Phragmites/papyrus highly developed aerenchyma
      aere_radius_z  = 8.0e-3_r8     ! ~8 mm stems
      aere_tillerC_z = 3.0_r8        ! ~3 gC/tiller (large emergent)
      aere_scale_z   = 1.0_r8
   ELSEIF (abs(dlat) <= 23.5_r8) THEN
      ! Zone 2: Tropical swamp (mineral soil, varzea grass / igapo forest mix).
      ! Smaller tillers + lower porosity since trees and grasses share patch;
      ! Pangala 2017 stem-conduit fluxes not yet modelled here.
      lai_fallback   = 3.0_r8
      anpp_tot       = 400._r8
      bg_ratio       = 0.3_r8        ! higher AG biomass for tree mix
      aere_poros_z   = 0.10_r8       ! reduced (trees + grass mix)
      aere_radius_z  = 3.0e-3_r8     ! grass scale
      aere_tillerC_z = 0.5_r8
      aere_scale_z   = 0.5_r8        ! half default tiller transport
   ELSEIF (abs(dlat) > 50._r8 .and. cellorg_top > peat_om_threshold) THEN
      ! Zone 5: Boreal/arctic Sphagnum bog (non-vascular moss).
      ! LAI fallback = 0 disables aerenchyma entirely; override params
      ! still set to 0 for safety in case lai_in trips the data branch.
      lai_fallback   = 0._r8
      anpp_tot       = 100._r8
      bg_ratio       = 0.2_r8
      aere_poros_z   = 0._r8
      aere_radius_z  = 1.e-6_r8      ! nonzero to avoid div-by-zero downstream
      aere_tillerC_z = 1._r8         ! nonzero to avoid div-by-zero
      aere_scale_z   = 0._r8
   ELSEIF (abs(dlat) > 50._r8) THEN
      ! Zone 4: Boreal/subarctic sedge fen (Carex-dominated, CTSM defaults)
      lai_fallback   = 2.0_r8
      anpp_tot       = 200._r8
      bg_ratio       = 0.6_r8
      aere_poros_z   = 0.30_r8
      aere_radius_z  = 2.9e-3_r8
      aere_tillerC_z = 0.3_r8
      aere_scale_z   = 1.0_r8
   ELSE
      ! Zone 3: Temperate marsh (Typha + Carex mix)
      lai_fallback   = 2.5_r8
      anpp_tot       = 300._r8
      bg_ratio       = 0.5_r8
      aere_poros_z   = 0.40_r8
      aere_radius_z  = 5.0e-3_r8
      aere_tillerC_z = 1.0_r8
      aere_scale_z   = 1.0_r8
   ENDIF

   ! Write per-patch aerenchyma override (consumed by methane_aere /
   ! SiteOxAere in Physics via getter functions).
   IF (allocated(wetland_aere_active) .and. &
       ipatch >= 1 .and. ipatch <= size(wetland_aere_active)) THEN
      wetland_aere_poros  (ipatch) = aere_poros_z
      wetland_aere_radius (ipatch) = aere_radius_z
      wetland_aere_tillerC(ipatch) = aere_tillerC_z
      wetland_aere_scale  (ipatch) = aere_scale_z
      wetland_aere_active (ipatch) = .true.
   ENDIF

   ! LAI: trust data when valid, else fall back to climate-zone peak.
   IF (lai_in > lai_in_min .and. lai_in < lai_in_max) THEN
      lai_out = lai_in
   ELSE
      lai_out = lai_fallback
   ENDIF

   annsum_npp_out = anpp_tot
   agnpp_out      = (1._r8 - bg_ratio) * anpp_tot / secsperyear
   bgnpp_out      =          bg_ratio  * anpp_tot / secsperyear

   ! Sedge-style shallow root profile, normalized.  ~80% in top 30 cm.
   rootfr_out(:) = 0._r8
   IF (nl_soil >= 1) rootfr_out(1) = 0.35_r8
   IF (nl_soil >= 2) rootfr_out(2) = 0.25_r8
   IF (nl_soil >= 3) rootfr_out(3) = 0.15_r8
   IF (nl_soil >= 4) rootfr_out(4) = 0.10_r8
   IF (nl_soil >= 5) rootfr_out(5) = 0.08_r8
   IF (nl_soil >= 6) rootfr_out(6) = 0.04_r8
   IF (nl_soil >= 7) rootfr_out(7) = 0.02_r8
   IF (nl_soil >= 8) rootfr_out(8) = 0.01_r8
   IF (sum(rootfr_out) > 0._r8) rootfr_out(:) = rootfr_out(:) / sum(rootfr_out)

   END SUBROUTINE get_wetland_veg_proxy

   !---------------------------------------------------------------------------
   SUBROUTINE get_rice_veg_proxy (lai_in, ipatch, rice_fraction)
   !
   ! Zone 6: Rice paddy aerenchyma (R4 extension; off unless DEF_METHANE%enable_rice_paddy).
   !
   ! Writes wetland_aere_* overrides for paddy-irrigated rice patches so
   ! methane_aere uses rice-tuned tiller geometry instead of the
   ! Wania-natural-wetland default (poros=0.30, radius=2.9mm).
   !
   ! Rice tiller params (literature):
   !   Justin & Armstrong 1991 (New Phytol)  : aerenchyma porosity 0.18-0.40
   !   Sasakawa & Yamamoto 1978 (Soil Sci PN): tiller base diameter ~1.5 mm
   !   Aulakh+2001 Plant&Soil               : tiller density 200-400/m2
   !   Holzapfel-Pschorn+1985 Plant&Soil    : plant-mediated CH4 ~30-50% total
   !
   ! Unlike get_wetland_veg_proxy this does NOT overwrite NPP / rootfr / LAI —
   ! CN already computes those correctly for rice CFTs (61/62).  Only the
   ! aerenchyma transport geometry needs override.
   !---------------------------------------------------------------------------
   real(r8), intent(in)  :: lai_in          ! input LAI from CoLM/CN [m2/m2]
   integer,  intent(in)  :: ipatch          ! patch index
   real(r8), intent(in)  :: rice_fraction   ! paddy rice area fraction in mixed soil patch

   real(r8) :: aere_poros_z, aere_radius_z, aere_tillerC_z, aere_scale_z

   ! Rice paddy fixed params.  aere_radius is treated as a RADIUS by
   ! SiteOxAere (area_tiller ~ PI * radius_eff**2).  Watanabe 2001
   ! / Aulakh 2001 report rice tiller cross-section ~1.5 mm in DIAMETER
   ! at the base, so the radius is ~0.75 mm; mis-reading diameter as
   ! radius would over-estimate aerenchyma cross-section by 4x.
   aere_poros_z   = 0.40_r8       ! Justin+Armstrong 1991 typical, more developed than sedge
   aere_radius_z  = 0.75e-3_r8    ! 0.75 mm radius = 1.5 mm tiller diameter (Watanabe 2001)
   aere_tillerC_z = 1.0_r8        ! denser tillers vs natural wetland 0.22 default

   ! Phenology proxy: scale aerenchyma area by relative LAI within the
   ! growing season.  LAI ~ 4 at peak tillering/heading for rice (IRRI
   ! N2/N3 cv).  Cap [0.5, 1.5] guards against early-season LAI<<1 from
   ! transplant transient or very dense double-crop LAI > 6 inflating
   ! the override.  The caller (methane_driver) only invokes this routine
   ! when is_paddy_rice_live(ipatch)==.true., so post-harvest stubble is
   ! NOT a concern here — the override simply isn't applied after harvest.
   IF (lai_in > 0._r8 .and. lai_in < 20._r8) THEN
      aere_scale_z = max(0.5_r8, min(1.5_r8, lai_in / 4._r8))
   ELSE
      aere_scale_z = 1.0_r8
   ENDIF
   ! Rice is represented as a CFT fraction inside a soil patch, not as a
   ! separate land patch.  Scale the rice aerenchyma override by the paddy
   ! fraction so a small rice tile cannot impose rice transport geometry on
   ! the whole mixed soil patch.
   aere_scale_z = aere_scale_z * min(max(rice_fraction, 0._r8), 1._r8)

   IF (allocated(wetland_aere_active) .and. &
       ipatch >= 1 .and. ipatch <= size(wetland_aere_active)) THEN
      wetland_aere_poros  (ipatch) = aere_poros_z
      wetland_aere_radius (ipatch) = aere_radius_z
      wetland_aere_tillerC(ipatch) = aere_tillerC_z
      wetland_aere_scale  (ipatch) = aere_scale_z
      wetland_aere_active (ipatch) = .true.
   ENDIF

   END SUBROUTINE get_rice_veg_proxy

   !---------------------------------------------------------------------------
   ! Returns .true. for rice CFTs that CH4 should treat as paddy rice.
   !
   ! Important: do not require irrig_method_paddy only.  The crop input stores
   ! rice paddy water management as irrig_method_flood (=3), and MOD_Irrigation
   ! converts flood -> paddy only when DEF_USE_IRRIGATION is true.  CH4 rice
   ! paddy mode is controlled by DEF_METHANE%enable_rice_paddy and must still
   ! recognise flood-managed rice when irrigation physics is disabled.
   !---------------------------------------------------------------------------
   logical function ch4_rice_pft_is_paddy(pft_class, irrig_method) result(yes)
#ifdef CROP
      USE MOD_Vars_Global, only: nrice, nirrig_rice, irrig_method_flood, irrig_method_paddy
#endif
      integer, intent(in) :: pft_class, irrig_method

      yes = .false.
#ifdef CROP
      IF (pft_class == nrice .or. pft_class == nirrig_rice) THEN
         yes = (irrig_method == irrig_method_paddy .or. irrig_method == irrig_method_flood)
      ENDIF
#endif
   END FUNCTION ch4_rice_pft_is_paddy

   !---------------------------------------------------------------------------
   ! Returns .true. if the patch contains at least one CH4 paddy rice tile
   ! (flood/paddy-managed CFT nrice / nirrig_rice) that CN reports as alive
   ! (croplive_p).
   ! This is the single place that touches BGC state from the methane chain,
   ! satisfying the architectural rule that MOD_MT_Methane_Driver and
   ! MOD_MT_Methane_Physics access BGC only through this module.
   ! Cross-year wrap is handled by CN's state machine (croplive_p stays
   ! .true. across the Dec->Jan boundary for southern-hemisphere rice).
   !---------------------------------------------------------------------------
   logical function is_paddy_rice_live(ipatch) result(yes)
#ifdef CROP
      USE MOD_BGC_Vars_PFTimeVariables, only: croplive_p
      USE MOD_LandPFT,             only: patch_pft_s, patch_pft_e
      USE MOD_Vars_PFTimeVariables,  only: irrig_method_p
      USE MOD_Vars_PFTimeInvariants, only: pftclass
      integer :: m_loc, ps_loc, pe_loc, n_pft_max
#endif
      integer, intent(in) :: ipatch

      yes = .false.
#ifdef CROP
      IF (.not. allocated(patch_pft_s) .or. .not. allocated(patch_pft_e)) RETURN
      IF (.not. allocated(pftclass) .or. .not. allocated(irrig_method_p)) RETURN
      IF (.not. allocated(croplive_p)) RETURN
      IF (ipatch < 1 .or. ipatch > size(patch_pft_s) .or. ipatch > size(patch_pft_e)) RETURN
      ps_loc = patch_pft_s(ipatch)
      pe_loc = patch_pft_e(ipatch)
      IF (ps_loc < 1 .or. pe_loc < ps_loc) RETURN
      n_pft_max = min(size(pftclass), size(irrig_method_p), size(croplive_p))
      IF (pe_loc > n_pft_max) RETURN
      DO m_loc = ps_loc, pe_loc
         IF (ch4_rice_pft_is_paddy(pftclass(m_loc), irrig_method_p(m_loc)) .and. &
             croplive_p(m_loc)) THEN
            yes = .true.
            RETURN
         ENDIF
      ENDDO
#endif
   END FUNCTION is_paddy_rice_live

   !---------------------------------------------------------------------------
   ! Returns .true. if the patch HOSTS a CH4 paddy rice (flood/paddy-managed
   ! CFT nrice / nirrig_rice) tile, regardless of CN croplive_p.  Mirrors the driver
   ! gate at CoLMDRIVER.F90:380 (rice_pft_frac > 0.01 -> run_methane=.true.)
   ! exactly: methane() runs on any patch with a paddy rice CFT, even
   ! before plant or after harvest, so the history/AccFlux mask must
   ! include those patches too.  Distinct from is_paddy_rice_live() which
   ! is used for the season-specific finundated override timing.
   !---------------------------------------------------------------------------
   !---------------------------------------------------------------------------
   ! Sum of pftfrac over CH4 paddy rice (flood/paddy-managed CFT nrice /
   ! nirrig_rice) tiles in this patch.  Returns 0.0 when patch is out of range, when
   ! BGC/CROP arrays are not allocated, or when no paddy rice tile is
   ! present.  Clamped to [0, 1] to guard against upstream pftfrac
   ! normalization slips.  THIS is the single source of truth for the
   ! "how much of this patch is paddy rice" question; everything else
   ! (driver gate, mask helper) goes through it.
   !---------------------------------------------------------------------------
   real(r8) function paddy_rice_fraction(ipatch) result(frac)
#ifdef CROP
      USE MOD_LandPFT,             only: patch_pft_s, patch_pft_e
      USE MOD_Vars_PFTimeVariables,  only: irrig_method_p
      USE MOD_Vars_PFTimeInvariants, only: pftclass, pftfrac
      integer :: m_loc, ps_loc, pe_loc, n_pft_max
#endif
      integer, intent(in) :: ipatch

      frac = 0._r8
#ifdef CROP
      IF (.not. allocated(patch_pft_s) .or. .not. allocated(patch_pft_e)) RETURN
      IF (.not. allocated(pftclass) .or. .not. allocated(irrig_method_p)) RETURN
      IF (.not. allocated(pftfrac)) RETURN
      IF (ipatch < 1 .or. ipatch > size(patch_pft_s) .or. ipatch > size(patch_pft_e)) RETURN
      ps_loc = patch_pft_s(ipatch)
      pe_loc = patch_pft_e(ipatch)
      IF (ps_loc < 1 .or. pe_loc < ps_loc) RETURN
      n_pft_max = min(size(pftclass), size(irrig_method_p), size(pftfrac))
      IF (pe_loc > n_pft_max) RETURN
      DO m_loc = ps_loc, pe_loc
         IF (ch4_rice_pft_is_paddy(pftclass(m_loc), irrig_method_p(m_loc))) THEN
            frac = frac + pftfrac(m_loc)
         ENDIF
      ENDDO
      frac = min(max(frac, 0._r8), 1._r8)
#endif
   END FUNCTION paddy_rice_fraction

   !---------------------------------------------------------------------------
   ! Convenience wrapper: patch hosts paddy rice tile above the gate
   ! threshold (matches driver gate at CoLMDRIVER.F90).
   !---------------------------------------------------------------------------
   logical function has_paddy_rice_tile(ipatch) result(yes)
      integer, intent(in) :: ipatch
      yes = paddy_rice_fraction(ipatch) > PADDY_RICE_FRAC_MIN
   END FUNCTION has_paddy_rice_tile

   !---------------------------------------------------------------------------
   ! Days past planting / since harvest for the first paddy rice CFT in the
   ! patch.  A wrapped event belongs to the previous calendar year, so its
   ! year length can differ from the current year across a leap-year boundary.
   ! Returns -1 when info is unavailable (no paddy rice,
   ! harvdate=NOT_Harvested, invalid day-of-year, etc).
   !---------------------------------------------------------------------------
   integer function rice_days_past_planting(ipatch, jday, year) result(idpp)
#ifdef CROP
      USE MOD_BGC_Vars_PFTimeVariables, only: croplive_p, idop_p
      USE MOD_LandPFT,             only: patch_pft_s, patch_pft_e
      USE MOD_Vars_PFTimeVariables,  only: irrig_method_p
      USE MOD_Vars_PFTimeInvariants, only: pftclass
      integer :: m_loc, ps_loc, pe_loc, n_pft_max
#endif
      integer, intent(in) :: ipatch, jday
      integer, intent(in), optional :: year
      integer :: dayspyr, event_dayspyr

      idpp = -1
      dayspyr = 365
      IF (present(year)) THEN
         IF (isleapyear(year)) dayspyr = 366
      ENDIF
      IF (jday < 1 .or. jday > dayspyr) RETURN
#ifdef CROP
      IF (.not. allocated(patch_pft_s) .or. .not. allocated(patch_pft_e)) RETURN
      IF (.not. allocated(pftclass) .or. .not. allocated(irrig_method_p)) RETURN
      IF (.not. allocated(croplive_p) .or. .not. allocated(idop_p)) RETURN
      IF (ipatch < 1 .or. ipatch > size(patch_pft_s) .or. ipatch > size(patch_pft_e)) RETURN
      ps_loc = patch_pft_s(ipatch)
      pe_loc = patch_pft_e(ipatch)
      IF (ps_loc < 1 .or. pe_loc < ps_loc) RETURN
      n_pft_max = min(size(pftclass), size(irrig_method_p), &
                      size(croplive_p), size(idop_p))
      IF (pe_loc > n_pft_max) RETURN
      DO m_loc = ps_loc, pe_loc
         IF (ch4_rice_pft_is_paddy(pftclass(m_loc), irrig_method_p(m_loc)) .and. &
             croplive_p(m_loc)) THEN
            IF (idop_p(m_loc) >= 1) THEN
               IF (jday >= idop_p(m_loc) .and. idop_p(m_loc) <= dayspyr) THEN
                  idpp = jday - idop_p(m_loc)
               ELSE
                  event_dayspyr = 365
                  IF (present(year)) THEN
                     IF (isleapyear(year - 1)) event_dayspyr = 366
                  ENDIF
                  IF (idop_p(m_loc) > event_dayspyr) CYCLE
                  idpp = event_dayspyr + jday - idop_p(m_loc)
               ENDIF
               RETURN
            ENDIF
         ENDIF
      ENDDO
#endif
   END FUNCTION rice_days_past_planting

   integer function rice_days_since_harvest(ipatch, jday, year) result(dsh)
#ifdef CROP
      USE MOD_BGC_Vars_PFTimeVariables, only: harvdate_p
      USE MOD_LandPFT,             only: patch_pft_s, patch_pft_e
      USE MOD_Vars_PFTimeVariables,  only: irrig_method_p
      USE MOD_Vars_PFTimeInvariants, only: pftclass
      integer :: m_loc, ps_loc, pe_loc, n_pft_max
      integer, parameter :: NOT_Harvested = 999
#endif
      integer, intent(in) :: ipatch, jday
      integer, intent(in), optional :: year
      integer :: dayspyr, event_dayspyr

      dsh = -1
      dayspyr = 365
      IF (present(year)) THEN
         IF (isleapyear(year)) dayspyr = 366
      ENDIF
      IF (jday < 1 .or. jday > dayspyr) RETURN
#ifdef CROP
      IF (.not. allocated(patch_pft_s) .or. .not. allocated(patch_pft_e)) RETURN
      IF (.not. allocated(pftclass) .or. .not. allocated(irrig_method_p)) RETURN
      IF (.not. allocated(harvdate_p)) RETURN
      IF (ipatch < 1 .or. ipatch > size(patch_pft_s) .or. ipatch > size(patch_pft_e)) RETURN
      ps_loc = patch_pft_s(ipatch)
      pe_loc = patch_pft_e(ipatch)
      IF (ps_loc < 1 .or. pe_loc < ps_loc) RETURN
      n_pft_max = min(size(pftclass), size(irrig_method_p), size(harvdate_p))
      IF (pe_loc > n_pft_max) RETURN
      DO m_loc = ps_loc, pe_loc
         IF (ch4_rice_pft_is_paddy(pftclass(m_loc), irrig_method_p(m_loc))) THEN
            IF (harvdate_p(m_loc) >= 1 .and. harvdate_p(m_loc) < NOT_Harvested) THEN
               IF (jday >= harvdate_p(m_loc) .and. harvdate_p(m_loc) <= dayspyr) THEN
                  dsh = jday - harvdate_p(m_loc)
               ELSE
                  event_dayspyr = 365
                  IF (present(year)) THEN
                     IF (isleapyear(year - 1)) event_dayspyr = 366
                  ENDIF
                  IF (harvdate_p(m_loc) > event_dayspyr) CYCLE
                  dsh = event_dayspyr + jday - harvdate_p(m_loc)
               ENDIF
               RETURN
            ENDIF
         ENDIF
      ENDDO
#endif
   END FUNCTION rice_days_since_harvest

   !---------------------------------------------------------------------------
   ! Build the per-patch logical mask of patches where methane() actually
   ! runs.  Same gate used by CoLMDRIVER.F90 (~line 380) so AccFlux and
   ! history filters stay in sync with the physics.
   !
   ! Active if any of:
   !   - patchtype == 2 (natural wetland)
   !   - patchtype == 0 (soil) and DEF_METHANE%only_wetland is false
   !   - patchtype == 0 (soil) and rice paddy mode enabled AND patch
   !     hosts a paddy rice tile (uses has_paddy_rice_tile NOT
   !     is_paddy_rice_live so that pre-plant and post-harvest paddy
   !     patches — where methane() still runs for drainage CH4 pulse —
   !     are included, matching the driver gate exactly).
   !---------------------------------------------------------------------------
   SUBROUTINE methane_patch_active_mask (numpatch_loc, only_wetland, &
                                          enable_rice_paddy, ptype, mask)
      integer, intent(in)  :: numpatch_loc
      logical, intent(in)  :: only_wetland, enable_rice_paddy
      integer, intent(in)  :: ptype(numpatch_loc)
      logical, intent(out) :: mask(numpatch_loc)
      integer :: i

      DO i = 1, numpatch_loc
         mask(i) = (ptype(i) == 2)
         IF (.not. only_wetland .and. ptype(i) == 0) mask(i) = .true.
#ifdef CROP
         IF (enable_rice_paddy .and. ptype(i) == 0) THEN
            IF (has_paddy_rice_tile(i)) mask(i) = .true.
         ENDIF
#endif
      ENDDO
   END SUBROUTINE methane_patch_active_mask

   !---------------------------------------------------------------------------
   real(r8) FUNCTION get_biome_f_methane (patchtype, dlat, cellorg_top, &
      is_rice_paddy, rice_fraction, rice_parameter_active, is_floodplain)
   !
   ! Biome-specific f_methane (CH4 yield ratio) lookup.
   ! Mirrors get_wetland_veg_proxy 5-zone classification so each climate
   ! zone gets a literature-informed CH4/HR ratio instead of CTSM's global
   ! 0.20 default.
   !
   ! IMPORTANT: The specific f_methane numbers in DEF_METHANE are author-
   ! selected midpoints within published literature RANGES, not direct
   ! quotations from any single paper.  Tuned against Pantanal in-situ
   ! (Marani & Alvalá 2007).  Other biomes use range midpoints.
   !
   ! Primary literature (full bibliographic detail; see also Const.F90 header):
   !   Bridgham, S. D., Cadillo-Quiroz, H., Keller, J. K., & Zhuang, Q. (2013).
   !     Methane emissions from wetlands: biogeochemical, microbial, and
   !     modeling perspectives from local to global scales.
   !     Global Change Biology, 19(5), 1325-1346.  doi:10.1111/gcb.12131
   !   Marani, L., & Alvalá, P. C. (2007).  Methane emissions from lakes and
   !     floodplains in Pantanal, Brazil.
   !     Atmospheric Environment, 41(8), 1627-1633.
   !     doi:10.1016/j.atmosenv.2006.10.046
   !   Wania, R., Ross, I., & Prentice, I. C. (2010).  Implementation and
   !     evaluation of a new methane model within a dynamic global vegetation
   !     model: LPJ-WHyMe v1.3.1.
   !     Geoscientific Model Development, 3(2), 565-584.
   !     doi:10.5194/gmd-3-565-2010
   !   Walter, B. P., Heimann, M., & Matthews, E. (2001).  Modeling modern
   !     methane emissions from natural wetlands: 1. Model description and
   !     results.  J. Geophys. Res. Atmos., 106(D24), 34189-34206.
   !     doi:10.1029/2001JD900165
   !   Holzapfel-Pschorn, A., Conrad, R., & Seiler, W. (1985).  Production,
   !     oxidation and emission of methane in rice paddies.
   !     FEMS Microbiology Ecology, 1(6), 343-351.
   !     doi:10.1111/j.1574-6968.1985.tb01605.x
   !   Le Mer, J., & Roger, P. (2001).  Production, oxidation, emission and
   !     consumption of methane by soils: A review.
   !     European Journal of Soil Biology, 37(1), 25-50.
   !     doi:10.1016/S1164-5563(01)01067-6
   !---------------------------------------------------------------------------
   USE MOD_MT_Methane_Const, only: DEF_METHANE
   integer,  intent(in) :: patchtype
   real(r8), intent(in) :: dlat
   real(r8), intent(in) :: cellorg_top
   logical,  intent(in) :: is_rice_paddy
   real(r8), intent(in) :: rice_fraction
   logical,  intent(in) :: rice_parameter_active
   logical,  intent(in), optional :: is_floodplain

   real(r8), parameter :: peat_om_threshold      = 150._r8
   real(r8), parameter :: tropical_peat_threshold = 80._r8
   logical :: floodplain_active
   real(r8) :: rice_frac
   real(r8) :: nonrice_f_methane

   ! Legacy path: use global scalar if lookup not enabled
   IF (.not. DEF_METHANE%use_biome_f_methane) THEN
      get_biome_f_methane = DEF_METHANE%f_methane
      RETURN
   ENDIF

   floodplain_active = .false.
   IF (present(is_floodplain)) floodplain_active = is_floodplain
   rice_frac = min(max(rice_fraction, 0._r8), 1._r8)

   ! First compute the non-rice background value for this patch.  Rice is a
   ! CFT fraction inside patchtype==0, so it must blend with the non-rice
   ! soil/floodplain value instead of replacing the whole mixed soil patch.
   IF (floodplain_active) THEN
      nonrice_f_methane = DEF_METHANE%f_methane_floodplain
   ELSE IF (patchtype /= 2) THEN
      nonrice_f_methane = DEF_METHANE%f_methane_upland_soil
   ELSE IF (abs(dlat) <= 23.5_r8 .and. cellorg_top >= tropical_peat_threshold) THEN
      nonrice_f_methane = DEF_METHANE%f_methane_tropical_peat
   ELSE IF (abs(dlat) <= 23.5_r8) THEN
      nonrice_f_methane = DEF_METHANE%f_methane_tropical_floodplain
   ELSE IF (abs(dlat) > 50._r8 .and. cellorg_top > peat_om_threshold) THEN
      nonrice_f_methane = DEF_METHANE%f_methane_boreal_bog
   ELSE IF (abs(dlat) > 50._r8) THEN
      nonrice_f_methane = DEF_METHANE%f_methane_boreal_fen
   ELSE
      nonrice_f_methane = DEF_METHANE%f_methane_temperate_marsh
   ENDIF

   IF (is_rice_paddy .and. rice_parameter_active .and. rice_frac > 0._r8) THEN
      get_biome_f_methane = (1._r8 - rice_frac) * nonrice_f_methane + &
                            rice_frac * DEF_METHANE%f_methane_rice_paddy
   ELSE
      get_biome_f_methane = nonrice_f_methane
   ENDIF

   END FUNCTION get_biome_f_methane

   !---------------------------------------------------------------------------
   real(r8) FUNCTION get_biome_redoxlag (patchtype, dlat, cellorg_top, &
      is_rice_paddy, rice_fraction, rice_parameter_active, is_floodplain)
   !
   ! Biome-specific redoxlag (microbial response time, days) lookup.
   ! Mirrors get_wetland_veg_proxy 5-zone classification.  Concept: faster
   ! response in warm tropical wetlands, slower in cold boreal peat.
   !
   ! IMPORTANT: Specific day values in DEF_METHANE are author-chosen
   ! order-of-magnitude estimates based on QUALITATIVE conclusions in the
   ! cited papers, NOT direct numerical quotes.  redoxlag is less well-
   ! constrained empirically than f_methane (wider uncertainty 50-100%).
   !
   ! Primary literature informing the qualitative ordering:
   !   Pangala, S. R., Enrich-Prast, A., Basso, L. S., Peixoto, R. B.,
   !     Bastviken, D., Hornibrook, E. R. C., et al. (2017).  Large
   !     emissions from floodplain trees close the Amazon methane budget.
   !     Nature, 552(7684), 230-234.  doi:10.1038/nature24639
   !     (Amazon tree stem rapid CH4 emission -> tropical rapid response)
   !   Walter, B. P., Heimann, M., & Matthews, E. (2001).  Modeling modern
   !     methane emissions from natural wetlands: 1. Model description and
   !     results.  J. Geophys. Res. Atmos., 106(D24), 34189-34206.
   !     doi:10.1029/2001JD900165
   !     (CTSM-baseline ~14-30 d wetland CH4 transition timescale)
   !   Whalen, S. C., & Reeburgh, W. S. (1990).  Consumption of atmospheric
   !     methane by tundra soils.
   !     Nature, 346(6280), 160-162.  doi:10.1038/346160a0
   !     (cold boreal microbe slow response)
   !   Holzapfel-Pschorn, A., Conrad, R., & Seiler, W. (1985).  Production,
   !     oxidation and emission of methane in rice paddies.
   !     FEMS Microbiology Ecology, 1(6), 343-351.
   !     doi:10.1111/j.1574-6968.1985.tb01605.x
   !     (rice paddy rapid response to managed flood/drain)
   !---------------------------------------------------------------------------
   USE MOD_MT_Methane_Const, only: DEF_METHANE
   integer,  intent(in) :: patchtype
   real(r8), intent(in) :: dlat
   real(r8), intent(in) :: cellorg_top
   logical,  intent(in) :: is_rice_paddy
   real(r8), intent(in) :: rice_fraction
   logical,  intent(in) :: rice_parameter_active
   logical,  intent(in), optional :: is_floodplain

   real(r8), parameter :: peat_om_threshold       = 150._r8
   real(r8), parameter :: tropical_peat_threshold = 80._r8
   logical :: floodplain_active
   real(r8) :: rice_frac
   real(r8) :: nonrice_redoxlag

   ! Legacy path: use global scalar if lookup not enabled
   IF (.not. DEF_METHANE%use_biome_redoxlag) THEN
      get_biome_redoxlag = DEF_METHANE%redoxlag
      RETURN
   ENDIF

   floodplain_active = .false.
   IF (present(is_floodplain)) floodplain_active = is_floodplain
   rice_frac = min(max(rice_fraction, 0._r8), 1._r8)

   ! Compute the non-rice background response time, then blend the paddy
   ! rapid-response parameter by rice area fraction for mixed soil patches.
   IF (floodplain_active) THEN
      nonrice_redoxlag = DEF_METHANE%redoxlag_tropical_floodplain
   ELSE IF (patchtype /= 2) THEN
      nonrice_redoxlag = DEF_METHANE%redoxlag_upland_soil
   ELSE IF (abs(dlat) <= 23.5_r8 .and. cellorg_top >= tropical_peat_threshold) THEN
      nonrice_redoxlag = DEF_METHANE%redoxlag_tropical_peat
   ELSE IF (abs(dlat) <= 23.5_r8) THEN
      nonrice_redoxlag = DEF_METHANE%redoxlag_tropical_floodplain
   ELSE IF (abs(dlat) > 50._r8 .and. cellorg_top > peat_om_threshold) THEN
      nonrice_redoxlag = DEF_METHANE%redoxlag_boreal_bog
   ELSE IF (abs(dlat) > 50._r8) THEN
      nonrice_redoxlag = DEF_METHANE%redoxlag_boreal_fen
   ELSE
      nonrice_redoxlag = DEF_METHANE%redoxlag_temperate_marsh
   ENDIF

   IF (is_rice_paddy .and. rice_parameter_active .and. rice_frac > 0._r8) THEN
      get_biome_redoxlag = (1._r8 - rice_frac) * nonrice_redoxlag + &
                           rice_frac * DEF_METHANE%redoxlag_rice_paddy
   ELSE
      get_biome_redoxlag = nonrice_redoxlag
   ENDIF

   END FUNCTION get_biome_redoxlag

END MODULE MOD_MT_Methane_BgcLink
#endif
