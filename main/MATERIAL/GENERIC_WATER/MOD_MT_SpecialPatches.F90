#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_SpecialPatches

   USE MOD_Precision
   USE MOD_SPMD_Task, only: CoLM_stop
   USE MOD_MT_Defs, only: nmts, materials, mt_tiny, mt_water_min_for_ratio, &
      mt_delta_sanity_max, material_init_water_ratio, material_can_use_fixed_signature, &
      mt_uses_generic_water, material_is_nonvolatile_solute, material_is_isotope
   USE MOD_MT_Solute, only: mt_has_dissolved_limit, mt_equilibrate_dissolved
   USE MOD_MT_Namelist, only: DEF_ISO_SUBL_SKIN_MM
   USE MOD_MT_Forcing, only: material_forcing_rain_value, material_forcing_snow_value, &
      material_forcing_vapor_value
   USE MOD_MT_EvapLimit, only: material_atmospheric_mt_loss
   USE MOD_MT_Frac, only: material_fractionation_active, material_surface_relhum, &
      material_alpha_kinetic_craig_gordon, material_craig_gordon_evap_affine, &
      material_alpha_kinetic_open_water, &
      material_equilibrium_deposition_ratio
   USE MOD_MT_Conservation, only: mt_save_storage, mt_balance_check, &
      mt_apply_decay
   USE MOD_MT_Hist, only: mt_hist_accumulate
   USE MOD_MT_Vars, only: mt_wliq_soisno, mt_wice_soisno, mt_scv, mt_wdsrf, &
      mt_ldew_rain, mt_ldew_snow, mt_rnof_step, a_mt_precip, material_book_evap_loss, &
      MT_EVAP_KIND_SOILEVAP, MT_EVAP_KIND_SUBL, &
      a_mt_rsur, a_mt_rnof, mt_wetwat, mt_waterstorage, mt_storage_beg, &
      mt_surface_residue, mt_subsurface_residue, mt_runtime_forced, &
      mt_solid_soisno, mt_canopy_solid, mt_surface_solid, &
      mt_subsurface_solid, mt_waterstorage_solid, sync_mt_patch_ratio, &
      mt_add_rnof_step
   USE MOD_MT_Vars, only: mt_aquifer_ref_water, mt_leaf_iso_storage, &
      mt_leaf_water_moles, mt_leaf_delta_e, mt_leaf_delta_b, mt_leaf_peclet

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: mt_glacier_patch
   PUBLIC :: mt_waterbody_patch

   logical, PUBLIC :: mt_waterbody_hist_sample = .true.

CONTAINS

   SUBROUTINE mt_glacier_patch(ipatch, maxsnl, nl_soil, deltim, &
      prc_rain, prl_rain, prc_snow, prl_snow, rnof, qseva, qsubl, qsdew, qfros, &
      endwb, totwb, glacier_overflow_mass, errorw, wdsrf, scv, t_grnd, forc_q, &
      forc_psrf, wliq_soisno, wice_soisno)

      IMPLICIT NONE
      integer,  intent(in) :: ipatch, maxsnl, nl_soil
      real(r8), intent(in) :: deltim
      real(r8), intent(in) :: prc_rain, prl_rain, prc_snow, prl_snow
      real(r8), intent(in) :: rnof, qseva, qsubl, qsdew, qfros
      real(r8), intent(in) :: endwb, totwb, glacier_overflow_mass, errorw
      real(r8), intent(in) :: wdsrf, scv, t_grnd, forc_q, forc_psrf
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil)

      integer  :: imt, j_mt, snl_mt
      real(r8) :: R_init, R_precip, R_rain, R_vapor, R_final
      real(r8) :: R_dew, R_frost, R_runoff
      real(r8) :: precip_mass, rain_mass, rnof_mass
      real(r8) :: evap_mass, dep_mass
      real(r8) :: evap_liq_mass, evap_ice_mass, dep_liq_mass, dep_ice_mass
      real(r8) :: water_dS, water_end, water_beg, water_input
      real(r8) :: water_before_output, water_after_evap, water_after_liq
      real(r8) :: mt_input, mt_evap, mt_evap_liq, mt_evap_ice, mt_after_liq
      real(r8) :: mt_rnof, mt_available, mt_final
      real(r8) :: mt_held_storage, surface_residue_beg
      real(r8) :: surface_liquid_end, surface_carrier, surface_residue_export
      real(r8) :: xerr_material
      logical  :: mixed_signature, fixed_signature, frac_active, nonvolatile_solute

      IF (nmts <= 0) RETURN
      IF (mt_aquifer_ref_water(ipatch) > 0._r8) &
         CALL CoLM_stop('glacier transition with isotope aquifer reference needs explicit water transfer')
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         ! Glacier physics exposes a single bulk snow carrier (scv).  Fold any
         ! layered state left by the preceding patch/step into that carrier
         ! before the balance snapshot instead of clearing it outside the
         ! accounting window.
         mt_scv(imt, ipatch) = mt_scv(imt, ipatch) + &
            sum(mt_wliq_soisno(imt, maxsnl+1:0, ipatch)) + &
            sum(mt_wice_soisno(imt, maxsnl+1:0, ipatch))
         mt_wliq_soisno(imt, maxsnl+1:0, ipatch) = 0._r8
         mt_wice_soisno(imt, maxsnl+1:0, ipatch) = 0._r8
         ! These carrier pools do not exist on glacier patches.  Quarantine
         ! conservative mass left by a patch-class transition before clearing
         ! the incompatible dissolved states; it will re-enter only through a
         ! resolved glacier surface carrier below.
         IF (material_is_nonvolatile_solute(imt)) THEN
            CALL move_special_subsurface_residue_to_surface(imt, ipatch)
            IF (mt_has_dissolved_limit(imt)) THEN
               mt_surface_solid(imt, ipatch) = mt_surface_solid(imt, ipatch) + &
                  mt_surface_residue(imt, ipatch) + mt_ldew_rain(imt, ipatch) + &
                  mt_ldew_snow(imt, ipatch) + mt_wetwat(imt, ipatch) + &
                  mt_canopy_solid(imt, ipatch)
               mt_surface_residue(imt, ipatch) = 0._r8
               mt_surface_solid(imt, ipatch) = mt_surface_solid(imt, ipatch) + &
                  sum(mt_solid_soisno(imt, maxsnl+1:0, ipatch))
               mt_solid_soisno(imt, maxsnl+1:0, ipatch) = 0._r8
               IF (allocated(mt_waterstorage)) THEN
                  mt_surface_solid(imt, ipatch) = mt_surface_solid(imt, ipatch) + &
                     mt_waterstorage(imt, ipatch) + mt_waterstorage_solid(imt, ipatch)
               ENDIF
            ELSE
               mt_surface_residue(imt, ipatch) = mt_surface_residue(imt, ipatch) + &
                  mt_ldew_rain(imt, ipatch) + mt_ldew_snow(imt, ipatch) + &
                  mt_wetwat(imt, ipatch)
               IF (allocated(mt_waterstorage)) THEN
                  mt_surface_residue(imt, ipatch) = mt_surface_residue(imt, ipatch) + &
                     mt_waterstorage(imt, ipatch)
               ENDIF
            ENDIF
         ENDIF
         mt_ldew_rain(imt, ipatch) = 0._r8
         mt_ldew_snow(imt, ipatch) = 0._r8
         mt_wetwat   (imt, ipatch) = 0._r8
         mt_canopy_solid(imt, ipatch) = 0._r8
         IF (allocated(mt_waterstorage)) mt_waterstorage(imt, ipatch) = 0._r8
         mt_waterstorage_solid(imt, ipatch) = 0._r8
      ENDDO

      snl_mt = 0
      CALL mt_save_storage(ipatch, snl_mt, nl_soil)
      precip_mass = (prc_rain + prl_rain + prc_snow + prl_snow) * deltim
      rain_mass   = (prc_rain + prl_rain) * deltim
      rnof_mass   = max(rnof, 0._r8) * deltim
      evap_liq_mass = max(qseva, 0._r8) * deltim
      evap_ice_mass = max(qsubl, 0._r8) * deltim
      dep_liq_mass  = max(qsdew, 0._r8) * deltim
      dep_ice_mass  = max(qfros, 0._r8) * deltim
      evap_mass = evap_liq_mass + evap_ice_mass
      dep_mass  = dep_liq_mass  + dep_ice_mass
      water_input = precip_mass + dep_mass
      water_dS = endwb - totwb - glacier_overflow_mass
      water_end = max(wdsrf, 0._r8) + max(scv, 0._r8)
      DO j_mt = 1, nl_soil
         water_end = water_end + max(wliq_soisno(j_mt), 0._r8) + max(wice_soisno(j_mt), 0._r8)
      ENDDO
      water_beg = water_end - water_dS

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         ! No per-substep reset: the lake/waterbody substeps all feed the
         ! same main-step river hand-off (mt_add_rnof_step below).
         R_init = material_init_water_ratio(imt)
         frac_active = material_fractionation_active(imt)
         nonvolatile_solute = material_is_nonvolatile_solute(imt)
         surface_residue_beg = mt_surface_residue(imt, ipatch)
         mt_held_storage = mt_subsurface_residue(imt, ipatch) + surface_residue_beg + &
            solid_inventory(imt, ipatch, 1, nl_soil)
         fixed_signature = material_can_use_fixed_signature(imt) .and. .not. frac_active
         IF (allocated(mt_runtime_forced)) THEN
            fixed_signature = fixed_signature .and. .not. mt_runtime_forced(imt)
         ENDIF
         ! The beginning inventory already includes the signed NSS anomaly.
         ! It must enter the receiving water once, never a fixed-signature reset.
         IF (material_is_isotope(imt)) THEN
            IF (mt_leaf_iso_storage(imt, ipatch) /= 0._r8) THEN
               IF (water_end <= mt_water_min_for_ratio .or. &
                   mt_storage_beg(imt, ipatch) - mt_held_storage < 0._r8) &
                  CALL CoLM_stop('special patch leaf isotope anomaly has no safe receiving water inventory')
               IF (fixed_signature) &
                  CALL CoLM_stop('special patch leaf isotope anomaly requires mixed-signature transport')
            ENDIF
         ENDIF
         mixed_signature = .not. fixed_signature

         IF (mixed_signature) THEN
            ! Snow ratio for all precipitation plus the rain correction: the
            ! correction is exactly zero where rain and snow share a ratio
            ! (see material_forcing_rain_value in MOD_MT_Forcing).
            R_precip = material_forcing_snow_value(imt, ipatch)
            R_rain   = material_forcing_rain_value(imt, ipatch)
            R_vapor  = material_forcing_vapor_value (imt, ipatch)
            IF (nonvolatile_solute) THEN
               R_dew = 0._r8
               R_frost = 0._r8
            ELSE
               R_dew = R_vapor
               R_frost = R_vapor
               IF (frac_active) THEN
                  R_dew = material_equilibrium_deposition_ratio(imt, R_vapor, t_grnd, .false.)
                  R_frost = material_equilibrium_deposition_ratio(imt, R_vapor, t_grnd, .true.)
               ENDIF
            ENDIF
            mt_input = precip_mass * R_precip + rain_mass * (R_rain - R_precip) &
               + dep_liq_mass * R_dew + dep_ice_mass * R_frost
            mt_available = max(mt_storage_beg(imt, ipatch) - &
               mt_held_storage + mt_input, 0._r8)
            water_before_output = water_beg + water_input
            IF (nonvolatile_solute) THEN
               mt_evap_liq = 0._r8
               mt_evap_ice = 0._r8
            ELSE
               mt_evap_liq = material_atmospheric_mt_loss(mt_available, water_before_output, &
                  evap_liq_mass, t_grnd, .false., glacier_evap_affine_for, mt_tiny, &
                  material_ratio_cap(imt), nonvolatile_solute)
               mt_after_liq = mt_available - mt_evap_liq
               water_after_liq = water_before_output - evap_liq_mass
               mt_evap_ice = material_atmospheric_mt_loss(mt_after_liq, water_after_liq, &
                  evap_ice_mass, t_grnd, .true., glacier_evap_affine_for, mt_tiny, &
                  material_ratio_cap(imt), nonvolatile_solute, skin_mass = DEF_ISO_SUBL_SKIN_MM)
            ENDIF
            mt_evap = mt_evap_liq + mt_evap_ice
            water_after_evap = water_before_output - evap_mass
            mt_final = max(mt_available - mt_evap, 0._r8)
            CALL mt_equilibrate_dissolved(imt, water_after_evap, mt_final, &
               mt_surface_solid(imt, ipatch))
            IF (water_after_evap > mt_water_min_for_ratio) THEN
               R_runoff = mt_final / water_after_evap
            ELSE
               R_runoff = 0._r8
            ENDIF
            mt_rnof = min(rnof_mass * R_runoff, mt_final)
            surface_liquid_end = max(wdsrf, 0._r8) + max(wliq_soisno(1), 0._r8)
            surface_carrier = surface_liquid_end + rnof_mass
            surface_residue_export = 0._r8
            IF (nonvolatile_solute .and. .not. mt_has_dissolved_limit(imt) .and. &
                surface_residue_beg > mt_tiny .and. &
                rnof_mass > mt_tiny .and. &
                surface_carrier > mt_water_min_for_ratio) THEN
               surface_residue_export = min(surface_residue_beg * &
                  rnof_mass / surface_carrier, surface_residue_beg)
               surface_residue_beg = surface_residue_beg - surface_residue_export
               mt_rnof = mt_rnof + surface_residue_export
            ENDIF
            ! surface_residue_export comes from held mass, not mt_available.
            mt_final = max(mt_final - mt_rnof + surface_residue_export, 0._r8)
            CALL mt_equilibrate_dissolved(imt, water_end, mt_final, &
               mt_surface_solid(imt, ipatch))
            IF (water_end > mt_water_min_for_ratio) THEN
               R_final = mt_final / water_end
               IF (nonvolatile_solute) mt_surface_residue(imt, ipatch) = surface_residue_beg
            ELSE
               R_final = 0._r8
               IF (nonvolatile_solute) THEN
                  mt_surface_residue(imt, ipatch) = surface_residue_beg + mt_final
               ENDIF
            ENDIF
         ELSE
            mt_input = water_input * R_init
            mt_evap_liq = evap_liq_mass * R_init
            mt_evap_ice = evap_ice_mass * R_init
            mt_evap = mt_evap_liq + mt_evap_ice
            mt_rnof  = rnof_mass * R_init
            R_final   = R_init
         ENDIF

         IF (mt_input > 0._r8) a_mt_precip(imt, ipatch) = a_mt_precip(imt, ipatch) + mt_input
         CALL material_book_evap_loss(imt, ipatch, mt_evap_liq, evap_liq_mass, MT_EVAP_KIND_SOILEVAP)
         CALL material_book_evap_loss(imt, ipatch, mt_evap_ice, evap_ice_mass, MT_EVAP_KIND_SUBL)
         IF (mt_rnof > 0._r8) THEN
            CALL mt_add_rnof_step (imt, ipatch, mt_rnof)
            a_mt_rsur(imt, ipatch) = a_mt_rsur(imt, ipatch) + mt_rnof
            a_mt_rnof(imt, ipatch) = a_mt_rnof(imt, ipatch) + mt_rnof
         ENDIF
         CALL sync_mt_patch_ratio(imt, ipatch, snl_mt, maxsnl, nl_soil, &
            wliq_soisno, wice_soisno, 0._r8, wdsrf, scv, R_final)
         IF (material_is_isotope(imt)) THEN
            ! sync has transferred the anomaly included in mt_storage_beg.
            mt_leaf_iso_storage(imt, ipatch) = 0._r8
            mt_leaf_water_moles(imt, ipatch) = 0._r8
            mt_leaf_delta_e(imt, ipatch) = materials(imt)%init_delta
            mt_leaf_delta_b(imt, ipatch) = materials(imt)%init_delta
            mt_leaf_peclet(imt, ipatch) = 1._r8
         ENDIF
         IF (nonvolatile_solute .and. .not. mt_has_dissolved_limit(imt) .and. &
             mt_surface_residue(imt, ipatch) > mt_tiny) THEN
            IF (wdsrf > mt_water_min_for_ratio) THEN
               mt_wdsrf(imt, ipatch) = mt_wdsrf(imt, ipatch) + &
                  mt_surface_residue(imt, ipatch)
               mt_surface_residue(imt, ipatch) = 0._r8
            ELSEIF (wliq_soisno(1) > mt_water_min_for_ratio) THEN
               mt_wliq_soisno(imt, 1, ipatch) = &
                  mt_wliq_soisno(imt, 1, ipatch) + &
                  mt_surface_residue(imt, ipatch)
               mt_surface_residue(imt, ipatch) = 0._r8
            ENDIF
         ENDIF
      ENDDO

      CALL mt_apply_decay(ipatch, snl_mt, nl_soil, deltim)
      CALL mt_balance_check(ipatch, snl_mt, nl_soil, deltim, xerr_material, &
         patchtype_in = 3, water_err_in = errorw, water_dS_in = water_dS, &
         water_input_in = precip_mass + dep_mass, water_output_in = evap_mass + rnof_mass, &
         water_evap_in = evap_mass, water_rnof_in = rnof_mass)
      CALL mt_hist_accumulate(ipatch, snl_mt, maxsnl, nl_soil, 0._r8, 0._r8, &
         wliq_soisno(snl_mt+1:nl_soil), wice_soisno(snl_mt+1:nl_soil), 0._r8, wdsrf, 0._r8, scv)

   CONTAINS

      SUBROUTINE glacier_evap_affine_for (temp_k, from_ice, slope, offset)
         real(r8), intent(in)  :: temp_k
         logical,  intent(in) :: from_ice
         real(r8), intent(out) :: slope, offset
         real(r8) :: relhum, alpha_k, alpha_k_liq, alpha_k_ice

         slope  = 1._r8
         offset = 0._r8
         IF (.not. material_fractionation_active(imt)) RETURN
         alpha_k_liq = material_alpha_kinetic_craig_gordon(imt, .false.)
         alpha_k_ice = material_alpha_kinetic_craig_gordon(imt, .true.)
         alpha_k = merge(alpha_k_ice, alpha_k_liq, from_ice)
         relhum = material_surface_relhum(forc_q, forc_psrf, temp_k, from_ice)
         CALL material_craig_gordon_evap_affine(imt, R_vapor, temp_k, relhum, &
            alpha_k, from_ice, slope, offset)
      END SUBROUTINE glacier_evap_affine_for

   END SUBROUTINE mt_glacier_patch

   SUBROUTINE mt_waterbody_patch(ipatch, maxsnl, nl_soil, snl, deltim, &
      forc_rain, forc_snow, lake_deficit, rnof, qseva, qsubl, qsdew, qfros, &
      endwb, totwb, errorw, wa, wdsrf, scv, t_grnd, forc_q, forc_psrf, &
      forc_us, forc_vs, &
      wliq_soisno, wice_soisno, use_dynamic_lake)

      IMPLICIT NONE
      integer,  intent(in) :: ipatch, maxsnl, nl_soil, snl
      real(r8), intent(in) :: deltim
      real(r8), intent(in) :: forc_rain, forc_snow, lake_deficit, rnof
      real(r8), intent(in) :: qseva, qsubl, qsdew, qfros
      real(r8), intent(in) :: endwb, totwb, errorw, wa, wdsrf, scv, t_grnd, forc_q, forc_psrf
      ! Wind is needed because open-water kinetic fractionation is
      ! wind-dependent (Merlivat & Jouzel 1979).
      real(r8), intent(in) :: forc_us, forc_vs
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil)
      logical,  intent(in) :: use_dynamic_lake

      integer  :: imt, j_mt
      real(r8) :: R_init, R_precip, R_rain, R_vapor, R_pool, R_final
      real(r8) :: R_dew, R_frost, R_runoff
      real(r8) :: atm_precip_mass, deficit_mass, precip_mass, rnof_mass
      real(r8) :: evap_mass, dep_mass
      real(r8) :: evap_liq_mass, evap_ice_mass, dep_liq_mass, dep_ice_mass
      real(r8) :: water_dS, water_end, water_beg, water_input
      real(r8) :: water_before_output, water_after_evap, water_after_liq
      real(r8) :: mt_input, mt_evap, mt_evap_liq, mt_evap_ice, mt_after_liq
      real(r8) :: mt_rnof, mt_available, mt_final
      real(r8) :: mt_held_storage, surface_residue_beg
      real(r8) :: surface_liquid_end, surface_carrier, surface_residue_export
      real(r8) :: xerr_material
      logical  :: mixed_signature, fixed_signature, frac_active, nonvolatile_solute

      IF (nmts <= 0) RETURN
      IF (mt_aquifer_ref_water(ipatch) > 0._r8) &
         CALL CoLM_stop('waterbody transition with isotope aquifer reference needs explicit water transfer')
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         ! Waterbody patches have no canopy, wetland, or irrigation-reservoir
         ! carrier in this branch. Preserve conservative transition mass in a
         ! surface quarantine instead of reporting a carrierless concentration.
         IF (material_is_nonvolatile_solute(imt)) THEN
            CALL move_special_subsurface_residue_to_surface(imt, ipatch)
            IF (mt_has_dissolved_limit(imt)) THEN
               mt_surface_solid(imt, ipatch) = mt_surface_solid(imt, ipatch) + &
                  mt_surface_residue(imt, ipatch) + mt_ldew_rain(imt, ipatch) + &
                  mt_ldew_snow(imt, ipatch) + mt_wetwat(imt, ipatch) + &
                  mt_canopy_solid(imt, ipatch)
               mt_surface_residue(imt, ipatch) = 0._r8
               mt_surface_solid(imt, ipatch) = mt_surface_solid(imt, ipatch) + &
                  sum(mt_solid_soisno(imt, maxsnl+1:0, ipatch))
               mt_solid_soisno(imt, maxsnl+1:0, ipatch) = 0._r8
               IF (allocated(mt_waterstorage)) THEN
                  mt_surface_solid(imt, ipatch) = mt_surface_solid(imt, ipatch) + &
                     mt_waterstorage(imt, ipatch) + mt_waterstorage_solid(imt, ipatch)
               ENDIF
            ELSE
               mt_surface_residue(imt, ipatch) = mt_surface_residue(imt, ipatch) + &
                  mt_ldew_rain(imt, ipatch) + mt_ldew_snow(imt, ipatch) + &
                  mt_wetwat(imt, ipatch)
               IF (allocated(mt_waterstorage)) THEN
                  mt_surface_residue(imt, ipatch) = mt_surface_residue(imt, ipatch) + &
                     mt_waterstorage(imt, ipatch)
               ENDIF
            ENDIF
         ENDIF
         mt_ldew_rain(imt, ipatch) = 0._r8
         mt_ldew_snow(imt, ipatch) = 0._r8
         mt_wetwat   (imt, ipatch) = 0._r8
         mt_canopy_solid(imt, ipatch) = 0._r8
         IF (allocated(mt_waterstorage)) mt_waterstorage(imt, ipatch) = 0._r8
         mt_waterstorage_solid(imt, ipatch) = 0._r8
      ENDDO

      CALL mt_save_storage(ipatch, maxsnl, nl_soil)
      atm_precip_mass = (forc_rain + forc_snow) * deltim
      IF (use_dynamic_lake) THEN
         precip_mass = atm_precip_mass
      ELSE
         precip_mass = (forc_rain + forc_snow + lake_deficit) * deltim
      ENDIF
      deficit_mass = precip_mass - atm_precip_mass
      rnof_mass   = max(rnof, 0._r8) * deltim
      evap_liq_mass = max(qseva, 0._r8) * deltim
      evap_ice_mass = max(qsubl, 0._r8) * deltim
      dep_liq_mass  = max(qsdew, 0._r8) * deltim
      dep_ice_mass  = max(qfros, 0._r8) * deltim
      evap_mass = evap_liq_mass + evap_ice_mass
      dep_mass  = dep_liq_mass  + dep_ice_mass
      water_input = precip_mass + dep_mass
      water_dS = endwb - totwb
      water_end = wa + max(wdsrf, 0._r8)
      DO j_mt = maxsnl + 1, nl_soil
         IF (j_mt >= snl + 1) THEN
            water_end = water_end + max(wliq_soisno(j_mt), 0._r8) + max(wice_soisno(j_mt), 0._r8)
         ENDIF
      ENDDO
      IF (snl >= 0) water_end = water_end + max(scv, 0._r8)
      water_beg = water_end - water_dS

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         ! No per-substep reset: the lake/waterbody substeps all feed the
         ! same main-step river hand-off (mt_add_rnof_step below).
         R_init = material_init_water_ratio(imt)
         frac_active = material_fractionation_active(imt)
         nonvolatile_solute = material_is_nonvolatile_solute(imt)
         surface_residue_beg = mt_surface_residue(imt, ipatch)
         mt_held_storage = mt_subsurface_residue(imt, ipatch) + surface_residue_beg + &
            solid_inventory(imt, ipatch, maxsnl + 1, nl_soil)
         fixed_signature = material_can_use_fixed_signature(imt) .and. .not. frac_active
         IF (allocated(mt_runtime_forced)) THEN
            fixed_signature = fixed_signature .and. .not. mt_runtime_forced(imt)
         ENDIF
         ! The beginning inventory already includes the signed NSS anomaly.
         ! It must enter the receiving water once, never a fixed-signature reset.
         IF (material_is_isotope(imt)) THEN
            IF (mt_leaf_iso_storage(imt, ipatch) /= 0._r8) THEN
               IF (water_end <= mt_water_min_for_ratio .or. &
                   mt_storage_beg(imt, ipatch) - mt_held_storage < 0._r8) &
                  CALL CoLM_stop('special patch leaf isotope anomaly has no safe receiving water inventory')
               IF (fixed_signature) &
                  CALL CoLM_stop('special patch leaf isotope anomaly requires mixed-signature transport')
            ENDIF
         ENDIF
         mixed_signature = .not. fixed_signature

         IF (mixed_signature) THEN
            ! As on glaciers: snow ratio for all precipitation plus the rain
            ! correction, exactly zero where the two ratios are the same.
            R_precip = material_forcing_snow_value(imt, ipatch)
            R_rain   = material_forcing_rain_value(imt, ipatch)
            R_vapor  = material_forcing_vapor_value (imt, ipatch)
            IF (water_beg > mt_water_min_for_ratio) THEN
               R_pool = max(mt_storage_beg(imt, ipatch) - &
                  mt_held_storage, 0._r8) / water_beg
            ELSE
               R_pool = R_rain
            ENDIF
            IF (nonvolatile_solute) THEN
               R_dew = 0._r8
               R_frost = 0._r8
            ELSE
               R_dew = R_vapor
               R_frost = R_vapor
               IF (frac_active) THEN
                  R_dew = material_equilibrium_deposition_ratio(imt, R_vapor, t_grnd, .false.)
                  R_frost = material_equilibrium_deposition_ratio(imt, R_vapor, t_grnd, .true.)
               ENDIF
            ENDIF
            mt_input = atm_precip_mass * R_precip + forc_rain * deltim * (R_rain - R_precip) &
               + dep_liq_mass * R_dew + dep_ice_mass * R_frost + deficit_mass * R_pool
            mt_available = max(mt_storage_beg(imt, ipatch) - &
               mt_held_storage + mt_input, 0._r8)
            water_before_output = water_beg + water_input
            IF (nonvolatile_solute) THEN
               mt_evap_liq = 0._r8
               mt_evap_ice = 0._r8
            ELSE
               mt_evap_liq = material_atmospheric_mt_loss(mt_available, water_before_output, &
                  evap_liq_mass, t_grnd, .false., waterbody_evap_affine_for, mt_tiny, &
                  material_ratio_cap(imt), nonvolatile_solute)
               mt_after_liq = mt_available - mt_evap_liq
               water_after_liq = water_before_output - evap_liq_mass
               mt_evap_ice = material_atmospheric_mt_loss(mt_after_liq, water_after_liq, &
                  evap_ice_mass, t_grnd, .true., waterbody_evap_affine_for, mt_tiny, &
                  material_ratio_cap(imt), nonvolatile_solute, skin_mass = DEF_ISO_SUBL_SKIN_MM)
            ENDIF
            mt_evap = mt_evap_liq + mt_evap_ice
            water_after_evap = water_before_output - evap_mass
            mt_final = max(mt_available - mt_evap, 0._r8)
            CALL mt_equilibrate_dissolved(imt, water_after_evap, mt_final, &
               mt_surface_solid(imt, ipatch))
            IF (water_after_evap > mt_water_min_for_ratio) THEN
               R_runoff = mt_final / water_after_evap
            ELSE
               R_runoff = 0._r8
            ENDIF
            mt_rnof = min(rnof_mass * R_runoff, mt_final)
            surface_liquid_end = max(wdsrf, 0._r8) + &
               max(wliq_soisno(snl+1), 0._r8)
            surface_carrier = surface_liquid_end + rnof_mass
            surface_residue_export = 0._r8
            IF (nonvolatile_solute .and. .not. mt_has_dissolved_limit(imt) .and. &
                surface_residue_beg > mt_tiny .and. &
                rnof_mass > mt_tiny .and. &
                surface_carrier > mt_water_min_for_ratio) THEN
               surface_residue_export = min(surface_residue_beg * &
                  rnof_mass / surface_carrier, surface_residue_beg)
               surface_residue_beg = surface_residue_beg - surface_residue_export
               mt_rnof = mt_rnof + surface_residue_export
            ENDIF
            ! surface_residue_export comes from held mass, not mt_available.
            mt_final = max(mt_final - mt_rnof + surface_residue_export, 0._r8)
            CALL mt_equilibrate_dissolved(imt, water_end, mt_final, &
               mt_surface_solid(imt, ipatch))
            IF (water_end > mt_water_min_for_ratio) THEN
               R_final = mt_final / water_end
               IF (nonvolatile_solute) mt_surface_residue(imt, ipatch) = surface_residue_beg
            ELSE
               R_final = 0._r8
               IF (nonvolatile_solute) THEN
                  mt_surface_residue(imt, ipatch) = surface_residue_beg + mt_final
               ENDIF
            ENDIF
         ELSE
            mt_input = water_input * R_init
            mt_evap_liq = evap_liq_mass * R_init
            mt_evap_ice = evap_ice_mass * R_init
            mt_evap = mt_evap_liq + mt_evap_ice
            mt_rnof  = rnof_mass * R_init
            R_final   = R_init
         ENDIF

         IF (mt_input > 0._r8) a_mt_precip(imt, ipatch) = a_mt_precip(imt, ipatch) + mt_input
         CALL material_book_evap_loss(imt, ipatch, mt_evap_liq, evap_liq_mass, MT_EVAP_KIND_SOILEVAP)
         CALL material_book_evap_loss(imt, ipatch, mt_evap_ice, evap_ice_mass, MT_EVAP_KIND_SUBL)
         IF (mt_rnof > 0._r8) THEN
            CALL mt_add_rnof_step (imt, ipatch, mt_rnof)
            a_mt_rsur(imt, ipatch) = a_mt_rsur(imt, ipatch) + mt_rnof
            a_mt_rnof(imt, ipatch) = a_mt_rnof(imt, ipatch) + mt_rnof
         ENDIF
         CALL sync_mt_patch_ratio(imt, ipatch, snl, maxsnl, nl_soil, &
            wliq_soisno, wice_soisno, wa, wdsrf, scv, R_final)
         IF (material_is_isotope(imt)) THEN
            ! sync has transferred the anomaly included in mt_storage_beg.
            mt_leaf_iso_storage(imt, ipatch) = 0._r8
            mt_leaf_water_moles(imt, ipatch) = 0._r8
            mt_leaf_delta_e(imt, ipatch) = materials(imt)%init_delta
            mt_leaf_delta_b(imt, ipatch) = materials(imt)%init_delta
            mt_leaf_peclet(imt, ipatch) = 1._r8
         ENDIF
         IF (nonvolatile_solute .and. .not. mt_has_dissolved_limit(imt) .and. &
             mt_surface_residue(imt, ipatch) > mt_tiny) THEN
            IF (wdsrf > mt_water_min_for_ratio) THEN
               mt_wdsrf(imt, ipatch) = mt_wdsrf(imt, ipatch) + &
                  mt_surface_residue(imt, ipatch)
               mt_surface_residue(imt, ipatch) = 0._r8
            ELSEIF (wliq_soisno(snl+1) > mt_water_min_for_ratio) THEN
               mt_wliq_soisno(imt, snl+1, ipatch) = &
                  mt_wliq_soisno(imt, snl+1, ipatch) + &
                  mt_surface_residue(imt, ipatch)
               mt_surface_residue(imt, ipatch) = 0._r8
            ENDIF
         ENDIF
      ENDDO

      CALL mt_apply_decay(ipatch, maxsnl, nl_soil, deltim)
      CALL mt_balance_check(ipatch, maxsnl, nl_soil, deltim, xerr_material, &
         patchtype_in = 4, water_err_in = errorw, water_dS_in = water_dS, &
         water_input_in = precip_mass + dep_mass, water_output_in = evap_mass + rnof_mass, &
         water_evap_in = evap_mass, water_rnof_in = rnof_mass)
      IF (mt_waterbody_hist_sample) CALL mt_hist_accumulate(ipatch, snl, maxsnl, nl_soil, 0._r8, 0._r8, &
         wliq_soisno(snl+1:nl_soil), wice_soisno(snl+1:nl_soil), wa, wdsrf, 0._r8, scv)

   CONTAINS

      SUBROUTINE waterbody_evap_affine_for (temp_k, from_ice, slope, offset)
         real(r8), intent(in)  :: temp_k
         logical,  intent(in) :: from_ice
         real(r8), intent(out) :: slope, offset
         real(r8) :: relhum, alpha_k, alpha_k_liq, alpha_k_ice

         slope  = 1._r8
         offset = 0._r8
         IF (.not. material_fractionation_active(imt)) RETURN
         alpha_k_liq = material_alpha_kinetic_open_water(imt, &
            sqrt(max(forc_us*forc_us + forc_vs*forc_vs, 0._r8)))
         alpha_k_ice = material_alpha_kinetic_craig_gordon(imt, .true.)
         alpha_k = merge(alpha_k_ice, alpha_k_liq, from_ice)
         relhum = material_surface_relhum(forc_q, forc_psrf, temp_k, from_ice)
         CALL material_craig_gordon_evap_affine(imt, R_vapor, temp_k, relhum, &
            alpha_k, from_ice, slope, offset)
      END SUBROUTINE waterbody_evap_affine_for

   END SUBROUTINE mt_waterbody_patch

   SUBROUTINE move_special_subsurface_residue_to_surface (imt, ipatch)
      integer, intent(in) :: imt, ipatch

      IF (mt_subsurface_residue(imt, ipatch) <= mt_tiny) RETURN
      IF (mt_has_dissolved_limit(imt)) THEN
         mt_surface_solid(imt, ipatch) = mt_surface_solid(imt, ipatch) + &
            mt_subsurface_residue(imt, ipatch)
      ELSE
         mt_surface_residue(imt, ipatch) = mt_surface_residue(imt, ipatch) + &
            mt_subsurface_residue(imt, ipatch)
      ENDIF
      mt_subsurface_residue(imt, ipatch) = 0._r8
   END SUBROUTINE move_special_subsurface_residue_to_surface

   real(r8) FUNCTION solid_inventory (imt, ipatch, lb, nl_soil)
      integer, intent(in) :: imt, ipatch, lb, nl_soil

      solid_inventory = mt_canopy_solid(imt, ipatch) + &
         mt_surface_solid(imt, ipatch) + mt_subsurface_solid(imt, ipatch) + &
         mt_waterstorage_solid(imt, ipatch)
      IF (lb <= nl_soil) solid_inventory = solid_inventory + &
         sum(mt_solid_soisno(imt, lb:nl_soil, ipatch))
   END FUNCTION solid_inventory


   real(r8) FUNCTION material_ratio_cap (imt)
      integer, intent(in) :: imt

      material_ratio_cap = merge(materials(imt)%ref_ratio * (1._r8 + mt_delta_sanity_max / 1000._r8), &
         0._r8, material_fractionation_active(imt))
   END FUNCTION material_ratio_cap

END MODULE MOD_MT_SpecialPatches
#endif
