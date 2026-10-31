#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_Evapo

   USE MOD_Precision
   USE MOD_MT_Defs, only: nmts, mt_tiny, mt_uses_generic_water, material_init_water_ratio, &
      material_is_nonvolatile_solute, materials, mt_delta_sanity_max, mt_water_min_for_ratio
   USE MOD_MT_Solute, only: mt_equilibrate_dissolved
   USE MOD_MT_Forcing, only: material_forcing_vapor_value
   USE MOD_MT_Frac, only: material_fractionation_active, material_alpha_kinetic_craig_gordon, &
      material_craig_gordon_evap_affine, material_equilibrium_deposition_ratio, &
      material_rayleigh_freezing_loss, material_surface_relhum, &
      material_equilibration_exchange, material_alpha_liq_vap, material_alpha_kinetic_open_water
   USE MOD_MT_Namelist, only: DEF_ISO_SUBL_SKIN_MM, DEF_ISO_CANOPY_EQUILIBRATION
   USE MOD_MT_EvapLimit, only: material_evaporative_mt_loss, material_atmospheric_mt_loss, &
      material_skin_limited_mt_loss
	   USE MOD_MT_Vars, only: mt_ldew_rain, mt_ldew_snow, &
	      mt_wliq_soisno, mt_wice_soisno, mt_solid_soisno, mt_canopy_solid, &
	      mt_numerical_residual_step, &
	      a_mt_precip, a_water_precip, a_mt_vapor_exchange, material_book_evap_loss, &
	      MT_EVAP_KIND_CANOPYEVAP, MT_EVAP_KIND_SOILEVAP, MT_EVAP_KIND_SUBL

   IMPLICIT NONE

CONTAINS

   SUBROUTINE mt_flood_evap_loss(ipatch, water_credit, material_credit, water_evap, &
      temp_k, forc_q, forc_psrf, forc_us, forc_vs, material_loss)
      ! Flood evaporation belongs to the routing pool, not land a_mt_evap.
      ! A negative Craig-Gordon loss is atmospheric isotope uptake.
      integer, intent(in) :: ipatch
      real(r8), intent(in) :: water_credit, material_credit(:), water_evap
      real(r8), intent(in) :: temp_k, forc_q, forc_psrf, forc_us, forc_vs
      real(r8), intent(out) :: material_loss(:)
      integer :: imt
      real(r8) :: vapor_ratio

      material_loss = 0._r8
      IF (water_credit <= mt_tiny .or. water_evap <= mt_tiny) RETURN
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         IF (material_is_nonvolatile_solute(imt)) CYCLE
         vapor_ratio = material_forcing_vapor_value(imt, ipatch)
         material_loss(imt) = material_atmospheric_mt_loss( &
            material_credit(imt), water_credit, water_evap, temp_k, .false., &
            evap_affine, mt_tiny, &
            merge(materials(imt)%ref_ratio * (1._r8+mt_delta_sanity_max/1000._r8), &
                  0._r8, material_fractionation_active(imt)), .false.)
      ENDDO

   CONTAINS
      SUBROUTINE evap_affine(temp, from_ice, slope, offset)
         real(r8), intent(in) :: temp
         logical, intent(in) :: from_ice
         real(r8), intent(out) :: slope, offset
         real(r8) :: relhum, alpha_k
         slope  = 1._r8
         offset = 0._r8
         IF (.not. material_fractionation_active(imt)) RETURN
         relhum = material_surface_relhum(forc_q, forc_psrf, temp, from_ice)
         alpha_k = material_alpha_kinetic_open_water(imt, &
            sqrt(max(forc_us*forc_us+forc_vs*forc_vs,0._r8)))
         CALL material_craig_gordon_evap_affine(imt, vapor_ratio, temp, relhum, &
            alpha_k, from_ice, slope, offset)
      END SUBROUTINE evap_affine
   END SUBROUTINE mt_flood_evap_loss

   !---------------------------------------------------------------
   ! Delta-based ET material update after THERMAL.
   !
   ! Compare post-THERMAL vs pre-THERMAL water states:
   !   Decrease → evaporation/sublimation/transpiration (output)
   !   Increase → dew/frost (input) or thaw/freeze (internal)
   !---------------------------------------------------------------
   SUBROUTINE mt_evapo (ipatch, deltim, snl, nl_soil, &
      ldew_rain, ldew_snow, ldew_rain_bef, ldew_snow_bef, &
      wliq_soisno, wice_soisno, &
      wliq_soisno_bef, wice_soisno_bef, &
      canopy_smelt_mass_th, canopy_frzc_mass_th, &
      soil_thaw_mass_th, soil_frzc_mass_th, &
      tleaf_frac, t_soisno_frac, forc_q_frac, forc_psrf_frac)

      IMPLICIT NONE
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: deltim
      integer,  intent(in) :: snl, nl_soil
      real(r8), intent(in) :: ldew_rain, ldew_snow
      real(r8), intent(in) :: ldew_rain_bef, ldew_snow_bef
      real(r8), intent(in) :: wliq_soisno(snl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno(snl+1:nl_soil)
      real(r8), intent(in) :: wliq_soisno_bef(snl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno_bef(snl+1:nl_soil)
      ! Explicit canopy snow→rain melt and rain→snow freeze
      ! masses (mm of water-equivalent over deltim) reported by THERMAL /
      ! LeafTemperature. When present, they replace the d_rain/d_snow
      ! sign-pattern heuristic so a simultaneous (qmelt + qevpl) or
      ! (qfrz + qsubl) episode no longer collapses one branch into the
      ! other. The heuristic is preserved as a fallback for callers that
      ! do not (yet) plumb these values — important for hot-start /
      ! third-party drivers that bind mt_evapo directly.
      real(r8), intent(in), optional :: canopy_smelt_mass_th
      real(r8), intent(in), optional :: canopy_frzc_mass_th
      real(r8), intent(in), optional :: soil_thaw_mass_th(snl+1:nl_soil)
      real(r8), intent(in), optional :: soil_frzc_mass_th(snl+1:nl_soil)
      real(r8), intent(in), optional :: tleaf_frac
      real(r8), intent(in), optional :: t_soisno_frac(snl+1:nl_soil)
      real(r8), intent(in), optional :: forc_q_frac
      real(r8), intent(in), optional :: forc_psrf_frac

      integer  :: imt, j, lb
      real(r8) :: ratio, mt_flux, R_atm, R_vapor
	      real(r8) :: d_rain, d_snow, d_wliq, d_wice, water_loss, mt_resid
      real(r8) :: thaw_amt, freeze_amt
      ! Post-internal-transfer pool sizes used for soil evaporation/sublimation.
      real(r8) :: wliq_post_phase, wice_post_phase
      ! Canopy phase-change detection variables.
      ! Mirrors the soil-layer thaw/freeze handling at L94-119: separate
      ! INTERNAL phase transfer (snow↔rain) from EXTERNAL atm I/O
      ! (evap / dew / sublim / frost). Without this split, LeafTemperature's
      ! canopy melt (which moves mass from ldew_snow to ldew_rain) was
      ! misclassified as simultaneous "dew on canopy" + "sublimation from
      ! canopy snow", inflating both a_mt_precip and a_mt_evap and
      ! polluting mt_ldew_rain with R_atm signature instead of the
      ! original snow signature R_canopy_snow.
      real(r8) :: ldew_rain_pre_phase, ldew_snow_pre_phase
      real(r8) :: d_rain_external, d_snow_external
      real(r8) :: equil_gain
      IF (nmts <= 0) RETURN
      lb = snl + 1

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE

         ! N04: a nonvolatile solute is retained in mt_wice while its ice
         ! carrier sublimates (it must not be lost with water vapour).  After
         ! the carrier is gone that inventory is invisible to the liquid/solid
         ! equilibrium and never re-dissolves on rewetting.  Move it once to
         ! the same layer's rewettable inventory as soon as the ice carrier is
         ! empty, then equilibrate against the current liquid water.
         IF (material_is_nonvolatile_solute(imt)) THEN
            DO j = lb, nl_soil
               IF (mt_wice_soisno(imt, j, ipatch) > mt_tiny .and. &
                   max(wice_soisno(j), 0._r8) <= mt_water_min_for_ratio) THEN
                  IF (max(wliq_soisno(j), 0._r8) > mt_water_min_for_ratio) THEN
                     mt_wliq_soisno(imt, j, ipatch) = &
                        mt_wliq_soisno(imt, j, ipatch) + mt_wice_soisno(imt, j, ipatch)
                  ELSE
                     mt_solid_soisno(imt, j, ipatch) = &
                        mt_solid_soisno(imt, j, ipatch) + mt_wice_soisno(imt, j, ipatch)
                  ENDIF
                  mt_wice_soisno(imt, j, ipatch) = 0._r8
               ENDIF
               CALL mt_equilibrate_dissolved(imt, max(wliq_soisno(j), 0._r8), &
                  mt_wliq_soisno(imt, j, ipatch), mt_solid_soisno(imt, j, ipatch))
            ENDDO
         ENDIF

         R_atm = material_forcing_vapor_value(imt, ipatch)
         R_vapor = R_atm

         ! LeafTemperature applies atmospheric canopy exchange first and only
         ! then melts/freezes the remaining carrier.  Preserve that ordering:
         ! fractionating evaporation/sublimation does not commute with a later
         ! internal phase transfer.
         d_rain = ldew_rain - ldew_rain_bef
         d_snow = ldew_snow - ldew_snow_bef

         ! Prefer the explicit qmelt / qfrz masses reported
         ! by THERMAL when the caller plumbs them. The d_rain/d_snow
         ! sign-pattern heuristic (legacy fallback below) silently
         ! collapses thaw_amt to min(|d_snow|, d_rain), which under-counts
         ! when canopy melt and rain-pool evaporation fire in the same
         ! step (e.g. tl ~ tfrz with qmelt > 0 AND qevpl > 0 → d_rain =
         ! qmelt - qevpl < qmelt; min picks d_rain instead of qmelt).
         ! Phase-1 invariant hides this because all R = R_atm = R_init,
         ! but Phase 2 (time-varying R_atm) would mis-charge the missing
         ! qmelt fraction to sublimation/dew.
         thaw_amt   = 0._r8   ! snow→rain (canopy melt)
         freeze_amt = 0._r8   ! rain→snow (canopy refreeze)
         IF (present(canopy_smelt_mass_th) .or. present(canopy_frzc_mass_th)) THEN
            IF (present(canopy_smelt_mass_th)) THEN
               thaw_amt = max(canopy_smelt_mass_th, 0._r8)
            ENDIF
            IF (present(canopy_frzc_mass_th)) THEN
               freeze_amt = max(canopy_frzc_mass_th, 0._r8)
            ENDIF
         ELSE
            ! Legacy heuristic — used only when the caller has not
            ! yet been updated to forward THERMAL's explicit phase-change
            ! mass. Misattributes mass when melt+evap (or freeze+frost)
            ! coexist in one step; see comment above.
            IF (d_snow < -mt_tiny .and. d_rain > mt_tiny) THEN
               thaw_amt = min(abs(d_snow), d_rain)
            ENDIF
            IF (d_rain < -mt_tiny .and. d_snow > mt_tiny) THEN
               freeze_amt = min(abs(d_rain), d_snow)
            ENDIF
         ENDIF

         ! Remove the explicit internal phase contribution from the total delta.
         d_rain_external = d_rain - thaw_amt + freeze_amt
         d_snow_external = d_snow + thaw_amt - freeze_amt

         ! --- Canopy rain external exchange: evap/dew (host operation 1) ---
         IF (d_rain_external < -mt_tiny) THEN
            ! Net rain loss = canopy evaporation
            IF (ldew_rain_bef > mt_tiny) THEN
	               mt_flux = evaporative_mt_loss(mt_ldew_rain(imt, ipatch), &
	                  ldew_rain_bef, abs(d_rain_external), canopy_temp(), .false.)
	               mt_ldew_rain(imt, ipatch) = mt_ldew_rain(imt, ipatch) - mt_flux
	               CALL material_book_evap_loss(imt, ipatch, mt_flux, abs(d_rain_external), &
	                  MT_EVAP_KIND_CANOPYEVAP)
            ENDIF
         ELSEIF (d_rain_external > mt_tiny) THEN
            ! Net rain gain = dew deposition
            ratio = deposition_ratio_for(canopy_temp(), .false.)
            mt_flux = d_rain_external * ratio
            mt_ldew_rain(imt, ipatch) = mt_ldew_rain(imt, ipatch) + mt_flux
            a_mt_precip(imt, ipatch) = a_mt_precip(imt, ipatch) + mt_flux
            a_water_precip(imt, ipatch) = a_water_precip(imt, ipatch) + d_rain_external
         ENDIF
         mt_ldew_rain(imt, ipatch) = max(mt_ldew_rain(imt, ipatch), 0._r8)

         ! --- Canopy snow external exchange: sublim/frost (host operation 1) ---
         IF (d_snow_external < -mt_tiny) THEN
            ! Net snow loss = canopy sublimation
            IF (ldew_snow_bef > mt_tiny) THEN
	               mt_flux = evaporative_mt_loss(mt_ldew_snow(imt, ipatch), &
	                  ldew_snow_bef, abs(d_snow_external), canopy_temp(), .true.)
	               mt_ldew_snow(imt, ipatch) = mt_ldew_snow(imt, ipatch) - mt_flux
	               CALL material_book_evap_loss(imt, ipatch, mt_flux, abs(d_snow_external), &
	                  MT_EVAP_KIND_SUBL)
            ENDIF
         ELSEIF (d_snow_external > mt_tiny) THEN
            ! Net snow gain = frost deposition
            ratio = deposition_ratio_for(canopy_temp(), .true.)
            mt_flux = d_snow_external * ratio
            mt_ldew_snow(imt, ipatch) = mt_ldew_snow(imt, ipatch) + mt_flux
            a_mt_precip(imt, ipatch) = a_mt_precip(imt, ipatch) + mt_flux
            a_water_precip(imt, ipatch) = a_water_precip(imt, ipatch) + d_snow_external
         ENDIF
         mt_ldew_snow(imt, ipatch) = max(mt_ldew_snow(imt, ipatch), 0._r8)

         ! Phase-change carrier sizes after external exchange.  The explicit
         ! THERMAL masses are already limited by these host pools; retain a
         ! defensive clamp for roundoff and third-party direct callers.
         ldew_rain_pre_phase = max(ldew_rain_bef + d_rain_external, 0._r8)
         ldew_snow_pre_phase = max(ldew_snow_bef + d_snow_external, 0._r8)
         thaw_amt = min(thaw_amt, ldew_snow_pre_phase)
         freeze_amt = min(freeze_amt, ldew_rain_pre_phase + thaw_amt)

         ! --- Internal MELT: snow → rain (host operation 2) ---
         IF (thaw_amt > mt_tiny .and. ldew_snow_pre_phase > mt_tiny) THEN
            ratio = mt_ldew_snow(imt, ipatch) / ldew_snow_pre_phase
            mt_flux = min(thaw_amt * ratio, max(mt_ldew_snow(imt, ipatch), 0._r8))
            mt_ldew_snow(imt, ipatch) = mt_ldew_snow(imt, ipatch) - mt_flux
            mt_ldew_rain(imt, ipatch) = mt_ldew_rain(imt, ipatch) + mt_flux
         ENDIF

         ! --- Internal FREEZE: rain → snow (host operation 2) ---
         IF (freeze_amt > mt_tiny .and. ldew_rain_pre_phase + thaw_amt > mt_tiny) THEN
            mt_flux = material_rayleigh_freezing_loss(imt, mt_ldew_rain(imt, ipatch), &
               ldew_rain_pre_phase + thaw_amt, freeze_amt, canopy_temp())
            mt_ldew_rain(imt, ipatch) = mt_ldew_rain(imt, ipatch) - mt_flux
            mt_ldew_snow(imt, ipatch) = mt_ldew_snow(imt, ipatch) + mt_flux
         ENDIF

         CALL mt_equilibrate_dissolved(imt, max(ldew_rain, 0._r8), &
            mt_ldew_rain(imt, ipatch), mt_canopy_solid(imt, ipatch))

         ! --- Wet-leaf two-way equilibrium exchange with ambient vapour ---
         ! Zero net water flux, non-zero net material flux, so this is booked
         ! into a dedicated boundary-flux accumulator rather than a_mt_precip
         ! (which is paired with a water flux).  Liquid only: exchange with a
         ! frozen canopy pool is orders of magnitude slower.
         IF (DEF_ISO_CANOPY_EQUILIBRATION > 0._r8 .and. &
             material_fractionation_active(imt) .and. &
             .not. material_is_nonvolatile_solute(imt) .and. &
             ldew_rain > mt_tiny) THEN
            equil_gain = material_equilibration_exchange(mt_ldew_rain(imt, ipatch), &
               max(ldew_rain, 0._r8), R_vapor, material_alpha_liq_vap(imt, canopy_temp()), &
               DEF_ISO_CANOPY_EQUILIBRATION)
            mt_ldew_rain(imt, ipatch) = max(mt_ldew_rain(imt, ipatch) + equil_gain, 0._r8)
            IF (allocated(a_mt_vapor_exchange)) THEN
               a_mt_vapor_exchange(imt, ipatch) = &
                  a_mt_vapor_exchange(imt, ipatch) + equil_gain
            ENDIF
         ENDIF

         ! --- Soil+snow layers: combined liquid+ice ---
         DO j = lb, nl_soil
            d_wliq = wliq_soisno(j) - wliq_soisno_bef(j)
            d_wice = wice_soisno(j) - wice_soisno_bef(j)

            ! Determine freeze/thaw amounts (internal transfer, not I/O)
            thaw_amt   = 0._r8  ! ice→liquid
            freeze_amt = 0._r8  ! liquid→ice

            IF (present(soil_thaw_mass_th) .or. present(soil_frzc_mass_th)) THEN
               IF (present(soil_thaw_mass_th)) THEN
                  thaw_amt = min(max(soil_thaw_mass_th(j), 0._r8), max(wice_soisno_bef(j), 0._r8))
               ENDIF
               IF (present(soil_frzc_mass_th)) THEN
                  ! Cap against the liquid that will actually be there: the
                  ! pre-THERMAL liquid plus whatever thaws first, matching the
                  ! Rayleigh denominator below. Capping on wliq_soisno_bef
                  ! alone silently truncated -- with a fully frozen layer it
                  ! dropped the freeze entirely, however much had just thawed.
                  ! CoLM's own MOD_PhaseChange reports thaw and freeze as the
                  ! positive and negative parts of one ice difference, so only
                  ! one is ever non-zero and this changes nothing today; the
                  ! interface accepts both independently and must be coherent
                  ! for any caller that supplies them that way.
                  freeze_amt = min(max(soil_frzc_mass_th(j), 0._r8), &
                                   max(wliq_soisno_bef(j), 0._r8) + thaw_amt)
               ENDIF
            ELSE
               IF (d_wice < -mt_tiny .and. d_wliq > mt_tiny) THEN
                  thaw_amt = min(abs(d_wice), d_wliq)
               ENDIF
               IF (d_wliq < -mt_tiny .and. d_wice > mt_tiny) THEN
                  freeze_amt = min(abs(d_wliq), d_wice)
               ENDIF
            ENDIF

            ! --- THAW: ice → liquid (internal transfer) ---
            IF (thaw_amt > mt_tiny) THEN
               IF (wice_soisno_bef(j) > mt_tiny) THEN
                  ratio = mt_wice_soisno(imt, j, ipatch) / wice_soisno_bef(j)
                  mt_flux = min(thaw_amt * ratio, max(mt_wice_soisno(imt, j, ipatch), 0._r8))
                  mt_wice_soisno(imt, j, ipatch) = mt_wice_soisno(imt, j, ipatch) - mt_flux
                  mt_wliq_soisno(imt, j, ipatch) = mt_wliq_soisno(imt, j, ipatch) + mt_flux
               ENDIF
            ENDIF

            ! --- FREEZE: liquid → ice (internal transfer) ---
            IF (freeze_amt > mt_tiny) THEN
               ! Denominator is the liquid pool the freeze actually draws from,
               ! which is the pre-THERMAL liquid PLUS whatever just thawed --
               ! the same quantity wliq_post_phase is built from a few lines
               ! below, and the same convention the canopy path uses
               ! (ldew_rain_pre_phase + thaw_amt). The numerator
               ! mt_wliq_soisno has already been incremented by the thaw flux
               ! above, so pairing it with the pre-thaw liquid overstated the
               ! ratio by (wliq_bef + thaw)/wliq_bef. Only reachable when
               ! THERMAL reports thaw and freeze for the same layer in one
               ! step, which the explicit soil_thaw_mass_th / soil_frzc_mass_th
               ! path allows; with thaw_amt = 0 this is exactly the old
               ! expression.
               IF (wliq_soisno_bef(j) + thaw_amt > mt_tiny) THEN
                  mt_flux = material_rayleigh_freezing_loss(imt, mt_wliq_soisno(imt, j, ipatch), &
                     wliq_soisno_bef(j) + thaw_amt, freeze_amt, layer_temp(j))
                  mt_wliq_soisno(imt, j, ipatch) = mt_wliq_soisno(imt, j, ipatch) - mt_flux
                  mt_wice_soisno(imt, j, ipatch) = mt_wice_soisno(imt, j, ipatch) + mt_flux
               ENDIF
            ENDIF

            ! Reconstruct the water pool sizes that participate in the external
            ! phase changes below: pre-evap liquid = pre-THERMAL liquid + thaw
            ! - freeze, pre-sublimation ice = pre-THERMAL ice - thaw + freeze.
	            wliq_post_phase = max(wliq_soisno_bef(j) + thaw_amt - freeze_amt, 0._r8)
	            wice_post_phase = max(wice_soisno_bef(j) - thaw_amt + freeze_amt, 0._r8)

	            ! THERMAL reports only internal phase changes for snow layers.
	            ! Atmospheric snow exchange (qsdew/qseva/qfros/qsubl) is applied
	            ! later by snowwater and mirrored in mt_soil_water. Treating
	            ! any remaining snow-layer residual here as dew/frost/evap/subl
	            ! would create material fluxes with no matching water-side flux.
	            IF (j < 1) THEN
               CALL mt_equilibrate_dissolved(imt, max(wliq_soisno(j), 0._r8), &
                  mt_wliq_soisno(imt, j, ipatch), mt_solid_soisno(imt, j, ipatch))
               CYCLE
            ENDIF

	            ! --- Net liquid change beyond thaw/freeze = evap/dew ---
            ! Net external liquid change = d_wliq - thaw + freeze
            !   (thaw adds liquid internally, freeze removes liquid internally)
	            mt_flux = d_wliq - thaw_amt + freeze_amt
	            IF (mt_flux < -mt_tiny) THEN
	               IF (j == 1) THEN
	                  ! Net exposed-soil liquid loss = evaporation.
	                  IF (wliq_post_phase > mt_tiny) THEN
	                     water_loss = abs(mt_flux)
	                     mt_flux = evaporative_mt_loss(mt_wliq_soisno(imt, j, ipatch), &
	                        wliq_post_phase, water_loss, layer_temp(j), .false.)
	                     mt_wliq_soisno(imt, j, ipatch) = mt_wliq_soisno(imt, j, ipatch) - mt_flux
	                     CALL material_book_evap_loss(imt, ipatch, mt_flux, water_loss, &
	                        MT_EVAP_KIND_SOILEVAP)
	                  ENDIF
	               ELSE ! DEEP_STORAGE_RESIDUAL liquid loss: not atmospheric soil evaporation.
	                  IF (wliq_post_phase > mt_tiny) THEN
	                     ratio = mt_wliq_soisno(imt, j, ipatch) / wliq_post_phase
	                  ELSE
	                     ratio = material_init_water_ratio(imt)
	                  ENDIF
	                  mt_resid = min(abs(mt_flux) * max(ratio, 0._r8), &
	                     max(mt_wliq_soisno(imt, j, ipatch), 0._r8))
	                  mt_wliq_soisno(imt, j, ipatch) = mt_wliq_soisno(imt, j, ipatch) - mt_resid
	                  IF (allocated(mt_numerical_residual_step)) THEN
	                     mt_numerical_residual_step(imt, ipatch) = &
	                        mt_numerical_residual_step(imt, ipatch) - mt_resid
	                  ENDIF
	               ENDIF
            ELSEIF (mt_flux > mt_tiny) THEN
               ! Net liquid gain = dew deposition, only on the exposed soil
               ! surface. Deeper positive residuals are internal storage
               ! redistribution and must not be recoloured as atmospheric dew.
               IF (j == 1) THEN
                  ratio = deposition_ratio_for(layer_temp(j), .false.)
                  mt_wliq_soisno(imt, j, ipatch) = mt_wliq_soisno(imt, j, ipatch) + mt_flux * ratio
                  a_mt_precip(imt, ipatch) = a_mt_precip(imt, ipatch) + mt_flux * ratio
                  a_water_precip(imt, ipatch) = a_water_precip(imt, ipatch) + mt_flux
               ELSE ! DEEP_STORAGE_RESIDUAL liquid gain: preserve local/fallback signature.
                  IF (wliq_post_phase > mt_tiny) THEN
                     ratio = mt_wliq_soisno(imt, j, ipatch) / wliq_post_phase
                  ELSE
                     ratio = material_init_water_ratio(imt)
                  ENDIF
                  mt_resid = mt_flux * max(ratio, 0._r8)
                  mt_wliq_soisno(imt, j, ipatch) = mt_wliq_soisno(imt, j, ipatch) + mt_resid
                  IF (allocated(mt_numerical_residual_step)) THEN
                     mt_numerical_residual_step(imt, ipatch) = &
                        mt_numerical_residual_step(imt, ipatch) + mt_resid
                  ENDIF
               ENDIF
            ENDIF

            ! --- Net ice change beyond thaw/freeze = sublimation/frost ---
	            mt_flux = d_wice + thaw_amt - freeze_amt
	            IF (mt_flux < -mt_tiny) THEN
	               IF (j == 1) THEN
	                  ! Net exposed-soil ice loss = sublimation.
	                  IF (wice_post_phase > mt_tiny) THEN
	                     water_loss = abs(mt_flux)
	                     mt_flux = evaporative_mt_loss(mt_wice_soisno(imt, j, ipatch), &
	                        wice_post_phase, water_loss, layer_temp(j), .true.)
	                     mt_wice_soisno(imt, j, ipatch) = mt_wice_soisno(imt, j, ipatch) - mt_flux
	                     CALL material_book_evap_loss(imt, ipatch, mt_flux, water_loss, &
	                        MT_EVAP_KIND_SUBL)
	                  ENDIF
	               ELSE ! DEEP_STORAGE_RESIDUAL ice loss: not atmospheric sublimation.
	                  IF (wice_post_phase > mt_tiny) THEN
	                     ratio = mt_wice_soisno(imt, j, ipatch) / wice_post_phase
	                  ELSE
	                     ratio = material_init_water_ratio(imt)
	                  ENDIF
	                  mt_resid = min(abs(mt_flux) * max(ratio, 0._r8), &
	                     max(mt_wice_soisno(imt, j, ipatch), 0._r8))
	                  mt_wice_soisno(imt, j, ipatch) = mt_wice_soisno(imt, j, ipatch) - mt_resid
	                  IF (allocated(mt_numerical_residual_step)) THEN
	                     mt_numerical_residual_step(imt, ipatch) = &
	                        mt_numerical_residual_step(imt, ipatch) - mt_resid
	                  ENDIF
	               ENDIF
            ELSEIF (mt_flux > mt_tiny) THEN
               ! Net ice gain = frost deposition, only on the exposed soil
               ! surface. Deeper positive residuals are internal storage
               ! redistribution and must not be recoloured as atmospheric frost.
               IF (j == 1) THEN
                  ratio = deposition_ratio_for(layer_temp(j), .true.)
                  mt_wice_soisno(imt, j, ipatch) = mt_wice_soisno(imt, j, ipatch) + mt_flux * ratio
                  a_mt_precip(imt, ipatch) = a_mt_precip(imt, ipatch) + mt_flux * ratio
                  a_water_precip(imt, ipatch) = a_water_precip(imt, ipatch) + mt_flux
               ELSE ! DEEP_STORAGE_RESIDUAL ice gain: preserve local/fallback signature.
                  IF (wice_post_phase > mt_tiny) THEN
                     ratio = mt_wice_soisno(imt, j, ipatch) / wice_post_phase
                  ELSE
                     ratio = material_init_water_ratio(imt)
                  ENDIF
                  mt_resid = mt_flux * max(ratio, 0._r8)
                  mt_wice_soisno(imt, j, ipatch) = mt_wice_soisno(imt, j, ipatch) + mt_resid
                  IF (allocated(mt_numerical_residual_step)) THEN
                     mt_numerical_residual_step(imt, ipatch) = &
                        mt_numerical_residual_step(imt, ipatch) + mt_resid
                  ENDIF
               ENDIF
            ENDIF
            CALL mt_equilibrate_dissolved(imt, max(wliq_soisno(j), 0._r8), &
               mt_wliq_soisno(imt, j, ipatch), mt_solid_soisno(imt, j, ipatch))
         ENDDO
      ENDDO

      CONTAINS

      real(r8) FUNCTION evaporative_mt_loss (pool_mt, pool_water, water_loss, temp_k, from_ice)
         real(r8), intent(in) :: pool_mt
         real(r8), intent(in) :: pool_water
         real(r8), intent(in) :: water_loss
         real(r8), intent(in) :: temp_k
         logical,  intent(in) :: from_ice

         IF (material_is_nonvolatile_solute(imt)) THEN
            evaporative_mt_loss = 0._r8
         ELSE
            ! r_max ceiling only for fractionating isotope materials; above the
            ! ceiling, evaporation becomes non-fractionating so the pool R stops
            ! increasing.  Nonvolatile solutes are handled above as zero
            ! evaporative material loss.
            IF (from_ice) THEN
               ! Sublimation: only a thin surface skin exchanges isotopically;
               ! deeper ice leaves bodily as the surface retreats.
               evaporative_mt_loss = material_skin_limited_mt_loss(pool_mt, pool_water, &
                  water_loss, DEF_ISO_SUBL_SKIN_MM, temp_k, from_ice, evap_affine_for, mt_tiny, &
                  merge(materials(imt)%ref_ratio * (1._r8 + mt_delta_sanity_max / 1000._r8), &
                        0._r8, material_fractionation_active(imt)))
            ELSE
               evaporative_mt_loss = material_evaporative_mt_loss(pool_mt, pool_water, &
                  water_loss, temp_k, from_ice, evap_affine_for, mt_tiny, &
                  merge(materials(imt)%ref_ratio * (1._r8 + mt_delta_sanity_max / 1000._r8), &
                        0._r8, material_fractionation_active(imt)))
            ENDIF
         ENDIF
      END FUNCTION evaporative_mt_loss

      real(r8) FUNCTION canopy_temp ()
         canopy_temp = 273.15_r8
         IF (present(tleaf_frac)) canopy_temp = tleaf_frac
      END FUNCTION canopy_temp

      real(r8) FUNCTION layer_temp (jlay)
         integer, intent(in) :: jlay

         layer_temp = 273.15_r8
         IF (present(t_soisno_frac)) layer_temp = t_soisno_frac(jlay)
      END FUNCTION layer_temp

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
         alpha_k = material_alpha_kinetic_craig_gordon(imt, from_ice)
         CALL material_craig_gordon_evap_affine(imt, R_vapor, temp_k, relhum, &
            alpha_k, from_ice, slope, offset)
      END SUBROUTINE evap_affine_for

      real(r8) FUNCTION deposition_ratio_for (temp_k, from_ice)
         real(r8), intent(in) :: temp_k
         logical,  intent(in) :: from_ice

         IF (material_is_nonvolatile_solute(imt)) THEN
            deposition_ratio_for = 0._r8
            RETURN
         ENDIF
         deposition_ratio_for = R_atm
         IF (.not. material_fractionation_active(imt)) RETURN
         deposition_ratio_for = material_equilibrium_deposition_ratio(imt, R_vapor, temp_k, from_ice)
      END FUNCTION deposition_ratio_for

   END SUBROUTINE mt_evapo

END MODULE MOD_MT_Evapo
#endif
