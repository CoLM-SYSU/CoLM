#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_Precip

   USE MOD_Precision
   USE MOD_MT_Defs, only: nmts, mt_tiny, mt_water_min_for_ratio, mt_uses_generic_water
   USE MOD_MT_Solute, only: mt_equilibrate_dissolved
   USE MOD_MT_Forcing, only: material_forcing_rain_value, material_forcing_snow_value
   USE MOD_MT_Vars, only: mt_ldew_rain, mt_ldew_snow, a_mt_precip, &
      mt_pg_rain_ground, mt_pg_snow_ground, mt_waterstorage, &
      mt_canopy_solid, mt_surface_solid, mt_waterstorage_solid
   USE MOD_Const_Physical, only: tfrz
   USE MOD_Namelist, only: DEF_VEG_SNOW

   IMPLICIT NONE


CONTAINS

   !---------------------------------------------------------------
   ! Canopy material update after LEAF_INTERCEPTION.
   !
   ! Physical process:
   !   1. Precipitation arrives: forc_rain, forc_snow
   !   2. Part intercepted (qintr): mixes with EXISTING canopy water
   !   3. Part passes through (throughfall): keeps precipitation R
   !   4. If canopy saturates, excess drips: carries MIXED canopy R
   !   5. Phase change inside LEAF_INTERCEPTION (NoahMP/MATSIRO/VIC/
   !      JULES) moves ldew_smelt_out from ldew_snow to ldew_rain and/or
   !      ldew_frzc_out from ldew_rain to ldew_snow. Material mass has to
   !      migrate coherently, otherwise the rain pool loses signal from
   !      the melt event (and vice versa for freeze) and the snow material
   !      gets mis-attributed as drip to mt_pg_snow_ground downstream.
   !
   ! Material to ground = throughfall*R_precip + drip*R_mixed
   !   (rain → mt_pg_rain_ground → surface pool via mt_soil_water;
   !    snow → mt_pg_snow_ground → mt_scv via mt_newsnow)
   !
   ! Canopy material: migrate phase change first (updates "effective"
   ! pre-interception anchor so the rain/snow stages below see a
   ! consistent (ldew_*_old, mt_ldew_*) pair), then for each phase
   ! remove xsc, mix interception with remainder, and finally strip drip.
   !---------------------------------------------------------------
   SUBROUTINE mt_precip (ipatch, deltim, &
      forc_rain, forc_snow, qintr, qintr_rain, qintr_snow, &
      pg_rain, pg_snow, ldew_rain, ldew_snow, &
      ldew_rain_old, ldew_snow_old, qflx_irrig_sprinkler, &
      gross_intr_rain, gross_intr_snow, &
      xsc_rain_out, xsc_snow_out, &
      ldew_smelt_mass, ldew_frzc_mass, waterstorage_before)

      IMPLICIT NONE
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: deltim
      real(r8), intent(in) :: forc_rain, forc_snow
      real(r8), intent(in) :: qintr, qintr_rain, qintr_snow
      real(r8), intent(in) :: pg_rain, pg_snow
      real(r8), intent(in) :: ldew_rain, ldew_snow         ! AFTER interception
      real(r8), intent(in) :: ldew_rain_old, ldew_snow_old ! BEFORE interception
      ! Sprinkler irrigation rate (mm/s). LEAF_INTERCEPTION merges this
      ! into the rain stream (MOD_LeafInterception.F90:187,232), so the
      ! "total rain arriving" at the patch is (forc_rain + this). Without
      ! counting it here the material input is under-reported and the
      ! throughfall-vs-interception partition uses the wrong denominator.
      real(r8), intent(in) :: qflx_irrig_sprinkler
      ! Gross rain/snow rate entering the canopy mixed pool [mm/s], >=0.
      ! Returned by LEAF_INTERCEPTION alongside qintr (net). Used here to
      ! correctly attribute isotope signatures on canopy-release events:
      ! schemes 3/5/6/7 can have qintr_rain < 0 when stored canopy water
      ! drips out while new rain simultaneously enters the pool. Using
      ! only qintr (net) would attribute the released canopy water to the
      ! fresh-precipitation signature R_snow_input / R_rain_input. Using
      ! gross_intr_* gives the right R_mixed for the drip and the fresh
      ! ratio for the pure-gap part.
      real(r8), intent(in) :: gross_intr_rain
      real(r8), intent(in) :: gross_intr_snow
      ! Pre-mix old-pool release rate [mm/s], >=0. In most schemes a
      ! portion of the OLD canopy pool is flushed BEFORE new precipitation
      ! mixes in (capacity overflow on LAI-driven satcap shrink, phase-
      ! change overflow, snow unloading in NoahMP, etc). That mass carries
      ! the pre-mix canopy signature, NOT R_mixed. Without separating it
      ! out, the material drip term would lump xsc with the post-mix tex
      ! and dilute old canopy signatures with fresh precipitation.
      real(r8), intent(in) :: xsc_rain_out
      real(r8), intent(in) :: xsc_snow_out
      ! Canopy rain<->snow phase change mass [grid-scale mm, >=0].
      ! ldew_smelt_mass is the canopy snow -> rain transfer this step
      ! (melt); ldew_frzc_mass is the canopy rain -> snow transfer
      ! (freeze). Each nonzero value triggers an equal-mass material
      ! migration between mt_ldew_snow and mt_ldew_rain before the
      ! rain/snow stages run, and the pre-interception anchors
      ! (ldew_*_old_eff, mt_ldew_*) are updated so R_canopy_pre in each
      ! stage reflects the post-migration pool composition. Without this
      ! step the snow-pool material would be spuriously classified as drip
      ! and sent to mt_pg_snow_ground (and the rain pool would
      ! under-report the melt signal).
      real(r8), intent(in) :: ldew_smelt_mass
      real(r8), intent(in) :: ldew_frzc_mass
      ! Reservoir water before sprinkler withdrawal. The bulk irrigation
      ! routine has already debited waterstorage by the time this routine runs.
      real(r8), intent(in), optional :: waterstorage_before

      integer  :: imt
      real(r8) :: R_snow_input, R_rain_atm, storage_ratio, R_rain_input
      real(r8) :: intercepted, drip, throughfall
      real(r8) :: mt_mixed, water_mixed, R_mixed
      real(r8) :: mt_throughfall, mt_drip
      real(r8) :: mt_rain_ground, mt_snow_ground
      real(r8) :: rain_total           ! forc_rain + sprinkler [mm/s]
      real(r8) :: sprinkler_rate       ! non-negative sprinkler water input [mm/s]
      real(r8) :: sprinkler_water      ! sprinkler water transferred this step [mm]
      real(r8) :: xsc_mass             ! pre-mix release in this step [mm]
      real(r8) :: mt_xsc              ! pre-mix release material mass
      real(r8) :: R_canopy_pre         ! canopy signature BEFORE mix
      real(r8) :: ldew_rain_pre_mix    ! post-xsc, pre-mix storage [mm]
      real(r8) :: ldew_snow_pre_mix
      ! Effective pre-interception anchors after phase-change migration.
      ! Equal to ldew_*_old when no phase change happened; otherwise
      ! reflect the post-migration pool so R_canopy_pre stays consistent
      ! with the (migrated) mt_ldew_* values that feed the rain/snow
      ! stages. Clamped to >=0 in case of float noise in the caller.
      real(r8) :: ldew_rain_old_eff
      real(r8) :: ldew_snow_old_eff
      real(r8) :: smelt_mass, frzc_mass ! effective phase-change masses
      real(r8) :: R_snow_pre, R_rain_pre
      real(r8) :: mt_smelt, mt_frzc   ! material mass migrated this step
      real(r8) :: canopy_mt_beg, canopy_mt_end, canopy_input_mt
      real(r8) :: canopy_resid, ground_mt_total, desired_ground_total
      real(r8) :: sprinkler_mt
      real(r8) :: canopy_old_water, canopy_release_mt
      real(r8) :: release_rain, release_snow, release_total

      IF (nmts <= 0) RETURN

      ! Interception / pg_rain / ldew already see the sprinkler-boosted
      ! rain stream from LEAF_INTERCEPTION. Use the boosted total here
      ! so the throughfall computation matches (intercepted + throughfall
      ! ≡ rain_total*deltim) and the system-input accounting covers
      ! both natural precipitation and sprinkler irrigation.
      sprinkler_rate = max(qflx_irrig_sprinkler, 0._r8)
      sprinkler_water = sprinkler_rate * deltim
      rain_total = max(forc_rain, 0._r8) + sprinkler_rate

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         ! Rain and snow take their own forcing ratio: they differ only where
         ! the forcing record is near-dry (see material_forcing_rain_value in
         ! MOD_MT_Forcing); elsewhere both are the same value as before.
         R_rain_atm   = material_forcing_rain_value(imt, ipatch)
         R_snow_input = material_forcing_snow_value(imt, ipatch)
         storage_ratio = R_rain_atm
         CALL mt_equilibrate_dissolved(imt, max(ldew_rain_old, 0._r8), &
            mt_ldew_rain(imt, ipatch), mt_canopy_solid(imt, ipatch))
         canopy_mt_beg = mt_ldew_rain(imt, ipatch) + mt_ldew_snow(imt, ipatch) &
            + mt_canopy_solid(imt, ipatch)

         ! Atmospheric precipitation is external input; count into the
         ! precip accumulator. Sprinkler irrigation also reaches the
         ! canopy/ground here via LEAF_INTERCEPTION, but water-side
         ! treats it as an internal transfer from waterstorage (not a
         ! term in forc_prc+forc_prl). Mirror that: debit
         ! mt_waterstorage for sprinkler and keep a_mt_precip purely
         ! atmospheric, so the material budget's "external input"
         ! matches the water-side definition.
         ! The second term is exactly zero when rain and snow share a ratio,
         ! which keeps this sum bit-identical to the single-ratio form.
         a_mt_precip(imt, ipatch) = a_mt_precip(imt, ipatch) &
            + (forc_rain + forc_snow) * R_snow_input * deltim &
            + forc_rain * (R_rain_atm - R_snow_input) * deltim

         IF (qflx_irrig_sprinkler > mt_tiny .and. allocated(mt_waterstorage)) THEN
            IF (present(waterstorage_before)) THEN
               CALL mt_equilibrate_dissolved(imt, max(waterstorage_before, 0._r8), &
                  mt_waterstorage(imt, ipatch), mt_waterstorage_solid(imt, ipatch))
               IF (waterstorage_before > mt_water_min_for_ratio) THEN
                  storage_ratio = mt_waterstorage(imt, ipatch) / waterstorage_before
               ENDIF
            ENDIF
            sprinkler_mt = sprinkler_water * storage_ratio
            sprinkler_mt = min(max(sprinkler_mt, 0._r8), max(mt_waterstorage(imt, ipatch), 0._r8))
            mt_waterstorage(imt, ipatch) = mt_waterstorage(imt, ipatch) &
               - sprinkler_mt
            ! The canopy must receive exactly the material actually removed,
            ! including near-dry and dissolved-inventory-limited withdrawals.
            IF (sprinkler_water > mt_tiny) storage_ratio = sprinkler_mt / sprinkler_water
            IF (present(waterstorage_before)) THEN
               CALL mt_equilibrate_dissolved(imt, &
                  max(waterstorage_before - sprinkler_water, 0._r8), &
                  mt_waterstorage(imt, ipatch), mt_waterstorage_solid(imt, ipatch))
            ENDIF
         ELSE
            sprinkler_mt = sprinkler_water * storage_ratio
         ENDIF

         R_rain_input = R_rain_atm
         IF (rain_total > mt_tiny) THEN
            R_rain_input = (max(forc_rain, 0._r8) * R_rain_atm + sprinkler_rate * storage_ratio) / rain_total
         ENDIF

         mt_rain_ground = 0._r8
         mt_snow_ground = 0._r8

         ! ---- Stage 0: canopy phase-change migration ----
         ! NoahMP/MATSIRO/VIC/JULES may have melted part of the canopy
         ! snow pool into the rain pool (ldew_smelt_mass > 0) or frozen
         ! part of the rain pool into the snow pool (ldew_frzc_mass > 0)
         ! BEFORE touching any new precipitation this step. Without
         ! migrating the material mass along with the bulk water, the snow
         ! block below would see d_ldew < 0 (for melt) and falsely
         ! classify ldew_smelt_mass as drip to mt_pg_snow_ground, while
         ! the rain pool would silently lose the meltwater's signature.
         ! Clamp each mass to the pool available when that transfer runs;
         ! melt precedes freeze in this aggregate, well-mixed canopy state.
         ldew_rain_old_eff = max(0._r8, ldew_rain_old)
         ldew_snow_old_eff = max(0._r8, ldew_snow_old)
         smelt_mass = min(max(0._r8, ldew_smelt_mass), ldew_snow_old_eff)
         frzc_mass  = max(0._r8, ldew_frzc_mass)

         IF (smelt_mass > 0._r8) THEN
            IF (ldew_snow_old_eff > mt_tiny) THEN
               R_snow_pre = mt_ldew_snow(imt, ipatch) / ldew_snow_old_eff
            ELSE
               R_snow_pre = R_snow_input
            ENDIF
            mt_smelt = smelt_mass * R_snow_pre
            mt_ldew_snow(imt, ipatch) = max(0._r8, &
               mt_ldew_snow(imt, ipatch) - mt_smelt)
            mt_ldew_rain(imt, ipatch) = mt_ldew_rain(imt, ipatch) + mt_smelt
            ldew_snow_old_eff = ldew_snow_old_eff - smelt_mass
            ldew_rain_old_eff = ldew_rain_old_eff + smelt_mass
         ENDIF

         frzc_mass = min(frzc_mass, ldew_rain_old_eff)
         IF (frzc_mass > 0._r8) THEN
            ! ponytail: interception exposes aggregate phase mass but no
            ! freeze-event temperature (and PFTs can differ); keep this
            ! transfer conservative until temperature-weighted diagnostics exist.
            IF (ldew_rain_old_eff > mt_tiny) THEN
               R_rain_pre = mt_ldew_rain(imt, ipatch) / ldew_rain_old_eff
            ELSE
               R_rain_pre = R_rain_atm
            ENDIF
            mt_frzc = frzc_mass * R_rain_pre
            mt_ldew_rain(imt, ipatch) = max(0._r8, &
               mt_ldew_rain(imt, ipatch) - mt_frzc)
            mt_ldew_snow(imt, ipatch) = mt_ldew_snow(imt, ipatch) + mt_frzc
            ldew_rain_old_eff = ldew_rain_old_eff - frzc_mass
            ldew_snow_old_eff = ldew_snow_old_eff + frzc_mass
         ENDIF

         ! LEAF_interception has a carrier-empty branch for lai+sai~=0 that
         ! releases the entire old single-bucket canopy store as rain or snow
         ! according to tleaf, then zeros both ldew phase pools.  Its xsc
         ! diagnostic therefore names the RELEASE phase, not necessarily the
         ! phase in which the old water/material entered this routine.  Repack the
         ! complete old canopy inventory into that host-selected phase before
         ! the ordinary xsc calculation.  Otherwise old snow released as rain
         ! remains stranded in mt_ldew_snow (and old rain released as snow is
         ! precipitated into mt_canopy_solid when the now-empty liquid carrier
         ! is equilibrated below).
         canopy_old_water = ldew_rain_old_eff + ldew_snow_old_eff
         IF (max(ldew_rain, 0._r8) + max(ldew_snow, 0._r8) <= mt_water_min_for_ratio .and. &
             canopy_old_water > mt_water_min_for_ratio .and. &
             (((max(xsc_rain_out, 0._r8) * deltim > mt_water_min_for_ratio) .and. &
                (max(xsc_snow_out, 0._r8) * deltim <= mt_water_min_for_ratio)) .or. &
              ((max(xsc_snow_out, 0._r8) * deltim > mt_water_min_for_ratio) .and. &
                (max(xsc_rain_out, 0._r8) * deltim <= mt_water_min_for_ratio))) .and. &
             (max(xsc_rain_out, 0._r8) + max(xsc_snow_out, 0._r8)) * deltim &
                >= canopy_old_water - mt_water_min_for_ratio) THEN
            canopy_release_mt = mt_ldew_rain(imt, ipatch) + &
               mt_ldew_snow(imt, ipatch) + mt_canopy_solid(imt, ipatch)
            mt_canopy_solid(imt, ipatch) = 0._r8
            IF (max(xsc_rain_out, 0._r8) * deltim > mt_water_min_for_ratio) THEN
               mt_ldew_rain(imt, ipatch) = canopy_release_mt
               mt_ldew_snow(imt, ipatch) = 0._r8
               ldew_rain_old_eff = canopy_old_water
               ldew_snow_old_eff = 0._r8
            ELSE
               mt_ldew_rain(imt, ipatch) = 0._r8
               mt_ldew_snow(imt, ipatch) = canopy_release_mt
               ldew_rain_old_eff = 0._r8
               ldew_snow_old_eff = canopy_old_water
            ENDIF
         ENDIF

         ! A single bucket can release both phases: PC/PFT interception uses
         ! each PFT's leaf temperature, not the patch mean used by canopy_lump.
         ! When a release lacks its nominal phase donor, reserve that water
         ! from the COMMON old inventory before the ordinary two-stage update.
         ! Keep the remainder in its original phase and leave already-matching
         ! donors untouched (including the homogeneous IGBP path).
         release_rain = max(xsc_rain_out, 0._r8) * deltim
         release_snow = max(xsc_snow_out, 0._r8) * deltim
         IF (.not. DEF_VEG_SNOW .and. canopy_old_water > mt_tiny .and. &
             (release_rain > ldew_rain_old_eff .or. release_snow > ldew_snow_old_eff)) THEN
            release_total = release_rain + release_snow
            IF (release_total > canopy_old_water) THEN
               release_rain = canopy_old_water * (release_rain / release_total)
               release_snow = canopy_old_water - release_rain
            ENDIF
            canopy_release_mt = mt_ldew_rain(imt, ipatch) + mt_ldew_snow(imt, ipatch)
            ! Retain the established empty-canopy flush of solid inventory,
            ! also when the complete release has BOTH phases.
            IF (max(ldew_rain, 0._r8) + max(ldew_snow, 0._r8) <= mt_water_min_for_ratio .and. &
                release_total >= canopy_old_water - mt_water_min_for_ratio) THEN
               canopy_release_mt = canopy_release_mt + mt_canopy_solid(imt, ipatch)
               mt_canopy_solid(imt, ipatch) = 0._r8
            ENDIF
            IF (release_snow > ldew_snow_old_eff) THEN
               ldew_snow_old_eff = release_snow
               ldew_rain_old_eff = canopy_old_water - ldew_snow_old_eff
            ELSE
               ldew_rain_old_eff = release_rain
               ldew_snow_old_eff = canopy_old_water - ldew_rain_old_eff
            ENDIF
            ! Compute the smaller share directly, avoiding cancellation for
            ! tiny releases from a much larger canopy reservoir.
            IF (ldew_rain_old_eff <= ldew_snow_old_eff) THEN
               mt_ldew_rain(imt, ipatch) = canopy_release_mt * (ldew_rain_old_eff / canopy_old_water)
               mt_ldew_snow(imt, ipatch) = canopy_release_mt - mt_ldew_rain(imt, ipatch)
            ELSE
               mt_ldew_snow(imt, ipatch) = canopy_release_mt * (ldew_snow_old_eff / canopy_old_water)
               mt_ldew_rain(imt, ipatch) = canopy_release_mt - mt_ldew_snow(imt, ipatch)
            ENDIF
         ENDIF

         ! ---- Rain component (2-stage isotope attribution) ----
         ! Stage 1: PRE-MIX release. xsc_rain_out is the fraction of the
         ! OLD canopy pool flushed out before any new rain mixes in
         ! (capacity overflow, phase-change overflow, unloading). It
         ! carries the canopy's pre-mix signature R_canopy_pre, derived
         ! from the entry material storage and entry bulk water.
         ! Use the post-phase-change effective old-pool so R_canopy_pre
         ! reflects the actual pool entering the xsc stage.
         xsc_mass = min(max(0._r8, xsc_rain_out) * deltim, ldew_rain_old_eff)
         IF (ldew_rain_old_eff > mt_tiny) THEN
            R_canopy_pre = mt_ldew_rain(imt, ipatch) / ldew_rain_old_eff
         ELSE
            R_canopy_pre = R_rain_input
         ENDIF
         mt_xsc = xsc_mass * R_canopy_pre
         ldew_rain_pre_mix = ldew_rain_old_eff - xsc_mass

         ! Stage 2: gross rain enters the (post-release) canopy pool at
         ! the mixed atmospheric+sprinkler rain signature. Blend to obtain R_mixed.
         intercepted = max(0._r8, gross_intr_rain) * deltim
         water_mixed = ldew_rain_pre_mix + intercepted
         mt_mixed = (mt_ldew_rain(imt, ipatch) - mt_xsc) &
                   + intercepted * R_rain_input
         IF (water_mixed > mt_tiny) THEN
            R_mixed = mt_mixed / water_mixed
         ELSE
            R_mixed = R_rain_input
         ENDIF

         ! Pure-gap throughfall (never touched canopy).
         throughfall = max(rain_total * deltim - intercepted, 0._r8)

         ! Actual post-mix drip is defined by the water-side ground flux.
         ! This avoids relying on ldew_rain as a phase-specific final store:
         ! scheme=1 with DEF_VEG_SNOW=false updates only the single ldew
         ! bucket during interception.
         drip = max(max(pg_rain, 0._r8) * deltim - xsc_mass - throughfall, 0._r8)
         drip = min(drip, max(intercepted, 0._r8))

         ! Rain material reaching the ground: pure-gap throughfall keeps the
         ! input rain signature, while drip carries the mixed canopy ratio.
         mt_throughfall = throughfall * R_rain_input
         mt_drip = drip * R_mixed

         ! Canopy material: blended pool minus drip
         mt_ldew_rain(imt, ipatch) = max(mt_mixed - drip * R_mixed, 0._r8)

         mt_rain_ground = mt_throughfall + mt_xsc + mt_drip

         ! ---- Snow component (same 2-stage logic) ----
         ! Use ldew_snow_old_eff (post-phase-change effective anchor) and
         ! the migrated mt_ldew_snow so R_canopy_pre matches the pool
         ! that actually feeds the xsc stage.
         xsc_mass = min(max(0._r8, xsc_snow_out) * deltim, ldew_snow_old_eff)
         IF (ldew_snow_old_eff > mt_tiny) THEN
            R_canopy_pre = mt_ldew_snow(imt, ipatch) / ldew_snow_old_eff
         ELSE
            R_canopy_pre = R_snow_input
         ENDIF
         mt_xsc = xsc_mass * R_canopy_pre
         ldew_snow_pre_mix = ldew_snow_old_eff - xsc_mass

         intercepted = max(0._r8, gross_intr_snow) * deltim
         water_mixed = ldew_snow_pre_mix + intercepted
         mt_mixed = (mt_ldew_snow(imt, ipatch) - mt_xsc) &
                   + intercepted * R_snow_input
         IF (water_mixed > mt_tiny) THEN
            R_mixed = mt_mixed / water_mixed
         ELSE
            R_mixed = R_snow_input
         ENDIF

         throughfall = max(forc_snow * deltim - intercepted, 0._r8)
         drip = max(max(pg_snow, 0._r8) * deltim - xsc_mass - throughfall, 0._r8)
         drip = min(drip, max(water_mixed, 0._r8))

         mt_ldew_snow(imt, ipatch) = max(mt_mixed - drip * R_mixed, 0._r8)

         mt_snow_ground = throughfall * R_snow_input + mt_xsc + drip * R_mixed

         ! Close the combined canopy interception budget exactly:
         !   Δ(canopy material) = rain/snow input material - ground material.
         ! The phase-specific partition above follows water-side release
         ! diagnostics, but several interception schemes expose only net
         ! pg_* plus clamped xsc/gross terms. That can leave O(1e-8)
         ! residuals which later show up as wetland/soil balance errors.
         ! Preserve the prognostic canopy pools and put the tiny correction
         ! on the internal ground flux consumed by mt_soil_water /
         ! mt_wetland; atmospheric a_mt_precip remains unchanged.
         CALL mt_equilibrate_dissolved(imt, max(ldew_rain, 0._r8), &
            mt_ldew_rain(imt, ipatch), mt_canopy_solid(imt, ipatch))
         canopy_mt_end = mt_ldew_rain(imt, ipatch) + mt_ldew_snow(imt, ipatch) &
            + mt_canopy_solid(imt, ipatch)
         canopy_input_mt = rain_total * deltim * R_rain_input + max(forc_snow, 0._r8) * deltim * R_snow_input
         canopy_resid = canopy_mt_end - canopy_mt_beg - canopy_input_mt &
                      + mt_rain_ground + mt_snow_ground
         IF (abs(canopy_resid) > mt_tiny) THEN
            ground_mt_total = mt_rain_ground + mt_snow_ground
            desired_ground_total = max(ground_mt_total - canopy_resid, 0._r8)
            IF (ground_mt_total > mt_tiny) THEN
               mt_rain_ground = mt_rain_ground * desired_ground_total / ground_mt_total
               mt_snow_ground = mt_snow_ground * desired_ground_total / ground_mt_total
            ELSEIF (rain_total > mt_tiny) THEN
               mt_rain_ground = desired_ground_total
            ELSE
               mt_snow_ground = desired_ground_total
            ENDIF
         ENDIF

         ! Ground rain is a liquid mobile pool. Excess solute precipitates
         ! directly into the surface solid inventory; existing surface solid
         ! may redissolve only up to the rainwater capacity.
         CALL mt_equilibrate_dissolved(imt, max(pg_rain, 0._r8) * deltim, &
            mt_rain_ground, mt_surface_solid(imt, ipatch))

        ! Store rain/snow material separately: rain feeds the surface mixed
        ! pool in mt_soil_water; snow feeds mt_scv in mt_newsnow.
         mt_pg_rain_ground(imt, ipatch) = mt_rain_ground
         mt_pg_snow_ground(imt, ipatch) = mt_snow_ground

      ENDDO

   END SUBROUTINE mt_precip

   !---------------------------------------------------------------
   ! Canopy rain/snow pools when the host keeps a single canopy store.
   !
   ! With DEF_VEG_SNOW = .false. the host interception works on ldew only
   ! (excess and new water go to the rain or the snow side by tleaf) and
   ! ldew_rain/ldew_snow are re-split afterwards in the leaf temperature
   ! code. The host values are left untouched; only the material pools are
   ! arranged here so that they match the store the host is acting on.
   !
   ! mt_canopy_lump (before interception): put all canopy material into the
   !   phase the host interception uses and return that phase's water.
   ! mt_canopy_split (before the energy balance): split ldew the way the
   !   leaf temperature code splits it, and the material in the same ratio.
   ! With DEF_VEG_SNOW = .true. both return the host components unchanged.
   !---------------------------------------------------------------
   SUBROUTINE mt_canopy_lump (ipatch, ldew, ldew_rain, ldew_snow, tleaf, rain_eff, snow_eff)

      IMPLICIT NONE
      integer,  intent(in)  :: ipatch
      real(r8), intent(in)  :: ldew, ldew_rain, ldew_snow, tleaf
      real(r8), intent(out) :: rain_eff, snow_eff

      integer :: imt

         IF (DEF_VEG_SNOW) THEN
            rain_eff = ldew_rain
            snow_eff = ldew_snow
            RETURN
         ENDIF

         IF (tleaf > tfrz) THEN
            rain_eff = ldew
            snow_eff = 0._r8
         ELSE
            rain_eff = 0._r8
            snow_eff = ldew
         ENDIF

         DO imt = 1, nmts
            IF (.not. mt_uses_generic_water(imt)) CYCLE
            IF (tleaf > tfrz) THEN
               mt_ldew_rain(imt, ipatch) = mt_ldew_rain(imt, ipatch) + mt_ldew_snow(imt, ipatch)
               mt_ldew_snow(imt, ipatch) = 0._r8
            ELSE
               mt_ldew_snow(imt, ipatch) = mt_ldew_snow(imt, ipatch) + mt_ldew_rain(imt, ipatch)
               mt_ldew_rain(imt, ipatch) = 0._r8
            ENDIF
         ENDDO

   END SUBROUTINE mt_canopy_lump

   SUBROUTINE mt_canopy_split (ipatch, ldew, ldew_rain, ldew_snow, tleaf, rain_eff, snow_eff)

      IMPLICIT NONE
      integer,  intent(in)  :: ipatch
      real(r8), intent(in)  :: ldew, ldew_rain, ldew_snow, tleaf
      real(r8), intent(out) :: rain_eff, snow_eff

      integer  :: imt
      real(r8) :: total

         CALL mt_canopy_parts (ldew, ldew_rain, ldew_snow, tleaf, rain_eff, snow_eff)
         IF (DEF_VEG_SNOW) RETURN

         DO imt = 1, nmts
            IF (.not. mt_uses_generic_water(imt)) CYCLE
            total = mt_ldew_rain(imt, ipatch) + mt_ldew_snow(imt, ipatch)
            IF (ldew > 0._r8) THEN
               mt_ldew_rain(imt, ipatch) = total * (rain_eff / ldew)
               mt_ldew_snow(imt, ipatch) = total - mt_ldew_rain(imt, ipatch)
            ENDIF
         ENDDO

   END SUBROUTINE mt_canopy_split

   !---------------------------------------------------------------
   ! Rain and snow parts of ldew as the DEF_VEG_SNOW = .false. branch of
   ! LeafTemperature splits it; no material pool is changed. Also used after
   ! the energy balance: LeafTemperaturePC has no such re-split, so its
   ! ldew_rain/ldew_snow need not add up to ldew there.
   !---------------------------------------------------------------
   SUBROUTINE mt_canopy_parts (ldew, ldew_rain, ldew_snow, tleaf, rain_eff, snow_eff)

      IMPLICIT NONE
      real(r8), intent(in)  :: ldew, ldew_rain, ldew_snow, tleaf
      real(r8), intent(out) :: rain_eff, snow_eff

         IF (DEF_VEG_SNOW) THEN
            rain_eff = ldew_rain
            snow_eff = ldew_snow
         ELSEIF (ldew_rain + ldew_snow > 1.e-10) THEN
            rain_eff = ldew * (ldew_rain / (ldew_rain + ldew_snow))
            snow_eff = ldew - rain_eff
         ELSEIF (ldew > 0.) THEN
            IF (tleaf > tfrz) THEN
               rain_eff = ldew
               snow_eff = 0._r8
            ELSE
               rain_eff = 0._r8
               snow_eff = ldew
            ENDIF
         ELSE
            rain_eff = 0._r8
            snow_eff = 0._r8
         ENDIF

   END SUBROUTINE mt_canopy_parts

END MODULE MOD_MT_Precip
#endif
