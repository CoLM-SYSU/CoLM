#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_EvapLimit
!=======================================================================
! Shared finite-pool evaporation/sublimation limiter for material losses.
! Callers provide the evaporate ratio of the pool as an affine function of
! the pool ratio (R_E = slope*R + offset, Craig-Gordon); the integration
! over the step's water loss lives here so canopy, soil and wetland paths
! cannot silently diverge.
!=======================================================================

   USE MOD_Precision

   IMPLICIT NONE
   PRIVATE

   abstract interface
      SUBROUTINE material_evap_affine_callback (temp_k, from_ice, slope, offset)
         IMPORT :: r8
         real(r8), intent(in)  :: temp_k
         logical,  intent(in)  :: from_ice
         real(r8), intent(out) :: slope, offset
      END SUBROUTINE material_evap_affine_callback
   end interface

   PUBLIC :: material_evaporative_mt_loss
   PUBLIC :: material_atmospheric_mt_loss
   PUBLIC :: material_skin_limited_mt_loss

CONTAINS

   real(r8) FUNCTION material_atmospheric_mt_loss (pool_mt, pool_water, water_loss, &
      temp_k, from_ice, evap_affine_fn, mt_tiny, r_max, is_nonvolatile, skin_mass)

      real(r8), intent(in) :: pool_mt
      real(r8), intent(in) :: pool_water
      real(r8), intent(in) :: water_loss
      real(r8), intent(in) :: temp_k
      logical,  intent(in) :: from_ice
      procedure(material_evap_affine_callback) :: evap_affine_fn
      real(r8), intent(in) :: mt_tiny
      real(r8), intent(in) :: r_max
      logical,  intent(in) :: is_nonvolatile
      ! Exchanging surface-skin mass for SUBLIMATION.  Supplied by callers that
      ! route snow/ice sublimation through the skin limiter; when omitted the
      ! whole pool is treated as well mixed (legacy behaviour).
      real(r8), intent(in), optional :: skin_mass

      IF (is_nonvolatile) THEN
         material_atmospheric_mt_loss = 0._r8
      ELSEIF (from_ice .and. present(skin_mass)) THEN
         material_atmospheric_mt_loss = material_skin_limited_mt_loss(pool_mt, pool_water, &
            water_loss, skin_mass, temp_k, from_ice, evap_affine_fn, mt_tiny, r_max)
      ELSE
         material_atmospheric_mt_loss = material_evaporative_mt_loss(pool_mt, pool_water, &
            water_loss, temp_k, from_ice, evap_affine_fn, mt_tiny, r_max)
      ENDIF
   END FUNCTION material_atmospheric_mt_loss

   !---------------------------------------------------------------
   ! Sublimation/evaporation where only a thin SURFACE SKIN of the pool can
   ! exchange isotopically within one step.
   !
   ! Applying Craig-Gordon to a whole snow layer treats tens of mm of ice as
   ! a single well-mixed reservoir and over-enriches the pack: physically
   ! only the top few mm exchange with the atmosphere, while deeper mass is
   ! removed bodily as the snow surface retreats, which does not fractionate.
   ! IsoGSM avoids the problem by not fractionating land/ice evaporation at
   ! all (gsml/moninp.F:384-389).
   !
   ! Split the step's water loss:
   !   w_frac = min(water_loss, skin_mass)  -> fractionates (Craig-Gordon)
   !   w_bulk = water_loss - w_frac         -> water-matched, no fractionation
   ! The skin's ratio is taken as the layer mean, i.e. the skin is assumed
   ! well mixed with its layer between steps.  That is deliberately
   ! conservative: resolving a skin that enriches independently would need a
   ! new prognostic sub-layer state (and restart/LULCC plumbing for it).
   !
   ! Limits: skin_mass -> infinity reproduces the layer-mixed behaviour
   ! exactly; skin_mass = 0 reproduces IsoGSM's non-fractionating limit.
   !---------------------------------------------------------------
   real(r8) FUNCTION material_skin_limited_mt_loss (pool_mt, pool_water, water_loss, &
      skin_mass, temp_k, from_ice, evap_affine_fn, mt_tiny, r_max)

      real(r8), intent(in) :: pool_mt
      real(r8), intent(in) :: pool_water
      real(r8), intent(in) :: water_loss
      real(r8), intent(in) :: skin_mass
      real(r8), intent(in) :: temp_k
      logical,  intent(in) :: from_ice
      procedure(material_evap_affine_callback) :: evap_affine_fn
      real(r8), intent(in) :: mt_tiny
      real(r8), intent(in) :: r_max

      real(r8) :: w_frac, w_bulk, frac_loss, bulk_loss
      real(r8) :: remaining_mt, remaining_water

      material_skin_limited_mt_loss = 0._r8
      IF (pool_water <= mt_tiny .or. water_loss <= mt_tiny) RETURN
      IF (water_loss >= pool_water * (1._r8 - 1.e-12_r8)) THEN
         ! Whole pool leaves: the integrated outgoing material is the whole
         ! inventory, independent of how much of it fractionated.
         material_skin_limited_mt_loss = max(pool_mt, 0._r8)
         RETURN
      ENDIF

      w_frac = min(water_loss, max(skin_mass, 0._r8))
      w_bulk = water_loss - w_frac

      frac_loss = 0._r8
      IF (w_frac > mt_tiny) THEN
         frac_loss = material_evaporative_mt_loss(pool_mt, pool_water, w_frac, &
            temp_k, from_ice, evap_affine_fn, mt_tiny, r_max)
      ENDIF

      remaining_mt = pool_mt - frac_loss
      remaining_water = pool_water - w_frac
      bulk_loss = 0._r8
      IF (w_bulk > mt_tiny) THEN
         IF (remaining_water > mt_tiny) THEN
            bulk_loss = min(w_bulk * max(remaining_mt, 0._r8) / remaining_water, &
               max(remaining_mt, 0._r8))
         ELSE
            bulk_loss = max(remaining_mt, 0._r8)
         ENDIF
      ENDIF

      material_skin_limited_mt_loss = min(frac_loss + bulk_loss, max(pool_mt, 0._r8))
   END FUNCTION material_skin_limited_mt_loss

   !---------------------------------------------------------------
   ! Material leaving a finite, well-mixed pool that loses water_loss of its
   ! pool_water to the atmosphere, the evaporate ratio being affine in the
   ! pool ratio R (R_E = slope*R + offset). Mass balance dM = R_E dW with
   ! M = R*W gives  W dR/dW = (slope-1) R + offset,  integrated exactly over
   ! the loss (the desiccating-pool solution, Gonfiantini 1986):
   !
   !   x  = ln(W1/W0),  y = (slope-1) x
   !   R1 = R0 exp(y) + offset x (exp(y)-1)/y
   !
   ! For slope > 1 the pool ratio relaxes to R* = -offset/(slope-1), the
   ! ratio in balance with the ambient vapour. Near saturation slope is of
   ! order 1/(1-h), so the relaxation takes a tiny part of the loss. The
   ! former explicit substeps (10% of the pool each) were unstable there:
   ! a 5e-3 mm snowmelt pool losing 79% was left at delta18O -969 permil at
   ! h = 0.99 and +505 permil at h = 0.97 (exact: -25 in both), the two
   ! isotopes overshooting differently, and the residue reached the river
   ! as runoff (NEC 2004-05 with IsoGSM forcing).
   !
   ! r_max is a NUMERICAL/PHYSICAL safety ceiling on the residual pool
   ! ratio -- not part of the fractionation. Over repeated wet/dry cycles a
   ! pool could otherwise enrich to non-physical values; at or above r_max
   ! evaporation is non-fractionating (water-matched, exactly conservative)
   ! and below it the residual ratio is capped at r_max. Pass r_max <= 0 to
   ! disable; nonvolatile solutes are handled by callers as zero loss.
   !
   ! The result is signed: when the ambient vapour is heavier than the
   ! pool's equilibrium vapour the pool takes up heavy isotopes.
   !---------------------------------------------------------------
   real(r8) FUNCTION material_evaporative_mt_loss (pool_mt, pool_water, water_loss, &
      temp_k, from_ice, evap_affine_fn, mt_tiny, r_max)

      real(r8), intent(in) :: pool_mt
      real(r8), intent(in) :: pool_water
      real(r8), intent(in) :: water_loss
      real(r8), intent(in) :: temp_k
      logical,  intent(in) :: from_ice
      procedure(material_evap_affine_callback) :: evap_affine_fn
      real(r8), intent(in) :: mt_tiny
      real(r8), intent(in) :: r_max

      real(r8) :: source_ratio, slope, offset
      real(r8) :: residual_water, residual_ratio, min_loss_for_cap
      real(r8) :: x, y, growth, offset_factor

      material_evaporative_mt_loss = 0._r8
      IF (pool_water <= mt_tiny .or. water_loss <= mt_tiny) RETURN
      IF (water_loss >= pool_water * (1._r8 - 1.e-12_r8)) THEN
         ! ALL_POOL_EVAPORATION: when the water-side flux removes the whole
         ! finite pool, the integrated outgoing material is the whole material
         ! inventory.  Applying a last instantaneous fractionation ratio would
         ! leave non-zero material in a zero-water pool and break conservation.
         material_evaporative_mt_loss = max(pool_mt, 0._r8)
         RETURN
      ENDIF

      source_ratio = max(pool_mt, 0._r8) / pool_water
      CALL evap_affine_fn(temp_k, from_ice, slope, offset)
      residual_water = pool_water - water_loss

      IF ((slope == 1._r8 .and. offset == 0._r8) .or. &
          (r_max > 0._r8 .and. source_ratio >= r_max)) THEN
         ! Non-fractionating (no fractionation, sublimation, or at the
         ! ceiling): the vapour leaves at the pool ratio.
         material_evaporative_mt_loss = min(water_loss * source_ratio, max(pool_mt, 0._r8))
         IF (r_max > 0._r8 .and. source_ratio <= r_max) THEN
            min_loss_for_cap = max(pool_mt, 0._r8) - r_max * residual_water
            material_evaporative_mt_loss = max(material_evaporative_mt_loss, min_loss_for_cap)
         ENDIF
         RETURN
      ENDIF

      x = log(residual_water / pool_water)
      y = (slope - 1._r8) * x
      growth = exp(y)
      IF (abs(y) < 1.e-8_r8) THEN
         offset_factor = x * (1._r8 + 0.5_r8 * y)
      ELSE
         offset_factor = x * (growth - 1._r8) / y
      ENDIF
      residual_ratio = max(source_ratio * growth + offset * offset_factor, 0._r8)
      IF (r_max > 0._r8) residual_ratio = min(residual_ratio, r_max)

      material_evaporative_mt_loss = min(max(pool_mt, 0._r8) - residual_ratio * residual_water, &
         max(pool_mt, 0._r8))

   END FUNCTION material_evaporative_mt_loss

END MODULE MOD_MT_EvapLimit
#endif
