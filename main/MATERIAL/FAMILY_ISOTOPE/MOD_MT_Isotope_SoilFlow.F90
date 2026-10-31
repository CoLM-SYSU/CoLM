#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_Isotope_SoilFlow
!-----------------------------------------------------------------------
! Isotope rule for the inter-layer (qlayer) sweep of MOD_MT_SoilWater.
!
! The sweep first moves every face unclamped and keeps the result unless a
! layer ends negative by more than rounding of what passed through it;
! otherwise it falls back to a clamped top-down sweep. That fallback moves
! too little isotope out of a layer that is refilled from below later in
! the same sweep, so the isotope stays behind in that layer while the water
! goes on.
!
! A frozen layer with no liquid that only passes water upward (plant
! hydraulic root return into it, flow from below, flow out above) ends the
! host step with a liquid water of a few 1e-15 mm below zero: the host's
! own rounding of qlayer*deltim. The unclamped sweep leaves the isotope
! there below zero by that water times the layer ratio, which is larger
! than 1e-12 of what passed through. Taking this as a failed sweep moved
! 1.7e-4 mm worth of isotope into the wrong layer in a NEC December run
! and broke the fixed-signature check for the rest of the day.
!
! For isotopes the stock is therefore accepted when the layer's water is
! itself below zero by no more than mt_water_min_for_ratio and the
! isotope deficit is within twice that water deficit at the layer ratio.
! The sweep then zeroes it as an explicit numerical source, as it does
! for any rounding dust. Other materials keep the original rule.
!
! Wetland vapour uptake.  The last leaf temperature iteration corrects
! etr linearly (etr + etr_dtl*dtl), which on cold nights can leave a
! slightly negative transpiration.  The host wetland takes it as water
! gained from the air (wetwat + (gwat - etr)*deltim) and fevpa contains
! it, while the material wetland pool only removes positive
! transpiration.  The isotope pool then misses that water, its fluxes no
! longer match the host water fluxes times the fixed signature, and a
! NEC December run stopped on that check (1.0e-9 mm against fevpa).
! isotope_wetland_vapor_uptake returns the water the wetland gains this
! way; the wetland branch adds it to the pool as vapour deposition, like
! dew.  Other materials keep the original rule.
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_MT_Defs, only: mt_water_min_for_ratio
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: isotope_sweep_deficit_from_host_water
   PUBLIC :: isotope_wetland_vapor_uptake

CONTAINS

   pure logical FUNCTION isotope_sweep_deficit_from_host_water (stock, water, ratio)
      real(r8), intent(in) :: stock   ! isotope in the layer after the unclamped sweep
      real(r8), intent(in) :: water   ! layer liquid water after the sweep [mm]
      real(r8), intent(in) :: ratio   ! isotope/water ratio the layer passed on

      isotope_sweep_deficit_from_host_water = water < 0._r8 .and. &
         water >= -mt_water_min_for_ratio .and. &
         stock >= 2._r8 * water * max(ratio, 0._r8)
   END FUNCTION isotope_sweep_deficit_from_host_water

   pure real(r8) FUNCTION isotope_wetland_vapor_uptake (etr, deltim)
      real(r8), intent(in) :: etr      ! host transpiration [mm/s], + = to the air
      real(r8), intent(in) :: deltim   ! time step [s]

      isotope_wetland_vapor_uptake = max(-etr, 0._r8) * deltim
   END FUNCTION isotope_wetland_vapor_uptake

END MODULE MOD_MT_Isotope_SoilFlow
#endif
