#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_Isotope_O18

   USE MOD_Precision
   USE MOD_MT_Namelist, only: DEF_ISO_KINETIC_SCHEME
   USE MOD_MT_Isotope_Registry, only: register_isotope_physics

   IMPLICIT NONE
   SAVE
   PRIVATE

   ! VSMOW reference ratio
   real(r8), parameter :: Rsmow_18O = 2.0052e-3_r8
   real(r8), parameter :: o18_diffusivity_ratio_cappa = 1.03189_r8
   real(r8), parameter :: o18_diffusivity_ratio_merlivat = 1.0285_r8

   PUBLIC :: register_o18_isotope_physics

CONTAINS

   SUBROUTINE register_o18_isotope_physics ()
      CALL register_isotope_physics(name='O18', name_patterns='18o,o18', &
         mj79_relative_factor=1.0_r8, &
         ref_ratio_hint=Rsmow_18O, &
         default_soil_init_varname='soilwat_O18', &
         alpha_liq_vap_fn=o18_alpha_liq_vap, &
         alpha_ice_vap_fn=o18_alpha_ice_vap, &
         diffusivity_ratio_air_fn=o18_diffusivity_ratio, &
         leaf_liquid_diffusivity_fn=o18_leaf_liquid_diffusivity)
   END SUBROUTINE register_o18_isotope_physics

   real(r8) FUNCTION o18_alpha_liq_vap (temp_k)
      real(r8), intent(in) :: temp_k
      real(r8) :: tk

      tk = max(temp_k, 150._r8)
      o18_alpha_liq_vap = exp(1137._r8 / (tk * tk) - 0.4156_r8 / tk - 0.0020667_r8)
   END FUNCTION o18_alpha_liq_vap

   real(r8) FUNCTION o18_alpha_ice_vap (temp_k)
      real(r8), intent(in) :: temp_k
      real(r8) :: tk

      tk = max(temp_k, 150._r8)
      o18_alpha_ice_vap = exp(11.839_r8 / tk - 0.028224_r8)
   END FUNCTION o18_alpha_ice_vap

   real(r8) FUNCTION o18_diffusivity_ratio ()
      IF (trim(DEF_ISO_KINETIC_SCHEME) == 'MERLIVAT1978') THEN
         o18_diffusivity_ratio = o18_diffusivity_ratio_merlivat
      ELSE
         o18_diffusivity_ratio = o18_diffusivity_ratio_cappa
      ENDIF
   END FUNCTION o18_diffusivity_ratio

   real(r8) FUNCTION o18_leaf_liquid_diffusivity (temp_k)
      real(r8), intent(in) :: temp_k
      real(r8) :: tk

      tk = max(temp_k, 150._r8)
      o18_leaf_liquid_diffusivity = 119.e-9_r8 * exp(-637._r8 / max(tk - 137._r8, 1._r8))
   END FUNCTION o18_leaf_liquid_diffusivity

END MODULE MOD_MT_Isotope_O18
#endif
