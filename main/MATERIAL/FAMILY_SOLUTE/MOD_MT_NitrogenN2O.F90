#include <define.h>
#ifdef MaterialTransport

MODULE MOD_MT_NitrogenN2O

!-----------------------------------------------------------------------
! Indirect N2O conversion at the land-to-river entry for the nitrogen
! solutes NH4 and NO3.  This is the FAMILY_SOLUTE-specific policy; the
! generic river code (GENERIC_WATER/MOD_MT_RiverLake) only asks this
! module for the factor and applies it to the incoming mass.
!
! The factor comes from the solute's own parameter file
! (DEF_SOLUTE%n2o_ef5).  A non-negative DEF_MT_N2O_EF5_NH4/NO3 in
! &nl_colm_material overrides it, which makes sensitivity runs possible
! without editing the parameter files.  Any other solute gets 0 and its
! transport is untouched, so configurations without NH4/NO3 (or with the
! factors left at 0) remain bit-identical.
!
! This is the inventory arithmetic (N_leached x EF5), not a mechanistic
! water-body model: it does not represent where in the river network the
! emission occurs or how long the water resides, and the factor must not
! be applied again downstream.  Reference values and sources are
! documented at the DEF_MT_N2O_EF5_* declarations in MOD_MT_Namelist and
! in run/standard_nh4_parameter.nml / run/standard_no3_parameter.nml
! (IPCC 2019 EF5g/EF5r/EF5e; Tian et al. 2019; Maavara et al. 2019).
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_MT_Defs, only: nmts, materials, material_lower
   USE MOD_MT_Namelist, only: DEF_MT_N2O_EF5_NH4, DEF_MT_N2O_EF5_NO3
   IMPLICIT NONE
   PUBLIC :: mt_n2o_ef5_for

CONTAINS

   real(r8) FUNCTION mt_n2o_ef5_for (imt)

      IMPLICIT NONE
      integer, intent(in) :: imt
      character(len=32) :: lname

      mt_n2o_ef5_for = 0._r8
      IF (imt < 1 .or. imt > nmts) RETURN

      lname = material_lower(trim(materials(imt)%name))
      SELECT CASE (lname)
      CASE ('nh4')
         mt_n2o_ef5_for = materials(imt)%n2o_ef5
         ! Any value >= 0 is a deliberate override, so 0 switches the
         ! conversion off; the -1 default means "use the parameter file".
         IF (DEF_MT_N2O_EF5_NH4 >= 0._r8) mt_n2o_ef5_for = DEF_MT_N2O_EF5_NH4
      CASE ('no3')
         mt_n2o_ef5_for = materials(imt)%n2o_ef5
         IF (DEF_MT_N2O_EF5_NO3 >= 0._r8) mt_n2o_ef5_for = DEF_MT_N2O_EF5_NO3
      END SELECT

   END FUNCTION mt_n2o_ef5_for

END MODULE MOD_MT_NitrogenN2O
#endif
