#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_Providers
!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!  The lists of compiled material species code:
!
!  - isotope fractionation physics (ensure_isotope_physics_registered);
!  - species that own state or phase callbacks (mt_register_providers,
!    in MT_RegisterProviders.F90, built after the providers it lists).
!
!  Parameter-only isotopes and solutes are not listed here. Add a line only
!  when the species code exists.
!-----------------------------------------------------------------------

   USE MOD_MT_Isotope_O18, only: register_o18_isotope_physics
   USE MOD_MT_Isotope_HDO, only: register_hdo_isotope_physics
   IMPLICIT NONE
   SAVE
   PRIVATE

   logical :: isotope_physics_registered = .false.

   PUBLIC :: ensure_isotope_physics_registered
   PUBLIC :: register_all_isotope_physics

CONTAINS

   SUBROUTINE ensure_isotope_physics_registered ()
      IF (isotope_physics_registered) RETURN
      CALL register_all_isotope_physics ()
      isotope_physics_registered = .true.
   END SUBROUTINE ensure_isotope_physics_registered

   SUBROUTINE register_all_isotope_physics ()
      CALL register_o18_isotope_physics ()
      CALL register_hdo_isotope_physics ()
   END SUBROUTINE register_all_isotope_physics

END MODULE MOD_MT_Providers
#endif
