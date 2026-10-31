#include <define.h>

#if (defined MaterialTransport) && (defined BGC)
MODULE MOD_MT_Methane_Registry
!=======================================================================
! Methane's material index is populated by lifecycle provider registration.
!
!=======================================================================

   IMPLICIT NONE
   SAVE
   PRIVATE

   ! Sentinel value: material absent (or non-reactive) in the registry.
   integer, parameter :: METHANE_GAS_ABSENT = -1

   integer :: igas_ch4 = METHANE_GAS_ABSENT

   PUBLIC :: methane_is_active
   PUBLIC :: igas_ch4, METHANE_GAS_ABSENT

CONTAINS

   logical FUNCTION methane_is_active ()
      methane_is_active = (igas_ch4 > 0)
   END FUNCTION methane_is_active


END MODULE MOD_MT_Methane_Registry
#endif
