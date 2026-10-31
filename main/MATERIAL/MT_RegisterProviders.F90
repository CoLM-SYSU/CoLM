#include <define.h>

#ifdef MaterialTransport
SUBROUTINE mt_register_providers ()
!-----------------------------------------------------------------------
!  Species with lifecycle callbacks. An external subroutine rather than a
!  module procedure because the providers USE MOD_MT_Lifecycle, which calls
!  this. A provider attaches only to a configured entry of the descriptor
!  table, so registering it with the species switched off does nothing.
!  Its own file, built with the main objects, because methane needs the
!  BGC soil modules there.
!-----------------------------------------------------------------------
#ifdef GridRiverLakeFlow
   USE MOD_MT_Sediment, only: register_sediment_mt_provider
#endif
#ifdef BGC
   USE MOD_MT_Methane, only: ch4_register_mt_provider
#endif
   IMPLICIT NONE
#ifdef GridRiverLakeFlow
   CALL register_sediment_mt_provider ()
#endif
#ifdef BGC
   CALL ch4_register_mt_provider ()
#endif
END SUBROUTINE mt_register_providers
#endif
