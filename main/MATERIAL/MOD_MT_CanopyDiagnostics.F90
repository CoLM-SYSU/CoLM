#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_CanopyDiagnostics
   USE MOD_Precision
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: mt_canopy_drainage, mt_canopy_drain

   ! Net interception cannot recover gross interception without actual drainage.
   TYPE mt_canopy_drainage
      real(r8) :: rain = 0._r8, snow = 0._r8
   CONTAINS
      PROCEDURE :: clear => clear_drainage
      PROCEDURE :: capture => capture_drainage
      PROCEDURE :: accumulate => accumulate_drainage
   END TYPE mt_canopy_drainage

   ! The existing last-interception-call ledger, consumed by the land handoff.
   TYPE(mt_canopy_drainage), SAVE :: mt_canopy_drain

CONTAINS

   SUBROUTINE clear_drainage(this)
      class(mt_canopy_drainage), intent(inout) :: this
      this%rain = 0._r8
      this%snow = 0._r8
   END SUBROUTINE clear_drainage

   SUBROUTINE capture_drainage(this, deltim, tex_rain, tex_snow)
      class(mt_canopy_drainage), intent(inout) :: this
      real(r8), intent(in) :: deltim, tex_rain, tex_snow
      IF (deltim > 0._r8) THEN
         this%rain = max(tex_rain, 0._r8) / deltim
         this%snow = max(tex_snow, 0._r8) / deltim
      ELSE
         CALL this%clear()
      ENDIF
   END SUBROUTINE capture_drainage

   SUBROUTINE accumulate_drainage(this, drainage, fraction)
      class(mt_canopy_drainage), intent(inout) :: this
      type(mt_canopy_drainage), intent(in) :: drainage
      real(r8), intent(in) :: fraction
      this%rain = this%rain + drainage%rain * fraction
      this%snow = this%snow + drainage%snow * fraction
   END SUBROUTINE accumulate_drainage

END MODULE MOD_MT_CanopyDiagnostics
#endif
