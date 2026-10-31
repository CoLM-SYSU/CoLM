#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_ThermalDiagnostics
   USE MOD_Precision
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: mt_canopy_diagnostics, mt_thermal_diagnostics
   PUBLIC :: mt_thermal_reset, mt_phase_reset, mt_phase_record, mt_canopy_average

   TYPE mt_canopy_diagnostics
      real(r8) :: melt = 0._r8, freeze = 0._r8 ! mm per step
      real(r8) :: raw = 0._r8                  ! moisture resistance [s/m]
   END TYPE mt_canopy_diagnostics

   ! Per-call observations; host energy and water solutions never read these.
   TYPE mt_thermal_diagnostics
      type(mt_canopy_diagnostics) :: canopy
      real(r8), allocatable :: soil_thaw(:), soil_freeze(:)
   END TYPE mt_thermal_diagnostics

CONTAINS

   SUBROUTINE mt_thermal_reset(diag, lb, nl_soil)
      type(mt_thermal_diagnostics), intent(inout), optional :: diag
      integer, intent(in) :: lb, nl_soil

      IF (.not. present(diag)) RETURN
      diag%canopy = mt_canopy_diagnostics()
      CALL mt_phase_reset(diag, lb, nl_soil)
   END SUBROUTINE mt_thermal_reset

   SUBROUTINE mt_phase_reset(diag, lb, nl_soil)
      type(mt_thermal_diagnostics), intent(inout), optional :: diag
      integer, intent(in) :: lb, nl_soil

      IF (.not. present(diag)) RETURN
      IF (allocated(diag%soil_thaw)) THEN
         IF (lbound(diag%soil_thaw,1) /= lb .or. ubound(diag%soil_thaw,1) /= nl_soil) &
            deallocate(diag%soil_thaw, diag%soil_freeze)
      ENDIF
      IF (.not. allocated(diag%soil_thaw)) &
         allocate(diag%soil_thaw(lb:nl_soil), diag%soil_freeze(lb:nl_soil))
      diag%soil_thaw = 0._r8
      diag%soil_freeze = 0._r8
   END SUBROUTINE mt_phase_reset

   SUBROUTINE mt_phase_record(diag, layer, ice_before, ice_after)
      type(mt_thermal_diagnostics), intent(inout), optional :: diag
      integer, intent(in) :: layer
      real(r8), intent(in) :: ice_before, ice_after

      IF (.not. present(diag)) RETURN
      diag%soil_thaw(layer) = max(ice_before - ice_after, 0._r8)
      diag%soil_freeze(layer) = max(ice_after - ice_before, 0._r8)
   END SUBROUTINE mt_phase_record

   SUBROUTINE mt_canopy_average(ps, pe, canopy_p, pftfrac, canopy)
      integer, intent(in) :: ps, pe
      type(mt_canopy_diagnostics), intent(in) :: canopy_p(ps:pe)
      real(r8), intent(in) :: pftfrac(ps:pe)
      type(mt_canopy_diagnostics), intent(out) :: canopy

      canopy%melt = sum(canopy_p(ps:pe)%melt * pftfrac(ps:pe))
      canopy%freeze = sum(canopy_p(ps:pe)%freeze * pftfrac(ps:pe))
      canopy%raw = sum(canopy_p(ps:pe)%raw * pftfrac(ps:pe))
   END SUBROUTINE mt_canopy_average
END MODULE MOD_MT_ThermalDiagnostics
#endif
