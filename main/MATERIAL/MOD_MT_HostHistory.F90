#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_HostHistory
   USE MOD_Precision
   USE MOD_Vars_Global, only: spval
   USE MOD_MT_Namelist, only: mt_any_enabled
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: frcsat_hist_valid, nac_frcsat
   PUBLIC :: mt_allocate_frcsat_valid, mt_deallocate_frcsat_valid
   PUBLIC :: mt_allocate_frcsat_count, mt_deallocate_frcsat_count, mt_flush_frcsat_count
   PUBLIC :: mt_accumulate_frcsat

   logical, allocatable :: frcsat_hist_valid(:)
   real(r8), allocatable :: nac_frcsat(:)

CONTAINS

   SUBROUTINE mt_allocate_frcsat_valid(numpatch)
      integer, intent(in) :: numpatch
      allocate(frcsat_hist_valid(numpatch)); frcsat_hist_valid(:) = .true.
   END SUBROUTINE mt_allocate_frcsat_valid

   SUBROUTINE mt_deallocate_frcsat_valid()
      deallocate(frcsat_hist_valid)
   END SUBROUTINE mt_deallocate_frcsat_valid

   SUBROUTINE mt_allocate_frcsat_count(numpatch)
      integer, intent(in) :: numpatch
      IF (mt_any_enabled()) allocate(nac_frcsat(numpatch))
   END SUBROUTINE mt_allocate_frcsat_count

   SUBROUTINE mt_deallocate_frcsat_count()
      IF (allocated(nac_frcsat)) deallocate(nac_frcsat)
   END SUBROUTINE mt_deallocate_frcsat_count

   SUBROUTINE mt_flush_frcsat_count()
      IF (allocated(nac_frcsat)) nac_frcsat(:) = 0._r8
   END SUBROUTINE mt_flush_frcsat_count

   SUBROUTINE mt_accumulate_frcsat(frcsat, a_frcsat, accumulate)
      real(r8), intent(in) :: frcsat(:)
      real(r8), intent(inout) :: a_frcsat(:)
      integer :: i
      ! Callback keeps the host's accumulator implementation without a module cycle.
      interface
         SUBROUTINE accumulate(var, s, filter)
            USE MOD_Precision
            real(r8), intent(in) :: var(:)
            real(r8), intent(inout) :: s(:)
            logical, intent(in), optional :: filter(:)
         END SUBROUTINE accumulate
      end interface

      IF (mt_any_enabled()) THEN
         CALL accumulate(frcsat, a_frcsat, frcsat_hist_valid)
         DO i = 1, size(frcsat)
            IF (.not. frcsat_hist_valid(i)) CYCLE
            ! A legacy partial window has no recoverable valid-sample count.
            IF (nac_frcsat(i) < 0._r8) CYCLE
            IF (frcsat(i) /= spval) nac_frcsat(i) = nac_frcsat(i) + 1._r8
         ENDDO
      ELSE
         CALL accumulate(frcsat, a_frcsat)
      ENDIF
   END SUBROUTINE mt_accumulate_frcsat


END MODULE MOD_MT_HostHistory
#endif
