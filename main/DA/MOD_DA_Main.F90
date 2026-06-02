#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Main
!-----------------------------------------------------------------------------
! DESCRIPTION:
!     Main procedures for data assimilation
!
! AUTHOR:
!     Lu Li, 12/2024
!     Zhilong Fan, Lu Li, 03/2024
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Spmd_task
   USE MOD_Namelist
   USE MOD_LandPatch
   USE MOD_DA_Assim_SM
   USE MOD_DA_Ens
   IMPLICIT NONE
   SAVE

   logical :: use_assim_sm = .false.

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_DA ()

!-----------------------------------------------------------------------------
   IMPLICIT NONE

   integer :: i

!-----------------------------------------------------------------------------

      use_assim_sm = .false.
      DO i = 1, size(DEF_DA_OBS_TARGET)
         IF (trim(DEF_DA_OBS_TARGET(i)) == 'SM') THEN
            use_assim_sm = .true.
            EXIT
         ENDIF
      ENDDO
      IF (use_assim_sm) CALL init_Assim_SM ()

   END SUBROUTINE init_DA

!-----------------------------------------------------------------------------

   SUBROUTINE run_DA (idate, deltim, dolai, doalb, dosst, oro)

!-----------------------------------------------------------------------------
   IMPLICIT NONE

   integer,  intent(in)    :: idate(3)
   real(r8), intent(in)    :: deltim
   logical,  intent(in)    :: dolai    
   logical,  intent(in)    :: doalb    
   logical,  intent(in)    :: dosst    
   real(r8), intent(inout) :: oro(numpatch)  

   integer :: i, n

!-----------------------------------------------------------------------------
      ! generate ensemble members
      IF (p_is_worker) THEN
         CALL model_to_ens_member(0)
         CALL ensemble (deltim)
         DO i = 1, DEF_DA_ENS_NUM
            CALL ens_member_to_model(i)
            CALL CoLMDRIVER (idate, deltim, dolai, doalb, dosst, oro)
            CALL model_to_ens_member(i)
         ENDDO
         CALL ens_member_to_model(0)
      ENDIF

      ! data assimilation
      IF (p_is_master) THEN
         print *, '[CoLM-DA] Start surface soil moisture data assimilation.'

         n = 0
         write(*,'(A)') '[CoLM-DA] Soil moisture assimilation observation sources:'
         write(*,'(A)') '[CoLM-DA] +------+----------+------------------+----------+------------+'
         write(*,'(A)') '[CoLM-DA] | No.  | Source   | Sensor           | Var      | Freq(GHz)  |'
         write(*,'(A)') '[CoLM-DA] +------+----------+------------------+----------+------------+'

         DO i = 1, size(DEF_DA_OBS_TARGET)
            IF (trim(DEF_DA_OBS_TARGET(i)) == 'SM') THEN
               n = n + 1
               write(*,'(A,I4,A,A8,A,A16,A,A8,A,F10.3,A)') &
                  '[CoLM-DA] | ', n, ' | ', &
                  DEF_DA_OBS_SOURCE(i), ' | ', &
                  DEF_DA_OBS_SENSOR(i), ' | ', &
                  DEF_DA_OBS_VAR(i),    ' | ', &
                  DEF_DA_OBS_FGHZ(i),         ' |'
            ENDIF
         ENDDO

         IF (n == 0) THEN
            write(*,'(A)') '[CoLM-DA] | none |          |                  |          |            |'
         ENDIF

         write(*,'(A)') '[CoLM-DA] +------+----------+------------------+----------+------------+'
      ENDIF


      CALL mpi_barrier (p_comm_glb, p_err)
      IF (use_assim_sm) THEN
         CALL run_Assim_SM (idate, deltim)
      ENDIF

   END SUBROUTINE run_DA

!-----------------------------------------------------------------------------

   SUBROUTINE end_DA ()

!-----------------------------------------------------------------------------
   IMPLICIT NONE

!-----------------------------------------------------------------------------

      IF (use_assim_sm) CALL end_Assim_SM ()
      use_assim_sm = .false.

   END SUBROUTINE end_DA

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Main
#endif
