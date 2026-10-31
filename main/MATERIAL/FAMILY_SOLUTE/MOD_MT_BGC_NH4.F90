#include <define.h>

#if (defined MaterialTransport) && (defined BGC)
MODULE MOD_MT_BGC_NH4
   USE MOD_Precision
   USE MOD_MT_Namelist, only: DEF_USE_NH4_LEACHING
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: smin_nh4_leached_vr, smin_nh4_leached
   PUBLIC :: mt_nh4_leach_layer, mt_nh4_update_layer, mt_nh4_sum_layer
   PUBLIC :: mt_nh4_balance_outputs, mt_nh4_report_loss

   ! BGC retains the allocation/reset hooks at its existing worker lifecycle sites.
   real(r8), allocatable :: smin_nh4_leached_vr(:,:) ! gN m-3 s-1
   real(r8), allocatable :: smin_nh4_leached(:)       ! gN m-2 s-1

CONTAINS

   SUBROUTINE mt_nh4_leach_layer(j, i, sf, nh4, dz, water, tot_water, drain_tot, deltim)
      integer, intent(in) :: j, i
      real(r8), intent(in) :: sf, nh4, dz, water, tot_water, drain_tot, deltim
      real(r8) :: disn_conc

      IF (DEF_USE_NH4_LEACHING) THEN
         disn_conc = 0._r8
         IF (water > 0._r8) THEN
            disn_conc = (sf * max(nh4,0._r8) * dz) / water
         ENDIF
         IF (tot_water > 0._r8) THEN
            smin_nh4_leached_vr(j,i) = disn_conc * drain_tot * water / (tot_water * dz)
         ELSE
            smin_nh4_leached_vr(j,i) = 0._r8
         ENDIF
         smin_nh4_leached_vr(j,i) = min(smin_nh4_leached_vr(j,i), (sf * max(nh4,0._r8))/deltim)
         smin_nh4_leached_vr(j,i) = max(smin_nh4_leached_vr(j,i), 0._r8)
      ELSE
         smin_nh4_leached_vr(j,i) = 0._r8
      ENDIF
   END SUBROUTINE mt_nh4_leach_layer

   SUBROUTINE mt_nh4_update_layer(j, i, deltim, nh4)
      integer, intent(in) :: j, i
      real(r8), intent(in) :: deltim
      real(r8), intent(inout) :: nh4

      IF (DEF_USE_NH4_LEACHING) THEN
         nh4 = max(nh4 - smin_nh4_leached_vr(j,i) * deltim, 0._r8)
      ENDIF
   END SUBROUTINE mt_nh4_update_layer

   SUBROUTINE mt_nh4_sum_layer(j, i, dz)
      integer, intent(in) :: j, i
      real(r8), intent(in) :: dz

      smin_nh4_leached(i) = smin_nh4_leached(i) + smin_nh4_leached_vr(j,i) * dz
   END SUBROUTINE mt_nh4_sum_layer

   SUBROUTINE mt_nh4_balance_outputs(i, col_noutputs, nit, no3_leached, no3_runoff)
      integer, intent(in) :: i
      real(r8), intent(inout) :: col_noutputs
      real(r8), intent(in) :: nit, no3_leached, no3_runoff

      ! Retain NH4 between nitrification and nitrate in the host's addition order.
      col_noutputs = col_noutputs + nit + smin_nh4_leached(i) + no3_leached + no3_runoff
   END SUBROUTINE mt_nh4_balance_outputs

   SUBROUTINE mt_nh4_report_loss(i, deltim)
      integer, intent(in) :: i
      real(r8), intent(in) :: deltim

      write(*,*)'NH4 leached mass         = ',smin_nh4_leached(i)*deltim
   END SUBROUTINE mt_nh4_report_loss

END MODULE MOD_MT_BGC_NH4
#endif
