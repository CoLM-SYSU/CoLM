#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_SnowHost
   USE MOD_Precision
   USE MOD_MT_Defs, only: nmts
   USE MOD_MT_Vars, only: mt_wliq_soisno, mt_wice_soisno, mt_solid_soisno, mt_scv
   USE MOD_SnowLayersCombineDivide, only: snowlayerscombine, snowlayerscombine_snicar, &
      snowlayersdivide, snowlayersdivide_snicar
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: mt_land_snow_combine, mt_land_snow_divide, mt_land_snow_normalize
CONTAINS
   SUBROUTINE mt_land_snow_normalize(ipatch, use_snicar, lb, snl, &
      z_soisno, dz_soisno, zi_soisno, wliq_soisno, wice_soisno, t_soisno, &
      scv, snowdp, &
      mss_bcpho,mss_bcphi,mss_ocpho,mss_ocphi,mss_dst1,mss_dst2,mss_dst3,mss_dst4)
      integer, intent(in) :: ipatch, lb
      logical, intent(in) :: use_snicar
      integer, intent(inout) :: snl
      real(r8), intent(inout) :: z_soisno(lb:1), dz_soisno(lb:1), zi_soisno(lb-1:1)
      real(r8), intent(inout) :: wliq_soisno(lb:1), wice_soisno(lb:1), t_soisno(lb:1)
      real(r8), intent(inout) :: scv, snowdp
      real(r8), intent(inout) :: mss_bcpho(lb:0), mss_bcphi(lb:0), mss_ocpho(lb:0), mss_ocphi(lb:0), &
         mss_dst1(lb:0), mss_dst2(lb:0), mss_dst3(lb:0), mss_dst4(lb:0)

      CALL mt_land_snow_combine(ipatch, use_snicar, lb, snl, &
         z_soisno, dz_soisno, zi_soisno, wliq_soisno, wice_soisno, t_soisno, scv, snowdp, &
         mss_bcpho,mss_bcphi,mss_ocpho,mss_ocphi,mss_dst1,mss_dst2,mss_dst3,mss_dst4)
      IF (snl < 0) CALL mt_land_snow_divide(ipatch, use_snicar, lb, snl, &
         z_soisno(lb:0), dz_soisno(lb:0), zi_soisno(lb-1:0), &
         wliq_soisno(lb:0), wice_soisno(lb:0), t_soisno(lb:0), &
         mss_bcpho,mss_bcphi,mss_ocpho,mss_ocphi,mss_dst1,mss_dst2,mss_dst3,mss_dst4)
   END SUBROUTINE mt_land_snow_normalize

   SUBROUTINE mt_land_snow_combine(ipatch, use_snicar, lb, snl, &
      z_soisno, dz_soisno, zi_soisno, wliq_soisno, wice_soisno, t_soisno, &
      scv, snowdp, &
      mss_bcpho,mss_bcphi,mss_ocpho,mss_ocphi,mss_dst1,mss_dst2,mss_dst3,mss_dst4)
      integer, intent(in) :: ipatch, lb
      logical, intent(in) :: use_snicar
      integer, intent(inout) :: snl
      real(r8), intent(inout) :: z_soisno(lb:1), dz_soisno(lb:1), zi_soisno(lb-1:1)
      real(r8), intent(inout) :: wliq_soisno(lb:1), wice_soisno(lb:1), t_soisno(lb:1)
      real(r8), intent(inout) :: scv, snowdp
      real(r8), intent(inout) :: mss_bcpho(lb:0), mss_bcphi(lb:0), mss_ocpho(lb:0), mss_ocphi(lb:0), &
         mss_dst1(lb:0), mss_dst2(lb:0), mss_dst3(lb:0), mss_dst4(lb:0)

      ! Materials are allocated only after registration; keep the absent-option path.
      IF (use_snicar) THEN
         IF (nmts > 0) THEN
            CALL snowlayerscombine_snicar (lb,snl,&
            z_soisno(lb:1),dz_soisno(lb:1),zi_soisno(lb-1:1),&
            wliq_soisno(lb:1),wice_soisno(lb:1),t_soisno(lb:1),scv,snowdp,&
            mss_bcpho(lb:0), mss_bcphi(lb:0), mss_ocpho(lb:0), mss_ocphi(lb:0),&
            mss_dst1(lb:0), mss_dst2(lb:0), mss_dst3(lb:0), mss_dst4(lb:0), &
            mt_wliq = mt_wliq_soisno(:, lb:1, ipatch), &
            mt_wice = mt_wice_soisno(:, lb:1, ipatch), &
            mt_solid = mt_solid_soisno(:, lb:1, ipatch), &
            mt_scv  = mt_scv(:, ipatch))
         ELSE
            CALL snowlayerscombine_snicar (lb,snl,&
            z_soisno(lb:1),dz_soisno(lb:1),zi_soisno(lb-1:1),&
            wliq_soisno(lb:1),wice_soisno(lb:1),t_soisno(lb:1),scv,snowdp,&
            mss_bcpho(lb:0), mss_bcphi(lb:0), mss_ocpho(lb:0), mss_ocphi(lb:0),&
            mss_dst1(lb:0), mss_dst2(lb:0), mss_dst3(lb:0), mss_dst4(lb:0) )
         ENDIF
      ELSE
         IF (nmts > 0) THEN
            CALL snowlayerscombine (lb,snl,&
            z_soisno(lb:1),dz_soisno(lb:1),zi_soisno(lb-1:1),&
            wliq_soisno(lb:1),wice_soisno(lb:1),t_soisno(lb:1),scv,snowdp, &
            mt_wliq = mt_wliq_soisno(:, lb:1, ipatch), &
            mt_wice = mt_wice_soisno(:, lb:1, ipatch), &
            mt_solid = mt_solid_soisno(:, lb:1, ipatch), &
            mt_scv  = mt_scv(:, ipatch))
         ELSE
            CALL snowlayerscombine (lb,snl,&
            z_soisno(lb:1),dz_soisno(lb:1),zi_soisno(lb-1:1),&
            wliq_soisno(lb:1),wice_soisno(lb:1),t_soisno(lb:1),scv,snowdp)
         ENDIF
      ENDIF
   END SUBROUTINE mt_land_snow_combine

   SUBROUTINE mt_land_snow_divide(ipatch, use_snicar, lb, snl, &
      z_soisno, dz_soisno, zi_soisno, wliq_soisno, wice_soisno, t_soisno, &
      mss_bcpho,mss_bcphi,mss_ocpho,mss_ocphi,mss_dst1,mss_dst2,mss_dst3,mss_dst4)
      integer, intent(in) :: ipatch, lb
      logical, intent(in) :: use_snicar
      integer, intent(inout) :: snl
      real(r8), intent(inout) :: z_soisno(lb:0), dz_soisno(lb:0), zi_soisno(lb-1:0)
      real(r8), intent(inout) :: wliq_soisno(lb:0), wice_soisno(lb:0), t_soisno(lb:0)
      real(r8), intent(inout) :: mss_bcpho(lb:0), mss_bcphi(lb:0), mss_ocpho(lb:0), mss_ocphi(lb:0), &
         mss_dst1(lb:0), mss_dst2(lb:0), mss_dst3(lb:0), mss_dst4(lb:0)

      ! Materials are allocated only after registration; keep the absent-option path.
      IF (use_snicar) THEN
         IF (nmts > 0) THEN
            CALL snowlayersdivide_snicar (lb,snl,&
            z_soisno(lb:0),dz_soisno(lb:0),zi_soisno(lb-1:0),&
            wliq_soisno(lb:0),wice_soisno(lb:0),t_soisno(lb:0),&
            mss_bcpho(lb:0),mss_bcphi(lb:0),mss_ocpho(lb:0),mss_ocphi(lb:0),&
            mss_dst1(lb:0),mss_dst2(lb:0),mss_dst3(lb:0),mss_dst4(lb:0), &
            mt_wliq = mt_wliq_soisno(:, lb:0, ipatch), &
            mt_wice = mt_wice_soisno(:, lb:0, ipatch), &
            mt_solid = mt_solid_soisno(:, lb:0, ipatch))
         ELSE
            CALL snowlayersdivide_snicar (lb,snl,&
            z_soisno(lb:0),dz_soisno(lb:0),zi_soisno(lb-1:0),&
            wliq_soisno(lb:0),wice_soisno(lb:0),t_soisno(lb:0),&
            mss_bcpho(lb:0),mss_bcphi(lb:0),mss_ocpho(lb:0),mss_ocphi(lb:0),&
            mss_dst1(lb:0),mss_dst2(lb:0),mss_dst3(lb:0),mss_dst4(lb:0) )
         ENDIF
      ELSE
         IF (nmts > 0) THEN
            CALL snowlayersdivide (lb,snl,&
            z_soisno(lb:0),dz_soisno(lb:0),zi_soisno(lb-1:0),&
            wliq_soisno(lb:0),wice_soisno(lb:0),t_soisno(lb:0), &
            mt_wliq = mt_wliq_soisno(:, lb:0, ipatch), &
            mt_wice = mt_wice_soisno(:, lb:0, ipatch), &
            mt_solid = mt_solid_soisno(:, lb:0, ipatch))
         ELSE
            CALL snowlayersdivide (lb,snl,&
            z_soisno(lb:0),dz_soisno(lb:0),zi_soisno(lb-1:0),&
            wliq_soisno(lb:0),wice_soisno(lb:0),t_soisno(lb:0))
         ENDIF
      ENDIF
   END SUBROUTINE mt_land_snow_divide

END MODULE MOD_MT_SnowHost
#endif
