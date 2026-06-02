#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Assim_SM
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Main soil-moisture assimilation driver for multiple passive-microwave
!    observation streams.
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_DA_Obs_SM
   USE MOD_DA_EnKF
   USE MOD_DA_Vars_TimeVariables
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_Global, ONLY: pi, nl_soil, dz_soi
   USE MOD_LandPatch
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Const_Physical, ONLY: denice, denh2o, tfrz
   IMPLICIT NONE

   type(SM_type) :: sm

   PUBLIC :: init_Assim_SM
   PUBLIC :: run_Assim_SM
   PUBLIC :: end_Assim_SM

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_Assim_SM ()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      integer :: is, i, nsensor
      type(SM_config_type), allocatable :: configs(:)

!-----------------------------------------------------------------------------
      nsensor = 0
      DO i = 1, size(DEF_DA_OBS_TARGET)
         IF (trim(DEF_DA_OBS_TARGET(i)) == 'SM') THEN
            nsensor = nsensor + 1
         ENDIF
      ENDDO
      IF (nsensor == 0) RETURN

      allocate(configs(nsensor))
      is = 0
      DO i = 1, size(DEF_DA_OBS_TARGET)
         IF (trim(DEF_DA_OBS_TARGET(i)) == 'SM') THEN
            is = is + 1
            configs(is)%source_name = DEF_DA_OBS_SOURCE(i)
            configs(is)%sensor_name = DEF_DA_OBS_SENSOR(i)
            configs(is)%var_name    = DEF_DA_OBS_VAR(i)
            configs(is)%fghz        = DEF_DA_OBS_FGHZ(i)
         ENDIF
      ENDDO

      CALL sm%init(landpatch, configs)

      deallocate(configs)

   END SUBROUTINE init_Assim_SM

!-----------------------------------------------------------------------------

   SUBROUTINE run_Assim_SM (idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      integer,  intent(in) :: idate(3)
      real(r8), intent(in) :: deltim

!------------------------ Local Variables ------------------------------------
      integer :: np
      integer :: iens
      integer :: il
      integer :: nobs_p
      integer :: npatch_assim
      integer :: nobs_p_max
      real(r8) :: eff_porsl
      real(r8), allocatable :: y_p(:)
      real(r8), allocatable :: hx_p(:,:)
      real(r8), allocatable :: r_p(:)
      real(r8), allocatable :: trans(:,:)
      real(r8), allocatable :: wliq_ens_f(:,:)
      real(r8), allocatable :: wliq_ens_a(:,:)

!-----------------------------------------------------------------------------
      ! calculate y and H(x) at y locations from different sources
      CALL sm%calcg(idate, deltim)
      CALL sm%calcp()
      CALL sm%save(.false.)

      IF (p_is_worker .and. numpatch > 0) THEN
         wliq_soisno_ol(:,:) = wliq_soisno_ens(:,0,:)
         wice_soisno_ol(:,:) = wice_soisno_ens(:,0,:)
         t_soisno_ol   (:,:) = t_soisno_ens   (:,0,:)
         wliq_soisno_f (:,:) = sum(wliq_soisno_ens(:,1:DEF_DA_ENS_NUM,:), dim=2) / DEF_DA_ENS_NUM
         wice_soisno_f (:,:) = sum(wice_soisno_ens(:,1:DEF_DA_ENS_NUM,:), dim=2) / DEF_DA_ENS_NUM
         t_soisno_f    (:,:) = sum(t_soisno_ens   (:,1:DEF_DA_ENS_NUM,:), dim=2) / DEF_DA_ENS_NUM
      ENDIF

      ! loop over patches and perform assimilation for each patch
      IF (p_is_worker) THEN
         npatch_assim = 0
         nobs_p_max = 0
         DO np = 1, numpatch
            IF (patchtype(np) >= 3) CYCLE

            ! prepare localized y, R, and H(x) for the target patch
            CALL sm%prepare(np, nobs_p, y_p, hx_p, r_p)
            IF (nobs_p == 0) CYCLE
            npatch_assim = npatch_assim + 1
            nobs_p_max = max(nobs_p_max, nobs_p)

            ! perform data assimilation
            allocate(trans(DEF_DA_ENS_NUM, DEF_DA_ENS_NUM))
            CALL letkf(DEF_DA_ENS_NUM, nobs_p, hx_p, y_p, r_p, sm%infl, trans)

            allocate(wliq_ens_f(nl_soil, DEF_DA_ENS_NUM))
            allocate(wliq_ens_a(nl_soil, DEF_DA_ENS_NUM))
            DO iens = 1, DEF_DA_ENS_NUM
               wliq_ens_f(:,iens) = wliq_soisno_ens(1:,iens,np)
            ENDDO
            CALL dgemm('N', 'N', nl_soil, DEF_DA_ENS_NUM, DEF_DA_ENS_NUM, &
                        1.0_r8, wliq_ens_f, nl_soil, trans, DEF_DA_ENS_NUM, &
                        0.0_r8, wliq_ens_a, nl_soil)
            wliq_ens_a = max(0.0_r8, wliq_ens_a)
            DO iens = 1, DEF_DA_ENS_NUM
               wliq_soisno_ens(1:2,iens,np) = wliq_ens_a(1:2,iens)
            ENDDO

            ! postprocess soil moisture to ensure physical consistency
            DO il = 1, nl_soil
               DO iens = 1, DEF_DA_ENS_NUM
                  wliq_soisno_ens(il,iens,np) = max(0.0_r8, wliq_soisno_ens(il,iens,np))
                  wice_soisno_ens(il,iens,np) = max(0.0_r8, wice_soisno_ens(il,iens,np))

                  IF (wliq_soisno_ens(il,iens,np) == 0.0_r8 .and. &
                      wice_soisno_ens(il,iens,np) == 0.0_r8) THEN
                     IF (t_soisno_ens(il,iens,np) - tfrz < -5.0_r8) THEN
                        wice_soisno_ens(il,iens,np) = 1.0e-10_r8
                     ELSE
                        wliq_soisno_ens(il,iens,np) = 1.0e-10_r8
                     ENDIF
                  ENDIF

                  wliq_soisno_ens(il,iens,np) = min( &
                     porsl(il,np)*(dz_soi(il)*denh2o), &
                     wliq_soisno_ens(il,iens,np))
                  eff_porsl = max(0.0_r8, porsl(il,np) - &
                     wliq_soisno_ens(il,iens,np)/(dz_soi(il)*denh2o))
                  wice_soisno_ens(il,iens,np) = min( &
                     eff_porsl*(dz_soi(il)*denice), &
                     wice_soisno_ens(il,iens,np))
               ENDDO
            ENDDO

            wa_ens(1:DEF_DA_ENS_NUM,np) = wa_ens(1:DEF_DA_ENS_NUM,np) - &
               sum(wliq_soisno_ens(1:,1:DEF_DA_ENS_NUM,np) - wliq_ens_f, dim=1)

            deallocate(y_p, r_p, hx_p, trans)
            deallocate(wliq_ens_f, wliq_ens_a)
         ENDDO
      ENDIF

      ! Recalculate H(x) after LETKF so history can compare prior and analysis
      ! values for each configured observation stream.
      IF (p_is_worker .and. numpatch > 0) THEN
         wliq_soisno_a(:,:) = sum(wliq_soisno_ens(:,1:DEF_DA_ENS_NUM,:), dim=2) / DEF_DA_ENS_NUM
         wice_soisno_a(:,:) = sum(wice_soisno_ens(:,1:DEF_DA_ENS_NUM,:), dim=2) / DEF_DA_ENS_NUM
         t_soisno_a   (:,:) = sum(t_soisno_ens   (:,1:DEF_DA_ENS_NUM,:), dim=2) / DEF_DA_ENS_NUM
      ENDIF
      CALL sm%calcp()
      CALL sm%save(.true.)
      CALL sm%clear()

   END SUBROUTINE run_Assim_SM

!-----------------------------------------------------------------------------

   SUBROUTINE end_Assim_SM()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      CALL sm%clear()
      IF (allocated(sm%cfgs)) deallocate(sm%cfgs)
      IF (allocated(sm%pm%sensors)) deallocate(sm%pm%sensors)
      IF (allocated(sm%pm%cfgs)) deallocate(sm%pm%cfgs)

      sm%pm%nsensor = 0
      sm%use_pm = .false.
      sm%use_gnos = .false.

   END SUBROUTINE end_Assim_SM

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Assim_SM
#endif
