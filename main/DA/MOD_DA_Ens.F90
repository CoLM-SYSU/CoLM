#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Ens
!-----------------------------------------------------------------------
! DESCRIPTION:
!    Provide functions to generate ensemble samples for data assimilation
!
! REFERENCES:
!    [1] Algorithm Theoretical Basis Document Level 4 Surface and Root
!        Zone Soil Moisture (L4_SM) Data Product
!
! AUTHOR:
!   Lu Li, 12/2024: Initial version
!   Lu Li, 10/2025: Consider correlation & AR(1) process
!-----------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_Vars_TimeVariables
   USE MOD_DA_Vars_TimeVariables
   USE MOD_Vars_1DForcing
   USE MOD_LandPatch
   USE MOD_SPMD_Task
   IMPLICIT NONE
   SAVE

   ! public functions
   PUBLIC :: ensemble
   PUBLIC :: model_to_ens_member
   PUBLIC :: ens_member_to_model

   ! local parameters
   ! forcing [parameters used here is consistent with SMAP L4 (Table 4 in [1])]
   integer,  parameter :: nvar = 4                                   ! number of pertutated forcing variables
   real(r8), parameter :: tau_ar = 24.0                              ! correlation time scale in hours

   real(r8) :: dt                                                    ! time step in hours
   real(r8) :: phi                                                   ! AR(1) autocorrelation coefficient consider time scale
   real(r8) :: sigma_eps                                             ! standard deviation of noise in AR(1) process
   real(r8), dimension(nvar) :: sigma = (/0.5, 0.3, 20.0, 1.0/)      ! standard deviation of perturbed forcing variables (prcp, sw, lw, t)
   real(r8), dimension(nvar, nvar) :: C = reshape([ &
         1.0, -0.8,  0.5, 0.0, &
        -0.8,  1.0, -0.5, 0.4, &
         0.5, -0.5,  1.0, 0.4, &
         0.0,  0.4,  0.4, 1.0], shape=[nvar, nvar])                  ! cross-correlation matrix between perturbed forcing variables
   real(r8), allocatable :: r_prev(:,:,:)                            ! previous perturbation (numpatch, nvar, DEF_DA_ENS_NUM)
   real(r8), allocatable :: r_curr(:,:,:)                            ! current perturbation (numpatch, nvar, DEF_DA_ENS_NUM)
   logical :: initialized = .false.                                  ! flag to indicate if is initialized

   ! soil moisture [default set 0.002 m3/m3 disterbulance]
   integer,  parameter :: nvar_sm = 2                                ! number of pertutated soil moisture layers
   real(r8), parameter :: tau_sm = 3.0                               ! correlation time scale in hours
   real(r8) :: phi_sm                                                ! AR(1) autocorrelation coefficient consider time scale
   real(r8) :: sigma_eps_sm                                          ! standard deviation of noise in AR(1) process
   real(r8), dimension(nvar_sm) :: sigma_sm = (/0.035, 0.0552/)      ! standard deviation of perturbed soil moisture (equal to 0.002 m3/m3)
   real(r8), allocatable :: r_prev_sm(:,:,:)                         ! previous perturbation (numpatch, nvar_sm, DEF_DA_ENS_NUM)
   real(r8), allocatable :: r_curr_sm(:,:,:)                         ! current perturbation (numpatch, nvar_sm, DEF_DA_ENS_NUM)

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE ensemble(deltim)

!-----------------------------------------------------------------------
      IMPLICIT NONE

      real(r8), intent(in) :: deltim

!------------------------ Local Variables ------------------------------
      integer  ::  np, i, j

      real(r8) :: cov_matrix(nvar, nvar)                      ! covariance matrix between perturbed forcing variables
      real(r8) :: L(nvar, nvar)                               ! Cholesky decomposition of correlation matrix
      integer  :: info                                        ! info flag for Cholesky decomposition
      real(r8) :: u1(DEF_DA_ENS_NUM/2), u2(DEF_DA_ENS_NUM/2)  ! uniform random variables
      real(r8) :: z(nvar, DEF_DA_ENS_NUM)                     ! standard normal random variables
      real(r8) :: mean_z(nvar)                                ! mean of perturbation (nvar)
      real(r8) :: std_z(nvar)                                 ! std of perturbation (nvar)
      real(r8) :: zxL(nvar, DEF_DA_ENS_NUM)                   ! correlated random variables (nvar, DEF_DA_ENS_NUM)

      real(r8) :: z_sm(DEF_DA_ENS_NUM)                        ! standard normal random variables
      real(r8) :: mean_z_sm                                   ! mean of perturbation
      real(r8) :: std_z_sm                                    ! std of perturbation
      real(r8) :: mean_r_sm(nvar_sm)                          ! mean of perturbation for soil moisture (nvar_sm)
      real(r8) :: std_r_sm(nvar_sm)                           ! std of perturbation for soil moisture (nvar_sm)
      real(r8) :: a1(DEF_DA_ENS_NUM)                          ! temporary disturbed variable for soil moisture layer 1
      real(r8) :: a2(DEF_DA_ENS_NUM)                          ! temporary disturbed variable for soil moisture layer 2

!-----------------------------------------------------------------------

      IF (mod(DEF_DA_ENS_NUM, 2) /= 0) THEN
         print *, 'Error: DEF_DA_ENS_NUM must be even for Box-Muller ensemble sampling.'
         CALL CoLM_stop()
      ENDIF

      ! initialize persistent variables
      IF (.not. initialized) THEN
         CALL random_seed()
         allocate(r_prev(numpatch, nvar, DEF_DA_ENS_NUM))
         allocate(r_curr(numpatch, nvar, DEF_DA_ENS_NUM))
         allocate(r_prev_sm(numpatch, nvar_sm, DEF_DA_ENS_NUM))
         allocate(r_curr_sm(numpatch, nvar_sm, DEF_DA_ENS_NUM))
         r_prev = 0.0_r8
         r_curr = 0.0_r8
         r_prev_sm = 0.0_r8
         r_curr_sm = 0.0_r8
         initialized = .true.
      ENDIF

      ! calculate AR(1) parameters
      dt = deltim/3600
      phi = exp(-dt/tau_ar)
      sigma_eps = sqrt(1.0 - phi**2)
      phi_sm = exp(-dt/tau_sm)
      sigma_eps_sm = sqrt(1.0 - phi_sm**2)

      ! calculate covariance matrix by cross-correlation matrix and standard deviation
      cov_matrix = 0.0_r8
      DO i = 1, nvar
         DO j = 1, nvar
            cov_matrix(i,j) = C(i,j) * sigma(i) * sigma(j)
         ENDDO
      ENDDO

      ! perform Cholesky decomposition of covariance matrix
      L = cov_matrix
      CALL dpotrf('L', nvar, L, nvar, info)
      DO i = 1, nvar
         DO j = i+1, nvar
            L(i,j) = 0.0_r8
         ENDDO
      ENDDO
      IF (info /= 0) THEN
         print *, 'Error: Cholesky decomposition failed'
         CALL CoLM_stop()
      ENDIF

      ! Generate ensemble samples for forcing variables
      DO np = 1, numpatch
         ! generate disturbance ensemble samples ~ N(0, I)
         CALL random_number(u1)
         CALL random_number(u2)
         DO i = 1, DEF_DA_ENS_NUM/2
            u1(i) = max(u1(i), 1e-10)  ! ensure u1 is not zero
            z(1,i*2-1) = sqrt(-2.0*log(u1(i))) * cos(2.0*pi*u2(i))
            z(1,i*2)   = sqrt(-2.0*log(u1(i))) * sin(2.0*pi*u2(i))
         ENDDO
         CALL random_number(u1)
         CALL random_number(u2)
         DO i = 1, DEF_DA_ENS_NUM/2
            u1(i) = max(u1(i), 1e-10)
            z(2,i*2-1) = sqrt(-2.0*log(u1(i))) * cos(2.0*pi*u2(i))
            z(2,i*2)   = sqrt(-2.0*log(u1(i))) * sin(2.0*pi*u2(i))
         ENDDO
         CALL random_number(u1)
         CALL random_number(u2)
         DO i = 1, DEF_DA_ENS_NUM/2
            u1(i) = max(u1(i), 1e-10)
            z(3,i*2-1) = sqrt(-2.0*log(u1(i))) * cos(2.0*pi*u2(i))
            z(3,i*2)   = sqrt(-2.0*log(u1(i))) * sin(2.0*pi*u2(i))
         ENDDO
         CALL random_number(u1)
         CALL random_number(u2)
         DO i = 1, DEF_DA_ENS_NUM/2
            u1(i) = max(u1(i), 1e-10)
            z(4,i*2-1) = sqrt(-2.0*log(u1(i))) * cos(2.0*pi*u2(i))
            z(4,i*2)   = sqrt(-2.0*log(u1(i))) * sin(2.0*pi*u2(i))
         ENDDO

         ! normalize z to mean 0 and std 1
         DO i = 1, nvar
            mean_z(i) = sum(z(i, :))/DEF_DA_ENS_NUM
            std_z(i) = sqrt(sum((z(i, :) - mean_z(i))**2)/(DEF_DA_ENS_NUM - 1))
         ENDDO
         DO i = 1, nvar
            z(i,:) = (z(i,:)-mean_z(i))/std_z(i)
         ENDDO

         ! multiply by Cholesky factor to introduce correlation (z*L)
         CALL dgemm('N', 'N', nvar, DEF_DA_ENS_NUM, nvar, 1.0_r8, L, nvar, z, nvar, 0.0_r8, zxL, nvar)

         ! introduce correlation using AR(1) process
         DO i = 1, nvar
            DO j = 1, DEF_DA_ENS_NUM
               r_curr(np,i,j) = phi * r_prev(np,i,j) + sigma_eps * zxL(i,j)
            ENDDO
         ENDDO
         ! no AR(1) process, directly use correlated random variables
         ! r_curr(np,:,:) = zxL

         ! normalize the disturbance ensemble samples to mean 0
         mean_z = sum(r_curr(np,:,:), dim=2)/DEF_DA_ENS_NUM
         DO i = 1, nvar
            DO j = 1, DEF_DA_ENS_NUM
               r_curr(np,i,j) = r_curr(np,i,j) - mean_z(i)
            ENDDO
         ENDDO

         ! save current perturbation as previous perturbation for next time step
         r_prev = r_curr

         ! generate ensemble samples according different types
         DO j = 1, DEF_DA_ENS_NUM
            forc_prc_ens(j,np) = forc_prc(np) * exp(r_curr(np,1,j) - 0.5 * sigma(1)**2)
            forc_prl_ens(j,np) = forc_prl(np) * exp(r_curr(np,1,j) - 0.5 * sigma(1)**2)
            forc_sols_ens(j,np) = forc_sols(np) * exp(r_curr(np,2,j) - 0.5 * sigma(2)**2)
            forc_soll_ens(j,np) = forc_soll(np) * exp(r_curr(np,2,j) - 0.5 * sigma(2)**2)
            forc_solsd_ens(j,np) = forc_solsd(np) * exp(r_curr(np,2,j) - 0.5 * sigma(2)**2)
            forc_solld_ens(j,np) = forc_solld(np) * exp(r_curr(np,2,j) - 0.5 * sigma(2)**2)
            forc_frl_ens(j,np) = forc_frl(np) + r_curr(np,3,j)
            forc_t_ens(j,np) = forc_t(np) + r_curr(np,4,j)
         ENDDO

         IF (DEF_DA_ENS_SM) THEN
            ! generate ensemble samples (0, I) for soil moisture
            CALL random_number(u1)
            CALL random_number(u2)
            DO i = 1, DEF_DA_ENS_NUM/2
               u1(i) = max(u1(i), 1e-10)
               z_sm(i*2-1) = sqrt(-2.0*log(u1(i))) * cos(2.0*pi*u2(i))
               z_sm(i*2)   = sqrt(-2.0*log(u1(i))) * sin(2.0*pi*u2(i))
            ENDDO
            mean_z_sm = sum(z_sm)/DEF_DA_ENS_NUM
            std_z_sm = sqrt(sum((z_sm - mean_z_sm)**2)/(DEF_DA_ENS_NUM - 1))
            z_sm = (z_sm - mean_z_sm)/std_z_sm

            ! introduce correlation using AR(1) process
            DO i = 1, nvar_sm
               DO j = 1, DEF_DA_ENS_NUM
                  r_curr_sm(np,i,j) = phi_sm * r_prev_sm(np,i,j) + sigma_eps_sm * sigma_sm(i) * z_sm(j)
               ENDDO
            ENDDO

            ! normalize the disturbance ensemble samples to mean 0
            mean_r_sm = sum(r_curr_sm(np,:,:), dim=2)/DEF_DA_ENS_NUM
            DO i = 1, nvar_sm
               DO j = 1, DEF_DA_ENS_NUM
                  r_curr_sm(np,i,j) = r_curr_sm(np,i,j) - mean_r_sm(i)
               ENDDO
            ENDDO

            ! save current perturbation as previous perturbation for next time step
            r_prev_sm = r_curr_sm

            ! generate ensemble samples according different types
            DO j = 1, DEF_DA_ENS_NUM
               a1(j) = wliq_soisno_ens(1,j,np) + r_curr_sm(np,1,j)
               a2(j) = wliq_soisno_ens(2,j,np) + r_curr_sm(np,2,j)
            ENDDO
            DO j = 1, DEF_DA_ENS_NUM
               a1(j) = max(1e-10, a1(j))
               a2(j) = max(1e-10, a2(j))
            ENDDO

            ! move residual water to water table
            DO j = 1, DEF_DA_ENS_NUM
               wa_ens(j, np) = wa_ens(j, np) - (sum(wliq_soisno_ens(1:2, j, np)) - a1(j) - a2(j))
               wliq_soisno_ens(1, j, np) = a1(j)
               wliq_soisno_ens(2, j, np) = a2(j)
            ENDDO
         ENDIF
      ENDDO

   END SUBROUTINE ensemble

!-----------------------------------------------------------------------

   SUBROUTINE model_to_ens_member(iens)

!-----------------------------------------------------------------------
      USE MOD_Vars_1DFluxes
      USE MOD_DA_Vars_1DFluxes
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------
      integer, intent(in) :: iens

!-----------------------------------------------------------------------
      z_sno_ens       (  :,iens,:) = z_sno
      dz_sno_ens      (  :,iens,:) = dz_sno
      t_soisno_ens    (  :,iens,:) = t_soisno
      wliq_soisno_ens (  :,iens,:) = wliq_soisno
      wice_soisno_ens (  :,iens,:) = wice_soisno
      smp_ens         (  :,iens,:) = smp
      hk_ens          (  :,iens,:) = hk
      t_grnd_ens      (    iens,:) = t_grnd
      tleaf_ens       (    iens,:) = tleaf
      ldew_ens        (    iens,:) = ldew
      ldew_rain_ens   (    iens,:) = ldew_rain
      ldew_snow_ens   (    iens,:) = ldew_snow
      fwet_snow_ens   (    iens,:) = fwet_snow
      sag_ens         (    iens,:) = sag
      scv_ens         (    iens,:) = scv
      snowdp_ens      (    iens,:) = snowdp
      fveg_ens        (    iens,:) = fveg
      fsno_ens        (    iens,:) = fsno
      sigf_ens        (    iens,:) = sigf
      green_ens       (    iens,:) = green
      tlai_ens        (    iens,:) = tlai
      lai_ens         (    iens,:) = lai
      tsai_ens        (    iens,:) = tsai
      sai_ens         (    iens,:) = sai
      alb_ens         (:,:,iens,:) = alb
      ssun_ens        (:,:,iens,:) = ssun
      ssha_ens        (:,:,iens,:) = ssha
      ssoi_ens        (:,:,iens,:) = ssoi
      ssno_ens        (:,:,iens,:) = ssno
      thermk_ens      (    iens,:) = thermk
      extkb_ens       (    iens,:) = extkb
      extkd_ens       (    iens,:) = extkd
      zwt_ens         (    iens,:) = zwt
      wdsrf_ens       (    iens,:) = wdsrf
      wa_ens          (    iens,:) = wa
      wetwat_ens      (    iens,:) = wetwat
      t_lake_ens      (  :,iens,:) = t_lake
      lake_icefrac_ens(  :,iens,:) = lake_icefrac
      savedtke1_ens   (    iens,:) = savedtke1

      trad_ens  (iens,:)   = trad
      tref_ens  (iens,:)   = tref
      qref_ens  (iens,:)   = qref
      ustar_ens (iens,:)   = ustar
      qstar_ens (iens,:)   = qstar
      tstar_ens (iens,:)   = tstar
      fm_ens    (iens,:)   = fm
      fh_ens    (iens,:)   = fh
      fq_ens    (iens,:)   = fq

      forc_t_ens    (iens,:) = forc_t
      forc_prc_ens  (iens,:) = forc_prc
      forc_prl_ens  (iens,:) = forc_prl
      forc_sols_ens (iens,:) = forc_sols
      forc_soll_ens (iens,:) = forc_soll
      forc_solsd_ens(iens,:) = forc_solsd
      forc_solld_ens(iens,:) = forc_solld
      forc_frl_ens  (iens,:) = forc_frl

      fsena_ens (iens,:) = fsena
      lfevpa_ens(iens,:) = lfevpa
      fevpa_ens (iens,:) = fevpa
      rsur_ens  (iens,:) = rsur

   END SUBROUTINE model_to_ens_member

!-----------------------------------------------------------------------

   SUBROUTINE ens_member_to_model(iens)

!-----------------------------------------------------------------------
      USE MOD_Vars_1DFluxes
      USE MOD_DA_Vars_1DFluxes
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------
      integer, intent(in) :: iens

!-----------------------------------------------------------------------
      forc_t     = forc_t_ens(iens,:)
      forc_prc   = forc_prc_ens(iens,:)
      forc_prl   = forc_prl_ens(iens,:)
      forc_sols  = forc_sols_ens(iens,:)
      forc_soll  = forc_soll_ens(iens,:)
      forc_solsd = forc_solsd_ens(iens,:)
      forc_solld = forc_solld_ens(iens,:)
      forc_frl   = forc_frl_ens(iens,:)

      z_sno        = z_sno_ens       (  :,iens,:)
      dz_sno       = dz_sno_ens      (  :,iens,:)
      t_soisno     = t_soisno_ens    (  :,iens,:)
      wliq_soisno  = wliq_soisno_ens (  :,iens,:)
      wice_soisno  = wice_soisno_ens (  :,iens,:)
      smp          = smp_ens         (  :,iens,:)
      hk           = hk_ens          (  :,iens,:)
      t_grnd       = t_grnd_ens      (    iens,:)
      tleaf        = tleaf_ens       (    iens,:)
      ldew         = ldew_ens        (    iens,:)
      ldew_rain    = ldew_rain_ens   (    iens,:)
      ldew_snow    = ldew_snow_ens   (    iens,:)
      fwet_snow    = fwet_snow_ens   (    iens,:)
      sag          = sag_ens         (    iens,:)
      scv          = scv_ens         (    iens,:)
      snowdp       = snowdp_ens      (    iens,:)
      fveg         = fveg_ens        (    iens,:)
      fsno         = fsno_ens        (    iens,:)
      sigf         = sigf_ens        (    iens,:)
      green        = green_ens       (    iens,:)
      tlai         = tlai_ens        (    iens,:)
      lai          = lai_ens         (    iens,:)
      tsai         = tsai_ens        (    iens,:)
      sai          = sai_ens         (    iens,:)
      alb          = alb_ens         (:,:,iens,:)
      ssun         = ssun_ens        (:,:,iens,:)
      ssha         = ssha_ens        (:,:,iens,:)
      ssoi         = ssoi_ens        (:,:,iens,:)
      ssno         = ssno_ens        (:,:,iens,:)
      thermk       = thermk_ens      (    iens,:)
      extkb        = extkb_ens       (    iens,:)
      extkd        = extkd_ens       (    iens,:)
      zwt          = zwt_ens         (    iens,:)
      wdsrf        = wdsrf_ens       (    iens,:)
      wa           = wa_ens          (    iens,:)
      wetwat       = wetwat_ens      (    iens,:)
      t_lake       = t_lake_ens      (  :,iens,:)
      lake_icefrac = lake_icefrac_ens(  :,iens,:)
      savedtke1    = savedtke1_ens   (    iens,:)

      trad   = trad_ens(iens,:)
      tref   = tref_ens(iens,:)
      qref   = qref_ens(iens,:)
      ustar  = ustar_ens(iens,:)
      qstar  = qstar_ens(iens,:)
      tstar  = tstar_ens(iens,:)
      fm     = fm_ens(iens,:)
      fh     = fh_ens(iens,:)
      fq     = fq_ens(iens,:)

      fsena  = fsena_ens(iens,:)
      lfevpa = lfevpa_ens(iens,:)
      fevpa  = fevpa_ens(iens,:)
      rsur   = rsur_ens(iens,:)

   END SUBROUTINE ens_member_to_model

!-----------------------------------------------------------------------
END MODULE MOD_DA_Ens
#endif
