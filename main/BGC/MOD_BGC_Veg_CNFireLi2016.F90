#include <define.h>
#ifdef BGC
MODULE MOD_BGC_Veg_CNFireLi2016

!-------------------------------------------------------------------------------------------------------------------
! !DESCRIPTION:
! This module calculate burned area of each fire. The burned area is used to calculate fire induced CN loss rates
! in bgc_veg_CNFireBaseMod.F90
!
! !REFERENCES:
! Li, F., Levis, S., and Ward, D. S. 2013a. Quantifying the role of fire in the Earth system - Part 1: Improved global fire
! modeling in the Community Earth System Model (CESM1). Biogeosciences 10:2293-2314.
! Li, F., and Lawrence, D. 2017. Role of fire in the global land water budget during the 20th century through changing
! ecosystems. J. Clim. 30: 1894-1908.
!
! !ORIGINAL:
! The Community Land Model version 5.0 (CLM5)
!
! !REVISION:
! Xingjie Lu, 2021, revised the CLM5 code to be compatible with CoLM code structure.

   USE MOD_Precision
   USE MOD_TimeManager
   USE MOD_Namelist, only: DEF_Runoff_SCHEME, DEF_TOPMOD_method
   USE MOD_IncompleteGamma, only: GRATIO
   USE MOD_Const_Physical, only: tfrz
   USE MOD_Vars_Global, only: spval
   USE MOD_BGC_Vars_1DFluxes, only: fire_btran2
   USE MOD_Vars_1DForcing, only: &
       forc_q, forc_t, forc_psrf, forc_us, forc_vs
   USE MOD_Const_PFT, only: &
       isshrub, isgrass, isbetr, isbdtr, isbare, iscrop, isnatveg, &
       fd_pft, fsr_pft, rootfr_p, rswf_min, rswf_max
   USE MOD_Vars_TimeInvariants, only: &
       i_cwd, gdp_lf, abm_lf, peatf_lf, &
       lfuel, ufuel, cropfire_a1, borealat, troplat, &
       non_boreal_peatfire_c, boreal_peatfire_c, rh_low, rh_hgh, &
       bt_min, bt_max, pot_hmn_ign_counts_alpha, g0_fire, &
       psi0, porsl, bsw, fsatmax, fsatdcf, topoweti, &
       alp_twi, chi_twi, mu_twi
#ifdef vanGenuchten_Mualem_SOIL_MODEL
   USE MOD_Vars_TimeInvariants, only: theta_r, alpha_vgm, n_vgm, L_vgm, sc_vgm, fc_vgm
   USE MOD_Hydro_SoilFunction, only: soil_psi_from_vliq
#endif
   USE MOD_Vars_TimeVariables, only: &
       decomp_cpools_vr , totlitc    , totvegc   ,  cropf      , lfwt     , fuelc     , fuelc_crop , fsr     , &
       fd               , rootc      , lgdp      , lgdp1       , lpop     , wtlf      , &
       trotr1           , trotr2     , hdm_lf    , lnfm        , baf_crop , baf_peatf , &
       farea_burned     , nfire      , prec60    , wf2        , &
       tsoi17           , rh30       , prec30    , t_soisno    , wliq_soisno, zwt
   USE MOD_Vars_1DFluxes, only: frcsat
   USE MOD_BGC_Vars_PFTimeVariables, only: &
       burndate_p
   USE MOD_Vars_PFTimeInvariants, only: pftclass, pftfrac
   USE MOD_BGC_Vars_PFTimeVariables, only:  leafc_p     , leafc_storage_p     , leafc_xfer_p     , &
                                       frootc_p    , frootc_storage_p    , frootc_xfer_p    , &
                                       deadcrootc_p, deadcrootc_storage_p, deadcrootc_xfer_p, &
                                       livecrootc_p, livecrootc_storage_p, livecrootc_xfer_p  
#ifdef CROP
   USE MOD_BGC_Vars_PFTimeVariables, only: croplive_p
#endif

   USE MOD_Eroot, only: eroot
   USE MOD_Qsadv, only: qsadv

   IMPLICIT NONE

   PUBLIC CNFireArea

CONTAINS

   SUBROUTINE CNFireArea(i,ps,pe,dlat,nl_soil,idate,dz_soi)

   integer ,intent(in) :: i                  ! patch index
   integer ,intent(in) :: ps                 ! start pft index
   integer ,intent(in) :: pe                 ! END pft index
   real(r8),intent(in) :: dlat               ! latitude (degree)
   integer ,intent(in) :: nl_soil            ! number of total soil layers
   integer ,intent(in) :: idate(3)           ! current date (year, day of the year, seconds of the day)
   real(r8),intent(in) :: dz_soi(1:nl_soil)  ! thicknesses of each soil layer

   integer  :: g,l,c,p,j,fc,fp,kyr, kmo, kda, mcsec   ! index variables
   integer  :: ivt
   real(r8) :: dayspyr  ! days per year
   real(r8) :: fb       ! availability of fuel for regs A and C
   real(r8) :: fhd      ! impact of hd on agricultural fire
   real(r8) :: fgdp     ! impact of gdp on agricultural fire
   real(r8) :: fire_m   ! combustability of fuel for fire occurrence
   real(r8) :: spread_m ! combustability of fuel for fire spread
   real(r8) :: Lb_lf    ! length-to-breadth ratio added by Lifang
   real(r8) :: lh       ! anthro. ignitions (count/km2/hr)
   real(r8) :: fs       ! hd-dependent fires suppression (0-1)
   real(r8) :: ig       ! total ignitions (count/km2/hr)
   real(r8) :: arh, arh30 !combustability of fuel related to RH and RH30
   real(r8) :: afuel    !weight for arh and arh30
   real(r8) :: eq
   real(r8) :: deqdT
   real(r8) :: qsatq
   real(r8) :: qsatqdT
   real(r8) :: forc_rh
   real(r8) :: rootr(nl_soil)
   real(r8) :: rresis(nl_soil)
   real(r8) :: smp_node
   real(r8) :: s_node
   real(r8) :: tmp1d(nl_soil)
   real(r8) :: tmp0d
   real(r8) :: btran2
   real(r8) :: btran2_p(ps:pe)
   real(r8) :: satfrac_fire
   real(r8) :: eta_fire, pgr0_fire, pgr1_fire, qgr_fire, gfun_fire
   integer  :: niter_fire
   logical  :: satfrac_input_ok

   real(r8),parameter :: secsphr = 3600._r8
   real(r8),parameter :: secspday = 86400._r8
   real(r8),parameter :: smpswc  = -1.5e5
   real(r8),parameter :: smpsfc  = -3.3e3
   real(r8),parameter :: occur_hi_gdp_tree = 0.33d00
   real(r8),parameter :: nonborpeat_fire_precip_denom = 6.5
   real(r8),parameter :: borpeat_fire_soilmoist_denom = 0.35
   real(r8),parameter :: max_rh30_affecting_fuel = 95.
   real(r8),parameter :: ignition_efficiency     = 0.22
   real(r8),parameter :: prh30 = 0.6
   real(r8),parameter :: PI = 4.*atan(1.)
   real(r8),parameter :: topmod_vdcf = 2.0_r8
   integer m

      tsoi17 = forc_t(i)  ! Temporarily use air temperature for tsoi17, need to revised later.
      wf2    = 0.5        ! Temporarily set up, need to revise later.

      CALL julian2monthday(idate(1),idate(2),kmo,kda)

!      DO m = ps, pe
!         CALL eroot(nl_soil,0._r8,porsl(1:,i),&
!#ifdef Campbell_SOIL_MODEL
!                    bsw(1:,i),&
!#endif
!#ifdef vanGenuchten_Mualem_SOIL_MODEL
!                    theta_r(1:,i), alpha_vgm(1:,i), n_vgm(1:,i), L_vgm(1:,i), sc_vgm(1:,i), fc_vgm(1:,i), &
!#endif
!                    psi0(1:,i),rootfr_p(1:,pftclass(m)),dz_soi(1:),&
!                    t_soisno(1:,i),wliq_soisno(1:,i),tmp1d,tmp0d,btran2_p(m))
!      ENDDO
    !
    ! Calculate fraction of crop (cropf_col) and non-crop and non-bare-soil
    ! vegetation (lfwt) in vegetated column
    !
      cropf(i) = 0._r8
      lfwt (i) = 0._r8

      ! For crop veg types
      DO m = ps, pe
         IF( iscrop(pftclass(m)) )THEN
            cropf(i) = cropf(i) + pftfrac(m)
         ENDIF
      ! For natural vegetation (non-crop and non-bare-soil)
         IF( isnatveg(pftclass(m)))THEN
            lfwt (i) = lfwt(i) + pftfrac(m)
         ENDIF
      ENDDO

      !
      ! Calculate crop fuel
      !
      fuelc_crop(i)=0._r8

      ! For crop PFTs, fuel load includes leaf and litter; only
      ! column-level litter carbon
      ! is available, so we use leaf carbon to estimate the
      ! litter carbon for crop PFTs
      DO m = ps, pe
         IF( iscrop(pftclass(m)) .and. pftfrac(m) .gt. 0 .and. &
                 sum(leafc_p(ps:pe)*pftfrac(ps:pe)) > 0._r8 )THEN
            fuelc_crop(i) = fuelc_crop(i) + (leafc_p(m) + leafc_storage_p(m) + leafc_xfer_p(m))&
                          * pftfrac(m)/cropf(i) + totlitc(i) * leafc_p(m) &
                          / sum(leafc_p(ps:pe) * pftfrac(ps:pe)) * pftfrac(m) / cropf(i)
         ENDIF
      ENDDO
      !
      ! Calculate noncrop column variables
      !
      fsr   (i) = 0._r8
      fd    (i) = 0._r8
      rootc (i) = 0._r8
      lgdp  (i) = 0._r8
      lgdp1 (i) = 0._r8
      lpop  (i) = 0._r8
      wtlf  (i) = 0._r8
      trotr1(i) = 0._r8
      trotr2(i) = 0._r8
      btran2    = 0._r8

      do m = ps, pe
         btran2_p(m) = 0
      end do

      do m = ps, pe
         if(btran2_p(m) > 1._r8) then
            btran2_p(m) = 1._r8
         end if
      end do

      do m = ps, pe
         ivt = pftclass(m)
         if(isnatveg(ivt) .and. cropf(i) .lt. 1._r8)then
            do j = 1, nl_soil
               s_node = max(wliq_soisno(j,i)/(1000.*dz_soi(j)*porsl(j,i)),0.001)
               s_node = min(1., s_node)
#ifdef Campbell_SOIL_MODEL
               smp_node = max(smpswc, psi0(j,i)*s_node**(-bsw(j,i)))
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
               smp_node = soil_psi_from_vliq ( s_node*(porsl(j,i)-theta_r(j,i)) + theta_r(j,i), &
                  porsl(j,i), theta_r(j,i), psi0(j,i), &
                  5, (/alpha_vgm(j,i), n_vgm(j,i), L_vgm(j,i), sc_vgm(j,i), fc_vgm(j,i)/))
               smp_node = max(smpswc, smp_node)
#endif
            
!               btran2_p(m) = btran2_p(m) + rootfr_p(j,ivt) * max(0._r8, min((smp_node - smpswc) / &
!                    (smpsfc - smpswc),1._r8))
               btran2_p(m) = btran2_p(m) + rootfr_p(j,ivt) * s_node
            end do
         end if
         btran2 = btran2  + max(0._r8, min(1._r8, &
             (btran2_p(m) -rswf_min(ivt))/(rswf_max(ivt) &
              -rswf_min(ivt)))) * pftfrac(m)   
         wtlf(i)  = wtlf(i) + pftfrac(m)
      end do

      ! Save the normalized wetness used below for fire diagnostics only.
      fire_btran2(i) = spval
      IF (wtlf(i) > 0._r8) fire_btran2(i) = btran2 / wtlf(i)

      ! Warning : ivt is not initialized.
      ! For non-crop -- natural vegetation and bare-soil
      do m = ps, pe
         ivt = pftclass(m)
         if(isnatveg(ivt) .and. cropf(i) .lt. 1._r8)then
            IF( isbetr(ivt) )THEN
               trotr1(i)= trotr1(i) + pftfrac(m)
            ENDIF
            IF( isbdtr(ivt) .and. abs(dlat) .lt. troplat)THEN
               trotr2(i)= trotr2(i) + pftfrac(m)
            ENDIF

            rootc(i) = rootc(i) + (frootc_p(m) + frootc_storage_p(m) + frootc_xfer_p(m) &
                 + deadcrootc_p(m) + deadcrootc_storage_p(m) + deadcrootc_xfer_p(m) &
                 + livecrootc_p(m) + livecrootc_storage_p(m) + livecrootc_xfer_p(m)) * pftfrac(m)

            fsr(i) = fsr(i) + fsr_pft(ivt) * pftfrac(m) / (1._r8-cropf(i))

         ! all these constants are in Li et al. BG (2012a,b;2013)

            IF( hdm_lf(i)  >  0.1_r8 )THEN
               ! For not bare-soil
               IF(.not. isbare(ivt) )THEN
                  ! For shrub and grass (crop already excluded above)
                  IF( isshrub(ivt) .or. isgrass(ivt) )THEN      !for shurb and grass
                     lgdp(i)  = lgdp(i) + (0.1_r8 + 0.9_r8*    &
                                 exp(-1._r8*PI* (gdp_lf(i)/8._r8)**0.5_r8)) &
                                 * pftfrac(m) / (1.0_r8-cropf(i))
                     lgdp1(i) = lgdp1(i) + (0.2_r8 + 0.8_r8*   &
                                 exp(-1._r8*PI* (gdp_lf(i)/7._r8))) &
                                 * pftfrac(m) / (1.0_r8-cropf(i))
                     lpop(i)  = lpop(i) + (0.2_r8 + 0.8_r8*    &
                                 exp(-1._r8*PI* (hdm_lf(i)/450._r8)**0.5_r8)) &
                                 * pftfrac(m) / (1.0_r8-cropf(i))
                  ELSE   ! for trees
                     IF( gdp_lf(i)  >  20._r8 )THEN
                        lgdp(i) =lgdp(i)+occur_hi_gdp_tree*pftfrac(m)/(1._r8-cropf(i))
                        lgdp1(i) =lgdp1(i)+0.62_r8*pftfrac(m)/(1._r8-cropf(i))
                     ELSE
                        IF( gdp_lf(i) > 8._r8 )THEN
                           lgdp(i)=lgdp(i)+0.79_r8*pftfrac(m)/(1._r8-cropf(i))
                           lgdp1(i)=lgdp1(i)+0.83_r8*pftfrac(m)/(1._r8-cropf(i))
                        ELSE
                           lgdp(i) = lgdp(i)+pftfrac(m)/(1._r8-cropf(i))
                           lgdp1(i)=lgdp1(i)+pftfrac(m)/(1._r8-cropf(i))
                        ENDIF
                     ENDIF
                     lpop(i) = lpop(i) + (0.4_r8 + 0.6_r8*    &
                                   exp(-1._r8*PI* & 
                                   (hdm_lf(i)/125._r8)))*pftfrac(m)/(1._r8-cropf(i))
                  ENDIF
               ENDIF
            ELSE
               lgdp(i)  = lgdp(i)  + pftfrac(m)/(1._r8-cropf(i))
               lgdp1(i) = lgdp1(i) + pftfrac(m)/(1._r8-cropf(i))
               lpop(i)  = lpop(i)  + pftfrac(m)/(1._r8-cropf(i))
            ENDIF

            fd(i) = fd(i) + fd_pft(ivt) *pftfrac(m) * secsphr / (1.0_r8-cropf(i))
         ENDIF
      ENDDO
      !
      ! calculate burned area fraction in cropland
      !
      ! Reset diagnostics also when the non-crop calculation is skipped.
      nfire(i) = 0._r8
      fuelc(i) = 0._r8
      baf_crop(i)=0._r8

      DO m = ps, pe
         IF( kmo == 1 .and. kda == 1 .and. idate(3) == 0 )THEN
            burndate_p(m) = 10000 ! init. value; actual range [0 365]
         ENDIF
      ENDDO

      ! For crop
      DO m = ps, pe
         ivt = pftclass(m)
         IF( forc_t(i)  >=  tfrz .and. iscrop(ivt) .and. kmo == abm_lf(i) &
             .and. burndate_p(m) >= 999 .and. pftfrac(m) .gt. 0._r8)THEN ! catch  crop burn time

         ! calculate human density impact on ag. fire
            fhd = 0.2_r8+0.8_r8*exp(-1._r8*PI*(hdm_lf(i)/400._r8))

         ! calculate impact of GDP on ag. fire
            fgdp = 0.05_r8+0.95_r8*exp(-1._r8*PI*(gdp_lf(i)/20._r8))

         ! calculate burned area
         ! crop fire only for generic crop types at this time
         ! managed crops are treated as grasses IF crop model is turned on
#ifdef CROP
            if(.not. croplive_p(m))then
               burndate_p(m)  =  kda
               baf_crop(i) = baf_crop(i) + cropfire_a1/secsphr*fhd*fgdp*pftfrac(m)
            endif
#else
            burndate_p(m)  =  kda
            baf_crop(i) = baf_crop(i) + cropfire_a1/secsphr*fhd*fgdp*pftfrac(m)
#endif
         ENDIF
      ENDDO

      !
      ! calculate saturated fraction used by boreal peat fire.
      ! For TOPMODEL, reproduce runoff fsat exactly. For other runoff schemes,
      ! use hydrology frcsat as the scheme-consistent saturated-area fraction.
      ! Fall back safely only when the hydrology value is unavailable.
      !
      satfrac_fire = 0._r8
      satfrac_input_ok = .false.
      IF (DEF_Runoff_SCHEME == 0) THEN
         IF ((zwt(i) /= spval) .AND. (zwt(i) == zwt(i)) .AND. &
             (ABS(zwt(i)) < 1.0e30_r8)) THEN
            IF ((DEF_TOPMOD_method == 0) .OR. (DEF_TOPMOD_method == 1)) THEN
               satfrac_input_ok = (fsatmax(i) /= spval) .AND. &
                    (fsatdcf(i) /= spval) .AND. &
                    (fsatmax(i) == fsatmax(i)) .AND. &
                    (fsatdcf(i) == fsatdcf(i)) .AND. &
                    (ABS(fsatmax(i)) < 1.0e30_r8) .AND. &
                    (ABS(fsatdcf(i)) < 1.0e30_r8)
               IF (satfrac_input_ok) THEN
                  satfrac_fire = fsatmax(i) * EXP(-fsatdcf(i) * topmod_vdcf * zwt(i))
               ENDIF
            ELSE
               satfrac_input_ok = (topoweti(i) /= spval) .AND. &
                    (alp_twi(i) /= spval) .AND. (chi_twi(i) /= spval) .AND. &
                    (mu_twi(i) /= spval) .AND. &
                    (topoweti(i) == topoweti(i)) .AND. &
                    (alp_twi(i) == alp_twi(i)) .AND. &
                    (chi_twi(i) == chi_twi(i)) .AND. &
                    (mu_twi(i) == mu_twi(i)) .AND. &
                    (ABS(topoweti(i)) < 1.0e30_r8) .AND. &
                    (ABS(alp_twi(i)) < 1.0e30_r8) .AND. &
                    (ABS(chi_twi(i)) < 1.0e30_r8) .AND. &
                    (ABS(mu_twi(i)) < 1.0e30_r8) .AND. &
                    (alp_twi(i) > 0._r8) .AND. (chi_twi(i) > 0._r8)
               IF (satfrac_input_ok) THEN
                  IF (zwt(i) <= 0._r8) THEN
                     satfrac_fire = 1._r8
                  ELSE
                     eta_fire = topoweti(i)
                     gfun_fire = 0._r8
                     DO niter_fire = 1, 20
                        CALL GRATIO(alp_twi(i)+1._r8, &
                             (eta_fire-mu_twi(i))/chi_twi(i), pgr1_fire, qgr_fire, 0)
                        CALL GRATIO(alp_twi(i), &
                             (eta_fire-mu_twi(i))/chi_twi(i), pgr0_fire, qgr_fire, 0)
                        gfun_fire = ((eta_fire-mu_twi(i))*pgr0_fire &
                             - chi_twi(i)*alp_twi(i)*pgr1_fire) / topmod_vdcf - zwt(i)
                        IF (ABS(gfun_fire) <= 1.e-6_r8 .OR. pgr0_fire <= 0._r8) EXIT
                        eta_fire = mu_twi(i) + &
                             (chi_twi(i)*alp_twi(i)*pgr1_fire + topmod_vdcf*zwt(i)) / pgr0_fire
                     ENDDO
                     CALL GRATIO(alp_twi(i), &
                          (eta_fire-mu_twi(i))/chi_twi(i), pgr0_fire, qgr_fire, 0)
                     satfrac_fire = qgr_fire
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      IF (.NOT. satfrac_input_ok) THEN
         IF ((frcsat(i) /= spval) .AND. (frcsat(i) == frcsat(i)) .AND. &
             (ABS(frcsat(i)) < 1.0e30_r8)) THEN
            satfrac_fire = frcsat(i)
         ELSE
            satfrac_fire = 0._r8
         ENDIF
      ENDIF

      IF ((satfrac_fire /= satfrac_fire) .OR. (ABS(satfrac_fire) >= 1.0e30_r8)) THEN
         IF ((frcsat(i) /= spval) .AND. (frcsat(i) == frcsat(i)) .AND. &
             (ABS(frcsat(i)) < 1.0e30_r8)) THEN
            satfrac_fire = frcsat(i)
         ELSE
            satfrac_fire = 0._r8
         ENDIF
      ENDIF
      satfrac_fire = MAX(0._r8, MIN(1._r8, satfrac_fire))

      !
      ! calculate peatland fire
      !
      IF(dlat < borealat )THEN
         if((trotr1(i)+trotr2(i)) .le. 0.8_r8 .and. trotr1(i)+trotr2(i) .gt. 0._r8)then
            baf_peatf(i) = non_boreal_peatfire_c/secsphr*max(0._r8, &
                 min(1._r8,(1.0_r8-prec30(i)*secspday/nonborpeat_fire_precip_denom)))*peatf_lf(i)
         else
            baf_peatf(i) = 0._r8
         end if
      ELSE
         baf_peatf(i) = boreal_peatfire_c/secsphr*exp(-PI*(max(wf2(i),0._r8)/borpeat_fire_soilmoist_denom))* &
                 max(0._r8,min(1._r8,(tsoi17(i)-tfrz)/10._r8))*peatf_lf(i)* &
                 (1._r8-satfrac_fire)
      ENDIF
      !
      ! calculate other fires
      !
      CALL qsadv(forc_t(i),forc_psrf(i),eq,deqdT,qsatq,qsatqdT)
      forc_rh = forc_q(i) / eq

      IF( cropf(i)  <  1._r8 )THEN
         fuelc(i) = totlitc(i)+totvegc(i)-rootc(i)-fuelc_crop(i)*cropf(i)
         DO j = 1, nl_soil
            fuelc(i) = fuelc(i)+decomp_cpools_vr(j,i_cwd,i) * dz_soi(j)
         ENDDO
         fuelc(i) = fuelc(i)/(1._r8-cropf(i))
         fb       = max(0.0_r8,min(1.0_r8,(fuelc(i)-lfuel)/(ufuel-lfuel)))
         afuel  =min(1._r8,max(0._r8,(fuelc(i)-2500._r8)/(5000._r8-2500._r8)))
         arh=1._r8-max(0._r8, min(1._r8,(forc_rh-rh_low)/(rh_hgh-rh_low)))
         arh30=1._r8-max(prh30, min(1._r8,rh30(i)/max_rh30_affecting_fuel))
         IF (forc_rh < rh_hgh .and. wtlf(i) > 0._r8 .and. tsoi17(i)> tfrz)THEN
            fire_m   = ((afuel*arh30+(1._r8-afuel)*arh)**1.5_r8) &
                    *((1._r8 - btran2/wtlf(i))**0.5_r8)
         ELSE
            fire_m   = 0._r8
         ENDIF
         lh       = pot_hmn_ign_counts_alpha*6.8_r8*hdm_lf(i)**(0.43_r8)/30._r8/24._r8
         fs       = 1._r8-(0.01_r8+0.98_r8*exp(-0.025_r8*hdm_lf(i)))
         IF (trotr1(i)+trotr2(i)<=0.6_r8) THEN
            ig    = (lh+lnfm(i)/(5.16_r8+2.16_r8*cos(PI/180._r8*3*min(60._r8,abs(dlat/PI*180))))*ignition_efficiency)  &
                    *(1._r8-fs)*(lfwt(i)**0.5)
         ELSE
            ig    = (lnfm(i)/(5.16_r8+2.16_r8*cos(PI/180._r8*3*min(60._r8,abs(dlat/PI*180))))*ignition_efficiency)  &
                    *(1._r8-fs)*(lfwt(i)**0.5)
         ENDIF
         nfire(i) = ig/secsphr*fb*fire_m*lgdp(i) !fire counts/km2/sec
         Lb_lf    = 1._r8+10._r8*(1._r8-EXP(-0.06_r8*sqrt(forc_us(i)*forc_us(i)+forc_vs(i)*forc_vs(i))))
         spread_m = fire_m**0.5_r8
         fd(i)    = (lfwt(i)*lgdp1(i)*lpop(i))**0.5_r8 * fd(i)
         farea_burned(i) = min(1._r8,(g0_fire*spread_m*fsr(i)* &
                    fd(i)/1000._r8)**2*nfire(i)*PI*Lb_lf+ &
                    baf_crop(i)+baf_peatf(i))  ! fraction (0-1) per sec
      ELSE
         farea_burned(i) = min(1._r8,baf_crop(i)+baf_peatf(i))
      ENDIF
   END SUBROUTINE CNFireArea

END MODULE MOD_BGC_Veg_CNFireLi2016
#endif


