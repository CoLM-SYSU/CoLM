#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Operator_RTM
!-----------------------------------------------------------------------
! DESCRIPTION:
!    Forward modeling of brightness temperature observations
!
! AUTHOR:
!   Lu Li, 12/2024: Initial version
!   Zhilong Fan, Lu Li, 03/2024: Debug and clean codes
!   Lu Li, 10/2025: Debug and clean codes
!-----------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Const_Physical
   USE MOD_Vars_1DForcing
   USE MOD_DA_Operator_Const
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only: nl_soil, nl_lake, N_land_classification, maxsnl, spval, dz_soi
   USE MOD_Namelist
   IMPLICIT NONE
   SAVE

! public functions
   PUBLIC   :: PM_RTM

! local variables (parameters depends on frequency and incidence angle of satellite)
   real(r8) :: fghz                       ! frequency of satellite (GHz)
   real(r8) :: theta                      ! incidence angle of satellite (rad)
   real(r8) :: f                          ! frequency (Hz)
   real(r8) :: omega                      ! radian frequency
   real(r8) :: lam                        ! wavelength (m)
   real(r8) :: k                          ! wave number (rad/m)
   real(r8) :: kcm                        ! wave number (rad/cm)
   real(r8) :: kr                         ! size parameter used in calcuate single-particle albedo

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE PM_RTM ( &
      patchtype, patchclass, dz_sno, &
      forc_topo, htop, &
      tref, t_soisno, tleaf, &
      wliq_soisno, wice_soisno, &
      snowdp, lai, sai, &
      wf_clay, wf_sand, wf_silt, BD_all, porsl, &
      sat_theta, sat_fghz, &
      tb_toa_h)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Forward modeling of brightness temperature observations
!-----------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      integer, intent(in)   :: patchtype                        ! land cover type
      integer, intent(in)   :: patchclass                       ! land cover class
      real(r8), intent(in)  :: dz_sno(maxsnl + 1:0)             ! layer thickness (m)
      real(r8), intent(in)  :: forc_topo                        ! topography [m]
      real(r8), intent(in)  :: htop                             ! upper height of vegetation [m]
      real(r8), intent(in)  :: tref                             ! 2 m height air temperature [kelvin]
      real(r8), intent(in)  :: tleaf                            ! leaf temperature [K]
      real(r8), intent(in)  :: t_soisno(maxsnl + 1:nl_soil)     ! soil temperature [K]
      real(r8), intent(in)  :: wliq_soisno(maxsnl + 1:nl_soil)  ! liquid water in layers [kg/m2]
      real(r8), intent(in)  :: wice_soisno(maxsnl + 1:nl_soil)  ! ice lens in layers [kg/m2]
      real(r8), intent(in)  :: snowdp                           ! snow depth [meter]
      real(r8), intent(in)  :: lai                              ! leaf area index
      real(r8), intent(in)  :: sai                              ! stem area index
      real(r8), intent(in)  :: wf_clay(nl_soil)                 ! gravimetric fraction of clay
      real(r8), intent(in)  :: wf_sand(nl_soil)                 ! gravimetric fraction of sand
      real(r8), intent(in)  :: wf_silt(nl_soil)                 ! gravimetric fraction of silt
      real(r8), intent(in)  :: BD_all(nl_soil)                  ! bulk density of soil (GRAVELS + ORGANIC MATTER + Mineral Soils,kg/m3)
      real(r8), intent(in)  :: porsl(nl_soil)                   ! fraction of soil that is voids [-]
      real(r8), intent(in)  :: sat_theta                        ! incidence angle of satellite (rad)
      real(r8), intent(in)  :: sat_fghz                         ! frequency of satellite (GHz)
      real(r8), intent(out) :: tb_toa_h                         ! brightness temperature of top-of-atmosphere for H- polarization
      !real(r8), intent(out) :: tb_toa_v                         ! brightness temperature of top-of-atmosphere for V- polarization

!----------------------- Local Variables -------------------------------
      logical  :: is_low_veg                  ! flag for low vegetation
      real(r8) :: dz_soisno(maxsnl+1:nl_soil) ! liquid water in layers [kg/m2]
      integer  :: lb                          ! lower bound of arrays
      real(r8) :: tau_atm                     ! atmospheric optical depth
      real(r8) :: r_r(2)                      ! rough surface reflectivity for H and V polarizations
      real(r8) :: r_sn(2)                     ! reflectivity between the snow and ground for H and V polarizations
      real(r8) :: r_snow(2)                   ! reflectivity of the snow for H and V polarizations
      real(r8) :: tb_soil(2)                  ! brightness temperature of soil for H and V polarizations
      real(r8) :: tb_tos(2)                   ! brightness temperature of snow-covered ground for H and V polarizations
      real(r8) :: tb_tov(2)                   ! brightness temperature of vegetation (consider snow) for H and V polarizations
      real(r8) :: tb_tov_noad(2)              ! brightness temperature of vegetation (no downwelling radiation) for H and V polarizations
      real(r8) :: tb_au(2)                    ! upwelling radiation (brightness temperature) of atmosphere
      real(r8) :: tb_ad(2)                    ! downwelling radiation (brightness temperature) of atmosphere
      real(r8) :: rho_snow                    ! snow density (g/cm3)
      real(r8) :: liq_snow                    ! snow liquid water content (cm3/cm3)
      real(r8) :: gamma_p(2)                  ! vegetation opacity for H- and V- polarization
      real(r8) :: tb_veg(2)                   ! brightness temperature of vegetation for H- and V- polarization
      real(r8) :: tb_2(2)                     ! the downwelling vegetation emission reflected by the soil and attenuated by the canopy layer
      real(r8) :: tb_3(2)                     ! upwelling soil emission attenuated by the canopy
      real(r8) :: tb_4(2)                     ! the downwelling cosmic ray reflected by the soil and attenuated by the canopy layer
      real(r8) :: tb_toa(2)                   ! brightness temperature of top-of-atmosphere for H- and V- polarization
      real(r8) :: wf_total(nl_soil)           ! total gravimetric
      real(r8) :: BD_all_surf                 ! bulk density of soil (g/m3) at surface
      real(r8) :: porsl_surf                  ! soil porosity at surface
      real(r8) :: t_surf                      ! soil temperature at surface (C)
      real(r8) :: t_deep                      ! soil temperature at deep layer (C)
      real(r8) :: liq_surf                    ! liquid volumetric water content at surface (m3/m3)
      real(r8) :: ice_surf                    ! ice volumetric water content at surface (m3/m3)
      real(r8) :: wf_clay_surf                ! gravimetric clay percent fraction(%) at surface
      real(r8) :: wf_sand_surf                ! gravimetric sand percent fraction(%) at surface
      integer  :: i

!-----------------------------------------------------------------------

!#############################################################################
! Prepare parameters & states used in the operator
!#############################################################################
      ! get depth of soil and snow layers
      dz_soisno(maxsnl+1:0) = dz_sno(maxsnl+1:0)
      dz_soisno(1:nl_soil) = dz_soi(1:nl_soil)

      ! calculate weighted parameters
      wf_total     = wf_clay + wf_sand + wf_silt
      wf_clay_surf = (wf_clay(1)/wf_total(1)*0.0175 + wf_clay(2)/wf_total(2)*0.0276)/0.0451*100
      wf_sand_surf = (wf_sand(1)/wf_total(1)*0.0175 + wf_sand(2)/wf_total(2)*0.0276)/0.0451*100
      BD_all_surf  = (BD_all(1)*0.0175 + BD_all(2)*0.0276)/0.0451/1000
      porsl_surf   = (porsl(1)*0.0175 + porsl(2)*0.0276)/0.0451

      ! caculate temperature (C) at surface and deep soil layers
      t_surf = ((t_soisno(1)*(0.0175) + t_soisno(2)*(0.0451 - 0.0175))/0.0451) - tfrz
      t_deep = ((t_soisno(7)*(0.8289-0.5) + t_soisno(8)*(1.0 - 0.8289))/0.5) - tfrz

      ! caculate liquid/ice volumetric water (first two layers)
      liq_surf = (wliq_soisno(1) + wliq_soisno(2))/(0.0451*denh2o)
      ice_surf = (wice_soisno(1) + wice_soisno(2))/(0.0451*denice)

      ! calculate lower bound of snow
      lb = 0
      DO i = maxsnl+1, 0
         IF (wliq_soisno(i) + wice_soisno(i) > 0.0) THEN
            lb = i
            EXIT
         ENDIF
      ENDDO

!#############################################################################
! Run the forward operator
!#############################################################################
      ! check the patch type
      IF (patchtype >= 3) THEN ! ocean, lake, ice
         tb_toa = spval
      ELSE
         ! calculate parameters used in operator varied with satellite
         CALL calc_parameters(sat_theta, sat_fghz)

!#############################################################################
! atmosphere module
!#############################################################################
         CALL atm(forc_topo, tref, tau_atm, tb_au, tb_ad)

!#############################################################################
! soil module
!#############################################################################
         CALL soil(&
            patchclass, &
            t_surf, t_deep, &
            liq_surf, ice_surf, &
            wf_sand_surf, wf_clay_surf, BD_all_surf, porsl_surf, &
            r_r, tb_soil)

!#############################################################################
! vegetation and snow module
!    We categorized four different cases for the calculations:
!    1) no vegetation and no snow
!    2) no vegetation with snow
!    3) vegetation without snow
!    4) vegetation with snow
!#############################################################################
         ! roughly judge low or high vegetation (only for IGBP)
         is_low_veg = .true.
         IF (patchclass >= 1 .and. patchclass <= 5) THEN
            is_low_veg = .false.
         END IF

         ! ensure snow density to <= 1 g/cm3
         IF (snowdp > 0.01) THEN
            rho_snow = (wliq_soisno(lb) + wice_soisno(lb))/(dz_soisno(lb)*1e3)
            liq_snow = wliq_soisno(lb)*rho_snow/(dz_soisno(lb)*1e3)

            IF (liq_snow > 1.0 .or. rho_snow > 1.0) then
               rho_snow = 1.0
               liq_snow = wliq_soisno(lb)*rho_snow/(dz_soisno(lb)*1e3)
            END IF
         END IF

         ! main procedures
         ! --------------------------------------------------------------------
         ! 1) no veg and no snow
         !    two components:
         !       (1) brightness temperature of soil
         !       (2) the downwelling radiation reflected by the soil
         ! --------------------------------------------------------------------
         IF ((lai + sai < 1e-6) .and. (snowdp < 0.01)) THEN
            tb_tov = tb_soil + tb_ad*r_r
            tb_tov_noad = tb_soil

         ! --------------------------------------------------------------------
         ! 2) no veg and has snow
         !    two components:
         !       (1) brightness temperature of snow
         !       (2) the downwelling radiation reflected by the snow
         ! --------------------------------------------------------------------
         ELSE IF ((lai + sai < 1e-6) .and. (snowdp > 0.01)) THEN
            ! calculate brightness temperature of snow-covered ground
            CALL snow(t_soisno(1), t_soisno(1), snowdp, rho_snow, liq_snow, r_r, r_snow, tb_tos)

            tb_tov = tb_tos + tb_ad*r_snow
            tb_tov_noad = tb_tos

         ! --------------------------------------------------------------------
         ! 3) has veg and no snow
         !    four components:
         !       (1) the direct upwelling vegetation emission,
         !       (2) the downwelling vegetation emission reflected by the soil and attenuated by the canopy layer
         !       (3) upwelling soil emission attenuated by the canopy
         !       (4) the downwelling reflected by the soil and attenuated by the canopy layer
         ! --------------------------------------------------------------------
         ELSE IF ((lai + sai > 1e-6) .and. (snowdp < 0.01)) THEN
            ! calculate brightness temperature of vegetation
            CALL veg(patchclass, lai, htop, 0.0, tleaf, tb_veg, gamma_p)

            DO i = 1, 2
               tb_2(i) = tb_veg(i)*gamma_p(i)*r_r(i)
               tb_3(i) = tb_soil(i)*gamma_p(i)
               tb_4(i) = tb_ad(i)*r_r(i)*(gamma_p(i)**2)
               tb_tov(i) = tb_veg(i) + tb_2(i) + tb_3(i) + tb_4(i)
               tb_tov_noad(i) = tb_veg(i) + tb_2(i) + tb_3(i)
            END DO

         ! --------------------------------------------------------------------
         ! 4) has veg and has snow
         !    We need to determine the positional relationship between vegetation and snow.
         !
         !    If vegetation is higher than snow,
         !    we first calculate brightness temperature of snow (soil boundary), then calculate
         !    four components to derive brightness temperature of top of vegetation:
         !       (1) the direct upwelling vegetation emission,
         !       (2) the downwelling vegetation emission reflected by the snow and attenuated by the canopy layer
         !       (3) upwelling snow emission attenuated by the canopy
         !       (4) the downwelling reflected by the snow and attenuated by the canopy layer
         !
         !    If vegetation is lower than snow
         !    we first calculate brightness temperature of top of vegetation (soil boundary), then calculate
         !    four components to derive brightness temperature of top of snow:
         !       (1) the direct upwelling vegetation emission,
         !       (2) the downwelling vegetation emission reflected by the soil and attenuated by the canopy layer
         !       (3) upwelling soil emission attenuated by the canopy
         !       (4) the downwelling reflected by the soil and attenuated by the canopy layer
         ! --------------------------------------------------------------------
         ELSE IF ((lai + sai > 1e-6) .and. (snowdp > 0.01)) THEN
            IF (htop < snowdp) THEN
               ! calculate brightness temperature of low vegetation
               CALL veg(patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)

               ! calculate brightness temperature of top of vegetation
               DO i = 1, 2
                  tb_2(i) = tb_veg(i)*gamma_p(i)*r_r(i)
                  tb_3(i) = tb_soil(i)*gamma_p(i)
                  tb_4(i) = tb_ad(i)*r_r(i)*(gamma_p(i)**2)
                  tb_tov(i) = tb_veg(i) + tb_2(i) + tb_3(i) + tb_4(i)
                  tb_tov_noad(i) = tb_veg(i) + tb_2(i) + tb_3(i)
               END DO

               ! calculate reflectivity between the snow and low veg (adopted from CMEM)
               r_sn(:) = 1.0 - tb_tov_noad(:)/t_soisno(1)

               ! calculate brightness temperature of snow-covered ground
               CALL snow(t_soisno(1), t_soisno(1), snowdp, rho_snow, liq_snow, r_sn, r_snow, tb_tos)

               ! calculate brightness temperature of top of snow
               tb_tov = tb_tos + tb_ad*r_snow
               tb_tov_noad = tb_tos

            ELSE
               ! calculate brightness temperature of snow-covered ground
               CALL snow(t_soisno(1), t_soisno(1), snowdp, rho_snow, liq_snow, r_r, r_snow, tb_tos)

               ! calculate brightness temperature of high vegetation
               CALL veg(patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)

               ! calculate brightness temperature of top of vegetation
               DO i = 1, 2
                  tb_2(i) = tb_veg(i)*gamma_p(i)*r_snow(i)
                  tb_3(i) = tb_tos(i)*gamma_p(i)
                  tb_4(i) = tb_ad(i)*r_snow(i)*(gamma_p(i)**2)
                  tb_tov(i) = tb_veg(i) + tb_2(i) + tb_3(i) + tb_4(i)
                  tb_tov_noad(i) = tb_veg(i) + tb_2(i) + tb_3(i)
               END DO
            END IF
         END IF

!#############################################################################
! Caculate brightness temperature of top-of-atmosphere
!#############################################################################
         tb_toa = tb_tov*exp(-tau_atm) + tb_au

      END IF

      tb_toa_h = tb_toa(1)
      !tb_toa_v = tb_toa(2)

   END SUBROUTINE PM_RTM


!-----------------------------------------------------------------------

   SUBROUTINE calc_parameters (sat_theta, sat_fghz)

!-----------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Argument -------------------------------
      real(r8), intent(in) :: sat_theta, sat_fghz

!----------------------- Local Variables -------------------------------

      theta = sat_theta                   ! incidence angle of satellite (rad)
      fghz = sat_fghz                     ! frequency of satellite (GHz)
      f = fghz*1e9                        ! frequency (Hz)
      omega = 2.0*pi*f                    ! radian frequency (rad/s)
      lam = C/f                           ! wavelength (m)
      k = 2*pi/lam                        ! wave number (rad/m)
      kcm = k/100.0                       ! wave number (rad/cm)
      kr = k*(0.5*1e-3)                   ! size parameter used in calcuate single-particle albedo

   END SUBROUTINE calc_parameters

!-----------------------------------------------------------------------

   SUBROUTINE atm(z, tref, tau_atm, tb_au, tb_ad)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the atmospheric opacity and up and downwelling brightness temperature
!
! REFERENCES:
!   [1] Pellarin, T., et al. (2003), Two-year global simulation of L-band brightness
!       temperature over land, IEEE Trans. Geosci. Remote Sens., 41, 2135–2139.
!-----------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      real(r8), intent(in)  :: z            ! altitude (m)
      real(r8), intent(in)  :: tref         ! 2m air temperature (K)
      real(r8), intent(out) :: tau_atm      ! atmospheric optical depth
      real(r8), intent(out) :: tb_au(2)     ! upwelling radiation (brightness temperature) of atmosphere
      real(r8), intent(out) :: tb_ad(2)     ! downwelling radiation (brightness temperature) of atmosphere

!----------------------- Local Variables -------------------------------
      real(r8) :: t_sky = 2.7    ! cosmic ray radiation (K)
      real(r8) :: t_eq           ! equivalent layer temperature
      real(r8) :: gossat

!-----------------------------------------------------------------------

      ! calculate optical depth of atmosphere [1] eq(A1)
      tau_atm = exp(-3.9262 - 0.2211*z/1000 - 0.00369*tref)/cos(theta)
      gossat = exp(-tau_atm)

      ! calculate equivalent layer temperature
      t_eq = exp(4.9274 + 0.002195*tref)

      ! upwelling radiation (brightness temperature) of atmosphere
      tb_au(:) = t_eq*(1.-gossat)

      ! downwelling radiation (brightness temperature) of atmosphere [1] eq(A2)
      tb_ad(:) = t_eq*(1.-gossat) + t_sky*gossat

   END SUBROUTINE atm

!-----------------------------------------------------------------------

   SUBROUTINE soil( &
      patchclass, &
      t_surf, t_deep, &
      liq_surf, ice_surf, &
      wf_sand_surf, wf_clay_surf, BD_all_surf, porsl_surf, &
      r_r, tb_soil)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate brightness temperature of soil surface
!
! REFERENCES:
!   [1] Wigneron et al., 2007, "L-band Microwave Emission of the Biosphere (L-MEB) Model:
!       Description and calibration against experimental
!       data sets over crop fields" Remote Sensing of Environment. Vol. 107, pp. 639-655k
!-----------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      integer, intent(in)   :: patchclass               ! land cover class
      real(r8), intent(in)  :: t_surf                   ! soil temperature at surface (C)
      real(r8), intent(in)  :: t_deep                   ! soil temperature at deep layer (C)
      real(r8), intent(in)  :: liq_surf                 ! liquid volumetric water content at surface (m3/m3)
      real(r8), intent(in)  :: ice_surf                 ! ice volumetric water content at surface (m3/m3)
      real(r8), intent(in)  :: wf_sand_surf             ! gravimetric sand percent fraction(%) at surface
      real(r8), intent(in)  :: wf_clay_surf             ! gravimetric clay percent fraction(%) at surface
      real(r8), intent(in)  :: BD_all_surf              ! bulk density of soil (g/m3) at surface
      real(r8), intent(in)  :: porsl_surf               ! soil porosity at surface
      real(r8), intent(out) :: r_r(2)                   ! rough surface reflectivity for H and V polarizations
      real(r8), intent(out) :: tb_soil(2)               ! brightness temperature of soil

!----------------------- Local Variables -------------------------------
      real(r8)    :: t_eff(2)                   ! effective temperature for H and V polarizations, [K]
      complex(r8) :: eps_soil                   ! dielectric constant of soil for H and V polarizations
      real(r8)    :: r_s(2)                     ! smooth surface reflectivity for H and V polarizations
      complex(r8) :: ew                         ! dielectric constant of water
      logical     :: is_desert                  ! flag for desert soil
      real(r8)    :: ffrz                       ! fraction of frozen soil
      complex(r8) :: eps_f = (5.0, 0.5)         ! dielectric constant of frozen soil
      real(r8)    :: sal_soil = 0.0             ! soil salinity (psu)

!-----------------------------------------------------------------------

      ! whether this patch is desert
      is_desert = .false.
      IF (liq_surf < 0.02 .and. wf_sand_surf > 90) THEN
         is_desert = .true.
      END IF

      ! calculate ratio of freezed soil
      IF (liq_surf + ice_surf <= 0.0d0) THEN
        ffrz = 0.0d0
      ELSE
        ffrz = ice_surf / (liq_surf + ice_surf)
      ENDIF

      ! caculate effective temperature
      CALL eff_soil_temp(liq_surf, t_surf, t_deep, t_eff)

      ! caculate dielectric constant of soil (mixture medium)
      IF (is_desert) THEN
         ! Microwave 1-10GHz permittivity of dry sand (matzler '98, eq.1)
         eps_soil = 2.53 + (2.79 - 2.53)/(1 - jj*(fghz/0.27)) + jj*0.002
      ELSE
         ! define bulk density and porosity (CMEM)
         ! BD_all_surf = (wf_sand_surf*1.60d0 + wf_clay_surf*1.10d0 + (100.0d0 - wf_sand_surf - wf_clay_surf)*1.20d0)/100.0d0
         ! porsl_surf = 1.0d0 - BD_all_surf/2.660d0

         ! caculate ice or water dielectric constant
         IF (ffrz > 0.95) THEN
            CALL diel_ice(t_surf, ew)
         ELSE
            CALL diel_water(2, liq_surf, t_surf, wf_sand_surf, wf_clay_surf, BD_all_surf, sal_soil, ew)
         END IF

         ! caculate dielectric constant in mixed soil
         IF (DEF_DA_OPERATOR_RTM_DIEL == 0) THEN
            CALL diel_soil_W80 (ew, t_surf, liq_surf, wf_sand_surf, wf_clay_surf, porsl_surf, eps_soil)
         ELSE IF (DEF_DA_OPERATOR_RTM_DIEL == 1) THEN
            CALL diel_soil_D85 (ew, liq_surf, wf_sand_surf, wf_clay_surf, BD_all_surf, eps_soil)
         ELSE IF (DEF_DA_OPERATOR_RTM_DIEL == 2) THEN
            CALL diel_soil_M04 (liq_surf, wf_clay_surf, eps_soil)
         ELSE IF (DEF_DA_OPERATOR_RTM_DIEL == 3) THEN
            CALL diel_soil_M09 (liq_surf, t_surf, wf_clay_surf, eps_soil)
         ENDIF
      END IF

      ! mix dielectric constant of frozen and non-frozen soil
      eps_soil = eps_soil*(1.-ffrz) + eps_f*ffrz

      ! caculate smooth surface reflectivity
      CALL smooth_reflectivity(eps_soil, r_s)

      ! caculate rough surface reflectivity
      CALL rough_reflectivity(is_desert, patchclass, r_s, r_r)

      ! calculate brightness temperature
      IF (is_desert) THEN
         CALL desert(t_eff, r_r, eps_soil, tb_soil)
      ELSE
         tb_soil = t_eff * (1 - r_r)
      END IF

   END SUBROUTINE soil

!-----------------------------------------------------------------------

   SUBROUTINE eff_soil_temp(wc_surf, t_surf, t_deep, t_eff)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the effective temperature of soil
!
! REFERENCES:
!   [1] An improved two-layer algorithm for estimating effective soil
!       temperature in microwave radiometry using in situ temperature
!       and soil moisture measurements
!-----------------------------------------------------------------------
      IMPLICIT NONE

! ------------------------ Dummy Argument ------------------------------
      real(r8), intent(in)  :: wc_surf        ! soil moisture at surface (m3/m3)
      real(r8), intent(in)  :: t_surf         ! soil temperature (C) at surface
      real(r8), intent(in)  :: t_deep         ! soil temperature (C) at deep layer
      real(r8), intent(out) :: t_eff(2)       ! effective temperature for H and V polarizations, [K]

!----------------------- Local Variables -------------------------------
      real(r8) :: C             ! parameter depending mainly on frequency
                                ! and soil moisture to describe the impact of
                                ! surface temperature on the effective temperature;
                                ! soil moisture increase, C large, teff close to tsurf
                                ! soil moisture decrease, C small, tdeep impact teff more
      real(r8) :: w0 = 0.30     ! parameter
      real(r8) :: bw = 0.30     ! parameter
!-----------------------------------------------------------------------

      IF (wc_surf < 0.0) THEN
        C = 0.001
      ELSE
        C = max(0.001, (wc_surf/w0)**bw)
      ENDIF
      t_eff(:) = t_deep + (t_surf - t_deep)*C + tfrz

   END SUBROUTINE eff_soil_temp

!-----------------------------------------------------------------------

   SUBROUTINE diel_ice(t, eps_i)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate dielectric constant of pure ice
!
! REFERENCES:
!   [1] Matzler, C. (2006). Thermal Microwave Radiation: Applications
!       for Remote Sensing p456-461
!-----------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      real(r8), intent(in)  :: t             ! temperature (C)
      complex(r8), intent(out) :: eps_i      ! dielectric constant of ice water

!----------------------- Local Variables -------------------------------
      real(r8) :: betam                ! beta parameter by Mishima et al. (1983)
      real(r8) :: dbeta                ! corrected delta beta parameter
      real(r8) :: beta                 ! beta parameter
      real(r8) :: t_inv                ! modified inverse temperature
      real(r8) :: tk                   ! temperature (K)
      real(r8) :: alpha                ! alpha parameter
      real(r8) :: eps_i_r              ! real part of pure ice dielectric constant
      real(r8) :: eps_i_i              ! imaginary part of pure ice dielectric constant

!-----------------------------------------------------------------------

      ! C to K
      tk = t + tfrz

      ! eq.(5.33): calculate beta parameter by Mishima et al. (1983)
      betam = (0.0207/tk)*(exp(335./tk)/((exp(335./tk) - 1.)**2.)) + 1.16e-11*(fghz**2.)      !  [1](5.33)

      ! eq.(5.35): calculate delta beta parameter
      dbeta = exp(-10.02 + 0.0364*t)                              ! [1](5.35)

      ! eq.(5.34): calculate beta parameter
      beta = betam + dbeta                                        ! [1](5.34)

      ! eq.(5.32): calculate alpha parameter
      t_inv = 300./tk - 1                                         ! [1](p.457)
      alpha = (0.00504 + 0.0062*t_inv)*exp(-22.1*t_inv) !(GHz)    ! [1](5.32)

      ! eq.(5.30): calculate real part of pure ice dielectric constant
      eps_i_r = 3.1884 + 9.1e-4*t                                 ! [1](5.30)

      ! eq.(5.31): calculate imaginary part of pure ice dielectric constant
      eps_i_i = alpha/fghz + beta*fghz                            ! [1](5.31)

      ! calculate dielectric constant of pure ice
      eps_i = eps_i_r - jj*eps_i_i                                ! [1](5.31)

   END SUBROUTINE diel_ice

!-----------------------------------------------------------------------

   SUBROUTINE diel_water(type, swc, t, wf_sand, wf_clay, BD_all, sal, eps_w)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate dielectric constant of water in water (saline water)
!
! REFERENCES:
!   [1] Ulaby FT, R. K. Moore, and A. K. Fung, Microwave Remote Sensing:
!       Active and Passive. Vol. III. From theory to applications. Artech House,
!       Norwood, MA., 1986
!   [2] Klein, L. A. and C. T. Swift (1977): An improved model
!       for the dielectric constant of sea water at microwave
!       frequencies, IEEE Transactions on  Antennas and Propagation,
!       Vol. AP-25, No. 1, 104-111.
!-----------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      integer, intent(in)      :: type       ! type of water, 0: pure water, 1: sea water, 2: soil water
      real(r8), intent(in)     :: swc        ! soil water content (m3/m3)
      real(r8), intent(in)     :: t          ! soil temperature (C)
      real(r8), intent(in)     :: wf_sand    ! gravimetric sand percent fraction(%)
      real(r8), intent(in)     :: wf_clay    ! gravimetric clay percent fraction(%)
      real(r8), intent(in)     :: BD_all     ! bulk density(g/cm3)
      real(r8), intent(in)     :: sal        ! water salinity (psu)
      complex(r8), intent(out) :: eps_w      ! dielectric constant of soil water

!----------------------- Local Variables -------------------------------
      real(r8) :: sigma               ! ionic conductivity (S/m)
      real(r8) :: a, b                ! parameters
      real(r8) :: tau_w               ! relaxation time of pure water
      real(r8) :: eps_w0              ! static dielectric constant of pure water
      real(r8) :: wc

!-----------------------------------------------------------------------

      ! [3] eq.16: tau(T, sal) = tau_w(T) * b(sal, T)
      ! calculate relaxation time of pure water (Stogryn)
      tau_w = 1.768e-11 - 6.068e-13*t + 1.104e-14*t**2 - 8.111e-17*t**3               ! [2](17)
      b = 1.000 + 2.282e-5*sal*t - 7.638e-4*sal - 7.760e-6*sal**2 + 1.105e-8*sal**3   ! [2](18)
      tau_w = tau_w*b                                                                 ! [2](16)

      ! [3] eq.13: eps_w0(sal, T) = eps_w0(T) * a(sal, T)
      ! static dielectric constant of pure water (Klein and Swift)
      eps_w0 = 87.134 - 1.949e-1*t - 1.276e-2*t**2 + 2.491e-4*t**3                    ! [2](14)
      a = 1.000 + 1.613e-5*sal*t - 3.656e-3*sal + 3.210e-5*sal**2 - 4.232e-7*sal**3   ! [2](15)
      eps_w0 = eps_w0*a                                                               ! [2](13)

      IF (type == 0) THEN  ! pure water
         ! [1] eq.19
         eps_w0 = 88.045 - 0.4147*t + 6.295e-4*t**2 + 1.075e-5*t**3
         eps_w = eps_w_inf + (eps_w0 - eps_w_inf)/(1 - jj*omega*tau_w)

      ELSEIF (type == 1) THEN  ! sea water
         ! calculate ionic conductivity [1] eq.27, eq.28
         sigma = sal*(0.182521 - 1.46192e-3*sal + 2.09324e-5*sal**2 - 1.28205e-7*sal**3) &
            *exp(-1.*(25 - t)* &
            (2.033e-2 + 1.266e-4*(25 - t) + 2.464e-6*(25 - t)**2 &
            - sal*(1.849e-5 - 2.551e-7*(25 - t) + 2.551e-8*(25 - t)**2)))

         ! diel constant of sea water [1] eq.21
         eps_w = eps_w_inf + (eps_w0 - eps_w_inf)/(1 - jj*omega*tau_w) + jj*sigma/(omega*eps_0)
      ELSEIF (type == 2) THEN
         ! calculate soil conductivity
         sigma = -1.645 + 1.939*BD_all - 0.02256*wf_sand + 0.01594*wf_clay
         IF (sigma < 0.) THEN
            sigma = 0. ! negative for very sandy soils with low bulk density
         END IF

         ! calculate dielectric constant of soil-water by modified Debye expression
         wc = max(0.001, swc)
         eps_w = eps_w_inf + (eps_w0 - eps_w_inf)/(1 - jj*omega*tau_w) &
            + jj*sigma/(omega*eps_0)*(rho_soil - BD_all)/(rho_soil*wc)
      END IF

   END SUBROUTINE diel_water

!-----------------------------------------------------------------------

   SUBROUTINE diel_soil_W80(ew, t, wc, wf_sand, wf_clay, porsl, eps)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the dielectric constant of a wet soil
!
! REFERENCES:
!   [1]  Matzler, C. (1998). Microwave permittivity of dry sand.
!   IEEE Transactions on Geoscience and Remote Sensing, 36(1), 317-319.
!
!   [2]  Wang and Schmugge, 1980: An empirical model for the
!   complex dielectric permittivity of soils as a function of water
!   content. IEEE Trans. Geosci. Rem. Sens., GE-18, No. 4, 288-295.
!-----------------------------------------------------------------------
      IMPLICIT NONE

! ------------------------ Dummy Argument ------------------------------
      real(r8), intent(in)  :: t                        ! soil temperature (C)
      real(r8), intent(in)  :: wc                       ! volumetric soil moisture (m3/m3)
      real(r8), intent(in)  :: wf_sand                  ! gravimetric sand percent fraction(%)
      real(r8), intent(in)  :: wf_clay                  ! gravimetric clay percent fraction(%)
      real(r8), intent(in)  :: porsl                    ! soil porosity at surface
      complex(r8), intent(in)  :: ew                    ! dielectric constant of water
      complex(r8), intent(out) :: eps

!----------------------- Local Variables -------------------------------
      real(r8) :: wp                       ! wilting point
      real(r8) :: wt                       ! transition moisture point (cm3/cm3)
      real(r8) :: gamma                    ! fitting parameter
      real(r8) :: ecl                      ! conductivity loss
      real(r8) :: alpha                    ! conductivity loss parameter
      real(r8) :: sal_sea = 32.5           ! sea water salinity (psu)
      real(r8) :: sal_soil = 0.0           ! soil salinity (psu)
      complex(r8) :: eps_x                 ! dielectric constant of the initially absorbed water
      complex(r8) :: eps_a = (1.0, 0.0)    ! dielectric constant of air,  [2]IV
      complex(r8) :: eps_r = (5.5, 0.2)    ! dielectric constant of rock, [2]IV
      complex(r8) :: eps_i = (3.2, 0.1)    ! dielectric constant of ice,  [2]IV
      complex(r8) :: eps_f = (5.0, 0.5)    ! dielectric constant of frozen soil

!-----------------------------------------------------------------------

      ! calculate wilting point at the soil layer
      wp = 0.06774 - 0.00064*wf_sand + 0.00478*wf_clay        ! [2](1)

      ! calculate fitting parameters
      gamma = -0.57*wp + 0.481                                ! [2](8)

      ! calculate transition moisture point
      wt = 0.49*wp + 0.165                                    ! [2](9)

      ! calculate dielectric constant of wet soil (when all soil freeze, eps_x = eps_i)
      IF (wc <= wt) THEN
        eps_x = eps_i + (ew - eps_i)*(wc/wt)*gamma                                ! [2](3)
        eps = wc*eps_x + (porsl - wc)*eps_a + (1.-porsl)*eps_r                    ! [2](2)
      ELSE
        eps_x = eps_i + (ew - eps_i)*gamma                                        ! [2](5)
        eps = wt*eps_x + (wc - wt)*ew + (porsl - wc)*eps_a + (1.-porsl)*eps_r     ! [2](4)
      END IF

      ! add conductivity loss for imaginary part
      alpha = min(100.*wp, 26.)
      ecl = alpha*wc**2                          ! [2](6)
      eps = eps + jj*ecl                         ! [2](6)

   END SUBROUTINE diel_soil_W80

!-----------------------------------------------------------------------

   SUBROUTINE diel_soil_M04(wc, wf_clay, eps)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the dielectric constant of a wet soil Developed and
!   validated from 1 to 10 GHz, adapted for a large range of soil moisture
!
! REFERENCES:
!   [1] Mironov et al, Generalized Refractive Mixing Dielectric Model for
!       moist soil. IEEE Trans. Geosc. Rem. Sens., vol 42 (4), 773-785. 2004.
!
!   [2] Mironov et al, Physically and Mineralogically Based Spectroscopic
!       Dielectric Model for Moist Soils. IEEE Trans. Geosc. Rem. Sens.,
!       vol 47 (7), 2059-2070. 2009.
!-----------------------------------------------------------------------
      IMPLICIT NONE

! ------------------------ Dummy Argument ------------------------------
      real(r8), intent(in)  :: wc          ! soil moisture (m3/m3)
      real(r8), intent(in)  :: wf_clay     ! gravimetric clay percent fraction(%)
      complex(r8), intent(out) :: eps

!----------------------- Local Variables -------------------------------
      real(r8) :: znd, zkd, zxmvt, zep0b, ztaub, zsigmab, zep0u, ztauu
      real(r8) :: zsigmau, zcxb, zepwbx, zepwby, zcxu, zepwux, zepwuy
      real(r8) :: znb, zkb, znu, zku, zxmvt2, znm, zkm, zepmx, zepmy
      integer  :: zflag

!-----------------------------------------------------------------------
!------------------------------------------------------------------------
!  Initializing the GRMDM spectroscopic parameters with clay (fraction)
!------------------------------------------------------------------------
      ! RI & NAC of dry soils
      znd = 1.634 - 0.539 * (wf_clay/100) + 0.2748 * (wf_clay/100) ** 2
      zkd = 0.03952 - 0.04038 * (wf_clay / 100)                                    ! [2](18)

      ! Maximum bound water fraction
      zxmvt = 0.02863 + 0.30673 * wf_clay / 100                                    ! [2](19)

      ! Bound water parameters
      zep0b   = 79.8 - 85.4  * (wf_clay / 100) + 32.7  * (wf_clay / 100)*(wf_clay / 100) ! [2](20)
      ztaub   = 1.062e-11 + 3.450e-12 * (wf_clay / 100)                            ! [2](21)
      zsigmab = 0.3112 + 0.467 * (wf_clay / 100)                                   ! [2](22)

      ! Unbound (free) water parameters
      zep0u   = 100                                                                ! [2](24)
      ztauu   = 8.5e-12                                                            ! [2](25)
      zsigmau = 0.3631 + 1.217 * (wf_clay / 100)

      ! Computation of epsilon water (bound & unbound)
      zcxb   = (zep0b - eps_w_inf) / (1. + (2.*pi*f*ztaub)**2)                     ! [2](16)
      zepwbx = eps_w_inf + zcxb                                                    ! [2](16)
      zepwby = zcxb * (2.*pi*f*ztaub) + zsigmab / (2.*pi*eps_0*f)                  ! [2](16)
      zcxu   = (zep0u - eps_w_inf) / (1 + (2*pi*f*ztauu)**2)                       ! [2](16)
      zepwux = eps_w_inf + zcxu                                                    ! [2](16)
      zepwuy = zcxu * (2.*pi*f*ztauu) + zsigmau/(2.*pi*eps_0*f)

      ! Computation of refractive index of water (bound & unbound)
      znb = sqrt( sqrt( zepwbx**2 + zepwby**2) + zepwbx ) / sqrt(2.0)              ! [2](14)
      zkb = sqrt( sqrt( zepwbx**2 + zepwby**2) - zepwbx ) / sqrt(2.0)              ! [2](15)
      znu = sqrt( sqrt( zepwux**2 + zepwuy**2) + zepwux ) / sqrt(2.0)              ! [2](14)
      zku = sqrt( sqrt( zepwux**2 + zepwuy**2) - zepwux ) / sqrt(2.0)              ! [2](15)

      ! Computation of soil refractive index (nm & km): xmv can be a vector
      zxmvt2 = min (wc, zxmvt)
      zflag  = 0
      IF ( wc >= zxmvt ) zflag = 1
      znm = znd + (znb - 1) * zxmvt2 + (znu - 1) * (wc-zxmvt) * zflag              ! [2](12)
      zkm = zkd + zkb * zxmvt2 + zku * (wc-zxmvt) * zflag                          ! [2](13)

      ! computation of soil dielectric constant:
      zepmx = znm ** 2 - zkm ** 2                                                  ! [2](11)
      zepmy = znm * zkm * 2                                                        ! [2](11)
      eps   = cmplx(zepmx, zepmy, kind=r8)

   END SUBROUTINE diel_soil_M04

!-----------------------------------------------------------------------

   SUBROUTINE diel_soil_M09(wc, t, wf_clay, eps)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the dielectric constant of a wet soil Developed and
!   validated from 1 to 10 GHz, adapted for a large range of soil moisture
!
! REFERENCES:
!   [1] V. L. Mironov, S. V. Fomin,
!       "Temperature and mineralogy dependable model for microwave dielectric
!       spectra of moist soils", PIERS Online, vol. 5, no. 5, pp. 411-415, 2009.
!
!   [2] Mironov et al, Physically and Mineralogically Based Spectroscopic Dielectric
!       Model for Moist Soils. IEEE Trans. Geosc. Rem. Sens., vol 47 (7), 2059-2070. 2009.
!-----------------------------------------------------------------------
      IMPLICIT NONE

! ------------------------ Dummy Argument ------------------------------
      real(r8), intent(in)  :: t               ! soil temperature
      real(r8), intent(in)  :: wc              ! soil moisture (m3/m3)
      real(r8), intent(in)  :: wf_clay         ! weighted fraction (%)
      complex(r8), intent(out) :: eps

!----------------------- Local Variables -------------------------------
      real(r8) :: nd, kd, mvt, ts, e0b, Bb, Bsgb, Fb, eb0
      real(r8) :: dSbR, taub, sigmabt, sigmab, e0u, Bu, Bsgu, Fu, eu0, dHuR, dSuR, tauu, dHbR
      real(r8) :: sigmau, sigmaut, cxb, eb_r, eb_i, cxu, eu_r, eu_i, nb, kb, nu, ku, nm, km, eps_r, eps_i

!-----------------------------------------------------------------------
!------------------------------------------------------------------------
!  Initializing the GRMDM spectroscopic parameters with clay (fraction)
!------------------------------------------------------------------------
      ! RI & NAC of dry soils
      nd = 1.634 - 0.539e-2 * wf_clay + 0.2748e-4 * (wf_clay ** 2)   ! [1](11)
      kd = 0.03952 - 0.04038e-2 * wf_clay                            ! [1](12)

      ! maximum bound water fraction
      mvt = 0.02863 + 0.30673e-2 * wf_clay                           ! [1](13)

      ! starting temperature for parameters' fit ([1] p.413)
      ts = 20.

      ! eb0 computation
      e0b  = 79.8 - 85.4e-2 * wf_clay + 32.7e-4 * (wf_clay **2)      ! [1](14)                                                               ! [1](14)
      Bb   = 8.67e-19 - 0.00126e-2 * wf_clay + 0.00184e-4 * (wf_clay ** 2)  - 9.77e-10*(wf_clay**3) - 1.39e-15 *(wf_clay**4)   ! [1](15)
      Bsgb = 0.0028  + 0.02094e-2*wf_clay - 0.01229e-4*(wf_clay**2) - 5.03e-22*(wf_clay**3) + 4.163e-24*(wf_clay**4)           ! [1](23)
      Fb   = log((e0b - 1)/(e0b + 2))                                                        ! [1](8)(ep0->e0p)
      eb0  = (1 + 2*exp(Fb-Bb*(t-ts))) / (1 - exp(Fb-Bb*(t-ts)))                             ! [1](7)(e0p->ep0)

      ! taub computation
      dHbR = 1467 + 2697e-2*wf_clay - 980e-4 *(wf_clay**2) + 1.368e-10*(wf_clay**3) - 8.61e-13 *(wf_clay**4)         ! [1](18)
      dSbR = 0.888 + 9.7e-2 *wf_clay - 4.262e-4*(wf_clay**2) + 6.79e-21 *(wf_clay**3) + 4.263e-22*(wf_clay**4)       ! [1](19)
      taub = 48e-12 * exp(dHbR/(t+tfrz)-dSbR)/(t+tfrz)                                      ! [1](9)

      ! sigmab computation
      sigmabt = 0.3112 + 0.467e-2*wf_clay                                                   ! [1](22)
      sigmab  = sigmabt + Bsgb*(t-ts)                                                       ! [1](10)
  
      ! unbound (free) water parameters
      !-------------------
      !  eu0 computation
      !-------------------
      e0u  = 100.                                                                                                    ! [1](16)
      Bu   = 1.11e-4 - 1.603e-7 *wf_clay + 1.239e-9 *(wf_clay**2) + 8.33e-13 *(wf_clay**3) - 1.007e-14*(wf_clay**4)  ! [1](17)
      Bsgu = 0.00108 + 0.1413e-2*wf_clay - 0.2555e-4*(wf_clay**2) + 0.2147e-6*(wf_clay**3) - 0.0711e-8*(wf_clay**4)  ! [1](25)
      Fu   = log((e0u - 1)/(e0u + 2))                                                                                ! [1](8)(ep0->e0p)
      eu0  = (1 + 2*exp(Fu-Bu*(t-ts))) / (1-exp(Fu-Bu*(t-ts)))                                                       ! [1](7))e0p->ep0)

      !--------------------
      !  tauu computation
      !--------------------
      dHuR = 2231 - 143.1e-2 *wf_clay + 223.2e-4*(wf_clay**2) - 142.1e-6*(wf_clay**3) + 27.14e-8 *(wf_clay**4)       ! [1](20)
      dSuR = 3.649 - 0.4894e-2*wf_clay + 0.763e-4*(wf_clay**2) - 0.4859e-6*(wf_clay**3) + 0.0928e-8*(wf_clay**4)     ! [1](21)
      tauu = 48e-12 * exp(dHuR/(t+tfrz)-dSuR)/(t+tfrz)                                                               ! [1](9)

      !----------------------
      !  sigmau computation
      !----------------------
      sigmaut = 0.05 + 1.4*(1.0 - (1.0 - wf_clay*1.e-2)**4.664)                                                      ! [1](24)
      sigmau  = sigmaut + Bsgu*(t-ts)                                                                                ! [1](10)

      !--------------------------------------------------
      !  computation of epsilon water (bound & unbound)
      !--------------------------------------------------
      cxb  = (eb0-eps_w_inf) / (1+(2*pi*f*taub)**2)           ! [1](6), [2](16)
      eb_r = eps_w_inf + cxb                                  ! [1](6), [2](16)
      eb_i = cxb*(2*pi*f*taub) + sigmab/(2*pi*eps_0*f)        ! [1](6), [2](16)
      cxu  = (eu0-eps_w_inf) / (1+(2*pi*f*tauu)**2)           ! [1](6), [2](16)
      eu_r = eps_w_inf + cxu                                  ! [1](6), [2](16)
      eu_i = cxu*(2*pi*f*tauu) + sigmau/(2*pi*eps_0*f)        ! [1](6), [2](16)

      !--------------------------------------------------------------
      !  computation of refractive index of water (bound & unbound)
      !--------------------------------------------------------------
      nb = sqrt(sqrt(eb_r**2+eb_i**2)+eb_r) / sqrt(2.0)       ! [1](5)
      kb = sqrt(sqrt(eb_r**2+eb_i**2)-eb_r) / sqrt(2.0)       ! [1](5)
      nu = sqrt(sqrt(eu_r**2+eu_i**2)+eu_r) / sqrt(2.0)       ! [1](5)
      ku = sqrt(sqrt(eu_r**2+eu_i**2)-eu_r) / sqrt(2.0)       ! [1](5)

      !--------------------------------------------------
      !  computation of soil refractive index (nm & km)
      !--------------------------------------------------
      IF (wc <= mvt) THEN
         nm = nd + (nb-1)*wc                                  ! [2](12)
         km = kd + kb*wc                                      ! [2](13)
      ELSE
         nm = nd + (nb-1)*mvt + (nu-1)*(wc-mvt)               ! [2](12)
         km = kd + kb*mvt + ku*(wc-mvt)                       ! [2](13)
      ENDIF

      !-------------------------------------------
      !  computation of soil dielectric constant
      !-------------------------------------------
      eps_r = nm**2 - km**2                                   ! [1](4)
      eps_i = 2* nm * km                                      ! [1](4)
      eps   = cmplx(eps_r, eps_i, kind=r8)

   END SUBROUTINE diel_soil_M09

!-----------------------------------------------------------------------

   SUBROUTINE diel_soil_D85(ew, swc, wf_sand, wf_clay, BD_all, eps)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the dielectric constant of a wet soil Developed and validated for 1.4 and 18 GHz.

! REFERENCES:
!   [1] Dobson et al., 1985: Microwave Dielectric behavior of
!       wet soil - part II: Dielectric mixing models,
!       IEEE Trans. Geosc. Rem. Sens., GE-23, No. 1, 35-46.

!   [2] N. R. Peplinski, F. T. Ulaby, and M. C. Dobson,
!       Dielectric Properties of Soils in the 0.3-1.3-GHz Range,
!       IEEE Trans. Geosc. Rem. Sens., vol. 33, pp. 803-807, May 1995

!   [3] N. R. Peplinski, F. T. Ulaby, and M. C. Dobson,
!       Corrections to “Dielectric Properties of Soils in the 0.3-1.3-GHz Range",
!       IEEE Trans. Geosc. Rem. Sens., vol. 33, p. 1340, November 1995

!-----------------------------------------------------------------------
      IMPLICIT NONE

! ------------------------ Dummy Argument ------------------------------
      complex(r8), intent(in)  :: ew                     ! dielectric constant of water
      real(r8), intent(in)  :: swc                       ! soil moisture
      real(r8), intent(in)  :: wf_sand                   ! gravimetric sand percent fraction(%)
      real(r8), intent(in)  :: wf_clay                   ! gravimetric clay percent fraction(%)
      real(r8), intent(in)  :: BD_all                    ! soil bulk density (g/cm3)
      complex(r8), intent(out) :: eps

!----------------------- Local Variables -------------------------------
      real(r8) :: alphas = 0.65_r8
      real(r8) :: beta, eaa, eps_s, epsi, epsr, wc

!-----------------------------------------------------------------------

      wc    = max(swc, 0.001_r8)
      eps_s = (1.01_r8 + 0.44_r8 * rho_soil)**2.0_r8 - 0.062_r8                 ! [1](22)
      beta  = (127.48_r8 - 0.519_r8 * wf_sand - 0.152_r8 * wf_clay) / 100.0_r8  ! [1](30)
      eaa   = 1.0_r8 + (BD_all / rho_soil) * (eps_s ** alphas - 1.0_r8) &
      &      + (wc ** beta) * (real(ew) ** alphas) - wc                         ! [1](28)
      epsr  = eaa ** (1.0_r8/alphas)                                            ! [2](2),[3]
      beta  = (133.797_r8 - 0.603_r8 * wf_sand - 0.166_r8 * wf_clay) / 100.0_r8 ! [1](31) 1.33797 -> 133.797, [2](5)
      eaa   = (wc ** beta) * (abs(aimag(ew)) ** alphas)                         ! [2](3)
      epsi  = eaa ** (1.0_r8/alphas)                                            ! [2](3)
      eps   = cmplx(epsr, epsi, kind=r8)

   END SUBROUTINE diel_soil_D85

!-----------------------------------------------------------------------

   SUBROUTINE smooth_reflectivity(eps, r_s)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the smooth surface reflectivity by Fresnel Law
!
! REFERENCES:
!   [1] Njoku and Kong, 1977: Theory for passive microwave remote sensing
!       of near-surface soil moisture. Journal of Geophysical Research,
!       Vol. 82, No. 20, 3108-3118.
!-----------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      complex(r8), intent(in) :: eps       ! dielectric constant of the surface
      real(r8), intent(out)   :: r_s(2)    ! reflectivities of flat surfaces for H and V polarizations

!----------------------- Local Variables -------------------------------
      complex(r8) :: g                     ! parameter in Fresnel Law

!-----------------------------------------------------------------------

      g = sqrt(eps - sin(theta)**2)
      r_s(1) = abs((cos(theta) - g)/(cos(theta) + g))**2.
      r_s(2) = abs((cos(theta)*eps - g)/(cos(theta)*eps + g))**2.

   END SUBROUTINE smooth_reflectivity

!-----------------------------------------------------------------------

   SUBROUTINE rough_reflectivity(is_desert, patchclass, r_s, r_r)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the rough surface reflectivity
!
! REFERENCES:
!   [1] Kerr and Njokui, 1990: A Semiempirical Model For Interpreting Microwave
!       Emission From Semiarid Land Surfaces as Seen From Space
!       IEEE Trans. Geosci. Rem. Sens., Vol.28, No.3, 384-393.
!
!   [2] Wigneron, J. P., Jackson, T. J., O'neill, P., De Lannoy, G., de Rosnay, P., Walker,
!       J. P., ... & Kerr, Y. (2017). Modelling the passive microwave signature from land surfaces:
!       A review of recent results and application to the L-band SMOS & SMAP soil moisture retrieval algorithms.
!       Remote Sensing of Environment, 192, 238-262.
!-----------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      logical, intent(in)   :: is_desert   ! flag for desert soil
      integer, intent(in)   :: patchclass  ! patch class
      real(r8), intent(in)  :: r_s(2)      ! reflectivities of flat surfaces for H and V polarizations
      real(r8), intent(out) :: r_r(2)      ! reflectivities of rough surfaces for H and V polarizations

!----------------------- Local Variables -------------------------------
      real(r8) :: Q                          ! parameter for polarization mixing
      real(r8) :: hr(N_land_classification)  ! roughness parameter
      real(r8) :: nrh(N_land_classification) ! parameter for H polarization
      real(r8) :: nrv(N_land_classification) ! parameter for V polarization

!-----------------------------------------------------------------------

      IF (is_desert) THEN
         r_r = r_s
      ELSE
         ! calculate parameter for polarization mixing due to surface roughness
         IF (fghz < 2.) THEN
            Q = 0. ! Q is assumed zero at low frequency
         ELSE
            Q = 0.35*(1.0 - exp(-0.6*rgh_surf**2*fghz))     !    [1](16)
         END IF

         ! calculate rough surface reflectivity (default settings used in [2])
         IF (DEF_DA_OPERATOR_RTM_ROUGH == 0) THEN
            hr(:) = (2.0*kcm*rgh_surf)**2.0
            nrh(:) = 0.0
            nrv(:) = 0.0
         ELSE IF (DEF_DA_OPERATOR_RTM_ROUGH == 1) THEN
            hr(:) = hr_SMOS
            nrh(:) = 2.0
            nrv(:) = 0.0
         ELSE IF (DEF_DA_OPERATOR_RTM_ROUGH == 2) THEN
            hr(:) = hr_SMAP
            nrh(:) = 2.0
            nrv(:) = 2.0
         ELSE IF (DEF_DA_OPERATOR_RTM_ROUGH == 3) THEN
            hr(:) = hr_P16
            nrh(:) = -1.0
            nrv(:) = -1.0
         END IF

         ! rough surface reflectivity for H and V polarizations
         r_r(1) = (Q*r_s(2) + (1.-Q)*r_s(1))*exp(-hr(patchclass)*cos(theta)**nrh(patchclass))
         r_r(2) = (Q*r_s(1) + (1.-Q)*r_s(2))*exp(-hr(patchclass)*cos(theta)**nrv(patchclass))
      END IF

   END SUBROUTINE rough_reflectivity

!-----------------------------------------------------------------------

   SUBROUTINE desert(t_soil, r_r, eps, tb_desert)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate desert emissivity using Grody and Weng, 2008
!
! REFERENCES:
!  [1] Grody, N. C., & Weng, F. (2008). Microwave emission and scattering from deserts:
!      Theory compared with satellite measurements.
!      IEEE Transactions on Geoscience and Remote Sensing, 46, 361–375.
!-----------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      real(r8), intent(in)    :: t_soil(2)     ! desert surface temperature (K) (h-pol and v-pol)
      real(r8), intent(in)    :: r_r(2)        ! reflectivity of rough surface  (h-pol and v-pol)
      complex(r8), intent(in) :: eps           ! diel. const of desert
      real(r8), intent(out)   :: tb_desert(2)  ! brightness temperature of soil for H- and V- polarization

!----------------------- Local Variables -------------------------------
      real(r8) :: f0 = 0.7                      ! the fractional volume of spherical particles
                                                ! (f = (4/3)*pi*r^3*n0),
                                                ! r : the particle radius = 0.5 (mm)
                                                ! n0 : the number of particles per unit volume.
      real(r8) :: w     ! single-particle albedo
      real(r8) :: g     ! asymmetry parameter
      real(r8) :: a     ! similarity parameter
      real(r8) :: em(2) ! desert soil emissivity
      real(r8) :: y_r   ! real part of y-parameters
      real(r8) :: y_i   ! imaginary part of y-parameters

!-----------------------------------------------------------------------
      ! calculate y-parameters (eq.A15)
      y_r = (real(eps) - 1)/(real(eps) + 2)                              ! [1](A15)
      y_i = 3*aimag(eps)/(real(eps) + 2)**2                              ! [1](A15)

      ! calculate single-particle albedo (eq.A16)
      w = (1 - f0)**4*kr**3*y_r**2/ &
         ((1 - f0)**4*kr**3*y_r**2 + 1.5*(1 + 2*f0)**2*y_i)              ! [1](A16)

      ! calculate asymmetry parameter (p.374)
      g = 0.23*kr**2                                                     ! [1]p.374

      ! calculate similarity parameter (eq.3b)
      a = sqrt((1 - w)/(1 - w*g))                                        ! [1](3b)

      ! calculate desert soil emissivity (eq.A13)
      em = (1 - r_r)*(2*a/((1 + a) - (1 - a)*r_r))                       ! [1](13)

      ! calculate brightness temperature of desert
      tb_desert = t_soil*em

   END SUBROUTINE desert

!-----------------------------------------------------------------------

   SUBROUTINE veg(patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the brightness temperature on the top of vegetation
!
! REFERENCES:
!   [1] Wigneron et al., 2007, "L-band Microwave Emission of the Biosphere (L-MEB) Model:
!       Description and calibration against experimental
!       data sets over crop fields" Remote Sensing of Environment. Vol. 107, pp. 639-655k
!
!-----------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      integer, intent(in)   :: patchclass     ! land cover class
      real(r8), intent(in)  :: lai            ! leaf area index
      real(r8), intent(in)  :: tleaf          ! leaf temperature (K)
      real(r8), intent(out) :: tb_veg(2)      ! brightness temperature of vegetation for H- and V- polarization
      real(r8), intent(out) :: gamma_p(2)     ! vegetation opacity for H- and V- polarization
      real(r8), intent(in)  :: htop, snowdp

!----------------------- Local Variables -------------------------------
      real(r8) :: tau_nadir    ! vegetation opacity at nadir
      real(r8) :: tau_veg(2)   ! vegetation opacity for H- and V- polarization
      integer  :: i

!-----------------------------------------------------------------------

      ! caculate vegetation opacity (optical depth) at nadir b*VWC
      IF (htop < snowdp) THEN
         tau_nadir = b1(patchclass)*lai + b2(patchclass) ! low veg              ! [1](22)
      ELSE
         tau_nadir = b3(patchclass) ! high veg
      END IF

      ! calculate vegetation optical depth at H- and V- polarizations
      tau_veg(1) = tau_nadir*(cos(theta)**2 + tth(patchclass)*sin(theta)**2)    ! [1](23)
      tau_veg(2) = tau_nadir*(cos(theta)**2 + ttv(patchclass)*sin(theta)**2)    ! [1](24)

      ! calculate brightness temperature of vegetation
      DO i = 1, 2
         gamma_p(i) = exp(-tau_veg(i)/cos(theta))                               !  [1](15)
         tb_veg(i) = (1.-w_CMEM(patchclass))*(1.-gamma_p(i))*tleaf
      END DO

   END SUBROUTINE veg

!-----------------------------------------------------------------------

   SUBROUTINE snow(t_snow, t, snowdp, rho_snow, liq_snow, r_sn, r_snow, tb_tos)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the brightness temperature of snow-covered ground
!
! REFERENCES:
!   [1] Christian Mätzler (1987) Applications of the interaction of
!       microwaves with the natural snow cover, Remote Sensing Reviews, 2:2, 259-387, DOI:
!       10.1080/02757258709532086
!
!   [2] Anderson, E. A., 1976: A point energy and mass balance model of a snow cover.
!       NOAA Tech. Rep. NWS 19, 150 pp. U.S. Dept. of Commer., Washington, D.C.(eq.5.1)
!
!   [3] Hallikainen, M. T., F. Ulaby, and T. Deventer. 1987. Extinction behavior of dry snow in the
!       18- to 90-GHz range. IEEE Trans. Geosci. Remote Sens., GE-25, 737–745.
!
!   [4] Microwave remote sensing : active and passive
!-----------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      real(r8), intent(in)  :: t_snow       ! average snow temperature (K)
      real(r8), intent(in)  :: t            ! temperature at bottom of snow (K), i.e., soil or leaf
      real(r8), intent(in)  :: snowdp       ! snow depth (m)
      real(r8), intent(in)  :: rho_snow     ! snow density (g/cm3)
      real(r8), intent(in)  :: liq_snow     ! snow liquid water content (cm3/cm3)
      real(r8), intent(in)  :: r_sn(2)      ! reflectivity between the snow and ground at (1, H-POL. 2, V.)
      real(r8), intent(out) :: r_snow(2)    ! reflectivity between the snow and air for H- and V- polarization
      real(r8), intent(out) :: tb_tos(2)    ! brightness temperature of snow-cover ground for H- and V- polarization

!----------------------- Local Variables -------------------------------
      real(r8) :: sal_snow = 0.0                   ! snow salinity (pmm)
      real(r8) :: eps_i_r                          ! real part of dielectric constant of ice
      real(r8) :: eps_i_i                          ! imaginary part of dielectric constant of ice
      real(r8) :: eps_i_is                         ! imaginary part of dielectric constant of impure ice -5(C)
      real(r8) :: eps_i_ip                         ! imaginary part of dielectric constant of pure ice -5(C)
      real(r8) :: eps_ds_r                         ! real part of dielectric constant of dry snow
      real(r8) :: eps_ds_i                         ! imaginary part of dielectric constant of dry snow
      real(r8) :: eps_ws_i                         ! imaginary part of dielectric constant of wet snow
      real(r8) :: eps_ws_r                         ! real part of dielectric constant of wet snow
      real(r8) :: eps_w_s = 88.                    ! dielectric constant of static water
      real(r8) :: eps_a_inf, eps_b_inf, eps_c_inf  ! infinite frequency dielectric constant of three parts
      real(r8) :: eps_a_s, eps_b_s, eps_c_s        ! static dielectric constant of three parts
      complex(r8) :: eps_i                         ! dielectric constant of ice
      complex(r8) :: eps_a
      complex(r8) :: eps_b
      complex(r8) :: eps_c                         ! dielectric constant of three parts
      complex(r8) :: eps                           ! dielectric constant of wet snow
      real(r8) :: rho_ds                           ! density of dry snow  (g/cm3)
      real(r8) :: rho_i = 0.916                    ! density of ice (g/cm3)
      real(r8) :: aa = 0.005
      real(r8) :: bb = 0.4975
      real(r8) :: cc = 0.4975                      ! fitting parameters
      real(r8) :: fa, fb, fc                       ! relaxation frequency of wet snow
      real(r8) :: d                                ! snow grain size (mm)
      real(r8) :: alpha, beta, pp, qq              ! parameter used to calculate propogation angle in snow
      real(r8) :: theta_s                          ! propogation angle in snow
      complex(r8) :: z_s                           ! wave impedance in snow
      real(r8) :: r_sa                             ! reflectivity between the snow and air for H- and V- polarization
      real(r8) :: tb_2                             ! the net apparent temperature contributions due to emission by layers 2 (snow)
      real(r8) :: tb_3                             ! the net apparent temperature contributions due to emission by layers 3 (soil)
      real(r8) :: l2_apu                           ! extinction coefficient of snow (Beer's Law)
      real(r8) :: q = 0.96                         ! parameter
      real(r8) :: ka_ws                            ! absorption coefficient of wet snow
      real(r8) :: ka_ds                            ! absorption coefficient of dry snow
      real(r8) :: ke                               ! extinction coefficient of wet snow
      real(r8) :: ke_ds                            ! extinction coefficient of dry snow
      real(r8) :: ks                               ! scattering coefficient of snow
      real(r8) :: b_ds, b_ws
      real(r8) :: wk_h                             ! [m], equal to cmem cmem_snow_set_var.F90:476
      integer  :: i

!-----------------------------------------------------------------------
      IF (snowdp > 0.01) THEN  ! > 1cm
         ! calculate dielectric constant of ice
         CALL diel_ice(t_snow - tfrz, eps_i)
         eps_i_r = real(eps_i)
         eps_i_i = aimag(eps_i)

         ! consider the effect of salinity on the dielectric constant of ice
         eps_i_is = 0.0026/fghz + 0.00023*(fghz**0.87) ! impure ice -5(C)
         eps_i_ip = 6.e-4/fghz + 6.5e-5*(fghz**1.07)   ! pure ice -5(C)
         eps_i_i = eps_i_i + (eps_i_is - eps_i_ip)*sal_snow/13.0d0 ! corrected imaginary part of diel cons of ice

         ! calculate dielectric constant of dry snow (mixed by air and ice) (Polder–van Santen mixing model)
         rho_ds = (rho_snow - liq_snow)/(1.0 - liq_snow)   ! caculate density of dry snow,
         eps_ds_r = 1.0 + 1.58*rho_ds/(1.0 - 0.365*rho_ds)
         eps_ds_i = 3.0*(rho_ds/rho_i)*eps_i_i*(eps_ds_r**2)*(2*eps_ds_r + 1)/ &
            ((eps_i_r + 2*eps_ds_r)*(eps_i_r + 2*eps_ds_r**2))   ! negative imaginary part of diel cons of dry snow

         ! calculate dielectric constant of wet snow (Matzler 1987)
         IF (liq_snow > 0.) THEN  ! wet snow
            ! caculate relaxation frequency of three parts (eq.2.26)
            fa = f0w*(1 + (aa*(eps_w_s - eps_w_inf)/(eps_ds_r + (aa*(eps_w_inf - eps_ds_r)))))
            fb = f0w*(1 + (bb*(eps_w_s - eps_w_inf)/(eps_ds_r + (bb*(eps_w_inf - eps_ds_r)))))
            fc = f0w*(1 + (cc*(eps_w_s - eps_w_inf)/(eps_ds_r + (cc*(eps_w_inf - eps_ds_r)))))

            ! caculate infinite frequency dielectric constant of three parts
            eps_a_inf = (liq_snow*(eps_w_inf - eps_ds_r)/3.)/ &
               (1.+aa*((eps_w_inf/eps_ds_r) - 1.))
            eps_b_inf = (liq_snow*(eps_w_inf - eps_ds_r)/3.)/ &
               (1.+bb*((eps_w_inf/eps_ds_r) - 1.))
            eps_c_inf = (liq_snow*(eps_w_inf - eps_ds_r)/3.)/ &
               (1.+cc*((eps_w_inf/eps_ds_r) - 1.))

            ! caculate static dielectric constant of three parts
            eps_a_s = (liq_snow/3.)*(eps_w_s - eps_ds_r)/ &
               (1.+aa*((eps_w_s/eps_ds_r) - 1.))
            eps_b_s = (liq_snow/3.)*(eps_w_s - eps_ds_r)/ &
               (1.+bb*((eps_w_s/eps_ds_r) - 1.))
            eps_c_s = (liq_snow/3.)*(eps_w_s - eps_ds_r)/ &
               (1.+cc*((eps_w_s/eps_ds_r) - 1.))

            ! Debye equations
            eps_a = eps_a_inf + (eps_a_s - eps_a_inf)/(1 + jj*fghz/fa)
            eps_b = eps_b_inf + (eps_b_s - eps_b_inf)/(1 + jj*fghz/fb)
            eps_c = eps_c_inf + (eps_c_s - eps_c_inf)/(1 + jj*fghz/fc)

            ! calculate dielectric constant of wet snow
            eps = eps_a + eps_b + eps_c + (eps_ds_r - jj*eps_ds_i)
            eps_ws_r = real(eps)
            eps_ws_i = -1.*aimag(eps)
         ELSE
            eps_ws_r = eps_ds_r
            eps_ws_i = eps_ds_i
         END IF

         ! caculate propogation angle in snow (change medium from air to snow)
         alpha = k*abs(aimag(sqrt(eps_ws_r - jj*eps_ws_i)))
         beta = k*real(sqrt(eps_ws_r - jj*eps_ws_i))
         pp = 2.*alpha*beta
         qq = beta**2 - alpha**2 - (k*k)*(sin(theta)**2)
         theta_s = atan(k*sin(theta)/((1./sqrt(2.)) &
            *sqrt(sqrt(pp**2.+qq**2.) + qq)))

         ! caclulate wave impedance in snow
         z_s = z0/sqrt(eps_ws_r - jj*eps_ws_i)

         ! calculate brightness temperature above snow for H- V- polarization
         DO i = 1, 2
            ! Fresnel reflection coefficient between snow and air
            IF (i == 1) THEN
               r_sa = abs((z_s*cos(theta) - z0*cos(theta_s))/ &
                  (z_s*cos(theta) + z0*cos(theta_s)))**2
            ELSE
               r_sa = abs((z0*cos(theta) - z_s*cos(theta_s))/ &
                  (z0*cos(theta) + z_s*cos(theta_s)))**2
            END IF

            ! calculate snow grain size (mm) (Anderson 1976, eq.5.1) 
            d = min(1000*(1.6e-4 + 1.1e-13*((rho_snow*1000.0)**4)), 3.0)

            ! extinction coefficient of dry snow
            !//TODO: the paper is not focus on L-band, thus the formula is not suitable
            ke_ds = 0.0018*(fghz**2.8)*(d**2)/4.3429    !  [3](14)

            ! absorption coefficient of dry snow
            b_ds = (eps_ds_i/eps_ds_r)**2
            ka_ds = 2.*omega*sqrt(mu0*eps_0*eps_ds_r)* &
               sqrt(b_ds/(2.*(sqrt(1.+b_ds) + 1.)))

            IF (ke_ds < ka_ds) ke_ds = ka_ds

            ! absorption coefficient of wet snow
            b_ws = (eps_ws_i/eps_ws_r)**2
            ka_ws = 2.*omega*sqrt(mu0*eps_0*eps_ws_r)* &
               sqrt(b_ws/(2.*(sqrt(1.+b_ws) + 1)))

            ! total extinction (assuming scattering is the same for dry and wet snow)
            ke = (ke_ds - ka_ds) + ka_ws

            ! scattering coefficient of dry and wet snow
            ks = ke_ds - ka_ds

            wk_h = snowdp/rho_snow
            IF (((ke - q*ks)*(1./cos(theta_s))*wk_h) > (log(HUGE(1.)) - 1.)) THEN
               l2_apu = sqrt(HUGE(1.))
            ELSE
               l2_apu = min(sqrt(HUGE(1.)), exp((ke - q*ks)*(1./cos(theta_s))*wk_h))
            END IF

            ! brightness temperature through snow-air layer ([4] pp.243, eq.4.161)
            tb_2 = (1.+r_sn(i)/l2_apu)*(1.-r_sa)*t_snow* &
               (ka_ws/(ke - q*ks))*(1.-1./l2_apu)/ &
               (1.-r_sn(i)*r_sa/l2_apu**2.)

            ! brightness temperature through soil-snow layer ([4] pp.243, eq.4.162)
            tb_3 = ((1.-r_sn(i))*(1.-r_sa)*t)/ &
               (l2_apu*(1.-r_sn(i)*r_sa/l2_apu**2.))

            ! brightness temperature of snow-cover ground
            tb_tos(i) = tb_2 + tb_3

            ! emission by layers 2 (snow) and 3 (soil)
            ! reflectivity by layers 2 (snow) and 3 (soil)
            r_snow(i) = 1 - (tb_2/t_snow + tb_3/t)
         END DO
      END IF

   END SUBROUTINE snow
!-----------------------------------------------------------------------
END MODULE MOD_DA_Operator_RTM
#endif
