#include <define.h>
MODULE MOD_LeafInterception
! -----------------------------------------------------------------
! !DESCRIPTION:
! For calculating vegetation canopy precipitation interception.
!
! This MODULE is the coupler for the colm and CaMa-Flood model.
! CoLM2024 interception: Zhongwang Wei @ SYSU.
! Li, Q., Jin, X., Wei, Z. et al. (2026), Three-dimensional canopy morphology
! and wind dynamics govern global rainfall interception.
! Commun. Earth Environ. 7, 699. doi:10.1038/s43247-026-03694-7

!ANCILLARY FUNCTIONS AND SUBROUTINES
!-------------------
   !* :SUBROUTINE:"LEAF_interception_CoLM2014" : Leaf interception and drainage schemes based on colm2014 version
   !* :SUBROUTINE:"LEAF_interception_CoLM2024" : Canopy-morphology and wind-dependent interception
   !* :SUBROUTINE:"LEAF_interception_pftwrap"  : wrapper for pft land use classification

!REVISION HISTORY:
!----------------
   ! 2026.01     Zhongwang Wei: Fully revise CLM4,5,Noah-MP,MATSIRO,VIC and JULES schemes.
   ! 2024.04     Hua Yuan: add option to account for vegetation snow process based on Niu et al., 2004
   ! 2023.07     Hua Yuan: remove wrapper PC by using PFT leaf interception
   ! 2023.06     Shupeng Zhang @ SYSU
   ! 2023.02.23  Zhongwang Wei @ SYSU
   ! 2021.12.12  Zhongwang Wei @ SYSU
   ! 2020.10.21  Zhongwang Wei @ SYSU
   ! 2019.06     Hua Yuan: 1) add wrapper for PFT and PC, and 2) remove sigf by using lai+sai
   ! 2014.04     Yongjiu Dai
   ! 2002.08.31  Yongjiu Dai
   USE MOD_Precision
   USE MOD_Const_Physical, only: tfrz, denh2o, denice, cpliq, cpice, hfus
   USE MOD_Namelist, only: DEF_Interception_scheme, DEF_VEG_SNOW
#ifdef MaterialTransport
   USE MOD_MT_CanopyDiagnostics, only: mt_canopy_drainage, mt_canopy_drain
#endif

   IMPLICIT NONE
#ifdef MaterialTransport
   PRIVATE :: mt_canopy_drainage, mt_canopy_drain
#endif

   real(r8), parameter ::  CICE        = 2.094E06  !specific heat capacity of ice (j/m3/k)
   real(r8), parameter ::  bp          = 20.
   real(r8), parameter ::  CWAT        = 4.188E06  !specific heat capacity of water (j/m3/k)
   real(r8), parameter ::  pcoefs(2,2) = reshape((/20.0_r8, 0.206e-8_r8, 0.0001_r8, 0.9999_r8/), (/2,2/))

   ! Minimum significant precipitation rate threshold [mm/s]
   ! Used across all schemes for numerical stability
   real(r8), parameter ::  PRECIP_THRESHOLD = 1.0e-8_r8

   ! Tolerance for interception water balance checks [mm]
   ! Used by check_interception_balance subroutine under CoLMDEBUG
   real(r8), parameter ::  INTERCEPTION_BALANCE_TOL = 1.0e-5_r8


   !----------------------- Dummy argument --------------------------------
   real(r8) :: satcap                     ! maximum allowed water on canopy [mm]
   real(r8) :: satcap_rain                ! maximum allowed rain on canopy [mm]
   real(r8) :: satcap_snow                ! maximum allowed snow on canopy [mm]
   real(r8) :: lsai                       ! sum of leaf area index and stem area index [-]
   real(r8) :: chiv                       ! leaf angle distribution factor
   real(r8) :: ppc                        ! convective precipitation in time-step [mm]
   real(r8) :: ppl                        ! large-scale precipitation in time-step [mm]
   real(r8) :: p0                         ! precipitation in time-step [mm]
   real(r8) :: fpi                        ! coefficient of interception
   real(r8) :: fpi_rain                   ! coefficient of interception of rain
   real(r8) :: fpi_snow                   ! coefficient of interception of snow
   real(r8) :: alpha_rain                 ! coefficient of interception of rain
   real(r8) :: alpha_snow                 ! coefficient of interception of snow
   real(r8) :: pinf                       ! interception of precipitation in time step [mm]
   real(r8) :: tti_rain                   ! direct rain throughfall in time step [mm]
   real(r8) :: tti_snow                   ! direct snow throughfall in time step [mm]
   real(r8) :: tex_rain                   ! canopy rain drainage in time step [mm]
   real(r8) :: tex_snow                   ! canopy snow drainage in time step [mm]
   real(r8) :: vegt                       ! sigf*lsai
   real(r8) :: xs                         ! proportion of the grid area where the intercepted rainfall
                                          ! plus the preexisting canopy water storage
   real(r8)  :: unl_snow_temp,U10,unl_snow_wind,unl_snow
   real(r8)  :: ap, cp, aa1, bb1, exrain, arg, w
   real(r8)  :: thru_rain, thru_snow
   real(r8)  :: xsc_rain, xsc_snow

   real(r8)  :: fvegc                     ! vegetation fraction
   real(r8)  :: FT                        ! the temperature factor for snow unloading
   real(r8)  :: FV                        ! the wind factor for snow unloading
   real(r8)  :: ICEDRIP                   ! snow unloading

   real(r8)  :: ldew_smelt
   real(r8)  :: ldew_frzc
   real(r8)  :: FP
   real(r8)  :: int_rain
   real(r8)  :: int_snow

CONTAINS

   PURE REAL(r8) FUNCTION canopy_storage_capacity_colm2024 (dewmx,lai,sai,forc_us,forc_vs, &
                                                             htop,ncd,ncw,bcw,veg_class,is_pft)
      USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
      real(r8), intent(in) :: dewmx, lai, sai, forc_us, forc_vs
      real(r8), intent(in) :: htop, ncd, ncw, bcw
      integer,  intent(in) :: veg_class
      logical,  intent(in) :: is_pft

      integer :: canopy_type
      real(r8) :: wind, needle_cap, broad_cap, crown_ratio
      logical :: needle_valid, broad_valid

      canopy_storage_capacity_colm2024 = dewmx * max(0._r8, lai+sai)
      canopy_type = 0

      IF (is_pft) THEN
         IF (veg_class >= 1 .and. veg_class <= 3) canopy_type = 1
         IF (veg_class >= 4 .and. veg_class <= 8) canopy_type = 2
         IF (veg_class >= 9 .and. veg_class <= 11) canopy_type = 3
      ELSE
#ifdef LULC_USGS
         SELECT CASE (veg_class)
         CASE (12,14)
            canopy_type = 1
         CASE (11,13)
            canopy_type = 2
         CASE (8)
            canopy_type = 3
         CASE (15)
            canopy_type = 4
         END SELECT
#elif defined LULC_IGBP
         SELECT CASE (veg_class)
         CASE (1,3)
            canopy_type = 1
         CASE (2,4)
            canopy_type = 2
         CASE (5)
            canopy_type = 4
         CASE (6,7)
            canopy_type = 3
         END SELECT
#endif
      ENDIF

      IF (canopy_type == 0) RETURN
      wind = sqrt(max(0._r8, forc_us*forc_us + forc_vs*forc_vs))
      needle_valid = .false.
      broad_valid = .false.
      IF (canopy_type == 1 .or. canopy_type == 4) THEN
         IF (all(ieee_is_finite([ncd, ncw]))) &
            needle_valid = ncd > 0._r8 .and. ncd < 1000._r8 .and. ncw > 0._r8 .and. ncw < 1000._r8
      ENDIF
      IF (canopy_type == 2 .or. canopy_type == 4) THEN
         IF (all(ieee_is_finite([bcw, htop]))) &
            broad_valid = bcw > 0._r8 .and. bcw < 1000._r8 .and. htop > 0._r8 .and. htop < 1000._r8
      ENDIF

      IF (needle_valid) needle_cap = &
         (min(11._r8,max(3._r8,ncd)) + min(7._r8,max(2.9_r8,ncw))) / &
         (4._r8 * (1._r8 + min(3.6_r8,max(1._r8,wind))))

      IF (broad_valid) THEN
         crown_ratio = min(7._r8,max(1._r8,htop/bcw))
         broad_cap = min(8._r8,max(2._r8,bcw)) / &
            (2._r8 * (min(4._r8,max(1.5_r8,wind)) + crown_ratio))
      ENDIF

      SELECT CASE (canopy_type)
      CASE (1)
         IF (needle_valid) canopy_storage_capacity_colm2024 = needle_cap
      CASE (2)
         IF (broad_valid) canopy_storage_capacity_colm2024 = broad_cap
      CASE (3)
         canopy_storage_capacity_colm2024 = &
            0.5_r8 * (1._r8 + 1._r8/(1._r8 + min(4._r8,max(1._r8,wind))))
      CASE (4)
         IF (needle_valid .and. broad_valid) &
            canopy_storage_capacity_colm2024 = 0.5_r8 * (needle_cap + broad_cap)
      END SELECT
   END FUNCTION canopy_storage_capacity_colm2024

   SUBROUTINE LEAF_interception_CoLM2014 (deltim,dewmx,forc_us,forc_vs,chil,sigf,lai,sai,tair,tleaf,&
                                          prc_rain,prc_snow,prl_rain,prl_snow,bifall,qflx_irrig_sprinkler,&
                                          ldew,ldew_rain,ldew_snow,z0m,hu,pg_rain,pg_snow,qintr,qintr_rain,qintr_snow,satcap_rain_override)
!DESCRIPTION
!===========
   ! Calculation of  interception and drainage of precipitation
   ! the treatment are based on Sellers et al. (1996)

!Original Author:
!-------------------
   !canopy interception scheme modified by Yongjiu Dai based on Sellers et al. (1996)

!References:
!-------------------
   !---Dai, Y., Zeng, X., Dickinson, R.E., Baker, I., Bonan, G.B., BosiloVICh,
   !   M.G., Denning, A.S., Dirmeyer, P.A., Houser, P.R., Niu, G. and Oleson,
   !   K.W., 2003.  The common land model. Bulletin of the American
   !   Meteorological Society, 84(8), pp.1013-1024.

   !---Lawrence, D.M., Thornton, P.E., Oleson, K.W. and Bonan, G.B., 2007.  The
   !   partitioning of evapotranspiration into transpiration, soil evaporation,
   !   and canopy evaporation in a GCM: Impacts on land-atmosphere interaction.
   !   Journal of Hydrometeorology, 8(4), pp.862-880.

   !---Oleson, K., Dai, Y., Bonan, B., BosiloVIChm, M., Dickinson, R.,
   !   Dirmeyer, P., Hoffman, F., Houser, P., Levis, S., Niu, G.Y. and
   !   Thornton, P., 2004.  Technical description of the community land model
   !   (CLM).

   !---Sellers, P.J., Randall, D.A., Collatz, G.J., Berry, J.A., Field, C.B.,
   !   Dazlich, D.A., Zhang, C., Collelo, G.D. and Bounoua, L., 1996. A revised
   !   land surface parameterization (SiB2) for atmospheric GCMs.  Part I:
   !   Model formulation. Journal of climate, 9(4), pp.676-705.

   !---Sellers, P.J., Tucker, C.J., Collatz, G.J., Los, S.O., Justice, C.O.,
   !   Dazlich, D.A. and Randall, D.A., 1996.  A revised land surface
   !   parameterization (SiB2) for atmospheric GCMs. Part II: The generation of
   !   global fields of terrestrial biophysical parameters from satellite data.
   !   Journal of climate, 9(4), pp.706-737.


!ANCILLARY FUNCTIONS AND SUBROUTINES
!-------------------

!REVISION HISTORY
!----------------
   !---2024.04.16  Hua Yuan: add option to account for vegetation snow process based on Niu et al., 2004
   !---2023.02.21  Zhongwang Wei @ SYSU : Snow and rain interception
   !---2021.12.08  Zhongwang Wei @ SYSU
   !---2019.06     Hua Yuan: remove sigf and USE lai+sai for judgement.
   !---2014.04     Yongjiu Dai
   !---2002.08.31  Yongjiu Dai
!=======================================================================

   IMPLICIT NONE

   real(r8), intent(in) :: deltim       !seconds in a time step [second]
   real(r8), intent(in) :: dewmx        !maximum dew [mm]
   real(r8), intent(in) :: forc_us      !wind speed
   real(r8), intent(in) :: forc_vs      !wind speed
   real(r8), intent(in) :: chil         !leaf angle distribution factor
   real(r8), intent(in) :: prc_rain     !convective rainfall [mm/s]
   real(r8), intent(in) :: prc_snow     !convective snowfall [mm/s]
   real(r8), intent(in) :: prl_rain     !large-scale rainfall [mm/s]
   real(r8), intent(in) :: prl_snow     !large-scale snowfall [mm/s]
   real(r8), intent(in) :: qflx_irrig_sprinkler ! irrigation and sprinkler water flux [mm/s]
   real(r8), intent(in) :: bifall       !bulk density of newly fallen dry snow [kg/m3]
   real(r8), intent(in) :: sigf         !fraction of veg cover, excluding snow-covered veg [-]
   real(r8), intent(in) :: lai          !leaf area index [-]
   real(r8), intent(in) :: sai          !stem area index [-]
   real(r8), intent(in) :: tair         !air temperature [K]
   real(r8), intent(in) :: tleaf        !sunlit canopy leaf temperature [K]

   real(r8), intent(inout) :: ldew      !depth of water on foliage [mm]
   real(r8), intent(inout) :: ldew_rain !depth of water on foliage [mm]
   real(r8), intent(inout) :: ldew_snow !depth of water on foliage [mm]
   real(r8), intent(in)    :: z0m       !roughness length
   real(r8), intent(in)    :: hu        !forcing height of U

   real(r8), intent(out) :: pg_rain     !rainfall onto ground including canopy runoff [kg/(m2 s)]
   real(r8), intent(out) :: pg_snow     !snowfall onto ground including canopy runoff [kg/(m2 s)]
   real(r8), intent(out) :: qintr       !interception [kg/(m2 s)]
   real(r8), intent(out) :: qintr_rain  !rainfall interception (mm h2o/s)
   real(r8), intent(out) :: qintr_snow  !snowfall interception (mm h2o/s)
   real(r8), intent(in), optional :: satcap_rain_override

!-----------------------------------------------------------------------

      IF (lai+sai > 1e-6) THEN
         lsai   = lai + sai
         vegt   = lsai
         satcap = dewmx*vegt
         satcap_rain = satcap
         IF (present(satcap_rain_override)) THEN
            satcap_rain = max(0._r8, satcap_rain_override)
            IF (.not. DEF_VEG_SNOW) satcap = satcap_rain
         ENDIF
         satcap_snow = 6.6*(0.27+46./bifall)*vegt  ! Niu et al., 2004
         satcap_snow = 48.*satcap                  ! Simple one without snow density input

         p0  = (prc_rain + prc_snow + prl_rain + prl_snow + qflx_irrig_sprinkler)*deltim
         ppc = (prc_rain + prc_snow)*deltim
         ppl = (prl_rain + prl_snow + qflx_irrig_sprinkler)*deltim

         w = ldew+p0
         IF (tleaf > tfrz) THEN
            xsc_rain = max(0., ldew-satcap)
            xsc_snow = 0.
         ELSE
            xsc_rain = 0.
            xsc_snow = max(0., ldew-satcap)
         ENDIF

         ldew = ldew - (xsc_rain + xsc_snow)

         !TODO-done: account for vegetation snow
         IF ( DEF_VEG_SNOW ) THEN
            xsc_rain  = max(0., ldew_rain-satcap_rain)
            xsc_snow  = max(0., ldew_snow-satcap_snow)
            ldew_rain = ldew_rain - xsc_rain
            ldew_snow = ldew_snow - xsc_snow
            ldew      = ldew_rain + ldew_snow
         ENDIF

         ap = pcoefs(2,1)
         cp = pcoefs(2,2)

         IF (p0 > 1.e-8) THEN
            ap = ppc/p0 * pcoefs(1,1) + ppl/p0 * pcoefs(2,1)
            cp = ppc/p0 * pcoefs(1,2) + ppl/p0 * pcoefs(2,2)

            !----------------------------------------------------------------------
            !      proportional saturated area (xs) and leaf drainage(tex)
            !-----------------------------------------------------------------------
            chiv = chil
            IF ( abs(chiv) .le. 0.01 ) chiv = 0.01
            aa1 = 0.5 - 0.633 * chiv - 0.33 * chiv * chiv
            bb1 = 0.877 * ( 1. - 2. * aa1 )
            exrain = aa1 + bb1

            ! coefficient of interception
            ! set fraction of potential interception to max 0.25 (Lawrence et al. 2007)
            ! assume alpha_rain = alpha_snow
            alpha_rain = 0.25
            fpi = alpha_rain * ( 1.-exp(-exrain*lsai) )
            tti_rain = (prc_rain+prl_rain+qflx_irrig_sprinkler)*deltim * ( 1.-fpi )
            tti_snow = (prc_snow+prl_snow)*deltim * ( 1.-fpi )

            xs = 1.
            IF (p0*fpi>1.e-9) THEN
               arg = (satcap-ldew)/(p0*fpi*ap) - cp/ap
               IF (arg>1.e-9) THEN
                  xs = -1./bp * log( arg )
                  xs = min( xs, 1. )
                  xs = max( xs, 0. )
               ENDIF
            ENDIF

            ! assume no fall down of the intercepted snowfall in a time step
            ! drainage
            tex_rain = (prc_rain+prl_rain+qflx_irrig_sprinkler)*deltim * fpi * (ap/bp*(1.-exp(-bp*xs))+cp*xs) &
                     - max(0., (satcap-ldew)) * xs
            tex_rain = max( tex_rain, 0. )
            ! Ensure physical constraint: tex_rain + tti_rain <= total rain input
            tex_rain = min( tex_rain, (prc_rain+prl_rain+qflx_irrig_sprinkler)*deltim - tti_rain )
            tex_snow = 0.

            ! 04/11/2024, yuan:
            !TODO-done: account for snow on vegetation,
            IF ( DEF_VEG_SNOW ) THEN

               ! re-calculate leaf rain drainage using ldew_rain

               xs = 1.
               IF (p0*fpi>1.e-9) THEN
                  arg = (satcap_rain-ldew_rain)/(p0*fpi*ap) - cp/ap
                  IF (arg>1.e-9) THEN
                     xs = -1./bp * log( arg )
                     xs = min( xs, 1. )
                     xs = max( xs, 0. )
                  ENDIF
               ENDIF

               tex_rain = (prc_rain+prl_rain+qflx_irrig_sprinkler)*deltim * fpi * (ap/bp*(1.-exp(-bp*xs))+cp*xs) &
                        - max(0., (satcap_rain-ldew_rain)) * xs
               tex_rain = max( tex_rain, 0. )
               ! Ensure physical constraint: tex_rain + tti_rain <= total rain input
               tex_rain = min( tex_rain, (prc_rain+prl_rain+qflx_irrig_sprinkler)*deltim - tti_rain )

               ! re-calculate the snow loading rate

               fvegc = 1. - exp(-0.52*lsai)
               FP    = (ppc + ppl) / (10.*ppc + ppl)
               qintr_snow = fvegc * (prc_snow+prl_snow) * FP
               qintr_snow = min (qintr_snow, (satcap_snow-ldew_snow)/deltim * (1.-exp(-(prc_snow+prl_snow)*deltim/satcap_snow)) )
               qintr_snow = max (qintr_snow, 0.)

               ! snow unloading rate

               FT = max(0.0, (tleaf - tfrz) / 1.87e5)
               FV = sqrt(forc_us*forc_us + forc_vs*forc_vs) / 1.56e5
               tex_snow = max(0., ldew_snow/deltim) * (FV+FT)
               tti_snow = (1.0-fvegc)*(prc_snow+prl_snow) + (fvegc*(prc_snow+prl_snow) - qintr_snow)

               ! rate -> mass

               tti_snow = tti_snow * deltim
               tex_snow = tex_snow * deltim
            ENDIF

#if (defined CoLMDEBUG)
            IF (tex_rain+tex_snow+tti_rain+tti_snow-p0 > 1.e-10 .and. .not.DEF_VEG_SNOW) THEN
               write(6,*) 'tex_ + tti_ > p0 in interception code : ',ldew,tex_rain,tex_snow,tti_rain,tti_snow,p0
            ENDIF
#endif

         ELSE
            ! all intercepted by canopy leaves for very small precipitation
            tti_rain = 0.
            tti_snow = 0.
            tex_rain = 0.
            tex_snow = 0.
         ENDIF

         !----------------------------------------------------------------------
         !   total throughfall (thru) and store augmentation
         !----------------------------------------------------------------------

         thru_rain = tti_rain + tex_rain
         thru_snow = tti_snow + tex_snow
         pinf = p0 - (thru_rain + thru_snow)
         ldew = ldew + pinf

         !TODO-done: IF DEF_VEG_SNOW, update ldew_rain, ldew_snow
         IF ( DEF_VEG_SNOW ) THEN
            ldew_rain = ldew_rain + (prc_rain+prl_rain+qflx_irrig_sprinkler)*deltim - thru_rain
            ldew_snow = ldew_snow + (prc_snow+prl_snow)*deltim - thru_snow
            ldew = ldew_rain + ldew_snow
         ENDIF

         pg_rain = (xsc_rain + thru_rain) / deltim
         pg_snow = (xsc_snow + thru_snow) / deltim
         qintr   = pinf / deltim

         qintr_rain = prc_rain + prl_rain + qflx_irrig_sprinkler - thru_rain / deltim
         qintr_snow = prc_snow + prl_snow - thru_snow / deltim

#ifdef MaterialTransport
         CALL mt_canopy_drain%capture(deltim, tex_rain, tex_snow)
#endif

#if (defined CoLMDEBUG)
         w = w - ldew - (pg_rain+pg_snow)*deltim
         IF (abs(w) > INTERCEPTION_BALANCE_TOL) THEN
            write(6,*) 'something wrong in interception code: '
            write(6,*) w, ldew, (pg_rain+pg_snow)*deltim, satcap
            CALL abort
         ENDIF

         IF (DEF_VEG_SNOW .and. abs(ldew-ldew_rain-ldew_snow) > INTERCEPTION_BALANCE_TOL) THEN
            write(6,*) 'something wrong in interception code when DEF_VEG_SNOW: '
            write(6,*) ldew, ldew_rain, ldew_snow
            CALL abort
         ENDIF
#endif

      ELSE
         ! 07/15/2023, yuan: #bug found for ldew value reset.
         !NOTE: this bug should exist in other interception schemes @Zhongwang.
         IF (ldew > 0.) THEN
            IF (tleaf > tfrz) THEN
               pg_rain = prc_rain + prl_rain + qflx_irrig_sprinkler + ldew/deltim
               pg_snow = prc_snow + prl_snow
            ELSE
               pg_rain = prc_rain + prl_rain + qflx_irrig_sprinkler
               pg_snow = prc_snow + prl_snow + ldew/deltim
            ENDIF
         ELSE
            pg_rain = prc_rain + prl_rain + qflx_irrig_sprinkler
            pg_snow = prc_snow + prl_snow
         ENDIF

         ldew       = 0.
         ldew_rain  = 0.
         ldew_snow  = 0.
         qintr      = 0.
         qintr_rain = 0.
         qintr_snow = 0.

#ifdef MaterialTransport
         CALL mt_canopy_drain%clear()
#endif

      ENDIF

   END SUBROUTINE LEAF_interception_CoLM2014

   SUBROUTINE LEAF_interception_CoLM2024 (deltim,dewmx,forc_us,forc_vs,chil,sigf,lai,sai,tair,tleaf,&
                                          prc_rain,prc_snow,prl_rain,prl_snow,qflx_irrig_sprinkler,&
                                          bifall,veg_class,is_pft,ncd,ncw,bcw,htop,&
                                          ldew,ldew_rain,ldew_snow,z0m,hu,pg_rain,pg_snow,&
                                          qintr,qintr_rain,qintr_snow)
!DESCRIPTION
!===========
   !CoLM2014 interception with the rain storage capacity set by canopy morphology and wind

!Original Author:
!-------------------
   !---Zhongwang Wei @ SYSU

!References:
!-------------------
   !---Li, Q., Jin, X., Wei, Z. et al., 2026. Three-dimensional canopy morphology and wind
   !   dynamics govern global rainfall interception. Commun. Earth Environ., 7, 699.
   !   doi:10.1038/s43247-026-03694-7

   IMPLICIT NONE
   real(r8), intent(in) :: deltim, dewmx, forc_us, forc_vs, chil, sigf, lai, sai, tair
   real(r8), intent(inout) :: tleaf
   real(r8), intent(in) :: prc_rain, prc_snow, prl_rain, prl_snow, qflx_irrig_sprinkler, bifall
   integer, intent(in) :: veg_class
   logical, intent(in) :: is_pft
   real(r8), intent(in) :: ncd, ncw, bcw, htop
   real(r8), intent(inout) :: ldew, ldew_rain, ldew_snow
   real(r8), intent(in) :: z0m, hu
   real(r8), intent(out) :: pg_rain, pg_snow, qintr, qintr_rain, qintr_snow
   real(r8) :: satcap_2024
      satcap_2024 = dewmx * max(0._r8, lai+sai)
      IF (lai+sai > 1.e-6_r8) THEN
         satcap_2024 = canopy_storage_capacity_colm2024 (dewmx,lai,sai,forc_us,forc_vs,&
                                                         htop,ncd,ncw,bcw,veg_class,is_pft)
      ENDIF
      CALL LEAF_interception_CoLM2014 (deltim,dewmx,forc_us,forc_vs,chil,sigf,lai,sai,tair,tleaf,&
         prc_rain,prc_snow,prl_rain,prl_snow,bifall,qflx_irrig_sprinkler,&
         ldew,ldew_rain,ldew_snow,z0m,hu,pg_rain,pg_snow,qintr,qintr_rain,qintr_snow,&
         satcap_rain_override=satcap_2024)
   END SUBROUTINE LEAF_interception_CoLM2024

   SUBROUTINE LEAF_interception_wrap(deltim,dewmx,forc_us,forc_vs,chil,sigf,lai,sai,tair,tleaf, &
                               prc_rain,prc_snow,prl_rain,prl_snow,bifall, &
                               patchclass,ncd,ncw,bcw,htop, &
                                                       ldew,ldew_rain,ldew_snow,z0m,hu,pg_rain, &
                                                            pg_snow,qintr,qintr_rain,qintr_snow,&
                                                                           qflx_irrig_sprinkler )
!DESCRIPTION
!===========
   !wrapper for calculation of canopy interception using USGS or IGBP land cover classification

!ANCILLARY FUNCTIONS AND SUBROUTINES
!-------------------

!Original Author:
!-------------------
   !---Shupeng Zhang

!References:


!REVISION HISTORY
!----------------

   IMPLICIT NONE

   real(r8), intent(in)    :: deltim     !seconds in a time step [second]
   real(r8), intent(in)    :: dewmx      !maximum dew [mm]
   real(r8), intent(in)    :: forc_us    !wind speed
   real(r8), intent(in)    :: forc_vs    !wind speed
   real(r8), intent(in)    :: chil       !leaf angle distribution factor
   real(r8), intent(in)    :: prc_rain   !convective rainfall [mm/s]
   real(r8), intent(in)    :: prc_snow   !convective snowfall [mm/s]
   real(r8), intent(in)    :: prl_rain   !large-scale rainfall [mm/s]
   real(r8), intent(in)    :: prl_snow   !large-scale snowfall [mm/s]
   real(r8), intent(in)    :: bifall     !bulk density of newly fallen dry snow [kg/m3]
   integer, intent(in)     :: patchclass
   real(r8), intent(in)    :: ncd, ncw, bcw, htop
   real(r8), intent(in)    :: sigf       !fraction of veg cover, excluding snow-covered veg [-]
   real(r8), intent(in)    :: lai        !leaf area index [-]
   real(r8), intent(in)    :: sai        !stem area index [-]
   real(r8), intent(in)    :: tair       !air temperature [K]
   real(r8), intent(inout) :: tleaf      !sunlit canopy leaf temperature [K]

   real(r8), intent(inout) :: ldew       !depth of water on foliage [mm]
   real(r8), intent(inout) :: ldew_rain  !depth of liquid on foliage [mm]
   real(r8), intent(inout) :: ldew_snow  !depth of liquid on foliage [mm]
   real(r8), intent(in)    :: z0m        !roughness length
   real(r8), intent(in)    :: hu         !forcing height of U


   real(r8), intent(out)   :: pg_rain    !rainfall onto ground including canopy runoff [kg/(m2 s)]
   real(r8), intent(out)   :: pg_snow    !snowfall onto ground including canopy runoff [kg/(m2 s)]
   real(r8), intent(out)   :: qintr      !interception [kg/(m2 s)]
   real(r8), intent(out)   :: qintr_rain !rainfall interception (mm h2o/s)
   real(r8), intent(out)   :: qintr_snow !snowfall interception (mm h2o/s)

   real(r8), intent(in), optional :: qflx_irrig_sprinkler

   real(r8) :: qflx_irrig_sprinkler_

      IF ( present(qflx_irrig_sprinkler) ) THEN
         qflx_irrig_sprinkler_ = qflx_irrig_sprinkler
      ELSE
         qflx_irrig_sprinkler_ = 0
      ENDIF

      IF (DEF_Interception_scheme==1) THEN
         CALL LEAF_interception_CoLM2014 (deltim,dewmx,forc_us,forc_vs,chil,sigf,lai,sai,tair,tleaf,&
                                             prc_rain,prc_snow,prl_rain,prl_snow,bifall,qflx_irrig_sprinkler_,&
                                             ldew,ldew_rain,ldew_snow,z0m,hu,pg_rain,&
                                             pg_snow,qintr,qintr_rain,qintr_snow)
      ELSEIF  (DEF_Interception_scheme==8) THEN
         CALL LEAF_interception_CoLM2024 (deltim,dewmx,forc_us,forc_vs,chil,sigf,lai,sai,tair,tleaf,&
                                             prc_rain,prc_snow,prl_rain,prl_snow,qflx_irrig_sprinkler_,&
                                             bifall,patchclass,.false.,ncd,ncw,bcw,htop,&
                                             ldew,ldew_rain,ldew_snow,z0m,hu,pg_rain,&
                                             pg_snow,qintr,qintr_rain,qintr_snow)
      ELSE
         write(6,*) 'LEAF_interception_wrap requires scheme 1 or 8'
         CALL abort
      ENDIF

   END SUBROUTINE LEAF_interception_wrap

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   SUBROUTINE LEAF_interception_pftwrap (ipatch,deltim,dewmx,forc_us,forc_vs,forc_t,&
                               prc_rain,prc_snow,prl_rain,prl_snow,bifall,&
                               ldew,ldew_rain,ldew_snow,z0m,hu,pg_rain,pg_snow,qintr,qintr_rain,qintr_snow,&
                               qflx_irrig_sprinkler)

! -----------------------------------------------------------------
! !DESCRIPTION:
! wrapper for calculation of canopy interception for PFTs within a land cover type.
!
! Created by Hua Yuan, 06/2019
!
! !REVISION HISTORY:
! 2023.02.21 Zhongwang Wei @ SYSU: add different options of canopy interception for PFTs
!
! -----------------------------------------------------------------

   USE MOD_Precision
   USE MOD_LandPFT
   USE MOD_Const_Physical, only: tfrz
   USE MOD_Vars_PFTimeInvariants
   USE MOD_Vars_PFTimeVariables
   USE MOD_Vars_1DPFTFluxes
   USE MOD_Const_PFT
   IMPLICIT NONE

   integer,  intent(in)    :: ipatch     !patch index
   real(r8), intent(in)    :: deltim     !seconds in a time step [second]
   real(r8), intent(in)    :: dewmx      !maximum dew [mm]
   real(r8), intent(in)    :: forc_us    !wind speed
   real(r8), intent(in)    :: forc_vs    !wind speed
   real(r8), intent(in)    :: forc_t     !air temperature
   real(r8), intent(in)    :: z0m        !roughness length
   real(r8), intent(in)    :: hu         !forcing height of U
   real(r8), intent(inout) :: ldew_rain  !depth of water on foliage [mm]
   real(r8), intent(inout) :: ldew_snow  !depth of water on foliage [mm]
   real(r8), intent(in)    :: prc_rain   !convective ranfall [mm/s]
   real(r8), intent(in)    :: prc_snow   !convective snowfall [mm/s]
   real(r8), intent(in)    :: prl_rain   !large-scale rainfall [mm/s]
   real(r8), intent(in)    :: prl_snow   !large-scale snowfall [mm/s]
   real(r8), intent(in)    :: bifall     ! bulk density of newly fallen dry snow [kg/m3]

   real(r8), intent(inout) :: ldew       !depth of water on foliage [mm]
   real(r8), intent(out)   :: pg_rain    !rainfall onto ground including canopy runoff [kg/(m2 s)]
   real(r8), intent(out)   :: pg_snow    !snowfall onto ground including canopy runoff [kg/(m2 s)]
   real(r8), intent(out)   :: qintr      !interception [kg/(m2 s)]
   real(r8), intent(out)   :: qintr_rain !rainfall interception (mm h2o/s)
   real(r8), intent(out)   :: qintr_snow !snowfall interception (mm h2o/s)

   real(r8), intent(in), optional :: qflx_irrig_sprinkler !irrigation and sprinkler water [mm/s]

   integer i, p, ps, pe
#ifdef CROP
   integer  :: irrig_flag  ! 1 if sprinker, 2 if others
#endif
#ifdef MaterialTransport
   type(mt_canopy_drainage) :: mt_drain_sum
#endif
   real(r8) pg_rain_tmp, pg_snow_tmp, qflx_irrig_sprinkler_

      pg_rain_tmp = 0.
      pg_snow_tmp = 0.

      IF ( present(qflx_irrig_sprinkler) ) THEN
         qflx_irrig_sprinkler_ = qflx_irrig_sprinkler
      ELSE
         qflx_irrig_sprinkler_ = 0
      ENDIF

      ps = patch_pft_s(ipatch)
      pe = patch_pft_e(ipatch)
#ifdef MaterialTransport
      CALL mt_drain_sum%clear()
#endif

      IF (DEF_Interception_scheme==1) THEN
         DO i = ps, pe
            p = pftclass(i)
            CALL LEAF_interception_CoLM2014 (deltim,dewmx,forc_us,forc_vs,chil_p(p),sigf_p(i),lai_p(i),sai_p(i),forc_t,tleaf_p(i),&
                                                prc_rain,prc_snow,prl_rain,prl_snow,bifall,qflx_irrig_sprinkler_,&
                                                ldew_p(i),ldew_rain_p(i),ldew_snow_p(i),z0m_p(i),hu,pg_rain,pg_snow,qintr_p(i),qintr_rain_p(i),qintr_snow_p(i))
            pg_rain_tmp = pg_rain_tmp + pg_rain*pftfrac(i)
            pg_snow_tmp = pg_snow_tmp + pg_snow*pftfrac(i)
#ifdef MaterialTransport
            CALL mt_drain_sum%accumulate(mt_canopy_drain, pftfrac(i))
#endif
         ENDDO
      ELSEIF (DEF_Interception_scheme==8) THEN
         DO i = ps, pe
            p = pftclass(i)
            CALL LEAF_interception_CoLM2024 (deltim,dewmx,forc_us,forc_vs,chil_p(p),sigf_p(i),lai_p(i),sai_p(i),forc_t,tleaf_p(i),&
                                             prc_rain,prc_snow,prl_rain,prl_snow,qflx_irrig_sprinkler_,&
                                             bifall,p,.true.,ncd_p(i),ncw_p(i),bcw_p(i),htop_p(i),&
                                             ldew_p(i),ldew_rain_p(i),ldew_snow_p(i),z0m_p(i),hu,pg_rain,pg_snow,qintr_p(i),qintr_rain_p(i),qintr_snow_p(i))
            pg_rain_tmp = pg_rain_tmp + pg_rain*pftfrac(i)
            pg_snow_tmp = pg_snow_tmp + pg_snow*pftfrac(i)
#ifdef MaterialTransport
            CALL mt_drain_sum%accumulate(mt_canopy_drain, pftfrac(i))
#endif
         ENDDO
      ELSE
         write(6,*) 'LEAF_interception_pftwrap requires scheme 1 or 8'
         CALL abort
      ENDIF

      pg_rain = pg_rain_tmp
      pg_snow = pg_snow_tmp
      ldew    = sum( ldew_p(ps:pe) * pftfrac(ps:pe))
      ldew_rain = sum( ldew_rain_p(ps:pe) * pftfrac(ps:pe))
      ldew_snow = sum( ldew_snow_p(ps:pe) * pftfrac(ps:pe))
      qintr   = sum(qintr_p(ps:pe) * pftfrac(ps:pe))
      qintr_rain = sum(qintr_rain_p(ps:pe) * pftfrac(ps:pe))
      qintr_snow = sum(qintr_snow_p(ps:pe) * pftfrac(ps:pe))
#ifdef MaterialTransport
      mt_canopy_drain = mt_drain_sum
#endif

   END SUBROUTINE LEAF_interception_pftwrap
#endif

   SUBROUTINE check_interception_balance(scheme_name, &
         ldew, ldew_rain, ldew_snow, pg_rain, pg_snow, &
         qintr, qintr_rain, qintr_snow)

   ! Validates interception water balance consistency.
   ! Called from CoLMDEBUG blocks after each scheme completes.

      character(len=*), intent(in) :: scheme_name
      real(r8), intent(in) :: ldew, ldew_rain, ldew_snow
      real(r8), intent(in) :: pg_rain, pg_snow
      real(r8), intent(in) :: qintr, qintr_rain, qintr_snow

      ! Check A: component consistency (ldew == ldew_rain + ldew_snow)
      IF (abs(ldew - (ldew_rain + ldew_snow)) > INTERCEPTION_BALANCE_TOL) THEN
         write(6,*) 'Component consistency error in ', scheme_name, ':'
         write(6,*) 'ldew=', ldew, ' ldew_rain+ldew_snow=', ldew_rain+ldew_snow
         write(6,*) 'diff=', ldew - (ldew_rain + ldew_snow)
         CALL abort
      ENDIF

      ! Check B: non-negativity
      IF (ldew < -INTERCEPTION_BALANCE_TOL .or. &
          ldew_rain < -INTERCEPTION_BALANCE_TOL .or. &
          ldew_snow < -INTERCEPTION_BALANCE_TOL .or. &
          pg_rain < -INTERCEPTION_BALANCE_TOL .or. &
          pg_snow < -INTERCEPTION_BALANCE_TOL) THEN
         write(6,*) 'Negative value error in ', scheme_name, ':'
         write(6,*) 'ldew=', ldew, ' ldew_rain=', ldew_rain, ' ldew_snow=', ldew_snow
         write(6,*) 'pg_rain=', pg_rain, ' pg_snow=', pg_snow
         CALL abort
      ENDIF

      ! Check C: flux consistency (qintr == qintr_rain + qintr_snow)
      IF (abs(qintr - (qintr_rain + qintr_snow)) > INTERCEPTION_BALANCE_TOL) THEN
         write(6,*) 'Flux consistency error in ', scheme_name, ':'
         write(6,*) 'qintr=', qintr, ' qintr_rain+qintr_snow=', qintr_rain+qintr_snow
         write(6,*) 'diff=', qintr - (qintr_rain + qintr_snow)
         CALL abort
      ENDIF

   END SUBROUTINE check_interception_balance

END MODULE MOD_LeafInterception
