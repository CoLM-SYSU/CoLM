#include <define.h>

MODULE MOD_ForcingDownscaling

!-----------------------------------------------------------------------------
! !DESCRIPTION:
!  Downscaling meteorological forcings
!
! !INITIAL:
!  The Community Land Model version 5.0 (CLM5.0)
!
! !REVISIONS:
!  Updated by Yongjiu Dai, January 16, 2023
!  Revised by Lu Li, January 24, 2024
!  Revised by Sisi Chen, Lu Li, June, 2024
!-----------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Qsadv
   USE MOD_Namelist
   USE MOD_Const_Physical
   USE MOD_Vars_Global
   USE MOD_UserDefFun

   IMPLICIT NONE

   real(r8), parameter :: SHR_CONST_MWDAIR = 28.966_r8       ! molecular weight dry air [kg/kmole]
   real(r8), parameter :: SHR_CONST_MWWV   = 18.016_r8       ! molecular weight water vapor
   real(r8), parameter :: SHR_CONST_AVOGAD = 6.02214e26_r8   ! Avogadro's number [molecules/kmole]
   real(r8), parameter :: SHR_CONST_BOLTZ  = 1.38065e-23_r8  ! Boltzmann's constant [J/K/molecule]
   ! Universal gas constant [J/K/kmole]
   real(r8), parameter :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ
   ! Dry air gas constant [J/K/kg]
   real(r8), parameter :: rair             = SHR_CONST_RGAS/SHR_CONST_MWDAIR

   ! On the windward side of the range, annual mean lapse rates of 3.9-5.2 (deg km-1), substantially
   ! smaller than the often-assumed 6.5 (deg km-1).  The data sets show similar seasonal and diurnal
   ! variability, with lapse rates smallest (2.5-3.5 deg km-1) in late-summer minimum temperatures,
   ! and largest (6.5-7.5 deg km-1) in spring maximum temperatures. Geographic (windward versus lee
   ! side) differences in lapse rates are substantial. [Minder et al., 2010, Surface temperature
   ! lapse rates over complex terrain: Lessons from the Cascade Mountains. JGR, 115,
   ! doi:10.1029/2009JD013493]
   !
   ! Kunkel, K. E., 1989: Simple procedures for extrapolation of humidity variables in the
   ! mountainous western United States.  J. Climate, 2, 656-669. lapse_rate = /Jan -
   ! Dec/4.4,5.9,7.1,7.8,8.1,8.2,8.1,8.1,7.7,6.8,5.5,4.7/ (deg km-1)
   real(r8), parameter :: lapse_rate = 0.006_r8         ! surface temperature lapse rate (deg m-1)

   SAVE

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: downscale_forcings     ! Downscale atmospheric forcing

! PRIVATE MEMBER FUNCTIONS:
   PRIVATE :: rhos                  ! calculate atmospheric density

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   PURE FUNCTION rhos(qbot, pbot, tbot)

!-----------------------------------------------------------------------------
! DESCRIPTION:
! Compute atmospheric density (kg/m**3)
!-----------------------------------------------------------------------------

   IMPLICIT NONE

   ! ARGUMENTS:
   real(r8) :: rhos              ! function result: atmospheric density (kg/m**3)
   real(r8), intent(in) :: qbot  ! atmospheric specific humidity (kg/kg)
   real(r8), intent(in) :: pbot  ! atmospheric pressure (Pa)
   real(r8), intent(in) :: tbot  ! atmospheric temperature (K)

   ! LOCAL VARIABLES:
   real(r8) :: egcm
   real(r8) :: wv_to_dair_weight_ratio ! ratio of molecular weight of water vapor to dry air [-]

      wv_to_dair_weight_ratio = SHR_CONST_MWWV/SHR_CONST_MWDAIR

      egcm = qbot*pbot / (wv_to_dair_weight_ratio + (1._r8 - wv_to_dair_weight_ratio)*qbot)
      rhos = (pbot - (1._r8 - wv_to_dair_weight_ratio)*egcm) / (rair*tbot)

   END FUNCTION rhos

!-----------------------------------------------------------------------------

   SUBROUTINE downscale_forcings (&
                  glaciers, &

                  ! non-adjusted forcing
                  forc_topo_g ,forc_maxelv_g ,forc_t_g   ,forc_th_g  ,forc_q_g     ,&
                  forc_pbot_g ,forc_rho_g    ,forc_prc_g ,forc_prl_g ,forc_lwrad_g ,&
                  forc_hgt_g  ,forc_swrad_g  ,forc_us_g  ,forc_vs_g  , &

                  ! topography-based factor on patch
                  slp_type_c, asp_type_c, cur_c, &

                  ! other factors
                  julian_day, coszen, cosazi, &

                  ! adjusted forcing
                  forc_topo_c ,forc_t_c   ,forc_th_c  ,forc_q_c     ,forc_pbot_c ,&
                  forc_rho_c  ,forc_prc_c ,forc_prl_c ,forc_lwrad_c, forc_swrad_c, &
                  forc_us_c   ,forc_vs_c, &

                  ! optional parameters for full downscaling
                  area_type_c, svf_c, alb, &
#ifdef SinglePoint
                  sf_lut_c &
#else
                  sf_curve_c &
#endif
                  )

!-----------------------------------------------------------------------------
! !DESCRIPTION:
!  Downscale atmospheric forcing fields.
!
!  Downscaling is done based on the difference between each land model
!  column's elevation and the atmosphere's surface elevation (which is
!  the elevation at which the atmospheric forcings are valid).
!
!  Note that the downscaling procedure can result in changes in grid
!  cell mean values compared to what was provided by the atmosphere. We
!  conserve fluxes of mass and energy, but allow states such as
!  temperature to differ.
!-----------------------------------------------------------------------------

   IMPLICIT NONE

   integer,  parameter :: S = 1370             ! solar constant (W/m**2)
   real(r8), parameter :: thr = 85*PI/180      ! threshold of zenith angle

   ! ARGUMENTS:
   logical,  intent(in) :: glaciers            ! true: glacier column (itypwat = 3)
   real(r8), intent(in) :: julian_day          ! day of year
   real(r8), intent(in) :: coszen              ! cosine of sun zenith angle at an hour
   real(r8), intent(in) :: cosazi              ! cosine of sun azimuth angle at an hour

   ! topography-based factor
   real(r8), intent(in) :: cur_c               ! curvature
   ! topographic aspect of each type of one patch(rad) - can be different dimensions
   real(r8), intent(in) :: asp_type_c (:)
   ! topographic slope of each character of one patch - can be different dimensions  
   real(r8), intent(in) :: slp_type_c (:)

   ! optional parameters for complex downscaling
   real(r8), intent(in), optional :: alb                 ! blue sky albedo
   real(r8), intent(in), optional :: svf_c               ! sky view factor
   ! area percentage of each character of one patch
   real(r8), intent(in), optional :: area_type_c(:)
#ifdef SinglePoint
   ! look up table of shadow mask of a patch
   real(r8), intent(in), optional :: sf_lut_c   (:,:)
#else
   ! curve of shadow mask of a patch
   real(r8), intent(in), optional :: sf_curve_c (:,:)
#endif

   ! non-downscaled fields:
   real(r8), intent(in) :: forc_topo_g   ! atmospheric surface height [m]
   real(r8), intent(in) :: forc_maxelv_g ! max atmospheric surface height [m]
   real(r8), intent(in) :: forc_t_g      ! atmospheric temperature [Kelvin]
   real(r8), intent(in) :: forc_th_g     ! atmospheric potential temperature [Kelvin]
   real(r8), intent(in) :: forc_q_g      ! atmospheric specific humidity [kg/kg]
   real(r8), intent(in) :: forc_pbot_g   ! atmospheric pressure [Pa]
   real(r8), intent(in) :: forc_rho_g    ! atmospheric density [kg/m**3]
   real(r8), intent(in) :: forc_prc_g    ! convective precipitation in grid [mm/s]
   real(r8), intent(in) :: forc_prl_g    ! large-scale precipitation in grid [mm/s]
   real(r8), intent(in) :: forc_lwrad_g  ! grid downward longwave [W/m**2]
   real(r8), intent(in) :: forc_swrad_g  ! grid downward shortwave [W/m**2]
   real(r8), intent(in) :: forc_hgt_g    ! atmospheric reference height [m]
   real(r8), intent(in) :: forc_us_g     ! eastward wind [m/s]
   real(r8), intent(in) :: forc_vs_g     ! northward wind [m/s]

   ! downscaled fields:
   real(r8), intent(in)  :: forc_topo_c  ! column surface height [m]
   real(r8), intent(out) :: forc_t_c     ! atmospheric temperature [Kelvin]
   real(r8), intent(out) :: forc_th_c    ! atmospheric potential temperature [Kelvin]
   real(r8), intent(out) :: forc_q_c     ! atmospheric specific humidity [kg/kg]
   real(r8), intent(out) :: forc_pbot_c  ! atmospheric pressure [Pa]
   real(r8), intent(out) :: forc_rho_c   ! atmospheric density [kg/m**3]
   real(r8), intent(out) :: forc_prc_c   ! column convective precipitation [mm/s]
   real(r8), intent(out) :: forc_prl_c   ! column large-scale precipitation [mm/s]
   real(r8), intent(out) :: forc_lwrad_c ! column downward longwave [W/m**2]
   real(r8), intent(out) :: forc_swrad_c ! column downward shortwave [W/m**2]
   real(r8), intent(out) :: forc_us_c    ! column eastward wind [m/s]
   real(r8), intent(out) :: forc_vs_c    ! column northward wind [m/s]

   ! Local variables for topo downscaling:
   real(r8) :: hsurf_g, hsurf_c
   real(r8) :: Hbot, zbot
   real(r8) :: tbot_g, pbot_g, thbot_g, qbot_g, qs_g, es_g, rhos_g
   real(r8) :: tbot_c, pbot_c, thbot_c, qbot_c, qs_c, es_c, rhos_c
   real(r8) :: rhos_c_estimate, rhos_g_estimate
   real(r8) :: dum1, dum2
   real(r8) :: max_elev_c    ! the maximum column level elevation value within the grid
   real(r8) :: delta_prc_c   ! deviation of the column convective prec. from the grid level prec.
   real(r8) :: delta_prl_c   ! deviation of the column large-scale prec. from the grid level prec.

!-----------------------------------------------------------------------------

      ! --------------------------------------------------------------------------------------
      ! 1. adjust air temperature, potential temperature, specific humidity, pressure, density
      ! --------------------------------------------------------------------------------------
      hsurf_g = forc_topo_g             ! gridcell sfc elevation
      tbot_g  = forc_t_g                ! atm sfc temp
      thbot_g = forc_th_g               ! atm sfc pot temp
      qbot_g  = forc_q_g                ! atm sfc spec humid
      pbot_g  = forc_pbot_g             ! atm sfc pressure
      rhos_g  = forc_rho_g              ! atm density
      zbot    = forc_hgt_g              ! atm ref height

      hsurf_c = forc_topo_c                           ! column sfc elevation
      tbot_c  = tbot_g-lapse_rate*(hsurf_c-hsurf_g)   ! adjust [temp] for column
      Hbot    = rair*0.5_r8*(tbot_g+tbot_c)/grav      ! scale ht at avg temp
      pbot_c  = pbot_g*exp(-(hsurf_c-hsurf_g)/Hbot)   ! adjust [press] for column

      ! Derivation of potential temperature calculation:
      !
      ! The textbook definition would be:
      ! thbot_c = tbot_c * (p0/pbot_c)^(rair/cpair)
      !
      ! Note that pressure is related to scale height as:
      ! pbot_c = p0 * exp(-zbot/Hbot)
      !
      ! Plugging this in to the textbook definition, then manipulating, we get:
      ! thbot_c = tbot_c * (p0/(p0*exp(-zbot/Hbot)))^(rair/cpair)
      !         = tbot_c * (1/exp(-zbot/Hbot))^(rair/cpair)
      !         = tbot_c * (exp(zbot/Hbot))^(rair/cpair)
      !         = tbot_c * exp((zbot/Hbot) * (rair/cpair))

      ! But we want everything expressed in delta form, resulting in:
      ! adjust [pot temp] for column
      thbot_c = thbot_g + (tbot_c - tbot_g)*exp((zbot/Hbot)*(rair/cpair))

      CALL Qsadv(tbot_g,pbot_g,es_g,dum1,qs_g,dum2) ! es, qs for gridcell
      CALL Qsadv(tbot_c,pbot_c,es_c,dum1,qs_c,dum2) ! es, qs for column
      qbot_c = qbot_g*(qs_c/qs_g) ! adjust [q] for column

      rhos_c_estimate = rhos(qbot=qbot_c, pbot=pbot_c, tbot=tbot_c)
      rhos_g_estimate = rhos(qbot=qbot_g, pbot=pbot_g, tbot=tbot_g)
      rhos_c = rhos_g * (rhos_c_estimate / rhos_g_estimate) ! adjust [density] for column

      ! save
      forc_t_c    = tbot_c
      forc_th_c   = thbot_c
      forc_q_c    = qbot_c
      forc_pbot_c = pbot_c
      forc_rho_c  = rhos_c

      ! --------------------------------------------------------------------------------------
      ! 2. adjust wind speed
      ! --------------------------------------------------------------------------------------
      forc_us_c = forc_us_g
      forc_vs_c = forc_vs_g

      ! --------------------------------------------------------------------------------------
      ! 3. adjust longwave radiation and shortwave radiation
      ! --------------------------------------------------------------------------------------
      CALL downscale_longwave (glaciers, &
                     forc_topo_g, forc_t_g, forc_q_g, forc_pbot_g, forc_lwrad_g, &
                     forc_topo_c, forc_t_c, forc_q_c, forc_pbot_c, forc_lwrad_c)

      ! Check if optional parameters are present for complex downscaling
      IF (present(area_type_c) .and. present(svf_c) .and. present(alb) .and. &
#ifdef SinglePoint
          present(sf_lut_c) &
#else
          present(sf_curve_c) &
#endif
          ) THEN
         ! Complex downscaling with topographic effects
         CALL downscale_shortwave(&
                        forc_topo_g, forc_pbot_g, forc_swrad_g, &
                        forc_topo_c, forc_pbot_c, forc_swrad_c, &
                        julian_day, coszen, cosazi, alb, &
                        slp_type_c, asp_type_c, svf_c,   &
#ifdef SinglePoint
                        sf_lut_c,   &
#else
                        sf_curve_c, &
#endif
                        area_type_c)
      ELSE
         ! Simple downscaling
         CALL downscale_shortwave_simple(&
                        forc_topo_g, forc_pbot_g, forc_swrad_g, &
                        forc_topo_c, forc_pbot_c, forc_swrad_c, &
                        julian_day, coszen, cosazi, &
                        slp_type_c, asp_type_c)
      ENDIF

      ! --------------------------------------------------------------------------------------
      ! 4. adjust precipitation
      ! --------------------------------------------------------------------------------------
      IF (trim(DEF_DS_precipitation_adjust_scheme) == 'I') THEN
         ! Tesfa et al, 2020: Exploring Topography-Based Methods for Downscaling
         ! Subgrid Precipitation for Use in Earth System Models. Equation (5)
         ! https://doi.org/ 10.1029/2019JD031456

         IF (forc_maxelv_g /= 0.) THEN
            delta_prc_c = forc_prc_g * (forc_topo_c - forc_topo_g) / forc_maxelv_g
         ELSE
            delta_prc_c = 0.
         ENDIF
         forc_prc_c  = forc_prc_g + delta_prc_c     ! convective precipitation [mm/s]
         IF (forc_prc_c<=0) forc_prc_c = forc_prc_g ! the limit value is non-negative
         IF (forc_prc_c==0) forc_prc_c = 1.0e-10_r8 ! avoid denominator being 0 when conserving water quantity

         IF (forc_maxelv_g /= 0.) THEN
            delta_prl_c = forc_prl_g * (forc_topo_c - forc_topo_g) / forc_maxelv_g
         ELSE
            delta_prl_c = 0.
         ENDIF
         forc_prl_c  = forc_prl_g + delta_prl_c    ! large scale precipitation [mm/s]
         IF (forc_prl_c<=0) forc_prl_c = forc_prl_g ! the limit value is non-negative
         IF (forc_prl_c==0) forc_prl_c = 1.0e-10_r8 ! avoid denominator being 0 when conserving water quantity

      ELSEIF (trim(DEF_DS_precipitation_adjust_scheme) == 'II') THEN
         ! Liston, G. E. and Elder, K.: A meteorological distribution system
         ! for high-resolution terrestrial modeling (MicroMet), J. Hydrometeorol., 7, 217-234, 2006.
         ! Equation (33) and Table 1: chi range from January to December:
         ! [0.35,0.35,0.35,0.30,0.25,0.20,0.20,0.20,0.20,0.25,0.30,0.35] (1/km)

         delta_prc_c = forc_prc_g *1.0*0.27*(forc_topo_c - forc_topo_g) &
            /(1.0 - 0.27*(forc_topo_c - forc_topo_g))
         forc_prc_c = forc_prc_g + delta_prc_c   ! large scale precipitation [mm/s]

         delta_prl_c = forc_prl_g *1.0*0.27*(forc_topo_c - forc_topo_g) &
            /(1.0 - 0.27*(forc_topo_c - forc_topo_g))
         forc_prl_c = forc_prl_g + delta_prl_c   ! large scale precipitation [mm/s]
      ENDIF

      IF (forc_prl_c < 0) THEN
         write(*,*) 'negative prl', forc_prl_g, forc_maxelv_g, forc_topo_c, forc_topo_g
         forc_prl_c = 0.
      ENDIF

      IF (forc_prc_c < 0) THEN
         write(*,*) 'negative prc', forc_prc_g, forc_maxelv_g, forc_topo_c, forc_topo_g
         forc_prc_c = 0.
      ENDIF

   END SUBROUTINE downscale_forcings
!-----------------------------------------------------------------------------

   SUBROUTINE downscale_wind(forc_us_g, forc_vs_g, &
                             slp_type_c, asp_type_c, area_type_c, cur_c)

!-----------------------------------------------------------------------------
! !DESCRIPTION:
!  Downscale wind speed
!
!  Liston, G. E. and Elder, K.: A meteorological distribution system for
!  high-resolution terrestrial modeling (MicroMet), J. Hydrometeorol.,
!  7, 217-234, 2006.
!-----------------------------------------------------------------------------

   IMPLICIT NONE

   ! ARGUMENTS:
   real(r8), intent(inout)  :: forc_us_g         ! eastward wind (m/s)
   real(r8), intent(inout)  :: forc_vs_g         ! northward wind (m/s)

   real(r8), intent(in) :: cur_c                 ! curvature
   ! topographic aspect of each character of one patch
   real(r8), intent(in) :: asp_type_c  (1:num_slope_type)
   ! topographic slope of each character of one patch
   real(r8), intent(in) :: slp_type_c  (1:num_slope_type)
   ! area percentage of each character of one patch
   real(r8), intent(in) :: area_type_c (1:num_slope_type)

   ! local variables
   real(r8) :: wind_dir                          ! wind direction
   real(r8) :: ws_g                              ! non-downscaled wind speed
   real(r8) :: wind_dir_slp (1:num_slope_type)   ! the slope in the direction of the wind
   real(r8) :: ws_c_type(1:num_slope_type)       ! downscaled wind speed of each type in each patch
   real(r8) :: ws_c                              ! downscaled wind speed
   real(r8) :: scale_factor                      ! Combined scaling factor for regulating wind speed
   integer :: g, c, i

!-----------------------------------------------------------------------------

      ! calculate wind direction
      IF (forc_us_g == 0.) THEN
         wind_dir  = PI/2
      ELSE
         wind_dir  = atan(forc_vs_g /forc_us_g)
      ENDIF

      ! non-adjusted wind speed
      ws_g  = sqrt(forc_vs_g *forc_vs_g +forc_us_g *forc_us_g )

      ! compute the slope in the direction of the wind
      DO i = 1, num_slope_type
         wind_dir_slp(i) = slp_type_c(i)*cos(wind_dir-asp_type_c(i))
      ENDDO

      ! compute wind speed adjustment
      DO i = 1, num_slope_type
         scale_factor = (1+(0.58*wind_dir_slp(i))+0.42*cur_c)
         ! Limiting the scope of proportionality adjustments
         IF (scale_factor>1.5) THEN
             scale_factor = 1.5
         ELSEIF (scale_factor<-1.5) THEN
             scale_factor = -1.5
         ENDIF
         ws_c_type(i) = ws_g *scale_factor*area_type_c(i)
      ENDDO

      ! adjusted wind speed
      ws_c = sum(ws_c_type(:))
      forc_us_g = ws_c*cos(wind_dir)
      forc_vs_g = ws_c*sin(wind_dir)

   END SUBROUTINE downscale_wind

!-----------------------------------------------------------------------------

   SUBROUTINE downscale_wind_simple(forc_us_g, forc_vs_g, &
                                    slp_type_c, area_type_c, cur_c)

!-----------------------------------------------------------------------------
! !DESCRIPTION:
!  Downscale wind speed using a simple method for land-atmosphere coupling models
!-----------------------------------------------------------------------------

   IMPLICIT NONE

   ! ARGUMENTS:
   real(r8), intent(inout)  :: forc_us_g         ! eastward wind (m/s)
   real(r8), intent(inout)  :: forc_vs_g         ! northward wind (m/s)

   real(r8), intent(in) :: cur_c                 ! curvature
   ! topographic slope of each character of one patch
   real(r8), intent(in) :: slp_type_c  (1:num_aspect_type)
   ! area percentage of each character of one patch
   real(r8), intent(in) :: area_type_c (1:num_aspect_type)

   ! local variables
   real(r8) :: asp_type_c(1:num_aspect_type)        ! topographic aspect of each type of one patch (rad)
   real(r8) :: slp_type_c_rad                       ! Convert tan slope value to slope angle value

   real(r8) :: wind_dir                             ! wind direction
   real(r8) :: ws_g                                 ! non-downscaled wind speed
   real(r8) :: wind_dir_slp (1:num_aspect_type)     ! the slope in the direction of the wind
   real(r8) :: ws_c_type(1:num_aspect_type)         ! downscaled wind speed of each type in each patch
   real(r8) :: ws_c                                 ! downscaled wind speed
   real(r8) :: scale_factor                         ! Combined scaling factor for regulating wind speed
   integer :: u_sign, v_sign, i, wind_dir_u, wind_dir_v

   ! Initialize aspect type
   asp_type_c(1) = 0.0_r8*PI/180     ! north
   asp_type_c(2) = 45.0_r8*PI/180    ! northeast
   asp_type_c(3) = 90.0_r8*PI/180    ! east
   asp_type_c(4) = 135.0_r8*PI/180   ! southeast
   asp_type_c(5) = 180.0_r8*PI/180   ! south
   asp_type_c(6) = 225.0_r8*PI/180   ! southwest
   asp_type_c(7) = 270.0_r8*PI/180   ! west
   asp_type_c(8) = 315.0_r8*PI/180   ! northwest
   asp_type_c(9) = -9999.0_r8        ! flat
!-----------------------------------------------------------------------------

      ! calculate wind direction
      IF (forc_us_g == 0.) THEN
         wind_dir  = PI/2
      ELSE
         wind_dir  = atan2(forc_vs_g, forc_us_g)
      ENDIF

      ! convert to 0-2*PI range
      !wind_dir = wind_dir + 3*PI/2 
      !wind_dir = mod(wind_dir, 2*PI)

      ! 0° is north 
      IF (wind_dir > PI/2) THEN
         wind_dir = 2*PI - wind_dir + PI/2
      ELSE
         wind_dir = PI/2 - wind_dir
      ENDIF

      ! non-adjusted wind speed
      ws_g  = sqrt(forc_vs_g *forc_vs_g +forc_us_g *forc_us_g )

      ! log the + - sign of the u,v direction
      IF (forc_us_g >= 0.0_r8) THEN
         u_sign = 1
      ELSE
         u_sign = -1
      ENDIF 
      IF (forc_vs_g >= 0.0_r8) THEN
         v_sign = 1
      ELSE
         v_sign = -1
      ENDIF


      ! compute the slope in the direction of the wind
      DO i = 1, num_aspect_type
         IF (slp_type_c(i) == -1.0e36) THEN
            ! no slope in the direction of the aspect
            wind_dir_slp(i) = -1.0e36
         ELSE
            ! slope in the direction of the wind
            slp_type_c_rad = atan(slp_type_c(i))
            wind_dir_slp(i) = slp_type_c_rad*cos(wind_dir-asp_type_c(i))
         ENDIF
      ENDDO

      ! compute wind speed adjustment
      DO i = 1, num_aspect_type
         ! For the flat area, we do not adjust the wind speed
         IF (asp_type_c(i) == -9999.0_r8) THEN
            scale_factor = 1.0_r8
         ELSE IF ((wind_dir_slp(i) == -1.0e36).or.(cur_c == -1.0e36)) THEN
            ! no slope in the direction of the aspect
            scale_factor = -1.0e36
         ELSE
            scale_factor = (1+(0.58*wind_dir_slp(i))+0.42*cur_c)
            !write(*,*) 'scale_factor', scale_factor, 'wind_dir_slp(i)', wind_dir_slp(i), 'cur_c', cur_c
         ENDIF

         ! Limiting the scope of proportionality adjustments
         IF (scale_factor>1.5) THEN
             scale_factor = 1.5
         ELSEIF (scale_factor<-1.5) THEN
             scale_factor = -1.5
         ENDIF

         ! Downscale wind speed for each type in each patch
         IF ((scale_factor == -1.0e36).or.(area_type_c(i) == -1.0e36)) THEN
            ws_c_type(i) = -1.0e36
         ELSE
            ws_c_type(i) = ws_g *scale_factor*area_type_c(i)
         ENDIF

      ENDDO

      ! adjusted wind speed
      ws_c = sum(ws_c_type(:),mask=ws_c_type(:) /= -1.0e36)

      ! caculate u and v components of the wind
      IF (wind_dir > PI/2) THEN
         wind_dir = 2*PI - wind_dir + PI/2
      ELSE
         wind_dir = PI/2 - wind_dir
      ENDIF
      
      forc_us_g = u_sign*sqrt((ws_c*cos(wind_dir))**2)
      forc_vs_g = v_sign*sqrt((ws_c*sin(wind_dir))**2)

   END SUBROUTINE downscale_wind_simple
!-----------------------------------------------------------------------------

   SUBROUTINE downscale_longwave (glaciers, &
      forc_topo_g, forc_t_g, forc_q_g, forc_pbot_g, forc_lwrad_g, &
      forc_topo_c, forc_t_c, forc_q_c, forc_pbot_c, forc_lwrad_c)

!-----------------------------------------------------------------------------
! !DESCRIPTION:
!  Downscale longwave radiation
!-----------------------------------------------------------------------------

   IMPLICIT NONE

   ! ARGUMENTS:
   logical,  intent(in) :: glaciers      ! true: glacier column

   real(r8), intent(in) :: forc_topo_g   ! atmospheric surface height (m)
   real(r8), intent(in) :: forc_t_g      ! atmospheric temperature [Kelvin]
   real(r8), intent(in) :: forc_q_g      ! atmospheric specific humidity [kg/kg]
   real(r8), intent(in) :: forc_pbot_g   ! atmospheric pressure [Pa]
   real(r8), intent(in) :: forc_lwrad_g  ! downward longwave (W/m**2)

   real(r8), intent(in) :: forc_topo_c   ! column surface height (m)
   real(r8), intent(in) :: forc_t_c      ! atmospheric temperature [Kelvin]
   real(r8), intent(in) :: forc_q_c      ! atmospheric specific humidity [kg/kg]
   real(r8), intent(in) :: forc_pbot_c   ! atmospheric pressure [Pa]
   real(r8), intent(out):: forc_lwrad_c  ! downward longwave (W/m**2)

   ! LOCAL VARIABLES:
   real(r8) :: hsurf_c  ! column-level elevation (m)
   real(r8) :: hsurf_g  ! gridcell-level elevation (m)

   real(r8) :: pv_g ! the water vapor pressure at grid cell (hPa)
   real(r8) :: pv_c ! the water vapor pressure at column (hPa)
   real(r8) :: emissivity_clearsky_g ! clear-sky emissivity at grid cell
   real(r8) :: emissivity_clearsky_c ! clear-sky emissivity at grid column
   real(r8) :: emissivity_allsky_g   ! all-sky emissivity at grid cell
   real(r8) :: es_g, es_c, dum1, dum2, dum3

   real(r8), parameter :: lapse_rate_longwave = 0.032_r8 ! longwave radiation lapse rate (W m-2 m-1)
   ! relative limit for how much longwave downscaling can be done (unitless)
   real(r8), parameter :: longwave_downscaling_limit = 0.5_r8

!--------------------------------------------------------------------------

      ! Initialize (needs to be done for ALL active columns)
      forc_lwrad_c = forc_lwrad_g
      hsurf_g = forc_topo_g
      hsurf_c = forc_topo_c

      IF (trim(DEF_DS_longwave_adjust_scheme) == 'I') THEN
         ! Fiddes and Gruber, 2014, TopoSCALE v.1.0: downscaling gridded climate data in complex
         ! terrain. Geosci. Model Dev., 7, 387-405. doi:10.5194/gmd-7-387-2014.  Equation (1) (2)
         ! (3); here, the empirical parameters x1 and x2 are different from Konzelmann et al. (1994)
         ! where x1 = 0.443 and x2 = 8 (optimal for measurements on the Greenland ice sheet)

         CALL Qsadv(forc_t_g, forc_pbot_g, es_g,dum1,dum2,dum3)
         CALL Qsadv(forc_t_c, forc_pbot_c, es_c,dum1,dum2,dum3)
         pv_g = forc_q_g*es_g/100._r8  ! (hPa)
         pv_c = forc_q_c*es_c/100._r8  ! (hPa)

         emissivity_clearsky_g = 0.23_r8 + 0.43_r8*(pv_g/forc_t_g)**(1._r8/5.7_r8)
         emissivity_clearsky_c = 0.23_r8 + 0.43_r8*(pv_c/forc_t_c)**(1._r8/5.7_r8)
         emissivity_allsky_g = forc_lwrad_g / (5.67e-8_r8*forc_t_g**4)

         forc_lwrad_c = &
            (emissivity_clearsky_c + (emissivity_allsky_g - emissivity_clearsky_g)) &
            * 5.67e-8_r8*forc_t_c**4
      ELSE
         ! Longwave radiation is downscaled by assuming a linear decrease in downwelling longwave
         ! radiation with increasing elevation (0.032 W m-2 m-1, limited to 0.5 - 1.5 times the
         ! gridcell mean value, then normalized to conserve gridcell total energy) (Van Tricht et
         ! al., 2016, TC) Figure 6, doi:10.5194/tc-10-2379-2016

         IF (glaciers) THEN
            forc_lwrad_c = forc_lwrad_g - lapse_rate_longwave * (hsurf_c-hsurf_g)

            ! Here we assume that deltaLW = (dLW/dT)*(dT/dz)*deltaz
            ! We get dLW/dT = 4*eps*sigma*T^3 = 4*LW/T from the Stefan-Boltzmann law,
            ! evaluated at the mean temp. We assume the same temperature lapse rate as above.

         ELSE
            forc_lwrad_c = forc_lwrad_g &
               - 4.0_r8 * forc_lwrad_g/(0.5_r8*(forc_t_c+forc_t_g)) &
               * lapse_rate * (hsurf_c - hsurf_g)
         ENDIF
      ENDIF

      ! But ensure that we don't depart too far from the atmospheric forcing value:
      ! negative values of lwrad are certainly bad, but small positive values might
      ! also be bad. We can especially run into trouble due to the normalization: a
      ! small lwrad value in one column can lead to a big normalization factor,
      ! leading to huge lwrad values in other columns.

      forc_lwrad_c = min(forc_lwrad_c, forc_lwrad_g * (1._r8 + longwave_downscaling_limit))
      forc_lwrad_c = max(forc_lwrad_c, forc_lwrad_g * (1._r8 - longwave_downscaling_limit))

   END SUBROUTINE downscale_longwave

!-----------------------------------------------------------------------------   
   SUBROUTINE downscale_shortwave( &
                        forc_topo_g, forc_pbot_g, forc_swrad_g, &
                        forc_topo_c, forc_pbot_c, forc_swrad_c, &
                        julian_day, coszen, cosazi, alb, &
                        slp_type_c, asp_type_c, svf_c,   &
#ifdef SinglePoint
                        sf_lut_c,   &
#else
                        sf_curve_c, &
#endif
                        area_type_c)

!-----------------------------------------------------------------------------
! !DESCRIPTION:
!
!  Rouf, T., Mei, Y., Maggioni, V., Houser, P., & Noonan, M. (2020). A
!  Physically Based Atmospheric Variables Downscaling Technique. Journal
!  of Hydrometeorology, 21(1), 93-108.
!  https://doi.org/10.1175/JHM-D-19-0109.1
!
!  Sisi Chen, Lu Li, Yongjiu Dai, et al. Exploring Topography
!  Downscaling Methods for Hyper-Resolution Land Surface Modeling.
!  Authorea. April 25, 2024.  DOI: 10.22541/au.171403656.68476353/v1
!
!  Must be done after downscaling of surface pressure
!-----------------------------------------------------------------------------

   IMPLICIT NONE

   integer,  parameter :: S = 1370                               ! solar constant (W/m**2)
   real(r8), parameter :: thr = 85*PI/180                        ! threshold of zenith
   ! relative limit for how much shortwave downscaling can be done (unitless)
   real(r8), parameter :: shortwave_downscaling_limit = 0.5_r8

   ! ARGUMENTS:
   real(r8), intent(in) :: julian_day          ! day of year
   real(r8), intent(in) :: coszen              ! zenith angle at an hour
   real(r8), intent(in) :: cosazi              ! azimuth angle at an hour
   real(r8), intent(in) :: alb                 ! blue sky albedo

   real(r8), intent(in) :: forc_topo_g         ! atmospheric surface height (m)
   real(r8), intent(in) :: forc_pbot_g         ! atmospheric pressure [Pa]
   real(r8), intent(in) :: forc_swrad_g        ! downward shortwave (W/m**2)

   real(r8), intent(in) :: forc_topo_c         ! column surface height (m)
   real(r8), intent(in) :: forc_pbot_c         ! atmospheric pressure [Pa]
   real(r8), intent(out):: forc_swrad_c        ! downward shortwave (W/m**2)

   real(r8), intent(in) :: svf_c               ! sky view factor
# ifdef SinglePoint
   ! look up table of shadow mask of a patch
   real(r8), intent(in) :: sf_lut_c   (1:num_azimuth,1:num_zenith)
# else
   ! curve of shadow mask of a patch
   real(r8), intent(in) :: sf_curve_c (1:num_azimuth,1:num_zenith_parameter)
# endif
   ! topographic aspect of each character of one patch (°)
   real(r8), intent(in) :: asp_type_c (1:num_slope_type)
   ! topographic slope of each character of one patch
   real(r8), intent(in) :: slp_type_c (1:num_slope_type)
   ! area percentage of each character of one patch
   real(r8), intent(in) :: area_type_c(1:num_slope_type)

   ! LOCAL VARIABLES:
   real(r8) :: zen_rad, zen_deg, azi_rad, azi_deg ! rad and deg of sun zenith and azimuth angle
   integer  :: idx_azi, idx_zen                   ! index used to cal shadow factor from
                                                  ! look up table
   real(r8) :: sf_c                               ! shadow factor
   real(r8) :: rt_R                               ! The ratio of the current distance between
                                                  ! the sun and the earth
   real(r8) :: toa_swrad                          ! top of atmosphere shortwave radiation
   real(r8) :: clr_idx                            ! atmospheric transparency
   real(r8) :: diff_wgt                           ! diffuse weight
   real(r8) :: k_c                                ! column broadband attenuation coefficient [Pa^-1]
   real(r8) :: opt_factor                         ! optical length factor
   real(r8) :: a_p
   real(r8) :: svf, balb

   real(r8) :: diff_swrad_g, beam_swrad_g              ! diffuse and beam radiation
   real(r8) :: diff_swrad_c, beam_swrad_c, refl_swrad_c! downscaled diffuse, beam radiation
                                                       ! and reflect radiation
   real(r8) :: beam_swrad_type (1:num_slope_type)      ! beam radiation of one characterized patch
   real(r8) :: refl_swrad_type (1:num_slope_type)      ! refl. radiation of one characterized patch
   real(r8) :: tcf_type        (1:num_slope_type)      ! terrain configure factor
   real(r8) :: cosill_type     (1:num_slope_type)      ! illumination angle (cos) at defined types

   real(r8) :: zenith_segment, a1, a2                  ! Segmented function segmentation
                                                       ! points (rad), parameter1, parameter2

   integer  :: i

!-----------------------------------------------------------------------------

      ! calculate shadow factor according to sun zenith and azimuth angle
      zen_rad = acos(coszen)
      azi_rad = acos(cosazi)
      azi_deg = azi_rad*180.0/PI ! turn deg

      idx_azi = INT(azi_deg*num_azimuth/360)

      IF (idx_azi==0) idx_azi = 1

#ifdef SinglePoint
      zen_deg = zen_rad*180/PI ! turn deg
      idx_zen = INT(zen_deg*num_zenith/90)
      IF (idx_zen==0) idx_zen = 1
      !constrain the upper boundary of zenith angle to 90 deg
      IF (idx_zen>num_zenith) idx_zen = num_zenith

      sf_c = sf_lut_c(idx_azi, idx_zen)
#else
      ! Constructing a shadow factor function from zenith angle parameters
      ! shadow factor = exp(-1*exp(a1*zenith+a2))
      zenith_segment = sf_curve_c(idx_azi, 1)       ! Segmented function segmentation points (rad)
      a1 = sf_curve_c(idx_azi, 2)                   ! parameter of function
      a2 = sf_curve_c(idx_azi, 3)                   ! parameter of function

      IF (zen_rad <= zenith_segment) THEN
         sf_c = 1.
      ELSEIF (a1<=1e-10) THEN
         sf_c = 1.
      ELSE
         sf_c = exp(-1*exp(min(a1*zen_rad+a2,3.5)))
      ENDIF
#endif

      IF (sf_c<0) sf_c = 0
      IF (sf_c>1) sf_c = 1

      ! calculate top-of-atmosphere incident shortwave radiation
      rt_R = 1-0.01672*cos(0.9856*(julian_day-4))
      toa_swrad = S*(rt_R**2)*coszen

      ! calculate clearness index
      IF (toa_swrad.eq.0) THEN
         clr_idx = 0
      ELSE
         clr_idx = forc_swrad_g/toa_swrad
      ENDIF
      IF (clr_idx>1) clr_idx = 1

      ! calculate diffuse weight
      ! Ruiz-Arias, J. A., Alsamamra, H., Tovar-Pescador, J., & Pozo-Vázquez, D. (2010).
      ! Proposal of a regressive model for the hourly diffuse solar radiation under all sky
      ! conditions. Energy Conversion and Management, 51(5), 881-893.
      ! https://doi.org/10.1016/j.enconman.2009.11.024
      diff_wgt = 0.952-1.041*exp(-1*exp(min(2.3-4.702*clr_idx,3.5)))
      IF (diff_wgt>1) diff_wgt = 1
      IF (diff_wgt<0) diff_wgt = 0

      ! calculate diffuse and beam radiation
      diff_swrad_g = forc_swrad_g*diff_wgt
      beam_swrad_g = forc_swrad_g*(1-diff_wgt)

      ! calculate broadband attenuation coefficient [Pa^-1]
      IF (clr_idx.le.0) THEN
         k_c = 0
      ELSE
         k_c = log(clr_idx)/forc_pbot_c
      ENDIF

      ! calculate factor to account for the difference of optical path length
      ! due to pressure difference
      opt_factor = exp(k_c*(forc_pbot_g-forc_pbot_c))
      ! Control the boundary of optical path length
      IF ((opt_factor>10000).or.(opt_factor<-10000)) opt_factor = 0

      ! Adjust the zenith angle so that the range of zenith angles is less than 85°
      IF (zen_rad>thr) zen_rad=thr

      ! loop for four defined types to downscale beam radiation
      DO i = 1, num_slope_type
         ! calculate the cosine of solar illumination angle, cos(θ),
         ! ranging between −1 and 1, indicates if the sun is below or
         ! above the local horizon (note that values lower than 0 are set to 0 indicate self shadow)
         cosill_type(i) = cos(slp_type_c(i))+tan(zen_rad)*sin(slp_type_c(i))*cos(asp_type_c(i))
         IF (cosill_type(i)>1) cosill_type(i) = 1
         IF (cosill_type(i)<0) cosill_type(i) = 0

         ! downscaling beam radiation
         a_p = area_type_c(i)
         IF (a_p.gt.1.0) a_p = 1
         IF (a_p.lt.0) a_p = 0
         beam_swrad_type(i) = sf_c*cosill_type(i)*opt_factor*a_p*beam_swrad_g
      ENDDO
      beam_swrad_c = sum(beam_swrad_type)

      ! downscaling diffuse radiation
      svf = svf_c
      IF (svf>1) svf = 1
      IF (svf<0) svf = 0
      diff_swrad_c = svf*diff_swrad_g

      ! downscaling reflected radiation
      balb = alb
      DO i = 1, num_slope_type
         tcf_type(i) = (1+cos(slp_type_c(i)))/2-svf
         IF (tcf_type(i)<0) tcf_type(i) = 0

         IF (isnan_ud(alb)) THEN
            refl_swrad_type(i) = -1.0e36
         ELSE
            IF ((balb<0).or.(balb>1)) balb = 0
            refl_swrad_type(i) = balb*tcf_type(i)*(beam_swrad_c*coszen+(1-svf)*diff_swrad_c)
         ENDIF
      ENDDO
      refl_swrad_c = sum(refl_swrad_type, mask = refl_swrad_type /= -1.0e36)
      forc_swrad_c = beam_swrad_c+diff_swrad_c+refl_swrad_c

      ! But ensure that we don't depart too far from the atmospheric forcing value:
      ! negative values of swrad are certainly bad, but small positive values might
      ! also be bad. We can especially run into trouble due to the normalization: a
      ! small swrad value in one column can lead to a big normalization factor,
      ! leading to huge swrad values in other columns.

      forc_swrad_c = min(forc_swrad_c, &
               forc_swrad_g * (1._r8 + shortwave_downscaling_limit))
      forc_swrad_c = max(forc_swrad_c, &
               forc_swrad_g * (1._r8 - shortwave_downscaling_limit))
      ! Ensure that the denominator is not 0 during shortwave normalization
      IF (forc_swrad_c==0.) forc_swrad_c = 0.0001

   END SUBROUTINE downscale_shortwave

!-----------------------------------------------------------------------------
   SUBROUTINE downscale_shortwave_simple( &
                        forc_topo_g, forc_pbot_g, forc_swrad_g, &
                        forc_topo_c, forc_pbot_c, forc_swrad_c, &
                        julian_day, coszen, cosazi, &
                        slp_type_c, area_type_c)
!-----------------------------------------------------------------------------
! !DESCRIPTION:
! This subroutine performs a simple downscaling of shortwave radiation for
! land-atmosphere coupling models. The adjustments are only made for direct
! radiation without considering the impact of shadow factor.
!
   IMPLICIT NONE

   integer,  parameter :: S = 1370                               ! solar constant (W/m**2)
   real(r8), parameter :: thr = 85*PI/180                        ! threshold of zenith
   ! relative limit for how much shortwave downscaling can be done (unitless)
   real(r8), parameter :: shortwave_downscaling_limit = 0.2_r8

   ! ARGUMENTS:
   real(r8), intent(in) :: julian_day          ! day of year
   real(r8), intent(in) :: coszen              ! zenith angle at an hour
   real(r8), intent(in) :: cosazi              ! azimuth angle at an hour

   real(r8), intent(in) :: forc_topo_g         ! atmospheric surface height (m)
   real(r8), intent(in) :: forc_pbot_g         ! atmospheric pressure [Pa]
   real(r8), intent(in) :: forc_swrad_g        ! downward shortwave (W/m**2)

   real(r8), intent(in) :: forc_topo_c         ! column surface height (m)
   real(r8), intent(in) :: forc_pbot_c         ! atmospheric pressure [Pa]
   real(r8), intent(out):: forc_swrad_c        ! downward shortwave (W/m**2)

   ! tan value of topographic slope of each direction of one patch
   real(r8), intent(in) :: slp_type_c (1:num_aspect_type)
   ! area percentage of each character of one patch
   real(r8), intent(in) :: area_type_c(1:num_aspect_type)

   ! LOCAL VARIABLES:

   real(r8) :: asp_type_c (1:num_aspect_type)     ! topographic aspect fraction of one patch (%100) num_aspect_type = 1:north, 2:northeast, 3:east, 
                                                  ! 4:southeast, 5:south, 6:southwest, 7:west, 8:northwest, 9:flat
   real(r8) :: slp_type_c_rad                     ! Convert tan slope value to slope angle value

   real(r8) :: zen_rad, azi_rad                   ! rad of sun zenith and azimuth angle
   integer  :: idx_azi, idx_zen                   ! index used to cal shadow factor from
                                                  ! look up table
   real(r8) :: sf_c                               ! shadow factor
   real(r8) :: rt_R                               ! The ratio of the current distance between
                                                  ! the sun and the earth
   real(r8) :: toa_swrad                          ! top of atmosphere shortwave radiation
   real(r8) :: clr_idx                            ! atmospheric transparency
   real(r8) :: diff_wgt                           ! diffuse weight
   real(r8) :: k_c                                ! column broadband attenuation coefficient [Pa^-1]
   real(r8) :: opt_factor                         ! optical length factor
   real(r8) :: a_p
   real(r8) :: diff_swrad_g, beam_swrad_g              ! diffuse and beam radiation
   real(r8) :: diff_swrad_c, beam_swrad_c              ! downscaled diffuse, beam radiation
                                                       
   real(r8) :: beam_swrad_type (1:num_aspect_type)      ! beam radiation of one characterized patch
   real(r8) :: cosill_type     (1:num_aspect_type)      ! illumination angle (cos) at defined types

   integer  :: i
!------------------------------------------------------------------------------

      ! Initialize aspect type
      asp_type_c(1) = 0.0_r8*PI/180     ! north
      asp_type_c(2) = 45.0_r8*PI/180    ! northeast
      asp_type_c(3) = 90.0_r8*PI/180    ! east
      asp_type_c(4) = 135.0_r8*PI/180   ! southeast
      asp_type_c(5) = 180.0_r8*PI/180   ! south
      asp_type_c(6) = 225.0_r8*PI/180   ! southwest
      asp_type_c(7) = 270.0_r8*PI/180   ! west
      asp_type_c(8) = 315.0_r8*PI/180   ! northwest
      asp_type_c(9) = -9999.0_r8        ! flat

      ! calculate shadow factor according to sun zenith and azimuth angle
      zen_rad = acos(coszen)
      azi_rad = acos(cosazi)
      
      ! calculate top-of-atmosphere incident shortwave radiation
      rt_R = 1-0.01672*cos(0.9856*(julian_day-4))
      toa_swrad = S*(rt_R**2)*coszen

      ! calculate clearness index
      IF (toa_swrad < 1.e-7) THEN
         clr_idx = 0
      ELSE
         clr_idx = forc_swrad_g/toa_swrad
      ENDIF
      IF (clr_idx>1) clr_idx = 1

      ! calculate diffuse weight
      ! Ruiz-Arias, J. A., Alsamamra, H., Tovar-Pescador, J., & Pozo-Vázquez, D. (2010).
      ! Proposal of a regressive model for the hourly diffuse solar radiation under all sky
      ! conditions. Energy Conversion and Management, 51(5), 881-893.
      ! https://doi.org/10.1016/j.enconman.2009.11.024
      diff_wgt = 0.952-1.041*exp(-1*exp(min(2.3-4.702*clr_idx,3.5)))
      IF (diff_wgt>1) diff_wgt = 1
      IF (diff_wgt<0) diff_wgt = 0

      ! calculate diffuse and beam radiation
      diff_swrad_g = forc_swrad_g*diff_wgt
      beam_swrad_g = forc_swrad_g*(1-diff_wgt)

      ! calculate broadband attenuation coefficient [Pa^-1]
      IF (clr_idx.le.0) THEN
         k_c = 0
      ELSE
         k_c = log(clr_idx)/forc_pbot_c
      ENDIF

      ! calculate factor to account for the difference of optical path length
      ! due to pressure difference
      opt_factor = exp(k_c*(forc_pbot_g-forc_pbot_c))
      ! Control the boundary of optical path length
      IF ((opt_factor>10000).or.(opt_factor<-10000)) opt_factor = 0

      ! Adjust the zenith angle so that the range of zenith angles is less than 85°
      IF (zen_rad>thr) zen_rad = thr

      ! loop for four defined types to downscale beam radiation
      DO i = 1, num_aspect_type
         ! calculate the cosine of solar illumination angle, cos(θ),
         ! ranging between −1 and 1, indicates if the sun is below or
         ! above the local horizon (note that values lower than 0 are set to 0 indicate self shadow)
         IF (i == 9) THEN
            ! flat area, no slope
            cosill_type(i) = 1.0_r8
         ELSE
            slp_type_c_rad = atan(slp_type_c(i))
            cosill_type(i) = cos(slp_type_c_rad)+tan(zen_rad)*sin(slp_type_c_rad)*cos(asp_type_c(i))
         ENDIF

         ! Ensure that the cosine of illumination angle is between 0 and 1
         IF (cosill_type(i)>1) cosill_type(i) = 1
         IF (cosill_type(i)<0) cosill_type(i) = 0

         ! downscaling beam radiation
         a_p = area_type_c(i)
         IF (a_p.gt.1.0) a_p = 1
         IF (a_p.lt.0) a_p = 0
         beam_swrad_type(i) = cosill_type(i)*opt_factor*a_p*beam_swrad_g
      ENDDO
      beam_swrad_c = sum(beam_swrad_type)

      ! do not downscale diffuse radiation
      diff_swrad_c = diff_swrad_g

      forc_swrad_c = beam_swrad_c+diff_swrad_c

      ! But ensure that we don't depart too far from the atmospheric forcing value:
      ! negative values of swrad are certainly bad, but small positive values might
      ! also be bad. We can especially run into trouble due to the normalization: a
      ! small swrad value in one column can lead to a big normalization factor,
      ! leading to huge swrad values in other columns.

      forc_swrad_c = min(forc_swrad_c, &
               forc_swrad_g * (1._r8 + shortwave_downscaling_limit)) 
      forc_swrad_c = max(forc_swrad_c, &
               forc_swrad_g * (1._r8 - shortwave_downscaling_limit))
      ! Ensure that the denominator is not 0 during shortwave normalization
      IF (forc_swrad_c < 1.e-4) forc_swrad_c = 0.0001

   END SUBROUTINE downscale_shortwave_simple



END MODULE MOD_ForcingDownscaling
