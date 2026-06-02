#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Operator_Const
!-----------------------------------------------------------------------
! DESCRIPTION:
! 1. Define constants (do not rely on satellite parameters) used in RTM.
! 2. Define parameters (IGBP) used in RTM.
!
! REFERENCES:
!   [1] L-band Microwave Emission of the Biosphere (L-MEB) Model:  Description 
!       and calibration against experimental data sets over crop fields.
!
!   [2] Wigneron, J. P., Jackson, T. J., O'neill, P., De Lannoy, G., de Rosnay, P., Walker,
!       J. P., ... & Kerr, Y. (2017). Modelling the passive microwave signature from land surfaces:
!       A review of recent results and application to the L-band SMOS & SMAP soil moisture retrieval algorithms.
!       Remote Sensing of Environment, 192, 238-262.
!
! AUTHOR:
!   Lu Li, 12/2024: Initial version
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Vars_Global, only: pi, N_land_classification

   IMPLICIT NONE
   SAVE

   ! Constant variables
   real(r8), parameter :: C = 2.998e8           ! speed of light (m/s)
   real(r8), parameter :: mu0 = 4.*pi*1e-7      ! vacuum permeability (H/m)
   real(r8), parameter :: eps_0 = 8.854e-12     ! dielectric constant of free space (Klein and Swift 1977) [Farads/meter]
   real(r8), parameter :: z0 = sqrt(mu0/eps_0)  ! impedance of free space (Ohm)
   real(r8), parameter :: eps_w_inf = 4.9       ! dielectric constant at infinite frequency (Stogryn 1971),
   real(r8), parameter :: rho_soil = 2.66       ! soil specific density (g/cm3)
   real(r8), parameter :: f0w = 9.              ! relaxation frequency of liquid water (GHz)
   real(r8), parameter :: rgh_surf = 2.2        ! soil surface roughness (cm)
   complex(r8), parameter :: jj = (0., 1.)      ! imaginary unit for complex number


#ifdef LULC_IGBP

   ! MODIS IGBP Land Use/Land Cover System Legend
   !---------------------------
   ! 0  Ocean                               !  海洋
   ! 1  Evergreen Needleleaf Forests        !  常绿针叶林
   ! 2  Evergreen Broadleaf Forests         !  常绿阔叶林
   ! 3  Deciduous Needleleaf Forests        !  落叶针叶林
   ! 4  Deciduous Broadleaf Forests         !  落叶阔叶林
   ! 5  Mixed Forests                       !  混交林
   ! 6  Closed Shrublands                   !  密闭灌丛
   ! 7  Open Shrublands                     !  稀疏灌丛
   ! 8  Woody Savannas                      !  木本稀树草原
   ! 9  Savannas                            !  稀树草原
   !10  Grasslands                          !  草地
   !11  Permanent Wetlands                  !  永久性湿地
   !12  Croplands                           !  农田
   !13  Urban and Built-up Lands            !  城市与建成区
   !14  Cropland/Natural Vegetation Mosaics !  农田-自然植被镶嵌区
   !15  Permanent Snow and Ice              !  永久冰雪
   !16  Barren                              !  裸地
   !17  Water Bodies                        !  水体

   ! empirical parameters to account for the dependence of optical depth on incidence angle [1]
   real(r8), parameter, dimension(N_land_classification) :: tth &
      = (/0.80, 1.00, 0.80, 0.49, 0.49, &
          1.00, 1.00, 1.00, 1.00, 1.00,  &
          1.00, 1.00, 1.00, 1.00, 1.00,  &
          1.00, 2.00/)
   real(r8), parameter, dimension(N_land_classification) :: ttv &
      = (/0.80, 1.00, 0.80, 0.46, 0.46, &
          1.00, 1.00, 1.00, 1.00, 1.00,  &
          1.00, 2.00, 1.00, 2.00, 1.00,  &
          1.00, 1.00/)

   ! empirical roughness parameters (Table 2 in [2])
   real(r8), parameter, dimension(N_land_classification) :: hr_SMAP &
      = (/0.160, 0.160, 0.160, 0.160, 0.160, &
          0.110, 0.110, 0.125, 0.156, 0.156, &
          0.100, 0.108, 0.000, 0.130, 0.000, &
          0.150, 0.000/)

   real(r8), parameter, dimension(N_land_classification) :: hr_SMOS &
      = (/0.300, 0.300, 0.300, 0.300, 0.300, &
          0.100, 0.100, 0.100, 0.100, 0.100, &
          0.100, 0.100, 0.100, 0.100, 0.000, &
          0.100, 0.000/)

   real(r8), parameter, dimension(N_land_classification) :: hr_P16  &
      = (/0.350, 0.460, 0.430, 0.450, 0.410, &
          0.260, 0.170, 0.350, 0.230, 0.130, &
          0.020, 0.170, 0.190, 0.220, 0.000, &
          0.020, 0.000/)

   ! b parameters for Wigneron vegetation model (ref?)
   real(r8), parameter, dimension(N_land_classification) :: b1 &
      = (/0.2600, 0.2260, 0.2600, 0.2260, 0.2260, &
          0.0375, 0.0375, 0.0375, 0.0375, 0.0375, &
          0.0000, 0.0500, 0.0000, 0.0500, 0.0000, &
          0.0000, 0.0500/)

   real(r8), parameter, dimension(N_land_classification) :: b2 &
      = (/0.0060, 0.0010, 0.0060, 0.0010, 0.0010, &
          0.0500, 0.0500, 0.0500, 0.0500, 0.0500, &
          0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
          0.0000, 0.0000/)

   real(r8), parameter, dimension(N_land_classification) :: b3 &
      = (/0.6900, 0.7000, 0.6900, 0.7000, 0.7000, &
          0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
          0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
          0.0000, 0.0000/)

   ! effective diffusion albedo (Table 3 in [2])
   real(r8), parameter, dimension(N_land_classification) :: w_SMAPL2 &
      = (/0.050, 0.050, 0.050, 0.050, 0.050, &
          0.050, 0.050, 0.050, 0.080, 0.050, &
          0.050, 0.000, 0.065, 0.000, 0.000, &
          0.000, 0.000/)

   real(r8), parameter, dimension(N_land_classification) :: w_CMEM &
      = (/0.080, 0.095, 0.080, 0.070, 0.070, &
          0.050, 0.050, 0.050, 0.050, 0.050, &
          0.000, 0.000, 0.000, 0.000, 0.000, &
          0.000, 0.000/)

   real(r8), parameter, dimension(N_land_classification) :: w_K16 &
      = (/0.050, 0.050, 0.060, 0.030, 0.050, &
          0.030, 0.050, 0.040, 0.020, 0.030, &
          0.000, 0.040, 0.000, 0.020, 0.000, &
          0.000, 0.000/)
   
   real(r8), parameter, dimension(N_land_classification) :: w_SMAPL4 &
      = (/0.120, 0.080, 0.120, 0.100, 0.120, &
          0.140, 0.110, 0.130, 0.120, 0.070, &
          0.000, 0.120, 0.000, 0.150, 0.000, &
          0.000, 0.000/)  
#endif

END MODULE MOD_DA_Operator_Const
!-----------------------------------------------------------------------------
#endif
