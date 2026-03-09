#include <define.h>

SUBROUTINE Aggregation_SoilParameters ( &
      gland, dir_rawdata, dir_model_landdata, lc_year)

!-----------------------------------------------------------------------
! !DESCRIPTION:
!  Create soil hydraulic and thermal parameters for the modeling reolustion
!
! !REFERENCES:
!  1) Dai, Y., Q. Xin, N. Wei, Y. Zhang, W. Shangguan, H. Yuan, S. Zhang, S.
!  Liu, X. Lu, 2019. A global high-resolution dataset of soil hydraulic and
!  thermal properties for land surface modeling. Journal of Advances in
!  Modeling Earth Systems,11, 2996-3023.
!
!  2) Dai, Y., N. Wei, H. Yuan, S. Zhang, W. Shangguan, S. Liu, and X. Lu, 2019.
!  Evaluation of soil thermal conductivity schemes for use in land surface
!  modelling, Journal of Advances in Modeling Earth Systems, 11, 3454-3473.
!
!  3) Dai, Y., W. Shangguan, Q. Duan, B. Liu, S. Fu, and G. Niu, 2013.
!  Development of a China dataset of soil hydraulic parameters using
!  pedotransfer functions for land surface modeling. Journal of
!  Hydrometeorology 14, 869-887
!
!  Original author: Yongjiu Dai and Wei Shangguan, 02/2014
!
! !REVISIONS:
!  06/2019, Nan Wei: add algorithms of fitting soil water retention curves to
!                    aggregate soil hydraulic parameters from pixels to a patch.
!
!  01/2022, Shupeng Zhang and Nan Wei: porting codes to MPI parallel version.
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_LandPatch
   USE MOD_Land2mWMO
   USE MOD_NetCDFBlock
   USE MOD_NetCDFVector
   USE MOD_AggregationRequestData
#ifdef RangeCheck
   USE MOD_RangeCheck
#endif
   USE MOD_Utils
   USE MOD_UserDefFun
#ifdef SrfdataDiag
   USE MOD_SrfdataDiag
#endif

   IMPLICIT NONE

   ! arguments:
   integer, intent(in) :: lc_year
   type(grid_type),  intent(in) :: gland
   character(len=*), intent(in) :: dir_rawdata
   character(len=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ---------------------------------------------------------------
   character(len=256) :: landdir, lndname, cyear
   character(len=256) :: c
   integer :: nsl, ipatch, L, np, LL, ipxstt, ipxend
   integer :: wmo_src

   type (block_data_real8_2d) :: vf_quartz_mineral_s_grid
   type (block_data_real8_2d) :: vf_gravels_s_grid
   type (block_data_real8_2d) :: vf_om_s_grid
   type (block_data_real8_2d) :: vf_sand_s_grid
   type (block_data_real8_2d) :: vf_clay_s_grid
   type (block_data_real8_2d) :: wf_gravels_s_grid
   type (block_data_real8_2d) :: wf_sand_s_grid
   type (block_data_real8_2d) :: wf_om_s_grid
   type (block_data_real8_2d) :: wf_clay_s_grid
   type (block_data_real8_2d) :: OM_density_s_grid
   type (block_data_real8_2d) :: BD_all_s_grid
   type (block_data_real8_2d) :: theta_s_grid
   type (block_data_real8_2d) :: psi_s_grid
   type (block_data_real8_2d) :: lambda_grid
#ifdef vanGenuchten_Mualem_SOIL_MODEL
   type (block_data_real8_2d) :: theta_r_grid
   type (block_data_real8_2d) :: alpha_vgm_grid
   type (block_data_real8_2d) :: L_vgm_grid
   type (block_data_real8_2d) :: n_vgm_grid
#endif
   type (block_data_real8_2d) :: k_s_grid
   type (block_data_real8_2d) :: csol_grid
   type (block_data_real8_2d) :: tksatu_grid
   type (block_data_real8_2d) :: tksatf_grid
   type (block_data_real8_2d) :: tkdry_grid
   type (block_data_real8_2d) :: k_solids_grid

   real(r8), allocatable :: vf_quartz_mineral_s_patches (:)
   real(r8), allocatable :: vf_gravels_s_patches (:)
   real(r8), allocatable :: vf_om_s_patches (:)
   real(r8), allocatable :: vf_sand_s_patches (:)
   real(r8), allocatable :: vf_clay_s_patches (:)
   real(r8), allocatable :: wf_gravels_s_patches (:)
   real(r8), allocatable :: wf_sand_s_patches (:)
   real(r8), allocatable :: wf_om_s_patches (:)
   real(r8), allocatable :: wf_clay_s_patches (:)
   real(r8), allocatable :: OM_density_s_patches (:)
   real(r8), allocatable :: BD_all_s_patches (:)
   real(r8), allocatable :: theta_s_patches (:)
   real(r8), allocatable :: psi_s_patches   (:)
   real(r8), allocatable :: lambda_patches  (:)
#ifdef vanGenuchten_Mualem_SOIL_MODEL
   real(r8), allocatable :: theta_r_patches (:)
   real(r8), allocatable :: alpha_vgm_patches  (:)
   real(r8), allocatable :: L_vgm_patches  (:)
   real(r8), allocatable :: n_vgm_patches  (:)
#endif
   real(r8), allocatable :: k_s_patches     (:)
   real(r8), allocatable :: csol_patches    (:)
   real(r8), allocatable :: tksatu_patches  (:)
   real(r8), allocatable :: tksatf_patches  (:)
   real(r8), allocatable :: tkdry_patches   (:)
   real(r8), allocatable :: k_solids_patches  (:)
   real(r8), allocatable :: BA_alpha_patches  (:)
   real(r8), allocatable :: BA_beta_patches  (:)

   real(r8), allocatable :: vf_quartz_mineral_s_one (:)
   real(r8), allocatable :: vf_gravels_s_one (:)
   real(r8), allocatable :: vf_om_s_one (:)
   real(r8), allocatable :: vf_sand_s_one (:)
   real(r8), allocatable :: vf_clay_s_one (:)
   real(r8), allocatable :: wf_gravels_s_one (:)
   real(r8), allocatable :: wf_sand_s_one (:)
   real(r8), allocatable :: wf_om_s_one (:)
   real(r8), allocatable :: wf_clay_s_one (:)
   real(r8), allocatable :: OM_density_s_one (:)
   real(r8), allocatable :: BD_all_s_one (:)
   real(r8), allocatable :: theta_s_one (:)
   real(r8), allocatable :: psi_s_one   (:)
   real(r8), allocatable :: lambda_one  (:)
#ifdef vanGenuchten_Mualem_SOIL_MODEL
   real(r8), allocatable :: theta_r_one  (:)
   real(r8), allocatable :: alpha_vgm_one  (:)
   real(r8), allocatable :: L_vgm_one  (:)
   real(r8), allocatable :: n_vgm_one  (:)
#endif
   real(r8), allocatable :: k_s_one     (:)
   real(r8), allocatable :: csol_one    (:)
   real(r8), allocatable :: tksatu_one  (:)
   real(r8), allocatable :: tksatf_one  (:)
   real(r8), allocatable :: tkdry_one   (:)
   real(r8), allocatable :: k_solids_one  (:)
   real(r8), allocatable :: area_one   (:)
   real(r8), allocatable :: BA_alpha_one  (:)
   real(r8), allocatable :: BA_beta_one  (:)

! local variables for estimating the upscaled soil parameters using the Levenberg-Marquardt fitting method
! ---------------------------------------------------------------
   integer, parameter   :: npointw  = 24
   integer, parameter   :: npointb  = 20
   real(r8),parameter   :: xdat(npointw) = (/1.,5.,10.,20.,30.,40.,50.,60.,70.,90.,110.,130.,150.,&
                                170.,210.,300.,345.,690.,1020.,5100.,15300.,20000.,100000.,1000000./)
                           !  points of soil pressure heads used for fitting SW retention curves

   real(r8),parameter   :: xdatsr(npointb)=(/0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,&
                                             0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0/)
                           !  points of soil saturation levels (Sr) used for fitting Ke-Sr relationship

   real(r8),allocatable :: ydatc  (:,:)   ! the Campbell SW retentions at fine grids
   real(r8),allocatable :: ydatcks(:,:)   ! the Campbell SW hydraulic conductivity at fine grids
   real(r8),allocatable :: ydatv  (:,:)   ! the van Genuchten & Mualem SW retentions at fine grids
   real(r8),allocatable :: ydatvks(:,:)   ! the van Genuchten & Mualem SW hydraulic conductivity at fine grids
   real(r8),allocatable :: THETA  (:)     ! the van Genuchten & Mualem relative SW at fine grids
   real(r8),allocatable :: ydatb(:,:)     ! the Balland and Arp (2005) Ke-Sr relationship at fine grids

   integer, parameter   :: nc = 3         ! the number of fitted parameters in Campbell SW retention curve (psi and lambda)
   integer, parameter   :: nv = 4         ! the number of fitted parameters in van Genuchten SW retention curve
                                          ! (theta_r, alpha and n)
   integer, parameter   :: nb = 2         ! the number of fitted parameters in Ke-Sr relationship (alpha and beta)

! Variables needed for Levenberg-Marquardt algorithm in MINPACK library
   real(r8),parameter   :: factor = 0.1
   real(r8),parameter   :: ftol = 1.0e-5
   real(r8),parameter   :: xtol = 1.0e-4
   real(r8),parameter   :: gtol = 0.0
   integer, parameter   :: mode = 1
   integer, parameter   :: nprint = 0
   integer              :: ldfjac,info,ipvtc(nc),ipvtv(nv),ipvtb(nb),maxfev,nfev,njev
   real(r8)             :: xc(nc),xv(nv),xb(nb),diagc(nc),diagv(nv),diagb(nb),qtfc(nc),qtfv(nv),qtfb(nb)
   real(r8),allocatable :: fjacc(:,:),fvecc(:),fjacv(:,:),fvecv(:),fjacb(:,:),fvecb(:)
   integer isiter                         ! flags to tell whether the iteration is completed, 1=Yes, 0=No

   ! Parameters to fill water body patches
   real(r8), parameter :: vf_quartz_mineral_fill_water(8) = 0.1
   real(r8), parameter :: vf_gravels_fill_water(8)        = 0.0
   real(r8), parameter :: vf_sand_fill_water(8)           = 0.09
   real(r8), parameter :: vf_clay_fill_water(8)           = 0.189
   real(r8), parameter :: vf_om_fill_water(8)             = 0.102
   real(r8), parameter :: wf_gravels_fill_water(8)        = 0.0
   real(r8), parameter :: wf_sand_fill_water(8)           = 0.1
   real(r8), parameter :: wf_om_fill_water(8)             = 0.1
   real(r8), parameter :: wf_clay_fill_water(8)           = 0.2
   real(r8), parameter :: theta_r_fill_water(8)           = 0.116
   real(r8), parameter :: alpha_vgm_fill_water(8)         = 0.01
   real(r8), parameter :: n_vgm_fill_water(8)             = 1.352
   real(r8), parameter :: theta_s_fill_water(8)           = 0.532
   real(r8), parameter :: k_s_fill_water(8)               = 11.616
   real(r8), parameter :: L_vgm_fill_water(8)             = 0.5
   real(r8), parameter :: psi_s_fill_water(8)             = -35.446
   real(r8), parameter :: lambda_fill_water(8)            = 0.108
   real(r8), parameter :: csol_fill_water(8)              = 1.102e6
   real(r8), parameter :: tksatu_fill_water(8)            = 1.145
   real(r8), parameter :: tksatf_fill_water(8)            = 2.401
   real(r8), parameter :: tkdry_fill_water(8)             = 0.136
   real(r8), parameter :: k_solids_fill_water(8)          = 1.545
   real(r8), parameter :: OM_density_fill_water(8)        = 62.064
   real(r8), parameter :: BD_all_fill_water(8)            = 1200
   real(r8), parameter :: BA_alpha_fill_water(8)          = 0.2
   real(r8), parameter :: BA_beta_fill_water (8)          = 10


#ifdef SrfdataDiag
   integer :: typpatch(N_land_classification+1), ityp
#endif

   external SW_CB_dist                    ! the objective function to be fitted for Campbell SW retention curve
   external SW_VG_dist                    ! the objective function to be fitted for van Genuchten SW retention curve
!   external Ke_Sr_dist                    ! the objective function to be fitted for Balland and Arp (2005) Ke-Sr relationship

      write(cyear,'(i4.4)') lc_year
      landdir = trim(dir_model_landdata) // '/soil/' // trim(cyear)

! ---------------------------------------------------------------------------------------------
! ... [2] aggregate the soil parameters from the resolution of raw data to modelling resolution
! ---------------------------------------------------------------------------------------------
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,'(/, A29)') 'Aggregate Soil Parameters ...'
         CALL system('mkdir -p ' // trim(adjustl(landdir)))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_worker) THEN

         allocate ( vf_quartz_mineral_s_patches(numpatch) )
         allocate ( vf_gravels_s_patches       (numpatch) )
         allocate ( vf_om_s_patches            (numpatch) )
         allocate ( vf_sand_s_patches          (numpatch) )
         allocate ( vf_clay_s_patches          (numpatch) )
         allocate ( wf_gravels_s_patches       (numpatch) )
         allocate ( wf_sand_s_patches          (numpatch) )
         allocate ( wf_om_s_patches            (numpatch) )
         allocate ( wf_clay_s_patches          (numpatch) )
         allocate ( OM_density_s_patches       (numpatch) )
         allocate ( BD_all_s_patches           (numpatch) )
         allocate ( theta_s_patches            (numpatch) )
         allocate ( psi_s_patches              (numpatch) )
         allocate ( lambda_patches             (numpatch) )
#ifdef vanGenuchten_Mualem_SOIL_MODEL
         allocate ( theta_r_patches   (numpatch) )
         allocate ( alpha_vgm_patches (numpatch) )
         allocate ( L_vgm_patches     (numpatch) )
         allocate ( n_vgm_patches     (numpatch) )
#endif
         allocate ( k_s_patches       (numpatch) )
         allocate ( csol_patches      (numpatch) )
         allocate ( tksatu_patches    (numpatch) )
         allocate ( tksatf_patches    (numpatch) )
         allocate ( tkdry_patches     (numpatch) )
         allocate ( k_solids_patches  (numpatch) )
         allocate ( BA_alpha_patches  (numpatch) )
         allocate ( BA_beta_patches   (numpatch) )

      ENDIF

      DO nsl = 1, 8

         write(c,'(i1)') nsl

         ! (1) volumetric fraction of quartz within mineral soil
         IF (p_is_io) THEN

            CALL allocate_block_data (gland, vf_quartz_mineral_s_grid)
            lndname = trim(dir_rawdata)//'/soil/vf_quartz_mineral_s.nc'
            CALL ncio_read_block (lndname, 'vf_quartz_mineral_s_l'//trim(c), gland, vf_quartz_mineral_s_grid)
#ifdef USEMPI
            CALL aggregation_data_daemon (gland, data_r8_2d_in1 = vf_quartz_mineral_s_grid)
#endif
         ENDIF

         IF (p_is_worker) THEN

            DO ipatch = 1, numpatch
               L = landpatch%settyp(ipatch)

               IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
                  wmo_src = wmo_source (landpatch%ielm(ipatch))
                  vf_quartz_mineral_s_patches (ipatch) = vf_quartz_mineral_s_patches (wmo_src)

                  CYCLE
               ENDIF

               IF (L /= 0) THEN
                  CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, area = area_one, &
                     data_r8_2d_in1 = vf_quartz_mineral_s_grid, data_r8_2d_out1 = vf_quartz_mineral_s_one)
                  CALL fillnan (vf_quartz_mineral_s_one, L == WATERBODY, vf_quartz_mineral_fill_water(nsl))
                  CALL fillnan (vf_quartz_mineral_s_one, L == GLACIERS , vf_quartz_mineral_fill_water(nsl))
                  vf_quartz_mineral_s_patches (ipatch) = sum (vf_quartz_mineral_s_one * (area_one/sum(area_one)))
               ELSE
                  vf_quartz_mineral_s_patches (ipatch) = -1.0e36_r8
               ENDIF

               IF (isnan_ud(vf_quartz_mineral_s_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in vf_quartz_mineral_s_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

            ENDDO

#ifdef USEMPI
            CALL aggregation_worker_done ()
#endif
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
         CALL check_vector_data ('vf_quartz_mineral_s lev '//trim(c), vf_quartz_mineral_s_patches)
#endif

         lndname = trim(landdir)//'/vf_quartz_mineral_s_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'vf_quartz_mineral_s_l'//trim(c)//'_patches', 'patch',&
                                 landpatch, vf_quartz_mineral_s_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (vf_quartz_mineral_s_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'vf_quartz_mineral_s_l'//trim(c), compress = 1, write_mode = 'one', &
            create_mode = (nsl==1))
#endif

         ! (2) volumetric fraction of gravels
         ! (3) volumetric fraction of sand
         ! (4) volumetric fraction of organic matter
         ! with the parameter alpha and beta in the Balland V. and P. A. Arp (2005) model
         IF (p_is_io) THEN

            CALL allocate_block_data (gland, vf_gravels_s_grid)
            lndname = trim(dir_rawdata)//'/soil/vf_gravels_s.nc'
            CALL ncio_read_block (lndname, 'vf_gravels_s_l'//trim(c), gland, vf_gravels_s_grid)

            CALL allocate_block_data (gland, vf_sand_s_grid)
            lndname = trim(dir_rawdata)//'/soil/vf_sand_s.nc'
            CALL ncio_read_block (lndname, 'vf_sand_s_l'//trim(c), gland, vf_sand_s_grid)

            CALL allocate_block_data (gland, vf_om_s_grid)
            lndname = trim(dir_rawdata)//'/soil/vf_om_s.nc'
            CALL ncio_read_block (lndname, 'vf_om_s_l'//trim(c), gland, vf_om_s_grid)

#ifdef USEMPI
            CALL aggregation_data_daemon (gland, data_r8_2d_in1 = vf_gravels_s_grid, &
               data_r8_2d_in2 = vf_sand_s_grid,  data_r8_2d_in3 = vf_om_s_grid)
#endif
         ENDIF

         IF (p_is_worker) THEN

            DO ipatch = 1, numpatch
               L = landpatch%settyp(ipatch)

               IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
                  wmo_src = wmo_source (landpatch%ielm(ipatch))

                  vf_gravels_s_patches (ipatch) = vf_gravels_s_patches (wmo_src)
                  vf_sand_s_patches    (ipatch) = vf_sand_s_patches    (wmo_src)
                  vf_om_s_patches      (ipatch) = vf_om_s_patches      (wmo_src)

                  CYCLE
               ENDIF

               IF (L /= 0) THEN
                  CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, area = area_one, &
                     data_r8_2d_in1 = vf_gravels_s_grid, data_r8_2d_out1 = vf_gravels_s_one, &
                     data_r8_2d_in2 = vf_sand_s_grid,    data_r8_2d_out2 = vf_sand_s_one, &
                     data_r8_2d_in3 = vf_om_s_grid,      data_r8_2d_out3 = vf_om_s_one)

                  CALL fillnan (vf_gravels_s_one, L == WATERBODY, vf_gravels_fill_water(nsl))
                  CALL fillnan (vf_sand_s_one   , L == WATERBODY, vf_sand_fill_water(nsl)   )
                  CALL fillnan (vf_om_s_one     , L == WATERBODY, vf_om_fill_water(nsl)     )

                  CALL fillnan (vf_gravels_s_one, L == GLACIERS , vf_gravels_fill_water(nsl))
                  CALL fillnan (vf_sand_s_one   , L == GLACIERS , vf_sand_fill_water(nsl)   )
                  CALL fillnan (vf_om_s_one     , L == GLACIERS , vf_om_fill_water(nsl)     )

                  vf_gravels_s_patches (ipatch) = sum (vf_gravels_s_one * (area_one/sum(area_one)))
                  vf_sand_s_patches (ipatch) = sum (vf_sand_s_one * (area_one/sum(area_one)))
                  vf_om_s_patches (ipatch) = sum (vf_om_s_one * (area_one/sum(area_one)))

                  ! the parameter values of Balland and Arp (2005) Ke-Sr relationship,
                  ! modified by Barry-Macaulay et al.(2015), Evaluation of soil thermal conductivity models

                  allocate(BA_alpha_one  (size(area_one)))
                  allocate(BA_beta_one   (size(area_one)))
                  WHERE ((vf_gravels_s_one + vf_sand_s_one) > 0.4)
                     BA_alpha_one = 0.38
                     BA_beta_one = 35.0
                  ELSEWHERE ((vf_gravels_s_one + vf_sand_s_one) > 0.25)
                     BA_alpha_one = 0.24
                     BA_beta_one = 26.0
                  ELSEWHERE
                     BA_alpha_one = 0.2
                     BA_beta_one = 10.0
                  END WHERE

                  BA_alpha_patches (ipatch) = median (BA_alpha_one, size(BA_alpha_one), spval)
                  BA_beta_patches (ipatch) = median (BA_beta_one, size(BA_beta_one), spval)

                  deallocate(BA_alpha_one)
                  deallocate(BA_beta_one)

               ELSE
                  vf_gravels_s_patches (ipatch) = -1.0e36_r8
                  vf_sand_s_patches (ipatch) = -1.0e36_r8
                  vf_om_s_patches (ipatch) = -1.0e36_r8
                  BA_alpha_patches (ipatch) = -1.0e36_r8
                  BA_beta_patches (ipatch) = -1.0e36_r8
               ENDIF

               IF (isnan_ud(vf_gravels_s_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in vf_gravels_s_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

               IF (isnan_ud(vf_sand_s_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in vf_sand_s_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

               IF (isnan_ud(vf_om_s_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in vf_om_s_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

               IF (isnan_ud(BA_alpha_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in BA_alpha_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

               IF (isnan_ud(BA_beta_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in BA_beta_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

            ENDDO

#ifdef USEMPI
            CALL aggregation_worker_done ()
#endif
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
         CALL check_vector_data ('vf_gravels_s lev '//trim(c), vf_gravels_s_patches)
         CALL check_vector_data ('vf_sand_s lev '//trim(c), vf_sand_s_patches)
         CALL check_vector_data ('vf_om_s lev '//trim(c), vf_om_s_patches)
         CALL check_vector_data ('BA_alpha lev '//trim(c), BA_alpha_patches)
         CALL check_vector_data ('BA_beta lev '//trim(c), BA_beta_patches)
#endif

         lndname = trim(landdir)//'/vf_gravels_s_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'vf_gravels_s_l'//trim(c)//'_patches', 'patch',&
                                 landpatch, vf_gravels_s_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (vf_gravels_s_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'vf_gravels_s_l'//trim(c), compress = 1, write_mode = 'one')
#endif

         lndname = trim(landdir)//'/vf_sand_s_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'vf_sand_s_l'//trim(c)//'_patches', 'patch',&
                                 landpatch, vf_sand_s_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (vf_sand_s_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'vf_sand_s_l'//trim(c), compress = 1, write_mode = 'one')
#endif

         lndname = trim(landdir)//'/vf_om_s_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'vf_om_s_l'//trim(c)//'_patches', 'patch',&
                                 landpatch, vf_om_s_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (vf_om_s_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'vf_om_s_l'//trim(c), compress = 1, write_mode = 'one')
#endif

         lndname = trim(landdir)//'/BA_alpha_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'BA_alpha_l'//trim(c)//'_patches', 'patch',&
                                 landpatch, BA_alpha_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (BA_alpha_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'BA_alpha_l'//trim(c), compress = 1, write_mode = 'one')
#endif

         lndname = trim(landdir)//'/BA_beta_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'BA_beta_l'//trim(c)//'_patches', 'patch',&
                                 landpatch, BA_beta_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (BA_beta_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'BA_beta_l'//trim(c), compress = 1, write_mode = 'one')
#endif

         ! (5) gravimetric fraction of gravels
         IF (p_is_io) THEN

            CALL allocate_block_data (gland, wf_gravels_s_grid)
            lndname = trim(dir_rawdata)//'/soil/wf_gravels_s.nc'
            CALL ncio_read_block (lndname, 'wf_gravels_s_l'//trim(c), gland, wf_gravels_s_grid)
#ifdef USEMPI
            CALL aggregation_data_daemon (gland, data_r8_2d_in1 = wf_gravels_s_grid)
#endif
         ENDIF

         IF (p_is_worker) THEN

            DO ipatch = 1, numpatch
               L = landpatch%settyp(ipatch)

               IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
                  wmo_src = wmo_source (landpatch%ielm(ipatch))

                  wf_gravels_s_patches (ipatch) = wf_gravels_s_patches (wmo_src)

                  CYCLE
               ENDIF

               IF (L /= 0) THEN
                  CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, area = area_one, &
                     data_r8_2d_in1 = wf_gravels_s_grid, data_r8_2d_out1 = wf_gravels_s_one)
                  CALL fillnan (wf_gravels_s_one, L == WATERBODY, wf_gravels_fill_water(nsl))
                  CALL fillnan (wf_gravels_s_one, L == GLACIERS , wf_gravels_fill_water(nsl))
                  wf_gravels_s_patches (ipatch) = sum (wf_gravels_s_one * (area_one/sum(area_one)))
               ELSE
                  wf_gravels_s_patches (ipatch) = -1.0e36_r8
               ENDIF

               IF (isnan_ud(wf_gravels_s_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in wf_gravels_s_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

            ENDDO

#ifdef USEMPI
            CALL aggregation_worker_done ()
#endif
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
         CALL check_vector_data ('wf_gravels_s lev '//trim(c), wf_gravels_s_patches)
#endif

         lndname = trim(landdir)//'/wf_gravels_s_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'wf_gravels_s_l'//trim(c)//'_patches', 'patch',&
                                 landpatch, wf_gravels_s_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (wf_gravels_s_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'wf_gravels_s_l'//trim(c), compress = 1, write_mode = 'one')
#endif

         ! (6) gravimetric fraction of sand
         IF (p_is_io) THEN

            CALL allocate_block_data (gland, wf_sand_s_grid)
            lndname = trim(dir_rawdata)//'/soil/wf_sand_s.nc'
            CALL ncio_read_block (lndname, 'wf_sand_s_l'//trim(c), gland, wf_sand_s_grid)
#ifdef USEMPI
            CALL aggregation_data_daemon (gland, data_r8_2d_in1 = wf_sand_s_grid)
#endif
         ENDIF

         IF (p_is_worker) THEN

            DO ipatch = 1, numpatch
               L = landpatch%settyp(ipatch)

               IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
                  wmo_src = wmo_source (landpatch%ielm(ipatch))

                  wf_sand_s_patches (ipatch) = wf_sand_s_patches (wmo_src)

                  CYCLE
               ENDIF

               IF (L /= 0) THEN
                  CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, area = area_one, &
                     data_r8_2d_in1 = wf_sand_s_grid, data_r8_2d_out1 = wf_sand_s_one)
                  CALL fillnan (wf_sand_s_one, L == WATERBODY, wf_sand_fill_water(nsl))
                  CALL fillnan (wf_sand_s_one, L == GLACIERS , wf_sand_fill_water(nsl))
                  wf_sand_s_patches (ipatch) = sum (wf_sand_s_one * (area_one/sum(area_one)))
               ELSE
                  wf_sand_s_patches (ipatch) = -1.0e36_r8
               ENDIF

               IF (isnan_ud(wf_sand_s_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in wf_sand_s_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

            ENDDO

#ifdef USEMPI
            CALL aggregation_worker_done ()
#endif
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
         CALL check_vector_data ('wf_sand_s lev '//trim(c), wf_sand_s_patches)
#endif

         lndname = trim(landdir)//'/wf_sand_s_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'wf_sand_s_l'//trim(c)//'_patches', 'patch',&
                                 landpatch, wf_sand_s_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (wf_sand_s_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'wf_sand_s_l'//trim(c), compress = 1, write_mode = 'one')
#endif


#ifdef vanGenuchten_Mualem_SOIL_MODEL

         ! (7) VGM's pore-connectivity parameter (L)
         ! (8) VGM's residual water content (theta_r) [cm3/cm3]
         ! (9) VGM's parameter corresponding approximately to the inverse of the air-entry value (alpha)
         ! (10) VGM's shape parameter (n)
         ! (11) saturated water content [cm3/cm3]
         ! (12) saturated hydraulic conductivity [cm/day]

         IF (p_is_io) THEN

            CALL allocate_block_data (gland, L_vgm_grid)
            lndname = trim(dir_rawdata)//'/soil/VGM_L.nc'
            CALL ncio_read_block (lndname, 'VGM_L_l'//trim(c), gland, L_vgm_grid)

            CALL allocate_block_data (gland, theta_r_grid)
            lndname = trim(dir_rawdata)//'/soil/VGM_theta_r.nc'
            CALL ncio_read_block (lndname, 'VGM_theta_r_l'//trim(c), gland, theta_r_grid)

            CALL allocate_block_data (gland, alpha_vgm_grid)
            lndname = trim(dir_rawdata)//'/soil/VGM_alpha.nc'
            CALL ncio_read_block (lndname, 'VGM_alpha_l'//trim(c), gland, alpha_vgm_grid)

            CALL allocate_block_data (gland, n_vgm_grid)
            lndname = trim(dir_rawdata)//'/soil/VGM_n.nc'
            CALL ncio_read_block (lndname, 'VGM_n_l'//trim(c), gland, n_vgm_grid)

            CALL allocate_block_data (gland, theta_s_grid)
            lndname = trim(dir_rawdata)//'/soil/theta_s.nc'
            CALL ncio_read_block (lndname, 'theta_s_l'//trim(c), gland, theta_s_grid)

            CALL allocate_block_data (gland, k_s_grid)
            lndname = trim(dir_rawdata)//'/soil/k_s.nc'
            CALL ncio_read_block (lndname, 'k_s_l'//trim(c), gland, k_s_grid)

#ifdef USEMPI
            CALL aggregation_data_daemon (gland, &
               data_r8_2d_in1 = theta_r_grid, data_r8_2d_in2 = alpha_vgm_grid, &
               data_r8_2d_in3 = n_vgm_grid,   data_r8_2d_in4 = theta_s_grid,   &
               data_r8_2d_in5 = k_s_grid,     data_r8_2d_in6 = L_vgm_grid    )
#endif
         ENDIF

         IF (p_is_worker) THEN

            DO ipatch = 1, numpatch
               L = landpatch%settyp(ipatch)

               IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
                  wmo_src = wmo_source (landpatch%ielm(ipatch))

                  theta_r_patches   (ipatch) = theta_r_patches   (wmo_src)
                  alpha_vgm_patches (ipatch) = alpha_vgm_patches (wmo_src)
                  n_vgm_patches     (ipatch) = n_vgm_patches     (wmo_src)
                  theta_s_patches   (ipatch) = theta_s_patches   (wmo_src)
                  k_s_patches       (ipatch) = k_s_patches       (wmo_src)
                  L_vgm_patches     (ipatch) = L_vgm_patches     (wmo_src)

                  CYCLE
               ENDIF

               IF (L /= 0) THEN
                  CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, area = area_one, &
                         data_r8_2d_in1 = theta_r_grid,   data_r8_2d_out1 = theta_r_one, &
                         data_r8_2d_in2 = alpha_vgm_grid, data_r8_2d_out2 = alpha_vgm_one, &
                         data_r8_2d_in3 = n_vgm_grid,     data_r8_2d_out3 = n_vgm_one, &
                         data_r8_2d_in4 = theta_s_grid,   data_r8_2d_out4 = theta_s_one, &
                         data_r8_2d_in5 = k_s_grid,       data_r8_2d_out5 = k_s_one, &
                         data_r8_2d_in6 = L_vgm_grid,     data_r8_2d_out6 = L_vgm_one  )

                  CALL fillnan (theta_r_one  , L == WATERBODY, theta_r_fill_water(nsl)  )
                  CALL fillnan (alpha_vgm_one, L == WATERBODY, alpha_vgm_fill_water(nsl))
                  CALL fillnan (n_vgm_one    , L == WATERBODY, n_vgm_fill_water(nsl)    )
                  CALL fillnan (theta_s_one  , L == WATERBODY, theta_s_fill_water(nsl)  )
                  CALL fillnan (k_s_one      , L == WATERBODY, k_s_fill_water(nsl)      )
                  CALL fillnan (L_vgm_one    , L == WATERBODY, L_vgm_fill_water(nsl)    )

                  CALL fillnan (theta_r_one  , L == GLACIERS , theta_r_fill_water(nsl)  )
                  CALL fillnan (alpha_vgm_one, L == GLACIERS , alpha_vgm_fill_water(nsl))
                  CALL fillnan (n_vgm_one    , L == GLACIERS , n_vgm_fill_water(nsl)    )
                  CALL fillnan (theta_s_one  , L == GLACIERS , theta_s_fill_water(nsl)  )
                  CALL fillnan (k_s_one      , L == GLACIERS , k_s_fill_water(nsl)      )
                  CALL fillnan (L_vgm_one    , L == GLACIERS , L_vgm_fill_water(nsl)    )

                  theta_r_patches (ipatch)   = sum (theta_r_one * (area_one/sum(area_one)))
                  alpha_vgm_patches (ipatch) = median (alpha_vgm_one, size(alpha_vgm_one), spval)
                  n_vgm_patches (ipatch)     = median (n_vgm_one, size(n_vgm_one), spval)
                  theta_s_patches (ipatch)   = sum (theta_s_one * (area_one/sum(area_one)))
                  k_s_patches (ipatch)       = product(k_s_one**(area_one/sum(area_one)))
                  L_vgm_patches (ipatch)     = median (L_vgm_one, size(L_vgm_one), spval)

                  IF (DEF_USE_SOILPAR_UPS_FIT) THEN
                     np = size(theta_r_one)

                     IF( np > 1 ) THEN
                        allocate ( ydatv  (1:np,npointw) )
                        allocate ( ydatvks(1:np,npointw) )
                        allocate ( THETA  (     npointw) )
! the jacobian matrix required in Levenberg-Marquardt fitting method
                        allocate ( fjacv  (npointw,nv) )           ! calculated in SW_VG_dist
! the values of objective functions to be fitted
                        allocate ( fvecv  (npointw)    )           ! calculated in SW_VG_dist

! SW VG retentions and hydraulic conductivity at fine grids for each patch
                        DO LL = 1,np
                           THETA         = (1 + (alpha_vgm_one(LL)*xdat)**n_vgm_one(LL))**(1.0/n_vgm_one(LL)-1)
                           ydatv(LL,:)   = theta_r_one(LL)+(theta_s_one(LL) - theta_r_one(LL)) * THETA
                           ydatvks(LL,:) = k_s_one(LL) * THETA**L_vgm_one(LL) &
                                         * (1-(1-THETA**(n_vgm_one(LL)/(n_vgm_one(LL)-1)))**(1.0-1.0/n_vgm_one(LL)))**2
                        ENDDO

! Fitting the van Genuchten SW retention parameters
                        ldfjac = npointw
                        xv(1) = theta_r_patches (ipatch)
                        xv(2) = alpha_vgm_patches (ipatch)
                        xv(3) = n_vgm_patches (ipatch)
                        xv(4) = k_s_patches (ipatch)
                        maxfev = 100 * ( nv + 1 )
                        isiter = 1

                        CALL lmder ( SW_VG_dist, npointw, nv, xv, fvecv, fjacv, ldfjac, ftol, xtol, gtol, maxfev, &
                              diagv, mode, factor, nprint, info, nfev, njev, ipvtv, qtfv,&
                              xdat, npointw, ydatv, ydatvks, np, theta_s_patches(ipatch),&
                              k_s_patches (ipatch), isiter, L_vgm_patches (ipatch))

                        IF ( xv(1) >= 0.0 .and. xv(1) <= theta_s_patches(ipatch) .and. xv(2) >= 1.0e-5 .and. xv(2) <= 1.0 .and. &
                             xv(3) >= 1.1 .and. xv(3) <= 10.0 .and. xv(4) > 0.0 .and. xv(4) <= 1.0e7 .and. isiter == 1) THEN
                           theta_r_patches(ipatch)   = xv(1)
                           alpha_vgm_patches(ipatch) = xv(2)
                           n_vgm_patches(ipatch)     = xv(3)
                           k_s_patches (ipatch)      = xv(4)
                        ENDIF

                        deallocate(ydatv)
                        deallocate(ydatvks)
                        deallocate(THETA)
                        deallocate(fjacv)
                        deallocate(fvecv)

                     ENDIF
                  ENDIF

               ELSE
                  theta_r_patches (ipatch)   = -1.0e36_r8
                  alpha_vgm_patches (ipatch) = -1.0e36_r8
                  n_vgm_patches (ipatch)     = -1.0e36_r8
                  theta_s_patches (ipatch)   = -1.0e36_r8
                  k_s_patches (ipatch)       = -1.0e36_r8
                  L_vgm_patches (ipatch)     = -1.0e36_r8
               ENDIF

               IF (isnan_ud(theta_r_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in theta_r_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

               IF (isnan_ud(alpha_vgm_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in alpha_vgm_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

               IF (isnan_ud(n_vgm_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in n_vgm_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

               IF (isnan_ud(theta_s_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in theta_s_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

               IF (isnan_ud(k_s_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in k_s_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

               IF (isnan_ud(L_vgm_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in L_vgm_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

            ENDDO

#ifdef USEMPI
            CALL aggregation_worker_done ()
#endif
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
         CALL check_vector_data ('theta_r lev '//trim(c), theta_r_patches)
         CALL check_vector_data ('alpha VGM lev '//trim(c), alpha_vgm_patches)
         CALL check_vector_data ('n VGM lev '//trim(c), n_vgm_patches)
         CALL check_vector_data ('theta_s lev '//trim(c), theta_s_patches)
         CALL check_vector_data ('k_s lev '//trim(c), k_s_patches)
         CALL check_vector_data ('L VGM lev '//trim(c), L_vgm_patches)
#endif

         lndname = trim(landdir)//'/theta_r_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'theta_r_l'//trim(c)//'_patches', 'patch', landpatch, &
                                 theta_r_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (theta_r_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'theta_r_l'//trim(c), compress = 1, write_mode = 'one')
#endif

         lndname = trim(landdir)//'/alpha_vgm_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'alpha_vgm_l'//trim(c)//'_patches', 'patch', landpatch, &
                                 alpha_vgm_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (alpha_vgm_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'alpha_vgm_l'//trim(c), compress = 1, write_mode = 'one')
#endif

         lndname = trim(landdir)//'/n_vgm_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'n_vgm_l'//trim(c)//'_patches', 'patch', landpatch, &
                                 n_vgm_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (n_vgm_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'n_vgm_l'//trim(c), compress = 1, write_mode = 'one')
#endif

         lndname = trim(landdir)//'/theta_s_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'theta_s_l'//trim(c)//'_patches', 'patch', landpatch, &
                                 theta_s_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (theta_s_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'theta_s_l'//trim(c), compress = 1, write_mode = 'one')
#endif

         lndname = trim(landdir)//'/k_s_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'k_s_l'//trim(c)//'_patches', 'patch', landpatch, &
                                 k_s_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (k_s_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'k_s_l'//trim(c), compress = 1, write_mode = 'one')
#endif

         lndname = trim(landdir)//'/L_vgm_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'L_vgm_l'//trim(c)//'_patches', 'patch', landpatch, &
                                 L_vgm_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (L_vgm_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'L_vgm_l'//trim(c), compress = 1, write_mode = 'one')
#endif

#endif

         ! (11) saturated water content [cm3/cm3]
         ! (12) saturated hydraulic conductivity [cm/day]
         ! (13) matric potential at saturation (psi_s) [cm]
         ! (14) pore size distribution index [dimensionless]

         IF (p_is_io) THEN

            CALL allocate_block_data (gland, theta_s_grid)
            lndname = trim(dir_rawdata)//'/soil/theta_s.nc'
            CALL ncio_read_block (lndname, 'theta_s_l'//trim(c), gland, theta_s_grid)

            CALL allocate_block_data (gland, k_s_grid)
            lndname = trim(dir_rawdata)//'/soil/k_s.nc'
            CALL ncio_read_block (lndname, 'k_s_l'//trim(c), gland, k_s_grid)

            CALL allocate_block_data (gland, psi_s_grid)
            lndname = trim(dir_rawdata)//'/soil/psi_s.nc'
            CALL ncio_read_block (lndname, 'psi_s_l'//trim(c), gland, psi_s_grid)

            CALL allocate_block_data (gland, lambda_grid)
            lndname = trim(dir_rawdata)//'/soil/lambda.nc'
            CALL ncio_read_block (lndname, 'lambda_l'//trim(c), gland, lambda_grid)

#ifdef USEMPI
            CALL aggregation_data_daemon (gland, data_r8_2d_in1 = theta_s_grid, data_r8_2d_in2 = k_s_grid, &
                                                 data_r8_2d_in3 = psi_s_grid,   data_r8_2d_in4 = lambda_grid)
#endif
         ENDIF

         IF (p_is_worker) THEN

            DO ipatch = 1, numpatch
               L = landpatch%settyp(ipatch)

               IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
                  wmo_src = wmo_source (landpatch%ielm(ipatch))

                  theta_s_patches (ipatch) = theta_s_patches (wmo_src)
                  k_s_patches     (ipatch) = k_s_patches     (wmo_src)
                  psi_s_patches   (ipatch) = psi_s_patches   (wmo_src)
                  lambda_patches  (ipatch) = lambda_patches  (wmo_src)

                  CYCLE
               ENDIF

               IF (L /= 0) THEN
                  CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, area = area_one, &
                               data_r8_2d_in1 = theta_s_grid, data_r8_2d_out1 = theta_s_one, &
                               data_r8_2d_in2 = k_s_grid,     data_r8_2d_out2 = k_s_one, &
                               data_r8_2d_in3 = psi_s_grid,   data_r8_2d_out3 = psi_s_one, &
                               data_r8_2d_in4 = lambda_grid,  data_r8_2d_out4 = lambda_one)

                  CALL fillnan (theta_s_one, L == WATERBODY, theta_s_fill_water(nsl))
                  CALL fillnan (k_s_one    , L == WATERBODY, k_s_fill_water(nsl)    )
                  CALL fillnan (psi_s_one  , L == WATERBODY, psi_s_fill_water(nsl)  )
                  CALL fillnan (lambda_one , L == WATERBODY, lambda_fill_water(nsl) )

                  CALL fillnan (theta_s_one, L == GLACIERS , theta_s_fill_water(nsl))
                  CALL fillnan (k_s_one    , L == GLACIERS , k_s_fill_water(nsl))
                  CALL fillnan (psi_s_one  , L == GLACIERS , psi_s_fill_water(nsl)  )
                  CALL fillnan (lambda_one , L == GLACIERS , lambda_fill_water(nsl) )

                  theta_s_patches (ipatch) = sum (theta_s_one * (area_one/sum(area_one)))
                  k_s_patches (ipatch)     = product(k_s_one**(area_one/sum(area_one)))
                  psi_s_patches (ipatch)   = median (psi_s_one, size(psi_s_one), spval)
                  lambda_patches (ipatch)  = median (lambda_one, size(lambda_one), spval)

                  IF (DEF_USE_SOILPAR_UPS_FIT) THEN
                     np = size(psi_s_one)

                     IF( np > 1 ) THEN
                        allocate ( ydatc  (1:np,npointw-7) )
                        allocate ( ydatcks(1:np,npointw-7) )
! the jacobian matrix required in Levenberg-Marquardt fitting method
                        allocate ( fjacc  (npointw-7,nc) )           ! calculated in SW_CB_dist
! the values of objective functions to be fitted
                        allocate ( fvecc  (npointw-7)    )           ! calculated in SW_CB_dist

! SW CB retentions at fine grids for each patch
                        DO LL = 1,np
                           ydatc(LL,:)   = (-1.0*xdat(8:)/psi_s_one(LL))**(-1.0*lambda_one(LL)) * theta_s_one(LL)
                           ydatcks(LL,:) = (-1.0*xdat(8:)/psi_s_one(LL))**(-3.0*lambda_one(LL)-2) * k_s_one(LL)
                        ENDDO

! Fitting the Campbell SW retention parameters
                        ldfjac = npointw-7
                        xc(1) = psi_s_patches (ipatch)
                        xc(2) = lambda_patches (ipatch)
                        xc(3) = k_s_patches (ipatch)
                        maxfev = 100 * ( nc + 1 )
                        isiter = 1

                        CALL lmder ( SW_CB_dist, npointw-7, nc, xc, fvecc, fjacc, ldfjac, ftol, xtol, gtol, maxfev, &
                              diagc, mode, factor, nprint, info, nfev, njev, ipvtc, qtfc,&
                              xdat(8:), npointw-7, ydatc, ydatcks, np, theta_s_patches(ipatch),k_s_patches(ipatch), isiter)

                        IF( xc(1) >= -300. .and. xc(1) < 0.0 .and. xc(2) > 0.0 .and. xc(2) <= 1.0 .and. &
                            xc(3) > 0.0 .and. xc(3) <= 1.0e7 .and. isiter == 1) THEN
                           psi_s_patches (ipatch) = xc(1)
                           lambda_patches(ipatch) = xc(2)
                           k_s_patches (ipatch)   = xc(3)
                        ENDIF

                        deallocate(ydatc)
                        deallocate(ydatcks)
                        deallocate(fjacc)
                        deallocate(fvecc)

                     ENDIF
                  ENDIF

               ELSE
                  theta_s_patches (ipatch) = -1.0e36_r8
                  k_s_patches (ipatch)     = -1.0e36_r8
                  psi_s_patches (ipatch)   = -1.0e36_r8
                  lambda_patches (ipatch)  = -1.0e36_r8
               ENDIF

               IF (isnan_ud(theta_s_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in theta_s_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

               IF (isnan_ud(k_s_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in k_s_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

               IF (isnan_ud(psi_s_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in psi_s_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

               IF (isnan_ud(lambda_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in lambda_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

            ENDDO

#ifdef USEMPI
            CALL aggregation_worker_done ()
#endif
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
         CALL check_vector_data ('theta_s lev '//trim(c), theta_s_patches)
         CALL check_vector_data ('k_s lev '//trim(c), k_s_patches)
         CALL check_vector_data ('psi_s lev '//trim(c), psi_s_patches)
         CALL check_vector_data ('lambda lev '//trim(c), lambda_patches)
#endif

#ifndef vanGenuchten_Mualem_SOIL_MODEL
         lndname = trim(landdir)//'/theta_s_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'theta_s_l'//trim(c)//'_patches', 'patch', landpatch, &
                                 theta_s_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (theta_s_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'theta_s_l'//trim(c), compress = 1, write_mode = 'one')
#endif

         lndname = trim(landdir)//'/k_s_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'k_s_l'//trim(c)//'_patches', 'patch', landpatch, &
                                 k_s_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (k_s_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'k_s_l'//trim(c), compress = 1, write_mode = 'one')
#endif
#endif

         lndname = trim(landdir)//'/psi_s_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'psi_s_l'//trim(c)//'_patches', 'patch', landpatch, &
                                 psi_s_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (psi_s_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'psi_s_l'//trim(c), compress = 1, write_mode = 'one')
#endif

         lndname = trim(landdir)//'/lambda_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'lambda_l'//trim(c)//'_patches', 'patch', landpatch, &
                                 lambda_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (lambda_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'lambda_l'//trim(c), compress = 1, write_mode = 'one')
#endif

         ! (15) heat capacity of soil solids [J/(m3 K)]
         IF (p_is_io) THEN
            CALL allocate_block_data (gland, csol_grid)
            lndname = trim(dir_rawdata)//'/soil/csol.nc'
            CALL ncio_read_block (lndname, 'csol_l'//trim(c), gland, csol_grid)
#ifdef USEMPI
            CALL aggregation_data_daemon (gland, data_r8_2d_in1 = csol_grid)
#endif
         ENDIF

         IF (p_is_worker) THEN

            DO ipatch = 1, numpatch
               L = landpatch%settyp(ipatch)

               IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
                  wmo_src = wmo_source (landpatch%ielm(ipatch))

                  csol_patches(ipatch) = csol_patches(wmo_src)

                  CYCLE
               ENDIF

               IF (L /= 0) THEN
                  CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, area = area_one, &
                     data_r8_2d_in1 = csol_grid, data_r8_2d_out1 = csol_one)
                  CALL fillnan (csol_one, L == WATERBODY, csol_fill_water(nsl))
                  CALL fillnan (csol_one, L == GLACIERS , csol_fill_water(nsl))
                  csol_patches (ipatch) = sum(csol_one*(area_one/sum(area_one)))
               ELSE
                  csol_patches (ipatch) = -1.0e36_r8
               ENDIF

               IF (isnan_ud(csol_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in csol_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

            ENDDO

#ifdef USEMPI
            CALL aggregation_worker_done ()
#endif
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
         CALL check_vector_data ('csol lev '//trim(c), csol_patches)
#endif

         lndname = trim(landdir)//'/csol_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'csol_l'//trim(c)//'_patches', 'patch', landpatch, &
                                 csol_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (csol_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'csol_l'//trim(c), compress = 1, write_mode = 'one')
#endif

         ! (16) thermal conductivity of unfrozen saturated soil [W/m-K]
         IF (p_is_io) THEN
            CALL allocate_block_data (gland, tksatu_grid)
            lndname = trim(dir_rawdata)//'/soil/tksatu.nc'
            CALL ncio_read_block (lndname, 'tksatu_l'//trim(c), gland, tksatu_grid)
#ifdef USEMPI
            CALL aggregation_data_daemon (gland, data_r8_2d_in1 = tksatu_grid)
#endif
         ENDIF

         IF (p_is_worker) THEN

            DO ipatch = 1, numpatch
               L = landpatch%settyp(ipatch)

               IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
                  wmo_src = wmo_source (landpatch%ielm(ipatch))

                  tksatu_patches(ipatch) = tksatu_patches(wmo_src)

                  CYCLE
               ENDIF

               IF (L /= 0) THEN
                  CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, area = area_one, &
                     data_r8_2d_in1 = tksatu_grid, data_r8_2d_out1 = tksatu_one)
                  CALL fillnan (tksatu_one, L == WATERBODY, tksatu_fill_water(nsl))
                  CALL fillnan (tksatu_one, L == GLACIERS , tksatu_fill_water(nsl))
                  tksatu_patches (ipatch) = product(tksatu_one**(area_one/sum(area_one)))
               ELSE
                  tksatu_patches (ipatch) = -1.0e36_r8
               ENDIF

               IF (isnan_ud(tksatu_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in tksatu_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

            ENDDO

#ifdef USEMPI
            CALL aggregation_worker_done ()
#endif
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
         CALL check_vector_data ('tksatu lev '//trim(c), tksatu_patches)
#endif

         lndname = trim(landdir)//'/tksatu_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'tksatu_l'//trim(c)//'_patches', 'patch', landpatch, &
                                 tksatu_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (tksatu_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'tksatu_l'//trim(c), compress = 1, write_mode = 'one')
#endif

         ! (17) thermal conductivity of frozen saturated soil [W/m-K]
         IF (p_is_io) THEN
            CALL allocate_block_data (gland, tksatf_grid)
            lndname = trim(dir_rawdata)//'/soil/tksatf.nc'
            CALL ncio_read_block (lndname, 'tksatf_l'//trim(c), gland, tksatf_grid)
#ifdef USEMPI
            CALL aggregation_data_daemon (gland, data_r8_2d_in1 = tksatf_grid)
#endif
         ENDIF

         IF (p_is_worker) THEN

            DO ipatch = 1, numpatch
               L = landpatch%settyp(ipatch)

               IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
                  wmo_src = wmo_source (landpatch%ielm(ipatch))

                  tksatf_patches(ipatch) = tksatf_patches(wmo_src)

                  CYCLE
               ENDIF

               IF (L /= 0) THEN
                  CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, area = area_one, &
                     data_r8_2d_in1 = tksatf_grid, data_r8_2d_out1 = tksatf_one)
                  CALL fillnan (tksatf_one, L == WATERBODY, tksatf_fill_water(nsl))
                  CALL fillnan (tksatf_one, L == GLACIERS , tksatf_fill_water(nsl))
                  tksatf_patches (ipatch) = product(tksatf_one**(area_one/sum(area_one)))
               ELSE
                  tksatf_patches (ipatch) = -1.0e36_r8
               ENDIF

               IF (isnan_ud(tksatf_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in tksatf_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

            ENDDO

#ifdef USEMPI
            CALL aggregation_worker_done ()
#endif
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
         CALL check_vector_data ('tksatf lev '//trim(c), tksatf_patches)
#endif

         lndname = trim(landdir)//'/tksatf_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'tksatf_l'//trim(c)//'_patches', 'patch', landpatch, &
                                 tksatf_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (tksatf_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'tksatf_l'//trim(c), compress = 1, write_mode = 'one')
#endif

         ! (18) thermal conductivity for dry soil [W/(m-K)]
         IF (p_is_io) THEN
            CALL allocate_block_data (gland, tkdry_grid)
            lndname = trim(dir_rawdata)//'/soil/tkdry.nc'
            CALL ncio_read_block (lndname, 'tkdry_l'//trim(c), gland, tkdry_grid)
#ifdef USEMPI
            CALL aggregation_data_daemon (gland, data_r8_2d_in1 = tkdry_grid)
#endif
         ENDIF

         IF (p_is_worker) THEN

            DO ipatch = 1, numpatch
               L = landpatch%settyp(ipatch)

               IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
                  wmo_src = wmo_source (landpatch%ielm(ipatch))

                  tkdry_patches(ipatch) = tkdry_patches(wmo_src)

                  CYCLE
               ENDIF

               IF (L /= 0) THEN
                  CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, area = area_one, &
                     data_r8_2d_in1 = tkdry_grid, data_r8_2d_out1 = tkdry_one)
                  CALL fillnan (tkdry_one, L == WATERBODY, tkdry_fill_water(nsl))
                  CALL fillnan (tkdry_one, L == GLACIERS , tkdry_fill_water(nsl))
                  tkdry_patches (ipatch) = product(tkdry_one**(area_one/sum(area_one)))
               ELSE
                  tkdry_patches (ipatch) = -1.0e36_r8
               ENDIF

               IF (isnan_ud(tkdry_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in tkdry_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

            ENDDO

#ifdef USEMPI
            CALL aggregation_worker_done ()
#endif
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
         CALL check_vector_data ('tkdry lev '//trim(c), tkdry_patches)
#endif

         lndname = trim(landdir)//'/tkdry_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'tkdry_l'//trim(c)//'_patches', 'patch', landpatch, &
                                 tkdry_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (tkdry_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'tkdry_l'//trim(c), compress = 1, write_mode = 'one')
#endif

         ! (19) thermal conductivity of soil solids [W/m-K]
         IF (p_is_io) THEN
            CALL allocate_block_data (gland, k_solids_grid)
            lndname = trim(dir_rawdata)//'/soil/k_solids.nc'
            CALL ncio_read_block (lndname, 'k_solids_l'//trim(c), gland, k_solids_grid)
#ifdef USEMPI
            CALL aggregation_data_daemon (gland, data_r8_2d_in1 = k_solids_grid)
#endif
         ENDIF

         IF (p_is_worker) THEN

            DO ipatch = 1, numpatch
               L = landpatch%settyp(ipatch)

               IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
                  wmo_src = wmo_source (landpatch%ielm(ipatch))

                  k_solids_patches(ipatch) = k_solids_patches(wmo_src)

                  CYCLE
               ENDIF

               IF (L /= 0) THEN
                  CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, area = area_one, &
                     data_r8_2d_in1 = k_solids_grid, data_r8_2d_out1 = k_solids_one)
                  CALL fillnan (k_solids_one, L == WATERBODY, k_solids_fill_water(nsl))
                  CALL fillnan (k_solids_one, L == GLACIERS , k_solids_fill_water(nsl))
                  k_solids_patches (ipatch) = product(k_solids_one**(area_one/sum(area_one)))
               ELSE
                  k_solids_patches (ipatch) = -1.0e36_r8
               ENDIF

               IF (isnan_ud(k_solids_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in k_solids_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

            ENDDO

#ifdef USEMPI
            CALL aggregation_worker_done ()
#endif
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
         CALL check_vector_data ('k_solids lev '//trim(c), k_solids_patches)
#endif

         lndname = trim(landdir)//'/k_solids_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'k_solids_l'//trim(c)//'_patches', 'patch', landpatch, &
                                 k_solids_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (k_solids_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'k_solids_l'//trim(c), compress = 1, write_mode = 'one')
#endif

         ! (20) OM_density [kg/m3]
         IF (p_is_io) THEN

            CALL allocate_block_data (gland, OM_density_s_grid)
            lndname = trim(dir_rawdata)//'/soil/OM_density_s.nc'
            CALL ncio_read_block (lndname, 'OM_density_s_l'//trim(c), gland, OM_density_s_grid)
#ifdef USEMPI
            CALL aggregation_data_daemon (gland, data_r8_2d_in1 = OM_density_s_grid)
#endif
         ENDIF

         IF (p_is_worker) THEN

            DO ipatch = 1, numpatch
               L = landpatch%settyp(ipatch)

               IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
                  wmo_src = wmo_source (landpatch%ielm(ipatch))

                  OM_density_s_patches(ipatch) = OM_density_s_patches(wmo_src)

                  CYCLE
               ENDIF

               IF (L /= 0) THEN
                  CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, area = area_one, &
                     data_r8_2d_in1 = OM_density_s_grid, data_r8_2d_out1 = OM_density_s_one)
                  CALL fillnan (OM_density_s_one, L == WATERBODY, OM_density_fill_water(nsl))
                  CALL fillnan (OM_density_s_one, L == GLACIERS , OM_density_fill_water(nsl))
                  OM_density_s_patches (ipatch) = sum (OM_density_s_one * (area_one/sum(area_one)))
               ELSE
                  OM_density_s_patches (ipatch) = -1.0e36_r8
               ENDIF

               IF (isnan_ud(OM_density_s_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in OM_density_s_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

            ENDDO

#ifdef USEMPI
            CALL aggregation_worker_done ()
#endif
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
         CALL check_vector_data ('OM_density_s lev '//trim(c), OM_density_s_patches)
#endif

         lndname = trim(landdir)//'/OM_density_s_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'OM_density_s_l'//trim(c)//'_patches', 'patch',&
                                 landpatch, OM_density_s_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (OM_density_s_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'OM_density_s_l'//trim(c), compress = 1, write_mode = 'one')
#endif

         ! (21) bulk density of soil (GRAVELS + OM + Mineral Soils)
         IF (p_is_io) THEN

            CALL allocate_block_data (gland, BD_all_s_grid)
            lndname = trim(dir_rawdata)//'/soil/BD_all_s.nc'
            CALL ncio_read_block (lndname, 'BD_all_s_l'//trim(c), gland, BD_all_s_grid)
#ifdef USEMPI
            CALL aggregation_data_daemon (gland, data_r8_2d_in1 = BD_all_s_grid)
#endif
         ENDIF

         IF (p_is_worker) THEN

            DO ipatch = 1, numpatch
               L = landpatch%settyp(ipatch)

               IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
                  wmo_src = wmo_source (landpatch%ielm(ipatch))

                  BD_all_s_patches(ipatch) = BD_all_s_patches(wmo_src)

                  CYCLE
               ENDIF


               IF (L /= 0) THEN
                  CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, area = area_one, &
                     data_r8_2d_in1 = BD_all_s_grid, data_r8_2d_out1 = BD_all_s_one)
                  CALL fillnan (BD_all_s_one, L == WATERBODY, BD_all_fill_water(nsl))
                  CALL fillnan (BD_all_s_one, L == GLACIERS , BD_all_fill_water(nsl))
                  BD_all_s_patches (ipatch) = sum (BD_all_s_one * (area_one/sum(area_one)))
               ELSE
                  BD_all_s_patches (ipatch) = -1.0e36_r8
               ENDIF

               IF (isnan_ud(BD_all_s_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in BD_all_s_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

            ENDDO

#ifdef USEMPI
            CALL aggregation_worker_done ()
#endif
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
         CALL check_vector_data ('BD_all_s lev '//trim(c), BD_all_s_patches)
#endif

         lndname = trim(landdir)//'/BD_all_s_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'BD_all_s_l'//trim(c)//'_patches', 'patch',&
                                 landpatch, BD_all_s_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (BD_all_s_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'BD_all_s_l'//trim(c), compress = 1, write_mode = 'one')
#endif

         ! (22) volumetric fraction of clay
         IF (p_is_io) THEN

            CALL allocate_block_data (gland, vf_clay_s_grid)
            lndname = trim(dir_rawdata)//'/soil/vf_clay_s.nc'
            CALL ncio_read_block (lndname, 'vf_clay_s_l'//trim(c), gland, vf_clay_s_grid)
#ifdef USEMPI
            CALL aggregation_data_daemon (gland, data_r8_2d_in1 = vf_clay_s_grid)
#endif
         ENDIF

         IF (p_is_worker) THEN

            DO ipatch = 1, numpatch
               L = landpatch%settyp(ipatch)

               IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
                  wmo_src = wmo_source (landpatch%ielm(ipatch))

                  vf_clay_s_patches(ipatch) = vf_clay_s_patches(wmo_src)

                  CYCLE
               ENDIF

               IF (L /= 0) THEN
                  CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, area = area_one, &
                     data_r8_2d_in1 = vf_clay_s_grid, data_r8_2d_out1 = vf_clay_s_one)
                  CALL fillnan (vf_clay_s_one, L == WATERBODY, vf_clay_fill_water(nsl))
                  CALL fillnan (vf_clay_s_one, L == GLACIERS , vf_clay_fill_water(nsl))
                  vf_clay_s_patches (ipatch) = sum (vf_clay_s_one * (area_one/sum(area_one)))
               ELSE
                  vf_clay_s_patches (ipatch) = -1.0e36_r8
               ENDIF

               IF (isnan_ud(vf_clay_s_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in vf_clay_s_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

            ENDDO

#ifdef USEMPI
            CALL aggregation_worker_done ()
#endif
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
         CALL check_vector_data ('vf_clay_s lev '//trim(c), vf_clay_s_patches)
#endif

         lndname = trim(landdir)//'/vf_clay_s_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'vf_clay_s_l'//trim(c)//'_patches', 'patch',&
                                 landpatch, vf_clay_s_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (vf_clay_s_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'vf_clay_s_l'//trim(c), compress = 1, write_mode = 'one')
#endif


         ! (23) gravimetric fraction of om
         IF (p_is_io) THEN

            CALL allocate_block_data (gland, wf_om_s_grid)
            lndname = trim(dir_rawdata)//'/soil/wf_om_s.nc'
            CALL ncio_read_block (lndname, 'wf_om_s_l'//trim(c), gland, wf_om_s_grid)
#ifdef USEMPI
            CALL aggregation_data_daemon (gland, data_r8_2d_in1 = wf_om_s_grid)
#endif
         ENDIF

         IF (p_is_worker) THEN

            DO ipatch = 1, numpatch
               L = landpatch%settyp(ipatch)

               IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
                  wmo_src = wmo_source (landpatch%ielm(ipatch))

                  wf_om_s_patches(ipatch) = wf_om_s_patches(wmo_src)

                  CYCLE
               ENDIF

               IF (L /= 0) THEN
                  CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, area = area_one, &
                     data_r8_2d_in1 = wf_om_s_grid, data_r8_2d_out1 = wf_om_s_one)
                  CALL fillnan (wf_om_s_one, L == WATERBODY, wf_om_fill_water(nsl))
                  CALL fillnan (wf_om_s_one, L == GLACIERS , wf_om_fill_water(nsl))
                  wf_om_s_patches (ipatch) = sum (wf_om_s_one * (area_one/sum(area_one)))
               ELSE
                  wf_om_s_patches (ipatch) = -1.0e36_r8
               ENDIF

               IF (isnan_ud(wf_om_s_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in wf_om_s_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

            ENDDO

#ifdef USEMPI
            CALL aggregation_worker_done ()
#endif
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
         CALL check_vector_data ('wf_om_s lev '//trim(c), wf_om_s_patches)
#endif

         lndname = trim(landdir)//'/wf_om_s_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'wf_om_s_l'//trim(c)//'_patches', 'patch',&
                                 landpatch, wf_om_s_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (wf_om_s_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'wf_om_s_l'//trim(c), compress = 1, write_mode = 'one')
#endif


         ! (24) gravimetric fraction of clay
         IF (p_is_io) THEN

            CALL allocate_block_data (gland, wf_clay_s_grid)
            lndname = trim(dir_rawdata)//'/soil/wf_clay_s.nc'
            CALL ncio_read_block (lndname, 'wf_clay_s_l'//trim(c), gland, wf_clay_s_grid)
#ifdef USEMPI
            CALL aggregation_data_daemon (gland, data_r8_2d_in1 = wf_clay_s_grid)
#endif
         ENDIF

         IF (p_is_worker) THEN

            DO ipatch = 1, numpatch
               L = landpatch%settyp(ipatch)

               IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
                  wmo_src = wmo_source (landpatch%ielm(ipatch))

                  wf_clay_s_patches(ipatch) = wf_clay_s_patches(wmo_src)

                  CYCLE
               ENDIF

               IF (L /= 0) THEN
                  CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, area = area_one, &
                     data_r8_2d_in1 = wf_clay_s_grid, data_r8_2d_out1 = wf_clay_s_one)
                  CALL fillnan (wf_clay_s_one, L == WATERBODY, wf_clay_fill_water(nsl))
                  CALL fillnan (wf_clay_s_one, L == GLACIERS , wf_clay_fill_water(nsl))
                  wf_clay_s_patches (ipatch) = sum (wf_clay_s_one * (area_one/sum(area_one)))
               ELSE
                  wf_clay_s_patches (ipatch) = -1.0e36_r8
               ENDIF

               IF (isnan_ud(wf_clay_s_patches(ipatch))) THEN
                  write(*,*) "Warning: NAN appears in wf_clay_s_patches."
                  write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
               ENDIF

            ENDDO

#ifdef USEMPI
            CALL aggregation_worker_done ()
#endif
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
         CALL check_vector_data ('wf_clay_s lev '//trim(c), wf_clay_s_patches)
#endif

         lndname = trim(landdir)//'/wf_clay_s_l'//trim(c)//'_patches.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'wf_clay_s_l'//trim(c)//'_patches', 'patch',&
                                 landpatch, wf_clay_s_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/soil_parameters_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (wf_clay_s_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'wf_clay_s_l'//trim(c), compress = 1, write_mode = 'one')
#endif
      ENDDO


      ! Deallocate the allocatable array
      ! --------------------------------

      IF (p_is_worker) THEN

         deallocate ( vf_quartz_mineral_s_patches )
         deallocate ( vf_gravels_s_patches )
         deallocate ( vf_om_s_patches )
         deallocate ( vf_sand_s_patches )
         deallocate ( vf_clay_s_patches )
         deallocate ( wf_gravels_s_patches )
         deallocate ( wf_sand_s_patches )
         deallocate ( wf_om_s_patches )
         deallocate ( wf_clay_s_patches )
         deallocate ( OM_density_s_patches )
         deallocate ( BD_all_s_patches )
         deallocate ( theta_s_patches )
         deallocate ( psi_s_patches   )
         deallocate ( lambda_patches  )
#ifdef vanGenuchten_Mualem_SOIL_MODEL
         deallocate ( theta_r_patches  )
         deallocate ( alpha_vgm_patches)
         deallocate ( L_vgm_patches    )
         deallocate ( n_vgm_patches    )
#endif
         deallocate ( k_s_patches     )
         deallocate ( csol_patches    )
         deallocate ( tksatu_patches  )
         deallocate ( tkdry_patches   )
         deallocate ( tksatf_patches  )
         deallocate ( k_solids_patches)
         deallocate ( BA_alpha_patches)
         deallocate ( BA_beta_patches )

         IF (allocated(vf_quartz_mineral_s_one)) deallocate (vf_quartz_mineral_s_one)
         IF (allocated(vf_gravels_s_one))        deallocate (vf_gravels_s_one)
         IF (allocated(vf_om_s_one))             deallocate (vf_om_s_one)
         IF (allocated(vf_sand_s_one))           deallocate (vf_sand_s_one)
         IF (allocated(vf_clay_s_one))           deallocate (vf_clay_s_one)
         IF (allocated(wf_gravels_s_one))        deallocate (wf_gravels_s_one)
         IF (allocated(wf_sand_s_one))           deallocate (wf_sand_s_one)
         IF (allocated(wf_om_s_one))             deallocate (wf_om_s_one)
         IF (allocated(wf_clay_s_one))           deallocate (wf_clay_s_one)
         IF (allocated(OM_density_s_one))        deallocate (OM_density_s_one)
         IF (allocated(BD_all_s_one))            deallocate (BD_all_s_one)
         IF (allocated(theta_s_one))             deallocate (theta_s_one)
         IF (allocated(psi_s_one  ))             deallocate (psi_s_one  )
         IF (allocated(lambda_one ))     deallocate (lambda_one )
#ifdef vanGenuchten_Mualem_SOIL_MODEL
         IF (allocated ( theta_r_one  )) deallocate ( theta_r_one  )
         IF (allocated ( alpha_vgm_one)) deallocate ( alpha_vgm_one)
         IF (allocated ( L_vgm_one    )) deallocate ( L_vgm_one    )
         IF (allocated ( n_vgm_one    )) deallocate ( n_vgm_one    )
#endif
         IF (allocated(k_s_one    ))     deallocate (k_s_one    )
         IF (allocated(csol_one   ))     deallocate (csol_one   )
         IF (allocated(tksatu_one ))     deallocate (tksatu_one )
         IF (allocated(tkdry_one  ))     deallocate (tkdry_one  )
         IF (allocated(tksatf_one ))     deallocate (tksatf_one)
         IF (allocated(k_solids_one))    deallocate (k_solids_one)
         IF (allocated(area_one))        deallocate (area_one)

      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

END SUBROUTINE Aggregation_SoilParameters


SUBROUTINE SW_CB_dist ( m, n, x, fvec, fjac, ldfjac, iflag, xdat, npoint, ydatc, ydatcks, nptf, phi, k_s, isiter)

!=================================================================
! DESCRIPTION:
! This is the subroutine for calculating the function/jacobian matrix
! of the distance between the fitted and prescribed SW retention and
! hydraulic conductivity curves for the Campbell model.
!
! Created by Nan Wei, 01/2019
! ----------------------------------------------------------------

   USE MOD_Precision
   IMPLICIT NONE

   integer m,n,ldfjac,iflag,i,nptf,isiter,npoint
   real(r8) x(n),fjac(ldfjac,n),fvec(m),xdat(npoint),ydatc(nptf,npoint),ydatcks(nptf,npoint),phi,k_s

      IF ( iflag == 0 ) THEN

         print*,x

      ELSEIF ( iflag == 1 ) THEN

         IF (x(1) >= 0.0 .or. abs(x(2)) >= 100. .or. x(3) <= 0.0) THEN
            isiter = 0
            RETURN
         ENDIF

         DO i = 1, m
            fvec(i) = sum((((-1.0*xdat(i)/x(1))**(-1.0*x(2)) * phi - ydatc(:,i))/phi)**2) &
                    + sum(((log10(-1.0*xdat(i)/x(1))*(-3.0*x(2)-2) + log10(x(3)) - log10(ydatcks(:,i)))/log10(k_s))**2)
         ENDDO

      ELSEIF ( iflag == 2 ) THEN

         IF (x(1) >= 0.0 .or. abs(x(2)) >= 100. .or. x(3) <= 0.0) THEN
            isiter = 0
            RETURN
         ENDIF

         DO i = 1, m
            fjac(i,1) = sum(2.0*(((-1.0*xdat(i)/x(1))**(-1.0*x(2)) * phi - ydatc(:,i))/phi)*&
                        x(2) * (-1.0*xdat(i)/x(1))**(-1.0*x(2)) / x(1)) &
                      + sum(2.0*((log10(-1.0*xdat(i)/x(1))*(-3.0*x(2)-2) + log10(x(3)) - log10(ydatcks(:,i)))/log10(k_s))*&
                        (3.0*x(2)+2)/(x(1)*log(10.))/log10(k_s))
            fjac(i,2) = sum(-2.0*(((-1.0*xdat(i)/x(1))**(-1.0*x(2)) * phi - ydatc(:,i))/phi)*&
                        (-1.0*xdat(i)/x(1))**(-1.0*x(2)) * log(-1.0*xdat(i)/x(1))) &
                      + sum(-6.0*((log10(-1.0*xdat(i)/x(1))*(-3.0*x(2)-2) + log10(x(3)) - log10(ydatcks(:,i)))/log10(k_s))*&
                        log10(-1.0*xdat(i)/x(1))/log10(k_s))
            fjac(i,3) = sum(2.0*((log10(-1.0*xdat(i)/x(1))*(-3.0*x(2)-2) + log10(x(3)) - log10(ydatcks(:,i)))/log10(k_s))/&
                        (x(3)*log(10.)*log10(k_s)))
         ENDDO

      ENDIF

END SUBROUTINE SW_CB_dist

SUBROUTINE SW_VG_dist ( m, n, x, fvec, fjac, ldfjac, iflag, xdat, npoint, ydatv, ydatvks, nptf, phi, k_s, isiter, L_vgm)

!=================================================================
! DESCRIPTION:
! This is the subroutine for calculating the function/jacobian matrix
! of the distance between the fitted and prescribed SW retention and
! hydraulic conductivity curves for the van Genuchten model.
!
! Created by Nan Wei, 01/2019
! ----------------------------------------------------------------

   USE MOD_Precision
   IMPLICIT NONE

   integer m,n,ldfjac,iflag,i,nptf,isiter,npoint
   real(r8) x(n),fjac(ldfjac,n),fvec(m),xdat(npoint),ydatv(nptf,npoint),ydatvks(nptf,npoint),phi,k_s,L_vgm

      IF ( iflag == 0 ) THEN

         print*,x

      ELSEIF ( iflag == 1 ) THEN

         IF (x(2) <= 0.0 .or. x(3) <= 0.1 .or. x(3) >= 100. .or. x(4) <= 0.0) THEN
            isiter = 0
            RETURN
         ENDIF

         DO i = 1, m
            fvec(i) = sum(((x(1) + (phi - x(1))*(1+(x(2)*xdat(i))**x(3))**(1.0/x(3)-1) - ydatv(:,i))/phi)**2) &
                    + sum(((log10(x(4)) + (1.0/x(3)-1)*L_vgm*log10(1+(x(2)*xdat(i))**x(3)) + &
                      log10((1.0-(1.0-1.0/(1+(x(2)*xdat(i))**x(3)))**(1-1.0/x(3)))**2) - log10(ydatvks(:,i)))/log10(k_s))**2)
         ENDDO

      ELSEIF ( iflag == 2 ) THEN

         IF (x(2) <= 0.0 .or. x(3) <= 0.1 .or. x(3) >= 100. .or. x(4) <= 0.0) THEN
            isiter = 0
            RETURN
         ENDIF

         DO i = 1, m
            fjac(i,1) = sum(2*((x(1) + (phi - x(1))*(1+(x(2)*xdat(i))**x(3))**(1.0/x(3)-1) - ydatv(:,i))/phi)*&
                        (1 - (1+(x(2)*xdat(i))**x(3))**(1.0/x(3)-1))/phi)
            fjac(i,2) = sum(2*((x(1) + (phi - x(1))*(1+(x(2)*xdat(i))**x(3))**(1.0/x(3)-1) - ydatv(:,i))/phi)/phi*&
                        (phi - x(1)) * (1 - x(3)) * (1+(x(2)*xdat(i))**x(3))**(1.0/x(3)-2) * x(2)**(x(3)-1) * xdat(i)**x(3)) &
                      + sum(2*((log10(x(4)) + (1.0/x(3)-1)*L_vgm*log10(1+(x(2)*xdat(i))**x(3)) + &
                        log10((1.0-(1.0-1.0/(1+(x(2)*xdat(i))**x(3)))**(1-1.0/x(3)))**2) - log10(ydatvks(:,i)))/log10(k_s)) * &
                        (L_vgm*(1.0-x(3))*x(2)**(x(3)-1)*xdat(i)**x(3)/((1+(x(2)*xdat(i))**x(3))*log(10.)) + &
                        2.0*(1.0-x(3))*(1.0-1.0/(1+(x(2)*xdat(i))**x(3)))**(-1.0/x(3))* &
                        x(2)**(x(3)-1)*xdat(i)**x(3)*(1+(x(2)*xdat(i))**x(3))**(-2)/ &
                        ((1.0-(1.0-1.0/(1+(x(2)*xdat(i))**x(3)))**(1-1.0/x(3)))*log(10.)))/log10(k_s))
            fjac(i,3) = sum(2*((x(1) + (phi - x(1))*(1+(x(2)*xdat(i))**x(3))**(1.0/x(3)-1) - ydatv(:,i))/phi)/phi*&
                        (phi - x(1)) * (1+(x(2)*xdat(i))**x(3))**(1.0/x(3)-1) *&
                        ((1.0-x(3))*(x(2)*xdat(i))**x(3)*log(x(2)*xdat(i))/(x(3)*(1+(x(2)*xdat(i))**x(3))) &
                        - log(1+(x(2)*xdat(i))**x(3))/x(3)**2)) &
                      + sum(2*((log10(x(4)) + (1.0/x(3)-1)*L_vgm*log10(1+(x(2)*xdat(i))**x(3)) + &
                        log10((1.0-(1.0-1.0/(1+(x(2)*xdat(i))**x(3)))**(1-1.0/x(3)))**2) - log10(ydatvks(:,i)))/log10(k_s)) * &
                        (-1.0*L_vgm*log10(1+(x(2)*xdat(i))**x(3))/x(3)**2 + &
                        (1.0/x(3)-1)*L_vgm*(x(2)*xdat(i))**x(3)*log10(x(2)*xdat(i))/(1+(x(2)*xdat(i))**x(3)) - &
                        2.0*(1.0-1.0/(1+(x(2)*xdat(i))**x(3)))**(1-1.0/x(3))/ &
                        (1.0-(1.0-1.0/(1+(x(2)*xdat(i))**x(3)))**(1-1.0/x(3))) * &
                        (log10(1.0-1.0/(1+(x(2)*xdat(i))**x(3)))/x(3)**2 + &
                        (1-1.0/x(3))*log10(x(2)*xdat(i))/(1+(x(2)*xdat(i))**x(3))))/log10(k_s))
            fjac(i,4) = sum(2*((log10(x(4)) + (1.0/x(3)-1)*L_vgm*log10(1+(x(2)*xdat(i))**x(3)) + &
                        log10((1.0-(1.0-1.0/(1+(x(2)*xdat(i))**x(3)))**(1-1.0/x(3)))**2) - log10(ydatvks(:,i)))/log10(k_s)) / &
                        (x(4)*log(10.))/log10(k_s))
         ENDDO

      ENDIF

END SUBROUTINE SW_VG_dist

!-----------------------------------------------------------------------
!EOP
