#include <define.h>

#ifdef SinglePoint
MODULE MOD_SingleSrfdata
!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!    This module includes subroutines to read or write surface data for
!    "SinglePoint".
!
!  Created by Shupeng Zhang, May 2023
!-----------------------------------------------------------------------

   USE MOD_Precision, only: r8
   USE MOD_Vars_Global
   USE MOD_Const_LC
   USE MOD_Namelist
   USE MOD_SPMD_Task
   IMPLICIT NONE
   SAVE

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   integer,  allocatable :: SITE_pfttyp  (:)
   real(r8), allocatable :: SITE_pctpfts (:)
#endif

#ifdef CROP
   integer,  allocatable :: SITE_croptyp (:)
   real(r8), allocatable :: SITE_pctcrop (:)
#endif

   real(r8) :: SITE_htop
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   real(r8), allocatable :: SITE_htop_pfts (:)
#endif

   real(r8), allocatable :: SITE_LAI_monthly (:,:)
   real(r8), allocatable :: SITE_SAI_monthly (:,:)
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   real(r8), allocatable :: SITE_LAI_pfts_monthly (:,:,:)
   real(r8), allocatable :: SITE_SAI_pfts_monthly (:,:,:)
#endif

   integer,  allocatable :: SITE_LAI_year (:)
   real(r8), allocatable :: SITE_LAI_8day (:,:)

   real(r8) :: SITE_lakedepth = 1.

   real(r8) :: SITE_soil_s_v_alb
   real(r8) :: SITE_soil_d_v_alb
   real(r8) :: SITE_soil_s_n_alb
   real(r8) :: SITE_soil_d_n_alb

   real(r8), allocatable :: SITE_soil_vf_quartz_mineral (:)
   real(r8), allocatable :: SITE_soil_vf_gravels        (:)
   real(r8), allocatable :: SITE_soil_vf_sand           (:)
   real(r8), allocatable :: SITE_soil_vf_clay           (:)
   real(r8), allocatable :: SITE_soil_vf_om             (:)
   real(r8), allocatable :: SITE_soil_wf_gravels        (:)
   real(r8), allocatable :: SITE_soil_wf_sand           (:)
   real(r8), allocatable :: SITE_soil_wf_clay           (:)
   real(r8), allocatable :: SITE_soil_wf_om             (:)
   real(r8), allocatable :: SITE_soil_OM_density        (:)
   real(r8), allocatable :: SITE_soil_BD_all            (:)
   real(r8), allocatable :: SITE_soil_theta_s           (:)
   real(r8), allocatable :: SITE_soil_k_s               (:)
   real(r8), allocatable :: SITE_soil_csol              (:)
   real(r8), allocatable :: SITE_soil_tksatu            (:)
   real(r8), allocatable :: SITE_soil_tksatf            (:)
   real(r8), allocatable :: SITE_soil_tkdry             (:)
   real(r8), allocatable :: SITE_soil_k_solids          (:)
   real(r8), allocatable :: SITE_soil_psi_s             (:)
   real(r8), allocatable :: SITE_soil_lambda            (:)
   real(r8), allocatable :: SITE_soil_theta_r           (:)
#ifdef vanGenuchten_Mualem_SOIL_MODEL
   real(r8), allocatable :: SITE_soil_alpha_vgm         (:)
   real(r8), allocatable :: SITE_soil_L_vgm             (:)
   real(r8), allocatable :: SITE_soil_n_vgm             (:)
#endif
   real(r8), allocatable :: SITE_soil_BA_alpha          (:)
   real(r8), allocatable :: SITE_soil_BA_beta           (:)

   integer  :: SITE_soil_texture

   real(r8) :: SITE_dbedrock  = 0.

   real(r8) :: SITE_elevation  = 0.
   real(r8) :: SITE_elvstd     = 0.
   real(r8) :: SITE_sloperatio = 0.

   ! topography factors used for downscaling
   real(r8) :: SITE_svf = 0.
   real(r8) :: SITE_cur = 0.
   real(r8), allocatable :: SITE_slp_type (:)
   real(r8), allocatable :: SITE_asp_type (:)
   real(r8), allocatable :: SITE_area_type(:)
   real(r8), allocatable :: SITE_sf_lut (:,:)

   logical :: u_site_landtype,   u_site_crop,         u_site_pfts,        u_site_htop,              &
              u_site_lai,        u_site_lakedepth,    u_site_soil_bright, u_site_vf_quartz_mineral, &
              u_site_vf_gravels, u_site_vf_sand,      u_site_vf_clay,     u_site_vf_om,             &
              u_site_wf_gravels, u_site_wf_sand,      u_site_wf_clay,     u_site_wf_om,             &
              u_site_OM_density, u_site_BD_all,       u_site_theta_s,     u_site_k_s,               &
              u_site_csol,       u_site_tksatu,       u_site_tksatf,      u_site_tkdry,             &
              u_site_k_solids,   u_site_psi_s,        u_site_lambda,      u_site_theta_r,           &
              u_site_alpha_vgm,  u_site_L_vgm,        u_site_n_vgm,       u_site_BA_alpha,          &
              u_site_BA_beta,    u_site_soil_texture, u_site_dbedrock,    u_site_elevation,         &
              u_site_elvstd,     u_site_svf,          u_site_cur,         u_site_slp_type,          &
              u_site_asp_type,   u_site_area_type,    u_site_sf_lut,      u_site_sloperatio


   integer  :: SITE_ncar_rid
   integer  :: SITE_urbtyp

   real(r8) :: SITE_lucyid

   real(r8) :: SITE_fveg_urb
   real(r8) :: SITE_htop_urb
   real(r8) :: SITE_flake_urb
   real(r8) :: SITE_froof
   real(r8) :: SITE_hroof
   real(r8) :: SITE_fgimp
   real(r8) :: SITE_fgper
   real(r8) :: SITE_hlr
   real(r8) :: SITE_lambdaw
   real(r8) :: SITE_popden

   real(r8) :: SITE_em_roof
   real(r8) :: SITE_em_wall
   real(r8) :: SITE_em_gimp
   real(r8) :: SITE_em_gper
   real(r8) :: SITE_t_roommax
   real(r8) :: SITE_t_roommin

   real(r8) :: SITE_thickroof
   real(r8) :: SITE_thickwall

   real(r8), allocatable :: SITE_cv_roof (:)
   real(r8), allocatable :: SITE_cv_wall (:)
   real(r8), allocatable :: SITE_cv_gimp (:)
   real(r8), allocatable :: SITE_tk_roof (:)
   real(r8), allocatable :: SITE_tk_wall (:)
   real(r8), allocatable :: SITE_tk_gimp (:)

   real(r8), allocatable :: SITE_alb_roof  (:,:)
   real(r8), allocatable :: SITE_alb_wall  (:,:)
   real(r8), allocatable :: SITE_alb_gimp  (:,:)
   real(r8), allocatable :: SITE_alb_gper  (:,:)

   logical :: u_site_froof, u_site_hroof, u_site_fgper  , u_site_hlr    , &
              u_site_fveg , u_site_htopu, u_site_urblai , u_site_urbsai , &
              u_site_flake, u_site_utype,                                 &
              u_site_albr , u_site_albw , u_site_albgimp, u_site_albgper, &
              u_site_emr  , u_site_emw  , u_site_emgimp , u_site_emgper , &
              u_site_cvr  , u_site_cvw  , u_site_cvgimp ,                 &
              u_site_tkr  , u_site_tkw  , u_site_tkgimp ,                 &
              u_site_tbmax, u_site_tbmin, u_site_thickr , u_site_thickw , &
              u_site_pop  , u_site_lucy

   ! -----------------------------------------------------------------------------------
   ! The soil color and reflectance is from the work:
   ! Peter J. Lawrence and Thomas N. Chase, 2007:
   ! Representing a MODIS consistent land surface in the Community Land Model (CLM 3.0):
   ! Part 1 generating MODIS consistent land surface parameters
   ! -----------------------------------------------------------------------------------
   real(r8), parameter :: soil_s_v_refl(20) = &  ! Saturated visible soil reflectance
      (/ 0.26, 0.24, 0.22, 0.20, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14, &
         0.13, 0.12, 0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04 /)
   real(r8), parameter :: soil_d_v_refl(20) = &  ! Dry visible soil reflectance
      (/ 0.37, 0.35, 0.33, 0.31, 0.30, 0.29, 0.28, 0.27, 0.26, 0.25, &
         0.24, 0.23, 0.22, 0.21, 0.20, 0.19, 0.18, 0.17, 0.16, 0.15 /)
   real(r8), parameter :: soil_s_n_refl(20) = &  ! Saturated near infrared soil reflectance
      (/ 0.52, 0.48, 0.44, 0.40, 0.38, 0.36, 0.34, 0.32, 0.30, 0.28, &
         0.26, 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.08 /)
   real(r8), parameter :: soil_d_n_refl(20) = &  ! Dry near infrared soil reflectance
      (/ 0.63, 0.59, 0.55, 0.51, 0.49, 0.47, 0.45, 0.43, 0.41, 0.39, &
         0.37, 0.35, 0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.19 /)

CONTAINS


!-----------------------------------------------------------------------
   SUBROUTINE read_surface_data_single (fsrfdata, mksrfdata)

   USE MOD_TimeManager
   USE MOD_Grid
   USE MOD_Block
   USE MOD_NetCDFSerial
   USE MOD_NetCDFPoint
   USE MOD_Namelist
   USE MOD_Utils
   USE MOD_Vars_Global, only: PI
   USE MOD_Const_LC
   USE MOD_Const_PFT
   USE MOD_SPMD_Task
   USE MOD_LandPatch
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_LandPFT
#endif
   USE MOD_Mesh, only: numelm
   USE MOD_LandElm
   IMPLICIT NONE

   character(len=*), intent(in) :: fsrfdata
   logical, intent(in) :: mksrfdata

   ! Local Variables
   real(r8) :: lat_in, lon_in
   real(r8) :: LAI, lakedepth, slp, asp, zenith_angle
   integer  :: i, isc, nsl, typ, a, z, arraysize
   integer  :: iyear, idate(3), simulation_lai_year_start, simulation_lai_year_end
   integer  :: start_year, end_year, ntime, itime

   character(len=256) :: filename, dir_5x5, fmt_str
   character(len=4)   :: cyear, c

   type(grid_type) :: gridpatch,  gridcrop, gridpft,  gridhtop, gridlai, gridlake,  &
                      gridbright, gridsoil, gridrock, gridtopo, grid_topo_factor

   integer,  allocatable :: croptyp(:), pfttyp (:)
   real(r8), allocatable :: pctcrop(:), pctpfts(:), pftLAI(:), pftSAI(:), tea_f(:), tea_b(:)

   integer, parameter :: N_PFT_modis = 16
   logical            :: readflag

      CALL Init_GlobalVars
      CALL Init_LC_Const

      IF (mksrfdata) THEN
         write(*,*)
         write(*,*) '  ----------------  Make Single Point Surface Data  ----------------  '
      ENDIF

      IF (ncio_var_exist(fsrfdata, 'latitude')) THEN
         CALL ncio_read_serial (fsrfdata, 'latitude',  lat_in)
         IF ((lat_in /= SITE_lat_location) .and. (SITE_lat_location /= -1.e36_r8)) THEN
            write(*,*) 'Warning: Latitude mismatch: ', &
               lat_in, ' in data file and ', SITE_lat_location, 'in namelist.'
         ENDIF
         SITE_lat_location = lat_in
      ENDIF

      IF (ncio_var_exist(fsrfdata, 'longitude')) THEN
         CALL ncio_read_serial (fsrfdata, 'longitude', lon_in)
         IF ((lon_in /= SITE_lon_location) .and. (SITE_lon_location /= -1.e36_r8)) THEN
            write(*,*) 'Warning: Longitude mismatch: ', &
               lon_in, ' in data file and ', SITE_lon_location, 'in namelist.'
         ENDIF
         SITE_lon_location = lon_in
      ENDIF

      CALL normalize_longitude (SITE_lon_location)

      IF (.not. isgreenwich) THEN
         LocalLongitude = SITE_lon_location
      ENDIF

      IF (mksrfdata) THEN
         write(*,'(A,F8.2)') 'Latitude  : ', SITE_lat_location
         write(*,'(A,F8.2)') 'Longitude : ', SITE_lon_location
      ENDIF


      DEF_domain%edges = floor(SITE_lat_location)
      DEF_domain%edgen = floor(SITE_lat_location) + 1.
      DEF_domain%edgew = floor(SITE_lon_location)
      DEF_domain%edgee = floor(SITE_lon_location) + 1.

      CALL gblock%set ()
      gblock%nblkme = 1
      allocate(gblock%xblkme(1))
      allocate(gblock%yblkme(1))
      gblock%xblkme(1) = find_nearest_west  (SITE_lon_location, gblock%nxblk, gblock%lon_w)
      gblock%yblkme(1) = find_nearest_south (SITE_lat_location, gblock%nyblk, gblock%lat_s)


      ! (1) build/read "land patch" by using land cover type data
      numpatch = 1

#ifdef LULC_USGS
      u_site_landtype = (.not. mksrfdata) .or. (SITE_landtype >= 0) &
         .or. (USE_SITE_landtype .and. ncio_var_exist(fsrfdata,'USGS_classification'))

      IF (u_site_landtype) THEN
         IF (SITE_landtype == -1) THEN
            CALL ncio_read_serial (fsrfdata, 'USGS_classification', SITE_landtype)
         ENDIF
      ELSE
         CALL gridpatch%define_by_name ('colm_1km')
         filename = trim(DEF_dir_rawdata)//'/landtypes/landtype-usgs-update.nc'
         CALL read_point_var_2d_int32 (gridpatch, filename, 'landtype', &
            SITE_lon_location, SITE_lat_location, SITE_landtype)
      ENDIF
#else
      u_site_landtype = (.not. mksrfdata) .or. (SITE_landtype >= 0) &
         .or. (USE_SITE_landtype .and. ncio_var_exist(fsrfdata,'IGBP_classification'))

      IF (u_site_landtype)  THEN
         IF (SITE_landtype == -1) THEN
            CALL ncio_read_serial (fsrfdata, 'IGBP_classification', SITE_landtype)
         ENDIF
      ELSE
         CALL gridpatch%define_by_name ('colm_500m')
         write(cyear,'(i4.4)') DEF_LC_YEAR
         filename = trim(DEF_dir_rawdata)//'landtypes/landtype-igbp-modis-'//trim(cyear)//'.nc'
         CALL read_point_var_2d_int32 (gridpatch, filename, 'landtype', &
            SITE_lon_location, SITE_lat_location, SITE_landtype)
      ENDIF
#endif

      IF (SITE_landtype < 0) THEN
         write(*,*) 'Error! Please set SITE_landtype in namelist file !'
         CALL CoLM_stop()
      ENDIF

#ifdef URBAN_MODEL
      IF (SITE_landtype /= URBAN) THEN
         write(*,*) 'Error! Please set SITE_landtype to URBAN in namelist file !'
         CALL CoLM_stop()
      ENDIF
#endif

      IF (mksrfdata) THEN
         write(*,'(A,A,3A)') 'Land cover type : ', trim(patchclassname(SITE_landtype)), &
            ' (from ',trim(datasource(u_site_landtype)),')'
      ENDIF

      ! (2) build/read "land crop" by using crop data
#ifdef CROP
      IF (SITE_landtype == CROPLAND) THEN

         readflag = ((.not. mksrfdata) .or. USE_SITE_pctcrop)
         u_site_crop = readflag &
            .and. ncio_var_exist(fsrfdata,'croptyp',readflag) .and. ncio_var_exist(fsrfdata,'pctcrop',readflag)

         IF (u_site_crop) THEN
            CALL ncio_read_serial (fsrfdata, 'croptyp', croptyp)
            CALL ncio_read_serial (fsrfdata, 'pctcrop', pctcrop)
         ELSE
            allocate (croptyp (N_CFT))
            croptyp = (/(i, i = 1, N_CFT)/)

            filename = trim(DEF_dir_rawdata) // '/global_CFT_surface_data.nc'
            CALL gridcrop%define_from_file (filename, 'lat', 'lon')
            CALL read_point_var_3d_first_real8 (gridcrop, filename, 'PCT_CFT', &
               SITE_lon_location, SITE_lat_location, N_CFT, pctcrop)
         ENDIF

         numpatch = count(pctcrop > 0.)

         IF (numpatch == 0) THEN
            write(*,*) 'There is no crop at this point!'
            CALL CoLM_stop()
         ENDIF

         allocate (SITE_croptyp (numpatch))
         allocate (SITE_pctcrop (numpatch))

         SITE_croptyp = pack(croptyp, pctcrop > 0.)
         SITE_pctcrop = pack(pctcrop, pctcrop > 0.) / sum(pctcrop)

         IF (mksrfdata) THEN
            write(c,'(I0)') numpatch
            write(*,'(A,'//trim(c)//'I5,3A)')   'crop type : ', SITE_croptyp, ' (from ',trim(datasource(u_site_crop)),')'
            write(*,'(A,'//trim(c)//'F5.2,3A)') 'crop frac : ', SITE_pctcrop, ' (from ',trim(datasource(u_site_crop)),')'
         ENDIF

      ENDIF
#endif


      ! (3) build/read "land pft" or "land pc" by using plant functional type data
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
#ifndef CROP
      IF (patchtypes(SITE_landtype) == 0) THEN
#else
      IF (patchtypes(SITE_landtype) == 0 .and. SITE_landtype /= CROPLAND) THEN
#endif
         readflag = ((.not. mksrfdata) .or. USE_SITE_pctpfts)
         u_site_pfts = readflag &
            .and. ncio_var_exist(fsrfdata,'pfttyp',readflag) .and. ncio_var_exist(fsrfdata,'pctpfts',readflag)

         IF (u_site_pfts) THEN
            CALL ncio_read_serial (fsrfdata, 'pfttyp',  pfttyp )
            CALL ncio_read_serial (fsrfdata, 'pctpfts', pctpfts)
         ELSE
            allocate (pfttyp (N_PFT_modis))
            pfttyp = (/(i, i = 0, N_PFT_modis-1)/)

            CALL gridpft%define_by_name ('colm_500m')

            dir_5x5 = trim(DEF_dir_rawdata) // '/plant_15s'
            write(cyear,'(i4.4)') DEF_LC_YEAR
            CALL read_point_5x5_var_3d_real8 (gridpft, dir_5x5, 'MOD'//trim(cyear), 'PCT_PFT', &
               SITE_lon_location, SITE_lat_location, N_PFT_modis, pctpfts)
         ENDIF

         numpft = count(pctpfts > 0.)

         allocate (SITE_pfttyp  (numpft))
         allocate (SITE_pctpfts (numpft))

         SITE_pfttyp  = pack(pfttyp,  pctpfts > 0.)
         SITE_pctpfts = pack(pctpfts, pctpfts > 0.) / sum(pctpfts)

#ifdef CROP
      ELSEIF (SITE_landtype == CROPLAND) THEN
         u_site_pfts = .false.
         numpft = numpatch
         allocate (SITE_pfttyp  (numpft))
         allocate (SITE_pctpfts (numpft))
         SITE_pfttyp  = SITE_croptyp + N_PFT - 1
         SITE_pctpfts = 1.
#endif
      ELSE
         numpft = 0
      ENDIF

      IF ((patchtypes(SITE_landtype) == 0) .and. (numpft == 0)) THEN
         write(*,*) 'Warning : There is no plant functional type at this site !    '
         CALL CoLM_stop()
      ENDIF

      IF (mksrfdata) THEN
         IF (numpft > 0) THEN
            write(*,'(4A)') 'PFT type and fraction : ', ' (from ',trim(datasource(u_site_pfts)),')'
            DO i = 1, numpft
               write(*,'(A,F5.2)') '  '//trim(pftclassname(SITE_pfttyp(i))), SITE_pctpfts(i)
            ENDDO
         ENDIF
      ENDIF

#endif


      ! (4) forest height
      readflag = (.not. mksrfdata) .or. USE_SITE_htop
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      IF (patchtypes(SITE_landtype) == 0) THEN
         u_site_htop = readflag .and. ncio_var_exist(fsrfdata,'canopy_height_pfts',readflag)
      ELSE
         u_site_htop = readflag .and. ncio_var_exist(fsrfdata,'canopy_height',readflag)
      ENDIF
#else
      u_site_htop = readflag .and. ncio_var_exist(fsrfdata,'canopy_height',readflag)
#endif

      IF (u_site_htop) THEN
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
         IF (patchtypes(SITE_landtype) == 0) THEN
            CALL ncio_read_serial (fsrfdata, 'canopy_height_pfts', SITE_htop_pfts)
         ELSE
            CALL ncio_read_serial (fsrfdata, 'canopy_height', SITE_htop)
         ENDIF
#else
         CALL ncio_read_serial (fsrfdata, 'canopy_height', SITE_htop)
#endif
      ELSE
#ifdef LULC_USGS
         CALL gridhtop%define_by_name ('colm_1km')
         filename = trim(DEF_dir_rawdata)//'/Forest_Height.nc'
         CALL read_point_var_2d_real8 (gridhtop, filename, 'forest_height', &
            SITE_lon_location, SITE_lat_location, SITE_htop)
#else

         CALL gridhtop%define_by_name ('colm_500m')

         dir_5x5 = trim(DEF_dir_rawdata) // '/plant_15s'
         write(cyear,'(i4.4)') DEF_LC_YEAR
         CALL read_point_5x5_var_2d_real8 (gridhtop, dir_5x5, 'MOD'//trim(cyear), 'HTOP', &
            SITE_lon_location, SITE_lat_location, SITE_htop)

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
         IF (numpft > 0) THEN
            allocate (SITE_htop_pfts (numpft))
            SITE_htop_pfts(:) = SITE_htop
         ENDIF
#endif
#endif
      ENDIF

      IF (mksrfdata) THEN
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
         IF (patchtypes(SITE_landtype) == 0) THEN
            arraysize = size(SITE_htop_pfts)
            write(fmt_str, '("(A,", I0, "F8.2,3A)")') arraysize
            write(*,fmt_str) 'Forest height : ', SITE_htop_pfts, ' (from ',trim(datasource(u_site_htop)),')'
         ELSE
            write(*,'(A,F8.2,3A)') 'Forest height : ', SITE_htop, ' (from ',trim(datasource(u_site_htop)),')'
         ENDIF
#else
         write(*,'(A,F8.2,3A)') 'Forest height : ', SITE_htop, ' (from ',trim(datasource(u_site_htop)),')'
#endif
      ENDIF


      ! (5) LAI
      readflag = ((.not. mksrfdata) .or. USE_SITE_LAI)
      readflag = readflag .and. ncio_var_exist(fsrfdata,'LAI_year',readflag)
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      IF (patchtypes(SITE_landtype) == 0) THEN
         u_site_lai = readflag .and. ncio_var_exist(fsrfdata,'LAI_pfts_monthly',readflag) &
            .and. ncio_var_exist(fsrfdata,'SAI_pfts_monthly',readflag)
      ELSE
         u_site_lai = readflag .and. ncio_var_exist(fsrfdata,'LAI_monthly',readflag) &
            .and. ncio_var_exist(fsrfdata,'SAI_monthly',readflag)
      ENDIF
#else
      IF (DEF_LAI_MONTHLY) THEN
         u_site_lai = readflag .and. ncio_var_exist(fsrfdata,'LAI_monthly',readflag) &
            .and. ncio_var_exist(fsrfdata,'SAI_monthly',readflag)
      ELSE
         u_site_lai = readflag .and. ncio_var_exist(fsrfdata,'LAI_8day',readflag)
      ENDIF
#endif

      IF (u_site_lai) THEN
         CALL ncio_read_serial (fsrfdata, 'LAI_year', SITE_LAI_year)
         start_year = 1
         end_year   = size(SITE_LAI_year)
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
         IF (patchtypes(SITE_landtype) == 0) THEN
            CALL ncio_read_serial (fsrfdata, 'LAI_pfts_monthly', SITE_LAI_pfts_monthly)
            CALL ncio_read_serial (fsrfdata, 'SAI_pfts_monthly', SITE_SAI_pfts_monthly)
            ntime = size(SITE_LAI_pfts_monthly,2)
         ELSE
            CALL ncio_read_serial (fsrfdata, 'LAI_monthly', SITE_LAI_monthly)
            CALL ncio_read_serial (fsrfdata, 'SAI_monthly', SITE_SAI_monthly)
            ntime = size(SITE_LAI_monthly,1)
         ENDIF
#else
         IF (DEF_LAI_MONTHLY) THEN
            CALL ncio_read_serial (fsrfdata, 'LAI_monthly', SITE_LAI_monthly)
            CALL ncio_read_serial (fsrfdata, 'SAI_monthly', SITE_SAI_monthly)
            ntime = size(SITE_LAI_monthly,1)
         ELSE
            CALL ncio_read_serial (fsrfdata, 'LAI_8day', SITE_LAI_8day)
            ntime = size(SITE_LAI_8day,1)
         ENDIF
#endif
      ELSE

         idate(1) = DEF_simulation_time%start_year
         IF (.not. isgreenwich) THEN
            idate(3) = DEF_simulation_time%start_sec
            CALL monthday2julian (idate(1), &
               DEF_simulation_time%start_month, DEF_simulation_time%start_day, idate(2))
            CALL localtime2gmt(idate)
         ENDIF

         simulation_lai_year_start = idate(1)

         idate(1) = DEF_simulation_time%end_year
         IF (.not. isgreenwich) THEN
            idate(3) = DEF_simulation_time%end_sec
            CALL monthday2julian (idate(1), &
               DEF_simulation_time%end_month, DEF_simulation_time%end_day, idate(2))
            CALL localtime2gmt(idate)
         ENDIF

         simulation_lai_year_end = idate(1)

         IF (DEF_LAI_CHANGE_YEARLY) THEN
            start_year = max(simulation_lai_year_start, DEF_LAI_START_YEAR )
            start_year = min(start_year,                DEF_LAI_END_YEAR   )
            end_year   = min(simulation_lai_year_end,   DEF_LAI_END_YEAR   )
            end_year   = max(end_year,                  DEF_LAI_START_YEAR )
         ELSE
            start_year = DEF_LC_YEAR
            end_year   = DEF_LC_YEAR
         ENDIF

         allocate (SITE_LAI_year (start_year:end_year))
         SITE_LAI_year = (/(iyear, iyear = start_year, end_year)/)

         IF (DEF_LAI_MONTHLY) THEN
            ntime = 12
            allocate (SITE_LAI_monthly (12,start_year:end_year))
            allocate (SITE_SAI_monthly (12,start_year:end_year))
         ELSE
            ntime = 46
            allocate (SITE_LAI_8day    (46,start_year:end_year))
         ENDIF

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
         IF (numpft > 0) THEN
            allocate (SITE_LAI_pfts_monthly (numpft,12,start_year:end_year))
            allocate (SITE_SAI_pfts_monthly (numpft,12,start_year:end_year))
         ENDIF
#endif

         CALL gridlai%define_by_name ('colm_500m')

         DO iyear = start_year, end_year

            write(cyear,'(i4.4)') iyear

            DO itime = 1, ntime
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)

               dir_5x5 = trim(DEF_dir_rawdata) // '/plant_15s'
               CALL read_point_5x5_var_3d_time_real8 (gridlai, dir_5x5, 'MOD'//trim(cyear), 'MONTHLY_PFT_LAI', &
                  SITE_lon_location, SITE_lat_location, N_PFT_modis, itime, pftLAI)
               CALL read_point_5x5_var_3d_time_real8 (gridlai, dir_5x5, 'MOD'//trim(cyear), 'MONTHLY_PFT_SAI', &
                  SITE_lon_location, SITE_lat_location, N_PFT_modis, itime, pftSAI)

#ifndef CROP
               IF (patchtypes(SITE_landtype) == 0) THEN
#else
               IF (patchtypes(SITE_landtype) == 0 .and. SITE_landtype /= CROPLAND) THEN
#endif
                  IF (allocated(pctpfts)) deallocate (pctpfts)
                  allocate(pctpfts (0:N_PFT_modis-1)); pctpfts(:) = 0.
                  pctpfts(SITE_pfttyp) = SITE_pctpfts

                  SITE_LAI_pfts_monthly(:,itime,iyear) = pack(pftLAI, pctpfts > 0.)
                  SITE_SAI_pfts_monthly(:,itime,iyear) = pack(pftSAI, pctpfts > 0.)
#ifdef CROP
               ELSEIF (SITE_landtype == CROPLAND) THEN
                  CALL read_point_5x5_var_3d_real8 (gridlai, dir_5x5, 'MOD'//trim(cyear), 'PCT_PFT', &
                     SITE_lon_location, SITE_lat_location, N_PFT_modis, pctpfts)
                  SITE_LAI_pfts_monthly(:,itime,iyear) = sum(pftLAI * pctpfts) / sum(pctpfts)
                  SITE_SAI_pfts_monthly(:,itime,iyear) = sum(pftSAI * pctpfts) / sum(pctpfts)
#endif
               ELSE
                  dir_5x5 = trim(DEF_dir_rawdata) // '/plant_15s'
                  CALL read_point_5x5_var_2d_time_real8 (gridlai, dir_5x5, 'MOD'//trim(cyear), &
                     'MONTHLY_LC_LAI', SITE_lon_location, SITE_lat_location, itime, &
                     SITE_LAI_monthly(itime,iyear))
                  CALL read_point_5x5_var_2d_time_real8 (gridlai, dir_5x5, 'MOD'//trim(cyear), &
                     'MONTHLY_LC_SAI', SITE_lon_location, SITE_lat_location, itime, &
                     SITE_SAI_monthly(itime,iyear))
               ENDIF

#else
               IF (DEF_LAI_MONTHLY) THEN
                  dir_5x5 = trim(DEF_dir_rawdata) // '/plant_15s'
                  CALL read_point_5x5_var_2d_time_real8 (gridlai, dir_5x5, 'MOD'//trim(cyear), &
                     'MONTHLY_LC_LAI', SITE_lon_location, SITE_lat_location, itime, &
                     SITE_LAI_monthly(itime,iyear))
                  CALL read_point_5x5_var_2d_time_real8 (gridlai, dir_5x5, 'MOD'//trim(cyear), &
                     'MONTHLY_LC_SAI', SITE_lon_location, SITE_lat_location, itime, &
                     SITE_SAI_monthly(itime,iyear))
               ELSE
                  filename = trim(DEF_dir_rawdata)//'/lai_15s_8day/lai_8-day_15s_'//trim(cyear)//'.nc'
                  CALL read_point_var_2d_time_real8 (gridlai, filename, 'lai', &
                     SITE_lon_location, SITE_lat_location, itime, LAI)
                  SITE_LAI_8day(itime,iyear) = LAI * 0.1
               ENDIF
#endif
            ENDDO
         ENDDO
      ENDIF

      IF (mksrfdata) THEN
         DO iyear = start_year, end_year
            write(c,'(i2)') ntime
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
            IF (patchtypes(SITE_landtype) == 0) THEN
               DO i = 1, numpft
                  write(*,'(A,I4,A,I2,A,'//trim(c)//'F8.2,4A)') 'LAI (year ', SITE_LAI_year(iyear), &
                     ', pft ', SITE_pfttyp(i),') : ', SITE_LAI_pfts_monthly(i,:,iyear), &
                     ' (from ',trim(datasource(u_site_lai)),')'
                  write(*,'(A,I4,A,I2,A,'//trim(c)//'F8.2,4A)') 'SAI (year ', SITE_LAI_year(iyear), &
                     ', pft ', SITE_pfttyp(i),') : ', SITE_SAI_pfts_monthly(i,:,iyear), &
                     ' (from ',trim(datasource(u_site_lai)),')'
               ENDDO
            ELSE
               write(*,'(A,I4,A,'//trim(c)//'F8.2,4A)') 'LAI (year ', SITE_LAI_year(iyear), ') : ', &
                  SITE_LAI_monthly(:,iyear), ' (from ',trim(datasource(u_site_lai)),')'
               write(*,'(A,I4,A,'//trim(c)//'F8.2,4A)') 'SAI (year ', SITE_LAI_year(iyear), ') : ', &
                  SITE_SAI_monthly(:,iyear), ' (from ',trim(datasource(u_site_lai)),')'
            ENDIF
#else
            IF (DEF_LAI_MONTHLY) THEN
               write(*,'(A,I4,A,'//trim(c)//'F8.2,4A)') 'LAI (year ', SITE_LAI_year(iyear), ') : ', &
                  SITE_LAI_monthly(:,iyear), ' (from ',trim(datasource(u_site_lai)),')'
               write(*,'(A,I4,A,'//trim(c)//'F8.2,4A)') 'SAI (year ', SITE_LAI_year(iyear), ') : ', &
                  SITE_SAI_monthly(:,iyear), ' (from ',trim(datasource(u_site_lai)),')'
            ELSE
               write(*,'(A,I4,A,'//trim(c)//'F8.2,4A)') 'LAI (year ', SITE_LAI_year(iyear), ') : ', &
                  SITE_LAI_8day(:,iyear), ' (from ',trim(datasource(u_site_lai)),')'
            ENDIF
#endif
         ENDDO
      ENDIF


      ! (6) lake depth
      readflag = ((.not. mksrfdata) .or. USE_SITE_lakedepth)
      u_site_lakedepth = readflag .and. ncio_var_exist(fsrfdata,'lakedepth',readflag)

      IF (u_site_lakedepth) THEN
         CALL ncio_read_serial (fsrfdata, 'lakedepth', SITE_lakedepth)
      ELSE
         CALL gridlake%define_by_name ('colm_500m')
         filename = trim(DEF_dir_rawdata)//'/lake_depth.nc'
         CALL read_point_var_2d_real8 (gridlake, filename, 'lake_depth', &
            SITE_lon_location, SITE_lat_location, lakedepth)
         SITE_lakedepth = lakedepth * 0.1
      ENDIF

      IF (mksrfdata) THEN
         write(*,'(A,F8.2,3A)') 'Lake depth : ', SITE_lakedepth, ' (from ',trim(datasource(u_site_lakedepth)),')'
      ENDIF


      ! (7) soil brightness parameters
      readflag = ((.not. mksrfdata) .or. USE_SITE_soilreflectance)
      u_site_soil_bright = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_s_v_alb',readflag) &
         .and. ncio_var_exist(fsrfdata,'soil_d_v_alb',readflag) &
         .and. ncio_var_exist(fsrfdata,'soil_s_n_alb',readflag) &
         .and. ncio_var_exist(fsrfdata,'soil_d_n_alb',readflag)

      IF (u_site_soil_bright) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_s_v_alb', SITE_soil_s_v_alb)
         CALL ncio_read_serial (fsrfdata, 'soil_d_v_alb', SITE_soil_d_v_alb)
         CALL ncio_read_serial (fsrfdata, 'soil_s_n_alb', SITE_soil_s_n_alb)
         CALL ncio_read_serial (fsrfdata, 'soil_d_n_alb', SITE_soil_d_n_alb)
      ELSE
         SITE_soil_s_v_alb = spval
         SITE_soil_d_v_alb = spval
         SITE_soil_s_n_alb = spval
         SITE_soil_d_n_alb = spval

         CALL gridbright%define_by_name ('colm_500m')
         filename = trim(DEF_dir_rawdata)//'/soil_brightness.nc'
         CALL read_point_var_2d_int32 (gridbright, filename, 'soil_brightness', &
            SITE_lon_location, SITE_lat_location, isc)

#ifdef LULC_USGS
         IF(SITE_landtype /= 16 .and. SITE_landtype /= 24)THEN  ! NOT WATER BODIES(16)/GLACIER and ICESHEET(24)
#else
         IF(SITE_landtype /= 17 .and. SITE_landtype /= 15)THEN  ! NOT WATER BODIES(17)/GLACIER and ICE SHEET(15)
#endif
            IF ((isc >= 1) .and. (isc <= 20)) THEN
               SITE_soil_s_v_alb = soil_s_v_refl( isc )
               SITE_soil_d_v_alb = soil_d_v_refl( isc )
               SITE_soil_s_n_alb = soil_s_n_refl( isc )
               SITE_soil_d_n_alb = soil_d_n_refl( isc )
            ENDIF
         ENDIF
      ENDIF

      IF (mksrfdata) THEN
         write(*,'(A,F8.2,3A)') 'Soil brightness s_v : ', SITE_soil_s_v_alb, ' (from ',trim(datasource(u_site_soil_bright)),')'
         write(*,'(A,F8.2,3A)') 'Soil brightness d_v : ', SITE_soil_d_v_alb, ' (from ',trim(datasource(u_site_soil_bright)),')'
         write(*,'(A,F8.2,3A)') 'Soil brightness s_n : ', SITE_soil_s_n_alb, ' (from ',trim(datasource(u_site_soil_bright)),')'
         write(*,'(A,F8.2,3A)') 'Soil brightness d_n : ', SITE_soil_d_n_alb, ' (from ',trim(datasource(u_site_soil_bright)),')'
      ENDIF


      ! (8) soil parameters

      CALL gridsoil%define_by_name ('colm_500m')

      readflag = ((.not. mksrfdata) .or. USE_SITE_soilparameters)
      u_site_vf_quartz_mineral = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_vf_quartz_mineral',readflag)
      IF (u_site_vf_quartz_mineral) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_vf_quartz_mineral', SITE_soil_vf_quartz_mineral)
      ELSE
         allocate (SITE_soil_vf_quartz_mineral (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/vf_quartz_mineral_s.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'vf_quartz_mineral_s_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_vf_quartz_mineral(nsl))
         ENDDO
      ENDIF

      u_site_vf_gravels = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_vf_gravels',readflag)
      IF (u_site_vf_gravels) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_vf_gravels', SITE_soil_vf_gravels)
      ELSE
         allocate (SITE_soil_vf_gravels (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/vf_gravels_s.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'vf_gravels_s_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_vf_gravels(nsl))
         ENDDO
      ENDIF

      u_site_vf_sand = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_vf_sand',readflag)
      IF (u_site_vf_sand) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_vf_sand', SITE_soil_vf_sand)
      ELSE
         allocate (SITE_soil_vf_sand (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/vf_sand_s.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'vf_sand_s_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_vf_sand(nsl))
         ENDDO
      ENDIF

      u_site_vf_clay = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_vf_clay',readflag)
      IF (u_site_vf_clay) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_vf_clay', SITE_soil_vf_clay)
      ELSE
         allocate (SITE_soil_vf_clay (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/vf_clay_s.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'vf_clay_s_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_vf_clay(nsl))
         ENDDO
      ENDIF

      u_site_vf_om = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_vf_om',readflag)
      IF (u_site_vf_om) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_vf_om', SITE_soil_vf_om)
      ELSE
         allocate (SITE_soil_vf_om (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/vf_om_s.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'vf_om_s_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_vf_om(nsl))
         ENDDO
      ENDIF

      u_site_wf_gravels = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_wf_gravels',readflag)
      IF (u_site_wf_gravels) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_wf_gravels', SITE_soil_wf_gravels)
      ELSE
         allocate (SITE_soil_wf_gravels (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/wf_gravels_s.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'wf_gravels_s_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_wf_gravels(nsl))
         ENDDO
      ENDIF

      u_site_wf_sand = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_wf_sand',readflag)
      IF (u_site_wf_sand) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_wf_sand', SITE_soil_wf_sand)
      ELSE
         allocate (SITE_soil_wf_sand (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/wf_sand_s.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'wf_sand_s_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_wf_sand(nsl))
         ENDDO
      ENDIF

      u_site_wf_clay = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_wf_clay',readflag)
      IF (u_site_wf_clay) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_wf_clay', SITE_soil_wf_clay)
      ELSE
         allocate (SITE_soil_wf_clay (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/wf_clay_s.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'wf_clay_s_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_wf_clay(nsl))
         ENDDO
      ENDIF

      u_site_wf_om = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_wf_om',readflag)
      IF (u_site_wf_om) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_wf_om', SITE_soil_wf_om)
      ELSE
         allocate (SITE_soil_wf_om (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/wf_om_s.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'wf_om_s_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_wf_om(nsl))
         ENDDO
      ENDIF

      u_site_OM_density = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_OM_density',readflag)
      IF (u_site_OM_density) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_OM_density', SITE_soil_OM_density)
      ELSE
         allocate (SITE_soil_OM_density (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/OM_density_s.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'OM_density_s_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_OM_density(nsl))
         ENDDO
      ENDIF

      u_site_BD_all = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_BD_all',readflag)
      IF (u_site_BD_all) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_BD_all', SITE_soil_BD_all)
      ELSE
         allocate (SITE_soil_BD_all (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/BD_all_s.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'BD_all_s_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_BD_all(nsl))
         ENDDO
      ENDIF

      u_site_theta_s = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_theta_s',readflag)
      IF (u_site_theta_s) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_theta_s', SITE_soil_theta_s)
      ELSE
         allocate (SITE_soil_theta_s (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/theta_s.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'theta_s_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_theta_s(nsl))
         ENDDO
      ENDIF

      u_site_k_s = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_k_s',readflag)
      IF (u_site_k_s) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_k_s', SITE_soil_k_s)
      ELSE
         allocate (SITE_soil_k_s (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/k_s.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'k_s_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_k_s(nsl))
         ENDDO
      ENDIF

      u_site_csol = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_csol',readflag)
      IF (u_site_csol) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_csol', SITE_soil_csol)
      ELSE
         allocate (SITE_soil_csol (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/csol.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'csol_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_csol(nsl))
         ENDDO
      ENDIF

      u_site_tksatu = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_tksatu',readflag)
      IF (u_site_tksatu) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_tksatu', SITE_soil_tksatu)
      ELSE
         allocate (SITE_soil_tksatu (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/tksatu.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'tksatu_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_tksatu(nsl))
         ENDDO
      ENDIF

      u_site_tksatf = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_tksatf',readflag)
      IF (u_site_tksatf) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_tksatf', SITE_soil_tksatf)
      ELSE
         allocate (SITE_soil_tksatf (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/tksatf.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'tksatf_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_tksatf(nsl))
         ENDDO
      ENDIF

      u_site_tkdry = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_tkdry',readflag)
      IF (u_site_tkdry) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_tkdry', SITE_soil_tkdry)
      ELSE
         allocate (SITE_soil_tkdry (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/tkdry.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'tkdry_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_tkdry(nsl))
         ENDDO
      ENDIF

      u_site_k_solids = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_k_solids',readflag)
      IF (u_site_k_solids) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_k_solids', SITE_soil_k_solids)
      ELSE
         allocate (SITE_soil_k_solids (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/k_solids.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'k_solids_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_k_solids(nsl))
         ENDDO
      ENDIF

      u_site_psi_s = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_psi_s',readflag)
      IF (u_site_psi_s) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_psi_s', SITE_soil_psi_s)
      ELSE
         allocate (SITE_soil_psi_s (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/psi_s.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'psi_s_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_psi_s(nsl))
         ENDDO
      ENDIF

      u_site_lambda = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_lambda',readflag)
      IF (u_site_lambda) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_lambda', SITE_soil_lambda)
      ELSE
         allocate (SITE_soil_lambda (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/lambda.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'lambda_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_lambda(nsl))
         ENDDO
      ENDIF

#ifdef vanGenuchten_Mualem_SOIL_MODEL
      u_site_theta_r = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_theta_r',readflag)
      IF (u_site_theta_r) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_theta_r', SITE_soil_theta_r)
      ELSE
         allocate (SITE_soil_theta_r (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/VGM_theta_r.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'VGM_theta_r_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_theta_r(nsl))
         ENDDO
      ENDIF

      u_site_alpha_vgm = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_alpha_vgm',readflag)
      IF (u_site_alpha_vgm) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_alpha_vgm', SITE_soil_alpha_vgm)
      ELSE
         allocate (SITE_soil_alpha_vgm (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/VGM_alpha.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'VGM_alpha_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_alpha_vgm(nsl))
         ENDDO
      ENDIF

      u_site_L_vgm = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_L_vgm',readflag)
      IF (u_site_L_vgm) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_L_vgm', SITE_soil_L_vgm)
      ELSE
         allocate (SITE_soil_L_vgm (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/VGM_L.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'VGM_L_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_L_vgm(nsl))
         ENDDO
      ENDIF

      u_site_n_vgm = readflag &
         .and. ncio_var_exist(fsrfdata,'soil_n_vgm',readflag)
      IF (u_site_n_vgm) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_n_vgm', SITE_soil_n_vgm)
      ELSE
         allocate (SITE_soil_n_vgm (8))
         DO nsl = 1, 8
            write(c,'(i1)') nsl
            filename = trim(DEF_dir_rawdata)//'/soil/VGM_n.nc'
            CALL read_point_var_2d_real8 (gridsoil, filename, 'VGM_n_l'//trim(c), &
               SITE_lon_location, SITE_lat_location, SITE_soil_n_vgm(nsl))
         ENDDO
      ENDIF
#endif

      u_site_BA_alpha = u_site_vf_gravels .and. u_site_vf_sand
      u_site_BA_beta  = u_site_vf_gravels .and. u_site_vf_sand
      allocate (SITE_soil_BA_alpha (8))
      allocate (SITE_soil_BA_beta  (8))
      DO nsl = 1, 8
         IF (SITE_soil_vf_gravels(nsl) + SITE_soil_vf_sand(nsl) > 0.4) THEN
            SITE_soil_BA_alpha(nsl) = 0.38
            SITE_soil_BA_beta (nsl) = 35.0
         ELSEIF (SITE_soil_vf_gravels(nsl) + SITE_soil_vf_sand(nsl) > 0.25) THEN
            SITE_soil_BA_alpha(nsl) = 0.24
            SITE_soil_BA_beta (nsl) = 26.0
         ELSE
            SITE_soil_BA_alpha(nsl) = 0.20
            SITE_soil_BA_beta (nsl) = 10.0
         ENDIF
      ENDDO


      IF (DEF_Runoff_SCHEME == 3) THEN ! for Simple VIC
         u_site_soil_texture = readflag &
            .and. ncio_var_exist(fsrfdata,'soil_texture',readflag)
         IF (u_site_soil_texture) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_texture', SITE_soil_texture)
         ELSE
            filename = trim(DEF_dir_rawdata)//'/soil/soiltexture_0cm-60cm_mean.nc'
            CALL read_point_var_2d_int32 (gridsoil, filename, 'soiltexture', &
               SITE_lon_location, SITE_lat_location, SITE_soil_texture)
         ENDIF
      ENDIF

      IF (mksrfdata) THEN
         write(*,'(A,8ES10.2,3A)') 'soil_vf_quartz_mineral : ', SITE_soil_vf_quartz_mineral(1:8), ' (from ',trim(datasource(u_site_vf_quartz_mineral)),')'
         write(*,'(A,8ES10.2,3A)') 'soil_vf_gravels        : ', SITE_soil_vf_gravels       (1:8), ' (from ',trim(datasource(u_site_vf_gravels       )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_vf_sand           : ', SITE_soil_vf_sand          (1:8), ' (from ',trim(datasource(u_site_vf_sand          )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_vf_clay           : ', SITE_soil_vf_clay          (1:8), ' (from ',trim(datasource(u_site_vf_clay          )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_vf_om             : ', SITE_soil_vf_om            (1:8), ' (from ',trim(datasource(u_site_vf_om            )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_wf_gravels        : ', SITE_soil_wf_gravels       (1:8), ' (from ',trim(datasource(u_site_wf_gravels       )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_wf_sand           : ', SITE_soil_wf_sand          (1:8), ' (from ',trim(datasource(u_site_wf_sand          )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_wf_clay           : ', SITE_soil_wf_clay          (1:8), ' (from ',trim(datasource(u_site_wf_clay          )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_wf_om             : ', SITE_soil_wf_om            (1:8), ' (from ',trim(datasource(u_site_wf_om            )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_OM_density        : ', SITE_soil_OM_density       (1:8), ' (from ',trim(datasource(u_site_OM_density       )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_BD_all            : ', SITE_soil_BD_all           (1:8), ' (from ',trim(datasource(u_site_BD_all           )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_theta_s           : ', SITE_soil_theta_s          (1:8), ' (from ',trim(datasource(u_site_theta_s          )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_k_s               : ', SITE_soil_k_s              (1:8), ' (from ',trim(datasource(u_site_k_s              )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_csol              : ', SITE_soil_csol             (1:8), ' (from ',trim(datasource(u_site_csol             )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_tksatu            : ', SITE_soil_tksatu           (1:8), ' (from ',trim(datasource(u_site_tksatu           )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_tksatf            : ', SITE_soil_tksatf           (1:8), ' (from ',trim(datasource(u_site_tksatf           )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_tkdry             : ', SITE_soil_tkdry            (1:8), ' (from ',trim(datasource(u_site_tkdry            )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_k_solids          : ', SITE_soil_k_solids         (1:8), ' (from ',trim(datasource(u_site_k_solids         )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_psi_s             : ', SITE_soil_psi_s            (1:8), ' (from ',trim(datasource(u_site_psi_s            )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_lambda            : ', SITE_soil_lambda           (1:8), ' (from ',trim(datasource(u_site_lambda           )),')'
#ifdef vanGenuchten_Mualem_SOIL_MODEL
         write(*,'(A,8ES10.2,3A)') 'soil_theta_r           : ', SITE_soil_theta_r          (1:8), ' (from ',trim(datasource(u_site_theta_r          )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_alpha_vgm         : ', SITE_soil_alpha_vgm        (1:8), ' (from ',trim(datasource(u_site_alpha_vgm        )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_L_vgm             : ', SITE_soil_L_vgm            (1:8), ' (from ',trim(datasource(u_site_L_vgm            )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_n_vgm             : ', SITE_soil_n_vgm            (1:8), ' (from ',trim(datasource(u_site_n_vgm            )),')'
#endif
         write(*,'(A,8ES10.2,3A)') 'soil_BA_alpha          : ', SITE_soil_BA_alpha         (1:8), ' (from ',trim(datasource(u_site_BA_alpha         )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_BA_beta           : ', SITE_soil_BA_beta          (1:8), ' (from ',trim(datasource(u_site_BA_beta          )),')'

         IF (DEF_Runoff_SCHEME == 3) THEN ! for Simple VIC
            write(*,'(A,I3,3A)') 'soil texture           : ', SITE_soil_texture, ' (from ',trim(datasource(u_site_soil_texture)),')'
         ENDIF
      ENDIF


      ! (9) depth to bedrock
      IF (DEF_USE_BEDROCK) THEN
         readflag = ((.not. mksrfdata) .or. USE_SITE_dbedrock)
         u_site_dbedrock = readflag &
            .and. ncio_var_exist (fsrfdata, 'depth_to_bedrock',readflag)
         IF (u_site_dbedrock) THEN
            CALL ncio_read_serial (fsrfdata, 'depth_to_bedrock', SITE_dbedrock)
         ELSE
            CALL gridrock%define_by_name ('colm_500m')
            filename = trim(DEF_dir_rawdata)//'/bedrock.nc'
            CALL read_point_var_2d_real8 (gridrock, filename, 'dbedrock', &
               SITE_lon_location, SITE_lat_location, SITE_dbedrock)
         ENDIF

         IF (mksrfdata) THEN
            write(*,'(A,F8.2,3A)') 'Depth to bedrock : ', SITE_dbedrock, ' (from ',trim(datasource(u_site_dbedrock)),')'
         ENDIF

      ENDIF

      ! (10) topography
      readflag = ((.not. mksrfdata) .or. USE_SITE_topography)

      u_site_elevation = readflag &
         .and. ncio_var_exist (fsrfdata, 'elevation',readflag)
      IF (u_site_elevation) THEN
         CALL ncio_read_serial (fsrfdata, 'elevation', SITE_elevation)
      ELSE
         CALL gridtopo%define_by_name ('colm_500m')
         filename = trim(DEF_dir_rawdata)//'/topography.nc'
         CALL read_point_var_2d_real8 (gridtopo, filename, 'elevation', &
            SITE_lon_location, SITE_lat_location, SITE_elevation)
      ENDIF

      u_site_elvstd = readflag &
         .and. ncio_var_exist (fsrfdata, 'elvstd',readflag)
      IF (u_site_elvstd) THEN
         CALL ncio_read_serial (fsrfdata, 'elvstd', SITE_elvstd)
      ELSE
         CALL gridtopo%define_by_name ('colm_500m')
         filename = trim(DEF_dir_rawdata)//'/topography.nc'
         CALL read_point_var_2d_real8 (gridtopo, filename, 'elvstd', &
            SITE_lon_location, SITE_lat_location, SITE_elvstd)
      ENDIF

      u_site_sloperatio = readflag &
         .and. ncio_var_exist (fsrfdata, 'sloperatio',readflag)
      IF (u_site_sloperatio) THEN
         CALL ncio_read_serial (fsrfdata, 'sloperatio', SITE_sloperatio)
      ELSE
         CALL gridtopo%define_by_name ('colm_500m')
         filename = trim(DEF_dir_rawdata)//'/topography.nc'
         CALL read_point_var_2d_real8 (gridtopo, filename, 'slope', &
            SITE_lon_location, SITE_lat_location, SITE_sloperatio)
      ENDIF

      IF (DEF_USE_Forcing_Downscaling) THEN

         filename = trim(DEF_DS_HiresTopographyDataDir)//"/slope.nc"
         IF (ncio_var_exist(filename,'lat') .and. ncio_var_exist(filename,'lon')) THEN
            CALL grid_topo_factor%define_from_file (filename, "lat", "lon")
         ENDIF

         u_site_svf = readflag &
            .and. ncio_var_exist (fsrfdata, 'SITE_svf',readflag)
         IF (u_site_svf) THEN
            CALL ncio_read_serial (fsrfdata, 'SITE_svf' , SITE_svf)
         ELSE
            filename = trim(DEF_DS_HiresTopographyDataDir)//"/sky_view_factor.nc"
            CALL read_point_var_2d_real8 (grid_topo_factor, filename, 'svf', &
               SITE_lon_location, SITE_lat_location, SITE_svf)
         ENDIF

         u_site_cur = readflag &
            .and. ncio_var_exist (fsrfdata, 'SITE_cur',readflag)
         IF (u_site_cur) THEN
            CALL ncio_read_serial (fsrfdata, 'SITE_cur' , SITE_cur)
         ELSE
            filename = trim(DEF_DS_HiresTopographyDataDir)//"/curvature.nc"
            CALL read_point_var_2d_real8 (grid_topo_factor, filename, 'curvature', &
               SITE_lon_location, SITE_lat_location, SITE_cur)
         ENDIF

         u_site_slp_type = readflag &
            .and. ncio_var_exist (fsrfdata, 'SITE_slp_type', readflag) &
            .and. ncio_var_exist (fsrfdata, 'SITE_asp_type', readflag) &
            .and. ncio_var_exist (fsrfdata, 'SITE_area_type',readflag)
         u_site_asp_type  = u_site_slp_type
         u_site_area_type = u_site_slp_type

         IF (u_site_slp_type) THEN
            CALL ncio_read_serial (fsrfdata, 'SITE_slp_type' , SITE_slp_type  )
            CALL ncio_read_serial (fsrfdata, 'SITE_asp_type' , SITE_asp_type  )
            CALL ncio_read_serial (fsrfdata, 'SITE_area_type', SITE_area_type )
         ELSE
            filename = trim(DEF_DS_HiresTopographyDataDir)//"/slope.nc"
            CALL read_point_var_2d_real8 (grid_topo_factor, filename, 'slope', &
               SITE_lon_location, SITE_lat_location, slp)

            filename = trim(DEF_DS_HiresTopographyDataDir)//"/aspect.nc"
            CALL read_point_var_2d_real8 (grid_topo_factor, filename, 'aspect', &
               SITE_lon_location, SITE_lat_location, asp)

            allocate (SITE_slp_type  (num_slope_type)); SITE_slp_type (:) = 0.
            allocate (SITE_asp_type  (num_slope_type)); SITE_asp_type (:) = 0.
            allocate (SITE_area_type (num_slope_type)); SITE_area_type(:) = 0.

            IF ((asp.ge.0 .and. asp.le.90*pi/180) .or. (asp.ge.270*pi/180 .and. asp.le.360*pi/180)) THEN
               IF ((slp.ge.15*pi/180)) THEN  ! north abrupt slope
                  typ = 1
               ELSE                          ! north gentle slope
                  typ = 2
               ENDIF
            ELSE
               IF ((slp.ge.15*pi/180)) THEN  ! south abrupt slope
                  typ = 3
               ELSE                          ! south gentle slope
                  typ = 4
               ENDIF
            ENDIF

            SITE_slp_type (typ) = slp
            SITE_asp_type (typ) = asp
            SITE_area_type(typ) = 1.

         ENDIF

         u_site_sf_lut = readflag &
            .and. ncio_var_exist (fsrfdata, 'SITE_sf_lut',readflag)
         IF (u_site_sf_lut) THEN
            CALL ncio_read_serial (fsrfdata, 'SITE_sf_lut', SITE_sf_lut)
         ELSE
            filename = trim(DEF_DS_HiresTopographyDataDir)//"/terrain_elev_angle_front.nc"
            CALL read_point_var_3d_first_real8 (grid_topo_factor, filename, 'tea_front', &
               SITE_lon_location, SITE_lat_location, num_azimuth, tea_f)

            filename = trim(DEF_DS_HiresTopographyDataDir)//"/terrain_elev_angle_back.nc"
            CALL read_point_var_3d_first_real8 (grid_topo_factor, filename, 'tea_back', &
               SITE_lon_location, SITE_lat_location, num_azimuth, tea_b)

            allocate (SITE_sf_lut (num_azimuth, num_zenith))

            DO a = 1, num_azimuth

               tea_f(a) = asin(max(min(tea_f(a),1.),-1.))
               tea_b(a) = asin(max(min(tea_b(a),1.),-1.))

               IF (tea_f(a) <= tea_b(a)) tea_f(a) = tea_b(a) + 0.001

               DO z = 1, num_zenith
                  zenith_angle = pi/(2*num_zenith)*(z-1)

                  IF (pi*0.5 - zenith_angle < tea_b(a)) THEN
                     SITE_sf_lut(a,z) = 0
                  ELSE IF (pi*0.5 - zenith_angle > tea_f(a)) THEN
                     SITE_sf_lut(a,z) = 1
                  ELSE
                     SITE_sf_lut(a,z) = (0.5*pi - zenith_angle - tea_b(a))/(tea_f(a) - tea_b(a))
                  ENDIF

               ENDDO
            ENDDO

         ENDIF
      ENDIF

      IF (mksrfdata) THEN
         write(*,'(A,F8.2,3A)') 'Elevation : ', SITE_elevation , ' (from ',trim(datasource(u_site_elevation)),')'
         write(*,'(A,F8.2,3A)') 'Elv std   : ', SITE_elvstd    , ' (from ',trim(datasource(u_site_elvstd)),')'
         write(*,'(A,F8.2,3A)') 'SlopeRatio: ', SITE_sloperatio, ' (from ',trim(datasource(u_site_sloperatio)),')'

         IF (DEF_USE_Forcing_Downscaling) THEN
            write(*,'(A,F8.2,3A)') 'Sky view factor : ', SITE_svf, ' (from ',trim(datasource(u_site_svf)),')'
            write(*,'(A,F8.2,3A)') 'Curvature       : ', SITE_cur, ' (from ',trim(datasource(u_site_cur)),')'
            write(c,'(I0)') num_slope_type
            write(*,'(A,'//trim(c)//'F8.2,3A)') 'Slope  type     : ', SITE_slp_type,  ' (from ',trim(datasource(u_site_slp_type)),')'
            write(*,'(A,'//trim(c)//'F8.2,3A)') 'Aspect type     : ', SITE_asp_type,  ' (from ',trim(datasource(u_site_slp_type)),')'
            write(*,'(A,'//trim(c)//'F8.2,3A)') 'Slope type area : ', SITE_area_type, ' (from ',trim(datasource(u_site_slp_type)),')'
            write(c,'(I0)') num_azimuth*num_zenith
            write(*,'(A,A,I3,A,I3,A,'//trim(c)//'F8.2,3A)') 'Shadow lookup table    : ', &
               '(', num_azimuth, ' in azimuth,', num_zenith, ' in zenith)', &
               SITE_sf_lut , ' (from ',trim(datasource(u_site_sf_lut)),')'
         ENDIF
      ENDIF


      IF (.not. mksrfdata) THEN

         landpatch%nset = numpatch

         allocate (landpatch%settyp (numpatch)); landpatch%settyp = SITE_landtype

         landpatch%nblkgrp = 1
         allocate (landpatch%xblkgrp(1));       landpatch%xblkgrp(1) = 1
         allocate (landpatch%yblkgrp(1));       landpatch%yblkgrp(1) = 1

         allocate (landpatch%vecgs%vlen(1,1));  landpatch%vecgs%vlen(1,1) = numpatch
         allocate (landpatch%vecgs%vstt(1,1));  landpatch%vecgs%vstt(1,1) = 1
         allocate (landpatch%vecgs%vend(1,1));  landpatch%vecgs%vend(1,1) = numpatch

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)

         IF (numpft > 0) THEN

            landpft%nset = numpft

            allocate (landpft%settyp (numpft)); landpft%settyp = SITE_pfttyp

            landpft%nblkgrp = 1
            allocate (landpft%xblkgrp(1));       landpft%xblkgrp(1) = 1
            allocate (landpft%yblkgrp(1));       landpft%yblkgrp(1) = 1

            allocate (landpft%vecgs%vlen(1,1));  landpft%vecgs%vlen(1,1) = numpft
            allocate (landpft%vecgs%vstt(1,1));  landpft%vecgs%vstt(1,1) = 1
            allocate (landpft%vecgs%vend(1,1));  landpft%vecgs%vend(1,1) = numpft

            allocate (patch_pft_s (numpatch))
            allocate (patch_pft_e (numpatch))
            allocate (pft2patch   (numpft  ))
#ifdef CROP
            IF (SITE_landtype == CROPLAND) THEN
               patch_pft_s = (/(i, i = 1, numpatch)/)
               patch_pft_e = (/(i, i = 1, numpatch)/)
               pft2patch   = (/(i, i = 1, numpatch)/)
            ELSE
               patch_pft_s = 1
               patch_pft_e = numpft
               pft2patch   = 1
            ENDIF
#else
            patch_pft_s = 1
            patch_pft_e = numpft
            pft2patch   = 1
#endif
         ENDIF
#endif

         numelm = 1
         allocate (landelm%settyp   (1));  landelm%settyp  (1) = 0
         allocate (elm_patch%substt (1));  elm_patch%substt(1) = 1
         allocate (elm_patch%subend (1));  elm_patch%subend(1) = numpatch

         allocate (elm_patch%subfrc (numpatch)); elm_patch%subfrc = 1./numpatch
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
#ifdef CROP
         elm_patch%subfrc = SITE_pctcrop
#endif
#endif

      ENDIF

   END SUBROUTINE read_surface_data_single

!-----------------------------------------------------------------------
   SUBROUTINE read_urban_surface_data_single (fsrfdata, mksrfdata, mkrun)

   USE MOD_TimeManager
   USE MOD_Grid
   USE MOD_Block
   USE MOD_NetCDFSerial
   USE MOD_NetCDFPoint
   USE MOD_Namelist
   USE MOD_Utils
   USE MOD_SPMD_Task
   USE MOD_LandPatch
   USE MOD_LandUrban
   USE MOD_Urban_Const_LCZ
   USE MOD_Vars_Global, only: PI, URBAN
   USE MOD_Mesh,        only: numelm
   USE MOD_LandElm
   IMPLICIT NONE

   character(len=*), intent(in) :: fsrfdata
   logical, intent(in) :: mksrfdata
   logical, intent(in), optional :: mkrun

   ! Local Variables
   real(r8), allocatable, dimension(:,:)     :: hwrbld_ncar , fgper_ncar  , htroof_ncar , wtroof_ncar
   real(r8), allocatable, dimension(:,:)     :: emroof_ncar , emwall_ncar , emgimp_ncar , emgper_ncar
   real(r8), allocatable, dimension(:,:)     :: thkroof_ncar, thkwall_ncar, tbldmin_ncar, tbldmax_ncar
   real(r8), allocatable, dimension(:,:,:)   :: cvroof_ncar , cvwall_ncar , cvgimp_ncar , &
                                                tkroof_ncar , tkwall_ncar , tkgimp_ncar
   real(r8), allocatable, dimension(:,:,:,:) :: albroof_ncar, albwall_ncar, albgimp_ncar, albgper_ncar

   real(r8) :: lat_in, lon_in
   real(r8) :: LAI, lakedepth, slp, asp, zenith_angle
   integer  :: i, isc, nsl, typ, a, z, rid, utyp
   integer  :: iyear, idate(3), simulation_lai_year_start, simulation_lai_year_end
   integer  :: start_year, end_year, ntime, itime, pop_i
   logical  :: readflag

   character(len=256) :: filename, dir_5x5
   character(len=4)   :: cyear, c, c5year

   type(grid_type) :: gridupatch, gridhroof, gridfroof , gridhtopu, gridfvegu, gridflakeu, gridlaiu, &
                      gridpopu  , gridlucy , gridbright, gridsoil , gridrock , gridtopo  , gridlake, &
                      grid_topo_factor

   real(r8), allocatable :: tea_f(:), tea_b(:)

      u_site_froof = .false.; u_site_hroof = .false.; u_site_fgper   = .false.; u_site_hlr     = .false.;
      u_site_fveg  = .false.; u_site_htopu = .false.; u_site_urblai  = .false.; u_site_urbsai  = .false.;
      u_site_flake = .false.;

      u_site_albr  = .false.; u_site_albw  = .false.; u_site_albgimp = .false.; u_site_albgper = .false.;
      u_site_emr   = .false.; u_site_emw   = .false.; u_site_emgimp  = .false.; u_site_emgper  = .false.;

      u_site_cvr   = .false.; u_site_cvw   = .false.; u_site_cvgimp  = .false.;
      u_site_tkr   = .false.; u_site_tkw   = .false.; u_site_tkgimp  = .false.;

      u_site_tbmax = .false.; u_site_tbmin = .false.; u_site_thickr  = .false.; u_site_thickw  = .false.;
      u_site_pop   = .false.; u_site_lucy  = .false.;

      IF (mksrfdata) THEN
         write(*,*)
         write(*,*) '  ----------------  Make Single Point Surface Data  ----------------  '
      ENDIF

      IF (ncio_var_exist(fsrfdata, 'latitude')) THEN
         CALL ncio_read_serial (fsrfdata, 'latitude',  lat_in)
         IF (lat_in /= SITE_lat_location) THEN
            write(*,*) 'Warning: Latitude mismatch: ', &
               lat_in, ' in data file and ', SITE_lat_location, 'in namelist.'
         ENDIF
         SITE_lat_location = lat_in
      ENDIF

      IF (ncio_var_exist(fsrfdata, 'longitude')) THEN
         CALL ncio_read_serial (fsrfdata, 'longitude', lon_in)
         IF (lon_in /= SITE_lon_location) THEN
            write(*,*) 'Warning: Longitude mismatch: ', &
               lon_in, ' in data file and ', SITE_lon_location, 'in namelist.'
         ENDIF
         SITE_lon_location = lon_in
      ENDIF

      CALL normalize_longitude (SITE_lon_location)

      IF (.not. isgreenwich) THEN
         LocalLongitude = SITE_lon_location
      ENDIF

      IF (mksrfdata) THEN
         write(*,'(A,F8.2)') 'Latitude  : ', SITE_lat_location
         write(*,'(A,F8.2)') 'Longitude : ', SITE_lon_location
      ENDIF

      DEF_domain%edges = floor(SITE_lat_location)
      DEF_domain%edgen = floor(SITE_lat_location) + 1.
      DEF_domain%edgew = floor(SITE_lon_location)
      DEF_domain%edgee = floor(SITE_lon_location) + 1.

      CALL gblock%set ()
      gblock%nblkme = 1
      allocate(gblock%xblkme(1))
      allocate(gblock%yblkme(1))
      gblock%xblkme(1) = find_nearest_west  (SITE_lon_location, gblock%nxblk, gblock%lon_w)
      gblock%yblkme(1) = find_nearest_south (SITE_lat_location, gblock%nyblk, gblock%lat_s)

      ! (1) build/read "land patch" by using land cover type data
      IF (trim(fsrfdata) /= 'null') THEN
         SITE_landtype = URBAN
      ELSEIF (SITE_landtype /= URBAN) THEN
         write(*,*) 'Error! Please set SITE_landtype to URBAN in namelist file !'
         CALL CoLM_stop()
      ENDIF

      numpatch = 1
      numurban = 1
      u_site_landtype = (SITE_landtype >= 0)
      IF (mksrfdata) THEN
         write(*,'(A,A,3A)') 'Land cover type : ', trim(patchclassname(SITE_landtype)), &
            ' (from ',trim(datasource(u_site_landtype)),')'
      ENDIF

      IF (mksrfdata) THEN
         ! (2) build/read "urban type" by using land cover type data
IF (DEF_URBAN_type_scheme == 1) THEN
         u_site_utype = ncio_var_exist(fsrfdata,'URBAN_DENSITY_CLASS')

         IF ( u_site_utype ) THEN
            CALL ncio_read_serial (fsrfdata, 'URBTYP'             , SITE_ncar_rid)
            CALL ncio_read_serial (fsrfdata, 'URBAN_DENSITY_CLASS', SITE_urbtyp  )
         ELSE
            CALL gridupatch%define_by_name ('colm_500m')

            dir_5x5 = trim(DEF_dir_rawdata) // '/urban_type/'
            CALL read_point_5x5_var_2d_int32 (gridupatch, dir_5x5, 'URBTYP', 'REGION_ID', &
               SITE_lon_location, SITE_lat_location, SITE_ncar_rid)

            CALL read_point_5x5_var_2d_int32 (gridupatch, dir_5x5, 'URBTYP', 'URBAN_DENSITY_CLASS', &
               SITE_lon_location, SITE_lat_location, SITE_urbtyp)

            write(*,'(A,I0,A,I0,3A)') 'Urban type : NCAR ', SITE_urbtyp, ' of Region ', SITE_ncar_rid, &
               ' (from ',trim(datasource(u_site_utype)),')'
         ENDIF
ELSE
         u_site_utype = ncio_var_exist(fsrfdata,'LCZ_DOM')
         IF ( u_site_utype ) THEN
            CALL ncio_read_serial (fsrfdata, 'LCZ_DOM', SITE_urbtyp  )
         ELSE
            CALL gridupatch%define_by_name ('colm_500m')

            dir_5x5 = trim(DEF_dir_rawdata) // '/urban_type/'
            CALL read_point_5x5_var_2d_int32 (gridupatch, dir_5x5, 'URBTYP', 'LCZ_DOM', &
               SITE_lon_location, SITE_lat_location, SITE_urbtyp)
         ENDIF
         write(*,'(A,I0,3A)') 'Urban type : LCZ ', SITE_urbtyp, &
            ' (from ',trim(datasource(u_site_utype)),')'
ENDIF


         ! (4) urban geometry
         readflag     = USE_SITE_urban_geometry
         u_site_hroof = ncio_var_exist(fsrfdata,'building_mean_height',readflag)
         IF ( u_site_hroof ) THEN
            CALL ncio_read_serial (fsrfdata, 'building_mean_height', SITE_hroof  )
         ELSE
            CALL gridhroof%define_by_name ('colm_500m')
            dir_5x5 = trim(DEF_dir_rawdata) // '/urban/'
            write(c5year, '(i4.4)') int(DEF_LC_YEAR/5)*5

IF (DEF_Urban_geom_data == 1) THEN
            CALL read_point_5x5_var_2d_real8 (gridhroof, dir_5x5, 'URBSRF'//trim(c5year), 'HT_ROOF_GHSL', &
               SITE_lon_location, SITE_lat_location, SITE_hroof)
ELSE
            CALL read_point_5x5_var_2d_real8 (gridhroof, dir_5x5, 'URBSRF'//trim(c5year), 'HT_ROOF_Li', &
               SITE_lon_location, SITE_lat_location, SITE_hroof)
ENDIF
         ENDIF

         u_site_froof = readflag .and. ncio_var_exist(fsrfdata,'roof_area_fraction',readflag)
         IF ( u_site_froof ) THEN
            CALL ncio_read_serial (fsrfdata, 'roof_area_fraction', SITE_froof  )
         ELSE
            CALL gridfroof%define_by_name ('colm_500m')
            dir_5x5 = trim(DEF_dir_rawdata) // '/urban/'
            write(c5year, '(i4.4)') int(DEF_LC_YEAR/5)*5
IF (DEF_Urban_geom_data == 1) THEN
            CALL read_point_5x5_var_2d_real8 (gridfroof, dir_5x5, 'URBSRF'//trim(c5year), 'PCT_ROOF_GHSL', &
               SITE_lon_location, SITE_lat_location, SITE_froof)
ELSE
            CALL read_point_5x5_var_2d_real8 (gridfroof, dir_5x5, 'URBSRF'//trim(c5year), 'PCT_ROOF_Li', &
               SITE_lon_location, SITE_lat_location, SITE_froof)
ENDIF
         ENDIF

         u_site_fgper  = readflag .and. ncio_var_exist(fsrfdata,'impervious_area_fraction',readflag)
         IF ( u_site_fgper ) THEN
            CALL ncio_read_serial (fsrfdata, 'impervious_area_fraction', SITE_fgimp  )
         ENDIF

         u_site_thickr  = readflag .and. ncio_var_exist(fsrfdata,'THICK_ROOF',readflag)
         IF ( u_site_thickr ) THEN
            CALL ncio_read_serial (fsrfdata, 'THICK_ROOF', SITE_thickroof  )
         ENDIF

         u_site_thickw  = readflag .and. ncio_var_exist(fsrfdata,'THICK_WALL',readflag)
         IF ( u_site_thickw ) THEN
            CALL ncio_read_serial (fsrfdata, 'THICK_WALL', SITE_thickwall  )
         ENDIF

IF (DEF_USE_CANYON_HWR) THEN
         u_site_hlr  = readflag .and. ncio_var_exist(fsrfdata,'canyon_height_width_ratio',readflag)
         IF ( u_site_hlr ) THEN
            CALL ncio_read_serial (fsrfdata, 'canyon_height_width_ratio', SITE_hlr  )
         ENDIF
ELSE
         u_site_hlr  = readflag .and. ncio_var_exist(fsrfdata,'wall_to_plan_area_ratio',readflag)
         IF ( u_site_hlr ) THEN
            CALL ncio_read_serial (fsrfdata, 'wall_to_plan_area_ratio', SITE_lambdaw  )
            SITE_hlr     = SITE_lambdaw/4/SITE_froof
         ENDIF
ENDIF
         ! (5) urban ecology
         readflag      = USE_SITE_urban_ecology
         u_site_htopu  = readflag .and. ncio_var_exist(fsrfdata,'tree_mean_height',readflag)
         IF ( u_site_htopu ) THEN
            CALL ncio_read_serial (fsrfdata, 'tree_mean_height', SITE_htop_urb  )
         ELSE
            CALL gridhtopu%define_by_name ('colm_500m')
            dir_5x5 = trim(DEF_dir_rawdata) // '/urban/'
            write(c5year, '(i4.4)') int(DEF_LC_YEAR/5)*5

            CALL read_point_5x5_var_2d_real8 (gridhtopu, dir_5x5, 'URBSRF'//trim(c5year), 'HTOP', &
               SITE_lon_location, SITE_lat_location, SITE_htop_urb)
         ENDIF

         u_site_flake  = readflag .and. ncio_var_exist(fsrfdata,'water_area_fraction',readflag)
         IF ( u_site_flake ) THEN
            CALL ncio_read_serial (fsrfdata, 'water_area_fraction', SITE_flake_urb  )
         ELSE
            CALL gridflakeu%define_by_name ('colm_500m')
            dir_5x5 = trim(DEF_dir_rawdata) // '/urban/'
            write(c5year, '(i4.4)') int(DEF_LC_YEAR/5)*5

            CALL read_point_5x5_var_2d_real8 (gridflakeu, dir_5x5, 'URBSRF'//trim(c5year), 'PCT_Water', &
               SITE_lon_location, SITE_lat_location, SITE_flake_urb)

            SITE_flake_urb = SITE_flake_urb/100
         ENDIF

         u_site_fveg = readflag .and. ncio_var_exist(fsrfdata,'tree_area_fraction',readflag)
         IF ( u_site_fveg ) THEN
            CALL ncio_read_serial (fsrfdata, 'tree_area_fraction', SITE_fveg_urb  )
         ELSE
            CALL gridfvegu%define_by_name ('colm_500m')
            dir_5x5 = trim(DEF_dir_rawdata) // '/urban/'
            write(c5year, '(i4.4)') int(DEF_LC_YEAR/5)*5

            CALL read_point_5x5_var_2d_real8 (gridfvegu, dir_5x5, 'URBSRF'//trim(c5year), 'PCT_Tree', &
               SITE_lon_location, SITE_lat_location, SITE_fveg_urb)

            SITE_fveg_urb = SITE_fveg_urb/100
         ENDIF

         u_site_urblai = readflag .and. ncio_var_exist(fsrfdata,'TREE_LAI',readflag)
         IF ( u_site_urblai) THEN
            CALL ncio_read_serial (fsrfdata, 'TREE_LAI', SITE_LAI_monthly  )
            CALL ncio_read_serial (fsrfdata, 'TREE_SAI', SITE_SAI_monthly  )
         ELSE

            idate(1) = DEF_simulation_time%start_year
            IF (.not. isgreenwich) THEN
               idate(3) = DEF_simulation_time%start_sec
               CALL monthday2julian (idate(1), &
                  DEF_simulation_time%start_month, DEF_simulation_time%start_day, idate(2))
               CALL localtime2gmt(idate)
            ENDIF
            simulation_lai_year_start = idate(1)

            idate(1) = DEF_simulation_time%end_year
            IF (.not. isgreenwich) THEN
               idate(3) = DEF_simulation_time%end_sec
               CALL monthday2julian (idate(1), &
                  DEF_simulation_time%end_month, DEF_simulation_time%end_day, idate(2))
               CALL localtime2gmt(idate)
            ENDIF
            simulation_lai_year_end = idate(1)

            IF (DEF_LAI_CHANGE_YEARLY) THEN
               start_year = max(simulation_lai_year_start, DEF_LAI_START_YEAR )
               start_year = min(start_year,                DEF_LAI_END_YEAR   )
               end_year   = min(simulation_lai_year_end,   DEF_LAI_END_YEAR   )
               end_year   = max(end_year,                  DEF_LAI_START_YEAR )
            ELSE
               start_year = DEF_LC_YEAR
               end_year   = DEF_LC_YEAR
            ENDIF

            allocate (SITE_LAI_year (start_year:end_year))
            SITE_LAI_year = (/(iyear, iyear = start_year, end_year)/)

            ntime = 12
            allocate (SITE_LAI_monthly (12,start_year:end_year))
            allocate (SITE_SAI_monthly (12,start_year:end_year))

            CALL gridlaiu%define_by_name ('colm_500m')
            dir_5x5 = trim(DEF_dir_rawdata) // '/urban_lai_500m/'

            DO iyear = start_year, end_year
               write(cyear,'(i4.4)') iyear
               DO itime = 1, ntime
                  CALL read_point_5x5_var_2d_time_real8 (gridlaiu, dir_5x5, 'URBLAI_'//trim(cyear), &
                        'URBAN_TREE_LAI', SITE_lon_location, SITE_lat_location, itime, &
                        SITE_LAI_monthly(itime,iyear))

                  CALL read_point_5x5_var_2d_time_real8 (gridlaiu, dir_5x5, 'URBLAI_'//trim(cyear), &
                        'URBAN_TREE_SAI', SITE_lon_location, SITE_lat_location, itime, &
                        SITE_SAI_monthly(itime,iyear))
               ENDDO
            ENDDO
         ENDIF

         ! (6) urban radiation
         readflag    = USE_SITE_urban_radiation
         u_site_albr = readflag .and. ncio_var_exist(fsrfdata,'ALB_ROOF',readflag)
         IF ( u_site_albr ) THEN
            CALL ncio_read_serial (fsrfdata, 'ALB_ROOF', SITE_alb_roof  )
         ENDIF

         u_site_albw = readflag .and. ncio_var_exist(fsrfdata,'ALB_WALL',readflag)
         IF ( u_site_albw ) THEN
            CALL ncio_read_serial (fsrfdata, 'ALB_WALL', SITE_alb_wall  )
         ENDIF

         u_site_albgper = readflag .and. ncio_var_exist(fsrfdata,'ALB_GPER',readflag)
         IF ( u_site_albgper ) THEN
            CALL ncio_read_serial (fsrfdata, 'ALB_GPER', SITE_alb_gper  )
         ENDIF

         u_site_albgimp = readflag .and. ncio_var_exist(fsrfdata,'ALB_GIMP',readflag)
         IF ( u_site_albgimp ) THEN
            CALL ncio_read_serial (fsrfdata, 'ALB_GIMP', SITE_alb_gimp  )
         ENDIF

         u_site_emr = readflag .and. ncio_var_exist(fsrfdata,'EM_ROOF',readflag)
         IF ( u_site_emr ) THEN
            CALL ncio_read_serial (fsrfdata, 'EM_ROOF', SITE_em_roof  )
         ENDIF

         u_site_emw = readflag .and. ncio_var_exist(fsrfdata,'EM_WALL',readflag)
         IF ( u_site_emw ) THEN
            CALL ncio_read_serial (fsrfdata, 'EM_WALL', SITE_em_wall  )
         ENDIF

         u_site_emgper = readflag .and. ncio_var_exist(fsrfdata,'EM_GPER',readflag)
         IF ( u_site_emgper ) THEN
            CALL ncio_read_serial (fsrfdata, 'EM_GPER', SITE_em_gper  )
         ENDIF

         u_site_emgimp = readflag .and. ncio_var_exist(fsrfdata,'EM_GIMP',readflag)
         IF ( u_site_emgimp ) THEN
            CALL ncio_read_serial (fsrfdata, 'EM_GIMP', SITE_em_gimp  )
         ENDIF

         ! (6) urban thermal
         readflag   = USE_SITE_urban_thermal
         u_site_cvr = readflag .and. ncio_var_exist(fsrfdata,'CV_ROOF',readflag)
         IF ( u_site_cvr ) THEN
            CALL ncio_read_serial (fsrfdata, 'CV_ROOF', SITE_cv_roof  )
         ENDIF

         u_site_cvw = readflag .and. ncio_var_exist(fsrfdata,'CV_WALL',readflag)
         IF ( u_site_cvw ) THEN
            CALL ncio_read_serial (fsrfdata, 'CV_WALL', SITE_cv_wall  )
         ENDIF

         u_site_cvgimp = readflag .and. ncio_var_exist(fsrfdata,'CV_GIMP',readflag)
         IF ( u_site_cvgimp ) THEN
            CALL ncio_read_serial (fsrfdata, 'CV_GIMP', SITE_cv_gimp  )
         ENDIF

         u_site_tkr = readflag .and. ncio_var_exist(fsrfdata,'TK_ROOF',readflag)
         IF ( u_site_tkr ) THEN
            CALL ncio_read_serial (fsrfdata, 'TK_ROOF', SITE_tk_roof  )
         ENDIF

         u_site_tkw = readflag .and. ncio_var_exist(fsrfdata,'TK_WALL',readflag)
         IF ( u_site_tkw ) THEN
            CALL ncio_read_serial (fsrfdata, 'TK_WALL', SITE_tk_wall  )
         ENDIF

         u_site_tkgimp = readflag .and. ncio_var_exist(fsrfdata,'TK_GIMP',readflag)
         IF ( u_site_tkgimp ) THEN
            CALL ncio_read_serial (fsrfdata, 'TK_GIMP', SITE_tk_gimp  )
         ENDIF

         ! (6) urban human
         readflag  = USE_SITE_urban_human
         u_site_pop= readflag .and. ncio_var_exist(fsrfdata,'resident_population_density',readflag)
         IF ( u_site_pop) THEN
            CALL ncio_read_serial (fsrfdata, 'resident_population_density', SITE_popden  )
         ELSE
            CALL gridpopu%define_by_name ('colm_500m')
            dir_5x5 = trim(DEF_dir_rawdata) // '/urban/'
            write(c5year, '(i4.4)') int(DEF_LC_YEAR/5)*5

            IF (mod(DEF_LC_YEAR,5) == 0) THEN
               pop_i = 1
            ELSE
               pop_i = 5 - (ceiling(DEF_LC_YEAR*1./5.)*5 - DEF_LC_YEAR) + 1
            ENDIF

            CALL read_point_5x5_var_2d_time_real8 (gridlaiu, dir_5x5, 'URBSRF'//trim(c5year), &
                  'POP_DEN', SITE_lon_location, SITE_lat_location, pop_i, &
                  SITE_popden)
         ENDIF

         u_site_lucy= readflag .and. ncio_var_exist(fsrfdata,'LUCY_ID',readflag)
         IF ( u_site_lucy) THEN
            CALL ncio_read_serial (fsrfdata, 'LUCY_ID', SITE_lucyid )
         ELSE
            CALL gridlucy%define_by_name ('colm_5km')
            filename = trim(DEF_dir_rawdata) // '/urban/LUCY_regionid.nc'

            CALL read_point_var_2d_real8 (gridlucy, filename, 'LUCY_REGION_ID', &
               SITE_lon_location, SITE_lat_location, SITE_lucyid)
         ENDIF

         u_site_tbmax= readflag .and. ncio_var_exist(fsrfdata,'T_BUILDING_MAX',readflag)
         IF ( u_site_tbmax) THEN
            CALL ncio_read_serial (fsrfdata, 'T_BUILDING_MAX', SITE_t_roommax  )
         ENDIF

         u_site_tbmin= readflag .and. ncio_var_exist(fsrfdata,'T_BUILDING_MIN',readflag)
         IF ( u_site_tbmin) THEN
            CALL ncio_read_serial (fsrfdata, 'T_BUILDING_MIN', SITE_t_roommin  )
         ENDIF

IF (DEF_URBAN_type_scheme == 1) THEN
         filename = trim(DEF_dir_rawdata)//'urban/NCAR_urban_properties.nc'

         CALL ncio_read_bcast_serial (filename,  "WTLUNIT_ROOF"  , wtroof_ncar )
         CALL ncio_read_bcast_serial (filename,  "HT_ROOF"       , htroof_ncar )
         CALL ncio_read_bcast_serial (filename,  "CANYON_HWR"    , hwrbld_ncar )
         CALL ncio_read_bcast_serial (filename,  "WTROAD_PERV"   , fgper_ncar  )
         CALL ncio_read_bcast_serial (filename,  "EM_ROOF"       , emroof_ncar )
         CALL ncio_read_bcast_serial (filename,  "EM_WALL"       , emwall_ncar )
         CALL ncio_read_bcast_serial (filename,  "EM_IMPROAD"    , emgimp_ncar )
         CALL ncio_read_bcast_serial (filename,  "EM_PERROAD"    , emgper_ncar )
         CALL ncio_read_bcast_serial (filename,  "ALB_ROOF"      , albroof_ncar)
         CALL ncio_read_bcast_serial (filename,  "ALB_WALL"      , albwall_ncar)
         CALL ncio_read_bcast_serial (filename,  "ALB_IMPROAD"   , albgimp_ncar)
         CALL ncio_read_bcast_serial (filename,  "ALB_PERROAD"   , albgper_ncar)
         CALL ncio_read_bcast_serial (filename,  "TK_ROOF"       , tkroof_ncar )
         CALL ncio_read_bcast_serial (filename,  "TK_WALL"       , tkwall_ncar )
         CALL ncio_read_bcast_serial (filename,  "TK_IMPROAD"    , tkgimp_ncar )
         CALL ncio_read_bcast_serial (filename,  "CV_ROOF"       , cvroof_ncar )
         CALL ncio_read_bcast_serial (filename,  "CV_WALL"       , cvwall_ncar )
         CALL ncio_read_bcast_serial (filename,  "CV_IMPROAD"    , cvgimp_ncar )
         CALL ncio_read_bcast_serial (filename,  "THICK_ROOF"    , thkroof_ncar)
         CALL ncio_read_bcast_serial (filename,  "THICK_WALL"    , thkwall_ncar)
         CALL ncio_read_bcast_serial (filename,  "T_BUILDING_MIN", tbldmin_ncar)
         CALL ncio_read_bcast_serial (filename,  "T_BUILDING_MAX", tbldmax_ncar)

         rid  = SITE_ncar_rid
         utyp = SITE_urbtyp

         IF (.not. u_site_emr    ) SITE_em_roof   = emroof_ncar(utyp, rid)
         IF (.not. u_site_emw    ) SITE_em_wall   = emwall_ncar(utyp, rid)
         IF (.not. u_site_emgimp ) SITE_em_gimp   = emgimp_ncar(utyp, rid)
         IF (.not. u_site_emgper ) SITE_em_gper   = emgper_ncar(utyp, rid)

         IF (.not. u_site_tbmax  ) SITE_t_roommax = tbldmax_ncar (utyp, rid)
         IF (.not. u_site_tbmin  ) SITE_t_roommin = tbldmin_ncar (utyp, rid)

         IF (.not. u_site_thickr ) SITE_thickroof = thkroof_ncar(utyp, rid)
         IF (.not. u_site_thickw ) SITE_thickwall = thkwall_ncar(utyp, rid)

         IF (.not. u_site_cvr   ) THEN
            allocate( SITE_cv_roof (nl_roof) )
            SITE_cv_roof(:) = cvroof_ncar(utyp, rid, :)
         ENDIF

         IF (.not. u_site_cvw   ) THEN
            allocate( SITE_cv_wall (nl_wall) )
            SITE_cv_wall(:) = cvwall_ncar(utyp, rid, :)
         ENDIF

         IF (.not. u_site_cvgimp) THEN
            allocate( SITE_cv_gimp (nl_soil) )
            SITE_cv_gimp(:) = cvgimp_ncar(utyp, rid, :)
         ENDIF

         IF (.not. u_site_tkr   ) THEN
            allocate( SITE_tk_roof (nl_roof) )
            SITE_tk_roof(:) = tkroof_ncar(utyp, rid, :)
         ENDIF

         IF (.not. u_site_tkw   ) THEN
            allocate( SITE_tk_wall (nl_wall) )
            SITE_tk_wall(:) = tkwall_ncar(utyp, rid, :)
         ENDIF

         IF (.not. u_site_tkgimp) THEN
            allocate( SITE_tk_gimp (nl_soil) )
            SITE_tk_gimp(:) = tkgimp_ncar(utyp, rid, :)
         ENDIF

         IF (.not. u_site_albr   ) THEN
            allocate( SITE_alb_roof (2, 2) )
            SITE_alb_roof(:,:) = albroof_ncar(utyp, rid, :, :)
         ENDIF

         IF (.not. u_site_albw   ) THEN
            allocate( SITE_alb_wall (2, 2) )
            SITE_alb_wall(:,:) = albwall_ncar(utyp, rid, :, :)
         ENDIF

         IF (.not. u_site_albgimp) THEN
            allocate( SITE_alb_gimp (2, 2) )
            SITE_alb_gimp(:,:) = albgimp_ncar(utyp, rid, :, :)
         ENDIF

         IF (.not. u_site_albgper) THEN
            allocate( SITE_alb_gper (2, 2) )
            SITE_alb_gper(:,:) = albgper_ncar(utyp, rid, :, :)
         ENDIF

         IF (.not. u_site_hlr  ) SITE_hlr   = hwrbld_ncar(utyp, rid)
         IF (.not. u_site_fgper) SITE_fgimp = 1-fgper_ncar(utyp, rid)
ELSE
         utyp = SITE_urbtyp
         IF (.not. u_site_emr   ) SITE_em_roof = emroof_lcz(utyp)
         IF (.not. u_site_emw   ) SITE_em_wall = emwall_lcz(utyp)
         IF (.not. u_site_emgimp) SITE_em_gimp = emgimp_lcz(utyp)
         IF (.not. u_site_emgper) SITE_em_gper = emgper_lcz(utyp)

         IF (.not. u_site_tbmax) SITE_t_roommax = tbldmax_lcz(utyp)
         IF (.not. u_site_tbmin) SITE_t_roommin = tbldmin_lcz(utyp)

         IF (.not. u_site_thickr) SITE_thickroof = thkroof_lcz(utyp)
         IF (.not. u_site_thickw) SITE_thickwall = thkwall_lcz(utyp)

         IF (.not. u_site_cvr   ) THEN
            allocate( SITE_cv_roof (nl_roof) )
            SITE_cv_roof(:) = cvroof_lcz(utyp)
         ENDIF

         IF (.not. u_site_cvw   ) THEN
            allocate( SITE_cv_wall (nl_wall) )
            SITE_cv_wall(:) = cvwall_lcz(utyp)
         ENDIF

         IF (.not. u_site_cvgimp) THEN
            allocate( SITE_cv_gimp (nl_soil) )
            SITE_cv_gimp(:) = cvgimp_lcz(utyp)
         ENDIF

         IF (.not. u_site_tkr   ) THEN
            allocate( SITE_tk_roof (nl_roof) )
            SITE_tk_roof(:) = tkroof_lcz(utyp)
         ENDIF

         IF (.not. u_site_tkw   ) THEN
            allocate( SITE_tk_wall (nl_wall) )
            SITE_tk_wall(:) = tkwall_lcz(utyp)
         ENDIF

         IF (.not. u_site_tkgimp) THEN
            allocate( SITE_tk_gimp (nl_soil) )
            SITE_tk_gimp(:) = tkgimp_lcz(utyp)
         ENDIF

         IF (.not. u_site_albr   ) THEN
            allocate( SITE_alb_roof (2, 2) )
            SITE_alb_roof(:,:) = albroof_lcz(utyp)
         ENDIF

         IF (.not. u_site_albw   ) THEN
            allocate( SITE_alb_wall (2, 2) )
            SITE_alb_wall(:,:) = albwall_lcz(utyp)
         ENDIF

         IF (.not. u_site_albgimp) THEN
            allocate( SITE_alb_gimp (2, 2) )
            SITE_alb_gimp(:,:) = albgimp_lcz(utyp)
         ENDIF

         IF (.not. u_site_albgper) THEN
            allocate( SITE_alb_gper (2, 2) )
            SITE_alb_gper(:,:) = albgper_lcz(utyp)
         ENDIF

         IF (.not. u_site_hlr  ) SITE_hlr   = hwrbld_lcz(utyp)

         IF (.not. u_site_fgper) SITE_fgimp = 1-fgper_lcz(utyp)/(1-SITE_froof)
ENDIF

IF (DEF_USE_CANYON_HWR) THEN
         SITE_hlr =SITE_hlr*(1-sqrt(SITE_froof))/sqrt(SITE_froof)
ENDIF

         ! (6) lake depth
         readflag         = u_site_lakedepth
         u_site_lakedepth = readflag .and. ncio_var_exist(fsrfdata,'lakedepth',readflag)
         IF (u_site_lakedepth) THEN
            CALL ncio_read_serial (fsrfdata, 'lakedepth', SITE_lakedepth)
         ELSE
            CALL gridlake%define_by_name ('colm_500m')
            filename = trim(DEF_dir_rawdata)//'/lake_depth.nc'
            CALL read_point_var_2d_real8 (gridlake, filename, 'lake_depth', &
               SITE_lon_location, SITE_lat_location, lakedepth)
            SITE_lakedepth = lakedepth * 0.1
         ENDIF

         write(*,'(A,F8.2,3A)') 'Lake depth : ', SITE_lakedepth, ' (from ',datasource(u_site_lakedepth),')'

         ! (7) soil brightness parameters
         readflag           = USE_SITE_soilreflectance
         u_site_soil_bright = readflag &
            .and. ncio_var_exist(fsrfdata,'soil_s_v_alb',readflag) &
            .and. ncio_var_exist(fsrfdata,'soil_d_v_alb',readflag) &
            .and. ncio_var_exist(fsrfdata,'soil_s_n_alb',readflag) &
            .and. ncio_var_exist(fsrfdata,'soil_d_n_alb',readflag)

         IF (u_site_soil_bright) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_s_v_alb', SITE_soil_s_v_alb)
            CALL ncio_read_serial (fsrfdata, 'soil_d_v_alb', SITE_soil_d_v_alb)
            CALL ncio_read_serial (fsrfdata, 'soil_s_n_alb', SITE_soil_s_n_alb)
            CALL ncio_read_serial (fsrfdata, 'soil_d_n_alb', SITE_soil_d_n_alb)
         ELSE
            SITE_soil_s_v_alb = spval
            SITE_soil_d_v_alb = spval
            SITE_soil_s_n_alb = spval
            SITE_soil_d_n_alb = spval

            CALL gridbright%define_by_name ('colm_500m')
            filename = trim(DEF_dir_rawdata)//'/soil_brightness.nc'
            CALL read_point_var_2d_int32 (gridbright, filename, 'soil_brightness', &
               SITE_lon_location, SITE_lat_location, isc)

#ifdef LULC_USGS
            IF(SITE_landtype /= 16 .and. SITE_landtype /= 24)THEN  ! NOT WATER BODIES(16)/GLACIER and ICESHEET(24)
#else
            IF(SITE_landtype /= 17 .and. SITE_landtype /= 15)THEN  ! NOT WATER BODIES(17)/GLACIER and ICE SHEET(15)
#endif
               IF ((isc >= 1) .and. (isc <= 20)) THEN
                  SITE_soil_s_v_alb = soil_s_v_refl( isc )
                  SITE_soil_d_v_alb = soil_d_v_refl( isc )
                  SITE_soil_s_n_alb = soil_s_n_refl( isc )
                  SITE_soil_d_n_alb = soil_d_n_refl( isc )
               ENDIF
            ENDIF
         ENDIF

         write(*,'(A,F8.2,3A)') 'Soil brightness s_v : ', SITE_soil_s_v_alb, ' (from ',trim(datasource(u_site_soil_bright)),')'
         write(*,'(A,F8.2,3A)') 'Soil brightness d_v : ', SITE_soil_d_v_alb, ' (from ',trim(datasource(u_site_soil_bright)),')'
         write(*,'(A,F8.2,3A)') 'Soil brightness s_n : ', SITE_soil_s_n_alb, ' (from ',trim(datasource(u_site_soil_bright)),')'
         write(*,'(A,F8.2,3A)') 'Soil brightness d_n : ', SITE_soil_d_n_alb, ' (from ',trim(datasource(u_site_soil_bright)),')'

         ! (8) soil parameters
         CALL gridsoil%define_by_name ('colm_500m')

         readflag = USE_SITE_soilparameters
         u_site_vf_quartz_mineral = readflag .and. ncio_var_exist(fsrfdata,'soil_vf_quartz_mineral',readflag)
         IF (u_site_vf_quartz_mineral) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_vf_quartz_mineral', SITE_soil_vf_quartz_mineral)
         ELSE
            allocate (SITE_soil_vf_quartz_mineral (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/vf_quartz_mineral_s.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'vf_quartz_mineral_s_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_vf_quartz_mineral(nsl))
            ENDDO
         ENDIF

         u_site_vf_gravels = readflag .and. ncio_var_exist(fsrfdata,'soil_vf_gravels',readflag)
         IF (u_site_vf_gravels) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_vf_gravels', SITE_soil_vf_gravels)
         ELSE
            allocate (SITE_soil_vf_gravels (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/vf_gravels_s.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'vf_gravels_s_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_vf_gravels(nsl))
            ENDDO
         ENDIF

         u_site_vf_sand = readflag .and. ncio_var_exist(fsrfdata,'soil_vf_sand',readflag)
         IF (u_site_vf_sand) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_vf_sand', SITE_soil_vf_sand)
         ELSE
            allocate (SITE_soil_vf_sand (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/vf_sand_s.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'vf_sand_s_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_vf_sand(nsl))
            ENDDO
         ENDIF

         u_site_vf_clay = readflag .and. ncio_var_exist(fsrfdata,'soil_vf_clay',readflag)
         IF (u_site_vf_clay) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_vf_clay', SITE_soil_vf_clay)
         ELSE
            allocate (SITE_soil_vf_clay (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/vf_clay_s.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'vf_clay_s_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_vf_clay(nsl))
            ENDDO
         ENDIF

         u_site_vf_om = readflag .and. ncio_var_exist(fsrfdata,'soil_vf_om',readflag)
         IF (u_site_vf_om) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_vf_om', SITE_soil_vf_om)
         ELSE
            allocate (SITE_soil_vf_om (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/vf_om_s.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'vf_om_s_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_vf_om(nsl))
            ENDDO
         ENDIF

         u_site_wf_gravels = readflag .and. ncio_var_exist(fsrfdata,'soil_wf_gravels',readflag)
         IF (u_site_wf_gravels) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_wf_gravels', SITE_soil_wf_gravels)
         ELSE
            allocate (SITE_soil_wf_gravels (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/wf_gravels_s.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'wf_gravels_s_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_wf_gravels(nsl))
            ENDDO
         ENDIF

         u_site_wf_sand = readflag .and. ncio_var_exist(fsrfdata,'soil_wf_sand',readflag)
         IF (u_site_wf_sand) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_wf_sand', SITE_soil_wf_sand)
         ELSE
            allocate (SITE_soil_wf_sand (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/wf_sand_s.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'wf_sand_s_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_wf_sand(nsl))
            ENDDO
         ENDIF

         u_site_wf_clay = readflag .and. ncio_var_exist(fsrfdata,'soil_wf_clay',readflag)
         IF (u_site_wf_clay) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_wf_clay', SITE_soil_wf_clay)
         ELSE
            allocate (SITE_soil_wf_clay (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/wf_clay_s.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'wf_clay_s_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_wf_clay(nsl))
            ENDDO
         ENDIF

         u_site_wf_om = readflag .and. ncio_var_exist(fsrfdata,'soil_wf_om',readflag)
         IF (u_site_wf_om) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_wf_om', SITE_soil_wf_om)
         ELSE
            allocate (SITE_soil_wf_om (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/wf_om_s.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'wf_om_s_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_wf_om(nsl))
            ENDDO
         ENDIF

         u_site_OM_density = readflag .and. ncio_var_exist(fsrfdata,'soil_OM_density',readflag)
         IF (u_site_OM_density) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_OM_density', SITE_soil_OM_density)
         ELSE
            allocate (SITE_soil_OM_density (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/OM_density_s.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'OM_density_s_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_OM_density(nsl))
            ENDDO
         ENDIF

         u_site_BD_all = readflag .and. ncio_var_exist(fsrfdata,'soil_BD_all',readflag)
         IF (u_site_BD_all) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_BD_all', SITE_soil_BD_all)
         ELSE
            allocate (SITE_soil_BD_all (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/BD_all_s.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'BD_all_s_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_BD_all(nsl))
            ENDDO
         ENDIF

         u_site_theta_s = readflag .and. ncio_var_exist(fsrfdata,'soil_theta_s',readflag)
         IF (u_site_theta_s) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_theta_s', SITE_soil_theta_s)
         ELSE
            allocate (SITE_soil_theta_s (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/theta_s.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'theta_s_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_theta_s(nsl))
            ENDDO
         ENDIF

         u_site_k_s = readflag .and. ncio_var_exist(fsrfdata,'soil_k_s',readflag)
         IF (u_site_k_s) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_k_s', SITE_soil_k_s)
         ELSE
            allocate (SITE_soil_k_s (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/k_s.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'k_s_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_k_s(nsl))
            ENDDO
         ENDIF

         u_site_csol = readflag .and. ncio_var_exist(fsrfdata,'soil_csol',readflag)
         IF (u_site_csol) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_csol', SITE_soil_csol)
         ELSE
            allocate (SITE_soil_csol (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/csol.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'csol_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_csol(nsl))
            ENDDO
         ENDIF

         u_site_tksatu = readflag .and. ncio_var_exist(fsrfdata,'soil_tksatu',readflag)
         IF (u_site_tksatu) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_tksatu', SITE_soil_tksatu)
         ELSE
            allocate (SITE_soil_tksatu (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/tksatu.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'tksatu_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_tksatu(nsl))
            ENDDO
         ENDIF

         u_site_tksatf = readflag .and. ncio_var_exist(fsrfdata,'soil_tksatf',readflag)
         IF (u_site_tksatf) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_tksatf', SITE_soil_tksatf)
         ELSE
            allocate (SITE_soil_tksatf (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/tksatf.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'tksatf_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_tksatf(nsl))
            ENDDO
         ENDIF

         u_site_tkdry = readflag .and. ncio_var_exist(fsrfdata,'soil_tkdry',readflag)
         IF (u_site_tkdry) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_tkdry', SITE_soil_tkdry)
         ELSE
            allocate (SITE_soil_tkdry (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/tkdry.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'tkdry_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_tkdry(nsl))
            ENDDO
         ENDIF

         u_site_k_solids = readflag .and. ncio_var_exist(fsrfdata,'soil_k_solids',readflag)
         IF (u_site_k_solids) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_k_solids', SITE_soil_k_solids)
         ELSE
            allocate (SITE_soil_k_solids (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/k_solids.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'k_solids_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_k_solids(nsl))
            ENDDO
         ENDIF

         u_site_psi_s = readflag .and. ncio_var_exist(fsrfdata,'soil_psi_s',readflag)
         IF (u_site_psi_s) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_psi_s', SITE_soil_psi_s)
         ELSE
            allocate (SITE_soil_psi_s (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/psi_s.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'psi_s_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_psi_s(nsl))
            ENDDO
         ENDIF

         u_site_lambda = readflag .and. ncio_var_exist(fsrfdata,'soil_lambda',readflag)
         IF (u_site_lambda) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_lambda', SITE_soil_lambda)
         ELSE
            allocate (SITE_soil_lambda (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/lambda.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'lambda_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_lambda(nsl))
            ENDDO
         ENDIF

#ifdef vanGenuchten_Mualem_SOIL_MODEL
         u_site_theta_r = readflag .and. ncio_var_exist(fsrfdata,'soil_theta_r',readflag)
         IF (u_site_theta_r) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_theta_r', SITE_soil_theta_r)
         ELSE
            allocate (SITE_soil_theta_r (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/VGM_theta_r.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'VGM_theta_r_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_theta_r(nsl))
            ENDDO
         ENDIF

         u_site_alpha_vgm = readflag .and. ncio_var_exist(fsrfdata,'soil_alpha_vgm',readflag)
         IF (u_site_alpha_vgm) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_alpha_vgm', SITE_soil_alpha_vgm)
         ELSE
            allocate (SITE_soil_alpha_vgm (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/VGM_alpha.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'VGM_alpha_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_alpha_vgm(nsl))
            ENDDO
         ENDIF

         u_site_L_vgm = readflag .and. ncio_var_exist(fsrfdata,'soil_L_vgm',readflag)
         IF (u_site_L_vgm) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_L_vgm', SITE_soil_L_vgm)
         ELSE
            allocate (SITE_soil_L_vgm (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/VGM_L.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'VGM_L_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_L_vgm(nsl))
            ENDDO
         ENDIF

         u_site_n_vgm = readflag .and. ncio_var_exist(fsrfdata,'soil_n_vgm',readflag)
         IF (u_site_n_vgm) THEN
            CALL ncio_read_serial (fsrfdata, 'soil_n_vgm', SITE_soil_n_vgm)
         ELSE
            allocate (SITE_soil_n_vgm (8))
            DO nsl = 1, 8
               write(c,'(i1)') nsl
               filename = trim(DEF_dir_rawdata)//'/soil/VGM_n.nc'
               CALL read_point_var_2d_real8 (gridsoil, filename, 'VGM_n_l'//trim(c), &
                  SITE_lon_location, SITE_lat_location, SITE_soil_n_vgm(nsl))
            ENDDO
         ENDIF
#endif

         u_site_BA_alpha = .false.
         u_site_BA_beta  = .false.
         allocate (SITE_soil_BA_alpha (8))
         allocate (SITE_soil_BA_beta  (8))
         DO nsl = 1, 8
            IF (SITE_soil_vf_gravels(nsl) + SITE_soil_vf_sand(nsl) > 0.4) THEN
               SITE_soil_BA_alpha(nsl) = 0.38
               SITE_soil_BA_beta (nsl) = 35.0
            ELSEIF (SITE_soil_vf_gravels(nsl) + SITE_soil_vf_sand(nsl) > 0.25) THEN
               SITE_soil_BA_alpha(nsl) = 0.24
               SITE_soil_BA_beta (nsl) = 26.0
            ELSE
               SITE_soil_BA_alpha(nsl) = 0.20
               SITE_soil_BA_beta (nsl) = 10.0
            ENDIF
         ENDDO

         IF (DEF_Runoff_SCHEME == 3) THEN ! for Simple VIC
            u_site_soil_texture = readflag .and. ncio_var_exist(fsrfdata,'soil_texture',readflag)
            IF (u_site_soil_texture) THEN
               CALL ncio_read_serial (fsrfdata, 'soil_texture', SITE_soil_texture)
            ELSE
               filename = trim(DEF_dir_rawdata)//'/soil/soiltexture_0cm-60cm_mean.nc'
               CALL read_point_var_2d_int32 (gridsoil, filename, 'soiltexture', &
                  SITE_lon_location, SITE_lat_location, SITE_soil_texture)
            ENDIF
         ENDIF

         write(*,'(A,8ES10.2,3A)') 'soil_vf_quartz_mineral : ', SITE_soil_vf_quartz_mineral(1:8), ' (from ',trim(datasource(u_site_vf_quartz_mineral)),')'
         write(*,'(A,8ES10.2,3A)') 'soil_vf_gravels        : ', SITE_soil_vf_gravels       (1:8), ' (from ',trim(datasource(u_site_vf_gravels       )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_vf_sand           : ', SITE_soil_vf_sand          (1:8), ' (from ',trim(datasource(u_site_vf_sand          )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_vf_clay           : ', SITE_soil_vf_clay          (1:8), ' (from ',trim(datasource(u_site_vf_clay          )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_vf_om             : ', SITE_soil_vf_om            (1:8), ' (from ',trim(datasource(u_site_vf_om            )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_wf_gravels        : ', SITE_soil_wf_gravels       (1:8), ' (from ',trim(datasource(u_site_wf_gravels       )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_wf_sand           : ', SITE_soil_wf_sand          (1:8), ' (from ',trim(datasource(u_site_wf_sand          )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_wf_clay           : ', SITE_soil_wf_clay          (1:8), ' (from ',trim(datasource(u_site_wf_clay          )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_wf_om             : ', SITE_soil_wf_om            (1:8), ' (from ',trim(datasource(u_site_wf_om            )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_OM_density        : ', SITE_soil_OM_density       (1:8), ' (from ',trim(datasource(u_site_OM_density       )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_BD_all            : ', SITE_soil_BD_all           (1:8), ' (from ',trim(datasource(u_site_BD_all           )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_theta_s           : ', SITE_soil_theta_s          (1:8), ' (from ',trim(datasource(u_site_theta_s          )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_k_s               : ', SITE_soil_k_s              (1:8), ' (from ',trim(datasource(u_site_k_s              )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_csol              : ', SITE_soil_csol             (1:8), ' (from ',trim(datasource(u_site_csol             )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_tksatu            : ', SITE_soil_tksatu           (1:8), ' (from ',trim(datasource(u_site_tksatu           )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_tksatf            : ', SITE_soil_tksatf           (1:8), ' (from ',trim(datasource(u_site_tksatf           )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_tkdry             : ', SITE_soil_tkdry            (1:8), ' (from ',trim(datasource(u_site_tkdry            )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_k_solids          : ', SITE_soil_k_solids         (1:8), ' (from ',trim(datasource(u_site_k_solids         )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_psi_s             : ', SITE_soil_psi_s            (1:8), ' (from ',trim(datasource(u_site_psi_s            )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_lambda            : ', SITE_soil_lambda           (1:8), ' (from ',trim(datasource(u_site_lambda           )),')'
#ifdef vanGenuchten_Mualem_SOIL_MODEL
         write(*,'(A,8ES10.2,3A)') 'soil_theta_r           : ', SITE_soil_theta_r          (1:8), ' (from ',trim(datasource(u_site_theta_r          )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_alpha_vgm         : ', SITE_soil_alpha_vgm        (1:8), ' (from ',trim(datasource(u_site_alpha_vgm        )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_L_vgm             : ', SITE_soil_L_vgm            (1:8), ' (from ',trim(datasource(u_site_L_vgm            )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_n_vgm             : ', SITE_soil_n_vgm            (1:8), ' (from ',trim(datasource(u_site_n_vgm            )),')'
#endif
         write(*,'(A,8ES10.2,3A)') 'soil_BA_alpha          : ', SITE_soil_BA_alpha         (1:8), ' (from ',trim(datasource(u_site_BA_alpha         )),')'
         write(*,'(A,8ES10.2,3A)') 'soil_BA_beta           : ', SITE_soil_BA_beta          (1:8), ' (from ',trim(datasource(u_site_BA_beta          )),')'

         IF (DEF_Runoff_SCHEME == 3) THEN ! for Simple VIC
            write(*,'(A,I3,3A)') 'soil texture           : ', SITE_soil_texture, ' (from ',trim(datasource(u_site_soil_texture)),')'
         ENDIF

         ! (9) depth to bedrock
         IF (DEF_USE_BEDROCK) THEN
            readflag = USE_SITE_dbedrock
            u_site_dbedrock = readflag .and. ncio_var_exist (fsrfdata, 'depth_to_bedrock',readflag)
            IF (u_site_dbedrock) THEN
               CALL ncio_read_serial (fsrfdata, 'depth_to_bedrock', SITE_dbedrock)
            ELSE
               CALL gridrock%define_by_name ('colm_500m')
               filename = trim(DEF_dir_rawdata)//'/bedrock.nc'
               CALL read_point_var_2d_real8 (gridrock, filename, 'dbedrock', &
                  SITE_lon_location, SITE_lat_location, SITE_dbedrock)
            ENDIF

            write(*,'(A,F8.2,3A)') 'Depth to bedrock : ', SITE_dbedrock, ' (from ',datasource(u_site_dbedrock),')'
         ENDIF

         ! (10) topography
         readflag = USE_SITE_topography

         u_site_elevation = readflag .and. ncio_var_exist (fsrfdata, 'elevation',readflag)
         IF (u_site_elevation) THEN
            CALL ncio_read_serial (fsrfdata, 'elevation', SITE_elevation)
         ELSE
            CALL gridtopo%define_by_name ('colm_500m')
            filename = trim(DEF_dir_rawdata)//'/elevation.nc'
            CALL read_point_var_2d_real8 (gridtopo, filename, 'elevation', &
               SITE_lon_location, SITE_lat_location, SITE_elevation)
         ENDIF

         u_site_elvstd = readflag &
            .and. ncio_var_exist (fsrfdata, 'elvstd',readflag)
         IF (u_site_elvstd) THEN
            CALL ncio_read_serial (fsrfdata, 'elvstd', SITE_elvstd)
         ELSE
            CALL gridtopo%define_by_name ('colm_500m')
            filename = trim(DEF_dir_rawdata)//'/topography.nc'
            CALL read_point_var_2d_real8 (gridtopo, filename, 'elvstd', &
               SITE_lon_location, SITE_lat_location, SITE_elvstd)
         ENDIF

         u_site_sloperatio = readflag &
            .and. ncio_var_exist (fsrfdata, 'sloperatio',readflag)
         IF (u_site_sloperatio) THEN
            CALL ncio_read_serial (fsrfdata, 'sloperatio', SITE_sloperatio)
         ELSE
            CALL gridtopo%define_by_name ('colm_500m')
            filename = trim(DEF_dir_rawdata)//'/topography.nc'
            CALL read_point_var_2d_real8 (gridtopo, filename, 'slope', &
               SITE_lon_location, SITE_lat_location, SITE_sloperatio)
         ENDIF

         IF (DEF_USE_Forcing_Downscaling) THEN

            filename = trim(DEF_DS_HiresTopographyDataDir)//"/slope.nc"
            IF (ncio_var_exist(filename,'lat') .and. ncio_var_exist(filename,'lon')) THEN
               CALL grid_topo_factor%define_from_file (filename, "lat", "lon")
            ENDIF

            u_site_svf = readflag .and. ncio_var_exist (fsrfdata, 'SITE_svf',readflag)
            IF (u_site_svf) THEN
               CALL ncio_read_serial (fsrfdata, 'SITE_svf' , SITE_svf)
            ELSE
               filename = trim(DEF_DS_HiresTopographyDataDir)//"/sky_view_factor.nc"
               CALL read_point_var_2d_real8 (grid_topo_factor, filename, 'svf', &
                  SITE_lon_location, SITE_lat_location, SITE_svf)
            ENDIF

            u_site_cur = readflag .and. ncio_var_exist (fsrfdata, 'SITE_cur',readflag)
            IF (u_site_cur) THEN
               CALL ncio_read_serial (fsrfdata, 'SITE_cur' , SITE_cur)
            ELSE
               filename = trim(DEF_DS_HiresTopographyDataDir)//"/curvature.nc"
               CALL read_point_var_2d_real8 (grid_topo_factor, filename, 'curvature', &
                  SITE_lon_location, SITE_lat_location, SITE_cur)
            ENDIF

            u_site_slp_type = readflag &
               .and. ncio_var_exist (fsrfdata, 'SITE_slp_type' ,readflag) &
               .and. ncio_var_exist (fsrfdata, 'SITE_asp_type' ,readflag) &
               .and. ncio_var_exist (fsrfdata, 'SITE_area_type',readflag)
            u_site_asp_type  = u_site_slp_type
            u_site_area_type = u_site_slp_type

            IF (u_site_slp_type) THEN
               CALL ncio_read_serial (fsrfdata, 'SITE_slp_type' , SITE_slp_type  )
               CALL ncio_read_serial (fsrfdata, 'SITE_asp_type' , SITE_asp_type  )
               CALL ncio_read_serial (fsrfdata, 'SITE_area_type', SITE_area_type )
            ELSE
               filename = trim(DEF_DS_HiresTopographyDataDir)//"/slope.nc"
               CALL read_point_var_2d_real8 (grid_topo_factor, filename, 'slope', &
                  SITE_lon_location, SITE_lat_location, slp)

               filename = trim(DEF_DS_HiresTopographyDataDir)//"/aspect.nc"
               CALL read_point_var_2d_real8 (grid_topo_factor, filename, 'aspect', &
                  SITE_lon_location, SITE_lat_location, asp)

               allocate (SITE_slp_type  (num_slope_type)); SITE_slp_type (:) = 0.
               allocate (SITE_asp_type  (num_slope_type)); SITE_asp_type (:) = 0.
               allocate (SITE_area_type (num_slope_type)); SITE_area_type(:) = 0.

               IF ((asp.ge.0 .and. asp.le.90*pi/180) .or. (asp.ge.270*pi/180 .and. asp.le.360*pi/180)) THEN
                  IF ((slp.ge.15*pi/180)) THEN  ! north abrupt slope
                     typ = 1
                  ELSE                          ! north gentle slope
                     typ = 2
                  ENDIF
               ELSE
                  IF ((slp.ge.15*pi/180)) THEN  ! south abrupt slope
                     typ = 3
                  ELSE                          ! south gentle slope
                     typ = 4
                  ENDIF
               ENDIF

               SITE_slp_type (typ) = slp
               SITE_asp_type (typ) = asp
               SITE_area_type(typ) = 1.

            ENDIF

            u_site_sf_lut = readflag .and. ncio_var_exist (fsrfdata, 'SITE_sf_lut',readflag)
            IF (u_site_sf_lut) THEN
               CALL ncio_read_serial (fsrfdata, 'SITE_sf_lut', SITE_sf_lut)
            ELSE
               filename = trim(DEF_DS_HiresTopographyDataDir)//"/terrain_elev_angle_front.nc"
               CALL read_point_var_3d_first_real8 (grid_topo_factor, filename, 'tea_front', &
                  SITE_lon_location, SITE_lat_location, num_azimuth, tea_f)

               filename = trim(DEF_DS_HiresTopographyDataDir)//"/terrain_elev_angle_back.nc"
               CALL read_point_var_3d_first_real8 (grid_topo_factor, filename, 'tea_back', &
                  SITE_lon_location, SITE_lat_location, num_azimuth, tea_b)

               allocate (SITE_sf_lut (num_azimuth, num_zenith))

               DO a = 1, num_azimuth

                  tea_f(a) = asin(max(min(tea_f(a),1.),-1.))
                  tea_b(a) = asin(max(min(tea_b(a),1.),-1.))

                  IF (tea_f(a) <= tea_b(a)) tea_f(a) = tea_b(a) + 0.001

                  DO z = 1, num_zenith
                     zenith_angle = pi/(2*num_zenith)*(z-1)

                     IF (pi*0.5 - zenith_angle < tea_b(a)) THEN
                        SITE_sf_lut(a,z) = 0
                     ELSE IF (pi*0.5 - zenith_angle > tea_f(a)) THEN
                        SITE_sf_lut(a,z) = 1
                     ELSE
                        SITE_sf_lut(a,z) = (0.5*pi - zenith_angle - tea_b(a))/(tea_f(a) - tea_b(a))
                     ENDIF

                  ENDDO
               ENDDO

            ENDIF
         ENDIF

         write(*,'(A,F8.2,3A)') 'Elevation : ', SITE_elevation, ' (from ',datasource(u_site_elevation),')'
         write(*,'(A,F8.2,3A)') 'Elv std   : ', SITE_elvstd,    ' (from ',datasource(u_site_elvstd),')'
         write(*,'(A,F8.2,3A)') 'SlopeRatio: ', SITE_sloperatio,' (from ',datasource(u_site_sloperatio),')'

         IF (DEF_USE_Forcing_Downscaling) THEN
            write(*,'(A,F8.2,3A)') 'Sky view factor : ', SITE_svf, ' (from ',datasource(u_site_svf),')'
            write(*,'(A,F8.2,3A)') 'Curvature       : ', SITE_cur, ' (from ',datasource(u_site_cur),')'
            write(c,'(I0)') num_slope_type
            write(*,'(A,'//trim(c)//'F8.2,3A)') 'Slope  type     : ', SITE_slp_type,  ' (from ',datasource(u_site_slp_type),')'
            write(*,'(A,'//trim(c)//'F8.2,3A)') 'Aspect type     : ', SITE_asp_type,  ' (from ',datasource(u_site_slp_type),')'
            write(*,'(A,'//trim(c)//'F8.2,3A)') 'Slope type area : ', SITE_area_type, ' (from ',datasource(u_site_slp_type),')'
            write(c,'(I0)') num_azimuth*num_zenith
            write(*,'(A,A,I3,A,I3,A,'//trim(c)//'F8.2,3A)') 'Shadow lookup table    : ', &
               '(', num_azimuth, ' in azimuth,', num_zenith, ' in zenith)', &
               SITE_sf_lut , ' (from ',datasource(u_site_sf_lut),')'
         ENDIF

      ELSE
         CALL ncio_read_serial (fsrfdata, 'LAI_year'      , SITE_LAI_year   )
         CALL ncio_read_serial (fsrfdata, 'TREE_LAI'      , SITE_LAI_monthly)
         CALL ncio_read_serial (fsrfdata, 'TREE_SAI'      , SITE_SAI_monthly)

         CALL ncio_read_serial (fsrfdata, 'URBAN_TYPE'    , SITE_urbtyp     )
         CALL ncio_read_serial (fsrfdata, 'LUCY_id'       , SITE_lucyid     )
         CALL ncio_read_serial (fsrfdata, 'PCT_Tree'      , SITE_fveg_urb   )
         CALL ncio_read_serial (fsrfdata, 'URBAN_TREE_TOP', SITE_htop_urb   )
         CALL ncio_read_serial (fsrfdata, 'PCT_Water'     , SITE_flake_urb  )
         CALL ncio_read_serial (fsrfdata, 'WT_ROOF'       , SITE_froof      )
         CALL ncio_read_serial (fsrfdata, 'HT_ROOF'       , SITE_hroof      )
         CALL ncio_read_serial (fsrfdata, 'WTROAD_PERV'   , SITE_fgper      )
         CALL ncio_read_serial (fsrfdata, 'BUILDING_HLR'  , SITE_hlr        )
         CALL ncio_read_serial (fsrfdata, 'POP_DEN'       , SITE_popden     )

         CALL ncio_read_serial (fsrfdata, 'EM_ROOF'       , SITE_em_roof    )
         CALL ncio_read_serial (fsrfdata, 'EM_WALL'       , SITE_em_wall    )
         CALL ncio_read_serial (fsrfdata, 'EM_IMPROAD'    , SITE_em_gimp    )
         CALL ncio_read_serial (fsrfdata, 'EM_PERROAD'    , SITE_em_gper    )
         CALL ncio_read_serial (fsrfdata, 'T_BUILDING_MAX', SITE_t_roommax  )
         CALL ncio_read_serial (fsrfdata, 'T_BUILDING_MIN', SITE_t_roommin  )
         CALL ncio_read_serial (fsrfdata, 'THICK_ROOF'    , SITE_thickroof  )
         CALL ncio_read_serial (fsrfdata, 'THICK_WALL'    , SITE_thickwall  )

         CALL ncio_read_serial (fsrfdata, 'ALB_ROOF'      , SITE_alb_roof   )
         CALL ncio_read_serial (fsrfdata, 'ALB_WALL'      , SITE_alb_wall   )
         CALL ncio_read_serial (fsrfdata, 'ALB_IMPROAD'   , SITE_alb_gimp   )
         CALL ncio_read_serial (fsrfdata, 'ALB_PERROAD'   , SITE_alb_gper   )

         CALL ncio_read_serial (fsrfdata, 'CV_ROOF'       , SITE_cv_roof    )
         CALL ncio_read_serial (fsrfdata, 'CV_WALL'       , SITE_cv_wall    )
         CALL ncio_read_serial (fsrfdata, 'CV_IMPROAD'    , SITE_cv_gimp    )
         CALL ncio_read_serial (fsrfdata, 'TK_ROOF'       , SITE_tk_roof    )
         CALL ncio_read_serial (fsrfdata, 'TK_WALL'       , SITE_tk_wall    )
         CALL ncio_read_serial (fsrfdata, 'TK_IMPROAD'    , SITE_tk_gimp    )

         CALL ncio_read_serial (fsrfdata, 'lakedepth'     , SITE_lakedepth  )

         CALL ncio_read_serial (fsrfdata, 'soil_s_v_alb'  , SITE_soil_s_v_alb)
         CALL ncio_read_serial (fsrfdata, 'soil_d_v_alb'  , SITE_soil_d_v_alb)
         CALL ncio_read_serial (fsrfdata, 'soil_s_n_alb'  , SITE_soil_s_n_alb)
         CALL ncio_read_serial (fsrfdata, 'soil_d_n_alb'  , SITE_soil_d_n_alb)

         CALL ncio_read_serial (fsrfdata, 'soil_vf_quartz_mineral', SITE_soil_vf_quartz_mineral)
         CALL ncio_read_serial (fsrfdata, 'soil_vf_gravels'       , SITE_soil_vf_gravels       )
         CALL ncio_read_serial (fsrfdata, 'soil_vf_sand'          , SITE_soil_vf_sand          )
         CALL ncio_read_serial (fsrfdata, 'soil_vf_clay'          , SITE_soil_vf_clay          )
         CALL ncio_read_serial (fsrfdata, 'soil_vf_om'            , SITE_soil_vf_om            )
         CALL ncio_read_serial (fsrfdata, 'soil_wf_gravels'       , SITE_soil_wf_gravels       )
         CALL ncio_read_serial (fsrfdata, 'soil_wf_sand'          , SITE_soil_wf_sand          )
         CALL ncio_read_serial (fsrfdata, 'soil_wf_clay'          , SITE_soil_wf_clay          )
         CALL ncio_read_serial (fsrfdata, 'soil_wf_om'            , SITE_soil_wf_om            )

         CALL ncio_read_serial (fsrfdata, 'soil_OM_density'       , SITE_soil_OM_density)
         CALL ncio_read_serial (fsrfdata, 'soil_BD_all'           , SITE_soil_BD_all    )
         CALL ncio_read_serial (fsrfdata, 'soil_theta_s'          , SITE_soil_theta_s   )
         CALL ncio_read_serial (fsrfdata, 'soil_k_s'              , SITE_soil_k_s       )
         CALL ncio_read_serial (fsrfdata, 'soil_csol'             , SITE_soil_csol      )
         CALL ncio_read_serial (fsrfdata, 'soil_tksatu'           , SITE_soil_tksatu    )
         CALL ncio_read_serial (fsrfdata, 'soil_tksatf'           , SITE_soil_tksatf    )
         CALL ncio_read_serial (fsrfdata, 'soil_tkdry'            , SITE_soil_tkdry     )
         CALL ncio_read_serial (fsrfdata, 'soil_k_solids'         , SITE_soil_k_solids  )
         CALL ncio_read_serial (fsrfdata, 'soil_lambda'           , SITE_soil_lambda    )
         CALL ncio_read_serial (fsrfdata, 'soil_psi_s'            , SITE_soil_psi_s     )

#ifdef vanGenuchten_Mualem_SOIL_MODEL
         CALL ncio_read_serial (fsrfdata, 'soil_theta_r'  , SITE_soil_theta_r  )
         CALL ncio_read_serial (fsrfdata, 'soil_alpha_vgm', SITE_soil_alpha_vgm)
         CALL ncio_read_serial (fsrfdata, 'soil_L_vgm'    , SITE_soil_L_vgm    )
         CALL ncio_read_serial (fsrfdata, 'soil_n_vgm'    , SITE_soil_n_vgm    )
#endif

         CALL ncio_read_serial (fsrfdata, 'soil_BA_alpha', SITE_soil_BA_alpha)
         CALL ncio_read_serial (fsrfdata, 'soil_BA_beta' , SITE_soil_BA_beta )

         IF (DEF_Runoff_SCHEME == 3) THEN ! for Simple VIC
            CALL ncio_read_serial (fsrfdata, 'soil_texture', SITE_soil_texture)
         ENDIF

         IF(DEF_USE_BEDROCK)THEN
            CALL ncio_read_serial (fsrfdata, 'depth_to_bedrock', SITE_dbedrock)
         ENDIF

         CALL ncio_read_serial (fsrfdata, 'elevation' , SITE_elevation )
         CALL ncio_read_serial (fsrfdata, 'elvstd'    , SITE_elvstd    )
         CALL ncio_read_serial (fsrfdata, 'sloperatio', SITE_sloperatio)

         ! used for downscaling
         IF (DEF_USE_Forcing_Downscaling) THEN
            CALL ncio_read_serial (fsrfdata, 'SITE_svf'      , SITE_svf      )
            CALL ncio_read_serial (fsrfdata, 'SITE_cur'      , SITE_cur      )
            CALL ncio_read_serial (fsrfdata, 'SITE_sf_lut'   , SITE_sf_lut   )
            CALL ncio_read_serial (fsrfdata, 'SITE_slp_type' , SITE_slp_type )
            CALL ncio_read_serial (fsrfdata, 'SITE_asp_type' , SITE_asp_type )
            CALL ncio_read_serial (fsrfdata, 'SITE_area_type', SITE_area_type)
         ENDIF
      ENDIF

      IF (.not. mksrfdata) THEN

         landpatch%nset = numpatch

         allocate (landpatch%settyp (numpatch)); landpatch%settyp = SITE_landtype

         landpatch%nblkgrp = 1
         allocate (landpatch%xblkgrp(1));       landpatch%xblkgrp(1) = 1
         allocate (landpatch%yblkgrp(1));       landpatch%yblkgrp(1) = 1

         allocate (landpatch%vecgs%vlen(1,1));  landpatch%vecgs%vlen(1,1) = numpatch
         allocate (landpatch%vecgs%vstt(1,1));  landpatch%vecgs%vstt(1,1) = 1
         allocate (landpatch%vecgs%vend(1,1));  landpatch%vecgs%vend(1,1) = numpatch

         landurban%nset = numurban

         allocate (landurban%settyp (numurban)); landurban%settyp = SITE_urbtyp

         landurban%nblkgrp = 1
         allocate (landurban%xblkgrp(1));       landurban%xblkgrp(1) = 1
         allocate (landurban%yblkgrp(1));       landurban%yblkgrp(1) = 1

         allocate (landurban%vecgs%vlen(1,1));  landurban%vecgs%vlen(1,1) = numurban
         allocate (landurban%vecgs%vstt(1,1));  landurban%vecgs%vstt(1,1) = 1
         allocate (landurban%vecgs%vend(1,1));  landurban%vecgs%vend(1,1) = numurban

         allocate (patch2urban (numurban));     patch2urban(1) = 1
         allocate (urban2patch (numurban));     urban2patch(1) = 1

         numelm = 1
         allocate (landelm%settyp   (1));  landelm%settyp  (1) = 0
         allocate (elm_patch%substt (1));  elm_patch%substt(1) = 1
         allocate (elm_patch%subend (1));  elm_patch%subend(1) = numpatch
         allocate (elm_patch%subfrc (numpatch)); elm_patch%subfrc = 1./numpatch

      ENDIF

   END SUBROUTINE read_urban_surface_data_single


!-----------------------------------------------------------------------
   SUBROUTINE write_surface_data_single (numpatch, numpft)

   USE MOD_NetCDFSerial
   USE MOD_Namelist
   USE MOD_Const_LC
   IMPLICIT NONE

   integer, intent(in) :: numpatch
   integer, intent(in), optional :: numpft

   ! Local Variables
   character(len=256) :: fsrfdata
   integer :: ipft, iyear, itime
   character(len=18)   :: source

      fsrfdata = trim(DEF_dir_landdata) // '/srfdata.nc'

      CALL ncio_create_file (fsrfdata)

      CALL ncio_define_dimension (fsrfdata, 'patch', numpatch)
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      IF (numpft > 0) THEN
         CALL ncio_define_dimension (fsrfdata, 'pft', numpft)
      ENDIF
#endif

      CALL ncio_define_dimension (fsrfdata, 'LAI_year', size(SITE_LAI_year))
      IF (DEF_LAI_MONTHLY) THEN
         CALL ncio_define_dimension (fsrfdata, 'month', 12)
      ELSE
         CALL ncio_define_dimension (fsrfdata, 'J8day', 46)
      ENDIF

      CALL ncio_define_dimension (fsrfdata, 'soil', 8)

      IF (DEF_USE_Forcing_Downscaling) THEN
         CALL ncio_define_dimension (fsrfdata, 'slope_type', num_slope_type)
         CALL ncio_define_dimension (fsrfdata, 'azi', num_azimuth)
         CALL ncio_define_dimension (fsrfdata, 'zen', num_zenith)
      ENDIF


      CALL ncio_write_serial (fsrfdata, 'latitude',  SITE_lat_location)
      CALL ncio_put_attr     (fsrfdata, 'latitude', 'units', 'degrees_north')

      CALL ncio_write_serial (fsrfdata, 'longitude', SITE_lon_location)
      CALL ncio_put_attr     (fsrfdata, 'longitude','units', 'degrees_east')

#ifdef LULC_USGS
      CALL ncio_write_serial (fsrfdata, 'USGS_classification', SITE_landtype)
      CALL ncio_put_attr     (fsrfdata, 'USGS_classification', 'source', trim(datasource(u_site_landtype)))
      CALL ncio_put_attr     (fsrfdata, 'USGS_classification', 'long_name', 'GLCC USGS Land Use/Land Cover')
#else
      CALL ncio_write_serial (fsrfdata, 'IGBP_classification', SITE_landtype)
      CALL ncio_put_attr     (fsrfdata, 'IGBP_classification', 'source', trim(datasource(u_site_landtype)))
      CALL ncio_put_attr     (fsrfdata, 'IGBP_classification', 'long_name', 'MODIS IGBP Land Use/Land Cover')
#endif

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      IF (numpft > 0) THEN
         CALL ncio_write_serial (fsrfdata, 'pfttyp',  SITE_pfttyp,  'pft')
         CALL ncio_put_attr     (fsrfdata, 'pfttyp',  'source', trim(datasource(u_site_pfts)))
         CALL ncio_put_attr     (fsrfdata, 'pfttyp',  'long_name', 'plant functional type')

         CALL ncio_write_serial (fsrfdata, 'pctpfts', SITE_pctpfts, 'pft')
         CALL ncio_put_attr     (fsrfdata, 'pctpfts', 'source', trim(datasource(u_site_pfts)))
         CALL ncio_put_attr     (fsrfdata, 'pctpfts', 'long_name', 'fraction of plant functional type')
      ENDIF
#endif
#if (defined CROP)
      IF (SITE_landtype == CROPLAND) THEN
         CALL ncio_write_serial (fsrfdata, 'croptyp', SITE_croptyp, 'patch')
         CALL ncio_put_attr     (fsrfdata, 'croptyp', 'source', trim(datasource(u_site_crop)))
         CALL ncio_put_attr     (fsrfdata, 'croptyp', 'long_name', 'crop type')

         CALL ncio_write_serial (fsrfdata, 'pctcrop', SITE_pctcrop, 'patch')
         CALL ncio_put_attr     (fsrfdata, 'pctcrop', 'source', trim(datasource(u_site_crop)))
         CALL ncio_put_attr     (fsrfdata, 'pctcrop', 'long_name', 'fraction of crop type')
      ENDIF
#endif

      CALL ncio_write_serial (fsrfdata, 'canopy_height', SITE_htop)
      CALL ncio_put_attr     (fsrfdata, 'canopy_height', 'source', trim(datasource(u_site_htop)))
      CALL ncio_put_attr     (fsrfdata, 'canopy_height', 'long_name', 'canopy height')
      CALL ncio_put_attr     (fsrfdata, 'canopy_height', 'units', 'm')
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      IF (numpft > 0) THEN
         CALL ncio_write_serial (fsrfdata, 'canopy_height_pfts', SITE_htop_pfts, 'pft')
         CALL ncio_put_attr     (fsrfdata, 'canopy_height_pfts', 'source', trim(datasource(u_site_htop)))
         CALL ncio_put_attr     (fsrfdata, 'canopy_height_pfts', 'long_name', 'canopy height')
         CALL ncio_put_attr     (fsrfdata, 'canopy_height_pfts', 'units', 'm')
      ENDIF
#endif

      source = trim(datasource(u_site_lai))
      CALL ncio_write_serial (fsrfdata, 'LAI_year', SITE_LAI_year, 'LAI_year')
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      IF (numpft > 0) THEN
         CALL ncio_write_serial (fsrfdata, 'LAI_pfts_monthly', SITE_LAI_pfts_monthly, 'pft', 'month', 'LAI_year')
         CALL ncio_put_attr     (fsrfdata, 'LAI_pfts_monthly', 'source', source)
         CALL ncio_put_attr     (fsrfdata, 'LAI_pfts_monthly', 'long_name', 'monthly leaf area index associated with PFT')

         CALL ncio_write_serial (fsrfdata, 'SAI_pfts_monthly', SITE_SAI_pfts_monthly, 'pft', 'month', 'LAI_year')
         CALL ncio_put_attr     (fsrfdata, 'SAI_pfts_monthly', 'source', source)
         CALL ncio_put_attr     (fsrfdata, 'SAI_pfts_monthly', 'long_name', 'monthly stem area index associated with PFT')
      ELSE
         CALL ncio_write_serial (fsrfdata, 'LAI_monthly', SITE_LAI_monthly, 'month', 'LAI_year')
         CALL ncio_put_attr     (fsrfdata, 'LAI_monthly', 'source', source)
         CALL ncio_put_attr     (fsrfdata, 'LAI_monthly', 'long_name', 'monthly leaf area index')

         CALL ncio_write_serial (fsrfdata, 'SAI_monthly', SITE_SAI_monthly, 'month', 'LAI_year')
         CALL ncio_put_attr     (fsrfdata, 'SAI_monthly', 'source', source)
         CALL ncio_put_attr     (fsrfdata, 'SAI_monthly', 'long_name', 'monthly stem area index')
      ENDIF
#else
      IF (DEF_LAI_MONTHLY) THEN
         CALL ncio_write_serial (fsrfdata, 'LAI_monthly', SITE_LAI_monthly, 'month', 'LAI_year')
         CALL ncio_put_attr     (fsrfdata, 'LAI_monthly', 'source', source)
         CALL ncio_put_attr     (fsrfdata, 'LAI_monthly', 'long_name', 'monthly leaf area index')

         CALL ncio_write_serial (fsrfdata, 'SAI_monthly', SITE_SAI_monthly, 'month', 'LAI_year')
         CALL ncio_put_attr     (fsrfdata, 'SAI_monthly', 'source', source)
         CALL ncio_put_attr     (fsrfdata, 'SAI_monthly', 'long_name', 'monthly stem area index')
      ELSE
         CALL ncio_write_serial (fsrfdata, 'LAI_8day', SITE_LAI_8day, 'J8day', 'LAI_year')
         CALL ncio_put_attr     (fsrfdata, 'LAI_8day', 'source', source)
         CALL ncio_put_attr     (fsrfdata, 'LAI_8day', 'long_name', '8-day leaf area index')
      ENDIF
#endif

      CALL ncio_write_serial (fsrfdata, 'lakedepth', SITE_lakedepth)
      CALL ncio_put_attr     (fsrfdata, 'lakedepth', 'source', trim(datasource(u_site_lakedepth)))
      CALL ncio_put_attr     (fsrfdata, 'lakedepth', 'long_name', 'lake depth')
      CALL ncio_put_attr     (fsrfdata, 'lakedepth', 'units', 'm')

      CALL ncio_write_serial (fsrfdata, 'soil_s_v_alb', SITE_soil_s_v_alb)
      CALL ncio_put_attr     (fsrfdata, 'soil_s_v_alb', 'source', trim(datasource(u_site_soil_bright)))
      CALL ncio_put_attr     (fsrfdata, 'soil_s_v_alb', 'long_name', 'albedo of visible of the saturated soil')

      CALL ncio_write_serial (fsrfdata, 'soil_d_v_alb', SITE_soil_d_v_alb)
      CALL ncio_put_attr     (fsrfdata, 'soil_d_v_alb', 'source', trim(datasource(u_site_soil_bright)))
      CALL ncio_put_attr     (fsrfdata, 'soil_d_v_alb', 'long_name', 'albedo of visible of the dry soil')

      CALL ncio_write_serial (fsrfdata, 'soil_s_n_alb', SITE_soil_s_n_alb)
      CALL ncio_put_attr     (fsrfdata, 'soil_s_n_alb', 'source', trim(datasource(u_site_soil_bright)))
      CALL ncio_put_attr     (fsrfdata, 'soil_s_n_alb', 'long_name', 'albedo of near infrared of the saturated soil')

      CALL ncio_write_serial (fsrfdata, 'soil_d_n_alb', SITE_soil_d_n_alb)
      CALL ncio_put_attr     (fsrfdata, 'soil_d_n_alb', 'source', trim(datasource(u_site_soil_bright)))
      CALL ncio_put_attr     (fsrfdata, 'soil_d_n_alb', 'long_name', 'albedo of near infrared of the dry soil')

      CALL ncio_write_serial (fsrfdata, 'soil_vf_quartz_mineral', SITE_soil_vf_quartz_mineral(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_quartz_mineral', 'source', trim(datasource(u_site_vf_quartz_mineral)))
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_quartz_mineral', 'long_name', 'volumetric fraction of quartz within mineral soil')

      CALL ncio_write_serial (fsrfdata, 'soil_vf_gravels', SITE_soil_vf_gravels(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_gravels', 'source', trim(datasource(u_site_vf_gravels)))
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_gravels', 'long_name', 'volumetric fraction of gravels')

      CALL ncio_write_serial (fsrfdata, 'soil_vf_sand', SITE_soil_vf_sand(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_sand', 'source', trim(datasource(u_site_vf_sand)))
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_sand', 'long_name', 'volumetric fraction of sand')

      CALL ncio_write_serial (fsrfdata, 'soil_vf_clay', SITE_soil_vf_clay(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_clay', 'source', trim(datasource(u_site_vf_clay)))
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_clay', 'long_name', 'volumetric fraction of clay')

      CALL ncio_write_serial (fsrfdata, 'soil_vf_om', SITE_soil_vf_om(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_om', 'source', trim(datasource(u_site_vf_om)))
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_om', 'long_name', 'volumetric fraction of organic matter')

      CALL ncio_write_serial (fsrfdata, 'soil_wf_gravels', SITE_soil_wf_gravels(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_wf_gravels', 'source', trim(datasource(u_site_wf_gravels)))
      CALL ncio_put_attr     (fsrfdata, 'soil_wf_gravels', 'long_name', 'gravimetric fraction of gravels')

      CALL ncio_write_serial (fsrfdata, 'soil_wf_sand', SITE_soil_wf_sand(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_wf_sand', 'source', trim(datasource(u_site_wf_sand)))
      CALL ncio_put_attr     (fsrfdata, 'soil_wf_sand', 'long_name', 'gravimetric fraction of sand')

      CALL ncio_write_serial (fsrfdata, 'soil_wf_clay', SITE_soil_wf_clay(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_wf_clay', 'source', trim(datasource(u_site_wf_clay)))
      CALL ncio_put_attr     (fsrfdata, 'soil_wf_clay', 'long_name', 'gravimetric fraction of clay')

      CALL ncio_write_serial (fsrfdata, 'soil_wf_om', SITE_soil_wf_om(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_wf_om', 'source', trim(datasource(u_site_wf_om)))
      CALL ncio_put_attr     (fsrfdata, 'soil_wf_om', 'long_name', 'gravimetric fraction of om')

      CALL ncio_write_serial (fsrfdata, 'soil_OM_density', SITE_soil_OM_density(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_OM_density', 'source', trim(datasource(u_site_OM_density)))
      CALL ncio_put_attr     (fsrfdata, 'soil_OM_density', 'long_name', 'OM density')
      CALL ncio_put_attr     (fsrfdata, 'soil_OM_density', 'units', 'kg/m3')

      CALL ncio_write_serial (fsrfdata, 'soil_BD_all', SITE_soil_BD_all(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_BD_all', 'source', trim(datasource(u_site_BD_all)))
      CALL ncio_put_attr     (fsrfdata, 'soil_BD_all', 'long_name', 'bulk density of soil (GRAVELS + OM + mineral soils)')
      CALL ncio_put_attr     (fsrfdata, 'soil_BD_all', 'units', 'kg/m3')

      CALL ncio_write_serial (fsrfdata, 'soil_theta_s', SITE_soil_theta_s(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_theta_s', 'source', trim(datasource(u_site_theta_s)))
      CALL ncio_put_attr     (fsrfdata, 'soil_theta_s', 'long_name', 'saturated water content')
      CALL ncio_put_attr     (fsrfdata, 'soil_theta_s', 'units', 'cm3/cm3')

      CALL ncio_write_serial (fsrfdata, 'soil_k_s', SITE_soil_k_s(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_k_s', 'source', trim(datasource(u_site_k_s)))
      CALL ncio_put_attr     (fsrfdata, 'soil_k_s', 'long_name', 'saturated hydraulic conductivity')
      CALL ncio_put_attr     (fsrfdata, 'soil_k_s', 'units', 'cm/day')

      CALL ncio_write_serial (fsrfdata, 'soil_csol', SITE_soil_csol(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_csol', 'source', trim(datasource(u_site_csol)))
      CALL ncio_put_attr     (fsrfdata, 'soil_csol', 'long_name', 'heat capacity of soil solids')
      CALL ncio_put_attr     (fsrfdata, 'soil_csol', 'units', 'J/(m3 K)')

      CALL ncio_write_serial (fsrfdata, 'soil_tksatu', SITE_soil_tksatu(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_tksatu', 'source', trim(datasource(u_site_tksatu)))
      CALL ncio_put_attr     (fsrfdata, 'soil_tksatu', 'long_name', 'thermal conductivity of saturated unfrozen soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_tksatu', 'units', 'W/m-K')

      CALL ncio_write_serial (fsrfdata, 'soil_tksatf', SITE_soil_tksatf(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_tksatf', 'source', trim(datasource(u_site_tksatf)))
      CALL ncio_put_attr     (fsrfdata, 'soil_tksatf', 'long_name', 'thermal conductivity of saturated frozen soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_tksatf', 'units', 'W/m-K')

      CALL ncio_write_serial (fsrfdata, 'soil_tkdry', SITE_soil_tkdry(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_tkdry', 'source', trim(datasource(u_site_tkdry)))
      CALL ncio_put_attr     (fsrfdata, 'soil_tkdry', 'long_name', 'thermal conductivity for dry soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_tkdry', 'units', 'W/(m-K)')

      CALL ncio_write_serial (fsrfdata, 'soil_k_solids', SITE_soil_k_solids(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_k_solids', 'source', trim(datasource(u_site_k_solids)))
      CALL ncio_put_attr     (fsrfdata, 'soil_k_solids', 'long_name', 'thermal conductivity of minerals soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_k_solids', 'units', 'W/m-K')

      CALL ncio_write_serial (fsrfdata, 'soil_lambda', SITE_soil_lambda(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_lambda', 'source', trim(datasource(u_site_lambda)))
      CALL ncio_put_attr     (fsrfdata, 'soil_lambda', 'long_name', 'pore size distribution index (dimensionless)')

      CALL ncio_write_serial (fsrfdata, 'soil_psi_s', SITE_soil_psi_s(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_psi_s ', 'source', trim(datasource(u_site_psi_s)))
      CALL ncio_put_attr     (fsrfdata, 'soil_psi_s ', 'long_name', 'matric potential at saturation')
      CALL ncio_put_attr     (fsrfdata, 'soil_psi_s ', 'units', 'cm')

#ifdef vanGenuchten_Mualem_SOIL_MODEL
      CALL ncio_write_serial (fsrfdata, 'soil_theta_r', SITE_soil_theta_r(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_theta_r', 'source', trim(datasource(u_site_theta_r)))
      CALL ncio_put_attr     (fsrfdata, 'soil_theta_r', 'long_name', 'residual water content')
      CALL ncio_put_attr     (fsrfdata, 'soil_theta_r', 'units', 'cm3/cm3')

      CALL ncio_write_serial (fsrfdata, 'soil_alpha_vgm', SITE_soil_alpha_vgm(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_alpha_vgm', 'source', trim(datasource(u_site_alpha_vgm)))
      CALL ncio_put_attr     (fsrfdata, 'soil_alpha_vgm', 'long_name', 'a parameter corresponding approximately to the inverse of the air-entry value')

      CALL ncio_write_serial (fsrfdata, 'soil_L_vgm', SITE_soil_L_vgm(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_L_vgm', 'source', trim(datasource(u_site_L_vgm)))
      CALL ncio_put_attr     (fsrfdata, 'soil_L_vgm', 'long_name', 'pore-connectivity parameter [dimensionless]')

      CALL ncio_write_serial (fsrfdata, 'soil_n_vgm', SITE_soil_n_vgm(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_n_vgm', 'source', trim(datasource(u_site_n_vgm)))
      CALL ncio_put_attr     (fsrfdata, 'soil_n_vgm', 'long_name', 'a shape parameter [dimensionless]')

#endif

      CALL ncio_write_serial (fsrfdata, 'soil_BA_alpha', SITE_soil_BA_alpha(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_BA_alpha', 'source', trim(datasource(u_site_BA_alpha)))
      CALL ncio_put_attr     (fsrfdata, 'soil_BA_alpha', 'long_name', 'alpha in Balland and Arp(2005) thermal conductivity scheme')

      CALL ncio_write_serial (fsrfdata, 'soil_BA_beta', SITE_soil_BA_beta(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_BA_beta', 'source', trim(datasource(u_site_BA_beta)))
      CALL ncio_put_attr     (fsrfdata, 'soil_BA_beta', 'long_name', 'beta in Balland and Arp(2005) thermal conductivity scheme')

      IF (DEF_Runoff_SCHEME == 3) THEN ! for Simple VIC
         CALL ncio_write_serial (fsrfdata, 'soil_texture ', SITE_soil_texture)
         CALL ncio_put_attr     (fsrfdata, 'soil_texture ', 'source', trim(datasource(u_site_soil_texture)))
         CALL ncio_put_attr     (fsrfdata, 'soil_texture ', 'long_name', 'USDA soil texture')
      ENDIF

      IF(DEF_USE_BEDROCK)THEN
         CALL ncio_write_serial (fsrfdata, 'depth_to_bedrock', SITE_dbedrock)
         CALL ncio_put_attr     (fsrfdata, 'depth_to_bedrock', 'source', trim(datasource(u_site_dbedrock)))
      ENDIF

      CALL ncio_write_serial (fsrfdata, 'elevation', SITE_elevation)
      CALL ncio_put_attr     (fsrfdata, 'elevation', 'source', trim(datasource(u_site_elevation)))

      CALL ncio_write_serial (fsrfdata, 'elvstd', SITE_elvstd)
      CALL ncio_put_attr     (fsrfdata, 'elvstd', 'source', trim(datasource(u_site_elvstd)))
      CALL ncio_put_attr     (fsrfdata, 'elvstd', 'long_name', 'standard deviation of elevation')

      CALL ncio_write_serial (fsrfdata, 'sloperatio', SITE_sloperatio)
      CALL ncio_put_attr     (fsrfdata, 'sloperatio', 'source', trim(datasource(u_site_sloperatio)))
      CALL ncio_put_attr     (fsrfdata, 'sloperatio', 'long_name', 'slope ratio')

      ! used for downscaling
      IF (DEF_USE_Forcing_Downscaling) THEN
         CALL ncio_write_serial (fsrfdata, 'SITE_svf', SITE_svf)
         CALL ncio_put_attr     (fsrfdata, 'SITE_svf','source', trim(datasource(u_site_svf)))
         CALL ncio_put_attr     (fsrfdata, 'SITE_svf','long_name', 'sky view factor')

         CALL ncio_write_serial (fsrfdata, 'SITE_cur', SITE_cur)
         CALL ncio_put_attr     (fsrfdata, 'SITE_cur','source', trim(datasource(u_site_cur)))
         CALL ncio_put_attr     (fsrfdata, 'SITE_cur','long_name', 'curvature')

         CALL ncio_write_serial (fsrfdata, 'SITE_sf_lut', SITE_sf_lut, 'azi', 'zen')
         CALL ncio_put_attr     (fsrfdata, 'SITE_sf_lut','source', trim(datasource(u_site_sf_lut)))
         CALL ncio_put_attr     (fsrfdata, 'SITE_sf_lut','long_name', 'look up table of shadow factor')

         CALL ncio_write_serial (fsrfdata, 'SITE_slp_type' , SITE_slp_type , 'type')
         CALL ncio_put_attr     (fsrfdata, 'SITE_slp_type','source', trim(datasource(u_site_slp_type)))
         CALL ncio_put_attr     (fsrfdata, 'SITE_slp_type','long_name', 'topographic slope of each character')

         CALL ncio_write_serial (fsrfdata, 'SITE_asp_type' , SITE_asp_type , 'type')
         CALL ncio_put_attr     (fsrfdata, 'SITE_asp_type','source', trim(datasource(u_site_asp_type)))
         CALL ncio_put_attr     (fsrfdata, 'SITE_asp_type','long_name', 'topographic aspect of each character')

         CALL ncio_write_serial (fsrfdata, 'SITE_area_type', SITE_area_type, 'type')
         CALL ncio_put_attr     (fsrfdata, 'SITE_area_type','source', trim(datasource(u_site_area_type)))
         CALL ncio_put_attr     (fsrfdata, 'SITE_area_type','long_name', 'area percentage of each character')
      ENDIF

   END SUBROUTINE write_surface_data_single


!-----------------------------------------------------------------------
   SUBROUTINE write_urban_surface_data_single (numurban)

   USE MOD_NetCDFSerial
   USE MOD_Namelist
   USE MOD_Const_LC
   IMPLICIT NONE

   integer, intent(in) :: numurban

   ! Local Variables
   character(len=256) :: fsrfdata
   integer :: iyear, itime
   character(len=8)   :: source

      fsrfdata = trim(DEF_dir_landdata) // '/srfdata.nc'

      SITE_fgper     = 1 - (SITE_fgimp-SITE_froof)/(1-SITE_froof-SITE_flake_urb)
      SITE_froof     = SITE_froof/(1-SITE_flake_urb)
      SITE_fveg_urb  = SITE_fveg_urb  * 100
      SITE_flake_urb = SITE_flake_urb * 100

      CALL ncio_create_file (fsrfdata)

      CALL ncio_define_dimension (fsrfdata, 'soil', 8)

      CALL ncio_define_dimension (fsrfdata, 'azi', num_azimuth)
      CALL ncio_define_dimension (fsrfdata, 'zen', num_zenith)
      CALL ncio_define_dimension (fsrfdata, 'slope_type', num_slope_type)

      CALL ncio_define_dimension (fsrfdata, 'patch', numurban)

      CALL ncio_define_dimension (fsrfdata, 'LAI_year', size(SITE_LAI_year))
      CALL ncio_define_dimension (fsrfdata, 'month', 12)

      CALL ncio_define_dimension (fsrfdata, 'ulev'    , 10)
      CALL ncio_define_dimension (fsrfdata, 'numsolar', 2 )
      CALL ncio_define_dimension (fsrfdata, 'numrad'  , 2 )

      CALL ncio_write_serial (fsrfdata, 'latitude' , SITE_lat_location)
      CALL ncio_write_serial (fsrfdata, 'longitude', SITE_lon_location)

      CALL ncio_write_serial (fsrfdata, 'LAI_year', SITE_LAI_year, 'LAI_year')
      CALL ncio_write_serial (fsrfdata, 'TREE_LAI', SITE_LAI_monthly, 'month', 'LAI_year')
      CALL ncio_write_serial (fsrfdata, 'TREE_SAI', SITE_SAI_monthly, 'month', 'LAI_year')
      CALL ncio_put_attr     (fsrfdata, 'TREE_LAI', 'source', trim(datasource(u_site_urblai)))
      CALL ncio_put_attr     (fsrfdata, 'TREE_SAI', 'source', trim(datasource(u_site_urbsai)))

      CALL ncio_write_serial (fsrfdata, 'URBAN_TYPE'    , SITE_urbtyp     )
      CALL ncio_write_serial (fsrfdata, 'LUCY_id'       , SITE_lucyid     )
      CALL ncio_write_serial (fsrfdata, 'PCT_Tree'      , SITE_fveg_urb   )
      CALL ncio_write_serial (fsrfdata, 'URBAN_TREE_TOP', SITE_htop_urb   )
      CALL ncio_write_serial (fsrfdata, 'PCT_Water'     , SITE_flake_urb  )
      CALL ncio_write_serial (fsrfdata, 'WT_ROOF'       , SITE_froof      )
      CALL ncio_write_serial (fsrfdata, 'HT_ROOF'       , SITE_hroof      )
      CALL ncio_write_serial (fsrfdata, 'WTROAD_PERV'   , SITE_fgper      )
      CALL ncio_write_serial (fsrfdata, 'BUILDING_HLR'  , SITE_hlr        )
      CALL ncio_write_serial (fsrfdata, 'POP_DEN'       , SITE_popden     )

      CALL ncio_put_attr     (fsrfdata, 'PCT_Tree'      , 'source', trim(datasource(u_site_fveg )))
      CALL ncio_put_attr     (fsrfdata, 'URBAN_TREE_TOP', 'source', trim(datasource(u_site_htopu)))
      CALL ncio_put_attr     (fsrfdata, 'PCT_Water'     , 'source', trim(datasource(u_site_flake)))
      CALL ncio_put_attr     (fsrfdata, 'WT_ROOF'       , 'source', trim(datasource(u_site_froof)))
      CALL ncio_put_attr     (fsrfdata, 'HT_ROOF'       , 'source', trim(datasource(u_site_hroof)))
      CALL ncio_put_attr     (fsrfdata, 'WTROAD_PERV'   , 'source', trim(datasource(u_site_fgper)))
      CALL ncio_put_attr     (fsrfdata, 'BUILDING_HLR'  , 'source', trim(datasource(u_site_hlr  )))
      CALL ncio_put_attr     (fsrfdata, 'POP_DEN'       , 'source', trim(datasource(u_site_pop  )))

      CALL ncio_write_serial (fsrfdata, 'EM_ROOF'       , SITE_em_roof   )
      CALL ncio_write_serial (fsrfdata, 'EM_WALL'       , SITE_em_wall   )
      CALL ncio_write_serial (fsrfdata, 'EM_IMPROAD'    , SITE_em_gimp   )
      CALL ncio_write_serial (fsrfdata, 'EM_PERROAD'    , SITE_em_gper   )
      CALL ncio_write_serial (fsrfdata, 'T_BUILDING_MAX', SITE_t_roommax )
      CALL ncio_write_serial (fsrfdata, 'T_BUILDING_MIN', SITE_t_roommin )
      CALL ncio_write_serial (fsrfdata, 'THICK_ROOF'    , SITE_thickroof )
      CALL ncio_write_serial (fsrfdata, 'THICK_WALL'    , SITE_thickwall )

      CALL ncio_put_attr     (fsrfdata, 'EM_ROOF'       , 'source', trim(datasource(u_site_emr   )))
      CALL ncio_put_attr     (fsrfdata, 'EM_WALL'       , 'source', trim(datasource(u_site_emw   )))
      CALL ncio_put_attr     (fsrfdata, 'EM_IMPROAD'    , 'source', trim(datasource(u_site_emgimp)))
      CALL ncio_put_attr     (fsrfdata, 'EM_PERROAD'    , 'source', trim(datasource(u_site_emgper)))
      CALL ncio_put_attr     (fsrfdata, 'T_BUILDING_MAX', 'source', trim(datasource(u_site_tbmax )))
      CALL ncio_put_attr     (fsrfdata, 'T_BUILDING_MIN', 'source', trim(datasource(u_site_tbmin )))
      CALL ncio_put_attr     (fsrfdata, 'THICK_ROOF'    , 'source', trim(datasource(u_site_thickr)))
      CALL ncio_put_attr     (fsrfdata, 'THICK_WALL'    , 'source', trim(datasource(u_site_thickw)))

      CALL ncio_write_serial (fsrfdata, 'ALB_ROOF'      , SITE_alb_roof   , 'numrad', 'numsolar')
      CALL ncio_write_serial (fsrfdata, 'ALB_WALL'      , SITE_alb_wall   , 'numrad', 'numsolar')
      CALL ncio_write_serial (fsrfdata, 'ALB_IMPROAD'   , SITE_alb_gimp   , 'numrad', 'numsolar')
      CALL ncio_write_serial (fsrfdata, 'ALB_PERROAD'   , SITE_alb_gper   , 'numrad', 'numsolar')

      CALL ncio_put_attr     (fsrfdata, 'ALB_ROOF'      , 'source', trim(datasource(u_site_albr   )))
      CALL ncio_put_attr     (fsrfdata, 'ALB_WALL'      , 'source', trim(datasource(u_site_albw   )))
      CALL ncio_put_attr     (fsrfdata, 'ALB_IMPROAD'   , 'source', trim(datasource(u_site_albgimp)))
      CALL ncio_put_attr     (fsrfdata, 'ALB_PERROAD'   , 'source', trim(datasource(u_site_albgper)))

      CALL ncio_write_serial (fsrfdata, 'CV_ROOF'       , SITE_cv_roof    , 'ulev')
      CALL ncio_write_serial (fsrfdata, 'CV_WALL'       , SITE_cv_wall    , 'ulev')
      CALL ncio_write_serial (fsrfdata, 'CV_IMPROAD'    , SITE_cv_gimp    , 'ulev')
      CALL ncio_write_serial (fsrfdata, 'TK_ROOF'       , SITE_tk_roof    , 'ulev')
      CALL ncio_write_serial (fsrfdata, 'TK_WALL'       , SITE_tk_wall    , 'ulev')
      CALL ncio_write_serial (fsrfdata, 'TK_IMPROAD'    , SITE_tk_gimp    , 'ulev')

      CALL ncio_put_attr     (fsrfdata, 'CV_ROOF'       , 'source', trim(datasource(u_site_cvr   )))
      CALL ncio_put_attr     (fsrfdata, 'CV_WALL'       , 'source', trim(datasource(u_site_cvw   )))
      CALL ncio_put_attr     (fsrfdata, 'CV_IMPROAD'    , 'source', trim(datasource(u_site_cvgimp)))
      CALL ncio_put_attr     (fsrfdata, 'TK_ROOF'       , 'source', trim(datasource(u_site_tkr   )))
      CALL ncio_put_attr     (fsrfdata, 'TK_WALL'       , 'source', trim(datasource(u_site_tkw   )))
      CALL ncio_put_attr     (fsrfdata, 'TK_IMPROAD'    , 'source', trim(datasource(u_site_tkgimp)))

      CALL ncio_write_serial (fsrfdata, 'lakedepth', SITE_lakedepth)
      CALL ncio_put_attr     (fsrfdata, 'lakedepth', 'source', trim(datasource(u_site_lakedepth)))
      CALL ncio_put_attr     (fsrfdata, 'lakedepth', 'long_name', 'lake depth')
      CALL ncio_put_attr     (fsrfdata, 'lakedepth', 'units', 'm')

      CALL ncio_write_serial (fsrfdata, 'soil_s_v_alb', SITE_soil_s_v_alb)
      CALL ncio_put_attr     (fsrfdata, 'soil_s_v_alb', 'source', trim(datasource(u_site_soil_bright)))
      CALL ncio_put_attr     (fsrfdata, 'soil_s_v_alb', 'long_name', 'albedo of visible of the saturated soil')

      CALL ncio_write_serial (fsrfdata, 'soil_d_v_alb', SITE_soil_d_v_alb)
      CALL ncio_put_attr     (fsrfdata, 'soil_d_v_alb', 'source', trim(datasource(u_site_soil_bright)))
      CALL ncio_put_attr     (fsrfdata, 'soil_d_v_alb', 'long_name', 'albedo of visible of the dry soil')

      CALL ncio_write_serial (fsrfdata, 'soil_s_n_alb', SITE_soil_s_n_alb)
      CALL ncio_put_attr     (fsrfdata, 'soil_s_n_alb', 'source', trim(datasource(u_site_soil_bright)))
      CALL ncio_put_attr     (fsrfdata, 'soil_s_n_alb', 'long_name', 'albedo of near infrared of the saturated soil')

      CALL ncio_write_serial (fsrfdata, 'soil_d_n_alb', SITE_soil_d_n_alb)
      CALL ncio_put_attr     (fsrfdata, 'soil_d_n_alb', 'source', trim(datasource(u_site_soil_bright)))
      CALL ncio_put_attr     (fsrfdata, 'soil_d_n_alb', 'long_name', 'albedo of near infrared of the dry soil')

      CALL ncio_write_serial (fsrfdata, 'soil_vf_quartz_mineral', SITE_soil_vf_quartz_mineral(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_quartz_mineral', 'source', trim(datasource(u_site_vf_quartz_mineral)))
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_quartz_mineral', 'long_name', 'volumetric fraction of quartz within mineral soil')

      CALL ncio_write_serial (fsrfdata, 'soil_vf_gravels', SITE_soil_vf_gravels(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_gravels', 'source', trim(datasource(u_site_vf_gravels)))
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_gravels', 'long_name', 'volumetric fraction of gravels')

      CALL ncio_write_serial (fsrfdata, 'soil_vf_sand', SITE_soil_vf_sand(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_sand', 'source', trim(datasource(u_site_vf_sand)))
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_sand', 'long_name', 'volumetric fraction of sand')

      CALL ncio_write_serial (fsrfdata, 'soil_vf_clay', SITE_soil_vf_clay(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_clay', 'source', trim(datasource(u_site_vf_clay)))
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_clay', 'long_name', 'volumetric fraction of clay')

      CALL ncio_write_serial (fsrfdata, 'soil_vf_om', SITE_soil_vf_om(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_om', 'source', trim(datasource(u_site_vf_om)))
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_om', 'long_name', 'volumetric fraction of organic matter')

      CALL ncio_write_serial (fsrfdata, 'soil_wf_gravels', SITE_soil_wf_gravels(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_wf_gravels', 'source', trim(datasource(u_site_wf_gravels)))
      CALL ncio_put_attr     (fsrfdata, 'soil_wf_gravels', 'long_name', 'gravimetric fraction of gravels')

      CALL ncio_write_serial (fsrfdata, 'soil_wf_sand', SITE_soil_wf_sand(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_wf_sand', 'source', trim(datasource(u_site_wf_sand)))
      CALL ncio_put_attr     (fsrfdata, 'soil_wf_sand', 'long_name', 'gravimetric fraction of sand')

      CALL ncio_write_serial (fsrfdata, 'soil_wf_clay', SITE_soil_wf_clay(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_wf_clay', 'source', trim(datasource(u_site_wf_clay)))
      CALL ncio_put_attr     (fsrfdata, 'soil_wf_clay', 'long_name', 'gravimetric fraction of clay')

      CALL ncio_write_serial (fsrfdata, 'soil_wf_om', SITE_soil_wf_om(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_wf_om', 'source', trim(datasource(u_site_wf_om)))
      CALL ncio_put_attr     (fsrfdata, 'soil_wf_om', 'long_name', 'gravimetric fraction of om')

      CALL ncio_write_serial (fsrfdata, 'soil_OM_density', SITE_soil_OM_density(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_OM_density', 'source', trim(datasource(u_site_OM_density)))
      CALL ncio_put_attr     (fsrfdata, 'soil_OM_density', 'long_name', 'OM density')
      CALL ncio_put_attr     (fsrfdata, 'soil_OM_density', 'units', 'kg/m3')

      CALL ncio_write_serial (fsrfdata, 'soil_BD_all', SITE_soil_BD_all(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_BD_all', 'source', trim(datasource(u_site_BD_all)))
      CALL ncio_put_attr     (fsrfdata, 'soil_BD_all', 'long_name', 'bulk density of soil (GRAVELS + OM + mineral soils)')
      CALL ncio_put_attr     (fsrfdata, 'soil_BD_all', 'units', 'kg/m3')

      CALL ncio_write_serial (fsrfdata, 'soil_theta_s', SITE_soil_theta_s(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_theta_s', 'source', trim(datasource(u_site_theta_s)))
      CALL ncio_put_attr     (fsrfdata, 'soil_theta_s', 'long_name', 'saturated water content')
      CALL ncio_put_attr     (fsrfdata, 'soil_theta_s', 'units', 'cm3/cm3')

      CALL ncio_write_serial (fsrfdata, 'soil_k_s', SITE_soil_k_s(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_k_s', 'source', trim(datasource(u_site_k_s)))
      CALL ncio_put_attr     (fsrfdata, 'soil_k_s', 'long_name', 'saturated hydraulic conductivity')
      CALL ncio_put_attr     (fsrfdata, 'soil_k_s', 'units', 'cm/day')

      CALL ncio_write_serial (fsrfdata, 'soil_csol', SITE_soil_csol(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_csol', 'source', trim(datasource(u_site_csol)))
      CALL ncio_put_attr     (fsrfdata, 'soil_csol', 'long_name', 'heat capacity of soil solids')
      CALL ncio_put_attr     (fsrfdata, 'soil_csol', 'units', 'J/(m3 K)')

      CALL ncio_write_serial (fsrfdata, 'soil_tksatu', SITE_soil_tksatu(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_tksatu', 'source', trim(datasource(u_site_tksatu)))
      CALL ncio_put_attr     (fsrfdata, 'soil_tksatu', 'long_name', 'thermal conductivity of saturated unfrozen soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_tksatu', 'units', 'W/m-K')

      CALL ncio_write_serial (fsrfdata, 'soil_tksatf', SITE_soil_tksatf(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_tksatf', 'source', trim(datasource(u_site_tksatf)))
      CALL ncio_put_attr     (fsrfdata, 'soil_tksatf', 'long_name', 'thermal conductivity of saturated frozen soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_tksatf', 'units', 'W/m-K')

      CALL ncio_write_serial (fsrfdata, 'soil_tkdry', SITE_soil_tkdry(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_tkdry', 'source', trim(datasource(u_site_tkdry)))
      CALL ncio_put_attr     (fsrfdata, 'soil_tkdry', 'long_name', 'thermal conductivity for dry soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_tkdry', 'units', 'W/(m-K)')

      CALL ncio_write_serial (fsrfdata, 'soil_k_solids', SITE_soil_k_solids(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_k_solids', 'source', trim(datasource(u_site_k_solids)))
      CALL ncio_put_attr     (fsrfdata, 'soil_k_solids', 'long_name', 'thermal conductivity of minerals soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_k_solids', 'units', 'W/m-K')

      CALL ncio_write_serial (fsrfdata, 'soil_lambda', SITE_soil_lambda(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_lambda', 'source', trim(datasource(u_site_lambda)))
      CALL ncio_put_attr     (fsrfdata, 'soil_lambda', 'long_name', 'pore size distribution index (dimensionless)')

      CALL ncio_write_serial (fsrfdata, 'soil_psi_s', SITE_soil_psi_s(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_psi_s ', 'source', trim(datasource(u_site_psi_s)))
      CALL ncio_put_attr     (fsrfdata, 'soil_psi_s ', 'long_name', 'matric potential at saturation')
      CALL ncio_put_attr     (fsrfdata, 'soil_psi_s ', 'units', 'cm')

#ifdef vanGenuchten_Mualem_SOIL_MODEL
      CALL ncio_write_serial (fsrfdata, 'soil_theta_r', SITE_soil_theta_r(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_theta_r', 'source', trim(datasource(u_site_theta_r)))
      CALL ncio_put_attr     (fsrfdata, 'soil_theta_r', 'long_name', 'residual water content')
      CALL ncio_put_attr     (fsrfdata, 'soil_theta_r', 'units', 'cm3/cm3')

      CALL ncio_write_serial (fsrfdata, 'soil_alpha_vgm', SITE_soil_alpha_vgm(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_alpha_vgm', 'source', trim(datasource(u_site_alpha_vgm)))
      CALL ncio_put_attr     (fsrfdata, 'soil_alpha_vgm', 'long_name', 'a parameter corresponding approximately to the inverse of the air-entry value')

      CALL ncio_write_serial (fsrfdata, 'soil_L_vgm', SITE_soil_L_vgm(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_L_vgm', 'source', trim(datasource(u_site_L_vgm)))
      CALL ncio_put_attr     (fsrfdata, 'soil_L_vgm', 'long_name', 'pore-connectivity parameter [dimensionless]')

      CALL ncio_write_serial (fsrfdata, 'soil_n_vgm', SITE_soil_n_vgm(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_n_vgm', 'source', trim(datasource(u_site_n_vgm)))
      CALL ncio_put_attr     (fsrfdata, 'soil_n_vgm', 'long_name', 'a shape parameter [dimensionless]')

#endif

      CALL ncio_write_serial (fsrfdata, 'soil_BA_alpha', SITE_soil_BA_alpha(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_BA_alpha', 'source', trim(datasource(u_site_BA_alpha)))
      CALL ncio_put_attr     (fsrfdata, 'soil_BA_alpha', 'long_name', 'alpha in Balland and Arp(2005) thermal conductivity scheme')

      CALL ncio_write_serial (fsrfdata, 'soil_BA_beta', SITE_soil_BA_beta(1:8), 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_BA_beta', 'source', trim(datasource(u_site_BA_beta)))
      CALL ncio_put_attr     (fsrfdata, 'soil_BA_beta', 'long_name', 'beta in Balland and Arp(2005) thermal conductivity scheme')

      IF (DEF_Runoff_SCHEME == 3) THEN ! for Simple VIC
         CALL ncio_write_serial (fsrfdata, 'soil_texture ', SITE_soil_texture)
         CALL ncio_put_attr     (fsrfdata, 'soil_texture ', 'source', trim(datasource(u_site_soil_texture)))
         CALL ncio_put_attr     (fsrfdata, 'soil_texture ', 'long_name', 'USDA soil texture')
      ENDIF

      IF(DEF_USE_BEDROCK)THEN
         CALL ncio_write_serial (fsrfdata, 'depth_to_bedrock', SITE_dbedrock)
         CALL ncio_put_attr     (fsrfdata, 'depth_to_bedrock', 'source', trim(datasource(u_site_dbedrock)))
      ENDIF

      CALL ncio_write_serial (fsrfdata, 'elevation', SITE_elevation)
      CALL ncio_put_attr     (fsrfdata, 'elevation', 'source', trim(datasource(u_site_elevation)))

      CALL ncio_write_serial (fsrfdata, 'elvstd', SITE_elvstd)
      CALL ncio_put_attr     (fsrfdata, 'elvstd', 'source', trim(datasource(u_site_elvstd)))
      CALL ncio_put_attr     (fsrfdata, 'elvstd', 'long_name', 'standard deviation of elevation')

      CALL ncio_write_serial (fsrfdata, 'sloperatio', SITE_sloperatio)
      CALL ncio_put_attr     (fsrfdata, 'sloperatio', 'source', trim(datasource(u_site_sloperatio)))
      CALL ncio_put_attr     (fsrfdata, 'sloperatio', 'long_name', 'slope ratio')

      ! used for downscaling
      IF (DEF_USE_Forcing_Downscaling) THEN
         CALL ncio_write_serial (fsrfdata, 'SITE_svf', SITE_svf)
         CALL ncio_put_attr     (fsrfdata, 'SITE_svf','source', trim(datasource(u_site_svf)))
         CALL ncio_put_attr     (fsrfdata, 'SITE_svf','long_name', 'sky view factor')

         CALL ncio_write_serial (fsrfdata, 'SITE_cur', SITE_cur)
         CALL ncio_put_attr     (fsrfdata, 'SITE_cur','source', trim(datasource(u_site_cur)))
         CALL ncio_put_attr     (fsrfdata, 'SITE_cur','long_name', 'curvature')

         CALL ncio_write_serial (fsrfdata, 'SITE_sf_lut', SITE_sf_lut, 'azi', 'zen')
         CALL ncio_put_attr     (fsrfdata, 'SITE_sf_lut','source', trim(datasource(u_site_sf_lut)))
         CALL ncio_put_attr     (fsrfdata, 'SITE_sf_lut','long_name', 'look up table of shadow factor')

         CALL ncio_write_serial (fsrfdata, 'SITE_slp_type' , SITE_slp_type , 'type')
         CALL ncio_put_attr     (fsrfdata, 'SITE_slp_type','source', trim(datasource(u_site_slp_type)))
         CALL ncio_put_attr     (fsrfdata, 'SITE_slp_type','long_name', 'topographic slope of each character')

         CALL ncio_write_serial (fsrfdata, 'SITE_asp_type' , SITE_asp_type , 'type')
         CALL ncio_put_attr     (fsrfdata, 'SITE_asp_type','source', trim(datasource(u_site_asp_type)))
         CALL ncio_put_attr     (fsrfdata, 'SITE_asp_type','long_name', 'topographic aspect of each character')

         CALL ncio_write_serial (fsrfdata, 'SITE_area_type', SITE_area_type, 'type')
         CALL ncio_put_attr     (fsrfdata, 'SITE_area_type','source', trim(datasource(u_site_area_type)))
         CALL ncio_put_attr     (fsrfdata, 'SITE_area_type','long_name', 'area percentage of each character')
      ENDIF

   END SUBROUTINE write_urban_surface_data_single


!-----------------------------------------------------------------------
   character(len=18) FUNCTION datasource (is_site)

   IMPLICIT NONE
   logical, intent(in) :: is_site

      IF (is_site) THEN
         datasource = 'SITE'
      ELSE
         datasource = 'CoLM 2024 raw data'
      ENDIF

   END FUNCTION datasource


!-----------------------------------------------------------------------
   SUBROUTINE single_srfdata_final ()

   IMPLICIT NONE

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      IF (allocated(SITE_pfttyp )) deallocate(SITE_pfttyp )
      IF (allocated(SITE_pctpfts)) deallocate(SITE_pctpfts)
#endif

#ifdef CROP
      IF (allocated(SITE_croptyp)) deallocate(SITE_croptyp)
      IF (allocated(SITE_pctcrop)) deallocate(SITE_pctcrop)
#endif

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      IF (allocated(SITE_htop_pfts)) deallocate(SITE_htop_pfts)
#endif

      IF (allocated(SITE_LAI_monthly)) deallocate(SITE_LAI_monthly)
      IF (allocated(SITE_SAI_monthly)) deallocate(SITE_SAI_monthly)
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      IF (allocated(SITE_LAI_pfts_monthly)) deallocate(SITE_LAI_pfts_monthly)
      IF (allocated(SITE_SAI_pfts_monthly)) deallocate(SITE_SAI_pfts_monthly)
#endif

      IF (allocated(SITE_LAI_year)) deallocate(SITE_LAI_year)
      IF (allocated(SITE_LAI_8day)) deallocate(SITE_LAI_8day)

      IF (allocated(SITE_soil_vf_quartz_mineral)) deallocate(SITE_soil_vf_quartz_mineral)
      IF (allocated(SITE_soil_vf_gravels       )) deallocate(SITE_soil_vf_gravels       )
      IF (allocated(SITE_soil_vf_sand          )) deallocate(SITE_soil_vf_sand          )
      IF (allocated(SITE_soil_vf_clay          )) deallocate(SITE_soil_vf_clay          )
      IF (allocated(SITE_soil_vf_om            )) deallocate(SITE_soil_vf_om            )
      IF (allocated(SITE_soil_wf_gravels       )) deallocate(SITE_soil_wf_gravels       )
      IF (allocated(SITE_soil_wf_sand          )) deallocate(SITE_soil_wf_sand          )
      IF (allocated(SITE_soil_wf_clay          )) deallocate(SITE_soil_wf_clay          )
      IF (allocated(SITE_soil_wf_om            )) deallocate(SITE_soil_wf_om            )
      IF (allocated(SITE_soil_OM_density       )) deallocate(SITE_soil_OM_density       )
      IF (allocated(SITE_soil_BD_all           )) deallocate(SITE_soil_BD_all           )
      IF (allocated(SITE_soil_theta_s          )) deallocate(SITE_soil_theta_s          )
      IF (allocated(SITE_soil_k_s              )) deallocate(SITE_soil_k_s              )
      IF (allocated(SITE_soil_csol             )) deallocate(SITE_soil_csol             )
      IF (allocated(SITE_soil_tksatu           )) deallocate(SITE_soil_tksatu           )
      IF (allocated(SITE_soil_tksatf           )) deallocate(SITE_soil_tksatf           )
      IF (allocated(SITE_soil_tkdry            )) deallocate(SITE_soil_tkdry            )
      IF (allocated(SITE_soil_k_solids         )) deallocate(SITE_soil_k_solids         )
      IF (allocated(SITE_soil_psi_s            )) deallocate(SITE_soil_psi_s            )
      IF (allocated(SITE_soil_lambda           )) deallocate(SITE_soil_lambda           )
      IF (allocated(SITE_soil_theta_r          )) deallocate(SITE_soil_theta_r          )
#ifdef vanGenuchten_Mualem_SOIL_MODEL
      IF (allocated(SITE_soil_alpha_vgm        )) deallocate(SITE_soil_alpha_vgm        )
      IF (allocated(SITE_soil_L_vgm            )) deallocate(SITE_soil_L_vgm            )
      IF (allocated(SITE_soil_n_vgm            )) deallocate(SITE_soil_n_vgm            )
#endif
      IF (allocated(SITE_soil_BA_alpha         )) deallocate(SITE_soil_BA_alpha         )
      IF (allocated(SITE_soil_BA_beta          )) deallocate(SITE_soil_BA_beta          )

      IF (allocated(SITE_sf_lut                )) deallocate(SITE_sf_lut                )
      IF (allocated(SITE_slp_type              )) deallocate(SITE_slp_type              )
      IF (allocated(SITE_asp_type              )) deallocate(SITE_asp_type              )
      IF (allocated(SITE_area_type             )) deallocate(SITE_area_type             )

   END SUBROUTINE single_srfdata_final

END MODULE MOD_SingleSrfdata
#endif
