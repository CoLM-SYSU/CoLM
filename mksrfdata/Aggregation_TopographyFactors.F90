#include <define.h>

SUBROUTINE Aggregation_TopographyFactors ( &
      grid_topo_factor , dir_topodata, dir_model_landdata, lc_year)
!-----------------------------------------------------------------------
!  Global topography-based factors data
!
!  Created by Sisi Chen, Lu Li, 06/2024
!-----------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_LandPatch
   USE MOD_Land2mWMO
   USE MOD_NetCDFVector
   USE MOD_NetCDFBlock
#ifdef RangeCheck
   USE MOD_RangeCheck
#endif
   USE MOD_AggregationRequestData
   USE MOD_Utils
#ifdef SrfdataDiag
   USE MOD_Mesh, only: numelm
   USE MOD_LandElm
   USE MOD_SrfdataDiag
#endif

   IMPLICIT NONE

   ! arguments:
   ! ---------------------------------------------------------------
   integer, intent(in) :: lc_year
   type(grid_type),  intent(in) :: grid_topo_factor    ! Grid structure for high resolution topography factors
   character(len=*), intent(in) :: dir_topodata        ! Direct of Rawdata
   character(len=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ---------------------------------------------------------------
   character(len=256) :: landdir, lndname, cyear
   character(len=3)   :: sdir, sdir1

   type (block_data_real8_2d) :: slp_grid    ! slope
   type (block_data_real8_2d) :: asp_grid    ! aspect
   type (block_data_real8_2d) :: svf_grid    ! sky view factor
   type (block_data_real8_2d) :: cur_grid    ! curvature
   type (block_data_real8_3d) :: tea_f_grid  ! sine of terrain elevation angle at front of grid
   type (block_data_real8_3d) :: tea_b_grid  ! sine of terrain elevation angle at back of grid

   ! patch
   real(r8), allocatable :: svf_patches (:)
   real(r8), allocatable :: cur_patches (:)
   real(r8), allocatable :: tea_f_azi_patches (:,:) ! shape as (azimuth, patches)
   real(r8), allocatable :: tea_b_azi_patches (:,:)
   real(r8), allocatable :: sf_lut_patches  (:,:,:) ! shape as (azimuth, zenith, patches)
   real(r8), allocatable :: sf_curve_patches(:,:,:) ! shape as (azimuth, parameters, patches)

   ! four defined types at all patches
   real(r8), allocatable :: asp_type_patches  (:,:) ! shape as (type, patches)
   real(r8), allocatable :: slp_type_patches  (:,:)
   real(r8), allocatable :: area_type_patches (:,:)

   ! pixelsets
   real(r8), allocatable :: slp_one         (:)
   real(r8), allocatable :: asp_one         (:)
   real(r8), allocatable :: svf_one         (:)
   real(r8), allocatable :: cur_one         (:)
   real(r8), allocatable :: area_one        (:)
   real(r8), allocatable :: sf_one          (:)
   real(r8), allocatable :: tea_f_azi_one (:,:)
   real(r8), allocatable :: tea_b_azi_one (:,:)
   real(r8), allocatable :: tea_f_one       (:)
   real(r8), allocatable :: tea_b_one       (:)
   logical , allocatable :: sf_mask_one     (:)
   logical , allocatable :: asp_mask_one    (:)
   logical , allocatable :: area_mask_one   (:)
   logical , allocatable :: slp_mask_one    (:)

   ! pixelsets of four defined types at each patch
   real(r8), allocatable :: asp_type_one  (:,:)
   real(r8), allocatable :: slp_type_one  (:,:)
   real(r8), allocatable :: area_type_one (:,:)

   real(r8) :: sum_area_one                ! sum of pixel area of a patch
   real(r8) :: zenith_angle(num_zenith)    ! sine of sun zenith angle (divided by num_zenith part)

   ! Intermediate variables used in reducing the dimensionality of masking factors
   real(r8) :: x2_sum
   real(r8) :: x_sum
   real(r8) :: y_sum
   real(r8) :: xy_sum
   real(r8) :: a1                          ! Function Parameters
   real(r8) :: a2                          ! Function Parameters
   real(r8) :: y(num_zenith)               ! shadow factor under a single azimuth and single patch
   real(r8) :: x(num_zenith)               ! Solar zenith angle used as a predictor
   real(r8), allocatable :: y_train(:)     ! The part of y used to fit the function
   real(r8), allocatable :: x_train(:)     ! The part of x used to fit the function
   real(r8), allocatable :: y_train_transform(:)       ! The transform function of y_train

   ! local variables
   integer :: ipatch, i, ps, pe, type, a, z, count_pixels, num_pixels, j, index, n
   integer :: wmo_src

#ifdef SrfdataDiag
   integer :: typpatch(N_land_classification+1), ityp  ! number of land classification
#endif
   write(cyear,'(i4.4)') lc_year
   landdir = trim(dir_model_landdata) // '/topography/' // trim(cyear)

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif
   IF (p_is_master) THEN
      write(*,'(/, A)') 'Aggregate topography factor ...'
      CALL system('mkdir -p ' // trim(adjustl(landdir)))
   ENDIF
#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

   ! -------------------------------------------------------------------
   ! read topography-based factor data
   ! -------------------------------------------------------------------
   IF (p_is_io) THEN
      lndname = trim(dir_topodata)//"/slope.nc"
      CALL allocate_block_data (grid_topo_factor, slp_grid)
      CALL ncio_read_block (lndname, 'slope', grid_topo_factor, slp_grid)

      lndname = trim(dir_topodata)//"/aspect.nc"
      CALL allocate_block_data (grid_topo_factor, asp_grid)
      CALL ncio_read_block (lndname, 'aspect', grid_topo_factor, asp_grid)

      lndname = trim(dir_topodata)//"/terrain_elev_angle_front.nc"
      CALL allocate_block_data (grid_topo_factor, tea_f_grid, num_azimuth)
      CALL ncio_read_block (lndname, 'tea_front', grid_topo_factor, num_azimuth, tea_f_grid)

      lndname = trim(dir_topodata)//"/terrain_elev_angle_back.nc"
      CALL allocate_block_data (grid_topo_factor, tea_b_grid, num_azimuth)
      CALL ncio_read_block (lndname, 'tea_back', grid_topo_factor, num_azimuth, tea_b_grid)

      lndname = trim(dir_topodata)//"/sky_view_factor.nc"
      CALL allocate_block_data (grid_topo_factor, svf_grid)
      CALL ncio_read_block (lndname, 'svf', grid_topo_factor, svf_grid)

      lndname = trim(dir_topodata)//"/curvature.nc"
      CALL allocate_block_data (grid_topo_factor, cur_grid)
      CALL ncio_read_block (lndname, 'curvature', grid_topo_factor, cur_grid)

   ! --------------------------------------------------------------------------
   ! aggregate the terrain factor data from the resolution of raw data to patch
   ! --------------------------------------------------------------------------
#ifdef USEMPI
      ! mpi send
      CALL aggregation_data_daemon (  grid_topo_factor, &
         data_r8_2d_in1 = slp_grid,   data_r8_2d_in2 = asp_grid,  &
         data_r8_2d_in3 = svf_grid,   data_r8_2d_in4 = cur_grid,  &
         data_r8_3d_in1 = tea_f_grid, n1_r8_3d_in1 = num_azimuth, &
         data_r8_3d_in2 = tea_b_grid, n1_r8_3d_in2 = num_azimuth)
#endif
   ENDIF


   IF (p_is_worker) THEN
      ! allocate for output variables at patches
      allocate (svf_patches      (numpatch))
      allocate (cur_patches      (numpatch))
      allocate (asp_type_patches (num_slope_type, numpatch))
      allocate (slp_type_patches (num_slope_type, numpatch))
      allocate (area_type_patches(num_slope_type, numpatch))
      allocate (sf_lut_patches   (num_azimuth, num_zenith, numpatch))
      ! generate sine of sun zenith angles at equal intervals
      DO i = 1, num_zenith
         zenith_angle(i) = pi/(2*num_zenith)*(i-1)
      ENDDO

      ! aggregate loop
      DO ipatch = 1, numpatch

         IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
           wmo_src = wmo_source (landpatch%ielm(ipatch))

           svf_patches (ipatch) = svf_patches (wmo_src)
           cur_patches (ipatch) = cur_patches (wmo_src)
           asp_type_patches (:, ipatch) = asp_type_patches (:, wmo_src)
           slp_type_patches (:, ipatch) = slp_type_patches (:, wmo_src)
           area_type_patches(:, ipatch) = area_type_patches(:, wmo_src)
           sf_lut_patches   (:, :, ipatch) = sf_lut_patches(:, :, wmo_src)
           CYCLE
         ENDIF

         CALL aggregation_request_data (landpatch, ipatch, grid_topo_factor, &
            zip = USE_zip_for_aggregation, area = area_one, &
            data_r8_2d_in1 = slp_grid,   data_r8_2d_out1 = slp_one, &
            data_r8_2d_in2 = asp_grid,   data_r8_2d_out2 = asp_one, &
            data_r8_2d_in3 = svf_grid,   data_r8_2d_out3 = svf_one, &
            data_r8_2d_in4 = cur_grid,   data_r8_2d_out4 = cur_one, &
            data_r8_3d_in1 = tea_f_grid, data_r8_3d_out1 = tea_f_azi_one, n1_r8_3d_in1 = num_azimuth, &
            data_r8_3d_in2 = tea_b_grid, data_r8_3d_out2 = tea_b_azi_one, n1_r8_3d_in2 = num_azimuth)

         ! ------------------------------------------------------------------
         ! aggregate sky view factor, curvature at patches
         ! ------------------------------------------------------------------
         IF (any(svf_one /= -9999.0)) THEN
            svf_patches (ipatch) = &
               sum(svf_one * area_one, mask = svf_one /= -9999.0) &
                  / sum(area_one, mask = svf_one /= -9999.0)
         ELSE
            svf_patches (ipatch) = -1.0e36
         ENDIF

         IF (any(cur_one /= -9999.0)) THEN
            cur_patches (ipatch) = &
               sum(cur_one * area_one, mask = cur_one /= -9999.0) &
                  / sum(area_one, mask = cur_one /= -9999.0)
         ELSE
            cur_patches (ipatch) = -1.0e36
         ENDIF

         ! ------------------------------------------------------------------------------
         ! aggregate look up table of shadow factor at patches
         ! ------------------------------------------------------------------------------
         ! number of pixels
         num_pixels = size(area_one)

         ! allocate pixel variables
         allocate(tea_f_one   (num_pixels))
         allocate(tea_b_one   (num_pixels))
         allocate(sf_one      (num_pixels))
         allocate(sf_mask_one (num_pixels))

         ! sum of areas of one patch
         sum_area_one = sum(area_one, mask = area_one>0)

         DO a = 1, num_azimuth
            DO z = 1, num_zenith
               ! terrain elevation angle at each azimuth
               tea_f_one(:) = tea_f_azi_one(a,:)
               tea_b_one(:) = tea_b_azi_one(a,:)

               ! count the pixels which are not missing value
               count_pixels = 0

               DO i = 1, num_pixels
                  IF ((isnan(tea_f_one(i))).or.(isnan(tea_b_one(i)))) THEN
                     sf_one(i) = 1    ! Not consider the effect of casting shadows
                  ELSE
                     IF (tea_f_one(i)>1)  tea_f_one(i) = 1
                     IF (tea_f_one(i)<-1) tea_f_one(i) = -1
                     IF (tea_b_one(i)>1)  tea_b_one(i) = 1
                     IF (tea_b_one(i)<-1) tea_b_one(i) = -1
                     tea_f_one(i) = asin(tea_f_one(i))
                     tea_b_one(i) = asin(tea_b_one(i))

                     ! Compare the sun's altitude angle to the terrain's altitude angle
                     ! to determine the value of the shadow factor.
                     ! -----------------------------------------------------------------
                     ! Sisi Chen, Lu Li, Yongjiu Dai, et al. Exploring Topography Downscaling
                     !     Methods for Hyper-Resolution Land Surface Modeling.
                     !     DOI: 10.22541/au.171403656.68476353/v1
                     ! -----------------------------------------------------------------
                     IF ((tea_b_one(i) /= -9999).and.(tea_f_one(i) /= -9999)) THEN
                        count_pixels = count_pixels+1

                        IF (pi*0.5 - zenith_angle(z) < tea_b_one(i)) THEN
                           sf_one(i) = 0
                        ELSEIF (pi*0.5 - zenith_angle(z) > tea_f_one(i)) THEN
                           sf_one(i) = 1
                        ELSE
                           IF (tea_f_one(i).eq.tea_b_one(i)) tea_f_one(i) = tea_b_one(i)+0.001
                           sf_one(i) = (0.5*pi - zenith_angle(z) - tea_b_one(i))/(tea_f_one(i) - tea_b_one(i))
                        ENDIF

                     ENDIF
                  ENDIF
               ENDDO
               sf_mask_one(:) = sf_one(:) /= -9999
               sf_lut_patches(a,z,ipatch) = sum(sf_one(:), mask = sf_mask_one)/count_pixels
            ENDDO
         ENDDO

         ! deallocate
         deallocate(tea_f_one)
         deallocate(tea_b_one)
         deallocate(sf_one)
         deallocate(sf_mask_one)

         ! -----------------------------------------------------------------------------------------------
         ! aggregate slope and aspect at four defined types at patches
         ! -----------------------------------------------------------------------------------------------
         ! allocate pixelsets variables
         allocate(asp_type_one(1:num_slope_type,1:num_pixels))
         allocate(slp_type_one(1:num_slope_type,1:num_pixels))
         allocate(area_type_one(1:num_slope_type,1:num_pixels))
         allocate(slp_mask_one(1:num_pixels))
         allocate(asp_mask_one(1:num_pixels))
         allocate(area_mask_one(1:num_pixels))

         DO i = 1, num_slope_type
            asp_type_one(i,:) = -9999
            slp_type_one(i,:) = -9999
            area_type_one(i,:) = -9999
         ENDDO

         DO i = 1, num_pixels
            ! Define the south slope, north slope, abrupt slope and gentle lope of target pixel
            IF (((asp_one(i).ge.0 .and. asp_one(i).le.90*pi/180) .or. (asp_one(i).ge.270*pi/180 .and. &
                 asp_one(i).le.360*pi/180)).and.(slp_one(i).ge.15*pi/180)) THEN    ! north abrupt slope
                 type = 1
            ELSEIF (((asp_one(i).ge.0 .and. asp_one(i).le.90*pi/180) .or. (asp_one(i).ge.270*pi/180 .and. &
                      asp_one(i).le.360*pi/180)).and.(slp_one(i)<15*pi/180)) THEN  ! north gentle slope
                 type = 2
            ELSEIF ((asp_one(i).gt.90*pi/180) .and. (asp_one(i).lt.270*pi/180) .and. &
                     (slp_one(i).ge.15*pi/180)) THEN  ! south abrupt slope
                 type = 3
            ELSEIF ((asp_one(i).gt.90*pi/180) .and. (asp_one(i).lt.270*pi/180) .and. &
                     (slp_one(i).lt.15*pi/180)) THEN  ! south gentle slope
                 type = 4
            ELSE ! missing value=-9999
                 CYCLE
            ENDIF

            IF ((area_one(i)>0).and.(area_one(i)<=sum_area_one)) THEN      ! quality control
                  area_type_one(type,i) = area_one(i)
                  asp_type_one (type,i) = asp_one(i)*area_one(i)
                  slp_type_one (type,i) = slp_one(i)*area_one(i)
            ENDIF
         ENDDO

         ! assign value to four types at patches
         DO i = 1, num_slope_type
            IF (sum_area_one.eq.0.0) THEN
               area_type_patches(i,ipatch) = 0
               asp_type_patches(i,ipatch) = 0
               slp_type_patches(i,ipatch) = 0
            ELSE
               area_mask_one(:) = area_type_one(i,:) /= -9999
               area_type_patches(i,ipatch) = sum(area_type_one(i,:), mask = area_mask_one(:))/sum_area_one
               asp_mask_one(:) = asp_type_one(i,:) /= -9999
               asp_type_patches(i,ipatch) = sum(asp_type_one(i,:), mask = asp_mask_one(:))/sum_area_one
               slp_mask_one(:) = slp_type_one(i,:) /= -9999
               slp_type_patches(i,ipatch) = sum(slp_type_one(i,:), mask = slp_mask_one(:))/sum_area_one
            ENDIF
         ENDDO

         ! deallocate
         deallocate(asp_type_one)
         deallocate(slp_type_one)
         deallocate(area_type_one)
         deallocate(slp_mask_one)
         deallocate(asp_mask_one)
         deallocate(area_mask_one)
      ENDDO


#ifdef USEMPI
      CALL aggregation_worker_done ()
#endif
   ENDIF

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
   CALL check_vector_data ('svf_patches       ', svf_patches      )
   CALL check_vector_data ('cur_patches       ', cur_patches      )
   CALL check_vector_data ('slp_type_patches  ', slp_type_patches )
   CALL check_vector_data ('asp_type_patches  ', asp_type_patches )
   CALL check_vector_data ('area_type_patches ', area_type_patches)
   CALL check_vector_data ('sf_lut_patches    ', sf_lut_patches   )
#endif

! Reduce the dimension of the shadow factor array
! Construct a new array with dimensions of sf_curve_patches(azimuth, shadow factor parameters, patches)
   IF (p_is_worker) THEN
      allocate(sf_curve_patches(num_azimuth,num_zenith_parameter,numpatch))
      DO a = 1, num_azimuth
         DO ipatch = 1, numpatch
            y(:) = sf_lut_patches(a,:,ipatch)
            x(:) = zenith_angle(:)

            ! Obtain the last position of y==1
            index = 1
            DO z = 1, num_zenith-1
               IF ((y(z)==1.).and.(y(z+1)<1.)) THEN
                  index = z+1
               ENDIF
            ENDDO

            ! allocate Allocate memory to dynamic arrays
            n = num_zenith - index +1
            allocate(y_train(n))
            allocate(x_train(n))
            allocate(y_train_transform(n))

            ! Obtain the predicted value y_train and prediction factor x_train for fitting,
            ! the form of the fitting function is
            ! ln(-ln(y_train)) = a1*x_train+a2
            y_train(:) = y(index:)
            x_train(:) = x(index:)

            ! Transform y_train to enable linear regression fitting
            DO i = 1, n
               IF (y_train(i) <= 0.) y_train(i) = 0.001
               IF (y_train(i) >= 1.) y_train(i) = 0.999
            ENDDO
            y_train_transform(:) = log(-1*log(y_train(:)))

            ! Obtain parameters a1 and a2 using the least squares method
            x_sum = 0.
            xy_sum = 0.
            y_sum = 0.
            x2_sum = 0.
            DO z = 1, n
               xy_sum = xy_sum + x_train(z)*y_train_transform(z)
               x_sum = x_sum + x_train(z)
               y_sum = y_sum + y_train_transform(z)
               x2_sum = x2_sum + x_train(z)*x_train(z)
            ENDDO
            IF (n*x2_sum - x_sum*x_sum == 0.) THEN
               a1 = 0
               a2 = 0
            ELSE
               a1 = (n*xy_sum - x_sum*y_sum)/(n*x2_sum - x_sum*x_sum)
               a2 = (y_sum - a1*x_sum)/n
            ENDIF
            sf_curve_patches(a,1,ipatch) = x(index-1) ! Minimum zenith angle at which occlusion begins
            sf_curve_patches(a,2,ipatch) = a1
            sf_curve_patches(a,3,ipatch) = a2

            ! deallocate
            deallocate(y_train)
            deallocate(x_train)
            deallocate(y_train_transform)
         ENDDO
      ENDDO
   ENDIF

   lndname = trim(landdir)//'/svf_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
   CALL ncio_write_vector (lndname, 'svf_patches', 'patch', landpatch, svf_patches, 1)

   lndname = trim(landdir)//'/cur_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
   CALL ncio_write_vector (lndname, 'cur_patches', 'patch', landpatch, cur_patches, 1)

   lndname = trim(landdir)//'/slp_type_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
   CALL ncio_define_dimension_vector (lndname, landpatch, 'slope_type', num_slope_type)
   CALL ncio_write_vector (lndname, 'slp_type_patches', 'slope_type', num_slope_type, 'patch', landpatch, slp_type_patches, 1)

   lndname = trim(landdir)//'/asp_type_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
   CALL ncio_define_dimension_vector (lndname, landpatch, 'slope_type', num_slope_type)
   CALL ncio_write_vector (lndname, 'asp_type_patches', 'slope_type', num_slope_type, 'patch', landpatch, asp_type_patches, 1)

   lndname = trim(landdir)//'/area_type_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
   CALL ncio_define_dimension_vector (lndname, landpatch, 'slope_type', num_slope_type)
   CALL ncio_write_vector (lndname, 'area_type_patches', 'slope_type', num_slope_type, 'patch', landpatch, area_type_patches, 1)

   lndname = trim(landdir)//'/sf_curve_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
   CALL ncio_define_dimension_vector (lndname, landpatch, 'azimuth', num_azimuth)
   CALL ncio_define_dimension_vector (lndname, landpatch, 'zenith_p', num_zenith_parameter)
   CALL ncio_write_vector (lndname, 'sf_curve_patches', 'azimuth', num_azimuth, 'zenith_p', num_zenith_parameter, 'patch', &
                           landpatch, sf_curve_patches, 1)

#ifdef SrfdataDiag
   typpatch = (/(ityp, ityp = 0, N_land_classification)/)

   ! only write the first type of slope and aspect at patches
   lndname  = trim(dir_model_landdata) // '/diag/topo_factor_slp_' // trim(cyear) // '.nc'
   DO i = 1, num_slope_type
      write(sdir,'(I0)') i
      CALL srfdata_map_and_write (slp_type_patches(i,:), landpatch%settyp, typpatch, m_patch2diag, &
         -1.0e36_r8, lndname, 'slp_'//trim(sdir), compress = 1, write_mode = 'one')
   ENDDO

   lndname  = trim(dir_model_landdata) // '/diag/topo_factor_asp_' // trim(cyear) // '.nc'
   DO i = 1, num_slope_type
      write(sdir,'(I0)') i
      CALL srfdata_map_and_write (asp_type_patches(i,:), landpatch%settyp, typpatch, m_patch2diag, &
         -1.0e36_r8, lndname, 'asp_'//trim(sdir), compress = 1, write_mode = 'one')
   ENDDO

   lndname  = trim(dir_model_landdata) // '/diag/topo_factor_svf_' // trim(cyear) // '.nc'
   CALL srfdata_map_and_write (svf_patches, landpatch%settyp, typpatch, m_patch2diag, &
      -1.0e36_r8, lndname, 'svf', compress = 1, write_mode = 'one')

   lndname  = trim(dir_model_landdata) // '/diag/topo_factor_cur_' // trim(cyear) // '.nc'
   CALL srfdata_map_and_write (cur_patches, landpatch%settyp, typpatch, m_patch2diag, &
      -1.0e36_r8, lndname, 'cur', compress = 1, write_mode = 'one')

   lndname  = trim(dir_model_landdata) // '/diag/topo_factor_sf_lut_' // trim(cyear) // '.nc'

   DO j = 1, num_azimuth
      DO i = 1, num_zenith
         write(sdir,'(I0)') j
         write(sdir1,'(I0)') i
         CALL srfdata_map_and_write (sf_lut_patches(j,i,:), landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, 'sf_'//trim(sdir)//'_'//trim(sdir1), compress = 1, write_mode = 'one')
      ENDDO
   ENDDO
#endif

   IF (p_is_worker) THEN
      IF (allocated(slp_type_patches)) deallocate ( slp_type_patches )
      IF (allocated(asp_type_patches)) deallocate ( asp_type_patches )
      IF (allocated(sf_lut_patches  )) deallocate ( sf_lut_patches   )
      IF (allocated(sf_curve_patches)) deallocate ( sf_curve_patches )
      IF (allocated(svf_patches     )) deallocate ( svf_patches      )
      IF (allocated(cur_patches     )) deallocate ( cur_patches      )
   ENDIF

END SUBROUTINE Aggregation_TopographyFactors
