#include <define.h>

SUBROUTINE Aggregation_SoilBrightness ( &
      gland, dir_rawdata, dir_model_landdata, lc_year)
!-----------------------------------------------------------------------
!  Creates land model surface dataset from original "raw" data files -
!      data with 30 arc seconds resolution
!
!  Created by Yongjiu Dai, 03/2014
!
! !REVISIONS:
!  Shupeng Zhang, 01/2022: porting codes to MPI parallel version.
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Block
   USE MOD_Grid
   USE MOD_LandPatch
   USE MOD_Land2mWMO
   USE MOD_NetCDFBlock
   USE MOD_NetCDFVector
#ifdef RangeCheck
   USE MOD_RangeCheck
#endif
   USE MOD_AggregationRequestData
   USE MOD_Utils
#ifdef SrfdataDiag
   USE MOD_SrfdataDiag
#endif
   USE MOD_Vars_Global, only: spval

   IMPLICIT NONE

   ! arguments:

   integer, intent(in) :: lc_year
   type(grid_type),  intent(in) :: gland
   character(len=*), intent(in) :: dir_rawdata
   character(len=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ---------------------------------------------------------------
   character(len=256) :: landdir, lndname, cyear

   type (block_data_int32_2d) :: isc
   type (block_data_real8_2d) :: a_s_v_refl
   type (block_data_real8_2d) :: a_d_v_refl
   type (block_data_real8_2d) :: a_s_n_refl
   type (block_data_real8_2d) :: a_d_n_refl

   real(r8), allocatable :: soil_s_v_alb (:)
   real(r8), allocatable :: soil_d_v_alb (:)
   real(r8), allocatable :: soil_s_n_alb (:)
   real(r8), allocatable :: soil_d_n_alb (:)

   integer :: ii, L
   integer :: ipatch, iblkme, iblk, jblk, ix, iy
   integer :: wmo_src
   real(r8), allocatable :: soil_one(:)

#ifdef SrfdataDiag
   integer :: typpatch(N_land_classification+1), ityp
#endif

   ! ----------------------------------------------------------------------
   ! The soil color and reflectance is from the work:
   ! Peter J. Lawrence and Thomas N. Chase, 2007:
   ! Representing a MODIS consistent land surface in the Community Land Model (CLM 3.0):
   ! Part 1 generating MODIS consistent land surface parameters
   real(r8) soil_s_v_refl(20) ! Saturated visible soil reflectance
   real(r8) soil_d_v_refl(20) ! Dry visible soil reflectance
   real(r8) soil_s_n_refl(20) ! Saturated near infrared soil reflectance
   real(r8) soil_d_n_refl(20) ! Dry near infrared soil reflectance

      soil_s_v_refl = (/ 0.26, 0.24, 0.22, 0.20, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14, &
         0.13, 0.12, 0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04 /)

      soil_d_v_refl = (/ 0.37, 0.35, 0.33, 0.31, 0.30, 0.29, 0.28, 0.27, 0.26, 0.25, &
         0.24, 0.23, 0.22, 0.21, 0.20, 0.19, 0.18, 0.17, 0.16, 0.15 /)

      soil_s_n_refl = (/ 0.52, 0.48, 0.44, 0.40, 0.38, 0.36, 0.34, 0.32, 0.30, 0.28, &
         0.26, 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.08 /)

      soil_d_n_refl = (/ 0.63, 0.59, 0.55, 0.51, 0.49, 0.47, 0.45, 0.43, 0.41, 0.39, &
         0.37, 0.35, 0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.19 /)

      write(cyear,'(i4.4)') lc_year
      landdir = trim(dir_model_landdata) // '/soil/' // trim(cyear)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,'(/, A29)') 'Aggregate Soil Brightness ...'
         CALL system('mkdir -p ' // trim(adjustl(landdir)))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

! -------------------------------------------------------------------------------------
! aggregate the soil parameters from the resolution of raw data to modelling resolution
! -------------------------------------------------------------------------------------

      lndname = trim(dir_rawdata)//'/soil_brightness.nc'

      IF (p_is_io) THEN

         CALL allocate_block_data (gland, isc)
         CALL allocate_block_data (gland, a_s_v_refl)
         CALL allocate_block_data (gland, a_d_v_refl)
         CALL allocate_block_data (gland, a_s_n_refl)
         CALL allocate_block_data (gland, a_d_n_refl)

         ! Read in the index of soil brightness (color)
         CALL ncio_read_block (lndname, 'soil_brightness', gland, isc)

         DO iblkme = 1, gblock%nblkme
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)

            DO iy = 1, gland%ycnt(jblk)
               DO ix = 1, gland%xcnt(iblk)

                  ii = isc%blk(iblk,jblk)%val(ix,iy)
                  IF ((ii >= 1) .and. (ii <= 20)) THEN
                     a_s_v_refl%blk(iblk,jblk)%val(ix,iy) = soil_s_v_refl( ii )
                     a_d_v_refl%blk(iblk,jblk)%val(ix,iy) = soil_d_v_refl( ii )
                     a_s_n_refl%blk(iblk,jblk)%val(ix,iy) = soil_s_n_refl( ii )
                     a_d_n_refl%blk(iblk,jblk)%val(ix,iy) = soil_d_n_refl( ii )
                  ELSE
                     a_s_v_refl%blk(iblk,jblk)%val(ix,iy) = spval
                     a_d_v_refl%blk(iblk,jblk)%val(ix,iy) = spval
                     a_s_n_refl%blk(iblk,jblk)%val(ix,iy) = spval
                     a_d_n_refl%blk(iblk,jblk)%val(ix,iy) = spval
                  ENDIF

               ENDDO
            ENDDO

         ENDDO

      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef USEMPI
      IF (p_is_io) THEN
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = a_s_v_refl)
      ENDIF
#endif

      IF (p_is_worker) THEN

         allocate ( soil_s_v_alb (numpatch) )

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)

            IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
               wmo_src = wmo_source (landpatch%ielm(ipatch))

               soil_s_v_alb(ipatch) = soil_s_v_alb(wmo_src)

               CYCLE
            ENDIF

#ifdef LULC_USGS
            IF(L/=16 .and. L/=24)THEN  ! NOT OCEAN(0)/WATER BODIES(16)/GLACIER and ICESHEET(24)
#else
            IF(L/=17 .and. L/=15)THEN  ! NOT OCEAN(0)/WATER BODIES(17)/GLACIER and ICE SHEET(15)
#endif
               CALL aggregation_request_data (landpatch, ipatch, gland, &
                  zip = USE_zip_for_aggregation, &
                  data_r8_2d_in1 = a_s_v_refl, data_r8_2d_out1 = soil_one)
               soil_s_v_alb (ipatch) = median (soil_one, size(soil_one), spval)

            ELSE
               soil_s_v_alb (ipatch) = spval
            ENDIF

         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef USEMPI
      IF (p_is_io) THEN
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = a_d_v_refl)
      ENDIF
#endif

      IF (p_is_worker) THEN

         allocate ( soil_d_v_alb (numpatch) )

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)

            IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
               wmo_src = wmo_source (landpatch%ielm(ipatch))

               soil_d_v_alb(ipatch) = soil_d_v_alb(wmo_src)

               CYCLE
            ENDIF

#ifdef LULC_USGS
            IF(L/=16 .and. L/=24)THEN  ! NOT OCEAN(0)/WATER BODIES(16)/GLACIER and ICESHEET(24)
#else
            IF(L/=17 .and. L/=15)THEN  ! NOT OCEAN(0)/WATER BODIES(17)/GLACIER and ICE SHEET(15)
#endif
               CALL aggregation_request_data (landpatch, ipatch, gland, &
                  zip = USE_zip_for_aggregation, &
                  data_r8_2d_in1 = a_d_v_refl, data_r8_2d_out1 = soil_one)
               soil_d_v_alb (ipatch) = median (soil_one, size(soil_one), spval)

            ELSE
               soil_d_v_alb (ipatch) = spval
            ENDIF

         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef USEMPI
      IF (p_is_io) THEN
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = a_s_n_refl)
      ENDIF
#endif

      IF (p_is_worker) THEN

         allocate ( soil_s_n_alb (numpatch) )

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)

            IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
               wmo_src = wmo_source (landpatch%ielm(ipatch))

               soil_s_n_alb(ipatch) = soil_s_n_alb(wmo_src)

               CYCLE
            ENDIF

#ifdef LULC_USGS
            IF(L/=16 .and. L/=24)THEN  ! NOT OCEAN(0)/WATER BODIES(16)/GLACIER and ICESHEET(24)
#else
            IF(L/=17 .and. L/=15)THEN  ! NOT OCEAN(0)/WATER BODIES(17)/GLACIER and ICE SHEET(15)
#endif
               CALL aggregation_request_data (landpatch, ipatch, gland, &
                  zip = USE_zip_for_aggregation, &
                  data_r8_2d_in1 = a_s_n_refl, data_r8_2d_out1 = soil_one)
               soil_s_n_alb (ipatch) = median (soil_one, size(soil_one), spval)

            ELSE
               soil_s_n_alb (ipatch) = spval
            ENDIF

         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef USEMPI
      IF (p_is_io) THEN
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = a_d_n_refl)
      ENDIF
#endif

      IF (p_is_worker) THEN

         allocate ( soil_d_n_alb (numpatch) )

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)

            IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
               wmo_src = wmo_source (landpatch%ielm(ipatch))

               soil_d_n_alb(ipatch) = soil_d_n_alb(wmo_src)

               CYCLE
            ENDIF

#ifdef LULC_USGS
            IF(L/=16 .and. L/=24)THEN  ! NOT OCEAN(0)/WATER BODIES(16)/GLACIER and ICESHEET(24)
#else
            IF(L/=17 .and. L/=15)THEN  ! NOT OCEAN(0)/WATER BODIES(17)/GLACIER and ICE SHEET(15)
#endif
               CALL aggregation_request_data (landpatch, ipatch, gland, &
                  zip = USE_zip_for_aggregation, &
                  data_r8_2d_in1 = a_d_n_refl, data_r8_2d_out1 = soil_one)
               soil_d_n_alb (ipatch) = median (soil_one, size(soil_one), spval)

            ELSE
               soil_d_n_alb (ipatch) = spval
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
      CALL check_vector_data ('s_v_alb ', soil_s_v_alb)
      CALL check_vector_data ('d_v_alb ', soil_d_v_alb)
      CALL check_vector_data ('s_n_alb ', soil_s_n_alb)
      CALL check_vector_data ('d_n_alb ', soil_d_n_alb)
#endif

      ! (1) Write-out the albedo of visible of the saturated soil
      lndname = trim(landdir)//'/soil_s_v_alb_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'soil_s_v_alb', 'patch', &
         landpatch, soil_s_v_alb, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typpatch = (/(ityp, ityp = 0, N_land_classification)/)
      lndname  = trim(dir_model_landdata) // '/diag/soil_brightness_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (soil_s_v_alb, landpatch%settyp, typpatch, m_patch2diag, &
         spval, lndname, 'soil_s_v_alb', compress = 1, write_mode = 'one', create_mode=.true.)
#endif

      ! (2) Write-out the albedo of visible of the dry soil
      lndname = trim(landdir)//'/soil_d_v_alb_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'soil_d_v_alb', 'patch', &
         landpatch, soil_d_v_alb, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typpatch = (/(ityp, ityp = 0, N_land_classification)/)
      lndname  = trim(dir_model_landdata) // '/diag/soil_brightness_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (soil_d_v_alb, landpatch%settyp, typpatch, m_patch2diag, &
         spval, lndname, 'soil_d_v_alb', compress = 1, write_mode = 'one')
#endif

      ! (3) Write-out the albedo of near infrared of the saturated soil
      lndname = trim(landdir)//'/soil_s_n_alb_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'soil_s_n_alb', 'patch', &
         landpatch, soil_s_n_alb, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typpatch = (/(ityp, ityp = 0, N_land_classification)/)
      lndname  = trim(dir_model_landdata) // '/diag/soil_brightness_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (soil_s_n_alb, landpatch%settyp, typpatch, m_patch2diag, &
         spval, lndname, 'soil_s_n_alb', compress = 1, write_mode = 'one')
#endif

      ! (4) Write-out the albedo of near infrared of the dry soil
      lndname = trim(landdir)//'/soil_d_n_alb_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'soil_d_n_alb', 'patch', &
         landpatch, soil_d_n_alb, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typpatch = (/(ityp, ityp = 0, N_land_classification)/)
      lndname  = trim(dir_model_landdata) // '/diag/soil_brightness_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (soil_d_n_alb, landpatch%settyp, typpatch, m_patch2diag, &
         spval, lndname, 'soil_d_n_alb', compress = 1, write_mode = 'one')
#endif

      ! Deallocate the allocatable array
      ! --------------------------------
      IF (p_is_worker) THEN
         deallocate ( soil_s_v_alb )
         deallocate ( soil_d_v_alb )
         deallocate ( soil_s_n_alb )
         deallocate ( soil_d_n_alb )
      ENDIF

END SUBROUTINE Aggregation_SoilBrightness
!-----------------------------------------------------------------------
!EOP
