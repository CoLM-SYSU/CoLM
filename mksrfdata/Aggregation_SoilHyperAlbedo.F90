#include <define.h>

SUBROUTINE Aggregation_SoilHyperAlbedo ( &
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

   IMPLICIT NONE

   ! arguments:

   integer, intent(in) :: lc_year
   type(grid_type),  intent(in) :: gland
   character(len=*), intent(in) :: dir_rawdata
   character(len=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ---------------------------------------------------------------
   character(len=256) :: landdir, lndname, cyear
   character(len=4) :: wavelength  ! wavelength in nm

   ! type (block_data_real8_2d) :: a_s_v_refl
   type (block_data_real8_2d) :: a_soil_hyper_alb

   ! real(r8), allocatable :: soil_s_v_alb (:)
   real(r8), allocatable :: soil_hyper_alb (:)

   integer :: ii, L
   integer :: ipatch, iblkme, iblk, jblk, ix, iy, wl
   real(r8), allocatable :: soil_one(:)

#ifdef SrfdataDiag
   integer :: typpatch(N_land_classification+1), ityp
#endif

   write(cyear,'(i4.4)') lc_year
   landdir = trim(dir_model_landdata) // '/HyperAlbedo/' // trim(cyear)

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif
   IF (p_is_master) THEN
      write(*,'(/, A29)') 'Aggregate Hyper Albedo ...'
      CALL system('mkdir -p ' // trim(adjustl(landdir)))
   ENDIF
#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

   ! 为 worker 在循环外分配一次，避免每次循环重复分配引发错误
   IF (p_is_worker) THEN
      IF (.NOT. ALLOCATED(soil_hyper_alb)) ALLOCATE(soil_hyper_alb(numpatch))
   ENDIF

   DO wl = 400, 2500, 10
      IF (p_is_master) write(*,'(A, i0, A)') '  Processing wavelength (nm): ', wl, ' ...'
      write(wavelength,'(i0)') wl

      ! Read in the hyper albedo
      IF (p_is_io) THEN

         CALL allocate_block_data (gland, a_soil_hyper_alb)
         lndname = trim(dir_rawdata)//'/colm_input_ghsad/colm_soil_albedo_'//trim(wavelength)//'nm.nc'
         ! Read in the soil albedo
         CALL ncio_read_block (lndname, 'albedo', gland, a_soil_hyper_alb)

         DO iblkme = 1, gblock%nblkme
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)
            DO iy = 1, gland%ycnt(jblk)
               DO ix = 1, gland%xcnt(iblk)
                  a_soil_hyper_alb%blk(iblk,jblk)%val(ix,iy) = a_soil_hyper_alb%blk(iblk,jblk)%val(ix,iy) / 10000.0
               END DO
            END DO
         END DO


#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = a_soil_hyper_alb)
#endif

      ENDIF

      IF (p_is_worker) THEN

         ! allocate 已在循环外完成，移除这里的重复分配
         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)
#ifdef LULC_USGS
            IF(L/=16 .and. L/=24)THEN  ! NOT OCEAN(0)/WATER BODIES(16)/GLACIER and ICESHEET(24)
#else
            IF(L/=17 .and. L/=15)THEN  ! NOT OCEAN(0)/WATER BODIES(17)/GLACIER and ICE SHEET(15)
#endif
               CALL aggregation_request_data (landpatch, ipatch, gland, &
                  zip = USE_zip_for_aggregation, &
                  data_r8_2d_in1 = a_soil_hyper_alb, data_r8_2d_out1 = soil_one)
               soil_hyper_alb (ipatch) = median (soil_one, size(soil_one))

            ELSE
               soil_hyper_alb (ipatch) = -1.0e36_r8
            ENDIF

         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF    ! IF (p_is_worker)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('soil_hyper_alb ', soil_hyper_alb)
#endif

      ! Write-out the hyper spectral soil albedo
      lndname = trim(landdir)//'/soil_hyper_alb_'//trim(wavelength)//'nm_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'soil_hyper_alb', 'patch', &
         landpatch, soil_hyper_alb, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typpatch = (/(ityp, ityp = 0, N_land_classification)/)
      lndname  = trim(dir_model_landdata) // '/diag/soil_hyper_alb_' // trim(wavelength) // 'nm_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (soil_hyper_alb, landpatch%settyp, typpatch, m_patch2diag, &
         -1.0e36_r8, lndname, 'soil_hyper_alb', compress = 1, write_mode = 'one')
#endif

   END DO   ! loop over "wl"

   ! Deallocate the allocatable array
   ! --------------------------------
   IF (p_is_worker) THEN
      deallocate ( soil_hyper_alb )
   ENDIF

END SUBROUTINE Aggregation_SoilHyperAlbedo
!-----------------------------------------------------------------------
!EOP
