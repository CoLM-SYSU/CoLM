#include <define.h>

SUBROUTINE Aggregation_SoilTexture ( &
      gland, dir_rawdata, dir_model_landdata, lc_year)

!-----------------------------------------------------------------------
! !DESCRIPTION:
!  Aggregate soil texture class within a patch.
!
!  Use the USDA soil texture triangle (using the amount of sand, clay, and
!  silt contents) to identify the soil texture in fine grid resolution and
!  then finding the major soil type in a patch by counting number of fine
!  grids with each type of soil and adopting the major one.
!
!  Created by Shupeng Zhang, 01/2025
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
   USE MOD_Utils, only: num_max_frequency
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
   integer :: ipatch
   integer :: wmo_src

   type(block_data_int32_2d) :: soiltext
   integer, allocatable :: soiltext_patches(:), soiltext_one(:)
#ifdef SrfdataDiag
   integer :: typpatch(N_land_classification+1), ityp
   real(r8), allocatable :: soiltext_r8 (:)
#endif

      write(cyear,'(i4.4)') lc_year
      landdir = trim(dir_model_landdata) // '/soil/' // trim(cyear)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,'(/, A)') 'Aggregate soil texture ...'
         CALL system('mkdir -p ' // trim(adjustl(landdir)))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      lndname = trim(dir_rawdata)//trim(DEF_rawdata%soil_property%dir)//&
                '/soiltexture_0cm-60cm_mean.nc'

      IF (p_is_io) THEN
         CALL allocate_block_data (gland, soiltext)
         CALL ncio_read_block (lndname, 'soiltexture', gland, soiltext)

#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_i4_2d_in1 = soiltext)
#endif
      ENDIF

      IF (p_is_worker) THEN

         IF (numpatch > 0) allocate (soiltext_patches (numpatch))

         DO ipatch = 1, numpatch

            IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
               wmo_src = wmo_source (landpatch%ielm(ipatch))

               soiltext_patches(ipatch) = soiltext_patches(wmo_src)

               CYCLE
            ENDIF

            CALL aggregation_request_data (landpatch, ipatch, gland, &
               zip = USE_zip_for_aggregation, &
               data_i4_2d_in1 = soiltext, data_i4_2d_out1 = soiltext_one)
            soiltext_patches(ipatch) = num_max_frequency (soiltext_one)
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('soiltext_patches ', soiltext_patches)
#endif

      lndname = trim(landdir)//'/soiltexture_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'soiltext_patches', 'patch', &
         landpatch, soiltext_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typpatch = (/(ityp, ityp = 0, N_land_classification)/)
      lndname = trim(dir_model_landdata)//'/diag/soiltexture_'//trim(cyear)//'.nc'
      IF (allocated(soiltext_patches)) THEN
         allocate (soiltext_r8 (size(soiltext_patches)))
         soiltext_r8 = real(soiltext_patches,r8)
      ENDIF

      CALL srfdata_map_and_write (soiltext_r8, landpatch%settyp, typpatch,  &
         m_patch2diag, -1., lndname, 'soiltexture', compress = 6, write_mode = 'one', &
         create_mode=.true.)

      IF (allocated(soiltext_r8)) deallocate(soiltext_r8)
#endif

      IF (p_is_worker) THEN
         IF (numpatch > 0) deallocate ( soiltext_patches )
      ENDIF

END SUBROUTINE Aggregation_SoilTexture
