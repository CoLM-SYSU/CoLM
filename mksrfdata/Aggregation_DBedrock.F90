#include <define.h>

SUBROUTINE Aggregation_DBedrock ( &
      gland, dir_rawdata, dir_model_landdata, lc_year)

!-----------------------------------------------------------------------
!  Depth to bedrock
!
!    Shangguan, W., Hengl, T., Mendes de Jesus, J., Yuan, H., Dai, Y. (2017).
!    Mapping the global depth to bedrock for land surface modeling.
!    Journal of Advances in Modeling Earth Systems, 9(1), 65-88.
!
!  Created by Shupeng Zhang, 05/2023
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_LandPatch
   USE MOD_Land2mWMO
   USE MOD_NetCDFVector
   USE MOD_NetCDFBlock
   USE MOD_RangeCheck
   USE MOD_AggregationRequestData

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
   character(len=256) :: landdir, lndname, cyear

   type (block_data_real8_2d) :: dbedrock
   real(r8), allocatable :: dbedrock_patches(:)
   real(r8), allocatable :: dbedrock_one(:), area_one(:)
   integer :: ipatch
   integer :: wmo_src

#ifdef SrfdataDiag
   integer :: typpatch(N_land_classification+1), ityp
#endif

      write(cyear,'(i4.4)') lc_year
      landdir = trim(dir_model_landdata) // '/dbedrock/' //trim(cyear)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,'(/, A)') 'Aggregate depth to bedrock ...'
         CALL system('mkdir -p ' // trim(adjustl(landdir)))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_io) THEN

         CALL allocate_block_data (gland, dbedrock)

         !TODO: add dir, file name and var name
         lndname = trim(dir_rawdata)//trim(DEF_rawdata%soil_bedrock%dir)//'/bedrock.nc'
         CALL ncio_read_block (lndname, 'dbedrock', gland, dbedrock)

#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = dbedrock)
#endif
      ENDIF

      IF (p_is_worker) THEN

         allocate (dbedrock_patches (numpatch))

         DO ipatch = 1, numpatch

            IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
               wmo_src = wmo_source (landpatch%ielm(ipatch))

               dbedrock_patches (ipatch) = dbedrock_patches (wmo_src)

               CYCLE
            ENDIF

            CALL aggregation_request_data (landpatch, ipatch, gland, &
               zip = USE_zip_for_aggregation, area = area_one, &
               data_r8_2d_in1 = dbedrock, data_r8_2d_out1 = dbedrock_one)
            dbedrock_patches (ipatch) = sum(dbedrock_one * area_one) / sum(area_one)
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('dbedrock_patches ', dbedrock_patches, -9999.0)
#endif

      ! Write-out the depth of the pacth in the gridcell
      lndname = trim(landdir)//'/dbedrock_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'dbedrock_patches', 'patch', &
           landpatch, dbedrock_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typpatch = (/(ityp, ityp = 0, N_land_classification)/)
      lndname  = trim(dir_model_landdata) // '/diag/dbedrock_patch_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (dbedrock_patches, landpatch%settyp, typpatch, m_patch2diag, &
         -1.0e36_r8, lndname, 'dbedrock', compress = 6, write_mode = 'one', create_mode=.true.)
#endif

      IF (p_is_worker) THEN
         deallocate ( dbedrock_patches )
      ENDIF

END SUBROUTINE Aggregation_DBedrock
