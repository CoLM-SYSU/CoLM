#include <define.h>

SUBROUTINE Aggregation_LakeDepth ( &
      gland, dir_rawdata, dir_model_landdata, lc_year)

!-----------------------------------------------------------------------
! !DESCRIPTION:
!  Aggregate lake depth of multiple pixels within a lake patch based on
!  Global land cover types (updated with the specific dataset)
!
!  Global Lake Coverage and Lake Depth (1km resolution)
!    (http://nwpi.krc.karelia.run/flake/)
!     Lake depth data legend
!     Value   Description
!     0       no lake indicated in this pixel
!     1       no any information about this lake and set the default value of 10 m
!     2       no information about depth for this lake and set the default value of 10 m
!     3       have the information about lake depth in this pixel
!     4       this is the river pixel according to our map, set the default value of 3 m
!
! !REFERENCES:
!  Kourzeneva, E., H. Asensio, E. Martin, and S. Faroux, 2012: Global gridded
!  dataset of lake coverage and lake depth for USE in numerical weather
!  prediction and climate modelling. Tellus A, 64, 15640.
!
!  Created by Yongjiu Dai, 02/2014
!
! !REVISIONS:
!  Shupeng Zhang, 01/2022: porting codes to MPI parallel version
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_LandPatch
   USE MOD_NetCDFVector
   USE MOD_NetCDFBlock
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
   integer :: L, ipatch

   type (block_data_real8_2d) :: lakedepth
   real(r8), allocatable :: lakedepth_patches(:), lakedepth_one(:)
#ifdef SrfdataDiag
   integer :: typlake(1) = (/17/)
#endif

      write(cyear,'(i4.4)') lc_year
      landdir = trim(dir_model_landdata) // '/lakedepth/' // trim(cyear)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,'(/, A)') 'Aggregate lake depth ...'
         CALL system('mkdir -p ' // trim(adjustl(landdir)))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

! ................................................
! global lake coverage and lake depth
! ................................................
      !TODO:
      lndname = trim(dir_rawdata)//trim(DEF_rawdata%hydro%dir)//'/lake_depth.nc'

      IF (p_is_io) THEN
         CALL allocate_block_data (gland, lakedepth)
         CALL ncio_read_block (lndname, 'lake_depth', gland, lakedepth)
         CALL block_data_linear_transform (lakedepth, scl = 0.1)

#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = lakedepth)
#endif
      ENDIF

! ----------------------------------------------------------------------------------
!   aggregate the lake depth from the resolution of raw data to modelling resolution
! ----------------------------------------------------------------------------------

      IF (p_is_worker) THEN

         allocate (lakedepth_patches (numpatch))

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)
            IF(L==WATERBODY)THEN  ! LAND WATER BODIES (17)
               CALL aggregation_request_data (landpatch, ipatch, gland, &
                  zip = USE_zip_for_aggregation, &
                  data_r8_2d_in1 = lakedepth, data_r8_2d_out1 = lakedepth_one)
               lakedepth_patches (ipatch) = median (lakedepth_one, size(lakedepth_one))
            ELSE
               lakedepth_patches (ipatch) = -1.0e36_r8
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
      CALL check_vector_data ('lakedepth_patches ', lakedepth_patches)
#endif

      lndname = trim(landdir)//'/lakedepth_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'lakedepth_patches', 'patch', &
         landpatch, lakedepth_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      lndname = trim(dir_model_landdata)//'/diag/lakedepth_'//trim(cyear)//'.nc'
      CALL srfdata_map_and_write (lakedepth_patches, landpatch%settyp, typlake, m_patch2diag, &
         -1.0e36_r8, lndname, 'lakedepth', compress = 6, write_mode = 'one', create_mode=.true.)
#endif

      IF (p_is_worker) THEN
         deallocate ( lakedepth_patches )
      ENDIF

END SUBROUTINE Aggregation_LakeDepth
