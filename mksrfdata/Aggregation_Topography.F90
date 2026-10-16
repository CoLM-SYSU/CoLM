#include <define.h>

SUBROUTINE Aggregation_Topography ( &
      gtopo, dir_rawdata, dir_model_landdata, lc_year)
!-----------------------------------------------------------------------
!  Global Topography data
!
!   Yamazaki, D., Ikeshima, D., Sosa, J.,Bates, P. D., Allen, G. H.,
!   Pavelsky, T. M. (2019).
!   MERIT Hydro: a high‐resolution global hydrographymap based on
!   latest topography dataset.Water Resources Research, 55, 5053-5073.
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
   integer, intent(in) :: lc_year
   type(grid_type),  intent(in) :: gtopo
   character(len=*), intent(in) :: dir_rawdata
   character(len=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ---------------------------------------------------------------
   character(len=256) :: landdir, lndname, cyear
   integer  :: ipatch, i, ps, pe
   integer  :: wmo_src
   real(r8) :: sumarea

   type (block_data_real8_2d) :: landarea
   type (block_data_real8_2d) :: elevation
   type (block_data_real8_2d) :: elvstd
   type (block_data_real8_2d) :: sloperatio

   real(r8), allocatable :: elevation_one    (:), elvstd_one    (:), sloperatio_one    (:)
   real(r8), allocatable :: elevation_patches(:), elvstd_patches(:), sloperatio_patches(:)
   real(r8), allocatable :: landarea_one     (:)

#ifdef SrfdataDiag
   integer :: typpatch(N_land_classification+1), ityp
#endif

      write(cyear,'(i4.4)') lc_year
      landdir = trim(dir_model_landdata) // '/topography/' // trim(cyear)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,'(/, A)') 'Aggregate topography ...'
         CALL system('mkdir -p ' // trim(adjustl(landdir)))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      lndname = trim(dir_rawdata)//trim(DEF_rawdata%topo%dir)//'/topography.nc'

! ---------------------------------------------------------------------------------
!   aggregate the elevation from the resolution of raw data to modelling resolution
! ---------------------------------------------------------------------------------

      IF (p_is_io) THEN

         CALL allocate_block_data (gtopo, landarea)
         CALL ncio_read_block (lndname, 'landarea', gtopo, landarea)

         CALL allocate_block_data (gtopo, elevation)
         CALL ncio_read_block (lndname, 'elevation', gtopo, elevation)

         CALL allocate_block_data (gtopo, elvstd)
         CALL ncio_read_block (lndname, 'elvstd', gtopo, elvstd)

         CALL allocate_block_data (gtopo, sloperatio)
         CALL ncio_read_block (lndname, 'slope', gtopo, sloperatio)

#ifdef USEMPI
         CALL aggregation_data_daemon (gtopo,                      &
            data_r8_2d_in1 = landarea, data_r8_2d_in2 = elevation, &
            data_r8_2d_in3 = elvstd  , data_r8_2d_in4 = sloperatio )
#endif
      ENDIF

      IF (p_is_worker) THEN

         allocate (elevation_patches  (numpatch))
         allocate (elvstd_patches     (numpatch))
         allocate (sloperatio_patches (numpatch))

         DO ipatch = 1, numpatch

            IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
               wmo_src = wmo_source (landpatch%ielm(ipatch))

               elevation_patches (ipatch) = elevation_patches (wmo_src)
               elvstd_patches    (ipatch) = elvstd_patches    (wmo_src)
               sloperatio_patches(ipatch) = sloperatio_patches(wmo_src)

               CYCLE
            ENDIF

            CALL aggregation_request_data (landpatch, ipatch, gtopo,          &
               zip = USE_zip_for_aggregation,                                 &
               data_r8_2d_in1 = landarea  , data_r8_2d_out1 = landarea_one  , &
               data_r8_2d_in2 = elevation , data_r8_2d_out2 = elevation_one , &
               data_r8_2d_in3 = elvstd    , data_r8_2d_out3 = elvstd_one    , &
               data_r8_2d_in4 = sloperatio, data_r8_2d_out4 = sloperatio_one  )

            IF (any(elevation_one /= -9999.0)) THEN

               sumarea = sum(landarea_one, mask = elevation_one /= -9999.0)

               elevation_patches(ipatch) = &
                  sum(elevation_one * landarea_one, mask = elevation_one /= -9999.0) / sumarea

               elvstd_patches(ipatch) = sqrt( &
                  sum(((elevation_one - elevation_patches(ipatch))**2 + elvstd_one**2) * landarea_one, &
                  mask = elevation_one /= -9999.0) / sumarea)

               sloperatio_patches(ipatch) = &
                  sum(sloperatio_one * landarea_one, mask = elevation_one /= -9999.0) / sumarea

            ELSE
               elevation_patches (ipatch) = 0.
               elvstd_patches    (ipatch) = 0.
               sloperatio_patches(ipatch) = 0.
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
      CALL check_vector_data ('elevation_patches ', elevation_patches )
      CALL check_vector_data ('elvstd_patches    ', elvstd_patches    )
      CALL check_vector_data ('sloperatio_patches', sloperatio_patches)
#endif

      lndname = trim(landdir)//'/elevation_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'elevation_patches', 'patch', landpatch, &
         elevation_patches, DEF_Srfdata_CompressLevel)

      lndname = trim(landdir)//'/elvstd_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'elvstd_patches', 'patch', landpatch, &
         elvstd_patches, DEF_Srfdata_CompressLevel)

      lndname = trim(landdir)//'/sloperatio_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'sloperatio_patches', 'patch', landpatch, &
         sloperatio_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typpatch = (/(ityp, ityp = 0, N_land_classification)/)
      lndname  = trim(dir_model_landdata) // '/diag/topography_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (elevation_patches, landpatch%settyp, typpatch, m_patch2diag, &
         -1.0e36_r8, lndname, 'elevation', compress = 6, write_mode = 'one', create_mode=.true.)

      typpatch = (/(ityp, ityp = 0, N_land_classification)/)
      lndname  = trim(dir_model_landdata) // '/diag/topography_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (elvstd_patches, landpatch%settyp, typpatch, m_patch2diag, &
         -1.0e36_r8, lndname, 'elvstd', compress = 6, write_mode = 'one')

      typpatch = (/(ityp, ityp = 0, N_land_classification)/)
      lndname  = trim(dir_model_landdata) // '/diag/topography_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (sloperatio_patches, landpatch%settyp, typpatch, m_patch2diag, &
         -1.0e36_r8, lndname, 'sloperatio', compress = 6, write_mode = 'one')
#endif

      IF (p_is_worker) THEN
         deallocate ( elevation_patches  )
         deallocate ( elvstd_patches     )
         deallocate ( sloperatio_patches )
      ENDIF

END SUBROUTINE Aggregation_Topography
