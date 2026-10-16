#include <define.h>

SUBROUTINE Aggregation_TopographyFactors_Simple ( &
      grid_topo_factor , dir_topodata, dir_model_landdata, lc_year)
!-----------------------------------------------------------------------
!  Global topography-based factors data
!
!  Created by Sisi Chen, 08/2025
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
   integer  :: wmo_src

   type (block_data_real8_3d) :: slp_grid    ! slope
   type (block_data_real8_3d) :: asp_grid    ! aspect
   type (block_data_real8_2d) :: cur_grid    ! curvature

   ! patch
   real(r8), allocatable :: cur_patches (:)

   ! nine defined types at all patches
   real(r8), allocatable :: asp_type_patches  (:,:) ! shape as (type, patches)
   real(r8), allocatable :: slp_type_patches  (:,:)

   ! pixelsets
   real(r8), allocatable :: cur_one         (:)
   real(r8), allocatable :: area_one        (:)

   ! pixelsets of nine defined types at each patch
   real(r8), allocatable :: slp_one         (:,:)
   real(r8), allocatable :: asp_one         (:,:)

   ! local variables
   integer :: ipatch, i

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
      lndname = trim(dir_topodata)//"/topography_MERITHydro.nc"
      CALL allocate_block_data (grid_topo_factor, slp_grid, num_aspect_type)
      CALL ncio_read_block (lndname, 'slp_aspect', grid_topo_factor, num_aspect_type, slp_grid)

      lndname = trim(dir_topodata)//"/topography_MERITHydro.nc"
      CALL allocate_block_data (grid_topo_factor, asp_grid, num_aspect_type)
      CALL ncio_read_block (lndname, 'pct_aspect', grid_topo_factor, num_aspect_type, asp_grid)

      lndname = trim(dir_topodata)//"/curvature_MERITHydro.nc"
      CALL allocate_block_data (grid_topo_factor, cur_grid)
      CALL ncio_read_block (lndname, 'curvature', grid_topo_factor, cur_grid)

   ! --------------------------------------------------------------------------
   ! aggregate the terrain factor data from the resolution of raw data to patch
   ! --------------------------------------------------------------------------
#ifdef USEMPI
      ! mpi send
      CALL aggregation_data_daemon (  grid_topo_factor, &
         data_r8_2d_in1 = cur_grid,  &
         data_r8_3d_in1 = slp_grid,   n1_r8_3d_in1 = num_aspect_type, &
         data_r8_3d_in2 = asp_grid,   n1_r8_3d_in2 = num_aspect_type)
#endif
   ENDIF


   IF (p_is_worker) THEN
      ! allocate for output variables at patches
      allocate (cur_patches      (numpatch))
      allocate (asp_type_patches (num_aspect_type, numpatch))
      allocate (slp_type_patches (num_aspect_type, numpatch))

      ! aggregate loop
      DO ipatch = 1, numpatch

         IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
            wmo_src = wmo_source (landpatch%ielm(ipatch))

            cur_patches     (  ipatch) = cur_patches     (  wmo_src)
            asp_type_patches(:,ipatch) = asp_type_patches(:,wmo_src)
            slp_type_patches(:,ipatch) = slp_type_patches(:,wmo_src)

            CYCLE
         ENDIF

         CALL aggregation_request_data (landpatch, ipatch, grid_topo_factor, &
            zip = USE_zip_for_aggregation, area = area_one, &
            data_r8_2d_in1 = cur_grid,   data_r8_2d_out1 = cur_one, &
            data_r8_3d_in1 = slp_grid,   data_r8_3d_out1 = slp_one,   n1_r8_3d_in1 = num_aspect_type, &
            data_r8_3d_in2 = asp_grid,   data_r8_3d_out2 = asp_one,   n1_r8_3d_in2 = num_aspect_type)

         ! ------------------------------------------------------------------
         ! aggregate curvature at patches
         ! ------------------------------------------------------------------

         IF (any(cur_one /= -9999.0)) THEN
            cur_patches (ipatch) = &
               sum(cur_one * area_one, mask = cur_one /= -9999.0) &
                  / sum(area_one, mask = cur_one /= -9999.0)
         ELSE
            cur_patches (ipatch) = -1.0e36
         ENDIF

         ! ------------------------------------------------------------------
         ! aggregate slope and aspect at each direction of aspect of patches.
         ! num_aspect_type = 1:north, 2:northeast, 3:east, 4:southeast, 5:south, 6:southwest, 7:west, 8:northwest, 9:flat
         ! ------------------------------------------------------------------
         DO i = 1, num_aspect_type
            IF (any(asp_one(i,:) /= -9999.0)) THEN
               asp_type_patches (i, ipatch) = &
                  sum(asp_one(i,:) * area_one(:), mask = asp_one(i,:) /= -9999.0) &
                     / sum(area_one, mask = asp_one(i,:) /= -9999.0)
            ELSE
               asp_type_patches (i, ipatch) = -1.0e36
            ENDIF
            IF (any(slp_one(i,:) /= -9999.0)) THEN
               slp_type_patches (i, ipatch) = &
                  sum(slp_one(i,:) * area_one(:), mask = slp_one(i,:) /= -9999.0) &
                     / sum(area_one, mask = slp_one(i,:) /= -9999.0)
            ELSE
               slp_type_patches (i, ipatch) = -1.0e36
            ENDIF

         ENDDO

      ENDDO


#ifdef USEMPI
      CALL aggregation_worker_done ()
#endif
   ENDIF

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
   CALL check_vector_data ('cur_patches       ', cur_patches      )
   CALL check_vector_data ('slp_type_patches  ', slp_type_patches )
   CALL check_vector_data ('asp_type_patches  ', asp_type_patches )
#endif

   ! --------------------------------------------------------------------------
   ! write aggregated data to netcdf files
   ! --------------------------------------------------------------------------

   lndname = trim(landdir)//'/cur_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
   CALL ncio_write_vector (lndname, 'cur_patches', 'patch', landpatch, cur_patches, 1)

   lndname = trim(landdir)//'/slp_type_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
   CALL ncio_define_dimension_vector (lndname, landpatch, 'slope_type', num_aspect_type)
   CALL ncio_write_vector (lndname, 'slp_type_patches', 'slope_type', num_aspect_type, 'patch', landpatch, slp_type_patches, 1)

   lndname = trim(landdir)//'/asp_type_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
   CALL ncio_define_dimension_vector (lndname, landpatch, 'slope_type', num_aspect_type)
   CALL ncio_write_vector (lndname, 'asp_type_patches', 'slope_type', num_aspect_type, 'patch', landpatch, asp_type_patches, 1)

   ! --------------------------------------------------------------------------
#ifdef SrfdataDiag
   typpatch = (/(ityp, ityp = 0, N_land_classification)/)

   ! only write the first type of slope and aspect at patches
   lndname  = trim(dir_model_landdata) // '/diag/topo_factor_slp_' // trim(cyear) // '.nc'
   DO i = 1, num_aspect_type
      write(sdir,'(I0)') i
      CALL srfdata_map_and_write (slp_type_patches(i,:), landpatch%settyp, typpatch, m_patch2diag, &
         -1.0e36_r8, lndname, 'slp_'//trim(sdir), compress = 6, write_mode = 'one')
   ENDDO

   lndname  = trim(dir_model_landdata) // '/diag/topo_factor_asp_' // trim(cyear) // '.nc'
   DO i = 1, num_aspect_type
      write(sdir,'(I0)') i
      CALL srfdata_map_and_write (asp_type_patches(i,:), landpatch%settyp, typpatch, m_patch2diag, &
         -1.0e36_r8, lndname, 'asp_'//trim(sdir), compress = 6, write_mode = 'one')
   ENDDO

   ! write curvature at patches
   lndname  = trim(dir_model_landdata) // '/diag/topo_factor_cur_' // trim(cyear) // '.nc'
   CALL srfdata_map_and_write (cur_patches, landpatch%settyp, typpatch, m_patch2diag, &
      -1.0e36_r8, lndname, 'cur', compress = 6, write_mode = 'one')

#endif

   IF (p_is_worker) THEN
      IF (allocated(slp_type_patches)) deallocate ( slp_type_patches )
      IF (allocated(asp_type_patches)) deallocate ( asp_type_patches )
      IF (allocated(cur_patches     )) deallocate ( cur_patches      )
   ENDIF

END SUBROUTINE Aggregation_TopographyFactors_Simple
