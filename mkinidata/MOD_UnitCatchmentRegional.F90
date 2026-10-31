#include <define.h>

MODULE MOD_UnitCatchmentRegional

!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!    Cut the unit-catchment routing network down to the river systems that receive
!    runoff from the land domain (DEF_UnitCatchment_regional = .true.).
!
!    GridRiverLakeFlow routes every unit catchment of the network file, wet or dry.
!    For a regional domain most of a global network never receives water, so the
!    run pays for it without using it.  This step finds the runoff-input grid cells
!    that the land patches overlap, in the same way GridRiverLakeFlow does at run
!    time (build_worker_remapdata on landpatch and the input grid), and lets
!    MOD_UnitCatchmentSubset write the regional network to
!    <DEF_dir_output>/<case>/riverlake/unitcatchment_regional.nc, which the model reads
!    instead of DEF_UnitCatchment_file (see get_unitcatchment_file).
!
!    Called by mkinidata before the land data of the initial year are loaded, so a
!    new network (or a new DEF_UnitCatchment_file) needs no new surface data.  The
!    file goes to the case directory, not to landdata, so landdata shared between
!    cases is never written.
!
!    With LULCC the model keeps one network for all years (its state is carried
!    across the yearly re-initialisation), so the cells of the land patches of every
!    land-cover year from the start year to the end year of the simulation are used.
!    A year whose land data are not there (for example an end year that the run only
!    reaches at its last instant) is skipped with a note; the model checks at start-up
!    and after each LULCC that the network covers the land patches of the year
!    (check_regional_coverage in MOD_Grid_RiverLakeNetwork).
!-------------------------------------------------------------------------------------

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: unitcatchment_regional_build

CONTAINS

   SUBROUTINE unitcatchment_regional_build (dir_landdata, year_first, year_last)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Grid
   USE MOD_Mesh,       only: mesh_free_mem, numelm
   USE MOD_LandElm,    only: landelm
   USE MOD_LandPatch,  only: landpatch, numpatch
   USE MOD_SrfdataRestart, only: mesh_load_from_file, pixelset_load_from_file
   USE MOD_WorkerPushData
   USE MOD_UnitCatchmentSubset
   IMPLICIT NONE

   character(len=*), intent(in) :: dir_landdata
   integer,          intent(in) :: year_first, year_last

   type(grid_type)             :: gridro
   type(worker_remapdata_type) :: remap
   logical, allocatable :: touched(:,:)
   integer, allocatable :: ids(:)
   integer :: nlon, nlat, ng, iworker, i, ix, iy, nkeep, nsystem, year, nyear
   logical :: have_year
   character(len=256) :: file_regional
   character(len=4)   :: cyear

      file_regional = regional_unitcatchment_file ()

      nlon = 0
      nlat = 0
      IF (p_is_master) THEN
         CALL unitcatchment_grid_size (trim(DEF_UnitCatchment_file), nlon, nlat)
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast (nlon, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (nlat, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif

      CALL gridro%define_by_ndims (nlon, nlat)

      allocate (touched (nlon, nlat))
      touched = .false.
      nyear = 0

      DO year = year_first, year_last

         ! The land data of the first year must be there (mkinidata loads them next);
         ! a later LULCC year may be missing.
         write(cyear,'(i4.4)') year
         IF (p_is_master) THEN
            inquire (file = trim(dir_landdata)//'/mesh/'//cyear//'/mesh.nc', exist = have_year)
         ENDIF
#ifdef USEMPI
         CALL mpi_bcast (have_year, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif
         IF (.not. have_year) THEN
            IF (year == year_first) CALL CoLM_stop ('Regional unit catchment: no land data for year ' &
               // cyear // ' in ' // trim(dir_landdata))
            IF (p_is_master) write(*,'(3A)') ' Regional unit-catchment network: no land data for ', &
               cyear, ', year skipped'
            CYCLE
         ENDIF

         ! Load the land patches of this year as LulccInitialize does, and free them
         ! again: mkinidata loads the initial year itself afterwards.
         CALL mesh_load_from_file     (dir_landdata, year)
         CALL pixelset_load_from_file (dir_landdata, 'landelm'  , landelm  , numelm  , year)
         CALL pixelset_load_from_file (dir_landdata, 'landpatch', landpatch, numpatch, year)

         ! Runoff-input cells overlapped by the land patches: the cells
         ! GridRiverLakeFlow will pull runoff from.
         CALL build_worker_remapdata (landpatch, gridro, remap)

#ifdef USEMPI
         IF (p_is_master) THEN
            DO iworker = 0, p_np_worker-1
               CALL mpi_recv (ng, 1, MPI_INTEGER, p_address_worker(iworker), mpi_tag_mesg, &
                  p_comm_glb, p_stat, p_err)
               IF (ng > 0) THEN
                  allocate (ids (ng))
                  CALL mpi_recv (ids, ng, MPI_INTEGER, p_address_worker(iworker), mpi_tag_data, &
                     p_comm_glb, p_stat, p_err)
                  DO i = 1, ng
                     ix = mod(ids(i)-1, nlon) + 1
                     iy = (ids(i)-1) / nlon + 1
                     touched(ix,iy) = .true.
                  ENDDO
                  deallocate (ids)
               ENDIF
            ENDDO
         ELSEIF (p_is_worker) THEN
            ng = remap%num_grid
            CALL mpi_send (ng, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
            IF (ng > 0) THEN
               CALL mpi_send (remap%ids_me, ng, MPI_INTEGER, p_address_master, mpi_tag_data, &
                  p_comm_glb, p_err)
            ENDIF
         ENDIF
#else
         ng = remap%num_grid
         DO i = 1, ng
            ix = mod(remap%ids_me(i)-1, nlon) + 1
            iy = (remap%ids_me(i)-1) / nlon + 1
            touched(ix,iy) = .true.
         ENDDO
#endif

         CALL worker_remapdata_free_mem (remap)
         CALL mesh_free_mem
         CALL landelm%forc_free_mem
         CALL landpatch%forc_free_mem
         nyear = nyear + 1

      ENDDO

      IF (p_is_master) THEN
         CALL execute_command_line ('mkdir -p ' // trim(file_regional(1:index(file_regional, '/', back=.true.))))
         CALL unitcatchment_subset_write (trim(DEF_UnitCatchment_file), trim(file_regional), nlon, nlat, &
            touched, .true., nkeep, nsystem)
         write(*,'(/,A,I0,A,I0,A,I0,A)') ' Regional unit-catchment network: ', nkeep, ' unit catchments in ', &
            nsystem, ' river systems (land patches of ', nyear, ' land-cover year(s))'
         write(*,'(2A)') '   written to ', trim(file_regional)
      ENDIF

      CALL grid_free_mem (gridro)
      deallocate (touched)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE unitcatchment_regional_build

END MODULE MOD_UnitCatchmentRegional
