#include <define.h>

MODULE MOD_MeshFilter

!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!    Mesh filter.
!    Mesh filter can be used to mask part of region or globe as needed.
!
!  Created by Shupeng Zhang, May 2023
!-----------------------------------------------------------------------

   USE MOD_Grid
   IMPLICIT NONE

   logical :: has_mesh_filter
   type(grid_type) :: grid_filter

CONTAINS

   logical FUNCTION inquire_mesh_filter ()

   USE MOD_SPMD_Task
   USE MOD_Namelist
   IMPLICIT NONE
   logical :: fexists

      IF (p_is_master) THEN

         inquire (file=trim(DEF_file_mesh_filter), exist=fexists)

         IF (.not. fexists) THEN
            write(*,'(/, 2A)') 'Note: Mesh Filter not used: file ', trim(DEF_file_mesh_filter)
         ELSE
            write(*,'(/, 2A)') 'Note: Mesh Filter from file ', trim(DEF_file_mesh_filter)
         ENDIF
      ENDIF

#ifdef USEMPI
      CALL mpi_bcast (fexists, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif

      inquire_mesh_filter = fexists

   END FUNCTION inquire_mesh_filter

   ! -------------
   SUBROUTINE mesh_filter (gridf, ffilter, fvname)

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_NetCDFBlock
   USE MOD_DataType
   USE MOD_LandElm
   USE MOD_Mesh
   USE MOD_AggregationRequestData
   USE MOD_Block
   IMPLICIT NONE

   type(grid_type),  intent(in) :: gridf
   character(len=*), intent(in) :: ffilter
   character(len=*), intent(in) :: fvname

   ! local variables:
   ! ---------------------------------------------------------------
   type (block_data_int32_2d) :: datafilter
   integer, allocatable :: ifilter(:), xtemp(:), ytemp(:)
   logical, allocatable :: filter(:)
   integer :: ielm, jelm, npxl, nelm_glb

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,'(/, A)') 'Filtering pixels ...'
      ENDIF

      IF (p_is_io) THEN
         CALL allocate_block_data (gridf, datafilter)
         CALL ncio_read_block (trim(ffilter), trim(fvname), gridf, datafilter)

#ifdef USEMPI
         CALL aggregation_data_daemon (gridf, data_i4_2d_in1 = datafilter)
#endif
      ENDIF

      IF (p_is_worker) THEN

         jelm = 0
         DO ielm = 1, numelm
            CALL aggregation_request_data (landelm, ielm, gridf, zip = .false., &
               data_i4_2d_in1 = datafilter, data_i4_2d_out1 = ifilter, &
               filledvalue_i4 = -1)

            allocate (filter (mesh(ielm)%npxl))
            filter = ifilter > 0

            IF (any(filter)) THEN
               jelm = jelm + 1
               IF (.not. all(filter)) THEN

                  npxl = count(filter)

                  allocate (xtemp(npxl))
                  allocate (ytemp(npxl))
                  xtemp = pack(mesh(ielm)%ilon, filter)
                  ytemp = pack(mesh(ielm)%ilat, filter)

                  deallocate(mesh(ielm)%ilon)
                  deallocate(mesh(ielm)%ilat)

                  mesh(ielm)%npxl = npxl

                  allocate(mesh(ielm)%ilon(npxl))
                  allocate(mesh(ielm)%ilat(npxl))
                  mesh(ielm)%ilon = xtemp
                  mesh(ielm)%ilat = ytemp

                  deallocate (xtemp)
                  deallocate (ytemp)
               ENDIF

               IF (jelm /= ielm) THEN
                  CALL copy_elm (mesh(ielm), mesh(jelm))
               ENDIF

            ENDIF

            deallocate (filter)
         ENDDO

         numelm = jelm

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

      IF (p_is_worker) THEN
         IF (allocated(landelm%eindex))   deallocate (landelm%eindex)
         IF (allocated(landelm%ipxstt))   deallocate (landelm%ipxstt)
         IF (allocated(landelm%ipxend))   deallocate (landelm%ipxend)
         IF (allocated(landelm%settyp))   deallocate (landelm%settyp)
         IF (allocated(landelm%ielm  ))   deallocate (landelm%ielm  )
      ENDIF

      CALL landelm_build ()

#ifdef USEMPI
      IF (p_is_worker) THEN
         CALL mpi_reduce (numelm, nelm_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', nelm_glb, ' elements after mesh filtering.'
         ENDIF
      ENDIF
#else
      write(*,'(A,I12,A)') 'Total: ', numelm, ' elements after mesh filtering.'
#endif

      ! Update nelm_blk
      nelm_blk(:,:) = 0
      IF (p_is_worker) THEN
         DO ielm = 1, numelm
            nelm_blk(mesh(ielm)%xblk,mesh(ielm)%yblk) = &
               nelm_blk(mesh(ielm)%xblk,mesh(ielm)%yblk) + 1
         ENDDO
      ENDIF
#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, nelm_blk, gblock%nxblk*gblock%nyblk, &
         MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE mesh_filter


   SUBROUTINE mesh_filter_5x5 (gridf, dir, fname, vname)

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_NetCDFBlock
   USE MOD_DataType
   USE MOD_LandElm
   USE MOD_Mesh
   USE MOD_5x5DataReadin
   USE MOD_AggregationRequestData
   USE MOD_Block
   IMPLICIT NONE

   type(grid_type),  intent(in) :: gridf
   character(len=*), intent(in) :: dir
   character(len=*), intent(in) :: fname
   character(len=*), intent(in) :: vname

   ! local variables:
   ! ---------------------------------------------------------------
   type (block_data_int32_2d) :: datafilter
   integer, allocatable :: ifilter(:), xtemp(:), ytemp(:)
   logical, allocatable :: filter(:)
   integer :: ielm, jelm, npxl, nelm_glb

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,'(/, A)') 'Filtering pixels ...'
      ENDIF

      IF (p_is_io) THEN
         CALL allocate_block_data (gridf, datafilter)
         CALL read_5x5_data (trim(dir), trim(fname), gridf, trim(vname), datafilter)

#ifdef USEMPI
         CALL aggregation_data_daemon (gridf, data_i4_2d_in1 = datafilter)
#endif
      ENDIF

      IF (p_is_worker) THEN

         jelm = 0
         DO ielm = 1, numelm
            CALL aggregation_request_data (landelm, ielm, gridf, zip = .false., &
               data_i4_2d_in1 = datafilter, data_i4_2d_out1 = ifilter, &
               filledvalue_i4 = -1)

            allocate (filter (mesh(ielm)%npxl))
            filter = ifilter > 0

            IF (any(filter)) THEN
               jelm = jelm + 1
               IF (.not. all(filter)) THEN

                  npxl = count(filter)

                  allocate (xtemp(npxl))
                  allocate (ytemp(npxl))
                  xtemp = pack(mesh(ielm)%ilon, filter)
                  ytemp = pack(mesh(ielm)%ilat, filter)

                  deallocate(mesh(ielm)%ilon)
                  deallocate(mesh(ielm)%ilat)

                  mesh(ielm)%npxl = npxl

                  allocate(mesh(ielm)%ilon(npxl))
                  allocate(mesh(ielm)%ilat(npxl))
                  mesh(ielm)%ilon = xtemp
                  mesh(ielm)%ilat = ytemp

                  deallocate (xtemp)
                  deallocate (ytemp)
               ENDIF

               IF (jelm /= ielm) THEN
                  CALL copy_elm (mesh(ielm), mesh(jelm))
               ENDIF

            ENDIF

            deallocate (filter)
         ENDDO

         numelm = jelm

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

      IF (p_is_worker) THEN
         IF (allocated(landelm%eindex))   deallocate (landelm%eindex)
         IF (allocated(landelm%ipxstt))   deallocate (landelm%ipxstt)
         IF (allocated(landelm%ipxend))   deallocate (landelm%ipxend)
         IF (allocated(landelm%settyp))   deallocate (landelm%settyp)
         IF (allocated(landelm%ielm  ))   deallocate (landelm%ielm  )
      ENDIF

      CALL landelm_build ()

#ifdef USEMPI
      IF (p_is_worker) THEN
         CALL mpi_reduce (numelm, nelm_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', nelm_glb, ' elements after mesh filtering.'
         ENDIF
      ENDIF
#else
      write(*,'(A,I12,A)') 'Total: ', numelm, ' elements after mesh filtering.'
#endif

      ! Update nelm_blk
      nelm_blk(:,:) = 0
      IF (p_is_worker) THEN
         DO ielm = 1, numelm
            nelm_blk(mesh(ielm)%xblk,mesh(ielm)%yblk) = &
               nelm_blk(mesh(ielm)%xblk,mesh(ielm)%yblk) + 1
         ENDDO
      ENDIF
#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, nelm_blk, gblock%nxblk*gblock%nyblk, &
         MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE mesh_filter_5x5


END MODULE MOD_MeshFilter
