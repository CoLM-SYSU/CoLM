#include <define.h>

#ifdef GridRiverLakeFlow
MODULE MOD_Grid_RiverLakeVectorIO
!-----------------------------------------------------------------------
! DESCRIPTION:
!
!    Gather and scatter of (nrow, ncol) matrices whose columns are distributed
!    over workers by a global column id (bifurcation pathways, restart identity
!    tables of unit catchments and reservoirs).
!    Kept here, next to the river-lake routing, so MOD_Vector_ReadWrite stays
!    as upstream; it complements vector_gather_to_master/vector_read_and_scatter.
!-----------------------------------------------------------------------

   PUBLIC :: vector_gather_matrix_to_master
   PUBLIC :: vector_read_matrix_and_scatter

CONTAINS

   SUBROUTINE vector_gather_matrix_to_master ( &
         matrix, nrow, ncol_local, ncol_global, global_id, wdata)

   USE MOD_Precision
   USE MOD_SPMD_Task
   IMPLICIT NONE

   real(r8), intent(in) :: matrix (:,:)
   integer,  intent(in) :: nrow
   integer,  intent(in) :: ncol_local
   integer,  intent(in) :: ncol_global
   integer,  intent(in) :: global_id (:)

   real(r8), allocatable, intent(inout) :: wdata (:,:)

   integer :: iwork, isrc, ndata, ipth, mesg(2)
   integer, allocatable :: icache(:)
   real(r8), allocatable :: rcache(:)

      IF (nrow <= 0 .or. ncol_global <= 0) RETURN

      IF (size(global_id) /= ncol_local) THEN
         CALL CoLM_stop ('vector_gather_matrix_to_master: global_id size mismatch')
      ENDIF

      IF (p_is_master) THEN
         allocate (wdata (nrow, ncol_global))
         wdata(:,:) = 0._r8
      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN
         mesg = (/p_iam_glb, ncol_local/)
         CALL mpi_send (mesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         IF (ncol_local > 0) THEN
            CALL mpi_send (global_id, ncol_local, MPI_INTEGER, &
               p_address_master, mpi_tag_mesg + 1, p_comm_glb, p_err)
            CALL mpi_send (matrix, nrow*ncol_local, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_err)
         ENDIF
      ENDIF

      IF (p_is_master) THEN
         DO iwork = 0, p_np_worker-1
            CALL mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               allocate (icache (ndata))
               allocate (rcache (nrow*ndata))

               CALL mpi_recv (icache, ndata, MPI_INTEGER, isrc, &
                  mpi_tag_mesg + 1, p_comm_glb, p_stat, p_err)
               CALL mpi_recv (rcache, nrow*ndata, MPI_REAL8, isrc, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)

               DO ipth = 1, ndata
                  IF (icache(ipth) < 1 .or. icache(ipth) > ncol_global) THEN
                     CALL CoLM_stop ('vector_gather_matrix_to_master: global_id out of range')
                  ENDIF
                  wdata(:, icache(ipth)) = rcache((ipth-1)*nrow+1:ipth*nrow)
               ENDDO

               deallocate (icache)
               deallocate (rcache)
            ENDIF
         ENDDO
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      DO ipth = 1, ncol_local
         IF (global_id(ipth) < 1 .or. global_id(ipth) > ncol_global) THEN
            CALL CoLM_stop ('vector_gather_matrix_to_master: global_id out of range')
         ENDIF
         wdata(:, global_id(ipth)) = matrix(:, ipth)
      ENDDO
#endif

   END SUBROUTINE vector_gather_matrix_to_master

   ! -------
   SUBROUTINE vector_read_matrix_and_scatter ( &
         filein, matrix, nrow, ncol_local, varname, global_id, ncol_global)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   IMPLICIT NONE

   character(len=*),       intent(in)    :: filein
   real(r8),  allocatable, intent(inout) :: matrix (:,:)
   integer,                intent(in)    :: nrow
   integer,                intent(in)    :: ncol_local
   character(len=*),       intent(in)    :: varname
   integer,                intent(in)    :: global_id (:)
   integer,                intent(in)    :: ncol_global

   integer :: iwork, isrc, ndata, ipth, mesg(2)
   integer, allocatable :: icache(:)
   real(r8), allocatable :: rdata(:,:), rcache(:)

      IF (nrow <= 0 .or. ncol_global <= 0) RETURN

      IF (size(global_id) /= ncol_local) THEN
         CALL CoLM_stop ('vector_read_matrix_and_scatter: global_id size mismatch')
      ENDIF

      IF (p_is_worker) THEN
         IF (allocated(matrix)) THEN
            IF (size(matrix,1) /= nrow .or. size(matrix,2) /= ncol_local) deallocate(matrix)
         ENDIF
         IF (.not. allocated(matrix)) allocate (matrix(nrow, ncol_local))
      ENDIF

      IF (p_is_master) THEN
         CALL ncio_read_serial (filein, varname, rdata)

         IF (size(rdata, 1) /= nrow .or. size(rdata, 2) /= ncol_global) THEN
            CALL CoLM_stop ('vector_read_matrix_and_scatter: restart matrix shape mismatch')
         ENDIF
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_worker) THEN
         mesg = (/p_iam_glb, ncol_local/)
         CALL mpi_send (mesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)

         IF (ncol_local > 0) THEN
            CALL mpi_send (global_id, ncol_local, MPI_INTEGER, &
               p_address_master, mpi_tag_mesg + 1, p_comm_glb, p_err)

            CALL mpi_recv (matrix, nrow*ncol_local, MPI_REAL8, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF
      ENDIF

      IF (p_is_master) THEN
         DO iwork = 0, p_np_worker-1
            CALL mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               allocate (icache (ndata))
               allocate (rcache (nrow*ndata))

               CALL mpi_recv (icache, ndata, MPI_INTEGER, isrc, &
                  mpi_tag_mesg + 1, p_comm_glb, p_stat, p_err)

               DO ipth = 1, ndata
                  IF (icache(ipth) < 1 .or. icache(ipth) > ncol_global) THEN
                     CALL CoLM_stop ('vector_read_matrix_and_scatter: global_id out of range')
                  ENDIF
                  rcache((ipth-1)*nrow+1:ipth*nrow) = rdata(:, icache(ipth))
               ENDDO

               CALL mpi_send (rcache, nrow*ndata, MPI_REAL8, isrc, &
                  mpi_tag_data, p_comm_glb, p_err)

               deallocate (icache)
               deallocate (rcache)
            ENDIF
         ENDDO
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      DO ipth = 1, ncol_local
         IF (global_id(ipth) < 1 .or. global_id(ipth) > ncol_global) THEN
            CALL CoLM_stop ('vector_read_matrix_and_scatter: global_id out of range')
         ENDIF
         matrix(:, ipth) = rdata(:, global_id(ipth))
      ENDDO
#endif

      IF (p_is_master) deallocate (rdata)

   END SUBROUTINE vector_read_matrix_and_scatter

END MODULE MOD_Grid_RiverLakeVectorIO
#endif
