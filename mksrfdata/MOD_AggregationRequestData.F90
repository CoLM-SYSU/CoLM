#include <define.h>

MODULE MOD_AggregationRequestData

!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!    Aggregation Utilities.
!
!    On IO processes, a data daemon is running to provide data
!       at fine resolutions for worker processes.
!    On worker processes, request is sent to IO processes and
!       data is returned from IO processes.
!
!  Created by Shupeng Zhang, May 2023
!-----------------------------------------------------------------------

   IMPLICIT NONE

   PUBLIC :: aggregation_request_data
   PUBLIC :: aggregation_request_data_multigrid

#ifdef USEMPI
   PUBLIC :: aggregation_data_daemon
   PUBLIC :: aggregation_data_daemon_multigrid
   PUBLIC :: aggregation_worker_done
   PUBLIC :: aggregation_worker_done_multigrid
#endif

! ---- subroutines ----
CONTAINS

#ifdef USEMPI
   SUBROUTINE aggregation_data_daemon (grid_in, &
         data_r8_2d_in1, data_r8_2d_in2, data_r8_2d_in3, data_r8_2d_in4, &
         data_r8_2d_in5, data_r8_2d_in6,                                 &
         data_r8_3d_in1, n1_r8_3d_in1  , data_r8_3d_in2, n1_r8_3d_in2,   &
         data_i4_2d_in1, data_i4_2d_in2)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_DataType

   IMPLICIT NONE

   type (grid_type), intent(in) :: grid_in

   ! 2D REAL data
   type (block_data_real8_2d), intent(in), optional :: data_r8_2d_in1
   type (block_data_real8_2d), intent(in), optional :: data_r8_2d_in2
   type (block_data_real8_2d), intent(in), optional :: data_r8_2d_in3
   type (block_data_real8_2d), intent(in), optional :: data_r8_2d_in4
   type (block_data_real8_2d), intent(in), optional :: data_r8_2d_in5
   type (block_data_real8_2d), intent(in), optional :: data_r8_2d_in6

   ! 3D REAL data
   integer, intent(in), optional :: n1_r8_3d_in1
   type (block_data_real8_3d), intent(in), optional :: data_r8_3d_in1

   integer, intent(in), optional :: n1_r8_3d_in2
   type (block_data_real8_3d), intent(in), optional :: data_r8_3d_in2

   ! 2D INTEGER data
   type (block_data_int32_2d), intent(in), optional :: data_i4_2d_in1
   type (block_data_int32_2d), intent(in), optional :: data_i4_2d_in2

   ! Local Variables
   integer :: nreq, ireq, rmesg(2), isrc, idest
   integer :: xblk, yblk, xloc, yloc
   integer,  allocatable :: ylist(:), xlist(:)

   real(r8), allocatable :: sbuf_r8_1d(:), sbuf_r8_2d(:,:)
   integer , allocatable :: sbuf_i4_1d(:)

   logical,  allocatable :: worker_done (:)

      IF (p_is_io) THEN

         allocate (worker_done (0:p_np_worker-1))

         worker_done(:) = .false.
         DO WHILE (any(.not. worker_done))

            CALL mpi_recv (rmesg, 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc = rmesg(1)
            nreq = rmesg(2)

            IF (nreq > 0) THEN

               allocate (xlist (nreq))
               allocate (ylist (nreq))

               CALL mpi_recv (xlist, nreq, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL mpi_recv (ylist, nreq, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               idest = isrc

               allocate (sbuf_r8_1d (nreq))

               IF (present(data_r8_2d_in1)) THEN
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuf_r8_1d(ireq) = data_r8_2d_in1%blk(xblk,yblk)%val(xloc,yloc)
                  ENDDO

                  CALL mpi_send (sbuf_r8_1d, nreq, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
               ENDIF

               IF (present(data_r8_2d_in2)) THEN
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuf_r8_1d(ireq) = data_r8_2d_in2%blk(xblk,yblk)%val(xloc,yloc)
                  ENDDO

                  CALL mpi_send (sbuf_r8_1d, nreq, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
               ENDIF

               IF (present(data_r8_2d_in3)) THEN
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuf_r8_1d(ireq) = data_r8_2d_in3%blk(xblk,yblk)%val(xloc,yloc)
                  ENDDO

                  CALL mpi_send (sbuf_r8_1d, nreq, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
               ENDIF

               IF (present(data_r8_2d_in4)) THEN
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuf_r8_1d(ireq) = data_r8_2d_in4%blk(xblk,yblk)%val(xloc,yloc)
                  ENDDO

                  CALL mpi_send (sbuf_r8_1d, nreq, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
               ENDIF

               IF (present(data_r8_2d_in5)) THEN
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuf_r8_1d(ireq) = data_r8_2d_in5%blk(xblk,yblk)%val(xloc,yloc)
                  ENDDO

                  CALL mpi_send (sbuf_r8_1d, nreq, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
               ENDIF

               IF (present(data_r8_2d_in6)) THEN
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuf_r8_1d(ireq) = data_r8_2d_in6%blk(xblk,yblk)%val(xloc,yloc)
                  ENDDO

                  CALL mpi_send (sbuf_r8_1d, nreq, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
               ENDIF


               deallocate (sbuf_r8_1d)

               IF (present(data_r8_3d_in1) .and. present(n1_r8_3d_in1)) THEN

                  allocate (sbuf_r8_2d (n1_r8_3d_in1,nreq))
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuf_r8_2d(:,ireq) = data_r8_3d_in1%blk(xblk,yblk)%val(:,xloc,yloc)
                  ENDDO

                  CALL mpi_send (sbuf_r8_2d, n1_r8_3d_in1*nreq, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)

                  deallocate (sbuf_r8_2d)
               ENDIF

               IF (present(data_r8_3d_in2) .and. present(n1_r8_3d_in2)) THEN

                  allocate (sbuf_r8_2d (n1_r8_3d_in2,nreq))
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuf_r8_2d(:,ireq) = data_r8_3d_in2%blk(xblk,yblk)%val(:,xloc,yloc)
                  ENDDO

                  CALL mpi_send (sbuf_r8_2d, n1_r8_3d_in2*nreq, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)

                  deallocate (sbuf_r8_2d)
               ENDIF

               allocate (sbuf_i4_1d (nreq))

               IF (present(data_i4_2d_in1)) THEN
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuf_i4_1d(ireq) = data_i4_2d_in1%blk(xblk,yblk)%val(xloc,yloc)
                  ENDDO

                  CALL mpi_send (sbuf_i4_1d, nreq, MPI_INTEGER, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
               ENDIF

               IF (present(data_i4_2d_in2)) THEN
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuf_i4_1d(ireq) = data_i4_2d_in2%blk(xblk,yblk)%val(xloc,yloc)
                  ENDDO

                  CALL mpi_send (sbuf_i4_1d, nreq, MPI_INTEGER, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
               ENDIF

               deallocate (sbuf_i4_1d)

               deallocate (ylist)
               deallocate (xlist)

            ELSE
               worker_done(p_itis_worker(isrc)) = .true.
            ENDIF

         ENDDO

         deallocate (worker_done)

      ENDIF

   END SUBROUTINE aggregation_data_daemon


   SUBROUTINE aggregation_data_daemon_multigrid ( &
         grid_in1, data_i4_2d_in1, data_r8_2d_in1, data_r8_3d_in1, n1_r8_3d_in1, &
         grid_in2, data_i4_2d_in2, data_r8_2d_in2, data_r8_3d_in2, n1_r8_3d_in2, &
         grid_in3, data_i4_2d_in3, data_r8_2d_in3, data_r8_3d_in3, n1_r8_3d_in3, &
         grid_in4, data_i4_2d_in4, data_r8_2d_in4, data_r8_3d_in4, n1_r8_3d_in4, &
         grid_in5, data_i4_2d_in5, data_r8_2d_in5, data_r8_3d_in5, n1_r8_3d_in5)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_DataType

   IMPLICIT NONE

   type (grid_type), intent(in) :: grid_in1
   type (grid_type), intent(in) :: grid_in2
   type (grid_type), intent(in), optional :: grid_in3
   type (grid_type), intent(in), optional :: grid_in4
   type (grid_type), intent(in), optional :: grid_in5

   ! 2D INTEGER data
   type (block_data_int32_2d), intent(in), optional :: data_i4_2d_in1
   type (block_data_int32_2d), intent(in), optional :: data_i4_2d_in2
   type (block_data_int32_2d), intent(in), optional :: data_i4_2d_in3
   type (block_data_int32_2d), intent(in), optional :: data_i4_2d_in4
   type (block_data_int32_2d), intent(in), optional :: data_i4_2d_in5

   ! 2D REAL data
   type (block_data_real8_2d), intent(in), optional :: data_r8_2d_in1
   type (block_data_real8_2d), intent(in), optional :: data_r8_2d_in2
   type (block_data_real8_2d), intent(in), optional :: data_r8_2d_in3
   type (block_data_real8_2d), intent(in), optional :: data_r8_2d_in4
   type (block_data_real8_2d), intent(in), optional :: data_r8_2d_in5

   ! 3D REAL data
   integer, intent(in), optional :: n1_r8_3d_in1
   type (block_data_real8_3d), intent(in), optional :: data_r8_3d_in1

   integer, intent(in), optional :: n1_r8_3d_in2
   type (block_data_real8_3d), intent(in), optional :: data_r8_3d_in2

   integer, intent(in), optional :: n1_r8_3d_in3
   type (block_data_real8_3d), intent(in), optional :: data_r8_3d_in3

   integer, intent(in), optional :: n1_r8_3d_in4
   type (block_data_real8_3d), intent(in), optional :: data_r8_3d_in4

   integer, intent(in), optional :: n1_r8_3d_in5
   type (block_data_real8_3d), intent(in), optional :: data_r8_3d_in5

   ! Local Variables
   integer :: nreq, ireq, rmesg(3), isrc, idest, grdn
   integer :: xblk, yblk, xloc, yloc
   integer,  allocatable :: ylist(:), xlist(:)

   real(r8), allocatable :: sbuf_r8_1d(:), sbuf_r8_2d(:,:)
   integer , allocatable :: sbuf_i4_1d(:)

   logical,  allocatable :: worker_done (:)

      IF (p_is_io) THEN

         allocate (worker_done (0:p_np_worker-1))

         worker_done(:) = .false.
         DO WHILE (any(.not. worker_done))

            CALL mpi_recv (rmesg, 3, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc = rmesg(1)
            nreq = rmesg(2)
            grdn = rmesg(3)

            IF (nreq > 0) THEN

               allocate (xlist (nreq))
               allocate (ylist (nreq))

               CALL mpi_recv (xlist, nreq, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL mpi_recv (ylist, nreq, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               idest = isrc

               allocate (sbuf_i4_1d (nreq))

               IF (present(data_i4_2d_in1)) THEN
                  IF (grdn == 1) THEN
                     DO ireq = 1, nreq
                        xblk = grid_in1%xblk(xlist(ireq))
                        yblk = grid_in1%yblk(ylist(ireq))
                        xloc = grid_in1%xloc(xlist(ireq))
                        yloc = grid_in1%yloc(ylist(ireq))

                        sbuf_i4_1d(ireq) = data_i4_2d_in1%blk(xblk,yblk)%val(xloc,yloc)
                     ENDDO

                     CALL mpi_send (sbuf_i4_1d, nreq, MPI_INTEGER, &
                        idest, mpi_tag_data, p_comm_glb, p_err)
                  ENDIF
               ENDIF

               IF (present(data_i4_2d_in2)) THEN
                  IF (grdn == 2) THEN
                     DO ireq = 1, nreq
                        xblk = grid_in2%xblk(xlist(ireq))
                        yblk = grid_in2%yblk(ylist(ireq))
                        xloc = grid_in2%xloc(xlist(ireq))
                        yloc = grid_in2%yloc(ylist(ireq))

                        sbuf_i4_1d(ireq) = data_i4_2d_in2%blk(xblk,yblk)%val(xloc,yloc)
                     ENDDO

                     CALL mpi_send (sbuf_i4_1d, nreq, MPI_INTEGER, &
                        idest, mpi_tag_data, p_comm_glb, p_err)
                  ENDIF
               ENDIF

               IF (present(data_i4_2d_in3)) THEN
                  IF (grdn == 3) THEN
                     DO ireq = 1, nreq
                        xblk = grid_in3%xblk(xlist(ireq))
                        yblk = grid_in3%yblk(ylist(ireq))
                        xloc = grid_in3%xloc(xlist(ireq))
                        yloc = grid_in3%yloc(ylist(ireq))

                        sbuf_i4_1d(ireq) = data_i4_2d_in3%blk(xblk,yblk)%val(xloc,yloc)
                     ENDDO

                     CALL mpi_send (sbuf_i4_1d, nreq, MPI_INTEGER, &
                        idest, mpi_tag_data, p_comm_glb, p_err)
                  ENDIF
               ENDIF

               IF (present(data_i4_2d_in4)) THEN
                  IF (grdn == 4) THEN
                     DO ireq = 1, nreq
                        xblk = grid_in4%xblk(xlist(ireq))
                        yblk = grid_in4%yblk(ylist(ireq))
                        xloc = grid_in4%xloc(xlist(ireq))
                        yloc = grid_in4%yloc(ylist(ireq))

                        sbuf_i4_1d(ireq) = data_i4_2d_in4%blk(xblk,yblk)%val(xloc,yloc)
                     ENDDO

                     CALL mpi_send (sbuf_i4_1d, nreq, MPI_INTEGER, &
                        idest, mpi_tag_data, p_comm_glb, p_err)
                  ENDIF
               ENDIF

               IF (present(data_i4_2d_in5)) THEN
                  IF (grdn == 5) THEN
                     DO ireq = 1, nreq
                        xblk = grid_in5%xblk(xlist(ireq))
                        yblk = grid_in5%yblk(ylist(ireq))
                        xloc = grid_in5%xloc(xlist(ireq))
                        yloc = grid_in5%yloc(ylist(ireq))

                        sbuf_i4_1d(ireq) = data_i4_2d_in5%blk(xblk,yblk)%val(xloc,yloc)
                     ENDDO

                     CALL mpi_send (sbuf_i4_1d, nreq, MPI_INTEGER, &
                        idest, mpi_tag_data, p_comm_glb, p_err)
                  ENDIF
               ENDIF

               deallocate (sbuf_i4_1d)

               allocate (sbuf_r8_1d (nreq))

               IF (present(data_r8_2d_in1)) THEN
                  IF (grdn == 1) THEN
                     DO ireq = 1, nreq
                        xblk = grid_in1%xblk(xlist(ireq))
                        yblk = grid_in1%yblk(ylist(ireq))
                        xloc = grid_in1%xloc(xlist(ireq))
                        yloc = grid_in1%yloc(ylist(ireq))

                        sbuf_r8_1d(ireq) = data_r8_2d_in1%blk(xblk,yblk)%val(xloc,yloc)

                     ENDDO
                     CALL mpi_send (sbuf_r8_1d, nreq, MPI_REAL8, &
                        idest, mpi_tag_data, p_comm_glb, p_err)

                  ENDIF
               ENDIF

               IF ( present(data_r8_2d_in2) ) THEN
                  IF (grdn == 2) THEN
                     DO ireq = 1, nreq
                        xblk = grid_in2%xblk(xlist(ireq))
                        yblk = grid_in2%yblk(ylist(ireq))
                        xloc = grid_in2%xloc(xlist(ireq))
                        yloc = grid_in2%yloc(ylist(ireq))

                        sbuf_r8_1d(ireq) = data_r8_2d_in2%blk(xblk,yblk)%val(xloc,yloc)

                     ENDDO

                     CALL mpi_send (sbuf_r8_1d, nreq, MPI_REAL8, &
                        idest, mpi_tag_data, p_comm_glb, p_err)

                  ENDIF
               ENDIF

               IF ( present(data_r8_2d_in3) ) THEN
                  IF (grdn == 3) THEN
                     DO ireq = 1, nreq
                        xblk = grid_in3%xblk(xlist(ireq))
                        yblk = grid_in3%yblk(ylist(ireq))
                        xloc = grid_in3%xloc(xlist(ireq))
                        yloc = grid_in3%yloc(ylist(ireq))

                        sbuf_r8_1d(ireq) = data_r8_2d_in3%blk(xblk,yblk)%val(xloc,yloc)

                     ENDDO

                     CALL mpi_send (sbuf_r8_1d, nreq, MPI_REAL8, &
                        idest, mpi_tag_data, p_comm_glb, p_err)

                  ENDIF
               ENDIF

               IF ( present(data_r8_2d_in4) ) THEN
                  IF (grdn == 4) THEN
                     DO ireq = 1, nreq
                        xblk = grid_in4%xblk(xlist(ireq))
                        yblk = grid_in4%yblk(ylist(ireq))
                        xloc = grid_in4%xloc(xlist(ireq))
                        yloc = grid_in4%yloc(ylist(ireq))

                        sbuf_r8_1d(ireq) = data_r8_2d_in4%blk(xblk,yblk)%val(xloc,yloc)

                     ENDDO

                     CALL mpi_send (sbuf_r8_1d, nreq, MPI_REAL8, &
                        idest, mpi_tag_data, p_comm_glb, p_err)

                  ENDIF
               ENDIF

               IF ( present(data_r8_2d_in5) ) THEN
                  IF (grdn == 5) THEN
                     DO ireq = 1, nreq
                        xblk = grid_in5%xblk(xlist(ireq))
                        yblk = grid_in5%yblk(ylist(ireq))
                        xloc = grid_in5%xloc(xlist(ireq))
                        yloc = grid_in5%yloc(ylist(ireq))

                        sbuf_r8_1d(ireq) = data_r8_2d_in5%blk(xblk,yblk)%val(xloc,yloc)

                     ENDDO

                     CALL mpi_send (sbuf_r8_1d, nreq, MPI_REAL8, &
                        idest, mpi_tag_data, p_comm_glb, p_err)

                  ENDIF
               ENDIF

               deallocate (sbuf_r8_1d)

               IF ( present(data_r8_3d_in1) .and. present(n1_r8_3d_in1) ) THEN

                  allocate (sbuf_r8_2d (n1_r8_3d_in1,nreq))
                  IF (grdn == 1) THEN
                     DO ireq = 1, nreq
                        xblk = grid_in1%xblk(xlist(ireq))
                        yblk = grid_in1%yblk(ylist(ireq))
                        xloc = grid_in1%xloc(xlist(ireq))
                        yloc = grid_in1%yloc(ylist(ireq))

                        sbuf_r8_2d(:,ireq) = data_r8_3d_in1%blk(xblk,yblk)%val(:,xloc,yloc)

                     ENDDO

                     CALL mpi_send (sbuf_r8_2d, n1_r8_3d_in1*nreq, MPI_REAL8, &
                        idest, mpi_tag_data, p_comm_glb, p_err)

                  ENDIF

                 deallocate (sbuf_r8_2d)
               ENDIF

               IF ( present(data_r8_3d_in2) .and. present(n1_r8_3d_in2) ) THEN

                  allocate (sbuf_r8_2d (n1_r8_3d_in2,nreq))
                  IF (grdn == 2) THEN
                     DO ireq = 1, nreq
                        xblk = grid_in2%xblk(xlist(ireq))
                        yblk = grid_in2%yblk(ylist(ireq))
                        xloc = grid_in2%xloc(xlist(ireq))
                        yloc = grid_in2%yloc(ylist(ireq))

                        sbuf_r8_2d(:,ireq) = data_r8_3d_in2%blk(xblk,yblk)%val(:,xloc,yloc)

                     ENDDO

                     CALL mpi_send (sbuf_r8_2d, n1_r8_3d_in2*nreq, MPI_REAL8, &
                        idest, mpi_tag_data, p_comm_glb, p_err)

                  ENDIF

                 deallocate (sbuf_r8_2d)
               ENDIF

               IF ( present(data_r8_3d_in3) .and. present(n1_r8_3d_in3) ) THEN

                  allocate (sbuf_r8_2d (n1_r8_3d_in3,nreq))
                  IF (grdn == 3) THEN
                     DO ireq = 1, nreq
                        xblk = grid_in3%xblk(xlist(ireq))
                        yblk = grid_in3%yblk(ylist(ireq))
                        xloc = grid_in3%xloc(xlist(ireq))
                        yloc = grid_in3%yloc(ylist(ireq))

                        sbuf_r8_2d(:,ireq) = data_r8_3d_in3%blk(xblk,yblk)%val(:,xloc,yloc)

                     ENDDO

                     CALL mpi_send (sbuf_r8_2d, n1_r8_3d_in3*nreq, MPI_REAL8, &
                        idest, mpi_tag_data, p_comm_glb, p_err)

                  ENDIF

                 deallocate (sbuf_r8_2d)
               ENDIF

               IF ( present(data_r8_3d_in4) .and. present(n1_r8_3d_in4) ) THEN

                  allocate (sbuf_r8_2d (n1_r8_3d_in4,nreq))
                  IF (grdn == 4) THEN
                     DO ireq = 1, nreq
                        xblk = grid_in4%xblk(xlist(ireq))
                        yblk = grid_in4%yblk(ylist(ireq))
                        xloc = grid_in4%xloc(xlist(ireq))
                        yloc = grid_in4%yloc(ylist(ireq))

                        sbuf_r8_2d(:,ireq) = data_r8_3d_in4%blk(xblk,yblk)%val(:,xloc,yloc)

                     ENDDO

                     CALL mpi_send (sbuf_r8_2d, n1_r8_3d_in4*nreq, MPI_REAL8, &
                        idest, mpi_tag_data, p_comm_glb, p_err)

                  ENDIF

                 deallocate (sbuf_r8_2d)
               ENDIF

               IF ( present(data_r8_3d_in5) .and. present(n1_r8_3d_in5) ) THEN

                  allocate (sbuf_r8_2d (n1_r8_3d_in5,nreq))
                  IF (grdn == 5) THEN
                     DO ireq = 1, nreq
                        xblk = grid_in5%xblk(xlist(ireq))
                        yblk = grid_in5%yblk(ylist(ireq))
                        xloc = grid_in5%xloc(xlist(ireq))
                        yloc = grid_in5%yloc(ylist(ireq))

                        sbuf_r8_2d(:,ireq) = data_r8_3d_in5%blk(xblk,yblk)%val(:,xloc,yloc)

                     ENDDO

                     CALL mpi_send (sbuf_r8_2d, n1_r8_3d_in5*nreq, MPI_REAL8, &
                        idest, mpi_tag_data, p_comm_glb, p_err)

                  ENDIF

                 deallocate (sbuf_r8_2d)
               ENDIF

               deallocate (ylist)
               deallocate (xlist)

            ELSE
               worker_done(p_itis_worker(isrc)) = .true.
            ENDIF

         ENDDO

         deallocate (worker_done)

      ENDIF

   END SUBROUTINE aggregation_data_daemon_multigrid
#endif

   !----------------------------------------------------
   SUBROUTINE aggregation_request_data (     &
         pixelset, iset, grid_in, zip, area, &
         data_r8_2d_in1, data_r8_2d_out1, &
         data_r8_2d_in2, data_r8_2d_out2, &
         data_r8_2d_in3, data_r8_2d_out3, &
         data_r8_2d_in4, data_r8_2d_out4, &
         data_r8_2d_in5, data_r8_2d_out5, &
         data_r8_2d_in6, data_r8_2d_out6, &
         data_r8_3d_in1, data_r8_3d_out1, n1_r8_3d_in1, lb1_r8_3d_in1, &
         data_r8_3d_in2, data_r8_3d_out2, n1_r8_3d_in2, lb1_r8_3d_in2, &
         data_i4_2d_in1, data_i4_2d_out1, &
         data_i4_2d_in2, data_i4_2d_out2, &
         filledvalue_i4)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Block
   USE MOD_Pixel
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_Mesh
   USE MOD_Pixelset
   USE MOD_Utils

   IMPLICIT NONE

   type (pixelset_type), intent(in) :: pixelset
   integer, intent(in)  :: iset

   type (grid_type), intent(in) :: grid_in
   logical, intent(in) :: zip

   real(r8), allocatable, intent(out), optional :: area(:)

   type (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in1
   real(r8), allocatable,      intent(out), optional :: data_r8_2d_out1 (:)

   type (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in2
   real(r8), allocatable,      intent(out), optional :: data_r8_2d_out2 (:)

   type (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in3
   real(r8), allocatable,      intent(out), optional :: data_r8_2d_out3 (:)

   type (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in4
   real(r8), allocatable,      intent(out), optional :: data_r8_2d_out4 (:)

   type (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in5
   real(r8), allocatable,      intent(out), optional :: data_r8_2d_out5 (:)

   type (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in6
   real(r8), allocatable,      intent(out), optional :: data_r8_2d_out6 (:)

   integer, intent(in), optional :: n1_r8_3d_in1, lb1_r8_3d_in1
   type (block_data_real8_3d), intent(in),  optional :: data_r8_3d_in1
   real(r8), allocatable,      intent(out), optional :: data_r8_3d_out1 (:,:)

   integer, intent(in), optional :: n1_r8_3d_in2, lb1_r8_3d_in2
   type (block_data_real8_3d), intent(in),  optional :: data_r8_3d_in2
   real(r8), allocatable,      intent(out), optional :: data_r8_3d_out2 (:,:)

   type (block_data_int32_2d), intent(in),  optional :: data_i4_2d_in1
   integer, allocatable,       intent(out), optional :: data_i4_2d_out1 (:)

   type (block_data_int32_2d), intent(in),  optional :: data_i4_2d_in2
   integer, allocatable,       intent(out), optional :: data_i4_2d_out2 (:)

   integer, intent(in), optional :: filledvalue_i4

   ! Local Variables
   integer :: totalreq, ireq, nreq, smesg(2), isrc, idest, iproc
   integer :: ilon, ilat, xblk, yblk, xloc, yloc, iloc, nx, ny, ix, iy, ig
   integer :: ie, ipxstt, ipxend, npxl, ipxl, lb1, xgrdthis, ygrdthis
   integer,  allocatable :: ylist(:), xlist(:), ipt(:), ibuf(:), rbuf_i4_1d(:)
   integer,  allocatable :: xsorted(:), ysorted(:), xy2d(:,:)
   real(r8), allocatable :: area2d(:,:), rbuf_r8_1d(:), rbuf_r8_2d(:,:)
   logical,  allocatable :: msk(:)


      ie     = pixelset%ielm  (iset)
      ipxstt = pixelset%ipxstt(iset)
      ipxend = pixelset%ipxend(iset)
      npxl   = ipxend - ipxstt + 1

      IF (zip) THEN

         allocate (xsorted(npxl))
         allocate (ysorted(npxl))

         nx = 0; ny = 0
         DO ipxl = ipxstt, ipxend
            xgrdthis = grid_in%xgrd(mesh(ie)%ilon(ipxl))
            ygrdthis = grid_in%ygrd(mesh(ie)%ilat(ipxl))
            CALL insert_into_sorted_list1 (xgrdthis, nx, xsorted, iloc)
            CALL insert_into_sorted_list1 (ygrdthis, ny, ysorted, iloc)
         ENDDO

         allocate (xy2d (nx,ny));     xy2d(:,:) = 0

         IF (present(area)) THEN
            allocate(area2d(nx,ny));  area2d(:,:) = 0.
         ENDIF

         DO ipxl = ipxstt, ipxend
            xgrdthis = grid_in%xgrd(mesh(ie)%ilon(ipxl))
            ygrdthis = grid_in%ygrd(mesh(ie)%ilat(ipxl))

            ix = find_in_sorted_list1(xgrdthis, nx, xsorted)
            iy = find_in_sorted_list1(ygrdthis, ny, ysorted)

            xy2d(ix,iy) = xy2d(ix,iy) + 1

            IF (present(area)) THEN
               area2d(ix,iy) = area2d(ix,iy) + areaquad (&
                  pixel%lat_s(mesh(ie)%ilat(ipxl)), pixel%lat_n(mesh(ie)%ilat(ipxl)), &
                  pixel%lon_w(mesh(ie)%ilon(ipxl)), pixel%lon_e(mesh(ie)%ilon(ipxl)) )
            ENDIF
         ENDDO

         totalreq = count(xy2d > 0)

         allocate (xlist (totalreq))
         allocate (ylist (totalreq))

         IF (present(area)) allocate(area(totalreq))

         ig = 0
         DO ix = 1, nx
            DO iy = 1, ny
               IF (xy2d(ix,iy) > 0) THEN
                  ig = ig + 1
                  xlist(ig) = xsorted(ix)
                  ylist(ig) = ysorted(iy)
                  IF (present(area)) area (ig) = area2d(ix,iy)
               ENDIF
            ENDDO
         ENDDO

         deallocate (xsorted, ysorted, xy2d)
         IF (present(area)) deallocate (area2d)

      ELSE

         allocate(xlist (npxl))
         allocate(ylist (npxl))

         IF (present(area)) allocate (area (npxl))

         totalreq = npxl
         DO ipxl = ipxstt, ipxend
            xlist(ipxl-ipxstt+1) = grid_in%xgrd(mesh(ie)%ilon(ipxl))
            ylist(ipxl-ipxstt+1) = grid_in%ygrd(mesh(ie)%ilat(ipxl))
            IF (present(area)) THEN
               area(ipxl-ipxstt+1) = areaquad (&
                  pixel%lat_s(mesh(ie)%ilat(ipxl)), pixel%lat_n(mesh(ie)%ilat(ipxl)), &
                  pixel%lon_w(mesh(ie)%ilon(ipxl)), pixel%lon_e(mesh(ie)%ilon(ipxl)) )
            ENDIF
         ENDDO

      ENDIF

      IF (present(data_r8_2d_in1) .and. present(data_r8_2d_out1))  allocate (data_r8_2d_out1 (totalreq))
      IF (present(data_r8_2d_in2) .and. present(data_r8_2d_out2))  allocate (data_r8_2d_out2 (totalreq))
      IF (present(data_r8_2d_in3) .and. present(data_r8_2d_out3))  allocate (data_r8_2d_out3 (totalreq))
      IF (present(data_r8_2d_in4) .and. present(data_r8_2d_out4))  allocate (data_r8_2d_out4 (totalreq))
      IF (present(data_r8_2d_in5) .and. present(data_r8_2d_out5))  allocate (data_r8_2d_out5 (totalreq))
      IF (present(data_r8_2d_in6) .and. present(data_r8_2d_out6))  allocate (data_r8_2d_out6 (totalreq))

      IF (present(data_r8_3d_in1) .and. present(data_r8_3d_out1) .and. present(n1_r8_3d_in1)) THEN
         IF (present(lb1_r8_3d_in1)) THEN
            lb1 = lb1_r8_3d_in1
         ELSE
            lb1 = 1
         ENDIF
         allocate (data_r8_3d_out1 (lb1:lb1-1+n1_r8_3d_in1,totalreq))
      ENDIF

      IF (present(data_r8_3d_in2) .and. present(data_r8_3d_out2) .and. present(n1_r8_3d_in2)) THEN
         IF (present(lb1_r8_3d_in2)) THEN
            lb1 = lb1_r8_3d_in2
         ELSE
            lb1 = 1
         ENDIF
         allocate (data_r8_3d_out2 (lb1:lb1-1+n1_r8_3d_in2,totalreq))
      ENDIF

      IF (present(data_i4_2d_in1) .and. present(data_i4_2d_out1)) THEN
         allocate (data_i4_2d_out1 (totalreq))
         IF (present(filledvalue_i4)) THEN
            data_i4_2d_out1 = filledvalue_i4
         ENDIF
      ENDIF

      IF (present(data_i4_2d_in2) .and. present(data_i4_2d_out2)) THEN
         allocate (data_i4_2d_out2 (totalreq))
         IF (present(filledvalue_i4)) THEN
            data_i4_2d_out2 = filledvalue_i4
         ENDIF
      ENDIF

#ifdef USEMPI

      allocate (ipt (totalreq))
      allocate (msk (totalreq))

      ipt(:) = -1

      DO ireq = 1, totalreq
         xblk = grid_in%xblk(xlist(ireq))
         yblk = grid_in%yblk(ylist(ireq))
         ipt(ireq) = gblock%pio(xblk,yblk)
      ENDDO

      DO iproc = 0, p_np_io-1
         msk = (ipt == p_address_io(iproc))
         nreq = count(msk)

         IF (nreq > 0) THEN

            smesg = (/p_iam_glb, nreq/)
            idest = p_address_io(iproc)
            CALL mpi_send (smesg, 2, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)

            allocate (ibuf (nreq))

            ibuf = pack(xlist(1:totalreq), msk)
            CALL mpi_send (ibuf, nreq, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

            ibuf = pack(ylist(1:totalreq), msk)
            CALL mpi_send (ibuf, nreq, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

            isrc = idest

            allocate (rbuf_r8_1d (nreq))

            IF (present(data_r8_2d_in1) .and. present(data_r8_2d_out1)) THEN
               CALL mpi_recv (rbuf_r8_1d, nreq, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_r8_1d, msk, data_r8_2d_out1)
            ENDIF

            IF (present(data_r8_2d_in2) .and. present(data_r8_2d_out2)) THEN
               CALL mpi_recv (rbuf_r8_1d, nreq, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_r8_1d, msk, data_r8_2d_out2)
            ENDIF

            IF (present(data_r8_2d_in3) .and. present(data_r8_2d_out3)) THEN
               CALL mpi_recv (rbuf_r8_1d, nreq, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_r8_1d, msk, data_r8_2d_out3)
            ENDIF

            IF (present(data_r8_2d_in4) .and. present(data_r8_2d_out4)) THEN
               CALL mpi_recv (rbuf_r8_1d, nreq, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_r8_1d, msk, data_r8_2d_out4)
            ENDIF

            IF (present(data_r8_2d_in5) .and. present(data_r8_2d_out5)) THEN
               CALL mpi_recv (rbuf_r8_1d, nreq, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_r8_1d, msk, data_r8_2d_out5)
            ENDIF

            IF (present(data_r8_2d_in6) .and. present(data_r8_2d_out6)) THEN
               CALL mpi_recv (rbuf_r8_1d, nreq, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_r8_1d, msk, data_r8_2d_out6)
            ENDIF

            deallocate (rbuf_r8_1d)

            IF (present(data_r8_3d_in1) .and. present(data_r8_3d_out1) .and. present(n1_r8_3d_in1)) THEN
               allocate (rbuf_r8_2d (n1_r8_3d_in1,nreq))
               CALL mpi_recv (rbuf_r8_2d, n1_r8_3d_in1*nreq, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_r8_2d, msk, data_r8_3d_out1)
               deallocate (rbuf_r8_2d)
            ENDIF

            IF (present(data_r8_3d_in2) .and. present(data_r8_3d_out2) .and. present(n1_r8_3d_in2)) THEN
               allocate (rbuf_r8_2d (n1_r8_3d_in2,nreq))
               CALL mpi_recv (rbuf_r8_2d, n1_r8_3d_in2*nreq, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_r8_2d, msk, data_r8_3d_out2)
               deallocate (rbuf_r8_2d)
            ENDIF

            allocate (rbuf_i4_1d (nreq))
            IF (present(data_i4_2d_in1) .and. present(data_i4_2d_out1)) THEN
               CALL mpi_recv (rbuf_i4_1d, nreq, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_i4_1d, msk, data_i4_2d_out1)
            ENDIF

            IF (present(data_i4_2d_in2) .and. present(data_i4_2d_out2)) THEN
               CALL mpi_recv (rbuf_i4_1d, nreq, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_i4_1d, msk, data_i4_2d_out2)
            ENDIF

            deallocate (rbuf_i4_1d)

            deallocate (ibuf)
         ENDIF
      ENDDO

      deallocate (xlist)
      deallocate (ylist)
      deallocate (ipt  )
      deallocate (msk  )

#else

      DO ireq = 1, totalreq

         xblk = grid_in%xblk(xlist(ireq))
         yblk = grid_in%yblk(ylist(ireq))
         xloc = grid_in%xloc(xlist(ireq))
         yloc = grid_in%yloc(ylist(ireq))

         IF (present(data_r8_2d_in1) .and. present(data_r8_2d_out1)) THEN
            data_r8_2d_out1(ireq) = data_r8_2d_in1%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_2d_in2) .and. present(data_r8_2d_out2)) THEN
            data_r8_2d_out2(ireq) = data_r8_2d_in2%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_2d_in3) .and. present(data_r8_2d_out3)) THEN
            data_r8_2d_out3(ireq) = data_r8_2d_in3%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_2d_in4) .and. present(data_r8_2d_out4)) THEN
            data_r8_2d_out4(ireq) = data_r8_2d_in4%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_2d_in5) .and. present(data_r8_2d_out5)) THEN
            data_r8_2d_out5(ireq) = data_r8_2d_in5%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_2d_in6) .and. present(data_r8_2d_out6)) THEN
            data_r8_2d_out6(ireq) = data_r8_2d_in6%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_3d_in1) .and. present(data_r8_3d_out1) .and. present(n1_r8_3d_in1)) THEN
            data_r8_3d_out1(:,ireq) = data_r8_3d_in1%blk(xblk,yblk)%val(:,xloc,yloc)
         ENDIF

         IF (present(data_r8_3d_in2) .and. present(data_r8_3d_out2) .and. present(n1_r8_3d_in2)) THEN
            data_r8_3d_out2(:,ireq) = data_r8_3d_in2%blk(xblk,yblk)%val(:,xloc,yloc)
         ENDIF

         IF (present(data_i4_2d_in1) .and. present(data_i4_2d_out1)) THEN
            data_i4_2d_out1(ireq) = data_i4_2d_in1%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_i4_2d_in2) .and. present(data_i4_2d_out2)) THEN
            data_i4_2d_out2(ireq) = data_i4_2d_in2%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

      ENDDO

#endif

   END SUBROUTINE aggregation_request_data

   SUBROUTINE aggregation_request_data_multigrid (     &
         pixelset, iset, grid_in1, area , &
         data_i4_2d_in1, data_i4_2d_out1, &
         data_r8_2d_in1, data_r8_2d_out1, &
         data_r8_3d_in1, data_r8_3d_out1, n1_r8_3d_in1, lb1_r8_3d_in1, &
         grid_in2, &
         data_i4_2d_in2, data_i4_2d_out2, &
         data_r8_2d_in2, data_r8_2d_out2, &
         data_r8_3d_in2, data_r8_3d_out2, n1_r8_3d_in2, lb1_r8_3d_in2, &
         grid_in3, &
         data_i4_2d_in3, data_i4_2d_out3, &
         data_r8_2d_in3, data_r8_2d_out3, &
         data_r8_3d_in3, data_r8_3d_out3, n1_r8_3d_in3, lb1_r8_3d_in3, &
         grid_in4, &
         data_i4_2d_in4, data_i4_2d_out4, &
         data_r8_2d_in4, data_r8_2d_out4, &
         data_r8_3d_in4, data_r8_3d_out4, n1_r8_3d_in4, lb1_r8_3d_in4, &
         grid_in5, &
         data_i4_2d_in5, data_i4_2d_out5, &
         data_r8_2d_in5, data_r8_2d_out5, &
         data_r8_3d_in5, data_r8_3d_out5, n1_r8_3d_in5, lb1_r8_3d_in5, &
         filledvalue_i4)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Block
   USE MOD_Pixel
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_Mesh
   USE MOD_Pixelset
   USE MOD_Utils

   IMPLICIT NONE

   type (pixelset_type), intent(in) :: pixelset
   integer, intent(in)  :: iset

   type (grid_type), intent(in) :: grid_in1
   type (grid_type), intent(in) :: grid_in2
   type (grid_type), intent(in), optional :: grid_in3
   type (grid_type), intent(in), optional :: grid_in4
   type (grid_type), intent(in), optional :: grid_in5

   real(r8), allocatable, intent(out), optional :: area(:)

   type (block_data_int32_2d), intent(in),  optional :: data_i4_2d_in1
   integer, allocatable,       intent(out), optional :: data_i4_2d_out1 (:)

   type (block_data_int32_2d), intent(in),  optional :: data_i4_2d_in2
   integer, allocatable,       intent(out), optional :: data_i4_2d_out2 (:)

   type (block_data_int32_2d), intent(in),  optional :: data_i4_2d_in3
   integer, allocatable,       intent(out), optional :: data_i4_2d_out3 (:)

   type (block_data_int32_2d), intent(in),  optional :: data_i4_2d_in4
   integer, allocatable,       intent(out), optional :: data_i4_2d_out4 (:)

   type (block_data_int32_2d), intent(in),  optional :: data_i4_2d_in5
   integer, allocatable,       intent(out), optional :: data_i4_2d_out5 (:)

   type (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in1
   real(r8), allocatable,      intent(out), optional :: data_r8_2d_out1 (:)

   type (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in2
   real(r8), allocatable,      intent(out), optional :: data_r8_2d_out2 (:)

   type (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in3
   real(r8), allocatable,      intent(out), optional :: data_r8_2d_out3 (:)

   type (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in4
   real(r8), allocatable,      intent(out), optional :: data_r8_2d_out4 (:)

   type (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in5
   real(r8), allocatable,      intent(out), optional :: data_r8_2d_out5 (:)

   integer, intent(in), optional :: n1_r8_3d_in1, lb1_r8_3d_in1
   type (block_data_real8_3d), intent(in),  optional :: data_r8_3d_in1
   real(r8), allocatable,      intent(out), optional :: data_r8_3d_out1 (:,:)

   integer, intent(in), optional :: n1_r8_3d_in2, lb1_r8_3d_in2
   type (block_data_real8_3d), intent(in),  optional :: data_r8_3d_in2
   real(r8), allocatable,      intent(out), optional :: data_r8_3d_out2 (:,:)

   integer, intent(in), optional :: n1_r8_3d_in3, lb1_r8_3d_in3
   type (block_data_real8_3d), intent(in),  optional :: data_r8_3d_in3
   real(r8), allocatable,      intent(out), optional :: data_r8_3d_out3 (:,:)

   integer, intent(in), optional :: n1_r8_3d_in4, lb1_r8_3d_in4
   type (block_data_real8_3d), intent(in),  optional :: data_r8_3d_in4
   real(r8), allocatable,      intent(out), optional :: data_r8_3d_out4 (:,:)

   integer, intent(in), optional :: n1_r8_3d_in5, lb1_r8_3d_in5
   type (block_data_real8_3d), intent(in),  optional :: data_r8_3d_in5
   real(r8), allocatable,      intent(out), optional :: data_r8_3d_out5 (:,:)

   integer, intent(in), optional :: filledvalue_i4

   ! Local Variables
   integer :: totalreq, ireq, nreq1, nreq2, nreq3, nreq4, nreq5, smesg(3), isrc, idest, iproc
   integer :: ilon, ilat, xblk, yblk, xloc, yloc, iloc, nx, ny, ix, iy, ig
   integer :: ie, ipxstt, ipxend, npxl, ipxl, lb1, xgrdthis, ygrdthis
   integer,  allocatable :: ylist1(:), xlist1(:), ipt1(:)
   integer,  allocatable :: ylist2(:), xlist2(:), ipt2(:)
   integer,  allocatable :: ylist3(:), xlist3(:), ipt3(:)
   integer,  allocatable :: ylist4(:), xlist4(:), ipt4(:)
   integer,  allocatable :: ylist5(:), xlist5(:), ipt5(:)
   integer,  allocatable :: ibuf(:), rbuf_i4_1d(:)
   integer,  allocatable :: xsorted(:), ysorted(:), xy2d(:,:)
   real(r8), allocatable :: area2d(:,:), rbuf_r8_1d(:), rbuf_r8_2d(:,:)
   logical,  allocatable :: msk1(:), msk2(:), msk3(:), msk4(:), msk5(:)


      ie     = pixelset%ielm  (iset)
      ipxstt = pixelset%ipxstt(iset)
      ipxend = pixelset%ipxend(iset)
      npxl   = ipxend - ipxstt + 1

      allocate(xlist1 (npxl))
      allocate(ylist1 (npxl))

      allocate(xlist2 (npxl))
      allocate(ylist2 (npxl))

      IF (present(grid_in3)) THEN
         allocate(xlist3 (npxl))
         allocate(ylist3 (npxl))
      ENDIF

      IF (present(grid_in4)) THEN
         allocate(xlist4 (npxl))
         allocate(ylist4 (npxl))
      ENDIF

      IF (present(grid_in5)) THEN
         allocate(xlist5 (npxl))
         allocate(ylist5 (npxl))
      ENDIF

      IF (present(area)) allocate (area (npxl))

      totalreq = npxl
      DO ipxl = ipxstt, ipxend
         xlist1(ipxl-ipxstt+1) = grid_in1%xgrd(mesh(ie)%ilon(ipxl))
         ylist1(ipxl-ipxstt+1) = grid_in1%ygrd(mesh(ie)%ilat(ipxl))

         xlist2(ipxl-ipxstt+1) = grid_in2%xgrd(mesh(ie)%ilon(ipxl))
         ylist2(ipxl-ipxstt+1) = grid_in2%ygrd(mesh(ie)%ilat(ipxl))

         IF (present(grid_in3)) THEN
            xlist3(ipxl-ipxstt+1) = grid_in3%xgrd(mesh(ie)%ilon(ipxl))
            ylist3(ipxl-ipxstt+1) = grid_in3%ygrd(mesh(ie)%ilat(ipxl))
         ENDIF

         IF (present(grid_in4)) THEN
            xlist4(ipxl-ipxstt+1) = grid_in4%xgrd(mesh(ie)%ilon(ipxl))
            ylist4(ipxl-ipxstt+1) = grid_in4%ygrd(mesh(ie)%ilat(ipxl))
         ENDIF

         IF (present(grid_in5)) THEN
            xlist5(ipxl-ipxstt+1) = grid_in5%xgrd(mesh(ie)%ilon(ipxl))
            ylist5(ipxl-ipxstt+1) = grid_in5%ygrd(mesh(ie)%ilat(ipxl))
         ENDIF

         IF (present(area)) THEN
            area(ipxl-ipxstt+1) = areaquad (&
               pixel%lat_s(mesh(ie)%ilat(ipxl)), pixel%lat_n(mesh(ie)%ilat(ipxl)), &
               pixel%lon_w(mesh(ie)%ilon(ipxl)), pixel%lon_e(mesh(ie)%ilon(ipxl)) )
         ENDIF
      ENDDO

      IF (present(data_i4_2d_in1) .and. present(data_i4_2d_out1)) THEN
         allocate (data_i4_2d_out1 (totalreq))
         IF (present(filledvalue_i4)) THEN
            data_i4_2d_out1 = filledvalue_i4
         ENDIF
      ENDIF

      IF (present(data_i4_2d_in2) .and. present(data_i4_2d_out2)) THEN
         allocate (data_i4_2d_out2 (totalreq))
         IF (present(filledvalue_i4)) THEN
            data_i4_2d_out2 = filledvalue_i4
         ENDIF
      ENDIF

      IF (present(data_i4_2d_in3) .and. present(data_i4_2d_out3)) THEN
         allocate (data_i4_2d_out3 (totalreq))
         IF (present(filledvalue_i4)) THEN
            data_i4_2d_out3 = filledvalue_i4
         ENDIF
      ENDIF

      IF (present(data_i4_2d_in4) .and. present(data_i4_2d_out4)) THEN
         allocate (data_i4_2d_out4 (totalreq))
         IF (present(filledvalue_i4)) THEN
            data_i4_2d_out4 = filledvalue_i4
         ENDIF
      ENDIF

      IF (present(data_i4_2d_in5) .and. present(data_i4_2d_out5)) THEN
         allocate (data_i4_2d_out5 (totalreq))
         IF (present(filledvalue_i4)) THEN
            data_i4_2d_out5 = filledvalue_i4
         ENDIF
      ENDIF

      IF (present(data_r8_2d_in1) .and. present(data_r8_2d_out1))  allocate (data_r8_2d_out1 (totalreq))
      IF (present(data_r8_2d_in2) .and. present(data_r8_2d_out2))  allocate (data_r8_2d_out2 (totalreq))
      IF (present(data_r8_2d_in3) .and. present(data_r8_2d_out3))  allocate (data_r8_2d_out3 (totalreq))
      IF (present(data_r8_2d_in4) .and. present(data_r8_2d_out4))  allocate (data_r8_2d_out4 (totalreq))
      IF (present(data_r8_2d_in5) .and. present(data_r8_2d_out5))  allocate (data_r8_2d_out5 (totalreq))

      IF (present(data_r8_3d_in1) .and. present(data_r8_3d_out1) .and. present(n1_r8_3d_in1)) THEN
         IF (present(lb1_r8_3d_in1)) THEN
            lb1 = lb1_r8_3d_in1
         ELSE
            lb1 = 1
         ENDIF
         allocate (data_r8_3d_out1 (lb1:lb1-1+n1_r8_3d_in1,totalreq))
      ENDIF

      IF (present(data_r8_3d_in2) .and. present(data_r8_3d_out2) .and. present(n1_r8_3d_in2)) THEN
         IF (present(lb1_r8_3d_in2)) THEN
            lb1 = lb1_r8_3d_in2
         ELSE
            lb1 = 1
         ENDIF
         allocate (data_r8_3d_out2 (lb1:lb1-1+n1_r8_3d_in2,totalreq))
      ENDIF

      IF (present(data_r8_3d_in3) .and. present(data_r8_3d_out3) .and. present(n1_r8_3d_in3)) THEN
         IF (present(lb1_r8_3d_in3)) THEN
            lb1 = lb1_r8_3d_in3
         ELSE
            lb1 = 1
         ENDIF
         allocate (data_r8_3d_out3 (lb1:lb1-1+n1_r8_3d_in3,totalreq))
      ENDIF

      IF (present(data_r8_3d_in4) .and. present(data_r8_3d_out4) .and. present(n1_r8_3d_in4)) THEN
         IF (present(lb1_r8_3d_in4)) THEN
            lb1 = lb1_r8_3d_in4
         ELSE
            lb1 = 1
         ENDIF
         allocate (data_r8_3d_out4 (lb1:lb1-1+n1_r8_3d_in4,totalreq))
      ENDIF

      IF (present(data_r8_3d_in5) .and. present(data_r8_3d_out5) .and. present(n1_r8_3d_in5)) THEN
         IF (present(lb1_r8_3d_in5)) THEN
            lb1 = lb1_r8_3d_in5
         ELSE
            lb1 = 1
         ENDIF
         allocate (data_r8_3d_out5 (lb1:lb1-1+n1_r8_3d_in5,totalreq))
      ENDIF

#ifdef USEMPI

      allocate (ipt1 (totalreq))
      allocate (msk1 (totalreq))

      allocate (ipt2 (totalreq))
      allocate (msk2 (totalreq))

      ipt1(:) = -1
      ipt2(:) = -1

      IF (present(grid_in3)) THEN
         allocate (ipt3 (totalreq))
         allocate (msk3 (totalreq))

         ipt3(:) = -1
      ENDIF

      IF (present(grid_in4)) THEN
         allocate (ipt4 (totalreq))
         allocate (msk4 (totalreq))

         ipt4(:) = -1
      ENDIF

      IF (present(grid_in5)) THEN
         allocate (ipt5 (totalreq))
         allocate (msk5 (totalreq))

         ipt5(:) = -1
      ENDIF

      DO ireq = 1, totalreq
         xblk = grid_in1%xblk(xlist1(ireq))
         yblk = grid_in1%yblk(ylist1(ireq))
         ipt1(ireq) = gblock%pio(xblk,yblk)

         xblk = grid_in2%xblk(xlist2(ireq))
         yblk = grid_in2%yblk(ylist2(ireq))
         ipt2(ireq) = gblock%pio(xblk,yblk)

         IF (present(grid_in3)) THEN
            xblk = grid_in3%xblk(xlist3(ireq))
            yblk = grid_in3%yblk(ylist3(ireq))
            ipt3(ireq) = gblock%pio(xblk,yblk)
         ENDIF

         IF (present(grid_in4)) THEN
            xblk = grid_in4%xblk(xlist4(ireq))
            yblk = grid_in4%yblk(ylist4(ireq))
            ipt4(ireq) = gblock%pio(xblk,yblk)
         ENDIF

         IF (present(grid_in5)) THEN
            xblk = grid_in5%xblk(xlist5(ireq))
            yblk = grid_in5%yblk(ylist5(ireq))
            ipt5(ireq) = gblock%pio(xblk,yblk)
         ENDIF
      ENDDO

      DO iproc = 0, p_np_io-1
         msk1 = (ipt1 == p_address_io(iproc))
         nreq1 = count(msk1)

         IF (nreq1 > 0) THEN

            smesg = (/p_iam_glb, nreq1, 1/)

            idest = p_address_io(iproc)
            CALL mpi_send (smesg, 3, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)

            allocate (ibuf (nreq1))

            ibuf = pack(xlist1(1:totalreq), msk1)
            CALL mpi_send (ibuf, nreq1, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

            ibuf = pack(ylist1(1:totalreq), msk1)
            CALL mpi_send (ibuf, nreq1, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

            isrc = idest

            allocate (rbuf_i4_1d (nreq1))

            IF (present(data_i4_2d_in1) .and. present(data_i4_2d_out1)) THEN
               CALL mpi_recv (rbuf_i4_1d, nreq1, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_i4_1d, msk1, data_i4_2d_out1)
            ENDIF

            deallocate (rbuf_i4_1d)

            allocate (rbuf_r8_1d (nreq1))

            IF (present(data_r8_2d_in1) .and. present(data_r8_2d_out1)) THEN
               CALL mpi_recv (rbuf_r8_1d, nreq1, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_r8_1d, msk1, data_r8_2d_out1)
            ENDIF

            deallocate (rbuf_r8_1d)

            IF (present(data_r8_3d_in1) .and. present(data_r8_3d_out1) .and. present(n1_r8_3d_in1)) THEN
               allocate (rbuf_r8_2d (n1_r8_3d_in1,nreq1))
               CALL mpi_recv (rbuf_r8_2d, n1_r8_3d_in1*nreq1, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_r8_2d, msk1, data_r8_3d_out1)
               deallocate (rbuf_r8_2d)
            ENDIF

            deallocate (ibuf)
         ENDIF

         msk2 = (ipt2 == p_address_io(iproc))
         nreq2 = count(msk2)

         IF (nreq2 > 0) THEN
            smesg = (/p_iam_glb, nreq2, 2/)

            idest = p_address_io(iproc)
            CALL mpi_send (smesg, 3, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)

            allocate (ibuf (nreq2))

            ibuf = pack(xlist2(1:totalreq), msk2)
            CALL mpi_send (ibuf, nreq2, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

            ibuf = pack(ylist2(1:totalreq), msk2)
            CALL mpi_send (ibuf, nreq2, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

            isrc = idest

            allocate (rbuf_i4_1d (nreq2))

            IF (present(data_i4_2d_in2) .and. present(data_i4_2d_out2)) THEN
               CALL mpi_recv (rbuf_i4_1d, nreq2, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_i4_1d, msk2, data_i4_2d_out2)
            ENDIF

            deallocate (rbuf_i4_1d)

            allocate (rbuf_r8_1d (nreq2))

            IF (present(data_r8_2d_in2) .and. present(data_r8_2d_out2)) THEN
               CALL mpi_recv (rbuf_r8_1d, nreq2, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_r8_1d, msk2, data_r8_2d_out2)
            ENDIF

            deallocate (rbuf_r8_1d)

            IF (present(data_r8_3d_in2) .and. present(data_r8_3d_out2)) THEN
               allocate (rbuf_r8_2d (n1_r8_3d_in2,nreq2))
               CALL mpi_recv (rbuf_r8_2d, n1_r8_3d_in2*nreq2, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_r8_2d, msk2, data_r8_3d_out2)
               deallocate (rbuf_r8_2d)
            ENDIF

            deallocate (ibuf)
         ENDIF

         IF (present(grid_in3)) THEN
            msk3 = (ipt3 == p_address_io(iproc))
            nreq3= count(msk3)

            IF (nreq3 > 0) THEN
               smesg = (/p_iam_glb, nreq3, 3/)

               idest = p_address_io(iproc)
               CALL mpi_send (smesg, 3, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)

               allocate (ibuf (nreq3))

               ibuf = pack(xlist3(1:totalreq), msk3)
               CALL mpi_send (ibuf, nreq3, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

               ibuf = pack(ylist3(1:totalreq), msk3)
               CALL mpi_send (ibuf, nreq3, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

               isrc = idest

               allocate (rbuf_i4_1d (nreq3))

               IF (present(data_i4_2d_in3) .and. present(data_i4_2d_out3)) THEN
                  CALL mpi_recv (rbuf_i4_1d, nreq3, MPI_INTEGER, &
                     isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
                  CALL unpack_inplace (rbuf_i4_1d, msk3, data_i4_2d_out3)
               ENDIF

               deallocate (rbuf_i4_1d)

               allocate (rbuf_r8_1d (nreq3))

               IF (present(data_r8_2d_in3) .and. present(data_r8_2d_out3)) THEN
                  CALL mpi_recv (rbuf_r8_1d, nreq3, MPI_REAL8, &
                     isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
                  CALL unpack_inplace (rbuf_r8_1d, msk3, data_r8_2d_out3)
               ENDIF

               deallocate (rbuf_r8_1d)

               IF (present(data_r8_3d_in3) .and. present(data_r8_3d_out3)) THEN
                  allocate (rbuf_r8_2d (n1_r8_3d_in3,nreq3))
                  CALL mpi_recv (rbuf_r8_2d, n1_r8_3d_in3*nreq3, MPI_REAL8, &
                     isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
                  CALL unpack_inplace (rbuf_r8_2d, msk3, data_r8_3d_out3)
                  deallocate (rbuf_r8_2d)
               ENDIF

               deallocate (ibuf)
            ENDIF
         ENDIF

         IF (present(grid_in4)) THEN
            msk4 = (ipt4 == p_address_io(iproc))
            nreq4= count(msk4)

            IF (nreq4 > 0) THEN
               smesg = (/p_iam_glb, nreq4, 4/)

               idest = p_address_io(iproc)
               CALL mpi_send (smesg, 3, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)

               allocate (ibuf (nreq4))

               ibuf = pack(xlist4(1:totalreq), msk4)
               CALL mpi_send (ibuf, nreq4, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

               ibuf = pack(ylist4(1:totalreq), msk4)
               CALL mpi_send (ibuf, nreq4, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

               isrc = idest

               allocate (rbuf_i4_1d (nreq4))

               IF (present(data_i4_2d_in4) .and. present(data_i4_2d_out4)) THEN
                  CALL mpi_recv (rbuf_i4_1d, nreq4, MPI_INTEGER, &
                     isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
                  CALL unpack_inplace (rbuf_i4_1d, msk4, data_i4_2d_out4)
               ENDIF

               deallocate (rbuf_i4_1d)

               allocate (rbuf_r8_1d (nreq4))

               IF (present(data_r8_2d_in4) .and. present(data_r8_2d_out4)) THEN
                  CALL mpi_recv (rbuf_r8_1d, nreq4, MPI_REAL8, &
                     isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
                  CALL unpack_inplace (rbuf_r8_1d, msk4, data_r8_2d_out4)
               ENDIF

               deallocate (rbuf_r8_1d)

               IF (present(data_r8_3d_in4) .and. present(data_r8_3d_out4)) THEN
                  allocate (rbuf_r8_2d (n1_r8_3d_in4,nreq4))
                  CALL mpi_recv (rbuf_r8_2d, n1_r8_3d_in4*nreq4, MPI_REAL8, &
                     isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
                  CALL unpack_inplace (rbuf_r8_2d, msk4, data_r8_3d_out4)
                  deallocate (rbuf_r8_2d)
               ENDIF

               deallocate (ibuf)
            ENDIF
         ENDIF

         IF (present(grid_in5)) THEN
            msk5 = (ipt5 == p_address_io(iproc))
            nreq5= count(msk5)

            IF (nreq5 > 0) THEN
               smesg = (/p_iam_glb, nreq5, 5/)

               idest = p_address_io(iproc)
               CALL mpi_send (smesg, 3, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)

               allocate (ibuf (nreq5))

               ibuf = pack(xlist5(1:totalreq), msk5)
               CALL mpi_send (ibuf, nreq5, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

               ibuf = pack(ylist5(1:totalreq), msk5)
               CALL mpi_send (ibuf, nreq5, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

               isrc = idest

               allocate (rbuf_i4_1d (nreq5))

               IF (present(data_i4_2d_in5) .and. present(data_i4_2d_out5)) THEN
                  CALL mpi_recv (rbuf_i4_1d, nreq5, MPI_INTEGER, &
                     isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
                  CALL unpack_inplace (rbuf_i4_1d, msk5, data_i4_2d_out5)
               ENDIF

               deallocate (rbuf_i4_1d)

               allocate (rbuf_r8_1d (nreq5))

               IF (present(data_r8_2d_in5) .and. present(data_r8_2d_out5)) THEN
                  CALL mpi_recv (rbuf_r8_1d, nreq5, MPI_REAL8, &
                     isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
                  CALL unpack_inplace (rbuf_r8_1d, msk5, data_r8_2d_out5)
               ENDIF

               deallocate (rbuf_r8_1d)

               IF (present(data_r8_3d_in5) .and. present(data_r8_3d_out5)) THEN
                  allocate (rbuf_r8_2d (n1_r8_3d_in5,nreq5))
                  CALL mpi_recv (rbuf_r8_2d, n1_r8_3d_in5*nreq5, MPI_REAL8, &
                     isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
                  CALL unpack_inplace (rbuf_r8_2d, msk5, data_r8_3d_out5)
                  deallocate (rbuf_r8_2d)
               ENDIF

               deallocate (ibuf)
            ENDIF
         ENDIF
      ENDDO

      deallocate (xlist1)
      deallocate (ylist1)
      deallocate (ipt1  )
      deallocate (msk1  )

      deallocate (xlist2)
      deallocate (ylist2)
      deallocate (ipt2  )
      deallocate (msk2  )

      IF (allocated(xlist3)) deallocate (xlist3)
      IF (allocated(ylist3)) deallocate (ylist3)
      IF (allocated(ipt3  )) deallocate (ipt3  )
      IF (allocated(msk3  )) deallocate (msk3  )

      IF (allocated(xlist4)) deallocate (xlist4)
      IF (allocated(ylist4)) deallocate (ylist4)
      IF (allocated(ipt4  )) deallocate (ipt4  )
      IF (allocated(msk4  )) deallocate (msk4  )

      IF (allocated(xlist5)) deallocate (xlist5)
      IF (allocated(ylist5)) deallocate (ylist5)
      IF (allocated(ipt5  )) deallocate (ipt5  )
      IF (allocated(msk5  )) deallocate (msk5  )
#else

      DO ireq = 1, totalreq

         IF (present(data_i4_2d_in1) .and. present(data_i4_2d_out1)) THEN
            xblk = grid_in1%xblk(xlist1(ireq))
            yblk = grid_in1%yblk(ylist1(ireq))
            xloc = grid_in1%xloc(xlist1(ireq))
            yloc = grid_in1%yloc(ylist1(ireq))

            data_i4_2d_out1(ireq) = data_i4_2d_in1%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_i4_2d_in2) .and. present(data_i4_2d_out2)) THEN
            xblk = grid_in2%xblk(xlist2(ireq))
            yblk = grid_in2%yblk(ylist2(ireq))
            xloc = grid_in2%xloc(xlist2(ireq))
            yloc = grid_in2%yloc(ylist2(ireq))

            data_i4_2d_out2(ireq) = data_i4_2d_in2%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_i4_2d_in3) .and. present(data_i4_2d_out3)) THEN
            xblk = grid_in3%xblk(xlist3(ireq))
            yblk = grid_in3%yblk(ylist3(ireq))
            xloc = grid_in3%xloc(xlist3(ireq))
            yloc = grid_in3%yloc(ylist3(ireq))

            data_i4_2d_out3(ireq) = data_i4_2d_in3%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_i4_2d_in4) .and. present(data_i4_2d_out4)) THEN
            xblk = grid_in4%xblk(xlist4(ireq))
            yblk = grid_in4%yblk(ylist4(ireq))
            xloc = grid_in4%xloc(xlist4(ireq))
            yloc = grid_in4%yloc(ylist4(ireq))

            data_i4_2d_out4(ireq) = data_i4_2d_in4%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_i4_2d_in5) .and. present(data_i4_2d_out5)) THEN
            xblk = grid_in5%xblk(xlist5(ireq))
            yblk = grid_in5%yblk(ylist5(ireq))
            xloc = grid_in5%xloc(xlist5(ireq))
            yloc = grid_in5%yloc(ylist5(ireq))

            data_i4_2d_out5(ireq) = data_i4_2d_in5%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_2d_in1) .and. present(data_r8_2d_out1)) THEN
            xblk = grid_in1%xblk(xlist1(ireq))
            yblk = grid_in1%yblk(ylist1(ireq))
            xloc = grid_in1%xloc(xlist1(ireq))
            yloc = grid_in1%yloc(ylist1(ireq))

            data_r8_2d_out1(ireq) = data_r8_2d_in1%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_2d_in2) .and. present(data_r8_2d_out2)) THEN
            xblk = grid_in2%xblk(xlist2(ireq))
            yblk = grid_in2%yblk(ylist2(ireq))
            xloc = grid_in2%xloc(xlist2(ireq))
            yloc = grid_in2%yloc(ylist2(ireq))

            data_r8_2d_out2(ireq) = data_r8_2d_in2%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_2d_in3) .and. present(data_r8_2d_out3)) THEN
            xblk = grid_in3%xblk(xlist3(ireq))
            yblk = grid_in3%yblk(ylist3(ireq))
            xloc = grid_in3%xloc(xlist3(ireq))
            yloc = grid_in3%yloc(ylist3(ireq))

            data_r8_2d_out3(ireq) = data_r8_2d_in3%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_2d_in4) .and. present(data_r8_2d_out4)) THEN
            xblk = grid_in4%xblk(xlist4(ireq))
            yblk = grid_in4%yblk(ylist4(ireq))
            xloc = grid_in4%xloc(xlist4(ireq))
            yloc = grid_in4%yloc(ylist4(ireq))

            data_r8_2d_out4(ireq) = data_r8_2d_in4%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_3d_in1) .and. present(data_r8_3d_out1) .and. present(n1_r8_3d_in1)) THEN
            xblk = grid_in1%xblk(xlist1(ireq))
            yblk = grid_in1%yblk(ylist1(ireq))
            xloc = grid_in1%xloc(xlist1(ireq))
            yloc = grid_in1%yloc(ylist1(ireq))

            data_r8_3d_out1(:,ireq) = data_r8_3d_in1%blk(xblk,yblk)%val(:,xloc,yloc)
         ENDIF

         IF (present(data_r8_3d_in2) .and. present(data_r8_3d_out2) .and. present(n1_r8_3d_in2)) THEN
            xblk = grid_in2%xblk(xlist2(ireq))
            yblk = grid_in2%yblk(ylist2(ireq))
            xloc = grid_in2%xloc(xlist2(ireq))
            yloc = grid_in2%yloc(ylist2(ireq))

            data_r8_3d_out2(:,ireq) = data_r8_3d_in2%blk(xblk,yblk)%val(:,xloc,yloc)
         ENDIF

         IF (present(data_r8_3d_in3) .and. present(data_r8_3d_out3) .and. present(n1_r8_3d_in3)) THEN
            xblk = grid_in3%xblk(xlist3(ireq))
            yblk = grid_in3%yblk(ylist3(ireq))
            xloc = grid_in3%xloc(xlist3(ireq))
            yloc = grid_in3%yloc(ylist3(ireq))

            data_r8_3d_out3(:,ireq) = data_r8_3d_in3%blk(xblk,yblk)%val(:,xloc,yloc)
         ENDIF

         IF (present(data_r8_3d_in4) .and. present(data_r8_3d_out4) .and. present(n1_r8_3d_in4)) THEN
            xblk = grid_in4%xblk(xlist4(ireq))
            yblk = grid_in4%yblk(ylist4(ireq))
            xloc = grid_in4%xloc(xlist4(ireq))
            yloc = grid_in4%yloc(ylist4(ireq))

            data_r8_3d_out4(:,ireq) = data_r8_3d_in4%blk(xblk,yblk)%val(:,xloc,yloc)
         ENDIF

      ENDDO

#endif

   END SUBROUTINE aggregation_request_data_multigrid

#ifdef USEMPI

   SUBROUTINE aggregation_worker_done ()

   USE MOD_SPMD_Task

   IMPLICIT NONE

   integer :: smesg(2), iproc, idest

      IF (p_is_worker) THEN
         DO iproc = 0, p_np_io-1
            smesg = (/p_iam_glb, -1/)
            idest = p_address_io(iproc)
            CALL mpi_send (smesg, 2, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)
         ENDDO
      ENDIF

   END SUBROUTINE aggregation_worker_done


   SUBROUTINE aggregation_worker_done_multigrid ()

   USE MOD_SPMD_Task

   IMPLICIT NONE

   integer :: smesg(3), iproc, idest

      IF (p_is_worker) THEN
         DO iproc = 0, p_np_io-1
            smesg = (/p_iam_glb, -1, 1/)
            idest = p_address_io(iproc)
            CALL mpi_send (smesg, 3, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)
          ENDDO
      ENDIF

   END SUBROUTINE

#endif

   SUBROUTINE fillnan (vec, fill, defval)

   USE MOD_Precision
   USE MOD_UserDefFun, only: isnan_ud
   IMPLICIT NONE

   real(r8), intent(inout) :: vec(:)
   logical,  intent(in)    :: fill
   real(r8), intent(in)    :: defval

   ! local variables
   integer  :: i, n
   real(r8) :: s

      n = 0
      s = 0.
      DO i = lbound(vec,1), ubound(vec,1)
         IF (.not. isnan_ud(vec(i))) THEN
            n = n + 1
            s = s + vec(i)
         ENDIF
      ENDDO

      IF ((n > 0) .and. (n < size(vec))) THEN
         s = s/n
         DO i = lbound(vec,1), ubound(vec,1)
            IF (isnan_ud(vec(i))) vec(i) = s
         ENDDO
      ENDIF

      IF ((n == 0) .and. fill) THEN
         vec(:) = defval
      ENDIF

   END SUBROUTINE fillnan

END MODULE MOD_AggregationRequestData
