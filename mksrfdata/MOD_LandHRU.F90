#include <define.h>

#ifdef CATCHMENT

MODULE MOD_LandHRU

!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!    Build pixelset "landhru".
!
!    In CoLM, the global/regional area is divided into a hierarchical structure:
!    1. If GRIDBASED or UNSTRUCTURED is defined, it is
!       ELEMENT >>> PATCH
!    2. If CATCHMENT is defined, it is
!       ELEMENT >>> HRU >>> PATCH
!    If Plant Function Type classification is used, PATCH is further divided into PFT.
!    If Plant Community classification is used,     PATCH is further divided into PC.
!
!    "landhru" refers to pixelset HRU.
!
!  Created by Shupeng Zhang, May 2023
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Pixelset
   USE MOD_Grid
   IMPLICIT NONE

   ! ---- Instance ----
   integer :: numhru
   type(grid_type)     :: grid_hru
   type(pixelset_type) :: landhru

   type(subset_type) :: elm_hru

CONTAINS

   ! -------------------------------
   SUBROUTINE landhru_build ()

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_NetcdfSerial
   USE MOD_Utils
   USE MOD_Block
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_Mesh
   USE MOD_LandElm
   USE MOD_CatchmentDataReadin
   USE MOD_Namelist
   USE MOD_AggregationRequestData

   IMPLICIT NONE

   ! Local Variables
   type (block_data_int32_2d) :: hrudata
   integer :: iwork, ncat, nhru, ie, typsgn, npxl, ipxl
   integer*8, allocatable :: catnum(:)
   integer,   allocatable :: numhru_all_g(:), lakeid(:)
   integer,   allocatable :: types(:), order(:), ibuff(:)
   integer :: nhru_glb

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN
         write(*,'(A)') 'Making land hydro units:'
      ENDIF

      IF (p_is_master) THEN
         CALL ncio_read_serial (DEF_CatchmentMesh_data, 'basin_numhru', numhru_all_g)
         CALL ncio_read_serial (DEF_CatchmentMesh_data, 'lake_id', lakeid)
      ENDIF

#ifdef USEMPI
      IF (p_is_master) THEN
         DO iwork = 0, p_np_worker-1

            CALL mpi_recv (ncat, 1, MPI_INTEGER, p_address_worker(iwork), mpi_tag_size, &
               p_comm_glb, p_stat, p_err)

            IF (ncat > 0) THEN
               allocate (catnum(ncat))
               allocate (ibuff (ncat))

               CALL mpi_recv (catnum, ncat, MPI_INTEGER8, p_address_worker(iwork), mpi_tag_data, &
                  p_comm_glb, p_stat, p_err)

               nhru = sum(numhru_all_g(catnum))
               CALL mpi_send (nhru, 1, MPI_INTEGER, &
                  p_address_worker(iwork), mpi_tag_size, p_comm_glb, p_err)

               ibuff = lakeid(catnum)
               CALL mpi_send (ibuff, ncat, MPI_INTEGER, &
                  p_address_worker(iwork), mpi_tag_data, p_comm_glb, p_err)

               deallocate(catnum)
               deallocate(ibuff )
            ENDIF
         ENDDO
      ENDIF

      IF (p_is_worker) THEN
         CALL mpi_send (numelm, 1, MPI_INTEGER, p_address_master, mpi_tag_size, p_comm_glb, p_err)
         IF (numelm > 0) THEN
            allocate (lakeid (numelm))
            CALL mpi_send (landelm%eindex, numelm, MPI_INTEGER8, p_address_master, mpi_tag_data, p_comm_glb, p_err)
            CALL mpi_recv (numhru, 1,      MPI_INTEGER, p_address_master, mpi_tag_size, p_comm_glb, p_stat, p_err)
            CALL mpi_recv (lakeid, numelm, MPI_INTEGER, p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)
         ELSE
            numhru = 0
         ENDIF
      ENDIF
#else
      numhru = sum(numhru_all_g)
#endif

      IF (p_is_master) THEN
         IF (allocated(numhru_all_g)) deallocate(numhru_all_g)
      ENDIF

      IF (p_is_io) CALL allocate_block_data (grid_hru, hrudata)
      CALL catchment_data_read (DEF_CatchmentMesh_data, 'ihydrounit2d', grid_hru, hrudata)

#ifdef USEMPI
      IF (p_is_io) THEN
         CALL aggregation_data_daemon (grid_hru, data_i4_2d_in1 = hrudata)
      ENDIF
#endif

      IF (p_is_worker) THEN

         IF (numhru > 0) THEN
            allocate (landhru%eindex (numhru))
            allocate (landhru%settyp (numhru))
            allocate (landhru%ipxstt (numhru))
            allocate (landhru%ipxend (numhru))
            allocate (landhru%ielm   (numhru))
         ENDIF

         numhru = 0

         DO ie = 1, numelm

            IF (lakeid(ie) > 0) THEN
               typsgn = -1
            ELSE
               typsgn = 1
            ENDIF

            npxl = mesh(ie)%npxl

            allocate (types (1:npxl))

            CALL aggregation_request_data (landelm, ie, grid_hru, zip = .false., &
               data_i4_2d_in1 = hrudata, data_i4_2d_out1 = ibuff)

            types = ibuff

            allocate (order (1:npxl))
            order = (/ (ipxl, ipxl = 1, npxl) /)

            CALL quicksort (npxl, types, order)

            mesh(ie)%ilon(1:npxl) = mesh(ie)%ilon(order)
            mesh(ie)%ilat(1:npxl) = mesh(ie)%ilat(order)

            DO ipxl = 1, npxl
               IF (ipxl == 1) THEN
                  numhru = numhru + 1
                  landhru%eindex (numhru) = mesh(ie)%indx
                  landhru%settyp (numhru) = types(ipxl) * typsgn
                  landhru%ipxstt (numhru) = ipxl
                  landhru%ielm   (numhru) = ie
               ELSEIF (types(ipxl) /= types(ipxl-1)) THEN
                  landhru%ipxend(numhru) = ipxl - 1

                  numhru = numhru + 1
                  landhru%eindex (numhru) = mesh(ie)%indx
                  landhru%settyp (numhru) = types(ipxl) * typsgn
                  landhru%ipxstt (numhru) = ipxl
                  landhru%ielm   (numhru) = ie
               ENDIF
            ENDDO
            landhru%ipxend(numhru) = npxl

            deallocate (ibuff)
            deallocate (types)
            deallocate (order)

         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

      landhru%nset = numhru
      CALL landhru%set_vecgs

#ifdef USEMPI
      IF (p_is_worker) THEN
         CALL mpi_reduce (numhru, nhru_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', nhru_glb, ' hydro units.'
         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      write(*,'(A,I12,A)') 'Total: ', numhru, ' hydro units.'
#endif

      IF (allocated(lakeid)) deallocate(lakeid)

   END SUBROUTINE landhru_build

END MODULE MOD_LandHRU
#endif
