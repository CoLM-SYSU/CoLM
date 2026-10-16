#include <define.h>

#if (defined UNSTRUCTURED || defined CATCHMENT)
MODULE MOD_HistVector

!----------------------------------------------------------------------------
! !DESCRIPTION:
!
!     Write out vectorized model results to history files.
!
!  Created by Shupeng Zhang, May 2023
!
!  TODO...(need complement)
!----------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Vars_Global, only: spval
   USE MOD_Mesh
   USE MOD_LandElm
#ifdef CATCHMENT
   USE MOD_LandHRU
#endif
   USE MOD_Pixelset
   USE MOD_NetCDFSerial
#ifdef CATCHMENT
   USE MOD_HRUVector
#else
   USE MOD_ElmVector
#endif

CONTAINS

   ! -- write history time --
   SUBROUTINE hist_vector_write_time (filename, filelast, dataname, time, itime_in_file)

   IMPLICIT NONE

   character (len=*), intent(in) :: filename
   character (len=*), intent(in) :: filelast
   character (len=*), intent(in) :: dataname
   integer, intent(in)  :: time(3)
   integer, intent(out) :: itime_in_file

   ! Local Variables
   logical :: fexists

      IF (p_is_master) THEN

         inquire (file=filename, exist=fexists)
         IF ((.not. fexists) .or. (trim(filename) /= trim(filelast))) THEN
            CALL ncio_create_file (trim(filename))
            CALL ncio_define_dimension(filename, 'time', 0)

#ifdef CATCHMENT
            CALL ncio_define_dimension(filename, 'hydrounit', totalnumhru)

            CALL ncio_write_serial (filename, 'bsn_hru', eindx_hru, 'hydrounit')
            CALL ncio_put_attr (filename, 'bsn_hru', 'long_name', &
               'basin index of hydrological units in mesh')

            CALL ncio_write_serial (filename, 'typ_hru' , htype_hru, 'hydrounit')
            CALL ncio_put_attr (filename, 'typ_hru' , 'long_name', &
               'index of hydrological units inside basin')
#else
            CALL ncio_define_dimension(filename, 'element', totalnumelm)
            CALL ncio_write_serial (filename, 'elmindex', eindex_glb, 'element')
            CALL ncio_put_attr (filename, 'elmindex', 'long_name', &
               'element index in mesh')
#endif

            CALL ncio_write_colm_dimension (filename)

         ENDIF

         CALL ncio_write_time (filename, dataname, time, itime_in_file, DEF_HIST_FREQ)

      ENDIF

   END SUBROUTINE hist_vector_write_time


   SUBROUTINE aggregate_to_vector_and_write_2d ( &
         acc_vec_patch, file_hist, varname, itime_in_file, filter, &
         longname, units, input_mode)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_LandPatch
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   real(r8), intent(in) :: acc_vec_patch (:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer,          intent(in) :: itime_in_file
   logical, intent(in) :: filter(:)

   character(len=*), intent(in) :: longname
   character(len=*), intent(in) :: units

   character(len=*), intent(in), optional :: input_mode

   ! Local variables
   integer :: numset, totalnumset, iset, istt, iend, iwork, mesg(2), isrc, ndata, compress
   logical,  allocatable :: mask(:)
   real(r8), allocatable :: frac(:)
   real(r8), allocatable :: acc_vec(:), rcache(:)
   real(r8) :: sumwt
   character(len=256) :: inmode

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_worker) THEN
#ifdef CATCHMENT
         numset = numhru
#else
         numset = numelm
#endif

         inmode = 'average'
         IF (present(input_mode)) inmode = trim(input_mode)

         IF (numset > 0) THEN

            allocate (acc_vec (numset))
            acc_vec(:) = spval

            DO iset = 1, numset
#ifdef CATCHMENT
               istt = hru_patch%substt(iset)
               iend = hru_patch%subend(iset)
#else
               istt = elm_patch%substt(iset)
               iend = elm_patch%subend(iset)
#endif

               IF ((istt > 0) .and. (iend >= istt)) THEN
                  allocate (mask(istt:iend))
                  allocate (frac(istt:iend))
                  mask = (acc_vec_patch(istt:iend) /= spval) .and. filter(istt:iend)
                  IF (any(mask)) THEN
                     IF (trim(inmode) == 'average') THEN
#ifdef CATCHMENT
                        frac = hru_patch%subfrc(istt:iend)
#else
                        frac = elm_patch%subfrc(istt:iend)
#endif
                        sumwt = sum(frac, mask = mask)
                        acc_vec(iset) = sum(frac * acc_vec_patch(istt:iend), mask = mask)
                        acc_vec(iset) = acc_vec(iset) / sumwt
                     ELSE
                        acc_vec(iset) = sum(acc_vec_patch(istt:iend), mask = mask)
                     ENDIF
                  ENDIF
                  deallocate(mask)
                  deallocate(frac)
               ENDIF
            ENDDO
         ENDIF

#ifdef USEMPI
         mesg = (/p_iam_glb, numset/)
         CALL mpi_send (mesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         IF (numset > 0) THEN
            CALL mpi_send (acc_vec, numset, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_err)
         ENDIF
#endif
      ENDIF

      IF (p_is_master) THEN

#ifdef CATCHMENT
         totalnumset = totalnumhru
#else
         totalnumset = totalnumelm
#endif

         IF (.not. allocated(acc_vec)) THEN
            allocate (acc_vec (totalnumset))
         ENDIF

#ifdef USEMPI
         DO iwork = 0, p_np_worker-1
            CALL mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               allocate(rcache (ndata))
               CALL mpi_recv (rcache, ndata, MPI_REAL8, isrc, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)

#ifdef CATCHMENT
               acc_vec(hru_data_address(p_itis_worker(isrc))%val) = rcache
#else
               acc_vec(elm_data_address(p_itis_worker(isrc))%val) = rcache
#endif

               deallocate (rcache)
            ENDIF
         ENDDO
#else
#ifdef CATCHMENT
         acc_vec(hru_data_address(0)%val) = acc_vec
#else
         acc_vec(elm_data_address(0)%val) = acc_vec
#endif
#endif
      ENDIF

      IF (p_is_master) THEN

         compress = DEF_HIST_CompressLevel

         IF (itime_in_file >= 1) THEN
#ifdef CATCHMENT
            CALL ncio_write_serial_time (file_hist, varname, itime_in_file, acc_vec, &
               'hydrounit', 'time', compress)
#else
            CALL ncio_write_serial_time (file_hist, varname, itime_in_file, acc_vec, &
               'element', 'time', compress)
#endif
         ELSE
#ifdef CATCHMENT
            CALL ncio_write_serial (file_hist, varname, acc_vec, 'hydrounit', compress)
#else
            CALL ncio_write_serial (file_hist, varname, acc_vec, 'element', compress)
#endif
         ENDIF

         IF (itime_in_file <= 1) THEN
            CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
            CALL ncio_put_attr (file_hist, varname, 'units', units)
            CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
         ENDIF

      ENDIF

      IF (allocated(acc_vec)) deallocate (acc_vec)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE aggregate_to_vector_and_write_2d


   SUBROUTINE aggregate_to_vector_and_write_3d ( &
         acc_vec_patch, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, filter, &
         longname, units)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_LandPatch
   USE MOD_Vars_1DAccFluxes,  only: nac
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   real(r8), intent(in) :: acc_vec_patch (lb1:,:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer,          intent(in) :: itime_in_file
   character(len=*), intent(in) :: dim1name
   integer,          intent(in) :: lb1, ndim1

   logical, intent(in) :: filter(:)

   character(len=*), intent(in) :: longname
   character(len=*), intent(in) :: units

   ! Local variables
   integer :: numset, totalnumset, iset, istt, iend, iwork, mesg(2), isrc, ndata, compress
   integer :: ub1, i1
   logical,  allocatable :: mask(:)
   real(r8), allocatable :: frac(:)
   real(r8), allocatable :: acc_vec(:,:), rcache(:,:)
   real(r8) :: sumwt

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      ub1 = lb1 + ndim1 - 1

      IF (p_is_worker) THEN
#ifdef CATCHMENT
         numset = numhru
#else
         numset = numelm
#endif

         IF (numset > 0) THEN

            allocate (acc_vec (lb1:ub1,numset))

            acc_vec(:,:) = spval

            DO iset = 1, numset
#ifdef CATCHMENT
               istt = hru_patch%substt(iset)
               iend = hru_patch%subend(iset)
#else
               istt = elm_patch%substt(iset)
               iend = elm_patch%subend(iset)
#endif

               IF ((istt > 0) .and. (iend >= istt)) THEN
                  allocate (mask(istt:iend))
                  allocate (frac(istt:iend))
                  DO i1 = lb1, ub1
                     mask = (acc_vec_patch(i1,istt:iend) /= spval) .and. filter(istt:iend)
                     IF (any(mask)) THEN
#ifdef CATCHMENT
                        frac = hru_patch%subfrc(istt:iend)
#else
                        frac = elm_patch%subfrc(istt:iend)
#endif
                        sumwt = sum(frac, mask = mask)
                        acc_vec(i1,iset) = sum(frac * acc_vec_patch(i1,istt:iend), mask = mask)
                        acc_vec(i1,iset) = acc_vec(i1,iset) / sumwt / nac
                     ENDIF
                  ENDDO
                  deallocate(mask)
                  deallocate(frac)
               ENDIF
            ENDDO
         ENDIF

#ifdef USEMPI
         mesg = (/p_iam_glb, numset/)
         CALL mpi_send (mesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         IF (numset > 0) THEN
            CALL mpi_send (acc_vec, ndim1 * numset, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_err)
         ENDIF
#endif
      ENDIF

      IF (p_is_master) THEN

#ifdef CATCHMENT
         totalnumset = totalnumhru
#else
         totalnumset = totalnumelm
#endif

         IF (.not. allocated(acc_vec)) THEN
            allocate (acc_vec (ndim1,totalnumset))
         ENDIF

#ifdef USEMPI
         DO iwork = 0, p_np_worker-1
            CALL mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               allocate(rcache (ndim1,ndata))
               CALL mpi_recv (rcache, ndim1*ndata, MPI_REAL8, isrc, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)

               DO i1 = 1, ndim1
#ifdef CATCHMENT
                  acc_vec(i1,hru_data_address(p_itis_worker(isrc))%val) = rcache(i1,:)
#else
                  acc_vec(i1,elm_data_address(p_itis_worker(isrc))%val) = rcache(i1,:)
#endif

               ENDDO

               deallocate (rcache)
            ENDIF
         ENDDO
#else
         DO i1 = lb1, ub1
#ifdef CATCHMENT
            acc_vec(i1,hru_data_address(0)%val) = acc_vec(i1,:)
#else
            acc_vec(i1,elm_data_address(0)%val) = acc_vec(i1,:)
#endif
         ENDDO
#endif
      ENDIF

      IF (p_is_master) THEN

         CALL ncio_define_dimension (file_hist, dim1name, ndim1)

         compress = DEF_HIST_CompressLevel

         IF (itime_in_file >= 1) THEN
#ifdef CATCHMENT
            CALL ncio_write_serial_time (file_hist, varname, itime_in_file, acc_vec, &
               dim1name, 'hydrounit', 'time', compress)
#else
            CALL ncio_write_serial_time (file_hist, varname, itime_in_file, acc_vec, &
               dim1name, 'element', 'time', compress)
#endif
         ELSE
#ifdef CATCHMENT
            CALL ncio_write_serial (file_hist, varname, acc_vec, &
               dim1name, 'hydrounit', compress)
#else
            CALL ncio_write_serial (file_hist, varname, acc_vec, &
               dim1name, 'element', compress)
#endif
         ENDIF

         IF (itime_in_file <= 1) THEN
            CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
            CALL ncio_put_attr (file_hist, varname, 'units', units)
            CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
         ENDIF

      ENDIF

      IF (allocated(acc_vec)) deallocate (acc_vec)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE aggregate_to_vector_and_write_3d


   SUBROUTINE aggregate_to_vector_and_write_4d ( &
         acc_vec_patch, file_hist, varname, itime_in_file,   &
         dim1name, lb1, ndim1, dim2name, lb2, ndim2, filter, &
         longname, units)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_LandPatch
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   real(r8), intent(in) :: acc_vec_patch (lb1:,lb2:,:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer,          intent(in) :: itime_in_file
   character(len=*), intent(in) :: dim1name, dim2name
   integer,          intent(in) :: lb1, ndim1, lb2, ndim2

   logical, intent(in) :: filter(:)

   character(len=*), intent(in) :: longname
   character(len=*), intent(in) :: units

   ! Local variables
   integer :: numset, totalnumset, iset, istt, iend, iwork, mesg(2), isrc, ndata, compress
   integer :: ub1, i1, ub2, i2
   logical,  allocatable :: mask(:)
   real(r8), allocatable :: frac(:)
   real(r8), allocatable :: acc_vec(:,:,:), rcache(:,:,:)
   real(r8) :: sumwt

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      ub1 = lb1 + ndim1 - 1
      ub2 = lb2 + ndim2 - 1

      IF (p_is_worker) THEN
#ifdef CATCHMENT
         numset = numhru
#else
         numset = numelm
#endif

         IF (numset > 0) THEN

            allocate (acc_vec (lb1:ub1,lb2:ub2,numset))

            acc_vec(:,:,:) = spval

            DO iset = 1, numset
#ifdef CATCHMENT
               istt = hru_patch%substt(iset)
               iend = hru_patch%subend(iset)
#else
               istt = elm_patch%substt(iset)
               iend = elm_patch%subend(iset)
#endif

               IF ((istt > 0) .and. (iend >= istt)) THEN
                  allocate (mask(istt:iend))
                  allocate (frac(istt:iend))
                  DO i1 = lb1, ub1
                     DO i2 = lb2, ub2
                        mask = (acc_vec_patch(i1,i2,istt:iend) /= spval) .and. filter(istt:iend)
                        IF (any(mask)) THEN
#ifdef CATCHMENT
                           frac = hru_patch%subfrc(istt:iend)
#else
                           frac = elm_patch%subfrc(istt:iend)
#endif
                           sumwt = sum(frac, mask = mask)
                           acc_vec(i1,i2,iset) = sum(frac * acc_vec_patch(i1,i2,istt:iend), mask = mask)
                           acc_vec(i1,i2,iset) = acc_vec(i1,i2,iset) / sumwt
                        ENDIF
                     ENDDO
                  ENDDO
                  deallocate(mask)
                  deallocate(frac)
               ENDIF
            ENDDO
         ENDIF

#ifdef USEMPI
         mesg = (/p_iam_glb, numset/)
         CALL mpi_send (mesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         IF (numset > 0) THEN
            CALL mpi_send (acc_vec, ndim1 * ndim2 * numset, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_err)
         ENDIF
#endif
      ENDIF

      IF (p_is_master) THEN

#ifdef CATCHMENT
         totalnumset = totalnumhru
#else
         totalnumset = totalnumelm
#endif

         IF (.not. allocated(acc_vec)) THEN
            allocate (acc_vec (ndim1,ndim2,totalnumset))
         ENDIF

#ifdef USEMPI
         DO iwork = 0, p_np_worker-1
            CALL mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               allocate(rcache (ndim1,ndim2,ndata))
               CALL mpi_recv (rcache, ndim1 * ndim2 * ndata, MPI_REAL8, isrc, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)

               DO i1 = 1, ndim1
                  DO i2 = 1, ndim2
#ifdef CATCHMENT
                     acc_vec(i1,i2,hru_data_address(p_itis_worker(isrc))%val) = rcache(i1,i2,:)
#else
                     acc_vec(i1,i2,elm_data_address(p_itis_worker(isrc))%val) = rcache(i1,i2,:)
#endif
                  ENDDO
               ENDDO

               deallocate (rcache)
            ENDIF
         ENDDO
#else
         DO i1 = lb1, ub1
            DO i2 = lb2, ub2
#ifdef CATCHMENT
               acc_vec(i1,i2,hru_data_address(0)%val) = acc_vec(i1,i2,:)
#else
               acc_vec(i1,i2,elm_data_address(0)%val) = acc_vec(i1,i2,:)
#endif
            ENDDO
         ENDDO
#endif
      ENDIF

      IF (p_is_master) THEN

         CALL ncio_define_dimension (file_hist, dim1name, ndim1)
         CALL ncio_define_dimension (file_hist, dim2name, ndim2)

         compress = DEF_HIST_CompressLevel

         IF (itime_in_file >= 1) THEN
#ifdef CATCHMENT
            CALL ncio_write_serial_time (file_hist, varname, itime_in_file, acc_vec, &
               dim1name, dim2name, 'hydrounit', 'time', compress)
#else
            CALL ncio_write_serial_time (file_hist, varname, itime_in_file, acc_vec, &
               dim1name, dim2name, 'element', 'time', compress)
#endif
         ELSE
#ifdef CATCHMENT
            CALL ncio_write_serial (file_hist, varname, acc_vec, &
               dim1name, dim2name, 'hydrounit', compress)
#else
            CALL ncio_write_serial (file_hist, varname, acc_vec, &
               dim1name, dim2name, 'element', compress)
#endif
         ENDIF

         IF (itime_in_file <= 1) THEN
            CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
            CALL ncio_put_attr (file_hist, varname, 'units', units)
            CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
         ENDIF

      ENDIF

      IF (allocated(acc_vec)) deallocate (acc_vec)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE aggregate_to_vector_and_write_4d


   SUBROUTINE aggregate_to_vector_and_write_urb_2d ( &
         acc_vec_patch, file_hist, varname, itime_in_file, filter, &
         longname, units, input_mode)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_LandPatch
   USE MOD_LandUrban
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   real(r8), intent(in) :: acc_vec_patch (:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer,          intent(in) :: itime_in_file
   logical, intent(in) :: filter(:)

   character(len=*), intent(in) :: longname
   character(len=*), intent(in) :: units

   character(len=*), intent(in), optional :: input_mode

   ! Local variables
   integer :: u, ipth, urb_s, urb_e, numpth
   integer :: numset, totalnumset, iset, iwork, mesg(2), isrc, ndata, compress
   logical,  allocatable :: mask(:)
   integer , allocatable :: locpth(:)
   real(r8), allocatable :: frac(:)
   real(r8), allocatable :: acc_vec(:), rcache(:)
   real(r8) :: sumwt
   character(len=256) :: inmode

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_worker) THEN
#ifdef CATCHMENT
         numset = numhru
#else
         numset = numelm
#endif

         inmode = 'average'
         IF (present(input_mode)) inmode = trim(input_mode)

         IF (numset > 0) THEN

            allocate (acc_vec (numset))
            acc_vec(:) = spval

            DO iset = 1, numset

               numpth = count(landurban%eindex==landelm%eindex(iset))
               IF (allocated(locpth)) deallocate(locpth)
               allocate(locpth(numpth))

               locpth = pack([(ipth, ipth=1, numurban)], &
                        landurban%eindex==landelm%eindex(iset))

               urb_s = minval(locpth)
               urb_e = maxval(locpth)

               IF ((urb_s > 0) .and. (urb_e >= urb_s)) THEN
                  allocate (mask(urb_s:urb_e))
                  allocate (frac(urb_s:urb_e))
                  mask = (acc_vec_patch(urb_s:urb_e) /= spval) .and. filter(urb_s:urb_e)
                  IF (any(mask)) THEN
                     IF (trim(inmode) == 'average') THEN
#ifdef CATCHMENT
                        frac = hru_patch%subfrc(urban2patch(urb_s:urb_e))
#else
                        frac = elm_patch%subfrc(urban2patch(urb_s:urb_e))
#endif
                        sumwt = sum(frac, mask = mask)
                        acc_vec(iset) = sum(frac * acc_vec_patch(urb_s:urb_e), mask = mask)
                        acc_vec(iset) = acc_vec(iset) / sumwt
                     ELSE
                        acc_vec(iset) = sum(acc_vec_patch(urb_s:urb_e), mask = mask)
                     ENDIF
                  ENDIF
                  deallocate(mask)
                  deallocate(frac)
               ENDIF
            ENDDO
         ENDIF

#ifdef USEMPI
         mesg = (/p_iam_glb, numset/)
         CALL mpi_send (mesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         IF (numset > 0) THEN
            CALL mpi_send (acc_vec, numset, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_err)
         ENDIF
#endif
      ENDIF

      IF (p_is_master) THEN

#ifdef CATCHMENT
         totalnumset = totalnumhru
#else
         totalnumset = totalnumelm
#endif

         IF (.not. allocated(acc_vec)) THEN
            allocate (acc_vec (totalnumset))
         ENDIF

#ifdef USEMPI
         DO iwork = 0, p_np_worker-1
            CALL mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               allocate(rcache (ndata))
               CALL mpi_recv (rcache, ndata, MPI_REAL8, isrc, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)

#ifdef CATCHMENT
               acc_vec(hru_data_address(p_itis_worker(isrc))%val) = rcache
#else
               acc_vec(elm_data_address(p_itis_worker(isrc))%val) = rcache
#endif

               deallocate (rcache)
            ENDIF
         ENDDO
#else
#ifdef CATCHMENT
         acc_vec(hru_data_address(0)%val) = acc_vec
#else
         acc_vec(elm_data_address(0)%val) = acc_vec
#endif
#endif
      ENDIF

      IF (p_is_master) THEN

         compress = DEF_HIST_CompressLevel

         IF (itime_in_file >= 1) THEN
#ifdef CATCHMENT
            CALL ncio_write_serial_time (file_hist, varname, itime_in_file, acc_vec, &
               'hydrounit', 'time', compress)
#else
            CALL ncio_write_serial_time (file_hist, varname, itime_in_file, acc_vec, &
               'element', 'time', compress)
#endif
         ELSE
#ifdef CATCHMENT
            CALL ncio_write_serial (file_hist, varname, acc_vec, 'hydrounit', compress)
#else
            CALL ncio_write_serial (file_hist, varname, acc_vec, 'element', compress)
#endif
         ENDIF

         IF (itime_in_file <= 1) THEN
            CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
            CALL ncio_put_attr (file_hist, varname, 'units', units)
            CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
         ENDIF

      ENDIF

      IF (allocated(acc_vec)) deallocate (acc_vec)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE aggregate_to_vector_and_write_urb_2d

END MODULE MOD_HistVector
#endif
