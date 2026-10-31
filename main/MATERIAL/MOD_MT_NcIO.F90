#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_NcIO
!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!  Restart presence checks for material (methane) state on vector files.
!  CoLM202X added these to share/MOD_NetCDFVector; upstream does not have
!  them, and the host file is left unchanged, so they live here and call
!  the upstream ncio_read_vector for the actual read.
!
!  - ncio_vector_var_present / ncio_vector_group_presence: is a variable
!    present in every block file of a vector restart (all IO processes
!    agree; present in some blocks only is an error);
!  - ncio_read_vector_complete: the same check (plus the shape), then the
!    upstream read; a missing variable is an error unless a default is
!    given and ncio_set_complete_require_present(.false.) is in effect.
!-----------------------------------------------------------------------

   USE MOD_Precision
   IMPLICIT NONE
   PRIVATE

   logical, save :: complete_require_present = .false.

   INTERFACE ncio_read_vector_complete
      MODULE procedure ncio_read_vector_complete_real8_1d
      MODULE procedure ncio_read_vector_complete_real8_2d
   END INTERFACE ncio_read_vector_complete

   PUBLIC :: ncio_read_vector_complete
   PUBLIC :: ncio_set_complete_require_present
   PUBLIC :: ncio_vector_group_presence
   PUBLIC :: ncio_vector_var_present

CONTAINS

   SUBROUTINE ncio_set_complete_require_present (require_present)

   IMPLICIT NONE
   logical, intent(in) :: require_present

      complete_require_present = require_present

   END SUBROUTINE ncio_set_complete_require_present

   SUBROUTINE ncio_require_complete_vector_var (filename, dataname, pixelset, &
         expected_rank, allow_missing, all_present, expected_dim1)

   USE MOD_NetCDFSerial, only: ncio_var_exist, ncio_inquire_varsize
   USE MOD_SPMD_Task
   USE MOD_Block, only: get_filename_block
   USE MOD_Pixelset, only: pixelset_type
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   type(pixelset_type), intent(in) :: pixelset
   integer, intent(in) :: expected_rank
   logical, intent(in) :: allow_missing
   logical, intent(out) :: all_present
   integer, intent(in), optional :: expected_dim1

   integer :: iblkgrp, iblk, jblk
   integer :: counts(3)
   integer, allocatable :: varsize(:)
   character(len=256) :: fileblock
   logical :: block_shape_ok

      counts(:) = 0
      all_present = .false.

      IF (.not. (p_is_io .or. p_is_worker)) RETURN

      IF (p_is_io) THEN
         counts(1) = pixelset%nblkgrp
         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)
            CALL get_filename_block(filename, iblk, jblk, fileblock)
            IF (ncio_var_exist(fileblock, dataname, readflag=.false.)) THEN
               counts(2) = counts(2) + 1
               CALL ncio_inquire_varsize(fileblock, dataname, varsize)

               block_shape_ok = .false.
               IF (allocated(varsize)) THEN
                  IF (size(varsize) == expected_rank) THEN
                     block_shape_ok = &
                        varsize(expected_rank) == pixelset%vecgs%vlen(iblk,jblk)
                     IF (block_shape_ok .and. present(expected_dim1)) &
                        block_shape_ok = varsize(1) == expected_dim1
                  ENDIF
                  deallocate(varsize)
               ENDIF
               IF (.not. block_shape_ok) counts(3) = counts(3) + 1
            ENDIF
         ENDDO
      ENDIF

#ifdef USEMPI
      IF (p_is_io) THEN
         CALL mpi_allreduce(MPI_IN_PLACE, counts, 3, MPI_INTEGER, MPI_SUM, p_comm_io, p_err)
      ENDIF
      CALL mpi_bcast(counts, 3, MPI_INTEGER, p_root, p_comm_group, p_err)
#endif

      IF (counts(2) > 0 .and. counts(2) < counts(1)) THEN
         IF (p_is_io .and. p_iam_io == p_root) WRITE(*,'(A)') &
            'ERROR: restart variable '//trim(dataname)//' is missing from some block files in '//trim(filename)//'.'
         CALL CoLM_stop()
      ENDIF
      IF (counts(2) == 0 .and. .not. allow_missing) THEN
         IF (p_is_io .and. p_iam_io == p_root) WRITE(*,'(A)') &
            'ERROR: required restart variable '//trim(dataname)//' is absent from all block files in '//trim(filename)//'.'
         CALL CoLM_stop()
      ENDIF
      IF (counts(3) > 0) THEN
         IF (p_is_io .and. p_iam_io == p_root) WRITE(*,'(A)') &
            'ERROR: restart variable '//trim(dataname)//' has an incompatible shape in '//trim(filename)//'.'
         CALL CoLM_stop()
      ENDIF
      all_present = counts(1) > 0 .and. counts(2) == counts(1)

   END SUBROUTINE ncio_require_complete_vector_var

   logical FUNCTION ncio_vector_var_present (filename, dataname, pixelset)

   USE MOD_NetCDFSerial, only: ncio_var_exist
   USE MOD_SPMD_Task
   USE MOD_Block, only: get_filename_block
   USE MOD_Pixelset, only: pixelset_type
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   type(pixelset_type), intent(in) :: pixelset

   integer :: iblkgrp, iblk, jblk
   integer :: counts(2), global_counts(2)
   character(len=256) :: fileblock

      counts(:) = 0
      IF (p_is_io) THEN
         counts(1) = pixelset%nblkgrp
         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)
            CALL get_filename_block(filename, iblk, jblk, fileblock)
            IF (ncio_var_exist(fileblock, dataname, readflag=.false.)) &
               counts(2) = counts(2) + 1
         ENDDO
      ENDIF

#ifdef USEMPI
      IF (p_is_io) THEN
         CALL mpi_allreduce(MPI_IN_PLACE, counts, 2, MPI_INTEGER, MPI_SUM, p_comm_io, p_err)
      ENDIF
      global_counts(:) = 0
      IF (p_is_io) THEN
         IF (p_iam_io == p_root) global_counts(:) = counts(:)
      ENDIF
      CALL mpi_allreduce(MPI_IN_PLACE, global_counts, 2, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#else
      global_counts(:) = counts(:)
#endif

      IF (global_counts(2) > 0 .and. global_counts(2) < global_counts(1)) THEN
         IF (p_is_master) WRITE(*,'(A)') &
            'ERROR: restart variable '//trim(dataname)//' is missing from some block files in '//trim(filename)//'.'
         CALL CoLM_stop()
      ENDIF
      ncio_vector_var_present = global_counts(1) > 0 .and. global_counts(2) == global_counts(1)

   END FUNCTION ncio_vector_var_present

   SUBROUTINE ncio_vector_group_presence (filename, datanames, pixelset, present_flags)

   USE MOD_NetCDFSerial, only: ncio_var_exist
   USE MOD_SPMD_Task
   USE MOD_Block, only: get_filename_block
   USE MOD_Pixelset, only: pixelset_type
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: datanames(:)
   type(pixelset_type), intent(in) :: pixelset
   logical, intent(out) :: present_flags(size(datanames))

   integer :: iblkgrp, iblk, jblk, ivar
   integer, allocatable :: counts(:), global_counts(:)
   character(len=256) :: fileblock

      allocate(counts(0:size(datanames)), global_counts(0:size(datanames)))
      counts(:) = 0
      IF (p_is_io) THEN
         counts(0) = pixelset%nblkgrp
         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)
            CALL get_filename_block(filename, iblk, jblk, fileblock)
            DO ivar = 1, size(datanames)
               IF (ncio_var_exist(fileblock, trim(datanames(ivar)), readflag=.false.)) &
                  counts(ivar) = counts(ivar) + 1
            ENDDO
         ENDDO
      ENDIF

#ifdef USEMPI
      IF (p_is_io) THEN
         CALL mpi_allreduce(MPI_IN_PLACE, counts, size(counts), MPI_INTEGER, MPI_SUM, p_comm_io, p_err)
      ENDIF
      global_counts(:) = 0
      IF (p_is_io) THEN
         IF (p_iam_io == p_root) global_counts(:) = counts(:)
      ENDIF
      CALL mpi_allreduce(MPI_IN_PLACE, global_counts, size(global_counts), &
         MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#else
      global_counts(:) = counts(:)
#endif

      DO ivar = 1, size(datanames)
         IF (global_counts(ivar) > 0 .and. global_counts(ivar) < global_counts(0)) THEN
            IF (p_is_master) WRITE(*,'(A,A,A)') 'ERROR: restart feature-group field ', &
               trim(datanames(ivar)), ' is missing from some vector blocks.'
            CALL CoLM_stop()
         ENDIF
         present_flags(ivar) = global_counts(0) > 0 .and. &
            global_counts(ivar) == global_counts(0)
      ENDDO

      deallocate(counts, global_counts)
   END SUBROUTINE ncio_vector_group_presence

   SUBROUTINE ncio_read_vector_complete_real8_1d (filename, dataname, pixelset, rdata, defval)

   USE MOD_Precision
   USE MOD_Pixelset, only: pixelset_type
   USE MOD_NetCDFVector, only: ncio_read_vector
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   type(pixelset_type), intent(in) :: pixelset
   real(r8), allocatable, intent(inout) :: rdata(:)
   real(r8), intent(in), optional :: defval
   logical :: field_present, allow_field_missing

      allow_field_missing = present(defval) .and. (.not. complete_require_present)
      CALL ncio_require_complete_vector_var(filename, dataname, pixelset, &
         expected_rank=1, allow_missing=allow_field_missing, all_present=field_present)
      IF (present(defval)) THEN
         CALL ncio_read_vector (filename, dataname, pixelset, rdata, defval)
      ELSE
         CALL ncio_read_vector (filename, dataname, pixelset, rdata)
      ENDIF

   END SUBROUTINE ncio_read_vector_complete_real8_1d

   SUBROUTINE ncio_read_vector_complete_real8_2d (filename, dataname, ndim1, pixelset, rdata, defval)

   USE MOD_Precision
   USE MOD_Pixelset, only: pixelset_type
   USE MOD_NetCDFVector, only: ncio_read_vector
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   integer, intent(in) :: ndim1
   type(pixelset_type), intent(in) :: pixelset
   real(r8), allocatable, intent(inout) :: rdata(:,:)
   real(r8), intent(in), optional :: defval
   logical :: field_present, allow_field_missing

      allow_field_missing = present(defval) .and. (.not. complete_require_present)
      CALL ncio_require_complete_vector_var(filename, dataname, pixelset, &
         expected_rank=2, allow_missing=allow_field_missing, all_present=field_present, &
         expected_dim1=ndim1)
      IF (present(defval)) THEN
         CALL ncio_read_vector (filename, dataname, ndim1, pixelset, rdata, defval)
      ELSE
         CALL ncio_read_vector (filename, dataname, ndim1, pixelset, rdata)
      ENDIF

   END SUBROUTINE ncio_read_vector_complete_real8_2d

END MODULE MOD_MT_NcIO
#endif
