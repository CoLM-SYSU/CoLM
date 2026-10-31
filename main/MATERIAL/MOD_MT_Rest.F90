#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_Rest

   USE MOD_Precision
   USE MOD_MT_Defs, only: nmts, material_init_water_ratio, mt_tiny, &
      mt_water_min_for_ratio, materials, material_uses_delta_diagnostics, mt_uses_generic_water, &
      material_is_nonvolatile_solute, material_is_isotope, material_aquifer_isotope_state_valid, &
      material_build_descriptor_identity, MATERIAL_DESCRIPTOR_IDENTITY_WIDTH
   USE MOD_MT_Solute, only: mt_equilibrate_dissolved
   USE MOD_MT_Vars
   USE MOD_LandPatch, only: landpatch
   USE MOD_Vars_TimeInvariants, only: patchtype
   USE MOD_Block, only: get_filename_block
   USE MOD_Namelist, only: DEF_REST_CompressLevel, DEF_USE_VariablySaturatedFlow
   USE MOD_MT_Namelist, only: DEF_MT_AQUIFER_MIXING_WATER_MM
   USE MOD_NetCDFSerial, only: ncio_var_exist, ncio_inquire_varsize, ncio_read_serial, &
      ncio_write_serial, ncio_define_dimension
   USE MOD_NetCDFVector, only: ncio_read_vector, ncio_write_vector
   USE MOD_MT_Lifecycle, only: material_lifecycle_land_write_restart
   USE MOD_MT_Forcing, only: material_forcing_write_restart
   USE MOD_SPMD_Task, only: CoLM_stop

   IMPLICIT NONE

   integer, parameter :: LAND_MT_RESTART_SCHEMA_VERSION = 5
   integer, parameter :: LAND_MT_DESCRIPTOR_FIELDS = MATERIAL_DESCRIPTOR_IDENTITY_WIDTH
   real(r8), parameter :: LAND_MT_RESTART_NEGATIVE_DUST = 1.0e-12_r8

   PRIVATE :: material_dim_matches, land_mt_descriptor_matches, &
      write_land_mt_descriptor_metadata, write_land_mt_transaction_marker, &
      land_transport_mt_count, read_transport_patch_field, &
      read_transport_soilsnow_field, pack_transport_patch, &
      pack_transport_soilsnow, validate_land_mt_restart_state

CONTAINS

   SUBROUTINE mt_write_host_restart (file_restart, maxsnl, nl_soil, numpatch, &
      ldew_rain, ldew_snow, wliq_soisno, wice_soisno, wa, wdsrf, wetwat, scv, &
      compress, waterstorage)
      USE MOD_SPMD_Task, only: p_is_worker
      character(len=*), intent(in) :: file_restart
      integer, intent(in) :: maxsnl, nl_soil, numpatch, compress
      real(r8), allocatable, intent(in) :: ldew_rain(:), ldew_snow(:), wa(:), wdsrf(:), wetwat(:), scv(:)
      real(r8), allocatable, intent(in) :: wliq_soisno(:,:), wice_soisno(:,:), waterstorage(:)
      real(r8) :: empty_patch(0), empty_soilsnow(maxsnl+1:nl_soil,0)

         IF (allocated(ldew_rain)) THEN
            IF (.not. p_is_worker) CALL CoLM_stop ('material restart patch state on non-worker')
            IF (numpatch < 0) CALL CoLM_stop ('invalid material restart patch count')
            IF (.not. allocated(ldew_snow) .or. .not. allocated(wliq_soisno) .or. &
                .not. allocated(wice_soisno) .or. .not. allocated(wa) .or. &
                .not. allocated(wdsrf) .or. .not. allocated(wetwat) .or. &
                .not. allocated(scv)) CALL CoLM_stop ('incomplete material restart water state')
            IF (size(ldew_rain) /= numpatch .or. size(ldew_snow) /= numpatch .or. &
                size(wa) /= numpatch .or. size(wdsrf) /= numpatch .or. &
                size(wetwat) /= numpatch .or. size(scv) /= numpatch .or. &
                size(wliq_soisno,1) /= nl_soil-maxsnl .or. size(wliq_soisno,2) /= numpatch .or. &
                size(wice_soisno,1) /= nl_soil-maxsnl .or. size(wice_soisno,2) /= numpatch) &
                CALL CoLM_stop ('material restart water shape mismatch')
            IF (allocated(waterstorage)) THEN
               IF (size(waterstorage) /= numpatch) CALL CoLM_stop ('material restart waterstorage shape mismatch')
               CALL write_mt_restart_all(file_restart, maxsnl, nl_soil, numpatch, &
                  ldew_rain, ldew_snow, wliq_soisno, wice_soisno, wa, wdsrf, wetwat, scv, &
                  compress, waterstorage)
            ELSE
               CALL write_mt_restart_all(file_restart, maxsnl, nl_soil, numpatch, &
                  ldew_rain, ldew_snow, wliq_soisno, wice_soisno, wa, wdsrf, wetwat, scv, compress)
            ENDIF
         ELSE
            IF (p_is_worker) THEN
               IF (numpatch /= 0) CALL CoLM_stop ('missing material restart patch state on worker')
            ENDIF
            CALL write_mt_restart_all(file_restart, maxsnl, nl_soil, 0, &
               empty_patch, empty_patch, empty_soilsnow, empty_soilsnow, &
               empty_patch, empty_patch, empty_patch, empty_patch, compress)
         ENDIF
   END SUBROUTINE mt_write_host_restart


   integer FUNCTION land_transport_mt_count ()
      IMPLICIT NONE
      integer :: imt

      land_transport_mt_count = 0
      DO imt = 1, nmts
         IF (mt_uses_generic_water(imt)) &
            land_transport_mt_count = land_transport_mt_count + 1
      ENDDO
   END FUNCTION land_transport_mt_count

   !-------------------------------------------------------------------
   ! Verify the fastest-varying extent of a restart variable matches the
   ! currently configured generic-transport row count. The land restart
   ! writer stores each pool as a compact (transport, [soilsnow,] patch)
   ! variable, and
   ! ncio_read_vector → ncio_read_serial auto-reallocates the per-block
   ! scatter buffer to the on-disk shape. If the file was written with a
   ! different transport descriptor, the downstream mpi_scatterv still uses
   ! the current compact row count as its element multiplier and will misalign or
   ! read past the buffer. Detect the mismatch up front so the caller
   ! can fall back to mt_init_from_water instead of crashing.
   !-------------------------------------------------------------------
   logical FUNCTION material_dim_matches (file_restart, varname, expect_soilsnow, patch_only)
#ifdef USEMPI
      USE MOD_SPMD_Task, only: p_is_io, p_comm_glb, p_err, &
         MPI_IN_PLACE, MPI_INTEGER, MPI_SUM
#else
      USE MOD_SPMD_Task, only: p_is_io
#endif
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart, varname
      ! Optional second-dimension size (soilsnow extent = nl_soil-maxsnl).
      ! When present, material_dim_matches additionally requires
      ! varsize(2) == expect_soilsnow; this guards against a restart
      ! written with different maxsnl / nl_soil values, which would
      ! otherwise pass the row-count-only check and then silently misalign
      ! per-layer data via ncio_read_vector's reshape.
      integer, intent(in), optional :: expect_soilsnow
      logical, intent(in), optional :: patch_only
      integer, allocatable :: varsize(:)
      integer :: iblkgrp, iblk, jblk, expected_rank, ntransport
      integer :: counts(3)
      character(len=256) :: fileblock
      logical :: block_shape_ok

      counts(:) = 0
      expected_rank = merge(3, 2, present(expect_soilsnow))
      IF (present(patch_only)) THEN
         IF (patch_only) expected_rank = 1
      ENDIF
      ntransport = land_transport_mt_count()

      ! Vector restart files are split by block via get_filename_block().
      ! Checking the unsuffixed base filename makes hot starts look like
      ! old/missing material restarts and silently reinitializes NSS state.
      IF (p_is_io) THEN
         counts(1) = landpatch%nblkgrp
         DO iblkgrp = 1, landpatch%nblkgrp
            iblk = landpatch%xblkgrp(iblkgrp)
            jblk = landpatch%yblkgrp(iblkgrp)
            CALL get_filename_block(file_restart, iblk, jblk, fileblock)

            IF (.not. ncio_var_exist(fileblock, varname, readflag = .false.)) THEN
               CYCLE
            ENDIF
            counts(2) = counts(2) + 1

            CALL ncio_inquire_varsize(fileblock, varname, varsize)
            block_shape_ok = .false.
            IF (allocated(varsize)) THEN
               IF (size(varsize) == expected_rank) THEN
                  block_shape_ok = varsize(expected_rank) == landpatch%vecgs%vlen(iblk,jblk)
                  IF (expected_rank > 1) block_shape_ok = block_shape_ok .and. varsize(1) == ntransport
                  IF (block_shape_ok .and. present(expect_soilsnow)) THEN
                     block_shape_ok = varsize(2) == expect_soilsnow
                  ENDIF
               ENDIF
            ENDIF
            IF (allocated(varsize)) deallocate(varsize)
            IF (.not. block_shape_ok) counts(3) = counts(3) + 1
         ENDDO
      ENDIF

#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, counts, 3, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif

      material_dim_matches = ntransport > 0 .and. counts(1) > 0 .and. &
         counts(2) == counts(1) .and. counts(3) == 0
   END FUNCTION material_dim_matches

   logical FUNCTION land_mt_descriptor_matches (file_restart, state_present, restart_schema)
#ifdef USEMPI
      USE MOD_SPMD_Task, only: p_is_io, p_iam_glb, p_comm_glb, p_err, &
         MPI_IN_PLACE, MPI_INTEGER, MPI_SUM, MPI_MIN
#else
      USE MOD_SPMD_Task, only: p_is_io
#endif
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      logical, intent(out), optional :: state_present
      integer, intent(out), optional :: restart_schema
      integer :: iblkgrp, iblk, jblk, schema, complete, transport_count, imt
      integer :: counts(7)
      integer, allocatable :: expected(:,:), ondisk(:,:)
      integer, allocatable :: local_reference_identity(:,:), reference_identity(:,:)
      integer, allocatable :: varsize(:)
      integer :: local_reference_schema, local_reference_count
      integer :: local_reference_rank, reference_rank(1), reference_metadata(2)
      integer :: generation_mismatch_count(1)
      real(r8) :: ondisk_mixing
      character(len=256) :: fileblock
      logical :: block_ok, has_schema, has_count, has_identity, has_commit, has_mixing
      logical :: has_local_generation, generation_mismatch, require_isotope_reference

      counts = 0
      has_local_generation = .false.
      generation_mismatch = .false.
      local_reference_schema = -1
      local_reference_count = -1
      reference_metadata = -1
      require_isotope_reference = .false.
      IF (DEF_USE_VariablySaturatedFlow) THEN
         DO imt = 1, nmts
            IF (mt_uses_generic_water(imt) .and. material_is_isotope(imt)) &
               require_isotope_reference = .true.
         ENDDO
      ENDIF
      CALL material_build_descriptor_identity(expected, transport_only=.true.)

      IF (p_is_io) THEN
         DO iblkgrp = 1, landpatch%nblkgrp
            iblk = landpatch%xblkgrp(iblkgrp)
            jblk = landpatch%yblkgrp(iblkgrp)
            CALL get_filename_block(file_restart, iblk, jblk, fileblock)
            counts(1) = counts(1) + 1
            has_schema = ncio_var_exist(fileblock, 'mt_land_restart_schema', readflag=.false.)
            has_count = ncio_var_exist(fileblock, 'mt_land_transport_count', readflag=.false.)
            has_identity = ncio_var_exist(fileblock, 'mt_land_descriptor_identity', readflag=.false.)
            has_commit = ncio_var_exist(fileblock, 'mt_land_restart_complete', readflag=.false.)
            has_mixing = ncio_var_exist(fileblock, 'mt_aquifer_mixing_water_mm', readflag=.false.)

            ! No marker is an old, non-transactional restart.  It is safe only
            ! when every block is old: the caller cold-starts the entire generic
            ! land domain before any vector scatter.
            IF (.not. has_commit) THEN
               ! A pre-transactional file is an old generic material restart
               ! only if it actually carries a generic prognostic pool.
               IF (ncio_var_exist(fileblock, 'mt_wa', readflag=.false.)) &
                  counts(6) = counts(6) + 1
               counts(2) = counts(2) + 1
               CYCLE
            ENDIF
            IF (.not. has_schema) THEN
               counts(5) = counts(5) + 1
               CYCLE
            ENDIF

            ! Validate ranks before dispatching to rank-specific readers.  Once
            ! a commit marker exists, partial/malformed metadata is corruption,
            ! not a legacy format that may be silently defaulted.
            block_ok = .true.
            schema = -1
            CALL ncio_inquire_varsize(fileblock, 'mt_land_restart_complete', varsize)
            IF (.not. allocated(varsize)) THEN
               block_ok = .false.
            ELSE
               block_ok = size(varsize) == 0
               deallocate(varsize)
            ENDIF
            CALL ncio_inquire_varsize(fileblock, 'mt_land_restart_schema', varsize)
            IF (.not. allocated(varsize)) THEN
               block_ok = .false.
            ELSE
               block_ok = size(varsize) == 0
               deallocate(varsize)
            ENDIF
            IF (block_ok) THEN
               CALL ncio_read_serial(fileblock, 'mt_land_restart_complete', complete)
               IF (complete /= 1) block_ok = .false.
            ENDIF
            IF (block_ok) THEN
               CALL ncio_read_serial(fileblock, 'mt_land_restart_schema', schema)
            ENDIF
            transport_count = -1
            IF (block_ok .and. &
                (schema == LAND_MT_RESTART_SCHEMA_VERSION .or. schema == 4)) THEN
               IF (.not. has_count) THEN
                  block_ok = .false.
               ELSE
                  CALL ncio_inquire_varsize(fileblock, 'mt_land_transport_count', varsize)
                  IF (.not. allocated(varsize)) THEN
                     block_ok = .false.
                  ELSE
                     block_ok = size(varsize) == 0
                     deallocate(varsize)
                  ENDIF
                  IF (block_ok) THEN
                     CALL ncio_read_serial(fileblock, 'mt_land_transport_count', transport_count)
                     block_ok = transport_count >= 0 .and. transport_count <= 1000
                  ENDIF
               ENDIF
            ENDIF
            IF (block_ok .and. schema == LAND_MT_RESTART_SCHEMA_VERSION) THEN
               IF (.not. has_mixing) THEN
                  block_ok = .false.
               ELSE
                  CALL ncio_inquire_varsize(fileblock, 'mt_aquifer_mixing_water_mm', varsize)
                  IF (.not. allocated(varsize)) THEN
                     block_ok = .false.
                  ELSE
                     block_ok = size(varsize) == 0
                     deallocate(varsize)
                  ENDIF
                  IF (block_ok) THEN
                     CALL ncio_read_serial(fileblock, 'mt_aquifer_mixing_water_mm', ondisk_mixing)
                     ! No silent reinterpretation of the persistent isotope inventory.
                     IF (require_isotope_reference .and. transport_count > 0) THEN
                        block_ok = ondisk_mixing == DEF_MT_AQUIFER_MIXING_WATER_MM
                        IF (.not. block_ok) counts(7) = counts(7) + 1
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF

            ! Schemas 4/5 represent a provider-only checkpoint explicitly with
            ! transport_count=0.  It intentionally carries no zero-length
            ! descriptor dimension; a stale descriptor in a reused file is
            ! ignored because the committed count is the source of truth.
            IF (block_ok .and. &
                (schema == LAND_MT_RESTART_SCHEMA_VERSION .or. schema == 4) .and. &
                transport_count > 0) THEN
               IF (.not. has_identity) THEN
                  block_ok = .false.
               ELSE
                  CALL ncio_inquire_varsize(fileblock, 'mt_land_descriptor_identity', varsize)
                  IF (.not. allocated(varsize)) THEN
                     block_ok = .false.
                  ELSE
                     IF (size(varsize) /= 2) THEN
                        block_ok = .false.
                     ELSE
                        block_ok = varsize(1) == LAND_MT_DESCRIPTOR_FIELDS .and. &
                           varsize(2) == transport_count
                     ENDIF
                     deallocate(varsize)
                  ENDIF
                  IF (block_ok) &
                     CALL ncio_read_serial(fileblock, 'mt_land_descriptor_identity', ondisk)
               ENDIF
            ELSEIF (block_ok .and. schema /= LAND_MT_RESTART_SCHEMA_VERSION .and. &
                    schema /= 4) THEN
               ! The previous transactional schema had no explicit count.  A
               ! committed descriptor is coherent but incompatible and must
               ! cold-start, while marker/schema-only fragments remain fatal.
               IF (.not. has_identity) THEN
                  block_ok = .false.
               ELSE
                  CALL ncio_inquire_varsize(fileblock, 'mt_land_descriptor_identity', varsize)
                  IF (.not. allocated(varsize)) THEN
                     block_ok = .false.
                  ELSE
                     block_ok = size(varsize) == 2
                     IF (block_ok) THEN
                        block_ok = varsize(1) == LAND_MT_DESCRIPTOR_FIELDS .and. &
                           varsize(2) > 0 .and. varsize(2) <= 1000
                        IF (block_ok) transport_count = varsize(2)
                     ENDIF
                     deallocate(varsize)
                  ENDIF
                  IF (block_ok) &
                     CALL ncio_read_serial(fileblock, 'mt_land_descriptor_identity', ondisk)
               ENDIF
            ENDIF

            ! Compare every individually coherent block with the first one
            ! owned by this IO rank.  The cross-rank reference below then
            ! compares those exact local representatives.  Count zero carries
            ! no identity by design, so stale identity bytes are never read or
            ! admitted into the generation key.
            IF (block_ok) THEN
               IF (transport_count > 0 .and. .not. allocated(ondisk)) block_ok = .false.
            ENDIF
            ! A committed count=0 is a provider-only/hydrology checkpoint,
            ! not old isotope state to reinterpret. The count dominates stale
            ! vector variables in a reused NetCDF target by design.
            IF (block_ok .and. transport_count > 0) counts(6) = counts(6) + 1
            IF (block_ok) THEN
               IF (.not. has_local_generation) THEN
                  has_local_generation = .true.
                  local_reference_schema = schema
                  local_reference_count = transport_count
                  IF (transport_count > 0) THEN
                     allocate(local_reference_identity(size(ondisk, 1), size(ondisk, 2)))
                     local_reference_identity = ondisk
                  ENDIF
               ELSEIF (schema /= local_reference_schema .or. &
                       transport_count /= local_reference_count) THEN
                  generation_mismatch = .true.
               ELSEIF (transport_count > 0) THEN
                  IF (.not. allocated(local_reference_identity)) THEN
                     generation_mismatch = .true.
                  ELSEIF (any(ondisk /= local_reference_identity)) THEN
                     generation_mismatch = .true.
                  ENDIF
               ENDIF
            ENDIF
            IF (.not. block_ok) THEN
               counts(5) = counts(5) + 1
            ELSEIF ((schema /= LAND_MT_RESTART_SCHEMA_VERSION .and. &
                     (schema /= 4 .or. require_isotope_reference)) .or. &
                    transport_count /= size(expected, 2)) THEN
               counts(4) = counts(4) + 1
            ELSEIF (transport_count == 0) THEN
               counts(3) = counts(3) + 1
            ELSEIF (&
                    size(ondisk, 1) /= size(expected, 1) .or. &
                    size(ondisk, 2) /= size(expected, 2)) THEN
               counts(4) = counts(4) + 1
            ELSEIF (.not. all(ondisk == expected)) THEN
               counts(4) = counts(4) + 1
            ELSE
               counts(3) = counts(3) + 1
            ENDIF
            IF (allocated(ondisk)) deallocate(ondisk)
         ENDDO
      ENDIF

#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, counts, 7, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)

      ! Select one rank's exact generation as the global reference.  A single
      ! broadcast plus a logical reduction is exact (unlike a checksum) while
      ! avoiding a maximum-sized identity reduction on every rank.
      local_reference_rank = huge(0)
      IF (has_local_generation) local_reference_rank = p_iam_glb
      reference_rank(1) = local_reference_rank
      CALL mpi_allreduce(MPI_IN_PLACE, reference_rank, 1, MPI_INTEGER, MPI_MIN, p_comm_glb, p_err)
      IF (reference_rank(1) /= huge(0)) THEN
         reference_metadata = 0
         IF (p_iam_glb == reference_rank(1)) THEN
            reference_metadata(1) = local_reference_schema
            reference_metadata(2) = local_reference_count
         ENDIF
         CALL mpi_bcast(reference_metadata, 2, MPI_INTEGER, reference_rank(1), p_comm_glb, p_err)
         IF (has_local_generation) THEN
            IF (local_reference_schema /= reference_metadata(1) .or. &
                local_reference_count /= reference_metadata(2)) generation_mismatch = .true.
         ENDIF
         IF (reference_metadata(2) > 0) THEN
            allocate(reference_identity(LAND_MT_DESCRIPTOR_FIELDS, reference_metadata(2)))
            reference_identity = 0
            IF (p_iam_glb == reference_rank(1)) THEN
               IF (allocated(local_reference_identity)) THEN
                  reference_identity = local_reference_identity
               ELSE
                  generation_mismatch = .true.
               ENDIF
            ENDIF
            CALL mpi_bcast(reference_identity, size(reference_identity), MPI_INTEGER, &
               reference_rank(1), p_comm_glb, p_err)
            IF (has_local_generation .and. &
                local_reference_schema == reference_metadata(1) .and. &
                local_reference_count == reference_metadata(2)) THEN
               IF (.not. allocated(local_reference_identity)) THEN
                  generation_mismatch = .true.
               ELSEIF (any(local_reference_identity /= reference_identity)) THEN
                  generation_mismatch = .true.
               ENDIF
            ENDIF
            deallocate(reference_identity)
         ENDIF
      ENDIF
      generation_mismatch_count(1) = merge(1, 0, generation_mismatch)
      CALL mpi_allreduce(MPI_IN_PLACE, generation_mismatch_count, 1, MPI_INTEGER, &
         MPI_SUM, p_comm_glb, p_err)
      generation_mismatch = generation_mismatch_count(1) > 0
#endif

      IF (counts(7) > 0) &
         CALL CoLM_stop('aquifer reference calibration differs from material restart')
      IF (generation_mismatch .or. counts(5) > 0 .or. &
          (counts(2) > 0 .and. counts(2) < counts(1)) .or. &
          (counts(3) > 0 .and. counts(4) > 0)) THEN
         CALL CoLM_stop('malformed/incomplete generic land material restart transaction')
      ENDIF

      ! A complete but different schema/descriptor is a coherent checkpoint
      ! for another configuration.  Cold-start the generic domain as one unit;
      ! never mix its rows with the current descriptor.
      land_mt_descriptor_matches = counts(1) > 0 .and. counts(3) == counts(1)
      IF (present(state_present)) state_present = counts(6) > 0
      IF (present(restart_schema)) THEN
#ifdef USEMPI
         restart_schema = reference_metadata(1)
#else
         restart_schema = local_reference_schema
#endif
      ENDIF
      IF (allocated(local_reference_identity)) deallocate(local_reference_identity)
      deallocate(expected)
   END FUNCTION land_mt_descriptor_matches

   SUBROUTINE write_land_mt_transaction_marker (file_restart, complete)
      USE MOD_SPMD_Task, only: p_is_io
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer, intent(in) :: complete
      integer :: iblkgrp, iblk, jblk
      character(len=256) :: fileblock

      IF (p_is_io) THEN
         DO iblkgrp = 1, landpatch%nblkgrp
            iblk = landpatch%xblkgrp(iblkgrp)
            jblk = landpatch%yblkgrp(iblkgrp)
            CALL get_filename_block(file_restart, iblk, jblk, fileblock)
            CALL ncio_write_serial(fileblock, 'mt_land_restart_complete', complete)
         ENDDO
      ENDIF
   END SUBROUTINE write_land_mt_transaction_marker

   SUBROUTINE write_land_mt_descriptor_metadata (file_restart)
      USE MOD_SPMD_Task, only: p_is_io
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer :: iblkgrp, iblk, jblk
      integer, allocatable :: identity(:,:)
      character(len=256) :: fileblock

      CALL material_build_descriptor_identity(identity, transport_only=.true.)
      IF (p_is_io) THEN
         DO iblkgrp = 1, landpatch%nblkgrp
            iblk = landpatch%xblkgrp(iblkgrp)
            jblk = landpatch%yblkgrp(iblkgrp)
            CALL get_filename_block(file_restart, iblk, jblk, fileblock)
            CALL ncio_write_serial(fileblock, 'mt_land_restart_schema', LAND_MT_RESTART_SCHEMA_VERSION)
            CALL ncio_write_serial(fileblock, 'mt_land_transport_count', size(identity, 2))
            CALL ncio_write_serial(fileblock, 'mt_aquifer_mixing_water_mm', &
               DEF_MT_AQUIFER_MIXING_WATER_MM)
            IF (size(identity, 2) > 0) THEN
               CALL ncio_define_dimension(fileblock, 'mt_land_descriptor_field', LAND_MT_DESCRIPTOR_FIELDS)
               CALL ncio_define_dimension(fileblock, 'mt_land_transport', size(identity, 2))
               CALL ncio_write_serial(fileblock, 'mt_land_descriptor_identity', identity, &
                  'mt_land_descriptor_field', 'mt_land_transport', DEF_REST_CompressLevel)
            ENDIF
         ENDDO
      ENDIF
      deallocate(identity)
   END SUBROUTINE write_land_mt_descriptor_metadata

   SUBROUTINE pack_transport_patch (source, packed)
      IMPLICIT NONE
      real(r8), intent(in) :: source(:,:)
      real(r8), intent(out) :: packed(:,:)
      integer :: imt, k

      packed = 0._r8
      k = 0
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         k = k + 1
         packed(k, :) = source(imt, :)
      ENDDO
   END SUBROUTINE pack_transport_patch

   SUBROUTINE pack_transport_soilsnow (source, packed)
      IMPLICIT NONE
      real(r8), intent(in) :: source(:,:,:)
      real(r8), intent(out) :: packed(:,:,:)
      integer :: imt, k

      packed = 0._r8
      k = 0
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         k = k + 1
         packed(k, :, :) = source(imt, :, :)
      ENDDO
   END SUBROUTINE pack_transport_soilsnow

   SUBROUTINE read_transport_patch_field (file_restart, varname, target)
      USE MOD_SPMD_Task, only: p_is_worker
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart, varname
      real(r8), allocatable, intent(inout) :: target(:,:)
      real(r8), allocatable :: packed(:,:)
      integer :: imt, k, ntransport

      ntransport = land_transport_mt_count()
      CALL ncio_read_vector(file_restart, varname, ntransport, landpatch, packed)
      IF (p_is_worker) THEN
         IF (landpatch%nset > 0 .and. .not. allocated(target)) &
            CALL CoLM_stop('generic land material restart target is not allocated')
         IF (allocated(target)) THEN
            target = 0._r8
            k = 0
            DO imt = 1, nmts
               IF (.not. mt_uses_generic_water(imt)) CYCLE
               k = k + 1
               target(imt, :) = packed(k, :)
            ENDDO
         ENDIF
      ENDIF
      IF (allocated(packed)) deallocate(packed)
   END SUBROUTINE read_transport_patch_field

   SUBROUTINE read_transport_soilsnow_field (file_restart, varname, nlayers, target)
      USE MOD_SPMD_Task, only: p_is_worker
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart, varname
      integer, intent(in) :: nlayers
      real(r8), allocatable, intent(inout) :: target(:,:,:)
      real(r8), allocatable :: packed(:,:,:)
      integer :: imt, k, ntransport

      ntransport = land_transport_mt_count()
      CALL ncio_read_vector(file_restart, varname, ntransport, nlayers, landpatch, packed)
      IF (p_is_worker) THEN
         IF (landpatch%nset > 0 .and. .not. allocated(target)) &
            CALL CoLM_stop('generic layered land material restart target is not allocated')
         IF (allocated(target)) THEN
            target = 0._r8
            k = 0
            DO imt = 1, nmts
               IF (.not. mt_uses_generic_water(imt)) CYCLE
               k = k + 1
               target(imt, :, :) = packed(k, :, :)
            ENDDO
         ENDIF
      ENDIF
      IF (allocated(packed)) deallocate(packed)
   END SUBROUTINE read_transport_soilsnow_field

   SUBROUTINE count_nonnegative_patch (field, invalid_counts)
      USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
      IMPLICIT NONE
      real(r8), intent(in) :: field(:,:)
      integer, intent(inout) :: invalid_counts(:)
      integer :: imt, ip

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         DO ip = 1, size(field, 2)
            IF (.not. ieee_is_finite(field(imt, ip))) THEN
               invalid_counts(1) = invalid_counts(1) + 1
            ELSEIF (field(imt, ip) < -LAND_MT_RESTART_NEGATIVE_DUST) THEN
               invalid_counts(2) = invalid_counts(2) + 1
            ENDIF
         ENDDO
      ENDDO
   END SUBROUTINE count_nonnegative_patch

   SUBROUTINE count_nonnegative_soilsnow (field, invalid_counts)
      USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
      IMPLICIT NONE
      real(r8), intent(in) :: field(:,:,:)
      integer, intent(inout) :: invalid_counts(:)
      integer :: imt, j, ip

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         DO ip = 1, size(field, 3)
            DO j = 1, size(field, 2)
               IF (.not. ieee_is_finite(field(imt, j, ip))) THEN
                  invalid_counts(1) = invalid_counts(1) + 1
               ELSEIF (field(imt, j, ip) < -LAND_MT_RESTART_NEGATIVE_DUST) THEN
                  invalid_counts(2) = invalid_counts(2) + 1
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   END SUBROUTINE count_nonnegative_soilsnow

   SUBROUTINE count_signed_patch (field, invalid_counts)
      USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
      IMPLICIT NONE
      real(r8), intent(in) :: field(:,:)
      integer, intent(inout) :: invalid_counts(:)
      integer :: imt, ip

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         DO ip = 1, size(field, 2)
            IF (.not. ieee_is_finite(field(imt, ip))) &
               invalid_counts(3) = invalid_counts(3) + 1
         ENDDO
      ENDDO
   END SUBROUTINE count_signed_patch

   SUBROUTINE count_peclet_patch (field, invalid_counts)
      USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
      IMPLICIT NONE
      real(r8), intent(in) :: field(:,:)
      integer, intent(inout) :: invalid_counts(:)
      integer :: imt, ip

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         DO ip = 1, size(field, 2)
            IF (.not. ieee_is_finite(field(imt, ip))) THEN
               invalid_counts(3) = invalid_counts(3) + 1
            ELSEIF (field(imt, ip) < 0._r8 .or. field(imt, ip) > 1._r8) THEN
               invalid_counts(4) = invalid_counts(4) + 1
            ENDIF
         ENDDO
      ENDDO
   END SUBROUTINE count_peclet_patch

   SUBROUTINE clamp_nonnegative_patch (field)
      IMPLICIT NONE
      real(r8), intent(inout) :: field(:,:)
      integer :: imt

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         WHERE (field(imt, :) < 0._r8) field(imt, :) = 0._r8
      ENDDO
   END SUBROUTINE clamp_nonnegative_patch

   SUBROUTINE clamp_nonnegative_soilsnow (field)
      IMPLICIT NONE
      real(r8), intent(inout) :: field(:,:,:)
      integer :: imt

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         WHERE (field(imt, :, :) < 0._r8) field(imt, :, :) = 0._r8
      ENDDO
   END SUBROUTINE clamp_nonnegative_soilsnow

   SUBROUTINE validate_land_mt_restart_state (wa)
      USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
#ifdef USEMPI
      USE MOD_SPMD_Task, only: p_is_worker, p_is_master, p_comm_glb, p_err, &
         MPI_IN_PLACE, MPI_INTEGER, MPI_SUM
#else
      USE MOD_SPMD_Task, only: p_is_worker, p_is_master
#endif
      IMPLICIT NONE
      real(r8), intent(in), optional :: wa(:)
      integer :: invalid_counts(7), imt, ip
      logical :: require_isotope_reference
      real(r8) :: expected_refmass

      ! [nonfinite physical amount, significantly negative physical amount,
      !  nonfinite signed/diagnostic value, Peclet outside [0,1],
      !  missing state allocation, invalid isotope aquifer state, reference]. mt_wa is
      ! deliberately signed because the wetland correction debt is prognostic;
      ! leaf isotope storage is likewise a signed NSS anomaly.
      invalid_counts = 0
      require_isotope_reference = .false.
      IF (DEF_USE_VariablySaturatedFlow) THEN
         DO imt = 1, nmts
            IF (mt_uses_generic_water(imt) .and. material_is_isotope(imt)) &
               require_isotope_reference = .true.
         ENDDO
      ENDIF
      IF (p_is_worker) THEN
#ifdef MATERIAL_SIGNED_SOIL_DEBT_EXPERIMENT
         IF (allocated(mt_wliq_soisno)) THEN
            IF (any(mt_wliq_soisno < 0._r8)) &
               CALL CoLM_stop('signed soil debt: rejected checkpoint; direct-kernel experiment only')
         ENDIF
#endif
         IF (landpatch%nset > 0) THEN
            IF (.not. (allocated(mt_ldew_rain) .and. allocated(mt_ldew_snow) .and. &
                       allocated(mt_wliq_soisno) .and. allocated(mt_wice_soisno) .and. &
                       allocated(mt_wa) .and. allocated(mt_aquifer_ref_water) .and. &
                       allocated(mt_aquifer_ref_mass) .and. &
                       allocated(mt_wdsrf) .and. allocated(mt_wetwat) .and. &
                       allocated(mt_surface_residue) .and. allocated(mt_subsurface_residue) .and. &
                       allocated(mt_solid_soisno) .and. allocated(mt_canopy_solid) .and. &
                       allocated(mt_surface_solid) .and. allocated(mt_subsurface_solid) .and. &
                       allocated(mt_waterstorage_solid) .and. allocated(mt_scv) .and. &
                       allocated(mt_waterstorage) .and. allocated(mt_leaf_delta_e) .and. &
                       allocated(mt_leaf_delta_b) .and. allocated(mt_leaf_peclet) .and. &
                       allocated(mt_leaf_water_moles) .and. allocated(mt_leaf_iso_storage))) &
               invalid_counts(5) = invalid_counts(5) + 1
         ENDIF
         IF (allocated(mt_ldew_rain)) CALL count_nonnegative_patch(mt_ldew_rain, invalid_counts)
         IF (allocated(mt_ldew_snow)) CALL count_nonnegative_patch(mt_ldew_snow, invalid_counts)
         IF (allocated(mt_wliq_soisno)) CALL count_nonnegative_soilsnow(mt_wliq_soisno, invalid_counts)
         IF (allocated(mt_wice_soisno)) CALL count_nonnegative_soilsnow(mt_wice_soisno, invalid_counts)
         IF (allocated(mt_wa)) CALL count_signed_patch(mt_wa, invalid_counts)
         IF (allocated(mt_aquifer_ref_mass)) CALL count_nonnegative_patch(mt_aquifer_ref_mass, invalid_counts)
         IF (allocated(mt_aquifer_ref_water)) THEN
            DO ip = 1, size(mt_aquifer_ref_water)
               IF (.not. ieee_is_finite(mt_aquifer_ref_water(ip))) THEN
                  invalid_counts(7) = invalid_counts(7) + 1
               ELSEIF (mt_aquifer_ref_water(ip) < 0._r8) THEN
                  invalid_counts(7) = invalid_counts(7) + 1
               ELSEIF (.not. require_isotope_reference .and. mt_aquifer_ref_water(ip) > 0._r8) THEN
                  invalid_counts(7) = invalid_counts(7) + 1
               ELSEIF (allocated(patchtype)) THEN
                  IF (ip <= size(patchtype)) THEN
                     IF (require_isotope_reference .and. &
                         (patchtype(ip) == 0 .or. patchtype(ip) == 2) .and. &
                         mt_aquifer_ref_water(ip) <= 0._r8) invalid_counts(7) = invalid_counts(7) + 1
                     IF (patchtype(ip) /= 0 .and. patchtype(ip) /= 2 .and. &
                         mt_aquifer_ref_water(ip) > 0._r8) invalid_counts(7) = invalid_counts(7) + 1
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
         IF (allocated(mt_aquifer_ref_water) .and. allocated(mt_aquifer_ref_mass)) THEN
            DO imt = 1, nmts
               IF (.not. mt_uses_generic_water(imt)) CYCLE
               DO ip = 1, min(size(mt_aquifer_ref_water), size(mt_aquifer_ref_mass, 2))
                  expected_refmass = 0._r8
                  IF (material_is_isotope(imt)) &
                     expected_refmass = mt_aquifer_ref_water(ip) * materials(imt)%ref_ratio
                  IF (.not. ieee_is_finite(expected_refmass)) THEN
                     invalid_counts(7) = invalid_counts(7) + 1
                  ELSEIF (abs(mt_aquifer_ref_mass(imt, ip) - expected_refmass) > &
                          max(1.e-12_r8, 1.e-10_r8 * abs(expected_refmass))) THEN
                     invalid_counts(7) = invalid_counts(7) + 1
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
         IF (present(wa)) THEN
            IF (allocated(mt_wa) .and. allocated(mt_aquifer_ref_water) .and. &
                allocated(mt_aquifer_ref_mass)) THEN
               DO imt = 1, nmts
                  IF (.not. material_is_isotope(imt)) CYCLE
                  IF (.not. mt_uses_generic_water(imt)) CYCLE
                  DO ip = 1, min(size(wa), size(mt_wa, 2))
                     IF (.not. material_aquifer_isotope_state_valid(wa(ip), mt_wa(imt, ip), &
                           materials(imt)%ref_ratio, mt_aquifer_ref_water(ip), &
                           mt_aquifer_ref_mass(imt, ip))) &
                        invalid_counts(6) = invalid_counts(6) + 1
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
         IF (allocated(mt_wdsrf)) CALL count_nonnegative_patch(mt_wdsrf, invalid_counts)
         IF (allocated(mt_wetwat)) CALL count_nonnegative_patch(mt_wetwat, invalid_counts)
         IF (allocated(mt_surface_residue)) CALL count_nonnegative_patch(mt_surface_residue, invalid_counts)
         IF (allocated(mt_subsurface_residue)) CALL count_nonnegative_patch(mt_subsurface_residue, invalid_counts)
         IF (allocated(mt_solid_soisno)) CALL count_nonnegative_soilsnow(mt_solid_soisno, invalid_counts)
         IF (allocated(mt_canopy_solid)) CALL count_nonnegative_patch(mt_canopy_solid, invalid_counts)
         IF (allocated(mt_surface_solid)) CALL count_nonnegative_patch(mt_surface_solid, invalid_counts)
         IF (allocated(mt_subsurface_solid)) CALL count_nonnegative_patch(mt_subsurface_solid, invalid_counts)
         IF (allocated(mt_waterstorage_solid)) CALL count_nonnegative_patch(mt_waterstorage_solid, invalid_counts)
         IF (allocated(mt_scv)) CALL count_nonnegative_patch(mt_scv, invalid_counts)
         IF (allocated(mt_waterstorage)) CALL count_nonnegative_patch(mt_waterstorage, invalid_counts)
         IF (allocated(mt_leaf_delta_e)) CALL count_signed_patch(mt_leaf_delta_e, invalid_counts)
         IF (allocated(mt_leaf_delta_b)) CALL count_signed_patch(mt_leaf_delta_b, invalid_counts)
         IF (allocated(mt_leaf_peclet)) CALL count_peclet_patch(mt_leaf_peclet, invalid_counts)
         IF (allocated(mt_leaf_water_moles)) CALL count_nonnegative_patch(mt_leaf_water_moles, invalid_counts)
         IF (allocated(mt_leaf_iso_storage)) CALL count_signed_patch(mt_leaf_iso_storage, invalid_counts)
      ENDIF
#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, invalid_counts, size(invalid_counts), &
         MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif

      IF (any(invalid_counts > 0)) THEN
         IF (p_is_master) WRITE(*,'(A,7(I0,1X))') &
            'ERROR generic land material restart state counts '// &
            '[amount_nan amount_neg signed_nan peclet_range missing_state aquifer_state reference]: ', invalid_counts
         CALL CoLM_stop('invalid generic land material restart state')
      ENDIF

      ! Only canonicalize the tiny negative dust admitted above.  Signed debt
      ! and anomaly fields are never clipped.
      IF (p_is_worker) THEN
         IF (allocated(mt_ldew_rain)) CALL clamp_nonnegative_patch(mt_ldew_rain)
         IF (allocated(mt_ldew_snow)) CALL clamp_nonnegative_patch(mt_ldew_snow)
         IF (allocated(mt_wliq_soisno)) CALL clamp_nonnegative_soilsnow(mt_wliq_soisno)
         IF (allocated(mt_wice_soisno)) CALL clamp_nonnegative_soilsnow(mt_wice_soisno)
         IF (allocated(mt_wdsrf)) CALL clamp_nonnegative_patch(mt_wdsrf)
         IF (allocated(mt_wetwat)) CALL clamp_nonnegative_patch(mt_wetwat)
         IF (allocated(mt_surface_residue)) CALL clamp_nonnegative_patch(mt_surface_residue)
         IF (allocated(mt_subsurface_residue)) CALL clamp_nonnegative_patch(mt_subsurface_residue)
         IF (allocated(mt_solid_soisno)) CALL clamp_nonnegative_soilsnow(mt_solid_soisno)
         IF (allocated(mt_canopy_solid)) CALL clamp_nonnegative_patch(mt_canopy_solid)
         IF (allocated(mt_surface_solid)) CALL clamp_nonnegative_patch(mt_surface_solid)
         IF (allocated(mt_subsurface_solid)) CALL clamp_nonnegative_patch(mt_subsurface_solid)
         IF (allocated(mt_waterstorage_solid)) CALL clamp_nonnegative_patch(mt_waterstorage_solid)
         IF (allocated(mt_scv)) CALL clamp_nonnegative_patch(mt_scv)
         IF (allocated(mt_waterstorage)) CALL clamp_nonnegative_patch(mt_waterstorage)
         IF (allocated(mt_leaf_water_moles)) CALL clamp_nonnegative_patch(mt_leaf_water_moles)
      ENDIF
   END SUBROUTINE validate_land_mt_restart_state

   SUBROUTINE mt_init_from_water (numpatch, maxsnl, nl_soil, &
      ldew_rain, ldew_snow, wliq_soisno, wice_soisno, &
      wa, wdsrf, wetwat, scv, waterstorage, patch_mask)

      IMPLICIT NONE
      integer,  intent(in) :: numpatch, maxsnl, nl_soil
      real(r8), intent(in) :: ldew_rain(numpatch)
      real(r8), intent(in) :: ldew_snow(numpatch)
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wa(numpatch)
      real(r8), intent(in) :: wdsrf(numpatch)
      real(r8), intent(in) :: wetwat(numpatch)
      real(r8), intent(in) :: scv(numpatch)
      ! Optional: per-patch irrigation reservoir from
      ! MOD_Vars_TimeVariables. Present when DEF_USE_IRRIGATION is on
      ! under CROP. When passed, mt_waterstorage is cold-started to
      ! waterstorage*R_init so the first mt_save_storage sees a
      ! correct starting inventory; without this the pool would be 0
      ! and the first step's balance check would under-count the
      ! reservoir mass that the water side already tracks.
      real(r8), intent(in), optional :: waterstorage(numpatch)
      ! Optional: cold-start only these patches (LULCC patches with no source
      ! patch, MOD_MT_Lulcc); the other patches and the provider-owned state
      ! are left as they are. Without it every patch is cold-started.
      logical,  intent(in), optional :: patch_mask(numpatch)

      integer  :: imt, ip, j, snl_local
      real(r8) :: R_init

      IF (present(patch_mask)) THEN
         DO ip = 1, numpatch
            IF (.not. patch_mask(ip)) CYCLE
            IF (allocated(mt_surface_residue)) mt_surface_residue(:,ip) = 0._r8
            IF (allocated(mt_subsurface_residue)) mt_subsurface_residue(:,ip) = 0._r8
            IF (allocated(mt_solid_soisno)) mt_solid_soisno(:,:,ip) = 0._r8
            IF (allocated(mt_canopy_solid)) mt_canopy_solid(:,ip) = 0._r8
            IF (allocated(mt_surface_solid)) mt_surface_solid(:,ip) = 0._r8
            IF (allocated(mt_subsurface_solid)) mt_subsurface_solid(:,ip) = 0._r8
            IF (allocated(mt_waterstorage_solid)) mt_waterstorage_solid(:,ip) = 0._r8
         ENDDO
      ELSE
         IF (allocated(mt_surface_residue)) mt_surface_residue = 0._r8
         IF (allocated(mt_subsurface_residue)) mt_subsurface_residue = 0._r8
         IF (allocated(mt_solid_soisno)) mt_solid_soisno = 0._r8
         IF (allocated(mt_canopy_solid)) mt_canopy_solid = 0._r8
         IF (allocated(mt_surface_solid)) mt_surface_solid = 0._r8
         IF (allocated(mt_subsurface_solid)) mt_subsurface_solid = 0._r8
         IF (allocated(mt_waterstorage_solid)) mt_waterstorage_solid = 0._r8
      ENDIF
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         R_init = material_init_water_ratio(imt)
         DO ip = 1, numpatch
            IF (present(patch_mask)) THEN
               IF (.not. patch_mask(ip)) CYCLE
            ENDIF
            mt_ldew_rain(imt, ip) = ldew_rain(ip) * R_init
            mt_ldew_snow(imt, ip) = ldew_snow(ip) * R_init
            CALL mt_equilibrate_dissolved(imt, max(ldew_rain(ip), 0._r8), &
               mt_ldew_rain(imt, ip), mt_canopy_solid(imt, ip))
            DO j = maxsnl + 1, nl_soil
               mt_wliq_soisno(imt, j, ip) = max(wliq_soisno(j, ip), 0._r8) * R_init
               mt_wice_soisno(imt, j, ip) = max(wice_soisno(j, ip), 0._r8) * R_init
               CALL mt_equilibrate_dissolved(imt, max(wliq_soisno(j, ip), 0._r8), &
                  mt_wliq_soisno(imt, j, ip), mt_solid_soisno(imt, j, ip))
            ENDDO
            ! wa uses the SIGNED water value so a hydrology restart with an
            ! aquifer-debt state (wa<0, recorded by the wetland branch at
            ! MOD_SoilSnowHydrology.F90:1200-1203) starts with a matching
            ! mt_wa<0. Without this the first wetland recovery step after a
            ! cold material-start would see pool_water = wa_bef + inputs but
            ! pool_material = 0 + inputs*R_atm (mt_wa_bef was clamped to 0),
            ! producing a 2× over-concentration in the recovered wetwat.
            ! wdsrf/wetwat stay with max(.,0) since WATER_VSF never leaves
            ! them negative by construction (overflow case L1196-1207).
            IF (material_is_isotope(imt) .and. DEF_USE_VariablySaturatedFlow .and. &
                DEF_MT_AQUIFER_MIXING_WATER_MM > 0._r8) THEN
               IF (.not. allocated(patchtype)) CALL CoLM_stop('missing patchtype for isotope aquifer reference')
               IF (ip > size(patchtype)) CALL CoLM_stop('missing patchtype for isotope aquifer reference')
               IF (patchtype(ip) == 0 .or. patchtype(ip) == 2) THEN
                  mt_aquifer_ref_water(ip) = DEF_MT_AQUIFER_MIXING_WATER_MM
                  mt_aquifer_ref_mass(imt, ip) = mt_aquifer_ref_water(ip) * materials(imt)%ref_ratio
                  IF (wa(ip) + mt_aquifer_ref_water(ip) < 0._r8) &
                     CALL CoLM_stop('negative isotope aquifer mixing water at cold start')
               ENDIF
               ! Initial actual composition follows INIT_DELTA, not the
               ! standard reference used only to anchor signed anomalies.
               mt_wa(imt, ip) = (wa(ip) + mt_aquifer_ref_water(ip)) * R_init &
                  - mt_aquifer_ref_mass(imt, ip)
            ELSE
               mt_wa(imt, ip) = wa(ip) * R_init
            ENDIF
            mt_wdsrf (imt, ip) = max(wdsrf(ip),  0._r8) * R_init
            mt_wetwat(imt, ip) = max(wetwat(ip), 0._r8) * R_init
            CALL mt_equilibrate_dissolved(imt, max(wdsrf(ip), 0._r8), &
               mt_wdsrf(imt, ip), mt_surface_solid(imt, ip))
            CALL mt_equilibrate_dissolved(imt, max(wetwat(ip), 0._r8), &
               mt_wetwat(imt, ip), mt_surface_solid(imt, ip))
            CALL mt_equilibrate_dissolved(imt, wa(ip), &
               mt_wa(imt, ip), mt_subsurface_solid(imt, ip))
            ! Reconstruct snl from the snow layer water content (same
            ! recipe as CoLMMAIN L770-774): snl counts non-empty snow
            ! layers from the top. snl is a runtime scalar, not a
            ! persistent per-patch array, so we infer it here.
            ! Relies on CoLM's snow-column contiguity convention: active
            ! snow layers occupy indices [snl+1 : 0] with no gaps, so
            ! iterating from j=0 downward and EXIT-ing at the first empty
            ! slot correctly recovers snl. If a future snow path violates
            ! this (gaps in the column), replace EXIT with a count of
            ! non-empty slots instead.
            snl_local = 0
            DO j = 0, maxsnl + 1, -1
               IF (wliq_soisno(j, ip) + wice_soisno(j, ip) > 0._r8) THEN
                  snl_local = snl_local - 1
               ELSE
                  EXIT
               ENDIF
            ENDDO
            ! mt_scv holds the pre-layer snow material pool. It is only
            ! populated when snl==0 (thin snow, no layer yet); once a
            ! snow layer is created, material lives in mt_wice/wliq and
            ! mt_scv stays zero.
            IF (snl_local == 0) THEN
               mt_scv(imt, ip) = max(scv(ip), 0._r8) * R_init
            ELSE
               mt_scv(imt, ip) = 0._r8
            ENDIF
            IF (present(waterstorage) .and. allocated(mt_waterstorage)) THEN
               mt_waterstorage(imt, ip) = max(waterstorage(ip), 0._r8) * R_init
               CALL mt_equilibrate_dissolved(imt, max(waterstorage(ip), 0._r8), &
                  mt_waterstorage(imt, ip), mt_waterstorage_solid(imt, ip))
            ENDIF
         ENDDO
      ENDDO
      IF (.not. present(patch_mask)) CALL zero_provider_owned_land_mt_state()
   END SUBROUTINE mt_init_from_water

   SUBROUTINE read_land_mt_restart (file_restart, maxsnl, nl_soil, found_restart, &
      scv_missing, waterstorage_missing, wa)
      USE MOD_SPMD_Task, only: p_is_master
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer, intent(in) :: maxsnl, nl_soil
      logical, intent(out) :: found_restart
      ! These outputs remain for caller ABI compatibility.  Transactional
      ! schema v5 never mixes an optional legacy pool into a hot generic load:
      ! an old checkpoint cold-starts the whole generic land domain instead.
      logical, optional, intent(out) :: scv_missing
      logical, optional, intent(out) :: waterstorage_missing
      real(r8), intent(in), optional :: wa(:)
      integer :: ntransport, imt, restart_schema
      logical :: descriptor_matches, restart_complete, field_matches, state_present, active_isotope

      found_restart = .false.
      IF (present(scv_missing)) scv_missing = .false.
      IF (present(waterstorage_missing)) waterstorage_missing = .false.
      CALL zero_provider_owned_land_mt_state()
      ntransport = land_transport_mt_count()
      descriptor_matches = land_mt_descriptor_matches(file_restart, state_present, restart_schema)
      IF (ntransport <= 0) THEN
         ! Provider-owned species have their own lifecycle restart callbacks.
         ! Validate any present transaction (marker=0/partial remains fatal),
         ! but do not create or consume zero-length generic material dimensions.
         found_restart = .true.
         RETURN
      ENDIF

      ! A missing commit marker is a legacy checkpoint and causes a coherent
      ! cold start.  Once a marker exists, the descriptor checker rejects
      ! incomplete metadata/transactions; a complete descriptor for another
      ! material configuration is safe but incompatible and also cold-starts.
      IF (.not. descriptor_matches) THEN
         active_isotope = .false.
         IF (DEF_USE_VariablySaturatedFlow) THEN
            DO imt = 1, nmts
               IF (mt_uses_generic_water(imt) .and. material_is_isotope(imt)) &
                  active_isotope = .true.
            ENDDO
         ENDIF
         IF (active_isotope .and. state_present) &
            CALL CoLM_stop('old/incompatible isotope material restart: explicit fresh material initialization required')
         IF (p_is_master) WRITE(*,'(3A)') &
            'Generic land material restart is legacy/incompatible in ', &
            TRIM(file_restart), '; using whole-domain cold start.'
         RETURN
      ENDIF

      ! A committed current transaction contains every generic prognostic and
      ! NSS field on every block with the exact compact transport-row extent.
      ! Finish all collective probes before the first vector read so no worker
      ! can observe a partially scattered checkpoint.
      restart_complete = .true.
      field_matches = material_dim_matches(file_restart, 'mt_ldew_rain')
      restart_complete = restart_complete .and. field_matches
      field_matches = material_dim_matches(file_restart, 'mt_ldew_snow')
      restart_complete = restart_complete .and. field_matches
      field_matches = material_dim_matches(file_restart, 'mt_wliq_soisno', &
                                         expect_soilsnow=nl_soil-maxsnl)
      restart_complete = restart_complete .and. field_matches
      field_matches = material_dim_matches(file_restart, 'mt_wice_soisno', &
                                         expect_soilsnow=nl_soil-maxsnl)
      restart_complete = restart_complete .and. field_matches
      field_matches = material_dim_matches(file_restart, 'mt_solid_soisno', &
                                         expect_soilsnow=nl_soil-maxsnl)
      restart_complete = restart_complete .and. field_matches
      field_matches = material_dim_matches(file_restart, 'mt_wa')
      restart_complete = restart_complete .and. field_matches
      IF (restart_schema == LAND_MT_RESTART_SCHEMA_VERSION) THEN
         field_matches = material_dim_matches(file_restart, 'mt_aquifer_ref_water', patch_only=.true.)
         restart_complete = restart_complete .and. field_matches
         field_matches = material_dim_matches(file_restart, 'mt_aquifer_ref_mass')
         restart_complete = restart_complete .and. field_matches
      ENDIF
      field_matches = material_dim_matches(file_restart, 'mt_wdsrf')
      restart_complete = restart_complete .and. field_matches
      field_matches = material_dim_matches(file_restart, 'mt_wetwat')
      restart_complete = restart_complete .and. field_matches
      field_matches = material_dim_matches(file_restart, 'mt_surface_residue')
      restart_complete = restart_complete .and. field_matches
      field_matches = material_dim_matches(file_restart, 'mt_subsurface_residue')
      restart_complete = restart_complete .and. field_matches
      field_matches = material_dim_matches(file_restart, 'mt_canopy_solid')
      restart_complete = restart_complete .and. field_matches
      field_matches = material_dim_matches(file_restart, 'mt_surface_solid')
      restart_complete = restart_complete .and. field_matches
      field_matches = material_dim_matches(file_restart, 'mt_subsurface_solid')
      restart_complete = restart_complete .and. field_matches
      field_matches = material_dim_matches(file_restart, 'mt_waterstorage_solid')
      restart_complete = restart_complete .and. field_matches
      field_matches = material_dim_matches(file_restart, 'mt_scv')
      restart_complete = restart_complete .and. field_matches
      field_matches = material_dim_matches(file_restart, 'mt_waterstorage')
      restart_complete = restart_complete .and. field_matches
      field_matches = material_dim_matches(file_restart, 'mt_leaf_delta_e')
      restart_complete = restart_complete .and. field_matches
      field_matches = material_dim_matches(file_restart, 'mt_leaf_delta_b')
      restart_complete = restart_complete .and. field_matches
      field_matches = material_dim_matches(file_restart, 'mt_leaf_peclet')
      restart_complete = restart_complete .and. field_matches
      field_matches = material_dim_matches(file_restart, 'mt_leaf_water_moles')
      restart_complete = restart_complete .and. field_matches
      field_matches = material_dim_matches(file_restart, 'mt_leaf_iso_storage')
      restart_complete = restart_complete .and. field_matches
      IF (.not. restart_complete) &
         CALL CoLM_stop('incomplete or malformed committed generic land material restart')

      CALL read_transport_patch_field(file_restart, 'mt_ldew_rain', mt_ldew_rain)
      CALL read_transport_patch_field(file_restart, 'mt_ldew_snow', mt_ldew_snow)
      CALL read_transport_soilsnow_field(file_restart, 'mt_wliq_soisno', &
         nl_soil-maxsnl, mt_wliq_soisno)
      CALL read_transport_soilsnow_field(file_restart, 'mt_wice_soisno', &
         nl_soil-maxsnl, mt_wice_soisno)
      CALL read_transport_soilsnow_field(file_restart, 'mt_solid_soisno', &
         nl_soil-maxsnl, mt_solid_soisno)
      CALL read_transport_patch_field(file_restart, 'mt_wa', mt_wa)
      IF (restart_schema == LAND_MT_RESTART_SCHEMA_VERSION) THEN
         CALL ncio_read_vector(file_restart, 'mt_aquifer_ref_water', landpatch, &
            mt_aquifer_ref_water)
         CALL read_transport_patch_field(file_restart, 'mt_aquifer_ref_mass', mt_aquifer_ref_mass)
      ELSE
         ! Schema 4 has no isotope reference; only non-VSF-isotope generic
         ! transport is admitted here, preserving its exact old hot state.
         mt_aquifer_ref_water = 0._r8
         mt_aquifer_ref_mass = 0._r8
      ENDIF
      CALL read_transport_patch_field(file_restart, 'mt_wdsrf', mt_wdsrf)
      CALL read_transport_patch_field(file_restart, 'mt_wetwat', mt_wetwat)
      CALL read_transport_patch_field(file_restart, 'mt_surface_residue', mt_surface_residue)
      CALL read_transport_patch_field(file_restart, 'mt_subsurface_residue', mt_subsurface_residue)
      CALL read_transport_patch_field(file_restart, 'mt_canopy_solid', mt_canopy_solid)
      CALL read_transport_patch_field(file_restart, 'mt_surface_solid', mt_surface_solid)
      CALL read_transport_patch_field(file_restart, 'mt_subsurface_solid', mt_subsurface_solid)
      CALL read_transport_patch_field(file_restart, 'mt_waterstorage_solid', mt_waterstorage_solid)
      CALL read_transport_patch_field(file_restart, 'mt_scv', mt_scv)
      CALL read_transport_patch_field(file_restart, 'mt_waterstorage', mt_waterstorage)
      CALL read_transport_patch_field(file_restart, 'mt_leaf_delta_e', mt_leaf_delta_e)
      CALL read_transport_patch_field(file_restart, 'mt_leaf_delta_b', mt_leaf_delta_b)
      CALL read_transport_patch_field(file_restart, 'mt_leaf_peclet', mt_leaf_peclet)
      CALL read_transport_patch_field(file_restart, 'mt_leaf_water_moles', mt_leaf_water_moles)
      CALL read_transport_patch_field(file_restart, 'mt_leaf_iso_storage', mt_leaf_iso_storage)

      CALL validate_land_mt_restart_state(wa)
      CALL zero_provider_owned_land_mt_state()
      found_restart = .true.
   END SUBROUTINE read_land_mt_restart

   !-------------------------------------------------------------------
   ! Rebuild only mt_scv from the current scv state, leaving every
   ! other material pool as-is. Used after read_land_mt_restart when
   ! the restart file was otherwise complete but lacked mt_scv (old
   ! format). Uses the same snl-from-water heuristic as
   ! mt_init_from_water to decide when mt_scv is meaningful.
   !-------------------------------------------------------------------
   SUBROUTINE material_init_scv_from_water (numpatch, maxsnl, nl_soil, &
      wliq_soisno, wice_soisno, scv)
      IMPLICIT NONE
      integer,  intent(in) :: numpatch, maxsnl, nl_soil
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: scv(numpatch)

      integer  :: imt, ip, j, snl_local
      real(r8) :: R_init

      IF (nmts <= 0) RETURN

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         R_init = material_init_water_ratio(imt)
         DO ip = 1, numpatch
            snl_local = 0
            DO j = 0, maxsnl + 1, -1
               IF (wliq_soisno(j, ip) + wice_soisno(j, ip) > 0._r8) THEN
                  snl_local = snl_local - 1
               ELSE
                  EXIT
               ENDIF
            ENDDO
            IF (snl_local == 0) THEN
               mt_scv(imt, ip) = max(scv(ip), 0._r8) * R_init
            ELSE
               mt_scv(imt, ip) = 0._r8
            ENDIF
         ENDDO
      ENDDO
   END SUBROUTINE material_init_scv_from_water

   SUBROUTINE material_init_waterstorage_from_ratio (numpatch, maxsnl, nl_soil, &
      ldew_rain, ldew_snow, wliq_soisno, wice_soisno, &
      wa, wdsrf, wetwat, scv, waterstorage)
      IMPLICIT NONE
      integer,  intent(in) :: numpatch, maxsnl, nl_soil
      real(r8), intent(in) :: ldew_rain(numpatch)
      real(r8), intent(in) :: ldew_snow(numpatch)
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wa(numpatch)
      real(r8), intent(in) :: wdsrf(numpatch)
      real(r8), intent(in) :: wetwat(numpatch)
      real(r8), intent(in) :: scv(numpatch)
      real(r8), intent(in) :: waterstorage(numpatch)

      integer  :: imt, ip, j, snl_local
      real(r8) :: R_init, ratio, water_ref, material_ref, aquifer_water_ref, aquifer_mass_ref

      IF (nmts <= 0) RETURN
      IF (.not. allocated(mt_waterstorage)) RETURN

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         R_init = material_init_water_ratio(imt)
         DO ip = 1, numpatch
            snl_local = 0
            DO j = 0, maxsnl + 1, -1
               IF (wliq_soisno(j, ip) + wice_soisno(j, ip) > 0._r8) THEN
                  snl_local = snl_local - 1
               ELSE
                  EXIT
               ENDIF
            ENDDO
            aquifer_water_ref = max(wa(ip), 0._r8)
            aquifer_mass_ref = max(mt_wa(imt, ip), 0._r8)
            IF (material_is_isotope(imt) .and. mt_aquifer_ref_water(ip) > 0._r8) THEN
               aquifer_water_ref = max(wa(ip) + mt_aquifer_ref_water(ip), 0._r8)
               aquifer_mass_ref = max(mt_wa(imt, ip) + mt_aquifer_ref_mass(imt, ip), 0._r8)
            ENDIF
            water_ref = max(ldew_rain(ip), 0._r8) + max(ldew_snow(ip), 0._r8) &
                      + aquifer_water_ref + max(wdsrf(ip), 0._r8) &
                      + max(wetwat(ip), 0._r8) + max(scv(ip), 0._r8)
            material_ref = max(mt_ldew_rain(imt, ip), 0._r8) &
                       + max(mt_ldew_snow(imt, ip), 0._r8) &
                       + aquifer_mass_ref &
                       + max(mt_wdsrf(imt, ip), 0._r8) &
                       + max(mt_wetwat(imt, ip), 0._r8) &
                       + max(mt_scv(imt, ip), 0._r8)
            DO j = snl_local + 1, nl_soil
               water_ref = water_ref + max(wliq_soisno(j, ip), 0._r8)
               material_ref = material_ref + max(mt_wliq_soisno(imt, j, ip), 0._r8)
               water_ref = water_ref + max(wice_soisno(j, ip), 0._r8)
               material_ref = material_ref + max(mt_wice_soisno(imt, j, ip), 0._r8)
            ENDDO
            IF (water_ref > mt_water_min_for_ratio) THEN
               ratio = material_ref / water_ref
            ELSE
               ratio = R_init
            ENDIF
            mt_waterstorage(imt, ip) = max(waterstorage(ip), 0._r8) * ratio
            CALL mt_equilibrate_dissolved(imt, max(waterstorage(ip), 0._r8), &
               mt_waterstorage(imt, ip), mt_waterstorage_solid(imt, ip))
         ENDDO
      ENDDO
   END SUBROUTINE material_init_waterstorage_from_ratio

   SUBROUTINE material_enforce_solubility_from_water (numpatch, maxsnl, nl_soil, &
      ldew_rain, wliq_soisno, wa, wdsrf, wetwat, waterstorage)
      integer, intent(in) :: numpatch, maxsnl, nl_soil
      real(r8), intent(in) :: ldew_rain(numpatch)
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wa(numpatch), wdsrf(numpatch), wetwat(numpatch)
      real(r8), intent(in), optional :: waterstorage(numpatch)

      integer :: imt, ip, j
      real(r8) :: surface_water, surface_material

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         DO ip = 1, numpatch
            CALL mt_equilibrate_dissolved(imt, max(ldew_rain(ip), 0._r8), &
               mt_ldew_rain(imt, ip), mt_canopy_solid(imt, ip))
            DO j = maxsnl + 1, nl_soil
               CALL mt_equilibrate_dissolved(imt, max(wliq_soisno(j, ip), 0._r8), &
                  mt_wliq_soisno(imt, j, ip), mt_solid_soisno(imt, j, ip))
            ENDDO
            surface_water = max(wdsrf(ip), 0._r8) + max(wetwat(ip), 0._r8)
            surface_material = mt_wdsrf(imt, ip) + mt_wetwat(imt, ip)
            CALL mt_equilibrate_dissolved(imt, surface_water, surface_material, &
               mt_surface_solid(imt, ip))
            IF (surface_water > mt_water_min_for_ratio) THEN
               mt_wdsrf(imt, ip) = surface_material * max(wdsrf(ip), 0._r8) / surface_water
               mt_wetwat(imt, ip) = surface_material * max(wetwat(ip), 0._r8) / surface_water
            ELSE
               mt_wdsrf(imt, ip) = 0._r8
               mt_wetwat(imt, ip) = 0._r8
            ENDIF
            CALL mt_equilibrate_dissolved(imt, wa(ip), mt_wa(imt, ip), &
               mt_subsurface_solid(imt, ip))
            IF (present(waterstorage) .and. allocated(mt_waterstorage)) THEN
               CALL mt_equilibrate_dissolved(imt, max(waterstorage(ip), 0._r8), &
                  mt_waterstorage(imt, ip), mt_waterstorage_solid(imt, ip))
            ENDIF
         ENDDO
      ENDDO
   END SUBROUTINE material_enforce_solubility_from_water

   SUBROUTINE write_land_mt_restart (file_restart, maxsnl, nl_soil, numpatch, &
      ldew_rain, ldew_snow, wliq_soisno, wice_soisno, wa, wdsrf, wetwat, scv, waterstorage)
      USE MOD_SPMD_Task, only: p_is_worker
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer, intent(in) :: maxsnl, nl_soil, numpatch
      real(r8), intent(in) :: ldew_rain(numpatch)
      real(r8), intent(in) :: ldew_snow(numpatch)
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wa(numpatch)
      real(r8), intent(in) :: wdsrf(numpatch)
      real(r8), intent(in) :: wetwat(numpatch)
      real(r8), intent(in) :: scv(numpatch)
      real(r8), intent(in), optional :: waterstorage(numpatch)

      integer :: ntransport
      real(r8), allocatable :: restart_patch(:,:)
      real(r8), allocatable :: restart_soilsnow(:,:,:)
      logical :: have_patch_data

      CALL zero_provider_owned_land_mt_state()
      ntransport = land_transport_mt_count()
      IF (ntransport <= 0) THEN
         ! Commit an explicit empty generic transaction.  The count, not any
         ! stale descriptor/vector left in a reused file, proves that this
         ! generation owns no generic state.  Never create a zero-length NetCDF
         ! dimension; provider-owned state is persisted by lifecycle callbacks.
         CALL write_land_mt_transaction_marker(file_restart, 0)
         CALL write_land_mt_descriptor_metadata(file_restart)
         CALL write_land_mt_transaction_marker(file_restart, 1)
         RETURN
      ENDIF

      ! Reject corrupt in-memory state before publishing the in-progress
      ! marker. A crash after marker=0 can never masquerade as a committed load.
      CALL validate_land_mt_restart_state(wa)
      CALL write_land_mt_transaction_marker(file_restart, 0)
      ! The packed vectors use mt_land_transport as a NetCDF dimension.
      ! Define their descriptor/dimension before the first vector write.
      CALL write_land_mt_descriptor_metadata(file_restart)

      have_patch_data = p_is_worker .and. numpatch > 0
      allocate(restart_patch(ntransport, numpatch))

      restart_patch = 0._r8
      IF (have_patch_data .and. allocated(mt_ldew_rain)) &
         CALL pack_transport_patch(mt_ldew_rain, restart_patch)
      CALL ncio_write_vector(file_restart, 'mt_ldew_rain', 'mt_land_transport', ntransport, &
         'patch', landpatch, restart_patch, DEF_REST_CompressLevel)

      restart_patch = 0._r8
      IF (have_patch_data .and. allocated(mt_ldew_snow)) &
         CALL pack_transport_patch(mt_ldew_snow, restart_patch)
      CALL ncio_write_vector(file_restart, 'mt_ldew_snow', 'mt_land_transport', ntransport, &
         'patch', landpatch, restart_patch, DEF_REST_CompressLevel)

      allocate(restart_soilsnow(ntransport, nl_soil-maxsnl, numpatch))
      restart_soilsnow = 0._r8
      IF (have_patch_data .and. allocated(mt_wliq_soisno)) &
         CALL pack_transport_soilsnow(mt_wliq_soisno, restart_soilsnow)
      CALL ncio_write_vector(file_restart, 'mt_wliq_soisno', 'mt_land_transport', ntransport, &
         'soilsnow', nl_soil-maxsnl, 'patch', landpatch, restart_soilsnow, DEF_REST_CompressLevel)

      restart_soilsnow = 0._r8
      IF (have_patch_data .and. allocated(mt_wice_soisno)) &
         CALL pack_transport_soilsnow(mt_wice_soisno, restart_soilsnow)
      CALL ncio_write_vector(file_restart, 'mt_wice_soisno', 'mt_land_transport', ntransport, &
         'soilsnow', nl_soil-maxsnl, 'patch', landpatch, restart_soilsnow, DEF_REST_CompressLevel)

      restart_soilsnow = 0._r8
      IF (have_patch_data .and. allocated(mt_solid_soisno)) &
         CALL pack_transport_soilsnow(mt_solid_soisno, restart_soilsnow)
      CALL ncio_write_vector(file_restart, 'mt_solid_soisno', 'mt_land_transport', ntransport, &
         'soilsnow', nl_soil-maxsnl, 'patch', landpatch, restart_soilsnow, DEF_REST_CompressLevel)
      deallocate(restart_soilsnow)

      restart_patch = 0._r8
      IF (have_patch_data .and. allocated(mt_wa)) CALL pack_transport_patch(mt_wa, restart_patch)
      CALL ncio_write_vector(file_restart, 'mt_wa', 'mt_land_transport', ntransport, &
         'patch', landpatch, restart_patch, DEF_REST_CompressLevel)

      CALL ncio_write_vector(file_restart, 'mt_aquifer_ref_water', 'patch', landpatch, &
         mt_aquifer_ref_water, DEF_REST_CompressLevel)
      restart_patch = 0._r8
      IF (have_patch_data .and. allocated(mt_aquifer_ref_mass)) &
         CALL pack_transport_patch(mt_aquifer_ref_mass, restart_patch)
      CALL ncio_write_vector(file_restart, 'mt_aquifer_ref_mass', 'mt_land_transport', ntransport, &
         'patch', landpatch, restart_patch, DEF_REST_CompressLevel)

      restart_patch = 0._r8
      IF (have_patch_data .and. allocated(mt_wdsrf)) CALL pack_transport_patch(mt_wdsrf, restart_patch)
      CALL ncio_write_vector(file_restart, 'mt_wdsrf', 'mt_land_transport', ntransport, &
         'patch', landpatch, restart_patch, DEF_REST_CompressLevel)

      restart_patch = 0._r8
      IF (have_patch_data .and. allocated(mt_wetwat)) CALL pack_transport_patch(mt_wetwat, restart_patch)
      CALL ncio_write_vector(file_restart, 'mt_wetwat', 'mt_land_transport', ntransport, &
         'patch', landpatch, restart_patch, DEF_REST_CompressLevel)

      restart_patch = 0._r8
      IF (have_patch_data .and. allocated(mt_surface_residue)) &
         CALL pack_transport_patch(mt_surface_residue, restart_patch)
      CALL ncio_write_vector(file_restart, 'mt_surface_residue', 'mt_land_transport', ntransport, &
         'patch', landpatch, restart_patch, DEF_REST_CompressLevel)

      restart_patch = 0._r8
      IF (have_patch_data .and. allocated(mt_subsurface_residue)) &
         CALL pack_transport_patch(mt_subsurface_residue, restart_patch)
      CALL ncio_write_vector(file_restart, 'mt_subsurface_residue', 'mt_land_transport', ntransport, &
         'patch', landpatch, restart_patch, DEF_REST_CompressLevel)

      restart_patch = 0._r8
      IF (have_patch_data .and. allocated(mt_canopy_solid)) &
         CALL pack_transport_patch(mt_canopy_solid, restart_patch)
      CALL ncio_write_vector(file_restart, 'mt_canopy_solid', 'mt_land_transport', ntransport, &
         'patch', landpatch, restart_patch, DEF_REST_CompressLevel)

      restart_patch = 0._r8
      IF (have_patch_data .and. allocated(mt_surface_solid)) &
         CALL pack_transport_patch(mt_surface_solid, restart_patch)
      CALL ncio_write_vector(file_restart, 'mt_surface_solid', 'mt_land_transport', ntransport, &
         'patch', landpatch, restart_patch, DEF_REST_CompressLevel)

      restart_patch = 0._r8
      IF (have_patch_data .and. allocated(mt_subsurface_solid)) &
         CALL pack_transport_patch(mt_subsurface_solid, restart_patch)
      CALL ncio_write_vector(file_restart, 'mt_subsurface_solid', 'mt_land_transport', ntransport, &
         'patch', landpatch, restart_patch, DEF_REST_CompressLevel)

      restart_patch = 0._r8
      IF (have_patch_data .and. allocated(mt_waterstorage_solid)) &
         CALL pack_transport_patch(mt_waterstorage_solid, restart_patch)
      CALL ncio_write_vector(file_restart, 'mt_waterstorage_solid', 'mt_land_transport', ntransport, &
         'patch', landpatch, restart_patch, DEF_REST_CompressLevel)

      restart_patch = 0._r8
      IF (have_patch_data .and. allocated(mt_scv)) CALL pack_transport_patch(mt_scv, restart_patch)
      CALL ncio_write_vector(file_restart, 'mt_scv', 'mt_land_transport', ntransport, &
         'patch', landpatch, restart_patch, DEF_REST_CompressLevel)

      restart_patch = 0._r8
      IF (have_patch_data .and. allocated(mt_waterstorage)) &
         CALL pack_transport_patch(mt_waterstorage, restart_patch)
      CALL ncio_write_vector(file_restart, 'mt_waterstorage', 'mt_land_transport', ntransport, &
         'patch', landpatch, restart_patch, DEF_REST_CompressLevel)

      restart_patch = 0._r8
      IF (have_patch_data .and. allocated(mt_leaf_delta_e)) &
         CALL pack_transport_patch(mt_leaf_delta_e, restart_patch)
      CALL ncio_write_vector(file_restart, 'mt_leaf_delta_e', 'mt_land_transport', ntransport, &
         'patch', landpatch, restart_patch, DEF_REST_CompressLevel)

      restart_patch = 0._r8
      IF (have_patch_data .and. allocated(mt_leaf_delta_b)) &
         CALL pack_transport_patch(mt_leaf_delta_b, restart_patch)
      CALL ncio_write_vector(file_restart, 'mt_leaf_delta_b', 'mt_land_transport', ntransport, &
         'patch', landpatch, restart_patch, DEF_REST_CompressLevel)

      restart_patch = 0._r8
      IF (have_patch_data .and. allocated(mt_leaf_peclet)) &
         CALL pack_transport_patch(mt_leaf_peclet, restart_patch)
      CALL ncio_write_vector(file_restart, 'mt_leaf_peclet', 'mt_land_transport', ntransport, &
         'patch', landpatch, restart_patch, DEF_REST_CompressLevel)

      restart_patch = 0._r8
      IF (have_patch_data .and. allocated(mt_leaf_water_moles)) &
         CALL pack_transport_patch(mt_leaf_water_moles, restart_patch)
      CALL ncio_write_vector(file_restart, 'mt_leaf_water_moles', 'mt_land_transport', ntransport, &
         'patch', landpatch, restart_patch, DEF_REST_CompressLevel)

      restart_patch = 0._r8
      IF (have_patch_data .and. allocated(mt_leaf_iso_storage)) &
         CALL pack_transport_patch(mt_leaf_iso_storage, restart_patch)
      CALL ncio_write_vector(file_restart, 'mt_leaf_iso_storage', 'mt_land_transport', ntransport, &
         'patch', landpatch, restart_patch, DEF_REST_CompressLevel)

      deallocate(restart_patch)
      CALL material_forcing_write_restart(file_restart)
      CALL write_land_mt_transaction_marker(file_restart, 1)
   END SUBROUTINE write_land_mt_restart

   SUBROUTINE write_mt_restart_all (file_restart, maxsnl, nl_soil, numpatch, &
      ldew_rain, ldew_snow, wliq_soisno, wice_soisno, wa, wdsrf, wetwat, scv, &
      compress, waterstorage)
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer, intent(in) :: maxsnl, nl_soil, numpatch, compress
      real(r8), intent(in) :: ldew_rain(numpatch)
      real(r8), intent(in) :: ldew_snow(numpatch)
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wa(numpatch)
      real(r8), intent(in) :: wdsrf(numpatch)
      real(r8), intent(in) :: wetwat(numpatch)
      real(r8), intent(in) :: scv(numpatch)
      real(r8), intent(in), optional :: waterstorage(numpatch)

      IF (present(waterstorage)) THEN
         CALL write_land_mt_restart(file_restart, maxsnl, nl_soil, numpatch, &
            ldew_rain, ldew_snow, wliq_soisno, wice_soisno, wa, wdsrf, wetwat, scv, waterstorage)
      ELSE
         CALL write_land_mt_restart(file_restart, maxsnl, nl_soil, numpatch, &
            ldew_rain, ldew_snow, wliq_soisno, wice_soisno, wa, wdsrf, wetwat, scv)
      ENDIF

      CALL material_lifecycle_land_write_restart(file_restart, compress)

   END SUBROUTINE write_mt_restart_all

END MODULE MOD_MT_Rest
#endif
