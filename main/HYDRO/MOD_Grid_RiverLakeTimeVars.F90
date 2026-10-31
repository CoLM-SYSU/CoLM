#include <define.h>

#ifdef GridRiverLakeFlow
MODULE MOD_Grid_RiverLakeTimeVars
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   Time Variables in gridded hydrological processes.
!
! Created by Shupeng Zhang, Oct 2025
!-------------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Grid_RiverLakeBifurcation, only: write_bifurcation_restart
   IMPLICIT NONE

   integer, parameter :: GRIDRIVER_RESTART_SCHEMA_VERSION = 2
   integer, parameter :: GRIDRIVER_UCATCH_IDENTITY_VERSION = 1
   integer, parameter :: GRIDRIVER_RESERVOIR_IDENTITY_VERSION = 1

   ! -- state variables --
   real(r8), allocatable, target :: wdsrf_ucat (:)
   real(r8), allocatable, target :: wdsrf_ucat_prev (:)
   logical :: wdsrf_ucat_prev_valid = .false.
   logical :: wdsrf_ucat_prev_restart_found = .false.
   logical :: restart_transaction_validated = .false.
   logical :: restart_feature_manifest_present = .false.
   logical :: restart_bifurcation_enabled = .false.
   logical :: restart_levee_enabled = .false.
   real(r8), allocatable, target :: veloc_riv  (:)
   real(r8), allocatable :: momen_riv  (:) ! unit river momentum       [m^2/s]
   real(r8), allocatable :: volresv    (:) ! reservoir water volume    [m^3]
   ! release of the reservoirs directly upstream that reached each reservoir in
   ! the previous routing sub-step [m^3/s]; part of the inflow seen by the
   ! reservoir operation (see grid_riverlake_flow)
   real(r8), allocatable :: qresv_ups_prev (:)
   real(r8), allocatable :: volwater_ucat (:)
   logical               :: volwater_ucat_valid = .false.

   real(r8), save       :: acctime_rnof = 0._r8
   real(r8), allocatable :: acc_rnof_uc (:)

   ! -- restart file path (saved for the deferred levee, bifurcation and material reads) --
   character(len=512) :: gridriver_restart_file = ''

   ! -- running state set aside while LULCC re-initializes (see hold_GridRiverLakeTimeVars_lulcc) --
   real(r8), allocatable :: held_wdsrf_ucat (:), held_wdsrf_ucat_prev (:), held_veloc_riv (:)
   real(r8), allocatable :: held_momen_riv (:), held_volresv (:), held_volwater_ucat (:)
   real(r8), allocatable :: held_acc_rnof_uc (:), held_qresv_ups_prev (:)
   real(r8)              :: held_acctime_rnof
   logical               :: held_flags (7)
   character(len=512)    :: held_restart_file

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_GridRiverLakeTimeVars
   PUBLIC :: deallocate_GridRiverLakeTimeVars
   PUBLIC :: hold_GridRiverLakeTimeVars_lulcc
   PUBLIC :: restore_GridRiverLakeTimeVars_lulcc

   PUBLIC :: read_GridRiverLakeTimeVars
   PUBLIC :: write_GridRiverLakeTimeVars
   PUBLIC :: commit_GridRiverLakeRestart

CONTAINS

   SUBROUTINE build_gridriver_ucatch_identity (identity)

   USE MOD_SPMD_Task
   USE MOD_Grid_RiverLakeNetwork, only: numucat, x_ucat, y_ucat, ucat_next
   IMPLICIT NONE

   real(r8), allocatable, intent(out) :: identity(:,:)

      allocate (identity(4, numucat))
      identity = 0._r8

      IF (p_is_worker .and. numucat > 0) THEN
         identity(1, :) = real(GRIDRIVER_UCATCH_IDENTITY_VERSION, r8)
         identity(2, :) = real(x_ucat, r8)
         identity(3, :) = real(y_ucat, r8)
         identity(4, :) = real(ucat_next, r8)
      ENDIF

   END SUBROUTINE build_gridriver_ucatch_identity


   SUBROUTINE write_gridriver_ucatch_identity (file_restart)

   USE MOD_SPMD_Task
   USE MOD_Namelist, only: DEF_REST_CompressLevel
   USE MOD_NetCDFSerial, only: ncio_define_dimension, ncio_write_serial
   USE MOD_Grid_RiverLakeVectorIO, only: vector_gather_matrix_to_master
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_ucid
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart
   integer, allocatable :: global_id_local(:)
   real(r8), allocatable :: identity(:,:), identity_global(:,:)

      IF (totalnumucat <= 0) RETURN

      CALL build_gridriver_ucatch_identity (identity)
      allocate (global_id_local(numucat))
      IF (p_is_worker .and. numucat > 0) global_id_local = ucat_ucid

      CALL vector_gather_matrix_to_master (identity, 4, numucat, totalnumucat, &
         global_id_local, identity_global)

      IF (p_is_master) THEN
         CALL ncio_define_dimension (file_restart, 'gridriver_ucatch_identity_field', 4)
         CALL ncio_write_serial (file_restart, 'gridriver_ucatch_identity', identity_global, &
            'gridriver_ucatch_identity_field', 'ucatch', DEF_REST_CompressLevel)
         deallocate (identity_global)
      ENDIF

      deallocate (identity)
      deallocate (global_id_local)

   END SUBROUTINE write_gridriver_ucatch_identity


   SUBROUTINE build_gridriver_reservoir_identity (identity)

   USE MOD_SPMD_Task
   USE MOD_Grid_Reservoir, only: numresv, resv_ucid
   IMPLICIT NONE

   real(r8), allocatable, intent(out) :: identity(:,:)

      allocate (identity(2, numresv))
      identity = 0._r8

      IF (p_is_worker .and. numresv > 0) THEN
         identity(1, :) = real(GRIDRIVER_RESERVOIR_IDENTITY_VERSION, r8)
         identity(2, :) = real(resv_ucid, r8)
      ENDIF

   END SUBROUTINE build_gridriver_reservoir_identity


   SUBROUTINE write_gridriver_reservoir_identity (file_restart)

   USE MOD_SPMD_Task
   USE MOD_Namelist, only: DEF_REST_CompressLevel
   USE MOD_NetCDFSerial, only: ncio_define_dimension, ncio_write_serial
   USE MOD_Grid_RiverLakeVectorIO, only: vector_gather_matrix_to_master
   USE MOD_Grid_Reservoir, only: numresv, totalnumresv, resv_global_id
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart
   integer, allocatable :: global_id_local(:)
   real(r8), allocatable :: identity(:,:), identity_global(:,:)

      IF (totalnumresv <= 0) RETURN

      CALL build_gridriver_reservoir_identity (identity)
      allocate (global_id_local(numresv))
      IF (p_is_worker .and. numresv > 0) global_id_local = resv_global_id

      CALL vector_gather_matrix_to_master (identity, 2, numresv, totalnumresv, &
         global_id_local, identity_global)

      IF (p_is_master) THEN
         CALL ncio_define_dimension (file_restart, 'gridriver_reservoir_identity_field', 2)
         CALL ncio_write_serial (file_restart, 'gridriver_reservoir_identity', identity_global, &
            'gridriver_reservoir_identity_field', 'reservoir', DEF_REST_CompressLevel)
         deallocate (identity_global)
      ENDIF

      deallocate (identity)
      deallocate (global_id_local)

   END SUBROUTINE write_gridriver_reservoir_identity


   SUBROUTINE validate_gridriver_restart_transaction (file_restart, legacy_restart)

   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial, only: ncio_var_exist, ncio_read_bcast_serial
   USE MOD_Grid_RiverLakeNetwork, only: totalnumucat
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart
   logical, intent(out) :: legacy_restart
   integer :: has_flags(5), restart_schema, restart_complete, identity_expected
   integer :: feature_bifurcation, feature_levee

      has_flags = 0
      restart_transaction_validated = .false.
      restart_feature_manifest_present = .false.
      restart_bifurcation_enabled = .false.
      restart_levee_enabled = .false.
      IF (p_is_master) THEN
         IF (ncio_var_exist(file_restart, 'gridriver_restart_schema', readflag = .false.)) has_flags(1) = 1
         IF (ncio_var_exist(file_restart, 'gridriver_restart_complete', readflag = .false.)) has_flags(2) = 1
         IF (ncio_var_exist(file_restart, 'gridriver_ucatch_identity', readflag = .false.)) has_flags(3) = 1
         IF (ncio_var_exist(file_restart, 'gridriver_restart_feature_bifurcation', &
            readflag = .false.)) has_flags(4) = 1
         IF (ncio_var_exist(file_restart, 'gridriver_restart_feature_levee', &
            readflag = .false.)) has_flags(5) = 1
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast (has_flags, 5, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif

      legacy_restart = all(has_flags == 0)
      IF (legacy_restart) THEN
         IF (p_is_master) THEN
            write(*,'(A)') 'WARNING: legacy GridRiverLake restart has no transaction or ucatch identity metadata.'
            call flush(6)
         ENDIF
         restart_transaction_validated = .true.
         RETURN
      ENDIF

      identity_expected = merge(1, 0, totalnumucat > 0)
      IF (has_flags(1) /= 1 .or. has_flags(2) /= 1 .or. has_flags(3) /= identity_expected) THEN
         CALL CoLM_stop ('incomplete GridRiverLake restart transaction metadata')
      ENDIF

      CALL ncio_read_bcast_serial (file_restart, 'gridriver_restart_schema', restart_schema)
      CALL ncio_read_bcast_serial (file_restart, 'gridriver_restart_complete', restart_complete)
      IF (restart_schema /= 1 .and. restart_schema /= GRIDRIVER_RESTART_SCHEMA_VERSION) THEN
         CALL CoLM_stop ('unsupported GridRiverLake restart schema')
      ENDIF
      IF (restart_complete /= 1) THEN
         CALL CoLM_stop ('GridRiverLake restart transaction is not complete')
      ENDIF

      IF (restart_schema == 1) THEN
         restart_transaction_validated = .true.
         RETURN
      ENDIF
      IF (restart_schema == GRIDRIVER_RESTART_SCHEMA_VERSION) THEN
         IF (has_flags(4) /= 1 .or. has_flags(5) /= 1) THEN
            CALL CoLM_stop ('GridRiverLake schema-v2 restart is missing feature manifest')
         ENDIF
         CALL ncio_read_bcast_serial (file_restart, &
            'gridriver_restart_feature_bifurcation', feature_bifurcation)
         CALL ncio_read_bcast_serial (file_restart, &
            'gridriver_restart_feature_levee', feature_levee)
         IF ((feature_bifurcation /= 0 .and. feature_bifurcation /= 1) .or. &
             (feature_levee /= 0 .and. feature_levee /= 1)) THEN
            CALL CoLM_stop ('invalid GridRiverLake restart feature manifest')
         ENDIF
         restart_feature_manifest_present = .true.
         restart_bifurcation_enabled = feature_bifurcation == 1
         restart_levee_enabled = feature_levee == 1
      ENDIF
      restart_transaction_validated = .true.

   END SUBROUTINE validate_gridriver_restart_transaction


   SUBROUTINE validate_gridriver_ucatch_identity (file_restart, legacy_restart)

   USE MOD_SPMD_Task
   USE MOD_Grid_RiverLakeVectorIO, only: vector_read_matrix_and_scatter
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_ucid
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart
   logical, intent(in) :: legacy_restart
   integer :: mismatch_count, first_mismatch_gid, i
   integer, allocatable :: global_id_local(:)
   real(r8), allocatable :: identity(:,:), identity_restart(:,:)

      IF (legacy_restart) RETURN
      IF (totalnumucat <= 0) RETURN

      CALL build_gridriver_ucatch_identity (identity)
      allocate (global_id_local(numucat))
      IF (p_is_worker .and. numucat > 0) global_id_local = ucat_ucid

      CALL vector_read_matrix_and_scatter (file_restart, identity_restart, 4, numucat, &
         'gridriver_ucatch_identity', global_id_local, totalnumucat)

      mismatch_count = 0
      first_mismatch_gid = huge(first_mismatch_gid)
      IF (p_is_worker) THEN
         DO i = 1, numucat
            IF (any(.not. ieee_is_finite(identity_restart(:, i)))) THEN
               mismatch_count = mismatch_count + 1
               first_mismatch_gid = min(first_mismatch_gid, global_id_local(i))
            ELSEIF (any(identity_restart(:, i) /= identity(:, i))) THEN
               mismatch_count = mismatch_count + 1
               first_mismatch_gid = min(first_mismatch_gid, global_id_local(i))
            ENDIF
         ENDDO
      ENDIF
#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, mismatch_count, 1, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, first_mismatch_gid, 1, MPI_INTEGER, MPI_MIN, p_comm_glb, p_err)
#endif

      IF (allocated(identity_restart)) deallocate (identity_restart)
      deallocate (identity)
      deallocate (global_id_local)

      IF (mismatch_count > 0) THEN
         IF (p_is_master) write(*,'(A,I0,A,I0)') &
            'ERROR: GridRiverLake restart ucatch identity mismatch; count=', &
            mismatch_count, ', first global ucatch=', first_mismatch_gid
         CALL CoLM_stop ('Refusing to load GridRiverLake state for a different ucatch network')
      ENDIF

   END SUBROUTINE validate_gridriver_ucatch_identity


   SUBROUTINE validate_gridriver_reservoir_identity (file_restart, legacy_restart)

   USE MOD_SPMD_Task
   USE MOD_Namelist, only: DEF_Reservoir_Method
   USE MOD_NetCDFSerial, only: ncio_var_exist
   USE MOD_Grid_RiverLakeVectorIO, only: vector_read_matrix_and_scatter
   USE MOD_Grid_Reservoir, only: numresv, totalnumresv, resv_global_id
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart
   logical, intent(in) :: legacy_restart
   logical :: has_identity
   integer :: mismatch_count, first_mismatch_gid, i
   integer, allocatable :: global_id_local(:)
   real(r8), allocatable :: identity(:,:), identity_restart(:,:)

      IF (legacy_restart) RETURN
      IF (DEF_Reservoir_Method <= 0 .or. totalnumresv <= 0) RETURN

      has_identity = .false.
      IF (p_is_master) has_identity = ncio_var_exist(&
         file_restart, 'gridriver_reservoir_identity', readflag = .false.)
#ifdef USEMPI
      CALL mpi_bcast (has_identity, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif
      IF (.not. has_identity) THEN
         IF (restart_feature_manifest_present) THEN
            CALL CoLM_stop ('GridRiverLake schema-v2 restart is missing required reservoir identity')
         ENDIF
         RETURN
      ENDIF

      CALL build_gridriver_reservoir_identity (identity)
      allocate (global_id_local(numresv))
      IF (p_is_worker .and. numresv > 0) global_id_local = resv_global_id

      CALL vector_read_matrix_and_scatter (file_restart, identity_restart, 2, numresv, &
         'gridriver_reservoir_identity', global_id_local, totalnumresv)

      mismatch_count = 0
      first_mismatch_gid = huge(first_mismatch_gid)
      IF (p_is_worker) THEN
         DO i = 1, numresv
            IF (any(.not. ieee_is_finite(identity_restart(:, i)))) THEN
               mismatch_count = mismatch_count + 1
               first_mismatch_gid = min(first_mismatch_gid, global_id_local(i))
            ELSEIF (any(identity_restart(:, i) /= identity(:, i))) THEN
               mismatch_count = mismatch_count + 1
               first_mismatch_gid = min(first_mismatch_gid, global_id_local(i))
            ENDIF
         ENDDO
      ENDIF
#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, mismatch_count, 1, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, first_mismatch_gid, 1, MPI_INTEGER, MPI_MIN, p_comm_glb, p_err)
#endif

      IF (allocated(identity_restart)) deallocate (identity_restart)
      deallocate (identity)
      deallocate (global_id_local)

      IF (mismatch_count > 0) THEN
         IF (p_is_master) write(*,'(A,I0,A,I0)') &
            'ERROR: GridRiverLake restart reservoir identity mismatch; count=', &
            mismatch_count, ', first reservoir ordinal=', first_mismatch_gid
         CALL CoLM_stop ('Refusing to load GridRiverLake state for a different reservoir configuration')
      ENDIF

   END SUBROUTINE validate_gridriver_reservoir_identity


   SUBROUTINE commit_GridRiverLakeRestart (file_restart)

   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial, only: ncio_write_serial
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) CALL ncio_write_serial (file_restart, 'gridriver_restart_complete', 1)

   END SUBROUTINE commit_GridRiverLakeRestart

   SUBROUTINE allocate_GridRiverLakeTimeVars

   USE MOD_SPMD_Task
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   USE MOD_Grid_Reservoir,        only: numresv
   IMPLICIT NONE

   integer :: ncell_state, nresv_state

      ncell_state = 0
      nresv_state = 0
      IF (p_is_worker) THEN
         ncell_state = numucat
         nresv_state = numresv
      ENDIF

      allocate (wdsrf_ucat (ncell_state))
      allocate (wdsrf_ucat_prev (ncell_state))
      allocate (veloc_riv  (ncell_state))
      allocate (momen_riv  (ncell_state))
      allocate (volresv (nresv_state))
      allocate (qresv_ups_prev (nresv_state))
      qresv_ups_prev = 0._r8
      allocate (acc_rnof_uc (ncell_state))
      acc_rnof_uc = 0._r8
      wdsrf_ucat_prev = 0._r8

      allocate (volwater_ucat(ncell_state))
      volwater_ucat = 0._r8

      acctime_rnof = 0._r8
      wdsrf_ucat_prev_valid = .false.
      wdsrf_ucat_prev_restart_found = .false.
      restart_transaction_validated = .false.
      restart_feature_manifest_present = .false.
      restart_bifurcation_enabled = .false.
      restart_levee_enabled = .false.
      volwater_ucat_valid = .false.
      gridriver_restart_file = ''

   END SUBROUTINE allocate_GridRiverLakeTimeVars


   SUBROUTINE READ_GridRiverLakeTimeVars (file_restart)

   USE MOD_SPMD_Task
   USE MOD_Namelist,              only: DEF_Reservoir_Method
   USE MOD_Vector_ReadWrite
   USE MOD_NetCDFSerial,          only: ncio_var_exist, ncio_read_bcast_serial
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_data_address, &
      totalnpthout, npthlev_bif
   USE MOD_Grid_Reservoir,        only: numresv, resv_data_address, totalnumresv
   USE MOD_Vars_Global,           only: spval
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart
   logical :: has_var, legacy_restart, restore_previous_depth
   integer :: base_var_flags(6), invalid_prev_count, invalid_base_count, i

      gridriver_restart_file = trim(file_restart)
      CALL validate_gridriver_restart_transaction (file_restart, legacy_restart)
      CALL validate_gridriver_ucatch_identity (file_restart, legacy_restart)
      CALL validate_gridriver_reservoir_identity (file_restart, legacy_restart)

      acctime_rnof = 0._r8
      IF (allocated(acc_rnof_uc)) acc_rnof_uc = 0._r8
      IF (allocated(wdsrf_ucat_prev)) wdsrf_ucat_prev = 0._r8
      wdsrf_ucat_prev_valid = .false.
      wdsrf_ucat_prev_restart_found = .false.
      IF (allocated(volwater_ucat)) volwater_ucat = 0._r8
      volwater_ucat_valid = .false.

      base_var_flags = 0
      IF (p_is_master) THEN
         IF (ncio_var_exist(file_restart, 'wdsrf_ucat', readflag = .false.)) base_var_flags(1) = 1
         IF (ncio_var_exist(file_restart, 'veloc_riv', readflag = .false.)) base_var_flags(2) = 1
         IF (ncio_var_exist(file_restart, 'acctime_rnof', readflag = .false.)) base_var_flags(3) = 1
         IF (ncio_var_exist(file_restart, 'acc_rnof_uc', readflag = .false.)) base_var_flags(4) = 1
         IF (ncio_var_exist(file_restart, 'volwater_ucat', readflag = .false.)) base_var_flags(5) = 1
         IF (ncio_var_exist(file_restart, 'volresv', readflag = .false.)) base_var_flags(6) = 1
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast (base_var_flags, 6, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif
      IF (restart_feature_manifest_present) THEN
         IF (base_var_flags(3) /= 1) THEN
            CALL CoLM_stop ('GridRiverLake schema-v2 restart is missing required base state')
         ENDIF
         IF (totalnumucat > 0) THEN
            IF (base_var_flags(1) /= 1 .or. base_var_flags(2) /= 1 .or. &
                base_var_flags(4) /= 1 .or. base_var_flags(5) /= 1) THEN
               CALL CoLM_stop ('GridRiverLake schema-v2 restart is missing required base state')
            ENDIF
         ENDIF
         IF (DEF_Reservoir_Method > 0 .and. totalnumresv > 0 .and. base_var_flags(6) /= 1) THEN
            CALL CoLM_stop ('GridRiverLake schema-v2 restart is missing required base state')
         ENDIF
      ENDIF

      IF (totalnumucat > 0) THEN
         CALL vector_read_and_scatter (file_restart, wdsrf_ucat, numucat, 'wdsrf_ucat', ucat_data_address)
         CALL vector_read_and_scatter (file_restart, veloc_riv,  numucat, 'veloc_riv',  ucat_data_address)
      ENDIF

      has_var = .false.
      IF (p_is_master .and. totalnumucat > 0) &
         has_var = ncio_var_exist(file_restart, 'wdsrf_ucat_prev', readflag = .false.)
#ifdef USEMPI
      CALL mpi_bcast (has_var, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif
      restore_previous_depth = restart_bifurcation_enabled .and. &
         totalnpthout > 0 .and. npthlev_bif > 0
      IF (restart_feature_manifest_present .and. .not. restore_previous_depth) has_var = .false.
      IF (.not. has_var .and. .not. legacy_restart .and. totalnumucat > 0 .and. &
          (.not. restart_feature_manifest_present .or. restore_previous_depth)) THEN
         CALL CoLM_stop ('GridRiverLake restart transaction is missing wdsrf_ucat_prev')
      ENDIF
      IF (has_var) THEN
         CALL vector_read_and_scatter (file_restart, wdsrf_ucat_prev, numucat, &
            'wdsrf_ucat_prev', ucat_data_address)

         invalid_prev_count = 0
         DO i = 1, size(wdsrf_ucat_prev)
            IF (.not. ieee_is_finite(wdsrf_ucat_prev(i))) THEN
               invalid_prev_count = invalid_prev_count + 1
            ELSEIF (wdsrf_ucat_prev(i) < 0._r8) THEN
               invalid_prev_count = invalid_prev_count + 1
            ENDIF
         ENDDO
#ifdef USEMPI
         CALL mpi_allreduce (MPI_IN_PLACE, invalid_prev_count, 1, MPI_INTEGER, &
            MPI_SUM, p_comm_glb, p_err)
#endif
         IF (invalid_prev_count > 0) THEN
            IF (.not. legacy_restart) THEN
               CALL CoLM_stop ('GridRiverLake restart has invalid wdsrf_ucat_prev')
            ENDIF
            wdsrf_ucat_prev = wdsrf_ucat
         ELSE
            wdsrf_ucat_prev_restart_found = .true.
         ENDIF
      ELSE
         wdsrf_ucat_prev = wdsrf_ucat
      ENDIF
      wdsrf_ucat_prev_valid = .true.

      IF (DEF_Reservoir_Method > 0) THEN
         IF (totalnumresv > 0) THEN
            CALL vector_read_and_scatter (file_restart, volresv, numresv, 'volresv', resv_data_address)
            ! absent in restarts written before it existed: start from no upstream release
            has_var = .false.
            IF (p_is_master) has_var = ncio_var_exist(file_restart, 'qresv_ups_prev', readflag = .false.)
#ifdef USEMPI
            CALL mpi_bcast (has_var, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif
            IF (has_var) THEN
               CALL vector_read_and_scatter (file_restart, qresv_ups_prev, numresv, &
                  'qresv_ups_prev', resv_data_address)
            ELSEIF (allocated(qresv_ups_prev)) THEN
               qresv_ups_prev = 0._r8
            ENDIF
         ENDIF
      ENDIF

      IF (base_var_flags(3) == 1) THEN
         CALL ncio_read_bcast_serial (file_restart, 'acctime_rnof', acctime_rnof)
      ENDIF

      IF (base_var_flags(4) == 1) THEN
         CALL vector_read_and_scatter (file_restart, acc_rnof_uc, numucat, 'acc_rnof_uc', ucat_data_address)
      ENDIF

      IF (base_var_flags(5) == 1) THEN
         CALL vector_read_and_scatter (file_restart, volwater_ucat, numucat, 'volwater_ucat', ucat_data_address)
      ENDIF

      invalid_base_count = 0
      DO i = 1, size(wdsrf_ucat)
         IF (.not. ieee_is_finite(wdsrf_ucat(i))) THEN
            invalid_base_count = invalid_base_count + 1
         ELSEIF (wdsrf_ucat(i) < 0._r8) THEN
            invalid_base_count = invalid_base_count + 1
         ENDIF
         IF (.not. ieee_is_finite(veloc_riv(i))) THEN
            invalid_base_count = invalid_base_count + 1
         ELSEIF (abs(veloc_riv(i)) > 50._r8) THEN
            invalid_base_count = invalid_base_count + 1
         ENDIF
         IF (base_var_flags(4) == 1) THEN
            IF (.not. ieee_is_finite(acc_rnof_uc(i))) &
               invalid_base_count = invalid_base_count + 1
         ENDIF
         IF (base_var_flags(5) == 1) THEN
            IF (.not. ieee_is_finite(volwater_ucat(i))) THEN
               invalid_base_count = invalid_base_count + 1
            ELSEIF (volwater_ucat(i) < 0._r8) THEN
               invalid_base_count = invalid_base_count + 1
            ENDIF
         ENDIF
      ENDDO
      IF (DEF_Reservoir_Method > 0 .and. totalnumresv > 0 .and. base_var_flags(6) == 1) THEN
         DO i = 1, size(volresv)
            IF (.not. ieee_is_finite(volresv(i))) THEN
               invalid_base_count = invalid_base_count + 1
            ELSEIF (volresv(i) < 0._r8 .and. volresv(i) /= spval) THEN
               invalid_base_count = invalid_base_count + 1
            ENDIF
         ENDDO
      ENDIF
      IF (p_is_master .and. base_var_flags(3) == 1) THEN
         IF (.not. ieee_is_finite(acctime_rnof)) THEN
            invalid_base_count = invalid_base_count + 1
         ELSEIF (acctime_rnof < 0._r8) THEN
            invalid_base_count = invalid_base_count + 1
         ENDIF
      ENDIF
#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, invalid_base_count, 1, MPI_INTEGER, &
         MPI_SUM, p_comm_glb, p_err)
#endif
      IF (invalid_base_count > 0) THEN
         CALL CoLM_stop ('invalid GridRiverLake restart base state')
      ENDIF
      IF (base_var_flags(5) == 1) volwater_ucat_valid = .true.

   END SUBROUTINE READ_GridRiverLakeTimeVars


   SUBROUTINE WRITE_GridRiverLakeTimeVars (file_restart)

   USE MOD_SPMD_Task
   USE MOD_Namelist,              only: DEF_Reservoir_Method, DEF_USE_LEVEE, DEF_USE_BIFURCATION
   USE MOD_NetCDFSerial
   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_data_address, &
      totalnpthout, npthlev_bif
   USE MOD_Grid_Reservoir,        only: numresv, totalnumresv, resv_data_address
   USE MOD_Grid_RiverLakeLevee,   only: levsto
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

      IF (p_is_master) THEN
         CALL ncio_create_file (trim(file_restart))
         CALL ncio_define_dimension(file_restart, 'ucatch', totalnumucat)
         CALL ncio_write_serial (file_restart, 'gridriver_restart_schema', GRIDRIVER_RESTART_SCHEMA_VERSION)
         CALL ncio_write_serial (file_restart, 'gridriver_restart_complete', 0)
         CALL ncio_write_serial (file_restart, 'gridriver_restart_feature_bifurcation', &
            merge(1, 0, DEF_USE_BIFURCATION))
         CALL ncio_write_serial (file_restart, 'gridriver_restart_feature_levee', &
            merge(1, 0, DEF_USE_LEVEE))
      ENDIF

      CALL write_gridriver_ucatch_identity (file_restart)

      CALL vector_gather_and_write (&
         wdsrf_ucat, numucat, totalnumucat, ucat_data_address, file_restart, 'wdsrf_ucat', 'ucatch')

      IF (DEF_USE_BIFURCATION .and. totalnpthout > 0 .and. npthlev_bif > 0) THEN
         IF (.not. wdsrf_ucat_prev_valid) THEN
            wdsrf_ucat_prev = wdsrf_ucat
            wdsrf_ucat_prev_valid = .true.
         ENDIF
         CALL vector_gather_and_write (&
            wdsrf_ucat_prev, numucat, totalnumucat, ucat_data_address, &
            file_restart, 'wdsrf_ucat_prev', 'ucatch')
      ENDIF

      CALL vector_gather_and_write (&
         veloc_riv, numucat, totalnumucat, ucat_data_address, file_restart, 'veloc_riv', 'ucatch')

      IF (p_is_master) CALL ncio_write_serial (file_restart, 'acctime_rnof', acctime_rnof)
      CALL vector_gather_and_write (&
         acc_rnof_uc, numucat, totalnumucat, ucat_data_address, file_restart, 'acc_rnof_uc', 'ucatch')

      IF (DEF_Reservoir_Method > 0) THEN
         IF (totalnumresv > 0) THEN

            IF (p_is_master) CALL ncio_define_dimension(file_restart, 'reservoir', totalnumresv)

            CALL write_gridriver_reservoir_identity (file_restart)

            CALL vector_gather_and_write (&
               volresv, numresv, totalnumresv, resv_data_address, file_restart, 'volresv', 'reservoir')
            CALL vector_gather_and_write (&
               qresv_ups_prev, numresv, totalnumresv, resv_data_address, file_restart, &
               'qresv_ups_prev', 'reservoir')
         ENDIF
      ENDIF

      CALL vector_gather_and_write (&
         volwater_ucat, numucat, totalnumucat, ucat_data_address, file_restart, 'volwater_ucat', 'ucatch')

      IF (DEF_USE_LEVEE) THEN
         CALL vector_gather_and_write (&
            levsto, numucat, totalnumucat, ucat_data_address, file_restart, 'levsto', 'ucatch')
      ENDIF

      IF (DEF_USE_BIFURCATION) THEN
         CALL write_bifurcation_restart(file_restart)
      ENDIF


   END SUBROUTINE WRITE_GridRiverLakeTimeVars

   SUBROUTINE deallocate_GridRiverLakeTimeVars

   IMPLICIT NONE

      IF (allocated (wdsrf_ucat )) deallocate (wdsrf_ucat )
      IF (allocated (wdsrf_ucat_prev)) deallocate (wdsrf_ucat_prev)
      IF (allocated (veloc_riv  )) deallocate (veloc_riv  )
      IF (allocated (momen_riv  )) deallocate (momen_riv  )
      IF (allocated (volresv    )) deallocate (volresv    )
      IF (allocated (qresv_ups_prev)) deallocate (qresv_ups_prev)
      IF (allocated (volwater_ucat)) deallocate (volwater_ucat)
      IF (allocated (acc_rnof_uc)) deallocate (acc_rnof_uc)
      acctime_rnof = 0._r8
      wdsrf_ucat_prev_valid = .false.
      wdsrf_ucat_prev_restart_found = .false.
      restart_transaction_validated = .false.
      restart_feature_manifest_present = .false.
      restart_bifurcation_enabled = .false.
      restart_levee_enabled = .false.
      volwater_ucat_valid = .false.
      gridriver_restart_file = ''

   END SUBROUTINE deallocate_GridRiverLakeTimeVars

   ! LULCC re-runs the whole initialization (deallocate_TimeVariables, then initialize), which
   ! frees and reallocates the river state.  The unit catchments do not change with land cover,
   ! so the running state is moved aside here and moved back by restore_GridRiverLakeTimeVars_lulcc.
   SUBROUTINE hold_GridRiverLakeTimeVars_lulcc ()

   IMPLICIT NONE

      CALL move_alloc (wdsrf_ucat     , held_wdsrf_ucat     )
      CALL move_alloc (wdsrf_ucat_prev, held_wdsrf_ucat_prev)
      CALL move_alloc (veloc_riv      , held_veloc_riv      )
      CALL move_alloc (momen_riv      , held_momen_riv      )
      CALL move_alloc (volresv        , held_volresv        )
      CALL move_alloc (volwater_ucat  , held_volwater_ucat  )
      CALL move_alloc (acc_rnof_uc    , held_acc_rnof_uc    )
      CALL move_alloc (qresv_ups_prev , held_qresv_ups_prev )
      held_acctime_rnof = acctime_rnof
      held_flags = (/ wdsrf_ucat_prev_valid, wdsrf_ucat_prev_restart_found, &
         restart_transaction_validated, restart_feature_manifest_present, &
         restart_bifurcation_enabled, restart_levee_enabled, volwater_ucat_valid /)
      held_restart_file = gridriver_restart_file

   END SUBROUTINE hold_GridRiverLakeTimeVars_lulcc

   SUBROUTINE restore_GridRiverLakeTimeVars_lulcc ()

   USE MOD_SPMD_Task, only: CoLM_stop
   IMPLICIT NONE

      ! The held state is indexed by local unit catchment and reservoir.  That holds
      ! because the rebuilt network distributes them as before (the distribution
      ! depends on the river network only, not on the land patches); stop if the
      ! rebuilt arrays nevertheless differ in size.
      CALL check_same_size (wdsrf_ucat, held_wdsrf_ucat, 'unit catchments')
      CALL check_same_size (volresv,    held_volresv,    'reservoirs')

      IF (allocated (wdsrf_ucat     )) deallocate (wdsrf_ucat     )
      IF (allocated (wdsrf_ucat_prev)) deallocate (wdsrf_ucat_prev)
      IF (allocated (veloc_riv      )) deallocate (veloc_riv      )
      IF (allocated (momen_riv      )) deallocate (momen_riv      )
      IF (allocated (volresv        )) deallocate (volresv        )
      IF (allocated (volwater_ucat  )) deallocate (volwater_ucat  )
      IF (allocated (acc_rnof_uc    )) deallocate (acc_rnof_uc    )
      IF (allocated (qresv_ups_prev )) deallocate (qresv_ups_prev )
      CALL move_alloc (held_wdsrf_ucat     , wdsrf_ucat     )
      CALL move_alloc (held_wdsrf_ucat_prev, wdsrf_ucat_prev)
      CALL move_alloc (held_veloc_riv      , veloc_riv      )
      CALL move_alloc (held_momen_riv      , momen_riv      )
      CALL move_alloc (held_volresv        , volresv        )
      CALL move_alloc (held_volwater_ucat  , volwater_ucat  )
      CALL move_alloc (held_acc_rnof_uc    , acc_rnof_uc    )
      CALL move_alloc (held_qresv_ups_prev , qresv_ups_prev )
      acctime_rnof                     = held_acctime_rnof
      wdsrf_ucat_prev_valid            = held_flags(1)
      wdsrf_ucat_prev_restart_found    = held_flags(2)
      restart_transaction_validated    = held_flags(3)
      restart_feature_manifest_present = held_flags(4)
      restart_bifurcation_enabled      = held_flags(5)
      restart_levee_enabled            = held_flags(6)
      volwater_ucat_valid              = held_flags(7)
      gridriver_restart_file           = held_restart_file


   CONTAINS

      SUBROUTINE check_same_size (rebuilt, held, what)
      real(r8), allocatable, intent(in) :: rebuilt(:), held(:)
      character(len=*), intent(in) :: what
         IF (allocated(rebuilt) .and. allocated(held)) THEN
            IF (size(rebuilt) /= size(held)) &
               CALL CoLM_stop ('LULCC: the rebuilt river network has a different number of ' // what)
         ENDIF
      END SUBROUTINE check_same_size

   END SUBROUTINE restore_GridRiverLakeTimeVars_lulcc

END MODULE MOD_Grid_RiverLakeTimeVars
#endif
