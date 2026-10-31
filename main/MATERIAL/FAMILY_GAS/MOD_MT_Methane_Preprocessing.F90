#include <define.h>

#if (defined MaterialTransport) && (defined BGC)
MODULE MOD_MT_Methane_Preprocessing

   USE MOD_SPMD_Task, only: CoLM_stop
   USE MOD_MT_Defs, only: material_defs_init, material_index_for_name, material_param_file_for_index, &
      material_is_gas
   USE MOD_MT_Methane_Const, only: DEF_METHANE, read_methane_namelist
   USE MOD_MT_Methane_Registry, only: igas_ch4, methane_is_active

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: methane_preprocessing_requirements

CONTAINS

   SUBROUTINE methane_preprocessing_requirements (requires_lake_soilc, requires_spatial_ph)
      IMPLICIT NONE
      logical, intent(out) :: requires_lake_soilc, requires_spatial_ph

      character(len=512) :: file_param
      logical :: found

      requires_lake_soilc = .false.
      requires_spatial_ph = .false.

      CALL material_defs_init ()
      igas_ch4 = material_index_for_name('CH4', 'METHANE')
      IF (.not. methane_is_active()) RETURN
      IF (.not. material_is_gas(igas_ch4)) THEN
         CALL CoLM_stop (' ***** ERROR: CH4/METHANE preprocessing descriptor must use family=gas')
      ENDIF

      ! Use the runtime parser so keyed aliases, positional files, separators,
      ! and null handling cannot drift between mksrfdata and colm.x.
      CALL material_param_file_for_index (igas_ch4, 'CH4,METHANE', file_param, found)
      IF (.not. found) THEN
         CALL CoLM_stop (' ***** ERROR: CH4 requires DEF_METHANE_PARAM_FILE to include a CH4 parameter file')
      ENDIF
      CALL read_methane_namelist (file_param)

      requires_lake_soilc = DEF_METHANE%allowlakeprod
      requires_spatial_ph = DEF_METHANE%use_spatial_ph

   END SUBROUTINE methane_preprocessing_requirements

END MODULE MOD_MT_Methane_Preprocessing
#endif
