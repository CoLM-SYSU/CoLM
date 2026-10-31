#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_Runtime
   USE MOD_Vars_Global, only: maxsnl, nl_soil
   USE MOD_LandPatch, only: numpatch
   USE MOD_Vars_TimeVariables, only: ldew_rain, ldew_snow, wliq_soisno, wice_soisno, &
      wa, wdsrf, wetwat, scv, waterstorage
   USE MOD_Forcing, only: gforc
   USE MOD_MT_LandPhase, only: mt_land_init, mt_land_final
   USE MOD_MT_Forcing, only: mt_forcing_init, mt_read_forcing, mt_forcing_reset, &
      mt_forcing_final, mt_forcing_read_restart, mt_forcing_lulcc_save, mt_forcing_lulcc_restore
   USE MOD_MT_Vars, only: mt_reset_rnof_step
   USE MOD_MT_Lifecycle, only: mt_lifecycle_reset
   USE MOD_MT_Defs, only: mt_defs_final
   USE MOD_MT_HostFixes, only: mt_bgc_data_init_year
   USE MOD_HistAccRestart, only: read_history_acc_restart, save_history_acc_restart, hist_acc_window_restored
#ifdef GridRiverLakeFlow
   USE MOD_MT_RiverLake, only: mt_riverlake_drop_hist_window
   USE MOD_MT_Sediment, only: mt_sediment_drop_hist_window
#endif
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: mt_runtime_land_init, mt_runtime_forcing_init, mt_runtime_history_restore
   PUBLIC :: mt_runtime_lulcc_forcing_save, mt_runtime_lulcc_forcing_restore
   PUBLIC :: mt_runtime_land_final, mt_runtime_registry_final
   PUBLIC :: mt_read_forcing, mt_forcing_reset, mt_reset_rnof_step
   PUBLIC :: mt_bgc_data_init_year, save_history_acc_restart

   logical :: material_loaded_restart
   character(len=256) :: material_restart_file

CONTAINS

   SUBROUTINE mt_runtime_land_init (s_month, lc_year, jdate, casename, dir_restart, dir_landdata)
      integer, intent(in) :: s_month, lc_year, jdate(3)
      character(len=*), intent(in) :: casename, dir_restart, dir_landdata

      CALL mt_land_init (numpatch, maxsnl, nl_soil, s_month, lc_year, jdate, &
         casename, dir_restart, dir_landdata, ldew_rain, ldew_snow, wliq_soisno, &
         wice_soisno, wa, wdsrf, wetwat, scv, waterstorage, &
         loaded_restart=material_loaded_restart, restart_file=material_restart_file)
   END SUBROUTINE mt_runtime_land_init

   SUBROUTINE mt_runtime_forcing_init ()
      CALL mt_forcing_init (gforc, numpatch)
      CALL mt_forcing_read_restart (material_restart_file, material_loaded_restart)
   END SUBROUTINE mt_runtime_forcing_init

   SUBROUTINE mt_runtime_history_restore (jdate)
      integer, intent(in) :: jdate(3)

      CALL read_history_acc_restart (jdate)
      ! Material rows must follow the host's restored partial-window decision.
      IF (.not. hist_acc_window_restored) THEN
#ifdef GridRiverLakeFlow
         CALL mt_riverlake_drop_hist_window ()
         CALL mt_sediment_drop_hist_window ()
#endif
      ENDIF
   END SUBROUTINE mt_runtime_history_restore

   SUBROUTINE mt_runtime_lulcc_forcing_save ()
      CALL mt_forcing_lulcc_save ()
      CALL mt_forcing_final ()
   END SUBROUTINE mt_runtime_lulcc_forcing_save

   SUBROUTINE mt_runtime_lulcc_forcing_restore ()
      CALL mt_forcing_init (gforc, numpatch)
      CALL mt_forcing_lulcc_restore ()
   END SUBROUTINE mt_runtime_lulcc_forcing_restore

   SUBROUTINE mt_runtime_land_final ()
      CALL mt_forcing_final ()
      CALL mt_land_final ()
   END SUBROUTINE mt_runtime_land_final

   SUBROUTINE mt_runtime_registry_final ()
      CALL mt_lifecycle_reset ()
      CALL mt_defs_final ()
   END SUBROUTINE mt_runtime_registry_final

END MODULE MOD_MT_Runtime
#endif
