#include <define.h>

MODULE MOD_HistAccRestart
!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!  Restart of the history accumulators: the a_* sums in MOD_Vars_1DAccFluxes
!  and, with GridRiverLakeFlow, the river-lake sums in
!  MOD_Grid_RiverLakeHistState.
!
!  Without this module a restart inside a history window (for example a
!  daily restart while history is written monthly) starts the window again
!  from zero, and the first history record after the restart averages only
!  the part of the window after the restart. Here the raw sums are saved
!  next to the model restart and read back, so a split run writes the same
!  history as a continuous run. The model state and the model restart files
!  are not touched.
!
!  Files, in the restart directory of the date:
!     <case>_restart_hist_<date>.nc        land sums, vector blocks
!     <case>_restart_hist_<date>.nc.river  river sums, one serial file
!
!  Save (save_history_acc_restart, called after WRITE_TimeVariables):
!     the sidecar is written with hacc_complete = 0, then hacc_complete is
!     set to 1. A sidecar with hacc_complete /= 1 was left by a run that
!     stopped while saving and is rejected. A restart directory without a
!     sidecar (older runs) starts a new history window, as before.
!
!  At the last time step hist_out writes the partial window and clears the
!  sums before the restart is saved, so hist_acc_final_snapshot keeps the raw
!  sums first.
!
!  The field list between the two marker lines in hacc_fields is checked
!  against allocate_acc_fluxes by .github/workflows/hist_acc_restart_fields.py.
!
!  Not covered: accumulators of an external lake model (EXTERNAL_LAKE).
!
! !HISTORY:
!  2026-10: ported from the CoLM202X land-history restart sidecar.
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Vars_1DAccFluxes
   USE MOD_MT_HistoryRestart, only: hacc_read_name => mt_history_read_name
#ifdef MaterialTransport
   USE MOD_MT_HistoryRestart, only: mt_history_restart_fields, mt_history_frcsat_count
#endif
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: read_history_acc_restart
   PUBLIC :: save_history_acc_restart
   ! .true. when the reader actually restored a PARTIAL window; the material
   ! river history follows the same decision so the two never disagree.
   PUBLIC :: hist_acc_window_restored
   PUBLIC :: hist_acc_final_snapshot

   integer, parameter :: HACC_WRITE = 1
   integer, parameter :: HACC_READ  = 2

   ! date of a sidecar already written by hist_acc_final_snapshot
   integer :: snapshot_date(3) = -1

   ! Set by read_history_acc_restart: .true. only when a partial window was
   ! restored.  A complete window is discarded, and then every material
   ! history window must start at the restart as well.
   logical :: hist_acc_window_restored = .false.

CONTAINS

   !-----------------------------------------------------------------------
   logical FUNCTION hist_acc_final_snapshot (idate, deltim, itstamp, etstamp)

   ! Called by hist_out after the sums of this step are accumulated. At the
   ! last step of the run, when the history window is not complete, hist_out
   ! writes the partial window and clears the sums; keep the raw sums for
   ! the restart that follows.

   USE MOD_TimeManager
   USE MOD_Namelist, only: DEF_HIST_FREQ
   IMPLICIT NONE

   integer,         intent(in) :: idate(3)
   real(r8),        intent(in) :: deltim
   type(timestamp), intent(in) :: itstamp, etstamp

   integer :: jdate(3)
   logical :: window_end

      hist_acc_final_snapshot = .false.
      IF (itstamp < etstamp) RETURN

      SELECT CASE (trim(adjustl(DEF_HIST_FREQ)))
      CASE ('HOURLY')
         window_end = isendofhour (idate, deltim)
      CASE ('DAILY')
         window_end = isendofday  (idate, deltim)
      CASE ('MONTHLY')
         window_end = isendofmonth(idate, deltim)
      CASE ('YEARLY')
         window_end = isendofyear (idate, deltim)
      CASE DEFAULT
         ! TIMESTEP: every step closes a window; unknown: never written
         window_end = .true.
      END SELECT

      IF (window_end) RETURN

      jdate = idate
      CALL adj2begin (jdate)
      CALL write_sidecar (jdate)
      snapshot_date = jdate
      hist_acc_final_snapshot = .true.

   END FUNCTION hist_acc_final_snapshot

   !-----------------------------------------------------------------------
   SUBROUTINE save_history_acc_restart (jdate)

   ! Called right after WRITE_TimeVariables with the same date.

   USE MOD_LandPatch, only: landpatch
   USE MOD_Namelist, only: DEF_simulation_time
   IMPLICIT NONE

   integer, intent(in) :: jdate(3)

   character(len=256) :: file
   real(r8) :: window_seconds, window_nac

      ! Nothing to carry over when the window is complete at this restart: the
      ! reader would discard it (see read_history_acc_restart), and the sidecar
      ! holds every accumulator at patch resolution, so writing it would cost
      ! hundreds of MB per restart for nothing.
      ! nac is maintained on workers; all ranks must take the same return.
      window_nac = 0._r8
      IF (p_is_worker) window_nac = nac
#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, window_nac, 1, MPI_REAL8, MPI_MAX, p_comm_glb, p_err)
#endif
      window_seconds = hacc_window_seconds (jdate)
      IF (window_seconds > 0._r8 .and. window_nac > 0._r8 .and. &
          window_nac * DEF_simulation_time%timestep >= window_seconds - &
          0.5_r8 * DEF_simulation_time%timestep) THEN
         snapshot_date = -1
         RETURN
      ENDIF

      IF (any(snapshot_date /= jdate)) CALL write_sidecar (jdate)
      snapshot_date = -1

      CALL sidecar_filename (jdate, file)
      CALL write_marker (file, 'hacc_complete', 1._r8)

   END SUBROUTINE save_history_acc_restart

   !-----------------------------------------------------------------------
   real(r8) FUNCTION hacc_window_seconds (jdate)

   ! Length in seconds of the history window that ended at jdate, or -1 when
   ! the frequency does not pin it down (TIMESTEP is handled by the caller's
   ! tolerance, MONTHLY/YEARLY only when the restart sits on the boundary).

   USE MOD_Namelist, only: DEF_HIST_FREQ, DEF_simulation_time
   USE MOD_TimeManager, only: julian2monthday
   IMPLICIT NONE

   integer, intent(in) :: jdate(3)

   integer :: month, mday, year, mprev, yprev
   integer :: mdays(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)

      hacc_window_seconds = -1._r8
      SELECT CASE (trim(adjustl(DEF_HIST_FREQ)))
      CASE ('TIMESTEP')
         hacc_window_seconds = DEF_simulation_time%timestep
      CASE ('HOURLY')
         hacc_window_seconds = 3600._r8
      CASE ('DAILY')
         hacc_window_seconds = 86400._r8
      CASE ('MONTHLY')
         CALL julian2monthday (jdate(1), jdate(2), month, mday)
         IF (mday == 1) THEN
            year = jdate(1)
            mprev = month - 1
            IF (mprev == 0) THEN
               mprev = 12
               year = year - 1
            ENDIF
            IF (mprev == 2) THEN
               mdays(2) = 28
               IF (mod(year,4) == 0 .and. (mod(year,100) /= 0 .or. mod(year,400) == 0)) mdays(2) = 29
            ENDIF
            hacc_window_seconds = real(mdays(mprev), r8) * 86400._r8
         ENDIF
      CASE ('YEARLY')
         IF (jdate(2) == 1) THEN
            yprev = jdate(1) - 1
            hacc_window_seconds = 365._r8 * 86400._r8
            IF (mod(yprev,4) == 0 .and. (mod(yprev,100) /= 0 .or. mod(yprev,400) == 0)) &
               hacc_window_seconds = 366._r8 * 86400._r8
         ENDIF
      END SELECT

   END FUNCTION hacc_window_seconds

   !-----------------------------------------------------------------------
   SUBROUTINE read_history_acc_restart (jdate)

   ! Called once at the start of a run, after hist_init and the river-lake
   ! initialization, with the date of the restart that was read.

   USE MOD_Namelist, only: DEF_HIST_FREQ, DEF_simulation_time
   USE MOD_LandPatch, only: landpatch, numpatch
   IMPLICIT NONE

   integer, intent(in) :: jdate(3)

   character(len=256) :: file
   real(r8) :: vmin, vmax, window_seconds

      hist_acc_window_restored = .false.
      CALL sidecar_filename (jdate, file)

      IF (.not. any_block_exists (file)) THEN
         IF (p_is_master) write(*,'(A)') &
            'Note: no history-accumulator restart; the history window starts at the model start.'
         RETURN
      ENDIF

      IF (.not. marker_present (file, 'hacc_complete', landpatch, 'patch')) THEN
         ! A sidecar left by a run killed between creating the blocks and
         ! writing the completion marker carries no trustworthy window.
         ! Ignore it and start a fresh history window instead of refusing to
         ! restart the model (the split run would otherwise be unusable).
         IF (p_is_master) write(*,'(A)') 'Warning: incomplete history-accumulator restart ' &
            //trim(file)//' (no completion marker); ignoring it and restarting the history window.'
         RETURN
      ENDIF
      CALL read_marker (file, 'hacc_complete', vmin, vmax)
      IF (vmin /= 1._r8 .or. vmax /= 1._r8) THEN
         IF (p_is_master) write(*,'(A)') 'Warning: incomplete history-accumulator restart ' &
            //trim(file)//' (hacc_complete is not 1); ignoring it and restarting the history window.'
         RETURN
      ENDIF

      CALL require_var (file, 'hacc_freq', landpatch, 'patch', 1, 0, 0)
      CALL read_marker (file, 'hacc_freq', vmin, vmax)
      IF (vmin /= vmax .or. vmin /= real(freq_code(DEF_HIST_FREQ), r8)) THEN
         IF (p_is_master) write(*,'(A)') &
            'Warning: DEF_HIST_FREQ changed across the restart; the history window starts at the model start.'
         RETURN
      ENDIF

      CALL require_var (file, 'hacc_nac', landpatch, 'patch', 1, 0, 0)
      CALL read_marker (file, 'hacc_nac', vmin, vmax)
      IF (vmin /= vmax .or. vmin < 0._r8 .or. vmin /= anint(vmin)) THEN
         IF (p_is_master) write(*,'(A)') 'ERROR: invalid sample count in '//trim(file)
         CALL CoLM_stop ()
      ENDIF

      ! A window that is already complete at the restart point has been, or
      ! will be, written out by the run that produced this restart: a daily
      ! restart at midnight always carries a full day.  Restoring it would make
      ! the first window after the restart cover two windows (a split run then
      ! reports mean(window1, window2) instead of window2).  Only a genuinely
      ! partial window -- a restart in the middle of one -- is restored.
      window_seconds = hacc_window_seconds (jdate)
      IF (window_seconds > 0._r8 .and. vmin > 0._r8 .and. &
          vmin * DEF_simulation_time%timestep >= window_seconds - &
          0.5_r8 * DEF_simulation_time%timestep) THEN
         IF (p_is_master) write(*,'(A)') &
            'Note: the restart closes a complete history window; the history window starts at the restart.'
         RETURN
      ENDIF

      IF (vmin > 0._r8) THEN
         IF (p_is_worker) nac = vmin
         CALL hacc_fields (file, HACC_READ)
         hist_acc_window_restored = .true.
#ifdef EXTERNAL_LAKE
         IF (p_is_master) write(*,'(A)') &
            'Warning: external-lake history accumulators are not restored; their first window is partial.'
#endif
      ENDIF

#ifdef GridRiverLakeFlow
      CALL require_var (file, 'hacc_river', landpatch, 'patch', 1, 0, 0)
      CALL read_marker (file, 'hacc_river', vmin, vmax)
      IF (vmax > 0._r8) THEN
         CALL river_sums (trim(file)//'.river', HACC_READ)
         hist_acc_window_restored = .true.
      ENDIF
#endif

      ! Every rank must reach the same verdict: the material river history
      ! drops its restored window when this is .false.
#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, hist_acc_window_restored, 1, MPI_LOGICAL, &
         MPI_LAND, p_comm_glb, p_err)
#endif

   END SUBROUTINE read_history_acc_restart

   !-----------------------------------------------------------------------
   SUBROUTINE write_sidecar (jdate)

   USE MOD_Namelist, only: DEF_HIST_FREQ, DEF_dir_restart
   USE MOD_LandPatch, only: landpatch
   USE MOD_LandUrban, only: landurban
   USE MOD_NetCDFVector, only: ncio_create_file_vector, ncio_define_dimension_vector
#ifdef GridRiverLakeFlow
   USE MOD_Grid_RiverLakeHistState, only: acctime_ucat
#endif
   IMPLICIT NONE

   integer, intent(in) :: jdate(3)

   character(len=256) :: file
   character(len=14)  :: cdate
   integer :: active

      CALL sidecar_filename (jdate, file)

      write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') jdate
      IF (p_is_master) CALL system ('mkdir -p '//trim(DEF_dir_restart)//'/'//trim(cdate))
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      CALL ncio_create_file_vector (file, landpatch)
      CALL ncio_define_dimension_vector (file, landpatch, 'patch')
#ifdef URBAN_MODEL
      CALL ncio_define_dimension_vector (file, landurban, 'urban')
#endif

      CALL write_marker (file, 'hacc_complete', 0._r8)
      CALL write_marker (file, 'hacc_freq', real(freq_code(DEF_HIST_FREQ), r8))
      CALL write_marker (file, 'hacc_nac', nac)

      active = 0
      IF (p_is_worker) THEN
         IF (nac > 0._r8) active = 1
      ENDIF
#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, active, 1, MPI_INTEGER, MPI_MAX, p_comm_glb, p_err)
#endif
      IF (active > 0) CALL hacc_fields (file, HACC_WRITE)

#ifdef GridRiverLakeFlow
      active = 0
      IF (p_is_worker .and. allocated(acctime_ucat)) THEN
         IF (any(acctime_ucat > 0._r8)) active = 1
      ENDIF
#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, active, 1, MPI_INTEGER, MPI_MAX, p_comm_glb, p_err)
#endif
      CALL write_marker (file, 'hacc_river', real(active, r8))
      IF (active > 0) CALL river_sums (trim(file)//'.river', HACC_WRITE)
#endif

   END SUBROUTINE write_sidecar

   !-----------------------------------------------------------------------
   SUBROUTINE hacc_fields (file, mode)

   USE MOD_LandPatch, only: landpatch
   USE MOD_LandUrban, only: landurban
   IMPLICIT NONE

   character(len=*), intent(in) :: file
   integer,          intent(in) :: mode

      ! Generated by .github/workflows/hist_acc_restart_fields.py --write;
      ! keep the marker lines.
      ! >>> fields from allocate_acc_fluxes
      CALL hacc_field1d (file, mode, 'a_us', a_us, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_vs', a_vs, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_t', a_t, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_q', a_q, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_prc', a_prc, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_prl', a_prl, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_pbot', a_pbot, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_frl', a_frl, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_solarin', a_solarin, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_hpbl', a_hpbl, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_taux', a_taux, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_tauy', a_tauy, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fsena', a_fsena, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_lfevpa', a_lfevpa, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fevpa', a_fevpa, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fsenl', a_fsenl, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fevpl', a_fevpl, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_etr', a_etr, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fseng', a_fseng, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fevpg', a_fevpg, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fgrnd', a_fgrnd, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_sabvsun', a_sabvsun, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_sabvsha', a_sabvsha, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_sabg', a_sabg, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_olrg', a_olrg, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_rnet', a_rnet, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_xerr', a_xerr, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_zerr', a_zerr, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_rsur', a_rsur, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_rsur_se', a_rsur_se, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_rsur_ie', a_rsur_ie, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_rsub', a_rsub, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_rnof', a_rnof, landpatch, 'patch')
#ifdef CatchLateralFlow
      CALL hacc_field1d (file, mode, 'a_xwsur', a_xwsur, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_xwsub', a_xwsub, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fldarea', a_fldarea, landpatch, 'patch')
#endif
      CALL hacc_field1d (file, mode, 'a_qintr', a_qintr, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_qinfl', a_qinfl, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_qdrip', a_qdrip, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_rstfacsun', a_rstfacsun, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_rstfacsha', a_rstfacsha, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_gssun', a_gssun, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_gssha', a_gssha, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_rss', a_rss, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_wdsrf', a_wdsrf, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_zwt', a_zwt, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_wa', a_wa, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_wat', a_wat, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_wetwat', a_wetwat, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_assim', a_assim, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_respc', a_respc, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_assimsun', a_assimsun, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_assimsha', a_assimsha, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_etrsun', a_etrsun, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_etrsha', a_etrsha, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_qcharge', a_qcharge, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_t_grnd', a_t_grnd, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_tleaf', a_tleaf, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_ldew_rain', a_ldew_rain, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_ldew_snow', a_ldew_snow, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_ldew', a_ldew, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_scv', a_scv, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_snowdp', a_snowdp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fsno', a_fsno, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_frcsat', a_frcsat, landpatch, 'patch')
#ifdef MaterialTransport
      CALL mt_history_frcsat_count(file, mode, HACC_WRITE, nac, a_frcsat, &
         hacc_field1d, write_marker, read_marker, require_var, global_max)
#endif
      CALL hacc_field1d (file, mode, 'a_sigf', a_sigf, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_green', a_green, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_lai', a_lai, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_laisun', a_laisun, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_laisha', a_laisha, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_sai', a_sai, landpatch, 'patch')
      CALL hacc_field3d (file, mode, 'a_alb', a_alb, landpatch, 'patch')
#ifdef HYPERSPECTRAL
      CALL hacc_field3d (file, mode, 'a_alb_hires', a_alb_hires, landpatch, 'patch')
      CALL hacc_field3d (file, mode, 'a_reflectance_out', a_reflectance_out, landpatch, 'patch')
      CALL hacc_field3d (file, mode, 'a_transmittance_out', a_transmittance_out, landpatch, 'patch')
#endif
      CALL hacc_field1d (file, mode, 'a_emis', a_emis, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_z0m', a_z0m, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_trad', a_trad, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_tref', a_tref, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_t2m_wmo', a_t2m_wmo, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_qref', a_qref, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_rain', a_rain, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_snow', a_snow, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_o3uptakesun', a_o3uptakesun, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_o3uptakesha', a_o3uptakesha, landpatch, 'patch')
#ifdef DataAssimilation
      CALL hacc_field3d (file, mode, 'a_h2osoi_ens', a_h2osoi_ens, landpatch, 'patch')
      CALL hacc_field3d (file, mode, 'a_t_brt_smap_ens', a_t_brt_smap_ens, landpatch, 'patch')
      CALL hacc_field3d (file, mode, 'a_t_brt_fy3d_ens', a_t_brt_fy3d_ens, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_t_brt_smap', a_t_brt_smap, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_t_brt_fy3d', a_t_brt_fy3d, landpatch, 'patch')
      CALL hacc_field3d (file, mode, 'a_wliq_soisno_ens', a_wliq_soisno_ens, landpatch, 'patch')
      CALL hacc_field3d (file, mode, 'a_wice_soisno_ens', a_wice_soisno_ens, landpatch, 'patch')
      CALL hacc_field3d (file, mode, 'a_t_soisno_ens', a_t_soisno_ens, landpatch, 'patch')
#endif
#ifdef URBAN_MODEL
      CALL hacc_field1d (file, mode, 'a_t_room', a_t_room, landurban, 'urban')
      CALL hacc_field1d (file, mode, 'a_tafu', a_tafu, landurban, 'urban')
      CALL hacc_field1d (file, mode, 'a_fhac', a_fhac, landurban, 'urban')
      CALL hacc_field1d (file, mode, 'a_fwst', a_fwst, landurban, 'urban')
      CALL hacc_field1d (file, mode, 'a_fach', a_fach, landurban, 'urban')
      CALL hacc_field1d (file, mode, 'a_fahe', a_fahe, landurban, 'urban')
      CALL hacc_field1d (file, mode, 'a_fhah', a_fhah, landurban, 'urban')
      CALL hacc_field1d (file, mode, 'a_vehc', a_vehc, landurban, 'urban')
      CALL hacc_field1d (file, mode, 'a_meta', a_meta, landurban, 'urban')
      CALL hacc_field1d (file, mode, 'a_senroof', a_senroof, landurban, 'urban')
      CALL hacc_field1d (file, mode, 'a_senwsun', a_senwsun, landurban, 'urban')
      CALL hacc_field1d (file, mode, 'a_senwsha', a_senwsha, landurban, 'urban')
      CALL hacc_field1d (file, mode, 'a_sengimp', a_sengimp, landurban, 'urban')
      CALL hacc_field1d (file, mode, 'a_sengper', a_sengper, landurban, 'urban')
      CALL hacc_field1d (file, mode, 'a_senurbl', a_senurbl, landurban, 'urban')
      CALL hacc_field1d (file, mode, 'a_lfevproof', a_lfevproof, landurban, 'urban')
      CALL hacc_field1d (file, mode, 'a_lfevpgimp', a_lfevpgimp, landurban, 'urban')
      CALL hacc_field1d (file, mode, 'a_lfevpgper', a_lfevpgper, landurban, 'urban')
      CALL hacc_field1d (file, mode, 'a_lfevpurbl', a_lfevpurbl, landurban, 'urban')
      CALL hacc_field1d (file, mode, 'a_troof', a_troof, landurban, 'urban')
      CALL hacc_field1d (file, mode, 'a_twall', a_twall, landurban, 'urban')
#endif
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      CALL hacc_field1d (file, mode, 'a_lai_enftemp', a_lai_enftemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_lai_enfboreal', a_lai_enfboreal, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_lai_dnfboreal', a_lai_dnfboreal, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_lai_ebftrop', a_lai_ebftrop, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_lai_ebftemp', a_lai_ebftemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_lai_dbftrop', a_lai_dbftrop, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_lai_dbftemp', a_lai_dbftemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_lai_dbfboreal', a_lai_dbfboreal, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_lai_ebstemp', a_lai_ebstemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_lai_dbstemp', a_lai_dbstemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_lai_dbsboreal', a_lai_dbsboreal, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_lai_c3arcgrass', a_lai_c3arcgrass, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_lai_c3grass', a_lai_c3grass, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_lai_c4grass', a_lai_c4grass, landpatch, 'patch')
#endif
#ifdef BGC
      CALL hacc_field1d (file, mode, 'a_farea_burned', a_farea_burned, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_baf_crop', a_baf_crop, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_baf_peatf', a_baf_peatf, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_nfire', a_nfire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fuelc', a_fuelc, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_btran2', a_btran2, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_pft_fire_closs', a_pft_fire_closs, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_pft_fire_nloss', a_pft_fire_nloss, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_col_fire_closs', a_col_fire_closs, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_col_fire_nloss', a_col_fire_nloss, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_somc_fire', a_somc_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_litfire', a_litfire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_somfire', a_somfire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_totfire', a_totfire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_leafc_to_fire', a_m_leafc_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_frootc_to_fire', a_m_frootc_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livestemc_to_fire', a_m_livestemc_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadstemc_to_fire', a_m_deadstemc_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livecrootc_to_fire', a_m_livecrootc_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadcrootc_to_fire', a_m_deadcrootc_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_leafc_storage_to_fire', a_m_leafc_storage_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_frootc_storage_to_fire', a_m_frootc_storage_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livestemc_storage_to_fire', a_m_livestemc_storage_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadstemc_storage_to_fire', a_m_deadstemc_storage_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livecrootc_storage_to_fire', a_m_livecrootc_storage_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadcrootc_storage_to_fire', a_m_deadcrootc_storage_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_gresp_storage_to_fire', a_m_gresp_storage_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_leafc_xfer_to_fire', a_m_leafc_xfer_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_frootc_xfer_to_fire', a_m_frootc_xfer_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livestemc_xfer_to_fire', a_m_livestemc_xfer_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadstemc_xfer_to_fire', a_m_deadstemc_xfer_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livecrootc_xfer_to_fire', a_m_livecrootc_xfer_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadcrootc_xfer_to_fire', a_m_deadcrootc_xfer_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_gresp_xfer_to_fire', a_m_gresp_xfer_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livestemc_to_deadstemc_fire', a_m_livestemc_to_deadstemc_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livecrootc_to_deadcrootc_fire', a_m_livecrootc_to_deadcrootc_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_leafc_to_litter_fire', a_m_leafc_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_frootc_to_litter_fire', a_m_frootc_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livestemc_to_litter_fire', a_m_livestemc_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadstemc_to_litter_fire', a_m_deadstemc_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livecrootc_to_litter_fire', a_m_livecrootc_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadcrootc_to_litter_fire', a_m_deadcrootc_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_leafc_storage_to_litter_fire', a_m_leafc_storage_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_frootc_storage_to_litter_fire', a_m_frootc_storage_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livestemc_storage_to_litter_fire', a_m_livestemc_storage_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadstemc_storage_to_litter_fire', a_m_deadstemc_storage_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livecrootc_storage_to_litter_fire', a_m_livecrootc_storage_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadcrootc_storage_to_litter_fire', a_m_deadcrootc_storage_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_gresp_storage_to_litter_fire', a_m_gresp_storage_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_leafc_xfer_to_litter_fire', a_m_leafc_xfer_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_frootc_xfer_to_litter_fire', a_m_frootc_xfer_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livestemc_xfer_to_litter_fire', a_m_livestemc_xfer_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadstemc_xfer_to_litter_fire', a_m_deadstemc_xfer_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livecrootc_xfer_to_litter_fire', a_m_livecrootc_xfer_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadcrootc_xfer_to_litter_fire', a_m_deadcrootc_xfer_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_gresp_xfer_to_litter_fire', a_m_gresp_xfer_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_leafn_to_fire', a_m_leafn_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_frootn_to_fire', a_m_frootn_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livestemn_to_fire', a_m_livestemn_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadstemn_to_fire', a_m_deadstemn_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livecrootn_to_fire', a_m_livecrootn_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadcrootn_to_fire', a_m_deadcrootn_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_leafn_storage_to_fire', a_m_leafn_storage_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_frootn_storage_to_fire', a_m_frootn_storage_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livestemn_storage_to_fire', a_m_livestemn_storage_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadstemn_storage_to_fire', a_m_deadstemn_storage_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livecrootn_storage_to_fire', a_m_livecrootn_storage_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadcrootn_storage_to_fire', a_m_deadcrootn_storage_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_leafn_xfer_to_fire', a_m_leafn_xfer_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_frootn_xfer_to_fire', a_m_frootn_xfer_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livestemn_xfer_to_fire', a_m_livestemn_xfer_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadstemn_xfer_to_fire', a_m_deadstemn_xfer_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livecrootn_xfer_to_fire', a_m_livecrootn_xfer_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadcrootn_xfer_to_fire', a_m_deadcrootn_xfer_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livestemn_to_deadstemn_fire', a_m_livestemn_to_deadstemn_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livecrootn_to_deadcrootn_fire', a_m_livecrootn_to_deadcrootn_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_retransn_to_fire', a_m_retransn_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_leafn_to_litter_fire', a_m_leafn_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_frootn_to_litter_fire', a_m_frootn_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livestemn_to_litter_fire', a_m_livestemn_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadstemn_to_litter_fire', a_m_deadstemn_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livecrootn_to_litter_fire', a_m_livecrootn_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadcrootn_to_litter_fire', a_m_deadcrootn_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_leafn_storage_to_litter_fire', a_m_leafn_storage_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_frootn_storage_to_litter_fire', a_m_frootn_storage_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livestemn_storage_to_litter_fire', a_m_livestemn_storage_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadstemn_storage_to_litter_fire', a_m_deadstemn_storage_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livecrootn_storage_to_litter_fire', a_m_livecrootn_storage_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadcrootn_storage_to_litter_fire', a_m_deadcrootn_storage_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_leafn_xfer_to_litter_fire', a_m_leafn_xfer_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_frootn_xfer_to_litter_fire', a_m_frootn_xfer_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livestemn_xfer_to_litter_fire', a_m_livestemn_xfer_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadstemn_xfer_to_litter_fire', a_m_deadstemn_xfer_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_livecrootn_xfer_to_litter_fire', a_m_livecrootn_xfer_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_deadcrootn_xfer_to_litter_fire', a_m_deadcrootn_xfer_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_retransn_to_litter_fire', a_m_retransn_to_litter_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_litr1_c_to_fire', a_m_litr1_c_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_litr1_n_to_fire', a_m_litr1_n_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_litr2_c_to_fire', a_m_litr2_c_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_litr2_n_to_fire', a_m_litr2_n_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_litr3_c_to_fire', a_m_litr3_c_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_litr3_n_to_fire', a_m_litr3_n_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_cwd_c_to_fire', a_m_cwd_c_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_m_cwd_n_to_fire', a_m_cwd_n_to_fire, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafc', a_leafc, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafc_storage', a_leafc_storage, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafc_xfer', a_leafc_xfer, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_frootc', a_frootc, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_frootc_storage', a_frootc_storage, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_frootc_xfer', a_frootc_xfer, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livestemc', a_livestemc, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livestemc_storage', a_livestemc_storage, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livestemc_xfer', a_livestemc_xfer, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadstemc', a_deadstemc, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadstemc_storage', a_deadstemc_storage, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadstemc_xfer', a_deadstemc_xfer, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livecrootc', a_livecrootc, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livecrootc_storage', a_livecrootc_storage, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livecrootc_xfer', a_livecrootc_xfer, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadcrootc', a_deadcrootc, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadcrootc_storage', a_deadcrootc_storage, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadcrootc_xfer', a_deadcrootc_xfer, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_grainc', a_grainc, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_grainc_storage', a_grainc_storage, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_grainc_xfer', a_grainc_xfer, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafn', a_leafn, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafn_storage', a_leafn_storage, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafn_xfer', a_leafn_xfer, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_frootn', a_frootn, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_frootn_storage', a_frootn_storage, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_frootn_xfer', a_frootn_xfer, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livestemn', a_livestemn, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livestemn_storage', a_livestemn_storage, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livestemn_xfer', a_livestemn_xfer, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadstemn', a_deadstemn, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadstemn_storage', a_deadstemn_storage, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadstemn_xfer', a_deadstemn_xfer, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livecrootn', a_livecrootn, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livecrootn_storage', a_livecrootn_storage, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livecrootn_xfer', a_livecrootn_xfer, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadcrootn', a_deadcrootn, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadcrootn_storage', a_deadcrootn_storage, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadcrootn_xfer', a_deadcrootn_xfer, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_grainn', a_grainn, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_grainn_storage', a_grainn_storage, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_grainn_xfer', a_grainn_xfer, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_retransn', a_retransn, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_gpp', a_gpp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_downreg', a_downreg, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_ar', a_ar, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_cwdprod', a_cwdprod, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_cwddecomp', a_cwddecomp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_hr', a_hr, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fpg', a_fpg, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fpi', a_fpi, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_totvegc', a_totvegc, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_totlitc', a_totlitc, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_totcwdc', a_totcwdc, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_totsomc', a_totsomc, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_totcolc', a_totcolc, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_totvegn', a_totvegn, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_totlitn', a_totlitn, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_totcwdn', a_totcwdn, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_totsomn', a_totsomn, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_totcoln', a_totcoln, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_gpp_enftemp', a_gpp_enftemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_gpp_enfboreal', a_gpp_enfboreal, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_gpp_dnfboreal', a_gpp_dnfboreal, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_gpp_ebftrop', a_gpp_ebftrop, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_gpp_ebftemp', a_gpp_ebftemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_gpp_dbftrop', a_gpp_dbftrop, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_gpp_dbftemp', a_gpp_dbftemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_gpp_dbfboreal', a_gpp_dbfboreal, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_gpp_ebstemp', a_gpp_ebstemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_gpp_dbstemp', a_gpp_dbstemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_gpp_dbsboreal', a_gpp_dbsboreal, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_gpp_c3arcgrass', a_gpp_c3arcgrass, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_gpp_c3grass', a_gpp_c3grass, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_gpp_c4grass', a_gpp_c4grass, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npp_enftemp', a_npp_enftemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npp_enfboreal', a_npp_enfboreal, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npp_dnfboreal', a_npp_dnfboreal, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npp_ebftrop', a_npp_ebftrop, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npp_ebftemp', a_npp_ebftemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npp_dbftrop', a_npp_dbftrop, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npp_dbftemp', a_npp_dbftemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npp_dbfboreal', a_npp_dbfboreal, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npp_ebstemp', a_npp_ebstemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npp_dbstemp', a_npp_dbstemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npp_dbsboreal', a_npp_dbsboreal, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npp_c3arcgrass', a_npp_c3arcgrass, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npp_c3grass', a_npp_c3grass, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npp_c4grass', a_npp_c4grass, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npptoleafc_enftemp', a_npptoleafc_enftemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npptoleafc_enfboreal', a_npptoleafc_enfboreal, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npptoleafc_dnfboreal', a_npptoleafc_dnfboreal, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npptoleafc_ebftrop', a_npptoleafc_ebftrop, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npptoleafc_ebftemp', a_npptoleafc_ebftemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npptoleafc_dbftrop', a_npptoleafc_dbftrop, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npptoleafc_dbftemp', a_npptoleafc_dbftemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npptoleafc_dbfboreal', a_npptoleafc_dbfboreal, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npptoleafc_ebstemp', a_npptoleafc_ebstemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npptoleafc_dbstemp', a_npptoleafc_dbstemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npptoleafc_dbsboreal', a_npptoleafc_dbsboreal, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npptoleafc_c3arcgrass', a_npptoleafc_c3arcgrass, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npptoleafc_c3grass', a_npptoleafc_c3grass, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_npptoleafc_c4grass', a_npptoleafc_c4grass, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafc_enftemp', a_leafc_enftemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafc_enfboreal', a_leafc_enfboreal, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafc_dnfboreal', a_leafc_dnfboreal, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafc_ebftrop', a_leafc_ebftrop, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafc_ebftemp', a_leafc_ebftemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafc_dbftrop', a_leafc_dbftrop, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafc_dbftemp', a_leafc_dbftemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafc_dbfboreal', a_leafc_dbfboreal, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafc_ebstemp', a_leafc_ebstemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafc_dbstemp', a_leafc_dbstemp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafc_dbsboreal', a_leafc_dbsboreal, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafc_c3arcgrass', a_leafc_c3arcgrass, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafc_c3grass', a_leafc_c3grass, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafc_c4grass', a_leafc_c4grass, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_O2_DECOMP_DEPTH_UNSAT', a_O2_DECOMP_DEPTH_UNSAT, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_CONC_O2_UNSAT', a_CONC_O2_UNSAT, landpatch, 'patch')
#ifdef CROP
      CALL hacc_field1d (file, mode, 'a_pdcorn', a_pdcorn, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_pdswheat', a_pdswheat, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_pdwwheat', a_pdwwheat, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_pdsoybean', a_pdsoybean, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_pdcotton', a_pdcotton, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_pdrice1', a_pdrice1, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_pdrice2', a_pdrice2, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_pdsugarcane', a_pdsugarcane, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_plantdate', a_plantdate, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_manunitro', a_manunitro, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fertnitro_corn', a_fertnitro_corn, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fertnitro_swheat', a_fertnitro_swheat, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fertnitro_wwheat', a_fertnitro_wwheat, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fertnitro_soybean', a_fertnitro_soybean, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fertnitro_cotton', a_fertnitro_cotton, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fertnitro_rice1', a_fertnitro_rice1, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fertnitro_rice2', a_fertnitro_rice2, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fertnitro_sugarcane', a_fertnitro_sugarcane, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_irrig_method_corn', a_irrig_method_corn, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_irrig_method_swheat', a_irrig_method_swheat, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_irrig_method_wwheat', a_irrig_method_wwheat, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_irrig_method_soybean', a_irrig_method_soybean, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_irrig_method_cotton', a_irrig_method_cotton, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_irrig_method_rice1', a_irrig_method_rice1, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_irrig_method_rice2', a_irrig_method_rice2, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_irrig_method_sugarcane', a_irrig_method_sugarcane, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_cphase', a_cphase, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_hui', a_hui, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_gddmaturity', a_gddmaturity, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_gddplant', a_gddplant, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_vf', a_vf, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_cropprod1c', a_cropprod1c, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_cropprod1c_loss', a_cropprod1c_loss, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_cropseedc_deficit', a_cropseedc_deficit, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_grainc_to_cropprodc', a_grainc_to_cropprodc, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_grainc_to_seed', a_grainc_to_seed, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fert_to_sminn', a_fert_to_sminn, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_sum_irrig', a_sum_irrig, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_sum_deficit_irrig', a_sum_deficit_irrig, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_sum_irrig_count', a_sum_irrig_count, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_waterstorage', a_waterstorage, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_groundwater_demand', a_groundwater_demand, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_groundwater_supply', a_groundwater_supply, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_reservoirriver_demand', a_reservoirriver_demand, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_reservoirriver_supply', a_reservoirriver_supply, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_reservoir_supply', a_reservoir_supply, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_river_supply', a_river_supply, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_runoff_supply', a_runoff_supply, landpatch, 'patch')
#endif
      CALL hacc_field1d (file, mode, 'a_ndep_to_sminn', a_ndep_to_sminn, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_abm', a_abm, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_gdp', a_gdp, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_peatf', a_peatf, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_hdm', a_hdm, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_lnfm', a_lnfm, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafcCap', a_leafcCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafc_storageCap', a_leafc_storageCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafc_xferCap', a_leafc_xferCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_frootcCap', a_frootcCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_frootc_storageCap', a_frootc_storageCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_frootc_xferCap', a_frootc_xferCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livestemcCap', a_livestemcCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livestemc_storageCap', a_livestemc_storageCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livestemc_xferCap', a_livestemc_xferCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadstemcCap', a_deadstemcCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadstemc_storageCap', a_deadstemc_storageCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadstemc_xferCap', a_deadstemc_xferCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livecrootcCap', a_livecrootcCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livecrootc_storageCap', a_livecrootc_storageCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livecrootc_xferCap', a_livecrootc_xferCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadcrootcCap', a_deadcrootcCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadcrootc_storageCap', a_deadcrootc_storageCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadcrootc_xferCap', a_deadcrootc_xferCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafnCap', a_leafnCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafn_storageCap', a_leafn_storageCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_leafn_xferCap', a_leafn_xferCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_frootnCap', a_frootnCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_frootn_storageCap', a_frootn_storageCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_frootn_xferCap', a_frootn_xferCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livestemnCap', a_livestemnCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livestemn_storageCap', a_livestemn_storageCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livestemn_xferCap', a_livestemn_xferCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadstemnCap', a_deadstemnCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadstemn_storageCap', a_deadstemn_storageCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadstemn_xferCap', a_deadstemn_xferCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livecrootnCap', a_livecrootnCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livecrootn_storageCap', a_livecrootn_storageCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_livecrootn_xferCap', a_livecrootn_xferCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadcrootnCap', a_deadcrootnCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadcrootn_storageCap', a_deadcrootn_storageCap, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_deadcrootn_xferCap', a_deadcrootn_xferCap, landpatch, 'patch')
#endif
      CALL hacc_field1d (file, mode, 'a_ozone', a_ozone, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_t_soisno', a_t_soisno, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_wliq_soisno', a_wliq_soisno, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_wice_soisno', a_wice_soisno, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_h2osoi', a_h2osoi, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_qlayer', a_qlayer, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_lake_deficit', a_lake_deficit, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_rootr', a_rootr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_BD_all', a_BD_all, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_wfc', a_wfc, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_OM_density', a_OM_density, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_vegwp', a_vegwp, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_dz_lake', a_dz_lake, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_t_lake', a_t_lake, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_lake_icefrac', a_lake_icefrac, landpatch, 'patch')
#ifdef BGC
      CALL hacc_field2d (file, mode, 'a_litr1c_vr', a_litr1c_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_litr2c_vr', a_litr2c_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_litr3c_vr', a_litr3c_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_soil1c_vr', a_soil1c_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_soil2c_vr', a_soil2c_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_soil3c_vr', a_soil3c_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_cwdc_vr', a_cwdc_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_litr1n_vr', a_litr1n_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_litr2n_vr', a_litr2n_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_litr3n_vr', a_litr3n_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_soil1n_vr', a_soil1n_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_soil2n_vr', a_soil2n_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_soil3n_vr', a_soil3n_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_cwdn_vr', a_cwdn_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_totsoiln_vr', a_totsoiln_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_sminn_vr', a_sminn_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_litr1cCap_vr', a_litr1cCap_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_litr2cCap_vr', a_litr2cCap_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_litr3cCap_vr', a_litr3cCap_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_soil1cCap_vr', a_soil1cCap_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_soil2cCap_vr', a_soil2cCap_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_soil3cCap_vr', a_soil3cCap_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_cwdcCap_vr', a_cwdcCap_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_litr1nCap_vr', a_litr1nCap_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_litr2nCap_vr', a_litr2nCap_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_litr3nCap_vr', a_litr3nCap_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_soil1nCap_vr', a_soil1nCap_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_soil2nCap_vr', a_soil2nCap_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_soil3nCap_vr', a_soil3nCap_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_cwdnCap_vr', a_cwdnCap_vr, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_t_scalar', a_t_scalar, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_w_scalar', a_w_scalar, landpatch, 'patch')
#endif
      CALL hacc_field1d (file, mode, 'a_ustar', a_ustar, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_ustar2', a_ustar2, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_tstar', a_tstar, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_qstar', a_qstar, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_zol', a_zol, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_rib', a_rib, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fm', a_fm, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fh', a_fh, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fq', a_fq, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_us10m', a_us10m, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_vs10m', a_vs10m, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_fm10m', a_fm10m, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_sr', a_sr, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_solvd', a_solvd, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_solvi', a_solvi, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_solnd', a_solnd, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_solni', a_solni, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_srvd', a_srvd, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_srvi', a_srvi, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_srnd', a_srnd, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_srni', a_srni, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_solvdln', a_solvdln, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_solviln', a_solviln, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_solndln', a_solndln, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_solniln', a_solniln, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_srvdln', a_srvdln, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_srviln', a_srviln, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_srndln', a_srndln, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'a_srniln', a_srniln, landpatch, 'patch')
#ifdef HYPERSPECTRAL
      CALL hacc_field2d (file, mode, 'a_sol_dir_ln_hires', a_sol_dir_ln_hires, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_sol_dif_ln_hires', a_sol_dif_ln_hires, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_sr_dir_ln_hires', a_sr_dir_ln_hires, landpatch, 'patch')
      CALL hacc_field2d (file, mode, 'a_sr_dif_ln_hires', a_sr_dif_ln_hires, landpatch, 'patch')
#endif
      CALL hacc_field2d (file, mode, 'a_sensors', a_sensors, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'nac_ln', nac_ln, landpatch, 'patch')
      CALL hacc_field1d (file, mode, 'nac_dt', nac_dt, landpatch, 'patch')
      ! <<< fields from allocate_acc_fluxes

#ifdef MaterialTransport
      CALL mt_history_restart_fields(file, mode, HACC_WRITE, hacc_field1d, hacc_field2d, hacc_field3d)
#endif

   END SUBROUTINE hacc_fields

   !-----------------------------------------------------------------------
   SUBROUTINE hacc_field1d (file, mode, name, var, pixelset, vecname)

   USE MOD_Namelist, only: DEF_REST_CompressLevel
   USE MOD_Pixelset, only: pixelset_type
   USE MOD_NetCDFVector, only: ncio_write_vector, ncio_read_vector
   IMPLICIT NONE

   character(len=*),      intent(in)    :: file, name, vecname
   integer,               intent(in)    :: mode
   real(r8), allocatable, intent(inout) :: var (:)
   type(pixelset_type),   intent(in)    :: pixelset

   integer :: n(1)
   real(r8), allocatable :: empty (:)

      n = 0
      IF (allocated(var)) n = 1
      CALL global_max (n)
      IF (n(1) == 0) RETURN

      IF (mode == HACC_WRITE) THEN
         IF (allocated(var)) THEN
            CALL ncio_write_vector (file, name, vecname, pixelset, var, DEF_REST_CompressLevel)
         ELSE
            allocate (empty (0))
            CALL ncio_write_vector (file, name, vecname, pixelset, empty, DEF_REST_CompressLevel)
            deallocate (empty)
         ENDIF
      ELSE
         CALL require_var (file, name, pixelset, vecname, 1, 0, 0)
         CALL ncio_read_vector (file, name, pixelset, var)
      ENDIF

   END SUBROUTINE hacc_field1d

   !-----------------------------------------------------------------------


   !-----------------------------------------------------------------------
   SUBROUTINE hacc_field2d (file, mode, name, var, pixelset, vecname)

   USE MOD_Namelist, only: DEF_REST_CompressLevel
   USE MOD_Pixelset, only: pixelset_type
   USE MOD_NetCDFVector, only: ncio_write_vector, ncio_read_vector, ncio_define_dimension_vector
   IMPLICIT NONE

   character(len=*),      intent(in)    :: file, name, vecname
   integer,               intent(in)    :: mode
   real(r8), allocatable, intent(inout) :: var (:,:)
   type(pixelset_type),   intent(in)    :: pixelset

   character(len=:), allocatable :: read_name
   integer :: n(1)
   real(r8), allocatable :: empty (:,:)

      n = 0
      IF (allocated(var)) n = size(var,1)
      CALL global_max (n)
      IF (n(1) == 0) RETURN

      IF (mode == HACC_WRITE) THEN
         CALL ncio_define_dimension_vector (file, pixelset, 'd1_'//trim(name), n(1))
         IF (allocated(var)) THEN
            CALL ncio_write_vector (file, name, 'd1_'//trim(name), n(1), vecname, pixelset, var, &
               DEF_REST_CompressLevel)
         ELSE
            allocate (empty (n(1),0))
            CALL ncio_write_vector (file, name, 'd1_'//trim(name), n(1), vecname, pixelset, empty, &
               DEF_REST_CompressLevel)
            deallocate (empty)
         ENDIF
      ELSE
         read_name = hacc_read_name (file, name, pixelset)
         CALL require_var (file, read_name, pixelset, vecname, 2, n(1), 0)
         CALL ncio_read_vector (file, read_name, n(1), pixelset, var)
      ENDIF

   END SUBROUTINE hacc_field2d

   !-----------------------------------------------------------------------
   SUBROUTINE hacc_field3d (file, mode, name, var, pixelset, vecname)

   USE MOD_Namelist, only: DEF_REST_CompressLevel
   USE MOD_Pixelset, only: pixelset_type
   USE MOD_NetCDFVector, only: ncio_write_vector, ncio_read_vector, ncio_define_dimension_vector
   IMPLICIT NONE

   character(len=*),      intent(in)    :: file, name, vecname
   integer,               intent(in)    :: mode
   real(r8), allocatable, intent(inout) :: var (:,:,:)
   type(pixelset_type),   intent(in)    :: pixelset

   character(len=:), allocatable :: read_name
   integer :: n(2)
   real(r8), allocatable :: empty (:,:,:)

      n = 0
      IF (allocated(var)) n = [size(var,1), size(var,2)]
      CALL global_max (n)
      IF (n(1) == 0 .or. n(2) == 0) RETURN

      IF (mode == HACC_WRITE) THEN
         CALL ncio_define_dimension_vector (file, pixelset, 'd1_'//trim(name), n(1))
         CALL ncio_define_dimension_vector (file, pixelset, 'd2_'//trim(name), n(2))
         IF (allocated(var)) THEN
            CALL ncio_write_vector (file, name, 'd1_'//trim(name), n(1), 'd2_'//trim(name), n(2), &
               vecname, pixelset, var, DEF_REST_CompressLevel)
         ELSE
            allocate (empty (n(1),n(2),0))
            CALL ncio_write_vector (file, name, 'd1_'//trim(name), n(1), 'd2_'//trim(name), n(2), &
               vecname, pixelset, empty, DEF_REST_CompressLevel)
            deallocate (empty)
         ENDIF
      ELSE
         read_name = hacc_read_name (file, name, pixelset)
         CALL require_var (file, read_name, pixelset, vecname, 3, n(1), n(2))
         CALL ncio_read_vector (file, read_name, n(1), n(2), pixelset, var)
      ENDIF

   END SUBROUTINE hacc_field3d

#ifdef GridRiverLakeFlow
   !-----------------------------------------------------------------------
   SUBROUTINE river_sums (file, mode)

   ! River-lake sums on unit catchments, reservoirs and bifurcation pathways,
   ! in one serial file written by the master. Pathway columns follow the
   ! global pathway order of the bifurcation restart, whose path signature
   ! has been checked in grid_riverlake_flow_init before this is read.

   USE MOD_Namelist, only: DEF_USE_LEVEE, DEF_USE_BIFURCATION, DEF_Reservoir_Method, &
      DEF_REST_CompressLevel
   USE MOD_NetCDFSerial, only: ncio_create_file, ncio_define_dimension, ncio_write_serial
   USE MOD_Grid_RiverLakeHistState
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_data_address, &
      npthlev_bif, npthout_local, totalnpthout, pth_global_id
   USE MOD_Grid_Reservoir, only: numresv, totalnumresv, resv_data_address
   USE MOD_Grid_RiverLakeVectorIO, only: vector_gather_matrix_to_master, vector_read_matrix_and_scatter
   IMPLICIT NONE

   character(len=*), intent(in) :: file
   integer,          intent(in) :: mode

   logical :: do_levee, do_bif, do_resv, do_path
   integer, allocatable  :: gid (:)
   real(r8), allocatable :: work (:,:), wout (:,:)

      do_levee = DEF_USE_LEVEE .and. allocated(a_levsto)
      do_bif   = DEF_USE_BIFURCATION .and. allocated(a_bifout)
      do_resv  = DEF_Reservoir_Method > 0 .and. totalnumresv > 0 .and. allocated(acctime_resv)
      do_path  = DEF_USE_BIFURCATION .and. totalnpthout > 0 .and. npthlev_bif > 0 &
         .and. allocated(a_bifflw_lev) .and. allocated(a_bifflw_acctime)

      IF (mode == HACC_WRITE) THEN
         IF (p_is_master) THEN
            CALL ncio_create_file (file)
            IF (totalnumucat > 0) CALL ncio_define_dimension (file, 'ucatch',    totalnumucat)
            IF (do_resv)          CALL ncio_define_dimension (file, 'reservoir', totalnumresv)
            IF (do_path) THEN
               CALL ncio_define_dimension (file, 'bifurcation_level',   npthlev_bif)
               CALL ncio_define_dimension (file, 'bifurcation_pathway', totalnpthout)
               CALL ncio_define_dimension (file, 'bifurcation_scalar',  1)
            ENDIF
         ENDIF
#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif
      ELSE
         IF (.not. serial_var_exists (file, 'hist_acctime_ucat') .and. totalnumucat > 0) THEN
            IF (p_is_master) write(*,'(A)') 'ERROR: missing river history-accumulator restart '//trim(file)
            CALL CoLM_stop ()
         ENDIF
      ENDIF

      IF (totalnumucat > 0) THEN
         CALL ucat_sum (file, mode, 'hist_acctime_ucat', acctime_ucat)
         CALL ucat_sum (file, mode, 'hist_wdsrf_ucat',   a_wdsrf_ucat)
         CALL ucat_sum (file, mode, 'hist_veloc_riv',    a_veloc_riv )
         CALL ucat_sum (file, mode, 'hist_discharge',    a_discharge )
         CALL ucat_sum (file, mode, 'hist_floodarea',    a_floodarea )
         CALL ucat_sum (file, mode, 'hist_rivsto',       a_rivsto    )
         CALL ucat_sum (file, mode, 'hist_fldsto',       a_fldsto    )
         CALL ucat_sum (file, mode, 'hist_flddph',       a_flddph    )
         CALL ucat_sum (file, mode, 'hist_storge',       a_storge    )
         CALL ucat_sum (file, mode, 'hist_sfcelv',       a_sfcelv    )
         IF (do_levee) THEN
            CALL ucat_sum (file, mode, 'hist_levsto', a_levsto)
            CALL ucat_sum (file, mode, 'hist_levdph', a_levdph)
         ENDIF
         IF (do_bif) CALL ucat_sum (file, mode, 'hist_bifout', a_bifout)
      ENDIF

      IF (do_resv) THEN
         CALL resv_sum (file, mode, 'hist_acctime_resv', acctime_resv)
         CALL resv_sum (file, mode, 'hist_volresv',      a_volresv   )
         CALL resv_sum (file, mode, 'hist_qresv_in',     a_qresv_in  )
         CALL resv_sum (file, mode, 'hist_qresv_out',    a_qresv_out )
      ENDIF

      IF (do_path) THEN
         IF (p_is_worker) THEN
            allocate (gid (npthout_local))
            gid = pth_global_id(1:npthout_local)
         ELSE
            allocate (gid (0))
         ENDIF

         IF (mode == HACC_WRITE) THEN
            CALL vector_gather_matrix_to_master (a_bifflw_lev, npthlev_bif, size(gid), &
               totalnpthout, gid, wout)
            IF (p_is_master) THEN
               CALL ncio_write_serial (file, 'hist_bifflw_lev', wout, &
                  'bifurcation_level', 'bifurcation_pathway', DEF_REST_CompressLevel)
               deallocate (wout)
            ENDIF

            allocate (work (1, size(gid)))
            IF (size(gid) > 0) work(1,:) = a_bifflw_acctime(1:size(gid))
            CALL vector_gather_matrix_to_master (work, 1, size(gid), totalnpthout, gid, wout)
            IF (p_is_master) THEN
               CALL ncio_write_serial (file, 'hist_bifflw_acctime', wout, &
                  'bifurcation_scalar', 'bifurcation_pathway', DEF_REST_CompressLevel)
               deallocate (wout)
            ENDIF
            deallocate (work)
         ELSE
            CALL require_serial (file, 'hist_bifflw_lev')
            CALL require_serial (file, 'hist_bifflw_acctime')
            CALL vector_read_matrix_and_scatter (file, a_bifflw_lev, npthlev_bif, size(gid), &
               'hist_bifflw_lev', gid, totalnpthout)
            allocate (work (1, size(gid)))
            work = 0._r8
            CALL vector_read_matrix_and_scatter (file, work, 1, size(gid), &
               'hist_bifflw_acctime', gid, totalnpthout)
            IF (size(gid) > 0) a_bifflw_acctime(1:size(gid)) = work(1,:)
            deallocate (work)
         ENDIF

         deallocate (gid)
      ENDIF

   CONTAINS

      SUBROUTINE ucat_sum (file, mode, name, var)
      USE MOD_Vector_ReadWrite, only: vector_gather_and_write, vector_read_and_scatter
      character(len=*),      intent(in)    :: file, name
      integer,               intent(in)    :: mode
      real(r8), allocatable, intent(inout) :: var (:)
         IF (mode == HACC_WRITE) THEN
            CALL vector_gather_and_write (var, numucat, totalnumucat, ucat_data_address, &
               file, name, 'ucatch')
         ELSE
            CALL require_serial (file, name)
            CALL vector_read_and_scatter (file, var, numucat, name, ucat_data_address)
         ENDIF
      END SUBROUTINE ucat_sum

      SUBROUTINE resv_sum (file, mode, name, var)
      USE MOD_Vector_ReadWrite, only: vector_gather_and_write, vector_read_and_scatter
      character(len=*),      intent(in)    :: file, name
      integer,               intent(in)    :: mode
      real(r8), allocatable, intent(inout) :: var (:)
         IF (mode == HACC_WRITE) THEN
            CALL vector_gather_and_write (var, numresv, totalnumresv, resv_data_address, &
               file, name, 'reservoir')
         ELSE
            CALL require_serial (file, name)
            CALL vector_read_and_scatter (file, var, numresv, name, resv_data_address)
         ENDIF
      END SUBROUTINE resv_sum

   END SUBROUTINE river_sums

   !-----------------------------------------------------------------------
   SUBROUTINE require_serial (file, name)

   IMPLICIT NONE
   character(len=*), intent(in) :: file, name

      IF (.not. serial_var_exists (file, name)) THEN
         IF (p_is_master) write(*,'(A)') 'ERROR: '//trim(name)//' missing in '//trim(file)
         CALL CoLM_stop ()
      ENDIF

   END SUBROUTINE require_serial

   !-----------------------------------------------------------------------
   logical FUNCTION serial_var_exists (file, name)

   ! Only the master opens the file (see the note in
   ! MOD_Grid_RiverLakeBifurcation on concurrent opens of HDF5 files).

   USE MOD_NetCDFSerial, only: ncio_var_exist
   IMPLICIT NONE
   character(len=*), intent(in) :: file, name
   logical :: fexists

      serial_var_exists = .false.
      IF (p_is_master) THEN
         inquire (file=trim(file), exist=fexists)
         IF (fexists) serial_var_exists = ncio_var_exist (file, name, readflag=.false.)
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast (serial_var_exists, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif

   END FUNCTION serial_var_exists
#endif

   !-----------------------------------------------------------------------
   SUBROUTINE sidecar_filename (jdate, file)

   USE MOD_Namelist, only: DEF_dir_restart, DEF_CASE_NAME
   IMPLICIT NONE

   integer,          intent(in)  :: jdate(3)
   character(len=*), intent(out) :: file
   character(len=14) :: cdate

      write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') jdate
      file = trim(DEF_dir_restart)//'/'//trim(cdate)//'/'//trim(DEF_CASE_NAME) &
         //'_restart_hist_'//trim(cdate)//'.nc'

   END SUBROUTINE sidecar_filename

   !-----------------------------------------------------------------------
   integer FUNCTION freq_code (freq)

   IMPLICIT NONE
   character(len=*), intent(in) :: freq

      SELECT CASE (trim(adjustl(freq)))
      CASE ('TIMESTEP')
         freq_code = 1
      CASE ('HOURLY')
         freq_code = 2
      CASE ('DAILY')
         freq_code = 3
      CASE ('MONTHLY')
         freq_code = 4
      CASE ('YEARLY')
         freq_code = 5
      CASE DEFAULT
         freq_code = 0
      END SELECT

   END FUNCTION freq_code

   !-----------------------------------------------------------------------
   SUBROUTINE write_marker (file, name, val)

   ! One value per patch, so that the marker is in every block file.

   USE MOD_LandPatch, only: landpatch, numpatch
   USE MOD_NetCDFVector, only: ncio_write_vector
   IMPLICIT NONE

   character(len=*), intent(in) :: file, name
   real(r8),         intent(in) :: val
   real(r8), allocatable :: marker (:)

      IF (p_is_worker) THEN
         allocate (marker (numpatch))
      ELSE
         allocate (marker (0))
      ENDIF
      marker(:) = val
      CALL ncio_write_vector (file, name, 'patch', landpatch, marker)
      deallocate (marker)

   END SUBROUTINE write_marker

   !-----------------------------------------------------------------------
   SUBROUTINE read_marker (file, name, vmin, vmax)

   ! Global minimum and maximum of a marker over all patches.

   USE, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   USE MOD_LandPatch, only: landpatch, numpatch
   USE MOD_NetCDFVector, only: ncio_read_vector
   IMPLICIT NONE

   character(len=*), intent(in)  :: file, name
   real(r8),         intent(out) :: vmin, vmax
   real(r8), allocatable :: marker (:)
   integer :: bad(1)

      CALL ncio_read_vector (file, name, landpatch, marker)
      ! MINVAL/MAXVAL can hide a NaN among otherwise valid marker values.
      bad = 0
      IF (allocated(marker)) THEN
         IF (any(.not. ieee_is_finite(marker))) bad = 1
      ENDIF
      CALL global_max (bad)
      IF (bad(1) /= 0) CALL CoLM_stop ('Non-finite history-accumulator marker '//trim(name))
      vmin =  huge(vmin)
      vmax = -huge(vmax)
      IF (p_is_worker .and. numpatch > 0) THEN
         vmin = minval(marker)
         vmax = maxval(marker)
      ENDIF
#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, vmin, 1, MPI_REAL8, MPI_MIN, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, vmax, 1, MPI_REAL8, MPI_MAX, p_comm_glb, p_err)
#endif
      IF (allocated(marker)) deallocate (marker)

   END SUBROUTINE read_marker

   !-----------------------------------------------------------------------
   SUBROUTINE global_max (n)

   IMPLICIT NONE
   integer, intent(inout) :: n(:)

#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, n, size(n), MPI_INTEGER, MPI_MAX, p_comm_glb, p_err)
#endif

   END SUBROUTINE global_max

   !-----------------------------------------------------------------------
   logical FUNCTION any_block_exists (file)

   USE MOD_Block, only: get_filename_block
   USE MOD_LandPatch, only: landpatch
   IMPLICIT NONE

   character(len=*), intent(in) :: file

   character(len=256) :: fileblock
   integer :: iblkgrp, found(1)
   logical :: fexists

      found = 0
      IF (p_is_io) THEN
         DO iblkgrp = 1, landpatch%nblkgrp
            CALL get_filename_block (file, landpatch%xblkgrp(iblkgrp), landpatch%yblkgrp(iblkgrp), fileblock)
            inquire (file=trim(fileblock), exist=fexists)
            IF (fexists) found = 1
         ENDDO
      ENDIF
      CALL global_max (found)
      any_block_exists = found(1) > 0

   END FUNCTION any_block_exists

   !-----------------------------------------------------------------------
   SUBROUTINE require_var (file, name, pixelset, vecname, rank, n1, n2)

   ! Stop unless every block file of the pixelset holds a double variable
   ! with the expected dimensions: (vecname), (d1_name, vecname) or
   ! (d1_name, d2_name, vecname), the vector dimension having the block's
   ! number of elements.

   USE netcdf
   USE MOD_Block, only: get_filename_block
   USE MOD_Pixelset, only: pixelset_type
   IMPLICIT NONE

   character(len=*),    intent(in) :: file, name, vecname
   type(pixelset_type), intent(in) :: pixelset
   integer,             intent(in) :: rank, n1, n2

   character(len=256) :: fileblock
   character(len=NF90_MAX_NAME) :: dname(3), expect(3)
   integer :: iblkgrp, iblk, jblk, ncid, varid, xtype, ndims, dimids(3), dlen(3), j, bad(1)

      bad = 0
      IF (p_is_io) THEN
         expect(1) = 'd1_'//trim(name)
         expect(2) = 'd2_'//trim(name)
         expect(rank) = vecname
         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)
            CALL get_filename_block (file, iblk, jblk, fileblock)
            IF (nf90_open (trim(fileblock), NF90_NOWRITE, ncid) /= NF90_NOERR) THEN
               bad = 1
               CYCLE
            ENDIF
            IF (nf90_inq_varid (ncid, trim(name), varid) /= NF90_NOERR) THEN
               bad = 1
            ELSEIF (nf90_inquire_variable (ncid, varid, xtype=xtype, ndims=ndims, dimids=dimids) &
                  /= NF90_NOERR) THEN
               bad = 1
            ELSEIF (xtype /= NF90_DOUBLE .or. ndims /= rank) THEN
               bad = 1
            ELSE
               DO j = 1, rank
                  IF (nf90_inquire_dimension (ncid, dimids(j), name=dname(j), len=dlen(j)) /= NF90_NOERR) &
                     bad = 1
               ENDDO
               IF (bad(1) == 0) THEN
                  IF (any(dname(1:rank) /= expect(1:rank))) bad = 1
                  IF (dlen(rank) /= pixelset%vecgs%vlen(iblk,jblk)) bad = 1
                  IF (rank >= 2) THEN
                     IF (dlen(1) /= n1) bad = 1
                  ENDIF
                  IF (rank == 3) THEN
                     IF (dlen(2) /= n2) bad = 1
                  ENDIF
               ENDIF
            ENDIF
            IF (nf90_close (ncid) /= NF90_NOERR) bad = 1
            IF (bad(1) /= 0) THEN
               write(*,'(A)') 'ERROR: '//trim(name)//' missing or of wrong shape in '//trim(fileblock)
               EXIT
            ENDIF
         ENDDO
      ENDIF

      CALL global_max (bad)
      IF (bad(1) /= 0) CALL CoLM_stop ()

   END SUBROUTINE require_var

   !-----------------------------------------------------------------------
   logical FUNCTION marker_present (file, name, pixelset, vecname)

   ! Non-fatal existence check: .true. when every block file of the pixel set
   ! holds `name`.  Used to detect a sidecar left incomplete by a killed run
   ! so it can be ignored instead of aborting the restart.

   USE netcdf
   USE MOD_Block, only: get_filename_block
   USE MOD_Pixelset, only: pixelset_type
   IMPLICIT NONE

   character(len=*),    intent(in) :: file, name, vecname
   type(pixelset_type), intent(in) :: pixelset

   character(len=256) :: fileblock
   integer :: iblkgrp, iblk, jblk, ncid, varid, bad(1)
   logical :: ok

      bad = 0
      IF (p_is_io) THEN
         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)
            CALL get_filename_block (file, iblk, jblk, fileblock)
            IF (nf90_open (trim(fileblock), NF90_NOWRITE, ncid) /= NF90_NOERR) THEN
               bad = 1
               CYCLE
            ENDIF
            ok = (nf90_inq_varid (ncid, trim(name), varid) == NF90_NOERR)
            IF (.not. ok) bad = 1
            IF (nf90_close (ncid) /= NF90_NOERR) bad = 1
            IF (bad(1) /= 0) EXIT
         ENDDO
      ENDIF
      CALL global_max (bad)
      marker_present = (bad(1) == 0)

   END FUNCTION marker_present

END MODULE MOD_HistAccRestart
