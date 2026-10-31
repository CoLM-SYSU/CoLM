#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_HostFixes

!-----------------------------------------------------------------------
! Fixes for two upstream problems that break restart reproducibility.
! Material runs need a split run (restart) to reproduce an unbroken
! run exactly, so both fixes are applied here, and ONLY when at least one
! material is switched on (mt_any_enabled). With materials off, or
! with MaterialTransport not defined, the host keeps its upstream
! behaviour bit for bit: the host calls below are inside
! #ifdef MaterialTransport and both functions return the upstream
! value when mt_any_enabled() is false.
!
! 1. Start year of the BGC nitrogen-deposition / fire data
!    (mt_bgc_data_init_year). CoLM.F90 passes sdate(1) after
!    adj2end(sdate). A run starting at 00:00 on 1 January has its start
!    date moved to the end of the previous year (day 365/366, 86400 s),
!    so the data of the previous year are read and kept for the whole
!    first year (the annual update in the time loop only fires when the
!    year changes). A run restarted later in the same year reads the
!    right year, so the split run and the unbroken run differ (ndep,
!    sminn_vr, ...). The year that the simulation actually starts in is
!    DEF_simulation_time%start_year (s_year), as CoLM202X uses.
!
! 2. History file continued after a restart (mt_hist_continue_file).
!    The host creates a new history file when the file name differs from
!    the last file written in this run; that name is 'null' at start-up,
!    so a run restarted inside a history group (e.g. on day 2 with
!    monthly files) creates the file again and loses the records written
!    before the restart. When the run starts inside the group of the file
!    it is about to write and the file exists, it is marked as the last
!    file, so the host appends to it (ncio_write_time finds the time
!    record or appends a new one). For gridded history in 'block' mode,
!    every per-block file of the IO processes must exist. CaMa history
!    files need nothing: they are only created when missing.
!-----------------------------------------------------------------------

   USE MOD_MT_Namelist, only: mt_any_enabled
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: mt_bgc_data_init_year
   PUBLIC :: mt_hist_continue_file

CONTAINS

   integer FUNCTION mt_bgc_data_init_year (sdate_year, start_year)
   ! sdate_year: sdate(1) after adj2end (upstream argument)
   ! start_year: DEF_simulation_time%start_year
   IMPLICIT NONE
   integer, intent(in) :: sdate_year, start_year

      IF (mt_any_enabled()) THEN
         mt_bgc_data_init_year = start_year
      ELSE
         mt_bgc_data_init_year = sdate_year
      ENDIF

   END FUNCTION mt_bgc_data_init_year


   SUBROUTINE mt_hist_continue_file (file_hist, idate, histform, file_last)
   ! file_hist: the history file about to be written
   ! idate:     current model time (end of the time step)
   ! histform:  'Gridded', 'Vector' or 'Single' (MOD_Hist HistForm)
   ! file_last: the host's last-written file name; set to file_hist when
   !            the existing file has to be continued
   USE MOD_SPMD_Task
   USE MOD_Namelist, only: DEF_simulation_time, DEF_HIST_groupby, DEF_HIST_mode
   USE MOD_TimeManager, only: julian2monthday, monthday2julian, minutes_since_1900
   USE MOD_Block, only: gblock, get_filename_block
   USE MOD_HistGridded, only: ghist
   IMPLICIT NONE
   character(len=*), intent(in)    :: file_hist
   integer,          intent(in)    :: idate(3)
   character(len=*), intent(in)    :: histform
   character(len=*), intent(inout) :: file_last

   integer :: month, mday, jday_group, jday_start, min_group, min_start
   integer :: iblkme, iblk, jblk
   character(len=256) :: fileblock
   logical :: lcontinue, fexists, per_block

      IF (.not. mt_any_enabled()) RETURN
      ! only the first write of this run
      IF (trim(file_last) /= 'null') RETURN
      ! only gridded history has per-block files
      per_block = (trim(histform) == 'Gridded') .and. (trim(DEF_HIST_mode) == 'block')

      ! beginning of the history group that idate belongs to
      ! (idate is end-of-step, so 86400 s still belongs to that day)
      SELECT CASE (trim(DEF_HIST_groupby))
      CASE ('YEAR')
         jday_group = 1
      CASE ('MONTH')
         CALL julian2monthday (idate(1), idate(2), month, mday)
         CALL monthday2julian (idate(1), month, 1, jday_group)
      CASE ('DAY')
         jday_group = idate(2)
      CASE DEFAULT
         RETURN
      END SELECT
      min_group = minutes_since_1900 (idate(1), jday_group, 0)

      CALL monthday2julian (DEF_simulation_time%start_year, DEF_simulation_time%start_month, &
         DEF_simulation_time%start_day, jday_start)
      min_start = minutes_since_1900 (DEF_simulation_time%start_year, jday_start, &
         DEF_simulation_time%start_sec)

      ! A run starting at the beginning of the group writes a new file, as upstream.
      IF (min_start <= min_group) RETURN

      IF (per_block) THEN
         ! Each IO process writes its own block files (same loop as
         ! hist_gridded_write_time); continue only if all of them exist.
         lcontinue = .true.
         IF (p_is_io) THEN
            DO iblkme = 1, gblock%nblkme
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)
               IF (ghist%ycnt(jblk) <= 0) CYCLE
               IF (ghist%xcnt(iblk) <= 0) CYCLE
               CALL get_filename_block (file_hist, iblk, jblk, fileblock)
               inquire (file=trim(fileblock), exist=fexists)
               lcontinue = lcontinue .and. fexists
            ENDDO
         ENDIF
#ifdef USEMPI
         CALL mpi_allreduce (MPI_IN_PLACE, lcontinue, 1, MPI_LOGICAL, MPI_LAND, p_comm_glb, p_err)
#endif
      ELSE
         ! The single file is opened by the master only; decide there so that
         ! all processes agree.
         lcontinue = .false.
         IF (p_is_master) inquire (file=trim(file_hist), exist=lcontinue)
#ifdef USEMPI
         CALL mpi_bcast (lcontinue, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif
      ENDIF

      IF (lcontinue) THEN
         file_last = file_hist
         IF (p_is_master) write(*,'(A)') 'Material restart: continuing history file ' // trim(file_hist)
      ENDIF

   END SUBROUTINE mt_hist_continue_file

END MODULE MOD_MT_HostFixes
#endif
