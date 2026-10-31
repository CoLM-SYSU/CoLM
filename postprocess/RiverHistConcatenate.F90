#include <define.h>

PROGRAM river_hist_concatenate

!-----------------------------------------------------------------------
! DESCRIPTION:
!
!    Reassemble the IO-group route-history shards written under
!    DEF_HIST_mode='block' into the single `_hist_unitcat_*.nc` file that
!    DEF_HIST_mode='one' produces, so existing analysis scripts need no
!    change.
!
!    Correctness is decided by the shard identity attributes, never by the
!    filename: schema version, case, target file, history period and grid
!    fingerprint must all agree before anything is merged, and the union of
!    global ids must cover the domain exactly once.  Shards from a different
!    restart segment ARE merged -- a run that restarts mid-period writes two
!    segments of the same period -- but a time record appearing twice with
!    different values is a conflict and stops the job.
!
!    Output is written to <target>.tmp, verified, closed and only then
!    renamed onto <target>, with a <target>.complete marker.  A half-written
!    aggregate must never be mistaken for the finished product.
!
! Usage:
!    river_hist_concatenate <colm_namelist> <target_hist_unitcat_file>
!
!    The shards are expected next to the target as <target minus .nc>_shardNNNNN.nc.
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE netcdf

   IMPLICIT NONE

   ! 2: shard filenames carry the run segment (<base>_<seg>_shardNNNNN.nc).
   !    A version-1 shard set has no segment in its name and is invisible to
   !    the discovery glob, so it must be refused by version rather than
   !    silently reported as "no shards found".
   integer, parameter :: SUPPORTED_SCHEMA_VERSION = 2
   real(r8), parameter :: SPVAL = -1.e36_r8

   character(len=256) :: nlfile, filetarget, filetmp, filedone
   character(len=256) :: fname, ref_case, ref_period, ref_grid, ref_target
   integer :: shard_count, ish, ncid, ierr, nlon, nlat, ntime, nvar, ivar
   integer :: totalnumucat, totalnpthout, npthlev
   integer, allocatable :: seen_ucat(:), seen_path(:)
   ! Run segments. A restart inside one history period writes a NEW complete
   ! shard set under a new segment id, so the period is the union of segments
   ! and each time record lives in exactly one of them.
   integer, parameter :: MAX_SEG = 64
   character(len=64) :: seg_id(MAX_SEG)
   integer :: seg_shards(MAX_SEG), seg_ntime(MAX_SEG), nseg
   real(r8), allocatable :: seg_time(:,:)          ! (record, segment)
   integer,  allocatable :: rec_seg(:), rec_local(:)  ! global record -> segment, local
   character(len=256), allocatable :: varnames(:)
   real(r8), allocatable :: tvals_global(:)
   integer,  allocatable :: vals_i(:)
   logical :: ok

#ifdef USEMPI
      CALL spmd_init
#endif

      CALL getarg (1, nlfile)
      CALL getarg (2, filetarget)

      IF (len_trim(nlfile) == 0 .or. len_trim(filetarget) == 0) THEN
         write(*,'(A)') 'Usage: river_hist_concatenate <colm_namelist> <target_hist_unitcat_file>'
         CALL CoLM_stop ('river_hist_concatenate: missing arguments')
      ENDIF

      ! BEFORE anything shells out. scan_shards builds an ls command with the
      ! target single-quoted, and a single quote is the one character single
      ! quotes cannot contain -- a target of  x'; rm -rf ~; #.nc  closes the
      ! quote and the rest runs as commands. Checking this inside rename_file
      ! was useless: that runs at the very end, long after the injected
      ! command has executed.
      CALL reject_unquotable (nlfile,     'namelist')
      CALL reject_unquotable (filetarget, 'target')

      CALL read_namelist (nlfile)

#ifdef USEMPI
      ! Bug fix: the merge runs on the master alone, but read_namelist
      ! broadcasts over p_comm_glb, so every process has to take part in it.
      ! The other processes used to leave before it, into the barrier of
      ! spmd_exit, while the master was in the broadcast: unmatched collectives
      ! on one communicator (a hang or an MPI error with more than one process).
      ! Now they leave after it; their spmd_exit barrier pairs with the one the
      ! master enters at the end.
      IF (.not. p_is_master) THEN
         CALL spmd_exit
         STOP
      ENDIF
#endif

      filetmp  = trim(filetarget) // '.tmp'
      filedone = trim(filetarget) // '.complete'

      write(*,'(A,A)') 'river_hist_concatenate: target ', trim(filetarget)

      CALL scan_shards ()
      CALL build_output ()
      CALL verify_and_promote ()

#ifdef USEMPI
      CALL spmd_exit
#endif

CONTAINS

   !> Refuse a path that cannot be safely single-quoted in a shell command.
   !!
   !! Every other metacharacter -- space, $, ;, &, backtick, newline -- is
   !! inert inside single quotes. A single quote is not: it ends the quoted
   !! string and hands the remainder to the shell. Rather than attempt to
   !! escape it, refuse the path. This must run before the first shell
   !! command, not before the last one.
   SUBROUTINE reject_unquotable (path, what)
   character(len=*), intent(in) :: path, what
      IF (index(path, "'") > 0) THEN
         write(*,'(4A)') 'ERROR: the ', trim(what), ' path contains a single quote, ', &
            'which cannot be passed safely to the shell:'
         write(*,'(2A)') '       ', trim(path)
         CALL CoLM_stop ('river_hist_concatenate: unsupported path')
      ENDIF
   END SUBROUTINE reject_unquotable

   !> Discover every run segment of this period and check identity agreement.
   !!
   !! A restart writes a new segment rather than reusing the previous shards,
   !! so the IO-group count may legitimately differ between segments. What must
   !! agree across all of them is the case, target file, period and grid; what
   !! must not overlap is the time records.
   SUBROUTINE scan_shards ()

   integer :: i, v, n, u, iseg, k, j
   character(len=256) :: c, listfile, line
   real(r8), allocatable :: tref(:), tcmp(:)
   logical :: fexists

      ! Enumerate segments from the shard-0 files on disk. Fortran has no
      ! portable directory listing, and the segment id is a timestamp that
      ! cannot be guessed, so shell out once.
      listfile = trim(filetarget) // '.segments.tmp'
      ! Single-quote both paths: they come from the command line, so a space
      ! would split the argument and a metacharacter would be interpreted.
      ! The glob must stay outside the quotes to remain a glob.
      CALL system ("ls '" // trim(stem()) // "'_seg*_shard00000.nc 2>/dev/null > '" &
         // trim(listfile) // "'")

      nseg = 0
      OPEN (newunit=u, file=trim(listfile), status='old', action='read', iostat=ierr)
      IF (ierr == 0) THEN
         DO
            READ (u,'(A)',iostat=ierr) line
            IF (ierr /= 0) EXIT
            IF (len_trim(line) == 0) CYCLE
            IF (nseg >= MAX_SEG) CALL CoLM_stop &
               ('river_hist_concatenate: too many run segments')
            nseg = nseg + 1
            seg_id(nseg) = segment_of (trim(line))
         ENDDO
         CLOSE (u, status='delete')
      ENDIF

      IF (nseg == 0) CALL CoLM_stop &
         ('river_hist_concatenate: no shards found for '//trim(filetarget))

      shard_count = -1
      DO iseg = 1, nseg
         ! count this segment's shards and check identity
         n = 0
         DO
            CALL shard_name (iseg, n, fname)
            inquire (file=trim(fname), exist=fexists)
            IF (.not. fexists) EXIT
            n = n + 1
         ENDDO
         seg_shards(iseg) = n

         DO ish = 0, n-1
            CALL shard_name (iseg, ish, fname)
            ierr = nf90_open (trim(fname), NF90_NOWRITE, ncid)
            IF (ierr /= NF90_NOERR) CALL CoLM_stop &
               ('river_hist_concatenate: cannot open '//trim(fname))

            CALL get_att_int (ncid, 'river_hist_shard_schema_version', v)
            IF (v /= SUPPORTED_SCHEMA_VERSION) THEN
               write(*,'(A,I0,A,I0,A,A)') 'ERROR: shard schema version ', v, &
                  ' is not supported (this build reads ', SUPPORTED_SCHEMA_VERSION, &
                  '): ', trim(fname)
               CALL CoLM_stop ('river_hist_concatenate: unsupported shard schema')
            ENDIF

            CALL get_att_int (ncid, 'shard_count', k)
            IF (k /= n) CALL identity_error (fname, 'shard_count')

            IF (shard_count < 0) THEN
               shard_count = n
               CALL get_att_str (ncid, 'case_name',               ref_case)
               CALL get_att_str (ncid, 'history_period_key',      ref_period)
               CALL get_att_str (ncid, 'grid_fingerprint',        ref_grid)
               CALL get_att_str (ncid, 'target_history_basename', ref_target)
            ELSE
               CALL get_att_str (ncid, 'case_name', c)
               IF (trim(c) /= trim(ref_case))   CALL identity_error (fname, 'case_name')
               CALL get_att_str (ncid, 'history_period_key', c)
               IF (trim(c) /= trim(ref_period)) CALL identity_error (fname, 'history_period_key')
               CALL get_att_str (ncid, 'grid_fingerprint', c)
               IF (trim(c) /= trim(ref_grid))   CALL identity_error (fname, 'grid_fingerprint')
               CALL get_att_str (ncid, 'target_history_basename', c)
               IF (trim(c) /= trim(ref_target)) CALL identity_error (fname, 'target_history_basename')
            ENDIF

            ! Identity is what decides membership, so check the attributes
            ! that name this file's place in the set -- not just the ones
            ! shared across the set. A shard whose recorded index disagrees
            ! with its filename means the set is not what it claims.
            CALL get_att_str (ncid, 'segment_id', c)
            IF (trim(c) /= trim(seg_id(iseg))) CALL identity_error (fname, 'segment_id')
            CALL get_att_int (ncid, 'shard_index', k)
            IF (k /= ish) CALL identity_error (fname, 'shard_index')

            ! Every shard of a segment must carry the SAME time axis: the
            ! aggregator addresses them all by one record index, so a shard
            ! whose axis differs would be silently combined at the wrong
            ! instant.
            IF (ish == 0) THEN
               CALL read_dim  (ncid, 'time', seg_ntime(iseg))
               CALL read_real (ncid, 'time', tref)
            ELSE
               CALL read_dim  (ncid, 'time', ntime)
               IF (ntime /= seg_ntime(iseg)) THEN
                  write(*,'(A,A,A,I0,A,I0)') 'ERROR: ', trim(fname), &
                     ' has ', ntime, ' time records but shard 0 has ', seg_ntime(iseg)
                  CALL CoLM_stop ('river_hist_concatenate: shard time axis length differs')
               ENDIF
               CALL read_real (ncid, 'time', tcmp)
               DO j = 1, seg_ntime(iseg)
                  IF (tcmp(j) /= tref(j)) THEN
                     write(*,'(A,A,A,I0,A,2ES16.8)') 'ERROR: ', trim(fname), &
                        ' time record ', j, ' differs from shard 0: ', tcmp(j), tref(j)
                     CALL CoLM_stop ('river_hist_concatenate: shard time axis differs')
                  ENDIF
               ENDDO
               IF (allocated(tcmp)) deallocate (tcmp)
            ENDIF
            ierr = nf90_close (ncid)
         ENDDO

         IF (allocated(tref)) deallocate (tref)
         write(*,'(A,A,A,I0,A,I0,A)') '  segment ', trim(seg_id(iseg)), ': ', &
            seg_shards(iseg), ' shards, ', seg_ntime(iseg), ' time records'
      ENDDO

      CALL build_time_map ()

   END SUBROUTINE scan_shards

   !> Global time axis = union of the segments, with each record owned by
   !! exactly one segment. The same time appearing twice is a conflict, not a
   !! merge: two runs wrote the same instant and we cannot know which is right.
   SUBROUTINE build_time_map ()

   integer :: iseg, it, g, j
   real(r8), allocatable :: tv(:)
   real(r8) :: tbest
   integer  :: nbest_s, nbest_r, total

      allocate (seg_time(maxval(seg_ntime(1:nseg)), nseg)); seg_time = -huge(1._r8)
      DO iseg = 1, nseg
         CALL shard_name (iseg, 0, fname)
         ierr = nf90_open (trim(fname), NF90_NOWRITE, ncid)
         CALL read_real (ncid, 'time', tv)
         DO it = 1, seg_ntime(iseg)
            seg_time(it, iseg) = tv(it)
         ENDDO
         IF (allocated(tv)) deallocate (tv)
         ierr = nf90_close (ncid)
      ENDDO

      total = sum(seg_ntime(1:nseg))
      allocate (rec_seg(total), rec_local(total))
      allocate (tvals_global(total))

      ! selection sort over the (segment, record) pairs by time value
      DO g = 1, total
         tbest = huge(1._r8); nbest_s = 0; nbest_r = 0
         DO iseg = 1, nseg
            DO it = 1, seg_ntime(iseg)
               IF (seg_time(it,iseg) < tbest) THEN
                  tbest = seg_time(it,iseg); nbest_s = iseg; nbest_r = it
               ENDIF
            ENDDO
         ENDDO
         IF (nbest_s == 0) CALL CoLM_stop ('river_hist_concatenate: time map failed')
         IF (g > 1) THEN
            IF (tbest == tvals_global(g-1)) THEN
               write(*,'(A,ES16.8,A)') 'ERROR: time ', tbest, &
                  ' appears in more than one run segment'
               CALL CoLM_stop ('river_hist_concatenate: conflicting time record')
            ENDIF
         ENDIF
         tvals_global(g) = tbest
         rec_seg(g)   = nbest_s
         rec_local(g) = nbest_r
         seg_time(nbest_r, nbest_s) = huge(1._r8)   ! consume
      ENDDO

      ntime = total
      write(*,'(A,I0,A,I0,A)') '  ', ntime, ' time records across ', nseg, ' segment(s)'

   END SUBROUTINE build_time_map

   !> filetarget without the trailing .nc, for glob construction.
   FUNCTION stem () RESULT (t)
   character(len=256) :: t
   integer :: i
      i = len_trim(filetarget)
      IF (i > 3 .and. filetarget(max(i-2,1):i) == '.nc') THEN
         t = filetarget(1:i-3)
      ELSE
         t = trim(filetarget)
      ENDIF
   END FUNCTION stem

   !> '<stem>_<seg>_shard00000.nc' -> '<seg>'
   FUNCTION segment_of (path) RESULT (sg)
   character(len=*), intent(in) :: path
   character(len=64) :: sg
   integer :: a, b
      sg = ''
      b = index(path, '_shard', back=.true.)
      IF (b <= 1) RETURN
      a = index(path(1:b-1), '_seg', back=.true.)
      IF (a <= 0) RETURN
      sg = path(a+1:b-1)
   END FUNCTION segment_of

   !> Rebuild every variable on the global grid / pathway axis.
   SUBROUTINE build_output ()

   integer :: k
   integer, allocatable :: ids(:), xs(:), ys(:)
   real(r8), allocatable :: lon(:), lat(:), tvals(:)

      ! geometry and time come from shard 0; identity checks above guarantee
      ! every shard agrees.
      CALL shard_name (1, 0, fname)
      ierr = nf90_open (trim(fname), NF90_NOWRITE, ncid)

      CALL read_dim  (ncid, 'lon_ucat', nlon)
      CALL read_dim  (ncid, 'lat_ucat', nlat)
      CALL read_real (ncid, 'lon_ucat', lon)
      CALL read_real (ncid, 'lat_ucat', lat)
      CALL collect_varnames (ncid, varnames, nvar)
      ierr = nf90_close (ncid)

      CALL ncio_create_file (trim(filetmp))
      CALL ncio_define_dimension (trim(filetmp), 'time', 0)
      CALL ncio_define_dimension (trim(filetmp), 'lat_ucat', nlat)
      CALL ncio_define_dimension (trim(filetmp), 'lon_ucat', nlon)
      CALL ncio_write_serial (trim(filetmp), 'lat_ucat', lat, 'lat_ucat')
      CALL ncio_write_serial (trim(filetmp), 'lon_ucat', lon, 'lon_ucat')
      CALL ncio_write_serial (trim(filetmp), 'time', tvals_global, 'time')
      ! ncio_write_serial writes values only; the units (the date origin) and
      ! long_name that ncio_write_time put on the shards are carried over
      ! here, as in the file DEF_HIST_mode='one' writes.
      CALL copy_var_attrs ('time')

      totalnumucat = 0
      DO ish = 0, seg_shards(rec_seg(1))-1
         CALL shard_name (rec_seg(1), ish, fname)
         ierr = nf90_open (trim(fname), NF90_NOWRITE, ncid)
         CALL read_int (ncid, 'ucat_ucid', ids)
         totalnumucat = totalnumucat + size(ids)
         IF (allocated(ids)) deallocate (ids)
         ierr = nf90_close (ncid)
      ENDDO
      allocate (seen_ucat(max(totalnumucat,1))); seen_ucat = 0

      DO ivar = 1, nvar
         CALL rebuild_ucat_variable (trim(varnames(ivar)))
      ENDDO

      ! Bug fix: the (time) window variables and the time-independent
      ! (unitcat_local) fields, e.g. mask_complete_upstream, are not
      ! (unitcat_local, time) fields, so collect_varnames skips them and they
      ! were missing from the aggregate although the 'one' file has them.
      CALL rebuild_window_variables ()
      CALL rebuild_static_ucat_variables ()

      CALL rebuild_resv_variables ()
      CALL rebuild_bif_matrix ()

      write(*,'(A,I0,A,I0,A)') '  rebuilt ', nvar, ' unit-catchment fields over ', &
         totalnumucat, ' catchments'

   END SUBROUTINE build_output

   !> One field: scatter each shard's values into the global grid by the
   !! (x,y) the shard stores for its own ids. Rank order is never consulted.
   SUBROUTINE rebuild_ucat_variable (varname)

   character(len=*), intent(in) :: varname
   real(r8), allocatable :: grid(:,:), vals(:)
   integer,  allocatable :: ids(:), xs(:), ys(:)
   integer :: it, k, vid, n

      allocate (grid(nlon,nlat))

      DO it = 1, ntime
         grid = SPVAL
         seen_ucat = 0

         DO ish = 0, seg_shards(rec_seg(it))-1
            CALL shard_name (rec_seg(it), ish, fname)
            ierr = nf90_open (trim(fname), NF90_NOWRITE, ncid)

            ierr = nf90_inq_varid (ncid, trim(varname), vid)
            IF (ierr /= NF90_NOERR) THEN
               ierr = nf90_close (ncid); CYCLE
            ENDIF

            CALL read_int (ncid, 'ucat_ucid', ids)
            CALL read_int (ncid, 'x_ucat',    xs)
            CALL read_int (ncid, 'y_ucat',    ys)
            n = size(ids)
            IF (n > 0) THEN
               CALL read_slice (ncid, vid, n, rec_local(it), vals)
               DO k = 1, n
                  IF (ids(k) < 1 .or. ids(k) > totalnumucat) THEN
                     write(*,'(A,I0,A,A)') 'ERROR: unit-catchment id ', ids(k), &
                        ' out of range in ', trim(fname)
                     CALL CoLM_stop ('river_hist_concatenate: id out of range')
                  ENDIF
                  IF (seen_ucat(ids(k)) /= 0) THEN
                     write(*,'(A,I0,A,A)') 'ERROR: unit-catchment id ', ids(k), &
                        ' appears in more than one shard, at ', trim(varname)
                     CALL CoLM_stop ('river_hist_concatenate: duplicate id')
                  ENDIF
                  seen_ucat(ids(k)) = 1
                  IF (xs(k) < 1 .or. xs(k) > nlon .or. ys(k) < 1 .or. ys(k) > nlat) THEN
                     CALL CoLM_stop ('river_hist_concatenate: grid index out of range')
                  ENDIF
                  grid(xs(k),ys(k)) = vals(k)
               ENDDO
               deallocate (vals)
            ENDIF
            IF (allocated(ids)) deallocate (ids)
            IF (allocated(xs))  deallocate (xs)
            IF (allocated(ys))  deallocate (ys)
            ierr = nf90_close (ncid)
         ENDDO

         IF (count(seen_ucat == 0) > 0) THEN
            write(*,'(A,I0,A,A)') 'ERROR: ', count(seen_ucat == 0), &
               ' unit catchments missing from the shard set for ', trim(varname)
            CALL CoLM_stop ('river_hist_concatenate: incomplete coverage')
         ENDIF

         CALL ncio_write_serial_time (trim(filetmp), trim(varname), it, grid, &
            'lon_ucat', 'lat_ucat', 'time', DEF_HIST_CompressLevel)
      ENDDO

      CALL copy_var_attrs (trim(varname))
      deallocate (grid)

   END SUBROUTINE rebuild_ucat_variable

   !> The per-record history window, (time) on every shard: each global
   !! record is taken from shard 0 of the segment that owns it.
   SUBROUTINE rebuild_window_variables ()

   character(len=*), parameter :: names(2) = &
      [character(len=32) :: 'history_window_seconds', 'history_window_end_minutes']
   integer  :: iv, it, vid, e
   real(r8) :: val(1)
   logical  :: present_on_shard

      DO iv = 1, size(names)
         CALL shard_name (1, 0, fname)
         ierr = nf90_open (trim(fname), NF90_NOWRITE, ncid)
         present_on_shard = nf90_inq_varid (ncid, trim(names(iv)), vid) == NF90_NOERR
         ierr = nf90_close (ncid)
         IF (.not. present_on_shard) CYCLE

         DO it = 1, ntime
            CALL shard_name (rec_seg(it), 0, fname)
            ierr = nf90_open (trim(fname), NF90_NOWRITE, ncid)
            e = nf90_inq_varid (ncid, trim(names(iv)), vid)
            IF (e == NF90_NOERR) e = nf90_get_var (ncid, vid, val, start=[rec_local(it)], count=[1])
            IF (e /= NF90_NOERR) THEN
               write(*,'(4A)') 'ERROR: cannot read ', trim(names(iv)), ' from ', trim(fname)
               CALL CoLM_stop ('river_hist_concatenate: missing history window record')
            ENDIF
            ierr = nf90_close (ncid)
            CALL ncio_write_serial_time (trim(filetmp), trim(names(iv)), it, val(1), 'time')
         ENDDO
         CALL copy_var_attrs (trim(names(iv)))
      ENDDO

   END SUBROUTINE rebuild_window_variables

   !> Time-independent real (unitcat_local) fields of the first segment, put on
   !! the grid by the (x,y) each shard stores for its own ids, as the 'one'
   !! writer does (missing cells SPVAL).  The id and index variables are
   !! integers and are not fields.
   SUBROUTINE rebuild_static_ucat_variables ()

   character(len=256) :: names(64), nm
   real(r8), allocatable :: grid(:,:), vals(:)
   integer,  allocatable :: ids(:), xs(:), ys(:)
   integer :: nv, v, nd, e, xtype, dimids(NF90_MAX_VAR_DIMS), udim
   integer :: nnames, iv, k, n, vid

      CALL shard_name (1, 0, fname)
      ierr = nf90_open (trim(fname), NF90_NOWRITE, ncid)
      nnames = 0
      e = nf90_inq_dimid (ncid, 'unitcat_local', udim)
      IF (e == NF90_NOERR) THEN
         ierr = nf90_inquire (ncid, nVariables=nv)
         DO v = 1, nv
            e = nf90_inquire_variable (ncid, v, name=nm, xtype=xtype, ndims=nd, dimids=dimids)
            IF (nd /= 1) CYCLE
            IF (dimids(1) /= udim) CYCLE
            IF (xtype /= NF90_DOUBLE .and. xtype /= NF90_FLOAT) CYCLE
            IF (nnames == size(names)) CALL CoLM_stop ('river_hist_concatenate: too many static fields')
            nnames = nnames + 1
            names(nnames) = nm
         ENDDO
      ENDIF
      ierr = nf90_close (ncid)
      IF (nnames == 0) RETURN

      allocate (grid(nlon,nlat))
      DO iv = 1, nnames
         grid = SPVAL
         seen_ucat = 0
         DO ish = 0, seg_shards(1)-1
            CALL shard_name (1, ish, fname)
            ierr = nf90_open (trim(fname), NF90_NOWRITE, ncid)
            IF (nf90_inq_varid (ncid, trim(names(iv)), vid) /= NF90_NOERR) THEN
               ierr = nf90_close (ncid); CYCLE
            ENDIF
            CALL read_int (ncid, 'ucat_ucid', ids)
            CALL read_int (ncid, 'x_ucat',    xs)
            CALL read_int (ncid, 'y_ucat',    ys)
            n = size(ids)
            IF (n > 0) THEN
               allocate (vals(n))
               e = nf90_get_var (ncid, vid, vals)
               IF (e /= NF90_NOERR) CALL CoLM_stop ('river_hist_concatenate: cannot read a static field')
               DO k = 1, n
                  IF (ids(k) < 1 .or. ids(k) > totalnumucat) &
                     CALL CoLM_stop ('river_hist_concatenate: id out of range')
                  IF (seen_ucat(ids(k)) /= 0) CALL CoLM_stop ('river_hist_concatenate: duplicate id')
                  seen_ucat(ids(k)) = 1
                  IF (xs(k) < 1 .or. xs(k) > nlon .or. ys(k) < 1 .or. ys(k) > nlat) &
                     CALL CoLM_stop ('river_hist_concatenate: grid index out of range')
                  grid(xs(k),ys(k)) = vals(k)
               ENDDO
               deallocate (vals)
            ENDIF
            IF (allocated(ids)) deallocate (ids)
            IF (allocated(xs))  deallocate (xs)
            IF (allocated(ys))  deallocate (ys)
            ierr = nf90_close (ncid)
         ENDDO
         IF (count(seen_ucat == 0) > 0) THEN
            write(*,'(A,I0,A,A)') 'ERROR: ', count(seen_ucat == 0), &
               ' unit catchments missing from the shard set for ', trim(names(iv))
            CALL CoLM_stop ('river_hist_concatenate: incomplete coverage')
         ENDIF
         CALL ncio_write_serial (trim(filetmp), trim(names(iv)), grid, &
            'lon_ucat', 'lat_ucat', DEF_REST_CompressLevel)
         CALL copy_var_attrs (trim(names(iv)))
         CALL ncio_put_attr (trim(filetmp), trim(names(iv)), 'missing_value', SPVAL)
      ENDDO
      deallocate (grid)

   END SUBROUTINE rebuild_static_ucat_variables

   !> Reservoir fields are (reservoir_local, time) -- the same shape as a
   !! unit-catchment field, which is why they must be selected by dimension
   !! name rather than by shape. Excluding them from the unit-catchment sweep
   !! is necessary but not sufficient: without this they were written into the
   !! shards and then silently absent from the aggregate, which is worse than
   !! being visibly wrong. Rebuilt on resv_global_index, never on rank order.
   SUBROUTINE rebuild_resv_variables ()

   character(len=256) :: names(256)
   real(r8), allocatable :: col(:), vals(:)
   integer,  allocatable :: ids(:), seen(:)
   integer :: nv, v, nd, e, dimids(NF90_MAX_VAR_DIMS), tdim, rdim
   integer :: nnames, it, k, n, i, vid, totalresv
   character(len=256) :: nm

      ! total reservoir count and the variable list, from shard 0
      nnames = 0
      totalresv = 0
      DO ish = 0, seg_shards(rec_seg(1))-1
         CALL shard_name (rec_seg(1), ish, fname)
         ierr = nf90_open (trim(fname), NF90_NOWRITE, ncid)
         CALL read_int (ncid, 'resv_global_index', ids)
         totalresv = totalresv + size(ids)
         IF (nnames == 0) THEN
            ierr = nf90_inq_dimid (ncid, 'time', tdim)
            e = nf90_inq_dimid (ncid, 'reservoir_local', rdim)
            IF (e /= NF90_NOERR) rdim = -1
            ierr = nf90_inquire (ncid, nVariables=nv)
            DO v = 1, nv
               e = nf90_inquire_variable (ncid, v, name=nm, ndims=nd, dimids=dimids)
               IF (nd /= 2) CYCLE
               IF (dimids(2) /= tdim) CYCLE
               IF (dimids(1) /= rdim) CYCLE
               nnames = nnames + 1
               names(nnames) = nm
            ENDDO
         ENDIF
         IF (allocated(ids)) deallocate (ids)
         ierr = nf90_close (ncid)
      ENDDO

      IF (totalresv == 0) RETURN

      CALL ncio_define_dimension (trim(filetmp), 'reservoir', totalresv)

      ! The reservoir identity itself. It has no time dimension, so it is not
      ! among the fields collected above, and without it the aggregate has a
      ! reservoir axis with no way to say which dam each entry is -- a schema
      ! the 'one' file does not share.
      CALL rebuild_resv_grand_id (totalresv)

      IF (nnames == 0) RETURN
      allocate (col(totalresv))
      allocate (seen(totalresv))

      DO i = 1, nnames
         DO it = 1, ntime
            col = SPVAL
            seen = 0
            DO ish = 0, seg_shards(rec_seg(it))-1
               CALL shard_name (rec_seg(it), ish, fname)
               ierr = nf90_open (trim(fname), NF90_NOWRITE, ncid)
               ierr = nf90_inq_varid (ncid, trim(names(i)), vid)
               IF (ierr /= NF90_NOERR) THEN
                  ierr = nf90_close (ncid); CYCLE
               ENDIF
               CALL read_int (ncid, 'resv_global_index', ids)
               n = size(ids)
               IF (n > 0) THEN
                  CALL read_slice (ncid, vid, n, rec_local(it), vals)
                  DO k = 1, n
                     IF (ids(k) < 1 .or. ids(k) > totalresv) CALL CoLM_stop &
                        ('river_hist_concatenate: reservoir id out of range')
                     IF (seen(ids(k)) /= 0) CALL CoLM_stop &
                        ('river_hist_concatenate: duplicate reservoir id')
                     seen(ids(k)) = 1
                     col(ids(k)) = vals(k)
                  ENDDO
                  deallocate (vals)
               ENDIF
               IF (allocated(ids)) deallocate (ids)
               ierr = nf90_close (ncid)
            ENDDO
            IF (count(seen == 0) > 0) THEN
               write(*,'(A,I0,A,A)') 'ERROR: ', count(seen == 0), &
                  ' reservoirs missing from the shard set for ', trim(names(i))
               CALL CoLM_stop ('river_hist_concatenate: incomplete reservoir coverage')
            ENDIF
            CALL ncio_write_serial_time (trim(filetmp), trim(names(i)), it, col, &
               'reservoir', 'time', DEF_HIST_CompressLevel)
         ENDDO
         CALL copy_var_attrs (trim(names(i)))
      ENDDO

      write(*,'(A,I0,A,I0,A)') '  rebuilt ', nnames, ' reservoir fields over ', &
         totalresv, ' reservoirs'
      deallocate (col, seen)

   END SUBROUTINE rebuild_resv_variables

   !> resv_GRAND_ID(reservoir), reassembled on resv_global_index like the
   !! time-varying reservoir fields.
   SUBROUTINE rebuild_resv_grand_id (totalresv)

   integer, intent(in) :: totalresv
   integer, allocatable :: gid(:), ids(:), seen(:)
   integer :: k, n, vid, e, nwith

      allocate (gid(totalresv)); gid = -9999
      allocate (seen(totalresv)); seen = 0
      nwith = 0

      DO ish = 0, seg_shards(1)-1
         CALL shard_name (1, ish, fname)
         ierr = nf90_open (trim(fname), NF90_NOWRITE, ncid)
         e = nf90_inq_varid (ncid, 'resv_GRAND_ID', vid)
         IF (e /= NF90_NOERR) THEN
            ierr = nf90_close (ncid); CYCLE
         ENDIF
         nwith = nwith + 1
         CALL read_int (ncid, 'resv_global_index', ids)
         n = size(ids)
         IF (n > 0) THEN
            CALL read_int (ncid, 'resv_GRAND_ID', vals_i)
            DO k = 1, n
               IF (ids(k) < 1 .or. ids(k) > totalresv) CALL CoLM_stop &
                  ('river_hist_concatenate: reservoir id out of range')
               gid(ids(k)) = vals_i(k)
               seen(ids(k)) = 1
            ENDDO
            IF (allocated(vals_i)) deallocate (vals_i)
         ENDIF
         IF (allocated(ids)) deallocate (ids)
         ierr = nf90_close (ncid)
      ENDDO

      ! No shard carries the variable at all -- a run without dams, or shards
      ! from a build that did not emit it. There is nothing to rebuild, and
      ! that is not an error. Partial coverage IS an error: it would write a
      ! GRAND ID array silently padded with -9999.
      IF (nwith == 0) THEN
         deallocate (gid, seen)
         RETURN
      ENDIF

      IF (count(seen == 0) == 0) THEN
         CALL ncio_write_serial (trim(filetmp), 'resv_GRAND_ID', gid, 'reservoir')
         CALL ncio_put_attr (trim(filetmp), 'resv_GRAND_ID', 'long_name', &
            'reservoir GRAND ID')
      ELSE
         write(*,'(A,I0,A)') 'ERROR: ', count(seen == 0), &
            ' reservoirs have no GRAND ID in the shard set'
         CALL CoLM_stop ('river_hist_concatenate: incomplete reservoir GRAND ID')
      ENDIF

      deallocate (gid, seen)

   END SUBROUTINE rebuild_resv_grand_id

   !> Bifurcation pathways are keyed strictly on pth_global_id.
   SUBROUTINE rebuild_bif_matrix ()

   real(r8), allocatable :: mat(:,:), part(:,:)
   integer,  allocatable :: ids(:)
   integer :: it, k, vid, n, l
   logical :: present_any

      present_any = .false.
      totalnpthout = 0
      npthlev = 0
      DO ish = 0, seg_shards(rec_seg(1))-1
         CALL shard_name (rec_seg(1), ish, fname)
         ierr = nf90_open (trim(fname), NF90_NOWRITE, ncid)
         ierr = nf90_inq_varid (ncid, 'f_bifflw_lev', vid)
         IF (ierr == NF90_NOERR) THEN
            present_any = .true.
            CALL read_dim (ncid, 'bifurcation_level', npthlev)
            CALL read_int (ncid, 'pth_global_id', ids)
            totalnpthout = totalnpthout + size(ids)
            IF (allocated(ids)) deallocate (ids)
         ENDIF
         ierr = nf90_close (ncid)
      ENDDO

      IF (.not. present_any .or. totalnpthout == 0) RETURN

      allocate (seen_path(totalnpthout))
      allocate (mat(npthlev, totalnpthout))

      CALL ncio_define_dimension (trim(filetmp), 'bifurcation_level', npthlev)
      CALL ncio_define_dimension (trim(filetmp), 'bifurcation_pathway', totalnpthout)

      DO it = 1, ntime
         mat = SPVAL
         seen_path = 0
         DO ish = 0, seg_shards(rec_seg(it))-1
            CALL shard_name (rec_seg(it), ish, fname)
            ierr = nf90_open (trim(fname), NF90_NOWRITE, ncid)
            ierr = nf90_inq_varid (ncid, 'f_bifflw_lev', vid)
            IF (ierr /= NF90_NOERR) THEN
               ierr = nf90_close (ncid); CYCLE
            ENDIF
            CALL read_int (ncid, 'pth_global_id', ids)
            n = size(ids)
            IF (n > 0) THEN
               allocate (part(npthlev, n))
               ierr = nf90_get_var (ncid, vid, part, start=[1,1,rec_local(it)], count=[npthlev,n,1])
               IF (ierr /= NF90_NOERR) CALL CoLM_stop &
                  ('river_hist_concatenate: cannot read f_bifflw_lev')
               DO k = 1, n
                  IF (ids(k) < 1 .or. ids(k) > totalnpthout) &
                     CALL CoLM_stop ('river_hist_concatenate: pathway id out of range')
                  IF (seen_path(ids(k)) /= 0) &
                     CALL CoLM_stop ('river_hist_concatenate: duplicate pathway id')
                  seen_path(ids(k)) = 1
                  DO l = 1, npthlev
                     mat(l, ids(k)) = part(l, k)
                  ENDDO
               ENDDO
               deallocate (part)
            ENDIF
            IF (allocated(ids)) deallocate (ids)
            ierr = nf90_close (ncid)
         ENDDO

         IF (count(seen_path == 0) > 0) THEN
            write(*,'(A,I0,A)') 'ERROR: ', count(seen_path == 0), &
               ' bifurcation pathways missing from the shard set'
            CALL CoLM_stop ('river_hist_concatenate: incomplete pathway coverage')
         ENDIF

         CALL ncio_write_serial_time (trim(filetmp), 'f_bifflw_lev', it, mat, &
            'bifurcation_level', 'bifurcation_pathway', 'time', DEF_HIST_CompressLevel)
      ENDDO

      CALL copy_var_attrs ('f_bifflw_lev')
      write(*,'(A,I0,A)') '  rebuilt bifurcation matrix over ', totalnpthout, ' pathways'

   END SUBROUTINE rebuild_bif_matrix

   !> Only a fully verified file is allowed to take the target name.
   SUBROUTINE verify_and_promote ()

   integer :: cmd, u
   logical :: fexists

      ierr = nf90_open (trim(filetmp), NF90_NOWRITE, ncid)
      IF (ierr /= NF90_NOERR) CALL CoLM_stop &
         ('river_hist_concatenate: the aggregate did not open for verification')
      ierr = nf90_close (ncid)

      inquire (file=trim(filedone), exist=fexists)
      IF (fexists) THEN
         OPEN (newunit=u, file=trim(filedone), status='old'); CLOSE (u, status='delete')
      ENDIF

      CALL rename_file (trim(filetmp), trim(filetarget))

      OPEN (newunit=u, file=trim(filedone), status='replace', action='write')
      write(u,'(A)')      'river_hist_concatenate complete'
      write(u,'(A,A)')    'target      : ', trim(filetarget)
      write(u,'(A,A)')    'case        : ', trim(ref_case)
      write(u,'(A,A)')    'period      : ', trim(ref_period)
      write(u,'(A,I0)')   'shards      : ', shard_count
      write(u,'(A,I0)')   'unitcat     : ', totalnumucat
      write(u,'(A,I0)')   'pathways    : ', totalnpthout
      write(u,'(A)')      'shards were NOT deleted; remove them only after checking this file'
      CLOSE (u)

      ! The pending note has served its purpose; leaving it would say the
      ! output still needs aggregating when it no longer does.
      inquire (file=trim(filetarget)//'.pending', exist=fexists)
      IF (fexists) THEN
         OPEN (newunit=u, file=trim(filetarget)//'.pending', status='old')
         CLOSE (u, status='delete')
      ENDIF

      write(*,'(A,A)') 'river_hist_concatenate: wrote ', trim(filetarget)
      write(*,'(A,A)') '                        marker ', trim(filedone)

   END SUBROUTINE verify_and_promote

   ! ================== small helpers ==================

   SUBROUTINE shard_name (iseg, idx, f)
   integer, intent(in) :: iseg, idx
   character(len=256), intent(out) :: f
   integer :: i
   character(len=8) :: c
   character(len=80) :: sg
      write(c,'(I5.5)') idx
      sg = ''
      IF (iseg >= 1 .and. iseg <= nseg) sg = '_' // trim(seg_id(iseg))
      i = len_trim(filetarget)
      IF (i > 3 .and. filetarget(max(i-2,1):i) == '.nc') THEN
         f = filetarget(1:i-3) // trim(sg) // '_shard' // trim(c) // '.nc'
      ELSE
         f = trim(filetarget) // trim(sg) // '_shard' // trim(c) // '.nc'
      ENDIF
   END SUBROUTINE shard_name

   SUBROUTINE identity_error (f, what)
   character(len=*), intent(in) :: f, what
      write(*,'(A,A,A,A)') 'ERROR: shard ', trim(f), ' disagrees on ', trim(what)
      CALL CoLM_stop ('river_hist_concatenate: shard identity mismatch')
   END SUBROUTINE identity_error

   SUBROUTINE get_att_int (nc, name, v)
   integer, intent(in) :: nc
   character(len=*), intent(in) :: name
   integer, intent(out) :: v
   integer :: e
      e = nf90_get_att (nc, NF90_GLOBAL, trim(name), v)
      IF (e /= NF90_NOERR) CALL CoLM_stop &
         ('river_hist_concatenate: missing shard attribute '//trim(name))
   END SUBROUTINE get_att_int

   SUBROUTINE get_att_str (nc, name, v)
   integer, intent(in) :: nc
   character(len=*), intent(in) :: name
   character(len=256), intent(out) :: v
   integer :: e
      v = ''
      e = nf90_get_att (nc, NF90_GLOBAL, trim(name), v)
      IF (e /= NF90_NOERR) CALL CoLM_stop &
         ('river_hist_concatenate: missing shard attribute '//trim(name))
   END SUBROUTINE get_att_str

   SUBROUTINE read_dim (nc, name, n)
   integer, intent(in) :: nc
   character(len=*), intent(in) :: name
   integer, intent(out) :: n
   integer :: d, e
      n = 0
      e = nf90_inq_dimid (nc, trim(name), d)
      IF (e == NF90_NOERR) e = nf90_inquire_dimension (nc, d, len=n)
   END SUBROUTINE read_dim

   SUBROUTINE read_int (nc, name, a)
   integer, intent(in) :: nc
   character(len=*), intent(in) :: name
   integer, allocatable, intent(out) :: a(:)
   integer :: v, e, n
      n = 0
      e = nf90_inq_varid (nc, trim(name), v)
      IF (e == NF90_NOERR) e = nf90_inquire_dimension (nc, &
         var_first_dim(nc, v), len=n)
      allocate (a(max(n,0)))
      IF (n > 0) e = nf90_get_var (nc, v, a)
   END SUBROUTINE read_int

   SUBROUTINE read_real (nc, name, a)
   integer, intent(in) :: nc
   character(len=*), intent(in) :: name
   real(r8), allocatable, intent(out) :: a(:)
   integer :: v, e, n
      n = 0
      e = nf90_inq_varid (nc, trim(name), v)
      IF (e == NF90_NOERR) e = nf90_inquire_dimension (nc, &
         var_first_dim(nc, v), len=n)
      allocate (a(max(n,1)))
      IF (n > 0) e = nf90_get_var (nc, v, a)
   END SUBROUTINE read_real

   integer FUNCTION var_first_dim (nc, v)
   integer, intent(in) :: nc, v
   integer :: dimids(NF90_MAX_VAR_DIMS), nd, e
      e = nf90_inquire_variable (nc, v, ndims=nd, dimids=dimids)
      var_first_dim = dimids(1)
   END FUNCTION var_first_dim

   SUBROUTINE read_slice (nc, vid, n, it, vals)
   integer, intent(in) :: nc, vid, n, it
   real(r8), allocatable, intent(out) :: vals(:)
   integer :: e
      allocate (vals(n))
      e = nf90_get_var (nc, vid, vals, start=[1,it], count=[n,1])
      IF (e /= NF90_NOERR) CALL CoLM_stop &
         ('river_hist_concatenate: cannot read a shard slice')
   END SUBROUTINE read_slice

   !> Carry long_name / units / missing_value across from shard 0.
   SUBROUTINE copy_var_attrs (varname)
   character(len=*), intent(in) :: varname
   integer :: nc, vid, e
   character(len=256) :: c
   real(r8) :: mv
      CALL shard_name (1, 0, fname)
      e = nf90_open (trim(fname), NF90_NOWRITE, nc)
      IF (e /= NF90_NOERR) RETURN
      e = nf90_inq_varid (nc, trim(varname), vid)
      IF (e == NF90_NOERR) THEN
         c = ''
         IF (nf90_get_att (nc, vid, 'long_name', c) == NF90_NOERR) &
            CALL ncio_put_attr (trim(filetmp), trim(varname), 'long_name', trim(c))
         c = ''
         IF (nf90_get_att (nc, vid, 'units', c) == NF90_NOERR) &
            CALL ncio_put_attr (trim(filetmp), trim(varname), 'units', trim(c))
         IF (nf90_get_att (nc, vid, 'missing_value', mv) == NF90_NOERR) &
            CALL ncio_put_attr (trim(filetmp), trim(varname), 'missing_value', mv)
      ENDIF
      e = nf90_close (nc)
   END SUBROUTINE copy_var_attrs

   !> Every time-varying field on the shard, i.e. what has to be rebuilt.
   SUBROUTINE collect_varnames (nc, names, n)
   integer, intent(in) :: nc
   character(len=256), allocatable, intent(out) :: names(:)
   integer, intent(out) :: n
   integer :: nv, v, nd, e, dimids(NF90_MAX_VAR_DIMS), tdim, udim
   character(len=256) :: nm
   character(len=256) :: tmp(512)
      n = 0
      e = nf90_inq_dimid (nc, 'time', tdim)
      e = nf90_inquire (nc, nVariables=nv)
      ! A unit-catchment field is (unitcat_local, time). Matching on "two
      ! dimensions, second is time" is NOT enough: reservoir fields are
      ! (reservoir_local, time) and would be swept up here, then rebuilt
      ! against ucat_ucid/x_ucat/y_ucat, whose length is unitcat_local. The
      ! first dimension must be unitcat_local explicitly.
      e = nf90_inq_dimid (nc, 'unitcat_local', udim)
      IF (e /= NF90_NOERR) udim = -1
      DO v = 1, nv
         e = nf90_inquire_variable (nc, v, name=nm, ndims=nd, dimids=dimids)
         IF (nd /= 2) CYCLE
         IF (dimids(2) /= tdim) CYCLE
         IF (dimids(1) /= udim) CYCLE
         IF (trim(nm) == 'time') CYCLE
         n = n + 1
         tmp(n) = nm
      ENDDO
      allocate (names(max(n,1)))
      IF (n > 0) names(1:n) = tmp(1:n)
   END SUBROUTINE collect_varnames

   !> Move src onto dst, or stop. Nothing downstream may treat a failed
   !! promotion as success: the .complete marker written afterwards is what
   !! tells an operator the aggregate is trustworthy.
   SUBROUTINE rename_file (src, dst)
   character(len=*), intent(in) :: src, dst
   integer :: u, estat, cstat
   logical :: fexists
   character(len=256) :: cmsg

      ! Both paths come from the command line. Unquoted, a path containing a
      ! space was split into separate mv arguments -- observed moving the file
      ! to a truncated name, or failing outright. Single quotes handle every
      ! shell metacharacter except a single quote itself, which cannot be
      ! escaped inside them, so that case is refused rather than mangled.
      ! Defence in depth; the real gate is reject_unquotable at startup.
      CALL reject_unquotable (src, 'source')
      CALL reject_unquotable (dst, 'destination')

      inquire (file=trim(dst), exist=fexists)
      IF (fexists) THEN
         OPEN (newunit=u, file=trim(dst), status='old'); CLOSE (u, status='delete')
      ENDIF

      ! execute_command_line, not system: it is standard Fortran and it reports
      ! the exit status. system() gave none, so mv could fail and the program
      ! carried on to write the success marker.
      cmsg = ''
      CALL execute_command_line ("mv '" // trim(src) // "' '" // trim(dst) // "'", &
         wait = .true., exitstat = estat, cmdstat = cstat, cmdmsg = cmsg)

      IF (cstat /= 0) THEN
         write(*,'(A,I0,2A)') 'ERROR: could not run mv (cmdstat=', cstat, '): ', trim(cmsg)
         CALL CoLM_stop ('river_hist_concatenate: rename failed')
      ENDIF
      IF (estat /= 0) THEN
         write(*,'(A,I0)') 'ERROR: mv exited with status ', estat
         write(*,'(4A)') '       ', trim(src), ' -> ', trim(dst)
         CALL CoLM_stop ('river_hist_concatenate: rename failed')
      ENDIF

      ! A zero exit status is not proof the file arrived. Check the filesystem.
      inquire (file=trim(dst), exist=fexists)
      IF (.not. fexists) THEN
         write(*,'(2A)') 'ERROR: rename reported success but the target is absent: ', trim(dst)
         CALL CoLM_stop ('river_hist_concatenate: rename failed')
      ENDIF
      inquire (file=trim(src), exist=fexists)
      IF (fexists) THEN
         write(*,'(2A)') 'ERROR: rename reported success but the source remains: ', trim(src)
         CALL CoLM_stop ('river_hist_concatenate: rename failed')
      ENDIF

   END SUBROUTINE rename_file

END PROGRAM river_hist_concatenate
! ----------------------------------------------------------------------
! EOP
