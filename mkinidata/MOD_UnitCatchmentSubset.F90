MODULE MOD_UnitCatchmentSubset

!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!    Serial NetCDF work of the regional unit-catchment network.
!
!    Given the runoff-input grid cells that the land domain touches, keep the river
!    systems (everything that drains to one river mouth or inland depression) that
!    receive runoff from those cells and write them to a smaller file with the layout
!    of the source file.  Kept cells retain their order, so seq_next still points
!    downstream to a larger index; indices held by seq, seq_next, seq_upst, the
!    bifurcation pathways and the dam table are renumbered.  Every other variable is
!    cut along the unit-catchment, pathway or dam dimension without interpretation,
!    so variables added to the network file later need no change here.
!
!    "Receives runoff" uses the input matrix of the network file (inpmat_x, inpmat_y),
!    the same matrix the model uses to route runoff into unit catchments.
!
!    Bifurcation pathways may join different river systems.  With bif_closure the
!    selection grows until every pathway that touches the selection lies wholly in it;
!    without it, pathways with one end outside are dropped.
!
!    The new file also carries seq_src_index (position of every kept unit catchment in
!    the source network) and the global attribute source_nseqmax; tools that must
!    translate indices of the source network, such as the reservoir catalogue, use them.
!
!    Errors stop the program: the only caller is a preprocessing step.
!-------------------------------------------------------------------------------------

   USE netcdf
   IMPLICIT NONE
   PRIVATE

   integer, parameter :: i8 = selected_int_kind(18)
   integer, parameter :: r8k = selected_real_kind(15)

   character(len=*), parameter :: SEQ_DIM  = 'nseqmax'
   character(len=*), parameter :: PATH_DIM = 'npthout'
   character(len=*), parameter :: DAM_DIM  = 'dam_ndams'

   PUBLIC :: unitcatchment_grid_size
   PUBLIC :: unitcatchment_subset_write

CONTAINS

   ! ---------
   SUBROUTINE fail (text)

   IMPLICIT NONE
   character(len=*), intent(in) :: text

      write(*,'(A,A)') 'ERROR (regional unit-catchment network): ', trim(text)
      error stop 1

   END SUBROUTINE fail

   ! ---------
   SUBROUTINE check (status, what)

   IMPLICIT NONE
   integer,          intent(in) :: status
   character(len=*), intent(in) :: what

      IF (status /= nf90_noerr) CALL fail (trim(what) // ': ' // trim(nf90_strerror(status)))

   END SUBROUTINE check

   ! ---------
   FUNCTION dim_length (ncid, name, required) RESULT (n)

   IMPLICIT NONE
   integer,          intent(in) :: ncid
   character(len=*), intent(in) :: name
   logical,          intent(in) :: required
   integer :: n

   integer :: dimid, status

      n = -1
      status = nf90_inq_dimid (ncid, name, dimid)
      IF (status == nf90_noerr) THEN
         CALL check (nf90_inquire_dimension (ncid, dimid, len = n), 'dimension ' // name)
      ELSEIF (required) THEN
         CALL fail ('dimension ' // name // ' not found')
      ENDIF

   END FUNCTION dim_length

   ! ---------
   SUBROUTINE unitcatchment_grid_size (file, nlon, nlat)

   ! Size of the runoff-input grid (the file's lon and lat coordinates).

   IMPLICIT NONE
   character(len=*), intent(in)  :: file
   integer,          intent(out) :: nlon, nlat

   integer :: ncid

      CALL check (nf90_open (trim(file), nf90_nowrite, ncid), 'open ' // trim(file))
      nlon = dim_length (ncid, 'nx', .true.)
      nlat = dim_length (ncid, 'ny', .true.)
      CALL check (nf90_close (ncid), 'close ' // trim(file))

   END SUBROUTINE unitcatchment_grid_size

   ! ---------
   SUBROUTINE read_int_1d (ncid, name, values)

   IMPLICIT NONE
   integer,          intent(in) :: ncid
   character(len=*), intent(in) :: name
   integer, allocatable, intent(out) :: values(:)

   integer :: varid, ndims, dimids(nf90_max_var_dims), n

      CALL check (nf90_inq_varid (ncid, name, varid), 'variable ' // name)
      CALL check (nf90_inquire_variable (ncid, varid, ndims = ndims, dimids = dimids), 'inquire ' // name)
      IF (ndims /= 1) CALL fail (name // ' must be one-dimensional')
      CALL check (nf90_inquire_dimension (ncid, dimids(1), len = n), 'length of ' // name)
      allocate (values (n))
      CALL check (nf90_get_var (ncid, varid, values), 'read ' // name)

   END SUBROUTINE read_int_1d

   ! ---------
   SUBROUTINE read_int_2d (ncid, name, values)

   IMPLICIT NONE
   integer,          intent(in) :: ncid
   character(len=*), intent(in) :: name
   integer, allocatable, intent(out) :: values(:,:)

   integer :: varid, ndims, dimids(nf90_max_var_dims), n1, n2

      CALL check (nf90_inq_varid (ncid, name, varid), 'variable ' // name)
      CALL check (nf90_inquire_variable (ncid, varid, ndims = ndims, dimids = dimids), 'inquire ' // name)
      IF (ndims /= 2) CALL fail (name // ' must be two-dimensional')
      CALL check (nf90_inquire_dimension (ncid, dimids(1), len = n1), 'length of ' // name)
      CALL check (nf90_inquire_dimension (ncid, dimids(2), len = n2), 'length of ' // name)
      allocate (values (n1,n2))
      CALL check (nf90_get_var (ncid, varid, values), 'read ' // name)

   END SUBROUTINE read_int_2d

   ! ---------
   SUBROUTINE copy_variable (ncin, ncout, varid, name, xtype, ndims, dimids, selected, isel, renumber, new_index)

   ! Copy one variable, keeping only the entries isel(:) along dimension `selected`
   ! (0: no selection) and renumbering positive entries when `renumber` is set.

   IMPLICIT NONE
   integer,          intent(in) :: ncin, ncout, varid, xtype, ndims, dimids(:), selected
   character(len=*), intent(in) :: name
   integer,          intent(in) :: isel(:)
   logical,          intent(in) :: renumber
   integer,          intent(in) :: new_index(:)

   integer :: n(2), k, i
   integer, allocatable :: idx1(:), idx2(:)
   integer,  allocatable :: i1(:), i2(:,:), j1(:), j2(:,:)
   real(r8k), allocatable :: d1(:), d2(:,:), e1(:), e2(:,:)
   integer :: i0
   real(r8k) :: d0

      IF (ndims > 2) CALL fail ('variable ' // name // ' has more than two dimensions')

      n = 1
      DO k = 1, ndims
         CALL check (nf90_inquire_dimension (ncin, dimids(k), len = n(k)), 'length of ' // name)
      ENDDO

      ! index vectors: identity except along the selected dimension
      allocate (idx1 (n(1)));  idx1 = (/ (i, i = 1, n(1)) /)
      allocate (idx2 (n(2)));  idx2 = (/ (i, i = 1, n(2)) /)
      IF (selected == 1) THEN
         deallocate (idx1);  allocate (idx1 (size(isel)));  idx1 = isel
      ELSEIF (selected == 2) THEN
         deallocate (idx2);  allocate (idx2 (size(isel)));  idx2 = isel
      ENDIF

      IF (renumber .and. xtype /= nf90_int) CALL fail ('index variable ' // name // ' must be integer')

      SELECT CASE (xtype)

      CASE (nf90_int)
         SELECT CASE (ndims)
         CASE (0)
            CALL check (nf90_get_var (ncin, varid, i0), 'read ' // name)
            CALL check (nf90_put_var (ncout, varid, i0), 'write ' // name)
         CASE (1)
            allocate (i1 (n(1)))
            CALL check (nf90_get_var (ncin, varid, i1), 'read ' // name)
            allocate (j1 (size(idx1)));  j1 = i1(idx1)
            IF (renumber) CALL renumber_flat (j1, size(j1), new_index, name)
            CALL check (nf90_put_var (ncout, varid, j1), 'write ' // name)
         CASE (2)
            allocate (i2 (n(1),n(2)))
            CALL check (nf90_get_var (ncin, varid, i2), 'read ' // name)
            allocate (j2 (size(idx1),size(idx2)));  j2 = i2(idx1,idx2)
            IF (renumber) CALL renumber_flat (j2, size(j2), new_index, name)
            CALL check (nf90_put_var (ncout, varid, j2), 'write ' // name)
         END SELECT

      CASE (nf90_double)
         SELECT CASE (ndims)
         CASE (0)
            CALL check (nf90_get_var (ncin, varid, d0), 'read ' // name)
            CALL check (nf90_put_var (ncout, varid, d0), 'write ' // name)
         CASE (1)
            allocate (d1 (n(1)))
            CALL check (nf90_get_var (ncin, varid, d1), 'read ' // name)
            allocate (e1 (size(idx1)));  e1 = d1(idx1)
            CALL check (nf90_put_var (ncout, varid, e1), 'write ' // name)
         CASE (2)
            allocate (d2 (n(1),n(2)))
            CALL check (nf90_get_var (ncin, varid, d2), 'read ' // name)
            allocate (e2 (size(idx1),size(idx2)));  e2 = d2(idx1,idx2)
            CALL check (nf90_put_var (ncout, varid, e2), 'write ' // name)
         END SELECT

      CASE (nf90_char)
         ! The first dimension of a character variable is the string length; the
         ! strings themselves can be cut, their characters cannot.
         IF (selected == 1) CALL fail ('character variable ' // name // ' is cut along its string length')
         SELECT CASE (ndims)
         CASE (1)
            CALL copy_strings (ncin, ncout, varid, name, n(1), 1, (/ 1 /))
         CASE (2)
            CALL copy_strings (ncin, ncout, varid, name, n(1), n(2), idx2)
         CASE DEFAULT
            CALL fail ('scalar character variable ' // name // ' is not supported')
         END SELECT

      CASE DEFAULT
         CALL fail ('variable ' // name // ' has an unsupported data type')
      END SELECT

   END SUBROUTINE copy_variable

   ! ---------
   SUBROUTINE copy_strings (ncin, ncout, varid, name, nlen, nstr, idx)

   IMPLICIT NONE
   integer,          intent(in) :: ncin, ncout, varid, nlen, nstr, idx(:)
   character(len=*), intent(in) :: name

   character(len=nlen), allocatable :: text(:), kept(:)

      allocate (text (nstr))
      CALL check (nf90_get_var (ncin, varid, text), 'read ' // name)
      allocate (kept (size(idx)))
      kept = text(idx)
      CALL check (nf90_put_var (ncout, varid, kept), 'write ' // name)

   END SUBROUTINE copy_strings

   ! ---------
   FUNCTION holds_unitcatchment_numbers (name) RESULT (yes)

   IMPLICIT NONE
   character(len=*), intent(in) :: name
   logical :: yes

      SELECT CASE (trim(name))
      CASE ('seq', 'seq_next', 'seq_upst', 'bifurcation_upst', 'bifurcation_down', 'dam_seq')
         yes = .true.
      CASE DEFAULT
         yes = .false.
      END SELECT

   END FUNCTION holds_unitcatchment_numbers

   ! ---------
   SUBROUTINE renumber_flat (values, n, new_index, name)

   IMPLICIT NONE
   integer,          intent(in)    :: n
   integer,          intent(inout) :: values(n)
   integer,          intent(in)    :: new_index(:)
   character(len=*), intent(in)    :: name

   integer :: i

      DO i = 1, n
         IF (values(i) > 0) THEN
            IF (values(i) > size(new_index)) CALL fail (name // ' refers beyond the network')
            IF (new_index(values(i)) <= 0) CALL fail (name // ' refers to a unit catchment outside the selection')
            values(i) = new_index(values(i))
         ENDIF
      ENDDO

   END SUBROUTINE renumber_flat

   ! ---------
   SUBROUTINE unitcatchment_subset_write (file_in, file_out, nlon, nlat, cell_touched, bif_closure, &
         nkeep, nsystem)

   IMPLICIT NONE
   character(len=*), intent(in)  :: file_in, file_out
   integer,          intent(in)  :: nlon, nlat
   logical,          intent(in)  :: cell_touched(nlon,nlat)   ! runoff-input grid cells of the land domain
   logical,          intent(in)  :: bif_closure
   integer,          intent(out) :: nkeep, nsystem

   integer :: ncin, ncout, ndim, nvar, nglobal, unlim
   integer :: nseq, npth, ndam
   integer :: dimid_seq, dimid_path, dimid_dam
   integer :: i, j, k, p, varid, varid_new, natt, xtype, ndims, iatt, status
   integer :: dimids(nf90_max_var_dims)
   integer, allocatable :: seq_next(:), inpmat_x(:,:), inpmat_y(:,:), bif_up(:), bif_dn(:), dam_seq(:)
   integer, allocatable :: mouth(:), new_index(:), isel_seq(:), isel_path(:), isel_dam(:)
   logical, allocatable :: receives(:), system_kept(:), keep(:), keep_path(:), keep_dam(:)
   logical :: changed
   integer :: nseqriv_old, nseqriv_new
   integer(i8) :: attr8 = 0
   character(len=nf90_max_name) :: name, attname, dimname
   logical :: renumber
   integer :: selected

      CALL check (nf90_open (trim(file_in), nf90_nowrite, ncin), 'open ' // trim(file_in))

      nseq = dim_length (ncin, SEQ_DIM, .true.)
      npth = dim_length (ncin, PATH_DIM, .false.)
      ndam = dim_length (ncin, DAM_DIM, .false.)

      ! ---- which unit catchments receive runoff ----
      CALL read_int_1d (ncin, 'seq_next', seq_next)
      CALL read_int_2d (ncin, 'inpmat_x', inpmat_x)
      CALL read_int_2d (ncin, 'inpmat_y', inpmat_y)
      IF (size(seq_next) /= nseq .or. size(inpmat_x,2) /= nseq .or. any(shape(inpmat_x) /= shape(inpmat_y))) &
         CALL fail ('seq_next and the input matrix do not agree on the number of unit catchments')

      allocate (mouth (nseq))
      DO i = 1, nseq
         mouth(i) = i
         IF (seq_next(i) > 0) THEN
            IF (seq_next(i) <= i) CALL fail ('seq_next is not ordered from upstream to downstream')
         ENDIF
      ENDDO
      DO i = nseq, 1, -1
         IF (seq_next(i) > 0) mouth(i) = mouth(seq_next(i))
      ENDDO

      allocate (receives (nseq));  receives = .false.
      DO i = 1, nseq
         DO j = 1, size(inpmat_x,1)
            IF (inpmat_x(j,i) >= 1 .and. inpmat_x(j,i) <= nlon .and. &
                inpmat_y(j,i) >= 1 .and. inpmat_y(j,i) <= nlat) THEN
               IF (cell_touched(inpmat_x(j,i), inpmat_y(j,i))) receives(i) = .true.
            ENDIF
         ENDDO
      ENDDO
      IF (.not. any(receives)) &
         CALL fail ('no unit catchment receives runoff from the land domain; check that the network matches the domain')

      allocate (system_kept (nseq));  system_kept = .false.
      DO i = 1, nseq
         IF (receives(i)) system_kept(mouth(i)) = .true.
      ENDDO

      IF (npth > 0) THEN
         CALL read_int_1d (ncin, 'bifurcation_upst', bif_up)
         CALL read_int_1d (ncin, 'bifurcation_down', bif_dn)
         IF (size(bif_up) /= npth .or. size(bif_dn) /= npth) CALL fail ('bifurcation arrays have the wrong length')
         IF (any(bif_up < 1 .or. bif_up > nseq .or. bif_dn < 1 .or. bif_dn > nseq)) &
            CALL fail ('a bifurcation pathway refers to a unit catchment outside the network')
         IF (bif_closure) THEN
            changed = .true.
            DO WHILE (changed)
               changed = .false.
               DO p = 1, npth
                  i = mouth(bif_up(p));  j = mouth(bif_dn(p))
                  IF (system_kept(i) .neqv. system_kept(j)) THEN
                     system_kept(i) = .true.;  system_kept(j) = .true.
                     changed = .true.
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
      ENDIF
      IF (ndam > 0) CALL read_int_1d (ncin, 'dam_seq', dam_seq)

      allocate (keep (nseq))
      DO i = 1, nseq
         keep(i) = system_kept(mouth(i))
      ENDDO
      nkeep = count(keep)
      nsystem = count(system_kept)

      allocate (new_index (nseq));  new_index = 0
      allocate (isel_seq (nkeep))
      k = 0
      DO i = 1, nseq
         IF (keep(i)) THEN
            k = k + 1
            new_index(i) = k
            isel_seq(k) = i
         ENDIF
      ENDDO

      IF (npth > 0) THEN
         allocate (keep_path (npth))
         keep_path = keep(bif_up) .and. keep(bif_dn)
         allocate (isel_path (count(keep_path)))
         k = 0
         DO p = 1, npth
            IF (keep_path(p)) THEN
               k = k + 1;  isel_path(k) = p
            ENDIF
         ENDDO
      ELSE
         allocate (isel_path (0))
      ENDIF

      IF (ndam > 0) THEN
         allocate (keep_dam (ndam))
         DO k = 1, ndam
            keep_dam(k) = .false.
            IF (dam_seq(k) >= 1 .and. dam_seq(k) <= nseq) keep_dam(k) = keep(dam_seq(k))
         ENDDO
         allocate (isel_dam (count(keep_dam)))
         j = 0
         DO k = 1, ndam
            IF (keep_dam(k)) THEN
               j = j + 1;  isel_dam(j) = k
            ENDIF
         ENDDO
      ELSE
         allocate (isel_dam (0))
      ENDIF

      ! ---- write the file ----
      CALL check (nf90_create (trim(file_out), ior(nf90_clobber, nf90_netcdf4), ncout), 'create ' // trim(file_out))
      CALL check (nf90_inquire (ncin, ndim, nvar, nglobal, unlim), 'inquire ' // trim(file_in))
      IF (unlim /= -1) CALL fail ('the network file has an unlimited dimension, which is not expected')

      CALL check (nf90_inq_dimid (ncin, SEQ_DIM, dimid_seq), 'dimension ' // SEQ_DIM)
      dimid_path = -1;  dimid_dam = -1
      IF (npth >= 0) status = nf90_inq_dimid (ncin, PATH_DIM, dimid_path)
      IF (ndam >= 0) status = nf90_inq_dimid (ncin, DAM_DIM,  dimid_dam)

      DO i = 1, ndim
         CALL check (nf90_inquire_dimension (ncin, i, name = dimname, len = k), 'dimension')
         IF (i == dimid_seq)  k = nkeep
         IF (i == dimid_path) k = size(isel_path)
         IF (i == dimid_dam)  k = size(isel_dam)
         CALL check (nf90_def_dim (ncout, trim(dimname), k, j), 'define dimension ' // trim(dimname))
         IF (j /= i) CALL fail ('dimension order changed while copying')
      ENDDO

      DO iatt = 1, nglobal
         CALL check (nf90_inq_attname (ncin, nf90_global, iatt, attname), 'global attribute')
         CALL check (nf90_copy_att (ncin, nf90_global, trim(attname), ncout, nf90_global), 'copy ' // trim(attname))
      ENDDO

      status = nf90_get_att (ncin, nf90_global, 'nseqriv', attr8)
      nseqriv_old = merge(int(attr8), nseq, status == nf90_noerr)
      nseqriv_new = count(keep(1:min(nseqriv_old, nseq)))

      CALL put_count_attr ('nseqall', nkeep)
      CALL put_count_attr ('nseqmax', nkeep)
      CALL put_count_attr ('nseqriv', nseqriv_new)
      IF (npth >= 0) CALL put_count_attr ('npthout', size(isel_path))
      IF (ndam >= 0) CALL put_count_attr ('dam_ndams', size(isel_dam))
      CALL check (nf90_put_att (ncout, nf90_global, 'source_nseqmax', int(nseq, i8)), 'source_nseqmax')
      CALL check (nf90_put_att (ncout, nf90_global, 'subset_bif_mode', &
         trim(merge('closure', 'drop   ', bif_closure))), 'subset_bif_mode')
      CALL check (nf90_put_att (ncout, nf90_global, 'subset_note', &
         'regional subset written by mksrfdata (DEF_UnitCatchment_regional)'), 'subset_note')

      DO varid = 1, nvar
         CALL check (nf90_inquire_variable (ncin, varid, name = name, xtype = xtype, ndims = ndims, &
            dimids = dimids, natts = natt), 'inquire variable')
         CALL check (nf90_def_var (ncout, trim(name), xtype, dimids(1:ndims), varid_new), 'define ' // trim(name))
         IF (varid_new /= varid) CALL fail ('variable order changed while copying')
         DO iatt = 1, natt
            CALL check (nf90_inq_attname (ncin, varid, iatt, attname), 'attribute of ' // trim(name))
            CALL check (nf90_copy_att (ncin, varid, trim(attname), ncout, varid), &
               'copy attribute of ' // trim(name))
         ENDDO
      ENDDO
      CALL check (nf90_def_var (ncout, 'seq_src_index', nf90_int, (/ dimid_seq /), varid_new), 'define seq_src_index')
      CALL check (nf90_put_att (ncout, varid_new, 'long_name', &
         'index of this unit catchment in the source network'), 'seq_src_index attribute')
      CALL check (nf90_enddef (ncout), 'end definition')

      DO varid = 1, nvar
         CALL check (nf90_inquire_variable (ncin, varid, name = name, xtype = xtype, ndims = ndims, &
            dimids = dimids), 'inquire variable')

         selected = 0
         DO k = 1, ndims
            IF (dimids(k) == dimid_seq .or. dimids(k) == dimid_path .or. dimids(k) == dimid_dam) THEN
               IF (selected /= 0) CALL fail ('variable ' // trim(name) // ' has two cut dimensions')
               selected = k
            ENDIF
         ENDDO

         renumber = holds_unitcatchment_numbers (name)
         IF (selected == 0) THEN
            CALL copy_variable (ncin, ncout, varid, trim(name), xtype, ndims, dimids, 0, (/ 0 /), renumber, new_index)
         ELSEIF (dimids(selected) == dimid_seq) THEN
            CALL copy_variable (ncin, ncout, varid, trim(name), xtype, ndims, dimids, selected, isel_seq, &
               renumber, new_index)
         ELSEIF (dimids(selected) == dimid_path) THEN
            CALL copy_variable (ncin, ncout, varid, trim(name), xtype, ndims, dimids, selected, isel_path, &
               renumber, new_index)
         ELSE
            CALL copy_variable (ncin, ncout, varid, trim(name), xtype, ndims, dimids, selected, isel_dam, &
               renumber, new_index)
         ENDIF
      ENDDO

      CALL check (nf90_inq_varid (ncout, 'seq_src_index', varid_new), 'seq_src_index')
      CALL check (nf90_put_var (ncout, varid_new, isel_seq), 'write seq_src_index')

      CALL check (nf90_close (ncin),  'close ' // trim(file_in))
      CALL check (nf90_close (ncout), 'close ' // trim(file_out))

   CONTAINS

      SUBROUTINE put_count_attr (attr, value)
      IMPLICIT NONE
      character(len=*), intent(in) :: attr
      integer,          intent(in) :: value
         CALL check (nf90_put_att (ncout, nf90_global, attr, int(value, i8)), 'attribute ' // attr)
      END SUBROUTINE put_count_attr

   END SUBROUTINE unitcatchment_subset_write

END MODULE MOD_UnitCatchmentSubset
