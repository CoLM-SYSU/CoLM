#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_SoilInit

   USE MOD_Precision
   USE MOD_Namelist, only: DEF_file_SoilInit
   USE MOD_MT_Namelist, only: DEF_MT_USE_SOIL_INIT, DEF_MT_SOIL_INIT_FILE, DEF_MT_SOIL_INIT_VARS
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_SpatialMapping
   USE MOD_LandPatch
   USE MOD_NetCDFSerial, only: ncio_read_bcast_serial, ncio_var_exist
   USE MOD_NetCDFBlock, only: ncio_read_block_time
   USE MOD_MT_Defs, only: nmts, materials, delta_to_R, material_is_isotope, &
      mt_delta_sanity_max, mt_uses_generic_water
   USE MOD_MT_Vars, only: mt_wliq_soisno, mt_wice_soisno
   USE MOD_MT_Isotope_Registry, only: isotope_default_soil_init_varname
   USE MOD_MT_Providers, only: ensure_isotope_physics_registered

   IMPLICIT NONE

   PUBLIC :: material_soil_init_from_file

CONTAINS

   SUBROUTINE material_soil_init_from_file (init_month, numpatch, maxsnl, nl_soil, &
      wliq_soisno, wice_soisno)

      IMPLICIT NONE
      integer,  intent(in) :: init_month, numpatch, maxsnl, nl_soil
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil, numpatch)

      type(grid_type) :: gsoil
      type(spatial_mapping_type) :: msoil2p
      type(block_data_real8_3d) :: delta_grid
      real(r8), allocatable :: soil_z(:)
      real(r8), allocatable :: delta_patch(:,:)
      character(len=256) :: soil_file, varname
      logical :: file_ok, has_soildepth, has_var
      integer :: month, nl_soil_ini, ncopy
      integer :: imt, ip, k
      real(r8) :: delta_layer, R_layer

      IF (.not. DEF_MT_USE_SOIL_INIT) RETURN
      IF (nmts <= 0) RETURN
      ! This routine contains MPI collectives inside NetCDF/grid mapping.
      ! IO/master ranks may not own patch material arrays, but they must
      ! still participate; only the worker-side storage overwrite is gated.

      soil_file = adjustl(DEF_MT_SOIL_INIT_FILE)
      IF (len_trim(soil_file) == 0 .or. trim(soil_file) == 'null') THEN
         soil_file = adjustl(DEF_file_SoilInit)
      ENDIF
      IF (len_trim(soil_file) == 0 .or. trim(soil_file) == 'null') RETURN

      file_ok = .false.
      IF (p_is_master) THEN
         inquire(file=trim(soil_file), exist=file_ok)
         IF (.not. file_ok) THEN
            WRITE(*,'(2A)') 'WARNING material soil init file not found: ', trim(soil_file)
         ENDIF
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast(file_ok, 1, mpi_logical, p_address_master, p_comm_glb, p_err)
#endif
      IF (.not. file_ok) RETURN

      has_soildepth = .false.
      IF (p_is_master) THEN
         has_soildepth = ncio_var_exist(trim(soil_file), 'soildepth', readflag=.false.)
         IF (.not. has_soildepth) THEN
            WRITE(*,'(2A)') 'WARNING material soil init missing soildepth: ', trim(soil_file)
         ENDIF
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast(has_soildepth, 1, mpi_logical, p_address_master, p_comm_glb, p_err)
#endif
      IF (.not. has_soildepth) RETURN

      CALL ncio_read_bcast_serial(trim(soil_file), 'soildepth', soil_z)
      nl_soil_ini = size(soil_z)
      IF (nl_soil_ini <= 0) RETURN

      month = min(12, max(1, init_month))
      ncopy = min(nl_soil_ini, nl_soil)

      CALL gsoil%define_from_file(trim(soil_file), latname='lat', lonname='lon')

      IF (p_is_io) THEN
         CALL allocate_block_data(gsoil, delta_grid, nl_soil_ini)
      ENDIF

      IF (p_is_worker .and. numpatch > 0) THEN
         allocate(delta_patch(nl_soil_ini, numpatch))
      ENDIF

      CALL msoil2p%build_arealweighted(gsoil, landpatch)

      IF (p_is_master) THEN
         WRITE(*,'(A,I0,2A)') 'Use material soil isotope initial field for month ', &
            month, ' from ', trim(soil_file)
      ENDIF

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         CALL material_soil_init_varname(imt, varname)
         IF (.not. material_is_isotope(imt)) THEN
            IF (len_trim(varname) > 0 .and. trim(varname) /= 'null' .and. p_is_master) THEN
               WRITE(*,'(A,I0,3A)') 'WARNING material soil init ignores non-isotope material ', &
                  imt, ' variable ', trim(varname), '.'
            ENDIF
            CYCLE
         ENDIF
         IF (len_trim(varname) == 0 .or. trim(varname) == 'null') CYCLE

         has_var = .false.
         IF (p_is_master) THEN
            has_var = ncio_var_exist(trim(soil_file), trim(varname), readflag=.false.)
            IF (.not. has_var) THEN
               WRITE(*,'(4A)') 'WARNING material soil init missing ', trim(varname), &
                  ' in ', trim(soil_file)
            ENDIF
         ENDIF
#ifdef USEMPI
         CALL mpi_bcast(has_var, 1, mpi_logical, p_address_master, p_comm_glb, p_err)
#endif
         IF (.not. has_var) CYCLE

         CALL ncio_read_block_time(trim(soil_file), trim(varname), gsoil, &
            nl_soil_ini, month, delta_grid)
         CALL msoil2p%grid2pset(delta_grid, nl_soil_ini, delta_patch)

         IF (p_is_worker .and. numpatch > 0 .and. allocated(mt_wliq_soisno)) THEN
            DO ip = 1, numpatch
               DO k = 1, ncopy
                  delta_layer = delta_patch(k, ip)
                  IF (material_soil_init_valid_delta(delta_layer)) THEN
                     R_layer = delta_to_R(delta_layer, materials(imt)%ref_ratio)
                     mt_wliq_soisno(imt, k, ip) = max(wliq_soisno(k, ip), 0._r8) * R_layer
                     mt_wice_soisno(imt, k, ip) = max(wice_soisno(k, ip), 0._r8) * R_layer
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
      ENDDO

      IF (allocated(delta_patch)) deallocate(delta_patch)
      IF (allocated(soil_z)) deallocate(soil_z)
      ! delta_grid carries a FINAL procedure, but release it explicitly so the
      ! IO ranks return the block memory before this routine returns instead
      ! of waiting for the finalizer at scope exit.  The call is idempotent:
      ! the later finalizer sees blk already deallocated.
      IF (p_is_io) CALL block_data_real8_3d_free_mem(delta_grid)

   END SUBROUTINE material_soil_init_from_file

   SUBROUTINE material_soil_init_varname (imt, varname)
      IMPLICIT NONE
      integer, intent(in) :: imt
      character(len=*), intent(out) :: varname
      character(len=256) :: token

      varname = 'null'
      CALL material_soil_init_csv_token(DEF_MT_SOIL_INIT_VARS, imt, token)
      IF (len_trim(token) > 0 .and. trim(token) /= 'null') THEN
         varname = trim(token)
         RETURN
      ENDIF

      CALL ensure_isotope_physics_registered ()
      CALL isotope_default_soil_init_varname(imt, varname)
   END SUBROUTINE material_soil_init_varname

   SUBROUTINE material_soil_init_csv_token (csv, idx, token)
      IMPLICIT NONE
      character(len=*), intent(in) :: csv
      integer, intent(in) :: idx
      character(len=*), intent(out) :: token
      integer :: i, n, start, finish

      token = ''
      IF (idx <= 0) RETURN

      n = 1
      start = 1
      DO i = 1, len_trim(csv) + 1
         IF (i > len_trim(csv) .or. csv(i:i) == ',') THEN
            finish = i - 1
            IF (n == idx) THEN
               IF (finish >= start) token = adjustl(csv(start:finish))
               RETURN
            ENDIF
            n = n + 1
            start = i + 1
         ENDIF
      ENDDO
   END SUBROUTINE material_soil_init_csv_token


   logical FUNCTION material_soil_init_valid_delta (delta)
      IMPLICIT NONE
      real(r8), intent(in) :: delta

      material_soil_init_valid_delta = (delta == delta) .and. &
         abs(delta) <= mt_delta_sanity_max .and. abs(delta) < huge(1._r8)
   END FUNCTION material_soil_init_valid_delta

END MODULE MOD_MT_SoilInit
#endif
