#include <define.h>

MODULE MOD_MT_HistoryRestart
   USE MOD_Precision
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: mt_history_read_name
#ifdef MaterialTransport
   PUBLIC :: mt_history_restart_fields, mt_history_frcsat_count
   ! Keep collective validation and I/O in the host's shared restart engine.
   abstract interface
      SUBROUTINE history_field1d(file, mode, name, var, pixelset, vecname)
         USE MOD_Precision
         USE MOD_Pixelset, only: pixelset_type
         character(len=*), intent(in) :: file, name, vecname
         integer, intent(in) :: mode
         real(r8), allocatable, intent(inout) :: var(:)
         type(pixelset_type), intent(in) :: pixelset
      END SUBROUTINE history_field1d
      SUBROUTINE history_field2d(file, mode, name, var, pixelset, vecname)
         USE MOD_Precision
         USE MOD_Pixelset, only: pixelset_type
         character(len=*), intent(in) :: file, name, vecname
         integer, intent(in) :: mode
         real(r8), allocatable, intent(inout) :: var(:,:)
         type(pixelset_type), intent(in) :: pixelset
      END SUBROUTINE history_field2d
      SUBROUTINE history_field3d(file, mode, name, var, pixelset, vecname)
         USE MOD_Precision
         USE MOD_Pixelset, only: pixelset_type
         character(len=*), intent(in) :: file, name, vecname
         integer, intent(in) :: mode
         real(r8), allocatable, intent(inout) :: var(:,:,:)
         type(pixelset_type), intent(in) :: pixelset
      END SUBROUTINE history_field3d
      SUBROUTINE history_write_marker(file, name, val)
         USE MOD_Precision
         character(len=*), intent(in) :: file, name
         real(r8), intent(in) :: val
      END SUBROUTINE history_write_marker
      SUBROUTINE history_read_marker(file, name, vmin, vmax)
         USE MOD_Precision
         character(len=*), intent(in) :: file, name
         real(r8), intent(out) :: vmin, vmax
      END SUBROUTINE history_read_marker
      SUBROUTINE history_require_var(file, name, pixelset, vecname, rank, n1, n2)
         USE MOD_Pixelset, only: pixelset_type
         character(len=*), intent(in) :: file, name, vecname
         type(pixelset_type), intent(in) :: pixelset
         integer, intent(in) :: rank, n1, n2
      END SUBROUTINE history_require_var
      SUBROUTINE history_global_max(n)
         integer, intent(inout) :: n(:)
      END SUBROUTINE history_global_max
   end interface
#endif

CONTAINS

   FUNCTION mt_history_read_name (file, name, pixelset) RESULT (read_name)

   USE MOD_Pixelset, only: pixelset_type
   USE MOD_NetCDFVector, only: ncio_vector_exist
   IMPLICIT NONE

   character(len=*), intent(in) :: file, name
   type(pixelset_type), intent(in) :: pixelset
   character(len=:), allocatable :: read_name

      read_name = name
      IF (index(name, 'a_mt_') /= 1) RETURN
      ! Prefer the new field globally, even when only some blocks contain it.
      ! require_var must reject an incomplete/corrupt new field, not read stale data.
      IF (.not. ncio_vector_exist(file, name, pixelset)) read_name = 'a_trc_'//name(6:)

   END FUNCTION mt_history_read_name

#ifdef MaterialTransport
   SUBROUTINE mt_history_restart_fields(file, mode, write_mode, field1d, field2d, field3d)
      USE MOD_MT_Vars
      USE MOD_LandPatch, only: landpatch
      USE MOD_SPMD_Task, only: p_is_master
      character(len=*), intent(in) :: file
      integer, intent(in) :: mode, write_mode
      procedure(history_field1d) :: field1d
      procedure(history_field2d) :: field2d
      procedure(history_field3d) :: field3d

      ! Material sums (MOD_MT_Vars, allocate_Material_Vars). After a cold
      ! start of the materials the saved sums belong to another state (or
      ! another species list) and are not read; their window starts again.
      IF (mode == write_mode .or. mt_land_restart_loaded) THEN
         CALL field2d (file, mode, 'a_mt_precip', a_mt_precip, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_vapor_exchange', a_mt_vapor_exchange, landpatch, 'patch')
         CALL field2d (file, mode, 'a_water_precip', a_water_precip, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_evap', a_mt_evap, landpatch, 'patch')
         CALL field2d (file, mode, 'a_water_evap_gross', a_water_evap_gross, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_transp', a_mt_transp, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_transp_src', a_mt_transp_src, landpatch, 'patch')
         CALL field2d (file, mode, 'a_water_transp', a_water_transp, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_soilevap', a_mt_soilevap, landpatch, 'patch')
         CALL field2d (file, mode, 'a_water_soilevap', a_water_soilevap, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_canopyevap', a_mt_canopyevap, landpatch, 'patch')
         CALL field2d (file, mode, 'a_water_canopyevap', a_water_canopyevap, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_subl', a_mt_subl, landpatch, 'patch')
         CALL field2d (file, mode, 'a_water_subl', a_water_subl, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_wetland_evap', a_mt_wetland_evap, landpatch, 'patch')
         CALL field2d (file, mode, 'a_water_wetland_evap', a_water_wetland_evap, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_rsur', a_mt_rsur, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_rsub', a_mt_rsub, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_rnof', a_mt_rnof, landpatch, 'patch')
         CALL field2d (file, mode, 'a_water_rnof', a_water_rnof, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_qinfl', a_mt_qinfl, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_qcharge', a_mt_qcharge, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_ldew_mass', a_mt_ldew_mass, landpatch, 'patch')
         CALL field1d (file, mode, 'a_water_ldew', a_water_ldew, landpatch, 'patch')
         CALL field3d (file, mode, 'a_mt_soil_mass', a_mt_soil_mass, landpatch, 'patch')
         CALL field2d (file, mode, 'a_water_soil', a_water_soil, landpatch, 'patch')
         CALL field3d (file, mode, 'a_mt_snow_mass', a_mt_snow_mass, landpatch, 'patch')
         CALL field2d (file, mode, 'a_water_snow', a_water_snow, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_wa_mass', a_mt_wa_mass, landpatch, 'patch')
         CALL field1d (file, mode, 'a_water_wa', a_water_wa, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_wa_debt_mass', a_mt_wa_debt_mass, landpatch, 'patch')
         CALL field1d (file, mode, 'a_water_wa_debt', a_water_wa_debt, landpatch, 'patch')
         CALL field1d (file, mode, 'a_water_aquifer_actual', a_water_aquifer_actual, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_aquifer_actual_mass', a_mt_aquifer_actual_mass, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_wdsrf_mass', a_mt_wdsrf_mass, landpatch, 'patch')
         CALL field1d (file, mode, 'a_water_wdsrf', a_water_wdsrf, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_wetwat_mass', a_mt_wetwat_mass, landpatch, 'patch')
         CALL field1d (file, mode, 'a_water_wetwat', a_water_wetwat, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_surface_residue_mass', a_mt_surface_residue_mass, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_subsurface_residue_mass', a_mt_subsurface_residue_mass, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_layer_dry_mass', a_mt_layer_dry_mass, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_solid_mass', a_mt_solid_mass, landpatch, 'patch')
         CALL field2d (file, mode, 'a_mt_scv_mass', a_mt_scv_mass, landpatch, 'patch')
         CALL field1d (file, mode, 'a_water_scv', a_water_scv, landpatch, 'patch')
      ELSEIF (p_is_master) THEN
         write(*,*) 'Note: materials were cold-started; their history window starts at the model start.'
      ENDIF
   END SUBROUTINE mt_history_restart_fields

   SUBROUTINE mt_history_frcsat_count(file, mode, write_mode, nac, a_frcsat, &
      hacc_field1d, write_marker, read_marker, require_var, global_max)

   USE, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   USE netcdf
   USE MOD_Block, only: get_filename_block
   USE MOD_LandPatch, only: landpatch
   USE MOD_MT_Namelist, only: mt_any_enabled
   USE MOD_Vars_Global, only: spval
   USE MOD_MT_HostHistory, only: nac_frcsat
   USE MOD_SPMD_Task, only: p_is_io, p_is_master, CoLM_stop
   IMPLICIT NONE

   character(len=*), intent(in) :: file
   integer, intent(in) :: mode, write_mode
   real(r8), intent(in) :: nac
   real(r8), allocatable, intent(in) :: a_frcsat(:)
   procedure(history_field1d) :: hacc_field1d
   procedure(history_write_marker) :: write_marker
   procedure(history_read_marker) :: read_marker
   procedure(history_require_var) :: require_var
   procedure(history_global_max) :: global_max
   character(len=256) :: fileblock
   integer :: iblkgrp, ncid, varid, status, found(2), bad(1)
   real(r8) :: vmin, vmax

      IF (.not. mt_any_enabled()) RETURN

      IF (mode == write_mode) THEN
         CALL write_marker (file, 'hacc_frcsat_count_version', 1._r8)
         CALL hacc_field1d (file, mode, 'nac_frcsat', nac_frcsat, landpatch, 'patch')
         RETURN
      ENDIF

      found = 0
      bad = 0
      IF (p_is_io) THEN
         DO iblkgrp = 1, landpatch%nblkgrp
            CALL get_filename_block (file, landpatch%xblkgrp(iblkgrp), landpatch%yblkgrp(iblkgrp), fileblock)
            IF (nf90_open (trim(fileblock), NF90_NOWRITE, ncid) /= NF90_NOERR) THEN
               bad = 1
               CYCLE
            ENDIF
            status = nf90_inq_varid (ncid, 'hacc_frcsat_count_version', varid)
            IF (status == NF90_NOERR) THEN
               found(1) = 1
            ELSEIF (status /= NF90_ENOTVAR) THEN
               bad = 1
            ENDIF
            status = nf90_inq_varid (ncid, 'nac_frcsat', varid)
            IF (status == NF90_NOERR) THEN
               found(2) = 1
            ELSEIF (status /= NF90_ENOTVAR) THEN
               bad = 1
            ENDIF
            IF (nf90_close (ncid) /= NF90_NOERR) bad = 1
         ENDDO
      ENDIF
      CALL global_max (bad)
      IF (bad(1) /= 0) CALL CoLM_stop ('Cannot inspect frcsat history sample-count metadata')
      CALL global_max (found)

      IF (found(1) == 0) THEN
         IF (found(2) /= 0) CALL CoLM_stop ('frcsat history count has no version marker')
         ! Keep the legacy sum, but never guess how many dry steps it contains.
         IF (allocated(nac_frcsat)) nac_frcsat = -1._r8
         IF (p_is_master) write(*,'(A)') &
            'Warning: legacy history restart has no frcsat valid-sample count; this window is missing for f_frcsat.'
         RETURN
      ENDIF

      CALL require_var (file, 'hacc_frcsat_count_version', landpatch, 'patch', 1, 0, 0)
      CALL read_marker (file, 'hacc_frcsat_count_version', vmin, vmax)
      IF (vmin /= 1._r8 .or. vmax /= 1._r8) &
         CALL CoLM_stop ('Unsupported frcsat history sample-count version')
      CALL hacc_field1d (file, mode, 'nac_frcsat', nac_frcsat, landpatch, 'patch')

      bad = 0
      IF (allocated(nac_frcsat)) THEN
         IF (any(.not. ieee_is_finite(nac_frcsat))) THEN
            bad = 1
         ELSEIF (any(nac_frcsat < -1._r8) .or. any(nac_frcsat > nac) .or. &
                 any(nac_frcsat /= anint(nac_frcsat))) THEN
            bad = 1
         ELSEIF (any((nac_frcsat >= 0._r8) .and. &
                     ((nac_frcsat == 0._r8) .neqv. (a_frcsat == spval)))) THEN
            bad = 1
         ENDIF
      ENDIF
      CALL global_max (bad)
      IF (bad(1) /= 0) CALL CoLM_stop ('Invalid frcsat history valid-sample count')

   END SUBROUTINE mt_history_frcsat_count
#endif

END MODULE MOD_MT_HistoryRestart
