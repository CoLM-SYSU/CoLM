#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_Forcing

   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   USE netcdf, only: NF90_NOERR, NF90_NOWRITE, NF90_MAX_NAME, NF90_MAX_VAR_DIMS, &
      nf90_open, nf90_close, nf90_inq_varid, nf90_inquire_variable, nf90_inquire_dimension
   USE MOD_Precision
   USE MOD_Namelist, only: DEF_forcing, DEF_Forcing_Interp_Method, DEF_REST_CompressLevel, &
      DEF_dir_forcing, DEF_USE_PFT, DEF_SOLO_PFT, DEF_FAST_PC
   USE MOD_Const_LC, only: patchtypes
   USE MOD_Vars_Global, only: CROPLAND, N_land_classification
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_NetCDFBlock, only: ncio_read_block_time
   USE MOD_NetCDFSerial, only: ncio_var_exist, ncio_inquire_varsize, ncio_read_serial, &
      ncio_write_serial, ncio_define_dimension
   USE MOD_NetCDFVector, only: ncio_read_vector, ncio_write_vector
   USE MOD_Block, only: get_filename_block
   USE MOD_SpatialMapping
   USE MOD_LandPatch
   USE MOD_TimeManager
   USE MOD_UserSpecifiedForcing, only: metfilename
   USE MOD_MT_Namelist, only: DEF_MT_FORCING_MIN_PRECIP, DEF_MT_FORCING_MEMORY_PRECIP
   USE MOD_MT_Defs, only: nmts, materials, material_precip_default_ratio, &
      material_vapor_default_ratio, &
      material_is_isotope, delta_to_R, mt_tiny, mt_delta_sanity_max, &
      mt_uses_generic_water, material_lower
   USE MOD_MT_Frac, only: material_fractionation_active, material_alpha_liq_vap, &
      material_alpha_ice_vap_deposition
   USE MOD_Vars_1DForcing, only: forc_t
   USE MOD_MT_Vars, only: mt_runtime_forced
   USE MOD_MT_Providers, only: ensure_isotope_physics_registered
   USE MOD_MT_ForcingInput, only: material_forcing_spec_type, material_forcing_input_load, &
      material_forcing_input_get, material_forcing_input_find, material_forcing_input_final

   IMPLICIT NONE
   SAVE

   integer, parameter :: STREAM_PRECIP = 1
   integer, parameter :: STREAM_VAPOR  = 2
   integer, parameter :: STREAM_TOTAL_PRECIP = 3
   integer, parameter :: STREAM_TOTAL_VAPOR  = 4

   integer, parameter :: DECODE_INVALID  = 0
   integer, parameter :: DECODE_NEAR_DRY = 1
   integer, parameter :: DECODE_VALID    = 2

   integer, parameter :: MODE_VALUE = 1
   integer, parameter :: MODE_DELTA = 2
   integer, parameter :: MODE_HEAVY_OVER_TOTAL = 3
   integer, parameter :: MODE_NORMALIZED_OVER_TOTAL = 4

   ! Heavy/total precipitation ratios are numerically unreliable for
   ! near-dry records; the threshold is DEF_MT_FORCING_MIN_PRECIP (see
   ! MOD_MT_Namelist for the IsoGSM evidence behind its default).
   real(r8), parameter :: mt_forc_min_q    = 1.0e-12_r8
   real(r8), parameter :: mt_forc_max_abs  = 1.0e10_r8
   ! Schema 2 adds mt_forcing_precip_reliable.  A schema-1 restart is still
   ! read: its remembered ratio was the last valid one, so it counts as reliable.
   integer, parameter :: MT_FORC_CACHE_SCHEMA = 2
   integer, parameter :: MT_FORC_ID_WIDTH = 8 + 6*256

   logical :: mt_runtime_forcing_enabled = .false.
   integer :: mt_forcing_log_count = 0
   integer :: n_mt_forc_vars = 0
   integer :: idx_total_precip = 0
   integer :: idx_total_vapor = 0

   integer, allocatable :: mt_var_stream(:)
   integer, allocatable :: mt_var_imt(:)
   integer, allocatable :: mt_var_mode(:)
   integer, allocatable :: mt_var_total(:)
   integer, allocatable :: mt_var_dtime(:)
   integer, allocatable :: mt_var_offset(:)
   character(len=256), allocatable :: mt_var_fprefix(:)
   character(len=256), allocatable :: mt_var_vname(:)
   character(len=256), allocatable :: mt_var_tintalgo(:)
   character(len=256), allocatable :: mt_var_timelog(:)
   type(timestamp), allocatable :: mt_tstamp_LB(:)
   type(timestamp), allocatable :: mt_tstamp_UB(:)

   type(spatial_mapping_type) :: mt_mg2p_forc
   type(grid_type) :: mt_gforc
   type(block_data_real8_2d) :: mt_metdata
   type(block_data_real8_2d), allocatable :: mt_forcn(:)
   type(block_data_real8_2d), allocatable :: mt_forcn_LB(:)
   type(block_data_real8_2d), allocatable :: mt_forcn_UB(:)
   real(r8), allocatable :: mt_forc_patch(:,:)

   real(r8), allocatable :: mt_forc_precip_value(:,:)
   real(r8), allocatable :: mt_forc_vapor_value(:,:)
   ! Precipitation ratio of the current step, valid where mt_forc_has_precip.
   ! mt_forc_precip_value is the remembered ratio used when the current
   ! record is near-dry (see material_forcing_snow_value).
   real(r8), allocatable :: mt_forc_precip_step(:,:)
   ! 1 where mt_forc_precip_value holds a ratio of a reliable record
   ! (total precipitation above DEF_MT_FORCING_MEMORY_PRECIP), 0 where it
   ! still holds the vapour-equilibrium estimate.  Kept as real so that it is
   ! restarted and LULCC-remapped like the ratios.
   real(r8), allocatable :: mt_forc_precip_reliable(:,:)
   ! Species whose below-threshold precipitation ratio is filled by phase
   ! (rain from the vapour equilibrium, snow from the remembered ratio):
   ! isotopes with both a precipitation and a vapour forcing.
   logical,  allocatable :: mt_forc_phase_fill(:)
   ! Last-valid compositions survive the annual patch rebuild; they are
   ! intensive forcing state, not material mass.
   real(r8), allocatable :: lulcc_precip_old(:,:), lulcc_vapor_old(:,:)
   real(r8), allocatable :: lulcc_precip_new(:,:), lulcc_vapor_new(:,:)
   real(r8), allocatable :: lulcc_reliable_old(:,:), lulcc_reliable_new(:,:)
   real(r8), allocatable :: lulcc_old_patch_area(:)
   logical, allocatable :: lulcc_phase_fill(:)
   logical,  allocatable :: mt_forc_has_precip(:,:)
   logical,  allocatable :: mt_forc_has_vapor(:,:)
   integer,  allocatable :: mt_forc_precip_status(:,:)
   integer,  allocatable :: mt_forc_vapor_status(:,:)
   logical,  allocatable :: mt_forc_precip_warned_dry(:)
   logical,  allocatable :: mt_forc_precip_warned_invalid(:)
   logical,  allocatable :: mt_forc_vapor_warned_dry(:)
   logical,  allocatable :: mt_forc_vapor_warned_invalid(:)

   PUBLIC :: mt_forcing_init
   PUBLIC :: mt_read_forcing
   PUBLIC :: mt_forcing_reset
   PUBLIC :: mt_forcing_final
   PUBLIC :: material_forcing_precip_value
   PUBLIC :: material_forcing_rain_value, material_forcing_snow_value
   PUBLIC :: material_forcing_vapor_value
   PUBLIC :: material_forcing_precip_ratio
   PUBLIC :: material_forcing_vapor_ratio
   PUBLIC :: material_forcing_has_precip
   PUBLIC :: material_forcing_has_vapor
   PUBLIC :: material_forcing_ratio_to_delta
   PUBLIC :: material_forcing_write_restart, mt_forcing_read_restart
   PUBLIC :: mt_forcing_lulcc_save, mt_forcing_lulcc_remap, mt_forcing_lulcc_restore

CONTAINS

   SUBROUTINE mt_forcing_init (gforc, numpatch)
      IMPLICIT NONE
      type(grid_type), intent(in) :: gforc
      integer, intent(in) :: numpatch
      integer :: iv

      CALL mt_forcing_final()
      IF (nmts <= 0) RETURN

      CALL material_forcing_configure()
      CALL material_forcing_allocate_state(numpatch)
      IF (.not. mt_runtime_forcing_enabled) RETURN
      IF (trim(DEF_forcing%dataset) == 'POINT') THEN
         IF (p_is_master) THEN
            WRITE(*,'(A)') 'ERROR material runtime forcing does not support DEF_forcing%dataset=POINT.'
            WRITE(*,'(A)') '      Disable material forcing variables or use gridded forcing for material fields.'
         ENDIF
         CALL CoLM_stop()
      ENDIF

      CALL mt_gforc%define_by_copy(gforc)
      IF (trim(DEF_Forcing_Interp_Method) == 'arealweight') THEN
         CALL mt_mg2p_forc%build_arealweighted(mt_gforc, landpatch)
      ELSEIF (trim(DEF_Forcing_Interp_Method) == 'bilinear') THEN
         CALL mt_mg2p_forc%build_bilinear(mt_gforc, landpatch)
      ELSE
         IF (p_is_master) WRITE(*,'(2A)') 'ERROR material forcing unsupported interp method: ', &
            trim(DEF_Forcing_Interp_Method)
         CALL CoLM_stop()
      ENDIF

      allocate(mt_forcn(n_mt_forc_vars))
      IF (p_is_io) THEN
         allocate(mt_forcn_LB(n_mt_forc_vars))
         allocate(mt_forcn_UB(n_mt_forc_vars))
         DO iv = 1, n_mt_forc_vars
            CALL allocate_block_data(mt_gforc, mt_forcn(iv))
            CALL allocate_block_data(mt_gforc, mt_forcn_LB(iv))
            CALL allocate_block_data(mt_gforc, mt_forcn_UB(iv))
         ENDDO
         CALL allocate_block_data(mt_gforc, mt_metdata)
      ENDIF

      allocate(mt_forc_patch(numpatch, n_mt_forc_vars))
      mt_forc_patch(:,:) = 0._r8

      IF (p_is_master) THEN
         WRITE(*,'(A,I0,A)') 'Material runtime forcing enabled: ', n_mt_forc_vars, &
            ' variable(s) from main forcing namelist.'
      ENDIF
   END SUBROUTINE mt_forcing_init

   SUBROUTINE mt_read_forcing (idate, dir_forcing)
      IMPLICIT NONE
      integer, intent(in) :: idate(3)
      character(len=*), intent(in) :: dir_forcing
      integer :: iv

      CALL material_forcing_prepare_step()
      IF (.not. mt_runtime_forcing_enabled) RETURN

      CALL material_forcing_read_LBUB(idate, dir_forcing)

      IF (p_is_io) THEN
         DO iv = 1, n_mt_forc_vars
            CALL material_forcing_interpolate_var(idate, iv)
         ENDDO
      ENDIF

      DO iv = 1, n_mt_forc_vars
         CALL mt_mg2p_forc%grid2pset(mt_forcn(iv), mt_forc_patch(:, iv))
      ENDDO

      CALL material_forcing_update_values()
      CALL material_forcing_log_ranges(idate)
   END SUBROUTINE mt_read_forcing

   SUBROUTINE mt_forcing_reset ()
      IMPLICIT NONE

      IF (allocated(mt_tstamp_LB)) mt_tstamp_LB(:) = timestamp(-1, -1, -1)
      IF (allocated(mt_tstamp_UB)) mt_tstamp_UB(:) = timestamp(-1, -1, -1)
      mt_forcing_log_count = 0
   END SUBROUTINE mt_forcing_reset

   SUBROUTINE mt_forcing_final ()
      IMPLICIT NONE

      CALL mt_mg2p_forc%forc_free_mem()
      IF (allocated(mt_forcn)) deallocate(mt_forcn)
      IF (allocated(mt_forcn_LB)) deallocate(mt_forcn_LB)
      IF (allocated(mt_forcn_UB)) deallocate(mt_forcn_UB)
      IF (allocated(mt_forc_patch)) deallocate(mt_forc_patch)
      CALL material_forcing_deallocate_config()
      CALL material_forcing_deallocate_state()
      CALL material_forcing_input_final()

      IF (allocated(mt_runtime_forced)) mt_runtime_forced(:) = .false.
      mt_runtime_forcing_enabled = .false.
      mt_forcing_log_count = 0
      n_mt_forc_vars = 0
      idx_total_precip = 0
      idx_total_vapor = 0
   END SUBROUTINE mt_forcing_final

   SUBROUTINE mt_forcing_lulcc_save ()
      IF (.not. p_is_worker .or. .not. mt_runtime_forcing_enabled) RETURN
      IF (.not. allocated(mt_forc_precip_value) .or. allocated(lulcc_precip_old)) &
         CALL CoLM_stop('material forcing LULCC snapshot is unavailable or already active')
      IF (size(mt_forc_precip_value,2) == 0) RETURN
      IF (.not. allocated(elm_patch%subfrc)) &
         CALL CoLM_stop('material forcing LULCC requires old element-patch area fractions')
      IF (size(elm_patch%subfrc) /= size(mt_forc_precip_value,2)) &
         CALL CoLM_stop('material forcing LULCC old element-patch area shape mismatch')
      lulcc_old_patch_area = elm_patch%subfrc
      CALL move_alloc(mt_forc_precip_value, lulcc_precip_old)
      CALL move_alloc(mt_forc_vapor_value, lulcc_vapor_old)
      CALL move_alloc(mt_forc_precip_reliable, lulcc_reliable_old)
      CALL move_alloc(mt_forc_phase_fill, lulcc_phase_fill)
   END SUBROUTINE mt_forcing_lulcc_save

   ! Called while LulccDriver still holds both patch maps and the transfer
   ! trace. Forcing init occurs later, so keep the mapped values until then.
   SUBROUTINE mt_forcing_lulcc_remap (patchclass_new, eindex_new, patchclass_old, eindex_old, &
      lccpct_patches)
      integer, intent(in) :: patchclass_new(:), patchclass_old(:)
      integer*8, intent(in) :: eindex_new(:), eindex_old(:)
      real(r8), intent(in), optional :: lccpct_patches(:,0:)
      integer :: imt, unsupported, ignored
      integer :: raw, source_class(0:N_land_classification)

      IF (.not. p_is_worker .or. .not. allocated(lulcc_precip_old)) RETURN
      IF (allocated(lulcc_precip_new)) CALL CoLM_stop('material forcing LULCC remap repeated')
      allocate(lulcc_precip_new(nmts, size(patchclass_new)))
      allocate(lulcc_vapor_new(nmts, size(patchclass_new)))
      ! A patch without a source starts unreliable: its remembered ratio is
      ! replaced by the vapour-equilibrium estimate at the next forcing update.
      allocate(lulcc_reliable_new(nmts, size(patchclass_new)))
      lulcc_reliable_new(:,:) = 0._r8
      DO imt = 1, nmts
         lulcc_precip_new(imt,:) = material_precip_default_ratio(imt)
         lulcc_vapor_new(imt,:) = material_vapor_default_ratio(imt)
      ENDDO
      source_class = [(raw, raw=0,N_land_classification)]
      IF ((DEF_USE_PFT .and. .not. DEF_SOLO_PFT) .or. DEF_FAST_PC) THEN
         DO raw = 1, N_land_classification
            IF (patchtypes(raw) /= 0) CYCLE
            source_class(raw) = 1
            IF (DEF_FAST_PC .and. (raw == CROPLAND .or. raw == 14)) &
               source_class(raw) = CROPLAND
         ENDDO
      ENDIF
      ! Only reliable records contribute to the inherited snow memory.
      DO imt = 1, nmts
         IF (lulcc_phase_fill(imt)) &
            lulcc_precip_old(imt,:) = lulcc_precip_old(imt,:) * lulcc_reliable_old(imt,:)
      ENDDO
      CALL material_forcing_lulcc_map(lulcc_precip_old, lulcc_precip_new, &
         patchclass_new, eindex_new, patchclass_old, eindex_old, unsupported, &
         lccpct_patches, lulcc_old_patch_area, source_class)
      CALL material_forcing_lulcc_map(lulcc_vapor_old, lulcc_vapor_new, &
         patchclass_new, eindex_new, patchclass_old, eindex_old, ignored, &
         lccpct_patches, lulcc_old_patch_area, source_class)
      ! Normalize by reliable source area, not by the area of all sources.
      CALL material_forcing_lulcc_map(lulcc_reliable_old, lulcc_reliable_new, &
         patchclass_new, eindex_new, patchclass_old, eindex_old, ignored, &
         lccpct_patches, lulcc_old_patch_area, source_class)
      DO imt = 1, nmts
         IF (.not. lulcc_phase_fill(imt)) CYCLE
         WHERE (lulcc_reliable_new(imt,:) > 0._r8)
            lulcc_precip_new(imt,:) = lulcc_precip_new(imt,:) / lulcc_reliable_new(imt,:)
            lulcc_reliable_new(imt,:) = 1._r8
         ELSEWHERE
            lulcc_precip_new(imt,:) = material_precip_default_ratio(imt)
         END WHERE
      ENDDO
      lulcc_reliable_new = merge(1._r8, 0._r8, lulcc_reliable_new >= 0.5_r8)
      IF (unsupported > 0) WRITE(*,'(A,I0,A)') 'Material forcing LULCC: ', unsupported, &
         ' patch(es) have no eligible transfer source; using configured forcing defaults.'
      deallocate(lulcc_precip_old, lulcc_vapor_old, lulcc_reliable_old, lulcc_old_patch_area, lulcc_phase_fill)
   END SUBROUTINE mt_forcing_lulcc_remap

   SUBROUTINE mt_forcing_lulcc_restore ()
      IF (.not. p_is_worker) RETURN
      IF (allocated(lulcc_precip_old)) THEN
         IF (.not. allocated(mt_forc_precip_value)) &
            CALL CoLM_stop('material forcing LULCC restore precedes forcing init')
         IF (size(mt_forc_precip_value,2) > 0) &
            CALL CoLM_stop('material forcing LULCC remap missing for rebuilt patches')
         deallocate(lulcc_precip_old, lulcc_vapor_old, lulcc_reliable_old, lulcc_old_patch_area, lulcc_phase_fill)
      ENDIF
      IF (.not. allocated(lulcc_precip_new)) RETURN
      IF (.not. allocated(mt_forc_precip_value)) &
         CALL CoLM_stop('material forcing LULCC restore precedes forcing init')
      IF (any(shape(mt_forc_precip_value) /= shape(lulcc_precip_new)) .or. &
          any(shape(mt_forc_vapor_value) /= shape(lulcc_vapor_new)) .or. &
          any(shape(mt_forc_precip_reliable) /= shape(lulcc_reliable_new))) &
         CALL CoLM_stop('material forcing LULCC patch/species count changed after remap')
      CALL move_alloc(lulcc_precip_new, mt_forc_precip_value)
      CALL move_alloc(lulcc_vapor_new, mt_forc_vapor_value)
      CALL move_alloc(lulcc_reliable_new, mt_forc_precip_reliable)
   END SUBROUTINE mt_forcing_lulcc_restore

   SUBROUTINE material_forcing_lulcc_map (old, mapped, patchclass_new, eindex_new, &
      patchclass_old, eindex_old, unsupported, lccpct_patches, old_patch_area, source_class)
      real(r8), intent(in) :: old(:,:)
      real(r8), intent(inout) :: mapped(:,:)
      integer, intent(in) :: patchclass_new(:), patchclass_old(:)
      integer*8, intent(in) :: eindex_new(:), eindex_old(:)
      integer, intent(out) :: unsupported
      real(r8), intent(in), optional :: lccpct_patches(:,0:), old_patch_area(:)
      integer, intent(in), optional :: source_class(0:)
      integer :: np, op, oq, c, raw, source_count, source_patch
      real(r8) :: weight, total_weight
      real(r8), allocatable :: class_area(:)
      logical :: same_class, old_sorted
      integer :: lo, hi

      IF (size(old,1) /= size(mapped,1) .or. size(old,2) /= size(patchclass_old) .or. &
          size(old,2) /= size(eindex_old) .or. size(mapped,2) /= size(patchclass_new) .or. &
          size(mapped,2) /= size(eindex_new)) CALL CoLM_stop('material forcing LULCC map shape mismatch')
      IF (present(lccpct_patches)) THEN
         IF (size(lccpct_patches,1) < size(patchclass_new)) &
            CALL CoLM_stop('material forcing LULCC trace is shorter than new patch map')
         IF (any(.not. ieee_is_finite(lccpct_patches))) &
            CALL CoLM_stop('material forcing LULCC trace has invalid source fractions')
         IF (any(lccpct_patches < 0._r8)) &
            CALL CoLM_stop('material forcing LULCC trace has invalid source fractions')
      ENDIF
      IF (present(old_patch_area)) THEN
         IF (size(old_patch_area) /= size(patchclass_old)) &
            CALL CoLM_stop('material forcing LULCC old-area shape mismatch')
         IF (any(.not. ieee_is_finite(old_patch_area))) &
            CALL CoLM_stop('material forcing LULCC old-area fractions are invalid')
         IF (any(old_patch_area < 0._r8)) &
            CALL CoLM_stop('material forcing LULCC old-area fractions are invalid')
      ENDIF
      IF (present(source_class) .and. present(lccpct_patches)) THEN
         IF (ubound(source_class,1) /= ubound(lccpct_patches,2)) &
            CALL CoLM_stop('material forcing LULCC source-class map shape mismatch')
      ENDIF
      IF (any(.not. ieee_is_finite(old))) &
         CALL CoLM_stop('material forcing LULCC old cache contains non-finite values')
      old_sorted = .true.
      IF (size(eindex_old) > 1) old_sorted = all(eindex_old(2:) >= eindex_old(:size(eindex_old)-1))
      IF (present(lccpct_patches) .and. present(old_patch_area)) THEN
         allocate(class_area(size(patchclass_old)))
         class_area = 0._r8
         DO op = 1, size(patchclass_old)
            CALL element_range (eindex_old(op), lo, hi)
            DO oq = lo, hi
               IF (eindex_old(oq) == eindex_old(op) .and. &
                   patchclass_old(oq) == patchclass_old(op)) &
                  class_area(op) = class_area(op) + old_patch_area(oq)
            ENDDO
         ENDDO
      ENDIF
      unsupported = 0
      DO np = 1, size(patchclass_new)
         total_weight = 0._r8
         source_count = 0
         source_patch = 0
         same_class = .false.
         ! SAT has no transfer trace: prefer the unchanged class, then the
         ! same element's area mean for a newly created class.
         CALL element_range (eindex_new(np), lo, hi)
         IF (.not. present(lccpct_patches)) THEN
            DO op = lo, hi
               IF (eindex_old(op) == eindex_new(np) .and. &
                   patchclass_old(op) == patchclass_new(np)) THEN
                  IF (present(old_patch_area)) THEN
                     IF (old_patch_area(op) <= tiny(1._r8)) CYCLE
                  ENDIF
                  same_class = .true.
               ENDIF
            ENDDO
         ENDIF
         DO op = lo, hi
            IF (eindex_old(op) /= eindex_new(np)) CYCLE
            IF (present(lccpct_patches)) THEN
               c = patchclass_old(op)
               weight = 0._r8
               IF (present(source_class)) THEN
                  DO raw = 0, ubound(source_class,1)
                     IF (source_class(raw) == c) &
                        weight = weight + lccpct_patches(np,raw)
                  ENDDO
               ELSE
                  IF (c < 0 .or. c > ubound(lccpct_patches,2)) CYCLE
                  weight = lccpct_patches(np,c)
               ENDIF
               IF (weight <= 0._r8) CYCLE
               IF (present(old_patch_area)) THEN
                  IF (class_area(op) <= tiny(1._r8)) CYCLE
                  weight = weight * old_patch_area(op) / class_area(op)
               ENDIF
            ELSE
               IF (same_class .and. patchclass_old(op) /= patchclass_new(np)) CYCLE
               weight = 1._r8
               IF (present(old_patch_area)) weight = max(0._r8, old_patch_area(op))
            ENDIF
            IF (weight <= 0._r8) CYCLE
            IF (total_weight == 0._r8) mapped(:,np) = 0._r8
            mapped(:,np) = mapped(:,np) + weight * old(:,op)
            total_weight = total_weight + weight
            source_count = source_count + 1
            source_patch = op
         ENDDO
         IF (total_weight > 0._r8) THEN
            IF (source_count == 1) THEN
               mapped(:,np) = old(:,source_patch)
            ELSE
               mapped(:,np) = mapped(:,np) / total_weight
            ENDIF
         ELSE
            unsupported = unsupported + 1
         ENDIF
      ENDDO

   CONTAINS

      SUBROUTINE element_range (eidx, lo, hi)
         integer*8, intent(in) :: eidx
         integer, intent(out) :: lo, hi
         integer :: a, b, m

         IF (.not. old_sorted) THEN
            lo = 1
            hi = size(eindex_old)
            RETURN
         ENDIF
         a = 1
         b = size(eindex_old) + 1
         DO WHILE (a < b)
            m = (a + b) / 2
            IF (eindex_old(m) < eidx) THEN
               a = m + 1
            ELSE
               b = m
            ENDIF
         ENDDO
         lo = a
         b = size(eindex_old) + 1
         DO WHILE (a < b)
            m = (a + b) / 2
            IF (eindex_old(m) <= eidx) THEN
               a = m + 1
            ELSE
               b = m
            ENDIF
         ENDDO
         hi = a - 1
      END SUBROUTINE element_range

   END SUBROUTINE material_forcing_lulcc_map

   SUBROUTINE forcing_identity_put (identity, column, slot, value)
      integer, intent(inout) :: identity(:,:)
      integer, intent(in) :: column, slot
      character(len=*), intent(in) :: value
      integer :: k, start

      start = 8 + (slot-1)*256
      DO k = 1, min(len_trim(value), 256)
         identity(start+k, column) = iachar(value(k:k))
      ENDDO
   END SUBROUTINE forcing_identity_put

   SUBROUTINE material_forcing_identity (identity)
      integer, allocatable, intent(out) :: identity(:,:)
      integer :: iv

      allocate(identity(MT_FORC_ID_WIDTH, n_mt_forc_vars+1))
      identity = 0
      identity(1:8,1) = [n_mt_forc_vars, nmts, DEF_forcing%startyr, &
         DEF_forcing%startmo, merge(1,0,DEF_forcing%leapyear), idx_total_precip, &
         idx_total_vapor, 0]
      CALL forcing_identity_put(identity, 1, 1, DEF_forcing%dataset)
      CALL forcing_identity_put(identity, 1, 2, DEF_forcing%groupby)
      CALL forcing_identity_put(identity, 1, 3, DEF_Forcing_Interp_Method)
      CALL forcing_identity_put(identity, 1, 4, DEF_forcing%fprefix(2))
      CALL forcing_identity_put(identity, 1, 5, DEF_forcing%fprefix(4))
      CALL forcing_identity_put(identity, 1, 6, DEF_dir_forcing)
      DO iv = 1, n_mt_forc_vars
         identity(1:8,iv+1) = [mt_var_stream(iv), mt_var_imt(iv), mt_var_mode(iv), &
            mt_var_total(iv), mt_var_dtime(iv), mt_var_offset(iv), 0, 0]
         CALL forcing_identity_put(identity, iv+1, 1, mt_var_fprefix(iv))
         CALL forcing_identity_put(identity, iv+1, 2, mt_var_vname(iv))
         CALL forcing_identity_put(identity, iv+1, 3, mt_var_tintalgo(iv))
         CALL forcing_identity_put(identity, iv+1, 4, mt_var_timelog(iv))
      ENDDO
   END SUBROUTINE material_forcing_identity

   logical FUNCTION material_forcing_dim_names_match (fileblock, varname, dim1, dim2)
      character(len=*), intent(in) :: fileblock, varname, dim1, dim2
      integer :: ncid, varid, ndims, dimids(NF90_MAX_VAR_DIMS), status
      character(len=NF90_MAX_NAME) :: name

      material_forcing_dim_names_match = .false.
      status = nf90_open(trim(fileblock), NF90_NOWRITE, ncid)
      IF (status /= NF90_NOERR) RETURN
      status = nf90_inq_varid(ncid, trim(varname), varid)
      IF (status == NF90_NOERR) status = nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids)
      IF (status == NF90_NOERR) THEN
         IF (ndims == 2) THEN
            status = nf90_inquire_dimension(ncid, dimids(1), name=name)
            IF (status == NF90_NOERR) THEN
               IF (trim(name) == dim1) THEN
                  status = nf90_inquire_dimension(ncid, dimids(2), name=name)
                  IF (status == NF90_NOERR) material_forcing_dim_names_match = trim(name) == dim2
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      status = nf90_close(ncid)
      IF (status /= NF90_NOERR) material_forcing_dim_names_match = .false.
   END FUNCTION material_forcing_dim_names_match

   ! Called while the generic land restart transaction marker is incomplete.
   SUBROUTINE material_forcing_write_restart (file_restart)
      character(len=*), intent(in) :: file_restart
      integer, allocatable :: identity(:,:)
      integer :: iblkgrp, iblk, jblk
      character(len=256) :: fileblock

      IF (nmts <= 0) RETURN
      IF (n_mt_forc_vars > 0 .and. .not. allocated(mt_forc_precip_value)) &
         CALL CoLM_stop('material forcing cache is not allocated at restart write')
      IF (n_mt_forc_vars > 0 .and. p_is_worker) THEN
         IF (size(mt_forc_precip_value,1) /= nmts .or. &
             size(mt_forc_precip_value,2) /= landpatch%nset .or. &
             any(shape(mt_forc_vapor_value) /= shape(mt_forc_precip_value))) &
            CALL CoLM_stop('material forcing cache shape mismatch at restart write')
         IF (any(.not. ieee_is_finite(mt_forc_precip_value)) .or. &
             any(.not. ieee_is_finite(mt_forc_vapor_value))) &
            CALL CoLM_stop('non-finite material forcing cache at restart write')
      ENDIF
      IF (p_is_io) THEN
         DO iblkgrp = 1, landpatch%nblkgrp
            iblk = landpatch%xblkgrp(iblkgrp)
            jblk = landpatch%yblkgrp(iblkgrp)
            CALL get_filename_block(file_restart, iblk, jblk, fileblock)
            CALL ncio_write_serial(fileblock, 'mt_forcing_cache_schema', MT_FORC_CACHE_SCHEMA)
            CALL ncio_write_serial(fileblock, 'mt_forcing_cache_count', n_mt_forc_vars)
         ENDDO
      ENDIF
      IF (n_mt_forc_vars <= 0) RETURN
      CALL material_forcing_identity(identity)
      IF (p_is_io) THEN
         DO iblkgrp = 1, landpatch%nblkgrp
            iblk = landpatch%xblkgrp(iblkgrp)
            jblk = landpatch%yblkgrp(iblkgrp)
            CALL get_filename_block(file_restart, iblk, jblk, fileblock)
            CALL ncio_define_dimension(fileblock, 'mt_forcing_id_field', MT_FORC_ID_WIDTH)
            CALL ncio_define_dimension(fileblock, 'mt_forcing_id_var', n_mt_forc_vars+1)
            CALL ncio_define_dimension(fileblock, 'mt_forcing_species', nmts)
            CALL ncio_write_serial(fileblock, 'mt_forcing_cache_identity', identity, &
               'mt_forcing_id_field', 'mt_forcing_id_var', DEF_REST_CompressLevel)
         ENDDO
      ENDIF
      deallocate(identity)
      CALL ncio_write_vector(file_restart, 'mt_forcing_precip_last', 'mt_forcing_species', nmts, &
         'patch', landpatch, mt_forc_precip_value, DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'mt_forcing_vapor_last', 'mt_forcing_species', nmts, &
         'patch', landpatch, mt_forc_vapor_value, DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'mt_forcing_precip_reliable', 'mt_forcing_species', nmts, &
         'patch', landpatch, mt_forc_precip_reliable, DEF_REST_CompressLevel)
   END SUBROUTINE material_forcing_write_restart

   ! Restore only after mt_forcing_init has created its patch-local arrays.
   ! A hydrology-only cold start has no previous cache and keeps the defaults.
   SUBROUTINE mt_forcing_read_restart (file_restart, loaded_restart)
#ifdef USEMPI
      USE MOD_SPMD_Task, only: p_comm_glb, p_err, MPI_IN_PLACE, MPI_INTEGER, MPI_SUM
#endif
      character(len=*), intent(in) :: file_restart
      logical, intent(in) :: loaded_restart
      integer, allocatable :: expected(:,:), ondisk(:,:), varsize(:)
      real(r8), allocatable :: precip(:,:), vapor(:,:), reliable(:,:)
      integer :: iblkgrp, iblk, jblk, schema, nvars, counts(6), imt
      logical :: have_schema, have_count, block_ok
      character(len=256) :: fileblock

      IF (.not. loaded_restart .or. nmts <= 0) RETURN
      IF (n_mt_forc_vars > 0) CALL material_forcing_identity(expected)
      counts = 0
      IF (p_is_io) THEN
         DO iblkgrp = 1, landpatch%nblkgrp
            iblk = landpatch%xblkgrp(iblkgrp)
            jblk = landpatch%yblkgrp(iblkgrp)
            CALL get_filename_block(file_restart, iblk, jblk, fileblock)
            have_schema = ncio_var_exist(fileblock, 'mt_forcing_cache_schema', readflag=.false.)
            have_count = ncio_var_exist(fileblock, 'mt_forcing_cache_count', readflag=.false.)
            IF (.not. have_schema .and. .not. have_count) THEN
               counts(1) = counts(1) + 1
               CYCLE
            ENDIF
            counts(4) = counts(4) + 1
            block_ok = have_schema .and. have_count
            IF (block_ok) THEN
               CALL ncio_inquire_varsize(fileblock, 'mt_forcing_cache_schema', varsize)
               block_ok = allocated(varsize)
               IF (block_ok) THEN
                  block_ok = size(varsize) == 0
                  deallocate(varsize)
               ENDIF
               CALL ncio_inquire_varsize(fileblock, 'mt_forcing_cache_count', varsize)
               block_ok = block_ok .and. allocated(varsize)
               IF (allocated(varsize)) THEN
                  block_ok = block_ok .and. size(varsize) == 0
                  deallocate(varsize)
               ENDIF
            ENDIF
            IF (.not. block_ok) THEN
               counts(2) = counts(2) + 1
               CYCLE
            ENDIF
            CALL ncio_read_serial(fileblock, 'mt_forcing_cache_schema', schema)
            CALL ncio_read_serial(fileblock, 'mt_forcing_cache_count', nvars)
            IF ((schema /= MT_FORC_CACHE_SCHEMA .and. schema /= 1) .or. &
                nvars /= n_mt_forc_vars) THEN
               counts(3) = counts(3) + 1
               CYCLE
            ENDIF
            IF (schema == 1) THEN
               counts(5) = counts(5) + 1
            ELSE
               counts(6) = counts(6) + 1
            ENDIF
            IF (nvars <= 0) CYCLE
            IF (schema /= 1) THEN
               block_ok = ncio_var_exist(fileblock, 'mt_forcing_precip_reliable', readflag=.false.)
               IF (block_ok) THEN
                  CALL ncio_inquire_varsize(fileblock, 'mt_forcing_precip_reliable', varsize)
                  block_ok = allocated(varsize)
                  IF (block_ok) THEN
                     block_ok = size(varsize) == 2
                     IF (block_ok) block_ok = varsize(1) == nmts .and. &
                        varsize(2) == landpatch%vecgs%vlen(iblk,jblk)
                     deallocate(varsize)
                  ENDIF
               ENDIF
               IF (block_ok) block_ok = material_forcing_dim_names_match(fileblock, &
                  'mt_forcing_precip_reliable', 'mt_forcing_species', 'patch')
               IF (.not. block_ok) THEN
                  counts(2) = counts(2) + 1
                  CYCLE
               ENDIF
            ENDIF
            IF (.not. ncio_var_exist(fileblock, 'mt_forcing_cache_identity', readflag=.false.) .or. &
                .not. ncio_var_exist(fileblock, 'mt_forcing_precip_last', readflag=.false.) .or. &
                .not. ncio_var_exist(fileblock, 'mt_forcing_vapor_last', readflag=.false.)) THEN
               counts(2) = counts(2) + 1
               CYCLE
            ENDIF
            CALL ncio_inquire_varsize(fileblock, 'mt_forcing_cache_identity', varsize)
            block_ok = allocated(varsize)
            IF (block_ok) THEN
               block_ok = size(varsize) == 2
               IF (block_ok) block_ok = all(varsize == shape(expected))
               deallocate(varsize)
            ENDIF
            CALL ncio_inquire_varsize(fileblock, 'mt_forcing_precip_last', varsize)
            block_ok = block_ok .and. allocated(varsize)
            IF (allocated(varsize)) THEN
               IF (size(varsize) == 2) THEN
                  block_ok = block_ok .and. varsize(1) == nmts .and. &
                     varsize(2) == landpatch%vecgs%vlen(iblk,jblk)
               ELSE
                  block_ok = .false.
               ENDIF
               deallocate(varsize)
            ENDIF
            CALL ncio_inquire_varsize(fileblock, 'mt_forcing_vapor_last', varsize)
            block_ok = block_ok .and. allocated(varsize)
            IF (allocated(varsize)) THEN
               IF (size(varsize) == 2) THEN
                  block_ok = block_ok .and. varsize(1) == nmts .and. &
                     varsize(2) == landpatch%vecgs%vlen(iblk,jblk)
               ELSE
                  block_ok = .false.
               ENDIF
               deallocate(varsize)
            ENDIF
            IF (block_ok) block_ok = material_forcing_dim_names_match(fileblock, &
               'mt_forcing_cache_identity', 'mt_forcing_id_field', 'mt_forcing_id_var') .and. &
               material_forcing_dim_names_match(fileblock, 'mt_forcing_precip_last', &
                  'mt_forcing_species', 'patch') .and. &
               material_forcing_dim_names_match(fileblock, 'mt_forcing_vapor_last', &
                  'mt_forcing_species', 'patch')
            IF (.not. block_ok) THEN
               counts(2) = counts(2) + 1
               CYCLE
            ENDIF
            CALL ncio_read_serial(fileblock, 'mt_forcing_cache_identity', ondisk)
            IF (any(ondisk /= expected)) counts(3) = counts(3) + 1
            deallocate(ondisk)
         ENDDO
      ENDIF
#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, counts, 6, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif
      IF (counts(1) > 0 .and. counts(4) > 0) &
         CALL CoLM_stop('mixed legacy/current material forcing cache blocks')
      IF (counts(5) > 0 .and. counts(6) > 0) &
         CALL CoLM_stop('mixed schema-1/schema-2 material forcing cache blocks')
      schema = MT_FORC_CACHE_SCHEMA
      IF (counts(5) > 0) schema = 1
      IF (counts(1) > 0 .and. n_mt_forc_vars > 0) &
         CALL CoLM_stop('old material restart lacks last-valid forcing cache; cannot reconstruct it')
      IF (counts(2) > 0) CALL CoLM_stop('incomplete or malformed material forcing cache restart')
      IF (counts(3) > 0) CALL CoLM_stop('material forcing cache configuration differs from restart')
      IF (n_mt_forc_vars <= 0) RETURN
      CALL ncio_read_vector(file_restart, 'mt_forcing_precip_last', nmts, landpatch, &
         precip)
      CALL ncio_read_vector(file_restart, 'mt_forcing_vapor_last', nmts, landpatch, &
         vapor)
      IF (schema /= 1) CALL ncio_read_vector(file_restart, 'mt_forcing_precip_reliable', nmts, &
         landpatch, reliable)
      counts = 0
      IF (p_is_worker) THEN
         IF (size(mt_forc_precip_value,2) > 0 .and. &
             (.not. allocated(precip) .or. .not. allocated(vapor))) counts(2) = counts(2) + 1
         IF (allocated(precip) .and. allocated(vapor)) THEN
            IF (any(.not. ieee_is_finite(precip)) .or. any(.not. ieee_is_finite(vapor))) THEN
               counts(2) = counts(2) + 1
            ELSE
               IF (any(abs(precip) >= mt_forc_max_abs) .or. &
                   any(abs(vapor) >= mt_forc_max_abs)) counts(2) = counts(2) + 1
               DO imt = 1, nmts
                  IF (.not. mt_uses_generic_water(imt)) CYCLE
                  IF (material_is_isotope(imt)) THEN
                     IF (any(precip(imt,:) <= mt_tiny) .or. any(vapor(imt,:) <= mt_tiny)) &
                        counts(2) = counts(2) + 1
                  ELSE
                     IF (any(precip(imt,:) < 0._r8) .or. any(vapor(imt,:) < 0._r8)) &
                        counts(2) = counts(2) + 1
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
         IF (schema /= 1 .and. size(mt_forc_precip_value,2) > 0) THEN
            IF (.not. allocated(reliable)) THEN
               counts(2) = counts(2) + 1
            ELSEIF (any(reliable /= 0._r8 .and. reliable /= 1._r8)) THEN
               counts(2) = counts(2) + 1
            ENDIF
         ENDIF
      ENDIF
#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, counts, 6, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif
      IF (counts(2) > 0) CALL CoLM_stop('non-finite or invalid material forcing cache restart')
      IF (p_is_worker .and. allocated(precip)) mt_forc_precip_value = precip
      IF (p_is_worker .and. allocated(vapor)) mt_forc_vapor_value = vapor
      ! A schema-1 cache remembered every valid record: take it as reliable.
      IF (p_is_worker) THEN
         IF (schema == 1) THEN
            mt_forc_precip_reliable(:,:) = 1._r8
         ELSEIF (allocated(reliable)) THEN
            mt_forc_precip_reliable = reliable
         ENDIF
      ENDIF
   END SUBROUTINE mt_forcing_read_restart

   SUBROUTINE material_forcing_configure ()
      IMPLICIT NONE
      integer :: imt, mode, dtime, offset, nold, k
      character(len=256) :: fprefix, vname, tintalgo, token
      logical, allocatable :: has_precip(:), has_vapor(:)
      type(material_forcing_spec_type) :: spec

      ! Load each material's OWN forcing specs from its parameter file
      ! (&nl_colm_mt_forcing). This replaces the former global
      ! DEF_forcing%material_* / legacy DEF_forcing%precipitation_O18_* path.
      CALL material_forcing_input_load()
      CALL material_forcing_deallocate_config()
      allocate(mt_var_stream(2*nmts + 6))
      allocate(mt_var_imt(2*nmts + 6))
      allocate(mt_var_mode(2*nmts + 6))
      allocate(mt_var_total(2*nmts + 6))
      allocate(mt_var_dtime(2*nmts + 6))
      allocate(mt_var_offset(2*nmts + 6))
      allocate(mt_var_fprefix(2*nmts + 6))
      allocate(mt_var_vname(2*nmts + 6))
      allocate(mt_var_tintalgo(2*nmts + 6))
      allocate(mt_var_timelog(2*nmts + 6))
      allocate(mt_tstamp_LB(2*nmts + 6))
      allocate(mt_tstamp_UB(2*nmts + 6))
      allocate(has_precip(nmts), has_vapor(nmts))

      n_mt_forc_vars = 0
      idx_total_precip = 0
      idx_total_vapor = 0
      has_precip(:) = .false.
      has_vapor(:) = .false.
      mt_tstamp_LB(:) = timestamp(-1, -1, -1)
      mt_tstamp_UB(:) = timestamp(-1, -1, -1)
      mt_var_timelog(:) = 'instant'

      ! Per-species forcing: each material declares its own precip/vapor inputs
      ! in &nl_colm_mt_forcing (loaded above by MOD_MT_ForcingInput).
      ! These two find() calls are the ONLY consumers of a forcing role in the
      ! whole model, and this loop skips provider-owned species. Anything else
      ! is refused by material_forcing_role_valid rather than stored and ignored.
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         k = material_forcing_input_find(imt, 'precip')
         IF (k > 0) THEN
            spec = material_forcing_input_get(imt, k)
            mode = material_forcing_parse_mode(spec%input_mode)
            nold = n_mt_forc_vars
            CALL material_forcing_add_var(STREAM_PRECIP, imt, mode, &
               spec%fprefix, spec%vname, spec%tintalgo, spec%dtime, spec%offset)
            has_precip(imt) = n_mt_forc_vars > nold
         ENDIF

         k = material_forcing_input_find(imt, 'vapor')
         IF (k > 0) THEN
            spec = material_forcing_input_get(imt, k)
            mode = material_forcing_parse_mode(spec%input_mode)
            nold = n_mt_forc_vars
            CALL material_forcing_add_var(STREAM_VAPOR, imt, mode, &
               spec%fprefix, spec%vname, spec%tintalgo, spec%dtime, spec%offset)
            has_vapor(imt) = n_mt_forc_vars > nold
         ENDIF

         IF (material_fractionation_active(imt) .and. .not. has_vapor(imt)) THEN
            IF (p_is_master) THEN
               WRITE(*,'(A,I0,3A)') &
                  'ERROR material_forcing_configure: fractionating material ', imt, &
                  ' (', trim(materials(imt)%name), &
                  ') must define vapor forcing with role=''vapor'' in &nl_colm_mt_forcing.'
            ENDIF
            CALL CoLM_stop()
         ENDIF
      ENDDO

      mt_runtime_forcing_enabled = n_mt_forc_vars > 0
      IF (allocated(mt_runtime_forced)) THEN
         mt_runtime_forced(:) = has_precip(:) .or. has_vapor(:)
      ENDIF
      deallocate(has_precip, has_vapor)
   END SUBROUTINE material_forcing_configure

   SUBROUTINE material_forcing_deallocate_config ()
      IMPLICIT NONE

      IF (allocated(mt_var_stream)) deallocate(mt_var_stream)
      IF (allocated(mt_var_imt)) deallocate(mt_var_imt)
      IF (allocated(mt_var_mode)) deallocate(mt_var_mode)
      IF (allocated(mt_var_total)) deallocate(mt_var_total)
      IF (allocated(mt_var_dtime)) deallocate(mt_var_dtime)
      IF (allocated(mt_var_offset)) deallocate(mt_var_offset)
      IF (allocated(mt_var_fprefix)) deallocate(mt_var_fprefix)
      IF (allocated(mt_var_vname)) deallocate(mt_var_vname)
      IF (allocated(mt_var_tintalgo)) deallocate(mt_var_tintalgo)
      IF (allocated(mt_var_timelog)) deallocate(mt_var_timelog)
      IF (allocated(mt_tstamp_LB)) deallocate(mt_tstamp_LB)
      IF (allocated(mt_tstamp_UB)) deallocate(mt_tstamp_UB)
   END SUBROUTINE material_forcing_deallocate_config

   SUBROUTINE material_forcing_add_var (stream, imt, mode, fprefix, vname, tintalgo, dtime, offset)
      IMPLICIT NONE
      integer, intent(in) :: stream, imt, mode, dtime, offset
      character(len=*), intent(in) :: fprefix, vname, tintalgo
      integer :: total_idx

      IF (dtime <= 0) THEN
         IF (p_is_master) WRITE(*,'(A,I0,A,I0,A)') &
            'ERROR material_forcing_add_var: forcing dtime for material ', imt, ' must be > 0, got ', dtime, '.'
         CALL CoLM_stop()
      ENDIF
      IF (.not. material_forcing_token_present(fprefix)) RETURN
      IF (.not. material_forcing_token_present(vname)) RETURN

      total_idx = 0
      IF (mode == MODE_HEAVY_OVER_TOTAL .or. mode == MODE_NORMALIZED_OVER_TOTAL) THEN
         total_idx = material_forcing_ensure_total(stream)
         IF (dtime /= mt_var_dtime(total_idx) .or. offset /= mt_var_offset(total_idx) .or. &
             trim(material_lower(tintalgo)) /= trim(mt_var_tintalgo(total_idx))) THEN
            IF (p_is_master) WRITE(*,'(A,I0,A)') &
               'ERROR material_forcing_add_var: raw/total forcing temporal config mismatch for material ', &
               imt, '.'
            CALL CoLM_stop()
         ENDIF
      ENDIF

      n_mt_forc_vars = n_mt_forc_vars + 1
      mt_var_stream(n_mt_forc_vars) = stream
      mt_var_imt(n_mt_forc_vars) = imt
      mt_var_mode(n_mt_forc_vars) = mode
      mt_var_total(n_mt_forc_vars) = total_idx
      mt_var_fprefix(n_mt_forc_vars) = trim(fprefix)
      mt_var_vname(n_mt_forc_vars) = trim(vname)
      mt_var_tintalgo(n_mt_forc_vars) = trim(tintalgo)
      mt_var_timelog(n_mt_forc_vars) = 'instant'
      IF (total_idx > 0) mt_var_timelog(n_mt_forc_vars) = mt_var_timelog(total_idx)
      mt_var_dtime(n_mt_forc_vars) = dtime
      mt_var_offset(n_mt_forc_vars) = offset
   END SUBROUTINE material_forcing_add_var

   integer FUNCTION material_forcing_ensure_total (stream)
      IMPLICIT NONE
      integer, intent(in) :: stream

      IF (stream == STREAM_PRECIP) THEN
         IF (idx_total_precip == 0) THEN
            IF (DEF_forcing%dtime(4) <= 0) THEN
               IF (p_is_master) WRITE(*,'(A,I0,A)') &
                  'ERROR material_forcing_ensure_total: total precipitation dtime must be > 0, got ', &
                  DEF_forcing%dtime(4), '.'
               CALL CoLM_stop()
            ENDIF
            n_mt_forc_vars = n_mt_forc_vars + 1
            idx_total_precip = n_mt_forc_vars
            mt_var_stream(idx_total_precip) = STREAM_TOTAL_PRECIP
            mt_var_imt(idx_total_precip) = 0
            mt_var_mode(idx_total_precip) = MODE_VALUE
            mt_var_total(idx_total_precip) = 0
            mt_var_fprefix(idx_total_precip) = trim(DEF_forcing%fprefix(4))
            mt_var_vname(idx_total_precip) = trim(DEF_forcing%vname(4))
            mt_var_tintalgo(idx_total_precip) = trim(material_lower(DEF_forcing%tintalgo(4)))
            mt_var_timelog(idx_total_precip) = trim(material_lower(DEF_forcing%timelog(4)))
            mt_var_dtime(idx_total_precip) = DEF_forcing%dtime(4)
            mt_var_offset(idx_total_precip) = DEF_forcing%offset(4)
         ENDIF
         material_forcing_ensure_total = idx_total_precip
      ELSE
         IF (idx_total_vapor == 0) THEN
            IF (DEF_forcing%dtime(2) <= 0) THEN
               IF (p_is_master) WRITE(*,'(A,I0,A)') &
                  'ERROR material_forcing_ensure_total: total vapor dtime must be > 0, got ', &
                  DEF_forcing%dtime(2), '.'
               CALL CoLM_stop()
            ENDIF
            n_mt_forc_vars = n_mt_forc_vars + 1
            idx_total_vapor = n_mt_forc_vars
            mt_var_stream(idx_total_vapor) = STREAM_TOTAL_VAPOR
            mt_var_imt(idx_total_vapor) = 0
            mt_var_mode(idx_total_vapor) = MODE_VALUE
            mt_var_total(idx_total_vapor) = 0
            mt_var_fprefix(idx_total_vapor) = trim(DEF_forcing%fprefix(2))
            mt_var_vname(idx_total_vapor) = trim(DEF_forcing%vname(2))
            mt_var_tintalgo(idx_total_vapor) = trim(material_lower(DEF_forcing%tintalgo(2)))
            mt_var_timelog(idx_total_vapor) = trim(material_lower(DEF_forcing%timelog(2)))
            mt_var_dtime(idx_total_vapor) = DEF_forcing%dtime(2)
            mt_var_offset(idx_total_vapor) = DEF_forcing%offset(2)
         ENDIF
         material_forcing_ensure_total = idx_total_vapor
      ENDIF
   END FUNCTION material_forcing_ensure_total

   SUBROUTINE material_forcing_allocate_state (numpatch)
      IMPLICIT NONE
      integer, intent(in) :: numpatch
      integer :: imt

      CALL material_forcing_deallocate_state()
      allocate(mt_forc_precip_value(nmts, numpatch))
      allocate(mt_forc_vapor_value(nmts, numpatch))
      allocate(mt_forc_precip_step(nmts, numpatch))
      allocate(mt_forc_precip_reliable(nmts, numpatch))
      allocate(mt_forc_phase_fill(nmts))
      allocate(mt_forc_has_precip(nmts, numpatch))
      allocate(mt_forc_has_vapor(nmts, numpatch))
      allocate(mt_forc_precip_status(nmts, numpatch))
      allocate(mt_forc_vapor_status(nmts, numpatch))
      allocate(mt_forc_precip_warned_dry(nmts))
      allocate(mt_forc_precip_warned_invalid(nmts))
      allocate(mt_forc_vapor_warned_dry(nmts))
      allocate(mt_forc_vapor_warned_invalid(nmts))
      mt_forc_precip_warned_dry(:) = .false.
      mt_forc_precip_warned_invalid(:) = .false.
      mt_forc_vapor_warned_dry(:) = .false.
      mt_forc_vapor_warned_invalid(:) = .false.
      DO imt = 1, nmts
         mt_forc_precip_value(imt, :) = material_precip_default_ratio(imt)
         mt_forc_vapor_value(imt, :) = material_vapor_default_ratio(imt)
         mt_forc_precip_step(imt, :) = material_precip_default_ratio(imt)
         mt_forc_phase_fill(imt) = material_is_isotope(imt) .and. &
            material_forcing_precip_configured(imt) .and. material_forcing_vapor_configured(imt)
      ENDDO
      IF (any(mt_forc_phase_fill) .and. &
          .not. (DEF_MT_FORCING_MEMORY_PRECIP >= DEF_MT_FORCING_MIN_PRECIP)) &
         CALL CoLM_stop('DEF_MT_FORCING_MEMORY_PRECIP must be >= DEF_MT_FORCING_MIN_PRECIP for phase fill.')
      mt_forc_precip_reliable(:,:) = 0._r8
      CALL material_forcing_prepare_step()
   END SUBROUTINE material_forcing_allocate_state

   SUBROUTINE material_forcing_deallocate_state ()
      IMPLICIT NONE

      IF (allocated(mt_forc_precip_value)) deallocate(mt_forc_precip_value)
      IF (allocated(mt_forc_vapor_value)) deallocate(mt_forc_vapor_value)
      IF (allocated(mt_forc_precip_step)) deallocate(mt_forc_precip_step)
      IF (allocated(mt_forc_precip_reliable)) deallocate(mt_forc_precip_reliable)
      IF (allocated(mt_forc_phase_fill)) deallocate(mt_forc_phase_fill)
      IF (allocated(mt_forc_has_precip)) deallocate(mt_forc_has_precip)
      IF (allocated(mt_forc_has_vapor)) deallocate(mt_forc_has_vapor)
      IF (allocated(mt_forc_precip_status)) deallocate(mt_forc_precip_status)
      IF (allocated(mt_forc_vapor_status)) deallocate(mt_forc_vapor_status)
      IF (allocated(mt_forc_precip_warned_dry)) deallocate(mt_forc_precip_warned_dry)
      IF (allocated(mt_forc_precip_warned_invalid)) deallocate(mt_forc_precip_warned_invalid)
      IF (allocated(mt_forc_vapor_warned_dry)) deallocate(mt_forc_vapor_warned_dry)
      IF (allocated(mt_forc_vapor_warned_invalid)) deallocate(mt_forc_vapor_warned_invalid)
   END SUBROUTINE material_forcing_deallocate_state

   SUBROUTINE material_forcing_prepare_step ()
      IMPLICIT NONE

      IF (.not. allocated(mt_forc_precip_value)) RETURN

      mt_forc_has_precip(:,:) = .false.
      mt_forc_has_vapor(:,:) = .false.
      mt_forc_precip_status(:,:) = DECODE_INVALID
      mt_forc_vapor_status(:,:) = DECODE_INVALID
   END SUBROUTINE material_forcing_prepare_step

   SUBROUTINE material_forcing_read_LBUB (idate, dir_forcing)
      IMPLICIT NONE
      integer, intent(in) :: idate(3)
      character(len=*), intent(in) :: dir_forcing
      integer :: iv, year, month, day, time_i
      type(timestamp) :: mtstamp
      character(len=256) :: filename

      mtstamp = idate

      DO iv = 1, n_mt_forc_vars
         IF (.not.(mt_tstamp_LB(iv) == 'NULL') .and. .not.(mt_tstamp_UB(iv) == 'NULL') .and. &
             mt_tstamp_LB(iv) <= mtstamp .and. mtstamp < mt_tstamp_UB(iv)) CYCLE

         IF (mt_tstamp_LB(iv) == 'NULL') THEN
            CALL material_forcing_setstamp_LB(mtstamp, iv, year, month, day, time_i)
            filename = trim(dir_forcing)//trim(material_forcing_filename(year, month, day, iv))
            IF (p_is_io) THEN
               CALL ncio_read_block_time(filename, trim(mt_var_vname(iv)), mt_gforc, time_i, mt_metdata)
               CALL block_data_copy(mt_metdata, mt_forcn_LB(iv))
            ENDIF
         ENDIF

         DO WHILE (mt_tstamp_UB(iv) == 'NULL' .or. mt_tstamp_UB(iv) <= mtstamp)
            IF (.not. (mt_tstamp_UB(iv) == 'NULL')) THEN
               IF (p_is_io) CALL block_data_copy(mt_forcn_UB(iv), mt_forcn_LB(iv))
            ENDIF

            CALL material_forcing_setstamp_UB(iv, year, month, day, time_i)
            filename = trim(dir_forcing)//trim(material_forcing_filename(year, month, day, iv))
            IF (p_is_io) THEN
               CALL ncio_read_block_time(filename, trim(mt_var_vname(iv)), mt_gforc, time_i, mt_metdata)
               CALL block_data_copy(mt_metdata, mt_forcn_UB(iv))
            ENDIF
         ENDDO
      ENDDO
   END SUBROUTINE material_forcing_read_LBUB

   SUBROUTINE material_forcing_interpolate_var (idate, iv)
      IMPLICIT NONE
      integer, intent(in) :: idate(3), iv
      type(timestamp) :: mtstamp
      integer :: dtLB, dtUB

      mtstamp = idate
      IF ((mtstamp < mt_tstamp_LB(iv)) .or. (mt_tstamp_UB(iv) < mtstamp)) THEN
         WRITE(6, *) 'material forcing data required is out of range! STOP!'
         CALL CoLM_stop()
      ENDIF

      dtLB = mtstamp - mt_tstamp_LB(iv)
      dtUB = mt_tstamp_UB(iv) - mtstamp

      IF (trim(mt_var_tintalgo(iv)) == 'linear') THEN
         IF ((dtLB + dtUB) > 0) THEN
            CALL block_data_linear_interp( &
               mt_forcn_LB(iv), real(dtUB, r8) / real(dtLB + dtUB, r8), &
               mt_forcn_UB(iv), real(dtLB, r8) / real(dtLB + dtUB, r8), &
               mt_forcn(iv))
         ELSE
            CALL block_data_copy(mt_forcn_LB(iv), mt_forcn(iv))
         ENDIF
      ELSEIF (trim(mt_var_tintalgo(iv)) == 'nearest') THEN
         IF (dtLB <= dtUB) THEN
            CALL block_data_copy(mt_forcn_LB(iv), mt_forcn(iv))
         ELSE
            CALL block_data_copy(mt_forcn_UB(iv), mt_forcn(iv))
         ENDIF
      ELSEIF (trim(mt_var_tintalgo(iv)) == 'uniform') THEN
         ! Match MOD_Forcing exactly: interval means logged forward belong
         ! to the lower timestamp; backward (and other non-forward logs)
         ! belong to the upper timestamp.
         IF (trim(mt_var_timelog(iv)) == 'forward') THEN
            CALL block_data_copy(mt_forcn_LB(iv), mt_forcn(iv))
         ELSE
            CALL block_data_copy(mt_forcn_UB(iv), mt_forcn(iv))
         ENDIF
      ELSE
         IF (p_is_master) WRITE(*,'(A,I0,3A)') &
            'ERROR material_forcing_interpolate_var: variable ', iv, &
            ' has invalid interpolation mode "', trim(mt_var_tintalgo(iv)), '".'
         CALL CoLM_stop()
      ENDIF
   END SUBROUTINE material_forcing_interpolate_var

   SUBROUTINE material_forcing_update_values ()
      IMPLICIT NONE
      integer :: iv, ip, imt, decode_status
      real(r8) :: value
      logical :: value_valid

      IF (.not. p_is_worker) RETURN
      IF (.not. allocated(mt_forc_patch)) RETURN
      IF (.not. allocated(mt_forc_precip_value)) RETURN

      DO iv = 1, n_mt_forc_vars
         imt = mt_var_imt(iv)
         IF (imt < 1 .or. imt > nmts) CYCLE
         DO ip = 1, size(mt_forc_patch, 1)
            CALL material_forcing_decode_value(iv, ip, value, value_valid, decode_status)
            IF (mt_var_stream(iv) == STREAM_PRECIP) mt_forc_precip_status(imt, ip) = decode_status
            IF (mt_var_stream(iv) == STREAM_VAPOR) mt_forc_vapor_status(imt, ip) = decode_status
            IF (.not. value_valid) CYCLE
            IF (mt_var_stream(iv) == STREAM_PRECIP) THEN
               mt_forc_precip_step(imt, ip) = value
               mt_forc_has_precip(imt, ip) = .true.
               IF (material_forcing_reliable_record(iv, ip)) THEN
                  mt_forc_precip_value(imt, ip) = value
                  mt_forc_precip_reliable(imt, ip) = 1._r8
               ENDIF
            ELSEIF (mt_var_stream(iv) == STREAM_VAPOR) THEN
               mt_forc_vapor_value(imt, ip) = value
               mt_forc_has_vapor(imt, ip) = .true.
            ENDIF
         ENDDO
      ENDDO

      ! Until a patch has seen a reliable record (cold start, or a new patch
      ! after LULCC), its remembered snow ratio is the ratio of ice deposited
      ! from the forcing vapour (Jouzel-Merlivat effective alpha at the air
      ! temperature), updated each step, instead of the fixed default ratio.
      DO imt = 1, nmts
         IF (.not. mt_forc_phase_fill(imt)) CYCLE
         DO ip = 1, size(mt_forc_precip_value, 2)
            IF (mt_forc_precip_reliable(imt, ip) > 0.5_r8) CYCLE
            value = material_alpha_ice_vap_deposition(imt, forc_t(ip)) * mt_forc_vapor_value(imt, ip)
            IF (material_forcing_valid_value(value) .and. value > mt_tiny) &
               mt_forc_precip_value(imt, ip) = value
         ENDDO
      ENDDO
   END SUBROUTINE material_forcing_update_values

   ! A valid precipitation record updates the remembered ratio only when it is
   ! reliable.  With phase fill (isotopes with precipitation and vapour
   ! forcing), a heavy/total ratio is reliable above
   ! DEF_MT_FORCING_MEMORY_PRECIP; without it every valid record is
   ! remembered, as before.
   logical FUNCTION material_forcing_reliable_record (iv, ip)
      IMPLICIT NONE
      integer, intent(in) :: iv, ip
      integer :: imt

      material_forcing_reliable_record = .true.
      imt = mt_var_imt(iv)
      IF (.not. mt_forc_phase_fill(imt)) RETURN
      IF (mt_var_mode(iv) /= MODE_HEAVY_OVER_TOTAL .and. &
          mt_var_mode(iv) /= MODE_NORMALIZED_OVER_TOTAL) RETURN
      IF (mt_var_total(iv) <= 0) RETURN
      material_forcing_reliable_record = &
         mt_forc_patch(ip, mt_var_total(iv)) > DEF_MT_FORCING_MEMORY_PRECIP
   END FUNCTION material_forcing_reliable_record

   SUBROUTINE material_forcing_decode_value (iv, ip, value, valid, status)
      IMPLICIT NONE
      integer, intent(in) :: iv, ip
      real(r8), intent(out) :: value
      logical, intent(out) :: valid
      integer, intent(out), optional :: status
      integer :: imt, total_idx, decode_status
      real(r8) :: raw, total, min_total, delta

      value = 0._r8
      valid = .false.
      decode_status = DECODE_INVALID
      raw = mt_forc_patch(ip, iv)

      imt = mt_var_imt(iv)
      SELECT CASE (mt_var_mode(iv))
      CASE (MODE_HEAVY_OVER_TOTAL, MODE_NORMALIZED_OVER_TOTAL)
         total_idx = mt_var_total(iv)
         IF (total_idx <= 0) THEN
            IF (present(status)) status = decode_status
            RETURN
         ENDIF
         total = mt_forc_patch(ip, total_idx)
         IF (.not. material_forcing_valid_value(total)) THEN
            IF (present(status)) status = decode_status
            RETURN
         ENDIF
         min_total = mt_forc_min_q
         IF (mt_var_stream(iv) == STREAM_PRECIP) min_total = DEF_MT_FORCING_MIN_PRECIP
         IF (total <= min_total) THEN
            decode_status = DECODE_NEAR_DRY
            IF (present(status)) status = decode_status
            RETURN
         ENDIF
         IF (.not. material_forcing_valid_value(raw)) THEN
            IF (present(status)) status = decode_status
            RETURN
         ENDIF
         IF (raw < 0._r8) THEN
            IF (present(status)) status = decode_status
            RETURN
         ENDIF
         value = raw / total
         IF (mt_var_mode(iv) == MODE_NORMALIZED_OVER_TOTAL) THEN
            ! `normalized_over_total` means the heavy/total stream is
            ! normalized by the material reference value. This is not limited
            ! to category=isotope: a no-fractionation HDO test
            ! still needs the IsoGSM normalized ratio converted back to the
            ! stored concentration per unit water. For solutes
            ! the default ref_ratio=1 keeps this identical to
            ! heavy_over_total.
            value = value * materials(imt)%ref_ratio
         ENDIF
      CASE (MODE_DELTA)
         IF (.not. material_forcing_valid_value(raw)) THEN
            IF (present(status)) status = decode_status
            RETURN
         ENDIF
         IF (material_is_isotope(imt)) THEN
            value = delta_to_R(raw, materials(imt)%ref_ratio)
         ELSE
            value = raw
         ENDIF
      CASE DEFAULT
         IF (.not. material_forcing_valid_value(raw)) THEN
            IF (present(status)) status = decode_status
            RETURN
         ENDIF
         value = raw
      END SELECT

      IF (material_is_isotope(imt)) THEN
         IF (value <= mt_tiny) THEN
            IF (present(status)) status = decode_status
            RETURN
         ENDIF
         delta = material_forcing_ratio_to_delta(value, materials(imt)%ref_ratio)
         IF (abs(delta) > mt_delta_sanity_max) THEN
            IF (present(status)) status = decode_status
            RETURN
         ENDIF
      ELSE
         ! Solute forcings are concentrations per unit water;
         ! zero is a valid concentration, negative values are not.
         IF (value < 0._r8) THEN
            IF (present(status)) status = decode_status
            RETURN
         ENDIF
      ENDIF
      valid = .true.
      decode_status = DECODE_VALID
      IF (present(status)) status = decode_status
   END SUBROUTINE material_forcing_decode_value

   logical FUNCTION material_forcing_valid_value (x)
      IMPLICIT NONE
      real(r8), intent(in) :: x

      material_forcing_valid_value = abs(x) < mt_forc_max_abs
   END FUNCTION material_forcing_valid_value

   SUBROUTINE material_forcing_setstamp_LB (mtstamp, iv, year, month, mday, time_i)
      IMPLICIT NONE
      type(timestamp), intent(in) :: mtstamp
      integer, intent(in) :: iv
      integer, intent(out) :: year, month, mday, time_i
      integer :: day, sec, sec_file
      integer :: months(0:12)

      CALL material_forcing_validate_var_dtime(iv)

      year = mtstamp%year
      day = mtstamp%day
      sec = mtstamp%sec

      mt_tstamp_LB(iv)%year = year
      mt_tstamp_LB(iv)%day = day

      IF (trim(DEF_forcing%groupby) == 'year') THEN
         sec_file = 86400 * (day - 1) + sec
         time_i = floor((sec_file - mt_var_offset(iv)) * 1.0_r8 / mt_var_dtime(iv)) + 1
         sec = (time_i - 1) * mt_var_dtime(iv) + mt_var_offset(iv) - 86400 * (day - 1)
         mt_tstamp_LB(iv)%sec = sec

         IF (sec < 0) THEN
            mt_tstamp_LB(iv)%sec = 86400 + sec
            mt_tstamp_LB(iv)%day = day - 1
            IF (mt_tstamp_LB(iv)%day == 0) THEN
               mt_tstamp_LB(iv)%year = year - 1
               IF (isleapyear(mt_tstamp_LB(iv)%year)) THEN
                  mt_tstamp_LB(iv)%day = 366
               ELSE
                  mt_tstamp_LB(iv)%day = 365
               ENDIF
            ENDIF
         ENDIF

         IF (sec < 0 .or. (sec == 0 .and. mt_var_offset(iv) /= 0)) THEN
            IF (year == DEF_forcing%startyr .and. DEF_forcing%startmo == 1 .and. day == 1) THEN
               sec = mt_var_offset(iv)
            ELSE
               sec = 86400 + sec
               day = day - 1
               IF (day == 0) THEN
                  year = year - 1
                  IF (isleapyear(year)) THEN
                     day = 366
                  ELSE
                     day = 365
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         IF (.not. DEF_forcing%leapyear .and. isleapyear(year) .and. day > 59) day = day - 1

         sec_file = 86400 * (day - 1) + sec
         time_i = floor((sec_file - mt_var_offset(iv)) * 1.0_r8 / mt_var_dtime(iv)) + 1
         CALL julian2monthday(year, day, month, mday)

      ELSEIF (trim(DEF_forcing%groupby) == 'month') THEN
         IF (isleapyear(year)) THEN
            months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
         ELSE
            months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
         ENDIF

         CALL julian2monthday(year, day, month, mday)
         sec_file = 86400 * (mday - 1) + sec
         time_i = floor((sec_file - mt_var_offset(iv)) * 1.0_r8 / mt_var_dtime(iv)) + 1
         sec = (time_i - 1) * mt_var_dtime(iv) + mt_var_offset(iv) - 86400 * (mday - 1)
         mt_tstamp_LB(iv)%sec = sec

         IF (sec < 0) THEN
            mt_tstamp_LB(iv)%sec = 86400 + sec
            mt_tstamp_LB(iv)%day = day - 1
            IF (mt_tstamp_LB(iv)%day == 0) THEN
               mt_tstamp_LB(iv)%year = year - 1
               IF (isleapyear(mt_tstamp_LB(iv)%year)) THEN
                  mt_tstamp_LB(iv)%day = 366
               ELSE
                  mt_tstamp_LB(iv)%day = 365
               ENDIF
            ENDIF
         ENDIF

         IF (sec < 0 .or. (sec == 0 .and. mt_var_offset(iv) /= 0)) THEN
            IF (year == DEF_forcing%startyr .and. month == DEF_forcing%startmo .and. mday == 1) THEN
               sec = mt_var_offset(iv)
            ELSE
               sec = 86400 + sec
               mday = mday - 1
               IF (mday == 0) THEN
                  month = month - 1
                  IF (month == 0) THEN
                     month = 12
                     year = year - 1
                     mday = 31
                  ELSE
                     mday = months(month) - months(month - 1)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         IF (.not. DEF_forcing%leapyear .and. isleapyear(year) .and. month == 2 .and. mday == 29) mday = 28

         sec_file = 86400 * (mday - 1) + sec
         time_i = floor((sec_file - mt_var_offset(iv)) * 1.0_r8 / mt_var_dtime(iv)) + 1

      ELSEIF (trim(DEF_forcing%groupby) == 'day') THEN
         CALL julian2monthday(year, day, month, mday)
         time_i = floor((sec - mt_var_offset(iv)) * 1.0_r8 / mt_var_dtime(iv)) + 1
         sec = (time_i - 1) * mt_var_dtime(iv) + mt_var_offset(iv)
         mt_tstamp_LB(iv)%sec = sec

         IF (sec < 0) THEN
            mt_tstamp_LB(iv)%sec = 86400 + sec
            mt_tstamp_LB(iv)%day = day - 1
            IF (mt_tstamp_LB(iv)%day == 0) THEN
               mt_tstamp_LB(iv)%year = year - 1
               IF (isleapyear(mt_tstamp_LB(iv)%year)) THEN
                  mt_tstamp_LB(iv)%day = 366
               ELSE
                  mt_tstamp_LB(iv)%day = 365
               ENDIF
            ENDIF

            IF (year == DEF_forcing%startyr .and. month == DEF_forcing%startmo .and. mday == 1) THEN
               sec = mt_var_offset(iv)
            ELSE
               sec = 86400 + sec
               year = mt_tstamp_LB(iv)%year
               CALL julian2monthday(mt_tstamp_LB(iv)%year, mt_tstamp_LB(iv)%day, month, mday)
            ENDIF
         ENDIF

         IF (.not. DEF_forcing%leapyear .and. isleapyear(year) .and. month == 2 .and. mday == 29) mday = 28
         time_i = floor((sec - mt_var_offset(iv)) * 1.0_r8 / mt_var_dtime(iv)) + 1

      ELSE
         IF (p_is_master) WRITE(*,'(2A)') 'ERROR material forcing unsupported groupby: ', trim(DEF_forcing%groupby)
         CALL CoLM_stop()
      ENDIF

      IF (time_i <= 0) THEN
         WRITE(6, *) 'got the wrong time record of material forcing! STOP!'
         CALL CoLM_stop()
      ENDIF
   END SUBROUTINE material_forcing_setstamp_LB

   SUBROUTINE material_forcing_setstamp_UB (iv, year, month, mday, time_i)
      IMPLICIT NONE
      integer, intent(in) :: iv
      integer, intent(out) :: year, month, mday, time_i
      integer :: day, sec, sec_file
      integer :: months(0:12)

      CALL material_forcing_validate_var_dtime(iv)

      IF (mt_tstamp_UB(iv) == 'NULL') THEN
         mt_tstamp_UB(iv) = mt_tstamp_LB(iv) + mt_var_dtime(iv)
      ELSE
         mt_tstamp_LB(iv) = mt_tstamp_UB(iv)
         mt_tstamp_UB(iv) = mt_tstamp_UB(iv) + mt_var_dtime(iv)
      ENDIF

      year = mt_tstamp_UB(iv)%year
      day = mt_tstamp_UB(iv)%day
      sec = mt_tstamp_UB(iv)%sec

      IF (trim(DEF_forcing%groupby) == 'year') THEN
         IF (sec == 86400 .and. mt_var_offset(iv) == 0) THEN
            sec = 0
            day = day + 1
            IF (isleapyear(year) .and. day == 367) THEN
               year = year + 1
               day = 1
            ENDIF
            IF ((.not. isleapyear(year)) .and. day == 366) THEN
               year = year + 1
               day = 1
            ENDIF
         ENDIF

         IF (.not. DEF_forcing%leapyear .and. isleapyear(year) .and. day > 59) day = day - 1

         sec_file = 86400 * (day - 1) + sec
         time_i = floor((sec_file - mt_var_offset(iv)) * 1.0_r8 / mt_var_dtime(iv)) + 1
         CALL julian2monthday(year, day, month, mday)

      ELSEIF (trim(DEF_forcing%groupby) == 'month') THEN
         IF (isleapyear(year)) THEN
            months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
         ELSE
            months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
         ENDIF

         CALL julian2monthday(year, day, month, mday)

         IF (sec == 86400 .and. mt_var_offset(iv) == 0) THEN
            sec = 0
            mday = mday + 1
            IF (mday > (months(month) - months(month - 1))) THEN
               mday = 1
               IF (month == 12) THEN
                  month = 1
                  year = year + 1
               ELSE
                  month = month + 1
               ENDIF
            ENDIF
         ENDIF

         IF (.not. DEF_forcing%leapyear .and. isleapyear(year) .and. month == 2 .and. mday == 29) mday = 28

         sec_file = 86400 * (mday - 1) + sec
         time_i = floor((sec_file - mt_var_offset(iv)) * 1.0_r8 / mt_var_dtime(iv)) + 1

      ELSEIF (trim(DEF_forcing%groupby) == 'day') THEN
         IF (isleapyear(year)) THEN
            months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
         ELSE
            months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
         ENDIF

         CALL julian2monthday(year, day, month, mday)

         IF (sec == 86400 .and. mt_var_offset(iv) == 0) THEN
            sec = 0
            mday = mday + 1
            IF (mday > (months(month) - months(month - 1))) THEN
               mday = 1
               IF (month == 12) THEN
                  month = 1
                  year = year + 1
               ELSE
                  month = month + 1
               ENDIF
            ENDIF
         ENDIF

         IF (.not. DEF_forcing%leapyear .and. isleapyear(year) .and. month == 2 .and. mday == 29) mday = 28
         time_i = floor((sec - mt_var_offset(iv)) * 1.0_r8 / mt_var_dtime(iv)) + 1

      ELSE
         IF (p_is_master) WRITE(*,'(2A)') 'ERROR material forcing unsupported groupby: ', trim(DEF_forcing%groupby)
         CALL CoLM_stop()
      ENDIF

      IF (time_i <= 0) THEN
         WRITE(6, *) 'got the wrong time record of material forcing! STOP!'
         CALL CoLM_stop()
      ENDIF
   END SUBROUTINE material_forcing_setstamp_UB

   SUBROUTINE material_forcing_validate_var_dtime (iv)
      IMPLICIT NONE
      integer, intent(in) :: iv

      IF (.not. allocated(mt_var_dtime)) THEN
         IF (p_is_master) WRITE(*,'(A)') 'ERROR material forcing dtime is not configured.'
         CALL CoLM_stop()
      ENDIF
      IF (iv < 1 .or. iv > n_mt_forc_vars) THEN
         IF (p_is_master) WRITE(*,'(A,I0,A)') 'ERROR material forcing variable index ', iv, ' is out of range.'
         CALL CoLM_stop()
      ENDIF
      IF (mt_var_dtime(iv) <= 0) THEN
         IF (p_is_master) WRITE(*,'(A,I0,A,I0,A)') 'ERROR material forcing dtime for variable ', iv, &
            ' must be > 0, got ', mt_var_dtime(iv), '.'
         CALL CoLM_stop()
      ENDIF
   END SUBROUTINE material_forcing_validate_var_dtime

   character(len=256) FUNCTION material_forcing_filename (year, month, day, iv)
      IMPLICIT NONE
      integer, intent(in) :: year, month, day, iv
      integer :: main_var
      character(len=16) :: yearstr, monthstr, daystr

      main_var = material_forcing_main_var_index(iv)
      IF (main_var > 0 .and. trim(mt_var_fprefix(iv)) == trim(DEF_forcing%fprefix(main_var))) THEN
         material_forcing_filename = metfilename(year, month, day, main_var)
         RETURN
      ENDIF

      write(yearstr, '(I4.4)') year
      write(monthstr, '(I2.2)') month
      write(daystr, '(I2.2)') day

      IF (trim(DEF_forcing%groupby) == 'year') THEN
         material_forcing_filename = '/'//trim(mt_var_fprefix(iv))//'_'//trim(yearstr)//'.nc'
      ELSEIF (trim(DEF_forcing%groupby) == 'month') THEN
         material_forcing_filename = '/'//trim(mt_var_fprefix(iv))//'_'//trim(yearstr)//'_'//trim(monthstr)//'.nc'
      ELSEIF (trim(DEF_forcing%groupby) == 'day') THEN
         material_forcing_filename = '/'//trim(mt_var_fprefix(iv))//'_'//trim(yearstr)//'_'//trim(monthstr)// &
            '_'//trim(daystr)//'.nc'
      ELSE
         material_forcing_filename = '/'//trim(mt_var_fprefix(iv))//'_'//trim(yearstr)//'.nc'
      ENDIF
   END FUNCTION material_forcing_filename

   integer FUNCTION material_forcing_main_var_index (iv)
      IMPLICIT NONE
      integer, intent(in) :: iv

      material_forcing_main_var_index = 0
      SELECT CASE (mt_var_stream(iv))
      CASE (STREAM_PRECIP, STREAM_TOTAL_PRECIP)
         material_forcing_main_var_index = 4
      CASE (STREAM_VAPOR, STREAM_TOTAL_VAPOR)
         material_forcing_main_var_index = 2
      END SELECT
   END FUNCTION material_forcing_main_var_index

   real(r8) FUNCTION material_forcing_precip_value (imt, ipatch)
      IMPLICIT NONE
      integer, intent(in) :: imt, ipatch

      material_forcing_precip_value = material_precip_default_ratio(imt)
      IF (.not. allocated(mt_forc_precip_value)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN
      IF (ipatch < 1 .or. ipatch > size(mt_forc_precip_value, 2)) RETURN
      material_forcing_precip_value = mt_forc_precip_value(imt, ipatch)
   END FUNCTION material_forcing_precip_value

   !---------------------------------------------------------------
   ! Precipitation ratio by phase.
   !
   ! Where the current record is valid (total above DEF_MT_FORCING_MIN_PRECIP)
   ! rain and snow both take it.  Where it is near-dry, the heavy/total ratio
   ! of IsoGSM is noise, and the two phases are filled differently
   ! (IsoGSM 2004, 6-hourly, predicting the ratio of records above 2e-5
   ! kg m-2 s-1, d18O bias / RMSE in permil):
   !
   ! - rain: in isotopic equilibrium with the forcing vapour at the air
   !   temperature (falling drops exchange with the surface air): NEC May
   !   -0.6 / 2.3, NEC Jul +0.5 / 1.5, AMZ +0.2..+0.6 / 0.9; the last valid
   !   ratio gives 0.9..5.4 RMSE.
   ! - snow: the ratio of the last reliable record (total above
   !   DEF_MT_FORCING_MEMORY_PRECIP): NEC Jan +0.1 / 4.5, 40-60N Jan +0.2 /
   !   4.5; ice in equilibrium with the vapour is +6.6..+8.7 too heavy (snow
   !   forms aloft from colder, more depleted vapour).  Before any reliable
   !   record the remembered ratio is the deposition-equilibrium estimate
   !   (material_forcing_update_values).
   !
   ! The fixed default ratio of the old fallback was the worst (NEC Jan
   ! +18, AMZ -7).  Species without phase fill (not an isotope, or no vapour
   ! forcing) keep the last valid ratio for both phases, as before.
   ! The alphas describe the atmosphere that made the forcing, so they do not
   ! depend on DEF_ISO_USE_FRACTIONATION.
   !---------------------------------------------------------------
   real(r8) FUNCTION material_forcing_rain_value (imt, ipatch)
      IMPLICIT NONE
      integer, intent(in) :: imt, ipatch
      real(r8) :: value

      material_forcing_rain_value = material_forcing_precip_value(imt, ipatch)
      IF (.not. material_forcing_phase_fill_at(imt, ipatch)) RETURN
      IF (mt_forc_has_precip(imt, ipatch)) THEN
         material_forcing_rain_value = mt_forc_precip_step(imt, ipatch)
      ELSEIF (mt_forc_precip_status(imt, ipatch) == DECODE_NEAR_DRY) THEN
         value = material_alpha_liq_vap(imt, forc_t(ipatch)) * mt_forc_vapor_value(imt, ipatch)
         IF (material_forcing_valid_value(value) .and. value > mt_tiny) &
            material_forcing_rain_value = value
      ENDIF
   END FUNCTION material_forcing_rain_value

   real(r8) FUNCTION material_forcing_snow_value (imt, ipatch)
      IMPLICIT NONE
      integer, intent(in) :: imt, ipatch

      material_forcing_snow_value = material_forcing_precip_value(imt, ipatch)
      IF (.not. material_forcing_phase_fill_at(imt, ipatch)) RETURN
      IF (mt_forc_has_precip(imt, ipatch)) &
         material_forcing_snow_value = mt_forc_precip_step(imt, ipatch)
   END FUNCTION material_forcing_snow_value

   logical FUNCTION material_forcing_phase_fill_at (imt, ipatch)
      IMPLICIT NONE
      integer, intent(in) :: imt, ipatch

      material_forcing_phase_fill_at = .false.
      IF (.not. allocated(mt_forc_phase_fill)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN
      IF (ipatch < 1 .or. ipatch > size(mt_forc_precip_value, 2)) RETURN
      material_forcing_phase_fill_at = mt_forc_phase_fill(imt)
   END FUNCTION material_forcing_phase_fill_at

   real(r8) FUNCTION material_forcing_vapor_value (imt, ipatch)
      IMPLICIT NONE
      integer, intent(in) :: imt, ipatch

      material_forcing_vapor_value = material_vapor_default_ratio(imt)
      IF (.not. allocated(mt_forc_vapor_value)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN
      IF (ipatch < 1 .or. ipatch > size(mt_forc_vapor_value, 2)) RETURN
      material_forcing_vapor_value = mt_forc_vapor_value(imt, ipatch)
   END FUNCTION material_forcing_vapor_value

   real(r8) FUNCTION material_forcing_precip_ratio (imt, ipatch)
      IMPLICIT NONE
      integer, intent(in) :: imt, ipatch
      material_forcing_precip_ratio = material_forcing_precip_value(imt, ipatch)
   END FUNCTION material_forcing_precip_ratio

   real(r8) FUNCTION material_forcing_vapor_ratio (imt, ipatch)
      IMPLICIT NONE
      integer, intent(in) :: imt, ipatch
      material_forcing_vapor_ratio = material_forcing_vapor_value(imt, ipatch)
   END FUNCTION material_forcing_vapor_ratio

   logical FUNCTION material_forcing_has_precip (imt, ipatch)
      IMPLICIT NONE
      integer, intent(in) :: imt, ipatch

      material_forcing_has_precip = .false.
      IF (.not. allocated(mt_forc_has_precip)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN
      IF (ipatch < 1 .or. ipatch > size(mt_forc_has_precip, 2)) RETURN
      material_forcing_has_precip = mt_forc_has_precip(imt, ipatch) .or. &
         material_forcing_precip_configured(imt)
   END FUNCTION material_forcing_has_precip

   logical FUNCTION material_forcing_has_vapor (imt, ipatch)
      IMPLICIT NONE
      integer, intent(in) :: imt, ipatch

      material_forcing_has_vapor = .false.
      IF (.not. allocated(mt_forc_has_vapor)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN
      IF (ipatch < 1 .or. ipatch > size(mt_forc_has_vapor, 2)) RETURN
      material_forcing_has_vapor = mt_forc_has_vapor(imt, ipatch) .or. &
         material_forcing_vapor_configured(imt)
   END FUNCTION material_forcing_has_vapor


   logical FUNCTION material_forcing_precip_configured (imt)
      IMPLICIT NONE
      integer, intent(in) :: imt
      integer :: iv

      material_forcing_precip_configured = .false.
      DO iv = 1, n_mt_forc_vars
         IF (mt_var_stream(iv) == STREAM_PRECIP .and. mt_var_imt(iv) == imt) THEN
            material_forcing_precip_configured = .true.
            RETURN
         ENDIF
      ENDDO
   END FUNCTION material_forcing_precip_configured

   logical FUNCTION material_forcing_vapor_configured (imt)
      IMPLICIT NONE
      integer, intent(in) :: imt
      integer :: iv

      material_forcing_vapor_configured = .false.
      DO iv = 1, n_mt_forc_vars
         IF (mt_var_stream(iv) == STREAM_VAPOR .and. mt_var_imt(iv) == imt) THEN
            material_forcing_vapor_configured = .true.
            RETURN
         ENDIF
      ENDDO
   END FUNCTION material_forcing_vapor_configured

   real(r8) FUNCTION material_forcing_ratio_to_delta (ratio, ref_ratio)
      IMPLICIT NONE
      real(r8), intent(in) :: ratio, ref_ratio

      IF (ratio > mt_tiny .and. ref_ratio > mt_tiny) THEN
         material_forcing_ratio_to_delta = (ratio / ref_ratio - 1._r8) * 1000._r8
      ELSE
         material_forcing_ratio_to_delta = 0._r8
      ENDIF
   END FUNCTION material_forcing_ratio_to_delta

   SUBROUTINE material_forcing_log_ranges (idate)
      IMPLICIT NONE
      integer, intent(in) :: idate(3)
      integer :: imt, ip
      real(r8), allocatable :: pmin(:), pmax(:), vmin(:), vmax(:)
      integer, allocatable :: pcnt(:), vcnt(:), pmiss(:), pdry(:), vmiss(:), vdry(:)
      real(r8) :: outval
      logical :: emit_header

      IF (.not. allocated(mt_forc_precip_value)) RETURN

      allocate(pmin(nmts), pmax(nmts), vmin(nmts), vmax(nmts))
      allocate(pcnt(nmts), vcnt(nmts), pmiss(nmts), pdry(nmts), vmiss(nmts), vdry(nmts))
      pmin(:) = huge(1._r8)
      vmin(:) = huge(1._r8)
      pmax(:) = -huge(1._r8)
      vmax(:) = -huge(1._r8)
      pcnt(:) = 0
      vcnt(:) = 0
      pmiss(:) = 0
      pdry(:) = 0
      vmiss(:) = 0
      vdry(:) = 0

      IF (p_is_worker) THEN
         DO imt = 1, nmts
            IF (.not. mt_uses_generic_water(imt)) CYCLE
            DO ip = 1, size(mt_forc_precip_value, 2)
               IF (mt_forc_has_precip(imt, ip)) THEN
                  outval = material_forcing_diag_value(imt, mt_forc_precip_step(imt, ip))
                  pmin(imt) = min(pmin(imt), outval)
                  pmax(imt) = max(pmax(imt), outval)
                  pcnt(imt) = pcnt(imt) + 1
               ELSEIF (material_is_isotope(imt) .and. material_forcing_precip_configured(imt)) THEN
                  IF (allocated(mt_forc_precip_status) .and. &
                      mt_forc_precip_status(imt, ip) == DECODE_NEAR_DRY) THEN
                     pdry(imt) = pdry(imt) + 1
                  ELSE
                     pmiss(imt) = pmiss(imt) + 1
                  ENDIF
               ENDIF
               IF (mt_forc_has_vapor(imt, ip)) THEN
                  outval = material_forcing_diag_value(imt, mt_forc_vapor_value(imt, ip))
                  vmin(imt) = min(vmin(imt), outval)
                  vmax(imt) = max(vmax(imt), outval)
                  vcnt(imt) = vcnt(imt) + 1
               ELSEIF (material_fractionation_active(imt)) THEN
                  IF (allocated(mt_forc_vapor_status) .and. &
                      mt_forc_vapor_status(imt, ip) == DECODE_NEAR_DRY) THEN
                     vdry(imt) = vdry(imt) + 1
                  ELSE
                     vmiss(imt) = vmiss(imt) + 1
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDIF

#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, pmin, nmts, MPI_REAL8, MPI_MIN, p_comm_glb, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, pmax, nmts, MPI_REAL8, MPI_MAX, p_comm_glb, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, vmin, nmts, MPI_REAL8, MPI_MIN, p_comm_glb, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, vmax, nmts, MPI_REAL8, MPI_MAX, p_comm_glb, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, pcnt, nmts, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, vcnt, nmts, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, pmiss, nmts, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, pdry, nmts, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, vmiss, nmts, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, vdry, nmts, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN
         emit_header = mt_forcing_log_count < 3
         DO imt = 1, nmts
            IF (.not. mt_uses_generic_water(imt)) CYCLE
            IF (pdry(imt) > 0 .and. &
                (mt_forcing_log_count < 3 .or. .not. mt_forc_precip_warned_dry(imt))) &
               emit_header = .true.
            IF (pmiss(imt) > 0 .and. &
                (mt_forcing_log_count < 3 .or. .not. mt_forc_precip_warned_invalid(imt))) &
               emit_header = .true.
            IF (vdry(imt) > 0 .and. &
                (mt_forcing_log_count < 3 .or. .not. mt_forc_vapor_warned_dry(imt))) &
               emit_header = .true.
            IF (vmiss(imt) > 0 .and. &
                (mt_forcing_log_count < 3 .or. .not. mt_forc_vapor_warned_invalid(imt))) &
               emit_header = .true.
         ENDDO
         IF (emit_header) WRITE(*,'(/,A,I4.4,A,I3.3,A,I5.5)') &
            'Checking material forcing at ', idate(1), '-', idate(2), '-', idate(3)
         DO imt = 1, nmts
            IF (.not. mt_uses_generic_water(imt)) CYCLE
            IF (mt_forcing_log_count < 3 .and. pcnt(imt) > 0) THEN
               WRITE(*,'(3A,I0,A,F12.5,A,F12.5)') '  precip ', trim(materials(imt)%name), &
                  ' n=', pcnt(imt), ' min=', pmin(imt), ' max=', pmax(imt)
            ENDIF
            IF (mt_forcing_log_count < 3 .and. vcnt(imt) > 0) THEN
               WRITE(*,'(3A,I0,A,F12.5,A,F12.5)') '  vapor   ', trim(materials(imt)%name), &
                  ' n=', vcnt(imt), ' min=', vmin(imt), ' max=', vmax(imt)
            ENDIF
            IF (pdry(imt) > 0 .and. &
                (mt_forcing_log_count < 3 .or. .not. mt_forc_precip_warned_dry(imt))) THEN
               IF (mt_forc_phase_fill(imt)) THEN
                  WRITE(*,'(3A,I0,A)') '  NOTE precip ', trim(materials(imt)%name), &
                     ' near-dry patches=', pdry(imt), &
                     ' (rain: vapour equilibrium; snow: last reliable precip ratio)'
               ELSE
                  WRITE(*,'(3A,I0,A)') '  WARNING precip ', trim(materials(imt)%name), &
                     ' near-dry fallback patches=', pdry(imt), &
                     ' (retaining the previous valid/default precip ratio)'
               ENDIF
               mt_forc_precip_warned_dry(imt) = .true.
            ENDIF
            IF (pmiss(imt) > 0 .and. &
                (mt_forcing_log_count < 3 .or. .not. mt_forc_precip_warned_invalid(imt))) THEN
               WRITE(*,'(3A,I0,A)') '  WARNING precip ', trim(materials(imt)%name), &
                  ' invalid/missing wet fallback patches=', pmiss(imt), &
                  ' (retaining the previous valid/default precip ratio)'
               mt_forc_precip_warned_invalid(imt) = .true.
            ENDIF
            IF (vdry(imt) > 0 .and. &
                (mt_forcing_log_count < 3 .or. .not. mt_forc_vapor_warned_dry(imt))) THEN
               WRITE(*,'(3A,I0,A)') '  WARNING vapor ', trim(materials(imt)%name), &
                  ' near-dry fallback patches=', vdry(imt), &
                  ' (retaining the previous valid/default vapor ratio)'
               mt_forc_vapor_warned_dry(imt) = .true.
            ENDIF
            IF (vmiss(imt) > 0 .and. &
                (mt_forcing_log_count < 3 .or. .not. mt_forc_vapor_warned_invalid(imt))) THEN
               WRITE(*,'(3A,I0,A)') '  WARNING vapor ', trim(materials(imt)%name), &
                  ' invalid/missing wet fallback patches=', vmiss(imt), &
                  ' (retaining the previous valid/default vapor ratio)'
               mt_forc_vapor_warned_invalid(imt) = .true.
            ENDIF
         ENDDO
      ENDIF

      mt_forcing_log_count = mt_forcing_log_count + 1
      deallocate(pmin, pmax, vmin, vmax, pcnt, vcnt, pmiss, pdry, vmiss, vdry)
   END SUBROUTINE material_forcing_log_ranges

   real(r8) FUNCTION material_forcing_diag_value (imt, value)
      IMPLICIT NONE
      integer, intent(in) :: imt
      real(r8), intent(in) :: value

      IF (material_is_isotope(imt)) THEN
         material_forcing_diag_value = material_forcing_ratio_to_delta(value, materials(imt)%ref_ratio)
      ELSE
         material_forcing_diag_value = value
      ENDIF
   END FUNCTION material_forcing_diag_value

   logical FUNCTION material_forcing_pair_present (fprefix, vname)
      IMPLICIT NONE
      character(len=*), intent(in) :: fprefix, vname
      material_forcing_pair_present = material_forcing_token_present(fprefix) .and. material_forcing_token_present(vname)
   END FUNCTION material_forcing_pair_present

   logical FUNCTION material_forcing_token_present (token)
      IMPLICIT NONE
      character(len=*), intent(in) :: token
      character(len=256) :: low

      low = material_lower(trim(token))
      material_forcing_token_present = len_trim(low) > 0 .and. trim(low) /= 'null' .and. trim(low) /= 'none'
   END FUNCTION material_forcing_token_present

   integer FUNCTION material_forcing_parse_mode (token)
      IMPLICIT NONE
      character(len=*), intent(in) :: token
      character(len=256) :: low

      low = material_lower(trim(token))
      IF (trim(low) == 'delta') THEN
         material_forcing_parse_mode = MODE_DELTA
      ELSEIF (trim(low) == 'heavy_over_total' .or. trim(low) == 'ratio_to_total' .or. &
              trim(low) == 'over_total') THEN
         material_forcing_parse_mode = MODE_HEAVY_OVER_TOTAL
      ELSEIF (trim(low) == 'normalized_over_total' .or. trim(low) == 'normalized_ratio' .or. &
              trim(low) == 'isogsm' .or. trim(low) == 'standard_over_total') THEN
         material_forcing_parse_mode = MODE_NORMALIZED_OVER_TOTAL
      ELSEIF (trim(low) == 'direct' .or. trim(low) == 'value' .or. trim(low) == 'raw') THEN
         material_forcing_parse_mode = MODE_VALUE
      ELSE
         IF (p_is_master) WRITE(*,'(3A)') &
            'ERROR material_forcing_parse_mode: unknown forcing_input_mode "', trim(token), '".'
         CALL CoLM_stop()
         material_forcing_parse_mode = MODE_VALUE
      ENDIF
   END FUNCTION material_forcing_parse_mode

   SUBROUTINE material_forcing_csv_token (csv, idx, token)
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
   END SUBROUTINE material_forcing_csv_token

   integer FUNCTION material_forcing_csv_int (csv, idx, default_value)
      IMPLICIT NONE
      character(len=*), intent(in) :: csv
      integer, intent(in) :: idx, default_value
      character(len=256) :: token
      integer :: ierr

      material_forcing_csv_int = default_value
      CALL material_forcing_csv_token(csv, idx, token)
      IF (.not. material_forcing_token_present(token)) RETURN
      read(token, *, iostat=ierr) material_forcing_csv_int
      IF (ierr /= 0) material_forcing_csv_int = default_value
   END FUNCTION material_forcing_csv_int

END MODULE MOD_MT_Forcing
#endif
