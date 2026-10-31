#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_Hist

   USE MOD_Precision
   USE MOD_Vars_Global, only: spval
   USE MOD_NetCDFSerial
   USE MOD_HistGridded
#ifdef UNSTRUCTURED
   USE MOD_HistVector
#endif
#ifdef SinglePoint
   USE MOD_HistSingle
#endif
   USE MOD_MT_Defs, only: nmts, materials, mass_to_delta, mt_tiny, &
      mt_delta_sanity_max, mt_flux_water_min_for_delta, mt_water_min_for_delta, &
      mt_water_min_for_ratio, &
      material_uses_delta_diagnostics, mt_uses_generic_water, &
      material_concentration_units, material_is_nonvolatile_solute, material_is_isotope, &
      material_aquifer_actual_water, material_aquifer_actual_mass
   USE MOD_MT_Vars
   USE MOD_MT_LandWork, only: mt_land_work
   USE MOD_MT_Conservation, only: mt_apply_decay, mt_balance_check
   USE MOD_MT_Namelist, only: mt_any_enabled
   USE MOD_MT_HostFixes, only: mt_hist_continue_file
   USE MOD_MT_HostHistory, only: nac_frcsat

   IMPLICIT NONE
   PRIVATE

   ! Host phases use the material facade; provider writers are explicitly exported.
   PUBLIC :: mt_hist_accumulate, mt_land_finish_step
   PUBLIC :: mt_hist_out
   PUBLIC :: HistForm
   PUBLIC :: write_history_mt_ratio_2d
   PUBLIC :: write_history_variable_2d
   PUBLIC :: write_history_variable_3d

   PUBLIC :: mt_history_record, mt_hist_initialize, mt_hist_flush, mt_hist_begin
   PUBLIC :: mt_hist_write_record, mt_hist_end, mt_hist_snapshot, mt_hist_write_frcsat

   type mt_history_record
      private
      character(len=256) :: filename
      integer :: itime
   end type mt_history_record
   character(len=256) :: file_last_material = 'null'

   character(len=10) :: HistForm = 'Gridded'

CONTAINS

   SUBROUTINE mt_hist_initialize(lulcc_call)
      USE MOD_Vars_1DAccFluxes, only: FLUSH_acc_fluxes
      logical, optional, intent(in) :: lulcc_call
      logical :: flush_mt

      ! Startup may already hold restarted material sums; only LULCC clears them.
      flush_mt = .false.
      IF (present(lulcc_call)) flush_mt = lulcc_call
      CALL FLUSH_acc_fluxes(flush_mt=flush_mt)
      CALL flush_Material_Acc()
   END SUBROUTINE mt_hist_initialize

   SUBROUTINE mt_hist_flush()
      USE MOD_Vars_1DAccFluxes, only: FLUSH_acc_fluxes
      CALL FLUSH_acc_fluxes()
      CALL flush_Material_Acc()
   END SUBROUTINE mt_hist_flush

   SUBROUTINE mt_hist_snapshot()
#ifdef GridRiverLakeFlow
      USE MOD_MT_RiverLake, only: mt_riverlake_snapshot_hist_acc
      USE MOD_MT_Sediment, only: mt_sediment_snapshot_hist_acc
      CALL mt_riverlake_snapshot_hist_acc()
      CALL mt_sediment_snapshot_hist_acc()
#endif
   END SUBROUTINE mt_hist_snapshot

   SUBROUTINE mt_hist_begin(record, dir_hist, casename, cdate, idate, HistForm, deltim, nac, write_time)
      USE MOD_SPMD_Task
      USE MOD_Namelist, only: DEF_HIST_mode, DEF_HIST_WriteBack
#ifdef SinglePoint
      USE MOD_Namelist, only: USE_SITE_HistWriteBack
#endif
      USE MOD_TimeManager, only: minutes_since_1900
      type(mt_history_record), intent(out) :: record
      character(len=*), intent(in) :: dir_hist, casename, cdate, HistForm
      integer, intent(in) :: idate(3)
      real(r8), intent(in) :: deltim, nac
      real(r8) :: history_window_seconds
#ifdef SinglePoint
      logical :: site_hist_writeback_saved
#endif
      ! Use the host dispatch without importing MOD_Hist back into its dependency.
      interface
         SUBROUTINE write_time(filename, filelast, dataname, time, itime)
            character(len=*), intent(in) :: filename, filelast, dataname
            integer, intent(in) :: time(3)
            integer, intent(out) :: itime
         END SUBROUTINE write_time
      end interface

         ! material history goes to its own file, with the averaging window;
         ! nothing is written when no material is enabled
         IF (mt_any_enabled()) THEN
            history_window_seconds=0._r8
            IF (p_is_worker) history_window_seconds=max(nac,0._r8)*deltim
#ifdef USEMPI
            CALL mpi_allreduce(MPI_IN_PLACE, history_window_seconds, 1, MPI_REAL8, MPI_MAX, p_comm_glb, p_err)
#endif
            record%filename = trim(dir_hist) // '/' // trim(casename) //'_hist_material_'//trim(cdate)//'.nc'
            CALL mt_hist_continue_file (record%filename, idate, HistForm, file_last_material)
#ifdef SinglePoint
            ! The single-site cache belongs to the host stream, not this file.
            site_hist_writeback_saved = USE_SITE_HistWriteBack
            USE_SITE_HistWriteBack = .false.
#endif
            CALL write_time (record%filename, file_last_material, 'time', idate, record%itime)
#ifdef SinglePoint
            USE_SITE_HistWriteBack = site_hist_writeback_saved
#endif
            ! Written by the master into the single file; in 'block' mode that file
            ! does not exist (per-block files), and with write-back it is created
            ! by the write-back process, so the window is not written there.
            IF (HistForm=='Gridded' .and. p_is_master .and. trim(DEF_HIST_mode)=='one' &
                .and. .not. DEF_HIST_WriteBack) THEN
               CALL ncio_write_serial_time (record%filename, 'history_window_seconds', &
                  record%itime, history_window_seconds, 'time')
               CALL ncio_write_serial_time (record%filename, 'history_window_end_minutes', &
                  record%itime, real(minutes_since_1900(idate(1), idate(2), idate(3)), r8) + &
                  real(mod(idate(3),60),r8)/60._r8, 'time')
               IF (record%itime==1) THEN
                  CALL ncio_put_attr(record%filename, 'history_window_seconds', 'units', 's')
                  CALL ncio_put_attr(record%filename, 'history_window_seconds', 'long_name', &
                     'elapsed window ending at history_window_end_minutes; terminal and resumed records can overlap')
                  CALL ncio_put_attr(record%filename, 'history_window_end_minutes', 'units', &
                     'minutes since 1900-1-1 0:0:0')
               ENDIF
            ENDIF
         ENDIF
   END SUBROUTINE mt_hist_begin

   SUBROUTINE mt_hist_write_record(record, hist_form_in, sumarea, filter, maxsnl, nl_soil, &
      forcing_has_missing_value, forcmask_pch)
      USE MOD_DataType, only: block_data_real8_2d
      type(mt_history_record), intent(in) :: record
      character(len=*), intent(in) :: hist_form_in
      type(block_data_real8_2d), intent(inout) :: sumarea
      logical, intent(inout) :: filter(:)
      integer, intent(in) :: maxsnl, nl_soil
      logical, intent(in) :: forcing_has_missing_value
      logical, allocatable, intent(in) :: forcmask_pch(:)
      logical :: forcmask_empty(0)

         IF (mt_any_enabled()) THEN
            IF (allocated(forcmask_pch)) THEN
               CALL mt_hist_out (record%filename, record%itime, hist_form_in, &
                  sumarea, filter, maxsnl, nl_soil, forcing_has_missing_value, forcmask_pch)
            ELSE
               CALL mt_hist_out (record%filename, record%itime, hist_form_in, &
                  sumarea, filter, maxsnl, nl_soil, forcing_has_missing_value, forcmask_empty)
            ENDIF
         ENDIF
   END SUBROUTINE mt_hist_write_record

   SUBROUTINE mt_hist_end(record)
      type(mt_history_record), intent(in) :: record
      IF (mt_any_enabled()) file_last_material = record%filename
   END SUBROUTINE mt_hist_end

   SUBROUTINE mt_hist_write_frcsat(is_hist, a_frcsat, file_hist, itime_in_file, sumarea, filter, write_field)
      USE MOD_DataType, only: block_data_real8_2d
      logical, intent(in) :: is_hist
      real(r8), allocatable, intent(inout) :: a_frcsat(:)
      character(len=*), intent(in) :: file_hist
      integer, intent(in) :: itime_in_file
      type(block_data_real8_2d), intent(in) :: sumarea
      logical, intent(in) :: filter(:)
      interface
         SUBROUTINE write_field(is_hist, acc_vec, file_hist, varname, itime_in_file, sumarea, filter, &
            longname, units, acc_num, input_mode)
            USE MOD_Precision
            USE MOD_DataType, only: block_data_real8_2d
            logical, intent(in) :: is_hist
            real(r8), intent(inout) :: acc_vec(:)
            character(len=*), intent(in) :: file_hist, varname, longname, units
            integer, intent(in) :: itime_in_file
            type(block_data_real8_2d), intent(in) :: sumarea
            logical, intent(in) :: filter(:)
            real(r8), intent(in), optional :: acc_num(:)
            character(len=*), intent(in), optional :: input_mode
         END SUBROUTINE write_field
      end interface

      IF (mt_any_enabled()) THEN
         IF (allocated(nac_frcsat)) THEN
            WHERE (nac_frcsat <= 0._r8) a_frcsat = spval
         ENDIF
         CALL write_field(is_hist, a_frcsat, file_hist, 'f_frcsat', itime_in_file, sumarea, filter, &
            'fraction of saturation area on ground', '-', acc_num=nac_frcsat)
      ELSE
         CALL write_field(is_hist, a_frcsat, file_hist, 'f_frcsat', itime_in_file, sumarea, filter, &
            'fraction of saturation area on ground', '-')
      ENDIF
   END SUBROUTINE mt_hist_write_frcsat


   SUBROUTINE mt_land_finish_step (work, ipatch, snl, maxsnl, nl_soil, deltim, patchtype, &
      errorw, endwb, totwb, forc_prc, forc_prl, fevpa, rnof, ldew_rain, ldew_snow, &
      wliq_soisno, wice_soisno, wa, wdsrf, wetwat, scv &
#if (defined CaMa_Flood) || (defined GridRiverLakeFlow)
      , flood_enabled, flood_input, flood_evaporation &
#endif
      )
      type(mt_land_work), intent(inout) :: work
      integer, intent(in) :: ipatch, snl, maxsnl, nl_soil, patchtype
      real(r8), intent(in) :: deltim, errorw, endwb, totwb, forc_prc, forc_prl, fevpa, rnof
      real(r8), intent(in) :: ldew_rain, ldew_snow, wa, wdsrf, wetwat, scv
      real(r8), intent(in) :: wliq_soisno(snl+1:nl_soil), wice_soisno(snl+1:nl_soil)
#if (defined CaMa_Flood) || (defined GridRiverLakeFlow)
      logical, intent(in) :: flood_enabled
      real(r8), intent(in) :: flood_input, flood_evaporation
#endif

      ! Water terms for the material balance. The host endwb has the flood
      ! infiltration taken out and fevpa includes flood evaporation; undo
      ! both so dS, input and output describe the same land water store.
      work%flood_input = 0._r8
      work%evaporation = fevpa
#if (defined CaMa_Flood) || (defined GridRiverLakeFlow)
      IF (flood_enabled .and. patchtype == 0) work%flood_input = flood_input
      IF (patchtype == 0) work%evaporation = fevpa - flood_evaporation
#endif
      CALL mt_apply_decay(ipatch, snl, nl_soil, deltim)
      CALL mt_balance_check(ipatch, snl, nl_soil, deltim, work%balance_error, &
         patchtype_in = patchtype, water_err_in = errorw, &
         water_dS_in = endwb - totwb + work%flood_input * deltim, &
         water_input_in = (forc_prc + forc_prl + work%flood_input) * deltim, &
         water_output_in = (work%evaporation + rnof) * deltim, &
         water_evap_in = work%evaporation * deltim, &
         water_rnof_in = rnof * deltim &
#ifdef GridRiverLakeFlow
         ! flood water may carry a river signature: no fixed-signature check
         , flood_heterogeneous_in = flood_enabled .and. patchtype == 0 &
#endif
         )

      CALL mt_hist_accumulate(ipatch, snl, maxsnl, nl_soil, ldew_rain, ldew_snow, &
         wliq_soisno(snl+1:nl_soil), wice_soisno(snl+1:nl_soil), &
         wa, wdsrf, wetwat, scv)
      CALL work%release()
   END SUBROUTINE mt_land_finish_step

   SUBROUTINE mt_hist_accumulate (ipatch, snl, maxsnl, nl_soil, &
      ldew_rain, ldew_snow, wliq_soisno, wice_soisno, &
      wa, wdsrf, wetwat, scv)
      IMPLICIT NONE
      integer,  intent(in) :: ipatch, snl, maxsnl, nl_soil
      real(r8), intent(in) :: ldew_rain, ldew_snow
      real(r8), intent(in) :: wliq_soisno(snl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno(snl+1:nl_soil)
      ! Post-step pool states paired with the existing canopy/soil/snow args.
      ! These feed four additional per-pool material-ratio diagnostics in
      ! MOD_Hist (f_mt_conc_wa / wdsrf / wetwat / scv). Without them the
      ! aquifer, surface water, wetland and pre-layer thin-snow δ values are
      ! only observable in the restart file — which rarely matches the
      ! history output cadence.
      real(r8), intent(in) :: wa, wdsrf, wetwat, scv

      integer :: imt, j, jsnow
      real(r8) :: layer_water, layer_material, actual_aquifer_water, actual_aquifer_mass

      IF (nmts <= 0) RETURN

      ! Canopy water
      a_water_ldew(ipatch) = a_water_ldew(ipatch) + (ldew_rain + ldew_snow)
      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         a_mt_ldew_mass(imt, ipatch) = a_mt_ldew_mass(imt, ipatch) &
            + mt_ldew_rain(imt, ipatch) + mt_ldew_snow(imt, ipatch)
      ENDDO

      ! Soil layers (1:nl_soil) — fixed dimension, safe for MPI gather
      DO j = 1, nl_soil
         layer_water = wliq_soisno(j) + wice_soisno(j)
         a_water_soil(j, ipatch) = a_water_soil(j, ipatch) + layer_water
         DO imt = 1, nmts
            IF (.not. mt_uses_generic_water(imt)) CYCLE
            layer_material = mt_wliq_soisno(imt, j, ipatch) + mt_wice_soisno(imt, j, ipatch)
            IF (material_is_nonvolatile_solute(imt) .and. &
                layer_water <= mt_water_min_for_ratio) THEN
               a_mt_layer_dry_mass(imt, ipatch) = a_mt_layer_dry_mass(imt, ipatch) + layer_material
            ELSE
               a_mt_soil_mass(imt, j, ipatch) = a_mt_soil_mass(imt, j, ipatch) + layer_material
            ENDIF
         ENDDO
      ENDDO

      ! Snow layers — separate accumulator with fixed dimension abs(maxsnl).
      ! Snow layer index j ranges from snl+1 to 0 (negative indices).
      ! Map to positive index: jsnow = j - maxsnl  (1-based, from top to bottom)
      ! e.g., maxsnl=-5: j=-4 → jsnow=1, j=-3 → jsnow=2, ..., j=0 → jsnow=5
      IF (snl < 0) THEN
         DO j = snl + 1, 0
            jsnow = j - maxsnl  ! positive index (1 to abs(maxsnl))
            a_water_snow(jsnow, ipatch) = a_water_snow(jsnow, ipatch) &
               + wliq_soisno(j) + wice_soisno(j)
            DO imt = 1, nmts
               IF (.not. mt_uses_generic_water(imt)) CYCLE
               layer_water = wliq_soisno(j) + wice_soisno(j)
               layer_material = mt_wliq_soisno(imt, j, ipatch) + mt_wice_soisno(imt, j, ipatch)
               IF (material_is_nonvolatile_solute(imt) .and. &
                   layer_water <= mt_water_min_for_ratio) THEN
                  a_mt_layer_dry_mass(imt, ipatch) = &
                     a_mt_layer_dry_mass(imt, ipatch) + layer_material
               ELSE
                  a_mt_snow_mass(imt, jsnow, ipatch) = &
                     a_mt_snow_mass(imt, jsnow, ipatch) + layer_material
               ENDIF
            ENDDO
         ENDDO
      ENDIF

            ! Aquifer: positive storage and signed debt are different diagnostics.
            ! wa<=0 is a wetland/aquifer debt state, not a physical concentration.
            ! Keep mt_wa state synchronized with wa, but filter sub-mm debt
            ! from concentration diagnostics to avoid tiny-denominator deltas.
            IF (wa > 1._r8) a_water_wa(ipatch) = a_water_wa(ipatch) + wa
            IF (wa < -1._r8) a_water_wa_debt(ipatch) = a_water_wa_debt(ipatch) - wa
            actual_aquifer_water = material_aquifer_actual_water(wa, mt_aquifer_ref_water(ipatch))
            IF (mt_aquifer_ref_water(ipatch) > 0._r8 .and. &
                actual_aquifer_water > mt_water_min_for_ratio) &
               a_water_aquifer_actual(ipatch) = a_water_aquifer_actual(ipatch) + actual_aquifer_water
            IF (mt_aquifer_ref_water(ipatch) <= 0._r8 .and. wa > 1._r8) &
               a_water_aquifer_actual(ipatch) = a_water_aquifer_actual(ipatch) + wa
         a_water_wdsrf (ipatch) = a_water_wdsrf (ipatch) + wdsrf
         a_water_wetwat(ipatch) = a_water_wetwat(ipatch) + wetwat
         ! CoLM's `scv` is total snow water equivalent even after layered
         ! snow exists, but `mt_scv` only mirrors the pre-layer thin-snow
         ! pool. Once snl<0, material mass lives in mt_wliq/mt_wice snow
         ! layers and is diagnosed via f_mt_conc_soisno_*; do not pair
         ! that total scv water with the thin-snow material denominator.
         IF (snl == 0 .and. scv > mt_tiny) THEN
            a_water_scv(ipatch) = a_water_scv(ipatch) + scv
         ENDIF
         DO imt = 1, nmts
            IF (.not. mt_uses_generic_water(imt)) CYCLE
            IF (material_is_isotope(imt) .and. mt_aquifer_ref_water(ipatch) > 0._r8 .and. &
                actual_aquifer_water > mt_water_min_for_ratio) THEN
               actual_aquifer_mass = material_aquifer_actual_mass(mt_wa(imt, ipatch), &
                  mt_aquifer_ref_mass(imt, ipatch))
               a_mt_aquifer_actual_mass(imt, ipatch) = &
                  a_mt_aquifer_actual_mass(imt, ipatch) + actual_aquifer_mass
               IF (wa < -1._r8) a_mt_wa_debt_mass(imt, ipatch) = &
                  a_mt_wa_debt_mass(imt, ipatch) - wa * actual_aquifer_mass / actual_aquifer_water
            ELSE
            IF (material_is_isotope(imt) .and. wa > 1._r8) &
               a_mt_aquifer_actual_mass(imt, ipatch) = &
                  a_mt_aquifer_actual_mass(imt, ipatch) + mt_wa(imt, ipatch)
            IF (wa > 1._r8) THEN
               a_mt_wa_mass(imt, ipatch) = a_mt_wa_mass(imt, ipatch) + mt_wa(imt, ipatch)
            ENDIF
               IF (wa < -1._r8) THEN
                  a_mt_wa_debt_mass(imt, ipatch) = a_mt_wa_debt_mass(imt, ipatch) &
                     + max(-mt_wa(imt, ipatch), 0._r8)
               ENDIF
            ENDIF
            a_mt_wdsrf_mass (imt, ipatch) = a_mt_wdsrf_mass (imt, ipatch) + mt_wdsrf (imt, ipatch)
         a_mt_wetwat_mass(imt, ipatch) = a_mt_wetwat_mass(imt, ipatch) + mt_wetwat(imt, ipatch)
         a_mt_surface_residue_mass(imt, ipatch) = a_mt_surface_residue_mass(imt, ipatch) + &
            mt_surface_residue(imt, ipatch)
         a_mt_subsurface_residue_mass(imt, ipatch) = a_mt_subsurface_residue_mass(imt, ipatch) + &
            mt_subsurface_residue(imt, ipatch)
         a_mt_solid_mass(imt, ipatch) = a_mt_solid_mass(imt, ipatch) + &
            mt_canopy_solid(imt, ipatch) + mt_surface_solid(imt, ipatch) + &
            mt_subsurface_solid(imt, ipatch) + mt_waterstorage_solid(imt, ipatch)
         DO j = max(snl + 1, lbound(mt_solid_soisno, 2)), nl_soil
            a_mt_solid_mass(imt, ipatch) = a_mt_solid_mass(imt, ipatch) + &
               mt_solid_soisno(imt, j, ipatch)
         ENDDO
         IF (snl == 0 .and. scv > mt_tiny) THEN
            a_mt_scv_mass(imt, ipatch) = a_mt_scv_mass(imt, ipatch) + mt_scv(imt, ipatch)
         ENDIF
      ENDDO
   END SUBROUTINE mt_hist_accumulate

   SUBROUTINE mt_hist_out (file_hist, itime_in_file, hist_form_in, &
      sumarea, filter, maxsnl, nl_soil, forcing_has_missing_value, forcmask_pch)
   USE MOD_Block
   USE MOD_DataType
   USE MOD_LandPatch, only: numpatch
   USE MOD_SPMD_Task, only: p_is_worker
   USE MOD_Vars_TimeInvariants, only: patchtype, patchmask
   USE MOD_Namelist, only: DEF_hist_vars
   USE MOD_MT_Lifecycle, only: material_lifecycle_land_history

   IMPLICIT NONE

   character(len=*), intent(in) :: file_hist
   integer, intent(in) :: itime_in_file
   character(len=*), intent(in) :: hist_form_in
   type(block_data_real8_2d), intent(inout) :: sumarea
   logical, intent(inout) :: filter(:)
   integer, intent(in) :: maxsnl, nl_soil
   logical, intent(in) :: forcing_has_missing_value
   logical, intent(in) :: forcmask_pch(:)
            real(r8), allocatable :: mt_mass_soisno (:,:)
            real(r8), allocatable :: water_soisno    (:,:)
            real(r8), allocatable :: mt_delta_flux  (:)
            real(r8), allocatable :: snowpack_mass_pat(:)
            real(r8), allocatable :: snowpack_water_pat(:)
            real(r8), allocatable :: ones_arr        (:)
            ! Per-patch (pseudo_mass, pseudo_water) pair for the
            ! leaf-water-weighted bulk delta. pseudo_water = leaf
            ! water moles, pseudo_mass = R_b * moles. Feeding these
            ! to write_history_mt_delta_2d makes the gridded
            ! result mass_to_delta(sum(area * R_b * moles),
            ! sum(area * moles)) -- a true leaf-water-weighted mean
            ! that obeys the sum-then-divide invariant. Reusing
            ! mt_delta_flux for this would mix two different units
            ! across the same buffer.
            real(r8), allocatable :: leaf_mass_pat   (:)
            real(r8), allocatable :: leaf_water_pat  (:)
            character(len=64)  :: mt_varname
            character(len=128) :: mt_longname
            character(len=32)  :: mt_ratio_word, mt_ratio_units
            real(r8) :: delta_loc, leaf_R_loc
            integer :: imt_loc, ip_loc, j_loc, jsnow_loc
#ifdef SinglePoint
            logical :: site_hist_writeback_saved
#endif

      HistForm = hist_form_in

#ifdef SinglePoint
      ! Match the material time writer: bypass the host's shared site cache,
      ! including provider output, then restore the host write-back policy.
      site_hist_writeback_saved = USE_SITE_HistWriteBack
      USE_SITE_HistWriteBack = .false.
#endif

      IF ((p_is_worker) .and. (numpatch > 0)) THEN
         filter(:) = patchtype < 99
         IF (forcing_has_missing_value) filter = filter .and. forcmask_pch
         filter = filter .and. patchmask
      ENDIF
      IF (HistForm == 'Gridded') THEN
         CALL mp2g_hist%get_sumarea (sumarea, filter)
      ENDIF

            IF (nmts > 0) THEN
               IF (p_is_worker) THEN
                  allocate (mt_mass_soisno (maxsnl+1:nl_soil, numpatch))
                  allocate (water_soisno    (maxsnl+1:nl_soil, numpatch))
                  allocate (mt_delta_flux  (numpatch))
                  allocate (snowpack_mass_pat (numpatch))
                  allocate (snowpack_water_pat(numpatch))
                  allocate (ones_arr        (numpatch));  ones_arr = 1._r8
                  allocate (leaf_mass_pat   (numpatch))
                  allocate (leaf_water_pat  (numpatch))
               ELSE
                  allocate (mt_mass_soisno (0,0))
                  allocate (water_soisno    (0,0))
                  allocate (mt_delta_flux  (0))
                  allocate (snowpack_mass_pat (0))
                  allocate (snowpack_water_pat(0))
                  allocate (ones_arr        (0))
                  allocate (leaf_mass_pat   (0))
                  allocate (leaf_water_pat  (0))
               ENDIF

               DO imt_loc = 1, nmts
                  IF (p_is_worker) THEN
                     IF (numpatch > 0) THEN
                        filter(:) = patchtype < 99
                        IF (forcing_has_missing_value) filter = filter .and. forcmask_pch
                        filter = filter .and. patchmask
                     ENDIF
                  ENDIF
                  IF (.not. mt_uses_generic_water(imt_loc)) CYCLE
                  ! The pool diagnostics below are named f_mt_conc_* for every
                  ! material, but for isotopes they hold a heavy/total ratio, not
                  ! a concentration -- unlike the flux diagnostics, which split
                  ! into f_mt_delta_* and f_mt_conc_* by category. The names
                  ! are kept because renaming them would silently drop these
                  ! fields from every existing DEF_hist_vars list and break
                  ! analysis scripts, so long_name has to carry the distinction:
                  ! say what the ratio is OF, not just that it is one.
                  IF (material_uses_delta_diagnostics(imt_loc)) THEN
                     mt_ratio_word  = 'heavy/total ratio'
                     mt_ratio_units = 'R'
                  ELSE
                     mt_ratio_word  = 'concentration'
                     mt_ratio_units = material_concentration_units(imt_loc)
                  ENDIF

                  IF (material_uses_delta_diagnostics(imt_loc)) THEN
                     write(mt_varname , '(A,A)') 'f_mt_delta_precip_', trim(materials(imt_loc)%name)
                     ! a_water_precip also books dew/frost deposition (MOD_MT_Evapo,
                     ! MOD_MT_SoilWater), so say so in the name.
                     write(mt_longname, '(A,A,A)') 'precipitation plus dew/frost deposition material delta (', &
                        trim(materials(imt_loc)%name), ')'
                     CALL write_history_mt_delta_2d (DEF_hist_vars%xy_prc .or. DEF_hist_vars%xy_prl, &
                        a_mt_precip(imt_loc, :), a_water_precip(imt_loc, :), &
                        materials(imt_loc)%ref_ratio, file_hist, trim(mt_varname), &
                        itime_in_file, filter, trim(mt_longname), 'permil', &
                        water_min_override = mt_flux_water_min_for_delta)

                     write(mt_varname , '(A,A)') 'f_mt_delta_runoff_', trim(materials(imt_loc)%name)
                     write(mt_longname, '(A,A,A)') 'total runoff material delta (', &
                        trim(materials(imt_loc)%name), ')'
                     CALL write_history_mt_delta_2d (DEF_hist_vars%rnof, &
                        a_mt_rnof(imt_loc, :), a_water_rnof(imt_loc, :), &
                        materials(imt_loc)%ref_ratio, file_hist, trim(mt_varname), &
                        itime_in_file, filter, trim(mt_longname), 'permil', &
                        water_min_override = mt_flux_water_min_for_delta)

                     ! --- Evapotranspiration flux delta ---
                     ! Pair a_mt_evap with gross evaporative water so the
                     ! denominator never goes near zero from dew/frost
                     ! cancellation (a_fevpa is the net flux; not valid as a
                     ! delta denominator). write_history_mt_delta_2d
                     ! maps mass and water separately to the grid and divides
                     ! once at the end -- preserving the
                     ! mass_to_delta(sum(mass), sum(water)) invariant
                     ! documented in MOD_MT_Defs.F90.
                     write(mt_varname , '(A,A)')   'f_mt_delta_evap_', trim(materials(imt_loc)%name)
                     write(mt_longname, '(A,A,A)') 'evapotranspiration material delta (', &
                        trim(materials(imt_loc)%name), ')'
                        CALL write_history_mt_delta_2d (DEF_hist_vars%fevpa, &
                           a_mt_evap(imt_loc, :), a_water_evap_gross(imt_loc, :), &
                           materials(imt_loc)%ref_ratio, file_hist, trim(mt_varname), &
                           itime_in_file, filter, trim(mt_longname), 'permil', &
                           water_min_override = mt_flux_water_min_for_delta)

                     write(mt_varname , '(A,A)')   'f_mt_evap_mass_', trim(materials(imt_loc)%name)
                     write(mt_longname, '(A,A,A)') 'signed evapotranspiration material mass (', &
                        trim(materials(imt_loc)%name), ')'
                     CALL write_history_variable_2d (DEF_hist_vars%fevpa, &
                        a_mt_evap(imt_loc, :), file_hist, trim(mt_varname), &
                        itime_in_file, sumarea, filter, trim(mt_longname), 'material amount/m2')

                     write(mt_varname , '(A,A)')   'f_mt_delta_soilevap_', trim(materials(imt_loc)%name)
                     write(mt_longname, '(A,A,A)') 'soil/surface evaporation material delta (', &
                        trim(materials(imt_loc)%name), ')'
                        CALL write_history_mt_delta_2d (DEF_hist_vars%fevpa, &
                           a_mt_soilevap(imt_loc, :), a_water_soilevap(imt_loc, :), &
                           materials(imt_loc)%ref_ratio, file_hist, trim(mt_varname), &
                           itime_in_file, filter, trim(mt_longname), 'permil', &
                           water_min_override = mt_flux_water_min_for_delta)

                     write(mt_varname , '(A,A)')   'f_mt_delta_canopyevap_', trim(materials(imt_loc)%name)
                     write(mt_longname, '(A,A,A)') 'canopy evaporation material delta (', &
                        trim(materials(imt_loc)%name), ')'
                        CALL write_history_mt_delta_2d (DEF_hist_vars%fevpa, &
                           a_mt_canopyevap(imt_loc, :), a_water_canopyevap(imt_loc, :), &
                           materials(imt_loc)%ref_ratio, file_hist, trim(mt_varname), &
                           itime_in_file, filter, trim(mt_longname), 'permil', &
                           water_min_override = mt_flux_water_min_for_delta)

                     write(mt_varname , '(A,A)')   'f_mt_delta_subl_', trim(materials(imt_loc)%name)
                     write(mt_longname, '(A,A,A)') 'sublimation material delta (', &
                        trim(materials(imt_loc)%name), ')'
                        CALL write_history_mt_delta_2d (DEF_hist_vars%fevpa, &
                           a_mt_subl(imt_loc, :), a_water_subl(imt_loc, :), &
                           materials(imt_loc)%ref_ratio, file_hist, trim(mt_varname), &
                           itime_in_file, filter, trim(mt_longname), 'permil', &
                           water_min_override = mt_flux_water_min_for_delta)

                     write(mt_varname , '(A,A)')   'f_mt_delta_wetland_evap_', trim(materials(imt_loc)%name)
                     write(mt_longname, '(A,A,A)') 'wetland evaporation/sublimation material delta (', &
                        trim(materials(imt_loc)%name), ')'
                        CALL write_history_mt_delta_2d (DEF_hist_vars%fevpa, &
                           a_mt_wetland_evap(imt_loc, :), a_water_wetland_evap(imt_loc, :), &
                           materials(imt_loc)%ref_ratio, file_hist, trim(mt_varname), &
                           itime_in_file, filter, trim(mt_longname), 'permil', &
                           water_min_override = mt_flux_water_min_for_delta)

                     ! --- Transpiration-only flux delta ---
                     ! write_history_mt_delta_2d does sum-then-divide on
                     ! the grid; no per-patch delta needed here.
                     write(mt_varname , '(A,A)')   'f_mt_delta_transp_', trim(materials(imt_loc)%name)
                     write(mt_longname, '(A,A,A)') 'transpiration material delta (', &
                        trim(materials(imt_loc)%name), ')'
                        CALL write_history_mt_delta_2d (DEF_hist_vars%etr, &
                           a_mt_transp(imt_loc, :), a_water_transp(imt_loc, :), &
                           materials(imt_loc)%ref_ratio, file_hist, trim(mt_varname), &
                           itime_in_file, filter, trim(mt_longname), 'permil', &
                           water_min_override = mt_flux_water_min_for_delta)

                     ! --- Transpiration source/xylem delta before NSS storage exchange ---
                     write(mt_varname , '(A,A)')   'f_mt_delta_transp_src_', trim(materials(imt_loc)%name)
                     write(mt_longname, '(A,A,A)') 'transpiration source/xylem material delta (', &
                        trim(materials(imt_loc)%name), ')'
                        CALL write_history_mt_delta_2d (DEF_hist_vars%etr, &
                           a_mt_transp_src(imt_loc, :), a_water_transp(imt_loc, :), &
                           materials(imt_loc)%ref_ratio, file_hist, trim(mt_varname), &
                           itime_in_file, filter, trim(mt_longname), 'permil', &
                           water_min_override = mt_flux_water_min_for_delta)

                     ! --- Leaf NSS state deltas ---
                  ! e-site: NSS evaporation-site state -- there is no
                  ! conserved per-patch water mass paired with this
                  ! delta (the Peclet number is a mixing factor, not
                  ! a mass), so area-weighted aggregation is the only
                  ! defensible choice. Longname documents that this
                  ! is a state diagnostic, not a mass-weighted bulk
                  ! pool delta.
                  IF (p_is_worker) THEN
                     DO ip_loc = 1, numpatch
                        delta_loc = mt_leaf_delta_e(imt_loc, ip_loc)
                        IF (delta_loc /= spval .and. &
                            abs(delta_loc) <= mt_delta_sanity_max) THEN
                           mt_delta_flux(ip_loc) = delta_loc
                        ELSE
                           mt_delta_flux(ip_loc) = spval
                        ENDIF
                     ENDDO
                  ENDIF
                  write(mt_varname , '(A,A)')   'f_mt_leaf_delta_e_', trim(materials(imt_loc)%name)
                  write(mt_longname, '(A,A,A)') 'leaf evaporation-site NSS delta, area-weighted state (', &
                     trim(materials(imt_loc)%name), ')'
                  CALL write_history_variable_2d (DEF_hist_vars%etr, &
                     mt_delta_flux, file_hist, trim(mt_varname), itime_in_file, &
                     sumarea, filter, trim(mt_longname), 'permil', acc_num = ones_arr)

                  ! b-site: bulk leaf water has a matching per-patch
                  ! mass scale (mt_leaf_water_moles, mol/m^2). Build
                  ! pseudo (mass, water) = (R_b * moles, moles) and
                  ! route through write_history_mt_delta_2d so the
                  ! grid result is mass_to_delta(sum(area*R_b*moles),
                  ! sum(area*moles)) -- a leaf-water-weighted bulk
                  ! delta consistent with the sum-then-divide
                  ! invariant. Patches with delta_b == spval, out-of-
                  ! sanity, or zero leaf water are spval'd in mass
                  ! and 0'd in water so the downstream filter at
                  ! water > mt_water_min_for_delta drops them.
                  IF (p_is_worker) THEN
                     DO ip_loc = 1, numpatch
                        delta_loc = mt_leaf_delta_b(imt_loc, ip_loc)
                        IF (delta_loc /= spval .and. &
                            abs(delta_loc) <= mt_delta_sanity_max .and. &
                            mt_leaf_water_moles(imt_loc, ip_loc) > 0._r8) THEN
                           leaf_R_loc = (1._r8 + delta_loc * 1.0e-3_r8) * &
                                        materials(imt_loc)%ref_ratio
                           leaf_mass_pat (ip_loc) = leaf_R_loc * &
                                                    mt_leaf_water_moles(imt_loc, ip_loc)
                           leaf_water_pat(ip_loc) = mt_leaf_water_moles(imt_loc, ip_loc)
                        ELSE
                           leaf_mass_pat (ip_loc) = spval
                           leaf_water_pat(ip_loc) = 0._r8
                        ENDIF
                     ENDDO
                  ENDIF
                  write(mt_varname , '(A,A)')   'f_mt_leaf_delta_b_', trim(materials(imt_loc)%name)
                  write(mt_longname, '(A,A,A)') 'bulk leaf-water NSS delta, leaf-water-mole weighted (', &
                     trim(materials(imt_loc)%name), ')'
                  CALL write_history_mt_delta_2d (DEF_hist_vars%etr, &
                     leaf_mass_pat, leaf_water_pat, &
                     materials(imt_loc)%ref_ratio, file_hist, trim(mt_varname), &
                     itime_in_file, filter, trim(mt_longname), 'permil')
                  ENDIF

                  ! Species-owned materials (for example CH4) were filtered by
                  ! mt_uses_generic_water above and write history
                  ! through callbacks. Solutes, decaying or not, keep these
                  ! water-pool diagnostics.
                  IF (.not. material_uses_delta_diagnostics(imt_loc)) THEN
                     write(mt_varname , '(A,A)') 'f_mt_conc_precip_', trim(materials(imt_loc)%name)
                     ! a_water_precip also books dew/frost deposition, which carries no
                     ! nonvolatile solute: dew-only periods read 0, not missing.
                     IF (material_is_nonvolatile_solute(imt_loc)) THEN
                        write(mt_longname, '(A,A,A)') &
                           'precipitation plus dew/frost deposition material concentration, dew/frost adding water only (', &
                           trim(materials(imt_loc)%name), ')'
                     ELSE
                        write(mt_longname, '(A,A,A)') &
                           'precipitation plus dew/frost deposition material concentration (', &
                           trim(materials(imt_loc)%name), ')'
                     ENDIF
                     CALL write_history_mt_ratio_2d (DEF_hist_vars%xy_prc .or. DEF_hist_vars%xy_prl, &
                        a_mt_precip(imt_loc, :), a_water_precip(imt_loc, :), &
                        file_hist, trim(mt_varname), itime_in_file, filter, &
                        trim(mt_longname), trim(mt_ratio_units))

                     write(mt_varname , '(A,A)') 'f_mt_conc_runoff_', trim(materials(imt_loc)%name)
                     write(mt_longname, '(A,A,A)') 'total runoff material concentration (', &
                        trim(materials(imt_loc)%name), ')'
                     CALL write_history_mt_ratio_2d (DEF_hist_vars%rnof, &
                        a_mt_rnof(imt_loc, :), a_water_rnof(imt_loc, :), &
                        file_hist, trim(mt_varname), itime_in_file, filter, &
                        trim(mt_longname), trim(mt_ratio_units))
                  ENDIF

                  ! Two-way equilibrium exchange with atmospheric vapour.
                  ! Written as a MASS, not a delta: it carries no water
                  ! flux by construction (zero net water, non-zero net
                  ! material), so there is nothing to normalise a delta by.
                  ! Signed -- positive means net uptake by the surface.
                  IF (material_uses_delta_diagnostics(imt_loc) .and. &
                      allocated(a_mt_vapor_exchange)) THEN
                     write(mt_varname , '(A,A)') 'f_mt_vapor_exchange_', &
                        trim(materials(imt_loc)%name)
                     write(mt_longname, '(A,A,A)') &
                        'zero-water-flux equilibrium vapour exchange, signed (', &
                        trim(materials(imt_loc)%name), ')'
                     CALL write_history_variable_2d (DEF_hist_vars%fevpa, &
                        a_mt_vapor_exchange(imt_loc, :), file_hist, trim(mt_varname), &
                        itime_in_file, sumarea, filter, trim(mt_longname), 'material amount/m2')
                  ENDIF

                  ! Waterless residue is an areal material inventory, not a
                  ! concentration. Tie its output to the existing surface-water
                  ! history switch until history namelist gains a dedicated key.
                  IF (material_is_nonvolatile_solute(imt_loc)) THEN
                     write(mt_varname , '(A,A)') 'f_mt_surface_residue_', trim(materials(imt_loc)%name)
                     write(mt_longname, '(A,A,A)') 'immobile surface material residue (', &
                        trim(materials(imt_loc)%name), ')'
                     CALL write_history_variable_2d (DEF_hist_vars%wdsrf, &
                        a_mt_surface_residue_mass(imt_loc, :), file_hist, trim(mt_varname), &
                        itime_in_file, sumarea, filter, trim(mt_longname), 'material amount/m2')

                     write(mt_varname , '(A,A)') 'f_mt_subsurface_residue_', trim(materials(imt_loc)%name)
                     write(mt_longname, '(A,A,A)') 'immobile subsurface material residue (', &
                        trim(materials(imt_loc)%name), ')'
                     CALL write_history_variable_2d (DEF_hist_vars%wa, &
                        a_mt_subsurface_residue_mass(imt_loc, :), file_hist, trim(mt_varname), &
                        itime_in_file, sumarea, filter, trim(mt_longname), 'material amount/m2')

                     write(mt_varname , '(A,A)') 'f_mt_layer_dry_inventory_', trim(materials(imt_loc)%name)
                     write(mt_longname, '(A,A,A)') 'dry snow/soil layer material inventory (', &
                        trim(materials(imt_loc)%name), ')'
                     CALL write_history_variable_2d (DEF_hist_vars%wliq_soisno, &
                        a_mt_layer_dry_mass(imt_loc, :), file_hist, trim(mt_varname), &
                        itime_in_file, sumarea, filter, trim(mt_longname), 'material amount/m2')

                     write(mt_varname , '(A,A)') 'f_mt_solid_inventory_', trim(materials(imt_loc)%name)
                     write(mt_longname, '(A,A,A)') 'precipitated solid material inventory (', &
                        trim(materials(imt_loc)%name), ')'
                     CALL write_history_variable_2d (DEF_hist_vars%wliq_soisno, &
                        a_mt_solid_mass(imt_loc, :), file_hist, trim(mt_varname), &
                        itime_in_file, sumarea, filter, trim(mt_longname), 'material amount/m2')
                  ENDIF

                  ! --- Canopy ratio ---
                  write(mt_varname , '(A,A)')   'f_mt_conc_ldew_', trim(materials(imt_loc)%name)
                  write(mt_longname, '(5A)') 'canopy material ', trim(mt_ratio_word), &
                     ' (', trim(materials(imt_loc)%name), ')'
                  CALL write_history_mt_ratio_2d (DEF_hist_vars%ldew, &
                     a_mt_ldew_mass(imt_loc, :), a_water_ldew, &
                     file_hist, trim(mt_varname), itime_in_file, &
                     filter, trim(mt_longname), trim(mt_ratio_units))

                  ! --- Soil + snow combined ratio ---
                  ! Snow accumulator is keyed by jsnow = j - maxsnl (positive);
                  ! soil accumulator is keyed by j (1..nl_soil). Map material
                  ! mass and water separately, then divide on the history grid.
                  ! This avoids area-fraction dilution when a grid cell has
                  ! multiple patch fractions with different material signatures.
                  IF (p_is_worker) THEN
                     mt_mass_soisno = spval
                     water_soisno    = spval
                     DO ip_loc = 1, numpatch
                        DO j_loc = maxsnl+1, 0
                           jsnow_loc = j_loc - maxsnl
                           IF (a_water_snow(jsnow_loc, ip_loc) > mt_tiny) THEN
                              mt_mass_soisno(j_loc, ip_loc) = &
                                 a_mt_snow_mass(imt_loc, jsnow_loc, ip_loc)
                              water_soisno(j_loc, ip_loc) = a_water_snow(jsnow_loc, ip_loc)
                           ENDIF
                        ENDDO
                        DO j_loc = 1, nl_soil
                           IF (a_water_soil(j_loc, ip_loc) > mt_tiny) THEN
                              mt_mass_soisno(j_loc, ip_loc) = &
                                 a_mt_soil_mass(imt_loc, j_loc, ip_loc)
                              water_soisno(j_loc, ip_loc) = a_water_soil(j_loc, ip_loc)
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDIF
                  write(mt_varname , '(A,A)')   'f_mt_conc_soisno_', trim(materials(imt_loc)%name)
                  write(mt_longname, '(5A)') 'soil/snow layer material ', trim(mt_ratio_word), &
                     ' (', trim(materials(imt_loc)%name), ')'
                  CALL write_history_mt_ratio_3d (DEF_hist_vars%wliq_soisno, &
                     mt_mass_soisno, water_soisno, file_hist, trim(mt_varname), itime_in_file, &
                     'soilsnow', maxsnl+1, nl_soil-maxsnl, filter, &
                     trim(mt_longname), trim(mt_ratio_units))

                  ! --- Total snowpack ratio/delta ---
                  ! CoLM's water-side scv is total SWE, while material_scv is
                  ! only the pre-layer thin-snow pool. Build a formal snowpack
                  ! diagnostic from the same pieces that own material mass:
                  ! layered snow accumulators plus thin-snow scv accumulators.
                  ! This keeps f_mt_conc_scv_* as a thin-snow diagnostic and
                  ! prevents the old thin-material / total-SWE denominator mix.
                  IF (p_is_worker) THEN
                     snowpack_mass_pat(:)  = a_mt_scv_mass(imt_loc, :)
                     snowpack_water_pat(:) = a_water_scv(:)
                     DO ip_loc = 1, numpatch
                        DO jsnow_loc = 1, abs(maxsnl)
                           snowpack_mass_pat(ip_loc) = snowpack_mass_pat(ip_loc) &
                              + a_mt_snow_mass(imt_loc, jsnow_loc, ip_loc)
                           snowpack_water_pat(ip_loc) = snowpack_water_pat(ip_loc) &
                              + a_water_snow(jsnow_loc, ip_loc)
                        ENDDO
                     ENDDO
                  ENDIF

                  write(mt_varname , '(A,A)')   'f_mt_conc_snowpack_', trim(materials(imt_loc)%name)
                  write(mt_longname, '(5A)') 'total snowpack material ', trim(mt_ratio_word), &
                     ' (', trim(materials(imt_loc)%name), ')'
                  CALL write_history_mt_ratio_2d (DEF_hist_vars%scv, &
                     snowpack_mass_pat, snowpack_water_pat, &
                     file_hist, trim(mt_varname), itime_in_file, &
                     filter, trim(mt_longname), trim(mt_ratio_units))

                  IF (material_uses_delta_diagnostics(imt_loc)) THEN
                     write(mt_varname , '(A,A)')   'f_mt_delta_snowpack_', trim(materials(imt_loc)%name)
                     write(mt_longname, '(A,A,A)') 'total snowpack material delta (', &
                        trim(materials(imt_loc)%name), ')'
                     CALL write_history_mt_delta_2d (DEF_hist_vars%scv, &
                        snowpack_mass_pat, snowpack_water_pat, &
                        materials(imt_loc)%ref_ratio, file_hist, trim(mt_varname), &
                        itime_in_file, filter, trim(mt_longname), 'permil')
                  ENDIF

                  ! --- Aquifer (wa) ratio ---
                  ! wa can be negative during wetland aquifer-debt episodes;
                  ! the mt_wetland redistribution keeps signed arithmetic.
                  write(mt_varname , '(A,A)')   'f_mt_conc_wa_', trim(materials(imt_loc)%name)
                  write(mt_longname, '(5A)') 'aquifer material ', trim(mt_ratio_word), &
                     ' (', trim(materials(imt_loc)%name), ')'
                     IF (material_is_isotope(imt_loc)) THEN
                        CALL write_history_mt_ratio_2d (DEF_hist_vars%wa, &
                           a_mt_aquifer_actual_mass(imt_loc, :), a_water_aquifer_actual, &
                           file_hist, trim(mt_varname), itime_in_file, &
                           filter, trim(mt_longname), trim(mt_ratio_units))
                     ELSE
                        CALL write_history_mt_ratio_2d (DEF_hist_vars%wa, &
                           a_mt_wa_mass(imt_loc, :), a_water_wa, &
                           file_hist, trim(mt_varname), itime_in_file, &
                           filter, trim(mt_longname), trim(mt_ratio_units))
                     ENDIF

                     write(mt_varname , '(A,A)')   'f_mt_conc_wa_debt_', trim(materials(imt_loc)%name)
                     write(mt_longname, '(5A)') 'aquifer debt material ', trim(mt_ratio_word), &
                        ' (', trim(materials(imt_loc)%name), ')'
                     CALL write_history_mt_ratio_2d (DEF_hist_vars%wa, &
                        a_mt_wa_debt_mass(imt_loc, :), a_water_wa_debt, &
                        file_hist, trim(mt_varname), itime_in_file, &
                        filter, trim(mt_longname), trim(mt_ratio_units))

                  ! --- Surface water (wdsrf) ratio ---
                  write(mt_varname , '(A,A)')   'f_mt_conc_wdsrf_', trim(materials(imt_loc)%name)
                  write(mt_longname, '(5A)') 'surface water material ', trim(mt_ratio_word), &
                     ' (', trim(materials(imt_loc)%name), ')'
                  CALL write_history_mt_ratio_2d (DEF_hist_vars%wdsrf, &
                     a_mt_wdsrf_mass(imt_loc, :), a_water_wdsrf, &
                     file_hist, trim(mt_varname), itime_in_file, &
                     filter, trim(mt_longname), trim(mt_ratio_units))

                  ! --- Wetland pool (wetwat) ratio ---
                  write(mt_varname , '(A,A)')   'f_mt_conc_wetwat_', trim(materials(imt_loc)%name)
                  write(mt_longname, '(5A)') 'wetland pool material ', trim(mt_ratio_word), &
                     ' (', trim(materials(imt_loc)%name), ')'
                  IF (p_is_worker) THEN
                     IF (numpatch > 0) THEN
                        filter(:) = patchtype == 2
                        IF (forcing_has_missing_value) THEN
                           filter = filter .and. forcmask_pch
                        ENDIF
                        filter = filter .and. patchmask
                     ENDIF
                  ENDIF
                  CALL write_history_mt_ratio_2d (DEF_hist_vars%wetwat, &
                     a_mt_wetwat_mass(imt_loc, :), a_water_wetwat, &
                     file_hist, trim(mt_varname), itime_in_file, &
                     filter, trim(mt_longname), trim(mt_ratio_units))

                  IF (p_is_worker) THEN
                     IF (numpatch > 0) THEN
                        filter(:) = patchtype < 99
                        IF (forcing_has_missing_value) filter = filter .and. forcmask_pch
                        filter = filter .and. patchmask
                     ENDIF
                  ENDIF

                  ! --- Pre-layer thin snow (scv) ratio ---
                  ! Only meaningful when snl==0 (accumulated snow below
                  ! layer-creation threshold). CoLM's water-side scv remains
                  ! total SWE after layered snow exists, but material mass then
                  ! lives in mt_wliq/mt_wice snow layers; the accumulator
                  ! therefore filters out snl<0 states so this ratio cannot
                  ! divide thin-snow material by total layered-snow SWE.
                  write(mt_varname , '(A,A)')   'f_mt_conc_scv_', trim(materials(imt_loc)%name)
                  write(mt_longname, '(5A)') 'thin-snow (scv) material ', trim(mt_ratio_word), &
                     ' (', trim(materials(imt_loc)%name), ')'
                  CALL write_history_mt_ratio_2d (DEF_hist_vars%scv, &
                     a_mt_scv_mass(imt_loc, :), a_water_scv, &
                     file_hist, trim(mt_varname), itime_in_file, &
                     filter, trim(mt_longname), trim(mt_ratio_units))

               ENDDO

               deallocate (mt_mass_soisno, water_soisno, mt_delta_flux, &
                           snowpack_mass_pat, snowpack_water_pat, ones_arr, &
                           leaf_mass_pat, leaf_water_pat)
            ENDIF

      CALL material_lifecycle_land_history (file_hist, itime_in_file, sumarea, filter, &
         nl_soil, forcing_has_missing_value, forcmask_pch)

#ifdef SinglePoint
      USE_SITE_HistWriteBack = site_hist_writeback_saved
#endif

   END SUBROUTINE mt_hist_out




#ifdef MaterialTransport
   ! Ratio writers (2d, 3d and their vector forms) drop a patch whose
   ! accumulated water is at most mt_water_min_for_delta, as the delta
   ! writers do. With mt_tiny they wrote the ratio of rounding dust: a
   ! NEC January IsoGSM run gave f_mt_conc_wdsrf = 0 (delta -1000) where
   ! the accumulated surface water was 1e-17..1e-14 mm.
   SUBROUTINE write_history_mt_ratio_2d ( is_hist, &
         material_mass_acc, water_acc, file_hist, varname, itime_in_file, &
         filter, longname, units)

   USE MOD_Block
   USE MOD_SPMD_Task, only: p_is_worker, p_is_io
#ifdef MaterialTransport
#endif
   USE MOD_Namelist, only: DEF_HIST_CompressLevel

   IMPLICIT NONE

   logical, intent(in) :: is_hist
   real(r8), intent(in) :: material_mass_acc(:)
   real(r8), intent(in) :: water_acc(:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer, intent(in) :: itime_in_file
   logical, intent(in) :: filter(:)
   character(len=*), intent(in) :: longname
   character(len=*), intent(in) :: units

   real(r8), allocatable :: material_map_vec(:)
   real(r8), allocatable :: water_map_vec(:)
   real(r8), allocatable :: ratio_vec(:)
   type(block_data_real8_2d) :: material_xy_2d
   type(block_data_real8_2d) :: water_xy_2d
   type(block_data_real8_2d) :: ratio_xy_2d
   integer :: ip, iblkme, xblk, yblk, xloc, yloc
   integer :: compress

      IF (.not. is_hist) RETURN
#ifdef SinglePoint
      IF (HistForm == 'Single') THEN
         allocate(ratio_vec(size(water_acc)))
         ratio_vec = spval
         DO ip = 1, size(water_acc)
            IF (abs(water_acc(ip)) > mt_water_min_for_delta .and. material_mass_acc(ip) /= spval) THEN
               ratio_vec(ip) = material_mass_acc(ip) / water_acc(ip)
            ENDIF
         ENDDO
         CALL single_write_2d (ratio_vec, file_hist, varname, itime_in_file, &
            longname, units)
         deallocate(ratio_vec)
         RETURN
      ENDIF
#endif
#ifdef UNSTRUCTURED
      IF (HistForm == 'Vector') THEN
         CALL write_history_mt_ratio_vector_2d (material_mass_acc, water_acc, &
            file_hist, varname, itime_in_file, filter, longname, units)
         RETURN
      ENDIF
#endif
      IF (HistForm /= 'Gridded') RETURN

      allocate(material_map_vec(size(material_mass_acc)))
      allocate(water_map_vec(size(water_acc)))
      material_map_vec = 0._r8
      water_map_vec = 0._r8

      IF (p_is_worker) THEN
         DO ip = 1, size(water_acc)
            IF (abs(water_acc(ip)) > mt_water_min_for_delta .and. material_mass_acc(ip) /= spval) THEN
               material_map_vec(ip) = material_mass_acc(ip)
               water_map_vec(ip) = water_acc(ip)
            ENDIF
         ENDDO
      ENDIF

      IF (p_is_io) THEN
         CALL allocate_block_data (ghist, material_xy_2d)
         CALL allocate_block_data (ghist, water_xy_2d)
         CALL allocate_block_data (ghist, ratio_xy_2d)
      ENDIF

      CALL mp2g_hist%pset2grid (material_map_vec, material_xy_2d, spv = spval, msk = filter)
      CALL mp2g_hist%pset2grid (water_map_vec, water_xy_2d, spv = spval, msk = filter)

      IF (p_is_io) THEN
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            DO yloc = 1, ghist%ycnt(yblk)
               DO xloc = 1, ghist%xcnt(xblk)
                  IF (water_xy_2d%blk(xblk,yblk)%val(xloc,yloc) /= spval .and. &
                      material_xy_2d%blk(xblk,yblk)%val(xloc,yloc) /= spval .and. &
                      abs(water_xy_2d%blk(xblk,yblk)%val(xloc,yloc)) > mt_tiny) THEN
                     ratio_xy_2d%blk(xblk,yblk)%val(xloc,yloc) = &
                        material_xy_2d%blk(xblk,yblk)%val(xloc,yloc) &
                        / water_xy_2d%blk(xblk,yblk)%val(xloc,yloc)
                  ELSE
                     ratio_xy_2d%blk(xblk,yblk)%val(xloc,yloc) = spval
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      compress = DEF_HIST_CompressLevel
      CALL hist_write_var_real8_2d (file_hist, varname, ghist, itime_in_file, &
         ratio_xy_2d, compress, longname, units)

      deallocate(material_map_vec)
      deallocate(water_map_vec)

   END SUBROUTINE write_history_mt_ratio_2d


   SUBROUTINE write_history_mt_ratio_3d ( is_hist, &
         material_mass_acc, water_acc, file_hist, varname, itime_in_file, &
         dim1name, lb1, ndim1, filter, longname, units)

   USE MOD_Block
   USE MOD_SPMD_Task, only: p_is_worker, p_is_io
   USE MOD_Vars_1DAccFluxes, only: nac
#ifdef MaterialTransport
#endif
   USE MOD_Namelist, only: DEF_HIST_CompressLevel

   IMPLICIT NONE

   logical, intent(in) :: is_hist
   real(r8), intent(in) :: material_mass_acc(:,:)
   real(r8), intent(in) :: water_acc(:,:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer, intent(in) :: itime_in_file
   character(len=*), intent(in) :: dim1name
   integer, intent(in) :: lb1, ndim1
   logical, intent(in) :: filter(:)
   character(len=*), intent(in) :: longname
   character(len=*), intent(in) :: units

   real(r8), allocatable :: material_map_vec(:,:)
   real(r8), allocatable :: water_map_vec(:,:)
   type(block_data_real8_3d) :: material_xy_3d
   type(block_data_real8_3d) :: water_xy_3d
   type(block_data_real8_3d) :: ratio_xy_3d
   integer :: ip, i1, iblkme, xblk, yblk, xloc, yloc
   integer :: lb_vec, ub_vec
   integer :: compress
   real(r8), allocatable :: ratio_vec(:,:)

      IF (.not. is_hist) RETURN
#ifdef SinglePoint
      IF (HistForm == 'Single') THEN
         lb_vec = lbound(material_mass_acc, 1)
         ub_vec = ubound(material_mass_acc, 1)
         allocate(ratio_vec(lb_vec:ub_vec, size(water_acc, 2)))
         ratio_vec = spval
         DO ip = 1, size(water_acc, 2)
            DO i1 = lb_vec, ub_vec
               IF (abs(water_acc(i1, ip)) > mt_water_min_for_delta .and. material_mass_acc(i1, ip) /= spval) THEN
                  ! single_write_3d divides by nac internally like ordinary
                  ! accumulated variables. Ratio is already a final diagnostic,
                  ! so pre-scale here to leave the written value unchanged.
                  ratio_vec(i1, ip) = material_mass_acc(i1, ip) / water_acc(i1, ip) * real(nac, r8)
               ENDIF
            ENDDO
         ENDDO
         CALL single_write_3d (ratio_vec, file_hist, varname, itime_in_file, &
            dim1name, ndim1, longname, units)
         deallocate(ratio_vec)
         RETURN
      ENDIF
#endif
#ifdef UNSTRUCTURED
      IF (HistForm == 'Vector') THEN
         CALL write_history_mt_ratio_vector_3d (material_mass_acc, water_acc, &
            file_hist, varname, itime_in_file, dim1name, lb1, ndim1, filter, &
            longname, units)
         RETURN
      ENDIF
#endif
      IF (HistForm /= 'Gridded') RETURN

      lb_vec = lbound(material_mass_acc, 1)
      ub_vec = ubound(material_mass_acc, 1)
      allocate(material_map_vec(lb_vec:ub_vec, size(material_mass_acc, 2)))
      allocate(water_map_vec  (lb_vec:ub_vec, size(water_acc, 2)))
      material_map_vec = spval
      water_map_vec  = spval

      IF (p_is_worker) THEN
         DO ip = 1, size(water_acc, 2)
            DO i1 = lb_vec, ub_vec
               IF (water_acc(i1, ip) > mt_water_min_for_delta .and. material_mass_acc(i1, ip) /= spval) THEN
                  material_map_vec(i1, ip) = material_mass_acc(i1, ip)
                  water_map_vec  (i1, ip) = water_acc(i1, ip)
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      IF (p_is_io) THEN
         CALL allocate_block_data (ghist, material_xy_3d, ndim1, lb1)
         CALL allocate_block_data (ghist, water_xy_3d,   ndim1, lb1)
         CALL allocate_block_data (ghist, ratio_xy_3d,   ndim1, lb1)
      ENDIF

      CALL mp2g_hist%pset2grid (material_map_vec, material_xy_3d, spv = spval, msk = filter)
      CALL mp2g_hist%pset2grid (water_map_vec,   water_xy_3d, spv = spval, msk = filter)

      IF (p_is_io) THEN
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            DO yloc = 1, ghist%ycnt(yblk)
               DO xloc = 1, ghist%xcnt(xblk)
                  DO i1 = ratio_xy_3d%lb1, ratio_xy_3d%ub1
                     IF (water_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc) /= spval .and. &
                         material_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc) /= spval .and. &
                         water_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc) > mt_tiny) THEN
                        ratio_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc) = &
                           material_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc) &
                           / water_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc)
                     ELSE
                        ratio_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc) = spval
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      compress = DEF_HIST_CompressLevel
      CALL hist_write_var_real8_3d (file_hist, varname, dim1name, ghist, &
         itime_in_file, ratio_xy_3d, compress, longname, units)

      deallocate(material_map_vec)
      deallocate(water_map_vec)

   END SUBROUTINE write_history_mt_ratio_3d


   SUBROUTINE write_history_mt_delta_2d ( is_hist, &
         material_mass_acc, water_acc, ref_ratio, file_hist, varname, itime_in_file, &
         filter, longname, units, water_min_override)

   USE MOD_Block
   USE MOD_SPMD_Task, only: p_is_worker, p_is_io
#ifdef MaterialTransport
#endif
   USE MOD_Namelist, only: DEF_HIST_CompressLevel

   IMPLICIT NONE

   logical, intent(in) :: is_hist
   real(r8), intent(in) :: material_mass_acc(:)
   real(r8), intent(in) :: water_acc(:)
   real(r8), intent(in) :: ref_ratio
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer, intent(in) :: itime_in_file
      logical, intent(in) :: filter(:)
      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units
      real(r8), intent(in), optional :: water_min_override

      real(r8), allocatable :: material_map_vec(:)
      real(r8), allocatable :: water_map_vec(:)
      real(r8), allocatable :: delta_vec(:)
      type(block_data_real8_2d) :: material_xy_2d
      type(block_data_real8_2d) :: water_xy_2d
      type(block_data_real8_2d) :: delta_xy_2d
      integer :: ip, iblkme, xblk, yblk, xloc, yloc
      integer :: compress
      real(r8) :: delta_loc, water_min

      IF (.not. is_hist) RETURN
      water_min = mt_water_min_for_delta
      IF (present(water_min_override)) water_min = water_min_override
#ifdef SinglePoint
      IF (HistForm == 'Single') THEN
         allocate(delta_vec(size(water_acc)))
         delta_vec = spval
         DO ip = 1, size(water_acc)
            IF (water_acc(ip) > water_min .and. material_mass_acc(ip) /= spval) THEN
               delta_loc = mass_to_delta(material_mass_acc(ip), water_acc(ip), ref_ratio)
               IF (delta_loc /= spval .and. abs(delta_loc) <= mt_delta_sanity_max) THEN
                  delta_vec(ip) = delta_loc
               ENDIF
            ENDIF
         ENDDO
         CALL single_write_2d (delta_vec, file_hist, varname, itime_in_file, &
            longname, units)
         deallocate(delta_vec)
         RETURN
      ENDIF
#endif
#ifdef UNSTRUCTURED
      IF (HistForm == 'Vector') THEN
         CALL write_history_mt_delta_vector_2d (material_mass_acc, water_acc, ref_ratio, &
            file_hist, varname, itime_in_file, filter, longname, units, water_min)
         RETURN
      ENDIF
#endif
      IF (HistForm /= 'Gridded') RETURN

      allocate(material_map_vec(size(material_mass_acc)))
         allocate(water_map_vec(size(water_acc)))
         material_map_vec = 0._r8
         water_map_vec = 0._r8

         IF (p_is_worker) THEN
            DO ip = 1, size(water_acc)
               IF (water_acc(ip) > water_min .and. &
                   material_mass_acc(ip) /= spval) THEN
               material_map_vec(ip) = material_mass_acc(ip)
               water_map_vec(ip) = water_acc(ip)
            ENDIF
         ENDDO
      ENDIF

      IF (p_is_io) THEN
         CALL allocate_block_data (ghist, material_xy_2d)
         CALL allocate_block_data (ghist, water_xy_2d)
         CALL allocate_block_data (ghist, delta_xy_2d)
      ENDIF

      CALL mp2g_hist%pset2grid (material_map_vec, material_xy_2d, spv = spval, msk = filter)
      CALL mp2g_hist%pset2grid (water_map_vec, water_xy_2d, spv = spval, msk = filter)

      IF (p_is_io) THEN
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            DO yloc = 1, ghist%ycnt(yblk)
               DO xloc = 1, ghist%xcnt(xblk)
                  IF (water_xy_2d%blk(xblk,yblk)%val(xloc,yloc) > mt_tiny) THEN
                     delta_loc = mass_to_delta( &
                        material_xy_2d%blk(xblk,yblk)%val(xloc,yloc), &
                        water_xy_2d%blk(xblk,yblk)%val(xloc,yloc), &
                        ref_ratio)
                     IF (delta_loc /= spval .and. abs(delta_loc) <= mt_delta_sanity_max) THEN
                        delta_xy_2d%blk(xblk,yblk)%val(xloc,yloc) = delta_loc
                     ELSE
                        delta_xy_2d%blk(xblk,yblk)%val(xloc,yloc) = spval
                     ENDIF
                  ELSE
                     delta_xy_2d%blk(xblk,yblk)%val(xloc,yloc) = spval
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      compress = DEF_HIST_CompressLevel
      CALL hist_write_var_real8_2d (file_hist, varname, ghist, itime_in_file, &
         delta_xy_2d, compress, longname, units)

      deallocate(material_map_vec)
      deallocate(water_map_vec)

   END SUBROUTINE write_history_mt_delta_2d

#ifdef UNSTRUCTURED
   SUBROUTINE write_history_mt_ratio_vector_2d (material_mass_acc, water_acc, &
      file_hist, varname, itime_in_file, filter, longname, units)

      IMPLICIT NONE
      real(r8), intent(in) :: material_mass_acc(:)
      real(r8), intent(in) :: water_acc(:)
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      integer, intent(in) :: itime_in_file
      logical, intent(in) :: filter(:)
      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units

      CALL write_history_mt_vector_2d (material_mass_acc, water_acc, 0._r8, &
         .false., mt_water_min_for_delta, file_hist, varname, itime_in_file, filter, &
         longname, units)
   END SUBROUTINE write_history_mt_ratio_vector_2d

   SUBROUTINE write_history_mt_delta_vector_2d (material_mass_acc, water_acc, &
      ref_ratio, file_hist, varname, itime_in_file, filter, longname, units, water_min)

      IMPLICIT NONE
      real(r8), intent(in) :: material_mass_acc(:)
      real(r8), intent(in) :: water_acc(:)
      real(r8), intent(in) :: ref_ratio
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      integer, intent(in) :: itime_in_file
      logical, intent(in) :: filter(:)
      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units
      real(r8), intent(in) :: water_min

      CALL write_history_mt_vector_2d (material_mass_acc, water_acc, ref_ratio, &
         .true., water_min, file_hist, varname, itime_in_file, filter, &
         longname, units)
   END SUBROUTINE write_history_mt_delta_vector_2d

   SUBROUTINE write_history_mt_vector_2d (material_mass_acc, water_acc, ref_ratio, &
      write_delta, water_min, file_hist, varname, itime_in_file, filter, longname, units)

      USE MOD_SPMD_Task
      USE MOD_Namelist, only: DEF_HIST_CompressLevel
      USE MOD_LandElm
      USE MOD_ElmVector
      USE MOD_LandPatch, only: elm_patch
      IMPLICIT NONE

      real(r8), intent(in) :: material_mass_acc(:)
      real(r8), intent(in) :: water_acc(:)
      real(r8), intent(in) :: ref_ratio
      logical, intent(in) :: write_delta
      real(r8), intent(in) :: water_min
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      integer, intent(in) :: itime_in_file
      logical, intent(in) :: filter(:)
      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units

      integer :: numset, totalnumset, iset, istt, iend, iwork, mesg(2)
      integer :: isrc, ndata, compress
      logical, allocatable :: mask(:)
      real(r8), allocatable :: frac(:)
      real(r8), allocatable :: mt_local(:), water_local(:), out_vec(:)
      real(r8), allocatable :: mt_global(:), water_global(:)
      real(r8), allocatable :: send_pair(:,:), recv_pair(:,:)
      real(r8) :: delta_loc

      IF (p_is_worker) THEN
         numset = numelm
         IF (numset > 0) THEN
            allocate(mt_local(numset), water_local(numset))
            mt_local = 0._r8
            water_local = 0._r8

            DO iset = 1, numset
               istt = elm_patch%substt(iset)
               iend = elm_patch%subend(iset)
               IF ((istt > 0) .and. (iend >= istt)) THEN
                  allocate(mask(istt:iend), frac(istt:iend))
                  IF (write_delta) THEN
                     mask = (water_acc(istt:iend) > water_min) .and. &
                        (material_mass_acc(istt:iend) /= spval) .and. filter(istt:iend)
                  ELSE
                     mask = (abs(water_acc(istt:iend)) > water_min) .and. &
                        (material_mass_acc(istt:iend) /= spval) .and. filter(istt:iend)
                  ENDIF
                  IF (any(mask)) THEN
                     frac = elm_patch%subfrc(istt:iend)
                     mt_local(iset) = sum(frac * material_mass_acc(istt:iend), mask = mask)
                     water_local(iset) = sum(frac * water_acc(istt:iend), mask = mask)
                  ENDIF
                  deallocate(mask, frac)
               ENDIF
            ENDDO
         ENDIF

#ifdef USEMPI
         mesg = (/p_iam_glb, numset/)
         CALL mpi_send (mesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         IF (numset > 0) THEN
            allocate(send_pair(2, numset))
            send_pair(1, :) = mt_local
            send_pair(2, :) = water_local
            CALL mpi_send (send_pair, 2*numset, MPI_REAL8, p_address_master, &
               mpi_tag_data, p_comm_glb, p_err)
            deallocate(send_pair)
         ENDIF
#endif
      ENDIF

      IF (p_is_master) THEN
         totalnumset = totalnumelm
         allocate(out_vec(totalnumset), mt_global(totalnumset), water_global(totalnumset))
         out_vec = spval
         mt_global = 0._r8
         water_global = 0._r8

#ifdef USEMPI
         DO iwork = 0, p_np_worker-1
            CALL mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)
            isrc = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               allocate(recv_pair(2, ndata))
               CALL mpi_recv (recv_pair, 2*ndata, MPI_REAL8, isrc, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)
               mt_global(elm_data_address(p_itis_worker(isrc))%val) = recv_pair(1, :)
               water_global(elm_data_address(p_itis_worker(isrc))%val) = recv_pair(2, :)
               deallocate(recv_pair)
            ENDIF
         ENDDO
#else
         IF (allocated(mt_local)) THEN
            mt_global(elm_data_address(0)%val) = mt_local
            water_global(elm_data_address(0)%val) = water_local
         ENDIF
#endif

         DO iset = 1, totalnumset
            IF (write_delta) THEN
               IF (water_global(iset) > water_min) THEN
                  delta_loc = mass_to_delta(mt_global(iset), water_global(iset), ref_ratio)
                  IF (delta_loc /= spval .and. abs(delta_loc) <= mt_delta_sanity_max) THEN
                     out_vec(iset) = delta_loc
                  ENDIF
               ENDIF
            ELSEIF (abs(water_global(iset)) > water_min) THEN
               out_vec(iset) = mt_global(iset) / water_global(iset)
            ENDIF
         ENDDO

         compress = DEF_HIST_CompressLevel
         IF (itime_in_file >= 1) THEN
            CALL ncio_write_serial_time (file_hist, varname, itime_in_file, out_vec, &
               'element', 'time', compress)
         ELSE
            CALL ncio_write_serial (file_hist, varname, out_vec, 'element', compress)
         ENDIF
         IF (itime_in_file <= 1) THEN
            CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
            CALL ncio_put_attr (file_hist, varname, 'units', units)
            CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
         ENDIF
      ENDIF

      IF (allocated(out_vec)) deallocate(out_vec)
      IF (allocated(mt_local)) deallocate(mt_local)
      IF (allocated(mt_global)) deallocate(mt_global)
      IF (allocated(water_local)) deallocate(water_local)
      IF (allocated(water_global)) deallocate(water_global)

#ifdef USEMPI
      ! Fixed tags are reused for every history field.  Keep a field-boundary
      ! synchronization so a fast worker cannot send the next field while the
      ! master is still receiving the current field from a slower worker.
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
   END SUBROUTINE write_history_mt_vector_2d

   SUBROUTINE write_history_mt_ratio_vector_3d (material_mass_acc, water_acc, &
      file_hist, varname, itime_in_file, dim1name, lb1, ndim1, filter, longname, units)

      USE MOD_SPMD_Task
      USE MOD_Namelist, only: DEF_HIST_CompressLevel
      USE MOD_LandElm
      USE MOD_ElmVector
      USE MOD_LandPatch, only: elm_patch
      IMPLICIT NONE

      real(r8), intent(in) :: material_mass_acc(lb1:,:)
      real(r8), intent(in) :: water_acc(lb1:,:)
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      integer, intent(in) :: itime_in_file
      character(len=*), intent(in) :: dim1name
      integer, intent(in) :: lb1, ndim1
      logical, intent(in) :: filter(:)
      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units

      integer :: numset, totalnumset, iset, istt, iend, iwork, mesg(2)
      integer :: isrc, ndata, compress, i1, k, ub1
      logical, allocatable :: mask(:)
      real(r8), allocatable :: frac(:)
      real(r8), allocatable :: mt_local(:,:), water_local(:,:), out_vec(:,:)
      real(r8), allocatable :: mt_global(:,:), water_global(:,:)
      real(r8), allocatable :: send_pair(:,:,:), recv_pair(:,:,:)

      ub1 = lb1 + ndim1 - 1
      IF (p_is_worker) THEN
         numset = numelm
         IF (numset > 0) THEN
            allocate(mt_local(ndim1, numset), water_local(ndim1, numset))
            mt_local = 0._r8
            water_local = 0._r8
            DO iset = 1, numset
               istt = elm_patch%substt(iset)
               iend = elm_patch%subend(iset)
               IF ((istt > 0) .and. (iend >= istt)) THEN
                  allocate(mask(istt:iend), frac(istt:iend))
                  frac = elm_patch%subfrc(istt:iend)
                  DO i1 = lb1, ub1
                     k = i1 - lb1 + 1
                     mask = (abs(water_acc(i1, istt:iend)) > mt_water_min_for_delta) .and. &
                        (material_mass_acc(i1, istt:iend) /= spval) .and. filter(istt:iend)
                     IF (any(mask)) THEN
                        mt_local(k, iset) = sum(frac * material_mass_acc(i1, istt:iend), mask = mask)
                        water_local(k, iset) = sum(frac * water_acc(i1, istt:iend), mask = mask)
                     ENDIF
                  ENDDO
                  deallocate(mask, frac)
               ENDIF
            ENDDO
         ENDIF

#ifdef USEMPI
         mesg = (/p_iam_glb, numset/)
         CALL mpi_send (mesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         IF (numset > 0) THEN
            allocate(send_pair(2, ndim1, numset))
            send_pair(1, :, :) = mt_local
            send_pair(2, :, :) = water_local
            CALL mpi_send (send_pair, 2*ndim1*numset, MPI_REAL8, p_address_master, &
               mpi_tag_data, p_comm_glb, p_err)
            deallocate(send_pair)
         ENDIF
#endif
      ENDIF

      IF (p_is_master) THEN
         totalnumset = totalnumelm
         allocate(out_vec(ndim1, totalnumset), mt_global(ndim1, totalnumset), &
            water_global(ndim1, totalnumset))
         out_vec = spval
         mt_global = 0._r8
         water_global = 0._r8

#ifdef USEMPI
         DO iwork = 0, p_np_worker-1
            CALL mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)
            isrc = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               allocate(recv_pair(2, ndim1, ndata))
               CALL mpi_recv (recv_pair, 2*ndim1*ndata, MPI_REAL8, isrc, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)
               DO k = 1, ndim1
                  mt_global(k, elm_data_address(p_itis_worker(isrc))%val) = recv_pair(1, k, :)
                  water_global(k, elm_data_address(p_itis_worker(isrc))%val) = recv_pair(2, k, :)
               ENDDO
               deallocate(recv_pair)
            ENDIF
         ENDDO
#else
         IF (allocated(mt_local)) THEN
            mt_global(:, elm_data_address(0)%val) = mt_local
            water_global(:, elm_data_address(0)%val) = water_local
         ENDIF
#endif

         DO iset = 1, totalnumset
            DO k = 1, ndim1
               IF (abs(water_global(k, iset)) > mt_tiny) THEN
                  out_vec(k, iset) = mt_global(k, iset) / water_global(k, iset)
               ENDIF
            ENDDO
         ENDDO

         CALL ncio_define_dimension (file_hist, dim1name, ndim1)
         compress = DEF_HIST_CompressLevel
         IF (itime_in_file >= 1) THEN
            CALL ncio_write_serial_time (file_hist, varname, itime_in_file, out_vec, &
               dim1name, 'element', 'time', compress)
         ELSE
            CALL ncio_write_serial (file_hist, varname, out_vec, dim1name, 'element', compress)
         ENDIF
         IF (itime_in_file <= 1) THEN
            CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
            CALL ncio_put_attr (file_hist, varname, 'units', units)
            CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
         ENDIF
      ENDIF

      IF (allocated(out_vec)) deallocate(out_vec)
      IF (allocated(mt_local)) deallocate(mt_local)
      IF (allocated(mt_global)) deallocate(mt_global)
      IF (allocated(water_local)) deallocate(water_local)
      IF (allocated(water_global)) deallocate(water_global)

#ifdef USEMPI
      ! See the 2-D gather above: the trailing barrier separates messages
      ! that intentionally reuse mpi_tag_mesg/mpi_tag_data across fields.
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
   END SUBROUTINE write_history_mt_ratio_vector_3d
#endif
#endif

   SUBROUTINE write_history_variable_2d ( is_hist, &
         acc_vec, file_hist, varname, itime_in_file, sumarea, filter, &
         longname, units, acc_num, input_mode)

   USE MOD_Vars_1DAccFluxes, only: nac

   IMPLICIT NONE

   logical, intent(in) :: is_hist

   real(r8),         intent(inout) :: acc_vec(:)
   character(len=*), intent(in)    :: file_hist
   character(len=*), intent(in)    :: varname
   integer,          intent(in)    :: itime_in_file
   character(len=*), intent(in)    :: longname
   character(len=*), intent(in)    :: units

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)
   real(r8), intent(in), optional  :: acc_num(:)

   character(len=*), intent(in), optional :: input_mode

      IF (.not. is_hist) RETURN

#ifndef SinglePoint
      IF ( .not. present(acc_num) ) THEN
         IF (p_is_worker) &
            WHERE (acc_vec /= spval) acc_vec = acc_vec / nac
      ELSE
         IF (p_is_worker) THEN
            WHERE (acc_vec/=spval .and. acc_num>0)
               acc_vec = acc_vec / acc_num
            ELSEWHERE (acc_vec/=spval .and. acc_num<=0)
               acc_vec = spval
            END WHERE
         ENDIF
      ENDIF
#else
      IF ( .not. present(acc_num) ) THEN
         WHERE (acc_vec /= spval)  acc_vec = acc_vec / nac
      ELSE
         WHERE (acc_vec/=spval .and. acc_num>0)
               acc_vec = acc_vec / acc_num
            ELSEWHERE (acc_vec/=spval .and. acc_num<=0)
               acc_vec = spval
            END WHERE
      ENDIF
#endif

      select CASE (HistForm)
      CASE ('Gridded')
         IF (present(input_mode)) THEN
            CALL flux_map_and_write_2d ( &
               acc_vec, file_hist, varname, itime_in_file, sumarea, filter, longname, units, input_mode)
         ELSE
            CALL flux_map_and_write_2d ( &
               acc_vec, file_hist, varname, itime_in_file, sumarea, filter, longname, units)
         ENDIF
#ifdef UNSTRUCTURED
      CASE ('Vector')
         IF (present(input_mode)) THEN
            CALL aggregate_to_vector_and_write_2d ( &
               acc_vec, file_hist, varname, itime_in_file, filter, longname, units, input_mode)
         ELSE
            CALL aggregate_to_vector_and_write_2d ( &
               acc_vec, file_hist, varname, itime_in_file, filter, longname, units)
         ENDIF
#endif
#ifdef SinglePoint
      CASE ('Single')
         CALL single_write_2d ( &
            acc_vec, file_hist, varname, itime_in_file, longname, units)
#endif
      END select

   END SUBROUTINE write_history_variable_2d


#ifdef URBAN_MODEL
   SUBROUTINE write_history_variable_urb_2d ( is_hist, &
         acc_vec, file_hist, varname, itime_in_file, sumarea, filter, &
         longname, units)

   IMPLICIT NONE

   logical, intent(in) :: is_hist

   real(r8), intent(inout) :: acc_vec(:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer,          intent(in) :: itime_in_file
   character(len=*), intent(in) :: longname
   character(len=*), intent(in) :: units

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)

      IF (.not. is_hist) RETURN

      select CASE (HistForm)
      CASE ('Gridded')
         CALL flux_map_and_write_urb_2d ( &
            acc_vec, file_hist, varname, itime_in_file, sumarea, filter, longname, units)
#ifdef UNSTRUCTURED
      CASE ('Vector')
         !TODO: currently, it is not applicable to urban variables
         CALL aggregate_to_vector_and_write_2d ( &
            acc_vec, file_hist, varname, itime_in_file, filter, longname, units)
#endif
#ifdef SinglePoint
      CASE ('Single')
         CALL single_write_2d ( &
            acc_vec, file_hist, varname, itime_in_file, longname, units)
#endif
      END select

   END SUBROUTINE write_history_variable_urb_2d
#endif


   SUBROUTINE write_history_variable_3d ( is_hist, &
         acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, &
         sumarea, filter, longname, units, acc_num)

   USE MOD_Vars_1DAccFluxes, only: nac

   IMPLICIT NONE

   logical, intent(in) :: is_hist

   real(r8), intent(inout) :: acc_vec(:,:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer, intent(in) :: itime_in_file
   character(len=*), intent(in) :: dim1name
   integer, intent(in) :: lb1, ndim1

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)
   character (len=*), intent(in) :: longname
   character (len=*), intent(in) :: units
   real(r8), intent(in), optional :: acc_num(:)

   ! Local variables
   integer :: iblkme, xblk, yblk, xloc, yloc, i1
   integer :: compress
   integer :: j

      IF (.not. is_hist) RETURN

      ! The lower-level 3D writers (flux_map_and_write_3d, aggregate_..., single_...)
      ! always divide acc_vec by `nac`. When a per-patch acc_num counter is supplied
      ! (e.g. methane accumulators that only tick on active patches), pre-scale by
      ! nac/acc_num so the lower-level division yields acc_vec/acc_num.
      IF (present(acc_num)) THEN
#ifndef SinglePoint
         IF (p_is_worker) THEN
            DO j = 1, size(acc_vec, 1)
               WHERE (acc_vec(j,:) /= spval .and. acc_num > 0._r8)
                  acc_vec(j,:) = acc_vec(j,:) * real(nac, r8) / acc_num
               ELSEWHERE (acc_vec(j,:) /= spval .and. acc_num <= 0._r8)
                  acc_vec(j,:) = spval
               END WHERE
            ENDDO
         ENDIF
#else
         DO j = 1, size(acc_vec, 1)
            WHERE (acc_vec(j,:) /= spval .and. acc_num > 0._r8)
               acc_vec(j,:) = acc_vec(j,:) * real(nac, r8) / acc_num
            ELSEWHERE (acc_vec(j,:) /= spval .and. acc_num <= 0._r8)
               acc_vec(j,:) = spval
            END WHERE
         ENDDO
#endif
      ENDIF

      select CASE (HistForm)
      CASE ('Gridded')
         CALL flux_map_and_write_3d ( &
            acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, &
            sumarea, filter, longname, units)
#ifdef UNSTRUCTURED
      CASE ('Vector')
         CALL aggregate_to_vector_and_write_3d ( &
            acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, &
            filter, longname, units)
#endif
#ifdef SinglePoint
      CASE ('Single')
         CALL single_write_3d (acc_vec, file_hist, varname, itime_in_file, &
            dim1name, ndim1, longname, units)
#endif
      END select

   END SUBROUTINE write_history_variable_3d

END MODULE MOD_MT_Hist
#endif
