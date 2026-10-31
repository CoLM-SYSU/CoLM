#include <define.h>

#if (defined MaterialTransport) && (defined BGC)
MODULE MOD_MT_Methane
!=======================================================================
! CH4 reactive-material implementation behind the generic reactive dispatch.
!=======================================================================

   USE MOD_Precision
   USE MOD_SPMD_Task, only: p_is_worker, p_is_master, CoLM_stop
#ifdef USEMPI
   USE MPI
   USE MOD_SPMD_Task, only: p_comm_glb, p_err
#endif
   USE MOD_MT_Defs, only: material_param_file_for_index, material_lower, &
      FAMILY_GAS, STATE_OWNER_PROVIDER, REACTION_PROVIDER
   USE MOD_MT_Lifecycle, only: material_lifecycle_hooks_type, register_mt_provider
   USE MOD_MT_Methane_Const, only: DEF_file_GIEMS, DEF_wetland_finundation_scheme
   USE MOD_Vars_TimeInvariants, only: patchtype, patchlatr, patchlonr
   USE MOD_MT_Methane_Registry, only: igas_ch4
   USE MOD_MT_Methane_Physics,  only: methane_host_water_reset, &
      methane_host_water_report
   USE MOD_MT_Methane_State,    only: allocate_methane_state, &
      init_methane_wetland_fraction_cache, deallocate_methane_state, &
      read_methane_restart, write_methane_restart, initialize_methane_lake_soilc_from_surface, &
      methane_lake_soilc_srf, load_methane_lake_soilc_srf, &
      save_methane_lulcc_state, remap_methane_lulcc_state, &
      publish_methane_levee_flood_patch, publish_methane_flood_patch
   USE MOD_MT_Methane_AccFlux,  only: allocate_methane_acc_fluxes, &
      deallocate_methane_acc_fluxes, flush_methane_acc_fluxes, accumulate_methane_fluxes, &
      read_methane_accflux_restart, write_methane_accflux_restart
   USE MOD_MT_Methane_Microbes, only: allocate_methane_microbes_state, &
      deallocate_methane_microbes_state, read_methane_microbes_restart, &
      validate_methane_microbes_restart_values, &
      write_methane_microbes_restart, save_methane_microbes_lulcc_state, &
      remap_methane_microbes_lulcc_state
   USE MOD_MT_Methane_Const,    only: read_methane_namelist, &
      configure_methane_inundation_mode, methane_history_accumulation_mode, DEF_METHANE
   USE MOD_MT_Methane_GIEMS,    only: allocate_methane_giems, &
      deallocate_methane_giems, read_methane_giems, giems_active
   USE MOD_MT_Methane_pH,       only: allocate_methane_ph, &
      deallocate_methane_ph, read_methane_ph_patch
   USE MOD_MT_Methane_VegOverride, only: allocate_wetland_aere_overrides, &
      deallocate_wetland_aere_overrides
   USE MOD_MT_Methane_Impl, only: ch4_impl_lake_step, &
      ch4_impl_wetland_decomp, ch4_impl_soil_step
   USE MOD_MT_Methane_Hist, only: methane_reactive_history
   USE MOD_MT_Methane_WetlandBGC, only: methane_wetland_bgc_allocate, &
      methane_wetland_bgc_deallocate, methane_wetland_bgc_cold_init, &
      methane_wetland_bgc_read_restart, methane_wetland_bgc_write_restart, &
      methane_wetland_bgc_save_lulcc, methane_wetland_bgc_remap_lulcc

   IMPLICIT NONE
   PRIVATE

   integer, parameter :: METHANE_RESTART_SCHEMA_VERSION = 5
   character(len=512), save :: last_methane_ph_patch_file = ''

   PUBLIC :: ch4_reactive_init, ch4_reactive_final
   PUBLIC :: ch4_reactive_lake_step, ch4_reactive_wetland_decomp
   PUBLIC :: ch4_reactive_soil_step, ch4_reactive_report
   PUBLIC :: ch4_reactive_write_restart, ch4_reactive_read_restart
   PUBLIC :: ch4_reactive_flush_acc_fluxes, ch4_reactive_accumulate_fluxes
   PUBLIC :: ch4_reactive_save_lulcc_state, ch4_reactive_remap_lulcc_state
   PUBLIC :: ch4_reactive_reload_lulcc_inputs
   PUBLIC :: ch4_reactive_publish_levee_flood, ch4_reactive_publish_flood
   PUBLIC :: ch4_register_mt_provider

CONTAINS

   SUBROUTINE ch4_register_mt_provider ()

      IMPLICIT NONE
      type(material_lifecycle_hooks_type) :: hooks

      hooks = material_lifecycle_hooks_type()
      hooks%land_init            => ch4_reactive_init
      hooks%land_read_restart    => ch4_reactive_read_restart
      hooks%land_write_restart   => ch4_reactive_write_restart
      hooks%lake_step            => ch4_reactive_lake_step
      hooks%wetland_decomp       => ch4_reactive_wetland_decomp
      hooks%soil_step            => ch4_reactive_soil_step
      hooks%report               => ch4_reactive_report
      hooks%flush_acc_fluxes     => ch4_reactive_flush_acc_fluxes
      hooks%accumulate_fluxes    => ch4_reactive_accumulate_fluxes
      hooks%land_history         => methane_reactive_history
      hooks%land_history_sidecar_write => ch4_reactive_write_history_sidecar
      hooks%land_history_sidecar_read  => ch4_reactive_read_history_sidecar
      hooks%land_save_lulcc      => ch4_reactive_save_lulcc_state
      hooks%land_remap_lulcc     => ch4_reactive_remap_lulcc_state
      hooks%land_reload_lulcc    => ch4_reactive_reload_lulcc_inputs
      hooks%publish_levee_flood  => ch4_reactive_publish_levee_flood
      hooks%publish_flood        => ch4_reactive_publish_flood
      hooks%land_final           => ch4_reactive_final

      CALL register_mt_provider ('CH4', 'METHANE', 'methane', &
         FAMILY_GAS, STATE_OWNER_PROVIDER, REACTION_PROVIDER, hooks, igas_ch4)

   END SUBROUTINE ch4_register_mt_provider

   logical FUNCTION ch4_reactive_has ()

      IMPLICIT NONE

      ch4_reactive_has = (igas_ch4 > 0)

   END FUNCTION ch4_reactive_has

   SUBROUTINE ch4_reactive_init (numpatch, lc_year, jdate, casename, dir_restart, dir_landdata, file_restart)

      USE MOD_LandPatch, only: landpatch
      USE MOD_MT_NcIO, only: ncio_vector_group_presence
      IMPLICIT NONE
      integer, intent(in) :: numpatch
      integer, intent(in) :: lc_year
      integer, intent(in) :: jdate(3)
      character(len=*), intent(in) :: casename
      character(len=*), intent(in) :: dir_restart
      character(len=*), intent(in) :: dir_landdata
      character(len=*), intent(in), optional :: file_restart

      character(len=256) :: cyear_restart
      character(len=256) :: file_param
      logical :: use_param
      logical :: lake_restart_present(2)
      character(len=32), parameter :: lake_restart_fields(2) = &
         [character(len=32) :: 'ch4_conc_methane', 'ch4_lake_soilc']
      real(r8), allocatable :: giems_dummy_patch(:)

      IF (.not. ch4_reactive_has()) RETURN

      CALL ch4_reactive_resolve_param_file (use_param, file_param)
      IF (use_param) THEN
         CALL read_methane_namelist (file_param)
      END IF
      CALL configure_methane_inundation_mode ()

#ifndef CROP
      IF (DEF_METHANE%enable_rice_paddy) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: DEF_METHANE%enable_rice_paddy requires compiling with CROP.'
         CALL CoLM_stop (' ***** ERROR: rice-paddy methane requested without CROP')
      ENDIF
#endif

      ! ch4_reactive_init can be re-entered by LULCC/restart workflows.
      ! Drop any previous CH4 allocatables before allocating the current
      ! numpatch layout; the deallocators are no-ops when not allocated.
      CALL deallocate_methane_acc_fluxes ()
      CALL deallocate_wetland_aere_overrides ()
      CALL deallocate_methane_ph ()
      CALL deallocate_methane_giems ()
      CALL deallocate_methane_microbes_state ()
      CALL deallocate_methane_state ()
      CALL methane_wetland_bgc_deallocate ()

      CALL allocate_methane_state (numpatch)
      ! Wetland shadow pools: cold-start from the host column now (land_init
      ! runs after READ_TimeVariables); land_read_restart replaces them when
      ! the methane checkpoint carries them.
      CALL methane_wetland_bgc_allocate (numpatch)
      IF (p_is_worker .and. allocated(patchtype)) THEN
         CALL methane_wetland_bgc_cold_init (patchtype, numpatch)
      ENDIF
      ! Per-rank diagnostic counters must not survive a re-initialisation in
      ! the same process, or the second run inherits the first's worst case
      ! and its already-warned flag.
      CALL methane_host_water_reset ()
      CALL init_methane_wetland_fraction_cache (numpatch)
      IF (DEF_METHANE%use_microbial_pools) THEN
         CALL allocate_methane_microbes_state (numpatch)
      ENDIF
      CALL allocate_methane_acc_fluxes (numpatch)

      write(cyear_restart,'(i4.4)') lc_year
      ! Probe on all ranks, but leave transaction/value validation and the
      ! actual read to land_read_restart. A valid zero (exhausted) inventory
      ! is still authoritative and must not be refilled from surface data.
      CALL load_methane_lake_soilc_srf (lc_year, dir_landdata, numpatch)
      lake_restart_present = .false.
      IF (present(file_restart)) THEN
         CALL ncio_vector_group_presence(file_restart, lake_restart_fields, landpatch, lake_restart_present)
      ENDIF
      IF (p_is_worker .and. allocated(patchtype) .and. .not. all(lake_restart_present)) THEN
         CALL initialize_methane_lake_soilc_from_surface (patchtype, methane_lake_soilc_srf, DEF_METHANE%allowlakeprod)
      ENDIF

      CALL allocate_methane_giems (numpatch)
      IF (DEF_wetland_finundation_scheme == 5) THEN
         IF (trim(DEF_file_GIEMS) == 'null') THEN
            CALL CoLM_stop (' ***** ERROR: DEF_wetland_finundation_scheme=5 requires DEF_file_GIEMS.')
         ENDIF
         CALL require_methane_giems_patch_coords (numpatch)
         IF (p_is_worker .and. numpatch > 0) THEN
            CALL read_methane_giems (DEF_file_GIEMS, patchlatr, patchlonr, numpatch)
         ELSE
            allocate(giems_dummy_patch(0))
            CALL read_methane_giems (DEF_file_GIEMS, giems_dummy_patch, giems_dummy_patch, 0)
            deallocate(giems_dummy_patch)
         ENDIF
         IF (.not. giems_active) THEN
            CALL CoLM_stop (' ***** ERROR: GIEMS file could not be loaded for methane scheme 5.')
         ENDIF
      ENDIF

      last_methane_ph_patch_file = ''
      CALL allocate_methane_ph (numpatch)
      IF (DEF_METHANE%use_spatial_ph) THEN
         last_methane_ph_patch_file = trim(dir_landdata)//'/soil/'//trim(cyear_restart)// &
            '/methane_ph_patches.nc'
         CALL read_methane_ph_patch (trim(last_methane_ph_patch_file), numpatch)
      ENDIF

      CALL allocate_wetland_aere_overrides (numpatch)

   END SUBROUTINE ch4_reactive_init


   SUBROUTINE ch4_reactive_resolve_param_file (use_param, file_param)

      IMPLICIT NONE
      logical, intent(out) :: use_param
      character(len=*), intent(out) :: file_param

      logical :: found

      CALL material_param_file_for_index (igas_ch4, 'CH4,METHANE', file_param, found)
      use_param = found .and. len_trim(file_param) > 0 .and. trim(material_lower(file_param)) /= 'null'
      IF (.not. use_param) THEN
         CALL CoLM_stop (' ***** ERROR: CH4 requires DEF_METHANE_PARAM_FILE to include a CH4 parameter file')
      ENDIF

   END SUBROUTINE ch4_reactive_resolve_param_file

   SUBROUTINE ch4_reactive_read_restart (file_restart)

      USE MOD_LandPatch, only: landpatch
      USE MOD_MT_NcIO, only: ncio_vector_group_presence, ncio_set_complete_require_present
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer, parameter :: nmicrobe_restart_fields = 25
      character(len=40), parameter :: provider_restart_probe_fields(5) = &
         [character(len=40) :: 'ch4_conc_methane', 'ch4_restart_schema', &
          'ch4_restart_complete', 'ch4_history_accumulation_mode', &
          'ch4_feature_microbe_pools']
      character(len=40), parameter :: microbe_restart_fields(nmicrobe_restart_fields) = &
         [character(len=40) :: &
          'ch4_B_methanogen', 'ch4_B_methanotroph', &
          'ch4_B_methanogen_dormant', 'ch4_B_methanotroph_dormant', &
          'ch4_B_methanogen_soil', 'ch4_B_methanogen_rice', &
          'ch4_B_methanotroph_soil', 'ch4_B_methanotroph_rice', &
          'ch4_B_methanogen_dormant_soil', 'ch4_B_methanogen_dormant_rice', &
          'ch4_B_methanotroph_dormant_soil', 'ch4_B_methanotroph_dormant_rice', &
          'ch4_a_B_methanogen', 'ch4_a_B_methanotroph', &
          'ch4_a_B_methanogen_dormant', 'ch4_a_B_methanotroph_dormant', &
          'ch4_a_f_T_methanogen', 'ch4_a_f_S_methanogen', &
          'ch4_a_f_O2_methanogen', 'ch4_a_f_T_methanotroph', &
          'ch4_a_methanogen_growth_rate', 'ch4_a_methanotroph_growth_rate', &
          'ch4_a_microbial_prod_potential', 'ch4_a_microbial_oxid_potential', &
          'ch4_a_methane_acc_num_microbe']
      logical :: file_has_pools, file_has_component_pools
      logical :: file_has_microbe_accumulators, strict_restart
      logical :: history_mode_changed, checkpoint_has_microbe_pools
      logical :: microbe_field_present(nmicrobe_restart_fields)
      logical :: wetland_pools_present
      logical :: provider_restart_probe_present(5)
      integer :: n_microbe_state_fields, n_microbe_component_fields
      integer :: n_microbe_accumulator_fields
      integer :: restart_schema

      IF (.not. ch4_reactive_has()) RETURN
      ! Generic land-water restart compatibility is independent of this
      ! provider-owned checkpoint. Conversely, a cold start may have no CH4
      ! fields at all. Probe the mandatory legacy/current field together with
      ! transaction markers: no provider fields is a clean cold start, while
      ! marker-only state is an interrupted/corrupt provider checkpoint.
      CALL ncio_vector_group_presence(file_restart, provider_restart_probe_fields, &
         landpatch, provider_restart_probe_present)
      IF (.not. any(provider_restart_probe_present)) RETURN
      IF (.not. provider_restart_probe_present(1)) THEN
         IF (p_is_master) WRITE(*,'(A)') &
            'ERROR: methane restart metadata exists without mandatory ch4_conc_methane state.'
         CALL CoLM_stop()
      ENDIF
      CALL validate_methane_restart_transaction (file_restart, strict_restart, restart_schema, &
         history_mode_changed, checkpoint_has_microbe_pools)
      ! C4: warn on a microbial-pool restart<->runtime-flag mismatch. The read
      ! below is flag-gated, so without this a restart written with pools ON but
      ! resumed with use_microbial_pools=.false. silently drops the prognostic
      ! biomass pools (and the reverse silently cold-starts them from B_init),
      ! breaking restart reproducibility with no message. Master-only probe.
      IF (restart_schema >= 5 .and. .not. checkpoint_has_microbe_pools) THEN
         ! Schema 5 makes feature ownership explicit.  When the marker is off,
         ! any microbial variables left in a reused NetCDF target belong to an
         ! older transaction.  Do not even probe those stale variables: an old
         ! partial group is outside the committed feature set as well.
         n_microbe_state_fields = 0
         n_microbe_component_fields = 0
         n_microbe_accumulator_fields = 0
         file_has_pools = .false.
         file_has_component_pools = .false.
         file_has_microbe_accumulators = .false.
      ELSE
         CALL ncio_vector_group_presence(file_restart, microbe_restart_fields, &
            landpatch, microbe_field_present)
         n_microbe_state_fields = count(microbe_field_present(1:4))
         n_microbe_component_fields = count(microbe_field_present(5:12))
         n_microbe_accumulator_fields = count(microbe_field_present(13:nmicrobe_restart_fields))
         IF (restart_schema >= 5) THEN
            IF (n_microbe_state_fields /= 4 .or. n_microbe_component_fields /= 8 .or. &
                n_microbe_accumulator_fields /= nmicrobe_restart_fields - 12) THEN
               IF (p_is_master) WRITE(*,'(A)') &
                  'ERROR: committed methane restart enables microbial pools but lacks the complete feature group.'
               CALL CoLM_stop()
            ENDIF
            file_has_pools = .true.
            file_has_component_pools = .true.
            file_has_microbe_accumulators = .true.
         ELSE
            ! Legacy schemas predate the authoritative feature marker.  Preserve
            ! their historical field-presence inference and migration behavior.
            file_has_pools = n_microbe_state_fields == 4
            file_has_component_pools = n_microbe_component_fields == 8
            file_has_microbe_accumulators = &
               n_microbe_accumulator_fields == nmicrobe_restart_fields - 12
            IF (strict_restart .and. &
                ((n_microbe_state_fields /= 0 .and. n_microbe_state_fields /= 4) .or. &
                 (n_microbe_component_fields /= 0 .and. n_microbe_component_fields /= 8) .or. &
                 (n_microbe_accumulator_fields /= 0 .and. &
                  n_microbe_accumulator_fields /= nmicrobe_restart_fields - 12) .or. &
                 (file_has_pools .neqv. file_has_microbe_accumulators) .or. &
                 (file_has_component_pools .and. .not. file_has_pools) .or. &
                 (restart_schema >= 3 .and. (file_has_pools .neqv. file_has_component_pools)))) THEN
               IF (p_is_master) WRITE(*,'(A)') &
                  'ERROR: committed methane restart contains a partial/inconsistent microbial feature group.'
               CALL CoLM_stop()
            ENDIF
            IF (.not. strict_restart .and. &
                ((n_microbe_state_fields /= 0 .and. n_microbe_state_fields /= 4) .or. &
                 (n_microbe_component_fields /= 0 .and. n_microbe_component_fields /= 8) .or. &
                 (n_microbe_accumulator_fields /= 0 .and. &
                  n_microbe_accumulator_fields /= nmicrobe_restart_fields - 12))) THEN
               IF (p_is_master) WRITE(*,'(A)') &
                  'WARNING: legacy methane restart has a partial microbial feature group; missing fields cold-start.'
            ENDIF
         ENDIF
      ENDIF
      IF (p_is_master) THEN
         IF (file_has_pools .and. (.not. DEF_METHANE%use_microbial_pools)) THEN
            WRITE(*,'(A)') 'WARNING: methane microbial-pool fields are present in the '// &
               'restart but use_microbial_pools=.false.; those pools are being ignored.'
         ELSEIF ((.not. file_has_pools) .and. DEF_METHANE%use_microbial_pools) THEN
            WRITE(*,'(A)') 'WARNING: use_microbial_pools=.true. but the restart has no '// &
               'microbial-pool fields; biomass cold-started from B_init (not reproducible).'
         ENDIF
      ENDIF

      ! Core schema fields are mandatory.  Microbial pools are a feature-
      ! conditional group: if the checkpoint contains the group, strict reads
      ! require every member; if it does not, enabling pools cold-starts from
      ! B_init as documented above rather than reclassifying a valid checkpoint
      ! as corrupt.
      CALL ncio_set_complete_require_present (strict_restart)
      ! Pass the validated schema so migration defaults are restricted to
      ! fields introduced after schema 1.
      CALL read_methane_restart (file_restart, strict_restart, restart_schema)
      ! Wetland shadow pools; without them (a restart from before the shadow
      ! pools existed) the cold start made in ch4_reactive_init is kept.
      CALL methane_wetland_bgc_read_restart (file_restart, wetland_pools_present)
      IF (p_is_master .and. .not. wetland_pools_present) WRITE(*,'(A)') &
         'WARNING: methane restart has no wetland BGC shadow pools; they are cold-started.'
      ! The accumulator reader uses the same schema boundary for the six new
      ! lake accumulators.
      CALL read_methane_accflux_restart (file_restart, file_has_pools, &
         file_has_microbe_accumulators, strict_restart, restart_schema)
      IF (restart_schema >= 2 .and. history_mode_changed) THEN
         IF (p_is_master) WRITE(*,'(A)') &
            'WARNING: methane history accumulation mode changed across restart; resetting the partial window.'
         CALL flush_methane_acc_fluxes ()
      ENDIF
      IF (file_has_pools) THEN
         IF (DEF_METHANE%use_microbial_pools) THEN
            CALL read_methane_microbes_restart (file_restart, strict_restart, restart_schema >= 3, &
               file_has_component_pools)
         ELSE
            CALL validate_methane_microbes_restart_values(file_restart, strict_restart, &
               file_has_component_pools)
         ENDIF
      ENDIF
      CALL ncio_set_complete_require_present (.false.)

   END SUBROUTINE ch4_reactive_read_restart

   SUBROUTINE ch4_reactive_write_restart (file_restart, compress)

      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer, intent(in) :: compress

      IF (.not. ch4_reactive_has()) RETURN
      ! Invalidate an older commit before overwriting any member field.  The
      ! final marker is written only after state, accumulators and optional
      ! microbial pools have all completed on every vector block.
      CALL write_methane_restart_marker(file_restart, 'ch4_restart_complete', 0._r8, compress)
      CALL write_methane_restart_marker(file_restart, 'ch4_restart_schema', &
         real(METHANE_RESTART_SCHEMA_VERSION, r8), compress)
      CALL write_methane_restart_marker(file_restart, 'ch4_history_accumulation_mode', &
         real(methane_history_accumulation_mode(), r8), compress)
      CALL write_methane_restart_marker(file_restart, 'ch4_feature_microbe_pools', &
         merge(1._r8, 0._r8, DEF_METHANE%use_microbial_pools), compress)
      CALL write_methane_restart(file_restart, compress)
      CALL write_methane_accflux_restart(file_restart, compress)
      CALL write_methane_microbes_restart(file_restart, compress)
      CALL methane_wetland_bgc_write_restart(file_restart, compress)
      CALL write_methane_restart_marker(file_restart, 'ch4_restart_complete', 1._r8, compress)

   END SUBROUTINE ch4_reactive_write_restart

   SUBROUTINE ch4_reactive_write_history_sidecar (file_hist_acc)

      USE MOD_LandPatch, only: landpatch
      USE MOD_NetCDFVector, only: ncio_define_dimension_vector
      USE MOD_Vars_Global, only: nl_soil
      IMPLICIT NONE
      character(len=*), intent(in) :: file_hist_acc

      IF (.not. ch4_reactive_has()) RETURN
      CALL ncio_define_dimension_vector (file_hist_acc, landpatch, 'soil', nl_soil)
      CALL write_methane_accflux_restart (file_hist_acc, 0)

   END SUBROUTINE ch4_reactive_write_history_sidecar

   SUBROUTINE ch4_reactive_read_history_sidecar (file_hist_acc)

      USE MOD_LandPatch, only: landpatch
      USE MOD_MT_NcIO, only: ncio_vector_var_present
      IMPLICIT NONE
      character(len=*), intent(in) :: file_hist_acc

      IF (.not. ch4_reactive_has()) RETURN
      IF (.not. ncio_vector_var_present(file_hist_acc, 'ch4_a_methane_acc_num', landpatch)) RETURN
      CALL read_methane_accflux_restart (file_hist_acc, &
         checkpoint_has_microbe_pools = DEF_METHANE%use_microbial_pools, &
         checkpoint_has_microbe_accumulators = DEF_METHANE%use_microbial_pools, &
         restart_schema = METHANE_RESTART_SCHEMA_VERSION)

   END SUBROUTINE ch4_reactive_read_history_sidecar

   SUBROUTINE write_methane_restart_marker (file_restart, varname, value, compress)

      USE MOD_LandPatch, only: landpatch
      USE MOD_NetCDFVector, only: ncio_write_vector
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart, varname
      real(r8), intent(in) :: value
      integer, intent(in) :: compress
      real(r8), allocatable :: marker(:)

      IF (p_is_worker) THEN
         allocate(marker(landpatch%nset))
         marker(:) = value
      ELSE
         allocate(marker(0))
      ENDIF
      CALL ncio_write_vector(file_restart, varname, 'patch', landpatch, marker, compress)
      deallocate(marker)

   END SUBROUTINE write_methane_restart_marker

   SUBROUTINE validate_methane_restart_transaction (file_restart, strict_restart, restart_schema, &
         history_mode_changed, checkpoint_has_microbe_pools)

      USE MOD_LandPatch, only: landpatch
      USE MOD_MT_NcIO, only: ncio_read_vector_complete, ncio_vector_group_presence
#ifdef USEMPI
      USE MOD_SPMD_Task, only: p_is_worker, p_is_master, p_comm_glb, p_err, &
         MPI_IN_PLACE, MPI_LOGICAL, MPI_LOR
#else
      USE MOD_SPMD_Task, only: p_is_worker, p_is_master
#endif
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      character(len=40), parameter :: restart_metadata_fields(4) = &
         [character(len=40) :: 'ch4_restart_schema', 'ch4_restart_complete', &
          'ch4_history_accumulation_mode', 'ch4_feature_microbe_pools']
      logical, intent(out) :: strict_restart
      integer, intent(out) :: restart_schema
      logical, intent(out) :: history_mode_changed, checkpoint_has_microbe_pools
      real(r8), allocatable :: schema(:), committed(:), history_mode(:), microbe_feature(:)
      logical :: has_schema, has_commit, has_history_mode, has_microbe_feature
      logical :: invalid_marker, invalid_microbe_feature
      logical :: microbe_feature_has_zero, microbe_feature_has_one
      logical :: restart_metadata_present(4)
      integer :: local_schema

      strict_restart = .false.
      restart_schema = 0
      history_mode_changed = .false.
      checkpoint_has_microbe_pools = .false.
      ! Only complete absence of transaction metadata identifies a legacy
      ! restart.  Orphaned optional metadata or one-sided schema/commit markers
      ! identify an interrupted/corrupt write and must never fall back.
      CALL ncio_vector_group_presence(file_restart, restart_metadata_fields, &
         landpatch, restart_metadata_present)
      has_schema = restart_metadata_present(1)
      has_commit = restart_metadata_present(2)
      has_history_mode = restart_metadata_present(3)
      has_microbe_feature = restart_metadata_present(4)
      IF (.not. has_schema .and. .not. has_commit) THEN
         IF (has_history_mode .or. has_microbe_feature) THEN
            IF (p_is_master) WRITE(*,'(A)') &
               'ERROR: orphaned methane restart transaction metadata; refusing legacy fallback.'
            CALL CoLM_stop()
         ENDIF
         RETURN
      ENDIF
      IF (.not. has_schema .or. .not. has_commit) THEN
         IF (p_is_master) WRITE(*,'(A)') &
            'ERROR: incomplete methane restart transaction metadata; refusing mixed checkpoint state.'
         CALL CoLM_stop()
      ENDIF

      CALL ncio_read_vector_complete(file_restart, 'ch4_restart_schema', landpatch, schema)
      CALL ncio_read_vector_complete(file_restart, 'ch4_restart_complete', landpatch, committed)
      local_schema = 0
      IF (p_is_worker .and. allocated(schema)) THEN
         IF (any(schema == real(METHANE_RESTART_SCHEMA_VERSION, r8))) THEN
            local_schema = METHANE_RESTART_SCHEMA_VERSION
         ELSEIF (any(schema == 4._r8)) THEN
            local_schema = 4
         ELSEIF (any(schema == 3._r8)) THEN
            local_schema = 3
         ELSEIF (any(schema == 2._r8)) THEN
            local_schema = 2
         ELSEIF (any(schema == 1._r8)) THEN
            local_schema = 1
         ENDIF
      ENDIF
      restart_schema = local_schema
#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, restart_schema, 1, MPI_INTEGER, MPI_MAX, p_comm_glb, p_err)
#endif
      invalid_marker = .false.
      IF (p_is_worker) THEN
         IF (allocated(schema)) invalid_marker = invalid_marker .or. &
            any(schema /= real(restart_schema, r8))
         IF (allocated(committed)) invalid_marker = invalid_marker .or. any(committed /= 1._r8)
      ENDIF
      invalid_marker = invalid_marker .or. &
         (restart_schema /= 1 .and. restart_schema /= 2 .and. restart_schema /= 3 .and. &
          restart_schema /= 4 .and. &
          restart_schema /= METHANE_RESTART_SCHEMA_VERSION)
#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, invalid_marker, 1, MPI_LOGICAL, MPI_LOR, p_comm_glb, p_err)
#endif
      IF (invalid_marker) THEN
         IF (p_is_master) WRITE(*,'(A)') &
            'ERROR: methane restart transaction is uncommitted or has an unsupported schema.'
         CALL CoLM_stop()
      ENDIF

      ! The coarse mode catches off/core/custom changes.  AccFlux separately
      ! validates the exact custom selector fingerprint because custom lists
      ! can accumulate different array families.
      IF (restart_schema >= 3 .and. .not. has_history_mode) THEN
         IF (p_is_master) WRITE(*,'(A)') &
            'ERROR: schema-3+ methane restart is missing the history accumulation mode.'
         CALL CoLM_stop()
      ENDIF
      history_mode_changed = .not. has_history_mode
      IF (has_history_mode) THEN
         CALL ncio_read_vector_complete(file_restart, 'ch4_history_accumulation_mode', &
            landpatch, history_mode)
         IF (p_is_worker .and. allocated(history_mode)) THEN
            history_mode_changed = any(history_mode /= &
               real(methane_history_accumulation_mode(), r8))
         ENDIF
#ifdef USEMPI
         CALL mpi_allreduce(MPI_IN_PLACE, history_mode_changed, 1, MPI_LOGICAL, &
            MPI_LOR, p_comm_glb, p_err)
#endif
         IF (allocated(history_mode)) deallocate(history_mode)
      ENDIF

      ! Schema 5 owns the microbial feature group through a transaction-local
      ! marker.  It must be complete, binary and identical on every patch so a
      ! reused target cannot resurrect stale ON-state after an OFF checkpoint.
      IF (restart_schema >= 5 .and. .not. has_microbe_feature) THEN
         IF (p_is_master) WRITE(*,'(A)') &
            'ERROR: schema-5 methane restart is missing the microbial-pool feature marker.'
         CALL CoLM_stop()
      ENDIF
      IF (restart_schema >= 5) THEN
         CALL ncio_read_vector_complete(file_restart, 'ch4_feature_microbe_pools', &
            landpatch, microbe_feature)
         invalid_microbe_feature = .false.
         microbe_feature_has_zero = .false.
         microbe_feature_has_one = .false.
         IF (p_is_worker .and. allocated(microbe_feature)) THEN
            invalid_microbe_feature = any((microbe_feature /= 0._r8) .and. &
               (microbe_feature /= 1._r8))
            microbe_feature_has_zero = any(microbe_feature == 0._r8)
            microbe_feature_has_one = any(microbe_feature == 1._r8)
         ENDIF
#ifdef USEMPI
         CALL mpi_allreduce(MPI_IN_PLACE, invalid_microbe_feature, 1, MPI_LOGICAL, &
            MPI_LOR, p_comm_glb, p_err)
         CALL mpi_allreduce(MPI_IN_PLACE, microbe_feature_has_zero, 1, MPI_LOGICAL, &
            MPI_LOR, p_comm_glb, p_err)
         CALL mpi_allreduce(MPI_IN_PLACE, microbe_feature_has_one, 1, MPI_LOGICAL, &
            MPI_LOR, p_comm_glb, p_err)
#endif
         invalid_microbe_feature = invalid_microbe_feature .or. &
            (microbe_feature_has_zero .eqv. microbe_feature_has_one)
         IF (invalid_microbe_feature) THEN
            IF (p_is_master) WRITE(*,'(A)') &
               'ERROR: schema-5 methane restart has an invalid or nonuniform microbial-pool feature marker.'
            CALL CoLM_stop()
         ENDIF
         checkpoint_has_microbe_pools = microbe_feature_has_one
         IF (allocated(microbe_feature)) deallocate(microbe_feature)
      ENDIF
      IF (allocated(schema)) deallocate(schema)
      IF (allocated(committed)) deallocate(committed)
      strict_restart = .true.

   END SUBROUTINE validate_methane_restart_transaction

   SUBROUTINE ch4_reactive_lake_step (istep_local, ipatch, idate, deltim_phy, isub, nsub)

      IMPLICIT NONE
      integer,  intent(in) :: istep_local
      integer,  intent(in) :: ipatch
      integer,  intent(in) :: idate(3)
      real(r8), intent(in) :: deltim_phy
      integer,  intent(in) :: isub
      integer,  intent(in) :: nsub

      IF (ch4_reactive_has()) THEN
         CALL ch4_impl_lake_step (istep_local, ipatch, idate, deltim_phy, isub, nsub)
      ENDIF

   END SUBROUTINE ch4_reactive_lake_step

   SUBROUTINE ch4_reactive_wetland_decomp (ipatch, deltim)

      IMPLICIT NONE
      integer, intent(in) :: ipatch
      real(r8), intent(in) :: deltim

      IF (ch4_reactive_has()) CALL ch4_impl_wetland_decomp (ipatch, deltim)

   END SUBROUTINE ch4_reactive_wetland_decomp

   SUBROUTINE ch4_reactive_soil_step (istep_local, ipatch, idate, deltim)

      IMPLICIT NONE
      integer,  intent(in) :: istep_local
      integer,  intent(in) :: ipatch
      integer,  intent(in) :: idate(3)
      real(r8), intent(in) :: deltim

      IF (ch4_reactive_has()) CALL ch4_impl_soil_step (istep_local, ipatch, idate, deltim)

   END SUBROUTINE ch4_reactive_soil_step

   SUBROUTINE ch4_reactive_report ()

      IMPLICIT NONE

      ! Generic material conservation is reported unconditionally from
      ! MOD_MT_LandPhase::mt_report so non-CH4 material runs are covered too.

      ! The host-water residual reduce stays on this per-timestep path on
      ! purpose. mt_report is reached only under CoLM.F90's IF (p_is_worker),
      ! so p_comm_worker is exactly the set of ranks that get here -- the one
      ! property that matters for a collective. The enclosing material_balance_report
      ! already performs four p_comm_worker collectives per call, so one more
      ! reduce of a single double adds no new synchronisation point. Moving it to
      ! ch4_reactive_final was tried and reverted: land_final is reachable from
      ! two call paths whose rank membership could not be established, and a
      ! collective on an unproven rank set is how the earlier deadlock happened.
      CALL methane_host_water_report ()

   END SUBROUTINE ch4_reactive_report

   SUBROUTINE ch4_reactive_flush_acc_fluxes ()

      IMPLICIT NONE

      IF (ch4_reactive_has()) CALL flush_methane_acc_fluxes ()

   END SUBROUTINE ch4_reactive_flush_acc_fluxes

   SUBROUTINE ch4_reactive_accumulate_fluxes ()

      IMPLICIT NONE

      IF (ch4_reactive_has()) CALL accumulate_methane_fluxes ()

   END SUBROUTINE ch4_reactive_accumulate_fluxes

   SUBROUTINE ch4_reactive_save_lulcc_state ()

      IMPLICIT NONE

      IF (.not. ch4_reactive_has()) RETURN
      CALL save_methane_lulcc_state ()
      IF (DEF_METHANE%use_microbial_pools) CALL save_methane_microbes_lulcc_state ()
      CALL methane_wetland_bgc_save_lulcc ()

   END SUBROUTINE ch4_reactive_save_lulcc_state

   SUBROUTINE ch4_reactive_remap_lulcc_state (patchclass_new, eindex_new, patchclass_old, eindex_old, &
      lccpct_patches, new_patch_area, old_patch_area)

      IMPLICIT NONE
      integer, intent(in) :: patchclass_new(:), patchclass_old(:)
      integer*8, intent(in) :: eindex_new(:), eindex_old(:)
      real(r8), intent(in), optional :: lccpct_patches(:,0:)
      real(r8), intent(in), optional :: new_patch_area(:)
      real(r8), intent(in), optional :: old_patch_area(:)
      integer :: nnew

      IF (.not. ch4_reactive_has()) RETURN
      CALL remap_methane_lulcc_state (patchclass_new, eindex_new, patchclass_old, eindex_old, &
         lccpct_patches, new_patch_area, old_patch_area)
      IF (DEF_METHANE%use_microbial_pools) THEN
         CALL remap_methane_microbes_lulcc_state (patchclass_new, eindex_new, patchclass_old, eindex_old, &
            lccpct_patches, new_patch_area, old_patch_area)
      ENDIF

      ! LULCC can change the worker-local patch count.  Methane state and
      ! microbial pools are remapped above; all other patch-sized methane
      ! buffers must be resized in the same transaction.  Leaving the old
      ! accumulator size would make acc1d(state, accumulator) write past the
      ! accumulator when numpatch grows, or keep stale/mis-aligned entries
      ! when it shrinks.
      nnew = size(patchclass_new)
      IF (allocated(patchtype)) THEN
         CALL methane_wetland_bgc_remap_lulcc (patchtype, patchclass_new, eindex_new, &
            patchclass_old, eindex_old)
      ENDIF
      CALL flush_methane_acc_fluxes ()
      CALL deallocate_methane_acc_fluxes ()
      CALL allocate_methane_acc_fluxes (nnew)

      CALL deallocate_wetland_aere_overrides ()
      CALL allocate_wetland_aere_overrides (nnew)

   END SUBROUTINE ch4_reactive_remap_lulcc_state

   SUBROUTINE ch4_reactive_reload_lulcc_inputs (lc_year, dir_landdata)

      IMPLICIT NONE
      integer, intent(in) :: lc_year
      character(len=*), intent(in) :: dir_landdata
      integer :: nnew
      character(len=4) :: cyear_lulcc
      real(r8), allocatable :: giems_dummy_patch(:)

      IF (.not. ch4_reactive_has()) RETURN
      IF (DEF_METHANE%use_spatial_ph .and. &
          (lc_year <= 0 .or. len_trim(dir_landdata) == 0)) THEN
         CALL CoLM_stop (' ***** ERROR: methane LULCC pH reload requires current landdata directory and year.')
      ENDIF

      ! GIEMS and spatial-pH loading are collective.  Derive the worker-local
      ! size here so non-worker ranks participate with zero-length vectors.
      nnew = 0
      IF (p_is_worker .and. allocated(patchtype)) nnew = size(patchtype)

      CALL deallocate_methane_giems ()
      CALL allocate_methane_giems (nnew)
      IF (DEF_wetland_finundation_scheme == 5) THEN
         IF (trim(DEF_file_GIEMS) == 'null') THEN
            CALL CoLM_stop (' ***** ERROR: DEF_wetland_finundation_scheme=5 requires DEF_file_GIEMS.')
         ENDIF
         CALL require_methane_giems_patch_coords (nnew)
         IF (p_is_worker .and. nnew > 0) THEN
            CALL read_methane_giems (DEF_file_GIEMS, patchlatr, patchlonr, nnew)
         ELSE
            allocate(giems_dummy_patch(0))
            CALL read_methane_giems (DEF_file_GIEMS, giems_dummy_patch, giems_dummy_patch, 0)
            deallocate(giems_dummy_patch)
         ENDIF
         IF (.not. giems_active) THEN
            CALL CoLM_stop (' ***** ERROR: GIEMS file could not be loaded for methane scheme 5.')
         ENDIF
      ENDIF

      CALL load_methane_lake_soilc_srf (lc_year, dir_landdata, nnew)

      CALL deallocate_methane_ph ()
      CALL allocate_methane_ph (nnew)
      IF (DEF_METHANE%use_spatial_ph) THEN
         write(cyear_lulcc,'(i4.4)') lc_year
         last_methane_ph_patch_file = trim(dir_landdata)//'/soil/'//trim(cyear_lulcc)// &
            '/methane_ph_patches.nc'
         CALL read_methane_ph_patch (trim(last_methane_ph_patch_file), nnew)
      ENDIF

   END SUBROUTINE ch4_reactive_reload_lulcc_inputs

   SUBROUTINE require_methane_giems_patch_coords (local_numpatch)

      IMPLICIT NONE
      integer, intent(in) :: local_numpatch
      logical :: invalid_coords

      invalid_coords = .false.
      IF (p_is_worker .and. local_numpatch > 0) THEN
         IF (.not. allocated(patchlatr) .or. .not. allocated(patchlonr)) THEN
            invalid_coords = .true.
         ELSE
            invalid_coords = size(patchlatr) < local_numpatch .or. &
                             size(patchlonr) < local_numpatch
         ENDIF
      ENDIF
#ifdef USEMPI
      CALL MPI_Allreduce (MPI_IN_PLACE, invalid_coords, 1, MPI_LOGICAL, MPI_LOR, p_comm_glb, p_err)
#endif
      IF (invalid_coords) THEN
         CALL CoLM_stop (' ***** ERROR: GIEMS methane input requires coordinates for every worker-local patch.')
      ENDIF

   END SUBROUTINE require_methane_giems_patch_coords

   SUBROUTINE ch4_reactive_publish_levee_flood (fldfrc_patch)

      IMPLICIT NONE
      real(r8), intent(in) :: fldfrc_patch(:)

      IF (.not. ch4_reactive_has()) RETURN
      CALL publish_methane_levee_flood_patch (fldfrc_patch)

   END SUBROUTINE ch4_reactive_publish_levee_flood

   SUBROUTINE ch4_reactive_publish_flood (fldfrc_patch, flddph_patch)

      IMPLICIT NONE
      real(r8), intent(in) :: fldfrc_patch(:)
      real(r8), intent(in) :: flddph_patch(:)

      IF (.not. ch4_reactive_has()) RETURN
      CALL publish_methane_flood_patch (fldfrc_patch, flddph_patch)

   END SUBROUTINE ch4_reactive_publish_flood

   SUBROUTINE ch4_reactive_final ()

      IMPLICIT NONE

      last_methane_ph_patch_file = ''

      CALL deallocate_methane_acc_fluxes ()
      CALL deallocate_wetland_aere_overrides ()
      CALL deallocate_methane_ph ()
      CALL deallocate_methane_giems ()
      CALL deallocate_methane_microbes_state ()
      CALL deallocate_methane_state ()
      CALL methane_wetland_bgc_deallocate ()

   END SUBROUTINE ch4_reactive_final


END MODULE MOD_MT_Methane
#endif
