#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_Lifecycle
!=======================================================================
! One lifecycle row per configured material. Provider names are resolved
! once during initialization; phase dispatch thereafter uses material index.
!=======================================================================

   USE MOD_Precision
   USE MOD_DataType, only: block_data_real8_2d
   USE MOD_SPMD_Task, only: CoLM_stop, p_is_master
   USE MOD_MT_Defs, only: nmts, materials, material_upper, material_index_for_name, &
      FAMILY_UNRESOLVED, FAMILY_ISOTOPE, FAMILY_SOLUTE, &
      FAMILY_PARTICLE, FAMILY_GAS, STATE_OWNER_UNKNOWN, &
      STATE_OWNER_GENERIC_WATER, STATE_OWNER_PROVIDER, &
      REACTION_NONE, REACTION_FIRST_ORDER, REACTION_PROVIDER

   IMPLICIT NONE
   PRIVATE

   abstract interface
      SUBROUTINE lifecycle_noarg_if ()
      END SUBROUTINE lifecycle_noarg_if

      SUBROUTINE lifecycle_land_init_if (numpatch, lc_year, jdate, casename, dir_restart, dir_landdata, file_restart)
         integer, intent(in) :: numpatch, lc_year, jdate(3)
         character(len=*), intent(in) :: casename, dir_restart, dir_landdata
         character(len=*), intent(in), optional :: file_restart
      END SUBROUTINE lifecycle_land_init_if

      SUBROUTINE lifecycle_read_restart_if (file_restart)
         character(len=*), intent(in) :: file_restart
      END SUBROUTINE lifecycle_read_restart_if

      SUBROUTINE lifecycle_land_write_restart_if (file_restart, compress)
         character(len=*), intent(in) :: file_restart
         integer, intent(in) :: compress
      END SUBROUTINE lifecycle_land_write_restart_if

      SUBROUTINE lifecycle_route_write_restart_if (file_restart)
         character(len=*), intent(in) :: file_restart
      END SUBROUTINE lifecycle_route_write_restart_if

      SUBROUTINE lifecycle_lake_step_if (istep_local, ipatch, idate, deltim_phy, isub, nsub)
         USE MOD_Precision
         integer, intent(in) :: istep_local, ipatch, idate(3), isub, nsub
         real(r8), intent(in) :: deltim_phy
      END SUBROUTINE lifecycle_lake_step_if

      SUBROUTINE lifecycle_wetland_decomp_if (ipatch, deltim)
         USE MOD_Precision
         integer, intent(in) :: ipatch
         real(r8), intent(in) :: deltim
      END SUBROUTINE lifecycle_wetland_decomp_if

      SUBROUTINE lifecycle_soil_step_if (istep_local, ipatch, idate, deltim)
         USE MOD_Precision
         integer, intent(in) :: istep_local, ipatch, idate(3)
         real(r8), intent(in) :: deltim
      END SUBROUTINE lifecycle_soil_step_if

      SUBROUTINE lifecycle_land_history_if (file_hist, itime_in_file, sumarea, filter, &
         nl_soil, forcing_has_missing_value, forcmask_pch)
         USE MOD_DataType, only: block_data_real8_2d
         character(len=*), intent(in) :: file_hist
         integer, intent(in) :: itime_in_file, nl_soil
         type(block_data_real8_2d), intent(inout) :: sumarea
         logical, intent(inout) :: filter(:)
         logical, intent(in) :: forcing_has_missing_value, forcmask_pch(:)
      END SUBROUTINE lifecycle_land_history_if

      SUBROUTINE lifecycle_land_remap_lulcc_if (patchclass_new, eindex_new, patchclass_old, eindex_old, &
         lccpct_patches, new_patch_area, old_patch_area)
         USE MOD_Precision
         integer, intent(in) :: patchclass_new(:), patchclass_old(:)
         integer*8, intent(in) :: eindex_new(:), eindex_old(:)
         real(r8), intent(in), optional :: lccpct_patches(:,0:), new_patch_area(:), old_patch_area(:)
      END SUBROUTINE lifecycle_land_remap_lulcc_if

      SUBROUTINE lifecycle_land_reload_lulcc_if (lc_year, dir_landdata)
         integer, intent(in) :: lc_year
         character(len=*), intent(in) :: dir_landdata
      END SUBROUTINE lifecycle_land_reload_lulcc_if

      SUBROUTINE lifecycle_publish_levee_flood_if (fldfrc_patch)
         USE MOD_Precision
         real(r8), intent(in) :: fldfrc_patch(:)
      END SUBROUTINE lifecycle_publish_levee_flood_if

      SUBROUTINE lifecycle_publish_flood_if (fldfrc_patch, flddph_patch)
         USE MOD_Precision
         real(r8), intent(in) :: fldfrc_patch(:), flddph_patch(:)
      END SUBROUTINE lifecycle_publish_flood_if

      SUBROUTINE lifecycle_route_forcing_put_if (precip, dt, valid_fraction)
         USE MOD_Precision
         real(r8), intent(in) :: precip(:), dt
         real(r8), optional, intent(in) :: valid_fraction(:)
      END SUBROUTINE lifecycle_route_forcing_put_if

      SUBROUTINE lifecycle_route_diag_if (dt_all, irivsys, ucatfilter, veloc, wdsrf, &
         rivsto_start, rivsto_input, rivout_fc, floodarea, protected_start, protected_end, protected_area)
         USE MOD_Precision
         real(r8), intent(in) :: dt_all(:), veloc(:), wdsrf(:)
         real(r8), intent(in) :: rivsto_start(:), rivsto_input(:), rivout_fc(:), floodarea(:)
         real(r8), optional, intent(in) :: protected_start(:), protected_end(:), protected_area(:)
         integer, intent(in) :: irivsys(:)
         logical, intent(in) :: ucatfilter(:)
      END SUBROUTINE lifecycle_route_diag_if

      SUBROUTINE lifecycle_route_calc_if (deltime)
         USE MOD_Precision
         real(r8), intent(in) :: deltime
      END SUBROUTINE lifecycle_route_calc_if

      SUBROUTINE lifecycle_route_sediment_bif_if (dt_all, irivsys, ucatfilter, bif_hflux_lev)
         USE MOD_Precision
         real(r8), intent(in) :: dt_all(:), bif_hflux_lev(:,:)
         integer, intent(in) :: irivsys(:)
         logical, intent(in) :: ucatfilter(:)
      END SUBROUTINE lifecycle_route_sediment_bif_if

      SUBROUTINE lifecycle_route_sediment_levee_if (i, visible_before, protected_before, &
         visible_after, protected_after)
         USE MOD_Precision
         integer, intent(in) :: i
         real(r8), intent(in) :: visible_before, protected_before, visible_after, protected_after
      END SUBROUTINE lifecycle_route_sediment_levee_if

      SUBROUTINE lifecycle_route_history_if (file_hist_ucat, itime_in_file_ucat)
         character(len=*), intent(in) :: file_hist_ucat
         integer, intent(in) :: itime_in_file_ucat
      END SUBROUTINE lifecycle_route_history_if
   end interface

   type, PUBLIC :: material_lifecycle_hooks_type
      procedure(lifecycle_land_init_if), pointer, nopass :: land_init => null()
      procedure(lifecycle_read_restart_if), pointer, nopass :: land_read_restart => null()
      procedure(lifecycle_land_write_restart_if), pointer, nopass :: land_write_restart => null()
      procedure(lifecycle_lake_step_if), pointer, nopass :: lake_step => null()
      procedure(lifecycle_wetland_decomp_if), pointer, nopass :: wetland_decomp => null()
      procedure(lifecycle_soil_step_if), pointer, nopass :: soil_step => null()
      procedure(lifecycle_noarg_if), pointer, nopass :: report => null()
      procedure(lifecycle_noarg_if), pointer, nopass :: flush_acc_fluxes => null()
      procedure(lifecycle_noarg_if), pointer, nopass :: accumulate_fluxes => null()
      procedure(lifecycle_land_history_if), pointer, nopass :: land_history => null()
      procedure(lifecycle_read_restart_if), pointer, nopass :: land_history_sidecar_write => null()
      procedure(lifecycle_read_restart_if), pointer, nopass :: land_history_sidecar_read => null()
      procedure(lifecycle_noarg_if), pointer, nopass :: land_save_lulcc => null()
      procedure(lifecycle_land_remap_lulcc_if), pointer, nopass :: land_remap_lulcc => null()
      procedure(lifecycle_land_reload_lulcc_if), pointer, nopass :: land_reload_lulcc => null()
      procedure(lifecycle_publish_levee_flood_if), pointer, nopass :: publish_levee_flood => null()
      procedure(lifecycle_publish_flood_if), pointer, nopass :: publish_flood => null()
      procedure(lifecycle_noarg_if), pointer, nopass :: land_final => null()
      procedure(lifecycle_noarg_if), pointer, nopass :: route_init => null()
      procedure(lifecycle_read_restart_if), pointer, nopass :: route_read_restart => null()
      procedure(lifecycle_route_write_restart_if), pointer, nopass :: route_write_restart => null()
      procedure(lifecycle_route_forcing_put_if), pointer, nopass :: route_forcing_put => null()
      procedure(lifecycle_route_diag_if), pointer, nopass :: route_diag_accumulate => null()
      procedure(lifecycle_route_calc_if), pointer, nopass :: route_calc => null()
      procedure(lifecycle_route_sediment_bif_if), pointer, nopass :: route_sediment_bif => null()
      procedure(lifecycle_route_sediment_levee_if), pointer, nopass :: route_sediment_levee => null()
      procedure(lifecycle_route_history_if), pointer, nopass :: route_history => null()
      procedure(lifecycle_noarg_if), pointer, nopass :: route_flush_history => null()
      procedure(lifecycle_noarg_if), pointer, nopass :: route_final => null()
   end type material_lifecycle_hooks_type

   type(material_lifecycle_hooks_type), allocatable, save :: lifecycle(:)
   character(len=32), allocatable, save :: provider_name(:)

   PUBLIC :: material_lifecycle_init, material_lifecycle_validate, mt_lifecycle_reset
   PUBLIC :: material_lifecycle_report_registry
   PUBLIC :: register_mt_provider
   PUBLIC :: material_lifecycle_size, material_lifecycle_has_provider, material_lifecycle_provider
   PUBLIC :: material_lifecycle_land_init, material_lifecycle_land_final
   PUBLIC :: material_lifecycle_lake_step, material_lifecycle_wetland_decomp, material_lifecycle_soil_step
   PUBLIC :: material_lifecycle_land_report, material_lifecycle_land_read_restart, material_lifecycle_land_write_restart
   PUBLIC :: material_lifecycle_land_history, material_lifecycle_land_flush_acc_fluxes
   PUBLIC :: material_lifecycle_land_history_sidecar
   PUBLIC :: material_lifecycle_land_accumulate_fluxes, mt_lifecycle_land_save_lulcc_state
   PUBLIC :: mt_lifecycle_land_remap_lulcc_state, mt_lifecycle_land_reload_lulcc_inputs
   PUBLIC :: material_lifecycle_has_levee_flood_publisher, material_lifecycle_has_flood_publisher
   PUBLIC :: mt_lifecycle_publish_levee_flood_patch, mt_lifecycle_publish_flood_patch
   PUBLIC :: mt_lifecycle_route_has_active, mt_lifecycle_route_init, mt_lifecycle_route_final
   PUBLIC :: mt_lifecycle_route_read_restart, mt_lifecycle_route_write_restart
   PUBLIC :: mt_lifecycle_route_forcing_put, mt_lifecycle_route_diag_accumulate, mt_lifecycle_route_calc
   PUBLIC :: mt_lifecycle_route_sediment_bif_accumulate, mt_lifecycle_route_sediment_levee_repartition
   PUBLIC :: mt_lifecycle_route_write_history, mt_lifecycle_route_flush_history

CONTAINS

   ! register_providers is mt_register_providers (MT_RegisterProviders.F90),
   ! passed in by the caller so that this module, linked into mkinidata too,
   ! does not reference the providers built with the main objects.
   SUBROUTINE material_lifecycle_init (register_providers)
      procedure(lifecycle_noarg_if) :: register_providers
      IF (allocated(lifecycle)) RETURN
      IF (nmts > 0 .and. .not. allocated(materials)) &
         CALL lifecycle_error ('descriptor table is not initialized')
      allocate(lifecycle(nmts), provider_name(nmts))
      provider_name = ''
      CALL register_providers ()
      CALL material_lifecycle_validate ()
      CALL material_lifecycle_report_registry ()
   END SUBROUTINE material_lifecycle_init

   SUBROUTINE mt_lifecycle_reset ()
      IF (allocated(lifecycle)) deallocate(lifecycle)
      IF (allocated(provider_name)) deallocate(provider_name)
   END SUBROUTINE mt_lifecycle_reset

   SUBROUTINE register_mt_provider (name, aliases, provider, declared_family, &
      declared_owner, declared_reaction, hooks, imt)
      character(len=*), intent(in) :: name, aliases, provider
      integer, intent(in) :: declared_family, declared_owner, declared_reaction
      type(material_lifecycle_hooks_type), intent(in) :: hooks
      integer, intent(out) :: imt

      integer :: i

      imt = 0
      IF (.not. allocated(lifecycle) .or. .not. allocated(provider_name)) &
         CALL lifecycle_error ('provider registration requires material_lifecycle_init')
      CALL validate_provider_declaration (name, aliases, provider, declared_family, &
         declared_owner, declared_reaction)
      IF (.not. lifecycle_hooks_present(hooks)) &
         CALL lifecycle_error ('compiled provider has no lifecycle hooks')

      imt = material_index_for_name(name, aliases)
      IF (imt == 0) RETURN
      IF (len_trim(provider_name(imt)) > 0) &
         CALL lifecycle_error ('duplicate provider attachment for configured material')

      DO i = 1, nmts
         IF (len_trim(provider_name(i)) > 0 .and. &
             trim(material_upper(provider_name(i))) == trim(material_upper(provider))) THEN
            CALL lifecycle_error ('duplicate provider identity')
         ENDIF
      ENDDO

      IF (materials(imt)%family_id == FAMILY_UNRESOLVED) THEN
         materials(imt)%family_id = declared_family
      ELSEIF (materials(imt)%family_id /= declared_family) THEN
         CALL lifecycle_error ('provider family does not match configured material family')
      ENDIF

      materials(imt)%state_owner = declared_owner
      materials(imt)%reaction_mode = declared_reaction
      materials(imt)%category = family_label(declared_family)
      IF (declared_owner == STATE_OWNER_PROVIDER .and. declared_family /= FAMILY_PARTICLE) &
         materials(imt)%unit_kind = 'species_owned'
      IF (declared_owner == STATE_OWNER_GENERIC_WATER .and. &
          trim(materials(imt)%unit_kind) == 'species_owned') &
         materials(imt)%unit_kind = 'material_per_water'

      provider_name(imt) = trim(provider)
      lifecycle(imt) = hooks
   END SUBROUTINE register_mt_provider

   SUBROUTINE material_lifecycle_validate ()
      integer :: i

      IF (.not. allocated(lifecycle) .or. .not. allocated(provider_name)) &
         CALL lifecycle_error ('lifecycle registry is not initialized')
      IF (size(lifecycle) /= nmts .or. size(provider_name) /= nmts) &
         CALL lifecycle_error ('lifecycle registry size does not match descriptor table')
      IF (nmts > 0 .and. .not. allocated(materials)) &
         CALL lifecycle_error ('descriptor table is not initialized')

      DO i = 1, nmts
         IF ((trim(material_upper(materials(i)%name)) == 'CH4' .or. &
              trim(material_upper(materials(i)%name)) == 'METHANE') .and. &
             len_trim(provider_name(i)) == 0) THEN
            CALL lifecycle_error ('CH4 requires its compiled gas provider')
         ENDIF

         IF (materials(i)%state_owner == STATE_OWNER_UNKNOWN) &
            CALL lifecycle_error ('material state owner remains unresolved')
         IF ((materials(i)%family_id == FAMILY_GAS .or. &
              materials(i)%family_id == FAMILY_PARTICLE .or. &
              materials(i)%state_owner == STATE_OWNER_PROVIDER .or. &
              materials(i)%reaction_mode == REACTION_PROVIDER) .and. &
             len_trim(provider_name(i)) == 0) THEN
            CALL lifecycle_error ('provider-owned material has no compiled provider')
         ENDIF
         IF ((materials(i)%family_id == FAMILY_GAS .or. &
              materials(i)%family_id == FAMILY_PARTICLE) .and. &
             materials(i)%state_owner /= STATE_OWNER_PROVIDER) THEN
            CALL lifecycle_error ('gas and particle families require provider-owned state')
         ENDIF
         IF ((materials(i)%family_id == FAMILY_ISOTOPE .or. &
              materials(i)%family_id == FAMILY_PARTICLE) .and. &
             materials(i)%reaction_mode /= REACTION_NONE) &
            CALL lifecycle_error ('isotope and particle reaction shape is not applicable')
         IF (materials(i)%reaction_mode == REACTION_FIRST_ORDER .and. &
             (materials(i)%family_id /= FAMILY_SOLUTE .or. &
              materials(i)%state_owner /= STATE_OWNER_GENERIC_WATER)) &
            CALL lifecycle_error ('first-order reaction requires generic-water solute state')
         IF (materials(i)%family_id == FAMILY_GAS .and. &
             .not. lifecycle_gas_hooks_complete(lifecycle(i))) &
            CALL lifecycle_error ('gas provider is missing required land lifecycle hooks')
         IF (materials(i)%family_id == FAMILY_PARTICLE .and. &
             .not. lifecycle_particle_hooks_complete(lifecycle(i))) &
            CALL lifecycle_error ('particle provider is missing required route lifecycle hooks')
         IF (materials(i)%reaction_mode == REACTION_PROVIDER .and. &
             .not. lifecycle_reaction_hooks_present(lifecycle(i))) &
            CALL lifecycle_error ('reaction provider is missing a reaction-step hook')
      ENDDO
   END SUBROUTINE material_lifecycle_validate

   SUBROUTINE material_lifecycle_report_registry ()
      integer :: i
      character(len=32) :: provider, hook_groups

      IF (.not. p_is_master) RETURN
      IF (.not. allocated(lifecycle) .or. .not. allocated(provider_name)) RETURN

      WRITE (*,'(A,I0)') 'MT lifecycle registry: count=', nmts
      DO i = 1, nmts
         provider = '-'
         IF (len_trim(provider_name(i)) > 0) provider = trim(provider_name(i))
         hook_groups = lifecycle_hook_groups(lifecycle(i))
         WRITE (*,'(A,I0,A,A,A,A,A,A,A,A,A,A,A,A)') &
            '  [', i, '] name=', trim(materials(i)%name), &
            ' family=', trim(family_label(materials(i)%family_id)), &
            ' owner=', trim(owner_label(materials(i)%state_owner)), &
            ' reaction=', trim(reaction_label(materials(i)%reaction_mode)), &
            ' provider=', trim(provider), ' hooks=', trim(hook_groups)
      ENDDO
   END SUBROUTINE material_lifecycle_report_registry

   integer FUNCTION material_lifecycle_size ()
      material_lifecycle_size = 0
      IF (allocated(lifecycle)) material_lifecycle_size = size(lifecycle)
   END FUNCTION material_lifecycle_size

   logical FUNCTION material_lifecycle_has_provider (imt)
      integer, intent(in) :: imt
      material_lifecycle_has_provider = lifecycle_row_registered(imt)
   END FUNCTION material_lifecycle_has_provider

   FUNCTION material_lifecycle_provider (imt) RESULT(provider)
      integer, intent(in) :: imt
      character(len=32) :: provider
      provider = ''
      IF (lifecycle_row_registered(imt)) provider = provider_name(imt)
   END FUNCTION material_lifecycle_provider

   SUBROUTINE material_lifecycle_land_init (numpatch, lc_year, jdate, casename, dir_restart, dir_landdata, file_restart)
      integer, intent(in) :: numpatch, lc_year, jdate(3)
      character(len=*), intent(in) :: casename, dir_restart, dir_landdata
      character(len=*), intent(in), optional :: file_restart
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%land_init)) &
            CALL lifecycle(i)%land_init(numpatch, lc_year, jdate, casename, dir_restart, dir_landdata, file_restart)
      ENDDO
   END SUBROUTINE material_lifecycle_land_init

   SUBROUTINE material_lifecycle_land_read_restart (file_restart)
      character(len=*), intent(in) :: file_restart
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%land_read_restart)) &
            CALL lifecycle(i)%land_read_restart(file_restart)
      ENDDO
   END SUBROUTINE material_lifecycle_land_read_restart

   SUBROUTINE material_lifecycle_land_write_restart (file_restart, compress)
      character(len=*), intent(in) :: file_restart
      integer, intent(in) :: compress
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%land_write_restart)) &
            CALL lifecycle(i)%land_write_restart(file_restart, compress)
      ENDDO
   END SUBROUTINE material_lifecycle_land_write_restart

   SUBROUTINE material_lifecycle_lake_step (istep_local, ipatch, idate, deltim_phy, isub, nsub)
      integer, intent(in) :: istep_local, ipatch, idate(3), isub, nsub
      real(r8), intent(in) :: deltim_phy
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%lake_step)) &
            CALL lifecycle(i)%lake_step(istep_local, ipatch, idate, deltim_phy, isub, nsub)
      ENDDO
   END SUBROUTINE material_lifecycle_lake_step

   SUBROUTINE material_lifecycle_wetland_decomp (ipatch, deltim)
      integer, intent(in) :: ipatch
      real(r8), intent(in) :: deltim
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%wetland_decomp)) &
            CALL lifecycle(i)%wetland_decomp(ipatch, deltim)
      ENDDO
   END SUBROUTINE material_lifecycle_wetland_decomp

   SUBROUTINE material_lifecycle_soil_step (istep_local, ipatch, idate, deltim)
      integer, intent(in) :: istep_local, ipatch, idate(3)
      real(r8), intent(in) :: deltim
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%soil_step)) &
            CALL lifecycle(i)%soil_step(istep_local, ipatch, idate, deltim)
      ENDDO
   END SUBROUTINE material_lifecycle_soil_step

   SUBROUTINE material_lifecycle_land_report ()
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%report)) CALL lifecycle(i)%report()
      ENDDO
   END SUBROUTINE material_lifecycle_land_report

   SUBROUTINE material_lifecycle_land_flush_acc_fluxes ()
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%flush_acc_fluxes)) &
            CALL lifecycle(i)%flush_acc_fluxes()
      ENDDO
   END SUBROUTINE material_lifecycle_land_flush_acc_fluxes

   SUBROUTINE material_lifecycle_land_history_sidecar (file_hist_acc, writing)
      character(len=*), intent(in) :: file_hist_acc
      logical, intent(in) :: writing
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (.not. lifecycle_row_registered(i)) CYCLE
         IF (writing) THEN
            IF (associated(lifecycle(i)%land_history_sidecar_write)) &
               CALL lifecycle(i)%land_history_sidecar_write(file_hist_acc)
         ELSE
            IF (associated(lifecycle(i)%land_history_sidecar_read)) &
               CALL lifecycle(i)%land_history_sidecar_read(file_hist_acc)
         ENDIF
      ENDDO
   END SUBROUTINE material_lifecycle_land_history_sidecar

   SUBROUTINE material_lifecycle_land_accumulate_fluxes ()
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%accumulate_fluxes)) &
            CALL lifecycle(i)%accumulate_fluxes()
      ENDDO
   END SUBROUTINE material_lifecycle_land_accumulate_fluxes

   SUBROUTINE material_lifecycle_land_history (file_hist, itime_in_file, sumarea, filter, &
      nl_soil, forcing_has_missing_value, forcmask_pch)
      character(len=*), intent(in) :: file_hist
      integer, intent(in) :: itime_in_file, nl_soil
      type(block_data_real8_2d), intent(inout) :: sumarea
      logical, intent(inout) :: filter(:)
      logical, intent(in) :: forcing_has_missing_value, forcmask_pch(:)
      integer :: i, active, sole
      type(block_data_real8_2d) :: callback_sumarea
      logical, allocatable :: callback_filter(:)

      IF (.not. allocated(lifecycle)) RETURN
      active = 0
      sole = 0
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%land_history)) THEN
            active = active + 1
            sole = i
         ENDIF
      ENDDO
      IF (active == 0) RETURN

      ! sumarea/filter are SCRATCH for the callee, not accumulators. They are
      ! intent(inout) because a provider overwrites filter outright and
      ! repeatedly -- methane_reactive_history does
      ! filter = (patchtype == 4) .and. patchmask -- to select the patches each
      ! of its variables applies to.
      !
      ! That is why the multi-provider path below hands each provider a private
      ! copy and never writes it back: writing back would give the next
      ! provider, or the caller, one provider's private mask. It is discarded
      ! on purpose, not by omission.
      !
      ! With one provider there is nobody to be isolated from, so it gets the
      ! caller's arrays directly and the deep copy of sumarea is skipped. The
      ! contract that makes this safe belongs to the caller: filter must be
      ! dead after this dispatch. Today it is -- mt_hist_out ends here and
      ! MOD_Hist deallocates filter on the next line. Anything reading filter
      ! after this call would see one provider's leftovers with a single
      ! provider registered and the original mask with two.
      IF (active == 1) THEN
         CALL lifecycle(sole)%land_history(file_hist, itime_in_file, sumarea, filter, &
            nl_soil, forcing_has_missing_value, forcmask_pch)
         RETURN
      ENDIF
      allocate(callback_filter(size(filter)))
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%land_history)) THEN
            callback_sumarea = sumarea
            callback_filter = filter
            CALL lifecycle(i)%land_history(file_hist, itime_in_file, callback_sumarea, callback_filter, &
               nl_soil, forcing_has_missing_value, forcmask_pch)
         ENDIF
      ENDDO
      deallocate(callback_filter)
   END SUBROUTINE material_lifecycle_land_history

   SUBROUTINE mt_lifecycle_land_save_lulcc_state ()
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%land_save_lulcc)) &
            CALL lifecycle(i)%land_save_lulcc()
      ENDDO
   END SUBROUTINE mt_lifecycle_land_save_lulcc_state

   SUBROUTINE mt_lifecycle_land_remap_lulcc_state (patchclass_new, eindex_new, patchclass_old, eindex_old, &
      lccpct_patches, new_patch_area, old_patch_area)
      integer, intent(in) :: patchclass_new(:), patchclass_old(:)
      integer*8, intent(in) :: eindex_new(:), eindex_old(:)
      real(r8), intent(in), optional :: lccpct_patches(:,0:), new_patch_area(:), old_patch_area(:)
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%land_remap_lulcc)) &
            CALL lifecycle(i)%land_remap_lulcc(patchclass_new, eindex_new, patchclass_old, eindex_old, &
               lccpct_patches, new_patch_area, old_patch_area)
      ENDDO
   END SUBROUTINE mt_lifecycle_land_remap_lulcc_state

   SUBROUTINE mt_lifecycle_land_reload_lulcc_inputs (lc_year, dir_landdata)
      integer, intent(in) :: lc_year
      character(len=*), intent(in) :: dir_landdata
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%land_reload_lulcc)) &
            CALL lifecycle(i)%land_reload_lulcc(lc_year, dir_landdata)
      ENDDO
   END SUBROUTINE mt_lifecycle_land_reload_lulcc_inputs

   logical FUNCTION material_lifecycle_has_levee_flood_publisher ()
      integer :: i
      material_lifecycle_has_levee_flood_publisher = .false.
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%publish_levee_flood)) THEN
            material_lifecycle_has_levee_flood_publisher = .true.
            RETURN
         ENDIF
      ENDDO
   END FUNCTION material_lifecycle_has_levee_flood_publisher

   logical FUNCTION material_lifecycle_has_flood_publisher ()
      integer :: i
      material_lifecycle_has_flood_publisher = .false.
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%publish_flood)) THEN
            material_lifecycle_has_flood_publisher = .true.
            RETURN
         ENDIF
      ENDDO
   END FUNCTION material_lifecycle_has_flood_publisher

   SUBROUTINE mt_lifecycle_publish_levee_flood_patch (fldfrc_patch)
      real(r8), intent(in) :: fldfrc_patch(:)
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%publish_levee_flood)) &
            CALL lifecycle(i)%publish_levee_flood(fldfrc_patch)
      ENDDO
   END SUBROUTINE mt_lifecycle_publish_levee_flood_patch

   SUBROUTINE mt_lifecycle_publish_flood_patch (fldfrc_patch, flddph_patch)
      real(r8), intent(in) :: fldfrc_patch(:), flddph_patch(:)
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%publish_flood)) &
            CALL lifecycle(i)%publish_flood(fldfrc_patch, flddph_patch)
      ENDDO
   END SUBROUTINE mt_lifecycle_publish_flood_patch

   SUBROUTINE material_lifecycle_land_final ()
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%land_final)) &
            CALL lifecycle(i)%land_final()
      ENDDO
   END SUBROUTINE material_lifecycle_land_final

   logical FUNCTION mt_lifecycle_route_has_active ()
      integer :: i
      mt_lifecycle_route_has_active = .false.
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. &
             (associated(lifecycle(i)%route_forcing_put) .or. &
              associated(lifecycle(i)%route_diag_accumulate) .or. &
              associated(lifecycle(i)%route_calc))) THEN
            mt_lifecycle_route_has_active = .true.
            RETURN
         ENDIF
      ENDDO
   END FUNCTION mt_lifecycle_route_has_active

   SUBROUTINE mt_lifecycle_route_init ()
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%route_init)) &
            CALL lifecycle(i)%route_init()
      ENDDO
   END SUBROUTINE mt_lifecycle_route_init

   SUBROUTINE mt_lifecycle_route_read_restart (file_restart)
      character(len=*), intent(in) :: file_restart
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%route_read_restart)) &
            CALL lifecycle(i)%route_read_restart(file_restart)
      ENDDO
   END SUBROUTINE mt_lifecycle_route_read_restart

   SUBROUTINE mt_lifecycle_route_write_restart (file_restart)
      character(len=*), intent(in) :: file_restart
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%route_write_restart)) &
            CALL lifecycle(i)%route_write_restart(file_restart)
      ENDDO
   END SUBROUTINE mt_lifecycle_route_write_restart

   SUBROUTINE mt_lifecycle_route_forcing_put (precip, dt, valid_fraction)
      real(r8), intent(in) :: precip(:), dt
      real(r8), optional, intent(in) :: valid_fraction(:)
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%route_forcing_put)) &
            CALL lifecycle(i)%route_forcing_put(precip, dt, valid_fraction)
      ENDDO
   END SUBROUTINE mt_lifecycle_route_forcing_put

   SUBROUTINE mt_lifecycle_route_diag_accumulate (dt_all, irivsys, ucatfilter, &
      veloc, wdsrf, rivsto_start, rivsto_input, rivout_fc, floodarea, protected_start, protected_end, protected_area)
      real(r8), intent(in) :: dt_all(:), veloc(:), wdsrf(:)
      real(r8), intent(in) :: rivsto_start(:), rivsto_input(:), rivout_fc(:), floodarea(:)
      real(r8), optional, intent(in) :: protected_start(:), protected_end(:), protected_area(:)
      integer, intent(in) :: irivsys(:)
      logical, intent(in) :: ucatfilter(:)
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%route_diag_accumulate)) &
            CALL lifecycle(i)%route_diag_accumulate(dt_all, irivsys, ucatfilter, &
               veloc, wdsrf, rivsto_start, rivsto_input, rivout_fc, floodarea, &
               protected_start, protected_end, protected_area)
      ENDDO
   END SUBROUTINE mt_lifecycle_route_diag_accumulate

   SUBROUTINE mt_lifecycle_route_sediment_bif_accumulate (dt_all, irivsys, ucatfilter, bif_hflux_lev)
      real(r8), intent(in) :: dt_all(:), bif_hflux_lev(:,:)
      integer, intent(in) :: irivsys(:)
      logical, intent(in) :: ucatfilter(:)
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%route_sediment_bif)) &
            CALL lifecycle(i)%route_sediment_bif(dt_all, irivsys, ucatfilter, bif_hflux_lev)
      ENDDO
   END SUBROUTINE mt_lifecycle_route_sediment_bif_accumulate

   SUBROUTINE mt_lifecycle_route_sediment_levee_repartition (i, visible_before, protected_before, &
      visible_after, protected_after)
      integer, intent(in) :: i
      real(r8), intent(in) :: visible_before, protected_before, visible_after, protected_after
      integer :: imt
      IF (.not. allocated(lifecycle)) RETURN
      DO imt = 1, size(lifecycle)
         IF (lifecycle_row_registered(imt) .and. associated(lifecycle(imt)%route_sediment_levee)) &
            CALL lifecycle(imt)%route_sediment_levee(i, visible_before, protected_before, &
               visible_after, protected_after)
      ENDDO
   END SUBROUTINE mt_lifecycle_route_sediment_levee_repartition

   SUBROUTINE mt_lifecycle_route_calc (deltime)
      real(r8), intent(in) :: deltime
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%route_calc)) &
            CALL lifecycle(i)%route_calc(deltime)
      ENDDO
   END SUBROUTINE mt_lifecycle_route_calc

   SUBROUTINE mt_lifecycle_route_write_history (file_hist_ucat, itime_in_file_ucat)
      character(len=*), intent(in) :: file_hist_ucat
      integer, intent(in) :: itime_in_file_ucat
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%route_history)) &
            CALL lifecycle(i)%route_history(file_hist_ucat, itime_in_file_ucat)
      ENDDO
   END SUBROUTINE mt_lifecycle_route_write_history

   SUBROUTINE mt_lifecycle_route_flush_history ()
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%route_flush_history)) &
            CALL lifecycle(i)%route_flush_history()
      ENDDO
   END SUBROUTINE mt_lifecycle_route_flush_history

   SUBROUTINE mt_lifecycle_route_final ()
      integer :: i
      IF (.not. allocated(lifecycle)) RETURN
      DO i = 1, size(lifecycle)
         IF (lifecycle_row_registered(i) .and. associated(lifecycle(i)%route_final)) &
            CALL lifecycle(i)%route_final()
      ENDDO
   END SUBROUTINE mt_lifecycle_route_final

   logical FUNCTION lifecycle_row_registered (imt)
      integer, intent(in) :: imt
      lifecycle_row_registered = .false.
      IF (.not. allocated(provider_name)) RETURN
      IF (imt < 1 .or. imt > size(provider_name)) RETURN
      lifecycle_row_registered = len_trim(provider_name(imt)) > 0
   END FUNCTION lifecycle_row_registered

   logical FUNCTION lifecycle_hooks_present (hooks)
      type(material_lifecycle_hooks_type), intent(in) :: hooks

      lifecycle_hooks_present = associated(hooks%land_init) .or. &
         associated(hooks%land_read_restart) .or. associated(hooks%land_write_restart) .or. &
         associated(hooks%lake_step) .or. associated(hooks%wetland_decomp) .or. &
         associated(hooks%soil_step) .or. associated(hooks%report) .or. &
         associated(hooks%flush_acc_fluxes) .or. associated(hooks%accumulate_fluxes) .or. &
         associated(hooks%land_history) .or. associated(hooks%land_save_lulcc) .or. &
         associated(hooks%land_history_sidecar_write) .or. associated(hooks%land_history_sidecar_read) .or. &
         associated(hooks%land_remap_lulcc) .or. associated(hooks%land_reload_lulcc) .or. &
         associated(hooks%publish_levee_flood) .or. associated(hooks%publish_flood) .or. &
         associated(hooks%land_final) .or. lifecycle_route_hooks_present(hooks)
   END FUNCTION lifecycle_hooks_present

   logical FUNCTION lifecycle_route_hooks_present (hooks)
      type(material_lifecycle_hooks_type), intent(in) :: hooks

      lifecycle_route_hooks_present = associated(hooks%route_init) .or. &
         associated(hooks%route_read_restart) .or. associated(hooks%route_write_restart) .or. &
         associated(hooks%route_forcing_put) .or. associated(hooks%route_diag_accumulate) .or. &
         associated(hooks%route_calc) .or. associated(hooks%route_history) .or. &
         associated(hooks%route_flush_history) .or. associated(hooks%route_final)
   END FUNCTION lifecycle_route_hooks_present

   logical FUNCTION lifecycle_gas_hooks_complete (hooks)
      type(material_lifecycle_hooks_type), intent(in) :: hooks
      lifecycle_gas_hooks_complete = associated(hooks%land_init) .and. &
         associated(hooks%land_read_restart) .and. associated(hooks%land_write_restart) .and. &
         associated(hooks%land_history) .and. associated(hooks%land_final)
   END FUNCTION lifecycle_gas_hooks_complete

   logical FUNCTION lifecycle_particle_hooks_complete (hooks)
      type(material_lifecycle_hooks_type), intent(in) :: hooks
      lifecycle_particle_hooks_complete = associated(hooks%route_init) .and. &
         associated(hooks%route_read_restart) .and. associated(hooks%route_write_restart) .and. &
         associated(hooks%route_calc) .and. associated(hooks%route_history) .and. &
         associated(hooks%route_flush_history) .and. associated(hooks%route_final)
   END FUNCTION lifecycle_particle_hooks_complete

   logical FUNCTION lifecycle_reaction_hooks_present (hooks)
      type(material_lifecycle_hooks_type), intent(in) :: hooks
      lifecycle_reaction_hooks_present = associated(hooks%lake_step) .or. &
         associated(hooks%wetland_decomp) .or. associated(hooks%soil_step)
   END FUNCTION lifecycle_reaction_hooks_present

   character(len=16) FUNCTION family_label (family)
      integer, intent(in) :: family

      SELECT CASE (family)
      CASE (FAMILY_ISOTOPE)
         family_label = 'isotope'
      CASE (FAMILY_SOLUTE)
         family_label = 'solute'
      CASE (FAMILY_PARTICLE)
         family_label = 'particle'
      CASE (FAMILY_GAS)
         family_label = 'gas'
      CASE DEFAULT
         family_label = 'unresolved'
      END SELECT
   END FUNCTION family_label

   character(len=16) FUNCTION owner_label (owner)
      integer, intent(in) :: owner

      SELECT CASE (owner)
      CASE (STATE_OWNER_GENERIC_WATER)
         owner_label = 'generic_water'
      CASE (STATE_OWNER_PROVIDER)
         owner_label = 'provider'
      CASE DEFAULT
         owner_label = 'unknown'
      END SELECT
   END FUNCTION owner_label

   character(len=16) FUNCTION reaction_label (reaction)
      integer, intent(in) :: reaction

      SELECT CASE (reaction)
      CASE (REACTION_NONE)
         reaction_label = 'none'
      CASE (REACTION_FIRST_ORDER)
         reaction_label = 'first_order'
      CASE (REACTION_PROVIDER)
         reaction_label = 'provider'
      CASE DEFAULT
         reaction_label = 'unknown'
      END SELECT
   END FUNCTION reaction_label

   character(len=32) FUNCTION lifecycle_hook_groups (hooks)
      type(material_lifecycle_hooks_type), intent(in) :: hooks
      logical :: has_land, has_reaction, has_route

      has_land = associated(hooks%land_init) .or. associated(hooks%land_read_restart) .or. &
         associated(hooks%land_write_restart) .or. associated(hooks%report) .or. &
         associated(hooks%flush_acc_fluxes) .or. associated(hooks%accumulate_fluxes) .or. &
         associated(hooks%land_history) .or. associated(hooks%land_save_lulcc) .or. &
         associated(hooks%land_history_sidecar_write) .or. associated(hooks%land_history_sidecar_read) .or. &
         associated(hooks%land_remap_lulcc) .or. associated(hooks%land_reload_lulcc) .or. &
         associated(hooks%publish_levee_flood) .or. associated(hooks%publish_flood) .or. &
         associated(hooks%land_final)
      has_reaction = lifecycle_reaction_hooks_present(hooks)
      has_route = lifecycle_route_hooks_present(hooks)

      lifecycle_hook_groups = '-'
      IF (has_land) lifecycle_hook_groups = 'land'
      IF (has_reaction) THEN
         IF (trim(lifecycle_hook_groups) == '-') THEN
            lifecycle_hook_groups = 'reaction'
         ELSE
            lifecycle_hook_groups = trim(lifecycle_hook_groups)//'+reaction'
         ENDIF
      ENDIF
      IF (has_route) THEN
         IF (trim(lifecycle_hook_groups) == '-') THEN
            lifecycle_hook_groups = 'route'
         ELSE
            lifecycle_hook_groups = trim(lifecycle_hook_groups)//'+route'
         ENDIF
      ENDIF
   END FUNCTION lifecycle_hook_groups

   SUBROUTINE validate_provider_declaration (name, aliases, provider, family, owner, reaction)
      character(len=*), intent(in) :: name, aliases, provider
      integer, intent(in) :: family, owner, reaction
      integer :: i, j, nitems, start_pos, comma_pos, end_pos, list_len, nseen
      character(len=32), allocatable :: seen(:)
      character(len=32) :: item

      IF (len_trim(name) == 0) CALL lifecycle_error ('provider registration name is empty')
      IF (len_trim(provider) == 0) CALL lifecycle_error ('provider identity is empty')
      IF (family < FAMILY_ISOTOPE .or. family > FAMILY_GAS) &
         CALL lifecycle_error ('provider declared an unknown family')
      IF (owner /= STATE_OWNER_GENERIC_WATER .and. owner /= STATE_OWNER_PROVIDER) &
         CALL lifecycle_error ('provider declared an unknown state owner')
      IF (reaction /= REACTION_NONE .and. reaction /= REACTION_FIRST_ORDER .and. &
          reaction /= REACTION_PROVIDER) &
         CALL lifecycle_error ('provider declared an unknown reaction mode')
      IF ((family == FAMILY_GAS .or. family == FAMILY_PARTICLE) .and. &
          owner /= STATE_OWNER_PROVIDER) &
         CALL lifecycle_error ('gas and particle providers require provider-owned state')
      IF ((family == FAMILY_ISOTOPE .or. family == FAMILY_PARTICLE) .and. &
          reaction /= REACTION_NONE) &
         CALL lifecycle_error ('isotope and particle reaction shape is not applicable')
      IF (reaction == REACTION_FIRST_ORDER .and. &
          (family /= FAMILY_SOLUTE .or. owner /= STATE_OWNER_GENERIC_WATER)) &
         CALL lifecycle_error ('first-order reaction requires generic-water solute state')

      list_len = len_trim(aliases)
      nitems = 1
      DO i = 1, list_len
         IF (aliases(i:i) == ',') nitems = nitems + 1
      ENDDO
      allocate(seen(nitems))
      seen = ''
      nseen = 0
      start_pos = 1
      DO WHILE (start_pos <= list_len)
         comma_pos = index(aliases(start_pos:list_len), ',')
         IF (comma_pos == 0) THEN
            end_pos = list_len
         ELSE
            end_pos = start_pos + comma_pos - 2
         ENDIF
         IF (end_pos >= start_pos) THEN
            IF (len_trim(adjustl(aliases(start_pos:end_pos))) > len(item)) &
               CALL lifecycle_error ('provider alias exceeds material-name length')
            item = material_upper(adjustl(trim(aliases(start_pos:end_pos))))
            IF (len_trim(item) > 0) THEN
               IF (trim(item) == trim(material_upper(name))) &
                  CALL lifecycle_error ('duplicate provider name in alias list')
               DO j = 1, nseen
                  IF (trim(item) == trim(seen(j))) CALL lifecycle_error ('duplicate provider alias')
               ENDDO
               nseen = nseen + 1
               seen(nseen) = item
            ENDIF
         ENDIF
         IF (comma_pos == 0) EXIT
         start_pos = end_pos + 2
      ENDDO
      deallocate(seen)
   END SUBROUTINE validate_provider_declaration

   SUBROUTINE lifecycle_error (message)
      character(len=*), intent(in) :: message
      CALL CoLM_stop ('MOD_MT_Lifecycle: '//trim(message))
   END SUBROUTINE lifecycle_error

END MODULE MOD_MT_Lifecycle
#endif
