#include <define.h>

#if (defined MaterialTransport) && (defined CATCHMENT)
#error "MaterialTransport does not support the catchment (unit-catchment) grid (CATCHMENT); use GRID or UNSTRUCTURED (vector)."
#endif

#ifdef MaterialTransport
MODULE MOD_MT_Defs

   USE MOD_Precision
   USE MOD_Namelist, only: DEF_USE_VariablySaturatedFlow
   USE MOD_MT_Namelist, only: mt_read_namelist, &
      DEF_USE_ISOTOPE, DEF_ISOTOPE_SPECIES, DEF_ISOTOPE_PARAM_FILE, &
      DEF_SOLUTE_NAMES, DEF_SOLUTE_PARAM_FILES, DEF_SEDIMENT_PARAM_FILE, &
      DEF_USE_METHANE, DEF_METHANE_PARAM_FILE, &
      DEF_MT_AQUIFER_MIXING_WATER_MM, &
      DEF_ISO_USE_FRACTIONATION, DEF_ISO_REF_RATIO_H2_18O, DEF_ISO_INIT_DELTA_H2_18O, &
      DEF_ISO_REF_RATIO_HDO, DEF_ISO_INIT_DELTA_HDO
   USE, INTRINSIC :: IEEE_ARITHMETIC, only: ieee_is_finite

   IMPLICIT NONE
   SAVE

   ! Each material is classified by three independent values:
   !
   !   family_id     : what the material physically is
   !     FAMILY_ISOTOPE  water isotopologue (H2_18O, HDO); moves with the land
   !                     and river water, fractionates on phase change and
   !                     leaves with the evaporated water.
   !     FAMILY_SOLUTE   dissolved, nonvolatile species (e.g. Cl); moves with
   !                     the water, stays behind on evaporation, may have a
   !                     dissolved limit (max_dissolved_conc) with the excess
   !                     kept as solid, and may decay at first order.
   !     FAMILY_PARTICLE suspended / bed material (sediment); exists only in
   !                     the river-lake flow, with its own state and transport.
   !     FAMILY_GAS      gas with dissolved and gaseous phases in the soil
   !                     column (methane); its own state and processes.
   !
   !   state_owner   : which code holds and moves the state
   !     STATE_OWNER_GENERIC_WATER  the shared water-carried transport
   !                                (isotopes and solutes)
   !     STATE_OWNER_PROVIDER       the species module, through the hooks it
   !                                registers (particle and gas)
   !
   !   reaction_mode : how the amount changes besides transport
   !     REACTION_NONE         none (also "not applicable" for isotopes and
   !                           particles)
   !     REACTION_FIRST_ORDER  first-order decay; generic-water solutes only
   !     REACTION_PROVIDER     computed by the species module (methane)
   integer, parameter :: FAMILY_UNRESOLVED = 0
   integer, parameter :: FAMILY_ISOTOPE    = 1
   integer, parameter :: FAMILY_SOLUTE     = 2
   integer, parameter :: FAMILY_PARTICLE   = 3
   integer, parameter :: FAMILY_GAS        = 4

   integer, parameter :: STATE_OWNER_UNKNOWN       = 0
   integer, parameter :: STATE_OWNER_GENERIC_WATER = 1
   integer, parameter :: STATE_OWNER_PROVIDER      = 2

   integer, parameter :: REACTION_NONE        = 0
   integer, parameter :: REACTION_FIRST_ORDER = 1
   integer, parameter :: REACTION_PROVIDER    = 2

   ! Fixed-width, byte-exact descriptor serialization shared by every
   ! restart domain.  Keeping the formatter here prevents land and routing
   ! checkpoints from drifting into different notions of material identity.
   integer, parameter :: MATERIAL_DESCRIPTOR_IDENTITY_WIDTH = 384

   type :: material_info_type
      character(len=32) :: name
      character(len=16) :: category
      integer           :: family_id = FAMILY_UNRESOLVED
      integer           :: state_owner = STATE_OWNER_UNKNOWN
      integer           :: reaction_mode = REACTION_NONE
      ! unit_kind controls diagnostic labels only; it does not rescale forcing
      ! or state values. Numeric inputs must already use the declared ratio.
      character(len=32) :: unit_kind = 'material_per_water'
      real(r8)          :: mol_weight
      real(r8)          :: ref_ratio
      real(r8)          :: init_delta
      ! Explicit non-isotope concentrations.  init_delta remains the legacy
      ! fallback so existing namelists retain their former signatures.
      real(r8)          :: init_conc
      real(r8)          :: precip_default_conc
      real(r8)          :: vapor_default_conc
      ! Positive dissolved concentration ceiling. huge() disables saturation
      ! for backward-compatible species that do not declare a solubility.
      real(r8)          :: max_dissolved_conc = huge(1.0_r8)
      real(r8)          :: decay_rate
      ! Indirect N2O emission factor at the land-to-river entry
      ! [kg N2O-N per kg N]; 0 disables it.  Set per solute in its
      ! &nl_colm_solute_parameter file.
      real(r8)          :: n2o_ef5 = 0._r8
      ! Ionic charge is registry metadata; no electroneutrality solver exists.
      integer           :: charge = 0
   end type material_info_type

   integer :: nmts = 0
   type(material_info_type), allocatable :: materials(:)

      real(r8), parameter :: mt_tiny = 1.0e-30_r8
      ! Numerical zero and physical denominator floors are deliberately
      ! separate.  mt_tiny is only an arithmetic nonzero guard.  Transport
      ! ratios from finite water pools need a larger physically resolved water
      ! floor, otherwise near-dry pools can export arbitrarily large R =
      ! material/water before downstream sanity caps see the problem.
      ! Units follow land water storage conventions (kg/m2 ~= mm).
      real(r8), parameter :: mt_water_min_for_ratio = 1.0e-12_r8
      ! ------------------------------------------------------------------
      ! DELTA DIAGNOSTIC INVARIANT
      !
      ! Any official isotope delta diagnostic must aggregate material material
      ! and water first, then convert once:
      !
      !    delta = mass_to_delta(sum(material_mass), sum(water_mass), ref_ratio)
      !
      ! Never average per-step/per-cell delta values, and do not average R
      ! unless that R is explicitly weighted by its matching water amount.
      ! Arithmetic means of delta or unweighted R bias low-water states and can
      ! create apparent isotope drift even when material mass is conserved.
      ! ------------------------------------------------------------------
      ! Physical minimum water mass (kg/m² ≈ mm) for delta diagnostics.
   ! Below this the water pool is floating-point noise (e.g. wa ≈ 1e-5
   ! mm in a patch that just went near zero), and R_sample = mt/water
   ! amplifies round-off into ±10000‰ deltas that are not physical.
   ! 1e-3 kg/m² = 1 μm water depth; this keeps legitimate trace water
   ! (dew/frost puddles, residual soil moisture) in the diagnostic while
   ! filtering pure FP dust. Previously 1e-6 (1 nm) was too permissive —
   ! flux deltas (qinfl/qcharge) with water_mass = rate*deltim near the
   ! FP floor produced arbitrarily large ratios that flowed into
   ! mt_delta_* history / check outputs.
   real(r8), parameter :: mt_water_min_for_delta = 1.0e-3_r8
   ! Flux isotope deltas need a stricter accumulated-water denominator than
   ! storage deltas. Below this history-window water loss, the flux signature is
   ! dominated by NSS/storage tendency or round-off rather than a resolved flux.
   real(r8), parameter :: mt_flux_water_min_for_delta = 1.0e-1_r8
   ! Hard sanity cap on |delta| (permil). Even extreme Rayleigh
   ! depletion in snowpack or evap residuals rarely exceeds a few
   ! hundred permil for meteoric isotope ratios; anything past this is
   ! a sign of material/water state divergence, not real physics. Used
   ! by delta diagnostics to clamp to spval so non-physical spikes do
   ! not poison min/max reports or downstream post-processing. 2000 ‰
   ! gives comfortable headroom for strongly fractionated pools.
   real(r8), parameter :: mt_delta_sanity_max = 2.0e3_r8

   PUBLIC :: material_defs_init, material_defs_validate, mt_defs_final
   PUBLIC :: mass_to_delta, delta_to_R, R_to_mass
   PUBLIC :: material_aquifer_isotope_state_valid
   PUBLIC :: material_aquifer_actual_water, material_aquifer_actual_mass, material_aquifer_isotope_ratio
   PUBLIC :: material_is_isotope, material_has_decay
   PUBLIC :: material_is_solute, material_is_particle, material_is_gas
   PUBLIC :: material_is_nonvolatile_solute, mt_uses_generic_water
   PUBLIC :: material_get_family, material_get_state_owner, material_get_reaction_mode
   PUBLIC :: material_uses_delta_diagnostics, material_can_use_fixed_signature
   PUBLIC :: material_concentration_units
   PUBLIC :: material_init_water_ratio, material_precip_default_ratio, material_vapor_default_ratio
   PUBLIC :: material_param_file_for_index, material_index_for_name
   PUBLIC :: material_build_descriptor_identity
   PUBLIC :: material_lower, material_upper
   PUBLIC :: nmts, materials
   PUBLIC :: material_info_type
   PUBLIC :: FAMILY_UNRESOLVED, FAMILY_ISOTOPE, FAMILY_SOLUTE, FAMILY_PARTICLE, FAMILY_GAS
   PUBLIC :: STATE_OWNER_UNKNOWN, STATE_OWNER_GENERIC_WATER, STATE_OWNER_PROVIDER
   PUBLIC :: REACTION_NONE, REACTION_FIRST_ORDER, REACTION_PROVIDER
   PUBLIC :: MATERIAL_DESCRIPTOR_IDENTITY_WIDTH
      PUBLIC :: mt_tiny, mt_water_min_for_ratio, &
                mt_water_min_for_delta, mt_flux_water_min_for_delta, &
                mt_delta_sanity_max

CONTAINS

   pure real(r8) FUNCTION material_aquifer_actual_water (wa, reference_water)
      real(r8), intent(in) :: wa, reference_water
      material_aquifer_actual_water = wa + reference_water
   END FUNCTION material_aquifer_actual_water

   pure real(r8) FUNCTION material_aquifer_actual_mass (relative_mass, reference_mass)
      real(r8), intent(in) :: relative_mass, reference_mass
      material_aquifer_actual_mass = relative_mass + reference_mass
   END FUNCTION material_aquifer_actual_mass

   pure real(r8) FUNCTION material_aquifer_isotope_ratio (wa, relative_mass, &
                                                        reference_water, reference_mass, fallback)
      real(r8), intent(in) :: wa, relative_mass, reference_water, reference_mass, fallback
      real(r8) :: actual_water

      actual_water = material_aquifer_actual_water(wa, reference_water)
      material_aquifer_isotope_ratio = fallback
      IF (abs(actual_water) > mt_water_min_for_ratio) &
         material_aquifer_isotope_ratio = &
            material_aquifer_actual_mass(relative_mass, reference_mass) / actual_water
   END FUNCTION material_aquifer_isotope_ratio

   pure logical FUNCTION material_aquifer_isotope_state_valid (water_mass, isotope_mass, ref_ratio, &
                                                             reference_water, reference_mass)
      real(r8), intent(in) :: water_mass, isotope_mass, ref_ratio
      real(r8), intent(in), optional :: reference_water, reference_mass
      real(r8) :: dust_mass, ratio, ref_water, ref_mass, actual_water, actual_mass

      material_aquifer_isotope_state_valid = .false.
      IF (.not. ieee_is_finite(water_mass) .or. .not. ieee_is_finite(isotope_mass)) RETURN
      IF (.not. ieee_is_finite(ref_ratio) .or. ref_ratio <= 0._r8) RETURN
      ref_water = 0._r8
      ref_mass = 0._r8
      IF (present(reference_water)) ref_water = reference_water
      IF (present(reference_mass)) ref_mass = reference_mass
      IF (.not. ieee_is_finite(ref_water) .or. .not. ieee_is_finite(ref_mass)) RETURN
      IF (ref_water < 0._r8 .or. ref_mass < 0._r8) RETURN
      IF (ref_water > huge(1._r8) - abs(water_mass)) RETURN
      IF (ref_mass > huge(1._r8) - abs(isotope_mass)) RETURN
      actual_water = material_aquifer_actual_water(water_mass, ref_water)
      actual_mass = material_aquifer_actual_mass(isotope_mass, ref_mass)
      ! Below the resolved water floor there is no aquifer carrier.  Allow
      ! only isotope mass at the corresponding round-off scale, not a
      ! material orphan that would otherwise be cleared as a residual.
      dust_mass = max(mt_water_min_for_ratio * ref_ratio, &
                      16._r8 * epsilon(1._r8) * max(abs(isotope_mass), ref_mass))
      IF (abs(actual_water) <= mt_water_min_for_ratio) THEN
         material_aquifer_isotope_state_valid = abs(actual_water) <= mt_water_min_for_ratio &
            .and. abs(actual_mass) <= dust_mass
      ELSE
         ! R=0 is the admissible isotope-free (delta=-1000 permil) limit.
         ! Check signs and the quotient range before division: debug builds
         ! trap floating-point overflow before ieee_is_finite can inspect it.
         IF (ref_water > 0._r8 .and. actual_water < 0._r8) RETURN
         IF (actual_water * sign(1._r8, actual_mass) < 0._r8 .and. actual_mass /= 0._r8) RETURN
         IF (abs(actual_water) < 1._r8) THEN
            IF (abs(actual_mass) > huge(1._r8) * abs(actual_water)) RETURN
         ENDIF
         ratio = actual_mass / actual_water
         material_aquifer_isotope_state_valid = ieee_is_finite(ratio) .and. ratio >= 0._r8
      ENDIF
   END FUNCTION material_aquifer_isotope_state_valid

   SUBROUTINE material_build_descriptor_identity (identity, transport_only)
      IMPLICIT NONE
      integer, allocatable, intent(out) :: identity(:,:)
      logical, intent(in), optional :: transport_only

      integer :: imt, k, nidentity, iidentity
      logical :: only_transport
      character(len=MATERIAL_DESCRIPTOR_IDENTITY_WIDTH) :: descriptor

      only_transport = .false.
      IF (present(transport_only)) only_transport = transport_only

      nidentity = 0
      DO imt = 1, nmts
         IF (only_transport .and. .not. mt_uses_generic_water(imt)) CYCLE
         nidentity = nidentity + 1
      ENDDO
      allocate(identity(MATERIAL_DESCRIPTOR_IDENTITY_WIDTH, nidentity))
      identity = 0

      iidentity = 0
      DO imt = 1, nmts
         IF (only_transport .and. .not. mt_uses_generic_water(imt)) CYCLE
         iidentity = iidentity + 1
         descriptor = ''
         WRITE(descriptor, '(A,"|",A,"|",A,"|",4(I0,"|"),9(ES24.16E3,"|"))') &
            trim(materials(imt)%name), trim(materials(imt)%category), &
            trim(materials(imt)%unit_kind), materials(imt)%family_id, &
            materials(imt)%state_owner, materials(imt)%reaction_mode, &
            materials(imt)%charge, materials(imt)%mol_weight, &
            materials(imt)%ref_ratio, materials(imt)%init_delta, &
            materials(imt)%init_conc, materials(imt)%precip_default_conc, &
            materials(imt)%vapor_default_conc, materials(imt)%max_dissolved_conc, &
            materials(imt)%decay_rate, materials(imt)%n2o_ef5
         DO k = 1, len(descriptor)
            identity(k, iidentity) = iachar(descriptor(k:k))
         ENDDO
      ENDDO
   END SUBROUTINE material_build_descriptor_identity

   SUBROUTINE material_defs_init ()
   !-----------------------------------------------------------------------
   !  Build the descriptor table from the switches in MOD_MT_Namelist:
   !  isotopes first (order of DEF_ISOTOPE_SPECIES), then solutes (order of
   !  DEF_SOLUTE_NAMES), then the river sediment (DEF_USE_SEDIMENT), then
   !  methane (DEF_USE_METHANE). The family of each entry follows from the
   !  switch that put it there.
   !-----------------------------------------------------------------------
      USE MOD_SPMD_Task, only: p_is_master, CoLM_stop
      USE MOD_Namelist, only: DEF_USE_SEDIMENT
      IMPLICIT NONE
      integer :: i, j, niso, nsol, nsed, ngas
      character(len=32), allocatable :: iso_tokens(:), sol_tokens(:)

      CALL mt_read_namelist ()

      IF (allocated(materials)) RETURN   ! already initialised (e.g. by earlier all-rank call)

      niso = 0
      nsol = 0
      IF (DEF_USE_ISOTOPE) CALL split_names (DEF_ISOTOPE_SPECIES, iso_tokens, niso)
      CALL split_names (DEF_SOLUTE_NAMES, sol_tokens, nsol)
      IF (DEF_USE_ISOTOPE .and. niso == 0) &
         CALL CoLM_stop ('DEF_USE_ISOTOPE is on but DEF_ISOTOPE_SPECIES is empty.')
      nsed = 0
      IF (DEF_USE_SEDIMENT) THEN
#ifdef GridRiverLakeFlow
         nsed = 1
#else
         CALL CoLM_stop ('DEF_USE_SEDIMENT needs the grid river-lake flow (GridRiverLakeFlow).')
#endif
      ENDIF

      ngas = 0
      IF (DEF_USE_METHANE) THEN
#ifdef BGC
         ngas = 1
#else
         CALL CoLM_stop ('DEF_USE_METHANE needs a BGC build.')
#endif
      ENDIF

      nmts = niso + nsol + nsed + ngas
      IF (nmts == 0) RETURN
      allocate(materials(nmts))

      DO i = 1, niso
         SELECT CASE (trim(material_upper(iso_tokens(i))))
         CASE ('H2_18O')
            materials(i)%name       = 'H2_18O'
            materials(i)%mol_weight = 20.0_r8
            materials(i)%ref_ratio  = DEF_ISO_REF_RATIO_H2_18O
            materials(i)%init_delta = DEF_ISO_INIT_DELTA_H2_18O
         CASE ('HDO')
            materials(i)%name       = 'HDO'
            materials(i)%mol_weight = 19.0_r8
            materials(i)%ref_ratio  = DEF_ISO_REF_RATIO_HDO
            materials(i)%init_delta = DEF_ISO_INIT_DELTA_HDO
         CASE DEFAULT
            CALL CoLM_stop ('Unknown isotope species "'//trim(iso_tokens(i))// &
               '" in DEF_ISOTOPE_SPECIES; available: H2_18O, HDO.')
         END SELECT
         materials(i)%category = 'isotope'
         CALL set_mt_category_defaults (i)
         ! Isotopes start from and fall back to init_delta; the
         ! concentration fields only enter the descriptor identity.
         materials(i)%init_conc           = materials(i)%init_delta
         materials(i)%precip_default_conc = materials(i)%init_delta
         materials(i)%vapor_default_conc  = materials(i)%init_delta
         materials(i)%decay_rate = 0._r8
      ENDDO

      DO j = 1, nsol
         i = niso + j
         materials(i)%name = sanitize_ncname (sol_tokens(j))
         materials(i)%category = 'solute'
         CALL set_mt_category_defaults (i)
         materials(i)%mol_weight          = 0._r8
         materials(i)%ref_ratio           = 1._r8
         materials(i)%init_delta          = 0._r8
         materials(i)%init_conc           = 0._r8
         materials(i)%precip_default_conc = 0._r8
         materials(i)%vapor_default_conc  = 0._r8
         materials(i)%decay_rate = 0._r8
      ENDDO

      IF (nsed > 0) THEN
         ! State and parameters are owned by MOD_MT_Sediment, which attaches
         ! itself to this entry (mt_register_providers).
         i = niso + nsol + 1
         materials(i)%name = 'SEDIMENT'
         materials(i)%category = 'particle'
         CALL set_mt_category_defaults (i)
         materials(i)%mol_weight          = 1._r8
         materials(i)%ref_ratio           = 1._r8
         materials(i)%init_delta          = 0._r8
         materials(i)%init_conc           = 0._r8
         materials(i)%precip_default_conc = 0._r8
         materials(i)%vapor_default_conc  = 0._r8
         materials(i)%decay_rate = 0._r8
      ENDIF

      IF (ngas > 0) THEN
         ! State, parameters and reactions are owned by MOD_MT_Methane, which
         ! attaches itself to this entry (mt_register_providers).
         i = niso + nsol + nsed + 1
         materials(i)%name = 'CH4'
         materials(i)%category = 'gas'
         CALL set_mt_category_defaults (i)
         materials(i)%mol_weight          = 16.04_r8
         materials(i)%ref_ratio           = 1._r8
         materials(i)%init_delta          = 0._r8
         materials(i)%init_conc           = 0._r8
         materials(i)%precip_default_conc = 0._r8
         materials(i)%vapor_default_conc  = 0._r8
         materials(i)%decay_rate = 0._r8
      ENDIF

      ! Names become NetCDF variable names in land and river files.
      DO i = 2, nmts
         DO j = 1, i - 1
            IF (material_param_equal(materials(i)%name, materials(j)%name)) &
               CALL CoLM_stop ('Material name "'//trim(materials(i)%name)//'" is used twice.')
         ENDDO
      ENDDO

      IF (allocated(iso_tokens)) deallocate(iso_tokens)
      IF (allocated(sol_tokens)) deallocate(sol_tokens)

   END SUBROUTINE material_defs_init

   SUBROUTINE material_defs_validate ()
   ! Checks of the descriptors, after the species parameter files
   ! (MOD_MT_Solute: mt_solute_read_parameters) have been applied.
      USE MOD_SPMD_Task, only: p_is_master, CoLM_stop
      IMPLICIT NONE
      integer :: i

      IF (.not. allocated(materials)) RETURN

      IF (DEF_ISO_USE_FRACTIONATION .and. DEF_USE_ISOTOPE .and. p_is_master) THEN
         WRITE(*,'(A)') 'WARNING material_defs_validate: isotope fractionation is experimental; enabled paths include ' // &
            'precipitation/evaporation, phase change, transpiration, wetland, glacier, and waterbody approximations.'
      ENDIF

      DO i = 1, nmts
         IF (.not. ieee_is_finite(materials(i)%decay_rate)) THEN
            CALL CoLM_stop ('Non-finite decay_rate for '//trim(materials(i)%name)//'.')
         ELSEIF (materials(i)%decay_rate < 0._r8) THEN
            CALL CoLM_stop ('Negative decay_rate for '//trim(materials(i)%name)//'.')
         ENDIF
         CALL derive_mt_taxonomy (i)
         CALL validate_mt_descriptor (i)
         IF (DEF_USE_VariablySaturatedFlow .and. material_is_isotope(i)) THEN
            IF (.not. ieee_is_finite(DEF_MT_AQUIFER_MIXING_WATER_MM) .or. &
                DEF_MT_AQUIFER_MIXING_WATER_MM <= 0._r8) THEN
               CALL CoLM_stop ('Isotopes with variably saturated flow need a positive '// &
                  'DEF_MT_AQUIFER_MIXING_WATER_MM.')
            ENDIF
         ENDIF
      ENDDO

   END SUBROUTINE material_defs_validate


   SUBROUTINE split_names (csvstr, tokens, ntokens)
   ! Comma-separated names; blank entries are skipped.
      USE MOD_SPMD_Task, only: CoLM_stop
      IMPLICIT NONE
      character(len=*), intent(in) :: csvstr
      character(len=32), allocatable, intent(out) :: tokens(:)
      integer, intent(out) :: ntokens
      integer :: i, j, slen

      slen = len_trim(csvstr)
      allocate(tokens(max(1, count([(csvstr(i:i) == ',', i = 1, slen)]) + 1)))
      ntokens = 0
      j = 1
      DO i = 1, slen + 1
         IF (i > slen) THEN
            CALL add_token (csvstr(j:slen))
         ELSEIF (csvstr(i:i) == ',') THEN
            CALL add_token (csvstr(j:i-1))
            j = i + 1
         ENDIF
      ENDDO

   CONTAINS

      SUBROUTINE add_token (raw)
         character(len=*), intent(in) :: raw
         IF (len_trim(adjustl(raw)) == 0) RETURN
         IF (len_trim(adjustl(raw)) > 32) &
            CALL CoLM_stop ('Material name longer than 32 characters: '//trim(adjustl(raw)))
         ntokens = ntokens + 1
         tokens(ntokens) = adjustl(raw)
      END SUBROUTINE add_token

   END SUBROUTINE split_names


   SUBROUTINE material_param_file_for_index (imt, aliases, file_param, found)
   !-----------------------------------------------------------------------
   !  Parameter file of one material: DEF_ISOTOPE_PARAM_FILE for every
   !  isotope, the DEF_SOLUTE_PARAM_FILES entry for a solute ('NAME:path',
   !  or a plain path at the solute's position), DEF_SEDIMENT_PARAM_FILE for
   !  the sediment, DEF_METHANE_PARAM_FILE for methane. 'null' means none.
   !-----------------------------------------------------------------------
      USE MOD_SPMD_Task, only: CoLM_stop
      IMPLICIT NONE
      integer, intent(in) :: imt
      character(len=*), intent(in) :: aliases
      character(len=*), intent(out) :: file_param
      logical, intent(out) :: found

      integer :: start_pos, end_pos, list_len, colon_pos, positional_index, isol, i
      character(len=len(DEF_SOLUTE_PARAM_FILES)) :: entry, key, value
      logical :: matched

      file_param = ''
      found = .false.
      matched = .false.
      IF (imt < 1 .or. .not. allocated(materials) .or. imt > nmts) RETURN

      IF (material_is_isotope(imt)) THEN
         file_param = DEF_ISOTOPE_PARAM_FILE
         found = len_trim(file_param) > 0 .and. trim(material_lower(file_param)) /= 'null'
         RETURN
      ENDIF
      IF (material_is_particle(imt)) THEN
         file_param = DEF_SEDIMENT_PARAM_FILE
         found = len_trim(file_param) > 0 .and. trim(material_lower(file_param)) /= 'null'
         RETURN
      ENDIF
      IF (material_is_gas(imt)) THEN
         file_param = DEF_METHANE_PARAM_FILE
         found = len_trim(file_param) > 0 .and. trim(material_lower(file_param)) /= 'null'
         RETURN
      ENDIF
      IF (.not. material_is_solute(imt)) RETURN

      ! Position of this solute among the solutes.
      isol = 0
      DO i = 1, imt
         IF (material_is_solute(i)) isol = isol + 1
      ENDDO

      list_len = len_trim(DEF_SOLUTE_PARAM_FILES)
      IF (list_len <= 0) RETURN
      IF (trim(material_lower(DEF_SOLUTE_PARAM_FILES)) == 'null') RETURN

      positional_index = 0
      start_pos = 1
      DO WHILE (start_pos <= list_len)
         end_pos = material_param_next_entry_end(DEF_SOLUTE_PARAM_FILES, start_pos, list_len)
         entry = adjustl(trim(DEF_SOLUTE_PARAM_FILES(start_pos:end_pos)))
         IF (len_trim(entry) > 0) THEN
            colon_pos = index(entry, ':')
            IF (colon_pos > 0) THEN
               key = adjustl(trim(entry(:colon_pos-1)))
               value = adjustl(trim(entry(colon_pos+1:)))
               IF (len_trim(key) <= 0 .or. len_trim(value) <= 0) THEN
                  CALL CoLM_stop('MOD_MT_Defs: empty solute parameter file mapping entry: '//trim(entry))
               ENDIF
               IF (len_trim(value) > len(file_param)) THEN
                  CALL CoLM_stop('MOD_MT_Defs: solute parameter file path exceeds supported length')
               ENDIF
               IF (.not. matched .and. material_param_key_matches(imt, key, aliases)) THEN
                  file_param = trim(value)
                  found = trim(material_lower(file_param)) /= 'null'
                  matched = .true.
               ENDIF
            ELSE
               IF (len_trim(entry) > len(file_param)) THEN
                  CALL CoLM_stop('MOD_MT_Defs: solute parameter file path exceeds supported length')
               ENDIF
               positional_index = positional_index + 1
               IF (.not. matched .and. positional_index == isol) THEN
                  file_param = trim(entry)
                  found = len_trim(file_param) > 0 .and. trim(material_lower(file_param)) /= 'null'
                  matched = .true.
               ENDIF
            ENDIF
         ENDIF
         start_pos = end_pos + 2
      ENDDO
   END SUBROUTINE material_param_file_for_index

   integer FUNCTION material_index_for_name (name, aliases)
      USE MOD_SPMD_Task, only: CoLM_stop
      IMPLICIT NONE
      character(len=*), intent(in) :: name, aliases
      integer :: imt, matches

      material_index_for_name = 0
      matches = 0
      IF (.not. allocated(materials)) RETURN

      DO imt = 1, nmts
         IF (material_name_matches(materials(imt)%name, name, aliases)) THEN
            matches = matches + 1
            material_index_for_name = imt
         ENDIF
      ENDDO
      IF (matches > 1) &
         CALL CoLM_stop('MOD_MT_Defs: material name/aliases match multiple configured materials')
   END FUNCTION material_index_for_name

   logical FUNCTION material_name_matches (candidate, name, aliases)
      IMPLICIT NONE
      character(len=*), intent(in) :: candidate, name, aliases
      integer :: start_pos, end_pos, list_len

      material_name_matches = material_param_equal(candidate, name)
      IF (material_name_matches) RETURN

      list_len = len_trim(aliases)
      start_pos = 1
      DO WHILE (start_pos <= list_len)
         end_pos = material_param_next_entry_end(aliases, start_pos, list_len)
         IF (end_pos >= start_pos) THEN
            IF (material_param_equal(candidate, aliases(start_pos:end_pos))) THEN
               material_name_matches = .true.
               RETURN
            ENDIF
         ENDIF
         start_pos = end_pos + 2
      ENDDO
   END FUNCTION material_name_matches

   logical FUNCTION material_param_key_matches (imt, raw_key, aliases)
      IMPLICIT NONE
      integer, intent(in) :: imt
      character(len=*), intent(in) :: raw_key, aliases

      integer :: start_pos, end_pos, list_len
      character(len=128) :: one_alias

      material_param_key_matches = .false.
      IF (imt < 1 .or. .not. allocated(materials) .or. imt > nmts) RETURN
      IF (material_param_equal(raw_key, materials(imt)%name)) THEN
         material_param_key_matches = .true.
         RETURN
      ENDIF

      list_len = len_trim(aliases)
      start_pos = 1
      DO WHILE (start_pos <= list_len)
         end_pos = material_param_next_entry_end(aliases, start_pos, list_len)
         one_alias = adjustl(trim(aliases(start_pos:end_pos)))
         IF (material_param_equal(raw_key, one_alias)) THEN
            material_param_key_matches = .true.
            RETURN
         ENDIF
         start_pos = end_pos + 2
      ENDDO
   END FUNCTION material_param_key_matches

   logical FUNCTION material_param_equal (a, b)
      IMPLICIT NONE
      character(len=*), intent(in) :: a, b
      material_param_equal = trim(material_lower(a)) == trim(material_lower(b))
   END FUNCTION material_param_equal

   integer FUNCTION material_param_next_entry_end (raw_list, start_pos, list_len)
      IMPLICIT NONE
      character(len=*), intent(in) :: raw_list
      integer, intent(in) :: start_pos, list_len
      integer :: i

      material_param_next_entry_end = list_len
      DO i = start_pos, list_len
         IF (raw_list(i:i) == ',' .or. raw_list(i:i) == ';') THEN
            material_param_next_entry_end = i - 1
            RETURN
         ENDIF
      ENDDO
   END FUNCTION material_param_next_entry_end

   FUNCTION material_lower (raw) RESULT(out)
      IMPLICIT NONE
      character(len=*), intent(in) :: raw
      character(len=max(1,len_trim(raw))) :: out
      integer :: i, ia

      out = adjustl(trim(raw))
      DO i = 1, len_trim(out)
         ia = iachar(out(i:i))
         IF (ia >= iachar('A') .and. ia <= iachar('Z')) THEN
            out(i:i) = achar(ia + iachar('a') - iachar('A'))
         ENDIF
      ENDDO
   END FUNCTION material_lower

   FUNCTION material_upper (raw) RESULT(out)
      IMPLICIT NONE
      character(len=*), intent(in) :: raw
      character(len=max(1,len_trim(raw))) :: out
      integer :: i, ia

      out = adjustl(trim(raw))
      DO i = 1, len_trim(out)
         ia = iachar(out(i:i))
         IF (ia >= iachar('a') .and. ia <= iachar('z')) THEN
            out(i:i) = achar(ia - iachar('a') + iachar('A'))
         ENDIF
      ENDDO
   END FUNCTION material_upper

   SUBROUTINE mt_defs_final ()
      IF (allocated(materials)) deallocate(materials)
      nmts = 0
   END SUBROUTINE mt_defs_final

   !-------------------------------------------------------------------
   ! Material taxonomy.
   !
   ! The transport math is intentionally expressed as "water flux times
   ! material per water" for all categories. For isotopes this is R_sample;
   ! for solutes it is concentration per unit water.
   !-------------------------------------------------------------------

   SUBROUTINE set_mt_category_defaults (imt)
      integer, intent(in) :: imt

      IF (.not. allocated(materials)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN

      materials(imt)%unit_kind = 'material_per_water'
      materials(imt)%charge = 0
      materials(imt)%max_dissolved_conc = huge(1.0_r8)
      materials(imt)%family_id = FAMILY_UNRESOLVED
      materials(imt)%state_owner = STATE_OWNER_UNKNOWN
      materials(imt)%reaction_mode = REACTION_NONE

      SELECT CASE (trim(materials(imt)%category))
      CASE ('isotope')
         materials(imt)%unit_kind = 'ratio'
         materials(imt)%family_id = FAMILY_ISOTOPE
         materials(imt)%state_owner = STATE_OWNER_GENERIC_WATER
      CASE ('solute')
         materials(imt)%family_id = FAMILY_SOLUTE
         materials(imt)%state_owner = STATE_OWNER_GENERIC_WATER
      CASE ('particle')
         materials(imt)%unit_kind = 'volume_fraction'
         materials(imt)%family_id = FAMILY_PARTICLE
         materials(imt)%state_owner = STATE_OWNER_PROVIDER
      CASE ('gas')
         materials(imt)%unit_kind = 'species_owned'
         materials(imt)%family_id = FAMILY_GAS
         materials(imt)%state_owner = STATE_OWNER_PROVIDER
      END SELECT
   END SUBROUTINE set_mt_category_defaults

   SUBROUTINE derive_mt_taxonomy (imt)
      integer, intent(in) :: imt

      IF (.not. allocated(materials)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN

      SELECT CASE (trim(materials(imt)%category))
      CASE ('isotope')
         materials(imt)%family_id = FAMILY_ISOTOPE
      CASE ('solute')
         materials(imt)%family_id = FAMILY_SOLUTE
      CASE ('particle')
         materials(imt)%family_id = FAMILY_PARTICLE
      CASE ('gas')
         materials(imt)%family_id = FAMILY_GAS
      CASE DEFAULT
         materials(imt)%family_id = FAMILY_UNRESOLVED
      END SELECT

      SELECT CASE (materials(imt)%family_id)
      CASE (FAMILY_ISOTOPE, FAMILY_SOLUTE)
         materials(imt)%state_owner = STATE_OWNER_GENERIC_WATER
      CASE (FAMILY_PARTICLE, FAMILY_GAS)
         materials(imt)%state_owner = STATE_OWNER_PROVIDER
      CASE DEFAULT
         materials(imt)%state_owner = STATE_OWNER_UNKNOWN
      END SELECT

      IF (materials(imt)%family_id == FAMILY_SOLUTE .and. &
          materials(imt)%decay_rate > 0._r8) THEN
         materials(imt)%reaction_mode = REACTION_FIRST_ORDER
      ELSE
         materials(imt)%reaction_mode = REACTION_NONE
      ENDIF
   END SUBROUTINE derive_mt_taxonomy

   SUBROUTINE validate_mt_descriptor (imt)
      USE MOD_SPMD_Task, only: p_is_master, CoLM_stop
      integer, intent(in) :: imt
      logical :: supported_unit

      IF (.not. allocated(materials)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN

      materials(imt)%unit_kind = material_lower(materials(imt)%unit_kind)
      SELECT CASE (materials(imt)%family_id)
      CASE (FAMILY_ISOTOPE, FAMILY_SOLUTE, FAMILY_PARTICLE, FAMILY_GAS)
         CONTINUE
      CASE DEFAULT
         CALL material_descriptor_error(imt, 'family_id', 'unsupported material family')
      END SELECT
      IF (materials(imt)%state_owner /= STATE_OWNER_GENERIC_WATER .and. &
          materials(imt)%state_owner /= STATE_OWNER_PROVIDER) THEN
         CALL material_descriptor_error(imt, 'state_owner', 'must resolve before transport initialization')
      ENDIF
      IF (materials(imt)%reaction_mode < REACTION_NONE .or. &
          materials(imt)%reaction_mode > REACTION_PROVIDER) THEN
         CALL material_descriptor_error(imt, 'reaction_mode', 'unsupported reaction capability')
      ENDIF
      IF ((material_is_isotope(imt) .or. material_is_particle(imt)) .and. &
          materials(imt)%reaction_mode /= REACTION_NONE) THEN
         CALL material_descriptor_error(imt, 'reaction_mode', &
            'isotope and particle families do not use the generic reaction capability')
      ENDIF
      supported_unit = trim(materials(imt)%unit_kind) == 'ratio' .or. &
         trim(materials(imt)%unit_kind) == 'material_per_water' .or. &
         trim(materials(imt)%unit_kind) == 'mass_fraction' .or. &
         trim(materials(imt)%unit_kind) == 'volume_fraction' .or. &
         trim(materials(imt)%unit_kind) == 'species_owned'
      IF (.not. supported_unit) THEN
         CALL material_descriptor_error(imt, 'unit_kind', 'unsupported unit convention')
      ENDIF
      IF (.not. ieee_is_finite(materials(imt)%mol_weight)) THEN
         CALL material_descriptor_error(imt, 'mol_weight', 'must be finite')
      ENDIF
      IF (.not. ieee_is_finite(materials(imt)%ref_ratio)) THEN
         CALL material_descriptor_error(imt, 'ref_ratio', 'must be finite')
      ENDIF
      IF (.not. ieee_is_finite(materials(imt)%init_delta)) THEN
         CALL material_descriptor_error(imt, 'init_delta', 'must be finite')
      ENDIF
      IF (.not. ieee_is_finite(materials(imt)%init_conc)) THEN
         CALL material_descriptor_error(imt, 'init_conc', 'must be finite')
      ENDIF
      IF (.not. ieee_is_finite(materials(imt)%precip_default_conc)) THEN
         CALL material_descriptor_error(imt, 'precip_default_conc', 'must be finite')
      ENDIF
      IF (.not. ieee_is_finite(materials(imt)%vapor_default_conc)) THEN
         CALL material_descriptor_error(imt, 'vapor_default_conc', 'must be finite')
      ENDIF
      IF (.not. ieee_is_finite(materials(imt)%max_dissolved_conc)) THEN
         CALL material_descriptor_error(imt, 'max_dissolved_conc', 'must be finite')
      ENDIF
      IF (.not. ieee_is_finite(materials(imt)%decay_rate)) THEN
         CALL material_descriptor_error(imt, 'decay_rate', 'must be finite')
      ENDIF
      IF (materials(imt)%mol_weight < 0._r8) THEN
         CALL material_descriptor_error(imt, 'mol_weight', 'must be non-negative')
      ENDIF
      IF (materials(imt)%ref_ratio <= 0._r8) THEN
         CALL material_descriptor_error(imt, 'ref_ratio', 'must be positive')
      ENDIF
      IF (material_is_isotope(imt)) THEN
         IF (materials(imt)%init_delta < -1000._r8) THEN
            CALL material_descriptor_error(imt, 'init_delta', 'must be at least -1000 permil')
         ENDIF
      ELSE
         IF (materials(imt)%init_conc < 0._r8) THEN
            CALL material_descriptor_error(imt, 'init_conc', 'must be non-negative')
         ENDIF
         IF (materials(imt)%precip_default_conc < 0._r8) THEN
            CALL material_descriptor_error(imt, 'precip_default_conc', 'must be non-negative')
         ENDIF
         IF (materials(imt)%vapor_default_conc < 0._r8) THEN
            CALL material_descriptor_error(imt, 'vapor_default_conc', 'must be non-negative')
         ENDIF
      ENDIF
      IF (materials(imt)%max_dissolved_conc <= 0._r8) THEN
         CALL material_descriptor_error(imt, 'max_dissolved_conc', 'must be positive')
      ENDIF
      IF (material_is_isotope(imt) .and. trim(materials(imt)%unit_kind) /= 'ratio') THEN
         CALL material_descriptor_error(imt, 'unit_kind', 'isotope materials require ratio')
      ENDIF
      IF (materials(imt)%charge /= 0 .and. .not. material_is_solute(imt)) THEN
         CALL material_descriptor_error(imt, 'charge', &
            'non-zero ionic charge is valid only for solute materials')
      ENDIF
      IF ((material_is_isotope(imt) .or. material_is_solute(imt)) .and. &
          materials(imt)%state_owner /= STATE_OWNER_GENERIC_WATER) THEN
         CALL material_descriptor_error(imt, 'state_owner', &
            'isotope and solute materials require generic land-water transport')
      ENDIF
      IF ((material_is_particle(imt) .or. material_is_gas(imt)) .and. &
          materials(imt)%state_owner /= STATE_OWNER_PROVIDER) THEN
         CALL material_descriptor_error(imt, 'state_owner', &
            'particle and gas materials require provider-owned state')
      ENDIF
      IF (materials(imt)%max_dissolved_conc < huge(1.0_r8) .and. &
          .not. material_is_solute(imt)) THEN
         CALL material_descriptor_error(imt, 'max_dissolved_conc', &
            'is valid only for solute materials')
      ENDIF
      IF (materials(imt)%decay_rate > 0._r8 .and. &
          .not. material_is_solute(imt)) THEN
         CALL material_descriptor_error(imt, 'decay_rate', &
            'first-order decay is supported only by solute materials')
      ENDIF

      IF (materials(imt)%state_owner == STATE_OWNER_PROVIDER .and. &
          .not. material_is_particle(imt) .and. &
          trim(materials(imt)%unit_kind) /= 'species_owned') THEN
         CALL material_descriptor_error(imt, 'unit_kind', &
            'species-owned state requires unit_kind=species_owned')
      ENDIF

#ifndef BGC
      IF ((trim(material_upper(materials(imt)%name)) == 'CH4' .or. &
           trim(material_upper(materials(imt)%name)) == 'METHANE') .and. &
          materials(imt)%state_owner == STATE_OWNER_PROVIDER) THEN
         CALL material_descriptor_error(imt, 'unit_kind', &
            'species-owned CH4 requires compiling with BGC; full soil/wetland/rice CH4 requires '// &
            'LULC_IGBP_PFT or LULC_IGBP_PC with real carbon, NPP, and root-respiration states')
      ENDIF
#endif

   CONTAINS

      SUBROUTINE material_descriptor_error (imt_local, field_name, reason)
         integer, intent(in) :: imt_local
         character(len=*), intent(in) :: field_name, reason

         IF (p_is_master) THEN
            write(*,'(5A)') 'ERROR material_defs_init: material "', trim(materials(imt_local)%name), &
               '" field ', trim(field_name), ': '//trim(reason)
         ENDIF
         CALL CoLM_stop()
      END SUBROUTINE material_descriptor_error

   END SUBROUTINE validate_mt_descriptor

   integer FUNCTION material_get_family (imt)
      integer, intent(in) :: imt
      material_get_family = FAMILY_UNRESOLVED
      IF (.not. allocated(materials)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN
      material_get_family = materials(imt)%family_id
   END FUNCTION material_get_family

   integer FUNCTION material_get_state_owner (imt)
      integer, intent(in) :: imt
      material_get_state_owner = STATE_OWNER_UNKNOWN
      IF (.not. allocated(materials)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN
      material_get_state_owner = materials(imt)%state_owner
   END FUNCTION material_get_state_owner

   integer FUNCTION material_get_reaction_mode (imt)
      integer, intent(in) :: imt
      material_get_reaction_mode = REACTION_NONE
      IF (.not. allocated(materials)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN
      material_get_reaction_mode = materials(imt)%reaction_mode
   END FUNCTION material_get_reaction_mode

   logical FUNCTION material_is_isotope (imt)
      integer, intent(in) :: imt
      material_is_isotope = .false.
      IF (.not. allocated(materials)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN
      material_is_isotope = materials(imt)%family_id == FAMILY_ISOTOPE
   END FUNCTION material_is_isotope

   logical FUNCTION material_is_solute (imt)
      integer, intent(in) :: imt
      material_is_solute = .false.
      IF (.not. allocated(materials)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN
      material_is_solute = materials(imt)%family_id == FAMILY_SOLUTE
   END FUNCTION material_is_solute

   logical FUNCTION material_has_decay (imt)
      integer, intent(in) :: imt
      material_has_decay = .false.
      IF (.not. allocated(materials)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN
      material_has_decay = materials(imt)%reaction_mode /= REACTION_NONE
   END FUNCTION material_has_decay

   logical FUNCTION material_is_nonvolatile_solute (imt)
      ! Nonvolatile solutes (decaying or not) are dissolved material,
      ! not water isotopologues. They move with liquid runoff/infiltration but
      ! atmospheric phase changes remove/add water only; evap/subl/transp
      ! leave solute behind and dew/frost do not import it unless a separate
      ! deposition source is provided.
      integer, intent(in) :: imt
      material_is_nonvolatile_solute = .false.
      IF (.not. allocated(materials)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN
      material_is_nonvolatile_solute = material_is_solute(imt) .and. &
         mt_uses_generic_water(imt)
   END FUNCTION material_is_nonvolatile_solute

   logical FUNCTION material_is_particle (imt)
      ! Particle materials (e.g. suspended sediment) carry concentration per
      ! unit water like solutes, but additionally settle and
      ! exchange mass with a non-water bed pool. Diagnostics are concentration
      ! / bed-change, never delta, so they share the non-isotope branches of
      ! material_uses_delta_diagnostics / material_can_use_fixed_signature.
      integer, intent(in) :: imt
      material_is_particle = .false.
      IF (.not. allocated(materials)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN
      material_is_particle = materials(imt)%family_id == FAMILY_PARTICLE
   END FUNCTION material_is_particle

   logical FUNCTION material_is_gas (imt)
      integer, intent(in) :: imt
      material_is_gas = .false.
      IF (.not. allocated(materials)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN
      material_is_gas = materials(imt)%family_id == FAMILY_GAS
   END FUNCTION material_is_gas

   logical FUNCTION mt_uses_generic_water (imt)
      ! Generic land-water material modules transport scalar signatures tied to
      ! canopy/soil/snow/aquifer water pools. Provider-owned particle and gas
      ! materials keep species-specific state and must not be initialized,
      ! transported, diagnosed, or restarted through those generic pools.
      integer, intent(in) :: imt
      mt_uses_generic_water = .false.
      IF (.not. allocated(materials)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN
      mt_uses_generic_water = materials(imt)%state_owner == STATE_OWNER_GENERIC_WATER
   END FUNCTION mt_uses_generic_water

   logical FUNCTION material_uses_delta_diagnostics (imt)
      integer, intent(in) :: imt
      ! Delta/NSS diagnostics are isotope-specific. Solutes
      ! store concentration per unit water and must not be converted
      ! with ref_ratio.
      material_uses_delta_diagnostics = material_is_isotope(imt)
   END FUNCTION material_uses_delta_diagnostics

   logical FUNCTION material_can_use_fixed_signature (imt)
      integer, intent(in) :: imt
      ! Fixed-signature resync/snap is a Phase-1 shortcut for isotope
      ! materials with inactive fractionation and no runtime forcing. Solutes
      ! may carry spatial gradients, restart
      ! memory, or source/sink history, so they must preserve the current
      ! pool-mixed concentration instead of being snapped to init_delta.
      material_can_use_fixed_signature = material_is_isotope(imt)
   END FUNCTION material_can_use_fixed_signature

   FUNCTION material_concentration_units (imt) RESULT(units)
      integer, intent(in) :: imt
      character(len=32) :: units

      units = 'material/water'
      IF (.not. allocated(materials)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN

      SELECT CASE (trim(materials(imt)%unit_kind))
      CASE ('ratio')
         units = 'R'
      CASE ('mass_fraction')
         units = 'kg/kg water'
      CASE ('volume_fraction')
         units = 'm3/m3 water'
      CASE ('species_owned')
         units = 'species-owned'
      CASE DEFAULT
         units = 'material/water'
      END SELECT
   END FUNCTION material_concentration_units

   real(r8) FUNCTION material_init_water_ratio (imt)
      integer, intent(in) :: imt

      material_init_water_ratio = 0._r8
      IF (.not. allocated(materials)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN

      IF (material_is_isotope(imt)) THEN
         material_init_water_ratio = delta_to_R(materials(imt)%init_delta, materials(imt)%ref_ratio)
      ELSE
         material_init_water_ratio = materials(imt)%init_conc
      ENDIF
   END FUNCTION material_init_water_ratio

   real(r8) FUNCTION material_precip_default_ratio (imt)
      integer, intent(in) :: imt

      material_precip_default_ratio = 0._r8
      IF (.not. allocated(materials)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN

      IF (material_is_isotope(imt)) THEN
         material_precip_default_ratio = delta_to_R(materials(imt)%init_delta, materials(imt)%ref_ratio)
      ELSE
         material_precip_default_ratio = materials(imt)%precip_default_conc
      ENDIF
   END FUNCTION material_precip_default_ratio

   real(r8) FUNCTION material_vapor_default_ratio (imt)
      integer, intent(in) :: imt

      material_vapor_default_ratio = 0._r8
      IF (.not. allocated(materials)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN

      IF (material_is_isotope(imt)) THEN
         material_vapor_default_ratio = delta_to_R(materials(imt)%init_delta, materials(imt)%ref_ratio)
      ELSE
         material_vapor_default_ratio = materials(imt)%vapor_default_conc
      ENDIF
   END FUNCTION material_vapor_default_ratio

      real(r8) FUNCTION mass_to_delta (mt_mass, water_mass, ref_ratio)
         USE MOD_Vars_Global, only: spval
         real(r8), intent(in) :: mt_mass, water_mass, ref_ratio
         real(r8) :: R_sample
         ! Public delta conversion point. Callers should pass already-aggregated
         ! material and water amounts for the requested history/check window.
         ! Do not call this on per-step values and then average the returned
         ! delta; aggregate mass and water first, convert last.
         ! abs(water_mass) guard: signed wa (aquifer-debt from wetland deficit)
      ! and its matching signed mt_wa still produce a well-defined ratio
      ! (both numerator and denominator share sign → R_sample positive).
      ! Physically non-negative pools (wliq/wice/wdsrf/wetwat/scv/ldew/flux)
      ! are unaffected since abs() equals value there.
      IF (abs(water_mass) > mt_tiny .and. ref_ratio > mt_tiny) THEN
         R_sample = mt_mass / water_mass
         IF (R_sample < -mt_tiny) THEN
            ! Sign mismatch between mt_mass and water_mass: upstream state
            ! machine is broken (e.g. mt_wa and wa drifted to opposite
            ! signs). Returning spval rather than silently clamping so
            ! history readers and assert tools can flag the cell instead
            ! of treating a bogus delta as physical.
#if (defined CoLMDEBUG)
            write(*,'(A,E12.5,A,E12.5)') &
               ' WARNING mass_to_delta: R_sample<0 mt=', mt_mass, ' water=', water_mass
#endif
            mass_to_delta = spval
         ELSE
            mass_to_delta = (R_sample / ref_ratio - 1.0_r8) * 1000.0_r8
         ENDIF
      ELSE
         mass_to_delta = 0.0_r8
      ENDIF
   END FUNCTION mass_to_delta

   pure real(r8) FUNCTION delta_to_R (delta, ref_ratio)
      real(r8), intent(in) :: delta, ref_ratio
      delta_to_R = ref_ratio * (1.0_r8 + delta / 1000.0_r8)
   END FUNCTION delta_to_R

   pure real(r8) FUNCTION R_to_mass (delta, water_mass, ref_ratio)
      real(r8), intent(in) :: delta, water_mass, ref_ratio
      R_to_mass = water_mass * delta_to_R(delta, ref_ratio)
   END FUNCTION R_to_mass

   !-------------------------------------------------------------------
   ! Sanitize a string for use as a NetCDF variable name component.
   ! Keeps alphanumerics, underscore and period; replaces anything else
   ! by dropping it. Land and river/lake history both splice materials(i)%name.
   !-------------------------------------------------------------------
   FUNCTION sanitize_ncname (raw) RESULT(clean)
      character(len=*), intent(in) :: raw
      character(len=32) :: clean
      integer :: i, j
      character(len=1) :: c
      clean = ''
      j = 0
      DO i = 1, len_trim(raw)
         c = raw(i:i)
         IF ((c >= 'A' .and. c <= 'Z') .or. (c >= 'a' .and. c <= 'z') .or. &
             (c >= '0' .and. c <= '9') .or. c == '_' .or. c == '.') THEN
            j = j + 1
            IF (j <= 32) clean(j:j) = c
         ENDIF
      ENDDO
      IF (len_trim(clean) == 0) clean = 'unnamed'
   END FUNCTION sanitize_ncname

END MODULE MOD_MT_Defs
#endif
