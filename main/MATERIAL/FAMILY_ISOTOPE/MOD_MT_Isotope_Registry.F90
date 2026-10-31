#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_Isotope_Registry

   USE MOD_Precision
   USE MOD_MT_Defs, only: materials, nmts, material_is_isotope, mt_tiny, material_lower
   USE MOD_SPMD_Task, only: CoLM_stop

   IMPLICIT NONE
   SAVE
   PRIVATE

   integer, parameter :: isotope_name_len = 32
   integer, parameter :: isotope_patterns_len = 160
   integer, parameter :: isotope_varname_len = 64

   ABSTRACT INTERFACE
      real(r8) FUNCTION isotope_alpha_temp_if (temp_k)
         IMPORT :: r8
         real(r8), intent(in) :: temp_k
      END FUNCTION isotope_alpha_temp_if

      real(r8) FUNCTION isotope_scalar_if ()
         IMPORT :: r8
      END FUNCTION isotope_scalar_if

   END INTERFACE

   TYPE :: isotope_physics_type
      character(len=isotope_name_len) :: name = ''
      character(len=isotope_patterns_len) :: name_patterns = ''
      real(r8) :: ref_ratio_hint = 0._r8
      real(r8) :: ref_ratio_tolerance = 0.1_r8
      character(len=isotope_varname_len) :: default_soil_init_varname = 'null'
      ! Species scaling of the Merlivat & Jouzel (1979) open-water kinetic
      ! factor, relative to 18O: 1.0 for 18O, 0.88 for HDO.  Published
      ! species property, so it lives beside the diffusivity ratio rather
      ! than as a name test inside the shared closure.
      real(r8) :: mj79_relative_factor = 1._r8
      procedure(isotope_alpha_temp_if), pointer, nopass :: alpha_liq_vap => null()
      procedure(isotope_alpha_temp_if), pointer, nopass :: alpha_ice_vap => null()
      procedure(isotope_scalar_if), pointer, nopass :: diffusivity_ratio_air => null()
      procedure(isotope_alpha_temp_if), pointer, nopass :: leaf_liquid_diffusivity => null()
   END TYPE isotope_physics_type

   TYPE(isotope_physics_type), allocatable :: isotope_physics(:)
   integer :: n_isotope_physics = 0

   PUBLIC :: isotope_physics_type
   PUBLIC :: register_isotope_physics
   PUBLIC :: find_isotope_physics
   PUBLIC :: isotope_fractionation_registered
   PUBLIC :: isotope_alpha_liq_vap
   PUBLIC :: isotope_alpha_ice_vap
   PUBLIC :: isotope_diffusivity_ratio_air
   PUBLIC :: isotope_mj79_relative_factor
   PUBLIC :: isotope_leaf_liquid_diffusivity
   PUBLIC :: isotope_default_soil_init_varname

CONTAINS

   SUBROUTINE register_isotope_physics (name, name_patterns, ref_ratio_hint, &
      default_soil_init_varname, ref_ratio_tolerance, mj79_relative_factor, &
      alpha_liq_vap_fn, alpha_ice_vap_fn, diffusivity_ratio_air_fn, &
      leaf_liquid_diffusivity_fn)
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: name_patterns
      real(r8), intent(in), optional :: ref_ratio_hint
      character(len=*), intent(in), optional :: default_soil_init_varname
      real(r8), intent(in), optional :: ref_ratio_tolerance
      real(r8), intent(in), optional :: mj79_relative_factor
      procedure(isotope_alpha_temp_if), optional :: alpha_liq_vap_fn
      procedure(isotope_alpha_temp_if), optional :: alpha_ice_vap_fn
      procedure(isotope_scalar_if), optional :: diffusivity_ratio_air_fn
      procedure(isotope_alpha_temp_if), optional :: leaf_liquid_diffusivity_fn
      integer :: idx

      IF (len_trim(name) <= 0) THEN
         CALL CoLM_stop ('MOD_MT_Isotope_Registry: cannot register isotope physics with empty name')
      ENDIF
      IF (find_registered_isotope_by_name(name) > 0) THEN
         CALL CoLM_stop ('MOD_MT_Isotope_Registry: duplicate isotope physics registration')
      ENDIF

      CALL ensure_isotope_physics_capacity(n_isotope_physics + 1)
      n_isotope_physics = n_isotope_physics + 1
      idx = n_isotope_physics

      isotope_physics(idx)%name = trim(name)
      isotope_physics(idx)%name_patterns = trim(name_patterns)
      IF (present(ref_ratio_hint)) isotope_physics(idx)%ref_ratio_hint = ref_ratio_hint
      IF (present(ref_ratio_tolerance)) isotope_physics(idx)%ref_ratio_tolerance = ref_ratio_tolerance
      IF (present(default_soil_init_varname)) &
         isotope_physics(idx)%default_soil_init_varname = trim(default_soil_init_varname)
      IF (present(mj79_relative_factor)) &
         isotope_physics(idx)%mj79_relative_factor = mj79_relative_factor
      IF (present(alpha_liq_vap_fn)) isotope_physics(idx)%alpha_liq_vap => alpha_liq_vap_fn
      IF (present(alpha_ice_vap_fn)) isotope_physics(idx)%alpha_ice_vap => alpha_ice_vap_fn
      IF (present(diffusivity_ratio_air_fn)) &
         isotope_physics(idx)%diffusivity_ratio_air => diffusivity_ratio_air_fn
      IF (present(leaf_liquid_diffusivity_fn)) &
         isotope_physics(idx)%leaf_liquid_diffusivity => leaf_liquid_diffusivity_fn
   END SUBROUTINE register_isotope_physics

   SUBROUTINE ensure_isotope_physics_capacity (need)
      integer, intent(in) :: need
      TYPE(isotope_physics_type), allocatable :: grown(:)
      integer :: old_size, new_size

      IF (.not. allocated(isotope_physics)) THEN
         new_size = max(need, 2)
         allocate(isotope_physics(new_size))
         RETURN
      ENDIF
      old_size = size(isotope_physics)
      IF (need <= old_size) RETURN

      new_size = max(need, old_size * 2)
      allocate(grown(new_size))
      IF (n_isotope_physics > 0) grown(1:n_isotope_physics) = isotope_physics(1:n_isotope_physics)
      CALL move_alloc(grown, isotope_physics)
   END SUBROUTINE ensure_isotope_physics_capacity

   integer FUNCTION find_registered_isotope_by_name (name)
      character(len=*), intent(in) :: name
      integer :: i

      find_registered_isotope_by_name = 0
      IF (.not. allocated(isotope_physics)) RETURN
      DO i = 1, n_isotope_physics
         IF (trim(material_lower(isotope_physics(i)%name)) == trim(material_lower(name))) THEN
            find_registered_isotope_by_name = i
            RETURN
         ENDIF
      ENDDO
   END FUNCTION find_registered_isotope_by_name

   integer FUNCTION find_isotope_physics (imt)
      integer, intent(in) :: imt
      character(len=32) :: lname
      integer :: i

      find_isotope_physics = 0
      IF (imt < 1 .or. imt > nmts) RETURN
      IF (.not. allocated(materials)) RETURN
      IF (.not. material_is_isotope(imt)) RETURN
      IF (.not. allocated(isotope_physics)) RETURN

      lname = material_lower(trim(materials(imt)%name))
      DO i = 1, n_isotope_physics
         IF (isotope_name_matches(lname, isotope_physics(i)%name_patterns)) THEN
            find_isotope_physics = i
            RETURN
         ENDIF
      ENDDO

      IF (materials(imt)%ref_ratio > 0._r8) THEN
         DO i = 1, n_isotope_physics
            IF (isotope_physics(i)%ref_ratio_hint <= mt_tiny) CYCLE
            IF (abs(materials(imt)%ref_ratio - isotope_physics(i)%ref_ratio_hint) / &
                isotope_physics(i)%ref_ratio_hint < isotope_physics(i)%ref_ratio_tolerance) THEN
               find_isotope_physics = i
               RETURN
            ENDIF
         ENDDO
      ENDIF
   END FUNCTION find_isotope_physics

   logical FUNCTION isotope_fractionation_registered (imt)
      integer, intent(in) :: imt

      isotope_fractionation_registered = find_isotope_physics(imt) > 0
   END FUNCTION isotope_fractionation_registered

   real(r8) FUNCTION isotope_alpha_liq_vap (imt, temp_k)
      integer, intent(in) :: imt
      real(r8), intent(in) :: temp_k
      integer :: idx

      isotope_alpha_liq_vap = 1._r8
      idx = find_isotope_physics(imt)
      IF (idx <= 0) RETURN
      IF (associated(isotope_physics(idx)%alpha_liq_vap)) &
         isotope_alpha_liq_vap = isotope_physics(idx)%alpha_liq_vap(temp_k)
   END FUNCTION isotope_alpha_liq_vap

   real(r8) FUNCTION isotope_alpha_ice_vap (imt, temp_k)
      integer, intent(in) :: imt
      real(r8), intent(in) :: temp_k
      integer :: idx

      isotope_alpha_ice_vap = 1._r8
      idx = find_isotope_physics(imt)
      IF (idx <= 0) RETURN
      IF (associated(isotope_physics(idx)%alpha_ice_vap)) &
         isotope_alpha_ice_vap = isotope_physics(idx)%alpha_ice_vap(temp_k)
   END FUNCTION isotope_alpha_ice_vap

   real(r8) FUNCTION isotope_diffusivity_ratio_air (imt)
      integer, intent(in) :: imt
      integer :: idx

      isotope_diffusivity_ratio_air = 1._r8
      idx = find_isotope_physics(imt)
      IF (idx <= 0) RETURN
      IF (associated(isotope_physics(idx)%diffusivity_ratio_air)) &
         isotope_diffusivity_ratio_air = isotope_physics(idx)%diffusivity_ratio_air()
   END FUNCTION isotope_diffusivity_ratio_air

   real(r8) FUNCTION isotope_mj79_relative_factor (imt)
      integer, intent(in) :: imt
      integer :: idx

      isotope_mj79_relative_factor = 1._r8
      idx = find_isotope_physics(imt)
      IF (idx > 0) isotope_mj79_relative_factor = isotope_physics(idx)%mj79_relative_factor
   END FUNCTION isotope_mj79_relative_factor

   real(r8) FUNCTION isotope_leaf_liquid_diffusivity (imt, temp_k)
      integer, intent(in) :: imt
      real(r8), intent(in) :: temp_k
      integer :: idx

      isotope_leaf_liquid_diffusivity = 0._r8
      idx = find_isotope_physics(imt)
      IF (idx <= 0) RETURN
      IF (associated(isotope_physics(idx)%leaf_liquid_diffusivity)) &
         isotope_leaf_liquid_diffusivity = isotope_physics(idx)%leaf_liquid_diffusivity(temp_k)
   END FUNCTION isotope_leaf_liquid_diffusivity

   SUBROUTINE isotope_default_soil_init_varname (imt, varname)
      integer, intent(in) :: imt
      character(len=*), intent(out) :: varname
      integer :: idx

      varname = 'null'
      idx = find_isotope_physics(imt)
      IF (idx > 0) varname = trim(isotope_physics(idx)%default_soil_init_varname)
   END SUBROUTINE isotope_default_soil_init_varname

   logical FUNCTION isotope_name_matches (lname, patterns)
      character(len=*), intent(in) :: lname
      character(len=*), intent(in) :: patterns
      character(len=isotope_patterns_len) :: rest
      character(len=32) :: token
      integer :: comma

      isotope_name_matches = .false.
      rest = material_lower(trim(patterns))
      DO WHILE (len_trim(rest) > 0)
         comma = index(rest, ',')
         IF (comma > 0) THEN
            token = adjustl(rest(1:comma-1))
            rest = adjustl(rest(comma+1:))
         ELSE
            token = adjustl(rest)
            rest = ''
         ENDIF
         IF (len_trim(token) > 0) THEN
            IF (token(1:1) == '=') THEN
               IF (trim(lname) == trim(token(2:))) THEN
                  isotope_name_matches = .true.
                  RETURN
               ENDIF
            ELSEIF (index(lname, trim(token)) > 0) THEN
               isotope_name_matches = .true.
               RETURN
            ENDIF
         ENDIF
      ENDDO
   END FUNCTION isotope_name_matches


END MODULE MOD_MT_Isotope_Registry
#endif
