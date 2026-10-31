#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_ForcingInput

   ! ------------------------------------------------------------------
   ! Per-species material forcing input.
   !
   ! Each material's parameter file (DEF_ISOTOPE_PARAM_FILE for isotopes,
   ! the DEF_SOLUTE_PARAM_FILES entry for a solute) may carry optional
   ! namelist groups describing that species' OWN forcing inputs:
   !
   !   &nl_colm_mt_forcing
   !      forcing_species    = 'H2_18O'
   !      forcing_num        = 2
   !      forcing_role       = 'precip', 'vapor'
   !      forcing_fprefix    = 'IsoGSM_prate', 'IsoGSM_Q'
   !      forcing_vname      = 'prate1sfc',    'spfh12m'
   !      forcing_tintalgo   = 'nearest',      'linear'
   !      forcing_dtime      = 21600, 21600
   !      forcing_offset     = 10800, 10800
   !      forcing_input_mode = 'normalized_over_total', 'normalized_over_total'
   !   /
   !
   ! The isotope file is shared by all isotopes, so it holds one group per
   ! species and each must name forcing_species. In a solute's own file
   ! forcing_species may be left out.
   !
   ! `role` names a CONSUMER, and only two exist: 'precip' and 'vapor', both
   ! read by MOD_MT_Forcing. A role is accepted only if something reads it.
   !
   ! 'inundation', 'atm' and 'erosion_yield' were once accepted here for gas and
   ! particle species. Nothing ever read them, so those species could be given
   ! only roles that did nothing: the forcing file was opened, validated, stored
   ! in the spec and then silently ignored for the whole run. They are now
   ! refused with a message naming where that field actually comes from --
   ! DEF_file_GIEMS for inundation, DEF_METHANE%atm_methane (or the transient
   ! table) for atmospheric CH4, and the internal precipitation-driven
   ! calculation in MOD_MT_Particle_Sediment for erosion yield.
   !
   ! This module replaces the former global DEF_forcing%material_* / legacy
   ! per-species DEF_forcing%precipitation_O18_* namelist path: material
   ! forcing now lives entirely in the per-species parameter files and is
   ! read by the MATERIAL subsystem itself.
   ! ------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_MT_Defs, only: nmts, materials, material_param_file_for_index, material_lower, &
                             mt_uses_generic_water, material_is_isotope

   IMPLICIT NONE
   SAVE
   PRIVATE

   integer, parameter, public :: MATERIAL_FORCING_MAX = 8

   type, public :: material_forcing_spec_type
      character(len=32)  :: role       = 'none'
      character(len=256) :: fprefix    = 'null'
      character(len=64)  :: vname      = 'null'
      character(len=32)  :: tintalgo   = 'linear'
      integer            :: dtime      = 21600
      integer            :: offset     = 0
      character(len=32)  :: input_mode = 'normalized_over_total'
   end type material_forcing_spec_type

   integer, allocatable, public :: material_forcing_n(:)                 ! (nmts)
   type(material_forcing_spec_type), allocatable, public :: &
        material_forcing_specs(:,:)                                      ! (MATERIAL_FORCING_MAX, nmts)

   PUBLIC :: material_forcing_input_load
   PUBLIC :: material_forcing_input_final
   PUBLIC :: material_forcing_input_count
   PUBLIC :: material_forcing_input_get
   PUBLIC :: material_forcing_input_find

CONTAINS

   SUBROUTINE material_forcing_input_load ()
      USE MOD_SPMD_Task, only: p_is_master, CoLM_stop
      IMPLICIT NONE
      integer :: imt, k, kdup, nf, ierr, unit_nml
      logical :: found, fexists
      character(len=256) :: nlfile
      character(len=256) :: role_hint
      character(len=512) :: iomsg

      ! namelist scratch (parallel arrays for cross-compiler robustness)
      character(len=32)  :: forcing_species
      integer            :: forcing_num
      character(len=32)  :: forcing_role(MATERIAL_FORCING_MAX)
      character(len=256) :: forcing_fprefix(MATERIAL_FORCING_MAX)
      character(len=64)  :: forcing_vname(MATERIAL_FORCING_MAX)
      character(len=32)  :: forcing_tintalgo(MATERIAL_FORCING_MAX)
      integer            :: forcing_dtime(MATERIAL_FORCING_MAX)
      integer            :: forcing_offset(MATERIAL_FORCING_MAX)
      character(len=32)  :: forcing_input_mode(MATERIAL_FORCING_MAX)
      integer            :: nmatch, ipass
      namelist /nl_colm_mt_forcing/ forcing_species, forcing_num, forcing_role, forcing_fprefix, &
         forcing_vname, forcing_tintalgo, forcing_dtime, forcing_offset, forcing_input_mode

      IF (nmts <= 0) RETURN
      IF (.not. allocated(material_forcing_n))     allocate(material_forcing_n(nmts))
      IF (.not. allocated(material_forcing_specs)) allocate(material_forcing_specs(MATERIAL_FORCING_MAX, nmts))
      material_forcing_n(:) = 0

      DO imt = 1, nmts
         CALL material_param_file_for_index (imt, '', nlfile, found)
         IF (.not. found) CYCLE
         INQUIRE (file=trim(nlfile), exist=fexists)
         IF (.not. fexists) CYCLE

         ! Read every &nl_colm_mt_forcing group of the file in turn. The
         ! first pass counts the groups of this species (more than one is an
         ! error: only one could ever be used); the second keeps it.
         nmatch = 0
         DO ipass = 1, 2
            IF (ipass == 2 .and. nmatch == 0) EXIT
            open (newunit=unit_nml, status='OLD', file=trim(nlfile), form='FORMATTED')
            DO
               CALL reset_scratch ()
               iomsg = ''
               read (unit_nml, nml=nl_colm_mt_forcing, iostat=ierr, iomsg=iomsg)
               IF (ierr < 0) EXIT
               IF (ierr > 0) THEN
                  IF (p_is_master) THEN
                     write(*,'(3A)') 'ERROR material_forcing_input_load: invalid &nl_colm_mt_forcing in ', &
                        trim(nlfile), ': '
                     write(*,'(A)') trim(iomsg)
                  ENDIF
                  CALL CoLM_stop()
               ENDIF
               IF (len_trim(forcing_species) == 0) THEN
                  IF (material_is_isotope(imt)) CALL CoLM_stop ( &
                     'forcing_species must be set in each &nl_colm_mt_forcing of '//trim(nlfile))
               ELSEIF (trim(material_lower(forcing_species)) /= trim(material_lower(materials(imt)%name))) THEN
                  CYCLE
               ENDIF
               IF (ipass == 2) EXIT
               nmatch = nmatch + 1
            ENDDO
            close (unit_nml)
         ENDDO
         IF (nmatch == 0) CYCLE
         IF (nmatch > 1) CALL CoLM_stop ('More than one &nl_colm_mt_forcing group for '// &
            trim(materials(imt)%name)//' in '//trim(nlfile))

         IF (forcing_num < 0 .or. forcing_num > MATERIAL_FORCING_MAX) THEN
            IF (p_is_master) THEN
               write(*,'(A,I0,A,I0,A,A,A)') 'ERROR material_forcing_input_load: forcing_num=', forcing_num, &
                  ' must be between 0 and ', MATERIAL_FORCING_MAX, ' in ', trim(nlfile), '.'
            ENDIF
            CALL CoLM_stop()
         ENDIF
         nf = forcing_num
         material_forcing_n(imt) = nf
         DO k = 1, nf
            IF (forcing_dtime(k) <= 0) THEN
               IF (p_is_master) THEN
                  write(*,'(A,I0,A,A,A,I0,A)') 'ERROR material_forcing_input_load: forcing_dtime(', k, &
                     ') for material "', trim(materials(imt)%name), '" must be > 0, got ', forcing_dtime(k), '.'
               ENDIF
               CALL CoLM_stop()
            ENDIF
            forcing_tintalgo(k) = material_lower(adjustl(forcing_tintalgo(k)))
            SELECT CASE (trim(forcing_tintalgo(k)))
            CASE ('linear', 'nearest')
            CASE DEFAULT
               IF (p_is_master) WRITE(*,'(A,I0,5A)') &
                  'ERROR material_forcing_input_load: forcing_tintalgo(', k, ') for material "', &
                  trim(materials(imt)%name), '" has invalid value "', trim(forcing_tintalgo(k)), '".'
               CALL CoLM_stop()
            END SELECT
            forcing_input_mode(k) = material_lower(adjustl(forcing_input_mode(k)))
            SELECT CASE (trim(forcing_input_mode(k)))
            CASE ('direct', 'value', 'raw', 'delta', &
                  'heavy_over_total', 'ratio_to_total', 'over_total', &
                  'normalized_over_total', 'normalized_ratio', 'isogsm', 'standard_over_total')
            CASE DEFAULT
               IF (p_is_master) WRITE(*,'(A,I0,5A)') &
                  'ERROR material_forcing_input_load: forcing_input_mode(', k, ') for material "', &
                  trim(materials(imt)%name), '" has invalid value "', trim(forcing_input_mode(k)), '".'
               CALL CoLM_stop()
            END SELECT
            forcing_role(k) = material_lower(adjustl(forcing_role(k)))
            IF (.not. material_forcing_role_valid(imt, forcing_role(k))) THEN
               role_hint = material_forcing_role_hint (imt, forcing_role(k))
               IF (p_is_master) THEN
                  WRITE(*,'(A,I0,5A)') &
                     'ERROR material_forcing_input_load: forcing_role(', k, ') for material "', &
                     trim(materials(imt)%name), '" has invalid value "', trim(forcing_role(k)), '".'
                  IF (len_trim(role_hint) > 0) THEN
                     ! Accepting it and reading nothing is what this replaces.
                     WRITE(*,'(2A)') '       no consumer reads this role: ', trim(role_hint)
                  ELSE
                     WRITE(*,'(A)')  '       supported roles are "precip" and "vapor".'
                  ENDIF
               ENDIF
               CALL CoLM_stop()
            ENDIF

            ! material_forcing_input_find returns the FIRST match, so a second
            ! entry with the same role is stored, passes every check above and
            ! is then never read -- the same silent no-op the role check just
            ! closed, arrived at from the other direction. A user who lists two
            ! precip files means one of them to be used.
            DO kdup = 1, k-1
               IF (trim(forcing_role(kdup)) == trim(forcing_role(k))) THEN
                  IF (p_is_master) THEN
                     WRITE(*,'(A,I0,A,I0,4A)') &
                        'ERROR material_forcing_input_load: forcing_role(', kdup, ') and (', k, &
                        ') for material "', trim(materials(imt)%name), '" are both "', &
                        trim(forcing_role(k))//'".'
                     WRITE(*,'(A)') '       only the first would ever be read; each role may appear once.'
                  ENDIF
                  CALL CoLM_stop()
               ENDIF
            ENDDO

            material_forcing_specs(k,imt)%role       = adjustl(forcing_role(k))
            material_forcing_specs(k,imt)%fprefix    = adjustl(forcing_fprefix(k))
            material_forcing_specs(k,imt)%vname      = adjustl(forcing_vname(k))
            material_forcing_specs(k,imt)%tintalgo   = forcing_tintalgo(k)
            material_forcing_specs(k,imt)%dtime      = forcing_dtime(k)
            material_forcing_specs(k,imt)%offset     = forcing_offset(k)
            material_forcing_specs(k,imt)%input_mode = forcing_input_mode(k)
         ENDDO

         IF (p_is_master) write(*,'(A,I0,A,A,A,I0,A)') &
            ' [material-forcing] material ', imt, ' (', trim(materials(imt)%name), '): ', nf, ' forcing input(s)'
      ENDDO
   CONTAINS

      SUBROUTINE reset_scratch ()
         forcing_species        = ''
         forcing_num            = 0
         forcing_role(:)        = 'none'
         forcing_fprefix(:)     = 'null'
         forcing_vname(:)       = 'null'
         forcing_tintalgo(:)    = 'linear'
         forcing_dtime(:)       = 21600
         forcing_offset(:)      = 0
         forcing_input_mode(:)  = 'normalized_over_total'
      END SUBROUTINE reset_scratch

   END SUBROUTINE material_forcing_input_load

   integer FUNCTION material_forcing_input_count (imt)
      IMPLICIT NONE
      integer, intent(in) :: imt
      material_forcing_input_count = 0
      IF (.not. allocated(material_forcing_n)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN
      material_forcing_input_count = material_forcing_n(imt)
   END FUNCTION material_forcing_input_count

   FUNCTION material_forcing_input_get (imt, k) RESULT(spec)
      IMPLICIT NONE
      integer, intent(in) :: imt, k
      type(material_forcing_spec_type) :: spec
      ! defaults from the type definition apply if out of range
      IF (.not. allocated(material_forcing_specs)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN
      IF (k < 1 .or. k > MATERIAL_FORCING_MAX) RETURN
      spec = material_forcing_specs(k,imt)
   END FUNCTION material_forcing_input_get

   integer FUNCTION material_forcing_input_find (imt, role)
      IMPLICIT NONE
      integer, intent(in) :: imt
      character(len=*), intent(in) :: role
      integer :: k
      material_forcing_input_find = 0
      IF (.not. allocated(material_forcing_n)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN
      DO k = 1, material_forcing_n(imt)
         IF (material_forcing_role_equal(material_forcing_specs(k,imt)%role, role)) THEN
            material_forcing_input_find = k
            RETURN
         ENDIF
      ENDDO
   END FUNCTION material_forcing_input_find

   logical FUNCTION material_forcing_role_equal (a, b)
      IMPLICIT NONE
      character(len=*), intent(in) :: a, b
      material_forcing_role_equal = (material_lower(a) == material_lower(b))
   END FUNCTION material_forcing_role_equal

   logical FUNCTION material_forcing_role_valid (imt, role)
      IMPLICIT NONE
      integer, intent(in) :: imt
      character(len=*), intent(in) :: role

      ! Validity means "some code reads this", not "this name is spelled
      ! correctly for this category". The only consumers are the two
      ! material_forcing_input_find calls in MOD_MT_Forcing, and that loop
      ! skips any material without generic land-water transport -- so the same
      ! role is live for one species and dead for another, and the category
      ! alone cannot tell them apart.
      material_forcing_role_valid = &
         (trim(role) == 'precip' .or. trim(role) == 'vapor') .and. &
         mt_uses_generic_water (imt)

   END FUNCTION material_forcing_role_valid

   !> Explains a role that names a real field but has no consumer here, so the
   !! error can point at the mechanism that does supply it. Empty when the role
   !! is simply unrecognised -- that is a typo, not a misdirected configuration.
   FUNCTION material_forcing_role_hint (imt, role) RESULT (hint)
      IMPLICIT NONE
      integer, intent(in) :: imt
      character(len=*), intent(in) :: role
      character(len=256) :: hint

      ! A supported role on a species the consumer loop skips. The name is
      ! right and the file would load, but MOD_MT_Forcing never reaches
      ! this material, so the data would go nowhere.
      IF ((trim(role) == 'precip' .or. trim(role) == 'vapor') .and. &
          .not. mt_uses_generic_water (imt)) THEN
         hint = 'this species is provider-owned and is skipped by the ' // &
                'precip/vapor loop in MOD_MT_Forcing'
         RETURN
      ENDIF

      SELECT CASE (trim(role))
      CASE ('inundation')
         hint = 'inundated fraction is read from DEF_file_GIEMS, not from material forcing'
      CASE ('atm')
         hint = 'atmospheric CH4 comes from DEF_METHANE%atm_methane or its transient table'
      CASE ('erosion_yield')
         hint = 'hillslope erosion is computed from precipitation in MOD_MT_Particle_Sediment'
      CASE DEFAULT
         hint = ''
      END SELECT

   END FUNCTION material_forcing_role_hint

   SUBROUTINE material_forcing_input_final ()
      IMPLICIT NONE
      IF (allocated(material_forcing_n))     deallocate(material_forcing_n)
      IF (allocated(material_forcing_specs)) deallocate(material_forcing_specs)
   END SUBROUTINE material_forcing_input_final

END MODULE MOD_MT_ForcingInput
#endif
