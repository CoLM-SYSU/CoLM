#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_Solute
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   Code of the solute family (FAMILY_SOLUTE): the per-solute parameter file
!   (&nl_colm_solute_parameter), the dissolved limit with the excess kept as
!   solid, and the first-order decay.  The transport itself is the shared
!   water-carried code (GENERIC_WATER).
!-------------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_MT_Namelist, only: mt_namelist_group_present
   USE MOD_MT_Defs, only: nmts, materials, material_is_solute, material_is_nonvolatile_solute, &
      material_param_file_for_index, material_lower, mt_water_min_for_ratio, REACTION_FIRST_ORDER
   IMPLICIT NONE
   SAVE
   PRIVATE

   ! Fields of &nl_colm_solute_parameter (DEF_SOLUTE%...). Defaults are
   ! filled in from the descriptor before the read.
   type :: solute_parameter_type
      character(len=32) :: unit_kind = 'material_per_water'
      real(r8) :: mol_weight = 0.0_r8
      real(r8) :: init_conc = 0.0_r8
      real(r8) :: precip_default_conc = 0.0_r8
      real(r8) :: vapor_default_conc = 0.0_r8
      real(r8) :: max_dissolved_conc = huge(1.0_r8)
      real(r8) :: decay_rate = 0.0_r8
      ! Indirect N2O emission factor at the land-to-river entry
      ! [kg N2O-N per kg N].  0 disables the conversion for this solute.
      real(r8) :: n2o_ef5 = 0.0_r8
      integer  :: charge = 0
   end type solute_parameter_type

   PUBLIC :: mt_solute_read_parameters
   PUBLIC :: mt_has_dissolved_limit, mt_equilibrate_dissolved
   PUBLIC :: material_decay_fraction

CONTAINS

   SUBROUTINE mt_solute_read_parameters ()
   ! &nl_colm_solute_parameter of each solute with a parameter file
   ! (DEF_SOLUTE_PARAM_FILES); called between material_defs_init and
   ! material_defs_validate.
      IMPLICIT NONE
      integer :: imt
      character(len=512) :: file_param
      logical :: found

      IF (.not. allocated(materials)) RETURN
      DO imt = 1, nmts
         IF (.not. material_is_solute(imt)) CYCLE
         CALL material_param_file_for_index (imt, '', file_param, found)
         IF (found) CALL read_solute_parameter_file (imt, file_param)
      ENDDO
   END SUBROUTINE mt_solute_read_parameters

   SUBROUTINE read_solute_parameter_file (imt, nlfile)
   !-----------------------------------------------------------------------
   !  &nl_colm_solute_parameter of one solute; fields left out keep the
   !  code defaults set in material_defs_init.
   !-----------------------------------------------------------------------
      USE MOD_SPMD_Task, only: p_is_master, CoLM_stop
      IMPLICIT NONE
      integer, intent(in) :: imt
      character(len=*), intent(in) :: nlfile

      type(solute_parameter_type) :: DEF_SOLUTE
      integer :: ierr, unit_nml
      character(len=512) :: iomsg
      namelist /nl_colm_solute_parameter/ DEF_SOLUTE

      IF (imt < 1 .or. imt > nmts) RETURN
      IF (.not. mt_namelist_group_present(nlfile, 'nl_colm_solute_parameter')) RETURN

      DEF_SOLUTE%unit_kind           = materials(imt)%unit_kind
      DEF_SOLUTE%mol_weight          = materials(imt)%mol_weight
      DEF_SOLUTE%init_conc           = materials(imt)%init_conc
      DEF_SOLUTE%precip_default_conc = materials(imt)%precip_default_conc
      DEF_SOLUTE%vapor_default_conc  = materials(imt)%vapor_default_conc
      DEF_SOLUTE%max_dissolved_conc  = materials(imt)%max_dissolved_conc
      DEF_SOLUTE%decay_rate          = materials(imt)%decay_rate
      DEF_SOLUTE%n2o_ef5             = materials(imt)%n2o_ef5
      DEF_SOLUTE%charge              = materials(imt)%charge

      open(newunit=unit_nml, status='OLD', file=trim(nlfile), form='FORMATTED')
      iomsg = ''
      read(unit_nml, nml=nl_colm_solute_parameter, iostat=ierr, iomsg=iomsg)
      close(unit_nml)
      IF (ierr /= 0) THEN
         IF (p_is_master) THEN
            write(*,'(3A)') 'ERROR read_solute_parameter_file: invalid &nl_colm_solute_parameter in ', &
               trim(nlfile), ': '
            write(*,'(A)') trim(iomsg)
         ENDIF
         CALL CoLM_stop()
      ENDIF

      materials(imt)%unit_kind           = material_lower(DEF_SOLUTE%unit_kind)
      materials(imt)%mol_weight          = DEF_SOLUTE%mol_weight
      materials(imt)%init_conc           = DEF_SOLUTE%init_conc
      materials(imt)%precip_default_conc = DEF_SOLUTE%precip_default_conc
      materials(imt)%vapor_default_conc  = DEF_SOLUTE%vapor_default_conc
      materials(imt)%max_dissolved_conc  = DEF_SOLUTE%max_dissolved_conc
      materials(imt)%decay_rate = DEF_SOLUTE%decay_rate
      materials(imt)%n2o_ef5            = DEF_SOLUTE%n2o_ef5
      materials(imt)%charge              = DEF_SOLUTE%charge

      ! The entry-point N2O factor is a fraction of the incoming mass; keep
      ! it inside [0,1) so a typo cannot create or destroy more nitrogen than
      ! enters the river.
      IF (materials(imt)%n2o_ef5 < 0._r8 .or. materials(imt)%n2o_ef5 >= 1._r8) &
         CALL CoLM_stop ('DEF_SOLUTE%n2o_ef5 must be in [0,1).')
   END SUBROUTINE read_solute_parameter_file

   logical FUNCTION mt_has_dissolved_limit (imt)
      integer, intent(in) :: imt

      mt_has_dissolved_limit = .false.
      IF (.not. allocated(materials)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN
      mt_has_dissolved_limit = material_is_nonvolatile_solute(imt) .and. &
         materials(imt)%max_dissolved_conc < huge(1.0_r8)
   END FUNCTION mt_has_dissolved_limit

   SUBROUTINE mt_equilibrate_dissolved (imt, water_mass, dissolved_mass, solid_mass)
      integer,  intent(in)    :: imt
      real(r8), intent(in)    :: water_mass
      real(r8), intent(inout) :: dissolved_mass, solid_mass
      real(r8) :: capacity, total_mass

      IF (.not. mt_has_dissolved_limit(imt)) RETURN
      ! Signed aquifer debt is accounting state, not a precipitating solution.
      IF (dissolved_mass < 0._r8 .or. solid_mass < 0._r8) RETURN

      total_mass = dissolved_mass + solid_mass
      capacity = 0._r8
      IF (water_mass > mt_water_min_for_ratio) THEN
         capacity = materials(imt)%max_dissolved_conc * water_mass
      ENDIF
      dissolved_mass = min(total_mass, capacity)
      solid_mass = total_mass - dissolved_mass
   END SUBROUTINE mt_equilibrate_dissolved

   real(r8) FUNCTION material_decay_fraction (imt, deltim)
      integer,  intent(in) :: imt
      real(r8), intent(in) :: deltim
      real(r8) :: kdt

      material_decay_fraction = 0._r8
      IF (.not. allocated(materials)) RETURN
      IF (imt < 1 .or. imt > nmts) RETURN
      IF (materials(imt)%reaction_mode /= REACTION_FIRST_ORDER) RETURN
      IF (deltim <= 0._r8) RETURN
      IF (materials(imt)%decay_rate <= 0._r8) RETURN

      kdt = min(materials(imt)%decay_rate * deltim, 700._r8)
      material_decay_fraction = 1._r8 - exp(-kdt)
   END FUNCTION material_decay_fraction

END MODULE MOD_MT_Solute
#endif
