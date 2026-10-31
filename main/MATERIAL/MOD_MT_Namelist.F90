#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_Namelist
!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!  Switches and settings of material transport.
!
!  They are kept out of MOD_Namelist on purpose: the main namelist file
!  carries one extra group, read here on the master and broadcast,
!
!     &nl_colm_material
!        DEF_USE_ISOTOPE        = .true.
!        DEF_ISOTOPE_SPECIES    = 'H2_18O,HDO'
!        DEF_ISOTOPE_PARAM_FILE = 'run/standard_isotope_parameter.nml'
!        DEF_SOLUTE_NAMES       = 'CL'
!        DEF_SOLUTE_PARAM_FILES = 'CL:run/standard_chloride_parameter.nml'
!        DEF_SEDIMENT_PARAM_FILE = 'run/standard_sediment_parameter.nml'
!        DEF_USE_METHANE        = .true.
!        DEF_METHANE_PARAM_FILE = 'run/standard_methane_parameter.nml'
!     /
!
!  The group may be left out; every setting then keeps its default and no
!  material is carried.  DEF_USE_SEDIMENT lives in &nl_colm_material as
!  well (it is owned by MOD_Namelist so the river code can read it), and
!  the river sediment needs MaterialTransport plus GridRiverLakeFlow; its
!  &nl_colm_sediment_parameter (DEF_SEDIMENT_PARAM_FILE, 'null' = code
!  defaults) is read by MOD_MT_Sediment.
!
!  Isotope physics settings live in &nl_colm_isotope_parameter of
!  DEF_ISOTOPE_PARAM_FILE ('null' = code defaults). Solute parameter files
!  are read per species by MOD_MT_Defs; per-species forcing by
!  MOD_MT_ForcingInput.  Methane (needs BGC and PFT/PC) reads its settings from
!  &nl_colm_methane_parameter of DEF_METHANE_PARAM_FILE (MOD_MT_Methane_Const).
!-----------------------------------------------------------------------

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

   ! ----- &nl_colm_material -----
   logical            :: DEF_USE_ISOTOPE        = .false.
   ! Only species with compiled physics: H2_18O, HDO.
   character(len=256) :: DEF_ISOTOPE_SPECIES    = 'H2_18O,HDO'
   character(len=256) :: DEF_ISOTOPE_PARAM_FILE = 'null'
   character(len=256) :: DEF_SOLUTE_NAMES       = ''
   ! Species-resolved NH4 leaching in BGC (MOD_BGC_Soil_BiogeochemNLeaching).
   ! 'auto' (default) turns it on exactly when DEF_SOLUTE_NAMES contains NH4;
   ! 'on'/'off' force the choice.  Only meaningful with BGC.  Off leaves the
   ! upstream BGC results bit-identical.
   character(len=16)  :: DEF_NH4_LEACHING = 'auto'
   logical            :: DEF_USE_NH4_LEACHING  = .false.
   ! 'NAME:path' entries separated by ',' or ';', or plain paths in the
   ! order of DEF_SOLUTE_NAMES. A solute without a file uses code defaults.
   character(len=2048):: DEF_SOLUTE_PARAM_FILES = 'null'
   character(len=256) :: DEF_SEDIMENT_PARAM_FILE = 'null'
   ! Soil, wetland and lake methane (FAMILY_GAS); needs BGC and a parameter file.
   logical            :: DEF_USE_METHANE        = .false.
   character(len=256) :: DEF_METHANE_PARAM_FILE = 'null'

   ! Stop when more than this many cells fail a check (default 0: any failure stops).
   integer            :: DEF_MT_BALANCE_ABORT_NBAD = 0
   integer            :: DEF_MT_RESID_ABORT_NBAD   = 0
   integer            :: DEF_MT_LULCC_ABORT_NBAD   = 0
   ! Water depth that carries the aquifer signature when the
   ! variably saturated flow scheme is on (required > 0 for isotopes then).
   real(r8)           :: DEF_MT_AQUIFER_MIXING_WATER_MM = -1._r8
   ! Total precipitation rate [kg m-2 s-1] at or below which a heavy/total
   ! forcing ratio is taken as near-dry and not used.
   ! IsoGSM 2004 global (January and July, 6-hourly): at 1e-7..1e-6 the
   ! ratio is noise (60-65% of records with |d-excess| > 100), 22-32% at
   ! 1e-6..2e-6, 3-8% at 2e-6..5e-6, at most 0.5% above 5e-6; 5e-6 leaves
   ! 1.2-1.5% of the precipitation mass to the fallback.  For isotopes with
   ! a vapour forcing the fallback is rain in equilibrium with the vapour and
   ! snow at the last reliable ratio (material_forcing_rain_value /
   ! material_forcing_snow_value in MOD_MT_Forcing); otherwise the last valid
   ! ratio (the default ratio before the first one).
   real(r8)           :: DEF_MT_FORCING_MIN_PRECIP = 5.e-6_r8
   ! Total precipitation rate [kg m-2 s-1] above which a record is remembered
   ! as the snow fallback ratio (isotopes with a vapour forcing).  Records just
   ! above DEF_MT_FORCING_MIN_PRECIP still carry noisy d-excess that a
   ! remembered ratio would repeat for days.  IsoGSM 2004, predicting the
   ! ratio of the next snowfall above 2e-5 (d18O RMSE / d-excess RMSE,
   ! permil): remembering records above 5e-6: NEC Jan 4.0 / 14.0, 40-60N Jan
   ! 4.4 / 12.2, 40-60N May 4.4 / 20.9; above 2e-5: 4.5 / 9.2, 4.5 / 8.4,
   ! 4.3 / 17.6; above 5e-5 the remembered ratio gets old (NEC Jan bias +4.5,
   ! RMSE 7.2).  Must be >= DEF_MT_FORCING_MIN_PRECIP when phase fill is used.
   real(r8)           :: DEF_MT_FORCING_MEMORY_PRECIP = 2.e-5_r8
   ! Indirect N2O emission factors for the nitrogen solutes at the land-to-
   ! river entry [kg N2O-N per kg N entering the river].  The default -1
   ! means "not set": the factor then comes from the solute parameter file
   ! (DEF_SOLUTE%n2o_ef5).  Any value >= 0 overrides it, so 0 disables the
   ! conversion for that species, which leaves every result bit-identical
   ! to a build without this option (a 0 default would be indistinguishable
   ! from "unset" and could never switch the conversion off).
   !
   ! Reference values (IPCC 2019 Refinement, Vol. 4 Ch. 11, Table 11.3):
   !   EF5g = 0.006  groundwater / surface drainage
   !   EF5r = 0.0025 rivers
   !   EF5e = 0.0025 estuaries            total EF5 = 0.011
   ! The 2006 Guidelines used a single 0.0075 (each component 0.0025).
   ! Tian, Cai & Akiyama (2019), Environ. Pollut. 245:300-306 (254
   ! measurements): groundwater/springs 0.0079, surface drainage 0.0040,
   ! lakes/ponds/reservoirs 0.0012, rivers 0.0030.  Maavara et al. (2019),
   ! Glob. Change Biol. 25:473-488 model nitrification/denitrification with
   ! residence time and conclude the IPCC factors may overestimate inland-
   ! water N2O by up to an order of magnitude; river N2O is mostly
   ! denitrification-derived while reservoirs/estuaries are nitrification-
   ! dominated, which a single factor cannot represent.
   !
   ! Method caveat: this is the inventory arithmetic (N_leached x EF5)
   ! applied once at the river entry, not a mechanistic water-body model.
   ! It does not represent where in the network the emission happens or how
   ! long the water resides, and it must not be applied again along the
   ! river.  Suggested use: start with the river component (0.0025) and use
   ! 0.006 / 0.011 as sensitivity bounds.
   real(r8)           :: DEF_MT_N2O_EF5_NH4 = -1._r8
   real(r8)           :: DEF_MT_N2O_EF5_NO3 = -1._r8
   logical            :: DEF_MT_USE_SOIL_INIT   = .false.
   character(len=256) :: DEF_MT_SOIL_INIT_FILE  = 'null'
   character(len=256) :: DEF_MT_SOIL_INIT_VARS  = ''

   ! ----- &nl_colm_isotope_parameter (in DEF_ISOTOPE_PARAM_FILE) -----
   logical            :: DEF_ISO_USE_FRACTIONATION     = .false.
   character(len=16)  :: DEF_ISO_KINETIC_SCHEME        = 'MERLIVAT1978'
   real(r8)           :: DEF_ISO_ICE_SUPERSAT_SLOPE    = 0.003_r8
   real(r8)           :: DEF_ISO_CG_RELHUM_MAX         = 0.99_r8
   character(len=16)  :: DEF_ISO_OPEN_WATER_KINETIC    = 'MJ79'
   real(r8)           :: DEF_ISO_SUBL_SKIN_MM          = 5.0_r8
   character(len=16)  :: DEF_ISO_SOIL_KINETIC          = 'RESISTANCE'
   logical            :: DEF_ISO_SOIL_DIFFUSION        = .true.
   logical            :: DEF_ISO_SOIL_VAPOR_DIFFUSION  = .true.
   real(r8)           :: DEF_ISO_CANOPY_EQUILIBRATION  = 0.0_r8
   real(r8)           :: DEF_ISO_SNOWMELT_EQUILIBRATION= 0.0_r8
   real(r8)           :: DEF_ISO_NSS_LEAF_WATER_PER_LAI= 0.12_r8
   real(r8)           :: DEF_ISO_NSS_LEAF_PATH_LENGTH  = 0.01_r8
   real(r8)           :: DEF_ISO_NSS_LEAF_RB           = 100._r8
   ! Per-species values; defaults are VSMOW and the former
   ! standard_O18/HDO_parameter.nml.
   real(r8)           :: DEF_ISO_REF_RATIO_H2_18O      = 2.0052e-3_r8
   real(r8)           :: DEF_ISO_INIT_DELTA_H2_18O     = -10._r8
   real(r8)           :: DEF_ISO_REF_RATIO_HDO         = 1.5576e-4_r8
   real(r8)           :: DEF_ISO_INIT_DELTA_HDO        = -70._r8

   logical, private :: loaded = .false.

   PUBLIC :: mt_read_namelist
   PUBLIC :: mt_any_enabled
   PUBLIC :: mt_namelist_group_present

CONTAINS

   SUBROUTINE mt_read_namelist (nlfile)
   !-----------------------------------------------------------------------
   !  Read the material groups once. nlfile is the main namelist file;
   !  without it the first command-line argument is used, as in CoLM.F90.
   !-----------------------------------------------------------------------
   USE MOD_SPMD_Task
   ! River-sediment switch: owned by MOD_Namelist (shared with the river
   ! code) but read here, from &nl_colm_material, not from &nl_colm.
   USE MOD_Namelist, only: DEF_USE_SEDIMENT
   IMPLICIT NONE

   character(len=*), intent(in), optional :: nlfile

   character(len=256) :: mainfile
   integer :: ierr, unit_nml
   character(len=512) :: iomsg

      namelist /nl_colm_material/ DEF_USE_ISOTOPE, DEF_ISOTOPE_SPECIES, &
         DEF_ISOTOPE_PARAM_FILE, DEF_SOLUTE_NAMES, DEF_SOLUTE_PARAM_FILES, &
         DEF_SEDIMENT_PARAM_FILE, DEF_USE_METHANE, DEF_METHANE_PARAM_FILE, &
         DEF_MT_BALANCE_ABORT_NBAD, DEF_MT_RESID_ABORT_NBAD, DEF_MT_LULCC_ABORT_NBAD, &
         DEF_MT_AQUIFER_MIXING_WATER_MM, DEF_MT_FORCING_MIN_PRECIP, &
         DEF_MT_FORCING_MEMORY_PRECIP, DEF_MT_USE_SOIL_INIT, DEF_MT_SOIL_INIT_FILE, DEF_MT_SOIL_INIT_VARS, &
         DEF_MT_N2O_EF5_NH4, DEF_MT_N2O_EF5_NO3, DEF_USE_SEDIMENT, DEF_NH4_LEACHING

      namelist /nl_colm_isotope_parameter/ DEF_ISO_USE_FRACTIONATION, &
         DEF_ISO_KINETIC_SCHEME, DEF_ISO_ICE_SUPERSAT_SLOPE, DEF_ISO_CG_RELHUM_MAX, &
         DEF_ISO_OPEN_WATER_KINETIC, DEF_ISO_SUBL_SKIN_MM, DEF_ISO_SOIL_KINETIC, &
         DEF_ISO_SOIL_DIFFUSION, DEF_ISO_SOIL_VAPOR_DIFFUSION, &
         DEF_ISO_CANOPY_EQUILIBRATION, DEF_ISO_SNOWMELT_EQUILIBRATION, &
         DEF_ISO_NSS_LEAF_WATER_PER_LAI, DEF_ISO_NSS_LEAF_PATH_LENGTH, DEF_ISO_NSS_LEAF_RB, &
         DEF_ISO_REF_RATIO_H2_18O, DEF_ISO_INIT_DELTA_H2_18O, &
         DEF_ISO_REF_RATIO_HDO, DEF_ISO_INIT_DELTA_HDO

      IF (loaded) RETURN

      IF (p_is_master) THEN

         IF (present(nlfile)) THEN
            mainfile = nlfile
         ELSE
            CALL get_command_argument (1, mainfile)
         ENDIF

         IF (mt_namelist_group_present (mainfile, 'nl_colm_material')) THEN
            open (newunit=unit_nml, status='OLD', file=trim(mainfile), form='FORMATTED')
            iomsg = ''
            read (unit_nml, nml=nl_colm_material, iostat=ierr, iomsg=iomsg)
            close (unit_nml)
            IF (ierr /= 0) THEN
               write(*,*) 'ERROR in ', trim(mainfile), ' : ', trim(iomsg)
               CALL CoLM_Stop (' ***** ERROR: Problem reading &nl_colm_material')
            ENDIF
         ENDIF

         IF (len_trim(DEF_SOLUTE_PARAM_FILES) == len(DEF_SOLUTE_PARAM_FILES)) &
            CALL CoLM_Stop ('DEF_SOLUTE_PARAM_FILES is too long.')

         IF (DEF_USE_ISOTOPE .and. mt_file_given (DEF_ISOTOPE_PARAM_FILE)) THEN
            IF (.not. mt_namelist_group_present (DEF_ISOTOPE_PARAM_FILE, &
                                                 'nl_colm_isotope_parameter')) THEN
               CALL CoLM_Stop ('No &nl_colm_isotope_parameter in '//trim(DEF_ISOTOPE_PARAM_FILE))
            ENDIF
            open (newunit=unit_nml, status='OLD', file=trim(DEF_ISOTOPE_PARAM_FILE), form='FORMATTED')
            iomsg = ''
            read (unit_nml, nml=nl_colm_isotope_parameter, iostat=ierr, iomsg=iomsg)
            close (unit_nml)
            IF (ierr /= 0) THEN
               write(*,*) 'ERROR in ', trim(DEF_ISOTOPE_PARAM_FILE), ' : ', trim(iomsg)
               CALL CoLM_Stop (' ***** ERROR: Problem reading &nl_colm_isotope_parameter')
            ENDIF
         ENDIF

         CALL resolve_nh4_leaching ()
         CALL check_settings ()

      ENDIF

#ifdef USEMPI
      CALL mpi_bcast (DEF_USE_ISOTOPE               ,1   ,mpi_logical   ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISOTOPE_SPECIES           ,256 ,mpi_character ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISOTOPE_PARAM_FILE        ,256 ,mpi_character ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_SOLUTE_NAMES              ,256 ,mpi_character ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_SOLUTE_PARAM_FILES        ,len(DEF_SOLUTE_PARAM_FILES) ,mpi_character ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_SEDIMENT_PARAM_FILE       ,256 ,mpi_character ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_USE_SEDIMENT              ,1   ,mpi_logical   ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_USE_NH4_LEACHING          ,1   ,mpi_logical   ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_NH4_LEACHING              ,16  ,mpi_character ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_USE_METHANE               ,1   ,mpi_logical   ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_METHANE_PARAM_FILE        ,256 ,mpi_character ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_MT_BALANCE_ABORT_NBAD     ,1   ,mpi_integer   ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_MT_RESID_ABORT_NBAD       ,1   ,mpi_integer   ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_MT_LULCC_ABORT_NBAD       ,1   ,mpi_integer   ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_MT_AQUIFER_MIXING_WATER_MM,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_MT_FORCING_MIN_PRECIP     ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_MT_FORCING_MEMORY_PRECIP  ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_MT_N2O_EF5_NH4            ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_MT_N2O_EF5_NO3            ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_MT_USE_SOIL_INIT          ,1   ,mpi_logical   ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_MT_SOIL_INIT_FILE         ,256 ,mpi_character ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_MT_SOIL_INIT_VARS         ,256 ,mpi_character ,p_address_master ,p_comm_glb ,p_err)

      CALL mpi_bcast (DEF_ISO_USE_FRACTIONATION     ,1   ,mpi_logical   ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_KINETIC_SCHEME        ,16  ,mpi_character ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_ICE_SUPERSAT_SLOPE    ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_CG_RELHUM_MAX         ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_OPEN_WATER_KINETIC    ,16  ,mpi_character ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_SUBL_SKIN_MM          ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_SOIL_KINETIC          ,16  ,mpi_character ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_SOIL_DIFFUSION        ,1   ,mpi_logical   ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_SOIL_VAPOR_DIFFUSION  ,1   ,mpi_logical   ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_CANOPY_EQUILIBRATION  ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_SNOWMELT_EQUILIBRATION,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_NSS_LEAF_WATER_PER_LAI,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_NSS_LEAF_PATH_LENGTH  ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_NSS_LEAF_RB           ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_REF_RATIO_H2_18O      ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_INIT_DELTA_H2_18O     ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_REF_RATIO_HDO         ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_INIT_DELTA_HDO        ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
#endif

      loaded = .true.

   END SUBROUTINE mt_read_namelist


   SUBROUTINE check_settings ()

   USE MOD_SPMD_Task, only: CoLM_Stop
   USE MOD_Namelist, only: DEF_USE_SEDIMENT
   IMPLICIT NONE

      ! The switch is read from &nl_colm_material; it is useless without the
      ! grid river-lake flow, so fail loudly here rather than silently
      ! ignoring it.
#ifndef GridRiverLakeFlow
      IF (DEF_USE_SEDIMENT) &
         CALL CoLM_Stop ('DEF_USE_SEDIMENT needs GridRiverLakeFlow in include/define.h.')
#endif

      IF (.not. (DEF_MT_FORCING_MIN_PRECIP >= 0._r8)) &
         CALL CoLM_Stop ('DEF_MT_FORCING_MIN_PRECIP must be >= 0.')
      IF (.not. (DEF_MT_FORCING_MEMORY_PRECIP >= 0._r8)) &
         CALL CoLM_Stop ('DEF_MT_FORCING_MEMORY_PRECIP must be >= 0.')
      ! -1 is the "not set" sentinel (use the solute parameter file); any
      ! other accepted value is a real emission factor in [0,1).
      IF ((DEF_MT_N2O_EF5_NH4 < 0._r8 .and. DEF_MT_N2O_EF5_NH4 /= -1._r8) .or. &
          DEF_MT_N2O_EF5_NH4 >= 1._r8 .or. &
          (DEF_MT_N2O_EF5_NO3 < 0._r8 .and. DEF_MT_N2O_EF5_NO3 /= -1._r8) .or. &
          DEF_MT_N2O_EF5_NO3 >= 1._r8) &
         CALL CoLM_Stop ('DEF_MT_N2O_EF5_NH4/NO3 must be -1 (unset) or in [0,1).')

      SELECT CASE (trim(adjustl(DEF_ISO_KINETIC_SCHEME)))
      CASE ('CAPPA2003', 'Cappa2003', 'cappa2003')
         DEF_ISO_KINETIC_SCHEME = 'CAPPA2003'
      CASE ('MERLIVAT1978', 'Merlivat1978', 'merlivat1978')
         DEF_ISO_KINETIC_SCHEME = 'MERLIVAT1978'
      CASE DEFAULT
         CALL CoLM_Stop ('Invalid DEF_ISO_KINETIC_SCHEME.')
      END SELECT

      SELECT CASE (trim(adjustl(DEF_ISO_OPEN_WATER_KINETIC)))
      CASE ('MJ79', 'mj79')
         DEF_ISO_OPEN_WATER_KINETIC = 'MJ79'
      CASE ('EXPONENT', 'exponent')
         DEF_ISO_OPEN_WATER_KINETIC = 'EXPONENT'
      CASE DEFAULT
         CALL CoLM_Stop ('Invalid DEF_ISO_OPEN_WATER_KINETIC.')
      END SELECT

      SELECT CASE (trim(adjustl(DEF_ISO_SOIL_KINETIC)))
      CASE ('RESISTANCE', 'resistance')
         DEF_ISO_SOIL_KINETIC = 'RESISTANCE'
      CASE ('EXPONENT', 'exponent')
         DEF_ISO_SOIL_KINETIC = 'EXPONENT'
      CASE DEFAULT
         CALL CoLM_Stop ('Invalid DEF_ISO_SOIL_KINETIC.')
      END SELECT

      IF (DEF_ISO_CG_RELHUM_MAX <= 0._r8 .or. DEF_ISO_CG_RELHUM_MAX >= 1._r8 .or. &
          DEF_ISO_CANOPY_EQUILIBRATION < 0._r8 .or. DEF_ISO_CANOPY_EQUILIBRATION > 1._r8 .or. &
          DEF_ISO_SNOWMELT_EQUILIBRATION < 0._r8 .or. DEF_ISO_SNOWMELT_EQUILIBRATION > 1._r8 .or. &
          DEF_ISO_SUBL_SKIN_MM < 0._r8 .or. DEF_ISO_ICE_SUPERSAT_SLOPE < 0._r8) &
         CALL CoLM_Stop ('Invalid isotope fractionation parameter.')

#ifndef BGC
      IF (DEF_USE_METHANE) CALL CoLM_Stop ('DEF_USE_METHANE needs a BGC build.')
#endif
      ! The fork checked this at compile time in define.h; methane is only
      ! built and tested on PFT/PC land cover, so keep the same requirement.
#if (!defined LULC_IGBP_PFT) && (!defined LULC_IGBP_PC)
      IF (DEF_USE_METHANE) CALL CoLM_Stop ('DEF_USE_METHANE needs LULC_IGBP_PFT or LULC_IGBP_PC.')
#endif
      IF (DEF_USE_METHANE .and. .not. mt_file_given (DEF_METHANE_PARAM_FILE)) &
         CALL CoLM_Stop ('DEF_USE_METHANE needs DEF_METHANE_PARAM_FILE.')

      IF (DEF_MT_BALANCE_ABORT_NBAD < 0 .or. DEF_MT_RESID_ABORT_NBAD < 0 .or. &
          DEF_MT_LULCC_ABORT_NBAD < 0) CALL CoLM_Stop ('Invalid DEF_MT_*_ABORT_NBAD.')

   END SUBROUTINE check_settings


   SUBROUTINE resolve_nh4_leaching ()
   ! Turn DEF_NH4_LEACHING into the logical read by BGC.  'auto'
   ! enables species-resolved NH4 leaching exactly when DEF_SOLUTE_NAMES
   ! contains a standalone NH4 token; 'on'/'off' force it.  The test is
   ! case-insensitive and tolerant of spaces around the commas.
   USE MOD_SPMD_Task, only: CoLM_Stop
   IMPLICIT NONE

   character(len=256) :: names
   character(len=16)  :: mode
   integer :: i, k, pos, n
   logical :: left_ok, right_ok

      mode = adjustl(DEF_NH4_LEACHING)
      CALL mt_lower_case (mode)

      SELECT CASE (trim(mode))
      CASE ('auto')
         names = trim(adjustl(DEF_SOLUTE_NAMES))
         CALL mt_lower_case (names)
         n = len_trim(names)
         DEF_USE_NH4_LEACHING = .false.
         IF (n > 0) THEN
            pos = 1
            DO
               k = index(names(pos:n), 'nh4')
               IF (k == 0) EXIT
               k = pos + k - 1
               left_ok  = (k == 1) .or. (names(k-1:k-1) == ' ') .or. &
                          (names(k-1:k-1) == ',')
               right_ok = (k+3 > n) .or. (names(k+3:k+3) == ' ') .or. &
                          (names(k+3:k+3) == ',')
               IF (left_ok .and. right_ok) THEN
                  DEF_USE_NH4_LEACHING = .true.
                  EXIT
               ENDIF
               pos = k + 3
               IF (pos > n) EXIT
            ENDDO
         ENDIF
      CASE ('on', 'true', 'yes')
         DEF_USE_NH4_LEACHING = .true.
      CASE ('off', 'false', 'no')
         DEF_USE_NH4_LEACHING = .false.
      CASE DEFAULT
         CALL CoLM_Stop ('DEF_NH4_LEACHING must be auto, on or off.')
      END SELECT

   END SUBROUTINE resolve_nh4_leaching


   SUBROUTINE mt_lower_case (s)
   ! ASCII lower-case in place; no locale/ICU dependency.
   IMPLICIT NONE
   character(len=*), intent(inout) :: s
   integer :: i, ic

      DO i = 1, len_trim(s)
         ic = iachar(s(i:i))
         IF (ic >= iachar('A') .and. ic <= iachar('Z')) &
            s(i:i) = achar(ic + 32)
      ENDDO

   END SUBROUTINE mt_lower_case


   logical FUNCTION mt_any_enabled ()
   ! True when at least one material is switched on.
   USE MOD_Namelist, only: DEF_USE_SEDIMENT
   IMPLICIT NONE
      mt_any_enabled = DEF_USE_ISOTOPE .or. len_trim(DEF_SOLUTE_NAMES) > 0 .or. DEF_USE_SEDIMENT &
         .or. DEF_USE_METHANE
   END FUNCTION mt_any_enabled


   logical FUNCTION mt_file_given (fname)
   IMPLICIT NONE
   character(len=*), intent(in) :: fname
      mt_file_given = len_trim(fname) > 0
      IF (mt_file_given) mt_file_given = trim(adjustl(fname)) /= 'null' &
                                   .and. trim(adjustl(fname)) /= 'NULL'
   END FUNCTION mt_file_given


   logical FUNCTION mt_namelist_group_present (fname, group)
   !-----------------------------------------------------------------------
   !  Whether a file holds "&group" at the start of a line, so that an
   !  optional group can be skipped rather than read to end of file.
   !-----------------------------------------------------------------------
   IMPLICIT NONE
   character(len=*), intent(in) :: fname, group

   character(len=1024) :: line
   character(len=128)  :: key
   logical :: fexists
   integer :: ierr, unit_nml, i, ia

      mt_namelist_group_present = .false.
      INQUIRE (file=trim(fname), exist=fexists)
      IF (.not. fexists) THEN
         IF (trim(group) == 'nl_colm_material') RETURN
         CALL stop_missing (fname)
      ENDIF

      key = '&'//trim(group)
      open (newunit=unit_nml, status='OLD', file=trim(fname), form='FORMATTED')
      DO
         read (unit_nml, '(A)', iostat=ierr) line
         IF (ierr /= 0) EXIT
         line = adjustl(line)
         DO i = 1, len_trim(line)
            ia = iachar(line(i:i))
            IF (ia >= iachar('A') .and. ia <= iachar('Z')) line(i:i) = achar(ia+32)
         ENDDO
         IF (line(1:1) == '$') line(1:1) = '&'
         IF (index(line, trim(key)) == 1) THEN
            i = len_trim(key) + 1
            IF (i > len_trim(line) .or. line(i:i) == ' ' .or. line(i:i) == '!') THEN
               mt_namelist_group_present = .true.
               EXIT
            ENDIF
         ENDIF
      ENDDO
      close (unit_nml)

   CONTAINS

      SUBROUTINE stop_missing (f)
      USE MOD_SPMD_Task, only: CoLM_Stop
      character(len=*), intent(in) :: f
         CALL CoLM_Stop ('Material parameter file not found: '//trim(f))
      END SUBROUTINE stop_missing

   END FUNCTION mt_namelist_group_present

END MODULE MOD_MT_Namelist
#endif
