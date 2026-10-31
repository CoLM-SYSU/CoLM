#include <define.h>
#if (defined MaterialTransport) && (defined BGC)
MODULE MOD_MT_Methane_Const
!=======================================================================
! methane constants
!=======================================================================
!
! PRIMARY REFERENCES (full bibliographic detail for citations used in
! parameter comments and biome lookup function headers below):
!
!   [Bridgham 2013] Bridgham, S. D., Cadillo-Quiroz, H., Keller, J. K., &
!     Zhuang, Q. (2013). Methane emissions from wetlands: biogeochemical,
!     microbial, and modeling perspectives from local to global scales.
!     Global Change Biology, 19(5), 1325-1346.  doi:10.1111/gcb.12131
!
!   [Hamilton 1996] Hamilton, S. K., Sippel, S. J., & Melack, J. M. (1996).
!     Inundation patterns in the Pantanal wetland of South America
!     determined from passive microwave remote sensing.
!     Archiv fur Hydrobiologie, 137(1), 1-23.
!
!   [Holzapfel-Pschorn 1985] Holzapfel-Pschorn, A., Conrad, R., & Seiler, W.
!     (1985). Production, oxidation and emission of methane in rice paddies.
!     FEMS Microbiology Ecology, 1(6), 343-351.
!     doi:10.1111/j.1574-6968.1985.tb01605.x
!
!   [Le Mer & Roger 2001] Le Mer, J., & Roger, P. (2001).  Production,
!     oxidation, emission and consumption of methane by soils: A review.
!     European Journal of Soil Biology, 37(1), 25-50.
!     doi:10.1016/S1164-5563(01)01067-6
!
!   [Marani & Alvalá 2007] Marani, L., & Alvalá, P. C. (2007).  Methane
!     emissions from lakes and floodplains in Pantanal, Brazil.
!     Atmospheric Environment, 41(8), 1627-1633.
!     doi:10.1016/j.atmosenv.2006.10.046
!
!   [Pangala 2017] Pangala, S. R., Enrich-Prast, A., Basso, L. S., Peixoto,
!     R. B., Bastviken, D., Hornibrook, E. R. C., Gatti, L. V., Marotta, H.,
!     Calazans, L. S. B., Sakuragui, C. M., Bastos, W. R., Malm, O., Gloor,
!     E., Miller, J. B., & Gauci, V. (2017).  Large emissions from
!     floodplain trees close the Amazon methane budget.
!     Nature, 552(7684), 230-234.  doi:10.1038/nature24639
!
!   [Walter 2001] Walter, B. P., Heimann, M., & Matthews, E. (2001).
!     Modeling modern methane emissions from natural wetlands: 1. Model
!     description and results.
!     J. Geophys. Res. Atmos., 106(D24), 34189-34206.
!     doi:10.1029/2001JD900165
!
!   [Wania 2010] Wania, R., Ross, I., & Prentice, I. C. (2010).
!     Implementation and evaluation of a new methane model within a
!     dynamic global vegetation model: LPJ-WHyMe v1.3.1.
!     Geoscientific Model Development, 3(2), 565-584.
!     doi:10.5194/gmd-3-565-2010
!
!   [Whalen & Reeburgh 1990] Whalen, S. C., & Reeburgh, W. S. (1990).
!     Consumption of atmospheric methane by tundra soils.
!     Nature, 346(6280), 160-162.  doi:10.1038/346160a0
!
! IMPORTANT: Specific numeric values in DEF_METHANE biome lookup arrays
! (f_methane_*, redoxlag_*, z0_methane_prod, hybrid_soil_threshold) are
! AUTHOR-SELECTED midpoints or order-of-magnitude estimates within the
! ranges discussed in these references, NOT direct numerical quotations.
! Default 'hybrid' mode is tuned against Pantanal in-situ flux
! [Marani & Alvalá 2007].
!=======================================================================
	USE MOD_Precision
	USE MOD_SPMD_Task, only: CoLM_Stop, p_is_master
	USE MOD_MT_Defs, only: material_lower
	USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite

	IMPLICIT NONE

	PRIVATE
	PUBLIC :: DEF_METHANE, DEF_METHANE_hydrology
	PUBLIC :: read_methane_namelist, configure_methane_inundation_mode
	PUBLIC :: methane_atm_mixing_ratio, methane_history_enabled
	PUBLIC :: methane_history_accumulation_mode
	PUBLIC :: DEF_wetland_finundation_scheme, DEF_file_GIEMS

	! CoLM202X kept these two in the host &nl_colm namelist; upstream has
	! neither.  The scheme is derived from DEF_METHANE%inundation_mode
	! (configure_methane_inundation_mode); the GIEMS file is read from
	! &nl_colm_methane_parameter.
	integer            :: DEF_wetland_finundation_scheme = 1
	character(len=256) :: DEF_file_GIEMS = 'null'

	!------------------------------------------------------------------
	! Constants for the Methane reactive material
	!------------------------------------------------------------------
	! Note some of these constants are also used in CNNitrifDenitrifMod

	integer :: iloop  ! loop index

	integer, public, parameter :: ngases      =   3     ! CH4, O2, & CO2
	integer, public, parameter :: METHANE_COMP_SOIL = 1
	integer, public, parameter :: METHANE_COMP_RICE = 2
	integer, public, parameter :: N_METHANE_COMP    = 2

	!------------------------------------------------------------------

	real(r8), public, parameter :: catomw = 12.011_r8 ! molar mass of C atoms (g/mol)
	! Carbon content of soil organic matter, gC per kg OM. Converts between
	! cellorg [kg OM m-3] and organic carbon: BgcLink divides by it going
	! C -> OM, the microbial substrate pools multiply by it going OM -> C.
	! Named because a bare 580 repeated across files is one edit away from
	! the two directions disagreeing, and nothing would catch that.
	real(r8), public, parameter :: gc_per_kg_om = 580._r8
	real(r8), public, parameter :: methane_atomw = 16.04_r8 ! molar mass of CH4 atoms (g/mol)

	real(r8), public :: s_con(ngases,4)    ! Schmidt # calculation constants (spp, #)
	data (s_con(1,iloop),iloop=1,4) /1898.0_r8, -110.1_r8, 2.834_r8, -0.02791_r8/ ! CH4
	data (s_con(2,iloop),iloop=1,4) /1801.0_r8, -120.1_r8, 3.7818_r8, -0.047608_r8/ ! O2
	data (s_con(3,iloop),iloop=1,4) /1911.0_r8, -113.7_r8, 2.967_r8, -0.02943_r8/ ! CO2

	real(r8), public :: d_con_w(ngases,3)    ! water diffusivity constants (spp, #)  (*10^-9 m2/s)
	data (d_con_w(1,iloop),iloop=1,3) /0.9798_r8, 0.02986_r8, 0.0004381_r8/ ! CH4
	data (d_con_w(2,iloop),iloop=1,3) /1.172_r8, 0.03443_r8, 0.0005048_r8/ ! O2
	data (d_con_w(3,iloop),iloop=1,3) /0.939_r8, 0.02671_r8, 0.0004095_r8/ ! CO2

	real(r8), public :: d_con_g(ngases,2)    ! gas diffusivity constants (spp, #) (*10^-4 m2/s)
	data (d_con_g(1,iloop),iloop=1,2) /0.1875_r8, 0.0013_r8/ ! CH4
	data (d_con_g(2,iloop),iloop=1,2) /0.1759_r8, 0.00117_r8/ ! O2
	data (d_con_g(3,iloop),iloop=1,2) /0.1325_r8, 0.0009_r8/ ! CO2

	real(r8), public :: c_h(ngases)    ! constant (K) for Henry's law (4.12, Wania)
	data c_h(1:3) /1600._r8, 1500._r8, 2400._r8/ ! CH4, O2, CO2

	real(r8), public :: kh_theta(ngases)    ! Henry's constant (mol/L/atm) at standard temperature (298K)
   data kh_theta(1:3) /1.4e-3_r8, 1.3e-3_r8, 3.4e-2_r8/ ! CH4, O2, CO2

	real(r8), public :: kh_tbase = 298.15_r8 ! base temperature for calculation of Henry's constant (K)

!------------------------------------------------------------------

	real(r8), public, parameter :: rgasm = 8.31446261815324_r8 ! Universal gas constant [J mol-1 K-1]
                                 ![J/K/mol]=[N*m/K/mol]=[Pa*m3/K/mol]
	!!! rgas Different from CoLM
	!!! rgas in CoLM is gas constant for dry air [J/kg/K]
	!!! not Universal gas constant
	real(r8), public, parameter :: rgasLatm = 0.08206_r8 ! L*atm/mol/K
   ! rgasLatm      = rgasm         / 101325 / 0.001
   ! [L*atm/mol/K] = [Pa*m3/K/mol] / [Pa/atm] / [m3/L]

	real(r8), public, parameter :: secspday = 86400._r8 ! Seconds per day

   type Methane_type
      ! ---------------------------------------------------- CH4 species controls ----------------------------------------------------
      ! These controls are private to the CH4 reactive material and are read from
      ! standard_ch4_parameter.nml with the rest of DEF_METHANE.  Keep the main
      ! nl_colm namelist species-neutral; it should only map reactive material
      ! names to parameter files via DEF_METHANE_PARAM_FILE.
      character(len=32) :: inundation_mode = 'hybrid'
      logical :: enable_rice_paddy = .false.
      logical :: only_wetland = .false.
      logical :: use_spatial_ph = .false.

      ! ---------------------------------------------------- Variable parameter ----------------------------------------------------
      ! methane production constants
      real(r8) :: q10methane =2._r8             ! additional Q10 for methane production ABOVE the soil decomposition temperature relationship (doc:Q10 Baseline:2 Range:1.5~4) (params:2.0)
      real(r8) :: f_methane = 0.2_r8            ! ratio of CH4 production to total C mineralization (Baseline:0.2 Range:NA params:0.2)

      ! Biome-specific f_methane lookup (CH4 yield = CH4 / total anaerobic decomp).
      ! When use_biome_f_methane=.true., per-patch f_methane is selected by
      ! climate zone (mirrors get_wetland_veg_proxy classification) instead
      ! of the global default 0.20.  Default .false. preserves backwards compat.
      !
      ! Values below are AUTHOR-SELECTED midpoints within published literature
      ! ranges, NOT direct numerical quotes from the cited papers.  The cited
      ! references provide the underlying f_methane / yield ranges (typically
      ! 0.03-0.50 across biomes per Bridgham 2013), from which the specific
      ! values below were chosen to match Pantanal in-situ flux (Marani &
      ! Alvalá 2007).  Re-tuning against other in-situ sites is encouraged.
      logical  :: use_biome_f_methane            = .false.
      real(r8) :: f_methane_tropical_peat        = 0.12_r8  ! Inferred from Bridgham 2013 GCB tropical range; Saunders et al. papyrus productivity work
      real(r8) :: f_methane_tropical_floodplain  = 0.07_r8  ! Tuned to match Marani & Alvalá 2007 Pantanal in-situ (200-300 / 30-80 mg/m²/d); within Bridgham 2013 0.03-0.15 tropical range
      real(r8) :: f_methane_floodplain           = 0.10_r8  ! Routing-activated floodplain / seasonally flooded non-rice soil CH4 yield
      real(r8) :: f_methane_temperate_marsh      = 0.15_r8  ! Inferred from Bridgham 2013 GCB Typha/Carex range 0.10-0.25
      real(r8) :: f_methane_boreal_fen           = 0.12_r8  ! Inferred from Wania et al. 2010 GMD LPJ-WHyMe fen range 0.10-0.20; lower than CTSM global 0.20
      real(r8) :: f_methane_boreal_bog           = 0.08_r8  ! Inferred from Bridgham 2013 Sphagnum range 0.05-0.15 (fermenter-dominated)
      real(r8) :: f_methane_rice_paddy           = 0.10_r8  ! Rice paddy CH4/CO2 ratio; conservative value within CLM4Me sensitivity range (0.1-0.3)
      real(r8) :: f_methane_upland_soil          = 0.05_r8  ! Inferred from Le Mer & Roger 2001 review (upland anaerobic microsites)

      ! Biome-specific redoxlag lookup (days).  When enabled, methane_prod
      ! consumes biome_redoxlag_patch(ipatch) instead of global DEF_METHANE%redoxlag.
      ! Concept: faster microbial response in warm tropical wetlands, slower
      ! in cold boreal peat.  Default .false. preserves backwards compatibility.
      !
      ! Values below are AUTHOR-CHOSEN order-of-magnitude estimates based on
      ! the qualitative conclusions in the cited papers, NOT direct numerical
      ! quotes.  Boreal slow-response from Whalen & Reeburgh 1990 (cold-microbe
      ! consumption); tropical rapid-response from Pangala et al. 2017 (Amazon
      ! tree stem CH4); CTSM legacy ~30d from Walter & Heimann 2001.
      ! Less constrained than f_methane lookup — wider uncertainty (50-100%).
      logical  :: use_biome_redoxlag             = .false.
      real(r8) :: redoxlag_tropical_peat         = 12._r8   ! Tropical mid-range
      real(r8) :: redoxlag_tropical_floodplain   = 7._r8    ! Tropical rapid (Pangala et al. 2017 Nature 552:230 qualitative)
      real(r8) :: redoxlag_temperate_marsh       = 20._r8   ! Temperate intermediate
      real(r8) :: redoxlag_boreal_fen            = 30._r8   ! CTSM legacy default (Walter & Heimann 2001 JGR ~14-30d transition)
      real(r8) :: redoxlag_boreal_bog            = 45._r8   ! Boreal slowest (cold Sphagnum, Whalen & Reeburgh 1990 Nature)
      real(r8) :: redoxlag_rice_paddy            = 25._r8   ! Managed flood/drain redox response; faster than wetland 30d but avoids unrealistically immediate production
      real(r8) :: redoxlag_upland_soil           = 30._r8   ! Default (only when Plan A activates upland CH4)

      ! Explicit methanogenesis depth attenuation.  CTSM BGC SOC profile
      ! decays with z0_BGC ~ 0.5 m, but boreal peatland observations
      ! (see Walter & Heimann 2001 JGR 106:34189) indicate CH4 production
      ! concentrates in the top 30 cm.  When z0_methane_prod > 0, partition_z
      ! is multiplied by exp(-z/z0_methane_prod) and the column re-normalized
      ! so total CH4 production is preserved (only redistributes vertically:
      ! top layers up, deep layers down).  Default 0 = disabled, falls back
      ! to CTSM BGC profile alone.
      !
      ! Recommended 0.30 m: chosen by author within the 0.2-0.5 m range
      ! discussed in Walter & Heimann 2001 (specific value 0.30 not directly
      ! quoted from paper).  Tropical wetlands may require shorter z0 (~0.15 m,
      ! Pangala 2017) — biome-specific z0 not yet implemented.
      real(r8) :: z0_methane_prod = 0._r8

      ! methane oxidation constants
	      ! Oxidation vmax is per aqueous/active-water volume; methane_oxid
	      ! multiplies by vol_aqu to produce bulk-soil mol m-3 s-1 rates.
	      real(r8) :: vmax_methane_oxid = 1.25e-5_r8       ! Unit: mol m-3-aqueous s-1
	      real(r8) :: vmax_oxid_unsat = 1.25e-6_r8         ! Unit: mol m-3-aqueous s-1
      real(r8) :: k_m = 5.e-3_r8                ! Michaelis-Menten oxidation rate constant for CH4 concentration (params:5e-3 code:5.e-6_r8 * 1000._r8) (doc:KCH4 Baseline:5e-3 Range:5e-4~5e-2)
      real(r8) :: k_m_unsat = 5.e-4_r8           ! Michaelis-Menten oxidation rate constant for CH4 concentration (params:5e-4 code:5.e-6_r8 * 1000._r8 / 10._r8) (doc:KCH4 Baseline:5e-3 Range:5e-4~5e-2)
      real(r8) :: k_m_o2 =2.e-2_r8             ! Michaelis-Menten oxidation rate constant for O2 concentration (params:2e-2 code:20.e-6_r8 * 1000._r8) (doc:KO2 Baseline:2e-2 Range:2e-3~2e-1)
      real(r8) :: q10_methane_oxid = 1.9_r8         ! Q10 oxidation constant (? params:1.9)
      ! Lake-only oxidation controls.  Defaults preserve the legacy/global
      ! methane oxidation kinetics; use these instead of changing k_m_o2 or
      ! vmax_methane_oxid globally when diagnosing lake CH4.
      real(r8) :: lake_oxid_scale = 1.0_r8          ! multiplicative factor on lake CH4 oxidation only
      real(r8) :: lake_k_m_o2 = -1.0_r8             ! >0 overrides k_m_o2 for lake oxidation only
      real(r8) :: lake_vmax_methane_oxid = -1.0_r8  ! >=0 overrides vmax_methane_oxid for lake oxidation only
      real(r8) :: lake_oxic_sediment_depth = -1.0_r8 ! m; >0 limits lake sediment oxidation to this top depth, -1 disables

      ! Experimental microbial-pool state.  Runtime validation rejects
      ! enabling it until biomass growth/loss is coupled to donor/sink carbon.
      logical  :: use_microbial_pools = .false.
      logical  :: use_microbial_flux_override = .false.
      logical  :: use_microbial_dormancy = .false.
      real(r8) :: B_init_methanogen = 1.0_r8
      real(r8) :: B_init_methanotroph = 1.0_r8
      real(r8) :: B_min_methanogen = 1.0e-3_r8
      real(r8) :: B_min_methanotroph = 1.0e-3_r8
      ! Biomass caps as fractions of local layer organic carbon
      ! (cellorg [kg OM m-3] * 580 gC kgOM-1).  Set <=0 to disable only
      ! this stricter fractional limit; the shared physical C ceiling remains.
      real(r8) :: B_max_fraction_methanogen = -1.0_r8
      real(r8) :: B_max_fraction_methanotroph = -1.0_r8
      real(r8) :: mu_max_methanogen = 0.2_r8
      real(r8) :: mu_max_methanotroph = 0.5_r8
      real(r8) :: gamma_methanogen = 0.05_r8
      real(r8) :: gamma_methanotroph = 0.10_r8
      real(r8) :: gamma_microbial_dormant = 0.005_r8
      real(r8) :: gamma_microbial_freeze = 0.001_r8
      real(r8) :: K_substrate_methanogen_pool = 0.04_r8   ! substrate-pool half-saturation [mol C m-3]
      real(r8) :: K_inh_O2_methanogen = 1.0e-3_r8
      real(r8) :: kappa_m_methanogen = 1.0e-2_r8
      real(r8) :: kappa_m_methanotroph = 1.0e-2_r8
      ! When microbial flux override is enabled, do not allow microbial
      ! production to exceed this multiple of the legacy production tendency.
      real(r8) :: max_microbe_prod_multiplier = 3.0_r8
      real(r8) :: q10_microbe_growth = 3.0_r8
      real(r8) :: T_ref_microbe = 298.15_r8
      real(r8) :: dormancy_rate_active = 0.1_r8
      real(r8) :: dormancy_rate_revive = 1.0_r8
      real(r8) :: dormancy_threshold_methanogen_fS = 0.1_r8
      real(r8) :: dormancy_threshold_methanogen_fO2 = 0.1_r8
      real(r8) :: dormancy_threshold_methanotroph_fS = 0.1_r8
      real(r8) :: dormancy_threshold_methanotroph_fO2 = 0.1_r8

      ! methane ebullition constants
      real(r8) :: vgc_max  =0.15_r8            ! gas-volume fraction scale for ebullition threshold/target (unitless)

      ! methane aerenchyma constants
      real(r8) :: nongrassporosratio = 1._r8/3._r8   ! Ratio of root porosity in non-grass to grass, used for aerenchyma transport (params:0.33)
      real(r8) :: poros_tiller = 0.3_r8    ! Porosity for grass tiller
      real(r8) :: unsat_aere_ratio = 0.05_r8/0.3_r8    ! Ratio to multiply upland vegetation aerenchyma porosity by compared to inundated systems (params:0.1666666667 code:0.05_r8 / 0.3_r8)
      real(r8) :: porosmin = 0.05_r8            ! minimum aerenchyma porosity (unitless)(params:0.05 code:0.05_r8)
      real(r8) :: aere_radius = 2.9e-3_r8 ! Aerenchyma radius
      real(r8) :: rob = 3._r8                 ! ratio of root length to vertical depth ("root obliquity") (params:3. code:3._r8)
      real(r8) :: scale_factor_aere = 1._r8   ! scale factor on the aerenchyma area for sensitivity tests (1 params:1.) (doc:Fa Baseline:1 Range:0.5~1.5)

      ! methane transport constants
      real(r8) :: scale_factor_gasdiff = 1._r8   ! For sensitivity tests; convection would allow this to be > 1(? params:1.) (doc:fD0 Basline:1 Range:1,10 Unit:m2 s-1)
      real(r8) :: scale_factor_liqdiff = 1._r8   ! For sensitivity tests; convection would allow this to be > 1(? params:1.) (doc:fD0 Basline:1 Range:1,10 Unit:m2 s-1)
      real(r8) :: lake_liqdiff_scale = 1._r8      ! Lake-only multiplier on saturated/liquid CH4 diffusivity
      real(r8) :: lake_o2_liqdiff_scale = 1._r8   ! Lake-only multiplier on saturated/liquid O2 diffusivity
      ! Fallback atmospheric CH4/O2 boundary conductance [m/s], used only
      ! when current CoLM Monin-Obukhov state is absent or invalid.
      real(r8) :: grnd_methane_cond_default = 1.e-6_r8

      ! -------------------------------------------------- Invariant parameter ----------------------------------------------------
      ! methane production constants
      real(r8) :: mino2lim = 0.2_r8         ! minimum anaerobic decomposition rate as a fraction of potential aerobic rate (0.2+ params:0.2)
      real(r8) :: q10methane_base = 295._r8 ! temperature at which the effective f_methane actually equals the constant f_methane (295+ params:295)
      real(r8) :: q10lakebase = 298._r8       ! (K) base temperature for lake CH4 production (params:298. code:298._r8)
      real(r8) :: cnscalefactor=1._r8        ! scale factor on CN decomposition for assigning methane flux (?- params:1.)

      real(r8) :: redoxlag =30._r8           ! Number of days to lag in the calculation of finundated_lag (30+ params:30.)
      real(r8) :: lake_decomp_fact =9.e-11_r8    ! Base decomposition rate (1/s) at 25C (1 params:9e-11)
      real(r8) :: redoxlag_vertical=30._r8   ! time lag (days) to inhibit production for newly unsaturated layers (30+ params:0.)
      real(r8) :: pHmax = 9._r8          ! maximum pH for methane production(params:9. code:9.)
      real(r8) :: pHmin = 2.2_r8         ! minimum pH for methane production(params:2.2 code:2.2)
      real(r8) :: oxinhib = 400._r8          ! inhibition of methane production by oxygen (m^3/mol) (400+? params:400.)

      ! methane oxidation constants
      real(r8) :: smp_crit =-2.4e5_r8            ! Critical soil moisture potential (mm) (params:-2.4e5)

      ! methane ebullition constants
      real(r8) :: bubble_f =0.57_r8            ! CH4 content in gas bubbles (Kellner et al. 2006)

      ! methane aerenchyma constants
      real(r8) :: aereoxid =0._r8            ! fraction of methane flux entering aerenchyma rhizosphere that will be(? params:0.)
      real(r8) :: tiller_C = 0.22_r8 ! Per tiller 0.22 g C [g C/tiller]

      ! methane transport constants
      real(r8) :: satpow  =2._r8             ! exponent on watsat for saturated soil solute diffusion (2? params:2.)
      real(r8) :: capthick = 100._r8         ! min thickness before assuming h2osfc is impermeable (mm) (params:100.code:100._r8)

      ! additional constants
      real(r8) :: atm_methane  = 1.7e-6_r8         ! Atmospheric CH4 mixing ratio fallback (mol/mol)
      logical :: use_transient_atm_methane = .false. ! true: read time-varying atmospheric CH4 from atm_methane_file
      character(len=256) :: atm_methane_file = 'null' ! ASCII table: year value, or year month value
      character(len=16) :: atm_methane_file_units = 'auto' ! auto, mol/mol, ppmv, or ppbv
      real(r8) :: om_frac_sf = 1._r8          ! Scale factor for organic matter fraction (unitless)(? params:NA)

      ! -------------------------------------------------- Switch ----------------------------------------------------
      logical :: use_aereoxid_prog = .true. ! if false then aereoxid is read off of
      ! the parameter file and may be modifed by the user (default aereoxid on the
      ! file is 0.0).

      logical :: transpirationloss = .true. ! switch for activating CH4 loss from transpiration
                                    ! Transpiration loss assumes that the methane concentration in dissolved soil
                                    ! water remains constant through the plant and is released when the water evaporates
                                    ! from the stomata.
                                    ! Currently hard-wired to true; impact is < 1 Tg CH4/yr

      logical :: allowlakeprod = .false. ! Switch to allow production under lakes based on soil carbon dataset
                                 ! (Methane can be produced, and CO2 produced from methane oxidation,
                                 ! which will slowly reduce the available carbon stock, if ! replenishlakec, but no other biogeochem is done.)
                                 ! Note: switching this off turns off ALL lake methane biogeochem. However, 0 values
                                 ! will still be averaged into the concentration _sat history fields.

	      logical :: usephfact = .true.  ! Switch to use pH factor in methane production; spatial pH input is required when use_spatial_ph is true.

      ! Smooth WTD->finundated transition for scheme 6 (logistic S-curve).
      ! Replaces the original step function (zwt<=0.30 -> 1, else 0).
      ! Reference: Walter & Heimann 2000 (sigmoid f(zwt)), Sundh 2000
      ! (exp(-zwt/0.3)), Bridgham 2013 GCB (review).
      !   finundated = 1 / (1 + exp((zwt - wtd_inflection)/wtd_steepness))
      ! Default wtd_inflection=0.30m preserves the original threshold as the
      ! S-curve's symmetric center (finundated = 0.5 at zwt = 30 cm).
      real(r8) :: wtd_inflection = 0.30_r8   ! [m] zwt at which finundated = 0.5 (wetland tile, patchtype==2)
      real(r8) :: wtd_steepness  = 0.05_r8   ! [m] S-curve width; smaller = sharper

      ! Soil-tile sigmoid params (patchtype==0/1): deeper inflection so seasonally
      ! wet upland (Pantanal floodplain etc.) can also produce CH4 when soil zwt
      ! rises into the 0.5-1.5 m range during wet season.  Disabled if =0
      ! (default 0 keeps backwards compat: soil patches use wetland's sigmoid).
      real(r8) :: wtd_inflection_soil = 0._r8 ! [m] >0 enables soil-tile dynamic flood extension
      real(r8) :: wtd_steepness_soil  = 0.3_r8

      ! Hybrid mode: on soil tiles (patchtype != 2) take finundated from the
      ! routing-published f_inund_flood_patch instead of sigmoid(zwt).  This
      ! captures tropical floodplain overland flooding (Pantanal, Amazon) that
      ! sigmoid(zwt) cannot resolve because soil zwt is too deep there.
      ! Wetland tiles still use sigmoid(zwt) to keep dyn_wtd seasonality.
      ! Requires GridRiverLakeFlow active + DEF_USE_Dynamic_Wetland=.true.
      logical  :: use_routing_for_soil = .false.

      ! Hybrid soil-tile threshold gate: when use_routing_for_soil=.true.,
      ! soil tile only produces CH4 when its routing fldfrc exceeds this
      ! threshold.  Below threshold, soil tile is treated as fully dry
      ! (finundated=0).  Conceptually represents hydrological connectivity:
      ! river overbank does not flood soil until stage exceeds bankfull.
      ! Default 0 disables the gate (back compat).
      !
      ! Recommended 0.05-0.10 for tropical 2-deg grid (empirically tuned
      ! against Pantanal P/T ratio match; not derived from any single paper).
      ! Pantanal routing fldfrc ranges 0.05-0.20 wet→dry, so threshold 0.05
      ! gates out only the driest months.
      real(r8) :: hybrid_soil_threshold    = 0._r8

      logical :: replenishlakec = .false. ! Finite lake-sediment carbon stock by default.
                                    ! Enabling replenishment imposes an external carbon source
                                    ! and is therefore suitable only for explicit sensitivity tests.

      logical :: methane_offline = .true.    ! Only offline land CH4 is implemented in this repository.
                                 ! Setting false is rejected during validation until a host atmosphere
                                 ! flux publisher and NEM-to-NEE coupling are available.

      logical :: methane_rmcnlim = .false.   ! Remove the N and low moisture limitations on SOM HR when calculating
                                 ! methanogenesis.
                                 ! Note: this option has not been extensively tested.
                                 ! Currently hardwired off.

      logical :: anoxicmicrosites = .false. ! Use Arah & Stephen 1998 expression to allow production above the water table
                                    ! Currently hardwired off; expression is crude.

      logical :: methane_frzout = .false.    ! Retired compatibility flag; must remain false.
                                 ! Ice is always excluded from mobile CH4 storage.

      ! public :: methane_conrd ! Read and initialize CH4 constants

      logical :: use_nitrif_denitrif = .true.

      logical :: anoxia  = .true. ! true => anoxia is applied to heterotrophic respiration also considered in CH4 model
                                    ! default value reset in controlMod
                                    ! Whether to enable the anoxia for the seasonally induatded zones
      logical :: use_vertical_redoxlag = .true. ! Whether to enable the vertical redox lag effect

      ! SIF is enabled only when coupled BGC does not already apply anoxia limits.
      logical :: bgc_anoxia_limits_decomp = .false.  ! true: BGC o_scalar already limits decomposition
      logical :: use_ch4_sif              = .true.   ! true: apply CH4 seasonal inundation factor
      ! If BGC later applies real o_scalar limits, set (true, false).

      ! Global-run diagnostics / guard rails.
      ! write_ch4_history=false suppresses all CH4 history variables.
      ! ch4_history_vars accepts:
      !   'core'       : compact global-run set (default): surface flux,
      !                  production, oxidation, and total column CH4;
      !   'diagnostic' : previous broad core set for debugging/global audits;
      !   'all'        : legacy behavior, emit every CH4 diagnostic;
      !   'none'       : suppress all CH4 diagnostics;
      !   comma-separated NetCDF variable names for an exact custom set.
      logical :: write_ch4_history = .true.
      character(len=4096) :: ch4_history_vars = 'core'
      ! A positive value aborts the whole MPI job when any single-timestep
      ! CH4/O2 numerical correction exceeds this mol/m2 threshold.  A
      ! non-positive value keeps the existing diagnostic-only behavior.
      real(r8) :: numerical_correction_fatal_threshold = -1._r8
      ! Tolerance, as a FRACTION OF PORE VOLUME, on the disagreement between the
      ! host soil water (wliq/denh2o + wice/denice) and what the inundated
      ! fraction implies, when methane partitions the column into its saturated
      ! and unsaturated halves.
      !
      ! The two quantities are computed independently -- vtot from the host soil
      ! state, finundated from a water-table S-curve -- so they cannot be
      ! compared at machine epsilon. The guard this replaces used
      ! 1.e-10*max(pore_volume,1), i.e. a flat 1e-10 m absolute tolerance
      ! (pore_volume is ~0.04 m, so the max() always selected 1.0), which is
      ! round-off level and aborted global runs on their first timestep.
      !
      ! Within the tolerance the host water is left exactly as it is; the
      ! saturated-column allocation is capped instead, so the area-weighted
      ! total still reproduces the host state and no water is created. The
      ! largest residual seen is reported once. Beyond the tolerance the run
      ! still stops: a disagreement that large is invalid forcing, not
      ! numerical drift.
      real(r8) :: host_water_tolerance = 0.05_r8
      ! By default a missing/non-positive lake depth emits one warning and
      ! continues.  Set true for production global runs to fail fast when
      ! lake CH4 is enabled but landdata lacks valid lakedepth.
      logical :: lake_zero_depth_fatal = .false.
      ! Targeted restart-continuity diagnostics for lake CH4/O2 state.
      ! Default is off.  When enabled, lake patches print one compact line
      ! for selected days/times so restart-vs-continuous divergence can be
      ! diagnosed without changing physics.
      logical :: lake_restart_debug = .false.
      integer :: lake_restart_debug_year = 1996
      integer :: lake_restart_debug_start_doy = 1
      integer :: lake_restart_debug_end_doy = 3
      integer :: lake_restart_debug_sec = 1800

      ! Water-balance finundated override for wetland tiles.
      ! Replaces scheme-computed finundated on patchtype==2 with
      ! clamp(wetwat / wetwatmax, 0, 1).  Soil patches keep scheme finundated.
      logical :: enable_wetwat_finundated_override = .false.

      ! Force the wetland unsaturated branch to use a dry soil-column state.
      ! This keeps saturated and unsaturated branches distinct when wetland
      ! finundated is supplied by a seasonal area signal.
      logical :: wetland_dry_unsat_branch = .false.

      ! Rice physiology/aerenchyma can remain active briefly after harvest,
      ! but inundation and water depth always come from host hydrology/routing.
      real(r8) :: rice_drain_window_days        = 30._r8

      ! R2 short-term SOC fix (methane-only): paddy soils accumulate SOC
      ! ~2-3x faster than upland soils under long flooding (Pan 2010 GCB,
      ! Inubushi 2003), but CoLM reads soil C profile as a grid-mean from
      ! cnsteadystate.nc with no paddy-specific value.  This systematically
      ! under-estimates the methanogenic substrate available on rice patches.
      ! At the same time CoLM assumes 100% straw return (stems all to litter
      ! at harvest), which slightly over-estimates substrate input; net of
      ! the two is ~ -30% to -50% on rice CH4 production.
      !
      ! Reserved compatibility parameter.  The former methane-only production
      ! multiplier was not carbon conservative because it did not debit the
      ! BGC substrate pool.  Validation therefore requires exactly 1 until a
      ! coupled carbon-debit path is implemented.
      real(r8) :: rice_substrate_boost          = 1.0_r8
   END type Methane_type

   type Methane_hydrology_type
      ! vdcf and pc remain in the namelist only for backward compatibility.
      ! The current CLM fill-and-spill formulation uses slopebeta/slopemax;
      ! non-default values for the retired knobs are rejected below.
      real(r8) :: vdcf = 2._r8
      real(r8) :: slopebeta = -3._r8
      real(r8) :: slopemax = 0.4_r8
      real(r8) :: pc = 0.4_r8
   END type Methane_hydrology_type

   type (Methane_type) :: DEF_METHANE
   type (Methane_hydrology_type) :: DEF_METHANE_hydrology

   real(r8), save :: atm_ch4_file_molmol(0:3000,12)
   logical,  save :: atm_ch4_file_loaded = .false.

CONTAINS

	   SUBROUTINE read_methane_namelist (nlfile)

   USE MOD_Namelist
   IMPLICIT NONE

   character(len=*), intent(in) :: nlfile
   ! Local variables
   logical :: fexists
   integer :: ierr
   character(len=512) :: iomsg
   integer :: unit_nml

   namelist /nl_colm_methane_parameter/ DEF_METHANE,DEF_METHANE_hydrology,DEF_file_GIEMS

      ! A model instance may be finalized and initialized again in the same
      ! process.  Reconstruct both parameter objects so fields omitted by the
      ! next namelist cannot inherit overrides from the previous instance.
      DEF_METHANE = Methane_type()
      DEF_METHANE_hydrology = Methane_hydrology_type()
      DEF_file_GIEMS = 'null'
      ! The forcing filename and contents are namelist state, so a new read
      ! must also invalidate the saved file cache even when the path is
      ! unchanged.
      atm_ch4_file_loaded = .false.
      atm_ch4_file_molmol(:,:) = -1._r8

      ! Read on every rank so all workers use the same methane constants.
      ! The previous master-only read left non-master ranks at default values
      ! whenever users overrode DEF_METHANE in the parameter namelist.
      INQUIRE (file=trim(nlfile), exist=fexists)
      IF (.not. fexists) THEN
         CALL CoLM_Stop (' ***** ERROR: methane parameter file does not exist: '// trim(nlfile))
      ENDIF
      open(newunit=unit_nml, status='OLD', file=trim(nlfile), form="FORMATTED")
      iomsg = ''
      read(unit_nml, nml=nl_colm_methane_parameter, iostat=ierr, iomsg=iomsg)
      IF (ierr /= 0) THEN
         close(unit_nml)
         IF (p_is_master) write(*,'(A,A,A)') &
            'ERROR read_methane_namelist: invalid &nl_colm_methane_parameter in ', trim(nlfile), ': '
         IF (p_is_master) write(*,'(A)') trim(iomsg)
         CALL CoLM_Stop (' ***** ERROR: Problem reading namelist: '// trim(nlfile))
      ENDIF
      close(unit_nml)

      CALL validate_methane_namelist ()

	   END SUBROUTINE read_methane_namelist

	   SUBROUTINE configure_methane_inundation_mode ()
	      ! Resolve the user-facing four-option CH4 inundation mode into the
	      ! internal scheme integer plus paired methane switches.  The old
	      ! integer scheme remains only as the physics dispatch key; users
	      ! should set DEF_METHANE%inundation_mode to one of:
	      !   wetwat, satellite/giems, routing, dynamic_wtd, hybrid.
	      USE MOD_Namelist, only: DEF_USE_Dynamic_Wetland
	      IMPLICIT NONE

	      character(len=32) :: mode

	      mode = material_lower(adjustl(trim(DEF_METHANE%inundation_mode)))

	      ! Reset only hydrology-mode-owned switches before dispatch so repeated
	      ! calls cannot retain a previous mode.  Biome production, redox,
	      ! vertical profile and numeric flood thresholds remain independent
	      ! namelist choices.
	      DEF_METHANE%enable_wetwat_finundated_override = .false.
	      DEF_METHANE%wetland_dry_unsat_branch = .true.
	      DEF_METHANE%use_routing_for_soil = .false.

	      SELECT CASE (trim(mode))
	      CASE ('wetwat')
	         DEF_wetland_finundation_scheme = 1
	         DEF_METHANE%enable_wetwat_finundated_override = .true.
	         DEF_METHANE%wetland_dry_unsat_branch = .true.
	         IF (DEF_USE_Dynamic_Wetland) THEN
	            IF (p_is_master) write(6,*) &
	               '***** ERROR: wetwat methane inundation mode requires DEF_USE_Dynamic_Wetland = .false.'
	            CALL CoLM_Stop (' ***** ERROR: invalid methane inundation mode / dynamic wetland combination')
	         ENDIF

	      CASE ('satellite','giems')
	         DEF_wetland_finundation_scheme = 5
	         DEF_METHANE%enable_wetwat_finundated_override = .false.
	         DEF_METHANE%wetland_dry_unsat_branch = .true.
	         IF (DEF_USE_Dynamic_Wetland) THEN
	            IF (p_is_master) write(6,*) &
	               '***** ERROR: satellite methane inundation mode requires DEF_USE_Dynamic_Wetland = .false.'
	            CALL CoLM_Stop (' ***** ERROR: invalid methane inundation mode / dynamic wetland combination')
	         ENDIF

	      CASE ('routing')
#ifndef GridRiverLakeFlow
	         ! The flood fraction comes only from MOD_Grid_RiverLakeFlow; without
	         ! it the soil column silently sees zero flooding.
	         IF (p_is_master) write(6,*) &
	            '***** ERROR: routing methane inundation mode requires a GridRiverLakeFlow-enabled kernel.'
	         CALL CoLM_Stop (' ***** ERROR: methane routing mode is unavailable in this kernel')
#endif
	         DEF_wetland_finundation_scheme = 7
	         DEF_METHANE%enable_wetwat_finundated_override = .false.
	         DEF_METHANE%wetland_dry_unsat_branch = .true.
	         IF (DEF_USE_Dynamic_Wetland) THEN
	            IF (p_is_master) write(6,*) &
	               '***** ERROR: routing methane inundation mode requires DEF_USE_Dynamic_Wetland = .false.'
	            CALL CoLM_Stop (' ***** ERROR: invalid methane inundation mode / dynamic wetland combination')
	         ENDIF

	      CASE ('dynamic_wtd','dynamic-wtd')
	         DEF_wetland_finundation_scheme = 6
	         DEF_METHANE%enable_wetwat_finundated_override = .false.
	         ! Dynamic wetland hydrology supplies the WTD forcing; keep the
	         ! dry unsaturated branch active for wetland tiles.
	         DEF_METHANE%wetland_dry_unsat_branch = .true.
	         DEF_METHANE%use_routing_for_soil = .false.
	         IF (.not. DEF_USE_Dynamic_Wetland) THEN
	            IF (p_is_master) write(6,*) &
	               '***** ERROR: dynamic_wtd requires DEF_USE_Dynamic_Wetland = .true.'
	            CALL CoLM_Stop (' ***** ERROR: invalid methane inundation mode / dynamic wetland combination')
	         ENDIF

	      CASE ('hybrid','dh_all_thr05','dyn_routing_hybrid')
#ifndef GridRiverLakeFlow
	         ! hybrid feeds the routing flood fraction to the soil column
	         ! (use_routing_for_soil); same requirement as 'routing'.
	         IF (p_is_master) write(6,*) &
	            '***** ERROR: hybrid methane inundation mode requires a GridRiverLakeFlow-enabled kernel.'
	         CALL CoLM_Stop (' ***** ERROR: methane hybrid mode is unavailable in this kernel')
#endif
		         ! Site-calibrated hybrid mode; not a globally validated default.
	         ! Combines routing and dynamic-WTD hydrology.  Biome yield,
	         ! redox lag and vertical source attenuation are independent
	         ! namelist controls; the standard parameter file enables their
	         ! recommended values explicitly.
	         ! Numeric values below are author-selected midpoints, not direct quotes:
	         !   - dyn_routing_hybrid: wetland sigmoid(zwt) + soil routing fldfrc
	         !   - biome f_methane lookup (range from Bridgham 2013 GCB review)
	         !   - biome redoxlag lookup (Pangala 2017 / Whalen 1990 qualitative)
	         !   - hybrid soil threshold 0.05 (empirical, tuned to Pantanal P/T)
	         !   - depth attenuation z0=0.30m (within Walter & Heimann 2001 range)
	         ! Tuned to match Pantanal in-situ flux (Marani & Alvalá 2007).
	         ! Names dyn_routing_hybrid/hybrid kept as backwards-compatible aliases.
	         DEF_wetland_finundation_scheme              = 6
	         DEF_METHANE%enable_wetwat_finundated_override = .false.
	         DEF_METHANE%wetland_dry_unsat_branch        = .true.
	         DEF_METHANE%use_routing_for_soil            = .true.
	         IF (.not. DEF_USE_Dynamic_Wetland) THEN
	            IF (p_is_master) write(6,*) &
	               '***** ERROR: hybrid mode requires DEF_USE_Dynamic_Wetland = .true.'
	            CALL CoLM_Stop (' ***** ERROR: invalid methane inundation mode / dynamic wetland combination')
	         ENDIF

	      CASE DEFAULT
	         IF (p_is_master) write(6,*) &
	            '***** ERROR: unsupported DEF_METHANE%inundation_mode = ', trim(DEF_METHANE%inundation_mode), &
	            '; expected wetwat, satellite, routing, dynamic_wtd, or hybrid.'
	         CALL CoLM_Stop (' ***** ERROR: unsupported methane inundation mode')
	      END SELECT

	      IF (p_is_master) write(6,'(A,A,A,I0,A,L1,A,L1)') &
	         ' CH4 inundation mode: ', trim(mode), &
	         ' -> scheme=', DEF_wetland_finundation_scheme, &
	         ' wetwat_override=', DEF_METHANE%enable_wetwat_finundated_override, &
	         ' dry_unsat=', DEF_METHANE%wetland_dry_unsat_branch
	   END SUBROUTINE configure_methane_inundation_mode

	   SUBROUTINE validate_methane_namelist ()
      ! Range-check user-overridable methane parameters.  Catches negative
      ! production rates, non-positive Q10 / Michaelis-Menten constants, and
      ! pH window inversions that would otherwise propagate as silent NaNs or
      ! negative fluxes.
	      IMPLICIT NONE
      logical :: bad

      bad = .false.

      ! Ordered comparisons do not reject NaN: every comparison with NaN is
      ! false.  Check the complete real-valued namelist surface once before
      ! the semantic range checks below so no user override can silently
      ! inject a non-finite rate, scale, threshold, or temperature.
      IF (.not. all(ieee_is_finite([ &
         DEF_METHANE%q10methane, DEF_METHANE%f_methane, DEF_METHANE%f_methane_tropical_peat, &
         DEF_METHANE%f_methane_tropical_floodplain, DEF_METHANE%f_methane_floodplain, &
         DEF_METHANE%f_methane_temperate_marsh, DEF_METHANE%f_methane_boreal_fen, &
         DEF_METHANE%f_methane_boreal_bog, DEF_METHANE%f_methane_rice_paddy, &
         DEF_METHANE%f_methane_upland_soil, DEF_METHANE%redoxlag_tropical_peat, &
         DEF_METHANE%redoxlag_tropical_floodplain, DEF_METHANE%redoxlag_temperate_marsh, &
         DEF_METHANE%redoxlag_boreal_fen, DEF_METHANE%redoxlag_boreal_bog, &
         DEF_METHANE%redoxlag_rice_paddy, DEF_METHANE%redoxlag_upland_soil, &
         DEF_METHANE%z0_methane_prod, DEF_METHANE%vmax_methane_oxid, &
         DEF_METHANE%vmax_oxid_unsat, DEF_METHANE%k_m, DEF_METHANE%k_m_unsat, &
         DEF_METHANE%k_m_o2, DEF_METHANE%q10_methane_oxid, DEF_METHANE%lake_oxid_scale, &
         DEF_METHANE%lake_k_m_o2, DEF_METHANE%lake_vmax_methane_oxid, &
         DEF_METHANE%lake_oxic_sediment_depth, DEF_METHANE%B_init_methanogen, &
         DEF_METHANE%B_init_methanotroph, DEF_METHANE%B_min_methanogen, &
         DEF_METHANE%B_min_methanotroph, DEF_METHANE%B_max_fraction_methanogen, &
         DEF_METHANE%B_max_fraction_methanotroph, DEF_METHANE%mu_max_methanogen, &
         DEF_METHANE%mu_max_methanotroph, DEF_METHANE%gamma_methanogen, &
         DEF_METHANE%gamma_methanotroph, DEF_METHANE%gamma_microbial_dormant, &
         DEF_METHANE%gamma_microbial_freeze, DEF_METHANE%K_substrate_methanogen_pool, &
         DEF_METHANE%K_inh_O2_methanogen, DEF_METHANE%kappa_m_methanogen, &
         DEF_METHANE%kappa_m_methanotroph, DEF_METHANE%max_microbe_prod_multiplier, &
         DEF_METHANE%q10_microbe_growth, DEF_METHANE%T_ref_microbe, &
         DEF_METHANE%dormancy_rate_active, DEF_METHANE%dormancy_rate_revive, &
         DEF_METHANE%dormancy_threshold_methanogen_fS, &
         DEF_METHANE%dormancy_threshold_methanogen_fO2, &
         DEF_METHANE%dormancy_threshold_methanotroph_fS, &
         DEF_METHANE%dormancy_threshold_methanotroph_fO2, DEF_METHANE%vgc_max, &
         DEF_METHANE%nongrassporosratio, DEF_METHANE%poros_tiller, &
         DEF_METHANE%unsat_aere_ratio, DEF_METHANE%porosmin, DEF_METHANE%aere_radius, &
         DEF_METHANE%rob, DEF_METHANE%scale_factor_aere, DEF_METHANE%scale_factor_gasdiff, &
         DEF_METHANE%scale_factor_liqdiff, DEF_METHANE%lake_liqdiff_scale, &
         DEF_METHANE%lake_o2_liqdiff_scale, DEF_METHANE%grnd_methane_cond_default, &
         DEF_METHANE%mino2lim, DEF_METHANE%q10methane_base, DEF_METHANE%q10lakebase, &
         DEF_METHANE%cnscalefactor, DEF_METHANE%redoxlag, DEF_METHANE%lake_decomp_fact, &
         DEF_METHANE%redoxlag_vertical, DEF_METHANE%pHmax, DEF_METHANE%pHmin, &
         DEF_METHANE%oxinhib, DEF_METHANE%smp_crit, DEF_METHANE%bubble_f, &
         DEF_METHANE%aereoxid, DEF_METHANE%tiller_C, DEF_METHANE%satpow, &
         DEF_METHANE%capthick, DEF_METHANE%atm_methane, DEF_METHANE%om_frac_sf, &
         DEF_METHANE%wtd_inflection, DEF_METHANE%wtd_steepness, &
         DEF_METHANE%wtd_inflection_soil, DEF_METHANE%wtd_steepness_soil, &
         DEF_METHANE%hybrid_soil_threshold, DEF_METHANE%rice_drain_window_days, &
         DEF_METHANE%rice_substrate_boost, DEF_METHANE%numerical_correction_fatal_threshold, &
         DEF_METHANE%host_water_tolerance, &
         DEF_METHANE_hydrology%vdcf, &
         DEF_METHANE_hydrology%slopebeta, DEF_METHANE_hydrology%slopemax, &
         DEF_METHANE_hydrology%pc]))) THEN
         IF (p_is_master) write(6,*) '***** ERROR: methane namelist contains a NaN or infinite real parameter.'
         bad = .true.
      ENDIF

      IF (.not. DEF_METHANE%methane_offline) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: methane_offline=.false. requires an online atmosphere/NEE coupling that is not implemented.'
         bad = .true.
      ENDIF
      IF (DEF_METHANE%methane_frzout) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: methane_frzout is retired; ice is always excluded from mobile CH4 storage.'
         bad = .true.
      ENDIF

      IF (DEF_METHANE%f_methane < 0._r8 .or. DEF_METHANE%f_methane > 0.5_r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: f_methane out of [0,0.5]: ', DEF_METHANE%f_methane
         bad = .true.
      ENDIF
      IF (DEF_METHANE%f_methane_tropical_peat < 0._r8 .or. &
          DEF_METHANE%f_methane_tropical_peat > 0.5_r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: f_methane_tropical_peat out of [0,0.5]: ', &
            DEF_METHANE%f_methane_tropical_peat
         bad = .true.
      ENDIF
      IF (DEF_METHANE%f_methane_tropical_floodplain < 0._r8 .or. &
          DEF_METHANE%f_methane_tropical_floodplain > 0.5_r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: f_methane_tropical_floodplain out of [0,0.5]: ', &
            DEF_METHANE%f_methane_tropical_floodplain
         bad = .true.
      ENDIF
      IF (DEF_METHANE%f_methane_floodplain < 0._r8 .or. &
          DEF_METHANE%f_methane_floodplain > 0.5_r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: f_methane_floodplain out of [0,0.5]: ', &
            DEF_METHANE%f_methane_floodplain
         bad = .true.
      ENDIF
      IF (DEF_METHANE%f_methane_temperate_marsh < 0._r8 .or. &
          DEF_METHANE%f_methane_temperate_marsh > 0.5_r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: f_methane_temperate_marsh out of [0,0.5]: ', &
            DEF_METHANE%f_methane_temperate_marsh
         bad = .true.
      ENDIF
      IF (DEF_METHANE%f_methane_boreal_fen < 0._r8 .or. &
          DEF_METHANE%f_methane_boreal_fen > 0.5_r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: f_methane_boreal_fen out of [0,0.5]: ', &
            DEF_METHANE%f_methane_boreal_fen
         bad = .true.
      ENDIF
      IF (DEF_METHANE%f_methane_boreal_bog < 0._r8 .or. &
          DEF_METHANE%f_methane_boreal_bog > 0.5_r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: f_methane_boreal_bog out of [0,0.5]: ', &
            DEF_METHANE%f_methane_boreal_bog
         bad = .true.
      ENDIF
      IF (DEF_METHANE%f_methane_rice_paddy < 0._r8 .or. &
          DEF_METHANE%f_methane_rice_paddy > 0.5_r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: f_methane_rice_paddy out of [0,0.5]: ', &
            DEF_METHANE%f_methane_rice_paddy
         bad = .true.
      ENDIF
      IF (DEF_METHANE%f_methane_upland_soil < 0._r8 .or. &
          DEF_METHANE%f_methane_upland_soil > 0.5_r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: f_methane_upland_soil out of [0,0.5]: ', &
            DEF_METHANE%f_methane_upland_soil
         bad = .true.
      ENDIF
      IF (DEF_METHANE%q10methane <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: q10methane must be > 0: ', DEF_METHANE%q10methane
         bad = .true.
      ENDIF
      IF (DEF_METHANE%q10_methane_oxid <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: q10_methane_oxid must be > 0: ', DEF_METHANE%q10_methane_oxid
         bad = .true.
      ENDIF
      IF (DEF_METHANE%vmax_methane_oxid < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: vmax_methane_oxid must be >= 0: ', DEF_METHANE%vmax_methane_oxid
         bad = .true.
      ENDIF
      IF (DEF_METHANE%vmax_oxid_unsat < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: vmax_oxid_unsat must be >= 0: ', DEF_METHANE%vmax_oxid_unsat
         bad = .true.
      ENDIF
      IF (DEF_METHANE%k_m       <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: k_m must be > 0: ', DEF_METHANE%k_m
         bad = .true.
      ENDIF
      IF (DEF_METHANE%k_m_unsat <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: k_m_unsat must be > 0: ', DEF_METHANE%k_m_unsat
         bad = .true.
      ENDIF
      IF (DEF_METHANE%k_m_o2    <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: k_m_o2 must be > 0: ', DEF_METHANE%k_m_o2
         bad = .true.
      ENDIF
      IF (DEF_METHANE%lake_oxid_scale < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: lake_oxid_scale must be >= 0: ', &
            DEF_METHANE%lake_oxid_scale
         bad = .true.
      ENDIF
      IF (DEF_METHANE%lake_k_m_o2 /= -1._r8 .and. DEF_METHANE%lake_k_m_o2 <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: lake_k_m_o2 must be -1 or > 0: ', &
            DEF_METHANE%lake_k_m_o2
         bad = .true.
      ENDIF
      IF (DEF_METHANE%lake_vmax_methane_oxid /= -1._r8 .and. &
          DEF_METHANE%lake_vmax_methane_oxid < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: lake_vmax_methane_oxid must be -1 or >= 0: ', &
            DEF_METHANE%lake_vmax_methane_oxid
         bad = .true.
      ENDIF
      IF (DEF_METHANE%lake_oxic_sediment_depth /= -1._r8 .and. &
          DEF_METHANE%lake_oxic_sediment_depth <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: lake_oxic_sediment_depth must be -1 or > 0: ', &
            DEF_METHANE%lake_oxic_sediment_depth
         bad = .true.
      ENDIF
      IF (DEF_METHANE%aereoxid < 0._r8 .or. DEF_METHANE%aereoxid > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: aereoxid out of [0,1]: ', DEF_METHANE%aereoxid
         bad = .true.
      ENDIF
      IF (DEF_METHANE%bubble_f <= 0._r8 .or. DEF_METHANE%bubble_f > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: bubble_f out of (0,1]: ', DEF_METHANE%bubble_f
         bad = .true.
      ENDIF
      IF (DEF_METHANE%pHmin >= DEF_METHANE%pHmax) THEN
         IF (p_is_master) write(6,*) '***** ERROR: pHmin >= pHmax: ', DEF_METHANE%pHmin, DEF_METHANE%pHmax
         bad = .true.
      ENDIF
      IF (DEF_METHANE%mino2lim < 0._r8 .or. DEF_METHANE%mino2lim > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: mino2lim out of [0,1]: ', DEF_METHANE%mino2lim
         bad = .true.
      ENDIF
      IF (DEF_METHANE%atm_methane < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: atm_methane must be >= 0: ', DEF_METHANE%atm_methane
         bad = .true.
      ENDIF
      ! A fraction of pore volume: 0 restores the machine-epsilon behaviour that
      ! aborted on its first timestep, and >=1 would accept any disagreement at
      ! all, which defeats the guard.
      IF (DEF_METHANE%host_water_tolerance < 0._r8 .or. &
          DEF_METHANE%host_water_tolerance >= 1._r8) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: host_water_tolerance out of [0,1): ', DEF_METHANE%host_water_tolerance
         bad = .true.
      ENDIF
      IF (DEF_METHANE%wtd_steepness <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wtd_steepness must be > 0: ', DEF_METHANE%wtd_steepness
         bad = .true.
      ENDIF
      IF (.not. ieee_is_finite(DEF_METHANE%wtd_inflection_soil) .or. &
          DEF_METHANE%wtd_inflection_soil < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wtd_inflection_soil must be finite and >= 0 m: ', &
            DEF_METHANE%wtd_inflection_soil
         bad = .true.
      ENDIF
      IF (.not. ieee_is_finite(DEF_METHANE%wtd_steepness_soil) .or. &
          DEF_METHANE%wtd_steepness_soil <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wtd_steepness_soil must be finite and > 0 m: ', &
            DEF_METHANE%wtd_steepness_soil
         bad = .true.
      ENDIF
      IF (.not. ieee_is_finite(DEF_METHANE%hybrid_soil_threshold) .or. &
          DEF_METHANE%hybrid_soil_threshold < 0._r8 .or. &
          DEF_METHANE%hybrid_soil_threshold > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: hybrid_soil_threshold must be finite and in [0,1]: ', &
            DEF_METHANE%hybrid_soil_threshold
         bad = .true.
      ENDIF
      IF (.not. ieee_is_finite(DEF_METHANE%z0_methane_prod) .or. &
          DEF_METHANE%z0_methane_prod < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: z0_methane_prod must be finite and >= 0 m: ', &
            DEF_METHANE%z0_methane_prod
         bad = .true.
      ENDIF
      IF (DEF_METHANE%vgc_max <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: vgc_max must be > 0: ', DEF_METHANE%vgc_max
         bad = .true.
      ENDIF
      IF (DEF_METHANE%poros_tiller < 0._r8 .or. DEF_METHANE%poros_tiller > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: poros_tiller out of [0,1]: ', DEF_METHANE%poros_tiller
         bad = .true.
      ENDIF
      IF (DEF_METHANE%nongrassporosratio < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: nongrassporosratio must be >= 0: ', DEF_METHANE%nongrassporosratio
         bad = .true.
      ENDIF
      IF (DEF_METHANE%unsat_aere_ratio < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: unsat_aere_ratio must be >= 0: ', DEF_METHANE%unsat_aere_ratio
         bad = .true.
      ENDIF
      IF (DEF_METHANE%porosmin < 0._r8 .or. DEF_METHANE%porosmin > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: porosmin out of [0,1]: ', DEF_METHANE%porosmin
         bad = .true.
      ENDIF
      IF (DEF_METHANE%aere_radius <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: aere_radius must be > 0: ', DEF_METHANE%aere_radius
         bad = .true.
      ENDIF
      IF (DEF_METHANE%rob <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: rob must be > 0: ', DEF_METHANE%rob
         bad = .true.
      ENDIF
      IF (DEF_METHANE%scale_factor_aere < 0._r8 .or. &
          DEF_METHANE%scale_factor_gasdiff < 0._r8 .or. &
          DEF_METHANE%scale_factor_liqdiff < 0._r8 .or. &
          DEF_METHANE%lake_liqdiff_scale < 0._r8 .or. &
          DEF_METHANE%lake_o2_liqdiff_scale < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: methane transport/aerenchyma scale factors must be >= 0: ', &
            DEF_METHANE%scale_factor_aere, DEF_METHANE%scale_factor_gasdiff, &
            DEF_METHANE%scale_factor_liqdiff, DEF_METHANE%lake_liqdiff_scale, &
            DEF_METHANE%lake_o2_liqdiff_scale
         bad = .true.
      ENDIF
      IF (DEF_METHANE%grnd_methane_cond_default <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: grnd_methane_cond_default must be > 0: ', &
            DEF_METHANE%grnd_methane_cond_default
         bad = .true.
      ENDIF
      IF (DEF_METHANE%q10methane_base <= 0._r8 .or. DEF_METHANE%q10lakebase <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: q10 base temperatures must be > 0 K: ', &
            DEF_METHANE%q10methane_base, DEF_METHANE%q10lakebase
         bad = .true.
      ENDIF
      IF (DEF_METHANE%cnscalefactor < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: cnscalefactor must be >= 0: ', DEF_METHANE%cnscalefactor
         bad = .true.
      ENDIF
      IF (DEF_METHANE%redoxlag < 0._r8 .or. DEF_METHANE%redoxlag_vertical < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: redox lags must be >= 0 days: ', &
            DEF_METHANE%redoxlag, DEF_METHANE%redoxlag_vertical
         bad = .true.
      ENDIF
      IF (min(DEF_METHANE%redoxlag_tropical_peat, &
              DEF_METHANE%redoxlag_tropical_floodplain, &
              DEF_METHANE%redoxlag_temperate_marsh, DEF_METHANE%redoxlag_boreal_fen, &
              DEF_METHANE%redoxlag_boreal_bog, DEF_METHANE%redoxlag_rice_paddy, &
              DEF_METHANE%redoxlag_upland_soil) < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: biome redox lags must be >= 0 days.'
         bad = .true.
      ENDIF
      IF (DEF_METHANE%oxinhib < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: oxinhib must be >= 0 m3 mol-1: ', DEF_METHANE%oxinhib
         bad = .true.
      ENDIF
      IF (DEF_METHANE%B_max_fraction_methanogen > 1._r8 .or. &
          DEF_METHANE%B_max_fraction_methanotroph > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: enabled microbial biomass fractions must be <= 1: ', &
            DEF_METHANE%B_max_fraction_methanogen, DEF_METHANE%B_max_fraction_methanotroph
         bad = .true.
      ENDIF
      IF (DEF_METHANE%lake_decomp_fact < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: lake_decomp_fact must be >= 0: ', DEF_METHANE%lake_decomp_fact
         bad = .true.
      ENDIF
      IF (DEF_METHANE%smp_crit >= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: smp_crit must be negative [mm]: ', DEF_METHANE%smp_crit
         bad = .true.
      ENDIF
      IF (DEF_METHANE%satpow <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: satpow must be > 0: ', DEF_METHANE%satpow
         bad = .true.
      ENDIF
      IF (DEF_METHANE%capthick < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: capthick must be >= 0: ', DEF_METHANE%capthick
         bad = .true.
      ENDIF
      ! Rice physiology and reserved substrate-control range checks.
      IF (DEF_METHANE%rice_drain_window_days <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: rice_drain_window_days must be > 0: ', &
            DEF_METHANE%rice_drain_window_days
         bad = .true.
      ENDIF
      IF (abs(DEF_METHANE%rice_substrate_boost - 1._r8) > 10._r8 * epsilon(1._r8)) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: rice_substrate_boost must remain 1 until methane production debits BGC carbon: ', &
            DEF_METHANE%rice_substrate_boost
         bad = .true.
      ENDIF
      IF (len_trim(DEF_METHANE%ch4_history_vars) == 0) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: ch4_history_vars must be core/diagnostic/all/none or a comma-separated list'
         bad = .true.
      ENDIF
      IF (DEF_METHANE%tiller_C <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: tiller_C must be > 0: ', DEF_METHANE%tiller_C
         bad = .true.
      ENDIF
      IF (DEF_METHANE%om_frac_sf < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: om_frac_sf must be >= 0: ', DEF_METHANE%om_frac_sf
         bad = .true.
      ENDIF
      IF (DEF_METHANE%K_substrate_methanogen_pool <= 0._r8 .or. &
          DEF_METHANE%K_inh_O2_methanogen <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: microbial half-saturation/inhibition constants must be > 0: ', &
            DEF_METHANE%K_substrate_methanogen_pool, DEF_METHANE%K_inh_O2_methanogen
         bad = .true.
      ENDIF
      IF (DEF_METHANE%B_init_methanogen < 0._r8 .or. DEF_METHANE%B_init_methanotroph < 0._r8 .or. &
          DEF_METHANE%B_min_methanogen < 0._r8 .or. DEF_METHANE%B_min_methanotroph < 0._r8 .or. &
          DEF_METHANE%mu_max_methanogen < 0._r8 .or. DEF_METHANE%mu_max_methanotroph < 0._r8 .or. &
          DEF_METHANE%gamma_methanogen < 0._r8 .or. DEF_METHANE%gamma_methanotroph < 0._r8 .or. &
          DEF_METHANE%gamma_microbial_dormant < 0._r8 .or. DEF_METHANE%gamma_microbial_freeze < 0._r8 .or. &
          DEF_METHANE%kappa_m_methanogen < 0._r8 .or. DEF_METHANE%kappa_m_methanotroph < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: microbial biomass/rate/loss parameters must be >= 0'
         bad = .true.
      ENDIF
      IF (DEF_METHANE%use_microbial_flux_override .and. &
          .not. DEF_METHANE%use_microbial_pools) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: use_microbial_flux_override requires use_microbial_pools.'
         bad = .true.
      ENDIF
      IF (DEF_METHANE%use_microbial_dormancy .and. &
          .not. DEF_METHANE%use_microbial_pools) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: use_microbial_dormancy requires use_microbial_pools.'
         bad = .true.
      ENDIF
      IF (DEF_METHANE%use_microbial_pools) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: microbial pools are disabled until biomass growth/loss has donor/sink carbon coupling.'
         bad = .true.
         IF (.not. ieee_is_finite(DEF_METHANE%B_max_fraction_methanogen) .or. &
             .not. ieee_is_finite(DEF_METHANE%B_max_fraction_methanotroph) .or. &
             DEF_METHANE%B_max_fraction_methanogen <= 0._r8 .or. &
             DEF_METHANE%B_max_fraction_methanotroph <= 0._r8 .or. &
             DEF_METHANE%B_max_fraction_methanogen > 1._r8 .or. &
             DEF_METHANE%B_max_fraction_methanotroph > 1._r8) THEN
            IF (p_is_master) write(6,*) &
               '***** ERROR: microbial pools require finite B_max fractions in (0,1]: ', &
               DEF_METHANE%B_max_fraction_methanogen, DEF_METHANE%B_max_fraction_methanotroph
            bad = .true.
         ENDIF
      ENDIF
      IF (DEF_METHANE%max_microbe_prod_multiplier <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: max_microbe_prod_multiplier must be > 0: ', &
            DEF_METHANE%max_microbe_prod_multiplier
         bad = .true.
      ENDIF
      IF (DEF_METHANE%q10_microbe_growth <= 0._r8 .or. DEF_METHANE%T_ref_microbe <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: microbial Q10 and reference temperature must be > 0: ', &
            DEF_METHANE%q10_microbe_growth, DEF_METHANE%T_ref_microbe
         bad = .true.
      ENDIF
      IF (DEF_METHANE%dormancy_rate_active < 0._r8 .or. DEF_METHANE%dormancy_rate_revive < 0._r8 .or. &
          DEF_METHANE%dormancy_threshold_methanogen_fS < 0._r8 .or. &
          DEF_METHANE%dormancy_threshold_methanogen_fO2 < 0._r8 .or. &
          DEF_METHANE%dormancy_threshold_methanotroph_fS < 0._r8 .or. &
          DEF_METHANE%dormancy_threshold_methanotroph_fO2 < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: microbial dormancy rates/thresholds must be >= 0'
         bad = .true.
      ENDIF
      IF (DEF_METHANE_hydrology%slopemax <= 0._r8 .or. &
          DEF_METHANE_hydrology%slopebeta >= 0._r8) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: methane hydrology requires slopemax > 0 and slopebeta < 0: ', &
            DEF_METHANE_hydrology%slopemax, DEF_METHANE_hydrology%slopebeta
         bad = .true.
      ENDIF
      IF (abs(DEF_METHANE_hydrology%vdcf - 2._r8) > 64._r8 * epsilon(1._r8) .or. &
          abs(DEF_METHANE_hydrology%pc - 0.4_r8) > 64._r8 * epsilon(1._r8)) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: retired methane hydrology knobs vdcf/pc must retain defaults 2.0/0.4: ', &
            DEF_METHANE_hydrology%vdcf, DEF_METHANE_hydrology%pc
         bad = .true.
      ENDIF
      IF (DEF_METHANE%use_transient_atm_methane .and. &
          (len_trim(DEF_METHANE%atm_methane_file) == 0 .or. &
           trim(material_lower(adjustl(DEF_METHANE%atm_methane_file))) == 'null')) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: transient atmospheric CH4 requires atm_methane_file.'
         bad = .true.
      ENDIF

      IF (bad) CALL CoLM_Stop (' ***** ERROR: methane namelist validation failed')
   END SUBROUTINE validate_methane_namelist

   real(r8) FUNCTION methane_atm_mixing_ratio (year, month)
      integer, intent(in) :: year, month
      integer :: iy, im

      methane_atm_mixing_ratio = DEF_METHANE%atm_methane
      IF (.not. DEF_METHANE%use_transient_atm_methane) RETURN

      IF (year < 0 .or. year > 3000 .or. month < 1 .or. month > 12) THEN
         write(6,*) 'ERROR: atmospheric CH4 request is outside file-table bounds: ', year, month
         CALL CoLM_Stop ('invalid transient atmospheric CH4 date')
      ENDIF

      CALL load_methane_atm_file ()
      iy = year
      im = month
      IF (.not. ieee_is_finite(atm_ch4_file_molmol(iy,im)) .or. &
          atm_ch4_file_molmol(iy,im) <= 0._r8) THEN
         write(6,*) 'ERROR: transient atmospheric CH4 table has no valid value for ', iy, im
         CALL CoLM_Stop ('incomplete transient atmospheric CH4 table')
      ENDIF
      methane_atm_mixing_ratio = atm_ch4_file_molmol(iy,im)
   END FUNCTION methane_atm_mixing_ratio

   integer FUNCTION methane_history_accumulation_mode ()
      ! AccFlux has three coarse storage modes.  The mode marker distinguishes
      ! off/core/custom, while AccFlux also stores an exact selector fingerprint
      ! because distinct custom selectors can accumulate different families.
      character(len=4096) :: selector

      methane_history_accumulation_mode = 0
      IF (.not. DEF_METHANE%write_ch4_history) RETURN

      selector = material_lower(adjustl(trim(DEF_METHANE%ch4_history_vars)))
      SELECT CASE (trim(selector))
      CASE ('none', 'off', 'false', '.false.')
         methane_history_accumulation_mode = 0
      CASE ('core', 'default', 'minimal', 'fast')
         methane_history_accumulation_mode = 1
      CASE DEFAULT
         methane_history_accumulation_mode = 2
      END SELECT
   END FUNCTION methane_history_accumulation_mode

   logical FUNCTION methane_history_enabled (varname)
      ! Runtime per-variable gate for methane history output.
      ! See DEF_METHANE%ch4_history_vars in Methane_type.
      character(len=*), intent(in) :: varname
      character(len=4096), save :: cached_raw = ' '
      character(len=4096), save :: cached_list = ' '
      integer, save :: cached_mode = -1
      character(len=4096) :: list
      character(len=256)  :: v

      methane_history_enabled = .false.
      IF (.not. DEF_METHANE%write_ch4_history) RETURN

      list = material_lower(adjustl(trim(DEF_METHANE%ch4_history_vars)))
      v    = material_lower(adjustl(trim(varname)))

      ! History gates are queried dozens of times per history write.  Cache
      ! the normalized selector/list so custom comma lists do not get reparsed
      ! for every variable.  Rebuild automatically if the namelist value is
      ! changed between calls.
      IF (trim(list) /= trim(cached_raw)) THEN
         cached_raw  = list
         cached_list = compact_commas(trim(list))
         SELECT CASE (trim(cached_list))
         CASE ('all','*')
            cached_mode = 1
         CASE ('none','off','false','.false.')
            cached_mode = 0
         CASE ('core','default','minimal','fast')
            cached_mode = 2
         CASE ('diagnostic','extended','debug')
            cached_mode = 3
         CASE DEFAULT
            cached_mode = 4
         END SELECT
      ENDIF

      SELECT CASE (cached_mode)
      CASE (1)
         methane_history_enabled = .true.
      CASE (2)
         methane_history_enabled = methane_history_is_core(trim(v))
      CASE (3)
         methane_history_enabled = methane_history_is_diagnostic(trim(v))
      CASE (4)
         methane_history_enabled = index(','//trim(cached_list)//',', ','//trim(v)//',') > 0
      CASE DEFAULT
         methane_history_enabled = .false.
      END SELECT
   END FUNCTION methane_history_enabled

   logical FUNCTION methane_history_is_core (v)
      character(len=*), intent(in) :: v

      ! Keep the default history set intentionally small for global CH4 runs.
      ! Use ch4_history_vars='diagnostic' to recover the previous broad core.
      SELECT CASE (trim(v))
      CASE ( &
         'f_methane_surf_flux_tot', &
         'f_methane_surf_flux_tot_active', &
         'f_methane_surf_flux_global_total_with_lake', &
	     'f_methane_surf_flux_global_phys_with_lake', &
	     'f_methane_balance_residual_global_with_lake', &
	     'f_methane_ch4_clip_credit_global_with_lake', &
         'f_methane_surf_flux_tot_phys', &
         'f_methane_balance_residual', &
         'f_methane_ch4_clip_credit', &
         'f_o2_cap_loss', &
         'f_o2_cap_gain', &
         'f_methane_surf_flux_wetland', &
         'f_methane_surf_flux_soil', &
         'f_methane_surf_flux_lake', &
         'f_methane_surf_flux_rice', &
         'f_methane_prod_tot', &
         'f_methane_oxid_tot', &
         'f_totcol_methane', &
         'f_co2_hr_methane')
         methane_history_is_core = .true.
      CASE DEFAULT
         methane_history_is_core = .false.
      END SELECT
   END FUNCTION methane_history_is_core

   logical FUNCTION methane_history_is_diagnostic (v)
      character(len=*), intent(in) :: v

      ! Broad diagnostic set preserved from the former 'core' selector.
      SELECT CASE (trim(v))
      CASE ( &
         'f_net_methane', &
         'f_co2_hr_methane', &
         'f_methane_surf_flux_tot', &
         'f_methane_surf_flux_tot_active', &
         'f_methane_surf_flux_active_total_without_lake', &
         'f_methane_surf_flux_global_total_with_lake', &
	     'f_methane_surf_flux_global_phys_with_lake', &
	     'f_methane_balance_residual_global_with_lake', &
	     'f_methane_ch4_clip_credit_global_with_lake', &
         'f_methane_surf_flux_tot_phys', &
         'f_methane_surf_aere', &
         'f_methane_surf_aere_soil', &
         'f_methane_surf_aere_rice', &
         'f_methane_surf_ebul', &
         'f_methane_surf_ebul_soil', &
         'f_methane_surf_ebul_rice', &
         'f_methane_surf_diff', &
         'f_methane_surf_diff_soil', &
         'f_methane_surf_diff_rice', &
         'f_methane_surf_diff_phys', &
         'f_methane_balance_residual', &
         'f_methane_ch4_clip_credit', &
         'f_o2_cap_loss', &
         'f_o2_cap_gain', &
         'f_methane_prod_tot', &
         'f_methane_prod_tot_soil', &
         'f_methane_prod_tot_rice', &
         'f_methane_oxid_tot', &
         'f_methane_oxid_tot_soil', &
         'f_methane_oxid_tot_rice', &
         'f_co2_decomp_tot', &
         'f_co2_oxid_tot', &
         'f_co2_net_tot', &
         'f_totcol_methane', &
         'f_grnd_methane_cond', &
         'f_methane_surf_flux_tot_lake', &
         'f_methane_surf_ebul_lake', &
         'f_methane_surf_diff_lake', &
         'f_methane_prod_tot_lake', &
         'f_methane_oxid_tot_lake', &
         'f_co2_net_tot_lake', &
         'f_totcol_methane_lake', &
         'f_lake_water_ch4_stock', &
         'f_lake_water_o2_stock', &
         'f_lake_water_ch4_oxid', &
         'f_lake_sed_ch4_flux', &
         'f_lake_sed_o2_flux', &
         'f_lake_air_o2_flux', &
         'f_forc_pmethanem', &
         'f_layer_sat_lag', &
         'f_annavg_finrw', &
         'f_methane_dfsat_tot', &
         'f_f_h2osfc', &
         'f_methane_finundated', &
         'f_methane_soil_finundated', &
         'f_methane_soil_zwt', &
         'f_inund_flood_patch', &
         'f_inund_flood_depth_patch', &
         'f_wetland_frac_patch', &
         'f_methane_surf_flux_wetland', &
         'f_methane_surf_flux_soil', &
         'f_methane_surf_flux_lake', &
         'f_methane_surf_flux_lake_intensive', &
         'f_methane_surf_flux_rice', &
         'f_methane_surf_flux_rice_intensive')
         methane_history_is_diagnostic = .true.
      CASE DEFAULT
         methane_history_is_diagnostic = .false.
      END SELECT
   END FUNCTION methane_history_is_diagnostic

   character(len=4096) FUNCTION compact_commas (s)
      ! Lower-level namelist convenience: allow users to write
      ! "f_a, f_b" with spaces after commas.
      character(len=*), intent(in) :: s
      integer :: i, n

      compact_commas = ' '
      n = 0
      DO i = 1, len_trim(s)
         IF (s(i:i) == ' ' .or. s(i:i) == char(9)) CYCLE
         n = n + 1
         IF (n <= len(compact_commas)) compact_commas(n:n) = s(i:i)
      ENDDO
   END FUNCTION compact_commas

   SUBROUTINE load_methane_atm_file ()
      character(len=512) :: line
      integer :: iu, ios, ios3, ios2, line_number
      integer :: yr, mon
      real(r8) :: val, val_molmol

      IF (atm_ch4_file_loaded) RETURN
      atm_ch4_file_loaded = .true.
      atm_ch4_file_molmol(:,:) = -1._r8

      open(newunit=iu, file=trim(DEF_METHANE%atm_methane_file), &
         status='old', action='read', iostat=ios)
      IF (ios /= 0) THEN
         write(6,*) 'ERROR: cannot open transient CH4 atmospheric file: ', &
            trim(DEF_METHANE%atm_methane_file)
         CALL CoLM_Stop ('cannot open transient atmospheric CH4 file')
      ENDIF

      line_number = 0
      DO
         read(iu,'(A)',iostat=ios) line
         IF (ios /= 0) EXIT
         line_number = line_number + 1
         line = adjustl(line)
         IF (len_trim(line) == 0) CYCLE
         IF (line(1:1) == '#' .or. line(1:1) == '!') CYCLE

         read(line,*,iostat=ios3) yr, mon, val
         IF (ios3 == 0) THEN
            IF (yr < 0 .or. yr > 3000 .or. mon < 1 .or. mon > 12) THEN
               write(6,*) 'ERROR: invalid atmospheric CH4 date at table line ', line_number
               CALL CoLM_Stop ('invalid transient atmospheric CH4 table date')
            ENDIF
            val_molmol = methane_atm_to_molmol(val)
            IF (.not. ieee_is_finite(val_molmol) .or. val_molmol <= 0._r8) THEN
               write(6,*) 'ERROR: invalid atmospheric CH4 value at table line ', line_number
               CALL CoLM_Stop ('invalid transient atmospheric CH4 value')
            ENDIF
            atm_ch4_file_molmol(yr,mon) = val_molmol
         ELSE
            read(line,*,iostat=ios2) yr, val
            IF (ios2 == 0 .and. yr >= 0 .and. yr <= 3000) THEN
               val_molmol = methane_atm_to_molmol(val)
               IF (.not. ieee_is_finite(val_molmol) .or. val_molmol <= 0._r8) THEN
                  write(6,*) 'ERROR: invalid atmospheric CH4 value at table line ', line_number
                  CALL CoLM_Stop ('invalid transient atmospheric CH4 value')
               ENDIF
               atm_ch4_file_molmol(yr,:) = val_molmol
            ELSE
               write(6,*) 'ERROR: malformed atmospheric CH4 table line ', line_number, ': ', trim(line)
               CALL CoLM_Stop ('malformed transient atmospheric CH4 table')
            ENDIF
         ENDIF
      ENDDO
      IF (ios > 0) THEN
         write(6,*) 'ERROR: failed while reading atmospheric CH4 table after line ', line_number
         CALL CoLM_Stop ('transient atmospheric CH4 table read failure')
      ENDIF
      close(iu)
   END SUBROUTINE load_methane_atm_file

   real(r8) FUNCTION methane_atm_to_molmol (val)
      real(r8), intent(in) :: val
      character(len=16) :: units

      units = adjustl(trim(DEF_METHANE%atm_methane_file_units))
      SELECT CASE (units)
      CASE ('mol/mol','molmol','vmr','MOL/MOL','MOLMOL','VMR')
         methane_atm_to_molmol = val
      CASE ('ppmv','ppm','PPMV','PPM')
         methane_atm_to_molmol = val * 1.e-6_r8
      CASE ('ppbv','ppb','PPBV','PPB')
         methane_atm_to_molmol = val * 1.e-9_r8
      CASE DEFAULT
         ! auto: accept common CH4 conventions: 1700 ppbv, 1.7 ppmv, or 1.7e-6 mol/mol.
         IF (val > 100._r8) THEN
            methane_atm_to_molmol = val * 1.e-9_r8
         ELSEIF (val > 0.1_r8) THEN
            methane_atm_to_molmol = val * 1.e-6_r8
         ELSE
            methane_atm_to_molmol = val
         ENDIF
      END SELECT
   END FUNCTION methane_atm_to_molmol
END MODULE MOD_MT_Methane_Const
#endif
