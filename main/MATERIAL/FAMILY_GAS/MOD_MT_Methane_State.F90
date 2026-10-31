#include <define.h>

#if (defined MaterialTransport) && (defined BGC)
MODULE MOD_MT_Methane_State
!=======================================================================
! Methane reactive-material state and restart fields.
!=======================================================================

   USE MOD_Precision
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   USE MOD_Vars_Global, only: nl_soil, spval, dz_soi, WATERBODY
   USE MOD_MT_Methane_Const, only: N_METHANE_COMP, &
      METHANE_COMP_SOIL, METHANE_COMP_RICE, DEF_METHANE

   IMPLICIT NONE
   SAVE
   PRIVATE

   ! Runtime patchtype convention used throughout methane physics.  This is
   ! distinct from the land-cover class code WATERBODY (=16/17), which is
   ! used only when remapping patchclass during LULCC.
   integer, parameter :: PATCHTYPE_LAKE = 4

   ! Explicit exports: keep module data private by default and expose
   ! only the current cross-module methane interface/state fields.
   PUBLIC :: accumulate_methane_lake_substep_diagnostics
   PUBLIC :: allocate_methane_state
   PUBLIC :: annavg_agnpp
   PUBLIC :: annavg_bgnpp
   PUBLIC :: annavg_finrw
   PUBLIC :: annavg_somhr
   PUBLIC :: annsum_counter
   PUBLIC :: biome_f_methane_patch
   PUBLIC :: biome_redoxlag_patch
   PUBLIC :: c_atm
   PUBLIC :: co2_aere_depth
   PUBLIC :: co2_aere_depth_sat
   PUBLIC :: co2_aere_depth_unsat
   PUBLIC :: co2_aere_tot
   PUBLIC :: co2_decomp_depth
   PUBLIC :: co2_decomp_depth_lake
   PUBLIC :: co2_decomp_depth_sat
   PUBLIC :: co2_decomp_depth_unsat
   PUBLIC :: co2_decomp_tot
   PUBLIC :: co2_decomp_tot_lake
   PUBLIC :: co2_decomp_tot_sat
   PUBLIC :: co2_decomp_tot_unsat
   PUBLIC :: co2_net_tot
   PUBLIC :: co2_net_tot_lake
   PUBLIC :: co2_net_tot_sat
   PUBLIC :: co2_net_tot_unsat
   PUBLIC :: co2_oxid_depth
   PUBLIC :: co2_oxid_depth_lake
   PUBLIC :: co2_oxid_depth_sat
   PUBLIC :: co2_oxid_depth_unsat
   PUBLIC :: co2_oxid_tot
   PUBLIC :: co2_oxid_tot_lake
   PUBLIC :: co2_oxid_tot_sat
   PUBLIC :: co2_oxid_tot_unsat
   PUBLIC :: compute_f_h2osfc
   PUBLIC :: conc_methane
   PUBLIC :: conc_methane_lake
   PUBLIC :: conc_methane_sat
   PUBLIC :: conc_methane_unsat
   PUBLIC :: conc_o2
   PUBLIC :: conc_o2_lake
   PUBLIC :: conc_o2_sat
   PUBLIC :: conc_o2_unsat
   PUBLIC :: deallocate_methane_state
   PUBLIC :: f_h2osfc
   PUBLIC :: f_inund_flood_depth_patch
   PUBLIC :: f_inund_flood_patch
   PUBLIC :: f_inund_levee_patch
   PUBLIC :: finundated_lag
   PUBLIC :: forc_pmethanem
   PUBLIC :: fsat_bef
   PUBLIC :: grnd_methane_cond
   PUBLIC :: grnd_methane_cond_lake
   PUBLIC :: grnd_methane_cond_sat
   PUBLIC :: grnd_methane_cond_unsat
   PUBLIC :: handle_methane_dry_lake_substep
   PUBLIC :: reset_methane_inactive_lake_diagnostics
   PUBLIC :: init_methane_wetland_fraction_cache
   PUBLIC :: initialize_methane_lake_soilc_from_surface
   PUBLIC :: methane_lake_soilc_srf, load_methane_lake_soilc_srf
   PUBLIC :: lake_air_o2_flux
   PUBLIC :: lake_sed_ch4_flux
   PUBLIC :: lake_sed_o2_flux
   PUBLIC :: lake_soilc
   PUBLIC :: lake_water_ch4_oxid
   PUBLIC :: lake_water_ch4_stock
   PUBLIC :: lake_water_o2_stock
   PUBLIC :: lake_frozen_ch4_stock
   PUBLIC :: lake_frozen_o2_stock
   PUBLIC :: lake_liquid_fraction_prev
   PUBLIC :: layer_sat_lag
   PUBLIC :: methane_aere_depth
   PUBLIC :: methane_aere_depth_sat
   PUBLIC :: methane_aere_depth_unsat
   PUBLIC :: methane_balance_residual
   PUBLIC :: methane_ch4_clip_credit
   PUBLIC :: methane_dfsat_tot
   PUBLIC :: methane_ebul_depth
   PUBLIC :: methane_ebul_depth_lake
   PUBLIC :: methane_ebul_depth_sat
   PUBLIC :: methane_ebul_depth_unsat
   PUBLIC :: methane_ebul_tot
   PUBLIC :: methane_ebul_tot_lake
   PUBLIC :: methane_ebul_tot_sat
   PUBLIC :: methane_ebul_tot_unsat
   PUBLIC :: methane_finundated
   PUBLIC :: methane_oxid_depth
   PUBLIC :: methane_oxid_depth_lake
   PUBLIC :: methane_oxid_depth_sat
   PUBLIC :: methane_oxid_depth_unsat
   PUBLIC :: methane_oxid_tot
   PUBLIC :: methane_oxid_tot_lake
   PUBLIC :: methane_oxid_tot_sat
   PUBLIC :: methane_oxid_tot_unsat
   PUBLIC :: methane_prod_depth
   PUBLIC :: methane_prod_depth_lake
   PUBLIC :: methane_prod_depth_sat
   PUBLIC :: methane_prod_depth_unsat
   PUBLIC :: methane_prod_tot
   PUBLIC :: methane_prod_tot_lake
   PUBLIC :: methane_prod_tot_sat
   PUBLIC :: methane_prod_tot_unsat
   PUBLIC :: methane_soil_finundated
   PUBLIC :: methane_soil_zwt
   PUBLIC :: methane_stress
   PUBLIC :: methane_stress_sat
   PUBLIC :: methane_stress_unsat
   PUBLIC :: methane_surf_aere
   PUBLIC :: methane_surf_aere_sat
   PUBLIC :: methane_surf_aere_unsat
   PUBLIC :: methane_surf_diff
   PUBLIC :: methane_surf_diff_lake
   PUBLIC :: methane_surf_diff_phys
   PUBLIC :: methane_surf_diff_phys_sat
   PUBLIC :: methane_surf_diff_phys_unsat
   PUBLIC :: methane_surf_diff_sat
   PUBLIC :: methane_surf_diff_unsat
   PUBLIC :: methane_surf_ebul
   PUBLIC :: methane_surf_ebul_lake
   PUBLIC :: methane_surf_ebul_sat
   PUBLIC :: methane_surf_ebul_unsat
   PUBLIC :: methane_surf_flux_lake
   PUBLIC :: methane_surf_flux_rice
   PUBLIC :: methane_surf_flux_soil
	PUBLIC :: methane_surf_aere_soil, methane_surf_aere_rice
	PUBLIC :: methane_surf_ebul_soil, methane_surf_ebul_rice
	PUBLIC :: methane_surf_diff_soil, methane_surf_diff_rice
	PUBLIC :: methane_prod_tot_soil, methane_prod_tot_rice
	PUBLIC :: methane_oxid_tot_soil, methane_oxid_tot_rice
   PUBLIC :: methane_surf_flux_tot
   PUBLIC :: methane_surf_flux_tot_lake
   PUBLIC :: methane_surf_flux_tot_phys
   PUBLIC :: methane_surf_flux_tot_sat
   PUBLIC :: methane_surf_flux_tot_unsat
   PUBLIC :: methane_surf_flux_wetland
   PUBLIC :: methane_tran_depth
   PUBLIC :: methane_tran_depth_sat
   PUBLIC :: methane_tran_depth_unsat
   PUBLIC :: net_methane
   PUBLIC :: co2_hr_methane
   PUBLIC :: net_methane_sat
   PUBLIC :: net_methane_unsat
   PUBLIC :: o2_aere_depth
   PUBLIC :: o2_aere_depth_sat
   PUBLIC :: o2_aere_depth_unsat
   PUBLIC :: o2_cap_gain
   PUBLIC :: o2_cap_loss
   PUBLIC :: o2_decomp_depth
   PUBLIC :: o2_decomp_depth_sat
   PUBLIC :: o2_decomp_depth_unsat
   PUBLIC :: o2_oxid_depth
   PUBLIC :: o2_oxid_depth_sat
   PUBLIC :: o2_oxid_depth_unsat
   PUBLIC :: o2stress
   PUBLIC :: o2stress_sat
   PUBLIC :: o2stress_unsat
   PUBLIC :: read_methane_restart
   PUBLIC :: remap_methane_lulcc_state
   PUBLIC :: publish_methane_levee_flood_patch
   PUBLIC :: publish_methane_flood_patch
   PUBLIC :: save_methane_lulcc_state
   PUBLIC :: tempavg_agnpp
   PUBLIC :: tempavg_bgnpp
   PUBLIC :: tempavg_finrw
   PUBLIC :: tempavg_somhr
   PUBLIC :: totcol_methane
   PUBLIC :: totcol_methane_lake
   PUBLIC :: totcol_methane_sat
   PUBLIC :: totcol_methane_unsat
   PUBLIC :: wetland_frac_per_patch
   PUBLIC :: write_methane_restart
   PUBLIC :: conc_o2_unsat_component, conc_o2_sat_component
   PUBLIC :: conc_methane_unsat_component, conc_methane_sat_component
   PUBLIC :: layer_sat_lag_component
   PUBLIC :: annavg_agnpp_component, annavg_bgnpp_component
   PUBLIC :: annavg_somhr_component, annavg_finrw_component
   PUBLIC :: tempavg_agnpp_component, tempavg_bgnpp_component
   PUBLIC :: annsum_counter_component, tempavg_somhr_component, tempavg_finrw_component
   PUBLIC :: fsat_bef_component, finundated_lag_component, rice_fraction_prev

   ! Public read-only data: external modules may inspect these arrays,
   ! but all writes stay inside this module through its APIs.
   PROTECTED :: f_inund_levee_patch, f_inund_flood_patch, f_inund_flood_depth_patch, &
      wetland_frac_per_patch, f_h2osfc

   ! -------------------- field declarations --------------------
   !!!! --------------------------------------------------------------------------------------------------------
   !!!!                                         sum data
   !!!! --------------------------------------------------------------------------------------------------------
	real(r8), allocatable :: net_methane           (:) ! CH4 oxidation - production, applied to BGC CO2-HR (mol/m2/s)
	real(r8), allocatable :: co2_hr_methane        (:) ! CO2-C part of the heterotrophic respiration, decomp_hr + catomw*net_methane (gC/m2/s)
	real(r8), allocatable :: methane_prod_depth      (:,:) ! production of CH4 in each soil layer  (mol/m3/s)
	real(r8), allocatable :: o2_decomp_depth     (:,:) ! O2 consumption during decomposition in each soil layer (mol/m3/s)
	real(r8), allocatable :: co2_decomp_depth    (:,:) ! diagnostic CO2 from decomposition/methanogenesis before O2-stress scaling (mol/m3/s)
	real(r8), allocatable :: methane_oxid_depth      (:,:) ! CH4 consumption rate via oxidation in each soil layer (mol/m3/s)
	real(r8), allocatable :: o2_oxid_depth       (:,:) ! O2 consumption rate via oxidation in each soil layer (mol/m3/s)
	real(r8), allocatable :: co2_oxid_depth      (:,:) ! CO2 production from CH4 oxidation (mol/m3/s)
	real(r8), allocatable :: methane_aere_depth      (:,:) ! CH4 loss rate via aerenchyma in each soil layer (mol/m3/s)
	real(r8), allocatable :: methane_tran_depth      (:,:) ! CH4 loss rate via transpiration in each soil layer (mol/m3/s)
	real(r8), allocatable :: o2_aere_depth       (:,:) ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s)
	real(r8), allocatable :: co2_aere_depth      (:,:) ! CO2 aerenchyma diagnostic flux (mol/m3/s)
	real(r8), allocatable :: methane_ebul_depth      (:,:) ! CH4 loss rate via ebullition in each soil layer (mol/m3/s)
	real(r8), allocatable :: o2stress            (:,:) ! Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs
	real(r8), allocatable :: methane_stress           (:,:) ! Ratio of methane available to the total per-timestep methane sinks
	real(r8), allocatable :: methane_surf_flux_tot     (:) ! CH4 flux to atm incl. numerical corrections (mol/m2/s)
	real(r8), allocatable :: methane_surf_flux_tot_phys(:) ! CH4 physical flux before CH4 clip/residual corrections (mol/m2/s)
	real(r8), allocatable :: methane_surf_aere         (:) ! Output: Total column CH4 aerenchyma (mol/m2/s)
	real(r8), allocatable :: methane_surf_ebul         (:) ! Output: CH4 ebullition to atmosphere (mol/m2/s)
	real(r8), allocatable :: methane_surf_diff         (:) ! Output: CH4 diffusion flux plus numerical closure (mol/m2/s)
	real(r8), allocatable :: methane_surf_diff_phys    (:) ! Output: CH4 pure-physics diffusion before clip/residual (mol/m2/s)
	real(r8), allocatable :: methane_balance_residual (:) ! numerical CH4 closure flux credited to methane_surf_diff (mol/m2/s)
	real(r8), allocatable :: methane_ch4_clip_credit(:) ! negative CH4 clip credited to methane_surf_diff (mol/m2/s)
	real(r8), allocatable :: o2_cap_loss              (:) ! O2 removed by post-solve physical cap (mol/m2/s)
	real(r8), allocatable :: o2_cap_gain              (:) ! O2 added by post-solve nonnegative floor (mol/m2/s)
	real(r8), allocatable :: methane_ebul_tot        (:) ! Output: Total column CH4 ebullition (mol/m2/s)
	real(r8), allocatable :: methane_prod_tot        (:) ! Output: Total column CH4 production (mol/m2/s)
	real(r8), allocatable :: methane_oxid_tot        (:) ! Output: Total column CH4 oxidation (mol/m2/s)
	real(r8), allocatable :: co2_decomp_tot       (:) ! total diagnostic CO2 from decomposition/methanogenesis before O2-stress scaling (mol/m2/s)
	real(r8), allocatable :: co2_oxid_tot         (:) ! total CO2 from CH4 oxidation (mol/m2/s)
	real(r8), allocatable :: co2_aere_tot         (:) ! total CO2 aerenchyma diagnostic (mol/m2/s)
	real(r8), allocatable :: co2_net_tot          (:) ! net diagnosed CO2 source from CH4 module (mol/m2/s)
   real(r8), allocatable :: totcol_methane             (:) ! total methane in soil column, start of timestep (mol/m2)
	real(r8), allocatable :: grnd_methane_cond         (:) ! effective material conductance incl. snow/pond resistance (m/s)
   real(r8), allocatable :: conc_o2             (:,:) ! O2 conc in each soil layer (mol/m3)
	real(r8), allocatable :: conc_methane            (:,:) ! CH4 conc in each soil layer (mol/m3)
   !!!! --------------------------------------------------------------------------------------------------------

   !!!! --------------------------------------------------------------------------------------------------------
   !!!!                                         sum data (unsaturated / saturated)
   !!!! --------------------------------------------------------------------------------------------------------
   real(r8), allocatable :: net_methane_unsat           (:)  ! CH4 oxidation - production, unsaturated subcolumn (mol/m2/s)
   real(r8), allocatable :: net_methane_sat             (:)  ! CH4 oxidation - production, saturated subcolumn (mol/m2/s)

   real(r8), allocatable :: methane_prod_depth_unsat      (:,:)  ! production of CH4 in each soil layer (unsaturated)  (mol/m3/s)
   real(r8), allocatable :: methane_prod_depth_sat        (:,:)  ! production of CH4 in each soil layer (saturated)    (mol/m3/s)

   real(r8), allocatable :: o2_decomp_depth_unsat     (:,:)  ! O2 consumption during decomposition (unsaturated)    (mol/m3/s)
   real(r8), allocatable :: o2_decomp_depth_sat       (:,:)  ! O2 consumption during decomposition (saturated)      (mol/m3/s)
   real(r8), allocatable :: co2_decomp_depth_unsat    (:,:)  ! diagnostic CO2 decomp/methanogenesis before O2-stress scaling (unsat) (mol/m3/s)
   real(r8), allocatable :: co2_decomp_depth_sat      (:,:)  ! diagnostic CO2 decomp/methanogenesis before O2-stress scaling (sat)   (mol/m3/s)

   real(r8), allocatable :: methane_oxid_depth_unsat      (:,:)  ! CH4 oxidation rate (unsaturated)                     (mol/m3/s)
   real(r8), allocatable :: methane_oxid_depth_sat        (:,:)  ! CH4 oxidation rate (saturated)                       (mol/m3/s)

   real(r8), allocatable :: o2_oxid_depth_unsat       (:,:)  ! O2 oxidation rate (unsaturated)                      (mol/m3/s)
   real(r8), allocatable :: o2_oxid_depth_sat         (:,:)  ! O2 oxidation rate (saturated)                        (mol/m3/s)
   real(r8), allocatable :: co2_oxid_depth_unsat      (:,:)  ! CO2 from CH4 oxidation (unsaturated)                 (mol/m3/s)
   real(r8), allocatable :: co2_oxid_depth_sat        (:,:)  ! CO2 from CH4 oxidation (saturated)                   (mol/m3/s)

   real(r8), allocatable :: methane_aere_depth_unsat      (:,:)  ! CH4 loss via aerenchyma (unsaturated)                (mol/m3/s)
   real(r8), allocatable :: methane_aere_depth_sat        (:,:)  ! CH4 loss via aerenchyma (saturated)                  (mol/m3/s)

   real(r8), allocatable :: methane_tran_depth_unsat      (:,:)  ! CH4 loss via transpiration (unsaturated)             (mol/m3/s)
   real(r8), allocatable :: methane_tran_depth_sat        (:,:)  ! CH4 loss via transpiration (saturated)               (mol/m3/s)

   real(r8), allocatable :: o2_aere_depth_unsat       (:,:)  ! O2 gain via aerenchyma (unsaturated)                 (mol/m3/s)
   real(r8), allocatable :: o2_aere_depth_sat         (:,:)  ! O2 gain via aerenchyma (saturated)                   (mol/m3/s)
   real(r8), allocatable :: co2_aere_depth_unsat      (:,:)  ! CO2 aerenchyma diagnostic (unsaturated)              (mol/m3/s)
   real(r8), allocatable :: co2_aere_depth_sat        (:,:)  ! CO2 aerenchyma diagnostic (saturated)                (mol/m3/s)

   real(r8), allocatable :: methane_ebul_depth_unsat      (:,:)  ! CH4 loss via ebullition (unsaturated)                (mol/m3/s)
   real(r8), allocatable :: methane_ebul_depth_sat        (:,:)  ! CH4 loss via ebullition (saturated)                  (mol/m3/s)

   real(r8), allocatable :: o2stress_unsat            (:,:)  ! O2 availability/demand ratio (unsaturated)           (-)
   real(r8), allocatable :: o2stress_sat              (:,:)  ! O2 availability/demand ratio (saturated)             (-)

   real(r8), allocatable :: methane_stress_unsat           (:,:)  ! CH4 availability/sinks ratio (unsaturated)           (-)
   real(r8), allocatable :: methane_stress_sat             (:,:)  ! CH4 availability/sinks ratio (saturated)             (-)

   real(r8), allocatable :: methane_surf_flux_tot_unsat     (:)  ! CH4 flux to atmosphere (unsaturated)                 (mol/m2/s)
   real(r8), allocatable :: methane_surf_flux_tot_sat       (:)  ! CH4 flux to atmosphere (saturated)                   (mol/m2/s)

   real(r8), allocatable :: methane_surf_aere_unsat         (:)  ! Total column CH4 aerenchyma (unsaturated)            (mol/m2/s)
   real(r8), allocatable :: methane_surf_aere_sat           (:)  ! Total column CH4 aerenchyma (saturated)              (mol/m2/s)

   real(r8), allocatable :: methane_surf_ebul_unsat         (:)  ! CH4 ebullition to atmosphere (unsaturated)           (mol/m2/s)
   real(r8), allocatable :: methane_surf_ebul_sat           (:)  ! CH4 ebullition to atmosphere (saturated)             (mol/m2/s)

   real(r8), allocatable :: methane_surf_diff_unsat         (:)  ! CH4 surface diffusive flux (unsaturated)             (mol/m2/s)
   real(r8), allocatable :: methane_surf_diff_sat           (:)  ! CH4 surface diffusive flux (saturated)               (mol/m2/s)
   real(r8), allocatable :: methane_surf_diff_phys_unsat    (:)  ! CH4 pure-physics diff (unsaturated)                  (mol/m2/s)
   real(r8), allocatable :: methane_surf_diff_phys_sat      (:)  ! CH4 pure-physics diff (saturated)                    (mol/m2/s)

   real(r8), allocatable :: methane_ebul_tot_unsat          (:)  ! Total column CH4 ebullition (unsaturated)            (mol/m2/s)
   real(r8), allocatable :: methane_ebul_tot_sat            (:)  ! Total column CH4 ebullition (saturated)              (mol/m2/s)

   real(r8), allocatable :: methane_prod_tot_unsat          (:)  ! Total column CH4 production (unsaturated)            (mol/m2/s)
   real(r8), allocatable :: methane_prod_tot_sat            (:)  ! Total column CH4 production (saturated)              (mol/m2/s)

   real(r8), allocatable :: methane_oxid_tot_unsat          (:)  ! Total column CH4 oxidation (unsaturated)             (mol/m2/s)
   real(r8), allocatable :: methane_oxid_tot_sat            (:)  ! Total column CH4 oxidation (saturated)               (mol/m2/s)
   real(r8), allocatable :: co2_decomp_tot_unsat       (:)  ! total CO2 decomp/methanogenesis (unsaturated)        (mol/m2/s)
   real(r8), allocatable :: co2_decomp_tot_sat         (:)  ! total CO2 decomp/methanogenesis (saturated)          (mol/m2/s)
   real(r8), allocatable :: co2_oxid_tot_unsat         (:)  ! total CO2 from CH4 oxidation (unsaturated)           (mol/m2/s)
   real(r8), allocatable :: co2_oxid_tot_sat           (:)  ! total CO2 from CH4 oxidation (saturated)             (mol/m2/s)
   real(r8), allocatable :: co2_net_tot_unsat          (:)  ! net diagnosed CO2 source (unsaturated)               (mol/m2/s)
   real(r8), allocatable :: co2_net_tot_sat            (:)  ! net diagnosed CO2 source (saturated)                 (mol/m2/s)

   real(r8), allocatable :: totcol_methane_unsat             (:)  ! total methane in soil column, start (unsaturated)     (mol/m2)
   real(r8), allocatable :: totcol_methane_sat               (:)  ! total methane in soil column, start (saturated)       (mol/m2)

   real(r8), allocatable :: grnd_methane_cond_unsat         (:)  ! effective material conductance incl. snow/pond (unsaturated) (m/s)
   real(r8), allocatable :: grnd_methane_cond_sat           (:)  ! effective material conductance incl. snow/pond (saturated) (m/s)

   real(r8), allocatable :: conc_o2_unsat             (:,:)  ! O2 concentration in each soil layer (unsaturated)    (mol/m3)
   real(r8), allocatable :: conc_o2_sat               (:,:)  ! O2 concentration in each soil layer (saturated)      (mol/m3)

	   real(r8), allocatable :: conc_methane_unsat            (:,:)  ! CH4 concentration in each soil layer (unsaturated)   (mol/m3)
	   real(r8), allocatable :: conc_methane_sat              (:,:)  ! CH4 concentration in each soil layer (saturated)     (mol/m3)
	   ! Independent soil/rice methane-process columns.  The legacy 2-D arrays
	   ! above remain patch-aggregate diagnostics and restart compatibility
	   ! fields; component state is the prognostic source for soil patches when
	   ! rice-paddy methane is enabled.
	   real(r8), allocatable :: conc_o2_unsat_component      (:,:,:)
	   real(r8), allocatable :: conc_o2_sat_component        (:,:,:)
	   real(r8), allocatable :: conc_methane_unsat_component (:,:,:)
	   real(r8), allocatable :: conc_methane_sat_component   (:,:,:)
	   !!!! --------------------------------------------------------------------------------------------------------

	   !!!! --------------------------------------------------------------------------------------------------------
	   !!!!                                         lake data (CTSM alignment)
	   !!!! --------------------------------------------------------------------------------------------------------
	   real(r8), allocatable :: methane_prod_depth_lake       (:,:)  ! lake CH4 production rate by layer (mol/m3/s)
	   real(r8), allocatable :: methane_oxid_depth_lake       (:,:)  ! lake CH4 oxidation rate by layer (mol/m3/s)
	   real(r8), allocatable :: methane_ebul_depth_lake       (:,:)  ! lake CH4 ebullition loss by layer (mol/m3/s)
	   real(r8), allocatable :: co2_decomp_depth_lake       (:,:)  ! lake CO2 from sediment decomposition/methanogenesis (mol/m3/s)
	   real(r8), allocatable :: co2_oxid_depth_lake         (:,:)  ! lake CO2 from CH4 oxidation (mol/m3/s)
	   real(r8), allocatable :: methane_surf_ebul_lake        (:)    ! lake surface ebullition flux (mol/m2/s)
	   real(r8), allocatable :: methane_surf_diff_lake        (:)    ! lake surface diffusive flux (mol/m2/s)
	   real(r8), allocatable :: methane_surf_flux_tot_lake    (:)    ! lake total surface CH4 flux = ebul + diff (mol/m2/s)
	   real(r8), allocatable :: methane_prod_tot_lake         (:)    ! lake total CH4 production (mol/m2/s)
	   real(r8), allocatable :: methane_oxid_tot_lake         (:)    ! lake total CH4 oxidation (mol/m2/s)
	   real(r8), allocatable :: methane_ebul_tot_lake         (:)    ! lake total ebullition (mol/m2/s)
	   real(r8), allocatable :: co2_decomp_tot_lake         (:)    ! lake total CO2 decomp/methanogenesis (mol/m2/s)
	   real(r8), allocatable :: co2_oxid_tot_lake           (:)    ! lake total CO2 oxidation product (mol/m2/s)
	   real(r8), allocatable :: co2_net_tot_lake            (:)    ! lake net diagnosed CO2 source (mol/m2/s)
	   real(r8), allocatable :: totcol_methane_lake           (:)    ! lake sediment CH4 column stock (mol/m2)
	   real(r8), allocatable :: grnd_methane_cond_lake        (:)    ! effective lake-atmosphere CH4 conductance (m/s)
	   real(r8), allocatable :: conc_o2_lake                  (:,:)  ! lake O2 concentration by layer (mol/m3)
	   real(r8), allocatable :: conc_methane_lake             (:,:)  ! lake CH4 concentration by layer (mol/m3)
	   real(r8), allocatable :: lake_water_ch4_stock          (:)    ! well-mixed lake-water CH4 inventory (mol/m2)
	   real(r8), allocatable :: lake_water_o2_stock           (:)    ! well-mixed lake-water O2 inventory (mol/m2)
	   real(r8), allocatable :: lake_frozen_ch4_stock         (:)    ! immobile CH4 retained during lake freeze (mol/m2)
	   real(r8), allocatable :: lake_frozen_o2_stock          (:)    ! immobile O2 retained during lake freeze (mol/m2)
	   real(r8), allocatable :: lake_liquid_fraction_prev     (:)    ! previous liquid fraction for conservative phase transfer
	   real(r8), allocatable :: lake_water_ch4_oxid           (:)    ! water-column CH4 oxidation (mol/m2/s)
	   real(r8), allocatable :: lake_sed_ch4_flux             (:)    ! sediment-to-water CH4 flux (mol/m2/s)
	   real(r8), allocatable :: lake_sed_o2_flux              (:)    ! sediment-to-water O2 flux (mol/m2/s)
	   real(r8), allocatable :: lake_air_o2_flux              (:)    ! water-to-atmosphere O2 flux (mol/m2/s)
	   !!!! --------------------------------------------------------------------------------------------------------

	   real(r8), allocatable :: c_atm               (:,:) ! CH4, O2, CO2 atmospheric conc  (mol/m3)
	real(r8), allocatable :: forc_pmethanem            (:) ! CH4 concentration in atmos. (pascals)
	real(r8), allocatable :: layer_sat_lag       (:,:)
	real(r8), allocatable :: layer_sat_lag_component (:,:,:)
	real(r8), allocatable :: lake_soilc          (:,:) ! total soil organic matter found in level (gC / m3)
   real(r8), allocatable :: annavg_agnpp          (:) ! annual average above-ground NPP (gC/m2/s)
	real(r8), allocatable :: annavg_bgnpp          (:) ! annual average below-ground NPP (gC/m2/s)
	real(r8), allocatable :: annavg_somhr          (:) ! annual average SOM heterotrophic resp. (gC/m2/s)
	real(r8), allocatable :: annavg_finrw          (:) ! respiration-weighted annual average of finundated
   real(r8), allocatable :: tempavg_agnpp         (:) ! temporary average above-ground NPP (gC/m2/s)
	real(r8), allocatable :: tempavg_bgnpp         (:) ! temporary average below-ground NPP (gC/m2/s)
	real(r8), allocatable :: annsum_counter        (:) ! seconds since last annual accumulator turnover
	real(r8), allocatable :: tempavg_somhr         (:) ! temporary average SOM heterotrophic resp. (gC/m2/s)
	real(r8), allocatable :: tempavg_finrw         (:) ! respiration-weighted annual average of finundated
	real(r8), allocatable :: annavg_agnpp_component   (:,:)
	real(r8), allocatable :: annavg_bgnpp_component   (:,:)
	real(r8), allocatable :: annavg_somhr_component   (:,:)
	real(r8), allocatable :: annavg_finrw_component   (:,:)
	real(r8), allocatable :: tempavg_agnpp_component  (:,:)
	real(r8), allocatable :: tempavg_bgnpp_component  (:,:)
	real(r8), allocatable :: annsum_counter_component (:,:)
	real(r8), allocatable :: tempavg_somhr_component  (:,:)
	real(r8), allocatable :: tempavg_finrw_component  (:,:)

   real(r8), allocatable :: fsat_bef              (:) ! finundated from previous timestep
   real(r8), allocatable :: finundated_lag        (:) ! time-lagged fractional inundated area
	real(r8), allocatable :: fsat_bef_component       (:,:)
	real(r8), allocatable :: finundated_lag_component (:,:)
	real(r8), allocatable :: rice_fraction_prev       (:)
   real(r8), allocatable :: methane_dfsat_tot         (:) ! CH4 flux to atm due to decreasing finundated [mol/m2/s]

   ! f_h2osfc: fractional area of surface water (0-1, dimensionless).
   ! Source-repo CLM5 microtopography-based prognostic h2osfc scheme.
   ! Maintained by compute_f_h2osfc (this module) before each methane_driver call.
	   real(r8), allocatable :: f_h2osfc              (:) ! fraction of surface water [-]
	   ! Diagnostics used to audit CH4 inundation choices and patchtype
	   ! contributions in history output.  These are not restart-critical state.
	   real(r8), allocatable :: methane_finundated        (:) ! actual CH4 finundated used by physics [-]
	   real(r8), allocatable :: methane_soil_finundated   (:) ! finundated on active soil/rice patches only [-]
	   real(r8), allocatable :: methane_soil_zwt          (:) ! zwt on active soil/rice patches only [m]
	   real(r8), allocatable :: methane_surf_flux_wetland (:) ! wetland contribution to CH4 surface flux [mol/m2/s]
	   real(r8), allocatable :: methane_surf_flux_soil    (:) ! non-rice soil contribution to CH4 surface flux [mol/m2/s]
	   real(r8), allocatable :: methane_surf_flux_lake    (:) ! lake contribution to CH4 surface flux [mol/m2/s]
	   real(r8), allocatable :: methane_surf_flux_rice    (:) ! rice-paddy contribution to CH4 surface flux [mol/m2/s]
	   real(r8), allocatable :: methane_surf_aere_soil(:), methane_surf_aere_rice(:)
	   real(r8), allocatable :: methane_surf_ebul_soil(:), methane_surf_ebul_rice(:)
	   real(r8), allocatable :: methane_surf_diff_soil(:), methane_surf_diff_rice(:)
	   real(r8), allocatable :: methane_prod_tot_soil(:), methane_prod_tot_rice(:)
	   real(r8), allocatable :: methane_oxid_tot_soil(:), methane_oxid_tot_rice(:)
	   ! Per-patch floodplain fraction (0-1) from GridRiverLakeFlow's levee
   ! diagnostic (levee_floodarea / topo_area), exposed to methane scheme 7.
   ! Default 0; populated by MOD_Grid_RiverLakeFlow via
   ! publish_levee_fldfrc_to_patches when GridRiverLakeFlow is compiled in and
   ! the catchment->patch mapping is available.
   real(r8), allocatable :: f_inund_levee_patch  (:) ! floodplain frac (-)

   ! General flood inundation fraction (levee+floodplain) for methane scheme 7.
   ! Published by MOD_Grid_RiverLakeFlow::publish_fldfrc_to_patches.
	   real(r8), allocatable :: f_inund_flood_patch  (:) ! flood frac (-)
	   real(r8), allocatable :: f_inund_flood_depth_patch (:) ! floodplain water depth [m]

	   ! Per-patch cache of the local wetland fraction in the parent
	   ! element/cell.  Used to convert grid/ucat-relative routing flood
	   ! fractions into wetland-relative finundated on patchtype==2:
	   !   fin_wet = min(1, flood_grid / max(wetland_frac_per_patch, 0.01)).
	   ! Built once after landpatch/patchtype are available, and rebuilt after
	   ! LULCC remaps.
	   real(r8), allocatable :: wetland_frac_per_patch (:) ! wetland area / active land area [-]

	   ! Biome-specific f_methane per-patch (set by Driver via get_biome_f_methane).
	   ! When DEF_METHANE%use_biome_f_methane=.true., methane_prod consumes this
	   ! array instead of the global DEF_METHANE%f_methane scalar.
	   real(r8), allocatable :: biome_f_methane_patch (:) ! [mol CH4 / mol CO2 anaerobic decomp]

	   ! Biome-specific redoxlag per-patch (set by Driver via get_biome_redoxlag).
	   ! When DEF_METHANE%use_biome_redoxlag=.true., methane uses this array
	   ! for finundated_lag time constant instead of global DEF_METHANE%redoxlag.
	   real(r8), allocatable :: biome_redoxlag_patch (:)  ! [days]

   ! LULCC remap snapshot. These arrays are private old-layout copies used
   ! to rebuild Methane state after landpatch/numpatch changes.
   logical :: methane_lulcc_snapshot_valid = .false.
   real(r8), allocatable :: lulcc_conc_o2_old(:,:), lulcc_conc_methane_old(:,:)
   real(r8), allocatable :: lulcc_co2_decomp_depth_old(:,:), lulcc_co2_oxid_depth_old(:,:), lulcc_co2_aere_depth_old(:,:)
   real(r8), allocatable :: lulcc_co2_decomp_depth_unsat_old(:,:), lulcc_co2_decomp_depth_sat_old(:,:)
   real(r8), allocatable :: lulcc_co2_oxid_depth_unsat_old(:,:), lulcc_co2_oxid_depth_sat_old(:,:)
   real(r8), allocatable :: lulcc_co2_aere_depth_unsat_old(:,:), lulcc_co2_aere_depth_sat_old(:,:)
   real(r8), allocatable :: lulcc_co2_decomp_depth_lake_old(:,:), lulcc_co2_oxid_depth_lake_old(:,:)
   real(r8), allocatable :: lulcc_conc_o2_unsat_old(:,:), lulcc_conc_o2_sat_old(:,:)
   real(r8), allocatable :: lulcc_conc_methane_unsat_old(:,:), lulcc_conc_methane_sat_old(:,:)
   real(r8), allocatable :: lulcc_conc_o2_lake_old(:,:), lulcc_conc_methane_lake_old(:,:)
   real(r8), allocatable :: lulcc_totcol_methane_old(:), lulcc_totcol_methane_unsat_old(:)
   real(r8), allocatable :: lulcc_co2_decomp_tot_old(:), lulcc_co2_oxid_tot_old(:), lulcc_co2_aere_tot_old(:), lulcc_co2_net_tot_old(:)
   real(r8), allocatable :: lulcc_co2_decomp_tot_unsat_old(:), lulcc_co2_decomp_tot_sat_old(:)
   real(r8), allocatable :: lulcc_co2_oxid_tot_unsat_old(:), lulcc_co2_oxid_tot_sat_old(:)
   real(r8), allocatable :: lulcc_co2_net_tot_unsat_old(:), lulcc_co2_net_tot_sat_old(:)
   real(r8), allocatable :: lulcc_co2_decomp_tot_lake_old(:), lulcc_co2_oxid_tot_lake_old(:), lulcc_co2_net_tot_lake_old(:)
   real(r8), allocatable :: lulcc_totcol_methane_sat_old(:), lulcc_totcol_methane_lake_old(:)
   real(r8), allocatable :: lulcc_lake_water_ch4_stock_old(:), lulcc_lake_water_o2_stock_old(:)
   real(r8), allocatable :: lulcc_lake_frozen_ch4_stock_old(:), lulcc_lake_frozen_o2_stock_old(:)
   real(r8), allocatable :: lulcc_lake_liquid_fraction_prev_old(:)
   real(r8), allocatable :: lulcc_grnd_methane_cond_old(:), lulcc_grnd_methane_cond_unsat_old(:)
   real(r8), allocatable :: lulcc_grnd_methane_cond_sat_old(:), lulcc_grnd_methane_cond_lake_old(:)
   real(r8), allocatable :: lulcc_layer_sat_lag_old(:,:), lulcc_lake_soilc_old(:,:)
   real(r8), allocatable :: lulcc_annavg_agnpp_old(:), lulcc_annavg_bgnpp_old(:)
   real(r8), allocatable :: lulcc_annavg_somhr_old(:), lulcc_annavg_finrw_old(:)
   real(r8), allocatable :: lulcc_tempavg_agnpp_old(:), lulcc_tempavg_bgnpp_old(:)
   real(r8), allocatable :: lulcc_annsum_counter_old(:), lulcc_tempavg_somhr_old(:)
   real(r8), allocatable :: lulcc_tempavg_finrw_old(:), lulcc_fsat_bef_old(:)
   real(r8), allocatable :: lulcc_finundated_lag_old(:), lulcc_methane_dfsat_tot_old(:)
   real(r8), allocatable :: lulcc_f_h2osfc_old(:), lulcc_forc_pmethanem_old(:)
	   real(r8), allocatable :: lulcc_f_inund_levee_patch_old(:)
	   real(r8), allocatable :: lulcc_f_inund_flood_patch_old(:)
	   real(r8), allocatable :: lulcc_f_inund_flood_depth_patch_old(:)
   real(r8), allocatable :: lulcc_c_atm_old(:,:)
	real(r8), allocatable :: lulcc_conc_o2_unsat_component_old(:,:,:)
	real(r8), allocatable :: lulcc_conc_o2_sat_component_old(:,:,:)
	real(r8), allocatable :: lulcc_conc_methane_unsat_component_old(:,:,:)
	real(r8), allocatable :: lulcc_conc_methane_sat_component_old(:,:,:)
	real(r8), allocatable :: lulcc_layer_sat_lag_component_old(:,:,:)
	real(r8), allocatable :: lulcc_annavg_agnpp_component_old(:,:)
	real(r8), allocatable :: lulcc_annavg_bgnpp_component_old(:,:)
	real(r8), allocatable :: lulcc_annavg_somhr_component_old(:,:)
	real(r8), allocatable :: lulcc_annavg_finrw_component_old(:,:)
	real(r8), allocatable :: lulcc_tempavg_agnpp_component_old(:,:)
	real(r8), allocatable :: lulcc_tempavg_bgnpp_component_old(:,:)
	real(r8), allocatable :: lulcc_annsum_counter_component_old(:,:)
	real(r8), allocatable :: lulcc_tempavg_somhr_component_old(:,:)
	real(r8), allocatable :: lulcc_tempavg_finrw_component_old(:,:)
	real(r8), allocatable :: lulcc_fsat_bef_component_old(:,:)
	real(r8), allocatable :: lulcc_finundated_lag_component_old(:,:)
	real(r8), allocatable :: lulcc_rice_fraction_prev_old(:)

	   ! Temporary lake-substep history buffers.  Lake methane runs on the
	   ! WATERBODY physics substep, while the normal history accumulator is
	   ! called once per land timestep.  These buffers keep time-weighted
	   ! diagnostic rates over the substeps and write the per-timestep mean
	   ! back to the module fields before accumulate_methane_fluxes samples
	   ! them.  Prognostic states such as concentrations and lake_soilc are
	   ! intentionally left at their final substep values.
	   integer, parameter :: methane_lake_substep_n2d = 44
	   integer, parameter :: methane_lake_substep_n1d = 57
	   real(r8), allocatable :: methane_lake_substep_acc2d(:,:)
	   real(r8), allocatable :: methane_lake_substep_acc1d(:)
	   integer :: methane_lake_substep_cached_ipatch = -1
	   integer :: methane_lake_substep_next_isub = 1

   ! Lake-bottom soil carbon of the surface data (gC/m3), the cold start of
   ! lake_soilc.  CoLM202X kept it as the time invariant lake_soilc_srf
   ! filled by mkinidata; upstream has no such field, so methane reads it
   ! itself (load_methane_lake_soilc_srf).
   real(r8), allocatable :: methane_lake_soilc_srf(:,:)

   ! -------------------- API --------------------

CONTAINS

   ! -------
   ! Lake-bottom soil carbon from <landdata>/soil/<year>/lake_soilc_patches.nc
   ! (mksrfdata Aggregation_LakeSoilC), as CoLM202X mkinidata read it; a
   ! missing file gives zero.  SinglePoint has no such file and uses the
   ! organic matter of the lake patch (580 gC per kg OM).  Collective: all
   ! ranks call it.
   SUBROUTINE load_methane_lake_soilc_srf (lc_year, dir_landdata, numpatch)

      USE MOD_SPMD_Task, only: p_is_worker
      USE MOD_Vars_TimeInvariants, only: patchtype, OM_density
#ifndef SinglePoint
      USE MOD_LandPatch, only: landpatch
      USE MOD_NetCDFVector, only: ncio_read_vector
#endif
      integer, intent(in) :: lc_year
      character(len=*), intent(in) :: dir_landdata
      integer, intent(in) :: numpatch

      real(r8), parameter :: carbon_per_kg_om = 580._r8
      character(len=4)   :: cyear
      character(len=512) :: lndname
      integer :: ipatch

      IF (allocated(methane_lake_soilc_srf)) deallocate (methane_lake_soilc_srf)
      IF (p_is_worker) THEN
         allocate (methane_lake_soilc_srf (nl_soil,numpatch))
         methane_lake_soilc_srf(:,:) = 0._r8
      ENDIF
      IF (.not. DEF_METHANE%allowlakeprod) RETURN

#ifndef SinglePoint
      write(cyear,'(i4.4)') lc_year
      lndname = trim(dir_landdata)//'/soil/'//trim(cyear)//'/lake_soilc_patches.nc'
      CALL ncio_read_vector (lndname, 'lake_soilc_patches', nl_soil, landpatch, &
         methane_lake_soilc_srf, defval = 0._r8)
#else
      IF (p_is_worker .and. numpatch > 0) THEN
         DO ipatch = 1, numpatch
            IF (patchtype(ipatch) == PATCHTYPE_LAKE) THEN
               methane_lake_soilc_srf(:,ipatch) = carbon_per_kg_om * max(OM_density(1:nl_soil,ipatch), 0._r8)
            ENDIF
         ENDDO
      ENDIF
#endif

   END SUBROUTINE load_methane_lake_soilc_srf

   SUBROUTINE allocate_methane_state (numpatch)
      USE MOD_MT_Methane_Const, only: DEF_METHANE
      integer, intent(in) :: numpatch
      ! Keep zero-length arrays allocated on IO/control ranks.  Methane
      ! restart reads use collective vector I/O; if numpatch==0 ranks return
      ! here, workers enter ncio_read_vector while IO ranks skip it and MPI
      ! can deadlock during global startup.

      allocate (net_methane                 (numpatch)); net_methane            (:) = 0._r8
      allocate (co2_hr_methane              (numpatch)); co2_hr_methane         (:) = 0._r8
      allocate (methane_prod_depth      (nl_soil,numpatch)); methane_prod_depth       (:,:) = 0._r8
      allocate (o2_decomp_depth     (nl_soil,numpatch)); o2_decomp_depth      (:,:) = 0._r8
      allocate (co2_decomp_depth    (nl_soil,numpatch)); co2_decomp_depth     (:,:) = 0._r8
      allocate (methane_oxid_depth      (nl_soil,numpatch)); methane_oxid_depth       (:,:) = 0._r8
      allocate (o2_oxid_depth       (nl_soil,numpatch)); o2_oxid_depth        (:,:) = 0._r8
      allocate (co2_oxid_depth      (nl_soil,numpatch)); co2_oxid_depth       (:,:) = 0._r8
      allocate (methane_aere_depth      (nl_soil,numpatch)); methane_aere_depth       (:,:) = 0._r8
      allocate (methane_tran_depth      (nl_soil,numpatch)); methane_tran_depth       (:,:) = 0._r8
      allocate (o2_aere_depth       (nl_soil,numpatch)); o2_aere_depth        (:,:) = 0._r8
      allocate (co2_aere_depth      (nl_soil,numpatch)); co2_aere_depth       (:,:) = 0._r8
      allocate (methane_ebul_depth      (nl_soil,numpatch)); methane_ebul_depth       (:,:) = 0._r8
      allocate (o2stress            (nl_soil,numpatch)); o2stress             (:,:) = 1.0_r8
      allocate (methane_stress           (nl_soil,numpatch)); methane_stress            (:,:) = 1.0_r8
      allocate (methane_surf_flux_tot           (numpatch)); methane_surf_flux_tot      (:) = 0._r8
      allocate (methane_surf_flux_tot_phys      (numpatch)); methane_surf_flux_tot_phys (:) = 0._r8
      allocate (methane_surf_aere               (numpatch)); methane_surf_aere          (:) = 0._r8
      allocate (methane_surf_ebul               (numpatch)); methane_surf_ebul          (:) = 0._r8
      allocate (methane_surf_diff               (numpatch)); methane_surf_diff          (:) = 0._r8
      allocate (methane_surf_diff_phys          (numpatch)); methane_surf_diff_phys     (:) = 0._r8
      allocate (methane_balance_residual       (numpatch)); methane_balance_residual  (:) = 0._r8
      allocate (methane_ch4_clip_credit        (numpatch)); methane_ch4_clip_credit   (:) = 0._r8
      allocate (o2_cap_loss                    (numpatch)); o2_cap_loss               (:) = 0._r8
      allocate (o2_cap_gain                    (numpatch)); o2_cap_gain               (:) = 0._r8
      allocate (methane_ebul_tot                (numpatch)); methane_ebul_tot           (:) = 0._r8
      allocate (methane_prod_tot                (numpatch)); methane_prod_tot           (:) = 0._r8
      allocate (methane_oxid_tot                (numpatch)); methane_oxid_tot           (:) = 0._r8
      allocate (co2_decomp_tot                  (numpatch)); co2_decomp_tot             (:) = 0._r8
      allocate (co2_oxid_tot                    (numpatch)); co2_oxid_tot               (:) = 0._r8
      allocate (co2_aere_tot                    (numpatch)); co2_aere_tot               (:) = 0._r8
      allocate (co2_net_tot                     (numpatch)); co2_net_tot                (:) = 0._r8

      allocate (totcol_methane                   (numpatch)); totcol_methane              (:) = 0._r8
      allocate (grnd_methane_cond               (numpatch)); grnd_methane_cond          (:) = DEF_METHANE%grnd_methane_cond_default
      ! Keep the combined cold-start O2 concentration consistent with the
      ! split unsat/sat pools.  Optional microbial pools may read conc_o2 on
      ! the first step before fsat_bef is valid.
      allocate (conc_o2             (nl_soil,numpatch)); conc_o2              (:,:) = 1.0_r8
      allocate (conc_methane            (nl_soil,numpatch)); conc_methane             (:,:) = 1.0e-6_r8
      !!!! --------------------------------------------------------------------------------------------------------

      !!!! --------------------------------------------------------------------------------------------------------
      !!!!                                         sum data (unsaturated / saturated)
      !!!! --------------------------------------------------------------------------------------------------------
      allocate (net_methane_unsat           (numpatch)); net_methane_unsat        (:)   = 0._r8
      allocate (net_methane_sat             (numpatch)); net_methane_sat          (:)   = 0._r8

      allocate (methane_prod_depth_unsat  (nl_soil,numpatch)); methane_prod_depth_unsat   (:,:) = 0._r8
      allocate (methane_prod_depth_sat    (nl_soil,numpatch)); methane_prod_depth_sat     (:,:) = 0._r8

      allocate (o2_decomp_depth_unsat (nl_soil,numpatch)); o2_decomp_depth_unsat  (:,:) = 0._r8
      allocate (o2_decomp_depth_sat   (nl_soil,numpatch)); o2_decomp_depth_sat    (:,:) = 0._r8
      allocate (co2_decomp_depth_unsat(nl_soil,numpatch)); co2_decomp_depth_unsat (:,:) = 0._r8
      allocate (co2_decomp_depth_sat  (nl_soil,numpatch)); co2_decomp_depth_sat   (:,:) = 0._r8

      allocate (methane_oxid_depth_unsat  (nl_soil,numpatch)); methane_oxid_depth_unsat   (:,:) = 0._r8
      allocate (methane_oxid_depth_sat    (nl_soil,numpatch)); methane_oxid_depth_sat     (:,:) = 0._r8

      allocate (o2_oxid_depth_unsat   (nl_soil,numpatch)); o2_oxid_depth_unsat    (:,:) = 0._r8
      allocate (o2_oxid_depth_sat     (nl_soil,numpatch)); o2_oxid_depth_sat      (:,:) = 0._r8
      allocate (co2_oxid_depth_unsat  (nl_soil,numpatch)); co2_oxid_depth_unsat   (:,:) = 0._r8
      allocate (co2_oxid_depth_sat    (nl_soil,numpatch)); co2_oxid_depth_sat     (:,:) = 0._r8

      allocate (methane_aere_depth_unsat  (nl_soil,numpatch)); methane_aere_depth_unsat   (:,:) = 0._r8
      allocate (methane_aere_depth_sat    (nl_soil,numpatch)); methane_aere_depth_sat     (:,:) = 0._r8

      allocate (methane_tran_depth_unsat  (nl_soil,numpatch)); methane_tran_depth_unsat   (:,:) = 0._r8
      allocate (methane_tran_depth_sat    (nl_soil,numpatch)); methane_tran_depth_sat     (:,:) = 0._r8

      allocate (o2_aere_depth_unsat   (nl_soil,numpatch)); o2_aere_depth_unsat    (:,:) = 0._r8
      allocate (o2_aere_depth_sat     (nl_soil,numpatch)); o2_aere_depth_sat      (:,:) = 0._r8
      allocate (co2_aere_depth_unsat  (nl_soil,numpatch)); co2_aere_depth_unsat   (:,:) = 0._r8
      allocate (co2_aere_depth_sat    (nl_soil,numpatch)); co2_aere_depth_sat     (:,:) = 0._r8

      allocate (methane_ebul_depth_unsat  (nl_soil,numpatch)); methane_ebul_depth_unsat   (:,:) = 0._r8
      allocate (methane_ebul_depth_sat    (nl_soil,numpatch)); methane_ebul_depth_sat     (:,:) = 0._r8

      allocate (o2stress_unsat        (nl_soil,numpatch)); o2stress_unsat         (:,:) = 1.0_r8
      allocate (o2stress_sat          (nl_soil,numpatch)); o2stress_sat           (:,:) = 1.0_r8

      allocate (methane_stress_unsat       (nl_soil,numpatch)); methane_stress_unsat        (:,:) = 1.0_r8
      allocate (methane_stress_sat         (nl_soil,numpatch)); methane_stress_sat          (:,:) = 1.0_r8

      allocate (methane_surf_flux_tot_unsat      (numpatch)); methane_surf_flux_tot_unsat  (:)   = 0._r8
      allocate (methane_surf_flux_tot_sat        (numpatch)); methane_surf_flux_tot_sat    (:)   = 0._r8

      allocate (methane_surf_aere_unsat          (numpatch)); methane_surf_aere_unsat      (:)   = 0._r8
      allocate (methane_surf_aere_sat            (numpatch)); methane_surf_aere_sat        (:)   = 0._r8

      allocate (methane_surf_ebul_unsat          (numpatch)); methane_surf_ebul_unsat      (:)   = 0._r8
      allocate (methane_surf_ebul_sat            (numpatch)); methane_surf_ebul_sat        (:)   = 0._r8

      allocate (methane_surf_diff_unsat          (numpatch)); methane_surf_diff_unsat      (:)   = 0._r8
      allocate (methane_surf_diff_sat            (numpatch)); methane_surf_diff_sat        (:)   = 0._r8
      allocate (methane_surf_diff_phys_unsat     (numpatch)); methane_surf_diff_phys_unsat (:)   = 0._r8
      allocate (methane_surf_diff_phys_sat       (numpatch)); methane_surf_diff_phys_sat   (:)   = 0._r8

      allocate (methane_ebul_tot_unsat           (numpatch)); methane_ebul_tot_unsat       (:)   = 0._r8
      allocate (methane_ebul_tot_sat             (numpatch)); methane_ebul_tot_sat         (:)   = 0._r8

      allocate (methane_prod_tot_unsat           (numpatch)); methane_prod_tot_unsat       (:)   = 0._r8
      allocate (methane_prod_tot_sat             (numpatch)); methane_prod_tot_sat         (:)   = 0._r8

      allocate (methane_oxid_tot_unsat           (numpatch)); methane_oxid_tot_unsat       (:)   = 0._r8
      allocate (methane_oxid_tot_sat             (numpatch)); methane_oxid_tot_sat         (:)   = 0._r8
      allocate (co2_decomp_tot_unsat             (numpatch)); co2_decomp_tot_unsat         (:)   = 0._r8
      allocate (co2_decomp_tot_sat               (numpatch)); co2_decomp_tot_sat           (:)   = 0._r8
      allocate (co2_oxid_tot_unsat               (numpatch)); co2_oxid_tot_unsat           (:)   = 0._r8
      allocate (co2_oxid_tot_sat                 (numpatch)); co2_oxid_tot_sat             (:)   = 0._r8
      allocate (co2_net_tot_unsat                (numpatch)); co2_net_tot_unsat            (:)   = 0._r8
      allocate (co2_net_tot_sat                  (numpatch)); co2_net_tot_sat              (:)   = 0._r8

      allocate (totcol_methane_unsat              (numpatch)); totcol_methane_unsat          (:)   = 0._r8
      allocate (totcol_methane_sat                (numpatch)); totcol_methane_sat            (:)   = 0._r8

      allocate (grnd_methane_cond_unsat          (numpatch)); grnd_methane_cond_unsat      (:)   = DEF_METHANE%grnd_methane_cond_default
      allocate (grnd_methane_cond_sat            (numpatch)); grnd_methane_cond_sat        (:)   = DEF_METHANE%grnd_methane_cond_default

      ! Physical cold-start concentrations for restart files without Methane state.
      allocate (conc_o2_unsat        (nl_soil,numpatch)); conc_o2_unsat           (:,:) = 1.0_r8
      allocate (conc_o2_sat          (nl_soil,numpatch)); conc_o2_sat             (:,:) = 1.0_r8

	      allocate (conc_methane_unsat       (nl_soil,numpatch)); conc_methane_unsat          (:,:) = 1.0e-6_r8
	      allocate (conc_methane_sat         (nl_soil,numpatch)); conc_methane_sat            (:,:) = 1.0e-6_r8
	      allocate (conc_o2_unsat_component      (nl_soil,N_METHANE_COMP,numpatch)); conc_o2_unsat_component      = 1.0_r8
	      allocate (conc_o2_sat_component        (nl_soil,N_METHANE_COMP,numpatch)); conc_o2_sat_component        = 1.0_r8
	      allocate (conc_methane_unsat_component (nl_soil,N_METHANE_COMP,numpatch)); conc_methane_unsat_component = 1.0e-6_r8
	      allocate (conc_methane_sat_component   (nl_soil,N_METHANE_COMP,numpatch)); conc_methane_sat_component   = 1.0e-6_r8
	      !!!! --------------------------------------------------------------------------------------------------------

	      !!!! --------------------------------------------------------------------------------------------------------
	      !!!!                                         lake data (CTSM alignment)
	      !!!! --------------------------------------------------------------------------------------------------------
	      allocate (methane_prod_depth_lake  (nl_soil,numpatch)); methane_prod_depth_lake     (:,:) = 0._r8
	      allocate (methane_oxid_depth_lake  (nl_soil,numpatch)); methane_oxid_depth_lake     (:,:) = 0._r8
	      allocate (methane_ebul_depth_lake  (nl_soil,numpatch)); methane_ebul_depth_lake     (:,:) = 0._r8
	      allocate (co2_decomp_depth_lake(nl_soil,numpatch)); co2_decomp_depth_lake     (:,:) = 0._r8
	      allocate (co2_oxid_depth_lake  (nl_soil,numpatch)); co2_oxid_depth_lake       (:,:) = 0._r8
	      allocate (methane_surf_ebul_lake          (numpatch)); methane_surf_ebul_lake       (:)   = 0._r8
	      allocate (methane_surf_diff_lake          (numpatch)); methane_surf_diff_lake       (:)   = 0._r8
	      allocate (methane_surf_flux_tot_lake      (numpatch)); methane_surf_flux_tot_lake   (:)   = 0._r8
	      allocate (methane_prod_tot_lake           (numpatch)); methane_prod_tot_lake        (:)   = 0._r8
	      allocate (methane_oxid_tot_lake           (numpatch)); methane_oxid_tot_lake        (:)   = 0._r8
	      allocate (methane_ebul_tot_lake           (numpatch)); methane_ebul_tot_lake        (:)   = 0._r8
	      allocate (co2_decomp_tot_lake           (numpatch)); co2_decomp_tot_lake        (:) = 0._r8
	      allocate (co2_oxid_tot_lake             (numpatch)); co2_oxid_tot_lake          (:) = 0._r8
	      allocate (co2_net_tot_lake              (numpatch)); co2_net_tot_lake           (:) = 0._r8
	      allocate (totcol_methane_lake             (numpatch)); totcol_methane_lake          (:)   = 0._r8
	      allocate (grnd_methane_cond_lake          (numpatch)); grnd_methane_cond_lake       (:)   = DEF_METHANE%grnd_methane_cond_default
	      allocate (conc_o2_lake             (nl_soil,numpatch)); conc_o2_lake                (:,:) = 1.0_r8
	      allocate (conc_methane_lake        (nl_soil,numpatch)); conc_methane_lake           (:,:) = 0._r8
	      allocate (lake_water_ch4_stock            (numpatch)); lake_water_ch4_stock         (:)   = 0._r8
	      allocate (lake_water_o2_stock             (numpatch)); lake_water_o2_stock          (:)   = 0._r8
	      allocate (lake_frozen_ch4_stock           (numpatch)); lake_frozen_ch4_stock        (:)   = 0._r8
	      allocate (lake_frozen_o2_stock            (numpatch)); lake_frozen_o2_stock         (:)   = 0._r8
	      allocate (lake_liquid_fraction_prev       (numpatch)); lake_liquid_fraction_prev    (:)   = spval
	      allocate (lake_water_ch4_oxid             (numpatch)); lake_water_ch4_oxid          (:)   = 0._r8
	      allocate (lake_sed_ch4_flux                (numpatch)); lake_sed_ch4_flux             (:)   = 0._r8
	      allocate (lake_sed_o2_flux                 (numpatch)); lake_sed_o2_flux              (:)   = 0._r8
	      allocate (lake_air_o2_flux                 (numpatch)); lake_air_o2_flux              (:)   = 0._r8
	      !!!! --------------------------------------------------------------------------------------------------------

	      allocate (c_atm                     (3,numpatch)); c_atm                (:,:) = 0._r8
      allocate (forc_pmethanem                  (numpatch)); forc_pmethanem             (:) = 0._r8
	      allocate (layer_sat_lag       (nl_soil,numpatch)); layer_sat_lag        (:,:) = spval
	      allocate (layer_sat_lag_component(nl_soil,N_METHANE_COMP,numpatch)); layer_sat_lag_component = spval
      allocate (lake_soilc          (nl_soil,numpatch)); lake_soilc           (:,:) = 0._r8
      allocate (annavg_agnpp                (numpatch)); annavg_agnpp           (:) = 0._r8
      allocate (annavg_bgnpp                (numpatch)); annavg_bgnpp           (:) = 0._r8
      allocate (annavg_somhr                (numpatch)); annavg_somhr           (:) = 0._r8
      ! annavg_finrw is a prerequisite for the seasonal inundation factor (SIF).
      ! Use spval until a full annual accumulator has been completed; otherwise a
      ! cold start / missing CH4 restart is interpreted as annavg_finrw=0 and SIF
      ! suppresses first-year saturated CH4 production by mino2lim.
      allocate (annavg_finrw                (numpatch)); annavg_finrw           (:) = spval
      allocate (tempavg_agnpp               (numpatch)); tempavg_agnpp          (:) = 0._r8
      allocate (tempavg_bgnpp               (numpatch)); tempavg_bgnpp          (:) = 0._r8
      allocate (annsum_counter              (numpatch)); annsum_counter         (:) = 0._r8
      allocate (tempavg_somhr               (numpatch)); tempavg_somhr          (:) = 0._r8
	      allocate (tempavg_finrw               (numpatch)); tempavg_finrw          (:) = 0._r8
	      allocate (annavg_agnpp_component  (N_METHANE_COMP,numpatch)); annavg_agnpp_component   = 0._r8
	      allocate (annavg_bgnpp_component  (N_METHANE_COMP,numpatch)); annavg_bgnpp_component   = 0._r8
	      allocate (annavg_somhr_component  (N_METHANE_COMP,numpatch)); annavg_somhr_component   = 0._r8
	      allocate (annavg_finrw_component  (N_METHANE_COMP,numpatch)); annavg_finrw_component   = spval
	      allocate (tempavg_agnpp_component (N_METHANE_COMP,numpatch)); tempavg_agnpp_component  = 0._r8
	      allocate (tempavg_bgnpp_component (N_METHANE_COMP,numpatch)); tempavg_bgnpp_component  = 0._r8
	      allocate (annsum_counter_component(N_METHANE_COMP,numpatch)); annsum_counter_component = 0._r8
	      allocate (tempavg_somhr_component (N_METHANE_COMP,numpatch)); tempavg_somhr_component  = 0._r8
	      allocate (tempavg_finrw_component (N_METHANE_COMP,numpatch)); tempavg_finrw_component  = 0._r8

	      allocate (fsat_bef                    (numpatch)); fsat_bef               (:) = spval
	      allocate (finundated_lag              (numpatch)); finundated_lag         (:) = spval
	      allocate (fsat_bef_component      (N_METHANE_COMP,numpatch)); fsat_bef_component       = spval
	      allocate (finundated_lag_component(N_METHANE_COMP,numpatch)); finundated_lag_component = spval
	      allocate (rice_fraction_prev      (numpatch)); rice_fraction_prev = 0._r8
      allocate (methane_dfsat_tot               (numpatch)); methane_dfsat_tot          (:) = 0._r8
      allocate (f_h2osfc                        (numpatch)); f_h2osfc                   (:) = 0._r8
      allocate (methane_finundated              (numpatch)); methane_finundated         (:) = 0._r8
      allocate (methane_soil_finundated         (numpatch)); methane_soil_finundated    (:) = 0._r8
      allocate (methane_soil_zwt                (numpatch)); methane_soil_zwt           (:) = spval
      allocate (methane_surf_flux_wetland       (numpatch)); methane_surf_flux_wetland  (:) = 0._r8
      allocate (methane_surf_flux_soil          (numpatch)); methane_surf_flux_soil     (:) = 0._r8
      allocate (methane_surf_flux_lake          (numpatch)); methane_surf_flux_lake     (:) = 0._r8
      allocate (methane_surf_flux_rice          (numpatch)); methane_surf_flux_rice     (:) = 0._r8
	  allocate (methane_surf_aere_soil(numpatch)); methane_surf_aere_soil = 0._r8
	  allocate (methane_surf_aere_rice(numpatch)); methane_surf_aere_rice = 0._r8
	  allocate (methane_surf_ebul_soil(numpatch)); methane_surf_ebul_soil = 0._r8
	  allocate (methane_surf_ebul_rice(numpatch)); methane_surf_ebul_rice = 0._r8
	  allocate (methane_surf_diff_soil(numpatch)); methane_surf_diff_soil = 0._r8
	  allocate (methane_surf_diff_rice(numpatch)); methane_surf_diff_rice = 0._r8
	  allocate (methane_prod_tot_soil(numpatch)); methane_prod_tot_soil = 0._r8
	  allocate (methane_prod_tot_rice(numpatch)); methane_prod_tot_rice = 0._r8
	  allocate (methane_oxid_tot_soil(numpatch)); methane_oxid_tot_soil = 0._r8
	  allocate (methane_oxid_tot_rice(numpatch)); methane_oxid_tot_rice = 0._r8
	      allocate (f_inund_levee_patch             (numpatch)); f_inund_levee_patch        (:) = 0._r8
	      allocate (f_inund_flood_patch             (numpatch)); f_inund_flood_patch        (:) = 0._r8
	      allocate (f_inund_flood_depth_patch       (numpatch)); f_inund_flood_depth_patch  (:) = 0._r8
	      allocate (wetland_frac_per_patch           (numpatch)); wetland_frac_per_patch     (:) = 1._r8
	      allocate (biome_f_methane_patch            (numpatch)); biome_f_methane_patch      (:) = 0.20_r8  ! default = legacy DEF_METHANE%f_methane
	      allocate (biome_redoxlag_patch             (numpatch)); biome_redoxlag_patch       (:) = 30._r8   ! default = legacy DEF_METHANE%redoxlag (days)

	   END SUBROUTINE allocate_methane_state


	   SUBROUTINE init_methane_wetland_fraction_cache (numpatch)
	      ! Build a patch-local lookup for the wetland fraction of each parent
	      ! element/cell.  This is intentionally a static landdata cache:
	      ! routing/GIEMS provide an absolute/grid flood fraction, while methane
	      ! wetland tiles need a patch-relative inundation fraction.
	      USE MOD_LandPatch,           only: landpatch
	      USE MOD_Mesh,                only: numelm, mesh
	      USE MOD_Pixel,               only: pixel
	      USE MOD_Utils,               only: areaquad
	      USE MOD_SPMD_Task,           only: p_is_worker
	      USE MOD_Vars_TimeInvariants, only: patchtype
	      IMPLICIT NONE

	      integer, intent(in) :: numpatch
	      real(r8), allocatable :: elm_wet_area(:), elm_act_area(:)
	      integer :: ipatch, ie
	      real(r8) :: area

	      IF (.not. allocated(wetland_frac_per_patch)) THEN
	         allocate(wetland_frac_per_patch(numpatch))
	      ELSEIF (size(wetland_frac_per_patch) /= numpatch) THEN
	         deallocate(wetland_frac_per_patch)
	         allocate(wetland_frac_per_patch(numpatch))
	      ENDIF
	      wetland_frac_per_patch(:) = 1._r8

	      IF (.not. p_is_worker) RETURN
	      IF (numpatch <= 0 .or. numelm <= 0) RETURN
	      IF (.not. allocated(patchtype)) RETURN
	      IF (.not. allocated(landpatch%ielm)) RETURN

	      allocate(elm_wet_area(numelm), elm_act_area(numelm))
	      elm_wet_area(:) = 0._r8
	      elm_act_area(:) = 0._r8

	      DO ipatch = 1, min(numpatch, size(patchtype))
	         IF (ipatch > size(landpatch%ielm)) CYCLE
	         ie = landpatch%ielm(ipatch)
	         IF (ie < 1 .or. ie > numelm) CYCLE
	         area = methane_patch_area(ipatch)
	         IF (area <= 0._r8) CYCLE
	         IF (patchtype(ipatch) == 0 .or. patchtype(ipatch) == 2) THEN
	            elm_act_area(ie) = elm_act_area(ie) + area
	         ENDIF
	         IF (patchtype(ipatch) == 2) THEN
	            elm_wet_area(ie) = elm_wet_area(ie) + area
	         ENDIF
	      ENDDO

	      DO ipatch = 1, min(numpatch, size(patchtype))
	         IF (ipatch > size(landpatch%ielm)) CYCLE
	         ie = landpatch%ielm(ipatch)
	         IF (ie < 1 .or. ie > numelm) CYCLE
	         IF (elm_act_area(ie) > 0._r8) THEN
	            wetland_frac_per_patch(ipatch) = max(0._r8, min(1._r8, &
	               elm_wet_area(ie) / elm_act_area(ie)))
	         ELSE
	            wetland_frac_per_patch(ipatch) = 0._r8
	         ENDIF
	      ENDDO

	      deallocate(elm_wet_area, elm_act_area)

	   CONTAINS

	      real(r8) FUNCTION methane_patch_area (ip)
	         integer, intent(in) :: ip
	         integer :: ipxstt, ipxend, ipxl, ie_local

	         methane_patch_area = 0._r8
	         IF (ip < 1) RETURN
	         IF (ip > size(landpatch%ielm)) RETURN
	         IF (.not. allocated(landpatch%ipxstt)) RETURN
	         IF (.not. allocated(landpatch%ipxend)) RETURN
	         IF (ip > size(landpatch%ipxstt) .or. ip > size(landpatch%ipxend)) RETURN

	         ie_local = landpatch%ielm(ip)
	         IF (ie_local < 1 .or. ie_local > numelm) RETURN
	         ipxstt = landpatch%ipxstt(ip)
	         ipxend = landpatch%ipxend(ip)
	         IF (ipxstt == -1 .and. ipxend == -1) THEN
	            ipxstt = 1
	            ipxend = mesh(ie_local)%npxl
	         ENDIF
	         IF (ipxstt < 1 .or. ipxend < ipxstt) RETURN
	         IF (ipxend > mesh(ie_local)%npxl) RETURN

	         DO ipxl = ipxstt, ipxend
	            methane_patch_area = methane_patch_area + areaquad ( &
	               pixel%lat_s(mesh(ie_local)%ilat(ipxl)), &
	               pixel%lat_n(mesh(ie_local)%ilat(ipxl)), &
	               pixel%lon_w(mesh(ie_local)%ilon(ipxl)), &
	               pixel%lon_e(mesh(ie_local)%ilon(ipxl)) )
	         ENDDO
	         IF (landpatch%has_shared .and. allocated(landpatch%pctshared)) THEN
	            IF (ip <= size(landpatch%pctshared)) THEN
	               methane_patch_area = methane_patch_area * max(0._r8, landpatch%pctshared(ip))
	            ENDIF
	         ENDIF
	      END FUNCTION methane_patch_area

	   END SUBROUTINE init_methane_wetland_fraction_cache


	   SUBROUTINE deallocate_methane_state ()

      IF (allocated(net_methane)) deallocate (net_methane)
      IF (allocated(co2_hr_methane)) deallocate (co2_hr_methane)
      IF (allocated(methane_prod_depth)) deallocate (methane_prod_depth)
      IF (allocated(o2_decomp_depth)) deallocate (o2_decomp_depth)
      IF (allocated(methane_oxid_depth)) deallocate (methane_oxid_depth)
      IF (allocated(o2_oxid_depth)) deallocate (o2_oxid_depth)
      IF (allocated(methane_aere_depth)) deallocate (methane_aere_depth)
      IF (allocated(methane_tran_depth)) deallocate (methane_tran_depth)
      IF (allocated(o2_aere_depth)) deallocate (o2_aere_depth)
      IF (allocated(methane_ebul_depth)) deallocate (methane_ebul_depth)
      IF (allocated(o2stress)) deallocate (o2stress)
      IF (allocated(methane_stress)) deallocate (methane_stress)
      IF (allocated(methane_surf_flux_tot)) deallocate (methane_surf_flux_tot)
      IF (allocated(methane_surf_flux_tot_phys)) deallocate (methane_surf_flux_tot_phys)
      IF (allocated(methane_surf_aere)) deallocate (methane_surf_aere)
      IF (allocated(methane_surf_ebul)) deallocate (methane_surf_ebul)
      IF (allocated(methane_surf_diff)) deallocate (methane_surf_diff)
      IF (allocated(methane_surf_diff_phys)) deallocate (methane_surf_diff_phys)
      IF (allocated(methane_balance_residual)) deallocate (methane_balance_residual)
      IF (allocated(methane_ch4_clip_credit)) deallocate (methane_ch4_clip_credit)
      IF (allocated(o2_cap_loss)) deallocate (o2_cap_loss)
      IF (allocated(o2_cap_gain)) deallocate (o2_cap_gain)
      IF (allocated(methane_ebul_tot)) deallocate (methane_ebul_tot)
      IF (allocated(methane_prod_tot)) deallocate (methane_prod_tot)
      IF (allocated(methane_oxid_tot)) deallocate (methane_oxid_tot)
      IF (allocated(totcol_methane)) deallocate (totcol_methane)
      IF (allocated(grnd_methane_cond)) deallocate (grnd_methane_cond)
      IF (allocated(conc_o2)) deallocate (conc_o2)
      IF (allocated(conc_methane)) deallocate (conc_methane)
      !!!! --------------------------------------------------------------------------------------------------------

      !!!! --------------------------------------------------------------------------------------------------------
      !!!!                                         sum data (unsaturated / saturated)
      !!!! --------------------------------------------------------------------------------------------------------
      IF (allocated(net_methane_unsat)) deallocate (net_methane_unsat)
      IF (allocated(net_methane_sat)) deallocate (net_methane_sat)
      IF (allocated(methane_prod_depth_unsat)) deallocate (methane_prod_depth_unsat)
      IF (allocated(methane_prod_depth_sat)) deallocate (methane_prod_depth_sat)
      IF (allocated(o2_decomp_depth_unsat)) deallocate (o2_decomp_depth_unsat)
      IF (allocated(o2_decomp_depth_sat)) deallocate (o2_decomp_depth_sat)
      IF (allocated(methane_oxid_depth_unsat)) deallocate (methane_oxid_depth_unsat)
      IF (allocated(methane_oxid_depth_sat)) deallocate (methane_oxid_depth_sat)
      IF (allocated(o2_oxid_depth_unsat)) deallocate (o2_oxid_depth_unsat)
      IF (allocated(o2_oxid_depth_sat)) deallocate (o2_oxid_depth_sat)
      IF (allocated(methane_aere_depth_unsat)) deallocate (methane_aere_depth_unsat)
      IF (allocated(methane_aere_depth_sat)) deallocate (methane_aere_depth_sat)
      IF (allocated(methane_tran_depth_unsat)) deallocate (methane_tran_depth_unsat)
      IF (allocated(methane_tran_depth_sat)) deallocate (methane_tran_depth_sat)
      IF (allocated(o2_aere_depth_unsat)) deallocate (o2_aere_depth_unsat)
      IF (allocated(o2_aere_depth_sat)) deallocate (o2_aere_depth_sat)
      IF (allocated(methane_ebul_depth_unsat)) deallocate (methane_ebul_depth_unsat)
      IF (allocated(methane_ebul_depth_sat)) deallocate (methane_ebul_depth_sat)
      IF (allocated(o2stress_unsat)) deallocate (o2stress_unsat)
      IF (allocated(o2stress_sat)) deallocate (o2stress_sat)
      IF (allocated(methane_stress_unsat)) deallocate (methane_stress_unsat)
      IF (allocated(methane_stress_sat)) deallocate (methane_stress_sat)
      IF (allocated(methane_surf_flux_tot_unsat)) deallocate (methane_surf_flux_tot_unsat)
      IF (allocated(methane_surf_flux_tot_sat)) deallocate (methane_surf_flux_tot_sat)
      IF (allocated(methane_surf_aere_unsat)) deallocate (methane_surf_aere_unsat)
      IF (allocated(methane_surf_aere_sat)) deallocate (methane_surf_aere_sat)
      IF (allocated(methane_surf_ebul_unsat)) deallocate (methane_surf_ebul_unsat)
      IF (allocated(methane_surf_ebul_sat)) deallocate (methane_surf_ebul_sat)
      IF (allocated(methane_surf_diff_unsat)) deallocate (methane_surf_diff_unsat)
      IF (allocated(methane_surf_diff_sat)) deallocate (methane_surf_diff_sat)
      IF (allocated(methane_surf_diff_phys_unsat)) deallocate (methane_surf_diff_phys_unsat)
      IF (allocated(methane_surf_diff_phys_sat)) deallocate (methane_surf_diff_phys_sat)
      IF (allocated(methane_ebul_tot_unsat)) deallocate (methane_ebul_tot_unsat)
      IF (allocated(methane_ebul_tot_sat)) deallocate (methane_ebul_tot_sat)
      IF (allocated(methane_prod_tot_unsat)) deallocate (methane_prod_tot_unsat)
      IF (allocated(methane_prod_tot_sat)) deallocate (methane_prod_tot_sat)
      IF (allocated(methane_oxid_tot_unsat)) deallocate (methane_oxid_tot_unsat)
      IF (allocated(methane_oxid_tot_sat)) deallocate (methane_oxid_tot_sat)
      IF (allocated(totcol_methane_unsat)) deallocate (totcol_methane_unsat)
      IF (allocated(totcol_methane_sat)) deallocate (totcol_methane_sat)
      IF (allocated(grnd_methane_cond_unsat)) deallocate (grnd_methane_cond_unsat)
      IF (allocated(grnd_methane_cond_sat)) deallocate (grnd_methane_cond_sat)
      IF (allocated(conc_o2_unsat)) deallocate (conc_o2_unsat)
	      IF (allocated(conc_o2_sat)) deallocate (conc_o2_sat)
	      IF (allocated(conc_methane_unsat)) deallocate (conc_methane_unsat)
	      IF (allocated(conc_methane_sat)) deallocate (conc_methane_sat)
	      IF (allocated(conc_o2_unsat_component)) deallocate (conc_o2_unsat_component)
	      IF (allocated(conc_o2_sat_component)) deallocate (conc_o2_sat_component)
	      IF (allocated(conc_methane_unsat_component)) deallocate (conc_methane_unsat_component)
	      IF (allocated(conc_methane_sat_component)) deallocate (conc_methane_sat_component)
	      !!!! --------------------------------------------------------------------------------------------------------

	      !!!! --------------------------------------------------------------------------------------------------------
	      !!!!                                         lake data (CTSM alignment)
	      !!!! --------------------------------------------------------------------------------------------------------
	      IF (allocated(methane_prod_depth_lake)) deallocate (methane_prod_depth_lake)
	      IF (allocated(methane_oxid_depth_lake)) deallocate (methane_oxid_depth_lake)
	      IF (allocated(methane_ebul_depth_lake)) deallocate (methane_ebul_depth_lake)
	      IF (allocated(methane_surf_ebul_lake)) deallocate (methane_surf_ebul_lake)
	      IF (allocated(methane_surf_diff_lake)) deallocate (methane_surf_diff_lake)
	      IF (allocated(methane_surf_flux_tot_lake)) deallocate (methane_surf_flux_tot_lake)
	      IF (allocated(methane_prod_tot_lake)) deallocate (methane_prod_tot_lake)
	      IF (allocated(methane_oxid_tot_lake)) deallocate (methane_oxid_tot_lake)
	      IF (allocated(methane_ebul_tot_lake)) deallocate (methane_ebul_tot_lake)
	      IF (allocated(totcol_methane_lake)) deallocate (totcol_methane_lake)
	      IF (allocated(grnd_methane_cond_lake)) deallocate (grnd_methane_cond_lake)
	      IF (allocated(conc_o2_lake)) deallocate (conc_o2_lake)
	      IF (allocated(conc_methane_lake)) deallocate (conc_methane_lake)
	      IF (allocated(lake_water_ch4_stock)) deallocate (lake_water_ch4_stock)
	      IF (allocated(lake_water_o2_stock)) deallocate (lake_water_o2_stock)
	      IF (allocated(lake_frozen_ch4_stock)) deallocate (lake_frozen_ch4_stock)
	      IF (allocated(lake_frozen_o2_stock)) deallocate (lake_frozen_o2_stock)
	      IF (allocated(lake_liquid_fraction_prev)) deallocate (lake_liquid_fraction_prev)
	      IF (allocated(lake_water_ch4_oxid)) deallocate (lake_water_ch4_oxid)
	      IF (allocated(lake_sed_ch4_flux)) deallocate (lake_sed_ch4_flux)
	      IF (allocated(lake_sed_o2_flux)) deallocate (lake_sed_o2_flux)
	      IF (allocated(lake_air_o2_flux)) deallocate (lake_air_o2_flux)
	      !!!! --------------------------------------------------------------------------------------------------------

	      IF (allocated(c_atm)) deallocate (c_atm)
      IF (allocated(forc_pmethanem)) deallocate (forc_pmethanem)
      IF (allocated(layer_sat_lag)) deallocate (layer_sat_lag)
	  IF (allocated(layer_sat_lag_component)) deallocate (layer_sat_lag_component)
      IF (allocated(lake_soilc)) deallocate (lake_soilc)
      IF (allocated(annavg_agnpp)) deallocate (annavg_agnpp)
      IF (allocated(annavg_bgnpp)) deallocate (annavg_bgnpp)
      IF (allocated(annavg_somhr)) deallocate (annavg_somhr)
      IF (allocated(annavg_finrw)) deallocate (annavg_finrw)
      IF (allocated(tempavg_agnpp)) deallocate (tempavg_agnpp)
      IF (allocated(tempavg_bgnpp)) deallocate (tempavg_bgnpp)
      IF (allocated(annsum_counter)) deallocate (annsum_counter)
      IF (allocated(tempavg_somhr)) deallocate (tempavg_somhr)
      IF (allocated(tempavg_finrw)) deallocate (tempavg_finrw)
	  IF (allocated(annavg_agnpp_component)) deallocate (annavg_agnpp_component)
	  IF (allocated(annavg_bgnpp_component)) deallocate (annavg_bgnpp_component)
	  IF (allocated(annavg_somhr_component)) deallocate (annavg_somhr_component)
	  IF (allocated(annavg_finrw_component)) deallocate (annavg_finrw_component)
	  IF (allocated(tempavg_agnpp_component)) deallocate (tempavg_agnpp_component)
	  IF (allocated(tempavg_bgnpp_component)) deallocate (tempavg_bgnpp_component)
	  IF (allocated(annsum_counter_component)) deallocate (annsum_counter_component)
	  IF (allocated(tempavg_somhr_component)) deallocate (tempavg_somhr_component)
	  IF (allocated(tempavg_finrw_component)) deallocate (tempavg_finrw_component)
      IF (allocated(fsat_bef)) deallocate (fsat_bef)
      IF (allocated(finundated_lag)) deallocate (finundated_lag)
	  IF (allocated(fsat_bef_component)) deallocate (fsat_bef_component)
	  IF (allocated(finundated_lag_component)) deallocate (finundated_lag_component)
	  IF (allocated(rice_fraction_prev)) deallocate (rice_fraction_prev)
      IF (allocated(methane_dfsat_tot)) deallocate (methane_dfsat_tot)
	      IF (allocated(f_h2osfc)) deallocate (f_h2osfc)
	      IF (allocated(methane_finundated)) deallocate (methane_finundated)
	      IF (allocated(methane_soil_finundated)) deallocate (methane_soil_finundated)
	      IF (allocated(methane_soil_zwt)) deallocate (methane_soil_zwt)
	      IF (allocated(methane_surf_flux_wetland)) deallocate (methane_surf_flux_wetland)
	      IF (allocated(methane_surf_flux_soil)) deallocate (methane_surf_flux_soil)
	      IF (allocated(methane_surf_flux_lake)) deallocate (methane_surf_flux_lake)
	      IF (allocated(methane_surf_flux_rice)) deallocate (methane_surf_flux_rice)
	      IF (allocated(methane_surf_aere_soil)) deallocate (methane_surf_aere_soil)
	      IF (allocated(methane_surf_aere_rice)) deallocate (methane_surf_aere_rice)
	      IF (allocated(methane_surf_ebul_soil)) deallocate (methane_surf_ebul_soil)
	      IF (allocated(methane_surf_ebul_rice)) deallocate (methane_surf_ebul_rice)
	      IF (allocated(methane_surf_diff_soil)) deallocate (methane_surf_diff_soil)
	      IF (allocated(methane_surf_diff_rice)) deallocate (methane_surf_diff_rice)
	      IF (allocated(methane_prod_tot_soil)) deallocate (methane_prod_tot_soil)
	      IF (allocated(methane_prod_tot_rice)) deallocate (methane_prod_tot_rice)
	      IF (allocated(methane_oxid_tot_soil)) deallocate (methane_oxid_tot_soil)
	      IF (allocated(methane_oxid_tot_rice)) deallocate (methane_oxid_tot_rice)
	      IF (allocated(f_inund_levee_patch)) deallocate (f_inund_levee_patch)
	      IF (allocated(f_inund_flood_patch)) deallocate (f_inund_flood_patch)
	      IF (allocated(f_inund_flood_depth_patch)) deallocate (f_inund_flood_depth_patch)
	      IF (allocated(wetland_frac_per_patch)) deallocate (wetland_frac_per_patch)
	      IF (allocated(biome_f_methane_patch)) deallocate (biome_f_methane_patch)
	      IF (allocated(biome_redoxlag_patch)) deallocate (biome_redoxlag_patch)
      IF (allocated(co2_decomp_depth       )) deallocate (co2_decomp_depth       )
      IF (allocated(co2_oxid_depth         )) deallocate (co2_oxid_depth         )
      IF (allocated(co2_aere_depth         )) deallocate (co2_aere_depth         )
      IF (allocated(co2_decomp_tot         )) deallocate (co2_decomp_tot         )
      IF (allocated(co2_oxid_tot           )) deallocate (co2_oxid_tot           )
      IF (allocated(co2_aere_tot           )) deallocate (co2_aere_tot           )
      IF (allocated(co2_net_tot            )) deallocate (co2_net_tot            )
      IF (allocated(co2_decomp_depth_unsat )) deallocate (co2_decomp_depth_unsat )
      IF (allocated(co2_decomp_depth_sat   )) deallocate (co2_decomp_depth_sat   )
      IF (allocated(co2_oxid_depth_unsat   )) deallocate (co2_oxid_depth_unsat   )
      IF (allocated(co2_oxid_depth_sat     )) deallocate (co2_oxid_depth_sat     )
      IF (allocated(co2_aere_depth_unsat   )) deallocate (co2_aere_depth_unsat   )
      IF (allocated(co2_aere_depth_sat     )) deallocate (co2_aere_depth_sat     )
      IF (allocated(co2_decomp_tot_unsat   )) deallocate (co2_decomp_tot_unsat   )
      IF (allocated(co2_decomp_tot_sat     )) deallocate (co2_decomp_tot_sat     )
      IF (allocated(co2_oxid_tot_unsat     )) deallocate (co2_oxid_tot_unsat     )
      IF (allocated(co2_oxid_tot_sat       )) deallocate (co2_oxid_tot_sat       )
      IF (allocated(co2_net_tot_unsat      )) deallocate (co2_net_tot_unsat      )
      IF (allocated(co2_net_tot_sat        )) deallocate (co2_net_tot_sat        )
      IF (allocated(co2_decomp_depth_lake  )) deallocate (co2_decomp_depth_lake  )
      IF (allocated(co2_oxid_depth_lake    )) deallocate (co2_oxid_depth_lake    )
      IF (allocated(co2_decomp_tot_lake    )) deallocate (co2_decomp_tot_lake    )
	      IF (allocated(co2_oxid_tot_lake      )) deallocate (co2_oxid_tot_lake      )
	      IF (allocated(co2_net_tot_lake       )) deallocate (co2_net_tot_lake       )
	      IF (allocated(methane_lake_substep_acc2d)) deallocate (methane_lake_substep_acc2d)
	      IF (allocated(methane_lake_substep_acc1d)) deallocate (methane_lake_substep_acc1d)
	      methane_lake_substep_cached_ipatch = -1
	      methane_lake_substep_next_isub = 1

	   END SUBROUTINE deallocate_methane_state


	   SUBROUTINE handle_methane_dry_lake_substep (ipatch, substep_dt)
	      integer,  intent(in) :: ipatch
	      real(r8), intent(in) :: substep_dt

	      real(r8) :: drydown_ch4_stock
	      real(r8) :: drydown_o2_stock
	      real(r8) :: drydown_ch4_flux

	      IF (.not. allocated(lake_water_ch4_stock)) RETURN
	      IF (ipatch < lbound(lake_water_ch4_stock,1) .or. &
	          ipatch > ubound(lake_water_ch4_stock,1)) RETURN

	      ! A dynamic dry lake follows the host soil-ground branch, so none of
	      ! the previous wet-lake process rates may survive into this substep.
	      ! Reset first; the one-time drydown exports are installed below.
	      CALL reset_methane_inactive_lake_diagnostics(ipatch)

	      ! Do not debit a prognostic inventory when no positive integration
	      ! interval exists to carry the compensating boundary flux.
	      IF (substep_dt <= 0._r8) RETURN

	      drydown_ch4_stock = max(lake_water_ch4_stock(ipatch), 0._r8) + &
	         max(lake_frozen_ch4_stock(ipatch), 0._r8)
	      drydown_o2_stock = max(lake_water_o2_stock(ipatch), 0._r8) + &
	         max(lake_frozen_o2_stock(ipatch), 0._r8)

	      lake_water_ch4_stock(ipatch)  = 0._r8
	      lake_frozen_ch4_stock(ipatch) = 0._r8
	      lake_water_o2_stock(ipatch)   = 0._r8
	      lake_frozen_o2_stock(ipatch)  = 0._r8

	      ! The sediment column retains a physical carrier.  The generic lake
	      ! inventory is redundant, so remove only the exported water/ice stock
	      ! and leave the sediment state and concentrations untouched.
	      totcol_methane(ipatch)       = totcol_methane_lake(ipatch)
	      totcol_methane_sat(ipatch)   = totcol_methane_lake(ipatch)
	      totcol_methane_unsat(ipatch) = 0._r8

	      IF (drydown_ch4_stock > 0._r8) THEN
	         drydown_ch4_flux = drydown_ch4_stock / substep_dt
	         methane_surf_diff(ipatch)          = drydown_ch4_flux
	         methane_surf_diff_phys(ipatch)     = drydown_ch4_flux
	         methane_surf_flux_tot(ipatch)      = drydown_ch4_flux
	         methane_surf_flux_tot_phys(ipatch) = drydown_ch4_flux
	         methane_surf_diff_lake(ipatch)     = drydown_ch4_flux
	         methane_surf_flux_tot_lake(ipatch) = drydown_ch4_flux
	         methane_surf_flux_lake(ipatch)     = drydown_ch4_flux
	      ENDIF

	      IF (drydown_o2_stock > 0._r8) &
	         lake_air_o2_flux(ipatch) = drydown_o2_stock / substep_dt

	      ! Force a deterministic water-column cold start after rewetting.  The
	      ! exported water/ice inventory must not be reconstructed from stale
	      ! inundation or phase-partition memory.
	      fsat_bef(ipatch)                   = spval
	      finundated_lag(ipatch)             = spval
	      layer_sat_lag(:,ipatch)            = spval
	      lake_liquid_fraction_prev(ipatch)  = spval

	   END SUBROUTINE handle_methane_dry_lake_substep


	   SUBROUTINE reset_methane_inactive_lake_diagnostics (ipatch)
	      integer, intent(in) :: ipatch

	      IF (.not. allocated(net_methane)) RETURN
	      IF (ipatch < lbound(net_methane,1) .or. &
	          ipatch > ubound(net_methane,1)) RETURN

	      ! Clear only instantaneous process and patch diagnostics.  Lake water,
	      ! frozen-water, and sediment inventories/concentrations are prognostic
	      ! and must survive a temporarily inactive lake-production pathway.

	      net_methane(ipatch) = 0._r8
	      co2_hr_methane(ipatch) = 0._r8
	      methane_prod_depth(:,ipatch) = 0._r8
	      o2_decomp_depth(:,ipatch) = 0._r8
	      co2_decomp_depth(:,ipatch) = 0._r8
	      methane_oxid_depth(:,ipatch) = 0._r8
	      o2_oxid_depth(:,ipatch) = 0._r8
	      co2_oxid_depth(:,ipatch) = 0._r8
	      methane_aere_depth(:,ipatch) = 0._r8
	      methane_tran_depth(:,ipatch) = 0._r8
	      o2_aere_depth(:,ipatch) = 0._r8
	      co2_aere_depth(:,ipatch) = 0._r8
	      methane_ebul_depth(:,ipatch) = 0._r8
	      o2stress(:,ipatch) = 0._r8
	      methane_stress(:,ipatch) = 0._r8

	      methane_surf_flux_tot(ipatch) = 0._r8
	      methane_surf_flux_tot_phys(ipatch) = 0._r8
	      methane_surf_aere(ipatch) = 0._r8
	      methane_surf_ebul(ipatch) = 0._r8
	      methane_surf_diff(ipatch) = 0._r8
	      methane_surf_diff_phys(ipatch) = 0._r8
	      methane_balance_residual(ipatch) = 0._r8
	      methane_ch4_clip_credit(ipatch) = 0._r8
	      o2_cap_loss(ipatch) = 0._r8
	      o2_cap_gain(ipatch) = 0._r8
	      methane_ebul_tot(ipatch) = 0._r8
	      methane_prod_tot(ipatch) = 0._r8
	      methane_oxid_tot(ipatch) = 0._r8
	      co2_decomp_tot(ipatch) = 0._r8
	      co2_oxid_tot(ipatch) = 0._r8
	      co2_aere_tot(ipatch) = 0._r8
	      co2_net_tot(ipatch) = 0._r8
	      ! Conductances are restart state, not process-rate diagnostics.  Keep
	      ! every inactive/dry reset inside the committed-restart write domain:
	      ! the strict reader requires all four values to remain finite and > 0.
	      grnd_methane_cond(ipatch) = DEF_METHANE%grnd_methane_cond_default

	      net_methane_unsat(ipatch) = 0._r8
	      net_methane_sat(ipatch) = 0._r8
	      methane_prod_depth_unsat(:,ipatch) = 0._r8
	      methane_prod_depth_sat(:,ipatch) = 0._r8
	      o2_decomp_depth_unsat(:,ipatch) = 0._r8
	      o2_decomp_depth_sat(:,ipatch) = 0._r8
	      co2_decomp_depth_unsat(:,ipatch) = 0._r8
	      co2_decomp_depth_sat(:,ipatch) = 0._r8
	      methane_oxid_depth_unsat(:,ipatch) = 0._r8
	      methane_oxid_depth_sat(:,ipatch) = 0._r8
	      o2_oxid_depth_unsat(:,ipatch) = 0._r8
	      o2_oxid_depth_sat(:,ipatch) = 0._r8
	      co2_oxid_depth_unsat(:,ipatch) = 0._r8
	      co2_oxid_depth_sat(:,ipatch) = 0._r8
	      methane_aere_depth_unsat(:,ipatch) = 0._r8
	      methane_aere_depth_sat(:,ipatch) = 0._r8
	      methane_tran_depth_unsat(:,ipatch) = 0._r8
	      methane_tran_depth_sat(:,ipatch) = 0._r8
	      o2_aere_depth_unsat(:,ipatch) = 0._r8
	      o2_aere_depth_sat(:,ipatch) = 0._r8
	      co2_aere_depth_unsat(:,ipatch) = 0._r8
	      co2_aere_depth_sat(:,ipatch) = 0._r8
	      methane_ebul_depth_unsat(:,ipatch) = 0._r8
	      methane_ebul_depth_sat(:,ipatch) = 0._r8
	      o2stress_unsat(:,ipatch) = 0._r8
	      o2stress_sat(:,ipatch) = 0._r8
	      methane_stress_unsat(:,ipatch) = 0._r8
	      methane_stress_sat(:,ipatch) = 0._r8
	      methane_surf_flux_tot_unsat(ipatch) = 0._r8
	      methane_surf_flux_tot_sat(ipatch) = 0._r8
	      methane_surf_aere_unsat(ipatch) = 0._r8
	      methane_surf_aere_sat(ipatch) = 0._r8
	      methane_surf_ebul_unsat(ipatch) = 0._r8
	      methane_surf_ebul_sat(ipatch) = 0._r8
	      methane_surf_diff_unsat(ipatch) = 0._r8
	      methane_surf_diff_sat(ipatch) = 0._r8
	      methane_surf_diff_phys_unsat(ipatch) = 0._r8
	      methane_surf_diff_phys_sat(ipatch) = 0._r8
	      methane_ebul_tot_unsat(ipatch) = 0._r8
	      methane_ebul_tot_sat(ipatch) = 0._r8
	      methane_prod_tot_unsat(ipatch) = 0._r8
	      methane_prod_tot_sat(ipatch) = 0._r8
	      methane_oxid_tot_unsat(ipatch) = 0._r8
	      methane_oxid_tot_sat(ipatch) = 0._r8
	      co2_decomp_tot_unsat(ipatch) = 0._r8
	      co2_decomp_tot_sat(ipatch) = 0._r8
	      co2_oxid_tot_unsat(ipatch) = 0._r8
	      co2_oxid_tot_sat(ipatch) = 0._r8
	      co2_net_tot_unsat(ipatch) = 0._r8
	      co2_net_tot_sat(ipatch) = 0._r8
	      grnd_methane_cond_unsat(ipatch) = DEF_METHANE%grnd_methane_cond_default
	      grnd_methane_cond_sat(ipatch) = DEF_METHANE%grnd_methane_cond_default

	      methane_prod_depth_lake(:,ipatch) = 0._r8
	      methane_oxid_depth_lake(:,ipatch) = 0._r8
	      methane_ebul_depth_lake(:,ipatch) = 0._r8
	      co2_decomp_depth_lake(:,ipatch) = 0._r8
	      co2_oxid_depth_lake(:,ipatch) = 0._r8
	      methane_surf_ebul_lake(ipatch) = 0._r8
	      methane_surf_diff_lake(ipatch) = 0._r8
	      methane_surf_flux_tot_lake(ipatch) = 0._r8
	      methane_prod_tot_lake(ipatch) = 0._r8
	      methane_oxid_tot_lake(ipatch) = 0._r8
	      methane_ebul_tot_lake(ipatch) = 0._r8
	      co2_decomp_tot_lake(ipatch) = 0._r8
	      co2_oxid_tot_lake(ipatch) = 0._r8
	      co2_net_tot_lake(ipatch) = 0._r8
	      grnd_methane_cond_lake(ipatch) = DEF_METHANE%grnd_methane_cond_default
	      lake_water_ch4_oxid(ipatch) = 0._r8
	      lake_sed_ch4_flux(ipatch) = 0._r8
	      lake_sed_o2_flux(ipatch) = 0._r8
	      lake_air_o2_flux(ipatch) = 0._r8

	      methane_dfsat_tot(ipatch) = 0._r8
	      f_h2osfc(ipatch) = 0._r8
	      methane_finundated(ipatch) = 0._r8
	      methane_soil_finundated(ipatch) = 0._r8
	      methane_soil_zwt(ipatch) = spval
	      methane_surf_flux_wetland(ipatch) = 0._r8
	      methane_surf_flux_soil(ipatch) = 0._r8
	      methane_surf_flux_lake(ipatch) = 0._r8
	      methane_surf_flux_rice(ipatch) = 0._r8
	      methane_surf_aere_soil(ipatch) = 0._r8
	      methane_surf_aere_rice(ipatch) = 0._r8
	      methane_surf_ebul_soil(ipatch) = 0._r8
	      methane_surf_ebul_rice(ipatch) = 0._r8
	      methane_surf_diff_soil(ipatch) = 0._r8
	      methane_surf_diff_rice(ipatch) = 0._r8
	      methane_prod_tot_soil(ipatch) = 0._r8
	      methane_prod_tot_rice(ipatch) = 0._r8
	      methane_oxid_tot_soil(ipatch) = 0._r8
	      methane_oxid_tot_rice(ipatch) = 0._r8

	   END SUBROUTINE reset_methane_inactive_lake_diagnostics


	   SUBROUTINE accumulate_methane_lake_substep_diagnostics (ipatch, substep_dt, isub, nsub)
	      integer,  intent(in) :: ipatch
	      real(r8), intent(in) :: substep_dt
	      integer,  intent(in) :: isub
	      integer,  intent(in) :: nsub

	      real(r8) :: total_dt

	      IF (nsub <= 1) RETURN
	      IF (substep_dt <= 0._r8) RETURN
	      IF (ipatch <= 0) RETURN

	      IF (.not. allocated(methane_lake_substep_acc2d)) THEN
	         allocate (methane_lake_substep_acc2d(nl_soil, methane_lake_substep_n2d))
	      ENDIF
	      IF (.not. allocated(methane_lake_substep_acc1d)) THEN
	         allocate (methane_lake_substep_acc1d(methane_lake_substep_n1d))
	      ENDIF

	      IF (isub <= 1 .or. ipatch /= methane_lake_substep_cached_ipatch .or. &
	          isub /= methane_lake_substep_next_isub) THEN
	         methane_lake_substep_acc2d(:,:) = 0._r8
	         methane_lake_substep_acc1d(:)   = 0._r8
	      ENDIF
	      methane_lake_substep_cached_ipatch = ipatch

	      CALL add2d( 1, methane_prod_depth(:,ipatch))
	      CALL add2d( 2, o2_decomp_depth(:,ipatch))
	      CALL add2d( 3, co2_decomp_depth(:,ipatch))
	      CALL add2d( 4, methane_oxid_depth(:,ipatch))
	      CALL add2d( 5, o2_oxid_depth(:,ipatch))
	      CALL add2d( 6, co2_oxid_depth(:,ipatch))
	      CALL add2d( 7, methane_aere_depth(:,ipatch))
	      CALL add2d( 8, methane_tran_depth(:,ipatch))
	      CALL add2d( 9, o2_aere_depth(:,ipatch))
	      CALL add2d(10, co2_aere_depth(:,ipatch))
	      CALL add2d(11, methane_ebul_depth(:,ipatch))
	      CALL add2d(12, o2stress(:,ipatch))
	      CALL add2d(13, methane_stress(:,ipatch))
	      CALL add2d(14, methane_prod_depth_unsat(:,ipatch))
	      CALL add2d(15, o2_decomp_depth_unsat(:,ipatch))
	      CALL add2d(16, co2_decomp_depth_unsat(:,ipatch))
	      CALL add2d(17, methane_oxid_depth_unsat(:,ipatch))
	      CALL add2d(18, o2_oxid_depth_unsat(:,ipatch))
	      CALL add2d(19, co2_oxid_depth_unsat(:,ipatch))
	      CALL add2d(20, methane_aere_depth_unsat(:,ipatch))
	      CALL add2d(21, methane_tran_depth_unsat(:,ipatch))
	      CALL add2d(22, o2_aere_depth_unsat(:,ipatch))
	      CALL add2d(23, co2_aere_depth_unsat(:,ipatch))
	      CALL add2d(24, methane_ebul_depth_unsat(:,ipatch))
	      CALL add2d(25, o2stress_unsat(:,ipatch))
	      CALL add2d(26, methane_stress_unsat(:,ipatch))
	      CALL add2d(27, methane_prod_depth_sat(:,ipatch))
	      CALL add2d(28, o2_decomp_depth_sat(:,ipatch))
	      CALL add2d(29, co2_decomp_depth_sat(:,ipatch))
	      CALL add2d(30, methane_oxid_depth_sat(:,ipatch))
	      CALL add2d(31, o2_oxid_depth_sat(:,ipatch))
	      CALL add2d(32, co2_oxid_depth_sat(:,ipatch))
	      CALL add2d(33, methane_aere_depth_sat(:,ipatch))
	      CALL add2d(34, methane_tran_depth_sat(:,ipatch))
	      CALL add2d(35, o2_aere_depth_sat(:,ipatch))
	      CALL add2d(36, co2_aere_depth_sat(:,ipatch))
	      CALL add2d(37, methane_ebul_depth_sat(:,ipatch))
	      CALL add2d(38, o2stress_sat(:,ipatch))
	      CALL add2d(39, methane_stress_sat(:,ipatch))
	      CALL add2d(40, methane_prod_depth_lake(:,ipatch))
	      CALL add2d(41, methane_oxid_depth_lake(:,ipatch))
	      CALL add2d(42, methane_ebul_depth_lake(:,ipatch))
	      CALL add2d(43, co2_decomp_depth_lake(:,ipatch))
	      CALL add2d(44, co2_oxid_depth_lake(:,ipatch))

	      CALL add1d( 1, net_methane(ipatch))
	      CALL add1d( 2, methane_surf_flux_tot(ipatch))
	      CALL add1d(46, methane_surf_flux_tot_phys(ipatch))
	      CALL add1d( 3, methane_surf_aere(ipatch))
	      CALL add1d( 4, methane_surf_ebul(ipatch))
	      CALL add1d( 5, methane_surf_diff(ipatch))
	      CALL add1d(43, methane_balance_residual(ipatch))
	      CALL add1d(47, methane_ch4_clip_credit(ipatch))
	      CALL add1d(44, o2_cap_loss(ipatch))
	      CALL add1d(45, o2_cap_gain(ipatch))
	      CALL add1d( 6, methane_ebul_tot(ipatch))
	      CALL add1d( 7, methane_prod_tot(ipatch))
	      CALL add1d( 8, methane_oxid_tot(ipatch))
	      CALL add1d( 9, co2_decomp_tot(ipatch))
	      CALL add1d(10, co2_oxid_tot(ipatch))
	      CALL add1d(11, co2_aere_tot(ipatch))
	      CALL add1d(12, co2_net_tot(ipatch))
	      CALL add1d(13, net_methane_unsat(ipatch))
	      CALL add1d(14, net_methane_sat(ipatch))
	      CALL add1d(15, methane_surf_flux_tot_unsat(ipatch))
	      CALL add1d(16, methane_surf_flux_tot_sat(ipatch))
	      CALL add1d(17, methane_surf_aere_unsat(ipatch))
	      CALL add1d(18, methane_surf_aere_sat(ipatch))
	      CALL add1d(19, methane_surf_ebul_unsat(ipatch))
	      CALL add1d(20, methane_surf_ebul_sat(ipatch))
	      CALL add1d(21, methane_surf_diff_unsat(ipatch))
	      CALL add1d(22, methane_surf_diff_sat(ipatch))
	      CALL add1d(23, methane_ebul_tot_unsat(ipatch))
	      CALL add1d(24, methane_ebul_tot_sat(ipatch))
	      CALL add1d(25, methane_prod_tot_unsat(ipatch))
	      CALL add1d(26, methane_prod_tot_sat(ipatch))
	      CALL add1d(27, methane_oxid_tot_unsat(ipatch))
	      CALL add1d(28, methane_oxid_tot_sat(ipatch))
	      CALL add1d(29, co2_decomp_tot_unsat(ipatch))
	      CALL add1d(30, co2_decomp_tot_sat(ipatch))
	      CALL add1d(31, co2_oxid_tot_unsat(ipatch))
	      CALL add1d(32, co2_oxid_tot_sat(ipatch))
	      CALL add1d(33, co2_net_tot_unsat(ipatch))
	      CALL add1d(34, co2_net_tot_sat(ipatch))
	      CALL add1d(35, methane_surf_ebul_lake(ipatch))
	      CALL add1d(36, methane_surf_diff_lake(ipatch))
	      CALL add1d(48, methane_surf_flux_tot_lake(ipatch))
	      CALL add1d(49, methane_surf_flux_lake(ipatch))
	      CALL add1d(37, methane_prod_tot_lake(ipatch))
	      CALL add1d(38, methane_oxid_tot_lake(ipatch))
	      CALL add1d(39, methane_ebul_tot_lake(ipatch))
	      CALL add1d(40, co2_decomp_tot_lake(ipatch))
	      CALL add1d(41, co2_oxid_tot_lake(ipatch))
	      CALL add1d(42, co2_net_tot_lake(ipatch))
	      CALL add1d(50, lake_water_ch4_oxid(ipatch))
	      CALL add1d(51, lake_sed_ch4_flux(ipatch))
	      CALL add1d(52, lake_sed_o2_flux(ipatch))
	      CALL add1d(53, lake_air_o2_flux(ipatch))
	      CALL add1d(54, grnd_methane_cond(ipatch))
	      CALL add1d(55, grnd_methane_cond_sat(ipatch))
	      CALL add1d(56, grnd_methane_cond_lake(ipatch))
	      CALL add1d(57, methane_surf_diff_phys(ipatch))

	      IF (isub == nsub) THEN
	         total_dt = substep_dt * real(nsub, r8)
	         IF (total_dt > 0._r8) THEN
	            CALL finish2d( 1, methane_prod_depth(:,ipatch))
	            CALL finish2d( 2, o2_decomp_depth(:,ipatch))
	            CALL finish2d( 3, co2_decomp_depth(:,ipatch))
	            CALL finish2d( 4, methane_oxid_depth(:,ipatch))
	            CALL finish2d( 5, o2_oxid_depth(:,ipatch))
	            CALL finish2d( 6, co2_oxid_depth(:,ipatch))
	            CALL finish2d( 7, methane_aere_depth(:,ipatch))
	            CALL finish2d( 8, methane_tran_depth(:,ipatch))
	            CALL finish2d( 9, o2_aere_depth(:,ipatch))
	            CALL finish2d(10, co2_aere_depth(:,ipatch))
	            CALL finish2d(11, methane_ebul_depth(:,ipatch))
	            CALL finish2d(12, o2stress(:,ipatch))
	            CALL finish2d(13, methane_stress(:,ipatch))
	            CALL finish2d(14, methane_prod_depth_unsat(:,ipatch))
	            CALL finish2d(15, o2_decomp_depth_unsat(:,ipatch))
	            CALL finish2d(16, co2_decomp_depth_unsat(:,ipatch))
	            CALL finish2d(17, methane_oxid_depth_unsat(:,ipatch))
	            CALL finish2d(18, o2_oxid_depth_unsat(:,ipatch))
	            CALL finish2d(19, co2_oxid_depth_unsat(:,ipatch))
	            CALL finish2d(20, methane_aere_depth_unsat(:,ipatch))
	            CALL finish2d(21, methane_tran_depth_unsat(:,ipatch))
	            CALL finish2d(22, o2_aere_depth_unsat(:,ipatch))
	            CALL finish2d(23, co2_aere_depth_unsat(:,ipatch))
	            CALL finish2d(24, methane_ebul_depth_unsat(:,ipatch))
	            CALL finish2d(25, o2stress_unsat(:,ipatch))
	            CALL finish2d(26, methane_stress_unsat(:,ipatch))
	            CALL finish2d(27, methane_prod_depth_sat(:,ipatch))
	            CALL finish2d(28, o2_decomp_depth_sat(:,ipatch))
	            CALL finish2d(29, co2_decomp_depth_sat(:,ipatch))
	            CALL finish2d(30, methane_oxid_depth_sat(:,ipatch))
	            CALL finish2d(31, o2_oxid_depth_sat(:,ipatch))
	            CALL finish2d(32, co2_oxid_depth_sat(:,ipatch))
	            CALL finish2d(33, methane_aere_depth_sat(:,ipatch))
	            CALL finish2d(34, methane_tran_depth_sat(:,ipatch))
	            CALL finish2d(35, o2_aere_depth_sat(:,ipatch))
	            CALL finish2d(36, co2_aere_depth_sat(:,ipatch))
	            CALL finish2d(37, methane_ebul_depth_sat(:,ipatch))
	            CALL finish2d(38, o2stress_sat(:,ipatch))
	            CALL finish2d(39, methane_stress_sat(:,ipatch))
	            CALL finish2d(40, methane_prod_depth_lake(:,ipatch))
	            CALL finish2d(41, methane_oxid_depth_lake(:,ipatch))
	            CALL finish2d(42, methane_ebul_depth_lake(:,ipatch))
	            CALL finish2d(43, co2_decomp_depth_lake(:,ipatch))
	            CALL finish2d(44, co2_oxid_depth_lake(:,ipatch))

	            CALL finish1d( 1, net_methane(ipatch))
	            CALL finish1d( 2, methane_surf_flux_tot(ipatch))
	            CALL finish1d(46, methane_surf_flux_tot_phys(ipatch))
	            CALL finish1d( 3, methane_surf_aere(ipatch))
	            CALL finish1d( 4, methane_surf_ebul(ipatch))
	            CALL finish1d( 5, methane_surf_diff(ipatch))
	            CALL finish1d(43, methane_balance_residual(ipatch))
	            CALL finish1d(47, methane_ch4_clip_credit(ipatch))
	            CALL finish1d(44, o2_cap_loss(ipatch))
	            CALL finish1d(45, o2_cap_gain(ipatch))
	            CALL finish1d( 6, methane_ebul_tot(ipatch))
	            CALL finish1d( 7, methane_prod_tot(ipatch))
	            CALL finish1d( 8, methane_oxid_tot(ipatch))
	            CALL finish1d( 9, co2_decomp_tot(ipatch))
	            CALL finish1d(10, co2_oxid_tot(ipatch))
	            CALL finish1d(11, co2_aere_tot(ipatch))
	            CALL finish1d(12, co2_net_tot(ipatch))
	            CALL finish1d(13, net_methane_unsat(ipatch))
	            CALL finish1d(14, net_methane_sat(ipatch))
	            CALL finish1d(15, methane_surf_flux_tot_unsat(ipatch))
	            CALL finish1d(16, methane_surf_flux_tot_sat(ipatch))
	            CALL finish1d(17, methane_surf_aere_unsat(ipatch))
	            CALL finish1d(18, methane_surf_aere_sat(ipatch))
	            CALL finish1d(19, methane_surf_ebul_unsat(ipatch))
	            CALL finish1d(20, methane_surf_ebul_sat(ipatch))
	            CALL finish1d(21, methane_surf_diff_unsat(ipatch))
	            CALL finish1d(22, methane_surf_diff_sat(ipatch))
	            CALL finish1d(23, methane_ebul_tot_unsat(ipatch))
	            CALL finish1d(24, methane_ebul_tot_sat(ipatch))
	            CALL finish1d(25, methane_prod_tot_unsat(ipatch))
	            CALL finish1d(26, methane_prod_tot_sat(ipatch))
	            CALL finish1d(27, methane_oxid_tot_unsat(ipatch))
	            CALL finish1d(28, methane_oxid_tot_sat(ipatch))
	            CALL finish1d(29, co2_decomp_tot_unsat(ipatch))
	            CALL finish1d(30, co2_decomp_tot_sat(ipatch))
	            CALL finish1d(31, co2_oxid_tot_unsat(ipatch))
	            CALL finish1d(32, co2_oxid_tot_sat(ipatch))
	            CALL finish1d(33, co2_net_tot_unsat(ipatch))
	            CALL finish1d(34, co2_net_tot_sat(ipatch))
	            CALL finish1d(35, methane_surf_ebul_lake(ipatch))
	            CALL finish1d(36, methane_surf_diff_lake(ipatch))
	            CALL finish1d(48, methane_surf_flux_tot_lake(ipatch))
	            CALL finish1d(49, methane_surf_flux_lake(ipatch))
	            CALL finish1d(37, methane_prod_tot_lake(ipatch))
	            CALL finish1d(38, methane_oxid_tot_lake(ipatch))
	            CALL finish1d(39, methane_ebul_tot_lake(ipatch))
	            CALL finish1d(40, co2_decomp_tot_lake(ipatch))
	            CALL finish1d(41, co2_oxid_tot_lake(ipatch))
	            CALL finish1d(42, co2_net_tot_lake(ipatch))
	            CALL finish1d(50, lake_water_ch4_oxid(ipatch))
	            CALL finish1d(51, lake_sed_ch4_flux(ipatch))
	            CALL finish1d(52, lake_sed_o2_flux(ipatch))
	            CALL finish1d(53, lake_air_o2_flux(ipatch))
	            CALL finish1d(54, grnd_methane_cond(ipatch))
	            CALL finish1d(55, grnd_methane_cond_sat(ipatch))
	            CALL finish1d(56, grnd_methane_cond_lake(ipatch))
	            CALL finish1d(57, methane_surf_diff_phys(ipatch))
	         ENDIF
	         methane_lake_substep_cached_ipatch = -1
	         methane_lake_substep_next_isub = 1
	      ELSE
	         methane_lake_substep_next_isub = isub + 1
	      ENDIF

	   CONTAINS
	      SUBROUTINE add2d (icol, var)
	         integer,  intent(in) :: icol
	         real(r8), intent(in) :: var(1:nl_soil)
	         methane_lake_substep_acc2d(:,icol) = methane_lake_substep_acc2d(:,icol) + var(:) * substep_dt
	      END SUBROUTINE add2d

	      SUBROUTINE add1d (icol, var)
	         integer,  intent(in) :: icol
	         real(r8), intent(in) :: var
	         methane_lake_substep_acc1d(icol) = methane_lake_substep_acc1d(icol) + var * substep_dt
	      END SUBROUTINE add1d

	      SUBROUTINE finish2d (icol, var)
	         integer,  intent(in)    :: icol
	         real(r8), intent(inout) :: var(1:nl_soil)
	         var(:) = methane_lake_substep_acc2d(:,icol) / total_dt
	      END SUBROUTINE finish2d

	      SUBROUTINE finish1d (icol, var)
	         integer,  intent(in)    :: icol
	         real(r8), intent(inout) :: var
	         var = methane_lake_substep_acc1d(icol) / total_dt
	      END SUBROUTINE finish1d
	   END SUBROUTINE accumulate_methane_lake_substep_diagnostics


	   !-------------------------------------------------------------------
	   ! compute_f_h2osfc — CLM5 microtopography-based prognostic surface water fraction
   !
   ! Ported from the Methane source restart flow.
   !   main/MOD_Runoff.F90 (line 332-462) / MOD_Runoff_h2osfc.F90 (line 247)
   !   main/MOD_SoilSnowHydrology.F90 (line 1307-1322)
   ! Maintains f_h2osfc(i) as patch-level diagnostic state, called by CoLMDRIVER
   ! immediately before methane_driver. Implementation is minimal-invasion:
   ! does NOT modify CoLM202X WATER_2014/WATER_VSF/Runoff physics. f_h2osfc
   ! is recomputed each step from current wdsrf + slpratio + microtopography
   ! params (DEF_METHANE_hydrology%slopemax/slopebeta).
   !
   ! Inputs:
   !   ipatch    — patch index
   !   slpratio  — slope ratio (from MOD_Vars_TimeInvariants)
   !   wdsrf     — surface water depth [mm] (from MOD_Vars_TimeVariables)
   ! Output:
   !   f_h2osfc(ipatch) is updated
   !-------------------------------------------------------------------
   SUBROUTINE compute_f_h2osfc (ipatch, slpratio_in, wdsrf_in)
      USE MOD_Vars_Global, only: PI
      USE MOD_MT_Methane_Const, only: DEF_METHANE_hydrology
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: slpratio_in      ! slope ratio [-]
      real(r8), intent(in) :: wdsrf_in         ! surface water depth [mm]

      real(r8) :: micro_sigma, sigma_mm, d, slope_angle, slope_arg
      real(r8) :: fd, dfdd, d_lo, d_hi, d_mid, f_mid
      integer  :: p
      logical  :: converged
      real(r8), parameter :: pondmin = 1.e-8_r8
      real(r8), parameter :: fd_tol  = 1.e-10_r8

      IF (.not. ieee_is_finite(wdsrf_in) .or. .not. ieee_is_finite(slpratio_in) .or. &
          wdsrf_in == spval .or. slpratio_in == spval .or. &
          wdsrf_in <= pondmin .or. abs(slpratio_in) >= 1.e30_r8 .or. &
          DEF_METHANE_hydrology%slopemax <= 0._r8 .or. &
          DEF_METHANE_hydrology%slopebeta >= 0._r8) THEN
         f_h2osfc(ipatch) = 0._r8
         RETURN
      END IF

      ! MOD_Vars_TimeInvariants defines slpratio as the dimensionless rise/run
      ! ratio, so the physical slope angle is unambiguously atan(slpratio).
      ! Guessing radians/degrees/percent from magnitude creates discontinuities
      ! for perfectly valid steep slopes (slpratio > 1).
      slope_arg = max(slpratio_in, 0._r8)
      slope_angle = atan(slope_arg)
      slope_angle = max(0._r8, min(0.5_r8*PI, slope_angle))
      micro_sigma = (slope_angle + &
                     DEF_METHANE_hydrology%slopemax**(1._r8/DEF_METHANE_hydrology%slopebeta) &
                    )**DEF_METHANE_hydrology%slopebeta
      micro_sigma = max(0._r8, min(DEF_METHANE_hydrology%slopemax, micro_sigma))
      sigma_mm = 1.0e3_r8 * micro_sigma   ! m -> mm

	      IF (sigma_mm > 1.e-3_r8) THEN
	         IF (wdsrf_in >= 10._r8 * sigma_mm) THEN
	            f_h2osfc(ipatch) = 1._r8
	            RETURN
	         ENDIF
	         ! Newton iteration for the CLM fill-and-spill microtopography
         ! relation:
         !   W(d) = 0.5*d*(1+erf(d/(sigma*sqrt(2))))
         !        + sigma/sqrt(2*pi)*exp(-d**2/(2*sigma**2))
         ! where W is grid-cell mean surface-water depth [mm].  The old
         ! implementation solved an unrelated pc threshold and then replaced
         ! d by wdsrf, which made wdsrf=0 diagnose f_h2osfc=0.5.
		         d = min(max(0._r8, wdsrf_in), 10._r8 * sigma_mm)
		         converged = .false.
	         DO p = 1, 20
	            fd = 0.5_r8 * d * (1.0_r8 + erf(d / (sigma_mm * sqrt(2.0_r8)))) &
	               + sigma_mm / sqrt(2.0_r8 * PI) &
	               * exp(-d**2 / (2.0_r8 * sigma_mm**2)) &
	               - wdsrf_in
	            dfdd = 0.5_r8 * (1.0_r8 + erf(d / (sigma_mm * sqrt(2.0_r8))))
	            IF (abs(fd) < fd_tol) THEN
	               converged = .true.
	               EXIT
	            ENDIF
	            IF (dfdd < 1.e-12_r8) EXIT
	            d = d - fd / dfdd
	            d = max(-10._r8 * sigma_mm, min(10._r8 * sigma_mm, d))
	         END DO
	         IF (.not. converged) THEN
	            ! Monotone fallback for rare Newton failures.  This avoids
	            ! silently using the last Newton iterate when the derivative is
	            ! tiny near the dry tail.
	            d_lo = -10._r8 * sigma_mm
	            d_hi =  10._r8 * sigma_mm
	            DO p = 1, 60
	               d_mid = 0.5_r8 * (d_lo + d_hi)
	               f_mid = 0.5_r8 * d_mid * (1.0_r8 + erf(d_mid / (sigma_mm * sqrt(2.0_r8)))) &
	                  + sigma_mm / sqrt(2.0_r8 * PI) &
	                  * exp(-d_mid**2 / (2.0_r8 * sigma_mm**2)) &
	                  - wdsrf_in
	               IF (abs(f_mid) < fd_tol) EXIT
	               IF (f_mid > 0._r8) THEN
	                  d_hi = d_mid
	               ELSE
	                  d_lo = d_mid
	               ENDIF
	            END DO
	            d = d_mid
	         ENDIF
         f_h2osfc(ipatch) = 0.5_r8 * (1.0_r8 + erf(d / (sigma_mm * sqrt(2.0_r8))))
         f_h2osfc(ipatch) = min(1._r8, max(0._r8, f_h2osfc(ipatch)))
      ELSE
         f_h2osfc(ipatch) = 0._r8
      END IF
   END SUBROUTINE compute_f_h2osfc


   SUBROUTINE write_methane_restart (file_restart, compress)
      USE MOD_LandPatch,     only: landpatch
      USE MOD_NetCDFVector, only: ncio_write_vector
      character(len=*), intent(in) :: file_restart
      integer,          intent(in) :: compress

      IF (.not. allocated(conc_methane)) RETURN

	      CALL ncio_write_vector (file_restart, 'ch4_conc_o2',          'soil', nl_soil, 'patch', landpatch, conc_o2,          compress)
	      CALL ncio_write_vector (file_restart, 'ch4_conc_methane',     'soil', nl_soil, 'patch', landpatch, conc_methane,     compress)
	      CALL ncio_write_vector (file_restart, 'ch4_totcol_methane',   'patch', landpatch, totcol_methane,       compress)
	      CALL ncio_write_vector (file_restart, 'ch4_grnd_methane_cond','patch', landpatch, grnd_methane_cond,    compress)
	      CALL ncio_write_vector (file_restart, 'ch4_conc_o2_unsat',    'soil', nl_soil, 'patch', landpatch, conc_o2_unsat,    compress)
	      CALL ncio_write_vector (file_restart, 'ch4_conc_o2_sat',      'soil', nl_soil, 'patch', landpatch, conc_o2_sat,      compress)
	      CALL ncio_write_vector (file_restart, 'ch4_conc_ch4_unsat',   'soil', nl_soil, &
	                              'patch', landpatch, conc_methane_unsat, compress)
	      CALL ncio_write_vector (file_restart, 'ch4_conc_ch4_sat',     'soil', nl_soil, &
	                              'patch', landpatch, conc_methane_sat,   compress)
	      CALL ncio_write_vector (file_restart, 'ch4_conc_o2_unsat_soil', 'soil', nl_soil, &
	                              'patch', landpatch, conc_o2_unsat_component(:,METHANE_COMP_SOIL,:), compress)
	      CALL ncio_write_vector (file_restart, 'ch4_conc_o2_unsat_rice', 'soil', nl_soil, &
	                              'patch', landpatch, conc_o2_unsat_component(:,METHANE_COMP_RICE,:), compress)
	      CALL ncio_write_vector (file_restart, 'ch4_conc_o2_sat_soil', 'soil', nl_soil, &
	                              'patch', landpatch, conc_o2_sat_component(:,METHANE_COMP_SOIL,:), compress)
	      CALL ncio_write_vector (file_restart, 'ch4_conc_o2_sat_rice', 'soil', nl_soil, &
	                              'patch', landpatch, conc_o2_sat_component(:,METHANE_COMP_RICE,:), compress)
	      CALL ncio_write_vector (file_restart, 'ch4_conc_ch4_unsat_soil', 'soil', nl_soil, &
	                              'patch', landpatch, conc_methane_unsat_component(:,METHANE_COMP_SOIL,:), compress)
	      CALL ncio_write_vector (file_restart, 'ch4_conc_ch4_unsat_rice', 'soil', nl_soil, &
	                              'patch', landpatch, conc_methane_unsat_component(:,METHANE_COMP_RICE,:), compress)
	      CALL ncio_write_vector (file_restart, 'ch4_conc_ch4_sat_soil', 'soil', nl_soil, &
	                              'patch', landpatch, conc_methane_sat_component(:,METHANE_COMP_SOIL,:), compress)
	      CALL ncio_write_vector (file_restart, 'ch4_conc_ch4_sat_rice', 'soil', nl_soil, &
	                              'patch', landpatch, conc_methane_sat_component(:,METHANE_COMP_RICE,:), compress)
	      CALL ncio_write_vector (file_restart, 'ch4_totcol_methane_unsat','patch', landpatch, totcol_methane_unsat, compress)
	      CALL ncio_write_vector (file_restart, 'ch4_totcol_methane_sat',  'patch', landpatch, totcol_methane_sat,   compress)
	      CALL ncio_write_vector (file_restart, 'ch4_grnd_methane_cond_unsat','patch', landpatch, grnd_methane_cond_unsat, compress)
	      CALL ncio_write_vector (file_restart, 'ch4_grnd_methane_cond_sat',  'patch', landpatch, grnd_methane_cond_sat,   compress)
		      CALL ncio_write_vector (file_restart, 'ch4_conc_o2_lake',     'soil', nl_soil, 'patch', landpatch, conc_o2_lake,     compress)
	      CALL ncio_write_vector (file_restart, 'ch4_conc_ch4_lake',    'soil', nl_soil, 'patch', landpatch, conc_methane_lake, compress)
	      CALL ncio_write_vector (file_restart, 'ch4_totcol_lake',      'patch', landpatch, totcol_methane_lake, compress)
	      CALL ncio_write_vector (file_restart, 'ch4_grnd_methane_cond_lake','patch', landpatch, grnd_methane_cond_lake, compress)
	      CALL ncio_write_vector (file_restart, 'ch4_lake_water_ch4_stock','patch', landpatch, lake_water_ch4_stock, compress)
	      CALL ncio_write_vector (file_restart, 'ch4_lake_water_o2_stock', 'patch', landpatch, lake_water_o2_stock,  compress)
	      CALL ncio_write_vector (file_restart, 'ch4_lake_frozen_ch4_stock','patch', landpatch, lake_frozen_ch4_stock, compress)
	      CALL ncio_write_vector (file_restart, 'ch4_lake_frozen_o2_stock', 'patch', landpatch, lake_frozen_o2_stock,  compress)
	      CALL ncio_write_vector (file_restart, 'ch4_lake_liquid_fraction_prev', 'patch', landpatch, &
	         lake_liquid_fraction_prev, compress)
	      CALL ncio_write_vector (file_restart, 'ch4_layer_sat_lag',    'soil', nl_soil, 'patch', landpatch, layer_sat_lag,    compress)
	      CALL ncio_write_vector (file_restart, 'ch4_layer_sat_lag_soil', 'soil', nl_soil, &
	                              'patch', landpatch, layer_sat_lag_component(:,METHANE_COMP_SOIL,:), compress)
	      CALL ncio_write_vector (file_restart, 'ch4_layer_sat_lag_rice', 'soil', nl_soil, &
	                              'patch', landpatch, layer_sat_lag_component(:,METHANE_COMP_RICE,:), compress)
	      CALL ncio_write_vector (file_restart, 'ch4_lake_soilc',       'soil', nl_soil, 'patch', landpatch, lake_soilc,       compress)
      CALL ncio_write_vector (file_restart, 'ch4_annavg_agnpp',     'patch', landpatch, annavg_agnpp,     compress)
      CALL ncio_write_vector (file_restart, 'ch4_annavg_bgnpp',     'patch', landpatch, annavg_bgnpp,     compress)
      CALL ncio_write_vector (file_restart, 'ch4_annavg_somhr',     'patch', landpatch, annavg_somhr,     compress)
      CALL ncio_write_vector (file_restart, 'ch4_annavg_finrw',     'patch', landpatch, annavg_finrw,     compress)
      CALL ncio_write_vector (file_restart, 'ch4_tempavg_agnpp',    'patch', landpatch, tempavg_agnpp,    compress)
      CALL ncio_write_vector (file_restart, 'ch4_tempavg_bgnpp',    'patch', landpatch, tempavg_bgnpp,    compress)
      CALL ncio_write_vector (file_restart, 'ch4_annsum_counter',   'patch', landpatch, annsum_counter,   compress)
      CALL ncio_write_vector (file_restart, 'ch4_tempavg_somhr',    'patch', landpatch, tempavg_somhr,    compress)
      CALL ncio_write_vector (file_restart, 'ch4_tempavg_finrw',    'patch', landpatch, tempavg_finrw,    compress)
	  CALL ncio_write_vector (file_restart, 'ch4_annavg_agnpp_soil', 'patch', landpatch, &
	     annavg_agnpp_component(METHANE_COMP_SOIL,:), compress)
	  CALL ncio_write_vector (file_restart, 'ch4_annavg_agnpp_rice', 'patch', landpatch, &
	     annavg_agnpp_component(METHANE_COMP_RICE,:), compress)
	  CALL ncio_write_vector (file_restart, 'ch4_annavg_bgnpp_soil', 'patch', landpatch, &
	     annavg_bgnpp_component(METHANE_COMP_SOIL,:), compress)
	  CALL ncio_write_vector (file_restart, 'ch4_annavg_bgnpp_rice', 'patch', landpatch, &
	     annavg_bgnpp_component(METHANE_COMP_RICE,:), compress)
	  CALL ncio_write_vector (file_restart, 'ch4_annavg_somhr_soil', 'patch', landpatch, &
	     annavg_somhr_component(METHANE_COMP_SOIL,:), compress)
	  CALL ncio_write_vector (file_restart, 'ch4_annavg_somhr_rice', 'patch', landpatch, &
	     annavg_somhr_component(METHANE_COMP_RICE,:), compress)
	  CALL ncio_write_vector (file_restart, 'ch4_annavg_finrw_soil', 'patch', landpatch, &
	     annavg_finrw_component(METHANE_COMP_SOIL,:), compress)
	  CALL ncio_write_vector (file_restart, 'ch4_annavg_finrw_rice', 'patch', landpatch, &
	     annavg_finrw_component(METHANE_COMP_RICE,:), compress)
	  CALL ncio_write_vector (file_restart, 'ch4_tempavg_agnpp_soil', 'patch', landpatch, &
	     tempavg_agnpp_component(METHANE_COMP_SOIL,:), compress)
	  CALL ncio_write_vector (file_restart, 'ch4_tempavg_agnpp_rice', 'patch', landpatch, &
	     tempavg_agnpp_component(METHANE_COMP_RICE,:), compress)
	  CALL ncio_write_vector (file_restart, 'ch4_tempavg_bgnpp_soil', 'patch', landpatch, &
	     tempavg_bgnpp_component(METHANE_COMP_SOIL,:), compress)
	  CALL ncio_write_vector (file_restart, 'ch4_tempavg_bgnpp_rice', 'patch', landpatch, &
	     tempavg_bgnpp_component(METHANE_COMP_RICE,:), compress)
	  CALL ncio_write_vector (file_restart, 'ch4_annsum_counter_soil', 'patch', landpatch, &
	     annsum_counter_component(METHANE_COMP_SOIL,:), compress)
	  CALL ncio_write_vector (file_restart, 'ch4_annsum_counter_rice', 'patch', landpatch, &
	     annsum_counter_component(METHANE_COMP_RICE,:), compress)
	  CALL ncio_write_vector (file_restart, 'ch4_tempavg_somhr_soil', 'patch', landpatch, &
	     tempavg_somhr_component(METHANE_COMP_SOIL,:), compress)
	  CALL ncio_write_vector (file_restart, 'ch4_tempavg_somhr_rice', 'patch', landpatch, &
	     tempavg_somhr_component(METHANE_COMP_RICE,:), compress)
	  CALL ncio_write_vector (file_restart, 'ch4_tempavg_finrw_soil', 'patch', landpatch, &
	     tempavg_finrw_component(METHANE_COMP_SOIL,:), compress)
	  CALL ncio_write_vector (file_restart, 'ch4_tempavg_finrw_rice', 'patch', landpatch, &
	     tempavg_finrw_component(METHANE_COMP_RICE,:), compress)
      CALL ncio_write_vector (file_restart, 'ch4_fsat_bef',         'patch', landpatch, fsat_bef,         compress)
	      CALL ncio_write_vector (file_restart, 'ch4_finundated_lag',   'patch', landpatch, finundated_lag,   compress)
	      CALL ncio_write_vector (file_restart, 'ch4_fsat_bef_soil', 'patch', landpatch, &
	         fsat_bef_component(METHANE_COMP_SOIL,:), compress)
	      CALL ncio_write_vector (file_restart, 'ch4_fsat_bef_rice', 'patch', landpatch, &
	         fsat_bef_component(METHANE_COMP_RICE,:), compress)
	      CALL ncio_write_vector (file_restart, 'ch4_finundated_lag_soil', 'patch', landpatch, &
	         finundated_lag_component(METHANE_COMP_SOIL,:), compress)
	      CALL ncio_write_vector (file_restart, 'ch4_finundated_lag_rice', 'patch', landpatch, &
	         finundated_lag_component(METHANE_COMP_RICE,:), compress)
	      CALL ncio_write_vector (file_restart, 'ch4_rice_fraction_prev', 'patch', landpatch, rice_fraction_prev, compress)
	      CALL ncio_write_vector (file_restart, 'ch4_methane_dfsat_tot','patch', landpatch, methane_dfsat_tot, compress)
	      CALL ncio_write_vector (file_restart, 'ch4_f_h2osfc',         'patch', landpatch, f_h2osfc,         compress)
	      CALL ncio_write_vector (file_restart, 'ch4_f_inund_levee_patch',       &
	                              'patch', landpatch, f_inund_levee_patch,       compress)
	      CALL ncio_write_vector (file_restart, 'ch4_f_inund_flood_patch',       &
	                              'patch', landpatch, f_inund_flood_patch,       compress)
	      CALL ncio_write_vector (file_restart, 'ch4_f_inund_flood_depth_patch', &
	                              'patch', landpatch, f_inund_flood_depth_patch, compress)
	   END SUBROUTINE write_methane_restart


	   SUBROUTINE read_methane_restart (file_restart, strict_restart, restart_schema)
	      USE MOD_LandPatch,     only: landpatch
	      USE MOD_MT_Methane_Const, only: DEF_METHANE
	      USE MOD_Vars_TimeInvariants, only: patchtype
	      USE MOD_MT_NcIO, only: ncio_read_vector => ncio_read_vector_complete, &
         ncio_set_complete_require_present, ncio_vector_var_present, ncio_vector_group_presence
#ifdef USEMPI
	      USE MOD_SPMD_Task, only: p_is_worker, p_is_master, p_comm_glb, p_err, &
	         MPI_IN_PLACE, MPI_LOGICAL, MPI_LOR, MPI_INTEGER, MPI_SUM, CoLM_stop
#else
	      USE MOD_SPMD_Task, only: p_is_worker, p_is_master, CoLM_stop
#endif
		      character(len=*), intent(in) :: file_restart
	      logical, intent(in), optional :: strict_restart
	      integer, intent(in), optional :: restart_schema
	      integer, parameter :: n_component_restart_fields = 33
	      character(len=40), parameter :: component_restart_fields(n_component_restart_fields) = &
	         [character(len=40) :: &
	          'ch4_conc_o2_unsat_soil', 'ch4_conc_o2_unsat_rice', &
	          'ch4_conc_o2_sat_soil', 'ch4_conc_o2_sat_rice', &
	          'ch4_conc_ch4_unsat_soil', 'ch4_conc_ch4_unsat_rice', &
	          'ch4_conc_ch4_sat_soil', 'ch4_conc_ch4_sat_rice', &
	          'ch4_layer_sat_lag_soil', 'ch4_layer_sat_lag_rice', &
	          'ch4_annavg_agnpp_soil', 'ch4_annavg_agnpp_rice', &
	          'ch4_annavg_bgnpp_soil', 'ch4_annavg_bgnpp_rice', &
	          'ch4_annavg_somhr_soil', 'ch4_annavg_somhr_rice', &
	          'ch4_annavg_finrw_soil', 'ch4_annavg_finrw_rice', &
	          'ch4_tempavg_agnpp_soil', 'ch4_tempavg_agnpp_rice', &
	          'ch4_tempavg_bgnpp_soil', 'ch4_tempavg_bgnpp_rice', &
	          'ch4_annsum_counter_soil', 'ch4_annsum_counter_rice', &
	          'ch4_tempavg_somhr_soil', 'ch4_tempavg_somhr_rice', &
	          'ch4_tempavg_finrw_soil', 'ch4_tempavg_finrw_rice', &
	          'ch4_fsat_bef_soil', 'ch4_fsat_bef_rice', &
	          'ch4_finundated_lag_soil', 'ch4_finundated_lag_rice', &
	          'ch4_rice_fraction_prev']
	      logical :: strict_restart_active
	      logical :: lake_water_ch4_present, lake_water_o2_present, lake_soilc_present
	      logical :: lake_phase_fields_present(3)
	      character(len=40), parameter :: lake_phase_fields(3) = [character(len=40) :: &
	         'ch4_lake_frozen_ch4_stock', 'ch4_lake_frozen_o2_stock', &
	         'ch4_lake_liquid_fraction_prev']
	      logical :: component_fields_present(n_component_restart_fields)
	      logical :: component_state_present, component_state_required
	      integer :: ipatch, npatch, component, restart_schema_active, corrupt_prognostic_values
	      integer :: invalid_lake_inventory
	      real(r8) :: lake_component_total, lake_inventory_tolerance
	      ! Round-off floor for lake stocks: a conservative exchange residual sits
	      ! a few 1e-26 below zero when the true stock is zero; reject only values
	      ! materially more negative than this.
	      real(r8), parameter :: restart_neg_floor = 1.e-18_r8

		      IF (.not. allocated(conc_methane)) RETURN
	      strict_restart_active = .false.
	      IF (present(strict_restart)) strict_restart_active = strict_restart
	      restart_schema_active = 0
	      IF (present(restart_schema)) restart_schema_active = restart_schema
	      component_state_required = strict_restart_active .and. restart_schema_active >= 3
	      CALL ncio_vector_group_presence(file_restart, component_restart_fields, &
	         landpatch, component_fields_present)
	      component_state_present = all(component_fields_present)
	      IF ((component_state_required .and. .not. component_state_present) .or. &
	          (any(component_fields_present) .and. .not. component_state_present)) THEN
	         IF (p_is_master) WRITE(*,'(A)') &
	            'ERROR: methane soil/rice component restart group is incomplete.'
	         CALL CoLM_stop()
	      ENDIF
	      ! Schema 1 predates these two stocks.  Presence probes distinguish a
	      ! legitimate missing-field sentinel from an explicitly stored corrupt
	      ! value while retaining mixed-block detection.  Schema 2's required
	      ! complete reads already perform the same presence/shape preflight.
	      IF (strict_restart_active .and. restart_schema_active >= 2) THEN
	         lake_water_ch4_present = .true.
	         lake_water_o2_present = .true.
	      ELSE
	         lake_water_ch4_present = ncio_vector_var_present(file_restart, &
	            'ch4_lake_water_ch4_stock', landpatch)
	         lake_water_o2_present = ncio_vector_var_present(file_restart, &
	            'ch4_lake_water_o2_stock', landpatch)
	      ENDIF
	      ! Committed schemas require lake_soilc.  Markerless legacy files may
	      ! predate it; in that case preserve the surface-data initialization
	      ! performed before restart loading instead of overwriting it with zero.
	      lake_soilc_present = strict_restart_active
	      IF (.not. strict_restart_active) lake_soilc_present = &
	         ncio_vector_var_present(file_restart, 'ch4_lake_soilc', landpatch)
	      IF (strict_restart_active .and. restart_schema_active >= 4) THEN
	         lake_phase_fields_present = .true.
	      ELSE
	         CALL ncio_vector_group_presence(file_restart, lake_phase_fields, &
	            landpatch, lake_phase_fields_present)
	         IF (any(lake_phase_fields_present) .and. .not. all(lake_phase_fields_present)) THEN
	            IF (p_is_master) WRITE(*,'(A)') &
	               'ERROR: methane lake frozen-phase restart group is incomplete.'
	            CALL CoLM_stop()
	         ENDIF
	      ENDIF
	      CALL ncio_read_vector (file_restart, 'ch4_conc_o2',          nl_soil, landpatch, conc_o2,          defval = 1._r8)
	      CALL ncio_read_vector (file_restart, 'ch4_conc_methane',     nl_soil, landpatch, conc_methane,     defval = 1.e-6_r8)
	      CALL ncio_read_vector (file_restart, 'ch4_totcol_methane',   landpatch, totcol_methane,            defval = spval)
		      CALL ncio_read_vector (file_restart, 'ch4_grnd_methane_cond',landpatch, grnd_methane_cond, &
		         defval = DEF_METHANE%grnd_methane_cond_default)
	      CALL ncio_read_vector (file_restart, 'ch4_conc_o2_unsat',    nl_soil, landpatch, conc_o2_unsat,    defval = 1._r8)
	      CALL ncio_read_vector (file_restart, 'ch4_conc_o2_sat',      nl_soil, landpatch, conc_o2_sat,      defval = 1._r8)
	      CALL ncio_read_vector (file_restart, 'ch4_conc_ch4_unsat',   nl_soil, landpatch, conc_methane_unsat, defval = 1.e-6_r8)
	      CALL ncio_read_vector (file_restart, 'ch4_conc_ch4_sat',     nl_soil, landpatch, conc_methane_sat,   defval = 1.e-6_r8)
	      IF (component_state_present) THEN
	         CALL read_component_2d('ch4_conc_o2_unsat_soil', conc_o2_unsat_component(:,METHANE_COMP_SOIL,:))
	         CALL read_component_2d('ch4_conc_o2_unsat_rice', conc_o2_unsat_component(:,METHANE_COMP_RICE,:))
	         CALL read_component_2d('ch4_conc_o2_sat_soil', conc_o2_sat_component(:,METHANE_COMP_SOIL,:))
	         CALL read_component_2d('ch4_conc_o2_sat_rice', conc_o2_sat_component(:,METHANE_COMP_RICE,:))
	         CALL read_component_2d('ch4_conc_ch4_unsat_soil', conc_methane_unsat_component(:,METHANE_COMP_SOIL,:))
	         CALL read_component_2d('ch4_conc_ch4_unsat_rice', conc_methane_unsat_component(:,METHANE_COMP_RICE,:))
	         CALL read_component_2d('ch4_conc_ch4_sat_soil', conc_methane_sat_component(:,METHANE_COMP_SOIL,:))
	         CALL read_component_2d('ch4_conc_ch4_sat_rice', conc_methane_sat_component(:,METHANE_COMP_RICE,:))
	      ELSE
	         conc_o2_unsat_component = spval
	         conc_o2_sat_component = spval
	         conc_methane_unsat_component = spval
	         conc_methane_sat_component = spval
	      ENDIF
	      CALL ncio_read_vector (file_restart, 'ch4_totcol_methane_unsat', landpatch, totcol_methane_unsat,    defval = spval)
	      CALL ncio_read_vector (file_restart, 'ch4_totcol_methane_sat',   landpatch, totcol_methane_sat,      defval = spval)
		      CALL ncio_read_vector (file_restart, 'ch4_grnd_methane_cond_unsat', landpatch, grnd_methane_cond_unsat, &
		         defval = DEF_METHANE%grnd_methane_cond_default)
		      CALL ncio_read_vector (file_restart, 'ch4_grnd_methane_cond_sat',   landpatch, grnd_methane_cond_sat, &
		         defval = DEF_METHANE%grnd_methane_cond_default)
	      CALL ncio_read_vector (file_restart, 'ch4_conc_o2_lake',     nl_soil, landpatch, conc_o2_lake,       defval = 1._r8)
	      CALL ncio_read_vector (file_restart, 'ch4_conc_ch4_lake',    nl_soil, landpatch, conc_methane_lake,  defval = 0._r8)
	      CALL ncio_read_vector (file_restart, 'ch4_totcol_lake',      landpatch, totcol_methane_lake,         defval = spval)
		      CALL ncio_read_vector (file_restart, 'ch4_grnd_methane_cond_lake', landpatch, grnd_methane_cond_lake, &
		         defval = DEF_METHANE%grnd_methane_cond_default)
	      IF (strict_restart_active .and. restart_schema_active == 1) THEN
	         CALL ncio_set_complete_require_present (.false.)
	      ENDIF
	      CALL ncio_read_vector (file_restart, 'ch4_lake_water_ch4_stock', landpatch, lake_water_ch4_stock, defval = spval)
	      CALL ncio_read_vector (file_restart, 'ch4_lake_water_o2_stock',  landpatch, lake_water_o2_stock,  defval = spval)
	      IF (strict_restart_active .and. restart_schema_active == 1) THEN
	         CALL ncio_set_complete_require_present (.true.)
	      ENDIF
	      IF (strict_restart_active .and. restart_schema_active < 4) &
	         CALL ncio_set_complete_require_present (.false.)
	      CALL ncio_read_vector (file_restart, 'ch4_lake_frozen_ch4_stock', &
	         landpatch, lake_frozen_ch4_stock, defval = 0._r8)
	      CALL ncio_read_vector (file_restart, 'ch4_lake_frozen_o2_stock', &
	         landpatch, lake_frozen_o2_stock, defval = 0._r8)
	      CALL ncio_read_vector (file_restart, 'ch4_lake_liquid_fraction_prev', &
	         landpatch, lake_liquid_fraction_prev, defval = spval)
	      IF (strict_restart_active .and. restart_schema_active < 4) &
	         CALL ncio_set_complete_require_present (.true.)
	      CALL ncio_read_vector (file_restart, 'ch4_layer_sat_lag',    nl_soil, landpatch, layer_sat_lag,    defval = spval)
	      IF (component_state_present) THEN
	         CALL read_component_2d('ch4_layer_sat_lag_soil', layer_sat_lag_component(:,METHANE_COMP_SOIL,:))
	         CALL read_component_2d('ch4_layer_sat_lag_rice', layer_sat_lag_component(:,METHANE_COMP_RICE,:))
	      ELSE
	         layer_sat_lag_component = spval
	      ENDIF
      IF (lake_soilc_present) CALL ncio_read_vector (file_restart, 'ch4_lake_soilc', &
         nl_soil, landpatch, lake_soilc, defval = 0._r8)
      CALL ncio_read_vector (file_restart, 'ch4_annavg_agnpp',     landpatch, annavg_agnpp,     defval = 0._r8)
      CALL ncio_read_vector (file_restart, 'ch4_annavg_bgnpp',     landpatch, annavg_bgnpp,     defval = 0._r8)
      CALL ncio_read_vector (file_restart, 'ch4_annavg_somhr',     landpatch, annavg_somhr,     defval = 0._r8)
      CALL ncio_read_vector (file_restart, 'ch4_annavg_finrw',     landpatch, annavg_finrw,     defval = spval)
      CALL ncio_read_vector (file_restart, 'ch4_tempavg_agnpp',    landpatch, tempavg_agnpp,    defval = 0._r8)
      CALL ncio_read_vector (file_restart, 'ch4_tempavg_bgnpp',    landpatch, tempavg_bgnpp,    defval = 0._r8)
      CALL ncio_read_vector (file_restart, 'ch4_annsum_counter',   landpatch, annsum_counter,   defval = 0._r8)
      CALL ncio_read_vector (file_restart, 'ch4_tempavg_somhr',    landpatch, tempavg_somhr,    defval = 0._r8)
      CALL ncio_read_vector (file_restart, 'ch4_tempavg_finrw',    landpatch, tempavg_finrw,    defval = 0._r8)
	  IF (component_state_present) THEN
	  CALL read_component_1d('ch4_annavg_agnpp_soil', annavg_agnpp_component(METHANE_COMP_SOIL,:))
	  CALL read_component_1d('ch4_annavg_agnpp_rice', annavg_agnpp_component(METHANE_COMP_RICE,:))
	  CALL read_component_1d('ch4_annavg_bgnpp_soil', annavg_bgnpp_component(METHANE_COMP_SOIL,:))
	  CALL read_component_1d('ch4_annavg_bgnpp_rice', annavg_bgnpp_component(METHANE_COMP_RICE,:))
	  CALL read_component_1d('ch4_annavg_somhr_soil', annavg_somhr_component(METHANE_COMP_SOIL,:))
	  CALL read_component_1d('ch4_annavg_somhr_rice', annavg_somhr_component(METHANE_COMP_RICE,:))
	  CALL read_component_1d('ch4_annavg_finrw_soil', annavg_finrw_component(METHANE_COMP_SOIL,:))
	  CALL read_component_1d('ch4_annavg_finrw_rice', annavg_finrw_component(METHANE_COMP_RICE,:))
	  CALL read_component_1d('ch4_tempavg_agnpp_soil', tempavg_agnpp_component(METHANE_COMP_SOIL,:))
	  CALL read_component_1d('ch4_tempavg_agnpp_rice', tempavg_agnpp_component(METHANE_COMP_RICE,:))
	  CALL read_component_1d('ch4_tempavg_bgnpp_soil', tempavg_bgnpp_component(METHANE_COMP_SOIL,:))
	  CALL read_component_1d('ch4_tempavg_bgnpp_rice', tempavg_bgnpp_component(METHANE_COMP_RICE,:))
	  CALL read_component_1d('ch4_annsum_counter_soil', annsum_counter_component(METHANE_COMP_SOIL,:))
	  CALL read_component_1d('ch4_annsum_counter_rice', annsum_counter_component(METHANE_COMP_RICE,:))
	  CALL read_component_1d('ch4_tempavg_somhr_soil', tempavg_somhr_component(METHANE_COMP_SOIL,:))
	  CALL read_component_1d('ch4_tempavg_somhr_rice', tempavg_somhr_component(METHANE_COMP_RICE,:))
	  CALL read_component_1d('ch4_tempavg_finrw_soil', tempavg_finrw_component(METHANE_COMP_SOIL,:))
	  CALL read_component_1d('ch4_tempavg_finrw_rice', tempavg_finrw_component(METHANE_COMP_RICE,:))
	  ELSE
	     annavg_agnpp_component = spval
	     annavg_bgnpp_component = spval
	     annavg_somhr_component = spval
	     annavg_finrw_component = spval
	     tempavg_agnpp_component = spval
	     tempavg_bgnpp_component = spval
	     annsum_counter_component = spval
	     tempavg_somhr_component = spval
	     tempavg_finrw_component = spval
	  ENDIF
	      CALL ncio_read_vector (file_restart, 'ch4_fsat_bef',         landpatch, fsat_bef,         defval = spval)
	      CALL ncio_read_vector (file_restart, 'ch4_finundated_lag',   landpatch, finundated_lag,   defval = spval)
	      IF (component_state_present) THEN
	         CALL read_component_1d('ch4_fsat_bef_soil', fsat_bef_component(METHANE_COMP_SOIL,:))
	         CALL read_component_1d('ch4_fsat_bef_rice', fsat_bef_component(METHANE_COMP_RICE,:))
	         CALL read_component_1d('ch4_finundated_lag_soil', finundated_lag_component(METHANE_COMP_SOIL,:))
	         CALL read_component_1d('ch4_finundated_lag_rice', finundated_lag_component(METHANE_COMP_RICE,:))
	         CALL ncio_read_vector (file_restart, 'ch4_rice_fraction_prev', landpatch, rice_fraction_prev, defval = spval)
	      ELSE
	         fsat_bef_component = spval
	         finundated_lag_component = spval
	         rice_fraction_prev = spval
	      ENDIF
	      CALL ncio_read_vector (file_restart, 'ch4_methane_dfsat_tot',landpatch, methane_dfsat_tot, defval = 0._r8)
	      CALL ncio_read_vector (file_restart, 'ch4_f_h2osfc',         landpatch, f_h2osfc,         defval = 0._r8)
	      CALL ncio_read_vector (file_restart, 'ch4_f_inund_levee_patch',       landpatch, &
	                             f_inund_levee_patch,       defval = 0._r8)
	      CALL ncio_read_vector (file_restart, 'ch4_f_inund_flood_patch',       landpatch, &
	                             f_inund_flood_patch,       defval = 0._r8)
	      CALL ncio_read_vector (file_restart, 'ch4_f_inund_flood_depth_patch', landpatch, &
	                             f_inund_flood_depth_patch, defval = 0._r8)

	      ! A committed transaction is reproducible only if its prognostic state
	      ! is used verbatim.  Reject corrupt/negative concentrations, inventories,
	      ! and lake substrate before any repair.  Schema 1 legitimately lacks the
	      ! two water-column stocks; those missing values are the sole strict-mode
	      ! migration exception.  Marker-less legacy files retain the documented
	      ! sanitation path, with a collective warning rather than silent mutation.
	      corrupt_prognostic_values = 0
	      IF (p_is_worker) THEN
	         corrupt_prognostic_values = &
	            count(invalid_restart_value(conc_o2) .or. conc_o2 < 0._r8) + &
	            count(invalid_restart_value(conc_o2_unsat) .or. conc_o2_unsat < 0._r8) + &
	            count(invalid_restart_value(conc_o2_sat) .or. conc_o2_sat < 0._r8) + &
	            count(invalid_restart_value(conc_o2_lake) .or. conc_o2_lake < 0._r8) + &
	            count(invalid_restart_value(conc_methane) .or. conc_methane < 0._r8) + &
	            count(invalid_restart_value(conc_methane_unsat) .or. conc_methane_unsat < 0._r8) + &
	            count(invalid_restart_value(conc_methane_sat) .or. conc_methane_sat < 0._r8) + &
	            count(invalid_restart_value(conc_methane_lake) .or. conc_methane_lake < 0._r8) + &
	            count(invalid_restart_value(totcol_methane) .or. totcol_methane < 0._r8) + &
	            count(invalid_restart_value(totcol_methane_unsat) .or. totcol_methane_unsat < 0._r8) + &
	            count(invalid_restart_value(totcol_methane_sat) .or. totcol_methane_sat < 0._r8) + &
	            count(invalid_restart_value(totcol_methane_lake) .or. totcol_methane_lake < 0._r8) + &
	            count(invalid_restart_value(lake_soilc) .or. lake_soilc < 0._r8) + &
	            ! Use < 0, not <= 0: a ground conductance of exactly zero is a
	            ! valid physical state, not corruption.  A component's conductance
	            ! is zero wherever it does not apply - the _lake term on non-lake
	            ! patches, the _sat term on unsaturated patches - and the paddy
	            ! aggregation writes it directly (Driver: wss*grnd_cond_sat +
	            ! wsr*grnd_cond_sat).  The sibling sanitize path below already
	            ! treats <= 0 as recoverable (it resets to the default rather than
	            ! aborting), so the strict path must not kill the run for the same
	            ! value; and grnd_methane_cond is recomputed every step, so the
	            ! restart value never enters the solve regardless.
	            ! NaN/spval and genuinely negative values are still rejected.
	            count(invalid_restart_value(grnd_methane_cond) .or. grnd_methane_cond < 0._r8) + &
	            count(invalid_restart_value(grnd_methane_cond_unsat) .or. grnd_methane_cond_unsat < 0._r8) + &
	            count(invalid_restart_value(grnd_methane_cond_sat) .or. grnd_methane_cond_sat < 0._r8) + &
	            count(invalid_restart_value(grnd_methane_cond_lake) .or. grnd_methane_cond_lake < 0._r8) + &
	            count(invalid_restart_value(annavg_agnpp) .or. annavg_agnpp < 0._r8) + &
	            count(invalid_restart_value(annavg_bgnpp) .or. annavg_bgnpp < 0._r8) + &
	            count(invalid_restart_value(annavg_somhr) .or. annavg_somhr < 0._r8) + &
	            count(invalid_restart_fraction_or_sentinel(annavg_finrw)) + &
	            count(invalid_restart_value(tempavg_agnpp) .or. tempavg_agnpp < 0._r8) + &
	            count(invalid_restart_value(tempavg_bgnpp) .or. tempavg_bgnpp < 0._r8) + &
	            count(invalid_restart_value(tempavg_somhr) .or. tempavg_somhr < 0._r8) + &
	            count(invalid_restart_value(tempavg_finrw) .or. tempavg_finrw < 0._r8) + &
	            count(invalid_restart_value(annsum_counter) .or. annsum_counter < 0._r8) + &
	            count(invalid_restart_fraction_or_sentinel(fsat_bef)) + &
	            count(invalid_restart_fraction_or_sentinel(finundated_lag)) + &
	            count(invalid_restart_fraction_or_sentinel(layer_sat_lag)) + &
	            count(invalid_restart_value(methane_dfsat_tot)) + &
	            count(invalid_restart_value(f_h2osfc) .or. f_h2osfc < 0._r8 .or. f_h2osfc > 1._r8) + &
	            count(invalid_restart_value(f_inund_levee_patch) .or. &
	                  f_inund_levee_patch < 0._r8 .or. f_inund_levee_patch > 1._r8) + &
	            count(invalid_restart_value(f_inund_flood_patch) .or. &
	                  f_inund_flood_patch < 0._r8 .or. f_inund_flood_patch > 1._r8) + &
	            count(invalid_restart_value(f_inund_flood_depth_patch) .or. &
	                  f_inund_flood_depth_patch < 0._r8)
	         ! Allow a round-off floor on the lake stocks.  These stocks are the
	         ! residual of conservative sediment/water/air exchanges, so a stock
	         ! that is physically zero lands a few 1e-26 below zero after the
	         ! solve.  Only a materially negative value signals corruption; the
	         ! sibling sanitize path already resets < 0 to zero rather than
	         ! aborting, so the strict path must not kill the run for round-off.
	         IF (lake_water_ch4_present) THEN
	            corrupt_prognostic_values = corrupt_prognostic_values + &
	               count(invalid_restart_value(lake_water_ch4_stock) .or. lake_water_ch4_stock < -restart_neg_floor)
	         ENDIF
	         IF (lake_water_o2_present) THEN
	            corrupt_prognostic_values = corrupt_prognostic_values + &
	               count(invalid_restart_value(lake_water_o2_stock) .or. lake_water_o2_stock < -restart_neg_floor)
	         ENDIF
	         IF (all(lake_phase_fields_present)) THEN
	            corrupt_prognostic_values = corrupt_prognostic_values + &
	               count(invalid_restart_value(lake_frozen_ch4_stock) .or. lake_frozen_ch4_stock < -restart_neg_floor) + &
	               count(invalid_restart_value(lake_frozen_o2_stock) .or. lake_frozen_o2_stock < -restart_neg_floor) + &
	               count(invalid_restart_fraction_or_sentinel(lake_liquid_fraction_prev))
	         ENDIF
	         IF (component_state_present) THEN
	            corrupt_prognostic_values = corrupt_prognostic_values + &
	               count(invalid_restart_value(conc_o2_unsat_component) .or. conc_o2_unsat_component < 0._r8) + &
	               count(invalid_restart_value(conc_o2_sat_component) .or. conc_o2_sat_component < 0._r8) + &
	               count(invalid_restart_value(conc_methane_unsat_component) .or. &
	                     conc_methane_unsat_component < 0._r8) + &
	               count(invalid_restart_value(conc_methane_sat_component) .or. &
	                     conc_methane_sat_component < 0._r8) + &
	               count(invalid_restart_fraction_or_sentinel(layer_sat_lag_component)) + &
	               count(invalid_restart_value(annavg_agnpp_component) .or. annavg_agnpp_component < 0._r8) + &
	               count(invalid_restart_value(annavg_bgnpp_component) .or. annavg_bgnpp_component < 0._r8) + &
	               count(invalid_restart_value(annavg_somhr_component) .or. annavg_somhr_component < 0._r8) + &
	               count(invalid_restart_fraction_or_sentinel(annavg_finrw_component)) + &
	               count(invalid_restart_value(tempavg_agnpp_component) .or. tempavg_agnpp_component < 0._r8) + &
	               count(invalid_restart_value(tempavg_bgnpp_component) .or. tempavg_bgnpp_component < 0._r8) + &
	               count(invalid_restart_value(annsum_counter_component) .or. annsum_counter_component < 0._r8) + &
	               count(invalid_restart_value(tempavg_somhr_component) .or. tempavg_somhr_component < 0._r8) + &
	               count(invalid_restart_value(tempavg_finrw_component) .or. tempavg_finrw_component < 0._r8) + &
	               count(invalid_restart_fraction_or_sentinel(fsat_bef_component)) + &
	               count(invalid_restart_fraction_or_sentinel(finundated_lag_component)) + &
	               count(invalid_restart_value(rice_fraction_prev) .or. rice_fraction_prev < 0._r8 .or. &
	                     rice_fraction_prev > 1._r8)
	         ENDIF
	      ENDIF
#ifdef USEMPI
	      CALL mpi_allreduce(MPI_IN_PLACE, corrupt_prognostic_values, 1, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif
	      IF (corrupt_prognostic_values > 0) THEN
	         IF (strict_restart_active) THEN
	            IF (p_is_master) WRITE(*,'(A,I0,A)') &
	               'ERROR: committed methane restart contains ', corrupt_prognostic_values, &
	               ' invalid or negative prognostic state values; refusing checkpoint.'
	            CALL CoLM_stop()
	         ELSEIF (p_is_master) THEN
	            WRITE(*,'(A,I0,A)') 'WARNING: legacy methane restart contains ', &
	               corrupt_prognostic_values, ' invalid or negative prognostic values; sanitizing once.'
	         ENDIF
	      ENDIF

	      ! Sanitize aggregate concentration sentinels before they are used as
	      ! the schema-1 fallback for either independent process column.
	      WHERE (invalid_restart_value(conc_o2))             conc_o2             = 1._r8
	      WHERE (invalid_restart_value(conc_o2_unsat))       conc_o2_unsat       = 1._r8
	      WHERE (invalid_restart_value(conc_o2_sat))         conc_o2_sat         = 1._r8
	      WHERE (invalid_restart_value(conc_methane))        conc_methane        = 1.e-6_r8
	      WHERE (invalid_restart_value(conc_methane_unsat))  conc_methane_unsat  = 1.e-6_r8
	      WHERE (invalid_restart_value(conc_methane_sat))    conc_methane_sat    = 1.e-6_r8

	      ! Schema-1 restarts contain only the aggregate methane column.  Seed
	      ! both independent process columns from that state; new restarts keep
	      ! their separate histories.  The aggregate fields remain the sole
	      ! restart clip-credit representation, so fallback copies are not
	      ! counted a second time.
	      DO component = 1, N_METHANE_COMP
	         WHERE (invalid_restart_value(conc_o2_unsat_component(:,component,:)))
	            conc_o2_unsat_component(:,component,:) = conc_o2_unsat
	         END WHERE
	         WHERE (invalid_restart_value(conc_o2_sat_component(:,component,:)))
	            conc_o2_sat_component(:,component,:) = conc_o2_sat
	         END WHERE
	         WHERE (invalid_restart_value(conc_methane_unsat_component(:,component,:)))
	            conc_methane_unsat_component(:,component,:) = conc_methane_unsat
	         END WHERE
	         WHERE (invalid_restart_value(conc_methane_sat_component(:,component,:)))
	            conc_methane_sat_component(:,component,:) = conc_methane_sat
	         END WHERE
	         IF (component_state_present) THEN
	            WHERE (invalid_restart_fraction_or_sentinel(layer_sat_lag_component(:,component,:)))
	               layer_sat_lag_component(:,component,:) = layer_sat_lag
	            END WHERE
	         ELSE
	            WHERE (invalid_restart_value(layer_sat_lag_component(:,component,:)))
	               layer_sat_lag_component(:,component,:) = layer_sat_lag
	            END WHERE
	         ENDIF
	         WHERE (invalid_restart_value(annavg_agnpp_component(component,:)))
	            annavg_agnpp_component(component,:) = annavg_agnpp
	         END WHERE
	         WHERE (invalid_restart_value(annavg_bgnpp_component(component,:)))
	            annavg_bgnpp_component(component,:) = annavg_bgnpp
	         END WHERE
	         WHERE (invalid_restart_value(annavg_somhr_component(component,:)))
	            annavg_somhr_component(component,:) = annavg_somhr
	         END WHERE
	         IF (component_state_present) THEN
	            WHERE (invalid_restart_fraction_or_sentinel(annavg_finrw_component(component,:)))
	               annavg_finrw_component(component,:) = annavg_finrw
	            END WHERE
	         ELSE
	            WHERE (invalid_restart_value(annavg_finrw_component(component,:)))
	               annavg_finrw_component(component,:) = annavg_finrw
	            END WHERE
	         ENDIF
	         WHERE (invalid_restart_value(tempavg_agnpp_component(component,:)))
	            tempavg_agnpp_component(component,:) = tempavg_agnpp
	         END WHERE
	         WHERE (invalid_restart_value(tempavg_bgnpp_component(component,:)))
	            tempavg_bgnpp_component(component,:) = tempavg_bgnpp
	         END WHERE
	         WHERE (invalid_restart_value(annsum_counter_component(component,:)))
	            annsum_counter_component(component,:) = annsum_counter
	         END WHERE
	         WHERE (invalid_restart_value(tempavg_somhr_component(component,:)))
	            tempavg_somhr_component(component,:) = tempavg_somhr
	         END WHERE
	         WHERE (invalid_restart_value(tempavg_finrw_component(component,:)))
	            tempavg_finrw_component(component,:) = tempavg_finrw
	         END WHERE
	         IF (component_state_present) THEN
	            WHERE (invalid_restart_fraction_or_sentinel(fsat_bef_component(component,:)))
	               fsat_bef_component(component,:) = fsat_bef
	            END WHERE
	         ELSE
	            WHERE (invalid_restart_value(fsat_bef_component(component,:)))
	               fsat_bef_component(component,:) = fsat_bef
	            END WHERE
	         ENDIF
	         IF (component_state_present) THEN
	            WHERE (invalid_restart_fraction_or_sentinel(finundated_lag_component(component,:)))
	               finundated_lag_component(component,:) = finundated_lag
	            END WHERE
	         ELSE
	            WHERE (invalid_restart_value(finundated_lag_component(component,:)))
	               finundated_lag_component(component,:) = finundated_lag
	            END WHERE
	         ENDIF
	      END DO
	      WHERE (invalid_restart_value(conc_o2_unsat_component)) conc_o2_unsat_component = 1._r8
	      WHERE (invalid_restart_value(conc_o2_sat_component)) conc_o2_sat_component = 1._r8
	      WHERE (invalid_restart_value(conc_methane_unsat_component)) &
	         conc_methane_unsat_component = 1.e-6_r8
	      WHERE (invalid_restart_value(conc_methane_sat_component)) &
	         conc_methane_sat_component = 1.e-6_r8
	      WHERE (invalid_restart_value(rice_fraction_prev)) rice_fraction_prev = 0._r8
	      rice_fraction_prev = min(max(rice_fraction_prev, 0._r8), 1._r8)

	      WHERE (invalid_restart_value(conc_o2))           conc_o2           = 1._r8
      WHERE (invalid_restart_value(conc_o2_unsat))     conc_o2_unsat     = 1._r8
      WHERE (invalid_restart_value(conc_o2_sat))       conc_o2_sat       = 1._r8
      WHERE (invalid_restart_value(conc_o2_lake))      conc_o2_lake      = 1._r8
      WHERE (invalid_restart_value(conc_methane))      conc_methane      = 1.e-6_r8
      WHERE (invalid_restart_value(conc_methane_unsat)) conc_methane_unsat = 1.e-6_r8
      WHERE (invalid_restart_value(conc_methane_sat))  conc_methane_sat  = 1.e-6_r8
      WHERE (invalid_restart_value(conc_methane_lake)) conc_methane_lake = 0._r8
	      WHERE (invalid_restart_value(grnd_methane_cond) .or. grnd_methane_cond <= 0._r8) &
	         grnd_methane_cond = DEF_METHANE%grnd_methane_cond_default
	      WHERE (invalid_restart_value(grnd_methane_cond_unsat) .or. grnd_methane_cond_unsat <= 0._r8) &
	         grnd_methane_cond_unsat = DEF_METHANE%grnd_methane_cond_default
	      WHERE (invalid_restart_value(grnd_methane_cond_sat) .or. grnd_methane_cond_sat <= 0._r8) &
	         grnd_methane_cond_sat = DEF_METHANE%grnd_methane_cond_default
		      WHERE (invalid_restart_value(grnd_methane_cond_lake) .or. grnd_methane_cond_lake <= 0._r8) &
		         grnd_methane_cond_lake = DEF_METHANE%grnd_methane_cond_default
	      WHERE (invalid_restart_value(f_inund_levee_patch) .or. f_inund_levee_patch < 0._r8) &
	         f_inund_levee_patch = 0._r8
	      WHERE (invalid_restart_value(f_inund_flood_patch) .or. f_inund_flood_patch < 0._r8) &
	         f_inund_flood_patch = 0._r8
	      WHERE (invalid_restart_value(f_inund_flood_depth_patch) .or. f_inund_flood_depth_patch < 0._r8) &
	         f_inund_flood_depth_patch = 0._r8
	      WHERE (f_inund_levee_patch > 1._r8) f_inund_levee_patch = 1._r8
	      WHERE (f_inund_flood_patch > 1._r8) f_inund_flood_patch = 1._r8
	      WHERE (invalid_restart_value(annavg_agnpp) .or. annavg_agnpp < 0._r8) annavg_agnpp = 0._r8
	      WHERE (invalid_restart_value(annavg_bgnpp) .or. annavg_bgnpp < 0._r8) annavg_bgnpp = 0._r8
	      WHERE (invalid_restart_value(annavg_somhr) .or. annavg_somhr < 0._r8) annavg_somhr = 0._r8
	      WHERE (invalid_restart_fraction_or_sentinel(annavg_finrw)) annavg_finrw = spval
	      WHERE (invalid_restart_value(tempavg_agnpp) .or. tempavg_agnpp < 0._r8) tempavg_agnpp = 0._r8
	      WHERE (invalid_restart_value(tempavg_bgnpp) .or. tempavg_bgnpp < 0._r8) tempavg_bgnpp = 0._r8
	      WHERE (invalid_restart_value(tempavg_somhr) .or. tempavg_somhr < 0._r8) tempavg_somhr = 0._r8
	      WHERE (invalid_restart_value(tempavg_finrw) .or. tempavg_finrw < 0._r8) tempavg_finrw = 0._r8
	      WHERE (invalid_restart_value(annsum_counter) .or. annsum_counter < 0._r8) annsum_counter = 0._r8
	      WHERE (invalid_restart_value(methane_dfsat_tot)) methane_dfsat_tot = 0._r8
	      WHERE (invalid_restart_value(f_h2osfc) .or. f_h2osfc < 0._r8) f_h2osfc = 0._r8
	      WHERE (f_h2osfc > 1._r8) f_h2osfc = 1._r8
	      ! These lag states feed differences and exponential memory updates before
	      ! any arithmetic can safely repair a NaN.  Preserve the normal Physics
	      ! cold-start path by converting corrupt restart values to its sentinel.
	      WHERE (invalid_restart_fraction_or_sentinel(fsat_bef)) fsat_bef = spval
	      WHERE (invalid_restart_fraction_or_sentinel(finundated_lag)) finundated_lag = spval
	      WHERE (invalid_restart_fraction_or_sentinel(layer_sat_lag)) layer_sat_lag = spval

	      ! invalid_restart_value catches NaN/spval but not negatives.  A stale
      ! restart with sub-zero CH4 or O2 would feed phase-partition and
      ! Michaelis-Menten kinetics, yielding negative oxidation / production.
	      ! Clip concentrations defensively.  They are diagnostic/split views of
	      ! the same physical inventory, so do not count them again in the restart
	      ! mass impulse below.
      WHERE (conc_o2           < 0._r8) conc_o2           = 0._r8
      WHERE (conc_o2_unsat     < 0._r8) conc_o2_unsat     = 0._r8
      WHERE (conc_o2_sat       < 0._r8) conc_o2_sat       = 0._r8
      WHERE (conc_o2_lake      < 0._r8) conc_o2_lake      = 0._r8
      WHERE (conc_methane      < 0._r8) conc_methane      = 0._r8
      WHERE (conc_methane_unsat < 0._r8) conc_methane_unsat = 0._r8
      WHERE (conc_methane_sat  < 0._r8) conc_methane_sat  = 0._r8
      WHERE (conc_methane_lake < 0._r8) conc_methane_lake = 0._r8
	  WHERE (conc_o2_unsat_component < 0._r8) conc_o2_unsat_component = 0._r8
	  WHERE (conc_o2_sat_component < 0._r8) conc_o2_sat_component = 0._r8
	  WHERE (conc_methane_unsat_component < 0._r8) conc_methane_unsat_component = 0._r8
	  WHERE (conc_methane_sat_component < 0._r8) conc_methane_sat_component = 0._r8

      ! Only force the cold-start sentinel when the *combined* column
      ! inventory (totcol_methane) is missing — that is the field the
      ! first-step balance check actually references.  A restart that
      ! lacks one of the sat/unsat/lake split inventories (e.g. an older
      ! file written before those were persisted) should NOT discard a
      ! valid fsat_bef; the per-variable defval at lines 1000-1003 will
      ! zero the missing split fields and physics will reconcile them on
      ! the next time step.
      WHERE (invalid_restart_value(totcol_methane))
         fsat_bef = spval
      END WHERE
	      ! Missing water-column fields identify an old restart.  Force one
	      ! cold-start step so Physics can initialize the dissolved inventories
	      ! at atmospheric Henry equilibrium without reporting that initialization
	      ! as a physical methane-budget residual.
	      IF (allocated(patchtype)) THEN
	         IF (size(patchtype) == size(fsat_bef)) THEN
	            WHERE (patchtype == PATCHTYPE_LAKE .and. &
	                   (invalid_restart_value(lake_water_ch4_stock) .or. &
	                    invalid_restart_value(lake_water_o2_stock)))
	               fsat_bef = spval
	            END WHERE
	         ENDIF
	      ENDIF
      WHERE (invalid_restart_value(totcol_methane))       totcol_methane       = 0._r8
      WHERE (invalid_restart_value(totcol_methane_unsat)) totcol_methane_unsat = 0._r8
      WHERE (invalid_restart_value(totcol_methane_sat))   totcol_methane_sat   = 0._r8
      WHERE (invalid_restart_value(totcol_methane_lake))  totcol_methane_lake  = 0._r8
	      WHERE (totcol_methane < 0._r8 .or. lake_water_ch4_stock < 0._r8)
	         fsat_bef = spval
	      END WHERE
      WHERE (totcol_methane       < 0._r8) totcol_methane       = 0._r8
      WHERE (totcol_methane_unsat < 0._r8) totcol_methane_unsat = 0._r8
      WHERE (totcol_methane_sat   < 0._r8) totcol_methane_sat   = 0._r8
      WHERE (totcol_methane_lake  < 0._r8) totcol_methane_lake  = 0._r8
	      WHERE (invalid_restart_value(lake_water_ch4_stock) .or. lake_water_ch4_stock < 0._r8) &
	         lake_water_ch4_stock = 0._r8
	      WHERE (invalid_restart_value(lake_water_o2_stock) .or. lake_water_o2_stock < 0._r8) &
	         lake_water_o2_stock = 0._r8
	      WHERE (invalid_restart_value(lake_frozen_ch4_stock) .or. lake_frozen_ch4_stock < 0._r8) &
	         lake_frozen_ch4_stock = 0._r8
	      WHERE (invalid_restart_value(lake_frozen_o2_stock) .or. lake_frozen_o2_stock < 0._r8) &
	         lake_frozen_o2_stock = 0._r8
	      WHERE (invalid_restart_value(lake_liquid_fraction_prev) .or. &
	             lake_liquid_fraction_prev < 0._r8 .or. lake_liquid_fraction_prev > 1._r8) &
	         lake_liquid_fraction_prev = spval

	      ! The combined lake inventory is redundant with sediment + water, but
	      ! it is the first-step budget reference.  A committed schema restart
	      ! must therefore be internally consistent; silently choosing either
	      ! side would create or destroy CH4.  Legacy files are repaired by
	      ! trusting the mechanistic component states and forcing a cold step.
	      invalid_lake_inventory = 0
	      IF (p_is_worker) THEN
	         IF (.not. allocated(patchtype)) THEN
	            invalid_lake_inventory = merge(1, 0, strict_restart_active)
	         ELSE
	            IF (strict_restart_active .and. &
	                (size(patchtype) /= size(totcol_methane) .or. &
	                 size(patchtype) /= size(totcol_methane_lake) .or. &
	                 size(patchtype) /= size(lake_water_ch4_stock) .or. &
	                 size(patchtype) /= size(lake_frozen_ch4_stock))) &
	               invalid_lake_inventory = 1
	            npatch = min(size(patchtype), size(totcol_methane), &
	               size(totcol_methane_lake), size(lake_water_ch4_stock), &
	               size(lake_frozen_ch4_stock))
	            DO ipatch = 1, npatch
	               IF (patchtype(ipatch) /= PATCHTYPE_LAKE) CYCLE
	               lake_component_total = totcol_methane_lake(ipatch) + &
	                  lake_water_ch4_stock(ipatch) + lake_frozen_ch4_stock(ipatch)
	               lake_inventory_tolerance = 1.e-12_r8 + 1.e-10_r8 * &
	                  max(abs(totcol_methane(ipatch)), abs(lake_component_total))
	               IF (abs(totcol_methane(ipatch) - lake_component_total) > lake_inventory_tolerance) THEN
	                  IF (strict_restart_active) THEN
	                     invalid_lake_inventory = 1
	                  ELSE
	                     totcol_methane(ipatch) = lake_component_total
	                     fsat_bef(ipatch) = spval
	                  ENDIF
	               ENDIF
	            ENDDO
	         ENDIF
	      ENDIF
#ifdef USEMPI
	      CALL mpi_allreduce(MPI_IN_PLACE, invalid_lake_inventory, 1, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif
	      IF (invalid_lake_inventory > 0) THEN
	         IF (p_is_master) WRITE(*,'(A)') &
	            'ERROR: methane restart lake inventory violates total = sediment + water; refusing corrupt checkpoint.'
	         CALL CoLM_stop()
	      ENDIF

      ! Clean spval/NaN/negative values from legacy files that contain this
      ! field.  A missing legacy field was left at its surface-data value above.
      WHERE (invalid_restart_value(lake_soilc) .or. lake_soilc < 0._r8) &
         lake_soilc = 0._r8

	   CONTAINS
	      SUBROUTINE read_component_2d(varname, target)
	         character(len=*), intent(in) :: varname
	         real(r8), intent(out) :: target(:,:)
	         real(r8), allocatable :: values(:,:)

	         allocate(values(nl_soil,size(target,2)))
	         CALL ncio_read_vector(file_restart, varname, nl_soil, landpatch, values, defval=spval)
	         target = values
	         deallocate(values)
	      END SUBROUTINE read_component_2d

	      SUBROUTINE read_component_1d(varname, target)
	         character(len=*), intent(in) :: varname
	         real(r8), intent(out) :: target(:)
	         real(r8), allocatable :: values(:)

	         allocate(values(size(target)))
	         CALL ncio_read_vector(file_restart, varname, landpatch, values, defval=spval)
	         target = values
	         deallocate(values)
	      END SUBROUTINE read_component_1d

	   END SUBROUTINE read_methane_restart

   SUBROUTINE save_methane_lulcc_state ()
      USE MOD_Vars_TimeInvariants, only: patchtype
	  integer :: ip, component, nsave

      IF (.not. allocated(conc_methane)) THEN
         methane_lulcc_snapshot_valid = .false.
         RETURN
      ENDIF

      CALL clear_methane_lulcc_snapshot ()

      allocate(lulcc_conc_o2_old(nl_soil,size(conc_o2,2))); lulcc_conc_o2_old = conc_o2
      allocate(lulcc_conc_methane_old(nl_soil,size(conc_methane,2))); lulcc_conc_methane_old = conc_methane
      allocate(lulcc_co2_decomp_depth_old(nl_soil,size(co2_decomp_depth,2))); lulcc_co2_decomp_depth_old = co2_decomp_depth
      allocate(lulcc_co2_oxid_depth_old(nl_soil,size(co2_oxid_depth,2))); lulcc_co2_oxid_depth_old = co2_oxid_depth
      allocate(lulcc_co2_aere_depth_old(nl_soil,size(co2_aere_depth,2))); lulcc_co2_aere_depth_old = co2_aere_depth
      allocate(lulcc_co2_decomp_depth_unsat_old(nl_soil,size(co2_decomp_depth_unsat,2))); lulcc_co2_decomp_depth_unsat_old = co2_decomp_depth_unsat
      allocate(lulcc_co2_decomp_depth_sat_old(nl_soil,size(co2_decomp_depth_sat,2))); lulcc_co2_decomp_depth_sat_old = co2_decomp_depth_sat
      allocate(lulcc_co2_oxid_depth_unsat_old(nl_soil,size(co2_oxid_depth_unsat,2))); lulcc_co2_oxid_depth_unsat_old = co2_oxid_depth_unsat
      allocate(lulcc_co2_oxid_depth_sat_old(nl_soil,size(co2_oxid_depth_sat,2))); lulcc_co2_oxid_depth_sat_old = co2_oxid_depth_sat
      allocate(lulcc_co2_aere_depth_unsat_old(nl_soil,size(co2_aere_depth_unsat,2))); lulcc_co2_aere_depth_unsat_old = co2_aere_depth_unsat
      allocate(lulcc_co2_aere_depth_sat_old(nl_soil,size(co2_aere_depth_sat,2))); lulcc_co2_aere_depth_sat_old = co2_aere_depth_sat
      allocate(lulcc_co2_decomp_depth_lake_old(nl_soil,size(co2_decomp_depth_lake,2))); lulcc_co2_decomp_depth_lake_old = co2_decomp_depth_lake
      allocate(lulcc_co2_oxid_depth_lake_old(nl_soil,size(co2_oxid_depth_lake,2))); lulcc_co2_oxid_depth_lake_old = co2_oxid_depth_lake
      allocate(lulcc_conc_o2_unsat_old(nl_soil,size(conc_o2_unsat,2))); lulcc_conc_o2_unsat_old = conc_o2_unsat
      allocate(lulcc_conc_o2_sat_old(nl_soil,size(conc_o2_sat,2))); lulcc_conc_o2_sat_old = conc_o2_sat
      allocate(lulcc_conc_methane_unsat_old(nl_soil,size(conc_methane_unsat,2))); lulcc_conc_methane_unsat_old = conc_methane_unsat
      allocate(lulcc_conc_methane_sat_old(nl_soil,size(conc_methane_sat,2))); lulcc_conc_methane_sat_old = conc_methane_sat
      allocate(lulcc_conc_o2_lake_old(nl_soil,size(conc_o2_lake,2))); lulcc_conc_o2_lake_old = conc_o2_lake
      allocate(lulcc_conc_methane_lake_old(nl_soil,size(conc_methane_lake,2))); lulcc_conc_methane_lake_old = conc_methane_lake
      allocate(lulcc_totcol_methane_old(size(totcol_methane))); lulcc_totcol_methane_old = totcol_methane
      allocate(lulcc_co2_decomp_tot_old(size(co2_decomp_tot))); lulcc_co2_decomp_tot_old = co2_decomp_tot
      allocate(lulcc_co2_oxid_tot_old(size(co2_oxid_tot))); lulcc_co2_oxid_tot_old = co2_oxid_tot
      allocate(lulcc_co2_aere_tot_old(size(co2_aere_tot))); lulcc_co2_aere_tot_old = co2_aere_tot
      allocate(lulcc_co2_net_tot_old(size(co2_net_tot))); lulcc_co2_net_tot_old = co2_net_tot
      allocate(lulcc_co2_decomp_tot_unsat_old(size(co2_decomp_tot_unsat))); lulcc_co2_decomp_tot_unsat_old = co2_decomp_tot_unsat
      allocate(lulcc_co2_decomp_tot_sat_old(size(co2_decomp_tot_sat))); lulcc_co2_decomp_tot_sat_old = co2_decomp_tot_sat
      allocate(lulcc_co2_oxid_tot_unsat_old(size(co2_oxid_tot_unsat))); lulcc_co2_oxid_tot_unsat_old = co2_oxid_tot_unsat
      allocate(lulcc_co2_oxid_tot_sat_old(size(co2_oxid_tot_sat))); lulcc_co2_oxid_tot_sat_old = co2_oxid_tot_sat
      allocate(lulcc_co2_net_tot_unsat_old(size(co2_net_tot_unsat))); lulcc_co2_net_tot_unsat_old = co2_net_tot_unsat
      allocate(lulcc_co2_net_tot_sat_old(size(co2_net_tot_sat))); lulcc_co2_net_tot_sat_old = co2_net_tot_sat
      allocate(lulcc_co2_decomp_tot_lake_old(size(co2_decomp_tot_lake))); lulcc_co2_decomp_tot_lake_old = co2_decomp_tot_lake
      allocate(lulcc_co2_oxid_tot_lake_old(size(co2_oxid_tot_lake))); lulcc_co2_oxid_tot_lake_old = co2_oxid_tot_lake
      allocate(lulcc_co2_net_tot_lake_old(size(co2_net_tot_lake))); lulcc_co2_net_tot_lake_old = co2_net_tot_lake
      allocate(lulcc_totcol_methane_unsat_old(size(totcol_methane_unsat))); lulcc_totcol_methane_unsat_old = totcol_methane_unsat
      allocate(lulcc_totcol_methane_sat_old(size(totcol_methane_sat))); lulcc_totcol_methane_sat_old = totcol_methane_sat
      allocate(lulcc_totcol_methane_lake_old(size(totcol_methane_lake))); lulcc_totcol_methane_lake_old = totcol_methane_lake
      allocate(lulcc_lake_water_ch4_stock_old(size(lake_water_ch4_stock))); &
         lulcc_lake_water_ch4_stock_old = lake_water_ch4_stock
      allocate(lulcc_lake_water_o2_stock_old(size(lake_water_o2_stock))); &
         lulcc_lake_water_o2_stock_old = lake_water_o2_stock
      allocate(lulcc_lake_frozen_ch4_stock_old(size(lake_frozen_ch4_stock))); &
         lulcc_lake_frozen_ch4_stock_old = lake_frozen_ch4_stock
      allocate(lulcc_lake_frozen_o2_stock_old(size(lake_frozen_o2_stock))); &
         lulcc_lake_frozen_o2_stock_old = lake_frozen_o2_stock
      allocate(lulcc_lake_liquid_fraction_prev_old(size(lake_liquid_fraction_prev))); &
         lulcc_lake_liquid_fraction_prev_old = lake_liquid_fraction_prev
      allocate(lulcc_grnd_methane_cond_old(size(grnd_methane_cond))); lulcc_grnd_methane_cond_old = grnd_methane_cond
      allocate(lulcc_grnd_methane_cond_unsat_old(size(grnd_methane_cond_unsat))); lulcc_grnd_methane_cond_unsat_old = grnd_methane_cond_unsat
      allocate(lulcc_grnd_methane_cond_sat_old(size(grnd_methane_cond_sat))); lulcc_grnd_methane_cond_sat_old = grnd_methane_cond_sat
      allocate(lulcc_grnd_methane_cond_lake_old(size(grnd_methane_cond_lake))); lulcc_grnd_methane_cond_lake_old = grnd_methane_cond_lake
      allocate(lulcc_layer_sat_lag_old(nl_soil,size(layer_sat_lag,2))); lulcc_layer_sat_lag_old = layer_sat_lag
      allocate(lulcc_lake_soilc_old(nl_soil,size(lake_soilc,2))); lulcc_lake_soilc_old = lake_soilc
      allocate(lulcc_annavg_agnpp_old(size(annavg_agnpp))); lulcc_annavg_agnpp_old = annavg_agnpp
      allocate(lulcc_annavg_bgnpp_old(size(annavg_bgnpp))); lulcc_annavg_bgnpp_old = annavg_bgnpp
      allocate(lulcc_annavg_somhr_old(size(annavg_somhr))); lulcc_annavg_somhr_old = annavg_somhr
      allocate(lulcc_annavg_finrw_old(size(annavg_finrw))); lulcc_annavg_finrw_old = annavg_finrw
      allocate(lulcc_tempavg_agnpp_old(size(tempavg_agnpp))); lulcc_tempavg_agnpp_old = tempavg_agnpp
      allocate(lulcc_tempavg_bgnpp_old(size(tempavg_bgnpp))); lulcc_tempavg_bgnpp_old = tempavg_bgnpp
      allocate(lulcc_annsum_counter_old(size(annsum_counter))); lulcc_annsum_counter_old = annsum_counter
      allocate(lulcc_tempavg_somhr_old(size(tempavg_somhr))); lulcc_tempavg_somhr_old = tempavg_somhr
      allocate(lulcc_tempavg_finrw_old(size(tempavg_finrw))); lulcc_tempavg_finrw_old = tempavg_finrw
      allocate(lulcc_fsat_bef_old(size(fsat_bef))); lulcc_fsat_bef_old = fsat_bef
      allocate(lulcc_finundated_lag_old(size(finundated_lag))); lulcc_finundated_lag_old = finundated_lag
      allocate(lulcc_methane_dfsat_tot_old(size(methane_dfsat_tot))); lulcc_methane_dfsat_tot_old = methane_dfsat_tot
      allocate(lulcc_f_h2osfc_old(size(f_h2osfc))); lulcc_f_h2osfc_old = f_h2osfc
      allocate(lulcc_forc_pmethanem_old(size(forc_pmethanem))); lulcc_forc_pmethanem_old = forc_pmethanem
	      allocate(lulcc_f_inund_levee_patch_old(size(f_inund_levee_patch)))
	      lulcc_f_inund_levee_patch_old = f_inund_levee_patch
	      allocate(lulcc_f_inund_flood_patch_old(size(f_inund_flood_patch)))
	      lulcc_f_inund_flood_patch_old = f_inund_flood_patch
	      allocate(lulcc_f_inund_flood_depth_patch_old(size(f_inund_flood_depth_patch)))
	      lulcc_f_inund_flood_depth_patch_old = f_inund_flood_depth_patch
      allocate(lulcc_c_atm_old(3,size(c_atm,2))); lulcc_c_atm_old = c_atm
	  allocate(lulcc_conc_o2_unsat_component_old(nl_soil,N_METHANE_COMP,size(conc_o2_unsat_component,3)))
	  lulcc_conc_o2_unsat_component_old = conc_o2_unsat_component
	  allocate(lulcc_conc_o2_sat_component_old(nl_soil,N_METHANE_COMP,size(conc_o2_sat_component,3)))
	  lulcc_conc_o2_sat_component_old = conc_o2_sat_component
	  allocate(lulcc_conc_methane_unsat_component_old(nl_soil,N_METHANE_COMP,size(conc_methane_unsat_component,3)))
	  lulcc_conc_methane_unsat_component_old = conc_methane_unsat_component
	  allocate(lulcc_conc_methane_sat_component_old(nl_soil,N_METHANE_COMP,size(conc_methane_sat_component,3)))
	  lulcc_conc_methane_sat_component_old = conc_methane_sat_component
	  allocate(lulcc_layer_sat_lag_component_old(nl_soil,N_METHANE_COMP,size(layer_sat_lag_component,3)))
	  lulcc_layer_sat_lag_component_old = layer_sat_lag_component
	  allocate(lulcc_annavg_agnpp_component_old(N_METHANE_COMP,size(annavg_agnpp_component,2)))
	  lulcc_annavg_agnpp_component_old = annavg_agnpp_component
	  allocate(lulcc_annavg_bgnpp_component_old(N_METHANE_COMP,size(annavg_bgnpp_component,2)))
	  lulcc_annavg_bgnpp_component_old = annavg_bgnpp_component
	  allocate(lulcc_annavg_somhr_component_old(N_METHANE_COMP,size(annavg_somhr_component,2)))
	  lulcc_annavg_somhr_component_old = annavg_somhr_component
	  allocate(lulcc_annavg_finrw_component_old(N_METHANE_COMP,size(annavg_finrw_component,2)))
	  lulcc_annavg_finrw_component_old = annavg_finrw_component
	  allocate(lulcc_tempavg_agnpp_component_old(N_METHANE_COMP,size(tempavg_agnpp_component,2)))
	  lulcc_tempavg_agnpp_component_old = tempavg_agnpp_component
	  allocate(lulcc_tempavg_bgnpp_component_old(N_METHANE_COMP,size(tempavg_bgnpp_component,2)))
	  lulcc_tempavg_bgnpp_component_old = tempavg_bgnpp_component
	  allocate(lulcc_annsum_counter_component_old(N_METHANE_COMP,size(annsum_counter_component,2)))
	  lulcc_annsum_counter_component_old = annsum_counter_component
	  allocate(lulcc_tempavg_somhr_component_old(N_METHANE_COMP,size(tempavg_somhr_component,2)))
	  lulcc_tempavg_somhr_component_old = tempavg_somhr_component
	  allocate(lulcc_tempavg_finrw_component_old(N_METHANE_COMP,size(tempavg_finrw_component,2)))
	  lulcc_tempavg_finrw_component_old = tempavg_finrw_component
	  allocate(lulcc_fsat_bef_component_old(N_METHANE_COMP,size(fsat_bef_component,2)))
	  lulcc_fsat_bef_component_old = fsat_bef_component
	  allocate(lulcc_finundated_lag_component_old(N_METHANE_COMP,size(finundated_lag_component,2)))
	  lulcc_finundated_lag_component_old = finundated_lag_component
	  allocate(lulcc_rice_fraction_prev_old(size(rice_fraction_prev)))
	  lulcc_rice_fraction_prev_old = rice_fraction_prev

	  ! Non-soil patches still use the legacy aggregate column.  Mirror that
	  ! live state into the snapshot's zero-rice component before LULCC so a
	  ! wetland/lake-to-soil transition cannot inherit stale hidden columns.
	  nsave = min(size(patchtype), size(lulcc_rice_fraction_prev_old))
	  DO ip = 1, nsave
	     IF (patchtype(ip) == 0) CYCLE
	     lulcc_rice_fraction_prev_old(ip) = 0._r8
	     DO component = 1, N_METHANE_COMP
	        lulcc_conc_o2_unsat_component_old(:,component,ip) = lulcc_conc_o2_unsat_old(:,ip)
	        lulcc_conc_o2_sat_component_old(:,component,ip) = lulcc_conc_o2_sat_old(:,ip)
	        lulcc_conc_methane_unsat_component_old(:,component,ip) = lulcc_conc_methane_unsat_old(:,ip)
	        lulcc_conc_methane_sat_component_old(:,component,ip) = lulcc_conc_methane_sat_old(:,ip)
	        lulcc_layer_sat_lag_component_old(:,component,ip) = lulcc_layer_sat_lag_old(:,ip)
	        lulcc_annavg_agnpp_component_old(component,ip) = lulcc_annavg_agnpp_old(ip)
	        lulcc_annavg_bgnpp_component_old(component,ip) = lulcc_annavg_bgnpp_old(ip)
	        lulcc_annavg_somhr_component_old(component,ip) = lulcc_annavg_somhr_old(ip)
	        lulcc_annavg_finrw_component_old(component,ip) = lulcc_annavg_finrw_old(ip)
	        lulcc_tempavg_agnpp_component_old(component,ip) = lulcc_tempavg_agnpp_old(ip)
	        lulcc_tempavg_bgnpp_component_old(component,ip) = lulcc_tempavg_bgnpp_old(ip)
	        lulcc_annsum_counter_component_old(component,ip) = lulcc_annsum_counter_old(ip)
	        lulcc_tempavg_somhr_component_old(component,ip) = lulcc_tempavg_somhr_old(ip)
	        lulcc_tempavg_finrw_component_old(component,ip) = lulcc_tempavg_finrw_old(ip)
	        lulcc_fsat_bef_component_old(component,ip) = lulcc_fsat_bef_old(ip)
	        lulcc_finundated_lag_component_old(component,ip) = lulcc_finundated_lag_old(ip)
	     ENDDO
	  ENDDO

      methane_lulcc_snapshot_valid = .true.
   END SUBROUTINE save_methane_lulcc_state

   SUBROUTINE publish_methane_levee_flood_patch (fldfrc_patch)

      IMPLICIT NONE
      real(r8), intent(in) :: fldfrc_patch(:)

      integer :: ncopy

      IF (.not. allocated(f_inund_levee_patch)) RETURN

      ncopy = min(size(f_inund_levee_patch), size(fldfrc_patch))
      IF (ncopy > 0) f_inund_levee_patch(1:ncopy) = fldfrc_patch(1:ncopy)
      IF (size(f_inund_levee_patch) > ncopy) f_inund_levee_patch(ncopy+1:) = 0._r8

   END SUBROUTINE publish_methane_levee_flood_patch

   SUBROUTINE publish_methane_flood_patch (fldfrc_patch, flddph_patch)

      IMPLICIT NONE
      real(r8), intent(in) :: fldfrc_patch(:)
      real(r8), intent(in) :: flddph_patch(:)

      integer :: ncopy

      IF (.not. allocated(f_inund_flood_patch)) RETURN

      ncopy = min(size(f_inund_flood_patch), size(fldfrc_patch))
      IF (ncopy > 0) f_inund_flood_patch(1:ncopy) = fldfrc_patch(1:ncopy)
      IF (size(f_inund_flood_patch) > ncopy) f_inund_flood_patch(ncopy+1:) = 0._r8

      IF (allocated(f_inund_flood_depth_patch)) THEN
         ncopy = min(size(f_inund_flood_depth_patch), size(flddph_patch))
         IF (ncopy > 0) f_inund_flood_depth_patch(1:ncopy) = flddph_patch(1:ncopy)
         IF (size(f_inund_flood_depth_patch) > ncopy) f_inund_flood_depth_patch(ncopy+1:) = 0._r8
      ENDIF

   END SUBROUTINE publish_methane_flood_patch


	   SUBROUTINE remap_methane_lulcc_state (patchclass_new, eindex_new, patchclass_old, eindex_old, &
	      lccpct_patches, new_patch_area, old_patch_area)
	      USE MOD_Vars_TimeInvariants, only: patchtype
	      USE MOD_SPMD_Task, only: CoLM_stop
	      integer, intent(in) :: patchclass_new(:), patchclass_old(:)
	      integer*8, intent(in) :: eindex_new(:), eindex_old(:)
      real(r8), intent(in), optional :: lccpct_patches(:,0:)
	      real(r8), intent(in), optional :: new_patch_area(:)
	      real(r8), intent(in), optional :: old_patch_area(:)
	      integer :: nnew, np, op, link
	      integer, allocatable :: map_start(:), map_old(:), fallback_map(:)
	      real(r8), allocatable :: map_source_weight(:), map_mass_weight(:)
	      logical, allocatable :: map_mass_available(:), initialize_lake_from_surface(:)

	      nnew = size(patchclass_new)
	      IF (allocated(conc_methane)) CALL deallocate_methane_state ()
	      CALL allocate_methane_state (nnew)
	      CALL init_methane_wetland_fraction_cache (nnew)

		      IF (.not. methane_lulcc_snapshot_valid) THEN
         IF (allocated(patchtype)) THEN
            CALL initialize_methane_lake_soilc_from_surface (patchtype, methane_lake_soilc_srf, &
               DEF_METHANE%allowlakeprod)
         ENDIF
         RETURN
      ENDIF

      CALL build_lulcc_remap_map ()

      ! LulccInitialize has already loaded this year's lake soil carbon.  Mark
      ! only water patches without a contributing old water patch so an
      ! exhausted existing lake is not refilled by the annual LULCC cycle.
      allocate(initialize_lake_from_surface(nnew))
      initialize_lake_from_surface = .false.
      DO np = 1, min(nnew, size(patchclass_new))
         IF (patchclass_new(np) /= WATERBODY) CYCLE
         initialize_lake_from_surface(np) = .true.
         DO link = map_start(np), map_start(np+1)-1
            op = map_old(link)
            IF (op > size(patchclass_old)) CYCLE
            IF (patchclass_old(op) /= WATERBODY) CYCLE
            IF (present(lccpct_patches)) THEN
               IF (map_source_weight(link) <= 0._r8) CYCLE
            ENDIF
            initialize_lake_from_surface(np) = .false.
            EXIT
         ENDDO
      ENDDO

      CALL remap2d(lulcc_conc_o2_old,              conc_o2)
      CALL remap2d(lulcc_conc_methane_old,         conc_methane)
      CALL remap2d(lulcc_co2_decomp_depth_old,     co2_decomp_depth)
      CALL remap2d(lulcc_co2_oxid_depth_old,       co2_oxid_depth)
      CALL remap2d(lulcc_co2_aere_depth_old,       co2_aere_depth)
      CALL remap2d(lulcc_co2_decomp_depth_unsat_old, co2_decomp_depth_unsat)
      CALL remap2d(lulcc_co2_decomp_depth_sat_old,   co2_decomp_depth_sat)
      CALL remap2d(lulcc_co2_oxid_depth_unsat_old, co2_oxid_depth_unsat)
      CALL remap2d(lulcc_co2_oxid_depth_sat_old,   co2_oxid_depth_sat)
      CALL remap2d(lulcc_co2_aere_depth_unsat_old, co2_aere_depth_unsat)
      CALL remap2d(lulcc_co2_aere_depth_sat_old,   co2_aere_depth_sat)
      CALL remap2d(lulcc_co2_decomp_depth_lake_old, co2_decomp_depth_lake)
      CALL remap2d(lulcc_co2_oxid_depth_lake_old,  co2_oxid_depth_lake)
      CALL remap2d(lulcc_conc_o2_unsat_old,        conc_o2_unsat)
      CALL remap2d(lulcc_conc_o2_sat_old,          conc_o2_sat)
      CALL remap2d(lulcc_conc_methane_unsat_old,   conc_methane_unsat)
      CALL remap2d(lulcc_conc_methane_sat_old,     conc_methane_sat)
      CALL remap2d(lulcc_conc_o2_lake_old,         conc_o2_lake)
      CALL remap2d(lulcc_conc_methane_lake_old,    conc_methane_lake)
      CALL remap2d(lulcc_layer_sat_lag_old,        layer_sat_lag)
      CALL remap2d(lulcc_lake_soilc_old,           lake_soilc)
      IF (allocated(patchtype)) THEN
         CALL initialize_methane_lake_soilc_from_surface (patchtype, methane_lake_soilc_srf, &
            DEF_METHANE%allowlakeprod, initialize_lake_from_surface)
      ENDIF
      CALL remap2d(lulcc_c_atm_old,                c_atm)
	  CALL remap_component_fraction(lulcc_rice_fraction_prev_old, rice_fraction_prev)
	  CALL remap_component_phase(lulcc_conc_o2_unsat_component_old, &
	     lulcc_conc_o2_sat_component_old, lulcc_fsat_bef_component_old, &
	     conc_o2_unsat_component, conc_o2_sat_component)
	  CALL remap_component_phase(lulcc_conc_methane_unsat_component_old, &
	     lulcc_conc_methane_sat_component_old, lulcc_fsat_bef_component_old, &
	     conc_methane_unsat_component, conc_methane_sat_component)
	  CALL remap_component_3d(lulcc_layer_sat_lag_component_old, layer_sat_lag_component)
	  CALL remap_component_2d(lulcc_annavg_agnpp_component_old, annavg_agnpp_component)
	  CALL remap_component_2d(lulcc_annavg_bgnpp_component_old, annavg_bgnpp_component)
	  CALL remap_component_2d(lulcc_annavg_somhr_component_old, annavg_somhr_component)
	  CALL remap_component_2d(lulcc_annavg_finrw_component_old, annavg_finrw_component)
	  CALL remap_component_2d(lulcc_tempavg_agnpp_component_old, tempavg_agnpp_component)
	  CALL remap_component_2d(lulcc_tempavg_bgnpp_component_old, tempavg_bgnpp_component)
	  CALL remap_component_2d(lulcc_annsum_counter_component_old, annsum_counter_component)
	  CALL remap_component_2d(lulcc_tempavg_somhr_component_old, tempavg_somhr_component)
	  CALL remap_component_2d(lulcc_tempavg_finrw_component_old, tempavg_finrw_component)
	  CALL remap_component_2d(lulcc_fsat_bef_component_old, fsat_bef_component)
	  CALL remap_component_2d(lulcc_finundated_lag_component_old, finundated_lag_component)
	  rice_fraction_prev = min(max(rice_fraction_prev, 0._r8), 1._r8)
	  WHERE (.not. invalid_restart_value(fsat_bef_component))
	     fsat_bef_component = min(max(fsat_bef_component, 0._r8), 1._r8)
	  ELSEWHERE
	     fsat_bef_component = spval
	  END WHERE
	  WHERE (.not. invalid_restart_value(finundated_lag_component))
	     finundated_lag_component = min(max(finundated_lag_component, 0._r8), 1._r8)
	  ELSEWHERE
	     finundated_lag_component = spval
	  END WHERE

      ! Concentrations, conductances, fractions, rates, and rolling diagnostics
      ! are remapped as intensive fields.  Column stocks and total per-area
      ! source/flux diagnostics are remapped as area-mass fields when both old
      ! and new patch areas are available, matching land-water material LULCC
      ! semantics for extensive per-area pools.
      CALL remap1d_mass(lulcc_totcol_methane_old,       totcol_methane)
      CALL remap1d_mass(lulcc_co2_decomp_tot_old,       co2_decomp_tot)
      CALL remap1d_mass(lulcc_co2_oxid_tot_old,         co2_oxid_tot)
      CALL remap1d_mass(lulcc_co2_aere_tot_old,         co2_aere_tot)
      CALL remap1d_mass(lulcc_co2_net_tot_old,          co2_net_tot)
      CALL remap1d_mass(lulcc_co2_decomp_tot_unsat_old, co2_decomp_tot_unsat)
      CALL remap1d_mass(lulcc_co2_decomp_tot_sat_old,   co2_decomp_tot_sat)
      CALL remap1d_mass(lulcc_co2_oxid_tot_unsat_old,   co2_oxid_tot_unsat)
      CALL remap1d_mass(lulcc_co2_oxid_tot_sat_old,     co2_oxid_tot_sat)
      CALL remap1d_mass(lulcc_co2_net_tot_unsat_old,    co2_net_tot_unsat)
      CALL remap1d_mass(lulcc_co2_net_tot_sat_old,      co2_net_tot_sat)
      CALL remap1d_mass(lulcc_co2_decomp_tot_lake_old,  co2_decomp_tot_lake)
      CALL remap1d_mass(lulcc_co2_oxid_tot_lake_old,    co2_oxid_tot_lake)
      CALL remap1d_mass(lulcc_co2_net_tot_lake_old,     co2_net_tot_lake)
      CALL remap1d_mass(lulcc_totcol_methane_unsat_old, totcol_methane_unsat)
      CALL remap1d_mass(lulcc_totcol_methane_sat_old,   totcol_methane_sat)
      CALL remap1d_mass(lulcc_totcol_methane_lake_old,  totcol_methane_lake)
      CALL remap1d_mass(lulcc_lake_water_ch4_stock_old, lake_water_ch4_stock)
      CALL remap1d_mass(lulcc_lake_water_o2_stock_old, lake_water_o2_stock)
      CALL remap1d_mass(lulcc_lake_frozen_ch4_stock_old, lake_frozen_ch4_stock)
      CALL remap1d_mass(lulcc_lake_frozen_o2_stock_old, lake_frozen_o2_stock)
      CALL remap1d(lulcc_lake_liquid_fraction_prev_old, lake_liquid_fraction_prev)
	  WHERE (invalid_restart_value(lake_liquid_fraction_prev) .or. &
	         lake_liquid_fraction_prev < 0._r8 .or. lake_liquid_fraction_prev > 1._r8) &
	     lake_liquid_fraction_prev = spval
	  WHERE (initialize_lake_from_surface)
	     lake_frozen_ch4_stock = 0._r8
	     lake_frozen_o2_stock = 0._r8
	     lake_liquid_fraction_prev = spval
	  END WHERE
      CALL remap1d(lulcc_fsat_bef_old,             fsat_bef)
      CALL repartition_ch4_totcol_after_lulcc()
      ! Keep remapped CH4 column stocks and layer concentrations consistent
      ! across land<->lake LULCC so the next balance step does not emit a
      ! one-time artificial residual flux.
      CALL sync_ch4_conc_to_totcol(conc_methane,       totcol_methane)
      CALL sync_ch4_conc_to_totcol(conc_methane_unsat, totcol_methane_unsat)
      CALL sync_ch4_conc_to_totcol(conc_methane_sat,   totcol_methane_sat)
      CALL sync_ch4_conc_to_totcol(conc_methane_lake,  totcol_methane_lake)
      CALL remap1d(lulcc_grnd_methane_cond_old,    grnd_methane_cond)
      CALL remap1d(lulcc_grnd_methane_cond_unsat_old, grnd_methane_cond_unsat)
      CALL remap1d(lulcc_grnd_methane_cond_sat_old,   grnd_methane_cond_sat)
      CALL remap1d(lulcc_grnd_methane_cond_lake_old,  grnd_methane_cond_lake)
      CALL remap1d(lulcc_annavg_agnpp_old,         annavg_agnpp)
      CALL remap1d(lulcc_annavg_bgnpp_old,         annavg_bgnpp)
      CALL remap1d(lulcc_annavg_somhr_old,         annavg_somhr)
      CALL remap1d(lulcc_annavg_finrw_old,         annavg_finrw)
      CALL remap1d(lulcc_tempavg_agnpp_old,        tempavg_agnpp)
      CALL remap1d(lulcc_tempavg_bgnpp_old,        tempavg_bgnpp)
      CALL remap1d(lulcc_annsum_counter_old,       annsum_counter)
      CALL remap1d(lulcc_tempavg_somhr_old,        tempavg_somhr)
      CALL remap1d(lulcc_tempavg_finrw_old,        tempavg_finrw)
      CALL remap1d(lulcc_finundated_lag_old,       finundated_lag)
      CALL remap1d_mass(lulcc_methane_dfsat_tot_old,    methane_dfsat_tot)
      CALL remap1d(lulcc_f_h2osfc_old,             f_h2osfc)
      CALL remap1d(lulcc_forc_pmethanem_old,       forc_pmethanem)
	      CALL remap1d(lulcc_f_inund_levee_patch_old, f_inund_levee_patch)
	      CALL remap1d(lulcc_f_inund_flood_patch_old, f_inund_flood_patch)
	      CALL remap1d(lulcc_f_inund_flood_depth_patch_old, f_inund_flood_depth_patch)
	      WHERE (ieee_is_nan(f_inund_levee_patch) .or. &
	             abs(f_inund_levee_patch) >= 0.5_r8 * abs(spval))
	         f_inund_levee_patch = 0._r8
	      END WHERE
	      WHERE (ieee_is_nan(f_inund_flood_patch) .or. &
	             abs(f_inund_flood_patch) >= 0.5_r8 * abs(spval))
	         f_inund_flood_patch = 0._r8
	      END WHERE
	      WHERE (ieee_is_nan(f_inund_flood_depth_patch) .or. &
	             abs(f_inund_flood_depth_patch) >= 0.5_r8 * abs(spval))
	         f_inund_flood_depth_patch = 0._r8
	      END WHERE
	      f_inund_levee_patch = min(max(f_inund_levee_patch, 0._r8), 1._r8)
	      f_inund_flood_patch = min(max(f_inund_flood_patch, 0._r8), 1._r8)
	      f_inund_flood_depth_patch = max(f_inund_flood_depth_patch, 0._r8)

	  CALL sync_component_aggregates_after_lulcc()

      CALL clear_methane_lulcc_snapshot ()

   CONTAINS
	  SUBROUTINE build_lulcc_remap_map ()
	     integer :: np, op, nold, nlink, link, c, rep, class_lo, class_hi
	     real(r8) :: base_weight, class_sum, target_area, denom
	     real(r8), allocatable :: old_group_class_area(:,:)
	     real(r8), allocatable :: new_target_class_area(:,:)
	     real(r8), allocatable :: target_group_class_area(:,:)
	     logical, allocatable :: old_group_area_ready(:)

	     nold = min(size(patchclass_old), size(eindex_old))
	     nlink = 0
	     DO np = 1, min(nnew, size(eindex_new))
	        DO op = 1, nold
	           IF (eindex_old(op) == eindex_new(np)) nlink = nlink + 1
	        ENDDO
	     ENDDO

	     allocate(map_start(nnew+1), fallback_map(nnew), map_mass_available(nnew))
	     allocate(map_old(nlink), map_source_weight(nlink), map_mass_weight(nlink))
	     map_start = 1
	     fallback_map = 0
	     map_mass_available = .false.
	     map_source_weight = 0._r8
	     map_mass_weight = 0._r8

	     link = 1
	     DO np = 1, nnew
	        map_start(np) = link
	        map_mass_available(np) = area_mass_remap_available(np)
	        IF (np <= size(eindex_new)) THEN
	           DO op = 1, nold
	              IF (eindex_old(op) /= eindex_new(np)) CYCLE
	              map_old(link) = op
	              IF (fallback_map(np) == 0 .and. np <= size(patchclass_new)) THEN
	                 IF (patchclass_old(op) == patchclass_new(np)) fallback_map(np) = op
	              ENDIF
	              link = link + 1
	           ENDDO
	        ENDIF
	     ENDDO
	     map_start(nnew+1) = link

	     IF (.not. present(lccpct_patches)) RETURN
	     class_lo = lbound(lccpct_patches,2)
	     class_hi = ubound(lccpct_patches,2)
	     allocate(old_group_class_area(nold,class_lo:class_hi))
	     allocate(new_target_class_area(nnew,class_lo:class_hi))
	     allocate(target_group_class_area(nold,class_lo:class_hi))
	     allocate(old_group_area_ready(nold))
	     old_group_class_area = 0._r8
	     new_target_class_area = 0._r8
	     target_group_class_area = 0._r8
	     old_group_area_ready = .false.

	     ! Precompute each new patch's target class areas once.  The mass
	     ! remap denominator is then accumulated by eindex group using the
	     ! first linked old patch as a stable group representative.
	     IF (present(new_patch_area)) THEN
	        DO np = 1, min(nnew, size(lccpct_patches,1), size(new_patch_area))
	           class_sum = sum(max(lccpct_patches(np,class_lo:class_hi), 0._r8))
	           IF (class_sum <= tiny(1._r8)) CYCLE
	           DO c = class_lo, class_hi
	              new_target_class_area(np,c) = max(0._r8, new_patch_area(np)) * &
	                 max(0._r8, lccpct_patches(np,c)) / class_sum
	           ENDDO
	        ENDDO
	     ENDIF
	     DO np = 1, nnew
	        IF (map_start(np) >= map_start(np+1)) CYCLE
	        rep = map_old(map_start(np))
	        target_group_class_area(rep,:) = target_group_class_area(rep,:) + &
	           new_target_class_area(np,:)
	     ENDDO

	     ! Sum old source area once per eindex/class group.  Every new patch in
	     ! the group has the same link set, so later link weights are O(1).
	     IF (present(old_patch_area)) THEN
	        DO np = 1, nnew
	           IF (map_start(np) >= map_start(np+1)) CYCLE
	           rep = map_old(map_start(np))
	           IF (old_group_area_ready(rep)) CYCLE
	           DO link = map_start(np), map_start(np+1)-1
	              op = map_old(link)
	              IF (op > size(old_patch_area)) CYCLE
	              c = patchclass_old(op)
	              IF (c < class_lo .or. c > class_hi) CYCLE
	              old_group_class_area(rep,c) = old_group_class_area(rep,c) + &
	                 max(0._r8, old_patch_area(op))
	           ENDDO
	           old_group_area_ready(rep) = .true.
	        ENDDO
	     ENDIF

	     DO np = 1, min(nnew, size(lccpct_patches,1))
	        IF (map_start(np) >= map_start(np+1)) CYCLE
	        rep = map_old(map_start(np))
	        DO link = map_start(np), map_start(np+1)-1
	           op = map_old(link)
	           c = patchclass_old(op)
	           IF (c < class_lo .or. c > class_hi) CYCLE
	           base_weight = max(0._r8, lccpct_patches(np,c))
	           map_source_weight(link) = base_weight
	           IF (base_weight > 0._r8 .and. present(old_patch_area)) THEN
	              IF (op <= size(old_patch_area)) THEN
	                 denom = old_group_class_area(rep,c)
	                 IF (denom > 0._r8) map_source_weight(link) = base_weight * &
	                    max(0._r8, old_patch_area(op)) / denom
	              ENDIF
	           ENDIF

	           IF (map_mass_available(np)) THEN
	              IF (op <= size(old_patch_area)) THEN
	                 target_area = new_target_class_area(np,c)
	                 denom = target_group_class_area(rep,c)
	                 IF (target_area > tiny(1._r8) .and. denom > tiny(1._r8)) THEN
	                    map_mass_weight(link) = max(0._r8, old_patch_area(op)) * &
	                       target_area / denom
	                 ENDIF
	              ENDIF
	           ENDIF
	        ENDDO
	     ENDDO
	  END SUBROUTINE build_lulcc_remap_map

	  SUBROUTINE remap_component_fraction(old_fraction, new_fraction)
	     real(r8), intent(in) :: old_fraction(:)
	     real(r8), intent(inout) :: new_fraction(:)
	     integer :: np, op, src, link
	     real(r8) :: w, wsum, rice_area

	     DO np = 1, min(size(new_fraction), nnew)
	        wsum = 0._r8
	        rice_area = 0._r8
	        IF (present(lccpct_patches)) THEN
	           DO link = map_start(np), map_start(np+1)-1
	              op = map_old(link)
	              IF (op > size(old_fraction)) CYCLE
	              w = component_base_weight(np, link)
	              IF (w <= 0._r8) CYCLE
	              rice_area = rice_area + w*min(max(old_fraction(op), 0._r8), 1._r8)
	              wsum = wsum + w
	           ENDDO
	        ENDIF
	        IF (wsum > tiny(1._r8)) THEN
	           new_fraction(np) = rice_area/wsum
	        ELSE
	           src = fallback_source(np, size(old_fraction))
	           IF (src > 0) new_fraction(np) = old_fraction(src)
	        ENDIF
	     ENDDO
	  END SUBROUTINE remap_component_fraction

	  SUBROUTINE remap_component_2d(old, new)
	     real(r8), intent(in) :: old(:,:)
	     real(r8), intent(inout) :: new(:,:)
	     integer :: np, op, component, src, link
	     real(r8) :: w, wsum, val

	     DO np = 1, min(size(new,2), nnew)
	        DO component = 1, min(size(new,1), size(old,1), N_METHANE_COMP)
	           wsum = 0._r8
	           val = 0._r8
	           IF (present(lccpct_patches)) THEN
	              DO link = map_start(np), map_start(np+1)-1
	                 op = map_old(link)
	                 IF (op > size(old,2)) CYCLE
	                 w = component_transfer_weight(np, link, component)
	                 IF (w <= 0._r8) CYCLE
	                 val = val + w*old(component,op)
	                 wsum = wsum + w
	              ENDDO
	           ENDIF
	           IF (wsum > tiny(1._r8)) THEN
	              new(component,np) = val/wsum
	           ELSE
	              src = fallback_source(np, size(old,2))
	              IF (src > 0) new(component,np) = old(component,src)
	           ENDIF
	        ENDDO
	     ENDDO
	  END SUBROUTINE remap_component_2d

	  SUBROUTINE remap_component_3d(old, new)
	     real(r8), intent(in) :: old(:,:,:)
	     real(r8), intent(inout) :: new(:,:,:)
	     integer :: np, op, component, src, nlev, link
	     real(r8) :: w, wsum
	     real(r8) :: val(size(new,1))

	     nlev = min(size(new,1), size(old,1))
	     DO np = 1, min(size(new,3), nnew)
	        DO component = 1, min(size(new,2), size(old,2), N_METHANE_COMP)
	           wsum = 0._r8
	           val = 0._r8
	           IF (present(lccpct_patches)) THEN
	              DO link = map_start(np), map_start(np+1)-1
	                 op = map_old(link)
	                 IF (op > size(old,3)) CYCLE
	                 w = component_transfer_weight(np, link, component)
	                 IF (w <= 0._r8) CYCLE
	                 val(1:nlev) = val(1:nlev) + w*old(1:nlev,component,op)
	                 wsum = wsum + w
	              ENDDO
	           ENDIF
	           IF (wsum > tiny(1._r8)) THEN
	              new(1:nlev,component,np) = val(1:nlev)/wsum
	           ELSE
	              src = fallback_source(np, size(old,3))
	              IF (src > 0) new(1:nlev,component,np) = old(1:nlev,component,src)
	           ENDIF
	        ENDDO
	     ENDDO
	  END SUBROUTINE remap_component_3d

	  SUBROUTINE remap_component_phase(old_unsat, old_sat, old_fsat, new_unsat, new_sat)
	     real(r8), intent(in) :: old_unsat(:,:,:), old_sat(:,:,:), old_fsat(:,:)
	     real(r8), intent(inout) :: new_unsat(:,:,:), new_sat(:,:,:)
	     integer :: np, op, component, src, nlev, link
	     real(r8) :: w, h, wunsat, wsat
	     real(r8) :: val_unsat(size(new_unsat,1)), val_sat(size(new_sat,1))

	     nlev = min(size(new_unsat,1), size(new_sat,1), size(old_unsat,1), size(old_sat,1))
	     DO np = 1, min(size(new_unsat,3), size(new_sat,3), nnew)
	        DO component = 1, min(size(new_unsat,2), size(new_sat,2), &
	           size(old_unsat,2), size(old_sat,2), N_METHANE_COMP)
	           wunsat = 0._r8
	           wsat = 0._r8
	           val_unsat = 0._r8
	           val_sat = 0._r8
	           IF (present(lccpct_patches)) THEN
	              DO link = map_start(np), map_start(np+1)-1
	                 op = map_old(link)
	                 IF (op > size(old_unsat,3) .or. op > size(old_sat,3) .or. &
	                     op > size(old_fsat,2)) CYCLE
	                 w = component_transfer_weight(np, link, component)
	                 IF (w <= 0._r8) CYCLE
	                 h = old_fsat(component,op)
	                 IF (invalid_restart_value(h) .or. h < 0._r8 .or. h > 1._r8) h = 0.5_r8
	                 val_unsat(1:nlev) = val_unsat(1:nlev) + &
	                    w*(1._r8-h)*old_unsat(1:nlev,component,op)
	                 val_sat(1:nlev) = val_sat(1:nlev) + w*h*old_sat(1:nlev,component,op)
	                 wunsat = wunsat + w*(1._r8-h)
	                 wsat = wsat + w*h
	              ENDDO
	           ENDIF
	           IF (wunsat > tiny(1._r8)) new_unsat(1:nlev,component,np) = val_unsat(1:nlev)/wunsat
	           IF (wsat > tiny(1._r8)) new_sat(1:nlev,component,np) = val_sat(1:nlev)/wsat
	           IF (wunsat <= tiny(1._r8) .and. wsat > tiny(1._r8)) &
	              new_unsat(1:nlev,component,np) = new_sat(1:nlev,component,np)
	           IF (wsat <= tiny(1._r8) .and. wunsat > tiny(1._r8)) &
	              new_sat(1:nlev,component,np) = new_unsat(1:nlev,component,np)
	           IF (wunsat <= tiny(1._r8) .and. wsat <= tiny(1._r8)) THEN
	              src = fallback_source(np, size(old_unsat,3))
	              IF (src > 0) THEN
	                 new_unsat(1:nlev,component,np) = old_unsat(1:nlev,component,src)
	                 new_sat(1:nlev,component,np) = old_sat(1:nlev,component,src)
	              ENDIF
	           ENDIF
	        ENDDO
	     ENDDO
	  END SUBROUTINE remap_component_phase

	  SUBROUTINE sync_component_aggregates_after_lulcc()
	     integer :: np, j, nlev
	     real(r8) :: r, ws, hs, hr, sat_area, unsat_area
	     real(r8) :: wss, wsr, wus, wur, dz, target, column, scale

	     nlev = min(nl_soil, size(conc_methane,1))
	     DO np = 1, min(nnew, size(patchclass_new), size(rice_fraction_prev))
	        IF (patchclass_new(np) == WATERBODY) THEN
	           rice_fraction_prev(np) = 0._r8
	           CYCLE
	        ENDIF
	        r = min(max(rice_fraction_prev(np), 0._r8), 1._r8)
	        ws = 1._r8-r
	        hs = fsat_bef_component(METHANE_COMP_SOIL,np)
	        hr = fsat_bef_component(METHANE_COMP_RICE,np)
	        IF (invalid_restart_value(hs) .or. hs < 0._r8 .or. hs > 1._r8) hs = 0.5_r8
	        IF (invalid_restart_value(hr) .or. hr < 0._r8 .or. hr > 1._r8) hr = 0.5_r8
	        sat_area = ws*hs+r*hr
	        unsat_area = ws*(1._r8-hs)+r*(1._r8-hr)
	        IF (sat_area > tiny(1._r8)) THEN
	           wss = ws*hs/sat_area
	           wsr = r*hr/sat_area
	        ELSE
	           wss = ws
	           wsr = r
	        ENDIF
	        IF (unsat_area > tiny(1._r8)) THEN
	           wus = ws*(1._r8-hs)/unsat_area
	           wur = r*(1._r8-hr)/unsat_area
	        ELSE
	           wus = ws
	           wur = r
	        ENDIF

	        ! The component-aware averages above preserve the split between
	        ! rice/non-rice and sat/unsat source areas.  Apply one common CH4
	        ! scale so their reconstructed patch inventory also matches the
	        ! area-mass-conservative legacy total remapped just above.
	        target = max(totcol_methane(np), 0._r8)
	        column = 0._r8
	        DO j = 1, nlev
	           dz = max(dz_soi(j), 0._r8)
	           column = column + dz*( &
	              ws*((1._r8-hs)*conc_methane_unsat_component(j,METHANE_COMP_SOIL,np) + &
	                  hs*conc_methane_sat_component(j,METHANE_COMP_SOIL,np)) + &
	              r*((1._r8-hr)*conc_methane_unsat_component(j,METHANE_COMP_RICE,np) + &
	                  hr*conc_methane_sat_component(j,METHANE_COMP_RICE,np)))
	        ENDDO
	        IF (target <= tiny(1._r8)) THEN
	           conc_methane_unsat_component(:,:,np) = 0._r8
	           conc_methane_sat_component(:,:,np) = 0._r8
	        ELSEIF (column > tiny(1._r8)) THEN
	           scale = target/column
	           conc_methane_unsat_component(:,:,np) = conc_methane_unsat_component(:,:,np)*scale
	           conc_methane_sat_component(:,:,np) = conc_methane_sat_component(:,:,np)*scale
	        ELSE
	           conc_methane_unsat_component(:,:,np) = target/max(sum(max(dz_soi(1:nlev),0._r8)),tiny(1._r8))
	           conc_methane_sat_component(:,:,np) = conc_methane_unsat_component(:,:,np)
	        ENDIF

	        conc_o2_unsat(1:nlev,np) = wus*conc_o2_unsat_component(1:nlev,METHANE_COMP_SOIL,np) + &
	           wur*conc_o2_unsat_component(1:nlev,METHANE_COMP_RICE,np)
	        conc_o2_sat(1:nlev,np) = wss*conc_o2_sat_component(1:nlev,METHANE_COMP_SOIL,np) + &
	           wsr*conc_o2_sat_component(1:nlev,METHANE_COMP_RICE,np)
	        conc_methane_unsat(1:nlev,np) = &
	           wus*conc_methane_unsat_component(1:nlev,METHANE_COMP_SOIL,np) + &
	           wur*conc_methane_unsat_component(1:nlev,METHANE_COMP_RICE,np)
	        conc_methane_sat(1:nlev,np) = &
	           wss*conc_methane_sat_component(1:nlev,METHANE_COMP_SOIL,np) + &
	           wsr*conc_methane_sat_component(1:nlev,METHANE_COMP_RICE,np)
	        conc_o2(1:nlev,np) = ws*((1._r8-hs)*conc_o2_unsat_component(1:nlev,METHANE_COMP_SOIL,np) + &
	           hs*conc_o2_sat_component(1:nlev,METHANE_COMP_SOIL,np)) + &
	           r*((1._r8-hr)*conc_o2_unsat_component(1:nlev,METHANE_COMP_RICE,np) + &
	           hr*conc_o2_sat_component(1:nlev,METHANE_COMP_RICE,np))
	        conc_methane(1:nlev,np) = &
	           ws*((1._r8-hs)*conc_methane_unsat_component(1:nlev,METHANE_COMP_SOIL,np) + &
	           hs*conc_methane_sat_component(1:nlev,METHANE_COMP_SOIL,np)) + &
	           r*((1._r8-hr)*conc_methane_unsat_component(1:nlev,METHANE_COMP_RICE,np) + &
	           hr*conc_methane_sat_component(1:nlev,METHANE_COMP_RICE,np))

	        totcol_methane(np) = 0._r8
	        totcol_methane_unsat(np) = 0._r8
	        totcol_methane_sat(np) = 0._r8
	        DO j = 1, nlev
	           dz = max(dz_soi(j), 0._r8)
	           totcol_methane(np) = totcol_methane(np) + conc_methane(j,np)*dz
	           totcol_methane_unsat(np) = totcol_methane_unsat(np) + conc_methane_unsat(j,np)*dz
	           totcol_methane_sat(np) = totcol_methane_sat(np) + conc_methane_sat(j,np)*dz
	        ENDDO
	        fsat_bef(np) = sat_area
	        layer_sat_lag(1:nlev,np) = ws*layer_sat_lag_component(1:nlev,METHANE_COMP_SOIL,np) + &
	           r*layer_sat_lag_component(1:nlev,METHANE_COMP_RICE,np)
	        annavg_agnpp(np) = ws*annavg_agnpp_component(METHANE_COMP_SOIL,np) + &
	           r*annavg_agnpp_component(METHANE_COMP_RICE,np)
	        annavg_bgnpp(np) = ws*annavg_bgnpp_component(METHANE_COMP_SOIL,np) + &
	           r*annavg_bgnpp_component(METHANE_COMP_RICE,np)
	        annavg_somhr(np) = ws*annavg_somhr_component(METHANE_COMP_SOIL,np) + &
	           r*annavg_somhr_component(METHANE_COMP_RICE,np)
	        annavg_finrw(np) = ws*annavg_finrw_component(METHANE_COMP_SOIL,np) + &
	           r*annavg_finrw_component(METHANE_COMP_RICE,np)
	        tempavg_agnpp(np) = ws*tempavg_agnpp_component(METHANE_COMP_SOIL,np) + &
	           r*tempavg_agnpp_component(METHANE_COMP_RICE,np)
	        tempavg_bgnpp(np) = ws*tempavg_bgnpp_component(METHANE_COMP_SOIL,np) + &
	           r*tempavg_bgnpp_component(METHANE_COMP_RICE,np)
	        annsum_counter(np) = ws*annsum_counter_component(METHANE_COMP_SOIL,np) + &
	           r*annsum_counter_component(METHANE_COMP_RICE,np)
	        tempavg_somhr(np) = ws*tempavg_somhr_component(METHANE_COMP_SOIL,np) + &
	           r*tempavg_somhr_component(METHANE_COMP_RICE,np)
	        tempavg_finrw(np) = ws*tempavg_finrw_component(METHANE_COMP_SOIL,np) + &
	           r*tempavg_finrw_component(METHANE_COMP_RICE,np)
	        finundated_lag(np) = ws*finundated_lag_component(METHANE_COMP_SOIL,np) + &
	           r*finundated_lag_component(METHANE_COMP_RICE,np)

	        IF (np <= size(patchtype) .and. patchtype(np) /= 0) THEN
	           rice_fraction_prev(np) = 0._r8
	           conc_o2_unsat_component(:,:,np) = spread(conc_o2_unsat(:,np), 2, N_METHANE_COMP)
	           conc_o2_sat_component(:,:,np) = spread(conc_o2_sat(:,np), 2, N_METHANE_COMP)
	           conc_methane_unsat_component(:,:,np) = spread(conc_methane_unsat(:,np), 2, N_METHANE_COMP)
	           conc_methane_sat_component(:,:,np) = spread(conc_methane_sat(:,np), 2, N_METHANE_COMP)
	           layer_sat_lag_component(:,:,np) = spread(layer_sat_lag(:,np), 2, N_METHANE_COMP)
	           annavg_agnpp_component(:,np) = annavg_agnpp(np)
	           annavg_bgnpp_component(:,np) = annavg_bgnpp(np)
	           annavg_somhr_component(:,np) = annavg_somhr(np)
	           annavg_finrw_component(:,np) = annavg_finrw(np)
	           tempavg_agnpp_component(:,np) = tempavg_agnpp(np)
	           tempavg_bgnpp_component(:,np) = tempavg_bgnpp(np)
	           annsum_counter_component(:,np) = annsum_counter(np)
	           tempavg_somhr_component(:,np) = tempavg_somhr(np)
	           tempavg_finrw_component(:,np) = tempavg_finrw(np)
	           fsat_bef_component(:,np) = fsat_bef(np)
	           finundated_lag_component(:,np) = finundated_lag(np)
	        ENDIF
	     ENDDO
	  END SUBROUTINE sync_component_aggregates_after_lulcc

      SUBROUTINE remap1d(old, new)
         real(r8), intent(in) :: old(:)
         real(r8), intent(inout) :: new(:)
         integer :: np, op, src, link
         real(r8) :: w, wsum, val
         DO np = 1, min(size(new), nnew)
            wsum = 0._r8
            val = 0._r8
            IF (present(lccpct_patches)) THEN
	               DO link = map_start(np), map_start(np+1)-1
	                  op = map_old(link)
	                  IF (op > size(old)) CYCLE
	                  w = map_source_weight(link)
	                  IF (w <= 0._r8) CYCLE
                  val = val + w * old(op)
                  wsum = wsum + w
               ENDDO
            ENDIF
            IF (wsum > 0._r8) THEN
               new(np) = val / wsum
            ELSE
               src = fallback_source(np, size(old))
               IF (src > 0) new(np) = old(src)
            ENDIF
         ENDDO
      END SUBROUTINE remap1d

      SUBROUTINE remap1d_mass(old, new)
         real(r8), intent(in) :: old(:)
         real(r8), intent(inout) :: new(:)
         integer :: np, op, src, link
         real(r8) :: w, wsum, val, denom
         logical :: conserve_area_mass

         DO np = 1, min(size(new), nnew)
            wsum = 0._r8
            val = 0._r8
	         conserve_area_mass = map_mass_available(np)
            IF (present(lccpct_patches)) THEN
	               DO link = map_start(np), map_start(np+1)-1
	                  op = map_old(link)
	                  IF (op > size(old)) CYCLE
                  IF (conserve_area_mass) THEN
	                     w = map_mass_weight(link)
                  ELSE
	                     w = map_source_weight(link)
                  ENDIF
                  IF (w <= 0._r8) CYCLE
                  val = val + w * old(op)
                  wsum = wsum + w
               ENDDO
            ENDIF
            IF (wsum > 0._r8) THEN
               denom = remap_denominator(np, wsum, conserve_area_mass)
               new(np) = val / denom
            ELSE
               src = fallback_source(np, size(old))
               IF (src > 0) new(np) = old(src)
            ENDIF
         ENDDO
      END SUBROUTINE remap1d_mass

	      SUBROUTINE remap2d(old, new)
         real(r8), intent(in) :: old(:,:)
         real(r8), intent(inout) :: new(:,:)
         integer :: np, op, src, link
         real(r8) :: w, wsum
         real(r8) :: default_vals(size(new,1))
         DO np = 1, min(size(new,2), nnew)
            wsum = 0._r8
            default_vals(:) = new(:,np)
            IF (present(lccpct_patches)) THEN
               new(:,np) = 0._r8
	               DO link = map_start(np), map_start(np+1)-1
	                  op = map_old(link)
	                  IF (op > size(old,2)) CYCLE
	                  w = map_source_weight(link)
	                  IF (w <= 0._r8) CYCLE
                  new(1:min(size(new,1),size(old,1)),np) = &
                     new(1:min(size(new,1),size(old,1)),np) + &
                     w * old(1:min(size(new,1),size(old,1)),op)
                  wsum = wsum + w
               ENDDO
            ENDIF
            IF (wsum > 0._r8) THEN
               new(:,np) = new(:,np) / wsum
            ELSE
               src = fallback_source(np, size(old,2))
               IF (src > 0) THEN
                  new(:,np) = 0._r8
                  new(1:min(size(new,1),size(old,1)),np) = old(1:min(size(new,1),size(old,1)),src)
               ELSE
                  new(:,np) = default_vals(:)
               ENDIF
            ENDIF
         ENDDO
	      END SUBROUTINE remap2d


         SUBROUTINE repartition_ch4_totcol_after_lulcc()
            integer :: np, n
            real(r8) :: total, water_total, land_eff, f, scale

            n = min(size(totcol_methane), size(totcol_methane_unsat), &
               size(totcol_methane_sat), size(totcol_methane_lake), &
               size(lake_water_ch4_stock), size(lake_water_o2_stock), &
               size(lake_frozen_ch4_stock), size(lake_frozen_o2_stock), &
               size(lake_liquid_fraction_prev), &
               size(patchclass_new), nnew)
            DO np = 1, n
               total = max(totcol_methane(np), 0._r8)
               totcol_methane(np) = total

               IF (patchclass_new(np) == WATERBODY) THEN
	                  ! The generic column is the canonical total inventory;
	                  ! keep the independent water stock once and assign only
	                  ! the remainder to the lake sediment representation.
	                  lake_water_ch4_stock(np) = max(lake_water_ch4_stock(np), 0._r8)
	                  lake_frozen_ch4_stock(np) = max(lake_frozen_ch4_stock(np), 0._r8)
	                  water_total = lake_water_ch4_stock(np) + lake_frozen_ch4_stock(np)
	                  IF (water_total > total .and. water_total > tiny(1._r8)) THEN
	                     scale = total / water_total
	                     lake_water_ch4_stock(np) = lake_water_ch4_stock(np) * scale
	                     lake_frozen_ch4_stock(np) = lake_frozen_ch4_stock(np) * scale
	                  ENDIF
	                  lake_water_o2_stock(np) = max(lake_water_o2_stock(np), 0._r8)
	                  lake_frozen_o2_stock(np) = max(lake_frozen_o2_stock(np), 0._r8)
	                  totcol_methane_lake(np) = total - lake_water_ch4_stock(np) - &
	                     lake_frozen_ch4_stock(np)
	                  totcol_methane_sat(np) = totcol_methane_lake(np)
                  totcol_methane_unsat(np) = 0._r8
               ELSE
	                  ! Any remapped lake-water CH4 is already present in the
	                  ! canonical total and is repartitioned into land branches
	                  ! below; clear the lake-only views without losing mass.
	                  lake_water_ch4_stock(np) = 0._r8
	                  lake_water_o2_stock(np) = 0._r8
	                  lake_frozen_ch4_stock(np) = 0._r8
	                  lake_frozen_o2_stock(np) = 0._r8
	                  lake_liquid_fraction_prev(np) = spval
                  totcol_methane_lake(np) = 0._r8
                  totcol_methane_sat(np) = max(totcol_methane_sat(np), 0._r8)
                  totcol_methane_unsat(np) = max(totcol_methane_unsat(np), 0._r8)
                  f = methane_lulcc_saturated_fraction(np)
                  land_eff = f * totcol_methane_sat(np) + &
                     (1._r8 - f) * totcol_methane_unsat(np)
                  IF (total <= tiny(1._r8)) THEN
                     totcol_methane_sat(np) = 0._r8
                     totcol_methane_unsat(np) = 0._r8
                  ELSEIF (land_eff > tiny(1._r8)) THEN
                     scale = total / land_eff
                     totcol_methane_sat(np) = totcol_methane_sat(np) * scale
                     totcol_methane_unsat(np) = totcol_methane_unsat(np) * scale
                  ELSE
                     ! ponytail: class-changed patches have no branch history;
                     ! seed both land branches equally so any finundated gives
                     ! the conserved remapped CH4 column without inventing a
                     ! saturated/unsaturated contrast.
                     totcol_methane_sat(np) = total
                     totcol_methane_unsat(np) = total
                  ENDIF
               ENDIF
            ENDDO
         END SUBROUTINE repartition_ch4_totcol_after_lulcc

         REAL(r8) FUNCTION methane_lulcc_saturated_fraction(np) RESULT(f)
            integer, intent(in) :: np

            f = 0._r8
            IF (np <= size(fsat_bef)) f = fsat_bef(np)
            IF (f < 0._r8 .or. f > 1._r8 .or. ieee_is_nan(f) .or. &
                abs(f) >= 0.5_r8 * abs(spval)) f = 0._r8
            f = min(max(f, 0._r8), 1._r8)
         END FUNCTION methane_lulcc_saturated_fraction

      SUBROUTINE sync_ch4_conc_to_totcol(conc, totcol)
	         real(r8), intent(inout) :: conc(:,:)
	         real(r8), intent(inout) :: totcol(:)
	         integer :: j, np, nlev
	         real(r8) :: col, target, dz, dzsum

	         nlev = min(size(conc,1), nl_soil)
	         dzsum = 0._r8
	         DO j = 1, nlev
	            dzsum = dzsum + max(dz_soi(j), 0._r8)
	         ENDDO
	         IF (dzsum <= tiny(1._r8)) RETURN

	         DO np = 1, min(size(conc,2), size(totcol), nnew)
	            target = max(totcol(np), 0._r8)
	            totcol(np) = target
	            col = 0._r8
	            DO j = 1, nlev
	               dz = max(dz_soi(j), 0._r8)
	               conc(j,np) = max(conc(j,np), 0._r8)
	               col = col + conc(j,np) * dz
	            ENDDO

	            IF (target <= tiny(1._r8)) THEN
	               conc(1:nlev,np) = 0._r8
	            ELSEIF (col > tiny(1._r8)) THEN
	               conc(1:nlev,np) = conc(1:nlev,np) * (target / col)
	            ELSE
	               conc(1:nlev,np) = target / dzsum
	            ENDIF
	         ENDDO
	      END SUBROUTINE sync_ch4_conc_to_totcol

		      INTEGER FUNCTION fallback_source(np, old_n) RESULT(src)
         integer, intent(in) :: np, old_n
         src = 0
	     IF (np < 1 .or. np > size(fallback_map)) RETURN
	     IF (fallback_map(np) <= old_n) src = fallback_map(np)
         ! Do not fall back across patch classes.  Methane lake/sediment
         ! inventories are class-specific and can be corrupted by copying from
         ! a same-eindex non-lake patch.
	      END FUNCTION fallback_source

	      REAL(r8) FUNCTION component_base_weight(np, link) RESULT(w)
	         integer, intent(in) :: np, link

	         IF (map_mass_available(np)) THEN
	            w = map_mass_weight(link)
	         ELSE
	            w = map_source_weight(link)
	         ENDIF
	      END FUNCTION component_base_weight

	      REAL(r8) FUNCTION component_transfer_weight(np, link, component) RESULT(w)
	         integer, intent(in) :: np, link, component
	         integer :: op
	         real(r8) :: r

	         w = component_base_weight(np, link)
	         op = map_old(link)
	         IF (w <= 0._r8 .or. op > size(lulcc_rice_fraction_prev_old)) RETURN
	         r = min(max(lulcc_rice_fraction_prev_old(op), 0._r8), 1._r8)
	         IF (component == METHANE_COMP_RICE) THEN
	            w = w*r
	         ELSE
	            w = w*(1._r8-r)
	         ENDIF
	      END FUNCTION component_transfer_weight

	      LOGICAL FUNCTION area_mass_remap_available(np) RESULT(ok)
	         integer, intent(in) :: np

	         ok = present(lccpct_patches) .and. present(old_patch_area) .and. present(new_patch_area)
	         IF (.not. ok) RETURN
	         ok = np <= size(new_patch_area)
	         IF (.not. ok) RETURN
	         ok = new_patch_area(np) > tiny(1._r8)
	      END FUNCTION area_mass_remap_available

	      REAL(r8) FUNCTION remap_denominator(np, wsum, conserve_mass) RESULT(denom)
	         integer, intent(in) :: np
	         real(r8), intent(in) :: wsum
	         logical, intent(in) :: conserve_mass

	         denom = max(wsum, tiny(1._r8))
	         IF (conserve_mass .and. present(new_patch_area)) THEN
	            IF (np <= size(new_patch_area) .and. new_patch_area(np) > tiny(1._r8)) THEN
	               denom = new_patch_area(np)
	            ENDIF
	         ENDIF
	      END FUNCTION remap_denominator
	   END SUBROUTINE remap_methane_lulcc_state


   SUBROUTINE clear_methane_lulcc_snapshot ()
      IF (allocated(lulcc_conc_o2_old)) deallocate(lulcc_conc_o2_old)
      IF (allocated(lulcc_conc_methane_old)) deallocate(lulcc_conc_methane_old)
      IF (allocated(lulcc_co2_decomp_depth_old)) deallocate(lulcc_co2_decomp_depth_old)
      IF (allocated(lulcc_co2_oxid_depth_old)) deallocate(lulcc_co2_oxid_depth_old)
      IF (allocated(lulcc_co2_aere_depth_old)) deallocate(lulcc_co2_aere_depth_old)
      IF (allocated(lulcc_co2_decomp_depth_unsat_old)) deallocate(lulcc_co2_decomp_depth_unsat_old)
      IF (allocated(lulcc_co2_decomp_depth_sat_old)) deallocate(lulcc_co2_decomp_depth_sat_old)
      IF (allocated(lulcc_co2_oxid_depth_unsat_old)) deallocate(lulcc_co2_oxid_depth_unsat_old)
      IF (allocated(lulcc_co2_oxid_depth_sat_old)) deallocate(lulcc_co2_oxid_depth_sat_old)
      IF (allocated(lulcc_co2_aere_depth_unsat_old)) deallocate(lulcc_co2_aere_depth_unsat_old)
      IF (allocated(lulcc_co2_aere_depth_sat_old)) deallocate(lulcc_co2_aere_depth_sat_old)
      IF (allocated(lulcc_co2_decomp_depth_lake_old)) deallocate(lulcc_co2_decomp_depth_lake_old)
      IF (allocated(lulcc_co2_oxid_depth_lake_old)) deallocate(lulcc_co2_oxid_depth_lake_old)
      IF (allocated(lulcc_conc_o2_unsat_old)) deallocate(lulcc_conc_o2_unsat_old)
      IF (allocated(lulcc_conc_o2_sat_old)) deallocate(lulcc_conc_o2_sat_old)
      IF (allocated(lulcc_conc_methane_unsat_old)) deallocate(lulcc_conc_methane_unsat_old)
      IF (allocated(lulcc_conc_methane_sat_old)) deallocate(lulcc_conc_methane_sat_old)
      IF (allocated(lulcc_conc_o2_lake_old)) deallocate(lulcc_conc_o2_lake_old)
      IF (allocated(lulcc_conc_methane_lake_old)) deallocate(lulcc_conc_methane_lake_old)
      IF (allocated(lulcc_totcol_methane_old)) deallocate(lulcc_totcol_methane_old)
      IF (allocated(lulcc_co2_decomp_tot_old)) deallocate(lulcc_co2_decomp_tot_old)
      IF (allocated(lulcc_co2_oxid_tot_old)) deallocate(lulcc_co2_oxid_tot_old)
      IF (allocated(lulcc_co2_aere_tot_old)) deallocate(lulcc_co2_aere_tot_old)
      IF (allocated(lulcc_co2_net_tot_old)) deallocate(lulcc_co2_net_tot_old)
      IF (allocated(lulcc_co2_decomp_tot_unsat_old)) deallocate(lulcc_co2_decomp_tot_unsat_old)
      IF (allocated(lulcc_co2_decomp_tot_sat_old)) deallocate(lulcc_co2_decomp_tot_sat_old)
      IF (allocated(lulcc_co2_oxid_tot_unsat_old)) deallocate(lulcc_co2_oxid_tot_unsat_old)
      IF (allocated(lulcc_co2_oxid_tot_sat_old)) deallocate(lulcc_co2_oxid_tot_sat_old)
      IF (allocated(lulcc_co2_net_tot_unsat_old)) deallocate(lulcc_co2_net_tot_unsat_old)
      IF (allocated(lulcc_co2_net_tot_sat_old)) deallocate(lulcc_co2_net_tot_sat_old)
      IF (allocated(lulcc_co2_decomp_tot_lake_old)) deallocate(lulcc_co2_decomp_tot_lake_old)
      IF (allocated(lulcc_co2_oxid_tot_lake_old)) deallocate(lulcc_co2_oxid_tot_lake_old)
      IF (allocated(lulcc_co2_net_tot_lake_old)) deallocate(lulcc_co2_net_tot_lake_old)
      IF (allocated(lulcc_totcol_methane_unsat_old)) deallocate(lulcc_totcol_methane_unsat_old)
      IF (allocated(lulcc_totcol_methane_sat_old)) deallocate(lulcc_totcol_methane_sat_old)
      IF (allocated(lulcc_totcol_methane_lake_old)) deallocate(lulcc_totcol_methane_lake_old)
      IF (allocated(lulcc_lake_water_ch4_stock_old)) deallocate(lulcc_lake_water_ch4_stock_old)
      IF (allocated(lulcc_lake_water_o2_stock_old)) deallocate(lulcc_lake_water_o2_stock_old)
      IF (allocated(lulcc_lake_frozen_ch4_stock_old)) deallocate(lulcc_lake_frozen_ch4_stock_old)
      IF (allocated(lulcc_lake_frozen_o2_stock_old)) deallocate(lulcc_lake_frozen_o2_stock_old)
      IF (allocated(lulcc_lake_liquid_fraction_prev_old)) deallocate(lulcc_lake_liquid_fraction_prev_old)
      IF (allocated(lulcc_grnd_methane_cond_old)) deallocate(lulcc_grnd_methane_cond_old)
      IF (allocated(lulcc_grnd_methane_cond_unsat_old)) deallocate(lulcc_grnd_methane_cond_unsat_old)
      IF (allocated(lulcc_grnd_methane_cond_sat_old)) deallocate(lulcc_grnd_methane_cond_sat_old)
      IF (allocated(lulcc_grnd_methane_cond_lake_old)) deallocate(lulcc_grnd_methane_cond_lake_old)
      IF (allocated(lulcc_layer_sat_lag_old)) deallocate(lulcc_layer_sat_lag_old)
      IF (allocated(lulcc_lake_soilc_old)) deallocate(lulcc_lake_soilc_old)
      IF (allocated(lulcc_annavg_agnpp_old)) deallocate(lulcc_annavg_agnpp_old)
      IF (allocated(lulcc_annavg_bgnpp_old)) deallocate(lulcc_annavg_bgnpp_old)
      IF (allocated(lulcc_annavg_somhr_old)) deallocate(lulcc_annavg_somhr_old)
      IF (allocated(lulcc_annavg_finrw_old)) deallocate(lulcc_annavg_finrw_old)
      IF (allocated(lulcc_tempavg_agnpp_old)) deallocate(lulcc_tempavg_agnpp_old)
      IF (allocated(lulcc_tempavg_bgnpp_old)) deallocate(lulcc_tempavg_bgnpp_old)
      IF (allocated(lulcc_annsum_counter_old)) deallocate(lulcc_annsum_counter_old)
      IF (allocated(lulcc_tempavg_somhr_old)) deallocate(lulcc_tempavg_somhr_old)
      IF (allocated(lulcc_tempavg_finrw_old)) deallocate(lulcc_tempavg_finrw_old)
      IF (allocated(lulcc_fsat_bef_old)) deallocate(lulcc_fsat_bef_old)
      IF (allocated(lulcc_finundated_lag_old)) deallocate(lulcc_finundated_lag_old)
      IF (allocated(lulcc_methane_dfsat_tot_old)) deallocate(lulcc_methane_dfsat_tot_old)
      IF (allocated(lulcc_f_h2osfc_old)) deallocate(lulcc_f_h2osfc_old)
      IF (allocated(lulcc_forc_pmethanem_old)) deallocate(lulcc_forc_pmethanem_old)
	   IF (allocated(lulcc_f_inund_levee_patch_old)) deallocate(lulcc_f_inund_levee_patch_old)
	   IF (allocated(lulcc_f_inund_flood_patch_old)) deallocate(lulcc_f_inund_flood_patch_old)
	   IF (allocated(lulcc_f_inund_flood_depth_patch_old)) deallocate(lulcc_f_inund_flood_depth_patch_old)
      IF (allocated(lulcc_c_atm_old)) deallocate(lulcc_c_atm_old)
	  IF (allocated(lulcc_conc_o2_unsat_component_old)) deallocate(lulcc_conc_o2_unsat_component_old)
	  IF (allocated(lulcc_conc_o2_sat_component_old)) deallocate(lulcc_conc_o2_sat_component_old)
	  IF (allocated(lulcc_conc_methane_unsat_component_old)) deallocate(lulcc_conc_methane_unsat_component_old)
	  IF (allocated(lulcc_conc_methane_sat_component_old)) deallocate(lulcc_conc_methane_sat_component_old)
	  IF (allocated(lulcc_layer_sat_lag_component_old)) deallocate(lulcc_layer_sat_lag_component_old)
	  IF (allocated(lulcc_annavg_agnpp_component_old)) deallocate(lulcc_annavg_agnpp_component_old)
	  IF (allocated(lulcc_annavg_bgnpp_component_old)) deallocate(lulcc_annavg_bgnpp_component_old)
	  IF (allocated(lulcc_annavg_somhr_component_old)) deallocate(lulcc_annavg_somhr_component_old)
	  IF (allocated(lulcc_annavg_finrw_component_old)) deallocate(lulcc_annavg_finrw_component_old)
	  IF (allocated(lulcc_tempavg_agnpp_component_old)) deallocate(lulcc_tempavg_agnpp_component_old)
	  IF (allocated(lulcc_tempavg_bgnpp_component_old)) deallocate(lulcc_tempavg_bgnpp_component_old)
	  IF (allocated(lulcc_annsum_counter_component_old)) deallocate(lulcc_annsum_counter_component_old)
	  IF (allocated(lulcc_tempavg_somhr_component_old)) deallocate(lulcc_tempavg_somhr_component_old)
	  IF (allocated(lulcc_tempavg_finrw_component_old)) deallocate(lulcc_tempavg_finrw_component_old)
	  IF (allocated(lulcc_fsat_bef_component_old)) deallocate(lulcc_fsat_bef_component_old)
	  IF (allocated(lulcc_finundated_lag_component_old)) deallocate(lulcc_finundated_lag_component_old)
	  IF (allocated(lulcc_rice_fraction_prev_old)) deallocate(lulcc_rice_fraction_prev_old)
      methane_lulcc_snapshot_valid = .false.
   END SUBROUTINE clear_methane_lulcc_snapshot


   SUBROUTINE initialize_methane_lake_soilc_from_surface (patchtype_in, lake_soilc_srf_in, allowlakeprod, &
      initialize_patch)
      USE MOD_SPMD_Task, only: CoLM_stop
      integer,  intent(in) :: patchtype_in(:)
      real(r8), intent(in), optional :: lake_soilc_srf_in(:,:)
      logical,  intent(in) :: allowlakeprod
      logical,  intent(in), optional :: initialize_patch(:)

      integer :: ipatch, npatch
      logical, allocatable :: needs_initialization(:)
      integer :: lake_soilc_nlake, lake_soilc_missing
      real(r8), parameter :: smallnumber = 1.e-12_r8

      IF (.not. allocated(lake_soilc)) RETURN
      IF (.not. allowlakeprod) RETURN

      lake_soilc_nlake = count(patchtype_in == PATCHTYPE_LAKE)
      IF (lake_soilc_nlake == 0) RETURN

      npatch = size(lake_soilc,2)
      IF (size(patchtype_in) /= npatch) THEN
         CALL CoLM_stop (' ***** ERROR: lake CH4 surface carbon dimensions do not match the active patch layout.')
      ENDIF
      IF (present(initialize_patch)) THEN
         IF (size(initialize_patch) /= npatch) THEN
            CALL CoLM_stop (' ***** ERROR: lake CH4 LULCC initialization mask has the wrong patch dimension.')
         ENDIF
      ENDIF

      ! Only cold/new lakes need surface carbon. Existing lakes may have
      ! legitimately exhausted their inventory and must never be refilled.
      allocate(needs_initialization(npatch))
      needs_initialization = patchtype_in == PATCHTYPE_LAKE
      IF (present(initialize_patch)) THEN
         needs_initialization = needs_initialization .and. initialize_patch
      ELSE
         DO ipatch = 1, npatch
            IF (sum(max(lake_soilc(:,ipatch), 0._r8)) > smallnumber) &
               needs_initialization(ipatch) = .false.
         ENDDO
      ENDIF
      IF (.not. any(needs_initialization)) RETURN
      IF (.not. present(lake_soilc_srf_in)) THEN
         CALL CoLM_stop (' ***** ERROR: lake CH4 initialization requires lake_soilc surface data.')
      ENDIF
      IF (size(lake_soilc_srf_in,1) < nl_soil .or. size(lake_soilc_srf_in,2) /= npatch) THEN
         CALL CoLM_stop (' ***** ERROR: lake CH4 surface carbon dimensions do not match the active patch layout.')
      ENDIF

      lake_soilc_missing = 0
      DO ipatch = 1, npatch
         IF (.not. needs_initialization(ipatch)) CYCLE
         IF (any(invalid_restart_value(lake_soilc_srf_in(1:nl_soil,ipatch))) .or. &
             any(lake_soilc_srf_in(1:nl_soil,ipatch) < 0._r8) .or. &
             sum(lake_soilc_srf_in(1:nl_soil,ipatch)) <= smallnumber) THEN
            lake_soilc_missing = lake_soilc_missing + 1
         ENDIF
      END DO

      IF (lake_soilc_missing > 0) THEN
         write(6,*) ' ERROR: lake CH4 production requires positive finite lake_soilc for every initialized lake patch; missing ', &
            lake_soilc_missing, ' of ', lake_soilc_nlake, ' local lake patches.'
         CALL CoLM_stop (' ***** ERROR: incomplete lake_soilc input while lake CH4 production is enabled.')
      ENDIF

      DO ipatch = 1, npatch
         IF (.not. needs_initialization(ipatch)) CYCLE
         lake_soilc(:,ipatch) = max(lake_soilc_srf_in(1:nl_soil,ipatch), 0._r8)
      END DO
   END SUBROUTINE initialize_methane_lake_soilc_from_surface


   ELEMENTAL LOGICAL FUNCTION invalid_restart_value (x)
      real(r8), intent(in) :: x

      invalid_restart_value = ieee_is_nan(x) .or. (abs(x) >= 0.5_r8 * abs(spval))
   END FUNCTION invalid_restart_value


   ELEMENTAL LOGICAL FUNCTION invalid_restart_fraction_or_sentinel (x)
      real(r8), intent(in) :: x

      IF (x == spval) THEN
         invalid_restart_fraction_or_sentinel = .false.
      ELSE
         invalid_restart_fraction_or_sentinel = invalid_restart_value(x) .or. &
            x < 0._r8 .or. x > 1._r8
      ENDIF
   END FUNCTION invalid_restart_fraction_or_sentinel

END MODULE MOD_MT_Methane_State
#endif
