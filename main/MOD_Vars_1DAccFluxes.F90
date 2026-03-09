#include <define.h>

MODULE MOD_Vars_1DAccFluxes

   USE MOD_Precision
#ifdef DataAssimilation
   USE MOD_DA_Vars_TimeVariables
   USE MOD_DA_Vars_1DFluxes
   USE MOD_Namelist
#endif
#ifdef EXTERNAL_LAKE
   USE MOD_Lake_1DAccVars
#endif

   real(r8) :: nac ! number of accumulation
   real(r8), allocatable :: nac_ln      (:)
   real(r8), allocatable :: nac_dt      (:)
   logical,  allocatable :: filter_dt   (:)

   real(r8), allocatable :: a_us        (:)
   real(r8), allocatable :: a_vs        (:)
   real(r8), allocatable :: a_t         (:)
   real(r8), allocatable :: a_q         (:)
   real(r8), allocatable :: a_prc       (:)
   real(r8), allocatable :: a_prl       (:)
   real(r8), allocatable :: a_pbot      (:)
   real(r8), allocatable :: a_frl       (:)
   real(r8), allocatable :: a_solarin   (:)
   real(r8), allocatable :: a_hpbl      (:)

   real(r8), allocatable :: a_taux      (:)
   real(r8), allocatable :: a_tauy      (:)
   real(r8), allocatable :: a_fsena     (:)
   real(r8), allocatable :: a_lfevpa    (:)
   real(r8), allocatable :: a_fevpa     (:)
   real(r8), allocatable :: a_fsenl     (:)
   real(r8), allocatable :: a_fevpl     (:)
   real(r8), allocatable :: a_etr       (:)
   real(r8), allocatable :: a_fseng     (:)
   real(r8), allocatable :: a_fevpg     (:)
   real(r8), allocatable :: a_fgrnd     (:)
   real(r8), allocatable :: a_sabvsun   (:)
   real(r8), allocatable :: a_sabvsha   (:)
   real(r8), allocatable :: a_sabg      (:)
   real(r8), allocatable :: a_olrg      (:)
   real(r8), allocatable :: a_rnet      (:)
   real(r8), allocatable :: a_xerr      (:)
   real(r8), allocatable :: a_zerr      (:)
   real(r8), allocatable :: a_rsur      (:)
   real(r8), allocatable :: a_rsur_se   (:)
   real(r8), allocatable :: a_rsur_ie   (:)
   real(r8), allocatable :: a_rsub      (:)
   real(r8), allocatable :: a_rnof      (:)
#ifdef CatchLateralFlow
   real(r8), allocatable :: a_xwsur     (:)
   real(r8), allocatable :: a_xwsub     (:)
   real(r8), allocatable :: a_fldarea   (:)
#endif
   real(r8), allocatable :: a_qintr     (:)
   real(r8), allocatable :: a_qinfl     (:)
   real(r8), allocatable :: a_qdrip     (:)
   real(r8), allocatable :: a_rstfacsun (:)
   real(r8), allocatable :: a_rstfacsha (:)
   real(r8), allocatable :: a_gssun     (:)
   real(r8), allocatable :: a_gssha     (:)
   real(r8), allocatable :: a_rss       (:)
   real(r8), allocatable :: a_wdsrf     (:)
   real(r8), allocatable :: a_zwt       (:)
   real(r8), allocatable :: a_wa        (:)
   real(r8), allocatable :: a_wat       (:)
   real(r8), allocatable :: a_wetwat    (:)
   real(r8), allocatable :: a_assim     (:)
   real(r8), allocatable :: a_respc     (:)
   real(r8), allocatable :: a_assimsun  (:)
   real(r8), allocatable :: a_assimsha  (:)
   real(r8), allocatable :: a_etrsun    (:)
   real(r8), allocatable :: a_etrsha    (:)

   real(r8), allocatable :: a_qcharge   (:)

   real(r8), allocatable :: a_t_grnd    (:)
   real(r8), allocatable :: a_tleaf     (:)
   real(r8), allocatable :: a_ldew      (:)
   real(r8), allocatable :: a_ldew_rain (:)
   real(r8), allocatable :: a_ldew_snow (:)
   real(r8), allocatable :: a_scv       (:)
   real(r8), allocatable :: a_snowdp    (:)
   real(r8), allocatable :: a_fsno      (:)
   real(r8), allocatable :: a_frcsat    (:)
   real(r8), allocatable :: a_sigf      (:)
   real(r8), allocatable :: a_green     (:)
   real(r8), allocatable :: a_lai       (:)
   real(r8), allocatable :: a_laisun    (:)
   real(r8), allocatable :: a_laisha    (:)
   real(r8), allocatable :: a_sai       (:)
   real(r8), allocatable :: a_alb       (:,:,:)

#ifdef HYPERSPECTRAL
   real(r8), allocatable :: a_alb_hires (:,:,:)
   real(r8), allocatable :: a_reflectance_out   (:,:,:)
   real(r8), allocatable :: a_transmittance_out (:,:,:)
#endif
   real(r8), allocatable :: a_emis      (:)
   real(r8), allocatable :: a_z0m       (:)
   real(r8), allocatable :: a_trad      (:)
   real(r8), allocatable :: a_tref      (:)
   real(r8), allocatable :: a_t2m_wmo   (:)
   real(r8), allocatable :: a_qref      (:)
   real(r8), allocatable :: a_rain      (:)
   real(r8), allocatable :: a_snow      (:)

   real(r8), allocatable :: a_o3uptakesun(:)
   real(r8), allocatable :: a_o3uptakesha(:)

#ifdef DataAssimilation
   real(r8), allocatable :: a_h2osoi_ens     (:,:,:)
   real(r8), allocatable :: a_t_brt_smap_ens (:,:,:)
   real(r8), allocatable :: a_t_brt_fy3d_ens (:,:,:)
   real(r8), allocatable :: a_t_brt_smap       (:,:)
   real(r8), allocatable :: a_t_brt_fy3d       (:,:)
   real(r8), allocatable :: a_wliq_soisno_ens(:,:,:)
   real(r8), allocatable :: a_wice_soisno_ens(:,:,:)
   real(r8), allocatable :: a_t_soisno_ens   (:,:,:)
#endif

#ifdef URBAN_MODEL
   real(r8), allocatable :: a_t_room    (:) !temperature of inner building [K]
   real(r8), allocatable :: a_tafu      (:) !temperature of outer building [K]
   real(r8), allocatable :: a_fhac      (:) !sensible flux from heat or cool AC [W/m2]
   real(r8), allocatable :: a_fwst      (:) !waste heat flux from heat or cool AC [W/m2]
   real(r8), allocatable :: a_fach      (:) !flux from inner and outer air exchange [W/m2]
   real(r8), allocatable :: a_fahe      (:) !flux from metabolic and vehicle [W/m2]
   real(r8), allocatable :: a_fhah      (:) !sensible flux from heating [W/m2]
   real(r8), allocatable :: a_vehc      (:) !flux from vehicle [W/m2]
   real(r8), allocatable :: a_meta      (:) !flux from metabolic [W/m2]

   real(r8), allocatable :: a_senroof   (:) !sensible heat flux from roof [W/m2]
   real(r8), allocatable :: a_senwsun   (:) !sensible heat flux from sunlit wall [W/m2]
   real(r8), allocatable :: a_senwsha   (:) !sensible heat flux from shaded wall [W/m2]
   real(r8), allocatable :: a_sengimp   (:) !sensible heat flux from impervious road [W/m2]
   real(r8), allocatable :: a_sengper   (:) !sensible heat flux from pervious road [W/m2]
   real(r8), allocatable :: a_senurbl   (:) !sensible heat flux from urban vegetation [W/m2]

   real(r8), allocatable :: a_lfevproof (:) !latent heat flux from roof [W/m2]
   real(r8), allocatable :: a_lfevpgimp (:) !latent heat flux from impervious road [W/m2]
   real(r8), allocatable :: a_lfevpgper (:) !latent heat flux from pervious road [W/m2]
   real(r8), allocatable :: a_lfevpurbl (:) !latent heat flux from urban vegetation [W/m2]

   real(r8), allocatable :: a_troof     (:) !temperature of roof [K]
   real(r8), allocatable :: a_twall     (:) !temperature of wall [K]
#endif

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   real(r8), allocatable :: a_lai_enftemp      (:) !1
   real(r8), allocatable :: a_lai_enfboreal    (:) !2
   real(r8), allocatable :: a_lai_dnfboreal    (:) !3
   real(r8), allocatable :: a_lai_ebftrop      (:) !4
   real(r8), allocatable :: a_lai_ebftemp      (:) !5
   real(r8), allocatable :: a_lai_dbftrop      (:) !6
   real(r8), allocatable :: a_lai_dbftemp      (:) !7
   real(r8), allocatable :: a_lai_dbfboreal    (:) !8
   real(r8), allocatable :: a_lai_ebstemp      (:) !9
   real(r8), allocatable :: a_lai_dbstemp      (:) !10
   real(r8), allocatable :: a_lai_dbsboreal    (:) !11
   real(r8), allocatable :: a_lai_c3arcgrass   (:) !12
   real(r8), allocatable :: a_lai_c3grass      (:) !13
   real(r8), allocatable :: a_lai_c4grass      (:) !14
#endif

#ifdef BGC
   real(r8), allocatable :: a_leafc              (:)
   real(r8), allocatable :: a_leafc_storage      (:)
   real(r8), allocatable :: a_leafc_xfer         (:)
   real(r8), allocatable :: a_frootc             (:)
   real(r8), allocatable :: a_frootc_storage     (:)
   real(r8), allocatable :: a_frootc_xfer        (:)
   real(r8), allocatable :: a_livestemc          (:)
   real(r8), allocatable :: a_livestemc_storage  (:)
   real(r8), allocatable :: a_livestemc_xfer     (:)
   real(r8), allocatable :: a_deadstemc          (:)
   real(r8), allocatable :: a_deadstemc_storage  (:)
   real(r8), allocatable :: a_deadstemc_xfer     (:)
   real(r8), allocatable :: a_livecrootc         (:)
   real(r8), allocatable :: a_livecrootc_storage (:)
   real(r8), allocatable :: a_livecrootc_xfer    (:)
   real(r8), allocatable :: a_deadcrootc         (:)
   real(r8), allocatable :: a_deadcrootc_storage (:)
   real(r8), allocatable :: a_deadcrootc_xfer    (:)
   real(r8), allocatable :: a_grainc             (:)
   real(r8), allocatable :: a_grainc_storage     (:)
   real(r8), allocatable :: a_grainc_xfer        (:)
   real(r8), allocatable :: a_leafn              (:)
   real(r8), allocatable :: a_leafn_storage      (:)
   real(r8), allocatable :: a_leafn_xfer         (:)
   real(r8), allocatable :: a_frootn             (:)
   real(r8), allocatable :: a_frootn_storage     (:)
   real(r8), allocatable :: a_frootn_xfer        (:)
   real(r8), allocatable :: a_livestemn          (:)
   real(r8), allocatable :: a_livestemn_storage  (:)
   real(r8), allocatable :: a_livestemn_xfer     (:)
   real(r8), allocatable :: a_deadstemn          (:)
   real(r8), allocatable :: a_deadstemn_storage  (:)
   real(r8), allocatable :: a_deadstemn_xfer     (:)
   real(r8), allocatable :: a_livecrootn         (:)
   real(r8), allocatable :: a_livecrootn_storage (:)
   real(r8), allocatable :: a_livecrootn_xfer    (:)
   real(r8), allocatable :: a_deadcrootn         (:)
   real(r8), allocatable :: a_deadcrootn_storage (:)
   real(r8), allocatable :: a_deadcrootn_xfer    (:)
   real(r8), allocatable :: a_grainn             (:)
   real(r8), allocatable :: a_grainn_storage     (:)
   real(r8), allocatable :: a_grainn_xfer        (:)
   real(r8), allocatable :: a_retransn           (:)
   real(r8), allocatable :: a_gpp                (:)
   real(r8), allocatable :: a_downreg            (:)
   real(r8), allocatable :: a_ar                 (:)
   real(r8), allocatable :: a_cwdprod            (:)
   real(r8), allocatable :: a_cwddecomp          (:)
   real(r8), allocatable :: a_hr                 (:)
   real(r8), allocatable :: a_fpg                (:)
   real(r8), allocatable :: a_fpi                (:)
   real(r8), allocatable :: a_totvegc            (:)
   real(r8), allocatable :: a_totlitc            (:)
   real(r8), allocatable :: a_totcwdc            (:)
   real(r8), allocatable :: a_totsomc            (:)
   real(r8), allocatable :: a_totcolc            (:)
   real(r8), allocatable :: a_totvegn            (:)
   real(r8), allocatable :: a_totlitn            (:)
   real(r8), allocatable :: a_totcwdn            (:)
   real(r8), allocatable :: a_totsomn            (:)
   real(r8), allocatable :: a_totcoln            (:)
   real(r8), allocatable :: a_gpp_enftemp        (:) !1
   real(r8), allocatable :: a_gpp_enfboreal      (:) !2
   real(r8), allocatable :: a_gpp_dnfboreal      (:) !3
   real(r8), allocatable :: a_gpp_ebftrop        (:) !4
   real(r8), allocatable :: a_gpp_ebftemp        (:) !5
   real(r8), allocatable :: a_gpp_dbftrop        (:) !6
   real(r8), allocatable :: a_gpp_dbftemp        (:) !7
   real(r8), allocatable :: a_gpp_dbfboreal      (:) !8
   real(r8), allocatable :: a_gpp_ebstemp        (:) !9
   real(r8), allocatable :: a_gpp_dbstemp        (:) !10
   real(r8), allocatable :: a_gpp_dbsboreal      (:) !11
   real(r8), allocatable :: a_gpp_c3arcgrass     (:) !12
   real(r8), allocatable :: a_gpp_c3grass        (:) !13
   real(r8), allocatable :: a_gpp_c4grass        (:) !14
   real(r8), allocatable :: a_npp_enftemp      (:) !1
   real(r8), allocatable :: a_npp_enfboreal    (:) !2
   real(r8), allocatable :: a_npp_dnfboreal    (:) !3
   real(r8), allocatable :: a_npp_ebftrop      (:) !4
   real(r8), allocatable :: a_npp_ebftemp      (:) !5
   real(r8), allocatable :: a_npp_dbftrop      (:) !6
   real(r8), allocatable :: a_npp_dbftemp      (:) !7
   real(r8), allocatable :: a_npp_dbfboreal    (:) !8
   real(r8), allocatable :: a_npp_ebstemp      (:) !9
   real(r8), allocatable :: a_npp_dbstemp      (:) !10
   real(r8), allocatable :: a_npp_dbsboreal    (:) !11
   real(r8), allocatable :: a_npp_c3arcgrass   (:) !12
   real(r8), allocatable :: a_npp_c3grass      (:) !13
   real(r8), allocatable :: a_npp_c4grass      (:) !14
   real(r8), allocatable :: a_npptoleafc_enftemp      (:) !1
   real(r8), allocatable :: a_npptoleafc_enfboreal    (:) !2
   real(r8), allocatable :: a_npptoleafc_dnfboreal    (:) !3
   real(r8), allocatable :: a_npptoleafc_ebftrop      (:) !4
   real(r8), allocatable :: a_npptoleafc_ebftemp      (:) !5
   real(r8), allocatable :: a_npptoleafc_dbftrop      (:) !6
   real(r8), allocatable :: a_npptoleafc_dbftemp      (:) !7
   real(r8), allocatable :: a_npptoleafc_dbfboreal    (:) !8
   real(r8), allocatable :: a_npptoleafc_ebstemp      (:) !9
   real(r8), allocatable :: a_npptoleafc_dbstemp      (:) !10
   real(r8), allocatable :: a_npptoleafc_dbsboreal    (:) !11
   real(r8), allocatable :: a_npptoleafc_c3arcgrass   (:) !12
   real(r8), allocatable :: a_npptoleafc_c3grass      (:) !13
   real(r8), allocatable :: a_npptoleafc_c4grass      (:) !14
   real(r8), allocatable :: a_leafc_enftemp      (:) !1
   real(r8), allocatable :: a_leafc_enfboreal    (:) !2
   real(r8), allocatable :: a_leafc_dnfboreal    (:) !3
   real(r8), allocatable :: a_leafc_ebftrop      (:) !4
   real(r8), allocatable :: a_leafc_ebftemp      (:) !5
   real(r8), allocatable :: a_leafc_dbftrop      (:) !6
   real(r8), allocatable :: a_leafc_dbftemp      (:) !7
   real(r8), allocatable :: a_leafc_dbfboreal    (:) !8
   real(r8), allocatable :: a_leafc_ebstemp      (:) !9
   real(r8), allocatable :: a_leafc_dbstemp      (:) !10
   real(r8), allocatable :: a_leafc_dbsboreal    (:) !11
   real(r8), allocatable :: a_leafc_c3arcgrass   (:) !12
   real(r8), allocatable :: a_leafc_c3grass      (:) !13
   real(r8), allocatable :: a_leafc_c4grass      (:) !14
   real(r8), allocatable :: a_O2_DECOMP_DEPTH_UNSAT (:,:)
   real(r8), allocatable :: a_CONC_O2_UNSAT         (:,:)
#ifdef CROP
   real(r8), allocatable :: a_pdcorn                (:)
   real(r8), allocatable :: a_pdswheat              (:)
   real(r8), allocatable :: a_pdwwheat              (:)
   real(r8), allocatable :: a_pdsoybean             (:)
   real(r8), allocatable :: a_pdcotton              (:)
   real(r8), allocatable :: a_pdrice1               (:)
   real(r8), allocatable :: a_pdrice2               (:)
   real(r8), allocatable :: a_pdsugarcane           (:)
   real(r8), allocatable :: a_plantdate             (:)
   real(r8), allocatable :: a_manunitro             (:)
   real(r8), allocatable :: a_fertnitro_corn        (:)
   real(r8), allocatable :: a_fertnitro_swheat      (:)
   real(r8), allocatable :: a_fertnitro_wwheat      (:)
   real(r8), allocatable :: a_fertnitro_soybean     (:)
   real(r8), allocatable :: a_fertnitro_cotton      (:)
   real(r8), allocatable :: a_fertnitro_rice1       (:)
   real(r8), allocatable :: a_fertnitro_rice2       (:)
   real(r8), allocatable :: a_fertnitro_sugarcane   (:)
   real(r8), allocatable :: a_irrig_method_corn     (:)
   real(r8), allocatable :: a_irrig_method_swheat   (:)
   real(r8), allocatable :: a_irrig_method_wwheat   (:)
   real(r8), allocatable :: a_irrig_method_soybean  (:)
   real(r8), allocatable :: a_irrig_method_cotton   (:)
   real(r8), allocatable :: a_irrig_method_rice1    (:)
   real(r8), allocatable :: a_irrig_method_rice2    (:)
   real(r8), allocatable :: a_irrig_method_sugarcane(:)

   real(r8), allocatable :: a_cphase                (:)
   real(r8), allocatable :: a_gddplant              (:)
   real(r8), allocatable :: a_gddmaturity           (:)
   real(r8), allocatable :: a_vf                    (:)
   real(r8), allocatable :: a_hui                   (:)
   real(r8), allocatable :: a_cropprod1c            (:)
   real(r8), allocatable :: a_cropprod1c_loss       (:)
   real(r8), allocatable :: a_cropseedc_deficit     (:)
   real(r8), allocatable :: a_grainc_to_cropprodc   (:)
   real(r8), allocatable :: a_grainc_to_seed        (:)
   real(r8), allocatable :: a_fert_to_sminn         (:)

   real(r8), allocatable :: a_sum_irrig          (:)
   real(r8), allocatable :: a_sum_deficit_irrig  (:)
   real(r8), allocatable :: a_sum_irrig_count    (:)
   real(r8), allocatable :: a_waterstorage       (:)
   real(r8), allocatable :: a_groundwater_demand (:)
   real(r8), allocatable :: a_groundwater_supply (:)
   real(r8), allocatable :: a_reservoirriver_demand(:)
   real(r8), allocatable :: a_reservoirriver_supply(:)
   real(r8), allocatable :: a_reservoir_supply     (:)
   real(r8), allocatable :: a_river_supply         (:)
   real(r8), allocatable :: a_runoff_supply        (:)
#endif
   real(r8), allocatable :: a_ndep_to_sminn         (:)
   real(r8), allocatable :: a_abm                   (:)
   real(r8), allocatable :: a_gdp                   (:)
   real(r8), allocatable :: a_peatf                 (:)
   real(r8), allocatable :: a_hdm                   (:)
   real(r8), allocatable :: a_lnfm                  (:)
   real(r8), allocatable :: a_leafcCap              (:)
   real(r8), allocatable :: a_leafc_storageCap      (:)
   real(r8), allocatable :: a_leafc_xferCap         (:)
   real(r8), allocatable :: a_frootcCap             (:)
   real(r8), allocatable :: a_frootc_storageCap     (:)
   real(r8), allocatable :: a_frootc_xferCap        (:)
   real(r8), allocatable :: a_livestemcCap          (:)
   real(r8), allocatable :: a_livestemc_storageCap  (:)
   real(r8), allocatable :: a_livestemc_xferCap     (:)
   real(r8), allocatable :: a_deadstemcCap          (:)
   real(r8), allocatable :: a_deadstemc_storageCap  (:)
   real(r8), allocatable :: a_deadstemc_xferCap     (:)
   real(r8), allocatable :: a_livecrootcCap         (:)
   real(r8), allocatable :: a_livecrootc_storageCap (:)
   real(r8), allocatable :: a_livecrootc_xferCap    (:)
   real(r8), allocatable :: a_deadcrootcCap         (:)
   real(r8), allocatable :: a_deadcrootc_storageCap (:)
   real(r8), allocatable :: a_deadcrootc_xferCap    (:)
   real(r8), allocatable :: a_leafnCap              (:)
   real(r8), allocatable :: a_leafn_storageCap      (:)
   real(r8), allocatable :: a_leafn_xferCap         (:)
   real(r8), allocatable :: a_frootnCap             (:)
   real(r8), allocatable :: a_frootn_storageCap     (:)
   real(r8), allocatable :: a_frootn_xferCap        (:)
   real(r8), allocatable :: a_livestemnCap          (:)
   real(r8), allocatable :: a_livestemn_storageCap  (:)
   real(r8), allocatable :: a_livestemn_xferCap     (:)
   real(r8), allocatable :: a_deadstemnCap          (:)
   real(r8), allocatable :: a_deadstemn_storageCap  (:)
   real(r8), allocatable :: a_deadstemn_xferCap     (:)
   real(r8), allocatable :: a_livecrootnCap         (:)
   real(r8), allocatable :: a_livecrootn_storageCap (:)
   real(r8), allocatable :: a_livecrootn_xferCap    (:)
   real(r8), allocatable :: a_deadcrootnCap         (:)
   real(r8), allocatable :: a_deadcrootn_storageCap (:)
   real(r8), allocatable :: a_deadcrootn_xferCap    (:)
#endif
! Ozone stress variables
   real(r8), allocatable :: a_ozone                 (:)
! End ozone stress variables

   real(r8), allocatable :: a_t_soisno    (:,:)
   real(r8), allocatable :: a_wliq_soisno (:,:)
   real(r8), allocatable :: a_wice_soisno (:,:)
   real(r8), allocatable :: a_h2osoi      (:,:)
   real(r8), allocatable :: a_qlayer      (:,:)
   real(r8), allocatable :: a_lake_deficit  (:)
   real(r8), allocatable :: a_rootr       (:,:)
   real(r8), allocatable :: a_BD_all      (:,:)
   real(r8), allocatable :: a_wfc         (:,:)
   real(r8), allocatable :: a_OM_density  (:,:)
!Plant Hydraulic parameters
   real(r8), allocatable :: a_vegwp       (:,:)
!END plant hydraulic parameters
   real(r8), allocatable :: a_dz_lake     (:,:)
   real(r8), allocatable :: a_t_lake      (:,:)
   real(r8), allocatable :: a_lake_icefrac(:,:)

#ifdef BGC
   real(r8), allocatable :: a_litr1c_vr   (:,:)
   real(r8), allocatable :: a_litr2c_vr   (:,:)
   real(r8), allocatable :: a_litr3c_vr   (:,:)
   real(r8), allocatable :: a_soil1c_vr   (:,:)
   real(r8), allocatable :: a_soil2c_vr   (:,:)
   real(r8), allocatable :: a_soil3c_vr   (:,:)
   real(r8), allocatable :: a_cwdc_vr     (:,:)
   real(r8), allocatable :: a_litr1n_vr   (:,:)
   real(r8), allocatable :: a_litr2n_vr   (:,:)
   real(r8), allocatable :: a_litr3n_vr   (:,:)
   real(r8), allocatable :: a_soil1n_vr   (:,:)
   real(r8), allocatable :: a_soil2n_vr   (:,:)
   real(r8), allocatable :: a_soil3n_vr   (:,:)
   real(r8), allocatable :: a_cwdn_vr     (:,:)
   real(r8), allocatable :: a_totsoiln_vr (:,:)
   real(r8), allocatable :: a_litr1cCap_vr(:,:)
   real(r8), allocatable :: a_litr2cCap_vr(:,:)
   real(r8), allocatable :: a_litr3cCap_vr(:,:)
   real(r8), allocatable :: a_soil1cCap_vr(:,:)
   real(r8), allocatable :: a_soil2cCap_vr(:,:)
   real(r8), allocatable :: a_soil3cCap_vr(:,:)
   real(r8), allocatable :: a_cwdcCap_vr  (:,:)
   real(r8), allocatable :: a_litr1nCap_vr(:,:)
   real(r8), allocatable :: a_litr2nCap_vr(:,:)
   real(r8), allocatable :: a_litr3nCap_vr(:,:)
   real(r8), allocatable :: a_soil1nCap_vr(:,:)
   real(r8), allocatable :: a_soil2nCap_vr(:,:)
   real(r8), allocatable :: a_soil3nCap_vr(:,:)
   real(r8), allocatable :: a_cwdnCap_vr  (:,:)
   real(r8), allocatable :: a_t_scalar    (:,:)
   real(r8), allocatable :: a_w_scalar    (:,:)
   real(r8), allocatable :: a_sminn_vr    (:,:)
   real(r8), allocatable :: decomp_vr_tmp (:,:)
#endif

   real(r8), allocatable :: a_ustar   (:)
   real(r8), allocatable :: a_ustar2  (:)
   real(r8), allocatable :: a_tstar   (:)
   real(r8), allocatable :: a_qstar   (:)
   real(r8), allocatable :: a_zol     (:)
   real(r8), allocatable :: a_rib     (:)
   real(r8), allocatable :: a_fm      (:)
   real(r8), allocatable :: a_fh      (:)
   real(r8), allocatable :: a_fq      (:)

   real(r8), allocatable :: a_us10m   (:)
   real(r8), allocatable :: a_vs10m   (:)
   real(r8), allocatable :: a_fm10m   (:)

   real(r8), allocatable :: a_sr      (:)
   real(r8), allocatable :: a_solvd   (:)
   real(r8), allocatable :: a_solvi   (:)
   real(r8), allocatable :: a_solnd   (:)
   real(r8), allocatable :: a_solni   (:)
   real(r8), allocatable :: a_srvd    (:)
   real(r8), allocatable :: a_srvi    (:)
   real(r8), allocatable :: a_srnd    (:)
   real(r8), allocatable :: a_srni    (:)
   real(r8), allocatable :: a_solvdln (:)
   real(r8), allocatable :: a_solviln (:)
   real(r8), allocatable :: a_solndln (:)
   real(r8), allocatable :: a_solniln (:)
   real(r8), allocatable :: a_srvdln  (:)
   real(r8), allocatable :: a_srviln  (:)
   real(r8), allocatable :: a_srndln  (:)
   real(r8), allocatable :: a_srniln  (:)
#ifdef HYPERSPECTRAL
   real(r8), allocatable :: a_sol_dir_ln_hires(:,:)
   real(r8), allocatable :: a_sol_dif_ln_hires(:,:)
   real(r8), allocatable :: a_sr_dir_ln_hires (:,:)
   real(r8), allocatable :: a_sr_dif_ln_hires (:,:)
#endif
   real(r8), allocatable :: a_sensors (:,:)

   PUBLIC :: allocate_acc_fluxes
   PUBLIC :: deallocate_acc_fluxes
   PUBLIC :: flush_acc_fluxes
   PUBLIC :: accumulate_fluxes

CONTAINS

   SUBROUTINE allocate_acc_fluxes

   USE MOD_SPMD_Task
   USE MOD_LandElm
   USE MOD_LandPatch
   USE MOD_LandUrban, only: numurban
   USE MOD_Vars_1DFluxes, only: nsensor
#ifdef CROP
   USE MOD_LandCrop
#endif
   USE MOD_Vars_Global
   IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN

            allocate (a_us        (numpatch))
            allocate (a_vs        (numpatch))
            allocate (a_t         (numpatch))
            allocate (a_q         (numpatch))
            allocate (a_prc       (numpatch))
            allocate (a_prl       (numpatch))
            allocate (a_pbot      (numpatch))
            allocate (a_frl       (numpatch))
            allocate (a_solarin   (numpatch))
            allocate (a_hpbl      (numpatch))

            allocate (a_taux      (numpatch))
            allocate (a_tauy      (numpatch))
            allocate (a_fsena     (numpatch))
            allocate (a_lfevpa    (numpatch))
            allocate (a_fevpa     (numpatch))
            allocate (a_fsenl     (numpatch))
            allocate (a_fevpl     (numpatch))
            allocate (a_etr       (numpatch))
            allocate (a_fseng     (numpatch))
            allocate (a_fevpg     (numpatch))
            allocate (a_fgrnd     (numpatch))
            allocate (a_sabvsun   (numpatch))
            allocate (a_sabvsha   (numpatch))
            allocate (a_sabg      (numpatch))
            allocate (a_olrg      (numpatch))
            allocate (a_rnet      (numpatch))
            allocate (a_xerr      (numpatch))
            allocate (a_zerr      (numpatch))
            allocate (a_rsur      (numpatch))
            allocate (a_rsur_se   (numpatch))
            allocate (a_rsur_ie   (numpatch))
            allocate (a_rsub      (numpatch))
            allocate (a_rnof      (numpatch))
#ifdef CatchLateralFlow
            allocate (a_xwsur     (numpatch))
            allocate (a_xwsub     (numpatch))
            allocate (a_fldarea   (numpatch))
#endif
            allocate (a_qintr     (numpatch))
            allocate (a_qinfl     (numpatch))
            allocate (a_qdrip     (numpatch))
            allocate (a_rstfacsun (numpatch))
            allocate (a_rstfacsha (numpatch))
            allocate (a_gssun     (numpatch))
            allocate (a_gssha     (numpatch))
            allocate (a_rss       (numpatch))
            allocate (a_wdsrf     (numpatch))

            allocate (a_zwt       (numpatch))
            allocate (a_wa        (numpatch))
            allocate (a_wat       (numpatch))
            allocate (a_wetwat    (numpatch))
            allocate (a_assim     (numpatch))
            allocate (a_respc     (numpatch))

            allocate (a_assimsun  (numpatch))
            allocate (a_assimsha  (numpatch))
            allocate (a_etrsun    (numpatch))
            allocate (a_etrsha    (numpatch))

            allocate (a_qcharge   (numpatch))

            allocate (a_t_grnd    (numpatch))
            allocate (a_tleaf     (numpatch))
            allocate (a_ldew_rain (numpatch))
            allocate (a_ldew_snow (numpatch))
            allocate (a_ldew      (numpatch))
            allocate (a_scv       (numpatch))
            allocate (a_snowdp    (numpatch))
            allocate (a_fsno      (numpatch))
            allocate (a_frcsat    (numpatch))
            allocate (a_sigf      (numpatch))
            allocate (a_green     (numpatch))
            allocate (a_lai       (numpatch))
            allocate (a_laisun    (numpatch))
            allocate (a_laisha    (numpatch))
            allocate (a_sai       (numpatch))

            allocate (a_alb       (2  ,2,numpatch))
#ifdef HYPERSPECTRAL
            allocate (a_alb_hires (211,2,numpatch))
            allocate (a_reflectance_out  (211,16,numpatch))
            allocate (a_transmittance_out(211,16,numpatch))
#endif
            allocate (a_emis      (numpatch))
            allocate (a_z0m       (numpatch))
            allocate (a_trad      (numpatch))
            allocate (a_tref      (numpatch))
            allocate (a_t2m_wmo   (numpatch))
            allocate (a_qref      (numpatch))
            allocate (a_rain      (numpatch))
            allocate (a_snow      (numpatch))

            allocate (a_o3uptakesun(numpatch))
            allocate (a_o3uptakesha(numpatch))

#ifdef DataAssimilation
            allocate (a_h2osoi_ens            (1:nl_soil,DEF_DA_ENS_NUM,numpatch))
            allocate (a_t_brt_smap_ens                (2,DEF_DA_ENS_NUM,numpatch))
            allocate (a_t_brt_fy3d_ens                (2,DEF_DA_ENS_NUM,numpatch))
            allocate (a_t_brt_smap                                   (2,numpatch))
            allocate (a_t_brt_fy3d                                   (2,numpatch))
            allocate (a_wliq_soisno_ens(maxsnl+1:nl_soil,DEF_DA_ENS_NUM,numpatch))
            allocate (a_wice_soisno_ens(maxsnl+1:nl_soil,DEF_DA_ENS_NUM,numpatch))
            allocate (a_t_soisno_ens   (maxsnl+1:nl_soil,DEF_DA_ENS_NUM,numpatch))
#endif

#ifdef URBAN_MODEL
            IF (numurban > 0) THEN
               allocate (a_t_room    (numurban))
               allocate (a_tafu      (numurban))
               allocate (a_fhac      (numurban))
               allocate (a_fwst      (numurban))
               allocate (a_fach      (numurban))
               allocate (a_fahe      (numurban))
               allocate (a_fhah      (numurban))
               allocate (a_vehc      (numurban))
               allocate (a_meta      (numurban))

               allocate (a_senroof   (numurban))
               allocate (a_senwsun   (numurban))
               allocate (a_senwsha   (numurban))
               allocate (a_sengimp   (numurban))
               allocate (a_sengper   (numurban))
               allocate (a_senurbl   (numurban))

               allocate (a_lfevproof (numurban))
               allocate (a_lfevpgimp (numurban))
               allocate (a_lfevpgper (numurban))
               allocate (a_lfevpurbl (numurban))

               allocate (a_troof     (numurban))
               allocate (a_twall     (numurban))
            ENDIF
#endif
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
            allocate (a_lai_enftemp        (numpatch))
            allocate (a_lai_enfboreal      (numpatch))
            allocate (a_lai_dnfboreal      (numpatch))
            allocate (a_lai_ebftrop        (numpatch))
            allocate (a_lai_ebftemp        (numpatch))
            allocate (a_lai_dbftrop        (numpatch))
            allocate (a_lai_dbftemp        (numpatch))
            allocate (a_lai_dbfboreal      (numpatch))
            allocate (a_lai_ebstemp        (numpatch))
            allocate (a_lai_dbstemp        (numpatch))
            allocate (a_lai_dbsboreal      (numpatch))
            allocate (a_lai_c3arcgrass     (numpatch))
            allocate (a_lai_c3grass        (numpatch))
            allocate (a_lai_c4grass        (numpatch))
#endif
#ifdef BGC
            allocate (a_leafc              (numpatch))
            allocate (a_leafc_storage      (numpatch))
            allocate (a_leafc_xfer         (numpatch))
            allocate (a_frootc             (numpatch))
            allocate (a_frootc_storage     (numpatch))
            allocate (a_frootc_xfer        (numpatch))
            allocate (a_livestemc          (numpatch))
            allocate (a_livestemc_storage  (numpatch))
            allocate (a_livestemc_xfer     (numpatch))
            allocate (a_deadstemc          (numpatch))
            allocate (a_deadstemc_storage  (numpatch))
            allocate (a_deadstemc_xfer     (numpatch))
            allocate (a_livecrootc         (numpatch))
            allocate (a_livecrootc_storage (numpatch))
            allocate (a_livecrootc_xfer    (numpatch))
            allocate (a_deadcrootc         (numpatch))
            allocate (a_deadcrootc_storage (numpatch))
            allocate (a_deadcrootc_xfer    (numpatch))
            allocate (a_grainc             (numpatch))
            allocate (a_grainc_storage     (numpatch))
            allocate (a_grainc_xfer        (numpatch))
            allocate (a_leafn              (numpatch))
            allocate (a_leafn_storage      (numpatch))
            allocate (a_leafn_xfer         (numpatch))
            allocate (a_frootn             (numpatch))
            allocate (a_frootn_storage     (numpatch))
            allocate (a_frootn_xfer        (numpatch))
            allocate (a_livestemn          (numpatch))
            allocate (a_livestemn_storage  (numpatch))
            allocate (a_livestemn_xfer     (numpatch))
            allocate (a_deadstemn          (numpatch))
            allocate (a_deadstemn_storage  (numpatch))
            allocate (a_deadstemn_xfer     (numpatch))
            allocate (a_livecrootn         (numpatch))
            allocate (a_livecrootn_storage (numpatch))
            allocate (a_livecrootn_xfer    (numpatch))
            allocate (a_deadcrootn         (numpatch))
            allocate (a_deadcrootn_storage (numpatch))
            allocate (a_deadcrootn_xfer    (numpatch))
            allocate (a_grainn             (numpatch))
            allocate (a_grainn_storage     (numpatch))
            allocate (a_grainn_xfer        (numpatch))
            allocate (a_retransn           (numpatch))
            allocate (a_gpp                (numpatch))
            allocate (a_downreg            (numpatch))
            allocate (a_ar                 (numpatch))
            allocate (a_cwdprod            (numpatch))
            allocate (a_cwddecomp          (numpatch))
            allocate (a_hr                 (numpatch))
            allocate (a_fpg                (numpatch))
            allocate (a_fpi                (numpatch))
            allocate (a_totvegc            (numpatch))
            allocate (a_totlitc            (numpatch))
            allocate (a_totcwdc            (numpatch))
            allocate (a_totsomc            (numpatch))
            allocate (a_totcolc            (numpatch))
            allocate (a_totvegn            (numpatch))
            allocate (a_totlitn            (numpatch))
            allocate (a_totcwdn            (numpatch))
            allocate (a_totsomn            (numpatch))
            allocate (a_totcoln            (numpatch))
            allocate (a_gpp_enftemp        (numpatch)) !1
            allocate (a_gpp_enfboreal      (numpatch)) !2
            allocate (a_gpp_dnfboreal      (numpatch)) !3
            allocate (a_gpp_ebftrop        (numpatch)) !4
            allocate (a_gpp_ebftemp        (numpatch)) !5
            allocate (a_gpp_dbftrop        (numpatch)) !6
            allocate (a_gpp_dbftemp        (numpatch)) !7
            allocate (a_gpp_dbfboreal      (numpatch)) !8
            allocate (a_gpp_ebstemp        (numpatch)) !9
            allocate (a_gpp_dbstemp        (numpatch)) !10
            allocate (a_gpp_dbsboreal      (numpatch)) !11
            allocate (a_gpp_c3arcgrass     (numpatch)) !12
            allocate (a_gpp_c3grass        (numpatch)) !13
            allocate (a_gpp_c4grass        (numpatch)) !14
            allocate (a_npp_enftemp        (numpatch)) !1
            allocate (a_npp_enfboreal      (numpatch)) !2
            allocate (a_npp_dnfboreal      (numpatch)) !3
            allocate (a_npp_ebftrop        (numpatch)) !4
            allocate (a_npp_ebftemp        (numpatch)) !5
            allocate (a_npp_dbftrop        (numpatch)) !6
            allocate (a_npp_dbftemp        (numpatch)) !7
            allocate (a_npp_dbfboreal      (numpatch)) !8
            allocate (a_npp_ebstemp        (numpatch)) !9
            allocate (a_npp_dbstemp        (numpatch)) !10
            allocate (a_npp_dbsboreal      (numpatch)) !11
            allocate (a_npp_c3arcgrass     (numpatch)) !12
            allocate (a_npp_c3grass        (numpatch)) !13
            allocate (a_npp_c4grass        (numpatch)) !14
            allocate (a_npptoleafc_enftemp        (numpatch)) !1
            allocate (a_npptoleafc_enfboreal      (numpatch)) !2
            allocate (a_npptoleafc_dnfboreal      (numpatch)) !3
            allocate (a_npptoleafc_ebftrop        (numpatch)) !4
            allocate (a_npptoleafc_ebftemp        (numpatch)) !5
            allocate (a_npptoleafc_dbftrop        (numpatch)) !6
            allocate (a_npptoleafc_dbftemp        (numpatch)) !7
            allocate (a_npptoleafc_dbfboreal      (numpatch)) !8
            allocate (a_npptoleafc_ebstemp        (numpatch)) !9
            allocate (a_npptoleafc_dbstemp        (numpatch)) !10
            allocate (a_npptoleafc_dbsboreal      (numpatch)) !11
            allocate (a_npptoleafc_c3arcgrass     (numpatch)) !12
            allocate (a_npptoleafc_c3grass        (numpatch)) !13
            allocate (a_npptoleafc_c4grass        (numpatch)) !14
            allocate (a_leafc_enftemp      (numpatch)) !1
            allocate (a_leafc_enfboreal    (numpatch)) !2
            allocate (a_leafc_dnfboreal    (numpatch)) !3
            allocate (a_leafc_ebftrop      (numpatch)) !4
            allocate (a_leafc_ebftemp      (numpatch)) !5
            allocate (a_leafc_dbftrop      (numpatch)) !6
            allocate (a_leafc_dbftemp      (numpatch)) !7
            allocate (a_leafc_dbfboreal    (numpatch)) !8
            allocate (a_leafc_ebstemp      (numpatch)) !9
            allocate (a_leafc_dbstemp      (numpatch)) !10
            allocate (a_leafc_dbsboreal    (numpatch)) !11
            allocate (a_leafc_c3arcgrass   (numpatch)) !12
            allocate (a_leafc_c3grass      (numpatch)) !13
            allocate (a_leafc_c4grass      (numpatch)) !14

            allocate (a_O2_DECOMP_DEPTH_UNSAT (1:nl_soil,numpatch))
            allocate (a_CONC_O2_UNSAT         (1:nl_soil,numpatch))

#ifdef CROP
            allocate (a_pdcorn             (numpatch))
            allocate (a_pdswheat           (numpatch))
            allocate (a_pdwwheat           (numpatch))
            allocate (a_pdsoybean          (numpatch))
            allocate (a_pdcotton           (numpatch))
            allocate (a_pdrice1            (numpatch))
            allocate (a_pdrice2            (numpatch))
            allocate (a_pdsugarcane        (numpatch))
            allocate (a_plantdate          (numpatch))
            allocate (a_manunitro          (numpatch))
            allocate (a_fertnitro_corn     (numpatch))
            allocate (a_fertnitro_swheat   (numpatch))
            allocate (a_fertnitro_wwheat   (numpatch))
            allocate (a_fertnitro_soybean  (numpatch))
            allocate (a_fertnitro_cotton   (numpatch))
            allocate (a_fertnitro_rice1    (numpatch))
            allocate (a_fertnitro_rice2    (numpatch))
            allocate (a_fertnitro_sugarcane(numpatch))
            allocate (a_irrig_method_corn     (numpatch))
            allocate (a_irrig_method_swheat   (numpatch))
            allocate (a_irrig_method_wwheat   (numpatch))
            allocate (a_irrig_method_soybean  (numpatch))
            allocate (a_irrig_method_cotton   (numpatch))
            allocate (a_irrig_method_rice1    (numpatch))
            allocate (a_irrig_method_rice2    (numpatch))
            allocate (a_irrig_method_sugarcane(numpatch))
            allocate (a_cphase             (numpatch))
            allocate (a_hui                (numpatch))
            allocate (a_gddmaturity        (numpatch))
            allocate (a_gddplant           (numpatch))
            allocate (a_vf                 (numpatch))
            allocate (a_cropprod1c         (numpatch))
            allocate (a_cropprod1c_loss    (numpatch))
            allocate (a_cropseedc_deficit  (numpatch))
            allocate (a_grainc_to_cropprodc(numpatch))
            allocate (a_grainc_to_seed     (numpatch))
            allocate (a_fert_to_sminn      (numpatch))

            allocate (a_sum_irrig          (numpatch))
            allocate (a_sum_deficit_irrig  (numpatch))
            allocate (a_sum_irrig_count    (numpatch))
            allocate (a_waterstorage       (numpatch))
            allocate (a_groundwater_demand (numpatch))
            allocate (a_groundwater_supply (numpatch))
            allocate (a_reservoirriver_demand(numpatch))
            allocate (a_reservoirriver_supply(numpatch))
            allocate (a_reservoir_supply     (numpatch))
            allocate (a_river_supply         (numpatch))
            allocate (a_runoff_supply        (numpatch))
#endif
            allocate (a_ndep_to_sminn      (numpatch))

            allocate (a_abm                (numpatch))
            allocate (a_gdp                (numpatch))
            allocate (a_peatf              (numpatch))
            allocate (a_hdm                (numpatch))
            allocate (a_lnfm               (numpatch))

            allocate (a_leafcCap              (numpatch))
            allocate (a_leafc_storageCap      (numpatch))
            allocate (a_leafc_xferCap         (numpatch))
            allocate (a_frootcCap             (numpatch))
            allocate (a_frootc_storageCap     (numpatch))
            allocate (a_frootc_xferCap        (numpatch))
            allocate (a_livestemcCap          (numpatch))
            allocate (a_livestemc_storageCap  (numpatch))
            allocate (a_livestemc_xferCap     (numpatch))
            allocate (a_deadstemcCap          (numpatch))
            allocate (a_deadstemc_storageCap  (numpatch))
            allocate (a_deadstemc_xferCap     (numpatch))
            allocate (a_livecrootcCap         (numpatch))
            allocate (a_livecrootc_storageCap (numpatch))
            allocate (a_livecrootc_xferCap    (numpatch))
            allocate (a_deadcrootcCap         (numpatch))
            allocate (a_deadcrootc_storageCap (numpatch))
            allocate (a_deadcrootc_xferCap    (numpatch))
            allocate (a_leafnCap              (numpatch))
            allocate (a_leafn_storageCap      (numpatch))
            allocate (a_leafn_xferCap         (numpatch))
            allocate (a_frootnCap             (numpatch))
            allocate (a_frootn_storageCap     (numpatch))
            allocate (a_frootn_xferCap        (numpatch))
            allocate (a_livestemnCap          (numpatch))
            allocate (a_livestemn_storageCap  (numpatch))
            allocate (a_livestemn_xferCap     (numpatch))
            allocate (a_deadstemnCap          (numpatch))
            allocate (a_deadstemn_storageCap  (numpatch))
            allocate (a_deadstemn_xferCap     (numpatch))
            allocate (a_livecrootnCap         (numpatch))
            allocate (a_livecrootn_storageCap (numpatch))
            allocate (a_livecrootn_xferCap    (numpatch))
            allocate (a_deadcrootnCap         (numpatch))
            allocate (a_deadcrootn_storageCap (numpatch))
            allocate (a_deadcrootn_xferCap    (numpatch))
#endif
! Ozone stress variables
            allocate (a_ozone              (numpatch))
! End ozone stress variables
            allocate (a_t_soisno    (maxsnl+1:nl_soil,numpatch))
            allocate (a_wliq_soisno (maxsnl+1:nl_soil,numpatch))
            allocate (a_wice_soisno (maxsnl+1:nl_soil,numpatch))
            allocate (a_h2osoi      (1:nl_soil,       numpatch))
            allocate (a_qlayer      (0:nl_soil,       numpatch))
            allocate (a_lake_deficit                 (numpatch))
            allocate (a_rootr       (1:nl_soil,       numpatch))
            allocate (a_BD_all      (1:nl_soil,       numpatch))
            allocate (a_wfc         (1:nl_soil,       numpatch))
            allocate (a_OM_density  (1:nl_soil,       numpatch))
!Plant Hydraulic parameters
            allocate (a_vegwp       (1:nvegwcs,       numpatch))
!End Plant Hydraulic parameters
            allocate (a_dz_lake     (nl_lake,         numpatch))
            allocate (a_t_lake      (nl_lake,         numpatch))
            allocate (a_lake_icefrac(nl_lake,         numpatch))

#ifdef BGC
            allocate (a_litr1c_vr   (1:nl_soil,       numpatch))
            allocate (a_litr2c_vr   (1:nl_soil,       numpatch))
            allocate (a_litr3c_vr   (1:nl_soil,       numpatch))
            allocate (a_soil1c_vr   (1:nl_soil,       numpatch))
            allocate (a_soil2c_vr   (1:nl_soil,       numpatch))
            allocate (a_soil3c_vr   (1:nl_soil,       numpatch))
            allocate (a_cwdc_vr     (1:nl_soil,       numpatch))
            allocate (a_litr1n_vr   (1:nl_soil,       numpatch))
            allocate (a_litr2n_vr   (1:nl_soil,       numpatch))
            allocate (a_litr3n_vr   (1:nl_soil,       numpatch))
            allocate (a_soil1n_vr   (1:nl_soil,       numpatch))
            allocate (a_soil2n_vr   (1:nl_soil,       numpatch))
            allocate (a_soil3n_vr   (1:nl_soil,       numpatch))
            allocate (a_cwdn_vr     (1:nl_soil,       numpatch))
            allocate (a_totsoiln_vr (1:nl_soil,       numpatch))
            allocate (a_sminn_vr    (1:nl_soil,       numpatch))
            allocate (decomp_vr_tmp (1:nl_soil,       numpatch))

            allocate (a_litr1cCap_vr(1:nl_soil,       numpatch))
            allocate (a_litr2cCap_vr(1:nl_soil,       numpatch))
            allocate (a_litr3cCap_vr(1:nl_soil,       numpatch))
            allocate (a_soil1cCap_vr(1:nl_soil,       numpatch))
            allocate (a_soil2cCap_vr(1:nl_soil,       numpatch))
            allocate (a_soil3cCap_vr(1:nl_soil,       numpatch))
            allocate (a_cwdcCap_vr  (1:nl_soil,       numpatch))
            allocate (a_litr1nCap_vr(1:nl_soil,       numpatch))
            allocate (a_litr2nCap_vr(1:nl_soil,       numpatch))
            allocate (a_litr3nCap_vr(1:nl_soil,       numpatch))
            allocate (a_soil1nCap_vr(1:nl_soil,       numpatch))
            allocate (a_soil2nCap_vr(1:nl_soil,       numpatch))
            allocate (a_soil3nCap_vr(1:nl_soil,       numpatch))
            allocate (a_cwdnCap_vr  (1:nl_soil,       numpatch))
            allocate (a_t_scalar    (1:nl_soil,       numpatch))
            allocate (a_w_scalar    (1:nl_soil,       numpatch))
#endif

            allocate (a_ustar     (numpatch))
            allocate (a_ustar2    (numpatch))
            allocate (a_tstar     (numpatch))
            allocate (a_qstar     (numpatch))
            allocate (a_zol       (numpatch))
            allocate (a_rib       (numpatch))
            allocate (a_fm        (numpatch))
            allocate (a_fh        (numpatch))
            allocate (a_fq        (numpatch))

            allocate (a_us10m     (numpatch))
            allocate (a_vs10m     (numpatch))
            allocate (a_fm10m     (numpatch))

            allocate (a_sr        (numpatch))
            allocate (a_solvd     (numpatch))
            allocate (a_solvi     (numpatch))
            allocate (a_solnd     (numpatch))
            allocate (a_solni     (numpatch))
            allocate (a_srvd      (numpatch))
            allocate (a_srvi      (numpatch))
            allocate (a_srnd      (numpatch))
            allocate (a_srni      (numpatch))
            allocate (a_solvdln   (numpatch))
            allocate (a_solviln   (numpatch))
            allocate (a_solndln   (numpatch))
            allocate (a_solniln   (numpatch))
            allocate (a_srvdln    (numpatch))
            allocate (a_srviln    (numpatch))
            allocate (a_srndln    (numpatch))
            allocate (a_srniln    (numpatch))
#ifdef HYPERSPECTRAL
            allocate (a_sol_dir_ln_hires(211, numpatch))
            allocate (a_sol_dif_ln_hires(211, numpatch))
            allocate (a_sr_dir_ln_hires (211, numpatch))
            allocate (a_sr_dif_ln_hires (211, numpatch))
#endif
            allocate (a_sensors (nsensor,numpatch))

            allocate (nac_ln      (numpatch))
            allocate (nac_dt      (numpatch))
            allocate (filter_dt   (numpatch))

         ENDIF
      ENDIF

#ifdef EXTERNAL_LAKE
      CALL allocate_LakeAccVars
#endif

      IF (p_is_worker) THEN
         CALL elm_patch%build (landelm, landpatch, use_frac = .true.)
      ENDIF

   END SUBROUTINE allocate_acc_fluxes

   SUBROUTINE deallocate_acc_fluxes ()

   USE MOD_SPMD_Task
   USE MOD_LandPatch, only: numpatch
   USE MOD_LandUrban, only: numurban
   IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN

            deallocate (a_us     )
            deallocate (a_vs     )
            deallocate (a_t      )
            deallocate (a_q      )
            deallocate (a_prc    )
            deallocate (a_prl    )
            deallocate (a_pbot   )
            deallocate (a_frl    )
            deallocate (a_solarin)
            deallocate (a_hpbl   )

            deallocate (a_taux      )
            deallocate (a_tauy      )
            deallocate (a_fsena     )
            deallocate (a_lfevpa    )
            deallocate (a_fevpa     )
            deallocate (a_fsenl     )
            deallocate (a_fevpl     )
            deallocate (a_etr       )
            deallocate (a_fseng     )
            deallocate (a_fevpg     )
            deallocate (a_fgrnd     )
            deallocate (a_sabvsun   )
            deallocate (a_sabvsha   )
            deallocate (a_sabg      )
            deallocate (a_olrg      )
            deallocate (a_rnet      )
            deallocate (a_xerr      )
            deallocate (a_zerr      )
            deallocate (a_rsur      )
            deallocate (a_rsur_se   )
            deallocate (a_rsur_ie   )
            deallocate (a_rsub      )
            deallocate (a_rnof      )
#ifdef CatchLateralFlow
            deallocate (a_xwsur     )
            deallocate (a_xwsub     )
            deallocate (a_fldarea   )
#endif
            deallocate (a_qintr     )
            deallocate (a_qinfl     )
            deallocate (a_qdrip     )
            deallocate (a_rstfacsun )
            deallocate (a_rstfacsha )
            deallocate (a_gssun     )
            deallocate (a_gssha     )
            deallocate (a_rss       )
            deallocate (a_wdsrf     )

            deallocate (a_zwt       )
            deallocate (a_wa        )
            deallocate (a_wat       )
            deallocate (a_wetwat    )
            deallocate (a_assim     )
            deallocate (a_respc     )

            deallocate (a_assimsun  ) !1
            deallocate (a_assimsha  ) !1
            deallocate (a_etrsun    ) !1
            deallocate (a_etrsha    ) !1

            deallocate (a_qcharge   )

            deallocate (a_t_grnd    )
            deallocate (a_tleaf     )
            deallocate (a_ldew_rain )
            deallocate (a_ldew_snow )
            deallocate (a_ldew      )
            deallocate (a_scv       )
            deallocate (a_snowdp    )
            deallocate (a_fsno      )
            deallocate (a_frcsat    )
            deallocate (a_sigf      )
            deallocate (a_green     )
            deallocate (a_lai       )
            deallocate (a_laisun    )
            deallocate (a_laisha    )
            deallocate (a_sai       )

            deallocate (a_alb  )
#ifdef HYPERSPECTRAL
            deallocate (a_alb_hires  )
            deallocate (a_reflectance_out  )
            deallocate (a_transmittance_out)
#endif

            deallocate (a_emis      )
            deallocate (a_z0m       )
            deallocate (a_trad      )
            deallocate (a_tref      )
            deallocate (a_t2m_wmo   )
            deallocate (a_qref      )
            deallocate (a_rain      )
            deallocate (a_snow      )

            deallocate (a_o3uptakesun)
            deallocate (a_o3uptakesha)

#ifdef DataAssimilation
            deallocate (a_h2osoi_ens     )
            deallocate (a_t_brt_smap_ens )
            deallocate (a_t_brt_fy3d_ens )
            deallocate (a_t_brt_fy3d     )
            deallocate (a_t_brt_smap     )
            deallocate (a_wliq_soisno_ens)
            deallocate (a_wice_soisno_ens)
            deallocate (a_t_soisno_ens   )
#endif

#ifdef URBAN_MODEL
            IF (numurban > 0) THEN
               deallocate (a_t_room    )
               deallocate (a_tafu      )
               deallocate (a_fhac      )
               deallocate (a_fwst      )
               deallocate (a_fach      )
               deallocate (a_fahe      )
               deallocate (a_fhah      )
               deallocate (a_vehc      )
               deallocate (a_meta      )

               deallocate (a_senroof   )
               deallocate (a_senwsun   )
               deallocate (a_senwsha   )
               deallocate (a_sengimp   )
               deallocate (a_sengper   )
               deallocate (a_senurbl   )

               deallocate (a_lfevproof )
               deallocate (a_lfevpgimp )
               deallocate (a_lfevpgper )
               deallocate (a_lfevpurbl )

               deallocate (a_troof     )
               deallocate (a_twall     )
            ENDIF
#endif

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
            deallocate (a_lai_enftemp        )
            deallocate (a_lai_enfboreal      )
            deallocate (a_lai_dnfboreal      )
            deallocate (a_lai_ebftrop        )
            deallocate (a_lai_ebftemp        )
            deallocate (a_lai_dbftrop        )
            deallocate (a_lai_dbftemp        )
            deallocate (a_lai_dbfboreal      )
            deallocate (a_lai_ebstemp        )
            deallocate (a_lai_dbstemp        )
            deallocate (a_lai_dbsboreal      )
            deallocate (a_lai_c3arcgrass     )
            deallocate (a_lai_c3grass        )
            deallocate (a_lai_c4grass        )
#endif

#ifdef BGC
            deallocate (a_leafc              )
            deallocate (a_leafc_storage      )
            deallocate (a_leafc_xfer         )
            deallocate (a_frootc             )
            deallocate (a_frootc_storage     )
            deallocate (a_frootc_xfer        )
            deallocate (a_livestemc          )
            deallocate (a_livestemc_storage  )
            deallocate (a_livestemc_xfer     )
            deallocate (a_deadstemc          )
            deallocate (a_deadstemc_storage  )
            deallocate (a_deadstemc_xfer     )
            deallocate (a_livecrootc         )
            deallocate (a_livecrootc_storage )
            deallocate (a_livecrootc_xfer    )
            deallocate (a_deadcrootc         )
            deallocate (a_deadcrootc_storage )
            deallocate (a_deadcrootc_xfer    )
            deallocate (a_grainc             )
            deallocate (a_grainc_storage     )
            deallocate (a_grainc_xfer        )
            deallocate (a_leafn              )
            deallocate (a_leafn_storage      )
            deallocate (a_leafn_xfer         )
            deallocate (a_frootn             )
            deallocate (a_frootn_storage     )
            deallocate (a_frootn_xfer        )
            deallocate (a_livestemn          )
            deallocate (a_livestemn_storage  )
            deallocate (a_livestemn_xfer     )
            deallocate (a_deadstemn          )
            deallocate (a_deadstemn_storage  )
            deallocate (a_deadstemn_xfer     )
            deallocate (a_livecrootn         )
            deallocate (a_livecrootn_storage )
            deallocate (a_livecrootn_xfer    )
            deallocate (a_deadcrootn         )
            deallocate (a_deadcrootn_storage )
            deallocate (a_deadcrootn_xfer    )
            deallocate (a_grainn             )
            deallocate (a_grainn_storage     )
            deallocate (a_grainn_xfer        )
            deallocate (a_retransn           )
            deallocate (a_gpp                )
            deallocate (a_downreg            )
            deallocate (a_ar                 )
            deallocate (a_cwdprod            )
            deallocate (a_cwddecomp          )
            deallocate (a_hr                 )
            deallocate (a_fpg                )
            deallocate (a_fpi                )
            deallocate (a_totvegc            )
            deallocate (a_totlitc            )
            deallocate (a_totcwdc            )
            deallocate (a_totsomc            )
            deallocate (a_totcolc            )
            deallocate (a_totvegn            )
            deallocate (a_totlitn            )
            deallocate (a_totcwdn            )
            deallocate (a_totsomn            )
            deallocate (a_totcoln            )
            deallocate (a_gpp_enftemp        ) !1
            deallocate (a_gpp_enfboreal      ) !2
            deallocate (a_gpp_dnfboreal      ) !3
            deallocate (a_gpp_ebftrop        ) !4
            deallocate (a_gpp_ebftemp        ) !5
            deallocate (a_gpp_dbftrop        ) !6
            deallocate (a_gpp_dbftemp        ) !7
            deallocate (a_gpp_dbfboreal      ) !8
            deallocate (a_gpp_ebstemp        ) !9
            deallocate (a_gpp_dbstemp        ) !10
            deallocate (a_gpp_dbsboreal      ) !11
            deallocate (a_gpp_c3arcgrass     ) !12
            deallocate (a_gpp_c3grass        ) !13
            deallocate (a_gpp_c4grass        ) !14
            deallocate (a_npp_enftemp        ) !1
            deallocate (a_npp_enfboreal      ) !2
            deallocate (a_npp_dnfboreal      ) !3
            deallocate (a_npp_ebftrop        ) !4
            deallocate (a_npp_ebftemp        ) !5
            deallocate (a_npp_dbftrop        ) !6
            deallocate (a_npp_dbftemp        ) !7
            deallocate (a_npp_dbfboreal      ) !8
            deallocate (a_npp_ebstemp        ) !9
            deallocate (a_npp_dbstemp        ) !10
            deallocate (a_npp_dbsboreal      ) !11
            deallocate (a_npp_c3arcgrass     ) !12
            deallocate (a_npp_c3grass        ) !13
            deallocate (a_npp_c4grass        ) !14
            deallocate (a_npptoleafc_enftemp        ) !1
            deallocate (a_npptoleafc_enfboreal      ) !2
            deallocate (a_npptoleafc_dnfboreal      ) !3
            deallocate (a_npptoleafc_ebftrop        ) !4
            deallocate (a_npptoleafc_ebftemp        ) !5
            deallocate (a_npptoleafc_dbftrop        ) !6
            deallocate (a_npptoleafc_dbftemp        ) !7
            deallocate (a_npptoleafc_dbfboreal      ) !8
            deallocate (a_npptoleafc_ebstemp        ) !9
            deallocate (a_npptoleafc_dbstemp        ) !10
            deallocate (a_npptoleafc_dbsboreal      ) !11
            deallocate (a_npptoleafc_c3arcgrass     ) !12
            deallocate (a_npptoleafc_c3grass        ) !13
            deallocate (a_npptoleafc_c4grass        ) !14
            deallocate (a_leafc_enftemp      ) !1
            deallocate (a_leafc_enfboreal    ) !2
            deallocate (a_leafc_dnfboreal    ) !3
            deallocate (a_leafc_ebftrop      ) !4
            deallocate (a_leafc_ebftemp      ) !5
            deallocate (a_leafc_dbftrop      ) !6
            deallocate (a_leafc_dbftemp      ) !7
            deallocate (a_leafc_dbfboreal    ) !8
            deallocate (a_leafc_ebstemp      ) !9
            deallocate (a_leafc_dbstemp      ) !10
            deallocate (a_leafc_dbsboreal    ) !11
            deallocate (a_leafc_c3arcgrass   ) !12
            deallocate (a_leafc_c3grass      ) !13
            deallocate (a_leafc_c4grass      ) !14

            deallocate (a_O2_DECOMP_DEPTH_UNSAT )
            deallocate (a_CONC_O2_UNSAT         )

#ifdef CROP
            deallocate (a_pdcorn             )
            deallocate (a_pdswheat           )
            deallocate (a_pdwwheat           )
            deallocate (a_pdsoybean          )
            deallocate (a_pdcotton           )
            deallocate (a_pdrice1            )
            deallocate (a_pdrice2            )
            deallocate (a_pdsugarcane        )
            deallocate (a_plantdate          )
            deallocate (a_manunitro          )
            deallocate (a_fertnitro_corn     )
            deallocate (a_fertnitro_swheat   )
            deallocate (a_fertnitro_wwheat   )
            deallocate (a_fertnitro_soybean  )
            deallocate (a_fertnitro_cotton   )
            deallocate (a_fertnitro_rice1    )
            deallocate (a_fertnitro_rice2    )
            deallocate (a_fertnitro_sugarcane)
            deallocate (a_irrig_method_corn     )
            deallocate (a_irrig_method_swheat   )
            deallocate (a_irrig_method_wwheat   )
            deallocate (a_irrig_method_soybean  )
            deallocate (a_irrig_method_cotton   )
            deallocate (a_irrig_method_rice1    )
            deallocate (a_irrig_method_rice2    )
            deallocate (a_irrig_method_sugarcane)
            deallocate (a_cphase             )
            deallocate (a_hui                )
            deallocate (a_vf                 )
            deallocate (a_gddmaturity        )
            deallocate (a_gddplant           )
            deallocate (a_cropprod1c         )
            deallocate (a_cropprod1c_loss    )
            deallocate (a_cropseedc_deficit  )
            deallocate (a_grainc_to_cropprodc)
            deallocate (a_grainc_to_seed     )
            deallocate (a_fert_to_sminn      )

            deallocate (a_sum_irrig          )
            deallocate (a_sum_deficit_irrig  )
            deallocate (a_sum_irrig_count    )
            deallocate (a_waterstorage       )
            deallocate (a_groundwater_demand )
            deallocate (a_groundwater_supply )
            deallocate (a_reservoirriver_demand)
            deallocate (a_reservoirriver_supply)
            deallocate (a_reservoir_supply     )
            deallocate (a_river_supply         )
            deallocate (a_runoff_supply        )
#endif
            deallocate (a_ndep_to_sminn      )

            deallocate (a_abm                )
            deallocate (a_gdp                )
            deallocate (a_peatf              )
            deallocate (a_hdm                )
            deallocate (a_lnfm               )

            deallocate (a_leafcCap             )
            deallocate (a_leafc_storageCap     )
            deallocate (a_leafc_xferCap        )
            deallocate (a_frootcCap            )
            deallocate (a_frootc_storageCap    )
            deallocate (a_frootc_xferCap       )
            deallocate (a_livestemcCap         )
            deallocate (a_livestemc_storageCap )
            deallocate (a_livestemc_xferCap    )
            deallocate (a_deadstemcCap         )
            deallocate (a_deadstemc_storageCap )
            deallocate (a_deadstemc_xferCap    )
            deallocate (a_livecrootcCap        )
            deallocate (a_livecrootc_storageCap)
            deallocate (a_livecrootc_xferCap   )
            deallocate (a_deadcrootcCap        )
            deallocate (a_deadcrootc_storageCap)
            deallocate (a_deadcrootc_xferCap   )
            deallocate (a_leafnCap             )
            deallocate (a_leafn_storageCap     )
            deallocate (a_leafn_xferCap        )
            deallocate (a_frootnCap            )
            deallocate (a_frootn_storageCap    )
            deallocate (a_frootn_xferCap       )
            deallocate (a_livestemnCap         )
            deallocate (a_livestemn_storageCap )
            deallocate (a_livestemn_xferCap    )
            deallocate (a_deadstemnCap         )
            deallocate (a_deadstemn_storageCap )
            deallocate (a_deadstemn_xferCap    )
            deallocate (a_livecrootnCap        )
            deallocate (a_livecrootn_storageCap)
            deallocate (a_livecrootn_xferCap   )
            deallocate (a_deadcrootnCap        )
            deallocate (a_deadcrootn_storageCap)
            deallocate (a_deadcrootn_xferCap   )
#endif
! Ozone stress variables
            deallocate (a_ozone              )
! END ozone stress variables

            deallocate (a_t_soisno    )
            deallocate (a_wliq_soisno )
            deallocate (a_wice_soisno )
            deallocate (a_h2osoi      )
            deallocate (a_qlayer      )
            deallocate (a_lake_deficit)
            deallocate (a_rootr       )
            deallocate (a_BD_all      )
            deallocate (a_wfc         )
            deallocate (a_OM_density  )
!Plant Hydraulic parameters
            deallocate (a_vegwp       )
!END plant hydraulic parameters
            deallocate (a_dz_lake     )
            deallocate (a_t_lake      )
            deallocate (a_lake_icefrac)
#ifdef BGC
            deallocate (a_litr1c_vr   )
            deallocate (a_litr2c_vr   )
            deallocate (a_litr3c_vr   )
            deallocate (a_soil1c_vr   )
            deallocate (a_soil2c_vr   )
            deallocate (a_soil3c_vr   )
            deallocate (a_cwdc_vr     )
            deallocate (a_litr1n_vr   )
            deallocate (a_litr2n_vr   )
            deallocate (a_litr3n_vr   )
            deallocate (a_soil1n_vr   )
            deallocate (a_soil2n_vr   )
            deallocate (a_soil3n_vr   )
            deallocate (a_cwdn_vr     )
            deallocate (a_totsoiln_vr )
            deallocate (a_sminn_vr    )
            deallocate (decomp_vr_tmp )
            deallocate (a_litr1cCap_vr)
            deallocate (a_litr2cCap_vr)
            deallocate (a_litr3cCap_vr)
            deallocate (a_soil1cCap_vr)
            deallocate (a_soil2cCap_vr)
            deallocate (a_soil3cCap_vr)
            deallocate (a_cwdcCap_vr  )
            deallocate (a_litr1nCap_vr)
            deallocate (a_litr2nCap_vr)
            deallocate (a_litr3nCap_vr)
            deallocate (a_soil1nCap_vr)
            deallocate (a_soil2nCap_vr)
            deallocate (a_soil3nCap_vr)
            deallocate (a_cwdnCap_vr  )
            deallocate (a_t_scalar    )
            deallocate (a_w_scalar    )
#endif

            deallocate (a_ustar     )
            deallocate (a_ustar2    )
            deallocate (a_tstar     )
            deallocate (a_qstar     )
            deallocate (a_zol       )
            deallocate (a_rib       )
            deallocate (a_fm        )
            deallocate (a_fh        )
            deallocate (a_fq        )

            deallocate (a_us10m     )
            deallocate (a_vs10m     )
            deallocate (a_fm10m     )

            deallocate (a_sr        )
            deallocate (a_solvd     )
            deallocate (a_solvi     )
            deallocate (a_solnd     )
            deallocate (a_solni     )
            deallocate (a_srvd      )
            deallocate (a_srvi      )
            deallocate (a_srnd      )
            deallocate (a_srni      )
            deallocate (a_solvdln   )
            deallocate (a_solviln   )
            deallocate (a_solndln   )
            deallocate (a_solniln   )
            deallocate (a_srvdln    )
            deallocate (a_srviln    )
            deallocate (a_srndln    )
            deallocate (a_srniln    )
#ifdef HYPERSPECTRAL
            deallocate (a_sol_dir_ln_hires)
            deallocate (a_sol_dif_ln_hires)
            deallocate (a_sr_dir_ln_hires )
            deallocate (a_sr_dif_ln_hires )
#endif
            deallocate (a_sensors   )

            deallocate (nac_ln      )
            deallocate (nac_dt      )
            deallocate (filter_dt   )

         ENDIF
      ENDIF

#ifdef EXTERNAL_LAKE
      CALL deallocate_LakeAccVars
#endif

   END SUBROUTINE deallocate_acc_fluxes

   !-----------------------
   SUBROUTINE FLUSH_acc_fluxes ()

      USE MOD_SPMD_Task
      USE MOD_LandPatch, only: numpatch
      USE MOD_LandUrban, only: numurban
      USE MOD_Vars_Global, only: spval
      IMPLICIT NONE

      IF (p_is_worker) THEN

         nac = 0

         IF (numpatch > 0) THEN

            ! flush the Fluxes for accumulation
            a_us        (:) = spval
            a_vs        (:) = spval
            a_t         (:) = spval
            a_q         (:) = spval
            a_prc       (:) = spval
            a_prl       (:) = spval
            a_pbot      (:) = spval
            a_frl       (:) = spval
            a_solarin   (:) = spval
            a_hpbl      (:) = spval

            a_taux      (:) = spval
            a_tauy      (:) = spval
            a_fsena     (:) = spval
            a_lfevpa    (:) = spval
            a_fevpa     (:) = spval
            a_fsenl     (:) = spval
            a_fevpl     (:) = spval
            a_etr       (:) = spval
            a_fseng     (:) = spval
            a_fevpg     (:) = spval
            a_fgrnd     (:) = spval
            a_sabvsun   (:) = spval
            a_sabvsha   (:) = spval
            a_sabg      (:) = spval
            a_olrg      (:) = spval
            a_rnet      (:) = spval
            a_xerr      (:) = spval
            a_zerr      (:) = spval
            a_rsur      (:) = spval
            a_rsur_se   (:) = spval
            a_rsur_ie   (:) = spval
            a_rsub      (:) = spval
            a_rnof      (:) = spval
#ifdef CatchLateralFlow
            a_xwsur     (:) = spval
            a_xwsub     (:) = spval
            a_fldarea   (:) = spval
#endif
            a_qintr     (:) = spval
            a_qinfl     (:) = spval
            a_qdrip     (:) = spval
            a_rstfacsun (:) = spval
            a_rstfacsha (:) = spval
            a_gssun     (:) = spval
            a_gssha     (:) = spval
            a_rss       (:) = spval

            a_wdsrf     (:) = spval
            a_zwt       (:) = spval
            a_wa        (:) = spval
            a_wat       (:) = spval
            a_wetwat    (:) = spval
            a_assim     (:) = spval
            a_respc     (:) = spval
            a_assimsun  (:) = spval
            a_assimsha  (:) = spval
            a_etrsun    (:) = spval
            a_etrsha    (:) = spval

            a_qcharge   (:) = spval

            a_t_grnd    (:) = spval
            a_tleaf     (:) = spval
            a_ldew_rain (:) = spval
            a_ldew_snow (:) = spval
            a_ldew      (:) = spval
            a_scv       (:) = spval
            a_snowdp    (:) = spval
            a_fsno      (:) = spval
            a_frcsat    (:) = spval
            a_sigf      (:) = spval
            a_green     (:) = spval
            a_lai       (:) = spval
            a_laisun    (:) = spval
            a_laisha    (:) = spval
            a_sai       (:) = spval

            a_alb   (:,:,:) = spval

#ifdef HYPERSPECTRAL
            a_alb_hires (:,:,:) = spval
#endif

            a_emis      (:) = spval
            a_z0m       (:) = spval
            a_trad      (:) = spval
            a_tref      (:) = spval
            a_t2m_wmo   (:) = spval
            a_qref      (:) = spval
            a_rain      (:) = spval
            a_snow      (:) = spval

            a_o3uptakesun(:) = spval
            a_o3uptakesha(:) = spval

#ifdef DataAssimilation
            a_h2osoi_ens     (:,:,:) = spval
            a_t_brt_smap_ens (:,:,:) = spval
            a_t_brt_fy3d_ens (:,:,:) = spval
            a_t_brt_fy3d       (:,:) = spval
            a_t_brt_smap       (:,:) = spval
            a_wliq_soisno_ens(:,:,:) = spval
            a_wice_soisno_ens(:,:,:) = spval
            a_t_soisno_ens   (:,:,:) = spval
#endif

#ifdef URBAN_MODEL
            IF (numurban > 0) THEN
               a_t_room   (:) = spval
               a_tafu     (:) = spval
               a_fhac     (:) = spval
               a_fwst     (:) = spval
               a_fach     (:) = spval
               a_fahe     (:) = spval
               a_fhah     (:) = spval
               a_vehc     (:) = spval
               a_meta     (:) = spval

               a_senroof  (:) = spval
               a_senwsun  (:) = spval
               a_senwsha  (:) = spval
               a_sengimp  (:) = spval
               a_sengper  (:) = spval
               a_senurbl  (:) = spval

               a_lfevproof(:) = spval
               a_lfevpgimp(:) = spval
               a_lfevpgper(:) = spval
               a_lfevpurbl(:) = spval

               a_troof    (:) = spval
               a_twall    (:) = spval
            ENDIF
#endif

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
            a_lai_enftemp      (:) = spval
            a_lai_enfboreal    (:) = spval
            a_lai_dnfboreal    (:) = spval
            a_lai_ebftrop      (:) = spval
            a_lai_ebftemp      (:) = spval
            a_lai_dbftrop      (:) = spval
            a_lai_dbftemp      (:) = spval
            a_lai_dbfboreal    (:) = spval
            a_lai_ebstemp      (:) = spval
            a_lai_dbstemp      (:) = spval
            a_lai_dbsboreal    (:) = spval
            a_lai_c3arcgrass   (:) = spval
            a_lai_c3grass      (:) = spval
            a_lai_c4grass      (:) = spval
#endif
#ifdef BGC
            a_leafc              (:) = spval
            a_leafc_storage      (:) = spval
            a_leafc_xfer         (:) = spval
            a_frootc             (:) = spval
            a_frootc_storage     (:) = spval
            a_frootc_xfer        (:) = spval
            a_livestemc          (:) = spval
            a_livestemc_storage  (:) = spval
            a_livestemc_xfer     (:) = spval
            a_deadstemc          (:) = spval
            a_deadstemc_storage  (:) = spval
            a_deadstemc_xfer     (:) = spval
            a_livecrootc         (:) = spval
            a_livecrootc_storage (:) = spval
            a_livecrootc_xfer    (:) = spval
            a_deadcrootc         (:) = spval
            a_deadcrootc_storage (:) = spval
            a_deadcrootc_xfer    (:) = spval
            a_grainc             (:) = spval
            a_grainc_storage     (:) = spval
            a_grainc_xfer        (:) = spval
            a_leafn              (:) = spval
            a_leafn_storage      (:) = spval
            a_leafn_xfer         (:) = spval
            a_frootn             (:) = spval
            a_frootn_storage     (:) = spval
            a_frootn_xfer        (:) = spval
            a_livestemn          (:) = spval
            a_livestemn_storage  (:) = spval
            a_livestemn_xfer     (:) = spval
            a_deadstemn          (:) = spval
            a_deadstemn_storage  (:) = spval
            a_deadstemn_xfer     (:) = spval
            a_livecrootn         (:) = spval
            a_livecrootn_storage (:) = spval
            a_livecrootn_xfer    (:) = spval
            a_deadcrootn         (:) = spval
            a_deadcrootn_storage (:) = spval
            a_deadcrootn_xfer    (:) = spval
            a_grainn             (:) = spval
            a_grainn_storage     (:) = spval
            a_grainn_xfer        (:) = spval
            a_retransn           (:) = spval
            a_gpp                (:) = spval
            a_downreg            (:) = spval
            a_ar                 (:) = spval
            a_cwdprod            (:) = spval
            a_cwddecomp          (:) = spval
            a_hr                 (:) = spval
            a_fpg                (:) = spval
            a_fpi                (:) = spval
            a_totvegc            (:) = spval
            a_totlitc            (:) = spval
            a_totcwdc            (:) = spval
            a_totsomc            (:) = spval
            a_totcolc            (:) = spval
            a_totvegn            (:) = spval
            a_totlitn            (:) = spval
            a_totcwdn            (:) = spval
            a_totsomn            (:) = spval
            a_totcoln            (:) = spval
            a_gpp_enftemp        (:) = spval
            a_gpp_enfboreal      (:) = spval
            a_gpp_dnfboreal      (:) = spval
            a_gpp_ebftrop        (:) = spval
            a_gpp_ebftemp        (:) = spval
            a_gpp_dbftrop        (:) = spval
            a_gpp_dbftemp        (:) = spval
            a_gpp_dbfboreal      (:) = spval
            a_gpp_ebstemp        (:) = spval
            a_gpp_dbstemp        (:) = spval
            a_gpp_dbsboreal      (:) = spval
            a_gpp_c3arcgrass     (:) = spval
            a_gpp_c3grass        (:) = spval
            a_gpp_c4grass        (:) = spval
            a_npp_enftemp        (:) = spval
            a_npp_enfboreal      (:) = spval
            a_npp_dnfboreal      (:) = spval
            a_npp_ebftrop        (:) = spval
            a_npp_ebftemp        (:) = spval
            a_npp_dbftrop        (:) = spval
            a_npp_dbftemp        (:) = spval
            a_npp_dbfboreal      (:) = spval
            a_npp_ebstemp        (:) = spval
            a_npp_dbstemp        (:) = spval
            a_npp_dbsboreal      (:) = spval
            a_npp_c3arcgrass     (:) = spval
            a_npp_c3grass        (:) = spval
            a_npp_c4grass        (:) = spval
            a_npptoleafc_enftemp        (:) = spval
            a_npptoleafc_enfboreal      (:) = spval
            a_npptoleafc_dnfboreal      (:) = spval
            a_npptoleafc_ebftrop        (:) = spval
            a_npptoleafc_ebftemp        (:) = spval
            a_npptoleafc_dbftrop        (:) = spval
            a_npptoleafc_dbftemp        (:) = spval
            a_npptoleafc_dbfboreal      (:) = spval
            a_npptoleafc_ebstemp        (:) = spval
            a_npptoleafc_dbstemp        (:) = spval
            a_npptoleafc_dbsboreal      (:) = spval
            a_npptoleafc_c3arcgrass     (:) = spval
            a_npptoleafc_c3grass        (:) = spval
            a_npptoleafc_c4grass        (:) = spval
            a_leafc_enftemp      (:) = spval
            a_leafc_enfboreal    (:) = spval
            a_leafc_dnfboreal    (:) = spval
            a_leafc_ebftrop      (:) = spval
            a_leafc_ebftemp      (:) = spval
            a_leafc_dbftrop      (:) = spval
            a_leafc_dbftemp      (:) = spval
            a_leafc_dbfboreal    (:) = spval
            a_leafc_ebstemp      (:) = spval
            a_leafc_dbstemp      (:) = spval
            a_leafc_dbsboreal    (:) = spval
            a_leafc_c3arcgrass   (:) = spval
            a_leafc_c3grass      (:) = spval
            a_leafc_c4grass      (:) = spval

            a_O2_DECOMP_DEPTH_UNSAT (:,:) = spval
            a_CONC_O2_UNSAT         (:,:) = spval
#ifdef CROP
            a_pdcorn             (:) = spval
            a_pdswheat           (:) = spval
            a_pdwwheat           (:) = spval
            a_pdsoybean          (:) = spval
            a_pdcotton           (:) = spval
            a_pdrice1            (:) = spval
            a_pdrice2            (:) = spval
            a_pdsugarcane        (:) = spval
            a_plantdate          (:) = spval
            a_manunitro          (:) = spval
            a_fertnitro_corn     (:) = spval
            a_fertnitro_swheat   (:) = spval
            a_fertnitro_wwheat   (:) = spval
            a_fertnitro_soybean  (:) = spval
            a_fertnitro_cotton   (:) = spval
            a_fertnitro_rice1    (:) = spval
            a_fertnitro_rice2    (:) = spval
            a_fertnitro_sugarcane(:) = spval
            a_irrig_method_corn     (:) = spval
            a_irrig_method_swheat   (:) = spval
            a_irrig_method_wwheat   (:) = spval
            a_irrig_method_soybean  (:) = spval
            a_irrig_method_cotton   (:) = spval
            a_irrig_method_rice1    (:) = spval
            a_irrig_method_rice2    (:) = spval
            a_irrig_method_sugarcane(:) = spval
            a_cphase             (:) = spval
            a_vf                 (:) = spval
            a_gddmaturity        (:) = spval
            a_gddplant           (:) = spval
            a_hui                (:) = spval
            a_cropprod1c         (:) = spval
            a_cropprod1c_loss    (:) = spval
            a_cropseedc_deficit  (:) = spval
            a_grainc_to_cropprodc(:) = spval
            a_grainc_to_seed     (:) = spval
            a_fert_to_sminn      (:) = spval

            a_sum_irrig          (:) = spval
            a_sum_deficit_irrig  (:) = spval
            a_sum_irrig_count    (:) = spval
            a_waterstorage       (:) = spval
            a_groundwater_demand (:) = spval
            a_groundwater_supply (:) = spval
            a_reservoirriver_demand(:) = spval
            a_reservoirriver_supply(:) = spval
            a_reservoir_supply     (:) = spval
            a_river_supply         (:) = spval
            a_runoff_supply        (:) = spval
#endif
            a_ndep_to_sminn      (:) = spval

            a_abm                (:) = spval
            a_gdp                (:) = spval
            a_peatf              (:) = spval
            a_hdm                (:) = spval
            a_lnfm               (:) = spval

            a_leafcCap             (:) = spval
            a_leafc_storageCap     (:) = spval
            a_leafc_xferCap        (:) = spval
            a_frootcCap            (:) = spval
            a_frootc_storageCap    (:) = spval
            a_frootc_xferCap       (:) = spval
            a_livestemcCap         (:) = spval
            a_livestemc_storageCap (:) = spval
            a_livestemc_xferCap    (:) = spval
            a_deadstemcCap         (:) = spval
            a_deadstemc_storageCap (:) = spval
            a_deadstemc_xferCap    (:) = spval
            a_livecrootcCap        (:) = spval
            a_livecrootc_storageCap(:) = spval
            a_livecrootc_xferCap   (:) = spval
            a_deadcrootcCap        (:) = spval
            a_deadcrootc_storageCap(:) = spval
            a_deadcrootc_xferCap   (:) = spval
            a_leafnCap             (:) = spval
            a_leafn_storageCap     (:) = spval
            a_leafn_xferCap        (:) = spval
            a_frootnCap            (:) = spval
            a_frootn_storageCap    (:) = spval
            a_frootn_xferCap       (:) = spval
            a_livestemnCap         (:) = spval
            a_livestemn_storageCap (:) = spval
            a_livestemn_xferCap    (:) = spval
            a_deadstemnCap         (:) = spval
            a_deadstemn_storageCap (:) = spval
            a_deadstemn_xferCap    (:) = spval
            a_livecrootnCap        (:) = spval
            a_livecrootn_storageCap(:) = spval
            a_livecrootn_xferCap   (:) = spval
            a_deadcrootnCap        (:) = spval
            a_deadcrootn_storageCap(:) = spval
            a_deadcrootn_xferCap   (:) = spval
#endif
            a_ozone                (:) = spval

            a_t_soisno     (:,:) = spval
            a_wliq_soisno  (:,:) = spval
            a_wice_soisno  (:,:) = spval
            a_h2osoi       (:,:) = spval
            a_qlayer       (:,:) = spval
            a_lake_deficit   (:) = spval
            a_rootr        (:,:) = spval
            a_BD_all       (:,:) = spval
            a_wfc          (:,:) = spval
            a_OM_density   (:,:) = spval
!Plant Hydraulic parameters
            a_vegwp        (:,:) = spval
!END plant hydraulic parameters
            a_dz_lake      (:,:) = spval
            a_t_lake       (:,:) = spval
            a_lake_icefrac (:,:) = spval
#ifdef BGC
            a_litr1c_vr    (:,:) = spval
            a_litr2c_vr    (:,:) = spval
            a_litr3c_vr    (:,:) = spval
            a_soil1c_vr    (:,:) = spval
            a_soil2c_vr    (:,:) = spval
            a_soil3c_vr    (:,:) = spval
            a_cwdc_vr      (:,:) = spval
            a_litr1n_vr    (:,:) = spval
            a_litr2n_vr    (:,:) = spval
            a_litr3n_vr    (:,:) = spval
            a_soil1n_vr    (:,:) = spval
            a_soil2n_vr    (:,:) = spval
            a_soil3n_vr    (:,:) = spval
            a_cwdn_vr      (:,:) = spval
            a_totsoiln_vr  (:,:) = spval

            a_litr1cCap_vr (:,:) = spval
            a_litr2cCap_vr (:,:) = spval
            a_litr3cCap_vr (:,:) = spval
            a_soil1cCap_vr (:,:) = spval
            a_soil2cCap_vr (:,:) = spval
            a_soil3cCap_vr (:,:) = spval
            a_cwdcCap_vr   (:,:) = spval
            a_litr1nCap_vr (:,:) = spval
            a_litr2nCap_vr (:,:) = spval
            a_litr3nCap_vr (:,:) = spval
            a_soil1nCap_vr (:,:) = spval
            a_soil2nCap_vr (:,:) = spval
            a_soil3nCap_vr (:,:) = spval
            a_cwdnCap_vr   (:,:) = spval

            a_t_scalar     (:,:) = spval
            a_w_scalar     (:,:) = spval

            a_sminn_vr     (:,:) = spval
#endif

            a_ustar (:) = spval
            a_ustar2(:) = spval
            a_tstar (:) = spval
            a_qstar (:) = spval
            a_zol   (:) = spval
            a_rib   (:) = spval
            a_fm    (:) = spval
            a_fh    (:) = spval
            a_fq    (:) = spval

            a_us10m (:) = spval
            a_vs10m (:) = spval
            a_fm10m (:) = spval

            a_sr       (:) = spval
            a_solvd    (:) = spval
            a_solvi    (:) = spval
            a_solnd    (:) = spval
            a_solni    (:) = spval
            a_srvd     (:) = spval
            a_srvi     (:) = spval
            a_srnd     (:) = spval
            a_srni     (:) = spval
            a_solvdln  (:) = spval
            a_solviln  (:) = spval
            a_solndln  (:) = spval
            a_solniln  (:) = spval
            a_srvdln   (:) = spval
            a_srviln   (:) = spval
            a_srndln   (:) = spval
            a_srniln   (:) = spval
#ifdef HYPERSPECTRAL
            a_sol_dir_ln_hires(:,:) = spval
            a_sol_dif_ln_hires(:,:) = spval
            a_sr_dir_ln_hires (:,:) = spval
            a_sr_dif_ln_hires (:,:) = spval
#endif
            a_sensors(:,:) = spval

            nac_ln     (:) = 0
            nac_dt     (:) = 0
            filter_dt  (:) = .true.

         ENDIF
      ENDIF

#ifdef EXTERNAL_LAKE
      CALL Flush_LakeAccVars
#endif

   END SUBROUTINE FLUSH_acc_fluxes

   SUBROUTINE accumulate_fluxes
! ----------------------------------------------------------------------
!  perform the grid average mapping: average a subgrid input 1d vector
!  of length numpatch to a output 2d array of length [ghist%xcnt,ghist%ycnt]
!
!  Created by Yongjiu Dai, 03/2014
!---------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE mod_forcing, only: forcmask_pch
   USE MOD_Mesh,    only: numelm
   USE MOD_LandElm
   USE MOD_LandPatch,      only: numpatch, elm_patch
   USE MOD_LandUrban,      only: numurban
   USE MOD_Const_Physical, only: vonkar, stefnc, cpair, rgas, grav
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_TimeVariables
   USE MOD_Vars_1DForcing
   USE MOD_Vars_1DFluxes
   USE MOD_FrictionVelocity
   USE MOD_Namelist, only: DEF_USE_CBL_HEIGHT, DEF_USE_OZONESTRESS, DEF_USE_PLANTHYDRAULICS, DEF_USE_NITRIF
   USE MOD_TurbulenceLEddy
   USE MOD_Vars_Global
#ifdef CatchLateralFlow
   USE MOD_Catch_Vars_1DFluxes
   USE MOD_Catch_Hist, only: accumulate_fluxes_basin
#endif

   IMPLICIT NONE

   ! Local Variables

   real(r8), allocatable :: r_trad   (:)
   real(r8), allocatable :: r_ustar  (:)
   real(r8), allocatable :: r_ustar2 (:) !define a temporary for estimating us10m only, output should be r_ustar. Shaofeng, 2023.05.20
   real(r8), allocatable :: r_tstar  (:)
   real(r8), allocatable :: r_qstar  (:)
   real(r8), allocatable :: r_zol    (:)
   real(r8), allocatable :: r_rib    (:)
   real(r8), allocatable :: r_fm     (:)
   real(r8), allocatable :: r_fh     (:)
   real(r8), allocatable :: r_fq     (:)

   real(r8), allocatable :: r_us10m  (:)
   real(r8), allocatable :: r_vs10m  (:)
   real(r8), allocatable :: r_fm10m  (:)

   logical,  allocatable :: filter   (:)

   !---------------------------------------------------------------------
   integer  ib, jb, i, j, ielm, istt, iend
   real(r8) sumwt
   real(r8) rhoair,thm,th,thv,ur,displa_av,zldis,hgt_u,hgt_t,hgt_q
   real(r8) hpbl ! atmospheric boundary layer height [m]
   real(r8) z0m_av,z0h_av,z0q_av,us,vs,tm,qm,psrf,taux_e,tauy_e,fsena_e,fevpa_e
   real(r8) r_ustar_e, r_tstar_e, r_qstar_e, r_zol_e, r_ustar2_e, r_fm10m_e
   real(r8) r_fm_e, r_fh_e, r_fq_e, r_rib_e, r_us10m_e, r_vs10m_e
   real(r8) obu,fh2m,fq2m
   real(r8) um,thvstar,beta,zii,wc,wc2

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN

            ! count for time steps
            nac = nac + 1

            ! count for local noon time steps
            DO i = 1, numpatch
               IF (solvdln(i) /= spval) THEN
                  nac_ln(i) = nac_ln(i) + 1
               ENDIF
            ENDDO

            ! set daytime filter
            filter_dt(:) = coszen(:) > 0

            ! count for daytime time steps
            WHERE ( filter_dt(:) ) nac_dt(:) = nac_dt(:) + 1

            CALL acc1d (forc_us    , a_us      )
            CALL acc1d (forc_vs    , a_vs      )
            CALL acc1d (forc_t     , a_t       )
            CALL acc1d (forc_q     , a_q       )
            CALL acc1d (forc_prc   , a_prc     )
            CALL acc1d (forc_prl   , a_prl     )
            CALL acc1d (forc_pbot  , a_pbot    )
            CALL acc1d (forc_frl   , a_frl     )

            CALL acc1d (forc_sols  , a_solarin )
            CALL acc1d (forc_soll  , a_solarin )
            CALL acc1d (forc_solsd , a_solarin )
            CALL acc1d (forc_solld , a_solarin )

            IF (DEF_USE_CBL_HEIGHT) THEN
               CALL acc1d (forc_hpbl , a_hpbl)
            ENDIF

            CALL acc1d (taux    , a_taux    )
            CALL acc1d (tauy    , a_tauy    )
            CALL acc1d (fsena   , a_fsena   )
            CALL acc1d (lfevpa  , a_lfevpa  )
            CALL acc1d (fevpa   , a_fevpa   )
            CALL acc1d (fsenl   , a_fsenl   )
            CALL acc1d (fevpl   , a_fevpl   )
            CALL acc1d (etr     , a_etr     )
            CALL acc1d (fseng   , a_fseng   )
            CALL acc1d (fevpg   , a_fevpg   )
            CALL acc1d (fgrnd   , a_fgrnd   )
            CALL acc1d (sabvsun , a_sabvsun )
            CALL acc1d (sabvsha , a_sabvsha )
            CALL acc1d (sabg    , a_sabg    )
            CALL acc1d (olrg    , a_olrg    )

            IF (DEF_forcing%has_missing_value) THEN
               WHERE (forcmask_pch)
                  rnet = sabg + sabvsun + sabvsha - olrg + forc_frl
               END WHERE
            ELSE
               WHERE(patchmask)
                 rnet = sabg + sabvsun + sabvsha - olrg + forc_frl
               END WHERE
            ENDIF
            CALL acc1d (rnet    , a_rnet    )

            CALL acc1d (xerr    , a_xerr    )
            CALL acc1d (zerr    , a_zerr    )
            CALL acc1d (rsur    , a_rsur    )
#ifndef CatchLateralFlow
            CALL acc1d (rsur_se , a_rsur_se )
            CALL acc1d (rsur_ie , a_rsur_ie )
#endif
            CALL acc1d (rsub          , a_rsub           )
            CALL acc1d (rnof          , a_rnof           )
#ifdef CatchLateralFlow
            CALL acc1d (xwsur         , a_xwsur          )
            CALL acc1d (xwsub         , a_xwsub          )
            CALL acc1d (fldarea       , a_fldarea        )
#endif
            CALL acc1d (qintr         , a_qintr          )
            CALL acc1d (qinfl         , a_qinfl          )
            CALL acc1d (qdrip         , a_qdrip          )

            CALL acc1d (rstfacsun_out , a_rstfacsun      )
            CALL acc1d (rstfacsha_out , a_rstfacsha      )

            CALL acc1d (gssun_out     , a_gssun          )
            CALL acc1d (gssha_out     , a_gssha          )

            CALL acc1d (rss           , a_rss            )
            CALL acc1d (wdsrf         , a_wdsrf          )
            CALL acc1d (zwt           , a_zwt            )
            CALL acc1d (wa            , a_wa             )
            CALL acc1d (wat           , a_wat            )
            CALL acc1d (wetwat        , a_wetwat         )
            CALL acc1d (assim         , a_assim          )
            CALL acc1d (respc         , a_respc          )
            CALL acc1d (assimsun_out  , a_assimsun       )
            CALL acc1d (assimsha_out  , a_assimsha       )
            CALL acc1d (etrsun_out    , a_etrsun         )
            CALL acc1d (etrsha_out    , a_etrsha         )

            CALL acc1d (qcharge       , a_qcharge        )

            CALL acc1d (t_grnd        , a_t_grnd         )
            CALL acc1d (tleaf         , a_tleaf          )
            CALL acc1d (ldew_rain     , a_ldew_rain      )
            CALL acc1d (ldew_snow     , a_ldew_snow      )
            CALL acc1d (ldew          , a_ldew           )
            CALL acc1d (scv           , a_scv            )
            CALL acc1d (snowdp        , a_snowdp         )
            CALL acc1d (fsno          , a_fsno           )
            CALL acc1d (frcsat        , a_frcsat         )
            CALL acc1d (sigf          , a_sigf           )
            CALL acc1d (green         , a_green          )
            CALL acc1d (lai           , a_lai            )
            CALL acc1d (laisun        , a_laisun         )
            CALL acc1d (laisha        , a_laisha         )
            CALL acc1d (sai           , a_sai            )

            ! only acc for daytime for albedo
            CALL acc3d (alb           , a_alb, filter_dt )
#ifdef HYPERSPECTRAL
            CALL acc3d (alb_hires     , a_alb_hires, filter_dt )
#endif

            CALL acc1d (emis          , a_emis           )
            CALL acc1d (z0m           , a_z0m            )

            allocate (r_trad (numpatch)) ; r_trad(:) = spval
            DO i = 1, numpatch
               IF (DEF_forcing%has_missing_value) THEN
                  IF (.not. forcmask_pch(i)) CYCLE
               ENDIF

               IF (.not. patchmask(i)) CYCLE
               r_trad(i) = (olrg(i)/stefnc)**0.25
            ENDDO
            CALL acc1d (r_trad , a_trad )
            deallocate (r_trad          )

            CALL acc1d (tref   , a_tref )
            CALL acc1d (qref   , a_qref )

            ! set 2m WMO temperature
            DO ielm = 1, numelm

               istt = elm_patch%substt(ielm)
               iend = elm_patch%subend(ielm)

               ! landelm%settyp=1 means 2m WMO patch exist,
               ! which is the last end patch in a element.
               IF (landelm%settyp(ielm)==1 .and. forcmask_pch(iend)) THEN
                  ! all set to the 2m WMO patch tref
                  t2m_wmo(istt:iend) = tref(iend)
               ELSE
                  ! if no 2m WMO patch, keep t2m_wmo to tref
                  t2m_wmo(istt:iend) = tref(istt:iend)
               ENDIF
            ENDDO

            CALL acc1d (t2m_wmo, a_t2m_wmo)

            CALL acc1d (forc_rain, a_rain )
            CALL acc1d (forc_snow, a_snow )

            IF (DEF_USE_OZONESTRESS)THEN
               CALL acc1d(o3uptakesun,a_o3uptakesun)
               CALL acc1d(o3uptakesha,a_o3uptakesha)
            ENDIF

#ifdef DataAssimilation
            CALL acc3d (h2osoi_ens     , a_h2osoi_ens     )
            CALL acc3d (t_brt_smap_ens , a_t_brt_smap_ens )
            CALL acc2d (t_brt_smap     , a_t_brt_smap     )
            CALL acc3d (t_brt_fy3d_ens , a_t_brt_fy3d_ens )
            CALL acc2d (t_brt_fy3d     , a_t_brt_fy3d     )
            CALL acc3d (wliq_soisno_ens, a_wliq_soisno_ens)
            CALL acc3d (wice_soisno_ens, a_wice_soisno_ens)
            CALL acc3d (t_soisno_ens   , a_t_soisno_ens   )
#endif

#ifdef URBAN_MODEL
            IF (numurban > 0) THEN
               CALL acc1d(t_room     , a_t_room    )
               CALL acc1d(tafu       , a_tafu      )
               CALL acc1d(fhac       , a_fhac      )
               CALL acc1d(fwst       , a_fwst      )
               CALL acc1d(fach       , a_fach      )
               CALL acc1d(fahe       , a_fahe      )
               CALL acc1d(fhah       , a_fhah      )
               CALL acc1d(vehc       , a_vehc      )
               CALL acc1d(meta       , a_meta      )

               CALL acc1d(fsen_roof  , a_senroof   )
               CALL acc1d(fsen_wsun  , a_senwsun   )
               CALL acc1d(fsen_wsha  , a_senwsha   )
               CALL acc1d(fsen_gimp  , a_sengimp   )
               CALL acc1d(fsen_gper  , a_sengper   )
               CALL acc1d(fsen_urbl  , a_senurbl   )

               CALL acc1d(lfevp_roof , a_lfevproof )
               CALL acc1d(lfevp_gimp , a_lfevpgimp )
               CALL acc1d(lfevp_gper , a_lfevpgper )
               CALL acc1d(lfevp_urbl , a_lfevpurbl )

               CALL acc1d(t_roof     , a_troof     )
               CALL acc1d(t_wall     , a_twall     )
            ENDIF
#endif

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
            CALL acc1d (lai_enftemp        , a_lai_enftemp       )
            CALL acc1d (lai_enfboreal      , a_lai_enfboreal     )
            CALL acc1d (lai_dnfboreal      , a_lai_dnfboreal     )
            CALL acc1d (lai_ebftrop        , a_lai_ebftrop       )
            CALL acc1d (lai_ebftemp        , a_lai_ebftemp       )
            CALL acc1d (lai_dbftrop        , a_lai_dbftrop       )
            CALL acc1d (lai_dbftemp        , a_lai_dbftemp       )
            CALL acc1d (lai_dbfboreal      , a_lai_dbfboreal     )
            CALL acc1d (lai_ebstemp        , a_lai_ebstemp       )
            CALL acc1d (lai_dbstemp        , a_lai_dbstemp       )
            CALL acc1d (lai_dbsboreal      , a_lai_dbsboreal     )
            CALL acc1d (lai_c3arcgrass     , a_lai_c3arcgrass    )
            CALL acc1d (lai_c3grass        , a_lai_c3grass       )
            CALL acc1d (lai_c4grass        , a_lai_c4grass       )
#endif
#ifdef BGC
            CALL acc1d (leafc              , a_leafc               )
            CALL acc1d (leafc_storage      , a_leafc_storage       )
            CALL acc1d (leafc_xfer         , a_leafc_xfer          )
            CALL acc1d (frootc             , a_frootc              )
            CALL acc1d (frootc_storage     , a_frootc_storage      )
            CALL acc1d (frootc_xfer        , a_frootc_xfer         )
            CALL acc1d (livestemc          , a_livestemc           )
            CALL acc1d (livestemc_storage  , a_livestemc_storage   )
            CALL acc1d (livestemc_xfer     , a_livestemc_xfer      )
            CALL acc1d (deadstemc          , a_deadstemc           )
            CALL acc1d (deadstemc_storage  , a_deadstemc_storage   )
            CALL acc1d (deadstemc_xfer     , a_deadstemc_xfer      )
            CALL acc1d (livecrootc         , a_livecrootc          )
            CALL acc1d (livecrootc_storage , a_livecrootc_storage  )
            CALL acc1d (livecrootc_xfer    , a_livecrootc_xfer     )
            CALL acc1d (deadcrootc         , a_deadcrootc          )
            CALL acc1d (deadcrootc_storage , a_deadcrootc_storage  )
            CALL acc1d (deadcrootc_xfer    , a_deadcrootc_xfer     )
            CALL acc1d (grainc             , a_grainc              )
            CALL acc1d (grainc_storage     , a_grainc_storage      )
            CALL acc1d (grainc_xfer        , a_grainc_xfer         )
            CALL acc1d (leafn              , a_leafn               )
            CALL acc1d (leafn_storage      , a_leafn_storage       )
            CALL acc1d (leafn_xfer         , a_leafn_xfer          )
            CALL acc1d (frootn             , a_frootn              )
            CALL acc1d (frootn_storage     , a_frootn_storage      )
            CALL acc1d (frootn_xfer        , a_frootn_xfer         )
            CALL acc1d (livestemn          , a_livestemn           )
            CALL acc1d (livestemn_storage  , a_livestemn_storage   )
            CALL acc1d (livestemn_xfer     , a_livestemn_xfer      )
            CALL acc1d (deadstemn          , a_deadstemn           )
            CALL acc1d (deadstemn_storage  , a_deadstemn_storage   )
            CALL acc1d (deadstemn_xfer     , a_deadstemn_xfer      )
            CALL acc1d (livecrootn         , a_livecrootn          )
            CALL acc1d (livecrootn_storage , a_livecrootn_storage  )
            CALL acc1d (livecrootn_xfer    , a_livecrootn_xfer     )
            CALL acc1d (deadcrootn         , a_deadcrootn          )
            CALL acc1d (deadcrootn_storage , a_deadcrootn_storage  )
            CALL acc1d (deadcrootn_xfer    , a_deadcrootn_xfer     )
            CALL acc1d (grainn             , a_grainn              )
            CALL acc1d (grainn_storage     , a_grainn_storage      )
            CALL acc1d (grainn_xfer        , a_grainn_xfer         )
            CALL acc1d (retransn           , a_retransn            )
            CALL acc1d (gpp                , a_gpp                 )
            CALL acc1d (downreg            , a_downreg             )
            CALL acc1d (ar                 , a_ar                  )
            CALL acc1d (cwdprod            , a_cwdprod             )
            CALL acc1d (cwddecomp          , a_cwddecomp           )
            CALL acc1d (decomp_hr          , a_hr                  )
            CALL acc1d (fpg                , a_fpg                 )
            CALL acc1d (fpi                , a_fpi                 )
            CALL acc1d (totvegc            , a_totvegc             )
            CALL acc1d (totlitc            , a_totlitc             )
            CALL acc1d (totcwdc            , a_totcwdc             )
            CALL acc1d (totsomc            , a_totsomc             )
            CALL acc1d (totcolc            , a_totcolc             )
            CALL acc1d (totvegn            , a_totvegn             )
            CALL acc1d (totlitn            , a_totlitn             )
            CALL acc1d (totcwdn            , a_totcwdn             )
            CALL acc1d (totsomn            , a_totsomn             )
            CALL acc1d (totcoln            , a_totcoln             )
            CALL acc1d (gpp_enftemp        , a_gpp_enftemp         )
            CALL acc1d (gpp_enfboreal      , a_gpp_enfboreal       )
            CALL acc1d (gpp_dnfboreal      , a_gpp_dnfboreal       )
            CALL acc1d (gpp_ebftrop        , a_gpp_ebftrop         )
            CALL acc1d (gpp_ebftemp        , a_gpp_ebftemp         )
            CALL acc1d (gpp_dbftrop        , a_gpp_dbftrop         )
            CALL acc1d (gpp_dbftemp        , a_gpp_dbftemp         )
            CALL acc1d (gpp_dbfboreal      , a_gpp_dbfboreal       )
            CALL acc1d (gpp_ebstemp        , a_gpp_ebstemp         )
            CALL acc1d (gpp_dbstemp        , a_gpp_dbstemp         )
            CALL acc1d (gpp_dbsboreal      , a_gpp_dbsboreal       )
            CALL acc1d (gpp_c3arcgrass     , a_gpp_c3arcgrass      )
            CALL acc1d (gpp_c3grass        , a_gpp_c3grass         )
            CALL acc1d (gpp_c4grass        , a_gpp_c4grass         )
            CALL acc1d (npp_enftemp        , a_npp_enftemp         )
            CALL acc1d (npp_enfboreal      , a_npp_enfboreal       )
            CALL acc1d (npp_dnfboreal      , a_npp_dnfboreal       )
            CALL acc1d (npp_ebftrop        , a_npp_ebftrop         )
            CALL acc1d (npp_ebftemp        , a_npp_ebftemp         )
            CALL acc1d (npp_dbftrop        , a_npp_dbftrop         )
            CALL acc1d (npp_dbftemp        , a_npp_dbftemp         )
            CALL acc1d (npp_dbfboreal      , a_npp_dbfboreal       )
            CALL acc1d (npp_ebstemp        , a_npp_ebstemp         )
            CALL acc1d (npp_dbstemp        , a_npp_dbstemp         )
            CALL acc1d (npp_dbsboreal      , a_npp_dbsboreal       )
            CALL acc1d (npp_c3arcgrass     , a_npp_c3arcgrass      )
            CALL acc1d (npp_c3grass        , a_npp_c3grass         )
            CALL acc1d (npp_c4grass        , a_npp_c4grass         )
            CALL acc1d (npptoleafc_enftemp        , a_npptoleafc_enftemp         )
            CALL acc1d (npptoleafc_enfboreal      , a_npptoleafc_enfboreal       )
            CALL acc1d (npptoleafc_dnfboreal      , a_npptoleafc_dnfboreal       )
            CALL acc1d (npptoleafc_ebftrop        , a_npptoleafc_ebftrop         )
            CALL acc1d (npptoleafc_ebftemp        , a_npptoleafc_ebftemp         )
            CALL acc1d (npptoleafc_dbftrop        , a_npptoleafc_dbftrop         )
            CALL acc1d (npptoleafc_dbftemp        , a_npptoleafc_dbftemp         )
            CALL acc1d (npptoleafc_dbfboreal      , a_npptoleafc_dbfboreal       )
            CALL acc1d (npptoleafc_ebstemp        , a_npptoleafc_ebstemp         )
            CALL acc1d (npptoleafc_dbstemp        , a_npptoleafc_dbstemp         )
            CALL acc1d (npptoleafc_dbsboreal      , a_npptoleafc_dbsboreal       )
            CALL acc1d (npptoleafc_c3arcgrass     , a_npptoleafc_c3arcgrass      )
            CALL acc1d (npptoleafc_c3grass        , a_npptoleafc_c3grass         )
            CALL acc1d (npptoleafc_c4grass        , a_npptoleafc_c4grass         )
            CALL acc1d (leafc_enftemp      , a_leafc_enftemp       )
            CALL acc1d (leafc_enfboreal    , a_leafc_enfboreal     )
            CALL acc1d (leafc_dnfboreal    , a_leafc_dnfboreal     )
            CALL acc1d (leafc_ebftrop      , a_leafc_ebftrop       )
            CALL acc1d (leafc_ebftemp      , a_leafc_ebftemp       )
            CALL acc1d (leafc_dbftrop      , a_leafc_dbftrop       )
            CALL acc1d (leafc_dbftemp      , a_leafc_dbftemp       )
            CALL acc1d (leafc_dbfboreal    , a_leafc_dbfboreal     )
            CALL acc1d (leafc_ebstemp      , a_leafc_ebstemp       )
            CALL acc1d (leafc_dbstemp      , a_leafc_dbstemp       )
            CALL acc1d (leafc_dbsboreal    , a_leafc_dbsboreal     )
            CALL acc1d (leafc_c3arcgrass   , a_leafc_c3arcgrass    )
            CALL acc1d (leafc_c3grass      , a_leafc_c3grass       )
            CALL acc1d (leafc_c4grass      , a_leafc_c4grass       )
            IF(DEF_USE_NITRIF)THEN
               CALL acc2d (to2_decomp_depth_unsat, a_O2_DECOMP_DEPTH_UNSAT)
               CALL acc2d (tconc_o2_unsat        , a_CONC_O2_UNSAT        )
            ENDIF
#ifdef CROP
            CALL acc1d (pdcorn             ,   a_pdcorn             )
            CALL acc1d (pdswheat           ,   a_pdswheat           )
            CALL acc1d (pdwwheat           ,   a_pdwwheat           )
            CALL acc1d (pdsoybean          ,   a_pdsoybean          )
            CALL acc1d (pdcotton           ,   a_pdcotton           )
            CALL acc1d (pdrice1            ,   a_pdrice1            )
            CALL acc1d (pdrice2            ,   a_pdrice2            )
            CALL acc1d (pdsugarcane        ,   a_pdsugarcane        )
            CALL acc1d (plantdate          ,   a_plantdate          )
            CALL acc1d (manunitro          ,   a_manunitro          )
            CALL acc1d (fertnitro_corn     ,   a_fertnitro_corn     )
            CALL acc1d (fertnitro_swheat   ,   a_fertnitro_swheat   )
            CALL acc1d (fertnitro_wwheat   ,   a_fertnitro_wwheat   )
            CALL acc1d (fertnitro_soybean  ,   a_fertnitro_soybean  )
            CALL acc1d (fertnitro_cotton   ,   a_fertnitro_cotton   )
            CALL acc1d (fertnitro_rice1    ,   a_fertnitro_rice1    )
            CALL acc1d (fertnitro_rice2    ,   a_fertnitro_rice2    )
            CALL acc1d (fertnitro_sugarcane,   a_fertnitro_sugarcane)
            CALL acc1d (real(irrig_method_corn     ,r8),   a_irrig_method_corn     )
            CALL acc1d (real(irrig_method_swheat   ,r8),   a_irrig_method_swheat   )
            CALL acc1d (real(irrig_method_wwheat   ,r8),   a_irrig_method_wwheat   )
            CALL acc1d (real(irrig_method_soybean  ,r8),   a_irrig_method_soybean  )
            CALL acc1d (real(irrig_method_cotton   ,r8),   a_irrig_method_cotton   )
            CALL acc1d (real(irrig_method_rice1    ,r8),   a_irrig_method_rice1    )
            CALL acc1d (real(irrig_method_rice2    ,r8),   a_irrig_method_rice2    )
            CALL acc1d (real(irrig_method_sugarcane,r8),   a_irrig_method_sugarcane)
            CALL acc1d (cphase             ,   a_cphase             )
            CALL acc1d (hui                ,   a_hui                )
            CALL acc1d (vf                 ,   a_vf                 )
            CALL acc1d (gddmaturity        ,   a_gddmaturity        )
            CALL acc1d (gddplant           ,   a_gddplant           )
            CALL acc1d (cropprod1c         ,   a_cropprod1c         )
            CALL acc1d (cropprod1c_loss    ,   a_cropprod1c_loss    )
            CALL acc1d (cropseedc_deficit  ,   a_cropseedc_deficit  )
            CALL acc1d (grainc_to_cropprodc,   a_grainc_to_cropprodc)
            CALL acc1d (grainc_to_seed     ,   a_grainc_to_seed     )
            CALL acc1d (fert_to_sminn      ,   a_fert_to_sminn      )

            a_sum_irrig = sum_irrig
            a_sum_irrig_count = sum_irrig_count
            a_waterstorage = waterstorage
            CALL acc1d (groundwater_demand   ,   a_groundwater_demand )
            CALL acc1d (groundwater_supply   ,   a_groundwater_supply )
            CALL acc1d (reservoirriver_demand,  a_reservoirriver_demand)
            CALL acc1d (reservoirriver_supply,  a_reservoirriver_supply)
            CALL acc1d (reservoir_supply     ,  a_reservoir_supply)
            CALL acc1d (river_supply         ,  a_river_supply)
            CALL acc1d (runoff_supply        ,  a_runoff_supply)
#endif
            CALL acc1d (ndep_to_sminn      ,   a_ndep_to_sminn      )
            IF(DEF_USE_FIRE)THEN
               CALL acc1d (abm_lf          ,   a_abm                )
               CALL acc1d (gdp_lf          ,   a_gdp                )
               CALL acc1d (peatf_lf        ,   a_peatf              )
               CALL acc1d (hdm_lf          ,   a_hdm                )
               CALL acc1d (lnfm            ,   a_lnfm               )
            ENDIF
            IF(DEF_USE_DiagMatrix)THEN
               CALL acc1d (leafcCap             ,a_leafcCap             )
               CALL acc1d (leafc_storageCap     ,a_leafc_storageCap     )
               CALL acc1d (leafc_xferCap        ,a_leafc_xferCap        )
               CALL acc1d (frootcCap            ,a_frootcCap            )
               CALL acc1d (frootc_storageCap    ,a_frootc_storageCap    )
               CALL acc1d (frootc_xferCap       ,a_frootc_xferCap       )
               CALL acc1d (livestemcCap         ,a_livestemcCap         )
               CALL acc1d (livestemc_storageCap ,a_livestemc_storageCap )
               CALL acc1d (livestemc_xferCap    ,a_livestemc_xferCap    )
               CALL acc1d (deadstemcCap         ,a_deadstemcCap         )
               CALL acc1d (deadstemc_storageCap ,a_deadstemc_storageCap )
               CALL acc1d (deadstemc_xferCap    ,a_deadstemc_xferCap    )
               CALL acc1d (livecrootcCap        ,a_livecrootcCap        )
               CALL acc1d (livecrootc_storageCap,a_livecrootc_storageCap)
               CALL acc1d (livecrootc_xferCap   ,a_livecrootc_xferCap   )
               CALL acc1d (deadcrootcCap        ,a_deadcrootcCap        )
               CALL acc1d (deadcrootc_storageCap,a_deadcrootc_storageCap)
               CALL acc1d (deadcrootc_xferCap   ,a_deadcrootc_xferCap   )
               CALL acc1d (leafnCap             ,a_leafnCap             )
               CALL acc1d (leafn_storageCap     ,a_leafn_storageCap     )
               CALL acc1d (leafn_xferCap        ,a_leafn_xferCap        )
               CALL acc1d (frootnCap            ,a_frootnCap            )
               CALL acc1d (frootn_storageCap    ,a_frootn_storageCap    )
               CALL acc1d (frootn_xferCap       ,a_frootn_xferCap       )
               CALL acc1d (livestemnCap         ,a_livestemnCap         )
               CALL acc1d (livestemn_storageCap ,a_livestemn_storageCap )
               CALL acc1d (livestemn_xferCap    ,a_livestemn_xferCap    )
               CALL acc1d (deadstemnCap         ,a_deadstemnCap         )
               CALL acc1d (deadstemn_storageCap ,a_deadstemn_storageCap )
               CALL acc1d (deadstemn_xferCap    ,a_deadstemn_xferCap    )
               CALL acc1d (livecrootnCap        ,a_livecrootnCap        )
               CALL acc1d (livecrootn_storageCap,a_livecrootn_storageCap)
               CALL acc1d (livecrootn_xferCap   ,a_livecrootn_xferCap   )
               CALL acc1d (deadcrootnCap        ,a_deadcrootnCap        )
               CALL acc1d (deadcrootn_storageCap,a_deadcrootn_storageCap)
               CALL acc1d (deadcrootn_xferCap   ,a_deadcrootn_xferCap   )
            ENDIF
#endif
            IF(DEF_USE_OZONESTRESS)THEN
               CALL acc1d (forc_ozone  , a_ozone       )
            ENDIF

            IF (.not. DEF_USE_Dynamic_Lake) THEN
               CALL acc1d (lake_deficit, a_lake_deficit)
            ENDIF

            CALL acc2d (t_soisno    , a_t_soisno     )
            CALL acc2d (wliq_soisno , a_wliq_soisno  )
            CALL acc2d (wice_soisno , a_wice_soisno  )

            CALL acc2d (h2osoi      , a_h2osoi       )
            CALL acc2d (qlayer      , a_qlayer       )
            CALL acc2d (rootr       , a_rootr        )
            CALL acc2d (BD_all      , a_BD_all       )
            CALL acc2d (wfc         , a_wfc          )
            CALL acc2d (OM_density  , a_OM_density   )
            IF(DEF_USE_PLANTHYDRAULICS)THEN
               CALL acc2d (vegwp    , a_vegwp        )
            ENDIF
            IF (DEF_USE_Dynamic_Lake) THEN
               CALL acc2d (dz_lake  , a_dz_lake      )
            ENDIF
            CALL acc2d (t_lake      , a_t_lake       )
            CALL acc2d (lake_icefrac, a_lake_icefrac )
#ifdef BGC
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_cpools_vr(j,i_met_lit,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_litr1c_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_cpools_vr(j,i_cel_lit,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_litr2c_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_cpools_vr(j,i_lig_lit,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_litr3c_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_cpools_vr(j,i_soil1,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_soil1c_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_cpools_vr(j,i_soil2,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_soil2c_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_cpools_vr(j,i_soil3,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_soil3c_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_cpools_vr(j,i_cwd,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_cwdc_vr     )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_npools_vr(j,i_met_lit,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_litr1n_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_npools_vr(j,i_cel_lit,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_litr2n_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_npools_vr(j,i_lig_lit,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_litr3n_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_npools_vr(j,i_soil1,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_soil1n_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_npools_vr(j,i_soil2,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_soil2n_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_npools_vr(j,i_soil3,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_soil3n_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_npools_vr(j,i_cwd,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_cwdn_vr     )
            CALL acc2d (totsoiln_vr  , a_totsoiln_vr )

            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_cpools_vr_Cap(j,i_met_lit,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_litr1cCap_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_cpools_vr_Cap(j,i_cel_lit,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_litr2cCap_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_cpools_vr_Cap(j,i_lig_lit,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_litr3cCap_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_cpools_vr_Cap(j,i_soil1,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_soil1cCap_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_cpools_vr_Cap(j,i_soil2,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_soil2cCap_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_cpools_vr_Cap(j,i_soil3,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_soil3cCap_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_cpools_vr_Cap(j,i_cwd,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_cwdcCap_vr     )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_npools_vr_Cap(j,i_met_lit,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_litr1nCap_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_npools_vr_Cap(j,i_cel_lit,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_litr2nCap_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_npools_vr_Cap(j,i_lig_lit,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_litr3nCap_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_npools_vr_Cap(j,i_soil1,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_soil1nCap_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_npools_vr_Cap(j,i_soil2,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_soil2nCap_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_npools_vr_Cap(j,i_soil3,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_soil3nCap_vr   )
            DO i = 1, numpatch
               DO j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_npools_vr_Cap(j,i_cwd,i)
               ENDDO
            ENDDO
            CALL acc2d (decomp_vr_tmp, a_cwdnCap_vr     )
            CALL acc2d (sminn_vr     , a_sminn_vr       )

            CALL acc2d (t_scalar     , a_t_scalar       )
            CALL acc2d (w_scalar     , a_w_scalar       )
#endif
            allocate (r_ustar  (numpatch));  r_ustar (:) = spval
            allocate (r_ustar2 (numpatch));  r_ustar2(:) = spval !Shaofeng, 2023.05.20
            allocate (r_tstar  (numpatch));  r_tstar (:) = spval
            allocate (r_qstar  (numpatch));  r_qstar (:) = spval
            allocate (r_zol    (numpatch));  r_zol   (:) = spval
            allocate (r_rib    (numpatch));  r_rib   (:) = spval
            allocate (r_fm     (numpatch));  r_fm    (:) = spval
            allocate (r_fh     (numpatch));  r_fh    (:) = spval
            allocate (r_fq     (numpatch));  r_fq    (:) = spval
            allocate (r_us10m  (numpatch));  r_us10m (:) = spval
            allocate (r_vs10m  (numpatch));  r_vs10m (:) = spval
            allocate (r_fm10m  (numpatch));  r_fm10m (:) = spval

            DO ielm = 1, numelm

               istt = elm_patch%substt(ielm)
               iend = elm_patch%subend(ielm)

               allocate (filter (istt:iend))
               filter(:) = .true.

               filter(:) = patchmask(istt:iend)

               IF (DEF_forcing%has_missing_value) THEN
                  WHERE (.not. forcmask_pch(istt:iend)) filter = .false.
                  filter = filter .and. forcmask_pch(istt:iend)
               ENDIF

               IF (.not. any(filter)) THEN
                  deallocate(filter)
                  CYCLE
               ENDIF

               sumwt = sum(elm_patch%subfrc(istt:iend), mask = filter)

               ! Aggregate variables from patches to element (gridcell in latitude-longitude mesh)
               z0m_av  = sum(z0m        (istt:iend) * elm_patch%subfrc(istt:iend), mask = filter) / sumwt
               hgt_u   = sum(forc_hgt_u (istt:iend) * elm_patch%subfrc(istt:iend), mask = filter) / sumwt
               hgt_t   = sum(forc_hgt_t (istt:iend) * elm_patch%subfrc(istt:iend), mask = filter) / sumwt
               hgt_q   = sum(forc_hgt_q (istt:iend) * elm_patch%subfrc(istt:iend), mask = filter) / sumwt
               us      = sum(forc_us    (istt:iend) * elm_patch%subfrc(istt:iend), mask = filter) / sumwt
               vs      = sum(forc_vs    (istt:iend) * elm_patch%subfrc(istt:iend), mask = filter) / sumwt
               tm      = sum(forc_t     (istt:iend) * elm_patch%subfrc(istt:iend), mask = filter) / sumwt
               qm      = sum(forc_q     (istt:iend) * elm_patch%subfrc(istt:iend), mask = filter) / sumwt
               psrf    = sum(forc_psrf  (istt:iend) * elm_patch%subfrc(istt:iend), mask = filter) / sumwt
               taux_e  = sum(taux       (istt:iend) * elm_patch%subfrc(istt:iend), mask = filter) / sumwt
               tauy_e  = sum(tauy       (istt:iend) * elm_patch%subfrc(istt:iend), mask = filter) / sumwt
               fsena_e = sum(fsena      (istt:iend) * elm_patch%subfrc(istt:iend), mask = filter) / sumwt
               fevpa_e = sum(fevpa      (istt:iend) * elm_patch%subfrc(istt:iend), mask = filter) / sumwt
               IF (DEF_USE_CBL_HEIGHT) THEN !//TODO: Shaofeng, 2023.05.18
                  hpbl = sum(forc_hpbl(istt:iend) * elm_patch%subfrc(istt:iend), mask = filter) / sumwt
               ENDIF

               z0h_av = z0m_av
               z0q_av = z0m_av

               displa_av = 2./3.*z0m_av/0.07

               hgt_u = max(hgt_u, 5.+displa_av)
               hgt_t = max(hgt_t, 5.+displa_av)
               hgt_q = max(hgt_q, 5.+displa_av)

               zldis = hgt_u-displa_av

               rhoair = (psrf - 0.378*qm*psrf/(0.622+0.378*qm)) / (rgas*tm)

               r_ustar_e = sqrt(max(1.e-6,sqrt(taux_e**2+tauy_e**2))/rhoair)
               r_tstar_e = -fsena_e/(rhoair*r_ustar_e)/cpair
               r_qstar_e = -fevpa_e/(rhoair*r_ustar_e)

               thm = tm + 0.0098*hgt_t
               th  = tm*(100000./psrf)**(rgas/cpair)
               thv = th*(1.+0.61*qm)

               r_zol_e = zldis*vonkar*grav * (r_tstar_e*(1.+0.61*qm)+0.61*th*r_qstar_e) &
                  / (r_ustar_e**2*thv)

               IF(r_zol_e >= 0.)THEN   !stable
                  r_zol_e = min(2.,max(r_zol_e,1.e-6))
               ELSE                       !unstable
                  r_zol_e = max(-100.,min(r_zol_e,-1.e-6))
               ENDIF

               beta = 1.
               zii = 1000.

               thvstar=r_tstar_e*(1.+0.61*qm)+0.61*th*r_qstar_e
               ur = sqrt(us*us+vs*vs)
               IF(r_zol_e >= 0.)THEN
                  um = max(ur,0.1)
               ELSE
                  IF (DEF_USE_CBL_HEIGHT) THEN !//TODO: Shaofeng, 2023.05.18
                     zii = max(5.*hgt_u,hpbl)
                  ENDIF !//TODO: Shaofeng, 2023.05.18
                  wc = (-grav*r_ustar_e*thvstar*zii/thv)**(1./3.)
                  wc2 = beta*beta*(wc*wc)
                  um = max(0.1,sqrt(ur*ur+wc2))
               ENDIF

               obu = zldis/r_zol_e
               IF (DEF_USE_CBL_HEIGHT) THEN
                  CALL moninobuk_leddy(hgt_u,hgt_t,hgt_q,displa_av,z0m_av,z0h_av,z0q_av,&
                     obu,um, hpbl, r_ustar2_e,fh2m,fq2m,r_fm10m_e,r_fm_e,r_fh_e,r_fq_e) !Shaofeng, 2023.05.20
               ELSE
                  CALL moninobuk(hgt_u,hgt_t,hgt_q,displa_av,z0m_av,z0h_av,z0q_av,&
                    obu,um,r_ustar2_e,fh2m,fq2m,r_fm10m_e,r_fm_e,r_fh_e,r_fq_e) !Shaofeng, 2023.05.20
               ENDIF

               ! bug found by chen qiying 2013/07/01
               r_rib_e = r_zol_e /vonkar * r_ustar2_e**2 / (vonkar/r_fh_e*um**2)
               r_rib_e = min(5.,r_rib_e)

               r_us10m_e = us/um * r_ustar2_e /vonkar * r_fm10m_e
               r_vs10m_e = vs/um * r_ustar2_e /vonkar * r_fm10m_e

               ! Assign values from element (gridcell in latitude-longitude mesh) to patches.
               ! Notice that all values on patches in an element are equal.
               r_ustar (istt:iend) = r_ustar_e
               r_ustar2(istt:iend) = r_ustar2_e
               r_tstar (istt:iend) = r_tstar_e
               r_qstar (istt:iend) = r_qstar_e
               r_zol   (istt:iend) = r_zol_e
               r_rib   (istt:iend) = r_rib_e
               r_fm    (istt:iend) = r_fm_e
               r_fh    (istt:iend) = r_fh_e
               r_fq    (istt:iend) = r_fq_e
               r_us10m (istt:iend) = r_us10m_e
               r_vs10m (istt:iend) = r_vs10m_e
               r_fm10m (istt:iend) = r_fm10m_e

               deallocate(filter)

            ENDDO

            CALL acc1d (r_ustar  , a_ustar  )
            CALL acc1d (r_ustar2 , a_ustar2 )
            CALL acc1d (r_tstar  , a_tstar  )
            CALL acc1d (r_qstar  , a_qstar  )
            CALL acc1d (r_zol    , a_zol    )
            CALL acc1d (r_rib    , a_rib    )
            CALL acc1d (r_fm     , a_fm     )
            CALL acc1d (r_fh     , a_fh     )
            CALL acc1d (r_fq     , a_fq     )

            CALL acc1d (r_us10m  , a_us10m  )
            CALL acc1d (r_vs10m  , a_vs10m  )
            CALL acc1d (r_fm10m  , a_fm10m  )

            deallocate (r_ustar )
            deallocate (r_ustar2) !Shaofeng, 2023.05.20
            deallocate (r_tstar )
            deallocate (r_qstar )
            deallocate (r_zol   )
            deallocate (r_rib   )
            deallocate (r_fm    )
            deallocate (r_fh    )
            deallocate (r_fq    )

            deallocate (r_us10m )
            deallocate (r_vs10m )
            deallocate (r_fm10m )

            CALL acc1d (sr      , a_sr      )
            CALL acc1d (solvd   , a_solvd   )
            CALL acc1d (solvi   , a_solvi   )
            CALL acc1d (solnd   , a_solnd   )
            CALL acc1d (solni   , a_solni   )
            CALL acc1d (srvd    , a_srvd    )
            CALL acc1d (srvi    , a_srvi    )
            CALL acc1d (srnd    , a_srnd    )
            CALL acc1d (srni    , a_srni    )
            CALL acc1d (solvdln , a_solvdln )
            CALL acc1d (solviln , a_solviln )
            CALL acc1d (solndln , a_solndln )
            CALL acc1d (solniln , a_solniln )
            CALL acc1d (srvdln  , a_srvdln  )
            CALL acc1d (srviln  , a_srviln  )
            CALL acc1d (srndln  , a_srndln  )
            CALL acc1d (srniln  , a_srniln  )

            CALL acc2d (sensors , a_sensors )
#ifdef HYPERSPECTRAL
            CALL acc2d (sol_dir_ln_hires, a_sol_dir_ln_hires)
            CALL acc2d (sol_dif_ln_hires, a_sol_dif_ln_hires)
            CALL acc2d (sr_dir_ln_hires , a_sr_dir_ln_hires )
            CALL acc2d (sr_dif_ln_hires , a_sr_dif_ln_hires )
#endif

         ENDIF
      ENDIF

#ifdef CatchLateralFlow
      CALL accumulate_fluxes_basin ()
#endif

#ifdef EXTERNAL_LAKE
      CALL accumulate_LakeTimeVars
#endif

   END SUBROUTINE accumulate_fluxes


   !------
   SUBROUTINE acc1d (var, s)

   USE MOD_Precision
   USE MOD_Vars_Global, only: spval

   IMPLICIT NONE

   real(r8), intent(in)    :: var(:)
   real(r8), intent(inout) :: s  (:)
   ! Local variables
   integer :: i

      DO i = lbound(var,1), ubound(var,1)
         IF (var(i) /= spval) THEN
            IF (s(i) /= spval) THEN
               s(i) = s(i) + var(i)
            ELSE
               s(i) = var(i)
            ENDIF
         ENDIF
      ENDDO

   END SUBROUTINE acc1d

   !------
   SUBROUTINE acc2d (var, s)

   USE MOD_Precision
   USE MOD_Vars_Global, only: spval

   IMPLICIT NONE

   real(r8), intent(in)    :: var(:,:)
   real(r8), intent(inout) :: s  (:,:)
   ! Local variables
   integer :: i1, i2

      DO i2 = lbound(var,2), ubound(var,2)
         DO i1 = lbound(var,1), ubound(var,1)
            IF (var(i1,i2) /= spval) THEN
               IF (s(i1,i2) /= spval) THEN
                  s(i1,i2) = s(i1,i2) + var(i1,i2)
               ELSE
                  s(i1,i2) = var(i1,i2)
               ENDIF
            ENDIF
         ENDDO
      ENDDO

   END SUBROUTINE acc2d

   !------
   SUBROUTINE acc3d (var, s, filter)

   USE MOD_Precision
   USE MOD_Vars_Global, only: spval

   IMPLICIT NONE

   real(r8), intent(in)    :: var(:,:,:)
   real(r8), intent(inout) :: s  (:,:,:)
   logical,  intent(in), optional :: filter(:)

   ! Local variables
   integer :: i1, i2, i3

      DO i3 = lbound(var,3), ubound(var,3)

         IF ( present(filter) ) THEN
            IF ( .not. filter(i3) ) CYCLE
         ENDIF

         DO i2 = lbound(var,2), ubound(var,2)
            DO i1 = lbound(var,1), ubound(var,1)
               IF (var(i1,i2,i3) /= spval) THEN
                  IF (s(i1,i2,i3) /= spval) THEN
                     s(i1,i2,i3) = s(i1,i2,i3) + var(i1,i2,i3)
                  ELSE
                     s(i1,i2,i3) = var(i1,i2,i3)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE acc3d

END MODULE MOD_Vars_1DAccFluxes
! ---------- EOP ------------
