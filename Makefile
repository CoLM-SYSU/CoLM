# Makefile for CoLM program

include include/Makeoptions
HEADER = include/define.h

INCLUDE_DIR = -Iinclude -I.bld/ -I${NETCDF_INC}
VPATH = include : share : mksrfdata : mkinidata \
	: main : main/HYDRO : main/BGC : main/URBAN : main/LULCC : main/DA \
	: main/ParaOpt : extends/CaMa/src : postprocess : .bld

# ********** Targets ALL **********
.PHONY: all
all : mkdir_build mksrfdata.x mkinidata.x colm.x postprocess.x lib
	@echo ''
	@echo '*******************************************************'
	@echo '*                                                     *'
	@echo '*        Making all CoLM programs successfully.       *'
	@echo '*                                                     *'
	@echo '*******************************************************'
# ******* End of Targets ALL ******

.PHONY: mkdir_build
mkdir_build :
	mkdir -p .bld

OBJS_SHARED =    \
				  MOD_Precision.o              \
				  MOD_SPMD_Task.o              \
				  MOD_Namelist.o               \
				  MOD_Vars_Global.o            \
				  MOD_Const_Physical.o         \
				  MOD_Const_LC.o               \
				  MOD_Utils.o                  \
				  MOD_IncompleteGamma.o        \
				  MOD_UserDefFun.o             \
				  MOD_TimeManager.o            \
				  MOD_Const_PFT.o              \
				  MOD_NetCDFSerial.o           \
				  MOD_Block.o                  \
				  MOD_Grid.o                   \
				  MOD_Pixel.o                  \
				  MOD_DataType.o               \
				  MOD_NetCDFPoint.o            \
				  MOD_NetCDFBlock.o            \
				  MOD_CatchmentDataReadin.o    \
				  MOD_5x5DataReadin.o          \
				  MOD_Mesh.o                   \
				  MOD_Pixelset.o               \
				  MOD_NetCDFVector.o           \
				  MOD_RangeCheck.o             \
				  MOD_SpatialMapping.o         \
				  MOD_WorkerPushData.o         \
				  MOD_AggregationRequestData.o \
				  MOD_PixelsetShared.o         \
				  MOD_LandElm.o                \
				  MOD_LandHRU.o                \
				  MOD_LandPatch.o              \
				  MOD_Land2mWMO.o              \
				  MOD_LandCrop.o               \
				  MOD_LandPFT.o                \
				  MOD_LandUrban.o              \
				  MOD_Urban_Const_LCZ.o        \
				  MOD_SingleSrfdata.o          \
				  MOD_SrfdataDiag.o            \
				  MOD_SrfdataRestart.o         \
				  MOD_ElmVector.o              \
				  MOD_HRUVector.o              \
				  MOD_MeshFilter.o             \
				  MOD_MT_SurfaceData.o            \
				  MOD_RegionClip.o

${OBJS_SHARED} : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o .bld/$@ $< ${MOD_CMD}.bld

OBJS_SHARED_T = $(addprefix .bld/,${OBJS_SHARED})

include main/MATERIAL/Makefile.material

OBJS_MKSRFDATA = \
				  Aggregation_PercentagesPFT.o           \
				  Aggregation_LAI.o                      \
				  Aggregation_SoilHyperAlbedo.o          \
				  Aggregation_SoilBrightness.o           \
				  Aggregation_LakeDepth.o                \
				  Aggregation_ForestHeight.o             \
				  Aggregation_CanopyStructure.o          \
				  Aggregation_SoilParameters.o           \
				  Aggregation_DBedrock.o                 \
				  Aggregation_Topography.o               \
				  Aggregation_TopoWetness.o              \
				  Aggregation_TopographyFactors.o        \
				  Aggregation_TopographyFactors_Simple.o \
				  Aggregation_Urban.o                    \
				  Aggregation_SoilTexture.o              \
				  MOD_Lulcc_TransferTrace.o              \
				  MOD_MT_Methane_Aggregation.o           \
				  MKSRFDATA.o

MKSRFDATA.o: MOD_MT_Methane_Preprocessing.o

$(OBJS_MKSRFDATA) : %.o : %.F90 ${HEADER} ${OBJS_SHARED} ${OBJS_MT_CONFIG}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o .bld/$@ $< ${MOD_CMD}.bld

OBJS_MKSRFDATA_T = $(addprefix .bld/,${OBJS_MKSRFDATA})

# ------- Target 1: mksrfdata --------
mksrfdata.x : mkdir_build ${HEADER} ${OBJS_SHARED} ${OBJS_MT_CONFIG} ${OBJS_MKSRFDATA}
	@echo ''
	@echo 'making CoLM surface data start >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	@echo ''
	${FF} ${FOPTS} ${OBJS_SHARED_T} ${OBJS_MT_CONFIG_T} ${OBJS_MKSRFDATA_T} -o run/mksrfdata.x ${LDFLAGS}
	@echo ''
	@echo '<<<<<<<<<<<<<<<<<<<<<<<<<< making CoLM surface data completed!'
	@echo ''
# ----- End of Target 1 mksrfdata ----

OBJS_BASIC =    \
				 MOD_Vector_ReadWrite.o          \
				 MOD_dataSpec_PDB.o              \
				 MOD_tav_abs.o                   \
				 MOD_prospect_DB.o               \
				 MOD_Catch_BasinNetwork.o        \
				 MOD_Catch_Vars_TimeVariables.o  \
				 MOD_Catch_Vars_1DFluxes.o       \
				 MOD_Grid_RiverLakePush.o        \
				 MOD_Grid_RiverLakeVectorIO.o    \
				 MOD_Grid_RiverLakeNetwork.o     \
				 MOD_Grid_Reservoir.o            \
				 MOD_Grid_ReservoirUpstream.o    \
				 MOD_Grid_RiverLakeLevee.o       \
				 MOD_Grid_RiverLakeBifurcation.o \
				 MOD_Grid_RiverLakeHistState.o   \
				 MOD_Grid_RiverLakeHistShard.o   \
				 MOD_Grid_RiverLakeHistRoute.o   \
				 MOD_MT_Solute.o                 \
				 MOD_MT_NitrogenN2O.o            \
				 MOD_MT_Lifecycle.o              \
				 MOD_Grid_RiverLakeTimeVars.o    \
				 MOD_MT_RiverLake.o              \
				 MOD_MT_Sediment.o               \
				 MOD_MT_HostHistory.o            \
				 MOD_MT_BGC_NH4.o                \
				 MOD_BGC_Vars_1DFluxes.o         \
				 MOD_BGC_Vars_1DPFTFluxes.o      \
				 MOD_BGC_Vars_PFTimeVariables.o  \
				 MOD_BGC_Vars_TimeInvariants.o   \
				 MOD_BGC_Vars_TimeVariables.o    \
				 MOD_Urban_Vars_1DFluxes.o       \
				 MOD_Urban_Vars_TimeVariables.o  \
				 MOD_Urban_Vars_TimeInvariants.o \
				 MOD_DA_Vars_1DFluxes.o          \
				 MOD_Vars_TimeInvariants.o       \
				 MOD_DA_Vars_TimeVariables.o     \
				 MOD_Vars_TimeVariables.o        \
				 MOD_Vars_1DPFTFluxes.o          \
				 MOD_Vars_1DFluxes.o             \
				 MOD_Vars_1DForcing.o            \
				 MOD_MT_HydroDiagnostics.o      \
				 MOD_MT_ThermalDiagnostics.o     \
				 MOD_MT_CanopyDiagnostics.o      \
				 MOD_MT_LandWork.o              \
				 MOD_MT_SnowRemap.o             \
				 MOD_Hydro_SoilFunction.o        \
				 MOD_Hydro_SoilWater.o           \
				 MOD_Eroot.o                     \
				 MOD_Qsadv.o                     \
				 MOD_UserSpecifiedForcing.o      \
				 MOD_MT_Isotope_Registry.o       \
				 MOD_MT_Isotope_O18.o            \
				 MOD_MT_Isotope_HDO.o            \
				 MOD_MT_Isotope_Lulcc.o          \
				 MOD_MT_Isotope_SoilFlow.o       \
				 MOD_MT_Providers.o              \
				 MOD_MT_Frac.o                   \
				 MOD_MT_EvapLimit.o              \
				 MOD_MT_Vars.o                   \
				 MOD_MT_Conservation.o           \
				 MOD_MT_SoilInit.o               \
				 MOD_MT_ForcingInput.o           \
				 MOD_MT_Forcing.o                \
				 MOD_MT_Precip.o                 \
				 MOD_MT_Evapo.o                  \
				 MOD_MT_SoilWater.o              \
				 MOD_MT_NitrogenLeaching.o       \
				 MOD_MT_Snow.o                   \
				 MOD_MT_HistoryRestart.o         \
				 MOD_MT_Rest.o                   \
				 MOD_LAIEmpirical.o              \
				 MOD_LAIReadin.o                 \
				 MOD_CropReadin.o                \
				 MOD_NitrifData.o                \
				 MOD_NdepData.o                  \
				 MOD_FireData.o                  \
				 MOD_OrbCoszen.o                 \
				 MOD_OrbCosazi.o                 \
				 MOD_HighRes_Parameters.o        \
				 MOD_3DCanopyRadiation.o         \
				 MOD_Aerosol.o                   \
				 MOD_SnowSnicar.o                \
				 MOD_Albedo.o                    \
				 MOD_SnowSnicar_HiRes.o          \
				 MOD_Albedo_HiRes.o              \
				 MOD_SnowFraction.o              \
				 MOD_Urban_LAIReadin.o           \
				 MOD_Urban_Shortwave.o           \
				 MOD_Urban_Albedo.o              \
				 MOD_MonthlyinSituCO2MaunaLoa.o  \
				 MOD_PercentagesPFTReadin.o      \
				 MOD_LakeDepthReadin.o           \
				 MOD_DBedrockReadin.o            \
				 MOD_SoilColorRefl.o             \
				 MOD_SoilParametersReadin.o      \
				 MOD_SoilTextureReadin.o         \
				 MOD_HtopReadin.o                \
				 MOD_UrbanReadin.o               \
				 MOD_BGC_CNSummary.o             \
				 MOD_IniTimeVariable.o           \
				 MOD_UrbanIniTimeVariable.o      \
				 MOD_ElementNeighbour.o          \
				 MOD_Catch_HillslopeNetwork.o    \
				 MOD_Catch_RiverLakeNetwork.o    \
				 MOD_Catch_Reservoir.o           \
				 MOD_VicParaReadin.o             \
				 MOD_Initialize.o


$(OBJS_BASIC) : %.o : %.F90 ${HEADER} ${OBJS_SHARED} ${OBJS_MT_CONFIG}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o .bld/$@ $< ${MOD_CMD}.bld

OBJS_BASIC_T = $(addprefix .bld/,${OBJS_BASIC})

OBJS_MKINIDATA = \
				  MOD_UnitCatchmentSubset.o              \
				  MOD_UnitCatchmentRegional.o            \
				  CoLMINI.o

MOD_UnitCatchmentRegional.o: MOD_UnitCatchmentSubset.o
CoLMINI.o: MOD_UnitCatchmentRegional.o

$(OBJS_MKINIDATA) : %.o : %.F90 ${HEADER} ${OBJS_SHARED} ${OBJS_MT_CONFIG} ${OBJS_BASIC}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o .bld/$@ $< ${MOD_CMD}.bld

OBJS_MKINIDATA_T = $(addprefix .bld/,${OBJS_MKINIDATA})

# -------- Target 2: mkinidata -------
mkinidata.x : mkdir_build ${HEADER} ${OBJS_SHARED} ${OBJS_MT_CONFIG} ${OBJS_BASIC} ${OBJS_MKINIDATA}
	@echo ''
	@echo 'making CoLM initial data start >>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	@echo ''
	${FF} ${FOPTS} ${OBJS_SHARED_T} ${OBJS_MT_CONFIG_T} ${OBJS_BASIC_T} ${OBJS_MKINIDATA_T} -o run/mkinidata.x ${LDFLAGS}
	@echo ''
	@echo '<<<<<<<<<<<<<<<<<<<<<<<<< making CoLM initial data completed!'
	@echo ''
# ----- End of Target 2 mkinidata ----

CaMa = $(shell \
	CAMA_INIT=$$(grep -E '^\#(define|undef)[[:space:]]+CaMa_Flood' $(HEADER) | head -1 | grep -q 'define' && echo 1 || echo 0); \
	SP=$$(grep -E '^\#(define|undef)[[:space:]]+SinglePoint' $(HEADER) | tail -1 | grep -q 'define' && echo 1 || echo 0); \
	MP=$$(if [ $$SP -eq 1 ]; then echo 0; else grep -E '^\#(define|undef)[[:space:]]+USEMPI' $(HEADER) | head -1 | grep -q 'define' && echo 1 || echo 0; fi); \
	if [ $$CAMA_INIT -eq 0 ]; then echo NO; elif [ $$SP -eq 1 ] || [ $$MP -eq 0 ]; then echo NO; else echo YES; fi)
ifeq (${CaMa},YES)# Compile CoLM decoupled with river routing scheme (CaMa-Flood)
#
OBJECTS_CAMA=\
				  parkind1.o              \
				  yos_cmf_input.o         \
				  yos_cmf_time.o          \
				  yos_cmf_map.o           \
				  yos_cmf_prog.o          \
				  yos_cmf_diag.o          \
				  cmf_utils_mod.o         \
				  cmf_coupling_budget_mod.o \
				  cmf_calc_outflw_mod.o   \
				  cmf_calc_pthout_mod.o   \
				  cmf_calc_fldstg_mod.o   \
				  cmf_calc_stonxt_mod.o   \
				  cmf_opt_outflw_mod.o    \
				  cmf_ctrl_tracer_mod.o   \
				  cmf_ctrl_mpi_mod.o      \
				  cmf_ctrl_damout_mod.o   \
				  cmf_ctrl_levee_mod.o    \
				  cmf_ctrl_forcing_mod.o  \
				  cmf_ctrl_boundary_mod.o \
				  cmf_ctrl_output_mod.o   \
				  cmf_ctrl_restart_mod.o  \
				  cmf_ctrl_sed_mod.o      \
				  cmf_calc_diag_mod.o     \
				  cmf_ctrl_physics_mod.o  \
				  cmf_ctrl_time_mod.o     \
				  cmf_ctrl_maps_mod.o     \
				  cmf_ctrl_vars_mod.o     \
				  cmf_ctrl_nmlist_mod.o   \
				  cmf_drv_control_mod.o   \
				  cmf_drv_advance_mod.o

cmf_coupling_budget_mod.o: parkind1.o

$(OBJECTS_CAMA) : %.o : %.F90 ${HEADER}
	$(FCMP)  -c ${FFLAGS} $(MODS) ${CFLAGS} $(INCLUDE_DIR) -o .bld/$@ $< ${MOD_CMD}.bld

OBJS_CAMA_T = $(addprefix .bld/,${OBJECTS_CAMA})

endif

OBJS_MAIN = \
				MOD_Catch_HillslopeFlow.o                 \
				MOD_Catch_SubsurfaceFlow.o                \
				MOD_Catch_RiverLakeFlow.o                 \
				MOD_Catch_Hist.o                          \
				MOD_Catch_WriteParameters.o               \
				MOD_BGC_CNCStateUpdate1.o                 \
				MOD_BGC_CNCStateUpdate2.o                 \
				MOD_BGC_CNCStateUpdate3.o                 \
				MOD_BGC_CNNStateUpdate1.o                 \
				MOD_BGC_CNNStateUpdate2.o                 \
				MOD_BGC_CNNStateUpdate3.o                 \
				MOD_BGC_Soil_BiogeochemNStateUpdate1.o    \
				MOD_BGC_Soil_BiogeochemNitrifDenitrif.o   \
				MOD_BGC_Soil_BiogeochemCompetition.o      \
				MOD_BGC_Soil_BiogeochemDecompCascadeBGC.o \
				MOD_BGC_Soil_BiogeochemDecomp.o           \
				MOD_BGC_Soil_BiogeochemLittVertTransp.o   \
				MOD_BGC_Soil_BiogeochemNLeaching.o        \
				MOD_BGC_Soil_BiogeochemPotential.o        \
				MOD_BGC_Soil_BiogeochemVerticalProfile.o  \
				MOD_BGC_Veg_CNGapMortality.o              \
				MOD_BGC_Veg_CNGResp.o                     \
				MOD_BGC_Veg_CNMResp.o                     \
				MOD_BGC_Daylength.o                       \
				MOD_BGC_Veg_CNPhenology.o                 \
				MOD_BGC_Veg_NutrientCompetition.o         \
				MOD_BGC_Veg_CNVegStructUpdate.o           \
				MOD_BGC_CNAnnualUpdate.o                  \
				MOD_BGC_CNZeroFluxes.o                    \
				MOD_BGC_CNBalanceCheck.o                  \
				MOD_BGC_CNSASU.o                          \
				MOD_BGC_Veg_CNNDynamics.o                 \
				MOD_BGC_Veg_CNFireBase.o                  \
				MOD_BGC_Veg_CNFireLi2016.o                \
				MOD_Vars_2DForcing.o                      \
				MOD_ForcingDownscaling.o                  \
				MOD_Forcing.o                             \
				MOD_DA_TWS.o                              \
				MOD_DA_Const.o                            \
				MOD_DA_RTM.o                              \
				MOD_DA_EnKF.o                             \
				MOD_DA_SM.o                               \
				MOD_DA_Ensemble.o                         \
				MOD_DA_Main.o                             \
				MOD_Opt_Baseflow.o                        \
				MOD_ParameterOptimization.o               \
				MOD_AssimStomataConductance.o             \
				MOD_PlantHydraulic.o                      \
				MOD_FrictionVelocity.o                    \
				MOD_TurbulenceLEddy.o                     \
				MOD_Ozone.o                               \
				MOD_CanopyLayerProfile.o                  \
				MOD_LeafTemperature.o                     \
				MOD_LeafTemperaturePC.o                   \
				MOD_SoilThermalParameters.o               \
				MOD_Hydro_VIC_Variables.o                 \
				MOD_Hydro_VIC.o                           \
				MOD_Runoff.o                              \
				MOD_SoilSnowHydrology.o                   \
				MOD_SnowLayersCombineDivide.o             \
				MOD_MT_SnowHost.o                         \
				MOD_PhaseChange.o                         \
				MOD_Glacier.o                             \
				MOD_Lake.o                                \
				MOD_SimpleOcean.o                         \
				MOD_GroundFluxes.o                        \
				MOD_GroundTemperature.o                   \
				MOD_LeafInterception.o                    \
				MOD_NetSolar.o                            \
				MOD_NetSolar_Hyper.o                      \
				MOD_WetBulb.o                             \
				MOD_RainSnowTemp.o                        \
				MOD_SoilSurfaceResistance.o               \
				MOD_NewSnow.o                             \
				MOD_CaMa_Vars.o                           \
				MOD_CaMa_colmCaMa.o                       \
				MOD_Thermal.o                             \
				MOD_MT_LandPhase.o                        \
				MOD_Vars_1DAccFluxes.o                    \
				MOD_Irrigation.o                          \
				MOD_BGC_driver.o                          \
				MOD_HistWriteBack.o                       \
				MOD_HistGridded.o                         \
				MOD_HistVector.o                          \
				MOD_HistSingle.o                          \
				MOD_MT_Hist.o                             \
				MOD_MT_Methane_VegOverride.o              \
				MOD_MT_Methane_GIEMS.o                    \
				MOD_MT_Methane_pH.o                       \
				MOD_MT_Methane_State.o                    \
				MOD_MT_Methane_Microbes.o                 \
				MOD_MT_Methane_WetlandBGC.o               \
				MOD_MT_Methane_BgcLink.o                  \
				MOD_MT_Methane_AccFlux.o                  \
				MOD_MT_Methane_Physics.o                  \
				MOD_MT_Methane_Driver.o                   \
				MOD_MT_Methane_Hist.o                     \
				MOD_MT_Methane_Impl.o                     \
				MOD_MT_Methane.o                          \
				MT_RegisterProviders.o                    \
				MOD_MT_HostFixes.o                        \
				MOD_MT_SpecialPatches.o                   \
				MOD_Grid_RiverLakeHist.o                  \
				MOD_HistAccRestart.o                      \
				MOD_Hist.o                                \
				MOD_CheckEquilibrium.o                    \
				MOD_LightningData.o                       \
				MOD_Catch_LateralFlow.o                   \
				MOD_MT_RiverLakeCoupling.o                \
				MOD_Grid_RiverLakeFlow.o                  \
				MOD_Urban_Longwave.o                      \
				MOD_Urban_NetSolar.o                      \
				MOD_Urban_Flux.o                          \
				MOD_Urban_GroundFlux.o                    \
				MOD_Urban_RoofFlux.o                      \
				MOD_Urban_RoofTemperature.o               \
				MOD_Urban_WallTemperature.o               \
				MOD_Urban_PerviousTemperature.o           \
				MOD_Urban_ImperviousTemperature.o         \
				MOD_Urban_Hydrology.o                     \
				MOD_Urban_BEM.o                           \
				MOD_Urban_LUCY.o                          \
				MOD_Urban_Thermal.o                       \
				CoLMMAIN_Urban.o                          \
				MOD_Lulcc_Vars_TimeInvariants.o           \
				MOD_Lulcc_Vars_TimeVariables.o            \
				MOD_Lulcc_TransferTraceReadin.o           \
				MOD_Lulcc_MassEnergyConserve.o            \
				MOD_Lulcc_Initialize.o                    \
				MOD_MT_Lulcc.o                            \
				MOD_Lulcc_Driver.o                        \
				MOD_MT_Runtime.o                          \
				CoLMDRIVER.o                              \
				CoLMMAIN.o                                \
				CoLM.o

$(OBJS_MAIN) : %.o : %.F90 ${HEADER} ${OBJS_SHARED} ${OBJS_BASIC}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o .bld/$@ $< ${MOD_CMD}.bld

MOD_Urban_Thermal.o: MOD_Urban_Flux.o
MOD_LeafTemperature.o MOD_LeafTemperaturePC.o: MOD_LeafInterception.o

OBJS_MAIN_T = $(addprefix .bld/,${OBJS_MAIN})

# ------ Target 3: main --------

ifneq (${CaMa},YES)# Compile CoLM decoupled without river routing scheme (CaMa-Flood)

colm.x : mkdir_build ${HEADER} ${OBJS_SHARED} ${OBJS_MT_CONFIG} ${OBJS_BASIC} ${OBJS_MAIN}
	@echo ''
	@echo 'making CoLM start >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	@echo ''
	${FF} ${FOPTS} ${OBJS_SHARED_T} ${OBJS_MT_CONFIG_T} ${OBJS_BASIC_T} ${OBJS_MAIN_T} -o run/colm.x ${LDFLAGS}
	@echo ''
	@echo '<<<<<<<<<<<<<<<<<<<<<<<<<<<<< making CoLM completed!'
	@echo ''

else
colm.x : mkdir_build  ${HEADER} ${OBJS_SHARED} ${OBJS_MT_CONFIG} ${OBJECTS_CAMA} ${OBJS_BASIC} ${OBJS_MAIN}
	@echo ''
	@echo 'making CoLM with CaMa start >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	@echo ''
	${FF} ${FOPTS} ${OBJS_SHARED_T} ${OBJS_MT_CONFIG_T} ${OBJS_BASIC_T} ${OBJS_CAMA_T} ${OBJS_MAIN_T} -o run/colm.x ${LDFLAGS}

	@echo ''
	@echo '<<<<<<<<<<<<<<<<<<<<<<<<<<<< making CoLM with CaMa completed!'
	@echo ''

endif

# ----- End of Target 3 main -----

OBJS_POST1 = MOD_Concatenate.o HistConcatenate.o
OBJS_POST2 = MOD_Vector2Grid.o POST_Vector2Grid.o
OBJS_POST3 = SrfDataConcatenate.o
OBJS_POST4 = RiverHistConcatenate.o
OBJS_POST1_T = $(addprefix .bld/,${OBJS_POST1})
OBJS_POST2_T = $(addprefix .bld/,${OBJS_POST2})
OBJS_POST3_T = $(addprefix .bld/,${OBJS_POST3})
OBJS_POST4_T = $(addprefix .bld/,${OBJS_POST4})

$(OBJS_POST1):%.o:%.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o .bld/$@ $< ${MOD_CMD}.bld

$(OBJS_POST2):%.o:%.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o .bld/$@ $< ${MOD_CMD}.bld

$(OBJS_POST3):%.o:%.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o .bld/$@ $< ${MOD_CMD}.bld

$(OBJS_POST4):%.o:%.F90 ${HEADER} ${OBJS_SHARED}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o .bld/$@ $< ${MOD_CMD}.bld

hist_concatenate.x : ${HEADER} ${OBJS_SHARED} ${OBJS_POST1}
	${FF} ${FOPTS} ${OBJS_SHARED_T} ${OBJS_POST1_T} -o run/$@ ${LDFLAGS}

post_vector2grid.x : ${HEADER} ${OBJS_SHARED} ${OBJS_POST2}
	${FF} ${FOPTS} ${OBJS_SHARED_T} ${OBJS_POST2_T} -o run/$@ ${LDFLAGS}

srfdata_concatenate.x : ${HEADER} ${OBJS_SHARED} ${OBJS_POST3}
	${FF} ${FOPTS} ${OBJS_SHARED_T} ${OBJS_POST3_T} -o run/$@ ${LDFLAGS}

river_hist_concatenate.x : ${HEADER} ${OBJS_SHARED} ${OBJS_POST4}
	${FF} ${FOPTS} ${OBJS_SHARED_T} ${OBJS_POST4_T} -o run/$@ ${LDFLAGS}

# ------ Target 4: postprocess --------
DEF = $(shell grep -i CATCHMENT include/define.h)
vector2grid = $(word 1, ${DEF})
ifneq (${vector2grid},\#define)
DEF = $(shell grep -i UNSTRUCTURED include/define.h)
vector2grid = $(word 1, ${DEF})
endif

.PHONY: postprocess.x
ifneq (${vector2grid},\#define)
postprocess.x : mkdir_build hist_concatenate.x river_hist_concatenate.x srfdata_concatenate.x
	@echo '<<<<<<<<<<<<<<<<<<<<<<<<< making CoLM postprocessing completed!'
else
postprocess.x : mkdir_build hist_concatenate.x river_hist_concatenate.x srfdata_concatenate.x post_vector2grid.x
	@echo '<<<<<<<<<<<<<<<<<<<<<<<<< making CoLM postprocessing completed!'
endif
# --- End of Target 4 postprocess ------

# ------ Target 5: static libs --------
.PHONY: lib
lib :
	@echo ''
	@echo 'making CoLM static library >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	mkdir -p lib
	cd lib && find ../.bld -name "*.o" ! -name "CoLM.o" ! -name "MKSRFDATA.o" ! -name "CoLMINI.o" -exec ln -sf {} ./ \;
	cd lib && ar rc libcolm.a *.o && ranlib libcolm.a
	ln -sf lib/libcolm.a ./libcolm.a
# ------End of Target 5: static libs --------

.PHONY: clean
clean :
	rm -rf .bld
	rm -rf lib libcolm.a
	rm -f run/mksrfdata.x run/mkinidata.x run/colm.x
	rm -f run/hist_concatenate.x run/river_hist_concatenate.x run/srfdata_concatenate.x run/post_vector2grid.x
	rm -f CaMa/src/*.o CaMa/src/*.mod CaMa/src/*.a

