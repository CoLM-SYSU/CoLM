#include <define.h>
#if (defined MaterialTransport) && (defined BGC)
MODULE MOD_MT_Methane_Driver
!=======================================================================
! Driver wrapper for Methane reactive-material physics.
!=======================================================================
	USE MOD_Precision, only: r8
	USE MOD_Vars_Global, only: nl_soil

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: methane_driver

	TYPE :: methane_column_result_type
		real(r8) :: net, surf_flux, surf_flux_phys, surf_aere, surf_ebul, surf_diff, surf_diff_phys
		real(r8) :: balance_residual, ch4_clip_credit, o2_cap_loss, o2_cap_gain
		real(r8) :: ebul_tot, prod_tot, oxid_tot, co2_decomp_tot, co2_oxid_tot, co2_aere_tot, co2_net_tot
		real(r8) :: totcol, grnd_cond, dfsat_tot, finundated, finundated_default
		real(r8) :: net_unsat, net_sat
		real(r8) :: surf_flux_unsat, surf_flux_sat, surf_aere_unsat, surf_aere_sat
		real(r8) :: surf_ebul_unsat, surf_ebul_sat, surf_diff_unsat, surf_diff_sat
		real(r8) :: surf_diff_phys_unsat, surf_diff_phys_sat
		real(r8) :: ebul_tot_unsat, ebul_tot_sat, prod_tot_unsat, prod_tot_sat
		real(r8) :: oxid_tot_unsat, oxid_tot_sat
		real(r8) :: co2_decomp_tot_unsat, co2_decomp_tot_sat
		real(r8) :: co2_oxid_tot_unsat, co2_oxid_tot_sat, co2_net_tot_unsat, co2_net_tot_sat
		real(r8) :: totcol_unsat, totcol_sat, grnd_cond_unsat, grnd_cond_sat
		real(r8) :: prod(nl_soil), o2_decomp(nl_soil), co2_decomp(nl_soil)
		real(r8) :: oxid(nl_soil), o2_oxid(nl_soil), co2_oxid(nl_soil)
		real(r8) :: aere(nl_soil), tran(nl_soil), o2_aere(nl_soil), co2_aere(nl_soil)
		real(r8) :: ebul(nl_soil), o2stress(nl_soil), ch4stress(nl_soil)
		real(r8) :: prod_unsat(nl_soil), prod_sat(nl_soil)
		real(r8) :: o2_decomp_unsat(nl_soil), o2_decomp_sat(nl_soil)
		real(r8) :: co2_decomp_unsat(nl_soil), co2_decomp_sat(nl_soil)
		real(r8) :: oxid_unsat(nl_soil), oxid_sat(nl_soil)
		real(r8) :: o2_oxid_unsat(nl_soil), o2_oxid_sat(nl_soil)
		real(r8) :: co2_oxid_unsat(nl_soil), co2_oxid_sat(nl_soil)
		real(r8) :: aere_unsat(nl_soil), aere_sat(nl_soil)
		real(r8) :: tran_unsat(nl_soil), tran_sat(nl_soil)
		real(r8) :: o2_aere_unsat(nl_soil), o2_aere_sat(nl_soil)
		real(r8) :: co2_aere_unsat(nl_soil), co2_aere_sat(nl_soil)
		real(r8) :: ebul_unsat(nl_soil), ebul_sat(nl_soil)
		real(r8) :: o2stress_unsat(nl_soil), o2stress_sat(nl_soil)
		real(r8) :: ch4stress_unsat(nl_soil), ch4stress_sat(nl_soil)
		real(r8) :: conc_o2(nl_soil), conc_ch4(nl_soil)
		real(r8) :: conc_o2_unsat(nl_soil), conc_o2_sat(nl_soil)
		real(r8) :: conc_ch4_unsat(nl_soil), conc_ch4_sat(nl_soil)
	END TYPE methane_column_result_type

CONTAINS

	SUBROUTINE methane_driver (istep,i,idate,patchclass,patchtype,deltim,lb,snl,dlon,dlat,&!input
		z_soisno,dz_soisno,zi_soisno,t_soisno,t_grnd,wliq_soisno,wice_soisno,&
			forc_t,forc_pbot,forc_po2m,forc_pco2m,forc_us,forc_vs,ustar,fq,&
		zwt,rootfr,snowdp,wat,rsur,etr,lakedepth,dz_lake,t_lake,lake_icefrac,wdsrf,wetwat,bsw,&
		smp,porsl,lai,sai,rootr,fsatmax,fsatdcf,frcsat,f_h2osfc, &
		is_rice_paddy_in, rice_pft_frac_in)

		use MOD_Precision
		use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
		use MOD_Const_Physical, only: rgas, denh2o, denice, tfrz, grav
		use MOD_MT_Methane_Const
		use MOD_Namelist, only : DEF_USE_VariablySaturatedFlow
		use MOD_Vars_Global, only : maxsnl,nl_soil,nl_lake,spval,PI,deg2rad,z_soi,zi_soi,dz_soi
		use MOD_MT_Methane_Physics
		use MOD_SPMD_Task
			USE MOD_MT_Methane_BgcLink, only: mt_ch4_bgc_patch_inputs, &
			     mt_ch4_bgc_component_veg_inputs, &
			     get_wetland_veg_proxy, get_rice_veg_proxy, is_paddy_rice_live, &
			     rice_days_since_harvest, organic_max, get_biome_f_methane, &
			     get_biome_redoxlag, mt_ch4_bgc_finalize_step
			USE MOD_MT_Methane_VegOverride, only: wetland_aere_active
		USE MOD_MT_Methane_Microbes, only: methane_microbes_step, &
		     aggregate_methane_microbes, repartition_methane_microbes, &
		     reset_methane_inactive_lake_microbe_diagnostics, &
		     microbial_prod_potential, microbial_oxid_potential, &
		     microbial_prod_potential_comp, microbial_oxid_potential_comp
		! ----- Methane material state (module-level allocatables) -----
		USE MOD_MT_Methane_State, only: net_methane, co2_hr_methane, methane_prod_depth, o2_decomp_depth, &
		     co2_decomp_depth, methane_oxid_depth, o2_oxid_depth, co2_oxid_depth, &
		     methane_aere_depth, methane_tran_depth, o2_aere_depth, co2_aere_depth, methane_ebul_depth, &
		     o2stress, methane_stress, &
		     methane_surf_flux_tot, methane_surf_flux_tot_phys, methane_surf_aere, methane_surf_ebul, methane_surf_diff, &
		     methane_surf_diff_phys, &
		     methane_balance_residual, methane_ch4_clip_credit, o2_cap_loss, o2_cap_gain, &
		     methane_ebul_tot, methane_prod_tot, methane_oxid_tot, &
		     co2_decomp_tot, co2_oxid_tot, co2_aere_tot, co2_net_tot, &
		     totcol_methane, grnd_methane_cond, conc_o2, conc_methane, &
		     net_methane_unsat, net_methane_sat, &
		     methane_prod_depth_unsat, methane_prod_depth_sat, &
		     o2_decomp_depth_unsat, o2_decomp_depth_sat, &
		     co2_decomp_depth_unsat, co2_decomp_depth_sat, &
		     methane_oxid_depth_unsat, methane_oxid_depth_sat, &
		     o2_oxid_depth_unsat, o2_oxid_depth_sat, &
		     co2_oxid_depth_unsat, co2_oxid_depth_sat, &
		     methane_aere_depth_unsat, methane_aere_depth_sat, &
		     methane_tran_depth_unsat, methane_tran_depth_sat, &
		     o2_aere_depth_unsat, o2_aere_depth_sat, &
		     co2_aere_depth_unsat, co2_aere_depth_sat, &
		     methane_ebul_depth_unsat, methane_ebul_depth_sat, &
		     o2stress_unsat, o2stress_sat, methane_stress_unsat, methane_stress_sat, &
		     methane_surf_flux_tot_unsat, methane_surf_flux_tot_sat, &
		     methane_surf_aere_unsat, methane_surf_aere_sat, &
		     methane_surf_ebul_unsat, methane_surf_ebul_sat, &
		     methane_surf_diff_unsat, methane_surf_diff_sat, &
		     methane_surf_diff_phys_unsat, methane_surf_diff_phys_sat, &
		     methane_ebul_tot_unsat, methane_ebul_tot_sat, &
		     methane_prod_tot_unsat, methane_prod_tot_sat, &
		     methane_oxid_tot_unsat, methane_oxid_tot_sat, &
		     co2_decomp_tot_unsat, co2_decomp_tot_sat, &
		     co2_oxid_tot_unsat, co2_oxid_tot_sat, &
		     co2_net_tot_unsat, co2_net_tot_sat, &
		     totcol_methane_unsat, totcol_methane_sat, &
		     grnd_methane_cond_unsat, grnd_methane_cond_sat, &
		     conc_o2_unsat, conc_o2_sat, conc_methane_unsat, conc_methane_sat, &
		     methane_prod_depth_lake, methane_oxid_depth_lake, methane_ebul_depth_lake, &
		     co2_decomp_depth_lake, co2_oxid_depth_lake, &
		     methane_surf_ebul_lake, methane_surf_diff_lake, methane_surf_flux_tot_lake, &
		     methane_prod_tot_lake, methane_oxid_tot_lake, methane_ebul_tot_lake, &
		     co2_decomp_tot_lake, co2_oxid_tot_lake, co2_net_tot_lake, &
		     totcol_methane_lake, grnd_methane_cond_lake, conc_o2_lake, conc_methane_lake, &
		     lake_water_ch4_stock, lake_water_o2_stock, lake_frozen_ch4_stock, lake_frozen_o2_stock, &
		     lake_liquid_fraction_prev, lake_water_ch4_oxid, &
		     lake_sed_ch4_flux, lake_sed_o2_flux, lake_air_o2_flux, &
		     c_atm, forc_pmethanem, layer_sat_lag, lake_soilc, &
		     annavg_agnpp, annavg_bgnpp, annavg_somhr, annavg_finrw, &
		     tempavg_agnpp, tempavg_bgnpp, annsum_counter, &
		     tempavg_somhr, tempavg_finrw, &
		     fsat_bef, finundated_lag, methane_dfsat_tot, &
		     biome_f_methane_patch, biome_redoxlag_patch, &
		     f_inund_flood_patch, wetland_frac_per_patch, &
		     conc_o2_unsat_component, conc_o2_sat_component, &
		     conc_methane_unsat_component, conc_methane_sat_component, &
		     layer_sat_lag_component, annavg_agnpp_component, annavg_bgnpp_component, &
		     annavg_somhr_component, annavg_finrw_component, tempavg_agnpp_component, &
		     tempavg_bgnpp_component, annsum_counter_component, tempavg_somhr_component, &
		     tempavg_finrw_component, fsat_bef_component, finundated_lag_component, &
		     rice_fraction_prev, methane_finundated, methane_soil_finundated, methane_soil_zwt, &
		     methane_surf_flux_wetland, methane_surf_flux_soil, methane_surf_flux_lake, &
		     methane_surf_flux_rice, methane_surf_aere_soil, methane_surf_aere_rice, &
		     methane_surf_ebul_soil, methane_surf_ebul_rice, methane_surf_diff_soil, &
		     methane_surf_diff_rice, methane_prod_tot_soil, methane_prod_tot_rice, &
		     methane_oxid_tot_soil, methane_oxid_tot_rice

		IMPLICIT NONE
		integer ,intent(in) :: istep

		integer ,intent(in) :: i         ! patch index
		integer ,intent(in) :: idate(1:3)  ! current date (year, day of the year, seconds of the day)
		integer, intent(in) :: &
							patchclass  ,&! land patch class of USGS classification or others
							patchtype     ! land patch type (0=soil, 1=urban and built-up,
								          ! 2=wetland, 3=land ice, 4=land water bodies, 99 = ocean)
		real(r8),intent(in) :: deltim    ! time step in seconds
		integer ,intent(in) :: &
				lb,&
				snl
		real(r8),intent(in) :: dlon     ! longitude (degrees)
		real(r8),intent(in) :: dlat     ! latitude (degrees)

		real(r8),intent(in) :: &
				z_soisno   (maxsnl+1:nl_soil) , &! layer depth (m)
				dz_soisno  (maxsnl+1:nl_soil) , &! layer thickness (m)
				zi_soisno  (maxsnl:nl_soil)   , &! interface level below a "z" level (m)
				t_soisno   (maxsnl+1:nl_soil) , &! soil + snow layer temperature [K]
				t_grnd                 		  , &! ground surface temperature [k]
				wliq_soisno(maxsnl+1:nl_soil) , &! liquid water (kg/m2)
				wice_soisno(maxsnl+1:nl_soil) , &! ice lens (kg/m2)
				forc_t                        , &! temperature at agcm reference height [kelvin]
					forc_pbot                     , &! atmosphere pressure at the bottom of the atmos. model level [pa]
					forc_po2m                     , &! partial pressure of O2 at observational height [pa]
					forc_pco2m                    , &! partial pressure of CO2 at observational height [Pa]
					forc_us                       , &! eastward wind speed [m/s], used as U10 proxy for lake gas exchange
					forc_vs                       , &! northward wind speed [m/s], used as U10 proxy for lake gas exchange
					ustar                         , &! current friction velocity [m/s]
					fq                            , &! Monin-Obukhov moisture profile integral [-]
					zwt                           , &! the depth to water table [m]
				rootfr     (1:nl_soil)        , &! fraction of roots in each soil layer
				snowdp                        , &! snow depth (m)
				wat                           , &! total water storage [mm] (reserved interface)
				rsur                          , &! surface runoff [mm H2O/s] (reserved interface)
				etr                           , &! transpiration rate [mm/s]
				lakedepth                     , &! static/capacity lake depth [m]
				dz_lake     (1:nl_lake)       , &! current lake layer thickness [m]
				t_lake      (1:nl_lake)       , &! current lake layer temperature [K]
				lake_icefrac(1:nl_lake)       , &! lake frozen mass fraction for lake CH4 exchange
				wdsrf                         , &! depth of surface water [mm]
		wetwat                        , &! water storage in wetland [mm]
				bsw         (1:nl_soil)       , &! clapp and hornbereger "b" parameter [-]
				smp         (1:nl_soil)       , &! soil matrix potential [mm]
				porsl       (1:nl_soil)       , &! fraction of soil that is voids [-]
				lai                           , &! leaf area index
				sai                           , &! stem area index
				rootr       (1:nl_soil)       , &! water exchange between soil and root. Positive: soil->root [?]

				fsatmax                       , &! maximum saturated area fraction [-]
				fsatdcf                       , &! decay factor in calculation of saturated area fraction [1/m]
		frcsat                        , &! fraction of saturation area
				f_h2osfc                        ! fraction of surface water [-], maintained by WATER_2014 / WATER_VSF

		! Optional R1 rice paddy inputs from CoLMDRIVER.  When absent the
		! routine behaves exactly as before (default-off).
		logical, intent(in), optional :: is_rice_paddy_in
		real(r8), intent(in), optional :: rice_pft_frac_in

		logical  :: is_rice_paddy
		logical  :: rice_parameter_active
		logical  :: is_floodplain_active
		real(r8) :: rice_pft_frac
		integer  :: rice_dsh

		real(r8):: &
				crootfr  (1:nl_soil)     , &! fraction of roots for carbon in each soil layer
				pH                       , &! soil water pH
				cellorg  (1:nl_soil)     , &! column 3D org (kg/m3 organic matter)
				t_h2osfc             	    ! surface water temperature

			real(r8) :: somhr_loc, lithr_loc, rr_loc, agnpp_loc, bgnpp_loc, annsum_npp_loc
			real(r8) :: hr_vr_loc(1:nl_soil), fphr_loc(1:nl_soil)
			real(r8) :: o_scalar_loc(1:nl_soil), pot_f_nit_vr_loc(1:nl_soil)
				real(r8) :: microbe_conc_o2(1:nl_soil), microbe_conc_ch4(1:nl_soil)
				real(r8) :: microbial_prod_potential_eff(1:nl_soil)
				real(r8) :: microbial_oxid_potential_eff(1:nl_soil)
			real(r8) :: lai_eff
			real(r8) :: rootfr_eff(1:nl_soil)
			real(r8) :: rootr_eff(1:nl_soil)
				real(r8) :: forc_t_eff, forc_pbot_eff, forc_po2m_eff, forc_pco2m_eff
				real(r8) :: forc_us_eff, forc_vs_eff
			real(r8) :: fprev
			real(r8) :: rice_weight
			real(r8) :: component_fraction(N_METHANE_COMP)
			real(r8) :: component_lai(N_METHANE_COMP)
			real(r8) :: component_crootfr(1:nl_soil,N_METHANE_COMP)
			real(r8) :: component_rr(N_METHANE_COMP)
			real(r8) :: component_agnpp(N_METHANE_COMP)
			real(r8) :: component_bgnpp(N_METHANE_COMP)
			real(r8) :: component_annsum_npp(N_METHANE_COMP)
			logical :: component_veg_ready(N_METHANE_COMP)
			logical  :: bgc_inputs_ready
			TYPE(methane_column_result_type) :: soil_column, rice_column
			logical, save :: warned_missing_bgc_inputs = .false.

		! Use raw t_grnd as proxy for surface-water temperature.
		! Do NOT clamp to tfrz: the ponddiff formula handles sub-freezing
		! T (parabolic in t_c) without issue, and the downstream ice-block
		! path (Physics.F90 line ~2274) needs t_h2osfc < tfrz to fire.
		! Clamping was the cause of the dead ice-block branch.
		t_h2osfc = t_grnd

		! Sanitize atmospheric forcing used by methane gas exchange.  The
		! subroutine arguments are intent(in), so keep local copies rather
		! than mutating caller state.  This prevents coast/domain-edge spval
		! or NaN forcing from poisoning c_atm and lake ebullition pressure.

		! R1 rice paddy: unpack optional inputs, default off.
		IF (present(is_rice_paddy_in)) THEN
		   is_rice_paddy = is_rice_paddy_in
		ELSE
		   is_rice_paddy = .false.
		ENDIF
		IF (present(rice_pft_frac_in)) THEN
		   rice_pft_frac = rice_pft_frac_in
		ELSE
		   rice_pft_frac = 0._r8
		ENDIF
		rice_pft_frac = min(max(rice_pft_frac, 0._r8), 1._r8)
		rice_parameter_active = is_rice_paddy .and. rice_pft_frac > 0._r8
		IF (rice_parameter_active .and. .not. is_paddy_rice_live(i)) THEN
			   rice_dsh = rice_days_since_harvest(i, idate(2), idate(1))
		   rice_parameter_active = rice_dsh >= 0 .and. &
		      real(rice_dsh, r8) < DEF_METHANE%rice_drain_window_days
		ENDIF

		forc_t_eff = forc_t
		if (ieee_is_nan(forc_t_eff) .or. abs(forc_t_eff) >= 0.5_r8*abs(spval) .or. &
		    forc_t_eff < 150._r8 .or. forc_t_eff > 350._r8) forc_t_eff = 288.15_r8

		forc_pbot_eff = forc_pbot
		if (ieee_is_nan(forc_pbot_eff) .or. abs(forc_pbot_eff) >= 0.5_r8*abs(spval) .or. &
		    forc_pbot_eff <= 0._r8) forc_pbot_eff = 101325._r8

		forc_po2m_eff = forc_po2m
		if (ieee_is_nan(forc_po2m_eff) .or. abs(forc_po2m_eff) >= 0.5_r8*abs(spval) .or. &
		    forc_po2m_eff <= 0._r8) forc_po2m_eff = 0.2095_r8 * forc_pbot_eff

			forc_pco2m_eff = forc_pco2m
			if (ieee_is_nan(forc_pco2m_eff) .or. abs(forc_pco2m_eff) >= 0.5_r8*abs(spval) .or. &
			    forc_pco2m_eff < 0._r8) forc_pco2m_eff = 415.e-6_r8 * forc_pbot_eff

			forc_us_eff = forc_us
			if (ieee_is_nan(forc_us_eff) .or. abs(forc_us_eff) >= 0.5_r8*abs(spval)) forc_us_eff = 0._r8

			forc_vs_eff = forc_vs
			if (ieee_is_nan(forc_vs_eff) .or. abs(forc_vs_eff) >= 0.5_r8*abs(spval)) forc_vs_eff = 0._r8

			CALL mt_ch4_bgc_patch_inputs (i, rootfr, crootfr, pH, cellorg, &
			     somhr_loc, lithr_loc, hr_vr_loc, rr_loc, agnpp_loc, bgnpp_loc, &
			     annsum_npp_loc, fphr_loc, o_scalar_loc, pot_f_nit_vr_loc, bgc_inputs_ready)
			component_fraction = 0._r8
			component_fraction(METHANE_COMP_SOIL) = 1._r8
			component_lai = 0._r8
			component_crootfr = 0._r8
			component_rr = 0._r8
			component_agnpp = 0._r8
			component_bgnpp = 0._r8
			component_annsum_npp = 0._r8
			component_veg_ready = .false.
			IF (patchtype == 0 .and. is_rice_paddy) THEN
				CALL mt_ch4_bgc_component_veg_inputs(i, rootfr, component_fraction, &
				   component_lai, component_crootfr, component_rr, component_agnpp, &
				   component_bgnpp, component_annsum_npp, component_veg_ready)
			ENDIF
			IF (patchtype == 0 .and. .not. bgc_inputs_ready) THEN
				IF (p_iam_worker == 0 .and. .not. warned_missing_bgc_inputs) THEN
					write(6,*) 'WARNING: Soil methane running with sanitized BGC defaults before BGC/PFT inputs are initialized.'
					warned_missing_bgc_inputs = .true.
				ENDIF
			ENDIF

			! ---- Wetland vegetation proxy (patchtype==2 only) -------------------
			! IGBP class 11 wetland patches carry no PFT, so NPP / rootfr arrive
			! as 0 and methane_aere (Wania 2010) cannot fire fully, dropping
			! the plant-mediated CH4 pathway that observations show carries
			! 50-90% of total wetland CH4 efflux (Bridgham 2013 GCB).
			!
			! LAI handling: CoLM mksrfdata aggregates Yuan+2011 LAI onto
			! patchtype==2 patches (fveg0_igbp(11)=1), so `lai` has real
			! data + DEF_LAI_MONTHLY seasonality.  get_wetland_veg_proxy
			! trusts it when valid (0.1 < lai < 20) and only falls back to
			! a climate-zone peak when missing.
			!
			! For patchtype /= 2: use the inputs as-is (vegetated patches
			! already have their own PFT-driven LAI/NPP/roots, which
			! correctly drives the aerenchyma path even under seasonal
			! inundation, e.g. floodplain pasture or seasonally flooded forest).
			lai_eff       = lai
			rootfr_eff(:) = rootfr(:)
			rootr_eff(:)  = rootr(:)
			! Clear any prior aerenchyma override (so non-wetland patches use
			! the DEF_METHANE%* defaults via the getter fallback path).
			IF (allocated(wetland_aere_active)) THEN
				IF (i >= 1 .and. i <= size(wetland_aere_active)) &
					wetland_aere_active(i) = .false.
			ENDIF
			IF (patchtype == 2) THEN
				CALL get_wetland_veg_proxy (dlat, cellorg(1), lai, i, &
				     lai_eff, annsum_npp_loc, agnpp_loc, bgnpp_loc, rootfr_eff)
				! BgcLink set crootfr from the original (zero) rootfr; replace it
				! with the proxy profile so methane_prod distributes root
				! respiration into the right layers.
				crootfr(1:nl_soil) = rootfr_eff(1:nl_soil)
				! Wetland patches have no PFT-level root water uptake or root
				! respiration diagnostics.  Use the same shallow proxy profile
				! for transpiration loss and a conservative belowground-respiration
				! proxy so root O2 demand is not silently zero.
				rootr_eff(1:nl_soil) = rootfr_eff(1:nl_soil)
				rr_loc = max(rr_loc, 0.5_r8 * bgnpp_loc)
			ENDIF
			! R4 rice paddy aerenchyma: rice patches (patchtype==0) get their own
			! Wania-style override (Zone 6) so methane_aere uses rice tiller
			! geometry instead of the natural-wetland default.  NPP / rootfr stay
			! whatever CN provides for the actual rice CFT — only the aerenchyma
			! geometry channel is overridden.
			!
			! Gate on croplive_p so winter stubble/fallow on rice patches drops
			! back to the wetland default (otherwise the LAI-scaled aere keeps
			! a residual 50% transport even after harvest).  rice_paddy_alive
			! flag is set by the same any_paddy_rice_live() helper used in
			! methane() for finundated override, keeping the two physics
			! decisions consistent.
			! All BGC state access goes through MOD_MT_Methane_BgcLink.
			! Rice aerenchyma is configured inside the pure-rice component call
			! with an intensive fraction of one.  The patch rice fraction is used
			! only once, when the two completed columns are aggregated.

			! Biome-specific f_methane lookup (Bridgham 2013 GCB).  Sets
			! biome_f_methane_patch(i) which methane_prod consumes when
			! DEF_METHANE%use_biome_f_methane is true; otherwise returns
			! the legacy DEF_METHANE%f_methane scalar (backwards compatible).
			is_floodplain_active = .false.
			IF (patchtype == 0 .and. DEF_METHANE%use_routing_for_soil .and. &
			    allocated(f_inund_flood_patch) .and. allocated(wetland_frac_per_patch)) THEN
				IF (i >= 1 .and. i <= size(f_inund_flood_patch) .and. &
				    i <= size(wetland_frac_per_patch)) &
					! Physics allocates grid flood to static wetland first.  Apply
					! floodplain biome parameters only when a positive soil residual
					! remains, and preserve the user threshold on the grid fraction.
					is_floodplain_active = &
						f_inund_flood_patch(i) > DEF_METHANE%hybrid_soil_threshold .and. &
						f_inund_flood_patch(i) > wetland_frac_per_patch(i)
			ENDIF
			IF (allocated(biome_f_methane_patch) .and. &
			    i >= 1 .and. i <= size(biome_f_methane_patch)) THEN
				biome_f_methane_patch(i) = get_biome_f_methane (patchtype, dlat, cellorg(1), &
				                                                is_rice_paddy, rice_pft_frac, &
				                                                rice_parameter_active, is_floodplain_active)
			ENDIF

			! Biome-specific redoxlag lookup (Pangala 2017, Whalen 1990).
			! Tropical warm wetlands respond faster (~5-15d); cold boreal
			! peat slower (~30-60d).  Consumed by methane when
			! DEF_METHANE%use_biome_redoxlag is true.
			IF (patchtype /= 0 .and. allocated(biome_redoxlag_patch) .and. &
			    i >= 1 .and. i <= size(biome_redoxlag_patch)) THEN
				biome_redoxlag_patch(i) = get_biome_redoxlag (patchtype, dlat, cellorg(1), &
				                                              is_rice_paddy, rice_pft_frac, &
				                                              rice_parameter_active, is_floodplain_active)
			ENDIF

			! Soil patches use two independent methane-process columns.  Each
			! column completes its nonlinear CH4/O2/microbial solve before the
			! patch rice fraction is applied exactly once to the outputs.
			IF (patchtype == 0) THEN
				rice_weight = merge(rice_pft_frac, 0._r8, is_rice_paddy)
				rice_weight = min(max(rice_weight, 0._r8), 1._r8)
				IF (rice_weight <= 1.e-14_r8) rice_weight = 0._r8
				IF (rice_weight >= 1._r8-1.e-14_r8) rice_weight = 1._r8
				IF (any(component_veg_ready) .and. &
				    abs(component_fraction(METHANE_COMP_RICE)-rice_weight) > 1.e-8_r8) THEN
					CALL CoLM_stop(' ***** ERROR: rice methane component fraction disagrees with PFT bridge')
				ENDIF
				CALL repartition_methane_column_state(i, rice_fraction_prev(i), rice_weight)
				CALL repartition_methane_microbes(i, rice_fraction_prev(i), rice_weight)
				IF (rice_weight <= 1.e-14_r8) THEN
					CALL run_methane_component(METHANE_COMP_SOIL, .false., soil_column)
					rice_column = soil_column
				ELSEIF (rice_weight >= 1._r8-1.e-14_r8) THEN
					CALL run_methane_component(METHANE_COMP_RICE, .true., rice_column)
					soil_column = rice_column
				ELSE
					CALL run_methane_component(METHANE_COMP_SOIL, .false., soil_column)
					CALL run_methane_component(METHANE_COMP_RICE, .true., rice_column)
				ENDIF
				CALL aggregate_methane_microbes(i, rice_weight)
				CALL aggregate_methane_columns(soil_column, rice_column, rice_weight)
				CALL mt_ch4_bgc_finalize_step(i, patchtype, deltim, net_methane(i), co2_hr_methane(i))
				rice_fraction_prev(i) = rice_weight
				RETURN
			ENDIF
					! f_h2osfc is maintained by the water module before methane_driver.
					! Lake CH4 uses the CTSM-style sediment-carbon pathway, not the optional
					! soil microbial-pool override; keep lake microbial pools from evolving as
					! diagnostics-only state.
					IF (patchtype /= 4) THEN
						IF (abs(fsat_bef(i)) < 1.e30_r8 .and. fsat_bef(i) >= 0._r8 .and. &
						    fsat_bef(i) <= 1._r8) THEN
							fprev = fsat_bef(i)
							microbe_conc_o2(:) = conc_o2_sat(1:nl_soil,i) * fprev + &
								conc_o2_unsat(1:nl_soil,i) * (1._r8 - fprev)
							microbe_conc_ch4(:) = conc_methane_sat(1:nl_soil,i) * fprev + &
								conc_methane_unsat(1:nl_soil,i) * (1._r8 - fprev)
						ELSE
							microbe_conc_o2(:) = conc_o2(1:nl_soil,i)
							microbe_conc_ch4(:) = conc_methane(1:nl_soil,i)
						ENDIF
							CALL methane_microbes_step (i, METHANE_COMP_SOIL, deltim, t_soisno(1:nl_soil), &
							     microbe_conc_o2, microbe_conc_ch4, hr_vr_loc, cellorg)
							CALL aggregate_methane_microbes(i, 0._r8)
						ELSE
							! Lake CH4 bypasses the optional soil microbial-pool override.
							! Reset every inactive per-step diagnostic through its owner;
							! prognostic biomass remains carried by the lake sediment.
							CALL reset_methane_inactive_lake_microbe_diagnostics(i)
						ENDIF

						microbial_prod_potential_eff(:) = 0._r8
						microbial_oxid_potential_eff(:) = 0._r8
						IF (allocated(microbial_prod_potential) .and. allocated(microbial_oxid_potential)) THEN
							IF (i >= lbound(microbial_prod_potential,2) .and. i <= ubound(microbial_prod_potential,2) .and. &
							    i >= lbound(microbial_oxid_potential,2) .and. i <= ubound(microbial_oxid_potential,2)) THEN
								microbial_prod_potential_eff(:) = microbial_prod_potential(1:nl_soil,i)
								microbial_oxid_potential_eff(:) = microbial_oxid_potential(1:nl_soil,i)
							ENDIF
						ENDIF

			CALL methane (istep,i,idate(1:3),patchclass,patchtype,lb,snl,dlon,dlat,deltim,&
		z_soisno(maxsnl+1:),dz_soisno(maxsnl+1:),zi_soisno(maxsnl:),t_soisno(maxsnl+1:),&
		t_grnd,wliq_soisno(maxsnl+1:),wice_soisno(maxsnl+1:),&
			forc_t_eff,forc_pbot_eff,forc_po2m_eff,forc_pco2m_eff,forc_us_eff,forc_vs_eff,&
		zwt,rootfr_eff,snowdp,wat,rsur,etr,lakedepth,dz_lake,t_lake,lake_icefrac,wdsrf,wetwat,bsw,&
		smp,porsl,lai_eff,sai,rootr_eff,&
		annsum_npp_loc, rr_loc,&
		fsatmax,fsatdcf,frcsat,&
		agnpp_loc, bgnpp_loc, somhr_loc,&
		crootfr(1:nl_soil), lithr_loc, hr_vr_loc(1:nl_soil), o_scalar_loc(1:nl_soil), &
		fphr_loc(1:nl_soil), pot_f_nit_vr_loc(1:nl_soil), pH,&
			cellorg(1:nl_soil),t_h2osfc,organic_max,&
				microbial_prod_potential_eff, microbial_oxid_potential_eff, &
		!!!! --------------------------------------------------------------------------------------------------------
		!!!!                                         sum data
		!!!! --------------------------------------------------------------------------------------------------------
		net_methane(i), &
		methane_prod_depth(1:nl_soil,i), o2_decomp_depth(1:nl_soil,i), co2_decomp_depth(1:nl_soil,i), &
		methane_oxid_depth(1:nl_soil,i), o2_oxid_depth(1:nl_soil,i), co2_oxid_depth(1:nl_soil,i), &
		methane_aere_depth(1:nl_soil,i), methane_tran_depth(1:nl_soil,i), &
		o2_aere_depth(1:nl_soil,i), co2_aere_depth(1:nl_soil,i), &
		methane_ebul_depth(1:nl_soil,i), &
		o2stress(1:nl_soil,i), methane_stress(1:nl_soil,i), &
		methane_surf_flux_tot(i), methane_surf_flux_tot_phys(i), methane_surf_aere(i), methane_surf_ebul(i), methane_surf_diff(i), &
		methane_surf_diff_phys(i), &
		methane_balance_residual(i), methane_ch4_clip_credit(i), o2_cap_loss(i), o2_cap_gain(i), &
		methane_ebul_tot(i), methane_prod_tot(i), methane_oxid_tot(i), &
		co2_decomp_tot(i), co2_oxid_tot(i), co2_aere_tot(i), co2_net_tot(i), &
		totcol_methane(i), grnd_methane_cond(i), conc_o2(1:nl_soil,i), conc_methane(1:nl_soil,i), &
		!!!! --------------------------------------------------------------------------------------------------------
		!!!! --------------------------------------------------------------------------------------------------------
		!!!!                                         sum data (unsaturated / saturated)
		!!!! --------------------------------------------------------------------------------------------------------
		net_methane_unsat(i), net_methane_sat(i), &
		methane_prod_depth_unsat(1:nl_soil,i), methane_prod_depth_sat(1:nl_soil,i), &
		o2_decomp_depth_unsat(1:nl_soil,i), o2_decomp_depth_sat(1:nl_soil,i), &
		co2_decomp_depth_unsat(1:nl_soil,i), co2_decomp_depth_sat(1:nl_soil,i), &
		methane_oxid_depth_unsat(1:nl_soil,i), methane_oxid_depth_sat(1:nl_soil,i), &
		o2_oxid_depth_unsat(1:nl_soil,i), o2_oxid_depth_sat(1:nl_soil,i), &
		co2_oxid_depth_unsat(1:nl_soil,i), co2_oxid_depth_sat(1:nl_soil,i), &
		methane_aere_depth_unsat(1:nl_soil,i), methane_aere_depth_sat(1:nl_soil,i), &
		methane_tran_depth_unsat(1:nl_soil,i), methane_tran_depth_sat(1:nl_soil,i), &
		o2_aere_depth_unsat(1:nl_soil,i), o2_aere_depth_sat(1:nl_soil,i), &
		co2_aere_depth_unsat(1:nl_soil,i), co2_aere_depth_sat(1:nl_soil,i), &
		methane_ebul_depth_unsat(1:nl_soil,i), methane_ebul_depth_sat(1:nl_soil,i), &
		o2stress_unsat(1:nl_soil,i), o2stress_sat(1:nl_soil,i), &
		methane_stress_unsat(1:nl_soil,i), methane_stress_sat(1:nl_soil,i), &
		methane_surf_flux_tot_unsat(i), methane_surf_flux_tot_sat(i), &
		methane_surf_aere_unsat(i), methane_surf_aere_sat(i), &
		methane_surf_ebul_unsat(i), methane_surf_ebul_sat(i), &
		methane_surf_diff_unsat(i), methane_surf_diff_sat(i), &
		methane_surf_diff_phys_unsat(i), methane_surf_diff_phys_sat(i), &
		methane_ebul_tot_unsat(i), methane_ebul_tot_sat(i), &
		methane_prod_tot_unsat(i), methane_prod_tot_sat(i), &
		methane_oxid_tot_unsat(i), methane_oxid_tot_sat(i), &
		co2_decomp_tot_unsat(i), co2_decomp_tot_sat(i), &
		co2_oxid_tot_unsat(i), co2_oxid_tot_sat(i), &
		co2_net_tot_unsat(i), co2_net_tot_sat(i), &
		totcol_methane_unsat(i), totcol_methane_sat(i), &
		grnd_methane_cond_unsat(i), grnd_methane_cond_sat(i), &
		conc_o2_unsat(1:nl_soil,i), conc_o2_sat(1:nl_soil,i), &
		conc_methane_unsat(1:nl_soil,i), conc_methane_sat(1:nl_soil,i), &
		methane_prod_depth_lake(1:nl_soil,i), methane_oxid_depth_lake(1:nl_soil,i), &
		methane_ebul_depth_lake(1:nl_soil,i), &
		co2_decomp_depth_lake(1:nl_soil,i), co2_oxid_depth_lake(1:nl_soil,i), &
		methane_surf_ebul_lake(i), methane_surf_diff_lake(i), methane_surf_flux_tot_lake(i), &
		methane_prod_tot_lake(i), methane_oxid_tot_lake(i), methane_ebul_tot_lake(i), &
		co2_decomp_tot_lake(i), co2_oxid_tot_lake(i), co2_net_tot_lake(i), &
		totcol_methane_lake(i), grnd_methane_cond_lake(i), conc_o2_lake(1:nl_soil,i), conc_methane_lake(1:nl_soil,i), &
		lake_water_ch4_stock(i), lake_water_o2_stock(i), lake_frozen_ch4_stock(i), lake_frozen_o2_stock(i), &
		lake_liquid_fraction_prev(i), lake_water_ch4_oxid(i), &
		lake_sed_ch4_flux(i), lake_sed_o2_flux(i), lake_air_o2_flux(i), &
		!!!! --------------------------------------------------------------------------------------------------------
		c_atm(1:3,i), forc_pmethanem(i), layer_sat_lag(1:nl_soil,i), lake_soilc(1:nl_soil,i), &
		annavg_agnpp(i), annavg_bgnpp(i), annavg_somhr(i), annavg_finrw(i), &
		tempavg_agnpp(i), tempavg_bgnpp(i), annsum_counter(i), tempavg_somhr(i), &
		tempavg_finrw(i), fsat_bef(i), finundated_lag(i), methane_dfsat_tot(i), f_h2osfc, &
			is_rice_paddy_in=is_rice_paddy, rice_pft_frac_in=rice_pft_frac, &
			ustar_in=ustar, fq_in=fq)
			CALL mt_ch4_bgc_finalize_step(i, patchtype, deltim, net_methane(i), co2_hr_methane(i))

	CONTAINS

		SUBROUTINE run_methane_component(component, rice_column_active, result)
			integer, intent(in) :: component
			logical, intent(in) :: rice_column_active
			TYPE(methane_column_result_type), intent(out) :: result
			real(r8) :: column_fraction, column_fsat
			real(r8) :: column_totcol, column_totcol_unsat, column_totcol_sat
			real(r8) :: column_grnd_cond, column_grnd_cond_unsat, column_grnd_cond_sat
			real(r8) :: column_dfsat_tot
			real(r8) :: column_conc_o2(nl_soil), column_conc_ch4(nl_soil)
			real(r8) :: finundated_used, finundated_default_used
			real(r8) :: column_lai, column_crootfr(nl_soil)
			real(r8) :: column_rr, column_agnpp, column_bgnpp, column_annsum_npp

			column_fraction = merge(1._r8, 0._r8, rice_column_active)
			column_lai = lai_eff
			column_crootfr = crootfr
			column_rr = rr_loc
			column_agnpp = agnpp_loc
			column_bgnpp = bgnpp_loc
			column_annsum_npp = annsum_npp_loc
			IF (component_veg_ready(component)) THEN
				column_lai = component_lai(component)
				column_crootfr = component_crootfr(:,component)
				column_rr = component_rr(component)
				column_agnpp = component_agnpp(component)
				column_bgnpp = component_bgnpp(component)
				column_annsum_npp = component_annsum_npp(component)
			ENDIF
			IF (allocated(wetland_aere_active)) wetland_aere_active(i) = .false.
			IF (rice_column_active .and. is_paddy_rice_live(i)) THEN
				CALL get_rice_veg_proxy(column_lai, i, 1._r8)
			ENDIF

			is_floodplain_active = .false.
			IF (.not. rice_column_active) THEN
				is_floodplain_active = DEF_METHANE%use_routing_for_soil .and. &
				   allocated(f_inund_flood_patch) .and. allocated(wetland_frac_per_patch) .and. &
				   i >= 1 .and. i <= size(f_inund_flood_patch) .and. &
				   i <= size(wetland_frac_per_patch) .and. &
				   f_inund_flood_patch(i) > DEF_METHANE%hybrid_soil_threshold .and. &
				   f_inund_flood_patch(i) > wetland_frac_per_patch(i)
			ENDIF
			IF (allocated(biome_f_methane_patch)) THEN
				biome_f_methane_patch(i) = get_biome_f_methane(patchtype, dlat, cellorg(1), &
				   rice_column_active, column_fraction, &
				   rice_column_active .and. rice_parameter_active, is_floodplain_active)
			ENDIF
			IF (allocated(biome_redoxlag_patch)) THEN
				biome_redoxlag_patch(i) = get_biome_redoxlag(patchtype, dlat, cellorg(1), &
				   rice_column_active, column_fraction, rice_column_active .and. rice_parameter_active, &
				   is_floodplain_active)
			ENDIF

			column_fsat = fsat_bef_component(component,i)
			IF (abs(column_fsat) < 1.e30_r8 .and. column_fsat >= 0._r8 .and. column_fsat <= 1._r8) THEN
				microbe_conc_o2 = conc_o2_sat_component(:,component,i) * column_fsat + &
				   conc_o2_unsat_component(:,component,i) * (1._r8 - column_fsat)
				microbe_conc_ch4 = conc_methane_sat_component(:,component,i) * column_fsat + &
				   conc_methane_unsat_component(:,component,i) * (1._r8 - column_fsat)
			ELSE
				microbe_conc_o2 = 0.5_r8 * (conc_o2_sat_component(:,component,i) + &
				   conc_o2_unsat_component(:,component,i))
				microbe_conc_ch4 = 0.5_r8 * (conc_methane_sat_component(:,component,i) + &
				   conc_methane_unsat_component(:,component,i))
			ENDIF
			CALL methane_microbes_step(i, component, deltim, t_soisno(1:nl_soil), &
			   microbe_conc_o2, microbe_conc_ch4, hr_vr_loc, cellorg)
			microbial_prod_potential_eff = 0._r8
			microbial_oxid_potential_eff = 0._r8
			IF (allocated(microbial_prod_potential_comp)) THEN
				microbial_prod_potential_eff = microbial_prod_potential_comp(:,component,i)
				microbial_oxid_potential_eff = microbial_oxid_potential_comp(:,component,i)
			ENDIF

			column_totcol_unsat = sum(conc_methane_unsat_component(:,component,i) * dz_soisno(1:nl_soil))
			column_totcol_sat = sum(conc_methane_sat_component(:,component,i) * dz_soisno(1:nl_soil))
			IF (column_fsat >= 0._r8 .and. column_fsat <= 1._r8) THEN
				column_totcol = column_fsat * column_totcol_sat + (1._r8-column_fsat) * column_totcol_unsat
				column_conc_o2 = column_fsat * conc_o2_sat_component(:,component,i) + &
				   (1._r8-column_fsat) * conc_o2_unsat_component(:,component,i)
				column_conc_ch4 = column_fsat * conc_methane_sat_component(:,component,i) + &
				   (1._r8-column_fsat) * conc_methane_unsat_component(:,component,i)
			ELSE
				column_totcol = 0.5_r8 * (column_totcol_sat + column_totcol_unsat)
				column_conc_o2 = microbe_conc_o2
				column_conc_ch4 = microbe_conc_ch4
			ENDIF
			column_grnd_cond = grnd_methane_cond(i)
			column_grnd_cond_unsat = grnd_methane_cond_unsat(i)
			column_grnd_cond_sat = grnd_methane_cond_sat(i)

			CALL methane (istep,i,idate(1:3),patchclass,patchtype,lb,snl,dlon,dlat,deltim,&
			z_soisno(maxsnl+1:),dz_soisno(maxsnl+1:),zi_soisno(maxsnl:),t_soisno(maxsnl+1:),&
			t_grnd,wliq_soisno(maxsnl+1:),wice_soisno(maxsnl+1:),&
			forc_t_eff,forc_pbot_eff,forc_po2m_eff,forc_pco2m_eff,forc_us_eff,forc_vs_eff,&
			zwt,rootfr_eff,snowdp,wat,rsur,etr,lakedepth,dz_lake,t_lake,lake_icefrac,wdsrf,wetwat,bsw,&
			smp,porsl,column_lai,sai,rootr_eff,column_annsum_npp,column_rr,fsatmax,fsatdcf,frcsat,&
			column_agnpp,column_bgnpp,somhr_loc,column_crootfr,lithr_loc,hr_vr_loc,o_scalar_loc,&
			fphr_loc,pot_f_nit_vr_loc,pH,cellorg,t_h2osfc,organic_max,&
			microbial_prod_potential_eff,microbial_oxid_potential_eff,&
			net_methane(i),methane_prod_depth(:,i),o2_decomp_depth(:,i),co2_decomp_depth(:,i),&
			methane_oxid_depth(:,i),o2_oxid_depth(:,i),co2_oxid_depth(:,i),&
			methane_aere_depth(:,i),methane_tran_depth(:,i),o2_aere_depth(:,i),co2_aere_depth(:,i),&
			methane_ebul_depth(:,i),o2stress(:,i),methane_stress(:,i),&
			methane_surf_flux_tot(i),methane_surf_flux_tot_phys(i),methane_surf_aere(i),&
			methane_surf_ebul(i),methane_surf_diff(i),methane_surf_diff_phys(i),&
			methane_balance_residual(i),methane_ch4_clip_credit(i),o2_cap_loss(i),o2_cap_gain(i),&
			methane_ebul_tot(i),methane_prod_tot(i),methane_oxid_tot(i),&
			co2_decomp_tot(i),co2_oxid_tot(i),co2_aere_tot(i),co2_net_tot(i),&
			column_totcol,column_grnd_cond,column_conc_o2,column_conc_ch4,&
			net_methane_unsat(i),net_methane_sat(i),&
			methane_prod_depth_unsat(:,i),methane_prod_depth_sat(:,i),&
			o2_decomp_depth_unsat(:,i),o2_decomp_depth_sat(:,i),&
			co2_decomp_depth_unsat(:,i),co2_decomp_depth_sat(:,i),&
			methane_oxid_depth_unsat(:,i),methane_oxid_depth_sat(:,i),&
			o2_oxid_depth_unsat(:,i),o2_oxid_depth_sat(:,i),&
			co2_oxid_depth_unsat(:,i),co2_oxid_depth_sat(:,i),&
			methane_aere_depth_unsat(:,i),methane_aere_depth_sat(:,i),&
			methane_tran_depth_unsat(:,i),methane_tran_depth_sat(:,i),&
			o2_aere_depth_unsat(:,i),o2_aere_depth_sat(:,i),&
			co2_aere_depth_unsat(:,i),co2_aere_depth_sat(:,i),&
			methane_ebul_depth_unsat(:,i),methane_ebul_depth_sat(:,i),&
			o2stress_unsat(:,i),o2stress_sat(:,i),methane_stress_unsat(:,i),methane_stress_sat(:,i),&
			methane_surf_flux_tot_unsat(i),methane_surf_flux_tot_sat(i),&
			methane_surf_aere_unsat(i),methane_surf_aere_sat(i),&
			methane_surf_ebul_unsat(i),methane_surf_ebul_sat(i),&
			methane_surf_diff_unsat(i),methane_surf_diff_sat(i),&
			methane_surf_diff_phys_unsat(i),methane_surf_diff_phys_sat(i),&
			methane_ebul_tot_unsat(i),methane_ebul_tot_sat(i),&
			methane_prod_tot_unsat(i),methane_prod_tot_sat(i),&
			methane_oxid_tot_unsat(i),methane_oxid_tot_sat(i),&
			co2_decomp_tot_unsat(i),co2_decomp_tot_sat(i),&
			co2_oxid_tot_unsat(i),co2_oxid_tot_sat(i),co2_net_tot_unsat(i),co2_net_tot_sat(i),&
			column_totcol_unsat,column_totcol_sat,column_grnd_cond_unsat,column_grnd_cond_sat,&
			conc_o2_unsat_component(:,component,i),conc_o2_sat_component(:,component,i),&
			conc_methane_unsat_component(:,component,i),conc_methane_sat_component(:,component,i),&
			methane_prod_depth_lake(:,i),methane_oxid_depth_lake(:,i),methane_ebul_depth_lake(:,i),&
			co2_decomp_depth_lake(:,i),co2_oxid_depth_lake(:,i),methane_surf_ebul_lake(i),&
			methane_surf_diff_lake(i),methane_surf_flux_tot_lake(i),methane_prod_tot_lake(i),&
			methane_oxid_tot_lake(i),methane_ebul_tot_lake(i),co2_decomp_tot_lake(i),&
			co2_oxid_tot_lake(i),co2_net_tot_lake(i),totcol_methane_lake(i),&
			grnd_methane_cond_lake(i),conc_o2_lake(:,i),conc_methane_lake(:,i),&
			lake_water_ch4_stock(i),lake_water_o2_stock(i),lake_frozen_ch4_stock(i),lake_frozen_o2_stock(i),&
			lake_liquid_fraction_prev(i),lake_water_ch4_oxid(i),&
			lake_sed_ch4_flux(i),lake_sed_o2_flux(i),lake_air_o2_flux(i),&
			c_atm(:,i),forc_pmethanem(i),layer_sat_lag_component(:,component,i),lake_soilc(:,i),&
			annavg_agnpp_component(component,i),annavg_bgnpp_component(component,i),&
			annavg_somhr_component(component,i),annavg_finrw_component(component,i),&
			tempavg_agnpp_component(component,i),tempavg_bgnpp_component(component,i),&
			annsum_counter_component(component,i),tempavg_somhr_component(component,i),&
			tempavg_finrw_component(component,i),fsat_bef_component(component,i),&
			finundated_lag_component(component,i),column_dfsat_tot,f_h2osfc,&
			is_rice_paddy_in=rice_column_active,rice_pft_frac_in=column_fraction,&
			ustar_in=ustar,fq_in=fq,&
			store_patch_diagnostics_in=.false.,finundated_used_out=finundated_used,&
			finundated_default_out=finundated_default_used)

			CALL capture_methane_column(component, result, column_totcol, column_grnd_cond, &
			   column_totcol_unsat, column_totcol_sat, column_grnd_cond_unsat, &
			   column_grnd_cond_sat, column_dfsat_tot, column_conc_o2, column_conc_ch4, &
			   finundated_used, finundated_default_used)
		END SUBROUTINE run_methane_component

		SUBROUTINE capture_methane_column(component, result, column_totcol, column_grnd_cond, &
		   column_totcol_unsat, column_totcol_sat, column_grnd_cond_unsat, &
		   column_grnd_cond_sat, column_dfsat_tot, column_conc_o2, column_conc_ch4, &
		   finundated_used, finundated_default_used)
			integer, intent(in) :: component
			TYPE(methane_column_result_type), intent(out) :: result
			real(r8), intent(in) :: column_totcol, column_grnd_cond
			real(r8), intent(in) :: column_totcol_unsat, column_totcol_sat
			real(r8), intent(in) :: column_grnd_cond_unsat, column_grnd_cond_sat
			real(r8), intent(in) :: column_dfsat_tot
			real(r8), intent(in) :: column_conc_o2(nl_soil), column_conc_ch4(nl_soil)
			real(r8), intent(in) :: finundated_used, finundated_default_used

			result%net = net_methane(i)
			result%surf_flux = methane_surf_flux_tot(i)
			result%surf_flux_phys = methane_surf_flux_tot_phys(i)
			result%surf_aere = methane_surf_aere(i)
			result%surf_ebul = methane_surf_ebul(i)
			result%surf_diff = methane_surf_diff(i)
			result%surf_diff_phys = methane_surf_diff_phys(i)
			result%balance_residual = methane_balance_residual(i)
			result%ch4_clip_credit = methane_ch4_clip_credit(i)
			result%o2_cap_loss = o2_cap_loss(i)
			result%o2_cap_gain = o2_cap_gain(i)
			result%ebul_tot = methane_ebul_tot(i)
			result%prod_tot = methane_prod_tot(i)
			result%oxid_tot = methane_oxid_tot(i)
			result%co2_decomp_tot = co2_decomp_tot(i)
			result%co2_oxid_tot = co2_oxid_tot(i)
			result%co2_aere_tot = co2_aere_tot(i)
			result%co2_net_tot = co2_net_tot(i)
			result%totcol = column_totcol
			result%grnd_cond = column_grnd_cond
			result%dfsat_tot = column_dfsat_tot
			result%finundated = finundated_used
			result%finundated_default = finundated_default_used
			result%net_unsat = net_methane_unsat(i)
			result%net_sat = net_methane_sat(i)
			result%surf_flux_unsat = methane_surf_flux_tot_unsat(i)
			result%surf_flux_sat = methane_surf_flux_tot_sat(i)
			result%surf_aere_unsat = methane_surf_aere_unsat(i)
			result%surf_aere_sat = methane_surf_aere_sat(i)
			result%surf_ebul_unsat = methane_surf_ebul_unsat(i)
			result%surf_ebul_sat = methane_surf_ebul_sat(i)
			result%surf_diff_unsat = methane_surf_diff_unsat(i)
			result%surf_diff_sat = methane_surf_diff_sat(i)
			result%surf_diff_phys_unsat = methane_surf_diff_phys_unsat(i)
			result%surf_diff_phys_sat = methane_surf_diff_phys_sat(i)
			result%ebul_tot_unsat = methane_ebul_tot_unsat(i)
			result%ebul_tot_sat = methane_ebul_tot_sat(i)
			result%prod_tot_unsat = methane_prod_tot_unsat(i)
			result%prod_tot_sat = methane_prod_tot_sat(i)
			result%oxid_tot_unsat = methane_oxid_tot_unsat(i)
			result%oxid_tot_sat = methane_oxid_tot_sat(i)
			result%co2_decomp_tot_unsat = co2_decomp_tot_unsat(i)
			result%co2_decomp_tot_sat = co2_decomp_tot_sat(i)
			result%co2_oxid_tot_unsat = co2_oxid_tot_unsat(i)
			result%co2_oxid_tot_sat = co2_oxid_tot_sat(i)
			result%co2_net_tot_unsat = co2_net_tot_unsat(i)
			result%co2_net_tot_sat = co2_net_tot_sat(i)
			result%totcol_unsat = column_totcol_unsat
			result%totcol_sat = column_totcol_sat
			result%grnd_cond_unsat = column_grnd_cond_unsat
			result%grnd_cond_sat = column_grnd_cond_sat
			result%prod = methane_prod_depth(:,i)
			result%o2_decomp = o2_decomp_depth(:,i)
			result%co2_decomp = co2_decomp_depth(:,i)
			result%oxid = methane_oxid_depth(:,i)
			result%o2_oxid = o2_oxid_depth(:,i)
			result%co2_oxid = co2_oxid_depth(:,i)
			result%aere = methane_aere_depth(:,i)
			result%tran = methane_tran_depth(:,i)
			result%o2_aere = o2_aere_depth(:,i)
			result%co2_aere = co2_aere_depth(:,i)
			result%ebul = methane_ebul_depth(:,i)
			result%o2stress = o2stress(:,i)
			result%ch4stress = methane_stress(:,i)
			result%prod_unsat = methane_prod_depth_unsat(:,i)
			result%prod_sat = methane_prod_depth_sat(:,i)
			result%o2_decomp_unsat = o2_decomp_depth_unsat(:,i)
			result%o2_decomp_sat = o2_decomp_depth_sat(:,i)
			result%co2_decomp_unsat = co2_decomp_depth_unsat(:,i)
			result%co2_decomp_sat = co2_decomp_depth_sat(:,i)
			result%oxid_unsat = methane_oxid_depth_unsat(:,i)
			result%oxid_sat = methane_oxid_depth_sat(:,i)
			result%o2_oxid_unsat = o2_oxid_depth_unsat(:,i)
			result%o2_oxid_sat = o2_oxid_depth_sat(:,i)
			result%co2_oxid_unsat = co2_oxid_depth_unsat(:,i)
			result%co2_oxid_sat = co2_oxid_depth_sat(:,i)
			result%aere_unsat = methane_aere_depth_unsat(:,i)
			result%aere_sat = methane_aere_depth_sat(:,i)
			result%tran_unsat = methane_tran_depth_unsat(:,i)
			result%tran_sat = methane_tran_depth_sat(:,i)
			result%o2_aere_unsat = o2_aere_depth_unsat(:,i)
			result%o2_aere_sat = o2_aere_depth_sat(:,i)
			result%co2_aere_unsat = co2_aere_depth_unsat(:,i)
			result%co2_aere_sat = co2_aere_depth_sat(:,i)
			result%ebul_unsat = methane_ebul_depth_unsat(:,i)
			result%ebul_sat = methane_ebul_depth_sat(:,i)
			result%o2stress_unsat = o2stress_unsat(:,i)
			result%o2stress_sat = o2stress_sat(:,i)
			result%ch4stress_unsat = methane_stress_unsat(:,i)
			result%ch4stress_sat = methane_stress_sat(:,i)
			result%conc_o2 = column_conc_o2
			result%conc_ch4 = column_conc_ch4
			result%conc_o2_unsat = conc_o2_unsat_component(:,component,i)
			result%conc_o2_sat = conc_o2_sat_component(:,component,i)
			result%conc_ch4_unsat = conc_methane_unsat_component(:,component,i)
			result%conc_ch4_sat = conc_methane_sat_component(:,component,i)
		END SUBROUTINE capture_methane_column

		SUBROUTINE aggregate_methane_columns(soil, rice, rice_fraction)
			TYPE(methane_column_result_type), intent(in) :: soil, rice
			real(r8), intent(in) :: rice_fraction
			real(r8) :: ws, wr, hs, hr, sat_area, unsat_area
			real(r8) :: wss, wsr, wus, wur

			wr = min(max(rice_fraction, 0._r8), 1._r8)
			ws = 1._r8 - wr
			hs = min(max(soil%finundated, 0._r8), 1._r8)
			hr = min(max(rice%finundated, 0._r8), 1._r8)
			sat_area = ws*hs + wr*hr
			unsat_area = ws*(1._r8-hs) + wr*(1._r8-hr)
			IF (sat_area > tiny(1._r8)) THEN
				wss = ws*hs/sat_area
				wsr = wr*hr/sat_area
			ELSE
				wss = ws
				wsr = wr
			ENDIF
			IF (unsat_area > tiny(1._r8)) THEN
				wus = ws*(1._r8-hs)/unsat_area
				wur = wr*(1._r8-hr)/unsat_area
			ELSE
				wus = ws
				wur = wr
			ENDIF
			net_methane(i) = ws*soil%net + wr*rice%net
			methane_surf_flux_tot(i) = ws*soil%surf_flux + wr*rice%surf_flux
			methane_surf_flux_tot_phys(i) = ws*soil%surf_flux_phys + wr*rice%surf_flux_phys
			methane_surf_aere(i) = ws*soil%surf_aere + wr*rice%surf_aere
			methane_surf_ebul(i) = ws*soil%surf_ebul + wr*rice%surf_ebul
			methane_surf_diff(i) = ws*soil%surf_diff + wr*rice%surf_diff
			methane_surf_diff_phys(i) = ws*soil%surf_diff_phys + wr*rice%surf_diff_phys
			methane_balance_residual(i) = ws*soil%balance_residual + wr*rice%balance_residual
			methane_ch4_clip_credit(i) = ws*soil%ch4_clip_credit + wr*rice%ch4_clip_credit
			o2_cap_loss(i) = ws*soil%o2_cap_loss + wr*rice%o2_cap_loss
			o2_cap_gain(i) = ws*soil%o2_cap_gain + wr*rice%o2_cap_gain
			methane_ebul_tot(i) = ws*soil%ebul_tot + wr*rice%ebul_tot
			methane_prod_tot(i) = ws*soil%prod_tot + wr*rice%prod_tot
			methane_oxid_tot(i) = ws*soil%oxid_tot + wr*rice%oxid_tot
			co2_decomp_tot(i) = ws*soil%co2_decomp_tot + wr*rice%co2_decomp_tot
			co2_oxid_tot(i) = ws*soil%co2_oxid_tot + wr*rice%co2_oxid_tot
			co2_aere_tot(i) = ws*soil%co2_aere_tot + wr*rice%co2_aere_tot
			co2_net_tot(i) = ws*soil%co2_net_tot + wr*rice%co2_net_tot
			totcol_methane(i) = ws*soil%totcol + wr*rice%totcol
			grnd_methane_cond(i) = ws*soil%grnd_cond + wr*rice%grnd_cond
			methane_dfsat_tot(i) = ws*soil%dfsat_tot + wr*rice%dfsat_tot

			net_methane_unsat(i) = wus*soil%net_unsat + wur*rice%net_unsat
			net_methane_sat(i) = wss*soil%net_sat + wsr*rice%net_sat
			methane_surf_flux_tot_unsat(i) = wus*soil%surf_flux_unsat + wur*rice%surf_flux_unsat
			methane_surf_flux_tot_sat(i) = wss*soil%surf_flux_sat + wsr*rice%surf_flux_sat
			methane_surf_aere_unsat(i) = wus*soil%surf_aere_unsat + wur*rice%surf_aere_unsat
			methane_surf_aere_sat(i) = wss*soil%surf_aere_sat + wsr*rice%surf_aere_sat
			methane_surf_ebul_unsat(i) = wus*soil%surf_ebul_unsat + wur*rice%surf_ebul_unsat
			methane_surf_ebul_sat(i) = wss*soil%surf_ebul_sat + wsr*rice%surf_ebul_sat
			methane_surf_diff_unsat(i) = wus*soil%surf_diff_unsat + wur*rice%surf_diff_unsat
			methane_surf_diff_sat(i) = wss*soil%surf_diff_sat + wsr*rice%surf_diff_sat
			methane_surf_diff_phys_unsat(i) = wus*soil%surf_diff_phys_unsat + wur*rice%surf_diff_phys_unsat
			methane_surf_diff_phys_sat(i) = wss*soil%surf_diff_phys_sat + wsr*rice%surf_diff_phys_sat
			methane_ebul_tot_unsat(i) = wus*soil%ebul_tot_unsat + wur*rice%ebul_tot_unsat
			methane_ebul_tot_sat(i) = wss*soil%ebul_tot_sat + wsr*rice%ebul_tot_sat
			methane_prod_tot_unsat(i) = wus*soil%prod_tot_unsat + wur*rice%prod_tot_unsat
			methane_prod_tot_sat(i) = wss*soil%prod_tot_sat + wsr*rice%prod_tot_sat
			methane_oxid_tot_unsat(i) = wus*soil%oxid_tot_unsat + wur*rice%oxid_tot_unsat
			methane_oxid_tot_sat(i) = wss*soil%oxid_tot_sat + wsr*rice%oxid_tot_sat
			co2_decomp_tot_unsat(i) = wus*soil%co2_decomp_tot_unsat + wur*rice%co2_decomp_tot_unsat
			co2_decomp_tot_sat(i) = wss*soil%co2_decomp_tot_sat + wsr*rice%co2_decomp_tot_sat
			co2_oxid_tot_unsat(i) = wus*soil%co2_oxid_tot_unsat + wur*rice%co2_oxid_tot_unsat
			co2_oxid_tot_sat(i) = wss*soil%co2_oxid_tot_sat + wsr*rice%co2_oxid_tot_sat
			co2_net_tot_unsat(i) = wus*soil%co2_net_tot_unsat + wur*rice%co2_net_tot_unsat
			co2_net_tot_sat(i) = wss*soil%co2_net_tot_sat + wsr*rice%co2_net_tot_sat
			totcol_methane_unsat(i) = wus*soil%totcol_unsat + wur*rice%totcol_unsat
			totcol_methane_sat(i) = wss*soil%totcol_sat + wsr*rice%totcol_sat
			grnd_methane_cond_unsat(i) = wus*soil%grnd_cond_unsat + wur*rice%grnd_cond_unsat
			grnd_methane_cond_sat(i) = wss*soil%grnd_cond_sat + wsr*rice%grnd_cond_sat

			methane_prod_depth(:,i) = ws*soil%prod + wr*rice%prod
			o2_decomp_depth(:,i) = ws*soil%o2_decomp + wr*rice%o2_decomp
			co2_decomp_depth(:,i) = ws*soil%co2_decomp + wr*rice%co2_decomp
			methane_oxid_depth(:,i) = ws*soil%oxid + wr*rice%oxid
			o2_oxid_depth(:,i) = ws*soil%o2_oxid + wr*rice%o2_oxid
			co2_oxid_depth(:,i) = ws*soil%co2_oxid + wr*rice%co2_oxid
			methane_aere_depth(:,i) = ws*soil%aere + wr*rice%aere
			methane_tran_depth(:,i) = ws*soil%tran + wr*rice%tran
			o2_aere_depth(:,i) = ws*soil%o2_aere + wr*rice%o2_aere
			co2_aere_depth(:,i) = ws*soil%co2_aere + wr*rice%co2_aere
			methane_ebul_depth(:,i) = ws*soil%ebul + wr*rice%ebul
			o2stress(:,i) = ws*soil%o2stress + wr*rice%o2stress
			methane_stress(:,i) = ws*soil%ch4stress + wr*rice%ch4stress
			methane_prod_depth_unsat(:,i) = wus*soil%prod_unsat + wur*rice%prod_unsat
			methane_prod_depth_sat(:,i) = wss*soil%prod_sat + wsr*rice%prod_sat
			o2_decomp_depth_unsat(:,i) = wus*soil%o2_decomp_unsat + wur*rice%o2_decomp_unsat
			o2_decomp_depth_sat(:,i) = wss*soil%o2_decomp_sat + wsr*rice%o2_decomp_sat
			co2_decomp_depth_unsat(:,i) = wus*soil%co2_decomp_unsat + wur*rice%co2_decomp_unsat
			co2_decomp_depth_sat(:,i) = wss*soil%co2_decomp_sat + wsr*rice%co2_decomp_sat
			methane_oxid_depth_unsat(:,i) = wus*soil%oxid_unsat + wur*rice%oxid_unsat
			methane_oxid_depth_sat(:,i) = wss*soil%oxid_sat + wsr*rice%oxid_sat
			o2_oxid_depth_unsat(:,i) = wus*soil%o2_oxid_unsat + wur*rice%o2_oxid_unsat
			o2_oxid_depth_sat(:,i) = wss*soil%o2_oxid_sat + wsr*rice%o2_oxid_sat
			co2_oxid_depth_unsat(:,i) = wus*soil%co2_oxid_unsat + wur*rice%co2_oxid_unsat
			co2_oxid_depth_sat(:,i) = wss*soil%co2_oxid_sat + wsr*rice%co2_oxid_sat
			methane_aere_depth_unsat(:,i) = wus*soil%aere_unsat + wur*rice%aere_unsat
			methane_aere_depth_sat(:,i) = wss*soil%aere_sat + wsr*rice%aere_sat
			methane_tran_depth_unsat(:,i) = wus*soil%tran_unsat + wur*rice%tran_unsat
			methane_tran_depth_sat(:,i) = wss*soil%tran_sat + wsr*rice%tran_sat
			o2_aere_depth_unsat(:,i) = wus*soil%o2_aere_unsat + wur*rice%o2_aere_unsat
			o2_aere_depth_sat(:,i) = wss*soil%o2_aere_sat + wsr*rice%o2_aere_sat
			co2_aere_depth_unsat(:,i) = wus*soil%co2_aere_unsat + wur*rice%co2_aere_unsat
			co2_aere_depth_sat(:,i) = wss*soil%co2_aere_sat + wsr*rice%co2_aere_sat
			methane_ebul_depth_unsat(:,i) = wus*soil%ebul_unsat + wur*rice%ebul_unsat
			methane_ebul_depth_sat(:,i) = wss*soil%ebul_sat + wsr*rice%ebul_sat
			o2stress_unsat(:,i) = wus*soil%o2stress_unsat + wur*rice%o2stress_unsat
			o2stress_sat(:,i) = wss*soil%o2stress_sat + wsr*rice%o2stress_sat
			methane_stress_unsat(:,i) = wus*soil%ch4stress_unsat + wur*rice%ch4stress_unsat
			methane_stress_sat(:,i) = wss*soil%ch4stress_sat + wsr*rice%ch4stress_sat

			conc_o2(:,i) = ws*soil%conc_o2 + wr*rice%conc_o2
			conc_methane(:,i) = ws*soil%conc_ch4 + wr*rice%conc_ch4
			conc_o2_unsat(:,i) = wus*soil%conc_o2_unsat + wur*rice%conc_o2_unsat
			conc_o2_sat(:,i) = wss*soil%conc_o2_sat + wsr*rice%conc_o2_sat
			conc_methane_unsat(:,i) = wus*soil%conc_ch4_unsat + wur*rice%conc_ch4_unsat
			conc_methane_sat(:,i) = wss*soil%conc_ch4_sat + wsr*rice%conc_ch4_sat
			layer_sat_lag(:,i) = ws*layer_sat_lag_component(:,METHANE_COMP_SOIL,i) + &
			   wr*layer_sat_lag_component(:,METHANE_COMP_RICE,i)
			annavg_agnpp(i) = ws*annavg_agnpp_component(METHANE_COMP_SOIL,i) + wr*annavg_agnpp_component(METHANE_COMP_RICE,i)
			annavg_bgnpp(i) = ws*annavg_bgnpp_component(METHANE_COMP_SOIL,i) + wr*annavg_bgnpp_component(METHANE_COMP_RICE,i)
			annavg_somhr(i) = ws*annavg_somhr_component(METHANE_COMP_SOIL,i) + wr*annavg_somhr_component(METHANE_COMP_RICE,i)
			annavg_finrw(i) = ws*annavg_finrw_component(METHANE_COMP_SOIL,i) + wr*annavg_finrw_component(METHANE_COMP_RICE,i)
			tempavg_agnpp(i) = ws*tempavg_agnpp_component(METHANE_COMP_SOIL,i) + wr*tempavg_agnpp_component(METHANE_COMP_RICE,i)
			tempavg_bgnpp(i) = ws*tempavg_bgnpp_component(METHANE_COMP_SOIL,i) + wr*tempavg_bgnpp_component(METHANE_COMP_RICE,i)
			annsum_counter(i) = ws*annsum_counter_component(METHANE_COMP_SOIL,i) + wr*annsum_counter_component(METHANE_COMP_RICE,i)
			tempavg_somhr(i) = ws*tempavg_somhr_component(METHANE_COMP_SOIL,i) + wr*tempavg_somhr_component(METHANE_COMP_RICE,i)
			tempavg_finrw(i) = ws*tempavg_finrw_component(METHANE_COMP_SOIL,i) + wr*tempavg_finrw_component(METHANE_COMP_RICE,i)
			fsat_bef(i) = ws*fsat_bef_component(METHANE_COMP_SOIL,i) + wr*fsat_bef_component(METHANE_COMP_RICE,i)
			finundated_lag(i) = ws*finundated_lag_component(METHANE_COMP_SOIL,i) + wr*finundated_lag_component(METHANE_COMP_RICE,i)

			methane_finundated(i) = ws*soil%finundated + wr*rice%finundated
			methane_soil_finundated(i) = soil%finundated_default
			methane_soil_zwt(i) = zwt
			methane_surf_flux_wetland(i) = 0._r8
			methane_surf_flux_lake(i) = 0._r8
			methane_surf_flux_soil(i) = ws*soil%surf_flux
			methane_surf_flux_rice(i) = wr*rice%surf_flux
			methane_surf_aere_soil(i) = ws*soil%surf_aere
			methane_surf_aere_rice(i) = wr*rice%surf_aere
			methane_surf_ebul_soil(i) = ws*soil%surf_ebul
			methane_surf_ebul_rice(i) = wr*rice%surf_ebul
			methane_surf_diff_soil(i) = ws*soil%surf_diff
			methane_surf_diff_rice(i) = wr*rice%surf_diff
			methane_prod_tot_soil(i) = ws*soil%prod_tot
			methane_prod_tot_rice(i) = wr*rice%prod_tot
			methane_oxid_tot_soil(i) = ws*soil%oxid_tot
			methane_oxid_tot_rice(i) = wr*rice%oxid_tot
		END SUBROUTINE aggregate_methane_columns

		SUBROUTINE repartition_methane_column_state(ipatch, old_fraction, new_fraction)
			integer, intent(in) :: ipatch
			real(r8), intent(in) :: old_fraction, new_fraction
			real(r8) :: rold, rnew

			rold = min(max(old_fraction, 0._r8), 1._r8)
			rnew = min(max(new_fraction, 0._r8), 1._r8)
			IF (abs(rnew-rold) <= 1.e-14_r8) RETURN
			CALL repartition_phase_state(conc_o2_unsat_component(:,METHANE_COMP_SOIL,ipatch), &
			   conc_o2_sat_component(:,METHANE_COMP_SOIL,ipatch), &
			   conc_o2_unsat_component(:,METHANE_COMP_RICE,ipatch), &
			   conc_o2_sat_component(:,METHANE_COMP_RICE,ipatch), &
			   fsat_bef_component(METHANE_COMP_SOIL,ipatch), &
			   fsat_bef_component(METHANE_COMP_RICE,ipatch), rold, rnew)
			CALL repartition_phase_state(conc_methane_unsat_component(:,METHANE_COMP_SOIL,ipatch), &
			   conc_methane_sat_component(:,METHANE_COMP_SOIL,ipatch), &
			   conc_methane_unsat_component(:,METHANE_COMP_RICE,ipatch), &
			   conc_methane_sat_component(:,METHANE_COMP_RICE,ipatch), &
			   fsat_bef_component(METHANE_COMP_SOIL,ipatch), &
			   fsat_bef_component(METHANE_COMP_RICE,ipatch), rold, rnew)
			CALL repartition_vector(layer_sat_lag_component(:,METHANE_COMP_SOIL,ipatch), &
			   layer_sat_lag_component(:,METHANE_COMP_RICE,ipatch), rold, rnew)
			CALL repartition_scalar(annavg_agnpp_component(METHANE_COMP_SOIL,ipatch), &
			   annavg_agnpp_component(METHANE_COMP_RICE,ipatch), rold, rnew)
			CALL repartition_scalar(annavg_bgnpp_component(METHANE_COMP_SOIL,ipatch), &
			   annavg_bgnpp_component(METHANE_COMP_RICE,ipatch), rold, rnew)
			CALL repartition_scalar(annavg_somhr_component(METHANE_COMP_SOIL,ipatch), &
			   annavg_somhr_component(METHANE_COMP_RICE,ipatch), rold, rnew)
			CALL repartition_scalar(annavg_finrw_component(METHANE_COMP_SOIL,ipatch), &
			   annavg_finrw_component(METHANE_COMP_RICE,ipatch), rold, rnew)
			CALL repartition_scalar(tempavg_agnpp_component(METHANE_COMP_SOIL,ipatch), &
			   tempavg_agnpp_component(METHANE_COMP_RICE,ipatch), rold, rnew)
			CALL repartition_scalar(tempavg_bgnpp_component(METHANE_COMP_SOIL,ipatch), &
			   tempavg_bgnpp_component(METHANE_COMP_RICE,ipatch), rold, rnew)
			CALL repartition_scalar(annsum_counter_component(METHANE_COMP_SOIL,ipatch), &
			   annsum_counter_component(METHANE_COMP_RICE,ipatch), rold, rnew)
			CALL repartition_scalar(tempavg_somhr_component(METHANE_COMP_SOIL,ipatch), &
			   tempavg_somhr_component(METHANE_COMP_RICE,ipatch), rold, rnew)
			CALL repartition_scalar(tempavg_finrw_component(METHANE_COMP_SOIL,ipatch), &
			   tempavg_finrw_component(METHANE_COMP_RICE,ipatch), rold, rnew)
			CALL repartition_scalar(fsat_bef_component(METHANE_COMP_SOIL,ipatch), &
			   fsat_bef_component(METHANE_COMP_RICE,ipatch), rold, rnew)
			CALL repartition_scalar(finundated_lag_component(METHANE_COMP_SOIL,ipatch), &
			   finundated_lag_component(METHANE_COMP_RICE,ipatch), rold, rnew)
		END SUBROUTINE repartition_methane_column_state

		SUBROUTINE repartition_phase_state(soil_unsat, soil_sat, rice_unsat, rice_sat, &
		   soil_fsat, rice_fsat, old_fraction, new_fraction)
			real(r8), intent(inout) :: soil_unsat(:), soil_sat(:), rice_unsat(:), rice_sat(:)
			real(r8), intent(in) :: soil_fsat, rice_fsat, old_fraction, new_fraction
			real(r8) :: hs, hr, delta, mixed
			integer :: j

			hs = soil_fsat
			IF (ieee_is_nan(hs) .or. abs(hs) >= 0.5_r8*abs(spval) .or. hs < 0._r8 .or. hs > 1._r8) hs = 0.5_r8
			hr = rice_fsat
			IF (ieee_is_nan(hr) .or. abs(hr) >= 0.5_r8*abs(spval) .or. hr < 0._r8 .or. hr > 1._r8) hr = 0.5_r8
			IF (new_fraction > old_fraction) THEN
				delta = new_fraction-old_fraction
				DO j = 1, size(soil_unsat)
					mixed = (old_fraction*(hr*rice_sat(j)+(1._r8-hr)*rice_unsat(j)) + &
					   delta*(hs*soil_sat(j)+(1._r8-hs)*soil_unsat(j))) / new_fraction
					rice_unsat(j) = mixed
					rice_sat(j) = mixed
				ENDDO
			ELSE
				delta = old_fraction-new_fraction
				DO j = 1, size(soil_unsat)
					mixed = ((1._r8-old_fraction)*(hs*soil_sat(j)+(1._r8-hs)*soil_unsat(j)) + &
					   delta*(hr*rice_sat(j)+(1._r8-hr)*rice_unsat(j))) / (1._r8-new_fraction)
					soil_unsat(j) = mixed
					soil_sat(j) = mixed
				ENDDO
			ENDIF
		END SUBROUTINE repartition_phase_state

		SUBROUTINE repartition_vector(soil, rice, old_fraction, new_fraction)
			real(r8), intent(inout) :: soil(:), rice(:)
			real(r8), intent(in) :: old_fraction, new_fraction
			integer :: j
			DO j = 1, size(soil)
				CALL repartition_scalar(soil(j), rice(j), old_fraction, new_fraction)
			ENDDO
		END SUBROUTINE repartition_vector

		SUBROUTINE repartition_scalar(soil, rice, old_fraction, new_fraction)
			real(r8), intent(inout) :: soil, rice
			real(r8), intent(in) :: old_fraction, new_fraction
			real(r8) :: delta
			logical :: soil_valid, rice_valid

			soil_valid = .not. ieee_is_nan(soil) .and. abs(soil) < 0.5_r8*abs(spval)
			rice_valid = .not. ieee_is_nan(rice) .and. abs(rice) < 0.5_r8*abs(spval)
			IF (new_fraction > old_fraction) THEN
				IF (.not. soil_valid) RETURN
				IF (.not. rice_valid .or. old_fraction <= 1.e-14_r8) THEN
					rice = soil
				ELSE
					delta = new_fraction - old_fraction
					rice = (old_fraction*rice + delta*soil) / new_fraction
				ENDIF
			ELSE
				IF (.not. rice_valid) RETURN
				IF (.not. soil_valid .or. old_fraction >= 1._r8-1.e-14_r8) THEN
					soil = rice
				ELSE
					delta = old_fraction - new_fraction
					soil = ((1._r8-old_fraction)*soil + delta*rice) / (1._r8-new_fraction)
				ENDIF
			ENDIF
		END SUBROUTINE repartition_scalar

	END SUBROUTINE methane_driver

END MODULE MOD_MT_Methane_Driver
#endif
