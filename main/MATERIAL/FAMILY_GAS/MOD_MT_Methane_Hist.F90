#include <define.h>

#if (defined MaterialTransport) && (defined BGC)
MODULE MOD_MT_Methane_Hist

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: methane_reactive_history

   integer, save :: methane_history_match_count = 0

CONTAINS

   SUBROUTINE methane_reactive_history (file_hist, itime_in_file, sumarea, filter, &
      nl_soil, forcing_has_missing_value, forcmask_pch)
   USE MOD_MT_Hist, only: HistForm, write_history_variable_2d, write_history_variable_3d, &
      write_history_mt_ratio_2d
   USE MOD_MT_Defs, only: material_lower
   USE MOD_HistGridded
   USE MOD_DataType
   USE MOD_LandPatch, only: numpatch
   USE MOD_SPMD_Task, only: p_is_worker
   USE MOD_Vars_TimeInvariants, only: patchtype, patchmask
   USE MOD_MT_Methane_BgcLink, only: methane_patch_active_mask
            USE MOD_MT_Methane_Registry, only: igas_ch4
	            USE MOD_MT_Methane_Const,    only: DEF_METHANE
            USE MOD_MT_Methane_AccFlux,  only: &
                 a_net_methane, a_co2_hr_methane, a_methane_prod_depth, a_o2_decomp_depth, a_co2_decomp_depth, &
                 a_methane_oxid_depth, a_o2_oxid_depth, a_co2_oxid_depth, &
                 a_methane_aere_depth, a_methane_tran_depth, a_o2_aere_depth, a_co2_aere_depth, &
                 a_methane_ebul_depth, a_o2stress, a_methane_stress, &
                 a_methane_surf_flux_tot, a_methane_surf_flux_tot_phys, a_methane_surf_aere, a_methane_surf_ebul, &
                    a_methane_surf_diff, a_methane_surf_diff_phys, a_methane_balance_residual, a_methane_ch4_clip_credit, &
                    a_methane_surf_aere_soil, a_methane_surf_aere_rice, &
                    a_methane_surf_ebul_soil, a_methane_surf_ebul_rice, &
                    a_methane_surf_diff_soil, a_methane_surf_diff_rice, &
                    a_o2_cap_loss, a_o2_cap_gain, &
                    a_methane_ebul_tot, a_methane_prod_tot, a_methane_oxid_tot, &
                    a_methane_prod_tot_soil, a_methane_prod_tot_rice, &
                    a_methane_oxid_tot_soil, a_methane_oxid_tot_rice, &
                 a_co2_decomp_tot, a_co2_oxid_tot, a_co2_aere_tot, a_co2_net_tot, &
                 a_totcol_methane, a_grnd_methane_cond, a_conc_o2, a_conc_methane, &
                    a_net_methane_unsat, a_net_methane_sat, &
                    a_methane_prod_depth_unsat, a_methane_prod_depth_sat, &
                       a_o2_decomp_depth_unsat, a_o2_decomp_depth_sat, &
                       a_co2_decomp_depth_unsat, a_co2_decomp_depth_sat, &
                       a_methane_oxid_depth_unsat, a_methane_oxid_depth_sat, &
                       a_o2_oxid_depth_unsat, a_o2_oxid_depth_sat, &
                       a_co2_oxid_depth_unsat, a_co2_oxid_depth_sat, &
                       a_methane_aere_depth_unsat, a_methane_aere_depth_sat, &
                       a_methane_tran_depth_unsat, a_methane_tran_depth_sat, &
                       a_o2_aere_depth_unsat, a_o2_aere_depth_sat, &
                       a_co2_aere_depth_unsat, a_co2_aere_depth_sat, &
                       a_methane_ebul_depth_unsat, a_methane_ebul_depth_sat, &
                       a_o2stress_unsat, a_o2stress_sat, &
                       a_methane_stress_unsat, a_methane_stress_sat, &
                       a_methane_surf_flux_tot_unsat, a_methane_surf_flux_tot_sat, &
                       a_methane_surf_aere_unsat, a_methane_surf_aere_sat, &
                       a_methane_surf_ebul_unsat, a_methane_surf_ebul_sat, &
                       a_methane_surf_diff_unsat, a_methane_surf_diff_sat, &
                       a_methane_ebul_tot_unsat, a_methane_ebul_tot_sat, &
                       a_methane_prod_tot_unsat, a_methane_prod_tot_sat, &
                       a_methane_oxid_tot_unsat, a_methane_oxid_tot_sat, &
                       a_co2_decomp_tot_unsat, a_co2_decomp_tot_sat, &
                       a_co2_oxid_tot_unsat, a_co2_oxid_tot_sat, &
                       a_co2_net_tot_unsat, a_co2_net_tot_sat, &
                       a_totcol_methane_unsat, a_totcol_methane_sat, &
                       a_grnd_methane_cond_unsat, a_grnd_methane_cond_sat, &
                    a_conc_o2_unsat, a_conc_o2_sat, &
                    a_conc_methane_unsat, a_conc_methane_sat, &
                 a_methane_prod_depth_lake, a_methane_oxid_depth_lake, a_methane_ebul_depth_lake, &
                 a_co2_decomp_depth_lake, a_co2_oxid_depth_lake, &
                 a_methane_surf_ebul_lake, a_methane_surf_diff_lake, a_methane_surf_flux_tot_lake, &
                 a_methane_prod_tot_lake, a_methane_oxid_tot_lake, a_methane_ebul_tot_lake, &
                 a_co2_decomp_tot_lake, a_co2_oxid_tot_lake, a_co2_net_tot_lake, &
                 a_totcol_methane_lake, a_grnd_methane_cond_lake, &
                 a_conc_o2_lake, a_conc_methane_lake, &
                 a_lake_water_ch4_stock, a_lake_water_o2_stock, &
                 a_lake_water_ch4_oxid, a_lake_sed_ch4_flux, &
                 a_lake_sed_o2_flux, a_lake_air_o2_flux, &
                    a_forc_pmethanem, a_layer_sat_lag, &
                    a_annavg_agnpp, a_annavg_bgnpp, a_annavg_somhr, a_annavg_finrw, &
                    a_methane_dfsat_tot, a_f_h2osfc, &
                    a_methane_finundated, a_methane_soil_finundated, a_methane_soil_zwt, &
                    a_f_inund_flood_patch, a_f_inund_flood_depth_patch, a_wetland_frac_per_patch, &
                    a_methane_surf_flux_wetland, a_methane_surf_flux_soil, &
                    a_methane_surf_flux_lake, a_methane_surf_flux_rice, &
                    a_methane_rice_fraction, &
                    a_B_methanogen, a_B_methanotroph, &
                 a_B_methanogen_dormant, a_B_methanotroph_dormant, &
                 a_f_T_methanogen, a_f_S_methanogen, a_f_O2_methanogen, a_f_T_methanotroph, &
                 a_methanogen_growth_rate, a_methanotroph_growth_rate, &
                 a_microbial_prod_potential, a_microbial_oxid_potential, &
                 a_methane_acc_num, a_methane_acc_num_unsat, a_methane_acc_num_sat, &
                 a_methane_acc_num_lake, a_methane_acc_num_extra, a_methane_acc_num_microbe, &
                 a_annavg_finrw_acc_num
   IMPLICIT NONE

   character(len=*), intent(in) :: file_hist
   integer, intent(in) :: itime_in_file
   type(block_data_real8_2d), intent(inout) :: sumarea
	   logical, intent(inout) :: filter(:)
	   integer, intent(in) :: nl_soil
	   logical, intent(in) :: forcing_has_missing_value
	   logical, intent(in) :: forcmask_pch(:)
		   real(r8), allocatable :: hist_ch4_active_without_lake(:)
		   real(r8), allocatable :: hist_ch4_global_with_lake(:)
		   real(r8), allocatable :: hist_ch4_global_phys_with_lake(:)
		   real(r8), allocatable :: hist_ch4_global_balance_with_lake(:)
		   real(r8), allocatable :: hist_ch4_global_clip_credit_with_lake(:)
		   real(r8), allocatable :: hist_ch4_land_total_surf(:)
		   real(r8), allocatable :: hist_ch4_land_flux_wetland(:)
		   real(r8), allocatable :: hist_ch4_land_flux_soil(:)
		   real(r8), allocatable :: hist_ch4_land_flux_lake(:)
		   real(r8), allocatable :: hist_ch4_land_flux_rice(:)
		   real(r8), allocatable :: hist_ch4_lake_flux_acc(:)
		   real(r8), allocatable :: hist_ch4_rice_flux_mean(:)
		   real(r8), allocatable :: hist_ch4_rice_area_frac(:)
		   real(r8), allocatable :: hist_ch4_acc_one(:)
		   logical, allocatable :: filter_active_without_lake(:)
		   logical, allocatable :: filter_all_land(:)
		   integer :: ipatch
		   logical :: need_active_without_lake, need_global_with_lake
		   logical :: need_global_phys_with_lake, need_global_balance_with_lake
		   logical :: need_global_clip_credit_with_lake
			   logical :: need_lake_intensive, need_rice_intensive, need_derived_ch4
			   logical :: need_land_flux_total, need_land_flux_wetland, need_land_flux_soil
			   logical :: need_land_flux_lake, need_land_flux_rice, need_land_flux_split
			   logical :: need_lake_history, need_extra_history, need_soil_history
			   logical :: need_active_history

	            IF (igas_ch4 <= 0) RETURN
	            IF (.not. DEF_METHANE%write_ch4_history) RETURN
	            IF (methane_history_selector_is_off(DEF_METHANE%ch4_history_vars)) RETURN
	            methane_history_match_count = 0
	            need_active_history = .not. methane_history_all_land_only(DEF_METHANE%ch4_history_vars)
	            IF (igas_ch4 > 0) THEN
               need_active_without_lake = mhist_on('f_methane_surf_flux_active_total_without_lake')
               need_global_with_lake    = mhist_on('f_methane_surf_flux_global_total_with_lake')
               need_global_phys_with_lake = mhist_on('f_methane_surf_flux_global_phys_with_lake')
               need_global_balance_with_lake = mhist_on('f_methane_balance_residual_global_with_lake')
               need_global_clip_credit_with_lake = mhist_on('f_methane_ch4_clip_credit_global_with_lake')
               need_lake_intensive      = mhist_on('f_methane_surf_flux_lake_intensive')
               need_rice_intensive      = mhist_on('f_methane_surf_flux_rice_intensive')
               need_land_flux_total     = mhist_on('f_methane_surf_flux_tot')
               need_land_flux_wetland   = mhist_on('f_methane_surf_flux_wetland')
               need_land_flux_soil      = mhist_on('f_methane_surf_flux_soil')
               need_land_flux_lake      = mhist_on('f_methane_surf_flux_lake')
               need_land_flux_rice      = mhist_on('f_methane_surf_flux_rice')
               need_land_flux_split = need_land_flux_total .or. need_land_flux_wetland .or. &
                                      need_land_flux_soil .or. need_land_flux_lake .or. &
                                      need_land_flux_rice
	               ! Keep one selector-function reference per statement: mhist_on
	               ! records whether the configured selector matched a real write.
	               need_lake_history = need_lake_intensive
	               need_lake_history = mhist_on('f_methane_prod_depth_lake') .or. need_lake_history
	               need_lake_history = mhist_on('f_methane_oxid_depth_lake') .or. need_lake_history
	               need_lake_history = mhist_on('f_co2_decomp_depth_lake') .or. need_lake_history
	               need_lake_history = mhist_on('f_co2_oxid_depth_lake') .or. need_lake_history
	               need_lake_history = mhist_on('f_methane_ebul_depth_lake') .or. need_lake_history
	               need_lake_history = mhist_on('f_methane_surf_ebul_lake') .or. need_lake_history
	               need_lake_history = mhist_on('f_methane_surf_diff_lake') .or. need_lake_history
	               need_lake_history = mhist_on('f_methane_surf_flux_tot_lake') .or. need_lake_history
	               need_lake_history = mhist_on('f_methane_prod_tot_lake') .or. need_lake_history
	               need_lake_history = mhist_on('f_methane_oxid_tot_lake') .or. need_lake_history
	               need_lake_history = mhist_on('f_methane_ebul_tot_lake') .or. need_lake_history
	               need_lake_history = mhist_on('f_co2_decomp_tot_lake') .or. need_lake_history
	               need_lake_history = mhist_on('f_co2_oxid_tot_lake') .or. need_lake_history
	               need_lake_history = mhist_on('f_co2_net_tot_lake') .or. need_lake_history
	               need_lake_history = mhist_on('f_totcol_methane_lake') .or. need_lake_history
	               need_lake_history = mhist_on('f_grnd_methane_cond_lake') .or. need_lake_history
	               need_lake_history = mhist_on('f_conc_o2_lake') .or. need_lake_history
	               need_lake_history = mhist_on('f_conc_methane_lake') .or. need_lake_history
	               need_lake_history = mhist_on('f_lake_water_ch4_stock') .or. need_lake_history
	               need_lake_history = mhist_on('f_lake_water_o2_stock') .or. need_lake_history
	               need_lake_history = mhist_on('f_lake_water_ch4_oxid') .or. need_lake_history
	               need_lake_history = mhist_on('f_lake_sed_ch4_flux') .or. need_lake_history
	               need_lake_history = mhist_on('f_lake_sed_o2_flux') .or. need_lake_history
	               need_lake_history = mhist_on('f_lake_air_o2_flux') .or. need_lake_history

	               need_extra_history = need_rice_intensive
	               need_extra_history = mhist_on('f_forc_pmethanem') .or. need_extra_history
	               need_extra_history = mhist_on('f_layer_sat_lag') .or. need_extra_history
	               need_extra_history = mhist_on('f_annavg_agnpp') .or. need_extra_history
	               need_extra_history = mhist_on('f_annavg_bgnpp') .or. need_extra_history
	               need_extra_history = mhist_on('f_annavg_somhr') .or. need_extra_history
	               need_extra_history = mhist_on('f_annavg_'//'finrw') .or. need_extra_history
	               need_extra_history = mhist_on('f_methane_dfsat_tot') .or. need_extra_history
	               need_extra_history = mhist_on('f_f_h2osfc') .or. need_extra_history
	               need_extra_history = mhist_on('f_methane_finundated') .or. need_extra_history
	               need_extra_history = mhist_on('f_inund_flood_patch') .or. need_extra_history
	               need_extra_history = mhist_on('f_inund_flood_depth_patch') .or. need_extra_history
	               need_extra_history = mhist_on('f_wetland_frac_patch') .or. need_extra_history

	               need_soil_history = mhist_on('f_methane_soil_zwt')
	               need_soil_history = mhist_on('f_methane_soil_finundated') .or. need_soil_history
               need_derived_ch4 = need_active_without_lake .or. need_global_with_lake .or. &
                                  need_global_phys_with_lake .or. need_global_balance_with_lake .or. &
                                  need_global_clip_credit_with_lake .or. &
                                  need_lake_intensive .or. need_rice_intensive .or. &
                                  need_land_flux_split
	               ! Outer filter (patchtype==0) excludes wetland (patchtype==2) and
	               ! lake (patchtype==4) — the only patches where the CH4 module
               ! actually runs.  Without overriding, every CH4 grid cell aggregates
               ! to zero because only soil patches are summed.  Re-derive a
               ! soil+wetland filter for the main/sat/unsat writes here, then
               ! swap to a lake-only filter before the lake section below.
                  ! Active CH4 mask mirrors the driver gate at
                  ! CoLMDRIVER.F90:381-384 — patchtype==2 and patchtype==0 are
                  ! active by default.  Urban (1) and ice (3) are NEVER active
                  ! even though "<=2" would include urban.
               IF ((p_is_worker) .and. (numpatch > 0)) THEN
                  ! Use BgcLink helper so the history filter matches the
                  ! driver gate + rice paddy override exactly.  Same call
                  ! site pattern as MOD_MT_Methane_AccFlux.F90.
                  CALL methane_patch_active_mask (numpatch, DEF_METHANE%only_wetland, &
                                                  DEF_METHANE%enable_rice_paddy, &
                                                  patchtype, filter)
                  filter = filter .and. patchmask
                  IF (forcing_has_missing_value) filter = filter .and. forcmask_pch
               ENDIF
		               IF (need_active_history .and. HistForm == 'Gridded') THEN
		                  CALL mp2g_hist%get_sumarea (sumarea, filter)
		               ENDIF

		               ! Allocate only the derived vectors selected for this write.
		               ! In particular, core history avoids lake-intensive and rice-
		               ! intensive temporaries while retaining its land-flux split.
               IF (need_derived_ch4) THEN
		               IF (need_active_without_lake) allocate (hist_ch4_active_without_lake(numpatch))
		               IF (need_global_with_lake) allocate (hist_ch4_global_with_lake(numpatch))
		               IF (need_global_phys_with_lake) allocate (hist_ch4_global_phys_with_lake(numpatch))
		               IF (need_global_balance_with_lake) allocate (hist_ch4_global_balance_with_lake(numpatch))
		               IF (need_global_clip_credit_with_lake) allocate (hist_ch4_global_clip_credit_with_lake(numpatch))
			               IF (need_land_flux_total) allocate (hist_ch4_land_total_surf(numpatch))
			               IF (need_land_flux_wetland) allocate (hist_ch4_land_flux_wetland(numpatch))
			               IF (need_land_flux_soil) allocate (hist_ch4_land_flux_soil(numpatch))
			               IF (need_land_flux_lake) allocate (hist_ch4_land_flux_lake(numpatch))
			               IF (need_land_flux_rice) allocate (hist_ch4_land_flux_rice(numpatch))
			               IF (need_lake_intensive) allocate (hist_ch4_lake_flux_acc(numpatch))
			               IF (need_rice_intensive) THEN
			                  allocate (hist_ch4_rice_flux_mean(numpatch))
			                  allocate (hist_ch4_rice_area_frac(numpatch))
			               ENDIF
		               IF (need_land_flux_split .or. need_global_with_lake .or. &
		                   need_global_phys_with_lake .or. need_global_balance_with_lake .or. &
		                   need_global_clip_credit_with_lake) &
		                  allocate (hist_ch4_acc_one(numpatch))
			               allocate (filter_active_without_lake(numpatch))
			               allocate (filter_all_land(numpatch))
		               IF (need_active_without_lake) hist_ch4_active_without_lake(:) = 0._r8
		               IF (need_global_with_lake) hist_ch4_global_with_lake(:) = 0._r8
		               IF (need_global_phys_with_lake) hist_ch4_global_phys_with_lake(:) = 0._r8
		               IF (need_global_balance_with_lake) hist_ch4_global_balance_with_lake(:) = 0._r8
		               IF (need_global_clip_credit_with_lake) hist_ch4_global_clip_credit_with_lake(:) = 0._r8
			               IF (need_land_flux_total) hist_ch4_land_total_surf(:) = 0._r8
			               IF (need_land_flux_wetland) hist_ch4_land_flux_wetland(:) = 0._r8
			               IF (need_land_flux_soil) hist_ch4_land_flux_soil(:) = 0._r8
			               IF (need_land_flux_lake) hist_ch4_land_flux_lake(:) = 0._r8
			               IF (need_land_flux_rice) hist_ch4_land_flux_rice(:) = 0._r8
			               IF (need_lake_intensive) hist_ch4_lake_flux_acc(:) = 0._r8
			               IF (need_rice_intensive) THEN
			                  hist_ch4_rice_flux_mean(:) = 0._r8
			                  hist_ch4_rice_area_frac(:) = 0._r8
			               ENDIF
		               IF (need_land_flux_split .or. need_global_with_lake .or. &
		                   need_global_phys_with_lake .or. need_global_balance_with_lake .or. &
		                   need_global_clip_credit_with_lake) &
		                  hist_ch4_acc_one(:) = 1._r8
	               filter_active_without_lake(:) = .false.
	               filter_all_land(:) = .false.
	               IF ((p_is_worker) .and. (numpatch > 0)) THEN
	                  CALL methane_patch_active_mask (numpatch, DEF_METHANE%only_wetland, &
	                                                  DEF_METHANE%enable_rice_paddy, &
	                                                  patchtype, filter_active_without_lake)
	                  filter_active_without_lake = filter_active_without_lake .and. patchmask
	                  filter_all_land = (patchtype < 99) .and. patchmask
	                  IF (forcing_has_missing_value) THEN
	                     filter_active_without_lake = filter_active_without_lake .and. forcmask_pch
	                     filter_all_land = filter_all_land .and. forcmask_pch
	                  ENDIF
			                  IF (need_active_without_lake) hist_ch4_active_without_lake = a_methane_surf_flux_tot
			                  IF (need_lake_intensive) hist_ch4_lake_flux_acc = a_methane_surf_flux_lake
			                  DO ipatch = 1, numpatch
			                     IF (filter_active_without_lake(ipatch) .and. &
			                         (a_methane_acc_num(ipatch) > 0._r8)) THEN
		                        IF (need_global_with_lake) hist_ch4_global_with_lake(ipatch) = &
		                           a_methane_surf_flux_tot(ipatch) / a_methane_acc_num(ipatch)
		                        IF (need_global_phys_with_lake) hist_ch4_global_phys_with_lake(ipatch) = &
		                           a_methane_surf_flux_tot_phys(ipatch) / a_methane_acc_num(ipatch)
		                        IF (need_global_balance_with_lake) hist_ch4_global_balance_with_lake(ipatch) = &
		                           a_methane_balance_residual(ipatch) / a_methane_acc_num(ipatch)
		                        IF (need_global_clip_credit_with_lake) hist_ch4_global_clip_credit_with_lake(ipatch) = &
		                           a_methane_ch4_clip_credit(ipatch) / a_methane_acc_num(ipatch)
			                        IF (need_land_flux_total) hist_ch4_land_total_surf(ipatch) = &
			                           a_methane_surf_flux_tot(ipatch) / a_methane_acc_num(ipatch)
			                        IF (need_land_flux_wetland) hist_ch4_land_flux_wetland(ipatch) = &
			                           a_methane_surf_flux_wetland(ipatch) / a_methane_acc_num(ipatch)
			                        IF (need_land_flux_soil) hist_ch4_land_flux_soil(ipatch) = &
			                           a_methane_surf_flux_soil(ipatch) / a_methane_acc_num(ipatch)
			                        IF (need_land_flux_rice) hist_ch4_land_flux_rice(ipatch) = &
			                           a_methane_surf_flux_rice(ipatch) / a_methane_acc_num(ipatch)
						IF (need_rice_intensive) THEN
						   ! Preserve the raw same-window sums.  The history ratio
						   ! writer maps numerator and denominator separately before
						   ! dividing, yielding a paddy-area-weighted grid mean.
						   hist_ch4_rice_flux_mean(ipatch) = a_methane_surf_flux_rice(ipatch)
						   hist_ch4_rice_area_frac(ipatch) = a_methane_rice_fraction(ipatch)
						ENDIF
			                     ENDIF
			                     IF ((patchtype(ipatch) == 4) .and. DEF_METHANE%allowlakeprod .and. &
			                         filter_all_land(ipatch) .and. &
			                         (a_methane_acc_num_lake(ipatch) > 0._r8)) THEN
		                        IF (need_global_with_lake) hist_ch4_global_with_lake(ipatch) = &
		                           a_methane_surf_flux_tot_lake(ipatch) / a_methane_acc_num_lake(ipatch)
		                        IF (need_global_phys_with_lake) hist_ch4_global_phys_with_lake(ipatch) = &
		                           a_methane_surf_flux_tot_phys(ipatch) / a_methane_acc_num_lake(ipatch)
		                        IF (need_global_balance_with_lake) hist_ch4_global_balance_with_lake(ipatch) = &
		                           a_methane_balance_residual(ipatch) / a_methane_acc_num_lake(ipatch)
		                        IF (need_global_clip_credit_with_lake) hist_ch4_global_clip_credit_with_lake(ipatch) = &
		                           a_methane_ch4_clip_credit(ipatch) / a_methane_acc_num_lake(ipatch)
			                        IF (need_land_flux_total) hist_ch4_land_total_surf(ipatch) = &
			                           a_methane_surf_flux_tot_lake(ipatch) / a_methane_acc_num_lake(ipatch)
			                        IF (need_land_flux_lake) hist_ch4_land_flux_lake(ipatch) = &
			                           a_methane_surf_flux_lake(ipatch) / a_methane_acc_num_lake(ipatch)
			                     ENDIF
			                  ENDDO
			               ENDIF
               ENDIF

	               CALL write_history_variable_2d (mhist_on('f_net_methane'), a_net_methane, file_hist, &
                  'f_net_methane', itime_in_file, sumarea, filter, &
                  'CH4 oxidation minus production; applied to soil/wetland BGC CO2 respiration', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_co2_hr_methane'), a_co2_hr_methane, file_hist, &
                  'f_co2_hr_methane', itime_in_file, sumarea, filter, &
                  'CO2-C heterotrophic respiration after the CH4 carbon split (host f_hr is not changed)', &
                  'gC/m2/s', acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_methane_prod_depth'), a_methane_prod_depth, file_hist, &
                  'f_methane_prod_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 production per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_o2_decomp_depth'), a_o2_decomp_depth, file_hist, &
                  'f_o2_decomp_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'O2 consumption during decomposition per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_co2_decomp_depth'), a_co2_decomp_depth, file_hist, &
                  'f_co2_decomp_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'diagnostic CO2 decomp/methanogenesis before O2-stress scaling per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_methane_oxid_depth'), a_methane_oxid_depth, file_hist, &
                  'f_methane_oxid_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 oxidation per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_o2_oxid_depth'), a_o2_oxid_depth, file_hist, &
                  'f_o2_oxid_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'O2 consumption via oxidation per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_co2_oxid_depth'), a_co2_oxid_depth, file_hist, &
                  'f_co2_oxid_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CO2 production from CH4 oxidation per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_methane_aere_depth'), a_methane_aere_depth, file_hist, &
                  'f_methane_aere_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 loss via aerenchyma per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_methane_tran_depth'), a_methane_tran_depth, file_hist, &
                  'f_methane_tran_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 loss via transpiration per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_o2_aere_depth'), a_o2_aere_depth, file_hist, &
                  'f_o2_aere_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'O2 gain via aerenchyma per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_co2_aere_depth'), a_co2_aere_depth, file_hist, &
                  'f_co2_aere_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'diagnosed CO2 aerenchyma flux per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_methane_ebul_depth'), a_methane_ebul_depth, file_hist, &
                  'f_methane_ebul_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 loss via ebullition per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_o2stress'), a_o2stress, file_hist, &
                  'f_o2stress', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'O2 availability/demand ratio', '-', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_methane_stress'), a_methane_stress, file_hist, &
                  'f_methane_stress', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 availability/sinks ratio', '-', &
	                  acc_num=a_methane_acc_num)
	                  CALL write_history_variable_2d (mhist_on('f_methane_surf_flux_tot_active'), &
	                     a_methane_surf_flux_tot, file_hist, &
	                     'f_methane_surf_flux_tot_active', itime_in_file, sumarea, filter, &
	                     'corrected CH4 total surface flux; active-CH4-area mean; lake excluded', 'mol/m2/s', &
	                     acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_surf_flux_tot_phys'), a_methane_surf_flux_tot_phys, file_hist, &
                     'f_methane_surf_flux_tot_phys', itime_in_file, sumarea, filter, &
                     'CH4 physical total surface flux before clip/residual; active-CH4-area mean; lake excluded', 'mol/m2/s', &
                     acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_surf_aere'), a_methane_surf_aere, file_hist, &
                  'f_methane_surf_aere', itime_in_file, sumarea, filter, &
                  'CH4 column aerenchyma flux', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_surf_ebul'), a_methane_surf_ebul, file_hist, &
                  'f_methane_surf_ebul', itime_in_file, sumarea, filter, &
                  'CH4 surface ebullition flux', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_surf_aere_soil'), a_methane_surf_aere_soil, file_hist, &
                  'f_methane_surf_aere_soil', itime_in_file, sumarea, filter, &
                  'non-rice CH4 aerenchyma flux contribution per active area', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_surf_aere_rice'), a_methane_surf_aere_rice, file_hist, &
                  'f_methane_surf_aere_rice', itime_in_file, sumarea, filter, &
                  'rice CH4 aerenchyma flux contribution per active area', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_surf_ebul_soil'), a_methane_surf_ebul_soil, file_hist, &
                  'f_methane_surf_ebul_soil', itime_in_file, sumarea, filter, &
                  'non-rice CH4 ebullition flux contribution per active area', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_surf_ebul_rice'), a_methane_surf_ebul_rice, file_hist, &
                  'f_methane_surf_ebul_rice', itime_in_file, sumarea, filter, &
                  'rice CH4 ebullition flux contribution per active area', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
                     CALL write_history_variable_2d (mhist_on('f_methane_surf_diff'), a_methane_surf_diff, file_hist, &
                        'f_methane_surf_diff', itime_in_file, sumarea, filter, &
                        'CH4 diffusive surface flux plus inundated-fraction and closure adjustments', 'mol/m2/s', &
                        acc_num=a_methane_acc_num)
                     CALL write_history_variable_2d (mhist_on('f_methane_surf_diff_soil'), &
                        a_methane_surf_diff_soil, file_hist, &
                        'f_methane_surf_diff_soil', itime_in_file, sumarea, filter, &
                        'non-rice CH4 diffusive flux contribution per active area', 'mol/m2/s', &
                        acc_num=a_methane_acc_num)
                     CALL write_history_variable_2d (mhist_on('f_methane_surf_diff_rice'), &
                        a_methane_surf_diff_rice, file_hist, &
                        'f_methane_surf_diff_rice', itime_in_file, sumarea, filter, &
                        'rice CH4 diffusive flux contribution per active area', 'mol/m2/s', &
                        acc_num=a_methane_acc_num)
                     CALL write_history_variable_2d (mhist_on('f_methane_surf_diff_phys'), a_methane_surf_diff_phys, file_hist, &
                        'f_methane_surf_diff_phys', itime_in_file, sumarea, filter, &
                        'CH4 pure-physics diffusive surface flux (before negative-conc clip and closure residual)', 'mol/m2/s', &
                        acc_num=a_methane_acc_num)
                  CALL write_history_variable_2d (mhist_on('f_methane_balance_residual'), a_methane_balance_residual, file_hist, &
                     'f_methane_balance_residual', itime_in_file, sumarea, filter, &
                     'CH4 numerical closure flux; active-CH4-area mean; lake excluded', 'mol/m2/s', &
                     acc_num=a_methane_acc_num)
                  CALL write_history_variable_2d (mhist_on('f_methane_ch4_clip_credit'), a_methane_ch4_clip_credit, file_hist, &
                     'f_methane_ch4_clip_credit', itime_in_file, sumarea, filter, &
                     'CH4 negative-concentration clip credit; active-CH4-area mean; lake excluded', 'mol/m2/s', &
                     acc_num=a_methane_acc_num)
                  CALL write_history_variable_2d (mhist_on('f_o2_cap_loss'), a_o2_cap_loss, file_hist, &
                     'f_o2_cap_loss', itime_in_file, sumarea, filter, &
                     'O2 column sink from post-solve physical concentration cap', 'mol/m2/s', &
                     acc_num=a_methane_acc_num)
                  CALL write_history_variable_2d (mhist_on('f_o2_cap_gain'), a_o2_cap_gain, file_hist, &
                     'f_o2_cap_gain', itime_in_file, sumarea, filter, &
                     'O2 column source from post-solve nonnegative concentration floor', 'mol/m2/s', &
                     acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_ebul_tot'), a_methane_ebul_tot, file_hist, &
                  'f_methane_ebul_tot', itime_in_file, sumarea, filter, &
                  'CH4 column ebullition', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_prod_tot'), a_methane_prod_tot, file_hist, &
                  'f_methane_prod_tot', itime_in_file, sumarea, filter, &
                  'CH4 column production', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_oxid_tot'), a_methane_oxid_tot, file_hist, &
                  'f_methane_oxid_tot', itime_in_file, sumarea, filter, &
                  'CH4 column oxidation', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_prod_tot_soil'), a_methane_prod_tot_soil, file_hist, &
                  'f_methane_prod_tot_soil', itime_in_file, sumarea, filter, &
                  'non-rice CH4 production contribution per active area', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_prod_tot_rice'), a_methane_prod_tot_rice, file_hist, &
                  'f_methane_prod_tot_rice', itime_in_file, sumarea, filter, &
                  'rice CH4 production contribution per active area', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_oxid_tot_soil'), a_methane_oxid_tot_soil, file_hist, &
                  'f_methane_oxid_tot_soil', itime_in_file, sumarea, filter, &
                  'non-rice CH4 oxidation contribution per active area', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_oxid_tot_rice'), a_methane_oxid_tot_rice, file_hist, &
                  'f_methane_oxid_tot_rice', itime_in_file, sumarea, filter, &
                  'rice CH4 oxidation contribution per active area', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_co2_decomp_tot'), a_co2_decomp_tot, file_hist, &
                  'f_co2_decomp_tot', itime_in_file, sumarea, filter, &
                     'diagnostic CO2 decomp/methanogenesis before O2-stress scaling', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_co2_oxid_tot'), a_co2_oxid_tot, file_hist, &
                  'f_co2_oxid_tot', itime_in_file, sumarea, filter, &
                  'CO2 column source from CH4 oxidation', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_co2_aere_tot'), a_co2_aere_tot, file_hist, &
                  'f_co2_aere_tot', itime_in_file, sumarea, filter, &
                  'diagnosed CO2 column aerenchyma flux', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_co2_net_tot'), a_co2_net_tot, file_hist, &
                  'f_co2_net_tot', itime_in_file, sumarea, filter, &
                  'net diagnosed CO2 source from Methane module', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_totcol_methane'), a_totcol_methane, file_hist, &
                  'f_totcol_methane', itime_in_file, sumarea, filter, &
                  'total CH4 in soil column', 'mol/m2', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_grnd_methane_cond'), a_grnd_methane_cond, file_hist, &
                  'f_grnd_methane_cond', itime_in_file, sumarea, filter, &
                  'boundary-layer CH4 conductance', 'm/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_conc_o2'), a_conc_o2, file_hist, &
                  'f_conc_o2', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'O2 concentration per soil layer', 'mol/m3', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_conc_methane'), a_conc_methane, file_hist, &
                  'f_conc_methane', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 concentration per soil layer', 'mol/m3', &
                  acc_num=a_methane_acc_num)

                  ! optional microbial-pool diagnostics
                  IF (allocated(a_B_methanogen)) THEN
                     CALL write_history_variable_3d (mhist_on('f_methane_B_methanogen'), a_B_methanogen, file_hist, &
                        'f_methane_B_methanogen', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'active methanogen biomass', 'gC/m3', &
                        acc_num=a_methane_acc_num_microbe)
                     CALL write_history_variable_3d (mhist_on('f_methane_B_methanotroph'), a_B_methanotroph, file_hist, &
                        'f_methane_B_methanotroph', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'active methanotroph biomass', 'gC/m3', &
                        acc_num=a_methane_acc_num_microbe)
                     CALL write_history_variable_3d (mhist_on('f_methane_B_methanogen_dormant'), &
                        a_B_methanogen_dormant, file_hist, &
                        'f_methane_B_methanogen_dormant', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'dormant methanogen biomass; active only when use_microbial_dormancy is true', 'gC/m3', &
                        acc_num=a_methane_acc_num_microbe)
                     CALL write_history_variable_3d (&
                        mhist_on('f_methane_B_methanotroph_dormant'), a_B_methanotroph_dormant, file_hist, &
                        'f_methane_B_methanotroph_dormant', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'dormant methanotroph biomass; active only when use_microbial_dormancy is true', 'gC/m3', &
                        acc_num=a_methane_acc_num_microbe)
                     CALL write_history_variable_3d (mhist_on('f_methane_microbe_f_T_methanogen'), a_f_T_methanogen, file_hist, &
                        'f_methane_microbe_f_T_methanogen', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'methanogen microbial temperature factor', '-', &
                        acc_num=a_methane_acc_num_microbe)
                     CALL write_history_variable_3d (mhist_on('f_methane_microbe_f_S_methanogen'), a_f_S_methanogen, file_hist, &
                        'f_methane_microbe_f_S_methanogen', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'methanogen substrate factor', '-', &
                        acc_num=a_methane_acc_num_microbe)
                     CALL write_history_variable_3d (mhist_on('f_methane_microbe_f_O2_methanogen'), a_f_O2_methanogen, file_hist, &
                        'f_methane_microbe_f_O2_methanogen', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'methanogen oxygen inhibition factor', '-', &
                        acc_num=a_methane_acc_num_microbe)
                     CALL write_history_variable_3d (mhist_on('f_methane_microbe_f_T_methanotroph'), &
                        a_f_T_methanotroph, file_hist, &
                        'f_methane_microbe_f_T_methanotroph', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'methanotroph microbial temperature factor', '-', &
                        acc_num=a_methane_acc_num_microbe)
                     CALL write_history_variable_3d (&
                        mhist_on('f_methane_methanogen_growth_rate'), a_methanogen_growth_rate, file_hist, &
                        'f_methane_methanogen_growth_rate', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'methanogen active-pool growth rate', 'day-1', &
                        acc_num=a_methane_acc_num_microbe)
                     CALL write_history_variable_3d (&
                        mhist_on('f_methane_methanotroph_growth_rate'), a_methanotroph_growth_rate, file_hist, &
                        'f_methane_methanotroph_growth_rate', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'methanotroph active-pool growth rate', 'day-1', &
                        acc_num=a_methane_acc_num_microbe)
                     CALL write_history_variable_3d (&
                        mhist_on('f_methane_microbial_prod_potential'), a_microbial_prod_potential, file_hist, &
                        'f_methane_microbial_prod_potential', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'carbon-capped microbial CH4 production potential', 'mol/m3/s', &
                        acc_num=a_methane_acc_num_microbe)
                     CALL write_history_variable_3d (&
                        mhist_on('f_methane_microbial_oxid_potential'), a_microbial_oxid_potential, file_hist, &
                        'f_methane_microbial_oxid_potential', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'CH4/O2-capped microbial CH4 oxidation potential', 'mol/m3/s', &
                        acc_num=a_methane_acc_num_microbe)
                  ENDIF

               ! unsat / sat
               CALL write_history_variable_2d (mhist_on('f_net_methane_unsat'), a_net_methane_unsat, file_hist, &
                  'f_net_methane_unsat', itime_in_file, sumarea, filter, &
                  'CH4 oxidation minus production (unsaturated subcolumn)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_unsat)
               CALL write_history_variable_2d (mhist_on('f_net_methane_sat'), a_net_methane_sat, file_hist, &
                  'f_net_methane_sat', itime_in_file, sumarea, filter, &
                  'CH4 oxidation minus production (saturated subcolumn)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_sat)
               CALL write_history_variable_3d (mhist_on('f_methane_prod_depth_unsat'), a_methane_prod_depth_unsat, file_hist, &
                  'f_methane_prod_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 production per layer (unsaturated)', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_unsat)
               CALL write_history_variable_3d (mhist_on('f_methane_prod_depth_sat'), a_methane_prod_depth_sat, file_hist, &
                  'f_methane_prod_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 production per layer (saturated)', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_sat)
               CALL write_history_variable_3d (mhist_on('f_methane_oxid_depth_unsat'), a_methane_oxid_depth_unsat, file_hist, &
                  'f_methane_oxid_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 oxidation per layer (unsaturated)', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_unsat)
                  CALL write_history_variable_3d (mhist_on('f_methane_oxid_depth_sat'), a_methane_oxid_depth_sat, file_hist, &
                     'f_methane_oxid_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'CH4 oxidation per layer (saturated)', 'mol/m3/s', &
                     acc_num=a_methane_acc_num_sat)
                  CALL write_history_variable_3d (mhist_on('f_o2_decomp_depth_unsat'), a_o2_decomp_depth_unsat, file_hist, &
                     'f_o2_decomp_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'O2 consumption from decomposition/root respiration per layer (unsaturated)', 'mol/m3/s', &
                     acc_num=a_methane_acc_num_unsat)
                     CALL write_history_variable_3d (mhist_on('f_o2_decomp_depth_sat'), a_o2_decomp_depth_sat, file_hist, &
                        'f_o2_decomp_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'O2 consumption from decomposition/root respiration per layer (saturated)', 'mol/m3/s', &
                        acc_num=a_methane_acc_num_sat)
                     CALL write_history_variable_3d (mhist_on('f_o2_oxid_depth_unsat'), a_o2_oxid_depth_unsat, file_hist, &
                        'f_o2_oxid_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'O2 consumption from CH4 oxidation per layer (unsaturated)', 'mol/m3/s', &
                        acc_num=a_methane_acc_num_unsat)
                     CALL write_history_variable_3d (mhist_on('f_o2_oxid_depth_sat'), a_o2_oxid_depth_sat, file_hist, &
                        'f_o2_oxid_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'O2 consumption from CH4 oxidation per layer (saturated)', 'mol/m3/s', &
                        acc_num=a_methane_acc_num_sat)
                     CALL write_history_variable_3d (mhist_on('f_o2_aere_depth_unsat'), a_o2_aere_depth_unsat, file_hist, &
                        'f_o2_aere_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'O2 aerenchyma supply per layer (unsaturated)', 'mol/m3/s', &
                        acc_num=a_methane_acc_num_unsat)
                     CALL write_history_variable_3d (mhist_on('f_o2_aere_depth_sat'), a_o2_aere_depth_sat, file_hist, &
                        'f_o2_aere_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'O2 aerenchyma supply per layer (saturated)', 'mol/m3/s', &
                        acc_num=a_methane_acc_num_sat)
                     CALL write_history_variable_3d (mhist_on('f_o2stress_unsat'), a_o2stress_unsat, file_hist, &
                        'f_o2stress_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'O2 availability stress factor (unsaturated)', '-', &
                        acc_num=a_methane_acc_num_unsat)
                     CALL write_history_variable_3d (mhist_on('f_o2stress_sat'), a_o2stress_sat, file_hist, &
                        'f_o2stress_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'O2 availability stress factor (saturated)', '-', &
                        acc_num=a_methane_acc_num_sat)
                     CALL write_history_variable_3d (mhist_on('f_methane_stress_unsat'), a_methane_stress_unsat, file_hist, &
                        'f_methane_stress_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'CH4 availability stress factor (unsaturated)', '-', &
                        acc_num=a_methane_acc_num_unsat)
                     CALL write_history_variable_3d (mhist_on('f_methane_stress_sat'), a_methane_stress_sat, file_hist, &
                        'f_methane_stress_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'CH4 availability stress factor (saturated)', '-', &
                        acc_num=a_methane_acc_num_sat)
                  CALL write_history_variable_3d (mhist_on('f_methane_aere_depth_unsat'), a_methane_aere_depth_unsat, file_hist, &
                     'f_methane_aere_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'CH4 aerenchyma loss per layer (unsaturated)', 'mol/m3/s', &
                     acc_num=a_methane_acc_num_unsat)
                  CALL write_history_variable_3d (mhist_on('f_methane_aere_depth_sat'), a_methane_aere_depth_sat, file_hist, &
                     'f_methane_aere_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'CH4 aerenchyma loss per layer (saturated)', 'mol/m3/s', &
                     acc_num=a_methane_acc_num_sat)
                  CALL write_history_variable_3d (mhist_on('f_methane_tran_depth_unsat'), a_methane_tran_depth_unsat, file_hist, &
                     'f_methane_tran_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'CH4 transpiration loss per layer (unsaturated)', 'mol/m3/s', &
                     acc_num=a_methane_acc_num_unsat)
                  CALL write_history_variable_3d (mhist_on('f_methane_tran_depth_sat'), a_methane_tran_depth_sat, file_hist, &
                     'f_methane_tran_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'CH4 transpiration loss per layer (saturated)', 'mol/m3/s', &
                     acc_num=a_methane_acc_num_sat)
                  CALL write_history_variable_3d (mhist_on('f_methane_ebul_depth_unsat'), a_methane_ebul_depth_unsat, file_hist, &
                     'f_methane_ebul_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'CH4 ebullition loss per layer (unsaturated)', 'mol/m3/s', &
                     acc_num=a_methane_acc_num_unsat)
                  CALL write_history_variable_3d (mhist_on('f_methane_ebul_depth_sat'), a_methane_ebul_depth_sat, file_hist, &
                     'f_methane_ebul_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'CH4 ebullition loss per layer (saturated)', 'mol/m3/s', &
                     acc_num=a_methane_acc_num_sat)
                  CALL write_history_variable_3d (mhist_on('f_co2_decomp_depth_unsat'), a_co2_decomp_depth_unsat, file_hist, &
                     'f_co2_decomp_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'diagnostic CO2 decomp/methanogenesis before O2-stress scaling (unsaturated)', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_unsat)
               CALL write_history_variable_3d (mhist_on('f_co2_decomp_depth_sat'), a_co2_decomp_depth_sat, file_hist, &
                  'f_co2_decomp_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'diagnostic CO2 decomp/methanogenesis before O2-stress scaling (saturated)', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_sat)
               CALL write_history_variable_3d (mhist_on('f_co2_oxid_depth_unsat'), a_co2_oxid_depth_unsat, file_hist, &
                  'f_co2_oxid_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CO2 from CH4 oxidation per layer (unsaturated)', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_unsat)
                  CALL write_history_variable_3d (mhist_on('f_co2_oxid_depth_sat'), a_co2_oxid_depth_sat, file_hist, &
                     'f_co2_oxid_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'CO2 from CH4 oxidation per layer (saturated)', 'mol/m3/s', &
                     acc_num=a_methane_acc_num_sat)
                  CALL write_history_variable_3d (mhist_on('f_co2_aere_depth_unsat'), a_co2_aere_depth_unsat, file_hist, &
                     'f_co2_aere_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'CO2 aerenchyma diagnostic per layer (unsaturated)', 'mol/m3/s', &
                     acc_num=a_methane_acc_num_unsat)
                  CALL write_history_variable_3d (mhist_on('f_co2_aere_depth_sat'), a_co2_aere_depth_sat, file_hist, &
                     'f_co2_aere_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'CO2 aerenchyma diagnostic per layer (saturated)', 'mol/m3/s', &
                     acc_num=a_methane_acc_num_sat)
               CALL write_history_variable_2d (mhist_on('f_co2_decomp_tot_unsat'), a_co2_decomp_tot_unsat, file_hist, &
                  'f_co2_decomp_tot_unsat', itime_in_file, sumarea, filter, &
                     'diagnostic CO2 decomp/methanogenesis before O2-stress scaling (unsaturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_unsat)
               CALL write_history_variable_2d (mhist_on('f_co2_decomp_tot_sat'), a_co2_decomp_tot_sat, file_hist, &
                  'f_co2_decomp_tot_sat', itime_in_file, sumarea, filter, &
                     'diagnostic CO2 decomp/methanogenesis before O2-stress scaling (saturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_sat)
               CALL write_history_variable_2d (mhist_on('f_co2_oxid_tot_unsat'), a_co2_oxid_tot_unsat, file_hist, &
                  'f_co2_oxid_tot_unsat', itime_in_file, sumarea, filter, &
                  'CO2 oxidation column source (unsaturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_unsat)
               CALL write_history_variable_2d (mhist_on('f_co2_oxid_tot_sat'), a_co2_oxid_tot_sat, file_hist, &
                  'f_co2_oxid_tot_sat', itime_in_file, sumarea, filter, &
                  'CO2 oxidation column source (saturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_sat)
               CALL write_history_variable_2d (mhist_on('f_co2_net_tot_unsat'), a_co2_net_tot_unsat, file_hist, &
                  'f_co2_net_tot_unsat', itime_in_file, sumarea, filter, &
                  'net diagnosed CO2 source (unsaturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_unsat)
               CALL write_history_variable_2d (mhist_on('f_co2_net_tot_sat'), a_co2_net_tot_sat, file_hist, &
                  'f_co2_net_tot_sat', itime_in_file, sumarea, filter, &
                  'net diagnosed CO2 source (saturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_sat)
               CALL write_history_variable_2d (&
                  mhist_on('f_methane_surf_flux_tot_unsat'), a_methane_surf_flux_tot_unsat, file_hist, &
                  'f_methane_surf_flux_tot_unsat', itime_in_file, sumarea, filter, &
                  'CH4 surface flux (unsaturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_unsat)
                  CALL write_history_variable_2d (mhist_on('f_methane_surf_flux_tot_sat'), a_methane_surf_flux_tot_sat, file_hist, &
                     'f_methane_surf_flux_tot_sat', itime_in_file, sumarea, filter, &
                     'CH4 surface flux (saturated)', 'mol/m2/s', &
                     acc_num=a_methane_acc_num_sat)
                  CALL write_history_variable_2d (mhist_on('f_methane_surf_aere_unsat'), a_methane_surf_aere_unsat, file_hist, &
                     'f_methane_surf_aere_unsat', itime_in_file, sumarea, filter, &
                     'CH4 column aerenchyma flux (unsaturated)', 'mol/m2/s', &
                     acc_num=a_methane_acc_num_unsat)
                  CALL write_history_variable_2d (mhist_on('f_methane_surf_aere_sat'), a_methane_surf_aere_sat, file_hist, &
                     'f_methane_surf_aere_sat', itime_in_file, sumarea, filter, &
                     'CH4 column aerenchyma flux (saturated)', 'mol/m2/s', &
                     acc_num=a_methane_acc_num_sat)
                  CALL write_history_variable_2d (mhist_on('f_methane_surf_ebul_unsat'), a_methane_surf_ebul_unsat, file_hist, &
                     'f_methane_surf_ebul_unsat', itime_in_file, sumarea, filter, &
                     'CH4 surface ebullition flux (unsaturated)', 'mol/m2/s', &
                     acc_num=a_methane_acc_num_unsat)
                  CALL write_history_variable_2d (mhist_on('f_methane_surf_ebul_sat'), a_methane_surf_ebul_sat, file_hist, &
                     'f_methane_surf_ebul_sat', itime_in_file, sumarea, filter, &
                     'CH4 surface ebullition flux (saturated)', 'mol/m2/s', &
                     acc_num=a_methane_acc_num_sat)
                  CALL write_history_variable_2d (mhist_on('f_methane_surf_diff_unsat'), a_methane_surf_diff_unsat, file_hist, &
                     'f_methane_surf_diff_unsat', itime_in_file, sumarea, filter, &
                     'CH4 diffusive surface flux (unsaturated)', 'mol/m2/s', &
                     acc_num=a_methane_acc_num_unsat)
                     CALL write_history_variable_2d (mhist_on('f_methane_surf_diff_sat'), a_methane_surf_diff_sat, file_hist, &
                        'f_methane_surf_diff_sat', itime_in_file, sumarea, filter, &
                        'CH4 diffusive surface flux (saturated)', 'mol/m2/s', &
                        acc_num=a_methane_acc_num_sat)
                     CALL write_history_variable_2d (mhist_on('f_methane_ebul_tot_unsat'), a_methane_ebul_tot_unsat, file_hist, &
                        'f_methane_ebul_tot_unsat', itime_in_file, sumarea, filter, &
                        'CH4 column ebullition (unsaturated)', 'mol/m2/s', &
                        acc_num=a_methane_acc_num_unsat)
                     CALL write_history_variable_2d (mhist_on('f_methane_ebul_tot_sat'), a_methane_ebul_tot_sat, file_hist, &
                        'f_methane_ebul_tot_sat', itime_in_file, sumarea, filter, &
                        'CH4 column ebullition (saturated)', 'mol/m2/s', &
                        acc_num=a_methane_acc_num_sat)
                     CALL write_history_variable_2d (mhist_on('f_methane_prod_tot_unsat'), a_methane_prod_tot_unsat, file_hist, &
                        'f_methane_prod_tot_unsat', itime_in_file, sumarea, filter, &
                        'CH4 column production (unsaturated)', 'mol/m2/s', &
                        acc_num=a_methane_acc_num_unsat)
                     CALL write_history_variable_2d (mhist_on('f_methane_prod_tot_sat'), a_methane_prod_tot_sat, file_hist, &
                        'f_methane_prod_tot_sat', itime_in_file, sumarea, filter, &
                        'CH4 column production (saturated)', 'mol/m2/s', &
                        acc_num=a_methane_acc_num_sat)
                     CALL write_history_variable_2d (mhist_on('f_methane_oxid_tot_unsat'), a_methane_oxid_tot_unsat, file_hist, &
                        'f_methane_oxid_tot_unsat', itime_in_file, sumarea, filter, &
                        'CH4 column oxidation (unsaturated)', 'mol/m2/s', &
                        acc_num=a_methane_acc_num_unsat)
                     CALL write_history_variable_2d (mhist_on('f_methane_oxid_tot_sat'), a_methane_oxid_tot_sat, file_hist, &
                        'f_methane_oxid_tot_sat', itime_in_file, sumarea, filter, &
                        'CH4 column oxidation (saturated)', 'mol/m2/s', &
                        acc_num=a_methane_acc_num_sat)
                     CALL write_history_variable_2d (mhist_on('f_totcol_methane_unsat'), a_totcol_methane_unsat, file_hist, &
                        'f_totcol_methane_unsat', itime_in_file, sumarea, filter, &
                        'total CH4 in unsaturated column', 'mol/m2', &
                        acc_num=a_methane_acc_num_unsat)
                     CALL write_history_variable_2d (mhist_on('f_totcol_methane_sat'), a_totcol_methane_sat, file_hist, &
                        'f_totcol_methane_sat', itime_in_file, sumarea, filter, &
                        'total CH4 in saturated column', 'mol/m2', &
                        acc_num=a_methane_acc_num_sat)
                     CALL write_history_variable_2d (mhist_on('f_grnd_methane_cond_unsat'), a_grnd_methane_cond_unsat, file_hist, &
                        'f_grnd_methane_cond_unsat', itime_in_file, sumarea, filter, &
                        'CH4 boundary-layer conductance (unsaturated)', 'm/s', &
                        acc_num=a_methane_acc_num_unsat)
                     CALL write_history_variable_2d (mhist_on('f_grnd_methane_cond_sat'), a_grnd_methane_cond_sat, file_hist, &
                        'f_grnd_methane_cond_sat', itime_in_file, sumarea, filter, &
                        'CH4 boundary-layer conductance (saturated)', 'm/s', &
                        acc_num=a_methane_acc_num_sat)
                  CALL write_history_variable_3d (mhist_on('f_conc_o2_unsat'), a_conc_o2_unsat, file_hist, &
                     'f_conc_o2_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'O2 concentration per layer (unsaturated)', 'mol/m3', &
                     acc_num=a_methane_acc_num_unsat)
                  CALL write_history_variable_3d (mhist_on('f_conc_o2_sat'), a_conc_o2_sat, file_hist, &
                     'f_conc_o2_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'O2 concentration per layer (saturated)', 'mol/m3', &
                     acc_num=a_methane_acc_num_sat)
                  CALL write_history_variable_3d (mhist_on('f_conc_methane_unsat'), a_conc_methane_unsat, file_hist, &
                     'f_conc_methane_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 concentration per layer (unsaturated)', 'mol/m3', &
                  acc_num=a_methane_acc_num_unsat)
	               CALL write_history_variable_3d (mhist_on('f_conc_methane_sat'), a_conc_methane_sat, file_hist, &
	                  'f_conc_methane_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                  'CH4 concentration per layer (saturated)', 'mol/m3', &
	                  acc_num=a_methane_acc_num_sat)

	               ! Land-area denominator diagnostics for global integrations.
	               ! active_total_without_lake = wetland + soil + rice, excluding lake.
	               ! global_total_with_lake    = wetland + soil + rice + lake.
		               IF ((need_active_without_lake .or. need_global_with_lake .or. &
		                       need_global_phys_with_lake .or. need_global_balance_with_lake .or. &
		                       need_global_clip_credit_with_lake .or. &
		                       need_land_flux_split) .and. &
                       HistForm == 'Gridded') THEN
	                  CALL mp2g_hist%get_sumarea (sumarea, filter_all_land)
	               ENDIF
	               IF (need_land_flux_total) THEN
	                  CALL write_history_variable_2d (.true., &
	                     hist_ch4_land_total_surf, file_hist, &
	                     'f_methane_surf_flux_tot', itime_in_file, &
	                     sumarea, filter_all_land, &
	                     'CH4 total surface flux including wetland/soil/rice/lake; land-area mean', &
	                     'mol/m2/s', acc_num=hist_ch4_acc_one)
	               ENDIF
	               IF (need_land_flux_wetland) THEN
	                  CALL write_history_variable_2d (.true., &
	                     hist_ch4_land_flux_wetland, file_hist, &
	                     'f_methane_surf_flux_wetland', itime_in_file, &
	                     sumarea, filter_all_land, &
	                     'wetland contribution to CH4 surface flux; land-area mean; multiply by landarea for global total', &
	                     'mol/m2/s', acc_num=hist_ch4_acc_one)
	               ENDIF
	               IF (need_land_flux_soil) THEN
	                  CALL write_history_variable_2d (.true., &
	                     hist_ch4_land_flux_soil, file_hist, &
	                     'f_methane_surf_flux_soil', itime_in_file, &
	                     sumarea, filter_all_land, &
	                     'non-rice soil contribution to CH4 surface flux; land-area mean; multiply by landarea for global total', &
	                     'mol/m2/s', acc_num=hist_ch4_acc_one)
	               ENDIF
	               IF (need_land_flux_lake) THEN
	                  CALL write_history_variable_2d (.true., &
	                     hist_ch4_land_flux_lake, file_hist, &
	                     'f_methane_surf_flux_lake', itime_in_file, &
	                     sumarea, filter_all_land, &
	                     'lake contribution to CH4 surface flux; land-area mean; multiply by landarea for global total', &
	                     'mol/m2/s', acc_num=hist_ch4_acc_one)
	               ENDIF
	               IF (need_land_flux_rice) THEN
	                  CALL write_history_variable_2d (.true., &
	                     hist_ch4_land_flux_rice, file_hist, &
	                     'f_methane_surf_flux_rice', itime_in_file, &
	                     sumarea, filter_all_land, &
	                     'rice-paddy contribution to CH4 surface flux; land-area mean; multiply by landarea for global total', &
	                     'mol/m2/s', acc_num=hist_ch4_acc_one)
	               ENDIF
	               IF (need_active_without_lake) THEN
	                  CALL write_history_variable_2d (.true., &
	                     hist_ch4_active_without_lake, file_hist, &
	                     'f_methane_surf_flux_active_total_without_lake', itime_in_file, &
	                     sumarea, filter_active_without_lake, &
	                     'active CH4 total surface flux excluding lake; land-area mean contribution', 'mol/m2/s', &
	                     acc_num=a_methane_acc_num)
	               ENDIF
		               IF (need_global_with_lake) THEN
	                  CALL write_history_variable_2d (.true., &
	                     hist_ch4_global_with_lake, file_hist, &
	                     'f_methane_surf_flux_global_total_with_lake', itime_in_file, &
	                     sumarea, filter_all_land, &
	                     'global CH4 total surface flux including lake; land-area mean contribution', 'mol/m2/s', &
	                     acc_num=hist_ch4_acc_one)
	               ENDIF
	               IF (need_global_phys_with_lake) THEN
	                  CALL write_history_variable_2d (.true., &
	                     hist_ch4_global_phys_with_lake, file_hist, &
	                     'f_methane_surf_flux_global_phys_with_lake', itime_in_file, &
	                     sumarea, filter_all_land, &
	                     'physical CH4 surface flux including lake; land-area mean contribution', 'mol/m2/s', &
	                     acc_num=hist_ch4_acc_one)
	               ENDIF
	               IF (need_global_balance_with_lake) THEN
	                  CALL write_history_variable_2d (.true., &
	                     hist_ch4_global_balance_with_lake, file_hist, &
	                     'f_methane_balance_residual_global_with_lake', itime_in_file, &
	                     sumarea, filter_all_land, &
	                     'CH4 column balance residual including lake; land-area mean contribution', 'mol/m2/s', &
	                     acc_num=hist_ch4_acc_one)
	               ENDIF
	               IF (need_global_clip_credit_with_lake) THEN
	                  CALL write_history_variable_2d (.true., &
	                     hist_ch4_global_clip_credit_with_lake, file_hist, &
	                     'f_methane_ch4_clip_credit_global_with_lake', itime_in_file, &
	                     sumarea, filter_all_land, &
	                     'CH4 nonnegative-storage correction including lake; land-area mean contribution', 'mol/m2/s', &
	                     acc_num=hist_ch4_acc_one)
	               ENDIF


	               ! lake diagnostics (CTSM lake CH4 path) — swap to lake-only filter
               IF ((p_is_worker) .and. (numpatch > 0)) THEN
                  filter = (patchtype == 4) .and. patchmask
                  IF (forcing_has_missing_value) filter = filter .and. forcmask_pch
               ENDIF
               IF (need_lake_history .and. HistForm == 'Gridded') THEN
                  CALL mp2g_hist%get_sumarea (sumarea, filter)
               ENDIF

               CALL write_history_variable_3d (mhist_on('f_methane_prod_depth_lake'), a_methane_prod_depth_lake, file_hist, &
                  'f_methane_prod_depth_lake', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'lake CH4 production per sediment layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_3d (mhist_on('f_methane_oxid_depth_lake'), a_methane_oxid_depth_lake, file_hist, &
                  'f_methane_oxid_depth_lake', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'lake CH4 oxidation per sediment layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_3d (mhist_on('f_co2_decomp_depth_lake'), a_co2_decomp_depth_lake, file_hist, &
                  'f_co2_decomp_depth_lake', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'lake CO2 decomp/methanogenesis per sediment layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_3d (mhist_on('f_co2_oxid_depth_lake'), a_co2_oxid_depth_lake, file_hist, &
                  'f_co2_oxid_depth_lake', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'lake CO2 from CH4 oxidation per sediment layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_3d (mhist_on('f_methane_ebul_depth_lake'), a_methane_ebul_depth_lake, file_hist, &
                  'f_methane_ebul_depth_lake', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'lake CH4 ebullition loss per sediment layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_methane_surf_ebul_lake'), a_methane_surf_ebul_lake, file_hist, &
                  'f_methane_surf_ebul_lake', itime_in_file, sumarea, filter, &
                  'lake CH4 surface ebullition flux', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_methane_surf_diff_lake'), a_methane_surf_diff_lake, file_hist, &
                  'f_methane_surf_diff_lake', itime_in_file, sumarea, filter, &
                  'lake CH4 surface diffusion flux', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
                  CALL write_history_variable_2d (mhist_on('f_methane_surf_flux_tot_lake'), &
                     a_methane_surf_flux_tot_lake, file_hist, &
                     'f_methane_surf_flux_tot_lake', itime_in_file, sumarea, filter, &
                     'lake CH4 total surface flux (ebullition + diffusion)', 'mol/m2/s', &
                     acc_num=a_methane_acc_num_lake)
		                  CALL write_history_variable_2d (.false., a_methane_surf_flux_lake, file_hist, &
		                     'f_methane_surf_flux_lake', itime_in_file, sumarea, filter, &
		                     'lake CH4 total surface flux; lake-area intensive, multiply by area_lake for global total', 'mol/m2/s', &
		                     acc_num=a_methane_acc_num_lake)
		                  IF (need_lake_intensive) THEN
		                     CALL write_history_variable_2d (.true., &
		                        hist_ch4_lake_flux_acc, file_hist, &
		                        'f_methane_surf_flux_lake_intensive', itime_in_file, sumarea, filter, &
		                        'lake CH4 total surface flux; explicit lake-area intensive alias', 'mol/m2/s', &
		                        acc_num=a_methane_acc_num_lake)
		                  ENDIF
                  CALL write_history_variable_2d (mhist_on('f_methane_prod_tot_lake'), a_methane_prod_tot_lake, file_hist, &
                  'f_methane_prod_tot_lake', itime_in_file, sumarea, filter, &
                  'lake CH4 column production', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_methane_oxid_tot_lake'), a_methane_oxid_tot_lake, file_hist, &
                  'f_methane_oxid_tot_lake', itime_in_file, sumarea, filter, &
                  'lake CH4 column oxidation', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_methane_ebul_tot_lake'), a_methane_ebul_tot_lake, file_hist, &
                  'f_methane_ebul_tot_lake', itime_in_file, sumarea, filter, &
                  'lake CH4 column ebullition', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_co2_decomp_tot_lake'), a_co2_decomp_tot_lake, file_hist, &
                  'f_co2_decomp_tot_lake', itime_in_file, sumarea, filter, &
                  'lake CO2 decomp/methanogenesis column source', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_co2_oxid_tot_lake'), a_co2_oxid_tot_lake, file_hist, &
                  'f_co2_oxid_tot_lake', itime_in_file, sumarea, filter, &
                  'lake CO2 oxidation column source', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_co2_net_tot_lake'), a_co2_net_tot_lake, file_hist, &
                  'f_co2_net_tot_lake', itime_in_file, sumarea, filter, &
                  'lake net diagnosed CO2 source', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_totcol_methane_lake'), a_totcol_methane_lake, file_hist, &
                  'f_totcol_methane_lake', itime_in_file, sumarea, filter, &
                  'total CH4 in lake sediment column', 'mol/m2', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_grnd_methane_cond_lake'), a_grnd_methane_cond_lake, file_hist, &
                  'f_grnd_methane_cond_lake', itime_in_file, sumarea, filter, &
                  'lake CH4 boundary-layer conductance', 'm/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_3d (mhist_on('f_conc_o2_lake'), a_conc_o2_lake, file_hist, &
                  'f_conc_o2_lake', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'lake O2 concentration per sediment layer', 'mol/m3', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_3d (mhist_on('f_conc_methane_lake'), a_conc_methane_lake, file_hist, &
                  'f_conc_methane_lake', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'lake CH4 concentration per sediment layer', 'mol/m3', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_lake_water_ch4_stock'), a_lake_water_ch4_stock, file_hist, &
                  'f_lake_water_ch4_stock', itime_in_file, sumarea, filter, &
                  'dissolved CH4 stock in the well-mixed lake water column', 'mol/m2', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_lake_water_o2_stock'), a_lake_water_o2_stock, file_hist, &
                  'f_lake_water_o2_stock', itime_in_file, sumarea, filter, &
                  'dissolved O2 stock in the well-mixed lake water column', 'mol/m2', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_lake_water_ch4_oxid'), a_lake_water_ch4_oxid, file_hist, &
                  'f_lake_water_ch4_oxid', itime_in_file, sumarea, filter, &
                  'lake-water CH4 oxidation sink', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_lake_sed_ch4_flux'), a_lake_sed_ch4_flux, file_hist, &
                  'f_lake_sed_ch4_flux', itime_in_file, sumarea, filter, &
                  'diffusive CH4 flux from sediment to lake water', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_lake_sed_o2_flux'), a_lake_sed_o2_flux, file_hist, &
                  'f_lake_sed_o2_flux', itime_in_file, sumarea, filter, &
                  'diffusive O2 flux from sediment to lake water (positive upward)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_lake_air_o2_flux'), a_lake_air_o2_flux, file_hist, &
                  'f_lake_air_o2_flux', itime_in_file, sumarea, filter, &
                  'lake-water to atmosphere O2 exchange flux (positive upward)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)

               ! extras + microbe diagnostics — swap back to active CH4 mask
               IF ((p_is_worker) .and. (numpatch > 0)) THEN
                  ! Use BgcLink helper so the history filter matches the
                  ! driver gate + rice paddy override exactly.  Same call
                  ! site pattern as MOD_MT_Methane_AccFlux.F90.
                  CALL methane_patch_active_mask (numpatch, DEF_METHANE%only_wetland, &
                                                  DEF_METHANE%enable_rice_paddy, &
                                                  patchtype, filter)
                  filter = filter .and. patchmask
                  IF (forcing_has_missing_value) filter = filter .and. forcmask_pch
               ENDIF
               IF (need_extra_history .and. HistForm == 'Gridded') THEN
                  CALL mp2g_hist%get_sumarea (sumarea, filter)
               ENDIF

               CALL write_history_variable_2d (mhist_on('f_forc_pmethanem'), a_forc_pmethanem, file_hist, &
                  'f_forc_pmethanem', itime_in_file, sumarea, filter, &
                  'atmospheric CH4 partial pressure forcing', 'Pa', &
                  acc_num=a_methane_acc_num_extra)
               CALL write_history_variable_3d (mhist_on('f_layer_sat_lag'), a_layer_sat_lag, file_hist, &
                  'f_layer_sat_lag', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'time-lagged saturation flag per layer', '-', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_annavg_agnpp'), a_annavg_agnpp, file_hist, &
                  'f_annavg_agnpp', itime_in_file, sumarea, filter, &
                  'annual avg above-ground NPP', 'gC/m2/s', &
                  acc_num=a_methane_acc_num_extra)
               CALL write_history_variable_2d (mhist_on('f_annavg_bgnpp'), a_annavg_bgnpp, file_hist, &
                  'f_annavg_bgnpp', itime_in_file, sumarea, filter, &
                  'annual avg below-ground NPP', 'gC/m2/s', &
                  acc_num=a_methane_acc_num_extra)
               CALL write_history_variable_2d (mhist_on('f_annavg_somhr'), a_annavg_somhr, file_hist, &
                  'f_annavg_somhr', itime_in_file, sumarea, filter, &
                  'annual avg SOM heterotrophic respiration', 'gC/m2/s', &
                  acc_num=a_methane_acc_num_extra)
               CALL write_history_variable_2d (mhist_on('f_annavg_finrw'), a_annavg_finrw, file_hist, &
                  'f_annavg_finrw', itime_in_file, sumarea, filter, &
                  'respiration-weighted annual avg inundated fraction', '-', &
                  acc_num=a_annavg_finrw_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_dfsat_tot'), a_methane_dfsat_tot, file_hist, &
                  'f_methane_dfsat_tot', itime_in_file, sumarea, filter, &
                  'CH4 flux from decreasing inundated fraction', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_extra)
                  CALL write_history_variable_2d (mhist_on('f_f_h2osfc'), a_f_h2osfc, file_hist, &
                     'f_f_h2osfc', itime_in_file, sumarea, filter, &
                     'fraction of surface water (CLM5 microtopography)', '-', &
                     acc_num=a_methane_acc_num_extra)
                  CALL write_history_variable_2d (mhist_on('f_methane_finundated'), a_methane_finundated, file_hist, &
                     'f_methane_finundated', itime_in_file, sumarea, filter, &
                     'actual inundated fraction used by CH4 physics', '-', &
                     acc_num=a_methane_acc_num)
                  ! f_methane_soil_finundated moved to soil-only filter block below (avoid wetland zero dilution)
                  CALL write_history_variable_2d (mhist_on('f_inund_flood_patch'), a_f_inund_flood_patch, file_hist, &
                     'f_inund_flood_patch', itime_in_file, sumarea, filter, &
                     'routing-published flood inundation fraction on methane-active patches', '-', &
                     acc_num=a_methane_acc_num)
                  CALL write_history_variable_2d (mhist_on('f_inund_flood_depth_patch'), a_f_inund_flood_depth_patch, file_hist, &
                     'f_inund_flood_depth_patch', itime_in_file, sumarea, filter, &
                     'routing-published floodplain water depth on methane-active patches', 'm', &
                     acc_num=a_methane_acc_num)
                  CALL write_history_variable_2d (mhist_on('f_wetland_frac_patch'), a_wetland_frac_per_patch, file_hist, &
                     'f_wetland_frac_patch', itime_in_file, sumarea, filter, &
                     'static wetland fraction of parent cell for methane patch', '-', &
                     acc_num=a_methane_acc_num)
                  ! Canonical wetland/soil split fluxes are written above as
                  ! all-land contributions.  Do not write those names again
                  ! with the active-CH4 denominator.
	                  CALL write_history_variable_2d (.false., a_methane_surf_flux_rice, file_hist, &
	                     'f_methane_surf_flux_rice', itime_in_file, sumarea, filter, &
	                     'rice-paddy contribution to total CH4 surface flux (active-mean)', 'mol/m2/s', &
	                     acc_num=a_methane_acc_num)
	                  ! Rice is a CFT fraction inside soil patches.  This diagnostic
	                  ! divides the accumulated rice flux contribution by the rice
	                  ! area fraction accumulated over the same samples.  The ratio
	                  ! therefore remains paddy-area intensive when that fraction
	                  ! changes during a history interval.
	                  IF (need_rice_intensive) THEN
	                     CALL write_history_mt_ratio_2d (.true., &
	                        hist_ch4_rice_flux_mean, hist_ch4_rice_area_frac, file_hist, &
	                        'f_methane_surf_flux_rice_intensive', itime_in_file, filter, &
	                        'rice-paddy CH4 surface flux; paddy-area intensive', 'mol/m2/s')
	                  ENDIF

	                  IF ((p_is_worker) .and. (numpatch > 0)) THEN
                     filter = (patchtype == 0) .and. patchmask
                     IF (forcing_has_missing_value) filter = filter .and. forcmask_pch
                  ENDIF
                  IF (need_soil_history .and. HistForm == 'Gridded') THEN
                     CALL mp2g_hist%get_sumarea (sumarea, filter)
                  ENDIF
                  CALL write_history_variable_2d (mhist_on('f_methane_soil_zwt'), a_methane_soil_zwt, file_hist, &
                     'f_methane_soil_zwt', itime_in_file, sumarea, filter, &
                     'water-table depth on methane-active soil/rice patches', 'm', &
                     acc_num=a_methane_acc_num)
                  CALL write_history_variable_2d (mhist_on('f_methane_soil_finundated'), a_methane_soil_finundated, file_hist, &
                     'f_methane_soil_finundated', itime_in_file, sumarea, filter, &
                     'CH4 inundated fraction on methane-active soil/rice patches (soil-only mean, no wetland-zero dilution)', '-', &
                     acc_num=a_methane_acc_num)
               ENDIF

	            IF (methane_history_match_count == 0) THEN
	               CALL CoLM_stop (' ***** ERROR: ch4_history_vars selects no known methane history variable: '// &
	                  trim(DEF_METHANE%ch4_history_vars))
	            ENDIF

   END SUBROUTINE methane_reactive_history

   logical FUNCTION mhist_on (varname)

      USE MOD_MT_Methane_Const, only: methane_history_enabled
      IMPLICIT NONE
      character(len=*), intent(in) :: varname

      mhist_on = methane_history_enabled(varname)
      IF (mhist_on) methane_history_match_count = methane_history_match_count + 1

   END FUNCTION mhist_on

   logical FUNCTION methane_history_selector_is_off (raw_selector)

      USE MOD_MT_Defs, only: material_lower
      IMPLICIT NONE
      character(len=*), intent(in) :: raw_selector
      character(len=len(raw_selector)) :: selector

      selector = material_lower(adjustl(trim(raw_selector)))
      SELECT CASE (trim(selector))
      CASE ('none', 'off', 'false', '.false.')
         methane_history_selector_is_off = .true.
      CASE DEFAULT
         methane_history_selector_is_off = .false.
      END SELECT

   END FUNCTION methane_history_selector_is_off

   logical FUNCTION methane_history_all_land_only (raw_selector)

      USE MOD_MT_Defs, only: material_lower
      IMPLICIT NONE
      character(len=*), intent(in) :: raw_selector
      character(len=len(raw_selector)) :: list, token
      integer :: i, n, comma

      list = ' '
      n = 0
      DO i = 1, len_trim(raw_selector)
         IF (raw_selector(i:i) == ' ' .or. raw_selector(i:i) == char(9)) CYCLE
         n = n + 1
         list(n:n) = raw_selector(i:i)
      ENDDO
      list = material_lower(list)
      methane_history_all_land_only = .false.
      IF (len_trim(list) == 0) RETURN

      DO
         comma = index(list, ',')
         IF (comma > 0) THEN
            IF (comma == 1) RETURN
            token = list(:comma-1)
            list = list(comma+1:)
         ELSE
            token = list
         ENDIF
         SELECT CASE (trim(token))
         CASE ('f_methane_surf_flux_tot', &
               'f_methane_surf_flux_wetland', &
               'f_methane_surf_flux_soil', &
               'f_methane_surf_flux_lake', &
               'f_methane_surf_flux_rice', &
               'f_methane_surf_flux_active_total_without_lake', &
               'f_methane_surf_flux_global_total_with_lake', &
               'f_methane_surf_flux_global_phys_with_lake', &
               'f_methane_balance_residual_global_with_lake', &
               'f_methane_ch4_clip_credit_global_with_lake')
         CASE DEFAULT
            RETURN
         END SELECT
         IF (comma == 0) EXIT
      ENDDO
      methane_history_all_land_only = .true.

   END FUNCTION methane_history_all_land_only

END MODULE MOD_MT_Methane_Hist
#endif
