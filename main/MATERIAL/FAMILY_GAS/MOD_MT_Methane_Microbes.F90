#include <define.h>

#if (defined MaterialTransport) && (defined BGC)
MODULE MOD_MT_Methane_Microbes
!=======================================================================
! Optional microbial-pool diagnostics for the Methane material.
!
! This module is deliberately off by default.  When enabled with
! DEF_METHANE%use_microbial_pools it advances simple active/dormant
! methanogen and methanotroph biomass pools and exposes guarded potential
! CH4 production / oxidation diagnostics.  The legacy CLM/CTSM-style
! methane fluxes remain unchanged unless DEF_METHANE%use_microbial_flux_override
! is explicitly enabled in the physics driver.
!
! Coupling note: methane_driver advances microbial pools before methane()
! physics and passes the resulting potentials into the same land step.  The
! potentials are therefore an explicit-Euler tendency based on start-of-step
! CH4/O2/carbon states, not an implicit solve against post-transport
! concentrations.
! Biomass here is an activity-state proxy, not an additional conserved BGC
! carbon pool.  A shared local organic-C ceiling prevents oversized proxy
! pools; optional flux overrides are separately bounded by host substrate,
! CH4, and O2 availability.
!=======================================================================

   USE MOD_Precision
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_nan
   USE MOD_Vars_Global,        only: nl_soil, spval
   USE MOD_Const_Physical,     only: tfrz
   USE MOD_MT_Methane_Const, only: DEF_METHANE, secspday, catomw, gc_per_kg_om, &
      METHANE_COMP_SOIL, METHANE_COMP_RICE, N_METHANE_COMP

   IMPLICIT NONE
   SAVE
   PRIVATE

   ! Prognostic biomass pools [gC biomass / m3 soil]
   real(r8), allocatable :: B_methanogen(:,:)
   real(r8), allocatable :: B_methanotroph(:,:)
   real(r8), allocatable :: B_methanogen_dormant(:,:)
   real(r8), allocatable :: B_methanotroph_dormant(:,:)

   ! Per-step diagnostics [-], [day-1], [mol m-3 s-1]
   real(r8), allocatable :: f_T_methanogen(:,:)
   real(r8), allocatable :: f_S_methanogen(:,:)
   real(r8), allocatable :: f_O2_methanogen(:,:)
   real(r8), allocatable :: f_T_methanotroph(:,:)
   real(r8), allocatable :: methanogen_growth_rate(:,:)
   real(r8), allocatable :: methanotroph_growth_rate(:,:)
   real(r8), allocatable :: microbial_prod_potential(:,:)
   real(r8), allocatable :: microbial_oxid_potential(:,:)

   ! Lake methane follows WATERBODY substeps, whereas history accumulation is
   ! once per land timestep.  Keep only diagnostic means here; prognostic
   ! biomass pools remain at their final substep state.
   integer, parameter :: methane_lake_microbe_substep_ndiag = 8
   real(r8), allocatable :: methane_lake_microbe_substep_acc(:,:)
   integer :: methane_lake_microbe_substep_cached_ipatch = -1
   integer :: methane_lake_microbe_substep_next_isub = 1

   ! Independent soil/rice microbial columns.  The legacy 2-D arrays above
   ! remain area-aggregated diagnostics and restart-compatibility fields.
   real(r8), allocatable :: B_methanogen_comp(:,:,:)
   real(r8), allocatable :: B_methanotroph_comp(:,:,:)
   real(r8), allocatable :: B_methanogen_dormant_comp(:,:,:)
   real(r8), allocatable :: B_methanotroph_dormant_comp(:,:,:)
   real(r8), allocatable :: f_T_methanogen_comp(:,:,:)
   real(r8), allocatable :: f_S_methanogen_comp(:,:,:)
   real(r8), allocatable :: f_O2_methanogen_comp(:,:,:)
   real(r8), allocatable :: f_T_methanotroph_comp(:,:,:)
   real(r8), allocatable :: methanogen_growth_rate_comp(:,:,:)
   real(r8), allocatable :: methanotroph_growth_rate_comp(:,:,:)
   real(r8), allocatable :: microbial_prod_potential_comp(:,:,:)
   real(r8), allocatable :: microbial_oxid_potential_comp(:,:,:)

   integer, parameter :: N_MICROBE_STATE_RESTART_FIELDS = 12
   character(len=40), parameter :: MICROBE_STATE_RESTART_FIELDS(N_MICROBE_STATE_RESTART_FIELDS) = &
      [character(len=40) :: &
       'ch4_B_methanogen', 'ch4_B_methanotroph', &
       'ch4_B_methanogen_dormant', 'ch4_B_methanotroph_dormant', &
       'ch4_B_methanogen_soil', 'ch4_B_methanogen_rice', &
       'ch4_B_methanotroph_soil', 'ch4_B_methanotroph_rice', &
       'ch4_B_methanogen_dormant_soil', 'ch4_B_methanogen_dormant_rice', &
       'ch4_B_methanotroph_dormant_soil', 'ch4_B_methanotroph_dormant_rice']

   logical :: methane_microbes_lulcc_snapshot_valid = .false.
   real(r8), allocatable :: lulcc_B_methanogen_old(:,:)
   real(r8), allocatable :: lulcc_B_methanotroph_old(:,:)
   real(r8), allocatable :: lulcc_B_methanogen_dormant_old(:,:)
   real(r8), allocatable :: lulcc_B_methanotroph_dormant_old(:,:)
   real(r8), allocatable :: lulcc_f_T_methanogen_old(:,:)
   real(r8), allocatable :: lulcc_f_S_methanogen_old(:,:)
   real(r8), allocatable :: lulcc_f_O2_methanogen_old(:,:)
   real(r8), allocatable :: lulcc_f_T_methanotroph_old(:,:)
   real(r8), allocatable :: lulcc_methanogen_growth_rate_old(:,:)
   real(r8), allocatable :: lulcc_methanotroph_growth_rate_old(:,:)
   real(r8), allocatable :: lulcc_microbial_prod_potential_old(:,:)
   real(r8), allocatable :: lulcc_microbial_oxid_potential_old(:,:)
   real(r8), allocatable :: lulcc_B_methanogen_comp_old(:,:,:)
   real(r8), allocatable :: lulcc_B_methanotroph_comp_old(:,:,:)
   real(r8), allocatable :: lulcc_B_methanogen_dormant_comp_old(:,:,:)
   real(r8), allocatable :: lulcc_B_methanotroph_dormant_comp_old(:,:,:)
   real(r8), allocatable :: lulcc_f_T_methanogen_comp_old(:,:,:)
   real(r8), allocatable :: lulcc_f_S_methanogen_comp_old(:,:,:)
   real(r8), allocatable :: lulcc_f_O2_methanogen_comp_old(:,:,:)
   real(r8), allocatable :: lulcc_f_T_methanotroph_comp_old(:,:,:)
   real(r8), allocatable :: lulcc_methanogen_growth_rate_comp_old(:,:,:)
   real(r8), allocatable :: lulcc_methanotroph_growth_rate_comp_old(:,:,:)
   real(r8), allocatable :: lulcc_microbial_prod_potential_comp_old(:,:,:)
   real(r8), allocatable :: lulcc_microbial_oxid_potential_comp_old(:,:,:)
   real(r8), allocatable :: lulcc_rice_fraction_prev_old(:)

   PUBLIC :: allocate_methane_microbes_state
   PUBLIC :: accumulate_methane_lake_microbe_substep_diagnostics
   PUBLIC :: deallocate_methane_microbes_state
   PUBLIC :: reset_methane_inactive_lake_microbe_diagnostics
   PUBLIC :: methane_microbes_step
   PUBLIC :: aggregate_methane_microbes
   PUBLIC :: repartition_methane_microbes
   PUBLIC :: read_methane_microbes_restart
   PUBLIC :: validate_methane_microbes_restart_values
   PUBLIC :: write_methane_microbes_restart
   PUBLIC :: save_methane_microbes_lulcc_state
   PUBLIC :: remap_methane_microbes_lulcc_state
   PUBLIC :: B_methanogen, B_methanotroph
   PUBLIC :: B_methanogen_dormant, B_methanotroph_dormant
   PUBLIC :: f_T_methanogen, f_S_methanogen, f_O2_methanogen, f_T_methanotroph
   PUBLIC :: methanogen_growth_rate, methanotroph_growth_rate
   PUBLIC :: microbial_prod_potential, microbial_oxid_potential
   PUBLIC :: B_methanogen_comp, B_methanotroph_comp
   PUBLIC :: B_methanogen_dormant_comp, B_methanotroph_dormant_comp
   PUBLIC :: f_T_methanogen_comp, f_S_methanogen_comp
   PUBLIC :: f_O2_methanogen_comp, f_T_methanotroph_comp
   PUBLIC :: methanogen_growth_rate_comp, methanotroph_growth_rate_comp
   PUBLIC :: microbial_prod_potential_comp, microbial_oxid_potential_comp

CONTAINS

	   SUBROUTINE allocate_methane_microbes_state(numpatch)
	      integer, intent(in) :: numpatch

	      IF (.not. DEF_METHANE%use_microbial_pools) RETURN
	      IF (allocated(B_methanogen)) RETURN
	      ! Keep zero-length arrays allocated on ranks with numpatch==0 when
	      ! microbial pools are enabled; restart/vector I/O paths are collective.

	      allocate(B_methanogen(nl_soil,numpatch))
      allocate(B_methanotroph(nl_soil,numpatch))
      allocate(B_methanogen_dormant(nl_soil,numpatch))
      allocate(B_methanotroph_dormant(nl_soil,numpatch))
      allocate(f_T_methanogen(nl_soil,numpatch))
      allocate(f_S_methanogen(nl_soil,numpatch))
      allocate(f_O2_methanogen(nl_soil,numpatch))
      allocate(f_T_methanotroph(nl_soil,numpatch))
      allocate(methanogen_growth_rate(nl_soil,numpatch))
      allocate(methanotroph_growth_rate(nl_soil,numpatch))
      allocate(microbial_prod_potential(nl_soil,numpatch))
      allocate(microbial_oxid_potential(nl_soil,numpatch))
      allocate(B_methanogen_comp(nl_soil,N_METHANE_COMP,numpatch))
      allocate(B_methanotroph_comp(nl_soil,N_METHANE_COMP,numpatch))
      allocate(B_methanogen_dormant_comp(nl_soil,N_METHANE_COMP,numpatch))
      allocate(B_methanotroph_dormant_comp(nl_soil,N_METHANE_COMP,numpatch))
      allocate(f_T_methanogen_comp(nl_soil,N_METHANE_COMP,numpatch))
      allocate(f_S_methanogen_comp(nl_soil,N_METHANE_COMP,numpatch))
      allocate(f_O2_methanogen_comp(nl_soil,N_METHANE_COMP,numpatch))
      allocate(f_T_methanotroph_comp(nl_soil,N_METHANE_COMP,numpatch))
      allocate(methanogen_growth_rate_comp(nl_soil,N_METHANE_COMP,numpatch))
      allocate(methanotroph_growth_rate_comp(nl_soil,N_METHANE_COMP,numpatch))
      allocate(microbial_prod_potential_comp(nl_soil,N_METHANE_COMP,numpatch))
      allocate(microbial_oxid_potential_comp(nl_soil,N_METHANE_COMP,numpatch))

      B_methanogen(:,:) = DEF_METHANE%B_init_methanogen
      B_methanotroph(:,:) = DEF_METHANE%B_init_methanotroph
      B_methanogen_dormant(:,:) = 0._r8
      B_methanotroph_dormant(:,:) = 0._r8
      f_T_methanogen(:,:) = spval
      f_S_methanogen(:,:) = spval
      f_O2_methanogen(:,:) = spval
      f_T_methanotroph(:,:) = spval
      methanogen_growth_rate(:,:) = spval
      methanotroph_growth_rate(:,:) = spval
      microbial_prod_potential(:,:) = spval
      microbial_oxid_potential(:,:) = spval
      B_methanogen_comp = DEF_METHANE%B_init_methanogen
      B_methanotroph_comp = DEF_METHANE%B_init_methanotroph
      B_methanogen_dormant_comp = 0._r8
      B_methanotroph_dormant_comp = 0._r8
      f_T_methanogen_comp = spval
      f_S_methanogen_comp = spval
      f_O2_methanogen_comp = spval
      f_T_methanotroph_comp = spval
      methanogen_growth_rate_comp = spval
      methanotroph_growth_rate_comp = spval
      microbial_prod_potential_comp = spval
      microbial_oxid_potential_comp = spval
   END SUBROUTINE allocate_methane_microbes_state


   SUBROUTINE deallocate_methane_microbes_state()
      IF (allocated(B_methanogen)) deallocate(B_methanogen)
      IF (allocated(B_methanotroph)) deallocate(B_methanotroph)
      IF (allocated(B_methanogen_dormant)) deallocate(B_methanogen_dormant)
      IF (allocated(B_methanotroph_dormant)) deallocate(B_methanotroph_dormant)
      IF (allocated(f_T_methanogen)) deallocate(f_T_methanogen)
      IF (allocated(f_S_methanogen)) deallocate(f_S_methanogen)
      IF (allocated(f_O2_methanogen)) deallocate(f_O2_methanogen)
      IF (allocated(f_T_methanotroph)) deallocate(f_T_methanotroph)
      IF (allocated(methanogen_growth_rate)) deallocate(methanogen_growth_rate)
      IF (allocated(methanotroph_growth_rate)) deallocate(methanotroph_growth_rate)
      IF (allocated(microbial_prod_potential)) deallocate(microbial_prod_potential)
      IF (allocated(microbial_oxid_potential)) deallocate(microbial_oxid_potential)
      IF (allocated(B_methanogen_comp)) deallocate(B_methanogen_comp)
      IF (allocated(B_methanotroph_comp)) deallocate(B_methanotroph_comp)
      IF (allocated(B_methanogen_dormant_comp)) deallocate(B_methanogen_dormant_comp)
      IF (allocated(B_methanotroph_dormant_comp)) deallocate(B_methanotroph_dormant_comp)
      IF (allocated(f_T_methanogen_comp)) deallocate(f_T_methanogen_comp)
      IF (allocated(f_S_methanogen_comp)) deallocate(f_S_methanogen_comp)
      IF (allocated(f_O2_methanogen_comp)) deallocate(f_O2_methanogen_comp)
      IF (allocated(f_T_methanotroph_comp)) deallocate(f_T_methanotroph_comp)
      IF (allocated(methanogen_growth_rate_comp)) deallocate(methanogen_growth_rate_comp)
      IF (allocated(methanotroph_growth_rate_comp)) deallocate(methanotroph_growth_rate_comp)
      IF (allocated(microbial_prod_potential_comp)) deallocate(microbial_prod_potential_comp)
      IF (allocated(microbial_oxid_potential_comp)) deallocate(microbial_oxid_potential_comp)
      IF (allocated(methane_lake_microbe_substep_acc)) &
         deallocate(methane_lake_microbe_substep_acc)
      methane_lake_microbe_substep_cached_ipatch = -1
      methane_lake_microbe_substep_next_isub = 1
      CALL clear_methane_microbes_lulcc_snapshot()
   END SUBROUTINE deallocate_methane_microbes_state


   SUBROUTINE reset_methane_inactive_lake_microbe_diagnostics (ipatch)
      integer, intent(in) :: ipatch

      IF (.not. allocated(f_T_methanogen)) RETURN
      IF (ipatch < lbound(f_T_methanogen,2) .or. &
          ipatch > ubound(f_T_methanogen,2)) RETURN

      ! Biomass remains with the sediment carrier; only per-substep microbial
      ! rates and response factors are inactive when lake microbial pools are
      ! bypassed (wet lake) or the dynamic lake is dry.
      f_T_methanogen(:,ipatch)           = 0._r8
      f_S_methanogen(:,ipatch)           = 0._r8
      f_O2_methanogen(:,ipatch)          = 0._r8
      f_T_methanotroph(:,ipatch)         = 0._r8
      methanogen_growth_rate(:,ipatch)   = 0._r8
      methanotroph_growth_rate(:,ipatch) = 0._r8
      microbial_prod_potential(:,ipatch) = 0._r8
      microbial_oxid_potential(:,ipatch) = 0._r8

   END SUBROUTINE reset_methane_inactive_lake_microbe_diagnostics


   SUBROUTINE accumulate_methane_lake_microbe_substep_diagnostics &
      (ipatch, substep_dt, isub, nsub)
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: substep_dt
      integer,  intent(in) :: isub
      integer,  intent(in) :: nsub

      real(r8) :: total_dt

      IF (nsub <= 1) RETURN
      IF (substep_dt <= 0._r8) RETURN
      IF (.not. allocated(f_T_methanogen)) RETURN
      IF (ipatch < lbound(f_T_methanogen,2) .or. &
          ipatch > ubound(f_T_methanogen,2)) RETURN

      IF (.not. allocated(methane_lake_microbe_substep_acc)) THEN
         allocate(methane_lake_microbe_substep_acc( &
            nl_soil, methane_lake_microbe_substep_ndiag))
      ENDIF

      IF (isub <= 1 .or. &
          ipatch /= methane_lake_microbe_substep_cached_ipatch .or. &
          isub /= methane_lake_microbe_substep_next_isub) THEN
         methane_lake_microbe_substep_acc(:,:) = 0._r8
      ENDIF
      methane_lake_microbe_substep_cached_ipatch = ipatch

      CALL add2d(1, f_T_methanogen(:,ipatch))
      CALL add2d(2, f_S_methanogen(:,ipatch))
      CALL add2d(3, f_O2_methanogen(:,ipatch))
      CALL add2d(4, f_T_methanotroph(:,ipatch))
      CALL add2d(5, methanogen_growth_rate(:,ipatch))
      CALL add2d(6, methanotroph_growth_rate(:,ipatch))
      CALL add2d(7, microbial_prod_potential(:,ipatch))
      CALL add2d(8, microbial_oxid_potential(:,ipatch))

      IF (isub == nsub) THEN
         total_dt = substep_dt * real(nsub, r8)
         IF (total_dt > 0._r8) THEN
            CALL finish2d(1, f_T_methanogen(:,ipatch))
            CALL finish2d(2, f_S_methanogen(:,ipatch))
            CALL finish2d(3, f_O2_methanogen(:,ipatch))
            CALL finish2d(4, f_T_methanotroph(:,ipatch))
            CALL finish2d(5, methanogen_growth_rate(:,ipatch))
            CALL finish2d(6, methanotroph_growth_rate(:,ipatch))
            CALL finish2d(7, microbial_prod_potential(:,ipatch))
            CALL finish2d(8, microbial_oxid_potential(:,ipatch))
         ENDIF
         methane_lake_microbe_substep_cached_ipatch = -1
         methane_lake_microbe_substep_next_isub = 1
      ELSE
         methane_lake_microbe_substep_next_isub = isub + 1
      ENDIF

   CONTAINS

      SUBROUTINE add2d (icol, var)
         integer,  intent(in) :: icol
         real(r8), intent(in) :: var(1:nl_soil)

         methane_lake_microbe_substep_acc(:,icol) = &
            methane_lake_microbe_substep_acc(:,icol) + var(:) * substep_dt
      END SUBROUTINE add2d

      SUBROUTINE finish2d (icol, var)
         integer,  intent(in)    :: icol
         real(r8), intent(inout) :: var(1:nl_soil)

         var(:) = methane_lake_microbe_substep_acc(:,icol) / total_dt
      END SUBROUTINE finish2d

   END SUBROUTINE accumulate_methane_lake_microbe_substep_diagnostics


   SUBROUTINE methane_microbes_step(ipatch, component, deltim, t_soisno, conc_o2, conc_ch4, hr_vr, cellorg)
      integer,  intent(in) :: ipatch, component
      real(r8), intent(in) :: deltim
      real(r8), intent(in) :: t_soisno(1:nl_soil)
      real(r8), intent(in) :: conc_o2(1:nl_soil)
      real(r8), intent(in) :: conc_ch4(1:nl_soil)
      real(r8), intent(in) :: hr_vr(1:nl_soil)
      real(r8), intent(in) :: cellorg(1:nl_soil)

      integer :: j
      real(r8) :: dt_day, tempfac, sub_pool, sub_rate, f_s, f_o2, ch4_mm, o2_mm
	      real(r8) :: mu_m, mu_o, loss_m, loss_o, to_dormant, from_dormant
	      real(r8) :: prod_pot, oxid_pot, carbon_cap, ch4_cap, o2_cap, freeze_loss
	      real(r8) :: growth_factor, dormant_loss, dormant_available
	      real(r8) :: organic_c_layer
      real(r8), parameter :: small = 1.e-30_r8

      IF (.not. DEF_METHANE%use_microbial_pools) RETURN
      IF (.not. allocated(B_methanogen_comp)) RETURN
      IF (component < 1 .or. component > N_METHANE_COMP) RETURN
      IF (ipatch < lbound(B_methanogen_comp,3) .or. ipatch > ubound(B_methanogen_comp,3)) RETURN
      IF (deltim <= 0._r8) RETURN

      dt_day = deltim / secspday

      DO j = 1, nl_soil
	         organic_c_layer = max(cellorg(j), 0._r8) * gc_per_kg_om
	         IF (t_soisno(j) <= tfrz) THEN
	            freeze_loss = min(1._r8, max(0._r8, DEF_METHANE%gamma_microbial_freeze * dt_day))
	            B_methanogen_comp(j,component,ipatch) = max(DEF_METHANE%B_min_methanogen, &
	               B_methanogen_comp(j,component,ipatch) * (1._r8 - freeze_loss))
	            B_methanotroph_comp(j,component,ipatch) = max(DEF_METHANE%B_min_methanotroph, &
	               B_methanotroph_comp(j,component,ipatch) * (1._r8 - freeze_loss))
	            B_methanogen_dormant_comp(j,component,ipatch) = max(0._r8, &
	               B_methanogen_dormant_comp(j,component,ipatch) * (1._r8 - freeze_loss))
	            B_methanotroph_dormant_comp(j,component,ipatch) = max(0._r8, &
	               B_methanotroph_dormant_comp(j,component,ipatch) * (1._r8 - freeze_loss))
	            CALL cap_microbe_biomass_layer(j, component, ipatch, organic_c_layer)

            f_T_methanogen_comp(j,component,ipatch) = 0._r8
            f_S_methanogen_comp(j,component,ipatch) = 0._r8
            f_O2_methanogen_comp(j,component,ipatch) = 0._r8
            f_T_methanotroph_comp(j,component,ipatch) = 0._r8
            methanogen_growth_rate_comp(j,component,ipatch) = 0._r8
            methanotroph_growth_rate_comp(j,component,ipatch) = 0._r8
            microbial_prod_potential_comp(j,component,ipatch) = 0._r8
            microbial_oxid_potential_comp(j,component,ipatch) = 0._r8
            CYCLE
         ELSE
            tempfac = DEF_METHANE%q10_microbe_growth ** &
               ((t_soisno(j) - DEF_METHANE%T_ref_microbe) / 10._r8)
         ENDIF

         ! cellorg is [kg OM m-3]; convert to [mol C m-3] using the
         ! same 580 gC kgOM-1 convention as the wetland/lake surface-data
         ! fallbacks and BgcLink.
         sub_pool = max(cellorg(j), 0._r8) * gc_per_kg_om / catomw
         sub_rate = max(hr_vr(j) / catomw, 0._r8)
         f_s = sub_pool / (DEF_METHANE%K_substrate_methanogen_pool + sub_pool + small)
         f_o2 = DEF_METHANE%K_inh_O2_methanogen / &
            (DEF_METHANE%K_inh_O2_methanogen + max(conc_o2(j), 0._r8) + small)

	         mu_m = DEF_METHANE%mu_max_methanogen * tempfac * f_s * f_o2
	         loss_m = DEF_METHANE%gamma_methanogen * tempfac
	         to_dormant = 0._r8
	         from_dormant = 0._r8
	         growth_factor = exp(max(-50._r8, min(50._r8, (mu_m - loss_m) * dt_day)))
	         B_methanogen_comp(j,component,ipatch) = max(DEF_METHANE%B_min_methanogen, &
            B_methanogen_comp(j,component,ipatch) * growth_factor)
	         IF (DEF_METHANE%use_microbial_dormancy) THEN
	            IF (f_s < DEF_METHANE%dormancy_threshold_methanogen_fS .or. &
	                f_o2 < DEF_METHANE%dormancy_threshold_methanogen_fO2) THEN
	               ! B_min is a numerical seed pool, not transferable biomass.
	               ! Restrict dormancy transfer to the active surplus so moving
	               ! biomass cannot trigger floor replacement and create carbon.
	               to_dormant = min(max(B_methanogen_comp(j,component,ipatch) - &
	                  DEF_METHANE%B_min_methanogen, 0._r8), &
	                  max(0._r8, DEF_METHANE%dormancy_rate_active * &
	                  B_methanogen_comp(j,component,ipatch) * dt_day))
	            ELSE
	               from_dormant = min(B_methanogen_dormant_comp(j,component,ipatch), &
	                  max(0._r8, DEF_METHANE%dormancy_rate_revive * &
	                  B_methanogen_dormant_comp(j,component,ipatch) * dt_day))
		            ENDIF
		         ENDIF

	         B_methanogen_comp(j,component,ipatch) = max(DEF_METHANE%B_min_methanogen, &
	            B_methanogen_comp(j,component,ipatch) - to_dormant + from_dormant)
	         dormant_available = max(B_methanogen_dormant_comp(j,component,ipatch) - from_dormant, 0._r8)
	         dormant_loss = min(dormant_available, &
	            max(0._r8, DEF_METHANE%gamma_microbial_dormant * tempfac * &
	            dormant_available * dt_day))
	         B_methanogen_dormant_comp(j,component,ipatch) = dormant_available + to_dormant - dormant_loss
	         B_methanogen_dormant_comp(j,component,ipatch) = &
	            max(B_methanogen_dormant_comp(j,component,ipatch), 0._r8)

         ch4_mm = max(conc_ch4(j), 0._r8) / (DEF_METHANE%k_m + max(conc_ch4(j), 0._r8) + small)
         o2_mm = max(conc_o2(j), 0._r8) / (DEF_METHANE%k_m_o2 + max(conc_o2(j), 0._r8) + small)
         mu_o = DEF_METHANE%mu_max_methanotroph * tempfac * ch4_mm * o2_mm
         loss_o = DEF_METHANE%gamma_methanotroph * tempfac

	         growth_factor = exp(max(-50._r8, min(50._r8, (mu_o - loss_o) * dt_day)))
	         B_methanotroph_comp(j,component,ipatch) = max(DEF_METHANE%B_min_methanotroph, &
            B_methanotroph_comp(j,component,ipatch) * growth_factor)
	         to_dormant = 0._r8
	         from_dormant = 0._r8
	         IF (DEF_METHANE%use_microbial_dormancy) THEN
	            IF (ch4_mm < DEF_METHANE%dormancy_threshold_methanotroph_fS .or. &
	                o2_mm < DEF_METHANE%dormancy_threshold_methanotroph_fO2) THEN
	               to_dormant = min(max(B_methanotroph_comp(j,component,ipatch) - &
	                  DEF_METHANE%B_min_methanotroph, 0._r8), &
	                  max(0._r8, DEF_METHANE%dormancy_rate_active * &
	                  B_methanotroph_comp(j,component,ipatch) * dt_day))
	            ELSE
	               from_dormant = min(B_methanotroph_dormant_comp(j,component,ipatch), &
	                  max(0._r8, DEF_METHANE%dormancy_rate_revive * &
	                  B_methanotroph_dormant_comp(j,component,ipatch) * dt_day))
	            ENDIF
	         ENDIF
	         B_methanotroph_comp(j,component,ipatch) = max(DEF_METHANE%B_min_methanotroph, &
	            B_methanotroph_comp(j,component,ipatch) - to_dormant + from_dormant)
	         dormant_available = max(B_methanotroph_dormant_comp(j,component,ipatch) - from_dormant, 0._r8)
	         dormant_loss = min(dormant_available, &
	            max(0._r8, DEF_METHANE%gamma_microbial_dormant * tempfac * &
	            dormant_available * dt_day))
	         B_methanotroph_dormant_comp(j,component,ipatch) = dormant_available + to_dormant - dormant_loss
	         B_methanotroph_comp(j,component,ipatch) = &
	            max(B_methanotroph_comp(j,component,ipatch), DEF_METHANE%B_min_methanotroph)
	         B_methanotroph_dormant_comp(j,component,ipatch) = &
	            max(B_methanotroph_dormant_comp(j,component,ipatch), 0._r8)

	         CALL cap_microbe_biomass_layer(j, component, ipatch, organic_c_layer)

	         ! B_* pools are stored as [gC biomass m-3 soil].  Treat kappa_m_*
	         ! as first-order biomass-C turnover [day-1] and convert biomass C
	         ! to mol C with catomw before exposing molar CH4 potentials.  This
	         ! keeps the optional microbial override in the same [mol m-3 s-1]
	         ! units as the legacy production/oxidation physics and the caps
	         ! below.
	         prod_pot = DEF_METHANE%kappa_m_methanogen * B_methanogen_comp(j,component,ipatch) / catomw * &
	            tempfac * f_s * f_o2 / secspday
         carbon_cap = sub_rate
         prod_pot = min(max(prod_pot, 0._r8), carbon_cap)

	         oxid_pot = DEF_METHANE%kappa_m_methanotroph * B_methanotroph_comp(j,component,ipatch) / catomw * &
	            tempfac * ch4_mm * o2_mm / secspday
         ch4_cap = max(conc_ch4(j), 0._r8) / deltim
         o2_cap = max(conc_o2(j), 0._r8) / (2._r8 * deltim)
         oxid_pot = min(max(oxid_pot, 0._r8), ch4_cap, o2_cap)

         f_T_methanogen_comp(j,component,ipatch) = tempfac
         f_S_methanogen_comp(j,component,ipatch) = f_s
         f_O2_methanogen_comp(j,component,ipatch) = f_o2
         f_T_methanotroph_comp(j,component,ipatch) = tempfac
         methanogen_growth_rate_comp(j,component,ipatch) = mu_m
         methanotroph_growth_rate_comp(j,component,ipatch) = mu_o
         microbial_prod_potential_comp(j,component,ipatch) = prod_pot
         microbial_oxid_potential_comp(j,component,ipatch) = oxid_pot
      END DO
   END SUBROUTINE methane_microbes_step


   SUBROUTINE cap_microbe_biomass_layer(j, component, ipatch, organic_c_layer)
      ! The biomass pools are diagnostic proxies and are not debited from
      ! cellorg.  Enforce both per-guild limits and one shared organic-C
      ! ceiling so two guilds cannot each claim the full represented pool.
      integer, intent(in) :: j, component, ipatch
      real(r8), intent(in) :: organic_c_layer
      real(r8) :: cap_m, cap_o, total_cap, seed_total
      real(r8) :: surplus, allowed_surplus, scale

      cap_m = max(DEF_METHANE%B_min_methanogen, organic_c_layer)
      IF (DEF_METHANE%B_max_fraction_methanogen > 0._r8) &
         cap_m = min(cap_m, max(DEF_METHANE%B_min_methanogen, &
            DEF_METHANE%B_max_fraction_methanogen * organic_c_layer))
	  B_methanogen_comp(j,component,ipatch) = min(max(B_methanogen_comp(j,component,ipatch), &
         DEF_METHANE%B_min_methanogen), cap_m)
	  B_methanogen_dormant_comp(j,component,ipatch) = &
	     min(max(B_methanogen_dormant_comp(j,component,ipatch), 0._r8), &
	         max(cap_m - B_methanogen_comp(j,component,ipatch), 0._r8))

      cap_o = max(DEF_METHANE%B_min_methanotroph, organic_c_layer)
      IF (DEF_METHANE%B_max_fraction_methanotroph > 0._r8) &
         cap_o = min(cap_o, max(DEF_METHANE%B_min_methanotroph, &
            DEF_METHANE%B_max_fraction_methanotroph * organic_c_layer))
	  B_methanotroph_comp(j,component,ipatch) = min(max(B_methanotroph_comp(j,component,ipatch), &
         DEF_METHANE%B_min_methanotroph), cap_o)
	  B_methanotroph_dormant_comp(j,component,ipatch) = &
	     min(max(B_methanotroph_dormant_comp(j,component,ipatch), 0._r8), &
	         max(cap_o - B_methanotroph_comp(j,component,ipatch), 0._r8))

      seed_total = DEF_METHANE%B_min_methanogen + DEF_METHANE%B_min_methanotroph
      total_cap = max(seed_total, organic_c_layer)
	  surplus = max(B_methanogen_comp(j,component,ipatch) - DEF_METHANE%B_min_methanogen, 0._r8) + &
	         B_methanogen_dormant_comp(j,component,ipatch) + &
	         max(B_methanotroph_comp(j,component,ipatch) - DEF_METHANE%B_min_methanotroph, 0._r8) + &
	         B_methanotroph_dormant_comp(j,component,ipatch)
      allowed_surplus = max(total_cap - seed_total, 0._r8)
      IF (surplus > allowed_surplus .and. surplus > 0._r8) THEN
         scale = allowed_surplus / surplus
	     B_methanogen_comp(j,component,ipatch) = DEF_METHANE%B_min_methanogen + &
	            (B_methanogen_comp(j,component,ipatch) - DEF_METHANE%B_min_methanogen) * scale
	     B_methanogen_dormant_comp(j,component,ipatch) = &
	        B_methanogen_dormant_comp(j,component,ipatch) * scale
	     B_methanotroph_comp(j,component,ipatch) = DEF_METHANE%B_min_methanotroph + &
	            (B_methanotroph_comp(j,component,ipatch) - DEF_METHANE%B_min_methanotroph) * scale
	     B_methanotroph_dormant_comp(j,component,ipatch) = &
	        B_methanotroph_dormant_comp(j,component,ipatch) * scale
      ENDIF
   END SUBROUTINE cap_microbe_biomass_layer


   SUBROUTINE aggregate_methane_microbes(ipatch, rice_fraction)
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: rice_fraction
      real(r8) :: rice_weight

      IF (.not. allocated(B_methanogen_comp)) RETURN
      IF (ipatch < lbound(B_methanogen_comp,3) .or. ipatch > ubound(B_methanogen_comp,3)) RETURN

      rice_weight = min(1._r8, max(0._r8, rice_fraction))
      CALL aggregate_field(B_methanogen_comp, B_methanogen(:,ipatch))
      CALL aggregate_field(B_methanotroph_comp, B_methanotroph(:,ipatch))
      CALL aggregate_field(B_methanogen_dormant_comp, B_methanogen_dormant(:,ipatch))
      CALL aggregate_field(B_methanotroph_dormant_comp, B_methanotroph_dormant(:,ipatch))
      CALL aggregate_field(f_T_methanogen_comp, f_T_methanogen(:,ipatch))
      CALL aggregate_field(f_S_methanogen_comp, f_S_methanogen(:,ipatch))
      CALL aggregate_field(f_O2_methanogen_comp, f_O2_methanogen(:,ipatch))
      CALL aggregate_field(f_T_methanotroph_comp, f_T_methanotroph(:,ipatch))
      CALL aggregate_field(methanogen_growth_rate_comp, methanogen_growth_rate(:,ipatch))
      CALL aggregate_field(methanotroph_growth_rate_comp, methanotroph_growth_rate(:,ipatch))
      CALL aggregate_field(microbial_prod_potential_comp, microbial_prod_potential(:,ipatch))
      CALL aggregate_field(microbial_oxid_potential_comp, microbial_oxid_potential(:,ipatch))

   CONTAINS
      SUBROUTINE aggregate_field(component_values, aggregate_values)
         real(r8), intent(in)  :: component_values(:,:,:)
         real(r8), intent(out) :: aggregate_values(:)

         IF (rice_weight <= 0._r8) THEN
            aggregate_values = component_values(:,METHANE_COMP_SOIL,ipatch)
         ELSEIF (rice_weight >= 1._r8) THEN
            aggregate_values = component_values(:,METHANE_COMP_RICE,ipatch)
         ELSE
            aggregate_values = (1._r8 - rice_weight) * &
               component_values(:,METHANE_COMP_SOIL,ipatch) + rice_weight * &
               component_values(:,METHANE_COMP_RICE,ipatch)
         ENDIF
      END SUBROUTINE aggregate_field
   END SUBROUTINE aggregate_methane_microbes


   SUBROUTINE repartition_methane_microbes(ipatch, old_rice_fraction, new_rice_fraction)
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: old_rice_fraction, new_rice_fraction
      real(r8) :: old_rice, new_rice

      IF (.not. allocated(B_methanogen_comp)) RETURN
      IF (ipatch < lbound(B_methanogen_comp,3) .or. ipatch > ubound(B_methanogen_comp,3)) RETURN

      old_rice = min(1._r8, max(0._r8, old_rice_fraction))
      new_rice = min(1._r8, max(0._r8, new_rice_fraction))
      IF (.not. (new_rice > old_rice .or. new_rice < old_rice)) RETURN

      CALL repartition_field(B_methanogen_comp)
      CALL repartition_field(B_methanotroph_comp)
      CALL repartition_field(B_methanogen_dormant_comp)
      CALL repartition_field(B_methanotroph_dormant_comp)
      CALL repartition_field(f_T_methanogen_comp)
      CALL repartition_field(f_S_methanogen_comp)
      CALL repartition_field(f_O2_methanogen_comp)
      CALL repartition_field(f_T_methanotroph_comp)
      CALL repartition_field(methanogen_growth_rate_comp)
      CALL repartition_field(methanotroph_growth_rate_comp)
      CALL repartition_field(microbial_prod_potential_comp)
      CALL repartition_field(microbial_oxid_potential_comp)

   CONTAINS
      SUBROUTINE repartition_field(component_values)
         real(r8), intent(inout) :: component_values(:,:,:)

         IF (new_rice > old_rice) THEN
            component_values(:,METHANE_COMP_RICE,ipatch) = &
               (old_rice * component_values(:,METHANE_COMP_RICE,ipatch) + &
                (new_rice - old_rice) * component_values(:,METHANE_COMP_SOIL,ipatch)) / new_rice
         ELSE
            component_values(:,METHANE_COMP_SOIL,ipatch) = &
               ((1._r8 - old_rice) * component_values(:,METHANE_COMP_SOIL,ipatch) + &
                (old_rice - new_rice) * component_values(:,METHANE_COMP_RICE,ipatch)) / &
               (1._r8 - new_rice)
         ENDIF
      END SUBROUTINE repartition_field
   END SUBROUTINE repartition_methane_microbes


   SUBROUTINE write_methane_microbes_restart(file_restart, compress)
      USE MOD_LandPatch,    only: landpatch
      USE MOD_NetCDFVector, only: ncio_write_vector
      character(len=*), intent(in) :: file_restart
      integer,          intent(in) :: compress

      IF (.not. allocated(B_methanogen)) RETURN
      IF (.not. DEF_METHANE%use_microbial_pools) RETURN

      CALL ncio_write_vector(file_restart, 'ch4_B_methanogen', 'soil', nl_soil, &
         'patch', landpatch, B_methanogen, compress)
      CALL ncio_write_vector(file_restart, 'ch4_B_methanotroph', 'soil', nl_soil, &
         'patch', landpatch, B_methanotroph, compress)
      CALL ncio_write_vector(file_restart, 'ch4_B_methanogen_dormant', 'soil', nl_soil, &
         'patch', landpatch, B_methanogen_dormant, compress)
      CALL ncio_write_vector(file_restart, 'ch4_B_methanotroph_dormant', 'soil', nl_soil, &
         'patch', landpatch, B_methanotroph_dormant, compress)
      CALL ncio_write_vector(file_restart, 'ch4_B_methanogen_soil', 'soil', nl_soil, &
         'patch', landpatch, B_methanogen_comp(:,METHANE_COMP_SOIL,:), compress)
      CALL ncio_write_vector(file_restart, 'ch4_B_methanogen_rice', 'soil', nl_soil, &
         'patch', landpatch, B_methanogen_comp(:,METHANE_COMP_RICE,:), compress)
      CALL ncio_write_vector(file_restart, 'ch4_B_methanotroph_soil', 'soil', nl_soil, &
         'patch', landpatch, B_methanotroph_comp(:,METHANE_COMP_SOIL,:), compress)
      CALL ncio_write_vector(file_restart, 'ch4_B_methanotroph_rice', 'soil', nl_soil, &
         'patch', landpatch, B_methanotroph_comp(:,METHANE_COMP_RICE,:), compress)
      CALL ncio_write_vector(file_restart, 'ch4_B_methanogen_dormant_soil', 'soil', nl_soil, &
         'patch', landpatch, B_methanogen_dormant_comp(:,METHANE_COMP_SOIL,:), compress)
      CALL ncio_write_vector(file_restart, 'ch4_B_methanogen_dormant_rice', 'soil', nl_soil, &
         'patch', landpatch, B_methanogen_dormant_comp(:,METHANE_COMP_RICE,:), compress)
      CALL ncio_write_vector(file_restart, 'ch4_B_methanotroph_dormant_soil', 'soil', nl_soil, &
         'patch', landpatch, B_methanotroph_dormant_comp(:,METHANE_COMP_SOIL,:), compress)
      CALL ncio_write_vector(file_restart, 'ch4_B_methanotroph_dormant_rice', 'soil', nl_soil, &
         'patch', landpatch, B_methanotroph_dormant_comp(:,METHANE_COMP_RICE,:), compress)
   END SUBROUTINE write_methane_microbes_restart


   SUBROUTINE read_methane_microbes_restart(file_restart, strict_restart, require_component_state, &
         component_state_present)
      USE MOD_LandPatch,    only: landpatch
      USE MOD_MT_NcIO, only: ncio_read_vector => ncio_read_vector_complete, ncio_vector_group_presence
#ifdef USEMPI
      USE MOD_SPMD_Task, only: p_is_worker, p_is_master, p_comm_glb, p_err, &
         MPI_IN_PLACE, MPI_INTEGER, MPI_SUM, CoLM_stop
#else
      USE MOD_SPMD_Task, only: p_is_worker, p_is_master, CoLM_stop
#endif
      character(len=*), intent(in) :: file_restart
      logical, intent(in), optional :: strict_restart
      logical, intent(in), optional :: require_component_state
      logical, intent(in), optional :: component_state_present
      logical :: strict_restart_active, require_components
      logical :: component_fields_present(8)
      integer :: invalid_biomass_values

      IF (.not. allocated(B_methanogen)) RETURN
      IF (.not. DEF_METHANE%use_microbial_pools) RETURN
      strict_restart_active = .false.
      IF (present(strict_restart)) strict_restart_active = strict_restart
      require_components = .false.
      IF (present(require_component_state)) require_components = require_component_state

      CALL ncio_read_vector(file_restart, 'ch4_B_methanogen', nl_soil, landpatch, &
         B_methanogen, defval=DEF_METHANE%B_init_methanogen)
      CALL ncio_read_vector(file_restart, 'ch4_B_methanotroph', nl_soil, landpatch, &
         B_methanotroph, defval=DEF_METHANE%B_init_methanotroph)
      CALL ncio_read_vector(file_restart, 'ch4_B_methanogen_dormant', nl_soil, landpatch, &
         B_methanogen_dormant, defval=0._r8)
      CALL ncio_read_vector(file_restart, 'ch4_B_methanotroph_dormant', nl_soil, landpatch, &
         B_methanotroph_dormant, defval=0._r8)

      invalid_biomass_values = 0
      IF (p_is_worker) THEN
         invalid_biomass_values = &
            count(ieee_is_nan(B_methanogen) .or. abs(B_methanogen) > 1.e30_r8 .or. B_methanogen < 0._r8) + &
            count(ieee_is_nan(B_methanotroph) .or. abs(B_methanotroph) > 1.e30_r8 .or. B_methanotroph < 0._r8) + &
            count(ieee_is_nan(B_methanogen_dormant) .or. abs(B_methanogen_dormant) > 1.e30_r8 .or. &
                  B_methanogen_dormant < 0._r8) + &
            count(ieee_is_nan(B_methanotroph_dormant) .or. abs(B_methanotroph_dormant) > 1.e30_r8 .or. &
                  B_methanotroph_dormant < 0._r8)
      ENDIF
#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, invalid_biomass_values, 1, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif
      IF (invalid_biomass_values > 0) THEN
         IF (strict_restart_active) THEN
            IF (p_is_master) WRITE(*,'(A,I0,A)') 'ERROR: committed methane restart contains ', &
               invalid_biomass_values, ' invalid microbial biomass values.'
            CALL CoLM_stop()
         ELSEIF (p_is_master) THEN
            WRITE(*,'(A,I0,A)') 'WARNING: legacy methane restart sanitizes ', &
               invalid_biomass_values, ' invalid microbial biomass values.'
         ENDIF
      ENDIF

      ! Clean spval/NaN/negative biomass from restart.  defval only fires
      ! on completely-missing fields, but a malformed restart could still
      ! contain sentinels which would propagate into methane prod/oxid.
      WHERE (ieee_is_nan(B_methanogen) .or. abs(B_methanogen) > 1.e30_r8 .or. B_methanogen < 0._r8) &
         B_methanogen = DEF_METHANE%B_init_methanogen
      WHERE (ieee_is_nan(B_methanotroph) .or. abs(B_methanotroph) > 1.e30_r8 .or. B_methanotroph < 0._r8) &
         B_methanotroph = DEF_METHANE%B_init_methanotroph
      WHERE (ieee_is_nan(B_methanogen_dormant) .or. abs(B_methanogen_dormant) > 1.e30_r8 .or. B_methanogen_dormant < 0._r8) &
         B_methanogen_dormant = 0._r8
      WHERE (ieee_is_nan(B_methanotroph_dormant) .or. abs(B_methanotroph_dormant) > 1.e30_r8 .or. B_methanotroph_dormant < 0._r8) &
         B_methanotroph_dormant = 0._r8

      ! Legacy restarts contain only the aggregate fields.  Seed both columns
      ! from those values, then replace them when the complete component set is
      ! available in a new-format restart.
      B_methanogen_comp(:,METHANE_COMP_SOIL,:) = B_methanogen
      B_methanogen_comp(:,METHANE_COMP_RICE,:) = B_methanogen
      B_methanotroph_comp(:,METHANE_COMP_SOIL,:) = B_methanotroph
      B_methanotroph_comp(:,METHANE_COMP_RICE,:) = B_methanotroph
      B_methanogen_dormant_comp(:,METHANE_COMP_SOIL,:) = B_methanogen_dormant
      B_methanogen_dormant_comp(:,METHANE_COMP_RICE,:) = B_methanogen_dormant
      B_methanotroph_dormant_comp(:,METHANE_COMP_SOIL,:) = B_methanotroph_dormant
      B_methanotroph_dormant_comp(:,METHANE_COMP_RICE,:) = B_methanotroph_dormant

      IF (present(component_state_present)) THEN
         component_fields_present(:) = component_state_present
      ELSE
         CALL ncio_vector_group_presence(file_restart, MICROBE_STATE_RESTART_FIELDS(5:12), &
            landpatch, component_fields_present)
      ENDIF

      IF ((require_components .and. .not. all(component_fields_present)) .or. &
          (any(component_fields_present) .and. .not. all(component_fields_present))) THEN
         IF (p_is_master) WRITE(*,'(A)') &
            'ERROR: methane microbial component restart fields are incomplete.'
         CALL CoLM_stop()
      ENDIF

      IF (all(component_fields_present)) THEN
         CALL read_component_field('ch4_B_methanogen_soil', &
            B_methanogen_comp(:,METHANE_COMP_SOIL,:))
         CALL read_component_field('ch4_B_methanogen_rice', &
            B_methanogen_comp(:,METHANE_COMP_RICE,:))
         CALL read_component_field('ch4_B_methanotroph_soil', &
            B_methanotroph_comp(:,METHANE_COMP_SOIL,:))
         CALL read_component_field('ch4_B_methanotroph_rice', &
            B_methanotroph_comp(:,METHANE_COMP_RICE,:))
         CALL read_component_field('ch4_B_methanogen_dormant_soil', &
            B_methanogen_dormant_comp(:,METHANE_COMP_SOIL,:))
         CALL read_component_field('ch4_B_methanogen_dormant_rice', &
            B_methanogen_dormant_comp(:,METHANE_COMP_RICE,:))
         CALL read_component_field('ch4_B_methanotroph_dormant_soil', &
            B_methanotroph_dormant_comp(:,METHANE_COMP_SOIL,:))
         CALL read_component_field('ch4_B_methanotroph_dormant_rice', &
            B_methanotroph_dormant_comp(:,METHANE_COMP_RICE,:))
      ENDIF

      invalid_biomass_values = 0
      IF (p_is_worker) THEN
         invalid_biomass_values = &
            count(ieee_is_nan(B_methanogen_comp) .or. abs(B_methanogen_comp) > 1.e30_r8 .or. &
                  B_methanogen_comp < 0._r8) + &
            count(ieee_is_nan(B_methanotroph_comp) .or. abs(B_methanotroph_comp) > 1.e30_r8 .or. &
                  B_methanotroph_comp < 0._r8) + &
            count(ieee_is_nan(B_methanogen_dormant_comp) .or. &
                  abs(B_methanogen_dormant_comp) > 1.e30_r8 .or. B_methanogen_dormant_comp < 0._r8) + &
            count(ieee_is_nan(B_methanotroph_dormant_comp) .or. &
                  abs(B_methanotroph_dormant_comp) > 1.e30_r8 .or. B_methanotroph_dormant_comp < 0._r8)
      ENDIF
#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, invalid_biomass_values, 1, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif
      IF (invalid_biomass_values > 0) THEN
         IF (strict_restart_active) THEN
            IF (p_is_master) WRITE(*,'(A,I0,A)') 'ERROR: committed methane restart contains ', &
               invalid_biomass_values, ' invalid component microbial biomass values.'
            CALL CoLM_stop()
         ELSEIF (p_is_master) THEN
            WRITE(*,'(A,I0,A)') 'WARNING: legacy methane restart sanitizes ', &
               invalid_biomass_values, ' invalid component microbial biomass values.'
         ENDIF
      ENDIF

      WHERE (ieee_is_nan(B_methanogen_comp) .or. abs(B_methanogen_comp) > 1.e30_r8 .or. &
             B_methanogen_comp < 0._r8) B_methanogen_comp = DEF_METHANE%B_init_methanogen
      WHERE (ieee_is_nan(B_methanotroph_comp) .or. abs(B_methanotroph_comp) > 1.e30_r8 .or. &
             B_methanotroph_comp < 0._r8) B_methanotroph_comp = DEF_METHANE%B_init_methanotroph
      WHERE (ieee_is_nan(B_methanogen_dormant_comp) .or. &
             abs(B_methanogen_dormant_comp) > 1.e30_r8 .or. B_methanogen_dormant_comp < 0._r8) &
         B_methanogen_dormant_comp = 0._r8
      WHERE (ieee_is_nan(B_methanotroph_dormant_comp) .or. &
             abs(B_methanotroph_dormant_comp) > 1.e30_r8 .or. B_methanotroph_dormant_comp < 0._r8) &
         B_methanotroph_dormant_comp = 0._r8

   CONTAINS
      SUBROUTINE read_component_field(varname, target)
         character(len=*), intent(in) :: varname
         real(r8), intent(out) :: target(:,:)
         real(r8), allocatable :: values(:,:)

         allocate(values(nl_soil,size(target,2)))
         CALL ncio_read_vector(file_restart, varname, nl_soil, landpatch, values)
         target = values
         deallocate(values)
      END SUBROUTINE read_component_field
   END SUBROUTINE read_methane_microbes_restart


   SUBROUTINE validate_methane_microbes_restart_values(file_restart, strict_restart, include_components)
      ! Validate persisted optional microbial state even when the resumed run
      ! disables microbial pools.  Runtime flags may choose not to retain the
      ! state, but cannot weaken a committed restart transaction.
      USE MOD_LandPatch,    only: landpatch
      USE MOD_MT_NcIO, only: ncio_read_vector => ncio_read_vector_complete
#ifdef USEMPI
      USE MOD_SPMD_Task, only: p_is_worker, p_is_master, p_comm_glb, p_err, &
         MPI_IN_PLACE, MPI_INTEGER, MPI_SUM, CoLM_stop
#else
      USE MOD_SPMD_Task, only: p_is_worker, p_is_master, CoLM_stop
#endif
      character(len=*), intent(in) :: file_restart
      logical, intent(in) :: strict_restart, include_components
      real(r8), allocatable :: values(:,:)
      integer :: nfield, ifield, invalid_values

      nfield = merge(N_MICROBE_STATE_RESTART_FIELDS, 4, include_components)
      IF (p_is_worker) THEN
         allocate(values(nl_soil, landpatch%nset))
      ELSE
         ! landpatch%nset is worker-owned metadata and is not guaranteed to be
         ! initialized on master-only/I/O ranks.  All ranks still participate
         ! in vector I/O, so give non-workers an explicit zero extent.
         allocate(values(nl_soil, 0))
      ENDIF
      invalid_values = 0
      DO ifield = 1, nfield
         CALL ncio_read_vector(file_restart, MICROBE_STATE_RESTART_FIELDS(ifield), &
            nl_soil, landpatch, values)
         IF (p_is_worker) invalid_values = invalid_values + &
            count(ieee_is_nan(values) .or. abs(values) > 1.e30_r8 .or. values < 0._r8)
      ENDDO
      deallocate(values)
#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, invalid_values, 1, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif
      IF (invalid_values > 0) THEN
         IF (strict_restart) THEN
            IF (p_is_master) WRITE(*,'(A,I0,A)') 'ERROR: committed methane restart contains ', &
               invalid_values, ' invalid ignored microbial biomass values.'
            CALL CoLM_stop()
         ELSEIF (p_is_master) THEN
            WRITE(*,'(A,I0,A)') 'WARNING: legacy methane restart ignores ', &
               invalid_values, ' invalid microbial biomass values.'
         ENDIF
      ENDIF
   END SUBROUTINE validate_methane_microbes_restart_values


   SUBROUTINE save_methane_microbes_lulcc_state()
      USE MOD_MT_Methane_State, only: rice_fraction_prev
	  USE MOD_Vars_TimeInvariants, only: patchtype
	  integer :: ip, component, nsave
      IF (.not. allocated(B_methanogen)) THEN
         methane_microbes_lulcc_snapshot_valid = .false.
         RETURN
      ENDIF

      CALL clear_methane_microbes_lulcc_snapshot()
      allocate(lulcc_B_methanogen_old(nl_soil,size(B_methanogen,2)))
      lulcc_B_methanogen_old = B_methanogen
      allocate(lulcc_B_methanotroph_old(nl_soil,size(B_methanotroph,2)))
      lulcc_B_methanotroph_old = B_methanotroph
      allocate(lulcc_B_methanogen_dormant_old(nl_soil,size(B_methanogen_dormant,2)))
      lulcc_B_methanogen_dormant_old = B_methanogen_dormant
      allocate(lulcc_B_methanotroph_dormant_old(nl_soil,size(B_methanotroph_dormant,2)))
      lulcc_B_methanotroph_dormant_old = B_methanotroph_dormant
      allocate(lulcc_f_T_methanogen_old(nl_soil,size(f_T_methanogen,2)))
      lulcc_f_T_methanogen_old = f_T_methanogen
      allocate(lulcc_f_S_methanogen_old(nl_soil,size(f_S_methanogen,2)))
      lulcc_f_S_methanogen_old = f_S_methanogen
      allocate(lulcc_f_O2_methanogen_old(nl_soil,size(f_O2_methanogen,2)))
      lulcc_f_O2_methanogen_old = f_O2_methanogen
      allocate(lulcc_f_T_methanotroph_old(nl_soil,size(f_T_methanotroph,2)))
      lulcc_f_T_methanotroph_old = f_T_methanotroph
      allocate(lulcc_methanogen_growth_rate_old(nl_soil,size(methanogen_growth_rate,2)))
      lulcc_methanogen_growth_rate_old = methanogen_growth_rate
      allocate(lulcc_methanotroph_growth_rate_old(nl_soil,size(methanotroph_growth_rate,2)))
      lulcc_methanotroph_growth_rate_old = methanotroph_growth_rate
      allocate(lulcc_microbial_prod_potential_old(nl_soil,size(microbial_prod_potential,2)))
      lulcc_microbial_prod_potential_old = microbial_prod_potential
      allocate(lulcc_microbial_oxid_potential_old(nl_soil,size(microbial_oxid_potential,2)))
      lulcc_microbial_oxid_potential_old = microbial_oxid_potential
      allocate(lulcc_B_methanogen_comp_old, source=B_methanogen_comp)
      allocate(lulcc_B_methanotroph_comp_old, source=B_methanotroph_comp)
      allocate(lulcc_B_methanogen_dormant_comp_old, source=B_methanogen_dormant_comp)
      allocate(lulcc_B_methanotroph_dormant_comp_old, source=B_methanotroph_dormant_comp)
      allocate(lulcc_f_T_methanogen_comp_old, source=f_T_methanogen_comp)
      allocate(lulcc_f_S_methanogen_comp_old, source=f_S_methanogen_comp)
      allocate(lulcc_f_O2_methanogen_comp_old, source=f_O2_methanogen_comp)
      allocate(lulcc_f_T_methanotroph_comp_old, source=f_T_methanotroph_comp)
      allocate(lulcc_methanogen_growth_rate_comp_old, source=methanogen_growth_rate_comp)
      allocate(lulcc_methanotroph_growth_rate_comp_old, source=methanotroph_growth_rate_comp)
      allocate(lulcc_microbial_prod_potential_comp_old, source=microbial_prod_potential_comp)
      allocate(lulcc_microbial_oxid_potential_comp_old, source=microbial_oxid_potential_comp)
      allocate(lulcc_rice_fraction_prev_old, source=rice_fraction_prev)
	  nsave = min(size(patchtype), size(lulcc_rice_fraction_prev_old))
	  DO ip = 1, nsave
	     IF (patchtype(ip) == 0) CYCLE
	     lulcc_rice_fraction_prev_old(ip) = 0._r8
	     DO component = 1, N_METHANE_COMP
	        lulcc_B_methanogen_comp_old(:,component,ip) = lulcc_B_methanogen_old(:,ip)
	        lulcc_B_methanotroph_comp_old(:,component,ip) = lulcc_B_methanotroph_old(:,ip)
	        lulcc_B_methanogen_dormant_comp_old(:,component,ip) = lulcc_B_methanogen_dormant_old(:,ip)
	        lulcc_B_methanotroph_dormant_comp_old(:,component,ip) = lulcc_B_methanotroph_dormant_old(:,ip)
	        lulcc_f_T_methanogen_comp_old(:,component,ip) = lulcc_f_T_methanogen_old(:,ip)
	        lulcc_f_S_methanogen_comp_old(:,component,ip) = lulcc_f_S_methanogen_old(:,ip)
	        lulcc_f_O2_methanogen_comp_old(:,component,ip) = lulcc_f_O2_methanogen_old(:,ip)
	        lulcc_f_T_methanotroph_comp_old(:,component,ip) = lulcc_f_T_methanotroph_old(:,ip)
	        lulcc_methanogen_growth_rate_comp_old(:,component,ip) = lulcc_methanogen_growth_rate_old(:,ip)
	        lulcc_methanotroph_growth_rate_comp_old(:,component,ip) = lulcc_methanotroph_growth_rate_old(:,ip)
	        lulcc_microbial_prod_potential_comp_old(:,component,ip) = lulcc_microbial_prod_potential_old(:,ip)
	        lulcc_microbial_oxid_potential_comp_old(:,component,ip) = lulcc_microbial_oxid_potential_old(:,ip)
	     ENDDO
	  ENDDO
      methane_microbes_lulcc_snapshot_valid = .true.
   END SUBROUTINE save_methane_microbes_lulcc_state


	   SUBROUTINE remap_methane_microbes_lulcc_state(patchclass_new, eindex_new, &
	      patchclass_old, eindex_old, lccpct_patches, new_patch_area, old_patch_area)
	      USE MOD_Vars_TimeInvariants, only: patchtype
	      integer, intent(in) :: patchclass_new(:), patchclass_old(:)
	      integer*8, intent(in) :: eindex_new(:), eindex_old(:)
	      real(r8), intent(in), optional :: lccpct_patches(:,0:)
	      real(r8), intent(in), optional :: new_patch_area(:)
	      real(r8), intent(in), optional :: old_patch_area(:)
      integer :: nnew, np
	  real(r8) :: remapped_rice_fraction(size(patchclass_new))
	  integer, allocatable :: map_start(:), map_old(:), fallback_map(:)
	  real(r8), allocatable :: map_source_weight(:), map_mass_weight(:)
	  logical, allocatable :: map_mass_available(:)

	      nnew = size(patchclass_new)
	      CALL deallocate_methane_microbes_arrays_only()
	      CALL allocate_methane_microbes_state(nnew)
	      IF (.not. allocated(B_methanogen)) RETURN
	      IF (.not. methane_microbes_lulcc_snapshot_valid) RETURN

	  CALL build_lulcc_remap_map()

      CALL remap2d(lulcc_B_methanogen_old, B_methanogen)
      CALL remap2d(lulcc_B_methanotroph_old, B_methanotroph)
      CALL remap2d(lulcc_B_methanogen_dormant_old, B_methanogen_dormant)
      CALL remap2d(lulcc_B_methanotroph_dormant_old, B_methanotroph_dormant)
      CALL remap2d(lulcc_f_T_methanogen_old, f_T_methanogen)
      CALL remap2d(lulcc_f_S_methanogen_old, f_S_methanogen)
      CALL remap2d(lulcc_f_O2_methanogen_old, f_O2_methanogen)
      CALL remap2d(lulcc_f_T_methanotroph_old, f_T_methanotroph)
      CALL remap2d(lulcc_methanogen_growth_rate_old, methanogen_growth_rate)
      CALL remap2d(lulcc_methanotroph_growth_rate_old, methanotroph_growth_rate)
      CALL remap2d(lulcc_microbial_prod_potential_old, microbial_prod_potential)
      CALL remap2d(lulcc_microbial_oxid_potential_old, microbial_oxid_potential)
	  CALL remap_component_fraction(remapped_rice_fraction)
	  CALL remap_component3d(lulcc_B_methanogen_comp_old, B_methanogen_comp)
	  CALL remap_component3d(lulcc_B_methanotroph_comp_old, B_methanotroph_comp)
	  CALL remap_component3d(lulcc_B_methanogen_dormant_comp_old, B_methanogen_dormant_comp)
	  CALL remap_component3d(lulcc_B_methanotroph_dormant_comp_old, B_methanotroph_dormant_comp)
	  CALL remap_component3d(lulcc_f_T_methanogen_comp_old, f_T_methanogen_comp)
	  CALL remap_component3d(lulcc_f_S_methanogen_comp_old, f_S_methanogen_comp)
	  CALL remap_component3d(lulcc_f_O2_methanogen_comp_old, f_O2_methanogen_comp)
	  CALL remap_component3d(lulcc_f_T_methanotroph_comp_old, f_T_methanotroph_comp)
	  CALL remap_component3d(lulcc_methanogen_growth_rate_comp_old, methanogen_growth_rate_comp)
	  CALL remap_component3d(lulcc_methanotroph_growth_rate_comp_old, methanotroph_growth_rate_comp)
	  CALL remap_component3d(lulcc_microbial_prod_potential_comp_old, microbial_prod_potential_comp)
	  CALL remap_component3d(lulcc_microbial_oxid_potential_comp_old, microbial_oxid_potential_comp)
	  DO np = 1, min(nnew, size(remapped_rice_fraction))
	     CALL aggregate_methane_microbes(np, remapped_rice_fraction(np))
	     IF (np <= size(patchtype) .and. patchtype(np) /= 0) &
	        CALL mirror_microbe_components_from_aggregate(np)
	  ENDDO
      CALL clear_methane_microbes_lulcc_snapshot()

   CONTAINS
	  SUBROUTINE build_lulcc_remap_map()
	     integer :: ipnew, op, nold, nlink, link, c, rep, class_lo, class_hi
	     real(r8) :: base_weight, class_sum, target_area, denom
	     real(r8), allocatable :: old_group_class_area(:,:)
	     real(r8), allocatable :: new_target_class_area(:,:)
	     real(r8), allocatable :: target_group_class_area(:,:)
	     logical, allocatable :: old_group_area_ready(:)

	     nold = min(size(patchclass_old), size(eindex_old))
	     nlink = 0
	     DO ipnew = 1, min(nnew, size(eindex_new))
	        DO op = 1, nold
	           IF (eindex_old(op) == eindex_new(ipnew)) nlink = nlink + 1
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
	     DO ipnew = 1, nnew
	        map_start(ipnew) = link
	        map_mass_available(ipnew) = area_mass_remap_available(ipnew)
	        IF (ipnew <= size(eindex_new)) THEN
	           DO op = 1, nold
	              IF (eindex_old(op) /= eindex_new(ipnew)) CYCLE
	              map_old(link) = op
	              IF (fallback_map(ipnew) == 0 .and. ipnew <= size(patchclass_new)) THEN
	                 IF (patchclass_old(op) == patchclass_new(ipnew)) fallback_map(ipnew) = op
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

	     IF (present(new_patch_area)) THEN
	        DO ipnew = 1, min(nnew, size(lccpct_patches,1), size(new_patch_area))
	           class_sum = sum(max(lccpct_patches(ipnew,class_lo:class_hi), 0._r8))
	           IF (class_sum <= tiny(1._r8)) CYCLE
	           DO c = class_lo, class_hi
	              new_target_class_area(ipnew,c) = max(0._r8, new_patch_area(ipnew)) * &
	                 max(0._r8, lccpct_patches(ipnew,c)) / class_sum
	           ENDDO
	        ENDDO
	     ENDIF
	     DO ipnew = 1, nnew
	        IF (map_start(ipnew) >= map_start(ipnew+1)) CYCLE
	        rep = map_old(map_start(ipnew))
	        target_group_class_area(rep,:) = target_group_class_area(rep,:) + &
	           new_target_class_area(ipnew,:)
	     ENDDO

	     IF (present(old_patch_area)) THEN
	        DO ipnew = 1, nnew
	           IF (map_start(ipnew) >= map_start(ipnew+1)) CYCLE
	           rep = map_old(map_start(ipnew))
	           IF (old_group_area_ready(rep)) CYCLE
	           DO link = map_start(ipnew), map_start(ipnew+1)-1
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

	     DO ipnew = 1, min(nnew, size(lccpct_patches,1))
	        IF (map_start(ipnew) >= map_start(ipnew+1)) CYCLE
	        rep = map_old(map_start(ipnew))
	        DO link = map_start(ipnew), map_start(ipnew+1)-1
	           op = map_old(link)
	           c = patchclass_old(op)
	           IF (c < class_lo .or. c > class_hi) CYCLE
	           base_weight = max(0._r8, lccpct_patches(ipnew,c))
	           map_source_weight(link) = base_weight
	           IF (base_weight > 0._r8 .and. present(old_patch_area)) THEN
	              IF (op <= size(old_patch_area)) THEN
	                 denom = old_group_class_area(rep,c)
	                 IF (denom > 0._r8) map_source_weight(link) = base_weight * &
	                    max(0._r8, old_patch_area(op)) / denom
	              ENDIF
	           ENDIF

	           IF (map_mass_available(ipnew)) THEN
	              IF (op <= size(old_patch_area)) THEN
	                 target_area = new_target_class_area(ipnew,c)
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

	  SUBROUTINE remap_component_fraction(new_fraction)
	     real(r8), intent(out) :: new_fraction(:)
	     integer :: ipnew, op, src, link
	     real(r8) :: w, wsum, rice_area

	     new_fraction = 0._r8
	     DO ipnew = 1, min(size(new_fraction), nnew)
	        wsum = 0._r8
	        rice_area = 0._r8
	        IF (present(lccpct_patches)) THEN
	           DO link = map_start(ipnew), map_start(ipnew+1)-1
	              op = map_old(link)
	              IF (op > size(lulcc_rice_fraction_prev_old)) CYCLE
		              w = component_base_weight(ipnew, link)
	              IF (w <= 0._r8) CYCLE
	              rice_area = rice_area + w*min(max(lulcc_rice_fraction_prev_old(op),0._r8),1._r8)
	              wsum = wsum + w
	           ENDDO
	        ENDIF
	        IF (wsum > tiny(1._r8)) THEN
	           new_fraction(ipnew) = rice_area/wsum
	        ELSE
	           src = fallback_source(ipnew, size(lulcc_rice_fraction_prev_old))
	           IF (src > 0) new_fraction(ipnew) = lulcc_rice_fraction_prev_old(src)
	        ENDIF
	     ENDDO
	  END SUBROUTINE remap_component_fraction

	  SUBROUTINE mirror_microbe_components_from_aggregate(ipatch)
	     integer, intent(in) :: ipatch
	     integer :: component

	     DO component = 1, N_METHANE_COMP
	        B_methanogen_comp(:,component,ipatch) = B_methanogen(:,ipatch)
	        B_methanotroph_comp(:,component,ipatch) = B_methanotroph(:,ipatch)
	        B_methanogen_dormant_comp(:,component,ipatch) = B_methanogen_dormant(:,ipatch)
	        B_methanotroph_dormant_comp(:,component,ipatch) = B_methanotroph_dormant(:,ipatch)
	        f_T_methanogen_comp(:,component,ipatch) = f_T_methanogen(:,ipatch)
	        f_S_methanogen_comp(:,component,ipatch) = f_S_methanogen(:,ipatch)
	        f_O2_methanogen_comp(:,component,ipatch) = f_O2_methanogen(:,ipatch)
	        f_T_methanotroph_comp(:,component,ipatch) = f_T_methanotroph(:,ipatch)
	        methanogen_growth_rate_comp(:,component,ipatch) = methanogen_growth_rate(:,ipatch)
	        methanotroph_growth_rate_comp(:,component,ipatch) = methanotroph_growth_rate(:,ipatch)
	        microbial_prod_potential_comp(:,component,ipatch) = microbial_prod_potential(:,ipatch)
	        microbial_oxid_potential_comp(:,component,ipatch) = microbial_oxid_potential(:,ipatch)
	     ENDDO
	  END SUBROUTINE mirror_microbe_components_from_aggregate

	  SUBROUTINE remap_component3d(old, new)
	     real(r8), intent(in) :: old(:,:,:)
	     real(r8), intent(inout) :: new(:,:,:)
	     integer :: ipnew, op, component, src, nlev, link
	     real(r8) :: w, wsum, r
	     real(r8) :: val(size(new,1))

	     nlev = min(size(new,1), size(old,1))
	     DO ipnew = 1, min(size(new,3), nnew)
	        DO component = 1, min(size(new,2), size(old,2), N_METHANE_COMP)
	           wsum = 0._r8
	           val = 0._r8
	           IF (present(lccpct_patches)) THEN
	              DO link = map_start(ipnew), map_start(ipnew+1)-1
	                 op = map_old(link)
	                 IF (op > size(old,3) .or. &
	                     op > size(lulcc_rice_fraction_prev_old)) CYCLE
		                 w = component_base_weight(ipnew, link)
	                 r = min(max(lulcc_rice_fraction_prev_old(op), 0._r8), 1._r8)
	                 IF (component == METHANE_COMP_RICE) THEN
	                    w = w*r
	                 ELSE
	                    w = w*(1._r8-r)
	                 ENDIF
	                 IF (w <= 0._r8) CYCLE
	                 val(1:nlev) = val(1:nlev) + w*old(1:nlev,component,op)
	                 wsum = wsum + w
	              ENDDO
	           ENDIF
	           IF (wsum > tiny(1._r8)) THEN
	              new(1:nlev,component,ipnew) = val(1:nlev)/wsum
	           ELSE
	              src = fallback_source(ipnew, size(old,3))
	              IF (src > 0) new(1:nlev,component,ipnew) = old(1:nlev,component,src)
	           ENDIF
	        ENDDO
	     ENDDO
	  END SUBROUTINE remap_component3d

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
               new(:,np) = default_vals(:)
               src = fallback_source(np, size(old,2))
               IF (src > 0) THEN
                  new(1:min(size(new,1),size(old,1)),np) = &
                     old(1:min(size(new,1),size(old,1)),src)
               ENDIF
            ENDIF
         ENDDO
      END SUBROUTINE remap2d
      integer FUNCTION fallback_source(np, nold) RESULT(src)
         integer, intent(in) :: np, nold
         src = 0
	     IF (np < 1 .or. np > size(fallback_map)) RETURN
	     IF (fallback_map(np) <= nold) src = fallback_map(np)
         ! Microbial pools are patch-class specific; avoid copying a same-eindex
         ! pool from a different land-cover class when no class match exists.
	      END FUNCTION fallback_source

		      REAL(r8) FUNCTION component_base_weight(np, link) RESULT(w)
		         integer, intent(in) :: np, link

		         IF (map_mass_available(np)) THEN
		            w = map_mass_weight(link)
		         ELSE
		            w = map_source_weight(link)
		         ENDIF
		      END FUNCTION component_base_weight

		      LOGICAL FUNCTION area_mass_remap_available(np) RESULT(ok)
		         integer, intent(in) :: np

		         ok = present(lccpct_patches) .and. present(old_patch_area) .and. present(new_patch_area)
		         IF (.not. ok) RETURN
		         ok = np <= size(new_patch_area)
		         IF (.not. ok) RETURN
		         ok = new_patch_area(np) > tiny(1._r8)
		      END FUNCTION area_mass_remap_available

		   END SUBROUTINE remap_methane_microbes_lulcc_state


   SUBROUTINE deallocate_methane_microbes_arrays_only()
      IF (allocated(B_methanogen)) deallocate(B_methanogen)
      IF (allocated(B_methanotroph)) deallocate(B_methanotroph)
      IF (allocated(B_methanogen_dormant)) deallocate(B_methanogen_dormant)
      IF (allocated(B_methanotroph_dormant)) deallocate(B_methanotroph_dormant)
      IF (allocated(f_T_methanogen)) deallocate(f_T_methanogen)
      IF (allocated(f_S_methanogen)) deallocate(f_S_methanogen)
      IF (allocated(f_O2_methanogen)) deallocate(f_O2_methanogen)
      IF (allocated(f_T_methanotroph)) deallocate(f_T_methanotroph)
      IF (allocated(methanogen_growth_rate)) deallocate(methanogen_growth_rate)
      IF (allocated(methanotroph_growth_rate)) deallocate(methanotroph_growth_rate)
      IF (allocated(microbial_prod_potential)) deallocate(microbial_prod_potential)
      IF (allocated(microbial_oxid_potential)) deallocate(microbial_oxid_potential)
      IF (allocated(B_methanogen_comp)) deallocate(B_methanogen_comp)
      IF (allocated(B_methanotroph_comp)) deallocate(B_methanotroph_comp)
      IF (allocated(B_methanogen_dormant_comp)) deallocate(B_methanogen_dormant_comp)
      IF (allocated(B_methanotroph_dormant_comp)) deallocate(B_methanotroph_dormant_comp)
      IF (allocated(f_T_methanogen_comp)) deallocate(f_T_methanogen_comp)
      IF (allocated(f_S_methanogen_comp)) deallocate(f_S_methanogen_comp)
      IF (allocated(f_O2_methanogen_comp)) deallocate(f_O2_methanogen_comp)
      IF (allocated(f_T_methanotroph_comp)) deallocate(f_T_methanotroph_comp)
      IF (allocated(methanogen_growth_rate_comp)) deallocate(methanogen_growth_rate_comp)
      IF (allocated(methanotroph_growth_rate_comp)) deallocate(methanotroph_growth_rate_comp)
      IF (allocated(microbial_prod_potential_comp)) deallocate(microbial_prod_potential_comp)
      IF (allocated(microbial_oxid_potential_comp)) deallocate(microbial_oxid_potential_comp)
   END SUBROUTINE deallocate_methane_microbes_arrays_only


   SUBROUTINE clear_methane_microbes_lulcc_snapshot()
      IF (allocated(lulcc_B_methanogen_old)) deallocate(lulcc_B_methanogen_old)
      IF (allocated(lulcc_B_methanotroph_old)) deallocate(lulcc_B_methanotroph_old)
      IF (allocated(lulcc_B_methanogen_dormant_old)) deallocate(lulcc_B_methanogen_dormant_old)
      IF (allocated(lulcc_B_methanotroph_dormant_old)) deallocate(lulcc_B_methanotroph_dormant_old)
      IF (allocated(lulcc_f_T_methanogen_old)) deallocate(lulcc_f_T_methanogen_old)
      IF (allocated(lulcc_f_S_methanogen_old)) deallocate(lulcc_f_S_methanogen_old)
      IF (allocated(lulcc_f_O2_methanogen_old)) deallocate(lulcc_f_O2_methanogen_old)
      IF (allocated(lulcc_f_T_methanotroph_old)) deallocate(lulcc_f_T_methanotroph_old)
      IF (allocated(lulcc_methanogen_growth_rate_old)) deallocate(lulcc_methanogen_growth_rate_old)
      IF (allocated(lulcc_methanotroph_growth_rate_old)) deallocate(lulcc_methanotroph_growth_rate_old)
      IF (allocated(lulcc_microbial_prod_potential_old)) deallocate(lulcc_microbial_prod_potential_old)
      IF (allocated(lulcc_microbial_oxid_potential_old)) deallocate(lulcc_microbial_oxid_potential_old)
      IF (allocated(lulcc_B_methanogen_comp_old)) deallocate(lulcc_B_methanogen_comp_old)
      IF (allocated(lulcc_B_methanotroph_comp_old)) deallocate(lulcc_B_methanotroph_comp_old)
      IF (allocated(lulcc_B_methanogen_dormant_comp_old)) deallocate(lulcc_B_methanogen_dormant_comp_old)
      IF (allocated(lulcc_B_methanotroph_dormant_comp_old)) deallocate(lulcc_B_methanotroph_dormant_comp_old)
      IF (allocated(lulcc_f_T_methanogen_comp_old)) deallocate(lulcc_f_T_methanogen_comp_old)
      IF (allocated(lulcc_f_S_methanogen_comp_old)) deallocate(lulcc_f_S_methanogen_comp_old)
      IF (allocated(lulcc_f_O2_methanogen_comp_old)) deallocate(lulcc_f_O2_methanogen_comp_old)
      IF (allocated(lulcc_f_T_methanotroph_comp_old)) deallocate(lulcc_f_T_methanotroph_comp_old)
      IF (allocated(lulcc_methanogen_growth_rate_comp_old)) deallocate(lulcc_methanogen_growth_rate_comp_old)
      IF (allocated(lulcc_methanotroph_growth_rate_comp_old)) deallocate(lulcc_methanotroph_growth_rate_comp_old)
      IF (allocated(lulcc_microbial_prod_potential_comp_old)) deallocate(lulcc_microbial_prod_potential_comp_old)
      IF (allocated(lulcc_microbial_oxid_potential_comp_old)) deallocate(lulcc_microbial_oxid_potential_comp_old)
      IF (allocated(lulcc_rice_fraction_prev_old)) deallocate(lulcc_rice_fraction_prev_old)
      methane_microbes_lulcc_snapshot_valid = .false.
   END SUBROUTINE clear_methane_microbes_lulcc_snapshot

END MODULE MOD_MT_Methane_Microbes
#endif
