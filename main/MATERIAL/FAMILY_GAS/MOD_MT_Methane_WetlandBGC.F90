#include <define.h>

#if (defined MaterialTransport) && (defined BGC)
MODULE MOD_MT_Methane_WetlandBGC
!=======================================================================
! !DESCRIPTION:
!
!  Soil decomposition of wetland patches (patchtype 2) for methane.
!
!  Upstream CoLM runs the BGC soil cascade on soil patches only, so a
!  wetland patch has no heterotrophic respiration for methane production.
!  CoLM202X (MOD_Material_Reactive_BgcShim) ran the cascade on the wetland
!  column of the host BGC arrays and updated those pools, which changes the
!  host BGC state, restart and history of every wetland patch.  Here the
!  host is left as upstream has it: the wetland carbon and nitrogen pools
!  are methane state (shadow pools, written to the methane restart), and
!  for each wetland patch and time step
!
!    1. the host column ipatch of every BGC array the cascade touches is
!       saved;
!    2. the shadow pools are put into the host column and the flux state
!       is cleared, as CoLM202X did;
!    3. the upstream routines decomp_rate_constants_bgc,
!       SoilBiogeochemPotential and SoilBiogeochemDecomp run, with the
!       no-plant competition below in between (CoLM202X added it to
!       MOD_BGC_Soil_BiogeochemCompetition);
!    4. what methane reads (respiration by transition, o_scalar,
!       pot_f_nit_vr) and what the pool update needs is copied out;
!    5. the host column is restored, so no host value changes.
!
!  The pool update is applied to the shadow pools when methane finalises
!  the step (methane_wetland_bgc_apply, from mt_ch4_bgc_finalize_step), as
!  CoLM202X did with CDecompStateUpdate / SoilBiogeochemNDecompStateUpdate
!  (same order of operations; the SASU accumulators are not updated, they
!  belong to the host soil patches).  Methane reads the pools of the step
!  before the update, as before.
!
!  Cold start: the shadow pools start from the host values of the wetland
!  column (zero after mkinidata unless the CN initial file is used); a
!  layer without carbon gets the organic matter of the soil data
!  (OM_density * 580 gC/kg, split 5/10/5/0/5/25/50 % into the litter, CWD
!  and soil pools, C:N 15).  CoLM202X mkinidata applied the same split,
!  but only inside its DEF_USE_CN_INIT branch, so with the default
!  (.false.) its wetland columns started without carbon and produced
!  almost no methane.  Here the organic matter start is applied in every
!  cold start, so that the default configuration has wetland substrate.
!  CoLM202X also filled the wetland column from the CN initial file; that
!  file is read by mkinidata only and is not used here.
!=======================================================================

   USE MOD_Precision
   USE MOD_SPMD_Task, only: CoLM_stop
   USE MOD_Vars_Global, only: nl_soil, nl_soil_full, z_soi, dz_soi, &
      ndecomp_pools, ndecomp_transitions, spval
   USE MOD_BGC_Soil_BiogeochemDecompCascadeBGC, only: decomp_rate_constants_bgc
   USE MOD_BGC_Soil_BiogeochemPotential,        only: SoilBiogeochemPotential
   USE MOD_BGC_Soil_BiogeochemDecomp,           only: SoilBiogeochemDecomp
   USE MOD_BGC_Vars_TimeInvariants, only: donor_pool, receiver_pool, &
      i_met_lit, i_cel_lit, i_lig_lit, i_cwd, i_soil1, i_soil2, i_soil3
   USE MOD_BGC_Vars_TimeVariables, only: decomp_cpools_vr, decomp_npools_vr, &
      sminn_vr, smin_nh4_vr, smin_no3_vr, decomp_k, t_scalar, w_scalar, o_scalar, &
      depth_scalar, cn_decomp_pools, fpi_vr, fpi, fpg
   USE MOD_BGC_Vars_1DFluxes, only: gross_nmin_vr, p_decomp_cpool_loss, phr_vr, &
      pmnf_decomp, potential_immob_vr, decomp_ctransfer_vr, decomp_hr_vr, &
      decomp_ntransfer_vr, decomp_sminn_flux_vr, gross_nmin, net_nmin, net_nmin_vr, &
      sminn_to_denit_decomp_vr, decomp_cpools_sourcesink, decomp_npools_sourcesink, &
      pot_f_nit_vr, decomp_hr, somc_fire, som_c_leached, som_n_leached, denit, &
      f_n2o_nit, smin_no3_leached, smin_no3_runoff, sminn_leached, sminn_to_plant, &
      sminn_to_plant_vr, smin_nh4_to_plant_vr, smin_no3_to_plant_vr, &
      supplement_to_sminn_vr, sminn_to_denit_excess_vr, f_nit_vr, f_denit_vr, &
      actual_immob_vr, actual_immob_nh4_vr, actual_immob_no3_vr

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: methane_wetland_bgc_allocate, methane_wetland_bgc_deallocate
   PUBLIC :: methane_wetland_bgc_cold_init
   PUBLIC :: methane_wetland_bgc_step, methane_wetland_bgc_apply
   PUBLIC :: methane_wetland_bgc_view_active
   PUBLIC :: methane_wetland_bgc_write_restart, methane_wetland_bgc_read_restart
   PUBLIC :: methane_wetland_bgc_save_lulcc, methane_wetland_bgc_remap_lulcc

   ! shadow pools of the wetland patches (methane state)
   real(r8), allocatable, PUBLIC :: wl_cpools   (:,:,:)  ! (nl_soil_full,ndecomp_pools,numpatch) gC/m3
   real(r8), allocatable, PUBLIC :: wl_npools   (:,:,:)  ! (nl_soil_full,ndecomp_pools,numpatch) gN/m3
   real(r8), allocatable, PUBLIC :: wl_sminn    (:,:)    ! (nl_soil,numpatch) gN/m3
   real(r8), allocatable, PUBLIC :: wl_smin_nh4 (:,:)
   real(r8), allocatable, PUBLIC :: wl_smin_no3 (:,:)

   ! what methane reads for the patch of the current step (view), and the
   ! fluxes of the pending pool update
   integer,  save :: view_patch = -1
   logical,  save :: pending = .false.
   real(r8), allocatable, PUBLIC :: wv_hr_vr       (:,:)  ! (nl_soil,ndecomp_transitions)
   real(r8), allocatable, PUBLIC :: wv_o_scalar    (:)
   real(r8), allocatable, PUBLIC :: wv_pot_f_nit_vr(:)
   real(r8), allocatable :: pd_ctransfer (:,:), pd_ntransfer (:,:), pd_sminn_flux (:,:), &
      pd_sminn_denit (:,:)
   real(r8), allocatable :: pd_gross_nmin (:), pd_immob_nh4 (:), pd_immob_no3 (:)

   ! saved host column
   real(r8), allocatable :: s3_cpools(:,:), s3_npools(:,:), s3_decomp_k(:,:), s3_cn_pools(:,:), &
      s3_p_loss(:,:), s3_pmnf(:,:), s3_ctransfer(:,:), s3_hr(:,:), s3_ntransfer(:,:), &
      s3_sminn_flux(:,:), s3_sminn_denit(:,:), s3_css(:,:), s3_nss(:,:)
   real(r8), allocatable :: s2_sminn(:), s2_nh4(:), s2_no3(:), s2_t_scalar(:), s2_w_scalar(:), &
      s2_o_scalar(:), s2_depth_scalar(:), s2_fpi_vr(:), s2_gross_nmin_vr(:), s2_phr_vr(:), &
      s2_potential_immob_vr(:), s2_net_nmin_vr(:), s2_pot_f_nit_vr(:), s2_sminn_to_plant_vr(:), &
      s2_nh4_to_plant_vr(:), s2_no3_to_plant_vr(:), s2_supplement(:), s2_denit_excess(:), &
      s2_f_nit_vr(:), s2_f_denit_vr(:), s2_actual_immob_vr(:), s2_immob_nh4(:), s2_immob_no3(:)
   ! shadow pools of the old patch layout across a LULCC remap
   real(r8), allocatable :: lulcc_cpools(:,:,:), lulcc_npools(:,:,:), lulcc_sminn(:,:), &
      lulcc_smin_nh4(:,:), lulcc_smin_no3(:,:)

   real(r8) :: s1_fpi, s1_fpg, s1_gross_nmin, s1_net_nmin, s1_decomp_hr, s1_somc_fire, &
      s1_som_c_leached, s1_som_n_leached, s1_denit, s1_f_n2o_nit, s1_no3_leached, &
      s1_no3_runoff, s1_sminn_leached, s1_sminn_to_plant

   real(r8), parameter :: wetland_cn_ratio = 15._r8     ! CoLM202X mkinidata
   real(r8), parameter :: carbon_per_kg_om = 580._r8

CONTAINS

   SUBROUTINE methane_wetland_bgc_allocate (numpatch)

      IMPLICIT NONE
      integer, intent(in) :: numpatch

      CALL methane_wetland_bgc_deallocate ()
      allocate (wl_cpools   (nl_soil_full,ndecomp_pools,numpatch)); wl_cpools   = 0._r8
      allocate (wl_npools   (nl_soil_full,ndecomp_pools,numpatch)); wl_npools   = 0._r8
      allocate (wl_sminn    (nl_soil,numpatch));                    wl_sminn    = 0._r8
      allocate (wl_smin_nh4 (nl_soil,numpatch));                    wl_smin_nh4 = 0._r8
      allocate (wl_smin_no3 (nl_soil,numpatch));                    wl_smin_no3 = 0._r8

      allocate (wv_hr_vr        (nl_soil,ndecomp_transitions))
      allocate (wv_o_scalar     (nl_soil))
      allocate (wv_pot_f_nit_vr (nl_soil))
      allocate (pd_ctransfer   (nl_soil,ndecomp_transitions))
      allocate (pd_ntransfer   (nl_soil,ndecomp_transitions))
      allocate (pd_sminn_flux  (nl_soil,ndecomp_transitions))
      allocate (pd_sminn_denit (nl_soil,ndecomp_transitions))
      allocate (pd_gross_nmin (nl_soil), pd_immob_nh4 (nl_soil), pd_immob_no3 (nl_soil))
      view_patch = -1
      pending = .false.

   END SUBROUTINE methane_wetland_bgc_allocate

   SUBROUTINE methane_wetland_bgc_deallocate ()

      IMPLICIT NONE

      IF (allocated(wl_cpools  )) deallocate (wl_cpools  )
      IF (allocated(wl_npools  )) deallocate (wl_npools  )
      IF (allocated(wl_sminn   )) deallocate (wl_sminn   )
      IF (allocated(wl_smin_nh4)) deallocate (wl_smin_nh4)
      IF (allocated(wl_smin_no3)) deallocate (wl_smin_no3)
      IF (allocated(wv_hr_vr       )) deallocate (wv_hr_vr       )
      IF (allocated(wv_o_scalar    )) deallocate (wv_o_scalar    )
      IF (allocated(wv_pot_f_nit_vr)) deallocate (wv_pot_f_nit_vr)
      IF (allocated(pd_ctransfer  )) deallocate (pd_ctransfer  )
      IF (allocated(pd_ntransfer  )) deallocate (pd_ntransfer  )
      IF (allocated(pd_sminn_flux )) deallocate (pd_sminn_flux )
      IF (allocated(pd_sminn_denit)) deallocate (pd_sminn_denit)
      IF (allocated(pd_gross_nmin )) deallocate (pd_gross_nmin )
      IF (allocated(pd_immob_nh4  )) deallocate (pd_immob_nh4  )
      IF (allocated(pd_immob_no3  )) deallocate (pd_immob_no3  )
      view_patch = -1
      pending = .false.

   END SUBROUTINE methane_wetland_bgc_deallocate

   ! -------
   ! Shadow pools of the wetland patches at a cold start (see the header).
   SUBROUTINE methane_wetland_bgc_cold_init (patchtype, numpatch)

      USE MOD_Vars_TimeInvariants, only: OM_density
      IMPLICIT NONE
      integer, intent(in) :: numpatch
      integer, intent(in) :: patchtype(:)

      integer  :: i, j, k
      real(r8) :: existing_c, tot_c

      IF (.not. allocated(wl_cpools)) RETURN

      DO i = 1, numpatch
         IF (patchtype(i) /= 2) CYCLE

         DO k = 1, ndecomp_pools
            DO j = 1, nl_soil_full
               wl_cpools(j,k,i) = valid_or (decomp_cpools_vr(j,k,i), 0._r8)
               wl_npools(j,k,i) = valid_or (decomp_npools_vr(j,k,i), 0._r8)
            ENDDO
         ENDDO
         DO j = 1, nl_soil
            ! upstream IniTimeVar: 5 + 5 gN/m3 mineral nitrogen
            wl_smin_nh4(j,i) = valid_or (smin_nh4_vr(j,i), 5._r8)
            wl_smin_no3(j,i) = valid_or (smin_no3_vr(j,i), 5._r8)
            wl_sminn   (j,i) = valid_or (sminn_vr   (j,i), 10._r8)
         ENDDO

         DO j = 1, nl_soil
            existing_c = sum(wl_cpools(j,:,i))
            IF (existing_c > 1.e-12_r8) CYCLE
            IF (.not. allocated(OM_density)) CYCLE
            IF (.not. (OM_density(j,i) > 0._r8 .and. OM_density(j,i) < 1.e30_r8)) CYCLE
            tot_c = OM_density(j,i) * carbon_per_kg_om
            wl_cpools(j,i_met_lit,i) = 0.05_r8 * tot_c
            wl_cpools(j,i_cel_lit,i) = 0.10_r8 * tot_c
            wl_cpools(j,i_lig_lit,i) = 0.05_r8 * tot_c
            wl_cpools(j,i_cwd    ,i) = 0._r8
            wl_cpools(j,i_soil1  ,i) = 0.05_r8 * tot_c
            wl_cpools(j,i_soil2  ,i) = 0.25_r8 * tot_c
            wl_cpools(j,i_soil3  ,i) = 0.50_r8 * tot_c
            wl_npools(j,i_met_lit,i) = 0.05_r8 * tot_c / wetland_cn_ratio
            wl_npools(j,i_cel_lit,i) = 0.10_r8 * tot_c / wetland_cn_ratio
            wl_npools(j,i_lig_lit,i) = 0.05_r8 * tot_c / wetland_cn_ratio
            wl_npools(j,i_cwd    ,i) = 0._r8
            wl_npools(j,i_soil1  ,i) = 0.05_r8 * tot_c / wetland_cn_ratio
            wl_npools(j,i_soil2  ,i) = 0.25_r8 * tot_c / wetland_cn_ratio
            wl_npools(j,i_soil3  ,i) = 0.50_r8 * tot_c / wetland_cn_ratio
         ENDDO
      ENDDO

   CONTAINS

      real(r8) FUNCTION valid_or (x, fallback)
         real(r8), intent(in) :: x, fallback
         IF (x >= 0._r8 .and. x < 1.e30_r8) THEN
            valid_or = x
         ELSE
            valid_or = fallback
         ENDIF
      END FUNCTION valid_or

   END SUBROUTINE methane_wetland_bgc_cold_init

   ! -------
   ! Decomposition cascade of wetland patch ipatch on its shadow pools.
   SUBROUTINE methane_wetland_bgc_step (ipatch, deltim)

      USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
      IMPLICIT NONE
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: deltim

      IF (.not. ieee_is_finite(deltim) .or. deltim <= 0._r8) THEN
         CALL CoLM_stop(' ***** ERROR: wetland CH4/BGC coupling requires a finite positive timestep')
      ENDIF
      IF (.not. allocated(wl_cpools)) &
         CALL CoLM_stop(' ***** ERROR: wetland CH4/BGC shadow pools are not allocated')

      CALL save_host_column (ipatch)

      ! Same clean per-patch flux state as CoLM202X (and the full BGC driver).
      decomp_cpools_sourcesink  (:,:,ipatch) = 0._r8
      decomp_npools_sourcesink  (:,:,ipatch) = 0._r8
      decomp_hr_vr              (:,:,ipatch) = 0._r8
      decomp_ctransfer_vr       (:,:,ipatch) = 0._r8
      decomp_ntransfer_vr       (:,:,ipatch) = 0._r8
      decomp_sminn_flux_vr      (:,:,ipatch) = 0._r8
      sminn_to_denit_decomp_vr  (:,:,ipatch) = 0._r8
      pmnf_decomp               (:,:,ipatch) = 0._r8
      p_decomp_cpool_loss       (:,:,ipatch) = 0._r8
      net_nmin_vr               (:,ipatch)   = 0._r8
      gross_nmin_vr             (:,ipatch)   = 0._r8
      potential_immob_vr        (:,ipatch)   = 0._r8
      phr_vr                    (:,ipatch)   = 0._r8
      pot_f_nit_vr              (:,ipatch)   = 0._r8
      o_scalar                  (:,ipatch)   = 1._r8
      fpi_vr                    (:,ipatch)   = 1._r8
      net_nmin                  (ipatch)     = 0._r8
      gross_nmin                (ipatch)     = 0._r8
      decomp_hr                 (ipatch)     = 0._r8
      somc_fire                 (ipatch)     = 0._r8
      som_c_leached             (ipatch)     = 0._r8
      som_n_leached             (ipatch)     = 0._r8
      denit                     (ipatch)     = 0._r8
      f_n2o_nit                 (ipatch)     = 0._r8
      smin_no3_leached          (ipatch)     = 0._r8
      smin_no3_runoff           (ipatch)     = 0._r8
      sminn_leached             (ipatch)     = 0._r8
      sminn_to_plant            (ipatch)     = 0._r8

      decomp_cpools_vr (:,:,ipatch) = wl_cpools   (:,:,ipatch)
      decomp_npools_vr (:,:,ipatch) = wl_npools   (:,:,ipatch)
      sminn_vr         (:,ipatch)   = wl_sminn    (:,ipatch)
      smin_nh4_vr      (:,ipatch)   = wl_smin_nh4 (:,ipatch)
      smin_no3_vr      (:,ipatch)   = wl_smin_no3 (:,ipatch)

      CALL decomp_rate_constants_bgc (ipatch, nl_soil, z_soi)
      CALL SoilBiogeochemPotential   (ipatch, nl_soil, ndecomp_pools, ndecomp_transitions)
      CALL competition_no_plant      (ipatch, deltim)
      CALL SoilBiogeochemDecomp      (ipatch, nl_soil, ndecomp_pools, ndecomp_transitions, dz_soi)

      wv_hr_vr        (:,:) = decomp_hr_vr        (1:nl_soil,:,ipatch)
      wv_o_scalar     (:)   = o_scalar            (1:nl_soil,ipatch)
      wv_pot_f_nit_vr (:)   = pot_f_nit_vr        (1:nl_soil,ipatch)
      pd_ctransfer    (:,:) = decomp_ctransfer_vr (1:nl_soil,:,ipatch)
      pd_ntransfer    (:,:) = decomp_ntransfer_vr (1:nl_soil,:,ipatch)
      pd_sminn_flux   (:,:) = decomp_sminn_flux_vr(1:nl_soil,:,ipatch)
      pd_sminn_denit  (:,:) = sminn_to_denit_decomp_vr(1:nl_soil,:,ipatch)
      pd_gross_nmin   (:)   = gross_nmin_vr       (1:nl_soil,ipatch)
      pd_immob_nh4    (:)   = actual_immob_nh4_vr (1:nl_soil,ipatch)
      pd_immob_no3    (:)   = actual_immob_no3_vr (1:nl_soil,ipatch)
      view_patch = ipatch
      pending    = .true.

      CALL restore_host_column (ipatch)

   END SUBROUTINE methane_wetland_bgc_step

   logical FUNCTION methane_wetland_bgc_view_active (ipatch)

      IMPLICIT NONE
      integer, intent(in) :: ipatch

      methane_wetland_bgc_view_active = allocated(wl_cpools) .and. (view_patch == ipatch)

   END FUNCTION methane_wetland_bgc_view_active

   ! -------
   ! Pool update of the wetland shadow pools with the fluxes of this step's
   ! cascade.  CoLM202X CDecompStateUpdate (apply_direct) and
   ! SoilBiogeochemNDecompStateUpdate (apply_direct), on a source/sink
   ! cleared at the start of the step, in the same order.
   SUBROUTINE methane_wetland_bgc_apply (ipatch, deltim)

      USE MOD_Namelist, only: DEF_USE_NITRIF
      IMPLICIT NONE
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: deltim

      real(r8) :: css (nl_soil,ndecomp_pools), nss (nl_soil,ndecomp_pools)
      integer  :: j, k

      IF (.not. pending .or. view_patch /= ipatch) RETURN

      css(:,:) = 0._r8
      DO k = 1, ndecomp_transitions
         DO j = 1, nl_soil
            css(j,donor_pool(k)) = css(j,donor_pool(k)) - (wv_hr_vr(j,k) + pd_ctransfer(j,k)) * deltim
         ENDDO
      ENDDO
      DO k = 1, ndecomp_transitions
         IF (receiver_pool(k) /= 0) THEN
            DO j = 1, nl_soil
               css(j,receiver_pool(k)) = css(j,receiver_pool(k)) + pd_ctransfer(j,k) * deltim
            ENDDO
         ENDIF
      ENDDO

      nss(:,:) = 0._r8
      DO k = 1, ndecomp_transitions
         DO j = 1, nl_soil
            nss(j,donor_pool(k)) = nss(j,donor_pool(k)) - pd_ntransfer(j,k) * deltim
         ENDDO
      ENDDO
      DO k = 1, ndecomp_transitions
         IF (receiver_pool(k) /= 0) THEN
            DO j = 1, nl_soil
               nss(j,receiver_pool(k)) = nss(j,receiver_pool(k)) + &
                  (pd_ntransfer(j,k) + pd_sminn_flux(j,k)) * deltim
            ENDDO
         ELSE
            DO j = 1, nl_soil
               nss(j,donor_pool(k)) = nss(j,donor_pool(k)) - pd_sminn_flux(j,k) * deltim
            ENDDO
         ENDIF
      ENDDO

      IF (.not. DEF_USE_NITRIF) THEN
         DO k = 1, ndecomp_transitions
            IF (receiver_pool(k) /= 0) THEN
               DO j = 1, nl_soil
                  wl_sminn(j,ipatch) = wl_sminn(j,ipatch) - &
                     (pd_sminn_denit(j,k) + pd_sminn_flux(j,k)) * deltim
               ENDDO
            ELSE
               DO j = 1, nl_soil
                  wl_sminn(j,ipatch) = wl_sminn(j,ipatch) + &
                     (pd_sminn_flux(j,k) - pd_sminn_denit(j,k)) * deltim
               ENDDO
            ENDIF
         ENDDO
      ELSE
         DO j = 1, nl_soil
            wl_smin_nh4(j,ipatch) = wl_smin_nh4(j,ipatch) + &
               (pd_gross_nmin(j) - pd_immob_nh4(j)) * deltim
            wl_smin_no3(j,ipatch) = wl_smin_no3(j,ipatch) - pd_immob_no3(j) * deltim
            wl_sminn(j,ipatch) = wl_smin_nh4(j,ipatch) + wl_smin_no3(j,ipatch)
         ENDDO
      ENDIF

      wl_cpools(1:nl_soil,:,ipatch) = wl_cpools(1:nl_soil,:,ipatch) + css(:,:)
      wl_npools(1:nl_soil,:,ipatch) = wl_npools(1:nl_soil,:,ipatch) + nss(:,:)

      pending = .false.

   END SUBROUTINE methane_wetland_bgc_apply

   ! -------
   ! CoLM202X SoilBiogeochemCompetitionNoPlant: mineral N competition of a
   ! column without plants (immobilisation only).  Runs on the host column
   ! while the shadow pools are swapped in.
   SUBROUTINE competition_no_plant (i, deltim)

      USE MOD_Namelist, only: DEF_USE_NITRIF
      IMPLICIT NONE
      integer,  intent(in) :: i
      real(r8), intent(in) :: deltim

      integer  :: j
      real(r8) :: available, remaining, actual_immob, potential_immob

      actual_immob = 0._r8
      potential_immob = 0._r8
      DO j = 1, nl_soil
         sminn_to_plant_vr(j,i) = 0._r8
         smin_nh4_to_plant_vr(j,i) = 0._r8
         smin_no3_to_plant_vr(j,i) = 0._r8
         supplement_to_sminn_vr(j,i) = 0._r8
         sminn_to_denit_excess_vr(j,i) = 0._r8
         f_nit_vr(j,i) = 0._r8
         f_denit_vr(j,i) = 0._r8

         IF (.not. DEF_USE_NITRIF) THEN
            available = max(sminn_vr(j,i), 0._r8) / deltim
            actual_immob_vr(j,i) = min(max(potential_immob_vr(j,i), 0._r8), available)
            actual_immob_nh4_vr(j,i) = 0._r8
            actual_immob_no3_vr(j,i) = 0._r8
         ELSE
            actual_immob_nh4_vr(j,i) = min(max(potential_immob_vr(j,i), 0._r8), &
               max(smin_nh4_vr(j,i), 0._r8) / deltim)
            remaining = max(potential_immob_vr(j,i) - actual_immob_nh4_vr(j,i), 0._r8)
            actual_immob_no3_vr(j,i) = min(remaining, max(smin_no3_vr(j,i), 0._r8) / deltim)
            actual_immob_vr(j,i) = actual_immob_nh4_vr(j,i) + actual_immob_no3_vr(j,i)
         ENDIF

         IF (potential_immob_vr(j,i) > 0._r8) THEN
            fpi_vr(j,i) = actual_immob_vr(j,i) / potential_immob_vr(j,i)
         ELSE
            fpi_vr(j,i) = 1._r8
         ENDIF
         actual_immob = actual_immob + actual_immob_vr(j,i) * dz_soi(j)
         potential_immob = potential_immob + max(potential_immob_vr(j,i), 0._r8) * dz_soi(j)
      ENDDO

      sminn_to_plant(i) = 0._r8
      fpg(i) = 1._r8
      IF (potential_immob > 0._r8) THEN
         fpi(i) = actual_immob / potential_immob
      ELSE
         fpi(i) = 1._r8
      ENDIF

   END SUBROUTINE competition_no_plant

   ! -------
   ! Every host BGC array the cascade above writes, column i.
   SUBROUTINE save_host_column (i)

      IMPLICIT NONE
      integer, intent(in) :: i

      s3_cpools      = decomp_cpools_vr        (:,:,i)
      s3_npools      = decomp_npools_vr        (:,:,i)
      s3_decomp_k    = decomp_k                (:,:,i)
      s3_cn_pools    = cn_decomp_pools         (:,:,i)
      s3_p_loss      = p_decomp_cpool_loss     (:,:,i)
      s3_pmnf        = pmnf_decomp             (:,:,i)
      s3_ctransfer   = decomp_ctransfer_vr     (:,:,i)
      s3_hr          = decomp_hr_vr            (:,:,i)
      s3_ntransfer   = decomp_ntransfer_vr     (:,:,i)
      s3_sminn_flux  = decomp_sminn_flux_vr    (:,:,i)
      s3_sminn_denit = sminn_to_denit_decomp_vr(:,:,i)
      s3_css         = decomp_cpools_sourcesink(:,:,i)
      s3_nss         = decomp_npools_sourcesink(:,:,i)

      s2_sminn              = sminn_vr                (:,i)
      s2_nh4                = smin_nh4_vr             (:,i)
      s2_no3                = smin_no3_vr             (:,i)
      s2_t_scalar           = t_scalar                (:,i)
      s2_w_scalar           = w_scalar                (:,i)
      s2_o_scalar           = o_scalar                (:,i)
      s2_depth_scalar       = depth_scalar            (:,i)
      s2_fpi_vr             = fpi_vr                  (:,i)
      s2_gross_nmin_vr      = gross_nmin_vr           (:,i)
      s2_phr_vr             = phr_vr                  (:,i)
      s2_potential_immob_vr = potential_immob_vr      (:,i)
      s2_net_nmin_vr        = net_nmin_vr             (:,i)
      s2_pot_f_nit_vr       = pot_f_nit_vr            (:,i)
      s2_sminn_to_plant_vr  = sminn_to_plant_vr       (:,i)
      s2_nh4_to_plant_vr    = smin_nh4_to_plant_vr    (:,i)
      s2_no3_to_plant_vr    = smin_no3_to_plant_vr    (:,i)
      s2_supplement         = supplement_to_sminn_vr  (:,i)
      s2_denit_excess       = sminn_to_denit_excess_vr(:,i)
      s2_f_nit_vr           = f_nit_vr                (:,i)
      s2_f_denit_vr         = f_denit_vr              (:,i)
      s2_actual_immob_vr    = actual_immob_vr         (:,i)
      s2_immob_nh4          = actual_immob_nh4_vr     (:,i)
      s2_immob_no3          = actual_immob_no3_vr     (:,i)

      s1_fpi            = fpi              (i)
      s1_fpg            = fpg              (i)
      s1_gross_nmin     = gross_nmin       (i)
      s1_net_nmin       = net_nmin         (i)
      s1_decomp_hr      = decomp_hr        (i)
      s1_somc_fire      = somc_fire        (i)
      s1_som_c_leached  = som_c_leached    (i)
      s1_som_n_leached  = som_n_leached    (i)
      s1_denit          = denit            (i)
      s1_f_n2o_nit      = f_n2o_nit        (i)
      s1_no3_leached    = smin_no3_leached (i)
      s1_no3_runoff     = smin_no3_runoff  (i)
      s1_sminn_leached  = sminn_leached    (i)
      s1_sminn_to_plant = sminn_to_plant   (i)

   END SUBROUTINE save_host_column

   SUBROUTINE restore_host_column (i)

      IMPLICIT NONE
      integer, intent(in) :: i

      decomp_cpools_vr        (:,:,i) = s3_cpools
      decomp_npools_vr        (:,:,i) = s3_npools
      decomp_k                (:,:,i) = s3_decomp_k
      cn_decomp_pools         (:,:,i) = s3_cn_pools
      p_decomp_cpool_loss     (:,:,i) = s3_p_loss
      pmnf_decomp             (:,:,i) = s3_pmnf
      decomp_ctransfer_vr     (:,:,i) = s3_ctransfer
      decomp_hr_vr            (:,:,i) = s3_hr
      decomp_ntransfer_vr     (:,:,i) = s3_ntransfer
      decomp_sminn_flux_vr    (:,:,i) = s3_sminn_flux
      sminn_to_denit_decomp_vr(:,:,i) = s3_sminn_denit
      decomp_cpools_sourcesink(:,:,i) = s3_css
      decomp_npools_sourcesink(:,:,i) = s3_nss

      sminn_vr                (:,i) = s2_sminn
      smin_nh4_vr             (:,i) = s2_nh4
      smin_no3_vr             (:,i) = s2_no3
      t_scalar                (:,i) = s2_t_scalar
      w_scalar                (:,i) = s2_w_scalar
      o_scalar                (:,i) = s2_o_scalar
      depth_scalar            (:,i) = s2_depth_scalar
      fpi_vr                  (:,i) = s2_fpi_vr
      gross_nmin_vr           (:,i) = s2_gross_nmin_vr
      phr_vr                  (:,i) = s2_phr_vr
      potential_immob_vr      (:,i) = s2_potential_immob_vr
      net_nmin_vr             (:,i) = s2_net_nmin_vr
      pot_f_nit_vr            (:,i) = s2_pot_f_nit_vr
      sminn_to_plant_vr       (:,i) = s2_sminn_to_plant_vr
      smin_nh4_to_plant_vr    (:,i) = s2_nh4_to_plant_vr
      smin_no3_to_plant_vr    (:,i) = s2_no3_to_plant_vr
      supplement_to_sminn_vr  (:,i) = s2_supplement
      sminn_to_denit_excess_vr(:,i) = s2_denit_excess
      f_nit_vr                (:,i) = s2_f_nit_vr
      f_denit_vr              (:,i) = s2_f_denit_vr
      actual_immob_vr         (:,i) = s2_actual_immob_vr
      actual_immob_nh4_vr     (:,i) = s2_immob_nh4
      actual_immob_no3_vr     (:,i) = s2_immob_no3

      fpi              (i) = s1_fpi
      fpg              (i) = s1_fpg
      gross_nmin       (i) = s1_gross_nmin
      net_nmin         (i) = s1_net_nmin
      decomp_hr        (i) = s1_decomp_hr
      somc_fire        (i) = s1_somc_fire
      som_c_leached    (i) = s1_som_c_leached
      som_n_leached    (i) = s1_som_n_leached
      denit            (i) = s1_denit
      f_n2o_nit        (i) = s1_f_n2o_nit
      smin_no3_leached (i) = s1_no3_leached
      smin_no3_runoff  (i) = s1_no3_runoff
      sminn_leached    (i) = s1_sminn_leached
      sminn_to_plant   (i) = s1_sminn_to_plant

   END SUBROUTINE restore_host_column

   ! -------
   SUBROUTINE methane_wetland_bgc_write_restart (file_restart, compress)

      USE MOD_LandPatch, only: landpatch
      USE MOD_NetCDFVector, only: ncio_write_vector, ncio_define_dimension_vector
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer, intent(in) :: compress

      CALL ncio_define_dimension_vector (file_restart, landpatch, 'patch')
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'soil_full', nl_soil_full)
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'soil', nl_soil)
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'ndecomp_pools', ndecomp_pools)
      CALL ncio_write_vector (file_restart, 'ch4_wetland_decomp_cpools_vr', 'soil_full', nl_soil_full, &
         'ndecomp_pools', ndecomp_pools, 'patch', landpatch, wl_cpools, compress)
      CALL ncio_write_vector (file_restart, 'ch4_wetland_decomp_npools_vr', 'soil_full', nl_soil_full, &
         'ndecomp_pools', ndecomp_pools, 'patch', landpatch, wl_npools, compress)
      CALL ncio_write_vector (file_restart, 'ch4_wetland_sminn_vr', 'soil', nl_soil, &
         'patch', landpatch, wl_sminn, compress)
      CALL ncio_write_vector (file_restart, 'ch4_wetland_smin_nh4_vr', 'soil', nl_soil, &
         'patch', landpatch, wl_smin_nh4, compress)
      CALL ncio_write_vector (file_restart, 'ch4_wetland_smin_no3_vr', 'soil', nl_soil, &
         'patch', landpatch, wl_smin_no3, compress)

   END SUBROUTINE methane_wetland_bgc_write_restart

   ! Returns present = .false. when the restart has no wetland pools (an
   ! older restart or a cold start); the caller then runs the cold start.
   SUBROUTINE methane_wetland_bgc_read_restart (file_restart, present_out)

      USE MOD_LandPatch, only: landpatch
      USE MOD_NetCDFVector, only: ncio_read_vector
      USE MOD_MT_NcIO, only: ncio_vector_group_presence
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      logical, intent(out) :: present_out

      character(len=40), parameter :: fields(5) = [character(len=40) :: &
         'ch4_wetland_decomp_cpools_vr', 'ch4_wetland_decomp_npools_vr', &
         'ch4_wetland_sminn_vr', 'ch4_wetland_smin_nh4_vr', 'ch4_wetland_smin_no3_vr']
      logical :: flags(5)

      CALL ncio_vector_group_presence (file_restart, fields, landpatch, flags)
      present_out = all(flags)
      IF (any(flags) .and. .not. present_out) &
         CALL CoLM_stop (' ***** ERROR: methane restart has only some of the wetland BGC pools')
      IF (.not. present_out) RETURN

      CALL ncio_read_vector (file_restart, 'ch4_wetland_decomp_cpools_vr', nl_soil_full, ndecomp_pools, &
         landpatch, wl_cpools)
      CALL ncio_read_vector (file_restart, 'ch4_wetland_decomp_npools_vr', nl_soil_full, ndecomp_pools, &
         landpatch, wl_npools)
      CALL ncio_read_vector (file_restart, 'ch4_wetland_sminn_vr',    nl_soil, landpatch, wl_sminn)
      CALL ncio_read_vector (file_restart, 'ch4_wetland_smin_nh4_vr', nl_soil, landpatch, wl_smin_nh4)
      CALL ncio_read_vector (file_restart, 'ch4_wetland_smin_no3_vr', nl_soil, landpatch, wl_smin_no3)

   END SUBROUTINE methane_wetland_bgc_read_restart

   ! -------
   ! LULCC: keep the shadow pools of the old patch layout.
   SUBROUTINE methane_wetland_bgc_save_lulcc ()

      IMPLICIT NONE

      IF (allocated(lulcc_cpools  )) deallocate (lulcc_cpools  )
      IF (allocated(lulcc_npools  )) deallocate (lulcc_npools  )
      IF (allocated(lulcc_sminn   )) deallocate (lulcc_sminn   )
      IF (allocated(lulcc_smin_nh4)) deallocate (lulcc_smin_nh4)
      IF (allocated(lulcc_smin_no3)) deallocate (lulcc_smin_no3)
      IF (.not. allocated(wl_cpools)) RETURN
      allocate (lulcc_cpools  , source = wl_cpools  )
      allocate (lulcc_npools  , source = wl_npools  )
      allocate (lulcc_sminn   , source = wl_sminn   )
      allocate (lulcc_smin_nh4, source = wl_smin_nh4)
      allocate (lulcc_smin_no3, source = wl_smin_no3)

   END SUBROUTINE methane_wetland_bgc_save_lulcc

   ! LULCC: shadow pools of the new patch layout.  A wetland patch that was a
   ! wetland patch of the same element before (same element index and patch
   ! class) keeps its pools; any other wetland patch gets the cold start.
   ! The area-weighted remap of the land-cover change itself is P6 work.
   SUBROUTINE methane_wetland_bgc_remap_lulcc (patchtype, patchclass_new, eindex_new, &
         patchclass_old, eindex_old)

      IMPLICIT NONE
      integer,   intent(in) :: patchtype(:)
      integer,   intent(in) :: patchclass_new(:), patchclass_old(:)
      integer*8, intent(in) :: eindex_new(:), eindex_old(:)

      integer :: nnew, np, op

      nnew = size(patchclass_new)
      CALL methane_wetland_bgc_allocate (nnew)
      CALL methane_wetland_bgc_cold_init (patchtype, nnew)
      IF (.not. allocated(lulcc_cpools)) RETURN

      DO np = 1, nnew
         IF (patchtype(np) /= 2) CYCLE
         DO op = 1, min(size(patchclass_old), size(lulcc_cpools,3))
            IF (eindex_old(op) /= eindex_new(np)) CYCLE
            IF (patchclass_old(op) /= patchclass_new(np)) CYCLE
            wl_cpools  (:,:,np) = lulcc_cpools  (:,:,op)
            wl_npools  (:,:,np) = lulcc_npools  (:,:,op)
            wl_sminn   (:,np)   = lulcc_sminn   (:,op)
            wl_smin_nh4(:,np)   = lulcc_smin_nh4(:,op)
            wl_smin_no3(:,np)   = lulcc_smin_no3(:,op)
            EXIT
         ENDDO
      ENDDO

      deallocate (lulcc_cpools, lulcc_npools, lulcc_sminn, lulcc_smin_nh4, lulcc_smin_no3)

   END SUBROUTINE methane_wetland_bgc_remap_lulcc

END MODULE MOD_MT_Methane_WetlandBGC
#endif
