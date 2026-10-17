#include <define.h>
#ifdef BGC
MODULE MOD_BGC_CNCStateUpdate3

!-------------------------------------------------------------------------------------------------------
! !DESCRIPTION
! First updates in vegetation and soil carbon. The major updates are included in bgc_CNCStateUpdate1Mod
!  1. Update fire-associated veg and soil(litter) C pool size changes
!  2. Record the accumulated C transfers associated to fire for semi-analytic spinup

! !ORIGINAL:
! The Community Land Model version 5.0 (CLM5.0)

! !REVISION:
! Xingjie Lu, 2022, 1) modify original CLM5 to be compatible with CoLM code structure.
!                   2) Record accumulated fire-associated C transfers for veg and soil C semi-analytic spinup

   USE MOD_Precision
   USE MOD_BGC_Vars_TimeInvariants, only: &
            i_met_lit,i_cel_lit,i_lig_lit ,i_cwd
   USE MOD_BGC_Vars_TimeVariables, only: &
     ! decomposition pools & fluxes variables (inout)
            decomp_cpools_vr

   USE MOD_BGC_Vars_1DFluxes, only: &
            m_decomp_cpools_to_fire_vr, &
            fire_mortality_to_met_c, fire_mortality_to_cel_c, &
            fire_mortality_to_lig_c, fire_mortality_to_cwdc

   USE MOD_BGC_Vars_PFTimeVariables, only: &
     ! vegetation carbon state variables (inout)
            leafc_p            , leafc_storage_p     , leafc_xfer_p     , &
            frootc_p           , frootc_storage_p    , frootc_xfer_p    , &
            livestemc_p        , livestemc_storage_p , livestemc_xfer_p , &
            deadstemc_p        , deadstemc_storage_p , deadstemc_xfer_p , &
            livecrootc_p       , livecrootc_storage_p, livecrootc_xfer_p, &
            deadcrootc_p       , deadcrootc_storage_p, deadcrootc_xfer_p, &
            gresp_storage_p    , gresp_xfer_p

   USE MOD_BGC_Vars_1DPFTFluxes, only: &
     ! vegetation carbon flux variables
            m_leafc_to_fire_p        , m_leafc_storage_to_fire_p     , m_leafc_xfer_to_fire_p     , &
            m_frootc_to_fire_p       , m_frootc_storage_to_fire_p    , m_frootc_xfer_to_fire_p    , &
            m_livestemc_to_fire_p    , m_livestemc_storage_to_fire_p , m_livestemc_xfer_to_fire_p , &
            m_deadstemc_to_fire_p    , m_deadstemc_storage_to_fire_p , m_deadstemc_xfer_to_fire_p , &
            m_livecrootc_to_fire_p   , m_livecrootc_storage_to_fire_p, m_livecrootc_xfer_to_fire_p, &
            m_deadcrootc_to_fire_p   , m_deadcrootc_storage_to_fire_p, m_deadcrootc_xfer_to_fire_p, &
            m_livestemc_to_deadstemc_fire_p , m_livecrootc_to_deadcrootc_fire_p , &
            m_gresp_storage_to_fire_p       , m_gresp_xfer_to_fire_p            , &

            m_leafc_to_litter_fire_p        , m_leafc_storage_to_litter_fire_p     , m_leafc_xfer_to_litter_fire_p     , &
            m_frootc_to_litter_fire_p       , m_frootc_storage_to_litter_fire_p    , m_frootc_xfer_to_litter_fire_p    , &
            m_livestemc_to_litter_fire_p    , m_livestemc_storage_to_litter_fire_p , m_livestemc_xfer_to_litter_fire_p , &
            m_deadstemc_to_litter_fire_p    , m_deadstemc_storage_to_litter_fire_p , m_deadstemc_xfer_to_litter_fire_p , &
            m_livecrootc_to_litter_fire_p   , m_livecrootc_storage_to_litter_fire_p, m_livecrootc_xfer_to_litter_fire_p, &
            m_deadcrootc_to_litter_fire_p   , m_deadcrootc_storage_to_litter_fire_p, m_deadcrootc_xfer_to_litter_fire_p, &
            m_gresp_storage_to_litter_fire_p, m_gresp_xfer_to_litter_fire_p

   USE MOD_Namelist, only: DEF_USE_SASU, DEF_USE_DiagMatrix

   USE MOD_BGC_Vars_PFTimeVariables, only: &
             AKX_leafc_exit_p_acc      , AKX_leafc_st_exit_p_acc     , AKX_leafc_xf_exit_p_acc     , &
             AKX_frootc_exit_p_acc     , AKX_frootc_st_exit_p_acc    , AKX_frootc_xf_exit_p_acc    , &
             AKX_livestemc_exit_p_acc  , AKX_livestemc_st_exit_p_acc , AKX_livestemc_xf_exit_p_acc , &
             AKX_deadstemc_exit_p_acc  , AKX_deadstemc_st_exit_p_acc , AKX_deadstemc_xf_exit_p_acc , &
             AKX_livecrootc_exit_p_acc , AKX_livecrootc_st_exit_p_acc, AKX_livecrootc_xf_exit_p_acc, &
             AKX_deadcrootc_exit_p_acc , AKX_deadcrootc_st_exit_p_acc, AKX_deadcrootc_xf_exit_p_acc, &
             AKX_livestemc_to_deadstemc_p_acc, AKX_livecrootc_to_deadcrootc_p_acc

   USE MOD_BGC_Vars_TimeVariables, only: &
             I_met_c_vr_acc, I_cel_c_vr_acc, I_lig_c_vr_acc, I_cwd_c_vr_acc, &
             AKX_met_exit_c_vr_acc, AKX_cel_exit_c_vr_acc, AKX_lig_exit_c_vr_acc, AKX_cwd_exit_c_vr_acc

   IMPLICIT NONE

   PUBLIC CStateUpdate3

CONTAINS

   SUBROUTINE CStateUpdate3(i, ps, pe, deltim, nl_soil, ndecomp_pools)

   integer ,intent(in) :: i             ! patch index
   integer ,intent(in) :: ps            ! start pft index
   integer ,intent(in) :: pe            ! END pft index
   real(r8),intent(in) :: deltim        ! time step in seconds
   integer ,intent(in) :: nl_soil       ! number of total soil number
   integer ,intent(in) :: ndecomp_pools ! number total litter & soil pools

   integer j,l,m

      DO j = 1, nl_soil
         ! patch-level wood to column-level CWD (uncombusted wood)
         decomp_cpools_vr(j,i_cwd,i)     = decomp_cpools_vr(j,i_cwd,i)     &
                                         + fire_mortality_to_cwdc (j,i) * deltim

         ! patch-level wood to column-level litter (uncombusted wood)
         decomp_cpools_vr(j,i_met_lit,i) = decomp_cpools_vr(j,i_met_lit,i) &
                                         + fire_mortality_to_met_c(j,i) * deltim
         decomp_cpools_vr(j,i_cel_lit,i) = decomp_cpools_vr(j,i_cel_lit,i) &
                                         + fire_mortality_to_cel_c(j,i) * deltim
         decomp_cpools_vr(j,i_lig_lit,i) = decomp_cpools_vr(j,i_lig_lit,i) &
                                          + fire_mortality_to_lig_c(j,i) * deltim

         IF(DEF_USE_SASU .or. DEF_USE_DiagMatrix)THEN
            I_cwd_c_vr_acc(j,i) = I_cwd_c_vr_acc(j,i) + fire_mortality_to_cwdc(j,i) * deltim
            I_met_c_vr_acc(j,i) = I_met_c_vr_acc(j,i) + fire_mortality_to_met_c(j,i) * deltim
            I_cel_c_vr_acc(j,i) = I_cel_c_vr_acc(j,i) + fire_mortality_to_cel_c(j,i) * deltim
            I_lig_c_vr_acc(j,i) = I_lig_c_vr_acc(j,i) + fire_mortality_to_lig_c(j,i) * deltim
         ENDIF
      ENDDO

          ! litter and CWD losses to fire
       DO l = 1, ndecomp_pools
          DO j = 1, nl_soil
             decomp_cpools_vr(j,l,i) = decomp_cpools_vr(j,l,i) &
                                     - m_decomp_cpools_to_fire_vr(j,l,i) * deltim

             IF(DEF_USE_SASU .or. DEF_USE_DiagMatrix)THEN
                IF(l == i_met_lit) AKX_met_exit_c_vr_acc(j,i) = AKX_met_exit_c_vr_acc(j,i) + m_decomp_cpools_to_fire_vr(j,l,i) * deltim
                IF(l == i_cel_lit) AKX_cel_exit_c_vr_acc(j,i) = AKX_cel_exit_c_vr_acc(j,i) + m_decomp_cpools_to_fire_vr(j,l,i) * deltim
                IF(l == i_lig_lit) AKX_lig_exit_c_vr_acc(j,i) = AKX_lig_exit_c_vr_acc(j,i) + m_decomp_cpools_to_fire_vr(j,l,i) * deltim
                IF(l == i_cwd    ) AKX_cwd_exit_c_vr_acc(j,i) = AKX_cwd_exit_c_vr_acc(j,i) + m_decomp_cpools_to_fire_vr(j,l,i) * deltim
             ENDIF
          ENDDO
       ENDDO

         ! patch-level carbon fluxes from fire
      DO m = ps , pe
         gresp_storage_p     (m) = gresp_storage_p     (m) &
                                 - m_gresp_storage_to_fire_p            (m) * deltim
         gresp_storage_p     (m) = gresp_storage_p     (m) &
                                 - m_gresp_storage_to_litter_fire_p     (m) * deltim
         gresp_xfer_p        (m) = gresp_xfer_p        (m) &
                                 - m_gresp_xfer_to_fire_p               (m) * deltim
         gresp_xfer_p        (m) = gresp_xfer_p        (m) &
                                 - m_gresp_xfer_to_litter_fire_p        (m) * deltim
         ! displayed pools
         leafc_p             (m) = leafc_p             (m) &
                                 - m_leafc_to_fire_p                    (m) * deltim
         leafc_p             (m) = leafc_p             (m) &
                                 - m_leafc_to_litter_fire_p             (m) * deltim
         frootc_p            (m) = frootc_p            (m) &
                                 - m_frootc_to_fire_p                   (m) * deltim
         frootc_p            (m) = frootc_p            (m) &
                                 - m_frootc_to_litter_fire_p            (m) * deltim
         livestemc_p         (m) = livestemc_p         (m) &
                                 - m_livestemc_to_fire_p                (m) * deltim
         livestemc_p         (m) = livestemc_p         (m) &
                                 - m_livestemc_to_litter_fire_p         (m) * deltim &
                                 - m_livestemc_to_deadstemc_fire_p      (m) * deltim
         deadstemc_p         (m) = deadstemc_p         (m) &
                                 - m_deadstemc_to_fire_p                (m) * deltim
         deadstemc_p         (m) = deadstemc_p         (m) &
                                 - m_deadstemc_to_litter_fire_p         (m) * deltim &
                                 + m_livestemc_to_deadstemc_fire_p      (m) * deltim
         livecrootc_p        (m) = livecrootc_p        (m) &
                                 - m_livecrootc_to_fire_p               (m) * deltim
         livecrootc_p        (m) = livecrootc_p        (m) &
                                 - m_livecrootc_to_litter_fire_p        (m) * deltim &
                                 - m_livecrootc_to_deadcrootc_fire_p    (m) * deltim
         deadcrootc_p        (m) = deadcrootc_p        (m) &
                                 - m_deadcrootc_to_fire_p               (m) * deltim
         deadcrootc_p        (m) = deadcrootc_p        (m) &
                                 - m_deadcrootc_to_litter_fire_p        (m) * deltim &
                                 + m_livecrootc_to_deadcrootc_fire_p    (m) * deltim

       ! storage pools
         leafc_storage_p     (m) = leafc_storage_p     (m) &
                                 - m_leafc_storage_to_fire_p            (m) * deltim
         leafc_storage_p     (m) = leafc_storage_p     (m) &
                                 - m_leafc_storage_to_litter_fire_p     (m) * deltim
         frootc_storage_p    (m) = frootc_storage_p    (m) &
                                 - m_frootc_storage_to_fire_p           (m) * deltim
         frootc_storage_p    (m) = frootc_storage_p    (m) &
                                 - m_frootc_storage_to_litter_fire_p    (m) * deltim
         livestemc_storage_p (m) = livestemc_storage_p (m) &
                                 - m_livestemc_storage_to_fire_p        (m) * deltim
         livestemc_storage_p (m) = livestemc_storage_p (m) &
                                 - m_livestemc_storage_to_litter_fire_p (m) * deltim
         deadstemc_storage_p (m) = deadstemc_storage_p (m) &
                                 - m_deadstemc_storage_to_fire_p        (m) * deltim
         deadstemc_storage_p (m) = deadstemc_storage_p (m) &
                                 - m_deadstemc_storage_to_litter_fire_p (m) * deltim
         livecrootc_storage_p(m) = livecrootc_storage_p(m) &
                                 - m_livecrootc_storage_to_fire_p       (m) * deltim
         livecrootc_storage_p(m) = livecrootc_storage_p(m) &
                                 - m_livecrootc_storage_to_litter_fire_p(m) * deltim
         deadcrootc_storage_p(m) = deadcrootc_storage_p(m) &
                                 - m_deadcrootc_storage_to_fire_p       (m) * deltim
         deadcrootc_storage_p(m) = deadcrootc_storage_p(m) &
                                 - m_deadcrootc_storage_to_litter_fire_p(m) * deltim

      ! transfer pools
         leafc_xfer_p        (m) = leafc_xfer_p        (m) &
                                 - m_leafc_xfer_to_fire_p               (m) * deltim
         leafc_xfer_p        (m) = leafc_xfer_p        (m) &
                                 - m_leafc_xfer_to_litter_fire_p        (m) * deltim
         frootc_xfer_p       (m) = frootc_xfer_p       (m) &
                                 - m_frootc_xfer_to_fire_p              (m) * deltim
         frootc_xfer_p       (m) = frootc_xfer_p       (m) &
                                 - m_frootc_xfer_to_litter_fire_p       (m) * deltim
         livestemc_xfer_p    (m) = livestemc_xfer_p    (m) &
                                 - m_livestemc_xfer_to_fire_p           (m) * deltim
         livestemc_xfer_p    (m) = livestemc_xfer_p    (m) &
                                 - m_livestemc_xfer_to_litter_fire_p    (m) * deltim
         deadstemc_xfer_p    (m) = deadstemc_xfer_p    (m) &
                                 - m_deadstemc_xfer_to_fire_p           (m) * deltim
         deadstemc_xfer_p    (m) = deadstemc_xfer_p    (m) &
                                 - m_deadstemc_xfer_to_litter_fire_p    (m) * deltim
         livecrootc_xfer_p   (m) = livecrootc_xfer_p   (m) &
                                 - m_livecrootc_xfer_to_fire_p          (m) * deltim
         livecrootc_xfer_p   (m) = livecrootc_xfer_p   (m) &
                                 - m_livecrootc_xfer_to_litter_fire_p   (m) * deltim
         deadcrootc_xfer_p   (m) = deadcrootc_xfer_p   (m) &
                                 - m_deadcrootc_xfer_to_fire_p          (m) * deltim
          deadcrootc_xfer_p   (m) = deadcrootc_xfer_p   (m) &
                                  - m_deadcrootc_xfer_to_litter_fire_p   (m) * deltim

         IF(DEF_USE_SASU .or. DEF_USE_DiagMatrix)THEN
            AKX_leafc_exit_p_acc      (m) = AKX_leafc_exit_p_acc      (m) + (m_leafc_to_fire_p      (m) + m_leafc_to_litter_fire_p      (m)) * deltim
            AKX_leafc_st_exit_p_acc   (m) = AKX_leafc_st_exit_p_acc   (m) + (m_leafc_storage_to_fire_p   (m) + m_leafc_storage_to_litter_fire_p   (m)) * deltim
            AKX_leafc_xf_exit_p_acc   (m) = AKX_leafc_xf_exit_p_acc   (m) + (m_leafc_xfer_to_fire_p   (m) + m_leafc_xfer_to_litter_fire_p   (m)) * deltim
            AKX_frootc_exit_p_acc     (m) = AKX_frootc_exit_p_acc     (m) + (m_frootc_to_fire_p     (m) + m_frootc_to_litter_fire_p     (m)) * deltim
            AKX_frootc_st_exit_p_acc  (m) = AKX_frootc_st_exit_p_acc  (m) + (m_frootc_storage_to_fire_p  (m) + m_frootc_storage_to_litter_fire_p  (m)) * deltim
            AKX_frootc_xf_exit_p_acc  (m) = AKX_frootc_xf_exit_p_acc  (m) + (m_frootc_xfer_to_fire_p  (m) + m_frootc_xfer_to_litter_fire_p  (m)) * deltim
            AKX_livestemc_exit_p_acc  (m) = AKX_livestemc_exit_p_acc  (m) + (m_livestemc_to_fire_p  (m) + m_livestemc_to_litter_fire_p  (m) + m_livestemc_to_deadstemc_fire_p(m)) * deltim
            AKX_livestemc_st_exit_p_acc(m) = AKX_livestemc_st_exit_p_acc(m) + (m_livestemc_storage_to_fire_p(m) + m_livestemc_storage_to_litter_fire_p(m)) * deltim
            AKX_livestemc_xf_exit_p_acc(m) = AKX_livestemc_xf_exit_p_acc(m) + (m_livestemc_xfer_to_fire_p(m) + m_livestemc_xfer_to_litter_fire_p(m)) * deltim
            AKX_deadstemc_exit_p_acc  (m) = AKX_deadstemc_exit_p_acc  (m) + (m_deadstemc_to_fire_p  (m) + m_deadstemc_to_litter_fire_p  (m)) * deltim
            AKX_deadstemc_st_exit_p_acc(m) = AKX_deadstemc_st_exit_p_acc(m) + (m_deadstemc_storage_to_fire_p(m) + m_deadstemc_storage_to_litter_fire_p(m)) * deltim
            AKX_deadstemc_xf_exit_p_acc(m) = AKX_deadstemc_xf_exit_p_acc(m) + (m_deadstemc_xfer_to_fire_p(m) + m_deadstemc_xfer_to_litter_fire_p(m)) * deltim
            AKX_livecrootc_exit_p_acc (m) = AKX_livecrootc_exit_p_acc (m) + (m_livecrootc_to_fire_p (m) + m_livecrootc_to_litter_fire_p (m) + m_livecrootc_to_deadcrootc_fire_p(m)) * deltim
            AKX_livecrootc_st_exit_p_acc(m) = AKX_livecrootc_st_exit_p_acc(m) + (m_livecrootc_storage_to_fire_p(m) + m_livecrootc_storage_to_litter_fire_p(m)) * deltim
            AKX_livecrootc_xf_exit_p_acc(m) = AKX_livecrootc_xf_exit_p_acc(m) + (m_livecrootc_xfer_to_fire_p(m) + m_livecrootc_xfer_to_litter_fire_p(m)) * deltim
            AKX_deadcrootc_exit_p_acc (m) = AKX_deadcrootc_exit_p_acc (m) + (m_deadcrootc_to_fire_p (m) + m_deadcrootc_to_litter_fire_p (m)) * deltim
            AKX_deadcrootc_st_exit_p_acc(m) = AKX_deadcrootc_st_exit_p_acc(m) + (m_deadcrootc_storage_to_fire_p(m) + m_deadcrootc_storage_to_litter_fire_p(m)) * deltim
            AKX_deadcrootc_xf_exit_p_acc(m) = AKX_deadcrootc_xf_exit_p_acc(m) + (m_deadcrootc_xfer_to_fire_p(m) + m_deadcrootc_xfer_to_litter_fire_p(m)) * deltim
            AKX_livestemc_to_deadstemc_p_acc(m)   = AKX_livestemc_to_deadstemc_p_acc(m)   + m_livestemc_to_deadstemc_fire_p(m) * deltim
            AKX_livecrootc_to_deadcrootc_p_acc(m) = AKX_livecrootc_to_deadcrootc_p_acc(m) + m_livecrootc_to_deadcrootc_fire_p(m) * deltim
         ENDIF
      ENDDO

   END SUBROUTINE CStateUpdate3

END MODULE MOD_BGC_CNCStateUpdate3
#endif
