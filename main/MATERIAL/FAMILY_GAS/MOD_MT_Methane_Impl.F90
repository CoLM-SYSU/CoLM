#include <define.h>

#if (defined MaterialTransport) && (defined BGC)
MODULE MOD_MT_Methane_Impl
! METHANE_IMPL_GLUE_POINT: this module is the single glue/adaptor point
! between CoLM global patch state and the argument-driven methane_driver.
! Keep direct USEs of MOD_Vars_* concentrated here; lower methane physics
! should remain parameterized through explicit arguments.

   USE MOD_Precision
   USE MOD_Const_Physical, only: denh2o, denice
   USE MOD_Namelist,       only: DEF_USE_Dynamic_Lake
   USE MOD_Const_LC,       only: rootfr_lc => rootfr
   USE MOD_Vars_Global, only: maxsnl, nl_soil, nl_lake, PI, &
      z_soi, dz_soi, zi_soi
   USE MOD_Vars_TimeInvariants, only: patchtype, patchclass, patchlonr, patchlatr, &
      slpratio, fsatmax, fsatdcf, lakedepth, dz_lake, bsw, porsl
   USE MOD_Vars_TimeVariables, only: t_soisno, t_grnd, wliq_soisno, wice_soisno, &
      zwt, snowdp, wat, t_lake, lake_icefrac, wdsrf, wetwat, smp, lai, sai, rootr, ustar, fq
   USE MOD_Vars_1DForcing, only: forc_t, forc_pbot, forc_po2m, forc_pco2m, &
      forc_us, forc_vs
   USE MOD_Vars_1DFluxes, only: rsur, etr, frcsat
   USE MOD_MT_Methane_Driver,   only: methane_driver
   USE MOD_MT_Methane_State,    only: f_h2osfc, compute_f_h2osfc, &
      handle_methane_dry_lake_substep, reset_methane_inactive_lake_diagnostics, &
      accumulate_methane_lake_substep_diagnostics
   USE MOD_MT_Methane_Registry, only: igas_ch4
   USE MOD_MT_Methane_Const,    only: DEF_METHANE
   USE MOD_MT_Methane_BgcLink,  only: paddy_rice_fraction, PADDY_RICE_FRAC_MIN
   USE MOD_MT_Methane_Microbes, only: &
      reset_methane_inactive_lake_microbe_diagnostics, &
      accumulate_methane_lake_microbe_substep_diagnostics
   USE MOD_MT_Methane_WetlandBGC,  only: methane_wetland_bgc_step

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: ch4_impl_lake_step
   PUBLIC :: ch4_impl_wetland_decomp
   PUBLIC :: ch4_impl_soil_step

CONTAINS

   SUBROUTINE ch4_impl_lake_step (istep_local, ipatch, idate, deltim_phy, isub, nsub)

      IMPLICIT NONE
      integer,  intent(in) :: istep_local
      integer,  intent(in) :: ipatch
      integer,  intent(in) :: idate(3)
      real(r8), intent(in) :: deltim_phy
      integer,  intent(in) :: isub
      integer,  intent(in) :: nsub

      integer  :: lb, snl_loc
      real(r8) :: z_soisno_m(maxsnl+1:nl_soil)
      real(r8) :: dz_soisno_m(maxsnl+1:nl_soil)
      real(r8) :: zi_soisno_m(maxsnl:nl_soil)

      ! Lake methane follows WATERBODY substepping so transport and
      ! sediment-carbon tendencies use the same deltim as lake physics.
      IF (igas_ch4 <= 0) RETURN
      IF (patchtype(ipatch) /= 4) RETURN
      IF (.not. DEF_METHANE%allowlakeprod) THEN
         ! The lake pathway is inactive, but AccFlux still samples the shared
         ! diagnostic arrays.  Replace restart/previous-step process values
         ! with zero on every WATERBODY substep without changing inventories.
         CALL reset_methane_inactive_lake_diagnostics(ipatch)
         CALL reset_methane_inactive_lake_microbe_diagnostics(ipatch)
         CALL accumulate_methane_lake_microbe_substep_diagnostics( &
            ipatch, deltim_phy, isub, nsub)
         CALL accumulate_methane_lake_substep_diagnostics(ipatch, deltim_phy, isub, nsub)
         RETURN
      ENDIF

      IF (DEF_USE_Dynamic_Lake .and. &
         (wdsrf(ipatch) < 100._r8 .or. zwt(ipatch) > 0._r8)) THEN
         CALL handle_methane_dry_lake_substep(ipatch, deltim_phy)
         CALL reset_methane_inactive_lake_microbe_diagnostics(ipatch)
         CALL accumulate_methane_lake_microbe_substep_diagnostics( &
            ipatch, deltim_phy, isub, nsub)
         CALL accumulate_methane_lake_substep_diagnostics(ipatch, deltim_phy, isub, nsub)
         RETURN
      ENDIF

      CALL compute_f_h2osfc (ipatch, slpratio(ipatch), wdsrf(ipatch))
      CALL methane_soisno_geometry (ipatch, lb, snl_loc, z_soisno_m, dz_soisno_m, zi_soisno_m)

      CALL methane_driver(istep_local, ipatch, idate(1:3), patchclass(ipatch), patchtype(ipatch), &
         deltim_phy, lb, snl_loc, &
         patchlonr(ipatch)*180._r8/PI, patchlatr(ipatch)*180._r8/PI, &
         z_soisno_m(maxsnl+1:nl_soil), dz_soisno_m(maxsnl+1:nl_soil), &
         zi_soisno_m(maxsnl:nl_soil), t_soisno(maxsnl+1:nl_soil,ipatch), &
         t_grnd(ipatch), wliq_soisno(maxsnl+1:nl_soil,ipatch), &
         wice_soisno(maxsnl+1:nl_soil,ipatch), forc_t(ipatch), forc_pbot(ipatch), &
         forc_po2m(ipatch), forc_pco2m(ipatch), forc_us(ipatch), forc_vs(ipatch), &
         ustar(ipatch), fq(ipatch), &
         zwt(ipatch), rootfr_lc(1:nl_soil,patchclass(ipatch)), snowdp(ipatch), wat(ipatch), &
         rsur(ipatch), etr(ipatch), lakedepth(ipatch), dz_lake(1:nl_lake,ipatch), &
         t_lake(1:nl_lake,ipatch), lake_icefrac(1:nl_lake,ipatch), &
         wdsrf(ipatch), wetwat(ipatch), bsw(1:nl_soil,ipatch), smp(1:nl_soil,ipatch), &
         porsl(1:nl_soil,ipatch), lai(ipatch), sai(ipatch), rootr(1:nl_soil,ipatch), &
         fsatmax(ipatch), fsatdcf(ipatch), frcsat(ipatch), f_h2osfc(ipatch))

      CALL accumulate_methane_lake_microbe_substep_diagnostics( &
         ipatch, deltim_phy, isub, nsub)
      CALL accumulate_methane_lake_substep_diagnostics(ipatch, deltim_phy, isub, nsub)

   END SUBROUTINE ch4_impl_lake_step

   SUBROUTINE ch4_impl_wetland_decomp (ipatch, deltim)

      IMPLICIT NONE
      integer, intent(in) :: ipatch
      real(r8), intent(in) :: deltim

      ! Wetland patches run the soil-decomposition cascade on methane-owned
      ! shadow pools (MOD_MT_Methane_WetlandBGC) so methane can read a
      ! heterotrophic respiration without PFTs; host wetland arrays stay as
      ! upstream left them.
      IF (igas_ch4 <= 0) RETURN
      IF (patchtype(ipatch) /= 2) RETURN

      CALL methane_wetland_bgc_step (ipatch, deltim)

   END SUBROUTINE ch4_impl_wetland_decomp

   SUBROUTINE ch4_impl_soil_step (istep_local, ipatch, idate, deltim)

      IMPLICIT NONE
      integer,  intent(in) :: istep_local
      integer,  intent(in) :: ipatch
      integer,  intent(in) :: idate(3)
      real(r8), intent(in) :: deltim

      logical  :: run_methane
      logical  :: is_rice_paddy
      real(r8) :: rice_pft_frac
      integer  :: lb, snl_loc
      real(r8) :: z_soisno_m(maxsnl+1:nl_soil)
      real(r8) :: dz_soisno_m(maxsnl+1:nl_soil)
      real(r8) :: zi_soisno_m(maxsnl:nl_soil)

      ! Runtime activation: only execute methane_driver when CH4 is
      ! registered as a reactive material in the namelist. Lake CH4 is
      ! handled inside WATERBODY substeps; this path is soil/wetland CH4.
      run_methane = .false.
      is_rice_paddy = .false.
      rice_pft_frac = 0._r8

      IF (igas_ch4 > 0) THEN
         IF (DEF_METHANE%only_wetland) THEN
            IF (patchtype(ipatch) == 2) run_methane = .true.
         ELSE
            IF ((patchtype(ipatch) == 2) .or. (patchtype(ipatch) == 0)) run_methane = .true.
         ENDIF
#ifdef CROP
         ! Detect rice paddy via BgcLink helper.  Single source of truth
         ! shared with AccFlux / Hist filter and the methane() finundated blend.
         IF (DEF_METHANE%enable_rice_paddy .and. patchtype(ipatch) == 0) THEN
            rice_pft_frac = paddy_rice_fraction(ipatch)
            IF (rice_pft_frac > PADDY_RICE_FRAC_MIN) THEN
               is_rice_paddy = .true.
               run_methane   = .true.
            ENDIF
         ENDIF
#endif
      ENDIF

      IF (.not. run_methane) RETURN

      CALL compute_f_h2osfc (ipatch, slpratio(ipatch), wdsrf(ipatch))
      CALL methane_soisno_geometry (ipatch, lb, snl_loc, z_soisno_m, dz_soisno_m, zi_soisno_m)

      CALL methane_driver(istep_local, ipatch, idate(1:3), patchclass(ipatch), patchtype(ipatch), &
         deltim, lb, snl_loc, &
         patchlonr(ipatch)*180._r8/PI, patchlatr(ipatch)*180._r8/PI, &
         z_soisno_m(maxsnl+1:nl_soil), dz_soisno_m(maxsnl+1:nl_soil), &
         zi_soisno_m(maxsnl:nl_soil), t_soisno(maxsnl+1:nl_soil,ipatch), &
         t_grnd(ipatch), wliq_soisno(maxsnl+1:nl_soil,ipatch), &
         wice_soisno(maxsnl+1:nl_soil,ipatch), forc_t(ipatch), forc_pbot(ipatch), &
         forc_po2m(ipatch), forc_pco2m(ipatch), forc_us(ipatch), forc_vs(ipatch), &
         ustar(ipatch), fq(ipatch), &
         zwt(ipatch), rootfr_lc(1:nl_soil,patchclass(ipatch)), snowdp(ipatch), wat(ipatch), &
         rsur(ipatch), etr(ipatch), lakedepth(ipatch), dz_lake(1:nl_lake,ipatch), &
         t_lake(1:nl_lake,ipatch), lake_icefrac(1:nl_lake,ipatch), &
         wdsrf(ipatch), wetwat(ipatch), bsw(1:nl_soil,ipatch), smp(1:nl_soil,ipatch), &
         porsl(1:nl_soil,ipatch), lai(ipatch), sai(ipatch), rootr(1:nl_soil,ipatch), &
         fsatmax(ipatch), fsatdcf(ipatch), frcsat(ipatch), f_h2osfc(ipatch), &
         is_rice_paddy_in=is_rice_paddy, rice_pft_frac_in=rice_pft_frac)

   END SUBROUTINE ch4_impl_soil_step

   SUBROUTINE methane_soisno_geometry (ipatch, lb, snl_loc, z_soisno_m, dz_soisno_m, zi_soisno_m)

      IMPLICIT NONE
      integer,  intent(in)  :: ipatch
      integer,  intent(out) :: lb
      integer,  intent(out) :: snl_loc
      real(r8), intent(out) :: z_soisno_m(maxsnl+1:nl_soil)
      real(r8), intent(out) :: dz_soisno_m(maxsnl+1:nl_soil)
      real(r8), intent(out) :: zi_soisno_m(maxsnl:nl_soil)

      integer  :: j_meth
      real(r8) :: snow_mass_total_m
      real(r8) :: snow_layer_mass_m
      real(r8) :: snow_layer_depth_m

      ! CoLMMAIN keeps soil+snow geometry as local work arrays, while
      ! CoLMDRIVER only has persistent soil temperatures/water masses.
      ! Reconstruct the geometry needed by the methane snow/soil diffusion
      ! wrapper from persistent snow water state plus the canonical soil grid.
      z_soisno_m(maxsnl+1:0)  = 0._r8
      dz_soisno_m(maxsnl+1:0) = 0._r8
      zi_soisno_m(maxsnl:0)   = 0._r8
      z_soisno_m(1:nl_soil)   = z_soi(1:nl_soil)
      dz_soisno_m(1:nl_soil)  = dz_soi(1:nl_soil)
      zi_soisno_m(0)          = 0._r8
      zi_soisno_m(1:nl_soil)  = zi_soi(1:nl_soil)

      snl_loc = 0
      DO j_meth = 0, maxsnl + 1, -1
         IF (wliq_soisno(j_meth,ipatch) + wice_soisno(j_meth,ipatch) > 0._r8) THEN
            snl_loc = snl_loc - 1
         ELSE
            EXIT
         ENDIF
      ENDDO
      lb = snl_loc + 1

      IF (snl_loc < 0) THEN
         snow_mass_total_m = 0._r8
         DO j_meth = 0, lb, -1
            snow_mass_total_m = snow_mass_total_m + &
               max(wliq_soisno(j_meth,ipatch) + wice_soisno(j_meth,ipatch), 0._r8)
         ENDDO

         DO j_meth = 0, lb, -1
            snow_layer_mass_m = max(wliq_soisno(j_meth,ipatch) + wice_soisno(j_meth,ipatch), 0._r8)
            IF (snowdp(ipatch) > 0._r8 .and. snow_mass_total_m > 0._r8) THEN
               snow_layer_depth_m = snowdp(ipatch) * snow_layer_mass_m / snow_mass_total_m
            ELSEIF (snow_layer_mass_m > 0._r8) THEN
               snow_layer_depth_m = wliq_soisno(j_meth,ipatch)/denh2o + &
                  wice_soisno(j_meth,ipatch)/denice
            ELSE
               snow_layer_depth_m = 0._r8
            ENDIF
            dz_soisno_m(j_meth) = max(snow_layer_depth_m, 1.e-12_r8)
            z_soisno_m(j_meth) = zi_soisno_m(j_meth) - 0.5_r8*dz_soisno_m(j_meth)
            zi_soisno_m(j_meth-1) = zi_soisno_m(j_meth) - dz_soisno_m(j_meth)
         ENDDO
      ENDIF

   END SUBROUTINE methane_soisno_geometry

END MODULE MOD_MT_Methane_Impl
#endif
