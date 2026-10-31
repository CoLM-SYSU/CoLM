#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_Snow

   USE MOD_Precision
   USE MOD_MT_Defs, only: nmts, mt_tiny, mt_uses_generic_water
   USE MOD_MT_Solute, only: mt_equilibrate_dissolved
	   USE MOD_MT_Vars, only: mt_wliq_soisno, mt_wice_soisno, mt_pg_snow_ground, &
	      mt_scv, mt_sm_carry, mt_wetwat, mt_solid_soisno, mt_numerical_residual_step
   USE MOD_Namelist, only: DEF_USE_VariablySaturatedFlow

   IMPLICIT NONE

CONTAINS

   SUBROUTINE mt_newsnow (ipatch, patchtype, snl, snl_old, pg_snow, deltim, &
      scv, scv_bef, wetwat_val, &
      wliq_soisno, wice_soisno, wice_soisno_bef)
      IMPLICIT NONE
      integer,  intent(in) :: ipatch, patchtype, snl, snl_old
      real(r8), intent(in) :: pg_snow                    ! snow throughfall [mm/s]
      real(r8), intent(in) :: deltim
      real(r8), intent(in) :: scv                        ! post-newsnow scv [mm]
      real(r8), intent(in) :: scv_bef                    ! pre-newsnow scv [mm]
      real(r8), intent(in) :: wetwat_val                 ! current wetwat [mm]
      ! These arrays are only valid when snl < 0
      real(r8), intent(in), optional :: wliq_soisno(:)
      real(r8), intent(in), optional :: wice_soisno(:)
      real(r8), intent(in), optional :: wice_soisno_bef(:)

      integer  :: imt, j
      real(r8) :: R_snow_step, R_layer, d_wice, snow_mass_step
      real(r8) :: total_material, layer_water, ice_water, liq_water
      real(r8) :: ice_material, liq_material, moved

      IF (nmts <= 0) RETURN

      DO imt = 1, nmts
         IF (.not. mt_uses_generic_water(imt)) CYCLE
         snow_mass_step = max(pg_snow, 0._r8) * deltim
         IF (snow_mass_step > mt_tiny) THEN
            R_snow_step = mt_pg_snow_ground(imt, ipatch) / snow_mass_step
         ELSE
            R_snow_step = 0._r8
         ENDIF

#if (defined CoLMDEBUG)
         ! Invariant: when the pre-newsnow state has layered snow
         ! (snl_old < 0), mt_scv must be zero — pre-layer accumulation
         ! stops once a real snow layer exists, and Cases B/C below will
         ! unconditionally reset mt_scv. A non-zero residual here means
         ! somebody else wrote mt_scv while snl_old<0 (a real bug
         ! upstream, not a simple rounding artefact). Flag it once per
	         ! material per occurrence so the root cause can be traced; Cases
	         ! B/C now fold any residual into the top snow layer before reset.
         IF (snl_old < 0 .and. abs(mt_scv(imt, ipatch)) > mt_tiny) THEN
            write(*,'(A,I8,A,I3,A,E12.5)') &
               ' WARNING mt_newsnow: mt_scv residual under layered snow ipatch=', &
               ipatch, ' imt=', imt, ' mt_scv_pre=', mt_scv(imt, ipatch)
         ENDIF
#endif

         IF (snl_old == 0 .and. snl < 0) THEN
            ! Case A: First snow layer just created from accumulated scv.
            ! Partition the accumulated snow material across the new layer's
            ! ice/liquid pools.  Do not put all mt_scv into ice and then add
            ! liquid material again; that duplicates material mass when the first
            ! layer is born with liquid water.
            total_material = mt_scv(imt, ipatch) + mt_pg_snow_ground(imt, ipatch)
            j = snl + 1  ! = 0
            ice_water = 0._r8
            liq_water = 0._r8
            IF (present(wice_soisno)) THEN
               IF (size(wice_soisno) > 0) ice_water = max(wice_soisno(1), 0._r8)
            ENDIF
            IF (present(wliq_soisno)) THEN
               IF (size(wliq_soisno) > 0) liq_water = max(wliq_soisno(1), 0._r8)
            ENDIF
            layer_water = ice_water + liq_water
            IF (layer_water <= mt_tiny) layer_water = max(scv, mt_tiny)
            R_layer = total_material / layer_water
            ice_material = min(max(ice_water * R_layer, 0._r8), max(total_material, 0._r8))
            liq_material = min(max(liq_water * R_layer, 0._r8), max(total_material - ice_material, 0._r8))
            IF (ice_water <= mt_tiny .and. liq_water <= mt_tiny) ice_material = max(total_material, 0._r8)
            mt_wice_soisno(imt, j, ipatch) = ice_material
            mt_wliq_soisno(imt, j, ipatch) = liq_material
            ! Clear mt_scv: all transferred to the snow layer
            mt_scv(imt, ipatch) = 0._r8

         ELSEIF (snl < 0 .and. snl < snl_old) THEN
            ! Case B: New layer added to existing snow pack (deeper snow).
            ! This is a less common case (snow pack growing a new layer).
            ! R_snow_step is correct ONLY when the new layer's mass came from
            ! fresh snowfall this step. SnowLayersCombineDivide handles
            ! split/combine via dedicated material hooks (see file-end
            ! comment), so reaching this branch via re-distribution would
            ! quietly use fresh-snow R on existing-pack water and erase the
            ! aged R signal. A CoLMDEBUG warning fires when the new
            ! layer mass disagrees with pg_snow*dt by more than FP noise
            ! — covers the regression risk if a future split path forgets
            ! to route through the hooks.
            j = snl + 1
#if (defined CoLMDEBUG)
            IF (present(wice_soisno)) THEN
               IF (abs(wice_soisno(1) - max(pg_snow, 0._r8) * deltim) > 1.e-6_r8) THEN
                  write(*,'(A,I8,A,I3,A,E12.5,A,E12.5)') &
                     ' WARNING mt_newsnow Case B: wice_soisno(1) /= pg_snow*dt ipatch=', &
                     ipatch, ' imt=', imt, &
                     ' wice=', wice_soisno(1), &
                     ' pg_snow*dt=', max(pg_snow, 0._r8) * deltim
               ENDIF
            ENDIF
#endif
            IF (present(wice_soisno)) THEN
               mt_wice_soisno(imt, j, ipatch) = wice_soisno(1) * R_snow_step
            ENDIF
	            IF (present(wliq_soisno)) THEN
	               IF (size(wliq_soisno) > 0) THEN
	                  mt_wliq_soisno(imt, j, ipatch) = wliq_soisno(1) * R_snow_step
	               ENDIF
	            ENDIF
	            IF (abs(mt_scv(imt, ipatch)) > mt_tiny) THEN
	               mt_wice_soisno(imt, j, ipatch) = mt_wice_soisno(imt, j, ipatch) &
	                  + mt_scv(imt, ipatch)
	            ENDIF
	            ! mt_scv should be 0 when layers already exist
	            mt_scv(imt, ipatch) = 0._r8

         ELSEIF (snl < 0 .and. snl == snl_old) THEN
            ! Case C: Snow added to existing top layer (no new node).
            ! Do not route through mt_scv here: once layered snow exists,
            ! the water-side pg_snow increment goes directly into the top
            ! snow ice layer. A previous scv -> d_wice-gated -> clear path
            ! could silently drop ground-snow material if mt_precip's
            ! canopy residual correction and the water increment differed by
            ! more than mt_tiny.
            j = snl + 1
            IF (present(wice_soisno) .and. present(wice_soisno_bef)) THEN
               d_wice = wice_soisno(1) - wice_soisno_bef(1)
	               IF (mt_pg_snow_ground(imt, ipatch) > mt_tiny) THEN
	                  mt_wice_soisno(imt, j, ipatch) = mt_wice_soisno(imt, j, ipatch) &
	                     + mt_pg_snow_ground(imt, ipatch)
	               ENDIF
#if (defined CoLMDEBUG)
               IF (abs(d_wice - snow_mass_step) > 1.e-6_r8 .and. &
                   mt_pg_snow_ground(imt, ipatch) > mt_tiny) THEN
                  write(*,'(A,I8,A,I3,A,E12.5,A,E12.5,A,E12.5)') &
                     ' WARNING mt_newsnow Case C: d_wice /= pg_snow*dt ipatch=', &
                     ipatch, ' imt=', imt, ' d_wice=', d_wice, &
                     ' pg_snow*dt=', snow_mass_step, &
                     ' mt_pg_snow_ground=', mt_pg_snow_ground(imt, ipatch)
	               ENDIF
#endif
	               IF (abs(mt_scv(imt, ipatch)) > mt_tiny) THEN
	                  mt_wice_soisno(imt, j, ipatch) = mt_wice_soisno(imt, j, ipatch) &
	                     + mt_scv(imt, ipatch)
	               ENDIF
            ENDIF
            ! mt_scv should be 0 when layers exist
            mt_scv(imt, ipatch) = 0._r8

         ELSE
            ! snl == 0: no snow layer yet. Accumulate in mt_scv until
            ! the first real snow layer is created.
            !
            ! scv sink handling (warm wetland melt, thin-snow melt) is done
            ! AFTER THERMAL in CoLMMAIN, not here. mt_newsnow only
            ! accumulates and transfers on layer creation.
            !
            ! One exception: warm wetland zeroing scv happens INSIDE newsnow
            ! (before THERMAL), so handle it here.
            mt_scv(imt, ipatch) = mt_scv(imt, ipatch) + mt_pg_snow_ground(imt, ipatch)
            IF (patchtype == 2 .and. scv < mt_tiny) THEN
               ! Warm wetland: newsnow moved scv to wetwat (variably saturated
               ! flow) or dropped it (otherwise), and set scv = 0. Follow it;
               ! a drop is booked as a numerical sink like the host's loss.
               IF (DEF_USE_VariablySaturatedFlow) THEN
                  mt_wetwat(imt, ipatch) = mt_wetwat(imt, ipatch) + mt_scv(imt, ipatch)
               ELSEIF (allocated(mt_numerical_residual_step)) THEN
                  mt_numerical_residual_step(imt, ipatch) = &
                     mt_numerical_residual_step(imt, ipatch) - mt_scv(imt, ipatch)
               ENDIF
               mt_scv(imt, ipatch) = 0._r8
            ENDIF
         ENDIF
         ! Warm wetland under layered snow: newsnow (as in the host) adds the
         ! whole scv, i.e. the water of all snow layers, to wetwat with
         ! variably saturated flow and sets scv = 0, but keeps the layers.
         ! The material follows the water: the layers keep their material and the
         ! same amount is added to wetwat, booked as a numerical source.
         IF (patchtype == 2 .and. snl < 0 .and. snl_old < 0 .and. scv < mt_tiny .and. &
             DEF_USE_VariablySaturatedFlow) THEN
            moved = sum(mt_wice_soisno(imt, snl+1:0, ipatch)) &
                  + sum(mt_wliq_soisno(imt, snl+1:0, ipatch))
            mt_wetwat(imt, ipatch) = mt_wetwat(imt, ipatch) + moved
            IF (allocated(mt_numerical_residual_step)) &
               mt_numerical_residual_step(imt, ipatch) = &
                  mt_numerical_residual_step(imt, ipatch) + moved
         ENDIF
         IF (snl < 0 .and. present(wliq_soisno)) THEN
            j = snl + 1
            IF (size(wliq_soisno) > 0) THEN
               CALL mt_equilibrate_dissolved(imt, max(wliq_soisno(1), 0._r8), &
                  mt_wliq_soisno(imt, j, ipatch), mt_solid_soisno(imt, j, ipatch))
            ENDIF
         ENDIF
      ENDDO
   END SUBROUTINE mt_newsnow

   SUBROUTINE mt_snowmelt_carry (ipatch, snl, scv, scv_before)
      integer, intent(in) :: ipatch, snl
      real(r8), intent(in) :: scv, scv_before
      integer :: imt_loc
      real(r8) :: ratio_loc

         ! snow melt carried from a vanishing or shrinking no-layer snow pack
         IF (snl == 0) THEN
            IF (scv < mt_tiny) THEN
               DO imt_loc = 1, nmts
                  IF (.not. mt_uses_generic_water(imt_loc)) CYCLE
                  mt_sm_carry(imt_loc, ipatch) = mt_scv(imt_loc, ipatch)
                  mt_scv(imt_loc, ipatch) = 0._r8
               ENDDO
            ELSEIF (scv < scv_before - mt_tiny) THEN
               ratio_loc = scv / max(scv_before, mt_tiny)
               ratio_loc = max(min(ratio_loc, 1._r8), 0._r8)
               DO imt_loc = 1, nmts
                  IF (.not. mt_uses_generic_water(imt_loc)) CYCLE
                  mt_sm_carry(imt_loc, ipatch) = mt_scv(imt_loc, ipatch) * (1._r8 - ratio_loc)
                  mt_scv(imt_loc, ipatch) = mt_scv(imt_loc, ipatch) * ratio_loc
               ENDDO
            ELSE
               DO imt_loc = 1, nmts
                  IF (.not. mt_uses_generic_water(imt_loc)) CYCLE
                  mt_sm_carry(imt_loc, ipatch) = 0._r8
               ENDDO
            ENDIF
         ELSE
            DO imt_loc = 1, nmts
               IF (.not. mt_uses_generic_water(imt_loc)) CYCLE
               mt_sm_carry(imt_loc, ipatch) = 0._r8
            ENDDO
         ENDIF

   END SUBROUTINE mt_snowmelt_carry

   ! material_snow_layer_adj was removed: the prior post-hoc redistribution
   ! (total-material / total-water share) homogenised the snow column every
   ! time snowlayerscombine / snowlayersdivide fired, erasing vertical
   ! isotope gradients that the model is otherwise supposed to preserve.
   ! The fix carries mt_wliq / mt_wice / mt_scv through the exact same
   ! per-layer topology as water via the new optional hooks on
   ! MOD_SnowLayersCombineDivide's combine/divide routines.

END MODULE MOD_MT_Snow
#endif
