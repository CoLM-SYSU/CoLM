#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_HydroDiagnostics

   USE MOD_Precision
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: mt_water_diagnostics

   ! Per-call water-path observations; the host water solution never reads them.
   TYPE mt_water_diagnostics
      ! Masses are mm per time step, except etroot and qgtop (mm/s).
      ! Actual root uptake is signed: negative values return water to the store.
      real(r8), allocatable :: etroot(:), etroot_actual(:), rsub_source_layer(:)
      real(r8), allocatable :: wblc_ice_sink(:) ! ice removed by the water balance correction
      real(r8), allocatable :: snow_qout_layer(:)
      logical, allocatable :: permeable_soil(:)
      real(r8) :: etroot_aquifer, etroot_surface
      real(r8) :: rsub_source_surface, rsub_source_aquifer
      ! WATER_VSF adds dew/frost directly to layer 1, so these transfers stay zero.
      real(r8) :: dew_overflow, frost_displaced, late_runoff
      real(r8) :: imperv_evap_wdsrf, imperv_evap_soil, imperv_subl_soil
      real(r8) :: qgtop
      real(r8), private :: root_water_before, surface_before, layer1_before
      real(r8), private :: aquifer_before, ponding_before
      real(r8), allocatable, private :: exchange_before(:)
   CONTAINS
      PROCEDURE :: reset => reset_water_diagnostics
      PROCEDURE :: reset_soil => reset_soil_diagnostics
      PROCEDURE :: begin_root_uptake
      PROCEDURE :: end_root_uptake
      PROCEDURE :: begin_aquifer_exchange
      PROCEDURE :: end_aquifer_exchange
      PROCEDURE :: begin_impervious_evaporation
      PROCEDURE :: end_impervious_evaporation
      PROCEDURE :: record_ice_sink
   END TYPE mt_water_diagnostics

CONTAINS

   SUBROUTINE reset_soil_diagnostics (this, nlev)
      class(mt_water_diagnostics), intent(inout) :: this
      integer, intent(in) :: nlev

      ! The soil solver must preserve diagnostics already sampled by WATER_VSF.
      IF (allocated(this%etroot)) THEN
         IF (size(this%etroot) /= nlev) &
            deallocate(this%etroot, this%etroot_actual, this%rsub_source_layer)
      ENDIF
      IF (.not. allocated(this%etroot)) &
         allocate(this%etroot(nlev), this%etroot_actual(nlev), this%rsub_source_layer(nlev))
      this%etroot = 0._r8
      this%etroot_actual = 0._r8
      this%etroot_aquifer = 0._r8
      this%etroot_surface = 0._r8
      this%rsub_source_layer = 0._r8
      this%rsub_source_surface = 0._r8
      this%rsub_source_aquifer = 0._r8
   END SUBROUTINE reset_soil_diagnostics

   SUBROUTINE reset_water_diagnostics (this, nl_soil, lb)
      class(mt_water_diagnostics), intent(inout) :: this
      integer, intent(in) :: nl_soil, lb

      CALL this%reset_soil(nl_soil)
      IF (allocated(this%wblc_ice_sink)) THEN
         IF (size(this%wblc_ice_sink) /= nl_soil) &
            deallocate(this%wblc_ice_sink, this%permeable_soil)
      ENDIF
      IF (.not. allocated(this%wblc_ice_sink)) &
         allocate(this%wblc_ice_sink(nl_soil), this%permeable_soil(nl_soil))
      IF (allocated(this%snow_qout_layer)) THEN
         IF (size(this%snow_qout_layer) /= max(1-lb, 0)) deallocate(this%snow_qout_layer)
      ENDIF
      IF (.not. allocated(this%snow_qout_layer)) allocate(this%snow_qout_layer(lb:0))
      this%wblc_ice_sink = 0._r8
      this%permeable_soil = .false.
      this%snow_qout_layer = 0._r8
      this%dew_overflow = 0._r8
      this%frost_displaced = 0._r8
      this%late_runoff = 0._r8
      this%imperv_evap_wdsrf = 0._r8
      this%imperv_evap_soil = 0._r8
      this%imperv_subl_soil = 0._r8
      this%qgtop = 0._r8
   END SUBROUTINE reset_water_diagnostics

   SUBROUTINE begin_root_uptake(this, liquid, thickness)
      class(mt_water_diagnostics), intent(inout) :: this
      real(r8), intent(in) :: liquid, thickness
      this%root_water_before = liquid * thickness
   END SUBROUTINE begin_root_uptake

   SUBROUTINE end_root_uptake(this, layer, liquid, thickness)
      class(mt_water_diagnostics), intent(inout) :: this
      integer, intent(in) :: layer
      real(r8), intent(in) :: liquid, thickness
      ! Capture actual uptake after the host clips liquid to [0, porosity].
      this%etroot_actual(layer) = this%root_water_before - liquid * thickness
   END SUBROUTINE end_root_uptake

   SUBROUTINE begin_aquifer_exchange(this, nlev, sp_zi, permeable, porosity, liquid, zwt, wa, ponding)
      class(mt_water_diagnostics), intent(inout) :: this
      integer, intent(in) :: nlev
      logical, intent(in) :: permeable(nlev)
      real(r8), intent(in) :: sp_zi(0:nlev), porosity(nlev), liquid(nlev), zwt, wa, ponding
      IF (allocated(this%exchange_before)) THEN
         IF (size(this%exchange_before) /= nlev) deallocate(this%exchange_before)
      ENDIF
      IF (.not. allocated(this%exchange_before)) allocate(this%exchange_before(nlev))
      this%ponding_before = ponding
      this%aquifer_before = wa
      CALL mt_exchange_layer_water(nlev, sp_zi, permeable, porosity, liquid, zwt, this%exchange_before)
   END SUBROUTINE begin_aquifer_exchange

   SUBROUTINE end_aquifer_exchange(this, nlev, sp_zi, permeable, porosity, liquid, zwt, wa, ponding, drain, et_left)
      class(mt_water_diagnostics), intent(inout) :: this
      integer, intent(in) :: nlev
      logical, intent(in) :: permeable(nlev)
      real(r8), intent(in) :: sp_zi(0:nlev), porosity(nlev), liquid(nlev), zwt, wa, ponding, drain, et_left
      real(r8) :: layer_after(nlev)
      ! Split the single host net exchange into root and runoff pathways.
      CALL mt_exchange_layer_water(nlev, sp_zi, permeable, porosity, liquid, zwt, layer_after)
      CALL mt_split_aquifer_exchange(nlev, drain, et_left, &
         this%exchange_before - layer_after, this%aquifer_before - wa, this%ponding_before - ponding, &
         this%etroot_actual, this%etroot_aquifer, this%etroot_surface, &
         this%rsub_source_layer, this%rsub_source_aquifer, this%rsub_source_surface)
   END SUBROUTINE end_aquifer_exchange

   SUBROUTINE begin_impervious_evaporation(this, surface, liquid)
      class(mt_water_diagnostics), intent(inout) :: this
      real(r8), intent(in) :: surface, liquid
      this%surface_before = surface
      this%layer1_before = liquid
   END SUBROUTINE begin_impervious_evaporation

   SUBROUTINE end_impervious_evaporation(this, surface, liquid)
      class(mt_water_diagnostics), intent(inout) :: this
      real(r8), intent(in) :: surface, liquid
      ! Host evaporation uses ponding first, then layer-1 liquid, clipping both.
      this%imperv_evap_wdsrf = this%surface_before - surface
      this%imperv_evap_soil = this%layer1_before - liquid
   END SUBROUTINE end_impervious_evaporation

   SUBROUTINE record_ice_sink(this, layer, ice, imbalance)
      class(mt_water_diagnostics), intent(inout) :: this
      integer, intent(in) :: layer
      real(r8), intent(in) :: ice, imbalance
      IF (ice > imbalance) THEN
         this%wblc_ice_sink(layer) = imbalance
      ELSE
         this%wblc_ice_sink(layer) = ice
      ENDIF
   END SUBROUTINE record_ice_sink

   ! --- material transport: water held by each layer ---
   ! Liquid water of each layer counting the saturated part below the water
   ! table at porosity, i.e. the storage soilwater_aquifer_exchange moves.
   SUBROUTINE mt_exchange_layer_water (nlev, sp_zi, is_permeable, porsl, ss_vliq, zwt, wlayer)

   IMPLICIT NONE

   integer,  intent(in)  :: nlev
   real(r8), intent(in)  :: sp_zi (0:nlev)
   logical,  intent(in)  :: is_permeable (1:nlev)
   real(r8), intent(in)  :: porsl (1:nlev)
   real(r8), intent(in)  :: ss_vliq (1:nlev)
   real(r8), intent(in)  :: zwt
   real(r8), intent(out) :: wlayer (1:nlev)

   integer :: ilev

      DO ilev = 1, nlev
         wlayer(ilev) = ss_vliq(ilev) * (sp_zi(ilev) - sp_zi(ilev-1))
         IF (is_permeable(ilev) .and. zwt < sp_zi(ilev)) THEN
            wlayer(ilev) = ss_vliq(ilev) * max(zwt - sp_zi(ilev-1), 0._r8) &
               + porsl(ilev) * min(sp_zi(ilev) - zwt, sp_zi(ilev) - sp_zi(ilev-1))
         ENDIF
      ENDDO

   END SUBROUTINE mt_exchange_layer_water

   ! --- material transport: split the net aquifer exchange ---
   ! The host moves the net amount drain + et_left (subsurface runoff plus the
   ! root demand that the unsaturated layers could not meet, negative when
   ! roots release water) in one soilwater_aquifer_exchange call. Given the
   ! loss of each store in that call (layers, aquifer, ponding water), split
   ! it into root water and subsurface runoff:
   !  * net removal, et_left >= 0: each store's loss is shared in proportion
   !    et_left : drain;
   !  * net removal, et_left < 0: the released root water is put into the
   !    losing stores in proportion to their loss, and leaves with the runoff;
   !  * net addition (et_left < 0 and larger than drain): the released water
   !    enters the gaining stores, the runoff leaves from them in proportion.
   ! In every case root water minus runoff of a store equals its net loss.
   SUBROUTINE mt_split_aquifer_exchange (nlev, drain, et_left, &
         dlayer, dwa, ddp, etroot_layer, etroot_wa, etroot_dp, &
         rsub_layer, rsub_wa, rsub_dp)

   IMPLICIT NONE

   integer,  intent(in)    :: nlev
   real(r8), intent(in)    :: drain, et_left    ! mm
   real(r8), intent(in)    :: dlayer(1:nlev), dwa, ddp  ! loss of each store (mm)
   real(r8), intent(inout) :: etroot_layer(1:nlev), etroot_wa, etroot_dp
   real(r8), intent(out)   :: rsub_layer(1:nlev), rsub_wa, rsub_dp

   real(r8) :: share(1:nlev+2), w(1:nlev+2), wsum, et_fraction, release

      rsub_layer = 0._r8
      rsub_wa    = 0._r8
      rsub_dp    = 0._r8

      IF (drain + et_left > 0._r8) THEN
         ! net removal; gains of single stores are not expected here
         share(1:nlev) = max(dlayer, 0._r8)
         share(nlev+1) = max(dwa, 0._r8)
         share(nlev+2) = max(ddp, 0._r8)
         IF (et_left >= 0._r8) THEN
            et_fraction = min(et_left / (drain + et_left), 1._r8)
            etroot_layer = etroot_layer + share(1:nlev) * et_fraction
            etroot_wa    = share(nlev+1) * et_fraction
            etroot_dp    = share(nlev+2) * et_fraction
            rsub_layer   = share(1:nlev) * (1._r8 - et_fraction)
            rsub_wa      = share(nlev+1) * (1._r8 - et_fraction)
            rsub_dp      = share(nlev+2) * (1._r8 - et_fraction)
         ELSE
            release = - et_left
            wsum = sum(share)
            IF (wsum > 0._r8) THEN
               w = share / wsum
            ELSE
               w = 0._r8;  w(nlev+1) = 1._r8
            ENDIF
            etroot_layer = etroot_layer - release * w(1:nlev)
            etroot_wa    = - release * w(nlev+1)
            etroot_dp    = - release * w(nlev+2)
            rsub_layer   = share(1:nlev) + release * w(1:nlev)
            rsub_wa      = share(nlev+1) + release * w(nlev+1)
            rsub_dp      = share(nlev+2) + release * w(nlev+2)
         ENDIF
      ELSEIF (et_left < 0._r8) THEN
         ! net addition of released root water
         share(1:nlev) = max(-dlayer, 0._r8)
         share(nlev+1) = max(-dwa, 0._r8)
         share(nlev+2) = max(-ddp, 0._r8)
         wsum = sum(share)
         IF (wsum > 0._r8) THEN
            w = share / wsum
         ELSE
            w = 0._r8;  w(nlev+1) = 1._r8
         ENDIF
         etroot_layer = etroot_layer - share(1:nlev) - drain * w(1:nlev)
         etroot_wa    = - share(nlev+1) - drain * w(nlev+1)
         etroot_dp    = - share(nlev+2) - drain * w(nlev+2)
         rsub_layer   = drain * w(1:nlev)
         rsub_wa      = drain * w(nlev+1)
         rsub_dp      = drain * w(nlev+2)
      ENDIF

#ifdef MaterialTransport
      ! Diagnostic for the aquifer-debt investigation: the split treats the
      ! host wa drop as an aquifer loss, so a host wa driven far below zero
      ! (or an unexpectedly large drain) shows up here.  Print-only; it does
      ! not change the split.
      IF (abs(dwa) > 1.0e3_r8 .or. abs(drain) > 1.0e3_r8) THEN
         write(*,'(A,6ES16.8)') 'MT_DIAG split_aquifer drain,dwa,ddp,et_left,rsub_wa,share_wa: ', &
            drain, dwa, ddp, et_left, rsub_wa, share(nlev+1)
      ENDIF
#endif

   END SUBROUTINE mt_split_aquifer_exchange

END MODULE MOD_MT_HydroDiagnostics
#endif
