#include <define.h>

#ifdef MaterialTransport
MODULE MOD_MT_LandWork

   USE MOD_Precision
   USE MOD_MT_HydroDiagnostics, only: mt_water_diagnostics
   USE MOD_MT_ThermalDiagnostics, only: mt_thermal_diagnostics

   IMPLICIT NONE
   PRIVATE

   ! Patch-local snapshots and scratch; no persistent state or restart fields.
   TYPE, PUBLIC :: mt_land_work
      real(r8), allocatable :: wliq_before(:), wice_before(:)
      type(mt_thermal_diagnostics) :: thermal
      real(r8), allocatable :: snow_ice_before(:)
      real(r8) :: balance_error, wa_before, wdsrf_before, wetwat_before
      real(r8) :: canopy_rain_before, canopy_snow_before
      real(r8) :: canopy_rain_intercepted, canopy_snow_intercepted
      real(r8) :: canopy_rain_prethermal, canopy_snow_prethermal
      real(r8) :: canopy_rain_postthermal, canopy_snow_postthermal
      real(r8) :: gross_rain, gross_snow, rain_release, snow_release
      real(r8) :: snow_before, glacier_overflow
      real(r8) :: recharge, infiltration, irrigation
      real(r8) :: irrigation_storage_before, irrigation_storage_ground
      real(r8) :: evaporation, flood_input
#ifdef GridRiverLakeFlow
      real(r8) :: flood_evap_temperature
#endif
   CONTAINS
      PROCEDURE :: prepare_snow => prepare_land_snow
      PROCEDURE :: prepare_thermal => prepare_land_thermal
      PROCEDURE :: save_water => save_land_water
      PROCEDURE :: release => release_land_work
   END TYPE mt_land_work

CONTAINS

   SUBROUTINE prepare_land_snow (this, maxsnl, nl_soil, snl, wice_soisno)
      class(mt_land_work), intent(inout) :: this
      integer, intent(in) :: maxsnl, nl_soil, snl
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil)

      allocate (this%snow_ice_before(maxsnl+1:nl_soil))
      this%snow_ice_before = 0._r8
      IF (snl < 0) this%snow_ice_before(snl+1:0) = wice_soisno(snl+1:0)
   END SUBROUTINE prepare_land_snow

   SUBROUTINE prepare_land_thermal (this, lb, nl_soil, wliq_soisno, wice_soisno, &
      wa, wdsrf, wetwat, water_diag, t_soisno, fsno, split_soilsnow)
      class(mt_land_work), intent(inout) :: this
      integer, intent(in) :: lb, nl_soil
      real(r8), intent(in) :: wliq_soisno(lb:nl_soil), wice_soisno(lb:nl_soil)
      real(r8), intent(in) :: wa, wdsrf, wetwat, t_soisno(lb:nl_soil), fsno
      logical, intent(in) :: split_soilsnow
      type(mt_water_diagnostics), intent(inout) :: water_diag

      allocate (this%wliq_before(lb:nl_soil), this%wice_before(lb:nl_soil))
      CALL this%save_water(lb, nl_soil, wliq_soisno, wice_soisno)
      this%wa_before = wa
      this%wdsrf_before = wdsrf
      this%wetwat_before = wetwat
      this%recharge = 0._r8
      CALL water_diag%reset(nl_soil, min(lb,0))
#ifdef GridRiverLakeFlow
      this%flood_evap_temperature = t_soisno(lb)
      IF (split_soilsnow) this%flood_evap_temperature = &
         fsno*t_soisno(lb) + (1._r8-fsno)*t_soisno(1)
#endif
   END SUBROUTINE prepare_land_thermal

   SUBROUTINE save_land_water (this, lb, nl_soil, wliq, wice)
      class(mt_land_work), intent(inout) :: this
      integer, intent(in) :: lb, nl_soil
      real(r8), intent(in) :: wliq(lb:nl_soil), wice(lb:nl_soil)

      this%wliq_before(lb:nl_soil) = wliq(lb:nl_soil)
      this%wice_before(lb:nl_soil) = wice(lb:nl_soil)
   END SUBROUTINE save_land_water

   SUBROUTINE release_land_work (this)
      class(mt_land_work), intent(inout) :: this

      deallocate (this%wliq_before, this%wice_before)
      deallocate (this%thermal%soil_thaw, this%thermal%soil_freeze)
      deallocate (this%snow_ice_before)
   END SUBROUTINE release_land_work

END MODULE MOD_MT_LandWork
#endif
