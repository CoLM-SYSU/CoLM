#include <define.h>
#if (defined MaterialTransport && defined BGC)

MODULE MOD_MT_NitrogenLeaching

!-----------------------------------------------------------------------
! Map the BGC soil mineral-N leaching loss onto the material-transport
! solute families NH4 and NO3.
!
! BGC computes the loss in MOD_BGC_Soil_BiogeochemNLeaching:
!   nitrification ON : smin_no3_leached_vr (sub-surface drainage) and
!                      smin_no3_runoff_vr  (surface runoff), for NO3, plus
!                      smin_nh4_leached_vr for NH4.  The NH4 row exists only
!                      in a MaterialTransport build and is zero unless
!                      DEF_USE_NH4_LEACHING is on.
!   nitrification OFF: sminn_leached_vr for total mineral N.  With no
!                      nitrification the mineral pool is NH4, so that flux
!                      is booked to NH4.
!
! This module deliberately does NOT change BGC physics: it reads the
! fluxes BGC already removed from its mineral-N pools and books the same
! mass as an external solute input on the receiving water pathway.
! Sub-surface leaching enters the base-flow bookkeeping (a_mt_rsub),
! surface runoff enters a_mt_rsur; both are added to the river hand-off
! (mt_rnof_step) and to a_mt_bgc_leach so the material balance check
! counts them as input without treating them as precipitation.
!
! Units: BGC fluxes are gN m-3 s-1 per soil layer, so the mass moved in
! one step is flux * dz_soi * deltim / 1000 [kg N m-2].  NH4 and NO3
! solute concentrations in the material transport are kg N per kg water;
! the tracked mass is kg N.
!
! Only soil patches (patchtype 0) run BGC; other patch types return
! immediately, and species other than NH4/NO3 are untouched.
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Vars_Global,        only: nl_soil, dz_soi
   USE MOD_Namelist,           only: DEF_USE_NITRIF
   USE MOD_MT_BGC_NH4,         only: smin_nh4_leached_vr
   USE MOD_BGC_Vars_1DFluxes,  only: sminn_leached_vr, &
                                     smin_no3_leached_vr, smin_no3_runoff_vr
   USE MOD_MT_Defs,            only: nmts, materials, material_lower, material_is_solute
   USE MOD_MT_Vars,            only: mt_rnof_step, a_mt_rsub, a_mt_rsur, &
                                     a_mt_rnof, a_mt_bgc_leach

   IMPLICIT NONE
   PUBLIC :: mt_bgc_n_leaching

CONTAINS

   SUBROUTINE mt_bgc_n_leaching (ipatch, patchtype, deltim)

      IMPLICIT NONE
      integer,  intent(in) :: ipatch, patchtype
      real(r8), intent(in) :: deltim

      integer :: imt, j
      real(r8) :: sub_mass, runoff_mass
      character(len=32) :: lname

      IF (patchtype /= 0) RETURN
      IF (nmts <= 0) RETURN

      DO imt = 1, nmts
         IF (.not. material_is_solute(imt)) CYCLE
         lname = material_lower(trim(materials(imt)%name))
         IF (lname /= 'nh4' .and. lname /= 'no3') CYCLE

         sub_mass    = 0._r8
         runoff_mass = 0._r8

         IF (DEF_USE_NITRIF) THEN
            ! Nitrification on: BGC computes species-resolved NH4 and NO3
            ! leaching (NH4 only when DEF_USE_NH4_LEACHING is on; otherwise
            ! that row is zero), plus NO3 surface runoff.
            IF (lname == 'nh4') THEN
               DO j = 1, nl_soil
                  sub_mass = sub_mass + smin_nh4_leached_vr(j,ipatch) &
                             * dz_soi(j) * deltim
               ENDDO
            ELSEIF (lname == 'no3') THEN
               DO j = 1, nl_soil
                  sub_mass = sub_mass + smin_no3_leached_vr(j,ipatch) &
                             * dz_soi(j) * deltim
                  runoff_mass = runoff_mass + smin_no3_runoff_vr(j,ipatch) &
                                * dz_soi(j) * deltim
               ENDDO
            ENDIF
         ELSE
            ! Nitrification off: the mineral-N pool is NH4, so the aggregate
            ! leaching flux is booked to NH4 (upstream BGC computes no
            ! species-resolved NH4 row in this branch).
            IF (lname == 'nh4') THEN
               DO j = 1, nl_soil
                  sub_mass = sub_mass + sminn_leached_vr(j,ipatch) &
                             * dz_soi(j) * deltim
               ENDDO
            ENDIF
         ENDIF

         ! gN m-2 -> kg N m-2
         sub_mass    = sub_mass    / 1000._r8
         runoff_mass = runoff_mass / 1000._r8
         IF (sub_mass + runoff_mass <= 0._r8) CYCLE

         a_mt_rsub     (imt,ipatch) = a_mt_rsub     (imt,ipatch) + sub_mass
         a_mt_rsur     (imt,ipatch) = a_mt_rsur     (imt,ipatch) + runoff_mass
         a_mt_rnof     (imt,ipatch) = a_mt_rnof     (imt,ipatch) &
                                     + sub_mass + runoff_mass
         mt_rnof_step   (imt,ipatch) = mt_rnof_step   (imt,ipatch) &
                                     + sub_mass + runoff_mass
         a_mt_bgc_leach(imt,ipatch) = a_mt_bgc_leach(imt,ipatch) &
                                     + sub_mass + runoff_mass
      ENDDO

   END SUBROUTINE mt_bgc_n_leaching

END MODULE MOD_MT_NitrogenLeaching
#endif
