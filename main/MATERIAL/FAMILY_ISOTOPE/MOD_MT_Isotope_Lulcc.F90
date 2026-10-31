#include <define.h>

#if (defined MaterialTransport) && (defined LULCC)
MODULE MOD_MT_Isotope_Lulcc
!-----------------------------------------------------------------------
! Isotope aquifer reference carrier across a land cover change.
!
! With variably saturated flow, every soil and wetland patch carries a
! reference aquifer water Vref = DEF_MT_AQUIFER_MIXING_WATER_MM with the
! anchor mass Mref = Vref * ref_ratio; the actual aquifer isotope is
! Mref + mt_wa in Vref + wa (see MOD_MT_Rest cold start). Lake, glacier,
! urban and other special patches carry no reference (Vref = Mref = 0).
!
! remap_land_mt_lulcc_state moves Vref and Mref by area like every
! other pool. Two outcomes cannot be valid on the new patches, and the run
! used to stop there:
!
!  - a special patch received reference water (soil area turned into, for
!    example, a lake). Its reference carrier is taken out: the isotope it
!    holds at the patch's actual aquifer ratio leaves, the rest stays in
!    mt_wa, and Vref = Mref = 0.
!
!  - a soil or wetland patch received no reference water (all its area
!    came from special patches, or with the SAT scheme no old patch of
!    the same class existed in the element). It gets the cold-start
!    carrier, Vref = DEF_MT_AQUIFER_MIXING_WATER_MM, at the patch's own
!    actual aquifer ratio (the initial ratio when its aquifer is empty).
!
! The isotope mass that leaves or enters with the carrier is returned per
! material (mass per patch area times patch area), so the MEC mass check
! can count it. Patches already valid after the area remap are not touched,
! so every land cover change that ran before gives the same result.
!-----------------------------------------------------------------------

   USE MOD_Precision
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: mt_isotope_lulcc_aquifer_reference

CONTAINS

   SUBROUTINE mt_isotope_lulcc_aquifer_reference (patchclass_new, wa_new, &
         ref_water, ref_mass, mt_wa, patch_area, carrier_change)

      USE MOD_SPMD_Task,   only: CoLM_stop, p_iam_glb
      USE MOD_Const_LC,    only: patchtypes
      USE MOD_Namelist,    only: DEF_USE_VariablySaturatedFlow
      USE MOD_MT_Namelist, only: DEF_MT_AQUIFER_MIXING_WATER_MM
      USE MOD_MT_Defs,     only: nmts, materials, mt_uses_generic_water, material_is_isotope, &
                                 material_init_water_ratio, material_aquifer_isotope_ratio

      integer,  intent(in)    :: patchclass_new(:)
      real(r8), intent(in)    :: wa_new(:)         ! host aquifer water on the new patches [mm]
      real(r8), intent(inout) :: ref_water(:)      ! Vref [mm]
      real(r8), intent(inout) :: ref_mass(:,:)     ! Mref (material, patch)
      real(r8), intent(inout) :: mt_wa(:,:)       ! relative aquifer mass (material, patch)
      real(r8), intent(in),    optional :: patch_area(:)
      real(r8), intent(inout), optional :: carrier_change(:)  ! mass entered - left, per material

      integer  :: np, imt, ptype, ncreated, nremoved
      real(r8) :: R_actual, area
      logical  :: has_isotope

      IF (.not. DEF_USE_VariablySaturatedFlow .or. DEF_MT_AQUIFER_MIXING_WATER_MM <= 0._r8) RETURN
      has_isotope = .false.
      DO imt = 1, nmts
         IF (mt_uses_generic_water(imt) .and. material_is_isotope(imt)) has_isotope = .true.
      ENDDO
      IF (.not. has_isotope) RETURN

      ncreated = 0
      nremoved = 0
      DO np = 1, size(patchclass_new)
         IF (patchclass_new(np) < lbound(patchtypes, 1) .or. &
             patchclass_new(np) > ubound(patchtypes, 1)) &
            CALL CoLM_stop('invalid LULCC patch class for aquifer isotope reference')
         IF (np > size(wa_new)) CALL CoLM_stop('MT LULCC isotope reference: missing host aquifer water')
         ptype = patchtypes(patchclass_new(np))
         area = 0._r8
         IF (present(patch_area)) area = max(patch_area(np), 0._r8)

         IF (ptype /= 0 .and. ptype /= 2 .and. ref_water(np) > 0._r8) THEN
            ! Special patch: take the reference carrier out at the actual ratio.
            DO imt = 1, nmts
               IF (.not. (mt_uses_generic_water(imt) .and. material_is_isotope(imt))) CYCLE
               R_actual = material_aquifer_isotope_ratio(wa_new(np), mt_wa(imt,np), &
                  ref_water(np), ref_mass(imt,np), ref_mass(imt,np) / ref_water(np))
               mt_wa(imt,np) = mt_wa(imt,np) + ref_mass(imt,np) - ref_water(np) * R_actual
               IF (present(carrier_change)) &
                  carrier_change(imt) = carrier_change(imt) - ref_water(np) * R_actual * area
               ref_mass(imt,np) = 0._r8
            ENDDO
            ref_water(np) = 0._r8
            nremoved = nremoved + 1

         ELSEIF ((ptype == 0 .or. ptype == 2) .and. ref_water(np) <= 0._r8) THEN
            ! Soil or wetland patch without reference: cold-start carrier at the
            ! patch's own actual aquifer ratio.
            DO imt = 1, nmts
               IF (.not. (mt_uses_generic_water(imt) .and. material_is_isotope(imt))) CYCLE
               R_actual = material_aquifer_isotope_ratio(wa_new(np), mt_wa(imt,np) + ref_mass(imt,np), &
                  0._r8, 0._r8, material_init_water_ratio(imt))
               mt_wa(imt,np) = mt_wa(imt,np) + ref_mass(imt,np) &
                  + DEF_MT_AQUIFER_MIXING_WATER_MM * (R_actual - materials(imt)%ref_ratio)
               ref_mass(imt,np) = DEF_MT_AQUIFER_MIXING_WATER_MM * materials(imt)%ref_ratio
               IF (present(carrier_change)) &
                  carrier_change(imt) = carrier_change(imt) &
                     + DEF_MT_AQUIFER_MIXING_WATER_MM * R_actual * area
            ENDDO
            ref_water(np) = DEF_MT_AQUIFER_MIXING_WATER_MM
            IF (wa_new(np) + ref_water(np) < 0._r8) &
               CALL CoLM_stop('negative isotope aquifer mixing water after LULCC')
            ncreated = ncreated + 1
         ENDIF
      ENDDO

      IF (ncreated + nremoved > 0) THEN
         write(*,'(A,I6,A,I8,A,I8)') ' MT LULCC isotope aquifer reference on process', p_iam_glb, &
            ': created on', ncreated, ' soil/wetland patches, removed from', nremoved
      ENDIF

   END SUBROUTINE mt_isotope_lulcc_aquifer_reference

END MODULE MOD_MT_Isotope_Lulcc
#endif
