#include <define.h>

#ifdef GridRiverLakeFlow
MODULE MOD_Grid_ReservoirUpstream
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   Outflow of the unit catchments directly upstream of a built reservoir, as in
!   CaMa-Flood v4 (cmf_ctrl_damout_mod.F90, UPDATE_INFLOW, cells with I1DAM = 10):
!   the outflow into the reservoir is a kinematic wave of the upstream cell's own
!   channel and floodplain, independent of the water level in the reservoir
!   ("to avoid storage buffer effect", Shin et al., 2019, WRR,
!   doi:10.1029/2019WR025265).
!
!   Why.  A built reservoir keeps its whole volume (GRanD, up to several km3) in
!   the one unit catchment of the dam, and its level is that cell's volume-depth
!   curve evaluated at this volume.  The real reservoir spreads over many cells
!   (Kentucky Lake: ~650 km2 and ~300 km long), the dam cell is 100-150 km2 at
!   0.1 degree and 500-1000 km2 at 0.25 degree, so the model level is far above
!   the real one: at the normal volume min(0.7*total, ConVol) the median dam cell
!   is 16 m deep and the 90th percentile 64 m at 0.1 degree (Kentucky 64 m, Toledo
!   Bend 85 m; the real lakes are ~20 m deep at the dam).  The face between the
!   upstream cell and the reservoir was an ordinary HLL face, so this level was the
!   downstream water height of that face: the outflow of the upstream cell stalled
!   or reversed, the cell filled until its level reached the reservoir level, and
!   the backwater went on to the next cells upstream.  At 0.1 degree 1666 of 5870
!   reservoirs have ordinary cells upstream whose bank top is below their normal
!   level (72150 cells in all; 153 above Kentucky Dam, up to 49 m over the bank);
!   at 0.25 degree 529 reservoirs and 10239 cells.  In a 0.1 degree run over the
!   lower Mississippi this showed up as river cells upstream of reservoirs whose
!   storage rose by about two orders of magnitude while their discharge fell, with
!   depths of tens of metres and flooded fractions near 100%.
!
!   What changes.  Only the outlet face of these cells: its water flux is the
!   kinematic outflow below, the cell's velocity is the kinematic channel velocity
!   (as D2RIVVEL in CaMa), and the momentum bound on dt does not apply to it, since
!   that momentum is not used.  The reservoir still receives all this water and
!   operates as before, and every other face is unchanged.  Runs without
!   reservoirs (DEF_Reservoir_Method = 0), and the years before a dam is built,
!   are unchanged bit for bit.
!
!   Kept from CaMa: minimum slope 1e-5 (PMINSLP), floodplain Manning 0.10
!   (PMANFLD), floodplain slope capped at 0.005, slope from the bank top
!   elevations over the distance to the downstream cell (D2ELEVTN, D2NXTDST;
!   here topo_rivelv + topo_rivhgt and topo_distance), and the floodplain flow
!   area fldsto/rivlen - flddph*rivwth with CaMa's split of the storage
!   (river storage rivstomax + rivlen*rivwth*flddph).  Not kept: the cap
!   outflow <= storage/dt (here dt is chosen after the fluxes, and the storage
!   bound on dt, or the outgoing-flux limiter with bifurcation, keeps the
!   storage non-negative as for every other face), and CaMa's closing of
!   bifurcation paths at cells upstream of a dam (here paths are closed at the
!   reservoir cells only, as before; a path from such a cell uses its own level,
!   which is no longer raised).
!-------------------------------------------------------------------------------------

   USE MOD_Precision
   IMPLICIT NONE

   real(r8), parameter :: KINE_MIN_SLOPE     = 1.e-5_r8  ! PMINSLP in CaMa-Flood
   real(r8), parameter :: KINE_FLD_MAX_SLOPE = 0.005_r8  ! floodplain slope cap in UPDATE_INFLOW
   real(r8), parameter :: KINE_FLD_MANNING   = 0.10_r8   ! PMANFLD in CaMa-Flood
   real(r8), parameter :: KINE_MIN_DEPTH     = 1.e-5_r8  ! = RIVERMIN of the routing

   ! slope to the downstream cell, (bank top - downstream bank top) / distance
   real(r8), allocatable :: kine_slope (:)

   PUBLIC :: reservoir_upstream_init
   PUBLIC :: reservoir_upstream_mark
   PUBLIC :: reservoir_upstream_outflow
   PUBLIC :: reservoir_upstream_final

CONTAINS

   ! -------
   SUBROUTINE reservoir_upstream_init ()

   USE MOD_SPMD_Task
   USE MOD_Namelist,              only: get_unitcatchment_file
   USE MOD_WorkerPushData
   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_next, push_next2ucat, &
      topo_rivelv, topo_rivhgt, topo_rivlen, readin_riverlake_parameter
   IMPLICIT NONE

   real(r8), allocatable :: nxtdst (:), elv (:), elv_next (:)
   integer :: i

      ! collective: the master reads and sends each worker its cells
      CALL readin_riverlake_parameter (get_unitcatchment_file (), 'topo_distance', rdata1d = nxtdst)

      IF (p_is_worker) THEN

         allocate (kine_slope (numucat))
         allocate (elv        (numucat))
         allocate (elv_next   (numucat))

         IF (numucat > 0) elv = topo_rivelv + topo_rivhgt

         CALL worker_push_data (push_next2ucat, elv, elv_next, fillvalue = 0._r8)

         DO i = 1, numucat
            kine_slope(i) = KINE_MIN_SLOPE
            IF (ucat_next(i) > 0) THEN
               IF (nxtdst(i) > 0._r8) THEN
                  kine_slope(i) = max((elv(i) - elv_next(i)) / nxtdst(i), KINE_MIN_SLOPE)
               ELSE
                  kine_slope(i) = max((elv(i) - elv_next(i)) / topo_rivlen(i), KINE_MIN_SLOPE)
               ENDIF
            ENDIF
         ENDDO

         deallocate (elv, elv_next)

      ENDIF

      IF (allocated(nxtdst)) deallocate (nxtdst)

   END SUBROUTINE reservoir_upstream_init

   ! -------
   ! Cells whose downstream cell is a built reservoir (I1DAM = 10 in CaMa).  Called
   ! once per routing call by all workers, as dams start in their build year.
   SUBROUTINE reservoir_upstream_mark (is_built_resv, is_resv_upstream)

   USE MOD_WorkerPushData
   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_next, push_next2ucat
   IMPLICIT NONE

   logical, intent(in)  :: is_built_resv    (:)
   logical, intent(out) :: is_resv_upstream (:)

   real(r8), allocatable :: built (:), built_next (:)
   integer :: i

      allocate (built      (numucat))
      allocate (built_next (numucat))

      DO i = 1, numucat
         built(i) = merge(1._r8, 0._r8, is_built_resv(i))
      ENDDO

      CALL worker_push_data (push_next2ucat, built, built_next, fillvalue = 0._r8)

      DO i = 1, numucat
         is_resv_upstream(i) = (.not. is_built_resv(i)) .and. (ucat_next(i) > 0) &
            .and. (built_next(i) > 0.5_r8)
      ENDDO

      deallocate (built, built_next)

   END SUBROUTINE reservoir_upstream_mark

   ! -------
   ! Kinematic outflow of cell i (CaMa UPDATE_INFLOW, step 1c), from its water depth
   ! wdsrf [m] and visible water volume volwater [m3]; hflux [m3/s] goes to the
   ! reservoir, veloc [m/s] is the channel velocity.
   SUBROUTINE reservoir_upstream_outflow (i, wdsrf, volwater, hflux, veloc)

   USE MOD_Grid_RiverLakeNetwork, only: topo_rivman, topo_rivwth, topo_rivhgt, &
      topo_rivlen, topo_rivstomax
   IMPLICIT NONE

   integer,  intent(in)  :: i
   real(r8), intent(in)  :: wdsrf, volwater
   real(r8), intent(out) :: hflux, veloc

   real(r8) :: flddph, rivsto, fldsto, area_fld, veloc_fld

      IF (wdsrf < KINE_MIN_DEPTH) THEN
         hflux = 0._r8
         veloc = 0._r8
         RETURN
      ENDIF

      ! river channel: Manning velocity with the channel depth, over width * depth
      veloc = topo_rivman(i)**(-1) * sqrt(kine_slope(i)) * wdsrf**(2._r8/3._r8)
      hflux = topo_rivwth(i) * wdsrf * veloc

      ! floodplain, above the bank
      flddph = max(wdsrf - topo_rivhgt(i), 0._r8)
      IF (flddph > 0._r8) THEN
         rivsto    = topo_rivstomax(i) + topo_rivlen(i) * topo_rivwth(i) * flddph
         fldsto    = max(volwater - rivsto, 0._r8)
         area_fld  = max(fldsto / topo_rivlen(i) - flddph * topo_rivwth(i), 0._r8)
         veloc_fld = KINE_FLD_MANNING**(-1) * sqrt(min(KINE_FLD_MAX_SLOPE, kine_slope(i))) &
            * flddph**(2._r8/3._r8)
         hflux     = hflux + area_fld * veloc_fld
      ENDIF

   END SUBROUTINE reservoir_upstream_outflow

   ! -------
   SUBROUTINE reservoir_upstream_final ()

   IMPLICIT NONE

      IF (allocated(kine_slope)) deallocate (kine_slope)

   END SUBROUTINE reservoir_upstream_final

END MODULE MOD_Grid_ReservoirUpstream
#endif
