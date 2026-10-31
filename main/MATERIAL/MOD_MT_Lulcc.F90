#include <define.h>

#if (defined MaterialTransport) && (defined LULCC)
MODULE MOD_MT_Lulcc
!-----------------------------------------------------------------------
! Material transport across a land cover change.
!
! LulccDriver calls two routines, so the host LULCC files keep their
! upstream layout:
!
!  mt_lulcc_save   before LulccInitialize, while landpatch/mesh/pixel still
!                  describe the old year. It snapshots the material pools
!                  and, for the MEC scheme, the physical area of every old
!                  patch.
!
!  mt_lulcc_remap  after the host variables are recovered (SAT or MEC) and
!                  before the transfer trace is freed. It moves the material
!                  pools, the provider pools and the forcing ratios onto the
!                  new patches and reloads provider inputs for the new year.
!
! With the MEC scheme the transfer trace says where each new patch's area
! came from. Before any pool is moved, the trace is checked against the old
! and new patch areas element by element, so a remap that could not conserve
! mass stops the run instead of silently creating or losing material.
!
! Moved here from the fork's MOD_Lulcc_Driver (#ifdef MATERIAL blocks).
!-----------------------------------------------------------------------

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

   PUBLIC :: mt_lulcc_save
   PUBLIC :: mt_lulcc_remap

   PRIVATE

   real(r8), allocatable :: old_patch_area(:)
   logical :: inventory_active = .false.

   ! Rows of the old host water passed through the remap to host_sync.
   integer, parameter :: HW_LDEW = 1, HW_LDEW_RAIN = 2, HW_LDEW_SNOW = 3, &
                         HW_THIN_SNOW = 4, HW_WDSRF = 5, HW_N2D = 5

CONTAINS

   SUBROUTINE mt_lulcc_save ()

   USE MOD_SPMD_Task
   USE MOD_Namelist,  only: DEF_LULCC_SCHEME
   USE MOD_LandPatch, only: landpatch, numpatch
   USE MOD_MT_Defs,   only: nmts
   USE MOD_MT_Vars,   only: save_land_mt_lulcc_state
   USE MOD_MT_Lifecycle, only: mt_lifecycle_land_save_lulcc_state
   IMPLICIT NONE

      CALL save_land_mt_lulcc_state ()
      CALL mt_lifecycle_land_save_lulcc_state ()

      IF (allocated(old_patch_area)) deallocate(old_patch_area)
      inventory_active = p_is_worker .and. nmts > 0 .and. DEF_LULCC_SCHEME == 2
      IF (inventory_active) THEN
         IF (numpatch > 0) THEN
            IF (.not. allocated(landpatch%eindex)) &
               CALL CoLM_stop('MT LULCC old patch map unavailable')
            IF (landpatch%nset /= numpatch) &
               CALL CoLM_stop('MT LULCC old patch count mismatch')
            CALL lulcc_patch_areas(landpatch, old_patch_area)
         ELSE
            allocate(old_patch_area(0))
         ENDIF
      ENDIF

   END SUBROUTINE mt_lulcc_save


   SUBROUTINE mt_lulcc_remap (lc_year, dir_landdata)

   USE MOD_SPMD_Task
   USE MOD_LandPatch, only: landpatch, numpatch
   USE MOD_Vars_TimeInvariants, only: patchclass
   USE MOD_Lulcc_Vars_TimeInvariants, only: landpatch_, patchclass_
   USE MOD_Lulcc_TransferTraceReadin, only: lccpct_patches
   USE MOD_Vars_TimeVariables, only: wa
   USE MOD_Vars_Global, only: maxsnl, nl_soil
   USE MOD_Lulcc_Vars_TimeVariables, only: ldew_, ldew_rain_, ldew_snow_, scv_, wdsrf_, &
                                           wliq_soisno_, wice_soisno_
   USE MOD_MT_Vars,   only: remap_land_mt_lulcc_state
   USE MOD_MT_Lifecycle, only: mt_lifecycle_land_remap_lulcc_state, &
                               mt_lifecycle_land_reload_lulcc_inputs
   USE MOD_MT_Conservation, only: deallocate_mt_conservation
   USE MOD_MT_Forcing, only: mt_forcing_lulcc_remap
   IMPLICIT NONE

   integer, intent(in) :: lc_year
   character(len=*), intent(in) :: dir_landdata

   real(r8), allocatable :: new_patch_area(:), inventory_trace(:,:)
   real(r8), allocatable :: host_water_old_2d(:,:), host_water_old_3d(:,:,:)
   integer :: nold

      IF (inventory_active) THEN
         IF ((size(old_patch_area) > 0) .neqv. (numpatch > 0)) &
            CALL CoLM_stop('MT LULCC worker element footprint changed')
         IF (size(old_patch_area) > 0) THEN
            IF (.not. allocated(patchclass_) .or. .not. allocated(landpatch_%eindex)) &
               CALL CoLM_stop('MT LULCC old patch map unavailable')
         ENDIF
         IF (.not. allocated(lccpct_patches)) &
            CALL CoLM_stop('MT LULCC MEC requires the transfer trace')
         IF (numpatch > 0) THEN
            IF (.not. allocated(landpatch%eindex)) &
               CALL CoLM_stop('MT LULCC new patch map unavailable')
            IF (landpatch%nset /= numpatch) &
               CALL CoLM_stop('MT LULCC new patch count mismatch')
            CALL lulcc_patch_areas(landpatch, new_patch_area)
         ELSE
            allocate(new_patch_area(0))
         ENDIF
         IF (allocated(patchclass_)) THEN
            IF (.not. allocated(patchclass)) &
               CALL CoLM_stop('MT LULCC MEC patch class map unavailable')
         ELSE
            IF (size(old_patch_area) > 0 .or. size(new_patch_area) > 0) &
               CALL CoLM_stop('MT LULCC MEC patch class map unavailable')
         ENDIF
         CALL lulcc_inventory_trace(lccpct_patches, inventory_trace)
         IF (size(old_patch_area) > 0 .or. size(new_patch_area) > 0) &
            CALL lulcc_check_inventory_transfer(patchclass_, landpatch_%eindex, old_patch_area, &
               patchclass, landpatch%eindex, new_patch_area, inventory_trace)
      ENDIF

      IF (p_is_worker .and. allocated(patchclass)) THEN
         IF (size(patchclass) > 0 .and. .not. allocated(patchclass_)) &
            CALL CoLM_stop('MT LULCC cannot remap a worker from zero to nonzero patches')
      ENDIF
      IF (p_is_worker .and. allocated(patchclass) .and. allocated(patchclass_) .and. &
          allocated(landpatch%eindex) .and. allocated(landpatch_%eindex)) THEN
         CALL deallocate_mt_conservation ()
         ! Old host water by compartment (see HW_* and mt_lulcc_host_sync);
         ! the host keeps the old arrays until LulccDriver frees them after
         ! this call.
         nold = size(patchclass_)
         allocate(host_water_old_2d(HW_N2D, nold))
         allocate(host_water_old_3d(2, maxsnl+1:nl_soil, nold))
         host_water_old_2d = 0._r8
         host_water_old_3d = 0._r8
         IF (nold > 0) THEN
            host_water_old_2d(HW_LDEW,      :) = ldew_     (1:nold)
            host_water_old_2d(HW_LDEW_RAIN, :) = ldew_rain_(1:nold)
            host_water_old_2d(HW_LDEW_SNOW, :) = ldew_snow_(1:nold)
            host_water_old_2d(HW_WDSRF,     :) = wdsrf_    (1:nold)
            ! Pre-layer snow: scv holds it while snow layer 0 is empty.
            WHERE (wliq_soisno_(0,1:nold) + wice_soisno_(0,1:nold) <= 0._r8)
               host_water_old_2d(HW_THIN_SNOW, :) = scv_(1:nold)
            END WHERE
            host_water_old_3d(1, :, :) = wliq_soisno_(:, 1:nold)
            host_water_old_3d(2, :, :) = wice_soisno_(:, 1:nold)
         ENDIF
         IF (inventory_active) THEN
            CALL remap_land_mt_lulcc_state (patchclass, landpatch%eindex, &
               patchclass_, landpatch_%eindex, wa, host_water_old_2d, host_water_old_3d, &
               mt_lulcc_host_sync, inventory_trace, new_patch_area, old_patch_area)
            CALL mt_lifecycle_land_remap_lulcc_state (patchclass, landpatch%eindex, &
               patchclass_, landpatch_%eindex, inventory_trace, new_patch_area, old_patch_area)
            CALL mt_forcing_lulcc_remap (patchclass, landpatch%eindex, &
               patchclass_, landpatch_%eindex, lccpct_patches)
         ELSE
            CALL remap_land_mt_lulcc_state (patchclass, landpatch%eindex, &
               patchclass_, landpatch_%eindex, wa, host_water_old_2d, host_water_old_3d, &
               mt_lulcc_host_sync)
            CALL mt_lifecycle_land_remap_lulcc_state (patchclass, landpatch%eindex, &
               patchclass_, landpatch_%eindex)
            CALL mt_forcing_lulcc_remap (patchclass, landpatch%eindex, &
               patchclass_, landpatch_%eindex)
         ENDIF
         deallocate(host_water_old_2d, host_water_old_3d)
      ENDIF

      IF (allocated(old_patch_area))  deallocate(old_patch_area)
      IF (allocated(new_patch_area))  deallocate(new_patch_area)
      IF (allocated(inventory_trace)) deallocate(inventory_trace)
      inventory_active = .false.

      CALL mt_lifecycle_land_reload_lulcc_inputs (lc_year, dir_landdata)

   END SUBROUTINE mt_lulcc_remap


   ! Called by remap_land_mt_lulcc_state after the pools are remapped.
   ! The host LULCC re-initialises every patch for the new year and then
   ! copies back only the variables it saved (REST_LulccTimeVariables,
   ! LulccMassEnergyConserve), so
   !  - a new patch with no source patch keeps the host cold-start water: its
   !    material pools get the material cold start (mt_init_from_water);
   !  - wetland water (wetwat) is not copied back and restarts from its cold
   !    value on every patch: the material in it restarts from the cold
   !    ratio as well;
   !  - elsewhere the host does not always carry its water over the way the
   !    material pools are moved (by area, compartment by compartment):
   !    MEC re-bins the snow layers to the dominant source patch and drops
   !    the pre-layer snow of the other sources, and with PFT or PC the
   !    patch canopy water is rebuilt from the new PFT fractions. The old
   !    host water, remapped with the same weights as the pools
   !    (src_water_2d/3d), shows where. In a compartment where it differs
   !    from the new host water, the material follows the host water at the
   !    ratio it had (material over remapped old water), or at the cold
   !    ratio where no old water was there. The snow column is handled as
   !    one: if any of its layers or its pre-layer snow differs, every snow
   !    pool gets the column ratio. Elsewhere nothing is touched, so the
   !    material profiles are kept wherever the host kept its water.
   ! The material mass this adds or removes (per unit area times patch area)
   ! is added to change for the LULCC mass check.
   SUBROUTINE mt_lulcc_host_sync (no_source, patch_area, change, src_water_2d, src_water_3d)

   USE MOD_SPMD_Task, only: p_iam_glb
   USE MOD_LandPatch, only: numpatch
   USE MOD_Namelist,  only: DEF_VEG_SNOW
   USE MOD_Vars_TimeInvariants, only: patchtype
   USE MOD_MT_Precip, only: mt_canopy_parts
   USE MOD_Vars_Global, only: maxsnl, nl_soil
   USE MOD_Vars_TimeVariables, only: ldew, ldew_rain, ldew_snow, wliq_soisno, wice_soisno, &
                                     wa, wdsrf, wetwat, scv, tleaf
   USE MOD_MT_Defs,   only: nmts, mt_uses_generic_water, material_init_water_ratio, &
                            mt_water_min_for_ratio
   USE MOD_MT_Vars,   only: mt_wetwat, mt_surface_solid, mt_patch_land_water_mass, &
                            mt_ldew_rain, mt_ldew_snow, mt_wliq_soisno, mt_wice_soisno, &
                            mt_scv, mt_wdsrf, mt_solid_soisno, mt_canopy_solid
   USE MOD_MT_Solute, only: mt_equilibrate_dissolved, mt_has_dissolved_limit
   USE MOD_MT_Rest,   only: mt_init_from_water
   IMPLICIT NONE

   logical,  intent(in)    :: no_source(:)
   real(r8), intent(in)    :: patch_area(:)
   real(r8), intent(inout) :: change(:)
   real(r8), intent(in)    :: src_water_2d(:,:)    ! remapped old host water, HW_* rows
   real(r8), intent(in)    :: src_water_3d(:,:,:)  ! (liquid/ice, maxsnl+1:nl_soil, patch)

   real(r8), allocatable :: before(:,:)
   real(r8) :: after(nmts), R_init
   integer  :: ip, imt, ncold, nfollow
   logical  :: followed

      IF (numpatch <= 0 .or. nmts <= 0) RETURN

      allocate(before(nmts, numpatch))
      DO ip = 1, numpatch
         CALL mt_patch_land_water_mass (ip, before(:,ip))
      ENDDO

      ncold = count(no_source(1:numpatch))
      IF (ncold > 0) THEN
         CALL mt_init_from_water (numpatch, maxsnl, nl_soil, &
            ldew_rain, ldew_snow, wliq_soisno, wice_soisno, wa, wdsrf, wetwat, scv, &
            patch_mask = no_source(1:numpatch))
      ENDIF

      nfollow = 0
      DO ip = 1, numpatch
         IF (no_source(ip)) CYCLE
         IF (ip > size(src_water_2d,2) .or. ip > size(src_water_3d,3)) CYCLE
         followed = .false.
         DO imt = 1, nmts
            IF (.not. mt_uses_generic_water(imt)) CYCLE
            R_init = material_init_water_ratio(imt)
            mt_wetwat(imt,ip) = 0._r8
            IF (patchtype(ip) == 2) mt_wetwat(imt,ip) = max(wetwat(ip), 0._r8) * R_init
            CALL follow_host_water (imt, ip, R_init, followed)
         ENDDO
         IF (followed) nfollow = nfollow + 1
      ENDDO

      DO ip = 1, numpatch
         CALL mt_patch_land_water_mass (ip, after)
         IF (ip <= size(patch_area)) change = change + (after - before(:,ip)) * patch_area(ip)
      ENDDO
      deallocate(before)

      IF (ncold > 0) write(*,'(A,I6,A,I8,A)') ' MT LULCC on process', p_iam_glb, ': ', ncold, &
         ' new patches without a source patch cold-started from the host water'
      IF (nfollow > 0) write(*,'(A,I6,A,I8,A)') ' MT LULCC on process', p_iam_glb, ': ', nfollow, &
         ' patches with material set to follow host water the host did not carry over by area'

   CONTAINS

      logical FUNCTION water_differs (a, b)
         real(r8), intent(in) :: a, b
         water_differs = abs(a - b) > mt_water_min_for_ratio + 1.e-12_r8 * max(abs(a), abs(b))
      END FUNCTION water_differs

      ! Material of one compartment at its own ratio, now on the host water.
      SUBROUTINE follow_pool (pool, w_src, w_new, R_init, followed)
         real(r8), intent(inout) :: pool
         real(r8), intent(in)    :: w_src, w_new, R_init
         logical,  intent(inout) :: followed
         IF (.not. water_differs(w_src, w_new)) RETURN
         IF (w_src > mt_water_min_for_ratio) THEN
            pool = pool / w_src * max(w_new, 0._r8)
         ELSE
            pool = R_init * max(w_new, 0._r8)
         ENDIF
         followed = .true.
      END SUBROUTINE follow_pool

      SUBROUTINE follow_host_water (imt, ip, R_init, followed)
         integer,  intent(in)    :: imt, ip
         real(r8), intent(in)    :: R_init
         logical,  intent(inout) :: followed
         integer  :: j, k
         logical  :: snow_differs
         real(r8) :: w_src, ratio, thin_new, solid, wsnow
         real(r8) :: canopy_liquid, canopy_snow, surface_water, surface_mass, wetland_water

         ! Canopy. Without DEF_VEG_SNOW the material splits the total
         ! ldew into its rain and snow pools itself (MOD_MT_Precip), so
         ! the total is what must match.
         IF (DEF_VEG_SNOW) THEN
            CALL follow_pool (mt_ldew_rain(imt,ip), src_water_2d(HW_LDEW_RAIN,ip), &
               ldew_rain(ip), R_init, followed)
            CALL follow_pool (mt_ldew_snow(imt,ip), src_water_2d(HW_LDEW_SNOW,ip), &
               ldew_snow(ip), R_init, followed)
         ELSE
            w_src = src_water_2d(HW_LDEW,ip)
            IF (water_differs(w_src, ldew(ip))) THEN
               IF (w_src > mt_water_min_for_ratio) THEN
                  mt_ldew_rain(imt,ip) = mt_ldew_rain(imt,ip) / w_src * max(ldew(ip), 0._r8)
                  mt_ldew_snow(imt,ip) = mt_ldew_snow(imt,ip) / w_src * max(ldew(ip), 0._r8)
               ELSE
                  mt_ldew_rain(imt,ip) = R_init * max(ldew(ip), 0._r8)
                  mt_ldew_snow(imt,ip) = 0._r8
               ENDIF
               followed = .true.
            ENDIF
         ENDIF

         ! Snow column: snow layers and the pre-layer snow (scv while snow
         ! layer 0 is empty), as in mt_init_from_water.
         thin_new = 0._r8
         IF (wliq_soisno(0,ip) + wice_soisno(0,ip) <= 0._r8) thin_new = scv(ip)
         snow_differs = water_differs(src_water_2d(HW_THIN_SNOW,ip), thin_new)
         DO j = maxsnl+1, 0
            k = j - maxsnl
            snow_differs = snow_differs .or. water_differs(src_water_3d(1,k,ip), wliq_soisno(j,ip)) &
                                        .or. water_differs(src_water_3d(2,k,ip), wice_soisno(j,ip))
         ENDDO
         IF (snow_differs) THEN
            w_src = src_water_2d(HW_THIN_SNOW,ip) + sum(src_water_3d(:,1:-maxsnl,ip))
            IF (w_src > mt_water_min_for_ratio) THEN
               ratio = (mt_scv(imt,ip) + sum(mt_wliq_soisno(imt,maxsnl+1:0,ip)) &
                  + sum(mt_wice_soisno(imt,maxsnl+1:0,ip))) / w_src
            ELSE
               ratio = R_init
            ENDIF
            mt_scv(imt,ip) = ratio * max(thin_new, 0._r8)
            DO j = maxsnl+1, 0
               mt_wliq_soisno(imt,j,ip) = ratio * max(wliq_soisno(j,ip), 0._r8)
               mt_wice_soisno(imt,j,ip) = ratio * max(wice_soisno(j,ip), 0._r8)
            ENDDO
            ! Undissolved material of the snow layers goes with the new snow
            ! layers, or to the top soil layer when no snow layer is left.
            solid = sum(mt_solid_soisno(imt,maxsnl+1:0,ip))
            wsnow = sum(max(wliq_soisno(maxsnl+1:0,ip), 0._r8) + max(wice_soisno(maxsnl+1:0,ip), 0._r8))
            IF (solid /= 0._r8) THEN
               IF (wsnow > 0._r8) THEN
                  DO j = maxsnl+1, 0
                     mt_solid_soisno(imt,j,ip) = solid * (max(wliq_soisno(j,ip), 0._r8) &
                        + max(wice_soisno(j,ip), 0._r8)) / wsnow
                  ENDDO
               ELSE
                  mt_solid_soisno(imt,maxsnl+1:0,ip) = 0._r8
                  mt_solid_soisno(imt,1,ip) = mt_solid_soisno(imt,1,ip) + solid
               ENDIF
            ENDIF
            followed = .true.
         ENDIF

         ! Soil layers, each phase on its own.
         DO j = 1, nl_soil
            k = j - maxsnl
            CALL follow_pool (mt_wliq_soisno(imt,j,ip), src_water_3d(1,k,ip), &
               wliq_soisno(j,ip), R_init, followed)
            CALL follow_pool (mt_wice_soisno(imt,j,ip), src_water_3d(2,k,ip), &
               wice_soisno(j,ip), R_init, followed)
         ENDDO

         ! Surface water.
         CALL follow_pool (mt_wdsrf(imt,ip), src_water_2d(HW_WDSRF,ip), wdsrf(ip), R_init, followed)
         ! Remapping can leave solid material beside an undersaturated liquid
         ! even when its water amount is unchanged. Rebalance each final pool.
         IF (mt_has_dissolved_limit(imt)) THEN
            CALL mt_canopy_parts(ldew(ip), ldew_rain(ip), ldew_snow(ip), tleaf(ip), &
               canopy_liquid, canopy_snow)
            CALL mt_equilibrate_dissolved(imt, max(canopy_liquid, 0._r8), &
               mt_ldew_rain(imt,ip), mt_canopy_solid(imt,ip))
            DO j = maxsnl+1, nl_soil
               CALL mt_equilibrate_dissolved(imt, max(wliq_soisno(j,ip), 0._r8), &
                  mt_wliq_soisno(imt,j,ip), mt_solid_soisno(imt,j,ip))
            ENDDO

            ! Wetland water and ponded water share one surface solid store;
            ! solve their combined capacity once, rather than in pool order.
            wetland_water = 0._r8
            IF (patchtype(ip) == 2) wetland_water = max(wetwat(ip), 0._r8)
            surface_water = max(wdsrf(ip), 0._r8) + wetland_water
            surface_mass = mt_wdsrf(imt,ip) + mt_wetwat(imt,ip)
            CALL mt_equilibrate_dissolved(imt, surface_water, surface_mass, mt_surface_solid(imt,ip))
            IF (surface_water > mt_water_min_for_ratio) THEN
               mt_wdsrf(imt,ip) = surface_mass * max(wdsrf(ip), 0._r8) / surface_water
               mt_wetwat(imt,ip) = surface_mass * wetland_water / surface_water
            ELSE
               mt_wdsrf(imt,ip) = 0._r8
               mt_wetwat(imt,ip) = 0._r8
            ENDIF
         ENDIF
      END SUBROUTINE follow_host_water

   END SUBROUTINE mt_lulcc_host_sync

   ! Physical area (m2) of every patch, from its pixels and shared fraction.
   SUBROUTINE lulcc_patch_areas(patches, areas)

   USE MOD_SPMD_Task, only: CoLM_stop
   USE MOD_Mesh,      only: mesh
   USE MOD_Pixel,     only: pixel
   USE MOD_Pixelset,  only: pixelset_type
   USE MOD_Utils,     only: areaquad
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   IMPLICIT NONE

   TYPE(pixelset_type), intent(in) :: patches
   real(r8), allocatable, intent(out) :: areas(:)
   integer :: ip, ie, ix, first, last

      allocate(areas(patches%nset))
      areas = 0._r8
      IF (patches%nset <= 0) RETURN
      IF (.not. allocated(mesh) .or. .not. allocated(patches%ielm) .or. &
          .not. allocated(patches%ipxstt) .or. .not. allocated(patches%ipxend) .or. &
          .not. allocated(pixel%lat_s) .or. .not. allocated(pixel%lat_n) .or. &
          .not. allocated(pixel%lon_w) .or. .not. allocated(pixel%lon_e)) &
         CALL CoLM_stop('MT LULCC patch geometry unavailable')
      IF (size(patches%ielm) /= patches%nset .or. size(patches%ipxstt) /= patches%nset .or. &
          size(patches%ipxend) /= patches%nset) &
         CALL CoLM_stop('MT LULCC patch geometry shape mismatch')
      IF (patches%has_shared) THEN
         IF (.not. allocated(patches%pctshared)) &
            CALL CoLM_stop('MT LULCC missing shared patch fractions')
         IF (size(patches%pctshared) /= patches%nset) &
            CALL CoLM_stop('MT LULCC shared patch area shape mismatch')
      ENDIF
      DO ip = 1, patches%nset
         ie = patches%ielm(ip)
         IF (ie < 1 .or. ie > size(mesh)) CALL CoLM_stop('MT LULCC invalid patch element')
         IF (.not. allocated(mesh(ie)%ilat) .or. .not. allocated(mesh(ie)%ilon)) &
            CALL CoLM_stop('MT LULCC element pixels unavailable')
         IF (size(mesh(ie)%ilat) < mesh(ie)%npxl .or. size(mesh(ie)%ilon) < mesh(ie)%npxl) &
            CALL CoLM_stop('MT LULCC element pixel shape mismatch')
         first = patches%ipxstt(ip)
         last = patches%ipxend(ip)
         IF (first == -1 .and. last == -1) CYCLE
         IF (first < 1 .or. last > mesh(ie)%npxl .or. last < first) &
            CALL CoLM_stop('MT LULCC invalid patch pixel range')
         DO ix = first, last
            IF (mesh(ie)%ilat(ix) < 1 .or. mesh(ie)%ilat(ix) > size(pixel%lat_s) .or. &
                mesh(ie)%ilat(ix) > size(pixel%lat_n) .or. mesh(ie)%ilon(ix) < 1 .or. &
                mesh(ie)%ilon(ix) > size(pixel%lon_w) .or. &
                mesh(ie)%ilon(ix) > size(pixel%lon_e)) &
               CALL CoLM_stop('MT LULCC invalid pixel coordinate index')
            areas(ip) = areas(ip) + 1.e6_r8 * areaquad( &
               pixel%lat_s(mesh(ie)%ilat(ix)), pixel%lat_n(mesh(ie)%ilat(ix)), &
               pixel%lon_w(mesh(ie)%ilon(ix)), pixel%lon_e(mesh(ie)%ilon(ix)))
         ENDDO
         IF (patches%has_shared) THEN
            IF (.not. ieee_is_finite(patches%pctshared(ip))) &
               CALL CoLM_stop('MT LULCC invalid shared patch fraction')
            IF (patches%pctshared(ip) < 0._r8 .or. patches%pctshared(ip) > 1._r8) &
               CALL CoLM_stop('MT LULCC invalid shared patch fraction')
            areas(ip) = areas(ip) * patches%pctshared(ip)
         ENDIF
         IF (.not. ieee_is_finite(areas(ip))) &
            CALL CoLM_stop('MT LULCC invalid physical patch area')
         IF (areas(ip) < 0._r8) &
            CALL CoLM_stop('MT LULCC invalid physical patch area')
      ENDDO

   END SUBROUTINE lulcc_patch_areas


   ! The transfer trace is read per land cover class. With PFT/PC all
   ! vegetated classes share one soil patch (class 1 here), and with the fast
   ! PC option cropland and class 14 share the crop patch, so fold the trace
   ! onto the classes the patches actually carry.
   SUBROUTINE lulcc_inventory_trace(raw, mapped)

   USE MOD_SPMD_Task, only: CoLM_stop
   USE MOD_Const_LC,  only: patchtypes
   USE MOD_Vars_Global, only: N_land_classification, CROPLAND
   USE MOD_Namelist,  only: DEF_USE_PFT, DEF_SOLO_PFT, DEF_FAST_PC
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   IMPLICIT NONE

   real(r8), intent(in) :: raw(:,0:)
   real(r8), allocatable, intent(out) :: mapped(:,:)
   integer :: c, dest

      IF (ubound(raw,2) /= N_land_classification) &
         CALL CoLM_stop('MT LULCC transfer trace class count mismatch')
      IF (any(.not. ieee_is_finite(raw))) &
         CALL CoLM_stop('MT LULCC invalid transfer trace')
      IF (any(raw < 0._r8)) &
         CALL CoLM_stop('MT LULCC invalid transfer trace')
      allocate(mapped(size(raw,1),0:N_land_classification))
      mapped = raw
      IF (.not. ((DEF_USE_PFT .and. .not. DEF_SOLO_PFT) .or. DEF_FAST_PC)) RETURN
      mapped = 0._r8
      mapped(:,0) = raw(:,0)
      DO c = 1, N_land_classification
         dest = c
         IF (patchtypes(c) == 0) dest = 1
         IF (DEF_FAST_PC .and. (c == CROPLAND .or. c == 14)) dest = CROPLAND
         mapped(:,dest) = mapped(:,dest) + raw(:,c)
      ENDDO

   END SUBROUTINE lulcc_inventory_trace


   ! Per element: the old area of each class must equal the new area that the
   ! trace says came from that class, and the element's total area must not
   ! change.
   SUBROUTINE lulcc_check_inventory_transfer(old_class, old_element, old_area, &
      new_class, new_element, new_area, trace)

   USE MOD_SPMD_Task, only: CoLM_stop
   USE MOD_Utils,     only: quicksort
   USE MOD_Vars_Global, only: N_land_classification
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   IMPLICIT NONE

   integer, intent(in) :: old_class(:), new_class(:)
   integer*8, intent(in) :: old_element(:), new_element(:)
   real(r8), intent(in) :: old_area(:), new_area(:), trace(:,0:)
   real(r8) :: source(0:N_land_classification), target(0:N_land_classification)
   real(r8) :: old_total, new_total, tol, rowsum
   integer*8, allocatable :: old_sorted(:), new_sorted(:)
   integer, allocatable :: old_order(:), new_order(:)
   integer :: os, oe, ns, ne, op, np, c, k

      IF (size(old_class) /= size(old_area) .or. size(old_class) /= size(old_element) .or. &
          size(new_class) /= size(new_area) .or. size(new_class) /= size(new_element) .or. &
          size(trace,1) /= size(new_class) .or. ubound(trace,2) /= N_land_classification) &
         CALL CoLM_stop('MT LULCC inventory map shape mismatch')
      IF (any(.not. ieee_is_finite(old_area)) .or. any(.not. ieee_is_finite(new_area)) .or. &
          any(.not. ieee_is_finite(trace))) &
         CALL CoLM_stop('MT LULCC invalid inventory area or transfer trace')
      IF (any(old_area < 0._r8) .or. any(new_area < 0._r8) .or. any(trace < 0._r8)) &
         CALL CoLM_stop('MT LULCC invalid inventory area or transfer trace')
      old_sorted = old_element
      new_sorted = new_element
      old_order = [(k, k=1,size(old_class))]
      new_order = [(k, k=1,size(new_class))]
      CALL quicksort(size(old_class), old_sorted, old_order)
      CALL quicksort(size(new_class), new_sorted, new_order)
      os = 1
      ns = 1
      DO WHILE (os <= size(old_class) .or. ns <= size(new_class))
         IF (os > size(old_class) .or. ns > size(new_class)) &
            CALL CoLM_stop('MT LULCC element footprint changed')
         IF (old_sorted(os) /= new_sorted(ns)) &
            CALL CoLM_stop('MT LULCC element footprint changed')
         oe = os
         DO WHILE (oe < size(old_class))
            IF (old_sorted(oe+1) /= old_sorted(os)) EXIT
            oe = oe + 1
         ENDDO
         ne = ns
         DO WHILE (ne < size(new_class))
            IF (new_sorted(ne+1) /= new_sorted(ns)) EXIT
            ne = ne + 1
         ENDDO
         source = 0._r8
         target = 0._r8
         DO k = os, oe
            op = old_order(k)
            IF (old_area(op) <= tiny(1._r8)) CYCLE
            c = old_class(op)
            IF (c < 0 .or. c > N_land_classification) &
               CALL CoLM_stop('MT LULCC invalid old patch class')
            source(c) = source(c) + old_area(op)
         ENDDO
         DO k = ns, ne
            np = new_order(k)
            IF (new_area(np) <= tiny(1._r8)) CYCLE
            c = new_class(np)
            IF (c < 0 .or. c > N_land_classification) &
               CALL CoLM_stop('MT LULCC invalid new patch class')
            rowsum = sum(trace(np,:))
            IF (abs(rowsum - 1._r8) > 1.e-10_r8) &
               CALL CoLM_stop('MT LULCC incomplete source trace')
            target = target + new_area(np) * trace(np,:) / rowsum
         ENDDO
         old_total = sum(source)
         new_total = sum(new_area(new_order(ns:ne)))
         tol = 1.e-10_r8 * max(old_total, new_total, 1._r8)
         IF (abs(old_total - new_total) > tol) THEN
            WRITE(*,'(A,I0,2(A,ES16.8))') 'MT LULCC element=', old_sorted(os), &
               ' old_area_m2=', old_total, ' new_area_m2=', new_total
            CALL CoLM_stop('MT LULCC element physical footprint changed')
         ENDIF
         DO c = 0, N_land_classification
            IF (((source(c) > 0._r8) .neqv. (target(c) > 0._r8)) .or. &
                abs(source(c) - target(c)) > 1.e-10_r8 * max(source(c), target(c))) THEN
               WRITE(*,'(A,I0,A,I0,2(A,ES16.8))') 'MT LULCC element=', old_sorted(os), &
                  ' class=', c, ' source_area_m2=', source(c), ' target_area_m2=', target(c)
               IF (c == 0 .and. source(c) <= 0._r8) &
                  WRITE(*,'(A)') '  hint: a transfer trace that puts all area in class 0 ' // &
                     'means the previous-year land cover was never read; ' // &
                     'mksrfdata MAKE_LulccTransferTrace needs the 5x5 RG blocks of the ' // &
                     'landcover dataset and silently fills zeros when they are absent.'
               CALL CoLM_stop('MT LULCC source-class physical area mismatch')
            ENDIF
         ENDDO
         os = oe + 1
         ns = ne + 1
      ENDDO

   END SUBROUTINE lulcc_check_inventory_transfer

END MODULE MOD_MT_Lulcc
#endif
