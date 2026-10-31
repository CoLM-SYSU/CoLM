#include <define.h>

#ifdef GridRiverLakeFlow
MODULE MOD_Grid_RiverLakeHist
!--------------------------------------------------------------------------------
! DESCRIPTION:
!
!     Write out model results in lateral hydrological processes to history files.
!
! Created by Shupeng Zhang, May 2023
!--------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Grid_RiverLakeHistState
#ifdef MaterialTransport
   USE MOD_MT_RiverLakeCoupling, only: mt_rlc_write_history, mt_rlc_flush_history
   PRIVATE :: mt_rlc_write_history, mt_rlc_flush_history
#endif

   ! -- PUBLIC SUBROUTINEs --
   PUBLIC :: hist_grid_riverlake_init
   PUBLIC :: hist_grid_riverlake_out
   PUBLIC :: hist_grid_riverlake_final
   PUBLIC :: flush_acc_fluxes_riverlake

!--------------------------------------------------------------------------
CONTAINS

   !---------------------------------------
   SUBROUTINE hist_grid_riverlake_init (histform)

   USE MOD_Block
   USE MOD_WorkerPushData
   USE MOD_HistGridded
   USE MOD_Grid_RiverLakeNetwork
   USE MOD_Grid_Reservoir,      only: numresv
   USE MOD_LandPatch,           only: numpatch
   USE MOD_Forcing,             only: forcmask_pch
   USE MOD_Vars_TimeInvariants, only: patchtype, patchmask
   USE MOD_Namelist

   IMPLICIT NONE

   character(len=*), intent(in) :: histform

   ! Local Variables
   logical,  allocatable :: filter_basic(:)
   real(r8), allocatable :: vec_ucat(:), vec_grid(:), vec_inpm(:), vec_patch(:)
   integer :: ilon, iblkme, iblk, jblk


      ! ----- allocate memory for accumulative variables -----

      allocate (acctime_ucat (numucat))
      allocate (a_wdsrf_ucat (numucat))
      allocate (a_veloc_riv  (numucat))
      allocate (a_discharge  (numucat))
      allocate (a_floodarea  (numucat))
      allocate (a_rivsto     (numucat))
      allocate (a_fldsto     (numucat))
      allocate (a_flddph     (numucat))
      allocate (a_storge     (numucat))
      allocate (a_sfcelv     (numucat))
      IF (DEF_USE_LEVEE) THEN
         allocate (a_levsto    (numucat))
         allocate (a_levdph    (numucat))
      ENDIF
      IF (DEF_USE_BIFURCATION) THEN
         allocate (a_bifout         (numucat))
         allocate (a_bifflw_lev     (npthlev_bif, npthout_local))
         allocate (a_bifflw_acctime (npthout_local))
      ENDIF

      IF (p_is_worker) THEN
         allocate (a_wdsrf_ucat_pch (numpatch))
         allocate (a_veloc_riv_pch  (numpatch))
         allocate (a_discharge_pch  (numpatch))
         allocate (a_dis_rmth_pch   (numpatch))
         allocate (a_floodfrc_pch   (numpatch))
      ENDIF

      allocate (acctime_resv (numresv))
      allocate (a_volresv    (numresv))
      allocate (a_qresv_in   (numresv))
      allocate (a_qresv_out  (numresv))

      CALL flush_acc_fluxes_riverlake ()

      ! ----- get longitude and latitude -----
      allocate (lon_ucat (griducat%nlon))
      allocate (lat_ucat (griducat%nlat))

      lat_ucat = (griducat%lat_s + griducat%lat_n) * 0.5

      DO ilon = 1, griducat%nlon
         IF (griducat%lon_w(ilon) > griducat%lon_e(ilon)) THEN
            lon_ucat(ilon) = (griducat%lon_w(ilon) + griducat%lon_e(ilon)+360.) * 0.5
            CALL normalize_longitude (lon_ucat(ilon))
         ELSE
            lon_ucat(ilon) = (griducat%lon_w(ilon) + griducat%lon_e(ilon)) * 0.5
         ENDIF
      ENDDO

      ! ----- for auxiliary data -----
      IF (p_is_worker) THEN
         allocate (vec_ucat  (numucat ))
         allocate (vec_grid  (numinpm ))
         allocate (vec_inpm  (numinpm ))
         allocate (vec_patch (numpatch))

         ! Patches excluding (type >= 99), virtual patches and thos forcing missed
         allocate (filter_basic (numpatch))
         IF (numpatch > 0) THEN
            filter_basic = patchtype < 99
            filter_basic = filter_basic .and. patchmask
            IF (DEF_forcing%has_missing_value) THEN
               filter_basic = filter_basic .and. forcmask_pch
            ENDIF
         ENDIF
      ENDIF

      ! ----- 1) area and filter covered by unit catchments -----
      IF (p_is_worker) THEN

         IF (numucat > 0) vec_ucat = 1.

         CALL worker_push_data (push_ucat2grid, vec_ucat, vec_grid, fillvalue = spval)
         CALL worker_remap_data_grid2pset ( &
            remap_patch2inpm, vec_grid, vec_patch, fillvalue = spval, mode = 'average')

         allocate (filter_ucat (numpatch))
         IF (numpatch > 0) THEN
            filter_ucat = filter_basic .and. (vec_patch /= spval)

            WHERE (filter_ucat)
               vec_patch = 1
            ELSE WHERE
               vec_patch = spval
            END WHERE
         ENDIF

         allocate (sum_grid_area (numinpm))
         CALL worker_remap_data_pset2grid (remap_patch2inpm, vec_patch, vec_grid, &
            fillvalue = spval, filter = filter_ucat)
         CALL worker_push_data (allreduce_inpm, vec_grid, sum_grid_area, fillvalue = spval)
      ENDIF

      IF (trim(histform) == 'Gridded') THEN
         IF (p_is_io)  CALL allocate_block_data (ghist, sumarea_ucat)
         CALL mp2g_hist%get_sumarea (sumarea_ucat, filter_ucat)
      ENDIF

      ! ----- 2) area and filter covered by river mouth -----
      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            WHERE (ucat_next == -9)
               vec_ucat = 1.
            ELSE WHERE
               vec_ucat = spval
            END WHERE
         ENDIF

         CALL worker_push_data (push_ucat2grid, vec_ucat, vec_grid, fillvalue = spval)
         CALL worker_remap_data_grid2pset (remap_patch2inpm, vec_grid, vec_patch, &
            fillvalue = spval, mode = 'average')

         allocate (filter_rivmth (numpatch))
         IF (numpatch > 0) THEN
            filter_rivmth = filter_ucat .and. (vec_patch /= spval)

            WHERE (filter_rivmth)
               vec_patch = 1
            ELSE WHERE
               vec_patch = spval
            END WHERE
         ENDIF

         allocate (sum_rmth_area (numinpm))

         CALL worker_remap_data_pset2grid (remap_patch2inpm, vec_patch, vec_grid, &
            fillvalue = spval, filter = filter_rivmth)
         CALL worker_push_data (allreduce_inpm, vec_grid, sum_rmth_area, fillvalue = spval)
      ENDIF

      ! ----- 3) area covered by input matrix -----
      IF (p_is_worker) THEN
         IF (numucat > 0) vec_ucat = 1.
         CALL worker_push_data (push_ucat2inpm, vec_ucat, vec_inpm, &
            fillvalue = spval, mode = 'average')
         CALL worker_remap_data_grid2pset (remap_patch2inpm, vec_inpm, vec_patch, &
            fillvalue = spval, mode = 'average')

         allocate (filter_inpm (numpatch))
         IF (numpatch > 0) THEN
            filter_inpm = filter_basic .and. (vec_patch /= spval)
         ENDIF
      ENDIF

      IF (trim(histform) == 'Gridded') THEN
         IF (p_is_io)  CALL allocate_block_data (ghist, sumarea_inpm)
         CALL mp2g_hist%get_sumarea (sumarea_inpm, filter_inpm)
      ENDIF

      ! ----- 4) mask of unit catchments with all upstreams in simulation region -----
      IF (p_is_worker) THEN
         allocate (allups_mask_pch (numpatch))
      ENDIF

      IF (p_is_worker) THEN
         CALL worker_push_data (push_ucat2grid, allups_mask_ucat, vec_grid, fillvalue = spval)
         CALL worker_remap_data_grid2pset ( &
            remap_patch2inpm, vec_grid, allups_mask_pch, fillvalue = spval, mode = 'average')
      ENDIF

      IF (trim(histform) == 'Gridded') THEN

         IF (p_is_io) CALL allocate_block_data (ghist, allups_mask_grid)

         CALL mp2g_hist%pset2grid (allups_mask_pch, allups_mask_grid, spv = spval, msk = filter_ucat)

         IF (p_is_io) THEN
            DO iblkme = 1, gblock%nblkme
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)
               WHERE (sumarea_ucat%blk(iblk,jblk)%val > 0.)
                  allups_mask_grid%blk(iblk,jblk)%val = &
                     allups_mask_grid%blk(iblk,jblk)%val / sumarea_ucat%blk(iblk,jblk)%val
               ELSE WHERE
                  allups_mask_grid%blk(iblk,jblk)%val = spval
               END WHERE
            ENDDO
         ENDIF

      ENDIF

      IF (allocated (vec_ucat    )) deallocate (vec_ucat    )
      IF (allocated (vec_grid    )) deallocate (vec_grid    )
      IF (allocated (vec_inpm    )) deallocate (vec_inpm    )
      IF (allocated (vec_patch   )) deallocate (vec_patch   )
      IF (allocated (filter_basic)) deallocate (filter_basic)

   END SUBROUTINE hist_grid_riverlake_init

   !---------------------------------------
   SUBROUTINE hist_grid_riverlake_out (file_hist, histform, idate, itime_in_file, is_first_in_file)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_WorkerPushData
   USE MOD_Grid_RiverLakeNetwork
   USE MOD_Grid_Reservoir
   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeHistRoute
   USE MOD_HistGridded
#ifdef UNSTRUCTURED
   USE MOD_HistVector
#endif
   USE MOD_LandPatch,   only: numpatch
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: histform
   integer, intent(in) :: idate(3)
   integer, intent(in) :: itime_in_file
   logical, intent(in) :: is_first_in_file

   ! Local variables
   character(len=256) :: file_hist_ucat
   integer :: itime_in_file_ucat, i, ncol_local_bif
   integer, allocatable :: pth_global_id_bif(:)

   real(r8), allocatable :: acc_vec_grid    (:)
   real(r8), allocatable :: bifflw_local    (:,:)
   real(r8), allocatable :: a_floodfrc_ucat (:)  ! flooded area fraction
   real(r8), allocatable :: levsto_local    (:)
   real(r8), allocatable :: levdph_local    (:)
   real(r8), allocatable :: bifout_local    (:)
   real(r8), allocatable :: volresv_local   (:)
   real(r8), allocatable :: qresv_in_local  (:)
   real(r8), allocatable :: qresv_out_local (:)
   real(r8), allocatable :: a_floodfrc_inpm (:)  ! flooded area fraction
   real(r8) :: window_seconds

      i = len_trim (file_hist)
      DO WHILE (file_hist(i:i) /= '_')
         i = i - 1
      ENDDO
      file_hist_ucat = file_hist(1:i) // 'unitcat_' // file_hist(i+1:)

      window_seconds=0._r8
      IF (p_is_worker .and. size(acctime_ucat)>0) window_seconds=maxval(acctime_ucat)
#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, window_seconds, 1, MPI_REAL8, MPI_MAX, p_comm_glb, p_err)
#endif
      CALL route_hist_begin (file_hist_ucat, idate, is_first_in_file, &
         lon_ucat, lat_ucat, itime_in_file_ucat, window_seconds)

      IF (is_first_in_file) THEN
         IF (trim(histform) == 'Gridded') THEN
            CALL hist_write_var_real8_2d ( &
               file_hist, 'mask_complete_upstream_regird', ghist, -1, allups_mask_grid, compress = 1,   &
               longname = 'Mask of grids with all upstream located in simulation region', units = '100%')
#ifdef UNSTRUCTURED
         ELSE
            CALL aggregate_to_vector_and_write_2d ( &
               allups_mask_pch, file_hist, 'mask_complete_upstream_regird', -1, filter_ucat,            &
               longname = 'Mask of grids with all upstream located in simulation region', units = '100%')
#endif
         ENDIF

         CALL route_hist_write_ucat (allups_mask_ucat, 'mask_complete_upstream', &
            longname = 'Mask of grids with all upstream located in simulation region', &
            units = '100%', &
            no_time = .true.)
      ENDIF

      IF (p_is_worker) THEN
         allocate (acc_vec_grid (numinpm ))
      ENDIF

      IF (DEF_hist_vars%riv_height) THEN
         IF (p_is_worker) THEN
            IF (numucat > 0) THEN
               WHERE (acctime_ucat > 0._r8)
                  a_wdsrf_ucat = a_wdsrf_ucat / acctime_ucat
               ELSE WHERE
                  a_wdsrf_ucat = spval
               END WHERE
            ENDIF
            CALL worker_push_data (push_ucat2grid, a_wdsrf_ucat, acc_vec_grid, fillvalue = spval)
            CALL worker_remap_data_grid2pset ( remap_patch2inpm, acc_vec_grid, a_wdsrf_ucat_pch, &
               fillvalue = spval, mode = 'average')
         ENDIF

         CALL route_hist_write_ucat (a_wdsrf_ucat, 'f_wdpth_ucat', &
            longname = 'deepest water depth in river and flood plain', &
            units = 'm')
      ENDIF

      IF (DEF_hist_vars%riv_veloct) THEN
         IF (p_is_worker) THEN
            IF (numucat > 0) THEN
               WHERE (acctime_ucat > 0._r8)
                  a_veloc_riv = a_veloc_riv / acctime_ucat
               ELSE WHERE
                  a_veloc_riv = spval
               END WHERE
            ENDIF
            CALL worker_push_data (push_ucat2grid, a_veloc_riv, acc_vec_grid, fillvalue = spval)
            CALL worker_remap_data_grid2pset (remap_patch2inpm, acc_vec_grid, a_veloc_riv_pch, &
               fillvalue = spval, mode = 'average')
         ENDIF

         CALL route_hist_write_ucat (a_veloc_riv, 'f_veloc_riv', &
            longname = 'water velocity in river', &
            units = 'm/s')
      ENDIF

      IF (DEF_hist_vars%discharge) THEN
         IF (p_is_worker) THEN
            IF (numucat > 0) THEN
               WHERE (acctime_ucat > 0._r8)
                  a_discharge = a_discharge / acctime_ucat
               ELSE WHERE
                  a_discharge = spval
               END WHERE
            ENDIF
            CALL worker_push_data (push_ucat2grid, a_discharge, acc_vec_grid, fillvalue = spval)

            IF (numinpm > 0)  THEN
               WHERE ((sum_grid_area /= spval) .and. (acc_vec_grid /= spval))
                  acc_vec_grid = acc_vec_grid / sum_grid_area
               ELSE WHERE
                  acc_vec_grid = spval
               END WHERE
            ENDIF

            CALL worker_remap_data_grid2pset (remap_patch2inpm, acc_vec_grid, a_discharge_pch, &
               fillvalue = spval, mode = 'sum')
         ENDIF

         CALL route_hist_write_ucat (a_discharge, 'f_discharge', &
            longname = 'discharge in river and flood plain', &
            units = 'm^3/s')

         IF (p_is_worker) THEN
            IF (numucat > 0)  THEN
               WHERE (ucat_next /= -9) a_discharge = spval
            ENDIF
            CALL worker_push_data (push_ucat2grid, a_discharge, acc_vec_grid, fillvalue = spval)

            IF (numinpm > 0)  THEN
               WHERE ((sum_rmth_area /= spval) .and. (acc_vec_grid /= spval))
                  acc_vec_grid = acc_vec_grid / sum_rmth_area
               ELSE WHERE
                  acc_vec_grid = spval
               END WHERE
            ENDIF

            CALL worker_remap_data_grid2pset (remap_patch2inpm, acc_vec_grid, a_dis_rmth_pch, &
               fillvalue = spval, mode = 'sum')
         ENDIF

         CALL route_hist_write_ucat (a_discharge, 'f_discharge_rivermouth', &
            longname = 'river mouth discharge into ocean', &
            units = 'm^3/s')
      ENDIF

      IF (DEF_hist_vars%floodfrc) THEN

         IF (p_is_worker) THEN

            IF (numucat > 0) THEN
               WHERE (acctime_ucat > 0._r8)
                  a_floodarea = a_floodarea / acctime_ucat
               ELSE WHERE
                  a_floodarea = spval
               END WHERE
            ENDIF

            IF (numucat > 0) THEN
               allocate (a_floodfrc_ucat (numucat))
               WHERE (topo_area > 0 .and. a_floodarea /= spval)
                  a_floodfrc_ucat = a_floodarea / topo_area
               ELSE WHERE
                  a_floodfrc_ucat = spval
               END WHERE
            ELSE
               allocate (a_floodfrc_ucat (0))
            ENDIF

            allocate (a_floodfrc_inpm (numinpm))

            CALL worker_push_data (push_ucat2inpm, a_floodfrc_ucat, a_floodfrc_inpm, &
               fillvalue = spval, mode = 'average')

            CALL worker_remap_data_grid2pset (remap_patch2inpm, a_floodfrc_inpm, a_floodfrc_pch, &
               fillvalue = spval, mode = 'average')
         ELSE
            allocate (a_floodfrc_ucat (0))
         ENDIF

         CALL route_hist_write_ucat (a_floodarea, 'f_floodarea', &
            longname = 'flooded area', &
            units = 'm^2')

         CALL route_hist_write_ucat (a_floodfrc_ucat, 'f_floodfrc', &
            longname = 'flooded area fraction', &
            units = '-')

      ENDIF

      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            WHERE (acctime_ucat > 0.)
               a_rivsto = a_rivsto / acctime_ucat
               a_fldsto = a_fldsto / acctime_ucat
               a_flddph = a_flddph / acctime_ucat
               a_storge = a_storge / acctime_ucat
               a_sfcelv = a_sfcelv / acctime_ucat
            ELSE WHERE
               a_rivsto = spval
               a_fldsto = spval
               a_flddph = spval
               a_storge = spval
               a_sfcelv = spval
            END WHERE
         ENDIF
      ENDIF

      CALL route_hist_write_ucat (a_rivsto, 'f_rivsto', &
         longname = 'below-bank river channel storage', &
         units = 'm^3')

      CALL route_hist_write_ucat (a_fldsto, 'f_fldsto', &
         longname = 'visible overbank storage excluding levee-protected storage', &
         units = 'm^3')

      CALL route_hist_write_ucat (a_flddph, 'f_flddph', &
         longname = 'visible river-side floodplain water depth excluding levee-protected depth', &
         units = 'm')

      CALL route_hist_write_ucat (a_storge, 'f_storge', &
         longname = 'total water storage (river+floodplain+levee)', &
         units = 'm^3')

      CALL route_hist_write_ucat (a_sfcelv, 'f_sfcelv', &
         longname = 'water surface elevation', &
         units = 'm')

      IF (DEF_USE_LEVEE) THEN
         IF (p_is_worker .and. numucat > 0 .and. allocated(a_levsto)) THEN
            allocate (levsto_local (numucat))
            allocate (levdph_local (numucat))
            WHERE (acctime_ucat > 0.)
               levsto_local = a_levsto / acctime_ucat
               levdph_local = a_levdph / acctime_ucat
            ELSE WHERE
               levsto_local = spval
               levdph_local = spval
            END WHERE
         ELSE
            allocate (levsto_local (0))
            allocate (levdph_local (0))
         ENDIF

         CALL route_hist_write_ucat (levsto_local, 'f_levsto', &
            longname = 'water storage in levee-protected area', &
            units = 'm^3')

         CALL route_hist_write_ucat (levdph_local, 'f_levdph', &
            longname = 'water depth in levee-protected area', &
            units = 'm')

         deallocate (levsto_local)
         deallocate (levdph_local)
      ENDIF

      IF (DEF_USE_BIFURCATION) THEN
         IF (p_is_worker .and. numucat > 0 .and. allocated(a_bifout)) THEN
            allocate (bifout_local (numucat))
            WHERE (acctime_ucat > 0.)
               bifout_local = a_bifout / acctime_ucat
            ELSE WHERE
               bifout_local = spval
            END WHERE
         ELSE
            allocate (bifout_local (0))
         ENDIF

         CALL route_hist_write_ucat (bifout_local, 'f_bifout', &
            longname = 'net bifurcation outflow', &
            units = 'm^3/s')

         deallocate (bifout_local)

         IF (npthlev_bif > 0 .and. totalnpthout > 0) THEN
            IF (p_is_worker) THEN
               ncol_local_bif = npthout_local
               IF (allocated(a_bifflw_lev) .and. allocated(a_bifflw_acctime)) THEN
                  allocate (bifflw_local (npthlev_bif, npthout_local))
                  DO i = 1, npthout_local
                     IF (a_bifflw_acctime(i) > 0._r8) THEN
                        bifflw_local(:, i) = a_bifflw_lev(:, i) / a_bifflw_acctime(i)
                     ELSE
                        bifflw_local(:, i) = spval
                     ENDIF
                  ENDDO
               ELSE
                  allocate (bifflw_local (npthlev_bif, npthout_local))
                  bifflw_local = spval
               ENDIF
               allocate (pth_global_id_bif (ncol_local_bif))
               pth_global_id_bif(:) = pth_global_id(:)
            ELSE
               ncol_local_bif = 0
               allocate (bifflw_local      (npthlev_bif, 0))
               allocate (pth_global_id_bif (0))
            ENDIF

            CALL route_hist_write_bif_matrix ( &
               bifflw_local, npthlev_bif, ncol_local_bif, pth_global_id_bif, totalnpthout, &
               'f_bifflw_lev', &
               longname = 'effective bifurcation pathway-layer flow', &
               units = 'm^3/s')

            deallocate (bifflw_local)
            deallocate (pth_global_id_bif)
         ENDIF
      ENDIF

      IF (allocated (a_floodfrc_ucat)) deallocate (a_floodfrc_ucat)
      IF (allocated (a_floodfrc_inpm)) deallocate (a_floodfrc_inpm)
      IF (allocated (acc_vec_grid   )) deallocate (acc_vec_grid   )

#ifdef MaterialTransport
      CALL mt_rlc_write_history (file_hist_ucat, itime_in_file_ucat, acctime_ucat)
#endif

      ! ----- reservoir variables -----
      IF (DEF_Reservoir_Method > 0) THEN
         IF (totalnumresv > 0) THEN

            IF (DEF_hist_vars%volresv) THEN
               IF (p_is_worker .and. numresv > 0) THEN
                  allocate (volresv_local (numresv))
                  WHERE (acctime_resv > 0)
                     volresv_local = a_volresv / acctime_resv
                  ELSEWHERE
                     volresv_local = spval
                  END WHERE
               ELSE
                  allocate (volresv_local (0))
               ENDIF

               CALL route_hist_write_resv (volresv_local, 'volresv', &
                  longname = 'reservoir water volume', units = 'm^3')
               deallocate (volresv_local)
            ENDIF

            IF (DEF_hist_vars%qresv_in) THEN
               IF (p_is_worker .and. numresv > 0) THEN
                  allocate (qresv_in_local (numresv))
                  WHERE (acctime_resv > 0)
                     qresv_in_local = a_qresv_in / acctime_resv
                  ELSEWHERE
                     qresv_in_local = spval
                  END WHERE
               ELSE
                  allocate (qresv_in_local (0))
               ENDIF

               CALL route_hist_write_resv (qresv_in_local, 'qresv_in', &
                  longname = 'reservoir inflow', units = 'm^3/s')
               deallocate (qresv_in_local)
            ENDIF

            IF (DEF_hist_vars%qresv_out) THEN
               IF (p_is_worker .and. numresv > 0) THEN
                  allocate (qresv_out_local (numresv))
                  WHERE (acctime_resv > 0)
                     qresv_out_local = a_qresv_out / acctime_resv
                  ELSEWHERE
                     qresv_out_local = spval
                  END WHERE
               ELSE
                  allocate (qresv_out_local (0))
               ENDIF

               CALL route_hist_write_resv (qresv_out_local, 'qresv_out', &
                  longname = 'reservoir outflow', units = 'm^3/s')
               deallocate (qresv_out_local)
            ENDIF

         ENDIF
      ENDIF

      CALL route_hist_end ()

      CALL flush_acc_fluxes_riverlake ()

   END SUBROUTINE hist_grid_riverlake_out

   !-----------------------
   SUBROUTINE flush_acc_fluxes_riverlake ()

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   USE MOD_Grid_Reservoir,        only: numresv
   IMPLICIT NONE

      IF (allocated(acctime_ucat)) acctime_ucat = 0.
      IF (allocated(a_wdsrf_ucat)) a_wdsrf_ucat = 0.
      IF (allocated(a_veloc_riv )) a_veloc_riv  = 0.
      IF (allocated(a_discharge )) a_discharge  = 0.
      IF (allocated(a_floodarea )) a_floodarea  = 0.
      IF (allocated(a_rivsto    )) a_rivsto     = 0.
      IF (allocated(a_fldsto    )) a_fldsto     = 0.
      IF (allocated(a_flddph    )) a_flddph     = 0.
      IF (allocated(a_storge    )) a_storge     = 0.
      IF (allocated(a_sfcelv    )) a_sfcelv     = 0.
      IF (allocated(a_levsto    )) a_levsto     = 0.
      IF (allocated(a_levdph    )) a_levdph     = 0.
      IF (allocated(a_bifout    )) a_bifout     = 0.
      IF (allocated(a_bifflw_lev)) a_bifflw_lev = 0.
      IF (allocated(a_bifflw_acctime)) a_bifflw_acctime = 0.

      IF (p_is_worker) THEN

         IF (numresv > 0) THEN
            acctime_resv (:) = 0.
            a_volresv    (:) = 0.
            a_qresv_in   (:) = 0.
            a_qresv_out  (:) = 0.
         ENDIF

#ifdef MaterialTransport
         CALL mt_rlc_flush_history ()
#endif

      ENDIF

   END SUBROUTINE flush_acc_fluxes_riverlake

   !---------------------------------------
   SUBROUTINE hist_grid_riverlake_final ()

   USE MOD_DataType, only: block_data_real8_2d_free_mem
   USE MOD_Grid_RiverLakeHistRoute, only: route_hist_final
   IMPLICIT NONE

      CALL route_hist_final ()

      CALL block_data_real8_2d_free_mem (sumarea_ucat)
      CALL block_data_real8_2d_free_mem (sumarea_inpm)
      CALL block_data_real8_2d_free_mem (allups_mask_grid)

      IF (allocated(acctime_ucat    )) deallocate (acctime_ucat    )
      IF (allocated(a_wdsrf_ucat    )) deallocate (a_wdsrf_ucat    )
      IF (allocated(a_veloc_riv     )) deallocate (a_veloc_riv     )
      IF (allocated(a_discharge     )) deallocate (a_discharge     )
      IF (allocated(a_floodarea     )) deallocate (a_floodarea     )
      IF (allocated(a_rivsto        )) deallocate (a_rivsto        )
      IF (allocated(a_fldsto        )) deallocate (a_fldsto        )
      IF (allocated(a_flddph        )) deallocate (a_flddph        )
      IF (allocated(a_storge        )) deallocate (a_storge        )
      IF (allocated(a_sfcelv        )) deallocate (a_sfcelv        )
      IF (allocated(a_levsto        )) deallocate (a_levsto        )
      IF (allocated(a_levdph        )) deallocate (a_levdph        )
      IF (allocated(a_bifout        )) deallocate (a_bifout        )
      IF (allocated(a_bifflw_lev    )) deallocate (a_bifflw_lev    )
      IF (allocated(a_bifflw_acctime)) deallocate (a_bifflw_acctime)

      IF (allocated(a_wdsrf_ucat_pch)) deallocate (a_wdsrf_ucat_pch)
      IF (allocated(a_veloc_riv_pch )) deallocate (a_veloc_riv_pch )
      IF (allocated(a_discharge_pch )) deallocate (a_discharge_pch )
      IF (allocated(a_dis_rmth_pch  )) deallocate (a_dis_rmth_pch  )
      IF (allocated(a_floodfrc_pch  )) deallocate (a_floodfrc_pch  )

      IF (allocated(acctime_resv    )) deallocate (acctime_resv    )
      IF (allocated(a_volresv       )) deallocate (a_volresv       )
      IF (allocated(a_qresv_in      )) deallocate (a_qresv_in      )
      IF (allocated(a_qresv_out     )) deallocate (a_qresv_out     )

      IF (allocated(lon_ucat        )) deallocate (lon_ucat        )
      IF (allocated(lat_ucat        )) deallocate (lat_ucat        )

      IF (allocated(filter_ucat     )) deallocate (filter_ucat     )
      IF (allocated(sum_grid_area   )) deallocate (sum_grid_area   )
      IF (allocated(filter_rivmth   )) deallocate (filter_rivmth   )
      IF (allocated(sum_rmth_area   )) deallocate (sum_rmth_area   )
      IF (allocated(filter_inpm     )) deallocate (filter_inpm     )

      IF (allocated(allups_mask_pch )) deallocate (allups_mask_pch )

   END SUBROUTINE hist_grid_riverlake_final

END MODULE MOD_Grid_RiverLakeHist

#endif
