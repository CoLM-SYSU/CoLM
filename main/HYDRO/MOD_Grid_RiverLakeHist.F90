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
   USE MOD_DataType
#ifdef GridRiverLakeSediment
   USE MOD_Grid_RiverLakeSediment, only: nsed, sed_hist_acctime, &
      a_sedcon, a_sedout, a_bedout, a_sedinp, a_netflw, a_layer, a_shearvel
#endif

   ! -- ACC Fluxes --
   real(r8), allocatable :: acctime_ucat     (:)

   real(r8), allocatable :: a_wdsrf_ucat     (:)
   real(r8), allocatable :: a_veloc_riv      (:)
   real(r8), allocatable :: a_discharge      (:)
   real(r8), allocatable :: a_floodarea      (:)  ! flooded area [m^2]

   real(r8), allocatable :: a_wdsrf_ucat_pch (:)
   real(r8), allocatable :: a_veloc_riv_pch  (:)
   real(r8), allocatable :: a_discharge_pch  (:)
   real(r8), allocatable :: a_dis_rmth_pch   (:)
   real(r8), allocatable :: a_floodfrc_pch   (:)  ! flooded area [m^2]

   ! for reservoirs
   real(r8), allocatable :: acctime_resv     (:)

   real(r8), allocatable :: a_volresv        (:)  ! reservoir water volume [m^3]
   real(r8), allocatable :: a_qresv_in       (:)  ! inflow to reservoir    [m^3/s]
   real(r8), allocatable :: a_qresv_out      (:)  ! outflow from reservoir [m^3/s]

   ! grid information
   real(r8), allocatable :: lon_ucat (:)
   real(r8), allocatable :: lat_ucat (:)

   ! auxiliary data
   type(block_data_real8_2d) :: sumarea_ucat        ! 1) area covered by unit catchments
   logical,  allocatable     :: filter_ucat     (:)
   real(r8), allocatable     :: sum_grid_area   (:) !    sum area of patches inside one grid

   logical,  allocatable     :: filter_rivmth   (:) ! 2) area covered by river mouths
   real(r8), allocatable     :: sum_rmth_area   (:)

   logical,  allocatable     :: filter_inpm     (:) ! 3) area covered by input matrix
   type(block_data_real8_2d) :: sumarea_inpm

   type(block_data_real8_2d) :: allups_mask_grid    ! 4) mask of unit catchments with all
   real(r8), allocatable     :: allups_mask_pch (:) !    upstreams in simulation region

   ! -- PUBLIC SUBROUTINEs --
   PUBLIC :: hist_grid_riverlake_init
   PUBLIC :: hist_grid_riverlake_out
   PUBLIC :: hist_grid_riverlake_final

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

   IMPLICIT NONE

   character(len=*), intent(in) :: histform

   ! Local Variables
   logical,  allocatable :: filter_basic(:)
   real(r8), allocatable :: vec_ucat(:), vec_grid(:), vec_inpm(:), vec_patch(:)
   integer :: ilon, iblkme, iblk, jblk


      ! ----- allocate memory for accumulative variables -----
      IF (p_is_worker) THEN

         IF (numucat > 0) THEN
            allocate (acctime_ucat (numucat))
            allocate (a_wdsrf_ucat (numucat))
            allocate (a_veloc_riv  (numucat))
            allocate (a_discharge  (numucat))
            allocate (a_floodarea  (numucat))
         ENDIF

         IF (numpatch > 0) THEN
            allocate (a_wdsrf_ucat_pch (numpatch))
            allocate (a_veloc_riv_pch  (numpatch))
            allocate (a_discharge_pch  (numpatch))
            allocate (a_dis_rmth_pch   (numpatch))
            allocate (a_floodfrc_pch   (numpatch))
         ENDIF

         IF (numresv > 0) THEN
            allocate (acctime_resv (numresv))
            allocate (a_volresv    (numresv))
            allocate (a_qresv_in   (numresv))
            allocate (a_qresv_out  (numresv))
         ENDIF

      ENDIF

      CALL flush_acc_fluxes_riverlake ()

      ! ----- get longitude and latitude -----
      IF (p_is_master) THEN
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
      ENDIF

      ! ----- for auxiliary data -----
      IF (p_is_worker) THEN
         IF (numucat  > 0) allocate (vec_ucat  (numucat ))
         IF (numinpm  > 0) allocate (vec_grid  (numinpm ))
         IF (numinpm  > 0) allocate (vec_inpm  (numinpm ))
         IF (numpatch > 0) allocate (vec_patch (numpatch))

         ! Patches excluding (type >= 99), virtual patches and thos forcing missed
         IF (numpatch > 0) THEN
            allocate (filter_basic (numpatch))
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

         IF (numpatch > 0) THEN
            allocate (filter_ucat (numpatch))
            filter_ucat = filter_basic .and. (vec_patch /= spval)

            WHERE (filter_ucat)
               vec_patch = 1
            ELSE WHERE
               vec_patch = spval
            END WHERE
         ENDIF

         IF (numinpm > 0) allocate (sum_grid_area (numinpm))
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

         IF (numpatch > 0) THEN
            allocate (filter_rivmth (numpatch))
            filter_rivmth = filter_ucat .and. (vec_patch /= spval)

            WHERE (filter_rivmth)
               vec_patch = 1
            ELSE WHERE
               vec_patch = spval
            END WHERE
         ENDIF

         IF (numinpm > 0) allocate (sum_rmth_area (numinpm))

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

         IF (numpatch > 0) THEN
            allocate (filter_inpm (numpatch))
            filter_inpm = filter_basic .and. (vec_patch /= spval)
         ENDIF
      ENDIF

      IF (trim(histform) == 'Gridded') THEN
         IF (p_is_io)  CALL allocate_block_data (ghist, sumarea_inpm)
         CALL mp2g_hist%get_sumarea (sumarea_inpm, filter_inpm)
      ENDIF

      ! ----- 4) mask of unit catchments with all upstreams in simulation region -----
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            allocate (allups_mask_pch (numpatch))
         ENDIF
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
   logical :: fexists
   integer :: itime_in_file_ucat, i

   real(r8), allocatable :: acc_vec_grid    (:)
   real(r8), allocatable :: a_floodfrc_ucat (:)  ! flooded area fraction
   real(r8), allocatable :: a_floodfrc_inpm (:)  ! flooded area fraction


      IF (p_is_master) THEN
         i = len_trim (file_hist)
         DO WHILE (file_hist(i:i) /= '_')
            i = i - 1
         ENDDO
         file_hist_ucat = file_hist(1:i) // 'unitcat_' // file_hist(i+1:)

         inquire (file=file_hist_ucat, exist=fexists)
         IF (.not. fexists) THEN

            CALL ncio_create_file (trim(file_hist_ucat))

            CALL ncio_define_dimension (file_hist_ucat, 'time', 0)
            CALL ncio_define_dimension (file_hist_ucat, 'lat_ucat', griducat%nlat)
            CALL ncio_define_dimension (file_hist_ucat, 'lon_ucat', griducat%nlon)

            CALL ncio_write_serial (file_hist_ucat, 'lat_ucat', lat_ucat, 'lat_ucat')
            CALL ncio_write_serial (file_hist_ucat, 'lon_ucat', lon_ucat, 'lon_ucat')
         ENDIF

         CALL ncio_write_time (file_hist_ucat, 'time', idate, itime_in_file_ucat, DEF_HIST_FREQ)

      ENDIF

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

         CALL vector_gather_map2grid_and_write ( &
            allups_mask_ucat, numucat, totalnumucat, ucat_data_address, griducat%nlon, x_ucat,       &
            griducat%nlat, y_ucat, file_hist_ucat, 'mask_complete_upstream', 'lon_ucat', 'lat_ucat', &
            longname = 'Mask of grids with all upstream located in simulation region', units = '100%')
      ENDIF

      IF (p_is_worker) THEN
         IF (numinpm  > 0) allocate (acc_vec_grid (numinpm ))
      ENDIF

      IF (DEF_hist_vars%riv_height) THEN
         IF (p_is_worker) THEN
            IF (numucat > 0)  a_wdsrf_ucat = a_wdsrf_ucat / acctime_ucat
            CALL worker_push_data (push_ucat2grid, a_wdsrf_ucat, acc_vec_grid, fillvalue = spval)
            CALL worker_remap_data_grid2pset ( remap_patch2inpm, acc_vec_grid, a_wdsrf_ucat_pch, &
               fillvalue = spval, mode = 'average')
         ENDIF

         CALL vector_gather_map2grid_and_write ( a_wdsrf_ucat, numucat,                    &
            totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
            file_hist_ucat, 'f_wdpth_ucat', 'lon_ucat', 'lat_ucat', itime_in_file_ucat,    &
            'deepest water depth in river and flood plain', 'm')
      ENDIF

      IF (DEF_hist_vars%riv_veloct) THEN
         IF (p_is_worker) THEN
            IF (numucat > 0)  a_veloc_riv = a_veloc_riv  / acctime_ucat
            CALL worker_push_data (push_ucat2grid, a_veloc_riv, acc_vec_grid, fillvalue = spval)
            CALL worker_remap_data_grid2pset (remap_patch2inpm, acc_vec_grid, a_veloc_riv_pch, &
               fillvalue = spval, mode = 'average')
         ENDIF

         CALL vector_gather_map2grid_and_write ( a_veloc_riv, numucat,                     &
            totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
            file_hist_ucat, 'f_veloc_riv', 'lon_ucat', 'lat_ucat', itime_in_file_ucat,     &
            'water velocity in river', 'm/s')
      ENDIF

      IF (DEF_hist_vars%discharge) THEN
         IF (p_is_worker) THEN
            IF (numucat > 0)  a_discharge = a_discharge  / acctime_ucat
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

         CALL vector_gather_map2grid_and_write ( a_discharge, numucat,                     &
            totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
            file_hist_ucat, 'f_discharge', 'lon_ucat', 'lat_ucat', itime_in_file_ucat,     &
            'discharge in river and flood plain', 'm^3/s')

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

         CALL vector_gather_map2grid_and_write ( a_discharge, numucat,                            &
            totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat,        &
            file_hist_ucat, 'f_discharge_rivermouth', 'lon_ucat', 'lat_ucat', itime_in_file_ucat, &
            'river mouth discharge into ocean', 'm^3/s')
      ENDIF

      IF (DEF_hist_vars%floodfrc) THEN

         IF (p_is_worker) THEN

            IF (numucat > 0) a_floodarea = a_floodarea / acctime_ucat

            IF (numucat > 0) THEN
               allocate (a_floodfrc_ucat (numucat))
               WHERE (topo_area > 0)
                  a_floodfrc_ucat = a_floodarea / topo_area
               ELSE WHERE
                  a_floodfrc_ucat = spval
               END WHERE
            ENDIF

            IF (numinpm > 0) allocate (a_floodfrc_inpm (numinpm))

            CALL worker_push_data (push_ucat2inpm, a_floodfrc_ucat, a_floodfrc_inpm, &
               fillvalue = spval, mode = 'average')

            CALL worker_remap_data_grid2pset (remap_patch2inpm, a_floodfrc_inpm, a_floodfrc_pch, &
               fillvalue = spval, mode = 'average')
         ENDIF

      ENDIF


      IF (allocated (a_floodfrc_ucat)) deallocate (a_floodfrc_ucat)
      IF (allocated (a_floodfrc_inpm)) deallocate (a_floodfrc_inpm)
      IF (allocated (acc_vec_grid   )) deallocate (acc_vec_grid   )

#ifdef GridRiverLakeSediment
      CALL write_sediment_history (file_hist_ucat, itime_in_file_ucat)
#endif

      ! ----- reservoir variables -----
      IF (DEF_Reservoir_Method > 0) THEN
         IF (totalnumresv > 0) THEN

            IF (p_is_master) THEN
               IF (.not. fexists) THEN
                  CALL ncio_define_dimension(file_hist_ucat, 'reservoir', totalnumresv)
                  CALL ncio_write_serial (file_hist_ucat, 'resv_GRAND_ID' , dam_GRAND_ID, 'reservoir')
                  CALL ncio_put_attr (file_hist_ucat, 'resv_GRAND_ID', 'long_name', 'reservoir GRAND ID')
               ENDIF
            ENDIF

            IF (DEF_hist_vars%volresv) THEN
               IF (p_is_worker) THEN
                  IF (numresv > 0) THEN
                     WHERE (acctime_resv > 0)
                        a_volresv = a_volresv / acctime_resv
                     ELSEWHERE
                        a_volresv = spval
                     END WHERE
                  ENDIF
               ENDIF

               CALL vector_gather_and_write ( a_volresv, numresv, totalnumresv, resv_data_address, &
                  file_hist_ucat, 'volresv', 'reservoir', itime_in_file_ucat, 'reservoir water volume', 'm^3')
            ENDIF

            IF (DEF_hist_vars%qresv_in) THEN
               IF (p_is_worker) THEN
                  IF (numresv > 0) THEN
                     WHERE (acctime_resv > 0)
                        a_qresv_in = a_qresv_in / acctime_resv
                     ELSEWHERE
                        a_qresv_in = spval
                     END WHERE
                  ENDIF
               ENDIF

               CALL vector_gather_and_write ( a_qresv_in, numresv, totalnumresv, resv_data_address, &
                  file_hist_ucat, 'qresv_in', 'reservoir', itime_in_file_ucat, 'reservoir inflow', 'm^3/s')
            ENDIF

            IF (DEF_hist_vars%qresv_out) THEN
               IF (p_is_worker) THEN
                  IF (numresv > 0) THEN
                     WHERE (acctime_resv > 0)
                        a_qresv_out = a_qresv_out / acctime_resv
                     ELSEWHERE
                        a_qresv_out = spval
                     END WHERE
                  ENDIF
               ENDIF

               CALL vector_gather_and_write ( a_qresv_out, numresv, totalnumresv, resv_data_address, &
                  file_hist_ucat, 'qresv_out', 'reservoir', itime_in_file_ucat, 'reservoir outflow', 'm^3/s')
            ENDIF

         ENDIF
      ENDIF

      CALL flush_acc_fluxes_riverlake ()

   END SUBROUTINE hist_grid_riverlake_out

   !-----------------------
   SUBROUTINE flush_acc_fluxes_riverlake ()

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   USE MOD_Grid_Reservoir,        only: numresv
   IMPLICIT NONE

      IF (p_is_worker) THEN

         IF (numucat > 0) THEN
            acctime_ucat (:) = 0.
            a_wdsrf_ucat (:) = 0.
            a_veloc_riv  (:) = 0.
            a_discharge  (:) = 0.
            a_floodarea  (:) = 0.
         ENDIF

         IF (numresv > 0) THEN
            acctime_resv (:) = 0.
            a_volresv    (:) = 0.
            a_qresv_in   (:) = 0.
            a_qresv_out  (:) = 0.
         ENDIF

#ifdef GridRiverLakeSediment
         IF (DEF_USE_SEDIMENT .and. allocated(a_sedcon)) THEN
            a_sedcon   = 0.
            a_sedout   = 0.
            a_bedout   = 0.
            a_sedinp   = 0.
            a_netflw   = 0.
            a_layer    = 0.
            a_shearvel = 0.
            sed_hist_acctime = 0.
         ENDIF
#endif

      ENDIF

   END SUBROUTINE flush_acc_fluxes_riverlake

   !---------------------------------------
   SUBROUTINE hist_grid_riverlake_final ()

   IMPLICIT NONE

      IF (allocated(acctime_ucat    )) deallocate (acctime_ucat    )
      IF (allocated(a_wdsrf_ucat    )) deallocate (a_wdsrf_ucat    )
      IF (allocated(a_veloc_riv     )) deallocate (a_veloc_riv     )
      IF (allocated(a_discharge     )) deallocate (a_discharge     )
      IF (allocated(a_floodarea     )) deallocate (a_floodarea     )

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

#ifdef GridRiverLakeSediment
   !---------------------------------------
   SUBROUTINE write_sediment_history (file_hist_ucat, itime_in_file_ucat)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_data_address, &
                                         x_ucat, y_ucat, griducat
   USE MOD_Vector_ReadWrite
   IMPLICIT NONE

   character(len=*), intent(in) :: file_hist_ucat
   integer, intent(in) :: itime_in_file_ucat

   real(r8), allocatable :: a_sedcon_avg(:,:)
   real(r8), allocatable :: a_sedout_avg(:,:)
   real(r8), allocatable :: a_bedout_avg(:,:)
   real(r8), allocatable :: a_sedinp_avg(:,:)
   real(r8), allocatable :: a_netflw_avg(:,:)
   real(r8), allocatable :: a_layer_avg(:,:)
   real(r8), allocatable :: a_shearvel_avg(:)
   integer :: ised
   character(len=16) :: cised

      IF (.not. DEF_USE_SEDIMENT) RETURN

      ! Allocate on ALL processes (zero-size on non-workers) to avoid
      ! passing unallocated arrays to vector_gather_map2grid_and_write
      IF (p_is_worker .and. numucat > 0) THEN
         allocate (a_sedcon_avg  (nsed, numucat))
         allocate (a_sedout_avg  (nsed, numucat))
         allocate (a_bedout_avg  (nsed, numucat))
         allocate (a_sedinp_avg  (nsed, numucat))
         allocate (a_netflw_avg  (nsed, numucat))
         allocate (a_layer_avg   (nsed, numucat))
         allocate (a_shearvel_avg(numucat))

         IF (sed_hist_acctime > 0._r8) THEN
            a_shearvel_avg = a_shearvel / sed_hist_acctime
            DO ised = 1, nsed
               a_sedcon_avg(ised,:) = a_sedcon(ised,:) / sed_hist_acctime
               a_sedout_avg(ised,:) = a_sedout(ised,:) / sed_hist_acctime
               a_bedout_avg(ised,:) = a_bedout(ised,:) / sed_hist_acctime
               a_sedinp_avg(ised,:) = a_sedinp(ised,:) / sed_hist_acctime
               a_netflw_avg(ised,:) = a_netflw(ised,:) / sed_hist_acctime
               a_layer_avg(ised,:)  = a_layer(ised,:)  / sed_hist_acctime
            ENDDO
         ELSE
            a_shearvel_avg = 0.
            a_sedcon_avg   = 0.
            a_sedout_avg   = 0.
            a_bedout_avg   = 0.
            a_sedinp_avg   = 0.
            a_netflw_avg   = 0.
            a_layer_avg    = 0.
         ENDIF
      ELSE
         ! Allocate with nsed in first dim so a_xxx_avg(ised,:) is a valid zero-length slice
         allocate (a_sedcon_avg  (nsed, 0))
         allocate (a_sedout_avg  (nsed, 0))
         allocate (a_bedout_avg  (nsed, 0))
         allocate (a_sedinp_avg  (nsed, 0))
         allocate (a_netflw_avg  (nsed, 0))
         allocate (a_layer_avg   (nsed, 0))
         allocate (a_shearvel_avg(0))
      ENDIF

      IF (DEF_hist_vars%sedcon) THEN
         DO ised = 1, nsed
            WRITE(cised, '(I0)') ised
            CALL vector_gather_map2grid_and_write ( a_sedcon_avg(ised,:), numucat,     &
               totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
               file_hist_ucat, 'f_sedcon_' // trim(cised), 'lon_ucat', 'lat_ucat', itime_in_file_ucat, &
               'suspended sediment concentration, size class ' // trim(cised), 'm^3/m^3')
         ENDDO
      ENDIF

      IF (DEF_hist_vars%sedout) THEN
         DO ised = 1, nsed
            WRITE(cised, '(I0)') ised
            CALL vector_gather_map2grid_and_write ( a_sedout_avg(ised,:), numucat,     &
               totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
               file_hist_ucat, 'f_sedout_' // trim(cised), 'lon_ucat', 'lat_ucat', itime_in_file_ucat, &
               'suspended sediment flux, size class ' // trim(cised), 'm^3/s')
         ENDDO
      ENDIF

      IF (DEF_hist_vars%bedout) THEN
         DO ised = 1, nsed
            WRITE(cised, '(I0)') ised
            CALL vector_gather_map2grid_and_write ( a_bedout_avg(ised,:), numucat,     &
               totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
               file_hist_ucat, 'f_bedout_' // trim(cised), 'lon_ucat', 'lat_ucat', itime_in_file_ucat, &
               'bedload flux, size class ' // trim(cised), 'm^3/s')
         ENDDO
      ENDIF

      IF (DEF_hist_vars%sedinp) THEN
         DO ised = 1, nsed
            WRITE(cised, '(I0)') ised
            CALL vector_gather_map2grid_and_write ( a_sedinp_avg(ised,:), numucat,     &
               totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
               file_hist_ucat, 'f_sedinp_' // trim(cised), 'lon_ucat', 'lat_ucat', itime_in_file_ucat, &
               'sediment erosion input, size class ' // trim(cised), 'm^3/s')
         ENDDO
      ENDIF

      IF (DEF_hist_vars%netflw) THEN
         DO ised = 1, nsed
            WRITE(cised, '(I0)') ised
            CALL vector_gather_map2grid_and_write ( a_netflw_avg(ised,:), numucat,     &
               totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
               file_hist_ucat, 'f_netflw_' // trim(cised), 'lon_ucat', 'lat_ucat', itime_in_file_ucat, &
               'net bed-water exchange flux (incl. shallow deposit), size class ' // trim(cised), 'm^3/s')
         ENDDO
      ENDIF

      IF (DEF_hist_vars%sedlayer) THEN
         DO ised = 1, nsed
            WRITE(cised, '(I0)') ised
            CALL vector_gather_map2grid_and_write ( a_layer_avg(ised,:), numucat,      &
               totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
               file_hist_ucat, 'f_layer_' // trim(cised), 'lon_ucat', 'lat_ucat', itime_in_file_ucat, &
               'active layer storage, size class ' // trim(cised), 'm^3')
         ENDDO
      ENDIF

      IF (DEF_hist_vars%shearvel) THEN
         CALL vector_gather_map2grid_and_write ( a_shearvel_avg, numucat,              &
            totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
            file_hist_ucat, 'f_shearvel', 'lon_ucat', 'lat_ucat', itime_in_file_ucat,  &
            'shear velocity', 'm/s')
      ENDIF

      IF (allocated(a_sedcon_avg  )) deallocate (a_sedcon_avg  )
      IF (allocated(a_sedout_avg  )) deallocate (a_sedout_avg  )
      IF (allocated(a_bedout_avg  )) deallocate (a_bedout_avg  )
      IF (allocated(a_sedinp_avg  )) deallocate (a_sedinp_avg  )
      IF (allocated(a_netflw_avg  )) deallocate (a_netflw_avg  )
      IF (allocated(a_layer_avg   )) deallocate (a_layer_avg   )
      IF (allocated(a_shearvel_avg)) deallocate (a_shearvel_avg)

   END SUBROUTINE write_sediment_history
#endif

END MODULE MOD_Grid_RiverLakeHist
#endif
