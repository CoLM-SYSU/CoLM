#include <define.h>

SUBROUTINE Aggregation_ForestHeight ( &
      gland, dir_rawdata, dir_model_landdata, lc_year)

!-----------------------------------------------------------------------
!  Global Forest Height
!     (http://lidarradar.jpl.nasa.gov/)
!      Simard, M., N. Pinto, J. B. Fisher, and A. Baccini, 2011: Mapping
!      forest canopy height globally with spaceborne lidar.
!      J. Geophys. Res., 116, G04021.
!
!  Created by Yongjiu Dai, 02/2014
!
! !REVISIONS:
!  Hua Yuan,      ?/2020 : for land cover land use classifications
!  Shupeng Zhang, 01/2022: porting codes to MPI parallel version
!  Jiayi Xiang & Yuan, 12/2025: add crown structure data for tree PFTs.
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_Vars_Global, only: spval
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_LandPatch
   USE MOD_Land2mWMO
   USE MOD_Block, only: get_filename_block
   USE MOD_NetCDFVector
   USE MOD_NetCDFBlock
   USE MOD_NetCDFSerial, only: ncio_put_attr
#ifdef RangeCheck
   USE MOD_RangeCheck
#endif
   USE MOD_AggregationRequestData
   USE MOD_Utils

   USE MOD_Const_LC
   USE MOD_Const_PFT, only: htop0_p, hbot0_p
   USE MOD_5x5DataReadin
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_LandPFT
#endif

#ifdef SrfdataDiag
   USE MOD_SrfdataDiag
#endif

   IMPLICIT NONE

   ! arguments:
   integer, intent(in) :: lc_year
   type(grid_type),  intent(in) :: gland
   character(len=*), intent(in) :: dir_rawdata
   character(len=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ---------------------------------------------------------------
   character(len=256) :: landdir, lndname, cyear
   integer :: L, ipatch, p
   integer :: wmo_src

   type (block_data_real8_2d) :: tree_height
   real(r8), allocatable :: tree_height_patches(:), tree_height_one(:)

   ! for IGBP data
   character(len=256) :: dir, fname
   type (block_data_real8_2d) :: htop
   type (block_data_real8_3d) :: pftPCT,cdepth
#ifdef LULC_IGBP_PC
   type (block_data_real8_3d) :: cratio
#endif
   real(r8), allocatable :: htop_patches(:), htop_pfts(:), htop_pcs(:,:)
   real(r8), allocatable :: hbot_pfts(:)
#ifdef LULC_IGBP_PC
   real(r8), allocatable :: cratio_pfts(:)
#endif
   real(r8), allocatable :: htop_one(:), area_one(:), pct_one(:,:)
   real(r8), allocatable :: cdepth_one(:,:)
#ifdef LULC_IGBP_PC
   real(r8), allocatable :: cratio_one(:,:)
#endif
   integer  :: ip, ipft, iblkgrp, iblk, jblk
   real(r8) :: sumarea
   character(len=256) :: fileblock

#ifdef SrfdataDiag
   integer :: typpatch(N_land_classification+1), ityp
#ifndef CROP
   integer :: typpft  (N_PFT)
#else
   integer :: typpft  (N_PFT+N_CFT)
#endif
   integer :: typpc   (N_land_classification+1)
#endif

      write(cyear,'(i4.4)') lc_year
      landdir = trim(dir_model_landdata) // '/htop/' //trim(cyear)
      ! print*,'landdir: ', trim(landdir)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,'(/, A)') 'Aggregate forest height ...'
         CALL system('mkdir -p ' // trim(adjustl(landdir)))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef LULC_USGS
      !TODO: add dir, file name and var name
      ! to check if they can be combined. -done
      lndname = trim(DEF_dir_rawdata) // trim(DEF_rawdata%htop%dir) &
              // '/Forest_Height.nc'

      IF (p_is_io) THEN
         CALL allocate_block_data (gland, tree_height)
         CALL ncio_read_block (lndname, 'forest_height', gland, tree_height)

#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = tree_height)
#endif
      ENDIF

      IF (p_is_worker) THEN

         allocate (tree_height_patches (numpatch))

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)

            IF(L/=0 .and. L/=1 .and. L/=16 .and. L/=24)THEN
               ! NOT OCEAN(0)/URBAN and BUILT-UP(1)/WATER BODIES(16)/ICE(24)
               CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, &
                  data_r8_2d_in1 = tree_height, data_r8_2d_out1 = tree_height_one)
               tree_height_patches (ipatch) = median (tree_height_one, size(tree_height_one))
            ELSE
               tree_height_patches (ipatch) = -1.0e36_r8
            ENDIF
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('htop_patches ', tree_height_patches)
#endif

      lndname = trim(landdir)//'/htop_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'htop_patches', 'patch', landpatch, tree_height_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typpatch = (/(ityp, ityp = 0, N_land_classification)/)
      lndname  = trim(dir_model_landdata) // '/diag/htop_patch_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (tree_height_patches, landpatch%settyp, typpatch, m_patch2diag, &
         -1.0e36_r8, lndname, 'htop', compress = 6, write_mode = 'one', defval=0._r8, create_mode=.true.)
#endif

      IF (p_is_worker) THEN
         deallocate ( tree_height_patches )
      ENDIF
#endif


#ifdef LULC_IGBP
      IF (p_is_io) THEN
         CALL allocate_block_data (gland, htop)
      ENDIF

      IF (p_is_io) THEN

         IF (DEF_rawdata_namelist == "colm2024.nml") THEN
            dir   = trim(dir_rawdata) // trim(DEF_rawdata%pft%dir)
            fname = trim(DEF_rawdata%pft%fname) //'.'// trim(cyear)

            CALL read_5x5_data (dir, fname, gland, 'HTOP', htop)
         ELSE
            dir   = trim(dir_rawdata) // trim(DEF_rawdata%htop%dir)
            fname = trim(DEF_rawdata%htop%fname)

            CALL read_5x5_data (dir, fname, gland, trim(DEF_rawdata%htop%vname), htop)
         ENDIF

#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = htop)
#endif
      ENDIF

      IF (p_is_worker) THEN

         allocate (htop_patches (numpatch))

         DO ipatch = 1, numpatch

            IF (landpatch%settyp(ipatch) /= 0) THEN

               CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, &
                  area = area_one, data_r8_2d_in1 = htop, data_r8_2d_out1 = htop_one)

               where (htop_one < 0.) htop_one = 0.

               htop_patches(ipatch) = sum(htop_one * area_one) / sum(area_one)
            ENDIF

         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('HTOP_patches ', htop_patches)
#endif

      lndname = trim(landdir)//'/htop_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'htop_patches', 'patch', landpatch, htop_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typpatch = (/(ityp, ityp = 0, N_land_classification)/)
      lndname  = trim(dir_model_landdata) // '/diag/htop_patch_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (htop_patches, landpatch%settyp, typpatch, m_patch2diag, &
         -1.0e36_r8, lndname, 'htop', compress = 6, write_mode = 'one', defval=0._r8, create_mode=.true.)
#endif

      IF (p_is_worker) THEN
         IF (allocated(htop_patches)) deallocate (htop_patches)
         IF (allocated(htop_one    )) deallocate (htop_one    )
         IF (allocated(area_one    )) deallocate (area_one    )
      ENDIF
#endif

! ================PFT/PC htop
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      IF (p_is_io) THEN
         CALL allocate_block_data (grid_pft, pftPCT, N_PFT_modis, lb1 = 0)
         IF (DEF_RS_CROWN_STRUCTURE) THEN
            CALL allocate_block_data (grid_pft, cdepth, N_PFT_modis, lb1 = 0)
#ifdef LULC_IGBP_PC
            CALL allocate_block_data (grid_pft, cratio, N_PFT_modis, lb1 = 0)
#endif
         ENDIF
         CALL allocate_block_data (gland   , htop)
      ENDIF

      IF (p_is_io) THEN

         IF (DEF_rawdata_namelist == "colm2024.nml") THEN
            dir   = trim(DEF_dir_rawdata) // trim(DEF_rawdata%pft%dir)
            fname = trim(DEF_rawdata%pft%fname) //'.'// trim(cyear)
            CALL read_5x5_data (dir, fname, gland, 'HTOP', htop)
         ELSE
            dir   = trim(DEF_dir_rawdata) // trim(DEF_rawdata%htop%dir)
            fname = trim(DEF_rawdata%htop%fname)
            CALL read_5x5_data (dir, fname, gland, trim(DEF_rawdata%htop%vname), htop)
         ENDIF

         dir   = trim(DEF_dir_rawdata) // trim(DEF_rawdata%pft%dir)
         fname = trim(DEF_rawdata%pft%fname) //'.'// trim(cyear)
         CALL read_5x5_data_pft (dir, fname, grid_pft, 'PCT_PFT', pftPCT)

         IF (DEF_RS_CROWN_STRUCTURE) THEN
            ! read crown depth from crown structure dataset, hbot = htop - cdepth
            dir = trim(DEF_dir_rawdata) // trim(DEF_rawdata%cdepth%dir)
            fname = trim(DEF_rawdata%cdepth%fname)
            CALL read_5x5_data_pft (dir, fname, grid_pft, trim(DEF_rawdata%cdepth%vname), &
               cdepth, lb=1, ub=8)
#ifdef LULC_IGBP_PC
            ! read crown aspect ratio from crown structure dataset
            dir = trim(DEF_dir_rawdata) // trim(DEF_rawdata%cratio%dir)
            fname = trim(DEF_rawdata%cratio%fname)
            CALL read_5x5_data_pft (dir, fname, grid_pft, trim(DEF_rawdata%cratio%vname), &
               cratio, lb=1, ub=8)
#endif
         ENDIF

#ifdef USEMPI
         IF (DEF_RS_CROWN_STRUCTURE) THEN
#ifdef LULC_IGBP_PC
            CALL aggregation_data_daemon_multigrid (grid_in1=grid_pft, data_r8_3d_in1=pftPCT, n1_r8_3d_in1=16, &
               grid_in2=gland, data_r8_2d_in2=htop, grid_in3=grid_pft, data_r8_3d_in3=cdepth, n1_r8_3d_in3=16, &
               grid_in4=grid_pft, data_r8_3d_in4=cratio, n1_r8_3d_in4=16)
#else
            CALL aggregation_data_daemon_multigrid (grid_in1=grid_pft, data_r8_3d_in1=pftPCT, n1_r8_3d_in1=16, &
               grid_in2=gland, data_r8_2d_in2=htop, grid_in3=grid_pft, data_r8_3d_in3=cdepth, n1_r8_3d_in3=16)
#endif
         ELSE
            CALL aggregation_data_daemon_multigrid (grid_in1=grid_pft, data_r8_3d_in1=pftPCT, n1_r8_3d_in1=16, &
               grid_in2=gland, data_r8_2d_in2=htop)
         ENDIF
#endif
      ENDIF

      IF (p_is_worker) THEN

         allocate (htop_patches   (numpatch))
         allocate (htop_pfts      (numpft  ))
         IF (DEF_RS_CROWN_STRUCTURE) THEN
            allocate (hbot_pfts   (numpft))
            hbot_pfts = spval
#ifdef LULC_IGBP_PC
            allocate (cratio_pfts (numpft))
            cratio_pfts = spval
#endif
         ENDIF

         DO ipatch = 1, numpatch

            IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
               wmo_src = wmo_source (landpatch%ielm(ipatch))

               ! set patch htop
               htop_patches(ipatch) = htop_patches(wmo_src)
               ! set pft htop at the same time
               ip = patch_pft_s(ipatch)
               htop_pfts(ip) = htop_patches(ipatch)

               CYCLE
            ENDIF
            IF (DEF_RS_CROWN_STRUCTURE) THEN
#ifdef LULC_IGBP_PC
               CALL aggregation_request_data_multigrid(landpatch, ipatch, grid_in1=grid_pft, area=area_one, &
                  data_r8_3d_in1=pftPCT, data_r8_3d_out1=pct_one, n1_r8_3d_in1=16, lb1_r8_3d_in1=0, &
                  grid_in2=gland, data_r8_2d_in2=htop, data_r8_2d_out2=htop_one, &
                  grid_in3=grid_pft, data_r8_3d_in3=cdepth, data_r8_3d_out3=cdepth_one, n1_r8_3d_in3=16, lb1_r8_3d_in3=0, &
                  grid_in4=grid_pft, data_r8_3d_in4=cratio, data_r8_3d_out4=cratio_one, n1_r8_3d_in4=16, lb1_r8_3d_in4=0)
#else
               CALL aggregation_request_data_multigrid(landpatch, ipatch, grid_in1=grid_pft, area=area_one, &
                  data_r8_3d_in1=pftPCT, data_r8_3d_out1=pct_one, n1_r8_3d_in1=16, lb1_r8_3d_in1=0, &
                  grid_in2=gland, data_r8_2d_in2=htop, data_r8_2d_out2=htop_one, &
                  grid_in3=grid_pft, data_r8_3d_in3=cdepth, data_r8_3d_out3=cdepth_one, n1_r8_3d_in3=16, lb1_r8_3d_in3=0)
#endif
            ELSE
               CALL aggregation_request_data_multigrid(landpatch, ipatch, grid_in1=grid_pft, area=area_one, &
                  data_r8_3d_in1=pftPCT, data_r8_3d_out1=pct_one, n1_r8_3d_in1=16, lb1_r8_3d_in1=0, &
                  grid_in2=gland, data_r8_2d_in2=htop, data_r8_2d_out2=htop_one)
            ENDIF
            where (htop_one < 0.) htop_one = 0.

            htop_patches(ipatch) = sum(htop_one * area_one) / sum(area_one)
            pct_one = max(pct_one , 0.0)

#ifndef CROP
            IF (patchtypes(landpatch%settyp(ipatch)) == 0) THEN
#else
            IF (patchtypes(landpatch%settyp(ipatch)) == 0 .and. landpatch%settyp(ipatch)/=CROPLAND) THEN
#endif
               DO ip = patch_pft_s(ipatch), patch_pft_e(ipatch)
                  p = landpft%settyp(ip)
                  sumarea = sum(pct_one(p,:) * area_one)
                  IF (sumarea > 0) THEN
                     htop_pfts(ip) = sum(htop_one * pct_one(p,:) * area_one) / sumarea
                  ELSE
                     htop_pfts(ip) = htop_patches(ipatch)
                  ENDIF
                  IF (DEF_RS_CROWN_STRUCTURE) THEN
                     sumarea = sum(area_one, &
                        mask=cdepth_one(p,:) > 0._r8 .and. cdepth_one(p,:) < htop_one)
                     IF (sumarea > 0._r8) THEN
                        hbot_pfts(ip) = sum((htop_one-cdepth_one(p,:)) * area_one, &
                           mask=cdepth_one(p,:) > 0._r8 .and. cdepth_one(p,:) < htop_one) / sumarea
                     ELSEIF (any(cdepth_one(p,:) > 0._r8 .and. htop_one > 0._r8 .and. &
                        cdepth_one(p,:) >= htop_one)) THEN
                        hbot_pfts(ip) = htop_pfts(ip) * hbot0_p(p) / htop0_p(p)
                        hbot_pfts(ip) = max(1._r8, hbot_pfts(ip))
                     ENDIF
#ifdef LULC_IGBP_PC
                     sumarea = sum(area_one, mask=cratio_one(p,:) > 0._r8)
                     IF (sumarea > 0._r8) THEN
                        cratio_pfts(ip) = sum(cratio_one(p,:) * area_one, &
                           mask=cratio_one(p,:) > 0._r8) / sumarea
                     ENDIF
#endif
                  ENDIF
               ENDDO
#ifdef CROP
            ELSEIF (landpatch%settyp(ipatch) == CROPLAND) THEN
               ip = patch_pft_s(ipatch)
               htop_pfts(ip) = htop_patches(ipatch)
#endif
            ENDIF
         ENDDO

#ifdef USEMPI
      CALL aggregation_worker_done_multigrid ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('HTOP_patches ', htop_patches)
      CALL check_vector_data ('HTOP_pfts    ', htop_pfts   )
      IF (DEF_RS_CROWN_STRUCTURE) THEN
         CALL check_vector_data ('HBOT_pfts    ', hbot_pfts)
#ifdef LULC_IGBP_PC
         CALL check_vector_data ('CRATIO_pfts  ', cratio_pfts)
#endif
      ENDIF
#endif

      lndname = trim(landdir)//'/htop_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'htop_patches', 'patch', landpatch, htop_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typpatch = (/(ityp, ityp = 0, N_land_classification)/)
      lndname  = trim(dir_model_landdata) // '/diag/htop_patch_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (htop_patches, landpatch%settyp, typpatch, m_patch2diag, &
         -1.0e36_r8, lndname, 'htop', compress = 6, write_mode = 'one', defval=0._r8, create_mode=.true.)
#endif

      lndname = trim(landdir)//'/htop_pfts.nc'
      CALL ncio_create_file_vector (lndname, landpft)
      CALL ncio_define_dimension_vector (lndname, landpft, 'pft')
      CALL ncio_write_vector (lndname, 'htop_pfts', 'pft', landpft, htop_pfts, DEF_Srfdata_CompressLevel)

      IF (DEF_RS_CROWN_STRUCTURE) THEN
         lndname = trim(landdir)//'/hbot_pfts.nc'
         CALL ncio_create_file_vector (lndname, landpft)
         CALL ncio_define_dimension_vector (lndname, landpft, 'pft')
         CALL ncio_write_vector (lndname, 'hbot_pfts', 'pft', landpft, hbot_pfts, DEF_Srfdata_CompressLevel)
         IF (p_is_io) THEN
            DO iblkgrp = 1, landpft%nblkgrp
               iblk = landpft%xblkgrp(iblkgrp)
               jblk = landpft%yblkgrp(iblkgrp)
               CALL get_filename_block (lndname, iblk, jblk, fileblock)
               CALL ncio_put_attr (fileblock, 'hbot_pfts', 'source', 'Derived from canopy height and crown depth')
               CALL ncio_put_attr (fileblock, 'hbot_pfts', 'references', &
                  'Crown depth: https://doi.org/10.1016/j.jag.2026.105401')
            ENDDO
         ENDIF
#ifdef LULC_IGBP_PC
         lndname = trim(landdir)//'/cratio_pfts.nc'
         CALL ncio_create_file_vector (lndname, landpft)
         CALL ncio_define_dimension_vector (lndname, landpft, 'pft')
         CALL ncio_write_vector (lndname, 'cratio_pfts', 'pft', landpft, cratio_pfts, DEF_Srfdata_CompressLevel)
         IF (p_is_io) THEN
            DO iblkgrp = 1, landpft%nblkgrp
               iblk = landpft%xblkgrp(iblkgrp)
               jblk = landpft%yblkgrp(iblkgrp)
               CALL get_filename_block (lndname, iblk, jblk, fileblock)
               CALL ncio_put_attr (fileblock, 'cratio_pfts', 'source', 'Remote-sensing crown structure data')
               CALL ncio_put_attr (fileblock, 'cratio_pfts', 'references', &
                  'https://doi.org/10.1016/j.jag.2026.105401')
            ENDDO
         ENDIF
#endif
      ENDIF
#ifdef SrfdataDiag
#ifndef CROP
      typpft  = (/(ityp, ityp = 0, N_PFT-1)/)
#else
      typpft  = (/(ityp, ityp = 0, N_PFT+N_CFT-1)/)
#endif
      lndname = trim(dir_model_landdata) // '/diag/htop_pft_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (htop_pfts, landpft%settyp, typpft, m_pft2diag, &
         -1.0e36_r8, lndname, 'htop_pft', compress = 6, write_mode = 'one', defval=0._r8, create_mode=.true.)
      IF (DEF_RS_CROWN_STRUCTURE) THEN
         lndname = trim(dir_model_landdata) // '/diag/hbot_pft_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (hbot_pfts, landpft%settyp, typpft, m_pft2diag, &
            spval, lndname, 'hbot_pft', compress=6, write_mode='one', defval=spval, create_mode=.true.)
         IF (p_is_master) THEN
            CALL ncio_put_attr (lndname, 'hbot_pft', 'source', 'Derived from canopy height and crown depth')
            CALL ncio_put_attr (lndname, 'hbot_pft', 'references', &
               'Crown depth: https://doi.org/10.1016/j.jag.2026.105401')
         ENDIF
#ifdef LULC_IGBP_PC
         lndname = trim(dir_model_landdata) // '/diag/cratio_pft_' // trim(cyear) // '.nc'
         CALL srfdata_map_and_write (cratio_pfts, landpft%settyp, typpft, m_pft2diag, &
            spval, lndname, 'cratio_pft', compress=6, write_mode='one', defval=spval, create_mode=.true.)
         IF (p_is_master) THEN
            CALL ncio_put_attr (lndname, 'cratio_pft', 'source', 'Remote-sensing crown structure data')
            CALL ncio_put_attr (lndname, 'cratio_pft', 'references', &
               'https://doi.org/10.1016/j.jag.2026.105401')
         ENDIF
#endif
      ENDIF
#endif

      IF (p_is_worker) THEN
         IF (allocated(htop_patches)) deallocate (htop_patches)
         IF (allocated(htop_pfts   )) deallocate (htop_pfts   )
         IF (allocated(htop_one    )) deallocate (htop_one    )
         IF (allocated(pct_one     )) deallocate (pct_one     )
         IF (allocated(area_one    )) deallocate (area_one    )
         IF (allocated(hbot_pfts   )) deallocate (hbot_pfts   )
         IF (allocated(cdepth_one  )) deallocate (cdepth_one  )
#ifdef LULC_IGBP_PC
         IF (allocated(cratio_pfts )) deallocate (cratio_pfts )
         IF (allocated(cratio_one  )) deallocate (cratio_one  )
#endif
      ENDIF
#endif

END SUBROUTINE Aggregation_ForestHeight
