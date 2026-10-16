#include <define.h>

SUBROUTINE Aggregation_PercentagesPFT (gland, dir_rawdata, dir_model_landdata, lc_year)

!-----------------------------------------------------------------------
!  Percentage of Plant Function Types
!
!  Original from Hua Yuan's OpenMP version.
!
! !REVISIONS:
!  Hua Yuan,      ?/2020 : for land cover land use classifications
!  Shupeng Zhang, 01/2022: porting codes to MPI parallel version
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_LandPatch
   USE MOD_Land2mWMO
#ifdef CROP
   USE MOD_LandCrop
#endif
   USE MOD_NetCDFBlock
   USE MOD_NetCDFVector
#ifdef RangeCheck
   USE MOD_RangeCheck
#endif
   USE MOD_AggregationRequestData

   USE MOD_Const_LC
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
   ! ----------------------------------------------------------------------
   character(len=256) :: landdir, lndname

   ! for IGBP data
   character(len=256) :: dir, fname, cyear
   ! for PFT
   type (block_data_real8_3d) :: pftPCT
   real(r8), allocatable :: pct_one(:), area_one(:)
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   real(r8), allocatable :: pct_pft_one(:,:)
   real(r8), allocatable :: pct_pfts(:)
#endif
   integer  :: ipatch, ipc, ipft, p
   integer  :: wmo_pth
   real(r8) :: sumarea, sum_pct_pfts
#ifdef SrfdataDiag
   integer :: typpft(0:N_PFT-1)
#ifdef CROP
   integer :: typcrop(N_CFT), ityp
#endif
#endif

      write(cyear,'(i4.4)') lc_year
      landdir = trim(dir_model_landdata) // '/pctpft/' // trim(cyear)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,'(/, A43)') 'Aggregate plant function type fractions ...'
         CALL system('mkdir -p ' // trim(adjustl(landdir)))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif


#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)

      dir   = trim(dir_rawdata) // trim(DEF_rawdata%pft%dir)
      fname = trim(DEF_rawdata%pft%fname)//'.'//trim(cyear)

      IF (p_is_io) THEN
         CALL allocate_block_data (gland, pftPCT, N_PFT_modis, lb1 = 0)
         CALL read_5x5_data_pft   (dir, fname, gland, 'PCT_PFT', pftPCT)
#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_3d_in1 = pftPCT, n1_r8_3d_in1 = N_PFT_modis)
#endif
      ENDIF

      IF (p_is_worker) THEN

         allocate(pct_pfts (numpft))
         pct_pfts(:) = 0.

         DO ipatch = 1, numpatch

            wmo_pth = wmo_patch(landpatch%ielm(ipatch))
            IF (ipatch == wmo_pth) THEN
               ipft = patch_pft_s(ipatch)

               pct_pfts(ipft) = 1.

               CYCLE
            ENDIF

#ifndef CROP
            IF (patchtypes(landpatch%settyp(ipatch)) == 0) THEN
#else
            IF (patchtypes(landpatch%settyp(ipatch)) == 0 &
               .and. landpatch%settyp(ipatch)/=CROPLAND) THEN
#endif
               CALL aggregation_request_data (landpatch, ipatch, gland, &
                  zip = USE_zip_for_aggregation, area = area_one, &
                  data_r8_3d_in1 = pftPCT, data_r8_3d_out1 = pct_pft_one, &
                  n1_r8_3d_in1 = N_PFT_modis, lb1_r8_3d_in1 = 0)

#ifdef CROP
               pct_pft_one(N_PFT_modis-1,:) = 0.
#endif

               pct_pft_one = max(pct_pft_one , 0.0)
               pct_one = sum(pct_pft_one, dim=1)
               pct_one = max(pct_one, 1.0e-6)
               sumarea = sum(area_one)

               DO ipft = patch_pft_s(ipatch), patch_pft_e(ipatch)
                  p = landpft%settyp(ipft)
                  pct_pfts(ipft) = sum(pct_pft_one(p,:) / pct_one * area_one) / sumarea
               ENDDO

               sum_pct_pfts = sum(pct_pfts(patch_pft_s(ipatch):patch_pft_e(ipatch)))

               IF (sum_pct_pfts > 0) THEN
                  pct_pfts(patch_pft_s(ipatch):patch_pft_e(ipatch)) = &
                     pct_pfts(patch_pft_s(ipatch):patch_pft_e(ipatch)) / sum_pct_pfts
               ELSE
                  ! in case of no PFT exist, but there is a patch type:
                  ! set bare soil 100%, be consistent with MOD_LandPFT.F90
                  pct_pfts(patch_pft_s(ipatch)) = 1.
               ENDIF

#ifdef CROP
            ELSEIF (landpatch%settyp(ipatch) == CROPLAND) THEN
               pct_pfts(patch_pft_s(ipatch):patch_pft_e(ipatch)) = 1.
#endif
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
      CALL check_vector_data ('PCT_PFTs ', pct_pfts)
#endif

      lndname = trim(landdir)//'/pct_pfts.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpft, 'pft')
      CALL ncio_write_vector (lndname, 'pct_pfts', 'pft', &
         landpft, pct_pfts, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typpft = (/(ipft, ipft = 0, N_PFT-1)/)
      lndname = trim(dir_model_landdata)//'/diag/pftfrac_elm_'//trim(cyear)//'.nc'
      CALL srfdata_map_and_write (pct_pfts, landpft%settyp, typpft, m_pft2diag, &
         -1.0e36_r8, lndname, 'pftfrac_elm', compress = 6, write_mode = 'one',  &
         defval=0._r8, stat_mode = 'fraction', pctshared = landpft%pctshared,   &
         create_mode=.true.)
#endif

      IF (p_is_worker) THEN
         IF (allocated(pct_pfts   )) deallocate(pct_pfts   )
         IF (allocated(pct_one    )) deallocate(pct_one    )
         IF (allocated(area_one   )) deallocate(area_one   )
         IF (allocated(pct_pft_one)) deallocate(pct_pft_one)
      ENDIF

#if (defined CROP)
      lndname = trim(landdir)//'/pct_crops.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'pct_crops', 'patch', &
         landpatch, cropfrac, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typcrop = (/(ityp, ityp = 1, N_CFT)/)
      lndname = trim(dir_model_landdata) // '/diag/cropfrac_elm_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (cropfrac, cropclass, typcrop, m_patch2diag,   &
         -1.0e36_r8, lndname, 'cropfrac_elm', compress = 6, write_mode = 'one', &
         defval=0._r8, stat_mode = 'fraction', pctshared = landpatch%pctshared, &
         create_mode=.true.)
#endif
#endif

#endif

END SUBROUTINE Aggregation_PercentagesPFT
