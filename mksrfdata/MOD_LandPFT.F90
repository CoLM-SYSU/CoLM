#include <define.h>

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)

MODULE MOD_LandPFT

!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!    Build pixelset "landpft" (Plant Function Type).
!
!    In CoLM, the global/regional area is divided into a hierarchical structure:
!    1. If GRIDBASED or UNSTRUCTURED is defined, it is
!       ELEMENT >>> PATCH
!    2. If CATCHMENT is defined, it is
!       ELEMENT >>> HRU >>> PATCH
!    If Plant Function Type classification is used, PATCH is further divided into PFT.
!    If Plant Community classification is used,     PATCH is further divided into PC.
!
!    "landpft" refers to pixelset PFT.
!
!  Created by Shupeng Zhang, May 2023
!    porting codes from Hua Yuan's OpenMP version to MPI parallel version.
!-----------------------------------------------------------------------

   USE MOD_Namelist
   USE MOD_Grid
   USE MOD_Pixelset
   USE MOD_Const_LC
   USE MOD_Vars_Global
   IMPLICIT NONE

   ! ---- Instance ----
   integer :: numpft
   type(grid_type)     :: grid_pft
   type(pixelset_type) :: landpft

   integer , allocatable :: pft2patch   (:)  !patch index of a PFT
   integer , allocatable :: patch_pft_s (:)  !start PFT index of a patch
   integer , allocatable :: patch_pft_e (:)  !end PFT index of a patch

   ! ---- PUBLIC routines ----
   PUBLIC :: landpft_build

CONTAINS

   ! -------------------------------
   SUBROUTINE landpft_build (lc_year)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_Namelist
   USE MOD_5x5DataReadin
   USE MOD_LandPatch
   USE MOD_Land2mWMO
   USE MOD_AggregationRequestData
   USE MOD_Const_LC
#ifdef CROP
   USE MOD_LandCrop
#endif

   IMPLICIT NONE

   integer, intent(in) :: lc_year

   ! Local Variables
   character(len=256) :: dir, fname, cyear
   type (block_data_real8_3d) :: pctpft
   real(r8), allocatable :: pctpft_patch(:,:), pctpft_one(:,:)
   real(r8), allocatable :: area_one  (:)
   logical,  allocatable :: patchmask (:)
   integer  :: ipatch, ipft, npatch, npft, npft_glb
   integer  :: wmo_src, ipft_grass
   real(r8) :: sumarea, maxgrass

      IF (p_is_master) THEN
         write(*,'(A)') 'Making land plant function type tiles:'
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      landpft%has_shared = .true.

      IF (p_is_io) THEN

         CALL allocate_block_data (grid_pft, pctpft, N_PFT_modis, lb1 = 0)
         CALL flush_block_data (pctpft, 1.0)

         dir = trim(DEF_dir_rawdata) // trim(DEF_rawdata%pft%dir)

         ! add parameter input for time year
         write(cyear,'(i4.4)') lc_year

         fname = trim(DEF_rawdata%pft%fname) //'.'// trim(cyear)

         CALL read_5x5_data_pft (dir, fname, grid_pft, 'PCT_PFT', pctpft)

#ifdef USEMPI
         CALL aggregation_data_daemon (grid_pft, data_r8_3d_in1 = pctpft, n1_r8_3d_in1 = N_PFT_modis)
#endif
      ENDIF


      IF (p_is_worker) THEN

         IF (numpatch > 0) THEN
            allocate (pctpft_patch (0:N_PFT-1,numpatch))
            allocate (patchmask (numpatch))

            pctpft_patch(:,:) = 0
            patchmask (:) = .true.
         ENDIF

         DO ipatch = 1, numpatch

            IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN

               wmo_src  = wmo_source(landpatch%ielm(ipatch))
               maxgrass = maxval(pctpft_patch(12:14,wmo_src))

               IF (maxgrass > 0) THEN
                  ipft_grass = maxloc(pctpft_patch(12:14,wmo_src), dim=1) + 11
                  pctpft_patch(:,ipatch) = 0
                  pctpft_patch(ipft_grass,ipatch) = 1.
               ELSE
                  pctpft_patch(0,ipatch) = 1.
               ENDIF

               CYCLE
            ENDIF
#ifndef CROP
            IF (patchtypes(landpatch%settyp(ipatch)) == 0) THEN
#else
            IF (patchtypes(landpatch%settyp(ipatch)) == 0 .and. landpatch%settyp(ipatch)/=CROPLAND) THEN
#endif
               CALL aggregation_request_data (landpatch, ipatch, grid_pft, zip = .false., area = area_one, &
                  data_r8_3d_in1 = pctpft, data_r8_3d_out1 = pctpft_one, n1_r8_3d_in1 = N_PFT_modis, lb1_r8_3d_in1 = 0)

               sumarea = sum(area_one * sum(pctpft_one(0:N_PFT-1,:),dim=1))

               ! in case of no PFT data, set to 100% bare when patchtype=0,
               ! be consistent with Aggregation_PercentagesPFT.F90.
               IF (sumarea <= 0.0) THEN
                  pctpft_patch(0,ipatch) = 1.
               ELSE
                  DO ipft = 0, N_PFT-1
                     pctpft_patch(ipft,ipatch) = sum(pctpft_one(ipft,:) * area_one) / sumarea
                  ENDDO
               ENDIF

            ENDIF
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif

         IF (numpatch > 0) THEN
            npatch = count(patchmask)
            numpft = count(pctpft_patch > 0.)
#ifdef CROP
            numpft = numpft + count(landpatch%settyp == CROPLAND)
#endif
            IF (npatch > 0) THEN
               allocate (patch_pft_s (npatch))
               allocate (patch_pft_e (npatch))
            ENDIF
         ELSE
            numpft = 0
         ENDIF

         IF (numpft > 0) THEN

            allocate (pft2patch      (numpft))

            allocate (landpft%eindex (numpft))
            allocate (landpft%settyp (numpft))
            allocate (landpft%ipxstt (numpft))
            allocate (landpft%ipxend (numpft))
            allocate (landpft%ielm   (numpft))

            allocate (landpft%pctshared (numpft))
            landpft%pctshared(:) = 1.

            npft = 0
            npatch = 0
            DO ipatch = 1, numpatch
               IF (patchmask(ipatch)) THEN
                  npatch = npatch + 1

#ifndef CROP
                  IF (patchtypes(landpatch%settyp(ipatch)) == 0) THEN
#else
                  IF (patchtypes(landpatch%settyp(ipatch)) == 0 .and. landpatch%settyp(ipatch)/=CROPLAND) THEN
#endif
                     patch_pft_s(npatch) = npft + 1
                     patch_pft_e(npatch) = npft + count(pctpft_patch(:,ipatch) > 0)

                     DO ipft = 0, N_PFT-1
                        IF (pctpft_patch(ipft,ipatch) > 0) THEN
                           npft = npft + 1

                           landpft%ielm  (npft) = landpatch%ielm  (ipatch)
                           landpft%eindex(npft) = landpatch%eindex(ipatch)
                           landpft%ipxstt(npft) = landpatch%ipxstt(ipatch)
                           landpft%ipxend(npft) = landpatch%ipxend(ipatch)
                           landpft%settyp(npft) = ipft

                           landpft%pctshared(npft) = pctpft_patch(ipft,ipatch)

                           pft2patch(npft) = npatch
                        ENDIF
                     ENDDO
#ifdef CROP
                  ELSEIF (landpatch%settyp(ipatch) == CROPLAND) THEN
                     npft = npft + 1
                     patch_pft_s(npatch) = npft
                     patch_pft_e(npatch) = npft

                     landpft%ielm  (npft) = landpatch%ielm  (ipatch)
                     landpft%eindex(npft) = landpatch%eindex(ipatch)
                     landpft%ipxstt(npft) = landpatch%ipxstt(ipatch)
                     landpft%ipxend(npft) = landpatch%ipxend(ipatch)
                     landpft%settyp(npft) = cropclass(ipatch) + N_PFT - 1

                     landpft%pctshared(npft) = landpatch%pctshared(ipatch)

                     pft2patch(npft) = npatch
#endif
                  ELSE
                     patch_pft_s(npatch) = -1
                     patch_pft_e(npatch) = -1
                  ENDIF
               ENDIF
            ENDDO

         ENDIF

      ENDIF

      CALL landpatch%pset_pack(patchmask, numpatch)

      landpft%nset = numpft
      CALL landpft%set_vecgs

#ifdef USEMPI
      IF (p_is_worker) THEN
         CALL mpi_reduce (numpft, npft_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', npft_glb, ' plant function type tiles.'
         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      write(*,'(A,I12,A)') 'Total: ', numpft, ' plant function type tiles.'
#endif

      IF (allocated(pctpft_patch)) deallocate (pctpft_patch)
      IF (allocated(pctpft_one  )) deallocate (pctpft_one  )
      IF (allocated(area_one    )) deallocate (area_one    )
      IF (allocated(patchmask   )) deallocate (patchmask   )

   END SUBROUTINE landpft_build

   ! ----------------------
   SUBROUTINE map_patch_to_pft

   USE MOD_SPMD_Task
   USE MOD_LandPatch
   USE MOD_Const_LC
   IMPLICIT NONE

   integer :: ipatch, ipft

      IF (p_is_worker) THEN

         IF ((numpatch <= 0) .or. (numpft <= 0)) RETURN

         IF (allocated(patch_pft_s)) deallocate(patch_pft_s)
         IF (allocated(patch_pft_e)) deallocate(patch_pft_e)
         IF (allocated(pft2patch  )) deallocate(pft2patch  )

         allocate (patch_pft_s (numpatch))
         allocate (patch_pft_e (numpatch))
         allocate (pft2patch   (numpft  ))

         ipft = 1
         DO ipatch = 1, numpatch
#ifndef CROP
            IF (patchtypes(landpatch%settyp(ipatch)) == 0) THEN
#else
            IF (patchtypes(landpatch%settyp(ipatch)) == 0 .and. landpatch%settyp(ipatch)/=CROPLAND) THEN
#endif

               patch_pft_s(ipatch) = ipft

               DO WHILE (ipft <= numpft)
                  IF ((landpft%eindex(ipft) == landpatch%eindex(ipatch))  &
                     .and. (landpft%ipxstt(ipft) == landpatch%ipxstt(ipatch))  &
                     .and. (landpft%settyp(ipft) < N_PFT)) THEN
                     pft2patch  (ipft  ) = ipatch
                     patch_pft_e(ipatch) = ipft
                     ipft = ipft + 1
                  ELSE
                     EXIT
                  ENDIF
               ENDDO
#ifdef CROP
            ELSEIF (landpatch%settyp(ipatch) == CROPLAND) THEN
               patch_pft_s(ipatch) = ipft
               patch_pft_e(ipatch) = ipft
               pft2patch  (ipft  ) = ipatch
               ipft = ipft + 1
#endif
            ELSE
               patch_pft_s(ipatch) = -1
               patch_pft_e(ipatch) = -1
            ENDIF

         ENDDO

      ENDIF

   END SUBROUTINE map_patch_to_pft

END MODULE MOD_LandPFT
#endif
