#include <define.h>

#ifdef CROP
MODULE MOD_LandCrop

!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!    Build crop patches.
!
!  Created by Shupeng Zhang, Sep 2023
!    porting codes from Hua Yuan's OpenMP version to MPI parallel version.
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Grid
   IMPLICIT NONE

   ! ---- Instance ----
   type(grid_type) :: grid_crop
   integer,  allocatable :: cropclass (:)
   real(r8), allocatable :: cropfrac  (:)

CONTAINS

   ! -------------------------------
   SUBROUTINE landcrop_build (lc_year)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Block
   USE MOD_DataType
   USE MOD_LandElm
#ifdef CATCHMENT
   USE MOD_LandHRU
#endif
   USE MOD_LandPatch
   USE MOD_NetCDFBlock
   USE MOD_PixelsetShared
   USE MOD_5x5DataReadin

   IMPLICIT NONE

   integer, intent(in) :: lc_year

   ! Local Variables
   character(len=255) :: cyear, file_patch, dir, fname
   integer :: npatch_glb
   type(block_data_real8_2d) :: pctcrop_xy
   type(block_data_real8_3d) :: pctshared_xy
   type(block_data_real8_3d) :: cropdata
   integer :: sharedfilter(1), cropfilter(1)
   integer :: iblkme, ib, jb
   real(r8), allocatable :: pctshared  (:)
   integer , allocatable :: classshared(:)

      write(cyear,'(i4.4)') lc_year
      IF (p_is_master) THEN
         write(*,'(A)') 'Making patches (crop shared):'
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_io) THEN

         dir = trim(DEF_dir_rawdata)//trim(DEF_rawdata%pft%dir)
         fname = trim(DEF_rawdata%pft%fname)//'.'//trim(cyear)

         CALL allocate_block_data (grid_patch, pctcrop_xy)
         CALL read_5x5_data (dir, fname, grid_patch, 'PCT_CROP', pctcrop_xy)

         CALL allocate_block_data (grid_patch, pctshared_xy, 2)
         DO iblkme = 1, gblock%nblkme
            ib = gblock%xblkme(iblkme)
            jb = gblock%yblkme(iblkme)
            pctshared_xy%blk(ib,jb)%val(1,:,:) = 1. - pctcrop_xy%blk(ib,jb)%val/100.
            pctshared_xy%blk(ib,jb)%val(2,:,:) = pctcrop_xy%blk(ib,jb)%val/100.
         ENDDO
      ENDIF

      sharedfilter = (/ 1 /)

      IF (landpatch%has_shared) then
         CALL pixelsetshared_build (landpatch, grid_patch, pctshared_xy, 2, sharedfilter, &
            pctshared, classshared, fracin = landpatch%pctshared)
      ELSE
         CALL pixelsetshared_build (landpatch, grid_patch, pctshared_xy, 2, sharedfilter, &
            pctshared, classshared)
      ENDIF

      IF (p_is_worker) THEN
         IF (landpatch%nset > 0) THEN
            WHERE (classshared == 2) landpatch%settyp = CROPLAND
         ENDIF
      ENDIF

      IF (p_is_io) THEN
         file_patch = trim(DEF_dir_rawdata)//trim(DEF_rawdata%crop%dir)//&
                      '/global_CFT_surface_data.nc'
         CALL allocate_block_data (grid_crop, cropdata, N_CFT)
         CALL ncio_read_block (file_patch, 'PCT_CFT', grid_crop, N_CFT, cropdata)
      ENDIF

      cropfilter = (/ CROPLAND /)

      CALL pixelsetshared_build (landpatch, grid_crop, cropdata, N_CFT, cropfilter, &
         cropfrac, cropclass, fracin = pctshared)

      numpatch = landpatch%nset

      landpatch%has_shared = .true.
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            IF (allocated(landpatch%pctshared)) THEN
               deallocate(landpatch%pctshared)
            ENDIF

            allocate(landpatch%pctshared(numpatch))
            landpatch%pctshared = cropfrac
         ENDIF
      ENDIF

      IF (allocated(pctshared  )) deallocate(pctshared  )
      IF (allocated(classshared)) deallocate(classshared)

#ifdef USEMPI
      IF (p_is_worker) THEN
         CALL mpi_reduce (numpatch, npatch_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', npatch_glb, ' patches (with crop).'
         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      write(*,'(A,I12,A)') 'Total: ', numpatch, ' patches.'
#endif

IF ( .not. DEF_Output_2mWMO ) THEN
      CALL elm_patch%build (landelm, landpatch, use_frac = .true.)
#ifdef CATCHMENT
      CALL hru_patch%build (landhru, landpatch, use_frac = .true.)
#endif

      CALL write_patchfrac (DEF_dir_landdata, lc_year)
ENDIF

   END SUBROUTINE landcrop_build

END MODULE MOD_LandCrop
#endif
