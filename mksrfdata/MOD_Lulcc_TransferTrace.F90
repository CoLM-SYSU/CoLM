#include <define.h>

MODULE MOD_Lulcc_TransferTrace

!------------------------------------------------------------------------
!
! !DESCRIPTION:
!  The transfer matrix and patch tracing vector were created using the
!  land cover type data of the adjacent two years. Based on next year's
!  patch, the pixels within the patch and last years' land cover type of
!  these pixels were obtained. Then the percent of source land cover
!  type of each patch was derived.
!
!  Created by Wanyi Lin, Shupeng Zhang and Hua Yuan, 07/2023
!
! !HISTORY:
!  05/2025, Wanyi Lin and Hua Yuan: code moved from main/LULCC/, now generate
!           the transfer matrix and patch tracing vector when making surface
!           data and SAVE it to a type (lulcc) of surface data.
!------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Vars_Global
   IMPLICIT NONE
   SAVE

   real(r8), allocatable, dimension(:,:) :: lccpct_patches(:,:) !frac of source patches in a patch
   real(r8), allocatable, dimension(:,:) :: lccpct_matrix (:,:) !frac of source patches in a grid

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_LulccTransferTrace
   PUBLIC :: deallocate_LulccTransferTrace
   PUBLIC :: MAKE_LulccTransferTrace

   ! PRIVATE MEMBER FUNCTIONS:


CONTAINS


   SUBROUTINE allocate_LulccTransferTrace
   ! --------------------------------------------------------------------
   ! Allocates memory for Lulcc time invariant variables
   ! --------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_LandPatch
   USE MOD_SPMD_Task

   IMPLICIT NONE

   integer :: nlc = N_land_classification

      IF (p_is_worker) THEN
         allocate (lccpct_patches (numpatch, 0:nlc))
         lccpct_patches (:,:) = 0
         allocate (lccpct_matrix  (numpatch, 0:nlc))
         lccpct_matrix  (:,:) = 0
      ENDIF

   END SUBROUTINE allocate_LulccTransferTrace


   SUBROUTINE MAKE_LulccTransferTrace (lc_year)

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_LandPatch
   USE MOD_Land2mWMO
   USE MOD_NetCDFVector
   USE MOD_NetCDFBlock
   USE MOD_AggregationRequestData
   USE MOD_Mesh
   USE MOD_MeshFilter
   USE MOD_LandElm
   USE MOD_DataType
   USE MOD_Block
   USE MOD_Pixel
   USE MOD_5x5DataReadin
   USE MOD_RegionClip
   USE MOD_Utils
#ifdef SrfdataDiag
   USE MOD_SrfdataDiag
#endif
#ifdef RangeCheck
   USE MOD_RangeCheck
#endif

   IMPLICIT NONE

   integer, intent(in) :: lc_year

!-------------------------- Local Variables ----------------------------
   character(len=256) :: dir_5x5, fname, lastyr, thisyr, dir_landdata, lndname
   character(len=4)   :: c2
   integer :: i,ipatch,ipxl,ipxstt,ipxend,numpxl,ilc
   integer, allocatable, dimension(:) :: locpxl
   type (block_data_int32_2d)         :: lcdatafr !land cover data of last year
   integer, allocatable, dimension(:) :: lcdatafr_one(:), lcfrbuff(:)
   real(r8),allocatable, dimension(:) :: area_one(:)    , areabuff(:)
   real(r8) :: sum_areabuff, gridarea
   integer, allocatable, dimension(:) :: grid_patch_s, grid_patch_e
   logical :: first_call
! for surface data diag
#ifdef SrfdataDiag
   integer  :: ityp
   integer, allocatable, dimension(:) :: typindex
#endif
!-----------------------------------------------------------------------

      IF ( (lc_year < 1990) .or. (lc_year < 2000 .and. MOD(lc_year, 5) /= 0) ) RETURN

      write(thisyr,'(i4.4)') lc_year
      IF (lc_year <= 2000 .and. MOD(lc_year, 5) == 0) THEN
         write(lastyr,'(i4.4)') MAX(1985, lc_year-5)
      ELSE
         write(lastyr,'(i4.4)') MAX(1985, lc_year-1)
      ENDIF

      first_call = .true.

#ifdef SrfdataDiag
      allocate( typindex(N_land_classification+1) )
#endif

      CALL allocate_LulccTransferTrace

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif


      IF (p_is_io) THEN

         CALL allocate_block_data (grid_patch, lcdatafr)
         dir_5x5 = trim(DEF_dir_rawdata) // trim(DEF_rawdata%landcover%dir)
         fname = trim(DEF_rawdata%landcover%fname)

         ! read the previous year land cover data
         CALL read_5x5_data (dir_5x5, fname, grid_patch, 'LC', lcdatafr)

#ifdef USEMPI
         CALL aggregation_data_daemon (grid_patch, data_i4_2d_in1 = lcdatafr)
#endif
      ENDIF

      ! -----------------------------------------------------------------
      ! extract the land cover type of pixels of last year for each patch
      ! -----------------------------------------------------------------
      IF (p_is_worker) THEN

         ! allocate with numelm
         allocate(grid_patch_s (numelm ))
         allocate(grid_patch_e (numelm ))

         grid_patch_e (:) = -1
         grid_patch_s (:) = -1

         DO i=1, numelm
            ! how many patches in ith element in this worker
            numpxl = count(landpatch%eindex==landelm%eindex(i))

            IF (allocated(locpxl)) deallocate(locpxl)
            allocate(locpxl(numpxl))

            ! get all patches' index that eindex is equal the i element
            locpxl = pack([(ipxl, ipxl=1, numpatch)], landpatch%eindex==landelm%eindex(i))
            ! the min index is the start of patch's index
            grid_patch_s(i) = minval(locpxl)
            ! the max index is the end of patch's index
            grid_patch_e(i) = maxval(locpxl)
         ENDDO

         DO i=1, numelm
            ipatch = grid_patch_s (i)

            IF (ipatch.le.0) CYCLE
            gridarea = 0

            DO WHILE (ipatch.le.grid_patch_e(i))

               !TODO-done: need to skip the 2m WMO patches
               IF (ipatch == wmo_patch(landpatch%ielm(ipatch)) ) THEN
                  ipatch = ipatch + 1
                  CYCLE
               ENDIF

               ! using this year patch mapping to aggregate the previous year land cover data
               CALL aggregation_request_data (landpatch, ipatch, grid_patch, zip = .true., &
                  area = area_one, data_i4_2d_in1 = lcdatafr, data_i4_2d_out1 = lcdatafr_one)

               ipxstt = landpatch%ipxstt(ipatch)
               ipxend = landpatch%ipxend(ipatch)

               IF (allocated(lcfrbuff)) deallocate(lcfrbuff)
               allocate(lcfrbuff(ipxstt:ipxend))
               lcfrbuff(:) = lcdatafr_one(:)

               IF (allocated(areabuff)) deallocate(areabuff)
               allocate(areabuff(ipxstt:ipxend))
               areabuff(:) = area_one(:)

               sum_areabuff = sum(areabuff)
               DO ipxl = ipxstt, ipxend
                  ! Transfer trace - the key codes to count for the source land cover types of LULCC
                  lccpct_patches(ipatch, lcfrbuff(ipxl)) = lccpct_patches(ipatch, lcfrbuff(ipxl)) &
                                                         + areabuff(ipxl) / sum_areabuff
                  lccpct_matrix (ipatch, lcfrbuff(ipxl)) = lccpct_matrix (ipatch, lcfrbuff(ipxl)) &
                                                         + areabuff(ipxl)
               ENDDO
               gridarea = gridarea + sum_areabuff
               ipatch = ipatch + 1
            ENDDO

            lccpct_matrix(grid_patch_s(i):grid_patch_e(i), :) = &
               lccpct_matrix (grid_patch_s(i):grid_patch_e(i), :) / gridarea

         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

      dir_landdata = DEF_dir_landdata
      CALL system('mkdir -p ' // trim(dir_landdata) // '/lulcc/' // trim(thisyr))
      DO ilc = 0, N_land_classification
         write(c2, '(i2.2)') ilc
         lndname = trim(dir_landdata)//'/lulcc/'//trim(thisyr)//&
            '/lccpct_patches_lc'//trim(c2)//'.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'lccpct_patches', 'patch', &
            landpatch, lccpct_patches(:,ilc), DEF_Srfdata_CompressLevel)
      ENDDO

#ifdef SrfdataDiag
      typindex = (/(ityp, ityp = 0, N_land_classification)/)
      lndname  = trim(dir_landdata) // &
                 '/diag/lccpct_matrix_' // trim(thisyr) // '.nc'
      DO ilc = 0, N_land_classification
         CALL srfdata_map_and_write (lccpct_matrix(:,ilc), landpatch%settyp, typindex, &
            m_patch2diag, -1.0e36_r8, lndname, 'lccpct_matrix', compress = 6, &
            write_mode = 'one', defval=0._r8, lastdimname = 'source_patch', lastdimvalue = ilc, &
            create_mode=first_call)
         IF ( first_call ) first_call = .false.
      ENDDO
      deallocate(typindex)
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('lccpct_patches', lccpct_patches)
      CALL check_vector_data ('lccpct_matrix' , lccpct_matrix )
#endif

      IF (p_is_worker) THEN
         IF (allocated(area_one))    deallocate (area_one)
      ENDIF

      CALL deallocate_LulccTransferTrace

   END SUBROUTINE MAKE_LulccTransferTrace



   SUBROUTINE deallocate_LulccTransferTrace
   ! --------------------------------------------------
   ! Deallocates memory for Lulcc time invariant variables
   ! --------------------------------------------------
   USE MOD_SPMD_Task

      IF (p_is_worker) THEN
         IF (allocated(lccpct_patches)) deallocate (lccpct_patches)
         IF (allocated(lccpct_matrix )) deallocate (lccpct_matrix )
      ENDIF

   END SUBROUTINE deallocate_LulccTransferTrace

END MODULE MOD_Lulcc_TransferTrace
! ---------- EOP ------------
