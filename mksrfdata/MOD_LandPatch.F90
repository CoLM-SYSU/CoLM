#include <define.h>

MODULE MOD_LandPatch

!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!    Build pixelset "landpatch".
!
!    In CoLM, the global/regional area is divided into a hierarchical structure:
!    1. If GRIDBASED or UNSTRUCTURED is defined, it is
!       ELEMENT >>> PATCH
!    2. If CATCHMENT is defined, it is
!       ELEMENT >>> HRU >>> PATCH
!    If Plant Function Type classification is used, PATCH is further divided into PFT.
!    If Plant Community classification is used,     PATCH is further divided into PC.
!
!    "landpatch" refers to pixelset PATCH.
!
!  Created by Shupeng Zhang, May 2023
!    porting codes from Hua Yuan's OpenMP version to MPI parallel version.
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_Pixelset
   USE MOD_Vars_Global
   USE MOD_Const_LC
   IMPLICIT NONE

   ! ---- Instance ----
   integer :: numpatch
   type(grid_type)     :: grid_patch
   type(pixelset_type) :: landpatch

   type(subset_type)   :: elm_patch
   type(superset_type) :: patch2elm

#ifdef CATCHMENT
   type(subset_type)   :: hru_patch
   type(superset_type) :: patch2hru
#endif


CONTAINS

   ! -------------------------------
   SUBROUTINE landpatch_build (lc_year)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Utils
   USE MOD_UserDefFun
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_Mesh
   USE MOD_LandElm
#ifdef CATCHMENT
   USE MOD_LandHRU
#endif
   USE MOD_Namelist
   USE MOD_NetCDFBlock
   USE MOD_5x5DataReadin
   USE MOD_AggregationRequestData

   IMPLICIT NONE

   integer, intent(in) :: lc_year
   ! Local Variables
   character(len=256) :: fname, dir
   character(len=255) :: cyear
   type (block_data_int32_2d) :: patchdata
   integer :: iloc, npxl, ipxl, numset
   integer :: ie, iset, ipxstt, ipxend
   integer,   allocatable :: types(:), order(:), ibuff(:)
   integer*8, allocatable :: eindex_tmp(:)
   integer,   allocatable :: settyp_tmp(:), ipxstt_tmp(:), ipxend_tmp(:), ielm_tmp(:)
   logical,   allocatable :: msk(:)
   integer :: npatch_glb
   integer :: dominant_type
   integer, allocatable :: npxl_types (:)

      write(cyear,'(i4.4)') lc_year
      IF (p_is_master) THEN
         write(*,'(A)') 'Making land patches:'
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_io) THEN

         CALL allocate_block_data (grid_patch, patchdata)

#ifndef LULC_USGS
         dir = trim(DEF_dir_rawdata) // trim(DEF_rawdata%landcover%dir)

         ! add parameter input for time year
         IF (DEF_rawdata%landcover%fname == "landtype-igbp-modis") THEN
            fname = trim(dir)//&
               trim(DEF_rawdata%landcover%fname)//'-'//trim(cyear)//'.nc'
            CALL ncio_read_block (fname, 'landtype', grid_patch, patchdata)
         ELSE
            fname = trim(DEF_rawdata%landcover%fname)//'.'//trim(cyear)
            CALL read_5x5_data (dir, fname, grid_patch, 'LC', patchdata)
         ENDIF
#else
         !TODO: need usgs land cover type data
         fname = trim(DEF_dir_rawdata)//trim(DEF_rawdata%landcover%dir)//&
                 '/landtype-usgs-update.nc'
         CALL ncio_read_block (fname, 'landtype', grid_patch, patchdata)
#endif

#ifdef USEMPI
         CALL aggregation_data_daemon (grid_patch, data_i4_2d_in1 = patchdata)
#endif
      ENDIF

      IF (p_is_worker) THEN

#ifdef CATCHMENT
         numset = numhru
#else
         numset = numelm
#endif

         IF (numset > 0) THEN
            allocate (eindex_tmp (numset*N_land_classification))
            allocate (settyp_tmp (numset*N_land_classification))
            allocate (ipxstt_tmp (numset*N_land_classification))
            allocate (ipxend_tmp (numset*N_land_classification))
            allocate (ielm_tmp   (numset*N_land_classification))
         ENDIF

         numpatch = 0

         DO iset = 1, numset
#ifdef CATCHMENT
            ie     = landhru%ielm  (iset)
            ipxstt = landhru%ipxstt(iset)
            ipxend = landhru%ipxend(iset)
#else
            ie     = landelm%ielm  (iset)
            ipxstt = landelm%ipxstt(iset)
            ipxend = landelm%ipxend(iset)
#endif

            npxl   = ipxend - ipxstt + 1

            allocate (types (ipxstt:ipxend))

#ifdef CATCHMENT
            CALL aggregation_request_data (landhru, iset, grid_patch, zip = .false., &
#else
            CALL aggregation_request_data (landelm, iset, grid_patch, zip = .false., &
#endif
               data_i4_2d_in1 = patchdata, data_i4_2d_out1 = ibuff)

            IF ( DEF_USE_GLC30 .or. DEF_USE_ESACCI ) THEN
               types(:) = LC_Map(ibuff)
            ELSE
               types(:) = ibuff
            ENDIF

            deallocate (ibuff)

#ifdef CATCHMENT
            IF (landhru%settyp(iset) <= 0) THEN
               types(ipxstt:ipxend) = WATERBODY
            ENDIF
            WHERE (types == 0)
               ! set land in MERITHydro while ocean in landtype data as water body
               types = WATERBODY
            END WHERE
            WHERE (types == 11)
               types = 10
            END WHERE
#endif

            IF ((DEF_USE_PFT .and. (.not. DEF_SOLO_PFT)) .or. DEF_FAST_PC) THEN
               ! For classification of plant function types or fast PC,
               ! merge all land types with soil ground
               DO ipxl = ipxstt, ipxend
                  IF (types(ipxl) > 0) THEN
                     IF (patchtypes(types(ipxl)) == 0) THEN
#ifndef CROP
                        ! Deal with cropland separately for fast PC
                        IF (DEF_FAST_PC .and. &
                           (types(ipxl)==CROPLAND .or. types(ipxl)==14)) THEN
                           types(ipxl) = CROPLAND
                        ELSE
                           types(ipxl) = 1
                        ENDIF
#else
                        types(ipxl) = 1
#endif
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF

            allocate (order (ipxstt:ipxend))
            order = (/ (ipxl, ipxl = ipxstt, ipxend) /)

            CALL quicksort (npxl, types, order)

            mesh(ie)%ilon(ipxstt:ipxend) = mesh(ie)%ilon(order)
            mesh(ie)%ilat(ipxstt:ipxend) = mesh(ie)%ilat(order)

            IF (DEF_USE_DOMINANT_PATCHTYPE) THEN
               allocate (npxl_types (0:maxval(types)))
               npxl_types(:) = 0
               DO ipxl = ipxstt, ipxend
                  npxl_types(types(ipxl)) = npxl_types(types(ipxl)) + 1
               ENDDO

               IF (any(types > 0)) THEN
                  iloc = findloc_ud(types > 0) + ipxstt - 1
                  dominant_type = maxloc(npxl_types(1:), dim=1)
                  types(iloc:ipxend) = dominant_type
               ENDIF

               deallocate(npxl_types)
            ENDIF

            DO ipxl = ipxstt, ipxend
               IF (ipxl == ipxstt) THEN
                  numpatch = numpatch + 1
                  eindex_tmp(numpatch) = mesh(ie)%indx
                  settyp_tmp(numpatch) = types(ipxl)
                  ipxstt_tmp(numpatch) = ipxl
                  ielm_tmp  (numpatch) = ie
               ELSEIF (types(ipxl) /= types(ipxl-1)) THEN
                  ipxend_tmp(numpatch) = ipxl - 1

                  numpatch = numpatch + 1
                  eindex_tmp(numpatch) = mesh(ie)%indx
                  settyp_tmp(numpatch) = types(ipxl)
                  ipxstt_tmp(numpatch) = ipxl
                  ielm_tmp  (numpatch) = ie
               ENDIF
            ENDDO
            ipxend_tmp(numpatch) = ipxend

            deallocate (types)
            deallocate (order)

         ENDDO

         IF (numpatch > 0) THEN
            allocate (landpatch%eindex (numpatch))
            allocate (landpatch%settyp (numpatch))
            allocate (landpatch%ipxstt (numpatch))
            allocate (landpatch%ipxend (numpatch))
            allocate (landpatch%ielm   (numpatch))

            landpatch%eindex = eindex_tmp(1:numpatch)
            landpatch%ipxstt = ipxstt_tmp(1:numpatch)
            landpatch%ipxend = ipxend_tmp(1:numpatch)
            landpatch%settyp = settyp_tmp(1:numpatch)
            landpatch%ielm   = ielm_tmp  (1:numpatch)
         ENDIF

         IF (numset > 0) THEN
            deallocate (eindex_tmp)
            deallocate (ipxstt_tmp)
            deallocate (ipxend_tmp)
            deallocate (settyp_tmp)
            deallocate (ielm_tmp  )
         ENDIF

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif

      ENDIF

      landpatch%nset = numpatch

      CALL landpatch%set_vecgs


#if (!defined(URBAN_MODEL) && !defined(CROP))
#ifdef USEMPI
      IF (p_is_worker) THEN
         CALL mpi_reduce (numpatch, npatch_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', npatch_glb, ' patches.'
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
#endif

   END SUBROUTINE landpatch_build

   ! -----
   SUBROUTINE write_patchfrac (dir_landdata, lc_year)

   USE MOD_Namelist
   USE MOD_NetCDFVector
   IMPLICIT NONE

   integer, intent(in) :: lc_year
   character(len=*), intent(in) :: dir_landdata
   character(len=256) :: lndname, cyear

      write(cyear,'(i4.4)') lc_year
      CALL system('mkdir -p ' // trim(dir_landdata) // '/landpatch/' // trim(cyear))

      lndname = trim(dir_landdata)//'/landpatch/'//trim(cyear)//'/patchfrac_elm.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'patchfrac_elm', 'patch', landpatch, elm_patch%subfrc, DEF_Srfdata_CompressLevel)

#ifdef CATCHMENT
      lndname = trim(dir_landdata)//'/landpatch/'//trim(cyear)//'/patchfrac_hru.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'patchfrac_hru', 'patch', landpatch, hru_patch%subfrc, DEF_Srfdata_CompressLevel)
#endif

   END SUBROUTINE write_patchfrac

END MODULE MOD_LandPatch
