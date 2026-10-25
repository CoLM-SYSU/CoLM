#include <define.h>

SUBROUTINE Aggregation_TopoWetness ( &
      gridtwi, dir_rawdata, dir_model_landdata, lc_year)
!-----------------------------------------------------------------------
!  Topographic Wetness Index data ln(a/tanB), calculated by using
!
!   Yamazaki, D., Ikeshima, D., Sosa, J.,Bates, P. D., Allen, G. H.,
!   Pavelsky, T. M. (2019).
!   MERIT Hydro: ahigh‐resolution global hydrographymap based on
!   latest topography dataset.Water Resources Research, 55, 5053-5073.
!
!  Created by Shupeng Zhang, 06/2025
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_Mesh
   USE MOD_LandElm
   USE MOD_LandPatch
   USE MOD_Land2mWMO
   USE MOD_NetCDFVector
   USE MOD_NetCDFBlock
   USE MOD_CatchmentDataReadin
#ifdef RangeCheck
   USE MOD_RangeCheck
#endif
   USE MOD_AggregationRequestData
   USE MOD_Utils
#ifdef SrfdataDiag
   USE MOD_SrfdataDiag
#endif

   IMPLICIT NONE
   ! arguments:
   integer, intent(in) :: lc_year
   type(grid_type),  intent(in) :: gridtwi
   character(len=*), intent(in) :: dir_rawdata
   character(len=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ---------------------------------------------------------------
   character(len=256) :: landdir, lndname, cyear
   integer  :: ipatch, npxl, i, im, ielm, istt, iend
   integer  :: wmo_src
   real(r8) :: mean_twi, sigma_twi, skew_twi, fsatmax, fsatdcf, alp_twi, chi_twi, mu_twi

   type (block_data_real8_3d) :: twi

   real(r8), allocatable :: twi_one (:,:)
   logical,  allocatable :: mask    (:,:)

   real(r8), allocatable :: mean_twi_patches(:)
   real(r8), allocatable :: fsatmax_patches (:), fsatdcf_patches (:)
   real(r8), allocatable :: alp_twi_patches (:), chi_twi_patches (:), mu_twi_patches (:)

   real(r8), allocatable :: twi_sort(:), xx(:), yy(:)
   integer,  allocatable :: order   (:)

#ifdef SrfdataDiag
   integer :: typpatch(N_land_classification+1), ityp
#endif

      write(cyear,'(i4.4)') lc_year
      landdir = trim(dir_model_landdata) // '/topography/' // trim(cyear)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,'(/, A)') 'Aggregate topographic wetness index ...'
         CALL system('mkdir -p ' // trim(adjustl(landdir)))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

! -----------------------------------------------------------------------------------------
!   aggregate the topographic index from the resolution of raw data to modelling resolution
! -----------------------------------------------------------------------------------------

      IF (p_is_io) THEN

         lndname = trim(dir_rawdata)//'/TWI.nc'
         CALL allocate_block_data (gridtwi, twi, 25)
         CALL ncio_read_block (lndname, 'twi', gridtwi, 25, twi)

#ifdef USEMPI
         CALL aggregation_data_daemon (gridtwi, data_r8_3d_in1 = twi, n1_r8_3d_in1 = 25)
#endif
      ENDIF

      IF (p_is_worker) THEN

         allocate (mean_twi_patches (numpatch));  mean_twi_patches(:) = spval
         allocate (fsatmax_patches  (numpatch));  fsatmax_patches (:) = spval
         allocate (fsatdcf_patches  (numpatch));  fsatdcf_patches (:) = spval
         allocate (alp_twi_patches  (numpatch));  alp_twi_patches (:) = spval
         allocate (chi_twi_patches  (numpatch));  chi_twi_patches (:) = spval
         allocate (mu_twi_patches   (numpatch));  mu_twi_patches  (:) = spval

         DO ipatch = 1, numpatch

            IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
               wmo_src = wmo_source (landpatch%ielm(ipatch))

               fsatmax_patches(ipatch) = fsatmax_patches(wmo_src)
               fsatdcf_patches(ipatch) = fsatdcf_patches(wmo_src)
               alp_twi_patches(ipatch) = alp_twi_patches(wmo_src)
               chi_twi_patches(ipatch) = chi_twi_patches(wmo_src)
               mu_twi_patches (ipatch) = mu_twi_patches (wmo_src)

               CYCLE
            ENDIF

            CALL aggregation_request_data (landpatch, ipatch, gridtwi, zip = USE_zip_for_aggregation, &
               data_r8_3d_in1 = twi, data_r8_3d_out1 = twi_one, n1_r8_3d_in1 = 25, lb1_r8_3d_in1 = 1  )

            allocate (mask (size(twi_one,1),size(twi_one,2)))
            mask = (twi_one > -1.e3)
            npxl = count(mask)

            IF (npxl >= 25) THEN

               allocate (twi_sort (npxl))
               allocate (order    (npxl))

               twi_sort = pack(twi_one, mask)

               CALL quicksort (npxl, twi_sort, order)

               mean_twi_patches(ipatch) = sum(twi_sort) / npxl

               im = 1
               DO WHILE ((twi_sort(im) < mean_twi_patches(ipatch)) .and. (im < npxl-1))
                  im = im + 1
               ENDDO

               fsatmax_patches(ipatch) = 1 - real(im-1)/npxl

               allocate (xx (npxl-im))
               allocate (yy (npxl-im))

               xx = -(twi_sort(im:(npxl-1)) - mean_twi_patches(ipatch))
               yy = (/(log((1-real(i)/npxl)/fsatmax_patches(ipatch)), i = im, npxl-1)/)

               fsatdcf_patches(ipatch) = sum(xx*yy)/sum(xx**2)

               mean_twi = mean_twi_patches(ipatch)
               sigma_twi = sqrt(sum((twi_sort-mean_twi)**2) / (npxl-1))

               IF (sigma_twi > 0) THEN
                  skew_twi  = real(npxl)/(real(npxl-1)*real(npxl-2)) &
                     * sum((twi_sort-mean_twi)**3) / sigma_twi**3
                  IF (skew_twi > 0) THEN
                     alp_twi_patches(ipatch) = (2./skew_twi)**2
                     chi_twi_patches(ipatch) = sigma_twi*skew_twi/2
                     mu_twi_patches (ipatch) = mean_twi - 2.*sigma_twi/skew_twi
                  ENDIF
               ENDIF

               fsatmax_patches(ipatch) = min(max(fsatmax_patches(ipatch),  0.1), 1.0)
               fsatdcf_patches(ipatch) = min(max(fsatdcf_patches(ipatch),  0.2), 2.0)
               alp_twi_patches(ipatch) = min(max(alp_twi_patches(ipatch),  0.1), 50.)
               chi_twi_patches(ipatch) = min(max(chi_twi_patches(ipatch), 0.01), 4.0)
               mu_twi_patches (ipatch) = min(max(mu_twi_patches (ipatch),   0.), 12.)

               deallocate (twi_sort)
               deallocate (order   )
               deallocate (xx      )
               deallocate (yy      )

            ENDIF

            deallocate (mask)

         ENDDO

         DO ielm = 1, numelm

            ! default value is set as global median at 0.1 degree resolution
            mean_twi = 9.27
            fsatmax  = 0.38
            fsatdcf  = 0.55
            alp_twi  = 1.34
            chi_twi  = 1.61
            mu_twi   = 6.95

            CALL aggregation_request_data (landelm, ielm, gridtwi, zip = USE_zip_for_aggregation,     &
               data_r8_3d_in1 = twi, data_r8_3d_out1 = twi_one, n1_r8_3d_in1 = 25, lb1_r8_3d_in1 = 1  )

            allocate (mask (size(twi_one,1),size(twi_one,2)))
            mask = (twi_one > -1.e3)
            npxl = count(mask)

            IF (npxl >= 25) THEN

               allocate (twi_sort (npxl))
               allocate (order    (npxl))

               twi_sort = pack(twi_one, mask)

               CALL quicksort (npxl, twi_sort, order)

               mean_twi = sum(twi_sort) / npxl

               im = 1
               DO WHILE ((twi_sort(im) < mean_twi) .and. (im < npxl-1))
                  im = im + 1
               ENDDO

               fsatmax = 1 - real(im-1)/npxl

               allocate (xx (npxl-im))
               allocate (yy (npxl-im))

               xx = -(twi_sort(im:(npxl-1)) - mean_twi)
               yy = (/(log((1-real(i)/npxl)/fsatmax), i = im, npxl-1)/)

               fsatdcf = sum(xx*yy)/sum(xx**2)

               sigma_twi = sqrt(sum((twi_sort-mean_twi)**2) / (npxl-1))

               IF (sigma_twi > 0) THEN
                  skew_twi  = real(npxl)/(real(npxl-1)*real(npxl-2)) &
                     * sum((twi_sort-mean_twi)**3) / sigma_twi**3
                  IF (skew_twi > 0) THEN
                     alp_twi = (2./skew_twi)**2
                     chi_twi = sigma_twi*skew_twi/2
                     mu_twi  = mean_twi - 2.*sigma_twi/skew_twi
                  ENDIF
               ENDIF

               fsatmax = min(max(fsatmax,  0.1), 1.0)
               fsatdcf = min(max(fsatdcf,  0.2), 2.0)
               alp_twi = min(max(alp_twi,  0.1), 50.)
               chi_twi = min(max(chi_twi, 0.01), 4.0)
               mu_twi  = min(max(mu_twi ,   0.), 12.)

               deallocate (twi_sort)
               deallocate (order   )
               deallocate (xx      )
               deallocate (yy      )

            ENDIF

            istt = elm_patch%substt(ielm)
            iend = elm_patch%subend(ielm)
            DO ipatch = istt, iend
               IF (mean_twi_patches(ipatch) == spval)   mean_twi_patches(ipatch) = mean_twi
               IF (fsatmax_patches (ipatch) == spval)   fsatmax_patches (ipatch) = fsatmax
               IF (fsatdcf_patches (ipatch) == spval)   fsatdcf_patches (ipatch) = fsatdcf
               IF (alp_twi_patches (ipatch) == spval)   alp_twi_patches (ipatch) = alp_twi
               IF (chi_twi_patches (ipatch) == spval)   chi_twi_patches (ipatch) = chi_twi
               IF (mu_twi_patches  (ipatch) == spval)   mu_twi_patches  (ipatch) = mu_twi
            ENDDO

            deallocate (mask)

         ENDDO


#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('mean_twi_patches  ', mean_twi_patches)
      CALL check_vector_data ('fsatmax_patches   ', fsatmax_patches )
      CALL check_vector_data ('fsatdcf_patches   ', fsatdcf_patches )
      CALL check_vector_data ('alp_twi_patches   ', alp_twi_patches )
      CALL check_vector_data ('chi_twi_patches   ', chi_twi_patches )
      CALL check_vector_data ('mu_twi_patches    ', mu_twi_patches  )
#endif

      lndname = trim(landdir)//'/mean_twi_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'mean_twi_patches', 'patch', landpatch, &
         mean_twi_patches, DEF_Srfdata_CompressLevel)

      lndname = trim(landdir)//'/fsatmax_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'fsatmax_patches', 'patch', landpatch, &
         fsatmax_patches, DEF_Srfdata_CompressLevel)

      lndname = trim(landdir)//'/fsatdcf_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'fsatdcf_patches', 'patch', landpatch, &
         fsatdcf_patches, DEF_Srfdata_CompressLevel)

      lndname = trim(landdir)//'/alp_twi_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'alp_twi_patches', 'patch', landpatch, &
         alp_twi_patches, DEF_Srfdata_CompressLevel)

      lndname = trim(landdir)//'/chi_twi_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'chi_twi_patches', 'patch', landpatch, &
         chi_twi_patches, DEF_Srfdata_CompressLevel)

      lndname = trim(landdir)//'/mu_twi_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'mu_twi_patches', 'patch', landpatch, &
         mu_twi_patches, DEF_Srfdata_CompressLevel)


#ifdef SrfdataDiag
      typpatch = (/(ityp, ityp = 0, N_land_classification)/)

      lndname  = trim(dir_model_landdata) // '/diag/twi_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (mean_twi_patches, landpatch%settyp, typpatch, m_patch2diag, &
         spval, lndname, 'mean_twi', compress = 1, write_mode = 'one', create_mode=.true.)
      CALL srfdata_map_and_write (fsatmax_patches,  landpatch%settyp, typpatch, m_patch2diag, &
         spval, lndname, 'fsatmax',  compress = 1, write_mode = 'one')
      CALL srfdata_map_and_write (fsatdcf_patches,  landpatch%settyp, typpatch, m_patch2diag, &
         spval, lndname, 'fsatdcf',  compress = 1, write_mode = 'one')
      CALL srfdata_map_and_write (alp_twi_patches,  landpatch%settyp, typpatch, m_patch2diag, &
         spval, lndname, 'alp_twi',  compress = 1, write_mode = 'one')
      CALL srfdata_map_and_write (chi_twi_patches,  landpatch%settyp, typpatch, m_patch2diag, &
         spval, lndname, 'chi_twi',  compress = 1, write_mode = 'one')
      CALL srfdata_map_and_write (mu_twi_patches,   landpatch%settyp, typpatch, m_patch2diag, &
         spval, lndname, 'mu_twi',   compress = 1, write_mode = 'one')
#endif

      IF (p_is_worker) THEN
         deallocate ( mean_twi_patches )
         deallocate ( fsatmax_patches  )
         deallocate ( fsatdcf_patches  )
         deallocate ( alp_twi_patches  )
         deallocate ( chi_twi_patches  )
         deallocate ( mu_twi_patches   )
      ENDIF

END SUBROUTINE Aggregation_TopoWetness
