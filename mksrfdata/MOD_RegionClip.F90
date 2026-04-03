#include <define.h>

MODULE MOD_RegionClip
!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!    This module includes subroutines to clip surface data from an existing
!    data in a larger region.
!
!    Please use namelist variable "USE_srfdata_from_larger_region" to call
!    these subroutines.
!
!  Created by Shupeng Zhang, May 2023
!-----------------------------------------------------------------------

CONTAINS

   ! ----- region clip -----
   SUBROUTINE srfdata_region_clip (dir_landdata_in, dir_landdata_out)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Utils
   USE MOD_Pixel
   USE MOD_Block
   USE MOD_NetCDFSerial

   IMPLICIT NONE

   character(len=*), intent(in) :: dir_landdata_in, dir_landdata_out

   ! Local Variables
   real(r8) :: edges, edgen, edgew, edgee
   character(len=256) :: file_in, file_out, fileblock
   integer :: iproc, iblk, jblk, ie, ipxl, ilon, ilat, i1
   integer :: nelm_in, nelm_out
   integer,   allocatable :: nelm_blk(:,:), IOproc(:,:)
   integer*8, allocatable :: elmindex(:)
   integer,   allocatable :: elmnpxl(:), elmpixels(:,:,:)

   logical, allocatable :: elmmask  (:)
   logical, allocatable :: patchmask(:)
#ifdef CATCHMENT
   logical, allocatable :: hrumask  (:)
#endif
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   logical, allocatable :: pftmask  (:)
#endif

   integer :: month, YY, itime, Julian_day, nsl
   integer :: start_year, end_year
   character(len=1) :: c1
   character(len=2) :: c2, cx
   character(len=3) :: c3
   character(len=4) :: cyear


      IF (p_is_master) THEN
         IF (trim(adjustl(dir_landdata_in)) == trim(adjustl(dir_landdata_out))) THEN
            write(*,*) 'Warning : Surface Data from ', trim(adjustl(dir_landdata_in)), &
               ' to ', trim(adjustl(dir_landdata_out))
            write(*,*) 'Data will be overwritten, please CHANGE directory.'
#ifdef USEMPI
            CALL mpi_abort (p_comm_glb, 1, p_err)
#else
            CALL EXIT(1)
#endif
         ELSE
            CALL system('mkdir -p ' // trim(dir_landdata_out))
         ENDIF
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      edges = DEF_domain%edges
      edgen = DEF_domain%edgen
      edgew = DEF_domain%edgew
      edgee = DEF_domain%edgee

      CALL normalize_longitude (edgew)
      CALL normalize_longitude (edgee)

      CALL pixel%load_from_file  (dir_landdata_in)

      file_in = trim(dir_landdata_in) // '/block.nc'
      CALL ncio_read_bcast_serial (file_in, 'lat_s', gblock%lat_s)
      CALL ncio_read_bcast_serial (file_in, 'lat_n', gblock%lat_n)
      CALL ncio_read_bcast_serial (file_in, 'lon_w', gblock%lon_w)
      CALL ncio_read_bcast_serial (file_in, 'lon_e', gblock%lon_e)

      gblock%nyblk = size(gblock%lat_s)
      gblock%nxblk = size(gblock%lon_w)

      IF (p_is_master) THEN

         file_in  = trim(dir_landdata_in ) // '/mesh/mesh.nc'
         file_out = trim(dir_landdata_out) // '/mesh/mesh.nc'
         CALL system('cp ' // trim(file_in) // ' ' // trim(file_out))

         file_in  = trim(dir_landdata_in ) // '/pixel.nc'
         file_out = trim(dir_landdata_out) // '/pixel.nc'
         CALL system('cp ' // trim(file_in) // ' ' // trim(file_out))

         file_in  = trim(dir_landdata_in ) // '/block.nc'
         file_out = trim(dir_landdata_out) // '/block.nc'
         CALL system('cp ' // trim(file_in) // ' ' // trim(file_out))

      ENDIF

      file_in = trim(dir_landdata_in) // '/mesh/mesh.nc'
      CALL ncio_read_bcast_serial (file_in, 'nelm_blk', nelm_blk)

      allocate (IOproc (gblock%nxblk, gblock%nyblk))
      IOproc (:,:) = -1

      iproc = -1
      DO jblk = 1, gblock%nyblk
         IF ((gblock%lat_s(jblk) >= edgen) .or. (gblock%lat_n(jblk) <= edges)) THEN
            CYCLE
         ENDIF

         DO iblk = 1, gblock%nxblk
            IF ((gblock%lon_w(iblk) == gblock%lon_e(iblk)) .or. (edgew == edgee)) THEN
               IF (nelm_blk(iblk,jblk) > 0) THEN
                  iproc  = mod(iproc + 1, p_np_glb)
                  IOproc(iblk,jblk) = iproc
               ENDIF
            ELSE
               IF (    lon_between_floor(gblock%lon_w(iblk), edgew, edgee) &
                  .or. lon_between_ceil (gblock%lon_e(iblk), edgew, edgee) &
                  .or. lon_between_floor(edgew, gblock%lon_w(iblk), gblock%lon_e(iblk)) &
                  .or. lon_between_ceil (edgee, gblock%lon_w(iblk), gblock%lon_e(iblk))) THEN
                  IF (nelm_blk(iblk,jblk) > 0) THEN
                     iproc  = mod(iproc + 1, p_np_glb)
                     IOproc(iblk,jblk) = iproc
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      WHERE (IOproc /= p_iam_glb)
         nelm_blk = 0
      END WHERE

      DO jblk = 1, gblock%nyblk
         DO iblk = 1, gblock%nxblk
            IF (IOproc(iblk,jblk) == p_iam_glb) THEN

               file_in = trim(dir_landdata_in) // '/mesh/mesh.nc'
               CALL get_filename_block (file_in, iblk, jblk, fileblock)

               CALL ncio_read_serial (fileblock, 'elmindex' , elmindex )
               CALL ncio_read_serial (fileblock, 'elmnpxl'  , elmnpxl  )
               CALL ncio_read_serial (fileblock, 'elmpixels', elmpixels)

               nelm_in = nelm_blk(iblk,jblk)

               allocate (elmmask (nelm_in))
               elmmask(:) = .false.

               DO ie = 1, nelm_in
                  DO ipxl = 1, elmnpxl(ie)
                     ilon = elmpixels(1,ipxl,ie)
                     ilat = elmpixels(2,ipxl,ie)

                     IF (((pixel%lat_s(ilat) < edgen) .and. (pixel%lat_n(ilat) > edges)) &
                        .and. (lon_between_floor(pixel%lon_w(ilon), edgew, edgee) &
                        .or. lon_between_ceil(pixel%lon_e(ilon), edgew, edgee))) THEN

                        elmmask(ie) = .true.

                        EXIT
                     ELSE
                        CYCLE
                     ENDIF
                  ENDDO
               ENDDO

               nelm_out = count(elmmask)
               nelm_blk(iblk,jblk) = nelm_out

               IF (nelm_out > 0) THEN

                  CALL clip_pixelset (dir_landdata_in, 'landpatch', iblk, jblk, elmmask, elmindex, patchmask)
#ifdef CATCHMENT
                  CALL clip_pixelset (dir_landdata_in, 'landhru'  , iblk, jblk, elmmask, elmindex, hrumask  )
#endif
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
                  CALL clip_pixelset (dir_landdata_in, 'landpft'  , iblk, jblk, elmmask, elmindex, pftmask  )
#endif
                  CALL system('mkdir -p ' // trim(dir_landdata_out) // '/mesh')
                  file_in  = trim(dir_landdata_in)  // '/mesh/mesh.nc'
                  file_out = trim(dir_landdata_out) // '/mesh/mesh.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'elmindex' , elmmask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'elmnpxl'  , elmmask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'elmpixels', elmmask)

                  CALL system('mkdir -p ' // trim(dir_landdata_out) // '/landelm')
                  file_in  = trim(dir_landdata_in)  // '/landelm/landelm.nc'
                  file_out = trim(dir_landdata_out) // '/landelm/landelm.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'eindex', elmmask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'ipxstt', elmmask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'ipxend', elmmask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'settyp', elmmask)

                  CALL system('mkdir -p ' // trim(dir_landdata_out) // '/landpatch')
                  file_in  = trim(dir_landdata_in)  // '/landpatch/landpatch.nc'
                  file_out = trim(dir_landdata_out) // '/landpatch/landpatch.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'eindex', patchmask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'ipxstt', patchmask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'ipxend', patchmask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'settyp', patchmask)

                  file_in  = trim(dir_landdata_in)  // '/landpatch/patchfrac_elm.nc'
                  file_out = trim(dir_landdata_out) // '/landpatch/patchfrac_elm.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'patchfrac_elm', patchmask)

#ifdef CATCHMENT
                  CALL system('mkdir -p ' // trim(dir_landdata_out) // '/landhru')
                  file_in  = trim(dir_landdata_in)  // '/landhru/landhru.nc'
                  file_out = trim(dir_landdata_out) // '/landhru/landhru.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'eindex', hrumask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'ipxstt', hrumask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'ipxend', hrumask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'settyp', hrumask)

                  file_in  = trim(dir_landdata_in)  // '/landpatch/patchfrac_hru.nc'
                  file_out = trim(dir_landdata_out) // '/landpatch/patchfrac_hru.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'patchfrac_hru', patchmask)
#endif

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
                  CALL system('mkdir -p ' // trim(dir_landdata_out) // '/landpft')
                  file_in  = trim(dir_landdata_in)  // '/landpft/landpft.nc'
                  file_out = trim(dir_landdata_out) // '/landpft/landpft.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'eindex', pftmask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'ipxstt', pftmask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'ipxend', pftmask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'settyp', pftmask)
#endif
               ENDIF

               ! Leaf Area Index
               CALL system('mkdir -p ' // trim(dir_landdata_out) // '/LAI')

               IF (DEF_LAI_CHANGE_YEARLY) THEN
                  start_year = DEF_simulation_time%start_year
                  end_year   = DEF_simulation_time%end_year
               ELSE
                  start_year = DEF_LC_YEAR
                  end_year   = DEF_LC_YEAR
               ENDIF

               DO YY = start_year, end_year

                  write(cyear,'(i4.4)') YY
                  CALL system('mkdir -p ' // trim(dir_landdata_out) // '/LAI/' // trim(cyear))

                  IF (DEF_LAI_MONTHLY) THEN
                     DO month = 1, 12
                        write(c2,'(i2.2)') month

                        file_in  = trim(dir_landdata_in)  // '/LAI/' // trim(cyear) // '/LAI_patches' // trim(c2) // '.nc'
                        file_out = trim(dir_landdata_out) // '/LAI/' // trim(cyear) // '/LAI_patches' // trim(c2) // '.nc'
                        CALL clip_vector (file_in, file_out, iblk, jblk, 'LAI_patches', patchmask)

                        file_in  = trim(dir_landdata_in)  // '/LAI/' // trim(cyear) // '/SAI_patches' // trim(c2) // '.nc'
                        file_out = trim(dir_landdata_out) // '/LAI/' // trim(cyear) // '/SAI_patches' // trim(c2) // '.nc'
                        CALL clip_vector (file_in, file_out, iblk, jblk, 'SAI_patches', patchmask)
                     ENDDO
                  ELSE
                     DO itime = 1, 46
                        Julian_day = 1 + (itime-1)*8
                        write(c3, '(i3.3)') Julian_day
                        file_in  = trim(dir_landdata_in) //'/LAI/'//trim(cyear)//'/LAI_patches'//trim(c3)// '.nc'
                        file_out = trim(dir_landdata_out)//'/LAI/'//trim(cyear)//'/LAI_patches'//trim(c3)// '.nc'

                        CALL clip_vector (file_in, file_out, iblk, jblk, 'LAI_patches', patchmask)
                     ENDDO
                  ENDIF

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
                  DO month = 1, 12
                     write(c2,'(i2.2)') month

                     file_in  = trim(dir_landdata_in) //'/LAI/'//trim(cyear)//'/LAI_pfts'//trim(c2)//'.nc'
                     file_out = trim(dir_landdata_out)//'/LAI/'//trim(cyear)//'/LAI_pfts'//trim(c2)//'.nc'
                     CALL clip_vector (file_in, file_out, iblk, jblk, 'LAI_pfts', pftmask)

                     file_in  = trim(dir_landdata_in) //'/LAI/'//trim(cyear)//'/SAI_pfts'//trim(c2)//'.nc'
                     file_out = trim(dir_landdata_out)//'/LAI/'//trim(cyear)//'/SAI_pfts'//trim(c2)//'.nc'
                     CALL clip_vector (file_in, file_out, iblk, jblk, 'SAI_pfts', pftmask)
                  ENDDO
#endif
               ENDDO

               ! depth to bedrock
               IF(DEF_USE_BEDROCK)THEN
                  CALL system('mkdir -p ' // trim(dir_landdata_out) // '/dbedrock')
                  file_in  = trim(dir_landdata_in)  // '/dbedrock/dbedrock_patches.nc'
                  file_out = trim(dir_landdata_out) // '/dbedrock/dbedrock_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'dbedrock_patches', patchmask)
               ENDIF

               ! forest height
               CALL system('mkdir -p ' // trim(dir_landdata_out) // '/htop')
               file_in  = trim(dir_landdata_in)  // '/htop/htop_patches.nc'
               file_out = trim(dir_landdata_out) // '/htop/htop_patches.nc'
               CALL clip_vector (file_in, file_out, iblk, jblk, 'htop_patches', patchmask)
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
               file_in  = trim(dir_landdata_in)  // '/htop/htop_pfts.nc'
               file_out = trim(dir_landdata_out) // '/htop/htop_pfts.nc'
               CALL clip_vector (file_in, file_out, iblk, jblk, 'htop_pfts', pftmask)
#endif

               ! lake depth
               CALL system('mkdir -p ' // trim(dir_landdata_out) // '/lakedepth')
               file_in  = trim(dir_landdata_in)  // '/lakedepth/lakedepth_patches.nc'
               file_out = trim(dir_landdata_out) // '/lakedepth/lakedepth_patches.nc'
               CALL clip_vector (file_in, file_out, iblk, jblk, 'lakedepth_patches', patchmask)

               ! plant function type percentage
               CALL system('mkdir -p ' // trim(dir_landdata_out) // '/pctpft')
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
               file_in  = trim(dir_landdata_in)  // '/pctpft/pct_pfts.nc'
               file_out = trim(dir_landdata_out) // '/pctpft/pct_pfts.nc'
               CALL clip_vector (file_in, file_out, iblk, jblk, 'pct_pfts', pftmask)
#ifdef CROP
               file_in  = trim(dir_landdata_in)  // '/pctpft/pct_crops.nc'
               file_out = trim(dir_landdata_out) // '/pctpft/pct_crops.nc'
               CALL clip_vector (file_in, file_out, iblk, jblk, 'pct_crops', patchmask)
#endif
#endif

               ! soil
               CALL system('mkdir -p ' // trim(dir_landdata_out) // '/soil')

               file_in  = trim(dir_landdata_in)  // '/soil/soil_s_v_alb_patches.nc'
               file_out = trim(dir_landdata_out) // '/soil/soil_s_v_alb_patches.nc'
               CALL clip_vector (file_in, file_out, iblk, jblk, 'soil_s_v_alb', patchmask)

               file_in  = trim(dir_landdata_in)  // '/soil/soil_d_v_alb_patches.nc'
               file_out = trim(dir_landdata_out) // '/soil/soil_d_v_alb_patches.nc'
               CALL clip_vector (file_in, file_out, iblk, jblk, 'soil_d_v_alb', patchmask)

               file_in  = trim(dir_landdata_in)  // '/soil/soil_s_n_alb_patches.nc'
               file_out = trim(dir_landdata_out) // '/soil/soil_s_n_alb_patches.nc'
               CALL clip_vector (file_in, file_out, iblk, jblk, 'soil_s_n_alb', patchmask)

               file_in  = trim(dir_landdata_in)  // '/soil/soil_d_n_alb_patches.nc'
               file_out = trim(dir_landdata_out) // '/soil/soil_d_n_alb_patches.nc'
               CALL clip_vector (file_in, file_out, iblk, jblk, 'soil_d_n_alb', patchmask)

               DO nsl = 1, 8

                  write(c1,'(i1)') nsl

                  ! (1) volumetric fraction of quartz within mineral soil
                  file_in  = trim(dir_landdata_in)  // '/soil/vf_quartz_mineral_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/vf_quartz_mineral_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'vf_quartz_mineral_s_l'//trim(c1)//'_patches', patchmask)

                  ! (2) volumetric fraction of gravels
                  file_in  = trim(dir_landdata_in)  // '/soil/vf_gravels_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/vf_gravels_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'vf_gravels_s_l'//trim(c1)//'_patches', patchmask)

                  ! (3) volumetric fraction of sand
                  file_in  = trim(dir_landdata_in)  // '/soil/vf_sand_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/vf_sand_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'vf_sand_s_l'//trim(c1)//'_patches', patchmask)

                  ! (4) volumetric fraction of organic matter with the parameter alpha and beta in the Balland V. and P. A. Arp (2005) model
                  file_in  = trim(dir_landdata_in)  // '/soil/vf_om_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/vf_om_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'vf_om_s_l'//trim(c1)//'_patches', patchmask)

                  file_in  = trim(dir_landdata_in)  // '/soil/BA_alpha_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/BA_alpha_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'BA_alpha_l'//trim(c1)//'_patches', patchmask)

                  file_in  = trim(dir_landdata_in)  // '/soil/BA_beta_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/BA_beta_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'BA_beta_l'//trim(c1)//'_patches', patchmask)

                  ! (5) gravimetric fraction of gravels
                  file_in  = trim(dir_landdata_in)  // '/soil/wf_gravels_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/wf_gravels_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'wf_gravels_s_l'//trim(c1)//'_patches', patchmask)

                  ! (6) gravimetric fraction of sand
                  file_in  = trim(dir_landdata_in)  // '/soil/wf_sand_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/wf_sand_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'wf_sand_s_l'//trim(c1)//'_patches', patchmask)

#ifdef vanGenuchten_Mualem_SOIL_MODEL
                  ! (7) VGM's pore-connectivity parameter (L)
                  file_in  = trim(dir_landdata_in)  // '/soil/L_vgm_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/L_vgm_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'L_vgm_l'//trim(c1)//'_patches', patchmask)

                  ! (8) VGM's residual water content (theta_r) [cm3/cm3]
                  file_in  = trim(dir_landdata_in)  // '/soil/theta_r_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/theta_r_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'theta_r_l'//trim(c1)//'_patches', patchmask)

                  ! (9) VGM's parameter corresponding approximately to the inverse of the air-entry value (alpha)
                  file_in  = trim(dir_landdata_in)  // '/soil/alpha_vgm_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/alpha_vgm_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'alpha_vgm_l'//trim(c1)//'_patches', patchmask)

                  ! (10) VGM's shape parameter (n)
                  file_in  = trim(dir_landdata_in)  // '/soil/n_vgm_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/n_vgm_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'n_vgm_l'//trim(c1)//'_patches', patchmask)
#endif
                  ! (11) saturated water content [cm3/cm3]
                  file_in  = trim(dir_landdata_in)  // '/soil/theta_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/theta_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'theta_s_l'//trim(c1)//'_patches', patchmask)

                  ! (12) matric potential at saturation (psi_s) [cm]
                  file_in  = trim(dir_landdata_in)  // '/soil/psi_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/psi_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'psi_s_l'//trim(c1)//'_patches', patchmask)

                  ! (13) pore size distribution index [dimensionless]
                  file_in  = trim(dir_landdata_in)  // '/soil/lambda_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/lambda_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'lambda_l'//trim(c1)//'_patches', patchmask)

                  ! (14) saturated hydraulic conductivity [cm/day]
                  file_in  = trim(dir_landdata_in)  // '/soil/k_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/k_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'k_s_l'//trim(c1)//'_patches', patchmask)

                  ! (15) heat capacity of soil solids [J/(m3 K)]
                  file_in  = trim(dir_landdata_in)  // '/soil/csol_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/csol_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'csol_l'//trim(c1)//'_patches', patchmask)

                  ! (16) thermal conductivity of unfrozen saturated soil [W/m-K]
                  file_in  = trim(dir_landdata_in)  // '/soil/tksatu_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/tksatu_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'tksatu_l'//trim(c1)//'_patches', patchmask)

                  ! (17) thermal conductivity of frozen saturated soil [W/m-K]
                  file_in  = trim(dir_landdata_in)  // '/soil/tksatf_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/tksatf_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'tksatf_l'//trim(c1)//'_patches', patchmask)

                  ! (18) thermal conductivity for dry soil [W/(m-K)]
                  file_in  = trim(dir_landdata_in)  // '/soil/tkdry_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/tkdry_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'tkdry_l'//trim(c1)//'_patches', patchmask)

                  ! (19) thermal conductivity of soil solids [W/m-K]
                  file_in  = trim(dir_landdata_in)  // '/soil/k_solids_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/k_solids_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'k_solids_l'//trim(c1)//'_patches', patchmask)

                  ! (20) OM_density [kg/m3]
                  file_in  = trim(dir_landdata_in)  // '/soil/OM_density_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/OM_density_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'OM_density_s_l'//trim(c1)//'_patches', patchmask)

                  ! (21) bulk density of soil (GRAVELS + OM + Mineral Soils)
                  file_in  = trim(dir_landdata_in)  // '/soil/BD_all_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/BD_all_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'BD_all_s_l'//trim(c1)//'_patches', patchmask)

                  ! (22) volumetric fraction of clay
                  file_in  = trim(dir_landdata_in)  // '/soil/vf_clay_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/vf_clay_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'vf_clay_s_l'//trim(c1)//'_patches', patchmask)

               ENDDO

               ! topography
               CALL system('mkdir -p ' // trim(dir_landdata_out) // '/topography')
               file_in  = trim(dir_landdata_in)  // '/topography/topography_patches.nc'
               file_out = trim(dir_landdata_out) // '/topography/topography_patches.nc'
               CALL clip_vector (file_in, file_out, iblk, jblk, 'topography_patches', patchmask)

            ENDIF
         ENDDO
      ENDDO

#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, nelm_blk, gblock%nxblk*gblock%nyblk, &
         MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN
         file_out = trim(dir_landdata_out) // '/mesh/mesh.nc'
         CALL ncio_write_serial (file_out, 'nelm_blk', nelm_blk, 'xblk', 'yblk')
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN
         write(*,*) 'Clip surface data from large region successful.'
      ENDIF

      IF (allocated(nelm_blk )) deallocate(nelm_blk )
      IF (allocated(IOproc   )) deallocate(IOproc   )
      IF (allocated(elmindex )) deallocate(elmindex )
      IF (allocated(elmnpxl  )) deallocate(elmnpxl  )
      IF (allocated(elmpixels)) deallocate(elmpixels)
      IF (allocated(elmmask  )) deallocate(elmmask  )
      IF (allocated(patchmask)) deallocate(patchmask)
#ifdef CATCHMENT
      IF (allocated(hrumask  )) deallocate(hrumask  )
#endif
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      IF (allocated(pftmask  )) deallocate(pftmask  )
#endif

   END SUBROUTINE srfdata_region_clip

   ! ----- pixelset clip -----
   SUBROUTINE clip_pixelset (dir_landdata_in, psetname, iblk, jblk, elmmask, elmindex, psetmask)

   USE MOD_NetCDFSerial
   USE MOD_Block

   IMPLICIT NONE

   character(len=*), intent(in) :: dir_landdata_in
   character(len=*), intent(in) :: psetname
   integer, intent(in) :: iblk, jblk
   logical, intent(in) :: elmmask (:)

   integer*8, intent(in) :: elmindex(:)

   logical, allocatable, intent(out) :: psetmask (:)

   ! Local Variables
   character(len=256) :: filename, fileblock
   logical :: fexists
   integer :: nset, ie, iset
   integer*8, allocatable :: eindex_p(:)

      filename = trim(dir_landdata_in) // '/' // trim(psetname) // '/' // trim(psetname) // '.nc'
      CALL get_filename_block (filename, iblk, jblk, fileblock)

      inquire (file=trim(fileblock), exist=fexists)
      IF (fexists) THEN
         CALL ncio_inquire_length (fileblock, 'eindex', nset)

         IF (nset > 0) THEN

            CALL ncio_read_serial (fileblock, 'eindex', eindex_p)

            allocate (psetmask (nset))
            psetmask(:) = .false.

            ie = 1
            DO iset = 1, nset
               DO WHILE ((eindex_p(iset) /= elmindex(ie)))
                  ie = ie + 1
               ENDDO
               psetmask(iset) = elmmask(ie)
            ENDDO

            deallocate(eindex_p)
         ENDIF
      ENDIF

   END SUBROUTINE clip_pixelset

   ! ----- vector clip -----
   SUBROUTINE clip_vector (file_in, file_out, iblk, jblk, varname, vecmask)

   USE netcdf
   USE MOD_NetCDFSerial
   USE MOD_Block
   IMPLICIT NONE

   character(len=*), intent(in) :: file_in, file_out
   character(len=*), intent(in) :: varname
   integer, intent(in) :: iblk, jblk
   logical, intent(in) :: vecmask(:)

   ! Local Variables
   character(len=256) :: input, output
   integer :: ncidin,  ncidout
   integer :: varidin, varidout
   integer :: vlen_in, vlen_out, xtype, ndims, id, i1, i2
   logical :: fexists, dim_exist
   integer, allocatable :: dimids (:), dimlens (:)
   character(len=256), allocatable :: dimnames (:)

   integer, allocatable :: data_i4_in1 (:)
   integer, allocatable :: data_i4_in2 (:,:)
   integer, allocatable :: data_i4_in3 (:,:,:)

   integer, allocatable :: data_i4_out1 (:)
   integer, allocatable :: data_i4_out2 (:,:)
   integer, allocatable :: data_i4_out3 (:,:,:)

   integer*8, allocatable :: data_i8_in1  (:)
   integer*8, allocatable :: data_i8_out1 (:)

   real(r8), allocatable :: data_r8_in1 (:)
   real(r8), allocatable :: data_r8_in2 (:,:)
   real(r8), allocatable :: data_r8_in3 (:,:,:)

   real(r8), allocatable :: data_r8_out1 (:)
   real(r8), allocatable :: data_r8_out2 (:,:)
   real(r8), allocatable :: data_r8_out3 (:,:,:)

      CALL get_filename_block (file_in,  iblk, jblk, input )
      CALL get_filename_block (file_out, iblk, jblk, output)

      write(*,*) 'Copy and clip ', trim(varname), ' from ', trim(input), ' to ', trim(output)

      IF (.not. ncio_var_exist(input, varname)) THEN
         write(*,*) 'Variable ', trim(varname), ' in ', trim(input), ' not found.'
      ENDIF

      vlen_in  = size (vecmask)
      vlen_out = count(vecmask)

      CALL nccheck( nf90_open (trim(input), NF90_NOWRITE, ncidin) )
      CALL nccheck( nf90_inq_varid (ncidin, trim(varname), varidin) )
      CALL nccheck( nf90_inquire_variable (ncidin, varidin, xtype = xtype, ndims = ndims) )
      allocate (dimids   (ndims))
      allocate (dimnames (ndims))
      allocate (dimlens  (ndims))
      CALL nccheck( nf90_inquire_variable (ncidin, varidin, dimids = dimids) )
      DO id = 1, ndims
         CALL nccheck( nf90_inquire_dimension (ncidin, dimids(id), dimnames(id), dimlens(id)) )
      ENDDO
      CALL nccheck( nf90_close (ncidin) )

      inquire (file=trim(output), exist=fexists)
      IF (.not. fexists) THEN
         CALL ncio_create_file (output)
      ENDIF

      DO id = 1, ndims
         CALL nccheck( nf90_open (trim(output), NF90_NOWRITE, ncidout) )
         dim_exist = (nf90_inq_dimid (ncidout, trim(dimnames(id)), dimids(id)) == NF90_NOERR)
         CALL nccheck( nf90_close (ncidout) )

         IF (.not. dim_exist) THEN
            IF (id < ndims) THEN
               CALL nccheck( nf90_open (trim(input), NF90_NOWRITE, ncidin) )
               CALL nccheck( nf90_inq_dimid (ncidin, trim(dimnames(id)), dimids(id)) )
               CALL nccheck( nf90_inquire_dimension (ncidin, dimids(id), len = dimlens(id)) )
               CALL nccheck( nf90_close (ncidin) )
            ELSE
               dimlens(ndims) = vlen_out
            ENDIF

            CALL nccheck( nf90_open (trim(output), NF90_WRITE, ncidout) )
            CALL nccheck( nf90_redef (ncidout) )
            CALL nccheck( nf90_def_dim (ncidout, trim(dimnames(id)), dimlens(id), dimids(id)) )
            CALL nccheck( nf90_enddef (ncidout) )
            CALL nccheck( nf90_close (ncidout) )
         ENDIF
      ENDDO

      CALL nccheck( nf90_open (trim(output), NF90_WRITE, ncidout) )
      DO id = 1, ndims
         CALL nccheck( nf90_inq_dimid (ncidout, dimnames(id), dimids(id)) )
      ENDDO

      IF (nf90_inq_varid(ncidout,trim(varname),varidout) /= nf90_noerr) THEN
         CALL nccheck( nf90_redef (ncidout) )
         CALL nccheck( nf90_def_var (ncidout, trim(varname), xtype, dimids, varidout, deflate_level = 1) )
         CALL nccheck( nf90_enddef (ncidout) )
      ENDIF

      CALL nccheck( nf90_open (trim(input),  NF90_NOWRITE, ncidin)  )
      CALL nccheck( nf90_inq_varid (ncidin,  trim(varname), varidin ) )

      SELECTCASE (ndims)
      CASE (1)
         IF (xtype == NF90_INT) THEN
            allocate (data_i4_in1 (vlen_in))
            CALL nccheck( nf90_get_var (ncidin,  varidin , data_i4_in1) )

            allocate (data_i4_out1 (vlen_out))
            data_i4_out1 = pack(data_i4_in1, vecmask)
            CALL nccheck( nf90_put_var (ncidout, varidout, data_i4_out1) )

            deallocate (data_i4_in1 )
            deallocate (data_i4_out1)
         ELSEIF (xtype == NF90_INT64) THEN
            allocate (data_i8_in1 (vlen_in))
            CALL nccheck( nf90_get_var (ncidin,  varidin , data_i8_in1) )

            allocate (data_i8_out1 (vlen_out))
            data_i8_out1 = pack(data_i8_in1, vecmask)
            CALL nccheck( nf90_put_var (ncidout, varidout, data_i8_out1) )

            deallocate (data_i8_in1 )
            deallocate (data_i8_out1)
         ELSEIF (xtype == NF90_DOUBLE) THEN
            allocate (data_r8_in1 (vlen_in))
            CALL nccheck( nf90_get_var (ncidin,  varidin , data_r8_in1) )

            allocate (data_r8_out1 (vlen_out))
            data_r8_out1 = pack(data_r8_in1, vecmask)
            CALL nccheck( nf90_put_var (ncidout, varidout, data_r8_out1) )

            deallocate (data_r8_in1 )
            deallocate (data_r8_out1)
         ENDIF

      CASE (2)
         IF (xtype == NF90_INT) THEN
            allocate (data_i4_in2 (dimlens(1),vlen_in))
            CALL nccheck( nf90_get_var (ncidin,  varidin , data_i4_in2) )

            allocate (data_i4_out2 (dimlens(1),vlen_out))
            DO i1 = 1, dimlens(1)
               data_i4_out2(i1,:) = pack(data_i4_in2(i1,:), vecmask)
            ENDDO
            CALL nccheck( nf90_put_var (ncidout, varidout, data_i4_out2) )

            deallocate (data_i4_in2 )
            deallocate (data_i4_out2)
         ELSEIF (xtype == NF90_DOUBLE) THEN
            allocate (data_r8_in2 (dimlens(1),vlen_in))
            CALL nccheck( nf90_get_var (ncidin,  varidin , data_r8_in2) )

            allocate (data_r8_out2 (dimlens(1),vlen_out))
            DO i1 = 1, dimlens(1)
               data_r8_out2(i1,:) = pack(data_r8_in2(i1,:), vecmask)
            ENDDO
            CALL nccheck( nf90_put_var (ncidout, varidout, data_r8_out2) )

            deallocate (data_r8_in2 )
            deallocate (data_r8_out2)
         ENDIF

      CASE (3)
         IF (xtype == NF90_INT) THEN
            allocate (data_i4_in3 (dimlens(1),dimlens(2),vlen_in))
            CALL nccheck( nf90_get_var (ncidin,  varidin , data_i4_in3) )

            allocate (data_i4_out3 (dimlens(1),dimlens(2),vlen_out))
            DO i1 = 1, dimlens(1)
               DO i2 = 1, dimlens(2)
                  data_i4_out3(i1,i2,:) = pack(data_i4_in3(i1,i2,:), vecmask)
               ENDDO
            ENDDO
            CALL nccheck( nf90_put_var (ncidout, varidout, data_i4_out3) )

            deallocate (data_i4_in3 )
            deallocate (data_i4_out3)
         ELSEIF (xtype == NF90_DOUBLE) THEN
            allocate (data_r8_in3 (dimlens(1),dimlens(2),vlen_in))
            CALL nccheck( nf90_get_var (ncidin,  varidin , data_r8_in3) )

            allocate (data_r8_out3 (dimlens(1),dimlens(2),vlen_out))
            DO i1 = 1, dimlens(1)
               DO i2 = 1, dimlens(2)
                  data_r8_out3(i1,i2,:) = pack(data_r8_in3(i1,i2,:), vecmask)
               ENDDO
            ENDDO
            CALL nccheck( nf90_put_var (ncidout, varidout, data_r8_out3) )

            deallocate (data_r8_in3 )
            deallocate (data_r8_out3)
         ENDIF

      CASE default
         write(*,*) 'Warning: there is no case for variable: ', trim(varname)
      ENDSELECT

      CALL nccheck( nf90_close (ncidin)  )
      CALL nccheck( nf90_close (ncidout) )

   END SUBROUTINE clip_vector

END MODULE MOD_RegionClip
