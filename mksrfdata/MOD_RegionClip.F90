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
   SUBROUTINE srfdata_region_clip (dir_landdata_in, dir_landdata_out, requires_lake_soilc, requires_spatial_ph)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Utils
#ifdef LULCC
   USE MOD_Const_LC, only: N_land_classification
#endif
   USE MOD_Pixel
   USE MOD_Block
   USE MOD_NetCDFSerial
#if (defined MaterialTransport) && (defined BGC)
   USE MOD_MT_SurfaceData, only: mt_clip_methane_surface_data
#endif

   IMPLICIT NONE

   character(len=*), intent(in) :: dir_landdata_in, dir_landdata_out
   logical, intent(in) :: requires_lake_soilc, requires_spatial_ph

   ! Local Variables
   real(r8) :: edges, edgen, edgew, edgee
   character(len=256) :: file_in, file_out, fileblock
   integer :: iproc, iblk, jblk, ie, ipxl, ilon, ilat, i1
   integer :: nelm_in, nelm_out, lc_year, pdsp
   integer,   allocatable :: nelm_blk(:,:), IOproc(:,:)
   integer*8, allocatable :: elmindex(:)
   integer,   allocatable :: elmnpxl(:), elmpixels(:,:,:), pixels(:,:), datasize(:)
   logical, allocatable :: pixelmask(:)
   logical :: fexists

   logical, allocatable :: elmmask  (:)
   logical, allocatable :: patchmask(:)
#ifdef CATCHMENT
   logical, allocatable :: hrumask  (:)
#endif
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   logical, allocatable :: pftmask  (:)
#endif

   integer :: month, YY, itime, Julian_day, nsl, ilc, ifield
   character(len=3), parameter :: canopy_fields(3) = ['ncd', 'ncw', 'bcw']
   character(len=8), parameter :: twi_fields(6) = &
      [character(len=8) :: 'mean_twi', 'fsatmax', 'fsatdcf', 'alp_twi', 'chi_twi', 'mu_twi']
   integer :: start_year, end_year
   character(len=1) :: c1
   character(len=2) :: c2, cx
   character(len=3) :: c3
   character(len=4) :: cyear, lc_cyear, input_year


      ! Like MKSRFDATA, clip one land-cover epoch per invocation. LULCC
      ! callers invoke mksrfdata for each epoch; LAI years are handled below.
      lc_year = DEF_LC_YEAR
#ifdef LULCC
      IF (lc_year < 2000) lc_year = MAX(1985, (lc_year / 5) * 5)
#endif
      write(lc_cyear,'(i4.4)') lc_year
      file_in = trim(dir_landdata_in)//'/mesh/'//lc_cyear//'/mesh.nc'
      inquire(file=trim(file_in), exist=fexists)
      input_year = lc_cyear
      IF (.not. fexists) THEN
         file_in = trim(dir_landdata_in)//'/mesh/mesh.nc'
         inquire(file=trim(file_in), exist=fexists)
         IF (.not. fexists) CALL clip_error('Missing mesh for requested land-cover year '//lc_cyear)
         ! Select the layout once, never fall back separately for individual fields.
         input_year = ''
         IF (p_is_master) write(*,*) 'RegionClip: reading legacy unyear surface layout for epoch ', lc_cyear
      ENDIF
#if defined URBAN_MODEL || defined HYPERSPECTRAL
      CALL clip_error('RegionClip does not support URBAN_MODEL or HYPERSPECTRAL surface data')
#endif
      IF (DEF_USE_Forcing_Downscaling .or. DEF_USE_Forcing_Downscaling_Simple) &
         CALL clip_error('RegionClip does not support downscaling surface data')

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

         CALL system('mkdir -p ' // trim(dir_landdata_out)//'/mesh/'//lc_cyear)
         file_in  = trim(dir_landdata_in ) // '/mesh/' // trim(input_year) // '/mesh.nc'
         file_out = trim(dir_landdata_out) // '/mesh/' // lc_cyear // '/mesh.nc'
         CALL system('cp ' // trim(file_in) // ' ' // trim(file_out))

         file_in  = trim(dir_landdata_in ) // '/pixel.nc'
         file_out = trim(dir_landdata_out) // '/pixel.nc'
         CALL system('cp ' // trim(file_in) // ' ' // trim(file_out))

         file_in  = trim(dir_landdata_in ) // '/block.nc'
         file_out = trim(dir_landdata_out) // '/block.nc'
         CALL system('cp ' // trim(file_in) // ' ' // trim(file_out))

      ENDIF

      file_in = trim(dir_landdata_in) // '/mesh/' // trim(input_year) // '/mesh.nc'
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

               file_in = trim(dir_landdata_in) // '/mesh/' // trim(input_year) // '/mesh.nc'
               CALL get_filename_block (file_in, iblk, jblk, fileblock)

               CALL ncio_read_serial (fileblock, 'elmindex' , elmindex )
               CALL ncio_read_serial (fileblock, 'elmnpxl'  , elmnpxl  )
               CALL ncio_inquire_varsize(fileblock, 'elmpixels', datasize)
               IF (size(datasize) == 3) THEN
                  CALL ncio_read_serial(fileblock, 'elmpixels', elmpixels)
               ELSEIF (size(datasize) == 2) THEN
                  CALL ncio_read_serial(fileblock, 'elmpixels', pixels)
                  allocate(pixelmask(size(pixels,2)))
                  pixelmask = .false.
               ELSE
                  CALL clip_error('Unsupported mesh pixel rank: '//trim(fileblock))
               ENDIF

               nelm_in = nelm_blk(iblk,jblk)

               IF (allocated(elmmask)) deallocate(elmmask)
               allocate (elmmask (nelm_in))
               elmmask(:) = .false.

               pdsp = 0
               DO ie = 1, nelm_in
                  DO ipxl = 1, elmnpxl(ie)
                     IF (size(datasize) == 3) THEN
                        ilon = elmpixels(1,ipxl,ie)
                        ilat = elmpixels(2,ipxl,ie)
                     ELSE
                        ilon = pixels(1,pdsp+ipxl)
                        ilat = pixels(2,pdsp+ipxl)
                     ENDIF

                     IF (((pixel%lat_s(ilat) < edgen) .and. (pixel%lat_n(ilat) > edges)) &
                        .and. (lon_between_floor(pixel%lon_w(ilon), edgew, edgee) &
                        .or. lon_between_ceil(pixel%lon_e(ilon), edgew, edgee))) THEN

                        elmmask(ie) = .true.

                        EXIT
                     ELSE
                        CYCLE
                     ENDIF
                  ENDDO
                  IF (allocated(pixelmask)) pixelmask(pdsp+1:pdsp+elmnpxl(ie)) = elmmask(ie)
                  pdsp = pdsp + elmnpxl(ie)
               ENDDO

               nelm_out = count(elmmask)
               nelm_blk(iblk,jblk) = nelm_out

               IF (nelm_out == 0) THEN
                  IF (allocated(pixelmask)) deallocate(pixelmask)
                  CYCLE
               ENDIF
               IF (nelm_out > 0) THEN

                  CALL clip_pixelset (dir_landdata_in, trim(input_year), 'landpatch', iblk, jblk, elmmask, elmindex, patchmask)
#ifdef CATCHMENT
                  CALL clip_pixelset (dir_landdata_in, trim(input_year), 'landhru'  , iblk, jblk, elmmask, elmindex, hrumask  )
#endif
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
                  CALL clip_pixelset (dir_landdata_in, trim(input_year), 'landpft'  , iblk, jblk, elmmask, elmindex, pftmask  )
#endif
                  CALL system('mkdir -p ' // trim(dir_landdata_out) // '/mesh/' // lc_cyear)
                  file_in  = trim(dir_landdata_in)  // '/mesh/' // trim(input_year) // '/mesh.nc'
                  file_out = trim(dir_landdata_out) // '/mesh/' // lc_cyear // '/mesh.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'elmindex' , elmmask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'elmnpxl'  , elmmask)
                  IF (size(datasize) == 3) THEN
                     CALL clip_vector(file_in, file_out, iblk, jblk, 'elmpixels', elmmask)
                  ELSE
                     CALL clip_vector(file_in, file_out, iblk, jblk, 'elmpixels', pixelmask)
                     deallocate(pixelmask)
                  ENDIF

                  CALL system('mkdir -p ' // trim(dir_landdata_out) // '/landelm/' // lc_cyear)
                  file_in  = trim(dir_landdata_in)  // '/landelm/' // trim(input_year) // '/landelm.nc'
                  file_out = trim(dir_landdata_out) // '/landelm/' // lc_cyear // '/landelm.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'eindex', elmmask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'ipxstt', elmmask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'ipxend', elmmask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'settyp', elmmask)
                  CALL get_filename_block(file_in, iblk, jblk, fileblock)
                  IF (ncio_var_exist(fileblock, 'pctshared', readflag=.false.)) &
                     CALL clip_vector(file_in, file_out, iblk, jblk, 'pctshared', elmmask)

                  CALL system('mkdir -p ' // trim(dir_landdata_out) // '/landpatch/' // lc_cyear)
                  file_in  = trim(dir_landdata_in)  // '/landpatch/' // trim(input_year) // '/landpatch.nc'
                  file_out = trim(dir_landdata_out) // '/landpatch/' // lc_cyear // '/landpatch.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'eindex', patchmask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'ipxstt', patchmask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'ipxend', patchmask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'settyp', patchmask)
                  CALL get_filename_block(file_in, iblk, jblk, fileblock)
                  IF (ncio_var_exist(fileblock, 'pctshared', readflag=.false.)) &
                     CALL clip_vector(file_in, file_out, iblk, jblk, 'pctshared', patchmask)

                  file_in  = trim(dir_landdata_in)  // '/landpatch/' // trim(input_year) // '/patchfrac_elm.nc'
                  file_out = trim(dir_landdata_out) // '/landpatch/' // lc_cyear // '/patchfrac_elm.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'patchfrac_elm', patchmask)

#ifdef CATCHMENT
                  CALL system('mkdir -p ' // trim(dir_landdata_out) // '/landhru/' // lc_cyear)
                  file_in  = trim(dir_landdata_in)  // '/landhru/' // trim(input_year) // '/landhru.nc'
                  file_out = trim(dir_landdata_out) // '/landhru/' // lc_cyear // '/landhru.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'eindex', hrumask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'ipxstt', hrumask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'ipxend', hrumask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'settyp', hrumask)
                  CALL get_filename_block(file_in, iblk, jblk, fileblock)
                  IF (ncio_var_exist(fileblock, 'pctshared', readflag=.false.)) &
                     CALL clip_vector(file_in, file_out, iblk, jblk, 'pctshared', hrumask)

                  file_in  = trim(dir_landdata_in)  // '/landpatch/' // trim(input_year) // '/patchfrac_hru.nc'
                  file_out = trim(dir_landdata_out) // '/landpatch/' // lc_cyear // '/patchfrac_hru.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'patchfrac_hru', patchmask)
#endif

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
                  CALL system('mkdir -p ' // trim(dir_landdata_out) // '/landpft/' // lc_cyear)
                  file_in  = trim(dir_landdata_in)  // '/landpft/' // trim(input_year) // '/landpft.nc'
                  file_out = trim(dir_landdata_out) // '/landpft/' // lc_cyear // '/landpft.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'eindex', pftmask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'ipxstt', pftmask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'ipxend', pftmask)
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'settyp', pftmask)
                  CALL get_filename_block(file_in, iblk, jblk, fileblock)
                  IF (ncio_var_exist(fileblock, 'pctshared', readflag=.false.)) &
                     CALL clip_vector(file_in, file_out, iblk, jblk, 'pctshared', pftmask)
#endif
               ENDIF

               ! Leaf Area Index
               CALL system('mkdir -p ' // trim(dir_landdata_out) // '/LAI')

#ifdef LULCC
               start_year = DEF_LC_YEAR
               end_year = DEF_LC_YEAR
#else
               IF (DEF_LAI_CHANGE_YEARLY) THEN
                  start_year = min(max(DEF_simulation_time%start_year, DEF_LAI_START_YEAR), DEF_LAI_END_YEAR)
                  end_year = max(min(DEF_simulation_time%end_year, DEF_LAI_END_YEAR), DEF_LAI_START_YEAR)
               ELSE
                  start_year = DEF_LC_YEAR
                  end_year   = DEF_LC_YEAR
               ENDIF

#endif
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

#ifdef LULCC
               IF (lc_year >= 1990) THEN
                  CALL system('mkdir -p '//trim(dir_landdata_out)//'/lulcc/'//lc_cyear)
                  DO ilc = 0, N_land_classification
                     write(c2,'(i2.2)') ilc
                     file_in = trim(dir_landdata_in)//'/lulcc/'//trim(input_year)//'/lccpct_patches_lc'//c2//'.nc'
                     file_out = trim(dir_landdata_out)//'/lulcc/'//lc_cyear//'/lccpct_patches_lc'//c2//'.nc'
                     CALL clip_vector(file_in, file_out, iblk, jblk, 'lccpct_patches', patchmask)
                  ENDDO
               ENDIF
#endif

               ! depth to bedrock
               IF(DEF_USE_BEDROCK)THEN
                  CALL system('mkdir -p ' // trim(dir_landdata_out) // '/dbedrock/' // lc_cyear)
                  file_in  = trim(dir_landdata_in)  // '/dbedrock/' // trim(input_year) // '/dbedrock_patches.nc'
                  file_out = trim(dir_landdata_out) // '/dbedrock/' // lc_cyear // '/dbedrock_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, 'dbedrock_patches', patchmask)
               ENDIF

               ! forest height
               CALL system('mkdir -p ' // trim(dir_landdata_out) // '/htop/' // lc_cyear)
               file_in  = trim(dir_landdata_in)  // '/htop/' // trim(input_year) // '/htop_patches.nc'
               file_out = trim(dir_landdata_out) // '/htop/' // lc_cyear // '/htop_patches.nc'
               CALL clip_vector (file_in, file_out, iblk, jblk, 'htop_patches', patchmask)
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
               file_in  = trim(dir_landdata_in)  // '/htop/' // trim(input_year) // '/htop_pfts.nc'
               file_out = trim(dir_landdata_out) // '/htop/' // lc_cyear // '/htop_pfts.nc'
               CALL clip_vector (file_in, file_out, iblk, jblk, 'htop_pfts', pftmask)
               IF (DEF_RS_CROWN_STRUCTURE) THEN
                  file_in = trim(dir_landdata_in)//'/htop/'//trim(input_year)//'/hbot_pfts.nc'
                  file_out = trim(dir_landdata_out)//'/htop/'//lc_cyear//'/hbot_pfts.nc'
                  CALL clip_vector(file_in, file_out, iblk, jblk, 'hbot_pfts', pftmask)
#ifdef LULC_IGBP_PC
                  file_in = trim(dir_landdata_in)//'/htop/'//trim(input_year)//'/cratio_pfts.nc'
                  file_out = trim(dir_landdata_out)//'/htop/'//lc_cyear//'/cratio_pfts.nc'
                  CALL clip_vector(file_in, file_out, iblk, jblk, 'cratio_pfts', pftmask)
#endif
               ENDIF
#endif

               IF (DEF_Interception_scheme == 8) THEN
                  CALL system('mkdir -p '//trim(dir_landdata_out)//'/cstructure/'//lc_cyear)
                  DO ifield = 1, size(canopy_fields)
                     file_in = trim(dir_landdata_in)//'/cstructure/'//trim(input_year)//'/'// &
                        canopy_fields(ifield)//'_patches.nc'
                     file_out = trim(dir_landdata_out)//'/cstructure/'//lc_cyear//'/'// &
                        canopy_fields(ifield)//'_patches.nc'
                     CALL clip_vector(file_in, file_out, iblk, jblk, canopy_fields(ifield)//'_patches', patchmask)
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
                     file_in = trim(dir_landdata_in)//'/cstructure/'//trim(input_year)//'/'// &
                        canopy_fields(ifield)//'_pfts.nc'
                     file_out = trim(dir_landdata_out)//'/cstructure/'//lc_cyear//'/'// &
                        canopy_fields(ifield)//'_pfts.nc'
                     CALL clip_vector(file_in, file_out, iblk, jblk, canopy_fields(ifield)//'_pfts', pftmask)
#endif
                  ENDDO
               ENDIF

               ! lake depth
               CALL system('mkdir -p ' // trim(dir_landdata_out) // '/lakedepth/' // lc_cyear)
               file_in  = trim(dir_landdata_in)  // '/lakedepth/' // trim(input_year) // '/lakedepth_patches.nc'
               file_out = trim(dir_landdata_out) // '/lakedepth/' // lc_cyear // '/lakedepth_patches.nc'
               CALL clip_vector (file_in, file_out, iblk, jblk, 'lakedepth_patches', patchmask)

               ! plant function type percentage
               CALL system('mkdir -p ' // trim(dir_landdata_out) // '/pctpft/' // lc_cyear)
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
               file_in  = trim(dir_landdata_in)  // '/pctpft/' // trim(input_year) // '/pct_pfts.nc'
               file_out = trim(dir_landdata_out) // '/pctpft/' // lc_cyear // '/pct_pfts.nc'
               CALL clip_vector (file_in, file_out, iblk, jblk, 'pct_pfts', pftmask)
#ifdef CROP
               file_in  = trim(dir_landdata_in)  // '/pctpft/' // trim(input_year) // '/pct_crops.nc'
               file_out = trim(dir_landdata_out) // '/pctpft/' // lc_cyear // '/pct_crops.nc'
               CALL clip_vector (file_in, file_out, iblk, jblk, 'pct_crops', patchmask)
#endif
#endif

               ! soil
               CALL system('mkdir -p ' // trim(dir_landdata_out) // '/soil/' // lc_cyear)

#if (defined MaterialTransport) && (defined BGC)
               CALL mt_clip_methane_surface_data(dir_landdata_in, dir_landdata_out, input_year, lc_cyear, &
                  iblk, jblk, patchmask, requires_lake_soilc, requires_spatial_ph, clip_vector, clip_error)
#endif

               file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/soil_s_v_alb_patches.nc'
               file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/soil_s_v_alb_patches.nc'
               CALL clip_vector (file_in, file_out, iblk, jblk, 'soil_s_v_alb', patchmask)

               file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/soil_d_v_alb_patches.nc'
               file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/soil_d_v_alb_patches.nc'
               CALL clip_vector (file_in, file_out, iblk, jblk, 'soil_d_v_alb', patchmask)

               file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/soil_s_n_alb_patches.nc'
               file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/soil_s_n_alb_patches.nc'
               CALL clip_vector (file_in, file_out, iblk, jblk, 'soil_s_n_alb', patchmask)

               file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/soil_d_n_alb_patches.nc'
               file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/soil_d_n_alb_patches.nc'
               CALL clip_vector (file_in, file_out, iblk, jblk, 'soil_d_n_alb', patchmask)

               DO nsl = 1, 8

                  write(c1,'(i1)') nsl

                  ! (1) volumetric fraction of quartz within mineral soil
                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/vf_quartz_mineral_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/vf_quartz_mineral_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'vf_quartz_mineral_s_l'//trim(c1)//'_patches', patchmask)

                  ! (2) volumetric fraction of gravels
                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/vf_gravels_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/vf_gravels_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'vf_gravels_s_l'//trim(c1)//'_patches', patchmask)

                  ! (3) volumetric fraction of sand
                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/vf_sand_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/vf_sand_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'vf_sand_s_l'//trim(c1)//'_patches', patchmask)

                  ! (4) volumetric fraction of organic matter with the parameter alpha and beta in the Balland V. and P. A. Arp (2005) model
                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/vf_om_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/vf_om_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'vf_om_s_l'//trim(c1)//'_patches', patchmask)

                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/BA_alpha_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/BA_alpha_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'BA_alpha_l'//trim(c1)//'_patches', patchmask)

                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/BA_beta_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/BA_beta_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'BA_beta_l'//trim(c1)//'_patches', patchmask)

                  ! (5) gravimetric fraction of gravels
                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/wf_gravels_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/wf_gravels_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'wf_gravels_s_l'//trim(c1)//'_patches', patchmask)

                  ! (6) gravimetric fraction of sand
                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/wf_sand_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/wf_sand_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'wf_sand_s_l'//trim(c1)//'_patches', patchmask)

#ifdef vanGenuchten_Mualem_SOIL_MODEL
                  ! (7) VGM's pore-connectivity parameter (L)
                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/L_vgm_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/L_vgm_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'L_vgm_l'//trim(c1)//'_patches', patchmask)

                  ! (8) VGM's residual water content (theta_r) [cm3/cm3]
                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/theta_r_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/theta_r_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'theta_r_l'//trim(c1)//'_patches', patchmask)

                  ! (9) VGM's parameter corresponding approximately to the inverse of the air-entry value (alpha)
                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/alpha_vgm_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/alpha_vgm_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'alpha_vgm_l'//trim(c1)//'_patches', patchmask)

                  ! (10) VGM's shape parameter (n)
                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/n_vgm_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/n_vgm_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'n_vgm_l'//trim(c1)//'_patches', patchmask)
#endif
                  ! (11) saturated water content [cm3/cm3]
                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/theta_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/theta_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'theta_s_l'//trim(c1)//'_patches', patchmask)

                  ! (12) matric potential at saturation (psi_s) [cm]
                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/psi_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/psi_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'psi_s_l'//trim(c1)//'_patches', patchmask)

                  ! (13) pore size distribution index [dimensionless]
                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/lambda_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/lambda_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'lambda_l'//trim(c1)//'_patches', patchmask)

                  ! (14) saturated hydraulic conductivity [cm/day]
                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/k_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/k_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'k_s_l'//trim(c1)//'_patches', patchmask)

                  ! (15) heat capacity of soil solids [J/(m3 K)]
                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/csol_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/csol_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'csol_l'//trim(c1)//'_patches', patchmask)

                  ! (16) thermal conductivity of unfrozen saturated soil [W/m-K]
                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/tksatu_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/tksatu_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'tksatu_l'//trim(c1)//'_patches', patchmask)

                  ! (17) thermal conductivity of frozen saturated soil [W/m-K]
                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/tksatf_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/tksatf_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'tksatf_l'//trim(c1)//'_patches', patchmask)

                  ! (18) thermal conductivity for dry soil [W/(m-K)]
                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/tkdry_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/tkdry_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'tkdry_l'//trim(c1)//'_patches', patchmask)

                  ! (19) thermal conductivity of soil solids [W/m-K]
                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/k_solids_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/k_solids_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'k_solids_l'//trim(c1)//'_patches', patchmask)

                  ! (20) OM_density [kg/m3]
                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/OM_density_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/OM_density_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'OM_density_s_l'//trim(c1)//'_patches', patchmask)

                  ! (21) bulk density of soil (GRAVELS + OM + Mineral Soils)
                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/BD_all_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/BD_all_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'BD_all_s_l'//trim(c1)//'_patches', patchmask)

                  file_in = trim(dir_landdata_in)//'/soil/'//trim(input_year)//'/wf_clay_s_l'//c1//'_patches.nc'
                  file_out = trim(dir_landdata_out)//'/soil/'//lc_cyear//'/wf_clay_s_l'//c1//'_patches.nc'
                  CALL clip_vector(file_in, file_out, iblk, jblk, 'wf_clay_s_l'//c1//'_patches', patchmask)
                  file_in = trim(dir_landdata_in)//'/soil/'//trim(input_year)//'/wf_om_s_l'//c1//'_patches.nc'
                  file_out = trim(dir_landdata_out)//'/soil/'//lc_cyear//'/wf_om_s_l'//c1//'_patches.nc'
                  CALL clip_vector(file_in, file_out, iblk, jblk, 'wf_om_s_l'//c1//'_patches', patchmask)
                  ! (22) volumetric fraction of clay
                  file_in  = trim(dir_landdata_in)  // '/soil/' // trim(input_year) // '/vf_clay_s_l'//trim(c1)//'_patches.nc'
                  file_out = trim(dir_landdata_out) // '/soil/' // lc_cyear // '/vf_clay_s_l'//trim(c1)//'_patches.nc'
                  CALL clip_vector (file_in, file_out, iblk, jblk, &
                     'vf_clay_s_l'//trim(c1)//'_patches', patchmask)

               ENDDO

               ! topography
               CALL system('mkdir -p ' // trim(dir_landdata_out) // '/topography/' // lc_cyear)
               file_in = trim(dir_landdata_in)//'/topography/'//trim(input_year)//'/elevation_patches.nc'
               file_out = trim(dir_landdata_out)//'/topography/'//lc_cyear//'/elevation_patches.nc'
               CALL clip_vector(file_in, file_out, iblk, jblk, 'elevation_patches', patchmask)
               file_in = trim(dir_landdata_in)//'/topography/'//trim(input_year)//'/elvstd_patches.nc'
               file_out = trim(dir_landdata_out)//'/topography/'//lc_cyear//'/elvstd_patches.nc'
               CALL clip_vector(file_in, file_out, iblk, jblk, 'elvstd_patches', patchmask)
               file_in = trim(dir_landdata_in)//'/topography/'//trim(input_year)//'/sloperatio_patches.nc'
               file_out = trim(dir_landdata_out)//'/topography/'//lc_cyear//'/sloperatio_patches.nc'
               CALL clip_vector(file_in, file_out, iblk, jblk, 'sloperatio_patches', patchmask)
               IF (DEF_Runoff_SCHEME == 0 .and. DEF_TOPMOD_method /= 0) THEN
                  DO ifield = 1, size(twi_fields)
                     file_in = trim(dir_landdata_in)//'/topography/'//trim(input_year)//'/'// &
                        trim(twi_fields(ifield))//'_patches.nc'
                     file_out = trim(dir_landdata_out)//'/topography/'//lc_cyear//'/'// &
                        trim(twi_fields(ifield))//'_patches.nc'
                     CALL clip_vector(file_in, file_out, iblk, jblk, trim(twi_fields(ifield))//'_patches', patchmask)
                  ENDDO
               ENDIF
               file_in = trim(dir_landdata_in)//'/soil/'//trim(input_year)//'/soiltexture_patches.nc'
               file_out = trim(dir_landdata_out)//'/soil/'//lc_cyear//'/soiltexture_patches.nc'
               CALL clip_vector(file_in, file_out, iblk, jblk, 'soiltext_patches', patchmask)


            ENDIF
         ENDDO
      ENDDO

#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, nelm_blk, gblock%nxblk*gblock%nyblk, &
         MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN
         file_out = trim(dir_landdata_out) // '/mesh/' // lc_cyear // '/mesh.nc'
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

   SUBROUTINE clip_error(message)
      USE MOD_SPMD_Task
      character(len=*), intent(in) :: message
      write(*,*) 'RegionClip: ', trim(message)
#ifdef USEMPI
      CALL mpi_abort(p_comm_glb, 1, p_err)
#else
      CALL EXIT(1)
#endif
   END SUBROUTINE clip_error

   ! ----- pixelset clip -----
   SUBROUTINE clip_pixelset (dir_landdata_in, lc_cyear, psetname, iblk, jblk, elmmask, elmindex, psetmask)

   USE MOD_NetCDFSerial
   USE MOD_Block
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_Const_LC, only: patchtypes_igbp
#endif

   IMPLICIT NONE

   character(len=*), intent(in) :: dir_landdata_in
   character(len=*), intent(in) :: lc_cyear, psetname
   integer, intent(in) :: iblk, jblk
   logical, intent(in) :: elmmask (:)

   integer*8, intent(in) :: elmindex(:)

   logical, allocatable, intent(out) :: psetmask (:)

   ! Local Variables
   character(len=256) :: filename, fileblock
   logical :: fexists
   integer :: nset, ie, iset
   integer*8, allocatable :: eindex_p(:)
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   integer, allocatable :: patchtypes_source(:)
#endif

      filename = trim(dir_landdata_in) // '/' // trim(psetname) // '/' // lc_cyear // '/' // trim(psetname) // '.nc'
      CALL get_filename_block (filename, iblk, jblk, fileblock)

      inquire (file=trim(fileblock), exist=fexists)
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      IF (.not. fexists .and. psetname == 'landpft') THEN
         ! Vector files have no block-count header: the producer omits zero-PFT
         ! blocks. Prove zero from ALL source patches, not just the clipped subset.
         ! Use the constant table: clip runs before Init_LC_Const initializes patchtypes.
         filename = trim(dir_landdata_in)//'/landpatch/'//lc_cyear//'/landpatch.nc'
         CALL get_filename_block(filename, iblk, jblk, fileblock)
         CALL ncio_read_serial(fileblock, 'settyp', patchtypes_source)
         IF (any(patchtypes_source < 1) .or. any(patchtypes_source > size(patchtypes_igbp))) &
            CALL clip_error('Invalid source patch types while checking absent PFT block: '//trim(fileblock))
         IF (size(patchtypes_source) > 0) THEN
            IF (all(patchtypes_igbp(patchtypes_source) /= 0)) THEN
               allocate(psetmask(0))
               RETURN
            ENDIF
         ENDIF
         CALL clip_error('Missing pixelset landpft for PFT-bearing source block: '//trim(fileblock))
      ENDIF
#endif
      IF (.not. fexists) CALL clip_error('Missing pixelset: '//trim(fileblock))
      IF (fexists) THEN
         CALL ncio_inquire_length (fileblock, 'eindex', nset)

         allocate(psetmask(nset))
         IF (nset > 0) THEN

            CALL ncio_read_serial (fileblock, 'eindex', eindex_p)

            psetmask(:) = .false.

            ie = 1
            DO iset = 1, nset
               DO WHILE ((eindex_p(iset) /= elmindex(ie)))
                  ie = ie + 1
                  IF (ie > size(elmindex)) CALL clip_error('Pixelset element absent from mesh: '//trim(fileblock))
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
   integer :: vlen_in, vlen_out, xtype, ndims, id, i1, i2, natts, iatt, dimlen_out, dimlen_expected
   character(len=256) :: attname
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
   integer*8, allocatable :: data_i8_in2  (:,:)
   integer*8, allocatable :: data_i8_in3  (:,:,:)
   integer*8, allocatable :: data_i8_out1 (:)
   integer*8, allocatable :: data_i8_out2 (:,:)
   integer*8, allocatable :: data_i8_out3 (:,:,:)

   real(r8), allocatable :: data_r8_in1 (:)
   real(r8), allocatable :: data_r8_in2 (:,:)
   real(r8), allocatable :: data_r8_in3 (:,:,:)

   real(r8), allocatable :: data_r8_out1 (:)
   real(r8), allocatable :: data_r8_out2 (:,:)
   real(r8), allocatable :: data_r8_out3 (:,:,:)

      ! Zero selected entries have no vector block file in the producer contract.
      CALL get_filename_block (file_out, iblk, jblk, output)
      IF (count(vecmask) == 0) THEN
         inquire (file=trim(output), exist=fexists)
         IF (fexists) CALL clip_error('Output dimension mismatch: '//trim(output)//' [empty selection]')
         RETURN
      ENDIF

      CALL get_filename_block (file_in,  iblk, jblk, input )

      write(*,*) 'Copy and clip ', trim(varname), ' from ', trim(input), ' to ', trim(output)

      inquire (file=trim(input), exist=fexists)
      IF (.not. fexists) CALL clip_error('Missing required block file: '//trim(input)//' [field '//trim(varname)//']')

      vlen_in  = size (vecmask)
      vlen_out = count(vecmask)

      CALL nccheck( nf90_open (trim(input), NF90_NOWRITE, ncidin), &
         'RegionClip required field '//trim(varname)//' in '//trim(input) )
      IF (nf90_inq_varid(ncidin, trim(varname), varidin) /= NF90_NOERR) THEN
         CALL nccheck(nf90_close(ncidin))
         CALL clip_error('Missing required field: '//trim(varname)//' in '//trim(input))
      ENDIF
      CALL nccheck( nf90_inquire_variable (ncidin, varidin, xtype = xtype, ndims = ndims) )
      IF (ndims < 1 .or. ndims > 3 .or. &
          (xtype /= NF90_INT .and. xtype /= NF90_INT64 .and. xtype /= NF90_DOUBLE)) &
         CALL clip_error('Unsupported vector type/rank: '//trim(input))
      allocate (dimids   (ndims))
      allocate (dimnames (ndims))
      allocate (dimlens  (ndims))
      CALL nccheck( nf90_inquire_variable (ncidin, varidin, dimids = dimids) )
      DO id = 1, ndims
         CALL nccheck( nf90_inquire_dimension (ncidin, dimids(id), dimnames(id), dimlens(id)) )
      ENDDO
      IF (dimlens(ndims) /= vlen_in) CALL clip_error('Vector/mask length mismatch: '//trim(input))
      CALL nccheck( nf90_close (ncidin) )

      inquire (file=trim(output), exist=fexists)
      IF (.not. fexists) THEN
         CALL ncio_create_file (output)
      ELSE
         ! Validate every existing dimension before changing a reused block file.
         CALL nccheck( nf90_open (trim(output), NF90_NOWRITE, ncidout) )
         DO id = 1, ndims
            IF (nf90_inq_dimid(ncidout, trim(dimnames(id)), dimids(id)) /= NF90_NOERR) CYCLE
            CALL nccheck(nf90_inquire_dimension(ncidout, dimids(id), len=dimlen_out))
            dimlen_expected = dimlens(id)
            IF (id == ndims) dimlen_expected = vlen_out
            IF (dimlen_out /= dimlen_expected) THEN
               CALL nccheck(nf90_close(ncidout))
               CALL clip_error('Output dimension mismatch: '//trim(output)//' ['//trim(dimnames(id))//']')
            ENDIF
         ENDDO
         CALL nccheck( nf90_close (ncidout) )
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
         CALL nccheck(nf90_open(trim(input), NF90_NOWRITE, ncidin))
         CALL nccheck(nf90_inq_varid(ncidin, trim(varname), varidin))
         CALL nccheck(nf90_inquire_variable(ncidin, varidin, natts=natts))
         DO iatt = 1, natts
            CALL nccheck(nf90_inq_attname(ncidin, varidin, iatt, attname))
            CALL nccheck(nf90_copy_att(ncidin, varidin, trim(attname), ncidout, varidout))
         ENDDO
         CALL nccheck(nf90_close(ncidin))
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
         ELSEIF (xtype == NF90_INT64) THEN
            allocate (data_i8_in2 (dimlens(1),vlen_in))
            CALL nccheck( nf90_get_var (ncidin, varidin, data_i8_in2) )

            allocate (data_i8_out2 (dimlens(1),vlen_out))
            DO i1 = 1, dimlens(1)
               data_i8_out2(i1,:) = pack(data_i8_in2(i1,:), vecmask)
            ENDDO
            CALL nccheck( nf90_put_var (ncidout, varidout, data_i8_out2) )

            deallocate (data_i8_in2 )
            deallocate (data_i8_out2)
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
         ELSEIF (xtype == NF90_INT64) THEN
            allocate (data_i8_in3 (dimlens(1),dimlens(2),vlen_in))
            CALL nccheck( nf90_get_var (ncidin, varidin, data_i8_in3) )

            allocate (data_i8_out3 (dimlens(1),dimlens(2),vlen_out))
            DO i1 = 1, dimlens(1)
               DO i2 = 1, dimlens(2)
                  data_i8_out3(i1,i2,:) = pack(data_i8_in3(i1,i2,:), vecmask)
               ENDDO
            ENDDO
            CALL nccheck( nf90_put_var (ncidout, varidout, data_i8_out3) )

            deallocate (data_i8_in3 )
            deallocate (data_i8_out3)
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
