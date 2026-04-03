#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_SM
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Data assimilation of surface soil moisture and temperature
!
! AUTHOR:
!   Lu Li, 12/2024: Initial version, based on SMAP L1C TB data
!   Zhilong Fan, Lu Li, 03/2024: Debug and clean codes
!   Lu Li, 07/2025: reframe codes
!   Lu Li, 10/2025: support SMAP/FY3D/SYNOP data
!-----------------------------------------------------------------------------
   USE MOD_DataType
   USE MOD_SpatialMapping
   USE MOD_DA_RTM
   USE MOD_DA_EnKF
   USE MOD_DA_Vars_TimeVariables
   USE MOD_Vars_Global, only: pi, nl_soil, maxsnl
   USE MOD_LandPatch
   USE MOD_Block
   USE MOD_Namelist
   USE MOD_DA_Vars_TimeVariables
   USE MOD_Pixelset
   USE MOD_Pixel
   USE MOD_Mesh
   IMPLICIT NONE
   SAVE

! public functions
   PUBLIC :: init_DA_SM
   PUBLIC :: run_DA_SM
   PUBLIC :: end_DA_SM

   PRIVATE

! local variables
!#############################################################################
! Universal variables
!#############################################################################
   ! time (UTC) at current step
   integer :: month, mday, hour                              ! month, day, hour of current step
   character(len=256) :: yearstr, monthstr, daystr, hourstr  ! string of year, month, day, hour
   integer :: idate_b, idate_e                               ! begin & end seconds since begin of current day (UTC)

   integer :: num_obs_p
   ! mask water
   logical, allocatable :: filter(:)                         ! to mask the water

   ! logical
   logical, allocatable :: has_DA                            ! whether has data assimilation

   ! setting of SMAP satellite
   real(r8), parameter :: smap_theta = 40.0*pi/180           ! incidence angle (rad)
   real(r8), parameter :: smap_fghz = 1.4                    ! frequency (GHz), L-band

   ! setting of FY satellite
   real(r8), parameter :: fy3d_theta = 40.0*pi/180           ! incidence angle (rad)
   real(r8), parameter :: fy3d_fghz = 1.4                    ! frequency (GHz), L-band

   ! parameters of LETKF
   real(r8), allocatable :: obs_err(:)                       ! observation error
   real(r8), parameter   :: dres = 0.4                       ! search localization radius (deg)
   real(r8), parameter   :: loc_r = 40.0                     ! localization radius (km)
   real(r8), parameter   :: infl = 1.2                       ! inflation factor

   ! data assimilation outputs
   real(r8), allocatable :: trans(:,:)                       ! transformation matrix on each patch
   real(r8), allocatable :: wliq_soi_ens(:,:)                ! soil liquid water content
   real(r8), allocatable :: wliq_soi_ens_da(:,:)             ! soil liquid water content after data assimilation
   real(r8), allocatable :: t_soi_ens(:,:)                   ! soil temperature
   real(r8), allocatable :: t_soi_ens_da(:,:)                ! soil temperature after data assimilation
   real(r8) :: eff_porsl                                     ! effective porosity of each soil layer
   integer :: end_idx(3)

   real(r8), allocatable :: pred_obs_p_ens(:,:)              ! predicted observations on each patch
   real(r8), allocatable :: obs_p(:)                         ! observations on each patch
   real(r8), allocatable :: obs_err_p(:)                     ! observation errors on each patch

   real(r8), allocatable :: smap_err(:)                      ! observation errors for SMAP TB
   real(r8), allocatable :: fy3d_err(:)                      ! observation errors for FY3D TB
   real(r8), allocatable :: synop_tref_err(:)                ! observation errors for SYNOP 2m temperature
   real(r8), allocatable :: synop_qref_err(:)                ! observation errors for SYNOP 2m humidity

!#############################################################################
! For SMAP L1C brightness temperature
!#############################################################################
   ! file path
   character(len=256) :: file_smap                           ! SMAP observation file path
   character(len=256) :: file_grid_smap                      ! SMAP world grid file path

   ! grid
   type(grid_type) :: grid_smap                              ! SMAP world grid
   type(spatial_mapping_type) :: mg2p_smap                   ! mapping between SMAP world grid to patch

   ! time variables used to determine whether has obs
   real(r8), allocatable :: smap_time(:)                     ! seconds of all obs since begin of current day (UTC)
   real(r8), allocatable :: dt_b_smap(:)                     ! delta time between SMAP and begin seconds of current day
   real(r8), allocatable :: dt_e_smap(:)                     ! delta time between SMAP and end seconds of current day

   ! logical variables
   logical :: has_smap_file                                  ! whether has file of SMAP at target hour
   logical :: has_smap_obs                                   ! whether has SMAP obs at current step

   ! observations (dimensions changes with time)
   integer :: num_smap_obs                                   ! number of all obs in current file
   integer :: num_smap_obs_domain                            ! number of all obs in simulation domain at current file
   real(r8), allocatable :: smap_lat(:)                      ! latitude of all obs
   real(r8), allocatable :: smap_lon(:)                      ! longitude of all obs
   real(r8), allocatable :: smap_tb_h(:)                     ! H- polarized brightness temperature of all obs ([K])
   real(r8), allocatable :: smap_tb_v(:)                     ! V- polarized brightness temperature of all obs ([K])
   integer,  allocatable :: smap_ii(:)                       ! i-th lat grid in world map of all obs
   integer,  allocatable :: smap_jj(:)                       ! j-th lon grid in world map of all obs

   ! history mean value of smap and colm
   real(r8), allocatable :: pred_smap_tb_h_mean(:)           ! history mean value of predicted H- polarized temp
   real(r8), allocatable :: pred_smap_tb_v_mean(:)           ! history mean value of predicted V- polarized temp
   real(r8), allocatable :: smap_tb_h_mean(:)                ! history mean value of SMAP H- polarized temp
   real(r8), allocatable :: smap_tb_v_mean(:)                ! history mean value of SMAP V- polarized temp

   ! ensemble predicted observations (at patch)
   real(r8), allocatable :: pred_smap_tb_h_pset_ens(:, :)    ! predicted H- polarized temp on patch
   real(r8), allocatable :: pred_smap_tb_v_pset_ens(:, :)    ! predicted V- polarized temp on patch
   real(r8), allocatable :: area_pset(:)                     ! area of each patch across world grid

   ! ensemble predicted observations (at world grid)
   type(block_data_real8_3d) :: pred_smap_tb_h_wgrid_ens     ! predicted H- polarized temp on world grid
   type(block_data_real8_3d) :: pred_smap_tb_v_wgrid_ens     ! predicted V- polarized temp on world grid
   type(block_data_real8_2d) :: area_smap_wgrid              ! area of each patch across world grid

   ! ensemble predicted observations (at obs grid)
   real(r8), allocatable :: pred_smap_tb_h_ogrid_ens(:, :)   ! predicted H- polarized temp on obs grid
   real(r8), allocatable :: pred_smap_tb_v_ogrid_ens(:, :)   ! predicted V- polarized temp on obs grid
   real(r8) :: area_smap_wgrid_obs                           ! area of each patch across world grid for each obs

   ! observations around patch (dimensions changes with patch)
   integer :: num_smap_p
   logical, allocatable :: index_smap_p(:)                   ! index of obs around each patch
   real(r8), allocatable :: smap_lat_p(:)                    ! latitude of obs around each patch
   real(r8), allocatable :: smap_lon_p(:)                    ! longitude of obs around each patch
   real(r8), allocatable :: smap_tb_h_p(:)                   ! H- polarized brightness temperature of obs around each patch ([K])
   real(r8), allocatable :: smap_tb_v_p(:)                   ! V- polarized brightness temperature of obs around each patch ([K])
   real(r8), allocatable :: pred_smap_tb_h_p_ens(:, :)       ! predicted H- polarized temp around patch
   real(r8), allocatable :: pred_smap_tb_v_p_ens(:, :)       ! predicted V- polarized temp around patch

   ! history mean value of brightness temperature around patch
   real(r8), allocatable :: pred_smap_tb_h_mean_p(:)         ! history mean value of predicted H- polarized temp around patch
   real(r8), allocatable :: pred_smap_tb_v_mean_p(:)         ! history mean value of predicted V- polarized temp around patch
   real(r8), allocatable :: smap_tb_h_mean_p(:)              ! history mean value of SMAP H- polarized temp around patch
   real(r8), allocatable :: smap_tb_v_mean_p(:)              ! history mean value of SMAP V- polarized temp around patch

   ! parameters for distance calculation
   real(r8), allocatable :: d_smap_p(:)                      ! distance between obs and patch center
   real(r8), parameter   :: static_smap_err = 0.1            ! static observation error

!#############################################################################
! For FY3D L1 brightness temperature
!#############################################################################
   ! file path
   character(len=256) :: file_fy3d                           ! FY3D observation file path
   character(len=256) :: file_grid_fy3d                      ! FY3D world grid file path

   ! grid
   type(grid_type) :: grid_fy3d                              ! FY3D world grid
   type(spatial_mapping_type) :: mg2p_fy3d                   ! mapping between FY3D world grid to patch

   ! time variables used to determine whether has obs
   real(r8), allocatable :: fy3d_time(:)                     ! seconds of all obs since begin of current day (UTC)
   real(r8), allocatable :: dt_b_fy3d(:)                     ! delta time between FY3D and begin seconds of current day
   real(r8), allocatable :: dt_e_fy3d(:)                     ! delta time between FY3D and end seconds of current day

   ! logical variables
   logical :: has_fy3d_file                                  ! whether has file of FY3D at target hour
   logical :: has_fy3d_obs                                   ! whether has FY3D obs at current step

   ! observations (dimensions changes with time)
   integer :: num_fy3d_obs                                   ! number of all obs in current file
   integer :: num_fy3d_obs_domain                            ! number of all obs in simulation domain at current file
   real(r8), allocatable :: fy3d_lat(:)                      ! latitude of all obs
   real(r8), allocatable :: fy3d_lon(:)                      ! longitude of all obs
   real(r8), allocatable :: fy3d_tb_h(:)                     ! H- polarized brightness temperature of all obs ([K])
   real(r8), allocatable :: fy3d_tb_v(:)                     ! V- polarized brightness temperature of all obs ([K])
   integer,  allocatable :: fy3d_ii(:)                       ! i-th lat grid in world map of all obs
   integer,  allocatable :: fy3d_jj(:)                       ! j-th lon grid in world map of all obs

   ! history mean value of fy3d and colm
   real(r8), allocatable :: pred_fy3d_tb_h_mean(:)           ! history mean value of predicted H- polarized temp
   real(r8), allocatable :: pred_fy3d_tb_v_mean(:)           ! history mean value of predicted V- polarized temp
   real(r8), allocatable :: fy3d_tb_h_mean(:)                ! history mean value of FY3D H- polarized temp
   real(r8), allocatable :: fy3d_tb_v_mean(:)                ! history mean value of FY3D V- polarized temp

   ! ensemble predicted observations (at patch)
   real(r8), allocatable :: pred_fy3d_tb_h_pset_ens(:, :)    ! predicted H- polarized temp on patch
   real(r8), allocatable :: pred_fy3d_tb_v_pset_ens(:, :)    ! predicted V- polarized temp on patch

   ! ensemble predicted observations (at world grid)
   type(block_data_real8_3d) :: pred_fy3d_tb_h_wgrid_ens     ! predicted H- polarized temp on world grid
   type(block_data_real8_3d) :: pred_fy3d_tb_v_wgrid_ens     ! predicted V- polarized temp on world grid
   type(block_data_real8_2d) :: area_fy3d_wgrid              ! area of each patch across world grid

   ! ensemble predicted observations (at obs grid)
   real(r8), allocatable :: pred_fy3d_tb_h_ogrid_ens(:, :)   ! predicted H- polarized temp on obs grid
   real(r8), allocatable :: pred_fy3d_tb_v_ogrid_ens(:, :)   ! predicted V- polarized temp on obs grid
   real(r8) :: area_fy3d_wgrid_obs                           ! area of each patch across world grid for each obs

   ! observations around patch (dimensions changes with patch)
   integer :: num_fy3d_p
   logical, allocatable :: index_fy3d_p(:)                   ! index of obs around each patch
   real(r8), allocatable :: fy3d_lat_p(:)                    ! latitude of obs around each patch
   real(r8), allocatable :: fy3d_lon_p(:)                    ! longitude of obs around each patch
   real(r8), allocatable :: fy3d_tb_h_p(:)                   ! H- polarized brightness temperature of obs around each patch ([K])
   real(r8), allocatable :: fy3d_tb_v_p(:)                   ! V- polarized brightness temperature of obs around each patch ([K])
   real(r8), allocatable :: pred_fy3d_tb_h_p_ens(:, :)       ! predicted H- polarized temp around patch
   real(r8), allocatable :: pred_fy3d_tb_v_p_ens(:, :)       ! predicted V- polarized temp around patch

   ! history mean value of brightness temperature around patch
   real(r8), allocatable :: pred_fy3d_tb_h_mean_p(:)         ! history mean value of predicted H- polarized temp around patch
   real(r8), allocatable :: pred_fy3d_tb_v_mean_p(:)         ! history mean value of predicted V- polarized temp around patch
   real(r8), allocatable :: fy3d_tb_h_mean_p(:)              ! history mean value of FY3D H- polarized temp around patch
   real(r8), allocatable :: fy3d_tb_v_mean_p(:)              ! history mean value of FY3D V- polarized temp around patch

   ! parameters for distance calculation
   real(r8), allocatable :: d_fy3d_p(:)                      ! distance between obs and patch center
   real(r8), parameter   :: static_fy3d_err = 10.0           ! static observation error

!#############################################################################
! For SYNOP 2m temperature and humidity observations
!#############################################################################
   ! derived types of pixel index in each worker
   ! firstly sorted by pixel latitude index (nlat) and record
   ! corresponding longitude index and patch id for each pixel
   type :: idx_type
      integer, allocatable :: ilon (:)
      integer, allocatable :: ipatch (:)
   END type idx_type
   type(idx_type), allocatable :: idx (:)                    ! derived types of pixel index at each worker
   integer, allocatable :: counter (:)                       ! counter of pixel longitude index of each latitude index at each worker

   ! file path
   character(len=256) :: file_synop                          ! SYNOP file path

   ! info of all SYNOP sites
   character(len=256) :: file_site                           ! file of location of all SYNOP sites
   real(r8), allocatable :: synop_lat_all(:)
   real(r8), allocatable :: synop_lon_all(:)                 ! latitude and longitude of all SYNOP sites
   integer :: nsite                                          ! number of all SYNOP sites
   integer, allocatable :: iloc_synop (:,:)                  ! global lat/lon index of pixel cover each site in all SYNOP sites
   integer :: counter_worker_nsite                           ! number of all SYNOP sites located at each worker
   integer, allocatable :: ip_worker (:,:)                   ! patch id of all SYNOP sites located at each worker
   integer, allocatable :: idx_lat(:), idx_lon(:), pos(:)    ! temporary array save global lat/lon index of each site
   integer, allocatable :: synop_lut (:,:)                   ! look-up-table of all SYNOP sites

   ! time variables used to determine whether has obs
   real(r8), allocatable :: synop_time(:)                    ! seconds of all obs since begin of current day (UTC)
   real(r8), allocatable :: dt_b_synop(:)                    ! delta time between SYNOP and begin seconds of current day
   real(r8), allocatable :: dt_e_synop(:)                    ! delta time between SYNOP and end seconds of current day

   ! logical variables
   logical :: has_synop_file                                 ! whether has file of SYNOP at target day
   logical :: has_synop_obs                                  ! whether has SYNOP obs at current step

   ! info of SYNOP sites at each step
   integer, allocatable :: synop_idx (:,:)                   ! index of SYNOP sites (worker/patch id) at each step
   integer, allocatable :: site_id_worker (:)                ! site id of SYNOP sites at each step at each worker
   real(r8), allocatable :: tref_ens_worker (:,:)            ! predicted 2m temperature of SYNOP sites at each step at each worker
   real(r8), allocatable :: qref_ens_worker (:,:)            ! predicted 2m humidity of SYNOP sites at each step at each worker
   real(r8), allocatable :: qref_ens_o (:,:)                 ! predicted 2m temperature of SYNOP sites at each step
   real(r8), allocatable :: tref_ens_o (:,:)                 ! predicted 2m humidity of SYNOP sites at each step

   ! observations (dimensions changes with time)
   integer :: num_synop_obs                                  ! number of all obs in current file
   integer :: num_synop_obs_domain                           ! number of all obs in simulation domain at current file
   real(r8), allocatable :: synop_lat(:)                     ! latitude of all obs
   real(r8), allocatable :: synop_lon(:)                     ! longitude of all obs
   integer,  allocatable :: synop_id(:)                      ! global id of all obs
   real(r8), allocatable :: synop_tref(:)                    ! 2m temperature of all obs ([K])
   integer,  allocatable :: synop_qref(:)                    ! 2m humidity of all obs ([K])

   ! observations around patch (dimensions changes with patch)
   integer :: num_synop_p
   logical, allocatable :: index_synop_p(:)                  ! index of obs around each patch
   real(r8), allocatable :: synop_lat_p(:)                   ! latitude of obs around each patch
   real(r8), allocatable :: synop_lon_p(:)                   ! longitude of obs around each patch
   real(r8), allocatable :: synop_qref_p(:)                  ! 2m temperature of obs around each patch ([K])
   real(r8), allocatable :: synop_tref_p(:)                  ! 2m humidity of obs around each patch ([K])
   real(r8), allocatable :: synop_p(:)                       ! concatenate 2m temperature and humidity of obs around each patch
   real(r8), allocatable :: qref_ens_p(:,:)                  ! predicted 2m temperature around patch
   real(r8), allocatable :: tref_ens_p(:,:)                  ! predicted 2m humidity around patch
   real(r8), allocatable :: pred_synop_ens_p(:,:)            ! predicted 2m temperature and humidity around patch

   ! parameters for distance calculation
   real(r8), allocatable :: d_synop_p(:)                     ! distance between obs and patch center
   real(r8), parameter   :: static_synop_tref_err = 1.0      ! static observation error (2m temperature)
   real(r8), parameter   :: static_synop_qref_err = 0.04     ! static observation error (2m humidity)

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_DA_SM()

!-----------------------------------------------------------------------------
      USE MOD_Spmd_Task
      USE MOD_Namelist, only: DEF_DA_obsdir
      USE MOD_Grid
      USE MOD_NetCDFSerial
      USE MOD_LandPatch
      USE MOD_Pixelset
      USE MOD_RangeCheck
      IMPLICIT NONE

      integer :: np, ie, ipxstt, ipxend, ipxl, i, ilat, isite, ilon, iwork, mesg(2), isrc, ndata, numpxl_lat, numpxl_lon, numpxl
      integer, allocatable :: temp (:,:)                        ! temporary array for receiving data from workers

!-----------------------------------------------------------------------------

#ifndef SinglePoint
      IF (DEF_DA_SM_SMAP) THEN
         ! grid file path of EASE v2.0, 36km world grid
         file_grid_smap = trim(DEF_DA_obsdir)//'/grid/'//'/SMAP_L1C_36km.nc'
         CALL grid_smap%define_from_file(file_grid_smap, 'latitude', 'longitude')

         ! map SMAP grid to patch
         CALL mg2p_smap%build_arealweighted(grid_smap, landpatch)
      ENDIF

      IF (DEF_DA_SM_FY) THEN
         ! grid file path of FY3D, 0.25 degree world grid
         file_grid_fy3d = trim(DEF_DA_obsdir)//'/grid/'//'/FY3D_L1_0p25.nc'

         CALL grid_fy3d%define_from_file(file_grid_fy3d, 'latitude', 'longitude')

         ! map FY3D grid to patch
         CALL mg2p_fy3d%build_arealweighted(grid_fy3d, landpatch)
      ENDIF

      IF (DEF_DA_SM_SYNOP) THEN
!#############################################################################
! Makeup derived types of pixel index for fast access at each worker
!#############################################################################
         IF (p_is_worker) THEN
            allocate (counter (pixel%nlat))
            counter(:) = 0

            ! count the number of pixel lon index for each pixel lat index
            DO np = 1, numpatch
               ie = landpatch%ielm(np)

               ipxstt = landpatch%ipxstt(np)
               ipxend = landpatch%ipxend(np)

               DO ipxl = ipxstt, ipxend
                  counter(mesh(ie)%ilat(ipxl)) = counter(mesh(ie)%ilat(ipxl)) + 1
               ENDDO
            ENDDO

            ! allocate derived types of index
            allocate (idx (pixel%nlat))
            DO i = 1, pixel%nlat
               IF (counter(i) > 0) THEN
                  allocate (idx(i)%ilon(counter(i)))
                  allocate (idx(i)%ipatch(counter(i)))

                  idx(i)%ilon(:) = 0
                  idx(i)%ipatch(:) = 0
               ENDIF
            ENDDO

            ! fill the index
            counter(:) = 0
            DO np = 1, numpatch
               ie = landpatch%ielm(np)

               ipxstt = landpatch%ipxstt(np)
               ipxend = landpatch%ipxend(np)

               DO ipxl = ipxstt, ipxend
                  ilat = mesh(ie)%ilat(ipxl)
                  counter(mesh(ie)%ilat(ipxl)) = counter(mesh(ie)%ilat(ipxl)) + 1
                  idx(ilat)%ilon(counter(mesh(ie)%ilat(ipxl))) = mesh(ie)%ilon(ipxl)
                  idx(ilat)%ipatch(counter(mesh(ie)%ilat(ipxl))) = np
               ENDDO
            ENDDO
         ENDIF


!#############################################################################
! Read the location of SYNOP sites and find the located pixel of each site
!#############################################################################
         ! file of SYNOP sites
         file_site = trim(DEF_DA_obsdir)//'/grid/'//'/SYNOP.nc'

         ! read latitude and longitude of sites
         IF (ncio_var_exist(file_site, 'latitude')) THEN
            CALL ncio_read_bcast_serial(file_site, 'latitude', synop_lat_all)
            CALL ncio_read_bcast_serial(file_site, 'longitude', synop_lon_all)
         ENDIF
         nsite = size(synop_lat_all)

         ! find the located pixel of each site & broadcast workers
         allocate (iloc_synop (2, nsite))
         iloc_synop(:,:) = -1
         IF (p_is_master) THEN
            DO i = 1, nsite
               numpxl_lat = count(pixel%lat_s <= synop_lat_all(i) .and. pixel%lat_n > synop_lat_all(i))
               numpxl_lon = count(pixel%lon_w <= synop_lon_all(i) .and. pixel%lon_e > synop_lon_all(i))

               IF (numpxl_lat ==1 .and. numpxl_lon ==1) THEN
                  IF (allocated (idx_lat)) deallocate (idx_lat)
                  IF (allocated (idx_lon)) deallocate (idx_lon)
                  allocate (idx_lat(numpxl_lat))
                  allocate (idx_lon(numpxl_lon))

                  idx_lat = pack([(isite, isite=1, pixel%nlat)], pixel%lat_s <= synop_lat_all(i) .and. pixel%lat_n > synop_lat_all(i))
                  idx_lon = pack([(isite, isite=1, pixel%nlon)], pixel%lon_w <= synop_lon_all(i) .and. pixel%lon_e > synop_lon_all(i))

                  iloc_synop(1,i) = idx_lat(1)
                  iloc_synop(2,i) = idx_lon(1)
               ELSE
                  iloc_synop(1,i) = -1
                  iloc_synop(2,i) = -1
               ENDIF
            ENDDO
         ENDIF
         CALL mpi_bcast(iloc_synop, 2*nsite, MPI_INTEGER, p_address_master, p_comm_glb, p_err)

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

!#############################################################################
! Assess the patch id of pixel that cover each observation at each worker
!#############################################################################
         IF (p_is_worker) THEN
            ! count the number of site that located in pixels of each worker
            counter_worker_nsite = 0
            DO i = 1, nsite
               IF (iloc_synop(1, i) > 0) THEN
                  IF (counter(iloc_synop(1, i)) > 0) THEN
                     IF (any(idx(iloc_synop(1, i))%ilon == iloc_synop(2, i))) THEN
                        counter_worker_nsite = counter_worker_nsite + 1
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO

            ! assess patch/site id of sites located at each worker
            IF (counter_worker_nsite > 0) THEN
               allocate (ip_worker (2, counter_worker_nsite))
               counter_worker_nsite = 0
               ip_worker(:,:) = -1

               DO i = 1, nsite
                  IF (iloc_synop(1, i) > 0) THEN
                     IF (counter(iloc_synop(1, i)) > 0) THEN
                        IF (any(idx(iloc_synop(1, i))%ilon == iloc_synop(2, i))) THEN
                           numpxl = count(idx(iloc_synop(1, i))%ilon == iloc_synop(2, i))
                           IF (allocated (pos)) deallocate (pos)
                           allocate (pos(numpxl))

                           pos = pack([(ilon, ilon=1, counter(iloc_synop(1, i)))], (idx(iloc_synop(1, i))%ilon == iloc_synop(2, i)))
                           counter_worker_nsite = counter_worker_nsite + 1
                           ip_worker(1, counter_worker_nsite) = i
                           ip_worker(2, counter_worker_nsite) = idx(iloc_synop(1, i))%ipatch(pos(1))
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF

            ! send the number of site and their patch id to master
#ifdef USEMPI
            mesg = (/p_iam_glb, counter_worker_nsite/)
            CALL mpi_send(mesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)

            IF (counter_worker_nsite > 0) THEN
               CALL mpi_send(ip_worker, 2*counter_worker_nsite, MPI_INTEGER, &
                  p_address_master, mpi_tag_data, p_comm_glb, p_err)
            ENDIF
#endif

            ! deallocate
            IF (allocated (counter)) deallocate (counter)
            DO i = 1, pixel%nlat
               IF (allocated(idx(i)%ilon)) deallocate (idx(i)%ilon)
               IF (allocated(idx(i)%ipatch)) deallocate (idx(i)%ipatch)
            ENDDO
            IF (allocated (idx)) deallocate (idx)
            IF (allocated (iloc_synop)) deallocate (iloc_synop)
            IF (allocated (ip_worker)) deallocate (ip_worker)
         ENDIF

!#############################################################################
! Generate look-up-table (contains the worker/patch id of each site) at master
!#############################################################################
         IF (p_is_master) THEN
            allocate (synop_lut (2, nsite))
            synop_lut(:,:) = -1

#ifdef USEMPI
            DO iwork = 0, p_np_worker-1
               CALL mpi_recv(mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
                  mpi_tag_mesg, p_comm_glb, p_stat, p_err)
               isrc = mesg(1)
               ndata = mesg(2)

               IF (ndata > 0) THEN
                  allocate (temp(2, ndata))
                  CALL mpi_recv(temp, 2*ndata, MPI_INTEGER, isrc, &
                     mpi_tag_data, p_comm_glb, p_stat, p_err)

                  DO i = 1, ndata
                     synop_lut(1, temp(1,i)) = isrc
                     synop_lut(2, temp(1,i)) = temp(2,i)
                  ENDDO

                  deallocate(temp)
               ENDIF
            ENDDO
#endif
         ENDIF
      ENDIF
#endif

   END SUBROUTINE init_DA_SM

!-----------------------------------------------------------------------------

   SUBROUTINE run_DA_SM(idate, deltim)

!-----------------------------------------------------------------------------
      USE MOD_Spmd_task
      USE MOD_TimeManager
      USE MOD_NetCDFBlock
      USE MOD_Mesh
      USE MOD_LandElm
      USE MOD_LandPatch
      USE MOD_Vars_Global
      USE MOD_Vars_1DFluxes
      USE MOD_Vars_1DForcing
      USE MOD_Vars_TimeVariables
      USE MOD_Vars_TimeInvariants
      USE MOD_DA_Vars_TimeVariables
      USE MOD_RangeCheck
      USE MOD_UserDefFun
      USE MOD_DA_EnKF
      USE MOD_Const_Physical, only: denice, denh2o
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      integer, intent(in) :: idate(3)
      real(r8), intent(in) :: deltim

!------------------------ Local Variables ------------------------------------
      real(r8) :: lat_p_n, lat_p_s, lon_p_w, lon_p_e
      real(r8) :: lat_p, lon_p
      real(r8) :: es, e
      integer  :: ib, jb, il, jl, ip, iens, iobs, np, i, n, ndata, isrc, j
      integer  :: sdate(3), smesg(2), rmesg(2)
      integer, allocatable  :: iloc(:)
      integer, allocatable  :: itemp(:)
      real(r8), allocatable :: dtemp(:,:)
      integer  :: iwork, mesg(2)

!-----------------------------------------------------------------------------

!#############################################################################
! Identify if there are observations at this time step
!#############################################################################
      ! Do not perform DA, only calcuate predict BRT for diagnostic
      IF (DEF_DA_ENS_NUM == 1) THEN
         has_smap_file = .false.
         has_fy3d_file = .false.
         has_synop_file = .false.
      ELSE
         ! covert local time to UTC for single point
         sdate = idate
         CALL adj2begin(sdate)
#ifdef SinglePoint
         IF (.not. DEF_simulation_time%greenwich) THEN
            CALL localtime2gmt(sdate)
         ENDIF
#endif
         ! calculate year/month/day/hour of current step
         CALL julian2monthday(sdate(1), sdate(2), month, mday)
         hour = int(sdate(3)/3600)

         ! whether has file at target day
         write (yearstr, '(I4.4)') sdate(1)
         write (monthstr, '(I2.2)') month
         write (daystr, '(I2.2)') mday
         write (hourstr, '(I2.2)') hour

         ! check SMAP brightness temperature observations
         IF (DEF_DA_SM_SMAP) THEN
            ! file path of SMAP
            file_smap = trim(DEF_DA_obsdir)//'/SMAP_L1C_D/'//'/SMAP_L1C_TB_'// &
               trim(yearstr)//'_'//trim(monthstr)//'_'//trim(daystr)//'_'//trim(hourstr)//'.nc'
            inquire (file=trim(file_smap), exist=has_smap_file)

            ! whether have obs at this time interval
            has_smap_obs = .false.
            IF (has_smap_file) THEN
               CALL ncio_read_bcast_serial(file_smap, 'time', smap_time)
               num_smap_obs = size(smap_time)
               idate_b = sdate(3)
               idate_e = sdate(3) + deltim
               allocate (dt_b_smap(num_smap_obs))
               allocate (dt_e_smap(num_smap_obs))
               dt_b_smap = smap_time - idate_b
               dt_e_smap = smap_time - idate_e
               IF (any(dt_b_smap >= 0 .and. dt_e_smap <= 0)) has_smap_obs = .true.
               deallocate (dt_b_smap)
               deallocate (dt_e_smap)
            ELSE
               has_smap_obs = .false.
            ENDIF

            ! check if there are obs in simulation domain
            IF (has_smap_obs) THEN
               CALL ncio_read_bcast_serial(file_smap, 'lat', smap_lat)
               CALL ncio_read_bcast_serial(file_smap, 'lon', smap_lon)

               num_smap_obs_domain = count(smap_lat(:) < DEF_domain%edgen .and. smap_lat(:) > DEF_domain%edges .and. &
                        smap_lon(:) > DEF_domain%edgew .and. smap_lon(:) < DEF_domain%edgee)
               IF (num_smap_obs_domain == 0) has_smap_obs = .false.
            ENDIF
         ELSE
            has_smap_file = .false.
            has_smap_obs = .false.
         ENDIF

         ! check FY brightness temperature observations
         IF (DEF_DA_SM_FY) THEN
            ! file path of FY3D
            file_fy3d = trim(DEF_DA_obsdir)//'/FY3D_L1/'//'/FY3D_L1_TB_'// &
               trim(yearstr)//'_'//trim(monthstr)//'_'//trim(daystr)//'_'//trim(hourstr)//'.nc'
            inquire (file=trim(file_fy3d), exist=has_fy3d_file)

            ! whether have obs at this time interval
            has_fy3d_obs = .false.
            IF (has_fy3d_file) THEN
               CALL ncio_read_bcast_serial(file_fy3d, 'time', fy3d_time)
               num_fy3d_obs = size(fy3d_time)
               idate_b = sdate(3)
               idate_e = sdate(3) + deltim
               allocate (dt_b_fy3d(num_fy3d_obs))
               allocate (dt_e_fy3d(num_fy3d_obs))
               dt_b_fy3d = fy3d_time - idate_b
               dt_e_fy3d = fy3d_time - idate_e
               IF (any(dt_b_fy3d >= 0 .and. dt_e_fy3d <= 0)) has_fy3d_obs = .true.
               deallocate (dt_b_fy3d)
               deallocate (dt_e_fy3d)
            ELSE
               has_fy3d_obs = .false.
            ENDIF

            ! check if there are obs in simulation domain
            IF (has_fy3d_obs) THEN
               CALL ncio_read_bcast_serial(file_fy3d, 'lat', fy3d_lat)
               CALL ncio_read_bcast_serial(file_fy3d, 'lon', fy3d_lon)

               num_fy3d_obs_domain = count(fy3d_lat(:) < DEF_domain%edgen .and. fy3d_lat(:) > DEF_domain%edges .and. &
                        fy3d_lon(:) > DEF_domain%edgew .and. fy3d_lon(:) < DEF_domain%edgee)
               IF (num_fy3d_obs_domain == 0) has_fy3d_obs = .false.
            ENDIF
         ELSE
            has_fy3d_file = .false.
            has_fy3d_obs = .false.
         ENDIF

         ! check SYNOP 2m observations
         IF (DEF_DA_SM_SYNOP) THEN
            file_synop = trim(DEF_DA_obsdir)//'/SYNOP/'//'/SYNOP_'// &
               trim(yearstr)//'_'//trim(monthstr)//'_'//trim(daystr)//'_'//trim(hourstr)//'.nc'
            inquire (file=trim(file_synop), exist=has_synop_file)

            ! whether have obs at this time interval
            has_synop_obs = .false.
            IF (has_synop_file) THEN
               CALL ncio_read_bcast_serial(file_synop, 'time', synop_time)
               num_synop_obs = size(synop_time)
               idate_b = sdate(3)
               idate_e = sdate(3) + deltim
               allocate (dt_b_synop(num_synop_obs))
               allocate (dt_e_synop(num_synop_obs))
               dt_b_synop = synop_time - idate_b
               dt_e_synop = synop_time - idate_e
               IF (any(dt_b_synop >= 0 .and. dt_e_synop <= 0)) has_synop_obs = .true.
               deallocate (dt_b_synop)
               deallocate (dt_e_synop)
            ELSE
               has_synop_obs = .false.
            ENDIF

            ! check if there are obs in simulation domain
            IF (has_synop_obs) THEN
               CALL ncio_read_bcast_serial(file_synop, 'lat', synop_lat)
               CALL ncio_read_bcast_serial(file_synop, 'lon', synop_lon)

               num_synop_obs_domain = count(synop_lat(:) < DEF_domain%edgen .and. synop_lat(:) > DEF_domain%edges .and. &
                        synop_lon(:) > DEF_domain%edgew .and. synop_lon(:) < DEF_domain%edgee)
               IF (num_synop_obs_domain == 0) has_synop_obs = .false.
            ENDIF
         ELSE
            has_synop_file = .false.
            has_synop_obs = .false.
         ENDIF
      ENDIF

      ! print info of observations
      IF (p_is_master) THEN
         IF (has_smap_obs) THEN
            print *, '[CoLM-DA] Have SMAP observations:', trim(file_smap)
         ELSE
            print *, '[CoLM-DA] No SMAP observations.'
         ENDIF
         IF (has_fy3d_obs) THEN
            print *, '[CoLM-DA] Have FY3D observations:', trim(file_fy3d)
         ELSE
            print *, '[CoLM-DA] No FY3D observations.'
         ENDIF
         IF (has_synop_obs) THEN
            print *, '[CoLM-DA] Have SYNOP observations:', trim(file_synop)
         ELSE
            print *, '[CoLM-DA] No SYNOP observations.'
         ENDIF
      ENDIF

!#############################################################################
! Allocate memory for variables
!#############################################################################
      IF (DEF_DA_SM_SMAP) THEN
         ! allocate memory for ensemble predicted observations at patch (for DA or no DA)
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               IF (allocated(pred_smap_tb_h_pset_ens)) deallocate (pred_smap_tb_h_pset_ens)
               IF (allocated(pred_smap_tb_v_pset_ens)) deallocate (pred_smap_tb_v_pset_ens)
               allocate (pred_smap_tb_h_pset_ens(DEF_DA_ENS_NUM, numpatch))
               allocate (pred_smap_tb_v_pset_ens(DEF_DA_ENS_NUM, numpatch))
            ENDIF
         ENDIF

         ! allocate memory for ensemble predicted observations in vectors (only for DA)
#ifndef SinglePoint
         IF (has_smap_obs) THEN
            IF (allocated(pred_smap_tb_h_ogrid_ens)) deallocate (pred_smap_tb_h_ogrid_ens)
            IF (allocated(pred_smap_tb_v_ogrid_ens)) deallocate (pred_smap_tb_v_ogrid_ens)
            allocate (pred_smap_tb_h_ogrid_ens(num_smap_obs, DEF_DA_ENS_NUM))
            allocate (pred_smap_tb_v_ogrid_ens(num_smap_obs, DEF_DA_ENS_NUM))
         ENDIF
#endif
      ENDIF

      IF (DEF_DA_SM_FY) THEN
         ! allocate memory for ensemble predicted observations at patch (for DA or no DA)
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               IF (allocated(pred_fy3d_tb_h_pset_ens)) deallocate (pred_fy3d_tb_h_pset_ens)
               IF (allocated(pred_fy3d_tb_v_pset_ens)) deallocate (pred_fy3d_tb_v_pset_ens)
               allocate (pred_fy3d_tb_h_pset_ens(DEF_DA_ENS_NUM, numpatch))
               allocate (pred_fy3d_tb_v_pset_ens(DEF_DA_ENS_NUM, numpatch))
            ENDIF
         ENDIF

         ! allocate memory for ensemble predicted observations in vectors (only for DA)
#ifndef SinglePoint
         IF (has_fy3d_obs) THEN
            IF (allocated(pred_fy3d_tb_h_ogrid_ens)) deallocate (pred_fy3d_tb_h_ogrid_ens)
            IF (allocated(pred_fy3d_tb_v_ogrid_ens)) deallocate (pred_fy3d_tb_v_ogrid_ens)
            allocate (pred_fy3d_tb_h_ogrid_ens(num_fy3d_obs, DEF_DA_ENS_NUM))
            allocate (pred_fy3d_tb_v_ogrid_ens(num_fy3d_obs, DEF_DA_ENS_NUM))
         ENDIF
#endif
      ENDIF

      IF (DEF_DA_SM_SYNOP) THEN
         ! allocate memory
         IF (has_synop_obs) THEN
            IF (allocated(synop_idx)) deallocate (synop_idx)
            IF (allocated(qref_ens_o)) deallocate (qref_ens_o)
            IF (allocated(tref_ens_o)) deallocate (tref_ens_o)
            allocate (synop_idx (2, num_synop_obs))
            allocate (qref_ens_o (num_synop_obs, DEF_DA_ENS_NUM))
            allocate (tref_ens_o (num_synop_obs, DEF_DA_ENS_NUM))
         ENDIF
      ENDIF

!#############################################################################
! Calculate predicted observations using observation operator
!#############################################################################
      ! forward model (for no DA)
      IF (p_is_worker) THEN
         ! SMAP forward model
         DO np = 1, numpatch
            lat_p   = patchlatr(np)*180/pi
            lon_p   = patchlonr(np)*180/pi
            CALL forward(&
               patchtype(np), patchclass(np), dz_sno(:,np), &
               forc_topo(np), htop(np), &
               tref(np), t_soisno(:,np), tleaf(np), &
               wliq_soisno(:,np), wice_soisno(:,np), h2osoi(:,np), &
               snowdp(np), lai(np), sai(np), &
               wf_clay(:, np), wf_sand(:, np), wf_silt(:, np), BD_all(:, np), porsl(:, np), &
               smap_theta, smap_fghz, &
               t_brt_smap(1,np), t_brt_smap(2,np))
         ENDDO

         ! FY3D forward model
         DO np = 1, numpatch
            CALL forward(&
               patchtype(np), patchclass(np), dz_sno(:,np), &
               forc_topo(np), htop(np), &
               tref(np), t_soisno(:,np), tleaf(np), &
               wliq_soisno(:,np), wice_soisno(:,np), h2osoi(:,np), &
               snowdp(np), lai(np), sai(np), &
               wf_clay(:, np), wf_sand(:, np), wf_silt(:, np), BD_all(:, np), porsl(:, np), &
               fy3d_theta, fy3d_fghz, &
               t_brt_fy3d(1,np), t_brt_fy3d(2,np))
         ENDDO
      ENDIF

      ! forward model for ensemble DA
      IF (p_is_worker) THEN
         IF (DEF_DA_SM_SMAP) THEN
            IF (DEF_DA_ENS_NUM > 1) THEN
               DO iens = 1, DEF_DA_ENS_NUM
                  DO np = 1, numpatch
                     CALL forward(&
                        patchtype(np), patchclass(np), dz_sno_ens(:, iens, np), &
                        forc_topo(np), htop(np), &
                        tref_ens(iens, np), t_soisno_ens(:,iens,np), tleaf_ens(iens, np), &
                        wliq_soisno_ens(:, iens, np), wice_soisno_ens(:, iens, np), h2osoi_ens(:, iens, np), &
                        snowdp_ens(iens, np), lai_ens(iens, np), sai_ens(iens, np), &
                        wf_clay(:, np), wf_sand(:, np), wf_silt(:, np), BD_all(:, np), porsl(:, np), &
                        smap_theta, smap_fghz, &
                        pred_smap_tb_h_pset_ens(iens, np), pred_smap_tb_v_pset_ens(iens, np))
                  ENDDO
               ENDDO
               t_brt_smap_ens(1,:,:) = pred_smap_tb_h_pset_ens
               t_brt_smap_ens(2,:,:) = pred_smap_tb_v_pset_ens
            ENDIF
         ENDIF

         IF (DEF_DA_SM_FY) THEN
            IF (DEF_DA_ENS_NUM > 1) THEN
               DO iens = 1, DEF_DA_ENS_NUM
                  DO np = 1, numpatch
                     CALL forward(&
                        patchtype(np), patchclass(np), dz_sno_ens(:, iens, np), &
                        forc_topo(np), htop(np), &
                        tref_ens(iens, np), t_soisno_ens(:,iens,np), tleaf_ens(iens, np), &
                        wliq_soisno_ens(:, iens, np), wice_soisno_ens(:, iens, np), h2osoi_ens(:, iens, np), &
                        snowdp_ens(iens, np), lai_ens(iens, np), sai_ens(iens, np), &
                        wf_clay(:, np), wf_sand(:, np), wf_silt(:, np), BD_all(:, np), porsl(:, np), &
                        fy3d_theta, fy3d_fghz, &
                        pred_fy3d_tb_h_pset_ens(iens, np), pred_fy3d_tb_v_pset_ens(iens, np))
                  ENDDO
               ENDDO
               t_brt_fy3d_ens(1,:,:) = pred_fy3d_tb_h_pset_ens
               t_brt_fy3d_ens(2,:,:) = pred_fy3d_tb_v_pset_ens
            ENDIF
         ENDIF

         IF (DEF_DA_SM_SYNOP) THEN
            IF (DEF_DA_ENS_NUM > 1) THEN
               ! calculate relative humidity from specific humidity
               DO ip = 1, numpatch
                  DO iens = 1, DEF_DA_ENS_NUM
                     es = 6.112 * exp((17.67 * (tref_ens(iens, ip) - 273.15))/ (tref_ens(iens, ip) - 273.15 + 243.5))
                     e = (qref_ens(iens, ip) * forc_psrf(ip) / 100) / (0.378 * qref_ens(iens, ip) + 0.622)
                     rhref_ens(iens, ip) = max(min(1.0d0, e / es), 0.0d0)
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
      ENDIF

!#############################################################################
! Reading observations data
!#############################################################################
      ! read observations data (only for DA)
      IF (has_smap_obs) THEN
         CALL ncio_read_bcast_serial(file_smap, 'ii'  , smap_ii  )
         CALL ncio_read_bcast_serial(file_smap, 'jj'  , smap_jj  )
         CALL ncio_read_bcast_serial(file_smap, 'tb_h', smap_tb_h)
         CALL ncio_read_bcast_serial(file_smap, 'tb_v', smap_tb_v)
         CALL ncio_read_bcast_serial(file_smap, 'tb_h_colm_mean', pred_smap_tb_h_mean)
         CALL ncio_read_bcast_serial(file_smap, 'tb_v_colm_mean', pred_smap_tb_v_mean)
         CALL ncio_read_bcast_serial(file_smap, 'tb_h_obs_mean',  smap_tb_h_mean)
         CALL ncio_read_bcast_serial(file_smap, 'tb_v_obs_mean',  smap_tb_v_mean)
      ENDIF
      IF (has_fy3d_obs) THEN
         CALL ncio_read_bcast_serial(file_fy3d, 'ii'  , fy3d_ii  )
         CALL ncio_read_bcast_serial(file_fy3d, 'jj'  , fy3d_jj  )
         CALL ncio_read_bcast_serial(file_fy3d, 'tb_h', fy3d_tb_h)
         CALL ncio_read_bcast_serial(file_fy3d, 'tb_v', fy3d_tb_v)
         CALL ncio_read_bcast_serial(file_fy3d, 'tb_h_colm_mean', pred_fy3d_tb_h_mean)
         CALL ncio_read_bcast_serial(file_fy3d, 'tb_v_colm_mean', pred_fy3d_tb_v_mean)
         CALL ncio_read_bcast_serial(file_fy3d, 'tb_h_obs_mean',  fy3d_tb_h_mean)
         CALL ncio_read_bcast_serial(file_fy3d, 'tb_v_obs_mean',  fy3d_tb_v_mean)
      ENDIF
      IF (has_synop_obs) THEN
         CALL ncio_read_bcast_serial(file_synop, 'lat',  synop_lat )
         CALL ncio_read_bcast_serial(file_synop, 'lon',  synop_lon )
         CALL ncio_read_bcast_serial(file_synop, 'id',   synop_id  )
         CALL ncio_read_bcast_serial(file_synop, 'tref', synop_tref)
         CALL ncio_read_bcast_serial(file_synop, 'qref', synop_qref)
      ENDIF

!#############################################################################
! Cropping predicted observations at observations space
!#############################################################################
#ifndef SinglePoint
      IF (DEF_DA_SM_SMAP) THEN
         IF (has_smap_obs) THEN
            ! prepare filter for land patches
            allocate (filter(numpatch))
            IF (p_is_worker) THEN
               IF (numpatch > 0) THEN
                  filter(:) = patchtype <= 2
               ENDIF
            ENDIF

            ! calculate area of each patch across world grid
            IF (p_is_io) CALL allocate_block_data(grid_smap, area_smap_wgrid)
            CALL mg2p_smap%get_sumarea(area_smap_wgrid, filter)

            ! mapping predicted observations from patch to world grid
            IF (p_is_io) THEN
               CALL allocate_block_data(grid_smap, pred_smap_tb_h_wgrid_ens, DEF_DA_ENS_NUM)
               CALL allocate_block_data(grid_smap, pred_smap_tb_v_wgrid_ens, DEF_DA_ENS_NUM)
            ENDIF
            CALL mg2p_smap%pset2grid(pred_smap_tb_h_pset_ens, pred_smap_tb_h_wgrid_ens, spv=spval, msk=filter)
            CALL mg2p_smap%pset2grid(pred_smap_tb_v_pset_ens, pred_smap_tb_v_wgrid_ens, spv=spval, msk=filter)
            deallocate (filter)

            ! crop the predicted observations
            IF (p_is_io) THEN
               allocate (iloc(num_smap_obs))
               pred_smap_tb_h_ogrid_ens = -9999.0
               pred_smap_tb_v_ogrid_ens = -9999.0

               ! crop obs from world grid to all obs grids
               ndata = 0
               DO i = 1, num_smap_obs
                  ib = grid_smap%xblk(smap_jj(i))
                  jb = grid_smap%yblk(smap_ii(i))
                  il = grid_smap%xloc(smap_jj(i))
                  jl = grid_smap%yloc(smap_ii(i))
                  IF (ib /= 0 .and. jb /= 0) THEN
                     IF (gblock%pio(ib, jb) == p_iam_glb) THEN
                        area_smap_wgrid_obs = area_smap_wgrid%blk(ib, jb)%val(il, jl)
                        IF (area_smap_wgrid_obs /= 0) THEN
#ifdef USEMPI
                           ndata = ndata + 1
                           iloc(ndata) = i
                           pred_smap_tb_h_ogrid_ens(ndata, :) = (pred_smap_tb_h_wgrid_ens%blk(ib, jb)%val(:, il, jl))/area_smap_wgrid_obs
                           pred_smap_tb_v_ogrid_ens(ndata, :) = (pred_smap_tb_v_wgrid_ens%blk(ib, jb)%val(:, il, jl))/area_smap_wgrid_obs
#else
                           pred_smap_tb_h_ogrid_ens(i, :) = (pred_smap_tb_h_wgrid_ens%blk(ib, jb)%val(:, il, jl))/area_smap_wgrid_obs
                           pred_smap_tb_v_ogrid_ens(i, :) = (pred_smap_tb_v_wgrid_ens%blk(ib, jb)%val(:, il, jl))/area_smap_wgrid_obs
#endif
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO

               ! send data from io to masters
#ifdef USEMPI
               smesg = (/p_iam_glb, ndata/)
               CALL mpi_send(smesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)

               IF (ndata > 0) THEN
                  CALL mpi_send(iloc(1:ndata), ndata, MPI_INTEGER, p_address_master, mpi_tag_data, p_comm_glb, p_err)

                  allocate (dtemp(ndata, DEF_DA_ENS_NUM))
                  dtemp = pred_smap_tb_h_ogrid_ens(1:ndata, :)
                  CALL mpi_send(dtemp, ndata*DEF_DA_ENS_NUM, MPI_REAL8, p_address_master, mpi_tag_data, p_comm_glb, p_err)
                  dtemp = pred_smap_tb_v_ogrid_ens(1:ndata, :)
                  CALL mpi_send(dtemp, ndata*DEF_DA_ENS_NUM, MPI_REAL8, p_address_master, mpi_tag_data + 1, p_comm_glb, p_err)
                  deallocate (dtemp)
               ENDIF
#endif
               deallocate (iloc)
               deallocate (pred_smap_tb_h_wgrid_ens%blk, pred_smap_tb_v_wgrid_ens%blk)
            ENDIF

            ! broadcast from master to all workers
#ifdef USEMPI
            IF (p_is_master) THEN
               pred_smap_tb_h_ogrid_ens = -9999.0
               pred_smap_tb_v_ogrid_ens = -9999.0
               DO ip = 0, p_np_io - 1
                  CALL mpi_recv(rmesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

                  ndata = rmesg(2)
                  IF (ndata > 0) THEN
                     allocate (itemp(ndata))
                     allocate (dtemp(ndata, DEF_DA_ENS_NUM))

                     isrc = rmesg(1)
                     CALL mpi_recv(itemp, ndata, MPI_INTEGER, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
                     CALL mpi_recv(dtemp, ndata*DEF_DA_ENS_NUM, MPI_REAL8, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
                     pred_smap_tb_h_ogrid_ens(itemp, :) = dtemp
                     CALL mpi_recv(dtemp, ndata*DEF_DA_ENS_NUM, MPI_REAL8, isrc, mpi_tag_data + 1, p_comm_glb, p_stat, p_err)
                     pred_smap_tb_v_ogrid_ens(itemp, :) = dtemp

                     deallocate (itemp)
                     deallocate (dtemp)
                  ENDIF
               ENDDO
            ENDIF
            CALL mpi_bcast(pred_smap_tb_h_ogrid_ens, num_smap_obs*DEF_DA_ENS_NUM, MPI_REAL8, p_address_master, p_comm_glb, p_err)
            CALL mpi_bcast(pred_smap_tb_v_ogrid_ens, num_smap_obs*DEF_DA_ENS_NUM, MPI_REAL8, p_address_master, p_comm_glb, p_err)
#endif
         ENDIF
      ENDIF

      IF (DEF_DA_SM_FY) THEN
         IF (has_fy3d_obs) THEN
            ! prepare filter for land patches
            allocate (filter(numpatch))
            IF (p_is_worker) THEN
               IF (numpatch > 0) THEN
                  filter(:) = patchtype <= 2
               ENDIF
            ENDIF

            ! calculate area of each patch across world grid
            IF (p_is_io) CALL allocate_block_data(grid_fy3d, area_fy3d_wgrid)
            CALL mg2p_fy3d%get_sumarea(area_fy3d_wgrid, filter)

            ! mapping predicted observations from patch to world grid
            IF (p_is_io) THEN
               CALL allocate_block_data(grid_fy3d, pred_fy3d_tb_h_wgrid_ens, DEF_DA_ENS_NUM)
               CALL allocate_block_data(grid_fy3d, pred_fy3d_tb_v_wgrid_ens, DEF_DA_ENS_NUM)
            ENDIF
            CALL mg2p_fy3d%pset2grid(pred_fy3d_tb_h_pset_ens, pred_fy3d_tb_h_wgrid_ens, spv=spval, msk=filter)
            CALL mg2p_fy3d%pset2grid(pred_fy3d_tb_v_pset_ens, pred_fy3d_tb_v_wgrid_ens, spv=spval, msk=filter)
            deallocate (filter)

            ! crop the predicted observations
            IF (p_is_io) THEN
               allocate (iloc(num_fy3d_obs))
               pred_fy3d_tb_h_ogrid_ens = -9999.0
               pred_fy3d_tb_v_ogrid_ens = -9999.0

               ! crop obs from world grid to all obs grids
               ndata = 0
               DO i = 1, num_fy3d_obs
                  ib = grid_fy3d%xblk(fy3d_jj(i))
                  jb = grid_fy3d%yblk(fy3d_ii(i))
                  il = grid_fy3d%xloc(fy3d_jj(i))
                  jl = grid_fy3d%yloc(fy3d_ii(i))
                  IF (ib /= 0 .and. jb /= 0) THEN
                     IF (gblock%pio(ib, jb) == p_iam_glb) THEN
                        area_fy3d_wgrid_obs = area_fy3d_wgrid%blk(ib, jb)%val(il, jl)
                        IF (area_fy3d_wgrid_obs /= 0) THEN
#ifdef USEMPI
                           ndata = ndata + 1
                           iloc(ndata) = i
                           pred_fy3d_tb_h_ogrid_ens(ndata, :) = (pred_fy3d_tb_h_wgrid_ens%blk(ib, jb)%val(:, il, jl))/area_fy3d_wgrid_obs
                           pred_fy3d_tb_v_ogrid_ens(ndata, :) = (pred_fy3d_tb_v_wgrid_ens%blk(ib, jb)%val(:, il, jl))/area_fy3d_wgrid_obs
#else
                           pred_fy3d_tb_h_ogrid_ens(i, :) = (pred_fy3d_tb_h_wgrid_ens%blk(ib, jb)%val(:, il, jl))/area_fy3d_wgrid_obs
                           pred_fy3d_tb_v_ogrid_ens(i, :) = (pred_fy3d_tb_v_wgrid_ens%blk(ib, jb)%val(:, il, jl))/area_fy3d_wgrid_obs
#endif
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO

               ! send data from io to masters
#ifdef USEMPI
               smesg = (/p_iam_glb, ndata/)
               CALL mpi_send(smesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)

               IF (ndata > 0) THEN
                  CALL mpi_send(iloc(1:ndata), ndata, MPI_INTEGER, p_address_master, mpi_tag_data, p_comm_glb, p_err)

                  allocate (dtemp(ndata, DEF_DA_ENS_NUM))
                  dtemp = pred_fy3d_tb_h_ogrid_ens(1:ndata, :)
                  CALL mpi_send(dtemp, ndata*DEF_DA_ENS_NUM, MPI_REAL8, p_address_master, mpi_tag_data, p_comm_glb, p_err)
                  dtemp = pred_fy3d_tb_v_ogrid_ens(1:ndata, :)
                  CALL mpi_send(dtemp, ndata*DEF_DA_ENS_NUM, MPI_REAL8, p_address_master, mpi_tag_data + 1, p_comm_glb, p_err)
                  deallocate (dtemp)
               ENDIF
#endif
               deallocate (iloc)
               deallocate (pred_fy3d_tb_h_wgrid_ens%blk, pred_fy3d_tb_v_wgrid_ens%blk)
            ENDIF

            ! broadcast from master to all workers
#ifdef USEMPI
            IF (p_is_master) THEN
               pred_fy3d_tb_h_ogrid_ens = -9999.0
               pred_fy3d_tb_v_ogrid_ens = -9999.0
               DO ip = 0, p_np_io - 1
                  CALL mpi_recv(rmesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

                  ndata = rmesg(2)
                  IF (ndata > 0) THEN
                     allocate (itemp(ndata))
                     allocate (dtemp(ndata, DEF_DA_ENS_NUM))

                     isrc = rmesg(1)
                     CALL mpi_recv(itemp, ndata, MPI_INTEGER, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
                     CALL mpi_recv(dtemp, ndata*DEF_DA_ENS_NUM, MPI_REAL8, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
                     pred_fy3d_tb_h_ogrid_ens(itemp, :) = dtemp
                     CALL mpi_recv(dtemp, ndata*DEF_DA_ENS_NUM, MPI_REAL8, isrc, mpi_tag_data + 1, p_comm_glb, p_stat, p_err)
                     pred_fy3d_tb_v_ogrid_ens(itemp, :) = dtemp

                     deallocate (itemp)
                     deallocate (dtemp)
                  ENDIF
               ENDDO
            ENDIF
            CALL mpi_bcast(pred_fy3d_tb_h_ogrid_ens, num_fy3d_obs*DEF_DA_ENS_NUM, MPI_REAL8, p_address_master, p_comm_glb, p_err)
            CALL mpi_bcast(pred_fy3d_tb_v_ogrid_ens, num_fy3d_obs*DEF_DA_ENS_NUM, MPI_REAL8, p_address_master, p_comm_glb, p_err)
#endif
         ENDIF
      ENDIF

      IF (DEF_DA_SM_SYNOP) THEN
         ! crop corresponding index (worker id and patch id) of each observation
         IF (has_synop_obs) THEN
            synop_idx(:,:) = -1
            IF (p_is_master) THEN
               synop_idx = synop_lut(:, synop_id)
            ENDIF
            CALL mpi_bcast(synop_idx, 2*num_synop_obs, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
         ENDIF

         ! recount the number of observations at each worker
         IF (has_synop_obs) THEN
            IF (p_is_worker) THEN
               counter_worker_nsite = 0
               DO i = 1, num_synop_obs
                  IF (synop_idx(1, i) == p_iam_glb) THEN
                     counter_worker_nsite = counter_worker_nsite + 1
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         ! allocate memory for worker
         IF (has_synop_obs) THEN
            IF (p_is_worker) THEN
               IF (counter_worker_nsite > 0) THEN
                  IF (allocated(tref_ens_worker)) deallocate (tref_ens_worker)
                  IF (allocated(qref_ens_worker)) deallocate (qref_ens_worker)
                  IF (allocated(site_id_worker)) deallocate (site_id_worker)
                  allocate (tref_ens_worker (counter_worker_nsite, DEF_DA_ENS_NUM))
                  allocate (qref_ens_worker (counter_worker_nsite, DEF_DA_ENS_NUM))
                  allocate (site_id_worker  (counter_worker_nsite))
               ENDIF
            ENDIF
         ENDIF

         ! crop observations at each worker according index & send to master
         IF (has_synop_obs) THEN
            IF (p_is_worker) THEN
               counter_worker_nsite = 0
               DO i = 1, num_synop_obs
                  IF (synop_idx(1, i) == p_iam_glb) THEN
                     counter_worker_nsite = counter_worker_nsite + 1
                     tref_ens_worker(counter_worker_nsite,:) = tref_ens(:, synop_idx(2, i))
                     qref_ens_worker(counter_worker_nsite,:) = rhref_ens(:, synop_idx(2, i))
                     site_id_worker (counter_worker_nsite)   = i
                  ENDIF
               ENDDO

#ifdef USEMPI
               ! send the number of site and their patch id to master
               mesg = (/p_iam_glb, counter_worker_nsite/)
               CALL mpi_send(mesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)

               IF (counter_worker_nsite > 0) THEN
                  CALL mpi_send(tref_ens_worker, DEF_DA_ENS_NUM*counter_worker_nsite, MPI_REAL8, p_address_master, mpi_tag_data, p_comm_glb, p_err)
                  CALL mpi_send(qref_ens_worker, DEF_DA_ENS_NUM*counter_worker_nsite, MPI_REAL8, p_address_master, mpi_tag_data, p_comm_glb, p_err)
                  CALL mpi_send(site_id_worker, counter_worker_nsite, MPI_INTEGER, p_address_master, mpi_tag_data, p_comm_glb, p_err)
               ENDIF
#endif

               IF (allocated(tref_ens_worker)) deallocate (tref_ens_worker)
               IF (allocated(qref_ens_worker)) deallocate (qref_ens_worker)
               IF (allocated(site_id_worker)) deallocate (site_id_worker)
            ENDIF
         ENDIF

         ! concatenate all predicted observations at master & broadcast to all workers
         IF (has_synop_obs) THEN
            qref_ens_o = spval
            tref_ens_o = spval
            IF (p_is_master) THEN
#ifdef USEMPI
               DO iwork = 0, p_np_worker-1
                  CALL mpi_recv(mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
                  isrc = mesg(1)
                  ndata = mesg(2)

                  IF (ndata > 0) THEN
                     allocate(tref_ens_worker(ndata, DEF_DA_ENS_NUM))
                     allocate(qref_ens_worker(ndata, DEF_DA_ENS_NUM))
                     allocate(site_id_worker (ndata))

                     CALL mpi_recv(tref_ens_worker, ndata*DEF_DA_ENS_NUM, MPI_REAL8, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
                     CALL mpi_recv(qref_ens_worker, ndata*DEF_DA_ENS_NUM, MPI_REAL8, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
                     CALL mpi_recv(site_id_worker, ndata, MPI_INTEGER, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

                     qref_ens_o(site_id_worker,:) = qref_ens_worker
                     tref_ens_o(site_id_worker,:) = tref_ens_worker
                  ENDIF

                  IF (allocated(tref_ens_worker)) deallocate(tref_ens_worker)
                  IF (allocated(qref_ens_worker)) deallocate(qref_ens_worker)
                  IF (allocated(site_id_worker)) deallocate(site_id_worker)
               ENDDO
#endif
            ENDIF
            CALL mpi_bcast(qref_ens_o, DEF_DA_ENS_NUM*num_synop_obs, MPI_REAL8, p_address_master, p_comm_glb, p_err)
            CALL mpi_bcast(tref_ens_o, DEF_DA_ENS_NUM*num_synop_obs, MPI_REAL8, p_address_master, p_comm_glb, p_err)
         ENDIF
      ENDIF
#endif

!#############################################################################
! Run data assimilation (for only DA)
!#############################################################################
      has_DA = .false.
      IF (has_smap_obs .or. has_fy3d_obs .or. has_synop_obs) THEN
         ! for grid data assimilation
         IF (p_is_worker) THEN
            DO np = 1, numpatch

!#############################################################################
! 1. Find observations around each patch
!#############################################################################
               ! regions info around target patch
               lat_p_n = patchlatr(np)*180/pi + dres
               lat_p_s = patchlatr(np)*180/pi - dres
               lon_p_w = patchlonr(np)*180/pi - dres
               lon_p_e = patchlonr(np)*180/pi + dres

               ! find observations around each patch
               num_obs_p = 0
               num_smap_p = 0
               num_fy3d_p = 0
               num_synop_p = 0
               IF (DEF_DA_SM_SMAP) THEN
                  IF (has_smap_obs) THEN
                     num_smap_p = count(smap_lat(:) < lat_p_n .and. smap_lat(:) > lat_p_s .and. &
                           smap_lon(:) > lon_p_w .and. smap_lon(:) < lon_p_e .and. &
                           smap_time(:) - idate_b >= 0 .and. smap_time(:) - idate_e <= 0 .and. &
                           smap_tb_h_mean(:) > 0 .and. pred_smap_tb_h_mean(:) > 0)
                  ENDIF
               ENDIF
               IF (DEF_DA_SM_FY) THEN
                  IF (has_fy3d_obs) THEN
                     num_fy3d_p = count(fy3d_lat(:) < lat_p_n .and. fy3d_lat(:) > lat_p_s .and. &
                           fy3d_lon(:) > lon_p_w .and. fy3d_lon(:) < lon_p_e .and. &
                           fy3d_time(:) - idate_b >= 0 .and. fy3d_time(:) - idate_e <= 0 .and. &
                           fy3d_tb_h_mean(:) > 0 .and. pred_fy3d_tb_h_mean(:) > 0)
                  ENDIF
               ENDIF
               IF (DEF_DA_SM_SYNOP) THEN
                  IF (has_synop_obs) THEN
                     num_synop_p = count( &
                        synop_lat(:) < lat_p_n .and. synop_lat(:) > lat_p_s .and. &
                        synop_lon(:) > lon_p_w .and. synop_lon(:) < lon_p_e .and. &
                        synop_time(:) - idate_b >= 0 .and. synop_time(:) - idate_e <= 0 .and. &
                        tref_ens_o(:,1) > 0)
                  ENDIF
               ENDIF
               num_obs_p = num_smap_p + num_fy3d_p + 2*num_synop_p

!#############################################################################
! 2. Perform data assimilation
!#############################################################################
               IF (num_obs_p > 0) THEN
                  ! allocate memory for data assimilation
                  allocate (obs_err        (num_obs_p                    ))
                  allocate (trans          (DEF_DA_ENS_NUM,DEF_DA_ENS_NUM))
                  allocate (wliq_soi_ens   (nl_soil,   DEF_DA_ENS_NUM    ))
                  allocate (wliq_soi_ens_da(nl_soil,   DEF_DA_ENS_NUM    ))
                  allocate (t_soi_ens      (nl_soil,   DEF_DA_ENS_NUM    ))
                  allocate (t_soi_ens_da   (nl_soil,   DEF_DA_ENS_NUM    ))
                  allocate (pred_obs_p_ens (num_obs_p, DEF_DA_ENS_NUM    ))
                  allocate (obs_p          (num_obs_p                    ))


!#############################################################################
! 2.1. Prepare observations around target patch for each observation dataset
!#############################################################################
                  IF (num_smap_p > 0) THEN
                     ! allocate memory
                     allocate (index_smap_p         (num_smap_p            ))
                     allocate (smap_lat_p           (num_smap_p            ))
                     allocate (smap_lon_p           (num_smap_p            ))
                     allocate (smap_tb_h_p          (num_smap_p            ))
                     allocate (smap_tb_v_p          (num_smap_p            ))
                     allocate (pred_smap_tb_h_p_ens (num_smap_p, DEF_DA_ENS_NUM))
                     allocate (pred_smap_tb_v_p_ens (num_smap_p, DEF_DA_ENS_NUM))
                     allocate (smap_tb_h_mean_p     (num_smap_p            ))
                     allocate (smap_tb_v_mean_p     (num_smap_p            ))
                     allocate (pred_smap_tb_h_mean_p(num_smap_p            ))
                     allocate (pred_smap_tb_v_mean_p(num_smap_p            ))
                     allocate (smap_err            (num_smap_p            ))
                     allocate (d_smap_p            (num_smap_p            ))

                     ! index of observations around target patch
                     index_smap_p = (smap_lat(:) < lat_p_n .and. smap_lat(:) > lat_p_s .and. &
                           smap_lon(:) > lon_p_w .and. smap_lon(:) < lon_p_e .and. &
                           smap_time(:) - idate_b >= 0 .and. smap_time(:) - idate_e <= 0 .and. &
                           smap_tb_h_mean(:) > 0 .and. pred_smap_tb_h_mean(:) > 0)

                     ! crop observations around target patch
                     smap_lat_p  = pack(smap_lat , index_smap_p)
                     smap_lon_p  = pack(smap_lon , index_smap_p)
                     smap_tb_h_p = pack(smap_tb_h, index_smap_p)
                     smap_tb_v_p = pack(smap_tb_v, index_smap_p)

                     ! crop ensemble predicted observations around target patch
                     DO i = 1, DEF_DA_ENS_NUM
                        pred_smap_tb_h_p_ens(:, i) = pack(pred_smap_tb_h_ogrid_ens(:, i), index_smap_p)
                        pred_smap_tb_v_p_ens(:, i) = pack(pred_smap_tb_v_ogrid_ens(:, i), index_smap_p)
                     ENDDO

                     ! crop mean values around target patch
                     smap_tb_h_mean_p  = pack(smap_tb_h_mean, index_smap_p)
                     smap_tb_v_mean_p  = pack(smap_tb_v_mean, index_smap_p)
                     pred_smap_tb_h_mean_p = pack(pred_smap_tb_h_mean, index_smap_p)
                     pred_smap_tb_v_mean_p = pack(pred_smap_tb_v_mean, index_smap_p)

                     ! scaling
                     smap_tb_h_p = smap_tb_h_p - smap_tb_h_mean_p
                     smap_tb_v_p = smap_tb_v_p - smap_tb_v_mean_p
                     DO i = 1, DEF_DA_ENS_NUM
                        pred_smap_tb_h_p_ens(:, i) = pred_smap_tb_h_p_ens(:, i) - pred_smap_tb_h_mean_p
                        pred_smap_tb_v_p_ens(:, i) = pred_smap_tb_v_p_ens(:, i) - pred_smap_tb_v_mean_p
                     ENDDO

                     ! calculate distance between observation and patch center
                     d_smap_p = 2*6.3781e3*asin(sqrt(sin((smap_lat_p*pi/180 - patchlatr(np))/2.0)**2 + &
                        cos(smap_lat_p*pi/180)*cos(patchlatr(np))*sin((smap_lon_p*pi/180 - patchlonr(np))/2.0)**2))

                     ! calculate weighted observation error
                     smap_err = static_smap_err/(exp((-d_smap_p**2)/(2*loc_r**2)))
                  ENDIF

                  IF (num_fy3d_p > 0) THEN
                     ! allocate memory
                     allocate (index_fy3d_p         (num_fy3d_p            ))
                     allocate (fy3d_lat_p           (num_fy3d_p            ))
                     allocate (fy3d_lon_p           (num_fy3d_p            ))
                     allocate (fy3d_tb_h_p          (num_fy3d_p            ))
                     allocate (fy3d_tb_v_p          (num_fy3d_p            ))
                     allocate (pred_fy3d_tb_h_p_ens (num_fy3d_p, DEF_DA_ENS_NUM))
                     allocate (pred_fy3d_tb_v_p_ens (num_fy3d_p, DEF_DA_ENS_NUM))
                     allocate (fy3d_tb_h_mean_p     (num_fy3d_p            ))
                     allocate (fy3d_tb_v_mean_p     (num_fy3d_p            ))
                     allocate (pred_fy3d_tb_h_mean_p(num_fy3d_p            ))
                     allocate (pred_fy3d_tb_v_mean_p(num_fy3d_p            ))
                     allocate (fy3d_err            (num_fy3d_p             ))
                     allocate (d_fy3d_p            (num_fy3d_p             ))

                     ! index of observations around target patch
                     index_fy3d_p = (fy3d_lat(:) < lat_p_n .and. fy3d_lat(:) > lat_p_s .and. &
                           fy3d_lon(:) > lon_p_w .and. fy3d_lon(:) < lon_p_e .and. &
                           fy3d_time(:) - idate_b >= 0 .and. fy3d_time(:) - idate_e <= 0 .and. &
                           fy3d_tb_h_mean(:) > 0 .and. pred_fy3d_tb_h_mean(:) > 0)

                     ! crop observations around target patch
                     fy3d_lat_p  = pack(fy3d_lat , index_fy3d_p)
                     fy3d_lon_p  = pack(fy3d_lon , index_fy3d_p)
                     fy3d_tb_h_p = pack(fy3d_tb_h, index_fy3d_p)
                     fy3d_tb_v_p = pack(fy3d_tb_v, index_fy3d_p)

                     ! crop ensemble predicted observations around target patch
                     DO i = 1, DEF_DA_ENS_NUM
                        pred_fy3d_tb_h_p_ens(:, i) = pack(pred_fy3d_tb_h_ogrid_ens(:, i), index_fy3d_p)
                        pred_fy3d_tb_v_p_ens(:, i) = pack(pred_fy3d_tb_v_ogrid_ens(:, i), index_fy3d_p)
                     ENDDO

                     ! crop mean values around target patch
                     fy3d_tb_h_mean_p  = pack(fy3d_tb_h_mean, index_fy3d_p)
                     fy3d_tb_v_mean_p  = pack(fy3d_tb_v_mean, index_fy3d_p)
                     pred_fy3d_tb_h_mean_p = pack(pred_fy3d_tb_h_mean, index_fy3d_p)
                     pred_fy3d_tb_v_mean_p = pack(pred_fy3d_tb_v_mean, index_fy3d_p)

                     ! scaling
                     fy3d_tb_h_p = fy3d_tb_h_p - fy3d_tb_h_mean_p
                     fy3d_tb_v_p = fy3d_tb_v_p - fy3d_tb_v_mean_p
                     DO i = 1, DEF_DA_ENS_NUM
                        pred_fy3d_tb_h_p_ens(:, i) = pred_fy3d_tb_h_p_ens(:, i) - pred_fy3d_tb_h_mean_p
                        pred_fy3d_tb_v_p_ens(:, i) = pred_fy3d_tb_v_p_ens(:, i) - pred_fy3d_tb_v_mean_p
                     ENDDO

                     ! calculate distance between observation and patch center
                     d_fy3d_p = 2*6.3781e3*asin(sqrt(sin((fy3d_lat_p*pi/180 - patchlatr(np))/2.0)**2 + &
                        cos(fy3d_lat_p*pi/180)*cos(patchlatr(np))*sin((fy3d_lon_p*pi/180 - patchlonr(np))/2.0)**2))

                     ! calculate weighted observation error
                     fy3d_err = static_fy3d_err/(exp((-d_fy3d_p**2)/(2*loc_r**2)))
                  ENDIF

                  IF (num_synop_p > 0) THEN
                     ! allocate memory
                     allocate (index_synop_p        (num_synop_p             ))
                     allocate (synop_lat_p          (num_synop_p             ))
                     allocate (synop_lon_p          (num_synop_p             ))
                     allocate (synop_qref_p         (num_synop_p             ))
                     allocate (synop_tref_p         (num_synop_p             ))
                     allocate (qref_ens_p           (num_synop_p, DEF_DA_ENS_NUM ))
                     allocate (tref_ens_p           (num_synop_p, DEF_DA_ENS_NUM ))
                     allocate (synop_tref_err       (num_synop_p             ))
                     allocate (synop_qref_err       (num_synop_p             ))
                     allocate (d_synop_p           (num_synop_p             ))

                     ! index of observations around target patch
                     index_synop_p = ( &
                        synop_lat(:) < lat_p_n .and. synop_lat(:) > lat_p_s .and. &
                        synop_lon(:) > lon_p_w .and. synop_lon(:) < lon_p_e .and. &
                        synop_time(:) - idate_b >= 0 .and. synop_time(:) - idate_e <= 0 .and. &
                        tref_ens_o(:,1) > 0)

                     ! crop observations around target patch
                     synop_lat_p  = pack(synop_lat, index_synop_p)
                     synop_lon_p  = pack(synop_lon, index_synop_p)
                     synop_qref_p = pack(synop_qref, index_synop_p)
                     synop_tref_p = pack(synop_tref, index_synop_p)

                     ! crop ensemble predicted observations around target patch
                     DO i = 1, DEF_DA_ENS_NUM
                       qref_ens_p(:, i) = pack(qref_ens_o(:, i), index_synop_p)
                       tref_ens_p(:, i) = pack(tref_ens_o(:, i), index_synop_p)
                     ENDDO

                     ! calculate distance between observation and patch center
                     d_synop_p = 2*6.3781e3*asin(sqrt(sin((synop_lat_p*pi/180 - patchlatr(np))/2.0)**2 + &
                        cos(synop_lat_p*pi/180)*cos(patchlatr(np))*sin((synop_lon_p*pi/180 - patchlonr(np))/2.0)**2))

                     ! calculate weighted observation error
                     synop_tref_err = static_synop_tref_err/(exp((-d_synop_p**2)/(2*loc_r**2)))
                     synop_qref_err = static_synop_qref_err/(exp((-d_synop_p**2)/(2*loc_r**2)))
                  ENDIF

!#############################################################################
! 2.2. Concat observations & ensemble predictions around target patch
!#############################################################################
                  ! index of end of each observation dataset
                  end_idx(1) = num_smap_p
                  end_idx(2) = num_smap_p + num_fy3d_p
                  end_idx(3) = num_smap_p + num_fy3d_p + num_synop_p

                  ! concatenate ensemble predictions
                  IF (num_smap_p > 0) THEN
                     pred_obs_p_ens(1:end_idx(1),:) = pred_smap_tb_h_p_ens
                     obs_p(1:end_idx(1)) = smap_tb_h_p
                     obs_err(1:end_idx(1)) = smap_err
                  ENDIF
                  IF (num_fy3d_p > 0) THEN
                     pred_obs_p_ens(end_idx(1)+1:end_idx(2),:) = pred_fy3d_tb_h_p_ens
                     obs_p(end_idx(1)+1:end_idx(2)) = fy3d_tb_h_p
                     obs_err(end_idx(1)+1:end_idx(2)) = fy3d_err
                  ENDIF
                  IF (num_synop_p > 0) THEN
                     pred_obs_p_ens(end_idx(2)+1:end_idx(3),:) = tref_ens_p
                     pred_obs_p_ens(end_idx(3)+1:,:) = qref_ens_p
                     obs_p(end_idx(2)+1:end_idx(3)) = synop_tref_p
                     obs_p(end_idx(3)+1:) = synop_qref_p
                     obs_err(end_idx(2)+1:end_idx(3)) = synop_tref_err
                     obs_err(end_idx(3)+1:) = synop_qref_err
                  ENDIF
                  !index_nearest = minloc(d_p, dim=1)

!#############################################################################
! Perform LETKF and postprocess analysis results
!#############################################################################
                  ! calculate transformation matrix
                  CALL letkf(DEF_DA_ENS_NUM, num_obs_p, &
                       pred_obs_p_ens, obs_p, obs_err, infl, &
                       trans)

                  ! calculate analysis value
                  IF (patchtype(np) < 3) THEN ! ocean, lake, ice
                     has_DA = .true.

                     ! soil layer
                     DO iens = 1, DEF_DA_ENS_NUM
                        wliq_soi_ens(:, iens) = wliq_soisno_ens(1:, iens, np)
                        t_soi_ens(:, iens)    = t_soisno_ens(1:, iens, np)
                     ENDDO

                     ! analysis
                     CALL dgemm('N', 'N', nl_soil, DEF_DA_ENS_NUM, DEF_DA_ENS_NUM, 1.0_8, wliq_soi_ens, &
                           nl_soil, trans, DEF_DA_ENS_NUM, 0.0_8, wliq_soi_ens_da, nl_soil)
                     CALL dgemm('N', 'N', nl_soil, DEF_DA_ENS_NUM, DEF_DA_ENS_NUM, 1.0_8, t_soi_ens, &
                           nl_soil, trans, DEF_DA_ENS_NUM, 0.0_8, t_soi_ens_da, nl_soil)

                     ! save analysis results
                     wliq_soi_ens_da = max(0.0, wliq_soi_ens_da)
                     DO iens = 1, DEF_DA_ENS_NUM
                        wliq_soisno_ens(1:2, iens, np) = wliq_soi_ens_da(1:2, iens)
                        ! limit the change of t_soisno in range of [-20.0, 20.0]
                        WHERE (t_soi_ens_da(1:2,iens) > t_soisno_ens(1:2,iens,np))
                           t_soisno_ens(1:2,iens,np) = t_soisno_ens(1:2,iens,np) + min(20.0, t_soi_ens_da(1:2,iens)-t_soisno_ens(1:2,iens,np))
                        ELSEWHERE
                           t_soisno_ens(1:2,iens,np) = t_soisno_ens(1:2,iens,np) + max(-20.0, t_soi_ens_da(1:2,iens)-t_soisno_ens(1:2,iens,np))
                        ENDWHERE
                     ENDDO

                     ! limit the soil liquid and ice water in a reasonable range
                     DO i = 1, nl_soil
                        DO iens = 1, DEF_DA_ENS_NUM
                           ! lower bound
                           wliq_soisno_ens(i, iens, np) = max(0.0d0, wliq_soisno_ens(i, iens, np))
                           wice_soisno_ens(i, iens, np) = max(0.0d0, wice_soisno_ens(i, iens, np))
                           IF (wliq_soisno_ens(i, iens, np) == 0.0 .and. wice_soisno_ens(i, iens, np) == 0.0) THEN
                              IF (t_soisno_ens(i, iens, np)-tfrz < -5.0) THEN
                                 wice_soisno_ens(i, iens, np) = 1e-10
                              ELSE
                                 wliq_soisno_ens(i, iens, np) = 1e-10
                              ENDIF
                           ENDIF

                           ! upper bound
                           wliq_soisno_ens(i, iens, np) = min(porsl(i, np)*(dz_soi(i)*denh2o), wliq_soisno_ens(i, iens, np))
                           eff_porsl = max(0.0d0, porsl(i, np) - wliq_soisno_ens(i, iens, np)/(dz_soi(i)*denh2o))
                           wice_soisno_ens(i, iens, np) = min(eff_porsl*(dz_soi(i)*denice), wice_soisno_ens(i, iens, np))
                        ENDDO
                     ENDDO

                     ! move residual water to water table
                     wa_ens(:, np) = wa_ens(:, np) - sum(wliq_soisno_ens(1:, :, np) - wliq_soi_ens, dim=1)

                     ! update volumetric water content for diagnostic
                     DO iens = 1, DEF_DA_ENS_NUM
                        h2osoi_ens(:, iens, np) = wliq_soisno_ens(1:, iens, np)/(dz_soi(:)*denh2o) + wice_soisno_ens(1:, iens, np)/(dz_soi(:)*denice)
                        h2osoi_ens(:, iens, np) = min(1.0d0, h2osoi_ens(:, iens, np))
                        h2osoi_ens(:, iens, np) = max(0.0d0, h2osoi_ens(:, iens, np))
                     ENDDO
                  ENDIF
               ENDIF


!#############################################################################
! deallocate memory changes with patch
!#############################################################################
               IF (allocated(trans                )) deallocate (trans                )
               IF (allocated(smap_err             )) deallocate (smap_err             )
               IF (allocated(d_smap_p             )) deallocate (d_smap_p             )
               IF (allocated(fy3d_err             )) deallocate (fy3d_err             )
               IF (allocated(d_fy3d_p             )) deallocate (d_fy3d_p             )
               IF (allocated(synop_tref_err       )) deallocate (synop_tref_err       )
               IF (allocated(synop_qref_err       )) deallocate (synop_qref_err       )
               IF (allocated(d_synop_p            )) deallocate (d_synop_p            )
               IF (allocated(pred_obs_p_ens       )) deallocate (pred_obs_p_ens       )
               IF (allocated(obs_p                )) deallocate (obs_p                )
               IF (allocated(obs_err              )) deallocate (obs_err              )
               IF (allocated(wliq_soi_ens         )) deallocate (wliq_soi_ens         )
               IF (allocated(wliq_soi_ens_da      )) deallocate (wliq_soi_ens_da      )
               IF (allocated(t_soi_ens            )) deallocate (t_soi_ens            )
               IF (allocated(t_soi_ens_da         )) deallocate (t_soi_ens_da         )

               IF (allocated(index_smap_p         )) deallocate (index_smap_p         )
               IF (allocated(smap_lat_p           )) deallocate (smap_lat_p           )
               IF (allocated(smap_lon_p           )) deallocate (smap_lon_p           )
               IF (allocated(smap_tb_h_p          )) deallocate (smap_tb_h_p          )
               IF (allocated(smap_tb_v_p          )) deallocate (smap_tb_v_p          )
               IF (allocated(pred_smap_tb_h_p_ens )) deallocate (pred_smap_tb_h_p_ens )
               IF (allocated(pred_smap_tb_v_p_ens )) deallocate (pred_smap_tb_v_p_ens )
               IF (allocated(smap_tb_h_mean_p     )) deallocate (smap_tb_h_mean_p     )
               IF (allocated(pred_smap_tb_h_mean_p)) deallocate (pred_smap_tb_h_mean_p)
               IF (allocated(smap_tb_v_mean_p     )) deallocate (smap_tb_v_mean_p     )
               IF (allocated(pred_smap_tb_v_mean_p)) deallocate (pred_smap_tb_v_mean_p)

               IF (allocated(index_fy3d_p         )) deallocate (index_fy3d_p         )
               IF (allocated(fy3d_lat_p           )) deallocate (fy3d_lat_p           )
               IF (allocated(fy3d_lon_p           )) deallocate (fy3d_lon_p           )
               IF (allocated(fy3d_tb_h_p          )) deallocate (fy3d_tb_h_p          )
               IF (allocated(fy3d_tb_v_p          )) deallocate (fy3d_tb_v_p          )
               IF (allocated(pred_fy3d_tb_h_p_ens )) deallocate (pred_fy3d_tb_h_p_ens )
               IF (allocated(pred_fy3d_tb_v_p_ens )) deallocate (pred_fy3d_tb_v_p_ens )
               IF (allocated(fy3d_tb_h_mean_p     )) deallocate (fy3d_tb_h_mean_p     )
               IF (allocated(pred_fy3d_tb_h_mean_p)) deallocate (pred_fy3d_tb_h_mean_p)
               IF (allocated(fy3d_tb_v_mean_p     )) deallocate (fy3d_tb_v_mean_p     )
               IF (allocated(pred_fy3d_tb_v_mean_p)) deallocate (pred_fy3d_tb_v_mean_p)

               IF (allocated(index_synop_p        )) deallocate (index_synop_p        )
               IF (allocated(synop_lat_p          )) deallocate (synop_lat_p          )
               IF (allocated(synop_lon_p          )) deallocate (synop_lon_p          )
               IF (allocated(synop_qref_p         )) deallocate (synop_qref_p         )
               IF (allocated(synop_tref_p         )) deallocate (synop_tref_p         )
               IF (allocated(qref_ens_p           )) deallocate (qref_ens_p           )
               IF (allocated(tref_ens_p           )) deallocate (tref_ens_p           )
            ENDDO
         ENDIF
      ENDIF


!#############################################################################
! Calculate ensemble brightness temperature after DA for diagnostic
!#############################################################################
      IF (has_DA) THEN
         IF (p_is_worker) THEN
            IF (DEF_DA_SM_SMAP) THEN
               IF (DEF_DA_ENS_NUM > 1) THEN
                  DO iens = 1, DEF_DA_ENS_NUM
                     DO np = 1, numpatch
                        CALL forward(&
                           patchtype(np), patchclass(np), dz_sno_ens(:, iens, np), &
                           forc_topo(np), htop(np), &
                           tref_ens(iens, np), t_soisno_ens(:,iens,np), tleaf_ens(iens, np), &
                           wliq_soisno_ens(:, iens, np), wice_soisno_ens(:, iens, np), h2osoi_ens(:, iens, np), &
                           snowdp_ens(iens, np), lai_ens(iens, np), sai_ens(iens, np), &
                           wf_clay(:, np), wf_sand(:, np), wf_silt(:, np), BD_all(:, np), porsl(:, np), &
                           smap_theta, smap_fghz, &
                           t_brt_smap_ens(1,iens,np), t_brt_smap_ens(2,iens,np))
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF

            IF (DEF_DA_SM_FY) THEN
               IF (DEF_DA_ENS_NUM > 1) THEN
                  DO iens = 1, DEF_DA_ENS_NUM
                     DO np = 1, numpatch
                        CALL forward(&
                           patchtype(np), patchclass(np), dz_sno_ens(:, iens, np), &
                           forc_topo(np), htop(np), &
                           tref_ens(iens, np), t_soisno_ens(:,iens,np), tleaf_ens(iens, np), &
                           wliq_soisno_ens(:, iens, np), wice_soisno_ens(:, iens, np), h2osoi_ens(:, iens, np), &
                           snowdp_ens(iens, np), lai_ens(iens, np), sai_ens(iens, np), &
                           wf_clay(:, np), wf_sand(:, np), wf_silt(:, np), BD_all(:, np), porsl(:, np), &
                           fy3d_theta, fy3d_fghz, &
                           t_brt_fy3d_ens(1,iens,np), t_brt_fy3d_ens(2,iens,np))
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
         ENDIF
      ENDIF

!#############################################################################
! Deallocate all memory
!#############################################################################
      IF (allocated (trans)) deallocate (trans)
      IF (allocated (obs_err)) deallocate (obs_err)
      IF (allocated (filter    )) deallocate (filter    )
      IF (allocated (wliq_soi_ens)) deallocate (wliq_soi_ens)
      IF (allocated (wliq_soi_ens_da)) deallocate (wliq_soi_ens_da)
      IF (allocated (t_soi_ens)) deallocate (t_soi_ens)
      IF (allocated (t_soi_ens_da)) deallocate (t_soi_ens_da)

      IF (allocated (smap_time)) deallocate (smap_time)
      IF (allocated (dt_b_smap)) deallocate (dt_b_smap)
      IF (allocated (dt_e_smap)) deallocate (dt_e_smap)
      IF (allocated (smap_lat)) deallocate (smap_lat)
      IF (allocated (smap_lon)) deallocate (smap_lon)
      IF (allocated (smap_tb_h)) deallocate (smap_tb_h)
      IF (allocated (smap_tb_v)) deallocate (smap_tb_v)
      IF (allocated (smap_ii)) deallocate (smap_ii)
      IF (allocated (smap_jj)) deallocate (smap_jj)
      IF (allocated (smap_tb_h_mean)) deallocate (smap_tb_h_mean)
      IF (allocated (smap_tb_v_mean)) deallocate (smap_tb_v_mean)
      IF (allocated (pred_smap_tb_h_mean)) deallocate (pred_smap_tb_h_mean)
      IF (allocated (pred_smap_tb_v_mean)) deallocate (pred_smap_tb_v_mean)
      IF (allocated (pred_smap_tb_h_pset_ens)) deallocate (pred_smap_tb_h_pset_ens)
      IF (allocated (pred_smap_tb_v_pset_ens)) deallocate (pred_smap_tb_v_pset_ens)
      IF (allocated (pred_smap_tb_h_ogrid_ens)) deallocate (pred_smap_tb_h_ogrid_ens)
      IF (allocated (pred_smap_tb_v_ogrid_ens)) deallocate (pred_smap_tb_v_ogrid_ens)

      IF (allocated (fy3d_time)) deallocate (fy3d_time)
      IF (allocated (dt_b_fy3d)) deallocate (dt_b_fy3d)
      IF (allocated (dt_e_fy3d)) deallocate (dt_e_fy3d)
      IF (allocated (fy3d_lat)) deallocate (fy3d_lat)
      IF (allocated (fy3d_lon)) deallocate (fy3d_lon)
      IF (allocated (fy3d_tb_h)) deallocate (fy3d_tb_h)
      IF (allocated (fy3d_tb_v)) deallocate (fy3d_tb_v)
      IF (allocated (fy3d_ii)) deallocate (fy3d_ii)
      IF (allocated (fy3d_jj)) deallocate (fy3d_jj)
      IF (allocated (fy3d_tb_h_mean)) deallocate (fy3d_tb_h_mean)
      IF (allocated (fy3d_tb_v_mean)) deallocate (fy3d_tb_v_mean)
      IF (allocated (pred_fy3d_tb_h_mean)) deallocate (pred_fy3d_tb_h_mean)
      IF (allocated (pred_fy3d_tb_v_mean)) deallocate (pred_fy3d_tb_v_mean)
      IF (allocated (pred_fy3d_tb_h_pset_ens)) deallocate (pred_fy3d_tb_h_pset_ens)
      IF (allocated (pred_fy3d_tb_v_pset_ens)) deallocate (pred_fy3d_tb_v_pset_ens)
      IF (allocated (pred_fy3d_tb_h_ogrid_ens)) deallocate (pred_fy3d_tb_h_ogrid_ens)
      IF (allocated (pred_fy3d_tb_v_ogrid_ens)) deallocate (pred_fy3d_tb_v_ogrid_ens)

      IF (allocated (synop_time)) deallocate (synop_time)
      IF (allocated (synop_lat)) deallocate (synop_lat)
      IF (allocated (synop_lon)) deallocate (synop_lon)
      IF (allocated (synop_qref)) deallocate (synop_qref)
      IF (allocated (synop_tref)) deallocate (synop_tref)
      IF (allocated (tref_ens_o)) deallocate (tref_ens_o)
      IF (allocated (qref_ens_o)) deallocate (qref_ens_o)

      IF (allocated (synop_idx)) deallocate (synop_idx)
      IF (allocated (site_id_worker)) deallocate (site_id_worker)
      IF (allocated (tref_ens_worker)) deallocate (tref_ens_worker)
      IF (allocated (qref_ens_worker)) deallocate (qref_ens_worker)
      IF (allocated (qref_ens_o)) deallocate (qref_ens_o)
      IF (allocated (tref_ens_o)) deallocate (tref_ens_o)



   END SUBROUTINE run_DA_SM

!-----------------------------------------------------------------------------

   SUBROUTINE end_DA_SM()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!-----------------------------------------------------------------------------
      IF (allocated (trans)) deallocate (trans)
      IF (allocated (obs_err)) deallocate (obs_err)
      IF (allocated (filter    )) deallocate (filter    )
      IF (allocated (wliq_soi_ens)) deallocate (wliq_soi_ens)
      IF (allocated (wliq_soi_ens_da)) deallocate (wliq_soi_ens_da)
      IF (allocated (t_soi_ens)) deallocate (t_soi_ens)
      IF (allocated (t_soi_ens_da)) deallocate (t_soi_ens_da)

      IF (allocated (smap_time)) deallocate (smap_time)
      IF (allocated (dt_b_smap)) deallocate (dt_b_smap)
      IF (allocated (dt_e_smap)) deallocate (dt_e_smap)
      IF (allocated (smap_lat)) deallocate (smap_lat)
      IF (allocated (smap_lon)) deallocate (smap_lon)
      IF (allocated (smap_tb_h)) deallocate (smap_tb_h)
      IF (allocated (smap_tb_v)) deallocate (smap_tb_v)
      IF (allocated (smap_ii)) deallocate (smap_ii)
      IF (allocated (smap_jj)) deallocate (smap_jj)
      IF (allocated (smap_tb_h_mean)) deallocate (smap_tb_h_mean)
      IF (allocated (smap_tb_v_mean)) deallocate (smap_tb_v_mean)
      IF (allocated (pred_smap_tb_h_mean)) deallocate (pred_smap_tb_h_mean)
      IF (allocated (pred_smap_tb_v_mean)) deallocate (pred_smap_tb_v_mean)
      IF (allocated (pred_smap_tb_h_pset_ens)) deallocate (pred_smap_tb_h_pset_ens)
      IF (allocated (pred_smap_tb_v_pset_ens)) deallocate (pred_smap_tb_v_pset_ens)
      IF (allocated (pred_smap_tb_h_ogrid_ens)) deallocate (pred_smap_tb_h_ogrid_ens)
      IF (allocated (pred_smap_tb_v_ogrid_ens)) deallocate (pred_smap_tb_v_ogrid_ens)

      IF (allocated (fy3d_time)) deallocate (fy3d_time)
      IF (allocated (dt_b_fy3d)) deallocate (dt_b_fy3d)
      IF (allocated (dt_e_fy3d)) deallocate (dt_e_fy3d)
      IF (allocated (fy3d_lat)) deallocate (fy3d_lat)
      IF (allocated (fy3d_lon)) deallocate (fy3d_lon)
      IF (allocated (fy3d_tb_h)) deallocate (fy3d_tb_h)
      IF (allocated (fy3d_tb_v)) deallocate (fy3d_tb_v)
      IF (allocated (fy3d_ii)) deallocate (fy3d_ii)
      IF (allocated (fy3d_jj)) deallocate (fy3d_jj)
      IF (allocated (fy3d_tb_h_mean)) deallocate (fy3d_tb_h_mean)
      IF (allocated (fy3d_tb_v_mean)) deallocate (fy3d_tb_v_mean)
      IF (allocated (pred_fy3d_tb_h_mean)) deallocate (pred_fy3d_tb_h_mean)
      IF (allocated (pred_fy3d_tb_v_mean)) deallocate (pred_fy3d_tb_v_mean)
      IF (allocated (pred_fy3d_tb_h_pset_ens)) deallocate (pred_fy3d_tb_h_pset_ens)
      IF (allocated (pred_fy3d_tb_v_pset_ens)) deallocate (pred_fy3d_tb_v_pset_ens)
      IF (allocated (pred_fy3d_tb_h_ogrid_ens)) deallocate (pred_fy3d_tb_h_ogrid_ens)
      IF (allocated (pred_fy3d_tb_v_ogrid_ens)) deallocate (pred_fy3d_tb_v_ogrid_ens)

      IF (allocated (synop_time)) deallocate (synop_time)
      IF (allocated (synop_lat)) deallocate (synop_lat)
      IF (allocated (synop_lon)) deallocate (synop_lon)
      IF (allocated (synop_qref)) deallocate (synop_qref)
      IF (allocated (synop_tref)) deallocate (synop_tref)
      IF (allocated (tref_ens_o)) deallocate (tref_ens_o)
      IF (allocated (qref_ens_o)) deallocate (qref_ens_o)

      IF (allocated (synop_idx)) deallocate (synop_idx)
      IF (allocated (site_id_worker)) deallocate (site_id_worker)
      IF (allocated (tref_ens_worker)) deallocate (tref_ens_worker)
      IF (allocated (qref_ens_worker)) deallocate (qref_ens_worker)
      IF (allocated (qref_ens_o)) deallocate (qref_ens_o)
      IF (allocated (tref_ens_o)) deallocate (tref_ens_o)


   END SUBROUTINE end_DA_SM

!-----------------------------------------------------------------------------
END MODULE MOD_DA_SM
#endif
