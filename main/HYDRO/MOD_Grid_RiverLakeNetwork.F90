#include <define.h>

#ifdef GridRiverLakeFlow
#ifndef USEMPI
! The serial path never builds the river-system labels (rivermouth) used below.
#error "GridRiverLakeFlow requires USEMPI (see include/define.h)"
#endif
MODULE MOD_Grid_RiverLakeNetwork
!--------------------------------------------------------------------------------
! DESCRIPTION:
!--------------------------------------------------------------------------------

   USE MOD_Grid
   USE MOD_WorkerPushData
   USE MOD_Grid_RiverLakePush
   USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
#ifdef USEMPI
   USE MOD_SPMD_Task, ONLY: MPI_COMM_NULL
#endif
   IMPLICIT NONE

   ! ----- River Lake network -----

   type(grid_type) :: griducat

   integer :: totalnumucat
   integer :: numucat
   integer, allocatable :: ucat_ucid (:)   ! index in unit catchment numbering
   integer, allocatable :: x_ucat    (:)   !
   integer, allocatable :: y_ucat    (:)   !
   integer, allocatable :: ucat_gdid (:)   !

   integer, allocatable :: numucat_wrk (:)
   type(pointer_int32_1d), allocatable :: ucat_data_address (:)

   ! ----- Part 1: between runoff input elements and unit catchments -----
   integer :: numinpm
   integer,  allocatable :: inpm_gdid (:)

   integer :: inpn
   integer,  allocatable :: idmap_gd2uc (:,:)
   real(r8), allocatable :: area_gd2uc  (:,:)

   integer :: nucpart
   integer,  allocatable :: idmap_uc2gd (:,:)
   real(r8), allocatable :: area_uc2gd  (:,:)

   type(worker_remapdata_type) :: remap_patch2inpm
   type(worker_pushdata_type)  :: push_inpm2ucat
   type(worker_pushdata_type)  :: push_ucat2inpm
   type(worker_pushdata_type)  :: push_ucat2grid
   type(worker_pushdata_type)  :: allreduce_inpm

   ! ----- Part 2: between upstream and downstream unit catchments -----
   integer,  allocatable :: ucat_next (:)  ! next unit catchment
   integer :: upnmax
   integer,  allocatable :: ucat_ups (:,:) ! upstream unit catchments
   real(r8), allocatable :: wts_ups  (:,:)

   type(worker_pushdata_type) :: push_next2ucat
   type(worker_pushdata_type) :: push_ups2ucat
   type(riverlake_pushcache_type) :: cache_ups2ucat   ! batched pushes (MOD_Grid_RiverLakePush)

   ! ----- Part 3: river systems -----
   integer :: numrivsys
   logical :: rivsys_by_multiple_procs
   integer, allocatable :: irivsys (:)
#ifdef USEMPI
   integer :: p_comm_rivsys = MPI_COMM_NULL
#endif

   ! Adaptive time step groups under bifurcation.
   ! Reason: with DEF_USE_BIFURCATION the routing sub-step used to be one global
   ! minimum dt over all river systems of all workers, so the whole globe
   ! advanced with the dt of its single most restrictive cell (about 1.5 times
   ! the unit-catchment sub-steps of the per-river-system dt, and about twice the
   ! wall time of a run without bifurcation).  A global dt is only needed where
   ! bifurcation pathways couple two river systems: the flux over a pathway must
   ! use the same dt on both sides to conserve water.  Here, in bifurcation mode,
   ! "river system" (irivsys) means a group of river systems that are connected
   ! by bifurcation pathways (or a single unconnected river system); groups
   ! advance with their own dt like river systems without bifurcation.
   ! Effect on results: this changes the time discretisation (each group takes
   ! its own adaptive sub-steps instead of the global minimum), so results differ
   ! at the level of the time-stepping error, as they would between two valid
   ! adaptive time steps.  Without bifurcation nothing changes.
   ! irivsys_shared(k) > 0 marks a group k that lies on more than one worker; its
   ! value indexes the vector that is reduced (MPI_MIN) over all workers so that
   ! every part of that group uses the same dt.  numrivsys_shared is the number of
   ! such groups, identical on all processes.
   integer :: numrivsys_shared = 0
   integer, allocatable :: irivsys_shared (:)


   ! ----- Parameters for River and Lake -----

   integer,  allocatable :: lake_type      (:)   ! 0: river; 2: reservoir.

   real(r8), allocatable :: topo_rivelv    (:)   ! river bed elevation [m]
   real(r8), allocatable :: topo_rivhgt    (:)   ! river channel depth [m]
   real(r8), allocatable :: topo_rivlen    (:)   ! river channel length [m]
   real(r8), allocatable :: topo_rivman    (:)   ! river manning coefficient [m]
   real(r8), allocatable :: topo_rivwth    (:)   ! river channel width [m]
   real(r8), allocatable :: topo_rivare    (:)   ! river channel area [m^2]
   real(r8), allocatable :: topo_rivstomax (:)   ! max river channel storage [m^3]

   real(r8), allocatable :: topo_area      (:)   ! floodplain area [m^2]
   real(r8), allocatable :: topo_fldhgt    (:,:) ! floodplain height profile [m]

   real(r8), allocatable :: levee_frc_data (:)
   real(r8), allocatable :: levee_hgt_data (:)

   integer  :: totalnpthout  = 0
   integer  :: npthout_local = 0
   integer  :: npthlev_bif   = 0

   integer,  allocatable :: pth_upst_local  (:)
   integer,  allocatable :: pth_down_local  (:)
   integer,  allocatable :: pth_down_ucid   (:)
   integer,  allocatable :: pth_global_id   (:)
   real(r8), allocatable :: pth_dst         (:)
   real(r8), allocatable :: pth_elv         (:,:)
   real(r8), allocatable :: pth_wth         (:,:)
   real(r8), allocatable :: pth_man         (:)

   integer  :: max_bif_incoming
   integer,  allocatable :: bif_incoming_pths (:,:)
   real(r8), allocatable :: bif_incoming_wts  (:,:)

   type(worker_pushdata_type) :: push_bif_dn2pth
   type(worker_pushdata_type) :: push_bif_influx
   type(riverlake_pushcache_type) :: cache_bif_dn2pth ! batched pushes (MOD_Grid_RiverLakePush)
   type(riverlake_pushcache_type) :: cache_bif_influx

   real(r8), allocatable :: bedelv_next    (:)   ! downstream river bed elevation [m]
   real(r8), allocatable :: outletwth      (:)   ! river outlet width [m]

   type :: vol_dep_curve_type
      integer  :: nlfp
      real(r8) :: rivhgt
      real(r8) :: rivare
      real(r8) :: rivstomax
      real(r8), allocatable :: flphgt    (:) ! floodplain height profile [m]
      real(r8), allocatable :: flparea   (:) ! flood plain area [m^2]
      real(r8), allocatable :: flpaccare (:) ! flood plain accumulated area [m^2]
      real(r8), allocatable :: flpstomax (:) ! max flood plain storage [m^3]
   CONTAINS
      procedure, PUBLIC :: depth     => retrieve_depth_from_volume
      procedure, PUBLIC :: volume    => retrieve_volume_from_depth
      procedure, PUBLIC :: floodarea => retrieve_area_from_depth
      final :: vol_depth_curve_free_mem
   END type vol_dep_curve_type

   type(vol_dep_curve_type), allocatable :: floodplain_curve (:)


   ! ----- Mask of Grids with all upstream area in the simulation region -----
   real(r8), allocatable :: allups_mask_ucat (:)

CONTAINS

   ! ----------
   SUBROUTINE build_riverlake_network ()

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_Mesh
   USE MOD_Utils
   USE MOD_LandPatch
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   ! Local Variables
   character(len=256)    :: parafile

   integer,  allocatable :: idmap_x(:,:), idmap_y(:,:)

   integer :: numrivmth
   integer,  allocatable :: rivermouth(:)
   integer,  allocatable :: shared_ucat(:)

   integer,  allocatable :: nups_nst  (:), iups_nst  (:), nups_all(:)
   integer,  allocatable :: uc_up2down(:), order_ucat(:)
   integer,  allocatable :: addr_ucat (:)

   integer , allocatable :: nuc_rs(:), iwrk_rs(:), nwrk_rs(:), nave_rs(:)
   real(r8), allocatable :: wt_uc (:), wt_rs  (:), wt_wrk (:), nuc_wrk(:)

   integer,  allocatable :: grdindex(:)


   integer,  allocatable :: idata1d(:), idata2d(:,:)
   real(r8), allocatable :: rdata1d(:), rdata2d(:,:)

   integer,  allocatable :: allgrd_in_inp (:), nucat_g2d(:,:), iucat_g(:)

   integer,  allocatable :: idmap_uc2gd_all(:,:)
   real(r8), allocatable :: area_uc2gd_all (:,:)

   real(r8), allocatable :: ucat_area_all (:)

   integer  :: nlat_ucat, nlon_ucat
   integer  :: nucat, iriv, ngrdall, igrd, ngrd
   integer  :: p_np_rivsys, color
   integer  :: iworker, iwrkdsp
   integer  :: iloc, i, j, ithis
   real(r8) :: sumwt
   logical  :: is_new


#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      ! read in parameters from file.
      IF (p_is_master) THEN

         parafile = get_unitcatchment_file ()

         CALL ncio_read_serial (parafile, 'seq_x', x_ucat)
         CALL ncio_read_serial (parafile, 'seq_y', y_ucat)

         CALL ncio_read_serial (parafile, 'seq_next', ucat_next)

         IF (DEF_UnitCatchment_regional) CALL verify_regional_network (parafile, x_ucat, y_ucat)

         CALL ncio_inquire_length (parafile, 'lon', nlon_ucat)
         CALL ncio_inquire_length (parafile, 'lat', nlat_ucat)

         CALL ncio_read_serial (parafile, 'inpmat_x', idmap_x)
         CALL ncio_read_serial (parafile, 'inpmat_y', idmap_y)
         CALL ncio_read_serial (parafile, 'inpmat_area', area_gd2uc)

      ENDIF

      IF (p_is_master) THEN

         totalnumucat = size(x_ucat)

         allocate (nups_nst (totalnumucat))
         allocate (iups_nst (totalnumucat))

         nups_nst(:) = 0
         DO i = 1, totalnumucat
            j = ucat_next(i)
            IF (j > 0) THEN
               nups_nst(j) = nups_nst(j) + 1
            ENDIF
         ENDDO

         ! sort unit catchment from upstream to downstream, recorded by "uc_up2down"
         allocate (uc_up2down (totalnumucat))

         ithis = 0
         iups_nst(:) = 0
         DO i = 1, totalnumucat
            IF (iups_nst(i) == nups_nst(i)) THEN

               ithis = ithis + 1
               uc_up2down(ithis) = i
               iups_nst(i) = -1

               j = ucat_next(i)
               DO WHILE (j > 0)

                  iups_nst(j) = iups_nst(j) + 1

                  IF (iups_nst(j) == nups_nst(j)) THEN
                     ithis = ithis + 1
                     uc_up2down(ithis) = j
                     iups_nst(j) = -1

                     j = ucat_next(j)
                  ELSE
                     EXIT
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

      ENDIF

#ifdef USEMPI
      ! divide unit catchments into groups and assign to workers
      IF (p_is_master) THEN

         allocate (wt_uc (totalnumucat));  wt_uc(:) = 1.

         allocate (rivermouth (totalnumucat))
         numrivmth = 0
         DO i = totalnumucat, 1, -1
            j = ucat_next(uc_up2down(i))
            IF (j <= 0) THEN
               numrivmth = numrivmth + 1
               rivermouth(uc_up2down(i)) = numrivmth
            ELSE
               rivermouth(uc_up2down(i)) = rivermouth(j)
            ENDIF
         ENDDO

         allocate (nuc_rs (numrivmth)); nuc_rs(:) = 0
         allocate (wt_rs  (numrivmth)); wt_rs (:) = 0.
         DO i = 1, totalnumucat
            nuc_rs(rivermouth(i)) = nuc_rs(rivermouth(i)) + 1
            wt_rs (rivermouth(i)) = wt_rs (rivermouth(i)) + wt_uc(i)
         ENDDO

         sumwt = sum(wt_rs)

         allocate (iwrk_rs (numrivmth))
         allocate (nwrk_rs (numrivmth))
         allocate (nave_rs (numrivmth))

         iwrkdsp = -1
         DO i = 1, numrivmth
            nwrk_rs(i) = floor(wt_rs(i)/sumwt * p_np_worker)
            IF (nwrk_rs(i) > 1) THEN

               nave_rs(i) = nuc_rs(i) / nwrk_rs(i)
               IF (mod(nuc_rs(i), nwrk_rs(i)) /= 0) THEN
                  nave_rs(i) = nave_rs(i) + 1
               ENDIF

               iwrk_rs(i) = iwrkdsp + 1
               iwrkdsp = iwrkdsp + nwrk_rs(i)
            ENDIF
         ENDDO

         allocate (nups_all (totalnumucat));  nups_all(:) = 1

         DO i = 1, totalnumucat
            j = ucat_next(uc_up2down(i))
            IF (j > 0) THEN
               nups_all(j) = nups_all(j) + nups_all(uc_up2down(i))
            ENDIF
         ENDDO

         allocate (addr_ucat (totalnumucat));  addr_ucat(:) = -1

         allocate (wt_wrk (0:p_np_worker-1));  wt_wrk (:) = 0
         allocate (nuc_wrk(0:p_np_worker-1));  nuc_wrk(:) = 0

         allocate (order_ucat (totalnumucat))
         order_ucat(uc_up2down) = (/(i, i = 1, totalnumucat)/)

         ithis = totalnumucat
         DO WHILE (ithis > 0)

            i = uc_up2down(ithis)

            IF (addr_ucat(i) >= 0) THEN
               ithis = ithis - 1
               CYCLE
            ENDIF

            j = ucat_next(i)
            IF (j > 0) THEN
               IF (addr_ucat(j) >= 0) THEN
                  addr_ucat(i) = addr_ucat(j)
                  ithis = ithis - 1
                  CYCLE
               ENDIF
            ENDIF

            iriv = rivermouth(i)
            IF (nwrk_rs(iriv) > 1) THEN
               iworker = iwrk_rs(iriv)
               IF (nups_all(i) <= nave_rs(iriv)-nuc_wrk(iworker)) THEN

                  addr_ucat(i) = p_address_worker(iworker)

                  nuc_wrk(iworker) = nuc_wrk(iworker) + nups_all(i)
                  IF (nuc_wrk(iworker) == nave_rs(iriv)) THEN
                     iwrk_rs(iriv) = iwrk_rs(iriv) + 1
                  ENDIF

                  j = ucat_next(i)
                  IF (j > 0) THEN
                     DO WHILE (j > 0)
                        nups_all(j) = nups_all(j) - nups_all(i)
                        ithis = order_ucat(j)
                        j = ucat_next(j)
                     ENDDO
                  ELSE
                     ithis = ithis - 1
                  ENDIF
               ELSE
                  ithis = ithis - 1
               ENDIF
            ELSE
               iworker = minloc(wt_wrk(iwrkdsp+1:p_np_worker-1), dim=1) + iwrkdsp

               addr_ucat(i) = p_address_worker(iworker)

               wt_wrk(iworker) = wt_wrk(iworker) + wt_rs(iriv)
               ithis = ithis - 1
            ENDIF

         ENDDO

         deallocate (order_ucat)
         deallocate (nups_all  )
         deallocate (nuc_rs    )
         deallocate (iwrk_rs   )
         deallocate (nwrk_rs   )
         deallocate (nave_rs   )
         deallocate (wt_uc     )
         deallocate (wt_rs     )
         deallocate (wt_wrk    )
         deallocate (nuc_wrk   )

      ENDIF

      IF (p_is_master) THEN

         allocate(ucat_ucid (totalnumucat))
         ucat_ucid = (/(i, i = 1, totalnumucat)/)

         allocate (numucat_wrk       (0:p_np_worker-1))
         allocate (ucat_data_address (0:p_np_worker-1))

         DO iworker = 0, p_np_worker-1
            nucat = count(addr_ucat == p_address_worker(iworker))
            numucat_wrk(iworker) = nucat
            IF (nucat > 0) THEN
               allocate (ucat_data_address(iworker)%val (nucat))
               ucat_data_address(iworker)%val = &
                  pack(ucat_ucid, mask = (addr_ucat == p_address_worker(iworker)))
            ENDIF
         ENDDO

         deallocate (ucat_ucid)
      ENDIF

      CALL mpi_bcast (totalnumucat, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)

      ! send unit catchment index to workers
      IF (p_is_master) THEN

         DO iworker = 0, p_np_worker-1

            CALL mpi_send (numucat_wrk(iworker), 1, MPI_INTEGER, p_address_worker(iworker), &
               mpi_tag_mesg, p_comm_glb, p_err)

            nucat = numucat_wrk(iworker)
            IF (nucat > 0) THEN

               CALL mpi_send (ucat_data_address(iworker)%val, nucat, MPI_INTEGER, &
                  p_address_worker(iworker), mpi_tag_data, p_comm_glb, p_err)

               allocate (idata1d (nucat))

               idata1d = x_ucat (ucat_data_address(iworker)%val)
               CALL mpi_send (idata1d, nucat, MPI_INTEGER, &
                  p_address_worker(iworker), mpi_tag_data, p_comm_glb, p_err)

               idata1d = y_ucat (ucat_data_address(iworker)%val)
               CALL mpi_send (idata1d, nucat, MPI_INTEGER, &
                  p_address_worker(iworker), mpi_tag_data, p_comm_glb, p_err)

               deallocate (idata1d)
            ENDIF
         ENDDO

      ELSEIF (p_is_worker) THEN

         CALL mpi_recv (numucat, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

         IF (numucat > 0) THEN
            allocate (ucat_ucid (numucat))
            allocate (x_ucat    (numucat))
            allocate (y_ucat    (numucat))
            CALL mpi_recv (ucat_ucid, numucat, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)
            CALL mpi_recv (x_ucat, numucat, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)
            CALL mpi_recv (y_ucat, numucat, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF

      ENDIF

      IF (.not. p_is_worker) THEN
         numucat = 0
         IF (.not. allocated(ucat_ucid)) allocate (ucat_ucid (0))
         IF (.not. allocated(x_ucat   )) allocate (x_ucat    (0))
         IF (.not. allocated(y_ucat   )) allocate (y_ucat    (0))
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      numucat = totalnumucat

      allocate(ucat_ucid (totalnumucat))
      ucat_ucid = (/(i, i = 1, totalnumucat)/)

      allocate (numucat_wrk (0:0))
      numucat_wrk(0) = numucat

      allocate (ucat_data_address (0:0))
      allocate (ucat_data_address(0)%val (numucat))
      ucat_data_address(0)%val = ucat_ucid
#endif

      IF (allocated(addr_ucat)) deallocate(addr_ucat)

      ! ----- Part 1: between runoff input elements and unit catchments -----

#ifdef USEMPI
      CALL mpi_bcast (nlon_ucat, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (nlat_ucat, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif

      CALL griducat%define_by_ndims (nlon_ucat, nlat_ucat)

      CALL build_worker_remapdata (landpatch, griducat, remap_patch2inpm)

      IF (DEF_UnitCatchment_regional) &
         CALL check_regional_coverage (nlon_ucat, nlat_ucat, idmap_x, idmap_y, remap_patch2inpm)

      IF (p_is_worker) THEN
         numinpm = remap_patch2inpm%num_grid
         IF (numinpm > 0) THEN
            allocate (inpm_gdid (numinpm))
            inpm_gdid = remap_patch2inpm%ids_me
         ENDIF
      ENDIF


      IF (p_is_master) THEN

         inpn = size(idmap_x,1)

         allocate(idmap_gd2uc (inpn,totalnumucat))

         idmap_gd2uc = (idmap_y-1)*nlon_ucat + idmap_x

         WHERE ((area_gd2uc <= 0) .or. (idmap_gd2uc <= 0))
            idmap_gd2uc = 0
            area_gd2uc  = 0.
         END WHERE

         allocate (nucat_g2d (nlon_ucat,nlat_ucat))
         nucat_g2d(:,:) = 0

         DO i = 1, totalnumucat
            DO j = 1, inpn
               IF (idmap_gd2uc(j,i) > 0) THEN
                  nucat_g2d(idmap_x(j,i),idmap_y(j,i)) = nucat_g2d(idmap_x(j,i),idmap_y(j,i)) + 1
               ENDIF
            ENDDO
         ENDDO

         nucpart = maxval(nucat_g2d)
         ngrdall = count(nucat_g2d > 0)

         allocate (allgrd_in_inp (ngrdall))

         igrd = 0
         DO i = 1, nlat_ucat
            DO j = 1, nlon_ucat
               IF (nucat_g2d(j,i) > 0) THEN
                  igrd = igrd + 1
                  allgrd_in_inp(igrd) = (i-1)*nlon_ucat + j
               ENDIF
            ENDDO
         ENDDO

         allocate (idmap_uc2gd_all (nucpart, ngrdall));  idmap_uc2gd_all(:,:) = 0
         allocate (area_uc2gd_all  (nucpart, ngrdall));  area_uc2gd_all (:,:) = 0.

         allocate (iucat_g (ngrdall)); iucat_g(:) = 0

         DO i = 1, totalnumucat
            DO j = 1, inpn
               IF (idmap_gd2uc(j,i) > 0) THEN
                  iloc = find_in_sorted_list1 (idmap_gd2uc(j,i), ngrdall, allgrd_in_inp(1:ngrdall))
                  iucat_g(iloc) = iucat_g(iloc) + 1
                  idmap_uc2gd_all(iucat_g(iloc),iloc) = i
                  area_uc2gd_all (iucat_g(iloc),iloc) = area_gd2uc(j,i)
               ENDIF
            ENDDO
         ENDDO

      ENDIF

#ifdef USEMPI
      CALL mpi_bcast (inpn, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)

      IF (p_is_master) THEN

         DO iworker = 0, p_np_worker-1

            nucat = numucat_wrk(iworker)

            IF (nucat > 0) THEN
               allocate (idata2d (inpn, nucat))
               DO i = 1, nucat
                  idata2d(:,i) = idmap_gd2uc(:,ucat_data_address(iworker)%val(i))
               ENDDO

               CALL mpi_send (idata2d, inpn*nucat, MPI_INTEGER, p_address_worker(iworker), &
                  mpi_tag_data, p_comm_glb, p_err)

               allocate (rdata2d (inpn, nucat))
               DO i = 1, nucat
                  rdata2d(:,i) = area_gd2uc(:,ucat_data_address(iworker)%val(i))
               ENDDO

               CALL mpi_send (rdata2d, inpn*nucat, MPI_REAL8, p_address_worker(iworker), &
                  mpi_tag_data, p_comm_glb, p_err)

               deallocate (idata2d)
               deallocate (rdata2d)
            ENDIF
         ENDDO

         deallocate (idmap_gd2uc)
         deallocate (area_gd2uc )

      ELSEIF (p_is_worker) THEN

         IF (numucat > 0) THEN

            allocate (idmap_gd2uc (inpn, numucat))
            CALL mpi_recv (idmap_gd2uc, inpn*numucat, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (area_gd2uc (inpn, numucat))
            CALL mpi_recv (area_gd2uc, inpn*numucat, MPI_REAL8, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)

         ENDIF

      ENDIF

      CALL mpi_bcast (nucpart, 1, mpi_integer, p_address_master, p_comm_glb, p_err)

      IF (p_is_master) THEN

         DO iworker = 0, p_np_worker-1

            CALL mpi_recv (ngrd, 1, MPI_INTEGER, &
               p_address_worker(iworker), mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            IF (ngrd > 0) THEN

               allocate (grdindex (ngrd))
               allocate (idata2d  (nucpart, ngrd));  idata2d(:,:) = 0
               allocate (rdata2d  (nucpart, ngrd));  rdata2d(:,:) = 0.

               CALL mpi_recv (grdindex, ngrd, MPI_INTEGER, &
                  p_address_worker(iworker), mpi_tag_data, p_comm_glb, p_stat, p_err)

               DO i = 1, ngrd
                  iloc = find_in_sorted_list1 (grdindex(i), ngrdall, allgrd_in_inp(1:ngrdall))
                  IF (iloc > 0) THEN
                     idata2d(:,i) = idmap_uc2gd_all(:,iloc)
                     rdata2d(:,i) = area_uc2gd_all (:,iloc)
                  ENDIF
               ENDDO

               CALL mpi_send (idata2d, nucpart*ngrd, MPI_INTEGER, p_address_worker(iworker), &
                  mpi_tag_data, p_comm_glb, p_err)
               CALL mpi_send (rdata2d, nucpart*ngrd, MPI_REAL8,   p_address_worker(iworker), &
                  mpi_tag_data, p_comm_glb, p_err)

               deallocate (grdindex)
               deallocate (idata2d )
               deallocate (rdata2d )
            ENDIF
         ENDDO

      ELSEIF (p_is_worker) THEN

         CALL mpi_send (numinpm, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)

         IF (numinpm > 0) THEN

            CALL mpi_send (inpm_gdid, numinpm, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_err)

            allocate (idmap_uc2gd (nucpart,numinpm))
            CALL mpi_recv (idmap_uc2gd, nucpart*numinpm, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (area_uc2gd (nucpart,numinpm))
            CALL mpi_recv (area_uc2gd, nucpart*numinpm, MPI_REAL8, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)

         ENDIF

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      allocate (idmap_uc2gd (nucpart,numinpm))
      allocate (area_uc2gd  (nucpart,numinpm))
      idmap_uc2gd = 0
      area_uc2gd  = 0.

      DO i = 1, numinpm
         iloc = find_in_sorted_list1 (inpm_gdid(i), ngrdall, allgrd_in_inp(1:ngrdall))
         IF (iloc > 0) THEN
            idmap_uc2gd(:,i) = idmap_uc2gd_all(:,iloc)
            area_uc2gd (:,i) = area_uc2gd_all (:,iloc)
         ENDIF
      ENDDO
#endif

      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            allocate (ucat_gdid (numucat))
            ucat_gdid = (y_ucat-1)*nlon_ucat + x_ucat
         ENDIF
      ENDIF

      CALL build_worker_pushdata (numinpm, inpm_gdid, numucat, idmap_gd2uc, area_gd2uc, push_inpm2ucat)
      CALL build_worker_pushdata (numucat, ucat_ucid, numinpm, idmap_uc2gd, area_uc2gd, push_ucat2inpm)
      CALL build_worker_pushdata (numucat, ucat_gdid, numinpm, inpm_gdid, push_ucat2grid)
      CALL build_worker_pushdata (numinpm, inpm_gdid, numinpm, inpm_gdid, allreduce_inpm)

      IF (p_is_master) THEN
         deallocate (idmap_x        )
         deallocate (idmap_y        )
         deallocate (allgrd_in_inp  )
         deallocate (nucat_g2d      )
         deallocate (iucat_g        )
         deallocate (idmap_uc2gd_all)
         deallocate (area_uc2gd_all )
      ENDIF

      ! ----- Part 2: between upstream and downstream unit catchments -----

      IF (p_is_master) THEN

         upnmax = maxval(nups_nst)
         allocate (ucat_ups (upnmax,totalnumucat))
         ucat_ups(:,:) = 0

         iups_nst(:) = 0
         DO i = 1, totalnumucat
            j = ucat_next(i)
            IF (j > 0) THEN
               iups_nst(j) = iups_nst(j) + 1
               ucat_ups(iups_nst(j),j) = i
            ENDIF
         ENDDO

      ENDIF


#ifdef USEMPI
      CALL mpi_bcast (upnmax, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)

      IF (p_is_master) THEN

         DO iworker = 0, p_np_worker-1

            nucat = numucat_wrk(iworker)

            IF (nucat > 0) THEN
               allocate (idata1d (nucat))
               idata1d = ucat_next(ucat_data_address(iworker)%val)
               CALL mpi_send (idata1d, nucat, MPI_INTEGER, p_address_worker(iworker), &
                  mpi_tag_data, p_comm_glb, p_err)

               allocate (idata2d (upnmax,nucat))
               DO i = 1, nucat
                  idata2d(:,i) = ucat_ups(:,ucat_data_address(iworker)%val(i))
               ENDDO
               CALL mpi_send (idata2d, upnmax*nucat, MPI_INTEGER, p_address_worker(iworker), &
                  mpi_tag_data, p_comm_glb, p_err)

               deallocate (idata1d)
               deallocate (idata2d)
            ENDIF
         ENDDO

         deallocate (ucat_ups )

      ELSEIF (p_is_worker) THEN

         IF (numucat > 0) THEN

            allocate (ucat_next (numucat))
            CALL mpi_recv (ucat_next, numucat, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (ucat_ups (upnmax, numucat))
            CALL mpi_recv (ucat_ups, upnmax*numucat, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)

         ENDIF

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)

      IF (.not. p_is_worker) THEN
         IF (.not. allocated(ucat_next)) allocate (ucat_next (0))
         IF (.not. allocated(ucat_ups )) allocate (ucat_ups  (upnmax, 0))
         IF (.not. allocated(wts_ups  )) allocate (wts_ups   (upnmax, 0))
      ENDIF
#endif

      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            allocate (wts_ups (upnmax,numucat))
            wts_ups(:,:) = 1.
         ENDIF
      ENDIF

      CALL build_worker_pushdata (numucat, ucat_ucid, numucat, ucat_next, push_next2ucat)
      CALL build_worker_pushdata (numucat, ucat_ucid, numucat, ucat_ups,  wts_ups, push_ups2ucat )
      CALL riverlake_pushcache_build (push_ups2ucat, cache_ups2ucat)

#ifdef CoLMDEBUG
      ! IF (p_is_worker) THEN
      !    write(*,'(A,I0,A,I0,A,I0,A)') 'worker ', p_iam_worker, ' has ', numucat, &
      !       ' unit catchment with ', sum(push_next2ucat%n_from_other), ' downstream to other workers'
      ! ENDIF
#endif

      ! ----- Part 3: river systems -----

#ifdef USEMPI
      IF (p_is_master) THEN
         DO iworker = 0, p_np_worker-1
            nucat = numucat_wrk(iworker)
            IF (nucat > 0) THEN
               allocate (idata1d (nucat))
               idata1d = rivermouth(ucat_data_address(iworker)%val)
               CALL mpi_send (idata1d, nucat, MPI_INTEGER, p_address_worker(iworker), &
                  mpi_tag_data, p_comm_glb, p_err)
               deallocate (idata1d)
            ENDIF
         ENDDO
      ELSEIF (p_is_worker) THEN
         IF (numucat > 0) THEN
            allocate (rivermouth (numucat))
            CALL mpi_recv (rivermouth, numucat, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF
      ENDIF

      rivsys_by_multiple_procs = .false.
      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            color = maxval(rivermouth)
            CALL mpi_comm_split (p_comm_worker, color, p_iam_worker, p_comm_rivsys, p_err)
         ELSE
            CALL mpi_comm_split (p_comm_worker, MPI_UNDEFINED, p_iam_worker, p_comm_rivsys, p_err)
         ENDIF

         IF (p_comm_rivsys /= MPI_COMM_NULL) THEN
            CALL mpi_comm_size (p_comm_rivsys, p_np_rivsys, p_err)
            IF (p_np_rivsys > 1) THEN
               rivsys_by_multiple_procs = .true.
            ENDIF
         ENDIF
      ENDIF

      ! With bifurcation the time step groups are the clusters of river systems
      ! connected by bifurcation pathways (see irivsys_shared above).  The
      ! river-mouth labels sent above have already set up p_comm_rivsys, which is
      ! not used in bifurcation mode; here the master relabels every river system
      ! with its cluster, and sends the labels and the index of clusters that are
      ! shared between workers.  The grouping below then works per cluster.
      IF (DEF_USE_BIFURCATION) THEN
         IF (p_is_master) THEN
            CALL bifurcation_dt_groups (parafile, numrivmth, rivermouth, shared_ucat)
            DO iworker = 0, p_np_worker-1
               nucat = numucat_wrk(iworker)
               IF (nucat > 0) THEN
                  allocate (idata1d (nucat))
                  idata1d = rivermouth(ucat_data_address(iworker)%val)
                  CALL mpi_send (idata1d, nucat, MPI_INTEGER, p_address_worker(iworker), &
                     mpi_tag_data, p_comm_glb, p_err)
                  idata1d = shared_ucat(ucat_data_address(iworker)%val)
                  CALL mpi_send (idata1d, nucat, MPI_INTEGER, p_address_worker(iworker), &
                     mpi_tag_data, p_comm_glb, p_err)
                  deallocate (idata1d)
               ENDIF
            ENDDO
         ELSEIF (p_is_worker) THEN
            IF (numucat > 0) THEN
               allocate (shared_ucat (numucat))
               CALL mpi_recv (rivermouth, numucat, MPI_INTEGER, p_address_master, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL mpi_recv (shared_ucat, numucat, MPI_INTEGER, p_address_master, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)
            ENDIF
         ENDIF
         CALL mpi_bcast (numrivsys_shared, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      ENDIF

#else
      rivsys_by_multiple_procs = .false.
#endif

      IF (p_is_worker) THEN

         numrivsys = 0
         allocate (irivsys (numucat))
         allocate (irivsys_shared (numucat))
         irivsys_shared(:) = 0

         ! In bifurcation mode the groups are always formed per cluster, also on
         ! workers that hold a part of a large river system; groups shared with
         ! other workers are synchronised through irivsys_shared instead.
         IF ((.not. rivsys_by_multiple_procs) .or. DEF_USE_BIFURCATION) THEN
            IF (numucat > 0) THEN

               allocate (order_ucat (numucat))
               order_ucat = (/(i, i = 1, numucat)/)

               CALL quicksort (numucat, rivermouth, order_ucat)

               numrivsys = 1
               irivsys(order_ucat(1)) = numrivsys
               DO i = 2, numucat
                  IF (rivermouth(i) /= rivermouth(i-1)) THEN
                     numrivsys = numrivsys + 1
                  ENDIF
                  irivsys(order_ucat(i)) = numrivsys
               ENDDO

               IF (allocated(shared_ucat)) THEN
                  DO i = 1, numucat
                     irivsys_shared(irivsys(i)) = shared_ucat(i)
                  ENDDO
               ENDIF

            ENDIF
         ELSE
            numrivsys  = 1
            irivsys(:) = 1
         ENDIF

      ENDIF

      IF (allocated(rivermouth )) deallocate(rivermouth )
      IF (allocated(shared_ucat)) deallocate(shared_ucat)
      IF (allocated(order_ucat )) deallocate(order_ucat )

      ! ----- Parameters for River and Lake -----

      CALL readin_riverlake_parameter (parafile, 'topo_rivelv',    rdata1d = topo_rivelv   )
      CALL readin_riverlake_parameter (parafile, 'topo_rivhgt',    rdata1d = topo_rivhgt   )
      CALL readin_riverlake_parameter (parafile, 'topo_rivlen',    rdata1d = topo_rivlen   )
      CALL readin_riverlake_parameter (parafile, 'topo_rivman',    rdata1d = topo_rivman   )
      CALL readin_riverlake_parameter (parafile, 'topo_rivwth',    rdata1d = topo_rivwth   )
      CALL readin_riverlake_parameter (parafile, 'topo_rivstomax', rdata1d = topo_rivstomax)
      CALL readin_riverlake_parameter (parafile, 'topo_area',      rdata1d = topo_area     )
      CALL readin_riverlake_parameter (parafile, 'topo_fldhgt',    rdata2d = topo_fldhgt   )

      IF (DEF_USE_LEVEE) THEN
         CALL readin_riverlake_parameter (parafile, 'levee_frc', rdata1d = levee_frc_data)
         CALL readin_riverlake_parameter (parafile, 'levee_hgt', rdata1d = levee_hgt_data)
      ENDIF

      IF (DEF_USE_BIFURCATION) THEN
         CALL read_and_distribute_bifurcation (parafile)
      ENDIF

      IF (p_is_worker) THEN
         IF (numucat > 0) THEN

            allocate (lake_type (numucat))
            lake_type(:) = 0

            allocate (topo_rivare (numucat))
            topo_rivare = topo_rivstomax / topo_rivhgt

            allocate (floodplain_curve (numucat))

            DO i = 1, numucat
               floodplain_curve(i)%nlfp      = size(topo_fldhgt,1)
               floodplain_curve(i)%rivhgt    = topo_rivhgt(i)
               floodplain_curve(i)%rivstomax = topo_rivstomax(i)
               floodplain_curve(i)%rivare    = topo_rivare(i)

               allocate (floodplain_curve(i)%flphgt    (0:floodplain_curve(i)%nlfp))
               allocate (floodplain_curve(i)%flparea   (0:floodplain_curve(i)%nlfp))
               allocate (floodplain_curve(i)%flpaccare (0:floodplain_curve(i)%nlfp))
               allocate (floodplain_curve(i)%flpstomax (0:floodplain_curve(i)%nlfp))

               floodplain_curve(i)%flphgt(0)  = 0.
               floodplain_curve(i)%flphgt(1:) = topo_fldhgt(:,i)

               floodplain_curve(i)%flparea(0)  = 0.
               floodplain_curve(i)%flparea(1:) = topo_area(i) / floodplain_curve(i)%nlfp

               floodplain_curve(i)%flpaccare(0) = 0.
               DO j = 1, floodplain_curve(i)%nlfp
                  floodplain_curve(i)%flpaccare(j) = &
                     floodplain_curve(i)%flpaccare(j-1) + floodplain_curve(i)%flparea(j)
               ENDDO

               ! [GridRiverLakeFlow] Bug fix: upstream integrates the storage of layer j
               ! with the per-layer area increments, 0.5*(flparea(j)+flparea(j-1)).  The
               ! flooded area at height h is the cumulative area flpaccare, so the storage
               ! between two heights is the trapezoid of the cumulative areas, as used by
               ! volume(), depth() and floodarea() of this type.  With the increments the
               ! layer thresholds disagree with those functions: the volume jumps down at
               ! each layer top (non-monotonic) and depth(volume(d)) /= d.  For example,
               ! with three 1 m layers over 3000 m2, volume(4 m) = 2200 m3 but
               ! volume(4.000001 m) = 1700 m3, and depth(2200 m3) = 4.236 m.
               ! This changes results of all runs whose water rises onto the floodplain
               ! (above the first layer).  The flood feedback always used this formula.
               floodplain_curve(i)%flpstomax(0) = 0.
               DO j = 1, floodplain_curve(i)%nlfp
                  floodplain_curve(i)%flpstomax(j) = floodplain_curve(i)%flpstomax(j-1)        &
                     + 0.5 * (floodplain_curve(i)%flpaccare(j) + floodplain_curve(i)%flpaccare(j-1)) &
                           * (floodplain_curve(i)%flphgt(j)  - floodplain_curve(i)%flphgt(j-1))
               ENDDO
            ENDDO

         ENDIF
      ENDIF

      ! Inactive ranks also pass these arrays to assumed-shape arguments.
      allocate (bedelv_next (numucat))
      allocate (outletwth   (numucat))

      CALL worker_push_data (push_next2ucat, topo_rivelv, bedelv_next, fillvalue = spval)
      CALL worker_push_data (push_next2ucat, topo_rivwth, outletwth  , fillvalue = spval)

      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            ! Outlet width is the width of this cell's own channel (it was the
            ! mean of this and the downstream channel), as in CaMa-Flood
            ! (cmf_calc_outflw_mod: flow area and outflow use D2RIVWTH of the
            ! upstream cell).  The momentum of a cell is divided by its own
            ! channel area topo_rivare = rivwth*rivlen, while the CFL limit
            ! uses rivlen.  Where a narrow tributary (~5 m) enters a wide river
            ! or lake (~2000 m) the mean width was hundreds of times the cell's
            ! own width, the effective length topo_rivare/outletwth only tens
            ! of metres, and the outflow flipped sign every other substep: a
            ! 1e-13 velocity perturbation grew to O(1) within a day (0.25 deg
            ! global network, 3 days: 73-89 cells with discharge >= 10 m3/s off
            ! by more than 1%, none now).  With the own width the effective
            ! length is rivlen and the CFL limit holds; the adaptive dt now
            ! stays at its 60 s cap (27% fewer substeps over 30 days).  The
            ! routing results change wherever the two widths differ (69% of
            ! the faces at 0.25 deg, 50% at 0.1 deg; most rivers widen
            ! downstream, so the outlet becomes narrower): on a cold start the
            ! discharge of the 10 largest river mouths averaged over 30 days is
            ! 0-13% lower, the river storage up to 1.6% higher.
            ! (min(w, w_next) was also tested: same stability and cost, about
            ! twice the change of the results.)
            WHERE (ucat_next > 0)
               outletwth = topo_rivwth
            ELSEWHERE
               outletwth = topo_rivwth
            END WHERE
         ENDIF
      ENDIF

      ! ----- Mask of Grids with all upstream area in the simulation region -----

      IF (p_is_master) allocate (ucat_area_all (totalnumucat))

#ifdef USEMPI
      IF (p_is_worker) THEN

         IF (numucat > 0) THEN
            CALL mpi_send (push_inpm2ucat%sum_area, numucat, MPI_REAL8, p_address_master, &
               mpi_tag_data, p_comm_glb, p_err)
         ENDIF

      ELSEIF (p_is_master) THEN

         DO iworker = 0, p_np_worker-1
            IF (numucat_wrk(iworker) > 0) THEN

               allocate (rdata1d (numucat_wrk(iworker)))
               CALL mpi_recv (rdata1d, numucat_wrk(iworker), MPI_REAL8, p_address_worker(iworker), &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)

               ucat_area_all(ucat_data_address(iworker)%val) = rdata1d

               deallocate (rdata1d)
            ENDIF

         ENDDO
      ENDIF
#else
      ucat_area_all = push_inpm2ucat%sum_area
#endif

      IF (p_is_master) THEN

         allocate (allups_mask_ucat (totalnumucat))
         allups_mask_ucat (:) = 0

         iups_nst(:) = 0
         DO i = 1, totalnumucat
            j = uc_up2down(i)
            IF (ucat_area_all(j) > 0.) THEN
               IF (iups_nst(j) == nups_nst(j)) THEN

                  allups_mask_ucat(j) = 1

                  IF (ucat_next(j) > 0) THEN
                     iups_nst(ucat_next(j)) = iups_nst(ucat_next(j)) + 1
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

      ENDIF

#ifdef USEMPI
      IF (p_is_master) THEN
         DO iworker = 0, p_np_worker-1
            IF (numucat_wrk(iworker) > 0) THEN
               allocate (rdata1d (numucat_wrk(iworker)))
               rdata1d = allups_mask_ucat(ucat_data_address(iworker)%val)

               CALL mpi_send (rdata1d, numucat_wrk(iworker), MPI_REAL8, p_address_worker(iworker), &
                  mpi_tag_data, p_comm_glb, p_err)

               deallocate (rdata1d)
            ENDIF
         ENDDO
      ELSEIF (p_is_worker) THEN
         IF (numucat > 0) THEN
            allocate (allups_mask_ucat (numucat))
            CALL mpi_recv (allups_mask_ucat, numucat, MPI_REAL8, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF
      ENDIF
#endif


      IF (.not. allocated(allups_mask_ucat)) allocate (allups_mask_ucat (0))

      IF (allocated (uc_up2down   )) deallocate (uc_up2down   )
      IF (allocated (nups_nst     )) deallocate (nups_nst     )
      IF (allocated (iups_nst     )) deallocate (iups_nst     )
      IF (allocated (ucat_area_all)) deallocate (ucat_area_all)

   END SUBROUTINE build_riverlake_network

   SUBROUTINE check_regional_coverage (nlon, nlat, idmap_x, idmap_y, remap)

   ! The regional network must take the runoff of every runoff-input cell that the
   ! land patches overlap and that the full network (DEF_UnitCatchment_file) routes;
   ! otherwise that runoff would be dropped without notice.  mkinidata cuts the
   ! network so, from the land patches of every land-cover year of the run, but the
   ! run may reach a year that mkinidata had no land data for, or use a file cut for
   ! another period or domain.  Called at start-up and after each LULCC.
   ! Collective over all processes.  idmap_x/idmap_y (the regional input matrix) are
   ! used on the master only.

   USE MOD_SPMD_Task
   USE MOD_Namelist,      only: DEF_UnitCatchment_file
   USE MOD_NetCDFSerial
   IMPLICIT NONE

   integer, intent(in) :: nlon, nlat
   integer, allocatable, intent(in) :: idmap_x(:,:), idmap_y(:,:)
   type(worker_remapdata_type), intent(in) :: remap

   logical, allocatable :: touched(:,:), routed_full(:,:), routed_regional(:,:)
   integer, allocatable :: src_x(:,:), src_y(:,:), ids(:)
   integer :: ng, ix, iy, nmiss
#ifdef USEMPI
   integer :: iworker
#endif

      IF (p_is_master) THEN
         allocate (touched (nlon,nlat));  touched = .false.
      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN
         ng = remap%num_grid
         CALL mpi_send (ng, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         IF (ng > 0) CALL mpi_send (remap%ids_me, ng, MPI_INTEGER, p_address_master, mpi_tag_data, &
            p_comm_glb, p_err)
      ELSEIF (p_is_master) THEN
         DO iworker = 0, p_np_worker-1
            CALL mpi_recv (ng, 1, MPI_INTEGER, p_address_worker(iworker), mpi_tag_mesg, p_comm_glb, &
               p_stat, p_err)
            IF (ng > 0) THEN
               allocate (ids (ng))
               CALL mpi_recv (ids, ng, MPI_INTEGER, p_address_worker(iworker), mpi_tag_data, &
                  p_comm_glb, p_stat, p_err)
               CALL mark (ids)
               deallocate (ids)
            ENDIF
         ENDDO
      ENDIF
#else
      IF (remap%num_grid > 0) CALL mark (remap%ids_me)
#endif

      nmiss = 0
      IF (p_is_master) THEN
         CALL ncio_read_serial (DEF_UnitCatchment_file, 'inpmat_x', src_x)
         CALL ncio_read_serial (DEF_UnitCatchment_file, 'inpmat_y', src_y)
         allocate (routed_full     (nlon,nlat));  CALL routed (src_x, src_y, routed_full)
         allocate (routed_regional (nlon,nlat));  CALL routed (idmap_x, idmap_y, routed_regional)

         nmiss = count(touched .and. routed_full .and. .not. routed_regional)
         IF (nmiss > 0) THEN
            write(*,'(A,I0,A)') 'ERROR: ', nmiss, ' runoff-input cells under the land patches are routed by'
            write(*,'(A)') '   DEF_UnitCatchment_file but not by the regional unit-catchment network.'
            write(*,'(A)') '   Run mkinidata again for this domain and simulation period'
            write(*,'(A)') '   (with LULCC, after mksrfdata has made the land data of every year).'
         ENDIF
         deallocate (touched, routed_full, routed_regional, src_x, src_y)
      ENDIF

#ifdef USEMPI
      CALL mpi_bcast (nmiss, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif
      IF (nmiss > 0) CALL CoLM_stop ('regional unit-catchment network does not cover the land patches')

   CONTAINS

      SUBROUTINE mark (cell_ids)
      integer, intent(in) :: cell_ids(:)
      integer :: k
         DO k = 1, size(cell_ids)
            ix = mod(cell_ids(k)-1, nlon) + 1
            iy = (cell_ids(k)-1) / nlon + 1
            touched(ix,iy) = .true.
         ENDDO
      END SUBROUTINE mark

      SUBROUTINE routed (mx, my, cell)
      integer, intent(in)  :: mx(:,:), my(:,:)
      logical, intent(out) :: cell(:,:)
      integer :: j, k
         cell = .false.
         DO k = 1, size(mx,2)
            DO j = 1, size(mx,1)
               IF (mx(j,k) >= 1 .and. mx(j,k) <= nlon .and. my(j,k) >= 1 .and. my(j,k) <= nlat) &
                  cell(mx(j,k), my(j,k)) = .true.
            ENDDO
         ENDDO
      END SUBROUTINE routed

   END SUBROUTINE check_regional_coverage

   SUBROUTINE verify_regional_network (file_regional, x_regional, y_regional)


   USE MOD_Namelist,      only: DEF_UnitCatchment_file
   USE MOD_NetCDFSerial
   IMPLICIT NONE

   character(len=*), intent(in) :: file_regional
   integer,          intent(in) :: x_regional(:), y_regional(:)

   integer, allocatable :: src_index(:), x_source(:), y_source(:)
   logical :: consistent

      CALL ncio_read_serial (file_regional, 'seq_src_index', src_index)
      CALL ncio_read_serial (DEF_UnitCatchment_file, 'seq_x', x_source)
      CALL ncio_read_serial (DEF_UnitCatchment_file, 'seq_y', y_source)

      consistent = (size(src_index) == size(x_regional))
      IF (consistent) consistent = all(src_index >= 1) .and. all(src_index <= size(x_source))
      IF (consistent) consistent = all(x_source(src_index) == x_regional) .and. &
                                   all(y_source(src_index) == y_regional)

      IF (.not. consistent) THEN
         write(*,'(A)') 'ERROR: the regional unit-catchment network does not belong to DEF_UnitCatchment_file.'
         write(*,'(2A)') '   regional file: ', trim(file_regional)
         write(*,'(A)') '   Run mkinidata again with the current DEF_UnitCatchment_file.'
         CALL CoLM_stop ()
      ENDIF

      deallocate (src_index, x_source, y_source)

   END SUBROUTINE verify_regional_network

   ! ---------
   SUBROUTINE readin_riverlake_parameter (parafile, varname, rdata1d, rdata2d, idata1d)

   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   IMPLICIT NONE

   character(len=*), intent(in) :: parafile
   character(len=*), intent(in) :: varname

   real(r8), allocatable, intent(inout), optional :: rdata1d (:)
   real(r8), allocatable, intent(inout), optional :: rdata2d (:,:)
   integer,  allocatable, intent(inout), optional :: idata1d (:)

   ! Local Variables
   integer :: iworker, nucat, ndim1, i
   real(r8), allocatable :: rsend1d (:)
   real(r8), allocatable :: rsend2d (:,:)
   integer,  allocatable :: isend1d (:)

      IF (p_is_master) THEN
         IF (present(rdata1d))  CALL ncio_read_serial (parafile, varname, rdata1d)
         IF (present(rdata2d))  CALL ncio_read_serial (parafile, varname, rdata2d)
         IF (present(idata1d))  CALL ncio_read_serial (parafile, varname, idata1d)
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (present(rdata2d)) THEN
         IF (p_is_master) ndim1 = size(rdata2d,1)
         CALL mpi_bcast (ndim1, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      ENDIF

      ! send unit catchment index to workers
      IF (p_is_master) THEN

         DO iworker = 0, p_np_worker-1

            nucat = numucat_wrk(iworker)

            IF (nucat > 0) THEN
               IF (present(rdata1d)) THEN
                  allocate (rsend1d (nucat))

                  rsend1d = rdata1d(ucat_data_address(iworker)%val)
                  CALL mpi_send (rsend1d, nucat, MPI_REAL8, p_address_worker(iworker), &
                     mpi_tag_data, p_comm_glb, p_err)

                  deallocate (rsend1d)
               ENDIF

               IF (present(rdata2d)) THEN
                  allocate (rsend2d (ndim1,nucat))

                  DO i = 1, nucat
                     rsend2d(:,i) = rdata2d(:,ucat_data_address(iworker)%val(i))
                  ENDDO
                  CALL mpi_send (rsend2d, ndim1*nucat, MPI_REAL8, p_address_worker(iworker), &
                     mpi_tag_data, p_comm_glb, p_err)

                  deallocate (rsend2d)
               ENDIF

               IF (present(idata1d)) THEN
                  allocate (isend1d (nucat))

                  isend1d = idata1d(ucat_data_address(iworker)%val)
                  CALL mpi_send (isend1d, nucat, MPI_INTEGER, p_address_worker(iworker), &
                     mpi_tag_data, p_comm_glb, p_err)

                  deallocate (isend1d)
               ENDIF
            ENDIF
         ENDDO

         IF (present(rdata1d))  deallocate (rdata1d)
         IF (present(rdata2d))  deallocate (rdata2d)
         IF (present(idata1d))  deallocate (idata1d)

      ELSEIF (p_is_worker) THEN

         IF (numucat > 0) THEN
            IF (present(rdata1d)) THEN
               allocate (rdata1d (numucat))
               CALL mpi_recv (rdata1d, numucat, MPI_REAL8, p_address_master, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)
            ENDIF

            IF (present(rdata2d)) THEN
               allocate (rdata2d (ndim1,numucat))
               CALL mpi_recv (rdata2d, ndim1*numucat, MPI_REAL8, p_address_master, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)
            ENDIF

            IF (present(idata1d)) THEN
               allocate (idata1d (numucat))
               CALL mpi_recv (idata1d, numucat, MPI_INTEGER, p_address_master, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)
            ENDIF
         ELSE
            IF (present(rdata1d)) allocate (rdata1d (0))
            IF (present(rdata2d)) allocate (rdata2d (ndim1,0))
            IF (present(idata1d)) allocate (idata1d (0))
         ENDIF

      ENDIF

      ! Master releases its global input above; IO ranks never read it.
      ! Keep legal actual arguments for the subsequent all-rank push calls.
      IF (.not. p_is_worker) THEN
         IF (present(rdata1d)) allocate (rdata1d (0))
         IF (present(rdata2d)) allocate (rdata2d (ndim1,0))
         IF (present(idata1d)) allocate (idata1d (0))
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE readin_riverlake_parameter

   !
   FUNCTION retrieve_depth_from_volume (this, volume) result(depth)

   IMPLICIT NONE

   class(vol_dep_curve_type) :: this
   real(r8), intent(in)      :: volume
   real(r8) :: depth

   ! Local Variables
   real(r8) :: v0, g
   integer  :: i

      v0 = volume - this%rivstomax
      IF (v0 <= 0) THEN
         depth = volume / this%rivare
      ELSE
         i = 1
         DO WHILE (i <= this%nlfp)
            IF (v0 > this%flpstomax(i)) THEN
               i = i + 1
            ELSE
               EXIT
            ENDIF
         ENDDO

         IF (i == this%nlfp+1) THEN
            depth = this%rivhgt + this%flphgt(this%nlfp) &
               + (v0-this%flpstomax(this%nlfp)) / this%flpaccare(this%nlfp)
         ELSE
            g = (this%flphgt(i)-this%flphgt(i-1))/this%flparea(i)
            depth = this%rivhgt + this%flphgt(i-1) &
               + g * (-this%flpaccare(i-1)+sqrt((this%flpaccare(i-1))**2+2*(v0-this%flpstomax(i-1))/g))
         ENDIF
      ENDIF

   END FUNCTION retrieve_depth_from_volume

   !
   FUNCTION retrieve_volume_from_depth (this, depth) result(volume)

   IMPLICIT NONE

   class(vol_dep_curve_type) :: this
   real(r8), intent(in)  :: depth
   real(r8) :: volume

   ! Local Variables
   real(r8) :: h, d
   integer  :: i

      IF (depth <= this%rivhgt) THEN
         volume = this%rivare * depth
      ELSE
         i = 1
         DO WHILE (i <= this%nlfp)
            IF (depth > this%rivhgt+this%flphgt(i)) THEN
               i = i + 1
            ELSE
               EXIT
            ENDIF
         ENDDO

         d = depth - this%rivhgt - this%flphgt(i-1)
         IF (i == this%nlfp+1) THEN
            volume = this%rivstomax + this%flpstomax(this%nlfp) + d * this%flpaccare(this%nlfp)
         ELSE
            h = this%flphgt(i)-this%flphgt(i-1)
            volume = this%rivstomax + this%flpstomax(i-1) &
               + (d/h*this%flparea(i)+2*this%flpaccare(i-1))*d*0.5
         ENDIF
      ENDIF

   END FUNCTION retrieve_volume_from_depth

   !
   FUNCTION retrieve_area_from_depth (this, depth) result(area)

   IMPLICIT NONE

   class(vol_dep_curve_type) :: this
   real(r8), intent(in)  :: depth
   real(r8) :: area

   ! Local Variables
   real(r8) :: h, d
   integer  :: i

      IF (depth <= this%rivhgt) THEN
         area = 0.
      ELSE
         i = 1
         DO WHILE (i <= this%nlfp)
            IF (depth > this%rivhgt+this%flphgt(i)) THEN
               i = i + 1
            ELSE
               EXIT
            ENDIF
         ENDDO

         IF (i == this%nlfp+1) THEN
            area = this%flpaccare(this%nlfp)
         ELSE
            h = this%flphgt(i)-this%flphgt(i-1)
            d = depth - this%rivhgt - this%flphgt(i-1)
            area = this%flpaccare(i-1) + d/h * this%flparea(i)
         ENDIF
      ENDIF

   END FUNCTION retrieve_area_from_depth

   ! ---
   SUBROUTINE vol_depth_curve_free_mem (this)

   IMPLICIT NONE
   type(vol_dep_curve_type) :: this

      IF (allocated(this%flphgt   )) deallocate (this%flphgt   )
      IF (allocated(this%flparea  )) deallocate (this%flparea  )
      IF (allocated(this%flpaccare)) deallocate (this%flpaccare)
      IF (allocated(this%flpstomax)) deallocate (this%flpstomax)

   END SUBROUTINE vol_depth_curve_free_mem

   ! ---------
   SUBROUTINE riverlake_network_final ()

#ifdef USEMPI
   USE MOD_SPMD_Task, ONLY: p_err
#endif
   IMPLICIT NONE

      CALL worker_pushdata_free_mem (push_inpm2ucat)
      CALL worker_pushdata_free_mem (push_ucat2inpm)
      CALL worker_pushdata_free_mem (push_ucat2grid)
      CALL worker_pushdata_free_mem (allreduce_inpm)
      CALL worker_pushdata_free_mem (push_next2ucat)
      CALL worker_pushdata_free_mem (push_ups2ucat)
      CALL worker_pushdata_free_mem (push_bif_dn2pth)
      CALL worker_pushdata_free_mem (push_bif_influx)
      CALL riverlake_pushcache_free_mem (cache_ups2ucat)
      CALL riverlake_pushcache_free_mem (cache_bif_dn2pth)
      CALL riverlake_pushcache_free_mem (cache_bif_influx)
      CALL worker_remapdata_free_mem (remap_patch2inpm)
      CALL grid_free_mem (griducat)

      IF (allocated(x_ucat           )) deallocate(x_ucat           )
      IF (allocated(y_ucat           )) deallocate(y_ucat           )

      IF (allocated(ucat_ucid        )) deallocate(ucat_ucid        )
      IF (allocated(ucat_gdid        )) deallocate(ucat_gdid        )

      IF (allocated(numucat_wrk      )) deallocate(numucat_wrk      )
      IF (allocated(ucat_data_address)) deallocate(ucat_data_address)

      IF (allocated(inpm_gdid        )) deallocate(inpm_gdid        )
      IF (allocated(idmap_gd2uc      )) deallocate(idmap_gd2uc      )
      IF (allocated(area_gd2uc       )) deallocate(area_gd2uc       )
      IF (allocated(idmap_uc2gd      )) deallocate(idmap_uc2gd      )
      IF (allocated(area_uc2gd       )) deallocate(area_uc2gd       )
      IF (allocated(ucat_next        )) deallocate(ucat_next        )
      IF (allocated(ucat_ups         )) deallocate(ucat_ups         )
      IF (allocated(wts_ups          )) deallocate(wts_ups          )
      IF (allocated(irivsys          )) deallocate(irivsys          )
      IF (allocated(irivsys_shared   )) deallocate(irivsys_shared   )
      IF (allocated(lake_type        )) deallocate(lake_type        )

      IF (allocated(topo_rivelv      )) deallocate(topo_rivelv      )
      IF (allocated(topo_rivhgt      )) deallocate(topo_rivhgt      )
      IF (allocated(topo_rivlen      )) deallocate(topo_rivlen      )
      IF (allocated(topo_rivman      )) deallocate(topo_rivman      )
      IF (allocated(topo_rivwth      )) deallocate(topo_rivwth      )
      IF (allocated(topo_rivare      )) deallocate(topo_rivare      )
      IF (allocated(topo_rivstomax   )) deallocate(topo_rivstomax   )
      IF (allocated(topo_area        )) deallocate(topo_area        )
      IF (allocated(topo_fldhgt      )) deallocate(topo_fldhgt      )
      IF (allocated(levee_frc_data   )) deallocate(levee_frc_data   )
      IF (allocated(levee_hgt_data   )) deallocate(levee_hgt_data   )
      IF (allocated(bedelv_next      )) deallocate(bedelv_next      )
      IF (allocated(outletwth        )) deallocate(outletwth        )

      IF (allocated(floodplain_curve )) deallocate(floodplain_curve )

      IF (allocated(allups_mask_ucat )) deallocate(allups_mask_ucat )

      IF (allocated(pth_upst_local    )) deallocate(pth_upst_local    )
      IF (allocated(pth_down_local    )) deallocate(pth_down_local    )
      IF (allocated(pth_down_ucid     )) deallocate(pth_down_ucid     )
      IF (allocated(pth_global_id     )) deallocate(pth_global_id     )
      IF (allocated(pth_dst           )) deallocate(pth_dst           )
      IF (allocated(pth_elv           )) deallocate(pth_elv           )
      IF (allocated(pth_wth           )) deallocate(pth_wth           )
      IF (allocated(pth_man           )) deallocate(pth_man           )
      IF (allocated(bif_incoming_pths )) deallocate(bif_incoming_pths )
      IF (allocated(bif_incoming_wts  )) deallocate(bif_incoming_wts  )

#ifdef USEMPI
      IF (p_comm_rivsys /= MPI_COMM_NULL) THEN
         CALL mpi_comm_free (p_comm_rivsys, p_err)
         p_comm_rivsys = MPI_COMM_NULL
      ENDIF
#endif

   END SUBROUTINE riverlake_network_final

   ! Time step groups for bifurcation (master only), see irivsys_shared.
   ! rivermouth(totalnumucat) holds the river system of each unit catchment on
   ! input and the cluster of river systems connected by bifurcation pathways on
   ! output; a cluster is labelled by one of its river systems.  shared_ucat gives,
   ! for each unit catchment, the index (1..numrivsys_shared) of its cluster when
   ! that cluster has unit catchments on more than one worker, and 0 otherwise.
   SUBROUTINE bifurcation_dt_groups (parafile, numrivmth, rivermouth, shared_ucat)

   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   IMPLICIT NONE

   character(len=*), intent(in) :: parafile
   integer, intent(in)    :: numrivmth
   integer, intent(inout) :: rivermouth (:)
   integer, allocatable, intent(out) :: shared_ucat (:)

   integer, allocatable :: bif_upst (:), bif_down (:)
   integer, allocatable :: parent (:), owner (:), ishared (:)
   integer :: npth, ip, i, a, b, iworker, ic

      CALL ncio_inquire_length (parafile, 'bifurcation_upst', npth)
      IF (npth > 0) THEN
         CALL ncio_read_serial (parafile, 'bifurcation_upst', bif_upst)
         CALL ncio_read_serial (parafile, 'bifurcation_down', bif_down)
      ENDIF

      ! union-find over river systems; a pathway into the ocean (down <= 0) or
      ! out of range joins nothing (read_bifurcation_global_arrays checks ranges)
      allocate (parent (numrivmth))
      parent = (/(i, i = 1, numrivmth)/)
      DO ip = 1, npth
         IF (bif_upst(ip) < 1 .or. bif_upst(ip) > size(rivermouth)) CYCLE
         IF (bif_down(ip) < 1 .or. bif_down(ip) > size(rivermouth)) CYCLE
         a = find_root (rivermouth(bif_upst(ip)))
         b = find_root (rivermouth(bif_down(ip)))
         IF (a /= b) parent(max(a,b)) = min(a,b)
      ENDDO

      DO i = 1, size(rivermouth)
         rivermouth(i) = find_root (rivermouth(i))
      ENDDO

      ! clusters with unit catchments on more than one worker
      allocate (owner (numrivmth));  owner(:) = -1
      allocate (ishared (numrivmth)); ishared(:) = 0
      DO iworker = 0, p_np_worker-1
         DO i = 1, numucat_wrk(iworker)
            ic = rivermouth(ucat_data_address(iworker)%val(i))
            IF (owner(ic) == -1) THEN
               owner(ic) = iworker
            ELSEIF (owner(ic) /= iworker) THEN
               ishared(ic) = 1
            ENDIF
         ENDDO
      ENDDO

      numrivsys_shared = 0
      DO ic = 1, numrivmth
         IF (ishared(ic) > 0) THEN
            numrivsys_shared = numrivsys_shared + 1
            ishared(ic) = numrivsys_shared
         ENDIF
      ENDDO

      allocate (shared_ucat (size(rivermouth)))
      shared_ucat = ishared(rivermouth)

      write(*,'(A,I0,A,I0,A,I0,A)') 'Bifurcation time step groups: ', &
         count(parent == (/(i, i = 1, numrivmth)/)), ' groups from ', numrivmth, &
         ' river systems, ', numrivsys_shared, ' shared by several workers.'

      deallocate (parent, owner, ishared)
      IF (allocated(bif_upst)) deallocate (bif_upst)
      IF (allocated(bif_down)) deallocate (bif_down)

   CONTAINS

      integer FUNCTION find_root (k)
      integer, intent(in) :: k
         find_root = k
         DO WHILE (parent(find_root) /= find_root)
            find_root = parent(find_root)
         ENDDO
      END FUNCTION find_root

   END SUBROUTINE bifurcation_dt_groups

   SUBROUTINE read_and_distribute_bifurcation (parafile)

   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   USE MOD_Utils
   IMPLICIT NONE

   character(len=*), intent(in) :: parafile

   integer,  allocatable :: bif_upst_all  (:)
   integer,  allocatable :: bif_down_all  (:)
   real(r8), allocatable :: bif_dist_all  (:)
   real(r8), allocatable :: bif_elev_all  (:,:)
   real(r8), allocatable :: bif_wdth_all  (:,:)
   real(r8), allocatable :: bif_mann_all  (:)

   integer,  allocatable :: iworker_of_ucat (:)
   integer,  allocatable :: pth_owner       (:)
   integer,  allocatable :: npth_wrk        (:)

   integer,  allocatable :: pth_upst_send (:)
   integer,  allocatable :: pth_down_send (:)
   integer,  allocatable :: pth_glid_send (:)
   real(r8), allocatable :: pth_dist_send (:)
   real(r8), allocatable :: pth_elev_send (:,:)
   real(r8), allocatable :: pth_wdth_send (:,:)

   integer,  allocatable :: bif_inc_cnt   (:)
   integer,  allocatable :: bif_inc_all   (:,:)
   integer,  allocatable :: bif_inc_send  (:,:)
   real(r8), allocatable :: bif_wt_send   (:,:)

   integer :: iworker, nucat, npth, ip, i, j, iloc
   integer :: max_bif_inc_global

#ifdef USEMPI

      IF (p_is_master) THEN

         CALL read_bifurcation_global_arrays (parafile, bif_upst_all, bif_down_all, &
            bif_dist_all, bif_elev_all, bif_wdth_all, bif_mann_all)

         allocate (iworker_of_ucat (totalnumucat))
         iworker_of_ucat(:) = -1
         DO iworker = 0, p_np_worker-1
            DO i = 1, numucat_wrk(iworker)
               iworker_of_ucat(ucat_data_address(iworker)%val(i)) = iworker
            ENDDO
         ENDDO

         allocate (pth_owner (totalnpthout))
         allocate (npth_wrk  (0:p_np_worker-1))
         npth_wrk(:) = 0
         DO ip = 1, totalnpthout
            pth_owner(ip) = iworker_of_ucat(bif_upst_all(ip))
            IF (pth_owner(ip) < 0) CALL CoLM_stop ('bifurcation upstream owner not found')
            npth_wrk(pth_owner(ip)) = npth_wrk(pth_owner(ip)) + 1
         ENDDO

         allocate (bif_inc_cnt (totalnumucat))
         bif_inc_cnt(:) = 0
         DO ip = 1, totalnpthout
            j = bif_down_all(ip)
            IF (j > 0 .and. j <= totalnumucat) THEN
               bif_inc_cnt(j) = bif_inc_cnt(j) + 1
            ENDIF
         ENDDO
         max_bif_inc_global = maxval(bif_inc_cnt)
         IF (max_bif_inc_global < 1) max_bif_inc_global = 1

         allocate (bif_inc_all (max_bif_inc_global, totalnumucat))
         bif_inc_all(:,:) = 0
         bif_inc_cnt(:) = 0
         DO ip = 1, totalnpthout
            j = bif_down_all(ip)
            IF (j > 0 .and. j <= totalnumucat) THEN
               bif_inc_cnt(j) = bif_inc_cnt(j) + 1
               bif_inc_all(bif_inc_cnt(j), j) = ip
            ENDIF
         ENDDO

      ENDIF

      CALL mpi_bcast (totalnpthout, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (npthlev_bif,  1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)

      IF (.not. p_is_master) allocate (pth_man (npthlev_bif))
      IF (p_is_master) THEN
         IF (.not. allocated(pth_man)) allocate (pth_man (npthlev_bif))
         pth_man(:) = bif_mann_all(:)
         deallocate (bif_mann_all)
      ENDIF
      CALL mpi_bcast (pth_man, npthlev_bif, MPI_REAL8, p_address_master, p_comm_glb, p_err)

      IF (p_is_master) max_bif_incoming = max_bif_inc_global
      CALL mpi_bcast (max_bif_incoming, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)

      IF (p_is_master) THEN

         DO iworker = 0, p_np_worker-1

            npth  = npth_wrk(iworker)
            nucat = numucat_wrk(iworker)

            CALL mpi_send (npth, 1, MPI_INTEGER, p_address_worker(iworker), &
               mpi_tag_mesg, p_comm_glb, p_err)

            IF (npth > 0) THEN
               allocate (pth_upst_send (npth))
               allocate (pth_down_send (npth))
               allocate (pth_glid_send (npth))
               allocate (pth_dist_send (npth))
               allocate (pth_elev_send (npthlev_bif, npth))
               allocate (pth_wdth_send (npthlev_bif, npth))

               j = 0
               DO ip = 1, totalnpthout
                  IF (pth_owner(ip) == iworker) THEN
                     j = j + 1
                     pth_upst_send(j) = bif_upst_all(ip)
                     pth_down_send(j) = bif_down_all(ip)
                     pth_glid_send(j) = ip
                     pth_dist_send(j) = bif_dist_all(ip)
                     pth_elev_send(:,j) = bif_elev_all(:,ip)
                     pth_wdth_send(:,j) = bif_wdth_all(:,ip)
                  ENDIF
               ENDDO

               CALL mpi_send (pth_upst_send, npth, MPI_INTEGER, &
                  p_address_worker(iworker), mpi_tag_data, p_comm_glb, p_err)
               CALL mpi_send (pth_down_send, npth, MPI_INTEGER, &
                  p_address_worker(iworker), mpi_tag_data, p_comm_glb, p_err)
               CALL mpi_send (pth_glid_send, npth, MPI_INTEGER, &
                  p_address_worker(iworker), mpi_tag_data, p_comm_glb, p_err)
               CALL mpi_send (pth_dist_send, npth, MPI_REAL8, &
                  p_address_worker(iworker), mpi_tag_data, p_comm_glb, p_err)
               CALL mpi_send (pth_elev_send, npthlev_bif*npth, MPI_REAL8, &
                  p_address_worker(iworker), mpi_tag_data, p_comm_glb, p_err)
               CALL mpi_send (pth_wdth_send, npthlev_bif*npth, MPI_REAL8, &
                  p_address_worker(iworker), mpi_tag_data, p_comm_glb, p_err)

               deallocate (pth_upst_send)
               deallocate (pth_down_send)
               deallocate (pth_glid_send)
               deallocate (pth_dist_send)
               deallocate (pth_elev_send)
               deallocate (pth_wdth_send)
            ENDIF

            IF (nucat > 0) THEN
               allocate (bif_inc_send (max_bif_inc_global, nucat))
               allocate (bif_wt_send  (max_bif_inc_global, nucat))
               DO i = 1, nucat
                  bif_inc_send(:,i) = bif_inc_all(:,ucat_data_address(iworker)%val(i))
               ENDDO
               bif_wt_send(:,:) = 0.
               WHERE (bif_inc_send > 0) bif_wt_send = 1.

               CALL mpi_send (bif_inc_send, max_bif_inc_global*nucat, MPI_INTEGER, &
                  p_address_worker(iworker), mpi_tag_data, p_comm_glb, p_err)
               CALL mpi_send (bif_wt_send,  max_bif_inc_global*nucat, MPI_REAL8, &
                  p_address_worker(iworker), mpi_tag_data, p_comm_glb, p_err)

               deallocate (bif_inc_send)
               deallocate (bif_wt_send)
            ENDIF

         ENDDO

         deallocate (bif_upst_all)
         deallocate (bif_down_all)
         deallocate (bif_dist_all)
         deallocate (bif_elev_all)
         deallocate (bif_wdth_all)
         deallocate (iworker_of_ucat)
         deallocate (pth_owner)
         deallocate (npth_wrk)
         deallocate (bif_inc_cnt)
         deallocate (bif_inc_all)

      ELSEIF (p_is_worker) THEN

         CALL mpi_recv (npthout_local, 1, MPI_INTEGER, p_address_master, &
            mpi_tag_mesg, p_comm_glb, p_stat, p_err)

         IF (npthout_local > 0) THEN
            allocate (pth_upst_local (npthout_local))
            allocate (pth_down_ucid  (npthout_local))
            allocate (pth_global_id  (npthout_local))
            allocate (pth_dst        (npthout_local))
            allocate (pth_elv        (npthlev_bif, npthout_local))
            allocate (pth_wth        (npthlev_bif, npthout_local))

            CALL mpi_recv (pth_upst_local, npthout_local, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)
            CALL mpi_recv (pth_down_ucid,  npthout_local, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)
            CALL mpi_recv (pth_global_id,  npthout_local, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)
            CALL mpi_recv (pth_dst,        npthout_local, MPI_REAL8,   p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)
            CALL mpi_recv (pth_elv,        npthlev_bif*npthout_local, MPI_REAL8, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)
            CALL mpi_recv (pth_wth,        npthlev_bif*npthout_local, MPI_REAL8, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (pth_down_local (npthout_local))
            CALL localize_bifurcation_path_indices (pth_upst_local, pth_down_ucid, pth_down_local)
         ELSE
            npthout_local = 0
            allocate (pth_upst_local (0))
            allocate (pth_down_ucid  (0))
            allocate (pth_down_local (0))
            allocate (pth_global_id  (0))
            allocate (pth_dst        (0))
            allocate (pth_elv        (npthlev_bif, 0))
            allocate (pth_wth        (npthlev_bif, 0))
         ENDIF

         IF (numucat > 0) THEN
            allocate (bif_incoming_pths (max_bif_incoming, numucat))
            allocate (bif_incoming_wts  (max_bif_incoming, numucat))
            CALL mpi_recv (bif_incoming_pths, max_bif_incoming*numucat, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)
            CALL mpi_recv (bif_incoming_wts,  max_bif_incoming*numucat, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)
         ELSE
            allocate (bif_incoming_pths (max_bif_incoming, 0))
            allocate (bif_incoming_wts  (max_bif_incoming, 0))
         ENDIF

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_io) THEN
         npthout_local = 0
         allocate (pth_upst_local (0))
         allocate (pth_down_ucid  (0))
         allocate (pth_down_local (0))
         allocate (pth_global_id  (0))
         allocate (pth_dst        (0))
         allocate (pth_elv        (npthlev_bif, 0))
         allocate (pth_wth        (npthlev_bif, 0))
         allocate (bif_incoming_pths (max_bif_incoming, 0))
         allocate (bif_incoming_wts  (max_bif_incoming, 0))
      ENDIF

#else

      CALL read_bifurcation_global_arrays (parafile, bif_upst_all, bif_down_all, &
         bif_dist_all, bif_elev_all, bif_wdth_all, bif_mann_all)

      npthout_local = totalnpthout

      allocate (pth_man (npthlev_bif))
      pth_man(:) = bif_mann_all(:)
      deallocate (bif_mann_all)

      IF (npthout_local > 0) THEN
         allocate (pth_upst_local (npthout_local))
         allocate (pth_down_ucid  (npthout_local))
         allocate (pth_down_local (npthout_local))
         allocate (pth_global_id  (npthout_local))
         allocate (pth_dst        (npthout_local))
         allocate (pth_elv        (npthlev_bif, npthout_local))
         allocate (pth_wth        (npthlev_bif, npthout_local))

         pth_down_ucid(:) = bif_down_all(:)
         pth_dst(:)       = bif_dist_all(:)
         pth_elv(:,:)     = bif_elev_all(:,:)
         pth_wth(:,:)     = bif_wdth_all(:,:)

         DO ip = 1, npthout_local
            pth_global_id(ip) = ip
         ENDDO

         pth_upst_local(:) = bif_upst_all(:)
         CALL localize_bifurcation_path_indices (pth_upst_local, pth_down_ucid, pth_down_local)
      ELSE
         npthout_local = 0
         allocate (pth_upst_local (0))
         allocate (pth_down_ucid  (0))
         allocate (pth_down_local (0))
         allocate (pth_global_id  (0))
         allocate (pth_dst        (0))
         allocate (pth_elv        (npthlev_bif, 0))
         allocate (pth_wth        (npthlev_bif, 0))
      ENDIF

      deallocate (bif_upst_all)
      deallocate (bif_down_all)
      deallocate (bif_dist_all)
      deallocate (bif_elev_all)
      deallocate (bif_wdth_all)

      allocate (bif_inc_cnt (totalnumucat))
      bif_inc_cnt(:) = 0
      DO ip = 1, totalnpthout
         j = pth_down_ucid(ip)
         IF (j > 0 .and. j <= totalnumucat) THEN
            bif_inc_cnt(j) = bif_inc_cnt(j) + 1
         ENDIF
      ENDDO
      max_bif_incoming = maxval(bif_inc_cnt)
      IF (max_bif_incoming < 1) max_bif_incoming = 1

      allocate (bif_incoming_pths (max_bif_incoming, numucat))
      allocate (bif_incoming_wts  (max_bif_incoming, numucat))
      bif_incoming_pths(:,:) = 0
      bif_incoming_wts(:,:)  = 0.

      bif_inc_cnt(:) = 0
      DO ip = 1, totalnpthout
         j = pth_down_ucid(ip)
         IF (j > 0 .and. j <= totalnumucat) THEN
            bif_inc_cnt(j) = bif_inc_cnt(j) + 1
            i = find_in_sorted_list1 (j, numucat, ucat_ucid)
            IF (i > 0) THEN
               bif_incoming_pths(bif_inc_cnt(j), i) = ip
               bif_incoming_wts (bif_inc_cnt(j), i) = 1.
            ENDIF
         ENDIF
      ENDDO

      deallocate (bif_inc_cnt)

#endif


      CALL build_worker_pushdata (numucat, ucat_ucid, npthout_local, pth_down_ucid, push_bif_dn2pth)

      CALL build_worker_pushdata (npthout_local, pth_global_id, numucat, &
         bif_incoming_pths, bif_incoming_wts, push_bif_influx)

      CALL riverlake_pushcache_build (push_bif_dn2pth, cache_bif_dn2pth)
      CALL riverlake_pushcache_build (push_bif_influx, cache_bif_influx)

   END SUBROUTINE read_and_distribute_bifurcation

   SUBROUTINE read_bifurcation_global_arrays (parafile, bif_upst_all, bif_down_all, &
      bif_dist_all, bif_elev_all, bif_wdth_all, bif_mann_all)

   USE MOD_NetCDFSerial
   USE MOD_Utils
   IMPLICIT NONE

   character(len=*), intent(in) :: parafile
   integer,  allocatable, intent(out) :: bif_upst_all  (:)
   integer,  allocatable, intent(out) :: bif_down_all  (:)
   real(r8), allocatable, intent(out) :: bif_dist_all  (:)
   real(r8), allocatable, intent(out) :: bif_elev_all  (:,:)
   real(r8), allocatable, intent(out) :: bif_wdth_all  (:,:)
   real(r8), allocatable, intent(out) :: bif_mann_all  (:)
   integer :: ip, ilev, prev_active_lev

      CALL ncio_inquire_length (parafile, 'bifurcation_upst',    totalnpthout)
      CALL ncio_inquire_length (parafile, 'bifurcation_manning', npthlev_bif)

      CALL ncio_read_serial (parafile, 'bifurcation_upst',      bif_upst_all)
      CALL ncio_read_serial (parafile, 'bifurcation_down',      bif_down_all)
      CALL ncio_read_serial (parafile, 'bifurcation_distance',  bif_dist_all)
      CALL ncio_read_serial (parafile, 'bifurcation_elevation', bif_elev_all)
      CALL ncio_read_serial (parafile, 'bifurcation_width',     bif_wdth_all)
      CALL ncio_read_serial (parafile, 'bifurcation_manning',   bif_mann_all)

      IF (size(bif_down_all) /= totalnpthout .or. &
          size(bif_dist_all) /= totalnpthout .or. &
          size(bif_elev_all, 1) /= npthlev_bif .or. &
          size(bif_elev_all, 2) /= totalnpthout .or. &
          size(bif_wdth_all, 1) /= npthlev_bif .or. &
          size(bif_wdth_all, 2) /= totalnpthout .or. &
          size(bif_mann_all) /= npthlev_bif) THEN
         CALL CoLM_stop ('bifurcation parameter dimensions are inconsistent')
      ENDIF

      IF (any(.not. ieee_is_finite(bif_dist_all))) THEN
         CALL CoLM_stop ('bifurcation distance must be finite and positive')
      ELSEIF (any(bif_dist_all <= 0._r8)) THEN
         CALL CoLM_stop ('bifurcation distance must be finite and positive')
      ENDIF
      IF (any(.not. ieee_is_finite(bif_elev_all))) THEN
         CALL CoLM_stop ('bifurcation elevation contains a non-finite value')
      ENDIF
      IF (any(.not. ieee_is_finite(bif_wdth_all))) THEN
         CALL CoLM_stop ('bifurcation width contains a non-finite value')
      ELSEIF (any(bif_wdth_all < 0._r8)) THEN
         CALL CoLM_stop ('bifurcation width must be non-negative')
      ENDIF
      IF (any(.not. ieee_is_finite(bif_mann_all))) THEN
         CALL CoLM_stop ('bifurcation Manning coefficient contains a non-finite value')
      ENDIF

      DO ilev = 1, npthlev_bif
         IF (any(bif_wdth_all(ilev, :) > 0._r8) .and. bif_mann_all(ilev) <= 0._r8) THEN
            CALL CoLM_stop ('active bifurcation layer requires a positive Manning coefficient')
         ENDIF
      ENDDO

      DO ip = 1, totalnpthout
         IF (bif_upst_all(ip) < 1 .or. bif_upst_all(ip) > totalnumucat) CALL CoLM_stop ( &
            'bifurcation upstream index out of range')
         IF (bif_down_all(ip) > totalnumucat) CALL CoLM_stop ('bifurcation downstream index out of range')
         IF (bif_down_all(ip) == bif_upst_all(ip)) CALL CoLM_stop ('bifurcation self-loop pathway is invalid')
         IF (.not. any(bif_wdth_all(:, ip) > 0._r8)) THEN
            CALL CoLM_stop ('bifurcation pathway has no active positive-width layer')
         ENDIF
         prev_active_lev = 0
         DO ilev = 1, npthlev_bif
            IF (bif_wdth_all(ilev, ip) <= 0._r8) CYCLE
            IF (prev_active_lev > 0) THEN
               IF (bif_elev_all(ilev, ip) < bif_elev_all(prev_active_lev, ip)) THEN
                  CALL CoLM_stop ( &
                     'active bifurcation layer elevation must be non-decreasing')
               ENDIF
            ENDIF
            prev_active_lev = ilev
         ENDDO
      ENDDO

   END SUBROUTINE read_bifurcation_global_arrays

   SUBROUTINE localize_bifurcation_path_indices (path_upst, path_down_ucid, path_down_local)

   USE MOD_Utils, only: find_in_sorted_list1
   IMPLICIT NONE

   integer, intent(inout) :: path_upst(:)
   integer, intent(in)    :: path_down_ucid(:)
   integer, intent(out)   :: path_down_local(:)
   integer :: ip, i

      DO ip = 1, size(path_upst)
         i = find_in_sorted_list1 (path_upst(ip), numucat, ucat_ucid)
         IF (i > 0) path_upst(ip) = i
      ENDDO

      path_down_local = -1
      DO ip = 1, size(path_down_ucid)
         i = find_in_sorted_list1 (path_down_ucid(ip), numucat, ucat_ucid)
         IF (i > 0) path_down_local(ip) = i
      ENDDO

   END SUBROUTINE localize_bifurcation_path_indices

END MODULE MOD_Grid_RiverLakeNetwork
#endif
