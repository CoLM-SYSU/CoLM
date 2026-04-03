#include <define.h>

#ifdef GridRiverLakeFlow
MODULE MOD_Grid_RiverLakeNetwork
!--------------------------------------------------------------------------------
! DESCRIPTION:
!--------------------------------------------------------------------------------

   USE MOD_Grid
   USE MOD_WorkerPushData
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

   ! ----- Part 3: river systems -----
   integer :: numrivsys
   logical :: rivsys_by_multiple_procs
   integer, allocatable :: irivsys (:)
#ifdef USEMPI
   integer :: p_comm_rivsys
#endif


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

         parafile = DEF_UnitCatchment_file

         CALL ncio_read_serial (parafile, 'seq_x', x_ucat)
         CALL ncio_read_serial (parafile, 'seq_y', y_ucat)

         CALL ncio_read_serial (parafile, 'seq_next', ucat_next)

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
#endif

      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            allocate (wts_ups (upnmax,numucat))
            wts_ups(:,:) = 1.
         ENDIF
      ENDIF

      CALL build_worker_pushdata (numucat, ucat_ucid, numucat, ucat_next, push_next2ucat)
      CALL build_worker_pushdata (numucat, ucat_ucid, numucat, ucat_ups,  wts_ups, push_ups2ucat )

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

      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            color = maxval(rivermouth)
            CALL mpi_comm_split (p_comm_worker, color, p_iam_worker, p_comm_rivsys, p_err)
         ELSE
            CALL mpi_comm_split (p_comm_worker, MPI_UNDEFINED, p_iam_worker, p_comm_rivsys, p_err)
         ENDIF

         rivsys_by_multiple_procs = .false.
         IF (p_comm_rivsys /= MPI_COMM_NULL) THEN
            CALL mpi_comm_size (p_comm_rivsys, p_np_rivsys, p_err)
            IF (p_np_rivsys > 1) THEN
               rivsys_by_multiple_procs = .true.
            ENDIF
         ENDIF
      ENDIF
#else
      rivsys_by_multiple_procs = .false.
#endif

      IF (p_is_worker) THEN

         IF (numucat > 0) allocate (irivsys (numucat))

         IF (.not. rivsys_by_multiple_procs) THEN
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

            ENDIF
         ELSE
            numrivsys  = 1
            irivsys(:) = 1
         ENDIF

      ENDIF

      IF (allocated(rivermouth)) deallocate(rivermouth)
      IF (allocated(order_ucat)) deallocate(order_ucat)

      ! ----- Parameters for River and Lake -----

      CALL readin_riverlake_parameter (parafile, 'topo_rivelv',    rdata1d = topo_rivelv   )
      CALL readin_riverlake_parameter (parafile, 'topo_rivhgt',    rdata1d = topo_rivhgt   )
      CALL readin_riverlake_parameter (parafile, 'topo_rivlen',    rdata1d = topo_rivlen   )
      CALL readin_riverlake_parameter (parafile, 'topo_rivman',    rdata1d = topo_rivman   )
      CALL readin_riverlake_parameter (parafile, 'topo_rivwth',    rdata1d = topo_rivwth   )
      CALL readin_riverlake_parameter (parafile, 'topo_rivstomax', rdata1d = topo_rivstomax)
      CALL readin_riverlake_parameter (parafile, 'topo_area',      rdata1d = topo_area     )
      CALL readin_riverlake_parameter (parafile, 'topo_fldhgt',    rdata2d = topo_fldhgt   )

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

               floodplain_curve(i)%flpstomax(0) = 0.
               DO j = 1, floodplain_curve(i)%nlfp
                  floodplain_curve(i)%flpstomax(j) = floodplain_curve(i)%flpstomax(j-1)        &
                     + 0.5 * (floodplain_curve(i)%flparea(j) + floodplain_curve(i)%flparea(j-1)) &
                           * (floodplain_curve(i)%flphgt(j)  - floodplain_curve(i)%flphgt(j-1))
               ENDDO
            ENDDO

            allocate (bedelv_next (numucat))
            allocate (outletwth   (numucat))

         ENDIF
      ENDIF

      CALL worker_push_data (push_next2ucat, topo_rivelv, bedelv_next, fillvalue = spval)
      CALL worker_push_data (push_next2ucat, topo_rivwth, outletwth  , fillvalue = spval)

      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            WHERE (ucat_next > 0)
               outletwth = (outletwth + topo_rivwth) * 0.5
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


      IF (allocated (uc_up2down   )) deallocate (uc_up2down   )
      IF (allocated (nups_nst     )) deallocate (nups_nst     )
      IF (allocated (iups_nst     )) deallocate (iups_nst     )
      IF (allocated (ucat_area_all)) deallocate (ucat_area_all)

   END SUBROUTINE build_riverlake_network

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
         ENDIF

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

   IMPLICIT NONE

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
      IF (allocated(irivsys          )) deallocate(irivsys          )

      IF (allocated(topo_rivelv      )) deallocate(topo_rivelv      )
      IF (allocated(topo_rivhgt      )) deallocate(topo_rivhgt      )
      IF (allocated(topo_rivlen      )) deallocate(topo_rivlen      )
      IF (allocated(topo_rivman      )) deallocate(topo_rivman      )
      IF (allocated(topo_rivwth      )) deallocate(topo_rivwth      )
      IF (allocated(topo_rivare      )) deallocate(topo_rivare      )
      IF (allocated(topo_rivstomax   )) deallocate(topo_rivstomax   )
      IF (allocated(topo_area        )) deallocate(topo_area        )
      IF (allocated(topo_fldhgt      )) deallocate(topo_fldhgt      )
      IF (allocated(bedelv_next      )) deallocate(bedelv_next      )
      IF (allocated(outletwth        )) deallocate(outletwth        )

      IF (allocated(floodplain_curve )) deallocate(floodplain_curve )

      IF (allocated(allups_mask_ucat )) deallocate(allups_mask_ucat )

   END SUBROUTINE riverlake_network_final

END MODULE MOD_Grid_RiverLakeNetwork
#endif
