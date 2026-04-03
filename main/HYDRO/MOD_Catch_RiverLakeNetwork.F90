#include <define.h>

#ifdef CatchLateralFlow
MODULE MOD_Catch_RiverLakeNetwork
!--------------------------------------------------------------------------------
! DESCRIPTION:
!
!    River networks: data and communication subroutines.
!
! Created by Shupeng Zhang, May 2023
!--------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Vars_Global, only: spval
   USE MOD_Pixelset
   USE MOD_Catch_BasinNetwork
   USE MOD_Catch_HillslopeNetwork
   USE MOD_Catch_Vars_TimeVariables
   IMPLICIT NONE

   ! -- river parameters --
   real(r8), allocatable :: riverlen  (:)
   real(r8), allocatable :: riverelv  (:)
   real(r8), allocatable :: riverwth  (:)
   real(r8), allocatable :: riverarea (:)
   real(r8), allocatable :: riverdpth (:)

   real(r8), allocatable :: basinelv  (:)
   real(r8), allocatable :: bedelv    (:)
   real(r8), allocatable :: handmin   (:)

   real(r8), allocatable :: wtsrfelv  (:)

   ! index of downstream river
   ! > 0 : other catchment;   0 : river mouth; -1 : inland depression
   integer, allocatable :: riverdown  (:)
   integer, allocatable :: ilocdown   (:)
   logical, allocatable :: to_lake    (:)

   real(r8), allocatable :: riverlen_ds (:)
   real(r8), allocatable :: wtsrfelv_ds (:)
   real(r8), allocatable :: riverwth_ds (:)
   real(r8), allocatable :: bedelv_ds   (:)

   real(r8), allocatable :: outletwth (:)

   integer :: riversystem

   integer :: numrivsys
   integer, allocatable :: irivsys (:)

#ifdef USEMPI
   integer :: p_comm_rivsys

   integer :: numbsnlink
   integer, allocatable :: linkbindex (:)

   integer :: nlink_me
   integer, allocatable :: linkpush (:)
   integer, allocatable :: linkpull (:)
#endif

   ! -- lake data type --
   type :: lake_info_type
      integer :: nsub
      real(r8), allocatable :: area0  (:) ! area data in HRU order
      real(r8), allocatable :: area   (:) ! area data in the order from deepest to shallowest HRU
      real(r8), allocatable :: depth0 (:) ! depth data in HRU order
      real(r8), allocatable :: depth  (:) ! depth data in the order from deepest to shallowest HRU
      ! a curve describing the relationship between depth of water from lake bottom and total water volume
      ! the i-th value corresponds to the volume when water depth is at i-th depth
      real(r8), allocatable :: dep_vol_curve (:)
   CONTAINS
      procedure, PUBLIC :: surface => retrieve_lake_surface_from_volume
      procedure, PUBLIC :: volume  => retrieve_lake_volume_from_surface
      final :: lake_info_free_mem
   END type lake_info_type

   ! -- lake information --
   type(lake_info_type), allocatable :: lakeinfo (:)

   ! -- information of HRU in basin --
   type(hillslope_network_type), pointer :: hillslope_basin (:)


   ! ----- subroutines and functions -----
   PUBLIC :: river_lake_network_init
   PUBLIC :: pull_from_downstream
   PUBLIC :: push_to_downstream
   PUBLIC :: calc_riverdepth_from_runoff
   PUBLIC :: retrieve_lake_surface_from_volume
   PUBLIC :: river_lake_network_final

CONTAINS

   ! ----------
   SUBROUTINE river_lake_network_init (patcharea)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_Mesh
   USE MOD_Pixel
   USE MOD_LandElm
   USE MOD_LandHRU
   USE MOD_LandPatch
   USE MOD_ElementNeighbour
   USE MOD_DataType
   USE MOD_Utils
   USE MOD_UserDefFun
   USE MOD_Vars_TimeInvariants, only: lakedepth
   IMPLICIT NONE

   real(r8), intent(in) :: patcharea (:)

   ! Local Variables
   character(len=256) :: river_file, rivdpt_file
   logical :: use_calc_rivdpt

   integer :: totalnumbasin, ibasin, jbasin, inb, iloc, ielm, i, j, ilink
   integer :: iworker, mesg(2), isrc, idest, nrecv

#ifdef USEMPI
   integer :: nblink_all
   integer,  allocatable :: linkbindex_all (:)
   integer,  allocatable :: linkrivmth_all (:)
#endif

   logical , allocatable :: is_link    (:)
   logical , allocatable :: link_on_me (:)

   integer :: nnode
   integer, allocatable :: route(:)

   integer , allocatable :: icache (:)
   real(r8), allocatable :: rcache (:)
   logical , allocatable :: lcache (:)

   integer, allocatable :: bindex(:), addrbasin(:), addrdown(:)
   integer, allocatable :: basin_sorted(:), basin_order(:), order (:)

   logical, allocatable :: bsnfilter(:)

   ! for lakes
   integer :: ps, pe, nsublake, hs, he, ihru
   integer,  allocatable :: all_lake_id    (:), lake_id_elm     (:)
   integer , allocatable :: lakedown_id_elm(:), lakedown_id_bsn (:)
   real(r8), allocatable :: lakedepth_hru  (:), lakedepth_bsnhru(:)
   real(r8), allocatable :: lakeoutlet_elm (:), lakeoutlet_bsn  (:)
   real(r8), allocatable :: unitarea_hru   (:), unitarea_bsnhru (:)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      use_calc_rivdpt = DEF_USE_EstimatedRiverDepth
      river_file      = DEF_CatchmentMesh_data

      ! step 1: read in parameters from file.
      IF (p_is_master) THEN

         CALL ncio_read_serial (river_file, 'lake_id', all_lake_id)

         CALL ncio_read_serial (river_file, 'basin_downstream', riverdown)
         CALL ncio_read_serial (river_file, 'river_length'    , riverlen )
         CALL ncio_read_serial (river_file, 'river_elevation' , riverelv )
         CALL ncio_read_serial (river_file, 'basin_elevation' , basinelv )

         IF (.not. use_calc_rivdpt) THEN
            CALL ncio_read_serial (river_file, 'river_depth' , riverdpth)
         ENDIF

         riverlen = riverlen * 1.e3 ! km to m

         totalnumbasin = size(riverdown)
         allocate (to_lake (totalnumbasin))
         to_lake = .false.
         DO i = 1, totalnumbasin
            IF (riverdown(i) > 0) THEN
               to_lake(i) = all_lake_id(riverdown(i)) > 0
            ENDIF
         ENDDO

      ENDIF

      ! step 2: Estimate river depth by using runoff data.
      IF (use_calc_rivdpt) THEN
         CALL calc_riverdepth_from_runoff (patcharea)
      ENDIF

#ifdef USEMPI
      IF (p_is_master) THEN

         DO iworker = 0, p_np_worker-1

            CALL mpi_recv (mesg(1:2), 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
            isrc  = mesg(1)
            nrecv = mesg(2)

            IF (nrecv > 0) THEN

               allocate (bindex (nrecv))
               allocate (icache (nrecv))
               allocate (rcache (nrecv))
               allocate (lcache (nrecv))

               CALL mpi_recv (bindex, nrecv, MPI_INTEGER, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               idest = isrc

               icache = riverdown(bindex)
               CALL mpi_send (icache, nrecv, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

               lcache = to_lake(bindex)
               CALL mpi_send (lcache, nrecv, MPI_LOGICAL, idest, mpi_tag_data, p_comm_glb, p_err)

               rcache = riverlen(bindex)
               CALL mpi_send (rcache, nrecv, MPI_REAL8, idest, mpi_tag_data, p_comm_glb, p_err)

               rcache = riverelv(bindex)
               CALL mpi_send (rcache, nrecv, MPI_REAL8, idest, mpi_tag_data, p_comm_glb, p_err)

               rcache = riverdpth(bindex)
               CALL mpi_send (rcache, nrecv, MPI_REAL8, idest, mpi_tag_data, p_comm_glb, p_err)

               rcache = basinelv(bindex)
               CALL mpi_send (rcache, nrecv, MPI_REAL8, idest, mpi_tag_data, p_comm_glb, p_err)

               deallocate (bindex)
               deallocate (icache)
               deallocate (rcache)
               deallocate (lcache)

            ENDIF
         ENDDO

         deallocate (all_lake_id)

      ENDIF

      IF (p_is_worker) THEN

         mesg(1:2) = (/p_iam_glb, numbasin/)
         CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)

         IF (numbasin > 0) THEN
            CALL mpi_send (basinindex, numbasin, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_err)

            allocate (riverdown (numbasin))
            CALL mpi_recv (riverdown, numbasin, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (to_lake (numbasin))
            CALL mpi_recv (to_lake, numbasin, MPI_LOGICAL, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (riverlen (numbasin))
            CALL mpi_recv (riverlen, numbasin, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (riverelv (numbasin))
            CALL mpi_recv (riverelv, numbasin, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (riverdpth (numbasin))
            CALL mpi_recv (riverdpth, numbasin, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (basinelv (numbasin))
            CALL mpi_recv (basinelv, numbasin, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      IF (numbasin > 0) THEN

         riverdown = riverdown(basinindex)
         to_lake   = to_lake  (basinindex)
         riverlen  = riverlen (basinindex)
         riverelv  = riverelv (basinindex)
         riverdpth = riverdpth(basinindex)
         basinelv  = basinelv (basinindex)

      ENDIF
#endif

#ifdef USEMPI
      ! get address of basins
      IF (p_is_master) THEN

         allocate (addrbasin (totalnumbasin)); addrbasin(:) = -1

         DO iworker = 0, p_np_worker-1

            CALL mpi_recv (mesg(1:2), 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
            isrc  = mesg(1)
            nrecv = mesg(2)

            IF (nrecv > 0) THEN
               allocate (bindex (nrecv))

               CALL mpi_recv (bindex, nrecv, MPI_INTEGER, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               addrbasin(bindex) = isrc

               deallocate(bindex)
            ENDIF

         ENDDO

      ELSEIF (p_is_worker) THEN

         mesg(1:2) = (/p_iam_glb, numbasin/)
         CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)

         IF (numbasin > 0) THEN
            CALL mpi_send (basinindex, numbasin, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_err)
         ENDIF

      ENDIF


      IF (p_is_master) THEN

         allocate (addrdown (totalnumbasin))
         allocate (is_link  (totalnumbasin)); is_link = .false.

         DO ibasin = 1, totalnumbasin
            IF (riverdown(ibasin) >= 1) THEN
               addrdown(ibasin) = addrbasin(riverdown(ibasin))
               IF (addrdown(ibasin) /= addrbasin(ibasin)) THEN
                  is_link(riverdown(ibasin)) = .true.
               ENDIF
            ELSE
               addrdown(ibasin) = addrbasin(ibasin)
            ENDIF
         ENDDO

         nblink_all = count(is_link)

         IF (nblink_all > 0) THEN
            allocate (linkbindex_all (nblink_all))
            allocate (linkrivmth_all (nblink_all))
            linkbindex_all = pack((/(ibasin, ibasin = 1, totalnumbasin)/), mask = is_link)
            linkrivmth_all = pack(rivermouth, mask = is_link)
         ENDIF

         deallocate (addrdown)
         deallocate (is_link )

         write(*,'(/,A,I5,A,/)') 'There are ', nblink_all, ' links between processors.'

      ENDIF

      CALL mpi_bcast (nblink_all, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      IF (nblink_all > 0) THEN
         IF (.not. allocated(linkbindex_all))  allocate (linkbindex_all (nblink_all))
         IF (.not. allocated(linkrivmth_all))  allocate (linkrivmth_all (nblink_all))
         CALL mpi_bcast (linkbindex_all, nblink_all, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
         CALL mpi_bcast (linkrivmth_all, nblink_all, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      ENDIF
#endif

      IF (p_is_worker) THEN

         IF (numbasin > 0) THEN
            allocate (basin_sorted (numbasin))
            allocate (basin_order  (numbasin))
            basin_sorted = basinindex
            basin_order  = (/(ibasin, ibasin = 1, numbasin)/)

            CALL quicksort (numbasin, basin_sorted, basin_order)
         ENDIF

         riversystem = -1

#ifdef USEMPI
         IF (nblink_all > 0) THEN
            DO ibasin = 1, numbasin

               iloc = find_in_sorted_list1 (riverdown(ibasin), nblink_all, linkbindex_all)
               IF (iloc <= 0) THEN
                  iloc = find_in_sorted_list1 (basinindex(ibasin), nblink_all, linkbindex_all)
               ENDIF

               IF (iloc > 0) THEN
                  IF (riversystem == -1) THEN
                     riversystem = linkrivmth_all(iloc)
                  ELSEIF (riversystem /= linkrivmth_all(iloc)) THEN
                     write(*,*) 'Warning: river system allocation error!'
                  ENDIF
               ENDIF
            ENDDO
         ENDIF

         IF (riversystem /= -1) THEN
            numbsnlink = count(linkrivmth_all == riversystem)
            allocate (linkbindex (numbsnlink))
            linkbindex = pack(linkbindex_all, linkrivmth_all == riversystem)

            allocate (link_on_me (numbsnlink)); link_on_me = .false.

            DO ibasin = 1, numbsnlink
               iloc = find_in_sorted_list1 (linkbindex(ibasin), numbasin, basin_sorted)
               IF (iloc > 0) THEN
                  link_on_me(ibasin) = .true.
               ENDIF
            ENDDO

            nlink_me = count(link_on_me)

            IF (nlink_me > 0) THEN
               allocate (linkpush (nlink_me))
               allocate (linkpull (nlink_me))
               ilink = 0
               DO ibasin = 1, numbsnlink
                  IF (link_on_me(ibasin)) THEN
                     ilink = ilink + 1
                     iloc = find_in_sorted_list1 (linkbindex(ibasin), numbasin, basin_sorted)
                     linkpush(ilink) = basin_order(iloc)
                     linkpull(ilink) = ibasin
                  ENDIF
               ENDDO
            ENDIF

            deallocate (link_on_me)
         ENDIF

         IF (riversystem /= -1) THEN
            CALL mpi_comm_split (p_comm_worker, riversystem, p_iam_worker, p_comm_rivsys, p_err)
         ELSE
            CALL mpi_comm_split (p_comm_worker, MPI_UNDEFINED, p_iam_worker, p_comm_rivsys, p_err)
         ENDIF
#endif


         IF (numbasin > 0) THEN

            allocate (ilocdown (numbasin)); ilocdown(:) = 0

            DO ibasin = 1, numbasin
               IF (riverdown(ibasin) > 0) THEN
                  iloc = find_in_sorted_list1 (riverdown(ibasin), numbasin, basin_sorted)
                  IF (iloc > 0) THEN
                     ilocdown(ibasin) = basin_order(iloc)
#ifdef USEMPI
                  ELSE
                     iloc = find_in_sorted_list1 (riverdown(ibasin), numbsnlink, linkbindex)
                     ilocdown(ibasin) = - iloc
#endif
                  ENDIF
               ENDIF
            ENDDO
         ENDIF


         IF (numbasin > 0) allocate (irivsys (numbasin))
         IF (numbasin > 0) allocate (route   (numbasin))

         IF (riversystem == -1) THEN
            irivsys(:) = -1
            numrivsys = 0
            DO ibasin = 1, numbasin
               IF (irivsys(ibasin) == -1) THEN

                  jbasin = ibasin
                  nnode = 1
                  route(nnode) = jbasin
                  DO WHILE ((riverdown(jbasin) > 0) .and. (irivsys(jbasin) == -1))
                     nnode = nnode + 1
                     route(nnode) = jbasin
                     jbasin = ilocdown(jbasin)
                  ENDDO

                  IF (irivsys(jbasin) == -1) THEN
                     numrivsys = numrivsys + 1
                     irivsys(jbasin) = numrivsys
                  ENDIF

                  irivsys(route(1:nnode)) = irivsys(jbasin)

               ENDIF
            ENDDO
         ELSE
            numrivsys  = 1
            irivsys(:) = 1
         ENDIF

         IF (numbasin > 0) deallocate (route)

      ENDIF

#ifdef USEMPI
      IF (allocated (linkbindex_all)) deallocate (linkbindex_all)
      IF (allocated (linkrivmth_all)) deallocate (linkrivmth_all)
#endif


      CALL hillslope_network_init (numbasin, basinindex, hillslope_basin)

      IF (p_is_worker) THEN

         IF (numelm > 0)    allocate (lake_id_elm     (numelm))
         IF (numhru > 0)    allocate (lakedepth_hru   (numhru))
         IF (numbsnhru > 0) allocate (lakedepth_bsnhru(numbsnhru))
         IF (numelm > 0)    allocate (lakedown_id_elm (numelm))
         IF (numbasin > 0)  allocate (lakedown_id_bsn (numbasin))
         IF (numelm > 0)    allocate (lakeoutlet_elm  (numelm))
         IF (numbasin > 0)  allocate (lakeoutlet_bsn  (numbasin))
         IF (numhru > 0)    allocate (unitarea_hru    (numhru))
         IF (numbsnhru > 0) allocate (unitarea_bsnhru (numbsnhru))

         DO ibasin = 1, numbasin

            lakedown_id_bsn(ibasin) = 0

            IF ((lake_id(ibasin) /= 0) .and. (to_lake(ibasin))) THEN
               ! lake to lake .or. lake catchment to lake
               lakedown_id_bsn(ibasin) = riverdown(ibasin)
            ENDIF
            IF ((lake_id(ibasin) > 0) .and. (riverdown(ibasin) == 0)) THEN
               ! lake to ocean
               lakedown_id_bsn(ibasin) = -9    ! -9 is ocean
            ENDIF
         ENDDO

      ENDIF

      CALL worker_push_data (push_bsn2elm, lake_id, lake_id_elm, -9999)
      CALL worker_push_data (push_bsn2elm, lakedown_id_bsn, lakedown_id_elm, -9999)

      IF (p_is_worker) THEN

         unitarea_hru   = 0.
         lakedepth_hru  = 0.
         lakeoutlet_elm = 0.

         DO ielm = 1, numelm
            hs = elm_hru%substt(ielm)
            he = elm_hru%subend(ielm)
            DO ihru = hs, he
               ps = hru_patch%substt(ihru)
               pe = hru_patch%subend(ihru)

               unitarea_hru(ihru) = sum(patcharea(ps:pe))

               IF (lake_id_elm(ielm) > 0) THEN
                  lakedepth_hru(ihru) = maxval(lakedepth(ps:pe))
               ENDIF
            ENDDO

            IF (lakedown_id_elm(ielm) /= 0) THEN
               inb = findloc_ud(elementneighbour(ielm)%glbindex == lakedown_id_elm(ielm))
               IF (inb > 0) lakeoutlet_elm(ielm) = elementneighbour(ielm)%lenbdr(inb)
            ENDIF
         ENDDO
      ENDIF

      CALL worker_push_data (push_elm2bsn, lakeoutlet_elm, lakeoutlet_bsn, spval)

      CALL worker_push_data (push_elmhru2bsnhru, unitarea_hru,  unitarea_bsnhru,  spval)
      CALL worker_push_data (push_elmhru2bsnhru, lakedepth_hru, lakedepth_bsnhru, spval)

      IF (allocated (lake_id_elm    )) deallocate (lake_id_elm    )
      IF (allocated (lakedepth_hru  )) deallocate (lakedepth_hru  )
      IF (allocated (lakedown_id_elm)) deallocate (lakedown_id_elm)
      IF (allocated (lakedown_id_bsn)) deallocate (lakedown_id_bsn)
      IF (allocated (lakeoutlet_elm )) deallocate (lakeoutlet_elm )
      IF (allocated (unitarea_hru   )) deallocate (unitarea_hru   )

      IF (p_is_worker) THEN

         IF (numbasin > 0) THEN

            allocate (lakeinfo    (numbasin))
            allocate (riverarea   (numbasin))
            allocate (riverwth    (numbasin))
            allocate (bedelv      (numbasin))
            allocate (handmin     (numbasin))
            allocate (wtsrfelv    (numbasin))
            allocate (riverlen_ds (numbasin))
            allocate (wtsrfelv_ds (numbasin))
            allocate (riverwth_ds (numbasin))
            allocate (bedelv_ds   (numbasin))
            allocate (outletwth   (numbasin))

            DO ibasin = 1, numbasin

               hs = basin_hru%substt(ibasin)
               he = basin_hru%subend(ibasin)

               IF (lake_id(ibasin) == 0) THEN

                  hillslope_basin(ibasin)%area = unitarea_bsnhru(hs:he)

                  riverarea(ibasin) = hillslope_basin(ibasin)%area(1)
                  riverwth (ibasin) = riverarea(ibasin) / riverlen(ibasin)

                  ! modify height above nearest drainage data to consider river depth
                  IF (hillslope_basin(ibasin)%nhru > 1) THEN
                     hillslope_basin(ibasin)%hand(2:) = &
                        hillslope_basin(ibasin)%hand(2:) + riverdpth(ibasin)
                  ENDIF

                  wtsrfelv(ibasin) = riverelv(ibasin)
                  bedelv  (ibasin) = riverelv(ibasin) - riverdpth(ibasin)

                  handmin(ibasin) = minval(hillslope_basin(ibasin)%hand)

               ELSEIF (lake_id(ibasin) > 0) THEN

                  wtsrfelv(ibasin) = basinelv(ibasin)

                  bedelv(ibasin) = basinelv(ibasin) - maxval(lakedepth_bsnhru(hs:he))

                  nsublake = he - hs + 1
                  lakeinfo(ibasin)%nsub = nsublake

                  allocate (lakeinfo(ibasin)%area0  (nsublake))
                  allocate (lakeinfo(ibasin)%area   (nsublake))
                  allocate (lakeinfo(ibasin)%depth0 (nsublake))
                  allocate (lakeinfo(ibasin)%depth  (nsublake))

                  lakeinfo(ibasin)%area  = unitarea_bsnhru (hs:he)
                  lakeinfo(ibasin)%depth = lakedepth_bsnhru(hs:he)

                  ! area data in HRU order
                  lakeinfo(ibasin)%area0 = lakeinfo(ibasin)%area

                  ! depth data in HRU order
                  lakeinfo(ibasin)%depth0 = lakeinfo(ibasin)%depth

                  allocate (order (1:nsublake))
                  order = (/(i, i = 1, nsublake)/)

                  CALL quicksort (nsublake, lakeinfo(ibasin)%depth, order)

                  ! area data in depth order
                  lakeinfo(ibasin)%area = lakeinfo(ibasin)%area(order)

                  ! adjust to be from deepest to shallowest
                  lakeinfo(ibasin)%depth = lakeinfo(ibasin)%depth(nsublake:1:-1)
                  lakeinfo(ibasin)%area  = lakeinfo(ibasin)%area (nsublake:1:-1)

                  allocate (lakeinfo(ibasin)%dep_vol_curve (nsublake))

                  lakeinfo(ibasin)%dep_vol_curve(1) = 0
                  DO i = 2, nsublake
                     lakeinfo(ibasin)%dep_vol_curve(i) = &
                        lakeinfo(ibasin)%dep_vol_curve(i-1) &
                        + sum(lakeinfo(ibasin)%area(1:i-1)) &
                        * (lakeinfo(ibasin)%depth(i-1)-lakeinfo(ibasin)%depth(i))
                  ENDDO

                  riverlen(ibasin) = 0.

                  deallocate (order)

               ELSEIF (lake_id(ibasin) < 0) THEN

                  hillslope_basin(ibasin)%area = unitarea_bsnhru(hs:he)
                  handmin(ibasin) = minval(hillslope_basin(ibasin)%hand)

               ENDIF

            ENDDO
         ENDIF

         IF (numbasin > 0) THEN
            allocate(bsnfilter (numbasin))
            bsnfilter(:) = .true.
         ENDIF

         CALL pull_from_downstream (riverlen, riverlen_ds, bsnfilter)
         CALL pull_from_downstream (wtsrfelv, wtsrfelv_ds, bsnfilter)
         CALL pull_from_downstream (riverwth, riverwth_ds, bsnfilter)
         CALL pull_from_downstream (bedelv  , bedelv_ds  , bsnfilter)

         IF (allocated(bsnfilter)) deallocate(bsnfilter)

         DO ibasin = 1, numbasin
            IF (lake_id(ibasin) < 0) THEN
               bedelv(ibasin) = wtsrfelv_ds(ibasin) + minval(hillslope_basin(ibasin)%hand)
            ENDIF
         ENDDO

         DO ibasin = 1, numbasin
            IF (lake_id(ibasin) == 0) THEN
               IF ((to_lake(ibasin)) .or. (riverdown(ibasin) <= 0)) THEN
                  ! river to lake, ocean, inland depression or out of region
                  outletwth(ibasin) = riverwth(ibasin)
               ELSE
                  ! river to river
                  outletwth(ibasin) = (riverwth(ibasin) + riverwth_ds(ibasin)) * 0.5
               ENDIF
            ELSEIF (lake_id(ibasin) /= 0) THEN
               IF ((.not. to_lake(ibasin)) .and. (riverdown(ibasin) /= 0)) THEN
                  IF (riverdown(ibasin) > 0) THEN
                     ! lake to river
                     outletwth(ibasin) = riverwth_ds(ibasin)
                  ELSEIF (riverdown(ibasin) == -1) THEN
                     ! lake is inland depression
                     outletwth(ibasin) = 0
                  ENDIF
               ELSEIF (to_lake(ibasin)) THEN
                  ! lake to lake .or. lake catchment to lake
                  outletwth(ibasin) = lakeoutlet_bsn(ibasin)
               ELSEIF (riverdown(ibasin) == 0) THEN
                  ! lake to ocean
                  outletwth(ibasin) = lakeoutlet_bsn(ibasin)
               ENDIF
            ENDIF
         ENDDO

      ENDIF

      IF (allocated (lakedepth_bsnhru)) deallocate (lakedepth_bsnhru)
      IF (allocated (lakeoutlet_bsn  )) deallocate (lakeoutlet_bsn  )
      IF (allocated (unitarea_bsnhru )) deallocate (unitarea_bsnhru )

      IF (allocated(basin_sorted  )) deallocate(basin_sorted  )
      IF (allocated(basin_order   )) deallocate(basin_order   )

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
      IF (p_is_master) write(*,'(A)') 'Building river network information done.'
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE river_lake_network_init

   ! ----- pull data from downstream basin -----
   SUBROUTINE pull_from_downstream (datain, dataout, bsnfilter)

   USE MOD_SPMD_Task

   IMPLICIT NONE

      real(r8), intent(in)    :: datain   (:)
      real(r8), intent(inout) :: dataout  (:)
      logical,  intent(in)    :: bsnfilter(:)

      ! local variables
      integer :: i, ibasin
      real(r8), allocatable :: datalink(:)

      IF (p_is_worker) THEN

#ifdef USEMPI
         IF (riversystem /= -1) THEN
            allocate (datalink (numbsnlink));  datalink(:) = 0.
            DO i = 1, nlink_me
               datalink(linkpull(i)) = datain(linkpush(i))
            ENDDO
            CALL mpi_allreduce (MPI_IN_PLACE, datalink, numbsnlink, MPI_REAL8, MPI_SUM, p_comm_rivsys, p_err)
         ENDIF
#endif

         DO ibasin = 1, numbasin
            IF (bsnfilter(ibasin)) THEN
               IF (riverdown(ibasin) > 0) THEN
                  IF (ilocdown(ibasin) > 0) THEN
                     dataout(ibasin) = datain(ilocdown(ibasin))
#ifdef USEMPI
                  ELSEIF (ilocdown(ibasin) < 0) THEN
                     dataout(ibasin) = datalink(-ilocdown(ibasin))
#endif
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

#ifdef USEMPI
         IF (allocated(datalink)) deallocate (datalink)
#endif
      ENDIF

   END SUBROUTINE pull_from_downstream

   ! ----- push data to downstream basin -----
   SUBROUTINE push_to_downstream (datain, dataout, bsnfilter)

   USE MOD_SPMD_Task

   IMPLICIT NONE

      real(r16), intent(in)    :: datain   (:)
      real(r16), intent(inout) :: dataout  (:)
      logical,   intent(in)    :: bsnfilter(:)

      ! local variables
      integer :: i, ibasin
      real(r16), allocatable :: datalink(:)


      IF (p_is_worker) THEN

#ifdef USEMPI
         IF (numbsnlink > 0) THEN
            allocate (datalink (numbsnlink));  datalink(:) = 0.
         ENDIF
#endif

         DO ibasin = 1, numbasin
            IF (bsnfilter(ibasin)) THEN
               IF (riverdown(ibasin) > 0) THEN
                  IF (ilocdown(ibasin) > 0) THEN
                     dataout(ilocdown(ibasin)) = dataout(ilocdown(ibasin)) + datain(ibasin)
#ifdef USEMPI
                  ELSEIF (ilocdown(ibasin) < 0) THEN
                     datalink(-ilocdown(ibasin)) = datalink(-ilocdown(ibasin)) + datain(ibasin)
#endif
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

#ifdef USEMPI
         IF (riversystem /= -1) THEN

            CALL mpi_allreduce (MPI_IN_PLACE, datalink, numbsnlink, MPI_REAL16, MPI_SUM, p_comm_rivsys, p_err)

            DO i = 1, nlink_me
               dataout(linkpush(i)) = dataout(linkpush(i)) + datalink(linkpull(i))
            ENDDO

            deallocate (datalink)
         ENDIF
#endif
      ENDIF

   END SUBROUTINE push_to_downstream

   ! ----- retrieve river depth from runoff -----
   SUBROUTINE calc_riverdepth_from_runoff (patcharea)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_DataType
   USE MOD_Utils
   USE MOD_NetCDFSerial
   USE MOD_NetCDFBlock
   USE MOD_Pixel
   USE MOD_Block
   USE MOD_Mesh
   USE MOD_Grid
   USE MOD_SpatialMapping
   USE MOD_LandElm
   USE MOD_LandPatch
   USE MOD_ElmVector
   IMPLICIT NONE

   real(r8), intent(in) :: patcharea (:)

   ! Local Variables
   character(len=256) :: file_rnof, file_rivdpt
   type(grid_type)    :: grid_rnof
   type(block_data_real8_2d)  :: f_rnof
   type(spatial_mapping_type) :: mg2p_rnof

   real(r8), allocatable :: bsnrnof(:) , bsndis(:)
   integer,  allocatable :: nups_riv(:), iups_riv(:), b_up2down(:)

   integer :: i, j, ithis, ib, jb, iblkme, ps, pe
   integer :: iwork, mesg(2), isrc, ndata
   real(r8), allocatable :: rcache(:)
   real(r8) :: myarea

   real(r8), parameter :: cH_rivdpt   = 0.1
   real(r8), parameter :: pH_rivdpt   = 0.5
   real(r8), parameter :: B0_rivdpt   = 0.0
   real(r8), parameter :: Bmin_rivdpt = 1.0


      file_rnof = trim(DEF_dir_runtime) // '/runoff_clim.nc'

      CALL grid_rnof%define_from_file (file_rnof, 'lat', 'lon')

      CALL mg2p_rnof%build_arealweighted (grid_rnof, landelm)

      IF (p_is_io) THEN
         CALL allocate_block_data (grid_rnof, f_rnof)
         CALL ncio_read_block (file_rnof, 'ro', grid_rnof, f_rnof)

         DO iblkme = 1, gblock%nblkme
            ib = gblock%xblkme(iblkme)
            jb = gblock%yblkme(iblkme)
            DO j = 1, grid_rnof%ycnt(jb)
               DO i = 1, grid_rnof%xcnt(ib)
                  f_rnof%blk(ib,jb)%val(i,j) = max(f_rnof%blk(ib,jb)%val(i,j), 0.)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      IF (p_is_worker) THEN
         IF (numelm > 0) allocate (bsnrnof (numelm))
      ENDIF

      CALL mg2p_rnof%grid2pset (f_rnof, bsnrnof)

      IF (p_is_worker) THEN
         IF (numelm > 0) THEN
            bsnrnof = bsnrnof /24.0/3600.0 ! from m/day to m/s
            DO i = 1, numelm
               ps = elm_patch%substt(i)
               pe = elm_patch%subend(i)
               myarea = sum(patcharea(ps:pe))
               ! total runoff in basin, from m/s to m3/s
               bsnrnof(i) = bsnrnof(i) * myarea
            ENDDO
         ENDIF
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_worker) THEN
         mesg = (/p_iam_glb, numelm/)
         CALL mpi_send (mesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         IF (numelm > 0) THEN
            CALL mpi_send (bsnrnof, numelm, MPI_REAL8, p_address_master, mpi_tag_data, p_comm_glb, p_err)
         ENDIF
      ENDIF

      IF (p_is_master) THEN

         allocate (bsnrnof (totalnumelm))

         DO iwork = 0, p_np_worker-1
            CALL mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               allocate(rcache (ndata))

               CALL mpi_recv (rcache, ndata, MPI_REAL8, isrc, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)

               bsnrnof(elm_data_address(p_itis_worker(isrc))%val) = rcache

               deallocate (rcache)
            ENDIF
         ENDDO
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      bsnrnof(elm_data_address(0)%val) = bsnrnof
#endif


      IF (p_is_master) THEN

         allocate (nups_riv (totalnumelm))
         allocate (iups_riv (totalnumelm))
         allocate (b_up2down(totalnumelm))

         allocate (bsndis   (totalnumelm))

         nups_riv(:) = 0
         DO i = 1, totalnumelm
            j = riverdown(i)
            IF (j > 0) THEN
               nups_riv(j) = nups_riv(j) + 1
            ENDIF
         ENDDO

         iups_riv(:) = 0
         ithis = 0
         DO i = 1, totalnumelm
            IF (iups_riv(i) == nups_riv(i)) THEN

               ithis = ithis + 1
               b_up2down(ithis) = i

               j = riverdown(i)
               DO WHILE (j > 0)

                  iups_riv(j) = iups_riv(j) + 1

                  IF (iups_riv(j) == nups_riv(j)) THEN
                     ithis = ithis + 1
                     b_up2down(ithis) = j
                     j = riverdown(j)
                  ELSE
                     EXIT
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

         bsndis(:) = 0.
         DO i = 1, totalnumelm
            j = b_up2down(i)
            bsndis(j) = bsndis(j) + bsnrnof(j)
            IF (riverdown(j) > 0) THEN
               bsndis(riverdown(j)) = bsndis(riverdown(j)) + bsndis(j)
            ENDIF
         ENDDO

         allocate (riverdpth (totalnumelm))
         DO i = 1, totalnumelm
            riverdpth(i) = max(cH_rivdpt * (bsndis(i)**pH_rivdpt) + B0_rivdpt, Bmin_rivdpt)
         ENDDO

      ENDIF

      IF (allocated (bsnrnof  )) deallocate(bsnrnof  )
      IF (allocated (bsndis   )) deallocate(bsndis   )
      IF (allocated (nups_riv )) deallocate(nups_riv )
      IF (allocated (iups_riv )) deallocate(iups_riv )
      IF (allocated (b_up2down)) deallocate(b_up2down)

   END SUBROUTINE calc_riverdepth_from_runoff

   !
   FUNCTION retrieve_lake_surface_from_volume (this, volume) result(surface)

   IMPLICIT NONE

   class(lake_info_type) :: this
   real(r8), intent(in)  :: volume
   real(r8) :: surface

   ! Local Variables
   integer :: i

      IF (volume <= 0) THEN
         surface = 0
         RETURN
      ENDIF

      IF (this%nsub == 1) THEN
         surface = volume / this%area(1)
      ELSE
         i = 1
         DO WHILE (i < this%nsub)
            IF (volume >= this%dep_vol_curve(i+1)) THEN
               i = i + 1
            ELSE
               EXIT
            ENDIF
         ENDDO
         surface = this%depth(1) - this%depth(i) + &
            (volume - this%dep_vol_curve(i)) / sum(this%area(1:i))
      ENDIF

   END FUNCTION retrieve_lake_surface_from_volume

   !
   FUNCTION retrieve_lake_volume_from_surface (this, surface) result(volume)

   IMPLICIT NONE

   class(lake_info_type) :: this
   real(r8), intent(in)  :: surface
   real(r8) :: volume

   ! Local Variables
   integer :: i

      IF (surface <= 0) THEN
         volume = 0
         RETURN
      ENDIF

      IF (this%nsub == 1) THEN
         volume = this%area(1) * surface
      ELSE
         i = 1
         DO WHILE (i < this%nsub)
            IF (surface >= this%depth(1)-this%depth(i+1)) THEN
               i = i + 1
            ELSE
               EXIT
            ENDIF
         ENDDO
         volume = this%dep_vol_curve(i) &
            + (surface - (this%depth(1) - this%depth(i))) * sum(this%area(1:i))
      ENDIF

   END FUNCTION retrieve_lake_volume_from_surface

   ! ----------
   SUBROUTINE river_lake_network_final ()

   IMPLICIT NONE

   ! Local Variables
   integer :: ilake

      IF (allocated (riverlen       ))  deallocate(riverlen       )
      IF (allocated (riverelv       ))  deallocate(riverelv       )
      IF (allocated (riverwth       ))  deallocate(riverwth       )
      IF (allocated (riverarea      ))  deallocate(riverarea      )
      IF (allocated (riverdpth      ))  deallocate(riverdpth      )
      IF (allocated (basinelv       ))  deallocate(basinelv       )
      IF (allocated (bedelv         ))  deallocate(bedelv         )
      IF (allocated (handmin        ))  deallocate(handmin        )
      IF (allocated (wtsrfelv       ))  deallocate(wtsrfelv       )
      IF (allocated (riverdown      ))  deallocate(riverdown      )
      IF (allocated (ilocdown       ))  deallocate(ilocdown       )
      IF (allocated (to_lake        ))  deallocate(to_lake        )
      IF (allocated (riverlen_ds    ))  deallocate(riverlen_ds    )
      IF (allocated (wtsrfelv_ds    ))  deallocate(wtsrfelv_ds    )
      IF (allocated (riverwth_ds    ))  deallocate(riverwth_ds    )
      IF (allocated (bedelv_ds      ))  deallocate(bedelv_ds      )
      IF (allocated (outletwth      ))  deallocate(outletwth      )
      IF (allocated (irivsys        ))  deallocate(irivsys        )
#ifdef USEMPI
      IF (allocated (linkbindex     ))  deallocate(linkbindex     )
      IF (allocated (linkpush       ))  deallocate(linkpush       )
      IF (allocated (linkpull       ))  deallocate(linkpull       )
#endif
      IF (allocated (lakeinfo       ))  deallocate(lakeinfo       )
      IF (associated(hillslope_basin))  deallocate(hillslope_basin)

   END SUBROUTINE river_lake_network_final

   ! ---------
   SUBROUTINE lake_info_free_mem (this)

   IMPLICIT NONE
   type(lake_info_type) :: this

      IF (allocated(this%area0 )) deallocate (this%area0 )
      IF (allocated(this%area  )) deallocate (this%area  )
      IF (allocated(this%depth0)) deallocate (this%depth0)
      IF (allocated(this%depth )) deallocate (this%depth )
      IF (allocated(this%dep_vol_curve)) deallocate (this%dep_vol_curve)

   END SUBROUTINE lake_info_free_mem

END MODULE MOD_Catch_RiverLakeNetwork
#endif
