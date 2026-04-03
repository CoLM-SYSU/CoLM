#include <define.h>

MODULE MOD_Mesh

!------------------------------------------------------------------------------------
! !DESCRIPTION:
!
!    MESH refers to the set of largest elements in CoLM.
!
!    In CoLM, the global/regional area is divided into a hierarchical structure:
!    1. If GRIDBASED or UNSTRUCTURED is defined, it is
!       ELEMENT >>> PATCH
!    2. If CATCHMENT is defined, it is
!       ELEMENT >>> HRU >>> PATCH
!    If Plant Function Type classification is used, PATCH is further divided into PFT.
!    If Plant Community classification is used,     PATCH is further divided into PC.
!
!    To represent ELEMENT in CoLM, the land surface is first divided into pixels,
!    which are rasterized points defined by fine-resolution data.
!
!    ELEMENT in MESH is set of pixels:
!    1. If GRIDBASED,    ELEMENT is set of pixels in a longitude-latitude rectangle.
!    2. If UNSTRUCTURED, ELEMENT is set of pixels in an irregular area (usually polygon).
!    3. If CATCHMENT,    ELEMENT is set of pixels in a catchment whose area is less than
!       a predefined value.
!
!    If GRIDBASED is defined, MESH is built by using input files containing mask of
!    land area or by defining the resolution of longitude-latitude grid.
!    If CATCHMENT or UNSTRUCTURED is defined, MESH is built by using input files
!    containing index of elements.
!
!  Created by Shupeng Zhang, May 2023
!------------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Grid
   IMPLICIT NONE

   ! ---- data types ----
   type :: irregular_elm_type

      integer*8 :: indx
      integer   :: xblk, yblk

      integer :: npxl
      integer, allocatable :: ilon(:)
      integer, allocatable :: ilat(:)

   END type irregular_elm_type

   ! ---- Instance ----
   type (grid_type) :: gridmesh

   integer :: numelm
   type (irregular_elm_type), allocatable :: mesh (:)

   integer, allocatable :: nelm_blk(:,:)

#ifdef GRIDBASED
   logical :: read_mesh_from_file = .true.
#endif

CONTAINS

   ! ------
#ifdef GRIDBASED
   SUBROUTINE init_gridbased_mesh_grid ()

   USE MOD_SPMD_Task
   USE MOD_Namelist
   IMPLICIT NONE

      IF (p_is_master) THEN
         inquire (file=trim(DEF_file_mesh), exist=read_mesh_from_file)
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast (read_mesh_from_file, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif
      IF (read_mesh_from_file) THEN
         CALL gridmesh%define_from_file (DEF_file_mesh)
      ELSE
         CALL gridmesh%define_by_res (DEF_GRIDBASED_lon_res, DEF_GRIDBASED_lat_res)
      ENDIF

   END SUBROUTINE init_gridbased_mesh_grid
#endif

   ! -------
   SUBROUTINE copy_elm (elm_from, elm_to)

   IMPLICIT NONE
   type (irregular_elm_type), intent(in)  :: elm_from
   type (irregular_elm_type), intent(out) :: elm_to

      elm_to%indx = elm_from%indx
      elm_to%npxl = elm_from%npxl
      elm_to%xblk = elm_from%xblk
      elm_to%yblk = elm_from%yblk

      IF (allocated(elm_to%ilat)) deallocate(elm_to%ilat)
      IF (allocated(elm_to%ilon)) deallocate(elm_to%ilon)

      allocate (elm_to%ilat (elm_to%npxl))
      allocate (elm_to%ilon (elm_to%npxl))
      elm_to%ilon = elm_from%ilon
      elm_to%ilat = elm_from%ilat

   END SUBROUTINE copy_elm

   ! --------------------------------
   SUBROUTINE mesh_build ()

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_NetCDFBlock
   USE MOD_Block
   USE MOD_Pixel
   USE MOD_Grid
   USE MOD_Utils
   USE MOD_DataType
   USE MOD_CatchmentDataReadin

   IMPLICIT NONE

   ! Local Variables
   type(block_data_int32_2d) :: datamesh

   integer  :: nelm, ie, je, dsp
   integer  :: iblkme, iblk, jblk, xloc, yloc, xg, yg, ixloc, iyloc
   integer  :: xp, yp, xblk, yblk, npxl, ipxl, ix, iy
   integer  :: iworker, iproc, idest, isrc, iloc, iloc_max(2)
   integer  :: ylg, yug, ysp, ynp, nyp
   integer  :: xlg, xug, xwp, xep, nxp
   real(r8) :: dlatp, dlonp
   logical  :: is_new
   integer  :: smesg(6), rmesg(6), blktag, nsend, nrecv, irecv
   integer  :: iblk_p, jblk_p, nelm_max_blk, nelm_glb

   integer, allocatable :: nelm_worker(:)
   type(pointer_int64_1d), allocatable :: elist_worker(:)

   integer*8 :: elmid
   integer*8, allocatable :: elist(:), elist2(:,:), sbuf64(:), elist_recv(:)

   integer, allocatable :: iaddr(:), elmindx  (:), order(:)
   integer, allocatable :: xlist(:), ylist    (:), npxl_(:), xlist_recv(:), ylist_recv(:), sbuf(:)
   logical, allocatable :: msk  (:), work_done(:)

   integer, allocatable :: xlist2(:,:), ylist2(:,:), ipt2(:,:), npxl_blk(:,:), blkdsp(:,:), blkcnt(:,:)
   logical, allocatable :: msk2  (:,:)

   type(irregular_elm_type), allocatable :: meshtmp (:)


      IF (p_is_io) THEN
         CALL allocate_block_data (gridmesh, datamesh)
      ENDIF

#ifdef GRIDBASED
      IF (read_mesh_from_file) THEN
         CALL ncio_read_block (DEF_file_mesh, 'landmask', gridmesh, datamesh)
      ELSE
         CALL flush_block_data (datamesh, 1)
      ENDIF
#endif

#ifdef CATCHMENT
      CALL catchment_data_read (DEF_CatchmentMesh_data, 'icatchment2d', gridmesh, datamesh, spv = -1)
#endif

#ifdef UNSTRUCTURED
      CALL ncio_read_block (DEF_file_mesh, 'elmindex', gridmesh, datamesh)
#endif

      ! Step 1: How many elms in each block?
      IF (p_is_io) THEN

         nelm = 0

         allocate (nelm_worker (0:p_np_worker-1))
         nelm_worker(:) = 0

         allocate (elist_worker (0:p_np_worker-1))
         DO iworker = 0, p_np_worker-1
            allocate (elist_worker(iworker)%val (1000))
         ENDDO

         DO iblkme = 1, gblock%nblkme
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)

            DO yloc = 1, gridmesh%ycnt(jblk)
               DO xloc = 1, gridmesh%xcnt(iblk)

#ifdef GRIDBASED
                  IF (datamesh%blk(iblk,jblk)%val(xloc,yloc) > 0) THEN
                     xg = gridmesh%xdsp(iblk) + xloc
                     IF (xg > gridmesh%nlon) xg = xg - gridmesh%nlon

                     yg = gridmesh%ydsp(jblk) + yloc

                     elmid = int(gridmesh%nlon,8) * (yg-1) + xg
                  ELSE
                     elmid = 0
                  ENDIF
#endif
#ifdef CATCHMENT
                  elmid = datamesh%blk(iblk,jblk)%val(xloc,yloc)
#endif
#ifdef UNSTRUCTURED
                  elmid = datamesh%blk(iblk,jblk)%val(xloc,yloc)
#endif

                  IF (elmid > 0) THEN

                     iworker = mod(elmid, p_np_worker)
                     CALL insert_into_sorted_list1 ( &
                        elmid, nelm_worker(iworker), elist_worker(iworker)%val, iloc)

                     IF (nelm_worker(iworker) == size(elist_worker(iworker)%val)) THEN
                        CALL expand_list (elist_worker(iworker)%val, 0.2_r8)
                     ENDIF

                  ENDIF

               ENDDO
            ENDDO

#ifdef USEMPI
            DO iworker = 0, p_np_worker-1
               IF (nelm_worker(iworker) > 0) THEN
                  idest = p_address_worker(iworker)
                  smesg(1:2) = (/p_iam_glb, nelm_worker(iworker)/)
                  ! send(01)
                  CALL mpi_send (smesg(1:2), 2, MPI_INTEGER, &
                     idest, mpi_tag_size, p_comm_glb, p_err)
               ENDIF
            ENDDO
#endif

            nelm = nelm + sum(nelm_worker)
            nelm_worker(:) = 0
         ENDDO

#ifdef USEMPI
         DO iworker = 0, p_np_worker-1
            idest = p_address_worker(iworker)
            ! send(02)
            smesg(1:2) = (/p_iam_glb, 0/)
            CALL mpi_send (smesg(1:2), 2, MPI_INTEGER, &
               idest, mpi_tag_size, p_comm_glb, p_err)
         ENDDO
#endif

         deallocate (nelm_worker)
         DO iworker = 0, p_np_worker-1
            deallocate (elist_worker(iworker)%val)
         ENDDO
         deallocate (elist_worker)

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN
         nelm = 0
         allocate(work_done(0:p_np_io-1))
         work_done(:) = .false.
         DO WHILE (.not. all(work_done))
            ! recv(01,02)
            CALL mpi_recv (rmesg(1:2), 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_size, p_comm_glb, p_stat, p_err)

            isrc  = rmesg(1)
            nrecv = rmesg(2)

            IF (nrecv > 0) THEN
               nelm = nelm + nrecv
            ELSE
               work_done(p_itis_io(isrc)) = .true.
            ENDIF
         ENDDO

         deallocate(work_done)
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      ! Step 2: Build pixel list for each elm.
      IF (p_is_worker) THEN
         IF (nelm > 0) THEN
            allocate (meshtmp (nelm))
            allocate (elist (nelm))
            allocate (iaddr (nelm))
         ENDIF
         nelm = 0
      ENDIF

      IF (p_is_io) THEN

         DO iblkme = 1, gblock%nblkme
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)
            IF (gridmesh%xcnt(iblk) <= 0) CYCLE
            IF (gridmesh%ycnt(jblk) <= 0) CYCLE

            ylg = gridmesh%ydsp(jblk) + 1
            yug = gridmesh%ydsp(jblk) + gridmesh%ycnt(jblk)
            IF (gridmesh%yinc == 1) THEN
               ysp = find_nearest_south (gridmesh%lat_s(ylg), pixel%nlat, pixel%lat_s)
               ynp = find_nearest_north (gridmesh%lat_n(yug), pixel%nlat, pixel%lat_n)
            ELSE
               ysp = find_nearest_south (gridmesh%lat_s(yug), pixel%nlat, pixel%lat_s)
               ynp = find_nearest_north (gridmesh%lat_n(ylg), pixel%nlat, pixel%lat_n)
            ENDIF

            nyp = ynp - ysp + 1

            xlg = gridmesh%xdsp(iblk) + 1
            xug = gridmesh%xdsp(iblk) + gridmesh%xcnt(iblk)
            IF (xug > gridmesh%nlon) xug = xug - gridmesh%nlon

            xwp = find_nearest_west (gridmesh%lon_w(xlg), pixel%nlon, pixel%lon_w)
            IF (.not. lon_between_floor(pixel%lon_w(xwp), gridmesh%lon_w(xlg), gridmesh%lon_e(xlg))) THEN
               xwp = mod(xwp,pixel%nlon) + 1
            ENDIF

            xep = find_nearest_east (gridmesh%lon_e(xug), pixel%nlon, pixel%lon_e)
            IF (.not. lon_between_ceil(pixel%lon_e(xep), gridmesh%lon_w(xug), gridmesh%lon_e(xug))) THEN
               xep = xep - 1
               IF (xep == 0) xep = pixel%nlon
            ENDIF

            nxp = xep - xwp + 1
            IF (nxp <= 0) nxp = nxp + pixel%nlon

            allocate (elist2 (nxp,nyp))
            allocate (xlist2 (nxp,nyp))
            allocate (ylist2 (nxp,nyp))
            allocate (msk2   (nxp,nyp))

            DO iy = ysp, ynp
               yg = gridmesh%ygrd(iy)
               yloc = gridmesh%yloc(yg)

               iyloc = iy - ysp + 1
               dlatp = pixel%lat_n(iy) - pixel%lat_s(iy)
               IF (dlatp < 1.0e-6_r8) THEN
                  elist2(:,iyloc) = 0
                  CYCLE
               ENDIF

               ix = xwp
               ixloc = 0
               DO WHILE (.true.)
                  ixloc = ixloc + 1
                  dlonp = pixel%lon_e(ix) - pixel%lon_w(ix)
                  IF (dlonp < 0) dlonp = dlonp + 360.0_r8

                  xg = gridmesh%xgrd(ix)
                  xloc = gridmesh%xloc(xg)

#ifdef GRIDBASED
                  IF (datamesh%blk(iblk,jblk)%val(xloc,yloc) > 0) THEN
                     elmid = int(gridmesh%nlon,8) * (yg-1) + xg
                  ELSE
                     elmid = 0
                  ENDIF
#endif
#ifdef CATCHMENT
                  elmid = datamesh%blk(iblk,jblk)%val(xloc,yloc)
#endif
#ifdef UNSTRUCTURED
                  elmid = datamesh%blk(iblk,jblk)%val(xloc,yloc)
#endif

                  xlist2(ixloc,iyloc) = ix
                  ylist2(ixloc,iyloc) = iy
                  elist2(ixloc,iyloc) = elmid

                  IF (dlonp < 1.0e-6_r8) THEN
                     elist2(ixloc,iyloc) = 0
                  ENDIF

                  IF (ix == xep) EXIT
                  ix = mod(ix,pixel%nlon) + 1
               ENDDO
            ENDDO

#ifdef USEMPI
            allocate (sbuf (nxp*nyp))
            allocate (ipt2 (nxp,nyp))

            allocate (sbuf64 (nxp*nyp))

            blktag = iblkme
            ipt2 = mod(elist2, p_np_worker)
            DO iproc = 0, p_np_worker-1
               msk2  = (ipt2 == iproc) .and. (elist2 > 0)
               nsend = count(msk2)
               IF (nsend > 0) THEN

                  idest = p_address_worker(iproc)

                  smesg(1:3) = (/p_iam_glb, nsend, blktag/)
                  ! send(03)
                  CALL mpi_send (smesg(1:3), 3, MPI_INTEGER, &
                     idest, mpi_tag_mesg, p_comm_glb, p_err)

                  sbuf64(1:nsend) = pack(elist2, msk2)
                  ! send(04)
                  CALL mpi_send (sbuf64(1:nsend), nsend, MPI_INTEGER8, &
                     idest, blktag, p_comm_glb, p_err)

                  sbuf(1:nsend) = pack(xlist2, msk2)
                  ! send(05)
                  CALL mpi_send (sbuf(1:nsend), nsend, MPI_INTEGER, &
                     idest, blktag, p_comm_glb, p_err)

                  sbuf(1:nsend) = pack(ylist2, msk2)
                  ! send(06)
                  CALL mpi_send (sbuf(1:nsend), nsend, MPI_INTEGER, &
                     idest, blktag, p_comm_glb, p_err)

               ENDIF
            ENDDO

            deallocate (sbuf  )
            deallocate (ipt2  )
            deallocate (sbuf64)
#else

            DO iy = 1, nyp
               DO ix = 1, nxp

                  elmid = elist2(ix,iy)
                  IF (elmid > 0) THEN

                     CALL insert_into_sorted_list1 (elmid, nelm, elist, iloc, is_new)

                     msk2 = (elist2 == elmid)
                     npxl = count(msk2)

                     IF (is_new) THEN
                        IF (iloc < nelm) THEN
                           iaddr(iloc+1:nelm) = iaddr(iloc:nelm-1)
                        ENDIF
                        iaddr(iloc) = nelm

                        meshtmp(iaddr(iloc))%indx = elmid
                        meshtmp(iaddr(iloc))%npxl = npxl
                     ELSE
                        meshtmp(iaddr(iloc))%npxl = meshtmp(iaddr(iloc))%npxl + npxl
                     ENDIF

                     allocate (xlist(npxl))
                     allocate (ylist(npxl))
                     xlist = pack(xlist2, msk2)
                     ylist = pack(ylist2, msk2)

                     CALL append_to_list (meshtmp(iaddr(iloc))%ilon, xlist)
                     CALL append_to_list (meshtmp(iaddr(iloc))%ilat, ylist)

                     WHERE(msk2) elist2 = -1

                     deallocate (xlist)
                     deallocate (ylist)
                  ENDIF

               ENDDO
            ENDDO
#endif

            deallocate (elist2)
            deallocate (xlist2)
            deallocate (ylist2)
            deallocate (msk2  )

         ENDDO

#ifdef USEMPI
         DO iworker = 0, p_np_worker-1
            idest = p_address_worker(iworker)
            ! send(07)
            smesg(1:3) = (/p_iam_glb, 0, 0/)
            CALL mpi_send (smesg(1:3), 3, MPI_INTEGER, &
               idest, mpi_tag_mesg, p_comm_glb, p_err)
         ENDDO
#endif

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN

         allocate(work_done(0:p_np_io-1))
         work_done(:) = .false.
         DO WHILE (.not. all(work_done))
            ! recv(03,07)
            CALL mpi_recv (rmesg(1:3), 3, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc   = rmesg(1)
            nrecv  = rmesg(2)
            blktag = rmesg(3)
            IF (nrecv > 0) THEN

               allocate (elist_recv (nrecv))
               ! recv(04)
               CALL mpi_recv (elist_recv, nrecv, MPI_INTEGER8, &
                  isrc, blktag, p_comm_glb, p_stat, p_err)

               allocate (xlist_recv (nrecv))
               ! recv(05)
               CALL mpi_recv (xlist_recv, nrecv, MPI_INTEGER, &
                  isrc, blktag, p_comm_glb, p_stat, p_err)

               allocate (ylist_recv (nrecv))
               ! recv(06)
               CALL mpi_recv (ylist_recv, nrecv, MPI_INTEGER, &
                  isrc, blktag, p_comm_glb, p_stat, p_err)

               allocate (msk(nrecv))

               DO irecv = 1, nrecv

                  elmid = elist_recv(irecv)

                  IF (elmid > 0) THEN

                     CALL insert_into_sorted_list1 (elmid, nelm, elist, iloc, is_new)

                     msk  = (elist_recv == elmid)
                     npxl = count(msk)

                     IF (is_new) THEN
                        IF (iloc < nelm) THEN
                           iaddr(iloc+1:nelm) = iaddr(iloc:nelm-1)
                        ENDIF
                        iaddr(iloc) = nelm

                        meshtmp(iaddr(iloc))%indx = elmid
                        meshtmp(iaddr(iloc))%npxl = npxl
                     ELSE
                        meshtmp(iaddr(iloc))%npxl = meshtmp(iaddr(iloc))%npxl + npxl
                     ENDIF

                     allocate (xlist(npxl))
                     allocate (ylist(npxl))
                     xlist = pack(xlist_recv, msk)
                     ylist = pack(ylist_recv, msk)

                     CALL append_to_list (meshtmp(iaddr(iloc))%ilon, xlist)
                     CALL append_to_list (meshtmp(iaddr(iloc))%ilat, ylist)

                     WHERE(msk) elist_recv = -1
                     deallocate (xlist)
                     deallocate (ylist)
                  ENDIF

               ENDDO

               deallocate (msk)
               deallocate (elist_recv)
               deallocate (xlist_recv)
               deallocate (ylist_recv)
            ELSE
               work_done(p_itis_io(isrc)) = .true.
            ENDIF
         ENDDO

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (allocated(elist)) deallocate (elist)
      IF (allocated(iaddr)) deallocate (iaddr)

      ! Step 3: Which block each elm locates at.
      IF (p_is_worker) THEN

         allocate (npxl_blk (gblock%nxblk,gblock%nyblk))
         allocate (nelm_blk (gblock%nxblk,gblock%nyblk))

         nelm_blk(:,:) = 0

         DO ie = 1, nelm

            npxl_blk (:,:) = 0

            DO ipxl = 1, meshtmp(ie)%npxl
               xp = meshtmp(ie)%ilon(ipxl)
               yp = meshtmp(ie)%ilat(ipxl)

               xg = gridmesh%xgrd(xp)
               yg = gridmesh%ygrd(yp)

               xblk = gridmesh%xblk(xg)
               yblk = gridmesh%yblk(yg)

               npxl_blk(xblk,yblk) = npxl_blk(xblk,yblk) + 1
            ENDDO

            iloc_max = maxloc(npxl_blk)
            meshtmp(ie)%xblk = iloc_max(1)
            meshtmp(ie)%yblk = iloc_max(2)

            nelm_blk(iloc_max(1), iloc_max(2)) = &
               nelm_blk(iloc_max(1), iloc_max(2)) + 1

         ENDDO

         deallocate (npxl_blk)
      ENDIF

#ifdef USEMPI
      IF (.not. p_is_worker) THEN
         allocate (nelm_blk (gblock%nxblk,gblock%nyblk))
         nelm_blk(:,:) = 0
      ENDIF

      CALL mpi_allreduce (MPI_IN_PLACE, nelm_blk, gblock%nxblk*gblock%nyblk, &
         MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif

      ! Step 4: IF MPI is used, sending elms from worker to their IO processes.

      IF (p_is_io) THEN

         allocate (blkdsp (gblock%nxblk, gblock%nyblk))
         blkdsp(1,1) = 0
         DO iblk = 1, gblock%nxblk
            DO jblk = 1, gblock%nyblk
               IF ((iblk /= 1) .or. (jblk /= 1)) THEN
                  IF (jblk == 1) THEN
                     iblk_p = iblk - 1
                     jblk_p = gblock%nyblk
                  ELSE
                     iblk_p = iblk
                     jblk_p = jblk - 1
                  ENDIF

                  IF (gblock%pio(iblk_p,jblk_p) == p_iam_glb) THEN
                     blkdsp(iblk,jblk) = blkdsp(iblk_p,jblk_p) + nelm_blk(iblk_p,jblk_p)
                  ELSE
                     blkdsp(iblk,jblk) = blkdsp(iblk_p,jblk_p)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN
         DO iblk = 1, gblock%nxblk
            DO jblk = 1, gblock%nyblk

               idest = gblock%pio(iblk,jblk)

               nsend = 0
               npxl  = 0
               DO ie = 1, nelm
                  IF ((meshtmp(ie)%xblk == iblk) .and. (meshtmp(ie)%yblk == jblk)) THEN
                     nsend = nsend + 1
                     npxl  = npxl  + meshtmp(ie)%npxl
                  ENDIF
               ENDDO

               IF (nsend > 0) THEN

                  allocate (elist (nsend))
                  allocate (npxl_ (nsend))
                  allocate (xlist (npxl ))
                  allocate (ylist (npxl ))

                  nsend = 0
                  npxl  = 0
                  DO ie = 1, nelm
                     IF ((meshtmp(ie)%xblk == iblk) .and. (meshtmp(ie)%yblk == jblk)) THEN

                        nsend = nsend + 1

                        elist(nsend) = meshtmp(ie)%indx
                        npxl_(nsend) = meshtmp(ie)%npxl

                        xlist(npxl+1:npxl+meshtmp(ie)%npxl) = meshtmp(ie)%ilon
                        ylist(npxl+1:npxl+meshtmp(ie)%npxl) = meshtmp(ie)%ilat

                        npxl = npxl + meshtmp(ie)%npxl
                     ENDIF
                  ENDDO

                  blktag = p_iam_glb + 1000

                  ! send(09)
                  smesg(1:6) = (/p_iam_glb, blktag, iblk, jblk, nsend, npxl/)
                  CALL mpi_send (smesg(1:6), 6, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)

                  ! send(10)
                  CALL mpi_send (elist, nsend, MPI_INTEGER8, idest, blktag, p_comm_glb, p_err)
                  CALL mpi_send (npxl_, nsend, MPI_INTEGER,  idest, blktag, p_comm_glb, p_err)

                  ! send(11)
                  CALL mpi_send (xlist, npxl, MPI_INTEGER, idest, blktag, p_comm_glb, p_err)
                  CALL mpi_send (ylist, npxl, MPI_INTEGER, idest, blktag, p_comm_glb, p_err)

                  deallocate (elist)
                  deallocate (npxl_)
                  deallocate (xlist)
                  deallocate (ylist)

               ENDIF
            ENDDO
         ENDDO
      ENDIF

      IF (p_is_io) THEN

         numelm = sum(nelm_blk, mask = gblock%pio == p_iam_glb)

         IF (numelm > 0) THEN

            allocate (mesh (numelm))

            allocate (blkcnt (gblock%nxblk, gblock%nyblk))
            blkcnt(:,:) = 0

            DO WHILE (sum(blkcnt) < numelm)

               ! recv(09)
               CALL mpi_recv (rmesg(1:6), 6, MPI_INTEGER, MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
               isrc   = rmesg(1)
               blktag = rmesg(2)
               xblk   = rmesg(3)
               yblk   = rmesg(4)
               nrecv  = rmesg(5)
               npxl   = rmesg(6)

               allocate (elist (nrecv))
               allocate (npxl_ (nrecv))
               allocate (xlist (npxl ))
               allocate (ylist (npxl ))

               ! recv(10)
               CALL mpi_recv (elist, nrecv, MPI_INTEGER8, isrc, blktag, p_comm_glb, p_stat, p_err)
               CALL mpi_recv (npxl_, nrecv, MPI_INTEGER,  isrc, blktag, p_comm_glb, p_stat, p_err)

               ! recv(11)
               CALL mpi_recv (xlist, npxl, MPI_INTEGER, isrc, blktag, p_comm_glb, p_stat, p_err)
               CALL mpi_recv (ylist, npxl, MPI_INTEGER, isrc, blktag, p_comm_glb, p_stat, p_err)

               dsp = 0
               DO ie = 1, nrecv

                  je = blkdsp(xblk,yblk) + blkcnt(xblk,yblk) + ie

                  mesh(je)%indx = elist(ie)
                  mesh(je)%xblk = xblk
                  mesh(je)%yblk = yblk
                  mesh(je)%npxl = npxl_(ie)

                  allocate (mesh(je)%ilon (npxl_(ie)))
                  allocate (mesh(je)%ilat (npxl_(ie)))

                  mesh(je)%ilon = xlist(dsp+1:dsp+npxl_(ie))
                  mesh(je)%ilat = ylist(dsp+1:dsp+npxl_(ie))

                  dsp = dsp + npxl_(ie)

               ENDDO

               blkcnt(xblk,yblk) = blkcnt(xblk,yblk) + nrecv

               deallocate (elist)
               deallocate (npxl_)
               deallocate (xlist)
               deallocate (ylist)

            ENDDO

         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)

#else
      numelm = nelm
      IF (numelm > 0) THEN

         allocate (mesh (numelm))

         allocate (blkcnt (gblock%nxblk, gblock%nyblk))
         blkcnt(:,:) = 0
         DO ie = 1, numelm

            xblk = meshtmp(ie)%xblk
            yblk = meshtmp(ie)%yblk

            blkcnt(xblk,yblk) = blkcnt(xblk,yblk) + 1
            je = blkdsp(xblk,yblk) + blkcnt(xblk,yblk)

            CALL copy_elm (meshtmp(ie), mesh(je))

         ENDDO

      ENDIF
#endif

      ! Step 4-2: sort elms.
      IF (p_is_io) THEN
         IF (allocated (meshtmp)) THEN
            DO ie = 1, size(meshtmp)
               IF (allocated(meshtmp(ie)%ilon))  deallocate (meshtmp(ie)%ilon)
               IF (allocated(meshtmp(ie)%ilat))  deallocate (meshtmp(ie)%ilat)
            ENDDO
            deallocate (meshtmp)
         ENDIF

         IF (numelm > 0) THEN
            allocate (meshtmp (numelm))
            DO ie = 1, numelm
               CALL copy_elm(mesh(ie), meshtmp(ie))
            ENDDO

            DO iblkme = 1, gblock%nblkme
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)

               IF (blkcnt(iblk,jblk) > 0) THEN
                  allocate (elmindx (blkcnt(iblk,jblk)))
                  allocate (order   (blkcnt(iblk,jblk)))

                  DO ie = blkdsp(iblk,jblk)+1, blkdsp(iblk,jblk)+blkcnt(iblk,jblk)
                     elmindx(ie-blkdsp(iblk,jblk)) = mesh(ie)%indx
                  ENDDO

                  order = (/ (ie, ie = 1, blkcnt(iblk,jblk)) /)
                  CALL quicksort (blkcnt(iblk,jblk), elmindx, order)

                  DO ie = 1, blkcnt(iblk,jblk)
                     CALL copy_elm (meshtmp(blkdsp(iblk,jblk)+order(ie)), &
                        mesh(blkdsp(iblk,jblk)+ie))
                  ENDDO

                  deallocate (elmindx)
                  deallocate (order  )
               ENDIF

            ENDDO
         ENDIF
      ENDIF

      IF (allocated(blkdsp)) deallocate(blkdsp)
      IF (allocated(blkcnt)) deallocate(blkcnt)

      IF (allocated (meshtmp)) THEN
         DO ie = 1, size(meshtmp)
            IF (allocated(meshtmp(ie)%ilon))  deallocate (meshtmp(ie)%ilon)
            IF (allocated(meshtmp(ie)%ilat))  deallocate (meshtmp(ie)%ilat)
         ENDDO

         deallocate (meshtmp )
      ENDIF

      ! Step 5: IF MPI is used, scatter elms from IO to workers.
#ifdef USEMPI
      CALL scatter_mesh_from_io_to_worker ()
#endif

      IF (p_is_master) THEN
         write(*,'(A)') 'Making mesh elements:'
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_io) THEN

         CALL mpi_reduce (numelm, nelm_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_io, p_err)
         IF (p_iam_io == p_root) THEN
            write(*,'(A,I12,A)') 'Total   : ', nelm_glb, ' elements.'
         ENDIF

         nelm_max_blk = maxval(nelm_blk, mask = gblock%pio == p_iam_glb)
         CALL mpi_allreduce (MPI_IN_PLACE, nelm_max_blk, 1, MPI_INTEGER, MPI_MAX, p_comm_io, p_err)
         IF (p_iam_io == p_root) THEN
            write(*,'(A,I12,A)') 'Maximum : ', nelm_max_blk, &
               ' elements in one block (More than 3600 is recommended).'
            write(*,'(/,A)') '   -----------------------------------------------------------------'
            write(*,'(A)')   '   |  Examples for setting of blocks and processor groupsize:      |'
            write(*,'(A)')   '   |  Resolution  DEF_nx_blocks  DEF_ny_blocks  DEF_PIO_groupsize  |'
            write(*,'(A)')   '   |         2x2             18              9                 15  |'
            write(*,'(A)')   '   |         1x1             18              9                 24  |'
            write(*,'(A)')   '   |     0.5x0.5             18              9                 36  |'
            write(*,'(A)')   '   |   0.25x0.25             30             15                 45  |'
            write(*,'(A)')   '   |     0.1x0.1             72             36                 64  |'
            write(*,'(A,/)') '   -----------------------------------------------------------------'
         ENDIF

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      write(*,'(A,I12,A)') 'Total: ', numelm, ' elements.'
#endif

   END SUBROUTINE mesh_build


#ifdef USEMPI
   ! --------------------------------
   SUBROUTINE scatter_mesh_from_io_to_worker

   USE MOD_SPMD_Task
   USE MOD_Block
   IMPLICIT NONE

   ! Local variables
   integer :: iblk, jblk, nave, nres, iproc, ndsp, nsend, idest, ie
   integer :: smesg(4), rmesg(4)
   integer, allocatable :: nelm_worker(:)
   integer :: iblkme
   character(len=20) :: wfmt

      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_io) THEN

         allocate (nelm_worker (1:p_np_group-1))
         nelm_worker(:) = 0

         DO iblkme = 1, gblock%nblkme
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)

            nave = nelm_blk(iblk,jblk) / (p_np_group-1)
            nres = mod(nelm_blk(iblk,jblk), p_np_group-1)
            DO iproc = 1, p_np_group-1
               nelm_worker(iproc) = nelm_worker(iproc) + nave
               IF (iproc <= nres)  nelm_worker(iproc) = nelm_worker(iproc) + 1
            ENDDO
         ENDDO

         IF (any(nelm_worker == 0)) THEN
            write(*,'(A,/,A)') 'Warning: there are idle workers, please use less processors ' // &
               'OR larger working group ', '  (set by DEF_PIO_groupsize in CoLM namelist).'
            write(wfmt,'(A,I0,A)') '(A,I6,A,', p_np_group-1, '(X,I0))'
            write(*,wfmt) 'Numbers of elements by workers in group ', p_iam_glb, ' are ', nelm_worker
         ENDIF

         DO iproc = 1, p_np_group-1
            CALL mpi_send (nelm_worker(iproc), 1, MPI_INTEGER, &
               iproc, mpi_tag_size, p_comm_group, p_err)
         ENDDO
         deallocate (nelm_worker)

         ndsp = 0
         DO iblkme = 1, gblock%nblkme
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)

            nave = nelm_blk(iblk,jblk) / (p_np_group-1)
            nres = mod(nelm_blk(iblk,jblk), p_np_group-1)
            DO iproc = 1, p_np_group-1
               nsend = nave
               IF (iproc <= nres)  nsend = nsend + 1

               DO ie = ndsp+1, ndsp+nsend
                  idest = iproc
                  CALL mpi_send (mesh(ie)%indx, 1, MPI_INTEGER8, &
                     idest, mpi_tag_mesg, p_comm_group, p_err)
                  smesg(1:3) = (/mesh(ie)%xblk, mesh(ie)%yblk, mesh(ie)%npxl/)
                  CALL mpi_send (smesg(1:3), 3, MPI_INTEGER, &
                     idest, mpi_tag_mesg, p_comm_group, p_err)
                  CALL mpi_send (mesh(ie)%ilon, mesh(ie)%npxl, &
                     MPI_INTEGER, idest, mpi_tag_data, p_comm_group, p_err)
                  CALL mpi_send (mesh(ie)%ilat, mesh(ie)%npxl, &
                     MPI_INTEGER, idest, mpi_tag_data, p_comm_group, p_err)
               ENDDO
               ndsp = ndsp + nsend
            ENDDO
         ENDDO

      ENDIF

      IF (p_is_worker) THEN

         CALL mpi_recv (numelm, 1, MPI_INTEGER, &
            p_root, mpi_tag_size, p_comm_group, p_stat, p_err)

         IF (numelm > 0) THEN
            allocate (mesh (numelm))

            DO ie = 1, numelm
               CALL mpi_recv (mesh(ie)%indx, 1, MPI_INTEGER8, &
                  p_root, mpi_tag_mesg, p_comm_group, p_stat, p_err)
               CALL mpi_recv (rmesg, 3, MPI_INTEGER, &
                  p_root, mpi_tag_mesg, p_comm_group, p_stat, p_err)

               mesh(ie)%xblk = rmesg(1)
               mesh(ie)%yblk = rmesg(2)
               mesh(ie)%npxl = rmesg(3)

               allocate (mesh(ie)%ilon (mesh(ie)%npxl))
               allocate (mesh(ie)%ilat (mesh(ie)%npxl))

               CALL mpi_recv (mesh(ie)%ilon, mesh(ie)%npxl, MPI_INTEGER, &
                  p_root, mpi_tag_data, p_comm_group, p_stat, p_err)
               CALL mpi_recv (mesh(ie)%ilat, mesh(ie)%npxl, MPI_INTEGER, &
                  p_root, mpi_tag_data, p_comm_group, p_stat, p_err)
            ENDDO

         ENDIF

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)

   END SUBROUTINE scatter_mesh_from_io_to_worker

#endif

   ! --------------------------------
   SUBROUTINE mesh_free_mem ()

   IMPLICIT NONE

   ! Local variables
   integer :: ie

      IF (allocated(mesh)) THEN
         DO ie = 1, numelm
            deallocate (mesh(ie)%ilon)
            deallocate (mesh(ie)%ilat)
         ENDDO

         deallocate (mesh)
      ENDIF

   END SUBROUTINE mesh_free_mem

END MODULE MOD_Mesh
