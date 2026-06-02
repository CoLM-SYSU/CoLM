#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Obs_PM
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Passive-microwave observation handling for data assimilation.
!
!    This module defines one passive-microwave observation stream (PM_type)
!    and a collection of streams (PM_set_type). A PM_type instance represents
!    one source/sensor/variable configuration and owns the observation reader,
!    part-space forward calculation, and mapping to observation-grid cells.
!
!    The grid-level workflow is:
!      1. clear time-step data and read current observations through
!         grid_obs_type;
!      2. map observation incidence angle from grid cells to patch-grid
!         intersection parts;
!      3. compute ensemble H(x) on valid parts with the passive-microwave
!         radiative-transfer operator, then area-weight valid parts back to
!         patch-space H(x);
!      4. map patch-space H(x) to the observation grid;
!      5. crop mapped H(x) to valid observation locations;
!      6. gather IO-task results on the master and broadcast this%hx so all
!         MPI tasks share the same observation-space H(x).
!
!    PM_set_type applies this workflow to all configured PM streams and
!    concatenates their observation values, errors, metadata, and H(x) into
!    set-level arrays for assimilation. Patch-space H(x), area-weighted from
!    valid parts, remains stored in each PM_type instance as hx_pset for
!    diagnostics.
!
! HISTORY:
!    Lu Li, 05/2026: First implementation
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_DataType
   USE MOD_DA_Obs_BasicType
   USE MOD_DA_Operator_RTM
   USE MOD_DA_Vars_TimeVariables
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_1DForcing
   USE MOD_Vars_Global
   USE MOD_LandPatch
   USE MOD_Block
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_SPMD_Task
   USE MOD_Pixelset
   USE MOD_TimeManager
   IMPLICIT NONE


   ! process one PM data
   type :: PM_type

      ! metadata of specific passive microwave observation stream
      real(r8)              :: fghz

      ! selected avaliable observation values and metadata at current timestep
      type(grid_obs_type)   :: y

      ! H(x) from forward model for diagnostics and data assimilation
      real(r8), allocatable :: hx_pset(:,:)  ! H(x) on patches for diagnostics [0:ens, patch]
      real(r8), allocatable :: hx(:,:)       ! H(x) mapped to y locations for assimilation 

   CONTAINS

      procedure, PUBLIC :: init     => PM_init         ! initialize the obs type (grid info, mapping, etc.)
      procedure, PUBLIC :: calcg    => PM_calc_on_grid ! calculate H(x) and y at y locations 
      procedure, PUBLIC :: calcp    => PM_calc_on_pset ! calculate H(x) on local patches only
      procedure, PUBLIC :: clear    => PM_clear        ! clear obs data for this time step (keep init info)

   END type PM_type


   type :: PM_config_type

      character(len=16)  :: source_name
      character(len=32)  :: sensor_name
      character(len=16)  :: var_name
      real(r8)           :: fghz

   END type PM_config_type


   ! process a set of PM data 
   type :: PM_set_type

      ! different sensors and configurations in the set
      integer :: nsensor
      type(PM_type), allocatable :: sensors(:)
      type(PM_config_type), allocatable :: cfgs(:)

      ! concatenated obs and H(x) from all sensors for assimilation
      integer :: nobs
      real(r8), allocatable :: lat(:)
      real(r8), allocatable :: lon(:)
      real(r8), allocatable :: y(:)
      real(r8), allocatable :: r(:)
      real(r8), allocatable :: hx(:,:)

   CONTAINS

      procedure, PUBLIC :: init     => PM_set_init
      procedure, PUBLIC :: concat   => PM_set_concat
      procedure, PUBLIC :: calcg    => PM_set_calc_on_grid
      procedure, PUBLIC :: calcp    => PM_set_calc_on_pset
      procedure, PUBLIC :: clear    => PM_set_clear

   END type PM_set_type


!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE PM_init (this, pixelset, source_name, sensor_name, var_name, fghz)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(PM_type),     intent(inout)  :: this
      type(pixelset_type), intent(in)    :: pixelset     ! target pixel/patch set
      character(len=*),    intent(in)    :: source_name  ! observation source class, e.g., "PM"
      character(len=*),    intent(in)    :: sensor_name  ! sensor/subdirectory name
      character(len=*),    intent(in)    :: var_name     ! observation variable/file prefix
      real(r8),            intent(in)    :: fghz         ! frequency [GHz]

!-----------------------------------------------------------------------------
      CALL this%y%init(pixelset, source_name, sensor_name, var_name)
      this%fghz = fghz

   END SUBROUTINE PM_init

!-----------------------------------------------------------------------------

   SUBROUTINE PM_calc_on_pset (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(PM_type), intent(inout) :: this

!------------------------ Local Variables ------------------------------------
      integer :: iobs                         ! observation index
      integer :: iens                         ! ensemble member index
      integer :: np                           ! patch index
      integer :: ipart                        ! patch-grid intersection part index
      integer :: ib, jb                       ! block indices of an observation grid cell
      integer :: il, jl                       ! local indices inside the block

      real(r8) :: hx_tmp                      ! H(x) from one part-level RTM call
      real(r8) :: sumarea(0:DEF_DA_ENS_NUM)  ! valid part area accumulated for each ensemble

      logical, allocatable :: filter(:)       ! local land-patch mask

      type(block_data_real8_2d) :: theta_wgrid                ! incidence angle on obs world grid
      type(pointer_real8_1d), allocatable :: theta_part(:)    ! incidence angle on patch-grid parts

!-----------------------------------------------------------------------------
      IF (allocated(this%hx_pset)) deallocate(this%hx_pset)
      allocate(this%hx_pset(0:DEF_DA_ENS_NUM, numpatch))
      this%hx_pset = spval

      IF (this%y%nobs == 0 .or. numpatch == 0) RETURN

      ! Map incidence angle from obs grid cells to patch-grid intersection parts.
      IF (p_is_io) THEN
         CALL allocate_block_data(this%y%grid, theta_wgrid)
         CALL flush_block_data(theta_wgrid, spval)

         DO iobs = 1, this%y%nobs
            ib = this%y%grid%xblk(this%y%jj(iobs))
            jb = this%y%grid%yblk(this%y%ii(iobs))
            il = this%y%grid%xloc(this%y%jj(iobs))
            jl = this%y%grid%yloc(this%y%ii(iobs))
            IF (ib /= 0 .and. jb /= 0) THEN
               IF (gblock%pio(ib,jb) == p_iam_glb) THEN
                  theta_wgrid%blk(ib,jb)%val(il,jl) = this%y%theta(iobs)
               ENDIF
            ENDIF
         ENDDO
      ENDIF
      CALL this%y%mg2p%allocate_part(theta_part)
      CALL this%y%mg2p%grid2part(theta_wgrid, theta_part)

      ! Compute part-level H(x), then area-average valid parts to patch space.
      IF (p_is_worker) THEN
         allocate(filter(numpatch))
         filter(:) = patchtype(:) <= 2

         DO np = 1, numpatch
            IF (this%y%mg2p%npart(np) == 0) CYCLE
            IF (.not. filter(np)) CYCLE

            sumarea(:) = 0.0_r8

            DO ipart = 1, this%y%mg2p%npart(np)
               IF (theta_part(np)%val(ipart) == spval) CYCLE

               SELECT CASE (trim(this%y%var_name))
               CASE ('TB')
                  DO iens = 0, DEF_DA_ENS_NUM
                     CALL PM_RTM( &
                        patchtype(np), patchclass(np), dz_sno_ens(:,iens,np), &
                        forc_topo(np), htop(np), &
                        tref_ens(iens,np), t_soisno_ens(:,iens,np), tleaf_ens(iens,np), &
                        wliq_soisno_ens(:,iens,np), wice_soisno_ens(:,iens,np), &
                        snowdp_ens(iens,np), &
                        lai_ens(iens,np), sai_ens(iens,np), &
                        wf_clay(:,np), wf_sand(:,np), wf_silt(:,np), &
                        BD_all(:,np), porsl(:,np), &
                        theta_part(np)%val(ipart), this%fghz, &
                        hx_tmp)

                     IF (hx_tmp == spval) CYCLE
                     IF (sumarea(iens) == 0) this%hx_pset(iens,np) = 0
                     sumarea(iens) = sumarea(iens) + this%y%mg2p%areapart(np)%val(ipart)
                     this%hx_pset(iens,np) = this%hx_pset(iens,np) + &
                        hx_tmp * this%y%mg2p%areapart(np)%val(ipart)
                  ENDDO
               END SELECT
            ENDDO

            DO iens = 0, DEF_DA_ENS_NUM
               IF (sumarea(iens) > 0) this%hx_pset(iens,np) = this%hx_pset(iens,np) / sumarea(iens)
            ENDDO
         ENDDO

         deallocate(filter)
      ENDIF

      CALL this%y%mg2p%deallocate_part(theta_part)
      IF (p_is_io) THEN
         IF (allocated(theta_wgrid%blk)) deallocate(theta_wgrid%blk)
      ENDIF

   END SUBROUTINE PM_calc_on_pset

!-----------------------------------------------------------------------------

   SUBROUTINE PM_calc_on_grid (this, idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(PM_type),      intent(inout) :: this
      integer,             intent(in)    :: idate(3)
      real(r8),            intent(in)    :: deltim

!------------------------ Local Variables ------------------------------------
      integer :: iobs                         ! observation index in this%y
      integer :: ip                           ! IO task loop index on master
      integer :: ib, jb                       ! block indices of an observation grid cell
      integer :: il, jl                       ! local indices inside the block

      integer :: ndata                        ! valid obs count handled by one IO task
      integer :: isrc                         ! MPI source rank for received obs data
      integer :: smesg(2)                     ! IO-to-master message header: rank, ndata
      integer :: rmesg(2)                     ! master receive message header: rank, ndata

      real(r8) :: area_obs                    ! mapped patch area in one observation grid cell
      type(block_data_real8_2d) :: area_wgrid ! mapped patch area on obs world grid (this%y%grid)
      type(block_data_real8_3d) :: hx_wgrid   ! ensemble H(x) mapped to obs world grid (this%y%grid)

      integer,  allocatable :: iloc(:)        ! obs indices handled by one IO task
      integer,  allocatable :: tmp_idx(:)     ! obs indices received by master
      real(r8), allocatable :: tmp_data(:,:)  ! H(x) values sent/received by MPI
      logical,  allocatable :: filter(:)      ! local patch mask used for mapping
      integer :: nvalid_patch                 ! debug count of patches with valid ensemble H(x)

!-----------------------------------------------------------------------------
      CALL this%clear()

      ! read obs
      IF (p_is_master) WRITE(*,'(A,A)') '[CoLM-DA-DEBUG] PM before obs read sensor=', trim(this%y%sensor_name)
      CALL this%y%read(idate, deltim)
      IF (p_is_master) WRITE(*,'(A,A,A,I0)') '[CoLM-DA-DEBUG] PM after obs read sensor=', &
         trim(this%y%sensor_name), ' nobs=', this%y%nobs
      IF (this%y%nobs == 0) RETURN

      ! forward operator and map H(x) to y locations
      IF (p_is_master) WRITE(*,'(A,A)') '[CoLM-DA-DEBUG] PM before calcp sensor=', trim(this%y%sensor_name)
      CALL this%calcp()
      IF (p_is_master) WRITE(*,'(A,A)') '[CoLM-DA-DEBUG] PM after calcp sensor=', trim(this%y%sensor_name)

      ! Keep the mapped numerator and area denominator on the same valid patches.
      allocate (filter(numpatch))
      IF (p_is_worker) filter(:) = patchtype(:) <= 2 .and. &
         all(this%hx_pset(1:DEF_DA_ENS_NUM,:) /= spval, dim=1)
      IF (p_is_worker) THEN
         nvalid_patch = count(filter)
         WRITE(*,'(A,I0,A,A,A,I0)') '[CoLM-DA-DEBUG] worker=', p_iam_glb, &
            ' PM valid hx patches sensor=', trim(this%y%sensor_name), ' count=', nvalid_patch
      ENDIF

      ! mapping H(x) from patch to y grid
      IF (p_is_io) THEN
         CALL allocate_block_data(this%y%grid, area_wgrid)
         CALL allocate_block_data(this%y%grid, hx_wgrid, DEF_DA_ENS_NUM)
      ENDIF
      CALL this%y%mg2p%get_sumarea(area_wgrid, filter)
      CALL this%y%mg2p%pset2grid(this%hx_pset(1:DEF_DA_ENS_NUM,:), hx_wgrid, spv=spval, msk=filter)
      deallocate (filter)

      ! crop H(x) corresponding to valid y locations
      allocate (this%hx(this%y%nobs, DEF_DA_ENS_NUM))

      IF (p_is_io) THEN
         allocate (iloc(this%y%nobs))
         this%hx = spval
         WRITE(*,'(A,I0,A,A)') '[CoLM-DA-DEBUG] io=', p_iam_glb, &
            ' PM crop hx start sensor=', trim(this%y%sensor_name)

         ndata = 0
         DO iobs = 1, this%y%nobs
            ib = this%y%grid%xblk(this%y%jj(iobs))
            jb = this%y%grid%yblk(this%y%ii(iobs))
            il = this%y%grid%xloc(this%y%jj(iobs))
            jl = this%y%grid%yloc(this%y%ii(iobs))
            IF (ib /= 0 .and. jb /= 0) THEN
               IF (gblock%pio(ib,jb) == p_iam_glb) THEN
                  area_obs = area_wgrid%blk(ib,jb)%val(il,jl)
                  IF (area_obs > 0 .and. &
                     all(hx_wgrid%blk(ib,jb)%val(:,il,jl) /= spval)) THEN
#ifdef USEMPI
                     ndata = ndata + 1
                     iloc(ndata) = iobs
                     this%hx(ndata,:) = hx_wgrid%blk(ib,jb)%val(:,il,jl) / area_obs
#endif
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

#ifdef USEMPI
         smesg = (/p_iam_glb, ndata/)
         WRITE(*,'(A,I0,A,I0,A,A)') '[CoLM-DA-DEBUG] io=', p_iam_glb, &
            ' PM send header ndata=', ndata, ' sensor=', trim(this%y%sensor_name)
         CALL mpi_send(smesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)

         IF (ndata > 0) THEN
            CALL mpi_send(iloc(1:ndata), ndata, MPI_INTEGER, p_address_master, mpi_tag_data, p_comm_glb, p_err)
            allocate (tmp_data(ndata, DEF_DA_ENS_NUM))
            tmp_data = this%hx(1:ndata,:)
            CALL mpi_send(tmp_data, ndata*DEF_DA_ENS_NUM, MPI_REAL8, p_address_master, mpi_tag_data+1, p_comm_glb, p_err)
            deallocate (tmp_data)
         ENDIF
#endif
         deallocate (iloc)
         IF (allocated(hx_wgrid%blk))   deallocate (hx_wgrid%blk)
         IF (allocated(area_wgrid%blk)) deallocate (area_wgrid%blk)
      ENDIF

#ifdef USEMPI
      IF (p_is_master) THEN
         this%hx = spval
         WRITE(*,'(A,I0,A,A)') '[CoLM-DA-DEBUG] master waiting PM io headers count=', &
            p_np_io, ' sensor=', trim(this%y%sensor_name)

         DO ip = 0, p_np_io - 1
            CALL mpi_recv(rmesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
            WRITE(*,'(A,I0,A,I0,A,A)') '[CoLM-DA-DEBUG] master got PM header source=', &
               rmesg(1), ' ndata=', rmesg(2), ' sensor=', trim(this%y%sensor_name)

            ndata = rmesg(2)
            IF (ndata > 0) THEN
               allocate (tmp_idx(ndata))
               allocate (tmp_data(ndata, DEF_DA_ENS_NUM))

               isrc = rmesg(1)
               CALL mpi_recv(tmp_idx, ndata, MPI_INTEGER, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL mpi_recv(tmp_data, ndata*DEF_DA_ENS_NUM, MPI_REAL8, isrc, mpi_tag_data+1, p_comm_glb, p_stat, p_err)
               this%hx(tmp_idx,:) = tmp_data

               deallocate (tmp_idx)
               deallocate (tmp_data)
            ENDIF
         ENDDO
         WRITE(*,'(A,A)') '[CoLM-DA-DEBUG] master before PM hx bcast sensor=', trim(this%y%sensor_name)
      ENDIF
      CALL mpi_bcast(this%hx, this%y%nobs*DEF_DA_ENS_NUM, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      IF (p_is_master) WRITE(*,'(A,A)') '[CoLM-DA-DEBUG] master after PM hx bcast sensor=', trim(this%y%sensor_name)
#endif

   END SUBROUTINE PM_calc_on_grid

!-----------------------------------------------------------------------------

   SUBROUTINE PM_clear (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(PM_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      CALL this%y%clear()

      IF (allocated(this%hx_pset)) deallocate(this%hx_pset)
      IF (allocated(this%hx))      deallocate(this%hx)

   END SUBROUTINE PM_clear

!-----------------------------------------------------------------------------

   SUBROUTINE PM_set_init (this, pixelset, configs)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(PM_set_type),   intent(inout) :: this
      type(pixelset_type),  intent(in)    :: pixelset
      type(PM_config_type), intent(in)    :: configs(:)

      integer :: i

      IF (allocated(this%sensors)) THEN
         DO i = 1, size(this%sensors)
            CALL this%sensors(i)%clear()
         ENDDO
         deallocate(this%sensors)
      ENDIF
      IF (allocated(this%cfgs)) deallocate(this%cfgs)

      this%nsensor = 0
      IF (size(configs) == 0) return

      this%nsensor = size(configs)

      allocate(this%cfgs(size(configs)))
      allocate(this%sensors(size(configs)))
      DO i = 1, size(configs)
         this%cfgs(i) = configs(i)
         CALL this%sensors(i)%init( &
            pixelset, &
            this%cfgs(i)%source_name, &
            this%cfgs(i)%sensor_name, &
            this%cfgs(i)%var_name, &
            this%cfgs(i)%fghz)
      ENDDO

   END SUBROUTINE PM_set_init

!-----------------------------------------------------------------------------

   SUBROUTINE PM_set_clear (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(PM_set_type), intent(inout) :: this
      integer :: is

!-----------------------------------------------------------------------------
      this%nobs = 0
      IF (allocated(this%lat))       deallocate(this%lat)
      IF (allocated(this%lon))       deallocate(this%lon)
      IF (allocated(this%y))         deallocate(this%y)
      IF (allocated(this%r))         deallocate(this%r)
      IF (allocated(this%hx))        deallocate(this%hx)

      IF (.not. allocated(this%sensors)) return

      DO is = 1, this%nsensor
         CALL this%sensors(is)%clear()
      ENDDO

   END SUBROUTINE PM_set_clear

!-----------------------------------------------------------------------------

   SUBROUTINE PM_set_concat (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(PM_set_type), intent(inout) :: this
      integer :: is, n, offset

!-----------------------------------------------------------------------------
      ! get total number of obs and allocate arrays
      DO is = 1, this%nsensor
         IF (this%sensors(is)%y%nobs == 0) CYCLE
         IF (.not. allocated(this%sensors(is)%hx)) CYCLE
         this%nobs = this%nobs + this%sensors(is)%y%nobs
      ENDDO
      IF (this%nobs == 0) return

      allocate(this%lat(this%nobs))
      allocate(this%lon(this%nobs))
      allocate(this%y(this%nobs))
      allocate(this%r(this%nobs))
      allocate(this%hx(this%nobs, DEF_DA_ENS_NUM))

      ! concatenate H(x) and y from each sensor to the set-level arrays
      offset = 0
      DO is = 1, this%nsensor
         IF (.not. allocated(this%sensors(is)%hx)) CYCLE

         n = this%sensors(is)%y%nobs
         IF (n == 0) CYCLE

         this%lat(offset+1:offset+n) = this%sensors(is)%y%lat
         this%lon(offset+1:offset+n) = this%sensors(is)%y%lon
         this%y(offset+1:offset+n) = this%sensors(is)%y%y
         this%r(offset+1:offset+n) = this%sensors(is)%y%r
         this%hx(offset+1:offset+n,:) = this%sensors(is)%hx

         offset = offset + n
      ENDDO

   END SUBROUTINE PM_set_concat

!-----------------------------------------------------------------------------

   SUBROUTINE PM_set_calc_on_grid(this, idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(PM_set_type),  intent(inout) :: this
      integer,             intent(in)    :: idate(3)
      real(r8),            intent(in)    :: deltim

      integer :: is

!-----------------------------------------------------------------------------
      CALL this%clear()

      DO is = 1, this%nsensor
         CALL this%sensors(is)%calcg(idate, deltim)
      ENDDO

      CALL this%concat()

   END SUBROUTINE PM_set_calc_on_grid

!-----------------------------------------------------------------------------

   SUBROUTINE PM_set_calc_on_pset(this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(PM_set_type), intent(inout) :: this
      integer :: is

!-----------------------------------------------------------------------------
      DO is = 1, this%nsensor
         CALL this%sensors(is)%calcp()
      ENDDO

   END SUBROUTINE PM_set_calc_on_pset

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Obs_PM
#endif
