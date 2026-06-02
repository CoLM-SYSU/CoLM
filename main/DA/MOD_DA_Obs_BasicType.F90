#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Obs_BasicType
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Basic observation derived types for data assimilation.
!
!    This module defines common observation containers and utilities shared by
!    source-specific observation modules. Observation sources are organized into
!    two broad classes:
!
!      1. Grid observations:
!         Satellite or gridded products with their own spatial grid. A
!         grid_obs_type stores the observation grid, the grid-to-patch spatial
!         mapping, and timestep-dependent observation vectors.
!
!      2. In-situ observations:
!         Point-based measurements such as stations or field networks. These
!         observations do not own a source grid, but still need common metadata,
!         location, time, value, and error containers.
!
!    Higher-level modules, such as PM, AM, GNOS, and soil-moisture
!    assimilation wrappers, should build on these basic types instead of
!    duplicating storage and search logic.
!
! SUPPORTED:
! +--------+----------+----------------------+-------------------------------+
! | Source | Variable | Assimilated state    | Example datasets              |
! +--------+----------+----------------------+-------------------------------+
! | PM     | TB       | Soil moisture, soil  | SMAP L1C, FY3D MWRI,          |
! |        |          | temperature, snow    | AMSR-E.                       |
! |        |          | variables            |                               |
! +--------+----------+----------------------+-------------------------------+
! | PM     | SM       | Soil moisture        | SMAP L2 SM, SMOS L2 SM,       |
! |        |          |                      | AMSR2 SM                      |
! +--------+----------+----------------------+-------------------------------+
! | AM     | TB       | Soil moisture, soil  | ATMS/MWHS/MWTS TB             |
! |        |          | temperature, snow    |                               |
! |        |          | variables            |                               |
! +--------+----------+----------------------+-------------------------------+
! | AM     | SM       | Soil moisture        | ASCAT, Sentinel-1             |
! +--------+----------+----------------------+-------------------------------+
! | GNOS   | DDM      | Soil moisture        | FY3E, CYGNSS                  |
! +--------+----------+----------------------+-------------------------------+
! | MR     | LST      | Land surface         | MODIS LST, VIIRS LST, FY3F    |
! |        |          | temperature, soil    |                               |
! |        |          | temperature          |                               |
! +--------+----------+----------------------+-------------------------------+
! | MR     | fSCA     | Snow variables       | MODIS, VIIRS, FY3F            |
! +--------+----------+----------------------+-------------------------------+
!   PM   = passive microwave remote sensing
!   AM   = active microwave remote sensing
!   GNOS = GNSS-R observation source
!   MR   = moderate-resolution optical/thermal remote sensing
!   TB   = brightness temperature
!   SM   = retrieved soil moisture
!   DDM  = delay-Doppler map
!   LST  = land surface temperature
!   fSCA = fractional snow-covered area
!
! HISTORY:
!    Lu Li, 05/2026: Initial version
!-----------------------------------------------------------------------------
    USE MOD_Precision
    USE MOD_SPMD_Task
    USE MOD_Grid
    USE MOD_Namelist
    USE MOD_NetCDFSerial
    USE MOD_Pixelset
    USE MOD_Vars_Global, ONLY: spval
    USE MOD_SpatialMapping
    USE MOD_TimeManager
    IMPLICIT NONE

    type :: grid_obs_type

        ! initialize once
        character(len=16) :: source_name = ''        ! "PM", "PM", "AM", "AM", "GNOS", "MR", "MR"
        character(len=32) :: sensor_name = ''        ! sensor/subdirectory name, e.g., "SMAP_L1C_D"
        character(len=16) :: var_name    = ''        ! "TB", "SM", "TB", "SM", "DDM", "LST", "fSCA"
        type(grid_type)   :: grid                    ! grid info for obs
        type(spatial_mapping_type) :: mg2p           ! mapping from obs grid to patch

        ! counting variables for this timestep
        integer           :: nobs                    ! total avaliable obs

        ! per-observation info (length = nobs)
        real(r8), allocatable :: lat(:)              ! latitude of each obs point
        real(r8), allocatable :: lon(:)              ! longitude of each obs point
        real(r8), allocatable :: time(:)             ! UTC time of each obs point (seconds since begin of current day)
        real(r8), allocatable :: y(:)                ! observed value
        real(r8), allocatable :: r(:)                ! observation error variance
        real(r8), allocatable :: theta(:)            ! incidence angle for each obs point [rad]
        integer,  allocatable :: ii(:)               ! the number of lat grids of this cell
        integer,  allocatable :: jj(:)               ! the number of lon grids of this cell

    CONTAINS

        procedure, PUBLIC :: init  => grid_obs_init  ! initialize the obs type (grid info, mapping, etc.)
        procedure, PUBLIC :: read  => grid_obs_read  ! read obs data for this time step (from file)
        procedure, PUBLIC :: clear => grid_obs_clear ! clear obs data for this time step (keep init info)

    END type grid_obs_type


!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

    SUBROUTINE grid_obs_init (this, pixelset, source_name, sensor_name, var_name)

!-----------------------------------------------------------------------------
        IMPLICIT NONE

        class(grid_obs_type),  intent(inout) :: this
        type(pixelset_type),   intent(in)    :: pixelset
        character(len=*),      intent(in)    :: source_name
        character(len=*),      intent(in)    :: sensor_name
        character(len=*),      intent(in)    :: var_name

        character(len=256) :: grid_file_name

        this%source_name = source_name
        this%sensor_name = sensor_name
        this%var_name = var_name
        this%nobs = 0

        grid_file_name = trim(DEF_DA_OBS_DIR)//'/'//trim(this%sensor_name)//'/grid.nc'
        CALL this%grid%define_from_file(grid_file_name, 'latitude', 'longitude')
        CALL this%mg2p%build_arealweighted(this%grid, pixelset)

    END SUBROUTINE grid_obs_init

!-----------------------------------------------------------------------------

    SUBROUTINE grid_obs_clear (this)

!-----------------------------------------------------------------------------
        IMPLICIT NONE

        class(grid_obs_type),   intent(inout) :: this

        IF (allocated(this%lat))  deallocate(this%lat )
        IF (allocated(this%lon))  deallocate(this%lon )
        IF (allocated(this%time)) deallocate(this%time)
        IF (allocated(this%y))    deallocate(this%y   )
        IF (allocated(this%r))    deallocate(this%r   )
        IF (allocated(this%theta)) deallocate(this%theta)
        IF (allocated(this%ii))   deallocate(this%ii  )
        IF (allocated(this%jj))   deallocate(this%jj  )

        this%nobs = 0

    END SUBROUTINE grid_obs_clear

!-----------------------------------------------------------------------------

    SUBROUTINE grid_obs_read (this, idate, deltim)

!-----------------------------------------------------------------------------
        IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
        class(grid_obs_type), intent(inout) :: this
        integer,              intent(in)    :: idate(3)
        real(r8),             intent(in)    :: deltim

!------------------------ Local Variables ------------------------------------
        logical :: exists_file                 ! whether the expected obs file exists
        integer :: sdate(3)                    ! timestep start time adjusted to the beginning of step
        integer :: month, mday, hour           ! calendar month, day of month, and UTC hour
        character(len=4)   :: yearstr          ! YYYY string for file name
        character(len=2)   :: monthstr         ! MM string for file name
        character(len=2)   :: daystr           ! DD string for file name
        character(len=2)   :: hourstr          ! HH string for file name
        character(len=256) :: obs_file_name    ! full path of current timestep obs file
        real(r8), allocatable :: lat(:)        ! raw obs latitudes read from file [deg]
        real(r8), allocatable :: lon(:)        ! raw obs longitudes read from file [deg]
        real(r8), allocatable :: time(:)       ! raw obs times read from file [s since current UTC day]
        real(r8), allocatable :: y(:)          ! raw obs values for var
        real(r8), allocatable :: r(:)          ! raw obs error variance for var
        real(r8), allocatable :: theta(:)      ! raw incidence angles read from file [rad]
        integer,  allocatable :: ii(:)         ! raw obs ii index read from file
        integer,  allocatable :: jj(:)         ! raw obs jj index read from file
        logical,  allocatable :: valid(:)      ! full-length mask for obs inside time window and domain
        integer :: nbad_r                      ! debug count of invalid observation error variances
        integer :: nbad_idx                    ! debug count of invalid grid indices in selected obs
        integer :: nbad_theta                  ! debug count of invalid incidence angles in selected obs

!-----------------------------------------------------------------------------

        CALL this%clear()

        ! UTC begin time of the current DA window
        sdate = idate
        CALL backtime(deltim, sdate)
        CALL adj2begin(sdate)
        CALL julian2monthday(sdate(1), sdate(2), month, mday)
        hour = int(sdate(3)/3600)
        write (yearstr, '(I4.4)') sdate(1)
        write (monthstr, '(I2.2)') month
        write (daystr, '(I2.2)') mday
        write (hourstr, '(I2.2)') hour

        ! whether exists obs for this timestep at domain
        obs_file_name = trim(DEF_DA_OBS_DIR)//'/'//trim(this%sensor_name)//'/'//trim(this%var_name)//'_'// &
            trim(yearstr)//'_'//trim(monthstr)//'_'//trim(daystr)//'_'//trim(hourstr)//'.nc'
        inquire (file=trim(obs_file_name), exist=exists_file)
        IF (p_is_master) WRITE(*,'(A,A,A,L1)') '[CoLM-DA-DEBUG] obs file=', &
            trim(obs_file_name), ' exists=', exists_file

        IF (exists_file) THEN
            CALL ncio_read_bcast_serial(obs_file_name, 'latitude', lat)
            CALL ncio_read_bcast_serial(obs_file_name, 'longitude', lon)
            CALL ncio_read_bcast_serial(obs_file_name, 'time', time)
            CALL ncio_read_bcast_serial(obs_file_name, trim(this%var_name), y)
            CALL ncio_read_bcast_serial(obs_file_name, 'r', r)
            CALL ncio_read_bcast_serial(obs_file_name, 'theta', theta)
            CALL ncio_read_bcast_serial(obs_file_name, 'ii', ii)
            CALL ncio_read_bcast_serial(obs_file_name, 'jj', jj)

            allocate (valid(size(time)))
            valid = time >= sdate(3)          .and. &
                    time <= sdate(3) + deltim .and. &
                    lat  >= DEF_domain%edges  .and. &
                    lat  <= DEF_domain%edgen  .and. &
                    lon  >= DEF_domain%edgew  .and. &
                    lon  <= DEF_domain%edgee
            this%nobs = count(valid)
            IF (p_is_master) WRITE(*,'(A,A,A,I0,A,I0)') '[CoLM-DA-DEBUG] obs read sensor=', &
                trim(this%sensor_name), ' raw_nobs=', size(time), ' selected_nobs=', this%nobs

            IF (this%nobs > 0) THEN
                allocate (this%lat(this%nobs))
                allocate (this%lon(this%nobs))
                allocate (this%time(this%nobs))
                allocate (this%y(this%nobs))
                allocate (this%r(this%nobs))
                allocate (this%theta(this%nobs))
                allocate (this%ii(this%nobs))
                allocate (this%jj(this%nobs))

                this%lat    = pack(lat, valid)
                this%lon    = pack(lon, valid)
                this%time   = pack(time, valid)
                this%y      = pack(y, valid)
                this%r      = pack(r, valid)
                this%theta  = pack(theta, valid)
                this%ii     = pack(ii, valid)
                this%jj     = pack(jj, valid)

                nbad_r = count(this%r <= 0.0_r8 .or. this%r == spval)
                nbad_idx = count(this%ii <= 0 .or. this%jj <= 0)
                nbad_theta = count(this%theta == spval)
                IF (p_is_master) WRITE(*,'(A,I0,A,ES12.4,A,ES12.4,A,I0,A,I0,A,I0)') &
                    '[CoLM-DA-DEBUG] obs selected stats nobs=', this%nobs, &
                    ' minR=', minval(this%r), ' maxR=', maxval(this%r), &
                    ' badR=', nbad_r, ' badIdx=', nbad_idx, ' badTheta=', nbad_theta
            ENDIF

            IF (allocated(valid)) deallocate(valid)
            IF (allocated(lat))   deallocate(lat)
            IF (allocated(lon))   deallocate(lon)
            IF (allocated(time))  deallocate(time)
            IF (allocated(y))     deallocate(y)
            IF (allocated(r))     deallocate(r)
            IF (allocated(theta)) deallocate(theta)
            IF (allocated(ii))    deallocate(ii)
            IF (allocated(jj))    deallocate(jj)
        ENDIF

    END SUBROUTINE grid_obs_read


!-----------------------------------------------------------------------------
END MODULE MOD_DA_Obs_BasicType
#endif
