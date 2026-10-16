MODULE CMF_CTRL_SED_MOD
  !==========================================================
  !* PURPOSE: CaMa-Flood sediment transport scheme
  ! (Mimicking CMF_CTRL_TRACER_MOD structure)
  !
  ! (C) M.Hatono (Hiroshima-U) May 2021 (original parts)
  !
  ! Licensed under the Apache License, Version 2.0 (the "License");
  !   You may not use this file except in compliance with the License.
  !   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
  !==========================================================
  USE PARKIND1,                ONLY: JPIM, JPRB, JPRM, JPRD
  USE CMF_UTILS_MOD,           ONLY: INQUIRE_FID, CMF_CheckNanB, mapR2vecD, vecD2mapR, CONV_END, &
                                     mapP2vecP, vecP2mapP, vecP2mapR
  USE YOS_CMF_INPUT,           ONLY: LOGNAM, TMPNAM, NX, NY, NXIN, NYIN, INPN, RMIS, DT, PGRV, PMANRIV
  USE YOS_CMF_INPUT,           ONLY: IFRQ_OUT, CSUFBIN, CSUFVEC, LSPAMAT, LOUTPUT, CSETFILE, NSETFILE
  USE YOS_CMF_MAP,             ONLY: NSEQMAX, NSEQALL, NSEQRIV, NPTHOUT, REGIONTHIS
  USE YOS_CMF_MAP,             ONLY: I1NEXT,  PTH_UPST,PTH_DOWN, I2MASK,D2GRAREA
  USE YOS_CMF_MAP,             ONLY: I1UPST,   I1UPN,    I1P_OUT,  I1P_OUTN, I1P_INF, I1P_INFN
  USE YOS_CMF_MAP,             ONLY: INPX, INPY, INPA, D2RIVLEN, D2RIVWTH
  USE YOS_CMF_PROG,            ONLY: P2RIVSTO, P2FLDSTO,D2RIVOUT
  USE YOS_CMF_DIAG,            ONLY: D2RIVDPH, D2RIVVEL, &
                                     D2OUTFLW_aAVG, D1PTHFLWSUM_aAVG
  USE YOS_CMF_TIME,            ONLY: KSTEP, NSTEPS, JYYYYMMDD, JHHMM, JDD, JHOUR, JMIN, &
                                     IYYYY, IMM, IDD, IHOUR, IMIN, &
                                     ISYYYY, ISMM, ISDD, ISHOUR, ISMIN ! For output/restart time naming
 !use CMF_CTRL_OUTPUT_MOD,     only: COUTDIR, IRECOUT, LOUTCDF, LOUTVEC, NDLEVEL
  use YOS_CMF_DIAG,            only: D2FLDFRC
  use YOS_CMF_INPUT,           only: TMPNAM,NXIN,NYIN,DTIN
  use YOS_CMF_INPUT,           only: NLFP
  use YOS_CMF_TIME,            only: KMINSTAIN, KMINSTART, KMINEND

#ifdef UseMPI_CMF
  USE CMF_CTRL_MPI_MOD,        ONLY: CMF_MPI_AllReduce_P2MAP, CMF_MPI_AllReduce_R2MAP ! From tracer
  ! Note: Original tracer also had "USE MPI", but CMF_CTRL_MPI_MOD might be a wrapper
#endif

  IMPLICIT NONE

  ! Declare CMF_ERROR_ABORT_STRING subroutine
  PRIVATE :: CMF_ERROR_ABORT_STRING
  PUBLIC  :: CMF_ERROR_ABORT
  PUBLIC  :: CMF_SED_DIAG_AVEMAX_ADPSTP
  PUBLIC  :: CMF_SED_RESTART_WRITE
  ! Define CMF_ERROR_ABORT interface
  INTERFACE CMF_ERROR_ABORT
    MODULE PROCEDURE CMF_ERROR_ABORT_STRING
  END INTERFACE CMF_ERROR_ABORT

  !----------------------------------------------------------------------
  ! Sediment Parameters (to be read from NAMSED or set as default)
  !----------------------------------------------------------------------
  INTEGER(KIND=JPIM)              :: nsed             ! number of sediment particle size classes
  INTEGER(KIND=JPIM)              :: totlyrnum        ! number of deposition layers under active layer
  LOGICAL                         :: revEgia         ! if true, use Egiazoroff equation
  REAL(KIND=JPRB)                 :: lambda           ! porosity
  REAL(KIND=JPRB)                 :: lyrdph          ! active layer depth (m)
  REAL(KIND=JPRB)                 :: psedD           ! density of sediment (e.g., 2.65 g/cm3)
  REAL(KIND=JPRB)                 :: pset       ! parameter for setting velocity
  REAL(KIND=JPRB)                 :: pwatD           ! density of water (e.g., 1.0 g/cm3)
  REAL(KIND=JPRB)                 :: visKin          ! kinematic viscosity of water (m2/s)
  REAL(KIND=JPRB)                 :: vonKar          ! von Karman coefficient

  CHARACTER(LEN=256)              :: crocdph    ! filename for bedrock depth data
  CHARACTER(LEN=256)              :: sedD     ! comma-separated string for sediment diameters
  CHARACTER(LEN=256)              :: csedfrc   ! filename for initial sediment layer thickness/fractions


  CHARACTER(LEN=256)              :: cslope            !  
  CHARACTER(LEN=256)              :: cinpmat_sed  ! input matrix for sediment
  REAL(KIND=JPRB)                 :: dsylunit         ! unit conversion for sediment [m3/km2] -> [m3/m2]
  REAL(KIND=JPRB)                 :: pyld, pyldc, pyldpc  ! parameters for sediment erosion calculation

  !*** namelist/sediment_output
  CHARACTER(LEN=256)              :: csedsout

  ! Restart controls (similar to tracer)
  LOGICAL                         :: LRESTDBL_SED     ! true: binary sediment restart in double precision
  INTEGER(KIND=JPIM)              :: IFRQ_RST_SED     ! sediment restart frequency (0: end, N: hourly, 30: monthly)

  ! Output controls (similar to tracer)
  CHARACTER(LEN=256)              :: COUTDIR_SED      ! sediment output directory
  CHARACTER(LEN=256)              :: COUTTAG_SED      ! sediment output tag name for each experiment
  LOGICAL                         :: LOUTVEC_SED      ! true for vectorial sediment output

  CHARACTER(LEN=256)              :: fName_sed

  
  !----------------------------------------------------------------------
  ! Sediment Variables (Arrays)
  !----------------------------------------------------------------------
  ! Input / Map data
  !REAL(KIND=JPRB), ALLOCATABLE     :: D2SEDFRC(:,:)    ! Sediment dist fraction in active layer [NSEQMAX, nsed]
  REAL(KIND=JPRB), ALLOCATABLE     :: sDiam(:)        ! Sediment grain diameter (m) [nsed]
  REAL(KIND=JPRB), ALLOCATABLE     :: SETVEL(:)        ! Settling velocity (m/s) [nsed]
  REAL(KIND=JPRB), ALLOCATABLE     :: D2ROCDPH(:)      ! Bedrock depth (m) [NSEQMAX] (placeholder if needed)

  ! Variables for time integration and diagnostics
  REAL(KIND=JPRB)                 :: SADD_RIV         ! Accumulator for time in DT for averaging river vars
  REAL(KIND=JPRB)                 :: SADD_OUT         ! Accumulator for time in IFRQ_OUT for averaging output vars (like NADD_out in tracer)

  REAL(KIND=JPRB), ALLOCATABLE     :: D2RIVOUT_SED(:)  ! River discharge averaged over DT [NSEQMAX]
  REAL(KIND=JPRB), ALLOCATABLE     :: D2RIVVEL_SED(:)  ! River velocity averaged over DT [NSEQMAX]
  REAL(KIND=JPRB), ALLOCATABLE     :: D2RIVSTO_PRE(:)  ! River storage from previous DT timestep [NSEQMAX]

  ! Core sediment state variables
  !REAL(KIND=JPRB), ALLOCATABLE     :: d2layer(:,:) ! Active layer sediment volume storage (m3) [NSEQMAX, nsed] (use JPRD for storage)
  !REAL(KIND=JPRD), ALLOCATABLE     :: P2SEDDEP_STO(:,:,:)! Sub-surface deposited sediment volume (m3) [NSEQMAX, totlyrnum, nsed] (use JPRD for storage)

  ! Fluxes and intermediate calculation variables (typically JPRB, as in tracer)
  !REAL(KIND=JPRB), ALLOCATABLE     :: d2sedout(:,:)    ! Suspended sediment outflow rate (m3/s) [NSEQMAX, nsed]
  !REAL(KIND=JPRB), ALLOCATABLE     :: d2sedcon(:,:)    ! Suspended sediment concentration (vol/vol) [NSEQMAX, nsed]
  !REAL(KIND=JPRB), ALLOCATABLE     :: d2sedinp(:,:) ! External sediment inflow source (mass/area/time or vol/area/time) [NSEQMAX, nsed]
  !REAL(KIND=JPRB), ALLOCATABLE     :: d2bedout(:,:)    ! Bedload sediment transport rate (m3/s) [NSEQMAX, nsed]
  !REAL(KIND=JPRB), ALLOCATABLE     :: d2netflw(:,:)    ! Net flux between suspension and bed (entrainment - deposition) (m3/s) [NSEQMAX, nsed]
  !REAL(KIND=JPRB), ALLOCATABLE       :: d2seddep(:,:,:)  ! deposition storage [NSEQMAX, totlyrnum, nsed]
  
  real(kind=JPRB),allocatable,target :: d2sedv(:,:,:)    ! storage array for sediment variables
  real(kind=JPRB),pointer            :: d2bedout(:,:)    ! bedflow (m3/s)
  real(kind=JPRB),pointer            :: d2layer(:,:)     ! exchange layer storage (m3)
  real(kind=JPRB),pointer            :: d2netflw(:,:)    ! suspension - deposition (m3/s)
  real(kind=JPRB),pointer            :: d2sedcon(:,:)    ! suspended sediment concentration (m3/m3)
  real(kind=JPRB),pointer            :: d2sedfrc(:,:)    ! sediment distribution fraction [-]
  real(kind=JPRB),pointer            :: d2sedout(:,:)    ! suspended sediment flow (m3/s)
  real(kind=JPRB),pointer            :: d2sedinp(:,:)    ! sediment inflow (m3/s)
  
  real(kind=JPRB),allocatable,target :: d2depv(:,:,:)    ! storage array for sediment variables
  real(kind=JPRB),pointer            :: d2seddep(:,:,:)  ! deposition storage
  
  real(KIND=JPRB),allocatable,target :: d2sedv_avg(:,:,:)    ! storage array for averaged sediment variables
  real(KIND=JPRB),pointer            :: d2bedout_avg(:,:)  ! bedflow (m3/s)
  real(KIND=JPRB),pointer            :: d2netflw_avg(:,:)    ! suspension - deposition (m3/s)
  real(KIND=JPRB),pointer            :: d2sedout_avg(:,:)    ! suspended sediment flow (m3/s)
  real(KIND=JPRB),pointer            :: d2sedinp_avg(:,:)    ! sediment inflow (m3/s)

  

  ! Buffer for forcing input (similar to TBUFF in tracer)
  REAL(KIND=JPRB), ALLOCATABLE       :: SBUFF(:,:,:)     ! Buffer to store forcing sediment input [NXIN, NYIN, nsed]
  real(KIND=JPRB),allocatable        :: d2slope(:,:)     ! floodplain slope [deg]
  integer(kind=JPIM)                 :: iseq


  ! Averaged variables for output (typically JPRB)
  ! For output file handling 
  INTEGER(KIND=JPIM)              :: NVARSOUT
  INTEGER(KIND=JPIM), PARAMETER   :: nvars=30
  CHARACTER(LEN=256) :: cvnames(nvars)

  TYPE TVAROUT_SED
    CHARACTER(LEN=256)              :: CVNAME             ! output variable name
    CHARACTER(LEN=256)              :: CVLNAME            ! output variable long name
    CHARACTER(LEN=256)              :: CVUNITS            ! output units
    CHARACTER(LEN=256)              :: CFILE              ! output full path file name 
    INTEGER(KIND=JPIM)              :: BINID              ! output binary output file ID
    INTEGER(KIND=JPIM)              :: NCID               ! output netCDF output file ID
    INTEGER(KIND=JPIM)              :: VARID              ! output netCDF output variable ID
    INTEGER(KIND=JPIM)              :: TIMID              ! output netCDF time   variable ID 
    INTEGER(KIND=JPIM)              :: IRECNC               ! Current time record for writting 
  END TYPE TVAROUT_SED 

  TYPE(TVAROUT_SED), ALLOCATABLE :: varout_sed(:)
 
  ! For restart file handling
  CHARACTER(LEN=256)              :: sedrest_infile, sedrest_outpre

  ! External function declarations
  PRIVATE :: CALC_SETTLINGVELOCITY

  NAMELIST /NAMSED/ nsed, totlyrnum, revEgia, lambda, lyrdph, psedD, pset, pwatD, visKin, vonKar, &
                    crocdph, sedD, csedfrc, &
                    cslope,dsylunit, pyld, pyldc, pyldpc,    &
                    cinpmat_sed,csedsout,&
                    sedrest_infile, sedrest_outpre, ifrq_rst_sed,&
                    COUTDIR_SED,COUTTAG_SED
CONTAINS

  ! Add CMF_ERROR_ABORT implementation
  SUBROUTINE CMF_ERROR_ABORT_STRING(CMESSAGE)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: CMESSAGE
    
    WRITE(LOGNAM,*) 'ERROR: ', TRIM(CMESSAGE)
    STOP 'CMF_ERROR_ABORT'
  END SUBROUTINE CMF_ERROR_ABORT_STRING

  SUBROUTINE CMF_SED_NMLIST
    !==========================================================
    !* PURPOSE: Read sediment namelist
    !==========================================================
    IMPLICIT NONE

    ! Local variables
    INTEGER :: IERR, IU
    LOGICAL :: LEXIST

    WRITE(LOGNAM,*) 'CMF_SED_NMLIST: READING SEDIMENT NAMELIST'

    ! Set Default Values
    lambda           = 0.4_JPRB
    lyrdph           = 0.00005d0    ! Default active layer depth 10cm
    nsed             = 3
    psedD            = 2.65d0       !  
    pset             = 1.d0     ! typical value for natural sand
    pwatD            = 1.d0
    totlyrnum        = 5
    revEgia          = .true.
    visKin           = 1.0E-6_JPRB  ! m2/s at 20degC
    vonKar           = 0.4_JPRB
 
    crocdph          = ""
    sedD             = "NONE"      ! Default to 1mm sand
    csedfrc          = ""

    !input
    
    cslope      ='./slope.nc'
    dsylunit    = 1.d-6
    pyld        = 0.01d0
    pyldc       = 2.d0
    pyldpc      = 2.d0
    cinpmat_sed = './inpmat.bin'

   ! output
    cvnames(:) = 'none'
    csedsout   ='./'
    
    ! restart
    sedrest_infile=''
    sedrest_outpre='./'
    ifrq_rst_sed=0

    NSETFILE=INQUIRE_FID()
    OPEN(NSETFILE,FILE=CSETFILE,STATUS="OLD")
    WRITE(LOGNAM,*) "CMF::FORCING_NMLIST: namelist OPEN in unit: ", TRIM(CSETFILE), NSETFILE 

!*** 3. read namelist
    REWIND(NSETFILE)
    READ(NSETFILE,NML=NAMSED)

    WRITE(LOGNAM,*) 'CMF_SED_NMLIST: --- Sediment Namelist Parameters ---'
    WRITE(LOGNAM,*) 'nsed            = ', nsed
    WRITE(LOGNAM,*) 'totlyrnum       = ', totlyrnum
    WRITE(LOGNAM,*) 'revEgia         = ', revEgia
    WRITE(LOGNAM,*) 'lambda          = ', lambda
    WRITE(LOGNAM,*) 'lyrdph          = ', lyrdph
    WRITE(LOGNAM,*) 'psedD           = ', psedD
    WRITE(LOGNAM,*) 'pset            = ', pset
    WRITE(LOGNAM,*) 'pwatD            = ', pwatD
    WRITE(LOGNAM,*) 'visKin           = ', visKin
    WRITE(LOGNAM,*) 'vonKar           = ', vonKar
    WRITE(LOGNAM,*) 'crocdph          = ', TRIM(crocdph)
    WRITE(LOGNAM,*) 'sedD             = ', TRIM(sedD)
    WRITE(LOGNAM,*) 'csedfrc          = ', TRIM(csedfrc)

    write(LOGNAM,*) 'cslope    :', trim(cslope)
    write(LOGNAM,*) 'dsylunit  :', dsylunit
    write(LOGNAM,*) 'pyld      :', pyld
    write(LOGNAM,*) 'pyldc     :', pyldc
    write(LOGNAM,*) 'pyldpc    :', pyldpc
    write(LOGNAM,*) 'cinpmat_sed:', trim(cinpmat_sed)


    CLOSE(NSETFILE)

    WRITE(LOGNAM,*) 'CMF_SED_NMLIST: --- End Sediment Namelist ---'
  END SUBROUTINE CMF_SED_NMLIST

  SUBROUTINE CMF_SED_INIT
    !==========================================================
    !* PURPOSE: Initialize sediment module
    !==========================================================
    USE NETCDF
    use CMF_UTILS_MOD,           only: NCERROR
    IMPLICIT NONE

    ! Local variables
    INTEGER :: I, J,tmpnam ,ised ,jf,irec,ilyr      ! Loop counters
    INTEGER :: IERR                     ! Error status
    CHARACTER(LEN=256)              :: ctmp(20),CDATE
    REAL(KIND=JPRM)                 :: r2temp(NX,NY)
    REAL(KIND=JPRM)                 :: r3temp(NX,NY,nsed)
    INTEGER(KIND=JPIM)              ::  NCID,VARID



    WRITE(LOGNAM,*) 'CMF_SED_INIT: INITIALIZING SEDIMENT MODULE'
    !sediment_vars_init
    ! Allocate core sediment arrays (dimensions: NSEQMAX for spatial, nsed for sediment classes)
    WRITE(LOGNAM,*) 'CMF_SED_INIT: Allocating sediment arrays for NSEQMAX = ', NSEQMAX, ', nsed = ', nsed
    allocate(d2sedv(NSEQMAX,nsed,6))
    d2sedv(:,:,:) = 0._JPRB
    d2sedout => d2sedv(:,:,1)
    d2sedcon => d2sedv(:,:,2)
    d2sedinp => d2sedv(:,:,3)
    d2bedout => d2sedv(:,:,4)
    d2netflw => d2sedv(:,:,5)
    d2layer  => d2sedv(:,:,6)
    
    allocate(d2depv(NSEQMAX,totlyrnum,nsed))
    d2depv(:,:,:) = 0._JPRB
    d2seddep => d2depv

    allocate(d2sedv_avg(NSEQMAX,nsed,4))
    d2sedv_avg(:,:,:) = 0._JPRB
    d2sedout_avg => d2sedv_avg(:,:,1)
    d2sedinp_avg => d2sedv_avg(:,:,2)
    d2bedout_avg => d2sedv_avg(:,:,3)
    d2netflw_avg => d2sedv_avg(:,:,4)

    !sediment_map_init
    ! 3. Parse sedD and determine nsed
    ctmp(:) = '-999'
    call splitchar(sedD,ctmp)
    ised = 0
    allocate(sDiam(nsed))
    do i = 1, nsed
      if ( ctmp(i) /= '-999' ) then
        ised = ised + 1
        read(ctmp(i),*) sDiam(ised)
      endif
    enddo
    if ( ised /= nsed ) then
      write(LOGNAM,*) 'nsed and sedD do not match',ised,nsed
      stop
    endif
    write(LOGNAM,*) ised,' grain sizes: ',sDiam(:)

    ! 3. Calculate Settling Velocities
    ALLOCATE(SETVEL(nsed), STAT=IERR)
    IF(IERR /= 0) CALL CMF_ERROR_ABORT('CMF_SED_INIT: ALLOCATE SETVEL FAILED')
    SETVEL(:) = calc_settlingVelocity()
    print *, "csedfrc:", csedfrc
    WRITE(LOGNAM,*) 'CMF_SED_INIT: open csedfrc ',  TRIM(csedfrc)
    ! 5. Read sediment fraction file 
    allocate(d2sedfrc(NSEQMAX,nsed))

    if ( REGIONTHIS == 1 ) then
      CALL NCERROR (NF90_OPEN(csedfrc,NF90_NOWRITE,NCID),'opening '//TRIM(csedfrc) )
      CALL NCERROR (NF90_INQ_VARID(NCID, 'sedfrc',VARID),'getting id' )
      CALL NCERROR (NF90_GET_VAR(NCID,VARID,r3temp),'reading data' )
      CALL NCERROR (NF90_CLOSE(NCID),'closing '//TRIM(csedfrc))
    endif

    do ised = 1, nsed
      if ( REGIONTHIS == 1 )  then
        ! NetCDF: (NX, NY, ised) -> Fortran: (NX, NY)
        r2temp(:,:) = r3temp(:,:,ised)
        call mapR2vecD(r2temp,d2sedfrc(:,ised))
      endif 
    enddo

    ! adjust if any fractions are negative or if sum is not equal to 1
    if ( nsed == 1 ) then
      d2sedfrc(:,:) = 1.d0
    else
      !$omp parallel do private(iseq)
      do iseq = 1, NSEQALL
        if ( minval(d2sedfrc(iseq,:)) < 0.d0 .or. sum(d2sedfrc(iseq,:)) == 0.d0 ) then
          d2sedfrc(iseq,:) = 1.d0 / dble(nsed)
        else if ( sum(d2sedfrc(iseq,:)) /= 1.d0 ) then
          d2sedfrc(iseq,:) = d2sedfrc(iseq,:) / sum(d2sedfrc(iseq,:))
        endif
      enddo
      !$omp end parallel do
    endif    
    print *, "d2sedfrc:", d2sedfrc(1,:)
    !sediment_input_init 
    ! read slope
    call read_slope
    print *, "d2slope:"

!sediment_restart_init
!!!restart
  if ( sedrest_infile == "" ) then  ! set layer/bedload if no restart file
    !$omp parallel do
    do iseq = 1, NSEQALL
      d2layer(iseq,:) = lyrdph * D2RIVWTH(iseq,1) * D2RIVLEN(iseq,1) * d2sedfrc(iseq,:)
      do ilyr = 1, totlyrnum-1
        d2seddep(iseq,ilyr,:) = d2layer(iseq,:)
      enddo
      d2seddep(iseq,totlyrnum,:) = ( max(10.d0-lyrdph*totlyrnum,0.d0) ) * D2RIVWTH(iseq,1) * D2RIVLEN(iseq,1) * d2sedfrc(iseq,:)
    enddo
    !$omp end parallel do
    print *, "d2layer:", d2layer(1,:)
    print *, "d2seddep:", d2seddep(1,:,:)

  else
    if ( REGIONTHIS == 1 ) then
      tmpnam = INQUIRE_FID()
      open(tmpnam,file=sedrest_infile,form='unformatted',access='direct',recl=4*NX*NY)
    endif
    do irec = 1, 2
      do ised = 1, nsed
        if ( REGIONTHIS == 1 ) read(tmpnam,rec=(irec-1)*nsed+ised) r2temp
#ifdef UseMPI_CMF
        call MPI_Bcast(r2temp(1,1),NX*NY,mpi_real4,0,MPI_COMM_CAMA,ierr)
#endif
        select case(irec)
          case (1)
            call mapR2vecD(r2temp,d2layer(:,ised))
          case (2)
            call mapR2vecD(r2temp,d2sedcon(:,ised))
        end select
      enddo
    enddo

    do irec = 1, totlyrnum
      do ised = 1, nsed
        if ( REGIONTHIS == 1 ) read(tmpnam,rec=(irec+1)*nsed+ised) r2temp
#ifdef UseMPI_CMF
        call MPI_Bcast(r2temp(1,1),NX*NY,mpi_real4,0,MPI_COMM_CAMA,ierr)
#endif
        call mapR2vecD(r2temp,d2seddep(:,irec,ised))
      enddo
    enddo
    if ( REGIONTHIS == 1 ) close(tmpnam)
    write(LOGNAM,*) 'read restart sediment',maxval(d2seddep(:,totlyrnum,:))
  endif

  allocate(d2rivsto_pre(NSEQMAX), d2rivout_sed(NSEQMAX), d2rivvel_sed(NSEQMAX))
  sadd_riv = 0.d0
  sadd_out = 0.d0
  d2rivsto_pre(:) = P2RIVSTO(:,1)
  d2rivout_sed(:) = 0.d0
  d2rivvel_sed(:) = 0.d0
  print *, "d2rivsto_pre:", d2rivsto_pre(1)
  print *, "d2rivout_sed:", d2rivout_sed(1)
  print *, "d2rivvel_sed:", d2rivvel_sed(1)
  END SUBROUTINE CMF_SED_INIT


  subroutine CMF_SED_RESTART_WRITE
    use YOS_CMF_TIME,            only: KSTEP, NSTEPS, JDD, JHHMM, JHOUR, JMIN, JYYYYMMDD
    use YOS_CMF_INPUT,           only: CSUFBIN, RMIS
    use CMF_CTRL_RESTART_MOD,    only: CRESTDIR
    use CMF_UTILS_MOD,           only: vecD2mapR
#ifdef UseMPI_CMF
    use CMF_CTRL_MPI_MOD,        only: CMF_MPI_AllReduce_R2MAP
#endif
    
    implicit none
    save
    integer(kind=JPIM)              :: irec, irest, ised, tmpnam
    real(kind=JPRM)                 :: r3final(NX,NY,nsed), r2temp(NX,NY)
    character(len=256)              :: cdate, cfile
  
    !irest = 0
  
    !if ( ifrq_rst_sed>=0 .and. KSTEP==NSTEPS ) then  !! end of run
    !  irest = 1
    !endif
  
    !if ( ifrq_rst_sed>=1 .and. ifrq_rst_sed<=24 ) then  !! at selected hour
    !  if ( mod(JHOUR,ifrq_rst_sed)==0 .and. JMIN==0 ) then
    !    irest = 1
    !  endif
    !endif
  
    !if ( ifrq_rst_sed==30 ) then  !! at end of month
    !  if ( JDD==1 .and. JHOUR==0 .and. JMIN==0 ) then
    !    irest = 1
    !   endif
    !endif
  
    !if ( irest==1 ) then
      write(LOGNAM,*) ""
      write(LOGNAM,*) "!---------------------!"
      write(LOGNAM,*) 'cmf::sediment_restart_write: write time: ' , JYYYYMMDD, JHHMM
      WRITE(CDATE,'(I8.8,I2.2)') JYYYYMMDD,JHOUR

      cfile=trim(CRESTDIR)//TRIM(sedrest_outpre)//TRIM(cdate)//TRIM(CSUFBIN)
      write(LOGNAM,*) 'wrte_rest_bin: restart file:',cfile
  
      !*** write restart data (2D map)
      if ( REGIONTHIS == 1 ) then
        tmpnam = INQUIRE_FID()
        open(TMPNAM,file=cfile,form='unformatted',access='direct',recl=4*NX*NY*nsed)
      endif
      do irec = 1, 2
       r3final(:,:,:) = RMIS
       do ised = 1, nsed
         select case(irec)
           case (1)
             call vecD2mapR(d2layer(:,ised),r2temp)
           case (2)
             call vecD2mapR(d2sedcon(:,ised),r2temp)
         end select
#ifdef UseMPI_CMF
         call CMF_MPI_AllReduce_R2MAP(r2temp)
#endif
         r3final(:,:,ised) = r2temp(:,:)
       enddo
       if ( REGIONTHIS == 1 ) write(tmpnam,rec=irec) r3final
      enddo
  
      do irec = 1, totlyrnum
        r3final(:,:,:) = RMIS
        do ised = 1, nsed
          call vecD2mapR(d2seddep(:,irec,ised),r2temp)
#ifdef UseMPI_CMF
          call CMF_MPI_AllReduce_R2MAP(r2temp)
#endif
          r3final(:,:,ised) = r2temp
        enddo
        if ( REGIONTHIS == 1 ) write(tmpnam,rec=irec+2) r3final
      enddo
  
      if ( REGIONTHIS == 1 ) close(tmpnam)
  
    !endif
  end subroutine CMF_SED_RESTART_WRITE
  !####################################################################




  subroutine CMF_SED_CALC_FLW
    use PARKIND1,                only: JPIM, JPRB
    use YOS_CMF_INPUT,           only: PGRV, DT
    use YOS_CMF_MAP,             only: D2RIVLEN, D2RIVWTH, NSEQALL
    use YOS_CMF_PROG,            only: P2RIVSTO
    use YOS_CMF_DIAG,            only: D2RIVDPH

    implicit none
    !$ SAVE
    save
    integer(kind=JPIM)              :: ISEQ
    real(kind=JPRB)                 :: sedsto(NSEQALL,nsed)
    real(kind=JPRB)                 :: shearVel(NSEQALL)
    real(kind=JPRB)                 :: critShearVel(NSEQALL,nsed), dMean(NSEQALL), susVel(NSEQALL,nsed)
    real(kind=JPRB), parameter      :: IGNORE_DPH = 0.05d0
    !================================================

    !! calculate average of river water variables within sediment time step
    call sed_diag_average
  
    !! sediment concentration
    !$omp parallel do
    do iseq = 1, NSEQALL 
      sedsto(iseq,:) = d2sedcon(iseq,:) * max(d2rivsto_pre(iseq), 0.d0)
    enddo
    !$omp end parallel do
  
    call calc_params         !! calculate sediment transfer parameters
    call calc_advection
    call calc_entrainment
    call calc_exchange
  
    !$omp parallel do
    do iseq = 1, NSEQALL
      if ( P2RIVSTO(iseq,1) < D2RIVWTH(iseq,1)*D2RIVLEN(iseq,1)*IGNORE_DPH ) cycle
      d2sedcon(iseq,:) = sedsto(iseq,:) / P2RIVSTO(iseq,1)
    enddo
    !$omp end parallel do
  
    call sed_diag_reset
  
  contains
  !==========================================================
  !+ calc_params       !! calculate sediment transfer parameters
  !+ calc_advection
  !+ calc_entrainment
  !+ calc_exchange
  !==========================================================
    subroutine calc_params
      implicit none
      save
      integer(kind=JPIM)            ::  ised, iseq
      real(kind=JPRB)               ::  csVel0, sTmp, sTmp1(nsed)
      !=====================================================
      
      do iseq = 1, NSEQALL
        !-------------------------!
        ! critical shear velocity !
        !-------------------------!
        
        if ( sum(d2layer(iseq,:)) <= 0.d0 ) then
          critShearVel(iseq,:) = 1e20
        else if ( revEgia ) then       !! use Egiazoroff equation for mixed-size sedimnt
          dMean(iseq) = 0.d0
          do ised = 1, nsed
            dMean(iseq) = dMean(iseq) + sDiam(ised)*d2layer(iseq,ised)/sum(d2layer(iseq,:))
          enddo
          csVel0 = calc_criticalShearVelocity(dMean(iseq))
          do ised = 1, nsed
            if ( sDiam(ised) / dMean(iseq) >= 0.4d0 ) then
              critShearVel(iseq,ised) = sqrt( csVel0*sDiam(ised)/dMean(iseq) ) * &
                & ( log10(19.d0)/log10(19.d0*sDiam(ised)/dMean(iseq)) ) * 0.01d0
            else
              critShearVel(iseq,ised) = sqrt( 0.85*csVel0 ) * 0.01d0
            endif
          enddo      
        else
          do ised = 1, nsed
            critShearVel(iseq,ised) = sqrt( calc_criticalShearVelocity(sDiam(ised)) ) * 0.01d0
          enddo
        endif
      
        !------------------------------------------------------!
        ! shear velocity, suspend velocity, Karman coefficient !
        !------------------------------------------------------!
        if ( d2rivvel_sed(iseq) == 0.d0 .or. D2RIVDPH(iseq,1) < IGNORE_DPH ) then
          shearVel(iseq) = 0.d0
          susVel(iseq,:) = 0.d0
        else
          shearVel(iseq) = calc_shearVelocity(d2rivvel_sed(iseq), D2RIVDPH(iseq,1))
          susVel(iseq,:) = calc_suspendVelocity(critShearVel(iseq,:), shearVel(iseq), setVel(:))
        endif
      enddo
    end subroutine calc_params
    !=====================================================
  
    subroutine calc_advection
      use YOS_CMF_MAP,           only:  I1NEXT
    !  use yos_cmf_sed,           only:  d2rivout_sed, d2bedout, d2sedout, &
     !                                   d2bedout_avg, d2sedout_avg, psedD, pwatD
      implicit none
      real(kind=JPRB)               ::  bOut(NSEQALL,nsed), brate(NSEQALL,nsed)
      real(kind=JPRB)               ::  sOut(NSEQALL,nsed), srate(NSEQALL,nsed)
      integer(kind=JPIM)            ::  ised, iseq
      ! save for omop
      real(kind=JPRB), save         ::  plusVel, minusVel
      integer(kind=JPIM), save      ::  iseq0, iseq1
      !$omp threadprivate ( plusVel, minusVel, iseq0, iseq1 )
      !========
  
      bOut(:,:) = 0.d0
      sOut(:,:) = 0.d0
      !$omp parallel do
      do iseq = 1, NSEQALL
        
        if ( d2rivout_sed(iseq) >= 0.d0 ) then
          iseq0 = iseq
          iseq1 = I1NEXT(iseq)
        else
          iseq0 = I1NEXT(iseq)
          iseq1 = iseq
        endif
  
        if ( d2rivout_sed(iseq) == 0.d0 ) then
          d2sedout(iseq,:) = 0.d0
          d2bedout(iseq,:) = 0.d0
          cycle
        endif
  
        !-------------------!
        ! calc suspend flow !
        !-------------------!
        if ( iseq0 < 0 ) then
          d2sedout(iseq,:) = d2sedcon(iseq1,:) * d2rivout_sed(iseq)
        else
          d2sedout(iseq,:) = d2sedcon(iseq0,:) * d2rivout_sed(iseq)
          sOut(iseq0,:) = sOut(iseq0,:) + abs(d2sedout(iseq,:))*DT
        endif
  
        !--------------!
        ! calc bedflow !
        !--------------!
        if ( minval(critShearVel(iseq,:)) >= shearVel(iseq) .or. sum(d2layer(iseq,:)) == 0.d0 .or. iseq0 < 0  ) then
          d2bedout(iseq,:) = 0.d0
        else 
          do ised = 1, nsed
            if ( critShearVel(iseq,ised) >= shearVel(iseq) .or. d2layer(iseq,ised) == 0.d0 ) then
              d2bedout(iseq,ised) = 0.d0
              cycle
            endif
            plusVel = shearVel(iseq) + critShearVel(iseq,ised)
            minusVel = shearVel(iseq) - critShearVel(iseq,ised)
            d2bedout(iseq,ised) = 17.d0 * D2RIVWTH(iseq,1) * plusVel * minusVel * minusVel & 
             & / ((psedD-pwatD)/pwatD) / PGRV * d2layer(iseq,ised) / sum(d2layer(iseq,:)) 
            bOut(iseq0,ised) = bOut(iseq0,ised) + d2bedout(iseq,ised)*DT
          enddo
        endif
      enddo
      !$omp end parallel do
  
      !--------------------------------------------!
      ! adjust outflow if larget than sedsto/layer !
      !--------------------------------------------!
      brate(:,:) = 1.d0
      srate(:,:) = 1.d0
      !$omp parallel do
      do iseq = 1, NSEQALL
        if ( minval(sOut(iseq,:)) <= 1e-8 ) then
          do ised = 1, nsed
            if ( sOut(iseq,ised) > 1e-8 ) then
              srate(iseq,ised) = min ( sedsto(iseq,ised) / sOut(iseq,ised), 1.d0 )
            endif
          enddo
        else
          srate(iseq,:) = min ( sedsto(iseq,:) / sOut(iseq,:), 1.d0 )
        endif
        if ( minval(bOut(iseq,:)) <= 1e-8 ) then
          do ised = 1, nsed
            if ( bOut(iseq,ised) > 1e-8 ) then
              brate(iseq,ised) = min( d2layer(iseq,ised) / bOut(iseq,ised), 1.d0 )
            endif
          enddo
        else
          brate(iseq,:) = min( d2layer(iseq,:) / bOut(iseq,:), 1.d0 )
        endif
      enddo
      !$omp end parallel do
      
      do iseq = 1, NSEQALL
        if ( d2rivout_sed(iseq) >= 0.d0 ) then
          iseq0 = iseq
          iseq1 = I1NEXT(iseq)
        else
          iseq0 = I1NEXT(iseq)
          iseq1 = iseq
        endif
  
        if ( iseq0 > 0 ) then
          d2sedout(iseq,:) = d2sedout(iseq,:) * srate(iseq0,:)
          sedsto(iseq0,:) = max( sedsto(iseq0,:)-abs(d2sedout(iseq,:))*DT, 0.d0 )
          d2bedout(iseq,:) = d2bedout(iseq,:) * brate(iseq0,:)
          d2layer(iseq0,:) = max( d2layer(iseq0,:)-abs(d2bedout(iseq,:))*DT, 0.d0 )
        endif
        if ( iseq1 > 0 ) then
          sedsto(iseq1,:) = max( sedsto(iseq1,:)+abs(d2sedout(iseq,:))*DT, 0.d0 )
          d2layer(iseq1,:) = max( d2layer(iseq1,:)+abs(d2bedout(iseq,:))*DT, 0.d0 ) 
        endif
  
        d2bedout_avg(iseq,:) = d2bedout_avg(iseq,:) + d2bedout(iseq,:)*DT
        d2sedout_avg(iseq,:) = d2sedout_avg(iseq,:) + d2sedout(iseq,:)*DT
      enddo
  
    end subroutine calc_advection
    !=====================================================
  
    subroutine calc_entrainment
    !  use yos_cmf_sed,           only:  vonKar, d2netflw, d2netflw_avg, d2sedinp, d2seddep, totlyrnum
      
      implicit none
      real(kind=JPRB)               ::  dTmp(NSEQALL,nsed), D(NSEQALL,nsed), Es(NSEQALL,nsed), Zd(NSEQALL,nsed)
      integer(kind=JPIM)            ::  ilyr, ised, iseq
      real(kind=JPRB),save          ::  dTmp1, layerP
      !$omp threadprivate  ( dTmp1, layerP )
      !========
  
      !$omp parallel do
      do iseq = 1, NSEQALL
        if ( D2RIVDPH(iseq,1) < IGNORE_DPH ) then
          d2netflw(iseq,:) = 0.d0
          cycle
        endif
  
        !----------------------!
        ! calculate suspension !
        !----------------------!
        if ( sum(d2layer(iseq,:)) == 0.d0 .or. all(susVel(iseq,:)==0.d0) ) then
          Es(iseq,:) = 0.d0
        else
          Es(iseq,:) = susVel(iseq,:) * (1.d0-lambda) * D2RIVWTH(iseq,1) * D2RIVLEN(iseq,1) * d2layer(iseq,:) / sum(d2layer(iseq,:))
          Es(iseq,:) = max( Es(iseq,:), 0.d0 )
        endif
  
        !----------------------!
        ! calculate deposition !
        !----------------------!
        if ( shearVel(iseq) == 0.d0 .or. all(setVel(:)==0.d0) ) then
          D(iseq,:) = 0.d0
        else
          Zd(iseq,:) = 6.d0 * setVel(:) / vonKar / shearVel(iseq)
          D(iseq,:) = setVel(:) * D2RIVWTH(iseq,1) * D2RIVLEN(iseq,1) * d2sedcon(iseq,:) * Zd(iseq,:) / (1.d0-exp(-Zd(iseq,:))) 
          D(iseq,:) = max( D(iseq,:), 0.d0 )
        endif
        d2netflw(iseq,:) = Es(iseq,:) - D(iseq,:)
     
        !-------------------------------------------!
        ! if >0, suspension ; if <0, deposition     !
        ! adjust netflw if larger than sedsto/layer !
        !-------------------------------------------!
        do ised = 1, nsed
          if ( d2netflw(iseq,ised) == 0.d0 ) then
            cycle
          else if ( d2netflw(iseq,ised) > 0.d0 ) then
            dTmp1 = d2netflw(iseq,ised)*DT/(1.d0-lambda)
            if ( dTmp1 < d2layer(iseq,ised) ) then
              d2layer(iseq,ised) = d2layer(iseq,ised) - dTmp1
            else
              d2netflw(iseq,ised) = d2layer(iseq,ised) * (1.d0-lambda) / DT
              d2layer(iseq,ised) = 0.d0
            endif
            sedsto(iseq,ised) = sedsto(iseq,ised) + d2netflw(iseq,ised) * DT
          else
            if ( abs(d2netflw(iseq,ised))*DT < sedsto(iseq,ised) ) then
              sedsto(iseq,ised) = max (sedsto(iseq,ised) - abs(d2netflw(iseq,ised))*DT, 0.d0 )
            else
              d2netflw(iseq,ised) = - sedsto(iseq,ised) / DT
              sedsto(iseq,ised) = 0.d0
            endif
            d2layer(iseq,ised) = d2layer(iseq,ised) + abs(d2netflw(iseq,ised))*DT/(1.d0-lambda)
          endif
        enddo
  
        sedsto(iseq,:) = sedsto(iseq,:) + d2sedinp(iseq,:)*DT    
        if ( sum(sedsto(iseq,:)) > P2RIVSTO(iseq,1) * 0.01d0 ) then
          dTmp(iseq,:) = ( sum(sedsto(iseq,:)) - P2RIVSTO(iseq,1)*0.01d0 ) * sedsto(iseq,:)/sum(sedsto(iseq,:))
          d2netflw(iseq,:) = d2netflw(iseq,:) - dTmp(iseq,:)/DT
          sedsto(iseq,:) = sedsto(iseq,:) - dTmp(iseq,:)
          d2layer(iseq,:) = d2layer(iseq,:) + dTmp(iseq,:)/(1.d0-lambda)
        endif
          
        d2netflw_avg(iseq,:) = d2netflw_avg(iseq,:) + d2netflw(iseq,:)*DT
      enddo
      !$omp end parallel do
    end subroutine calc_entrainment
    !=====================================================
  
    subroutine calc_exchange
      ! redistribute into vertical bed layers
     ! use yos_cmf_sed,           only:  d2seddep, lyrdph, totlyrnum
      
      implicit none
      integer(kind=JPIM)            ::  ilyr, ised, iseq, jlyr, slyr
      real(kind=JPRB)               ::  diff, lyrvol, layerP(nsed), seddepP(totlyrnum+1,nsed), tmp(nsed)
  
      do iseq = 1, NSEQALL
        lyrvol = lyrdph * D2RIVWTH(iseq,1) * D2RIVLEN(iseq,1)
  
        if ( minval(d2layer(iseq,:)) < 0.d0 ) d2layer(iseq,:) = max( d2layer(iseq,:), 0.d0 )
        if ( minval(d2seddep(iseq,:,:)) < 0.d0 ) d2seddep(iseq,:,:) = max( d2seddep(iseq,:,:), 0.d0 )
  
        !---------------------------------------!
        ! if bed storage less than layer volume !
        !---------------------------------------!
        if ( sum(d2layer(iseq,:)) + sum(d2seddep(iseq,:,:)) <= lyrvol ) then
          d2layer(iseq,:) = d2layer(iseq,:) + sum(d2seddep(iseq,:,:),dim=1)
          d2seddep(iseq,:,:) = 0.d0
          cycle
        endif
  
        !------------------------------------!
        ! distribute into top exchange layer !
        !------------------------------------!
        layerP(:) = d2layer(iseq,:)
        if ( sum(layerP(:)) >= lyrvol ) then
          d2layer(iseq,:) = layerP(:) * min( lyrvol/sum(layerP(:)), 1.d0 )
          layerP(:) = max( layerP(:) - d2layer(iseq,:), 0.d0 )
          slyr = 0
        else if ( sum(d2seddep(iseq,:,:)) > 0.d0 ) then
          layerP(:) = 0.d0
          do ilyr = 1, totlyrnum
            diff = lyrvol - sum(d2layer(iseq,:))
            if ( diff <= 0.d0 ) exit
            if ( sum(d2seddep(iseq,ilyr,:)) <= diff ) then
              d2layer(iseq,:) = d2layer(iseq,:) + d2seddep(iseq,ilyr,:)
              d2seddep(iseq,ilyr,:) = 0.d0
              slyr = ilyr + 1
            else
              tmp(:) = diff * d2seddep(iseq,ilyr,:) / sum(d2seddep(iseq,ilyr,:))
              d2layer(iseq,:) = d2layer(iseq,:) + tmp(:)
              d2seddep(iseq,ilyr,:) = max( d2seddep(iseq,ilyr,:) - tmp(:), 0.d0 )
              slyr = ilyr
              exit
            endif
          enddo
        else
          d2seddep(iseq,:,:) = 0.d0
          cycle
        endif
        if ( sum(d2seddep(iseq,:,:)) == 0.d0 ) cycle
  
        !-----------------------------------!
        ! distribute remaining bedload into !
        ! vertical deposition layers        !
        !-----------------------------------!
        seddepP(1,:) = layerP(:)
        seddepP(2:,:) = d2seddep(iseq,:,:)
        d2seddep(iseq,:,:) = 0.d0
        do ilyr = 1, totlyrnum - 1
          if ( sum(d2seddep(iseq,ilyr,:)) == lyrvol ) cycle
          do jlyr = slyr+1, totlyrnum + 1
            diff = lyrvol - sum(d2seddep(iseq,ilyr,:))
            if ( diff <= 0.d0 ) exit
            if ( sum(seddepP(jlyr,:)) <= diff ) then
              d2seddep(iseq,ilyr,:) = d2seddep(iseq,ilyr,:) + seddepP(jlyr,:)
              seddepP(jlyr,:) = 0.d0
            else
              tmp(:) = diff * seddepP(jlyr,:) / sum(seddepP(jlyr,:))
              d2seddep(iseq,ilyr,:) = d2seddep(iseq,ilyr,:) + tmp(:)
              seddepP(jlyr,:) = max(seddepP(jlyr,:) - tmp(:), 0.d0)
              exit
            endif
          enddo
        enddo
        
        if ( sum(seddepP) > 0.d0 ) then
          d2seddep(iseq,totlyrnum,:) = sum(seddepP, dim=1)
        endif
      enddo
  
    end subroutine calc_exchange
  end subroutine CMF_SED_CALC_FLW

!==========================================================
  subroutine calc_sedyld(pbuffin)
    use PARKIND1,                only: JPIM, JPRB
    use YOS_CMF_INPUT,           only: DTIN
    
    implicit none
    save
    real(kind=JPRB), intent(in)     :: pbuffin(:)
    real(kind=JPRB)                 :: sbuff(NSEQMAX)

    !================================================
  
    
    call prcp_convert_sed(pbuffin, sbuff) ! convert precipitation to sediment yield based on Sunada&Hasegawa(1993)
   
    !$omp parallel do 
    do iseq = 1, NSEQALL
      d2sedinp(iseq,:) = sbuff(iseq) * d2sedfrc(iseq,:)  ! distribute sediment yield to proportionate to sediment grain fraction
      d2sedinp_avg(iseq,:) = d2sedinp_avg(iseq,:) + d2sedinp(iseq,:) * DTIN
    enddo
    !$omp end parallel do
  
  contains
  !=============================
  !+ prcp_convert_sed
  !=============================
    subroutine prcp_convert_sed(pbuffin,pbuffout)
      use YOS_CMF_DIAG,          only: D2FLDFRC
  
      implicit none
      save
      real(kind=JPRB), intent(in)   :: pbuffin(:)     !! kg/m2/s
      real(kind=JPRB), intent(out)  :: pbuffout(:)  !! m3/s
      integer(kind=JPIM)            :: i
  
      !$omp parallel do
      do iseq = 1, NSEQALL
        pbuffout(iseq) = 0.d0
        if ( pbuffin(iseq) * 86400.d0 <= 10.d0 ) cycle
  
        do i = 1, NLFP
          if ( D2FLDFRC(iseq,1) * NLFP > dble(i) ) cycle  ! no erosion if submerged
          pbuffout(iseq) = pbuffout(iseq) + pyld * (pbuffin(iseq)*3600.d0)**pyldpc * d2slope(iseq,i)**pyldc / 3600.d0 & 
            & * D2GRAREA(iseq,1) * min(dble(i)/dble(NLFP)-D2FLDFRC(iseq,1), 1.d0/dble(NLFP)) * dsylunit
        enddo
      enddo
      !$omp end parallel do
    end subroutine prcp_convert_sed
      
  end subroutine calc_sedyld
  !==========================================================
  !+
  !==========================================================
  subroutine sedinp_interp(pbuffin,pbuffout)
  ! interporlate sediment forcing data using "input matrix"
    use YOS_CMF_INPUT,           only: NXIN, NYIN, INPN, RMIS
    
    implicit none
    real(kind=JPRM),intent(in)      :: pbuffin(:,:)     !! default for prcp[kg/m2/s]
    real(kind=JPRB),intent(out)     :: pbuffout(:)    !! kg/m2/s
    ! save for omp
    integer(kind=jpim),save  ::   ixin, iyin, inpi  !! for output
    !$omp threadprivate    (ixin, iyin)
    !============================
    !$omp parallel do
    do iseq=1, NSEQALL
      pbuffout(iseq)=0._JPRB
      do inpi=1, INPN
        ixin=INPX(iseq,inpi)
        iyin=INPY(iseq,inpi)
        if( ixin>0 )then
          if( ixin > NXIN .or. iyin > NYIN ) then
            write(LOGNAM,*)  "error"
            write(LOGNAM,*)  'xxx',iseq,inpi,ixin,iyin
            cycle
          endif
          if( pbuffin(ixin,iyin).ne.RMIS )then
            pbuffout(iseq) = pbuffout(iseq) + pbuffin(ixin,iyin) * INPA(iseq,inpi) / D2GRAREA(iseq,1)
          endif
        endif
      end do
      pbuffout(iseq)=max(pbuffout(iseq), 0._JPRB)
    end do
    !$omp end parallel do
  end subroutine sedinp_interp



  subroutine splitchar(allvars,vnames)
    ! same function as splitting characters in CaMa
    use PARKIND1,                only: JPIM
    implicit none
    save
    character(len=256), intent(in)  :: allvars
    character(len=256), intent(out) :: vnames(:)
    integer(kind=JPIM)              :: nvarsout, j0, j
    character(len=256)              :: ctmp
  
    nvarsout = 0
    j0 = 1
    do j = 1, len(trim(allvars))
      if ( (j>j0) .and. (allvars(j:j).eq.',') ) then
        ctmp = trim(adjustl(allvars(j0:j-1)))
        if ( len(ctmp) > 0 ) then
          nvarsout = nvarsout + 1
          vnames(nvarsout) = ctmp
        endif
        j0 = j + 1
      endif
    enddo
  
    ! last one
    if ( j0 < len(trim(allvars)) ) then
      j = len(trim(allvars))
      ctmp = trim(adjustl(allvars(j0:j)))
      if ( len(ctmp) > 0 ) then
        nvarsout = nvarsout + 1
        vnames(nvarsout) = ctmp
      endif
    endif
  end subroutine splitchar

subroutine sed_diag_average
  implicit none
  !! calculate average of river water variables within sediment time step
  d2rivout_sed(:) = d2rivout_sed(:) /dble(sadd_riv)
  d2rivvel_sed(:) = d2rivvel_sed(:) /dble(sadd_riv)
end subroutine sed_diag_average

!==========================================================
!+
!==========================================================
subroutine sed_diag_reset
  use PARKIND1,                only: JPRB
  use YOS_CMF_PROG,            only: P2RIVSTO

  implicit none

  sadd_riv = 0
  d2rivout_sed(:) = 0._JPRB
  d2rivvel_sed(:) = 0._JPRB
  d2rivsto_pre(:) = P2RIVSTO(:,1)

  sadd_out = sadd_out + DT
end subroutine sed_diag_reset


function calc_settlingVelocity() result(setVel)
  use YOS_CMF_INPUT,           only:  PGRV
  implicit none
  save
  real(kind=JPRB)                 ::  setVel(nsed) ! settling velocity [m/s]
  real(kind=JPRB)                 ::  sTmp(nsed)

  sTmp(:) = 6.d0 * visKin / sDiam(:)
  setVel(:) = pset * ( sqrt( 2.d0/3.d0*(psedD-pwatD)/pwatD*PGRV*sDiam(:) &
                             + sTmp(:)*sTmp(:) ) - sTmp(:) )
end function calc_settlingVelocity

function calc_criticalShearVelocity(diam) result(csVel)
  implicit none
  save
  real(kind=JPRB)                 ::  csVel ! critical shear velocity[(cm/s)^2]
  real(kind=JPRB), intent(in)     ::  diam ![m]
  real(kind=JPRB)                 ::  cA, cB
  !========
  cB = 1.d0
  if ( diam >= 0.00303d0 ) then
    cA = 80.9d0
  else if ( diam >= 0.00118d0 ) then
    cA = 134.6d0
    cB = 31.d0 / 32.d0
  else if ( diam >= 0.000565d0 ) then
    cA = 55.d0
  else if ( diam >= 0.000065d0 ) then
    cA = 8.41d0
    cB = 11.d0 / 32.d0
  else
    cA = 226.d0
  endif
  
  csVel = cA * ( diam*100.d0 ) ** cB
  return
end function calc_criticalShearVelocity

function calc_shearVelocity(rivvel,rivdph) result(sVel)
  use YOS_CMF_INPUT,           only:  PGRV, PMANRIV
  implicit none
  save
  real(kind=JPRB)                 ::  sVel ! shear velocity[m/s]
  real(kind=JPRB), intent(in)     ::  rivvel, rivdph 
  !========

  sVel = sqrt ( PGRV * PMANRIV**2.d0 * rivvel**2.d0 * rivdph**(-1.d0/3.d0) )  !bug fix 2022/11/22
  return
end function calc_shearVelocity

function calc_suspendVelocity(csVel,sVel,setVel) result(susVel) ! Uchida and Fukuoka (2019) Eq.44
  implicit none
  save
  real(kind=JPRB)                 ::  susVel(nsed)   ! suspend velocity [m/s]
!!  real(kind=JPRB), intent(in)     ::  csVel(nsed), sVel, setVel(nsed)    ! crit shear, shear velocity, setting velocity [m/s]
  real(kind=JPRB), intent(in)     ::  csVel(:), sVel, setVel(:)    ! crit shear, shear velocity, settling velocity [m/s]
  integer(kind=JPIM)              ::  ised
  real(kind=JPRB)                 ::  alpha, a, cB, sTmp
  !========
  
  !--------------!
  ! set constant !
  !--------------!
  alpha = vonKar / 6.d0
  a = 0.08d0
  cB = 1.d0 - lambda
  
  !-----------------------!
  ! calc suspend velocity !
  !-----------------------!
  susVel(:) = 0.d0
  do ised = 1, nsed
    if ( csVel(ised) > sVel ) cycle
    sTmp = setVel(ised) / alpha / sVel
    susVel(ised) = max( setVel(ised) * cB / (1.d0+sTmp) * (1.d0-a*sTmp) / (1.d0+(1.d0-a)*sTmp), 0.d0 )
  enddo
end function calc_suspendVelocity

subroutine read_slope
  use YOS_CMF_INPUT,         only: NX,NY, NLFP
  use YOS_CMF_MAP,           only: REGIONTHIS
  use CMF_UTILS_MOD,           only: NCERROR
  use NETCDF

  implicit none
  integer                       ::  tmpnam, i
  real(kind=jprm)               :: r2temp(nx,ny)
  real(kind=jprm)               :: r3temp(NX,NY,NLFP)
  INTEGER(KIND=JPIM)              ::  NCID,VARID

  allocate(d2slope(NSEQMAX,NLFP))
  CALL NCERROR (NF90_OPEN(cslope,NF90_NOWRITE,NCID),'opening '//TRIM(cslope) )
  CALL NCERROR (NF90_INQ_VARID(NCID, 'slope',VARID),'getting id' )
  CALL NCERROR (NF90_GET_VAR(NCID,VARID,r3temp),'reading data' )
  CALL NCERROR (NF90_CLOSE(NCID),'closing '//TRIM(cslope))
  do i = 1, NLFP
    r2temp(:,:) = r3temp(:,:,i)
    call mapR2vecD(r2temp,d2slope(:,i))
  enddo

end subroutine read_slope


!####################################################################
SUBROUTINE CMF_SED_FORCING_PUT(PBUFF)
  ! Simplified sediment forcing - use pre-allocated sediment buffer
  ! -- called from "Main Program / Coupler" or CMF_DRV_ADVANCE
  IMPLICIT NONE
  ! Declaration of arguments
  real(KIND=JPRB), intent(in)     :: PBUFF(:,:)  ! 2D sediment buffer [NX, NY]
  REAL(KIND=JPRB)                 :: D2TEMP(NSEQMAX,1)  ! Changed to 2D array

    CALL SEDINP_INTERP(REAL(PBUFF, KIND=JPRM),D2TEMP(:,1))
    CALL CALC_SEDYLD(D2TEMP(:,1))  ! Pass 1D slice

END SUBROUTINE CMF_SED_FORCING_PUT

SUBROUTINE CMF_SED_DIAG_AVEMAX_ADPSTP
  IMPLICIT NONE
  !====================
   !calculate average rivout and rivvel for sediment timestep
  sadd_riv = sadd_riv + DT
  !$OMP PARALLEL DO SIMD
  DO ISEQ=1, NSEQMAX
    d2rivout_sed(ISEQ) = d2rivout_sed(ISEQ)+D2RIVOUT(ISEQ,1)*DT
    d2rivvel_sed(ISEQ) = d2rivvel_sed(ISEQ)+D2RIVVEL(ISEQ,1)*DT
  END DO
  !$OMP END PARALLEL DO SIMD
  END SUBROUTINE CMF_SED_DIAG_AVEMAX_ADPSTP



  !####################################################################
!SUBROUTINE CMF_DIAG_GETAVE_OUTPUT
!  USE YOS_CMF_TIME,       ONLY: JYYYYMMDD, JHHMM
!  IMPLICIT NONE
!  INTEGER(KIND=JPIM),SAVE  ::  ISEQ, IPTH
!  INTEGER(KIND=JPIM),SAVE  ::  ISED
  !================================================
  
  ! Average sediment variables
!   !$OMP PARALLEL DO SIMD
!    DO ISEQ=1, NSEQMAX
!      DO ISED=1, nsed
!        d2sedout_avg(ISEQ,ISED) = d2sedout_avg(ISEQ,ISED) / REAL(NADD_out,KIND=JPRB)
!        d2sedinp_avg(ISEQ,ISED) = d2sedinp_avg(ISEQ,ISED) / REAL(NADD_out,KIND=JPRB)
!        d2bedout_avg(ISEQ,ISED) = d2bedout_avg(ISEQ,ISED) / REAL(NADD_out,KIND=JPRB)
!        d2netflw_avg(ISEQ,ISED) = d2netflw_avg(ISEQ,ISED) / REAL(NADD_out,KIND=JPRB)
!      END DO
!      ! Note: d2sedcon and d2layer are instantaneous values, not accumulated
!      ! They don't need averaging like the flux variables
!      ! d2seddep is also not accumulated in the current implementation
!   END DO
!  !$OMP END PARALLEL DO SIMD
  
!  END SUBROUTINE CMF_DIAG_GETAVE_OUTPUT

END MODULE CMF_CTRL_SED_MOD