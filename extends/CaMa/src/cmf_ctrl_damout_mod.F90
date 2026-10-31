MODULE CMF_CTRL_DAMOUT_MOD
!==========================================================
!* PURPOSE: CaMa-Flood reservoir operation scheme (under development)
!
! (C) R. Hanazaki & D.Yamazaki (U-Tokyo)  Feb 2020
!
!* CONTAINS:
! -- CMF_DEM_NMLIST  : Read setting from namelist
! -- CMF_DAM_INIT    : Initialize dam data
! -- CMF_CALC_DAMOUT : Calculate inflow and outflow at dam
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
USE PARKIND1,                ONLY: JPIM, JPRB, JPRM, JPRD
USE CMF_UTILS_MOD,           ONLY: INQUIRE_FID
USE YOS_CMF_INPUT,           ONLY: LOGNAM, IMIS, LDAMOUT, LPTHOUT, LRESTART, LDAMIRR
USE YOS_CMF_INPUT,           ONLY: NX, NY, DT
USE YOS_CMF_MAP,             ONLY: I2VECTOR, I1NEXT, NSEQMAX, NSEQALL, NSEQRIV
USE YOS_CMF_MAP,             ONLY: NPTHOUT,  NPTHLEV, PTH_UPST, PTH_DOWN, PTH_ELV, I2MASK!! bifurcation pass
USE YOS_CMF_PROG,            ONLY: D2RIVOUT, D2FLDOUT, P2RIVSTO, P2FLDSTO, P2DAMSTO, P2DAMINF, D2RUNOFF
USE YOS_CMF_PROG,            ONLY: dirrig_cama, release_cama_riv, release_cama_dam, dirrig_cama_unmt
USE YOS_CMF_DIAG,            ONLY: D2RIVINF, D2FLDINF
USE YOS_CMF_TIME,            ONLY: IYYYYMMDD,ISYYYY
!============================
IMPLICIT NONE
SAVE
!*** NAMELIST/NDAM/
CHARACTER(LEN=256)              :: CDAMFILE    !! dam paramter file
LOGICAL                         :: LDAMTXT     !! true: dam inflow-outflw txt output
LOGICAL                         :: LDAMH22     !! true: Use Hanazaki 2022 scheme
LOGICAL                         :: LDAMYBY     !! true: Use Year-By-Year dam activation
LOGICAL                         :: LiVnorm     !! true: initialize dam storage with Normal Volume
CHARACTER(LEN=3)                :: LDAMOPT     !! dam scheme
NAMELIST/NDAMOUT/   CDAMFILE, LDAMTXT, LDAMH22, LDAMYBY, LiVnorm, LDAMOPT

!*** dam map
INTEGER(KIND=JPIM),ALLOCATABLE  :: DamSeq(:)   !! coresponding ISEQ of each dam
INTEGER(KIND=JPIM),ALLOCATABLE  :: I1DAM(:)    !! dam map: 1=dam, 10=upstream of dam, 11: dam grid & downstream is also dam, 0=other

!*** dam parameters
INTEGER(KIND=JPIM)              :: IDAM, NDAM  !! number of dams
INTEGER(KIND=JPIM)              :: NDAMX       !! exclude dams

INTEGER(KIND=JPIM),ALLOCATABLE  :: DamID(:) !! Dam ID
CHARACTER(LEN=256),ALLOCATABLE  :: DamName(:)  !! 
INTEGER(KIND=JPIM),ALLOCATABLE  :: DamIX(:), DamIY(:)  !! IX,IY of dam grid
REAL(KIND=JPRB),ALLOCATABLE     :: DamLon(:), DamLat(:)  !! longitude, latitude of dam body
REAL(KIND=JPRB),ALLOCATABLE     :: upreal(:)   !! observed drainage area of reservoir
REAL(KIND=JPRB),ALLOCATABLE     :: R_VolUpa(:) !! ratio: flood storage capacity / drainage area
REAL(KIND=JPRB),ALLOCATABLE     :: Qf(:), Qn(:), Qe(:)!! Qf: flood discharge, Qn: normal discharge
CHARACTER(LEN=256),ALLOCATABLE  :: MainUse(:)             !! main use of dam
INTEGER(KIND=JPIM),ALLOCATABLE  :: DamYear(:)  !! Dam activation year
INTEGER(KIND=JPIM),ALLOCATABLE  :: DamStat(:)  !! Dam Status, 2=old, 1=new, -1=not_yet, IMIS=out_of_domain

REAL(KIND=JPRB),ALLOCATABLE     :: TotVol(:)            !! total storage capacity of reservoir (mcm)
REAL(KIND=JPRB),ALLOCATABLE     :: EmeVol(:)   !! storage volume to start emergency operation
REAL(KIND=JPRB),ALLOCATABLE     :: FldVol(:)   !! flood control volume: exclusive for flood control
REAL(KIND=JPRB),ALLOCATABLE     :: ConVol(:)   !! conservative volume: mainly for water supply
REAL(KIND=JPRB),ALLOCATABLE     :: NorVol(:)   !! normal storage volume: impoundment

! internal dam param for stability
REAL(KIND=JPRB),ALLOCATABLE     :: AdjVol(:)   !! Dam storage for stabilization 
REAL(KIND=JPRB),ALLOCATABLE     :: Qa(:)       !! Dam outflow for stabilization

!*** water use data 
REAL(KIND=JPRB),ALLOCATABLE     :: WUSE_DD(:,:)         !! daily water demand (m3/s)
REAL(KIND=JPRB),ALLOCATABLE     :: WUSE_AD(:)           !! average daily demand (m3/s)

!!!!----- for IRR
!*** dam parameters for different schemes
!*** H22
REAL(KIND=JPRB),ALLOCATABLE     :: H22_EmeVol(:)        !! storage volume to start emergency operation (m3)
REAL(KIND=JPRB),ALLOCATABLE     :: H22_k(:)             !! release coefficient (-)
!*** H06
REAL(KIND=JPRB),ALLOCATABLE     :: H06_c(:)             !! the ratio between capacity and mean annual inflow (-)
REAL(KIND=JPRB),ALLOCATABLE     :: H06_DPI(:)           !! the ratio between annual mean demand and annual mean inflow (-) 
!*** V13
INTEGER(KIND=JPIM),ALLOCATABLE  :: StFC_Mth(:)          !! the start of the flood control period
INTEGER(KIND=JPIM),ALLOCATABLE  :: NdFC_Mth(:)          !! the end of the flood control period
INTEGER(KIND=JPIM),ALLOCATABLE  :: StOP_Mth(:)          !! the start of the operational year

!*** water use grids
INTEGER(KIND=JPIM),ALLOCATABLE     :: serial_num(:),grid_num(:)
INTEGER(KIND=JPIM),ALLOCATABLE     :: grids_x(:,:),grids_y(:,:)
CHARACTER(LEN=16), ALLOCATABLE     :: dam_id(:)
REAL(KIND=JPRB),   ALLOCATABLE     :: area_1(:),area_2(:)
REAL(KIND=JPRB),   ALLOCATABLE     :: grids_share(:,:)
INTEGER(KIND=JPIM)                 :: max_gridnum

!*** dam water use data
REAL(KIND=JPRB), ALLOCATABLE       :: dam_tot_demand(:), dam_tot_demand_save(:), dam_tot_demand_unmt(:)
REAL(KIND=JPRB), ALLOCATABLE       :: dam_grid_demand(:,:), dam_grid_demand_unmt(:,:)
REAL(KIND=JPRB), ALLOCATABLE       :: save_damwithdraw(:)

!*** river water use
INTEGER(KIND=JPIM)                 :: NRIV           !! number of irrig grids
INTEGER(KIND=JPIM)                 :: NRIVX   
INTEGER(KIND=JPIM),ALLOCATABLE     :: IX_RIV(:), IY_RIV(:)
REAL(KIND=JPRB), ALLOCATABLE       :: riv_tot_demand(:), riv_tot_demand_unmt(:)
REAL(KIND=JPRB), ALLOCATABLE       :: save_rivwithdraw(:)

CONTAINS
!####################################################################
!* CONTAINS:
! -- CMF_DAMOUT_NMLIST_ORG  : Read setting from namelist
! -- CMF_DAMOUT_INIT_ORG    : Initialize dam data
! -- CMF_DAMOUT_CALC_ORG    : Calculate inflow and outflow at dam
! -- CMF_DAMOUT_WATBAL_ORG  : Calculate water balance at dam
! -- CMF_DAMOUT_WRTE_ORG   : Write dam-related variables in text file
!####################################################################
SUBROUTINE CMF_DAMOUT_NMLIST_ORG
! reed setting from namelist
! -- Called from CMF_DRV_NMLIST
USE YOS_CMF_INPUT,      ONLY: CSETFILE,NSETFILE,LDAMOUT
IMPLICIT NONE
!================================================
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!---------------------!"

!*** 1. open namelist
NSETFILE=INQUIRE_FID()
OPEN(NSETFILE,FILE=CSETFILE,STATUS="OLD")
WRITE(LOGNAM,*) "CMF::DAMOUT_NMLIST: namelist OPEN in unit: ", TRIM(CSETFILE), NSETFILE 

!*** 2. default value
CDAMFILE="./dam_params.csv"
LDAMTXT=.TRUE.
LDAMH22=.FALSE.
LDAMYBY=.FALSE.
LiVnorm=.FALSE.

!*** 3. read namelist
REWIND(NSETFILE)
READ(NSETFILE,NML=NDAMOUT)

IF( LDAMOUT )THEN
  WRITE(LOGNAM,*)   "=== NAMELIST, NDAMOUT ==="
  WRITE(LOGNAM,*)   "CDAMFILE: " , CDAMFILE
  WRITE(LOGNAM,*)   "LDAMTXT:  " , LDAMTXT
  WRITE(LOGNAM,*)   "LDAMH22:  " , LDAMH22
  WRITE(LOGNAM,*)   "LDAMYBY:  " , LDAMYBY
  WRITE(LOGNAM,*)   "LiVnorm: " , LiVnorm
ENDIF

CLOSE(NSETFILE)

WRITE(LOGNAM,*) "CMF::DAMOUT_NMLIST: end" 

END SUBROUTINE CMF_DAMOUT_NMLIST_ORG
!####################################################################





!####################################################################
SUBROUTINE CMF_DAMOUT_INIT_ORG
! reed setting from CDAMFILE
IMPLICIT NONE
INTEGER(KIND=JPIM)         :: NDAMFILE
INTEGER(KIND=JPIM)         :: ISEQ, JSEQ
INTEGER(KIND=JPIM)         :: IX, IY
REAL(KIND=JPRB)            :: FldVol_mcm, ConVol_mcm, TotVol_mcm !! from file in Million Cubic Metter
REAL(KIND=JPRB)            :: Qsto, Vyr

INTEGER(KIND=JPIM)         :: IPTH, ILEV, ISEQP, JSEQP
!####################################################################
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!---------------------!"
WRITE(LOGNAM,*) "CMF::DAMOUT_INIT: initialize dam", trim(CDAMFILE) 

!==========
NDAMFILE=INQUIRE_FID()
OPEN(NDAMFILE,FILE=CDAMFILE,STATUS="OLD")
READ(NDAMFILE,*) NDAM
READ(NDAMFILE,*)        !! skip header

WRITE(LOGNAM,*) "CMF::DAMOUT_INIT: number of dams", NDAM

!! === ALLOCATE ===
!! from CDAMFILE
ALLOCATE(DamID(NDAM),DamName(NDAM))
ALLOCATE(DamIX(NDAM),DamIY(NDAM),DamLon(NDAM),DamLat(NDAM))
ALLOCATE(upreal(NDAM))
ALLOCATE(Qf(NDAM),Qn(NDAM))
ALLOCATE(DamYear(NDAM),DamStat(NDAM))

!! calculate from CDAMFILE
ALLOCATE(DamSeq(NDAM))
ALLOCATE(FldVol(NDAM),ConVol(NDAM),EmeVol(NDAM),NorVol(NDAM))

!! for outflw stability
ALLOCATE(AdjVol(NDAM),Qa(NDAM))

!! H22scheme parameter (FldVol/Upreal)
ALLOCATE(R_VolUpa(NDAM))

!! dam map, dam variable
ALLOCATE(I1DAM(NSEQMAX))
!! =================
DamSeq(:) =IMIS
DamStat(:)=IMIS
I1DAM(:)=0
NDAMX=0
!! read dam parameters
DO IDAM = 1, NDAM
  IF( LDAMYBY) THEN
    READ(NDAMFILE,*) DamID(IDAM), DamName(IDAM), DamLat(IDAM), DamLon(IDAM), upreal(IDAM), &
     DamIX(IDAM), DamIY(IDAM), FldVol_mcm, ConVol_mcm, TotVol_mcm, Qn(IDAM), Qf(IDAM), DamYear(IDAM)
  ELSE
    READ(NDAMFILE,*) DamID(IDAM), DamName(IDAM), DamLat(IDAM), DamLon(IDAM), upreal(IDAM), &
     DamIX(IDAM), DamIY(IDAM), FldVol_mcm, ConVol_mcm, TotVol_mcm, Qn(IDAM), Qf(IDAM)
  ENDIF

  !! storage parameter --- from Million Cubic Meter to m3
  FldVol(IDAM) = FldVol_mcm * 1.E6_JPRB  ! Flood control storage capacity: exclusive for flood control
  ConVol(IDAM) = ConVol_mcm * 1.E6_JPRB

  EmeVol(IDAM) = ConVol(IDAM) + FldVol(IDAM) * 0.95_JPRB  ! storage to start emergency operation

  IX=DamIX(IDAM)
  IY=DamIY(IDAM)
  IF (IX<=0 .or. IX > NX .or. IY<=0 .or. IY > NY ) cycle

  ISEQ=I2VECTOR(IX,IY)
  IF( ISEQ<=0 ) cycle
  IF( I1NEXT(ISEQ)==-9999 ) cycle
  NDAMX=NDAMX+1

  DamSeq(IDAM) =ISEQ
  DamStat(IDAM)=2

  I1DAM(ISEQ)=1
  I2MASK(ISEQ,1)=2   !! reservoir grid. skipped for adaptive time step

  IF( LDAMH22 )THEN    !! Hanazaki 2022 scheme 
    NorVol(IDAM)   = ConVol(IDAM) * 0.5_JPRB    ! normal storage
    R_VolUpa(IDAM) = FldVol(IDAM) * 1.E-6 / upreal(IDAM)

  ELSE  !! Yamazaki&Funato scheme (paper in prep)
    Vyr =Qn(IDAM)*(365.*24.*60.*60.)                    !! Annual inflow -> assume dry period inflow is 1/8 of annual flow 
    Qsto=(ConVol(IDAM)*0.7+Vyr/4.)/(180.*24.*60.*60.)   !! possible mean outflow in dry period (6month, ConVol*0.7 + Inflow)
    Qn(IDAM)=min(Qn(IDAM),Qsto)*1.5                     !! Outflow at normal volume (*1.5 is parameter to decide outflw balance)

    AdjVol(IDAM)=ConVol(IDAM)  + FldVol(IDAM)*0.1       !! AdjVol is for outflow stability (result is not so sensitive)
    Qa(IDAM)=( Qn(IDAM)+Qf(IDAM) )*0.5                  !! Qa is also for stability
  ENDIF

  !! Year-by-Year scheme. If dam is not yet constructed
  IF( LDAMYBY )THEN
    IF( ISYYYY==DamYear(IDAM) )THEN
      DamStat(IDAM)=1   !! new this year
    ELSEIF( ISYYYY<DamYear(IDAM) .and. DamYear(IDAM)>0 )THEN
      DamStat(IDAM)=-1  !! not yet activated
      I1DAM(ISEQ)=-1
      FldVol(IDAM)=0._JPRB
      ConVol(IDAM)=0._JPRB
    ENDIF
  ENDIF

END DO
CLOSE(NDAMFILE)

WRITE(LOGNAM,*) "CMF::DAMOUT_INIT: allocated dams:", NDAMX
IF (NDAMX == 0) THEN
  WRITE(LOGNAM, *) "CMF::DAMOUT_INIT: None of the given dams allocated."
  STOP 9
ENDIF
!==========

!! mark upstream of dam grid, for applying kinematic wave routine to suppress storage buffer effect.
DO ISEQ=1, NSEQALL
  IF( I1DAM(ISEQ)<=0 .and. I1NEXT(ISEQ)>0 )THEN !! if target is non-dam grid
    JSEQ=I1NEXT(ISEQ)
    IF( I1DAM(JSEQ)==1 .or. I1DAM(JSEQ)==11 )THEN !! if downstream is dam
      I1DAM(ISEQ)=10            !! mark upstream of dam grid by "10"
      I2MASK(ISEQ,1)=1   !! reservoir upstream grid. skipped for adaptive time step
    ENDIF
  ENDIF

  IF( I1DAM(ISEQ)==1 .and. I1NEXT(ISEQ)>0 )THEN !! if target is dam grid
    JSEQ=I1NEXT(ISEQ)
    IF( I1DAM(JSEQ)==1 .or. I1DAM(JSEQ)==11 )THEN !! if downstream is dam
      I1DAM(ISEQ)=11            !! mark upstream of dam grid by "11"
      I2MASK(ISEQ,1)=2   !! reservoir grid (cascading). skipped for adaptive time step
    ENDIF
  ENDIF
END DO

!! Initialize dam storage
IF( .not. LRESTART )THEN  !! Initialize without restart data
  P2DAMSTO(:,1)=0._JPRD
  DO IDAM=1, NDAM
    IF( DamStat(IDAM)==IMIS )CYCLE !
    ISEQ=DamSeq(IDAM)
    IF( DamStat(IDAM)==-1 )THEN !! Dam not yet constructed
      P2DAMSTO(ISEQ,1)= P2RIVSTO(ISEQ,1)+P2FLDSTO(ISEQ,1)
    ELSE
      P2DAMSTO(ISEQ,1)=P2RIVSTO(ISEQ,1)+P2FLDSTO(ISEQ,1)
      IF( P2DAMSTO(ISEQ,1)<ConVol(IDAM) )THEN  !! If Normal Volume > initial storage, replace
        P2DAMSTO(ISEQ,1)=ConVol(IDAM)
        P2RIVSTO(ISEQ,1)=ConVol(IDAM)  
        P2FLDSTO(ISEQ,1)=0._JPRD
      ENDIF
    ENDIF
  END DO
ELSE       !! if from restart file
  IF( LDAMYBY )THEN   !! for restart with year-by-year option, set damsto for newly constructed dam 
    DO IDAM=1, NDAM
      IF( DamStat(IDAM)==1 )THEN !! Dam newly activated from this year
        ISEQ=DamSeq(IDAM)
        P2DAMSTO(ISEQ,1)=P2RIVSTO(ISEQ,1)+P2FLDSTO(ISEQ,1)
        IF( LiVnorm .and. P2DAMSTO(ISEQ,1)<ConVol(IDAM) )THEN  !! If Initialize Vnormal option & Vnor>Riv+Fld sto, replace
          P2DAMSTO(ISEQ,1)=ConVol(IDAM)
          P2RIVSTO(ISEQ,1)=ConVol(IDAM)  
          P2FLDSTO(ISEQ,1)=0._JPRD
        ENDIF
      ENDIF
    END DO
  ENDIF
ENDIF

!! Initialize dam inflow
DO ISEQ=1, NSEQALL
  P2DAMINF(ISEQ,1)=0._JPRD
END DO

!! Stop bifurcation at dam & dam-upstream grids
IF( LPTHOUT )THEN
  DO IPTH=1, NPTHOUT
    ISEQP=PTH_UPST(IPTH)
    JSEQP=PTH_DOWN(IPTH)
    IF( ISEQP<=0 .or. JSEQP<=0) CYCLE
    IF( I1DAM(ISEQP)>0 .or. I1DAM(JSEQP)>0 )THEN
      DO ILEV=1, NPTHLEV
        PTH_ELV(IPTH,ILEV)=1.E20  !! no bifurcation
      END DO
    ENDIF
  END DO
ENDIF

END SUBROUTINE CMF_DAMOUT_INIT_ORG
!####################################################################





!####################################################################
SUBROUTINE CMF_DAMOUT_CALC_ORG
! local
IMPLICIT NONE
! SAVE for OMP
INTEGER(KIND=JPIM),SAVE    :: ISEQD
!** dam variables
REAL(KIND=JPRB),SAVE       :: DamVol
REAL(KIND=JPRB),SAVE       :: DamInflow
REAL(KIND=JPRB),SAVE       :: DamOutflw           !! Total outflw 
REAL(KIND=JPRB),SAVE       :: DamOutTmp           !! Total outflw 
!$OMP THREADPRIVATE    (ISEQD,DamVol,DamInflow,DamOutflw,DamOutTmp)
!====================
!CONTAINS
!+ UPDATE_INFLOW: replace dam upstream with kinamatic wave, calculate inflow to dam
!+ MODIFY_OUTFLW: modify outflw to avoid negative storage
!+
!==========================================================

!* (1) Replace discharge in upstream grids with kinematic outflow
!     to avoid storage buffer effect (Shin et al., 2019, WRR)
! ------  rivout at upstream grids of dam, rivinf to dam grids are updated.
CALL UPDATE_INFLOW


!* (2) Reservoir Operation
!====================================
!     -- compare DamVol against storage level (NorVol, ConVol, EmeVol) & DamInflow against Qf
!$OMP PARALLEL DO
DO IDAM=1, NDAM
  IF( DamStat(IDAM)<=0 ) CYCLE  !! no calculation for dams not activated

  !! *** 2a update dam volume and inflow -----------------------------------
  ISEQD=DamSeq(IDAM)
  DamVol    = REAL(P2DAMSTO(ISEQD,1),KIND=JPRB)    
  DamInflow = REAL(P2DAMINF(ISEQD,1),KIND=JPRB)

  !! *** 2b Reservoir Operation          ------------------------------
  !===========================
  IF( LDAMH22 )THEN !! Hanazaki 2022 scheme
    !! case1: Water
    IF( DamVol <= NorVol(IDAM) )THEN
      DamOutflw = Qn(IDAM) * (DamVol / ConVol(IDAM) )
    !! case2: water supply
    ELSEIF( NorVol(IDAM)<DamVol .and. DamVol<=ConVol(IDAM) )THEN
      IF( Qf(IDAM)<=DamInflow )THEN
        DamOutflw = Qn(IDAM)*0.5 +   (DamVol-NorVol(IDAM))/( ConVol(IDAM)-NorVol(IDAM))      * (Qf(IDAM) - Qn(IDAM))
      ELSE
        DamOutflw = Qn(IDAM)*0.5 + (((DamVol-NorVol(IDAM))/( EmeVol(IDAM)-NorVol(IDAM)))**2) * (Qf(IDAM) - Qn(IDAM))
      ENDIF  
    !! case3: flood control
    ELSEIF( ConVol(IDAM)<DamVol .and. DamVol<EmeVol(IDAM) ) THEN
      IF( Qf(IDAM) <= DamInflow ) THEN
        DamOutflw = Qf(IDAM) + max((1. - R_VolUpa(IDAM)/0.2),0._JPRB) &
          * (DamVol-ConVol(IDAM))/(EmeVol(IDAM)-ConVol(IDAM)) * (DamInflow-Qf(IDAM))
      !! pre- and after flood control
      ELSE
        DamOutflw = Qn(IDAM)*0.5 + (((DamVol-NorVol(IDAM))/(EmeVol(IDAM)-NorVol(IDAM)))**2)* (Qf(IDAM) - Qn(IDAM))
      ENDIF
    !! case4: emergency operation
    ELSE
      DamOutflw = max(DamInflow, Qf(IDAM))
    ENDIF

  !===========================
  ELSE  !! (not LDAMH22) improved reservoir operation 'Yamazaki & Funato'

    !! Case 1: water use
    IF( DamVol<=ConVol(IDAM) )THEN
      DamOutflw = Qn(IDAM) * (DamVol/ConVol(IDAM))**0.5
    !! case 2: water excess (just avobe ConVol, for outflow stability)
    ELSEIF( DamVol>ConVol(IDAM) .and. DamVol<=AdjVol(IDAM) ) THEN
      DamOutflw = Qn(IDAM) + ( (DamVol-ConVol(IDAM)) / (AdjVol(IDAM)-ConVol(IDAM)) )**3.0 * (Qa(IDAM) - Qn(IDAM))
    !! case 3: water excess
    ELSEIF( DamVol>AdjVol(IDAM) .and. DamVol<=EmeVol(IDAM) ) THEN
      !! (flood period)
      IF( DamInflow >= Qf(IDAM) ) THEN
        !!figure left side No.2
        DamOutflw = Qn(IDAM) + ( (DamVol-ConVol(IDAM)) / (EmeVol(IDAM)-ConVol(IDAM)) )**1.0 * (DamInflow- Qn(IDAM))
        DamOutTmp = Qa(IDAM) + ( (DamVol-AdjVol(IDAM)) / (EmeVol(IDAM)-AdjVol(IDAM)) )**0.1 * (Qf(IDAM) - Qa(IDAM))
        DamOutflw = max( DamOutflw,DamOutTmp )
      !! (non-flood period)
      ELSE
        DamOutflw = Qa(IDAM) + ( (DamVol-AdjVol(IDAM)) / (EmeVol(IDAM)-AdjVol(IDAM)) )**0.1 * (Qf(IDAM) - Qa(IDAM))
      ENDIF
    !! case 4: emergency operation(no.1)
    ELSEIF( DamVol>EmeVol(IDAM) )THEN
      !! (flood period)
      IF( DamInflow >= Qf(IDAM) ) THEN
        DamOutflw = DamInflow
      !! (non-flood period)
      ELSE
        DamOutflw = Qf(IDAM)
      ENDIF
    ENDIF
  ENDIF

  !! *** 2c flow limitter
  DamOutflw = min( DamOutflw, DamVol/DT, real(P2RIVSTO(ISEQD,1)+P2FLDSTO(ISEQD,1),JPRB)/DT )
  DamOutflw = max( DamOutflw, 0._JPRB )

  !! update CaMa variables  (treat all outflow as RIVOUT in dam grid, no fldout)
  D2RIVOUT(ISEQD,1) = DamOutflw
  D2FLDOUT(ISEQD,1) = 0._JPRB
END DO
!$OMP END PARALLEL DO
!====================================

CONTAINS
!==========================================================
!+ UPDATE_INFLOW: replace dam upstream with kinamatic wave, calculate inflow to dam
!+ MODIFY_OUTFLW: modify outflw to avoid negative storage
!+
!==========================================================
SUBROUTINE UPDATE_INFLOW
USE YOS_CMF_INPUT,      ONLY: PMINSLP, PMANFLD
USE YOS_CMF_MAP,        ONLY: D2RIVLEN, D2RIVMAN, D2ELEVTN, D2NXTDST, D2RIVWTH
USE YOS_CMF_PROG,       ONLY: D2RIVOUT_PRE, D2FLDOUT_PRE
USE YOS_CMF_DIAG,       ONLY: D2RIVDPH, D2RIVVEL, D2FLDDPH
IMPLICIT NONE
! SAVE for OpenMP
INTEGER(KIND=JPIM),SAVE    :: ISEQ, JSEQ
REAL(KIND=JPRB),SAVE       :: DSLOPE,DAREA,DVEL,DSLOPE_F,DARE_F,DVEL_F
!$OMP THREADPRIVATE     (JSEQ,DSLOPE,DAREA,DVEL,DSLOPE_F,DARE_F,DVEL_F)
!============================

!*** 1a. reset outflw & dam inflow
!$OMP PARALLEL DO SIMD
DO ISEQ=1, NSEQALL
  IF( I1DAM(ISEQ)>0 )THEN  !! if dam grid or upstream of dam, reset variables
    D2RIVOUT(ISEQ,1) = 0._JPRB
    D2FLDOUT(ISEQ,1) = 0._JPRB
    P2DAMINF(ISEQ,1) = D2RUNOFF(ISEQ,1)  !! consider local runoff as inflow to dam
  ENDIF
END DO
!$OMP END PARALLEL DO SIMD

!*** 1b. calculate dam inflow, using previous tstep discharge
#ifndef NoAtom_CMF
!$OMP PARALLEL DO  !! No OMP Atomic for bit-identical simulation (set in Mkinclude)
#endif
DO ISEQ=1, NSEQALL
  IF( I1DAM(ISEQ)==10 .or. I1DAM(ISEQ)==11 )THEN  !! if dam grid or upstream of dam
    JSEQ=I1NEXT(ISEQ)
#ifndef NoAtom_CMF
!$OMP ATOMIC
#endif
    P2DAMINF(JSEQ,1) = P2DAMINF(JSEQ,1) + D2RIVOUT_PRE(ISEQ,1) + D2FLDOUT_PRE(ISEQ,1) 
  ENDIF
END DO
#ifndef NoAtom_CMF
!$OMP END PARALLEL DO
#endif

!*** 1c. discharge for upstream grids of dams
!$OMP PARALLEL DO  !! No OMP Atomic for bit-identical simulation (set in Mkinclude)
DO ISEQ=1, NSEQRIV
  IF( I1DAM(ISEQ)==10 )THEN  !! if downstream is DAM
    JSEQ   = I1NEXT(ISEQ)
    ! === river flow
    DSLOPE = (D2ELEVTN(ISEQ,1)-D2ELEVTN(JSEQ,1)) * D2NXTDST(ISEQ,1)**(-1.)
    DSLOPE = max(DSLOPE,PMINSLP)

    DVEL   = D2RIVMAN(ISEQ,1)**(-1.) * DSLOPE**0.5 * D2RIVDPH(ISEQ,1)**(2./3.)
    DAREA  = D2RIVWTH(ISEQ,1) * D2RIVDPH(ISEQ,1)

    D2RIVVEL(ISEQ,1) = DVEL
    D2RIVOUT(ISEQ,1) = DAREA * DVEL
    D2RIVOUT(ISEQ,1) = MIN( D2RIVOUT(ISEQ,1), real(P2RIVSTO(ISEQ,1),JPRB)/DT )
    !=== floodplain flow
    DSLOPE_F = min( 0.005_JPRB,DSLOPE )    !! set min [instead of using weirequation for efficiency]
    DVEL_F   = PMANFLD**(-1.) * DSLOPE_F**0.5 * D2FLDDPH(ISEQ,1)**(2./3.)
    DARE_F   = REAL(P2FLDSTO(ISEQ,1),KIND=JPRB) * D2RIVLEN(ISEQ,1)**(-1.)
    DARE_F   = MAX( DARE_F - D2FLDDPH(ISEQ,1)*D2RIVWTH(ISEQ,1), 0._JPRB )   !!remove above river channel     area

    D2FLDOUT(ISEQ,1) = DARE_F * DVEL_F
    D2FLDOUT(ISEQ,1) = MIN(  D2FLDOUT(ISEQ,1), REAL(P2FLDSTO(ISEQ,1),KIND=JPRB)/DT )
  ENDIF
END DO
!$OMP END PARALLEL DO

END SUBROUTINE UPDATE_INFLOW
!==========================================================
END SUBROUTINE CMF_DAMOUT_CALC_ORG
!####################################################################
!
!
!
!####################################################################
SUBROUTINE CMF_DAMOUT_WATBAL_ORG
IMPLICIT NONE
! SAVE for OMP
INTEGER(KIND=JPIM),SAVE    :: ISEQD
!*** water balance
REAL(KIND=JPRB),SAVE       :: DamInflow
REAL(KIND=JPRB),SAVE       :: DamOutflw           !! Total outflw 
REAL(KIND=JPRD),SAVE       :: GlbDAMSTO, GlbDAMSTONXT, GlbDAMINF, GlbDAMOUT, DamMiss
!$OMP THREADPRIVATE    (ISEQD,DamInflow,DamOutflw)
! ==========================================
!* 4) update reservoir storage and check water DamMiss --------------------------
GlbDAMSTO    = 0._JPRB
GlbDAMSTONXT = 0._JPRB
GlbDAMINF    = 0._JPRB
GlbDAMOUT    = 0._JPRB

!$OMP PARALLEL DO REDUCTION(+:GlbDAMSTO, GlbDAMSTONXT, GlbDAMINF, GlbDAMOUT)
DO IDAM=1, NDAM
  IF( DamStat(IDAM)==IMIS ) CYCLE  !! do not calculate for dams outside of calculation domain
  ISEQD = DamSeq(IDAM)

  DamInflow = D2RIVINF(ISEQD,1) + D2FLDINF(ISEQD,1) + D2RUNOFF(ISEQD,1)
  DamOutflw = D2RIVOUT(ISEQD,1) + D2FLDOUT(ISEQD,1)
!!P2DAMINF(ISEQD,1)=DamInflow   !! if water balance needs to be checked in the output file, P2DAMINF should be updated.

  GlbDAMSTO = GlbDAMSTO + P2DAMSTO(ISEQD,1)
  GlbDAMINF = GlbDAMINF + DamInflow*DT
  GlbDAMOUT = GlbDAMOUT + DamOutflw*DT

  P2DAMSTO(ISEQD,1) = P2DAMSTO(ISEQD,1) + DamInflow * DT - DamOutflw * DT

  GlbDAMSTONXT = GlbDAMSTONXT + P2DAMSTO(ISEQD,1)
END DO
!$OMP END PARALLEL DO

DamMiss = GlbDAMSTO-GlbDAMSTONXT+GlbDAMINF-GlbDAMOUT
WRITE(LOGNAM,*) "CMF::DAM_CALC: DamMiss at all dams:", DamMiss*1.D-9

END SUBROUTINE CMF_DAMOUT_WATBAL_ORG
!####################################################################
!
!
!
!####################################################################
SUBROUTINE CMF_DAMOUT_WRTE_ORG
! local
CHARACTER(len=36)          :: WriteTXT(NDAMX), WriteTXT2(NDAMX)
REAL(KIND=JPRB)            :: DDamInf, DDamOut
! File IO
INTEGER(KIND=JPIM),SAVE    :: ISEQD, JDAM
INTEGER(KIND=JPIM),SAVE    :: LOGDAM
CHARACTER(len=4),SAVE      :: cYYYY
CHARACTER(len=256),SAVE    :: CLEN, CFMT
CHARACTER(len=256),SAVE    :: DAMTXT
LOGICAL,SAVE               :: IsOpen
DATA IsOpen       /.FALSE./
! ==========================================

IF( LDAMTXT )THEN

  IF( .not. IsOpen )THEN
    IsOpen=.TRUE.
    WRITE(CYYYY,'(i4.4)') ISYYYY
    DAMTXT='./damtxt-'//trim(cYYYY)//'.txt'

    LOGDAM=INQUIRE_FID()
    OPEN(LOGDAM,FILE=DAMTXT,FORM='formatted')

    WRITE(CLEN,'(i0)') NDAMX
    CFMT="(i10,"//TRIM(CLEN)//"(a36))"

    JDAM=0
    DO IDAM=1, NDAM
      IF( DamStat(IDAM)==IMIS )CYCLE
      JDAM=JDAM+1
      IF( DamStat(IDAM)==-1 ) THEN   !! dam not activated yet
        WRITE(WriteTxt(JDAM), '(i12,2f12.2)') DamID(IDAM), -9., -9.
        WRITE(WriteTxt2(JDAM),'(3f12.2)') upreal(IDAM),   Qf(IDAM), Qn(IDAM)
      ELSE
        WRITE(WriteTxt(JDAM), '(i12,2f12.2)') DamID(IDAM), (FldVol(IDAM)+ConVol(IDAM))*1.E-9, ConVol(IDAM)*1.E-9
        WRITE(WriteTxt2(JDAM),'(3f12.2)') upreal(IDAM),   Qf(IDAM), Qn(IDAM)
      ENDIF
    END DO

    WRITE(LOGDAM,CFMT) NDAMX, (WriteTXT(JDAM) ,JDAM=1, NDAMX)

    CFMT="(a10,"//TRIM(CLEN)//"(a36))"
    WRITE(LOGDAM,CFMT)  "Date", (WriteTXT2(JDAM),JDAM=1, NDAMX)
  ENDIF

  JDAM=0
  DO IDAM=1, NDAM
    IF( DamStat(IDAM)==IMIS ) CYCLE
    JDAM=JDAM+1
    ISEQD=DamSeq(IDAM)
    DDamInf=REAL( P2DAMINF(ISEQD,1),KIND=JPRB)
    DDamOut=D2RIVOUT(ISEQD,1) + D2FLDOUT(ISEQD,1)
    WRITE(WriteTxt(JDAM), '(3f12.2)') P2DAMSTO(ISEQD,1)*1.E-9, DDamInf, DDamOut
  END DO

  CFMT="(i10,"//TRIM(CLEN)//"(a36))"  
  WRITE(LOGDAM,CFMT) IYYYYMMDD, (WriteTXT(JDAM),JDAM=1, NDAMX)

ENDIF

END SUBROUTINE CMF_DAMOUT_WRTE_ORG
!####################################################################






!####################################################################
!* CONTAINS:
! -- CMF_DAMOUT_NMLIST_IRR  : Read setting from namelist
! -- CMF_DAMOUT_INIT_IRR    : Initialize dam data
! -- CMF_DAMOUT_CALC_IRR    : Calculate inflow and outflow at dam
! -- CMF_DAMOUT_WATBAL_IRR  : Calculate water balance at dam
! -- CMF_DAMOUT_WRTE_IRR   : Write dam-related variables in text file
!####################################################################
SUBROUTINE CMF_DAMOUT_NMLIST_IRR
! reed setting from namelist
! -- Called from CMF_DRV_NMLIST
USE YOS_CMF_INPUT,      ONLY: CSETFILE,NSETFILE,LDAMOUT
IMPLICIT NONE
!================================================
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!---------------------!"

!*** 1. open namelist
NSETFILE=INQUIRE_FID()
OPEN(NSETFILE,FILE=CSETFILE,STATUS="OLD")
WRITE(LOGNAM,*) "CMF::DAMOUT_NMLIST: namelist OPEN in unit: ", TRIM(CSETFILE), NSETFILE 

!*** 2. default value
! CDAMFILE="./dam_params.csv"
LDAMTXT=.TRUE.
LDAMYBY=.FALSE.
LiVnorm=.FALSE.

!*** 3. read namelist
REWIND(NSETFILE)
READ(NSETFILE,NML=NDAMOUT)

IF( LDAMOUT )THEN
  WRITE(LOGNAM,*)   "=== NAMELIST, NDAMOUT ==="
  WRITE(LOGNAM,*)   "CDAMFILE: " , trim(CDAMFILE)
  WRITE(LOGNAM,*)   "LDAMOPT:  " , LDAMOPT
  WRITE(LOGNAM,*)   "LDAMYBY:  " , LDAMYBY
  WRITE(LOGNAM,*)   "LiVnorm:  " , LiVnorm
  WRITE(LOGNAM,*)   "LDAMTXT:  " , LDAMTXT
ENDIF

CLOSE(NSETFILE)

WRITE(LOGNAM,*) "CMF::DAMOUT_NMLIST: end" 

END SUBROUTINE CMF_DAMOUT_NMLIST_IRR
!####################################################################

!####################################################################
SUBROUTINE CMF_DAMOUT_INIT_IRR
! reed setting from CDAMFILE
IMPLICIT NONE
INTEGER(KIND=JPIM)         :: NDAMFILE,NRIVFILE
INTEGER(KIND=JPIM)         :: ISEQ, JSEQ
INTEGER(KIND=JPIM)         :: IX, IY
INTEGER(KIND=JPIM)         :: IPTH, ILEV, ISEQP, JSEQP
CHARACTER(LEN=256)         :: CDAMFILE_tmp,CRIVFILE_tmp
REAL(KIND=JPRB)            :: Qsto, Vyr
!####################################################################
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!---------------------!"
CDAMFILE_tmp = trim(CDAMFILE)//'dam_params_us_15min.csv'
WRITE(LOGNAM,*) "CMF::DAMOUT_INIT: initialize dam", trim(CDAMFILE_tmp) 
!==========
NDAMFILE=INQUIRE_FID()
OPEN(NDAMFILE,FILE=CDAMFILE_tmp,STATUS="OLD")
READ(NDAMFILE,*) NDAM
READ(NDAMFILE,*)        !! skip header

WRITE(LOGNAM,*) "CMF::DAMOUT_INIT: number of dams", NDAM

NRIV = 0
IF (LDAMIRR) THEN
  !! ================ READ irrig-river grid ================
  WRITE(LOGNAM,*) "!---------------------!"
  CRIVFILE_tmp = trim(CDAMFILE)//"WaterUse_grids/"//'river_grid_ixiy_us_15min.csv'
  WRITE(LOGNAM,*) "CMF::DAMOUT_INIT: irrig-river grid: ", trim(CRIVFILE_tmp)

  NRIVFILE=INQUIRE_FID()
  OPEN(NRIVFILE,FILE=CRIVFILE_tmp,STATUS="OLD")
  READ(NRIVFILE,*) NRIV
  READ(NRIVFILE,*)

  WRITE(LOGNAM,*) "CMF::DAMOUT_INIT: number of irrig-river grids: ", NRIV

  ALLOCATE(IX_RIV(NRIV), IY_RIV(NRIV))  
  ALLOCATE(riv_tot_demand_unmt(NRIV))
  riv_tot_demand_unmt = 0._JPRB
ENDIF

!!! update NDAM
NDAM = NDAM + NRIV

!! === ALLOCATE ===
!! from CDAMFILE
ALLOCATE(DamID(NDAM),DamName(NDAM))
ALLOCATE(DamIX(NDAM),DamIY(NDAM),DamLon(NDAM),DamLat(NDAM))
ALLOCATE(MainUse(NDAM))
ALLOCATE(upreal(NDAM))

!! calculate from CDAMFILE
ALLOCATE(DamSeq(NDAM))
ALLOCATE(Qf(NDAM),Qn(NDAM),Qe(NDAM))
ALLOCATE(TotVol(NDAM),FldVol(NDAM),NorVol(NDAM),ConVol(NDAM),EmeVol(NDAM))

!! dam map, dam variable
ALLOCATE(I1DAM(NSEQMAX))

!! for outflw stability
ALLOCATE(AdjVol(NDAM),Qa(NDAM))
ALLOCATE(DamYear(NDAM),DamStat(NDAM))
!! =================
DamSeq(:) =IMIS
DamStat(:)=IMIS
I1DAM(:)=0
NDAMX=0
NRIVX= 0

!! read dam parameters
DO IDAM = 1, NDAM

  IF (LDAMIRR .and. (IDAM <= NRIV) ) then
    !! ================ READ irrig-river grid ================
    ! write(LOGNAM,*) "CMF::DAMOUT_INIT: READ irrig-river grid: IDAM = ", IDAM
    READ(NRIVFILE,*) IX_RIV(IDAM), IY_RIV(IDAM)
    IX = IX_RIV(IDAM)
    IY = IY_RIV(IDAM)
    NorVol(IDAM) = 0._JPRB

    IF (IX<=0 .or. IX > NX .or. IY<=0 .or. IY > NY ) cycle
    ISEQ=I2VECTOR(IX,IY)  
    IF( I1NEXT(ISEQ)==-9999 .or. ISEQ<=0 ) cycle
    NRIVX=NRIVX+1
  ELSE
    !! ================ READ dam basic information ================
    ! write(LOGNAM,*) "CMF::DAMOUT_INIT: READ dam basic information: IDAM = ", IDAM - NRIV
    READ(NDAMFILE,*) DamID(IDAM),DamName(IDAM),DamLon(IDAM),DamLat(IDAM),MainUse(IDAM),DamYear(IDAM),DamIX(IDAM),DamIY(IDAM), &
      upreal(IDAM),Qn(IDAM),Qf(IDAM),TotVol(IDAM),FldVol(IDAM),NorVol(IDAM),ConVol(IDAM)

    Qe(IDAM) = Qn(IDAM) * 0.1_JPRB  
    !! --------------

    IX=DamIX(IDAM)
    IY=DamIY(IDAM)
    IF (IX<=0 .or. IX > NX .or. IY<=0 .or. IY > NY ) cycle
    ISEQ=I2VECTOR(IX,IY)  
    IF( I1NEXT(ISEQ)==-9999 .or. ISEQ<=0 ) cycle
    NDAMX=NDAMX+1
  ENDIF

  DamSeq(IDAM) =ISEQ
  DamStat(IDAM)=2
  I1DAM(ISEQ)=1
  I2MASK(ISEQ,1)=2   !! reservoir grid. skipped for adaptive time step
ENDDO
CLOSE(NDAMFILE)
CLOSE(NRIVFILE)

WRITE(LOGNAM,*) "CMF::DAMOUT_INIT: allocated riv grids:", NRIVX 
WRITE(LOGNAM,*) "CMF::DAMOUT_INIT: allocated dams:", NDAMX 

IF (NDAMX == 0) THEN
  WRITE(LOGNAM, *) "CMF::DAMOUT_INIT: None of the given dams allocated."
  STOP 9
ENDIF

  !! storage parameter --- from Million Cubic Meter to m3
  TotVol = TotVol * 1.E6_JPRB
  FldVol = FldVol * 1.E6_JPRB  ! Flood control storage capacity: exclusive for flood control
  NorVol = NorVol * 1.E6_JPRB
  ConVol = ConVol * 1.E6_JPRB
DO IDAM = 1, NDAM
  EmeVol(IDAM) = ConVol(IDAM) + FldVol(IDAM) * 0.95_JPRB  ! storage to start emergency operation
ENDDO

!! ================ READ water use data ================
IF (LDAMIRR) THEN
  WRITE(LOGNAM,*) "CMF::DAMOUT_INIT: READ_WUSE_GRID "
  CALL READ_WUSE_GRID
END IF

IF(LDAMOPT == "H06" .OR. LDAMOPT == "V13")THEN  
  WRITE(LOGNAM,*) "CMF::DAMOUT_INIT: READ_WATER_USE "
  CALL READ_WATER_USE
ENDIF


!! ================ parameters for different schemes ================
IF(LDAMOPT == "H06" .OR. LDAMOPT == "V13")THEN  
  !! --- ALLOCATE ---
  ALLOCATE(H06_DPI(NDAM), H06_c(NDAM))
  H06_DPI = WUSE_AD / Qn    
  H06_c = TotVol / (Qn * 86400. * 365.) !! Qn m3/s to m3
ENDIF

IF(LDAMOPT == "H22")THEN 
  !! --- ALLOCATE ---
  ALLOCATE(H22_EmeVol(NDAM), H22_k(NDAM))
  H22_EmeVol = FldVol + (TotVol - FldVol)*0.2 
  H22_k = max((1. - (TotVol - FldVol) / upreal /0.2),0._JPRB)
ENDIF

IF(LDAMOPT == "V13")THEN  
  WRITE(LOGNAM,*) "!---------------------!"
  CDAMFILE_tmp = trim(CDAMFILE)//'damfcperiod.csv'
  WRITE(LOGNAM,*) "CMF::DAMOUT_INIT: dam flow parameters:", trim(CDAMFILE_tmp) 
  NDAMFILE=INQUIRE_FID()
  OPEN(NDAMFILE,FILE=CDAMFILE_tmp,STATUS="OLD")
  READ(NDAMFILE,*) 
  !! --- ALLOCATE ---
  ALLOCATE(StFC_Mth(NDAM),NdFC_Mth(NDAM),StOP_Mth(NDAM))
  DO IDAM = 1, NDAM
    READ(NDAMFILE,*) DamID(IDAM),StFC_Mth(IDAM),NdFC_Mth(IDAM),StOP_Mth(IDAM)
  ENDDO
  CLOSE(NDAMFILE)
ENDIF

IF(LDAMOPT == "YF")THEN
  DO IDAM = 1, NDAM
    Vyr =Qn(IDAM)*(365.*24.*60.*60.)                    !! Annual inflow -> assume dry period inflow is 1/8 of annual flow 
    Qsto=(ConVol(IDAM)*0.7+Vyr/4.)/(180.*24.*60.*60.)   !! possible mean outflow in dry period (6month, ConVol*0.7 + Inflow)
    Qn(IDAM)=min(Qn(IDAM),Qsto)*1.5                     !! Outflow at normal volume (*1.5 is parameter to decide outflw balance)
    AdjVol(IDAM)=ConVol(IDAM)  + FldVol(IDAM)*0.1       !! AdjVol is for outflow stability (result is not so sensitive)
    Qa(IDAM)=( Qn(IDAM)+Qf(IDAM) )*0.5                  !! Qa is also for stability
  ENDDO                
ENDIF

!! Year-by-Year scheme. If dam is not yet constructed
IF( LDAMYBY )THEN
  DO IDAM = 1, NDAM
    IF( ISYYYY==DamYear(IDAM) )THEN
      DamStat(IDAM)=1   !! new this year
    ELSEIF( ISYYYY<DamYear(IDAM) .and. DamYear(IDAM)>0 )THEN
      DamStat(IDAM)=-1  !! not yet activated
      I1DAM(ISEQ)=-1
      FldVol(IDAM)=0._JPRB
      ConVol(IDAM)=0._JPRB
    ENDIF
  ENDDO
ENDIF
!==========

!! mark upstream of dam grid, for applying kinematic wave routine to suppress storage buffer effect.
DO ISEQ=1, NSEQALL
  IF( I1DAM(ISEQ)<=0 .and. I1NEXT(ISEQ)>0 )THEN !! if target is non-dam grid
    JSEQ=I1NEXT(ISEQ)
    IF( I1DAM(JSEQ)==1 .or. I1DAM(JSEQ)==11 )THEN !! if downstream is dam
      I1DAM(ISEQ)=10            !! mark upstream of dam grid by "10"
      I2MASK(ISEQ,1)=1   !! reservoir upstream grid. skipped for adaptive time step
    ENDIF
  ENDIF

  IF( I1DAM(ISEQ)==1 .and. I1NEXT(ISEQ)>0 )THEN !! if target is dam grid
    JSEQ=I1NEXT(ISEQ)
    IF( I1DAM(JSEQ)==1 .or. I1DAM(JSEQ)==11 )THEN !! if downstream is dam
      I1DAM(ISEQ)=11            !! mark upstream of dam grid by "11"
      I2MASK(ISEQ,1)=2          !! reservoir grid (cascading). skipped for adaptive time step
    ENDIF
  ENDIF
ENDDO

!! Initialize dam storage
IF( .not. LRESTART )THEN  !! Initialize without restart data
  P2DAMSTO(:,1)=0._JPRD
  DO IDAM=1, NDAM
    IF( DamStat(IDAM)==IMIS )CYCLE
    IF( DamSeq(IDAM)>0 )THEN
      ISEQ=DamSeq(IDAM)
      P2DAMSTO(ISEQ,1)=NorVol(IDAM)*0.5  !! set initial storage to Normal Storage Volume
      P2RIVSTO(ISEQ,1)=max(P2RIVSTO(ISEQ,1),NorVol(IDAM)*0.5_JPRD) !! also set initial river storage, in order to keep consistency
    ENDIF
  ENDDO
ELSE
  IF( LDAMYBY )THEN
    DO IDAM=1, NDAM
      IF( DamStat(IDAM)==1 )THEN !! Dam newly activated from this year
        ISEQ=DamSeq(IDAM)
        P2DAMSTO(ISEQ,1)=P2RIVSTO(ISEQ,1)+P2FLDSTO(ISEQ,1)
        IF( LiVnorm .AND. P2DAMSTO(ISEQ,1)<ConVol(IDAM) )THEN  !! If Initialize Vnormal option & Vnor>Riv+Fld sto, replace
          P2DAMSTO(ISEQ,1)=ConVol(IDAM)
          P2RIVSTO(ISEQ,1)=ConVol(IDAM)  
          P2FLDSTO(ISEQ,1)=0._JPRD
        ENDIF
      ENDIF
    ENDDO
  ENDIF
ENDIF

!! Initialize dam inflow
DO ISEQ=1, NSEQALL
  P2DAMINF(ISEQ,1)=0._JPRD
ENDDO

!! Stop bifurcation at dam & dam-upstream grids
IF( LPTHOUT )THEN
  DO IPTH=1, NPTHOUT
    ISEQP=PTH_UPST(IPTH)
    JSEQP=PTH_DOWN(IPTH)
    IF( ISEQP<=0 .or. JSEQP<=0) CYCLE
    IF( I1DAM(ISEQP)>0 .or. I1DAM(JSEQP)>0 )THEN
      DO ILEV=1, NPTHLEV
        PTH_ELV(IPTH,ILEV)=1.E20  !! no bifurcation
      ENDDO
    ENDIF
  ENDDO
ENDIF

!####################################################################
! !print*,"LHB debug line209 CaMa dam init error : dam start finish"
CONTAINS

SUBROUTINE READ_WUSE_GRID
  IMPLICIT NONE
  ! local variables
  CHARACTER(LEN=256)                         :: CDAMFILE_WUSE_GRID
  ! read ix
  CDAMFILE_WUSE_GRID = TRIM(CDAMFILE)//"WaterUse_grids/"//'irrig_ix_us_15min.txt'
  NDAMFILE=INQUIRE_FID()
  OPEN(NDAMFILE, file=TRIM(CDAMFILE_WUSE_GRID), status='old', form='formatted')
  READ(NDAMFILE, *)
  READ(NDAMFILE, *) max_gridnum

  ALLOCATE(serial_num(NDAM),dam_id(NDAM),grid_num(NDAM),area_1(NDAM),area_2(NDAM))
  ALLOCATE(grids_x(NDAM,max_gridnum))
  ALLOCATE(grids_y(NDAM,max_gridnum))
  ALLOCATE(grids_share(NDAM,max_gridnum))

  grids_x = 0
  grids_y = 0 
  grids_share = 0._JPRB

  DO IDAM = NRIV+1, NDAM
    READ(NDAMFILE, *) serial_num(IDAM), dam_id(IDAM), area_1(IDAM), area_2(IDAM), grid_num(IDAM), grids_x(IDAM,1:grid_num(IDAM))
  ENDDO
  CLOSE(NDAMFILE)

  ! read iy
  CDAMFILE_WUSE_GRID = TRIM(CDAMFILE)//"WaterUse_grids/"//'irrig_iy_us_15min.txt'
  NDAMFILE=INQUIRE_FID()
  OPEN(NDAMFILE, file=TRIM(CDAMFILE_WUSE_GRID), status='old', form='formatted')
  READ(NDAMFILE, *)
  READ(NDAMFILE, *) 
  DO IDAM = NRIV+1, NDAM
    READ(NDAMFILE, *) serial_num(IDAM), dam_id(IDAM), area_1(IDAM), area_2(IDAM), grid_num(IDAM), grids_y(IDAM,1:grid_num(IDAM))
  ENDDO
  CLOSE(NDAMFILE)
  
  ! read share
  CDAMFILE_WUSE_GRID = TRIM(CDAMFILE)//"WaterUse_grids/"//'irrig_grid_share_us_15min.txt'
  NDAMFILE=INQUIRE_FID()
  OPEN(NDAMFILE, file=TRIM(CDAMFILE_WUSE_GRID), status='old', form='formatted')
  READ(NDAMFILE, *)
  READ(NDAMFILE, *) 
  DO IDAM = NRIV+1, NDAM
    READ(NDAMFILE, *) serial_num(IDAM), dam_id(IDAM), area_1(IDAM), area_2(IDAM), grid_num(IDAM), grids_share(IDAM,1:grid_num(IDAM))
  ENDDO
  CLOSE(NDAMFILE)


  !! ----- add allocation for irrigation water withdraw -----
  ALLOCATE(dam_tot_demand(NDAM))
  ALLOCATE(dam_tot_demand_save(NDAM))
  ALLOCATE(dam_grid_demand(NDAM, max_gridnum))
  ALLOCATE(save_damwithdraw(NDAM))
  !! ----- add allocation for irrigation water withdraw -----

END SUBROUTINE READ_WUSE_GRID

!####################################################################
SUBROUTINE READ_WATER_USE
  USE YOS_CMF_TIME,       ONLY: NSTEPS   

  IMPLICIT NONE
  CHARACTER(LEN=256)         :: CDAMFILE_WUSE_YEAR
  CHARACTER(len=4)           :: CYYYY
  CHARACTER(LEN=16)          :: tmp_i, tmp_name  
  INTEGER(KIND=JPIM)         :: NDAMFILE

  !=======================================  
  WRITE(CYYYY,'(I4)') ISYYYY
  CDAMFILE_WUSE_YEAR = TRIM(CDAMFILE)//'water_use/'//TRIM(adjustl(CYYYY))//'.txt'
  WRITE(LOGNAM,*) "CMF::DAMOUT_INIT: dam water use file:", TRIM(CDAMFILE_WUSE_YEAR)

  !! --- ALLOCATE ---
  !! from dam water use data
  ALLOCATE(WUSE_DD(NDAM, NSTEPS))    
  ALLOCATE(WUSE_AD(NDAM))          

  !! read dam water use
  NDAMFILE=INQUIRE_FID()
  OPEN(NDAMFILE,FILE=TRIM(CDAMFILE_WUSE_YEAR),STATUS="OLD")

  DO IDAM = 1, NDAM
    READ(NDAMFILE,*) tmp_i, tmp_name, WUSE_DD(IDAM,:)  
    WUSE_AD(IDAM) = SUM(WUSE_DD(IDAM,:))/NSTEPS     
  END DO
  CLOSE(NDAMFILE)

END SUBROUTINE READ_WATER_USE
END SUBROUTINE CMF_DAMOUT_INIT_IRR
!####################################################################


!####################################################################
SUBROUTINE CMF_DAMOUT_WATBAL_IRR
  IMPLICIT NONE

  ! SAVE for OMP
  INTEGER(KIND=JPIM),SAVE    :: ISEQD
  !*** water balance
  REAL(KIND=JPRB),SAVE       :: DamInflow
  REAL(KIND=JPRB),SAVE       :: DamOutflw           !! Total outflw 
  REAL(KIND=JPRD),SAVE       :: GlbDAMSTO, GlbDAMSTONXT, GlbDAMINF, GlbDAMOUT, DamMiss
  
  !$OMP THREADPRIVATE    (ISEQD,DamInflow,DamOutflw)
  ! ==========================================
  !* 4) update reservoir storage and check water DamMiss --------------------------
  GlbDAMSTO    = 0._JPRB
  GlbDAMSTONXT = 0._JPRB
  GlbDAMINF    = 0._JPRB
  GlbDAMOUT    = 0._JPRB
  
  !$OMP PARALLEL DO REDUCTION(+:GlbDAMSTO, GlbDAMSTONXT, GlbDAMINF, GlbDAMOUT)
  DO IDAM=1, NDAM
    IF( DamSeq(IDAM)<=0 ) CYCLE
    ISEQD = DamSeq(IDAM)
  
    DamInflow = D2RIVINF(ISEQD,1) + D2FLDINF(ISEQD,1) + D2RUNOFF(ISEQD,1)
    DamOutflw = D2RIVOUT(ISEQD,1) + D2FLDOUT(ISEQD,1)
  
    GlbDAMSTO = GlbDAMSTO + P2DAMSTO(ISEQD,1)
    GlbDAMINF = GlbDAMINF + DamInflow*DT
    GlbDAMOUT = GlbDAMOUT + DamOutflw*DT
  
    P2DAMSTO(ISEQD,1) = P2DAMSTO(ISEQD,1) + DamInflow * DT - DamOutflw * DT

    GlbDAMSTONXT = GlbDAMSTONXT + P2DAMSTO(ISEQD,1)
  ENDDO
  !$OMP END PARALLEL DO
  
  DamMiss = GlbDAMSTO-GlbDAMSTONXT+GlbDAMINF-GlbDAMOUT
END SUBROUTINE CMF_DAMOUT_WATBAL_IRR
!####################################################################

!####################################################################
SUBROUTINE CMF_DAM_WUSE_UPDATE
  IMPLICIT NONE
  ! local variables
  INTEGER(KIND=JPIM) :: i_idam
  INTEGER(KIND=JPIM) :: i
  !##############################################

  !###### gridded daily irr_demand to dam-scale daily irr_demand
  dam_grid_demand(:,:) = 0._JPRB
  dam_tot_demand(:) = 0._JPRB

  DO i_idam = NRIV+1, NDAM
    DO i = 1, grid_num(i_idam)
        dam_grid_demand(i_idam, i) = grids_share(i_idam, i) * dirrig_cama(grids_x(i_idam, i), grids_y(i_idam, i))
    ENDDO    
    dam_tot_demand(i_idam) = sum(dam_grid_demand(i_idam, :))
  ENDDO

  ! intialize 
  dam_tot_demand_save(:) = 0._JPRB
  save_damwithdraw(:)    = 0._JPRB
END SUBROUTINE CMF_DAM_WUSE_UPDATE

!#################################################################### added!!!
SUBROUTINE CMF_DAM_WUSE_ALLOC
  IMPLICIT NONE
  ! local
  REAL(KIND=JPRB), ALLOCATABLE           :: release_water_temp(:,:)
  REAL(KIND=JPRB), ALLOCATABLE           :: release_cama_dam_temp(:,:)
  INTEGER(KIND=JPIM)                     :: i,i_idam

  !####################################################################
  ALLOCATE(release_water_temp(NX,NY))
  ALLOCATE(release_cama_dam_temp(NX,NY))
  release_cama_dam_temp(:,:) = 0._JPRB

  DO i_idam = NRIV+1, NDAM
    IF (dam_tot_demand_save(i_idam) .GT. 0.0) then
      release_water_temp(:,:) = 0._JPRB
      DO i = 1, grid_num(i_idam)        
        release_water_temp (grids_x(i_idam,i), grids_y(i_idam,i)) = (dam_grid_demand(i_idam,i) / & 
            dam_tot_demand_save(i_idam)) * save_damwithdraw(i_idam)
      ENDDO

      release_cama_dam_temp = release_cama_dam_temp + release_water_temp   ! unit: m3/day
    ENDIF
  ENDDO
  
  DEALLOCATE(release_water_temp)

  release_cama_dam_temp = MIN(dirrig_cama, release_cama_dam_temp)
  ! update dirrig_cama
  dirrig_cama = dirrig_cama - release_cama_dam_temp
  ! update release_cama_dam
  release_cama_dam = release_cama_dam + release_cama_dam_temp

END SUBROUTINE CMF_DAM_WUSE_ALLOC

!####################################################################
SUBROUTINE CMF_DAM_DEMAND_UPDATE
  USE YOS_CMF_TIME,       only: IHHMM, JMM, JDD
  IMPLICIT NONE
  ! 每天600时刻重新初始化dirrig_cama
  IF ((JMM == 12) .AND. (JDD == 31) .AND. (IHHMM == 1800)) THEN
    dirrig_cama_unmt = 0._JPRB
  ELSE
    dirrig_cama_unmt = dirrig_cama
  ENDIF
END SUBROUTINE CMF_DAM_DEMAND_UPDATE

!####################################################################
SUBROUTINE CMF_DAMOUT_CALC_IRR
USE YOS_CMF_TIME,       ONLY: KSTEP, ISMM 
! local
IMPLICIT NONE
! SAVE for OMP
INTEGER(KIND=JPIM),SAVE    :: ISEQD
!** dam variables
REAL(KIND=JPRB),SAVE       :: DamVol
REAL(KIND=JPRB),SAVE       :: DamInflow
REAL(KIND=JPRB),SAVE       :: DamOutflw           !! Total outflw 
REAL(KIND=JPRB),SAVE       :: DamWithdraw_temp  !  dam_demand
REAL(KIND=JPRB),SAVE       :: riv_demand_temp, RivWithdraw_temp  !  riv_demand
!$OMP THREADPRIVATE    (ISEQD,DamVol,DamInflow,DamOutflw)
!====================
!CONTAINS
!+ UPDATE_INFLOW: replace dam upstream with kinamatic wave, calculate inflow to dam
!+ MODIFY_OUTFLW: modify outflw to avoid negative storage
!==========================================================

!* (1) Replace discharge in upstream grids with kinematic outflow
!     to avoid storage buffer effect (Shin et al., 2019, WRR)
! ------  rivout at upstream grids of dam, rivinf to dam grids are updated.
CALL UPDATE_INFLOW

!* (2) Reservoir Operation
!====================================
!     -- compare DamVol against storage level (NorVol, ConVol, EmeVol) & DamInflow against Qf
!! $OMP PARALLEL DO
DO IDAM=1, NDAM
  IF( DamSeq(IDAM)<=0 ) CYCLE  !! no calculation for dams not activated

  !! *** 2a update dam volume and inflow -----------------------------------
  ISEQD     = DamSeq(IDAM)
  DamVol    = REAL(P2DAMSTO(ISEQD,1),KIND=JPRB)    
  DamInflow = REAL(P2DAMINF(ISEQD,1),KIND=JPRB)

  IF( LDAMIRR .and. (IDAM <= NRIV)) THEN

    !!! =========================== river water use ===========================    
    DamOutflw = DamInflow
    riv_demand_temp = dirrig_cama(IX_RIV(IDAM), IY_RIV(IDAM)) / DT   ! unit: m3 to m3/s
    IF (riv_demand_temp .GT. 0.0) THEN
      RivWithdraw_temp = MIN(MAX(0._JPRB, DamOutflw * 0.8), riv_demand_temp)
      DamOutflw = DamOutflw - RivWithdraw_temp 
      !! add release and reduce demand
      release_cama_riv(IX_RIV(IDAM), IY_RIV(IDAM)) = release_cama_riv(IX_RIV(IDAM), IY_RIV(IDAM)) + RivWithdraw_temp * DT
      dirrig_cama(IX_RIV(IDAM), IY_RIV(IDAM)) = MAX(dirrig_cama(IX_RIV(IDAM), IY_RIV(IDAM)) - RivWithdraw_temp * DT, 0._JPRB)
    ENDIF

  ELSE
    IF (LDAMIRR) THEN
      !!! ============================== dam w ater use =============================      
      !!! ============================== step1: update water demand for each dam
      IF (IDAM == NRIV+1) THEN         
        CALL CMF_DAM_WUSE_UPDATE
      ENDIF

      !!! ============================== step2: dam water use
      IF (dam_tot_demand(IDAM) .GT. 0.0) THEN
        DamWithdraw_temp = MIN(MAX(0._JPRB, DamVol - ConVol(IDAM) * 0.5), dam_tot_demand(IDAM))  ! ConVol: conservation storage mcm m3    
        DamVol = DamVol - DamWithdraw_temp
        P2DAMSTO(ISEQD,1) = DamVol
        save_damwithdraw(IDAM) = MIN(DamWithdraw_temp, dam_tot_demand(IDAM)) !+ save_damwithdraw(IDAM)
        dam_tot_demand_save(IDAM) = dam_tot_demand(IDAM)
        dam_tot_demand(IDAM) = dam_tot_demand(IDAM) - DamWithdraw_temp
      ENDIF
    ENDIF

  !! *** 2b Reservoir Operation          ------------------------------
  !===========================
!! option: Voisin 2013 scheme
    IF( LDAMOPT == "V13" )THEN
      CALL DAM_OPERATION_V13(DamVol, DamInflow, DamOutflw, &
                            MainUse(IDAM), TotVol(IDAM), Qn(IDAM),H06_DPI(IDAM), H06_c(IDAM), &
                            WUSE_DD(IDAM, KSTEP), WUSE_AD(IDAM),&
                            StFC_Mth(IDAM), NdFC_Mth(IDAM), StOP_Mth(IDAM))
    ENDIF

    !! option: LISFLOOD scheme
    IF( LDAMOPT == "LIS" )THEN
      CALL DAM_OPERATION_LISFLOOD(DamVol, DamInflow, DamOutflw, &
                                  ConVol(IDAM), NorVol(IDAM), FldVol(IDAM), TotVol(IDAM), Qn(IDAM), Qf(IDAM))
    ENDIF

    !! option: Hanazaki 2022 scheme
    IF( LDAMOPT == "H22" )THEN
      CALL DAM_OPERATION_H22(DamVol, DamInflow, DamOutflw, &
                            FldVol(IDAM)*0.5, FldVol(IDAM), H22_EmeVol(IDAM), Qn(IDAM), Qf(IDAM),H22_k(IDAM))
    ENDIF

    IF( LDAMOPT == "YF" )THEN
      CALL DAM_OPERATION_YF(DamVol, DamInflow, DamOutflw, &
                            ConVol(IDAM), EmeVol(IDAM), AdjVol(IDAM), Qn(IDAM), Qf(IDAM), Qa(IDAM))
    ENDIF

    !! *** 2c flow limitter
    DamOutflw = MIN( DamOutflw, DamVol/DT, real(P2RIVSTO(ISEQD,1)+P2FLDSTO(ISEQD,1),JPRB)/DT )
    DamOutflw = MAX( DamOutflw, 0._JPRB )
  ENDIF

  !! update CaMa variables  (treat all outflow as RIVOUT in dam grid, no fldout)
  !print*,"LHB debug line771 camarun error : update CaMa variables"
  D2RIVOUT(ISEQD,1) = DamOutflw
  D2FLDOUT(ISEQD,1) = 0._JPRB
ENDDO
!!$OMP END PARALLEL DO

!====================================
if (LDAMIRR) then
  CALL CMF_DAM_WUSE_ALLOC
end if

!* 3) modify outflow to suppless negative discharge, update RIVOUT,FLDOUT,RIVINF,FLDINF
CALL MODIFY_OUTFLW

!* 4) update reservoir storage and check water DamMiss --------------------------
CALL CMF_DAMOUT_WATBAL_IRR

CONTAINS
!==========================================================
!+ UPDATE_INFLOW: replace dam upstream with kinamatic wave, calculate inflow to dam
!+ MODIFY_OUTFLW: modify outflw to avoid negative storage
!+ DAM_OPERATION_H06
!+ DAM_OPERATION_V13
!+ DAM_OPERATION_LIS
!+ DAM_OPERATION_H22
!==========================================================
SUBROUTINE DAM_OPERATION_H06(Vol_dam, inflw, outflw, &
                             dam_purpose, Vol_tot, Q_n, DPI, c, &
                             dam_wuse, dam_wuse_avg)
  !! reference: Hanasaki, N., Kanae, S., & Oki, T. (2006).
  !! A reservoir operation scheme for global river routing models. 
  !! Journal of Hydrology, 327(1-2), 22-41.

  IMPLICIT NONE
  REAL(KIND=JPRB), INTENT(IN)             :: Vol_dam       ! Current volume of the dam (m3)
  REAL(KIND=JPRB), INTENT(IN)             :: inflw         ! Inflow rate to the dam (m3/s)
  REAL(KIND=JPRB), INTENT(OUT)            :: outflw        ! Outflow rate from the dam (m3/s)
  !*** parameter
  CHARACTER(LEN=256), INTENT(IN)          :: dam_purpose   ! main use of dam
  REAL(KIND=JPRB), INTENT(IN)             :: Vol_tot       ! total storage capacity (m3)
  REAL(KIND=JPRB), INTENT(IN)             :: Q_n           ! normal discharge (m3/s)
  REAL(KIND=JPRB), INTENT(IN)             :: DPI           ! the ratio between annual mean demand and annual mean inflo (-)
  REAL(KIND=JPRB), INTENT(IN)             :: c             ! the ratio between capacity and mean annual inflow (-)
  REAL(KIND=JPRB), INTENT(IN)             :: dam_wuse      ! daily water demand (m3/s)
  REAL(KIND=JPRB), INTENT(IN)             :: dam_wuse_avg  ! average daily demand (m3/s)
  !*** local
  REAL(KIND=JPRB)                         :: a=0.85            !! adjustment factor for Krls
  REAL(KIND=JPRB)                         :: M=0.5             !! the ratio between minimum release and long‐term annual mean inflow
  REAL(KIND=JPRB)                         :: R             !! demand‐controlled release ratio = min(1, a*c)
  REAL(KIND=JPRB)                         :: Krls          !! the ratio between initial storage and the long‐term target storage = S/(a*C)
  REAL(KIND=JPRB)                         :: OutTmp        !! temporary outflow
  !=====================================================
  ! !! parameter
  ! M = 0.5
  ! a = 0.85 
  !! parameter calculation
  Krls =  Vol_dam / (Vol_tot * a)
  R    =  MIN(1., a * c)

  !! calculate DamOutTmp
  IF( dam_purpose == "Irrigation" )THEN
    !! irrigation reservoir
    IF ( DPI < (1. - M) ) THEN
      OutTmp = Q_n + dam_wuse - dam_wuse_avg
    ELSE
      OutTmp = Q_n * (M + (1. - M) * dam_wuse / dam_wuse_avg)
    ENDIF
  ELSE
    !! not irrigation reservoir
    OutTmp = Q_n 
  ENDIF

  !! calculate DamOutflw
  IF ( c >= 0.5 ) THEN
    outflw = Krls * OutTmp 
  ELSE
    outflw = R * Krls * OutTmp + (1. - R) * inflw
  ENDIF
END SUBROUTINE DAM_OPERATION_H06

SUBROUTINE DAM_OPERATION_V13(Vol_dam, inflw, outflw, &
                             dam_purpose, Vol_tot, Q_n, DPI, c, &
                             dam_wuse, dam_wuse_avg, &
                             MthStFC,MthNdFC,MthStOP)
  !! reference: Voisin, N., Li, H., Ward, D., Huang, M., Wigmosta, M., & Leung, L. R. (2013). 
  !! On an improved sub-regional water resources management representation for integration into earth system models. 
  !! Hydrology and Earth System Sciences, 17(9), 3605-3622.
  !! reference: https://github.com/IMMM-SFA/mosartwmpy

  IMPLICIT NONE
  REAL(KIND=JPRB), INTENT(IN)             :: Vol_dam       ! Current volume of the dam (m3)
  REAL(KIND=JPRB), INTENT(IN)             :: inflw         ! Inflow rate to the dam (m3/s)
  REAL(KIND=JPRB), INTENT(OUT)            :: outflw        ! Outflow rate from the dam (m3/s)
  !*** parameter
  CHARACTER(LEN=256), INTENT(IN)          :: dam_purpose   ! main use of dam
  REAL(KIND=JPRB), INTENT(IN)             :: Vol_tot       ! total storage capacity (m3)
  REAL(KIND=JPRB), INTENT(IN)             :: Q_n           ! normal discharge (m3/s)
  REAL(KIND=JPRB), INTENT(IN)             :: DPI           ! the ratio between annual mean demand and annual mean inflo (-)
  REAL(KIND=JPRB), INTENT(IN)             :: c             ! the ratio between capacity and mean annual inflow (-)
  REAL(KIND=JPRB), INTENT(IN)             :: dam_wuse      ! daily water demand (m3/s)
  REAL(KIND=JPRB), INTENT(IN)             :: dam_wuse_avg  ! average daily demand (m3/s)
  INTEGER(KIND=JPIM), INTENT(IN)          :: MthStFC       ! the start of the flood control period
  INTEGER(KIND=JPIM), INTENT(IN)          :: MthNdFC       ! the end of the flood control period
  INTEGER(KIND=JPIM), INTENT(IN)          :: MthStOP       ! the start of the operational year
  !*** local
  REAL(KIND=JPRB)                         :: a=0.85            !! adjustment factor for Krls
  REAL(KIND=JPRB)                         :: M=0.5             !! the ratio between minimum release and long‐term annual mean inflow
  REAL(KIND=JPRB)                         :: R             !! demand‐controlled release ratio = min(1, a*c)
  REAL(KIND=JPRB)                         :: Krls          !! the ratio between initial storage and the long‐term target storage = S/(a*C)
  REAL(KIND=JPRB)                         :: OutTmp        !! temporary outflow
  REAL(KIND=JPRB)                         :: drop          !! temporary drop flow
  !=====================================================
  !! parameter calculation
  Krls =  Vol_dam / (Vol_tot * a)
  R    =  MIN(1., a * c)

  IF( dam_purpose == "Irrigation" .or. dam_purpose == "Flood-control") THEN
    !! irrigation reservoir
    !! calculate OutTmp
    IF ( DPI < (1. - M) ) THEN
      OutTmp = Q_n + dam_wuse - dam_wuse_avg
    ELSE
      OutTmp = Q_n * (M + (1. - M) * dam_wuse / dam_wuse_avg)
    ENDIF

    !! conmbined irrigation & flood control
    !*** stage-1: from MthStFC to MthNdFC
    IF ( MthStFC <= MthNdFC ) THEN
      IF ( ISMM >= MthStFC .AND. ISMM < MthNdFC) THEN
        IF ( inflw < Q_n ) THEN
          drop = abs(inflw - Q_n)
          OutTmp = OutTmp + drop
        ENDIF
      ENDIF
    ELSEIF ( MthStFC > MthNdFC ) THEN
      IF ( ISMM >= MthStFC .OR. ISMM < MthNdFC) THEN
        IF ( inflw < Q_n ) THEN
          drop = abs(inflw - Q_n)
          OutTmp = OutTmp + drop
        ENDIF
      ENDIF
    ENDIF

    !*** stage-2: from MthNdFC to MthStOP
    IF ( MthNdFC <= MthStOP ) THEN
      IF ( ISMM >= MthNdFC .AND. ISMM < MthStOP ) THEN
          ! IF ( inflw > Q_n ) THEN
          !   fill = abs(inflw - Q_n)
          ! END
          IF (OutTmp > Q_n) THEN
            OutTmp = Q_n
          ENDIF
      ENDIF
    elseif ( MthNdFC > MthStOP ) THEN
      IF ( ISMM >= MthNdFC .or. ISMM < MthStOP ) THEN
          ! IF ( inflw > Q_n ) THEN
          !   fill = abs(inflw - Q_n)
          ! END
          IF (OutTmp > Q_n) THEN
            OutTmp = Q_n
          ENDIF
      ENDIF
    ENDIF
  !! not irrigation reservoir
  ELSE    
    OutTmp = Q_n   
  ENDIF

  !! calculate DamOutflw
  IF ( c >= 0.5 ) THEN
    outflw = Krls * OutTmp 
  ELSE
    outflw = R * Krls * OutTmp + (1. - R) * inflw
  ENDIF
  ! the end of the flood control period (NDFC) is defined as the first month of the wet period preceding the start of the operational year
  ! The start of the flood control period (STFC) is defined as the month with the lowest flow within the dry period preceding the start of the operational period.
  ! the first month at which the long-term mean monthly flow falls below the long- term mean annual flow
END SUBROUTINE DAM_OPERATION_V13


SUBROUTINE DAM_OPERATION_LISFLOOD(Vol_dam, inflw, outflw, &
                                  Vol_con, Vol_nor, Vol_fld, Vol_tot, Q_n, Q_f)
  !! reference: https://ec-jrc.github.io/lisflood-model/3_03_optLISFLOOD_reservoirs/
  !! reference: https://github.com/ec-jrc/lisflood-code

  IMPLICIT NONE
  ! Input variables
  REAL(KIND=JPRB), INTENT(IN)                :: Vol_dam     ! Current volume of the dam (m3)
  REAL(KIND=JPRB), INTENT(IN)                :: inflw       ! Inflow rate to the dam (m3/s)
  REAL(KIND=JPRB), INTENT(OUT)               :: outflw      ! Outflow rate from the dam (m3/s)
  ! Parameters
  REAL(KIND=JPRB), INTENT(IN)  :: Vol_con          ! volume of water at the conservation level of the dam (m3)
  REAL(KIND=JPRB), INTENT(IN)  :: Vol_nor          ! volume of water at the normal level of the dam (m3)
  REAL(KIND=JPRB), INTENT(IN)  :: Vol_fld          ! volume of water at the flood controld level of the dam (m3)
  REAL(KIND=JPRB), INTENT(IN)  :: Vol_tot          ! total storage capacity (m3)
  REAL(KIND=JPRB), INTENT(IN)  :: Q_n              ! normal discharge (m3/s)
  REAL(KIND=JPRB), INTENT(IN)  :: Q_f              ! flood control discharge  (m3/s)

  ! Local variables
  REAL(KIND=JPRB)              :: F                ! fractional level of dam
  REAL(KIND=JPRB)              :: L_c              ! fractional level at conservation level
  REAL(KIND=JPRB)              :: L_n              ! fractional level at normal operating level
  REAL(KIND=JPRB)              :: L_f              ! fractional level at flood control level
  REAL(KIND=JPRB)              :: Q_min            ! minimum outflow rate
  REAL(KIND=JPRB)              :: Q_max            ! maximum outflow rate
  REAL(KIND=JPRB)              :: AdjLn            ! adjustment factor for L_n,in range 0.01 to 0.99
  REAL(KIND=JPRB)              :: AdjQnor          ! adjustment factor for Q_nor, in range 0.25 to 2
  REAL(KIND=JPRB)              :: Ladj_f           ! adjusted fractional level for maximum outflow
  REAL(KIND=JPRB)              :: Qadj_nor         ! adjusted normal outflow rate

  !=====================================================
  AdjLn = 0.5 ! adjustment factor for L_n,in range 0.01 to 0.99
  AdjQnor = 1 ! adjustment factor for Q_n, in range 0.25 to 2
  
  !*** Parameters calculation  
  F = Vol_dam / Vol_tot
  L_c = Vol_con / Vol_tot
  L_n = Vol_nor / Vol_tot
  L_f = Vol_fld / Vol_tot
  Ladj_f  = L_n + AdjLn * (L_f - L_n)

  Q_min = Q_n * 0.1
  Qadj_nor = MAX(Q_min, MIN(AdjQnor * Q_n, Q_f)) 
  
  !*** operation scheme
  IF (F <= 2* L_c) THEN
    outflw = MIN(Q_min, Vol_dam/86400.)
  
  ELSEIF (F > 2* L_c .and. F <= L_n) THEN
    outflw = Q_min + (Qadj_nor - Q_min) * ((F - 2*L_c)/(L_n - 2*L_c))

  ELSEIF (F > L_n .and. F <= Ladj_f) THEN
    outflw = Qadj_nor 

  ELSEIF (F > Ladj_f .and. F <= L_f) THEN
    outflw = Qadj_nor + (Q_f - Qadj_nor) * ((F - Ladj_f)/(L_f - Ladj_f))
  
  ELSEIF (F > L_f) THEN
    Q_max = MIN(Q_f, MAX(1.2 * inflw, Qadj_nor))   
    outflw = MAX(Q_max, (F - L_f - 0.01) * (Vol_tot/86400.))
  ENDIF
  
  ! the condition described below is applied in order to prevent outflow values that are too large compared to the inflow value.
  ! reference: https://github.com/ec-jrc/lisflood-code
  IF (F < L_f .and. outflw > MIN(1.2 * inflw, Qadj_nor)) THEN
    outflw = MIN(outflw, MAX(inflw, Qadj_nor))
  ENDIF
END SUBROUTINE DAM_OPERATION_LISFLOOD

SUBROUTINE DAM_OPERATION_H22(Vol_dam, inflw, outflw, &
                             Vol_con, Vol_fld, Vol_eme, Q_n, Q_f, k)
  !! reference: Hanazaki, R., Yamazaki, D., & Yoshimura, K. (2022). 
  !! Development of a reservoir flood control scheme for global flood models. 
  !! Journal of Advances in Modeling Earth Systems, 14(3), e2021MS002944.
  
  IMPLICIT NONE
  REAL(KIND=JPRB), INTENT(IN)                :: Vol_dam     ! Current volume of the dam (m3)
  REAL(KIND=JPRB), INTENT(IN)                :: inflw       ! Inflow rate to the dam (m3/s)
  REAL(KIND=JPRB), INTENT(OUT)               :: outflw      ! Outflow rate from the dam (m3/s)
  !*** parameter
  REAL(KIND=JPRB), INTENT(IN)                :: Vol_con     ! conservative storage (m3) 
  REAL(KIND=JPRB), INTENT(IN)                :: Vol_fld     ! flood control storage (m3) 
  REAL(KIND=JPRB), INTENT(IN)                :: Vol_eme     ! emergency storage (m3)   
  REAL(KIND=JPRB), INTENT(IN)                :: Q_n         ! normal discharge (m3/s)
  REAL(KIND=JPRB), INTENT(IN)                :: Q_f         ! flood discharge (m3/s)
  REAL(KIND=JPRB), INTENT(IN)                :: k           ! release coefficient
  !=====================================================

  !! case1: impoundment
  IF( Vol_dam <=  Vol_con)THEN
    outflw = Q_n * (Vol_dam / Vol_fld )

  !! case2: water supply
  ELSEIF( Vol_con < Vol_dam .and. Vol_dam <= Vol_fld )THEN
    IF( Q_f <= inflw )THEN
      outflw = Q_n * 0.5 +   (Vol_dam-Vol_con)/( Vol_fld-Vol_con)      * (Q_f - Q_n)
    ELSE
      outflw = Q_n * 0.5 + (((Vol_dam-Vol_con)/( Vol_eme-Vol_con))**2) * (Q_f - Q_n)
    ENDIF  

  !! case3: flood control
  ELSEIF( Vol_fld < Vol_dam .and. Vol_dam <= Vol_eme ) THEN
    IF( Q_f <= inflw ) THEN
      outflw = Q_f +  k * (Vol_dam-Vol_fld)/(Vol_eme-Vol_fld) * (inflw-Q_f)
    ELSE
      outflw = Q_n*0.5 + (((Vol_dam-Vol_con)/(Vol_eme-Vol_con))**2)* (Q_f - Q_n)
    ENDIF

  !! case4: emergency operation
  ELSE
    outflw = MAX(inflw, Q_f)
  ENDIF
END SUBROUTINE DAM_OPERATION_H22

SUBROUTINE DAM_OPERATION_YF(Vol_dam, inflw, outflw, &
                             Vol_con, Vol_eme, Vol_adj, Q_n, Q_f, Q_a)

!improved reservoir operation 'Yamazaki & Funato'
  REAL(KIND=JPRB), INTENT(IN)                :: Vol_dam     ! Current volume of the dam (m3)
  REAL(KIND=JPRB), INTENT(IN)                :: inflw       ! Inflow rate to the dam (m3/s)
  REAL(KIND=JPRB), INTENT(OUT)               :: outflw      ! Outflow rate from the dam (m3/s)
  !*** parameter
  REAL(KIND=JPRB), INTENT(IN)                :: Vol_con     ! conservative storage (m3) 
  REAL(KIND=JPRB), INTENT(IN)                :: Vol_eme     ! emergency storage (m3) 
  REAL(KIND=JPRB), INTENT(IN)                :: Vol_adj     ! stabilization storage (m3)
    
  REAL(KIND=JPRB), INTENT(IN)                :: Q_n         ! normal discharge (m3/s)
  REAL(KIND=JPRB), INTENT(IN)                :: Q_f         ! flood discharge (m3/s)
  REAL(KIND=JPRB), INTENT(IN)                :: Q_a         ! stabilization discharge (m3/s)

  REAL(KIND=JPRB)                            :: DamOutTmp
!=====================================================
  
!! Case 1: water use
  IF( Vol_dam<=Vol_con )THEN
    outflw = Q_n * (Vol_dam/Vol_con)**0.5
  !! case 2: water excess (just avobe ConVol, for outflow stability)
  ELSEIF( Vol_dam>Vol_con .and. Vol_dam<=Vol_adj ) THEN
    outflw = Q_n + ( (Vol_dam-Vol_con) / (Vol_adj-Vol_con) )**3.0 * (Q_a - Q_n)
  !! case 3: water excess
  ELSEIF( Vol_dam>Vol_adj .and. Vol_dam<=Vol_eme ) THEN
    !! (flood period)
    IF( inflw >= Q_f ) THEN
      !!figure left side No.2
      outflw = Q_n + ( (Vol_dam-Vol_con) / (Vol_eme-Vol_con) )**1.0 * (inflw- Q_n)
      DamOutTmp = Q_a + ( (Vol_dam-Vol_adj) / (Vol_eme-Vol_adj) )**0.1 * (Q_f - Q_a)
      outflw = max( outflw,DamOutTmp )
    !! (non-flood period)
    ELSE
      outflw = Q_a + ( (Vol_dam-Vol_adj) / (Vol_eme-Vol_adj) )**0.1 * (Q_f - Q_a)
    ENDIF
  !! case 4: emergency operation(no.1)
  ELSEIF( Vol_dam>Vol_eme )THEN
    !! (flood period)
    IF( inflw >= Q_f ) THEN
      outflw = inflw
    !! (non-flood period)
    ELSE
      outflw = Q_f
    ENDIF
  ENDIF
END SUBROUTINE DAM_OPERATION_YF

SUBROUTINE UPDATE_INFLOW
USE YOS_CMF_INPUT,      ONLY: PMINSLP, PMANFLD
USE YOS_CMF_MAP,        ONLY: D2RIVLEN, D2RIVMAN, D2ELEVTN, D2NXTDST, D2RIVWTH
USE YOS_CMF_PROG,       ONLY: D2RIVOUT_PRE, D2FLDOUT_PRE
USE YOS_CMF_DIAG,       ONLY: D2RIVDPH, D2RIVVEL, D2FLDDPH
IMPLICIT NONE
! SAVE for OpenMP
INTEGER(KIND=JPIM),SAVE    :: ISEQ, JSEQ
REAL(KIND=JPRB),SAVE       :: DSLOPE,DAREA,DVEL,DSLOPE_F,DARE_F,DVEL_F
!$OMP THREADPRIVATE     (JSEQ,DSLOPE,DAREA,DVEL,DSLOPE_F,DARE_F,DVEL_F)
!============================

!*** 1a. reset outflw & dam inflow
!$OMP PARALLEL DO SIMD
DO ISEQ=1, NSEQALL
  IF( I1DAM(ISEQ)>0 )THEN  !! if dam grid or upstream of dam, reset variables
    D2RIVOUT(ISEQ,1) = 0._JPRB
    D2FLDOUT(ISEQ,1) = 0._JPRB
    P2DAMINF(ISEQ,1) = 0._JPRD
  ENDIF
ENDDO
!$OMP END PARALLEL DO SIMD

!*** 1b. calculate dam inflow, using previous tstep discharge
#ifndef NoAtom_CMF
!$OMP PARALLEL DO  !! No OMP Atomic for bit-identical simulation (set in Mkinclude)
#endif
DO ISEQ=1, NSEQALL
  IF( I1DAM(ISEQ)==10 .or. I1DAM(ISEQ)==11 )THEN  !! if dam grid or upstream of dam
    JSEQ=I1NEXT(ISEQ)
#ifndef NoAtom_CMF
!$OMP ATOMIC
#endif
    P2DAMINF(JSEQ,1) = P2DAMINF(JSEQ,1) + D2RIVOUT_PRE(ISEQ,1) + D2FLDOUT_PRE(ISEQ,1) 
  ENDIF
END DO
#ifndef NoAtom_CMF
!$OMP END PARALLEL DO
#endif

!*** 1c. discharge for upstream grids of dams
!$OMP PARALLEL DO  !! No OMP Atomic for bit-identical simulation (set in Mkinclude)
DO ISEQ=1, NSEQRIV
  IF( I1DAM(ISEQ)==10 )THEN  !! if downstream is DAM
    JSEQ   = I1NEXT(ISEQ)
    ! === river flow
    DSLOPE = (D2ELEVTN(ISEQ,1)-D2ELEVTN(JSEQ,1)) * D2NXTDST(ISEQ,1)**(-1.)
    DSLOPE = max(DSLOPE,PMINSLP)

    DVEL   = D2RIVMAN(ISEQ,1)**(-1.) * DSLOPE**0.5 * D2RIVDPH(ISEQ,1)**(2./3.)
    DAREA  = D2RIVWTH(ISEQ,1) * D2RIVDPH(ISEQ,1)

    D2RIVVEL(ISEQ,1) = DVEL
    D2RIVOUT(ISEQ,1) = DAREA * DVEL
    D2RIVOUT(ISEQ,1) = MIN( D2RIVOUT(ISEQ,1), real(P2RIVSTO(ISEQ,1),JPRB)/DT )
    !=== floodplain flow
    DSLOPE_F = min( 0.005_JPRB,DSLOPE )    !! set min [instead of using weirequation for efficiency]
    DVEL_F   = PMANFLD**(-1.) * DSLOPE_F**0.5 * D2FLDDPH(ISEQ,1)**(2./3.)
    DARE_F   = REAL(P2FLDSTO(ISEQ,1),KIND=JPRB) * D2RIVLEN(ISEQ,1)**(-1.)
    DARE_F   = MAX( DARE_F - D2FLDDPH(ISEQ,1)*D2RIVWTH(ISEQ,1), 0._JPRB )   !!remove above river channel     area

    D2FLDOUT(ISEQ,1) = DARE_F * DVEL_F
    D2FLDOUT(ISEQ,1) = MIN(  D2FLDOUT(ISEQ,1), REAL(P2FLDSTO(ISEQ,1),KIND=JPRB)/DT )
  ENDIF
ENDDO
!$OMP END PARALLEL DO
END SUBROUTINE UPDATE_INFLOW

!==========================================================
SUBROUTINE MODIFY_OUTFLW
! modify outflow in order to avoid negative storage
IMPLICIT NONE
REAL(KIND=JPRD)            :: P2STOOUT(NSEQMAX,1)                      !! total outflow from a grid     [m3]
REAL(KIND=JPRD)            :: P2RIVINF(NSEQMAX,1)                      !! 
REAL(KIND=JPRD)            :: P2FLDINF(NSEQMAX,1)                      !! 
REAL(KIND=JPRB)            :: D2RATE(NSEQMAX,1)                        !! outflow correction
! SAVE for OpenMP
INTEGER(KIND=JPIM),SAVE    :: ISEQ, JSEQ
REAL(KIND=JPRB),SAVE       :: OUT_R1, OUT_R2, OUT_F1, OUT_F2, DIUP, DIDW
!$OMP THREADPRIVATE     (JSEQ,OUT_R1, OUT_R2, OUT_F1, OUT_F2, DIUP, DIDW)
!================================================
  
!*** 1. initialize & calculate P2STOOUT for normal cells
!$OMP PARALLEL DO
DO ISEQ=1, NSEQALL
  P2RIVINF(ISEQ,1) = 0._JPRD
  P2FLDINF(ISEQ,1) = 0._JPRD
  P2STOOUT(ISEQ,1) = 0._JPRD
  D2RATE(ISEQ,1) = 1._JPRB
ENDDO
!$OMP END PARALLEL DO

!! for normal cells ---------
#ifndef NoAtom_CMF
!$OMP PARALLEL DO
#endif
DO ISEQ=1, NSEQRIV                                                    !! for normalcells
  JSEQ=I1NEXT(ISEQ) ! next cell's pixel
  OUT_R1 = max(  D2RIVOUT(ISEQ,1),0._JPRB )
  OUT_R2 = max( -D2RIVOUT(ISEQ,1),0._JPRB )
  OUT_F1 = max(  D2FLDOUT(ISEQ,1),0._JPRB )
  OUT_F2 = max( -D2FLDOUT(ISEQ,1),0._JPRB )
  DIUP=(OUT_R1+OUT_F1)*DT
  DIDW=(OUT_R2+OUT_F2)*DT
#ifndef NoAtom_CMF
!$OMP ATOMIC
#endif
  P2STOOUT(ISEQ,1) = P2STOOUT(ISEQ,1) + DIUP 
#ifndef NoAtom_CMF
!$OMP ATOMIC
#endif
  P2STOOUT(JSEQ,1) = P2STOOUT(JSEQ,1) + DIDW 
END DO
#ifndef NoAtom_CMF
!$OMP END PARALLEL DO
#endif

!! for river mouth grids ------------
!$OMP PARALLEL DO
DO ISEQ=NSEQRIV+1, NSEQALL
  OUT_R1 = max( D2RIVOUT(ISEQ,1), 0._JPRB )
  OUT_F1 = max( D2FLDOUT(ISEQ,1), 0._JPRB )
  P2STOOUT(ISEQ,1) = P2STOOUT(ISEQ,1) + OUT_R1*DT + OUT_F1*DT
END DO
!$OMP END PARALLEL DO

!============================
!*** 2. modify outflow
!$OMP PARALLEL DO
DO ISEQ=1, NSEQALL
  IF ( P2STOOUT(ISEQ,1) > 1.E-8 ) THEN
    D2RATE(ISEQ,1) = MIN( (P2RIVSTO(ISEQ,1)+P2FLDSTO(ISEQ,1)) * P2STOOUT(ISEQ,1)**(-1.), 1._JPRD )
  ENDIF
ENDDO
!$OMP END PARALLEL DO

!! normal pixels------
#ifndef NoAtom_CMF
!$OMP PARALLEL DO  !! No OMP Atomic for bit-identical simulation (set in Mkinclude)
#endif
DO ISEQ=1, NSEQRIV ! for normal pixels
  JSEQ=I1NEXT(ISEQ)
  IF( D2RIVOUT(ISEQ,1) >= 0._JPRB )THEN
    D2RIVOUT(ISEQ,1) = D2RIVOUT(ISEQ,1)*D2RATE(ISEQ,1)
    D2FLDOUT(ISEQ,1) = D2FLDOUT(ISEQ,1)*D2RATE(ISEQ,1)
  ELSE
    D2RIVOUT(ISEQ,1) = D2RIVOUT(ISEQ,1)*D2RATE(JSEQ,1)
    D2FLDOUT(ISEQ,1) = D2FLDOUT(ISEQ,1)*D2RATE(JSEQ,1)
  ENDIF
#ifndef NoAtom_CMF
!$OMP ATOMIC
#endif
  P2RIVINF(JSEQ,1) = P2RIVINF(JSEQ,1) + D2RIVOUT(ISEQ,1)             !! total inflow to a grid (from upstream)
#ifndef NoAtom_CMF
!$OMP ATOMIC
#endif
  P2FLDINF(JSEQ,1) = P2FLDINF(JSEQ,1) + D2FLDOUT(ISEQ,1)
ENDDO
#ifndef NoAtom_CMF
!$OMP END PARALLEL DO
#endif

D2RIVINF(:,:)=P2RIVINF(:,:)  !! needed for SinglePrecisionMode
D2FLDINF(:,:)=P2FLDINF(:,:)

!! river mouth-----------------
!$OMP PARALLEL DO
DO ISEQ=NSEQRIV+1, NSEQALL
  D2RIVOUT(ISEQ,1) = D2RIVOUT(ISEQ,1)*D2RATE(ISEQ,1)
  D2FLDOUT(ISEQ,1) = D2FLDOUT(ISEQ,1)*D2RATE(ISEQ,1)
ENDDO
!$OMP END PARALLEL DO
END SUBROUTINE MODIFY_OUTFLW
!==========================================================
END SUBROUTINE CMF_DAMOUT_CALC_IRR

!####################################################################
SUBROUTINE CMF_DAMOUT_WRTE_IRR
! local
CHARACTER(len=36)          :: WriteTXT(NDAMX), WriteTXT2(NDAMX)
REAL(KIND=JPRB)            :: DDamInf, DDamOut
! File IO
INTEGER(KIND=JPIM),SAVE    :: ISEQD, JDAM
INTEGER(KIND=JPIM),SAVE    :: LOGDAM
CHARACTER(len=4),SAVE      :: cYYYY
CHARACTER(len=256),SAVE    :: CLEN, CFMT
CHARACTER(len=256),SAVE    :: DAMTXT
LOGICAL,SAVE               :: IsOpen
DATA IsOpen       /.FALSE./
! ==========================================

IF( LDAMTXT )THEN
  IF( .not. IsOpen )THEN
    IsOpen=.TRUE.
    WRITE(CYYYY,'(i4.4)') ISYYYY
    DAMTXT='./damtxt-'//trim(cYYYY)//'.txt'
    LOGDAM=INQUIRE_FID()
    OPEN(LOGDAM,FILE=DAMTXT,FORM='formatted')

    WRITE(CLEN,'(i0)') NDAMX
    CFMT="(i10,"//TRIM(CLEN)//"(a36))"

    JDAM=0
    DO IDAM=1, NDAM
      IF( DamSeq(IDAM)<=0 )CYCLE
      JDAM=JDAM+1
      ISEQD=DamSeq(IDAM)
      WRITE(WriteTxt(JDAM), '(i12)') DamID(IDAM)
    ENDDO

    WRITE(LOGDAM,CFMT) NDAMX, (WriteTXT(JDAM) ,JDAM=1, NDAMX)

    CFMT="(a10,"//TRIM(CLEN)//"(a36))"
    WRITE(LOGDAM,CFMT)  "Date", (WriteTXT2(JDAM),JDAM=1, NDAMX)
  ENDIF

  JDAM=0
  DO IDAM=1, NDAM
    IF( DamSeq(IDAM)<=0 ) CYCLE
    JDAM=JDAM+1
    ISEQD=DamSeq(IDAM)
    DDamInf=REAL( P2DAMINF(ISEQD,1),KIND=JPRB)
    DDamOut=D2RIVOUT(ISEQD,1) + D2FLDOUT(ISEQD,1)
    WRITE(WriteTxt(JDAM), '(3f12.2)') P2DAMSTO(ISEQD,1)*1.E-9, DDamInf, DDamOut
  ENDDO

  CFMT="(i10,"//TRIM(CLEN)//"(a36))"  
  WRITE(LOGDAM,CFMT) IYYYYMMDD, (WriteTXT(JDAM),JDAM=1, NDAMX)

ENDIF

END SUBROUTINE CMF_DAMOUT_WRTE_IRR

END MODULE CMF_CTRL_DAMOUT_MOD
