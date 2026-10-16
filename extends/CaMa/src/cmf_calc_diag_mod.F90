MODULE CMF_CALC_DIAG_MOD
!==========================================================
!* PURPOSE: Manage average and max diagnostic vars for output in CaMa-Flood
!
!* CONTAINS:
! -- CMF_DIAG_AVEMAX_OUTPUT   : Add / Max of diagnostic variables at time step
! -- CMF_DIAG_GETAVE_OUTPUT   : Calculate time-average of Diagnostic Variables
! -- CMF_DIAG_RESET_OUTPUT    : Reset Diagnostic Variables (Average & Maximum )
!!
! (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  Aug 2019
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
USE PARKIND1,           ONLY: JPIM, JPRM, JPRB
USE YOS_CMF_INPUT,      ONLY: LOGNAM
USE YOS_CMF_INPUT,      ONLY: DT, LPTHOUT,  LDAMOUT,  LWEVAP,LWINFILT,LSEDIMENT
USE YOS_CMF_MAP,        ONLY: NSEQMAX,      NPTHOUT,      NPTHLEV
USE YOS_CMF_PROG,       ONLY: D2RIVOUT,     D2FLDOUT,     D1PTHFLW,     D2GDWRTN, &
                            & D2RUNOFF,     D2ROFSUB,     P2DAMINF
USE YOS_CMF_DIAG,       ONLY: D2OUTFLW,     D2RIVVEL,     D2PTHOUT,     D2PTHINF, &
                            & D2RIVDPH,     D2STORGE,     D2WEVAPEX,   D2WINFILTEX
USE YOS_CMF_DIAG,       ONLY: D2RIVOUT_aAVG, D2FLDOUT_aAVG, D1PTHFLW_aAVG, D2GDWRTN_aAVG, D2RUNOFF_aAVG, D2ROFSUB_aAVG, &
                            & D2OUTFLW_aAVG, D2RIVVEL_aAVG, D2PTHOUT_aAVG, D2DAMINF_aAVG, D2WEVAPEX_aAVG,  D2WINFILTEX_aAVG,&
                            & D2OUTFLW_aMAX, D2RIVDPH_aMAX, D2STORGE_aMAX, NADD_adp,      D1PTHFLWSUM_aAVG
USE YOS_CMF_DIAG,       ONLY: D2RIVOUT_oAVG, D2FLDOUT_oAVG, D1PTHFLW_oAVG, D2GDWRTN_oAVG, D2RUNOFF_oAVG, D2ROFSUB_oAVG, &
                            & D2OUTFLW_oAVG, D2RIVVEL_oAVG, D2PTHOUT_oAVG, D2DAMINF_oAVG, D2WEVAPEX_oAVG,  D2WINFILTEX_oAVG, &
                            & D2OUTFLW_oMAX, D2RIVDPH_oMAX, D2STORGE_oMAX, NADD_out

USE YOS_CMF_INPUT,      ONLY: LSEDIMENT
USE CMF_CTRL_SED_MOD,   ONLY: nsed, d2sedout_avg, d2sedinp_avg, d2bedout_avg, d2netflw_avg
IMPLICIT NONE
CONTAINS 
!####################################################################
! -- CMF_DIAG_AVEMAX_ADPSTP   : Add / Max of diagnostic variables within adaptive time step
! -- CMF_DIAG_GETAVE_ADPSTP   : Calculate time-average of Diagnostic Variables for adaptive steps
! -- CMF_DIAG_RESET_ADPSTP    : Reset Diagnostic Variables (Average & Maximum ) for adaptive steps
!
! -- CMF_DIAG_AVEMAX_OUTPUT   : Add / Max of diagnostic variables at time step for output time step
! -- CMF_DIAG_GETAVE_OUTPUT   : Calculate time-average of Diagnostic Variables for output
! -- CMF_DIAG_RESET_OUTPUT    : Reset Diagnostic Variables (Average & Maximum ) for output
!
!####################################################################
SUBROUTINE CMF_DIAG_RESET_ADPSTP
USE YOS_CMF_TIME,       ONLY: JYYYYMMDD, JHHMM
IMPLICIT NONE
INTEGER(KIND=JPIM),SAVE  ::  ISEQ, IPTH
!================================================
WRITE(LOGNAM,*) "CMF::DIAG_AVERAGE: reset", JYYYYMMDD, JHHMM
NADD_adp=0

!$OMP PARALLEL DO SIMD
DO ISEQ=1, NSEQMAX
  D2RIVOUT_aAVG(ISEQ,1) = 0._JPRB
  D2FLDOUT_aAVG(ISEQ,1) = 0._JPRB
  D2OUTFLW_aAVG(ISEQ,1) = 0._JPRB
  D2RIVVEL_aAVG(ISEQ,1) = 0._JPRB
  D2PTHOUT_aAVG(ISEQ,1) = 0._JPRB
  D2GDWRTN_aAVG(ISEQ,1) = 0._JPRB
  D2RUNOFF_aAVG(ISEQ,1) = 0._JPRB
  D2ROFSUB_aAVG(ISEQ,1) = 0._JPRB

  D2STORGE_aMAX(ISEQ,1)=0._JPRB
  D2OUTFLW_aMAX(ISEQ,1)=0._JPRB
  D2RIVDPH_aMAX(ISEQ,1)=0._JPRB
ENDDO
!$OMP END PARALLEL DO SIMD

IF ( LDAMOUT ) THEN
  !$OMP PARALLEL DO SIMD
  DO ISEQ=1, NSEQMAX
    D2DAMINF_aAVG(ISEQ,1)  = 0._JPRB
  ENDDO
  !$OMP END PARALLEL DO SIMD
ENDIF
IF ( LWEVAP ) THEN
  !$OMP PARALLEL DO SIMD
  DO ISEQ=1, NSEQMAX
    D2WEVAPEX_aAVG(ISEQ,1) = 0._JPRB
  ENDDO
  !$OMP END PARALLEL DO SIMD
ENDIF

!$OMP PARALLEL DO SIMD
DO IPTH=1, NPTHOUT
  D1PTHFLW_aAVG(IPTH,:) = 0._JPRB 
  D1PTHFLWSUM_aAVG(IPTH)= 0._JPRB
ENDDO
!$OMP END PARALLEL DO SIMD

!reset sediment variables
!not here, need new variables

!IF (LSEDIMENT) THEN
!  !$OMP PARALLEL DO SIMD
!  DO ISEQ=1, NSEQMAX
!    DO ISED=1, nsed
!      d2sedout_avg(ISEQ,ISED) = 0._JPRB
!      d2sedinp_avg(ISEQ,ISED) = 0._JPRB
!      d2bedout_avg(ISEQ,ISED) = 0._JPRB
!      d2netflw_avg(ISEQ,ISED) = 0._JPRB
!    END DO
!  ENDDO
!  !$OMP END PARALLEL DO SIMD
!ENDIF


END SUBROUTINE CMF_DIAG_RESET_ADPSTP
!####################################################################
!
!
!
!####################################################################
SUBROUTINE CMF_DIAG_AVEMAX_ADPSTP
IMPLICIT NONE
INTEGER(KIND=JPIM),SAVE  ::  ISEQ, IPTH
!====================
NADD_adp=NADD_adp+DT
!$OMP PARALLEL DO SIMD
DO ISEQ=1, NSEQMAX
  D2RIVOUT_aAVG(ISEQ,1)=D2RIVOUT_aAVG(ISEQ,1)+D2RIVOUT(ISEQ,1)*DT
  D2FLDOUT_aAVG(ISEQ,1)=D2FLDOUT_aAVG(ISEQ,1)+D2FLDOUT(ISEQ,1)*DT
  D2RIVVEL_aAVG(ISEQ,1)=D2RIVVEL_aAVG(ISEQ,1)+D2RIVVEL(ISEQ,1)*DT
  D2OUTFLW_aAVG(ISEQ,1)=D2OUTFLW_aAVG(ISEQ,1)+D2OUTFLW(ISEQ,1)*DT

  D2PTHOUT_aAVG(ISEQ,1)=D2PTHOUT_aAVG(ISEQ,1)+D2PTHOUT(ISEQ,1)*DT-D2PTHINF(ISEQ,1)*DT

  D2GDWRTN_aAVG(ISEQ,1)=D2GDWRTN_aAVG(ISEQ,1)+D2GDWRTN(ISEQ,1)*DT
  D2RUNOFF_aAVG(ISEQ,1)=D2RUNOFF_aAVG(ISEQ,1)+D2RUNOFF(ISEQ,1)*DT
  D2ROFSUB_aAVG(ISEQ,1)=D2ROFSUB_aAVG(ISEQ,1)+D2ROFSUB(ISEQ,1)*DT

  D2OUTFLW_aMAX(ISEQ,1)=max( D2OUTFLW_aMAX(ISEQ,1), abs(D2OUTFLW(ISEQ,1)) )
  D2RIVDPH_aMAX(ISEQ,1)=max( D2RIVDPH_aMAX(ISEQ,1),     D2RIVDPH(ISEQ,1)  )
  D2STORGE_aMAX(ISEQ,1)=max( D2STORGE_aMAX(ISEQ,1),     D2STORGE(ISEQ,1)  )
         !recheck here zhongwang@Apr15,2024
         !IF( LWEVAP )THEN
         !   D2WEVAPEX_AVG(ISEQ,1)= (D2WEVAPEX_AVG(ISEQ,1) +D2WEVAPEX(ISEQ,1)*DT) /D2GRAREA(ISEQ,1)
         !ENDIF
         !IF( LWINFILT )THEN
         !   D2WINFILTEX_AVG(ISEQ,1)= (D2WINFILTEX_AVG(ISEQ,1) +D2WINFILTEX(ISEQ,1)*DT)/D2GRAREA(ISEQ,1)
         !ENDIF
ENDDO
!$OMP END PARALLEL DO SIMD

!! loop for optional variable (separated for computational efficiency)
IF( LDAMOUT )THEN
  !$OMP PARALLEL DO SIMD
  DO ISEQ=1, NSEQMAX
    D2DAMINF_aAVG(ISEQ,1)=D2DAMINF_aAVG(ISEQ,1)+REAL( P2DAMINF(ISEQ,1)*DT, KIND=JPRB)
  ENDDO
  !$OMP END PARALLEL DO SIMD
ENDIF
IF( LWEVAP )THEN
  !$OMP PARALLEL DO SIMD
  DO ISEQ=1, NSEQMAX
    D2WEVAPEX_aAVG(ISEQ,1)= D2WEVAPEX_aAVG(ISEQ,1) +D2WEVAPEX(ISEQ,1)*DT
  ENDDO
  !$OMP END PARALLEL DO SIMD
ENDIF

IF( LWINFILT )THEN
  !$OMP PARALLEL DO SIMD
  DO ISEQ=1, NSEQMAX
    D2WINFILTEX_aAVG(ISEQ,1)= D2WINFILTEX_aAVG(ISEQ,1) +D2WINFILTEX(ISEQ,1)*DT
  ENDDO
  !$OMP END PARALLEL DO SIMD
ENDIF


!$OMP PARALLEL DO SIMD
DO IPTH=1, NPTHOUT
  D1PTHFLW_aAVG(IPTH,:)=D1PTHFLW_aAVG(IPTH,:)+D1PTHFLW(IPTH,:)*DT
ENDDO
!$OMP END PARALLEL DO SIMD


!calculate average rivout and rivvel for sediment timestep
!IF( LSEDIMENT )THEN
!  sadd_riv = sadd_riv + DT
!  !$OMP PARALLEL DO SIMD
!  DO ISEQ=1, NSEQMAX
!    d2rivout_sed(ISEQ) = d2rivout_sed(ISEQ)+D2RIVOUT(ISEQ,1)*DT
!    d2rivvel_sed(ISEQ) = d2rivvel_sed(ISEQ)+D2RIVVEL(ISEQ,1)*DT
!  ENDDO
!  !$OMP END PARALLEL DO SIMD
!ENDIF

END SUBROUTINE CMF_DIAG_AVEMAX_ADPSTP
!####################################################################
!
!
!
!####################################################################
SUBROUTINE CMF_DIAG_GETAVE_ADPSTP
USE YOS_CMF_TIME,       ONLY: JYYYYMMDD, JHHMM
IMPLICIT NONE
INTEGER(KIND=JPIM),SAVE  ::  ISEQ, IPTH, ILEV
!================================================
WRITE(LOGNAM,*) "CMF::DIAG_AVERAGE: time-average", NADD_adp, JYYYYMMDD, JHHMM

!$OMP PARALLEL DO SIMD
DO ISEQ=1, NSEQMAX
  D2RIVOUT_aAVG(ISEQ,1) = D2RIVOUT_aAVG(ISEQ,1) / REAL(NADD_adp,KIND=JPRB)
  D2FLDOUT_aAVG(ISEQ,1) = D2FLDOUT_aAVG(ISEQ,1) / REAL(NADD_adp,KIND=JPRB)
  D2OUTFLW_aAVG(ISEQ,1) = D2OUTFLW_aAVG(ISEQ,1) / REAL(NADD_adp,KIND=JPRB)
  D2RIVVEL_aAVG(ISEQ,1) = D2RIVVEL_aAVG(ISEQ,1) / REAL(NADD_adp,KIND=JPRB)
  D2PTHOUT_aAVG(ISEQ,1) = D2PTHOUT_aAVG(ISEQ,1) / REAL(NADD_adp,KIND=JPRB)
  D2GDWRTN_aAVG(ISEQ,1) = D2GDWRTN_aAVG(ISEQ,1) / REAL(NADD_adp,KIND=JPRB)
  D2RUNOFF_aAVG(ISEQ,1) = D2RUNOFF_aAVG(ISEQ,1) / REAL(NADD_adp,KIND=JPRB)
  D2ROFSUB_aAVG(ISEQ,1) = D2ROFSUB_aAVG(ISEQ,1) / REAL(NADD_adp,KIND=JPRB)
ENDDO
!$OMP END PARALLEL DO SIMD

IF ( LDAMOUT ) THEN
!$OMP PARALLEL DO SIMD
  DO ISEQ=1, NSEQMAX
    D2DAMINF_aAVG(ISEQ,1)  = D2DAMINF_aAVG(ISEQ,1)  / REAL(NADD_adp,KIND=JPRB)
  ENDDO
!$OMP END PARALLEL DO SIMD
ENDIF

IF ( LWEVAP ) THEN
!$OMP PARALLEL DO SIMD
  DO ISEQ=1, NSEQMAX
    D2WEVAPEX_aAVG(ISEQ,1) = D2WEVAPEX_aAVG(ISEQ,1) / REAL(NADD_adp,KIND=JPRB)
  ENDDO
!$OMP END PARALLEL DO SIMD
ENDIF

IF ( LWINFILT ) THEN
!$OMP PARALLEL DO SIMD
  DO ISEQ=1, NSEQMAX
    D2WINFILTEX_oAVG(ISEQ,1) = D2WINFILTEX_aAVG(ISEQ,1) / REAL(NADD_adp,KIND=JPRB)
  ENDDO
!$OMP END PARALLEL DO SIMD
ENDIF

DO ILEV=1, NPTHLEV
!$OMP PARALLEL DO SIMD
  DO IPTH=1, NPTHOUT
    D1PTHFLW_aAVG(IPTH,ILEV) = D1PTHFLW_aAVG(IPTH,ILEV) / REAL(NADD_adp,KIND=JPRB)
    D1PTHFLWSUM_aAVG(IPTH)   = D1PTHFLWSUM_aAVG(IPTH) + D1PTHFLW_aAVG(IPTH,ILEV)  !! bifurcation height layer summation
  ENDDO
!$OMP END PARALLEL DO SIMD
ENDDO

END SUBROUTINE CMF_DIAG_GETAVE_ADPSTP
!####################################################################
!++
!++
!++
!++
!####################################################################
SUBROUTINE CMF_DIAG_RESET_OUTPUT
USE YOS_CMF_TIME,       ONLY: JYYYYMMDD, JHHMM
!USE CMF_CTRL_SED_MOD,        ONLY: d2sedout_avg, d2sedinp_avg, d2bedout_avg, d2netflw_avg, &
!                              d2sedv_avg, sadd_out, nsed, totlyrnum
IMPLICIT NONE
INTEGER(KIND=JPIM),SAVE  ::  ISEQ, IPTH
INTEGER(KIND=JPIM),SAVE  ::  ISED
!================================================
WRITE(LOGNAM,*) "CMF::DIAG_AVERAGE: reset", JYYYYMMDD, JHHMM
NADD_out=0
!$OMP PARALLEL DO SIMD
DO ISEQ=1, NSEQMAX
  D2RIVOUT_oAVG(ISEQ,1) = 0._JPRB
  D2FLDOUT_oAVG(ISEQ,1) = 0._JPRB
  D2OUTFLW_oAVG(ISEQ,1) = 0._JPRB
  D2RIVVEL_oAVG(ISEQ,1) = 0._JPRB
  D2PTHOUT_oAVG(ISEQ,1) = 0._JPRB
  D2GDWRTN_oAVG(ISEQ,1) = 0._JPRB
  D2RUNOFF_oAVG(ISEQ,1) = 0._JPRB
  D2ROFSUB_oAVG(ISEQ,1) = 0._JPRB

  D2STORGE_oMAX(ISEQ,1)=0._JPRB
  D2OUTFLW_oMAX(ISEQ,1)=0._JPRB
  D2RIVDPH_oMAX(ISEQ,1)=0._JPRB
ENDDO
!$OMP END PARALLEL DO SIMD

IF ( LDAMOUT ) THEN
!$OMP PARALLEL DO SIMD
  DO ISEQ=1, NSEQMAX
    D2DAMINF_oAVG(ISEQ,1)  = 0._JPRB
  ENDDO
!$OMP END PARALLEL DO SIMD
ENDIF
IF ( LWEVAP ) THEN
!$OMP PARALLEL DO SIMD
  DO ISEQ=1, NSEQMAX
    D2WEVAPEX_oAVG(ISEQ,1) = 0._JPRB
  ENDDO
!$OMP END PARALLEL DO SIMD
ENDIF

IF ( LWINFILT ) THEN
!$OMP PARALLEL DO SIMD
  DO ISEQ=1, NSEQMAX
    D2WINFILTEX_oAVG(ISEQ,1) = 0._JPRB
  ENDDO
!$OMP END PARALLEL DO SIMD
ENDIF


!$OMP PARALLEL DO SIMD
DO IPTH=1, NPTHOUT
  D1PTHFLW_oAVG(IPTH,:) = 0._JPRB
ENDDO
!$OMP END PARALLEL DO SIMD

! Reset sediment variables
IF (LSEDIMENT) THEN
  !$OMP PARALLEL DO SIMD
  DO ISEQ=1, NSEQMAX
    DO ISED=1, nsed
      d2sedout_avg(ISEQ,ISED) = 0._JPRB
      d2sedinp_avg(ISEQ,ISED) = 0._JPRB
      d2bedout_avg(ISEQ,ISED) = 0._JPRB
      d2netflw_avg(ISEQ,ISED) = 0._JPRB
    END DO
  ENDDO
  !$OMP END PARALLEL DO SIMD
ENDIF

END SUBROUTINE CMF_DIAG_RESET_OUTPUT
!####################################################################
!
!
!
!####################################################################
SUBROUTINE CMF_DIAG_AVEMAX_OUTPUT
!USE CMF_CTRL_SED_MOD,        ONLY: d2sedout_avg, d2sedinp_avg, d2bedout_avg, d2netflw_avg, &
!                              d2sedout, d2sedinp, d2bedout, d2netflw, d2sedcon, d2layer, &
!                              d2seddep, d2sedv_avg, sadd_out, sedDT, nsed, lsedflw, totlyrnum
IMPLICIT NONE
INTEGER(KIND=JPIM),SAVE  ::  ISEQ, IPTH
INTEGER(KIND=JPIM),SAVE  ::  ISED
!====================
NADD_out=NADD_out+DT


!$OMP PARALLEL DO SIMD
DO ISEQ=1, NSEQMAX
  D2RIVOUT_oAVG(ISEQ,1)=D2RIVOUT_oAVG(ISEQ,1)+D2RIVOUT_aAVG(ISEQ,1)*DT
  D2FLDOUT_oAVG(ISEQ,1)=D2FLDOUT_oAVG(ISEQ,1)+D2FLDOUT_aAVG(ISEQ,1)*DT
  D2RIVVEL_oAVG(ISEQ,1)=D2RIVVEL_oAVG(ISEQ,1)+D2RIVVEL_aAVG(ISEQ,1)*DT
  D2OUTFLW_oAVG(ISEQ,1)=D2OUTFLW_oAVG(ISEQ,1)+D2OUTFLW_aAVG(ISEQ,1)*DT

  D2PTHOUT_oAVG(ISEQ,1)=D2PTHOUT_oAVG(ISEQ,1)+D2PTHOUT_aAVG(ISEQ,1)*DT

  D2GDWRTN_oAVG(ISEQ,1)=D2GDWRTN_oAVG(ISEQ,1)+D2GDWRTN_aAVG(ISEQ,1)*DT
  D2RUNOFF_oAVG(ISEQ,1)=D2RUNOFF_oAVG(ISEQ,1)+D2RUNOFF_aAVG(ISEQ,1)*DT
  D2ROFSUB_oAVG(ISEQ,1)=D2ROFSUB_oAVG(ISEQ,1)+D2ROFSUB_aAVG(ISEQ,1)*DT

  D2OUTFLW_oMAX(ISEQ,1)=max( D2OUTFLW_oMAX(ISEQ,1), abs(D2OUTFLW_aMAX(ISEQ,1)) )
  D2RIVDPH_oMAX(ISEQ,1)=max( D2RIVDPH_oMAX(ISEQ,1),     D2RIVDPH_aMAX(ISEQ,1)  )
  D2STORGE_oMAX(ISEQ,1)=max( D2STORGE_oMAX(ISEQ,1),     D2STORGE_aMAX(ISEQ,1)  )

  IF( LWEVAP )THEN
    D2WEVAPEX_oAVG(ISEQ,1)= D2WEVAPEX_oAVG(ISEQ,1) +D2WEVAPEX_aAVG(ISEQ,1)*DT
  ENDIF

  IF( LWINFILT )THEN
    D2WINFILTEX_oAVG(ISEQ,1)= D2WINFILTEX_oAVG(ISEQ,1) +D2WINFILTEX_aAVG(ISEQ,1)*DT
  ENDIF
!need to add sediment variables oavg omax e.g.,
!  IF (LSEDIMENT) THEN
!    DO ISED=1, nsed
!      d2sedout_avg(ISEQ,ISED) = d2sedout_avg(ISEQ,ISED) + d2sedout(ISEQ,ISED) * DT
!      d2sedinp_avg(ISEQ,ISED) = d2sedinp_avg(ISEQ,ISED) + d2sedinp(ISEQ,ISED) * DT
!      d2bedout_avg(ISEQ,ISED) = d2bedout_avg(ISEQ,ISED) + d2bedout(ISEQ,ISED) * DT
!      d2netflw_avg(ISEQ,ISED) = d2netflw_avg(ISEQ,ISED) + d2netflw(ISEQ,ISED) * DT
!    END DO
!  ENDIF
ENDDO
!$OMP END PARALLEL DO SIMD

!! loop for optional variable (separated for computational efficiency)
IF( LDAMOUT )THEN
!$OMP PARALLEL DO SIMD
  DO ISEQ=1, NSEQMAX
    D2DAMINF_oAVG(ISEQ,1)=D2DAMINF_oAVG(ISEQ,1)+D2DAMINF_aAVG(ISEQ,1)*DT
  ENDDO
!$OMP END PARALLEL DO SIMD
ENDIF

!$OMP PARALLEL DO SIMD
DO IPTH=1, NPTHOUT
  D1PTHFLW_oAVG(IPTH,:)=D1PTHFLW_oAVG(IPTH,:)+D1PTHFLW_aAVG(IPTH,:)*DT
ENDDO
!$OMP END PARALLEL DO SIMD

END SUBROUTINE CMF_DIAG_AVEMAX_OUTPUT
!####################################################################
!
!
!
!####################################################################
SUBROUTINE CMF_DIAG_GETAVE_OUTPUT
USE YOS_CMF_TIME,       ONLY: JYYYYMMDD, JHHMM
!USE CMF_CTRL_SED_MOD,        ONLY: d2sedout_avg, d2sedcon, d2sedinp_avg, d2bedout_avg, d2layer, d2netflw_avg, &
!                             d2seddep, d2sedv_avg, nsed, totlyrnum
IMPLICIT NONE
INTEGER(KIND=JPIM),SAVE  ::  ISEQ, IPTH
INTEGER(KIND=JPIM),SAVE  ::  ISED
!================================================
WRITE(LOGNAM,*) "CMF::DIAG_AVERAGE: time-average", NADD_out, JYYYYMMDD, JHHMM
IF (LSEDIMENT) THEN
  WRITE(LOGNAM,*) "CMF::SED_DIAG_AVERAGE: sediment time-average", NADD_out
ENDIF

!$OMP PARALLEL DO SIMD
DO ISEQ=1, NSEQMAX
  D2RIVOUT_oAVG(ISEQ,1) = D2RIVOUT_oAVG(ISEQ,1) / REAL(NADD_out,KIND=JPRB)
  D2FLDOUT_oAVG(ISEQ,1) = D2FLDOUT_oAVG(ISEQ,1) / REAL(NADD_out,KIND=JPRB)
  D2OUTFLW_oAVG(ISEQ,1) = D2OUTFLW_oAVG(ISEQ,1) / REAL(NADD_out,KIND=JPRB)
  D2RIVVEL_oAVG(ISEQ,1) = D2RIVVEL_oAVG(ISEQ,1) / REAL(NADD_out,KIND=JPRB)
  D2PTHOUT_oAVG(ISEQ,1) = D2PTHOUT_oAVG(ISEQ,1) / REAL(NADD_out,KIND=JPRB)
  D2GDWRTN_oAVG(ISEQ,1) = D2GDWRTN_oAVG(ISEQ,1) / REAL(NADD_out,KIND=JPRB)
  D2RUNOFF_oAVG(ISEQ,1) = D2RUNOFF_oAVG(ISEQ,1) / REAL(NADD_out,KIND=JPRB)
  D2ROFSUB_oAVG(ISEQ,1) = D2ROFSUB_oAVG(ISEQ,1) / REAL(NADD_out,KIND=JPRB)
  
  IF ( LDAMOUT ) THEN
    D2DAMINF_oAVG(ISEQ,1)  = D2DAMINF_oAVG(ISEQ,1)  / REAL(NADD_out,KIND=JPRB)
  ENDIF
  IF ( LWEVAP ) THEN
    D2WEVAPEX_oAVG(ISEQ,1) = D2WEVAPEX_oAVG(ISEQ,1) / REAL(NADD_out,KIND=JPRB)
  ENDIF
  IF ( LWINFILT ) THEN
    D2WINFILTEX_oAVG(ISEQ,1) = D2WINFILTEX_oAVG(ISEQ,1) / REAL(NADD_out,KIND=JPRB)
  ENDIF
ENDDO
!$OMP END PARALLEL DO SIMD

!$OMP PARALLEL DO SIMD
DO IPTH=1, NPTHOUT
  D1PTHFLW_oAVG(IPTH,:) = D1PTHFLW_oAVG(IPTH,:) / REAL(NADD_out,KIND=JPRB)
ENDDO
!$OMP END PARALLEL DO SIMD

! Average sediment variables
IF (LSEDIMENT) THEN
  !$OMP PARALLEL DO SIMD
  DO ISEQ=1, NSEQMAX
    DO ISED=1, nsed
      ! Average sediment variables
      d2sedout_avg(ISEQ,ISED) = d2sedout_avg(ISEQ,ISED) / REAL(NADD_out,KIND=JPRB)
      d2sedinp_avg(ISEQ,ISED) = d2sedinp_avg(ISEQ,ISED) / REAL(NADD_out,KIND=JPRB)
      d2bedout_avg(ISEQ,ISED) = d2bedout_avg(ISEQ,ISED) / REAL(NADD_out,KIND=JPRB)
      d2netflw_avg(ISEQ,ISED) = d2netflw_avg(ISEQ,ISED) / REAL(NADD_out,KIND=JPRB)
    END DO
    ! Note: d2sedcon and d2layer are instantaneous values, not accumulated
    ! They don't need averaging like the flux variables
    ! d2seddep is also not accumulated in the current implementation 
    !!!maybe at CMF_DIAG_AVEMAX_OUTPUT? or ??? somewhere else?
END DO
 !$OMP END PARALLEL DO SIMD
ENDIF

END SUBROUTINE CMF_DIAG_GETAVE_OUTPUT
!####################################################################



END MODULE CMF_CALC_DIAG_MOD
