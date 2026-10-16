#include <define.h>

#ifdef LULCC
MODULE MOD_Lulcc_MassEnergyConserve

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: LulccMassEnergyConserve

CONTAINS


   SUBROUTINE LulccMassEnergyConserve

!-----------------------------------------------------------------------
!
! !DESCRIPTION
!  This is the main subroutine to execute the calculation of the restart
!  variables for the begin of next year.  There are mainly three ways to
!  adjust restart variables:
!
!  1) variable related to mass: area weighted mean of the source
!  patches, e.g., ldew, wliq_soisno.  variable related to energy: keep
!  energy conserve after the change of temperature, e.g., t_soisno.
!
!  2) recalculate according to physical process, e.g., dz_sno, scv,
!  fsno.
!
!  Created by Wanyi Lin and Hua Yuan, 07/2023
!
! !REVISIONS:
!
!  10/2023, Wanyi Lin: share the codes with REST_LulccTimeVariables(),
!           and simplify the codes in this subroutine.
!
!  01/2024, Wanyi Lin: use "enthalpy conservation" for snow layer
!           temperature calculation.
!
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_LandPatch
   USE MOD_LandElm
   USE MOD_Mesh
   USE MOD_SPMD_Task
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_TimeVariables
   USE MOD_Lulcc_Vars_TimeInvariants
   USE MOD_Lulcc_Vars_TimeVariables
   USE MOD_Lulcc_TransferTraceReadin
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_LandPFT
   USE MOD_Vars_PFTimeInvariants
   USE MOD_Vars_PFTimeVariables
#endif
#ifdef URBAN_MODEL
   USE MOD_LandUrban
   USE MOD_Urban_Vars_TimeVariables
   USE MOD_Urban_Vars_TimeInvariants
#endif
   USE MOD_Const_Physical, only: cpice, cpliq, hfus, tfrz, denh2o, denice
   USE MOD_GroundTemperature
   USE MOD_SnowFraction
   USE MOD_Albedo
   USE MOD_Namelist
   USE MOD_Hydro_SoilWater

   IMPLICIT NONE

!-------------------------- Local Variables ----------------------------
   integer, allocatable, dimension(:) :: grid_patch_s , grid_patch_e
   integer, allocatable, dimension(:) :: grid_patch_s_, grid_patch_e_
   integer, allocatable, dimension(:) :: locpxl
   integer :: numpxl,ipxl

   integer, allocatable :: frnp_(:)     !index of source patches
   integer, allocatable :: gu_(:)       !index of urban patches in last year's grid
   real(r8),allocatable :: cvsoil_(:,:) !heat capacity [J/(m2 K)]
   real(r8),allocatable :: h_(:,:)      !enthalpy of source patches [J/m2]

   integer :: k, ilc, num, inp_
   integer :: i, j, np, np_, selfnp_, l, ipft, ip, ip_, pc, pc_
   integer :: nsl                       !number of snow layer of the source patch with maximum area
   integer :: nsl_max                   !maximum number of snow layer considering all source patches
   integer :: u, u_, iu, selfu_, nurb, duclass
   integer :: nlc = N_land_classification
   integer :: ps, pe, ps_, pe_          !start and end index of patch pft

   real(r8), dimension(1:N_land_classification) :: lccpct_np
   real(r8):: sum_lccpct_np, wgt(maxsnl+1:nl_soil), hc(maxsnl+1:0)
   real(r8):: zi_sno(maxsnl+1:0)        !local variable for snow node and depth calculation
   real(r8):: vf_water                  !volumetric fraction liquid water within soil
   real(r8):: vf_ice                    !volumetric fraction ice lens within soil
   real(r8):: hcap                      !J/(m3 K)
   real(r8):: c_water                   !Specific heat of water * density of liquid water
   real(r8):: c_ice                     !Specific heat of ice   * density of ice
   real(r8):: denice_np(maxsnl+1:0), denh2o_np(maxsnl+1:0), rhosnow_np(maxsnl+1:0)
   real(r8):: wbef,wpre                 !water before and water present for water calculation heck
   ! real(r8):: fmelt                   !dimensionless melting factor
   real(r8):: wt                        !fraction of vegetation covered with snow [-]
   real(r8), parameter :: m = 1.0       !the value of m used in CLM4.5 is 1.0.
   ! real(r8) :: deltim = 1800.         !time step (seconds) TODO: be intent in
   logical :: FROM_SOIL, FOUND

   ! update for zwt, wa, ldew
   real(r8) :: tolerance, tol_z, tol_v, zi_soisno(0:nl_soil), sp_zi(0:nl_soil), sp_dz(1:nl_soil)
#ifdef Campbell_SOIL_MODEL
   integer, parameter :: nprms = 1
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
   integer, parameter :: nprms = 5
#endif
   real(r8) :: prms(nprms, 1:nl_soil)
   real(r8) :: ldew_tmp
!-----------------------------------------------------------------------

      IF (p_is_worker) THEN

         ! allocate with numelm
         allocate(grid_patch_s (numelm ))
         allocate(grid_patch_e (numelm ))
         allocate(grid_patch_s_(numelm_))
         allocate(grid_patch_e_(numelm_))

         grid_patch_e (:) = -1
         grid_patch_s (:) = -1
         grid_patch_e_(:) = -1
         grid_patch_s_(:) = -1

         ! loop for numelm of next year, patches at the beginning and end of
         ! the element were recorded landpatch%eindex is arranged in order,
         ! and the not land element is skipped so, if element is missing, the
         ! recorder is -1.
         DO i=1, numelm
            ! how many patches in ith element in this worker
            numpxl = count(landpatch%eindex==landelm%eindex(i))

            IF (allocated(locpxl)) deallocate(locpxl)
            allocate(locpxl(numpxl))

            ! get all patches' index that eindex is equal the i element
            locpxl = pack([(ipxl, ipxl=1, numpatch)], landpatch%eindex==landelm%eindex(i))
            ! the min index is the start of patch's index
            grid_patch_s(i) = minval(locpxl)
            ! the max index is the end of patch's index
            grid_patch_e(i) = maxval(locpxl)
         ENDDO

         ! same as above, loop for numelm of previous year
         ! patches at the beginning and end of the element were recorded
         DO i=1, numelm_
            numpxl = count(landpatch_%eindex==landelm_%eindex(i))

            IF (allocated(locpxl)) deallocate(locpxl)
            allocate(locpxl(numpxl))

            locpxl = pack([(ipxl, ipxl=1, numpatch_)], landpatch_%eindex==landelm_%eindex(i))

            grid_patch_s_(i) = minval(locpxl)
            grid_patch_e_(i) = maxval(locpxl)
         ENDDO

         DO i=1, numelm
            DO j=1,numelm_
               IF (landelm%eindex(i) == landelm_%eindex(j)) THEN
                  np = grid_patch_s (i)
                  np_= grid_patch_s_(j)

                  IF (np.le.0) CYCLE

                  DO WHILE (np.le.grid_patch_e(i))

IF (patchtype(np) .ne. 3) THEN !not a glacier patch

IF (DEF_USE_PFT .or. DEF_FAST_PC) THEN
                     lccpct_np(:) = 0
                     lccpct_np(1) = sum(lccpct_patches(np,1:), mask=patchtypes(:)==0)
                     lccpct_np(URBAN  )   = lccpct_patches(np,URBAN  )
                     lccpct_np(WETLAND)   = lccpct_patches(np,WETLAND)
                     lccpct_np(WATERBODY) = lccpct_patches(np,WATERBODY)
                     lccpct_np(GLACIERS)  = lccpct_patches(np,GLACIERS)
#ifndef CROP
                     !TODO: treat cropland separately for FAST_PC
                     IF ( DEF_FAST_PC ) THEN
                        lccpct_np(CROPLAND) = lccpct_patches(np,CROPLAND) &
                                            + lccpct_patches(np,14)
                        lccpct_np(1) = lccpct_np(1) - lccpct_np(CROPLAND)
                     ENDIF
#endif
ELSE
                     lccpct_np(:) = lccpct_patches(np,1:nlc)
ENDIF

                     num = count(lccpct_np .gt. 0)
                     sum_lccpct_np = sum(lccpct_np)
                     allocate ( frnp_  (                 num))
                     allocate ( cvsoil_(maxsnl+1:nl_soil,num))
                     allocate ( h_     (maxsnl+1:0,      num))

                     ! Source patch type which differs from np's type exists
                     IF ( (sum_lccpct_np - lccpct_np(patchclass(np))) .gt. 0 ) THEN

                        ! Get the index of source patches, and stored as frnp_
                        k = 0
                        DO ilc = 1, nlc

                           IF (lccpct_np(ilc) .gt. 0) THEN
                              FOUND = .FALSE.
                              k = k + 1
                              inp_ = np_

                              DO WHILE (inp_ .le. grid_patch_e_(j))

                                 ! Get the index of source patch that has the
                                 ! same LC, and stored as selfnp_
                                 IF (patchclass_(inp_) .eq. patchclass(np)) THEN
                                    selfnp_ = inp_
                                 ENDIF

                                 IF (patchclass_(inp_) .eq. ilc) THEN
                                    frnp_(k) = inp_
                                    FOUND = .TRUE.
                                    EXIT
                                 ENDIF
                                 inp_ = inp_ + 1
                              ENDDO
                              IF (.not.(FOUND)) THEN
                                 PRINT*, 'source patch not found, np', np, 'patchclass_', &
                                    patchclass_(grid_patch_s_(j):grid_patch_e_(j))
                              ENDIF
                           ELSE
                             CYCLE
                           ENDIF
                        ENDDO

                        ! Initialize
                        wliq_soisno (:,np)    = 0  !liquid water in layers [kg/m2]
                        wice_soisno (:,np)    = 0  !ice lens in layers [kg/m2]
                        t_soisno    (:,np)    = 0  !soil + snow layer temperature [K]
                        z_sno       (:,np)    = 0  !node depth [m]
                        dz_sno      (:,np)    = 0  !interface depth [m]
                        t_grnd        (np)    = 0  !ground surface temperature [K]
                        tleaf         (np)    = 0  !leaf temperature [K]
                        ldew          (np)    = 0  !depth of water on foliage [mm]
                        ldew_rain     (np)    = 0  !depth of rain on foliage [mm]
                        ldew_snow     (np)    = 0  !depth of snow on foliage [mm]
                        sag           (np)    = 0  !non dimensional snow age [-]
                        scv           (np)    = 0  !snow cover, water equivalent [mm]
                        snowdp        (np)    = 0  !snow depth [meter]
                        fsno          (np)    = 0  !frac of snow cover on ground
                        sigf          (np)    = 0  !frac of veg cover, exclud snow-covered veg [-]
                        zwt           (np)    = 0  !the depth to water table [m]
                        wa            (np)    = 0  !water storage in aquifer [mm]
                        wdsrf         (np)    = 0  !depth of surface water [mm]
                        smp         (:,np)    = 0  !soil matrix potential [mm]
                        hk          (:,np)    = 0  !hydraulic conductivity [mm h2o/s]

                        IF(DEF_USE_PLANTHYDRAULICS)THEN
                           vegwp    (:,np)    = 0  !vegetation water potential [mm]
                           gs0sun     (np)    = 0  !working copy of sunlit stomata conductance
                           gs0sha     (np)    = 0  !working copy of shaded stomata conductance
                        ENDIF

                        IF(DEF_USE_OZONESTRESS)THEN
                           lai_old    (np)    = 0  !lai in last time step
                        ENDIF

                        snw_rds     (:,np)    = 0  !effective grain radius (col,lyr) [microns, m-6]
                        mss_bcpho   (:,np)    = 0  !mass of hydrophobic BC in snow  (col,lyr) [kg]
                        mss_bcphi   (:,np)    = 0  !mass of hydrophillic BC in snow (col,lyr) [kg]
                        mss_ocpho   (:,np)    = 0  !mass of hydrophobic OC in snow  (col,lyr) [kg]
                        mss_ocphi   (:,np)    = 0  !mass of hydrophillic OC in snow (col,lyr) [kg]
                        mss_dst1    (:,np)    = 0  !mass of dust species 1 in snow  (col,lyr) [kg]
                        mss_dst2    (:,np)    = 0  !mass of dust species 2 in snow  (col,lyr) [kg]
                        mss_dst3    (:,np)    = 0  !mass of dust species 3 in snow  (col,lyr) [kg]
                        mss_dst4    (:,np)    = 0  !mass of dust species 4 in snow  (col,lyr) [kg]
                        ssno_lyr(:,:,:,np)    = 0  !snow layer absorption [-]

                        trad          (np)    = 0  !radiative temperature of surface [K]
                        tref          (np)    = 0  !2 m height air temperature [kelvin]
                        qref          (np)    = 0  !2 m height air specific humidity
                        rst           (np)    = 0  !canopy stomatal resistance (s/m)
                        emis          (np)    = 0  !averaged bulk surface emissivity
                        z0m           (np)    = 0  !effective roughness [m]
                        zol           (np)    = 0  !dimensionless height (z/L) used in M-O theory
                        rib           (np)    = 0  !bulk Richardson number in surface layer
                        ustar         (np)    = 0  !u* in similarity theory [m/s]
                        qstar         (np)    = 0  !q* in similarity theory [kg/kg]
                        tstar         (np)    = 0  !t* in similarity theory [K]
                        fm            (np)    = 0  !integral of profile function for momentum
                        fh            (np)    = 0  !integral of profile function for heat
                        fq            (np)    = 0  !integral of profile function for moisture


                        ! =============================================================
                        ! 1) Mass and Energy conserve adjustment (except for dz_sno).
                        ! =============================================================

                        ! Calculate the weight of temperature adjustment
                        c_water = cpliq * denh2o ! J/(m3 K) = 4188   [J/(kg K)]*1000(kg/m3)
                        c_ice   = cpice * denice ! J/(m3 K) = 2117.27[J/(kg K)]*917 (kg/m3)
                        cvsoil_(:,:) = 0
                        h_     (:,:) = 0
                        wgt(maxsnl+1:nl_soil) = 0
                        hc (maxsnl+1:0      ) = 0

                        DO k = 1, num

                           ! Soil ground and wetland heat capacity from: MOD_GroundTemperature.F90
                           DO l = 1, nl_soil
                              vf_water = wliq_soisno_(l,frnp_(k))/(dz_soi(l)*denh2o)
                              vf_ice   = wice_soisno_(l,frnp_(k))/(dz_soi(l)*denice)
                              hcap     = csol_(l,frnp_(k)) + vf_water*c_water + vf_ice*c_ice
                              cvsoil_(l,k) = hcap*dz_soi(l)
                           ENDDO

                           ! no snow layer exist
                           IF( dz_sno_(0,frnp_(k))<1.e-6 .and. scv_(frnp_(k))>0.) THEN
                              cvsoil_(1,k) = cvsoil_(1,k) + cpice*scv_(frnp_(k))
                           ENDIF

                           ! Snow heat capacity
                           IF( z_sno_(0,frnp_(k)) < 0 ) THEN
                              cvsoil_(:0,k) = cpliq*wliq_soisno_(:0,frnp_(k)) &
                                            + cpice*wice_soisno_(:0,frnp_(k))
                              h_(:0,k)      = (cpliq*wliq_soisno_(:0,frnp_(k)) + &
                                               cpice*wice_soisno_(:0,frnp_(k))) &
                                            * (t_soisno_(:0,frnp_(k)) - tfrz) &
                                            + hfus*wliq_soisno_(:0,frnp_(k))
                           ENDIF

                           wgt(maxsnl+1:nl_soil) = wgt(maxsnl+1:nl_soil) &
                              + cvsoil_(maxsnl+1:nl_soil,k) * lccpct_np(patchclass_(frnp_(k)))
                           hc(:0) = hc(:0) &
                                  + h_(:0,k) * lccpct_np(patchclass_(frnp_(k))) / sum_lccpct_np
                        ENDDO

                        ! Get the maximum lccpct for snow layers assignment
                        inp_ = frnp_(1)
                        k    = 2
                        DO WHILE (k .le. num)
                           IF ( lccpct_np(patchclass_(frnp_(k))) &
                                .gt. lccpct_np(patchclass_(inp_)) ) THEN
                              inp_ = frnp_(k)
                           ENDIF
                           k = k + 1
                        ENDDO

                        ! check if snow layer exist in patch inp_
                        nsl     = count(z_sno_(:,inp_) .lt. 0)
                        nsl_max = count(wgt(:0)        .gt. 0)
                        ! denh2o_np(maxsnl+1:0) = 0
                        ! denice_np(maxsnl+1:0) = 0
                        rhosnow_np(maxsnl+1:0) = 0 ! partial density of water/snow (ice + liquid)

                        IF (nsl > 0) THEN
                           ! move wgt above nsl to nsl
                           IF ( nsl_max > nsl) THEN
                              DO l = nsl+1, nsl_max
                                 wgt(-nsl+1) = wgt(-nsl+1) + wgt(-l+1)
                                 hc (-nsl+1) = hc (-nsl+1) + hc (-l+1)
                              ENDDO
                           ENDIF

                           DO k = 1, num

                              t_soisno (-nsl+1:0,np) = t_soisno (-nsl+1:0,np) &
                                 + t_soisno_(-nsl+1:0,frnp_(k))*cvsoil_(-nsl+1:0,k) &
                                 * lccpct_np(patchclass_(frnp_(k)))/wgt(-nsl+1:0)
                              wliq_soisno (-nsl+1:0,np) = wliq_soisno (-nsl+1:0,np) &
                                 + wliq_soisno_(-nsl+1:0,frnp_(k)) &
                                 * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                              wice_soisno (-nsl+1:0,np) = wice_soisno (-nsl+1:0,np) &
                                 + wice_soisno_(-nsl+1:0,frnp_(k)) &
                                 * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np

                              l = 1
                              DO WHILE ( (l .le. nsl) .and. (dz_sno_(-l+1,frnp_(k)) .gt. 0) )
                                 ! denh2o_np (-l+1) = denh2o_np(-l+1) &
                                 !    + wliq_soisno_(-l+1,frnp_(k))/dz_sno_(-l+1,frnp_(k)) &
                                 !    * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                                 ! denice_np (-l+1) = denice_np(-l+1) &
                                 !    + wice_soisno_(-l+1,frnp_(k))/dz_sno_(-l+1,frnp_(k)) &
                                 !    * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np

                                 ! Reference: MOD_GroundTemperature.F90 Line205
                                 rhosnow_np  (-l+1) = rhosnow_np(-l+1) &
                                    + (wliq_soisno_(-l+1,frnp_(k)) + wice_soisno_(-l+1,frnp_(k))) &
                                    / dz_sno_(-l+1,frnp_(k)) &
                                    * (lccpct_np(patchclass_(frnp_(k))) / sum_lccpct_np)
                                 l = l + 1
                                 IF (l .gt. -maxsnl) EXIT
                              ENDDO

                              ! if source patch has more snow layer than the main patch
                              IF (nsl .lt. -maxsnl) THEN
                                 l = nsl+1
                                 DO WHILE ( (l .le. -maxsnl) .and. (dz_sno_(-l+1,frnp_(k)) .gt. 0) )

                                    wliq_soisno(-nsl+1,np) = wliq_soisno (-nsl+1,np) &
                                       + wliq_soisno_(-l+1,frnp_(k)) &
                                       * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                                    wice_soisno(-nsl+1,np) = wice_soisno (-nsl+1,np) &
                                       + wice_soisno_(-l+1,frnp_(k)) &
                                       * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np

                                    t_soisno (-nsl+1,np) = t_soisno (-nsl+1,np) &
                                       + t_soisno_(-l+1,frnp_(k))*cvsoil_(-l+1,k) &
                                       * lccpct_np(patchclass_(frnp_(k)))/wgt(-nsl+1)

                                    ! denh2o_np (-nsl+1) = denh2o_np(-nsl+1) &
                                    !    + wliq_soisno_(-l+1,frnp_(k))/dz_sno_(-l+1,frnp_(k)) &
                                    !    * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                                    ! denice_np (-nsl+1) = denice_np(-nsl+1) &
                                    !    + wice_soisno_(-l+1,frnp_(k))/dz_sno_(-l+1,frnp_(k)) &
                                    !    * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                                    rhosnow_np  (-nsl+1) = rhosnow_np(-nsl+1) &
                                       + (wliq_soisno_(-l+1,frnp_(k))+wice_soisno_(-l+1,frnp_(k))) &
                                       / dz_sno_(-l+1,frnp_(k)) &
                                       * (lccpct_np(patchclass_(frnp_(k))) / sum_lccpct_np)
                                    l = l + 1
                                    IF (l .gt. -maxsnl) EXIT
                                 ENDDO
                              ENDIF
                           ENDDO



                           ! snow layer node and depth calculation according to new mass and density
                           zi_sno(0) = 0._r8
                           DO l = 0, -nsl+1, -1

                              ! IF (denice_np(l)>0 .and. denh2o_np(l)>0) THEN
                              !    dz_sno (l,np) = wice_soisno(l,np)/denice_np(l) &
                              !                  + wliq_soisno(l,np)/denh2o_np(l)

                              ! ELSEIF (denice_np(l)==0 .and. denh2o_np(l)>0) THEN
                              !    dz_sno (l,np) = wliq_soisno(l,np)/denh2o_np(l)
                              !    ! print*, 'denice=0! np=',np,'igbp=',patchclass(np),&
                              !              'nsl=',nsl,'frnp_=',frnp_
                              !    ! DO k = 1,num
                              !    !    print*,'frnp_=',frnp_(k),'wice=',wice_soisno(:0,frnp_(k))
                              !    ! ENDDO

                              ! ELSEIF (denh2o_np(l)==0 .and. denice_np(l)>0) THEN
                              !    dz_sno (l,np) = wice_soisno(l,np)/denice_np(l)
                              !    ! print*, 'denh2o=0! np=',np,'igbp=',patchclass(np),&
                              !              'nsl=',nsl,'frnp_=',frnp_
                              !    ! DO k = 1,num
                              !    !    print*,'frnp_=',frnp_(k),'wliq=',wliq_soisno(:0,frnp_(k))
                              !    ! ENDDO

                              ! ELSE
                              !    print*, 'denh2o and denice == 0! np=',np,'igbp=',patchclass(np),&
                              !            'nsl=',nsl,'frnp_=',frnp_
                              !    DO k = 1,num
                              !       print*,'frnp_=',frnp_(k),'wliq=',wliq_soisno(:0,frnp_(k))
                              !       print*,'frnp_=',frnp_(k),'wice=',wice_soisno(:0,frnp_(k))
                              !    ENDDO
                              !    CALL CoLM_stop()
                              ! ENDIF


                              ! Reference: MOD_SnowLayersCombineDivide.F90's subroutine combo
                              IF (hc(l) < 0.) THEN
                                 t_soisno (l,np) = tfrz + hc(l) / (cpice*wice_soisno(l,np) &
                                                 + cpliq*wliq_soisno(l,np))
                              ELSEIF (hc(l) .le. hfus*wliq_soisno(l,np)) THEN
                                 t_soisno (l,np) = tfrz
                              ELSE
                                 t_soisno (l,np) = tfrz + (hc(l) - hfus*wliq_soisno(l,np)) &
                                                 / (cpice*wice_soisno(l,np)+cpliq*wliq_soisno(l,np))
                              ENDIF


                              dz_sno (l,np) = (wice_soisno(l,np) + wliq_soisno(l,np))/rhosnow_np(l)

                              z_sno  (l,np) = zi_sno(l) - 0.5_r8*dz_sno(l,np)
                              IF (l-1 .lt. maxsnl+1) EXIT
                              zi_sno (l-1)  = zi_sno(l) - dz_sno(l,np)

                           ENDDO

                        ELSE
                           ! no snow layer exist in the main patch, add a layer
                           ! move wgt above soil to layer 0
                           IF ( nsl_max > nsl) THEN
                              DO l = nsl+1, nsl_max
                                 wgt(0) = wgt(0) + wgt(-l+1)
                              ENDDO
                           ENDIF

                           l = 0
                           DO WHILE (wgt(l) .gt. 0)
                              DO k = 1, num

                                 wliq_soisno(0,np) = wliq_soisno(0,np) &
                                    + wliq_soisno_(l,frnp_(k)) &
                                    * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                                 wice_soisno(0,np) = wice_soisno(0,np) &
                                    + wice_soisno_(l,frnp_(k)) &
                                    * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                                 t_soisno   (0,np) = t_soisno   (0,np) &
                                    + t_soisno_(l,frnp_(k))*cvsoil_(l,k) &
                                    * lccpct_np(patchclass_(frnp_(k)))/wgt(0)

                                 IF (dz_sno_(l,frnp_(k)) .gt. 0) THEN
                                    ! denh2o_np(0) = denh2o_np(0) &
                                    !    + wliq_soisno_(l,frnp_(k))/dz_sno_(l,frnp_(k)) &
                                    !    * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                                    ! denice_np(0) = denice_np(0) &
                                    !    + wice_soisno_(l,frnp_(k))/dz_sno_(l,frnp_(k)) &
                                    !    * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                                    rhosnow_np(0) = rhosnow_np(0) &
                                       + (wliq_soisno_(l,frnp_(k)) + wice_soisno_(l,frnp_(k))) &
                                       / dz_sno_(l,frnp_(k)) &
                                       * (lccpct_np(patchclass_(frnp_(k))) / sum_lccpct_np)
                                 ENDIF
                              ENDDO

                              l = l-1
                              IF (l .lt. maxsnl+1) EXIT
                           ENDDO

                           IF (wgt(0) .gt. 0) THEN

                              ! snow layer node and depth calculation according to
                              ! new mass and density
                              ! IF (denh2o_np(0)>0 .and. denh2o_np(0)>0) THEN
                              !    dz_sno (0,np) = wice_soisno(0,np)/denice_np(0) &
                              !                  + wliq_soisno(0,np)/denh2o_np(0)
                              ! ELSEIF (denice_np(0)==0 .and. denh2o_np(0)>0) THEN
                              !    dz_sno (0,np) = wliq_soisno(0,np)/denh2o_np(0)
                              !    ! print*, 'denice=0! stop! np=',np,'igbp=',patchclass(np),&
                              !              'nsl=',nsl,'frnp_=',frnp_
                              !    ! DO k = 1,num
                              !    !    print*,'frnp_=',frnp_(k),'wice=',wice_soisno(:0,frnp_(k))
                              !    ! ENDDO
                              ! ELSEIF (denice_np(0)>0 .and. denh2o_np(0)==0) THEN
                              !    dz_sno (0,np) = wice_soisno(0,np)/denice_np(0)
                              !    ! print*, 'denh2o=0! stop! np=',np,'igbp=',patchclass(np),&
                              !              'nsl=',nsl,'frnp_=',frnp_
                              !    ! DO k = 1,num
                              !    !    print*,'frnp_=',frnp_(k),'wliq=',wliq_soisno(:0,frnp_(k))
                              !    ! ENDDO
                              ! ELSE
                              !    print*, 'denh2o and denice == 0! np=',np,'igbp=',patchclass(np),&
                              !            'nsl=',nsl,'frnp_=',frnp_
                              !    DO k = 1,num
                              !       print*,'frnp_=',frnp_(k),'wliq=',wliq_soisno(:0,frnp_(k))
                              !       print*,'frnp_=',frnp_(k),'wice=',wice_soisno(:0,frnp_(k))
                              !    ENDDO
                              !    CALL CoLM_stop()
                              ! ENDIF


                              ! Reference: MOD_SnowLayersCombineDivide.F90's subroutine combo
                              IF (hc(0) < 0.) THEN
                                 t_soisno (0,np) = tfrz &
                                    + hc(0) / (cpice*wice_soisno(0,np) + cpliq*wliq_soisno(0,np))
                              ELSEIF (hc(0) .le. hfus*wliq_soisno(0,np)) THEN
                                 t_soisno (0,np) = tfrz
                              ELSE
                                 t_soisno (0,np) = tfrz + (hc(0) - hfus*wliq_soisno(0,np)) &
                                                 / (cpice*wice_soisno(0,np)+cpliq*wliq_soisno(0,np))
                              ENDIF

                              dz_sno (0,np) = (wice_soisno(0,np) + wliq_soisno(0,np))/rhosnow_np(0)

                              zi_sno (0)    = 0._r8
                              z_sno  (0,np) = zi_sno(0) - 0.5_r8*dz_sno(0,np)
                           ENDIF

                        ENDIF


                        ! Variable adjustment
                        DO k = 1, num

                           wliq_soisno (1:nl_soil,np) = wliq_soisno (1:nl_soil,np) &
                              + wliq_soisno_(1:nl_soil,frnp_(k)) &
                              * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           wice_soisno (1:nl_soil,np) = wice_soisno (1:nl_soil,np) &
                              + wice_soisno_(1:nl_soil,frnp_(k)) &
                              * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           t_soisno (1:nl_soil,np) = t_soisno (1:nl_soil,np) &
                              + t_soisno_(1:nl_soil,frnp_(k))*cvsoil_(1:nl_soil,k) &
                              * lccpct_np(patchclass_(frnp_(k)))/wgt(1:nl_soil)

                           tleaf     (np) = tleaf     (np) + tleaf_     (frnp_(k)) &
                                          * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           ldew      (np) = ldew      (np) + ldew_      (frnp_(k)) &
                                          * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           ldew_rain (np) = ldew_rain (np) + ldew_rain_ (frnp_(k)) &
                                          * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           ldew_snow (np) = ldew_snow (np) + ldew_snow_ (frnp_(k)) &
                                          * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           sag       (np) = sag       (np) + sag_       (frnp_(k)) &
                                          * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np

                           ! TODO: use MOD_SnowFraction.F90 to calculate sigf later - DONE
                           ! sigf    (np) = sigf      (np) + sigf_      (frnp_(k)) &
                           !              * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           wa        (np) = wa        (np) + wa_        (frnp_(k)) &
                                          * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           wdsrf     (np) = wdsrf     (np) + wdsrf_     (frnp_(k)) &
                                          * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           IF (.not. DEF_USE_VariablySaturatedFlow) THEN
                              zwt    (np) = zwt       (np) + zwt_       (frnp_(k)) &
                                          * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           ENDIF

                           snw_rds   (:,np) = snw_rds   (:,np) + snw_rds_   (:,frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           mss_bcpho (:,np) = mss_bcpho (:,np) + mss_bcpho_ (:,frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           mss_bcphi (:,np) = mss_bcphi (:,np) + mss_bcphi_ (:,frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           mss_ocpho (:,np) = mss_ocpho (:,np) + mss_ocpho_ (:,frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           mss_ocphi (:,np) = mss_ocphi (:,np) + mss_ocphi_ (:,frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           mss_dst1  (:,np) = mss_dst1  (:,np) + mss_dst1_  (:,frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           mss_dst2  (:,np) = mss_dst2  (:,np) + mss_dst2_  (:,frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           mss_dst3  (:,np) = mss_dst3  (:,np) + mss_dst3_  (:,frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           mss_dst4  (:,np) = mss_dst4  (:,np) + mss_dst4_  (:,frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           ssno_lyr  (:,:,:,np) = ssno_lyr (:,:,:,np) + ssno_lyr_ (:,:,:,frnp_(k)) &
                                                * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np

                           ! TODO: or use same type assignment
                           smp       (:,np) = smp    (:,np) + smp_   (:,frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           hk        (:,np) = hk     (:,np) + hk_    (:,frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np

                           IF(DEF_USE_PLANTHYDRAULICS)THEN
                              vegwp  (:,np) = vegwp  (:,np) + vegwp_ (:,frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                              gs0sun   (np) = gs0sun   (np) + gs0sun_  (frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                              gs0sha   (np) = gs0sha   (np) + gs0sha_  (frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           ENDIF

                           IF(DEF_USE_OZONESTRESS)THEN
                              lai_old  (np) = lai_old  (np) + lai_old_ (frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           ENDIF

                           trad        (np) = trad     (np) + trad_    (frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           tref        (np) = tref     (np) + tref_    (frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           qref        (np) = qref     (np) + qref_    (frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           rst         (np) = rst      (np) + rst_     (frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           emis        (np) = emis     (np) + emis_    (frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           z0m         (np) = z0m      (np) + z0m_     (frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           zol         (np) = zol      (np) + zol_     (frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           rib         (np) = rib      (np) + rib_     (frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           ustar       (np) = ustar    (np) + ustar_   (frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           qstar       (np) = qstar    (np) + qstar_   (frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           tstar       (np) = tstar    (np) + tstar_   (frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           fm          (np) = fm       (np) + fm_      (frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           fh          (np) = fh       (np) + fh_      (frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           fq          (np) = fq       (np) + fq_      (frnp_(k)) &
                                            * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                        ENDDO

                        ! water balance check
                        wbef = 0
                        wpre = 0
                        DO k = 1, num
                           wbef = wbef + ldew_(frnp_(k)) &
                                       * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           wbef = wbef + scv_ (frnp_(k)) &
                                       * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           wbef = wbef + wa_  (frnp_(k)) &
                                       * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           wbef = wbef + sum(wliq_soisno_(maxsnl+1:nl_soil,frnp_(k))) &
                                       * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           wbef = wbef + sum(wice_soisno_(maxsnl+1:nl_soil,frnp_(k))) &
                                       * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                        ENDDO

                        wpre = ldew(np) + scv(np) + wa(np) + sum(wliq_soisno(maxsnl+1:nl_soil,np)) &
                             + sum(wice_soisno(maxsnl+1:nl_soil,np))
                        IF (wpre-wbef > 1.e-6) THEN
                           print*,'np=',np,'total err=',wpre-wbef
                        ENDIF

                        wbef = 0
                        wpre = 0
                        DO k = 1, num
                           wbef = wbef + sum(wliq_soisno_(maxsnl+1:nl_soil,frnp_(k))) &
                                       * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           wbef = wbef + sum(wice_soisno_(maxsnl+1:nl_soil,frnp_(k))) &
                                       * lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                        ENDDO

                        wpre = sum(wliq_soisno(maxsnl+1:nl_soil,np)) &
                             + sum(wice_soisno(maxsnl+1:nl_soil,np))
                        IF (wpre-wbef > 1.e-6) THEN
                           print*,'np=',np,'wice+wliq err=',wpre-wbef
                        ENDIF


                        ! =============================================================
                        ! 2) adjusted based on code of physical process.
                        ! =============================================================

                        DO l = maxsnl+1, 0
                           IF ( z_sno(l,np) .lt. 0 ) THEN
                              scv(np) = scv(np) + wice_soisno(l,np) + wliq_soisno(l,np)
                              snowdp(np) = snowdp(np) + dz_sno(l,np)
                           ENDIF
                        ENDDO

                        ! Use restart value from the same type of
                        ! source patch or remain initialized
                        ! IF (lccpct_np(patchclass(np)) .gt. 0) THEN
                        !    tleaf         (np) = tleaf_         (selfnp_)
                        !    lake_icefrac(:,np) = lake_icefrac_(:,selfnp_)
                        !    t_lake      (:,np) = t_lake_      (:,selfnp_)
                        ! ENDIF

                        ! Fraction of soil covered by snow
                        zlnd     = 0.01    !Roughness length for soil [m]
                        ! fsno(np) = 0.0
                        ! IF (snowdp(np) > 0.) THEN
                        !    fmelt = (scv(np)/snowdp(np)/100.) ** m
                        !    fsno(np)  = tanh(snowdp(np)/(2.5 * zlnd * fmelt))
                        ! ENDIF

                        ! Sigf, fsno
                        CALL snowfraction (tlai(np),tsai(np),z0m(np),zlnd,&
                                           scv(np),snowdp(np),wt,sigf(np),fsno(np))
                        sai(np) = tsai(np) * sigf(np)

                        ! account for vegetation snow
                        IF ( DEF_VEG_SNOW ) THEN
                           lai(np) = tlai(np) * sigf(np)
                        ENDIF

                        ! ! In case lai+sai come into existence this year, set sigf to 1;
                        ! Update: won't happen if CALL snowfraction
                        ! IF ( (lai(np) + sai(np)).gt.0 .and. sigf(np).eq.0 ) THEN
                        !    sigf(np) = 1
                        ! ENDIF

                        ! Set Groud temperature
                        IF ( sum( z_sno(:,np) ) .eq. 0 )  THEN
                           t_grnd(np) = t_soisno(1,np)
                        ELSE
                           DO k = maxsnl+1, 0
                              IF ( z_sno(k,np) .lt. 0 ) THEN
                                 t_grnd(np) = t_soisno(k,np)
                                 EXIT
                              ENDIF
                           ENDDO
                        ENDIF


                        ! Zwt calculation suggested by Shupeng Zhang
                        IF (DEF_USE_VariablySaturatedFlow) THEN
                           IF (wa(np) >= 0) THEN
                              DO l = nl_soil, 1, -1
                                 vf_water = (wliq_soisno(l,np)/denh2o) / dz_soi(l)
                                 vf_ice   = (wice_soisno(l,np)/denice) / dz_soi(l)
                                 IF ((vf_water + vf_ice) < porsl(l,np)) THEN ! fraction of soil that is voids [-]
                                    zwt(np) = zi_soi(l) ! interface depths (lower bound) of soil layers
                                    EXIT
                                 ELSEIF (l==1) THEN
                                    zwt(np) = 0
                                 ENDIF
                              ENDDO
                           ELSE
                              PRINT*, 'Water table is below the soil column.'
#ifdef Campbell_SOIL_MODEL
                              prms(1,1:nl_soil) = bsw(1:nl_soil,np)
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
                              prms(1,1:nl_soil) = alpha_vgm(1:nl_soil,np)
                              prms(2,1:nl_soil) = n_vgm    (1:nl_soil,np)
                              prms(3,1:nl_soil) = L_vgm    (1:nl_soil,np)
                              prms(4,1:nl_soil) = sc_vgm   (1:nl_soil,np)
                              prms(5,1:nl_soil) = fc_vgm   (1:nl_soil,np)
#endif

                              tolerance = 1.e-3
                              tol_z = tolerance / real(nl_soil,r8) /2.0
                              zi_soisno = (/0., zi_soi/)
                              sp_zi(0:nl_soil) = zi_soisno(0:nl_soil) * 1000.0   ! from meter to mm
                              sp_dz(1:nl_soil) = sp_zi(1:nl_soil) - sp_zi(0:nl_soil-1)
                              tol_v = tol_z / maxval(sp_dz)
                              CALL get_zwt_from_wa(porsl(nl_soil,np), theta_r(nl_soil,np), psi0(nl_soil,np), hksati(nl_soil,np), &
                                                   nprms, prms(:,nl_soil), tol_v, tol_z, wa(np), sp_zi(nl_soil), zwt(np)         )
                           ENDIF
                        ENDIF


                     ! ELSE
                     ! ! Patch area stay unchanged or decrease,
                     !   use restart value or remain initialized
                     ! ! TODO: CALL REST - DONE
                     !    inp_ = np_
                     !    DO WHILE (inp_ .le. grid_patch_e_(j))
                     !       IF (patchclass_(inp_) .eq. patchclass(np)) THEN
                     !          selfnp_            = inp_
                     !          frnp_(1)           = inp_
                     !          wliq_soisno (:,np) = wliq_soisno_ (:,inp_)
                     !          wice_soisno (:,np) = wice_soisno_ (:,inp_)
                     !          t_soisno    (:,np) = t_soisno_    (:,inp_)
                     !          z_sno       (:,np) = z_sno_       (:,inp_)
                     !          dz_sno      (:,np) = dz_sno_      (:,inp_)
                     !          t_grnd        (np) = t_grnd_        (inp_)
                     !          tleaf         (np) = tleaf_         (inp_)
                     !          ldew          (np) = ldew_          (inp_)
                     !          sag           (np) = sag_           (inp_)
                     !          scv           (np) = scv_           (inp_)
                     !          snowdp        (np) = snowdp_        (inp_)
                     !          fsno          (np) = fsno_          (inp_)
                     !          sigf          (np) = sigf_          (inp_)
                     !          ! In case lai+sai come into existence this year, set sigf to 1
                     !          IF ( (lai(np) + sai(np)).gt.0 .and. sigf(np).eq.0 ) THEN
                     !             sigf(np) = 1
                     !          ENDIF
                     !          zwt           (np) = zwt_           (inp_)
                     !          wa            (np) = wa_            (inp_)
                     !          EXIT
                     !       ENDIF
                     !
                     !       inp_ = inp_ + 1
                     !    ENDDO

                     ENDIF

! ELSEIF (patchtype(np)==3) THEN !glacier patch
!                    ! Used restart value for GLACIERS patches if patchclass exists last year,
!                    ! or remain initialized
!                    ! TODO: CALL REST - DONE
!                    inp_ = np_
!                    DO WHILE (inp_ .le. grid_patch_e_(j))
!                       IF (patchclass_(inp_) .eq. patchclass(np)) THEN
!                          wliq_soisno (:,np) = wliq_soisno_ (:,inp_)
!                          wice_soisno (:,np) = wice_soisno_ (:,inp_)
!                          t_soisno    (:,np) = t_soisno_    (:,inp_)
!                          z_sno       (:,np) = z_sno_       (:,inp_)
!                          dz_sno      (:,np) = dz_sno_      (:,inp_)
!                          t_grnd        (np) = t_grnd_        (inp_)
!                          tleaf         (np) = tleaf_         (inp_)
!                          ldew          (np) = ldew_          (inp_)
!                          sag           (np) = sag_           (inp_)
!                          scv           (np) = scv_           (inp_)
!                          snowdp        (np) = snowdp_        (inp_)
!                          fsno          (np) = fsno_          (inp_)
!                          sigf          (np) = sigf_          (inp_)
!                          zwt           (np) = zwt_           (inp_)
!                          wa            (np) = wa_            (inp_)
!                          EXIT
!                       ENDIF
!                       inp_ = inp_ + 1
!                    ENDDO
ENDIF

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)

                     IF (patchtype(np)==0) THEN
                        ps  = patch_pft_s(np)
                        pe  = patch_pft_e(np)
                        ! ps_ = patch_pft_s_(selfnp_)
                        ! pe_ = patch_pft_e_(selfnp_)

                        ! if totally come from other types,ldew set to zero since ldew_p(:)=0
                        ! using MEC replace STA
                        ldew_tmp = sum( ldew_p(ps:pe)*pftfrac(ps:pe) )
                        IF ((ldew_tmp) > 0) THEN
                           ldew_p(ps:pe) = ldew_p(ps:pe) * (ldew(np) / ldew_tmp)
                        ELSE
                           ldew(np) = 0
                        ENDIF

                        ! z0m_p was same-type assigned, then here we update sigf_p, sigf, fsno
                        CALL snowfraction_pftwrap (np,zlnd,scv(np),snowdp(np),wt,sigf(np),fsno(np))

                        sai_p(ps:pe) = tsai_p(ps:pe) * sigf_p(ps:pe)
                        sai(np) = sum(sai_p(ps:pe)*pftfrac(ps:pe))

                        ! account for vegetation snow
                        IF ( DEF_VEG_SNOW ) THEN
                           lai_p(ps:pe) = tlai_p(ps:pe) * sigf_p(ps:pe)
                           lai(np) = sum(lai_p(ps:pe)*pftfrac(ps:pe))
                        ENDIF

                     ENDIF

                     ! ! TODO: CALL REST - DONE
                     ! IF (patchtype(np)==0 .and. lccpct_np(patchclass(np)) .gt. 0) THEN
                     !    ! Used restart value of the same pftclass for pft-specific variables
                     !    ! Note: For ip-specific variables, remain initialized value
                     !      for new soil patch or pftclass
                     !    ip = ps
                     !    ip_= ps_
                     !
                     !    IF (ip.le.0 .or. ip_.le.0) THEN
                     !       print *, "Error in LuLccMassEnergyConserve LULC_IGBP_PFT|LULC_IGBP_PC!"
                     !       STOP
                     !    ENDIF
                     !
                     !    DO WHILE (ip.le.pe .and. ip_.le.pe_)
                     !       ! if a PFT is missing, CYCLE
                     !       IF (pftclass(ip) > pftclass_(ip_)) THEN
                     !         ip_= ip_+ 1
                     !         CYCLE
                     !       ENDIF
                     !
                     !       ! if a PFT is added, CYCLE
                     !       IF (pftclass(ip) < pftclass_(ip_)) THEN
                     !         ip = ip + 1
                     !         CYCLE
                     !       ENDIF
                     !
                     !       ! for the same PFT, set PFT value
                     !       tleaf_p    (ip) = tleaf_p_    (ip_)
                     !       ldew_p     (ip) = ldew_p_     (ip_)
                     !       ! use MOD_SnowFraction.F90 later
                     !       sigf_p     (ip) = sigf_p_     (ip_)
                     !       IF ( (lai_p(ip) + sai_p(ip)).gt.0 .and. sigf_p(ip).eq.0 ) THEN
                     !          sigf_p(ip) = 1
                     !       ENDIF
                     !
                     !       ip = ip + 1
                     !       ip_= ip_+ 1
                     !    ENDDO
                     !    ldew(np) = sum( ldew_p(ps:pe)*pftfrac(ps:pe) )
                     ! ENDIF

#endif

#ifdef URBAN_MODEL
                     IF (patchclass(np)==URBAN) THEN

                        ! If there isn't any urban patch in last year's grid,initialized value was
                        ! remained.  Though the first source soil patch would be used for pervious
                        ! ground related variables.
                        u = patch2urban (np)
                        nurb = count( patchclass_(grid_patch_s_(j):grid_patch_e_(j)) == URBAN )

                        ! Get the index of urban patches in last year's grid, and index of urban
                        ! patch with the same urbclass
                        IF (nurb > 0) THEN

                           allocate(gu_(nurb)) ! ind of urban patches in last year's grid
                           selfu_   = -1       ! ind of urban patch with same class in last year
                           inp_     = np_      ! for loop to record the index of urban patch
                           iu       = 0

                           DO WHILE (inp_ .le. grid_patch_e_(j))
                              IF (patchclass_(inp_) == URBAN) THEN
                                 iu = iu + 1
                                 gu_(iu) = patch2urban_(inp_)
                                 IF (landurban%settyp(u) == urbclass_(gu_(iu))) THEN
                                    selfu_ = gu_(iu)
                                 ENDIF
                              ENDIF
                              inp_ = inp_ + 1
                           ENDDO
                        ENDIF

                        ! Index of the same urbclass or the nearest class would be used for new
                        ! year's assignment
                        IF (selfu_ > 0) THEN
                           u_ = selfu_

                        ELSEIF (nurb > 0) THEN
                           duclass = abs ( landurban%settyp(u) - urbclass_(gu_(1)) )
                           u_ = gu_(1)
                           iu = 2
                           DO WHILE (iu .le. nurb)
                              IF (duclass .gt. abs( landurban%settyp(u) - urbclass_(gu_(iu)) )) THEN
                                 u_ = gu_(iu)
                                 duclass = abs( landurban%settyp(u) - urbclass_(u_) )
                              ENDIF
                              iu = iu + 1
                           ENDDO
                        ENDIF

                        IF (u.le.0 .or. u_.le.0) THEN
                           print *, "Error in LuLccMassEnergyConserve URBAN_MODEL!"
                           STOP
                        ENDIF

                        fwsun          (u) = fwsun_          (u_)
                        dfwsun         (u) = dfwsun_         (u_)

                        sroof      (:,:,u) = sroof_      (:,:,u_)
                        swsun      (:,:,u) = swsun_      (:,:,u_)
                        swsha      (:,:,u) = swsha_      (:,:,u_)
                        sgimp      (:,:,u) = sgimp_      (:,:,u_)
                        sgper      (:,:,u) = sgper_      (:,:,u_)
                        slake      (:,:,u) = slake_      (:,:,u_)

                        z_sno_roof   (:,u) = z_sno_roof_   (:,u_)
                        z_sno_gimp   (:,u) = z_sno_gimp_   (:,u_)
                        z_sno_gper   (:,u) = z_sno_gper_   (:,u_)
                        z_sno_lake   (:,u) = z_sno_lake_   (:,u_)

                        dz_sno_roof  (:,u) = dz_sno_roof_  (:,u_)
                        dz_sno_gimp  (:,u) = dz_sno_gimp_  (:,u_)
                        dz_sno_gper  (:,u) = dz_sno_gper_  (:,u_)
                        dz_sno_lake  (:,u) = dz_sno_lake_  (:,u_)

                        lwsun          (u) = lwsun_          (u_)
                        lwsha          (u) = lwsha_          (u_)
                        lgimp          (u) = lgimp_          (u_)
                        lgper          (u) = lgper_          (u_)
                        lveg           (u) = lveg_           (u_)

                        t_roofsno    (:,u) = t_roofsno_    (:,u_)
                        t_wallsun    (:,u) = t_wallsun_    (:,u_)
                        t_wallsha    (:,u) = t_wallsha_    (:,u_)
                        t_gimpsno    (:,u) = t_gimpsno_    (:,u_)
                        t_gpersno    (:,u) = t_gpersno_    (:,u_)
                        t_lakesno    (:,u) = t_lakesno_    (:,u_)

                        troof_inner    (u) = troof_inner_    (u_)
                        twsun_inner    (u) = twsun_inner_    (u_)
                        twsha_inner    (u) = twsha_inner_    (u_)

                        wliq_roofsno (:,u) = wliq_roofsno_ (:,u_)
                        wice_roofsno (:,u) = wice_roofsno_ (:,u_)
                        wliq_gimpsno (:,u) = wliq_gimpsno_ (:,u_)
                        wice_gimpsno (:,u) = wice_gimpsno_ (:,u_)
                        wliq_gpersno (:,u) = wliq_gpersno_ (:,u_)
                        wice_gpersno (:,u) = wice_gpersno_ (:,u_)
                        wliq_lakesno (:,u) = wliq_lakesno_ (:,u_)
                        wice_lakesno (:,u) = wice_lakesno_ (:,u_)

                        sag_roof       (u) = sag_roof_       (u_)
                        sag_gimp       (u) = sag_gimp_       (u_)
                        sag_gper       (u) = sag_gper_       (u_)
                        sag_lake       (u) = sag_lake_       (u_)
                        scv_roof       (u) = scv_roof_       (u_)
                        scv_gimp       (u) = scv_gimp_       (u_)
                        scv_gper       (u) = scv_gper_       (u_)
                        scv_lake       (u) = scv_lake_       (u_)
                        fsno_roof      (u) = fsno_roof_      (u_)
                        fsno_gimp      (u) = fsno_gimp_      (u_)
                        fsno_gper      (u) = fsno_gper_      (u_)
                        fsno_lake      (u) = fsno_lake_      (u_)
                        snowdp_roof    (u) = snowdp_roof_    (u_)
                        snowdp_gimp    (u) = snowdp_gimp_    (u_)
                        snowdp_gper    (u) = snowdp_gper_    (u_)
                        snowdp_lake    (u) = snowdp_lake_    (u_)

                        Fhac           (u) = Fhac_           (u_)
                        Fwst           (u) = Fwst_           (u_)
                        Fach           (u) = Fach_           (u_)
                        Fahe           (u) = Fahe_           (u_)
                        Fhah           (u) = Fhah_           (u_)
                        vehc           (u) = vehc_           (u_)
                        meta           (u) = meta_           (u_)
                        t_room         (u) = t_room_         (u_)
                        t_roof         (u) = t_roof_         (u_)
                        t_wall         (u) = t_wall_         (u_)
                        tafu           (u) = tafu_           (u_)
                        urb_green      (u) = urb_green_      (u_)

                        ! used soil patch value for variable on pervious ground
                        FROM_SOIL = .false.
                        IF (selfu_ < 0) THEN
                           DO k = 1, num
                              IF (patchtype_(frnp_(k)) == 0) THEN
                                 FROM_SOIL = .true.
                              ENDIF
                           ENDDO
                        ENDIF

                        ! Use the first source soil patch temporarily
                        IF (FROM_SOIL) THEN
                           z_sno_gper   (:,u) = z_sno_    (:,frnp_(1))
                           sag_gper       (u) = sag_        (frnp_(1))
                           scv_gper       (u) = scv_        (frnp_(1))
                           fsno_gper      (u) = fsno_       (frnp_(1))
                           snowdp_gper    (u) = snowdp_     (frnp_(1))
                        ENDIF

                        !TODO: need to recalculate wliq_soisno, wice_soisno and scv value - DONE
                        wliq_soisno(: ,np) = 0.
                        wliq_soisno(:1,np) = wliq_roofsno(:1,u )*froof(u)
                        wliq_soisno(: ,np) = wliq_soisno (: ,np) &
                                           + wliq_gpersno(: ,u )*(1-froof(u))*fgper(u)
                        wliq_soisno(:1,np) = wliq_soisno (:1,np) &
                                           + wliq_gimpsno(:1,u )*(1-froof(u))*(1-fgper(u))

                        wice_soisno(: ,np) = 0.
                        wice_soisno(:1,np) = wice_roofsno(:1,u )*froof(u)
                        wice_soisno(: ,np) = wice_soisno (: ,np) &
                                           + wice_gpersno(: ,u )*(1-froof(u))*fgper(u)
                        wice_soisno(:1,np) = wice_soisno (:1,np) &
                                           + wice_gimpsno(:1,u )*(1-froof(u))*(1-fgper(u))

                        scv(np) = scv_roof(u)*froof(u) + scv_gper(u)*(1-froof(u))*fgper(u) &
                                + scv_gimp(u)*(1-froof(u))*(1-fgper(u))

                     ENDIF
#endif

                     ! CALL albland (np, patchtype(np),deltim,&
                     !    soil_s_v_alb(np),soil_d_v_alb(np),soil_s_n_alb(np),soil_d_n_alb(np),&
                     !    chil(patchclass(np)),rho(1:,1:,patchclass(np)),tau(1:,1:,patchclass(np)),&
                     !    fveg(np),green(np),lai(np),sai(np),coszen(np),&
                     !    wt,fsno(np),scv(np),scvold(np),sag(np),ssw,pg_snow(np),&
                     !    forc_t(np),t_grnd(np),t_soisno_(maxsnl+1:,np),&
                     !    dz_soisno_(maxsnl+1:,np),snl,wliq_soisno(maxsnl+1:,np),&
                     !    wice_soisno(maxsnl+1:,np),snw_rds(maxsnl+1:0,np),snofrz,&
                     !    mss_bcpho(maxsnl+1:0,np),mss_bcphi(maxsnl+1:0,np),&
                     !    mss_ocpho(maxsnl+1:0,np),mss_ocphi(maxsnl+1:0,np),&
                     !    mss_dst1(maxsnl+1:0,np),mss_dst2(maxsnl+1:0,np),&
                     !    mss_dst3(maxsnl+1:0,np),mss_dst4(maxsnl+1:0,np),&
                     !    alb(1:,1:,np),ssun(1:,1:,np),ssha(1:,1:,np),ssno(:,:,:,np),&
                     !    thermk(np),extkb(np),extkd(np))


                     IF (allocated(frnp_  )) deallocate(frnp_  )
                     IF (allocated(gu_    )) deallocate(gu_    )
                     IF (allocated(cvsoil_)) deallocate(cvsoil_)
                     IF (allocated(h_     )) deallocate(h_     )
                     np = np + 1
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      IF (p_is_worker) THEN
         IF (allocated(grid_patch_s )) deallocate(grid_patch_s )
         IF (allocated(grid_patch_e )) deallocate(grid_patch_e )
         IF (allocated(grid_patch_s_)) deallocate(grid_patch_s_)
         IF (allocated(grid_patch_e_)) deallocate(grid_patch_e_)
         IF (allocated(locpxl       )) deallocate(locpxl       )
      ENDIF

   END SUBROUTINE LulccMassEnergyConserve

END MODULE MOD_Lulcc_MassEnergyConserve
#endif
! ---------- EOP ------------
