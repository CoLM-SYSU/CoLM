#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Assim_TWS
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Data assimilation of terrestrial water storge from GRACE satellite
!
! AUTHOR:
!   Shupeng Zhang: Initial version
!-----------------------------------------------------------------------------
   USE MOD_DataType
   USE MOD_SpatialMapping

   IMPLICIT NONE

   PUBLIC :: init_DA_GRACE
   PUBLIC :: run_DA_GRACE
   PUBLIC :: end_DA_GRACE

   real(r8), allocatable, PUBLIC :: fslp_k_mon (:,:) ! slope factor of runoff
   real(r8), allocatable, PUBLIC :: fslp_k (:) ! slope factor of runoff

   PRIVATE

   character(len=256) :: file_grace
   type(grid_type)    :: grid_grace

   real(r8), allocatable :: longrace(:)
   real(r8), allocatable :: latgrace(:)

   integer :: nobstime
   integer,  allocatable :: obsyear  (:)
   integer,  allocatable :: obsmonth (:)

   type (spatial_mapping_type) :: mg2p_grace

   real(r8), allocatable :: lwe_obs_this (:)
   real(r8), allocatable :: err_obs_this (:)

   real(r8), allocatable :: lwe_obs_prev (:)
   real(r8), allocatable :: err_obs_prev (:)

   real(r8), allocatable :: wat_prev_m (:)
   real(r8), allocatable :: wat_this_m (:)

   real(r8), allocatable :: rnof_acc_prev_m (:)
   real(r8), allocatable :: rnof_acc_this_m (:)
   real(r8), allocatable :: zwt_acc_prev_m  (:)
   real(r8), allocatable :: zwt_acc_this_m  (:)

   real(r8), allocatable :: rnof_prev_m0 (:)
   real(r8), allocatable :: rnof_prev_m1 (:)
   real(r8), allocatable :: rnof_this_m  (:)

   logical,  allocatable :: rnofmask (:)

   logical :: has_prev_grace_obs
   integer :: nac_grace_this, nac_grace_prev

   integer :: year_prev, month_prev

CONTAINS

   ! ----------
   SUBROUTINE init_DA_GRACE ()

   USE MOD_Spmd_Task
   USE MOD_Namelist, only : DEF_DA_OBS_DIR
   USE MOD_Grid
   USE MOD_NetCDFSerial
   USE MOD_Mesh,     only : numelm
   USE MOD_LandElm,  only : landelm
   USE MOD_LandPatch
#ifdef CROP
   USE MOD_LandCrop
#endif
   USE MOD_Pixelset
   USE MOD_Vars_TimeInvariants, only : patchtype
   USE MOD_Forcing, only : forcmask_pch
   USE MOD_RangeCheck
   IMPLICIT NONE

   ! Local Variables

   real(r8), allocatable :: time_real8(:)
   integer :: itime

      file_grace = trim(DEF_DA_OBS_DIR) &
         // '/GRACE_JPL/GRCTellus.JPL.200204_202207.GLO.RL06M.MSCNv02CRI.nc'

      CALL ncio_read_bcast_serial (file_grace, 'time', time_real8)

      nobstime = size(time_real8)
      allocate (obsyear (nobstime))
      allocate (obsmonth(nobstime))

      DO itime = 1, nobstime
         CALL retrieve_yymm_from_days (time_real8(itime), obsyear(itime), obsmonth(itime))
      ENDDO

      IF (p_is_master) THEN
         write(*,*) 'Assimilate GRACE data at'
         DO itime = 1, nobstime
            write(*,*) obsyear(itime), obsmonth(itime)
         ENDDO
      ENDIF

      CALL ncio_read_bcast_serial (file_grace, 'lon', longrace)
      CALL ncio_read_bcast_serial (file_grace, 'lat', latgrace)

      CALL grid_grace%define_by_center (latgrace,longrace)

      CALL mg2p_grace%build_arealweighted (grid_grace, landelm)

      IF (p_is_worker) THEN
         IF (numelm > 0) THEN
            allocate (lwe_obs_this (numelm))
            allocate (err_obs_this (numelm))
            allocate (lwe_obs_prev (numelm))
            allocate (err_obs_prev (numelm))
         ENDIF

         IF (numpatch > 0) THEN
            allocate (wat_prev_m      (numpatch))
            allocate (wat_this_m      (numpatch))
            allocate (rnof_acc_prev_m (numpatch))
            allocate (rnof_acc_this_m (numpatch))
            allocate (zwt_acc_prev_m  (numpatch))
            allocate (zwt_acc_this_m  (numpatch))
            allocate (rnof_prev_m0    (numpatch))
            allocate (rnof_prev_m1    (numpatch))
            allocate (rnof_this_m     (numpatch))
            allocate (rnofmask        (numpatch))

            allocate (fslp_k_mon   (12,numpatch))
            allocate (fslp_k          (numpatch))
         ENDIF
      ENDIF

      IF (p_is_worker) THEN
         CALL elm_patch%build (landelm, landpatch, use_frac = .true.)
      ENDIF

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            rnofmask = patchtype == 0
            IF (DEF_forcing%has_missing_value) THEN
               rnofmask = rnofmask .and. forcmask_pch
            ENDIF
         ENDIF
      ENDIF


      has_prev_grace_obs = .false.

      nac_grace_this = 0
      nac_grace_prev = 0

      IF (p_is_worker) THEN
         wat_this_m     (:) = 0.
         rnof_acc_this_m(:) = 0.
         rnof_this_m    (:) = 0.
         zwt_acc_this_m (:) = 0.
         fslp_k_mon (:,:) = 1.0
         fslp_k (:) = 1.0
      ENDIF

      deallocate (time_real8)

   END SUBROUTINE init_DA_GRACE

   ! ----------
   SUBROUTINE run_DA_GRACE (idate, deltim)

   USE MOD_Spmd_task
   USE MOD_TimeManager
   USE MOD_NetCDFBlock
   USE MOD_Mesh
   USE MOD_LandElm
   USE MOD_LandPatch
   USE MOD_Vars_1DFluxes,       only : rnof, rsur
   USE MOD_Vars_TimeVariables,  only : wat, wa, wdsrf, zwt
   USE MOD_RangeCheck
   USE MOD_UserDefFun
   IMPLICIT NONE

   integer,  intent(in) :: idate(3)
   real(r8), intent(in) :: deltim

   ! Local Variables
   logical :: is_obs_time
   integer :: month, mday, itime, ielm, istt, iend, nextmonth

   real(r8) :: sumpct
   real(r8) :: w1, w0, r1, r0, var_o, var_m, dw_f, dw_o, dw_a, rr, zwt_ave
   real(r8) :: fscal, fprev, fthis

   type(block_data_real8_2d) :: f_grace_lwe ! unit: cm
   type(block_data_real8_2d) :: f_grace_err ! unit: cm

   character(len=256) :: sid, logfile

      CALL julian2monthday (idate(1), idate(2), month, mday)

      is_obs_time = any((obsyear == idate(1)) .and. (obsmonth == month))

      IF (p_is_master) THEN
         IF (is_obs_time) THEN
            write(*,*) 'GRACE at this time.'
         ENDIF
      ENDIF

      IF (p_is_worker) THEN

         IF (has_prev_grace_obs) THEN

            nac_grace_prev = nac_grace_prev + 1

            rnof_acc_prev_m = rnof_acc_prev_m + rnof * deltim
            IF (is_obs_time) THEN
               rnof_prev_m1 = rnof_prev_m1 + rnof_acc_prev_m
            ENDIF

            zwt_acc_prev_m = zwt_acc_prev_m + zwt
         ENDIF

         IF (is_obs_time) THEN

            nac_grace_this = nac_grace_this + 1

            wat_this_m = wat_this_m + wat + wa + wdsrf

            rnof_acc_this_m = rnof_acc_this_m + rnof * deltim
            rnof_this_m = rnof_this_m + rnof_acc_this_m

            zwt_acc_this_m = zwt_acc_this_m + zwt
         ENDIF

      ENDIF

      IF (is_obs_time .and. (isendofmonth(idate, deltim))) THEN

         itime = findloc_ud((obsyear == idate(1)) .and. (obsmonth == month))

         IF (p_is_io) THEN
            CALL allocate_block_data (grid_grace, f_grace_lwe)
            CALL allocate_block_data (grid_grace, f_grace_err)
            CALL ncio_read_block_time (file_grace, 'lwe_thickness', grid_grace, itime, f_grace_lwe)
            CALL ncio_read_block_time (file_grace, 'uncertainty'  , grid_grace, itime, f_grace_err)
         ENDIF

         CALL mg2p_grace%grid2pset (f_grace_lwe, lwe_obs_this)
         CALL mg2p_grace%grid2pset (f_grace_err, err_obs_this)

         IF (p_is_worker) THEN

            lwe_obs_this = lwe_obs_this * 10.0 ! from cm to mm
            err_obs_this = err_obs_this * 10.0 ! from cm to mm

            wat_this_m  = wat_this_m  / nac_grace_this

            zwt_acc_prev_m = zwt_acc_prev_m / nac_grace_prev

            IF (has_prev_grace_obs) THEN
               rnof_prev_m1 = rnof_prev_m1 / nac_grace_this
            ENDIF

            IF (has_prev_grace_obs .and. &
               (((idate(1) == year_prev)  .and. (month_prev == month-1)) &
               .or. ((idate(1) == year_prev+1) .and. (month_prev == 12) .and. (month == 1)))) &
               THEN

               ! write(sid,'(I0)') p_iam_worker
               ! logfile = 'log/grace_log_' // trim(sid) // '.txt'
               ! open(12, file = trim(logfile), position = 'append')

               DO ielm = 1, numelm
                  istt = elm_patch%substt(ielm)
                  iend = elm_patch%subend(ielm)

                  sumpct = sum(elm_patch%subfrc(istt:iend), mask = rnofmask(istt:iend))

                  IF (sumpct <= 0) THEN
                     CYCLE
                  ENDIF

                  w1 = sum(wat_this_m  (istt:iend) * elm_patch%subfrc(istt:iend), &
                     mask = rnofmask(istt:iend)) / sumpct
                  w0 = sum(wat_prev_m  (istt:iend) * elm_patch%subfrc(istt:iend), &
                     mask = rnofmask(istt:iend)) / sumpct
                  r1 = sum(rnof_prev_m1(istt:iend) * elm_patch%subfrc(istt:iend), &
                     mask = rnofmask(istt:iend)) / sumpct
                  r0 = sum(rnof_prev_m0(istt:iend) * elm_patch%subfrc(istt:iend), &
                     mask = rnofmask(istt:iend)) / sumpct

                  zwt_ave = sum(zwt_acc_prev_m(istt:iend) * elm_patch%subfrc(istt:iend), &
                     mask = rnofmask(istt:iend)) / sumpct


                  var_o = err_obs_this(ielm)**2 + err_obs_prev(ielm)**2

                  dw_f = w1 - w0
                  dw_o = lwe_obs_this(ielm) - lwe_obs_prev(ielm)
                  var_m = (dw_f-dw_o)**2 - var_o

                  IF (var_m > 0) THEN

                     dw_a = (var_o * dw_f + var_m * dw_o) / (var_m+var_o)

                     rr = r1 - r0

                     IF (rr > 0) THEN

                        fscal = (1-(dw_a-dw_f)/rr)

                        ! (2) method 2: one parameters adjusted
                        fprev = fslp_k_mon(month,istt)
                        fthis = fprev * fscal
                        fthis = min(max(fthis, fprev*0.5), fprev*2.0)
                        fslp_k_mon(month,istt:iend) = fthis

                        fprev = fslp_k_mon(month_prev,istt)
                        fthis = fprev * fscal
                        fthis = min(max(fthis, fprev*0.5), fprev*2.0)
                        fslp_k_mon(month_prev,istt:iend) = fthis

                        ! write(12,'(I4,I3,I8,8ES11.2)') idate(1), month, landelm%eindex(ielm), &
                        !    dw_o, sqrt(var_o), dw_f, sqrt(var_m), dw_a, &
                        !    rr, zwt_ave, fscal

                     ENDIF

                  ENDIF

               ENDDO

               ! close(12)
            ENDIF

            lwe_obs_prev = lwe_obs_this
            err_obs_prev = err_obs_this

            wat_prev_m = wat_this_m
            wat_this_m = 0.

            rnof_acc_prev_m = rnof_acc_this_m
            rnof_acc_this_m = 0.

            zwt_acc_prev_m = zwt_acc_this_m
            zwt_acc_this_m = 0.

            rnof_prev_m0 = rnof_this_m / nac_grace_this
            rnof_prev_m1 = 0.
            rnof_this_m  = 0.

            nac_grace_prev = nac_grace_this
            nac_grace_this = 0

         ENDIF

         has_prev_grace_obs = .true.
         year_prev  = idate(1)
         month_prev = month

      ENDIF

      IF (isendofmonth(idate, deltim)) THEN
         IF (p_is_worker .and. (numpatch > 0)) THEN
            nextmonth = mod(month+1,12)+1
            fslp_k = fslp_k_mon(nextmonth,:)
         ENDIF
      ENDIF


   END SUBROUTINE run_DA_GRACE

   ! ---------
   SUBROUTINE end_DA_GRACE ()

   IMPLICIT NONE

      IF (allocated(lwe_obs_this))    deallocate(lwe_obs_this)
      IF (allocated(err_obs_this))    deallocate(err_obs_this)
      IF (allocated(lwe_obs_prev))    deallocate(lwe_obs_prev)
      IF (allocated(err_obs_prev))    deallocate(err_obs_prev)
      IF (allocated(wat_prev_m     )) deallocate(wat_prev_m     )
      IF (allocated(wat_this_m     )) deallocate(wat_this_m     )
      IF (allocated(rnof_acc_prev_m)) deallocate(rnof_acc_prev_m)
      IF (allocated(rnof_acc_this_m)) deallocate(rnof_acc_this_m)
      IF (allocated(rnof_prev_m0   )) deallocate(rnof_prev_m0   )
      IF (allocated(rnof_prev_m1   )) deallocate(rnof_prev_m1   )
      IF (allocated(rnof_this_m    )) deallocate(rnof_this_m    )
      IF (allocated(rnofmask       )) deallocate(rnofmask       )

      IF (allocated(fslp_k_mon)) deallocate(fslp_k_mon)
      IF (allocated(fslp_k)) deallocate(fslp_k)

      IF (allocated(longrace)) deallocate(longrace)
      IF (allocated(latgrace)) deallocate(latgrace)

   END SUBROUTINE end_DA_GRACE

   ! ---------
   SUBROUTINE retrieve_yymm_from_days (days, yy, mm)

   IMPLICIT NONE
   real(r8), intent(in)  :: days
   integer,  intent(out) :: yy, mm

   ! Local Variables
   real(r8) :: resday
   integer  :: mdays(12)

      yy = 2002
      mm = 1
      mdays = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

      resday = days
      DO WHILE (resday > mdays(mm))

         resday = resday - mdays(mm)

         mm = mm + 1
         IF (mm > 12) THEN
            yy = yy + 1
            mm = 1
            IF( (mod(yy,4)==0 .and. mod(yy,100)/=0) .or. mod(yy,400)==0 ) THEN
               mdays(2) = 29
            ELSE
               mdays(2) = 28
            ENDIF
         ENDIF

      ENDDO

   END SUBROUTINE retrieve_yymm_from_days

END MODULE MOD_DA_Assim_TWS
#endif
