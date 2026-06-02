#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_EnKF
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    ensemble Kalman filter (EnKF) types
!
! AUTHOR:
!   Lu Li, 12/2024: Initial version
!   Zhilong Fan, Lu Li, 03/2024: Debug and clean codes
!-----------------------------------------------------------------------------
    USE MOD_Precision
    USE MOD_SPMD_Task
    IMPLICIT NONE
    SAVE

! public functions
    PUBLIC :: letkf


!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

    SUBROUTINE letkf (&
        num_ens, num_obs, &
        hx, y, R, infl, &
        trans)

!-----------------------------------------------------------------------------
! Description:
!   local transform ensemble Kalman filter
!
! Original author :
!   Lu Li, 12/2024
!-----------------------------------------------------------------------------
    USE MOD_Precision
    IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
    integer, intent(in)     :: num_ens                              ! ensemble size
    integer, intent(in)     :: num_obs                              ! number of observations
    real(r8), intent(in)    :: hx(num_obs, num_ens)                 ! ensemble H(x) matrix
    real(r8), intent(in)    :: y(num_obs)                           ! observation vector
    real(r8), intent(in)    :: R(num_obs)                           ! observation error variance
    real(r8), intent(in)    :: infl                                 ! inflation factor
    real(r8), intent(out)   :: trans(num_ens, num_ens)              ! transform matrix (k x k)

!------------------------ Local Variables ------------------------------------
    real(r8) :: hx_mean(num_obs)                    ! mean of ensemble H(x) (l)
    real(r8) :: dhx(num_obs, num_ens)               ! hx - mean(hx) (l x k)
    real(r8) :: dhx_t(num_ens, num_obs)             ! transpose of dhx (k x l)
    real(r8) :: C(num_ens, num_obs)                 ! C = (dhx)^T * (R)^-1 (k x l)
    real(r8) :: M1(num_ens, num_ens)                ! C * dhx (k x k)
    real(r8) :: pa_inv(num_ens, num_ens)            ! inverse of background error covariance matrix (k x k)
    real(r8) :: eigval(num_ens)                     ! eigenvalues of pa_inv (k)
    real(r8) :: eigvec(num_ens, num_ens)            ! eigenvectors of pa_inv (k x k)
    integer  :: lwork
    real(r8), allocatable :: work(:)
    integer  :: err
    real(r8) :: M2(num_ens, num_ens)                ! M2 = eigvec * eigval^-1 (k x k)
    real(r8) :: pa(num_ens, num_ens)                ! background error covariance matrix (k x k)
    real(r8) :: M3(num_ens, num_obs)                ! M3 = pa * C (k x l)
    real(r8) :: delta(num_obs)                      ! increment of observation (l)
    real(r8) :: w_avg(num_ens)                      ! weight (k)
    real(r8) :: M4(num_ens, num_ens)                ! M4 = eigvec * sqrt((k-1)/eigval) (k x k)
    real(r8) :: trans_pert(num_ens, num_ens)        ! perturbation transform matrix (k x k)
    real(r8) :: I0(num_ens, num_ens)                ! identity matrix (k x k)
    real(r8) :: one_div_ens(num_ens)                ! 1/num_ens
    integer  :: i, j

!-----------------------------------------------------------------------------

        ! calculate observation space perturbation
        hx_mean = sum(hx, dim=2) / size(hx, dim=2)   !(lx1)
        DO j = 1, num_ens
            dhx(:,j) = hx(:,j) - hx_mean !(lxk)
        ENDDO

        ! calculate C, intermediate matrix in localized observation
        dhx_t = transpose(dhx) !(kxl)
        DO j = 1, num_obs
            C(:,j) = dhx_t(:,j) / (R(j)) !(kxl)
        ENDDO

        ! calculate C*dhx, intermediate matrix in background error M1
        CALL dgemm('N', 'N', num_ens, num_ens, num_obs, 1.0_8, C, num_ens, dhx, num_obs, 0.0_8, M1, num_ens)

        ! calculate inverse of background error
        pa_inv = M1
        do i=1, num_ens
            pa_inv(i,i) = M1(i,i) + (num_ens-1)*1.0d0/infl
        end do

        ! eigenvalues and eigenvectors of inverse of background error
        lwork = 4 * num_ens
        allocate( work(lwork) )
        CALL dsyev('V', 'U', num_ens, pa_inv, num_ens, eigval, work, lwork, err)

        IF (err < 0) THEN
            WRITE(*,*) 'LETKF error: DSYEV argument ', -err, ' had an illegal value.'
            CALL CoLM_stop()
        ELSE IF (err > 0) THEN
            WRITE(*,*) 'LETKF error: DSYEV failed to converge, INFO = ', err
            CALL CoLM_stop()
        END IF

        eigvec = pa_inv !(kxk)

        ! calculate background error covariance matrix pa = eigvec (eigval)^-1 eigvec^T
        DO i = 1, num_ens
            M2(:,i) = eigvec(:,i) / eigval(i) !(kxk)
        ENDDO
        CALL dgemm('N', 'T', num_ens, num_ens, num_ens, 1.0_8, M2, num_ens, eigvec, num_ens, 0.0_8, pa, num_ens) !(kxk)

        ! caculate pa * C, intermediate matrix in Kalman gain M3
        CALL dgemm('N', 'N', num_ens, num_obs, num_ens, 1.0_8, pa, num_ens, C, num_ens, 0.0_8, M3, num_ens) !(kxl)

        ! calculate weight
        delta = y - hx_mean
        CALL dgemm('N', 'N', num_ens, 1, num_obs, 1.0_8, M3, num_ens, delta, num_obs, 0.0_8, w_avg, num_ens) !(kx1)

        ! calculate pertubation transform matrix
        DO j = 1, num_ens
            M4(:,j) = eigvec(:,j) * sqrt((num_ens-1) / eigval(j))  !(kxk)
        ENDDO
        CALL dgemm('N', 'T', num_ens, num_ens, num_ens, 1.0_8, M4, num_ens, eigvec, num_ens, 0.0_8, trans_pert, num_ens) !(kxk)

        ! calculate transform matrix
        one_div_ens(:) = 1./num_ens
        I0 = -1.0 / num_ens
        DO i = 1, num_ens
            I0(i, i) = 1.0 - 1./num_ens
        ENDDO
        DO j = 1, num_ens
            trans(:, j) = matmul(trans_pert(:, j) + w_avg(:), I0) + one_div_ens(:)
        ENDDO

    END SUBROUTINE letkf

!-----------------------------------------------------------------------------
END MODULE MOD_DA_EnKF
#endif
