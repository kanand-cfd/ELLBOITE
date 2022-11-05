subroutine ELLIPSOID_CONTACT_OPTIMIZATION(alpha, IFLAG, &
                                      IXP, IYP, IZP, IQUAT, &
                                      JXP, JYP, JZP, JQUAT)

use param_phys
use ellipsoid_particle
use mod_quaternion

implicit none

!===================================================!
real(kind=8), parameter :: TOL = 1.0E-6

real(kind=8), dimension(2), intent(inout) :: alpha

integer, intent(inout) :: IFLAG
! Particle 1
real(kind=8), intent(in) :: IXP, IYP, IZP
type(quaternion), intent(in) :: IQUAT

! Particle 2
real(kind=8), intent(in) :: JXP, JYP, JZP
type(quaternion), intent(in) :: JQUAT

!===================================================!

real(kind=8) :: func

real(kind=8), dimension(2) :: grad_func

real(kind=8), dimension(2,2) :: hessian_func, inv_hess_func, approx_hessian

real(kind=8) :: ETA, GRAD_NORM

real(kind=8), dimension(2,1) :: pk, grad_func_mat

real(kind=8) :: mean, det, eign1, eign2
!===================================================!

ETA = 10.0

call eval_function2(alpha, &
                    func, grad_func, hessian_func, &
                    IXP, IYP, IZP, IQUAT, &
                    JXP, JYP, JZP, JQUAT)

!! Calculate norm of Gradient
GRAD_NORM = sqrt(grad_func(1)*grad_func(1) + grad_func(2)*grad_func(2))
!write(*,*) GRAD_NORM

if(GRAD_NORM < TOL) then
    write(*,*) 'Minimum Norm of Grad', GRAD_NORM
    IFLAG = 0
    RETURN
end if

grad_func_mat(1,1) = grad_func(1)
grad_func_mat(2,1) = grad_func(2)

!! Finding out if hessian matrix is positive definite !!
! Calculate Determinant of hessian 
det = hessian_func(1,1)*hessian_func(2,2) - hessian_func(1,2)*hessian_func(2,1)

approx_hessian = matmul(grad_func_mat, transpose(grad_func_mat))

if((det > 0) .and. hessian_func(1,1) > 0) then ! Hessian is positive definite

    inv_hess_func(1,1) =   hessian_func(2,2)
    inv_hess_func(1,2) = - hessian_func(1,2)
    inv_hess_func(2,1) = - hessian_func(2,1)
    inv_hess_func(2,2) =   hessian_func(1,1)

    inv_hess_func = (1.0/det)*inv_hess_func

    ! Solve the linear system [Hess_func] * pk = - [Grad_func]
    pk = - matmul(inv_hess_func, grad_func_mat)

else

    inv_hess_func(1,1) =   approx_hessian(2,2)
    inv_hess_func(1,2) = - approx_hessian(1,2)
    inv_hess_func(2,1) = - approx_hessian(2,1)
    inv_hess_func(2,2) =   approx_hessian(1,1)

    inv_hess_func = (1.0/det)*inv_hess_func

    ! Solve the linear system [Hess_func] * pk = - [Grad_func]
    pk = - matmul(inv_hess_func, grad_func_mat)



!    STOP 'Hessian not positive definite'
end if


alpha = alpha + ETA*pk(:,1)

IFLAG = 1

!stop
RETURN
end subroutine ELLIPSOID_CONTACT_OPTIMIZATION