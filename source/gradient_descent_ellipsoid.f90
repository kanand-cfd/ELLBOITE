subroutine gradient_descent_ellipsoid(alpha, IFLAG, &
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

real(kind=8), dimension(2):: grad_func

real(kind=8), dimension(2) :: pk!,alpha_k, grad_func_dummy

real(kind=8) :: func!, func_dummy, func_dummy2

real(kind=8) :: ETA, GRAD_NORM

!real(kind=8) :: rho, c1

!===================================================!

call eval_function(alpha, func, grad_func, &
                   IXP, IYP, IZP, IQUAT, &
                   JXP, JYP, JZP, JQUAT)


!! Calculate norm of Gradient
GRAD_NORM = sqrt(grad_func(1)*grad_func(1) + grad_func(2)*grad_func(2))

!write(*,*) 'GRAD_NORM', GRAD_NORM
!write(*,*) 'ETA', ETA
!write(*,*) 'alpha', alpha

if(GRAD_NORM < TOL) then
    write(*,*) 'Minimum Norm of Grad', GRAD_NORM
    IFLAG = 0
    RETURN
end if

!write(*,*) 'Calling Gradient Descent'!

pk = - grad_func

ETA = 1.0

alpha =  alpha + ETA*pk

IFLAG = 1

RETURN
end subroutine gradient_descent_ellipsoid


!! Implement Armijo's Condition for calculating step length 
!write(*,*) 'Implement Armijo Condition for calculating step length'
!ETA = 1.0

!c1 = 0.01

!rho = 0.8

!alpha_k = alpha + ETA*pk

!call eval_function(alpha_k, func_dummy, grad_func_dummy)

!func_dummy2 = func + c1*ETA*(grad_func(1)*pk(1) + grad_func(2)*pk(2))

!do while (func_dummy > func_dummy2)
!    ETA = rho*ETA

!    alpha_k = alpha_k + ETA*pk

!    call eval_function(alpha_k, func_dummy, grad_func_dummy)

!    func_dummy2 = func_dummy + c1*ETA*(grad_func_dummy(1)*pk(1) + grad_func_dummy(2)*pk(2))

!enddo

!! End of Armijo Condition
!write(*,*) 'End of Armijo Condition'