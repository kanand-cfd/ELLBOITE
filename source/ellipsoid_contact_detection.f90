subroutine ELLIPSOID_CONTACT_DETECTION(IXP, IYP, IZP, IQUAT, &
                                       JXP, JYP, JZP, JQUAT, &
                                       flag, EPT1, EPT2, cnormal)

use param_phys
use mod_quaternion
use ellipsoid_particle

implicit none

!=================== Input Parameters =================!
! Particle 1
real(kind=8), intent(in) :: IXP, IYP, IZP
type(quaternion), intent(in) :: IQUAT

! Particle 2
real(kind=8), intent(in) :: JXP, JYP, JZP
type(quaternion), intent(in) :: JQUAT

logical, intent(inout) :: flag

real(kind=8), dimension(ndim, 1), intent(inout) :: EPT1, EPT2, cnormal
!======================================================!

!=================== Local  Variables =================!
real(kind=8) :: ti , tf

real(kind=8), dimension(ndim, 1) :: depth, PT1, PT2

real(kind=8), dimension(ndim, 1) :: normal1, normal2

real(kind=8), dimension(ndim, 1) :: global_min

real(kind=8) :: contact, norm_depth, dist_center

real(kind=8) :: global_min_norm, XRAND

! =====================Used for minimizing function======================== !

! Arbitrary unit vector at common normal
real(kind=8), dimension(ndim, 1):: center_normal

real(kind=8), dimension(ndim, 1):: cvec1, cvec2

! Ellipsoid centres
real(kind=8), dimension(ndim, 1):: b1, b2

! Ellipsoid matrix after rotation
real(kind=8), dimension(ndim, ndim) :: E1_matrix, E2_matrix

! Scalar normal parameter
real(kind=8):: lambda1, lambda2

! Temporary vectors
real(kind=8), dimension(ndim, 1) :: tmp1, tmp2

! ====================== Gradient Descent ======================= !

! Optimization parameters
real(kind=8), dimension(2):: alpha

! Optimization function and gradient
real(kind=8) :: c_norm

! Parameters for gradient descent
integer :: ITER, IFLAG, J

!======================================================!
! Call the CPU time
call CPU_TIME(ti)


dist_center = (IXP - JXP)**2 + (IYP - JYP)**2 + (IZP - JZP)**2
dist_center = sqrt(dist_center)


if(dist_center > (5.0*ELL_A)) then
    !STOP 'Too far'
    flag = .FALSE.
    return
end if



!call minimize_function(depth, PT1, PT2, normal1, normal2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Call Gradient Descent

! Initial guess
call random_number(XRAND)
center_normal(1,1) = XRAND !IXP - JXP !JXP - IXP

call random_number(XRAND)
center_normal(2,1) = XRAND !IYP - JYP !JYP - IYP

call random_number(XRAND)
center_normal(3,1) = XRAND !IZP - JZP !JZP - IZP

c_norm = sqrt(center_normal(1,1)**2 + center_normal(2,1)**2 + center_normal(3,1)**2)

!write(*,*) c_norm

if (c_norm /= 0.0) center_normal = center_normal/ c_norm


!100 
alpha(2) = asin(center_normal(3,1))
alpha(1) = acos(center_normal(1,1))/cos(alpha(2))


write(*,*) 'Initial Guess', alpha, 'based on normal unit vector', center_normal

IFLAG = 0

do J = 1, 1000000

    call gradient_descent_ellipsoid(alpha, IFLAG, IXP, IYP, IZP, IQUAT, JXP, JYP, JZP, JQUAT)

    !call ELLIPSOID_CONTACT_OPTIMIZATION(alpha, IFLAG, IXP, IYP, IZP, IQUAT, JXP, JYP, JZP, JQUAT)


    if(IFLAG .eq. 0) then
        !call eval_function(alpha, func, grad_func)
    !else
        ITER = J
        GO TO 50
    end if
    !write(*,*) 'ITERATIONS: ', J
    !write(*,*) 'alpha:', alpha
    !write(*,*) ' '
end do

STOP "TOO MANY ITERATIONS"

50 write(*,*) 'Minimum found at ITERATION =', ITER, 'for alpha = ', alpha

!===================== Fin de Gradient Descent ========================!
! Unit vector at Common Normal
cvec1(1,1) = cos(alpha(1))*cos(alpha(2))
cvec1(2,1) = sin(alpha(1))*cos(alpha(2))
cvec1(3,1) = sin(alpha(2))

!write(*,*) 'Minimum found at ITERATION =', ITER, 'for alpha = ', alpha, 'based on normal unit vector', cvec1

100 cvec2 = - cvec1

! Centre of Ellipsoid 1
b1(1,1) = IXP
b1(2,1) = IYP
b1(3,1) = IZP

! Centre of Ellipsoid 2
b2(1,1) = JXP
b2(2,1) = JYP
b2(3,1) = JZP

! Ellipsoid matrix 1 after rotation
call calculate_Ellipsoid_matrix(E1_matrix, IQUAT)

! Ellipsoid matrix 2 after rotation
call calculate_Ellipsoid_matrix(E2_matrix, JQUAT)

tmp1 = matmul(E1_matrix, cvec1)
tmp2 = matmul(E2_matrix, cvec2)

lambda1 = 0.25*(cvec1(1,1)*tmp1(1,1) + cvec1(2,1)*tmp1(2,1) + cvec1(3,1)*tmp1(3,1))
lambda2 = 0.25*(cvec2(1,1)*tmp2(1,1) + cvec2(2,1)*tmp2(2,1) + cvec2(3,1)*tmp2(3,1))

lambda1 = sqrt(lambda1)
lambda2 = sqrt(lambda2)

! Contact Point 1
PT1 = (1.0/(2.0*lambda1))*matmul(E1_matrix, cvec1) + b1

! Contact Point 2
PT2 = (1.0/(2.0*lambda2))*matmul(E2_matrix, cvec2) + b2

! Penetration Depth
depth = PT1 - PT2

! Common normal
normal1 = cvec1
normal2 = cvec2

!if((depth(1,1)*nrm1(1,1) + depth(2,1)*nrm1(2,1) + depth(3,1)*nrm1(3,1)) < 0.0 .and. J > 500) RETURN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contact = depth(1,1)*normal1(1,1) + depth(2,1)*normal1(2,1) + depth(3,1)*normal1(3,1)

!write(*,*) 'Contact Condition', contact
norm_depth = sqrt(depth(1,1)*depth(1,1) + depth(2,1)*depth(2,1) + depth(3,1)*depth(3,1))

global_min(1,1) =  depth(2,1)*normal1(3,1) - normal1(2,1)*depth(3,1)
global_min(2,1) = -depth(1,1)*normal1(3,1) + normal1(1,1)*depth(3,1)
global_min(3,1) =  depth(1,1)*normal1(2,1) - normal1(1,1)*depth(2,1)

global_min_norm = sqrt(global_min(1,1)**2 + global_min(2,1)**2 + global_min(3,1)**2)

if (contact .ge. 0.0) then

    !write(*,*) "Contact", contact
    write(*,*) "Ellipsoids are in contact"
    !write(*,*) "Point on Ellipsoid 1", PT1
    !write(*,*) "Point on Ellipsoid 2", PT2
    write(*,*) "Penetration Distance", norm_depth
    !write(*,*) 'Minimum found at ITERATION =', ITER
    !write(*,*) "Common Normal", normal1
    !write(*,*) 'Minimum found at ITERATION =', ITER, 'for alpha = ', alpha, 'based on normal unit vector', cvec1
    write(*,*) "Condition for global_min", global_min_norm

    !!=== Calculate margin function ===!!
    !if(fmargin(PT2, b1, IQUAT) < 1.0) 
    write(*,*) 'Margin Ellipsoid 1 point2',fmargin(PT2, b1, IQUAT)  
    !if(fmargin(PT1, b2, JQUAT) < 1.0) 
    write(*,*) 'Margin Ellipsoid 2 point1',fmargin(PT1, b2, JQUAT)
    write(*,*) ' '

    ! Output for collision subroutine
    flag = .TRUE.    

    EPT1 = PT1
    EPT2 = PT2
    cnormal = normal1

    CONTTACT = CONTTACT + 1


else

    !write(*,*) "Contact", contact
    write(*,*) "Ellipsoid do not touch each other"
    !write(*,*) "Distance", norm_depth
    !write(*,*) "Common Normal", normal1
    !write(*,*) "Point on Ellipsoid 1", PT1
    !write(*,*) "Point on Ellipsoid 2", PT2
    !write(*,*) 'Initial Guess', asin(center_normal(3,1)), acos(center_normal(1,1))/cos(alpha(2)), 'based on normal unit vector', center_normal
    !write(*,*) 'Minimum found at ITERATION =', ITER, 'for alpha = ', alpha, 'based on normal unit vector', cvec1
    !write(*,*) "Condition for global_min", global_min_norm
    write(*,*) 'Margin Ellipsoid 1 point2',fmargin(PT2, b1, IQUAT)  
    write(*,*) 'Margin Ellipsoid 2 point1',fmargin(PT1, b2, JQUAT)
    !write(*,*) ' '
    !flag = .FALSE.

    !NO_CONTACT = NO_CONTACT + 1

    !EPT1 = PT1
    !EPT2 = PT2
    !cnormal = normal1
    
    cvec1 = - cvec1
    !center_normal = - cvec1

    GO TO 100

end if


call CPU_TIME(tf)
!write(*,*) '  '
!print*, ' Time Elapsed', tf - ti
!write(*,*) '  '  
stop
end subroutine ELLIPSOID_CONTACT_DETECTION