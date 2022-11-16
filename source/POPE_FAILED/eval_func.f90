subroutine eval_function(alpha, func, Grad_func, &
                         IXP, IYP, IZP, IQUAT, &
                         JXP, JYP, JZP, JQUAT)

use param_phys
use ellipsoid_particle
use mod_quaternion

implicit none

! Optimization Parameters
real(kind=8), dimension(2), intent(in) :: alpha

! Optimization Function
real(kind=8), intent(inout) :: func

! Gradient of Optimization Function 
real(kind=8), dimension(2), intent(inout) :: Grad_func

! Particle 1
real(kind=8), intent(in) :: IXP, IYP, IZP
type(quaternion), intent(in) :: IQUAT
! Particle 2
real(kind=8), intent(in) :: JXP, JYP, JZP
type(quaternion), intent(in) :: JQUAT
!! ======================================================== !!

! Penetration depth and contact points
real(kind=8), dimension(ndim, 1):: depth, pt1, pt2
 
! Arbitrary unit vector at common normal
real(kind=8), dimension(ndim, 1):: cvec1, cvec2

! Ellipsoid centres
real(kind=8), dimension(ndim, 1):: b1, b2

! Ellipsoid matrix after rotation
real(kind=8), dimension(ndim, ndim) :: E1_matrix, E2_matrix

! Scalar normal parameter
real(kind=8):: lambda1, lambda2

! Temporary vectors
real(kind=8), dimension(ndim, 1) :: tmp1, tmp2

! Gradient of unit vector at common normal, contact points, depth of penetration  
real(kind=8), dimension(ndim, 2) :: Grad_cvec, Grad_pt1, Grad_pt2, Grad_depth

! matrix coefficient for calculation of contact point gradient
real(kind=8), dimension(ndim, ndim) :: Coeff_grad_pt1, Coeff_grad_pt2

!! ======================================================== !!

! Unit vector at Common Normal
cvec1(1,1) = cos(alpha(1))*cos(alpha(2))
cvec1(2,1) = sin(alpha(1))*cos(alpha(2))
cvec1(3,1) = sin(alpha(2))

cvec2 = - cvec1

!write(*,*) cvec

! Centre of Ellipsoid 1
b1(1,1) = IXP
b1(2,1) = IYP
b1(3,1) = IZP

!write(*,*) 'b1', b1

! Centre of Ellipsoid 2
b2(1,1) = JXP
b2(2,1) = JYP
b2(3,1) = JZP

!write(*,*) 'b2', b2

! Ellipsoid matrix 1 after rotation
call calculate_Ellipsoid_matrix(E1_matrix, IQUAT)
!write(*,*) 'Ellipsoid matrix 1'
!write(*,*) E1_matrix(1,:)
!write(*,*) E1_matrix(2,:)
!write(*,*) E1_matrix(3,:)

! Ellipsoid matrix 2 after rotation
call calculate_Ellipsoid_matrix(E2_matrix, JQUAT)
!write(*,*) 'Ellipsoid matrix 2'
!write(*,*) E2_matrix(1,:)
!write(*,*) E2_matrix(2,:)
!write(*,*) E2_matrix(3,:)

tmp1 = matmul(E1_matrix, cvec1)
tmp2 = matmul(E2_matrix, cvec2)

lambda1 = 0.25*(cvec1(1,1)*tmp1(1,1) + cvec1(2,1)*tmp1(2,1) + cvec1(3,1)*tmp1(3,1))

lambda2 = 0.25*(cvec2(1,1)*tmp2(1,1) + cvec2(2,1)*tmp2(2,1) + cvec2(3,1)*tmp2(3,1))

lambda1 = sqrt(lambda1)
!write(*,*) 'lambda1', lambda1

lambda2 = sqrt(lambda2)
!write(*,*) 'lambda2', lambda2

! Contact Point 1
pt1 = (1.0/(2.0*lambda1))*matmul(E1_matrix, cvec1) + b1
!write(*,*) 'pt1', pt1

! Contact Point 2
pt2 = (1.0/(2.0*lambda2))*matmul(E2_matrix, cvec2) + b2
!write(*,*) 'pt2', pt2

! Penetration Depth
depth = pt1 - pt2
!write(*,*) 'depth', depth

! Optimizing Function
func = depth(1,1)*depth(1,1) + depth(2,1)*depth(2,1) + depth(3,1)*depth(3,1)


!write(*,*) 'func', func

!!====================== Gradient calculation ==========================!!
! Gradient of Unit vector wrt alpha(1)
Grad_cvec(1,1) = -sin(alpha(1))*cos(alpha(2))
Grad_cvec(2,1) = cos(alpha(1))*cos(alpha(2))
Grad_cvec(3,1) = 0.0

! Gradient of Unit vector wrt alpha(2)
Grad_cvec(1,2) = -cos(alpha(1))*sin(alpha(2))
Grad_cvec(2,2) = -sin(alpha(1))*sin(alpha(2))
Grad_cvec(3,2) = cos(alpha(2))

! Gradient of Contact Point 1
Coeff_grad_pt1 = (1.0/(2.0*lambda1))*E1_matrix - &
                 (1.0/(8.0 * lambda1**3.0))*(matmul(E1_matrix,matmul(cvec1, matmul(transpose(cvec1), E1_matrix))))

Grad_pt1 = matmul(Coeff_grad_pt1, Grad_cvec)

!write(*,*) 'Grad_pt1'
!write(*,*) Grad_pt1(1,:)
!write(*,*) Grad_pt1(2,:)
!write(*,*) Grad_pt1(3,:)

! Gradient of Contact Point 2
Coeff_grad_pt2 = (1.0/(2.0*lambda2))*E2_matrix - &
                 (1.0/(8.0 * lambda2**3.0))*(matmul(E2_matrix,matmul(cvec2, matmul(transpose(cvec2), E2_matrix))))

Grad_pt2 = matmul(Coeff_grad_pt2, -Grad_cvec)

!write(*,*) 'Grad_pt2'
!write(*,*) Grad_pt2(1,:)
!write(*,*) Grad_pt2(2,:)
!write(*,*) Grad_pt2(3,:)

! Gradient of Penetration Depth
Grad_depth = Grad_pt1 - Grad_pt2

!write(*,*) 'Grad_depth'
!write(*,*) Grad_depth(1,:)
!write(*,*) Grad_depth(2,:)
!write(*,*) Grad_depth(3,:)

! Gradient of Optimizing Function
Grad_func(1) = 2.0*(depth(1,1)*Grad_depth(1,1) + depth(2,1)*Grad_depth(2,1) + depth(3,1)*Grad_depth(3,1))
Grad_func(2) = 2.0*(depth(1,1)*Grad_depth(1,2) + depth(2,1)*Grad_depth(2,2) + depth(3,1)*Grad_depth(3,2))

!write(*,*) 'Grad_func', Grad_func




end subroutine eval_function