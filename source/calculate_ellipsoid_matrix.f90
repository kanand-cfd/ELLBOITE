subroutine calculate_ELLIPSOID_matrix(e_matrix, q_n)

use param_phys
use ellipsoid_particle
use mod_quaternion

implicit none

! The distance of the impact point from the particle center
real(kind=8), dimension(ndim, ndim), intent(inout) :: e_matrix

! Input Orientation at nth time step
type(quaternion), intent(in) :: q_n

! Euler parameters
real(kind=8) :: q0, q1, q2, q3

! Global Rotation Matrix and its transpose
real(kind=8), dimension(ndim, ndim) :: R, Trn_R !, Temp1, Temp2

! Matrix to represent the basic and general rotated Ellipsoid
real(kind=8), dimension(ndim, ndim) :: A_base, A_gen

!===============================================================================!

! Initializing the Matrices
A_base(:,:) = 0.0
A_gen(:,:) = 0.0


A_base(1,1) = ELL_A*ELL_A
A_base(2,2) = ELL_B*ELL_B
A_base(3,3) = ELL_C*ELL_C

! Initialise the Euler Parameters
q0 = q_n%a
q1 = q_n%b
q2 = q_n%c
q3 = q_n%d

! Initialise Rotation Matrix
R(1,1) = 1.0 - 2.0*(q2*q2 + q3*q3)
R(1,2) = 2.0*(q1*q2 - q0*q3)
R(1,3) = 2.0*(q1*q3 + q0*q2)

R(2,1) = 2.0*(q1*q2 + q0*q3)
R(2,2) = 1.0 - 2.0*(q1*q1 + q3*q3)
R(2,3) = 2.0*(q2*q3 - q0*q1)

R(3,1) = 2.0*(q1*q3 - q0*q2)
R(3,2) = 2.0*(q2*q3 + q0*q1)
R(3,3) = 1.0 - 2.0*(q1*q1 + q2*q2)

! Transpose of R 
Trn_R = TRANSPOSE(R) ! R is orthogonal so Trn_R = Inverse(R)

! General Rotated Matrix
A_gen = MATMUL(R, MATMUL(A_base, Trn_R))

e_matrix = A_gen

end subroutine calculate_ELLIPSOID_matrix