module ellipsoid_particle

use param_phys
use mod_quaternion

implicit none

type ELL_PART

! Particle position
real(kind=8):: XP
real(kind=8):: YP
real(kind=8):: ZP

! Particle velocity
real(kind=8):: UP
real(kind=8):: VP
real(kind=8):: WP

! Particle Orientation Quaternions
real(kind=8):: QUAT_A
real(kind=8):: QUAT_B
real(kind=8):: QUAT_C
real(kind=8):: QUAT_D

type(quaternion):: ELLQUAT 

! Particle angular velocity 
real(kind=8):: OMEGAX
real(kind=8):: OMEGAY
real(kind=8):: OMEGAZ

! Global Angular Velocity
real(kind=8) :: GOMEGAX
real(kind=8) :: GOMEGAY
real(kind=8) :: GOMEGAZ

integer :: IDP

! Color flag for Collision
real :: COLOR

end type ELL_PART


public :: fmargin

private :: margin_function

interface fmargin 
    module procedure margin_function 
end interface

type(ELL_PART), dimension(:), allocatable :: ELLP 


contains

function margin_function(point, center, quat) result(m)

    real(kind=8) :: m

    real(kind=8), dimension(1, 1) :: margin_vec

    real(kind=8) , dimension(ndim, 1), intent(in) :: point, center

    type(quaternion), intent(in) :: quat

    real(kind=8), dimension(ndim, ndim) :: ELL, R, Trn_R

    real(kind=8), dimension(ndim, 1) :: ELL_, vec

    ! Initialize the ellipsoid matrix
    ELL(:,:) = 0.0
    ELL(1,1) = 1.0/(ELL_A*ELL_A)
    ELL(2,2) = 1.0/(ELL_B*ELL_B)
    ELL(3,3) = 1.0/(ELL_C*ELL_C)

    !write(*,*) "Quaternion : ", quat%a, quat%b, quat%c, quat%d
    !write(*,*) "Norm :", norm_q(q_p)

    ! Initialise Rotation Matrix
    R(1,1) = 1.0 - 2.0*(quat%c*quat%c + quat%d*quat%d)
    R(1,2) = 2.0*(quat%b*quat%c - quat%a*quat%d)
    R(1,3) = 2.0*(quat%b*quat%d + quat%a*quat%c)

    R(2,1) = 2.0*(quat%b*quat%c + quat%a*quat%d)
    R(2,2) = 1.0 - 2.0*(quat%b*quat%b + quat%d*quat%d)
    R(2,3) = 2.0*(quat%c*quat%d - quat%a*quat%b)

    R(3,1) = 2.0*(quat%b*quat%d - quat%a*quat%c)
    R(3,2) = 2.0*(quat%c*quat%d + quat%a*quat%b)
    R(3,3) = 1.0 - 2.0*(quat%b*quat%b + quat%c*quat%c)

    ! Transpose of R 
    Trn_R = TRANSPOSE(R) ! R is orthogonal so Trn_R = Inverse(R)

    ELL = matmul(R, matmul(ELL, Trn_R))

    vec = point - center

    ELL_ = matmul(ELL, vec)

    margin_vec = matmul(transpose(vec), ELL_)

    m = margin_vec(1,1)

end function margin_function


end module ellipsoid_particle