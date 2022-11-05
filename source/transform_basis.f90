subroutine transform_basis(A_b, q_p, dims)

use param_phys
use mod_quaternion

implicit none

!===========================================================!
!===========================================================!
! Dimensions
!integer, parameter :: ndim = 3

integer, dimension(2), intent(in) :: dims

! Matrix to be transformed
real(kind=8), dimension(dims(1), dims(2)), intent(inout) :: A_b

! Orientation in space
type(quaternion), intent(in) :: q_p 

! Temp Quaternion
type(quaternion) :: qp_temp, qp_temp1

! Temp Matrix
real(kind=8), dimension(dims(1), dims(2)) :: A_Tr

! Global Rotation Matrix and its transpose
real(kind=8), dimension(ndim, ndim) :: R, Trn_R 

real(kind=8) :: q0, q1, q2, q3
!===========================================================!
!===========================================================!
! Initialise the Euler Parameters
q0 = q_p%a
q1 = q_p%b
q2 = q_p%c
q3 = q_p%d

!write(*,*) "Quaternion : ", q0, q1, q2, q3
!write(*,*) "Norm :", norm_q(q_p)

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

!===========================================================!
!===========================================================!
! case 1 : Vector (3x1)
if(dims(2) == 1) then

    qp_temp%a = 0.0
    qp_temp%b = A_b(1,1)
    qp_temp%c = A_b(2,1)
    qp_temp%d = A_b(3,1)

    qp_temp1 = mult_quat(q_p, mult_quat(qp_temp, conj_q(q_p)))

    !if (qp_temp1%a /= 0.0) write(*,*) "Error in quaternion to vector transform", qp_temp1%a

    A_Tr(1,1) = qp_temp1%b
    A_Tr(2,1) = qp_temp1%c
    A_Tr(3,1) = qp_temp1%d

! case 2 : Tensor (3x3)
else

    A_Tr = matmul(R, matmul(A_b, Trn_R))

end if



A_b = A_Tr

!write(*,*) 'Rotation Matrix'
!write(*,*) R(1,:)
!write(*,*) R(2,:)
!write(*,*) R(3,:)

!write(*,*) 'transform_basis of Rotation Matrix'
!write(*,*) Trn_R(1,:)
!write(*,*) Trn_R(2,:)
!write(*,*) Trn_R(3,:)


!write(*,*) A_b(1,:)
!write(*,*) A_b(2,:)
!write(*,*) A_b(3,:)

end subroutine transform_basis