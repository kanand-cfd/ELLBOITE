subroutine QUATERNION_INTEGRATION(qp_n, omgx, omgy, omgz)

!========================================================!
! This integration method is based on Direct Multiplica- !
! -tion Method introduced in the patent by Whitmore.     !
! Refer to Zhang, van Wachem, Acta Mech. 2013			 !
!========================================================!

use param_phys
use mod_quaternion

implicit none

!========================================================!
! 					Input Variables 					 !
!========================================================!
real(kind=8), intent(in) :: omgx, omgy, omgz

type(quaternion), intent(inout) :: qp_n

!========================================================!
! 					Local Variables 					 !
!========================================================!
type(quaternion) :: qp_map, qp_next

real(kind=8) :: norm_omg

real(kind=8):: dt

real(kind=8) :: q0, q1, q2, q3

real(kind=8) :: opx, opy, opz

integer :: METHOD ! 1 = Euler, 2 = Direct Multiplication Method 


dt = DT_INIT

METHOD = 1
!========================================================!
!               Direct Multiplication Method             !
!========================================================!
if (METHOD == 1) then
    ! Calculate the norm of Angular velocity
    norm_omg = sqrt(omgx*omgx + omgy*omgy + omgz*omgz)
    !write(*,*) norm_omg

    if (norm_omg /= 0) then
        ! Calculation of exponential map converted to quaternion
        q0 = cos((norm_omg*dt)/2.0)

        q1 = sin((norm_omg*dt)/2.0)*(omgx/norm_omg) 

        q2 = sin((norm_omg*dt)/2.0)*(omgy/norm_omg)

        q3 = sin((norm_omg*dt)/2.0)*(omgz/norm_omg)

    else

        q0 = 1.0
        q1 = 0.0
        q2 = 0.0
        q3 = 0.0

    end if ! if (norm_omg /= 0) then

    ! Update the quaternion map  
    qp_map%a = q0
    qp_map%b = q1
    qp_map%c = q2
    qp_map%d = q3

    ! Update the Quaternion at next time-step by direct multiplication
    qp_next = mult_quat(qp_map, qp_n)

    qp_n = qp_next



else
!========================================================!
!                      Euler-Method                      !
!========================================================!
    q0 = qp_n%a
    q1 = qp_n%b
    q2 = qp_n%c
    q3 = qp_n%d

    opx = omgx
    opy = omgy
    opz = omgy

    qp_map%a = q0 + dt*(-q1*opx - q2*opy - q3*opz)
    qp_map%b = q1 + dt*(q0*opx - q3*opy + q2*opz)
    qp_map%c = q2 + dt*(q3*opx + q0*opy - q1*opz)
    qp_map%d = q3 + dt*(-q2*opx + q1*opy + q0*opz)

    qp_map = unit_quat(qp_map)

    qp_next = qp_map

    qp_n = qp_next


end if ! if (METHOD == 1) then

!write(*,*) qp_next
!write(*,*) norm_q(qp_next)
!stop

end subroutine QUATERNION_INTEGRATION