subroutine EULER_INTEGRATION(omgx, omgy, omgz)

use param_phys
use ellipsoid_particle

implicit none

!=========== INPUT VARIABLES =================!

real(kind=8), intent(inout) :: omgx, omgy, omgz

!=============================================!
real(kind=8) :: dt

!!========== Newton Raphson Method ==========!!
! Number of dimensions and total iterations 
integer, parameter:: ITER = 10000

integer :: k, METHOD

! Solution vector OMEGA
real(kind=8), dimension(ndim, 1) :: OMEGA, OMEGA_n_1, OMG

! Function vector F_x
real(kind=8), dimension(ndim, 1) :: F_x

! Jacobian Matrix and its inverse
real(kind=8), dimension(ndim, ndim) :: J_omega, invJ_omega

! Tolerance 
real(kind=8), parameter:: TOL = 1.0e-06

!!============== RK4 Method =================!!
real(kind=8) :: k1_x, k2_x, k3_x, k4_x
real(kind=8) :: k1_y, k2_y, k3_y, k4_y
real(kind=8) :: k1_z, k2_z, k3_z, k4_z

real(kind=8) :: w0_x, w0_y, w0_z
real(kind=8) :: w1_x, w1_y, w1_z
!=============================================!

dt = DT_INIT

METHOD = 2

if (METHOD == 1) then
    k = 1

    ! Initial guess for the solution vector
    OMEGA(1,1) = omgx
    OMEGA(2,1) = omgy
    OMEGA(3,1) = omgz


    do while(k .le. ITER)

        k = k + 1

        ! Initialise the solution from previous iteration
        OMEGA_n_1 = OMEGA

        ! Calculate the function vector
        F_x(1,1) = OMEGA(1,1) - dt*OMEGA(2,1)*OMEGA(3,1)*((IPYY - IPZZ)/IPXX) - omgx
        F_x(2,1) = OMEGA(2,1) - dt*OMEGA(3,1)*OMEGA(1,1)*((IPZZ - IPXX)/IPYY) - omgy
        F_x(3,1) = OMEGA(3,1) - dt*OMEGA(1,1)*OMEGA(2,1)*((IPXX - IPYY)/IPZZ) - omgz

        ! Calculate the Jacobian Matrix
        J_omega(1,1) = 1.0
        J_omega(1,2) = - dt*OMEGA(3,1)*((IPYY - IPZZ)/IPXX)
        J_omega(1,3) = - dt*OMEGA(2,1)*((IPYY - IPZZ)/IPXX)

        J_omega(2,1) = - dt*OMEGA(3,1)*((IPZZ - IPXX)/IPYY)
        J_omega(2,2) = 1.0
        J_omega(2,3) = - dt*OMEGA(1,1)*((IPZZ - IPXX)/IPYY)

        J_omega(3,1) = - dt*OMEGA(2,1)*((IPXX - IPYY)/IPZZ)
        J_omega(3,2) = - dt*OMEGA(1,1)*((IPXX - IPYY)/IPZZ)
        J_omega(3,3) = 1.0

        ! Invert the Jacobian Matrix
        call invert_ndim3_matrix(J_omega, invJ_omega)

        ! Solve the system J.Omg = F_x
        OMG = - matmul(invJ_omega, F_x)

        OMEGA = OMEGA_n_1 + OMG

        if (max(abs(OMG(1,1)), abs(OMG(2,1)), abs(OMG(3,1))) < TOL) then

            GO TO 100

        else if (k == ITER) then

            write(*,*) 'Maximum Iterations Reached'
            write(*,*) 'No Solution found - ERROR'
            stop

        end if

    end do

    100 CONTINUE 


    omgx = OMEGA(1,1)
    omgy = OMEGA(2,1)
    omgz = OMEGA(3,1)


else if (METHOD == 2) then

    ! Initial 
    w0_x = omgx
    w0_y = omgy
    w0_z = omgz

    ! First step
    k1_x = dt*(w0_y*w0_z*((IPYY - IPZZ)/IPXX)) ! omega_x
    k1_y = dt*(w0_z*w0_x*((IPZZ - IPXX)/IPYY)) ! omega_y
    k1_z = dt*(w0_x*w0_y*((IPXX - IPYY)/IPZZ)) ! omega_z

    !Second step
    k2_x = dt*((w0_y + 0.5*k1_y)*(w0_z + 0.5*k1_z)*((IPYY - IPZZ)/IPXX)) ! omega_x
    k2_y = dt*((w0_z + 0.5*k1_z)*(w0_x + 0.5*k1_x)*((IPZZ - IPXX)/IPYY)) ! omega_y
    k2_z = dt*((w0_x + 0.5*k1_x)*(w0_y + 0.5*k1_y)*((IPXX - IPYY)/IPZZ)) ! omega_z

    ! Third Step
    k3_x = dt*((w0_y + 0.5*k2_y)*(w0_z + 0.5*k2_z)*((IPYY - IPZZ)/IPXX)) ! omega_x
    k3_y = dt*((w0_z + 0.5*k2_z)*(w0_x + 0.5*k2_x)*((IPZZ - IPXX)/IPYY)) ! omega_y
    k3_z = dt*((w0_x + 0.5*k2_x)*(w0_y + 0.5*k2_y)*((IPXX - IPYY)/IPZZ)) ! omega_z

    ! Fourth Step
    k4_x = dt*((w0_y + k3_y)*(w0_z + k3_z)*((IPYY - IPZZ)/IPXX)) ! omega_x
    k4_y = dt*((w0_z + k3_z)*(w0_x + k3_x)*((IPZZ - IPXX)/IPYY)) ! omega_y
    k4_z = dt*((w0_x + k3_x)*(w0_y + k3_y)*((IPXX - IPYY)/IPZZ)) ! omega_z


    ! Final Step
    w1_x = w0_x + (k1_x + 2.0*k2_x + 2.0*k3_x + k4_x)/6.0
    w1_y = w0_y + (k1_y + 2.0*k2_y + 2.0*k3_y + k4_y)/6.0
    w1_z = w0_z + (k1_z + 2.0*k2_z + 2.0*k3_z + k4_z)/6.0


    omgx = w1_x
    omgy = w1_y
    omgz = w1_z


end if


!if (impact == 1) then 

!    write(*,*) 'Solution :', omgxn, omgyn, omgzn
!    write(*,*) 'Final Increment :', max(abs(OMG(1,1)), abs(OMG(2,1)), abs(OMG(3,1)))
!    write(*,*) 'Iteration :', k
!    write(*,*) 'Rotational Kinetic Energy: ', 0.5*(IPXX*omgxn*omgxn &
!                                                 + IPYY*omgyn*omgyn &
!                                                 + IPZZ*omgzn*omgzn )
!    write(*,*) ' '

!end if

    
end subroutine EULER_INTEGRATION