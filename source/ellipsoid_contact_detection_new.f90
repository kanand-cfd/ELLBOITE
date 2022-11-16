subroutine ellipsoid_contact_detection(IXP, IYP, IZP, IQUAT, &
                                       JXP, JYP, JZP, JQUAT, &
                                       flag, EPT1, EPT2, cnormal)

use param_phys
use ellipsoid_particle
use mod_quaternion
use minpack_module, only: wp, enorm, lmder1

implicit none

!=================== Input Arguments =================!
! Particle 1
real(kind=8), intent(in) :: IXP, IYP, IZP
type(quaternion), intent(in) :: IQUAT

! Particle 2
real(kind=8), intent(in) :: JXP, JYP, JZP
type(quaternion), intent(in) :: JQUAT

!=================== Output Arguments ================!
logical, intent(inout) :: flag

real(kind=8), dimension(ndim, 1), intent(inout) :: EPT1, EPT2, cnormal
!=====================================================!

!=================== Local  Variables =================!
real(kind=8) :: ti , tf

real(kind=8), dimension(ndim, 1) :: depth, PT1, PT2

real(kind=8), dimension(ndim, 1) :: global_min

real(kind=8) :: contact, norm_depth, dist_center, c_norm

real(kind=8) :: global_min_norm !, XRAND

! Arbitrary unit vector at common normal
real(kind=8), dimension(ndim, 1):: center_normal

real(kind=8), dimension(ndim, 1):: cvec1 !, cvec2

! Ellipsoid centres
real(kind=8), dimension(ndim, 1):: b1, b2

integer :: COUNT

!real(wp) :: x_loc(2)

! ====================== LMDER ======================= !
integer, parameter :: n =2
integer, parameter :: m = 3
integer, parameter :: lwa = 5*n + m

integer :: info
real(wp) :: tol, x(n), fvec(m), fjac(m,n)
integer :: ipvt(n)
real(wp) :: wa(lwa) 
!======================================================!

! Call the CPU time
call CPU_TIME(ti)

COUNT = 0

tol = sqrt(epsilon(1.0))

dist_center = (IXP - JXP)**2 + (IYP - JYP)**2 + (IZP - JZP)**2
dist_center = sqrt(dist_center)


if(dist_center > (2.5*ELL_A)) then
    !STOP 'Too far'
    flag = .FALSE.
    return
end if

! Initial guess 
!!call random_number(XRAND)
center_normal(1,1) = IXP - JXP
center_normal(2,1) = IYP - JYP
center_normal(3,1) = IZP - JZP

c_norm = sqrt(center_normal(1,1)**2 + center_normal(2,1)**2 + center_normal(3,1)**2)
!write(*,*) c_norm

if (c_norm /= 0.0) center_normal = center_normal/ c_norm

x(2) = asin(center_normal(3,1))
x(1) = acos(center_normal(1,1))/cos(x(2))

!100 write(*,*) 'Initial Guess', x


call lmder1(fcn, m, n, x, fvec, fjac, m, tol, info, ipvt, wa, lwa, &
            IXP, IYP, IZP, IQUAT, &
            JXP, JYP, JZP, JQUAT)


!write(*,*) "Solution ", x
!write(*,*) enorm(m, fvec)
!write(*,*) "EXIT parameter:", info


!! ======================================================== !!
! Centre of Ellipsoid 1
b1(1,1) = IXP
b1(2,1) = IYP
b1(3,1) = IZP

! Centre of Ellipsoid 1
b2(1,1) = JXP
b2(2,1) = JYP
b2(3,1) = JZP

call pair_depth(x, b1, IQUAT, &
                b2, JQUAT, &
                cvec1, PT1, PT2)


if((fmargin(PT2, b1, IQUAT) < 1.0) .and. (fmargin(PT1, b2, JQUAT) < 1.0)) then

    ! Penetration Depth
    depth = pt1 - pt2


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    contact = depth(1,1)*cvec1(1,1) + depth(2,1)*cvec1(2,1) + depth(3,1)*cvec1(3,1)

    !write(*,*) 'Contact Condition', contact
    norm_depth = sqrt(depth(1,1)*depth(1,1) + depth(2,1)*depth(2,1) + depth(3,1)*depth(3,1))

    global_min(1,1) =  depth(2,1)*cvec1(3,1) - cvec1(2,1)*depth(3,1)
    global_min(2,1) = -depth(1,1)*cvec1(3,1) + cvec1(1,1)*depth(3,1)
    global_min(3,1) =  depth(1,1)*cvec1(2,1) - cvec1(1,1)*depth(2,1)

    global_min_norm = sqrt(global_min(1,1)**2 + global_min(2,1)**2 + global_min(3,1)**2)

    !write(*,*) "Contact", contact
    !write(*,*) "Ellipsoids are in contact"
    !write(*,*) "Point on Ellipsoid 1", PT1
    !write(*,*) "Point on Ellipsoid 2", PT2
    !write(*,*) ' '
    !write(*,*) "Center of ellipsoid 1", IXP, IYP, IZP
    !write(*,*) "Orientation of Ell 1", IQUAT
    !write(*,*) ' ' 
    !write(*,*) "Center of ellipsoid 2", JXP, JYP, JZP
    !write(*,*) "Orientation of Ell 2", JQUAT
    !write(*,*) ' '
    !write(*,*) "Common normal", cvec1
    !write(*,*) "Penetration Distance", norm_depth
    !write(*,*) 'Minimum found for alpha = ', x, 'based on normal unit vector', cvec1
    !write(*,*) "Condition for global_min", global_min_norm
    !write(*,*) ' '
    !!=== Calculate margin function ===!!
    !write(*,*) 'Margin Ellipsoid 1 point2',fmargin(PT2, b1, IQUAT)  
    !write(*,*) 'Margin Ellipsoid 2 point1',fmargin(PT1, b2, JQUAT)
    !write(*,*) ' '


    EPT1 = PT1
    EPT2 = PT2
    cnormal = cvec1

    flag = .true.


else

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !contact = depth(1,1)*cvec1(1,1) + depth(2,1)*cvec1(2,1) + depth(3,1)*cvec1(3,1)    

    !write(*,*) 'Contact Condition', contact
    !write(*,*) "Point on Ellipsoid 1", PT1
    !write(*,*) "Point on Ellipsoid 2", PT2
    !write(*,*) ' '
    !write(*,*) "Center of ellipsoid 1", IXP, IYP, IZP
    !write(*,*) "Orientation of Ell 1", IQUAT
    !write(*,*) ' '
    !write(*,*) "Center of ellipsoid 2", JXP, JYP, JZP
    !write(*,*) "Orientation of Ell 2", JQUAT
    !write(*,*) ' '
    !write(*,*) "Common normal", cvec1
    !write(*,*) ' '
    !write(*,*) 'Margin Ellipsoid 1 point2',fmargin(PT2, b1, IQUAT)
    !write(*,*) 'Margin Ellipsoid 2 point1',fmargin(PT1, b2, JQUAT)
    !stop

    EPT1 = PT1
    EPT2 = PT2
    cnormal = cvec1

    flag = .false.

end if


call CPU_TIME(tf)
!write(*,*) '  '
!print*, ' Time Elapsed', tf - ti
!write(*,*) '  '  
return


contains


subroutine fcn(m, n, x, fvec, fjac, ldfjac, iflag, &
               IXP, IYP, IZP, IQUAT, &
               JXP, JYP, JZP, JQUAT)
    
    use param_phys
    use ellipsoid_particle
    use mod_quaternion
    use minpack_module, only: wp

    implicit none

    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: ldfjac
    real(wp), intent(in) :: x(n)
    real(wp), intent(inout) :: fvec(m)
    real(wp), intent(inout) :: fjac(ldfjac, n)
    integer, intent(inout) :: iflag

    ! Particle 1
    real(kind=8), intent(in) :: IXP, IYP, IZP
    type(quaternion), intent(in) :: IQUAT
    ! Particle 2
    real(kind=8), intent(in) :: JXP, JYP, JZP
    type(quaternion), intent(in) :: JQUAT
    !! ======================================================== !!
    ! Arbitrary unit vector at common normal
    real(kind=8), dimension(ndim, 1):: cvec1, cvec2

    ! Penetration depth and contact points
    real(kind=8), dimension(ndim, 1):: depth, pt1, pt2

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
    ! Centre of Ellipsoid 1
    b1(1,1) = IXP
    b1(2,1) = IYP
    b1(3,1) = IZP

    ! Centre of Ellipsoid 2
    b2(1,1) = JXP
    b2(2,1) = JYP
    b2(3,1) = JZP

    cvec1(1,1) = cos(x(1))*cos(x(2))
    cvec1(2,1) = sin(x(1))*cos(x(2))
    cvec1(3,1) = sin(x(2))

    !tmp1 = cvec1

    !call transform_basis(cvec1, IQUAT, shape(cvec1))

    !call transform_basis(tmp1, JQUAT, shape(tmp1))
    
    cvec2 = - cvec1!- tmp1

    ! Ellipsoid matrix 1 after rotation
    call calculate_Ellipsoid_matrix(E1_matrix, IQUAT)

    ! Ellipsoid matrix 2 after rotation
    call calculate_Ellipsoid_matrix(E2_matrix, JQUAT)

    tmp1 = matmul(E1_matrix, cvec1)
    tmp2 = matmul(E2_matrix, cvec2)

    lambda1 = 0.25*(cvec1(1,1)*tmp1(1,1) + cvec1(2,1)*tmp1(2,1) + cvec1(3,1)*tmp1(3,1))
    lambda1 = sqrt(lambda1)

    lambda2 = 0.25*(cvec2(1,1)*tmp2(1,1) + cvec2(2,1)*tmp2(2,1) + cvec2(3,1)*tmp2(3,1))
    lambda2 = sqrt(lambda2)

    ! Contact Point 1
    pt1 = (1.0/(2.0*lambda1))*matmul(E1_matrix, cvec1) + b1

    ! Contact Point 2
    pt2 = (1.0/(2.0*lambda2))*matmul(E2_matrix, cvec2) + b2

    ! Penetration Depth
    depth = pt1 - pt2

    ! Optimizing Function
    if(iflag == 1) then

        fvec(1) = depth(1,1)*depth(1,1)
        fvec(2) = depth(2,1)*depth(2,1)
        fvec(3) = depth(3,1)*depth(3,1)

    else
    !!====================== Gradient calculation ==========================!!
    ! Gradient of Unit vector wrt alpha(1)
        Grad_cvec(1,1) = -sin(x(1))*cos(x(2))
        Grad_cvec(2,1) =  cos(x(1))*cos(x(2))
        Grad_cvec(3,1) =  0.0

    ! Gradient of Unit vector wrt alpha(2)
        Grad_cvec(1,2) = -cos(x(1))*sin(x(2))
        Grad_cvec(2,2) = -sin(x(1))*sin(x(2))
        Grad_cvec(3,2) =  cos(x(2))

    ! Gradient of Contact Point 1
        Coeff_grad_pt1 = (1.0/(2.0*lambda1))*E1_matrix - &
                     (1.0/(8.0 * lambda1**3.0))*(matmul(E1_matrix,matmul(cvec1, matmul(transpose(cvec1), E1_matrix))))

        Grad_pt1 = matmul(Coeff_grad_pt1, Grad_cvec)

    ! Gradient of Contact Point 2
        Coeff_grad_pt2 = (1.0/(2.0*lambda2))*E2_matrix - &
                     (1.0/(8.0 * lambda2**3.0))*(matmul(E2_matrix,matmul(cvec2, matmul(transpose(cvec2), E2_matrix))))

        Grad_pt2 = matmul(Coeff_grad_pt2, -Grad_cvec)

    ! Gradient of Penetration Depth
        Grad_depth = Grad_pt1 - Grad_pt2

    ! Jacobian of Optimizing Function
        fjac(1,1) = 2.0*(depth(1,1)*Grad_depth(1,1))
        fjac(2,1) = 2.0*(depth(2,1)*Grad_depth(2,1))
        fjac(3,1) = 2.0*(depth(3,1)*Grad_depth(3,1))

        fjac(1,2) = 2.0*(depth(1,1)*Grad_depth(1,2))
        fjac(2,2) = 2.0*(depth(2,1)*Grad_depth(2,2))
        fjac(3,2) = 2.0*(depth(3,1)*Grad_depth(3,2))

    end if
    
end subroutine fcn

end subroutine ellipsoid_contact_detection