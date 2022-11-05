subroutine COLLISION_RESPONSE(part1, part2, point1, point2, normal, I, J)

use param_phys
use ellipsoid_particle
use mod_quaternion

implicit none

!=================================================!
!=================================================!
type(ELL_PART), intent(inout) :: part1, part2

real(kind=8), dimension(ndim, 1), intent(in) :: point1, point2, normal

integer, intent(in) :: I, J
!=================================================!
!=================================================!
! Impact Arm
real(kind=8), dimension(ndim, 1) :: R_impact1, R_impact2

! Velocity at contact point
real(kind=8), dimension(ndim, 1) :: VPC1, VPC2

! Relative Velocity of Contact Points
real(kind=8), dimension(ndim, 1) :: VREL

real(kind=8) :: VRN

! Velocity Difference 
real(kind=8), dimension(ndim, 1) :: DELV
!=================================================!
!=================================================!
! Inverse Mass Identity Matrix
real(kind=8), dimension(ndim, ndim) :: m_I

! Inertia Matrix and its Inverse
real(kind=8), dimension(ndim, ndim) :: I_p, Inv_Ip!, I_p1, I_p2, Inv_Ip1, Inv_Ip2

! Skew Symmetric Matrix of the Impact Arm vector and its transpose
real(kind=8), dimension(ndim, ndim) :: R_X1, R_X_Tr1, R_X2, R_X_Tr2

! Final mass-Inetria Symmetric Matrix and its Inverse
real(kind=8), dimension(ndim, ndim) :: K1, K2, K, invK

!=================================================!
!=================================================!

! Impulse vector with Inertia
real(kind=8), dimension(ndim, 1) :: PC, K_N, V_IMP
!real(kind=8) :: PCX, PCY, PCZ

! Tangential and Normal Component of Impulse Vector
real(kind=8) :: PCT, PCN, KN_dot_N, PCT_X, PCT_Y, PCT_Z

real(kind=8) :: NRM_TW

real(kind=8), dimension(ndim, 1) :: TW

! Velcity Vector
!real(kind=8), dimension(ndim, 1) :: U_P

! Angular Velocity vector
!real(kind=8), dimension(ndim, 1) :: OMG

real(kind=8) :: E_initial, E_final, E_final_trans, E_final_rot

!=================================================!
!=================================================!
!write(*,*) 'Normal', normal(:,1)

!write(*,*) 'PT1 - PT2', (point1 - point2)/(sqrt((point1(1,1)-point2(1,1))**2 + (point1(2,1)-point2(2,1))**2 + (point1(3,1)-point2(3,1))**2))

! Calculate the velocity at Contact Point
! Particle 1
R_impact1(1,1) = abs(point1(1,1) - part1%XP)
R_impact1(2,1) = abs(point1(2,1) - part1%YP)
R_impact1(3,1) = abs(point1(3,1) - part1%ZP)

VPC1(1,1) = part1%UP + (part1%OMEGAY*R_impact1(3,1) - part1%OMEGAZ*R_impact1(2,1))
VPC1(2,1) = part1%VP + (part1%OMEGAZ*R_impact1(1,1) - part1%OMEGAX*R_impact1(3,1))
VPC1(3,1) = part1%WP + (part1%OMEGAX*R_impact1(2,1) - part1%OMEGAY*R_impact1(1,1))

! Particle 2
R_impact2(1,1) = abs(point2(1,1) - part2%XP)
R_impact2(2,1) = abs(point2(2,1) - part2%YP)
R_impact2(3,1) = abs(point2(3,1) - part2%ZP)              

VPC2(1,1) = part2%UP + (part2%OMEGAY*R_impact2(3,1) - part2%OMEGAZ*R_impact2(2,1))
VPC2(2,1) = part2%VP + (part2%OMEGAZ*R_impact2(1,1) - part2%OMEGAX*R_impact2(3,1))
VPC2(3,1) = part2%WP + (part2%OMEGAX*R_impact2(2,1) - part2%OMEGAY*R_impact2(1,1))


! Relative Velocity of Contact Points
VREL = VPC2 - VPC1

! Dot Product of Relative Velocity and normal vector
VRN = VREL(1,1)*normal(1,1) + VREL(2,1)*normal(2,1) + VREL(3,1)*normal(3,1)

! Velocity Difference DelV
DELV(:,1) = VREL(:,1) + EW*VRN*normal(:,1)

!=================================================!
!================================================================!
! Calculation of mass-Inertia Positive Definite Symmetric Matrix !
!================================================================!
! Inverse Mass Identity Matrix
m_I(:,:) = 0.0

m_I(1,1) = 1.0
m_I(2,2) = 1.0
m_I(3,3) = 1.0

! Inertia matrix 
I_p(:,:) = 0.0
Inv_Ip(:,:) = 0.0

I_p(1,1) = IPXX
I_p(2,2) = IPYY
I_p(3,3) = IPZZ

!I_p1 = I_p
!I_p2 = I_p 
!write(*,*) 'I_p' I_p(1,:)

! Transform Inertia Matrix to the global coordinate system
! We don't do this at the moment
!call transform_basis(I_p1, part1%ELLQUAT, shape(I_p1))

!call transform_basis(I_p2, part1%ELLQUAT, shape(I_p2))

!write(*,*) 'I_p after transform'

! Invert the Transformed Inertia Matrix
!call invert_ndim3_matrix(I_p1, Inv_Ip1)

!call invert_ndim3_matrix(I_p2, Inv_Ip2)

call invert_ndim3_matrix(I_p, Inv_Ip)
!write(*,*) 'Inverse of I_p'

!!==== Skew Symmetric Matrix of the Impact Arm vector ====!!
! Particle 1
R_X1(1,1) =  0.0
R_X1(1,2) = -R_impact1(3,1)
R_X1(1,3) =  R_impact1(2,1)

R_X1(2,1) =  R_impact1(3,1) 
R_X1(2,2) =  0.0
R_X1(2,3) = -R_impact1(1,1)

R_X1(3,1) = -R_impact1(2,1)
R_X1(3,2) =  R_impact1(1,1)
R_X1(3,3) =  0.0

R_X_Tr1 = transpose(R_X1)

!write(*,*) 'R_X'

!write(*,*) 'transpose of R_X'

K1 = m_I + matmul(R_X_Tr1, matmul(Inv_Ip, R_X1))

!write(*,*) 'K1'

!write(*,*) 'Inverse of K'

! Particle 2
R_X2(1,1) =  0.0
R_X2(1,2) = -R_impact2(3,1)
R_X2(1,3) =  R_impact2(2,1)

R_X2(2,1) =  R_impact2(3,1) 
R_X2(2,2) =  0.0
R_X2(2,3) = -R_impact2(1,1)

R_X2(3,1) = -R_impact2(2,1)
R_X2(3,2) =  R_impact2(1,1)
R_X2(3,3) =  0.0

R_X_Tr2 = transpose(R_X2)

!write(*,*) 'R_X'

!write(*,*) 'transpose of R_X'

K2 = m_I + matmul(R_X_Tr2, matmul(Inv_Ip, R_X2))

!write(*,*) 'K2'

K = K1 + K2

call invert_ndim3_matrix(K, invK)

!write(*,*) 'K'

!write(*,*) 'Inverse of K'

!======================================================================!
!                       Calculation of Impulse                         !
!======================================================================!
!! Impulse for sticking
PC = - matmul(invK, DelV) 

!! Check Friction Cone
! Normal Impulse 
PCN = PC(1,1)*normal(1,1) + PC(2,1)*normal(2,1) + PC(3,1)*normal(3,1)
!write(*,*) PCN

! Tangential Impulse
PCT_X = PC(1,1) - PCN*normal(1,1)
PCT_Y = PC(2,1) - PCN*normal(2,1)
PCT_Z = PC(3,1) - PCN*normal(3,1)

PCT = sqrt(PCT_X*PCT_X + PCT_Y*PCT_Y + PCT_Z*PCT_Z)
!write(*,*) PCT


if (PCT > MUW*PCN) then

    ! Tangential unit vector
    TW(1,1) = VREL(1,1) - VRN*normal(1,1)
    TW(2,1) = VREL(2,1) - VRN*normal(2,1)
    TW(3,1) = VREL(3,1) - VRN*normal(3,1)

    NRM_TW = sqrt(TW(1,1)*TW(1,1) + TW(2,1)*TW(2,1) + TW(3,1)*TW(3,1))

    if (NRM_TW /= 0.0) then

        TW(1,1) = -TW(1,1)/NRM_TW
        TW(2,1) = -TW(2,1)/NRM_TW
        TW(3,1) = -TW(3,1)/NRM_TW

    end if

    !! Impulse for sliding
    ! Calculate (n + muw*t)
    V_IMP = normal + MUW*TW 

    ! Calculate K.(n + muw*t)
    K_N = matmul(K, V_IMP)

    ! Dot product of K.(n + muw*t) and n
    KN_dot_N = K_N(1,1)*normal(1,1) + K_N(2,1)*normal(2,1) + K_N(3,1)*normal(3,1)

    ! Dot product of DelV and n
    PCN = -1.0*(DelV(1,1)*normal(1,1) + DelV(2,1)*normal(2,1) + DelV(3,1)*normal(3,1))

    ! Final expression of Impulse
    PC = (PCN/KN_dot_N)*(V_IMP)

    !write(*,*) '       '
    !write(*,*) 'Checked slip condition'
    !write(*,*) 'Norm of Tangential Velocity :', NRM_TW
    !write(*,*) 'Tangential Vector :', TW
    !write(*,*) '       '

end if
!write(*,*) 'Impulse vector', PC

!----------------------------------------------------------------!
!############################################################!
!################# ENERGY CALCULATION #######################!
!############################################################!
E_initial = 0.5*(part1%UP*part1%UP + part1%VP*part1%VP + part1%WP*part1%WP &
               + IPXX*(part1%OMEGAX*part1%OMEGAX)  &
               + IPYY*(part1%OMEGAY*part1%OMEGAY)  &
               + IPZZ*(part1%OMEGAZ*part1%OMEGAZ)  &
               + part2%UP*part2%UP + part2%VP*part2%VP + part2%WP*part2%WP &
               + IPXX*(part2%OMEGAX*part2%OMEGAX)  &
               + IPYY*(part2%OMEGAY*part2%OMEGAY)  &
               + IPZZ*(part2%OMEGAZ*part2%OMEGAZ))
!----------------------------------------------------------------!

!============================================================!
!=================== Update the Velocity ====================!
!============================================================!
! Particle 1
part1%UP = part1%UP - PC(1,1)
part1%VP = part1%VP - PC(2,1)
part1%WP = part1%WP - PC(3,1)

part1%OMEGAX = part1%OMEGAX + (R_impact1(2,1)*(-PC(3,1)) - R_impact1(3,1)*(-PC(2,1)))/IPXX
part1%OMEGAY = part1%OMEGAY + (R_impact1(3,1)*(-PC(1,1)) - R_impact1(1,1)*(-PC(3,1)))/IPYY
part1%OMEGAZ = part1%OMEGAZ + (R_impact1(1,1)*(-PC(2,1)) - R_impact1(2,1)*(-PC(1,1)))/IPZZ


! Particle 2
part2%UP = part2%UP + PC(1,1)
part2%VP = part2%VP + PC(2,1)
part2%WP = part2%WP + PC(3,1)

part2%OMEGAX = part2%OMEGAX + (R_impact2(2,1)*PC(3,1) - R_impact2(3,1)*PC(2,1))/IPXX
part2%OMEGAY = part2%OMEGAY + (R_impact2(3,1)*PC(1,1) - R_impact2(1,1)*PC(3,1))/IPYY
part2%OMEGAZ = part2%OMEGAZ + (R_impact2(1,1)*PC(2,1) - R_impact2(2,1)*PC(1,1))/IPZZ

!============================================================!
!============================================================!
!============================================================!

!############################################################!
!################# ENERGY CALCULATION #######################!
!############################################################!
E_final = 0.5*(part1%UP*part1%UP + part1%VP*part1%VP + part1%WP*part1%WP &
             + IPXX*(part1%OMEGAX*part1%OMEGAX)  &
             + IPYY*(part1%OMEGAY*part1%OMEGAY)  &
             + IPZZ*(part1%OMEGAZ*part1%OMEGAZ)) &
             + 0.5*(part2%UP*part2%UP + part2%VP*part2%VP + part2%WP*part2%WP &
             + IPXX*(part2%OMEGAX*part2%OMEGAX)  &
             + IPYY*(part2%OMEGAY*part2%OMEGAY)  &
             + IPZZ*(part2%OMEGAZ*part2%OMEGAZ))


E_final_trans = 0.5*(part1%UP*part1%UP + part1%VP*part1%VP + part1%WP*part1%WP &
                   + part2%UP*part2%UP + part2%VP*part2%VP + part2%WP*part2%WP) 

E_final_rot = 0.5*(IPXX*(part1%OMEGAX*part1%OMEGAX)  &
                 + IPYY*(part1%OMEGAY*part1%OMEGAY)  &
                 + IPZZ*(part1%OMEGAZ*part1%OMEGAZ)  &
                 + IPXX*(part2%OMEGAX*part2%OMEGAX)  &
                 + IPYY*(part2%OMEGAY*part2%OMEGAY)  &
                 + IPZZ*(part2%OMEGAZ*part2%OMEGAZ))


if(abs((E_initial- E_final)/E_initial) > 1.0e-9) then
    ERROR_COUNT = ERROR_COUNT + 1

    AVG_ERROR = AVG_ERROR + abs((E_initial- E_final)/E_initial)
    !write(*,*) '   '
    !write(*,*) ' Error in Particle Collision '
    !write(*,*) ' BALANCE OF KINETIC ENERGY'
    !write(*,*) '-----------------------------------------------------------'    
    !print*, ' KE Initial', E_initial
    !print*, ' KE Final  ', E_final
    !print*, ' Error in KE (Total)%', abs((E_initial- E_final)/E_initial)*100.0
    !write(*,*) '-----------------------------------------------------------'
    !print*, ' KE Final Translational', E_final_trans
    !print*, ' KE Final Rotational', E_final_rot
    !print*, ' KE Final Total', E_final_trans + E_final_rot
    !write(*,*) '-----------------------------------------------------------'
    !write(*,*) 'Particle 1 information'
    !write(*,*) I, part1%OVERLAP
    !write(*,*) 'X = ', part1%XP, 'Y = ', part1%YP, ' Z =', part1%ZP
    !write(*,*) 'Particle 2 information'
    !write(*,*) J, part2%OVERLAP
    !write(*,*) 'X = ', part2%XP, 'Y = ', part2%YP, ' Z =', part2%ZP
    !write(*,*) 'Penetration =', sqrt((point1(1,1) - point2(1,1))**2 + (point1(2,1) - point2(2,1))**2 + (point1(3,1) - point2(3,1))**2)
    !write(*,*) 'Common normal', normal
    !write(*,*) 'Point 1', point1
    !write(*,*) 'Point 2', point2
    !write(*,*) 'Point1 - Point2', point1 - point2

    !write(*, *) 'Collision count =',  COLLISION_COUNT 
    !stop
end if 
!stop
end subroutine COLLISION_RESPONSE