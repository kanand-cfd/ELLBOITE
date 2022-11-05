subroutine WALL_ELLIPSOID_REBOUND(u, v, w, &
                                  omegax, omegay, omegaz, &
                                  q_p, R_impact, norm_W)
use ellipsoid_particle
use param_phys
use mod_quaternion

implicit none

!========================================================!
! 					Input Variables 					 !
!========================================================!

! Translational velocity 
real(kind=8), intent(inout) :: u, v, w 

! Angular velocity
real(kind=8), intent(inout) :: omegax, omegay, omegaz

! Orientation of the particle
type(quaternion), intent(in) :: q_p

! Distance beTWeen particle center and contact point
real(kind=8), dimension(ndim,1), intent(in) :: R_Impact

! Wall Normal Vector
real(kind=8), dimension(ndim,1), intent(in) :: norm_W

!========================================================!
! 					Local Variables 					 !
!========================================================!
! Impact arm vector
!real(kind=8) :: RX, RY, RZ

! Normal vector at Contact Point
!real(kind=8) :: NWX, NWY, NWZ
 
!real(kind=8), dimension(ndim, 1) :: NW

! Relative Velocity
real(kind=8) :: VRX, VRY, VRZ

! Velocity at Contact Point
real(kind=8) :: VPCX, VPCY, VPCZ

! Tangential unit vector
!real(kind=8) :: TWX, TWY, TWZ
real(kind=8) :: NRM_TW

real(kind=8), dimension(ndim, 1) :: TW

! Dot product of relative velocity and normal vector
real(kind=8) :: VRN

! Dot product of relative velocity and tangential vector 
!real(kind=8) :: VRT

!========================================================!
! Inverse Mass Identity Matrix
real(kind=8), dimension(ndim, ndim) :: m_I

! Inertia Matrix and its Inverse
real(kind=8), dimension(ndim, ndim) :: I_p, Inv_Ip

! Skew Symmetric Matrix of the Impact Arm vector and its transpose
real(kind=8), dimension(ndim, ndim) :: R_X, R_X_Tr

! Final mass-Inetria Symmetric Matrix and its Inverse
real(kind=8), dimension(ndim, ndim) :: K, K_1

!========================================================!
! Impulse vector without considering Inertia
!real(kind=8), dimension(ndim, 1) :: PC_Ip

! Impulse vector with Inertia
real(kind=8), dimension(ndim, 1) :: PC, K_N, V_IMP
!real(kind=8) :: PCX, PCY, PCZ

! Tangential and Normal Component of Impulse Vector
real(kind=8) :: PCT, PCN, KN_dot_N, PCT_X, PCT_Y, PCT_Z

! Velcity Vector
real(kind=8), dimension(ndim, 1) :: U_P, DelV

! Angular Velocity vector
real(kind=8), dimension(ndim, 1) :: OMG

real(kind=8) :: E_initial, E_final, E_final_trans, E_final_rot
!========================================================!
!========================================================!
! Initialize the velocity and Angular velocity vectors
U_P(1,1) = u
U_P(2,1) = v
U_P(3,1) = w 

OMG(1,1) = omegax
OMG(2,1) = omegay
OMG(3,1) = omegaz


! Velocity of contact point 
VPCX = U_P(1,1) + (OMG(2,1)*R_Impact(3,1) - OMG(3,1)*R_Impact(2,1))
VPCY = U_P(2,1) + (OMG(3,1)*R_Impact(1,1) - OMG(1,1)*R_Impact(3,1))
VPCZ = U_P(3,1) + (OMG(1,1)*R_Impact(2,1) - OMG(2,1)*R_Impact(1,1))

!write(*,*) 'VPC :', VPCX, VPCY, VPCZ

! Relative velocity
VRX = VPCX
VRY = VPCY
VRZ = VPCZ

! Dot product of relative velocity and normal vector
VRN = VRX*norm_W(1,1) + VRY*norm_W(2,1) + VRZ*norm_W(3,1)

!write(*,*) 'VRN', VRN

!! Relative Velocity before and after collision
DelV(1,1) = VRX + EW*VRN*norm_W(1,1)
DelV(2,1) = VRY + EW*VRN*norm_W(2,1)
DelV(3,1) = VRZ + EW*VRN*norm_W(3,1)


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

!write(*,*) 'I_p' I_p(1,:)

! Transform Inertia Matrix to the global coordinate system
! We don't do this at the moment
!call transform_basis(I_p, q_p, shape(I_p))

!write(*,*) 'I_p after transform'

! Invert the Transformed Inertia Matrix
call invert_ndim3_matrix(I_p, Inv_Ip)

!write(*,*) 'Inverse of I_p'

! Skew Symmetric Matrix of the Impact Arm vector
R_X(1,1) =  0.0
R_X(1,2) = -R_Impact(3,1)
R_X(1,3) =  R_Impact(2,1)

R_X(2,1) =  R_Impact(3,1) 
R_X(2,2) =  0.0
R_X(2,3) = -R_Impact(1,1)

R_X(3,1) = -R_Impact(2,1)
R_X(3,2) =  R_Impact(1,1)
R_X(3,3) =  0.0

R_X_Tr = transpose(R_X)

!write(*,*) 'R_X'

!write(*,*) 'transpose of R_X'

K = m_I + matmul(R_X_Tr, matmul(Inv_Ip, R_X))
call invert_ndim3_matrix(K, K_1)

!write(*,*) 'K'

!write(*,*) 'Inverse of K'
!======================================================================!
!======================================================================!

!======================================================================!
!                       Calculation of Impulse                         !
!======================================================================!
!! Impulse for sticking
PC = - matmul(K_1, DelV) 

!! Check Friction Cone
! Normal Impulse 
PCN = PC(1,1)*norm_W(1,1) + PC(2,1)*norm_W(2,1) + PC(3,1)*norm_W(3,1)
!write(*,*) PCN

! Tangential Impulse
PCT_X = PC(1,1) - PCN*norm_W(1,1)
PCT_Y = PC(2,1) - PCN*norm_W(2,1)
PCT_Z = PC(3,1) - PCN*norm_W(3,1)

PCT = sqrt(PCT_X*PCT_X + PCT_Y*PCT_Y + PCT_Z*PCT_Z)
!write(*,*) PCT


if (PCT > MUW*PCN) then

    ! Tangential unit vector
    TW(1,1) = VRX - VRN*norm_W(1,1)
    TW(2,1) = VRY - VRN*norm_W(2,1)
    TW(3,1) = VRZ - VRN*norm_W(3,1)

    NRM_TW = sqrt(TW(1,1)*TW(1,1) + TW(2,1)*TW(2,1) + TW(3,1)*TW(3,1))

    if (NRM_TW /= 0.0) then

        TW(1,1) = -TW(1,1)/NRM_TW
        TW(2,1) = -TW(2,1)/NRM_TW
        TW(3,1) = -TW(3,1)/NRM_TW

    end if

    !! Impulse for sliding
    ! Calculate (n + muw*t)
    V_IMP = norm_W + MUW*TW 

    ! Calculate K.(n + muw*t)
    K_N = matmul(K, V_IMP)

    ! Dot product of K.(n + muw*t) and n
    KN_dot_N = K_N(1,1)*norm_W(1,1) + K_N(2,1)*norm_W(2,1) + K_N(3,1)*norm_W(3,1)


    ! Dot product of DelV and n
    PCN = -1.0*(DelV(1,1)*norm_W(1,1) + DelV(2,1)*norm_W(2,1) + DelV(3,1)*norm_W(3,1))


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
   E_initial = 0.5*(u*u + v*v + w*w &
               + IPXX*omegax*omegax &
               + IPYY*omegay*omegay &
               + IPZZ*omegaz*omegaz)

!============================================================!
!=================== Update the Velocity ====================!
!============================================================!
u = u + PC(1,1)
v = v + PC(2,1)
w = w + PC(3,1)

omegax = omegax + (R_Impact(2,1)*PC(3,1) - R_Impact(3,1)*PC(2,1))/IPXX
omegay = omegay + (R_Impact(3,1)*PC(1,1) - R_Impact(1,1)*PC(3,1))/IPYY
omegaz = omegaz + (R_Impact(1,1)*PC(2,1) - R_Impact(2,1)*PC(1,1))/IPZZ

!############################################################!
!################# ENERGY CALCULATION #######################!
!############################################################!
E_final = 0.5*( u*u + v*v + w*w &
               + IPXX*omegax*omegax &
               + IPYY*omegay*omegay &
               + IPZZ*omegaz*omegaz)


E_final_trans = 0.5*(u*u + v*v + w*w)

E_final_rot = 0.5*(IPXX*omegax*omegax &
                 + IPYY*omegay*omegay & 
                 + IPZZ*omegaz*omegaz)

WALL_COUNT = WALL_COUNT + 1

if(abs((E_initial- E_final)/E_initial)*100.0 > 1.0e-9) then
    write(*,*) ' Wall Collision '
    write(*,*) ' BALANCE OF KINETIC ENERGY'
    write(*,*) '-----------------------------------------------------------'    
    print*, ' KE Initial', E_initial
    print*, ' KE Final  ', E_final
    print*, ' Error in KE (Total)%', abs((E_initial- E_final)/E_initial)*100.0
    write(*,*) '-----------------------------------------------------------'
    print*, ' KE Final Translational', E_final_trans
    print*, ' KE Final Rotational', E_final_rot
    print*, ' KE Final Total', E_final_trans + E_final_rot
    write(*,*) '-----------------------------------------------------------'
    stop
end if 

end subroutine WALL_ELLIPSOID_REBOUND