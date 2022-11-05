subroutine EULER_INTEGRATION_TEST(omgx, omgy, omgz, qt)

use param_phys
use mod_quaternion
use ellipsoid_particle

implicit none

!=========== INPUT VARIABLES =================!

real(kind=8), intent(inout) :: omgx, omgy, omgz

type(quaternion), intent(in) :: qt

!=============================================!

real(kind=8) :: dt

!!============== RK4 Method =================!!
real(kind=8) :: k1_x, k2_x, k3_x, k4_x
real(kind=8) :: k1_y, k2_y, k3_y, k4_y
real(kind=8) :: k1_z, k2_z, k3_z, k4_z

real(kind=8) :: w0_x, w0_y, w0_z
real(kind=8) :: w1_x, w1_y, w1_z
!=============================================!
real(kind=8), dimension(ndim,1) :: K1, K2, K3, K4

real(kind=8), dimension(ndim,1) :: OMEGA_0
real(kind=8), dimension(ndim,1) :: OMEGA_1

! Moment of inertia and its inverse
real(kind=8), dimension(ndim, ndim) :: MOI, inv_MOI

! Anti-symmetric matrix of omega
real(kind=8), dimension(ndim, ndim) :: MAT_OMG

real(kind=8), dimension(ndim, ndim) :: I_OMG            
!=============================================!

dt = DT_INIT

! Preliminary steps
MOI(1,1) = IPXX
MOI(2,2) = IPYY
MOI(3,3) = IPZZ

! Transform the Inertia Matrix acc. to Orientation
call transform_basis(MOI, qt, shape(MOI))

! Invert the inertia matrix
call invert_ndim3_matrix(MOI, inv_MOI)

! Initialise for the RK4 step

OMEGA_0(1,1) = omgx
OMEGA_0(2,1) = omgy
OMEGA_0(3,1) = omgz

!!====================== 1st step ==========================!!
! Initially Populate the Anti-symmetric matrix of omega
MAT_OMG(:,:) = 0.0

MAT_OMG(1,2) = - OMEGA_0(3,1)
MAT_OMG(1,3) =   OMEGA_0(2,1)

MAT_OMG(2,1) =   OMEGA_0(3,1)
MAT_OMG(2,3) = - OMEGA_0(1,1)

MAT_OMG(3,1) = - OMEGA_0(2,1)
MAT_OMG(3,2) =   OMEGA_0(1,1)

I_OMG = matmul(inv_MOI, matmul(MAT_OMG, MOI))

K1 = dt*(matmul(I_OMG, OMEGA_0))

!!====================== 2nd step ==========================!!
MAT_OMG(:,:) = 0.0

MAT_OMG(1,2) = - OMEGA_0(3,1) + 0.5*K1(3,1)
MAT_OMG(1,3) =   OMEGA_0(2,1) + 0.5*K1(2,1)

MAT_OMG(2,1) =   OMEGA_0(3,1) + 0.5*K1(3,1)
MAT_OMG(2,3) = - OMEGA_0(1,1) + 0.5*K1(1,1)

MAT_OMG(3,1) = - OMEGA_0(2,1) + 0.5*K1(2,1)
MAT_OMG(3,2) =   OMEGA_0(1,1) + 0.5*K1(1,1)


I_OMG = matmul(inv_MOI, matmul(MAT_OMG, MOI))

K2 = dt*(matmul(I_OMG, OMEGA_0 + 0.5*K1))

!!====================== 3rd step ==========================!!
MAT_OMG(:,:) = 0.0

MAT_OMG(1,2) = - OMEGA_0(3,1) + 0.5*K2(3,1)
MAT_OMG(1,3) =   OMEGA_0(2,1) + 0.5*K2(2,1)

MAT_OMG(2,1) =   OMEGA_0(3,1) + 0.5*K2(3,1)
MAT_OMG(2,3) = - OMEGA_0(1,1) + 0.5*K2(1,1)

MAT_OMG(3,1) = - OMEGA_0(2,1) + 0.5*K2(2,1)
MAT_OMG(3,2) =   OMEGA_0(1,1) + 0.5*K2(1,1)


I_OMG = matmul(inv_MOI, matmul(MAT_OMG, MOI))

K3 = dt*(matmul(I_OMG, OMEGA_0 + 0.5*K2))

!!====================== 4th step ==========================!!
MAT_OMG(:,:) = 0.0

MAT_OMG(1,2) = - OMEGA_0(3,1) + K3(3,1)
MAT_OMG(1,3) =   OMEGA_0(2,1) + K3(2,1)

MAT_OMG(2,1) =   OMEGA_0(3,1) + K3(3,1)
MAT_OMG(2,3) = - OMEGA_0(1,1) + K3(1,1)

MAT_OMG(3,1) = - OMEGA_0(2,1) + K3(2,1)
MAT_OMG(3,2) =   OMEGA_0(1,1) + K3(1,1)


I_OMG = matmul(inv_MOI, matmul(MAT_OMG, MOI))

K4 = dt*(matmul(I_OMG, OMEGA_0 + K3))

!!==========================================================!!

OMEGA_1 = OMEGA_0 + (K1 + 2.0*K2 + 2.0*K3 + K4)/6.0


omgx = OMEGA_1(1,1)
omgy = OMEGA_1(2,1)
omgz = OMEGA_1(3,1)

end subroutine EULER_INTEGRATION_TEST