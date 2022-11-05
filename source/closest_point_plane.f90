subroutine closest_point_plane(XP, YP, ZP, &
                              quat, WPOS, WNORMAL, &
                              plane_vec, margin)

use param_phys
use mod_quaternion
use ellipsoid_particle

implicit none

!==============================================================!    
real(kind=8), intent(in) :: XP, YP, ZP

type(quaternion), intent(in) :: quat

real(kind=8), dimension(ndim, 1), intent(in) :: WPOS, WNORMAL

real(kind=8), dimension(ndim, 1), intent(inout) :: plane_vec

real(kind=8), intent(inout) :: margin
!==============================================================!
! Ellipsoid Center Position Vector
real(kind=8), dimension(ndim, 1) :: ELLC 

! Ellipsoid Matrix and it inverse
real(kind=8), dimension(ndim, ndim) :: ELL, inv_ELL

! Transpose of wall normal vector
real(kind=8), dimension(1, ndim) :: trWNORMAL

! Product of inversed ELL and wall normal
real(kind=8), dimension(ndim, 1) :: inv_ELL_N

! Numerator & Denominator of Lagrange multiplier
real(kind=8), dimension(1,1) :: trN_wp_C, trN_invELL_N

! Lagrange Multiplier 
real(kind=8) :: lagrange

! Margin function
real(kind=8), dimension(ndim,1) :: vec

real(kind=8), dimension(ndim, 1) :: ELL_

real(kind=8), dimension(1, 1) :: margin_vec
!==============================================================!
! Initialize the Ellipsoid Center
ELLC(1,1) = XP
ELLC(2,1) = YP
ELLC(3,1) = ZP

! Initialize the ellipsoid matrix
ELL(:,:) = 0.0
ELL(1,1) = 1.0/(ELL_A*ELL_A)
ELL(2,2) = 1.0/(ELL_B*ELL_B)
ELL(3,3) = 1.0/(ELL_C*ELL_C)

! Transform the Ellipsoid Matrix acc. to Orientation
call transform_basis(ELL, quat, shape(ELL))

! Invert ELL
call invert_ndim3_matrix(ELL, inv_ELL)

! Transpose of wall normal vector
trWNORMAL = transpose(WNORMAL)

! [ELL^(-1)]N
inv_ELL_N = matmul(inv_ELL, WNORMAL)

! trN [ELL^(-1)] N
trN_invELL_N = matmul(trWNORMAL, inv_ELL_N)

! N^T * (WPOS - C)
trN_wp_C = matmul(trWNORMAL, (WPOS - ELLC))

lagrange = trN_wp_C(1,1)/trN_invELL_N(1,1)

plane_vec = ELLC + lagrange*inv_ELL_N  


!!=== Calculate margin function ===!!
vec = plane_vec - ELLC

ELL_ = matmul(ELL, vec)


margin_vec = matmul(transpose(vec), ELL_)


margin = margin_vec(1,1)

!if(margin < 1.0) then
!    write(*,*) 'Margin',margin
!    write(*,*) 'Closest Point', plane_vec
!    write(*,*) 'Ellipsoid Center Position', ELLC
!    write(*,*) ELL_(:,1)
!    write(*,*) 'Difference between closest point & center', vec
!    stop
!end if



end subroutine closest_point_plane