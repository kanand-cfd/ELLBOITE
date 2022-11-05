subroutine ELLIPSOID_OVERLAP_DETECTION(IXP, IYP, IZP, IQUAT, &
                                       JXP, JYP, JZP, JQUAT, &
                                       flag)

use param_phys
use ellipsoid_particle
use mod_quaternion

implicit none

integer, parameter :: n=3, nl = (n*(n+1))/2

real(kind=8), intent(in) :: IXP, IYP, IZP
real(kind=8), intent(in) :: JXP, JYP, JZP

type(quaternion), intent(in) :: IQUAT, JQUAT

logical, intent(inout) :: flag

!-------------------------------------------------------------!

integer    :: i, lu=1, max_it
real(kind=8) :: a1(n,n), a2(n,n), c1(n), c2(n), u1(n,n), u2(n,n), lam1(n), lam2(n), &
              L1(nl), L2(nl), g1(nl), g2(nl)
              
real(kind=8) :: u(n,n), ut(n,n), eye(n,n), xh(n), v(n), qual, dist_center
!logical    :: intersect 

!open( lu, file = 'f.dat' )
!eye = 0.d0
!do i = 1, n
!   eye(i,i) = 1.d0
!end do
dist_center = (IXP - JXP)**2 + (IYP - JYP)**2 + (IZP - JZP)**2
dist_center = sqrt(dist_center)
!write(*,*) 'center Distance', dist_center


if(dist_center > (2.5*ELL_A)) then
    !STOP 'Too far'
    flag = .FALSE.
    return
end if


! Ellipsoid 1 center
c1(1) = IXP
c1(2) = IYP
c1(3) = IZP

! Ellipsoid 2 center
c2(1) = JXP
c2(2) = JYP
c2(3) = JZP

!call random_number( u )
!u  = u + eye -0.5d0
!ut = transpose( u )
!a1 = matmul( u, ut )

call calculate_Ellipsoid_matrix(a1, IQUAT)

!call random_number( u )
!u  = u + eye -0.5d0
!ut = transpose( u )
!a2 = matmul( u, ut )

call calculate_Ellipsoid_matrix(a2, JQUAT)

call ell_full2eig( n, a1, u1, lam1 )
call ell_full2eig( n, a2, u2, lam2 )

!if( n== 2 ) write(lu,100) 4.d0, 1.d0, c1, lam1, u1
!if( n== 2 ) write(lu,100) 4.d0, 2.d0, c2, lam2, u2
!100   format((1p,2e20.10))

!  Cholesky representation

call ell_full2low( n, a1, L1 )  ! g contains lower triangle of aa in packed format
call ell_low2chol( n, L1 ,g1 )  ! g contains lower Cholesky triangle in packed format

call ell_full2low( n, a2, L2 )  ! g contains lower triangle of aa in packed format
call ell_low2chol( n, L2 ,g2 )  ! g contains lower Cholesky triangle in packed format

qual   = 0.999
max_it = 50
call ell_pair_separate( n, c1, g1, c2, g2, qual, max_it, xh, v, flag )
!write(0,*)'intersect = ', intersect
!write(0,*)'xh = ', xh
!write(0,*)'v  = ', v
!flag = intersect

return    
end subroutine ELLIPSOID_OVERLAP_DETECTION