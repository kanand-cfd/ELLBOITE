subroutine ell_full2eig( n, A, U, lam )

!  A is an n x n PSD matrix with eigen-decomposition: 
!  A = U * lam^2 * U^T.  This routine returns U and lam.

implicit none
!integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n
real(kind=8), intent(in)  :: a(n,n)
real(kind=8), intent(out) :: u(n,n), lam(n)

real(kind=8) :: vt(n,n), work(5*n*n+20*n)
integer    :: lwork, info, i

lwork = 5*n*n+20*n

vt = a
call dgesvd( 'A', 'N', n, n, vt, n, lam, u, n, vt, n, work, lwork, info)

if( info /= 0 ) then
   write(0,*)'ell_full2eig: info, lwork, work(1) = ', info, lwork, work(1)
   stop
endif

do i = 1, n
   lam(i) = sqrt( lam(i) )
end do

return
end subroutine ell_full2eig
