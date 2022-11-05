subroutine invert_ndim3_matrix(M, Inv_M)

implicit none

integer, parameter :: ndim = 3

real(kind=8), dimension(ndim, ndim), intent(in) :: M
real(kind=8), dimension(ndim, ndim), intent(out) :: Inv_M

real(kind=8) :: detinv

!========================================================!
! Calculate the Determinant of the Matrix
detinv = 1/(M(1,1)*M(2,2)*M(3,3) - M(1,1)*M(2,3)*M(3,2)&
          - M(1,2)*M(2,1)*M(3,3) + M(1,2)*M(2,3)*M(3,1)&
          + M(1,3)*M(2,1)*M(3,2) - M(1,3)*M(2,2)*M(3,1))


if (detinv == 0) then
    stop 'Matrix is numerically singular!'
end if


! Calculate the inverse of the matrix
    Inv_M(1,1) = +detinv * (M(2,2)*M(3,3) - M(2,3)*M(3,2))
    Inv_M(2,1) = -detinv * (M(2,1)*M(3,3) - M(2,3)*M(3,1))
    Inv_M(3,1) = +detinv * (M(2,1)*M(3,2) - M(2,2)*M(3,1))
    Inv_M(1,2) = -detinv * (M(1,2)*M(3,3) - M(1,3)*M(3,2))
    Inv_M(2,2) = +detinv * (M(1,1)*M(3,3) - M(1,3)*M(3,1))
    Inv_M(3,2) = -detinv * (M(1,1)*M(3,2) - M(1,2)*M(3,1))
    Inv_M(1,3) = +detinv * (M(1,2)*M(2,3) - M(1,3)*M(2,2))
    Inv_M(2,3) = -detinv * (M(1,1)*M(2,3) - M(1,3)*M(2,1))
    Inv_M(3,3) = +detinv * (M(1,1)*M(2,2) - M(1,2)*M(2,1))

!write(*,*) Inv_M(1, :)
!write(*,*) Inv_M(2, :)
!write(*,*) Inv_M(3, :)

!stop
end subroutine invert_ndim3_matrix