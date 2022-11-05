module mod_quaternion

implicit none


! Define the type quaternion 
type quaternion
    real(kind=8) :: a, b, c, d 
end type quaternion

! Initialize the public interface for operations
public :: norm_q, neg_q, conj_q, unit_quat
public :: add_quat
public :: mult_quat
public :: vec_mult_quat

! Initialize an instance of quaternion
!type(quaternion), dimension(:), allocatable :: Quat 

! Quaternion opertations (private to the module)
private :: q_plus_q, q_plus_r, &
           q_mult_q, q_mult_r, &
           q_plus_v, quat_prod_vec, &
           vec_prod_quat, &
           norm, neg, conj, &
           normalize_quat

! Link the public operations to private operations
!! Norm of the quaternion
interface norm_q
    module procedure norm
end interface

!! Negatvie of a quaternion
interface neg_q
    module procedure neg
end interface

!! Conjugate of a quaternion
interface conj_q
    module procedure conj
end interface

!! Addition of Quaternion
interface add_quat
    module procedure q_plus_q, q_plus_r, q_plus_v
end interface

!! Multiplication of Quaternion
interface mult_quat
    module procedure q_mult_q, q_mult_r
end interface

!! Normalization of a Quaternion to get unit quaternion
interface unit_quat
    module procedure normalize_quat
end interface

!! Multiplication of a quaternion and a vector 
interface vec_mult_quat
    module procedure quat_prod_vec, vec_prod_quat
end interface


contains

!!=========================================================!!
!!              Definition of Private Functions            !! 
!!=========================================================!!

!! Function def: Takes in Quaternion returns Norm
function norm(x) result(res)
    real(kind=8) :: res
    type (quaternion), intent(in) :: x

    res = sqrt(x%a*x%a + x%b*x%b + x%c*x%c + x%d*x%d)

end function norm

!! Function def: Takes in Quaternion returns Negative of the quaternion
function neg(x) result(res)
    type (quaternion) :: res
    type (quaternion), intent(in) :: x

    res%a = -x%a
    res%b = -x%b
    res%c = -x%c
    res%d = -x%d

end function neg

!! Function def: Takes in Quaternion returns Conjugate of the quaternion
function conj(x) result(res)
    type (quaternion) :: res
    type (quaternion), intent(in) :: x

    res%a = x%a
    res%b = -x%b
    res%c = -x%c
    res%d = -x%d

end function conj

!! Function def: Takes in two Quaternions returns the sum of the two
function q_plus_q(x,y) result(res)
    type (quaternion) :: res
    type (quaternion), intent(in) :: x, y

    res%a = x%a + y%a
    res%b = x%b + y%b
    res%c = x%c + y%c
    res%d = x%d + y%d

end function q_plus_q

!! Function def: Takes in a Quaternion and a scalar returns the sum 
!! of the two in the form of a Quaternion
function q_plus_r(x,r) result(res)
    type (quaternion) :: res
    type (quaternion), intent(in) :: x
    real(kind=8), intent(in) :: r

    res%a = x%a + r
    res%b = x%b
    res%c = x%c
    res%d = x%d    

end function q_plus_r

!! Function def: Takes in two Quaternions returns the product of the two
function q_mult_q(x,y) result(res)
    type (quaternion) :: res
    type (quaternion), intent(in) :: x, y

   res%a = x%a*y%a - x%b*y%b - x%c*y%c - x%d*y%d
   res%b = x%a*y%b + x%b*y%a + x%c*y%d - x%d*y%c
   res%c = x%a*y%c - x%b*y%d + x%c*y%a + x%d*y%b
   res%d = x%a*y%d + x%b*y%c - x%c*y%b + x%d*y%a

end function q_mult_q

!! Function def: Takes in a Quaternion and a scalar returns the product 
!! of the two in the form of a Quaternion
function q_mult_r(x,r) result(res)
    type (quaternion) :: res
    type (quaternion), intent(in) :: x
    real(kind=8), intent(in) :: r

    res%a = x%a*r
    res%b = x%b*r
    res%c = x%c*r
    res%d = x%d*r

end function q_mult_r


!! Function def: Takes in Quaternion returns 
!! the normalized Quaternion with norm =1
function normalize_quat(x) result(res)
    type (quaternion) :: res
    type (quaternion), intent(in) :: x

    real :: norm_

    norm_ = norm(x)

    res%a = x%a/norm_
    res%b = x%b/norm_
    res%c = x%c/norm_
    res%d = x%d/norm_

end function normalize_quat


!! Function def: Takes in a Quaternion and a scalar returns the sum 
!! of the two in the form of a Quaternion
function q_plus_v(x,v) result(res)
    type (quaternion) :: res
    type (quaternion), intent(in) :: x
    real(kind=8), dimension(3,1), intent(in) :: v

    res%a = x%a
    res%b = x%b + v(1,1)
    res%c = x%c + v(2,1)
    res%d = x%d + v(3,1)

end function q_plus_v

!! Function def: Takes in a vector and a Quaternion returns the product 
!! of the two in the form of a Quaternion
function vec_prod_quat(vec,q) result(res)
    type (quaternion) :: res
    real(kind=8), dimension(3,1), intent(in) :: vec
    type (quaternion), intent(in) :: q
    type (quaternion) :: y

    y%a = 0.0
    y%b = vec(1,1)
    y%c = vec(2,1)
    y%d = vec(3,1)
 
   res%a = q%a*y%a - q%b*y%b - q%c*y%c - q%d*y%d
   res%b = q%a*y%b + q%b*y%a + q%c*y%d - q%d*y%c
   res%c = q%a*y%c - q%b*y%d + q%c*y%a + q%d*y%b
   res%d = q%a*y%d + q%b*y%c - q%c*y%b + q%d*y%a

end function vec_prod_quat

!! Function def: Takes in a Quaternion and a vector returns the product 
!! of the two in the form of a vector
function quat_prod_vec(x,v) result(res)
    real(kind=8), dimension(3,1) :: res
    type (quaternion), intent(in) :: x
    real(kind=8), dimension(3,1), intent(in) :: v
    type (quaternion) :: y, z


    y%a = 0.0
    y%b = v(1,1)
    y%c = v(2,1)
    y%d = v(3,1)


    z%a = x%a*y%a - x%b*y%b - x%c*y%c - x%d*y%d
    z%b = x%a*y%b + x%b*y%a + x%c*y%d - x%d*y%c
    z%c = x%a*y%c - x%b*y%d + x%c*y%a + x%d*y%b
    z%d = x%a*y%d + x%b*y%c - x%c*y%b + x%d*y%a

    if (z%a == 0.0) then

        res(1,1) = z%b
        res(2,1) = z%c
        res(3,1) = z%d

    else

        write(*,*) "Error"
        stop

    end if

end function quat_prod_vec


end module mod_quaternion