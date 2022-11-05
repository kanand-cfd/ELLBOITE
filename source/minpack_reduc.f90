module minpack_module

    use mod_quaternion
    use param_phys
    use ellipsoid_particle
    use iso_fortran_env, only: wp => real64

    implicit none

    real(wp), dimension(3), parameter :: dpmpar = [epsilon(1.0_wp), &
                                                   tiny(1.0_wp), &
                                                   huge(1.0_wp)] !! machine constants

    real(wp), parameter, private :: epsmch = dpmpar(1) !! the machine precision
    real(wp), parameter, private :: one = 1.0_wp
    real(wp), parameter, private :: zero = 0.0_wp

    abstract interface

        subroutine fcn_lmder(m, n, x, fvec, fjac, ldfjac, iflag, &
                                           IXP, IYP, IZP, IQUAT, &
                                           JXP, JYP, JZP, JQUAT)

            !! user-supplied subroutine for [[lmder]] and [[lmder1]]
            import :: wp
            import :: quaternion
            implicit none
            integer, intent(in) :: m !! the number of functions.
            integer, intent(in) :: n !! the number of variables.
            integer, intent(in) :: ldfjac !! leading dimension of the array fjac.
            integer, intent(inout) :: iflag !! if iflag = 1 calculate the functions at x and
                                           !! return this vector in fvec. do not alter fjac.
                                           !! if iflag = 2 calculate the jacobian at x and
                                           !! return this matrix in fjac. do not alter fvec.
                                           !!
                                           !! the value of iflag should not be changed by fcn unless
                                           !! the user wants to terminate execution of lmder.
                                           !! in this case set iflag to a negative integer.
            real(wp), intent(in) :: x(n) !! independent variable vector
            real(wp), intent(inout) :: fvec(m) !! value of function at `x`
            real(wp), intent(inout) :: fjac(ldfjac, n) !! jacobian matrix at `x`

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Particle 1
            real(kind=8), intent(in) :: IXP, IYP, IZP
            type(quaternion), intent(in) :: IQUAT
            ! Particle 2
            real(kind=8), intent(in) :: JXP, JYP, JZP
            type(quaternion), intent(in) :: JQUAT
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        end subroutine fcn_lmder

    end interface

contains
!*****************************************************************************************


!*****************************************************************************************
!>
!  this subroutine checks the gradients of m nonlinear functions
!  in n variables, evaluated at a point x, for consistency with
!  the functions themselves.
!
!  the subroutine does not perform reliably if cancellation or
!  rounding errors cause a severe loss of significance in the
!  evaluation of a function. therefore, none of the components
!  of x should be unusually small (in particular, zero) or any
!  other value which may cause loss of significance.

    subroutine chkder(m, n, x, Fvec, Fjac, Ldfjac, Xp, Fvecp, Mode, Err)

        implicit none

        integer, intent(in) :: m !! a positive integer input variable set to the number
                            !! of functions.
        integer, intent(in) :: n !! a positive integer input variable set to the number
                            !! of variables.
        integer, intent(in) :: Ldfjac !! a positive integer input parameter not less than m
                                    !! which specifies the leading dimension of the array fjac.
        integer, intent(in) :: Mode !! an integer input variable set to 1 on the first call
                                !! and 2 on the second. other values of mode are equivalent
                                !! to mode = 1.
                                !!
                                !! the user must call chkder twice,
                                !! first with mode = 1 and then with mode = 2.
                                !!
                                !!  * mode = 1. **on input**, x must contain the point of evaluation.
                                !!    **on output**, xp is set to a neighboring point.
                                !!
                                !!  * mode = 2. **on input**, fvec must contain the functions and the
                                !!    rows of fjac must contain the gradients
                                !!    of the respective functions each evaluated
                                !!    at x, and fvecp must contain the functions
                                !!    evaluated at xp.
                                !!    **on output**, err contains measures of correctness of
                                !!    the respective gradients.
        real(wp), intent(in) :: x(n) !! input array
        real(wp), intent(in) :: Fvec(m) !! an array of length m. on input when mode = 2,
                                    !! fvec must contain the functions evaluated at x.
        real(wp), intent(in) :: Fjac(Ldfjac, n) !! an m by n array. on input when mode = 2,
                                            !! the rows of fjac must contain the gradients of
                                            !! the respective functions evaluated at x.
        real(wp), intent(out) :: Xp(n) !! an array of length n. on output when mode = 1,
                                    !! xp is set to a neighboring point of x.
        real(wp), intent(in) :: Fvecp(m) !! an array of length m. on input when mode = 2,
                                    !! fvecp must contain the functions evaluated at xp.
        real(wp), intent(out) :: Err(m) !! an array of length m. on output when mode = 2,
                                    !! err contains measures of correctness of the respective
                                    !! gradients. if there is no severe loss of significance,
                                    !! then if err(i) is 1.0 the i-th gradient is correct,
                                    !! while if err(i) is 0.0 the i-th gradient is incorrect.
                                    !! for values of err between 0.0 and 1.0, the categorization
                                    !! is less certain. in general, a value of err(i) greater
                                    !! than 0.5 indicates that the i-th gradient is probably
                                    !! correct, while a value of err(i) less than 0.5 indicates
                                    !! that the i-th gradient is probably incorrect.

        integer :: i, j
        real(wp) :: temp

        real(wp), parameter :: eps = sqrt(epsmch)
        real(wp), parameter :: factor = 100.0_wp
        real(wp), parameter :: epsf = factor*epsmch
        real(wp), parameter :: epslog = log10(eps)

        select case (Mode)
        case (2)
            Err = zero
            do j = 1, n
                temp = abs(x(j))
                if (temp == zero) temp = one
                do i = 1, m
                    Err(i) = Err(i) + temp*Fjac(i, j)
                end do
            end do
            do i = 1, m
                temp = one
                if (Fvec(i) /= zero .and. Fvecp(i) /= zero .and. abs(Fvecp(i) - Fvec(i)) >= epsf*abs(Fvec(i))) &
                    temp = eps*abs((Fvecp(i) - Fvec(i))/eps - Err(i))/(abs(Fvec(i)) + abs(Fvecp(i)))
                Err(i) = one
                if (temp > epsmch .and. temp < eps) Err(i) = (log10(temp) - epslog)/epslog
                if (temp >= eps) Err(i) = zero
            end do
        case (1)
            do j = 1, n
                temp = eps*abs(x(j))
                if (temp == zero) temp = eps
                Xp(j) = x(j) + temp
            end do
        case default
            error stop 'invalid mode in chkder'
        end select

    end subroutine chkder
!*****************************************************************************************


!*****************************************************************************************
!>
!  given an n-vector x, this function calculates the
!  euclidean norm of x.
!
!  the euclidean norm is computed by accumulating the sum of
!  squares in three different sums. the sums of squares for the
!  small and large components are scaled so that no overflows
!  occur. non-destructive underflows are permitted. underflows
!  and overflows do not occur in the computation of the unscaled
!  sum of squares for the intermediate components.
!  the definitions of small, intermediate and large components
!  depend on two constants, rdwarf and rgiant. the main
!  restrictions on these constants are that rdwarf**2 not
!  underflow and rgiant**2 not overflow. the constants
!  given here are suitable for every known computer.

    pure real(wp) function enorm(n, x)

        implicit none

        integer, intent(in) :: n !! a positive integer input variable.
        real(wp), intent(in) :: x(n) !! an input array of length n.

        integer :: i
        real(wp) :: agiant, s1, s2, s3, xabs, x1max, x3max

        real(wp), parameter :: rdwarf = 3.834e-20_wp
        real(wp), parameter :: rgiant = 1.304e19_wp

        s1 = zero
        s2 = zero
        s3 = zero
        x1max = zero
        x3max = zero
        agiant = rgiant/real(n, wp)
        do i = 1, n
            xabs = abs(x(i))
            if (xabs > rdwarf .and. xabs < agiant) then
                ! sum for intermediate components.
                s2 = s2 + xabs**2
            elseif (xabs <= rdwarf) then
                ! sum for small components.
                if (xabs <= x3max) then
                    if (xabs /= zero) s3 = s3 + (xabs/x3max)**2
                else
                    s3 = one + s3*(x3max/xabs)**2
                    x3max = xabs
                end if
                ! sum for large components.
            elseif (xabs <= x1max) then
                s1 = s1 + (xabs/x1max)**2
            else
                s1 = one + s1*(x1max/xabs)**2
                x1max = xabs
            end if
        end do

        ! calculation of norm.

        if (s1 /= zero) then
            enorm = x1max*sqrt(s1 + (s2/x1max)/x1max)
        elseif (s2 == zero) then
            enorm = x3max*sqrt(s3)
        else
            if (s2 >= x3max) enorm = sqrt(s2*(one + (x3max/s2)*(x3max*s3)))
            if (s2 < x3max) enorm = sqrt(x3max*((s2/x3max) + (x3max*s3)))
        end if

    end function enorm
!*****************************************************************************************


!*****************************************************************************************
!>
!  the purpose of lmder is to minimize the sum of the squares of
!  m nonlinear functions in n variables by a modification of
!  the levenberg-marquardt algorithm. the user must provide a
!  subroutine which calculates the functions and the jacobian.

    subroutine lmder(fcn, m, n, x, Fvec, Fjac, Ldfjac, Ftol, Xtol, Gtol, Maxfev, &
                     Diag, Mode, Factor, Nprint, Info, Nfev, Njev, Ipvt, Qtf, &
                     Wa1, Wa2, Wa3, Wa4, &
                     IXP, IYP, IZP, IQUAT, &
                     JXP, JYP, JZP, JQUAT)

        implicit none

        procedure(fcn_lmder) :: fcn !! the user-supplied subroutine which
                                !! calculates the functions and the jacobian
        integer, intent(in) :: m !! a positive integer input variable set to the number
                            !! of functions.
        integer, intent(in) :: n !! a positive integer input variable set to the number
                            !! of variables. n must not exceed m.
        integer, intent(in) :: Ldfjac !! a positive integer input variable not less than m
                                 !! which specifies the leading dimension of the array fjac.
        integer, intent(in) :: Maxfev !! a positive integer input variable. termination
                                 !! occurs when the number of calls to fcn with iflag = 1
                                 !! has reached maxfev.
        integer, intent(in) :: Mode !! an integer input variable. if mode = 1, the
                                !! variables will be scaled internally. if mode = 2,
                                !! the scaling is specified by the input diag. other
                                !! values of mode are equivalent to mode = 1.
        integer, intent(in) :: Nprint !! an integer input variable that enables controlled
                                 !! printing of iterates if it is positive. in this case,
                                 !! fcn is called with iflag = 0 at the beginning of the first
                                 !! iteration and every nprint iterations thereafter and
                                 !! immediately prior to return, with x, fvec, and fjac
                                 !! available for printing. fvec and fjac should not be
                                 !! altered. if nprint is not positive, no special calls
                                 !! of fcn with iflag = 0 are made.
        integer, intent(out) :: Info !! an integer output variable. if the user has
                                !! terminated execution, info is set to the (negative)
                                !! value of iflag. see description of fcn. otherwise,
                                !! info is set as follows:
                                !!
                                !!  * ***info = 0***  improper input parameters.
                                !!  * ***info = 1***  both actual and predicted relative reductions
                                !!    in the sum of squares are at most ftol.
                                !!  * ***info = 2***  relative error between two consecutive iterates
                                !!    is at most xtol.
                                !!  * ***info = 3***  conditions for info = 1 and info = 2 both hold.
                                !!  * ***info = 4***  the cosine of the angle between fvec and any
                                !!    column of the jacobian is at most gtol in
                                !!    absolute value.
                                !!  * ***info = 5***  number of calls to fcn with iflag = 1 has
                                !!    reached maxfev.
                                !!  * ***info = 6***  ftol is too small. no further reduction in
                                !!    the sum of squares is possible.
                                !!  * ***info = 7***  xtol is too small. no further improvement in
                                !!    the approximate solution x is possible.
                                !!  * ***info = 8***  gtol is too small. fvec is orthogonal to the
                                !!    columns of the jacobian to machine precision.
        integer, intent(out) :: Nfev !! an integer output variable set to the number of
                                !! calls to fcn with iflag = 1.
        integer, intent(out) :: Njev !! an integer output variable set to the number of
                                !! calls to fcn with iflag = 2.
        integer, intent(out) :: Ipvt(n) !! an integer output array of length n. ipvt
                                   !! defines a permutation matrix p such that jac*p = q*r,
                                   !! where jac is the final calculated jacobian, q is
                                   !! orthogonal (not stored), and r is upper triangular
                                   !! with diagonal elements of nonincreasing magnitude.
                                   !! column j of p is column ipvt(j) of the identity matrix.
        real(wp), intent(in) :: Ftol !! a nonnegative input variable. termination
                                !! occurs when both the actual and predicted relative
                                !! reductions in the sum of squares are at most ftol.
                                !! therefore, ftol measures the relative error desired
                                !! in the sum of squares.
        real(wp), intent(in) :: Xtol !! a nonnegative input variable. termination
                                !! occurs when the relative error between two consecutive
                                !! iterates is at most xtol. therefore, xtol measures the
                                !! relative error desired in the approximate solution.
        real(wp), intent(in) :: Gtol !! a nonnegative input variable. termination
                                !! occurs when the cosine of the angle between fvec and
                                !! any column of the jacobian is at most gtol in absolute
                                !! value. therefore, gtol measures the orthogonality
                                !! desired between the function vector and the columns
                                !! of the jacobian.
        real(wp), intent(in) :: Factor !! a positive input variable used in determining the
                                  !! initial step bound. this bound is set to the product of
                                  !! factor and the euclidean norm of diag*x if nonzero, or else
                                  !! to factor itself. in most cases factor should lie in the
                                  !! interval (.1,100.).100. is a generally recommended value.
        real(wp), intent(inout) :: x(n) !! an array of length n. on input x must contain
                                   !! an initial estimate of the solution vector. on output x
                                   !! contains the final estimate of the solution vector.
        real(wp), intent(out) :: Fvec(m) !! an output array of length m which contains
                                    !! the functions evaluated at the output x.
        real(wp), intent(out) :: Fjac(Ldfjac, n) !! an output m by n array. the upper n by n submatrix
                                            !! of fjac contains an upper triangular matrix r with
                                            !! diagonal elements of nonincreasing magnitude such that
                                            !!```
                                            !!        t     t           t
                                            !!       p *(jac *jac)*p = r *r,
                                            !!```
                                            !! where p is a permutation matrix and jac is the final
                                            !! calculated jacobian. column j of p is column ipvt(j)
                                            !! (see below) of the identity matrix. the lower trapezoidal
                                            !! part of fjac contains information generated during
                                            !! the computation of r.
        real(wp), intent(inout) :: Diag(n) !! an array of length n. if mode = 1 (see
                                      !! below), diag is internally set. if mode = 2, diag
                                      !! must contain positive entries that serve as
                                      !! multiplicative scale factors for the variables.
        real(wp), intent(out) :: Qtf(n) !! an output array of length n which contains
                                   !! the first n elements of the vector (q transpose)*fvec.
        real(wp), intent(inout) :: Wa1(n) !! work array of length n.
        real(wp), intent(inout) :: Wa2(n) !! work array of length n.
        real(wp), intent(inout) :: Wa3(n) !! work array of length n.
        real(wp), intent(inout) :: Wa4(m) !! work array of length m.


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Particle 1
        real(kind=8), intent(in) :: IXP, IYP, IZP
        type(quaternion), intent(in) :: IQUAT
        ! Particle 2
        real(kind=8), intent(in) :: JXP, JYP, JZP
        type(quaternion), intent(in) :: JQUAT
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        integer :: i, iflag, iter, j, l
        real(wp) :: actred, delta, dirder, fnorm, fnorm1, gnorm, par, &
                    pnorm, prered, ratio, sum, temp, temp1, temp2, xnorm

        real(wp), parameter :: p1 = 1.0e-1_wp
        real(wp), parameter :: p5 = 5.0e-1_wp
        real(wp), parameter :: p25 = 2.5e-1_wp
        real(wp), parameter :: p75 = 7.5e-1_wp
        real(wp), parameter :: p0001 = 1.0e-4_wp

        Info = 0
        iflag = 0
        Nfev = 0
        Njev = 0

        main : block

            ! check the input parameters for errors.

            if (n > 0 .and. m >= n .and. Ldfjac >= m .and. Ftol >= zero .and. &
                Xtol >= zero .and. Gtol >= zero .and. Maxfev > 0 .and. &
                Factor > zero) then
                if (Mode == 2) then
                    do j = 1, n
                        if (Diag(j) <= zero) exit main
                    end do
                end if
            else
                exit main
            end if

            ! evaluate the function at the starting point
            ! and calculate its norm.

            iflag = 1
            call fcn(m, n, x, Fvec, Fjac, Ldfjac, iflag, &
                    IXP, IYP, IZP, IQUAT, &
                    JXP, JYP, JZP, JQUAT)

            Nfev = 1
            if (iflag < 0) exit main
            fnorm = enorm(m, Fvec)

            ! initialize levenberg-marquardt parameter and iteration counter.

            par = zero
            iter = 1

            ! beginning of the outer loop.

            outer : do

                ! calculate the jacobian matrix.

                iflag = 2
                call fcn(m, n, x, Fvec, Fjac, Ldfjac, iflag, &
                        IXP, IYP, IZP, IQUAT, &
                        JXP, JYP, JZP, JQUAT)

                Njev = Njev + 1
                if (iflag < 0) exit main

                ! if requested, call fcn to enable printing of iterates.

                if (Nprint > 0) then
                    iflag = 0
                    if (mod(iter - 1, Nprint) == 0) &
                        call fcn(m, n, x, Fvec, Fjac, Ldfjac, iflag, &
                                 IXP, IYP, IZP, IQUAT, &
                                 JXP, JYP, JZP, JQUAT)

                    if (iflag < 0) exit main
                end if

                ! compute the qr factorization of the jacobian.

                call qrfac(m, n, Fjac, Ldfjac, .true., Ipvt, n, Wa1, Wa2, Wa3)

                ! on the first iteration and if mode is 1, scale according
                ! to the norms of the columns of the initial jacobian.

                if (iter == 1) then
                    if (Mode /= 2) then
                        do j = 1, n
                            Diag(j) = Wa2(j)
                            if (Wa2(j) == zero) Diag(j) = one
                        end do
                    end if

                    ! on the first iteration, calculate the norm of the scaled x
                    ! and initialize the step bound delta.

                    do j = 1, n
                        Wa3(j) = Diag(j)*x(j)
                    end do
                    xnorm = enorm(n, Wa3)
                    delta = Factor*xnorm
                    if (delta == zero) delta = Factor
                end if

                ! form (q transpose)*fvec and store the first n components in
                ! qtf.

                do i = 1, m
                    Wa4(i) = Fvec(i)
                end do
                do j = 1, n
                    if (Fjac(j, j) /= zero) then
                        sum = zero
                        do i = j, m
                            sum = sum + Fjac(i, j)*Wa4(i)
                        end do
                        temp = -sum/Fjac(j, j)
                        do i = j, m
                            Wa4(i) = Wa4(i) + Fjac(i, j)*temp
                        end do
                    end if
                    Fjac(j, j) = Wa1(j)
                    Qtf(j) = Wa4(j)
                end do

                ! compute the norm of the scaled gradient.

                gnorm = zero
                if (fnorm /= zero) then
                    do j = 1, n
                        l = Ipvt(j)
                        if (Wa2(l) /= zero) then
                            sum = zero
                            do i = 1, j
                                sum = sum + Fjac(i, j)*(Qtf(i)/fnorm)
                            end do
                            gnorm = max(gnorm, abs(sum/Wa2(l)))
                        end if
                    end do
                end if

                ! test for convergence of the gradient norm.

                if (gnorm <= Gtol) Info = 4
                if (Info /= 0) exit main

                ! rescale if necessary.

                if (Mode /= 2) then
                    do j = 1, n
                        Diag(j) = max(Diag(j), Wa2(j))
                    end do
                end if

                ! beginning of the inner loop.
                inner : do

                    ! determine the levenberg-marquardt parameter.

                    call lmpar(n, Fjac, Ldfjac, Ipvt, Diag, Qtf, delta, par, Wa1, Wa2, Wa3, Wa4)

                    ! store the direction p and x + p. calculate the norm of p.

                    do j = 1, n
                        Wa1(j) = -Wa1(j)
                        Wa2(j) = x(j) + Wa1(j)
                        Wa3(j) = Diag(j)*Wa1(j)
                    end do
                    pnorm = enorm(n, Wa3)

                    ! on the first iteration, adjust the initial step bound.

                    if (iter == 1) delta = min(delta, pnorm)

                    ! evaluate the function at x + p and calculate its norm.

                    iflag = 1
                    call fcn(m, n, Wa2, Wa4, Fjac, Ldfjac, iflag, &
                            IXP, IYP, IZP, IQUAT, &
                            JXP, JYP, JZP, JQUAT)

                    Nfev = Nfev + 1
                    if (iflag < 0) exit main
                    fnorm1 = enorm(m, Wa4)

                    ! compute the scaled actual reduction.

                    actred = -one
                    if (p1*fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2

                    ! compute the scaled predicted reduction and
                    ! the scaled directional derivative.

                    do j = 1, n
                        Wa3(j) = zero
                        l = Ipvt(j)
                        temp = Wa1(l)
                        do i = 1, j
                            Wa3(i) = Wa3(i) + Fjac(i, j)*temp
                        end do
                    end do
                    temp1 = enorm(n, Wa3)/fnorm
                    temp2 = (sqrt(par)*pnorm)/fnorm
                    prered = temp1**2 + temp2**2/p5
                    dirder = -(temp1**2 + temp2**2)

                    ! compute the ratio of the actual to the predicted
                    ! reduction.

                    ratio = zero
                    if (prered /= zero) ratio = actred/prered

                    ! update the step bound.

                    if (ratio <= p25) then
                        if (actred >= zero) temp = p5
                        if (actred < zero) temp = p5*dirder/(dirder + p5*actred)
                        if (p1*fnorm1 >= fnorm .or. temp < p1) temp = p1
                        delta = temp*min(delta, pnorm/p1)
                        par = par/temp
                    elseif (par == zero .or. ratio >= p75) then
                        delta = pnorm/p5
                        par = p5*par
                    end if

                    ! test for successful iteration.

                    if (ratio >= p0001) then
                        ! successful iteration. update x, fvec, and their norms.
                        do j = 1, n
                            x(j) = Wa2(j)
                            Wa2(j) = Diag(j)*x(j)
                        end do
                        do i = 1, m
                            Fvec(i) = Wa4(i)
                        end do
                        xnorm = enorm(n, Wa2)
                        fnorm = fnorm1
                        iter = iter + 1
                    end if

                    ! tests for convergence.
                    if (abs(actred) <= Ftol .and. prered <= Ftol .and. p5*ratio <= one) Info = 1
                    if (delta <= Xtol*xnorm) Info = 2
                    if (abs(actred) <= Ftol .and. prered <= Ftol .and. p5*ratio <= one .and. Info == 2) Info = 3
                    if (Info /= 0) exit main

                    ! tests for termination and stringent tolerances.
                    if (Nfev >= Maxfev) Info = 5
                    if (abs(actred) <= epsmch .and. prered <= epsmch .and. p5*ratio <= one) Info = 6
                    if (delta <= epsmch*xnorm) Info = 7
                    if (gnorm <= epsmch) Info = 8
                    if (Info /= 0) exit main

                    if (ratio >= p0001) exit inner

                end do inner ! end of the inner loop. repeat if iteration unsuccessful.

            end do outer ! end of the outer loop

        end block main

        ! termination, either normal or user imposed.

        if (iflag < 0) Info = iflag
        iflag = 0
        if (Nprint > 0) call fcn(m, n, x, Fvec, Fjac, Ldfjac, iflag, &
                                IXP, IYP, IZP, IQUAT, &
                                JXP, JYP, JZP, JQUAT)


    end subroutine lmder
!*****************************************************************************************


!*****************************************************************************************
!>
!  the purpose of lmder1 is to minimize the sum of the squares of
!  m nonlinear functions in n variables by a modification of the
!  levenberg-marquardt algorithm. this is done by using the more
!  general least-squares solver lmder. the user must provide a
!  subroutine which calculates the functions and the jacobian.

    subroutine lmder1(fcn, m, n, x, Fvec, Fjac, Ldfjac, Tol, Info, Ipvt, Wa, Lwa, &
                      IXP, IYP, IZP, IQUAT, &
                      JXP, JYP, JZP, JQUAT)

        implicit none

        procedure(fcn_lmder) :: fcn !! user-supplied subroutine which
                                    !! calculates the functions and the jacobian.
        integer, intent(in) :: m !! a positive integer input variable set to the number
                                !! of functions.
        integer, intent(in) :: n !! a positive integer input variable set to the number
                                !! of variables. n must not exceed m.
        integer, intent(in) :: Ldfjac !! a positive integer input variable not less than m
                                     !! which specifies the leading dimension of the array fjac.
        integer, intent(out) :: Info !! an integer output variable. if the user has
                                    !! terminated execution, info is set to the (negative)
                                    !! value of iflag. see description of fcn. otherwise,
                                    !! info is set as follows.
                                    !!
                                    !!  * ***info = 0***  improper input parameters.
                                    !!  * ***info = 1***  algorithm estimates that the relative error
                                    !!    in the sum of squares is at most tol.
                                    !!  * ***info = 2***  algorithm estimates that the relative error
                                    !!    between x and the solution is at most tol.
                                    !!  * ***info = 3***  conditions for info = 1 and info = 2 both hold.
                                    !!  * ***info = 4***  fvec is orthogonal to the columns of the
                                    !!    jacobian to machine precision.
                                    !!  * ***info = 5***  number of calls to fcn with iflag = 1 has
                                    !!    reached 100*(n+1).
                                    !!  * ***info = 6***  tol is too small. no further reduction in
                                    !!    the sum of squares is possible.
                                    !!  * ***info = 7***  tol is too small. no further improvement in
                                    !!    the approximate solution x is possible.
        integer, intent(in) :: Lwa !! a positive integer input variable not less than 5*n+m.
        integer, intent(out) :: Ipvt(n) !! an integer output array of length n. ipvt
                                       !! defines a permutation matrix p such that jac*p = q*r,
                                       !! where jac is the final calculated jacobian, q is
                                       !! orthogonal (not stored), and r is upper triangular
                                       !! with diagonal elements of nonincreasing magnitude.
                                       !! column j of p is column ipvt(j) of the identity matrix.
        real(wp), intent(in) :: Tol !! a nonnegative input variable. termination occurs
                                   !! when the algorithm estimates either that the relative
                                   !! error in the sum of squares is at most tol or that
                                   !! the relative error between x and the solution is at
                                   !! most tol.
        real(wp), intent(inout) :: x(n) !! an array of length n. on input x must contain
                                       !! an initial estimate of the solution vector. on output x
                                       !! contains the final estimate of the solution vector.
        real(wp), intent(out) :: Fvec(m) !! an output array of length m which contains
                                        !! the functions evaluated at the output x.
        real(wp), intent(out) :: Fjac(Ldfjac, n) !! an output m by n array. the upper n by n submatrix
                                                !! of fjac contains an upper triangular matrix r with
                                                !! diagonal elements of nonincreasing magnitude such that
                                                !!```
                                                !!        t     t           t
                                                !!       p *(jac *jac)*p = r *r,
                                                !!```
                                                !! where p is a permutation matrix and jac is the final
                                                !! calculated jacobian. column j of p is column ipvt(j)
                                                !! (see below) of the identity matrix. the lower trapezoidal
                                                !! part of fjac contains information generated during
                                                !! the computation of r.
        real(wp), intent(inout) :: Wa(Lwa) !! a work array of length lwa.

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Particle 1
        real(kind=8), intent(in) :: IXP, IYP, IZP
        type(quaternion), intent(in) :: IQUAT
        ! Particle 2
        real(kind=8), intent(in) :: JXP, JYP, JZP
        type(quaternion), intent(in) :: JQUAT
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        integer :: maxfev, mode, nfev, njev, nprint
        real(wp) :: ftol, gtol, xtol

        real(wp), parameter :: factor = 100.0_wp

        Info = 0

        ! check the input parameters for errors.

        if (n > 0 .and. m >= n .and. Ldfjac >= m .and. Tol >= zero .and. &
            Lwa >= 5*n + m) then
            ! call lmder.
            maxfev = 100*(n + 1)
            ftol = Tol
            xtol = Tol
            gtol = zero
            mode = 1
            nprint = 0
            call lmder(fcn, m, n, x, Fvec, Fjac, Ldfjac, ftol, xtol, gtol, maxfev,   &
                     & Wa(1), mode, factor, nprint, Info, nfev, njev, Ipvt, Wa(n + 1)&
                     & , Wa(2*n + 1), Wa(3*n + 1), Wa(4*n + 1), Wa(5*n + 1), &
                     IXP, IYP, IZP, IQUAT, &
                     JXP, JYP, JZP, JQUAT)

            if (Info == 8) Info = 4
        end if

    end subroutine lmder1
!*****************************************************************************************

!*****************************************************************************************
!>
!  given an m by n matrix a, an n by n nonsingular diagonal
!  matrix d, an m-vector b, and a positive number delta,
!  the problem is to determine a value for the parameter
!  par such that if x solves the system
!```
!        a*x = b ,     sqrt(par)*d*x = 0 ,
!```
!  in the least squares sense, and dxnorm is the euclidean
!  norm of d*x, then either par is zero and
!```
!        (dxnorm-delta) <= 0.1*delta ,
!```
!  or par is positive and
!```
!        abs(dxnorm-delta) <= 0.1*delta .
!```
!  this subroutine completes the solution of the problem
!  if it is provided with the necessary information from the
!  qr factorization, with column pivoting, of a. that is, if
!  a*p = q*r, where p is a permutation matrix, q has orthogonal
!  columns, and r is an upper triangular matrix with diagonal
!  elements of nonincreasing magnitude, then lmpar expects
!  the full upper triangle of r, the permutation matrix p,
!  and the first n components of (q transpose)*b. on output
!  lmpar also provides an upper triangular matrix s such that
!```
!         t   t                   t
!        p *(a *a + par*d*d)*p = s *s .
!```
!  s is employed within lmpar and may be of separate interest.
!
!  only a few iterations are generally needed for convergence
!  of the algorithm. if, however, the limit of 10 iterations
!  is reached, then the output par will contain the best
!  value obtained so far.

    subroutine lmpar(n, r, Ldr, Ipvt, Diag, Qtb, Delta, Par, x, Sdiag, Wa1, Wa2)
        implicit none

        integer, intent(in) :: n !! a positive integer input variable set to the order of r.
        integer, intent(in) :: Ldr !! a positive integer input variable not less than n
                                  !! which specifies the leading dimension of the array r.
        integer, intent(in) :: Ipvt(n) !! an integer input array of length n which defines the
                                      !! permutation matrix p such that a*p = q*r. column j of p
                                      !! is column ipvt(j) of the identity matrix.
        real(wp) :: Delta !! a positive input variable which specifies an upper
                          !! bound on the euclidean norm of d*x.
        real(wp), intent(inout) :: Par !! a nonnegative variable. on input par contains an
                                      !! initial estimate of the levenberg-marquardt parameter.
                                      !! on output par contains the final estimate.
        real(wp), intent(inout) :: r(Ldr, n) !! an n by n array. on input the full upper triangle
                                            !! must contain the full upper triangle of the matrix r.
                                            !! on output the full upper triangle is unaltered, and the
                                            !! strict lower triangle contains the strict upper triangle
                                            !! (transposed) of the upper triangular matrix s.
        real(wp), intent(in) :: Diag(n) !! an input array of length n which must contain the
                                       !! diagonal elements of the matrix d.
        real(wp), intent(in) :: Qtb(n) !! an input array of length n which must contain the first
                                      !! n elements of the vector (q transpose)*b.
        real(wp), intent(out) :: x(n) !! an output array of length n which contains the least
                                     !! squares solution of the system a*x = b, sqrt(par)*d*x = 0,
                                     !! for the output par.
        real(wp), intent(out) :: Sdiag(n) !! an output array of length n which contains the
                                         !! diagonal elements of the upper triangular matrix s.
        real(wp), intent(inout) :: Wa1(n) !! work array of length n.
        real(wp), intent(inout) :: Wa2(n) !! work array of length n.

        integer :: i, iter, j, jm1, jp1, k, l, nsing
        real(wp) :: dxnorm, fp, gnorm, parc, parl, paru, sum, temp

        real(wp), parameter :: p1 = 1.0e-1_wp
        real(wp), parameter :: p001 = 1.0e-3_wp
        real(wp), parameter :: dwarf = dpmpar(2) !! the smallest positive magnitude

        ! compute and store in x the gauss-newton direction. if the
        ! jacobian is rank-deficient, obtain a least squares solution.

        nsing = n
        do j = 1, n
            Wa1(j) = Qtb(j)
            if (r(j, j) == zero .and. nsing == n) nsing = j - 1
            if (nsing < n) Wa1(j) = zero
        end do
        if (nsing >= 1) then
            do k = 1, nsing
                j = nsing - k + 1
                Wa1(j) = Wa1(j)/r(j, j)
                temp = Wa1(j)
                jm1 = j - 1
                if (jm1 >= 1) then
                    do i = 1, jm1
                        Wa1(i) = Wa1(i) - r(i, j)*temp
                    end do
                end if
            end do
        end if
        do j = 1, n
            l = Ipvt(j)
            x(l) = Wa1(j)
        end do

        ! initialize the iteration counter.
        ! evaluate the function at the origin, and test
        ! for acceptance of the gauss-newton direction.

        iter = 0
        do j = 1, n
            Wa2(j) = Diag(j)*x(j)
        end do
        dxnorm = enorm(n, Wa2)
        fp = dxnorm - Delta
        if (fp <= p1*Delta) then
            ! termination.
            if (iter == 0) Par = zero
        else

            ! if the jacobian is not rank deficient, the newton
            ! step provides a lower bound, parl, for the zero of
            ! the function. otherwise set this bound to zero.

            parl = zero
            if (nsing >= n) then
                do j = 1, n
                    l = Ipvt(j)
                    Wa1(j) = Diag(l)*(Wa2(l)/dxnorm)
                end do
                do j = 1, n
                    sum = zero
                    jm1 = j - 1
                    if (jm1 >= 1) then
                        do i = 1, jm1
                            sum = sum + r(i, j)*Wa1(i)
                        end do
                    end if
                    Wa1(j) = (Wa1(j) - sum)/r(j, j)
                end do
                temp = enorm(n, Wa1)
                parl = ((fp/Delta)/temp)/temp
            end if

            ! calculate an upper bound, paru, for the zero of the function.

            do j = 1, n
                sum = zero
                do i = 1, j
                    sum = sum + r(i, j)*Qtb(i)
                end do
                l = Ipvt(j)
                Wa1(j) = sum/Diag(l)
            end do
            gnorm = enorm(n, Wa1)
            paru = gnorm/Delta
            if (paru == zero) paru = dwarf/min(Delta, p1)

            ! if the input par lies outside of the interval (parl,paru),
            ! set par to the closer endpoint.

            Par = max(Par, parl)
            Par = min(Par, paru)
            if (Par == zero) Par = gnorm/dxnorm

            ! beginning of an iteration.
            do

                iter = iter + 1

                ! evaluate the function at the current value of par.

                if (Par == zero) Par = max(dwarf, p001*paru)
                temp = sqrt(Par)
                do j = 1, n
                    Wa1(j) = temp*Diag(j)
                end do
                call qrsolv(n, r, Ldr, Ipvt, Wa1, Qtb, x, Sdiag, Wa2)
                do j = 1, n
                    Wa2(j) = Diag(j)*x(j)
                end do
                dxnorm = enorm(n, Wa2)
                temp = fp
                fp = dxnorm - Delta

                ! if the function is small enough, accept the current value
                ! of par. also test for the exceptional cases where parl
                ! is zero or the number of iterations has reached 10.

                if (abs(fp) <= p1*Delta .or. parl == zero .and. fp <= temp .and. &
                    temp < zero .or. iter == 10) then
                    if (iter == 0) Par = zero
                    exit
                else

                    ! compute the newton correction.

                    do j = 1, n
                        l = Ipvt(j)
                        Wa1(j) = Diag(l)*(Wa2(l)/dxnorm)
                    end do
                    do j = 1, n
                        Wa1(j) = Wa1(j)/Sdiag(j)
                        temp = Wa1(j)
                        jp1 = j + 1
                        if (n >= jp1) then
                            do i = jp1, n
                                Wa1(i) = Wa1(i) - r(i, j)*temp
                            end do
                        end if
                    end do
                    temp = enorm(n, Wa1)
                    parc = ((fp/Delta)/temp)/temp

                    ! depending on the sign of the function, update parl or paru.

                    if (fp > zero) parl = max(parl, Par)
                    if (fp < zero) paru = min(paru, Par)

                    ! compute an improved estimate for par.

                    Par = max(parl, Par + parc)

                end if

            end do ! end of an iteration.

        end if

    end subroutine lmpar
!*****************************************************************************************

!*****************************************************************************************
!>
!  this subroutine uses householder transformations with column
!  pivoting (optional) to compute a qr factorization of the
!  m by n matrix a. that is, qrfac determines an orthogonal
!  matrix q, a permutation matrix p, and an upper trapezoidal
!  matrix r with diagonal elements of nonincreasing magnitude,
!  such that a*p = q*r. the householder transformation for
!  column k, k = 1,2,...,min(m,n), is of the form
!```
!                        t
!        i - (1/u(k))*u*u
!```
!  where u has zeros in the first k-1 positions. the form of
!  this transformation and the method of pivoting first
!  appeared in the corresponding linpack subroutine.

    subroutine qrfac(m, n, a, Lda, Pivot, Ipvt, Lipvt, Rdiag, Acnorm, Wa)
        implicit none

        integer, intent(in) :: m !! a positive integer input variable set to the number
                                !! of rows of a.
        integer, intent(in) :: n !! a positive integer input variable set to the number
                                !! of columns of a.
        integer, intent(in) :: Lda !! a positive integer input variable not less than m
                                  !! which specifies the leading dimension of the array a.
        integer, intent(in) :: Lipvt !! a positive integer input variable. if pivot is false,
                                    !! then lipvt may be as small as 1. if pivot is true, then
                                    !! lipvt must be at least n.
        integer, intent(out) :: Ipvt(Lipvt) !! an integer output array of length lipvt. ipvt
                                           !! defines the permutation matrix p such that a*p = q*r.
                                           !! column j of p is column ipvt(j) of the identity matrix.
                                           !! if pivot is false, ipvt is not referenced.
        logical, intent(in) :: Pivot !! a logical input variable. if pivot is set true,
                                    !! then column pivoting is enforced. if pivot is set false,
                                    !! then no column pivoting is done.
        real(wp), intent(inout) :: a(Lda, n) !! an m by n array. on input a contains the matrix for
                                            !! which the qr factorization is to be computed. on output
                                            !! the strict upper trapezoidal part of a contains the strict
                                            !! upper trapezoidal part of r, and the lower trapezoidal
                                            !! part of a contains a factored form of q (the non-trivial
                                            !! elements of the u vectors described above).
        real(wp), intent(out) :: Rdiag(n) !! an output array of length n which contains the
                                         !! diagonal elements of r.
        real(wp), intent(out) :: Acnorm(n) !! an output array of length n which contains the
                                          !! norms of the corresponding columns of the input matrix a.
                                          !! if this information is not needed, then acnorm can coincide
                                          !! with rdiag.
        real(wp), intent(inout) :: Wa(n) !! a work array of length n. if pivot is false, then wa
                                        !! can coincide with rdiag.

        integer :: i, j, jp1, k, kmax, minmn
        real(wp) :: ajnorm, sum, temp

        real(wp), parameter :: p05 = 5.0e-2_wp

        ! compute the initial column norms and initialize several arrays.

        do j = 1, n
            Acnorm(j) = enorm(m, a(1, j))
            Rdiag(j) = Acnorm(j)
            Wa(j) = Rdiag(j)
            if (Pivot) Ipvt(j) = j
        end do

        ! reduce a to r with householder transformations.

        minmn = min(m, n)
        do j = 1, minmn
            if (Pivot) then

                ! bring the column of largest norm into the pivot position.

                kmax = j
                do k = j, n
                    if (Rdiag(k) > Rdiag(kmax)) kmax = k
                end do
                if (kmax /= j) then
                    do i = 1, m
                        temp = a(i, j)
                        a(i, j) = a(i, kmax)
                        a(i, kmax) = temp
                    end do
                    Rdiag(kmax) = Rdiag(j)
                    Wa(kmax) = Wa(j)
                    k = Ipvt(j)
                    Ipvt(j) = Ipvt(kmax)
                    Ipvt(kmax) = k
                end if
            end if

            ! compute the householder transformation to reduce the
            ! j-th column of a to a multiple of the j-th unit vector.

            ajnorm = enorm(m - j + 1, a(j, j))
            if (ajnorm /= zero) then
                if (a(j, j) < zero) ajnorm = -ajnorm
                do i = j, m
                    a(i, j) = a(i, j)/ajnorm
                end do
                a(j, j) = a(j, j) + one

                ! apply the transformation to the remaining columns
                ! and update the norms.

                jp1 = j + 1
                if (n >= jp1) then
                    do k = jp1, n
                        sum = zero
                        do i = j, m
                            sum = sum + a(i, j)*a(i, k)
                        end do
                        temp = sum/a(j, j)
                        do i = j, m
                            a(i, k) = a(i, k) - temp*a(i, j)
                        end do
                        if (.not. (.not. Pivot .or. Rdiag(k) == zero)) then
                            temp = a(j, k)/Rdiag(k)
                            Rdiag(k) = Rdiag(k)*sqrt(max(zero, one - temp**2))
                            if (p05*(Rdiag(k)/Wa(k))**2 <= epsmch) then
                                Rdiag(k) = enorm(m - j, a(jp1, k))
                                Wa(k) = Rdiag(k)
                            end if
                        end if
                    end do
                end if
            end if
            Rdiag(j) = -ajnorm
        end do

    end subroutine qrfac
!*****************************************************************************************

!*****************************************************************************************
!>
!  given an m by n matrix a, an n by n diagonal matrix d,
!  and an m-vector b, the problem is to determine an x which
!  solves the system
!```
!        a*x = b ,     d*x = 0 ,
!```
!  in the least squares sense.
!
!  this subroutine completes the solution of the problem
!  if it is provided with the necessary information from the
!  qr factorization, with column pivoting, of a. that is, if
!  a*p = q*r, where p is a permutation matrix, q has orthogonal
!  columns, and r is an upper triangular matrix with diagonal
!  elements of nonincreasing magnitude, then qrsolv expects
!  the full upper triangle of r, the permutation matrix p,
!  and the first n components of (q transpose)*b. the system
!  a*x = b, d*x = 0, is then equivalent to
!```
!               t       t
!        r*z = q *b ,  p *d*p*z = 0 ,
!```
!  where x = p*z. if this system does not have full rank,
!  then a least squares solution is obtained. on output qrsolv
!  also provides an upper triangular matrix s such that
!```
!         t   t               t
!        p *(a *a + d*d)*p = s *s .
!```
!  s is computed within qrsolv and may be of separate interest.

    subroutine qrsolv(n, r, Ldr, Ipvt, Diag, Qtb, x, Sdiag, Wa)
        implicit none

        integer, intent(in) :: n !! a positive integer input variable set to the order of r.
        integer, intent(in) :: Ldr !! a positive integer input variable not less than n
                                  !! which specifies the leading dimension of the array r.
        integer, intent(in) :: Ipvt(n) !! an integer input array of length n which defines the
                                      !! permutation matrix p such that a*p = q*r. column j of p
                                      !! is column ipvt(j) of the identity matrix.
        real(wp), intent(inout) :: r(Ldr, n) !! an n by n array. on input the full upper triangle
                                            !! must contain the full upper triangle of the matrix r.
                                            !! on output the full upper triangle is unaltered, and the
                                            !! strict lower triangle contains the strict upper triangle
                                            !! (transposed) of the upper triangular matrix s.
        real(wp), intent(in) :: Diag(n) !! an input array of length n which must contain the
                                       !! diagonal elements of the matrix d.
        real(wp), intent(in) :: Qtb(n) !! an input array of length n which must contain the first
                                      !! n elements of the vector (q transpose)*b.
        real(wp), intent(out) :: x(n) !! an output array of length n which contains the least
                                     !! squares solution of the system a*x = b, d*x = 0.
        real(wp), intent(out) :: Sdiag(n) !! an output array of length n which contains the
                                         !! diagonal elements of the upper triangular matrix s.
        real(wp), intent(inout) :: Wa(n) !! a work array of length n.

        integer :: i, j, jp1, k, kp1, l, nsing
        real(wp) :: cos, cotan, qtbpj, sin, sum, tan, temp

        real(wp), parameter :: p5 = 5.0e-1_wp
        real(wp), parameter :: p25 = 2.5e-1_wp

        ! copy r and (q transpose)*b to preserve input and initialize s.
        ! in particular, save the diagonal elements of r in x.

        do j = 1, n
            do i = j, n
                r(i, j) = r(j, i)
            end do
            x(j) = r(j, j)
            Wa(j) = Qtb(j)
        end do

        ! eliminate the diagonal matrix d using a givens rotation.

        do j = 1, n

            ! prepare the row of d to be eliminated, locating the
            ! diagonal element using p from the qr factorization.

            l = Ipvt(j)
            if (Diag(l) /= zero) then
                do k = j, n
                    Sdiag(k) = zero
                end do
                Sdiag(j) = Diag(l)

                ! the transformations to eliminate the row of d
                ! modify only a single element of (q transpose)*b
                ! beyond the first n, which is initially zero.

                qtbpj = zero
                do k = j, n

                    ! determine a givens rotation which eliminates the
                    ! appropriate element in the current row of d.

                    if (Sdiag(k) /= zero) then
                        if (abs(r(k, k)) >= abs(Sdiag(k))) then
                            tan = Sdiag(k)/r(k, k)
                            cos = p5/sqrt(p25 + p25*tan**2)
                            sin = cos*tan
                        else
                            cotan = r(k, k)/Sdiag(k)
                            sin = p5/sqrt(p25 + p25*cotan**2)
                            cos = sin*cotan
                        end if

                        ! compute the modified diagonal element of r and
                        ! the modified element of ((q transpose)*b,0).

                        r(k, k) = cos*r(k, k) + sin*Sdiag(k)
                        temp = cos*Wa(k) + sin*qtbpj
                        qtbpj = -sin*Wa(k) + cos*qtbpj
                        Wa(k) = temp

                        ! accumulate the tranformation in the row of s.

                        kp1 = k + 1
                        if (n >= kp1) then
                            do i = kp1, n
                                temp = cos*r(i, k) + sin*Sdiag(i)
                                Sdiag(i) = -sin*r(i, k) + cos*Sdiag(i)
                                r(i, k) = temp
                            end do
                        end if
                    end if
                end do
            end if

            ! store the diagonal element of s and restore
            ! the corresponding diagonal element of r.

            Sdiag(j) = r(j, j)
            r(j, j) = x(j)
        end do

        ! solve the triangular system for z. if the system is
        ! singular, then obtain a least squares solution.

        nsing = n
        do j = 1, n
            if (Sdiag(j) == zero .and. nsing == n) nsing = j - 1
            if (nsing < n) Wa(j) = zero
        end do
        if (nsing >= 1) then
            do k = 1, nsing
                j = nsing - k + 1
                sum = zero
                jp1 = j + 1
                if (nsing >= jp1) then
                    do i = jp1, nsing
                        sum = sum + r(i, j)*Wa(i)
                    end do
                end if
                Wa(j) = (Wa(j) - sum)/Sdiag(j)
            end do
        end if

        ! permute the components of z back to components of x.

        do j = 1, n
            l = Ipvt(j)
            x(l) = Wa(j)
        end do

    end subroutine qrsolv
!*****************************************************************************************

!*****************************************************************************************
end module minpack_module
!*****************************************************************************************