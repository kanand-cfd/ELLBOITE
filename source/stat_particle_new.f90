subroutine STAT_PARTICLE(NCYCLE, TIME)

use param_phys
use mod_quaternion
use ellipsoid_particle

implicit none

!=================================!
integer, intent(in) :: NCYCLE
real(kind=8), intent(in) :: TIME
!=================================!
integer :: I

real(kind=8), dimension(100) :: MEAN_PART
!=================================!

real(kind=8) :: UPM, VPM, WPM

real(kind=8) :: UPFLC, VPFLC, WPFLC
!=================================!

real(kind=8) :: OXM, OYM, OZM

real(kind=8) :: OXFLC, OYFLC, OZFLC 

real(kind=8), dimension(ndim, 1) :: global_angular_velocity, omega_fluct, Amom

real(kind=8), dimension(ndim, ndim) :: I_p
!=================================!

MEAN_PART(:) = ZEROP

! Initialization of Inertia Matrix
I_p(:,:) = 0.0
I_p(1,1) = IPXX
I_p(2,2) = IPYY
I_p(3,3) = IPZZ

!!=============================================!!
!!			Statistics for Translation 		   !!
!!=============================================!!
do I = 1, NPART_MAX
!!----------- Particle Velocities -----------!!
!! <Up>
    MEAN_PART(1) = MEAN_PART(1) + ELLP(I)%UP

!! <Vp>
    MEAN_PART(2) = MEAN_PART(2) + ELLP(I)%VP

!! <Wp>
    MEAN_PART(3) = MEAN_PART(3) + ELLP(I)%WP
end do

UPM = MEAN_PART(1)/NPART_MAX
VPM = MEAN_PART(2)/NPART_MAX
WPM = MEAN_PART(3)/NPART_MAX


do I = 1, NPART_MAX

!!- u'p = up-<up>
    UPFLC = ELLP(I)%UP - UPM

!!- v'p = vp-<vp>
    VPFLC = ELLP(I)%VP - VPM

!!- w'p = wp-<wp>
    WPFLC = ELLP(I)%WP - WPM

!!- <upup>
    MEAN_PART(4) = MEAN_PART(4) + UPFLC*UPFLC

!!- <vpvp>
    MEAN_PART(5) = MEAN_PART(5) + VPFLC*VPFLC

!!- <wpwp>
    MEAN_PART(6) = MEAN_PART(6) + WPFLC*WPFLC

!!- <up*vp>
    MEAN_PART(7) = MEAN_PART(7) + UPFLC*VPFLC

!!- <up*wp>
    MEAN_PART(8) = MEAN_PART(8) + UPFLC*WPFLC

!!- <vp*wp>
    MEAN_PART(9) = MEAN_PART(9) + VPFLC*WPFLC    
 
end do

!! - q_p - Translation Kinetic Energy
MEAN_PART(10) = 0.5*(MEAN_PART(4) &
                   + MEAN_PART(5) &
                   + MEAN_PART(6))

!!=============================================!!
!!			Statistics for Rotation 		   !!
!!=============================================!!
do I = 1, NPART_MAX


    MEAN_PART(11) = MEAN_PART(11) + ELLP(I)%OMEGAX

    MEAN_PART(12) = MEAN_PART(12) + ELLP(I)%OMEGAY

    MEAN_PART(13) = MEAN_PART(13) + ELLP(I)%OMEGAZ

end do

OXM = MEAN_PART(11)/NPART_MAX
OYM = MEAN_PART(12)/NPART_MAX
OZM = MEAN_PART(13)/NPART_MAX

do I = 1, NPART_MAX

    OXFLC = ELLP(I)%OMEGAX - OXM

    OYFLC = ELLP(I)%OMEGAY - OYM

    OZFLC = ELLP(I)%OMEGAZ - OZM

    MEAN_PART(14) = MEAN_PART(14) + IPXX*OXFLC*OXFLC

    MEAN_PART(15) = MEAN_PART(15) + IPYY*OYFLC*OYFLC

    MEAN_PART(16) = MEAN_PART(16) + IPZZ*OZFLC*OZFLC

    !MEAN_PART(17) = MEAN_PART(17) + OXFLC*OYFLC

    !MEAN_PART(18) = MEAN_PART(18) + OXFLC*OZFLC

    !MEAN_PART(19) = MEAN_PART(19) + OYFLC*OZFLC


end do

MEAN_PART(21) = 0.5*(IPXX*MEAN_PART(14) + &
                     IPYY*MEAN_PART(15) + &
                     IPZZ*MEAN_PART(16))


MEAN_PART(:) = MEAN_PART(:)/NPART_MAX



if(mod(NCYCLE, STAT_WRITE)==0)then

    write(300, 10000) &
    TIME, &
    MEAN_PART(1), &   !! <Up>
    MEAN_PART(2), &   !! <Vp>
    MEAN_PART(3), &   !! <Wp>
    MEAN_PART(4), &   !! <upup>
    MEAN_PART(5), &   !! <vpvp>
    MEAN_PART(6), &   !! <wpwp>
    MEAN_PART(7), &   !! <upvp>
    MEAN_PART(8), &   !! <upwp>
    MEAN_PART(9), &   !! <wpwp>
    MEAN_PART(10)     !! <qp>

    write(400, 10000) &
    TIME, &
    MEAN_PART(11), &  !! <omegaxp>
    MEAN_PART(12), &  !! <omegayp>
    MEAN_PART(13), &  !! <omegazp>
    MEAN_PART(14), &  !! <oxp.oxp>
    MEAN_PART(15), &  !! <oyp.oyp>
    MEAN_PART(16), &  !! <ozp.ozp>
    MEAN_PART(17), &  !! <oxp.oyp>
    MEAN_PART(18), &  !! <oxp.ozp>
    MEAN_PART(19), &  !! <oyp.ozp>
    MEAN_PART(21)     !! <qp> - Rotational  

end if

!!=============================================!!
!!            Calculation of PDFs              !!
!!=============================================!!




10000 format (30(e17.7))   
end subroutine STAT_PARTICLE