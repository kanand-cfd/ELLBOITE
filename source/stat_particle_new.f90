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

!!--------- Calculate Global Angular Velocity ----------!!
    global_angular_velocity(1,1) = ELLP(I)%OMEGAX
    global_angular_velocity(2,1) = ELLP(I)%OMEGAY
    global_angular_velocity(3,1) = ELLP(I)%OMEGAZ

    call transform_basis(global_angular_velocity, ELLP(I)%ELLQUAT, shape(global_angular_velocity))

    ELLP(I)%GOMEGAX = global_angular_velocity(1,1)
    ELLP(I)%GOMEGAY = global_angular_velocity(2,1)
    ELLP(I)%GOMEGAZ = global_angular_velocity(3,1)

    MEAN_PART(11) = MEAN_PART(11) + ELLP(I)%GOMEGAX

    MEAN_PART(12) = MEAN_PART(12) + ELLP(I)%GOMEGAY

    MEAN_PART(13) = MEAN_PART(13) + ELLP(I)%GOMEGAZ

! global Inertia Matrix
    call transform_basis(I_p, ELLP(I)%ELLQUAT, shape(I_p))
!    write(*,*) I_p(1,:)
!    write(*,*) I_p(2,:)
!    write(*,*) I_p(3,:)
!    stop

end do

OXM = MEAN_PART(11)/NPART_MAX
OYM = MEAN_PART(12)/NPART_MAX
OZM = MEAN_PART(13)/NPART_MAX

do I = 1, NPART_MAX

    OXFLC = ELLP(I)%GOMEGAX !- OXM

    OYFLC = ELLP(I)%GOMEGAY !- OYM

    OZFLC = ELLP(I)%GOMEGAZ !- OZM

    MEAN_PART(14) = MEAN_PART(14) + OXFLC*OXFLC

    MEAN_PART(15) = MEAN_PART(15) + OYFLC*OYFLC

    MEAN_PART(16) = MEAN_PART(16) + OZFLC*OZFLC

    MEAN_PART(17) = MEAN_PART(17) + OXFLC*OYFLC

    MEAN_PART(18) = MEAN_PART(18) + OXFLC*OZFLC

    MEAN_PART(19) = MEAN_PART(19) + OYFLC*OZFLC
 
    omega_fluct(1,1) = OXFLC
    omega_fluct(2,1) = OYFLC
    omega_fluct(3,1) = OZFLC

    Amom = matmul(I_p, omega_fluct)

    ! Rotational Kinetic Energy
    MEAN_PART(20) = MEAN_PART(20) + 0.5*(OXFLC*Amom(1,1) + OYFLC*Amom(2,1) + OZFLC*Amom(3,1))

end do

MEAN_PART(21) = 0.5*(MEAN_PART(14) + &
                     MEAN_PART(15) + &
                     MEAN_PART(16))


MEAN_PART(:) = MEAN_PART(:)/NPART_MAX


write(300, 10000) &
        TIME, MEAN_PART(1), &   !! <Up>
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
        TIME, MEAN_PART(11), &  !! <omegaxp>
              MEAN_PART(12), &  !! <omegayp>
              MEAN_PART(13), &  !! <omegazp>
              MEAN_PART(14), &  !! <oxp.oxp>
              MEAN_PART(15), &  !! <oyp.oyp>
              MEAN_PART(16), &  !! <ozp.ozp>
              MEAN_PART(17), &  !! <oxp.oyp>
              MEAN_PART(18), &  !! <oxp.ozp>
              MEAN_PART(19), &  !! <oyp.ozp>
              MEAN_PART(21)     !! <qp> - Rotational  

!!=============================================!!
!!            Calculation of PDFs              !!
!!=============================================!!




10000 format (30(e17.7))   
end subroutine STAT_PARTICLE