subroutine STAT_PARTICLE_OLD(NCYCLE, TIME)

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

!real(kind=8) :: NPM

real(kind=8) :: UPM, VPM, WPM

real(kind=8) :: UPFLC, VPFLC, WPFLC

real(kind=8) :: OXM, OYM, OZM

!real(kind=8) :: gl_OXM, gl_OYM, gl_OZM

real(kind=8) :: OXFLC, OYFLC, OZFLC

!real(kind=8) :: gl_OXFLC, gl_OYFLC, gl_OZFLC

!real(kind=8), dimension(ndim, 1) :: global_angular_velocity


MEAN_PART(:) = ZEROP
!global_angular_velocity(:,1) = ZEROP

do I = 1, NPART_MAX

!!----------- Particle Velocities -----------!!
!! <Up>
    MEAN_PART(1) = MEAN_PART(1) + ELLP(I)%UP

!! <Vp>
    MEAN_PART(2) = MEAN_PART(2) + ELLP(I)%VP

!! <Wp>
    MEAN_PART(3) = MEAN_PART(3) + ELLP(I)%WP


!!------ Particle  Angular Velocities ------!!
!! <OMEGAX>
    MEAN_PART(4) = MEAN_PART(4) + ELLP(I)%OMEGAX

!! <OMEGAY>
    MEAN_PART(5) = MEAN_PART(5) + ELLP(I)%OMEGAY

!! <OMEGAZ>
    MEAN_PART(6) = MEAN_PART(6) + ELLP(I)%OMEGAZ


end do

!!- Compute the mean velocities
UPM = MEAN_PART(1)/NPART_MAX
VPM = MEAN_PART(2)/NPART_MAX
WPM = MEAN_PART(3)/NPART_MAX

OXM = MEAN_PART(4)/NPART_MAX
OYM = MEAN_PART(5)/NPART_MAX
OZM = MEAN_PART(6)/NPART_MAX


do I = 1, NPART_MAX

!!- u'p = up-<up>
    UPFLC = ELLP(I)%UP !- UPM

!!- v'p = vp-<vp>
    VPFLC = ELLP(I)%VP !- VPM

!!- w'p = wp-<wp>
    WPFLC = ELLP(I)%WP !- WPM

!!- ox'p = oxp - <oxp>
    OXFLC = ELLP(I)%OMEGAX !- OXM

!!- oy'p = oyp - <oyp>
    OYFLC = ELLP(I)%OMEGAY !- OYM

!!- oz'p = ozp - <ozp>
    OZFLC = ELLP(I)%OMEGAZ !- OZM


!!- <upup>
    MEAN_PART(7) = MEAN_PART(7) + UPFLC*UPFLC

!!- <vpvp>
    MEAN_PART(8) = MEAN_PART(8) + VPFLC*VPFLC

!!- <wpwp>
    MEAN_PART(9) = MEAN_PART(9) + WPFLC*WPFLC

!!- <oxp.oxp>
    MEAN_PART(10) = MEAN_PART(10) + OXFLC*OXFLC

!!- <oyp.oyp>
    MEAN_PART(11) = MEAN_PART(11) + OYFLC*OYFLC

!!- <ozp.ozp>
    MEAN_PART(12) = MEAN_PART(12) + OZFLC*OZFLC


!!- <up*vp>
    MEAN_PART(13) = MEAN_PART(13) + UPFLC*VPFLC

!!- <up*wp>
    MEAN_PART(14) = MEAN_PART(14) + UPFLC*WPFLC

!!- <vp*wp>
    MEAN_PART(15) = MEAN_PART(15) + VPFLC*WPFLC


!!- <oxp.oxp>
    MEAN_PART(16) = MEAN_PART(16) + OXFLC*OYFLC

!!- <oyp.oyp>
    MEAN_PART(17) = MEAN_PART(17) + OXFLC*OZFLC

!!- <ozp.ozp>
    MEAN_PART(18) = MEAN_PART(18) + OYFLC*OZFLC

end do


!!- qp - Linear Velocity
MEAN_PART(19) = 0.5*(MEAN_PART(7) + &
                     MEAN_PART(8) + &
                     MEAN_PART(9))

!!- qp - Angular Velocity
MEAN_PART(20) = 0.5*(IPXX*MEAN_PART(10) + &
                     IPYY*MEAN_PART(11) + &
                     IPZZ*MEAN_PART(12))


MEAN_PART(:) = MEAN_PART(:)/NPART_MAX


if(mod(NCYCLE, STAT_WRITE)==0) then 

    write(500, 10000) &
    TIME, &
    MEAN_PART(1), &   !! <Up>
    MEAN_PART(2), &   !! <Vp>
    MEAN_PART(3), &   !! <Wp>
    MEAN_PART(4), &   !! <OMEGAX>
    MEAN_PART(5), &   !! <OMEGAY>
    MEAN_PART(6), &   !! <OMEGAZ>
    MEAN_PART(7), &   !!- <upup>
    MEAN_PART(8), &   !!- <vpvp>
    MEAN_PART(9), &   !!- <wpwp>
    MEAN_PART(10), &  !!- <oxp.oxp>
    MEAN_PART(11), &  !!- <oyp.oyp>
    MEAN_PART(12), &  !!- <ozp.ozp>
    MEAN_PART(13), &  !!- <upvp>
    MEAN_PART(14), &  !!- <upwp>
    MEAN_PART(15), &  !!- <vpwp>
    MEAN_PART(16), &  !!- <oxp.oyp>
    MEAN_PART(17), &  !!- <oxp.ozp>
    MEAN_PART(18), &  !!- <oyp.ozp>  
    MEAN_PART(19), &  !!- qp - Linear Velocity
    MEAN_PART(20)!, &  !!- qp - Angular Velocity      

end if 

 
10000 format (30(e17.7))   
end subroutine STAT_PARTICLE_OLD