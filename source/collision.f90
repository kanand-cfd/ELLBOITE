subroutine COLLISION(NCYCLE)

use param_phys
use mod_quaternion
use ellipsoid_particle

implicit none

integer, intent(in) :: NCYCLE

integer :: I, J, NO

logical :: flag

! Vector linking particle centers
real(kind=8), dimension(ndim, 1) :: VPC

real(kind=8) :: NRM_VPC

! Ellipsoid 1
real(kind=8) :: UP1, VP1, WP1 ! Velocity of contact points
real(kind=8) :: RX1, RY1, RZ1 ! Dist. of ct Pt from center

! Ellipsoid 2
real(kind=8) :: UP2, VP2, WP2 ! Velocity of contact points
real(kind=8) :: RX2, RY2, RZ2 ! Dist. of ct Pt from center

! Relative velocity of particles
real(kind=8) :: VRELX, VRELY, VRELZ


real(kind=8) :: NVREL

! Output of the contact detection algorithm
! Point on Ellipsoid 1
real(kind=8), dimension(ndim, 1) :: EPT1

! Point on Ellipsoid 2
real(kind=8), dimension(ndim, 1) :: EPT2

! Common normal
real(kind=8), dimension(ndim, 1) :: cnormal
!====================================!

ELLP(:)%COLOR = 1.0
COLLISION_COUNT = 0
ERROR_COUNT = 0
AVG_ERROR = 0.0

CONTTACT = 0
NO_CONTACT = 0

NO = 0

do I = 1, NPART_MAX
    do J = I+1, NPART_MAX

        
        flag = .FALSE.

        !!- Allow only one collision per particle
        if(ELLP(I)%COLOR == 1.0 .and. ELLP(J)%COLOR == 1.0) then

            call ELLIPSOID_CONTACT_DETECTION(ELLP(I)%XP, ELLP(I)%YP, ELLP(I)%ZP, ELLP(I)%ELLQUAT, &
                                                 ELLP(J)%XP, ELLP(J)%YP, ELLP(J)%ZP, ELLP(J)%ELLQUAT, & 
                                                 flag, EPT1, EPT2, cnormal)

            ! Apply Pope's algorithm to detect whether we have intersection
            !call ELLIPSOID_OVERLAP_DETECTION(ELLP(I)%XP, ELLP(I)%YP, ELLP(I)%ZP, ELLP(I)%ELLQUAT, &
            !                                 ELLP(J)%XP, ELLP(J)%YP, ELLP(J)%ZP, ELLP(J)%ELLQUAT, & 
            !                                 flag)


            ! if true -->  call collision detection
            if (flag) then
                !write(*,*) flag
                

                ! For 2 overlapped ellipsoids we need to figure out whether they are approaching
                ! or going away, this can be accurately done by considering the point on the  
                ! ellipsoid surface having the maximum encroachment into the other ellipsoid

                ! Calculate the vector linking the max. encroachment points                
                VPC(1,1) = EPT1(1,1) - EPT2(1,1)
                VPC(2,1) = EPT1(2,1) - EPT2(2,1)
                VPC(3,1) = EPT1(3,1) - EPT2(3,1)

                NRM_VPC = sqrt(VPC(1,1)*VPC(1,1) + VPC(2,1)*VPC(2,1) + VPC(3,1)*VPC(3,1)) 

                if (NRM_VPC /= 0.0) VPC = VPC / NRM_VPC

                ! Calculate the velocity of the max. encroachment points
                ! Particle 1
                RX1 = abs(EPT1(1,1) - ELLP(I)%XP)
                RY1 = abs(EPT1(2,1) - ELLP(I)%YP)
                RZ1 = abs(EPT1(3,1) - ELLP(I)%ZP)

                UP1 = ELLP(I)%UP + (ELLP(I)%OMEGAY*RZ1 - ELLP(I)%OMEGAZ*RY1)
                VP1 = ELLP(I)%VP + (ELLP(I)%OMEGAZ*RX1 - ELLP(I)%OMEGAX*RZ1)
                WP1 = ELLP(I)%WP + (ELLP(I)%OMEGAX*RY1 - ELLP(I)%OMEGAY*RX1)

                ! Particle 2
                RX2 = abs(EPT2(1,1) - ELLP(J)%XP)
                RY2 = abs(EPT2(2,1) - ELLP(J)%YP)
                RZ2 = abs(EPT2(3,1) - ELLP(J)%ZP)              

                UP2 = ELLP(J)%UP + (ELLP(J)%OMEGAY*RZ2 - ELLP(J)%OMEGAZ*RY2)
                VP2 = ELLP(J)%VP + (ELLP(J)%OMEGAZ*RX2 - ELLP(J)%OMEGAX*RZ2)
                WP2 = ELLP(J)%WP + (ELLP(J)%OMEGAX*RY2 - ELLP(J)%OMEGAY*RX2)

                ! Calculate the relative velocity
                VRELX = UP2 - UP1
                VRELY = VP2 - VP1
                VRELZ = WP2 - WP1

                ! Dot Product of Relative Velocity and Vector linking encroachment points)
                NVREL = VRELX*VPC(1,1) + VRELY*VPC(2,1) + VRELZ*VPC(3,1)

                ! if particles are approaching
                if(NVREL < ZEROP) then

                    !write(*,*) 'particles are approaching'
                    !write(*,*) 'Before Collision '
                    !write(*,*) ELLP(I)%UP, ELLP(I)%VP, ELLP(I)%WP
                    !write(*,*) ELLP(J)%UP, ELLP(J)%VP, ELLP(J)%WP

                    ! Calculate the impulse on the particles and Update the velocities
                    if(COLL_FLAG) call  COLLISION_RESPONSE(ELLP(I), ELLP(J), EPT1, EPT2, cnormal, I, J)

                    !write(*,*) 'After Collision '
                    !write(*,*) ELLP(I)%UP, ELLP(I)%VP, ELLP(I)%WP
                    !write(*,*) ELLP(J)%UP, ELLP(J)%VP, ELLP(J)%WP

                    !write(*,*) COLLISION_COUNT
                    COLLISION_COUNT = COLLISION_COUNT + 1 

                    ! Update COLOR
                    ELLP(I)%COLOR = -1.0
                    ELLP(J)%COLOR = -1.0


                end if ! NVREL < ZEROP

            end if !if overlapping == TRUE

        end if ! if(ELLP(I)%COLOR == 1 .and. ELLP(J)%COLOR == 1)

    end do ! J = 1, NPART_MAX

end do ! I = 1, NPART_MAX

!write(*, *) ' ' 
write(*, *) 'Collision count in cycle ', NCYCLE, ' =',  COLLISION_COUNT
if(COLL_FLAG) then
    write(*, *) 'Collision erros in cycle ', NCYCLE, ' =',  ERROR_COUNT
    if(ERROR_COUNT > 0) AVG_ERROR = AVG_ERROR / ERROR_COUNT
    write(*, *) 'Average ERROR   in cycle ', NCYCLE, ' =', AVG_ERROR
end if 

!write(*, *) '               Contact                = ', CONTTACT
!write(*, *) '               No Contact             = ', NO_CONTACT
!write(*, *) '               No Contacts collided   = ', NO
!write(*, *) ' '
!write(*,*) 'End of Collision step'
!stop    
end subroutine COLLISION