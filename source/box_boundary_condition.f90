subroutine BOX_BOUNDARY_CONDITION(part_ellp)

use param_phys
use mod_quaternion
use ellipsoid_particle

implicit none

!=================================================!
! 		            Input Variables               !
!=================================================!
type(ELL_PART), intent(inout) :: part_ellp
!=================================================!

! Wall Parameters
real(kind=8), dimension(ndim,1) :: wall_normal ! normal vector
real(kind=8), dimension(ndim,1) :: wall_pos ! Position vector

! Closest point on plane from Ellipsoid
real(kind=8), dimension(ndim,1) :: plane_pt

! Maximum Encroachment point on Ellipsoid
real(kind=8), dimension(ndim,1) :: ellip_pt

! Depth of penetration
real(kind=8) :: depth

! Margin function for Ellipsoid
real(kind=8) :: margin

! Impact Arm
real(kind=8), dimension(ndim,1) :: impact_arm
!=================================================!
!=================================================!

! Call subroutine for closest point between ellipsoid and plane !
! by passing the plane location and particle location. Check if !
! the particle center  has position greater or  lesser than the !
! closest point between the ellipsoid and the plane.  Apply the !
! collision force accordingly. Repeat the procedure for all the !
! 6 walls.														!

!=====================WALL#1======================!
! Wall at (0,0,0) with normal (1,0,0)
wall_pos(1,1) = 0.0
wall_pos(2,1) = 0.0
wall_pos(3,1) = 0.0

wall_normal(1,1) = 1.0
wall_normal(2,1) = 0.0
wall_normal(3,1) = 0.0

! call closest point
call closest_point_plane(part_ellp%XP,part_ellp%YP,part_ellp%ZP, &
                         part_ellp%ELLQUAT,wall_pos,wall_normal, &
                         plane_pt,margin) 

if(margin .LE. 1.0) then

    call max_encroachment_ellipsoid(part_ellp%XP,part_ellp%YP,part_ellp%ZP, &
                                    part_ellp%ELLQUAT,wall_normal,ellip_pt)

    ! Depth
    depth = abs(plane_pt(1,1) - ellip_pt(1,1))

    ! Trace the particle back to state of contact
    !part_ellp%XP = part_ellp%XP - depth

    !ellip_pt(1,1) = ellip_pt(1,1) - depth

    ! Calculate impact arm
    impact_arm(1,1) = abs(part_ellp%XP - ellip_pt(1,1))
    impact_arm(2,1) = abs(part_ellp%YP - ellip_pt(2,1))
    impact_arm(3,1) = abs(part_ellp%ZP - ellip_pt(3,1))

    ! Calculate number of collisions
    ! Calculate PDF of velocities before rebound
    !call PDF_BEFORE(part_ellp%UP, part_ellp%VP, part_ellp%WP, &
    !                part_ellp%OMEGAX, part_ellp%OMEGAY, part_ellp%OMEGAZ)

    ! calculate impulse
    call wall_ellipsoid_rebound(part_ellp%UP,part_ellp%VP,part_ellp%WP, &
                                part_ellp%OMEGAX,part_ellp%OMEGAY,part_ellp%OMEGAZ, &
                                part_ellp%ELLQUAT,impact_arm,wall_normal)

    
    !stop
    ! Calculate PDF of velocities after rebound
    !call PDF_AFTER(part_ellp%)

end if     
!=================================================!



!=====================WALL#2======================!
! Wall at (LX,0,0) with normal (-1,0,0)
wall_pos(1,1) = LXMAX
wall_pos(2,1) = 0.0
wall_pos(3,1) = 0.0

wall_normal(1,1) = -1.0
wall_normal(2,1) = 0.0
wall_normal(3,1) = 0.0

call closest_point_plane(part_ellp%XP,part_ellp%YP,part_ellp%ZP, &
                         part_ellp%ELLQUAT,wall_pos,wall_normal, &
                         plane_pt,margin) 

if(margin .LE. 1) then

    call max_encroachment_ellipsoid(part_ellp%XP,part_ellp%YP,part_ellp%ZP, &
                                    part_ellp%ELLQUAT,wall_normal,ellip_pt)

    ! Depth
    depth = abs(plane_pt(1,1) - ellip_pt(1,1))

    ! Trace the particle back to state of contact
    !part_ellp%XP = part_ellp%XP - depth

    !ellip_pt(1,1) = ellip_pt(1,1) - depth

    ! Calculate impact arm
    impact_arm(1,1) = abs(part_ellp%XP - ellip_pt(1,1))
    impact_arm(2,1) = abs(part_ellp%YP - ellip_pt(2,1))
    impact_arm(3,1) = abs(part_ellp%ZP - ellip_pt(3,1))

    ! Calculate number of collisions
    ! Calculate PDF of velocities before rebound

    ! calculate impulse
    call wall_ellipsoid_rebound(part_ellp%UP,part_ellp%VP,part_ellp%WP, &
                                part_ellp%OMEGAX,part_ellp%OMEGAY,part_ellp%OMEGAZ, &
                                part_ellp%ELLQUAT,impact_arm,wall_normal)
    !stop
    ! Calculate PDF of velocities after rebound
    
end if
!=================================================!     



!=====================WALL#3======================!
! Wall at (0,0,0) with normal (0,1,0)
wall_pos(1,1) = 0.0
wall_pos(2,1) = 0.0
wall_pos(3,1) = 0.0

wall_normal(1,1) = 0.0
wall_normal(2,1) = 1.0
wall_normal(3,1) = 0.0

! call subroutine

call closest_point_plane(part_ellp%XP,part_ellp%YP,part_ellp%ZP, &
                         part_ellp%ELLQUAT,wall_pos,wall_normal, &
                         plane_pt, margin) 

if(margin .LE. 1) then

    call max_encroachment_ellipsoid(part_ellp%XP,part_ellp%YP,part_ellp%ZP, &
                                    part_ellp%ELLQUAT,wall_normal,ellip_pt)

    ! Depth
    depth = abs(plane_pt(2,1) - ellip_pt(2,1))

    ! Trace the particle back to state of contact
    !part_ellp%YP = part_ellp%YP - depth

    !ellip_pt(2,1) = ellip_pt(2,1) - depth

    ! Calculate impact arm
    impact_arm(1,1) = abs(part_ellp%XP - ellip_pt(1,1))
    impact_arm(2,1) = abs(part_ellp%YP - ellip_pt(2,1))
    impact_arm(3,1) = abs(part_ellp%ZP - ellip_pt(3,1))

    ! Calculate number of collisions
    ! Calculate PDF of velocities before rebound

    ! calculate impulse
    call wall_ellipsoid_rebound(part_ellp%UP,part_ellp%VP,part_ellp%WP, &
                                part_ellp%OMEGAX,part_ellp%OMEGAY,part_ellp%OMEGAZ, &
                                part_ellp%ELLQUAT,impact_arm,wall_normal)
    !stop
    ! Calculate PDF of velocities after rebound

end if   
!=================================================!  



!=====================WALL#4======================!
! Wall at (0,LY,0) with normal (0,-1,0)
wall_pos(1,1) = 0.0
wall_pos(2,1) = LYMAX
wall_pos(3,1) = 0.0

wall_normal(1,1) = 0.0
wall_normal(2,1) = -1.0
wall_normal(3,1) = 0.0

! call subroutine
call closest_point_plane(part_ellp%XP,part_ellp%YP,part_ellp%ZP, &
                         part_ellp%ELLQUAT,wall_pos,wall_normal, &
                         plane_pt,margin) 

if(margin < 1) then

    call max_encroachment_ellipsoid(part_ellp%XP,part_ellp%YP,part_ellp%ZP, &
                                    part_ellp%ELLQUAT,wall_normal,ellip_pt)

    ! Depth
    depth = abs(plane_pt(2,1) - ellip_pt(2,1))

    ! Trace the particle back to state of contact
    !part_ellp%YP = part_ellp%YP - depth

    !ellip_pt(2,1) = ellip_pt(2,1) - depth

    ! Calculate impact arm
    impact_arm(1,1) = abs(part_ellp%XP - ellip_pt(1,1))
    impact_arm(2,1) = abs(part_ellp%YP - ellip_pt(2,1))
    impact_arm(3,1) = abs(part_ellp%ZP - ellip_pt(3,1))

    ! calculate impulse
    !write(*,*) 'depth ', depth
    !write(*,*) 'Updated max encorachment point', ellip_pt
    !write(*,*) 'Ellipsoid Center ', part_ellp%XP, part_ellp%YP, part_ellp%ZP
    !write(*,*) 'impact_arm ', impact_arm
    
    ! Calculate number of collisions
    ! Calculate PDF of velocities before rebound

    ! calculate impulse
    call wall_ellipsoid_rebound(part_ellp%UP,part_ellp%VP,part_ellp%WP, &
                                part_ellp%OMEGAX,part_ellp%OMEGAY,part_ellp%OMEGAZ, &
                                part_ellp%ELLQUAT,impact_arm,wall_normal)
    !stop
    ! Calculate PDF of velocities after rebound
    
end if
!=================================================!     


!=====================WALL#5======================!
! Wall at (0,0,0) with normal (0,0,1)
wall_pos(1,1) = 0.0
wall_pos(2,1) = 0.0
wall_pos(3,1) = 0.0

wall_normal(1,1) = 0.0
wall_normal(2,1) = 0.0
wall_normal(3,1) = 1.0

! call subroutine
call closest_point_plane(part_ellp%XP,part_ellp%YP,part_ellp%ZP, &
                         part_ellp%ELLQUAT,wall_pos,wall_normal, &
                         plane_pt,margin) 

if(margin .LE. 1) then

    call max_encroachment_ellipsoid(part_ellp%XP,part_ellp%YP,part_ellp%ZP, &
                                    part_ellp%ELLQUAT,wall_normal,ellip_pt)

    ! Depth
    depth = abs(plane_pt(3,1) - ellip_pt(3,1))

    ! Trace the particle back to state of contact
    !part_ellp%ZP = part_ellp%ZP - depth

    !ellip_pt(3,1) = ellip_pt(3,1) - depth

    ! Calculate impact arm
    impact_arm(1,1) = abs(part_ellp%XP - ellip_pt(1,1))
    impact_arm(2,1) = abs(part_ellp%YP - ellip_pt(2,1))
    impact_arm(3,1) = abs(part_ellp%ZP - ellip_pt(3,1))

    ! Calculate number of collisions
    ! Calculate PDF of velocities before rebound

    ! calculate impulse
    call wall_ellipsoid_rebound(part_ellp%UP,part_ellp%VP,part_ellp%WP, &
                                part_ellp%OMEGAX,part_ellp%OMEGAY,part_ellp%OMEGAZ, &
                                part_ellp%ELLQUAT,impact_arm,wall_normal)
    !stop
    ! Calculate PDF of velocities after rebound

end if     
!=================================================!

!=====================WALL#6======================!
! Wall at (0,0,LZ) with normal (0,0,-1)
wall_pos(1,1) = 0.0
wall_pos(2,1) = 0.0
wall_pos(3,1) = LZMAX

wall_normal(1,1) = 0.0
wall_normal(2,1) = 0.0
wall_normal(3,1) = -1.0

! call subroutine
call closest_point_plane(part_ellp%XP,part_ellp%YP,part_ellp%ZP, &
                         part_ellp%ELLQUAT,wall_pos,wall_normal, &
                         plane_pt,margin) 

if(margin .LE. 1) then

    call max_encroachment_ellipsoid(part_ellp%XP,part_ellp%YP,part_ellp%ZP, &
                                    part_ellp%ELLQUAT,wall_normal,ellip_pt)

    ! Depth
    depth = abs(plane_pt(3,1) - ellip_pt(3,1))

    ! Trace the particle back to state of contact
    !part_ellp%ZP = part_ellp%ZP - depth

    !ellip_pt(3,1) = ellip_pt(3,1) - depth

    ! Calculate impact arm
    impact_arm(1,1) = abs(part_ellp%XP - ellip_pt(1,1))
    impact_arm(2,1) = abs(part_ellp%YP - ellip_pt(2,1))
    impact_arm(3,1) = abs(part_ellp%ZP - ellip_pt(3,1))

    ! Calculate number of collisions
    ! Calculate PDF of velocities before rebound

    ! calculate impulse
    call wall_ellipsoid_rebound(part_ellp%UP,part_ellp%VP,part_ellp%WP, &
                                part_ellp%OMEGAX,part_ellp%OMEGAY,part_ellp%OMEGAZ, &
                                part_ellp%ELLQUAT,impact_arm,wall_normal)
    !stop
    ! Calculate PDF of velocities after rebound

end if     
!=================================================!

end subroutine BOX_BOUNDARY_CONDITION