subroutine ADV_PARTICLE(NCYCLE)

use param_phys
use ellipsoid_particle
use mod_quaternion
    
implicit none

integer, intent(in) :: NCYCLE

integer :: I

WALL_COUNT = 0

do I = 1, NPART_MAX

    call BOX_BOUNDARY_CONDITION(ELLP(I))

    ELLP(I)%XP = ELLP(I)%XP + DT_INIT*ELLP(I)%UP
    ELLP(I)%YP = ELLP(I)%YP + DT_INIT*ELLP(I)%VP
    ELLP(I)%ZP = ELLP(I)%ZP + DT_INIT*ELLP(I)%WP

    call QUATERNION_INTEGRATION(ELLP(I)%ELLQUAT, ELLP(I)%OMEGAX, ELLP(I)%OMEGAY, ELLP(I)%OMEGAZ)

    call EULER_INTEGRATION(ELLP(I)%OMEGAX, ELLP(I)%OMEGAY, ELLP(I)%OMEGAZ)

end do

if(WALL .eqv. .FALSE.) call BOUNDARY_PARTICLE

write(*, *) 'Wall  Rebound  in  cycle ', NCYCLE, ' =',  WALL_COUNT
    
end subroutine ADV_PARTICLE