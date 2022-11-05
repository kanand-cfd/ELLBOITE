subroutine save_particle

use param_phys
use ellipsoid_particle
use mod_quaternion
    
implicit none

character(len=40) :: FILENAME

integer :: I, J

J = 1

write(FILENAME, 10100) 'part_p', J,'.end'

open(unit=150, file= trim(FILENAME), status='replace', form='unformatted')

write(150) NPART_MAX
write(150)(ELLP(I)%XP, I=1, NPART_MAX)
write(150)(ELLP(I)%YP, I=1, NPART_MAX)
write(150)(ELLP(I)%ZP, I=1, NPART_MAX)

write(150)(ELLP(I)%UP, I=1, NPART_MAX)
write(150)(ELLP(I)%VP, I=1, NPART_MAX)
write(150)(ELLP(I)%WP, I=1, NPART_MAX)

write(150)(ELLP(I)%OMEGAX, I=1, NPART_MAX)
write(150)(ELLP(I)%OMEGAY, I=1, NPART_MAX)
write(150)(ELLP(I)%OMEGAZ, I=1, NPART_MAX)

write(150)(ELLP(I)%ELLQUAT%a, I =1, NPART_MAX)
write(150)(ELLP(I)%ELLQUAT%b, I =1, NPART_MAX)
write(150)(ELLP(I)%ELLQUAT%c, I =1, NPART_MAX)
write(150)(ELLP(I)%ELLQUAT%d, I =1, NPART_MAX)

close(150)

write(*,*) 'Final particle position and velocity --> Binary File'

10100 format(A,I2.2,A)    

end subroutine save_particle