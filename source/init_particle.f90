subroutine init_particle

use param_phys
use ellipsoid_particle
use mod_quaternion

implicit none

integer :: I, J

real(kind=8):: XRAND
real(kind=8), dimension(3) :: PART_MEAN

real(kind=8) :: ti, tf

character(len=40) :: FILENAME

integer :: NP_READ

call CPU_TIME(ti)

PART_MEAN(:) = ZEROP


!=======================================!
!      Particle Position Initiation     !
!=======================================!
if(INIT == 0) then

    do I = 1, NPART_MAX

        call random_number(XRAND)
        ELLP(I)%XP = XRAND*(LXMAX)! + DX)

        call random_number(XRAND)
        ELLP(I)%YP = XRAND*(LYMAX)! + DY)

        call random_number(XRAND)
        ELLP(I)%ZP = XRAND*(LZMAX)! + DZ)


    end do

    write(*,*)'Particle position initiation --> Random'
    !stop

!=======================================!
!      Particle Orientation Initiation  !
!=======================================!
    do I = 1, NPART_MAX
    
        call random_number(XRAND)
        ELLP(I)%ELLQUAT%a = XRAND

        call random_number(XRAND)
        ELLP(I)%ELLQUAT%b = XRAND 

        call random_number(XRAND)
        ELLP(I)%ELLQUAT%c = XRAND
 
        call random_number(XRAND)
        ELLP(I)%ELLQUAT%d = XRAND

        ELLP(I)%ELLQUAT = unit_quat(ELLP(I)%ELLQUAT)

    end do

    write(*,*)'Particle Orientation initiation --> Random'


!=======================================!
!      Particle Velocity Initiation     !
!=======================================!
    do I = 1, NPART_MAX
    
        call random_number(XRAND)
        ELLP(I)%UP = XRAND !- 0.5
    
        call random_number(XRAND)
        ELLP(I)%VP = XRAND !- 0.5

        call random_number(XRAND)
        ELLP(I)%WP = XRAND !- 0.5

        PART_MEAN(1) = PART_MEAN(1) + ELLP(I)%UP
        PART_MEAN(2) = PART_MEAN(2) + ELLP(I)%VP
        PART_MEAN(3) = PART_MEAN(3) + ELLP(I)%WP

    end do

    PART_MEAN(:) = PART_MEAN(:)/NPART_MAX

    do I = 1, NPART_MAX

        ELLP(I)%UP = ELLP(I)%UP - PART_MEAN(1)

        ELLP(I)%VP = ELLP(I)%VP - PART_MEAN(2)

        ELLP(I)%WP = ELLP(I)%WP - PART_MEAN(3)

        !write(*,*) I, ELLP(I)%XP, ELLP(I)%YP, ELLP(I)%ZP

    end do

    write(*,*)'Particle velocity initiation: --> Random'

!=======================================!
! Particle Angular Velocity Initiation  !
!=======================================!
    do I = 1, NPART_MAX

    !call random_number(XRAND)
        ELLP(I)%OMEGAX = 0.0!1.0e04*XRAND

    !call random_number(XRAND)
        ELLP(I)%OMEGAY = 0.0!1.0e04*XRAND    

    !call random_number(XRAND)
        ELLP(I)%OMEGAZ = 0.0!1.0e04*XRAND

    end do

    write(*,*)'Particle angular velocity initiation: --> Uniform (0.0)'

!==================================================================!
!               Particle Initiation from Binary File               !
!==================================================================!
elseif(INIT == 1) then

    write(*,*)'Particle position initiation: Read from single binary file'

    J = 1

    write(FILENAME, 10100)'part_p',J,'.ini'

    open(unit=150, file=trim(FILENAME), status='old', form = 'unformatted')

    read(150)NP_READ

    if((NP_READ==NPART_MAX) .eqv. .FALSE.) then
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)'!!                     ERROR                    !!'
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)'!! file     : initiation_particle_position.f90'
        write(*,*)'!!'
        write(*,*)'!!   NP_READ=',NP_READ
        write(*,*)'!!   NPMAX=', NPART_MAX
        write(*,*)'!!'
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        stop
    end if


    NPART_MAX = NP_READ

    read(150)(ELLP(I)%XP, I=1, NP_READ)
    read(150)(ELLP(I)%YP, I=1, NP_READ)
    read(150)(ELLP(I)%ZP, I=1, NP_READ)

    read(150)(ELLP(I)%UP, I=1, NP_READ)
    read(150)(ELLP(I)%VP, I=1, NP_READ)
    read(150)(ELLP(I)%WP, I=1, NP_READ)

    read(150)(ELLP(I)%OMEGAX, I=1, NP_READ)
    read(150)(ELLP(I)%OMEGAY, I=1, NP_READ)
    read(150)(ELLP(I)%OMEGAZ, I=1, NP_READ)

    read(150)(ELLP(I)%ELLQUAT%a, I =1, NP_READ)
    read(150)(ELLP(I)%ELLQUAT%b, I =1, NP_READ)
    read(150)(ELLP(I)%ELLQUAT%c, I =1, NP_READ)
    read(150)(ELLP(I)%ELLQUAT%d, I =1, NP_READ)


    close(150)

    write(*,*)'Particle position initiation: Read from file --> OK'

end if

call CPU_TIME(tf)

    write(*,*) 'Time elapsed in Initialization', tf - ti 


10100 format(A,I2.2,A)
!stop
end subroutine init_particle