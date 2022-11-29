program Ellipsoid_Box

use param_phys
use ellipsoid_particle

implicit none

real(kind=8) :: ti, tf, TIME

real(kind=8) :: CPU_INIT, CPU_ELAPSED

logical :: CONT

integer :: NCYCLE

call CPU_TIME(ti)

call READPARAMS

call INIT_PARAM

call INIT_PARTICLE

COLL_FLAG = .false.

CPU_INIT = 0.0
CPU_ELAPSED = 0.0

CONT = .true.
NCYCLE = 1
TIME = 0.0

call CPU_TIME(tf)

CPU_INIT = tf - ti

CPU_ELAPSED = CPU_ELAPSED + CPU_INIT

write(*,10700)100.0*NCYCLE/NCYCLEMAX, CPU_ELAPSED


write(*,*) 'Start the time-loop !!'

do while(CONT)

    call CPU_TIME(ti)

    if(mod(NCYCLE, STAT_CALC) == 0) then

        call STAT_PARTICLE_OLD(NCYCLE, TIME)

        call STAT_PARTICLE(NCYCLE, TIME)

    end if

    call ADV_PARTICLE(NCYCLE)

    if(COLL_FLAG) call COLLISION(NCYCLE)

    if(NCYCLE == NCYCLEMAX) CONT = .false.
    
    TIME = TIME + DT_INIT

    call CPU_TIME(tf)

    CPU_ELAPSED = CPU_ELAPSED + tf - ti

!!- Print percentage of simulation accomplished
    if((NCYCLEMAX>=10).and.(mod(NCYCLE,NCYCLEMAX/10)==0)) write(*,10700)100.*NCYCLE/NCYCLEMAX,CPU_ELAPSED
    NCYCLE  = NCYCLE + 1

end do

write(*,*) ' Time loop ended !!!'


call save_particle

write(*,*) '**************************************'
write(*,*) '       END COMPUTATION FOR NCYCLE'
write(*,*) '**************************************'
write(*,*) '           Cycle =', NCYCLEMAX
write(*,*) '      Ending time=', TIME
write(*,*) '**************************************'
write(*,*)

10700 format (2x,' Computation at ',f6.2,' %, Elapsed time:',f12.3,' s')

end program Ellipsoid_Box