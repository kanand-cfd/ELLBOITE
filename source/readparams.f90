subroutine READPARAMS

use param_phys

implicit none

open(unit=200,file='params.in',status='old')

read(200,*) ! comment line
read(200,*) ! comment line
read(200,*) ! comment line
read(200,*) ! comment line
read(200,*) ! comment line

read(200,*) DT_INIT

read(200,*) ! comment line

read(200,*) NCYCLEMAX

read(200,*) ! comment line

read(200,*) LXMAX
read(200,*) LYMAX
read(200,*) LZMAX

read(200,*) ! comment line

read(200,*) NPART_MAX
read(200,*) INIT
read(200,*) ELL_A
read(200,*) lambda

read(200,*) ! comment line

read(200,*) WALL  ! Wall Boundary Condition
read(200,*) EW    ! Coefficient of Normal Restituion
read(200,*) MUW   ! Coefficient of Friction
read(200,*) BETAW ! Coefficient of Tangential Restitution

read(200,*) ! comment line

write(*,*) 'Parameter reading ---> DONE'
    
end subroutine READPARAMS