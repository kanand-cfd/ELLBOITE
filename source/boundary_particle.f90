!!====================================================================
!!
!!   3-dimensional periodical boundary for particles
!!
!!====================================================================

subroutine BOUNDARY_PARTICLE

!!====================================================================
!!
!!
!!====================================================================

use param_phys
use ellipsoid_particle

implicit none

!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------

!- Index
integer :: I
!---------------------------------------------------------------------

do I = 1, NPART_MAX

!!- x-component
  if(ELLP(I)%XP > LXMAX) ELLP(I)%XP = ELLP(I)%XP - LXMAX
  if(ELLP(I)%XP < 0.   ) ELLP(I)%XP = ELLP(I)%XP + LXMAX

!!- y-component
  if(ELLP(I)%YP > LYMAX) ELLP(I)%YP = ELLP(I)%YP - LYMAX
  if(ELLP(I)%YP < 0.   ) ELLP(I)%YP = ELLP(I)%YP + LYMAX

!!- z-component
  if(ELLP(I)%ZP > LZMAX) ELLP(I)%ZP = ELLP(I)%ZP - LZMAX
  if(ELLP(I)%ZP < 0.   ) ELLP(I)%ZP = ELLP(I)%ZP + LZMAX


end do


!------------------------------------------------------------------
end subroutine BOUNDARY_PARTICLE
