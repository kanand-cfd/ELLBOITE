module param_phys

implicit none

integer, parameter :: ndim=3

real(kind=8), parameter :: PPI = 4.0*ATAN(1.0)
real(kind=8), parameter :: TWOPI = 8.0*ATAN(1.0)

real(kind=8), parameter :: ZEROP = 1.0e-33

! Dimensions of the box
real(kind=8):: LXMAX
real(kind=8):: LYMAX
real(kind=8):: LZMAX

! Total number of particles
integer :: NPART_MAX

! Particle Initiation
integer :: INIT 

! Dimensions of the ellipsoidal particles
real(kind=8):: ELL_A
real(kind=8):: ELL_B
real(kind=8):: ELL_C

! Aspect Ratio
real(kind=8):: lambda

! Moment of Inertia
real(kind=8):: IPXX
real(kind=8):: IPYY
real(kind=8):: IPZZ


! Simulation parameters
integer :: NCYCLEMAX
real(kind=8) :: DT_INIT

! Wall Properties
logical :: WALL
real(kind=8) :: EW
real(kind=8) :: MUW
real(kind=8) :: BETAW

! Statistics
logical :: COLL_FLAG
integer :: NREBOUND, COLLISION_COUNT, ERROR_COUNT, WALL_COUNT, CONTTACT, NO_CONTACT
real(kind=8) :: AVG_ERROR
integer, parameter :: NPDF = 1000
! PDF translational velocity
real(kind=8), dimension(NPDF) :: PDFUP_BEFORE, PDFVP_BEFORE, PDFWP_BEFORE
real(kind=8), dimension(NPDF) :: PDFUP_AFTER, PDFVP_AFTER, PDFWP_AFTER
! PDF global rotational velocity
real(kind=8), dimension(NPDF) :: PDFOXP_BEFORE, PDFOYP_BEFORE, PDFOZP_BEFORE
real(kind=8), dimension(NPDF) :: PDFOXP_AFTER, PDFOYP_AFTER, PDFOZP_AFTER

end module param_phys