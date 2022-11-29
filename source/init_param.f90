subroutine init_param
    
use param_phys
use ellipsoid_particle

implicit none

!integer :: NUMFILE

character(len=40) :: FILENAME

real(kind=8) :: Volume_Box, Volume_Particle, alpha
    
allocate(ELLP(NPART_MAX))

! Semi Minor axes
ELL_B = ELL_A/lambda
ELL_C = ELL_B

! Moment of Inertia of the Ellipsoid
IPXX = (ELL_B**2 + ELL_C**2)/5.0
IPYY = (ELL_A**2 + ELL_C**2)/5.0
IPZZ = (ELL_A**2 + ELL_B**2)/5.0

! Volume of Ellipsoidal Particle = 4/3*pi*a*b*c
Volume_Particle = (4.0/3.0)*ppi*ELL_A*ELL_B*ELL_C

! Volume of Box
Volume_Box = LXMAX*LYMAX*LZMAX

! Volume fraction
alpha = (NPART_MAX * Volume_Particle)/Volume_Box

! Parameters for statistics caluclation and writing data
STAT_CALC = 1000
STAT_WRITE = 100



write(*,*) '  '
print*, '#####################################'
write(*,10000) ' Dimension of the cubic box :', LXMAX
write(*,10000) ' Volume of the box :', Volume_Box 
write(*,*) 'Total number of Particles :', NPART_MAX
write(*,10000) ' Volume fraction of Particles : ',alpha 
print*, '#####################################'


write(*,*) '  '
print*, '#####################################'
write(*,10000) ' Semi Major Axis, a =', ELL_A
write(*,10000) ' Semi Minor Axes, b = c =', ELL_B 
write(*,10000) ' Aspect Ratio = ', lambda
print*, '#####################################'

print*, ' Moment of Inertia'
write(*,10000) ' I_pxx = ', IPXX
write(*,10000) ' I_pyy = ', IPYY
write(*,10000) ' I_pzz = ', IPZZ
print*, '#####################################'

print*,' Wall Properties'
write(*,10000) ' e_w =', EW
!write(*,10000) ' beta_w =', BETAW
write(*,10000) ' mu_w =', MUW
print*, '#####################################'

! Initialise filename
FILENAME = 'ell_part_translation.stat'

open(unit=300, file=trim(FILENAME), status='replace')
write(300,20000)
write(*,10700)trim(FILENAME)

FILENAME = 'ell_part_angular.stat'

open(unit=400, file=trim(FILENAME), status='replace')
write(400,30000)
write(*,10700)trim(FILENAME)

FILENAME = 'ell_part.stat'

open (unit=500, file=trim(FILENAME), status='replace')
write(500,40000)
write(*,10700)trim(FILENAME)
! Fin d'initialisation

10000 format(A, ES14.7)
10700 format (1x,'   +   ',A)


20000 format('# t, <up>, <vp>, <wp>, <up.up>, <vp.vp>, <wp.wp>, <up.vp>, <up.wp>, <vp.wp>, qp')
30000 format('# t, <omegax>, <omegay>, <omegaz>, <oxp.oxp>, <oyp.oyp>, <ozp.ozp>, <oxp.oyp>, <oxp.ozp>, <oyp.ozp>, qp_omg') 
40000 format('# t, <up>, <vp>, <wp>, <up.up>, <vp.vp>, <wp.wp>, <up.vp>, <up.wp>, <vp.wp>, qp, <omegax>, <omegay>, <omegaz>, <oxp.oxp>, <oyp.oyp>, <ozp.ozp>, <oxp.oyp>, <oxp.ozp>, <oyp.ozp>, qp_omg')   
end subroutine init_param