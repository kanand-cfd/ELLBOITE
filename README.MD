# ELLBOITE - Ellipsoids in a Box

## Introduction
ELLBOX is a Lagrangian Particle Simulation code to track the movement of ellipsoidal particles in a cubical box. The ellipsoidal particles can collide with each other and the walls. The linear motion is computed using a simple forward-Euler scheme. The rotational motion is computed using the Euler equation of rotational motion using an RK4 scheme [1]. The orientation of the ellipsoids are tracked using Quaternions. Quaternion integration is performed by directly solving the governing ODE which relates the orientation with angular velocity [3].

Apart from tracking the position and motion of every particle, this code also allows them to collide with each other. At every time-step binary collisions are tracked, i.e. only a pair of particles are considered. The collision pairs are searched in a brute force manner where each particle pair in the box is considered (there are more efficient ways to look for collision pairs, so room for improvement). To detect whether two ellipsoids are in contact a robust contact-detection algorithm has been developed, based on unconstrained optimization algorithm [4, 5] solved using Levenberg Marquardt Method [6] (The details regarding the algorithm will be updated shortly). The code also tracks collisions between the walls and the ellipsoids. Wall-Ellipsoid collision detection is solved using Lagrange Multiplier for a constrained optimization formulation.

Once it is detected that a pair of ellipsoids are in overlap a proper collision impulse is imparted on the ellipsoids using a collision model based on hard-sphere model extended to ellipsiods [2].

## Code Overview

### 
> main.f90          : Main file of the program

> readparams.f90    : Reads the 'params.in' file to initialise key parameters of the simulation

> init_params.f90   : Initialisation of phyiscal parameters, particle arrays to be used in the simulation

> init_particle.f90 : Random initialisation of the particle position, velocity, orientation (by virtue of quaternions) and angular velocity

### In the time loop

> adv_particle.f90  : Advance the position, orientation, linear and angular velocity of the particles for the next time interval

> collision.f90     : Handling collision between particles

> stat_particle.f90 : Calculation of time-averaged statistics (under progress)

### Apart from the time integration of linear velocity and position of ellipsoidal particles, 'adv_particle.f90' subroutine calls 

> box_boundary_condtions.f90 : Detects and calculates collision impulse of ellipsoids and the box walls

> quaternion_integration.f90 : Time integration of quaternion of ellipsoids to update the orientation of ellipsoids

> euler_integration.f90 : Time integration of angular velocities using Euler Equations (uses an RK4 scheme)

> boundary_particle.f90: Periodic Boundary Conditions (if one decides not to use walls : possibly a thing of future)
 
### The subroutine "collision.f90" detects whether a particle pair is overlapping. If an overlap is detected an impulse is imparted on the particle pair. It calls the following subroutines :

> ellipsoid_contact_detection.f90 : Detects whether a particle pair is overlapping or not

> collision_response.f90 : Imparts appropriate impulse on the particle pair once overlap is confirmed

### After the time loop

> save_particle.f90 : At the end of the simulation this subroutine saves the last state (position, linear & angular velocities and orientation) of the particles

### Miscellaneous subroutines

> transform_basis.f90 : Transforms the frame of reference of the ellipsoid from local to global by taking into account the orientation. 

> invert_ndim3_matrix.f90 : invert a 3x3 matrix



## References

1. Goldstein, Herbert, Charles Poole, and John Safko. "Classical mechanics." (2002): 782-783.

2. Tschisgale, Silvio, Louis Thiry, and Jochen Fröhlich. "A constraint-based collision model for Cosserat rods." Archive of Applied Mechanics 89.2 (2019): 167-193.

3. Zhao, F., and B. G. M. Van Wachem. "A novel Quaternion integration approach for describing the behaviour of non-spherical particles." Acta Mechanica 224.12 (2013): 3091-3109.

4. Wellmann, Christian, Claudia Lillie, and Peter Wriggers. “A contact detection algorithm for superellipsoids based on the common-normal concept.” Engineering Computations (2008).

5. Pazouki, A., H. Mazhar, and D. Negrut. “Parallel Contact Detection between Ellipsoids with Applications in Granular Dynamics.” Mathematics and Computers in Simulation, Submit- ted (2010).

6. Minpack library for modern fortran available in the repository: https://github.com/fortran-lang/minpack. 
