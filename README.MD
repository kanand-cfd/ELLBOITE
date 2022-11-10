## Ellipsoids in a Box - ELLBOITE 

ELLBOX is a Lagrangian Particle Simulation code to track the movement of ellipsoidal particles in a cubical box. The ellipsoidal particles can collide with each other and the walls. The linear motion is computed using a simple forward-Euler scheme. The rotational motion is computed using the Euler equation of rotational motion using an RK4 scheme [1]. The orientation of the ellipsoids are tracked using Quaternions. Quaternion integration is performed by directly solving the governing ODE which relates the orientation with angular velocity [3].

Apart from tracking the position and motion of every particle, this code also allows them to collide with each other. At every time-step binary collisions are tracked. To detect whether two ellipsoids are in contact a robust contact-detection algorithm has been developed, based on unconstrained optimization algorithm [4, 5] solved using Levenberg Marquardt Method [6] (The details regarding the algorithm will be updated shortly). The code also tracks collisions between the walls and the ellipsoids. Wall-Ellipsoid collision detection is solved using Lagrange Multiplier for a constrained optimization formulation.

Once it is detected that a pair of ellipsoids are in overlap a proper collision impulse is imparted on the ellipsoids [2].


## References

1. Goldstein, Herbert, Charles Poole, and John Safko. "Classical mechanics." (2002): 782-783.

2. Tschisgale, Silvio, Louis Thiry, and Jochen Fröhlich. "A constraint-based collision model for Cosserat rods." Archive of Applied Mechanics 89.2 (2019): 167-193.

3. Zhao, F., and B. G. M. Van Wachem. "A novel Quaternion integration approach for describing the behaviour of non-spherical particles." Acta Mechanica 224.12 (2013): 3091-3109.

4. Wellmann, Christian, Claudia Lillie, and Peter Wriggers. “A contact detection algorithm for superellipsoids based on the common-normal concept.” Engineering Computations (2008).

5. Pazouki, A., H. Mazhar, and D. Negrut. “Parallel Contact Detection between Ellipsoids with Applications in Granular Dynamics.” Mathematics and Computers in Simulation, Submit- ted (2010).

6. Minpack library for modern fortran available in the repository: https://github.com/fortran-lang/minpack. 
