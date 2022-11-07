## Abstract 

ELLBOX is a Lagrangian Particle Simulation code to track the movement of ellipsoidal particles in a cubical box. The ellipsoidal particles can collide with each other and the walls. The linear motion is computed using a simple forward-Euler scheme. The rotational motion is computed using the Euler equation of rotational motion using an RK4 scheme. The orientation of the ellipsoids are tracked using Quaternions. Apart from tracking the position and motion of every particle, this code also allows them to collide with each other. At every time-step binary collisions are tracked. To detect whether two ellipsoids are in contact a robust contact-detection algorithm has been developed, based on unconstrained optimization algorithm solved using Levenberg Marquardt Method (The details regarding the algorithm will be updated shortly). The code also tracks collisions between the walls and the ellipsoids. Wall-Ellipsoid collision detection is solved using Lagrange Multiplier for a constrained optimization formulation.

This code uses the Minpack library for modern fortran available under the repository: https://github.com/fortran-lang/minpack. 

