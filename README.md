# The meridian CFD solver
## A 2D, unstructured, finite volume solver of the Euler flow equations
Meridian is a 2D finite volume CFD solver of the Euler equations,
intended to solve low-to-medium complexity problems with a basic
accuracy level.

This code features a series of classes providing:
* Unstructured mesh representation, generated from SU2 format
files (compatible with SU2 version 6.2)
* Space integration, using first order discretization
(constant reconstruction) and both centered / upwinding
schemes. Boundary conditions include solid walls, farfield
state, stagnation inflow and backpressure outflow.
* Time integration, using the explicit backward Euler scheme.
* Input / Output functionalities, allowing to dump the mesh
and the convergence history in tecplot format.

Usage:
* Compile the Main.cpp script filled-in with the use case 
settings (tested with C++ 11 / GCC 4.8.5).
* Launch the executable.

Upcoming improvements:
* Use case definition from input file.
* Second-order accurate spatial discretization.
* Multi-stage time-integration schemes.
