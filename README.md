# Meridian CFD solver.
Meridian is a personal project aimed at developing a 2D unstructured
finite volume CFD solver of the Euler equations for solveing 
low-to-medium complexity problems with a basic accuracy level.

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

How to launch a simulation:
* Compile the Main.cpp script (tested with C++ 11 / GCC 4.8.5).
* Get your CFD mesh into SU2 format. Only quad and tri elements
are currently supported.
* Fill-in the configuration file with your settings.
* Launch the executable.

Upcoming improvements:
* Code architecture refactoring (source & header separation).
* Doxygen-compatible comments.
* Cmake-based compilation.
* Mesh metric pre-computation and storage.
* Gradient computation for second order accuracy schemes.
* Multi-stage time-integration schemes.
