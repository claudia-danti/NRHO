The main code is in the file main.cpp
There are several .h that contain different routines:

- stm.h/stm.cpp: is the library that computes the state transition matrix of a given 3BP with state vector x. In particular this routine gives the differential equations that the state vector and the STM have to satisfy. So, for example, for a 3D system, this library returns a vector "derivatives" that has, in order of rows, the equations of the state vector (with jacobi constant) and each component of the stm.

derivative[0]=dx/dt;
	.
	.
derivative[3]=dvx/dt
	.
derivative[5]=dvz/dt
derivative[6]=jacobi constant
derivative[7]=dphi/dt(element 1,1)
	.
	.
	.
derivative[43]=dphi/dt(element 6,6)
 
The output files are organized in folders:
- corrections: it contains all the file needed to plot the orbits (the initial conditions and the two corrected conditions) toghether with the appropriate python scripts to plot the orbits
- jacobi: the folder contains the files with the values of the Jacobi contant over one period and the related python script to plot it.

The files period_correction contain the values of the final initial vector and the period for the two approximation methods. (N.B.: every time the program runs they add lines)

In threeBody and threeBodyPhi are implemented the differential eqations of the problem.
