CODE EXPLANATION
The main code is in the file main.cpp
There are several .h that contain different routines:

- stm.h/stm.cpp: is the library that computes the state transition matrix of a given 3BP with state vector x. In particular this routine gives the differential equations of the system (the ones of the state vector) and the differential equations of the STM. So, for example, for a 3D system, this library returns a vector "derivatives" that has, in order of rows, the equations of the state vector and then the equations of each component of the stm.

derivative[0]=dx/dt
derivative[1]=dy/dt
	.
	.
	.
derivative[5]=dvz/dt
derivative[6]=dphi/dt(element 1,1)
	.
	.
	.
derivative[42]=dphi/dt(element 6,6)

- jacobian.h/jacobian.cpp: this library contains the routine to compute the jacobian of the problem, that is given by dK/dx, where by definition K is the vector of the difference between the final (after one period) and initial state vector. It can be computed using the following relation:
dK/dx=(dxT/dx0-I)
	(dC/dx0)
where I is the identity matrix and C the jacobi constant. Further manipulations of this expression lead to a calculation that involves the STM, that will be the one used to write this routine. In particular we address the case of initial conditions that lead to Halo orbits (this means with state vector x=(x0,0,z0,0,vy0,0)).
 
The output files are organized in folders:
- corrections: it contains all the file needed to plot the orbits (the initial conditions and the two corrected conditions) toghether with the appropriate python scripts to plot the orbits
- jacobi: the folder contains the files with the values of the Jacobi contant over one period and the related python script to plot it.

The files period_correction contain the values of the final initial vector and the period for the two approximation methods. (N.B.: every time the program runs they add lines)

In threeBody and threeBodyPhi are implemented the differential eqations of the problem.
