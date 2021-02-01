The main code is in the file main.cpp
The output files are organized in folders:
- corrections: it contains all the file needed to plot the orbits (the initial conditions and the two corrected conditions) toghether with the appropriate python scripts to plot the orbits
- jacobi: the folder contains the files with the values of the Jacobi contant over one period and the related python script to plot it.

The files period_correction contain the values of the final initial vector and the period for the two approximation methods. (N.B.: every time the program runs they add lines)

In threeBody and threeBodyPhi are implemented the differential eqations of the problem.