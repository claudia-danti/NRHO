#ifndef _RK4_h_
#define _RK4_h_

#include <cmath>
#include <vector>
#include <iostream>
#include "in-out.h"
#include "vector2D.h"
#include "algebra.h"	//this library contains the overloadings of the + - * / operators for vectro/matrixes
						//it allows to simply deal with operations between matrixes without having to do the cicles
						//in order to fill all the [i][j] elements of the matrix.

using namespace std;
//the object function is a vectorial function that must have the evaluate method (evaluating the function in a certain point x)
template <typename function>																			//parameters: initial time t0, final time tfin
vector2D<double> RK4(double t0, double tfin, double h, function & f, vector<double> y0){	//initial condition vector y0, function f, step h
	
	//number of steps of RK method (it's also the dimension of the matrix)
	int NrIterations = ((tfin-t0)/h);

	//definition of the 4 increasing vectors k
	vector<double> k1(y0.size());
	vector<double> k2(y0.size());
	vector<double> k3(y0.size());
	vector<double> k4(y0.size());

	//vector that contains all the times 
	vector<double> t (NrIterations);
	for(int i = 0; i < NrIterations; i++){
		t[i] = t0 + i*h;
	}
	
	//matrix of the solution
	//NrIteration vectors that contain the values of the solutions at a time t
	/* considering n=NrIteration time steps (so that t=[t0;tn]) and a problem with dimension d
	(so that the vector y=[y0;yd]), the matrix can be imagined as y[i][j] (i line index, j column index)

	y0(t0), y1(t0), ..., yd(t0)
	y0(t1), y1(t1), ..., yd(t1)
	...,				, ...
	...,				, ...	
	...,				, ...
	y0(tn), y1(tn), ..., yd(tn)

	*/
	vector2D<double> y (NrIterations,y0.size());
	//setting in the first raw the initial conditions
	y.setRow(0, y0);

	//RK method
		for(int i = 1; i < NrIterations; i++){

			k1 = h*f.evaluate(t[i-1],y[i-1]);							//the values of the column vector k are always
			k2 = h*f.evaluate(t[i-1] + 0.5*h, y[i-1] + 0.5*k1);		//overwritten, since we have no interest in keeping
			k3 = h*f.evaluate(t[i-1] + 0.5*h, y[i-1] + 0.5*k2);		//track of them (so there is no point in having a matrix)
			k4 = h*f.evaluate(t[i-1] + h, y[i-1] + k3);

		//thanks to the overloading of the operators (algebra.h library) and to the vector2d library we can write in a simplier way the
		//values of the matrix y (in particular the [j] component of y authomatically corresponds to the one of k)
			y.setRow( i, y[i-1] + (1./6.)*(k1 + 2.*k2 + 2.*k3 + k4) );
			cout << i*100/NrIterations << "%\r" << flush;	//prints the percentual of the computation
		}


	//insert the vector t at the beginning of the matrix, so that the uotput of the function is a matrix
	//that contains the vector of times and the matrix of the solution of the ODE
	//the matrix contains the vector of the solution evaluated for times times t in [t0,tfin]
	//y = (t0,...,tfin; y1(t0),...,y1(tfin);...;yp(t0),...yp(tfin))
	//where p is the dimesion of the vector y0.

		y.insertCol(0, t);


	return y;

	}


//OVERLOAD with control
//the object function is a vectorial function that must have the evaluate method (evaluating the function in a certain point x)
template <typename function, typename control>																			//parameters: initial time t0, final time tfin
vector2D<double> RK4(double t0, control &c, double h, function & f, vector<double> y0){	//initial condition vector y0, function f, step h

	//definition of the 4 increasing vectors k
	vector<double> k1(y0.size());
	vector<double> k2(y0.size());
	vector<double> k3(y0.size());
	vector<double> k4(y0.size());

		
	//matrix of the solution
	//NrIteration vectors that contain the values of the solutions at a time t
	/* considering n=NrIteration time steps (so that t=[t0;tn]) and a problem with dimension d
	(so that the vector y=[y0;yd]), the matrix can be imagined as y[i][j] (i line index, j column index)

	y0(t0), y1(t0), ..., yd(t0)
	y0(t1), y1(t1), ..., yd(t1)
	...,				, ...
	...,				, ...	
	...,				, ...
	y0(tn), y1(tn), ..., yd(tn)

	*/
//output.open("debug.dat");
	//vector that contains all the times 
	vector<double> t;
	t.push_back(t0);

	vector2D<double> y;
	//setting in the first raw the initial conditions
	y.pushBackRow(y0);


	int j = 0;

		while(c.control(y[j])){

			j++;

			k1 = h*f.evaluate(t[j-1],y[j-1]);							//the values of the column vector k are always
			k2 = h*f.evaluate(t[j-1] + 0.5*h, y[j-1] + 0.5*k1);		//overwritten, since we have no interest in keeping
			k3 = h*f.evaluate(t[j-1] + 0.5*h, y[j-1] + 0.5*k2);		//track of them (so there is no point in having a matrix)
			k4 = h*f.evaluate(t[j-1] + h, y[j-1] + k3);

		//thanks to the overloading of the operators (algebra.h library) and to the vector2d library we can write in a simplier way the
		//values of the matrix y (in particular the [j] component of y authomatically corresponds to the one of k)
			y.pushBackRow(y[j-1] + (1./6.)*(k1 + 2.*k2 + 2.*k3 + k4) );
			t.push_back(t0 + j*h);

			//output<<y[j]<<endl;
				
		}

	//insert the vector t at the beginning of the matrix, so that the uotput of the function is a matrix
	//that contains the vector of times and the matrix of the solution of the ODE
	//the matrix contains the vector of the solution evaluated for times times t in [t0,tfin]
	//y = (t0,...,tfin; y1(t0),...,y1(tfin);...;yp(t0),...yp(tfin))
	//where p is the dimesion of the vector y0.

		y.insertCol(0, t);

	return y;

	}

#endif
