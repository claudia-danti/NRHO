#ifndef _stm_h_
#define _stm_h_

#include <vector>
#include <cmath>
#include"vector2D.h"


using namespace std;

//computation of the jacobian matrix for the Halo orbits correction routine, the input is the stm at a certain time t
vector2D<double> jacobian2DHalo (vector2D<double> &phi,double vy, vector<double> xT);
vector2D<double> jacobian3DHalo (vector2D<double> &phi, double vy, vector<double> xT);

class Stm {

public:

	//constructor
	Stm(double mu): m_mu(mu), x1(-mu), x2(1-mu){}

	//methods
	//evaluation of the differential equations for state vector and stm
	//the input is a vector x(state vector, C, stm) and therefore all the differential equations of the system are outputted
	vector<double> evaluate(double t, const vector<double> & x);
	//evaluation of the Jacobi constant and potential
	//these are needed to insert at the end of the equation for the state vector the evaluation of the Jacobi cosntant
	double jacobi(const vector<double> &x);
	double potential (const vector<double> & x);
	//computation of the derivative of the jacobian constant with respect to the initial conditions
	vector<double> dC_dx(const vector<double> &x);


protected:	
	double m_mu;
	double x1;
	double x2;
};


#endif