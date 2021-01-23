#include <iostream>
#include <armadillo>
#include "algebra.h"
#include "stm.h"

using namespace std;
using namespace arma;


//uses the vector x t compute the differential equations
//N.B.: in the main x will be the state vector appended with the jacobi constant and with the phi matrix (so for a 2D problem is a vector with 4+1+16 rows)
vector<double> Stm::evaluate(double t, const vector<double> & x){
//create a generic definition for r1 and r2, regardless of the dimension of the problem
	double r1, r2;
	int dimension;
//this is the case of a 2D problem (dim_state_vect=4, phi=4x4)
//N.B.: the dimension is one more because there is the jacobi constant at the end of the state vector
	if(x.size()==21){

		r1 = sqrt(pow(x[1]-x1,2)+pow(x[2],2));
		r2 = sqrt(pow(x[1]-x2,2)+pow(x[2],2));
		dimension = 4;

	}
//this is the case of a 3D problem (dim_state_vect=6, phi=6x6)
//N.B.: the dimension is one more because there is the jacobi constant at the end of the state vector
	if(x.size()==43){

		r1 = sqrt(pow(x[0]-x1,2)+pow(x[1],2)+pow(x[2],2));
		r2 = sqrt(pow(x[0]-x2,2)+pow(x[1],2)+pow(x[2],2));
		dimension = 6;

	}

	vector<double> state(x.begin(), x.begin()+dimension+1);	//store the state vector
	vector<double> v_phi (x.begin()+dimension+1, x.end());	//store the stm into a vector v_phi
	mat phi (const & v_phi, dimension, dimension)			//store the stm into a matrix phi

	//differential equations for sensitivity matrix
	double Uxx = 1 - (1-m_mu)/pow(r1,3) -m_mu/pow(r2,3) + 3*(1-m_mu)*pow(state[0]-x1, 2)/pow(r1,5) + 3*m_mu*pow(state[0]-x2, 2)/pow(r2,5);
	double Uxy = 3*state[1]*( (1-m_mu)*(state[0]-x1)/pow(r1,5) + m_mu*(state[0]-x2)/pow(r2,5) );
	double Uyy = 1 - (1-m_mu)/pow(r1,3) -m_mu/pow(r2,3) +3*pow(state[1],2)*( ((1-m_mu)/pow(r1,5)) + (m_mu/pow(r2,5)) );
	double Uxz = 3*((1-m_mu*(state[0]-x1)/pow(r1,5)) + m_mu*(state[0]-x2)/pow(r2,5) )*state[2];
	double Uyz = 3*((1-m_mu/pow(r1,5)) + m_mu/pow(r2,5) )*state[1]*state[2];
	double Uzz = - (1-m_mu)/pow(r1,3) - (m_mu)/pow(r2,3) + 3*((1-m_mu)/pow(r1,5) + m_mu/pow(r2,5))*pow(state[2],2);
	//define the A matrix that contains the differential equations for the stm
	mat A;

	//again we make a distinction if the problem is 2D or 3D
	if (x.size() == 21){
		
		A = {{0.,0.,1.,0.}, {0.,0.,0.,1.}, {Uxx, Uxy, 0.,2.}, {Uxy, Uyy, -2., 0.}};
	}

	if(x.size() == 43){

		A = {{0.,0.,0.,1.,0.,0.}, {0.,0.,0.,0.,1.,0.},{0.,0.,0.,0.,0.,1.}, {Uxx, Uxy,Uxz, 0.,2.,0.}, {Uxy, Uyy, Uyz, -2., 0.,0.},{Uxz, Uyz, Uzz, 0.,0.,0.}};
	}

	mat phi_dot = A*phi;

	//the return are the differential equations for state vector and matrix phi, in a vectorial form
	vector<double> derivative(x.size());

	//sets in the derivative vector the velocities
	for(int i = 0; i< dimension/2; i++){
		derivative[i] = state[i+2];
	}

	//again we have a difference for the 2D and 3D case
	if (dimension == 4){

		derivative[2] = state[0] - (1-m_mu)*(state[0] - x1)/pow(r1,3) - m_mu*(state[0]-x2)/pow(r2,3) + 2*state[3];
		derivative[3] = state[1] - (1-m_mu)*state[1]/pow(r1,3) - m_mu*state[1]/pow(r2,3) - 2*state[2];
		derivative[4] = jacobi(x);
	}

	if(dimension == 6){

		derivative[3] = state[0] - (1-m_mu)*(state[0] - x1)/pow(r1,3) - m_mu*(state[0]-x2)/pow(r2,3) + 2*state[3];
		derivative[4] = state[1] - (1-m_mu)*state[1]/pow(r1,3) - m_mu*state[1]/pow(r2,3) - 2*state[2];
		derivative[5] = -(1-m_mu)*state[2]/pow(r1,3) - m_mu*state[2]/pow(r2,3);
		derivative[6] = jacobi(x);
	}

	//vectorised form of the dphi/dt matrix
	vector<double> v_phi_dot = vectorise(phi_dot, dim=1);
	//add this vector at the end of derivative (that for now contained only the equations for the state vector)
	derivative.insert(derivative.end(), v_phi_dot.begin(), v_phi_dot.end());

return derivative;

}

//evaluation of the Jacobi constant
double Stm::jacobi(const vector<double> &x){

	int dimension;
	double C;
	double U = potential(x);

	if(x.size() == 21){
		C = 2*U - (pow(x[2],2)+pow(x[3],2));

	}

	if(x.size() == 43){
		C = 2*U - (pow(x[3],2)+pow(x[4],2)+pow(x[5],2));
	}

	return C;
}

//evaluation of the potential of the problem
double Stm::potential (const vector<double> &x){
//we differentiate the 2D from the 3D case
double r1, r2;

if(x.size() == 21){


	r1 = sqrt(pow(x[1]-x1,2)+pow(x[2],2));
	r2 = sqrt(pow(x[1]-x2,2)+pow(x[2],2));
}

if(x.size() == 43){

	r1 = sqrt(pow(x[0]-x1,2)+pow(x[1],2)+pow(x[2],2));
	r2 = sqrt(pow(x[0]-x2,2)+pow(x[1],2)+pow(x[2],2));

}

	return (1-m_mu)/r1 +(m_mu)/r2 + (1./2)*(pow(x[1],2)+pow(x[2],2));
}