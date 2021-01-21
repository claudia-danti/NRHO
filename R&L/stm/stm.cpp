#include "stm.h"
#include <iostream>
#include "algebra.h"
#include <armadillo>
#include "threeBody.h"
//#include "in-out.h"
using namespace std;
using namespace arma;



//takes as input the time t at which to compute the stm and the matrix phi(0)(matrix of initial conditions)
vector<double> Stm::evaluate(double t, const vecotr<double> & x){

	double r1 = sqrt(pow(x[0]+m_mu,2)+pow(x[1],2)+pow(x[2],2));
	double r2 = sqrt(pow(x[0]-1.+m_mu,2)+pow(x[1],2)+pow(x[2],2));

	vector<double> state(x.begin(), x.begin()+6);	//store the state vector (x,y,z,vx,vy,vz)
	vector2D phi(6, vector<double>(x.begin()+6, x.end()));	//matrix 6x6

	//differential equations for sensitivity matrix
	double Uxx = 1 - (1-m_mu)/pow(r1,3) -m_mu/pow(r2,3) + 3*(1-m_mu)*pow(state[0]-x1, 2)/pow(r1,5) + 3*m_mu*pow(state[0]-x2, 2)/pow(r2,5);
	double Uxy = 3*state[1]*( (1-m_mu)*(state[0]-x1)/pow(r1,5) + m_mu*(state[0]-x2)/pow(r2,5) );
	double Uyy = 1 - (1-m_mu)/pow(r1,3) -m_mu/pow(r2,3) +3*pow(state[1],2)*( ((1-m_mu)/pow(r1,5)) + (m_mu/pow(r2,5)) );
	double Uxz = 3*((1-m_mu*(state[0]-x1)/pow(r1,5)) + m_mu*(state[0]-x2)/pow(r2,5) )*state[2];
	double Uyz = 3*((1-m_mu/pow(r1,5)) + m_mu/pow(r2,5) )*state[1]*state[2];
	double Uzz = - (1-m_mu)/pow(r1,3) - (m_mu)/pow(r2,3) + 3*((1-m_mu)/pow(r1,5) + m_mu/pow(r2,5))*pow(state[2],2);

	vector2D A = {{0.,0.,0.,1.,0.,0.}, {0.,0.,0.,0.,1.,0.}, {Uxx, Uxy,Uxz, 0.,2.,0.}, {Uxy, Uyy, Uyz,-2., 0.,0.}, {Uxz, Uyz, Uzz,0.,0.,0.}};
	vector2D phi_dot = A*phi;


	//create a vector that contains all the differential equations that need to be solved (42 if the system is 3D)
	//differential equations for state vector
	vector<double> derivative(42);
	derivative[0] = state[3];
	derivative[1] = state[4];
	derivative[2] = state[5];
	derivative[3] = state[0] - (1-m_mu)*(state[0] - x1)/pow(r1,3) - m_mu*(state[0]-x2)/pow(r2,3) + 2*state[3];
	derivative[4] = state[1] - (1-m_mu)*state[1]/pow(r1,3) - m_mu*state[1]/pow(r2,3) - 2*state[2];
	derivative[5] = -(1-m_mu)*state[2]/pow(r1,3) - m_mu*state[2]/pow(r2,3);

	for(int i = 6; i < derivative.size(); i++){
		derivative[i] = phi_dot.element(i/6-1, 6%6);
	}

return derivative;

}