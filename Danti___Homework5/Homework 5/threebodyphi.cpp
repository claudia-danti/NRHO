#include "threebodyphi.h"
#include <iostream>
#include "algebra.h"
#include "vector2D.h"
//#include "in-out.h"
using namespace std;



//uses the state vector of initial conditions to compute the equations of motion
vector<double> ThreeBodyPhi::evaluate(double t, const vector<double> & x){

	double r1 = sqrt(pow(x[0]+m_mu,2)+pow(x[1],2));
	double r2 = sqrt(pow(x[0]-1.+m_mu,2)+pow(x[1],2));

	vector<double> state(x.begin(), x.begin()+4);	//store the state vector
	vector2D<double> phi(4, vector<double>(x.begin()+4, x.end()));	//matrix 4x4

	//differential equations for sensitivity matrix
	double Uxx = 1 - (1-m_mu)/pow(r1,3) -m_mu/pow(r2,3) + 3*(1-m_mu)*pow(state[0]-x1, 2)/pow(r1,5) + 3*m_mu*pow(state[0]-x2, 2)/pow(r2,5);
	double Uxy = 3*state[1]*( (1-m_mu)*(state[0]-x1)/pow(r1,5) + m_mu*(state[0]-x2)/pow(r2,5) );
	double Uyy = 1 - (1-m_mu)/pow(r1,3) -m_mu/pow(r2,3) +3*pow(state[1],2)*( ((1-m_mu)/pow(r1,5)) + (m_mu/pow(r2,5)) );

	vector2D<double> A = {{0.,0.,1.,0.}, {0.,0.,0.,1.}, {Uxx, Uxy, 0.,2.}, {Uxy, Uyy, -2., 0.}};
	vector2D<double> phi_dot = A*phi;

	//differential equations for state vector
	vector<double> derivative(20);
	derivative[0] = state[2];
	derivative[1] = state[3];
	derivative[2] = state[0] - (1-m_mu)*(state[0] - x1)/pow(r1,3) - m_mu*(state[0]-x2)/pow(r2,3) + 2*state[3];
	derivative[3] = state[1] - (1-m_mu)*state[1]/pow(r1,3) - m_mu*state[1]/pow(r2,3) - 2*state[2];

	for(int i = 4; i < derivative.size(); i++){
		derivative[i] = phi_dot.element(i/4-1, 4%4);
	}

return derivative;

}