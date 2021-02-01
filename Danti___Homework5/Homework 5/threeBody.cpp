#include "threeBody.h"
#include <iostream>
#include "algebra.h"
//#include "in-out.h"
using namespace std;


//uses the state vector of initial conditions to compute the equations of motion
vector<double> ThreeBody::evaluate(double t, const vector<double> & x){

	vector<double> derivative(x.size());

	int a = x.size()/2;

	vector<double> acc = acceleration(x);

	for(int i = 0; i < a; i++){
		derivative[i] = x[i+a];	//inserts the velocities in the first 2 positions of the vector derivative
		derivative[i+a] = acc[i];
	}
	return derivative;
}

void ThreeBody::operator()(const vector<double> & x, vector<double> &derivative, double t){

	int a = x.size()/2;

	vector<double> acc = acceleration(x);

	for(int i = 0; i < a; i++){
		derivative[i] = x[i+a];	//inserts the velocities in the first 2 positions of the vector derivative
		derivative[i+a] = acc[i];
	}

	derivative[0] = x[2];
	derivative[1] = x[3];
	derivative[2] = acc[0];
	derivative[3] = acc[1];

}


vector<double> ThreeBody::acceleration(const vector<double> & x){

	vector<double> acc(2);
	//position of the two massive bodies in the synodic frame

	double r1 = sqrt(pow(x[0]-x1,2)+pow(x[1],2));
	double r2 = sqrt(pow(x[0]-x2,2)+pow(x[1],2));

	acc[0] = x[0] - (1-m_mu)*(x[0] - x1)/pow(r1,3) - m_mu*(x[0]-x2)/pow(r2,3) + 2*x[3];
	acc[1] = x[1] - (1-m_mu)*x[1]/pow(r1,3) - m_mu*x[1]/pow(r2,3) - 2*x[2];

	return acc;
}


double ThreeBody::potential (const vector<double> &x){

	double r1 = sqrt(pow(x[0]-x1,2)+pow(x[1],2));
	double r2 = sqrt(pow(x[0]-x2,2)+pow(x[1],2));

	return (1-m_mu)/r1 +(m_mu)/r2 + (1./2)*(pow(x[0],2)+pow(x[1],2));
}