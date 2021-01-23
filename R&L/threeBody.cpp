#include <iostream>
#include "threeBody.h"
#include "algebra.h"
//#include "in-out.h"
using namespace std;


//uses the state vector of initial conditions to compute the equations of motion
vector<double> ThreeBody::evaluate(double t, const vector<double> & x){

	vector<double> derivative(x.size());

	int a = x.size()/2;

	vector<double> acc = acceleration(x);

	for(int i = 0; i < a; i++){
		derivative[i] = x[i+a];		//inserts the velocities in the first half of the vector derivative
		derivative[i+a] = acc[i];	//inserts accelerations in the last half of the vector
	}

	//add the jacobi constant to the state vector
	derivative.push_back(potential(x));

	return derivative;
}

//overload of the operator ()
void ThreeBody::operator()(const vector<double> & x, vector<double> & derivative, double t){

	int a = x.size()/2;

	vector<double> acc = acceleration(x);

	for(int i = 0; i < a; i++){
		derivative[i] = x[i+a];		//inserts the velocities in the first 2 positions of the vector derivative
		derivative[i+a] = acc[i];	//inserts the acceleration in the last positions
	}

}

//evaluation of the acceleration of the problem
vector<double> ThreeBody::acceleration(const vector<double> & x){

//create a vector for the accelerations that has the size of half of x 
//so that it's generic regardless of the dimension of the vector x
vector<double> acc(x.size()/2);
double r1, r2;

	if (x.size()/2 == 2){

	//position of the two primaries in the synodic frame
	r1 = sqrt(pow(x[0]-x1,2)+pow(x[1],2));
	r2 = sqrt(pow(x[0]-x2,2)+pow(x[1],2));

	acc[0] = x[0] - (1-m_mu)*(x[0] - x1)/pow(r1,3) - m_mu*(x[0]-x2)/pow(r2,3) + 2*x[3];
	acc[1] = x[1] - (1-m_mu)*x[1]/pow(r1,3) - m_mu*x[1]/pow(r2,3) - 2*x[2];
}

if (x.size()/2 ==3){

	//position of the two primaries in the synodic frame
	r1 = sqrt(pow(x[0]-x1,2)+pow(x[1],2)+pow(x[2],2));
	r2 = sqrt(pow(x[0]-x2,2)+pow(x[1],2)+pow(x[2],2));

	acc[0] = x[0] - (1-m_mu)*(x[0] - x1)/pow(r1,3) - m_mu*(x[0]-x2)/pow(r2,3) + 2*x[3];
	acc[1] = x[1] - (1-m_mu)*x[1]/pow(r1,3) - m_mu*x[1]/pow(r2,3) - 2*x[2];
	acc[2] = -(1-m_mu)*x[2]/pow(r1,3) - m_mu*x[2]/pow(r2,3);

}


	return acc;

}

//evaluation of the potential of the problem
double ThreeBody::potential (const vector<double> &x){
//we differentiate the 2D from the 3D case
int dimension = x.size()/2;
double r1, r2;

if(dimension==2){


	r1 = sqrt(pow(x[1]-x1,2)+pow(x[2],2));
	r2 = sqrt(pow(x[1]-x2,2)+pow(x[2],2));

}

if(dimension==3){

	r1 = sqrt(pow(x[0]-x1,2)+pow(x[1],2)+pow(x[2],2));
	r2 = sqrt(pow(x[0]-x2,2)+pow(x[1],2)+pow(x[2],2));

}

	return (1-m_mu)/r1 +(m_mu)/r2 + (1./2)*(pow(x[1],2)+pow(x[2],2));
}

//evaluation of the Jacobi constant
double ThreeBody::jacobi(const vector<double> &x){

	int dimension = x.size()/2;
	double C;
	double U = potential(x);

	if(dimension == 2){
		C = 2*U - (pow(x[2],2)+pow(x[3],2));

	}

	if(dimension == 3){
		C = 2*U - (pow(x[3],2)+pow(x[4],2)+pow(x[5],2));
	}

	return C;
}