#include <iostream>
#include "algebra.h"
#include "stm.h"

using namespace std;



//uses the vector x t compute the differential equations
//N.B.: in the main x will be the state vector appended with the jacobi constant and with the phi matrix (so for a 2D problem is a vector with 4+1+16 rows)
vector<double> Stm::evaluate(double t, const vector<double> & x){
//create a generic definition for r1 and r2, regardless of the dimension of the problem
	double r1, r2;
	int dimension;
//this is the case of a 2D problem (dim_state_vect=4, phi=4x4)
//N.B.: the dimension is one more because there is the jacobi constant at the end of the state vector
	if(x.size()==20){

		r1 = sqrt(pow(x[1]-x1,2)+pow(x[2],2));
		r2 = sqrt(pow(x[1]-x2,2)+pow(x[2],2));
		dimension = 4;

	}
//this is the case of a 3D problem (dim_state_vect=6, phi=6x6)
//N.B.: the dimension is one more because there is the jacobi constant at the end of the state vector
	if(x.size()==42){

		r1 = sqrt(pow(x[0]-x1,2)+pow(x[1],2)+pow(x[2],2));
		r2 = sqrt(pow(x[0]-x2,2)+pow(x[1],2)+pow(x[2],2));
		dimension = 6;

	}

	vector<double> state(x.begin(), x.begin()+dimension);	//store the state vector
	vector<double> v_phi (x.begin()+dimension, x.end());	//store the stm into a vector v_phi
	vector2D<double> phi(dimension, vector<double>(x.begin()+dimension, x.end()));	//matrix 4x4

	//differential equations for sensitivity matrix
	double Uxx = 1 - (1-m_mu)/pow(r1,3) -m_mu/pow(r2,3) + 3*(1-m_mu)*pow(state[0]-x1, 2)/pow(r1,5) + 3*m_mu*pow(state[0]-x2, 2)/pow(r2,5);
	double Uxy = 3*state[1]*( (1-m_mu)*(state[0]-x1)/pow(r1,5) + m_mu*(state[0]-x2)/pow(r2,5) );
	double Uyy = 1 - (1-m_mu)/pow(r1,3) -m_mu/pow(r2,3) +3*pow(state[1],2)*( ((1-m_mu)/pow(r1,5)) + (m_mu/pow(r2,5)) );
	double Uxz = 3*((1-m_mu*(state[0]-x1)/pow(r1,5)) + m_mu*(state[0]-x2)/pow(r2,5) )*state[2];
	double Uyz = 3*((1-m_mu/pow(r1,5)) + m_mu/pow(r2,5) )*state[1]*state[2];
	double Uzz = - (1-m_mu)/pow(r1,3) - (m_mu)/pow(r2,3) + 3*((1-m_mu)/pow(r1,5) + m_mu/pow(r2,5))*pow(state[2],2);
	//define the A matrix that contains the differential equations for the stm
	vector2D<double> A;

	//again we make a distinction if the problem is 2D or 3D
	if (x.size() == 20){
		
		A = {{0.,0.,1.,0.}, {0.,0.,0.,1.}, {Uxx, Uxy, 0.,2.}, {Uxy, Uyy, -2., 0.}};
	}

	if(x.size() == 42){

		A = {{0.,0.,0.,1.,0.,0.}, {0.,0.,0.,0.,1.,0.},{0.,0.,0.,0.,0.,1.}, {Uxx, Uxy,Uxz, 0.,2.,0.}, {Uxy, Uyy, Uyz, -2., 0.,0.},{Uxz, Uyz, Uzz, 0.,0.,0.}};
	}

	vector2D<double> phi_dot = A*phi;

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
	}

	if(dimension == 6){

		derivative[3] = state[0] - (1-m_mu)*(state[0] - x1)/pow(r1,3) - m_mu*(state[0]-x2)/pow(r2,3) + 2*state[3];
		derivative[4] = state[1] - (1-m_mu)*state[1]/pow(r1,3) - m_mu*state[1]/pow(r2,3) - 2*state[2];
		derivative[5] = -(1-m_mu)*state[2]/pow(r1,3) - m_mu*state[2]/pow(r2,3);
	}

	//insert the differential equations for phi into the derivative vector
	for(int i = dimension; i < derivative.size(); i++){
		derivative[i] = phi_dot.element(i/dimension-1, dimension%dimension);
	}

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

	return (1-m_mu)/r1 +(m_mu)/r2 + (1./2)*(pow(x[0],2)+pow(x[1],2));
}


vector2D<double> jacobian2DHalo (vector2D<double> &phi, double vy, vector<double> xT){
//we write the jacobian using the stm: jacobian = phi (ereased) - 1/vy*(xT*vy)
	//store the vector of the y direction
	vector<double> phi_y = phi.getRow(1);
	//erease the y component of the vector
	phi_y.erase(phi_y.begin()+1);
	//erase the vx component of the vector
	phi_y.erase(phi_y.begin()+1);

	//erase rows and cols of the stm that correspond to y and vx dimensions
	//erase the y row
	phi.eraseRow(1);
	//erase the vx row (that after the first erase has become the first row)
	phi.eraseRow(1);

	//erase the y column
	phi.eraseCol(1);
	//erase the vx column
	phi.eraseCol(1);
	//stores the modified phi into a matrix a
	vector2D<double> a = phi;

//second memeber
	//computation of b/vy
	cout<<"prima di puter product"<<endl;
	vector2D<double> b = outerproduct(xT,phi_y);
	cout<<"outer product"<<endl;
	double cols = b.col();
	//we need to linearize the matrix in order to be able to divide every component of the matrix by vy
	vector<double> b_lin = b.getLinearVector();
	vector<double> b_lin_mod= b_lin/vy;
	//this is the new matrix
	vector2D<double> b_mod (cols, b_lin_mod);

	return a - b_mod;
}

vector2D<double> jacobian3DHalo (vector2D<double> &phi, double vy, vector<double> xT){
	//we write the jacobian using the stm: jacobian = phi (ereased) - 1/vy*(xT*vy)
	//store the vector of the y direction
	vector<double> phi_y = phi.getRow(1);
	//erease the y component of the vector
	phi_y.erase(phi_y.begin()+1);
	//erase the vx component of the vector
	phi_y.erase(phi_y.begin()+2);
	//erase the vz component of the vector
	phi_y.erase(phi_y.begin()+3);

	//erase rows and cols of the stm that correspond to y and vx dimensions
	//erase the y row
	phi.eraseRow(1);
	//erase the vx row (that after the first erase has become the first row)
	phi.eraseRow(2);
	//erease the vz row
	phi.eraseRow(3);

	//erase the y column
	phi.eraseCol(1);
	//erase the vx column
	phi.eraseCol(2);
	//erase the vz column
	phi.eraseCol(3);
	//stores the modified phi into a matrix a
	vector2D<double> a = phi;

	//computation of b/vy
	vector2D<double> b = outerproduct(xT,phi_y);
	double cols = b.col();
	//we need to linearize the matrix in order to be able to divide every component of the matrix by vy
	vector<double> b_lin = b.getLinearVector();
	vector<double> b_lin_mod= b_lin/vy;
	//this is the new matrix
	vector2D<double> b_mod (cols, b_lin_mod);

	return a - b_mod;

}


vector<double> Stm::dC_dx (const vector<double> & x){

	if(x.size() == 4){

		vector<double> dC;
		double r1 = sqrt(pow(x[1]-x1,2)+pow(x[2],2));
		double r2 = sqrt(pow(x[1]-x2,2)+pow(x[2],2));

		dC.push_back(2*(x[0]-(1-m_mu)*(x[0]-x1)/pow(r1,3) - m_mu/pow(r2,3)*(x[0]-x2)));
		dC.push_back(2*(x[1] - (1-m_mu)*x[1]/pow(r1,3) - m_mu*x[1]/pow(r2,3) - 2*x[2]));
		dC.push_back(-2*x[2]);
		dC.push_back(-2*x[3]);

		return dC;

	}

	else if(x.size() == 6){

		vector<double> dC;
		double r1 = sqrt(pow(x[0]-x1,2)+pow(x[1],2)+pow(x[2],2));
		double r2 = sqrt(pow(x[0]-x2,2)+pow(x[1],2)+pow(x[2],2));

		dC.push_back(2*(x[0]-(1-m_mu)*(x[0]-x1)/pow(r1,3) - m_mu/pow(r2,3)*(x[0]-x2)));
		dC.push_back(2*(x[1] - (1-m_mu)*x[1]/pow(r1,3) - m_mu*x[1]/pow(r2,3) - 2*x[2]));
		dC.push_back(2*(-(1-m_mu)*x[2]/pow(r1,3) - m_mu*x[2]/pow(r2,3)));
		dC.push_back(-2*x[3]);
		dC.push_back(-2*x[4]);
		dC.push_back(-2*x[5]);

		return dC;

	}

	else{
		cerr<<"The problem is neither 2 nor 3 dimensional."<<endl;

		vector<double> e (5);

		return e;
	}

}