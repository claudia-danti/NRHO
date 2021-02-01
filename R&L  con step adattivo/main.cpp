#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <boost/numeric/odeint.hpp>
#include <armadillo>
#include "vector2D.h"
#include "in-out.h"
#include "stm.h"
#include "threeBody.h"
#include "RK4.h"
#include "pca.h"

using namespace std;
using namespace boost::numeric::odeint;
using namespace arma;

ofstream output; 


//Needed to see at which time step the body reaches the x-axis again in x direction
class Check{

public:
	//constructor
	Check(double y, double vy): m_y_old(y), m_vy_in(vy){}

	//controls if the y variable in the state vector changes sign
	bool control(const vector<double> & y){

		//checks when y=0
		bool y_check = m_y_old*y[1] >= 0;
		m_y_old = y[1];

		bool vy_check = y[3] * m_vy_in < 0;

		return y_check || vy_check;
	}

protected:
	double m_y_old;
	double m_vy_in;
};

////////////////////////////NEEDED FOR ODEINT INGEGRATION////////////////////////////////////

class observer{

public:

	//constructor (assigns to the os member a value (ostream &os))
	observer(ostream & os, vector2D<double> & result) : os(os), storage(result){};	
	
	template<typename State>		
	void operator()(const State & psi, double t){//overloading of operator ()
		vector<double> q = psi; //state of integration is stored in q
		q.insert(q.begin(),t);
		
		for(auto & i : q){
			os << "\t" << i;	//output of the components of the state vector separated by a tab
		}
		os << "\n";
		storage.pushBackRow(psi);
	}	

protected:
	ostream &os;
	vector2D<double> & storage;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////





int main(){


//OBJECTS AND INITIAL CONDITIONS

double mu = 0.01227741;	//define the gravitational parameter using the units of the CR3BP
ThreeBody tb (mu);	//this object contains the differental equations
Stm Phi (mu);	//object Stm (will give us the differential equations for the whole system (state vector and state transition matrix))

/*
vector<vector<double>> xin({{0.3,0.,0.,1.8},{0.5,0.,0.,1.2},{0.7,0.,0.,0.53},{0.8,0.,0.,0.35,3.0777}});	//matrix that contains the different initial conditions (x,y,vx,vy) that will lead to different orbits
vector<vector<double>> xinC({{0.3,0.,0.,1.8,2.968},{0.5,0.,0.,1.2,2.7066},{0.7,0.,0.,0.53,3.0572},{0.8,0.,0.,0.35,3.0777}});	//initial conditions with jacobi constant
vector<double> Cin = {2.968,2.7066,3.0572,3.0777};	//jacobi constants*/
vector2D<double> phi0 = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};	//defining the sensitivity matrix initial conditions (identity 4x4 matrix)
vector<double> x0 = {0.7,0.,0.,0.5};
vector2D<double> I = {{1,0},{0,1}};

vector2D<double> x;	//here the result of the integration is stored
double t0 = 0.;	//starting time of integration
double tf = 20.;
double h = 0.001;	//step of integration

//objects needed in the cycle
vector<double> xTphiT;	//stm after one period
vector<double> K;	//constraint vector
double T;	//period
double vyT;	//velocity in y after one period
double d;	//constraint violation
double threshold = pow(10,-4);	//threshold under which the diagonal values of the jacobian matrix are set to 0
double tolerance = pow(10,-11); //tolerance that regulates the iteration of the correction
								//the corrector algorithm goes on until the difference between final and initial condition is <tolerance

//ERRORS AND ADAPTIVE STEP FOR THE ODEINT INTEGRATOR
double rel_err = pow(10,-15);	//10 digits
double abs_err = pow(10,-20);
//initializing the integration stepper of odeint
auto thestepper = make_controlled<runge_kutta_fehlberg78<vector<double>>> (rel_err, abs_err);
Check c (x0[1], x0[3]);	//the object check is needed to check if the body reaches the x-axis with vx=0

/////////////INTEGRATION OF THE ORBIT WITH THE UNMODIFIED INITIAL CONDITIONS////////////
	vector<double> state = x0;
	//creating a matrix
	vector2D<double> initial;
	//open the output file and start the interation with the adaptive step
	output.open("orbits/initial.dat");
	observer obs(output, initial);	//the results are written directly in the file
	integrate_adaptive( thestepper , tb , state , t0 , tf , h , obs );
	output.close();

//////////////CYCLE THAT CORRECTS THE INITIAL CONDITIONS ////////////////////////
int cont = 0;
do{
		cout<<"*"<<cont<<endl;
		//we want to integrate the eqations of motion of the problem using one vector (for both x and phi)
		vector<double> x0phi0(x0);
		x0phi0.insert(x0phi0.end(), phi0.getLinearVector().begin(), phi0.getLinearVector().end());	//now x0phi0 contains all the initial conditions

		double C0 = tb.jacobi(x0);	//initial jacobi constant

		c = Check(x0[1], x0[3]);



		x = RK4(t0, c, h, Phi, x0phi0);

		output.open("orbits/rk4.dat");
		output<<x<<endl;
		output.close();


		cout<<"period: "<<x[x.row()-1][0]<<endl;

		////////////////////////CYCLE TO INTEGRATE WITH A SMALLER STEP AROUND THE Y=0 POINT/////////////
		//When we are near the end of the period, to achieve better precision, we reduce the integration 
		//step, in order to get better "after one period" values and therefore a more precise correction
		double H = h/100;	//step for the integration at the end of the period
		x.eraseRow(x.row()-1);	//this is the row where we already crossed the y axis again (so we must erase it)
		double tol = pow(10,-20);	//tolerance with which we accept the y value

		//the cycle keeps on going until y(T) > tol
		while( abs(x.element(x.row()-1, 2)) > tol ){

			H = H / 10;	//we reduce the integration step
			double t_start = x.element(x.row()-1, 0);	//store the initial time of integration
			vector<double> x0phi0_new = x.getRow(x.row()-1);	//store the vector of initial conditions
			x0phi0_new.erase(x0phi0_new.begin());	//we erase the time element (so that we have a vector (x,y,vx,vy))

			x = RK4(t_start, c, H, Phi, x0phi0_new);	//integrate with the new step and initial conditions
														//stops again when crossing the x-axis
			while( x.row() == 1 ){	//this cycle here prevents the fact that the rRK4 stops after one row
				H = H / 10;
				x = RK4(t_start, c, H, Phi, x0phi0_new);
			}

//			cout<<"y(T)= "<<x.element(x.row()-1,2)<<endl;
			x.eraseRow(x.row()-1);	//we erase the last row and restart the procedure to reduce the step H again
//			cout<<"Y(T)= "<<x.element(x.row()-1,2)<<endl;
		}


		//STM AND STATE VECTOR AFTER ONE PERIOD
		T = x[x.row()-1][0];	//period
		xTphiT = x[x.row()-1];	//state vector and stm after one period (contains also the time)
		vector<double> xT(xTphiT.begin()+1, xTphiT.begin()+5); //sv after one period (no time included)
		vector2D<double> phiT(4, vector<double>(xTphiT.begin()+5, xTphiT.end())); //state vector after one period
		vyT = xT[3];	//velocity in y direction after one period
		double CT = tb.jacobi(xT);

		///////////////////////////////CORRECTION ROUTINE///////////////////////////////////////

		K = xT-x0;	//definition of the constraint vector
		vector<double> xdotT = {xT[2],tb.acceleration(xT)[1]};	//vector that contains (vx,ay)//used to compute the jacobian matrix
		vector2D<double> J = jacobian2DHalo(phiT,vyT, xdotT)-I;	//jacobian matrix
		J.pushBackRow(Phi.dC_dx(x0));	//we add the row with the jacobi constant

		//CONVERSION IN ARMADILLO OBJECT
		//now it must be converted in an armadillo matrix:
		vector<double> j = J.getLinearVector();	//first we get a vector out of the matrix (so we use the armadillo constructor that takes a std:vector)
		mat jac ( & j.front(), J.row(), J.col());	//this is the jacobian matrix as an armadillo object

		//SINGULAR VALUE DECOMPOSITION
		//now we use the svd on jac, and set to 0 all the diagonal element below the treshold
		//and get as a result the matrix to be applied to the constraint vector to get the correction
		mat svd_jac = pca(jac, threshold);
		vec k (K);	//convert constraint vector K into armadillo vector
		k.shed_row(1);	//from the constraint vector K we need to delete 
		k.shed_row(1);	//the lines containing y and vx
		k.insert_rows(k.n_elem,1);	//insert in the constraint vector K the jacobi constant
		k(k.n_elem-1) = CT-C0;


		vec delta = -svd_jac*k;	//apply the matrix to k to get the correction vector
		vector<double> Delta = conv_to<vector<double>>::from(delta);	//convert the correction vector from armadillo vector to std: vector

		//applying the correction
		x0[0] += Delta[0];
		x0[3] += Delta[1];


		//defining the tolerance to the constraint violation
		vector<double> pos0 = {x0[0],x0[1]};
		vector<double> posT = {xT[0],xT[1]};
		vector<double> vel0 = {x0[2],x0[3]};
		vector<double> velT = {xT[2],xT[3]};

		cout<<"C0: "<<C0<<endl;
		cout<<"CT: "<<CT<<endl;
		d = norm(posT-pos0)/norm(pos0) + norm(velT-vel0)/norm(vel0) + abs((CT-C0)/C0);

		cout<<"d = "<<d<<endl;
		cout<<"Delta = "<<Delta<<endl;
		cout<<"XT= "<<xT<<endl;
		cout<<"x0= "<<x0<<endl;

		cont++;
}while(abs(d) > tolerance && cont < 50);	//the condition to stop the routine is to have a sufficiently small constrint vector K;
cout << "stop after " << cont << " iterations" << endl;

	//print the corrected vector for 6 periods in the file orbit.dat
	vector2D<double> orbit;
	output.open("orbits/orbit.dat");
	observer obs2 (output, orbit);
	integrate_adaptive( thestepper , tb , x0 , t0 , 9*T , h , obs2 );
	output.close();



	return 0;
}