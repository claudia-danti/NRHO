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


//Needed to see at which time step the body reaches the x-axis again with 0 velocity in x direction
class Check{

public:
	//constructor
	Check(double y, double vx): m_y_old(y), m_vx_old(vx){}

	//controls if the y variable in the state vector changes sign
	bool control(const vector<double> & y){

		//checks when y=0
		bool c = m_y_old*y[1] >= 0;
		m_y_old = y[1];

		//checks when vx=0
		bool d = m_vx_old*y[2] >= 0;
		m_vx_old = y[2];

		//checks if both conditions are satisfied
		bool a = c && d;

		return a;
	}

protected:
	double m_y_old;
	double m_vx_old;
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
//define the gravitational parameter using the units of the CR3BP
double mu = 0.01227741;
//this object contains the differental equations
ThreeBody tb (mu);
//matrix that contains the initial conditions (x,y,vx,vy)
vector<vector<double>> xin({{0.3,0.,0.,1.8},{0.5,0.,0.,1.2},{0.7,0.,0.,0.53},{0.8,0.,0.,0.35,3.0777}});
//matrix that contains the different initial conditions (x,y,vx,vy,C) that will lead to different orbits
vector<vector<double>> xinC({{0.3,0.,0.,1.8,2.968},{0.5,0.,0.,1.2,2.7066},{0.7,0.,0.,0.53,3.0572},{0.8,0.,0.,0.35,3.0777}});
vector<double> Cin = {2.968,2.7066,3.0572,3.0777};
//object Stm (will give us the differential equations for the whole system (state vector and state transition matrix))
Stm Phi (mu);
//defining the sensitivity matrix initial conditions (identity 4x4 matrix)
vector2D<double> phi0 = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};


//EXTREMES OF INTEGRATION
double t0 = 0.;
double tf = 20.;
double h = 0.01;
double H = 0.0001;

//AUXILIARY VARIABLES USED IN THE WHILE CYCLE
vector2D<double> res;	//matrix of the result of integration
vector2D<double> phi;	//stm
vector2D<double> res_jac;	//matrix fo rjacobi constant
vector<double> accT;	//acceleration vector
vector<double> K;	//constraint vector
vector<double> C;	//jacoi constant vector
vector<double> xT;	//define the state vector after one period
double C0;
double CT;
double j;
double k;
double T;
double d;
double treshold = pow(10,-6);	//treshold under which the diagonal values of the jacobian matrix are set to 0
double tolerance = pow(10,-11); //tolerance that regulates the iteration of the correction
								//the corrector algorithm goes on until the difference between final and initial condition is <tolerance

//ERRORS AND ADAPTIVE STEP FOR THE ODEINT INTEGRATOR
double rel_err = pow(10,-15);	//10 digits
double abs_err = pow(10,-10);
//initializing the integration stepper of odeint
auto thestepper = make_controlled<runge_kutta_fehlberg78<vector<double>>> (rel_err, abs_err);

int i = 1;


//this for cycle repeats the process for every initial condition vector proposed
for(auto x0 : xin){

	cout<<"CALCOLO DEL VETTORE "+to_string(i)<<"\t"<<flush;

	//initializing state to the vector of initial conditions
	vector<double> state = x0;
	//creating a matrix
	vector2D<double> initial;

	C0=tb.jacobi(x0);
	//open the output file and start the interation with the adaptive step
	output.open("corrections/initial"+to_string(i)+".dat");
	observer obs(output, initial);	//the results are written directly in the file
	integrate_adaptive( thestepper , tb , state , t0 , tf , h , obs );
	output.close();

	//the object check is the one that controls at which tikìme t the body reaches the x-axis again woth vx=0
	//it is mandatory to have it if we want to use the overloaded RK4 integrator that stops when the body reaches the x-axis
	Check c (x0[1], x0[2]);

	H = 1./10000;
	double v_old = 0.;
	bool flag = true;



	do{
		cout<<"*"<<flush;
		//initializing x0phi0 to the vector x0
		vector<double> x0phi0(x0);
		//add as vector the state transition matrix phi0 after the x0 vector, so that we get a unique vector that contains both the state vector and the stm
		x0phi0.insert(x0phi0.end(), phi0.getLinearVector().begin(), phi0.getLinearVector().end());

		//using the overloaded RK4 that stops the integration when we reach x-axis again, with vx=0
		double tol = pow(10,-20);
		res = RK4(t0, c, H, Phi, x0phi0);
		res.eraseRow(res.row()-1);	//we erase the last row of the res vector (it takes track of the time, so we get only the state vector plus phi)

		//the while has a conrol over the velocity vx, until it stays over the tolerance, the cycle continues
/*			//the integration step is reduced in order to get a better precision
			H = H / 10;
			//we store the initial time of integration (that is in the element res(last row, first column))
			double t_start = res.element(res.row()-1, 0);
			//we store the state vector at final t obtained by the rk integration into x0phi0
			vector<double> x0phi0 = res.getRow(res.row()-1);
			x0phi0.erase(x0phi0.begin());
			res = RK4(t_start, c, H, Phi, x0phi0);
			while( res.row() == 1 ){
				H = H / 10;
				res = RK4(t_start, c, H, Phi, x0phi0);
			}
			res.eraseRow(res.row()-1);
		}
*/
		///////////////////////PERIOD PARAMETERS/////////////////////////////
		j = res.row();	//reads how many iterations we have done to arrive at T
		T = res[j-1][0];	//stores the value of the period
		H = T / 10000;
		cout<<" T= "<<T<<endl;
		//store in xTphiT the values of the vector and stm after a period
		vector<double> xTphiT(res[j-1]);
		//store the value of the state vector at the period
		vector<double> xT(xTphiT.begin()+1, xTphiT.begin()+5);
		cout<<"state vector a metà periodo"<<xT<<endl;
		//store the value of the phi matrix at one period
		vector2D<double> phiT(4, vector<double>(xTphiT.begin()+5, xTphiT.end()));
		double vT = xT[3];	//velocity in y after one period
		CT = tb.jacobi(xT);

		flag = v_old*xT[2]>=0;
		v_old = xT[2];

		///////////////////////////////CORRECTION ROUTINE///////////////////////////////////////
		//constraint vector K
		K = xT-x0;
		//create the jacobian matrix
		vector2D<double> J = jacobian2DHalo(phiT,vT, xT);
		//the row of dC/dx must be added
		J.pushBackRow(Phi.dC_dx(x0));
		//now it must be converted in an armadillo matrix:
		//first we get a vector out of the matrix (so we use the armadillo constructor that takes a std:vector)
		vector<double> j = J.getLinearVector();
		//this is the jacobian matrix as an armadillo object
		mat jac ( & j.front(), J.row(), J.col());
		//now we use the svd on jac, and set to 0 all the diagonal element below the treshold
		//and get as a result the matrix to be applied to the constraint vector to get the correction
		mat svd_jac = pca(j, treshold);
		//convert constraint vector K into armadillo vector
		vec k (K);
		//apply the matrix to k to get the correction vector
		vec delta = - svd_jac*k;
		//convert the correction vector from armadillo vector to std: vector
		vector<double> Delta = conv_to<vector<double>>::from(delta);


		//finally we get the new initial conditions vector
		x0 = x0 + Delta;


		//tolerance to the constraint violation
		vector<double> pos0 = {x0[0],x0[1]};
		vector<double> posT = {xT[0],xT[1]};
		vector<double> vel0 = {x0[2],x0[3]};
		vector<double> velT = {xT[2],xT[3]};
		vector<double> deltapos = posT-pos0;
		vector<double> deltavel = velT-vel0;

		d = norm(deltapos)/norm(pos0) + norm(deltavel)/norm(vel0) + (CT-C0)/C0;



	}while(abs(d) > tolerance);	//the do-while stops when we have that the constraint violation is less than a tolerance

	cout<<endl;

	output.open("period_correction.dat", fstream::app);
	output << "VECTOR " << i << "\n";
	output << "period:\t" << setprecision(15)<<T << "\n";
	output << "corrected initial conditions:"<< setprecision(15)<<x0 << "\n\n";
	output.close();


	//print the corrected orbit into a file using the new orbit

	vector2D<double> correction;
	output.open("corrections/correction"+to_string(i)+".dat");
	observer obs2 (output, correction);
	integrate_adaptive( thestepper , tb , x0 , t0 , 6*T , h , obs2 );
	output.close();

	//state vector over one period
	res_jac = RK4(t0, T, T/100000, tb, x0);
	cout<<"state vector over one period: "<<"\n"<<res_jac<<endl;


	i++;

}


	return 0;
}