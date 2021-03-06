#include "threeBody.h"
#include "RK4.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include "vector2D.h"
#include "in-out.h"
#include "threebodyphi.h"
#include <iomanip>
#include <string>
#include <boost/numeric/odeint.hpp>


using namespace std;
using namespace boost::numeric::odeint;

ofstream output; 

class Check{

public:
	//constructor
	Check(double y): m_y_old(y){}

	//controls if the y variable in the state vector changes sign
	bool control(const vector<double> & y){

		bool c = m_y_old*y[1] >= 0;
		m_y_old = y[1];

		return c;
	}

protected:
	double m_y_old;
};

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

int main(){

//define the gravitational parameter using the units of the CR3BP
double mu = 0.01227741;

vector<vector<double>> xin({{0.3,0.,0.,1.8},{0.5,0.,0.,1.2},{0.7,0.,0.,0.53},{0.8,0.,0.,0.35}});

ThreeBody tb (mu);

//EXTREMES OF INTEGRATION
double t0 = 0.;
double tf = 20.;
double h = 0.01;
double H = 0.0001;

//defining all the variable that will be used in the while
vector2D<double> res;
double vxT2;
vector<double> accT2;
vector2D<double> phi;


double j;
double k;
double T2;

double Bxdot;
double Axdot;
double Ax;
double Bx;
double Aydot;
double Bydot;
double dvy0;
vector<double> corr;
vector2D<double> res_jac;

vector<double> C;

//defining the sensitivity matrix
vector2D<double> phi0 = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
ThreeBodyPhi Phi (mu);

//defining absolute and relative erro we want
double rel_err = pow(10,-15);	//10 digits
double abs_err = pow(10,-10);
//qua stanziamo lo stepper
auto thestepper = make_controlled<runge_kutta_fehlberg78<vector<double>>> (rel_err, abs_err);

int i = 1;

///////////////////STRATEGY 1: VINCOLATED INITIAL POSITION/////////////////////

for(auto x0 : xin){

	cout<<"CALCOLO DEL VETTORE "+to_string(i)<<"\t"<<flush;

	vector<double> state = x0;
	vector2D<double> initial;
	output.open("corrections/initial"+to_string(i)+".dat");
	observer obs(output, initial);
	integrate_adaptive( thestepper , tb , state , t0 , tf , h , obs );
	output.close();

	Check c (x0[1]);

	H = 1./10000;
	double v_old = 0.;
	bool flag = true;

	do{
		cout<<"*"<<flush;
		vector<double> x0phi0(x0);
		x0phi0.insert(x0phi0.end(), phi0.getLinearVector().begin(), phi0.getLinearVector().end());

		//using the overloaded RK4 that stops the integration when we reach x-axis again
		double tolerance = pow(10,-20);
		res = RK4(t0, c, H, Phi, x0phi0);
		
		res.eraseRow(res.row()-1);

		while( res.element(res.row()-1, 2) > tolerance ){
			H = H / 10;
			double t_start = res.element(res.row()-1, 0);
			vector<double> x0phi0 = res.getRow(res.row()-1);
			x0phi0.erase(x0phi0.begin());
			res = RK4(t_start, c, H, Phi, x0phi0);
			while( res.row() == 1 ){
				H = H / 10;
				res = RK4(t_start, c, H, Phi, x0phi0);
			}
			res.eraseRow(res.row()-1);
		}
		///////////////////////HALF OF THE PERIOD PARAMETERS/////////////////////////////
		j = res.row();	//reads how many iterations we have done to arrive at T/2
		T2 = res[j-1][0];	//stores the value of half of the period
		H = T2 / 10000;
//		cout<<" T/2= "<<T2<<endl;

		vector<double> xT2phiT2(res[j-1]);

		//store the value of the state vector at half of the period
		vector<double> xT2(xT2phiT2.begin()+1, xT2phiT2.begin()+5);
//		cout<<"state vector a metà periodo"<<xT2<<endl;
		vxT2 = xT2[2];

		//store the value of the phi matrix at half of the period
		vector2D<double> phiT2(4, vector<double>(xT2phiT2.begin()+5, xT2phiT2.end()));

		flag = v_old*xT2[2]>=0;
		v_old = xT2[2];

		//gets the value of the acceleration at half of the period
		accT2 = tb.acceleration(xT2);

		//defining the parameters
		Bxdot = phiT2.element(2, 3) * xT2[3] - phiT2.element(1, 3) * accT2[0];

		///////////////////STRATEGY 1: VINCOLATED INITIAL POSITION/////////////////////

		corr = {0., 0., 0., -xT2[2]*xT2[3]/Bxdot};

//		cout<<"correzione sulla velocità iniziale "<<corr[3]<<endl;

		//correct the initial conditions
		x0 = x0 + corr;

	}while(abs(corr[3])>pow(10,-15));

	cout<<endl;

	output.open("period_correction.dat", fstream::app);
	output << "VECTOR " << i << "\n";
	output << "period:\t" << setprecision(15)<<2*T2 << "\n";
	output << "corrected initial conditions:"<< setprecision(15)<<x0 << "\n\n";
	output.close();


	//print the corrected orbit into a file using the new orbit

	vector2D<double> correction;
	output.open("corrections/correction"+to_string(i)+".dat");
	observer obs2 (output, correction);
	integrate_adaptive( thestepper , tb , x0 , t0 , 6*T2 , h , obs2 );
	output.close();

	//state vector over one period
	res_jac = RK4(t0, 2*T2, T2/100000, tb, x0);

	output.open("jacobi/jacobi1"+to_string(i)+".dat");

	//vector of the angular momentum (has only the z component)
	for(int k = 0; k < res_jac.row(); k++ ){
		output << res_jac.element(k, 0) << "\t" <<setprecision(15)<< 2*tb.potential(res_jac[k]) -
		(pow(res_jac.element(k, 3), 2) + pow(res_jac.element(k, 4), 2)) << "\n";
	}


	output.close();

	i++;

}

///////////////////STRATEGY 2: VINCOLATED INITIAL VELOCITY/////////////////////
i = 1;

for(auto x0 : xin){

	cout<<"CALCOLO DEL VETTORE "+to_string(i)<<"\t"<<flush;

	vector<double> state = x0;
	vector2D<double> initial;
	output.open("corrections/initial2"+to_string(i)+".dat");
	observer obs(output, initial);
	integrate_adaptive( thestepper , tb , state , t0 , tf , h , obs );
	output.close();

	Check c (x0[1]);

	T2 = 1;
	double v_old = 0.;
	bool flag = true;

	do{
		cout<<"*"<<flush;
		vector<double> x0phi0(x0);
		x0phi0.insert(x0phi0.end(), phi0.getLinearVector().begin(), phi0.getLinearVector().end());

		//using the overloaded RK4 that stops the integration when we reach x-axis again
		double tolerance = pow(10,-20);
		res = RK4(t0, c, H, Phi, x0phi0);
		res.eraseRow(res.row()-1);
		while( res.element(res.row()-1, 2) > tolerance ){
			H = H / 10;
			double t_start = res.element(res.row()-1, 0);
			vector<double> x0phi0 = res.getRow(res.row()-1);
			x0phi0.erase(x0phi0.begin());
			res = RK4(t_start, c, H, Phi, x0phi0);
			while( res.row() == 1 ){
				H = H / 10;
				res = RK4(t_start, c, H, Phi, x0phi0);
			}
			res.eraseRow(res.row()-1);
		}

		///////////////////////HALF OF THE PERIOD PARAMETERS/////////////////////////////
		j = res.row();	//reads how many iterations we have done to arrive at T/2
		T2 = res[j-1][0];	//stores the value of half of the period
		H = T2 / 10000;
//		cout<<" T/2= "<<T2<<endl;

		vector<double> xT2phiT2(res[j-1]);

		//store the value of the state vector at half of the period
		vector<double> xT2(xT2phiT2.begin()+1, xT2phiT2.begin()+5);
//		cout<<"state vector a metà periodo"<<xT2<<endl;
		vxT2 = xT2[2];

		//store the value of the phi matrix at half of the period
		vector2D<double> phiT2(4, vector<double>(xT2phiT2.begin()+5, xT2phiT2.end()));

		flag = v_old*xT2[2]>=0;
		v_old = xT2[2];

		//gets the value of the acceleration at half of the period
		accT2 = tb.acceleration(xT2);

		//defining the parameters
		Axdot = phiT2.element(2, 0) * xT2[3] - phiT2.element(1, 0) * accT2[0];

		corr = {-xT2[2]*xT2[3]/Axdot, 0., 0., 0.};

//		cout<<"correzione sulla velocità iniziale "<<corr[0]<<endl;

		//correct the initial conditions
		x0 = x0 + corr;

	}while(abs(corr[0])>pow(10,-15));

	cout<<endl;

	output.open("period_correction2.dat", fstream::app);
	output << "VECTOR " << i << "\n";
	output << "period:\t" << setprecision(15)<<2*T2 << "\n";
	output << "corrected initial conditions:" << setprecision(15)<< x0 << "\n\n";
	output.close();


	//print the corrected orbit into a file using the new orbit

	vector2D<double> correction;
	output.open("corrections/correction2"+to_string(i)+".dat");
	observer obs2 (output, correction);
	integrate_adaptive( thestepper , tb , x0 , t0 , 6*T2 , h , obs2 );
	output.close();


	//state vector over one period
	res_jac = RK4(t0, 2*T2, T2/100000, tb, x0);

	output.open("jacobi/jacobi2"+to_string(i)+".dat");

	//vector of the angular momentum (has only the z component)
	for(int k = 0; k < res_jac.row(); k++ ){
		output << res_jac.element(k, 0) << "\t" <<setprecision(15)<< 2*tb.potential(res_jac[k]) -
		(pow(res_jac.element(k, 3), 2) + pow(res_jac.element(k, 4), 2)) << "\n";
	}


	output.close();


	i++;

}


/*
///////////////////STRATEGY 3: VINCOLATED FINAL POSITION/////////////////////
i = 1;

for(auto x0 : xin){

	cout<<"CALCOLO DEL VETTORE "+to_string(i)<<"\t"<<flush;

	vector<double> state = x0;
	vector2D<double> initial;
	output.open("initial3"+to_string(i)+".dat");
	observer obs(output, initial);
	integrate_adaptive( thestepper , tb , state , t0 , tf , h , obs );
	output.close();

	Check c (x0[1]);

	T2 = 1;
	double v_old = 0.;
	bool flag = true;

	do{
		cout<<"*"<<flush;
		vector<double> x0phi0(x0);
		x0phi0.insert(x0phi0.end(), phi0.getLinearVector().begin(), phi0.getLinearVector().end());

		//using the overloaded RK4 that stops the integration when we reach x-axis again
		double tolerance = pow(10,-20);
		res = RK4(t0, c, H, Phi, x0phi0);
		res.eraseRow(res.row()-1);
		while( res.element(res.row()-1, 2) > tolerance ){
			H = H / 10;
			double t_start = res.element(res.row()-1, 0);
			vector<double> x0phi0 = res.getRow(res.row()-1);
			x0phi0.erase(x0phi0.begin());
			res = RK4(t_start, c, H, Phi, x0phi0);
			while( res.row() == 1 ){
				H = H / 10;
				res = RK4(t_start, c, H, Phi, x0phi0);
			}
			res.eraseRow(res.row()-1);
		}

		///////////////////////HALF OF THE PERIOD PARAMETERS/////////////////////////////
		j = res.row();	//reads how many iterations we have done to arrive at T/2
		T2 = res[j-1][0];	//stores the value of half of the period
		H = T2 / 10000;

		vector<double> xT2phiT2(res[j-1]);

		//store the value of the state vector at half of the period
		vector<double> xT2(xT2phiT2.begin()+1, xT2phiT2.begin()+5);
//		cout<<"state vector a metà periodo"<<xT2<<endl;
		vxT2 = xT2[2];

		//store the value of the phi matrix at half of the period
		vector2D<double> phiT2(4, vector<double>(xT2phiT2.begin()+5, xT2phiT2.end()));

		flag = v_old*xT2[2]>=0;
		v_old = xT2[2];

		//gets the value of the acceleration at half of the period
		accT2 = tb.acceleration(xT2);

		//defining the parameters
		Axdot = phiT2.element(2, 0) * xT2[3] - phiT2.element(1, 0) * accT2[0];
		Bxdot = phiT2.element(2, 3) * xT2[3] - phiT2.element(1, 3) * accT2[0];
		Bx = phiT2.element(0,3) * xT2[3] - phiT2.element(1,3) * xT2[2];
		Ax = phiT2.element(0,0) * xT2[3] - phiT2.element(1,0) * xT2[2];

		corr = {xT2[2]*xT2[3]*Bx/(Ax*Bxdot - Axdot*Bx), 0., 0., xT2[2]*xT2[3]*Ax/(Axdot*Bx - Ax*Bxdot)};

		//correct the initial conditions
		x0 = x0 + corr;

	}while(abs(corr[0])>pow(10,-8) && abs(corr[3])>pow(10,-8));

	cout<<endl;

	output.open("period_correction3.dat", fstream::app);
	output << "VECTOR " << i << "\n";
	output << "period:\t" << 2*T2 << "\n";
	output << "corrected initial conditions:" << x0 << "\n\n";
	output.close();


	//print the corrected orbit into a file using the new orbit

	vector2D<double> correction;
	output.open("correction3"+to_string(i)+".dat");
	observer obs2 (output, correction);
	integrate_adaptive( thestepper , tb , x0 , t0 , 8*T2 , h , obs2 );
	output.close();

/*	output.open("period_initial_SV"+to_string(i)+".dat");
	output<<"Period of the orbit"<<endl;
	output<<setprecision(13)<<T2*2<<endl;
	output<<"Final initial state vector"<<endl;
	for(int a = 1; a<5; a++){
		output<<setprecision(13)<<res[0][a]<<"\t";
	}
	output.close();

	//computing the value of the Jacobi constant over one period
	output.open("onePeriod"+to_string(i)+".dat");
	//state vector over one period
	res_jac = RK4(t0, 2*T2, h, tb, x0);
	output<<res_jac<<endl;
	output.close();

	double r = res_jac.row();

	output.open("jacobi1"+to_string(i)+".dat");

	//vector of the angular momentum (has only the z component)
	for(int a = 0; a < r; a++ ){

		C.push_back(2*(tb.potential(res_jac[a])) - pow(res_jac[a][3],2) - pow(res_jac[a][4],2));
		output<<res_jac[a][0]<<"\t"<<C[a]<<endl;
	}


	output.close();

	i++;

}


///////////////////STRATEGY 4: VINCOLATED FINAL VELOCITY/////////////////////
i = 1;

for(auto x0 : xin){

	cout<<"CALCOLO DEL VETTORE "+to_string(i)<<"\t"<<flush;

	vector<double> state = x0;
	vector2D<double> initial;
	output.open("initial4"+to_string(i)+".dat");
	observer obs(output, initial);
	integrate_adaptive( thestepper , tb , state , t0 , tf , h , obs );
	output.close();

	Check c (x0[1]);

	T2 = 1;
	double v_old = 0.;
	bool flag = true;

	do{
		cout<<"*"<<flush;
		vector<double> x0phi0(x0);
		x0phi0.insert(x0phi0.end(), phi0.getLinearVector().begin(), phi0.getLinearVector().end());

		//using the overloaded RK4 that stops the integration when we reach x-axis again
		double tolerance = pow(10,-15);
		res = RK4(t0, c, H, Phi, x0phi0);
		res.eraseRow(res.row()-1);
		while( res.element(res.row()-1, 2) > tolerance ){
			H = H / 10;
			double t_start = res.element(res.row()-1, 0);
			vector<double> x0phi0 = res.getRow(res.row()-1);
			x0phi0.erase(x0phi0.begin());
			res = RK4(t_start, c, H, Phi, x0phi0);
			while( res.row() == 1 ){
				H = H / 10;
				res = RK4(t_start, c, H, Phi, x0phi0);
			}
			res.eraseRow(res.row()-1);
		}

		///////////////////////HALF OF THE PERIOD PARAMETERS/////////////////////////////
		j = res.row();	//reads how many iterations we have done to arrive at T/2
		T2 = res[j-1][0];	//stores the value of half of the period
		H = T2 / 10000;

		vector<double> xT2phiT2(res[j-1]);

		//store the value of the state vector at half of the period
		vector<double> xT2(xT2phiT2.begin()+1, xT2phiT2.begin()+5);
//		cout<<"state vector a metà periodo"<<xT2<<endl;
		vxT2 = xT2[2];

		//store the value of the phi matrix at half of the period
		vector2D<double> phiT2(4, vector<double>(xT2phiT2.begin()+5, xT2phiT2.end()));

		flag = v_old*xT2[2]>=0;
		v_old = xT2[2];

		//gets the value of the acceleration at half of the period
		accT2 = tb.acceleration(xT2);

		//defining the parameters
		Axdot = phiT2.element(2, 0) * xT2[3] - phiT2.element(1, 0) * accT2[0];
		Bxdot = phiT2.element(2, 3) * xT2[3] - phiT2.element(1, 3) * accT2[0];
		Aydot = phiT2.element(3,0) * xT2[3] - phiT2.element(1,0) * accT2[1];
		Bydot = phiT2.element(3,3) * xT2[3] - phiT2.element(1,3) * accT2[1];

		corr = {xT2[2]*xT2[3]*Bydot/(Aydot*Bxdot - Axdot*Bydot), 0., 0., xT2[2]*xT2[3]*Aydot/(Axdot*Bydot - Aydot*Bxdot)};

		//correct the initial conditions
		x0 = x0 + corr;

	}while(abs(corr[0])>pow(10,-10) && abs(corr[3])>pow(10,-10));

	cout<<endl;

	output.open("period_correction4.dat", fstream::app);
	output << "VECTOR " << i << "\n";
	output << "period:\t" << 2*T2 << "\n";
	output << "corrected initial conditions:" << x0 << "\n\n";
	output.close();


	//print the corrected orbit into a file using the new orbit

	vector2D<double> correction;
	output.open("correction4"+to_string(i)+".dat");
	observer obs2 (output, correction);
	integrate_adaptive( thestepper , tb , x0 , t0 , 8*T2 , h , obs2 );
	output.close();

/*	output.open("period_initial_SV"+to_string(i)+".dat");
	output<<"Period of the orbit"<<endl;
	output<<setprecision(13)<<T2*2<<endl;
	output<<"Final initial state vector"<<endl;
	for(int a = 1; a<5; a++){
		output<<setprecision(13)<<res[0][a]<<"\t";
	}
	output.close();

	//computing the value of the Jacobi constant over one period
	output.open("onePeriod"+to_string(i)+".dat");
	//state vector over one period
	res_jac = RK4(t0, 2*T2, h, tb, x0);
	output<<res_jac<<endl;
	output.close();

	double r = res_jac.row();

	output.open("jacobi1"+to_string(i)+".dat");

	//vector of the angular momentum (has only the z component)
	for(int a = 0; a < r; a++ ){

		C.push_back(2*(tb.potential(res_jac[a])) - pow(res_jac[a][3],2) - pow(res_jac[a][4],2));
		output<<res_jac[a][0]<<"\t"<<C[a]<<endl;
	}


	output.close();

	i++;

}
*/

	return 0;
}