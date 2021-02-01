#ifndef _ThreeBody_h_
#define _ThreeBody_h_

#include <vector>
#include <cmath>

using namespace std;


class ThreeBody{

public:

	//constructor
	ThreeBody(double mu): m_mu(mu), x1(-mu), x2(1-mu) {};

	//methods
	//evaluate the differential equations
	vector<double> evaluate(double t, const vector<double> & x);
	//overload of the operator ()
	void operator()(const vector<double> & x, vector<double> & derivative, double t);
	//evaluate the acceleration
	vector<double> acceleration(const vector<double> & x);
	//evaluates the potential of the problem
	double potential (const vector<double> &x);
	//evaluates the Jacobi constant
	double jacobi(const vector<double> & x);


protected:
	double m_mu;
	double x1;
	double x2;
};

#endif