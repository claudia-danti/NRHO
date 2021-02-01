#ifndef _ThreeBody_h_
#define _ThreeBody_h_

#include <vector>
#include <cmath>
#include "vector2D.h"
using namespace std;


class ThreeBody{

public:

	//constructor
	ThreeBody(double mu): m_mu(mu), x1(-mu), x2(1-mu) {};

	//methods
	vector<double> evaluate(double t, const vector<double> & x);
	void operator()(const vector<double> & x, vector<double> &derivative, double t);
	vector<double> acceleration(const vector<double> & x);
	double potential (const vector<double> &x);


protected:
	double m_mu;
	double x1;
	double x2;
};

#endif