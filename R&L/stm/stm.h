#ifndef _stm_h_
#define _stm_h_

#include <vector>
#include <cmath>
#include <armadillo>
#include "threeBody.h"

using namespace std;
using namespace arma;


class Stm {

public:

	//constructor
	Stm(double mu): m_mu(mu), x1(-mu), x2(1-mu){}

	//methods
	//evaluates the differential equations of the stm dphi/dt=A*phi
	vector<double> evaluate(double t, const vector<double> & x);

protected:	
	double m_mu;
	double x1;
	double x2;
};



















#endif