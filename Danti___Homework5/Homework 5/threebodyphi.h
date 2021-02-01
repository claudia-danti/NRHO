#ifndef _threebodyphi_h_
#define _threebodyphi_h_

#include <vector>
#include <cmath>
using namespace std;



class ThreeBodyPhi {

public:

	//constructor
	ThreeBodyPhi(double mu): m_mu(mu), x1(-mu), x2(1-mu){}

	//methods
	vector<double> evaluate(double t, const vector<double> & x);

protected:	
	double m_mu;
	double x1;
	double x2;
};


















#endif