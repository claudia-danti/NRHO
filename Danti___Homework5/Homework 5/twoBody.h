#ifndef _twoBody_h_
#define _twoBody_h_
#include <vector>
#include <cmath>

using namespace std;
const double G = 6.67408*pow(10,-20);

class TwoBody{

public:

	//creators
	TwoBody();
	TwoBody(double m1, double m2);

	//destrutor
	~TwoBody();

	//set methods
	void set_M1(double m);
	void set_M2(double m);

	//get methods
	double get_M1();
	double get_M2();

	//evaluates the function in a point x at a time t
	vector<double> evaluate(double t, vector<double> x);


protected:
	double m_M1;
	double m_M2;

};

#endif
