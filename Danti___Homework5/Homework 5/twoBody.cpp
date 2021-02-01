#include "twoBody.h"
#include <cmath>
#include <iostream>
#include "algebra.h"
using namespace std;

//creator
TwoBody::TwoBody(double m1, double m2){
	m_M1 = m1;
	m_M2 = m2;
}

//distructor
TwoBody::~TwoBody(){};

//set methods
void TwoBody::set_M1(double m){
	m_M1 = m;
}
void TwoBody::set_M2(double m){
	m_M2 = m;
}

//get methods
double TwoBody::get_M1(){
	return m_M1;
}

double TwoBody::get_M2(){
	return m_M2;
}

//evaluates the derivative of positons and velocities of both bodies
//the vector x is a column vector that contains tha spatial coordinates of thr first body,
//the spatial coordinates of the second body and velocities of the two bodies:
//v=(x1,y1,...,x2,y2,...,vx1,vy1,...,vx2,vy2,...)
vector<double> TwoBody::evaluate(double t, vector<double> x){

//we know the vector x contains 4 informations for each mass (2 positions, 2 velocities)
double a = x.size()/4;

//we create the vectors that contain the positions and velocities of the two masses
vector<double> pos1(x.begin(), x.begin()+a);
vector<double> pos2(x.begin()+a, x.begin()+2*a);
vector<double> vel1(x.begin()+2*a, x.begin()+3*a);
vector<double> vel2(x.begin()+3*a, x.begin()+4*a);

//we compute the distance between the two masses
double distance = sqrt( (pos1-pos2)*(pos1-pos2) );

//create the vector that will contain the "differential equation"
vector<double> derivative;

//cicles that assignes the values of the derivative of the spatial coordinates of the 2 bodies
//(the spatial coordinates occupy the first 2*n places of the vector, where n is the dimension of the problem)
//dx/dt=v (dx1/dt=vx1,dy1/dt=vy1,...,dx2/dt=vx2,dy2/dt=vy2,...).
for(auto i : vel1){
	derivative.push_back(i);
}
for(auto i : vel2){
	derivative.push_back(i);
}

//cicle that assigns the values of the velocities
//(velocities occupy the last 2*n places of the vector)
//dv/dt=-(Gm/r^3)*(r1-r2) (dvx1/dt=-(Gm2/r^3)*(x1-x2), dvy1/dt=-(Gm2/r^3)*(y1-y2),...,dvx/dt=-(Gm1/r^3)*(x2-x1),...)   )

for(int i = 0; i < pos1.size(); ++i){
	derivative.push_back( -G*get_M2()/pow(distance, 3) * (pos1[i]-pos2[i]) );
}
for(int i = 0; i < pos2.size(); ++i){
	derivative.push_back( -G*get_M1()/pow(distance, 3) * (pos2[i]-pos1[i]) );
}

return derivative;

}
