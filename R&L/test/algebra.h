#ifndef algebra_h
#define algebra_h

#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>


using namespace std;


template<typename T>
short sign(T val){
	if( val > 0 ){
		return 1;
	} else if( val == 0 ){
		return 0;
	}
	return -1;
}


template<typename pointType>
void normalize(pointType & v, double radius = 1.){
	double norm = sqrt(inner_product(v.begin(), v.end(), v.begin(), 0.));
	for(auto & i : v){
		i = i * radius / norm;
	}
	return;
}


template<typename pointType>
double norm(const pointType & v){
	return sqrt(inner_product(v.begin(), v.end(), v.begin(), 0.));
}


template<typename pointType>
pointType operator+(pointType a, const pointType & b){
	for(int i = 0; i < a.size(); ++i){
		a[i] = a[i] + b[i];
	}
	return a;
}


template<typename pointType>
pointType operator-(pointType a, const pointType & b){
	for(int i = 0; i < a.size(); ++i){
		a[i] = a[i] - b[i];
	}
	return a;
}


template<typename pointType>
pointType operator*(pointType a, double b){
	for(int i = 0; i < a.size(); ++i){
		a[i] = a[i] * b;
	}
	return a;
}
template<typename pointType>
pointType operator*(double b, pointType a){
	return a * b;
}


template<typename pointType>
pointType operator/(pointType a, double b){
	for(int i = 0; i < a.size(); ++i){
		a[i] = a[i] / b;
	}
	return a;
}


template<typename pointType>
double operator*(const pointType & a, const pointType & b){
	return inner_product(a.begin(), a.end(), b.begin(), 0.);
}


#endif
