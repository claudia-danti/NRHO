#ifndef _pca_h_
#define _pca_h_


#include <armadillo>
#include<cmath>

using namespace std;
using namespace arma;

/*
	This function takes as argument the Jacobian and returns the pseudoinverse after svd
	it is based upon the relation that occurs between the svd of a matrix and its inverse:
	If X is a matrix with svd:
		X=UDV^T

	then its pseudoinverse has svd:
		X^+=VSU^T

	where s_i,i=1/d_i,i (the eigenvalues of the S diagonal matrix are the inverse of the one of D).
*/

template<typename matrix>
matrix pca(matrix X, double threshold){

	//define the matrices used in the svd
	matrix U, V;
	vec s;

	//applying armadillo's svd
	svd(U,s,V,X);
cout<<"autovalori"<<1/s<<endl;
	//create a diagonal matrix with the eigenvalues 1/s
	int S_rows = V.n_cols;
	int S_cols = U.n_cols;
	matrix S (S_rows, S_cols, fill::zeros);

	//sets to 0 the eigenvalues below the treshold
	for(int i = 0; i < s.n_elem; i++){
		if (s(i) < threshold){
			S(i,i) = 0.;
		}

		else{
			S(i,i) = 1/s(i);
		}

	}

	return V*S*trans(U);

}


#endif
