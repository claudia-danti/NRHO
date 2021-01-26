#ifndef _pca_h_
#define _pca_h_


#include <armadillo>
#include<cmath>

using namespace std;
using namespace arma;

/*	Given a matrix X, its svd is:
	X = VSU^T
	S is a diagonal matrix with eigenvalues s_i.
	Armadillo's svd returns the vector s:i, this function sets to 0 the s_i<treshold,
	then it returns the product of the new matrix S with U and V: USV^T.

*/


//N:B.: this function is to be used with the pseudoinverse of the Jacobian

mat pca_inv(mat X, double treshold){
	//store the number of rows and cols of the matrix X
	double rows = X.n_rows;
	double cols = X.n_cols;

	mat U, V;
	vec s;

	//applying armadillo's svd
	svd(V,s,U,X);

	//sets to 0 the eigenvalues below the treshold
	for(int i = 0; i < s.n_elem; i++){
		if (s(i) < treshold){
			s(i) = 0.;
		}

	}

	mat Ut = trans(U);

	//create a diagonal matrix with the eigenvalues s
	mat S;
	//separate the underconstrained and overconstrained cases
	//overconstrained problem
	if(rows > cols){
			S = diagmat(s);
			S.insert_rows( rows, 1 );	//inserts at the end (number of row = rows = end of the matrix) one row of 0
	}

	//underconstrained problem

	else if(rows < cols){
		S = diagmat(s);
		S.insert_cols( cols, 1);	//adds a column of 0 at the lest column
	}


	return V*S*Ut;

}

//ALTERNATIVE

/*
	This function takes as argument the Jacobian (not its pseudoinverse)
	it is based upon the relation that occurs between the svd of a matrix and its inverse:
	If X is a matrix with svd:
		X=UDV^T

	then its pseudoinverse has svd:
		X^+=VSU^T

	where s_i,i=1/d_i,i (the eigenvalues of the S diagonal matrix are the inverse of the one of D).
*/


mat pca(mat X, double treshold){
	//store the number of rows and cols of the matrix X
	double rows = X.n_rows;
	double cols = X.n_cols;

	mat U, V;
	vec d;

	//applying armadillo's svd
	svd(U,d,V,X);

	//sets to 0 the eigenvalues below the treshold
	for(int i = 0; i < d.n_elem; i++){
		if (d(i) < treshold){
			d(i) = 0.;
		}

		else{
			d(i) = 1/d(i);
		}

	}

	mat Vt = trans(V);
	//create a diagonal matrix with the eigenvalues 1/d
	mat D;
	//separate the underconstrained and overconstrained cases
	//overconstrained problem
	if(rows > cols){
			D = diagmat(d);
			D.insert_rows( rows, 1 );	//inserts at the end (number of row = rows = end of the matrix) one row of 0
	}

	//underconstrained problem

	else if(rows < cols){
		D = diagmat(d);
		D.insert_cols( cols, 1);	//adds a column of 0 at the lest column
	}


	return U*D*Vt;

}


#endif
