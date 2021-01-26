#ifndef vector2D_h
#define vector2D_h

#include <iostream>
#include <algorithm>
#include <vector>
#include <initializer_list>
#include <cassert>
#include <utility>
#include <numeric>
using namespace std;

template<typename T = double>
class vector2D{
public:
	vector2D() : rows(0), cols(0), v() {};
	vector2D(int r, int c) : rows(r), cols(c), v(rows*cols) {};
	vector2D(int r, int c, const T & initial_val) : rows(r), cols(c), v(rows*cols, initial_val) {};
	vector2D(int c, const vector<T> vec) : rows(vec.size() / c), cols(c), v(vec) {
		assert(rows * cols == vec.size() && "The construction cannot be done properly, wrong splitting of input vector");
		v.resize(rows * cols);
	}
	vector2D(const vector<vector<T>> & vecofvecs) : rows(vecofvecs.size()), cols(vecofvecs[0].size()), v(rows*cols) {
		for(auto & i : vecofvecs){
			assert( cols == i.size() && "vector<vector<T>> with inhomogeneous length");
		}
		for(int i = 0; i < v.size(); ++i){
			v[i] = vecofvecs[i / cols][i % cols];
		}
	}
	vector2D(initializer_list<initializer_list<T>> list) : rows(list.size()), cols((*list.begin()).size()), v(rows*cols) {
		int i = 0;
		for(auto it_row = list.begin(); it_row != list.end(); ++it_row){
			assert( cols == (*it_row).size() && "initializer_list with inhomogeneous length" );
			for(auto it_col = (*it_row).begin(); it_col != (*it_row).end(); ++it_col){
				v[i] = *it_col;
				++i;
			}
		}
	}
	int row() const { return rows; }
	int col() const { return cols; }
	T & element(int r, int c){
		return v[r*cols + c];
	}
	T element(int r, int c) const{
		return v[r*cols + c];
	}
	vector<T> & getLinearVector(){ return v; };
	vector<T> getLinearVector() const{ return v; };
	vector<T> operator[](int r) const{
		return vector<T>(v.begin() + r*cols, v.begin() + (r+1)*cols);
	}
	vector<T> getRow(int r) const{
		return vector<T>(v.begin() + r*cols, v.begin() + (r+1)*cols);
	}
	vector<T> getCol(int c) const{
		vector<T> to_return(rows);
		for(int i = 0; i < to_return.size(); ++i){
			to_return[i] = v[i*cols + c];
		}
		return to_return;
	}
	void setRow(int r, const vector<T> & input){
		if( checkColSize(input.size()) ){
			for(int i = 0; i < input.size(); ++i){
				v[r*cols + i] = input[i];
			}
		}
	}
	void setCol(int c, const vector<T> & input){
		if( checkRowSize(input.size()) ){
			for(int i = 0; i < input.size(); ++i){
				v[i*cols + c] = input[i];
			}
		}
	}
	void insertRow(int r, const vector<T> & input){
		if( checkColSize(input.size()) ){
			++rows;
			v.insert(v.begin() + r*cols, input.begin(), input.end());
		}
	}
	void insertCol(int c, const vector<T> & input){
		if( checkRowSize(input.size()) ){
			++cols;
			int index_old = v.size()-1;
			int index_input = input.size()-1;
			v.resize(rows*cols);
			for(int i = 0; i < rows*cols; ++i){
				int index = rows*cols - 1 - i;
				if( index % cols == c ){
					v[index] = input[index_input];
					--index_input;
				} else{
					v[index] = v[index_old];
					--index_old;
				}
			}
		}
	}
	void eraseRow(int r){
		v.erase(v.begin() + r*cols, v.begin() + (r+1)*cols);
		--rows;
	}
	void eraseCol(int c){
		--cols;
		int erased = 0;
		for(int i = c; i < v.size(); ++i){
			if( i % cols == c ){
				++erased;
			}
			v[i] = v[i+erased];
		}
		v.resize(rows*cols);
	}
	void pushBackRow(const vector<T> & input){
		if( checkColSize(input.size()) ){
			++rows;
			v.insert(v.end(), input.begin(), input.end());
		}
	}
	void pushBackCol(const vector<T> & input){
		if( checkRowSize(input.size()) ){
			++cols;
			int index_old = v.size()-1;
			int index_input = input.size()-1;
			v.resize(rows*cols);
			for(int i = 0; i < rows*cols; ++i){
				int index = rows*cols - 1 - i;
				if( index % cols == cols - 1 ){
					v[index] = input[index_input];
					--index_input;
				} else{
					v[index] = v[index_old];
					--index_old;
				}
			}
		}
	}
	vector2D<T> transposed() const{
		vector2D to_return(cols, rows);
		for(int i = 0; i < to_return.row(); ++i){
			for(int j = 0; j < to_return.col(); ++j){
				to_return.element(i,j) = element(j, i);
			}
		}
		return to_return;
	}
	void print(ostream & os) const{
		for(int i = 0; i < rows; ++i){
			for(int j = 0; j < cols; ++j){
				os << v[i*cols + j] << "\t";
			}
			os << endl;
		}
	}
	void printTransposed(ostream & os) const{
		for(int j = 0; j < cols; ++j){
			for(int i = 0; i < rows; ++i){
				os << v[i*cols + j] << "\t";
			}
			os << endl;
		}
	}
protected:
	int rows;
	int cols;
	vector<T> v;
private:
	bool checkColSize(int inputsize){
		if( v.size() == 0 ){
			cols = inputsize;
			return true;
		}
		if( cols == inputsize ){
			return true;
		} else{
			cerr << "The insertion cannot be done, wrong number of components in the inserting vector" << endl;
			return false;
		}
	}
	bool checkRowSize(int inputsize){
		if( v.size() == 0 ){
			rows = inputsize;
			return true;
		}
		if( rows == inputsize ){
			return true;
		} else{
			cerr << "The insertion cannot be done, wrong number of components in the inserting vector" << endl;
			return false;
		}
	}
};

template<typename T>
vector2D<T> operator+(const vector2D<T> & a , const vector2D<T> & b){


	if(a.col() == b.col() && a.row() == b.row()){

		vector<T> a_lin = a.getLinearVector();
		vector<T> b_lin = b.getLinearVector();
		vector<T> sum; 

		for(int i = 0; i < a_lin.size(); i++){
			sum[i] = a_lin [i] + b_lin[i];
		}
		vector2D<T> sum_mat (a.col(), sum);

		return sum_mat;
	}

	else{
		cerr<<"You cannot sum two matrices that have different dimensions."<<endl;
		return vector2D<T> ();
	}

}


template<typename T>
vector2D<T> operator-(const vector2D<T> & a , const vector2D<T> & b){


	if(a.col() == b.col() && a.row() == b.row()){

		vector<T> a_lin = a.getLinearVector();
		vector<T> b_lin = b.getLinearVector();
		vector<T> diff;

		for(int i = 0; i < a_lin.size(); i++){
			diff[i] = a_lin [i] - b_lin[i];
		}

		vector2D<T> diff_mat (a.col(), diff);

		return diff_mat;
	}

	else{
		cerr<<"You cannot subract two matrices that have different dimensions."<<endl;

		return vector2D<T> ();
	}

}

template<typename T>
vector2D<T> operator/(const vector2D<T> & a , T b){

	 vector<T> a_lin = a.getLinearVector();
	 vector<T> a_lin_b(a_lin.size());

	 for(int i = 0; i < a_lin.size(); i++){
	 	a_lin_b[i] = a_lin[i]/b;
	 }

	 vector2D<T> c (a.col(), a_lin_b);

	 return c;

}

template<typename T>
vector<T> operator*(const vector2D<T> & a, T b){

	vector<T> a_lin = a.getLinearVector();
	vector<T> a_lin_b(a_lin.size());

	for(int i = 0; i < a_lin.size(); i++){
	 	a_lin_b[i] = a_lin[i]*b;
	}

	vector2D<T> c (a.col(), a_lin_b);

	return c;
}



template<typename T>
vector<T> operator*(const vector2D<T> & m, const vector<T> & v){
	assert( m.col() == v.size() && "incompatibility betweeen vector2D cols and vector rows");
	vector<T> to_return(m.row());
	for(int i = 0; i < m.row(); ++i){
		to_return[i] = 0;
		for(int j = 0; j < m.col(); ++j){
			to_return[i] += m.element(i,j) * v[j];
		}
	}
	return to_return;
}


template<typename T>
vector<T> operator*(const vector<T> & v, const vector2D<T> & m){
	return m.transposed() * v;
}


template<typename T>
vector2D<T> operator*(const vector2D<T> & A, const vector2D<T> & B){
	assert( A.row() == B.col() && A.col() == B.row() && "impossible to multiply matrix with different size");
	vector2D<T> to_return(A.row(), A.col());
	vector<double> temp1, temp2;
	for(int i = 0; i < to_return.row(); ++i){
		for(int j = 0; j < to_return.col(); ++j){
			temp1 = A.getRow(i);
			temp2 = B.getCol(j);
			to_return.element(i,j) = inner_product(temp1.begin(), temp1.end(), temp2.begin(), 0.);
		}
	}
	return to_return;
}


template<typename T>
ostream & operator<<(ostream & os, const vector2D<T> & v2d){
	os << endl;
	for(int i = 0; i < v2d.row(); ++i){
		for(int j = 0; j < v2d.col(); ++j){
			if( j != 0 ){
				os << "\t";
			}
			os << v2d.element(i,j);
		}
		os << endl;
	}
	return os;
}


//prodotto colonna per riga
template<typename pointType>
vector2D<double> outerproduct (pointType a, pointType b){
	
	int dim_a = a.size();

	vector2D<double> c (a.size(), b.size());

	for(int i = 0; i < dim_a; i++){
		c.setRow(i, b*a[i]);
	}

	return c;

}







#endif