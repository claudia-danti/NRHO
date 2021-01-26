#include <iostream>
#include <armadillo>
#include "stm.h"
#include "pca.h"

using namespace std;
using namespace arma;


int main (){


mat Birm ={{1,2},{3,4},{5,6}};
 Birm.print("matrice B:");

  vec s = svd(Birm);
  s.print("svd:");

vec a;
mat U,V;

svd(U,a,V,Birm);

U.print("U:");
V.print("V");
a.print("singular values");





	return 0;
}