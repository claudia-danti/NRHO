#include <iostream>
#include <armadillo>
#include "stm.h"
#include "pca.h"

using namespace std;
using namespace arma;


int main (){


/*mat Birm ={{1,2},{3,4},{5,6}};
Birm.print("matrice B:");

vec s = svd(Birm);
s.print("svd autovalori:");

vec a;
mat U,V;

svd(U,a,V,Birm);

U.print("U:");
V.print("V");
a.print("singular values");
Birm.print("Birm after svd");
*/
double treshold = pow(10,-4);


mat B = {{1,22,1,2},{35,4,2,2},{51,4,50,2},{1,23,2,4},{7,4,9,9}};

B.print("matrix B");

mat p = pca(B, treshold);

p.print("matrice post pca:");

	return 0;
}