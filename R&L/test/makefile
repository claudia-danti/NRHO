all: testsvd.x 
testsvd.x: testsvd.o stm.o
	g++ testsvd.o stm.o -O3 -o testsvd.x
testsvd.o: testsvd.cpp pca.h stm.h
	g++ -c testsvd.cpp -O3 -o testsvd.o -std=c++11 -O2 -larmadillo
stm.o: stm.cpp stm.h vector2D.h
	g++ -c stm.cpp -O3 -o stm.o
clean:
	rm -f *.o *.x
