all: main.x #odeint.x
main.x: main.o threeBody.o stm.o
	g++ main.o threeBody.o stm.o -O3 -std=c++11 -larmadillo -o main.x
main.o: main.cpp algebra.h threeBody.h vector2D.h in-out.h RK4.h pca.h
	g++ -c main.cpp -O3 -o main.o
threeBody.o: threeBody.cpp threeBody.h
	g++ -c threeBody.cpp -O3 -o threeBody.o
stm.o: stm.cpp stm.h
	g++ -c stm.cpp -O3 -o stm.o
clean:
	rm -f *.o *.x
