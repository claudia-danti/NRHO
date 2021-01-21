all: main.x #odeint.x
main.x: main.o threeBody.o threebodyphi.o
	g++ main.o threeBody.o threebodyphi.o -O3 -o main.x
main.o: main.cpp algebra.h threeBody.h vector2D.h in-out.h RK4.h
	g++ -c main.cpp -O3 -o main.o
threeBody.o: threeBody.cpp threeBody.h twoBody.h
	g++ -c threeBody.cpp -O3 -o threeBody.o
twoBody.o: twoBody.cpp twoBody.h
	g++ -c twoBody.cpp -O3 -o twoBody.o
threebodyphi.o: threebodyphi.cpp threebodyphi.h
	g++ -c threebodyphi.cpp -O3 -o threebodyphi.o
clean:
	rm -f *.o *.x
