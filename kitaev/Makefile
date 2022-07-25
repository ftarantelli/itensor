CC = g++ -g -Wall

all: sim.out

sim.out: kitaevcy.cpp
	$(CC) kitaevcy.cpp -osim.out -larmadillo -llapack -lblas
	./sim.out

show: plot
	evince ./data/graph.eps &
plot: data/
	cd data/
	gnuplot plot.plt
	cd ..

check: 
	$(CC) checkskitaevcy.cpp -osim.out -larmadillo -llapack -lblas
0check:
	$(CC) checkskitaevcy0.cpp -osim.out -larmadillo -llapack -lblas

clean:
	rm sim.out, ./data/graph.eps
	
	
fortran: out-exact
out-exact: disskit.f90
	touch ddissip.dat
	gfortran disskit.f90 lapack_tebd.lib.f -o out-exact
	./out-exact
	sudo rm out-exact ddissip.dat parameters.mod*

