output = "../bin/solver.out"
objects = types.f90 solvers.f90 main.f90
olevel = -g
compile:
	gfortran $(objects) $(olevel) -o $(output)
check:
	gfortran -g -fbounds-check -Wall -fbacktrace -finit-real=nan $(objects) -o $(output)
clean:
	rm -fv *.out *.mod *.MOD
