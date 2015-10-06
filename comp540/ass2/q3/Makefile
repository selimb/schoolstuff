objects = types.f90 globals.f90 setup.f90 utils.f90 matops.f90 linsolve.f90 main.f90
compile:
	gfortran $(objects)
check:
	gfortran -g -fbounds-check -Wall -fbacktrace -finit-real=nan $(objects)
clean:
	rm -fv *.out *.mod *.MOD
