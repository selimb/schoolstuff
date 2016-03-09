output = "../bin/solver.out"
objects = modules.f90 solver.f90 post.f90 main.f90
# olevel = -g
olevel = -O3
opt =
compile:
	gfortran $(objects) $(opt) $(olevel) -o $(output)
grid:
	gfortran generategrid.f90 $(opt) $(olevel) -o "../bin/gridder.out"
test:
	gfortran -g -C modules.f90 solver.f90 testslvtridiag.f90 -O3 -o "test.out"
check:
	gfortran -g -fbounds-check -Wall -fbacktrace -finit-real=nan $(objects) -o $(output)
clean:
	rm -fv *.out *.mod *.MOD
