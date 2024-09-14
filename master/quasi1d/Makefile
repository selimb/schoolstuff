output = "q1d.out"
objects = types.f constants.f inputs.f setup.f utils.f common_calcs.f flx_schemes.f bc.f timestepping.f main.f
compile:
	gfortran $(objects) -o $(output)
check:
	gfortran -g -fbounds-check -Wall -fbacktrace -finit-real=nan $(objects)
clean:
	rm -fv *.out *.mod *.MOD
