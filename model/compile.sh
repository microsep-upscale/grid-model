gfortran -c coeff_io.f90
gfortran -c grid.f90
gfortran -o grid grid.o coeff_io.o -lm