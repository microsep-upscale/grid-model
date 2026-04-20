#!/bin/bash

set -e

mkdir -p output

# clean old build artifacts
rm -f *.o *.mod grid

gfortran -c coeff_io.f90
gfortran -c poly_fit.f90
gfortran -c io_profiles.f90
gfortran -c init_profiles.f90
gfortran -c tables_io.f90
gfortran -c grid.f90

gfortran coeff_io.o poly_fit.o io_profiles.o init_profiles.o tables_io.o grid.o -o grid

echo "Build complete"