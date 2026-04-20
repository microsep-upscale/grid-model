#!/bin/bash

mkdir -p output

gfortran -c coeff_io.f90
gfortran -c poly_fit.f90
gfortran -c grid.f90

gfortran coeff_io.o poly_fit.o grid.o -o grid

echo "Build complete"