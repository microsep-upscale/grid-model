#!/bin/bash

set -e

mkdir -p output

# clean old build artifacts
rm -f *.o *.mod grid

# compiler + flags
FC=gfortran
FCFLAGS="-g" # "-g -O0 -Wall -Wextra -Wunused -Wimplicit-interface -Wsurprising -fcheck=all -fbacktrace"

# compile sources
$FC $FCFLAGS -c coeff_io.f90
$FC $FCFLAGS -c spline_data.f90
$FC $FCFLAGS -c spline_io.f90
$FC $FCFLAGS -c spline_eval.f90
$FC $FCFLAGS -c poly_fit.f90
$FC $FCFLAGS -c io_profiles.f90
$FC $FCFLAGS -c init_profiles.f90
$FC $FCFLAGS -c tables_io.f90
$FC $FCFLAGS -c grid.f90

# link
$FC coeff_io.o poly_fit.o io_profiles.o init_profiles.o tables_io.o spline_io.o spline_eval.o spline_data.o grid.o -o grid
