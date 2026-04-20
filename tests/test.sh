#!/bin/bash

gfortran -c ../model/coeff_io.f90 -J .
gfortran test_load_coeffs.f90 coeff_io.o -o test_load

./test_load