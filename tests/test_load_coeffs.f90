program test_load_coeffs

  use coeff_io
  implicit none

  real(8), allocatable :: coeffs(:)
  integer :: degree, i

  call load_coeffs("../model/input/nemd_k_mu_coeffs.dat", coeffs, degree)

  print *, "degree =", degree
  print *, "coefficients:"

  do i = 1, size(coeffs)
     print *, i, coeffs(i)
  end do

end program test_load_coeffs