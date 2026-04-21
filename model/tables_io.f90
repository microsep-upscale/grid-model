module tables_io
  use poly_fit_mod
implicit none
contains

function invert_poly2(y, coeffs) result(mu)
  real(8), intent(in) :: y
  real(8), intent(in) :: coeffs(:)
  real(8) :: mu, a, b, c, discriminant
  a = coeffs(1); b = coeffs(2); c = coeffs(3)
  discriminant = b**2 - 4*a*(c - y)
  mu = (-b + sqrt(discriminant)) / (2*a)
end function invert_poly2

subroutine generate_tables(mu_min, mu_max, npts, rho_vs_mu, deg_rho, M_vs_mu, deg_M)
implicit none
real(8), intent(in) :: mu_min, mu_max
integer, intent(in) :: npts
real(8), intent(in) :: rho_vs_mu(:), M_vs_mu(:)
integer, intent(in) :: deg_rho, deg_M
real(8) :: mu, dmu
real(8) :: rho, M
real(8) :: mu_from_rho, mu_from_M
real(8) :: rho_back, M_back
real(8) :: rho_min, rho_max, drho
real(8) :: M_min, M_max, dM
integer :: i
integer :: u1, u2, u3, u4

  dmu = (mu_max - mu_min) / real(npts-1, 8)
  rho_min =  1d300;  rho_max = -1d300
  M_min   =  1d300;  M_max   = -1d300

  do i = 0, npts-1
    mu  = mu_min + i*dmu
    rho = poly_fit(mu, rho_vs_mu, deg_rho)
    M   = poly_fit(mu, M_vs_mu,   deg_M)
    if (rho < rho_min) rho_min = rho
    if (rho > rho_max) rho_max = rho
    if (M   < M_min)   M_min   = M
    if (M   > M_max)   M_max   = M
  end do

  drho = (rho_max - rho_min) / real(npts-1, 8)
  dM   = (M_max   - M_min)   / real(npts-1, 8)

  open(newunit=u1, file="output/table_rho_vs_mu.dat", status="replace")
  open(newunit=u2, file="output/table_M_vs_mu.dat",   status="replace")
  open(newunit=u3, file="output/table_mu_vs_rho.dat", status="replace")
  open(newunit=u4, file="output/table_mu_vs_M.dat",   status="replace")

  write(u1,*) "# mu[J/mol]    rho[1/m3]"
  write(u2,*) "# mu[J/mol]    M[s/kg/m]"
  write(u3,*) "# mu_from_rho[J/mol] rho_back[1/m3]"
  write(u4,*) "# mu_from_M[J/mol]   M_back[s/kg/m]"

  do i = 0, npts-1
    mu  = mu_min + i*dmu
    rho = poly_fit(mu, rho_vs_mu, deg_rho)
    M   = poly_fit(mu, M_vs_mu,   deg_M)
    write(u1,'(2ES20.10)') mu, rho
    write(u2,'(2ES20.10)') mu, M
  end do

  ! rho consistency (analytical inverse)
  do i = 0, npts-1
    rho         = rho_min + i*drho
    mu_from_rho = invert_poly2(rho, rho_vs_mu)
    rho_back    = poly_fit(mu_from_rho, rho_vs_mu, deg_rho)
    write(u3,'(2ES20.10)') mu_from_rho, rho_back
  end do

  ! M consistency (analytical inverse)
  do i = 0, npts-1
    M         = M_min + i*dM
    mu_from_M = invert_poly2(M, M_vs_mu)
    M_back    = poly_fit(mu_from_M, M_vs_mu, deg_M)
    write(u4,'(2ES20.10)') mu_from_M, M_back
  end do

  close(u1); close(u2); close(u3); close(u4)
end subroutine generate_tables
end module tables_io