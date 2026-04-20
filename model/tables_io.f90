module tables_io
  use poly_fit_mod
  implicit none
contains

  subroutine generate_tables(mu_min, mu_max, npts, &
                             rho_coeffs, deg_rho, &
                             k_coeffs, deg_k)

    implicit none

    real(8), intent(in) :: mu_min, mu_max
    integer, intent(in) :: npts
    real(8), intent(in) :: rho_coeffs(:), k_coeffs(:)
    integer, intent(in) :: deg_rho, deg_k

    real(8) :: mu, dmu, rho, K
    integer :: i, u1, u2

    dmu = (mu_max - mu_min) / real(npts-1,8)

    open(newunit=u1, file="output/rho_vs_mu_md.dat", status="replace")
    open(newunit=u2, file="output/k_vs_mu_md.dat",   status="replace")

    write(u1,*) "# mu[J/mol] rho[kg/m3]"
    write(u2,*) "# mu[J/mol] K[mol/J.m.s]"

    do i = 0, npts-1

      mu = mu_min + i*dmu

      rho = poly_fit(mu, rho_coeffs, deg_rho)
      K   = poly_fit(mu, k_coeffs, deg_k)

      write(u1,'(2ES20.10)') mu, rho
      write(u2,'(2ES20.10)') mu, K

    end do

    close(u1)
    close(u2)

  end subroutine generate_tables

end module tables_io