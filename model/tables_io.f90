module tables_io
  use poly_fit_mod
  implicit none

contains

  subroutine generate_tables(mu_min, mu_max, npts, &
                             rho_coeffs, deg_rho, &
                             k_coeffs, deg_k, &
                             mu_vs_rho, deg_mu_rho, &
                             mu_vs_k, deg_mu_k)

    implicit none

    real(8), intent(in) :: mu_min, mu_max
    integer, intent(in) :: npts

    real(8), intent(in) :: rho_coeffs(:), k_coeffs(:)
    real(8), intent(in) :: mu_vs_rho(:), mu_vs_k(:)

    integer, intent(in) :: deg_rho, deg_k
    integer, intent(in) :: deg_mu_rho, deg_mu_k

    real(8) :: mu, dmu
    real(8) :: rho, K
    real(8) :: mu_from_rho, mu_from_k

    integer :: i
    integer :: u1, u2, u3, u4

    ! ----------------------------
    ! step size
    ! ----------------------------
    dmu = (mu_max - mu_min) / real(npts-1,8)

    ! ----------------------------
    ! open files
    ! ----------------------------
    open(newunit=u1, file="output/rho_vs_mu.dat", status="replace")
    open(newunit=u2, file="output/k_vs_mu.dat", status="replace")
    open(newunit=u3, file="output/mu_vs_rho.dat", status="replace")
    open(newunit=u4, file="output/mu_vs_k.dat", status="replace")

    write(u1,*) "# mu[J/mol] rho[kg/m3]"
    write(u2,*) "# mu[J/mol] K"
    write(u3,*) "# rho[kg/m3] mu[J/mol]"
    write(u4,*) "# K mu[J/mol]"

    ! ----------------------------
    ! forward sweep in mu
    ! ----------------------------
    do i = 0, npts-1

      mu = mu_min + i*dmu

      rho = poly_fit(mu, rho_coeffs, deg_rho)
      K   = poly_fit(mu, k_coeffs, deg_k)

      write(u1,'(2ES20.10)') mu, rho
      write(u2,'(2ES20.10)') mu, K

    end do

    ! ----------------------------
    ! inverse consistency sweep
    ! ----------------------------
    do i = 0, npts-1

      ! sample rho-space and K-space consistently
      mu = mu_min + i*dmu   ! reuse spacing just for sampling

      ! invert using polynomials
      mu_from_rho = poly_fit(mu, mu_vs_rho, deg_mu_rho)
      mu_from_k   = poly_fit(mu, mu_vs_k, deg_mu_k)

      write(u3,'(2ES20.10)') mu, mu_from_rho
      write(u4,'(2ES20.10)') mu, mu_from_k

    end do

    close(u1)
    close(u2)
    close(u3)
    close(u4)

  end subroutine generate_tables

end module tables_io