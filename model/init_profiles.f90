module init_profiles

    use poly_fit_mod

    implicit none

contains

  subroutine init_mu_profile(mu, x_axis, n, left_mu, right_mu, mode)
    implicit none

    real(8), intent(out) :: mu(:)
    real(8), intent(in)  :: x_axis(:)
    integer, intent(in)   :: n
    real(8), intent(in)   :: left_mu, right_mu
    integer, intent(in)   :: mode   ! 1=linear, 2=left, 3=right

    integer :: i

    mu(1) = left_mu
    mu(n) = right_mu

    do i = 2, n-1

      select case (mode)

      case (1)
        ! linear gradient
        mu(i) = left_mu + (right_mu - left_mu) * &
                (i-1) / real(n-1,8)

      case (2)
        ! left reservoir
        mu(i) = left_mu

      case (3)
        ! right reservoir
        mu(i) = right_mu

      case default
        stop "Unknown mu initialization mode"

      end select

    end do

  end subroutine init_mu_profile

  subroutine init_rho_from_mu(rho, mu, n, coeffs, degree)
    implicit none

    real(8), intent(out) :: rho(:)
    real(8), intent(in)  :: mu(:)
    real(8), intent(in)  :: coeffs(:)
    integer, intent(in)  :: n, degree

    integer :: i

    rho(1) = 0
    rho(n) = 0
 
    do i = 2, n-1
      rho(i) = poly_fit(mu(i), coeffs, degree)
    end do

  end subroutine init_rho_from_mu

end module init_profiles