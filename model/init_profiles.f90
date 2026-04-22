module init_profiles

    use poly_fit_mod
    use spline_eval

    implicit none

contains

    ! Initialize chemical potential profile across the grid.
    ! Supports linear gradient or uniform values set to left/right reservoir.
    subroutine init_mu_profile(mu, n, left_mu, right_mu, mode)
        implicit none

        real(8), intent(out) :: mu(:)
        integer, intent(in) :: n
        real(8), intent(in) :: left_mu, right_mu
        integer, intent(in) :: mode   ! 1=linear, 2=left, 3=right

        integer :: i

        ! Assign the reservoir chemical potential to the edges
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

    ! ! Compute density profile from chemical potential using polynomial fit.
    ! ! Applies relation rho(mu) on interior nodes (edges handled separately).
    ! subroutine compute_rho_from_mu(rho, mu, n, coeffs, degree)
    !     implicit none

    !     real(8), intent(out) :: rho(:)
    !     real(8), intent(in)  :: mu(:)
    !     real(8), intent(in)  :: coeffs(:)
    !     integer, intent(in)  :: n, degree

    !     integer :: block

    !     ! Important note: the value of rho for block=1 and block=n is wrong,
    !     ! in practice, it should be calculated from a theory accounting for entrance effects
    !     ! This will be done someday
    
    !     ! Estimate rho from mu
    !     do block = 1, n
    !         rho(block) = poly_fit(mu(block), coeffs, degree)
    !     end do

    ! end subroutine compute_rho_from_mu

    ! Compute density profile from chemical potential using spline fit.
    ! Applies relation rho(mu) on all nodes.
    subroutine compute_rho_from_mu(rho, mu, n, spl)
    use spline_data
    use spline_eval
    implicit none
    real(8),        intent(out) :: rho(:)
    real(8),        intent(in)  :: mu(:)
    integer,        intent(in)  :: n
    type(spline_t), intent(in)  :: spl
    integer :: block

    ! Important note: the value of rho for block=1 and block=n is wrong,
    ! in practice, it should be calculated from a theory accounting for entrance effects
    ! This will be done someday

    do block = 1, n
        rho(block) = eval_spline(mu(block), spl)
    end do

    end subroutine compute_rho_from_mu

end module init_profiles