module mu_rho_inversion

    use density_interpolation
    use io_profiles

    implicit none

    private
    public :: invert_mu_from_rho, invert_mu_grid_validate

contains

    ! rho(muA,muB) + Jacobian via central finite differences on the
    ! log-transformed bicubic surface (matches how rhoA_fine/rhoB_fine
    ! were built in interpolate_density_grids).
    subroutine eval_rho_and_jacobian(muA_values, muB_values, rhoA_log_grid, rhoB_log_grid, &
                                      muA_q, muB_q, rhoA_q, rhoB_q, &
                                      dRhoA_dMuA, dRhoA_dMuB, dRhoB_dMuA, dRhoB_dMuB)
        real(8), intent(in)  :: muA_values(:), muB_values(:)
        real(8), intent(in)  :: rhoA_log_grid(:,:), rhoB_log_grid(:,:)
        real(8), intent(in)  :: muA_q, muB_q
        real(8), intent(out) :: rhoA_q, rhoB_q
        real(8), intent(out) :: dRhoA_dMuA, dRhoA_dMuB, dRhoB_dMuA, dRhoB_dMuB

        real(8) :: hA, hB
        real(8) :: rhoA_pA, rhoA_mA, rhoA_pB, rhoA_mB
        real(8) :: rhoB_pA, rhoB_mA, rhoB_pB, rhoB_mB

        hA = 1.0d-4 * (muA_values(size(muA_values)) - muA_values(1))
        hB = 1.0d-4 * (muB_values(size(muB_values)) - muB_values(1))

        rhoA_q = exp(bicubic_eval_point(muA_values, muB_values, rhoA_log_grid, muA_q, muB_q))
        rhoB_q = exp(bicubic_eval_point(muA_values, muB_values, rhoB_log_grid, muA_q, muB_q))

        rhoA_pA = exp(bicubic_eval_point(muA_values, muB_values, rhoA_log_grid, muA_q+hA, muB_q))
        rhoA_mA = exp(bicubic_eval_point(muA_values, muB_values, rhoA_log_grid, muA_q-hA, muB_q))
        rhoA_pB = exp(bicubic_eval_point(muA_values, muB_values, rhoA_log_grid, muA_q, muB_q+hB))
        rhoA_mB = exp(bicubic_eval_point(muA_values, muB_values, rhoA_log_grid, muA_q, muB_q-hB))

        rhoB_pA = exp(bicubic_eval_point(muA_values, muB_values, rhoB_log_grid, muA_q+hA, muB_q))
        rhoB_mA = exp(bicubic_eval_point(muA_values, muB_values, rhoB_log_grid, muA_q-hA, muB_q))
        rhoB_pB = exp(bicubic_eval_point(muA_values, muB_values, rhoB_log_grid, muA_q, muB_q+hB))
        rhoB_mB = exp(bicubic_eval_point(muA_values, muB_values, rhoB_log_grid, muA_q, muB_q-hB))

        dRhoA_dMuA = (rhoA_pA - rhoA_mA) / (2d0*hA)
        dRhoA_dMuB = (rhoA_pB - rhoA_mB) / (2d0*hB)
        dRhoB_dMuA = (rhoB_pA - rhoB_mA) / (2d0*hA)
        dRhoB_dMuB = (rhoB_pB - rhoB_mB) / (2d0*hB)

    end subroutine eval_rho_and_jacobian

    ! Newton solve for (muA,muB) given target (rhoA,rhoB). Clamps to the
    ! table domain each step to prevent divergence outside the known range.
    subroutine invert_mu_from_rho(muA_values, muB_values, rhoA_log_grid, rhoB_log_grid, &
                                   rhoA_target, rhoB_target, muA0, muB0, &
                                   muA, muB, converged, n_iter_used, tol, maxiter)
        real(8), intent(in)  :: muA_values(:), muB_values(:)
        real(8), intent(in)  :: rhoA_log_grid(:,:), rhoB_log_grid(:,:)
        real(8), intent(in)  :: rhoA_target, rhoB_target
        real(8), intent(in)  :: muA0, muB0
        real(8), intent(out) :: muA, muB
        logical, intent(out) :: converged
        integer, intent(out) :: n_iter_used
        real(8), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter

        real(8) :: tol_, rhoA_q, rhoB_q
        integer :: maxiter_, it
        real(8) :: resA, resB, J11, J12, J21, J22, det, dmuA, dmuB
        real(8) :: muA_min, muA_max, muB_min, muB_max

        tol_ = 1d-10
        if (present(tol)) tol_ = tol
        maxiter_ = 50
        if (present(maxiter)) maxiter_ = maxiter

        muA_min = muA_values(1);  muA_max = muA_values(size(muA_values))
        muB_min = muB_values(1);  muB_max = muB_values(size(muB_values))

        muA = muA0
        muB = muB0
        converged = .false.
        n_iter_used = 0

        do it = 1, maxiter_
            call eval_rho_and_jacobian(muA_values, muB_values, rhoA_log_grid, rhoB_log_grid, &
                                        muA, muB, rhoA_q, rhoB_q, J11, J12, J21, J22)
            resA = rhoA_q - rhoA_target
            resB = rhoB_q - rhoB_target

            if (abs(resA) < tol_ .and. abs(resB) < tol_) then
                converged = .true.
                n_iter_used = it
                return
            end if

            det = J11*J22 - J12*J21
            if (abs(det) < 1.0d-300) then
                n_iter_used = it
                return
            end if

            dmuA = -( J22*resA - J12*resB) / det
            dmuB = -(-J21*resA + J11*resB) / det

            muA = min(max(muA + dmuA, muA_min), muA_max)
            muB = min(max(muB + dmuB, muB_min), muB_max)
        end do

        n_iter_used = maxiter_
    end subroutine invert_mu_from_rho

    ! Validation pass: invert every point of the fine density grid back to
    ! (muA,muB), compare with the known MU_A_fine/MU_B_fine, print and
    ! write out error statistics -- Fortran analogue of the Python check.
    subroutine invert_mu_grid_validate(muA_values, muB_values, rhoA_log_grid, rhoB_log_grid, &
                                        muA_fine, muB_fine, rhoA_fine, rhoB_fine, output_dir)
        real(8), intent(in) :: muA_values(:), muB_values(:)
        real(8), intent(in) :: rhoA_log_grid(:,:), rhoB_log_grid(:,:)
        real(8), intent(in) :: muA_fine(:), muB_fine(:)
        real(8), intent(in) :: rhoA_fine(:,:), rhoB_fine(:,:)
        character(len=*), intent(in) :: output_dir

        integer :: i, j, ni, nj, n_iter_used, n_conv
        real(8) :: muA_r, muB_r
        logical :: conv
        real(8), allocatable :: muA_error(:,:), muB_error(:,:)
        real(8) :: max_errA, max_errB, rms_errA, rms_errB
        integer :: n_valid

        ni = size(muA_fine)
        nj = size(muB_fine)
        allocate(muA_error(ni,nj), muB_error(ni,nj))
        muA_error = 0d0
        muB_error = 0d0
        n_conv = 0
        n_valid = 0
        max_errA = 0d0; max_errB = 0d0
        rms_errA = 0d0; rms_errB = 0d0

        do j = 1, nj
            do i = 1, ni
                ! perturb initial guess away from the true answer (real test)
                call invert_mu_from_rho(muA_values, muB_values, rhoA_log_grid, rhoB_log_grid, &
                                         rhoA_fine(i,j), rhoB_fine(i,j), &
                                         muA_fine(i) + 0.05d0, muB_fine(j) - 0.05d0, &
                                         muA_r, muB_r, conv, n_iter_used)
                if (conv) then
                    n_conv = n_conv + 1
                    n_valid = n_valid + 1
                    muA_error(i,j) = muA_r - muA_fine(i)
                    muB_error(i,j) = muB_r - muB_fine(j)
                    max_errA = max(max_errA, abs(muA_error(i,j)))
                    max_errB = max(max_errB, abs(muB_error(i,j)))
                    rms_errA = rms_errA + muA_error(i,j)**2
                    rms_errB = rms_errB + muB_error(i,j)**2
                end if
            end do
        end do

        if (n_valid > 0) then
            rms_errA = sqrt(rms_errA / real(n_valid,8))
            rms_errB = sqrt(rms_errB / real(n_valid,8))
        end if

        ! write(*,*) "Converged fraction:", real(n_conv,8) / real(ni*nj,8)
        ! write(*,*) "Max absolute error in muA:", max_errA
        ! write(*,*) "RMS error in muA:", rms_errA
        ! write(*,*) "Max absolute error in muB:", max_errB
        ! write(*,*) "RMS error in muB:", rms_errB

        call write_grid2d(trim(output_dir)//"/muA_error.dat", muA_error, "# muA reconstruction error")
        call write_grid2d(trim(output_dir)//"/muB_error.dat", muB_error, "# muB reconstruction error")

        deallocate(muA_error, muB_error)
    end subroutine invert_mu_grid_validate

end module mu_rho_inversion