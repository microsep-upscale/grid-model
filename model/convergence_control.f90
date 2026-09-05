module convergence_control

  implicit none
  public :: check_flux_conservation, check_steady_state

contains

  ! Stopping criterion for mode 1 (linear gradient / driven flux).
  ! At steady state all interior edge fluxes should be equal, so the
  ! spread (std/mean) of the interior fluxes is used as the score.
  ! Not meaningful for mode 2, where flux_mean -> 0 by symmetry and
  ! the ratio diverges instead of vanishing.
  subroutine check_flux_conservation(iter, check_interval, time_step, &
                                      flux_edges, number_block, conv_tol, &
                                      conv_unit, flux_mean, flux_std, &
                                      flux_conservation, converged)
    implicit none
    integer(kind=8), intent(in)  :: iter
    integer,         intent(in)  :: check_interval
    real(8),         intent(in)  :: time_step
    real(8),         intent(in)  :: flux_edges(:)
    integer,         intent(in)  :: number_block
    real(8),         intent(in)  :: conv_tol
    integer,         intent(in)  :: conv_unit
    real(8),         intent(out) :: flux_mean, flux_std, flux_conservation
    logical,         intent(out) :: converged

    integer :: block, edge

    flux_mean = 0d0
    do block = 2, number_block-1
      edge = block - 1
      flux_mean = flux_mean + flux_edges(edge)
    end do
    flux_mean = flux_mean / real(number_block-2, 8)

    flux_std = 0d0
    do block = 2, number_block-1
      edge = block - 1
      flux_std = flux_std + (flux_edges(edge) - flux_mean)**2
    end do
    flux_std = sqrt(flux_std / real(number_block-2, 8))

    ! normalized score: 0 = perfectly conserved, 1 = large variation
    if (abs(flux_mean) > 0d0) then
      flux_conservation = flux_std / abs(flux_mean)
    else
      flux_conservation = 0d0
    end if

    converged = .false.

    if (flux_conservation > 10d0) then
      write(*,*) "time_step =", time_step
      write(*,*) flux_conservation, flux_std, flux_mean
      write(*,*) flux_edges
      stop "Conservation score diverged — simulation aborted"
    end if

    if (flux_conservation < conv_tol .and. iter > check_interval) then
      converged = .true.
    end if

    write(conv_unit,'(I10,4ES20.10)') iter, flux_conservation, flux_mean, flux_std, time_step

  end subroutine check_flux_conservation

  ! Stopping criterion for mode 2 (symmetric reservoirs) and, in
  ! general, any case where flux_mean is not a reliable score. Tracks
  ! the mean interior density and stops once its fractional drift
  ! between checks stays below steady_tol for n_steady consecutive
  ! checks. steady_count and prev_mean_density persist across calls
  ! (inout); initialize steady_count = 0 and prev_mean_density = -1d0
  ! before the first call. Set steady_tol <= 0  disable.
  subroutine check_steady_state(iter, fluid_density, &
                                 number_block, steady_tol, n_steady, &
                                 time_step, conv_unit, &
                                 steady_count, prev_mean_density, converged)
    implicit none
    integer(kind=8), intent(in)    :: iter
    real(8),         intent(in)    :: fluid_density(:)
    integer,         intent(in)    :: number_block
    real(8),         intent(in)    :: steady_tol
    integer,         intent(in)    :: n_steady
    real(8),         intent(in)    :: time_step
    integer,         intent(in)    :: conv_unit
    integer,         intent(inout) :: steady_count
    real(8),         intent(inout) :: prev_mean_density
    logical,         intent(out)   :: converged

    real(8) :: mean_density, density_drift
    integer :: block

    converged = .false.

    if (steady_tol <= 0d0) return

    mean_density = 0d0
    do block = 2, number_block-1
      mean_density = mean_density + fluid_density(block)
    end do
    mean_density = mean_density / real(number_block-2, 8)

    density_drift = 0d0
    if (prev_mean_density > 0d0) then
      density_drift = abs(mean_density - prev_mean_density) / abs(prev_mean_density)
      if (density_drift < steady_tol) then
        steady_count = steady_count + 1
      else
        steady_count = 0
      end if
    end if
    prev_mean_density = mean_density

    if (steady_count >= n_steady) then
      converged = .true.
    end if

    write(conv_unit,'(I10,2ES20.10)') iter, density_drift, time_step

  end subroutine check_steady_state

end module convergence_control