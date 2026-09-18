module density_interpolation
    use spline_data
    use spline_eval
    use io_profiles

    implicit none

    private
    public :: interpolate_density_grids

contains

    subroutine interpolate_density_grids(muA_values, muB_values, rhoA_grid, rhoB_grid, &
                                         multiplication, muA_fine, muB_fine, &
                                         rhoA_fine, rhoB_fine, output_dir)
        implicit none

        ! Inputs
        real(8), intent(in) :: muA_values(:), muB_values(:)
        real(8), intent(in) :: rhoA_grid(:,:), rhoB_grid(:,:)
        integer, intent(in) :: multiplication
        character(len=*), intent(in) :: output_dir

        ! Outputs
        real(8), allocatable, intent(out) :: muA_fine(:), muB_fine(:)
        real(8), allocatable, intent(out) :: rhoA_fine(:,:), rhoB_fine(:,:)

        ! Local variables
        integer :: i, j, k, n_muA, n_muB, grid_size
        real(8), allocatable :: rhoA_log_grid(:,:), rhoB_log_grid(:,:)
        real(8), allocatable :: muA_fine_1d(:), muB_fine_1d(:)
        real(8), allocatable :: rhoA_log_fine(:,:), rhoB_log_fine(:,:)

        ! Step 1: Define finer grid
        n_muA = size(muA_values)
        n_muB = size(muB_values)
        grid_size = n_muA * multiplication

        allocate(muA_fine(grid_size), muB_fine(grid_size))
        call linspace(muA_values(1), muA_values(n_muA), muA_fine)
        call linspace(muB_values(1), muB_values(n_muB), muB_fine)

        ! Step 2: Log-transform the input grids
        allocate(rhoA_log_grid(n_muA, n_muB), rhoB_log_grid(n_muA, n_muB))
        rhoA_log_grid = log(rhoA_grid)
        rhoB_log_grid = log(rhoB_grid)

        ! Step 3: Interpolate onto the finer grid
        allocate(rhoA_log_fine(grid_size, grid_size), rhoB_log_fine(grid_size, grid_size))

        ! Loop over each point in the fine grid and interpolate
        do j = 1, grid_size
            do i = 1, grid_size
                rhoA_log_fine(i,j) = cubic_interp_2d(muA_values, muB_values, rhoA_log_grid, &
                                                    muA_fine(i), muB_fine(j))
                rhoB_log_fine(i,j) = cubic_interp_2d(muA_values, muB_values, rhoB_log_grid, &
                                                    muA_fine(i), muB_fine(j))
            end do
        end do

        ! Step 4: Exponentiate to get back to linear scale
        allocate(rhoA_fine(grid_size, grid_size), rhoB_fine(grid_size, grid_size))
        rhoA_fine = exp(rhoA_log_fine)
        rhoB_fine = exp(rhoB_log_fine)

        ! Cleanup
        deallocate(rhoA_log_grid, rhoB_log_grid, rhoA_log_fine, rhoB_log_fine)

        call write_grid2d(trim(output_dir)//"/interpolated_rhoA_grid.dat", rhoA_fine, "# Interpolated rhoA grid")
        call write_grid2d(trim(output_dir)//"/interpolated_rhoB_grid.dat", rhoB_fine, "# Interpolated rhoB grid")

    end subroutine interpolate_density_grids

    ! Helper: Linear spacing (like np.linspace)
    subroutine linspace(x_min, x_max, x_out)
        real(8), intent(in) :: x_min, x_max
        real(8), intent(out) :: x_out(:)
        integer :: i, n
        real(8) :: dx

        n = size(x_out)
        dx = (x_max - x_min) / (n - 1)
        do i = 1, n
            x_out(i) = x_min + (i-1) * dx
        end do
    end subroutine linspace

    ! Helper: 2D cubic interpolation (simplified)
    function cubic_interp_2d(x, y, z, x_target, y_target) result(z_interp)
        real(8), intent(in) :: x(:), y(:), z(:,:)
        real(8), intent(in) :: x_target, y_target
        real(8) :: z_interp

        ! For simplicity, use bilinear interpolation here.
        ! Replace with a proper 2D spline if available in your codebase.
        z_interp = bilinear_interp(x, y, z, x_target, y_target)
    end function cubic_interp_2d

    ! Helper: Bilinear interpolation (fallback)
    function bilinear_interp(x, y, z, x_target, y_target) result(z_interp)
        real(8), intent(in) :: x(:), y(:), z(:,:)
        real(8), intent(in) :: x_target, y_target
        real(8) :: z_interp

        ! Find the nearest indices (simplified)
        integer :: i, j, n_x, n_y
        real(8) :: x1, x2, y1, y2, z11, z12, z21, z22

        n_x = size(x)
        n_y = size(y)

        ! Find the interval for x
        do i = 1, n_x - 1
            if (x(i) <= x_target .and. x_target <= x(i+1)) exit
        end do
        x1 = x(i)
        x2 = x(i+1)

        ! Find the interval for y
        do j = 1, n_y - 1
            if (y(j) <= y_target .and. y_target <= y(j+1)) exit
        end do
        y1 = y(j)
        y2 = y(j+1)

        ! Get the 4 surrounding points
        z11 = z(i, j)
        z12 = z(i, j+1)
        z21 = z(i+1, j)
        z22 = z(i+1, j+1)

        ! Bilinear interpolation
        z_interp = (z11 * (x2 - x_target) * (y2 - y_target) + &
                    z21 * (x_target - x1) * (y2 - y_target) + &
                    z12 * (x2 - x_target) * (y_target - y1) + &
                    z22 * (x_target - x1) * (y_target - y1)) / &
                   ((x2 - x1) * (y2 - y1))
    end function bilinear_interp

end module density_interpolation