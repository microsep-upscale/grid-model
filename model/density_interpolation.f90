module density_interpolation
    use spline_data
    use spline_eval
    use io_profiles

    implicit none

    private
    public :: interpolate_density_grids, bicubic_eval_point

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
        integer :: n_muA, n_muB, grid_size
        real(8), allocatable :: rhoA_log_grid(:,:), rhoB_log_grid(:,:)
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

        ! Step 3: Bicubic-spline interpolation onto the finer grid
        allocate(rhoA_log_fine(grid_size, grid_size), rhoB_log_fine(grid_size, grid_size))

        call bicubic_interp_2d(muA_values, muB_values, rhoA_log_grid, muA_fine, muB_fine, rhoA_log_fine)
        call bicubic_interp_2d(muA_values, muB_values, rhoB_log_grid, muA_fine, muB_fine, rhoB_log_fine)

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

    ! ------------------------------------------------------------------
    ! Bicubic interpolation via "spline of splines":
    !   1) For each row i, precompute a natural cubic spline in y
    !      (i.e. through z(i,:) as a function of y).
    !   2) For each target y_target, evaluate all n_muA row-splines at
    !      y_target to get a 1D profile tmp(1:n_muA) as a function of x.
    !   3) Spline tmp(:) in x and evaluate at each x_target.
    ! This gives a C2-smooth surface in x and in y (no bilinear kinks).
    ! ------------------------------------------------------------------
    subroutine bicubic_interp_2d(x, y, z, x_fine, y_fine, z_fine)
        implicit none
        real(8), intent(in)  :: x(:), y(:), z(:,:)
        real(8), intent(in)  :: x_fine(:), y_fine(:)
        real(8), intent(out) :: z_fine(:,:)

        integer :: n_x, n_y, nxf, nyf, i, j
        real(8), allocatable :: y2_rows(:,:)   ! spline coeffs of each row in y, (n_x, n_y)
        real(8), allocatable :: tmp(:), tmp_y2(:)

        n_x  = size(x)
        n_y  = size(y)
        nxf  = size(x_fine)
        nyf  = size(y_fine)

        allocate(y2_rows(n_x, n_y))
        do i = 1, n_x
            call spline_coeffs(y, z(i,:), n_y, y2_rows(i,:))
        end do

        allocate(tmp(n_x), tmp_y2(n_x))

        do j = 1, nyf
            ! Evaluate each row's y-spline at y_fine(j) to build the x-profile
            do i = 1, n_x
                call spline_eval1(y, z(i,:), y2_rows(i,:), n_y, y_fine(j), tmp(i))
            end do
            ! Spline that x-profile and evaluate at every x_fine value
            call spline_coeffs(x, tmp, n_x, tmp_y2)
            do i = 1, nxf
                call spline_eval1(x, tmp, tmp_y2, n_x, x_fine(i), z_fine(i,j))
            end do
        end do

        deallocate(y2_rows, tmp, tmp_y2)
    end subroutine bicubic_interp_2d

    ! ------------------------------------------------------------------
    ! Natural cubic spline: compute second derivatives y2(1:n)
    ! (standard tridiagonal algorithm, y2(1) = y2(n) = 0)
    ! ------------------------------------------------------------------
    subroutine spline_coeffs(x, y, n, y2)
        implicit none
        integer, intent(in)  :: n
        real(8), intent(in)  :: x(n), y(n)
        real(8), intent(out) :: y2(n)
        real(8), allocatable :: u(:)
        real(8) :: sig, p
        integer :: i, k

        allocate(u(n))
        y2(1) = 0d0
        u(1)  = 0d0

        do i = 2, n-1
            sig   = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
            p     = sig * y2(i-1) + 2d0
            y2(i) = (sig - 1d0) / p
            u(i)  = (6d0 * ((y(i+1)-y(i))/(x(i+1)-x(i)) - (y(i)-y(i-1))/(x(i)-x(i-1))) &
                     / (x(i+1)-x(i-1)) - sig*u(i-1)) / p
        end do

        y2(n) = 0d0
        do k = n-1, 1, -1
            y2(k) = y2(k)*y2(k+1) + u(k)
        end do

        deallocate(u)
    end subroutine spline_coeffs

    ! ------------------------------------------------------------------
    ! Evaluate a natural cubic spline (given precomputed y2) at xq.
    ! Uses binary search to locate the interval, then clamps xq to the
    ! table range to avoid extrapolation past the edges.
    ! ------------------------------------------------------------------
    subroutine spline_eval1(xa, ya, y2a, n, xq_in, yq)
        implicit none
        integer, intent(in)  :: n
        real(8), intent(in)  :: xa(n), ya(n), y2a(n)
        real(8), intent(in)  :: xq_in
        real(8), intent(out) :: yq
        integer :: klo, khi, k
        real(8) :: h, a, b, xq

        xq = min(max(xq_in, xa(1)), xa(n))  ! clamp to avoid extrapolation

        klo = 1
        khi = n
        do while (khi - klo > 1)
            k = (khi + klo) / 2
            if (xa(k) > xq) then
                khi = k
            else
                klo = k
            end if
        end do

        h = xa(khi) - xa(klo)
        if (h <= 0d0) then
            write(*,*) "Error: spline_eval1 - duplicate x values in table"
            stop 1
        end if

        a = (xa(khi) - xq) / h
        b = (xq - xa(klo)) / h

        yq = a*ya(klo) + b*ya(khi) + &
             ((a**3 - a)*y2a(klo) + (b**3 - b)*y2a(khi)) * (h**2) / 6d0
    end subroutine spline_eval1

    ! Evaluate the bicubic-spline surface at a single arbitrary point.
    ! x,y: coarse grid (muA_values/muB_values); z: coarse grid values (e.g. log(rho)).
    function bicubic_eval_point(x, y, z, xq, yq) result(zq)
        implicit none
        real(8), intent(in) :: x(:), y(:), z(:,:)
        real(8), intent(in) :: xq, yq
        real(8) :: zq

        integer :: n_x, n_y, i
        real(8), allocatable :: tmp(:), tmp_y2(:), row_y2(:)

        n_x = size(x)
        n_y = size(y)

        allocate(tmp(n_x), tmp_y2(n_x), row_y2(n_y))

        do i = 1, n_x
            call spline_coeffs(y, z(i,:), n_y, row_y2)
            call spline_eval1(y, z(i,:), row_y2, n_y, yq, tmp(i))
        end do

        call spline_coeffs(x, tmp, n_x, tmp_y2)
        call spline_eval1(x, tmp, tmp_y2, n_x, xq, zq)

        deallocate(tmp, tmp_y2, row_y2)
    end function bicubic_eval_point

end module density_interpolation