module numerical_gradient
    implicit none
    private
    public :: gradient_2d

contains

    ! 1D non-uniform gradient, 2nd-order accurate everywhere (matches 
    ! the python implementatoin of np.gradient with edge_order=2)
    subroutine gradient_1d(f, x, n, df)
        implicit none
        integer, intent(in)  :: n
        real(8), intent(in)  :: f(n), x(n)
        real(8), intent(out) :: df(n)
        integer :: i
        real(8) :: dx1, dx2, a, b, c

        if (n < 3) then
            write(*,*) "Error: gradient_1d requires at least 3 points for edge_order=2"
            stop 1
        end if

        ! Left edge (one-sided, 2nd order)
        dx1 = x(2) - x(1)
        dx2 = x(3) - x(2)
        a = -(2d0*dx1 + dx2) / (dx1*(dx1+dx2))
        b = (dx1 + dx2) / (dx1*dx2)
        c = -dx1 / (dx2*(dx1+dx2))
        df(1) = a*f(1) + b*f(2) + c*f(3)

        ! Interior (central, non-uniform)
        do i = 2, n-1
            dx1 = x(i)   - x(i-1)
            dx2 = x(i+1) - x(i)
            a = -dx2 / (dx1*(dx1+dx2))
            b = (dx2 - dx1) / (dx1*dx2)
            c = dx1 / (dx2*(dx1+dx2))
            df(i) = a*f(i-1) + b*f(i) + c*f(i+1)
        end do

        ! Right edge (one-sided, 2nd order)
        dx1 = x(n-1) - x(n-2)
        dx2 = x(n)   - x(n-1)
        a = dx2 / (dx1*(dx1+dx2))
        b = -(dx1 + dx2) / (dx1*dx2)
        c = (2d0*dx2 + dx1) / (dx2*(dx1+dx2))
        df(n) = a*f(n-2) + b*f(n-1) + c*f(n)
    end subroutine gradient_1d

    ! 2D gradient: dfdx = d f/d x (axis 1), dfdy = d f/d y (axis 2)
    subroutine gradient_2d(f, x, y, dfdx, dfdy)
        implicit none
        real(8), intent(in)  :: f(:,:), x(:), y(:)
        real(8), intent(out) :: dfdx(:,:), dfdy(:,:)
        integer :: i, j, nx, ny

        nx = size(f,1)
        ny = size(f,2)

        do j = 1, ny
            call gradient_1d(f(:,j), x, nx, dfdx(:,j))
        end do

        do i = 1, nx
            call gradient_1d(f(i,:), y, ny, dfdy(i,:))
        end do
    end subroutine gradient_2d

end module numerical_gradient