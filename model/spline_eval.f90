module spline_eval
  use spline_data
  implicit none
  public :: eval_spline, invert_spline

contains

  ! Forward: given mu, return rho or M
  function eval_spline(x, spl) result(y)
    implicit none
    real(8),        intent(in) :: x
    type(spline_t), intent(in) :: spl
    real(8) :: y, dx
    integer :: i

    ! Bisection search for the interval
    i = find_interval(x, spl)
    dx = x - spl%x_left(i)
    y  = spl%c0(i) + dx*(spl%c1(i) + dx*(spl%c2(i) + dx*spl%c3(i)))

  end function eval_spline

  ! Inverse: given rho, return mu (bisection on spline)
  function invert_spline(y, spl) result(x)
    implicit none
    real(8),        intent(in) :: y
    type(spline_t), intent(in) :: spl
    real(8) :: x
    real(8) :: xa, xb, fa, fb, xm, fm
    real(8), parameter :: tol = 1d-10
    integer :: i, iter

    ! Find the interval where the sign changes
    x = spl%x_left(1)  ! default fallback
    do i = 1, spl%n
      xa = spl%x_left(i)
      xb = spl%x_right(i)
      fa = eval_spline(xa, spl) - y
      fb = eval_spline(xb, spl) - y
      if (fa * fb <= 0d0) then
        ! Brentq-style bisection within this interval
        do iter = 1, 100
          xm = (xa + xb) / 2d0
          fm = eval_spline(xm, spl) - y
          if (abs(fm) < tol .or. (xb - xa) < tol) exit
          if (fa * fm < 0d0) then
            xb = xm; fb = fm
          else
            xa = xm; fa = fm
          end if
        end do
        x = xm
        return
      end if
    end do

  end function invert_spline

  ! Helper: find interval index for x using bisection search
  function find_interval(x, spl) result(i)
    implicit none
    real(8),        intent(in) :: x
    type(spline_t), intent(in) :: spl
    integer :: i, lo, hi, mid

    lo = 1
    hi = spl%n
    do while (hi - lo > 1)
      mid = (lo + hi) / 2
      if (x >= spl%x_left(mid)) then
        lo = mid
      else
        hi = mid
      end if
    end do
    i = lo

  end function find_interval

end module spline_eval