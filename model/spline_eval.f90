module spline_eval
  use spline_data
  implicit none
  public :: eval_spline, invert_spline, invert_spline_scan

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

  ! Inverse: given rho, return mu.
  ! Fast path (monotone splines): locate the interval by bisecting on the
  ! cached knot values y_left (no polynomial evaluation), then solve the
  ! cubic with safeguarded Newton (falls back to bisection if a Newton
  ! step would leave the bracket). ~5 evaluations per call.
  ! Falls back to the original full interval scan + bisection for any
  ! spline that isn't monotone at its knots.
  function invert_spline(y, spl) result(x)
    implicit none
    real(8),        intent(in) :: y
    type(spline_t), intent(in) :: spl
    real(8) :: x
    integer :: i

    if (.not. spl%monotone) then
      x = invert_spline_scan(y, spl)
      return
    end if

    i = find_interval_by_value(y, spl)
    x = invert_cubic_newton(y, spl, i)

  end function invert_spline

  ! Helper: find interval index for x using bisection search on x_left
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

  ! Helper: find interval index for target y using bisection search on
  ! the cached y_left knot values. Requires spl%monotone.
  function find_interval_by_value(y, spl) result(i)
    implicit none
    real(8),        intent(in) :: y
    type(spline_t), intent(in) :: spl
    integer :: i, lo, hi, mid
    logical :: increasing

    increasing = spl%y_left(spl%n) > spl%y_left(1)

    lo = 1
    hi = spl%n
    do while (hi - lo > 1)
      mid = (lo + hi) / 2
      if (increasing) then
        if (y >= spl%y_left(mid)) then
          lo = mid
        else
          hi = mid
        end if
      else
        if (y <= spl%y_left(mid)) then
          lo = mid
        else
          hi = mid
        end if
      end if
    end do
    i = lo

  end function find_interval_by_value

  ! Helper: solve c0 + c1*dx + c2*dx^2 + c3*dx^3 = y for dx in
  ! [0, x_right(i)-x_left(i)] via Newton's method safeguarded by
  ! bisection (bracket always maintained, so it cannot diverge or
  ! leave the interval).
  function invert_cubic_newton(y, spl, i) result(x)
    implicit none
    real(8),        intent(in) :: y
    type(spline_t), intent(in) :: spl
    integer,        intent(in) :: i
    real(8) :: x
    real(8) :: dxa, dxb, dx, dx_new, fa, f, fp
    real(8), parameter :: tol = 1d-12
    integer :: iter
    integer, parameter :: max_iter = 50

    dxa = 0d0
    dxb = spl%x_right(i) - spl%x_left(i)
    fa  = spl%y_left(i) - y

    dx = 0.5d0 * (dxa + dxb)

    do iter = 1, max_iter
      f  = spl%c0(i) + dx*(spl%c1(i) + dx*(spl%c2(i) + dx*spl%c3(i))) - y
      fp = spl%c1(i) + dx*(2d0*spl%c2(i) + 3d0*spl%c3(i)*dx)

      ! keep the bracket valid
      if (fa * f <= 0d0) then
        dxb = dx
      else
        dxa = dx
        fa  = f
      end if

      if (abs(f) < tol .or. (dxb - dxa) < tol) exit

      if (abs(fp) > 1d-14) then
        dx_new = dx - f / fp
      else
        dx_new = 0.5d0 * (dxa + dxb)
      end if

      ! reject Newton step if it leaves the bracket
      if (dx_new <= dxa .or. dx_new >= dxb) then
        dx_new = 0.5d0 * (dxa + dxb)
      end if

      dx = dx_new
    end do

    x = spl%x_left(i) + dx

  end function invert_cubic_newton

  ! Original method: full scan over all intervals to find a sign change,
  ! then plain bisection. Used automatically when a spline is not
  ! monotone at its knots. ~100 polynomial evaluations per call.
  function invert_spline_scan(y, spl) result(x)
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

  end function invert_spline_scan

end module spline_eval