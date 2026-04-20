module poly_fit_mod
  implicit none
  private
  public :: poly_fit

contains

  function poly_fit(x, coeffs, degree) result(y)
    implicit none

    real(8), intent(in) :: x
    real(8), intent(in) :: coeffs(:)
    integer, intent(in) :: degree
    real(8) :: y
    integer :: i

    y = 0d0

    do i = 0, degree
       y = y + coeffs(i+1) * x**(degree - i)
    end do

  end function poly_fit

end module poly_fit_mod
