module coeff_io
  implicit none
contains

  subroutine load_coeffs(filename, coeffs, degree)
    implicit none

    character(len=*), intent(in) :: filename
    real(8), allocatable, intent(out) :: coeffs(:)
    integer, intent(out) :: degree

    integer :: n, i, unit

    open(newunit=unit, file=filename, status='old', action='read')

    read(unit,*) n
    allocate(coeffs(n))

    do i = 1, n
       read(unit,*) coeffs(i)
    end do

    close(unit)

    degree = n - 1
  end subroutine

end module coeff_io