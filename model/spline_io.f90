module spline_io
  use spline_data
  implicit none
  public :: load_spline

contains

  subroutine load_spline(filename, spl)
    implicit none
    character(len=*), intent(in)  :: filename
    type(spline_t),   intent(out) :: spl
    integer :: unit, i
    character(len=256) :: line  ! for skipping header lines

    open(newunit=unit, file=filename, status='old', action='read')

    ! Skip 2 header lines (written by numpy savetxt with header=...)
    read(unit, '(A)') line
    read(unit, '(A)') line

    ! Count intervals — read all lines, rewind, allocate
    spl%n = 0
    do
      read(unit, *, end=10) line
      spl%n = spl%n + 1
    end do
    10 continue
    rewind(unit)
    read(unit, '(A)') line  ! skip headers again
    read(unit, '(A)') line

    allocate(spl%x_left(spl%n), spl%x_right(spl%n))
    allocate(spl%c0(spl%n), spl%c1(spl%n), spl%c2(spl%n), spl%c3(spl%n))

    do i = 1, spl%n
      read(unit, *) spl%x_left(i), spl%x_right(i), &
                    spl%c0(i), spl%c1(i), spl%c2(i), spl%c3(i)
    end do

    close(unit)
  end subroutine load_spline

end module spline_io