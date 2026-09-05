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
    real(8) :: dx
    logical :: increasing
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
    allocate(spl%y_left(spl%n), spl%y_right(spl%n))
 
    do i = 1, spl%n
      read(unit, *) spl%x_left(i), spl%x_right(i), &
                    spl%c0(i), spl%c1(i), spl%c2(i), spl%c3(i)
    end do

    close(unit)

    ! Cache y at each knot so invert_spline can locate the interval
    ! without evaluating the cubic. y_left(i) = poly at dx=0 = c0(i).
    do i = 1, spl%n
      spl%y_left(i) = spl%c0(i)
      dx = spl%x_right(i) - spl%x_left(i)
      spl%y_right(i) = spl%c0(i) + dx*(spl%c1(i) + dx*(spl%c2(i) + dx*spl%c3(i)))
    end do

    ! Determine monotonicity from y_left alone (see note: y_right(i) and
    ! y_left(i+1) can disagree in the last bits and must not be mixed
    ! into this test, or the fast path silently never triggers)
    spl%monotone = .true.
    if (spl%n >= 2) then
      increasing = spl%y_left(spl%n) > spl%y_left(1)
      do i = 2, spl%n
        if (increasing) then
          if (spl%y_left(i) <= spl%y_left(i-1)) then
            spl%monotone = .false.
            exit
          end if
        else
          if (spl%y_left(i) >= spl%y_left(i-1)) then
            spl%monotone = .false.
            exit
          end if
        end if
      end do
    end if

  end subroutine load_spline

end module spline_io