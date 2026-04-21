module spline_data
  implicit none
  public

  type spline_t
    integer :: n                        ! number of intervals
    real(8), allocatable :: x_left(:)  ! left knot of each interval
    real(8), allocatable :: x_right(:) ! right knot
    real(8), allocatable :: c0(:)      ! constant term
    real(8), allocatable :: c1(:)      ! linear term
    real(8), allocatable :: c2(:)      ! quadratic term
    real(8), allocatable :: c3(:)      ! cubic term
  end type spline_t

end module spline_data
