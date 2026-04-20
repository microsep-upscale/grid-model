module io_profiles
  implicit none
contains

  subroutine write_profile(filename, x_axis, field, number_block, scale_x, header)
    implicit none

    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: header
    integer, intent(in) :: number_block
    real(8), intent(in) :: x_axis(:)
    real(8), intent(in) :: field(:)
    real(8), intent(in) :: scale_x

    integer :: i, unit

    open(newunit=unit, file=filename, status="replace", action="write")

    write(unit,*) header

    do i = 1, number_block
        write(unit,'(2ES20.10)') x_axis(i)*scale_x, field(i)
    end do

    close(unit)

  end subroutine write_profile

end module io_profiles