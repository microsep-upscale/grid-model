module io_profiles
  implicit none
contains

  subroutine write_profile(filename, x_axis, field, number_block, scale_x, header, edge)
    implicit none

    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: header
    integer, intent(in) :: number_block
    real(8), intent(in) :: x_axis(:)
    real(8), intent(in) :: field(:)
    real(8), intent(in) :: scale_x
    logical, intent(in) :: edge
    
    integer :: i, unit

    open(newunit=unit, file=filename, status="replace", action="write")

    write(unit,*) header

    if (edge) then
        do i = 1, number_block
            write(unit,'(2ES20.10)') x_axis(i)*scale_x, field(i)
        end do
    else
        do i = 2, number_block-1
            write(unit,'(2ES20.10)') x_axis(i)*scale_x, field(i)
        end do
    end if

    close(unit)

  end subroutine write_profile

end module io_profiles