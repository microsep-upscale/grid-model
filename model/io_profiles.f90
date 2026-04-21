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

    subroutine write_profiles(iter, block_centers, chemical_potential, fluid_density, flux, number_block)
        implicit none
        integer, intent(in) :: iter
        integer, intent(in) :: number_block
        real(8), intent(in) :: block_centers(:), chemical_potential(:), fluid_density(:), flux(:)
        character(len=4) :: iter_str

        write(iter_str,'(I4.4)') iter

        call write_profile("output/mu_"//iter_str//".dat", &
                        block_centers, chemical_potential, &
                        number_block, 1d9, "# x[nm] mu[J/mol]", .true.)
        call write_profile("output/rho_"//iter_str//".dat", &
                        block_centers, fluid_density, &
                        number_block, 1d9, "# x[nm] rho[1/m3]", .false.)
        call write_profile("output/flux_"//iter_str//".dat", &
                        block_centers, flux, &
                        number_block, 1d9, "# x[nm] flux[1/s]", .false.)
    end subroutine write_profiles

end module io_profiles

