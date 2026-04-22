module io_profiles
  implicit none
contains

  subroutine write_profile(filename, x_axis, field, n, scale_x, header)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: header
    integer,          intent(in) :: n
    real(8),          intent(in) :: x_axis(n)
    real(8),          intent(in) :: field(n)
    real(8),          intent(in) :: scale_x
    integer :: i, unit

    open(newunit=unit, file=filename, status="replace", action="write")
    write(unit,*) header
    do i = 1, n
      write(unit,'(2ES20.10)') x_axis(i)*scale_x, field(i)
    end do
    close(unit)
  end subroutine write_profile

    subroutine write_profiles(iter, time, block_centers, block_edges, chemical_potential, fluid_density, &
                            permeability, flux_edges, grad_mu, number_block, number_edge, label)
        implicit none
        integer, intent(in)                    :: iter, number_block, number_edge
        real(8), intent(in)                    :: time
        real(8), intent(in)                    :: block_centers(number_block)
        real(8), intent(in)                    :: block_edges(number_edge)
        real(8), intent(in)                    :: chemical_potential(number_block)
        real(8), intent(in)                    :: fluid_density(number_block)
        real(8), intent(in)                    :: permeability(number_block)
        real(8), intent(in)                    :: flux_edges(number_edge)
        real(8), intent(in)                    :: grad_mu(number_edge)
        character(len=*), intent(in), optional :: label

        character(len=32)  :: iter_str
        character(len=128) :: meta

        if (present(label)) then
            iter_str = trim(label)
        else
            write(iter_str,'(I8.8)') iter
        end if

        write(meta,'(A,I10,A,ES14.6)') "# iter= ", iter, "  time[s]= ", time

        call write_profile("output/mu_"//trim(iter_str)//".dat", &
                            block_centers, chemical_potential, &
                            number_block, 1d9, trim(meta)//"  x[nm]  mu[J/mol]")

        call write_profile("output/rho_"//trim(iter_str)//".dat", &
                            block_centers, fluid_density, &
                            number_block, 1d9, trim(meta)//"  x[nm]  rho[1/m3]")

        call write_profile("output/perm_"//trim(iter_str)//".dat", &
                            block_centers, permeability, &
                            number_block, 1d9, trim(meta)//"  x[nm]  M[s/kg/m]")

        call write_profile("output/flux_"//trim(iter_str)//".dat", &
                            block_edges, flux_edges, &
                            number_edge, 1d9, trim(meta)//"  x[nm]  flux[1/s]")

        call write_profile("output/grad_mu_"//trim(iter_str)//".dat", &
                            block_edges, grad_mu, &
                            number_edge, 1d9, trim(meta)//"  x[nm]  grad_mu[J/mol/m]")

    end subroutine write_profiles

end module io_profiles