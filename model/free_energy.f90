module free_energy

    use io_profiles

    implicit none
    
    private
    public :: compute_free_energy_field

contains

    subroutine compute_free_energy_field(muA_fine, muB_fine, rhoA_fine, rhoB_fine, &
                                          p_fine, f_fine, output_dir)
        implicit none
        real(8), intent(in)  :: muA_fine(:), muB_fine(:)
        real(8), intent(in)  :: rhoA_fine(:,:), rhoB_fine(:,:)
        real(8), intent(in)  :: p_fine(:,:)
        real(8), allocatable, intent(out) :: f_fine(:,:)
        character(len=*), intent(in) :: output_dir

        integer :: i, j, grid_size
        character(len=64) :: size_str

        grid_size = size(muA_fine)

        allocate(f_fine(grid_size, grid_size))

        ! indexing="ij": MU_A_fine(i,j) = muA_fine(i), MU_B_fine(i,j) = muB_fine(j)
        do j = 1, grid_size
            do i = 1, grid_size
                f_fine(i,j) = muA_fine(i)*rhoA_fine(i,j) + muB_fine(j)*rhoB_fine(i,j) - p_fine(i,j)
            end do
        end do

        write(size_str,'(I0,A,I0)') grid_size, "x", grid_size
        call write_grid2d(trim(output_dir)//"/f_fine_"//trim(size_str)//".dat", &
                           f_fine, "# free energy field: f = muA*rhoA + muB*rhoB - p")

    end subroutine compute_free_energy_field

end module free_energy