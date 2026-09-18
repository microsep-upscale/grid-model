module run_binary_fluid_mod

    use table2d_data
    use table2d_io
    use table2d_eval
    use io_profiles
    use init_profiles
    use free_energy
    use timestep_control
    use convergence_control
    use numerical_gradient
    use pressure_integration
    use density_interpolation

    implicit none
    private
    public :: run_binary_fluid

contains

    subroutine run_binary_fluid(block_size_x, block_size_y, block_size_z, &
                                 system_size_x, time_step_in, max_time_step, &
                                 left_muA_in, right_muA_in, inside_muA_in, &
                                 left_muB_in, right_muB_in, inside_muB_in, mu_mode, &
                                 n_iter, n_jump, check_interval, conv_tol, &
                                 tol_min, tol_max, growth_factor, &
                                 rhoA_matrix, rhoB_matrix, muA_array, muB_array, &
                                 steady_tol, n_steady, output_dir)

        implicit none

        ! ---- arguments (values as read from the namelist) ----
        real(8), intent(in) :: block_size_x, block_size_y, block_size_z
        real(8), intent(in) :: system_size_x
        real(8), intent(in) :: time_step_in, max_time_step
        real(8), intent(in) :: left_muA_in, right_muA_in, inside_muA_in
        real(8), intent(in) :: left_muB_in, right_muB_in, inside_muB_in
        integer, intent(in) :: mu_mode
        integer(kind=8), intent(in) :: n_iter
        integer, intent(in) :: n_jump, check_interval
        real(8), intent(in) :: conv_tol, tol_min, tol_max, growth_factor
        character(len=*), intent(in) :: rhoA_matrix, rhoB_matrix, muA_array, muB_array
        real(8), intent(in) :: steady_tol
        integer, intent(in) :: n_steady
        character(len=*), intent(in) :: output_dir

        integer :: i, j

        ! Inside run_binary_fluid:
        integer, parameter :: multiplication = 5
        real(8), allocatable :: muA_fine(:), muB_fine(:)
        real(8), allocatable :: rhoA_fine(:,:), rhoB_fine(:,:)
        real(8), allocatable :: muA_values(:), muB_values(:)
        real(8), allocatable :: rhoA_grid(:,:), rhoB_grid(:,:)
        real(8), allocatable :: p_fine(:,:), p_alt(:,:), p_avg(:,:)
        real(8), allocatable :: f_fine(:,:)
        real(8), allocatable :: rhoA_log_fine(:,:), rhoB_log_fine(:,:)
        real(8), allocatable :: dlogrhoA_dmuA(:,:), dlogrhoA_dmuB(:,:)
        real(8), allocatable :: dlogrhoB_dmuA(:,:), dlogrhoB_dmuB(:,:)
        real(8), allocatable :: drhoA_dmuA(:,:), drhoA_dmuB(:,:)
        real(8), allocatable :: drhoB_dmuA(:,:), drhoB_dmuB(:,:)

        ! Read muA_values, muB_values, rhoA_grid, rhoB_grid from files
        call read_2d_grid(rhoA_matrix, rhoA_grid)
        call read_2d_grid(rhoB_matrix, rhoB_grid)
        call read_1d_array(muA_array, muA_values)
        call read_1d_array(muB_array, muB_values)

        ! Print rhoA_grid (9x9)
        if (.not. allocated(rhoA_grid) .or. size(rhoA_grid) == 0) then
            write(*,*) "Error: rhoA_grid is empty"
            stop 1
        end if

        ! Print rhoB_grid (9x9)
        if (.not. allocated(rhoB_grid) .or. size(rhoB_grid) == 0) then
            write(*,*) "Error: rhoB_grid is empty"
            stop 1
        end if

        ! Interpolate
        call interpolate_density_grids(muA_values, muB_values, rhoA_grid, rhoB_grid, &
                                    multiplication, muA_fine, muB_fine, rhoA_fine, &
                                    rhoB_fine, output_dir)

        allocate(rhoA_log_fine, source=log(rhoA_fine))
        allocate(rhoB_log_fine, source=log(rhoB_fine))
        allocate(dlogrhoA_dmuA(size(muA_fine), size(muB_fine)))
        allocate(dlogrhoA_dmuB(size(muA_fine), size(muB_fine)))
        allocate(dlogrhoB_dmuA(size(muA_fine), size(muB_fine)))
        allocate(dlogrhoB_dmuB(size(muA_fine), size(muB_fine)))

        call gradient_2d(rhoA_log_fine, muA_fine, muB_fine, dlogrhoA_dmuA, dlogrhoA_dmuB)
        call gradient_2d(rhoB_log_fine, muA_fine, muB_fine, dlogrhoB_dmuA, dlogrhoB_dmuB)

        allocate(drhoA_dmuA, source=rhoA_fine*dlogrhoA_dmuA)
        allocate(drhoA_dmuB, source=rhoA_fine*dlogrhoA_dmuB)
        allocate(drhoB_dmuA, source=rhoB_fine*dlogrhoB_dmuA)
        allocate(drhoB_dmuB, source=rhoB_fine*dlogrhoB_dmuB)

        call write_grid2d(trim(output_dir)//"/drhoA_dmuA_grid.dat", drhoA_dmuA, "# Calculated drhoA_dmuA")
        call write_grid2d(trim(output_dir)//"/drhoB_dmuA_grid.dat", drhoB_dmuA, "# Calculated drhoB_dmuA")
        call write_grid2d(trim(output_dir)//"/drhoA_dmuB_grid.dat", drhoA_dmuB, "# Calculated drhoA_dmuB")
        call write_grid2d(trim(output_dir)//"/drhoB_dmuB_grid.dat", drhoB_dmuB, "# Calculated drhoB_dmuB")

        call compute_pressure_fields(muA_fine, muB_fine, rhoA_fine, rhoB_fine, &
            p_fine, p_alt, p_avg, output_dir)

        call compute_free_energy_field(muA_fine, muB_fine, rhoA_fine, rhoB_fine, &
                                            p_fine, f_fine, output_dir)

    end subroutine run_binary_fluid

end module run_binary_fluid_mod