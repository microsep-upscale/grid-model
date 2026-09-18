program grid_model

    use spline_data
    use spline_io
    use spline_eval
    use coeff_io
    use poly_fit_mod
    use io_profiles
    use init_profiles
    use tables_io
    use timestep_control
    use convergence_control
    use run_single_fluid_mod
    use run_binary_fluid_mod

    implicit none

    integer(kind=8) :: n_iter
    integer :: mu_mode
    integer ::n_jump
    integer :: input_unit
    integer :: ios
    integer :: check_interval
    integer :: fluid_mode
    integer :: n_steady

    real(8) :: block_size_x, block_size_y, block_size_z
    real(8) :: system_size_x

    real(8) :: time_step, max_time_step
    real(8) :: tol_min, tol_max, growth_factor
    real(8) :: conv_tol, steady_tol
    character(len=256) :: output_dir

    ! --- single-fluid boundary conditions ---
    real(8) :: left_mu, right_mu, inside_mu
    character(len=256) :: rho_spline_file, M_spline_file

    ! --- binary-fluid boundary conditions ---
    real(8) :: left_muA, right_muA, inside_muA
    real(8) :: left_muB, right_muB, inside_muB
    character(len=256) :: rhoA_matrix_path, rhoB_matrix_path, muA_array_path, muB_array_path

    ! Declare a namelist and give the variables defaults
    namelist /params/ fluid_mode, &
                    block_size_x, block_size_y, block_size_z, system_size_x, &
                    time_step, max_time_step, mu_mode, n_iter, &
                    n_jump, check_interval, conv_tol, tol_min, tol_max, &
                    growth_factor, steady_tol, n_steady, output_dir, &
                    left_mu, right_mu, inside_mu, rho_spline_file, M_spline_file, &
                    left_muA, right_muA, inside_muA, &
                    left_muB, right_muB, inside_muB, &
                    rhoA_matrix_path, rhoB_matrix_path, muA_array_path, muB_array_path

    ! Set defaults (fallback if grid.in doesn't define them)

    ! Choose between single-phase fluid and two-phase fluid 
    fluid_mode = 1

    ! --- Geometry ---
    block_size_x = 1d-9         ! m
    block_size_y = 1d-9         ! m
    block_size_z = 1d-9         ! m
    system_size_x = 100d-9      ! m

    ! --- Time integration ---
    time_step = 1e-15           ! s
    max_time_step = 1e-12       ! s

    ! --- Boundary conditions ---
    left_mu = -9.5d0            ! kcal/mol
    inside_mu = -11d0           ! kcal/mol
    right_mu = -13.0d0          ! kcal/mol
    mu_mode = 1

    ! --- Run control ---
    n_iter = 500000000_8
    n_jump = 500000
    check_interval = 10000
    conv_tol = 1e-3
    tol_min = 1e-7
    tol_max = 1e-5
    growth_factor = 1.1
    steady_tol = 1d-6
    n_steady = 5

    ! --- Output ---
    output_dir = "output"

    ! --- Single-fluid defaults ---
    left_mu = -9.5d0
    inside_mu = -11d0
    right_mu = -13.0d0
    mu_mode = 1
    rho_spline_file = "../data/single-phase/lj-T300/h1.0/spline_rho_vs_mu.txt"
    M_spline_file   = "../data/single-phase/lj-T300/h1.0/spline_M_vs_mu.txt"

    ! --- Binary-fluid defaults ---
    left_muA   = -9.5d0
    inside_muA = -11d0
    right_muA  = -13.0d0
    left_muB   = -9.5d0
    inside_muB = -11d0
    right_muB  = -13.0d0
    rhoA_matrix_path = "../data/two-phase/lj-slit/rhoA_9x9.dat"
    rhoB_matrix_path = "../data/two-phase/lj-slit/rhoB_9x9.dat"
    muA_array_path = "../data/two-phase/lj-slit/muA_unique.dat"
    muB_array_path = "../data/two-phase/lj-slit/muA_unique.dat"

    open(newunit=input_unit, file="grid.in", status="old", action="read")
    read(input_unit, nml=params, iostat=ios)
    if (ios /= 0) then
        write(*,*) "Namelist read failed, iostat=", ios
        stop 1
    end if
    close(input_unit)
    call execute_command_line("mkdir -p " // trim(output_dir))

    select case (fluid_mode)
    case (1)
        write(*,*) "Starting single fluid"
        call run_single_fluid(block_size_x, block_size_y, block_size_z, &
                               system_size_x, time_step, max_time_step, &
                               left_mu, right_mu, inside_mu, mu_mode, &
                               n_iter, n_jump, check_interval, conv_tol, &
                               tol_min, tol_max, growth_factor, &
                               rho_spline_file, M_spline_file, &
                               steady_tol, n_steady, output_dir)
    case (2)
        write(*,*) "Starting binary fluid"
        call run_binary_fluid(block_size_x, block_size_y, block_size_z, &
                               system_size_x, time_step, max_time_step, &
                               left_muA, right_muA, inside_muA, &
                               left_muB, right_muB, inside_muB, mu_mode, &
                               n_iter, n_jump, check_interval, conv_tol, &
                               tol_min, tol_max, growth_factor, &
                               rhoA_matrix_path, rhoB_matrix_path, muA_array_path, muB_array_path, &
                               steady_tol, n_steady, output_dir)
    case default
        write(*,*) "Unknown fluid_mode =", fluid_mode, " (must be 1 or 2)"
        stop 1
    end select

end program grid_model
