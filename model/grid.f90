program grid_model

    use coeff_io
    use poly_fit_mod
    use io_profiles
    use init_profiles
    use tables_io

    implicit none

    integer :: number_block, i, degree, gradient
    integer :: deg1, deg2, deg3, deg4
    integer :: n_iter, iter
    integer :: log_unit
    character(len=256) :: filename
    integer :: data_unit
    character(len=4) :: iter_str
    integer :: mu_mode

    real(8) :: rho, mu, rho_rec, K, mu_rec
    real(8) :: block_size, system_size, time_step, init_mu
    real(8) :: left_mu, right_mu

    real(8), allocatable :: x_axis(:)
    real(8), allocatable :: chemical_potential(:)
    real(8), allocatable :: fluid_density(:)
    real(8), allocatable :: permeability(:)
    real(8), allocatable :: flux(:)
    real(8), allocatable :: grad_mu(:)

    real(8), allocatable :: maniac_mu_rho_coeffs(:)
    real(8), allocatable :: maniac_rho_mu_coeffs(:)
    real(8), allocatable :: nemd_k_mu_coeffs(:)
    real(8), allocatable :: nemd_mu_k_coeffs(:)

    ! ===============================
    ! Load coefficients
    ! ===============================
    call load_coeffs("input/maniac_mu_rho_coeffs.dat", maniac_mu_rho_coeffs, deg1)
    call load_coeffs("input/maniac_rho_mu_coeffs.dat", maniac_rho_mu_coeffs, deg2)
    call load_coeffs("input/nemd_k_mu_coeffs.dat", nemd_k_mu_coeffs, deg3)
    call load_coeffs("input/nemd_mu_k_coeffs.dat", nemd_mu_k_coeffs, deg4)

    ! ===============================
    ! System definition
    ! ===============================
    block_size = 20d-9
    system_size = 200d-9
    number_block = int(system_size / block_size)

    time_step = 1

    ! Allocations
    allocate(x_axis(number_block))
    allocate(chemical_potential(number_block))
    allocate(fluid_density(number_block))
    allocate(permeability(number_block))
    allocate(flux(number_block))
    allocate(grad_mu(number_block))

    do i = 1, number_block
        x_axis(i) = (i-1) * block_size
    end do

    ! ===============================
    ! Initial conditions
    ! ===============================
    left_mu = 2d0 * 4184d0
    right_mu = 3d0 * 4184d0
    mu_mode = 2 ! 1=linear, 2=left, 3=right

    ! For sanity check, generate the rho/k vs mu curve that have been imported from MD
    ! call generate_tables(left_mu, right_mu, 200, maniac_rho_mu_coeffs, deg2, nemd_k_mu_coeffs, deg3)

    ! Initialise the chemical potential profile within the pore
    call init_mu_profile(chemical_potential, x_axis, number_block, left_mu, right_mu, mu_mode)

    ! Initialise the density profile within the pore
    call init_rho_from_mu(fluid_density, chemical_potential, number_block, maniac_rho_mu_coeffs, deg2)

    flux(1) = 0.0d0

    fluid_density(1) = poly_fit(chemical_potential(1), maniac_rho_mu_coeffs, deg2)
    fluid_density(number_block) = poly_fit(chemical_potential(number_block), maniac_rho_mu_coeffs, deg2)

    ! ===============================
    ! Iteration
    ! ===============================
    n_iter = 10
    log_unit = 99
    data_unit = 98

    open(newunit=log_unit, file="output/grid.log", status="replace", action="write")
    write(log_unit,*) "Simulation started"
    write(log_unit,*) "Number of iterations =", n_iter

    do iter = 1, n_iter

        ! -------------------------
        ! write full profile file
        ! -------------------------
        write(iter_str,'(I4.4)') iter
        call write_profile("output/mu_"//iter_str//".dat", &
                        x_axis, chemical_potential, &
                        number_block, 1d9, "# x[nm] mu[J/mol]", .true.)

        call write_profile("output/rho_"//iter_str//".dat", &
                        x_axis, fluid_density, &
                        number_block, 1d9, "# x[nm] rho[kg/m3]", .false.)

        ! --- evaluate grad mu before any update ---
        do i = 2, number_block-1
            grad_mu(i) = (chemical_potential(i+1) - chemical_potential(i-1)) / (2 * block_size)
        end do

        ! --- update system ---
        do i = 2, number_block-1

            permeability(i) = poly_fit(chemical_potential(i), nemd_k_mu_coeffs, deg3)

            flux(i) = permeability(i) * grad_mu(i)

            fluid_density(i) = fluid_density(i) + (flux(i-1) - flux(i+1)) * time_step

            chemical_potential(i) = poly_fit(fluid_density(i), maniac_mu_rho_coeffs, deg1)

        end do

        ! log per iteration
        write(log_unit,*) "iter =", iter

    end do

    ! -------------------------
    ! write full profile file
    ! -------------------------
    write(iter_str,'(I4.4)') iter
    call write_profile("output/mu_"//iter_str//".dat", &
                    x_axis, chemical_potential, &
                    number_block, 1d9, "# x[nm] mu[J/mol]", .true.)

    call write_profile("output/rho_"//iter_str//".dat", &
                    x_axis, fluid_density, &
                    number_block, 1d9, "# x[nm] rho[kg/m3]", .false.)

    write(log_unit,*) "Simulation finished"
    close(log_unit)

    ! ===============================
    ! cleanup
    ! ===============================
    deallocate(x_axis, chemical_potential, fluid_density, permeability, flux)

end program grid_model