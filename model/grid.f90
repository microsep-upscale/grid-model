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
    real(8) :: block_size, system_size, time_step, init_mu, area, Na, m_fluid
    real(8) :: left_mu, right_mu

    real(8), allocatable :: block_centers(:)
    real(8), allocatable :: chemical_potential(:)
    real(8), allocatable :: fluid_density(:)
    real(8), allocatable :: permeability(:)
    real(8), allocatable :: flux(:)
    real(8), allocatable :: grad_mu(:)

    real(8), allocatable :: mu_vs_rho(:)
    real(8), allocatable :: rho_vs_mu(:)
    real(8), allocatable :: mu_vs_k(:)
    real(8), allocatable :: k_vs_mu(:)

    real(8), parameter :: kcal_to_j = 4184.0d0

    ! ===============================
    ! Load coefficients
    ! ===============================
    call load_coeffs("../data/T300/mu_vs_rho_fit.txt", mu_vs_rho, deg1)
    call load_coeffs("../data/T300/rho_vs_mu_fit.txt", rho_vs_mu, deg2)
    call load_coeffs("../data/T300/mu_vs_K_fit.txt", mu_vs_k, deg3)
    call load_coeffs("../data/T300/K_vs_mu_fit.txt", k_vs_mu, deg4)

    ! ===============================
    ! System definition
    ! ===============================
    block_size = 5d-9 ! m
    system_size = 200d-9 ! m
    number_block = int(system_size / block_size)

    time_step = 1e-12 ! s
    area = 1e-18 ! m**2
    Na = 6.022e23 ! mol-1
    m_fluid = 10d-3  ! kg/mol

    ! Allocations
    allocate(block_centers(number_block))
    allocate(chemical_potential(number_block))
    allocate(fluid_density(number_block))
    allocate(permeability(number_block))
    allocate(flux(number_block))
    allocate(grad_mu(number_block))

    do i = 1, number_block
        block_centers(i) = (i-1) * block_size + block_size/2
    end do

    ! ===============================
    ! Initial conditions
    ! ===============================
    left_mu = -6d0 * kcal_to_j
    right_mu = -1d0 * kcal_to_j
    mu_mode = 2 ! 1=linear, 2=left, 3=right

    ! For sanity check, generate the rho/k vs mu curve that have been imported from MD
    ! call generate_tables(left_mu, right_mu, 10, &
    !                     rho_vs_mu, deg2, &
    !                     k_vs_mu, deg4, &
    !                     mu_vs_rho, deg1, &
    !                     mu_vs_k, deg3)

    ! Initialise the chemical potential profile within the pore
    call init_mu_profile(chemical_potential, block_centers, number_block, left_mu, right_mu, mu_mode)

    ! Initialise the density profile within the pore
    call init_rho_from_mu(fluid_density, chemical_potential, number_block, rho_vs_mu, deg2)

    flux = 0.0d0

    ! write(*,*) fluid_density

    ! ===============================
    ! Iteration
    ! ===============================
    n_iter = 1000
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
                        block_centers, chemical_potential, &
                        number_block, 1d9, "# x[nm] mu[J/mol]", .true.)

        call write_profile("output/rho_"//iter_str//".dat", &
                        block_centers, fluid_density, &
                        number_block, 1d9, "# x[nm] rho[kg/m3]", .false.)

        ! --- evaluate grad mu before any update ---
        do i = 2, number_block-1
            permeability(i) = poly_fit(chemical_potential(i), k_vs_mu, deg4) ! [mol/(J·m·s)]
            grad_mu(i) = (chemical_potential(i+1) - chemical_potential(i-1)) / (2 * block_size) ! [J/mol] / [m] = [J/(mol·m)]
            flux(i) = permeability(i) * grad_mu(i) ! [mol/(J·m·s)] * [J/(mol·m)] = [1/(m²·s)]
        end do

        ! --- update system ---
        do i = 2, number_block-1

            fluid_density(i) = fluid_density(i) + (flux(i-1) - flux(i+1)) * time_step / block_size ! 1/(m²·s) * s * / m = 1/m*3
                        
            chemical_potential(i) = poly_fit(fluid_density(i), mu_vs_rho, deg1) ! [kg/m³] -> [J/mol]

        end do

        ! log per iteration
        write(log_unit,*) "iter =", iter

    end do

    ! -------------------------
    ! write full profile file
    ! -------------------------
    write(iter_str,'(I4.4)') iter
    call write_profile("output/mu_"//iter_str//".dat", &
                    block_centers, chemical_potential, &
                    number_block, 1d9, "# x[nm] mu[J/mol]", .true.)

    call write_profile("output/rho_"//iter_str//".dat", &
                    block_centers, fluid_density, &
                    number_block, 1d9, "# x[nm] rho[kg/m3]", .false.)

    write(log_unit,*) "Simulation finished"
    close(log_unit)

end program grid_model