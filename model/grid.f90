program grid_model

    use coeff_io
    use poly_fit_mod
    use io_profiles
    use init_profiles
    use tables_io

    implicit none

    integer :: number_block, number_edge, i, degree, gradient, block, edge
    integer :: deg1, deg2, deg3, deg4
    integer :: n_iter, iter
    integer :: log_unit
    character(len=256) :: filename
    integer :: data_unit
    character(len=4) :: iter_str
    integer :: mu_mode
    integer :: n_jump, check_interval

    real(8) :: rho, mu, rho_rec, K, mu_rec, flux
    real(8) :: block_size, system_size, time_step, init_mu, block_area, block_volume, Na, m_fluid, max_time_step
    real(8) :: left_mu, right_mu
    real(8) :: max_rel_change, tol_min, tol_max, growth_factor, grad_mu_edge, force_edge, net_flux
    real(8) :: force_right, force_left, flux_left, flux_right

    real(8), allocatable :: block_centers(:)
    real(8), allocatable :: block_edges(:)
    real(8), allocatable :: chemical_potential(:)
    real(8), allocatable :: fluid_density(:)
    real(8), allocatable :: delta_density(:)
    real(8), allocatable :: permeability(:)
    real(8), allocatable :: flux_edges(:)
    real(8), allocatable :: grad_mu(:)
    real(8), allocatable :: total_force(:)

    real(8), allocatable :: rho_vs_mu(:)
    real(8), allocatable :: M_vs_mu(:)

    real(8), parameter :: kcal_to_j = 4184.0d0

    ! Load coefficients
    call load_coeffs("../data/T84/rho_vs_mu_fit.txt", rho_vs_mu, deg2)
    call load_coeffs("../data/T84/M_vs_mu_fit.txt", M_vs_mu, deg4)

    ! System definition
    block_size = 5d-9 ! m
    system_size = 200d-9 ! m
    number_block = int(system_size / block_size)
    number_edge = number_block - 1 ! number of edge = number of block  + 1 - number of reservoirs
    block_area = block_size*block_size ! m**2
    block_volume = block_size*block_size*block_size ! m**3

    ! Parameters
    time_step = 1e-15 ! s (initial timestep)
    max_time_step = 1e-9 ! s (max timestep)
    Na = 6.022e23 ! mol-1
    m_fluid = 40d-3  ! kg/mol

    ! Quantities that are define within cell
    allocate(block_centers(number_block))
    allocate(chemical_potential(number_block)) ! J/mol
    allocate(fluid_density(number_block)) ! m-3
    allocate(delta_density(number_block)) ! m-3
    allocate(permeability(number_block)) ! s/kg/m
    allocate(total_force(number_block)) ! N

    ! Quantities that are defined between cells (at the edges)
    allocate(block_edges(number_edge))
    allocate(grad_mu(number_edge)) ! J/(mol·m)
    allocate(flux_edges(number_edge)) ! s-1

    block_edges(1) = 0
    do block = 1, number_block
        block_centers(block) = (block-1) * block_size + block_size/2
        block_edges(block+1) = block * block_size
    end do

    ! Initial conditions (make sure this corresponds to the range simulated in MD)
    left_mu = -3.0d0 * kcal_to_j ! J/mol
    right_mu = -2.0d0 * kcal_to_j ! J/mol

    mu_mode = 1 ! Pick the initial chemical potential profile
    ! 1: linear increase from left to right
    ! 2: use the value from the left reservoir
    ! 3: use the value from the right reservoir

    ! For sanity check, generate the curves that have been imported from MD
    call generate_tables(left_mu, right_mu, 10, rho_vs_mu, deg2, M_vs_mu, deg4)

    ! Initialise the chemical potential profile within the pore
    call init_mu_profile(chemical_potential, block_centers, number_block, left_mu, right_mu, mu_mode)

    ! Initialise the density profile within the pore
    call compute_rho_from_mu(fluid_density, chemical_potential, number_block, rho_vs_mu, deg2)

    ! Initialise to 0
    flux_edges = 0.0d0
    permeability = 0.0d0
    grad_mu = 0.0d0
    delta_density = 0.0d0

    ! ===============================
    ! Iteration
    ! ===============================
    n_iter = 10000
    n_jump = 10
    check_interval = 10
    log_unit = 99
    data_unit = 98

    ! timestep reevaluation
    tol_min=1e-7
    tol_max=1e-6
    growth_factor=2.0

    open(newunit=log_unit, file="output/grid.log", status="replace", action="write")
    write(log_unit,*) "Simulation started"
    write(log_unit,*) "Number of iterations =", n_iter

    iter = 0

    do while (iter < n_iter) ! .and. time_step <= max_time_step) (iter < n_iter .and. time_step <= max_time_step)

        iter = iter + 1

        if (mod(iter, n_jump) == 0) then
            ! write profiles to files
            call write_profiles(int(iter/n_jump), block_centers, chemical_potential, fluid_density, flux_edges, number_block)
            ! log per iteration
            write(log_unit,*) "iter =", iter
        end if

        ! evaluate grad mu before any update
        do block = 2, number_block
            edge = block-1
            ! gradient at the edge (difference between neighboring block centers)
            grad_mu(edge) = (chemical_potential(block) - chemical_potential(block-1)) / block_size ! [J/(mol·m)]
        end do

        ! evaluate permeability before any update
        do block = 2, number_block-1
            permeability(block) = poly_fit(chemical_potential(block), M_vs_mu, deg4) ! [s/kg/m]
        end do

        ! update density with safety check
        do block = 2, number_block-1
            edge = block-1
            force_left = -grad_mu(edge) * fluid_density(block) * block_volume / Na ! N
            force_right = -grad_mu(edge+1) * fluid_density(block) * block_volume / Na ! N
            flux_left = permeability(block) * force_left ! 1/s
            flux_right = permeability(block) * force_right ! 1/s
            net_flux = flux_left - flux_right ! 1/s
            delta_density(block) = net_flux * time_step / block_volume
            fluid_density(block) = fluid_density(block) + delta_density(block)
            chemical_potential(block) = invert_poly2(fluid_density(block), rho_vs_mu)

            flux_edges(edge) = net_flux

            ! guard against out-of-range density before inversion
            if (fluid_density(block) <= 0d0) then
                write(*,*) "Negative density at iter =", iter, " block =", i
                exit
            end if

        end do

        ! check for NaN values
        do i = 2, number_block-1
            if (isnan(fluid_density(i)) .or. isnan(chemical_potential(i))) then
                write(*,*) "NaN detected at iter =", iter, " block =", i
                write(*,*)
                write(*,*) "  delta_density  =", delta_density
                write(*,*)
                write(*,*) "  fluid_density  =", fluid_density
                write(*,*)
                write(*,*) "  chemical_potential =", chemical_potential
                write(*,*)
                write(*,*) "  flux           =", flux
                write(*,*)
                write(*,*) "  time_step      =", time_step
                write(*,*)
                write(log_unit,*) "NaN detected at iter =", iter, " block =", i, " — stopping"
                stop "NaN detected — simulation aborted"
            end if
        end do

        ! every check_interval iterations, adapt timestep
        if (mod(iter, check_interval) == 0) then

            ! evaluate the max relative change in density
            max_rel_change = 0d0
            do i = 2, number_block-1
                max_rel_change = max(max_rel_change, abs(delta_density(i) / fluid_density(i)))
            end do

            ! adjust timestep to make sure the typical change in density is within desired windows
            if (max_rel_change < tol_min) then
                time_step = time_step * growth_factor
                if (time_step > max_time_step) then ! do not let timestep diverge
                    time_step = max_time_step
                end if
            else if (max_rel_change > tol_max) then
                time_step = time_step / growth_factor
            end if
        end if

    end do

    ! write final profiles to files
    call write_profiles(int(iter/n_jump), block_centers, chemical_potential, fluid_density, flux_edges, number_block)

end program grid_model