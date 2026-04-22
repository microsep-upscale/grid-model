program grid_model

    use spline_data
    use spline_io
    use spline_eval
    use coeff_io
    use spline_io
    use poly_fit_mod
    use io_profiles
    use init_profiles
    use tables_io

    implicit none

    type(spline_t) :: spl_rho, spl_M

    integer(kind=8) :: n_iter, iter
    integer :: number_block, number_edge, i, block, edge, edge1, edge2
    integer :: log_unit, data_unit, mu_mode, n_jump, check_interval, conv_unit
    real(8) :: flux, block_size, system_size, time_step, block_area, block_volume, Na, m_fluid, max_time_step
    real(8) :: left_mu, right_mu, density_edge, permeability_edge
    real(8) :: max_rel_change, tol_min, tol_max, growth_factor, net_flux
    real(8) :: force_edge, flux_edge, flux_mean, flux_std, flux_conservation
    real(8) :: time, conv_tol

    real(8), allocatable :: block_centers(:)
    real(8), allocatable :: block_edges(:)
    real(8), allocatable :: chemical_potential(:)
    real(8), allocatable :: fluid_density(:)
    real(8), allocatable :: delta_density(:)
    real(8), allocatable :: permeability(:)
    real(8), allocatable :: flux_edges(:)
    real(8), allocatable :: grad_mu(:)

    real(8), parameter :: kcal_to_j = 4184.0d0

    ! ! Load coefficients
    ! call load_coeffs("../data/T84/rho_vs_mu_fit.txt", rho_vs_mu, deg2)
    ! call load_coeffs("../data/T84/M_vs_mu_fit.txt", M_vs_mu, deg4)

    ! Load splines (replaces load_coeffs)
    call load_spline("../data/T300/spline_rho_vs_mu.txt", spl_rho)
    call load_spline("../data/T300/spline_M_vs_mu.txt",   spl_M)

    ! For simulation stoping
    conv_tol = 1e-3

    ! System definition
    block_size = 2d-9 ! m
    system_size = 200d-9 ! m
    number_block = int(system_size / block_size)
    number_edge = number_block - 1 ! number of edge = number of block  + 1 - number of reservoirs
    block_area = block_size*block_size ! m**2
    block_volume = block_size*block_size*block_size ! m**3

    ! Parameters
    time_step = 1e-15 ! s (initial timestep)
    max_time_step = 1e-12 ! s (max timestep)
    Na = 6.022e23 ! mol-1
    m_fluid = 40d-3  ! kg/mol

    ! Quantities that are define within cell
    allocate(block_centers(number_block))
    allocate(chemical_potential(number_block)) ! J/mol
    allocate(fluid_density(number_block)) ! m-3
    allocate(delta_density(number_block)) ! m-3
    allocate(permeability(number_block)) ! s/kg/m

    ! Quantities that are defined between cells (at the edges)
    allocate(block_edges(number_edge))
    allocate(grad_mu(number_edge)) ! J/(mol·m)
    allocate(flux_edges(number_edge)) ! s-1

    ! Vector "position" along the pore
    do block = 1, number_block
        block_centers(block) = (block-1) * block_size + block_size/2 ! in m
    end do
    do edge = 1, number_edge
        block_edges(edge) = (edge-1) * block_size ! in m
    end do

    ! Initial conditions (make sure this corresponds to the range simulated in MD)
    left_mu = -3.0d0 * kcal_to_j ! J/mol
    right_mu = -2.0d0 * kcal_to_j ! J/mol

    mu_mode = 1 ! Pick the initial chemical potential profile
    ! 1: linear increase from left to right
    ! 2: use the value from the left reservoir
    ! 3: use the value from the right reservoir

    ! ! For sanity check, generate the curves that have been imported from MD
    ! call generate_tables(left_mu, right_mu, 10, rho_vs_mu, deg2, M_vs_mu, deg4)

    ! Initialise the chemical potential profile within the pore
    call init_mu_profile(chemical_potential, number_block, left_mu, right_mu, mu_mode)

    ! Initialise the density profile within the pore
    ! call compute_rho_from_mu(fluid_density, chemical_potential, number_block, rho_vs_mu, deg2)
    call compute_rho_from_mu(fluid_density, chemical_potential, number_block, spl_rho)

    ! Initialise to 0
    flux_edges = 0.0d0
    permeability = 0.0d0
    grad_mu = 0.0d0
    delta_density = 0.0d0

    ! Iteration
    n_iter = 500000000_8 ! simulation max duration (note that the simulation stop if convergence reached before)
    n_jump = 50000 ! interval for data printing
    check_interval = 10000 ! interval for timestep reevaluation

    ! For output file
    log_unit = 99
    data_unit = 98
    conv_unit = 97

    ! Timestep reevaluation
    tol_min=1e-7
    tol_max=1e-5
    growth_factor=1.1

    open(newunit=log_unit, file="output/grid.log", status="replace", action="write")
    write(log_unit,*) "Simulation started"
    write(log_unit,*) "Number of iterations =", n_iter

    open(newunit=conv_unit, file="output/conservation.dat", status="replace", action="write")
    write(conv_unit,*) "# iter    flux_conservation    flux_mean[1/s]    flux_std[1/s]    time_step[s]"

    iter = 0
    time = 0d0 ! in seconde

    do while (iter < n_iter)

        iter = iter + 1
        time = time + time_step ! in second

        ! outputs
        if (mod(iter, n_jump) == 0) then
            ! write profiles to files
            call write_profiles(int(iter/n_jump), time, &
                                block_centers, block_edges, &
                                chemical_potential, fluid_density, &
                                permeability, flux_edges, grad_mu, &
                                number_block, number_edge)

            ! log per iteration
            write(log_unit,*) "iter =", iter
            ! terminal progress
            write(*,'(A,I10,A,F6.1,A,ES10.3,A,ES10.3)') &
                "  iter=", iter, &
                "  progress=", 100.0*int(iter/n_iter), &
                "  dt=", time_step, &
                "  flux_mean=", flux_mean
        end if

        ! evaluate grad mu *before* any update of the chemical potential
        do block = 2, number_block
            edge = block-1
            grad_mu(edge) = (chemical_potential(block) - chemical_potential(block-1)) / block_size ! [J/(mol·m)]
        end do

        ! evaluate permeability before any update
        ! Important note: the value of permeability for block=1 and block=n is wrong,
        ! in practice, it should be calculated from a theory accounting for entrance effects
        ! This will be done someday

        ! Polynomial
        ! do block = 1, number_block
        !     permeability(block) = poly_fit(chemical_potential(block), M_vs_mu, deg4) ! [s/kg/m]
        ! end do
        ! Spline
        do block = 1, number_block
            permeability(block) = eval_spline(chemical_potential(block), spl_M)
        end do

        ! update density with safety check
        do edge = 1, number_edge
            block = edge + 1 
            density_edge = (fluid_density(block) + fluid_density(block-1))/2 ! m-3
            permeability_edge = (permeability(block) + permeability(block-1))/2 ! m-3
            force_edge = -grad_mu(edge) * density_edge * block_volume / Na ! N
            flux_edge = permeability_edge * force_edge ! 1/s
            flux_edges(edge) = flux_edge
        end do

        ! update density with safety check
        do block = 2, number_block-1

            edge1 = block-1
            edge2 = block

            ! Net flux in the block is the sum over the two edges
            net_flux = flux_edges(edge1) - flux_edges(edge2)

            ! Updated density based on the flux
            delta_density(block) = net_flux * time_step / block_volume
            fluid_density(block) = fluid_density(block) + delta_density(block)
            
            ! Update the chemical potential
            ! chemical_potential(block) = invert_poly2(fluid_density(block), rho_vs_mu)
            chemical_potential(block) = invert_spline(fluid_density(block), spl_rho)

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

        ! every check_interval iterations, check flux conservation
        if (mod(iter, check_interval) == 0) then

            ! flux conservation score: std dev of interior edge fluxes
            ! at steady state all interior fluxes should be equal
            flux_mean = 0d0
            do block = 2, number_block-1
                edge = block - 1
                flux_mean = flux_mean + flux_edges(edge)
            end do
            flux_mean = flux_mean / real(number_block-2, 8)

            flux_std = 0d0
            do block = 2, number_block-1
                edge = block - 1
                flux_std = flux_std + (flux_edges(edge) - flux_mean)**2
            end do
            flux_std = sqrt(flux_std / real(number_block-2, 8))

            ! normalized score: 0 = perfectly conserved, 1 = large variation
            if (abs(flux_mean) > 0d0) then
                flux_conservation = flux_std / abs(flux_mean)
            else
                flux_conservation = 0d0
            end if

            if (flux_conservation > 10d0) then
                write(*,*) "time_step =", time_step
                write(*,*) flux_conservation, flux_std, flux_mean
                write(*,*) flux_edges
                stop "Conservation score diverged — simulation aborted"
            end if

            ! check for convergence
            if (flux_conservation < conv_tol .and. iter > check_interval) then
                write(*,*) "Converged at iter =", iter
                write(*,*) "  flux_conservation =", flux_conservation
                write(*,*) "  flux_mean =", flux_mean
                write(log_unit,*) "Converged at iter =", iter, " flux_conservation =", flux_conservation

                call write_profiles(int(iter/n_jump), time, &
                                    block_centers, block_edges, &
                                    chemical_potential, fluid_density, &
                                    permeability, flux_edges, grad_mu, &
                                    number_block, number_edge, label="final")

                write(conv_unit,'(I10,4ES20.10)') iter, flux_conservation, flux_mean, flux_std, time_step
                stop "The simulation has converged successfully"
            end if

            write(conv_unit,'(I10,4ES20.10)') iter, flux_conservation, flux_mean, flux_std, time_step

       end if

    end do


end program grid_model
