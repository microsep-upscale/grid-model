module run_single_fluid_mod

    use spline_data
    use spline_io
    use spline_eval
    use io_profiles
    use init_profiles
    use timestep_control
    use convergence_control

    implicit none
    private
    public :: run_single_fluid

contains

    subroutine run_single_fluid(block_size_x, block_size_y, block_size_z, &
                                 system_size_x, time_step_in, max_time_step, &
                                 left_mu_in, right_mu_in, inside_mu_in, mu_mode, &
                                 n_iter, n_jump, check_interval, conv_tol, &
                                 tol_min, tol_max, growth_factor, &
                                 rho_spline_file, M_spline_file, &
                                 steady_tol, n_steady, output_dir)

        implicit none

        ! ---- arguments (values as read from the namelist) ----
        real(8), intent(in) :: block_size_x, block_size_y, block_size_z
        real(8), intent(in) :: system_size_x
        real(8), intent(in) :: time_step_in, max_time_step
        real(8), intent(in) :: left_mu_in, right_mu_in, inside_mu_in
        integer, intent(in) :: mu_mode
        integer(kind=8), intent(in) :: n_iter
        integer, intent(in) :: n_jump, check_interval
        real(8), intent(in) :: conv_tol, tol_min, tol_max, growth_factor
        character(len=*), intent(in) :: rho_spline_file, M_spline_file
        real(8), intent(in) :: steady_tol
        integer, intent(in) :: n_steady
        character(len=*), intent(in) :: output_dir

        ! ---- local state (was module/program-level before) ----
        type(spline_t) :: spl_rho, spl_M

        integer(kind=8) :: iter
        integer :: number_block, number_edge, i, block, edge, edge1, edge2
        integer :: log_unit, conv_unit, steady_count

        real(8) :: block_area_yz, block_volume_xyz
        real(8) :: time_step
        real(8) :: left_mu, right_mu, inside_mu
        real(8) :: density_edge, permeability_edge
        real(8) :: net_flux, force_edge, flux_edge
        real(8) :: flux_mean, flux_std, flux_conservation
        real(8) :: time, prev_mean_density
        logical :: converged, use_steady_state

        real(8), allocatable :: block_centers(:)
        real(8), allocatable :: block_edges(:)
        real(8), allocatable :: chemical_potential(:)
        real(8), allocatable :: fluid_density(:)
        real(8), allocatable :: delta_density(:)
        real(8), allocatable :: permeability(:)
        real(8), allocatable :: flux_edges(:)
        real(8), allocatable :: grad_mu(:)

        real(8), parameter :: Na = 6.022e23             ! mol-1
        real(8), parameter :: kcal_to_j = 4184.0d0       ! J/kcal
        real(8), parameter :: mu_equal_tol = 1d-6        ! J/mol

        ! ---- unit conversion (unchanged from original) ----
        left_mu   = left_mu_in   * kcal_to_j   ! J/mol
        right_mu  = right_mu_in  * kcal_to_j   ! J/mol
        inside_mu = inside_mu_in * kcal_to_j   ! J/mol
        time_step = time_step_in

        ! Decide once which stopping criterion applies: symmetric reservoirs
        ! give flux_mean -> 0, so flux-conservation can never converge there.
        use_steady_state = (abs(right_mu - left_mu) < mu_equal_tol)

        ! Load splines
        call load_spline(trim(rho_spline_file), spl_rho)
        call load_spline(trim(M_spline_file), spl_M)

        ! System definition
        number_block = nint(system_size_x / block_size_x) ! nint, not int to avoid truncating

        number_edge = number_block - 1                  ! number of edge = number of block + 1 - number of reservoirs
        block_area_yz = block_size_y*block_size_z        ! m^2
        block_volume_xyz = block_size_x*block_area_yz     ! m^3

        ! Quantities that are defined within cell
        allocate(block_centers(number_block))            ! m
        allocate(chemical_potential(number_block))       ! J/mol
        allocate(fluid_density(number_block))             ! m^-3
        allocate(delta_density(number_block))              ! m^-3
        allocate(permeability(number_block))                ! s/kg

        ! Quantities that are defined between cells (at the edges)
        allocate(block_edges(number_edge))               ! m
        allocate(grad_mu(number_edge))                   ! J/(mol·m)
        allocate(flux_edges(number_edge))                ! s^-1

        ! Vector "position" along the pore
        do block = 1, number_block
            block_centers(block) = (block-1) * block_size_x + block_size_x/2 ! in m
        end do
        do edge = 1, number_edge
            block_edges(edge) = (edge-1) * block_size_x ! in m
        end do

        ! Pick the initial chemical potential profile
        ! 1: linear increase from left to right
        ! 2: constant value inside the pore
        call init_mu_profile(chemical_potential, number_block, left_mu, right_mu, inside_mu, mu_mode)

        ! Initialise the density profile within the pore
        call compute_rho_from_mu(fluid_density, chemical_potential, number_block, spl_rho)

        ! Initialise to 0
        flux_edges = 0.0d0
        permeability = 0.0d0
        grad_mu = 0.0d0
        delta_density = 0.0d0

        ! For output file
        log_unit = 99
        conv_unit = 97

        open(newunit=log_unit, file=trim(output_dir)//"/grid.log", status="replace", action="write")
        write(log_unit,*) "Simulation started"
        write(log_unit,*) "Number of iterations =", n_iter

        open(newunit=conv_unit, file=trim(output_dir)//"/conservation.dat", status="replace", action="write")
        if (use_steady_state) then
            write(conv_unit,*) "# iter    density_drift    time_step[s]"
        else
            write(conv_unit,*) "# iter    flux_conservation    flux_mean[1/s]    flux_std[1/s]    time_step[s]"
        end if

        iter = 0
        time = 0d0

        steady_count = 0
        prev_mean_density = -1d0
        flux_mean = 0d0

        do while (iter < n_iter)

            iter = iter + 1
            time = time + time_step         ! s

            ! outputs
            if (mod(iter, n_jump) == 0) then
                ! write profiles to files
                call write_profiles(int(iter/n_jump), time, &
                                    block_centers, block_edges, &
                                    chemical_potential, fluid_density, &
                                    permeability, flux_edges, grad_mu, &
                                    number_block, number_edge, output_dir)

                ! log per iteration
                write(log_unit,*) "iter =", iter
                ! terminal progress
                write(*,'(A,I10,A,ES10.3,A,ES10.3)') &
                    "  iter=", iter, &
                    "  dt=", time_step, &
                    "  flux_mean=", flux_mean
            end if

            ! evaluate grad mu *before* any update of the chemical potential
            do block = 2, number_block
                edge = block-1
                grad_mu(edge) = (chemical_potential(block) - chemical_potential(block-1)) / block_size_x ! J/(mol·m)
            end do

            ! evaluate permeability before any update
            ! Important note: the value of permeability for block=1 and block=n is wrong,
            ! in practice, it should be calculated from a theory accounting for entrance effects
            do block = 1, number_block
                permeability(block) = eval_spline(chemical_potential(block), spl_M)
            end do

            do edge = 1, number_edge
                block = edge + 1
                density_edge = (fluid_density(block) + fluid_density(block-1))/2        ! m^-3
                permeability_edge = (permeability(block) + permeability(block-1))/2     ! s/kg
                force_edge = -grad_mu(edge) / Na                                        ! N (per-particle force)
                flux_edge = permeability_edge * force_edge * density_edge * block_area_yz ! s^-1
                flux_edges(edge) = flux_edge
            end do

            ! update density with safety check
            do block = 2, number_block-1

                edge1 = block-1
                edge2 = block

                ! Net flux in the block is the sum over the two edges
                net_flux = flux_edges(edge1) - flux_edges(edge2)                            ! s^-1

                ! Updated density based on the flux
                delta_density(block) = net_flux * time_step / block_volume_xyz              ! m^-3
                fluid_density(block) = fluid_density(block) + delta_density(block)          ! m^-3

                ! Update the chemical potential
                chemical_potential(block) = invert_spline(fluid_density(block), spl_rho)    ! J/mol

            end do

            ! check for NaN values
            if (mod(iter, check_interval) == 0) then

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
                        write(*,*) "  time_step      =", time_step
                        write(*,*)
                        write(log_unit,*) "NaN detected at iter =", iter, " block =", i, " — stopping"
                        stop "NaN detected — simulation aborted"
                    end if
                end do

                ! check for negative density
                do i = 2, number_block-1
                    if (fluid_density(i) <= 0d0) then
                        write(*,*) "Negative density at iter =", iter, " block =", i
                        write(*,*) "  delta_density  =", delta_density(i)
                        write(*,*) "  time_step      =", time_step
                        write(log_unit,*) "Negative density at iter =", iter, " block =", i, " — stopping"
                        stop "Negative density — simulation aborted"
                    end if
                end do

                ! Adapt timestep every "check_interval"
                call adapt_timestep(time_step, delta_density, fluid_density, number_block, &
                                    tol_min, tol_max, growth_factor, max_time_step)

            end if

            ! every check_interval iterations, check for convergence
            ! (criterion depends on mu_mode: mode 1 is driven, flux-based;
            ! mode 2 is symmetric, flux_mean -> 0, so density-based instead)
            if (mod(iter, check_interval) == 0) then

                ! Criterion depends on whether the reservoirs drive a net
                ! flux: if they're equal, flux_mean -> 0 by symmetry and the
                ! flux-conservation ratio can never converge, so density
                ! steady-state is used instead.
                if (use_steady_state) then
                    call check_steady_state(iter, fluid_density, &
                                             number_block, steady_tol, n_steady, &
                                             time_step, conv_unit, &
                                             steady_count, prev_mean_density, converged)
                else
                    call check_flux_conservation(iter, check_interval, time_step, &
                                                  flux_edges, number_block, conv_tol, &
                                                  conv_unit, flux_mean, flux_std, &
                                                  flux_conservation, converged)
                end if

                if (converged) then
                    write(*,*) "Converged at iter =", iter
                    write(log_unit,*) "Converged at iter =", iter

                    call write_profiles(int(iter/n_jump), time, &
                                        block_centers, block_edges, &
                                        chemical_potential, fluid_density, &
                                        permeability, flux_edges, grad_mu, &
                                        number_block, number_edge, output_dir, label="final")

                    close(log_unit)
                    close(conv_unit)
                    return

                end if

           end if

        end do

        close(log_unit)
        close(conv_unit)

    end subroutine run_single_fluid

end module run_single_fluid_mod