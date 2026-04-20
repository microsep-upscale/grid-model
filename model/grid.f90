program grid_model

    use coeff_io
    use poly_fit_mod
    use io_profiles
    use init_profiles

    implicit none

    integer :: number_block, dt, i, degree, gradient
    integer :: deg1, deg2, deg3, deg4
    integer :: n_iter, iter
    integer :: log_unit
    character(len=256) :: filename
    integer :: data_unit
    character(len=4) :: iter_str
    integer :: mu_mode

    real(8) :: rho, mu, rho_rec, K, mu_rec, rho_last
    real(8) :: block_size, system_size, time_step, init_mu
    real(8) :: left_mu, right_mu

    real(8), allocatable :: x_axis(:)
    real(8), allocatable :: chemical_potential(:)
    real(8), allocatable :: fluid_density(:)
    real(8), allocatable :: permeability(:)
    real(8), allocatable :: flux(:)
    real(8), allocatable :: new_chemical_potential(:)
    real(8), allocatable :: new_fluid_density(:)

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
    ! Test section
    ! ===============================
    rho = 200.0d0
    print *, "rho =", rho

    mu = poly_fit(rho, maniac_mu_rho_coeffs, deg1)
    rho_rec = poly_fit(mu, maniac_rho_mu_coeffs, deg2)
    K = poly_fit(mu, nemd_k_mu_coeffs, deg3)
    mu_rec = poly_fit(K, nemd_mu_k_coeffs, deg4)
    rho_last = poly_fit(mu_rec, maniac_rho_mu_coeffs, deg2)

    ! ===============================
    ! System definition
    ! ===============================
    block_size = 5d-9
    system_size = 200d-9
    number_block = int(system_size / block_size)

    time_step = 100d-15
    dt = 100000000

    allocate(x_axis(number_block))

    do i = 1, number_block
        x_axis(i) = (i-1) * block_size
    end do

    ! ===============================
    ! Initial conditions
    ! ===============================
    left_mu = 2d0 * 4184d0
    right_mu = 3d0 * 4184d0
    mu_mode = 1 ! 1=linear, 2=left, 3=right

    allocate(chemical_potential(number_block))
    allocate(fluid_density(number_block))

    call init_mu_profile(chemical_potential, x_axis, number_block, &
                        left_mu, right_mu, mu_mode)

    call init_rho_from_mu(fluid_density, chemical_potential, &
                        number_block, maniac_rho_mu_coeffs, deg2)

    ! ===============================
    ! Allocations
    ! ===============================
    allocate(permeability(number_block))
    allocate(flux(number_block))
    allocate(new_chemical_potential(number_block))
    allocate(new_fluid_density(number_block))

    do i = 1, number_block
        new_chemical_potential(i) = chemical_potential(i)
        new_fluid_density(i) = fluid_density(i)
    end do

    flux(1) = 0.0d0

    chemical_potential(1) = left_mu
    chemical_potential(number_block) = right_mu

    fluid_density(1) = poly_fit(chemical_potential(1), maniac_rho_mu_coeffs, deg2)
    fluid_density(number_block) = poly_fit(chemical_potential(number_block), maniac_rho_mu_coeffs, deg2)

    ! ===============================
    ! Iteration
    ! ===============================
    n_iter = 2
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
                        x_axis, new_chemical_potential, &
                        number_block, 1d9, &
                        "# x[nm] mu[J/mol]")

        call write_profile("output/rho_"//iter_str//".dat", &
                        x_axis, new_fluid_density, &
                        number_block, 1d9, &
                        "# x[nm] rho[kg/m3]")

        ! --- update system ---
        do i = 2, number_block-1

            permeability(i) = poly_fit(new_chemical_potential(i), nemd_k_mu_coeffs, deg3)

            flux(i) = permeability(i) * &
                    (new_chemical_potential(i+1) - new_chemical_potential(i)) * &
                    (1d0 / block_size) * (-1d0)

            new_fluid_density(i) = new_fluid_density(i) + &
                (flux(i-1) - flux(i)) * time_step * dt * iter

            new_chemical_potential(i) = poly_fit(new_fluid_density(i), maniac_mu_rho_coeffs, deg1)

        end do

        ! log per iteration
        write(log_unit,*) "iter =", iter

    end do

    ! -------------------------
    ! write full profile file
    ! -------------------------
    write(iter_str,'(I4.4)') iter
    call write_profile("output/mu_"//iter_str//".dat", &
                    x_axis, new_chemical_potential, &
                    number_block, 1d9, &
                    "# x[nm] mu[J/mol]")

    call write_profile("output/rho_"//iter_str//".dat", &
                    x_axis, new_fluid_density, &
                    number_block, 1d9, &
                    "# x[nm] rho[kg/m3]")

    write(log_unit,*) "Simulation finished"
    close(log_unit)

    ! ===============================
    ! cleanup
    ! ===============================
    deallocate(x_axis, chemical_potential, fluid_density)
    deallocate(permeability, flux, new_chemical_potential, new_fluid_density)

end program grid_model