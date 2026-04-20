program grid_model

    use coeff_io

    implicit none

    integer :: number_block, dt, i, degree, gradient
    integer :: deg1, deg2, deg3, deg4

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
    print *, "mu =", mu

    rho_rec = poly_fit(mu, maniac_rho_mu_coeffs, deg2)
    print *, "rho_rec =", rho_rec

    K = poly_fit(mu, nemd_k_mu_coeffs, deg3)
    print *, "K =", K*1d6

    mu_rec = poly_fit(K, nemd_mu_k_coeffs, deg4)
    print *, "mu_rec =", mu_rec

    rho_last = poly_fit(mu_rec, maniac_rho_mu_coeffs, deg2)
    print *, "rho_last =", rho_last

    ! ===============================
    ! System definition
    ! ===============================
    block_size = 5d-9
    system_size = 200d-9
    number_block = system_size / block_size

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
    gradient = 1
    init_mu = left_mu

    allocate(chemical_potential(number_block))
    allocate(fluid_density(number_block))

    open(10, file="output/initial_pore.dat")

    do i = 2, number_block-1

        if (gradient == 1) then
            chemical_potential(i) = left_mu + (right_mu-left_mu) * &
                                    (i-1) / real(number_block-1)

        else
            chemical_potential(i) = init_mu
        end if

        fluid_density(i) = poly_fit(chemical_potential(i), maniac_rho_mu_coeffs, deg2)

        write(10,*) x_axis(i)*1d9, chemical_potential(i), fluid_density(i)

    end do

    ! ===============================
    ! Allocations
    ! ===============================
    allocate(permeability(number_block))
    allocate(flux(number_block))
    allocate(new_chemical_potential(number_block))
    allocate(new_fluid_density(number_block))

    flux(1) = 0.0d0

    chemical_potential(1) = left_mu
    chemical_potential(number_block) = right_mu

    fluid_density(1) = poly_fit(chemical_potential(1), maniac_rho_mu_coeffs, deg2)
    fluid_density(number_block) = poly_fit(chemical_potential(number_block), maniac_rho_mu_coeffs, deg2)

    ! ===============================
    ! First iteration
    ! ===============================
    open(20, file="output/first_iteration.dat")

    do i = 2, number_block-1

        permeability(i) = poly_fit(chemical_potential(i), nemd_k_mu_coeffs, deg3)

        flux(i) = permeability(i) * &
                (chemical_potential(i+1)-chemical_potential(i)) * &
                (1d0/block_size) * (-1d0)

        new_fluid_density(i) = fluid_density(i) + &
            (flux(i-1)-flux(i)) * time_step * dt

        new_chemical_potential(i) = poly_fit(new_fluid_density(i), maniac_mu_rho_coeffs, deg1)

        write(20,*) x_axis(i)*1d9, new_chemical_potential(i), new_fluid_density(i)

        chemical_potential(i) = new_chemical_potential(i)
        fluid_density(i) = new_fluid_density(i)

    end do

    ! ===============================
    ! cleanup
    ! ===============================
    deallocate(x_axis, chemical_potential, fluid_density)
    deallocate(permeability, flux, new_chemical_potential, new_fluid_density)

contains

    function poly_fit(x, coeffs, degree) result(y)
        implicit none
        real(8), intent(in) :: x
        real(8), intent(in) :: coeffs(:)
        integer, intent(in) :: degree
        real(8) :: y
        integer :: i

        y = 0d0
        do i = 0, degree
            y = y + coeffs(i+1) * x**(degree-i)
        end do
    end function

end program grid_model