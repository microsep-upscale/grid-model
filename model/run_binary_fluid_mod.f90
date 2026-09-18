module run_binary_fluid_mod

    use table2d_data
    use table2d_io
    use table2d_eval
    use init_profiles
    use timestep_control
    use convergence_control

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





        ! ! ---- local state ----
        ! type(table2d_t) :: tbl

        ! integer(kind=8) :: iter
        ! integer :: number_block, number_edge, i, block, edge, edge1, edge2
        ! integer :: log_unit, conv_unit, steady_count

        ! real(8) :: block_area_yz, block_volume_xyz
        ! real(8) :: time_step, time_step_A, time_step_B
        ! real(8) :: left_muA, right_muA, inside_muA
        ! real(8) :: left_muB, right_muB, inside_muB
        ! real(8) :: net_flux_A, net_flux_B
        ! real(8) :: force_edge_A, force_edge_B, flux_edge_A, flux_edge_B
        ! real(8) :: mAA_edge, mAB_edge, mBB_edge
        ! real(8) :: flux_mean_A, flux_std_A, flux_conservation_A
        ! real(8) :: flux_mean_B, flux_std_B, flux_conservation_B
        ! real(8) :: time, prev_mean_density_A, prev_mean_density_B
        ! logical :: converged, converged_A, converged_B, use_steady_state
        ! logical :: inversion_ok

        ! real(8), allocatable :: block_centers(:)
        ! real(8), allocatable :: block_edges(:)
        ! real(8), allocatable :: muA_field(:), muB_field(:)
        ! real(8), allocatable :: rhoA_field(:), rhoB_field(:)
        ! real(8), allocatable :: delta_rhoA(:), delta_rhoB(:)
        ! real(8), allocatable :: mAA_block(:), mAB_block(:), mBB_block(:)
        ! real(8), allocatable :: flux_edges_A(:), flux_edges_B(:)
        ! real(8), allocatable :: grad_muA(:), grad_muB(:)
        ! real(8), allocatable :: p_dummy(:)   ! forward-eval also returns p; not used in dynamics

        ! real(8), parameter :: Na = 6.022e23             ! mol-1
        ! real(8), parameter :: kcal_to_j = 4184.0d0       ! J/kcal
        ! real(8), parameter :: mu_equal_tol = 1d-6        ! J/mol

        ! ! ---- unit conversion ----
        ! left_muA   = left_muA_in   * kcal_to_j
        ! right_muA  = right_muA_in  * kcal_to_j
        ! inside_muA = inside_muA_in * kcal_to_j

        ! left_muB   = left_muB_in   * kcal_to_j
        ! right_muB  = right_muB_in  * kcal_to_j
        ! inside_muB = inside_muB_in * kcal_to_j

        ! time_step = time_step_in

        ! ! Same rationale as single-fluid: if BOTH species have symmetric
        ! ! reservoirs, flux_mean -> 0 for both and flux-conservation can
        ! ! never converge, so density steady-state is used instead.
        ! use_steady_state = (abs(right_muA - left_muA) < mu_equal_tol) .and. &
        !                     (abs(right_muB - left_muB) < mu_equal_tol)

        ! ! Load 2D equation-of-state and Onsager mobility-matrix table
        ! call load_table2d(trim(table2d_dir), tbl)

        ! ! System definition
        ! number_block = nint(system_size_x / block_size_x)
        ! number_edge = number_block - 1
        ! block_area_yz = block_size_y*block_size_z
        ! block_volume_xyz = block_size_x*block_area_yz

        ! allocate(block_centers(number_block))
        ! allocate(muA_field(number_block), muB_field(number_block))
        ! allocate(rhoA_field(number_block), rhoB_field(number_block))
        ! allocate(delta_rhoA(number_block), delta_rhoB(number_block))
        ! allocate(mAA_block(number_block), mAB_block(number_block), mBB_block(number_block))
        ! allocate(p_dummy(number_block))

        ! allocate(block_edges(number_edge))
        ! allocate(grad_muA(number_edge), grad_muB(number_edge))
        ! allocate(flux_edges_A(number_edge), flux_edges_B(number_edge))

        ! do block = 1, number_block
        !     block_centers(block) = (block-1) * block_size_x + block_size_x/2
        ! end do
        ! do edge = 1, number_edge
        !     block_edges(edge) = (edge-1) * block_size_x
        ! end do

        ! ! Initial chemical-potential profiles, one per species
        ! call init_mu_profile(muA_field, number_block, left_muA, right_muA, inside_muA, mu_mode)
        ! call init_mu_profile(muB_field, number_block, left_muB, right_muB, inside_muB, mu_mode)

        ! ! Initial densities from the 2D table (forward map)
        ! do block = 1, number_block
        !     call eval_forward(muA_field(block), muB_field(block), tbl, &
        !                        rhoA_field(block), rhoB_field(block), p_dummy(block))
        ! end do

        ! flux_edges_A = 0.0d0; flux_edges_B = 0.0d0
        ! mAA_block = 0.0d0; mAB_block = 0.0d0; mBB_block = 0.0d0
        ! grad_muA = 0.0d0; grad_muB = 0.0d0
        ! delta_rhoA = 0.0d0; delta_rhoB = 0.0d0

        ! log_unit = 99
        ! conv_unit = 97

        ! open(newunit=log_unit, file=trim(output_dir)//"/grid.log", status="replace", action="write")
        ! write(log_unit,*) "Binary-fluid simulation started"
        ! write(log_unit,*) "Number of iterations =", n_iter

        ! open(newunit=conv_unit, file=trim(output_dir)//"/conservation.dat", status="replace", action="write")
        ! if (use_steady_state) then
        !     write(conv_unit,*) "# iter    density_drift_A    density_drift_B    time_step[s]"
        ! else
        !     write(conv_unit,*) "# iter    flux_conservation_A    flux_conservation_B    &
        !                         &flux_mean_A[1/s]    flux_mean_B[1/s]    time_step[s]"
        ! end if

        ! iter = 0
        ! time = 0d0
        ! steady_count = 0
        ! prev_mean_density_A = -1d0
        ! prev_mean_density_B = -1d0
        ! flux_mean_A = 0d0
        ! flux_mean_B = 0d0

        ! do while (iter < n_iter)

        !     iter = iter + 1
        !     time = time + time_step

        !     ! ---- outputs ----
        !     if (mod(iter, n_jump) == 0) then
        !         call write_profiles_binary(int(iter/n_jump), time, &
        !                                     block_centers, block_edges, &
        !                                     muA_field, muB_field, &
        !                                     rhoA_field, rhoB_field, &
        !                                     mAA_block, mAB_block, mBB_block, &
        !                                     flux_edges_A, flux_edges_B, &
        !                                     grad_muA, grad_muB, &
        !                                     number_block, number_edge, output_dir)

        !         write(log_unit,*) "iter =", iter
        !         write(*,'(A,I10,A,ES10.3,A,ES10.3,A,ES10.3)') &
        !             "  iter=", iter, &
        !             "  dt=", time_step, &
        !             "  fluxA_mean=", flux_mean_A, &
        !             "  fluxB_mean=", flux_mean_B
        !     end if

        !     ! ---- gradients, evaluated before any update ----
        !     do block = 2, number_block
        !         edge = block-1
        !         grad_muA(edge) = (muA_field(block) - muA_field(block-1)) / block_size_x
        !         grad_muB(edge) = (muB_field(block) - muB_field(block-1)) / block_size_x
        !     end do

        !     ! ---- Onsager mobility matrix, evaluated before any update ----
        !     ! NOTE: as in the single-fluid model, block=1 and block=number_block
        !     ! are wrong (no entrance-effect theory applied).
        !     do block = 1, number_block
        !         call eval_mobility(muA_field(block), muB_field(block), tbl, &
        !                             mAA_block(block), mAB_block(block), mBB_block(block))
        !     end do

        !     do edge = 1, number_edge
        !         block = edge + 1

        !         mAA_edge = (mAA_block(block) + mAA_block(block-1))/2
        !         mAB_edge = (mAB_block(block) + mAB_block(block-1))/2
        !         mBB_edge = (mBB_block(block) + mBB_block(block-1))/2

        !         force_edge_A = -grad_muA(edge) / Na
        !         force_edge_B = -grad_muB(edge) / Na

        !         flux_edge_A = (mAA_edge * force_edge_A + mAB_edge * force_edge_B) * block_area_yz
        !         flux_edge_B = (mAB_edge * force_edge_A + mBB_edge * force_edge_B) * block_area_yz

        !         flux_edges_A(edge) = flux_edge_A
        !         flux_edges_B(edge) = flux_edge_B
        !     end do

        !     ! ---- update densities, then invert back to mu (2D Newton) ----
        !     do block = 2, number_block-1

        !         edge1 = block-1
        !         edge2 = block

        !         net_flux_A = flux_edges_A(edge1) - flux_edges_A(edge2)
        !         net_flux_B = flux_edges_B(edge1) - flux_edges_B(edge2)

        !         delta_rhoA(block) = net_flux_A * time_step / block_volume_xyz
        !         delta_rhoB(block) = net_flux_B * time_step / block_volume_xyz

        !         rhoA_field(block) = rhoA_field(block) + delta_rhoA(block)
        !         rhoB_field(block) = rhoB_field(block) + delta_rhoB(block)

        !         ! warm-start the inversion from this block's previous mu
        !         call mu_from_rho2d(rhoA_field(block), rhoB_field(block), &
        !                             muA_field(block), muB_field(block), tbl, &
        !                             muA_field(block), muB_field(block), inversion_ok)

        !         if (.not. inversion_ok) then
        !             write(*,*) "mu_from_rho2d failed to converge at iter =", iter, &
        !                         " block =", block
        !             write(*,*) "  rhoA =", rhoA_field(block), "  rhoB =", rhoB_field(block)
        !             write(log_unit,*) "Inversion failure at iter =", iter, " block =", block, &
        !                                " — stopping"
        !             stop "Binary-fluid inversion failed — simulation aborted"
        !         end if

        !     end do

        !     ! ---- NaN / negative-density checks + timestep adaptation ----
        !     if (mod(iter, check_interval) == 0) then

        !         do i = 2, number_block-1
        !             if (isnan(rhoA_field(i)) .or. isnan(rhoB_field(i)) .or. &
        !                 isnan(muA_field(i))  .or. isnan(muB_field(i))) then
        !                 write(*,*) "NaN detected at iter =", iter, " block =", i
        !                 write(log_unit,*) "NaN detected at iter =", iter, " block =", i, " — stopping"
        !                 stop "NaN detected — simulation aborted"
        !             end if
        !         end do

        !         do i = 2, number_block-1
        !             if (rhoA_field(i) <= 0d0 .or. rhoB_field(i) <= 0d0) then
        !                 write(*,*) "Negative density at iter =", iter, " block =", i
        !                 write(log_unit,*) "Negative density at iter =", iter, " block =", i, " — stopping"
        !                 stop "Negative density — simulation aborted"
        !             end if
        !         end do

        !         ! Adapt timestep independently for each species, then take
        !         ! the more restrictive (smaller) of the two proposed steps.
        !         time_step_A = time_step
        !         time_step_B = time_step
        !         call adapt_timestep(time_step_A, delta_rhoA, rhoA_field, number_block, &
        !                             tol_min, tol_max, growth_factor, max_time_step)
        !         call adapt_timestep(time_step_B, delta_rhoB, rhoB_field, number_block, &
        !                             tol_min, tol_max, growth_factor, max_time_step)
        !         time_step = min(time_step_A, time_step_B)

        !     end if

        !     ! ---- convergence check ----
        !     if (mod(iter, check_interval) == 0) then

        !         if (use_steady_state) then
        !             call check_steady_state(iter, rhoA_field, &
        !                                      number_block, steady_tol, n_steady, &
        !                                      time_step, conv_unit, &
        !                                      steady_count, prev_mean_density_A, converged_A)
        !             call check_steady_state(iter, rhoB_field, &
        !                                      number_block, steady_tol, n_steady, &
        !                                      time_step, conv_unit, &
        !                                      steady_count, prev_mean_density_B, converged_B)
        !         else
        !             call check_flux_conservation(iter, check_interval, time_step, &
        !                                           flux_edges_A, number_block, conv_tol, &
        !                                           conv_unit, flux_mean_A, flux_std_A, &
        !                                           flux_conservation_A, converged_A)
        !             call check_flux_conservation(iter, check_interval, time_step, &
        !                                           flux_edges_B, number_block, conv_tol, &
        !                                           conv_unit, flux_mean_B, flux_std_B, &
        !                                           flux_conservation_B, converged_B)
        !         end if

        !         ! Require BOTH species to satisfy the convergence criterion.
        !         ! NOTE: check_steady_state/check_flux_conservation each write
        !         ! their own line to conv_unit, so calling them twice per
        !         ! check doubles the row count in conservation.dat (once for
        !         ! A, once for B) — fine for post-processing but worth knowing.
        !         converged = converged_A .and. converged_B

        !         if (converged) then
        !             write(*,*) "Converged at iter =", iter
        !             write(log_unit,*) "Converged at iter =", iter

        !             call write_profiles_binary(int(iter/n_jump), time, &
        !                                         block_centers, block_edges, &
        !                                         muA_field, muB_field, &
        !                                         rhoA_field, rhoB_field, &
        !                                         mAA_block, mAB_block, mBB_block, &
        !                                         flux_edges_A, flux_edges_B, &
        !                                         grad_muA, grad_muB, &
        !                                         number_block, number_edge, output_dir, &
        !                                         label="final")

        !             close(log_unit)
        !             close(conv_unit)
        !             return
        !         end if

        !     end if

        ! end do

        ! close(log_unit)
        ! close(conv_unit)

    end subroutine run_binary_fluid

    ! ------------------------------------------------------------
    ! Minimal profile writer for the two-species case. Mirrors the
    ! layout of io_profiles::write_profiles but with muA/muB and
    ! rhoA/rhoB side by side instead of a single mu/rho column.
    ! Replace with a proper io_profiles_binary module if this needs
    ! to match existing post-processing scripts.
    ! ------------------------------------------------------------
    subroutine write_profiles_binary(snap_index, time, block_centers, block_edges, &
                                      muA_field, muB_field, rhoA_field, rhoB_field, &
                                      mAA_block, mAB_block, mBB_block, &
                                      flux_edges_A, flux_edges_B, &
                                      grad_muA, grad_muB, &
                                      number_block, number_edge, output_dir, label)
        implicit none
        integer, intent(in) :: snap_index, number_block, number_edge
        real(8), intent(in) :: time
        real(8), intent(in) :: block_centers(number_block), block_edges(number_edge)
        real(8), intent(in) :: muA_field(number_block), muB_field(number_block)
        real(8), intent(in) :: rhoA_field(number_block), rhoB_field(number_block)
        real(8), intent(in) :: mAA_block(number_block), mAB_block(number_block), mBB_block(number_block)
        real(8), intent(in) :: flux_edges_A(number_edge), flux_edges_B(number_edge)
        real(8), intent(in) :: grad_muA(number_edge), grad_muB(number_edge)
        character(len=*), intent(in) :: output_dir
        character(len=*), intent(in), optional :: label

        integer :: unit, i
        character(len=256) :: filename
        character(len=32)  :: tag

        if (present(label)) then
            tag = trim(label)
        else
            write(tag,'(I0)') snap_index
        end if

        filename = trim(output_dir)//"/profile_block_"//trim(tag)//".dat"
        open(newunit=unit, file=trim(filename), status="replace", action="write")
        write(unit,*) "# time =", time, " s"
        write(unit,*) "# x[m]  muA[J/mol]  muB[J/mol]  rhoA[m^-3]  rhoB[m^-3]  M_AA  M_AB  M_BB"
        do i = 1, number_block
            write(unit,'(8ES16.8)') block_centers(i), muA_field(i), muB_field(i), &
                                     rhoA_field(i), rhoB_field(i), &
                                     mAA_block(i), mAB_block(i), mBB_block(i)
        end do
        close(unit)

        filename = trim(output_dir)//"/profile_edge_"//trim(tag)//".dat"
        open(newunit=unit, file=trim(filename), status="replace", action="write")
        write(unit,*) "# time =", time, " s"
        write(unit,*) "# x[m]  grad_muA  grad_muB  flux_A[1/s]  flux_B[1/s]"
        do i = 1, number_edge
            write(unit,'(5ES16.8)') block_edges(i), grad_muA(i), grad_muB(i), &
                                     flux_edges_A(i), flux_edges_B(i)
        end do
        close(unit)

    end subroutine write_profiles_binary

end module run_binary_fluid_mod