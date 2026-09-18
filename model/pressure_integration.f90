module pressure_integration

    use io_profiles

    implicit none
    private

    public :: compute_pressure_fields

contains

    subroutine compute_pressure_fields(muA_fine, muB_fine, rhoA_fine, rhoB_fine, &
                                        p_fine, p_alt, p_avg, output_dir)
        implicit none
        real(8), intent(in)  :: muA_fine(:), muB_fine(:)
        real(8), intent(in)  :: rhoA_fine(:,:), rhoB_fine(:,:)
        real(8), allocatable, intent(out) :: p_fine(:,:), p_alt(:,:), p_avg(:,:)
        character(len=*), intent(in) :: output_dir

        integer :: i, j, grid_size, ip, jp
        real(8) :: dmuA, dmuB
        real(8), allocatable :: p_path_error(:,:)
        real(8) :: max_err, rms_err
        character(len=64) :: size_str

        grid_size = size(muA_fine)

        allocate(p_fine(grid_size, grid_size))
        allocate(p_alt(grid_size, grid_size))
        allocate(p_avg(grid_size, grid_size))
        allocate(p_path_error(grid_size, grid_size))

        p_fine = 0d0
        p_alt  = 0d0

        ! ---- p_fine: integrate along muA first (at j=1), then along muB ----
        do i = 2, grid_size
            dmuA = muA_fine(i) - muA_fine(i-1)
            p_fine(i,1) = p_fine(i-1,1) + 0.5d0*(rhoA_fine(i-1,1) + rhoA_fine(i,1))*dmuA
        end do

        do i = 1, grid_size
            do j = 2, grid_size
                dmuB = muB_fine(j) - muB_fine(j-1)
                p_fine(i,j) = p_fine(i,j-1) + 0.5d0*(rhoB_fine(i,j-1) + rhoB_fine(i,j))*dmuB
            end do
        end do

        ! ---- p_alt: integrate along muB first (at i=1), then along muA ----
        do j = 2, grid_size
            dmuB = muB_fine(j) - muB_fine(j-1)
            p_alt(1,j) = p_alt(1,j-1) + 0.5d0*(rhoB_fine(1,j-1) + rhoB_fine(1,j))*dmuB
        end do

        do j = 1, grid_size
            do i = 2, grid_size
                dmuA = muA_fine(i) - muA_fine(i-1)
                p_alt(i,j) = p_alt(i-1,j) + 0.5d0*(rhoA_fine(i-1,j) + rhoA_fine(i,j))*dmuA
            end do
        end do

        ! ---- Path-dependence diagnostics ----
        p_path_error = p_fine - p_alt
        max_err = maxval(abs(p_path_error))
        rms_err = sqrt(sum(p_path_error**2) / real(grid_size*grid_size, 8))

        ip = 1
        jp = 1
        do i = 1, grid_size
            do j = 1, grid_size
                if (abs(p_path_error(i,j)) > abs(p_path_error(ip,jp))) then
                    ip = i
                    jp = j
                end if
            end do
        end do

        ! write(*,*) "max path difference = ", max_err
        ! write(*,*) "RMS path difference = ", rms_err
        ! write(*,*) "worst path difference at: muA=", muA_fine(ip), " muB=", muB_fine(jp), &
        !             " = ", p_path_error(ip,jp)

        ! ---- Symmetrized pressure field ----
        p_avg = 0.5d0 * (p_fine + p_alt)

        ! ---- Write outputs ----
        write(size_str,'(I0,A,I0)') grid_size, "x", grid_size

        call write_grid2d(trim(output_dir)//"/p_fine_"//trim(size_str)//".dat", &
                           p_fine, "# pressure field, muA-first path")
        call write_grid2d(trim(output_dir)//"/p_alt_"//trim(size_str)//".dat", &
                           p_alt, "# pressure field, muB-first path")
        call write_grid2d(trim(output_dir)//"/p_avg_"//trim(size_str)//".dat", &
                           p_avg, "# pressure field, path-averaged")
        call write_grid2d(trim(output_dir)//"/p_error_"//trim(size_str)//".dat", &
                           p_path_error, "# path integration error (p_fine - p_alt)")

        deallocate(p_path_error)

    end subroutine compute_pressure_fields

end module pressure_integration