module timestep_control
    implicit none
contains
    subroutine adapt_timestep(time_step, delta_density, fluid_density, &
                             number_block, tol_min, tol_max, growth_factor, max_time_step)
        implicit none
        real(8), intent(inout) :: time_step
        real(8), intent(in) :: delta_density(:), fluid_density(:)
        integer, intent(in) :: number_block
        real(8), intent(in) :: tol_min, tol_max, growth_factor, max_time_step
        real(8) :: max_rel_change
        integer :: i

        max_rel_change = 0d0
        do i = 2, number_block-1
            max_rel_change = max(max_rel_change, abs(delta_density(i)) / max(fluid_density(i), 1d-30))
        end do

        if (max_rel_change < tol_min) then
            time_step = time_step * growth_factor
            if (time_step > max_time_step) time_step = max_time_step
        else if (max_rel_change > tol_max) then
            time_step = time_step / growth_factor
        end if
    end subroutine adapt_timestep
end module timestep_control