module table2d_eval
  use table2d_data
  implicit none
  public :: eval_forward, eval_jacobian, mu_from_rho2d

contains

  ! Locate the lower-left grid index of the cell containing (muA_q,
  ! muB_q), and the fractional position (tA,tB) in [0,1] within it.
  ! Uniform-grid assumption gives O(1) lookup. Queries outside the
  ! table are clamped to the boundary cell rather than extrapolated.
  subroutine locate_cell(muA_q, muB_q, tbl, iA, iB, tA, tB)
    implicit none
    real(8),         intent(in)  :: muA_q, muB_q
    type(table2d_t), intent(in)  :: tbl
    integer,         intent(out) :: iA, iB
    real(8),         intent(out) :: tA, tB

    iA = floor((muA_q - tbl%muA(1)) / tbl%dmuA) + 1
    iB = floor((muB_q - tbl%muB(1)) / tbl%dmuB) + 1

    iA = max(1, min(tbl%nA - 1, iA))
    iB = max(1, min(tbl%nB - 1, iB))

    tA = (muA_q - tbl%muA(iA)) / tbl%dmuA
    tB = (muB_q - tbl%muB(iB)) / tbl%dmuB

    tA = max(0d0, min(1d0, tA))
    tB = max(0d0, min(1d0, tB))

  end subroutine locate_cell

  real(8) function bilinear(field, tbl, iA, iB, tA, tB) result(val)
    implicit none
    real(8),         intent(in) :: field(:,:)
    type(table2d_t), intent(in) :: tbl
    integer,         intent(in) :: iA, iB
    real(8),         intent(in) :: tA, tB

    val = (1d0-tA)*(1d0-tB) * field(iA,   iB)   &
        +      tA *(1d0-tB) * field(iA+1, iB)   &
        + (1d0-tA)*     tB  * field(iA,   iB+1) &
        +      tA *     tB  * field(iA+1, iB+1)

  end function bilinear

  ! Forward map: (muA,muB) -> (rhoA,rhoB,p)
  subroutine eval_forward(muA_q, muB_q, tbl, rhoA_q, rhoB_q, p_q)
    implicit none
    real(8),         intent(in)  :: muA_q, muB_q
    type(table2d_t), intent(in)  :: tbl
    real(8),         intent(out) :: rhoA_q, rhoB_q, p_q
    integer :: iA, iB
    real(8) :: tA, tB

    call locate_cell(muA_q, muB_q, tbl, iA, iB, tA, tB)

    rhoA_q = bilinear(tbl%rhoA, tbl, iA, iB, tA, tB)
    rhoB_q = bilinear(tbl%rhoB, tbl, iA, iB, tA, tB)
    p_q    = bilinear(tbl%p,    tbl, iA, iB, tA, tB)

  end subroutine eval_forward

  ! Local susceptibility matrix, interpolated rather than
  ! re-differentiated — an approximate Jacobian is enough for Newton.
  subroutine eval_jacobian(muA_q, muB_q, tbl, a, b, c, d)
    implicit none
    real(8),         intent(in)  :: muA_q, muB_q
    type(table2d_t), intent(in)  :: tbl
    real(8),         intent(out) :: a, b, c, d
    integer :: iA, iB
    real(8) :: tA, tB

    call locate_cell(muA_q, muB_q, tbl, iA, iB, tA, tB)

    a = bilinear(tbl%dA_dA, tbl, iA, iB, tA, tB)
    b = bilinear(tbl%dA_dB, tbl, iA, iB, tA, tB)
    c = bilinear(tbl%dB_dA, tbl, iA, iB, tA, tB)
    d = bilinear(tbl%dB_dB, tbl, iA, iB, tA, tB)

  end subroutine eval_jacobian

  ! bilinear lookup of the mobility matrix at a query point
  subroutine eval_mobility(muA_q, muB_q, tbl, mAA, mAB, mBB)
    implicit none
    real(8),         intent(in)  :: muA_q, muB_q
    type(table2d_t), intent(in)  :: tbl
    real(8),         intent(out) :: mAA, mAB, mBB
    integer :: iA, iB
    real(8) :: tA, tB

    call locate_cell(muA_q, muB_q, tbl, iA, iB, tA, tB)

    mAA = bilinear(tbl%M_AA, tbl, iA, iB, tA, tB)
    mAB = bilinear(tbl%M_AB, tbl, iA, iB, tA, tB)
    mBB = bilinear(tbl%M_BB, tbl, iA, iB, tA, tB)

  end subroutine eval_mobility

  ! Inverse map: (rhoA,rhoB) -> (muA,muB) via damped 2D Newton,
  ! warm-started from (muA_guess,muB_guess). Fortran analogue of
  ! mu_from_rho() in the Python reference implementation.
  subroutine mu_from_rho2d(rhoA_t, rhoB_t, muA_guess, muB_guess, tbl, &
                            muA_sol, muB_sol, converged)
    implicit none
    real(8),         intent(in)  :: rhoA_t, rhoB_t
    real(8),         intent(in)  :: muA_guess, muB_guess
    type(table2d_t), intent(in)  :: tbl
    real(8),         intent(out) :: muA_sol, muB_sol
    logical,         intent(out) :: converged

    real(8), parameter :: tol      = 1d-10
    integer, parameter :: max_iter = 50
    real(8), parameter :: max_step = 0.5d0   ! per-iteration damping cap

    real(8) :: muA, muB, rhoA_c, rhoB_c, p_c
    real(8) :: resA, resB, a, b, c, d, det
    real(8) :: dmuA_step, dmuB_step, step_norm
    integer :: iter, iA, iB
    real(8) :: tA, tB

    muA = muA_guess
    muB = muB_guess
    converged = .false.

    do iter = 1, max_iter

       call eval_forward(muA, muB, tbl, rhoA_c, rhoB_c, p_c)

       resA = rhoA_c - rhoA_t
       resB = rhoB_c - rhoB_t

       if (sqrt(resA**2 + resB**2) < tol) then
          converged = .true.
          exit
       end if

       call eval_jacobian(muA, muB, tbl, a, b, c, d)

       det = a*d - b*c
       if (abs(det) < 1d-14) exit   ! near-singular susceptibility, bail out

       dmuA_step = -( d*resA - b*resB) / det
       dmuB_step = -(-c*resA + a*resB) / det

       step_norm = sqrt(dmuA_step**2 + dmuB_step**2)
       if (step_norm > max_step) then
          dmuA_step = dmuA_step * max_step / step_norm
          dmuB_step = dmuB_step * max_step / step_norm
       end if

       muA = muA + dmuA_step
       muB = muB + dmuB_step

       muA = max(tbl%muA(1), min(tbl%muA(tbl%nA), muA))
       muB = max(tbl%muB(1), min(tbl%muB(tbl%nB), muB))

    end do

    muA_sol = muA
    muB_sol = muB
    if (.not. converged) return

    ! reject solutions landing in an unsafe (unstable/spinodal) cell
    call locate_cell(muA_sol, muB_sol, tbl, iA, iB, tA, tB)
    if (.not. (tbl%safe(iA,iB)   .and. tbl%safe(iA+1,iB) .and. &
               tbl%safe(iA,iB+1) .and. tbl%safe(iA+1,iB+1))) then
       converged = .false.
    end if

  end subroutine mu_from_rho2d

end module table2d_eval