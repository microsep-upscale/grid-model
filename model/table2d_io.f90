module table2d_io
  use table2d_data
  implicit none
  public :: load_table2d

contains

  subroutine load_table2d(dir, tbl)
    implicit none
    character(len=*), intent(in)  :: dir
    type(table2d_t),  intent(out) :: tbl
    integer :: unit, i, j
    integer, allocatable :: safe_int(:,:)

    ! ---- axes ----
    tbl%nA = count_lines(trim(dir)//"/muA_45x45.dat")
    tbl%nB = count_lines(trim(dir)//"/muB_45x45.dat")

    allocate(tbl%muA(tbl%nA), tbl%muB(tbl%nB))

    open(newunit=unit, file=trim(dir)//"/muA_45x45.dat", status="old", action="read")
    read(unit,*) (tbl%muA(i), i=1,tbl%nA)
    close(unit)

    open(newunit=unit, file=trim(dir)//"/muB_45x45.dat", status="old", action="read")
    read(unit,*) (tbl%muB(j), j=1,tbl%nB)
    close(unit)

    tbl%dmuA = tbl%muA(2) - tbl%muA(1)
    tbl%dmuB = tbl%muB(2) - tbl%muB(1)

    ! ---- fields ----
    allocate(tbl%rhoA(tbl%nA,tbl%nB), tbl%rhoB(tbl%nA,tbl%nB), tbl%p(tbl%nA,tbl%nB))
    allocate(tbl%dA_dA(tbl%nA,tbl%nB), tbl%dA_dB(tbl%nA,tbl%nB))
    allocate(tbl%dB_dA(tbl%nA,tbl%nB), tbl%dB_dB(tbl%nA,tbl%nB))
    allocate(safe_int(tbl%nA,tbl%nB), tbl%safe(tbl%nA,tbl%nB))

    call read_matrix(trim(dir)//"/rhoA_reconstructed_LS_45x45.dat", tbl%rhoA, tbl%nA, tbl%nB)
    call read_matrix(trim(dir)//"/rhoB_reconstructed_LS_45x45.dat", tbl%rhoB, tbl%nA, tbl%nB)
    call read_matrix(trim(dir)//"/p_reconstructed_LS_45x45.dat",    tbl%p,    tbl%nA, tbl%nB)

    allocate(tbl%M_AA(tbl%nA,tbl%nB), tbl%M_AB(tbl%nA,tbl%nB), tbl%M_BB(tbl%nA,tbl%nB))

    call read_matrix(trim(dir)//"/M_AA_45x45.dat", tbl%M_AA, tbl%nA, tbl%nB)
    call read_matrix(trim(dir)//"/M_AB_45x45.dat", tbl%M_AB, tbl%nA, tbl%nB)
    call read_matrix(trim(dir)//"/M_BB_45x45.dat", tbl%M_BB, tbl%nA, tbl%nB)

    open(newunit=unit, file=trim(dir)//"/safe_mask_45x45.dat", status="old", action="read")
    do i = 1, tbl%nA
       read(unit,*) (safe_int(i,j), j=1,tbl%nB)
    end do
    close(unit)
    tbl%safe = (safe_int == 1)

    ! ---- local Jacobian by central differences (interior),
    !      one-sided at the boundary ----
    call central_diff_axis1(tbl%rhoA, tbl%muA, tbl%nA, tbl%nB, tbl%dA_dA)
    call central_diff_axis2(tbl%rhoA, tbl%muB, tbl%nA, tbl%nB, tbl%dA_dB)
    call central_diff_axis1(tbl%rhoB, tbl%muA, tbl%nA, tbl%nB, tbl%dB_dA)
    call central_diff_axis2(tbl%rhoB, tbl%muB, tbl%nA, tbl%nB, tbl%dB_dB)

  end subroutine load_table2d

  ! ---------------------------------------------------------------

  integer function count_lines(filename) result(n)
    implicit none
    character(len=*), intent(in) :: filename
    integer :: unit, ios
    real(8) :: dummy
    n = 0
    open(newunit=unit, file=filename, status="old", action="read")
    do
       read(unit,*,iostat=ios) dummy
       if (ios /= 0) exit
       n = n + 1
    end do
    close(unit)
  end function count_lines

  subroutine read_matrix(filename, mat, nA, nB)
    implicit none
    character(len=*), intent(in)  :: filename
    integer,          intent(in)  :: nA, nB
    real(8),          intent(out) :: mat(nA,nB)
    integer :: unit, i, j
    open(newunit=unit, file=filename, status="old", action="read")
    do i = 1, nA
       read(unit,*) (mat(i,j), j=1,nB)
    end do
    close(unit)
  end subroutine read_matrix

  subroutine central_diff_axis1(f, x, nA, nB, df)
    ! d f / d x, x varies along first (row) index
    implicit none
    integer, intent(in)  :: nA, nB
    real(8), intent(in)  :: f(nA,nB), x(nA)
    real(8), intent(out) :: df(nA,nB)
    integer :: i, j
    do j = 1, nB
       df(1,j)  = (f(2,j)  - f(1,j))    / (x(2)  - x(1))
       df(nA,j) = (f(nA,j) - f(nA-1,j)) / (x(nA) - x(nA-1))
       do i = 2, nA-1
          df(i,j) = (f(i+1,j) - f(i-1,j)) / (x(i+1) - x(i-1))
       end do
    end do
  end subroutine central_diff_axis1

  subroutine central_diff_axis2(f, y, nA, nB, df)
    ! d f / d y, y varies along second (column) index
    implicit none
    integer, intent(in)  :: nA, nB
    real(8), intent(in)  :: f(nA,nB), y(nB)
    real(8), intent(out) :: df(nA,nB)
    integer :: i, j
    do i = 1, nA
       df(i,1)  = (f(i,2)  - f(i,1))    / (y(2)  - y(1))
       df(i,nB) = (f(i,nB) - f(i,nB-1)) / (y(nB) - y(nB-1))
       do j = 2, nB-1
          df(i,j) = (f(i,j+1) - f(i,j-1)) / (y(j+1) - y(j-1))
       end do
    end do
  end subroutine central_diff_axis2

end module table2d_io