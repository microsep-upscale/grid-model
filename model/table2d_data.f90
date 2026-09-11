module table2d_data
  implicit none

  type :: table2d_t
     integer :: nA, nB
     real(8), allocatable :: muA(:)          ! (nA)
     real(8), allocatable :: muB(:)          ! (nB)
     real(8) :: dmuA, dmuB                   ! uniform grid spacing

     real(8), allocatable :: rhoA(:,:)       ! (nA,nB)
     real(8), allocatable :: rhoB(:,:)
     real(8), allocatable :: p(:,:)

     ! local susceptibility (Jacobian of rho wrt mu), central differences
     real(8), allocatable :: dA_dA(:,:)      ! d(rhoA)/d(muA)
     real(8), allocatable :: dA_dB(:,:)      ! d(rhoA)/d(muB)
     real(8), allocatable :: dB_dA(:,:)      ! d(rhoB)/d(muA)
     real(8), allocatable :: dB_dB(:,:)      ! d(rhoB)/d(muB)

     ! Onsager mobility matrix fields, symmetric: MAB == MBA
     real(8), allocatable :: M_AA(:,:)
     real(8), allocatable :: M_AB(:,:)   ! = M_BA(:,:)
     real(8), allocatable :: M_BB(:,:)

     logical, allocatable :: safe(:,:)
  end type table2d_t

end module table2d_data