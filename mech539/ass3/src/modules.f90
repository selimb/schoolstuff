module types
    implicit none
    public
    integer, parameter :: dp=kind(1.d0)
end module

module grid
    use types, only: dp
    implicit none
    public
    real(dp), parameter :: tc = 0.08_dp
    integer :: nx, ny
    real(dp), dimension(:), allocatable :: x, y
end module

module params
   use types, only: dp
   implicit none
   public
   real(dp), parameter :: gamma = 1.4_dp
   real(dp) :: Minf
   integer :: solverID
end module
