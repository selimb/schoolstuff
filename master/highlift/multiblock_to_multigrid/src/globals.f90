module globals
    public
    ! Dimensions
    integer :: nblocks
    integer :: imax, jmax, lmax
    integer, parameter :: kmax = 9
!   integer, parameter :: lmax = max(imax, jmax)
    integer, dimension(:), allocatable :: ni, nj, nk
    ! Grid points
    real*8, dimension(:,:,:,:), allocatable :: X2D, X3D
!   real*8, dimension(2,imax,jmax,nbmax) :: X2D
!   real*8, dimension(3,imax,jmax,kmax) :: X3D
    ! Connectivity
    integer, parameter :: num_axes = 2
    integer, parameter :: num_dirs = 2
    integer, dimension(:, :, :), allocatable :: ITYPE, BTYPE
    integer, dimension(:, :, :, :), allocatable :: IA
!   integer, dimension(num_dirs, num_axes, nbmax) :: ITYPE
!   integer, dimension(3, num_dirs, num_axes, nbmax) :: IA
!   integer, dimension(num_dirs, num_axes, nbmax) :: BTYPE
    integer, parameter :: KTYPE = 0
    integer, parameter :: IJK = 1
    integer, parameter :: ROT = 0
    ! Bookkeeping
!   integer :: num_undone
!   integer, dimension(nbmax) :: UNDONE
    contains
subroutine set_sizes
    imax = maxval(ni(1:nblocks)) + 6
    jmax = maxval(nj(1:nblocks)) + 6
    write(*,*) imax, jmax
    lmax = max(imax, jmax)
    allocate(X2D(2,imax,jmax,nblocks))
    allocate(X3D(3,imax,jmax,kmax))
    allocate(ITYPE(num_dirs, num_axes, nblocks))
    allocate(IA(3, num_dirs, num_axes, nblocks))
    allocate(BTYPE(num_dirs, num_axes, nblocks))
    ITYPE = 0
    IA = 0
    BTYPE = -1
    X2D = 0
    X3D = 0
end subroutine

pure function NBL(ax, m)
    integer, intent(in) :: ax, m
    integer :: NBL
    if (ax .eq. 1) NBL = ni(m)
    if (ax .eq. 2) NBL = nj(m)
end function

end module
