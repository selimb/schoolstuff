module globals
    public
    ! Dimensions
    integer, parameter :: imax = 1684
    integer, parameter :: jmax = 302
    integer, parameter :: kmax = 5
    integer, parameter :: lmax = max(imax, jmax)
    integer, parameter :: nbmax = 200
    integer, dimension(nbmax) :: ni, nj, nk
    ! Grid points
    real*8, dimension(2,imax,jmax,nbmax) :: X2D
    real*8, dimension(3,imax,jmax,kmax) :: X3D
    ! Connectivity
    integer, parameter :: num_axes = 2
    integer, parameter :: num_dirs = 2
    integer, dimension(num_dirs, num_axes, nbmax) :: ITYPE
    integer, dimension(3, num_dirs, num_axes, nbmax) :: IA
    integer, dimension(num_dirs, num_axes, nbmax) :: BTYPE
    integer, parameter :: KTYPE = 0
    integer, parameter :: ROT = 0
    ! Blocks
    integer :: nblocks
    integer, dimension(-10:nbmax) :: IDX
end module
