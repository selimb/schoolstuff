module solver
use types, only: dp
use grid, only: nx, ny, x, y
use params, only: Minf, solverID, gamma, Ttot_in, Rgas
implicit none
real(dp), public, parameter :: tol = 1e-15
integer, public, parameter  :: itermax = 1e7
real(dp), public :: Uinf, ainf
private
real(dp) :: dy0, mMinf2, gp1m2u
! xlim1 and xlim2 correspond to the indices enclosing:
!   20 <= x <= 21
integer :: xlim1, xlim2
public solve, slvtridiag
contains

subroutine init()
   real(dp) :: T
   T = Ttot_in*(1 + 0.5*(gamma - 1)*Minf**2)**(-1)
   ainf = sqrt(gamma*Rgas*T)
   Uinf = Minf*ainf
   dy0 = y(2) - y(1)
   mMinf2 = (1 - Minf**2)
   gp1m2u = (gamma + 1)*Minf**2/Uinf
   call calc_xlims
end subroutine

subroutine calc_xlims()
   logical, dimension(nx) :: mask
   mask = x .ge. 20 .and. x .le. 21
   xlim1 = minloc(x, 1, mask)
   xlim2 = maxloc(x, 1, mask)
end subroutine

subroutine apply_neumann_bc(phi)
   use grid, only: tc
   real(dp), dimension(nx,ny), intent(inout) :: phi
   real(dp) :: dydx
   integer i
   ! Symmetry BCs
   !   Upstream
   do i = 1, xlim1 - 1
       phi(i,1) = phi(i,2)
   end do
   !   Downstream
   do i = xlim2 + 1, nx
       phi(i,1) = phi(i,2)
   end do
   ! Wall BC
   ! y = tc*(-2*x**2 + 82*x - 840)
   ! dy/dx = tc*(-4*x + 82)
   do i = xlim1, xlim2
       dydx = tc*(-4*x(i) + 82)
       phi(i,1) = phi(i,2) - Uinf*dy0*dydx
   end do
end subroutine

subroutine calc_coeffs(phi, a, b, c, d, e, g)
   real(dp), dimension(nx, ny), intent(in) :: phi
   real(dp), dimension(3:nx-1, 2:ny-1), intent(out) :: a, b, c, d, e, g
   real(dp), dimension(2:nx-1, 2:ny-1) :: AA
   integer,  dimension(2:nx-1, 2:ny-1) :: mu
   real(dp) :: dx, dphi
   real(dp) :: dxo, dxm, dxp, dyo, dyp, dym
   real(dp) :: dxmo, dxmm
   real(dp) :: muAA, muAAM
   integer :: i, j
   ! Calculate A and mu
   mMinf2 = (1 - Minf**2)
   gp1m2u = (gamma + 1)*(Minf**2)/Uinf
   do j = 2, ny-1
   do i = 2, nx-1
       dx = x(i+1) - x(i-1)
       dphi = phi(i+1,j) - phi(i-1,j)
       AA(i,j) = mMinf2 - gp1m2u*(dphi/dx)
       if (AA(i,j) .ge. 0) then
           mu(i,j) = 0
       else
           mu(i,j) = 1
       end if
   end do
   end do
   ! Calculate coefficients
   do j = 2, ny-1
       dyo = y(j+1) - y(j-1)
       dym = y(j)   - y(j-1)
       dyp = y(j+1) - y(j)
       do i = 3, nx-1
           dxo = x(i+1) - x(i-1)
           dxm = x(i)   - x(i-1)
           dxp = x(i+1) - x(i)
           dxmo = x(i) - x(i-2)
           dxmm = x(i-1) - x(i-2)
           muAA = AA(i,j)*(1 - mu(i,j))
           muAAm = AA(i-1,j)*mu(i-1,j)
           a(i,j) = -2*( &
               muAA*(1/dxp + 1/dxm)/dxo &
             + (1/dyp + 1/dym)/dyo &
             - muAAm/(dxm*dxmo) )
           b(i,j) = 2/(dyo*dyp)
           c(i,j) = 2/(dyo*dym)
           d(i,j) = 2*( &
               muAA/(dxo*dxm) &
             - muAAm*(1/dxm + 1/dxmm)/dxmo )
           e(i,j) = 2*muAA/(dxo*dxp)
           g(i,j) = 2*muAAm/(dxmo*dxmm)
       end do
   end do
end subroutine

subroutine gauss_update(phi, err)
   real(dp), dimension(nx,ny), intent(inout) :: phi
   real(dp)                  , intent(out)   :: err
   real(dp), dimension(3:nx-1, 2:ny-1) :: a, b, c, d, e, g
   real(dp) :: maxphi, maxerr, corr
   integer :: i, j
   maxphi = maxval(abs(phi))
   maxerr = 0
   call calc_coeffs(phi, a, b, c, d, e, g)
   do j = 2, ny-1
       do i = 3, nx-1
           ! Calculate correction
           corr = -(1/a(i,j))*( &
               c(i,j)*phi(i,j-1) + b(i,j)*phi(i,j+1) &
             + g(i,j)*phi(i-2,j) + d(i,j)*phi(i-1,j) + e(i,j)*phi(i+1,j) )
           ! Calculate error
           err = abs(corr - phi(i,j))
           maxerr = max(maxerr, err)
           ! Assign correction
           phi(i,j) = corr
       end do
   end do
   err = maxerr!/maxphi
end subroutine

subroutine gauss_line_update(phi, err)
   real(dp), dimension(nx,ny), intent(inout) :: phi
   real(dp)                  , intent(out)   :: err
   real(dp), dimension(3:nx-1, 2:ny-1) :: a, b, c, d, e, g
   real(dp), dimension(2:ny-1) :: diag, rhs, corr
   real(dp), dimension(3:ny-1) :: l
   real(dp), dimension(2:ny-2) :: u
   real(dp) :: maxphi, maxerr
   integer :: i, j
   maxphi = maxval(phi)
   call calc_coeffs(phi, a, b, c, d, e, g)
   maxerr = 0
   do i = 3, nx-1
       ! Calculate correction for whole j line
       !   Construct diagonal and sub-diagonals
       diag = a(i,:)
       l = c(i,3:ny-1)
       u = b(i,2:ny-2)
       !   Construct RHS
       do j = 2, ny-1
           rhs(j) = -g(i,j)*phi(i-2,j) - d(i,j)*phi(i-1,j) - e(i,j)*phi(i+1,j)
       end do
       j = 2
       rhs(j) = rhs(j) - c(i,j)*phi(i,j-1)
       j = ny - 1
       rhs(j) = rhs(j) - b(i,j)*phi(i,j+1)
       !   Solve system.
       call slvtridiag(l, diag, u, rhs, corr)
       do j = 2, ny-1
           ! Calculate error
           err = abs(corr(j) - phi(i,j))
           if (i .lt. 10) write(*,*) i,j, err
           maxerr = max(maxerr, err)
           ! Assign correction
           phi(i,j) = corr(j)
       end do
   end do
end subroutine

subroutine solve(phi, r, t, iters)
   real(dp), dimension(nx,ny), intent(out) :: phi
   real(dp), dimension(:)    , intent(out) :: r, t
   integer                   , intent(out) :: iters
   integer :: k
   real(dp) :: tic, toc, err
   ! Initialize solution field
   phi = 0
   ! Set constants
   call init
   ! Start timer and solve
   call cpu_time(tic)
   do k = 1, itermax
       call apply_neumann_bc(phi)
       if (solverID .eq. 1) then
           call gauss_update(phi, err)
       else
           call gauss_line_update(phi, err)
       end if
       r(k) = err
       call cpu_time(toc)
       t(k) = toc - tic
       write(*,*) k, err
       if (err .lt. tol) then
           exit
       end if
       if (mod(k, 1000) == 0) write(*,*) k, err
   end do
   iters = k
   write(*,*) "Converged after ", iters, " iterations in ", t(iters), " seconds."
end subroutine

! Algorithm from
! http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
subroutine slvtridiag(a, b, c, d, x)
    use types, only: dp
    implicit none
    real(dp), dimension(:)        , intent(inout) :: b
    real(dp), dimension(2:size(b)), intent(in)    :: a
    real(dp), dimension(size(b)-1), intent(in)    :: c
    real(dp), dimension(size(b))  , intent(inout) :: d
    real(dp), dimension(size(b))  , intent(out)   :: x
    integer :: k, n
    real(dp) :: m
    n = size(d)
    ! Forward elimination
    do k = 2, n
        m = a(k)/b(k-1)
        b(k) = b(k) - m*c(k-1)
        d(k) = d(k) - m*d(k-1)
    end do
    ! Backward substitution
    x = d/b
    do k = n-1, 1, -1
        write(*,*) k
        x(k) = (d(k) - c(k)*x(k+1))/b(k)
    end do
end subroutine
end module
