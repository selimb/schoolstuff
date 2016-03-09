program main
    use types, only: dp
    use grid, only: nx, ny, x, y
    use params, only: Minf, solverID
    use solver, only: itermax, solve
    use params
    implicit none
    real(dp), dimension(:,:), allocatable :: phi, u, up, v, M, P, Cp
    real(dp), dimension(:), allocatable :: r, t
    integer :: i, j, k, iters
    ! Read input
    open(100, file='input.prm')
    read(100,*) Minf
    read(100,*) solverID
    open(101, file='grid')
    read(101,*) nx, ny
    allocate(x(nx))
    allocate(y(ny))
    read(101,*) (x(i), i=1,nx)
    read(101,*) (y(j), j=1,ny)
    ! Allocate
    allocate(phi(nx,ny))
    allocate(u(nx,ny))
    allocate(up(nx,ny))
    allocate(v(nx,ny))
    allocate(M(nx,ny))
    allocate(P(nx,ny))
    allocate(Cp(nx,ny))
    allocate(r(itermax))
    allocate(t(itermax))
    ! Run solver
    write(*,*) 'Solving'
    call solve(phi, r, t, iters)
    ! Post-Process
    call post(phi, u, up, v, M, P, Cp)
    ! Output
    open(10, file='phi.dat')
    open(11, file='M.dat')
    open(12, file='P.dat')
    open(13, file='Cp.dat')
    open(14, file='u.dat')
    open(15, file='up.dat')
    open(16, file='v.dat')
    do i = 1, nx
        write(10,*) (phi(i,j), j=1,ny)
        write(11,*) (M(i,j), j=1,ny)
        write(12,*) (P(i,j), j=1,ny)
        write(13,*) (Cp(i,j), j=1,ny)
        write(14,*) (u(i,j), j=1,ny)
        write(15,*) (up(i,j), j=1,ny)
        write(16,*) (v(i,j), j=1,ny)
    end do
    open(20, file='r.dat')
    open(30, file='t.dat')
    do k = 1, iters
        write(20,*) r(k)
        write(30,*) t(k)
    end do
end program
