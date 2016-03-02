program main
    use types, only: dp
    use grid, only: nx, ny, x, y
    use params, only: Minf, solverID
    use solver, only: itermax, solve
    use params
    implicit none
    real(dp), dimension(:,:), allocatable :: phi
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
    allocate(r(itermax))
    allocate(t(itermax))
    ! Run solver
    write(*,*) 'Solving'
    call solve(phi, r, t, iters)
    ! Output
    open(10, file='phi.dat')
!   do i = 1, size(phi)
!       write(10,*) phi(i)
!   end do
    open(20, file='r.dat')
    open(30, file='t.dat')
    do k = 1, iters
        write(20,*) r(k)
        write(30,*) t(k)
    end do
end program
