program main
    use types, only: sp
    use solvers
    implicit none
    real(sp), dimension(:), allocatable :: u, r, t
    integer :: i, k, iters, l
    ! Read input
    open(100, file='input.prm')
    read(100,*) nx
    read(100,*) solverID
    read(100,*) relax
    ! Set other parameters
    l = nx*nx
    itermax = 1e9
    tol = 1e-6
    ! Allocate
    allocate(u(nx*nx))
    allocate(r(itermax))
    allocate(t(itermax))
    ! Run solver
    call solve(u, r, t, iters)
    ! Output
    open(10, file='u.dat')
    do i = 1, size(u)
        write(10,*) u(i)
    end do
    open(20, file='r.dat')
    open(30, file='t.dat')
    do k = 1, iters
        write(20,*) r(k)
        write(30,*) t(k)
    end do

end program
