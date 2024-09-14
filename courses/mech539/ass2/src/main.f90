program main
    use types, only: dp
    use solvers
    implicit none
    real(dp), dimension(:), allocatable :: u, r, t, z
    integer :: i, k, iters, l
    integer :: do_cond
    real(dp) :: cond
    ! Read input
    open(100, file='input.prm')
    read(100,*) nx
    read(100,*) solverID
    read(100,*) relax
    read(100,*) do_cond
    ! Set other parameters
    l = nx*nx
    itermax = 1e9
    if (do_cond .eq. 1) then
        tol = 1e-15
    else
        tol = 1e-6
    end if
    ! Allocate
    allocate(u(nx*nx))
    allocate(r(itermax))
    allocate(t(itermax))
    if (do_cond .eq. 1) allocate(z(size(u)))
    ! Run solver
    write(*,*) 'Solving'
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
    ! Condition Number if applicable
    if (do_cond .eq. 1) then
        call solvez(u, z)
        cond = norm2(z)/norm2(u)/EPSILON(1.0_dp)
        open(50, file='cond.dat')
        write(50,*) cond
        write(*,*) "Conditoning number:", cond
    end if

end program
