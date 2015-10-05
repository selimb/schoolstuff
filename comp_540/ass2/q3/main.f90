program main
    use globals, only: dotest
    use types
    use setup
    use utils
    use linsolve
    implicit none
    real(dp), dimension(:, :), allocatable :: A, L, U
    real(dp), dimension(:), allocatable :: b, x, r, d, e, y
    integer, dimension(:), allocatable :: p
    real(dp), parameter :: tol = 10E-15_dp
    real(dp) :: growth, err
    integer :: iter, n
    call init(A)
    n = size(A,1)
    allocate(L(n, n), U(n, n))
    allocate(b(n), x(n), r(n), d(n), e(n), y(n))
    allocate(p(n - 1))
    e(:) = 1.0_dp
    b = matmul(A, e)
    call lu(A, L, U, p)
    growth = maxval(abs(U))/maxval(abs(A))
    x = solve(L, U, p, b)
    if (dotest == 1) then
        write(*,*) 'L'
        call write2d(L)
        write(*,*) 'U'
        call write2d(U)
        write(*,*) 'p'
        write(*,*) p
        write(*,*) 'x'
        call write1d(x)
        stop
    end if
    write(*,*) ''
    write(*,*) 'Question c)'
    write(*,*) '-----------'
    write(*,*) 'Growth Factor'
    write(*,*) growth
    write(*,*) '||x - xc|| / ||xc||'
    write(*,*) maxval(abs(x-e))/maxval(abs(x))
    write(*,*) ''
    write(*,*) 'Question d)'
    write(*,*) '-----------'
    write(*,*) 'Beginning Iterative Refinement'
    iter = 1
    err = 1
    do while (err > tol)
        write(*,*) ''
        write(*, '(A,I3)') 'Refinement ', iter
        r = b - matmul(A, x)
        d = solve(L, U, p, r)
        x = x + d
        err = norm2(d)/norm2(x)
        write(*,*) '||x - xc|| / ||xc||'
        write(*,*) maxval(abs(x-e))/maxval(abs(x))
        write(*,*) '||d|| / ||x||'
        write(*,*) err
        iter = iter + 1
        if (iter > 100) then
            write(*,*) 'Maximum number of iterations reached'
            stop
        end if
    end do
    write(*,*) 'Final x'
    write(*, '(F20.16)') x
end program
