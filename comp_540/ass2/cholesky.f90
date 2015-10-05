module types
    integer, parameter :: sp = kind(1.0)
end module
module globals
    use types
    real(sp), dimension(:, :), allocatable :: A, L
end module
module setup
    use types
    implicit none
    private
    public init
    contains
    subroutine read_test()
        use globals, only: A
        open(20, file='test')
        allocate(A(3, 3))
        read(20, *) A
        A = transpose(A)
    end subroutine
    subroutine hilbert(n)
        use globals, only: A
        integer :: n
        integer :: i, j
        allocate(A(n, n))
        forall(i=1:n, j=1:n) &
            A(i, j) = 1.0_sp/(i + j - 1.0_sp)
    end subroutine
    subroutine init()
        integer :: n
        write(*,*) 'Enter integer value for n -- -1 for test'
        read(*, *) n
        if (n == -1) then
            call read_test
        else
            call hilbert(n)
        end if
    end subroutine
end module
module utils
    use types
    implicit none
    contains
    subroutine write1d(a)
        real(sp), dimension(:) :: a
        integer :: n
        character(len=10) :: fmt
        n = size(a, 1)
        write (fmt, '(A1,I1,A6)') '(', n, 'F13.9)'
        write (*, fmt) a
    end subroutine
    subroutine write2d(M)
        real(sp), dimension(:, :) :: M
        integer :: i, n
        character(len=10) :: fmt
        n = size(M,1)
        do i = 1, n
            call write1d(M(i, :))
        end do
    end subroutine
    function mysqrt(x)
        use globals, only: L
        real(sp), intent(in) :: x
        real(sp) :: mysqrt
        if (x < 0 ) then
            write(*,*) 'Error square-rooting: ', x
            write(*,*) 'L:'
            call write2d(L)
            stop
        end if
        mysqrt = sqrt(x)
    end function
    function mydiv(a, b)
        use globals, only: L
        real(sp), intent(in) :: a, b
        real(sp) :: mydiv
        if (abs(b) < tiny(0.0_sp)) then
            write(*,*) 'Divide by zero'
            call write2d(L)
            stop
        end if
        mydiv = a/b
    end function
!   Mostly used for testing purposes
    function normF(M)
        real(sp), dimension(:, :) :: M
        real(sp) :: normF
        integer :: i, n
        n = size(M,1)
        normF = 0
        do i = 1, n
            normF = normF + norm2(M(i, :))
        end do
        normF = sqrt(normF)
    end function normF
end module utils
module linearsolve
    use types
    implicit none
    contains
    function forward_sub(L, b) result (x)
        real(sp), dimension(:), intent(in) :: b
        real(sp), dimension(size(b), size(b)), intent(in) :: L
        real(sp), dimension(size(b)) :: x
        real(sp) :: s
        integer :: i, j, n
        n = size(b)
        do i = 1, n
            s = b(i)
            do j = 1, i - 1
                s = s - L(i, j)*x(j)
            end do
            x(i) = s/L(i, i)
        end do
    end function
    function back_sub(U, b) result (x)
        real(sp), dimension(:, :), intent(in) :: U
        real(sp), dimension(:), intent(in) :: b
        real(sp), dimension(size(b)) :: x
        real(sp) :: s
        integer :: i, j, n
        n = size(b)
        do i = n, 1, -1
            s = b(i)
            do j = i + 1, n
                s = s - U(i, j)*x(j)
            end do
            x(i) = s/U(i, i)
        end do
    end function
end module linearsolve
! Main Program
program main
    use types
    use setup
    use utils
    use globals
    use linearsolve
    implicit none
    integer :: i, j, n
    real(sp), dimension(:, :), allocatable :: A2
    real(sp), dimension(:), allocatable :: x, b, e, y
    real(sp) :: rel_err, rel_res, rel_mat_res, frob_a, frob_adiff
    logical :: res
    call init
    write(*,*) 'A'
    call write2d(A)
    n = size(A,1)
    allocate(L(n, n))
    allocate(A2(n, n))
    allocate(x(n))
    allocate(b(n))
    allocate(e(n))
    allocate(y(n))
    e(:) = 1
    b = matmul(A, e)
    L(:,:) = 0
    do j = 1, n
        L(j, j) = mysqrt( &
            A(j, j) &
            - dot_product( L(j,1:j-1), L(j, 1:j-1) ) &
        )
        do i = j + 1, n
            L(i, j) = mydiv( &
                A(i, j) - dot_product ( L(i, 1:j-1), L(j, 1:j-1) ), &
                L(j, j) &
            )
        end do
    end do
    A2 = matmul(L, transpose(L))
    write(*,*) 'L'
    call write2d(L)
    ! write(*,*) 'L x LT'
    ! call write2d(A2)
    ! write(*,*) 'L x LT = A'
    ! res = all(A2 == A)
    ! write(*,*) res
    y = forward_sub(L, b)
    x = back_sub(transpose(L), y)
    write(*,*) 'x'
    call write1d(x)
    write(*,*) 'Ax'
    call write1d(matmul(A, x))
    ! Compute Errors
    rel_err = norm2(x - e)/norm2(x)
    frob_a = normF(A)
    rel_res = norm2(b - matmul(A, x))/(frob_a*norm2(x))
    frob_adiff = normF(A - A2)
    rel_mat_res = frob_adiff/frob_a
    write(*,*) 'Relative Error:'
    write(*,*) rel_err
    write(*,*) 'Relative Residual:'
    write(*,*) rel_res
    write(*,*) 'Relative Matrix Residual'
    write(*,*) rel_mat_res
end program
