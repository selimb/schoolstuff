module setup
    use types
    implicit none
    private
    public init
    contains
    subroutine read_input(A)
        real(dp), dimension(:, :), allocatable, intent(out) :: A
        integer :: n = 4
        open(20, file='test')
        allocate(A(n, n))
        read(20, *) A
        A = transpose(A)
    end subroutine
    subroutine mkA0(A, n)
        real(dp), dimension(:, :), allocatable, intent(out) :: A
        integer, intent(in) :: n
        integer :: i, j
        allocate(A(n, n))
        A(:, :) = 0
        do i = 1, n
        do j = 1, n
            if (i == j) then
                A(i, j) = 1
            else if (i > j) then
                A(i, j) = -1
            else if (j == n) then
                A(i, j) = 1
            end if
        end do
        end do
    end subroutine
    subroutine init(A)
        use globals, only: dotest
        real(dp), dimension(:, :), allocatable :: A
        integer :: i, n = 30
        write(*,*) 'Do verification?'
        read(*, *) dotest
        if (dotest == 1) then
            call read_input(A)
        else
            call mkA0(A, 30)
            do i = 1, n
                A(i, i) = A(i, i) + 1E-8
            end do
        end if
    end subroutine
end module
