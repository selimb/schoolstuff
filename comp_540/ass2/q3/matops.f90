module matops
    use types
    implicit none
    contains
    pure subroutine swaprows(A, n, m)
        real(dp), dimension(:, :), intent(inout) :: A
        integer, intent(in) :: n, m
        real(dp), dimension(size(A, 1)) :: temp
        temp = A(n, :)
        A(n, :) = A(m, :)
        A(m, :) = temp
    end subroutine
!   Recursively swap rows given a permutation vector
    pure subroutine swaprows_recurs(A, p)
        real(dp), dimension(:, :), intent(inout) :: A
        integer, dimension(size(A,1) - 1), intent(in) :: p
        integer :: i, n
        n = size(p)
        do i = 1, n
            call swaprows(A, i, p(i))
        end do
    end subroutine
    pure function swaprows_recurs_1d(a, p) result(b)
        real(dp), dimension(:), intent(in) :: a
        integer, dimension(size(a)-1), intent(in) :: p
        real(dp), dimension(size(a)) :: b
        real(dp) :: temp
        integer :: i, n, k
        n = size(a)
        b(:) = a(:)
        do i = 1, n - 1
            k = p(i)
            temp = b(i)
            b(i) = b(k)
            b(k) = temp
        end do
    end function
!   Find row below diagonal with maximum value in given column.
    pure function find_max(A, j) result(k)
        real(dp), dimension(:, :), intent(in) :: A
        integer, intent(in) :: j
        integer :: k
        integer :: t(1)
        t = maxloc(abs(A(j:, j)))
        k = t(1) + j - 1
    end function
end module
