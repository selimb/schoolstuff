module linsolve
    use types
    use matops
    implicit none
    private
    public lu, solve
    contains
    pure function forward_sub(L, b) result (x)
        real(dp), dimension(:), intent(in) :: b
        real(dp), dimension(size(b), size(b)), intent(in) :: L
        real(dp), dimension(size(b)) :: x
        real(dp) :: s
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
    pure function back_sub(U, b) result (x)
        real(dp), dimension(:, :), intent(in) :: U
        real(dp), dimension(:), intent(in) :: b
        real(dp), dimension(size(b)) :: x
        real(dp) :: s
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
    pure subroutine lu(A, L, U, p)
        real(dp), dimension(:,:), intent(in) :: A
        real(dp), dimension(size(A,1),size(A,1)), intent(out) :: L, U
        integer, dimension(size(A,1)), intent(out) :: p
        integer :: maxrow, i, j, k, n
        n = size(A,1)
        U(:, :) = A(:, :)
        L(:, :) = 0
        do j = 1, n - 1
            maxrow = find_max(U, j)
            p(j) = maxrow
            call swaprows(U, j, maxrow)
            call swaprows(L, j, maxrow)
            do i = j + 1, n
                L(i, j) = U(i, j)/U(j, j)
                U(i, j) = 0
                do k = j + 1, n
                    U(i, k) = U(i, k) - L(i, j)*U(j, k)
                end do
            end do
        end do
        do i = 1, n
            L(i, i) = 1
        end do
    end subroutine
    function solve(L, U, p, b) result (x)
        real(dp), dimension(:, :), intent(in) :: L
        real(dp), dimension(size(L,1), size(L,1)), intent(in) :: U
        real(dp), dimension(size(L,1)), intent(in) :: b
        integer, dimension(size(L,1)-1), intent(in) :: p
        real(dp), dimension(size(L,1)) :: x
        real(dp), dimension(size(L,1)) :: pb, y
        pb = swaprows_recurs_1d(b, p)
        y = forward_sub(L, pb)
        x = back_sub(U, y)
    end function
end module linsolve
