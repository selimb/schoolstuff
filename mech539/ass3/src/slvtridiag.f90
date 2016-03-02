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
    do k = 1, n-1
        x(k) = (d(k) - c(k)*x(k+1))/b(k)
    end do
end subroutine
