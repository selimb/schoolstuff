module utils
    use types
    implicit none
    contains
    subroutine write1d(a)
        real(dp), dimension(:) :: a
        integer :: n
        character(len=12) :: fmt
        n = size(a, 1)
        write (fmt, '(A1,I2,A6)') '(', n, 'F13.9)'
        write (*, fmt) a
    end subroutine
    subroutine write2d(M)
        real(dp), dimension(:, :) :: M
        integer :: i, n
        character(len=10) :: fmt
        n = size(M,1)
        do i = 1, n
            call write1d(M(i, :))
        end do
    end subroutine
end module
