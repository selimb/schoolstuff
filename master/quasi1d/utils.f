C       ===============================================================
C       Utility functions
C       ===============================================================
        module utils
        use types, only: dp
        implicit none
        contains

        subroutine write1d(x)
            real(dp), dimension(:), intent(in) :: x
            integer :: i, n
            n = size(x)
            do i = 1, n
                write(*,*) x(i)
            end do
        end subroutine

        subroutine write2d(X)
            real(dp), dimension(:, :), intent(in) :: X
            integer :: n, k
            n = size(X, 2)
            do k = 1, 3
                write(*,*) k
                call write1d(X(k, :))
            end do
            write(*, *) ''
            write(*, *) ''
            write(*, *) ''
            write (*, *) 'Iteration complete'
            write(*, *) ''
            write(*, *) ''
            write(*, *) ''
        end subroutine
        end module
