module solvers
    use types, only: sp
    implicit none
    real(sp), private, parameter :: AII = -4
    real(sp), private, parameter :: AIJ = 1
    real(sp), private, parameter :: BI = 0
    real(sp), public :: relax, tol
    integer, public :: nx, solverID, itermax
    public
    contains

    pure function get_row(i, j) result (row)
        integer, intent(in) :: i, j
        integer :: row
        row = j + (i-1)*nx
    end function

    pure function is_bc(i, j) result (ret)
        integer, intent(in) :: i, j
        integer :: ret
        if (i .eq. 1 .or. j .eq. 1 .or. j .eq. nx .or. i .eq. nx) then
            ret = 1
        else
            ret = 0
        end if
    end function

    subroutine update(un, maxerr)
        real(sp), dimension(:), intent(inout) :: un
        real(sp), intent(out) :: maxerr
        real(sp), allocatable, dimension(:) :: u
        real(sp) :: s, corr, err
        integer, dimension(4) :: stencil
        integer :: i, j, row
        maxerr = 0
        if (solverID .eq. 1) then
            allocate(u(size(un)))
            u = un
        end if
        do i = 1, nx
            do j = 1, nx
                row = get_row(i, j)
                if (is_bc(i,j) .eq. 1) then
                    cycle
                end if
                ! Calculate correction
                stencil = (/ row-nx, row-1, row+1, row+nx /)
                if (solverID .eq. 1) then
                    s = sum(u(stencil))
                else
                    s = sum(un(stencil))
                end if
                corr = (BI - AIJ*s)/AII
                if (solverID .eq. 3) then
                    corr = (1.0 - relax)*un(row) + relax*corr
                end if
                ! Calculate error
                err = abs(corr - un(row))
                maxerr = max(maxerr, err)
                ! Assign correction
                un(row) = corr
            end do
        end do
    end subroutine

    subroutine printu(u)
        real(sp), dimension(:) :: u
        integer :: i, j, row
        do i = 1, nx
            do j = 1, nx
                row = get_row(i, j)
                write(*,'(F6.2)',advance='no') u(row)
            end do
            write(*,*) ''
        end do
    end subroutine

    subroutine solve(u, r, t, iters)
        real(sp), dimension(:), intent(out) :: u, r, t
        integer, intent(out) :: iters
        integer :: i, top_bgn, k
        real(sp) :: tic, toc, err
        u = 0
        top_bgn = nx*(nx-1) + 1
        do i = top_bgn, size(u)
            u(i) = 1
        end do
        call cpu_time(tic)
        do k = 1, itermax
            call update(u, err)
            r(k) = err
            call cpu_time(toc)
            t(k) = toc - tic
            if (err .lt. tol) then
                write(*,*) "Converged after ", k, " iterations in ", t(k), " seconds."
                exit
            end if
        end do
        iters = k
    end subroutine
end module
