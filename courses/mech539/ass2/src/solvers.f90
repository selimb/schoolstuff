module solvers
    use types, only: dp
    implicit none
    real(dp), private, parameter :: AII = -4
    real(dp), private, parameter :: AIJ = 1
    real(dp), private, parameter :: BI = 0
    real(dp), public :: relax, tol
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
        logical :: ret
        if (i .eq. 1 .or. j .eq. 1 .or. j .eq. nx .or. i .eq. nx) then
            ret = .true.
        else
            ret = .false.
        end if
    end function

    subroutine update(un, maxerr)
        real(dp), dimension(:), intent(inout) :: un
        real(dp), intent(out) :: maxerr
        real(dp), allocatable, dimension(:) :: u
        real(dp) :: s, corr, err
        integer, dimension(4) :: stencil
        integer :: i, j, row
        maxerr = 0
        if (solverID .eq. 1) then
            allocate(u(size(un)))
            u = un
        end if
        do i = 2, nx-1
            do j = 2, nx-1
                row = get_row(i, j)
                ! Calculate correction
                stencil = (/ row-nx, row-1, row+1, row+nx /)
                if (solverID .eq. 1) then
                    s = sum(u(stencil))
                else
                    s = sum(un(stencil))
                end if
                corr = (BI - AIJ*s)/AII
                if (solverID .eq. 3) then
                    corr = (1.0_dp - relax)*un(row) + relax*corr
                end if
                ! Calculate error
                err = abs(corr - un(row))
                maxerr = max(maxerr, err)
                ! Assign correction
                un(row) = corr
            end do
        end do
    end subroutine

    subroutine solvez(u, z)
        real(dp), dimension(:), intent(in) :: u
        real(dp), dimension(size(u)), intent(out) :: z
        real(dp) :: maxerr
        real(dp) :: r, s, corr, err
        integer, dimension(4) :: stencil
        integer :: i, j, k, row
        maxerr = 0
        ! Initialize z to 1's everywhere except on boundary
        z = 1
        do i = 1, nx
            do j = 1, nx
                row  = get_row(i,j)
                if (is_bc(i,j)) z(row) = 0
            end do
        end do
        do k = 1, itermax
            maxerr = 0
            do i = 2, nx-1
                do j = 2, nx-1
                    row = get_row(i,j)
                    ! Calculate correction
                    stencil = (/ row-nx, row-1, row+1, row+nx /)
                    r = BI - (AIJ*sum(u(stencil)) + AII*u(row))
                    s = sum(z(stencil))
                    corr = (r - AIJ*s)/AII
                    corr = (1.0_dp - relax)*z(row) + relax*corr
                    ! Calculate error
                    err = abs(corr - z(row))
                    maxerr = max(maxerr, err)
                    ! Assign correction
                    z(row) = corr
                end do
            end do
            if (maxerr .lt. tol) then
                write(*,*) "Converged conditioning after ", k, " iterations"
                exit
            end if
        end do
!       call printu(z)
    end subroutine

    subroutine printu(u)
        real(dp), dimension(:) :: u
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
        real(dp), dimension(:), intent(out) :: u, r, t
        integer, intent(out) :: iters
        integer :: i, top_bgn, k
        real(dp) :: tic, toc, err
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
