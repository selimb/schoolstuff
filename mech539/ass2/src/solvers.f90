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

    pure function get_stencil(row) result (stencil)
        integer, intent(in) :: row
        integer, dimension(4) :: stencil
        stencil(1) = row - nx
        stencil(2) = row - 1
        stencil(3) = row + 1
        stencil(4) = row + nx
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

    pure function jacobi(u, un, row) result (ui)
        real(sp), dimension(:), intent(in) :: u, un
        integer, intent(in) :: row
        real(sp) :: ui
        integer, dimension(4) :: stencil
        real(sp) :: s
        integer :: i
        stencil = get_stencil(row)
        s = 0
        do i = 1, size(stencil)
            s = s + u(stencil(i))
        end do
        ui = (BI - AIJ*s)/AII
    end function

    pure function gauss(u, un, row) result (ui)
        real(sp), dimension(:), intent(in) :: u, un
        integer, intent(in) :: row
        real(sp) :: ui
        integer, dimension(4) :: stencil
        real(sp) :: s
        integer :: i
        stencil = get_stencil(row)
        s = 0
        do i = 1, 2
            s = s + un(stencil(i))
            s = s + u(stencil(i+2))
        end do
        ui = (BI - AIJ*s)/AII
    end function

    pure function sor(u, un, row) result (ui)
        real(sp), dimension(:), intent(in) :: u, un
        integer, intent(in) :: row
        real(sp) :: ui
        real(sp) :: ugs
        ugs = gauss(u, un, row)
        ui = (1.0 - relax)*u(row) + relax*ugs
    end function

    function solver_step(u, un, row) result (ui)
        real(sp), dimension(:), intent(in) :: u, un
        integer, intent(in) :: row
        real(sp) :: ui
        select case (solverID)
            case(1)
                ui = jacobi(u, un, row)
            case(2)
                ui = gauss(u, un, row)
            case(3)
                ui = sor(u, un, row)
        end select
    end function

    function update(u) result (un)
        real(sp), dimension(:), intent(in) :: u
        real(sp), dimension(size(u)) :: un
        integer :: i, j, row
        un = u
        do i = 1, nx
            do j = 1, nx
                row = get_row(i, j)
                if (is_bc(i,j) .eq. 1) then
                    cycle
                end if
                un(row) = solver_step(u, un, row)
            end do
        end do
    end function

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

    subroutine solve(un, r, t, iters)
        real(sp), dimension(:), intent(out) :: un, r, t
        integer, intent(out) :: iters
        real(sp), dimension(size(un)) :: u
        real(sp) :: err
        integer :: i, top_bgn, k
        real(sp) :: tic, toc
        u = 0
        top_bgn = nx*(nx-1) + 1
        do i = top_bgn, size(u)
            u(i) = 1
        end do
        call cpu_time(tic)
        do k = 1, itermax
            un = update(u)
            err = norm2(u - un)
            r(k) = err
            call cpu_time(toc)
            t(k) = toc - tic
            if (err .lt. tol) then
                write(*,*) "Converged after ", k, " iterations."
                exit
            end if
            u = un
        end do
        iters = k
    end subroutine
end module
