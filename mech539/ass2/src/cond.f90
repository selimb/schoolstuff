module cond
    use types, only: sp
    integer, parameter :: sp=kind(1.0_d0)
    implicit none
    real(sp), private, parameter :: AII = -4
    real(sp), private, parameter :: AIJ = 1
    real(sp), private, parameter :: BI = 0
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

    subroutine doit()
        real(sp) :: cond
    end subroutine
end module
