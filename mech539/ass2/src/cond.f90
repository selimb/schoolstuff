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

    function update(u) result (un)
        real(sp), dimension(:), intent(in) :: u
        real(sp), dimension(size(u)) :: un
        integer :: i, j, row
        un = u
        do i = 1, nx
            do j = 1, nx
                row = get_row(i,j)
                ! Compute r. b = 0 unless at top
                ! and Ax0 is 0 unless last 2 rows

                r = b - Ax
                ! Compute ui
            end do
        end do
    subroutine doit()
        real(sp) :: cond

        ! Make b

    end subroutine
end module
