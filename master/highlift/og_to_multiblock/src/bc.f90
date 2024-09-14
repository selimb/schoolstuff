subroutine bc(block, ax, dir, to_block)
    use conn, only: fetch_face
    use globals, only: lmax
    implicit none
    integer, intent(in) :: block, ax, dir
    integer, intent(out) :: to_block
    integer :: length
    real*8, dimension(2, lmax) :: FACE
    real*8, dimension(2) :: point
    real*8 :: x, y
    call fetch_face(block, ax, dir, 0, FACE, length)
    point = FACE(:, 1)
    x = point(1)
    y = point(2)
    if (x .gt. -5 .and. x .lt. 25 .and. &
        y .gt. -5 .and. y .lt. 5) then
        to_block = -7
    else
        to_block = -5
    end if
end subroutine
