module utils
    use globals, only: num_axes, num_dirs, nblocks
    implicit none
    public
    integer, parameter :: GET = 2, &
                          SET = 3
    contains

subroutine read_conn()
    use globals, only: ITYPE, IA, BTYPE, num_dirs, num_axes
    integer :: m, block, a, d, L
    integer, parameter :: iin = 11
    integer :: ktype
    open(unit=iin, form='formatted', file='modded.CONN')
    read(iin,*) ! nblocks
    do block = 1, nblocks
        m = block
        read(iin,501) m, &
            ((ITYPE(d, a, m), BTYPE(d, a, m), &
             (IA(L, d, a, m),L=1,3), ktype, &
                 d=1,num_dirs), a=1,num_axes)
    end do
    501  FORMAT(I5,1X,4(I3,5I2,2X))
end subroutine

function conn_intact()
    use globals, only: ITYPE, nblocks, lmax, num_axes, num_dirs, lmax
    logical :: conn_intact
    integer :: m, nb, ax, dir, to_ax, to_dir, from_l, to_l, reversed
    real*8, dimension(2, lmax) :: from_face, to_face
    real*8 :: diff
    conn_intact = .true.
    do m = 1, nblocks
        do ax = 1, num_axes
            do dir = 1, num_dirs
                nb = ITYPE(dir, ax, m)
                if (nb.lt.0) cycle
                call get_nb_csys(m, ax, dir, to_ax, to_dir, reversed)
                call fetch_face(m, ax, dir, 0, from_face, from_l, GET)
                call fetch_face(nb, to_ax, to_dir, reversed, to_face, to_l, GET)
                if (from_l .ne. to_l) then
                    write(*,*) m, ax, dir
                    write(*,*) nb, to_ax, to_dir, reversed
                    write(*,*) "length mismatch", from_l, to_l
                    conn_intact = .false.
                    cycle
                end if
                diff = maxval(abs(from_face(:,1:from_l) - to_face(:,1:to_l)))
                if (diff > 1E-14) then
                    write(*,*) m, ax, dir
                    write(*,*) nb, to_ax, to_dir, reversed
                    write(*,*) "face not same", diff
                    conn_intact = .false.
                end if
            end do
        end do
    end do
end function

subroutine output(suff)
    use globals, only: ni, nj, nblocks
    character(1) :: suff
    integer :: m
    integer, dimension(2) :: maxlev
    integer, parameter :: iout = 14
    open(unit=iout, form='formatted', file='output'//suff)
    write(iout, '(A)', advance='no') ""
    do m = 1, nblocks
        call cmp_maxlev(m, maxlev)
        if (maxlev(1).lt.3 .or. maxlev(2).lt.3) then
            write(iout, *) "Block:", m
            write(iout, *) "Dims:", ni(m), nj(m)
            write(iout, *) "Maxlev:", maxlev(:)
        end if
    end do
    close(iout)
end subroutine

subroutine cmp_maxlev(block, maxlev)
    use globals, only: ni, nj
    integer, intent(in) :: block
    integer, dimension(2), intent(out) :: maxlev
    integer :: ilev, jlev, sub
    ilev = 1
    sub = ni(block)
    do while (mod(sub - 1, 2).eq.0 .and. sub.gt.2)
        sub = (sub - 1)/2 + 1
        ilev = ilev + 1
    end do
    jlev = 1
    sub = nj(block)
    do while (mod(sub - 1, 2).eq.0 .and. sub.gt.2)
        sub = (sub - 1)/2 + 1
        jlev = jlev + 1
    end do
    maxlev = (/ ilev, jlev /)
end subroutine

subroutine extend
    use globals, only: ITYPE, nblocks, num_axes, num_dirs, lmax
    integer :: m, ax, dir, length, n
    real*8 :: xavg
    real*8, dimension(2, lmax) :: face
    do m = 1, nblocks
        do ax = 1, num_axes
            do dir = 1, num_dirs
                if (ITYPE(dir, ax, m).ne.-5) cycle
                call fetch_face(m, ax, dir, 0, face, length, GET)
                xavg = sum(face(1,1:length))/length
                if (xavg.gt.148) cycle
                do n = 1, 3
                    call extend_1(m, ax, dir)
                end do
            end do
        end do
    end do
end subroutine

subroutine extend_1(m, ax, dir)
    use globals, only: ni, nj, X2D
    integer, intent(in) :: m, ax, dir
    integer :: i, j, l
    real*8, dimension(2) :: offset
    if (ax .eq. 1) then
        do j = 1, nj(m)
            if (dir .eq. 1) then
                offset = X2D(:, 1, j, m) - X2D(:, 2, j, m)
                X2D(:, 2:ni(m)+1, j, m) = X2D(:, 1:ni(m), j, m)
                do l = 1, 2
                    X2D(l, 1, j, m) = &
                        X2D(l, 2, j, m) + offset(l)
                end do
            else
                offset = X2D(:, ni(m), j, m) - X2D(:, ni(m)-1, j, m)
                do l = 1, 2
                    X2D(l, ni(m)+1, j, m) = &
                        X2D(l, ni(m), j, m) + offset(l)
                end do
            end if
        end do
        ni(m) = ni(m) + 1
    else
        do i = 1, ni(m)
            if (dir .eq. 1) then
                offset = X2D(:, i, 1, m) - X2D(:, i, 2, m)
                X2D(:, i, 2:nj(m)+1, m) = X2D(:, i, 1:nj(m), m)
                do l = 1, 2
                    X2D(l, i, 1, m) = &
                        X2D(l, i, 2, m) + offset(l)
                end do
            else
                offset = X2D(:, i, nj(m), m) - X2D(:, i, nj(m)-1, m)
                do l = 1, 2
                    X2D(l, i, nj(m)+1, m) = &
                        X2D(l, i, nj(m), m) + offset(l)
                end do
            end if
        end do
        nj(m) = nj(m) + 1
    end if
end subroutine

recursive subroutine dup_propagate(m_in, ax_in, dir_in)
    use globals, only: ITYPE, NBL
    integer, intent(in) :: m_in, ax_in, dir_in
    integer :: from_ax, to_ax, from_dir, to_dir, reversed
    integer :: nb, nb_ax, nb_dir
    integer :: num0, num1, nb_num
    num0 = NBL(ax_in, m_in)
    call dup_1(m_in, ax_in, dir_in)
    num1 = NBL(ax_in, m_in)
    if (num1 - num0 .ne. 1) then
        write(*,*) "ASSERT"
        write(*,*) num1, num0
        stop
    end if
    from_ax = perpendicular_ax(ax_in)
    do from_dir = 1, 2
        nb = ITYPE(from_dir, from_ax, m_in)
        if (nb .lt. 0) cycle
        call get_nb_csys(m_in, from_ax, from_dir, to_ax, to_dir, reversed)
        if (reversed .eq. 1) then
            nb_dir = perpendicular_ax(dir_in)
        else
            nb_dir = dir_in
        end if
        nb_ax = perpendicular_ax(to_ax)
        nb_num = NBL(nb_ax, nb)
        if (nb_num .eq. num1) then
            ! write(*,*) "Already done", nb
            cycle
        else if (nb_num .ne. num0) then
            write(*,*) "Matched wrong ax"
            write(*,*) m_in, from_ax, from_dir
            write(*,*) nb, to_ax, to_dir, reversed
            write(*,*) nb_ax
            stop
        end if
        call dup_propagate(nb, nb_ax, nb_dir)
    end do
end subroutine

subroutine dup_1(m, ax, dir)
    use globals, only: ni, nj, X2D, ITYPE
    integer, intent(in) :: m, ax, dir
    integer :: i, j, to_ax, to_dir, reversed, nb, from_dir, l
    real*8, dimension(2) :: offset
    if (ax .eq. 1) then
        do j = 1, nj(m)
            if (dir .eq. 1) then
                offset = X2D(:, 2, j, m) - X2D(:, 1, j, m)
                do l = 1, 2
                    X2D(l, 3:ni(m)+1, j, m) = &
                        X2D(l, 2:ni(m), j, m) + offset(l)
                end do
            else
                offset = X2D(:, ni(m)-1, j, m) - X2D(:, ni(m), j, m)
                X2D(:, 2:ni(m)+1, j, m) = X2D(:, 1:ni(m), j, m)
                do l = 1, 2
                    X2D(l, ni(m)-1:1:-1, j, m) = &
                        X2D(l, ni(m):2:-1, j, m) + offset(l)
                end do
            end if
        end do
        ni(m) = ni(m) + 1
    else
        do i = 1, ni(m)
            if (dir .eq. 1) then
                offset = X2D(:, i, 2, m) - X2D(:, i, 1, m)
                do l = 1, 2
                    X2D(l, i, 3:nj(m)+1, m) = &
                        X2D(l, i, 2:nj(m), m) + offset(l)
                end do
            else
                offset = X2D(:, i, nj(m)-1, m) - X2D(:, i, nj(m), m)
                X2D(:, i, 2:nj(m)+1, m) = X2D(:, i, 1:nj(m), m)
                do l = 1, 2
                    X2D(l, i, nj(m)-1:1:-1, m) = &
                        X2D(l, i, nj(m):2:-1, m) + offset(l)
                end do
            end if
        end do
        nj(m) = nj(m) + 1
    end if
    from_dir = perpendicular_ax(dir)
    nb = ITYPE(from_dir, ax, m)
    if (nb .lt. 0) then
        write(*,*) "CANNOT EDIT WALL"
        stop
    end if
    call sync(m, nb)
    call get_nb_csys(m, ax, from_dir, to_ax, to_dir, reversed)
    ! call fetch_face(m, ax, from_dir, 0, from_face, from_l, GET)
    ! call fetch_face(nb, to_ax, to_dir, reversed, from_face, from_l, SET)
    call merge_2(nb, to_ax, to_dir)
end subroutine

subroutine sync(m, nb)
    use globals, only: num_axes, num_dirs, ITYPE, lmax
    integer, intent(in) :: m, nb
    integer :: ax, dir, to_ax, to_dir, reversed, from_l
    real*8, dimension(2, lmax) :: from_face
    logical :: found
    found = .false.
    do ax = 1, num_axes
        do dir = 1, num_dirs
            if (ITYPE(dir, ax, m) .ne. nb) cycle
            if (found) then
                write(*,*) "ASSERTION. FOUND MULTIPLE SYNC"
                write(*,*) m, nb
                stop
            end if
            found = .true.
            call get_nb_csys(m, ax, dir, to_ax, to_dir, reversed)
            call fetch_face(m, ax, dir, 0, from_face, from_l, GET)
            call fetch_face(nb, to_ax, to_dir, reversed, from_face, from_l, SET)
        end do
    end do
    if (.not. found) then
        write(*,*) "COULD NOT SYNC"
        stop
    end if
end subroutine

subroutine merge_2(m, ax, dir)
    use globals, only: ni, nj, X2D
    integer, intent(in) :: m, ax, dir
    real*8, dimension(2) :: avg
    integer :: i, j
    if (ax .eq. 1) then
        do j = 1, nj(m)
            if (dir .eq. 1) then
                avg = 0.5*(X2D(:, 3, j, m) + X2D(:, 1, j, m))
                X2D(:, 2, j, m) = avg
            else
                avg = 0.5*(X2D(:, ni(m)-2, j, m) + X2D(:, ni(m), j, m))
                X2D(:, ni(m)-1, j, m) = avg
            end if
        end do
    else
        do i = 1, ni(m)
            if (dir .eq. 1) then
                avg = 0.5*(X2D(:, i, 3, m) + X2D(:, i, 1, m))
                X2D(:, i, 2, m) = avg
            else
                avg = 0.5*(X2D(:, i, nj(m)-2, m) + X2D(:, i, nj(m), m))
                X2D(:, i, nj(m)-1, m) = avg
            end if
        end do
    end if
end subroutine

subroutine rem_1(m, ax, idx)
    use globals, only: ni, nj, X2D
    integer, intent(in) :: m, ax, idx
    if (ax .eq. 1) then
        X2D(:, idx:ni(m)-1, :, m) = X2D(:, idx+1:ni(m), :, m)
        ni(m) = ni(m) - 1
    else
        X2D(:, :, idx:nj(m)-1, m) = X2D(:, :, idx+1:nj(m), m)
        nj(m) = nj(m) - 1
    end if
end subroutine

subroutine resize_blocks()
    use globals, only: nblocks, num_axes, NBL
    integer :: m, ax
    do m = 1, nblocks
        do ax = 1, num_axes
            if (NBL(ax, m) .eq. 102) then
                call resize_block(m, 101, ax)
            end if
        end do
    end do
end subroutine

subroutine resize_block(m, num, ax)
    use globals, only: X2D, ni, nj
    integer, intent(in) :: m, num, ax
    integer :: i, j
    real*8, dimension(2) :: bgn, fin
    if (ax .eq. 1) then
        do j = 1, nj(m)
            bgn = X2D(:, 1, j, m)
            fin = X2D(:, ni(m), j, m)
            X2D(:, 1:num, j, m) = linspace_2d(bgn, fin, num)
        end do
        ni(m) = num
    else
        do i = 1, ni(m)
            bgn = X2D(:, i, 1, m)
            fin = X2D(:, i, nj(m), m)
            X2D(:, i, 1:num, m) = linspace_2d(bgn, fin, num)
        end do
        nj(m) = num
    end if
end subroutine

function linspace_2d(bgn, fin, num)
    real*8, dimension(2), intent(in) :: bgn, fin
    integer, intent(in) :: num
    real*8, dimension(2, num) :: linspace_2d
    real*8, dimension(2) :: step, val
    integer :: n
    step = (fin - bgn)/(num - 1)
    val = bgn
    linspace_2d(:, 1) = bgn
    linspace_2d(:, num) = fin
    do n = 1, num - 2
        val = val + step
        linspace_2d(:, n+1) = val
    end do
    val = val + step
    if (maxval(abs(val - linspace_2d(:, num))/abs(val)) .gt. 1E-10) then
        write(*,*) "PROBLEM", val, linspace_2d(:, num)
        stop
    end if
end function

subroutine grow_blocks(ax_in, dir_in, idx_in, inc, blocks)
    use globals, only: ITYPE
    integer, intent(in) :: ax_in, dir_in, idx_in, inc
    integer, dimension(:), intent(in) :: blocks
    integer :: ax0, dir0, idx0, ax, dir, idx, m, m0, reversed
    integer :: n, num_blocks
    integer :: nb
    call grow_block(blocks(1), ax_in, inc, idx_in)
    ax0 = ax_in
    idx0 = idx_in
    dir0 = dir_in
    num_blocks = size(blocks)
    do n = 2, num_blocks
        m0 = blocks(n-1)
        m = blocks(n)
        call get_nb_csys(m0, ax0, dir0, ax, dir, reversed)
        nb = ITYPE(dir0, ax0, m0)
        if (nb .ne. m) then
            write(*,*) "WRONG NEIGHBOR", m0, ax0, dir0
            write(*,*) m, nb
            stop
        end if
        if (reversed .eq. 1) then
            idx = -idx0
        else
            idx = idx0
        end if
        call grow_block(m, ax, inc, idx)
        ax0 = ax
        dir0 = perpendicular_ax(dir)
        idx0 = idx
    end do
end subroutine

subroutine grow_block(m, ax, inc, idx)
    use globals, only: X2D, ni, nj
    integer, intent(in) :: m, ax, inc, idx
    real*8, dimension(:, :), allocatable :: line_old, line_new
    integer :: num0, num, i, j, s, pos
    if (inc.le.0) then
        write(*,*) "must be positive growth", inc
        stop
    end if
    if (ax.eq.2) then
        num0 = ni(m)
        ni(m) = num0 + inc
    else
        num0 = nj(m)
        nj(m) = num0 + inc
    end if
    num = num0 + inc
    pos = idx
    allocate(line_old(2, num))
    allocate(line_new(2, num))
    if (ax.eq.2) then
        do j = 1, nj(m)
            do s = num0, num - 1
                if (idx .lt. 0) pos = s - abs(idx) + 1
                line_old(:, 1:s) = X2D(:, 1:s, j, m)
                call add_one(s, line_old, line_new(:, 1:s+1), pos)
                X2D(:, 1:s+1, j, m) = line_new(:, 1:s+1)
            end do
        end do
    else
        do i = 1, ni(m)
            do s = num0, num - 1
                if (idx .lt. 0) pos = s - abs(idx) + 1
                line_old(:, 1:s) = X2D(:, i, 1:s, m)
                call add_one(s, line_old, line_new(:, 1:s+1), pos)
                X2D(:, i, 1:s+1, m) = line_new(:, 1:s+1)
            end do
        end do
    end if
end subroutine

subroutine tmp
    real*8, dimension(2, 4) :: line
    real*8, dimension(2, 5) :: line_new
    line(1,:) = (/0.0, 1.0, 2.0, 3.0/)
    line(2,:) = (/0.0, 1.0, 2.0, 3.0/)
    write(*,*) "og", line(1,:)
    call add_one(4, line, line_new, 1)
    write(*,*) "new", line_new(1,:)
    call add_one(4, line, line_new, 2)
    write(*,*) "new", line_new(1,:)
    call add_one(4, line, line_new, 4 - 1 + 1)
    write(*,*) "new", line_new(1,:)

    write(*,*) "og", line(1,4:1:-1)
    call add_one(4, line(:,4:1:-1), line_new, 4 - 1 + 1)
    write(*,*) "new", line_new(1,5:1:-1)
end subroutine

subroutine add_one(length, line_old, line_new, pos)
    integer, intent(in) :: length, pos
    real*8, dimension(2, length), intent(in) :: line_old
    real*8, dimension(2, length+1), intent(out) :: line_new
    real*8 :: ratio
    real*8, dimension(2) :: diff, remain
    integer :: n, l
    ratio = 1.0d0*length/(length + 1.0d0)
    line_new = 0
    line_new(:, 1:length) = line_old(:, :)
    do n = 2, length
        diff = line_old(:, n) - line_old(:, n-1)
        line_new(:, n) = line_new(:, n-1) + diff*ratio
    end do
    remain = line_old(:, length) - line_new(:, length)
    do l = 1, 2
        line_new(l, pos+1:length+1) = line_new(l, pos:length) + remain(l)
    end do
!   line_new(:, 1) = line_old(:, 1)
    line_new(:, length+1) = line_old(:, length)
!   do n = 2, length
!       diff = line_old(:, n) - line_old(:, n-1)
!       if (n .gt. pos) then
!           new_loc = n+pos
!       else
!           new_loc = n
!       end if
!       line_new(:, new_loc) = line_new(:, new_loc-1) + diff*ratio
!   end do
end subroutine

subroutine extrude(block, bgn, fin, num)
    use globals, only: X2D, X3D, nk, ni
    integer, intent(in) :: block
    real*8, intent(in) :: bgn, fin
    integer, intent(in) :: num
    real*8, dimension(num) :: z
    integer :: k, m, n
    real*8 :: step
    step = (fin - bgn)/(num - 1)
    z = (/ ((n*step), n = 0, num - 1) /)
    m = block
    nk(m) = num
    do k = 1, num
        X3D(1:2,1:ni(m),:,k) = X2D(1:2,1:ni(m),:,m)
        X3D(3,:,:,k) = z(k)
    end do
end subroutine

subroutine fetch_face(block, ax, dir, reversed, FACE, length, stat)
    use globals, only: X2D, ni, nj, lmax
    integer, intent(in) :: block, ax, dir, reversed, stat
    real*8, dimension(2, lmax), intent(inout) :: FACE
    integer, intent(inout) :: length
    integer :: m
    integer :: i, j
    m = block
    if (ax .eq. 1) then
        if (stat .eq. GET) length = nj(m)
        if (dir .eq. 1) then
            i = 1
        else
            i = ni(m)
        end if
        if (reversed .eq. 0) then
            if (stat .eq. GET) then
                FACE(:, 1:length) = X2D(:, i, 1:nj(m), m)
            else
                X2D(:, i, 1:nj(m), m) = FACE(:, 1:length)
            end if
        else
            if (stat .eq. GET) then
                FACE(:, 1:length) = X2D(:, i, nj(m):1:-1, m)
            else
                X2D(:, i, nj(m):1:-1, m) = FACE(:, 1:length)
            end if
        end if
    else if (ax .eq. 2) then
        if (stat .eq. GET) length = ni(m)
        if (dir .eq. 1) then
            j = 1
        else
            j = nj(m)
        end if
        if (reversed .eq. 0) then
            if (stat .eq. GET) then
                FACE(1:2, 1:length) = X2D(1:2, 1:ni(m), j, m)
            else
                X2D(:, 1:ni(m), j, m) = FACE(:, 1:length)
            end if
        else
            if (stat .eq. GET) then
                FACE(1:2, 1:length) = X2D(1:2, ni(m):1:-1, j, m)
            else
                X2D(:, ni(m):1:-1, j, m) = FACE(:, 1:length)
            end if
        end if
    end if
end subroutine

function cmp_face(FROM_FACE, TO_FACE) result(same)
    real*8, dimension(:, :), intent(in) :: FROM_FACE
    real*8, dimension(size(FROM_FACE, 1), size(FROM_FACE, 2)), intent(in) :: TO_FACE
    logical :: same
    real*8, dimension(size(FROM_FACE, 1), size(FROM_FACE, 2)) :: diff
    real*8, parameter :: tol = 1e-15
    diff = abs(FROM_FACE - TO_FACE)
    same = .false.
    if (maxval(diff) < tol) same = .true.
end function

subroutine get_nb_csys(m, ax, dir, to_ax, to_dir, reversed)
    use globals, only: IA
    integer, intent(in) :: m, ax, dir
    integer, intent(out) :: to_ax, to_dir, reversed
    integer, dimension(3) :: IJ
    IJ = IA(:, dir, ax, m)
    if (abs(IJ(1)) .eq. 1) then
        to_ax = ax
    else
        to_ax = perpendicular_ax(ax)
    end if
    if (IJ(ax) .lt. 0) then
        to_dir = dir
    else
        to_dir = perpendicular_ax(dir)
    end if
    if (IJ(perpendicular_ax(ax)).lt.0) then
        reversed = 1
    else
        reversed = 0
    end if
end subroutine

function perpendicular_ax(ax) result(ax1)
    integer, intent(in) :: ax
    integer :: ax1
    if (ax .eq. 1) then
        ax1 = 2
    else
        ax1 = 1
    end if
end function

end module
