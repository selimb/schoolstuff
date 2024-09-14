module conn
    use globals, only: num_axes, num_dirs, nblocks
    implicit none
    public
    contains

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

subroutine write_conn(iout, m)
    use globals, only: ITYPE, IA, BTYPE, KTYPE, ROT, num_dirs, num_axes
    integer, intent(in) :: iout
    integer :: m, a, d, L, ijk
    call validatehand(m, ijk)
    write (iout,501) m, &
        ((ITYPE(d, a, m), BTYPE(d, a, m), &
         (IA(L, d, a, m),L=1,3), KTYPE, &
             d=1,num_dirs), a=1,num_axes), &
        (-4, 0, 1, 2, 3, 0, d=1,2), &
        ijk, ROT
    501  FORMAT(I5,1X,6(I3,5I2,2X),I2,2X,I2)
end subroutine

function fix_conn() result(do_again)
    use globals, only: ITYPE, nblocks, lmax
    logical :: do_again
    integer :: m
    ITYPE = 0
    do_again = .false.
    do m = 1, nblocks
        do_again = (fix_block(m) .or. do_again)
    end do
end function

!
function fix_block(block) result(did_split)
    use globals, only: num_axes, num_dirs, lmax
    logical :: did_split
    integer :: block, ax, dir
    integer :: to_block, to_ax, to_dir, reversed, to_bgn, length
    integer :: num_created
    did_split = .false.
    do ax = 1, num_axes
        do dir = 1, num_dirs
            call find_exact_match(block, ax, dir, to_block, to_ax, to_dir, reversed)
            if (to_block .ne. 0) then
                call set_face_conn(block, ax, dir, to_block, to_ax, to_dir, reversed)
                call set_face_conn(to_block, to_ax, to_dir, block, ax, dir, reversed)
                cycle
            end if
            call find_longest_match(block, ax, dir, &
                to_block, to_ax, to_dir, reversed, to_bgn, length)
            if (to_block .eq. 0) then
                call bc(block, ax, dir, to_block)
                call set_face_conn(block, ax, dir, to_block, 0, 0, 0)
                cycle
            end if
            call split_block(block, ax, (/1, 1+length/), 0, num_created)
            call split_block(to_block, to_ax, (/to_bgn, to_bgn+length/), reversed, num_created)
            did_split = .true.
            return
        end do
    end do
end function

! Finds `to_block, to_ax, to_dir, reversed, to_bgn` and maximizes
! `length` such that FROM_FACE(:, 1:1+length) == TO_FACE(:, to_bgn:to_bgn+length)
! where FROM_FACE is the result of `fetch_face(block, ax, dir, 0)
! and TO_FACE is the result of `fetch_face(to_block, to_ax, to_dir, reversed)
subroutine find_longest_match(block, ax, dir, to_block, to_ax, to_dir, reversed, to_bgn, length)
    use globals, only: lmax
    integer, intent(in) :: block, ax, dir
    integer, intent(out) :: to_block, to_ax, to_dir, reversed, to_bgn, length
    real*8, dimension(2, lmax) :: FROM_FACE, TO_FACE
    integer :: i, from_length, to_length, max_length
    call fetch_face(block, ax, dir, 0, FROM_FACE, from_length)
    call find_two_point_match(block, FROM_FACE, to_block, to_ax, to_dir, reversed, to_bgn)
    if (to_block .eq. 0) return
    call fetch_face(to_block, to_ax, to_dir, reversed, TO_FACE, to_length)
    max_length = max(to_length - to_bgn, from_length - 1)
    do i = max_length, 1, -1
        if (cmp_face(FROM_FACE(:, 1:1+i), TO_FACE(:, to_bgn:to_bgn+i))) then
            length = i
            return
        end if
    end do
end subroutine

subroutine find_two_point_match(block, FACE, to_block, to_ax, to_dir, reversed, idx)
    use globals, only: lmax
    integer, intent(in) :: block
    integer, intent(out) :: to_block, to_ax, to_dir, reversed, idx
    real*8, dimension(2, lmax) :: FACE, TO_FACE
    integer :: to_length
    integer :: m, a, d, r, i
    to_block = 0
    do m = 1, nblocks
        if (m .eq. block) cycle
        do a = 1, num_axes
            do d = 1, num_dirs
                do r = 0, 1
                    if (m .gt. 200) stop
                    call fetch_face(m, a, d, r, TO_FACE, to_length)
                    do i = 1, to_length-1
                        if (cmp_face(FACE(:, 1:2), TO_FACE(:, i:i+1))) then
                            to_block = m
                            to_ax = a
                            to_dir = d
                            reversed = r
                            idx = i
                            return
                        end if
                    end do
                end do
            end do
        end do
    end do
end subroutine

subroutine fetch_face(block, ax, dir, reversed, FACE, length)
    use globals, only: X2D, ni, nj, lmax
    integer, intent(in) :: block, ax, dir, reversed
    real*8, dimension(2, lmax), intent(out) :: FACE
    integer, intent(out) :: length
    integer :: m
    integer :: i, j
    m = block
    if (ax .eq. 1) then
        length = nj(m)
        if (dir .eq. 1) then
            i = 1
        else
            i = ni(m)
        end if
        if (reversed .eq. 0) then
                FACE(:, 1:length) = X2D(:, i, 1:nj(m), m)
        else
                FACE(:, 1:length) = X2D(:, i, nj(m):1:-1, m)
        end if
    else if (ax .eq. 2) then
        length = ni(m)
        if (dir .eq. 1) then
            j = 1
        else
            j = nj(m)
        end if
        if (reversed .eq. 0) then
                FACE(1:2, 1:length) = X2D(1:2, 1:ni(m), j, m)
        else
                FACE(1:2, 1:length) = X2D(1:2, ni(m):1:-1, j, m)
        end if
    end if
end subroutine

subroutine find_exact_match(block, ax, dir, to_block, to_ax, to_dir, reversed)
    use globals, only: ni, nj, lmax, nblocks
    implicit none
    integer, intent(in) :: block, ax, dir
    integer, intent(out) :: to_block, to_ax, to_dir, reversed
    real*8, dimension(2, lmax) :: FROM_FACE, TO_FACE
    integer :: length
    integer, dimension(2) :: axes
    integer :: m, a, d, r, dummy
    call fetch_face(block, ax, dir, 0, FROM_FACE, length)
    to_block = 0
    to_ax = 0
    to_dir = 0
    reversed = 0
    do m = 1, nblocks
        axes = 0
        if (m .eq. block) cycle
        if (nj(m) .eq. length) axes(1) = 1
        if (ni(m) .eq. length) axes(2) = 1
        if (sum(axes) .eq. 0) cycle
        do a = 1, num_axes
            if (axes(a) .eq. 0) cycle
            do d = 1, num_dirs
                do r = 0, 1
                    call fetch_face(m, a, d, r, TO_FACE, dummy)
                    if (cmp_face(FROM_FACE(:, 1:length), TO_FACE(:, 1:length))) then
                        to_block = m
                        to_ax = a
                        to_dir = d
                        reversed = r
                        return
                    end if
                end do
            end do
        end do
    end do
end subroutine

subroutine set_face_conn(from_block, from_ax, from_dir, to_block, to_ax, to_dir, reversed)
    use globals, only: ITYPE, IA, BTYPE
    integer, intent(in) :: from_block, from_ax, from_dir, to_block, to_ax, to_dir, reversed
    integer, dimension(3) :: IJ
    integer :: sb, from_tangent_ax
    sb = 0
    IJ(3) = 3
    if (to_block .lt. 0) then
        if (to_block .eq. -7) sb = 1
        IJ(:) = 0
    else
        from_tangent_ax = perpendicular_ax(from_ax)
        if (from_ax .eq. to_ax) then
            IJ(1) = 1
            IJ(2) = 2
        else
            IJ(1) = 2
            IJ(2) = 1
        end if
        if (from_dir .eq. to_dir) IJ(from_ax) = -1*IJ(from_ax)
        if (reversed .eq. 1) IJ(from_tangent_ax) = -1*IJ(from_tangent_ax)
    end if
    ITYPE(from_dir, from_ax, from_block) = to_block
    IA(1:3, from_dir, from_ax, from_block) = IJ(:)
    BTYPE(from_dir, from_ax, from_block) = sb
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

! Split block `m` at indices `LOCS(1), LOCS(2)` for faces with normal axis `ax`.
! Increments nblocks by 0, 1 or 2
! Examples: Let ni(m) = 30, nj(m) = 50, nblocks0=nblocks before call to split_block.
! 1) Args: LOCS=(10, 25), ax=1, reversed=0
!    -> Split at j=(10, 25)
!       nj(nblocks0 + 1) = 10
!           X2D(:, :, 1:10, nblocks0+1) = X2D(:, :, 1:10, m)
!       nj(nblocks0 + 2) = 26 = 50 - 25 + 1
!           X2D(:, :, 1:26, nblocks0+2) = X2D(:, :, 25:50, m)
!       nj(m) = 16 = 25 - 10 + 1
!           X2D(:, :, 1:16, m) = X2D(:, :, 10:25, m)
! 2) Args: LOCS=(10, 25), ax=1, reversed=1
!    -> Split at j=(10, 25) starting from end == split at j=(26, 41)
! 3) Args: LOCS=(1, 10), ax=2, reversed=1
!    -> Split at j=(41, 50)
!       Only one block created
subroutine split_block(m, ax, in_LOCS, reversed, num_created)
    use globals, only: X2D, ni, nj
    implicit none
    integer, intent(in) :: m, ax, reversed
    integer, dimension(2), intent(in) :: in_LOCS
    integer, intent(out) :: num_created
    integer, dimension(size(in_LOCS)) :: LOCS
    integer :: lim
    if (ax .eq. 1) lim = nj(m)
    if (ax .eq. 2) lim = ni(m)
    if (reversed .ne. 0) then
        LOCS(1) = lim - in_LOCS(2) + 1
        LOCS(2) = lim - in_LOCS(1) + 1
    else
        LOCS = in_LOCS
    end if
    num_created = 0
    if (LOCS(1) .ge. LOCS(2)) then
        write(*,*) "Invalid locations for split block", LOCS
        stop
    end if
    if (LOCS(2) .gt. lim) then
        write(*,*) "Invalid split_block", m, LOCS, ax, reversed
        write(*,*) "Block limit: ", lim
        stop
    end if

    if (ax .eq. 1) then
        if (LOCS(1) .gt. 1) then
            nblocks = nblocks + 1
            nj(nblocks) = LOCS(1)
            ni(nblocks) = ni(m)
            X2D(:, :, 1:nj(nblocks), nblocks) = X2D(:, :, 1:nj(nblocks), m)
            call insert_block(nblocks, m)
            num_created = num_created + 1
        end if
        if (LOCS(2) .lt. nj(m)) then
            nblocks = nblocks + 1
            nj(nblocks) = nj(m) - LOCS(2) + 1
            ni(nblocks) = ni(m)
            X2D(:, :, 1:nj(nblocks), nblocks) = X2D(:, :, LOCS(2):nj(m), m)
            call insert_block(nblocks, m + num_created + 1)
            num_created = num_created + 1
        end if
        nj(m) = LOCS(2) - LOCS(1) + 1
        X2D(:, :, 1:nj(m), m) = X2D(:, :, LOCS(1):LOCS(2), m)
    else
        if (LOCS(1) .gt. 1) then
            nblocks = nblocks + 1
            ni(nblocks) = LOCS(1)
            nj(nblocks) = nj(m)
            X2D(:, 1:ni(nblocks), :, nblocks) = X2D(:, 1:ni(nblocks), :, m)
            call insert_block(nblocks, m)
            num_created = num_created + 1
        end if
        if (LOCS(2) .lt. ni(m)) then
            nblocks = nblocks + 1
            ni(nblocks) = ni(m) - LOCS(2) + 1
            nj(nblocks) = nj(m)
            X2D(:, 1:ni(nblocks), :, nblocks) = X2D(:, LOCS(2):ni(m), :, m)
            call insert_block(nblocks, m + num_created + 1)
            num_created = num_created + 1
        end if
        ni(m) = LOCS(2) - LOCS(1) + 1
        X2D(:, 1:ni(m), :, m) = X2D(:, LOCS(1):LOCS(2), :, m)
    end if
end subroutine

subroutine insert_block(block, pos)
    use globals, only: IDX, nbmax
    implicit none
    integer, intent(in) :: block, pos
    IDX(pos+1 : nbmax) = IDX(pos : nbmax-1)
    IDX(pos) = block
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
