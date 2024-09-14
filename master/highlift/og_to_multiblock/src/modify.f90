program main
    use globals, only: ni, nj, X2D, ITYPE, IA, BTYPE, nblocks, IDX ,nk, X3D
    use conn
    implicit none
    integer :: nblocks0
    integer :: i, j, m, block, n
    integer :: k, num
    logical :: do_fix
    real*8 :: bgn, fin
    ! Initialize
    IDX = 0
    ITYPE = 0
    IA = 0
    BTYPE = -1

    ! Input
    open(unit=10, form='formatted', file='nasa30p30n-multiblock.plot3d.x')
    read(10,*) nblocks0
    nblocks = nblocks0
    read(10,*) ( ni(m), nj(m), m = 1, nblocks )
    do  m = 1, nblocks
      read(10,*) &
          (( X2D(1,i,j,m), i=1,ni(m)), j=1,nj(m)), &
          (( X2D(2,i,j,m), i=1,ni(m)), j=1,nj(m))
    enddo
    do m = -10, nblocks0
        IDX(m) = m
    end do

    ! Split blocks to fix connectivity
    do_fix = .true.
    do while (do_fix)
        do_fix = fix_conn()
    end do

    ! Split blocks into smaller chunks
    call split_block(15, 2, (/ 101, 201 /), 0, m)
    write(*,*) m
    do_fix = .true.
    do while (do_fix)
        do_fix = fix_conn()
    end do
    call split_block(18, 1, (/ 101, 201 /), 0, m)
    write(*,*) m
    do_fix = .true.
    do while (do_fix)
        do_fix = fix_conn()
    end do
    call split_block(70, 1, (/ 1, 101 /), 0, m)
    do_fix = .true.
    do while (do_fix)
        do_fix = fix_conn()
    end do
    call split_block(8, 1, (/ 1, 137 /), 0, m)
    do_fix = .true.
    do while (do_fix)
        do_fix = fix_conn()
    end do

    ! Output
    open(unit=20, form='formatted', file='modded2D.x')
    write(20,*) nblocks
    do block = 1, nblocks
!      m = IDX(block)
       m = block
       write(20,*) ni(m), nj(m)
    end do
    do block = 1, nblocks
!       m = IDX(block)
        m = block
        do n = 1, 2
            write(20,*) (( X2D(n,i,j,m), i=1,ni(m)), j=1,nj(m))
        end do
    end do

    bgn = 0
    fin = 100
    num = 3
    open(unit=30, form='formatted', file='modded3D.x')
    open(unit=31, form='formatted', file='modded.CONN')
    write(30,*) nblocks
    write(31,*) nblocks
    do block = 1, nblocks
        m = IDX(block)
        m = block
        write(30,*) ni(m), nj(m), num
    end do
    do block = 1, nblocks
        m = IDX(block)
        m = block
        call extrude(m, bgn, fin, num)
        do n = 1, 3
            write(30,*) ((( X3D(n,i,j,k), i=1,ni(m)), j=1,nj(m)), k=1,nk(m))
        end do
        call write_conn(31, m)
    end do

end program
