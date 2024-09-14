program main
    implicit none
    integer, parameter :: imax = 2000
    integer, parameter :: jmax = 400
    integer, parameter :: nbmax = 200
    integer :: nblocks
    integer :: i, j, k, m, l, num
    integer, dimension(nbmax) :: ni, nj, nk
    integer, dimension(nbmax) :: IDX
    real*8, dimension(3,imax,jmax,nbmax) :: X2D
    character(100) :: meshname, conname, datname
    ! Reading MBL1.CONN
    integer :: nn, ijk, irot
    integer, dimension(6, nbmax) :: itype
    integer, dimension(6) :: btype, ktype
    integer, dimension(3, 6) :: ia
    IDX = 0
    read(*,*) meshname, conname, datname
    ! Input
    open(unit=10, form='formatted', file=meshname)
    read(10,*) nblocks
    read(10,*) ( ni(m), nj(m), m = 1, nblocks )
    do  m = 1, nblocks
      read(10,*) &
        (( X2D(1,i,j,m), i=1,ni(m)), j=1,nj(m)), &
        (( X2D(2,i,j,m), i=1,ni(m)), j=1,nj(m))
    enddo

    open(unit=11, form='formatted', file=conname)
    read(11,*)
    do m = 1, nblocks
        READ  (11,501) NN,(ITYPE(K,m),BTYPE(K),(IA(L,K),L=1,3), KTYPE(K),K=1,6),IJK,IROT
    end do

    open(unit=20, form='formatted', file=datname)
    write(20,*) ' TITLE = "NASA30P30N_GRID"'
    write(20,*) ' FILETYPE = FULL'
    write(20,*) ' VARIABLES = "X" "Y" "I" "J" "BLOCK" "I1" "I2" "I3" "I4" "NI" "NJ"'
    do m = 1, nblocks
    write(20,*) ' ZONE'
    write(20,*) 'I = ', ni(m), ' J = ', nj(m), ' F=POINT'
    num = 0
    do j = 1, nj(m)
      do i = 1, ni(m)
        num = num + 1
        write(20,*) X2D(1,i,j,m), X2D(2,i,j,m), i, j, m, ITYPE(1:4,m), ni(m), nj(m)
      end do
    end do
    end do

501  FORMAT(I5,1X,6(I3,5I2,2X),I2,2X,I2)
end program

