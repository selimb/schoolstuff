program main
    use types, only: dp
    use solver, only: slvtridiag
    implicit none
    integer, parameter :: n = 9
    real(dp), dimension(n) :: diag, rhs, corr
    real(dp), dimension(n-1) :: l, u
    diag = 20
    rhs = 0.1
    l = -10
    u = -10
!   write(*,*) diag
!   write(*,*) l
!   write(*,*) u
!   write(*,*) rhs
    call slvtridiag(l, diag, u, rhs, corr)
    write(*,*) corr

end program
