subroutine post(phi, u, up, v, M, P, Cp)
    use types, only: dp
    use params, only: Minf, Pinf, gamma
    use solver, only: Uinf, ainf
    use grid, only: nx, ny, x, y
    implicit none
    real(dp), dimension(nx, ny), intent(in) :: phi
    real(dp), dimension(nx, ny), intent(out) :: u, up, v, M, P, Cp
    real(dp), dimension(nx, ny) :: a, umag
    real(dp) :: a0inf2
    integer :: i, j
    ! Calculate u
    up(1,:) = (phi(2,:) - phi(1,:))/(x(2) - x(1))
    up(nx, :) = (phi(nx,:) - phi(nx-1,:))/(x(nx) - x(nx-1))
    do i = 2, nx-1
        up(i,:) = (phi(i+1,:) - phi(i-1,:))/(x(i+1) - x(i-1))
    end do
    u = up + Uinf
    ! Calculate v = vp
    v(:,1) = (phi(:,2) - phi(:,1))/(y(2) - y(1))
    v(:,ny) = (phi(:,ny) - phi(:,ny-1))/(y(ny) - y(ny-1))
    do j = 2, ny-1
        v(:,j) = (phi(:,j+1) - phi(:,j-1))/(y(j+1) - y(j-1))
    end do
    ! Calculate velocity magnitude
    umag = sqrt(u**2 + v**2)
    ! Calculate speed of sound
    a0inf2 = ainf**2 + 0.5*(gamma-1)*Uinf**2
    a = sqrt(a0inf2 - 0.5*(gamma-1)*umag**2)
    ! Mach Number
    M = umag/a
    ! Pressure
    P = Pinf*( 1 + 0.5*(gamma-1)*(Minf**2)*(1 - (umag/Uinf)**2) )**(gamma/(gamma-1))
    ! Coefficient of Pressure
    Cp = (P/Pinf - 1)/(0.5*gamma*(Uinf/ainf)**2)
end subroutine
