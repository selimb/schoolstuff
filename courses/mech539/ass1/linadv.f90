program main
    implicit real (A-H, O-Z)
    real, dimension(:), allocatable :: u, x
    TMAX = 10
    A = 0.5
    GRID_LENGTH = 40

    dx = 1e-5
    cfl = 0.5
    dt = CFL*dx/A
    nx = (GRID_LENGTH/dx) + 1
    allocate(u(nx))
    allocate(x(nx))
    do i = 1, nx
        x(i) = (i-1)*dx
        t = tanh(250.0*(x(i) - 20.0))
        u(i) = 0.5*(1.0 + t)
    end do
    nt = TMAX/dt
    t = 0
    icnt = 1
    ievery = nt/100
    write(*,*) nt
    write(*,*) ievery
    do j = 1, nt
        t = t + dt
        u = timestep(u)
        if (mod(j,ievery) .eq. 0) then
            write(*,*) icnt
            call flush
            icnt = icnt + 1
        end if
    end do
    write(*,*) t
    open(unit=2, file='q5')
    do i = 1, nx
        write(2,*) u(i)
    end do
    contains

        function timestep(u) result (unew)
            real, dimension(:) :: u
            real, dimension(size(u)) :: unew
            unew = u
            n = size(u)
            do i = 2, n-1
                first = 0.5*A*dt/dx*(u(i+1) - u(i-1))
                second = 0.5*(A*dt/dx)**2*(u(i+1) - 2*u(i) + u(i-1))
                unew(i) = u(i) - first + second
            end do
        end function
end program
