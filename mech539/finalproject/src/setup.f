C       ===============================================================
C       Initialization
C       ===============================================================
        module setup
        use types, only: dp
        implicit none
        public mkgrid, init_state
        contains

C       Make initial grid.
        subroutine mkgrid(x, s)
            use inputs, only: params
            use constants, only: h, pi, t1, t2
            real(dp), dimension(:), intent(out) :: x, s
            integer :: n, i
            n = size(x)
            do i=1, n
                x(i) = (i-0.5)*params%dx
                s(i) = 1.0_dp - h*(sin(pi*x(i)**t1))**t2
            end do
        end
C       Initialize field
        subroutine init_state(prim)
            use constants, only: gam, M_in, ptot_in, ttot_in, Rgas
            use constants, only: cv
            use inputs, only: params
            real(dp), dimension(:, :), intent(out) :: prim
            real(dp) :: rho, p, t, u, c, e, M_term
            integer :: i, n
            n = size(prim, 2)
            M_term = ( 1.0 + 0.5*(gam - 1.0)*M_in**2 )
            t = ttot_in/M_term
            p = ptot_in/( M_term**(gam/(gam - 1.0)) )
            rho = p/(Rgas*t)
            c = sqrt(gam*p/rho)
            u = M_in*c
            e = rho*(0.5*u**2 + p/(rho*(gam - 1)))
            do i = 1, n
                prim(1, i) = rho
                prim(2, i) = u
                prim(3, i) = p
                prim(4, i) = e
                prim(5, i) = c
            end do
C           Impose static pressure
            do i = 2, n
                prim(3, i) = params%p_exit
                prim(1, i) = params%p_exit/(Rgas*t)
                prim(5, i) = sqrt(gam*params%p_exit/prim(1, i))
                prim(2, i) = M_in*prim(5, i)
                prim(4, i) = prim(1, i)*(cv*t + 0.5*u**2)
            end do
        end
        end module setup
