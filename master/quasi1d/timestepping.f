C       ===============================================================
C       Time stepping
C       ===============================================================
        module timestepping
        use types, only: dp
        use inputs, only: params
        use common_calcs
        use flx_schemes, only: flx_eval
        use bc, only: update_bc
        implicit none
        private
        public timestep
        contains

        pure function calc_dt(prim) result(dt)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(size(prim, 2)) :: dt
            real(dp) :: u, c, lambda_max
            integer :: i, n
            n = size(prim, 2)
            do i = 1, n
                u = prim(2, i)
                c = prim(5, i)
                lambda_max = calc_lambda_max(u, c)
                dt(i) = params%dx*params%cfl/lambda_max
            end do
        end function
C       Euler explicit scheme
C       ---------------------------------------------------------------
        subroutine euler_xp(prim, s, r)
            real(dp), dimension(:, :), intent(inout) :: prim
            real(dp), dimension(size(prim, 2)), intent(in) :: s
            real(dp), dimension(3, size(prim, 2)), intent(out) :: r
            real(dp), dimension(3, size(prim, 2)) :: w, f, q, w_n
            real(dp), dimension(3, size(prim, 2) - 1) :: f_edge
            real(dp), dimension(size(prim, 2)) :: dt
            real(dp) :: dt_v
            integer :: i, k, n
C           Compute time step
            dt = calc_dt(prim)
            w = calc_w(prim)
            f = calc_f(prim)
            q = calc_q(prim, s)
C           Compute flux across edges
            f_edge = flx_eval(prim, w, f)
C           Compute residual
            r = calc_r(s, f_edge, q)
C           Update W
            n = size(prim, 2)
            do i = 2, n - 1
                dt_v = dt(i)/(s(i)*params%dx)
                do k = 1, 3
                    w_n(k, i) = w(k, i) - dt_v*r(k, i)
                end do
            end do
C           Update BCs
            call update_bc(prim, dt)
C           Update state vector
            prim(:, 2:n-1) = calc_prim(w_n(:, 2:n-1))
        end subroutine

C       Jameson Runge Kutta
C       ---------------------------------------------------------------
        subroutine jameson(prim, s, r_new)
            real(dp), dimension(:, :), intent(inout) :: prim
            real(dp), dimension(size(prim, 2)), intent(in) :: s
            real(dp), dimension(3, size(prim, 2)), intent(out) :: r_new
            real(dp), dimension(3, size(prim, 2)) :: w0
            real(dp), dimension(size(prim, 1), size(prim, 2)) ::
     &          prim_new
            real(dp), dimension(size(prim, 2)) :: dt, dt_dx
            integer ::  i, k, n
            n = size(prim, 2)
            dt = calc_dt(prim)
            dt_dx = dt/params%dx
            w0 = calc_w(prim)
            prim_new(:, :) = prim(:, :)
            do k = 1, 4
                call jameson_step(prim_new, k*1.0_dp)
            end do
C           Calculate residual from difference of densities
            do i = 1, n
                r_new(1, i) = (s(i)/dt_dx(i))*(
     &              prim_new(1, i) - prim(1, i)
     &          )
            end do
            call update_bc(prim, dt)
            prim(:, 2:n-1) = prim_new(:, 2:n-1)

            contains

            subroutine jameson_step(prim, k)
                real(dp), dimension(:, :), intent(inout) :: prim
                real(dp), intent(in) :: k
                real(dp), dimension(3, size(prim, 2)) :: w, f, q, r
                real(dp), dimension(3, size(prim, 2) - 1) :: f_edge
                real(dp) :: alpha
                integer :: i, j, n
                n = size(prim, 2)
                w = calc_w(prim)
                f = calc_f(prim)
                q = calc_q(prim, s)
                f_edge = flx_eval(prim, w, f)
                r = calc_r(s, f_edge, q)
                alpha = (1.0_dp)/(5.0_dp - k)
                do i = 2, n - 1
                    do j = 1, 3
                        w(j, i) = w0(j, i) - alpha*dt_dx(i)*r(j, i)
                    end do
                end do
                prim = calc_prim(w)
            end subroutine
        end subroutine

C       ---------------------------------------------------------------
C       Choose a timestepping scheme based on input
C       1 : Euler Explicit
C       2 : Jameson
        subroutine timestep(prim, s, r)
            real(dp), dimension(:, :), intent(inout) :: prim
            real(dp), dimension(size(prim, 2)), intent(in) :: s
            real(dp), dimension(3, size(prim, 2)), intent(out) :: r
            select case (params%timestep_scheme)
                case (1)
                    call euler_xp(prim, s, r)
                case (2)
                    call jameson(prim, s, r)
                case default
                    call euler_xp(prim, s, r)
            end select
        end subroutine
        end module
