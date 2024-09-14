C       ===============================================================
C       Helper functions to calculate commonly required quantities
C       ===============================================================
        module common_calcs
        use types, only: dp
        implicit none
        contains

C       Calculate prim from W
        pure function calc_prim(w) result(prim)
            use constants, only: gam
            real(dp), dimension(:, :), intent(in) :: w
            real(dp), dimension(5, size(w, 2)) :: prim
            real(dp) :: rho, u, p, e
            integer :: i, n
            n = size(w, 2)
            do i = 1, n
                rho = w(1, i)
                u = w(2, i)/rho
                e = w(3, i)
                p = (gam - 1.0)*rho*(e/rho - 0.5*u**2)
                prim(1, i) = rho
                prim(2, i) = u
                prim(3, i) = p
                prim(4, i) = e
                prim(5, i) = sqrt(gam*p/rho)
            end do
        end function
C       Calculate W from prim
        pure function calc_w(prim) result(w)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)) :: w
            integer :: i, n
            n = size(prim, 2)
            do i = 1, n
                w(1, i) = prim(1, i)
                w(2, i) = prim(1, i)*prim(2, i)
                w(3, i) = prim(4, i)
            end do
        end function
C       Calculate F
        pure function calc_f(prim) result(f)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)) :: f
C           f(1) = rho*u
C           f(2) = rho*u**2 + p
C           f(3) = (e + p)*u
            integer :: i, n
            n = size(prim, 2)
            do i=1, n
                f(1, i) = prim(1, i)*prim(2, i)
                f(2, i) = prim(1, i)*prim(2, i)**2 + prim(3, i)
                f(3, i) = prim(2, i)*(prim(4, i) + prim(3, i))
            end do
        end function
C       Calculate Q
        pure function calc_q(prim, s) result(q)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(size(prim, 2)), intent(in) :: s
            real(dp), dimension(3, size(prim, 2)) :: q
            integer :: i, n
            n = size(prim, 2)
            do i=2, n
                q(1, i) = 0
                q(2, i) = prim(3, i)*(s(i) - s(i-1))
                q(3, i) = 0
            end do
        end function
C       Calculate maximum eigenvalue
        elemental function calc_lambda_max(u, c) result(lambda_max)
            real(dp), intent(in) :: u, c
            real(dp) :: lambda_max
            lambda_max = max(u, u + c, u - c)
        end function
C       Calculate residuals
        pure function calc_r(s, f_edge, q) result(r)
            real(dp), dimension(:), intent(in) :: s
            real(dp), dimension(3, size(s)-1), intent(in) :: f_edge
            real(dp), dimension(3, size(s)), intent(in) :: q
            real(dp), dimension(3, size(s)) :: r
            real(dp), dimension(size(s)-1) :: s_edge
            integer :: i, k, n
            n = size(s)
            do i = 1, n - 1
                s_edge(i) = 0.5*(s(i) + s(i + 1))
            end do
            do i = 2, n - 1
            do k = 1, 3
                r(k, i) = f_edge(k, i)*s_edge(i)
     &                    - f_edge(k, i - 1)*s_edge(i - 1)
     &                    - q(k, i)
            end do
            end do
        end function
C       Calculate error over all residuals
        pure function calc_err(r) result(err)
            real(dp), dimension(:, :), intent(in) :: r
            real(dp) :: err
            err = maxval(abs(r(1, :)))
        end function
        end module common_calcs
