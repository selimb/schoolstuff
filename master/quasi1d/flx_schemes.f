C       ===============================================================
C       Picture of grid indices:
C             node        edge
C           |---o---|---o---|---o ... ---|---o---|
C               1   1   2   2   3       J-1  J
C       ===============================================================
C       ===============================================================
C       Schemes for flux evaluation.
C       ===============================================================
        module flx_schemes
        use types, only: dp
        use constants, only: gam
        use inputs, only: params
        use common_calcs
        implicit none
        contains

        pure function calc_lambdas(u, c) result (lambdas)
            real(dp), intent(in) :: u, c
            real(dp), dimension(3) :: lambdas
            lambdas(1) = u
            lambdas(2) = u + c
            lambdas(3) = u - c
        end function

C       Calculate the S diagonalizer and its inverse
        pure subroutine calc_S(rho, u, c, S, Sinv)
            real(dp), intent(in) :: rho, u, c
            real(dp), dimension(3, 3), intent(out) :: S, Sinv
            real(dp) :: beta, alpha
            beta = gam - 1
            alpha = 0.5*u**2
            S(:, :) = 0
            Sinv(:, :) = 0
            S(1, 1) = 1
            S(2, 1) = -u/rho
            S(3, 1) = alpha*beta
            S(2, 2) = 1/rho
            S(3, 2) = -u*beta
            S(3, 3) = beta
            Sinv(1, 1) = 1
            Sinv(2, 1) = u
            Sinv(3, 1) = alpha
            Sinv(2, 2) = rho
            Sinv(3, 2) = rho*u
            Sinv(3, 3) = 1/beta
        end subroutine

C       Calculate the C diagonalizer and its inverse
        pure subroutine calc_C(rho, u, c, CC, Cinv)
            real(dp), intent(in) :: rho, u, c
            real(dp), dimension(3, 3), intent(out) :: CC, Cinv
            real(dp) :: c2
            c2 = c**2
            CC(:, :) = 0
            Cinv(:, :) = 0
            CC(1, 1) = 1
            CC(2, 2) = rho*c
            CC(3, 2) = -rho*c
            CC(1, 3) = -1/(c2)
            CC(2, 3) = 1
            CC(3, 3) = 1
            Cinv(1, 1) = 1
            Cinv(1, 2) = 1/(2*c2)
            Cinv(2, 2) = 1/(2*rho*c)
            Cinv(3, 2) = 0.5
            Cinv(1, 3) = 1/(2*c2)
            Cinv(2, 3) = -1/(2*rho*c)
            Cinv(3, 3) = 0.5
        end subroutine

C       Calculate the diagonalizers S^-1*C^-1 and CS
        pure subroutine calc_diagonalizers(rho, u, c, SCinv, CS)
            real(dp), intent(in) :: rho, u, c
            real(dp), dimension(3, 3), intent(out) :: SCinv, CS
            real(dp), dimension(3, 3) :: S, CC, Sinv, Cinv
            call calc_C(rho, u, c, CC, Cinv)
            call calc_S(rho, u, c, S, Sinv)
            SCinv = matmul(Sinv, Cinv)
            CS = matmul(CC, S)
        end subroutine calc_diagonalizers

C       Calculate the diagonal matrices Yp and Ym according to the
C       Steger Warming scheme
        pure subroutine calc_diag_sw(u, c, Yp, Ym)
            real(dp), intent(in) :: u, c
            real(dp), dimension(3, 3), intent(out) :: Yp, Ym
            real(dp), dimension(3) :: lambdas
            real(dp) :: sq
            integer :: i
            Yp(:, :) = 0
            Ym(:, :) = 0
            lambdas = calc_lambdas(u, c)
            do i = 1, 3
                sq = sqrt(lambdas(i)**2 + params%eps**2)
                Yp(i, i) = 0.5*(lambdas(i) + sq)
                Ym(i, i) = 0.5*(lambdas(i) - sq)
            end do
        end subroutine

C       Calculate Jacobian from components
        pure function calc_jacob(SCinv, Y, CS) result (A)
            real(dp), dimension(3, 3), intent(in) :: SCinv, Y, CS
            real(dp), dimension(3, 3) :: A
            A = matmul( matmul(SCinv, Y), CS )
        end function

C       Calculate the jacobian for Steger-Warming schemes
        pure subroutine calc_jacob_sw(prim, Ap, Am)
            real(dp), dimension(:), intent(in) :: prim
            real(dp), dimension(3, 3), intent(out) :: Ap, Am
            real(dp), dimension(3, 3) :: Yp, Ym, CS, SCinv
            real(dp) :: rho, u, c
            rho = prim(1)
            u = prim(2)
            c = prim(5)
            call calc_diagonalizers(rho, u, c, SCinv, CS)
            call calc_diag_sw(u, c, Yp, Ym)
            Ap = calc_jacob(SCinv, Yp, CS)
            Am = calc_jacob(SCinv, Ym, CS)
        end subroutine
C       ---------------------------------------------------------------
C       Steger-Warming
C       ---------------------------------------------------------------
        pure function flx_sw(prim, w, f) result (f_edge)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)), intent(in) :: w, f
            real(dp), dimension(3, size(prim, 2) - 1) :: f_edge
            real(dp), dimension(3, 3, size(prim, 2)) :: Ap, Am
            integer :: i, n
            n = size(prim, 2)
            do i = 1, n
                call calc_jacob_sw(prim(:, i), Ap(:, :, i), Am(:, :, i))
            end do
            do i = 1, n - 1
                f_edge(:, i) = matmul(Ap(:, :, i), w(:, i))
     &                         + matmul(Am(:, :, i + 1), w(:, i + 1))
            end do
        end function flx_sw

C       ---------------------------------------------------------------
C       Modified Steger-Warming
C       ---------------------------------------------------------------
        pure function flx_msw(prim, w, f) result (f_edge)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)), intent(in) :: w, f
            real(dp), dimension(3, size(prim, 2) - 1) :: f_edge
            real(dp), dimension(3, 3) :: Ap_edge, Am_edge
            real(dp), dimension(5) :: prim_edge
            integer :: i, k, n, m
            n = size(prim, 2)
            m = size(prim, 1)
            do i = 1, n - 1
                do k = 1 , m
                    prim_edge(k) = 0.5*(prim(k, i) + prim(k, i + 1))
                end do
                call calc_jacob_sw(prim_edge, Ap_edge, Am_edge)
                f_edge(:, i) = matmul(Ap_edge, w(:, i))
     &                         + matmul(Am_edge, w(:, i + 1))
            end do
        end function flx_msw

C       ---------------------------------------------------------------
C       Corrected-Modified Steger-Warming
C       ---------------------------------------------------------------
        pure function flx_cmsw(prim, w, f) result(f_edge)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)), intent(in) :: w, f
            real(dp), dimension(3, size(prim, 2) - 1) :: f_edge
            real(dp), dimension(3, size(prim, 2) - 1) :: f_msw, f_sw
            real(dp) :: omega, p, pp, dp_dx
            integer :: i, k, n
            f_msw = flx_msw(prim, w, f)
            f_sw = flx_sw(prim, w, f)
            n = size(prim, 2)
            do i = 1, n - 1
                p = prim(3, i)
                pp = prim(3, i + 1)
                dp_dx = (pp - p)/min(p, pp)
                omega = 1/(1 + dp_dx**2)
                do k = 1, 3
                    f_edge(k, i) = omega*f_msw(k, i)
     &                             + (1 - omega)*f_sw(k, i)
                end do
            end do
        end function flx_cmsw

C       ---------------------------------------------------------------
C       Roe
C       ---------------------------------------------------------------
        function flx_roe(prim, w, f) result (f_edge)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)), intent(in) :: w, f
            real(dp), dimension(3, size(prim, 2) - 1) :: f_edge
            real(dp), dimension(3, 3) :: Y, CS, SCinv, A
            real(dp) :: rho, rhoh, rhop, sqrho, sqrhop, u, up,
     &          e, ep, p, pp, c, cp, ch, uh, hh, eps
            real(dp), dimension(3) :: lam, lamp, lamh
            integer :: i, k, n
            n = size(prim, 2)
            do i = 1, n - 1
                rho = prim(1, i)
                rhop = prim(1, i + 1)
                sqrho = sqrt(rho)
                sqrhop = sqrt(rhop)
                u = prim(2, i)
                up = prim(2, i + 1)
                p = prim(3, i)
                pp = prim(3, i + 1)
                e = prim(4, i)
                ep = prim(4, i + 1)
                c = prim(5, i)
                cp = prim(5, i + 1)
                rhoh = sqrho*sqrhop
                uh = (sqrho*u + sqrhop*up)/(sqrho + sqrhop)
                hh = ( sqrho*((e + p)/rho) + sqrhop*((ep + pp)/rhop))
     &               /(sqrho + sqrhop)
C               TODO Maybe this is not squared
                ch = sqrt((gam - 1)*(hh - 0.5*uh**2))
                call calc_diagonalizers(rhoh, uh, ch, SCinv, CS)
C               Calculate eigenvalues
                Y(:, :) = 0
                lam = calc_lambdas(u, c)
                lamp = calc_lambdas(up, cp)
                lamh = calc_lambdas(uh, ch)
                do k = 1, 3
                    eps = max(0.0_dp, lamh(k)-lam(k), lamp(k)-lamh(k))
                    if (abs(lamh(k)) <= eps) then
                        lamh(k) = 0.5*(lamh(k)**2/eps + eps)
                    end if
                    Y(k, k) = lamh(k)
                end do
C               Calculate Jacobian and Flux
C               A = abs(calc_jacob(SCinv, Y, CS))
                A = Y
                f_edge(:, i) = 0.5*( f(:, i) + f(:, i + 1)
     &                               - matmul(A, w(:, i+1) - w(:, i))
     &          )
C               if (i < 5) then
C                   write(*,*) 'i = '
C                   write(*,*) i
C                   write(*,*) 'rhoH'
C                   write(*,*) rhoh
C                   write(*,*) 'uH'
C                   write(*,*) uh
C                   write(*,*) 'hh'
C                   write(*,*) hh
C                   write(*,*) 'ch'
C                   write(*,*) ch
C                   write(*,*) 'lamh'
C                   write(*,*) lamh
C                   write(*,*) 'Y'
C                   write(*,*) Y
C                   write(*,*) 'A'
C                   write(*,*) A
C                   write(*,*) 'w'
C                   write(*,*) w(:, i)
C                   write(*,*) 'f'
C                   write(*,*) f(:, i)
C                   write(*,*) 'f_edge'
C                   write(*,*) f_edge(:, i)
C               end if
            end do
        end function

C       ---------------------------------------------------------------
C       Scalar Dissipation
C       ---------------------------------------------------------------
        pure function flx_scalar(prim, w, f) result (f_edge)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)), intent(in) :: w, f
            real(dp), dimension(3, size(prim, 2) - 1) :: f_edge
            real(dp) :: u_avg, c_avg, lambda_max
            integer :: i, k, n
            n = size(prim, 2)
            do i = 1, n - 1
                u_avg = 0.5*(prim(2, i) + prim(2, i+1))
                c_avg = 0.5*(prim(5, i) + prim(5, i+1))
                lambda_max = calc_lambda_max(u_avg, c_avg)
                do k = 1, 3
                    f_edge(k, i) = 0.5*(
     &                  f(k, i) + f(k, i+1)
     &                  - params%eps*lambda_max*(w(k, i+1) - w(k, i))
     &              )
                end do
            end do
        end function flx_scalar
C       Choose a flux evaluation scheme based on input
C       1 : Scalar dissipation
C       2 : Steger Warming
C       3 : Modified Steger Warming
C       4 : Roe
        function flx_eval(prim, w, f) result (f_edge)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)), intent(in) :: w, f
            real(dp), dimension(3, size(prim, 2) - 1) :: f_edge
            select case (params%flx_scheme)
                case (1)
                    f_edge = flx_scalar(prim, w, f)
                case (2)
                    f_edge = flx_sw(prim, w, f)
                case (3)
                    f_edge = flx_msw(prim, w, f)
                case (4)
                    f_edge = flx_cmsw(prim, w, f)
                case (5)
                    f_edge = flx_roe(prim, w, f)
                case default
                    write(*,*) 'Invalid Scheme'
                    f_edge = flx_scalar(prim, w, f)
            end select
        end function
        end module flx_schemes
