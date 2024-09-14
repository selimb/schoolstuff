C       ===============================================================
C       Constants
C       ===============================================================
        module constants
            use types, only: dp
            implicit none
C           Math
            real(dp), public, parameter :: pi = 4.0_dp*atan(1.0_dp)
C           Grid constants
            real(dp), public, parameter :: h = 0.15_dp
            real(dp), public, parameter :: t1 = 0.8_dp
            real(dp), public, parameter :: t2 = 3
C           Fluid constants
            real(dp), public, parameter :: gam = 1.4_dp
            real(dp), public, parameter :: Rgas = 1716
            real(dp), public, parameter :: cv = Rgas/(gam - 1)
C           BCs
            real(dp), public, parameter :: ttot_in = 531.2_dp
            real(dp), public, parameter :: ptot_in = 2117
            real(dp), public, parameter :: M_in = 0.2_dp
C           Used in inlet BC
            real(dp), public, parameter ::
     &          astar = 2*gam*((gam - 1)/(gam + 1))*cv*ttot_in
        end module
