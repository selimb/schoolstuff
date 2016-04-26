C       ===============================================================
C       Read input file
C       ===============================================================
        module inputs
            use types, only: dp
            use constants, only: ptot_in
            implicit none
            type Parameters
                integer :: flx_scheme
                integer :: timestep_scheme
                integer :: nx
                integer :: max_iter
                real(dp) :: p_exit_ratio
                real(dp) :: eps
                real(dp) :: tol
                real(dp) :: cfl
C               Calculated values.
                real(dp) :: dx
                real(dp) :: p_exit
            end type
            Type(Parameters), public :: params
            namelist /dat/ params
            contains

            subroutine read_input_file(filename)
                character(len=*), intent(in) :: filename
                open(10, file=filename)
                read(10, nml=dat)
                params%dx = 1.0_dp/params%nx
                params%p_exit = params%p_exit_ratio*ptot_in
            end subroutine
        end module
