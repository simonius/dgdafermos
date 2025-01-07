! The external parameters of the solver
module param
        use iso_fortran_env
        implicit none
        logical, parameter :: debug = .false.                   ! Print debug information
        integer, parameter :: hp=real128                        ! The high precision floating points
        integer, parameter :: np=kind(1.0d0)                    ! The normal precision floating points
        real(hp), parameter :: pihp = 4.0_hp*ATAN(1.0_hp)       ! a PI constant with high precision
        real(np), parameter :: pinp = 4.0_np*ATAN(1.0_np)       ! a PI constant with normal precision
end module param
