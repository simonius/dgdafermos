! This module contains initial conditions and (one) reference solution for our testcases. 
! All tests are initial conditions for the Euler equations.

module testcases
        use param
        use trianggrid
        use trianghp
        use euler

        implicit none

        contains

        
        ! This function returns the one-dimensional Lax-Shock tube initial condition
        ! depending on a given x coordinate. The y moment is set to zero.
        function Laxu0(x) result (u0)
                real(np), intent(in) :: x
                real(np) :: u0(4)

                if (x .LE. 0.5_np) then
                        u0 = ConsVar(1.0_np, 0.0_np, 0.0_np, 1.0_np)
                else
                        u0 = ConsVar(0.125_np, 0.0_np, 0.0_np, 0.1_np)
                end if
        end function Laxu0

        ! A two-dimensional version of the Lax Shock Tube, 
        ! assuming symmetry in the y direction.
        ! The complete tensor u will be initialized according to grid gr,
        ! if this function is called.
        ! gr: the used grid
        ! rt: the reference triangle
        ! u: the array of conserved variables
        subroutine LaxInit(gr, rt, u)
                type(Grid) :: gr
                type(nprt) :: rt
                real(np), intent(out) :: u(:, :, :)

                integer :: tr, nd
                real(np) :: x(2)

                do tr=1,size(u, 1)
                        do nd=1, size(u, 2)
                                x = GetNodeKord(gr, RT, tr, nd)
                                u(tr, nd, :) = Laxu0(x(1)) 
                        end do
                end do
        end subroutine LaxInit

        ! A radial two-dimensional version of the Lax Shock Tube
        ! The Sedov/Kurganov blast wave
        ! gr: the used grid
        ! u: the array of conserved variables
        subroutine KSInit(gr, rt, u)
                type(Grid) :: gr
                type(nprt) :: rt
                real(np), intent(out) :: u(:, :, :)

                integer :: tr, nd
                real(np) :: x(2)
                real(np) :: xcent(2) = [0.0D0, 0.0D0]

                do tr=1,size(u, 1)
                        do nd=1, size(u, 2)
                                x = GetNodeKord(gr, RT, tr, nd)
                                if (norm2(x-xcent) < 0.08_np) then
                                        u(tr, nd, :) = ConsVar(1.0_np, 0.0_np, 0.0_np, 1.0_np)
                                else
                                        u(tr, nd, :) = ConsVar(0.125_np, 0.0_np, 0.0_np, 0.1_np)
                                end if 
                        end do
                end do
        end subroutine KSInit

        ! Friedrichs smooth bump function, has value 1 at r=0
        ! and 0 at r .GE. 1
        function FB(r)
                real(np) :: FB
                real(np) :: r
                if (r .GE. 1.0_np) then
                        FB = 0.0_np
                else
                        FB = exp(-1/(1-r**2))*exp(1.0_np)
                end if
                FB = FB**6
                !FB = exp(-10*r**2)
        end function FB

        ! A moving bump, used for exactness tests.
        ! gr: the used grid
        ! rt: the reference triangle
        ! u: the array of conserved variables
        subroutine MBInit(gr, rt, u)
                type(Grid) :: gr
                type(nprt) :: rt
                double precision, intent(out) :: u(:, :, :)

                integer :: tr, nd
                double precision, dimension(2) :: x
                double precision, dimension(2) :: xcent = [-1.0_np, -0.0_np]
                double precision :: rholoc

                do tr=1,size(u, 1)
                        do nd=1, size(u, 2)
                                x = GetNodeKord(gr, RT, tr, nd)
                                rholoc = FB(norm2(x-xcent)*3.0) + 1.0D0
                                u(tr, nd, :) = ConsVar(rholoc, 1.0D0, 0.0D0, 1.0D0)
                        end do
                end do
        end subroutine MBInit

        ! This subroutine writes the reference solution to MBInit at time t
        ! into the array u.
        subroutine MBref(gr, rt, u, t)
                type(Grid) :: gr
                type(nprt) :: rt
                real(np), intent(out) :: u(:, :, :)
                real(np) :: t

                integer :: tr, nd
                real(np) :: x(2)
                real(np) :: xcent(2) = [-1.0_np, -0.0_np]
                real(np) :: rholoc

                xcent(1) = xcent(1) + t

                do tr=1,size(u, 1)
                        do nd=1, size(u, 2)
                                x = GetNodeKord(gr, RT, tr, nd)
                                rholoc = FB(norm2(x-xcent)*3.0_np) + 1.0_np
                                u(tr, nd, :) = ConsVar(rholoc, 1.0_np, 0.0_np, 1.0_np)
                        end do
                end do
        end subroutine MBref

        subroutine KHI(gr, rt, u, t)
                type(Grid) :: gr
                type(nprt) :: rt
                real(np), intent(out) :: u(:, :, :)
                real(np) :: t
        end subroutine KHI

        ! A homogeneous initial condition with the value from the boundary function bf.
        ! Used for the forward facing step and airfoil testcases.
        subroutine HomInit(gr, rt, u, bf)
                type(Grid) :: gr
                type(nprt) :: rt
                real(np), intent(out) :: u(:, :, :)
                external :: bf
                integer :: tr, nd
                real(np) :: x(2)
                real(np) :: uval(size(u, 3))

                do tr=1, size(u, 1)
                        do nd=1, size(u, 2)
                                x = GetNodeKord(gr, rt, tr, nd)
                                call bf(uval, uval, x, 0.0_np)
                                u(tr, nd, :) = uval
                        end do
                end do
        end subroutine HomInit
end module testcases
