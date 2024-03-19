! This module contains initial conditions and (one) reference solution for our testcases. 
! All tests are initial conditions for the Euler equations.

module testcases
        use triangles
        use euler

        implicit none

        double precision, parameter :: pie = 4*ATAN(1.0D0)
        contains

        
        ! This function returns the one-dimensional Lax-Shock tube initial condition
        ! depending on a given x coordinate. The y moment is set to zero.
        function Laxu0(x) result (u0)
                double precision, intent(in) :: x
                double precision :: u0(4)

                if (x .LE. 0.5D0) then
                        u0 = ConsVar(1.0D0, 0.0D0, 0.0D0, 1.0D0)
                else
                        u0 = ConsVar(0.125D0, 0.0D0, 0.0D0, 0.1D0)
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
                type(RefTriangle) :: rt
                double precision, intent(out) :: u(:, :, :)

                integer :: tr, nd
                double precision, dimension(2) :: x

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
                type(RefTriangle) :: rt
                double precision, intent(out) :: u(:, :, :)

                integer :: tr, nd
                double precision, dimension(2) :: x
                double precision, dimension(2) :: xcent = [0.0D0, 0.0D0]

                do tr=1,size(u, 1)
                        do nd=1, size(u, 2)
                                x = GetNodeKord(gr, RT, tr, nd)
                                if (norm2(x-xcent) < 0.8D-1) then
                                        u(tr, nd, :) = ConsVar(1.0D0, 0.0D0, 0.0D0, 1.0D0)
                                else
                                        u(tr, nd, :) = ConsVar(0.125D0, 0.0D0, 0.0D0, 0.1D0)
                                end if 
                        end do
                end do
        end subroutine KSInit

        ! Friedrichs smooth bump function, has value 1 at r=0
        ! and 0 at r .GE. 1
        function FB(r)
                double precision :: FB
                double precision :: r
                if (r .GE. 1.0) then
                        FB = 0.D0
                else
                        FB = exp(-1/(1-r**2))*exp(1.0D0)
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
                type(RefTriangle) :: rt
                double precision, intent(out) :: u(:, :, :)

                integer :: tr, nd
                double precision, dimension(2) :: x
                double precision, dimension(2) :: xcent = [-1.0D0, -0.0D0]
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
                type(RefTriangle) :: rt
                double precision, intent(out) :: u(:, :, :)
                double precision :: t

                integer :: tr, nd
                double precision, dimension(2) :: x
                double precision, dimension(2) :: xcent = [-1.0D0, -0.0D0]
                double precision :: rholoc

                xcent(1) = xcent(1) + t

                do tr=1,size(u, 1)
                        do nd=1, size(u, 2)
                                x = GetNodeKord(gr, RT, tr, nd)
                                rholoc = FB(norm2(x-xcent)*3.0) + 1.0D0
                                u(tr, nd, :) = ConsVar(rholoc, 1.0D0, 0.0D0, 1.0D0)
                        end do
                end do
        end subroutine MBref

        ! A homogeneous initial condition with the value from the boundary function bf.
        ! Used for the forward facing step and airfoil testcases.
        subroutine HomInit(gr, rt, u, bf)
                type(Grid) :: gr
                type(RefTriangle) :: rt
                double precision, intent(out) :: u(:, :, :)
                external :: bf
                integer :: tr, nd
                double precision, dimension(2) :: x
                double precision, dimension(size(u, 3)) :: uval

                do tr=1, size(u, 1)
                        do nd=1, size(u, 2)
                                x = GetNodeKord(gr, rt, tr, nd)
                                call bf(uval, uval, x, 0.0D0)
                                u(tr, nd, :) = uval
                        end do
                end do
        end subroutine HomInit
end module testcases
