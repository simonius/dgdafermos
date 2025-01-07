! This program is a proof of concept implementation of a unstructured DG
! solver for the Euler equations implementing the entropy rate criterion.
! It loads a triangular grid in the form of two filenames corresponding to 
! a grid in the file format of Skewshuks Triangle program and solves the 
! Euler equations on it up to the given time tend. Afterwards it
! saves the grid and solution and exits.
!       ./twdrp nodefile elefile outputname tend 

! we advide tend = 0.8 on a [1 x 1] domain

! This version is specialized on tests that starts with a two-dimensional riemann problem


! The boundary function for the free flow around the domain.
! uo: resulting outer value of the conserved variables
! u: inner value of the conserved variables
! x: position
! t: time.
subroutine twdrpbf(uo, u, x, t)
        use euler
        real(np), intent(in) :: u(4)
        real(np), intent(in) :: x(2)
        real(np), intent(in) :: t
        real(np), intent(out) :: uo(4)
        real(np), parameter :: x0 = 0.3, y0 = 0.3 ! The initial discontinuity position.
        real(np) :: u0NE(4), u0NW(4), u0SW(4), u0SE(4), fd(4)
        real(np) :: St, Sl, Sb, Sr              ! The speeds of the corresponding shocks.

        u0NE = consvar(1.5_np, 0.0_np, 0.0_np, 1.5_np)
        u0NW = ConsVar(0.5322580645_np, 1.2060453_np, 0.0_np, 0.3_np)
        u0SW = ConsVar(0.13799283_np, 1.206045378_np, 1.2060453783_np, 0.029032258_np)
        u0SE = ConsVar(0.5322580645_np, 0.0_np, 1.2060453_np, 0.3_np)
        

        fd = fEuler(u0NW) - fEuler(u0NE)
        St = fd(1) / (u0NW(1) - u0NE(1))
        fd = gEuler(u0SW) - gEuler(u0NW)
        Sl = fd(1) / (u0SW(1) - u0NW(1))
        fd = fEuler(u0SE) - fEuler(u0SW)
        Sb = fd(1) / (u0SE(1) - u0SW(1))
        fd = gEuler(u0NE) - gEuler(u0SE)
        Sr = fd(1) / (u0NE(1) - u0SE(1))
        ! North-East
        if (x(1) > x0 + St*t .AND. x(2) > y0 + Sr*t) then
                uo = consvar(1.5_np, 0.0_np, 0.0_np, 1.5_np)
        elseif (x(1) < x0 + St*t .AND. x(2) > y0 + Sl*t) then ! North-West
                uo = ConsVar(0.5322580645_np, 1.2060453_np, 0.0_np, 0.3_np)
        elseif (x(1) < x0 + Sb*t .AND. x(2) < y0 + Sl*t) then ! Sout-West
                uo = ConsVar(0.13799283_np, 1.206045378_np, 1.2060453783_np, 0.029032258_np) 
        else    ! Sout-East
                uo = ConsVar(0.5322580645_np, 0.0_np, 1.2060453_np, 0.3_np)
        end if
end subroutine twdrpbf

program twdrp
        use sdgsolver
        use trianggrid
        use trianghp
        use testcases
        use opshp
        implicit none

        type(Grid) :: gr
        type(nprt) :: rt
        double precision, allocatable :: u0(:, :, :)
        double precision, allocatable :: ures(:, :, :)
        double precision :: tend, gs, cfl, fs
        integer :: Nt, ord
        external :: twdrpbf
        character(len=256) :: nfname, efname, ofname, ordstring, tmaxstring

        if (command_argument_count() .NE. 5) then
                write (*,*) "Missing input filenames!"
                stop
        end if

        ! Parsing of the command line arguments
        call get_command_argument(1, nfname)
        call get_command_argument(2, efname)
        call get_command_argument(3, ofname)
        call get_command_argument(4, ordstring)
        call get_command_argument(5, tmaxstring)

        read (ordstring, *) ord
        read (tmaxstring, *) tend
        
        cfl = 0.3
        cfl = cfl / (dble(ord)**2 + dble(ord))

        ! Initialization of the reference triangle
        call DGinit(rt, ord)
        
        ! Load the grid from the HDD 
        call GridTriangleImp(gr, nfname, efname)

        ! Allocation of arrays for the unknowns
        Nt = size(gr%Verts, 1)
        allocate(u0(size(gr%verts, 1), size(rt%cpts, 1), 4))
        allocate(ures(size(gr%verts, 1), size(rt%cpts, 1), 4))

        ! We use a constant state of the boundary value as initial condition.
        call HomInit(gr, rt, u0, twdrpbf)
        
        ! Determine the stiffness of the grid and an initial speed estimate.
        gs = FindStiff(gr)
        fs = SpeedEstSplit(gr, rt, u0, twdrpbf, 0.0D0)
        write (*,*) "System speed estimate is", fs, " Grid stiffness is", gs
        write (*,*) "using as CFL ", cfl

        ! Solve the Euler equations
        call solDSDGSSPRK33(gr, rt, ures, u0, twdrpbf, tend, cfl)

        ! Export the results
        write (*,*) "writing ", size(ures, 1), " x ", size(ures, 2), " nodal values"
        call GridHDF5Export(gr, rt, "output/" // trim(ofname) // "gridexport.h5")
        call SolHDF5export(ures, tend, "output/" // trim(ofname) // "solution.h5")
end program twdrp
