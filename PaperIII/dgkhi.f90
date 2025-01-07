! This program is a proof of concept implementation of a unstructured DG
! solver for the Euler equations implementing the entropy rate criterion.
! It loads a triangular grid in the form of two filenames corresponding to 
! a grid in the file format of Skewshuks Triangle program and solves the 
! Euler equations on it up to the given time tend. Afterwards it
! saves the grid and solution and exits.
!       ./dgkhi nodefile elefile outputname order tend 

! we advise tend = 3.5 on a 1.0 x 3.5 domain.

! This version is specialized to test a Kelvin-Helmholtz instability. 


! The boundary function for the free flow around the domain.
! uo: resulting outer value of the conserved variables
! u: inner value of the conserved variables
! x: position
! t: time.
subroutine KHIbf(uo, u, x, t)
        use euler
        real(np), intent(in) :: u(4)
        real(np), intent(in) :: x(2)
        real(np), intent(in) :: t
        real(np), intent(out) :: uo(4)
        real(np)  :: vy
        if (x(1) < 10.0_np) then
                vy = 0.0_np
                if (0.33_np .LE. x(2) .AND. x(2) .LE. 0.66_np) then
                        uo = ConsVar(2.0_np, 0.5_np, vy, 1.0_np)
                else
                        uo = ConsVar(1.0_np, 0.0_np, vy, 1.0_np)
                end if                        
        else
                uo = ConsVar(1.0_np, 0.0_np, 0.0_np, 1.0_np) 
        end if
end subroutine KHIbf

program dgkhi
        use sdgsolver
        use trianggrid
        use trianghp
        use testcases
        use opshp
        implicit none

        type(Grid) :: gr
        type(nprt) :: rt
        real(np), allocatable :: u0(:, :, :)
        real(np), allocatable :: ures(:, :, :)
        real(np) :: tend, gs, cfl, fs
        integer :: Nt, ord
        external :: KHIbf
        character(len=256) :: nfname, efname, ofname, ordstring, tmaxstring

        cfl = 0.15_np
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
        
        ! Initialization of the reference triangle
        call DGinit(rt, ord)
        
        ! Load the grid from the HDD 
        call GridTriangleImp(gr, nfname, efname)

        ! Allocation of arrays for the unknowns
        Nt = size(gr%Verts, 1)
        allocate(u0(size(gr%verts, 1), size(rt%cpts, 1), 4))
        allocate(ures(size(gr%verts, 1), size(rt%cpts, 1), 4))

        ! We use a constant state of the boundary value as initial condition.
        call HomInit(gr, rt, u0, KHIbf)
        
        ! Determine the stiffness of the grid and an initial speed estimate.
        gs = FindStiff(gr)
        fs = SpeedEstSplit(gr, rt, u0, KHIbf, 0.0D0)
        cfl = cfl / (dble(ord)**2 + dble(ord))
        write (*,*) "System speed estimate is", fs, " Grid stiffness is", gs
        write (*,*) "using as CFL ", cfl

        ! Solve the Euler equations
        call solDSDGSSPRK33(gr, rt, ures, u0, KHIbf, tend, cfl)

        ! Export the results
        write (*,*) "writing ", size(ures, 1), " x ", size(ures, 2), " nodal values"
        call GridHDF5Export(gr, rt, "output/" // trim(ofname) // "gridexport.h5")
        call SolHDF5export(ures, tend, "output/" // trim(ofname) // "solution.h5")
end program  dgkhi
