! This program is a proof of concept implementation of a unstructured DG
! solver for the Euler equations implementing the entropy rate criterion.
! It loads a triangular grid in the form of two filenames corresponding to 
! a grid in the file format of Skewshuks Triangle program and solves the 
! Euler equations on it up to the given time tend. Afterwards it
! saves the grid and solution and exits.
!       ./udsd nodefile elefile outputname tend 


! The boundary function for the free flow around the domain.
! uo: resulting outer value of the conserved variables
! u: inner value of the conserved variables
! x: position
! t: time.
subroutine FFbf(uo, u, x, t)
        use euler
        double precision, intent(in) :: u(4)
        double precision, intent(in) :: x(2)
        double precision, intent(in) :: t
        double precision, intent(out) :: uo(4)
        double precision, parameter :: M = 2.0, alp = 1.25
        double precision, parameter :: pie = 4*ATAN(1.0D0)

        ! Value for stream with angle alp and mach number M        
        uo = ConsVar(1.00D0, M*sqrt(1.4D0), M*sqrt(1.4D0)*sin(alp/360*2*pie), 1.0D0)
        
        ! Value for the forward facing step problem.
        !uo = ConsVar(1.4D0, 3.0D0, 0.0D0, 1.0D0) 
end subroutine FFbf

program udsdg
        use sdgsolver
        use triangles
        use testcases
        use ops
        implicit none

        type(Grid) :: gr
        type(RefTriangle) :: rt
        double precision, allocatable :: u0(:, :, :)
        double precision, allocatable :: ures(:, :, :)
        double precision :: tend, gs, cfl, fs
        integer :: Nt
        external :: FFbf
        character(len=256) :: nfname, efname, ofname, tmaxstring

        if (command_argument_count() .NE. 4) then
                write (*,*) "Missing input filenames!"
                stop
        end if

        ! Parsing of the command line arguments
        call get_command_argument(1, nfname)
        call get_command_argument(2, efname)
        call get_command_argument(3, ofname)
        call get_command_argument(4, tmaxstring)

        read (tmaxstring, *) tend
        
        ! Initialization of the reference triangle
        call DGinit(rt)
        
        ! Load the grid from the HDD 
        call GridTriangleImp(gr, nfname, efname)

        ! Allocation of arrays for the unknowns
        Nt = size(gr%Verts, 1)
        allocate(u0(size(gr%verts, 1), size(rt%cpts, 1), 4))
        allocate(ures(size(gr%verts, 1), size(rt%cpts, 1), 4))

        ! We use a constant state of the boundary value as initial condition.
        call HomInit(gr, rt, u0, FFbf)
        
        ! Determine the stiffness of the grid and an initial speed estimate.
        gs = FindStiff(gr)
        fs = SpeedEstSplit(gr, rt, u0, FFbf, 0.0D0)
        cfl = 0.1
        write (*,*) "System speed estimate is", fs, " Grid stiffness is", gs
        write (*,*) "using as CFL ", cfl

        ! Solve the Euler equations
        call solDSDGSSPRK33(gr, rt, ures, u0, FFbf, tend, cfl)

        ! Export the results
        write (*,*) "writing ", size(ures, 1), " x ", size(ures, 2), " nodal values"
        call GridHDF5Export(gr, rt, "output/" // trim(ofname) // "gridexport.h5")
        call SolHDF5export(ures, tend, "output/" // trim(ofname) // "solution.h5")
end program udsdg
