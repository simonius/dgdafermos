! This file carries out a shock tube style test for the stabilized DG scheme.
! called as ./dsdgtest gridnodes gridtriangles savename tmax

! The corresponding boundary function for free flow
subroutine FFbf(uo, u, x, t)
        use euler
        double precision, intent(in) :: u(4)
        double precision, intent(in) :: x(2)
        double precision, intent(in) :: t
        double precision, intent(out) :: uo(4)
        uo = ConsVar(0.125D0, 0.0D0, 0.0D0, 0.1D0)
end subroutine FFbf


program dgtest
        use sdgsolver   ! import needed components
        use triangles
        use testcases
        use ops
        implicit none

        type(Grid) :: gr
        type(RefTriangle) :: rt
        external :: FFbf
        double precision, allocatable :: u0(:, :, :)
        double precision, allocatable :: ures(:, :, :)
        double precision :: tend = 0.4, cfl = 0.1
        character(len=256) :: nfname, efname, ofname, tmaxstring

        ! Parse filenames
        if (command_argument_count() .NE. 4) then
                write (*,*) "Missing input filenames!"
                stop
        end if
        call get_command_argument(1, nfname)
        call get_command_argument(2, efname)
        call get_command_argument(3, ofname)
        call get_command_argument(4, tmaxstring)

        read (tmaxstring, *) tend

        ! Initialize the reference triangle
        call DGinit(rt)
        ! Import the grid
        call GridTriangleImp(gr, nfname, efname)
        
        ! Export both
        call GridHDF5export(gr, rt, "output/" // trim(ofname) // "gridexport.h5")

        ! Allocate space for conserved variables
        allocate(u0(size(gr%verts, 1), size(rt%cpts, 1), 4))
        allocate(ures(size(gr%verts, 1), size(rt%cpts, 1), 4))

        ! Diagnostic output
        write(*,*) "conserved variavles have size", size(u0, 1), size(u0, 2), size(u0, 3)
        write(*,*) "Every boundary has ", rt%Nbdt, " nodes"
        
        ! Initialize the conserved variables. Lax shocktube is a 2d version 
        ! of Lax shocktube using translation invariance.
        ! The Kurganov Sedov test uses a rotational symmetry.

        !call LaxInit(gr, rt, u0)
        call KSinit(gr, rt, u0)
        
        ! Carry out timesteps
        call solDSDGSSPRK33(gr, rt, ures, u0, FFbf, tend, cfl)

        ! Export the results
        call SolHDF5export(ures, tend, "output/" // trim(ofname) // "solution.h5")
end program dgtest
