! This file carries out a shock tube style test for the stabilized DG scheme.
! called as 
!       >> ./dsdgtest gridnodes gridtriangles savename ord tmax

! The corresponding boundary function for free flow
subroutine FFbf(uo, u, x, t)
        use euler
        real(np), intent(in) :: u(4)
        real(np), intent(in) :: x(2)
        real(np), intent(in) :: t
        real(np), intent(out) :: uo(4)
        uo = ConsVar(0.125_np, 0.0_np, 0.0_np, 0.1_np)
end subroutine FFbf


program dgtest
        use sdgsolver   ! import needed components
        use trianggrid
        use trianghp
        use testcases
        use opshp
        implicit none

        type(Grid) :: gr
        type(nprt) :: rt
        integer :: ord
        external :: FFbf
        real(np), allocatable :: u0(:, :, :)
        real(np), allocatable :: ures(:, :, :)
        real(np) :: tend = 0.4, cfl = 0.3
        character(len=256) :: nfname, efname, ofname, tmaxstring, ordstring

        ! Parse filenames
        if (command_argument_count() .NE. 5) then
                write (*,*) "Missing input arguments!"
                stop
        end if
        call get_command_argument(1, nfname)
        call get_command_argument(2, efname)
        call get_command_argument(3, ofname)
        call get_command_argument(4, ordstring) 
        call get_command_argument(5, tmaxstring)

        read (ordstring, *) ord
        read (tmaxstring, *) tend
        cfl = cfl / (dble(ord)**2 + dble(ord))


        ! Initialize the reference triangle
        call DGinit(rt, ord)
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
