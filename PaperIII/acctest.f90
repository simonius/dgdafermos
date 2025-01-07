! This file carries out a convergence analysis using a moving bump
! that wanders from the left of the domain to the right.
! The resulting program is called as

!       ./acctest grids/acc1.1.node grids/acc1.1.ele acc1 ord 2.0

! to carry out the test on grid acc1, save the output with the
! keyword acc1 and solve up to t=2.0

! A boundary function for free-flow
subroutine FFbf(uo, u, x, t)
        use euler
        real(np), intent(in) :: u(4)
        real(np), intent(in) :: x(2)
        real(np), intent(in) :: t
        real(np), intent(out) :: uo(4)
        uo = ConsVar(1.0_np, 1.0_np, 0.0_np, 1.0_np)
end subroutine FFbf


program acctest
        use sdgsolver           ! Include components needed
        use trianggrid
        use trianghp
        use testcases
        use opshp
        implicit none

        type(Grid) :: gr
        type(nprt) :: rt
        external :: FFbf
        real(np), allocatable :: u0(:, :, :)
        real(np), allocatable :: ures(:, :, :)
        real(np) :: tend = 0.1, cfl = 0.3
        real(np) :: dnodes
        integer :: ord
        character(len=256) :: nfname, efname, ofname, ordstring, tmaxstring

        ! Parse filenames
        if (command_argument_count() .NE. 5) then
                write (*,*) "wrong number of arguments!"
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
        
        ! Load the grid in the two filenames from memory
        call GridTriangleImp(gr, nfname, efname)
        
        ! Export the grid combined with the reference triangle
        call GridHDF5export(gr, rt, "output/" // trim(ofname) // "gridexport.h5")

        ! Allocate memory for the conserved variables
        allocate(u0(size(gr%verts, 1), size(rt%cpts, 1), 4))
        allocate(ures(size(gr%verts, 1), size(rt%cpts, 1), 4))

        ! Print diagnostic messages
        write(*,*) "conserved variavles have size", size(u0, 1), size(u0, 2), size(u0, 3)
        write(*,*) "Every boundary has ", rt%Nbdt, " nodes"

        ! Initialize the conserved variables with the initial condition
        call MBinit(gr, rt, u0)
        
        ! Carry out timesteps into ures
        call solDSDGSSPRK33(gr, rt, ures, u0, FFbf, tend, cfl)

        ! Overwrite initial condition in u0 with reference solution at tend.
        call MBref(gr, rt, u0, tend)
        dnodes = dble(size(gr%verts, 1)*size(rt%cpts, 1))
        
        ! Calculate and print the error norm
        write(*,*) "Error is", norm2((u0(:, :, 1) - ures(:, :, 1)) / dnodes)

        ! Export the difference for inspection
        call SolHDF5export(u0 - ures, tend, "output/" // trim(ofname) // "solution.h5")
end program acctest
