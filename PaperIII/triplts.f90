! This program exports the coordinates of the reference triangle
! to plot them in the publication

program triplts
        use param
        use csv
        use trianghp
        use opshp

        implicit none
        integer :: p
        character(len=256) :: Norder
        character(len=256) :: Savename
        type(nprt) :: rt
        
        if (command_argument_count() .NE. 2) then
                write (*,*) "Wrong arguments. Call >> ./triplts [p] [savename]"
                stop
        end if
        ! Parsing of the command line arguments
        call get_command_argument(1, Norder)
        call get_command_argument(2, Savename)
        
        read (Norder, *) p
        
        ! Initalization of the reference triangle
        call DGinit(rt, p)
        call writecsv(rt%cpts(:, 1), rt%cpts(:, 2), savename)
        
end program triplts
