module csv
        use param
        implicit none
        contains

        ! writes the values (x_k, u_k) as a comma separated list
        ! under the filename name 
        subroutine writecsv(x, u, name)
                real(np) :: x(:), u(:)
                character(Len=*) :: name
                ! auxillary
                integer :: io, i
                open(newunit=io, file=name)
                write(io, *) "x ; u"
                do i=1, size(x, 1)
                        write(io, *) x(i), " ; ", u(i)
                end do
                close(io)
        end subroutine writecsv
end module csv

