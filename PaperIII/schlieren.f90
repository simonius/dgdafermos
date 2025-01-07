! This module is used to calculates numerical schlieren images
! used as a f2py module
module schlieren
    contains
    function ptschl(Pset, u, comp) result (sar)
        double precision :: Pset(:, :, :)
        double precision :: u(:, :, :)
        integer :: comp
        double precision :: sar(size(Pset, 1)*size(Pset, 2))
        ! aux
        double precision :: s, u1, u2
        integer :: i, j, k, Ntr, Ncp
        Ntr = size(Pset, 1)
        Ncp = size(Pset, 2)
        write (*,*) "We got an array of ", Ntr, " triangles"
        write (*,*) "These consists of", Ncp, " points per triangle"
        write (*,*) "The corresponding value array has size", size(u, 1), size(u, 2)
        do i=1, Ntr
            do j = 1, Ncp
                s = 0.0d0
                do k=1, Ncp
                    if (j .NE. k) then
                        s = max(s, abs(u(i, j, comp) - u(i, k, comp)) / norm2(Pset(i, j, :) - Pset(i, k, :)))      
                    end if
                end do
                sar((i-1)*Ncp + j) = s
            end do
        end do
        write (*,*) "Calculated schlieren images"
    end function ptschl
end module schlieren
