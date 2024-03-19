! This file contains the unit tests for the double precision linear algebra
! implemented in dpla.f90 - a QR decomposition is carried out for a random
! and backsubstitution and inversion are tested.

program latest
        use dpla
        implicit none
        double precision, dimension(4, 4) ::  A, Q, R, D, Ai
        double precision, dimension(4) :: x, b
        real, dimension(4, 4) :: C
        integer :: i, n
        n = 4
        call random_number(C)
        A = dble(C)
        
        call QRH(Q, R, A)
        
        write (*,*) "A" 
        do i=1, n
                write (*,*) A(i, :)
        end do
        write(*,*) "Q"
        do i=1, n
                write (*,*) Q(i, :)
        end do
        write(*,*) "R"
        do i=1, n
                write (*,*) R(i, :)
        end do
        
        D = matmul(Q, R) - A
        write (*,*) "Reassembly error"
        do i=1, n
                write (*,*) D(i, :)
        end do

        x = BS(R, b)
        write (*,*) "Backsubstitution error"
        write (*,*) matmul(R, x)
        
        Ai = invQR(A)
        D = matmul(Ai, A)
        write (*,*) "inverse error"
        do i=1, n
                write (*,*) D(i, :)
        end do

               
end program latest
