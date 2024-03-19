! This program contains unit tests for the module ndgmats, in particular the creation
! of discrete orthogonal polynomial bases and Gaussian quadratures.

program tests
        use ndgmats
        
        integer, parameter :: n = 3
        integer, parameter :: nquad = 4
        integer :: i, j, k
        double precision, dimension(nquad) :: glpts, glws
        double precision, dimension(10, 10) :: M, Q
        double precision, dimension(10, 2) :: cpts
        double precision :: r
        
        glpts = GL(nquad)
        glws = GLW(nquad, glpts)
        write (*, *) glpts
        write (*, *) glws
        do k = 0, 10
                write (*, *) mod(k, 3), modulo(k, 3)
        end do
       
        cpts = pts(n)
        do k = 1, 10
                write (*, '(F4.3, xxx, F4.3)') cpts(k, :)
        end do

        write (*,*) " vandermondes"
        do i=1, 10
                do j=1,10
                        Q(i, j) = MonBF(j-1, cpts(i, :))
                end do
        end do
        
        do i=1, 10
                do j=1, 10
                        write (*, '(F4.2, xxx)', advance="no") Q(i, j)
                end do
                write (*,*) ""
        end do

        write (*, *) "Monomial Integrals up to order 10"
        do i=0, 65
                write (*, *) MonBFI(i)
        end do

        write (*,*) "The orthonormal description"
        call DOBM(M, 10, cpts)
        do i=1, 10
                do j=1, 10
                        write (*, '(F8.1, xxx)', advance="no") M(i, j)
                end do
                write (*,*) ""
        end do
        write (*,*) "Now the resulting vandermondes" 
        do i=1, 10
                do j=1,10
                        Q(i, j) = DOBF(M, 10, j, cpts(i, :))
                end do
        end do
        
        do i=1, 10
                do j=1, 10
                        write (*, '(F4.2, xxx)', advance="no") Q(i, j)
                end do
                write (*,*) ""
        end do

        write (*,*) "And their orthogonality"
        do i=1, 10
                do j=1, 10
                        write (*, '(F8.2, xxx)', advance="no") dot_product(Q(:, i), Q(:, j))
                end do
                write (*,*) ""
        end do
end program tests
