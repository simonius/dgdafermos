! A fortran programm to test the quadratures
program qtests
        use ndgmats
        use quads
        implicit none
        integer, parameter :: n = 10
        integer, parameter :: q = 3
        integer, parameter :: gls = 5
        integer :: i, j, k
        double precision, dimension(n, n) :: M
        double precision, dimension(n, 2) :: cpts
        double precision, dimension(n) :: w

        double precision, dimension(gls*gls, 2) :: glT
        double precision, dimension(gls*gls) :: glTw

        double precision :: qr1, qr2, qr3

        cpts = pts(q)
        call DOBM(M, n, cpts)

        call FindQuad(M, n, cpts, w)
        write (*, *) w

        write (*, *) "The GL quadrature on the triangle is given by"
        call GLQuadTri(gls, glT, glTw)
        do i=1, gls*gls
                write(*,"(F6.3, x, F6.3, x, F6.3)") glT(i, :), glTw(i)
        end do

        write (*,*) "now testing all quadratures"

        do i=1, n
                ! calculate first value
                qr1 = 0.D0
                do j=1,n
                        qr1 = qr1 + w(j)*DOBF(M, n, i, cpts(j, :))
                end do
                qr2 = 0.D0
                do j=1,gls*gls
                        qr2 = qr2 + glTw(j)*DOBF(M, n, i, glT(j, :))
                end do
                qr3 = DOBI(M, n, i)
                write (*,"(D15.6, x, D15.6, x, D15.6)") qr1, qr2, qr3
        end do

end program qtests
