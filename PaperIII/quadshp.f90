! A high precision version of the quads module

module quadshp
        use param
        use bf
        use qpla
        implicit none
        contains

        ! finds the quadrature weights w.
        ! npts: number of collocation points.
        ! pts: the array of points
        ! w: the resulting quadrature.
        subroutine FindQuad(pts, w)
                real(hp), intent(out) :: w(:)
                real(hp), intent(in) :: pts(:, :)
                ! Aux
                integer :: i, j
                real(hp) :: A(size(pts, 1), size(pts, 1))
                real(hp) :: b(size(pts, 1))
                A = 0.0_hp
                b = 0.0_hp
                w = 0.0_hp
                do i=1, size(pts, 1)
                        do j=1, size(pts, 1)
                                A(i, j) = MonBF(i-1, pts(j, :))
                        end do
                        b(i) = MonBFI(i-1)
                end do
                w = matmul(invQR(A), b)
                do i=1, size(pts, 1)
                        if (w(i) < -epsilon(1000.0_hp)) then
                                write(*,*) "Quadrature node ", i, " is negative ", w(i)
                        end if
                end do
        end subroutine FindQuad

        ! This routine Calculates a quadrature for a triangle using
        ! Fubini and a collapsed coordinate transform.
        ! the resulting quadrature has far more points than needed, but is
        ! guaranteed to exist and be exact. This is used to integrate
        ! basis functions etc.
        subroutine GLQuadTri(p, pts, w)
                integer, intent(in) :: p
                real(hp), intent(out), dimension(p*p, 2) :: pts
                real(hp), intent(out), dimension(p*p) :: w
                ! aux variables
                real(hp), dimension(p) :: glpts
                real(hp), dimension(p) :: glws
                real(hp) :: x, y, wx, wy
                integer :: i, j, k

                ! Use a one-dimensional Gauß-Lobatto quadrature

                glpts = GLE(p)
                glws = GLW(p, glpts)

                if (debug) then
                        write(*,*) "Using one-dimensionsal GLE points", glpts
                        write(*,*) "Using one-dimensional GLE weights", glws
                end if
                ! Use the one-dimensional quadrature for the double integral
                ! \Int_0^1 \Int_0^(1-x) f(x, y) dy dx
                do i=1, p
                        do j=1, p
                                k = j+(i-1)*p
                                x = 0.5_hp*(glpts(i) + 1.0_hp)
                                y = 0.5_hp*(glpts(j) + 1.0_hp)*(1-x)
                                pts(k, 1) = x
                                pts(k, 2) = y
                                wx = 0.5_hp*glws(i)
                                wy = 0.5_hp*glws(j)*(1-x)
                                w(k) = wx*wy
                        end do
                end do
        end subroutine GLQuadTri
end module quadshp
