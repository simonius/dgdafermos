! This module contains routines to construct positive quadratures on scattered points
! These allow a mimetic discretisation of the per cell entropy

module quads
        ! The number of iterations for the POCS algorithm
        use ndgmats
        implicit none
        integer, parameter :: Niter = 1.0E5
        double precision, parameter :: eps = 1.0D-10
        

        contains

        ! Solves Aw = b with positive w using the POCS 
        ! - Projection onto convex sets algorithm
        ! every row of the linear system is thought as describing one hyperplane
        ! i.e. a convex set, and every positivity constraint is a convex set
        ! (i.e. a halfspace). We find a point in the intersection, if it exists
        ! by iterating the least squares projection onto both.
        ! A, w, b: The linear system under consideration
        ! n, m: the size of A
        ! maxrow: the last row considered. Used to iteratively search for maximal order
        ! returns the residual.
        function PosLS(A, w, b, n, m, maxrow)
                integer, intent(in) :: n, m, maxrow
                double precision, dimension(n, m), intent(in) :: A
                double precision, dimension(n), intent(in) :: b
                double precision, dimension(m), intent(inout) :: w
                double precision :: PosLS

                integer :: i, j
                double precision, dimension(m) :: wold
                do i=1,Niter
                        wold = w
                        do j=1,maxrow
                                ! Project onto the hyperplane given by the row
                                w = w + A(j, :)*(b(j) - dot_product(A(j, :), w))/(dot_product(A(j, :), A(j, :)) +1.0D-10)
                                ! Project onto the positivity constraint
                                w = max(0.0D0, w)
                        end do
                end do
                PosLS = norm2(matmul(A, w) - b)
                write (*,"(A, D8.3, A, D8.3)") "Last iterion moved by ", norm2(w - wold), " Residual is ", PosLS
        end function PosLS

        ! Tries to find a maximal accurate quadrature on the given points
        ! M: should describe a DOP basis,
        ! npts: number of collocation points.
        ! pts: the array of points
        ! w: the resulting quadrature. 
        subroutine FindQuad(M, npts, pts, w)
                integer, intent(in) :: npts
                double precision, intent(in), dimension(npts, npts) :: M
                double precision, intent(out), dimension(npts) :: w
                double precision, intent(in), dimension(npts, 2) :: pts
                ! Aux
                integer :: i, j, nBFmax, maxmon
                double precision :: Pstep
                double precision, dimension(npts) :: wold
                double precision, dimension(npts, npts) :: A
                double precision, dimension(npts) :: b

                ! We initialize everything with zero
                A = 0.0D0
                b = 0.0D0
                w = 0.0D0
                ! And we count the monomials for which we are exact.
                maxmon = 0
                ! Lets assume more than 2 times the point number is not realistic
                nBFmax = npts
                do i=1,nBFmax
                        ! Fill in the exactness matrix
                        do j=1,npts
                                A(i, j) = DOBF(M, npts, i, pts(j, :))
                        end do
                        ! And the integral of the basis functions
                        b(i) = DOBI(M, npts, i)
                        ! We search for a positive quadrature up to basis functio i
                        Pstep = PosLS(A, w, b, nBFmax, npts, i)
                        ! If the residual jumps above our prescribed epsilon we used one monomial
                        ! to much.
                        if (Pstep > eps) then
                                w = wold
                                write (*,*) "quadrature is exact for up to ", i -1, " monomials"
                                exit
                        else            ! Otherwise we retry with one basis function more
                                wold = w
                                maxmon = maxmon + 1
                        end if
                end do
                write (*,*) "quadrature is exact for up to ", maxmon, " monomials"
        end subroutine FindQuad

        ! This routine Calculates a quadrature for a triangle using 
        ! Fubini and a collapsed coordinate transform.
        ! the resulting quadrature has far more points than needed, but is
        ! guaranteed to exist and be exact. This is used to integrate
        ! basis functions etc.
        subroutine GLQuadTri(p, pts, w)
                integer, intent(in) :: p
                double precision, intent(out), dimension(p*p, 2) :: pts
                double precision, intent(out), dimension(p*p) :: w
                ! aux variables
                double precision, dimension(p) :: glpts
                double precision, dimension(p) :: glws
                double precision :: x, y, wx, wy
                integer :: i, j, k
                
                ! Use a one-dimensional Gauß-Lobatto quadrature
                glpts = GL(p)
                glws = GLW(p, glpts)

                ! Use the one-dimensional quadrature for the double integral
                ! \Int_0^1 \Int_0^(1-x) f(x, y) dy dx
                do i=1, p
                        do j=1, p
                                k = j+(i-1)*p
                                x = 0.5D0*(glpts(i) + 1.0D0)
                                y = 0.5D0*(glpts(j) + 1.0D0)*(1-x)
                                pts(k, 1) = x
                                pts(k, 2) = y
                                wx = 0.5D0*glws(i)
                                wy = 0.5D0*glws(j)*(1-x)
                                w(k) = wx*wy
                        end do
                end do
        end subroutine GLQuadTri
end module quads
