! This module contains functions that calculate the 
! matrices needed for our nodal dg solver.
! A second part can be found in the file "ops.f90"
module NDGMats
        implicit none

!       the 3 edges of the reference triangle
        double precision, dimension(2) :: v1 = [0.0D0, 0.0D0]
        double precision, dimension(2) :: v2 = [1.0D0, 0.0D0]
        double precision, dimension(2) :: v3 = [0.0D0, 1.0D0] 
        double precision :: piconst = 4*ATAN(1.0D0)
        contains
! returns the Legendere Polynomial
! checked
        pure recursive function Leg(n, x) result (L)
                integer, intent(in) :: n
                double precision, intent(in) :: x
                double precision :: L

                if (n == 0) then
                        L = 1.0
                else if (n == 1) then
                        L = x
                else
                        L = ((2*n-1)*x*Leg(n-1, x) - (n-1)*Leg(n-2, x))/n
                end if
        end function Leg

! recursively gives the derivative of
! Lagrange polynomial $n$, counting from zero,
! at position x
        pure recursive function dLeg(n, x) result (dL)
                integer, intent(in) :: n
                double precision, intent(in) :: x
                double precision :: dL
                if (n == 0) then
                        dL = 0.0
                else if (n == 1) then
                        dL = 1.0
                else
                        dL = ((2*n-1)*(Leg(n-1, x) + x*dLeg(n-1, x)) - (n-1)*dLeg(n-2, x))/n
                end if
        end function dLeg

! Finds a zero of the n-th Legendre polynomials using
! bisection
! checked
        pure function LegZero(n, xl, xr)
                integer, intent(in) :: n
                double precision, intent(in) :: xl, xr
                double precision :: LegZero
                integer :: i
                double precision :: fl, fm, fr, lb, ub
                lb = xl
                ub = xr
                fl = Leg(n, lb)
                fr = Leg(n, ub)

                do i=1, 64
                        fm = Leg(n, 0.5*(lb + ub))
                        if (fm*(fr - fl) > 0) then
                                ub = 0.5*(lb + ub)
                        else
                                lb = 0.5*(lb + ub)
                        end if
                end do
                LegZero = 0.5*(lb + ub)
        end function LegZero

! returns the n-th chebyhow points
        function CP(n) result (pts)
                integer, intent(in) :: n
                double precision, dimension(n) :: pts
                integer :: k
                do k=1, n
                        pts(k) = cos((2*k-1)/(2*n)*piconst)
                end do 
        end function CP

! returns the n-th Gauss-Legendere points
        recursive function GL(n) result (pts)
                integer, intent(in) :: n
                double precision, dimension(n) ::  pts

                integer :: i
                double precision, dimension(n-1) :: lopts
                ! For linear, it is clear. 
                if (n == 1) then
                        pts(1) = 0.0
                ! Else they lie outside the assumed domain.
                else
                        lopts = GL(n-1)
                        pts(1) = LegZero(n, -1.0D0, lopts(1))
                        pts(n) = LegZero(n, lopts(n-1), 1.0D0)
                        do i = 2, n-1
                                pts(i) = LegZero(n, lopts(i-1), lopts(i))
                        end do
                end if
        end function GL

! returns the n-th Gauss-Legendre weights
        function GLW(n, pts) result (w)
                integer, intent(in) :: n
                double precision, dimension(n), intent(in) :: pts
                double precision, dimension(n) :: w

                integer :: i
                do i=1, n
                        w(i) = 2.0D0/((1.0D0-pts(i)**2)*dLeg(n, pts(i))**2)
                end do
        end function GLW

! Enters points, at the moment only order 3
        function pts(q)
                integer, intent(in) :: q
                double precision, dimension ((q+2)*(q+1)/2, 2) :: pts
                double precision, dimension (q) :: glpts
                integer :: i
                double precision :: r

                glpts = GL(q)

                do i=1, q
                        r = (glpts(i) + 1)*0.5D0
                        pts(i, :) = v2*r + (1-r)*v1
                        pts((q)+i, :) = v3*r + (1-r)*v2
                        pts(2*(q)+i, :) = v1*r + (1-r)*v3
                end do
                if ((q+2)*(q+1)/2 > 3*q) then
                        pts((q+2)*(q+1)/2, :) = (v1 + v2 + v3)/3
                end if
        end function pts


        function factorial(n) result(fac)
                integer, intent(in) :: n
                integer :: fac
                integer :: i
                fac = 1
                do i=2, n
                        fac = fac*i
                end do
        end function factorial

        ! returns the n-th element of the pascal triangle of monomials
        ! n is zero indexed !
        function MonBF(n, a)
                integer, intent(in) :: n
                double precision, intent(in), dimension(2) :: a
                double precision :: MonBF
                integer :: i, j, num
                j = 0
                num = 0
                iloop: do i=0,n
                        do j=0,i
                                if (num == n) then
                                        exit iloop
                                else
                                        num = num + 1
                                end if
                        end do
                end do iloop
                MonBF = factorial(i-j)*a(1)**(i-j)*factorial(j)*a(2)**(j) 
        end function MonBF

        ! returns the derivative of the
        ! n-th element of the pascal triangle of monomials
        ! n is zero indexed !
        function MonBFx(n, a)
                integer, intent(in) :: n
                double precision, intent(in), dimension(2) :: a
                double precision :: MonBFx
                ! row (and column in the triangle
                integer :: i, j, num
                i=0
                j=0
                num = 0
                iloop: do i=0,n
                        do j=0,i
                                if (num == n) then
                                        exit iloop
                                else
                                        num = num + 1
                                end if
                        end do
                end do iloop
                if (i == j) then
                        MonBFx = 0.0D0
                else
                        MonBFx = factorial(i-j)*(i-j)*a(1)**(i-j-1)*factorial(j)*a(2)**(j)
                end if
        end function MonBFx

        ! returns the n-th element of the pascal triangle of monomials
        ! n is zero indexed !
        function MonBFy(n, a)
                integer, intent(in) :: n
                double precision, intent(in), dimension(2) :: a
                double precision :: MonBFy
                ! row (and column in the triangle
                integer :: i, j, num
                j = 0 
                num = 0
                iloop: do i=0,n
                        do j=0,i
                                if (num == n) then
                                        exit iloop
                                else
                                        num = num + 1
                                end if
                        end do
                end do iloop
                !write (*,*)  "selected monomial ", i-j, j
                if (0 == j) then
                        MonBFy = 0.0D0
                else    
                        MonBFy = factorial(i-j)*a(1)**(i-j)*factorial(j)*j*a(2)**(j-1)
                end if
        end function MonBFy


        ! Solves II x^k y^l d y d x  
        ! = _I x^k [(1-x)^(l+1)/(l+1)] d x
        recursive function PI(k, l) result (val)
                integer, intent(in) :: k, l
                double precision :: val 
                if (l > 0) then
                ! I x^k (1-x)^(l+1)/(l+1) d x = [- x^k (1-x)^l] + I x^(k+1)*(1-x)^l/(k+1)
                        val =  l/dble(k+1)*PI(k+1, l-1)
                else
                        ! I x^k (1-x) d x = [x^(k+1)/(k+1)]_0^1 - [x^(k+2)/(k+2)]
                        val = 1/dble(k+1) - 1/dble(k+2)
                end if
                
        end function PI

        ! Integrates a monomial over the reference triangle using a quadrature rule.
        function II(k, l) result (val)
                integer, intent(in) :: k, l
                double precision :: val
                integer, parameter :: nquad = 5
                double precision, dimension(nquad) :: qpts, qw
                integer :: i
                double precision :: f

                qpts = GL(nquad)
                qw = GLW(nquad, qpts)
                qpts = (qpts + 1.D00) * 0.5D0
                qw = qw * 0.5D0
                val = 0.0D0
                do i=1,nquad
                        f = qpts(i)**k*(1-qpts(i))**(l+1)/dble(l+1)
                        val = val + qw(i)*f
                end do
        end function II

        function MonBFI(n)
                integer, intent(in) :: n
                double precision :: MonBFI
                ! row (and column in the triangle
                integer :: i, j, num, k, l
                num = 0
                j = 0
                iloop: do i=0,n
                        do j=0,i
                                if (num == n) then
                                        exit iloop
                                else
                                        num = num + 1
                                end if
                        end do
                end do iloop
                k = i-j
                l = j 
                
                ! The integration of monomials over the reference triangle can be carried out in two ways
                ! either by exact integration (function PI) or by a quadrature rule of sufficient order (II)
                MonBFI = dble(factorial(i-j)*factorial(j))*PI(k, l)
                !MonBFI = dble(factorial(i-j)*factorial(j))*II(k, l) 
       end function MonBFI

        ! Calculates a discrete orthogonal basis from the monomials, saved in a matrix M
        ! Every row column in M contains the coefficients for the monomial basis.
        subroutine DOBM(M, n, pts)
                integer, intent(in) :: n
                double precision, dimension(n, n), intent(out) :: M
                double precision, dimension(n, 2), intent(in) :: pts
                
                integer :: i, j
                double precision, dimension(n) :: a
                double precision :: r
                ! space for the new orthogonal basis
                double precision, dimension(n, n) :: O

                M = 0.0D0
                ! loop over all base vectors
                do i=1, n
                        M(i, i) = 1.0D0
                        ! evaluate the function that should be evaluated
                        do j=1, n
                                a(j) = MonBF(i-1, pts(j, :))
                        end do
                        ! loop over all previous new base vectors
                        do j=1,i-1
                                M(:, i) = M(:, i) -dot_product(a, O(:, j))/(dot_product(O(:, j), O(:, j)) + 1.0D-10)*M(:, j)
                        end do
                        ! save the new base vectors
                        do j=1,n
                                O(j, i) = DOBF(M,n, i, pts(j, :))
                        end do
                        ! Normalize
                        r = norm2(O(:, i)) + 1.0D-10
                        O(:, i) = O(:, i) / r
                        M(:, i) = M(:, i) / r

                end do               
        end subroutine DOBM

        ! Evaluates the Discrete Orthogonal basis function k at x
        ! k is one-indexed
        function DOBF(M, n, k, x)
                integer, intent(in) :: n, k
                double precision, dimension(n, n), intent(in) :: M
                double precision, dimension(2), intent(in) :: x
                double precision :: DOBF
                integer :: i

                DOBF = 0.0D0
                do i=1,k
                        DOBF = DOBF + MonBF(i-1, x)*M(i, k)
                end do
        end function DOBF

        ! Evaluates the derivative in x direction of the
        ! Discrete Orthogonal basis function k at x
        ! k is one-indexed
        function DOBFx(M, n, k, x)
                integer, intent(in) :: n, k
                double precision, dimension(n, n), intent(in) :: M
                double precision, dimension(2), intent(in) :: x
                double precision :: DOBFx
                integer :: i

                DOBFx = 0.0D0
                do i=1,k
                        DOBFx = DOBFx + MonBFx(i-1, x)*M(i, k)
                end do
        end function DOBFx

        ! Evaluates the derivative in the y direction of the
        ! Discrete Orthogonal basis function k at x
        ! k is one-indexed
        function DOBFy(M, n, k, x)
                integer, intent(in) :: n, k
                double precision, dimension(n, n), intent(in) :: M
                double precision, dimension(2), intent(in) :: x
                double precision :: DOBFy
                integer :: i

                DOBFy = 0.0D0
                do i=1,k
                        DOBFy = DOBFy + MonBFy(i-1, x)*M(i, k)
                end do
        end function DOBFy


        ! Calculates the Vandermonde Matrix with respect to this discrete orthonormal basis
        subroutine DOBV(V, M, n, pts)
                integer, intent(in) :: n
                double precision, dimension(n, n), intent(out) :: V
                double precision, dimension(n, n), intent(in) :: M
                double precision, dimension(n, 2), intent(in) :: pts
                integer :: i, j
                do i=1,n
                        do j=1,n
                                V(j, i) = DOBF(M,n, i, pts(j, :))
                        end do
                end do
        end subroutine DOBV

        ! Calculates the der. Vandermonde Matrix in the x direction
        ! with respect to this discrete orthonormal basis
        subroutine DOBVx(Vx, M, n, pts)
                integer, intent(in) :: n
                double precision, dimension(n, n), intent(out) :: Vx
                double precision, dimension(n, n), intent(in) :: M
                double precision, dimension(n, 2), intent(in) :: pts
                integer :: i, j
                do i=1,n
                        do j=1,n
                                Vx(j, i) = DOBFx(M,n, i, pts(j, :))
                        end do
                end do
        end subroutine DOBVx

        ! Calculates the der. Vandermonde Matrix in the y direction
        ! with respect to this discrete orthonormal basis
        subroutine DOBVy(Vy, M, n, pts)
                integer, intent(in) :: n
                double precision, dimension(n, n), intent(out) :: Vy
                double precision, dimension(n, n), intent(in) :: M
                double precision, dimension(n, 2), intent(in) :: pts
                integer :: i, j
                do i=1,n
                        do j=1,n
                                Vy(j, i) = DOBFy(M,n, i, pts(j, :))
                        end do
                end do
        end subroutine DOBVy

        ! Calculates the integral of the DOBF k over the reference triangle
        function DOBI(M, n, k)
                integer, intent(in) :: n, k
                double precision, dimension(n, n) :: M
                double precision :: DOBI
                integer :: i
                DOBI = 0.0D0
                do i=1,k
                        DOBI = DOBI + MonBFI(i-1)*M(i, k)
                end do
        end function DOBI

        ! Monte Carlo approximation of DOBI, used for verification
        function DOBImc(M, n, k) result (res)
                integer, intent(in) :: n, k
                double precision, dimension(n, n) :: M
                double precision :: res
                double precision, dimension(2) :: r
                integer :: i, Niter = 10E6
                res = 0.0D0
                do i=1,Niter
                        call random_number(r)
                        if (r(1) + r(2) .LE. 1.0D0) then
                                res = res + DOBF(M, n, k, r)
                        end if
                end do
                res = res / Niter
        end function DOBImc
end module NDGMats
