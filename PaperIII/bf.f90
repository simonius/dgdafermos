! Basis functions and 1D quadrature nodes for the rewritten matrix calculation
! using high precision floating point.

module bf
        use param
        implicit none
        
        contains

        ! evaluates Legendre polynomial n at x
        pure function Legf(n, x) result (f)
                integer, intent(in) :: n
                real(hp), intent(in) :: x
                real(hp) :: f
                ! auxillary
                real(hp), dimension(n+1) :: p
                integer :: k

                if (n .EQ. 0) then
                        f = 1.0_hp                 !*sqrt((2.0D0*n+1) / 2)
                else if (n .EQ. 1) then
                        f = x                   !*sqrt((2.0D0*n+1) / 2)
                else
                        p(1) = 1
                        p(2) = x
                        do k=2, n
                                p(k+1) = ((2*k-1)*x*p(k) - (k-1)*p(k-1))/k
                        end do
                        f = p(n+1)              !*sqrt((2.0D0*n+1) / 2)
                end if
        end function Legf


        ! recursively gives the derivative of
        ! Lagrange polynomial $n$, counting from zero,
        ! at position x
        pure recursive function dLeg(n, x) result (dL)
                integer, intent(in) :: n
                real(hp), intent(in) :: x
                real(hp) :: dL
                if (n == 0) then
                        dL = 0.0_hp
                else if (n == 1) then
                        dL = 1.0_hp
                else
                        dL = ((2*n-1)*(Legf(n-1, x) + x*dLeg(n-1, x)) - (n-1)*dLeg(n-2, x))/real(n, hp)
                end if
        end function dLeg

        ! Finds a zero of the n-th Legendre polynomial
        ! in the interval [xl, xr] using bisection
        ! checked
        pure function LegZero(n, xl, xr)
                integer, intent(in) :: n
                real(hp), intent(in) :: xl, xr
                real(hp) :: LegZero
                integer :: i
                real(hp) :: fl, fm, fr, lb, ub
                lb = xl
                ub = xr
                fl = Legf(n, lb)
                fr = Legf(n, ub)

                do i=1, 128
                        fm = Legf(n, 0.5_hp*(lb + ub))
                        if (fm*(fr - fl) > 0.0_hp) then
                                ub = 0.5_hp*(lb + ub)
                        else
                                lb = 0.5_hp*(lb + ub)
                        end if
                end do
                LegZero = 0.5_hp*(lb + ub)
        end function LegZero

        ! Finds a zero of the first derivaative of the n-th Legendre polynomial
        ! in the interval [xl, xr] using bisection
        ! checked
        pure function dLegZero(n, xl, xr)
                integer, intent(in) :: n
                real(hp), intent(in) :: xl, xr
                real(hp) :: dLegZero
                integer :: i
                real(hp) :: fl, fm, fr, lb, ub
                lb = xl
                ub = xr
                fl = dLeg(n, lb)
                fr = dLeg(n, ub)

                do i=1, 128
                        fm = dLeg(n, 0.5_hp*(lb + ub))
                        if (fm*(fr - fl) > 0.0_hp) then
                                ub = 0.5_hp*(lb + ub)
                        else
                                lb = 0.5_hp*(lb + ub)
                        end if
                end do
                dLegZero = 0.5_hp*(lb + ub)
        end function dLegZero

        ! returns the n-th chebyhow points
        function CP(n) result (pts)
                integer, intent(in) :: n
                real(hp), dimension(n) :: pts
                integer :: k
                do k=1, n
                        pts(k) = - cos(real(2*k-1, hp)/real(2*n, hp)*pihp)
                end do
        end function CP

        ! returns the chebyshow points of the second kind
        function CP2(n) result (pts)
                integer, intent(in) :: n
                real(hp), dimension(n) :: pts
                integer :: k
                do k=1, n
                        pts(k) = -cos(pihp*real(k-1, hp)/(real(n, hp) - 1.0_hp))
                end do
        end function CP2

        ! returns the n-th Gauss-Legendere points
        recursive function GLE(n) result (pts)
                integer, intent(in) :: n
                real(hp), dimension(n) ::  pts

                integer :: i
                real(hp), dimension(n-1) :: lopts
                ! For linear, it is clear. 
                if (n == 1) then
                        pts(1) = 0.0_hp
                ! Else they lie outside the assumed domain.
                else
                        lopts = GLE(n-1)
                        pts(1) = LegZero(n, -1.0_hp, lopts(1))
                        pts(n) = LegZero(n, lopts(n-1), 1.0_hp)
                        do i = 2, n-1
                                pts(i) = LegZero(n, lopts(i-1), lopts(i))
                        end do
                end if
        end function GLE

        ! returns the n-th Gauss-Lobatto points
        recursive function GLO(n) result (pts)
                integer, intent(in) :: n
                real(hp), dimension(n) ::  pts

                integer :: i
                real(hp), dimension(n-1) :: lopts
                ! For linear, it is clear. 
                if (n == 2) then
                        pts(1) = -1.0_hp
                        pts(2) = 1.0_hp
                ! Else they lie outside the assumed domain.
                else
                        lopts = GLO(n-1)
                        pts(1) = -1.0_hp
                        pts(n) = 1.0_hp
                        do i = 2, n-1
                                pts(i) = dLegZero(n-1, lopts(i-1), lopts(i))
                        end do
                end if
        end function GLO

        ! returns the n-th Gauss-Legendre weights
        function GLW(n, pts) result (w)
                integer, intent(in) :: n
                real(hp), dimension(n), intent(in) :: pts
                real(hp), dimension(n) :: w

                integer :: i
                do i=1, n
                        w(i) = 2.0_hp/((1.0_hp-pts(i)**2)*dLeg(n, pts(i))**2)
                end do
        end function GLW

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
                real(hp), intent(in), dimension(2) :: a
                real(Hp) :: MonBF
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
                MonBF = a(1)**(i-j)*a(2)**(j)
        end function MonBF

        ! returns the derivative of the
        ! n-th element of the pascal triangle of monomials
        ! n is zero indexed !
        function MonBFx(n, a)
                integer, intent(in) :: n
                real(hp), intent(in), dimension(2) :: a
                real(hp) :: MonBFx
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
                        MonBFx = 0.0_hp
                else
                        MonBFx = (i-j)*a(1)**(i-j-1)*a(2)**(j)
                end if
        end function MonBFx

        ! returns the n-th element of the pascal triangle of monomials
        ! n is zero indexed !
        function MonBFy(n, a)
                integer, intent(in) :: n
                real(hp), intent(in), dimension(2) :: a
                real(hp) :: MonBFy
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
                        MonBFy = 0.0_hp
                else
                        MonBFy = a(1)**(i-j)*j*a(2)**(j-1)
                end if
        end function MonBFy

        ! Solves II x^k y^l d y d x  
        ! = _I x^k [(1-x)^(l+1)/(l+1)] d x
        recursive function PI(k, l) result (val)
                integer, intent(in) :: k, l
                real(hp) :: val
                if (l > 0) then
                ! I x^k (1-x)^(l+1)/(l+1) d x = [- x^k (1-x)^l] + I x^(k+1)*(1-x)^l/(k+1)
                        val =  l/real(k+1, hp)*PI(k+1, l-1)
                else
                        ! I x^k (1-x) d x = [x^(k+1)/(k+1)]_0^1 - [x^(k+2)/(k+2)]
                        val = 1/real(k+1, hp) - 1/real(k+2, hp)
                end if

        end function PI

        function MonBFI(n)
                integer, intent(in) :: n
                real(hp) :: MonBFI
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

                ! The integration of monomials over the reference triangle is carried out
                ! by exact integration (function PI)
                MonBFI = PI(k, l)
       end function MonBFI

end module bf
