! This module constructs the two-dimensional matrices for our nodal DG method
module ops
        use ndgmats             ! inclusion of needed modules
        use quads
        use triangles 
        use dpla                ! double precision linear algebra
        implicit none

        contains

        ! Calculates the Mass matrix and the stiffness matrices
        ! in x and y directions using numerical quadrature.
        ! M: Description of the DOP basis. 
        ! Mn: The Mass matrix in the nodal basis
        ! Sxn, Syn: The Stiffness matrices in x and y direction
        ! cpts: The used collocation points
        ! n: the number of collocation points.
        subroutine MassMat(M, Mn, Sxn, Syn, cpts, n)
                implicit none
                double precision, dimension(n, n) :: M
                double precision, dimension(n, n) :: Mn, Sxn, Syn
                double precision, dimension(n, 2) :: cpts
                integer, intent(in) :: n
                ! Aux variables
                double precision, dimension(n*n, 2) :: glT
                double precision, dimension(n*n) :: glTw
                double precision, dimension(n, n) :: V, Mm, Sxm, Sym
                double precision :: qr, sxqr, syqr
                integer :: i, j, k
                ! Find a suitble quadrature
                call GLQuadTri(n, glT, glTw)
                ! Get the Vandermone matrix
                call DOBV(V, M, n, cpts)
                
                ! Calculation of the Modal Mass matrix using the quadrature
                do i=1, n
                        do j=1, n
                                qr = 0.0D0
                                sxqr = 0.0D0
                                syqr = 0.0D0
                                do k=1,n*n
                                        qr = qr + glTw(k)*DOBF(M, n, i, glT(k, :))*DOBF(M, n, j, glT(k, :))
                                        sxqr = sxqr + glTw(k)*DOBFx(M, n, i, glT(k, :))*DOBF(M, n, j, glT(k, :))
                                        syqr = syqr + glTw(k)*DOBFy(M, n, i, glT(k, :))*DOBF(M, n, j, glT(k, :))
                                end do
                                Mm(i, j) = qr 
                                Sxm(i, j) = sxqr
                                Sym(i, j) = syqr
                        end do
                end do
                ! Info: V^-1 = V^T, b.c. discrete orthogonal polynomials (DOPs) are used
                ! We can therefore use V * M * V^T to convert from DOP to nodal basis. 
                Mn = matmul(V, matmul(Mm, transpose(V)))
                Sxn = matmul(V, matmul(Sxm, transpose(V)))
                Syn = matmul(V, matmul(Sym, transpose(V)))
        end subroutine MassMat


        ! Calculate the surface matrix elements
        ! in the modal basis and the Side-Splitting framework
        ! M: defines the DOP basis
        ! B: is a stack of 3 surface matrices.
        ! rt: defines the reference triangle
        subroutine SurfMatSplit(M, B, rt)
                implicit none
                type(RefTriangle) :: rt
                double precision, dimension(rt%ncol, rt%ncol) :: M
                double precision, dimension(3, rt%ncol, rt%ncol) :: B
                

                double precision, dimension(rt%ncol) :: gln, glws
                double precision :: lam
                double precision, dimension(2) :: x
                integer :: i, j, k
                ! Calculation of a suitable quadrature
                gln = GL(rt%ncol)
                glws = GLW(rt%ncol, gln)

                ! calculate the matrix elements
                B = 0.0D0
                do i=1, rt%ncol
                        do j=1, rt%ncol
                                ! Sum the quadrature
                                do k=1,rt%ncol
                                        lam = 0.5D0*(gln(k) + 1.0D0)
                                        x = (1-lam) * rt%x(1, :) + lam*rt%x(2, :)
                                        B(1, i, j) = B(1, i, j) &
                                                + norm2(rt%n(1, :))*0.5*glws(k)*DOBF(M, rt%ncol, i, x)*DOBF(M, rt%ncol, j, x)
                                        x = (1-lam) * rt%x(2, :) + lam*rt%x(3, :)
                                        B(2, i, j) = B(2, i, j) + &
                                                norm2(rt%n(2, :))*0.5*glws(k)*DOBF(M, rt%ncol, i, x)*DOBF(M, rt%ncol, j, x)
                                        x = (1-lam) * rt%x(3, :) + lam*rt%x(1, :)
                                        B(3, i, j) = B(3, i, j) + &
                                                norm2(rt%n(3, :))*0.5*glws(k)*DOBF(M, rt%ncol, i, x)*DOBF(M, rt%ncol, j, x)
                                end do
                        end do
                end do
                
        end subroutine SurfMatSplit


        ! Calculates the surface matrix elements
        ! In the modal basis and the MFSBP setting
        ! M: Describes the DOP basis
        ! Bx, By: The surface operators in x and y direction
        ! rt: the reference triangle
        ! n: number of nodes
        subroutine SurfMat(M, Bx, By, rt, n)
                implicit none
                double precision, dimension(n, n) :: M
                type(RefTriangle) :: rt
                double precision, intent(out), dimension(n, n) :: Bx, By
                integer, intent(in) :: n
                integer :: nsurf
                
                double precision, dimension(n) :: gln, glws
                double precision, dimension(3*n, 2) :: sqn
                double precision, dimension(3*n, 2) :: sqw
                double precision, dimension(2) :: sqr
                double precision :: lam
                integer i, j, k

                nsurf = 3*rt%Nbdt

                ! Calculate a surface quadrature of high enough order from a GL quadrature
                gln = GL(n)
                glws = GLW(n, gln)
                do i=1, n
                        lam = 0.5D0*(gln(i) + 1.0D0)
                        ! The nodes of the first, second and third side
                        do j=1, 3
                                sqn((j-1)*n + i, :) = RT%x(j, :) + lam*(RT%x(mod(j, 3)+1, :) - RT%x(j, :))
                        end do
                        ! Now the weights
                        do j=1, 3
                                sqw((j-1)*n + i, :) = RT%n(j, :)*glws(i)*0.5D0
                        end do
                end do

                ! Using this quadrature the matrix elements can be calculated.
                do i=1, n
                        do j=1, n
                                sqr = 0.0D0
                                do k=1,3*n
                                       sqr = sqr + sqw(k, :)*DOBF(M, n, i, sqn(k, :))*DOBF(M, n, j, sqn(k, :))
                                end do
                                Bx(i, j) = sqr(1)
                                By(i, j) = sqr(2)
                        end do
                end do

        end subroutine SurfMat

        ! Calculates the complete set of needed matrices/ functionals
        ! q: order      n: number of collocation points
        ! Mn: nodal mass matrix, Mni: its inverse
        ! Sn: the stiffness matrix
        ! w: a positive quadrature
        ! G: the dissipation matrix
        subroutine NodalMats(rt)
                type(RefTriangle) :: rt
                ! aux
                double precision, dimension(rt%ncol, rt%ncol) :: Mn, Mni, Sxn, Syn
                double precision, dimension(rt%ncol, rt%ncol) :: O, Bxm, Bym, V, G
                double precision, dimension(3, rt%ncol, rt%ncol) :: Bsp
                double precision, dimension(rt%ncol) :: w
                integer :: i, j

                ! Orthonormal basis 
                call DOBM(O, rt%ncol, rt%cpts)
                
                ! Positive quadrature on the collocation points
                call FindQuad(O, rt%ncol, rt%cpts, w)

                ! Mass and Stiffness matrix
                call MassMat(O, Mn, Sxn, Syn, rt%cpts, rt%ncol)
                
                ! Inversion of the Mass matrix
                Mni = invQR(Mn)

                allocate(rt%Br(rt%ncol, rt%ncol))
                allocate(rt%Bs(rt%ncol, rt%ncol))
                ! Calculation of the surface matrices
                call SurfMat(O, Bxm, Bym, rt, rt%ncol)
                call DOBV(V, O, rt%ncol, rt%cpts)
                rt%Br = matmul(V, matmul(Bxm, transpose(V))) !(:, 1:3*rt%Nbdt)
                rt%Bs = matmul(V, matmul(Bym, transpose(V))) !(:, 1:3*rt%Nbdt)


                ! Calcualtion of the split surface matrices
                call SurfMatSplit(O, Bsp, rt)
                allocate(rt%B(3, rt%ncol, rt%Nbdt))

                ! Translate into the nodal basis
                do i=1, 3
                        Bsp(i, :,:) =  matmul(V, matmul(Bsp(i, :, :), transpose(V)))
                end do
                ! Only save needed columns
                do i=1, 3
                        do j=1, rt%Nbdt
                                rt%B(i, :, j) = Bsp(i, :, rt%Bi(j, i))
                        end do
                end do

                ! We need a surface quadrature to integrate the entropy inequality predictors
                ! summing over the surface matrices suffices to get one. 
                allocate(rt%wsurf(rt%Nbdt))
                rt%wsurf = 0.0D0
                do i=1, rt%Nbdt
                       rt%wsurf(i) = sum(rt%B(1, :, i))
                end do
                write(*,*) "using the surface quadrature ", rt%wsurf


                allocate(rt%M(rt%ncol, rt%ncol))
                rt%M = Mn
                allocate(rt%Mi(rt%ncol, rt%ncol))
                rt%Mi = Mni
                allocate(rt%Sr(rt%ncol, rt%ncol))
                rt%Sr = Sxn
                allocate(rt%Ss(rt%ncol, rt%ncol))
                rt%Ss = Syn
                allocate(rt%w(rt%ncol))
                rt%w = w

                ! As a last element we calculate a suitable dissipation operator.
                call FindDispG(rt, G) 
                allocate(rt%G(rt%ncol, rt%ncol))
                rt%G = G
        end subroutine NodalMats


        ! A simple matrix exponential, uses the identity exp(z) = lim_k (I + A/k)^k
        function MatExp(Q) result (A)
                double precision, dimension(:, :), intent(in) :: Q
                double precision, dimension(size(Q, 1), size(Q, 2)) :: A
                double precision, dimension(size(Q, 1), size(Q, 2)) :: B
                integer :: i, maxit
                maxit = 40*int(norm2(Q))
                B = Q / maxit
                A = 0.0D0
                do i=1, size(Q, 1)
                        A(i, i) = 1.0D0
                end do
                do i=1, maxit
                        A = A + matmul(B, A) 
                end do

        end function MatExp

        ! Find a G using the approach given in part II
        subroutine FindDispG(rt, G)
                type(RefTriangle) :: rt
                double precision, dimension(size(rt%M, 1), size(rt%M, 2)), intent(out) :: G
                double precision, dimension(size(rt%M, 1), size(rt%M, 2)) :: Q
                double precision :: tlow = 0.0D0, tup = 1.0D0
                integer :: i, itmax = 32
                
                ! Calculate the heat equation discretisation on the element
                Q = matmul(rt%Sr, matmul(rt%Mi,transpose(rt%Sr))) &
                        + matmul(rt%Ss, matmul(rt%Mi, transpose(rt%Ss)))
                Q = matmul(rt%Mi, Q)

                ! Search for the minimal time until the matrix elements 
                ! in the resulting oeprator are all non-negative
                ! The bisection is carried out up to itmax times giving an
                ! answer exact up to 2^(-32).
                do i=1, itmax
                        G = MatExp(-0.5*(tlow + tup)*Q)

                        if (all(G .GE. -1.0D-8)) then
                                tup = 0.5*(tlow + tup)
                        else
                                tlow = 0.5*(tlow + tup)
                        end if
                end do
                write (*,*) "Heat equation positivity was reached in intervall", tlow, tup

                ! In the case of, for example, linear elements, the matrix is always positive as
                ! Q is always a positive generator itself. We use Q directly in this case.
                if ((tlow + tup) < 1.0E-6) then
                        write (*,*) "G Matrix is always positive. Writing out discretisation itself instead"
                        tup = 0.0
                        do i=1, size(Q, 1)
                                tup = max(tup, abs(Q(i, i)))
                        end do
                        G = -Q/max(norm2(Q), tup)
                        do i=1, size(Q, 1) 
                                G(i, i) = G(i, i) + 1.0
                        end do
                        G = MatExp(-Q)
                else
                        G = MatExp(-0.5*(tlow + tup)*Q)
                end if
                write (*,*) "G Matrix is"
                do i=1, size(G, 1)
                        write (*,"(10f8.3)") G(i, :)
                end do
                write (*,*) " "
                ! To convert between a generator and the time-step operator
                ! the identity must be subtracted
                do i=1, size(G, 1)
                        G(i, i) = G(i, i) - 1.0D0
                end do
        end subroutine

        ! This subroutine is used to calcualte the projection to lower order 
        ! polynomials on the reference triangle rt.
        subroutine FindLO(rt)
                type(RefTriangle) :: rt
                ! aux
                double precision, dimension(size(rt%M, 1),size(rt%M, 2)) :: ONB, DONB, LO, LOspec
                double precision, dimension(size(rt%M, 1)) :: a
                integer :: i, j, n
                n = size(rt%M, 1)
                
                ! First we find an DONB
                call DOBM(DONB, n, rt%cpts)

                ! Then we find an ONB, a matrix with a M orthogonal basis in its rows
                ! Here, M shall be the Gramian of the _continuous_ inner product on the rt.
                ONB = 0.0D0
                do i=1, size(rt%M, 1)
                        ! Calc nodal values
                        do j=1, size(rt%M, 1)
                                a(j) = DOBF(DONB, n, i, rt%cpts(j, :))
                        end do
                        ! Project out previous base vectors
                        do j=1, i-1
                                a = a - dot_product(ONB(:, j), matmul(rt%M, a))*ONB(:, j)      
                        end do
                        ! Normalize
                        a = a / sqrt(dot_product(a, matmul(rt%M, a)))
                        ! Save
                        ONB(:, j) = a
                end do

                ! Using this orthonormal basis we can write the projection
                ! in the ONB, i.e. the first basis vectors are kept,
                ! the higher modes are set to zero.
                LOspec = 0.0D0
                do i=1, min(6, n)
                        LOspec(i, i) = 1.0D0
                end do

                !  We check the result
                write (*,*) "Checking if ONB construction was succesfull"
                LO = matmul(transpose(ONB), matmul(rt%M, ONB))
                do i=1, size(rt%M, 1)
                        write (*,"(10f8.3)") LO(i, :)
                end do
                write (*,*) "Now resulting matrix"
                LO = matmul(ONB, matmul(LOspec, matmul(transpose(ONB), rt%M)))
                do i=1, size(rt%M, 1)
                        write (*,"(10f8.3)") LO(i, :)
                end do

                do i=1, size(rt%M, 1)
                        do j=1, size(rt%M, 1)
                                a(j) = DOBF(DONB, n, i, rt%cpts(j, :))
                        end do
                        write (*,*) "action onto basis function ", i, norm2(a - matmul(LO, a)), norm2(matmul(LO, a))
                end do

                ! And save it in the reference triangle
                allocate(rt%LO(n, n))
                rt%LO = LO
        end subroutine FindLO

        ! This subroutine calls all needed subroutines to initialize a reference triangle for our DG method.
        subroutine DGinit(rt)
                type(RefTriangle), intent(inout) :: rt

                call RefInitDG(rt)
                call NodalMats(rt)
                call FindLO(rt)
        end subroutine DGinit
end module ops
