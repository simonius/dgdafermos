! A high-precision version of the quadrature module.
! This module constructs the two-dimensional matrices for our nodal DG method
module opshp
        use quadshp
        use qpla
        use trianghp

        implicit none
        

        contains

        ! Calculates the Volume matrices, also allocates them in rt
        subroutine VolMat(rt)
                implicit none
                type(hprt) :: rt
                ! Aux variables
                real(hp), dimension(4*rt%q*rt%q, 2) :: glT  ! The auxillary volume quadrature
                real(hp), dimension(4*rt%q*rt%q) :: glTw
                real(hp), dimension(rt%ncol, rt%ncol) :: Mm, Sxm, Sym, Qm
                real(hp) :: qr, sxqr, syqr, diffqr
                integer :: i, j, k
                
                
                ! Find a suitble quadrature
                write (*,*) "Calculating High Order quadrature"
                call GLQuadTri(2*rt%q, glT, glTw)

                write (*,*) "Calculating the Vandermonde Matrix"
                ! Get the Vandermone matrix
                allocate(rt%V(rt%Ncol, rt%Ncol))
                do i=1, rt%ncol
                        do j=1, rt%ncol
                                rt%V(i, j) = MonBF(j-1, rt%cpts(i, :))
                        end do
                end do

                allocate(rt%Vi(rt%ncol, rt%ncol))
                rt%Vi = invQR(rt%V)

                write (*,*) "Calculating Mass, Derivative and Stiffness matrices"
                ! Calculation of the Modal Mass matrix using the quadrature
                do i=1, rt%ncol
                        do j=1, rt%ncol
                                qr = 0.0_hp
                                sxqr = 0.0_hp
                                syqr = 0.0_hp
                                diffqr = 0.0_hp
                                do k=1,4*rt%q*rt%q
                                        qr = qr + glTw(k)*MonBF(i-1, glt(k, :))*MonBF(j-1, glT(k, :))
                                        sxqr = sxqr + glTw(k)*MonBFx(i-1, glT(k, :))*MonBF(j-1, glT(k, :))
                                        syqr = syqr + glTw(k)*MonBFy(i-1, glT(k, :))*MonBF(j-1, glT(k, :))
                                        diffqr = diffqr + gltw(k)*(MonBFx(i-1, glT(k, :))*MonBFx(j-1, glT(k, :)) +&
                                                                MonBFy(i-1, glT(k, :))*MonBFy(j-1, glT(k, :)))
                                end do
                                Mm(i, j) = qr
                                Sxm(i, j) = sxqr
                                Sym(i, j) = syqr
                                Qm(i, j) = diffqr
                        end do
                end do
                
                allocate(rt%M(rt%ncol, rt%ncol))
                allocate(rt%Mi(rt%ncol, rt%ncol))
                allocate(rt%Sr(rt%ncol, rt%ncol))
                allocate(rt%Ss(rt%ncol, rt%ncol))
                allocate(rt%S2(rt%ncol, rt%ncol))
                rt%M = matmul(transpose(rt%Vi), matmul(Mm, rt%Vi))
                rt%Mi = invqr(rt%M)
                rt%Sr = matmul(transpose(rt%Vi), matmul(Sxm, rt%Vi))
                rt%Ss = matmul(transpose(rt%Vi), matmul(Sym, rt%Vi))
                rt%S2 = matmul(transpose(rt%Vi), matmul(Qm, rt%Vi))
        end subroutine VolMat

        ! Calculate the surface matrix elements
        ! B: is a stack of 3 surface matrices.
        ! rt: defines the reference triangle
        subroutine SurfMatSplit(rt)
                implicit none
                type(hprt) :: rt
                ! aux
                real(hp), dimension(3, rt%ncol, rt%ncol) :: Bmod, Bnod
                real(hp), dimension(2*rt%nbdt) :: gln, glws
                real(hp) :: lam
                real(hp), dimension(2) :: x
                integer :: i, j, k
                ! Calculation of a suitable quadrature
                gln = GLE(2*rt%nbdt)
                glws = GLW(2*rt%nbdt, gln)

                ! calculate the matrix elements
                Bmod = 0.0_hp
                do i=1, rt%ncol
                        do j=1, rt%ncol
                                ! Sum the quadrature
                                do k=1,2*rt%nbdt
                                        lam = 0.5_hp*(gln(k) + 1.0_hp)
                                        x = (1-lam) * rt%x(1, :) + lam*rt%x(2, :)
                                        Bmod(1, i, j) = Bmod(1, i, j) &
                                                + norm2(rt%n(1, :))*0.5_hp*glws(k)*MonBF(i-1, x)*MonBF(j-1, x)
                                        x = (1-lam) * rt%x(2, :) + lam*rt%x(3, :)
                                        Bmod(2, i, j) = Bmod(2, i, j) + &
                                                norm2(rt%n(2, :))*0.5_hp*glws(k)*MonBF(i-1, x)*MonBF(j-1, x)
                                        x = (1-lam) * rt%x(3, :) + lam*rt%x(1, :)
                                        Bmod(3, i, j) = Bmod(3, i, j) + &
                                                norm2(rt%n(3, :))*0.5_hp*glws(k)*MonBF(i-1, x)*MonBF(j-1, x)
                                end do
                        end do
                end do

                if(debug) write(*,*) "Split Surface Matrix 1 is", Bmod(1, :, :)
                if(debug) write(*,*) "Split Surface Matrix 2 is", Bmod(2, :, :)
                if(debug) write(*,*) "Split Surface Matrix 3 is", Bmod(3, :, :)

                allocate(rt%B(3, rt%ncol, rt%Nbdt))
                ! Translate into the nodal basis
                do i=1, 3
                        Bnod(i, :,:) =  matmul(transpose(rt%Vi), matmul(Bmod(i, :, :), rt%Vi))
                end do

                if (debug) write(*,*) "Split surface Matrix 1 is", Bnod(1, :, :)

                ! Only save needed columns
                do i=1, 3
                        do j=1, rt%Nbdt
                                rt%B(i, :, j) = Bnod(i, :, rt%Bi(j, i))
                        end do
                end do
        
                ! We need a surface quadrature to integrate the entropy inequality predictors
                ! summing over the surface matrices suffices to get one.
                allocate(rt%wsurf(rt%Nbdt))
                rt%wsurf = 0.0_hp
                do i=1, rt%Nbdt
                       rt%wsurf(i) = sum(rt%B(1, :, i))
                end do

                if (debug) write(*,*) "using the surface quadrature ", rt%wsurf
                
        end subroutine SurfMatSplit

         ! Calculates the surface matrix elements
        ! In the modal basis and the MFSBP setting
        ! M: Describes the DOP basis
        ! Bx, By: The surface operators in x and y direction
        ! rt: the reference triangle
        ! n: number of nodes
        subroutine SurfMat(rt)
                implicit none
                type(hprt) :: rt
                ! Aux
                real(hp), dimension(rt%ncol, rt%ncol) :: Bxmod, Bymod
                real(hp), dimension(rt%nbdt) :: gln, glws             ! A one-dimensional GLE quadrature
                real(hp), dimension(3*rt%nbdt, 2) :: sqn                ! A surface quadrature, derived from the GLE quad.
                real(hp), dimension(3*rt%nbdt, 2) :: sqw                ! And the corresponding weights.
                real(hp), dimension(2) :: sqr
                real(hp) :: lam
                integer i, j, k

                ! Calculate a surface quadrature of high enough order from a GL quadrature
                gln = GLE(rt%nbdt)
                glws = GLW(rt%nbdt, gln)
                do i=1, rt%nbdt
                        lam = 0.5_hp*(gln(i) + 1.0_hp)
                        ! The nodes of the first, second and third side
                        do j=1, 3
                                sqn((j-1)*rt%nbdt + i, :) = RT%x(j, :) + lam*(RT%x(mod(j, 3)+1, :) - RT%x(j, :))
                        end do
                        ! Now the weights
                        do j=1, 3
                                sqw((j-1)*rt%nbdt + i, :) = RT%n(j, :)*glws(i)*0.5_hp
                        end do
                end do

                ! Using this quadrature the matrix elements can be calculated.
                do i=1, rt%ncol
                        do j=1, rt%ncol
                                sqr = 0.0_hp
                                do k=1,3*rt%nbdt
                                       sqr = sqr + sqw(k, :)*MonBF(i-1, sqn(k, :))*MonBF(j-1, sqn(k, :))
                                end do
                                Bxmod(i, j) = sqr(1)
                                Bymod(i, j) = sqr(2)
                        end do
                end do
                allocate(rt%Br(rt%ncol, rt%ncol))
                allocate(rt%Bs(rt%ncol, rt%ncol))
                rt%Br = matmul(transpose(rt%Vi), matmul(Bxmod, rt%Vi))
                rt%Bs = matmul(transpose(rt%Vi), matmul(Bymod, rt%Vi))

        end subroutine SurfMat

        ! Find a G using the approach given in part II
        subroutine FindDispG(rt)
                type(hprt) :: rt
                ! Aux
                real(hp), dimension(size(rt%M, 1), size(rt%M, 2)) :: G
                real(hp), dimension(size(rt%M, 1), size(rt%M, 2)) :: Q
                real(hp) :: tlow = 0.0_hp, tup = 10.0_hp
                integer :: i, itmax = 32

                ! Calculate the heat equation discretisation on the element
                !Q = matmul(rt%Sr, matmul(rt%Mi,transpose(rt%Sr))) &
                 !       + matmul(rt%Ss, matmul(rt%Mi, transpose(rt%Ss)))
                Q = matmul(rt%Mi, rt%S2)

                if (debug) then
                        write (*,*) "Q Matrix is"
                        do i=1, size(Q, 1)
                                write (*,"(10f8.2)") Q(i, :)
                        end do
                end if


                ! Search for the minimal time until the matrix elements
                ! in the resulting oeprator are all non-negative
                ! The bisection is carried out up to itmax times giving an
                ! answer exact up to 2^(-32).
                do i=1, itmax
                        G = MatExp(-0.5_hp*(tlow + tup)*Q)

                        if (all(G .GE. real(-1.0E-8_hp, hp))) then
                                tup = 0.5_hp*(tlow + tup)
                        else
                                tlow = 0.5_hp*(tlow + tup)
                        end if
                end do
                write (*,'(A, 2F10.5)') "Heat equation positivity was reached in intervall", tlow, tup



                ! In the case of, for example, linear elements, the matrix is always positive as
                ! Q is always a positive generator itself. We use Q directly in this case.
                if ((tlow + tup) < real(1.0E-6, hp)) then
                        write (*,*) "G Matrix is always positive. Writing out discretisation itself instead"
                        tup = 0.0_hp
                        do i=1, size(Q, 1)
                                tup = max(tup, abs(Q(i, i)))
                        end do
                        G = -Q/max(norm2(Q), tup)
                        do i=1, size(Q, 1)
                                G(i, i) = G(i, i) + 1.0_hp
                        end do
                        G = MatExp(-Q)
                else
                        G = MatExp(-0.5_hp*(tlow + tup)*Q)
                end if

                if (debug) then
                        write (*,*) "G Matrix is"
                        do i=1, size(G, 1)
                                write (*,"(10f8.3)") G(i, :)
                        end do
                        write (*,*) " "
                end if 
                ! To convert between a generator and the time-step operator
                ! the identity must be subtracted
                do i=1, size(G, 1)
                        G(i, i) = G(i, i) - 1.0_hp
                end do
                allocate(rt%G(size(G, 1), size(G, 2)))
                rt%G = G
        end subroutine

        ! A new version - uses least squares interpolation.
        subroutine FindLO(rt)
                type(hprt) :: rt
                ! aux
                integer, parameter :: nlo(9) = (/3, 6, 6, 10, 15, 21, 28, 36, 45/)
                integer :: i, j, n, k
                real(hp) :: V(size(rt%M, 1), size(rt%M, 2))
                real(hp) :: VTV(size(rt%M, 1), size(rt%M, 2))
                real(hp) :: VTVI(size(rt%M, 1), size(rt%M, 2))
                
                n = size(rt%M, 1)

                k = nlo(rt%Nbdt-1)
                write (*,*) "We are using ", k, " basis functions for the projection hardening"
                do i=1, n
                        do j=1, k
                                V(i, j) = MonBF(j-1, rt%cpts(i, :))
                        end do
                end do

                VTV(1:k, 1:k) = matmul(transpose(V(1:n, 1:k)), V(1:n, 1:k))
                VTVI(1:k, 1:k) = invQR(VTV(1:k, 1:k))
                allocate(rt%LO(n, n))
                rt%LO = matmul(V(1:n, 1:k), matmul(VTVI(1:k, 1:k), transpose(V(1:n, 1:k))))

                if (debug) then
                        write (*,*) "The resulting matrix is"
                        do i=1, size(rt%M, 1)
                                write (*,"(10f8.3)") rt%LO(i, :)
                        end do
                end if

        end subroutine

        ! The main driver routine for the generation of the DG operators.
        ! rt is the reference triangle, q the order of the operators
        subroutine DGinit(rt, q)
                type(nprt) :: rt
                integer :: q
                ! Aux
                type(hprt) :: rthp
                call RefInitDG(rthp, q)
                call VolMat(rthp)
                call SurfMat(rthp)
                call SurfMatsplit(rthp)
                allocate(rthp%w(size(rthp%cpts, 1)))
                call FindQuad(rthp%cpts, rthp%w)
                call FindLO(rthp)
                call FindDispG(rthp)

                if (debug) then
                        write(*,*) "Operator construction complete. Sanity checks...."
                        call TestOp(rthp)
                end if
                call ConvertDown(rt, rthp)
        end subroutine DGinit
end module opshp
