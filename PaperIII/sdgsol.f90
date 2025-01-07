! This file constains the polygonal surface split DG solver
! that is stabilized using the Entropy Rate Criterion.

module sdgsolver
        use param
        use euler       ! Inclusion of needed modules
        use trianggrid
        use trianghp

        implicit none
        !logical, parameter :: projdown = .false.
        logical, parameter :: projdown = .true.

        contains

        ! Estimates the highest signal speed in the domain to control
        ! time stepping
        ! gr: Description of the Grid
        ! rt: Description of the Reference Triangle
        ! u: Array of conserved variable nodal values
        !       first dim: triangle
        !       second dim: collocation node
        !       last dim: conserved variable
        ! bf(uo, u, x, t) : the boundary function, 
        !       taking an internal and external state and a position
        ! t: the time at which $u$ was captured
        ! return value: the speed estimate
        function SpeedEstSplit(gr, rt, u, bf, t)
                type(Grid) :: gr
                type(nprt) :: rt
                double precision, intent(in) :: u(:, :, :)
                external :: bf
                real(np) :: t
                ! aux
                real(np) :: SpeedEstSplit
                real(np) :: ui(size(u, 3)), uo(size(u, 3))
                real(np) :: tse(size(u, 1))
                real(np) :: v(2), n(2), x(2)
                integer :: i, j, k, neib, AdjE

                !$OMP PARALLEL DO PRIVATE(j, k, v, n, x, neib, ui, uo, AdjE) SHARED(tse)
                do i=1, size(u, 1)
                        ! the speed bound is first calculated per triangle (triangle speed estimate)
                        tse(i) = 0.0_np
                        ! by looking at every surface and internal node
                        ! surface nodes
                        do j=1,3
                                neib = gr%adj(i, j)
                                do k=1, rt%Nbdt
                                        ui = u(i, rt%Bi(k, j), :)
                                        if (neib > 0) then ! Internal Edge
                                                AdjE = gr%AdjO(i, j)
                                                uo = u(neib, rt%Bi(rt%Nbdt - k + 1, AdjE), :)
                                        else if (neib .EQ. -1) then ! Reflection
                                                n = gr%nn(i, j, :)
                                                uo = MirrorConsVar(ui, n)
                                        else    ! Outer edge
                                                x = GetNodeKord(gr, rt, i, rt%Bi(k, j))
                                                call bf(uo, ui, x, t)
                                        end if
                                        tse(i) = max(tse(i), cmax(ui, uo))
                                end do

                        end do

                        ! internal nodes
                        do j=(rt%q+1)*3+1, size(u, 2)
                                v = u(i, j, 2:3)/u(i, j, 1)
                                tse(i) = max(tse(i), a(u(i, j, :)) + norm2(v))
                        end do
                end do
                !$OMP END PARALLEL DO
                ! The global estimate is then calculated via a standard library routine
                SpeedEstSplit = maxval(tse)
        end function SpeedEstSplit

        ! stable division with rounding error eps.
        function sdiv(a, b, eps)
                real(np) :: a, b, eps
                real(np) :: sdiv
                sdiv = a*b/(b**2 + eps**2)
        end function sdiv


        ! Main discrectisation routine of the stabilized DG method.
        ! gr: description of the grid
        ! rt: the reference triangle
        ! du: calculated time derivative of the conserved variable u
        ! u: tensor of conserved variables
        !       first dim: triangle
        !       second dim: node in the triangle
        !       last dim: conserved quantity
        ! bf: boundary function that should be enforced
        ! t: time
        ! dt: timestep size, needed for stability of the correction direction added.
        subroutine DSDG(gr, rt, du, u, bf, t, dt)
                ! in and aoutput
                type(Grid) :: gr
                type(nprt) :: rt
                real(np), intent(out) :: du(:, :, :)
                real(np), intent(in) :: u(:, :, :)
                external :: bf
                real(np), intent(in) :: t, dt
                real(np), parameter :: eps1 = 1.0D-7                            ! The assumed precision, sqrt of machine epsilon
                ! auxillary variables
                real(np), dimension(size(gr%edges, 1)) :: sigma, sigmaLO        ! Total entropy dissipation pred. per edge
                real(np), dimension(size(u, 1)):: lam                           ! disp scaling factor
                real(np), dimension(size(u, 1), size(u, 2), size(u, 3)) :: Gu   ! Disp directions
                real(np), dimension(size(u, 1)) :: dE, dEcorr                   ! Entropy dissipation needed and entropy dissipation of the direction
                real(np), dimension(size(u, 1)) :: Eflow                        ! Entropy change due to entropy flow
                real(np), dimension(size(u, 3)) :: ui, uo, uiLO, uoLO           ! inner and outer value of conserved variables
                real(np), dimension(size(u, 1), size(u, 2), size(u, 3)):: uLO   ! low order projection of the conserved variables
                real(np), dimension(size(u, 1)) :: Eb, Ec                       ! Entropy of the basis solution and of the correction.
                real(np), dimension(size(u, 1), size(u, 2), size(u, 3)) :: Upsu ! Solution after application of the positive conservative filter. 
                real(np), dimension(2) :: n, x                                  ! normal and local position
                real(np) :: ced, bsed, edgelam                                  ! The combined entropy dissipation of two cells
                real(np) :: eps                                                 ! the scaled machine precision
                real(np), parameter :: delta = 0.5_np                              ! the stepsize for the FD disp. calc.
                integer :: i, j, k                              ! counting indices
                integer ::  Nt, edge, adjE, Ne                  ! Number of triangles, local edge number (1-3), corresponding edge in the other triangle, # of edges.
                Nt = size(u, 1)
                Ne = size(gr%edges, 1)

                ! Calculation of the uncorrected time derivative using a DG scheme 
                call SDG(gr, rt, du, u, bf, t, dt)

                ! The first loop below calculates the flow of entropy, the correction direction,
                ! the low order projection, the effectivity of the correction direction per 
                ! triangle and the correction needed to be entropy dissipative.

                !$OMP PARALLEL DO PRIVATE(n, j, k, ui, uo, x) SHARED(dEcorr, Gu, dE)
                do i=1, Nt
                        ! Calculation of a scaled machine precision.
                        !eps = eps1*gr%Vol(i)
                        eps = eps1*gr%bl(i, 1)
                        !eps = eps1
                        ! Projection onto lower order polynomials
                        do k =1, size(u, 3)
                                uLO(i, :, k) = matmul(rt%LO, u(i, :, k))
                        end do
                        ! Calculation of the correction direction Gu and the dissipation speed of Gu
                        Gu(i, :, :) = matmul(rt%G, u(i, :, :))

                        ! A finite difference approximation is used
                        upsu(i, :, :) = u(i, :, :) + delta*Gu(i, :, :)
                        Eb(i) = 0.0_np
                        Ec(i) = 0.0_np
                        ! Eb holds the original entropy, Ec after the application of 0.1 G
                        do k=1, size(u, 2)
                                Eb(i) = Eb(i) + rt%w(k)*PUEuler(u(i, k, :))
                                Ec(i) = Ec(i) + rt%w(k)*PUEuler(upsu(i, k, :))
                        end do
                        Eb(i) = Eb(i)*gr%Vol(i)*2.0_np
                        Ec(i) = Ec(i)*gr%Vol(i)*2.0_np

                        ! This yields the following entropy derivative bound from the application of G 
                        dEcorr(i) = (Ec(i) - Eb(i)) / delta
                        ! Calculation of the entropy equality violation by the uncorrected time derivative
                        dE(i) = 0.0_np
                        ! Volume integral of the total entropy change in cell i.
                        do j=1, size(rt%w)
                                dE(i) = dE(i) + rt%w(j)*dot_product(dPUEuler(u(i, j, :)), du(i, j, :))
                        end do
                        dE(i) = dE(i)*gr%Vol(i)*2
                        ! Surface flux integral
                        Eflow(i) = 0.0_np
                        do j=1, 3
                                n = gr%nn(i, j, :)
                                do k=1, rt%Nbdt
                                        ui = u(i, rt%Bi(k ,j), :)
                                        if (gr%Adj(i, j) > 0) then ! Internal Edge
                                                uo = u(gr%Adj(i, j), rt%Bi(rt%Nbdt -k + 1, gr%AdjO(i, j)), :)
                                        else if (gr%Adj(i, j) .EQ. -1) then ! Reflection
                                                uo = MirrorConsVar(ui, n)
                                        else    ! Outer edge
                                                x = GetNodeKord(gr, rt, i, rt%Bi(k ,j))
                                                call bf(uo, ui, x, t)
                                        end if
                                        ! The entropy change from entropy flow can be calculated via integration over the cell boundary
                                        Eflow(i) = Eflow(i) - rt%wsurf(k)*gr%bl(i, j)*PEulerHLL(ui, uo, n)
                                end do
                        end do
                        ! Calculation of the needed dissipation to get entropy dissipative, the difference between change from entropy
                        ! flow and the total entropy change has to be taken into account.
                        lam(i) = max(0.0_np, dE(i) - Eflow(i))/(eps - dEcorr(i))
                        ! We update dE as the negative part of the entropy equality violation. 
                        dE(i) = min(0.0_np, dE(i) - Eflow(i))
                end do
                !$OMP END PARALLEL DO
                bsed = 0.0_np
                bsedt = 0.0_np

                ! This second loop is used to calculate the entropy inequality predictor output per edge

                !$OMP PARALLEL DO PRIVATE(edge, AdjE, n, ui, uo, uiLO, uoLO, k, ced, EdgeLam, bsed) SHARED(sigma, sigmaLO, lam)
                do i=1, Ne
                        ! This if clause detects if the edge is an inner, reflecting or outer edge and integrates the EIP
                        ! accordingly.
                        if (gr%edges(i, 2) > 0) then            ! An inner edge
                                ! Use the scaling of the first triangle on the edge and the normal pointing to the neighbor
                                ! to calculate the corresponding normal
                                edge = gr%edgeNumb(i,1)
                                AdjE = gr%edgeNumb(i,2)
                                n = gr%nn(gr%edges(i, 1), edge, :)
                                sigma(i) = 0.0_np
                                sigmaLO(i) = 0.0_np
                                ! Integrate the entropy inequality predictor over the edge, also in the low order version
                                do k=1, rt%Nbdt
                                        ui = u(gr%edges(i, 1), rt%Bi(k, edge), :)
                                        uiLO = uLO(gr%edges(i, 1), rt%Bi(k, edge), :)
                                        uo = u(gr%edges(i, 2), rt%Bi(rt%Nbdt - k + 1, adjE), :)
                                        uoLO = uLO(gr%edges(i, 2), rt%Bi(rt%Nbdt - k + 1, adjE), :)
                                        sigma(i) = sigma(i) + rt%wsurf(k)&
                                                *EulerDispOpHLL(ui, uo, n)*gr%bl(gr%edges(i, 1),edge)
                                        sigmaLO(i) = sigmaLO(i) + rt%wsurf(k)*&
                                                EulerDispOpHLL(uiLO, uoLO, n)*gr%bl(gr%edges(i, 1),edge)
                                end do
                        else if (gr%edges(i, 2) .EQ. -1) then ! A reflecting boundary. Everything as before, but with reflected states.
                                edge = gr%edgeNumb(i, 1)
                                n = gr%nn(gr%edges(i, 1), edge, :)
                                sigma(i) = 0.0_np
                                sigmaLO(i) = 0.0_np
                                do k=1, rt%Nbdt
                                        ui = u(gr%edges(i, 1), rt%Bi(k, edge), :)
                                        uo = MirrorConsVar(ui, n)
                                        sigma(i) = sigma(i) + 0.5_np*rt%wsurf(k)&
                                                *EulerDispOpHLL(ui, uo, n)*gr%bl(gr%edges(i, 1),edge)
                                        uiLO = uLO(gr%edges(i, 1), rt%Bi(k, edge), :)
                                        uoLO = MirrorConsVar(uiLO, n)
                                        sigmaLO(i) = sigmaLO(i)+0.5_np*rt%wsurf(k)*&
                                                EulerDispOpHLL(uiLO, uoLO, n)*gr%bl(gr%edges(i, 1),edge)
                                end do
                        else if (gr%edges(i, 2) .EQ. -2) then ! A controled boundary. I.e. with the external state from the boundary function.
                                edge = gr%edgeNumb(i, 1)
                                n = gr%nn(gr%edges(i, 1), edge, :)
                                sigma(i) = 0.0_np
                                sigmaLO(i) = 0.0_np
                                do k=1, rt%Nbdt
                                        ui = u(gr%edges(i, 1), rt%Bi(k, edge), :)
                                        x = GetNodeKord(gr, rt, gr%edges(i, 1), rt%Bi(k ,edge))
                                        call bf(uo, ui, x, t)
                                        sigma(i) = sigma(i) + 0.5_np*rt%wsurf(k)&
                                                *EulerDispOpHLL(ui, uo, n)*gr%bl(gr%edges(i, 1),edge)
                                        uiLO = uLO(gr%edges(i, 1), rt%Bi(k, edge), :)
                                        call bf(uoLO, uiLO, x, t)
                                        sigmaLO(i) = sigmaLO(i) + 0.5_np*rt%wsurf(k)*EulerDispOpHLL(uiLO, uoLO, n)&
                                                                        *gr%bl(gr%edges(i, 1),edge)
                                end do

                        end if

                        ! The prediction per edge is the minimum of the low order and high order version.
                        if (projdown) then
                                sigma(i) = min(sigma(i), sigmaLO(i))
                        end if

                        
                        
                        ! This if clause splits the entropy dissipation mandated by edge i onto the up to two triangles
                        ! adjacent to it.
                        if (gr%edges(i, 2) > 0) then                    ! An inner edge
                                eps = eps1*min(gr%bl(gr%edges(i, 1), 1), gr%bl(gr%edges(i, 2), 1))
                                ! Calculation with the first entropy
                                ced = dEcorr(gr%edges(i, 1))
                                bsed = 0.5_np*(dE(gr%edges(i, 1)) + dE(gr%edges(i, 2)))
                                EdgeLam = min(0.0_np, (sigma(i) - bsed))/(ced - eps)
                                lam(gr%edges(i, 1)) = lam(gr%edges(i,1)) + EdgeLam
                                lam(gr%edges(i, 2)) = lam(gr%edges(i,2)) + EdgeLam

                        else if (gr%edges(i, 2) .EQ. -1) then           ! Reflection
                                eps = eps1*gr%bl(gr%edges(i, 1), 1)
                                ced = dEcorr(gr%edges(i, 1))
                                bsed = 0.5_np*dE(gr%edges(i, 1))
                                EdgeLam = min(0.0_np, (sigma(i) - bsed))/(ced - eps)
                                lam(gr%edges(i, 1)) = lam(gr%edges(i,1)) + EdgeLam
                        else if (gr%edges(i, 2) .EQ. -2) then           ! Controled boundary
                                eps = eps1*gr%bl(gr%edges(i, 1), 1)
                                ced = dEcorr(gr%edges(i, 1))
                                bsed = 0.5_np*dE(gr%edges(i, 1))
                                EdgeLam = min(0.0_np, (sigma(i) - bsed))/(ced - eps)
                                lam(gr%edges(i, 1)) = lam(gr%edges(i,1)) + EdgeLam
                        end if

                end do
                !$OMP END PARALLEL DO

                ! The last loop of the function adds the calculated amounts of dissipation onto the time derivative.
                !$OMP PARALLEL DO SHARED(du, lam, GU)
                do i=1, Nt
                        lam(i) = min(lam(i), 0.5_np/dt)
                        du(i, :, :) = du(i, :, :) + lam(i)*Gu(i, :, :)
                end do
                !$OMP END PARALLEL DO
        end subroutine DSDG
                
        ! A unstructured, edge split implementation of the Nodal Discontinuous Galerkin Method+
        ! on triangular grids.
        ! gr: describes the grid
        ! rt: holds the reference triangle
        ! du: time derivative of the conserved variables
        ! u: conserved variables, 3-Tensor
        !       First dim: triangle index
        !       second dim: node index in the triangle
        !       last dim: number of conserbed variable
        ! bf(x, t, ui, uo): the boundary function
        ! t: time of u
        ! dt: timestep. 
        subroutine SDG(gr, rt, du, u, bf, t, dt)
                type(Grid) :: gr
                type(nprt) :: rt
                real(np), intent(out) :: du(:, :, :)
                real(np), intent(in) :: u(:, :, :)
                external :: bf
                real(np) :: t, dt
                ! auxillary variables
                real(np) :: n(2), x(2)                          ! Normal and position
                real(np), dimension(size(u, 3)) :: uo                   ! outer value of conserved variables
                real(np), dimension(rt%Nbdt, size(u, 3)) :: fi          ! Intercell flux
                real(np), dimension(size(u, 2), size(u, 3)) :: f, g    ! the nodal values of the flux at collocation points
                real(np), dimension(size(u, 2), size(u, 3)) :: h        ! time derivative before multiplication with inverse mass matrix
                integer :: Nt, Nc, k, l, m, neib, adjE
                Nt = size(u, 1)
                Nc = size(u, 2)

                ! This loop iterates over all triangles and calculates the time derivative of their local degrees of freedom.
                !$OMP PARALLEL DO PRIVATE(n, x, uo, fi, f, g, h, k, l, m, neib, adjE) SHARED(du, u, gr, rt, t) 
                do k=1, Nt
                        ! First, we calculate the flux values at the nodes.
                        do l=1, Nc
                                f(l, :) = fEuler(u(k, l, :))
                                g(l, :) = gEuler(u(k, l, :))
                        end do
                        ! And multiply with the stiffness matrices for the volume parts of the time derivative
                        h = matmul(rt%Sr, f)*gr%ToRef(k, 1, 1) + matmul(rt%Ss, f)*gr%ToRef(k, 2, 1) &
                         + matmul(rt%Sr, g)*gr%ToRef(k, 1, 2) + matmul(rt%Ss, g)*gr%ToRef(k, 2, 2)

                        ! Then we add the boundary terms from the surface fluxes by iteraing over the 3 edges of the triangle.
                        do m=1, 3
                                ! The grid data structure holds the normals and adjacent triangle indices for the 3 edges.
                                n = gr%nn(k, m, :)
                                neib = gr%adj(k, m)
                                ! The calculation of the surface fluxes depends on the boundary condition.
                                ! Positive neib values denote other triangles, negative values a boundary condition.
                                if (neib > 0) then              ! It is an internal edge between two triangles
                                        adjE = gr%adjO(k, m)
                                        do l=1, rt%Nbdt                 
                                                fi(l, :) = EulerHLLnf(u(k, rt%Bi(l,m), :), u(neib, rt%Bi(rt%Nbdt-l+1, adjE), :), n)
                                        end do
                                else if (neib .EQ. -1) then     !Reflection
                                        do l=1, rt%Nbdt
                                                fi(l, :) = EulerReflHLL(u(k, rt%Bi(l,m), :), n)
                                        end do  
                                else    ! Far Field boundary
                                        do l=1, rt%Nbdt
                                                x = GetNodeKord(gr, rt, k, rt%Bi(l, m))
                                                call bf(uo, u(k, rt%Bi(l, m), :), x, t)
                                                fi(l, :) = EulerHLLnf(u(k, rt%Bi(l, m), :), uo, n)
                                        end do
                                end if
                                ! The calculated values of fi are used in the discretisation of the boundary operators by
                                ! the three B matrices of the reference triangle.
                                h = h  - gr%ToRefDet(k)*matmul(rt%B(m, :, :), fi)*gr%bl(k, m)/norm2(rt%n(m, :))
                        end do
                        ! The resulting time derivative of the inner product description of u is converted to 
                        ! a time derivative of the nodal values by multiplication with the inverse mass matrix.
                        du(k, :, :) = matmul(rt%Mi, h)
                end do
        end subroutine SDG

        ! This limiter can be called on the vector of conserved variables to cure possible violations of
        ! positivity of pressure and density. It works by checking for violations and applies the 
        ! positive conservative filter if a violation was detected. If the average values of the entire
        ! cell is still positive results a correction to a positive value, as negative pressure and 
        ! negative density are concave. cf Hillebrand, Klein, Öffner:
        !        Applications of Limiters, Neural Networks and Polynomial Annihilation in Higher-Order FD/FV
        ! gr: Grid description
        ! rt: Reference triangle description
        ! u: tensor of conserved variable nodal values.
        subroutine PosLimiter(gr, rt, u)
                implicit none
                type(Grid) :: gr
                type(nprt) :: rt
                real(np), intent(inout) :: u(:, :, :)
                ! Auxillary variables
                real(np) :: avu(size(u, 3))
                real(np) :: lam, minp, minrho
                integer :: i, j, k
                real(np), parameter :: rhobound = 1.0E-12_np
                real(np), parameter :: pbound = 1.0E-12_np

                !$OMP PARALLEL DO PRIVATE(avu, minp, minrho, lam, j, k) SHARED(u)
                do i=1, size(u, 1)
                        do k=1,size(u, 3)
                                avu(k) = dot_product(rt%w, u(i, :, k))*2.0_np
                        end do
                        minrho = minval(u(i, :, 1))
                        if (minrho < rhobound) then
                                lam = min(1.0_np, (rhobound-minrho)/(avu(1)-minrho))
                                do j=1,size(u, 2)
                                        u(i, j, :) = lam*avu + (1-lam)*u(i, j, :)
                                end do
                        end if

                        minp = p(u(i, 1, :))
                        do j=2,size(u, 2)
                                minp = min(minp, p(u(i, j, :)))
                        end do
                        if (minp < pbound) then
                                lam = min(1.0, (pbound-minp)/(p(avu)-minp))
                                do j=1,size(u, 2)
                                        u(i, j, :) = lam*avu + (1-lam)*u(i, j, :)
                                end do
                        end if
                end do
                !$OMP END PARALLEL DO
        end subroutine PosLimiter


        ! Calculates SSPRK33 step
        ! Using the Nodal DG scheme
        ! gr: Grid description
        ! rt: Reference Triangle
        ! unew: the resulting values for the conserved variables at t+dt
        ! u: the initial condition of the conserved variables at t
        ! t: time at IC
        ! dt: timestep 
        subroutine SSPRK33SDG(gr,rt, unew, u, bf, t, dt)
                implicit none
                type(Grid), intent(in) :: gr
                type(nprt), intent(in) :: rt
                real(np), allocatable, intent(out) :: unew(:, :, :)
                real(np), allocatable, intent(inout) :: u(:, :, :)
                external :: bf
                real(np) :: t, dt
                real(np), allocatable :: du(:, :, :), u1(:, :, :), u2(:, :, :)
                integer :: n, k, l
                n = size(u, 1)
                k = size(u, 2)
                l = size(u, 3)
                allocate(du(n, k, l))
                allocate(u1(n, k, l))
                allocate(u2(n, k, l))
                call SDG(gr, rt, du, u, bf, t, dt)
                u1 = u + dt*du
                call PosLimiter(gr, rt, u1)
                call SDG(gr, rt, du, u1, bf, t, dt)
                u2 = 0.75_np*u + 0.25_np*u1 + 0.25_np*dt*du
                call PosLimiter(gr, rt, u2)
                call SDG(gr, rt, du, u2, bf, t, dt)
                unew = u/3.0_np + 2.0_np/3.0_np*u2 + 2.0_np/3.0_np*dt*du
                call PosLimiter(gr, rt, unew)
        end subroutine SSPRK33SDG

        ! Calculates SSPRK33 step
        ! Using the entropy rate criterion stabilized DG scheme
        ! gr: Grid description
        ! rt: Reference Triangle
        ! unew: the resulting values for the conserved variables at t+dt
        ! u: the initial condition of the conserved variables at t
        ! t: time at IC
        ! dt: timestep 
        subroutine SSPRK33DSDG(gr,rt, unew, u, bf, t, dt)
                implicit none
                type(Grid), intent(in) :: gr
                type(nprt), intent(in) :: rt
                real(np), allocatable, intent(out) :: unew(:, :, :)
                real(np), allocatable, intent(inout) :: u(:, :, :)
                external :: bf
                real(np) :: t, dt
                real(np), allocatable :: du(:, :, :), u1(:, :, :), u2(:, :, :)
                integer :: n, k, l
                n = size(u, 1)
                k = size(u, 2)
                l = size(u, 3)
                allocate(du(n, k, l))
                allocate(u1(n, k, l))
                allocate(u2(n, k, l))
                call DSDG(gr, rt, du, u, bf, t, dt)
                u1 = u + dt*du
                call PosLimiter(gr, rt, u1)
                call DSDG(gr, rt, du, u1, bf, t, dt)
                u2 = 0.75_np*u + 0.25_np*u1 + 0.25_np*dt*du
                call PosLimiter(gr, rt, u2)
                call DSDG(gr, rt, du, u2, bf, t, dt)
                unew = u/3.0_np + 2.0_np/3.0_np*u2 + 2.0_np/3.0_np*dt*du
                call PosLimiter(gr, rt, unew)
        end subroutine SSPRK33DSDG

        ! Solves the Euler equations using a Nodal DG method up to tend with
        ! SSPRK33 time integration. Timesteping is done with respect to a grid
        ! stiffness measure and a speed estimate for the highest signal speeds.
        ! gr: Grid description
        ! rt: Reference triangle
        ! ures: resulting conserved variables at tend
        ! u0: initial condition
        ! bf: boundary function
        ! tend: final time of the simulation
        ! cfl: Courant-Friedrichs-Lewy number 
        subroutine solSDGSSPRK33(gr, rt, ures, u0, bf, tend, cfl)
                implicit none
                type(Grid), intent(in) :: gr
                type(nprt), intent(in) :: rt
                real(np), allocatable, intent(out) :: ures(:, :, :)
                real(np), allocatable, intent(inout) :: u0(:, :, :)
                external :: bf
                real(np), intent(in) :: tend, cfl
                real(np) :: dt, t = 0.0D0
                real(np) :: gs, fs
                integer :: k


                do while (t < tend)
                        gs = FindStiff(gr)
                        fs = SpeedEstSplit(gr, rt, u0, bf, t)
                        dt = cfl / (gs*fs)
                        do k=1,50
                                call SSPRK33SDG(gr, rt, ures, u0, bf, t, dt)
                                t = t + dt
                                call SSPRK33SDG(gr, rt, u0, ures, bf, t, dt)
                                t = t + dt
                        end do

                        write (*, '(f8.4, A, f8.4, X, I3, A, f8.4, A, D12.4)') &
                                t, " /", tend, int((t*100)/tend), " % Wavespeed:", fs, " dt: ", dt
                end do
                ures(:, :, :) = u0(:, :, :)
        end subroutine solSDGSSPRK33

        ! Solves the Euler equations using a stabilized Nodal DG method 
        ! up to tend with SSPRK33 time integration. 
        ! Timesteping is done with respect to a grid stiffness
        ! measure and a speed estimate for the highest signal speeds.
        ! gr: Grid description
        ! rt: Reference triangle
        ! ures: resulting conserved variables at tend
        ! u0: initial condition
        ! bf: boundary function
        ! tend: final time of the simulation
        ! cfl: Courant-Friedrichs-Lewy number 
        subroutine solDSDGSSPRK33(gr, rt, ures, u0, bf, tend, cfl)
                implicit none
                type(Grid), intent(in) :: gr
                type(nprt), intent(in) :: rt
                real(np), allocatable, intent(out) :: ures(:, :, :)
                real(np), allocatable, intent(inout) :: u0(:, :, :)
                external :: bf
                real(np), intent(in) :: tend, cfl
                real(np) :: dt = 0.0_np, t = 0.0_np
                real(np) :: gs, fs
                integer :: k = 0


                do while (t < tend-2.0_np*dt)
                        if (k < 1) then
                                gs = FindStiff(gr)
                                fs = SpeedEstSplit(gr, rt, u0, bf, t)
                                dt = cfl / (gs*fs)
                                write (*, '(f8.4, A, f8.4, X, I3, A, f8.4, A, D12.4)') &
                                 t, " /", tend, int((t*100)/tend), " % Wavespeed:", fs, " dt: ", dt
                                k = 100
                        end if 
                        call SSPRK33DSDG(gr, rt, ures, u0, bf, t, dt)
                        t = t + dt
                        call SSPRK33DSDG(gr, rt, u0, ures, bf, t, dt)
                        t = t + dt
                        k = k-1   
                end do
                dt = 0.5_np*(tend - t)
                call SSPRK33DSDG(gr, rt, ures, u0, bf, t, dt)
                t = t + dt
                call SSPRK33DSDG(gr, rt, u0, ures, bf, t, dt)
                t = t + dt

                ures(:, :, :) = u0(:, :, :)
        end subroutine solDSDGSSPRK33


end module
