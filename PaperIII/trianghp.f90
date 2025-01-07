! A reimplementation of the triangle module in high precision.
module trianghp
        use param
        use bf

        implicit none
        type hprt       ! The reference triangle in high precision arithmetic
                real(hp), allocatable :: V(:, :), Vi(:, :)                        ! The Vandermonde and inverse Vandermonde matrix
                real(hp), allocatable :: M(:, :), Sr(:, :), Ss(:, :), Mi(:, :)    ! Mass, Stiffness and inverse Mass matrix
                real(hp), allocatable :: S2(:, :)                                 ! Second derivative stiffness matrix.
                real(hp), allocatable :: Br(:, :), Bs(:, :)                       ! MFSBP type boundary matrices
                real(hp), allocatable :: B(:, :, :)                               ! Split boundary matrix
                real(hp), allocatable :: G(:, :)                                  ! Dissipation Operator
                real(hp), allocatable :: Lo(:, :)                                 ! Lower order projection operator
                real(hp), allocatable :: w(:)                                     ! Positive volume quadrature
                real(hp), allocatable :: cpts(:, :)                               ! Collocation points
                real(hp), dimension(3, 2) :: x                                    ! Vertices of the reference triangle
                real(hp), dimension(3, 2) :: n                                    ! normal vectors
                real(hp), allocatable :: wsurf(:)                                 ! A quadrature for the edge
                integer :: Nbdt, Ncol, q                                          ! Number of collocation points per edge, of collocation points, polynomial degree.
                integer, allocatable :: Bi(:, :)                                  ! Bi(i, j) is Node/Vector index of i-th node on j-th boundary

        end type hprt
                
        type nprt       ! The reference triangle in normal precision arithmetic
                real(np), allocatable :: V(:, :), Vi(:, :)                        ! The Vandermonde and inverse Vandermonde matrix
                real(np), allocatable :: M(:, :), Sr(:, :), Ss(:, :), Mi(:, :)    ! Mass, Stiffness and inverse Mass matrix
                real(np), allocatable :: S2(:, :)                                 ! Second derivative stiffness matrix.
                real(np), allocatable :: Br(:, :), Bs(:, :)                       ! MFSBP type boundary matrices
                real(np), allocatable :: B(:, :, :)                               ! Split boundary matrix
                real(np), allocatable :: G(:, :)                                  ! Dissipation oeprator
                real(np), allocatable :: Lo(:, :)                                 ! Lower order projection operator
                real(np), allocatable :: w(:)                                     ! Positive volume quadrature
                real(np), allocatable :: cpts(:, :)                               ! Collocation points
                real(np), dimension(3, 2) :: x                                    ! Vertices of the reference triangle
                real(np), dimension(3, 2) :: n                                    ! normal vectors
                real(np), allocatable :: wsurf(:)                                 ! A quadrature for the edge
                integer :: Nbdt, Ncol, q                                          ! Number of collocation points per edge, of collocation points, polynomial degree.
                integer, allocatable :: Bi(:, :)                                  ! Bi(i, j) is Node/Vector index of i-th node on j-th boundary

        end type nprt

        

        contains

        ! Tests the high precision operators for the reference triangle rt
        subroutine TestOp(rt)
                implicit none
                type(hprt) :: rt
                ! Auxillary
                real(hp) :: V(rt%ncol, rt%ncol), Vx(rt%ncol, rt%ncol), Vy(rt%ncol, rt%ncol)
                real(hp) :: Gx(rt%ncol, rt%ncol), Gy(rt%ncol, rt%ncol)
                real(hp) :: Dx(rt%ncol, rt%ncol), Dy(rt%ncol, rt%ncol)
                integer :: i, j

                write(*,*) "Mass Matrix is"
                do i=1, rt%ncol
                        write(*, '(10F10.4)') rt%M(i, :)
                end do

                write(*,*) "Inverse Mass Matrix is"
                do i=1, rt%ncol
                        write(*, '(10F10.4)') rt%Mi(i, :)
                end do

                write(*,*) "First Boundary Operator"
                do i=1, size(rt%B, 2)
                        write(*, '(10F10.4)') rt%B(1, i, :)
                end do
                
                write(*,*) "Second Boundary Operator"
                do i=1, size(rt%B, 2)
                        write(*, '(10F10.4)') rt%B(2, i, :)
                end do

                write(*,*) "Third Boundary Operator"
                do i=1, size(rt%B, 2)
                        write(*, '(10F10.4)') rt%B(3, i, :)
                end do 

                write(*,*) "Dissipation Generator"
                do i=1, rt%ncol
                        write(*, '(10F10.4)') rt%G(i, :)
                end do

                write(*,*) "Low order Projection"
                do i=1, rt%ncol
                        write(*, '(10F10.4)') rt%Lo(i, :)
                end do
                
                write(*,*) "Volume quadrature"
                write(*, '(10F10.4)') rt%w

                write(*,*) "Surface Quadrature"
                write(*,'(10F10.4)') rt%wsurf

                do i=1, rt%ncol
                        do j=1, rt%ncol
                                V(i, j) = MonBF(j-1, rt%cpts(i, :))
                                Vx(i, j) = MonBFx(j-1, rt%cpts(i, :))
                                Vy(i, j) = MonBFy(j-1, rt%cpts(i, :))
                        end do        
                end do
                
                Dx = matmul(rt%Mi, transpose(rt%Sr))
                Dy = matmul(rt%Mi, transpose(rt%Ss))

                write (*,*) "Testing SBP property in x"
                do i=1,rt%ncol
                        write (*,"(10e10.2)") (rt%Sr(i,:) + rt%Sr(:, i) - rt%Br(i, :))
                end do
                write (*,*) "Testing SBP property in y"
                do i=1,rt%ncol
                        write (*, "(10e10.2)") (rt%Ss(i,:) + rt%Ss(:, i) - rt%Bs(i, :))
                end do
        
                Gx = matmul(Dx, V)
                Gy = matmul(Dy, V)

                write (*,*) "Testing Derivative exactness in x"
                do i=1,rt%ncol
                        write (*,"(10e10.2)") (Vx(i, :) - Gx(i, :))
                end do
                write (*,*) "Testing Derivative exactness in y"
                do i=1, rt%ncol
                        write (*,"(10e10.2)") (Vy(i, :) - Gy(i, :))
                end do
                write (*,*) "Operatornorm of the derivative error is", norm2(Vx - Gx) + norm2(Vy - Gy)

        end subroutine TestOp

        ! This subroutine converts a reference triangle in high precision to normal precision
        subroutine ConvertDown(rto, rti)
                type(nprt) :: rto
                type(hprt) :: rti
                
                ! copy fixed size parts
                rto%x = real(rti%x, np)
                rto%n = real(rti%n, np)
                rto%nbdt = rti%nbdt
                rto%ncol = rti%ncol
                rto%q = rti%q

                ! allocate everything allocatable
                allocate(rto%V(rti%ncol, rti%ncol))
                allocate(rto%Vi(rti%ncol, rti%ncol))

                allocate(rto%M(rti%Ncol, rti%Ncol))
                allocate(rto%Sr(rti%ncol, rti%ncol))
                allocate(rto%Ss(rti%ncol, rti%ncol))
                allocate(rto%Mi(rti%ncol, rti%ncol))
                allocate(rto%S2(rti%ncol, rti%ncol))
                allocate(rto%Br(rti%ncol, rti%ncol))
                allocate(rto%Bs(rti%ncol, rti%ncol))
                allocate(rto%B(3, rti%ncol, rti%nbdt))
                allocate(rto%Lo(rti%ncol, rti%ncol))
                allocate(rto%G(rti%ncol, rti%ncol))
                allocate(rto%w(rti%ncol))
                allocate(rto%wsurf(rti%nbdt))
                allocate(rto%cpts(rti%ncol, 2))
                allocate(rto%Bi(rti%nbdt, 3))

                ! cast arrays down
                rto%V = real(rti%V, np)
                rto%Vi = real(rti%Vi, np)
                rto%M = real(rti%M, np)
                rto%Mi = real(rti%Mi, np)
                rto%Sr = real(rti%Sr, np)
                rto%Ss = real(rti%Ss, np)
                rto%Mi = real(rti%Mi, np)
                rto%S2 = real(rti%S2, np)
                rto%Br = real(rti%Br, np)
                rto%Bs = real(rti%Bs, np)
                rto%B = real(rti%B, np)
                rto%Lo = real(rti%Lo, np)
                rto%G = real(rti%G, np)
                rto%w = real(rti%w, np)
                rto%wsurf = real(rti%wsurf, np)
                rto%cpts = real(rti%cpts, np)
                rto%Bi = rti%Bi
        end subroutine ConvertDown

        ! Initializes the reference triangle, only sets up basic geometry
        subroutine RefInitBasic(RT)
                type(hprt) :: RT
                RT%x(1, :) = [0.0_hp, 0.0_hp]
                RT%x(2, :) = [1.0_hp, 0.0_hp]
                RT%x(3, :) = [0.0_hp, 1.0_hp]
                RT%n(1, :) = [0.0_hp, -1.0_hp]
                RT%n(2, :) = [1.0_hp, 1.0_hp]
                RT%n(3, :) = [-1.0_hp, 0.0_hp]
        end subroutine RefInitBasic

        ! Initializes a complete reference triangle for a nodal DG method
        subroutine RefInitDG(RT, q)
                type(hprt) :: RT
                integer :: q
                ! aux
                !integer :: q
                integer :: n, i, j, k, l
                real(hp) :: lam
                real(hp) :: lam1, lam2, lam3
                ! aux
                integer, parameter :: Nmax = 20
                real(hp) :: glopts(Nmax)
                ! Call the geometry init
                call RefInitBasic(RT)
                if (OpDesign == "CG1") then             ! The linear CG triangle, one node per vertex.
                        q = 1                           ! polynomial degree one.
                        n = 3                           ! 3 collocation nodes.
                        RT%Nbdt = 2                     ! 2 nodes per edge
                        allocate(RT%cpts(n, 2))
                        RT%cpts(:, :) = RT%x            ! the vertices are the nodes.
                        allocate(RT%Bi(RT%Nbdt, 3))
                        do i=1, rt%Nbdt
                                rt%Bi(i, 1) = i
                                rt%Bi(i, 2) = i + 1
                                rt%Bi(i, 3) = mod(i+1, n)+1
                        end do
                else if (OpDesign == "CG2") then        ! Quadratic CG triangles
                        q = 2                           ! polynomial degree two.
                        n = 6                           ! 6 collocation nodes.
                        rt%Nbdt = 3                     ! 3 nodes per triangle edge
                        allocate(RT%cpts(n, 2))
                        do i=1, 3
                                rt%cpts(2*(i-1) + 1, :) = RT%x(i,:)
                                rt%cpts(2*i, :) = 0.5_hp*(RT%x(i, :) + RT%x(mod(i-1, 3)+2, :))
                        end do
                        allocate(RT%Bi(RT%Nbdt, 3))
                        do i=1, rt%Nbdt
                                rt%Bi(i, 1) = i
                                rt%Bi(i, 2) = i + rt%Nbdt-1
                                rt%Bi(i, 3) = mod(i+rt%Nbdt, n)+1

                        end do
                else if (OpDesign == "CGgl") then       ! GL triangle nodes (orig order 3)
                        glopts(1:q+1) = GLO(q+1)
                        glopts = 0.5_hp*(glopts + 1.0_hp)
                        
                        n = (q+1)*(q+2)/2               ! 15 Collocation nodes
                        rt%nbdt = q+1                     ! 5 Nodes per edge, 3 Nodes in the interrior
                        allocate(Rt%cpts(n, 2))
                        do i=1, 3
                                do j=1, q               ! Add the surface nodes on all three edges
                                        lam = glopts(j)
                                        rt%cpts((rt%Nbdt-1)*(i-1) + j, :) = lam*RT%x(mod(i, 3)+1,:) + (1-lam)*RT%x(i, :)
                                end do
                        end do
                        k = 3*(rt%nbdt-1) + 1


                        do i=1, rt%nbdt                       ! The inner nodes lie in a "triangular" pattern
                                do j=1, rt%nbdt+1-i
                                        l = rt%nbdt + 2 - i - j
                                        lam1 = (1.0_hp + 2.0_hp*glopts(i) - glopts(j) - glopts(l)) / 3.0_hp
                                        lam2 = (1.0_hp - glopts(i) +  2.0d0*glopts(j) - glopts(l)) / 3.0_hp
                                        lam3 = 1.0_hp - lam1 - lam2
                                        write(*, '(A, 2I3, A, 2F10.5)') "Collocation point with index", i, j,&
                                                         " and coordinates ", lam1, lam2
                                        if (i .NE. 1 .AND. j .NE. 1 .AND. i+j .NE. rt%nbdt+1) then  ! We are not on the surface
                                                write (*,*) "Its is an inner point"
                                                rt%cpts(k, :) = lam1*rt%x(2, :) + lam2*rt%x(3, :) + lam3*rt%x(1, :)
                                                k = k + 1
                                        end if
                                end do
                        end do
                        allocate(RT%Bi(rt%Nbdt, 3))
                        do i=1, rt%Nbdt
                                rt%Bi(i, 1) = i
                                rt%Bi(i, 2) = i + rt%Nbdt-1
                                rt%Bi(i, 3) = mod(i+2*(rt%Nbdt-1)-1, 3*(rt%nbdt-1))+1
                                write (*,'(A, I3, A, I3, A, I3, A, I3)') "Node", i, " on side 1 has index", rt%Bi(i, 1), &
                                        " on side 2 has index", rt%Bi(i, 2), " on side 3 has index", rt%Bi(i, 3)

                        end do

                else
                        write(*,*) "implement DG reference triangle", OpDesign, " !"
                end if
                rt%ncol = n
                rt%q = q
        end subroutine RefInitDG


end module trianghp
