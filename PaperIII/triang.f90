! The module triangles contains the routines and data structures to handle the 
! unstructured conforming triangle grids of the solver.
! Especially the creation of the reference triangles and the loading of
! grids in the output format of Jonathan Skewchuks Triangle library. 

module triangles
        use HDF5
        use ndgmats
        use param
        implicit none

        ! The reference triangle struct encodes a reference triangle with associated operators and collocation points.
        type RefTriangle
                double precision, allocatable :: M(:, :), Sr(:, :), Ss(:, :), Mi(:, :)    ! Mass, Stiffness and inverse Mass matrix
                double precision, allocatable :: Br(:, :), Bs(:, :)                       ! Boundary matrices in the MFSBP Language
                double precision, allocatable :: B(:,:,:)                                 ! The three side boundary matrices for side splitting framework
                double precision, allocatable :: Lo(:, :)                                 ! Projection operator onto the next lower order
                double precision, allocatable :: G(:, :)                                  ! Dissipation Operator
                double precision, allocatable :: w(:)                                     ! Positive volume quadrature
                double precision, allocatable :: cpts(:, :)                               ! Collocation points
                double precision, dimension(3, 2) :: x                                    ! Vertices of the reference triangle
                double precision, dimension(3, 2) :: n                                    ! normal vectors
                double precision, allocatable :: wsurf(:)                                 ! A quadrature for the edge
                integer :: Nbdt, Ncol, q                                                  ! Number of collocation points per edge, of collocation points, polynomial degree.
                integer, allocatable :: Bi(:, :)                                          ! Bi(i, j) is Node/Vector index of i-th node on j-th boundary
        end type

        ! A struct of arrays saving the Grid
        type Grid
                integer, allocatable :: verts(:, :)                     ! Vertice coordinates
                integer, allocatable :: adj(:, :)                       ! Adjacent triangles
                integer, allocatable :: adjO(:, :)                      ! The Orientation of adjacent triangles, i.e. if this triangle is the first, 
                                                                        !       second or third triangle for the other one
                integer, allocatable :: Mark(:)                         ! Markers for nodes

                integer, allocatable :: edges(:, :)                     ! Edges in the grid, (e, 1) is first triangle
                integer, allocatable :: edgenumb(:, :)                  ! Number of the edge in the triangle
                integer, allocatable :: edgeInd(:, :)                   ! A Nt x 3 array, the index into the edges array for a given side of a triangle 


                double precision, allocatable :: Pts(:, :)              ! Used points
                double precision, allocatable :: Cents(:, :)            ! Centers
                double precision, allocatable :: RotScals(:, :, :)      ! Rotation and scaling used to map Reference triangle onto triangle
                double precision, allocatable :: ToRef(:, :, :)         ! Its inverse
                double precision, allocatable :: Vol(:)                 ! Volume scaling
                double precision, allocatable :: ToRefDet(:)            ! Determinant of the mapping to the reference triangle
                double precision, allocatable :: Voli(:)                ! Inverted volume - removes division
                double precision, allocatable :: bl(:, :)               ! boundary length - removes two norm
                double precision, allocatable :: nn(:, :, :)            ! normed normal vectors
        end type


        contains

        ! Invert a 2x2 matrix in closed form
        function Invert2d(A) result (Ai)
                double precision, dimension(2, 2), intent(in) :: A
                double precision, dimension(2, 2) :: Ai
                double precision :: Adet

                Adet = A(1, 1)*A(2, 2) - A(1, 2)*A(2, 1)
                Ai(1, 1) = A(2, 2)/Adet
                Ai(1, 2) = -A(1, 2)/Adet
                Ai(2, 1) = -A(2, 1)/Adet
                Ai(2, 2) = A(1, 1)/Adet
        end function Invert2d

        ! Initializes the reference triangle, only sets up basic geometry
        subroutine RefInitBasic(RT)
                type(RefTriangle) :: RT
                RT%x(1, :) = [0.0D0, 0.0D0]
                RT%x(2, :) = [1.0D0, 0.0D0]
                RT%x(3, :) = [0.0D0, 1.0D0]
                RT%n(1, :) = [0.0D0, -1.0D0]
                RT%n(2, :) = [1.0D0, 1.0D0]
                RT%n(3, :) = [-1.0D0, 0.0D0]
        end subroutine RefInitBasic

        ! Initializes a complete reference triangle for a nodal DG method
        subroutine RefInitDG(RT)
                type(RefTriangle) :: RT
                integer :: q
                integer :: n, i, j
                double precision :: lam
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
                                rt%cpts(2*i, :) = 0.5*(RT%x(i, :) + RT%x(mod(i-1, 3)+2, :))
                        end do
                        allocate(RT%Bi(RT%Nbdt, 3))
                        do i=1, rt%Nbdt
                                rt%Bi(i, 1) = i
                                rt%Bi(i, 2) = i + rt%Nbdt-1
                                rt%Bi(i, 3) = mod(i+rt%Nbdt, n)+1

                        end do
                else if (OpDesign == "CG3") then        ! cubic CG triangles
                        q = 3                           ! polynomial degree three
                        n = 10                          ! 10 collocation nodes
                        rt%Nbdt = 4                     ! 4 nodes per edge, one in the interrior.
                        allocate(RT%cpts(n, 2))
                        do i=1, 3
                                do j=1, q
                                        lam = dble(j-1)/dble(q) 
                                        rt%cpts((rt%Nbdt-1)*(i-1) + j, :) = lam*RT%x(mod(i, 3)+1,:) + (1-lam)*RT%x(i, :)
                                end do
                        end do
                        rt%cpts(n, :) = [sum(rt%cpts(:, 1)), sum(rt%cpts(:, 2))] / 9.0D0
                        allocate(RT%Bi(rt%Nbdt, 3))
                        do i=1, rt%Nbdt
                                rt%Bi(i, 1) = i
                                rt%Bi(i, 2) = i + rt%Nbdt-1
                                rt%Bi(i, 3) = mod(i+rt%Nbdt+1, n-1)+1
                        end do
                else
                        write(*,*) "implement DG reference triangle", OpDesign, " !"
                end if
                rt%ncol = n
                rt%q = q
        end subroutine RefInitDG


        ! Calculates the coordinates of node nd in triangle tr
        ! using the grid gr and the reference triangle rt
        function GetNodeKord(gr, RT, tr, nd) result (x)
                type(Grid), intent(in) :: gr
                type(RefTriangle), intent(in) :: RT
                integer, intent(in) :: tr, nd
                double precision, dimension(2) :: x

                double precision, dimension(2) :: xloc

                ! Scale and rotate the point from the reference triangle to the local coordinates
                xloc = matmul(gr%RotScals(tr, :, :), RT%cpts(nd, :))
                ! sum with the 0 of the local coordinate system
                x = gr%Pts(gr%verts(tr, 1), :) + xloc 
        end function GetNodeKord


        ! Finds the Adjacent triangles and saves their connectivity in the grid structure.
        subroutine FindAdj(gr)
                implicit none
                type(Grid), intent(inout) :: gr
                ! aux
                integer :: i, j, k, l, Nt
                integer :: v1, v2, w1, w2
                Nt = size(gr%verts, 1)
                allocate(gr%Adj(Nt, 3))

                ! We iterate over all triangles. 
                do i=1, Nt
                        ! We iterate over the 3 edges.
 sideloop:              do j=1, 3
                                v1 = gr%verts(i, j)                     
                                v2 = gr%verts(i, modulo(j, 3) + 1)
                                ! We compare the markers of the node. If they are not zero it is an outer edge.
                                if (gr%Mark(v1) .NE. 0 .AND. gr%Mark(v2) .NE. 0) then   ! In this case it is an outer edge
                                        gr%Adj(i, j) = - gr%mark(v1)
                                        cycle
                                end if
                                ! If it is an inner edge we iterate over all other triangles an search for the same point
                                ! in another triangle.
                searchloop:     do k=1, Nt
                                        if (i .EQ. k) then              ! Skip, if it is the same triangle.
                                                cycle searchloop
                                        end if
                                        ! If two points of another edge match we save the found adjacent triangle
                                        do l=1, 3
                                                w1 = gr%verts(k, l)
                                                w2 = gr%verts(k, modulo(l, 3)+1)
                                                if (v1 .EQ. w1 .AND. v2 .EQ. w2) then
                                                        gr%Adj(i, j) = k
                                                        cycle sideloop
                                                end if
                                                if (v1 .EQ. w2 .AND. v2 .EQ. w1) then
                                                        gr%Adj(i, j) = k
                                                        cycle sideloop
                                                end if
                                        end do
                                end do searchloop
                                write(*,*) "No adj triangle found!", i, j, gr%Mark(v1), gr%Mark(v2)
                        end do  sideloop
                end do
        end subroutine FindAdj


        ! Sets the geometric fields of the grid based solely on the vertex coordinates
        subroutine FindGeo(gr)
                implicit none
                type(Grid), intent(inout) :: gr
                ! aux
                integer :: i, Nt
                double precision, dimension(2) :: v1, v2        ! two vectors spanning the triangle
                Nt = size(gr%verts, 1)
                ! First, we allocate memory
                allocate(gr%RotScals(Nt, 2, 2))
                allocate(gr%ToRef(Nt, 2, 2))
                allocate(gr%Vol(Nt))
                allocate(gr%ToRefDet(Nt))
                 
                do i=1, Nt
                        v1 = gr%Pts(gr%verts(i, 2), :) - gr%Pts(gr%verts(i, 1), :)
                        v2 = gr%Pts(gr%verts(i, 3), :) - gr%Pts(gr%verts(i, 1), :)

                        gr%RotScals(i, :, 1) = v1                               ! Maps the reference triangle to
                        gr%RotScals(i, :, 2) = v2                               ! the triangle in the grid
                        gr%ToRef(i, :, :) = Invert2D(gr%RotScals(i, :, :))      ! its inverse
                        gr%Vol(i) = abs(v1(1)*v2(2) - v1(2)*v2(1))/2.0          ! the volume of the triangle in the grid
                        gr%ToRefDet(i) = 0.5D0/gr%Vol(i)                ! the determinant of the mapping to the ref. triangle
                end do
        end subroutine FindGeo
        
        ! This subroutine is used to find the orientation of an adjacent triangle, i.e. if
        ! triangle i is adjacent to triangle nb on side j, on which side of triangle nb
        ! lies triangle i? 
        subroutine FindAdjO(gr)
                implicit none
                type(Grid), intent(inout) :: gr
                ! aux
                integer :: i, j, k, nb, Nt

                Nt = size(gr%verts, 1)
                allocate(gr%AdjO(Nt, 3))
                do i=1,Nt
                        do j=1,3
                                ! We iterate over every edge of every triangle. 
                                ! If the selected edge has a neighboring triangle
                                nb = Gr%Adj(i, j)
                                if (nb < 1) then
                                        cycle
                                end if
                                ! the corresponding edge index of the original triangle
                                ! is saved.
                                do k=1,3
                                        if (i .EQ. gr%Adj(nb, k)) then
                                                gr%AdjO(i, j) = k
                                        end if
                                end do
                        end do
                end do
        end subroutine FindAdjO

         
        ! Initializes the edge information in gr from the triangle information
        ! a list of all edges, i.e. the indices of the (up to two) triangles
        ! adjacent to every edge is kept in gr.
        subroutine EdgeListInit(gr)
                implicit none
                type(Grid), intent(inout) :: gr
                ! aux. variables
                integer, allocatable :: Edges(:, :), Edgenumb(:, :)
                integer :: i, j, k, Nt, Ne = 0, Noe = 0
                integer :: neib 
                Nt = size(gr%Verts, 1)
                ! We allocate buffers with a sufficient size.
                ! The memory in the struct is only allocated if we know
                ! the exact number of edges
                allocate(Edges(3*Nt, 2))
                allocate(Edgenumb(3*Nt, 2))
                do i=1, Nt
edgeloop:               do j = 1, 3
                                neib = gr%Adj(i, j)
                                ! Test if it is an outer triangle
                                if (neib < 1) then
                                        ! we mark outer edge by setting one of the triangle indices
                                        ! to a negative value
                                        Ne = Ne + 1
                                        Noe = Noe + 1
                                        Edges(Ne, 1) = i
                                        Edges(Ne, 2) = neib
                                        Edgenumb(Ne, 1) = j
                                        cycle edgeloop
                                end if
                                ! Test if the edge is already known
                                ! leave out self references
                                do k=1, Ne
                                        if ((Edges(k, 1) .EQ. neib) .AND. (Edges(k, 2) .EQ. i)) then
                                                cycle edgeloop
                                        end if
                                end do
                                ! We only arrive here if the edge is unknown.
                                ! We add it to the list.
                                Ne = Ne + 1
                                Edges(Ne, 1) = i
                                Edges(Ne, 2) = neib
                                Edgenumb(Ne, 1) = j
                                Edgenumb(Ne, 2)= gr%AdjO(i, j)
                        end do edgeloop

                end do
                write (*,*) "found a total of ", Ne, "edges and ", Noe, " outer edges"
                ! Allocate the memory in the structs, copy the array.
                allocate(gr%edges(Ne, 2))
                allocate(gr%edgenumb(Ne, 2))
                gr%edges(:, :) = Edges(1:Ne, :)
                gr%Edgenumb(:, :) = edgenumb(1:Ne, :)
                allocate(gr%edgeInd(Nt, 3))
                gr%edgeInd = 0
                do i=1, Ne
                        gr%edgeInd(gr%edges(i,1), gr%edgenumb(i, 1)) = i
                        if (gr%edges(i, 2) > 0) then     ! It is an internal edge
                               gr%edgeInd(gr%edges(i, 2), gr%edgenumb(i, 2)) = i
                        end if
                end do  
        end subroutine EdgeListInit

        ! Estimates the influence of the grid on the stiffnes of
        ! the spacial discretisation by finding the maximum of the length of a side
        ! divided by the triangle Volume
        function FindStiff(gr)
                type(Grid), intent(in) :: gr
                double precision :: FindStiff
                integer :: i
                double precision, dimension(2) :: p1, p2, p3

                FindStiff = 0.0D0
                do i=1, size(gr%Verts, 1)
                        p1 = gr%Pts(gr%Verts(i, 1), :)
                        p2 = gr%Pts(gr%Verts(i, 2), :)
                        p3 = gr%Pts(gr%Verts(i, 3), :)
                        FindStiff = max(FindStiff, norm2(p1-p2)/gr%Vol(i), norm2(p2-p3)/gr%Vol(i), norm2(p3-p1)/gr%Vol(i))
                end do
        end function FindStiff

        ! Ensures that the normal n points _outwards_
        ! by testing if the normal, when footed at fp, 
        ! categorizes the center point c as a point on the inside
        function Outerize(c, fp, n) result (on)
                implicit none
                double precision, dimension(2) :: c, fp, n
                double precision, dimension(2) :: on
                ! auxiallry
                if (dot_product(fp - c, n) .GE. 0.0D0) then
                        on = n
                else
                        on = -n
                end if
        end function Outerize

        ! errects a normal on v (rotation by 90 degrees)
        function errect(v)
                implicit none
                double precision, dimension(2) :: v
                double precision, dimension(2) :: errect
                errect(2) = v(1)
                errect(1) = -v(2)
        end function errect

        ! Initializes the normal information in the grid, i.e. edge
        ! lengths and corresponding _normalized_ normal vectors.
        subroutine GridPP(gr)
                implicit none
                type(Grid), intent(inout) :: gr
                ! auxillary vars
                double precision, dimension(2) :: cent, x1, x2, x3
                double precision, dimension(3, 2) :: n
                integer :: i,j, Ntr

                Ntr = size(gr%verts, 1)
                allocate(gr%Voli(Ntr))
                allocate(gr%bl(Ntr, 3))
                allocate(gr%nn(Ntr, 3, 2))
                do i=1, Ntr
                        gr%Voli(i) = 1.0D0/gr%vol(i)
                        x1 = gr%Pts(gr%Verts(i, 1), :)
                        x2 = gr%Pts(gr%Verts(i, 2), :)
                        x3 = gr%Pts(gr%Verts(i, 3), :)
                        cent = (x1 + x2 + x3) / 3.0D0
                        n(1, :) = outerize(cent, 0.5D0*(x1 + x2), errect(x2 - x1))
                        n(2, :) = outerize(cent, 0.5D0*(x2 + x3), errect(x3 - x2))
                        n(3, :) = outerize(cent, 0.5D0*(x3 + x1), errect(x1 - x3))
                        do j=1, 3
                                gr%bl(i, j) = norm2(n(j, :))
                                gr%nn(i, j, :) = n(j, :) / gr%bl(i, j)
                        end do
                end do
        end subroutine GridPP

        ! Prints some information in the grid struct gr for debugging purposes.
        subroutine GridDiag(gr)
                type(Grid), intent(in) :: gr
                integer :: i, j

                j = size(gr%verts, 1)

                888 FORMAT (A, I3, A, F6.3, F6.3)
                do i=1, size(gr%Pts, 1)
                        write (*, 888) "Point #", i, " has coordinates", gr%Pts(i, :)
                end do

                999 FORMAT (A, I5, A, I5, I5, I5, A, I5, I5, I5, A, F6.3)
                do i=1, j
                        write (*, 999) "Triangle #", i," has pts ", gr%verts(i, :), &
                                &", Adj. ", gr%adj(i, :), " and Volume ", gr%Vol(i)
                end do

        end subroutine GridDiag


        ! Saves the grid and the corresponding reference triangle
        ! in a HDF5 File. The triangles are exportes as an array
        ! of 3 integers per row denoting a triangle by its 3 points.
        ! the 3 integers refer to indices in an array of 2 doubles per
        ! row denoting a point each. The adjacent triangle array is
        ! exported as is.
        ! gr: the grid struct
        ! rt: the reference triangle struct
        ! filename: the filename for the HDF5 file.
        subroutine GridHDF5export(gr, rt, filename)
                ! Arguments
                type(Grid), intent(in) :: gr
                type(RefTriangle), intent(in) :: rt
                Character(LEN=*), intent(in) :: filename
                ! Parameters
                Character(LEN=*), parameter :: Tdsetname = "Triangles"
                Character(LEN=*), parameter :: Pdsetname = "Points"
                Character(LEN=*), parameter :: Adsetname = "Adjacents"
                Character(LEN=*), parameter :: CPdsetname = "ColPoints"
                Character(LEN=*), parameter :: groupname = "Grid"
                ! Aux. vars.
                Integer(HSIZE_T) :: file_id, Tdset_id, Pdset_id, Adset_id, CPdset_id
                Integer(HSIZE_T) :: Tdspace_id, Pdspace_id, Adspace_id, CPdspace_id
                Integer :: error

                Integer, parameter :: Trank = 2
                Integer, parameter :: Prank = 2
                Integer, parameter :: Arank = 2
                Integer, parameter :: CPrank = 2

                Integer(HSIZE_T), dimension(Trank) :: Tdims
                Integer(HSIZE_T), dimension(Prank) :: Pdims
                Integer(HSIZE_T), dimension(Arank) :: Adims
                Integer(HSIZE_T), dimension(CPrank) :: CPdims


                ! Set the dimensions
                Tdims = shape(gr%verts)
                Pdims = shape(gr%Pts)
                Adims = shape(gr%adj)
                CPdims = shape(rt%cpts)
                write (*,*) "Exporting ", Tdims, " triangle nodes"

                ! Open the HDF5 Fortran Interface
                call h5open_f(error)

                ! Create the File
                call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

                ! Create Dataspaces for Triangles, Points and Adjacencies
                call h5screate_simple_f(Trank, Tdims, Tdspace_id, error)
                call h5screate_simple_f(Prank, Pdims, Pdspace_id, error)
                call h5screate_simple_f(Arank, Adims, Adspace_id, error)
                call h5screate_simple_f(CPrank, CPdims, CPdspace_id, error)

                ! Create The Datasets
                call h5dcreate_f(file_id, Tdsetname, H5T_NATIVE_INTEGER, Tdspace_id, Tdset_id, error)
                call h5dcreate_f(file_id, Pdsetname, H5T_NATIVE_DOUBLE, Pdspace_id, Pdset_id, error)
                call h5dcreate_f(file_id, Adsetname, H5T_NATIVE_INTEGER, Adspace_id, Adset_id, error)
                call h5dcreate_f(file_id, CPdsetname, H5T_NATIVE_DOUBLE, CPdspace_id, CPdset_id, error)

                ! Write the Datasets
                CALL h5dwrite_f(Tdset_id, H5T_NATIVE_INTEGER, gr%verts, Tdims, error)
                CALL h5dwrite_f(Pdset_id, H5T_NATIVE_DOUBLE, gr%Pts, Pdims, error)
                CALL h5dwrite_f(Adset_id, H5T_NATIVE_INTEGER, gr%adj, Adims, error)
                CALL h5dwrite_f(CPdset_id, H5T_NATIVE_DOUBLE, rt%cpts, CPdims, error)

                ! Close everything down
                CALL h5dclose_f(Tdset_id, error)
                CALL h5dclose_f(Pdset_id, error)
                CALL h5dclose_f(Adset_id, error)
                CALL h5dclose_f(CPdset_id, error)

                CALL h5sclose_f(Tdspace_id, error)
                CALL h5sclose_f(Pdspace_id, error)
                CALL h5sclose_f(Adspace_id, error)
                CALL h5sclose_f(CPdspace_id, error)

                CALL h5fclose_f(file_id, error)
                CALL h5close_f(error)
        end subroutine GridHDF5export

        ! This subroutine saves the tensor of conserved variables in an HDF5 file
        ! u: 3-Tensor, first dim for the triangle index, second dim for collocation pt.
        !       third dim for the number of the conserved variable.
        ! t: time
        subroutine SolHDF5export(u, t, filename)
                ! Arguments
                double precision, intent(in) :: u(:, :, :)
                double precision :: t
                Character(LEN=*), intent(in) :: filename
                ! Parameters
                Character(LEN=*), parameter :: udsetname = "u"
                Character(LEN=*), parameter :: tdsetname = "t"
                Character(LEN=*), parameter :: groupname = "Solution"
                ! Aux. variables.
                Integer(HSIZE_T) :: file_id, udset_id, tdset_id
                Integer(HSIZE_T) :: udspace_id, tdspace_id
                Integer :: error

                Integer, parameter :: urank = 3         ! The triangle number, the node number, the conserved variables
                Integer, parameter :: trank = 1

                Integer(HSIZE_T), dimension(urank) :: udims
                Integer(HSIZE_T), dimension(trank) :: tdims


                ! Set the dimensions
                udims = shape(u)
                tdims(1) = 1

                ! Open the HDF5 Fortran Interface
                call h5open_f(error)

                ! Create the File
                call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

                ! Create Dataspaces for Triangles, Points and Adjacencies
                call h5screate_simple_f(urank, udims, udspace_id, error)
                call h5screate_simple_f(trank, tdims, tdspace_id, error)

                ! Create The Datasets
                call h5dcreate_f(file_id, udsetname, H5T_NATIVE_DOUBLE, udspace_id, udset_id, error)
                call h5dcreate_f(file_id, tdsetname, H5T_NATIVE_DOUBLE, tdspace_id, tdset_id, error)

                ! Write the Datasets
                CALL h5dwrite_f(udset_id, H5T_NATIVE_DOUBLE, u, udims, error)
                CALL h5dwrite_f(tdset_id, H5T_NATIVE_DOUBLE, [t], tdims, error)

                ! Close everything down
                CALL h5dclose_f(udset_id, error)
                CALL h5dclose_f(tdset_id, error)
                CALL h5sclose_f(udspace_id, error)
                CALL h5sclose_f(tdspace_id, error)
                CALL h5fclose_f(file_id, error)
                CALL h5close_f(error)
        end subroutine SolHDF5export

        ! Imports a unstructured triangle grid in the format of Skewchuks Triangle
        ! grid generator. 
        ! gr: resulting grid struct
        ! nfname: Filename of the nodes/points
        ! efname: Filename of the .ele file holding the triangles.
        subroutine GridTriangleImp(gr, nfname, efname)
                implicit none
                type(Grid), intent(inout) :: gr
                Character(LEN=*), intent(in) :: nfname, efname
                ! auxillaries
                integer, parameter :: nfunit = 10, efunit = 11
                integer :: Nnodes, Ntriangles, dims, nattr, nmark
                integer :: i, j

                open (nfunit, file=nfname)
                open (efunit, file=efname)

                read (nfunit, *) Nnodes, dims, nattr, nmark
                if (dims .NE. 2) then
                        write (*,*) "wrong number of node dimension"
                        stop
                end if
                if (nattr .NE. 0) then
                        write (*,*) "wrong number of attributes"
                        stop
                end if
                if (nmark .NE. 1) then
                        write (*,*) "wrong number of markers"
                        stop
                end if

                allocate(gr%Pts(Nnodes, 2))
                allocate(gr%Mark(Nnodes))

                do i=1, Nnodes
                        read (nfunit, *) j, gr%Pts(i, 1), gr%Pts(i, 2), gr%Mark(i)
                end do

                read (efunit, *) Ntriangles, dims, nattr
                allocate(gr%verts(Ntriangles, 3))
                do i=1, Ntriangles
                        read (efunit, *) j, gr%verts(i, 1), gr%verts(i, 2), gr%verts(i, 3)
                end do
                close (nfunit)
                close (efunit)
                write (*,*) "Loaded ", Nnodes, " Nodes and ", Ntriangles, " triangles from files ", trim(nfname), " ", trim(efname)
    
                ! Find the adjoint triangles
                call FindAdj(gr)
                ! Initialize the geometry fields
                call FindGeo(gr)
                ! Initialize the Adjoint Orientation fields
                call FindAdjO(gr)
                ! Initialize the list of edges
                call EdgeListInit(gr)
                call gridpp(gr) 
        end subroutine GridTriangleImp

end module triangles

