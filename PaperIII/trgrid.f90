! The module triangles contains the routines and data structures to handle the 
! unstructured conforming triangle grids of the solver.
! Especially the loading of
! grids in the output format of Jonathan Skewchuks Triangle library. 

module trianggrid
        use HDF5
        use param
        use trianghp
        implicit none


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


                real(np), allocatable :: Pts(:, :)              ! Used points
                real(np), allocatable :: Cents(:, :)            ! Centers
                real(np), allocatable :: RotScals(:, :, :)      ! Rotation and scaling used to map Reference triangle onto triangle
                real(np), allocatable :: ToRef(:, :, :)         ! Its inverse
                real(np), allocatable :: Vol(:)                 ! Volume scaling
                real(np), allocatable :: ToRefDet(:)            ! Determinant of the mapping to the reference triangle
                real(np), allocatable :: Voli(:)                ! Inverted volume - removes division
                real(np), allocatable :: bl(:, :)               ! boundary length - removes two norm
                real(np), allocatable :: nn(:, :, :)            ! normed normal vectors
        end type


        contains

        ! Invert a 2x2 matrix in closed form
        function Invert2d(A) result (Ai)
                real(np), dimension(2, 2), intent(in) :: A
                real(np), dimension(2, 2) :: Ai
                real(np) :: Adet

                Adet = A(1, 1)*A(2, 2) - A(1, 2)*A(2, 1)
                Ai(1, 1) = A(2, 2)/Adet
                Ai(1, 2) = -A(1, 2)/Adet
                Ai(2, 1) = -A(2, 1)/Adet
                Ai(2, 2) = A(1, 1)/Adet
        end function Invert2d


        ! Calculates the coordinates of node nd in triangle tr
        ! using the grid gr and the reference triangle rt
        function GetNodeKord(gr, RT, tr, nd) result (x)
                type(Grid), intent(in) :: gr
                type(nprt), intent(in) :: RT
                integer, intent(in) :: tr, nd
                real(np), dimension(2) :: x

                real(np), dimension(2) :: xloc

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
                real(np), dimension(2) :: v1, v2        ! two vectors spanning the triangle
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
                        gr%ToRefDet(i) = 0.5_np/gr%Vol(i)                ! the determinant of the mapping to the ref. triangle
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
                real(np), dimension(2) :: p1, p2, p3

                FindStiff = 0.0_np
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
                real(np), dimension(2) :: c, fp, n
                real(np), dimension(2) :: on
                ! auxiallry
                if (dot_product(fp - c, n) .GE. 0.0_np) then
                        on = n
                else
                        on = -n
                end if
        end function Outerize

        ! errects a normal on v (rotation by 90 degrees)
        function errect(v)
                implicit none
                real(np), dimension(2) :: v
                real(np), dimension(2) :: errect
                errect(2) = v(1)
                errect(1) = -v(2)
        end function errect

        ! Initializes the normal information in the grid, i.e. edge
        ! lengths and corresponding _normalized_ normal vectors.
        subroutine GridPP(gr)
                implicit none
                type(Grid), intent(inout) :: gr
                ! auxillary vars
                real(np), dimension(2) :: cent, x1, x2, x3
                real(np), dimension(3, 2) :: n
                integer :: i,j, Ntr

                Ntr = size(gr%verts, 1)
                allocate(gr%Voli(Ntr))
                allocate(gr%bl(Ntr, 3))
                allocate(gr%nn(Ntr, 3, 2))
                do i=1, Ntr
                        gr%Voli(i) = 0.0_np/gr%vol(i)
                        x1 = gr%Pts(gr%Verts(i, 1), :)
                        x2 = gr%Pts(gr%Verts(i, 2), :)
                        x3 = gr%Pts(gr%Verts(i, 3), :)
                        cent = (x1 + x2 + x3) / 3.0_np
                        n(1, :) = outerize(cent, 0.5_np*(x1 + x2), errect(x2 - x1))
                        n(2, :) = outerize(cent, 0.5_np*(x2 + x3), errect(x3 - x2))
                        n(3, :) = outerize(cent, 0.5_np*(x3 + x1), errect(x1 - x3))
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
                type(nprt), intent(in) :: rt
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

end module trianggrid

