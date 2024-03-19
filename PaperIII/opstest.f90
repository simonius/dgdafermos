! This program serves as a unit test to guarantee that the DG operator
! construction routines work as intended

program opstest
        use ops
        use triangles

        type(RefTriangle) :: rt
        double precision, allocatable :: O(:,:), V(:,:), Vx(:,:), Vy(:,:), Dx(:,:), Dy(:,:), Gx(:,:), Gy(:,:)

        integer :: i, n

        call DGinit(rt)
        allocate(O(rt%ncol, rt%ncol))
        allocate(V(rt%ncol, rt%ncol))
        allocate(Vx(rt%ncol, rt%ncol))
        allocate(Vy(rt%ncol, rt%ncol))
        allocate(Dx(rt%ncol, rt%ncol))
        allocate(Dy(rt%ncol, rt%ncol))
        allocate(Gx(rt%ncol, rt%ncol))
        allocate(Gy(rt%ncol, rt%ncol))
        
        n = rt%ncol

        write(*,*) "The nodal Mass matrix is"
        do i=1,n
                write (*,"(10f8.3)") rt%M(i,:)
        end do

        write(*,*) "The inverse nodal Mass matrix is"
        do i=1,n
                write (*,"(10f8.3)") rt%Mi(i,:)
        end do

        write(*,*) "The nodal Stiffness matrix in x is"
        do i=1,n
                write (*,"(10f8.3)") rt%Sr(i,:)
        end do

        write(*,*) "The nodal Stiffness matrix in y is"
        do i=1,n
                write (*,"(10f8.3)") rt%Ss(i,:)
        end do

        write(*,*) "The nodal boundary matrix in x is"
        do i=1,n
                write (*,"(10f8.3)") rt%Br(i,:)
        end do

        write(*,*) "The nodal boundary matrix in y is"
        do i=1,n
                write (*,"(10f8.3)") rt%Bs(i,:)
        end do

        write (*,*) "Testing SBP property in x and y"
        do i=1,n
                write (*,"(10f8.3)") (rt%Sr(i,:) + rt%Sr(:, i) - rt%Br(i, :))
                write (*, "(10f8.3)") (rt%Ss(i,:) + rt%Ss(:, i) - rt%Bs(i, :))
        end do

        

        Dx = matmul(rt%Mi, transpose(rt%Sr))
        Dy = matmul(rt%Mi, transpose(rt%Ss))
        
        ! Orthonormal basis 
        call DOBM(O, rt%ncol, rt%cpts)
        call DOBV(V, O, rt%ncol, rt%cpts)
        call DOBVx(Vx, O, rt%ncol, rt%cpts)
        call DOBVy(Vy, O, rt%ncol, rt%cpts)
        Gx = matmul(Dx, V)
        Gy = matmul(Dy, V)
        
        write (*,*) "Testing Derivative exactness in x"
        do i=1,n
                write (*,"(10f10.3)") (Vx(i, :) - Gx(i, :))
                write (*,"(10f8.3)") (Vy(i, :) - Gy(i, :))
        end do




end program opstest
