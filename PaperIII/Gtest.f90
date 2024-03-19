! This program tests the construction of a positive, conservative filter
! in two space dimensions on a triangle, as defined in the publication.

program Gtest
        use ops
        use Triangles
        
        type(RefTriangle) :: rt
        double precision, dimension(10, 10) :: G
        integer :: i
        call DGinit(rt)
        call FindDispG(rt, G)
        
        write(*,*) "G Matrix is"
        do i=1,10
                write(*,"(10f10.3)") G(i, :)
        end do 
        
        write(*,*) " Now checking row sums"
        do i=1, 10
                write (*,*) "row ", i, " Has sum 1 + ", sum(G(i, :))-1.0D0
        end do

        write(*, *) " Now checking column integral"
        do i=1, 10
                write (*,*) "column ", i, " Has integral error", sum(matmul(rt%M, G(:, i))) - sum(rt%M(:, i))
        end do
end program Gtest
