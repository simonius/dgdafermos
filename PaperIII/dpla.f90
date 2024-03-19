! This module contains a quick and dirty implementation of linear algebra routines
! in double precision to remove LAPACK and BLAS from the list of dependencies.
! These versions are _not_ optimized for speed as they only see limited use, for example to
! invert mass matrices.

module dpla
        implicit none
        contains

        ! carries out a QR factorization using Householder reflections
        subroutine QRH(Q, R, A)
                implicit none
                double precision, intent(in) :: A(:, :)
                double precision, intent(out) :: Q(:, :), R(:, :)
                ! auxillary
                double precision, dimension(size(A, 1), 1) :: u
                double precision :: alp
                integer :: i, n, k
                ! copy the initial matrix into R
                R = A
                n = size(A, 1)
                k = size(A, 2)
                ! initialize Q as identity
                Q = 0.0D0
                do i=1, n
                        Q(i, i) = 1.0D0
                end do
                ! Reflect until R is upper triangular, concatenate
                ! the inverse reflections in Q.
                do i=1, n-1
                        alp = -sign(norm2(R(i:n, i)), R(i, i))
                        u(i:n, 1) = R(i:n, i)
                        u(i, 1) = u(i, 1) - alp
                        u(i:n, 1) = u(i:n, 1) / norm2(u(i:n, 1))
                        R(i:n, i:k) = R(i:n, i:k) - 2.0D0*matmul(u(i:n, :), matmul(transpose(u(i:n, :)), R(i:n, i:k)))
                        Q(:, i:k) = Q(:, i:k) - 2.0D0*matmul(matmul(Q(:, i:k), u(i:n, :)), transpose(u(i:n, :)))
                end do
        end subroutine QRH

        ! carries out a backsubstitution for an upper triangular R
        function BS(R, b) result (x)
                double precision, intent(in) :: R(:, :)
                double precision, intent(in) :: b(:)
                double precision, dimension(size(R, 2)) :: x
                !  axillary
                integer :: i, n
                n = size(R, 1)
                do i=n, 1, -1
                        x(i) = (b(i) - dot_product(R(i, i+1:n), x(i+1:n))) / R(i, i)
                end do
        end function BS

        ! inverts A using a QR decomposition
        function invQR(A) result (Ai)
                double precision, intent(in) :: A(:, :)
                double precision :: Ai(size(A, 1), size(A, 2))
                ! auxillary
                double precision, dimension(size(A, 1), size(A, 2)) :: Q, R
                integer :: i
                
                call QRH(Q, R, A)
                do i=1, size(Ai, 1)
                        Ai(:, i) = BS(R, Q(i, :))
                end do
        end function invQR
end module dpla
