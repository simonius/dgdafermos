! This module contains a quick and dirty implementation of linear algebra routines
! in quad precision to remove LAPACK and BLAS from the list of dependencies.
! These versions are _not_ optimized for speed as they only see limited use, for example to
! invert mass matrices.

module qpla
        use param
        implicit none

        contains

        ! carries out a QR factorization using Householder reflections
        subroutine QRH(Q, R, A)
                implicit none
                real(hp), intent(in) :: A(:, :)
                real(hp), intent(out) :: Q(:, :), R(:, :)
                ! auxillary
                real(hp), dimension(size(A, 1), 1) :: u
                real(hp) :: alp
                integer :: i, n, k
                ! copy the initial matrix into R
                R = A
                n = size(A, 1)
                k = size(A, 2)
                ! initialize Q as identity
                Q = 0.0_hp
                do i=1, n
                        Q(i, i) = 1.0_hp
                end do
                ! Reflect until R is upper triangular, concatenate
                ! the inverse reflections in Q.
                do i=1, n-1
                        alp = -sign(norm2(R(i:n, i)), R(i, i))
                        u(i:n, 1) = R(i:n, i)
                        u(i, 1) = u(i, 1) - alp
                        u(i:n, 1) = u(i:n, 1) / norm2(u(i:n, 1))
                        R(i:n, i:k) = R(i:n, i:k) - 2.0_hp*matmul(u(i:n, :), matmul(transpose(u(i:n, :)), R(i:n, i:k)))
                        Q(:, i:k) = Q(:, i:k) - 2.0_hp*matmul(matmul(Q(:, i:k), u(i:n, :)), transpose(u(i:n, :)))
                end do
        end subroutine QRH

        ! carries out a backsubstitution for an upper triangular R
        function BS(R, b) result (x)
                real(hp), intent(in) :: R(:, :)
                real(hp), intent(in) :: b(:)
                real(hp), dimension(size(R, 2)) :: x
                !  axillary
                integer :: i, n
                n = size(R, 1)
                do i=n, 1, -1
                        x(i) = (b(i) - dot_product(R(i, i+1:n), x(i+1:n))) / R(i, i)
                end do
        end function BS

        ! inverts A using a QR decomposition
        function invQR(A) result (Ai)
                real(hp), intent(in) :: A(:, :)
                real(hp) :: Ai(size(A, 1), size(A, 2))
                ! auxillary
                real(hp), dimension(size(A, 1), size(A, 2)) :: Q, R
                integer :: i
                
                call QRH(Q, R, A)
                do i=1, size(Ai, 1)
                        Ai(:, i) = BS(R, Q(i, :))
                end do
        end function invQR

        function MatExpSeries(Q) result(A)
                real(hp), dimension(:, :), intent(in) :: Q
                real(hp), dimension(size(Q, 1), size(Q, 2)) :: A
                ! Aux
                real(hp) :: B(size(Q, 1), size(Q, 2))
                integer :: i
                B = Q
                A = 0.0_hp
                do i=1, size(Q, 1)
                        A(i, i) = 1.0_hp
                end do
                do i=1, 30
                        A = A + B
                        B = matmul(B,Q)/real(i+1, hp)
                        !write(*,"(10f10.5)") B(i, :)
                end do
        end function MatExpSeries

        ! We calculate the matrix exponential using the division and squaring approach.
        function MatExp(Q) result (A)
                real(hp), dimension(:, :), intent(in) :: Q
                real(hp), dimension(size(Q, 1), size(Q, 2)) :: A
                ! Aux
                real(hp) :: matnorm
                real(hp) :: S(size(Q, 1), size(Q, 2)), R(size(Q, 1), size(Q, 2))
                integer :: i, divexp
                matnorm = norm2(Q)
                divexp = ceiling(log(matnorm) / log(2.0_hp))
                S = Q / 2.0_hp**divexp
                R = MatExpSeries(S)
                do i=1, divexp
                        R = matmul(R, R)
                end do
                A = R
        end function MatExp

end module qpla
