module matrixinv
    use types, only: sp, dp
    implicit none
    private

    public :: MatrixInverse

contains
    ! Returns the inverse of a matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.
    !
    ! Arguments
    !   A: [double compelx NXN matrix, In/Out] orginal matrix on input, and the inverse of it on output
    subroutine MatrixInverse( A )
        complex(kind=dp), dimension(:,:), intent(inout) :: A

        complex(kind=dp), dimension(size(A,1)) :: work  ! work array for LAPACK
        integer, dimension(size(A,1)) :: ipiv   ! pivot indices
        integer :: n, info

        ! External procedures defined in LAPACK
        external ZGETRF
        external ZGETRI

        ! dimension of A
        n = size(A,1)

        ! ZGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call ZGETRF(n, n, A, n, ipiv, info)

        if (info /= 0) then
            stop 'Matrix is numerically singular!'
        end if

        ! ZGETRI computes the inverse of a matrix using the LU factorization
        ! computed by ZGETRF.
        call ZGETRI(n, A, n, ipiv, work, n, info)

        if (info /= 0) then
            stop 'Matrix inversion failed!'
        end if
    end subroutine MatrixInverse

end module matrixinv