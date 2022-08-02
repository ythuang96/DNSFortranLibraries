module wallnormalgrid
    use types, only: sp, dp
    implicit none
    private
    include 'parameters'

    public :: prepare_grid, d_dy_complexvector

    ! derivative matrices
    real(kind=cp) :: dt11(5,my),dt12(5,my)


contains
    ! subroutine prepare_grid
    ! Computes the variables for the y grid
    ! The subroutine will compute:
    ! y               : the y grid
    ! dt11, dt12, fmap: matrices for the first derivative in y
    ! trp             : trapezoid integration coefficients in y
    ! All of these will be set as global variables that can be accessed anywhere
    !
    ! call prepare_grid once when initializing the code
    !
    ! To use dt11 and dt12, put "use matrices"
    ! To use y, fmap, trp, put:
    ! "real(kind=cp) :: fmap, y, trp"
    ! "common /ygrid/ y(my), fmap(my), trp(my)"
    ! "save   /ygrid/"
    subroutine prepare_grid
        real(kind=cp) :: fmap, y, trp
        common /ygrid/ y(my), fmap(my), trp(my)
        save   /ygrid/

        real(kind=cp), dimension(my) :: b
        integer :: i, j


        ! Compute y grid
        y = y_grid( gamma )

        ! Compute first derivate matrix dt11 and dt12
        call first_derivative_matrix(dt11,dt12)
        ! Process dt11 for the 5 diagonal solver
        call bandec5(dt11,my)

        ! Compute fmap
        ! b = D12 * y
        b = dt12_times_realvector( y )
        ! b = D11\b = D11\(D12*y) = derivative of y w.r.t uniform grid
        call banbks5(dt11,my,b)
        ! fmap = 1/derivative of y
        fmap = 1.0_cp / b

        ! Compute trapezoid integration coefficient
        trp(1) = 0.5_cp*( y(2)-y(1) )
        do j=2,my-1
            trp(j) = 0.5_cp*( y(j+1)-y(j-1) )
        enddo
        trp(my)= 0.5_cp*( y(my)-y(my-1) )
    end subroutine prepare_grid


    ! function dudy = d_dy_complexvector( u )
    ! Computes the first y derivative for a complex vector
    !
    ! Arguments:
    !   u   : [complex double/single, size my, Input]
    !         Input vector
    ! Output:
    !   dudy: [complex double/single, size my, Output]
    !         first y derivative of input vector
    function d_dy_complexvector( u ) result( dudy )
        complex(kind=cp), dimension(my), intent(in) :: u
        complex(kind=cp), dimension(my) :: dudy

        real(kind=cp) :: fmap, y, trp
        common /ygrid/ y(my), fmap(my), trp(my)
        save   /ygrid/

        real(kind=cp), dimension(my) :: real_buffer, imag_buffer


        real_buffer = dt12_times_realvector( real( u, cp ) )
        call banbks5( dt11,my,real_buffer )
        real_buffer = real_buffer * fmap

        imag_buffer = dt12_times_realvector( aimag( u ) )
        call banbks5( dt11,my,imag_buffer )
        imag_buffer = imag_buffer * fmap

        dudy = CMPLX( real_buffer, imag_buffer, cp)
    end function d_dy_complexvector


    ! function y = y_grid( g )
    ! Generates a stretched sin grid
    !
    ! Arguments:
    !   g: [double, Input]
    !      grid stretching parameter
    ! Output:
    !   y: [double, size my, Output]
    !      y coordinates
    function y_grid( g ) result( y )
        real(kind=cp), intent(in) :: g
        real(kind=cp), dimension(my) :: y

        integer :: j
        real(kind=cp) :: pi

        pi = acos(-1.0_cp)

        do j = 1,my
            y(j) = sin( g*pi/2.0_cp * ( 2.0_cp*real(j-1, cp)/real(my-1, cp)-1.0_cp ) ) / &
                   sin( g*pi/2.0_cp )
        enddo
    end function y_grid


    ! subroutine first_derivative_matrix( d11,d12 )
    ! Generates d11 and d12 matrices for the first derivative in y
    ! d11, d12 are the five diagonals for matrix D11, D12, the first row of d11
    ! is the lowest diagonal of D11.
    ! D11 * du/dy = D12 * u
    !
    ! Arguments:
    !   d11: [double, size (5,my), Output]
    !        the five diagonals for matrix D11
    !   d12: [double, size (5,my), Output]
    !        the five diagonals for matrix D12
    subroutine first_derivative_matrix( d11,d12 )
        real(kind=cp), intent(out), dimension(5, my) :: d11, d12

        real(kind=cp) :: b(9), h
        integer :: j


        h = 2.0_cp/real(my-1, cp)

        b(1) = real( 40941015625d-12, cp)
        b(2) = real(  4905703125d-10, cp)
        b(3) = ( - 16.0_cp*b(1) + 2.0_cp*b(2) + 4.0_cp ) / (  6.0_cp*h )
        b(4) = (   22.0_cp*b(1) + 4.0_cp*b(2) - 1.0_cp ) / ( 12.0_cp*h )

        do j=3,my-2
            d11(1:5,j) = (/ b(1), b(2), 1.0_cp, b(2), b(1) /)

            d12(1:5,j) = (/ -b(4), -b(3), 0.0_cp, b(3), b(4) /)
        enddo

        d11(1:5,1   ) = (/ 0.0_cp,  0.0_cp, 1.0_cp,  2.0_cp, 0.0_cp /)
        d11(1:5,2   ) = (/ 0.0_cp, 0.25_cp, 1.0_cp, 0.25_cp, 0.0_cp /)
        d11(1:5,my-1) = (/ 0.0_cp, 0.25_cp, 1.0_cp, 0.25_cp, 0.0_cp /)
        d11(1:5,my  ) = (/ 0.0_cp,  2.0_cp, 1.0_cp,  0.0_cp, 0.0_cp /)

        d12(1:5,1   ) = (/ 0.0_cp, 0.0_cp, -15.0_cp/(6.0_cp*h),  2.0_cp/h,  1.0_cp/(2.0_cp*h) /)
        d12(1:5,2   ) = (/ 0.0_cp, -3.0_cp/(4.0_cp*h), 0.0_cp, 3.0_cp/(4.0_cp*h), 0.0_cp /)
        d12(1:5,my-1) = (/ 0.0_cp, -3.0_cp/(4.0_cp*h), 0.0_cp, 3.0_cp/(4.0_cp*h), 0.0_cp /)
        d12(1:5,my  ) = (/ -1.0_cp/(2.0_cp*h), -2.0_cp/h,  15.0_cp/(6.0_cp*h), 0.0_cp, 0.0_cp /)
    end subroutine first_derivative_matrix


    ! function vector_out = dt12_times_realvector( vector_in )
    ! Applies dt12 to a real vector, where d12 is the five diagonals for matrix
    ! D12, and the first row of dt12 is the lowest diagonal of D12.
    !
    ! Arguments:
    !   vector_in : [double, size my, Input]
    !               input vector
    ! Output:
    !   vector_out: [double, size my, Output]
    !               output vector = D12 * vector_in
    function dt12_times_realvector( vector_in ) result( vector_out )
        real(kind=cp), dimension(myf), intent(in) :: vector_in
        real(kind=cp), dimension(myf) :: vector_out

        integer :: i, j


        vector_out(1) = dt12(3,1)*vector_in(1) &
                      + dt12(4,1)*vector_in(2) &
                      + dt12(5,1)*vector_in(3)
        vector_out(2) = dt12(2,2)*vector_in(1) &
                      + dt12(3,2)*vector_in(2) &
                      + dt12(4,2)*vector_in(3) &
                      + dt12(5,2)*vector_in(4)
        do j = 3,my-2
            vector_out(j) = 0.0_cp
            do i = 1,5
                vector_out(j) = vector_out(j) + dt12(i,j)*vector_in(i+j-3)
            enddo
        enddo
        vector_out(my-1) = dt12(1,my-1)*vector_in(my-3) &
                         + dt12(2,my-1)*vector_in(my-2) &
                         + dt12(3,my-1)*vector_in(my-1) &
                         + dt12(4,my-1)*vector_in(my)
        vector_out(my)   = dt12(1,my  )*vector_in(my-2) &
                         + dt12(2,my  )*vector_in(my-1) &
                         + dt12(3,my  )*vector_in(my)
    end function dt12_times_realvector


    ! subroutine banbks5(a,n,b)
    ! subroutine bandec5(a,n)
    !
    ! Solve a 5 diagonal system b = A\b, where A is a 5 diagonal n X n matrix
    ! with diagonals given by the rows of a (first row is the lowest diagonal)
    ! and b is a size n column vector
    !
    ! first run bandec5 to update a (modify a in preperation for banbks5)
    ! then run banbks5 with the updated a to update b to the solution
    SUBROUTINE banbks5(a,n,b)
        INTEGER, intent(in) :: n
        REAL(kind=cp), intent(in) :: a(5,n)
        REAL(kind=cp), intent(inout) :: b(n)

        INTEGER :: i,k

        do k=1,n-2
            b(k+1) = b(k+1)-a(4,k)*b(k)
            b(k+2) = b(k+2)-a(5,k)*b(k)
        enddo
        b(n) = b(n)- a(4,n-1)*b(n-1)

        ! back substitution
        b(n) = b(n)*a(1,n)
        b(n-1) = (b(n-1)-a(2,n-1)*b(n))*a(1,n-1)
        do i=n-2,1,-1
            b(i) = (b(i)-a(2,i)*b(1+i)-a(3,i)*b(2+i))*a(1,i)
        enddo
    END SUBROUTINE banbks5

    SUBROUTINE bandec5(a,n)
        INTEGER, intent(in) :: n
        REAL(kind=cp), intent(inout) :: a(5,n)

        INTEGER :: j,k

        do j = 1,3
            a(j,1)=a(j+2,1)
        enddo
        do j = 1,4
            a(j,2)=a(j+1,2)
        enddo

        do k = 1,n-2
            a(1,k)   = 1.0_cp/a(1,k)

            a(4,k)   = a(1,k+1)*a(1,k)

            a(1,k+1) = a(2,k+1)-a(4,k)*a(2,k)
            a(2,k+1) = a(3,k+1)-a(4,k)*a(3,k)
            a(3,k+1) = a(4,k+1)

            a(5,k)   = a(1,k+2)*a(1,k)

            a(1,k+2) = a(2,k+2)-a(5,k)*a(2,k)
            a(2,k+2) = a(3,k+2)-a(5,k)*a(3,k)
            a(3,k+2) = a(4,k+2)
            a(4,k+2) = a(5,k+2)
        enddo

        a(1,n-1) = 1.0_cp/a(1,n-1)

        a(4,n-1)=a(1,n)*a(1,n-1)

        a(1,n) = a(2,n)-a(4,n-1)*a(2,n-1)
        a(2,n) = a(3,n)-a(4,n-1)*a(3,n-1)
        a(3,n) = a(4,n)

        a(1,n)=1.0_cp/a(1,n)
    END SUBROUTINE bandec5





end module wallnormalgrid