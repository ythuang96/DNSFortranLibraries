module wavenumbers
    use types, only: sp, dp
    implicit none
    private
#   include "parameters"

    real(kind=cp) :: kx_full, kx_pos, kz_full, kz_pos
    common /wave/ kx_full(nkx_full), kx_pos(nkx_pos), &
                  kz_full(nkz_full), kz_pos(nkz_pos)
    save /wave/

    public :: generate_wavenumbers, kx_derivative, kz_derivative


contains
    ! subroutine generate_wavenumbers
    ! This function generates the wavenumbers
    ! kx_full, kx_pos, kz_full, kz_pos
    ! full indicates it includes both the positive and negative wavenumbers
    ! pos  indicates it only includes the 0 and positive wavenumbers
    !
    ! These four variables are set as global variables
    subroutine generate_wavenumbers
        integer :: ii


        ! kx >= 0 part
        DO ii = 1, nkx_pos
            kx_full(ii) = real(ii-1, cp) * dkx
        ENDDO
        ! kx <  0 part
        DO ii = nkx_pos+1, nkx_full
            kx_full(ii) = real(ii-1-nkx_full, cp) * dkx
        ENDDO

        ! kz >= 0 part
        DO ii = 1, nkz_pos
            kz_full(ii) = real(ii-1, cp) * dkz
        ENDDO
        ! kz <  0 part
        DO ii = nkz_pos+1, nkz_full
            kz_full(ii) = real(ii-1-nkz_full, cp) * dkz
        ENDDO

        kx_pos(:) = kx_full(1:nkx_pos)
        kz_pos(:) = kz_full(1:nkz_pos)

    end subroutine generate_wavenumbers


    ! function matrix_out = kx_derivative( matrix_in )
    ! This function multiplies the matrix_in by i kx, which is differenetiating
    ! in x in the Fourier domain
    !
    ! Arguments:
    !   matrix_in : [double/single complex, size (mxf,mzf), Input]
    !               velocity field in the Fourier domain at a single y plane
    ! Output:
    !   matrix_out: [double/single complex, size (mxf,mzf), Output]
    !               x derivative of the velocity field in the Fourier domain
    !               at a single y plane
    function kx_derivative( matrix_in ) result( matrix_out )
        complex(kind=cp), dimension(mxf,mzf), intent(in) :: matrix_in
        complex(kind=cp), dimension(mxf,mzf) :: matrix_out

        complex(kind=cp) :: complex_i
        integer :: ii


        complex_i = ( 0.0_cp, 1.0_cp )

        DO ii = 1, mzf
            matrix_out(:, ii) = matrix_in(:, ii) * kx_pos * complex_i
        ENDDO
    end function kx_derivative


    ! function matrix_out = kz_derivative( matrix_in )
    ! This function multiplies the matrix_in by i kz, which is differenetiating
    ! in z in the Fourier domain
    !
    ! Arguments:
    !   matrix_in : [double/single complex, size (mxf,mzf), Input]
    !               velocity field in the Fourier domain at a single y plane
    ! Output:
    !   matrix_out: [double/single complex, size (mxf,mzf), Output]
    !               z derivative of the velocity field in the Fourier domain
    !               at a single y plane
    function kz_derivative( matrix_in ) result( matrix_out )
        complex(kind=cp), dimension(mxf,mzf), intent(in) :: matrix_in
        complex(kind=cp), dimension(mxf,mzf) :: matrix_out

        complex(kind=cp) :: complex_i
        integer :: ii


        complex_i = ( 0.0_cp, 1.0_cp )

        DO ii = 1, mzf
            matrix_out(:, ii) = matrix_in(:, ii) * kz_full(ii) * complex_i
        ENDDO
    end function kz_derivative




end module wavenumbers