module wavenumbers
    use types, only: sp, dp
    implicit none
    private
#   include "parameters"

    real(kind=cp) :: kx_full, kx_pos, kz_full, kz_pos
    common /wave/ kx_full(nkx_full), kx_pos(nkx_pos), &
                  kz_full(nkz_full), kz_pos(nkz_pos)
    save /wave/

    public :: generate_wavenumbers, kx_derivative, kz_derivative, truncate_kx_2D, truncate_kx_3D


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


    ! subroutine truncate_kx_2D( matrix )
    ! This function truncate the data in kx by setting all data with kx wavenumber
    ! greater than kx_max_truncate to 0 (kx = kx_max_truncate is retained)
    ! If kx_max_truncate < max(kx_pos) then this function does nothing
    !
    ! Arguments:
    !   matrix : [double/single complex, size (mxf,mzf), Input/Output]
    !            Any Forier domain data at a single y plane
    !   onoff  : [logical]
    !            turn on or turn off truncation
    !            if set to false, no truncation is performed regardless of the parameter kx_max_truncate
    !            if set to true , truncation is performed if kx_max_truncate < kx_max
    subroutine truncate_kx_2D( matrix, onoff)
        complex(kind=cp), intent(inout), dimension(mxf, mzf) :: matrix
        logical, intent(in) :: onoff

        real(kind=cp) :: kx_max
        integer :: loc(1), ind


        if ( onoff .eqv. .true. ) then
            ! Max value in the kx wavenumber vector
            kx_max = maxval(kx_pos)

            ! Truncate data if kx_max_truncate is smaller than kx_max
            ! Otherwise, no truncation is needed
            if (kx_max_truncate .lt. kx_max) then
                ! find the index where kx_pos if first larger than kx_max_truncate
                loc = findloc( (kx_pos .gt. kx_max_truncate), .true. )
                ind = loc(1)

                ! set data to 0 for kx greater than the set truncation point
                matrix( ind:mxf, : ) = (0.0_cp, 0.0_cp)
            endif
        endif
    end subroutine truncate_kx_2D


    ! subroutine truncate_kx_3D( matrix )
    ! This function truncate the data in kx by setting all data with kx wavenumber
    ! greater than kx_max_truncate to 0 (kx = kx_max_truncate is retained)
    ! If kx_max_truncate < max(kx_pos) then this function does nothing
    !
    ! Arguments:
    !   matrix : [double/single complex, size (mxf,mzf,myf), Input/Output]
    !            Any Forier domain data at a all y planes
    !   onoff  : [logical]
    !            turn on or turn off truncation
    !            if set to false, no truncation is performed regardless of the parameter kx_max_truncate
    !            if set to true , truncation is performed if kx_max_truncate < kx_max
    subroutine truncate_kx_3D( matrix , onoff)
        complex(kind=cp), intent(inout), dimension(mxf, mzf, myf) :: matrix
        logical, intent(in) :: onoff

        real(kind=cp) :: kx_max
        integer :: loc(1), ind


        if ( onoff .eqv. .true. ) then
            ! Max value in the kx wavenumber vector
            kx_max = maxval(kx_pos)

            ! Truncate data if kx_max_truncate is smaller than kx_max
            ! Otherwise, no truncation is needed
            if (kx_max_truncate .lt. kx_max) then
                ! find the index where kx_pos if first larger than kx_max_truncate
                loc = findloc( (kx_pos .gt. kx_max_truncate), .true. )
                ind = loc(1)

                ! set data to 0 for kx greater than the set truncation point
                matrix( ind:mxf, :, : ) = (0.0_cp, 0.0_cp)
            endif
        endif
    end subroutine truncate_kx_3D


end module wavenumbers