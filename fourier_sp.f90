module fourier_sp
    use, intrinsic :: iso_c_binding
    use types, only: sp, dp
    implicit none
    private
#   include "parameters"
    include 'fftw3.f03'

    type(C_PTR) :: plan_ifftx, plan_ifftz, plan_fft2, plan_ifft2, plan_fftt

    public :: fft_plan, fft2, ifft2, ifftx, ifftz, fftt


contains
    ! subroutine fft_plan
    ! fft_plan_dp for double precision computation
    ! fft_plan_dp for single precision computation
    !
    ! This function generates 4 fftw plans:
    ! plan_ifftx: inverse fft for x         complex vector size mgalx -> complex vector size mgalx
    ! plan_ifftz: inverse fft for z         complex vector size mgalz -> complex vector size mgalz
    ! plan_fft2 : fft in both x z           real matrix size mgalx,mgalz -> complex matrix size mgalx/2+1,mgalz
    ! plan_ifft2: inverse fft in both x z   complex matrix size mgalx/2+1,mgalz -> real matrix size mgalx,mgalz
    subroutine fft_plan
        complex(C_FLOAT_COMPLEX), dimension(mgalx) :: vector_x
        complex(C_FLOAT_COMPLEX), dimension(mgalz) :: vector_z
        complex(C_FLOAT_COMPLEX), dimension(mgalx/2+1,mgalz) :: matrix_2d_comp
        real   (C_FLOAT        ), dimension(mgalx    ,mgalz) :: matrix_2d_real
        complex(C_FLOAT_COMPLEX), dimension(nt) :: vector_t


        ! Generate plans for ifft in x and z seperately
        ! These two are inplace transforms
        plan_ifftx = fftwf_plan_dft_1d(mgalx, vector_x,vector_x, FFTW_BACKWARD, FFTW_PATIENT )
        plan_ifftz = fftwf_plan_dft_1d(mgalz, vector_z,vector_z, FFTW_BACKWARD, FFTW_PATIENT )

        ! Note: the r2c and c2r versions of the transform uses the hermitian symmetry
        ! The symmetric part of the kx wavenumber is removed
        ! Therefore the real data has size (mgalx, mgalz)
        ! the complex data has size (mgalx/2+1, mgalz)
        ! See FFTW documentation for more details:
        ! https://www.fftw.org/fftw3_doc/One_002dDimensional-DFTs-of-Real-Data.html
        ! https://www.fftw.org/fftw3_doc/Multi_002dDimensional-DFTs-of-Real-Data.html#Multi_002dDimensional-DFTs-of-Real-Data
        ! https://www.fftw.org/fftw3_doc/Reversing-array-dimensions.html
        plan_fft2  = fftwf_plan_dft_r2c_2d(mgalz,mgalx, matrix_2d_real,matrix_2d_comp, FFTW_PATIENT )
        plan_ifft2 = fftwf_plan_dft_c2r_2d(mgalz,mgalx, matrix_2d_comp,matrix_2d_real, FFTW_PATIENT )

        ! Generate plan for time FFT
        plan_fftt = fftwf_plan_dft_1d(nt, vector_t,vector_t, FFTW_FORWARD, FFTW_PATIENT )
    end subroutine fft_plan


    ! subroutine zeropad_2d( matrix_in, matrix_out )
    ! This function zeros pads the high frequencies of kx and kz
    ! Note that it does not fill in the negative kx part using hermitian
    ! symmetry as the c2r fftw routines do not include the symmetric part
    !
    ! Arguments:
    !   matrix_in:  [double/single complex, size (mxf,mzf), Input]
    !               the fourier coefficients at a single y plane
    !   matrix_out: [double/single complex, size (mgalx/2+1,mgalz), Output]
    !               the fourier coefficients at a single y plane padded with
    !               zeros for the high frequencies
    subroutine zeropad_2d( matrix_in, matrix_out )
        complex(kind=sp), intent( in), dimension(mxf,mzf) :: matrix_in
        complex(kind=sp), intent(out), dimension(mgalx/2+1,mgalz) :: matrix_out

        integer :: nkx, nkz, nposkz, nnegkz, nx, nz


        ! initialize output matrix to 0
        matrix_out = (0.0_sp, 0.0_sp)

        nkx = mxf
        nkz = mzf
        nposkz = ( nkz - 1 )/2 + 1 ! number of positive kz (including 0)
        nnegkz = ( nkz - 1 )/2     ! number of negative kz
        nx  = mgalx
        nz  = mgalz

        ! kx >= 0, kz >= 0
        matrix_out(1:nkx, 1:nposkz) = matrix_in(1:nkx, 1:nposkz)
        ! kx >= 0, kz <  0
        matrix_out(1:nkx, nz-nnegkz+1:nz) = matrix_in(1:nkx, nkz-nnegkz+1:nkz)
    end subroutine zeropad_2d


    ! subroutine removezeropad_2d( matrix_in, matrix_out )
    ! This function removes the zeros pad at the high frequencies of kx and kz
    !
    ! Arguments:
    !   matrix_in:  [double/single complex, size (mgalx/2+1,mgalz), Input]
    !               the fourier coefficients at a single y plane padded with
    !               zeros for the high frequencies without negative kx
    !               wavenumbers
    !   matrix_out: [double/single complex, size (mxf,mzf), Output]
    !               the fourier coefficients at a single y plane with zero pad
    !               and kx < 0 wavenumbers removed
    subroutine removezeropad_2d( matrix_in, matrix_out )
        complex(kind=sp), intent( in), dimension(mgalx/2+1,mgalz) :: matrix_in
        complex(kind=sp), intent(out), dimension(mxf,mzf) :: matrix_out

        integer :: nkx, nkz, nposkz, nnegkz, nx, nz


        nkx = mxf
        nkz = mzf
        nposkz = ( nkz - 1 )/2 + 1 ! number of positive kz (including 0)
        nnegkz = ( nkz - 1 )/2     ! number of negative kz
        nx  = mgalx
        nz  = mgalz

        ! kx >= 0, kz >= 0
        matrix_out(1:nkx, 1:nposkz) = matrix_in(1:nkx, 1:nposkz)
        ! kx >= 0, kz <  0
        matrix_out(1:nkx, nposkz+1:nkz) = matrix_in(1:nkx, nz-nnegkz+1:nz)
    end subroutine removezeropad_2d


    ! subroutine zeropad_x( matrix_in, matrix_out )
    ! This function rearranges the wavenumbers in preperation for a fourier
    ! transform in the x direction only.
    ! It rearranges the wavenumbers using hermitian symmetry so that kx now
    ! contains the full range with both positive and negative and kz only
    ! contains the non-negative wavenumbers. Then zero pad the high kx
    ! wavenumbers.
    !
    ! Arguments:
    !   matrix_in:  [double/single complex, size (mxf,mzf), Input]
    !               the fourier coefficients at a single y plane
    !   matrix_out: [double/single complex, size (mgalx, nkz_pos), Output]
    !               the fourier coefficients at a single y plane padded with 0
    !               for the high kx wavenumbers, and kx now contains the full
    !               range while kz only contains the non-negative wavenumbers
    subroutine zeropad_x( matrix_in, matrix_out )
        complex(kind=sp), intent( in), dimension(mxf,mzf) :: matrix_in
        complex(kind=sp), intent(out), dimension(mgalx, nkz_pos) :: matrix_out


        ! initialize to 0
        matrix_out = (0.0_sp, 0.0_sp)

        ! kx >= 0, kz >= 0
        ! exact 1 to 1 mapping
        matrix_out(1:nkx_pos, 1:nkz_pos) = matrix_in(1:nkx_pos, 1:nkz_pos)
        ! kx <  0, kz =  0
        ! this is the complex conjugate of kx > 0, kz = 0
        matrix_out(mgalx-nkx_pos+2:mgalx, 1) = CONJG( matrix_in(nkx_pos:2:-1,1) )
        ! kx <  0, kz >  0
        ! map from  kx > 0, kz < 0 quadrant using hermitian symmetry
        matrix_out(mgalx-nkx_pos+2:mgalx, 2:nkz_pos) = CONJG( matrix_in(nkx_pos:2:-1, nkz_full:nkz_pos+1:-1) )

    end subroutine zeropad_x


    ! subroutine zeropad_z( matrix_in, matrix_out )
    ! This function rearranges the wavenumbers in preperation for a fourier
    ! transform in the z direction only.
    ! It zero pad the high kz wavenumbers, while leaving the kx wavenumbers
    ! unchanged
    !
    ! Arguments:
    !   matrix_in:  [double/single complex, size (mxf,mzf), Input]
    !               the fourier coefficients at a single y plane
    !   matrix_out: [double/single complex, size (mgalz,nkx_pos), Output]
    !               the fourier coefficients at a single y plane padded with 0
    !               for the high kz wavenumbers, and with a transpose built in
    !               to exchange the kx and kz dimesions for improved indexing
    !               effeciency
    subroutine zeropad_z( matrix_in, matrix_out )
        complex(kind=sp), intent( in), dimension(mxf,mzf) :: matrix_in
        complex(kind=sp), intent(out), dimension(mgalz,nkx_pos) :: matrix_out

        integer :: ii


        ! initialize to 0
        matrix_out = (0.0_sp, 0.0_sp)

        DO ii = 1,nkx_pos
            ! kx >= 0, kz >= 0
            ! exact 1 to 1 mapping
            matrix_out(1:nkz_pos, ii) = matrix_in(ii, 1:nkz_pos)
            ! kx >= 0, kz <  0
            ! shifted 1 to 1 mapping
            matrix_out(mgalz-nkz_pos+2:mgalz, ii) = matrix_in(ii, nkz_pos+1:nkz_full)
        ENDDO

    end subroutine zeropad_z


    ! subroutine fft2( matrix_in, matrix_out )
    ! Computes the 2d fft for a 2d matrix (single y location)
    !
    ! Arguments:
    !   matrix_in:  [double/single, size (mgalx,mgalz), Input]
    !               physical domain data
    !   matrix_out: [double/single complex, size (mxf,mzf), Output]
    !               fourier domain data, with only positive kx wavenumbers
    !               and high wavenumber zero padding removed, and with
    !               normalization factor applied
    subroutine fft2( matrix_in, matrix_out )
        real   (kind=sp), intent( in), dimension(mgalx,mgalz) :: matrix_in
        complex(kind=sp), intent(out), dimension(mxf,mzf) :: matrix_out

        real   (kind=sp), dimension(mgalx,mgalz) :: temp1
        complex(kind=sp), dimension(mgalx/2+1,mgalz) :: temp2


        temp1 = matrix_in
        call fftwf_execute_dft_r2c(plan_fft2, temp1, temp2)
        call removezeropad_2d( temp2, matrix_out )
        ! normalization factor
        matrix_out = matrix_out/real(mgalx,sp)/real(mgalz,sp)
    end subroutine fft2


    ! subroutine ifft2( matrix_in, matrix_out )
    ! Computes the 2d ifft for a 2d matrix (single y location)
    !
    ! Arguments:
    !   matrix_in:  [double/single complex, size (mxf,mzf)), Input]
    !               fourier domain data
    !   matrix_out: [double/single, size (mgalx,mgalz), Output]
    !               physical domain data
    subroutine ifft2( matrix_in, matrix_out )
        complex(kind=sp), intent( in), dimension(mxf,mzf) :: matrix_in
        real   (kind=sp), intent(out), dimension(mgalx,mgalz) :: matrix_out

        complex(kind=sp), dimension(mgalx/2+1,mgalz) :: temp


        call zeropad_2d( matrix_in, temp )
        call fftwf_execute_dft_c2r(plan_ifft2, temp, matrix_out)
    end subroutine ifft2


    ! subroutine ifftx( matrix_in, matrix_out )
    ! Computes the ifft in x only for a 2d matrix (single y location) for only
    ! the non-negative kz wavenumbers
    !
    ! Arguments:
    !   matrix_in:  [double/single complex, size (mxf,mzf), Input]
    !               fourier domain data
    !   matrix_out: [double/single complex, size (mgalx, nkz_pos), Output]
    !               ifft in x only, for the non-negative kz wavenumbers
    subroutine ifftx( matrix_in, matrix_out )
        complex(kind=sp), intent( in), dimension(mxf,mzf) :: matrix_in
        complex(kind=sp), intent(out), dimension(mgalx, nkz_pos) :: matrix_out

        integer :: ii


        ! zero pad in x
        call zeropad_x( matrix_in, matrix_out )
        ! for each kz, ifft in x
        DO ii = 1, nkz_pos
            call fftwf_execute_dft(plan_ifftx, matrix_out(:,ii), matrix_out(:,ii) )
        ENDDO
    end subroutine ifftx


    ! subroutine ifftz( matrix_in, matrix_out )
    ! Computes the ifft in z only for a 2d matrix (single y location) for only
    ! the non-negative kx wavenumbers
    !
    ! Arguments:
    !   matrix_in:  [double/single complex, size (mxf,mzf), Input]
    !               fourier domain data
    !   matrix_out: [double/single complex, size (mgalz, nkx_pos), Output]
    !               ifft in z only, for the non-negative kx wavenumbers
    !               Note that this has a built in transpose that exchanges the
    !               kx dimension with the z dimension for improved indexing
    !               effeciency
    subroutine ifftz( matrix_in, matrix_out )
        complex(kind=sp), intent( in), dimension(mxf,mzf) :: matrix_in
        complex(kind=sp), intent(out), dimension(mgalz,nkx_pos) :: matrix_out

        integer :: ii


        ! zero pad in z
        call zeropad_z( matrix_in, matrix_out )
        ! for each kx, ifft in z
        DO ii = 1, nkx_pos
            call fftwf_execute_dft(plan_ifftz, matrix_out(:,ii), matrix_out(:,ii) )
        ENDDO
    end subroutine ifftz


    ! subroutine fftt( matrix )
    ! Computes the fft in time for a 3d matrix
    !
    ! Arguments
    !   matrix: [single complex, size (:,:,nt), Input/Output]
    !           Time domain data as the input, and updated to become frequency domain data
    !           The third dimension needs to be time, and updated to become frequency
    subroutine fftt( matrix )
        complex(kind=sp), intent(inout), dimension(:,:,:) :: matrix

        complex(kind=sp), dimension(nt) :: vec

        integer :: ii, jj, kk
        integer :: size_matrix(3), size_dim1, size_dim2

        ! get matrix dimensions
        size_matrix = shape(matrix)
        size_dim1   = size_matrix(1)
        size_dim2   = size_matrix(2)

        ! loop over dimension 1 and 2 and transform dimension 3 (which is time)
        DO jj = 1, size_dim2
            DO ii = 1, size_dim1
                DO kk = 1,nt
                    vec(kk) = matrix(ii,jj,kk)
                ENDDO

                ! perform a inplace transform
                call fftwf_execute_dft( plan_fftt, vec, vec)

                DO kk = 1,nt
                    matrix(ii,jj,kk) = vec(kk)/real(nt,sp)
                ENDDO
            ENDDO
        ENDDO
    end subroutine fftt


end module fourier_sp