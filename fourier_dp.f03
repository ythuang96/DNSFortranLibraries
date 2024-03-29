module fourier_dp
    use, intrinsic :: iso_c_binding
    use types, only: sp, dp
    implicit none
    private
#   include "parameters"
    include 'fftw3.f03'

    type(C_PTR) :: plan_ifftx, plan_fftx, plan_ifftz, plan_fftz
    type(C_PTR) :: plan_fft2, plan_ifft2, plan_fftt, plan_fftt_c2c

    complex(C_DOUBLE_COMPLEX), dimension(mgalx) :: vector_ifftx
    complex(C_DOUBLE_COMPLEX), dimension(mgalz) :: vector_ifftz
    complex(C_DOUBLE_COMPLEX), dimension(mgalx/2+1,mgalz) :: matrix_2d_comp
    real   (C_DOUBLE        ), dimension(mgalx    ,mgalz) :: matrix_2d_real
    real   (C_DOUBLE        ), dimension(ntseg    ) :: vector_t
    complex(C_DOUBLE_COMPLEX), dimension(ntseg/2+1) :: vector_om

    complex(C_DOUBLE_COMPLEX), dimension(ntseg) :: vector_tom_c2c

    public :: fft_plan, fft2, ifft2, ifftx, ifftz, fft_tseg_r2c, fft_tseg_c2c, fft_tseg_r2c_halfom
    public :: ifftz_h_vec, ifftx_h_vec


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
    ! plan_fftt : fft in time               real vector size nt -> complex vector size nt/2+1 (nt = nom)
    subroutine fft_plan
        ! Generate plans for ifft in x and z seperately
        ! These two are inplace transforms (the two pointers can be equal, indicating an in-place transform, ignore compiler warning )
        plan_ifftx = fftw_plan_dft_1d(mgalx, vector_ifftx,vector_ifftx, FFTW_BACKWARD, FFTW_PATIENT )
        plan_fftx  = fftw_plan_dft_1d(mgalx, vector_ifftx,vector_ifftx, FFTW_FORWARD , FFTW_PATIENT )
        plan_ifftz = fftw_plan_dft_1d(mgalz, vector_ifftz,vector_ifftz, FFTW_BACKWARD, FFTW_PATIENT )
        plan_fftz  = fftw_plan_dft_1d(mgalz, vector_ifftz,vector_ifftz, FFTW_FORWARD , FFTW_PATIENT )

        ! Note: the r2c and c2r versions of the transform uses the hermitian symmetry
        ! The symmetric part of the kx wavenumber is removed
        ! Therefore the real data has size (mgalx, mgalz)
        ! the complex data has size (mgalx/2+1, mgalz)
        ! See FFTW documentation for more details:
        ! https://www.fftw.org/fftw3_doc/One_002dDimensional-DFTs-of-Real-Data.html
        ! https://www.fftw.org/fftw3_doc/Multi_002dDimensional-DFTs-of-Real-Data.html#Multi_002dDimensional-DFTs-of-Real-Data
        ! https://www.fftw.org/fftw3_doc/Reversing-array-dimensions.html
        plan_fft2  = fftw_plan_dft_r2c_2d(mgalz,mgalx, matrix_2d_real,matrix_2d_comp, FFTW_PATIENT )
        plan_ifft2 = fftw_plan_dft_c2r_2d(mgalz,mgalx, matrix_2d_comp,matrix_2d_real, FFTW_PATIENT )

        ! Generate plan for time FFT
        ! that is NOT consistent with the reolvent definition with a negative sign
        ! need to flip the omega in later code
        plan_fftt = fftw_plan_dft_r2c_1d( ntseg, vector_t,vector_om, FFTW_PATIENT )

        ! Generate plan for time FFT of a complex data
        ! use FFTW_BACKWARD because that will automatically result in an omega
        ! that is consistent with the reolvent definition with a negative sign
        plan_fftt_c2c = fftw_plan_dft_1d( ntseg, vector_tom_c2c, vector_tom_c2c, FFTW_BACKWARD, FFTW_PATIENT)
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
        complex(kind=dp), intent( in), dimension(mxf,mzf) :: matrix_in
        complex(kind=dp), intent(out), dimension(mgalx/2+1,mgalz) :: matrix_out

        integer :: nkx, nkz, nposkz, nnegkz, nx, nz


        ! initialize output matrix to 0
        matrix_out = (0.0_dp, 0.0_dp)

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
        complex(kind=dp), intent( in), dimension(mgalx/2+1,mgalz) :: matrix_in
        complex(kind=dp), intent(out), dimension(mxf,mzf) :: matrix_out

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
        complex(kind=dp), intent( in), dimension(mxf,mzf) :: matrix_in
        complex(kind=dp), intent(out), dimension(mgalx, nkz_pos) :: matrix_out


        ! initialize to 0
        matrix_out = (0.0_dp, 0.0_dp)

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


    ! subroutine zeropad_x_vec( vec_in, vec_out )
    ! This function rearranges the wavenumbers in preperation for a fourier
    ! transform in the x direction only.
    ! It zero pad the high kx wavenumber
    !
    ! Arguments:
    !   vec_in:  [double complex, size (nkx_full), Input]
    !            the fourier coefficients, with both positive and negative kx
    !   vec_out: [double complex, size (mgalx), Output]
    !            the fourier coefficients padded with 0 for the high kx wavenumbers
    subroutine zeropad_x_vec( vec_in, vec_out )
        complex(kind=dp), intent( in), dimension(nkx_full) :: vec_in
        complex(kind=dp), intent(out), dimension(mgalx) :: vec_out

        ! initialize to 0
        vec_out = (0.0_dp, 0.0_dp)

        ! kx >= 0
        ! exact 1 to 1 mapping
        vec_out(1:nkx_pos) = vec_in(1:nkx_pos)
        ! kx <  0
        ! shifted 1 to 1 mapping
        vec_out(mgalx-nkx_pos+2:mgalx) = vec_in(nkx_pos+1:nkx_full)
    end subroutine zeropad_x_vec


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
        complex(kind=dp), intent( in), dimension(mxf,mzf) :: matrix_in
        complex(kind=dp), intent(out), dimension(mgalz,nkx_pos) :: matrix_out

        integer :: ii


        ! initialize to 0
        matrix_out = (0.0_dp, 0.0_dp)

        DO ii = 1,nkx_pos
            ! kx >= 0, kz >= 0
            ! exact 1 to 1 mapping
            matrix_out(1:nkz_pos, ii) = matrix_in(ii, 1:nkz_pos)
            ! kx >= 0, kz <  0
            ! shifted 1 to 1 mapping
            matrix_out(mgalz-nkz_pos+2:mgalz, ii) = matrix_in(ii, nkz_pos+1:nkz_full)
        ENDDO

    end subroutine zeropad_z


    ! subroutine zeropad_z_vec( vec_in, vec_out )
    ! This function rearranges the wavenumbers in preperation for a fourier
    ! transform in the z direction only.
    ! It zero pad the high kz wavenumber
    !
    ! Arguments:
    !   vec_in:  [double complex, size (mzf), Input]
    !            the fourier coefficients
    !   vec_out: [double complex, size (mgalz), Output]
    !            the fourier coefficients padded with 0 for the high kz wavenumbers
    subroutine zeropad_z_vec( vec_in, vec_out )
        complex(kind=dp), intent( in), dimension(mzf) :: vec_in
        complex(kind=dp), intent(out), dimension(mgalz) :: vec_out

        ! initialize to 0
        vec_out = (0.0_dp, 0.0_dp)

        ! kx >= 0, kz >= 0
        ! exact 1 to 1 mapping
        vec_out(1:nkz_pos) = vec_in(1:nkz_pos)
        ! kx >= 0, kz <  0
        ! shifted 1 to 1 mapping
        vec_out(mgalz-nkz_pos+2:mgalz) = vec_in(nkz_pos+1:nkz_full)
    end subroutine zeropad_z_vec


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
        real   (kind=dp), intent( in), dimension(mgalx,mgalz) :: matrix_in
        complex(kind=dp), intent(out), dimension(mxf,mzf) :: matrix_out

        matrix_2d_real = matrix_in
        call fftw_execute_dft_r2c(plan_fft2, matrix_2d_real, matrix_2d_comp)
        call removezeropad_2d( matrix_2d_comp, matrix_out )
        ! normalization factor
        matrix_out = matrix_out/real(mgalx,dp)/real(mgalz,dp)
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
        complex(kind=dp), intent( in), dimension(mxf,mzf) :: matrix_in
        real   (kind=dp), intent(out), dimension(mgalx,mgalz) :: matrix_out

        call zeropad_2d( matrix_in, matrix_2d_comp )
        call fftw_execute_dft_c2r(plan_ifft2, matrix_2d_comp, matrix_2d_real)
        matrix_out = matrix_2d_real
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
        complex(kind=dp), intent( in), dimension(mxf,mzf) :: matrix_in
        complex(kind=dp), intent(out), dimension(mgalx, nkz_pos) :: matrix_out

        integer :: ii


        ! zero pad in x
        call zeropad_x( matrix_in, matrix_out )
        ! for each kz, ifft in x
        DO ii = 1, nkz_pos
            vector_ifftx = matrix_out(:,ii)
            call fftw_execute_dft(plan_ifftx, vector_ifftx, vector_ifftx )
            matrix_out(:,ii) = vector_ifftx
        ENDDO
    end subroutine ifftx


    ! subroutine ifftx_h_vec( vec_in, vec_out )
    ! Computes the special ifft for h in x only for a 1d matrix
    ! The special ifft for h is in fact just a fft
    !
    ! Arguments:
    !   vec_in:  [double complex, size (nkx_full), Input]
    !            fourier domain data, with both positive and negative kx
    !   vec_out: [double complex, size (mgalx), Output]
    !            special ifft (fft) in x only
    subroutine ifftx_h_vec( vec_in, vec_out )
        complex(kind=dp), intent( in), dimension(nkx_full) :: vec_in
        complex(kind=dp), intent(out), dimension(mgalx) :: vec_out

        ! zero pad in x
        call zeropad_x_vec( vec_in, vec_out )
        ! special ifft in z for h ( the special ifft for h is simply a fft)
        vector_ifftx = vec_out
        call fftw_execute_dft(plan_fftx, vector_ifftx, vector_ifftx )
        vec_out = vector_ifftx
    end subroutine ifftx_h_vec


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
        complex(kind=dp), intent( in), dimension(mxf,mzf) :: matrix_in
        complex(kind=dp), intent(out), dimension(mgalz,nkx_pos) :: matrix_out

        integer :: ii


        ! zero pad in z
        call zeropad_z( matrix_in, matrix_out )
        ! for each kx, ifft in z
        DO ii = 1, nkx_pos
            vector_ifftz = matrix_out(:,ii)
            call fftw_execute_dft(plan_ifftz, vector_ifftz, vector_ifftz )
            matrix_out(:,ii) = vector_ifftz
        ENDDO
    end subroutine ifftz


    ! subroutine ifftz_h_vec( vec_in, vec_out )
    ! Computes the special ifft for h in z only for a 1d matrix
    ! The special ifft for h is in fact just a fft
    !
    ! Arguments:
    !   vec_in:  [double complex, size (mzf), Input]
    !            fourier domain data
    !   vec_out: [double complex, size (mgalz), Output]
    !            special ifft (fft) in z only
    subroutine ifftz_h_vec( vec_in, vec_out )
        complex(kind=dp), intent( in), dimension(mzf) :: vec_in
        complex(kind=dp), intent(out), dimension(mgalz) :: vec_out

        ! zero pad in z
        call zeropad_z_vec( vec_in, vec_out )
        ! special ifft in z for h ( the special ifft for h is simply a fft)
        vector_ifftz = vec_out
        call fftw_execute_dft(plan_fftz, vector_ifftz, vector_ifftz )
        vec_out = vector_ifftz
    end subroutine ifftz_h_vec


    ! subroutine fft_tseg_r2c( vector_in, matrix_out, WindowFunction )
    ! Computes the segmented and windowed fft in time for a temporal real vector
    !
    ! Arguments
    !   vector_in    : [double real, size (nt), Input]
    !                  Time domain data as the input
    !   matrix_out   : [double complex, size (nom,nseg), Output]
    !                  Frequency domain data. The first dimension is omega,
    !                  with the redundant (conjugate) data included, organized as ... -om2, -om1, 0, om1, om2 ...
    !                  Note that omega is consistent with the resolvent definition with a negative sign
    !                  The second dimension is the segment number
    !   WindowFuncion: [double real, size (ntseg), Input]
    !                  A real vector with same length as the temporal segment, for the window function
    subroutine fft_tseg_r2c( vector_in, matrix_out, WindowFunction )
        real   (kind=dp), intent( in), dimension(nt)       :: vector_in
        complex(kind=dp), intent(out), dimension(nom,nseg) :: matrix_out
        real   (kind=dp), intent( in), dimension(ntseg)    :: WindowFunction

        integer :: ii, jj, ind_start

        ! loop over all the welch segments
        DO ii = 1, nseg

            ! starting index of this time segment
            ind_start = (ii-1)*( ntseg - noverlap )
            ! The data of this time segment multiplied with the window function
            vector_t = vector_in( ind_start+1:ind_start+ntseg ) * WindowFunction

            ! perform transform
            ! Note that omega is NOT consistent with the resolvent definition with a negative sign
            ! as the fft plan uses forward fft
            call fftw_execute_dft_r2c( plan_fftt, vector_t, vector_om)

            ! Perform the truncation in omega and copy the transformed vector to the output data matrix
            ! Since the vector_om has a missing negative sign
            ! it is organized as: 0, -om1, -om2, ...
            ! The output should be organized as: ... -om2, -om1, 0, om1, om2 ...

            ! omega < 0
            DO jj = 1, (nom-1)/2
                matrix_out(jj, ii) = vector_om( (nom-1)/2+2-jj ) /real(ntseg,dp)
            ENDDO
            ! omega = 0
            matrix_out( (nom-1)/2+1, ii ) = vector_om(1) /real(ntseg,dp)
            ! omega > 0
            DO jj = 1,(nom-1)/2
                matrix_out((nom-1)/2+1+jj,ii) = conjg( vector_om(jj+1) )/real(ntseg,dp)
                ! conjg to conform to the resolvent standart with negative sign for omega
            ENDDO

        ENDDO

    end subroutine fft_tseg_r2c


    ! subroutine fft_tseg_r2c_halfom( vector_in, matrix_out, WindowFunction )
    ! Computes the segmented and windowed fft in time for a temporal real vector
    !
    ! Arguments
    !   vector_in    : [double real, size (nt), Input]
    !                  Time domain data as the input
    !   matrix_out   : [double complex, size ((nom+1)/2,nseg), Output]
    !                  Frequency domain data. The first dimension is omega,
    !                  with the redundant (conjugate) data excluded, organized as 0, om1, om2 ...
    !                  Note that omega is consistent with the resolvent definition with a negative sign
    !                  The second dimension is the segment number
    !   WindowFuncion: [double real, size (ntseg), Input]
    !                  A real vector with same length as the temporal segment, for the window function
    subroutine fft_tseg_r2c_halfom( vector_in, matrix_out, WindowFunction )
        real   (kind=dp), intent( in), dimension(nt)             :: vector_in
        complex(kind=dp), intent(out), dimension((nom+1)/2,nseg) :: matrix_out
        real   (kind=dp), intent( in), dimension(ntseg)          :: WindowFunction

        integer :: ii, ind_start

        ! loop over all the welch segments
        DO ii = 1, nseg

            ! starting index of this time segment
            ind_start = (ii-1)*( ntseg - noverlap )
            ! The data of this time segment multiplied with the window function
            vector_t = vector_in( ind_start+1:ind_start+ntseg ) * WindowFunction

            ! perform transform
            ! Note that omega is NOT consistent with the resolvent definition with a negative sign
            ! as the fft plan uses forward fft
            call fftw_execute_dft_r2c( plan_fftt, vector_t, vector_om)

            ! Perform the truncation in omega and copy the transformed vector to the output data matrix
            ! Since the vector_om has a missing negative sign
            ! it is organized as: 0, -om1, -om2, ...
            ! The output should be organized as: 0, om1, om2 ...

            ! omega = 0
            matrix_out( 1, ii ) = vector_om(1) /real(ntseg,dp)
            ! omega > 0
            matrix_out( 2:(nom+1)/2, ii) = conjg( vector_om(2:(nom+1)/2) )/real(ntseg,dp)
            ! conjg to conform to the resolvent standart with negative sign for omega

        ENDDO

    end subroutine fft_tseg_r2c_halfom


    ! subroutine fft_tseg_c2c( vector_in, matrix_out, WindowFuncion)
    ! Computes the segmented and windowed fft in time for a temporal complex vector
    !
    ! Arguments
    !   vector_in    : [double complex, size (nt), Input]
    !                  Time domain data as the input
    !   matrix_out   : [double complex, size (nom,nseg), Output]
    !                  Frequency domain data. The first dimension is omega,
    !                  organized as ... -om2, -om1, 0, om1, om2
    !                  Note that omega is consistent with the resolvent definition with a negative sign
    !                  The second dimension is the segment number
    !   WindowFuncion: [double real, size (ntseg), Input]
    !                  A real vector with same length as the temporal segment, for the window function
    subroutine fft_tseg_c2c( vector_in, matrix_out, WindowFunction)
        complex(kind=dp), intent( in), dimension(nt)       :: vector_in
        complex(kind=dp), intent(out), dimension(nom,nseg) :: matrix_out
        real   (kind=dp), intent( in), dimension(ntseg)    :: WindowFunction

        integer :: ii, ind_start

        ! loop over all the welch segments
        DO ii = 1, nseg

            ! starting index of this time segment
            ind_start = (ii-1)*( ntseg - noverlap )
            ! The data of this time segment multiplied with the window function
            vector_tom_c2c = vector_in( ind_start+1:ind_start+ntseg ) * WindowFunction

            ! perform transform
            ! Note that omega is consistent with the resolvent definition with a negative sign
            ! as the fft plan uses backwards fft
            call fftw_execute_dft( plan_fftt_c2c, vector_tom_c2c, vector_tom_c2c )
            ! The resulting omega is organized as 0, om1, om2, ... -om2, -om1

            ! Perform the truncation in omega and copy the transformed vector to the output data matrix
            ! omega < 0
            matrix_out(     1:(nom-1)/2, ii) = vector_tom_c2c( ntseg-(nom-1)/2+1:ntseg ) /real(ntseg,dp)
            ! omega >= 0
            matrix_out( (nom-1)/2+1:nom, ii) = vector_tom_c2c(      1:(nom-1)/2+1      ) /real(ntseg,dp)
            ! divide by ntseg as fftw does not perform normalization

        ENDDO

    end subroutine fft_tseg_c2c


end module fourier_dp