module ComputeTriadicProjection
    use types, only: sp, dp
    use time, only: timer_continue, timer_stop

    implicit none
    private
#   include "parameters"

    ! Global Timer Variables
    real(kind=dp) :: timer_total, timer_loop, timer_read, timer_distribute, timer_dy, &
                     timer_change, timer_f, timer_pkx, timer_pkz, &
                     timer_pkx1, timer_pkx2, timer_pkz1, timer_pkz2
    common/timers/   timer_total, timer_loop, timer_read, timer_distribute, timer_dy, &
                     timer_change, timer_f, timer_pkx, timer_pkz, &
                     timer_pkx1, timer_pkx2, timer_pkz1, timer_pkz2
    save  /timers/

    public :: compute_Pkx_Pkz

contains
    subroutine compute_Pkx_Pkz(uf,vf,wf,dudyf,dvdyf,dwdyf, yplane, Px_kx, Py_kx, Pz_kx, Px_kz, Py_kz, Pz_kz)
        use ComputeForcing, only: compute_forcing_meansubtract
        ! Input/Output
        complex(kind=cp), intent( in), dimension(mxf,mzf) :: uf,vf,wf, dudyf,dvdyf,dwdyf
        integer, intent(in) :: yplane
        complex(kind=cp), intent(inout), dimension(nkx_pos,nkx_full) :: Px_kx, Py_kx, Pz_kx
        complex(kind=cp), intent(inout), dimension(nkz_pos,nkz_full) :: Px_kz, Py_kz, Pz_kz
        ! Forcing
        complex(kind=cp), dimension(mxf,mzf) :: fxf,fyf,fzf


        ! Compute forcing
        call timer_continue( timer_f )
        call compute_forcing_meansubtract(uf,vf,wf,dudyf,dvdyf,dwdyf, yplane, fxf,fyf,fzf)
        call timer_stop( timer_f )

        ! Compute P_kx
        call timer_continue( timer_pkx )
        call compute_Pkx(uf,vf,wf,dudyf,dvdyf,dwdyf,fxf,fyf,fzf, Px_kx, Py_kx, Pz_kx)
        call timer_stop( timer_pkx )

        ! Compute P_kz
        call timer_continue( timer_pkz )
        call compute_Pkz(uf,vf,wf,dudyf,dvdyf,dwdyf,fxf,fyf,fzf, Px_kz, Py_kz, Pz_kz)
        call timer_stop( timer_pkz )

    end subroutine compute_Pkx_Pkz


    ! subroutine compute_Pkx(uf,vf,wf,dudyf,dvdyf,dwdyf,fxf,fyf,fzf, Px_kx, Py_kx, Pz_kx)
    ! This function computes the projection coefficient P_kx at a given yplane
    !
    ! Arguments:
    !   uf,vf,wf:            [double/single complex, Size (mxf,mzf), Input]
    !                        velocity fields at a single y plane
    !   dudyf,dvdyf,dwdyf:   [double/single complex, Size (mxf,mzf), Input]
    !                        y derivatives of velocity fields at a single y plane
    !   fxf,fyf,fzf:         [double/single complex, Size (mxf,mzf), Input]
    !                        mean subtracted non-linear forcing at a single y plane
    !   Px_kx, Py_kx, Pz_kx: [double/single complex, Size (nkx_pos,nkx_full), Output]
    !                        The computed projection coefficients
    subroutine compute_Pkx(uf,vf,wf,dudyf,dvdyf,dwdyf,fxf,fyf,fzf, Px_kx, Py_kx, Pz_kx)
        use FFT_PRECISION_CONTROL, only: ifftz
        use wavenumbers, only: kx_derivative, kz_derivative
        ! Input/Outputs
        complex(kind=cp), intent( in), dimension(mxf,mzf) :: uf,vf,wf, dudyf,dvdyf,dwdyf, fxf,fyf,fzf
        complex(kind=cp), intent(inout), dimension(nkx_pos,nkx_full) :: Px_kx, Py_kx, Pz_kx

        ! Results after ifft in z
        ! Note that these all have a built in transpose that exchanges the
        ! kx dimension with the z dimension for improved indexing effeciency
        complex(kind=cp), dimension(mgalz,nkx_pos ) :: uf_t, vf_t, wf_t
        complex(kind=cp), dimension(mgalz,nkx_full) :: fx_t, fy_t, fz_t
        complex(kind=cp), dimension(mgalz,nkx_full) :: dudx_t, dvdx_t, dwdx_t
        complex(kind=cp), dimension(mgalz,nkx_full) :: dudy_t, dvdy_t, dwdy_t
        complex(kind=cp), dimension(mgalz,nkx_full) :: dudz_t, dvdz_t, dwdz_t
        ! Temp variables for the loop
        complex(kind=cp), dimension(mgalz) :: uf_temp, vf_temp, wf_temp
        integer :: jj, kk


        call timer_continue( timer_pkx1 )
        ! ifft in z, only kx >= 0, size (mgalz, nkx_pos)
        call ifftz( uf, uf_t )
        call ifftz( vf, vf_t )
        call ifftz( wf, wf_t )

        ! ifft in z, and fill in kx < 0, size (mgalz, nkx_full,)
        call ifftz_fillnegkx( fxf, fx_t )
        call ifftz_fillnegkx( fyf, fy_t )
        call ifftz_fillnegkx( fzf, fz_t )

        call ifftz_fillnegkx( kx_derivative(uf), dudx_t )
        call ifftz_fillnegkx( kx_derivative(vf), dvdx_t )
        call ifftz_fillnegkx( kx_derivative(wf), dwdx_t )

        call ifftz_fillnegkx( dudyf, dudy_t )
        call ifftz_fillnegkx( dvdyf, dvdy_t )
        call ifftz_fillnegkx( dwdyf, dwdy_t )

        call ifftz_fillnegkx( kz_derivative(uf), dudz_t )
        call ifftz_fillnegkx( kz_derivative(vf), dvdz_t )
        call ifftz_fillnegkx( kz_derivative(wf), dwdz_t )

        call timer_stop( timer_pkx1 )

        call timer_continue( timer_pkx2 )
        ! loop over nkx_pos
        DO jj = 1, nkx_pos, DS_kx
            ! kx1 is non-negative, with size nkx_pos
            ! kx2 is the full range, with size nkx_full
            !
            ! The larger kx1 is, the smaller kx2 has to be to keep kx3 in range.
            ! If kx1 = 0, then kx2, kx3 will both be the full range.
            !
            ! Therefore, we will remove (jj-1) points from the highest kx2,
            ! and also remove (jj-1) points from the lowest kx3

            ! u v w for the current kx1, size (mgalz)
            uf_temp = uf_t(:,jj)
            vf_temp = vf_t(:,jj)
            wf_temp = wf_t(:,jj)


            ! Loop over the active range of kx2, which has (jj-1) points from the highest kx2 removed
            DO kk = 1, nkx_full-(jj-1)
                ! --------------------------- Px_kx ---------------------------
                Px_kx(jj,kk) = SUM( &
                    ! compute f for the current kx1 + kx2
                    (-uf_temp*dudx_t(:,kk) -vf_temp*dudy_t(:,kk) -wf_temp*dudz_t(:,kk)) &
                    ! project onto this kx3, which is shifted by (jj-1) points compared to kx2
                    *CONJG(fx_t(:,kk+jj-1)) &
                    ! average in z and time
                    )/real(mgalz, cp)/real(nt,cp) &
                    ! add to the cumulative time sum
                    + Px_kx(jj,kk)

                ! --------------------------- Py_kx ---------------------------
                Py_kx(jj,kk) = SUM( &
                    ! compute f for the current kx1 + kx2
                    (-uf_temp*dvdx_t(:,kk) -vf_temp*dvdy_t(:,kk) -wf_temp*dvdz_t(:,kk)) &
                    ! project onto this kx3, which is shifted by (jj-1) points compared to kx2
                    *CONJG(fy_t(:,kk+jj-1)) &
                    ! average in z and time
                    )/real(mgalz, cp)/real(nt,cp) &
                    ! add to the cumulative time sum
                    + Py_kx(jj,kk)

                ! --------------------------- Pz_kx ---------------------------
                Pz_kx(jj,kk) = SUM( &
                    ! compute f for the current kx1 + kx2
                    (-uf_temp*dwdx_t(:,kk) -vf_temp*dwdy_t(:,kk) -wf_temp*dwdz_t(:,kk)) &
                    ! project onto this kx3, which is shifted by (jj-1) points compared to kx2
                    *CONJG(fz_t(:,kk+jj-1)) &
                    ! average in z and time
                    )/real(mgalz, cp)/real(nt,cp) &
                    ! add to the cumulative time sum
                    + Pz_kx(jj,kk)
            ENDDO

        ENDDO
        call timer_stop( timer_pkx2 )

    end subroutine compute_Pkx


    ! subroutine compute_Pkz(uf,vf,wf,dudyf,dvdyf,dwdyf,fxf,fyf,fzf, Px_kz, Py_kz, Pz_kz)
    ! This function computes the projection coefficient P_kz at a given yplane
    !
    ! Arguments:
    !   uf,vf,wf:            [double complex, Size (mxf,mzf), Input]
    !                        velocity fields at a single y plane
    !   dudyf,dvdyf,dwdyf:   [double complex, Size (mxf,mzf), Input]
    !                        y derivatives of velocity fields at a single y plane
    !   fxf,fyf,fzf:         [double complex, Size (mxf,mzf), Input]
    !                        mean subtracted non-linear forcing at a single y plane
    !   Px_kz, Py_kz, Pz_kz: [double complex, Size (nkz_pos,nkz_full), Output]
    !                        The computed projection coefficients
    subroutine compute_Pkz(uf,vf,wf,dudyf,dvdyf,dwdyf,fxf,fyf,fzf, Px_kz, Py_kz, Pz_kz)
        use FFT_PRECISION_CONTROL, only: ifftx
        use wavenumbers, only: kx_derivative, kz_derivative
        ! Input/Outputs
        complex(kind=cp), intent( in), dimension(mxf,mzf) :: uf,vf,wf, dudyf,dvdyf,dwdyf, fxf,fyf,fzf
        complex(kind=cp), intent(inout), dimension(nkz_pos,nkz_full) :: Px_kz, Py_kz, Pz_kz

        ! Results after ifft in x
        complex(kind=cp), dimension(mgalx, nkz_pos ) :: uf_t, vf_t, wf_t
        complex(kind=cp), dimension(mgalx, nkz_full) :: fx_t, fy_t, fz_t
        complex(kind=cp), dimension(mgalx, nkz_full) :: dudx_t, dvdx_t, dwdx_t
        complex(kind=cp), dimension(mgalx, nkz_full) :: dudy_t, dvdy_t, dwdy_t
        complex(kind=cp), dimension(mgalx, nkz_full) :: dudz_t, dvdz_t, dwdz_t
        ! Temp variables for the loop
        complex(kind=cp), dimension(mgalx) :: uf_temp, vf_temp, wf_temp
        integer :: jj, kk


        call timer_continue( timer_pkz1 )
        ! ifft in x, only kz >= 0, size (mgalx, nkz_pos )
        call ifftx( uf, uf_t )
        call ifftx( vf, vf_t )
        call ifftx( wf, wf_t )

        ! ifft in x, and fill in kz < 0, size (mgalx, nkz_full)
        call ifftx_fillnegkz( fxf, fx_t )
        call ifftx_fillnegkz( fyf, fy_t )
        call ifftx_fillnegkz( fzf, fz_t )

        call ifftx_fillnegkz( kx_derivative(uf), dudx_t )
        call ifftx_fillnegkz( kx_derivative(vf), dvdx_t )
        call ifftx_fillnegkz( kx_derivative(wf), dwdx_t )

        call ifftx_fillnegkz( dudyf, dudy_t )
        call ifftx_fillnegkz( dvdyf, dvdy_t )
        call ifftx_fillnegkz( dwdyf, dwdy_t )

        call ifftx_fillnegkz( kz_derivative(uf), dudz_t )
        call ifftx_fillnegkz( kz_derivative(vf), dvdz_t )
        call ifftx_fillnegkz( kz_derivative(wf), dwdz_t )

        call timer_stop( timer_pkz1 )

        call timer_continue( timer_pkz2 )
        ! loop over nkz_pos
        DO jj = 1, nkz_pos, DS_kz
            ! kz1 is non-negative, with size nkz_pos
            ! kz2 is the full range, with size nkz_full
            !
            ! The larger kz1 is, the smaller kz2 has to be to keep kz3 in range.
            ! If kz1 = 0, then kz2, kz3 will both be the full range.
            !
            ! Therefore, we will remove (jj-1) points from the highest kz2,
            ! and also remove (jj-1) points from the lowest kz3

            ! u v w for the current kz1, size (mgalx)
            uf_temp = uf_t(:,jj)
            vf_temp = vf_t(:,jj)
            wf_temp = wf_t(:,jj)


            ! Loop over the active range of kz2, which has (jj-1) points from the highest kz2 removed
            DO kk = 1, nkz_full-(jj-1)
                ! --------------------------- Px_kz ---------------------------
                Px_kz(jj,kk) = SUM( &
                    ! compute f for the current kz1 + kz2
                    (-uf_temp*dudx_t(:,kk) -vf_temp*dudy_t(:,kk) -wf_temp*dudz_t(:,kk)) &
                    ! project onto this kz3, which is shifted by (jj-1) points compared to kz2
                    *CONJG(fx_t(:,kk+jj-1)) &
                    ! average in x and time
                    )/real(mgalz, cp)/real(nt,cp) &
                    ! add to the cumulative time sum
                    + Px_kz(jj,kk)

                ! --------------------------- Py_kz ---------------------------
                Py_kz(jj,kk) = SUM( &
                    ! compute f for the current kz1 + kz2
                    (-uf_temp*dvdx_t(:,kk) -vf_temp*dvdy_t(:,kk) -wf_temp*dvdz_t(:,kk)) &
                    ! project onto this kz3, which is shifted by (jj-1) points compared to kz2
                    *CONJG(fy_t(:,kk+jj-1)) &
                    ! average in x and time
                    )/real(mgalz, cp)/real(nt,cp) &
                    ! add to the cumulative time sum
                    + Py_kz(jj,kk)

                ! --------------------------- Pz_kz ---------------------------
                Pz_kz(jj,kk) = SUM( &
                    ! compute f for the current kz1 + kz2
                    (-uf_temp*dwdx_t(:,kk) -vf_temp*dwdy_t(:,kk) -wf_temp*dwdz_t(:,kk)) &
                    ! project onto this kz3, which is shifted by (jj-1) points compared to kz2
                    *CONJG(fz_t(:,kk+jj-1)) &
                    ! average in x and time
                    )/real(mgalz, cp)/real(nt,cp) &
                    ! add to the cumulative time sum
                    + Pz_kz(jj,kk)
            ENDDO

        ENDDO
        call timer_stop( timer_pkz2 )

    end subroutine compute_Pkz


    ! subroutine ifftx_fillnegkz( matrix_in, matrix_out )
    ! Computes the ifft in x only for a 2d matrix (single y location)
    ! and then fill the negative kz parts using hermitian symmetry
    !
    ! Arguments:
    !   matrix_in:  [double complex, size (mxf,mzf), Input]
    !               fourier domain data
    !   matrix_out: [double complex, size (mgalx, nkz_full), Output]
    !               ifft in x only, for the both positive and negative kz wavenumbers
    subroutine ifftx_fillnegkz( matrix_in, matrix_out )
        use FFT_PRECISION_CONTROL, only: ifftx
        complex(kind=cp), intent( in), dimension(mxf,mzf) :: matrix_in
        complex(kind=cp), intent(out), dimension(mgalx, nkz_full) :: matrix_out

        ! kz = 0 is now in the middle, which is position nkz_pos
        call ifftx( matrix_in, matrix_out(:,nkz_pos:nkz_full))
        ! kz < 0 using hermitian symmetry, copy data from the kz > 0 part
        matrix_out(:,1:nkz_pos-1) =  CONJG( matrix_out(:,nkz_full:nkz_pos+1:-1) )
    end subroutine ifftx_fillnegkz


    ! subroutine ifftz_fillnegkx( matrix_in, matrix_out )
    ! Computes the ifft in z only for a 2d matrix (single y location)
    ! and then fill the negative kx parts using hermitian symmetry
    !
    ! Arguments:
    !   matrix_in:  [double complex, size (mxf,mzf), Input]
    !               fourier domain data
    !   matrix_out: [double complex, size (mgalz, nkx_full), Output]
    !               ifft in z only, for the both positive and negative kx wavenumbers
    !               Note that this has a built in transpose that exchanges the
    !               kx dimension with the z dimension for improved indexing
    !               effeciency
    subroutine ifftz_fillnegkx( matrix_in, matrix_out )
        use FFT_PRECISION_CONTROL, only: ifftz
        complex(kind=cp), intent( in), dimension(mxf,mzf) :: matrix_in
        complex(kind=cp), intent(out), dimension(mgalz, nkx_full) :: matrix_out


        ! kx = 0 is now in the middle, which is position nkx_pos
        call ifftz( matrix_in, matrix_out(:,nkx_pos:nkx_full))
        ! kx < 0 using hermitian symmetry, copy data from the kx > 0 part
        matrix_out(:,1:nkx_pos-1) =  CONJG( matrix_out(:,nkx_full:nkx_pos+1:-1) )
    end subroutine ifftz_fillnegkx


end module ComputeTriadicProjection