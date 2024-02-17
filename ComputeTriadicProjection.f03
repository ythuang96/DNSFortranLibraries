module ComputeTriadicProjection
    use types, only: sp, dp

    implicit none
    private
#   include "parameters"

    public :: compute_Pkx, compute_Pkz, compute_Pom

contains
    ! subroutine compute_Pkx(uf,vf,wf,dudyf,dvdyf,dwdyf,fxf,fyf,fzf, Px_kx, Py_kx, Pz_kx)
    ! This function computes the projection coefficient P_kx at a given yplane
    ! and for ONE SPECIFIC TEMPORAL SNAPSHOT
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
        complex(kind=cp), intent(out), dimension(nkx_pos,nkx_full) :: Px_kx, Py_kx, Pz_kx

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

        ! loop over nkx_pos
        DO jj = 1, nkx_pos
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
                    ! average in z
                    )/real(mgalz, cp)

                ! --------------------------- Py_kx ---------------------------
                Py_kx(jj,kk) = SUM( &
                    ! compute f for the current kx1 + kx2
                    (-uf_temp*dvdx_t(:,kk) -vf_temp*dvdy_t(:,kk) -wf_temp*dvdz_t(:,kk)) &
                    ! project onto this kx3, which is shifted by (jj-1) points compared to kx2
                    *CONJG(fy_t(:,kk+jj-1)) &
                    ! average in z
                    )/real(mgalz, cp)

                ! --------------------------- Pz_kx ---------------------------
                Pz_kx(jj,kk) = SUM( &
                    ! compute f for the current kx1 + kx2
                    (-uf_temp*dwdx_t(:,kk) -vf_temp*dwdy_t(:,kk) -wf_temp*dwdz_t(:,kk)) &
                    ! project onto this kx3, which is shifted by (jj-1) points compared to kx2
                    *CONJG(fz_t(:,kk+jj-1)) &
                    ! average in z
                    )/real(mgalz, cp)
            ENDDO

        ENDDO

    end subroutine compute_Pkx


    ! subroutine compute_Pkz(uf,vf,wf,dudyf,dvdyf,dwdyf,fxf,fyf,fzf, Px_kz, Py_kz, Pz_kz)
    ! This function computes the projection coefficient P_kz at a given yplane
    ! and for ONE SPECIFIC TEMPORAL SNAPSHOT
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
        complex(kind=cp), intent(out), dimension(nkz_pos,nkz_full) :: Px_kz, Py_kz, Pz_kz

        ! Results after ifft in x
        complex(kind=cp), dimension(mgalx, nkz_pos ) :: uf_t, vf_t, wf_t
        complex(kind=cp), dimension(mgalx, nkz_full) :: fx_t, fy_t, fz_t
        complex(kind=cp), dimension(mgalx, nkz_full) :: dudx_t, dvdx_t, dwdx_t
        complex(kind=cp), dimension(mgalx, nkz_full) :: dudy_t, dvdy_t, dwdy_t
        complex(kind=cp), dimension(mgalx, nkz_full) :: dudz_t, dvdz_t, dwdz_t
        ! Temp variables for the loop
        complex(kind=cp), dimension(mgalx) :: uf_temp, vf_temp, wf_temp
        integer :: jj, kk


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

        ! loop over nkz_pos
        DO jj = 1, nkz_pos
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
                    ! average in x
                    )/real(mgalz, cp)

                ! --------------------------- Py_kz ---------------------------
                Py_kz(jj,kk) = SUM( &
                    ! compute f for the current kz1 + kz2
                    (-uf_temp*dvdx_t(:,kk) -vf_temp*dvdy_t(:,kk) -wf_temp*dvdz_t(:,kk)) &
                    ! project onto this kz3, which is shifted by (jj-1) points compared to kz2
                    *CONJG(fy_t(:,kk+jj-1)) &
                    ! average in x
                    )/real(mgalz, cp)

                ! --------------------------- Pz_kz ---------------------------
                Pz_kz(jj,kk) = SUM( &
                    ! compute f for the current kz1 + kz2
                    (-uf_temp*dwdx_t(:,kk) -vf_temp*dwdy_t(:,kk) -wf_temp*dwdz_t(:,kk)) &
                    ! project onto this kz3, which is shifted by (jj-1) points compared to kz2
                    *CONJG(fz_t(:,kk+jj-1)) &
                    ! average in x
                    )/real(mgalz, cp)
            ENDDO

        ENDDO

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


    ! subroutine compute_Pom(u_xom,v_xom,w_xom, ddx_xom,ddy_xom,ddz_xom, f_xom, P_om)
    ! This function computes the projection coefficient P_om at a given yplane
    !
    ! Arguments:
    !   u_xom, v_xom, w_xom      : [double/single complex, Size (mgalx, b_z:e_z, nom_pos), Input]
    !                              velocity fields at a single y plane, in physical x, z space and Fourier omega space, distributed in z
    !   ddx_xom, ddy_xom, ddz_xom: [double/single complex, Size (mgalx, b_z:e_z, nom_pos), Input]
    !                              x, y, z derivatives
    !   f_xom                    : [double/single complex, Size (mgalx, b_z:e_z, nom_pos), Input]
    !                              non-linear forcing
    !   P_om                     : [double/single complex, Size (nom_pos,nom), Output]
    !                              The computed projection coefficients
    !
    ! Note:
    !   To compute Px, input the u derivatives dudx, dudy, dudz, and fx
    !   Similarly, to compute Py, input the v derivatives and fy
    !   to compute Pz, input the w derivatives and fz
    !
    !   All input data are in physical x, z space and Fourier omega space, distributed in z
    !   And with all necessary mean subtraction
    subroutine compute_Pom(u_xom,v_xom,w_xom, ddx_xom,ddy_xom,ddz_xom, f_xom, P_om)
        use mpi
        use parallelization, only: allreduce_C2
        ! Input/Outputs
        complex(kind=cp), intent( in), dimension(:,:,:) :: u_xom,v_xom,w_xom, ddx_xom,ddy_xom,ddz_xom, f_xom
        complex(kind=cp), intent(out), dimension(nom_pos,nom) :: P_om

        ! Variables with full omega range
        complex(kind=cp), dimension(:,:,:), allocatable :: ddx_fullom, ddy_fullom, ddz_fullom, f_fullom

        ! Matrix dimensions
        integer, dimension(3) :: data_dim
        ! Temp variables for the loop
        complex(kind=cp), dimension(:,:), allocatable :: u_temp, v_temp, w_temp
        integer :: jj, kk

        ! ------------------------- Allocate Vaiables -------------------------
        ! Input matrix dimensions
        data_dim = shape(u_xom)
        ! Allocate temp variables with full omega range
        allocate( ddx_fullom(data_dim(1), data_dim(2), nom) )
        allocate( ddy_fullom(data_dim(1), data_dim(2), nom) )
        allocate( ddz_fullom(data_dim(1), data_dim(2), nom) )
        allocate(   f_fullom(data_dim(1), data_dim(2), nom) )
        ! Allocate temp variables for u, v, w at a single omega for the loop
        allocate( u_temp(data_dim(1), data_dim(2)) )
        allocate( v_temp(data_dim(1), data_dim(2)) )
        allocate( w_temp(data_dim(1), data_dim(2)) )

        ! Zero results
        P_om = (0.0_cp, 0.0_cp)

        ! ----------------------- Fillin negative omega -----------------------
        ! Fill the other half of omega using hermitian symmetry
        call fillnegomega( ddx_xom, ddx_fullom )
        call fillnegomega( ddy_xom, ddy_fullom )
        call fillnegomega( ddz_xom, ddz_fullom )
        call fillnegomega(   f_xom,   f_fullom )

        ! ------------------------- Loop over omega_1 -------------------------
        DO jj = 1, nom_pos
            ! om1 is non-negative omega, with size nom_pos
            ! om2 is the full range omega, with size nom
            !
            ! The larger om1 is, the smaller om2 has to be to keep om3 in range.
            ! If om1 = 0, then om2, om3 will both be the full range.
            !
            ! Therefore, we will remove (jj-1) points from the highest om2,
            ! and also remove (jj-1) points from the lowest om3

            ! u v w for the current om1, size (mgalx, b_z:e_z)
            u_temp = u_xom(:,:,jj)
            v_temp = v_xom(:,:,jj)
            w_temp = w_xom(:,:,jj)

            ! Loop over the active range of om2, which has (jj-1) points from the highest kx2 removed
            DO kk = 1, nom-(jj-1)
                P_om(jj,kk) = SUM( &
                    ! compute f for the current om1 + om2
                    (-u_temp*ddx_fullom(:,:,kk) -v_temp*ddy_fullom(:,:,kk) -w_temp*ddz_fullom(:,:,kk)) &
                    ! project onto this om3, which is shifted by (jj-1) points compared to om2
                    *CONJG(f_fullom(:,:,kk+jj-1)) &
                    ! average in x and z
                    )/real(mgalx, cp)/real(mgalz,cp)
            ENDDO
        ENDDO

        ! ------------------- Sum result over all processors -------------------
        ! data is distributed in z, use all reduce to sum in z
        call allreduce_C2( P_om, MPI_SUM )

        ! ------------------------------ Clean Up ------------------------------
        deallocate( ddx_fullom )
        deallocate( ddy_fullom )
        deallocate( ddz_fullom )
        deallocate(   f_fullom )
        deallocate(   u_temp   )
        deallocate(   v_temp   )
        deallocate(   w_temp   )

    end subroutine compute_Pom


    ! subroutine fillnegomega(matrix_in, matrix_out)
    ! This subroutine fills the negative omega parts for a matrix
    !
    ! Arguements
    !   matrix_in : [complex matrix size (:,:,nom_pos), Input ]
    !               matrix with omega as dimension 3, only the non-negative omegas
    !   matrix_out: [complex matrix size (:,:,nom    ), Output]
    !               matrix with both postivie and negative omega
    !               dimension 3 corresponds to omegas that are montonically increasing
    !
    ! Data Arrangement
    !   dimension 3 of input:
    !     matrix index         1     ...   nom_pos
    !     corresponding om:  0(DC)   ...   om_max
    !   dimension 3 of output:
    !     matrix index          1      ...   (nom+1)/2 = nom_pos   ...     nom
    !     corresponding om:  -om_max   ...          0(DC)          ...   om_max
    subroutine fillnegomega(matrix_in, matrix_out)
        ! Input/Outputs
        complex(kind=cp), intent( in), dimension(:,:,:) :: matrix_in
        complex(kind=cp), intent(out), dimension(:,:,:) :: matrix_out

        ! Loop index
        integer :: kk

        ! negative omega range, take the conjugate and flip
        DO kk = 1, nom_pos-1
            matrix_out(:,:,kk) = conjg( matrix_in(:,:,nom_pos+1-kk))
        ENDDO
        ! omega >= 0 component (just shifted)
        matrix_out(:,:,nom_pos:nom) = matrix_in
    end subroutine fillnegomega

end module ComputeTriadicProjection