module ComputeForcing
    use types, only: sp, dp

    implicit none
    private
#   include "parameters"

    public :: compute_forcing

contains
    ! subroutine compute_forcing(uf,vf,wf,dudyf,dvdyf,dwdyf, fxf,fyf,fzf, yplane)
    ! This function computes the non-linear forcing fxf, fyf, fzf at a given y
    ! plan then subtracts the mean (mean subtraction in controlled by the logical flag meansubtract1D_forcing )
    !
    ! Arguments:
    !   uf,vf,wf:          [double/single complex, Size (mxf,mzf), Input]
    !                      velocity fields at a single y plane
    !   dudyf,dvdyf,dwdyf: [double/single complex, Size (mxf,mzf), Input]
    !                      y derivatives of velocity fields at a single y plane
    !   yplane :           [Integer, Input]
    !                      index of this current y plane
    !   fxf,fyf,fzf:       [double/single complex, Size (mxf,mzf), Output]
    !                      (mean subtracted) non-linear forcing at a single y plane
    subroutine compute_forcing(uf,vf,wf,dudyf,dvdyf,dwdyf, yplane, fxf,fyf,fzf)
        use FFT_PRECISION_CONTROL, only: ifft2, fft2
        use wavenumbers, only: kx_derivative, kz_derivative, truncate_kx_2D
        complex(kind=cp), intent( in), dimension(mxf,mzf) :: uf,vf,wf, dudyf,dvdyf,dwdyf
        complex(kind=cp), intent(out), dimension(mxf,mzf) :: fxf,fyf,fzf
        integer, intent(in) :: yplane

        ! Global stats variables
        real(kind=cp) :: U_cl, u00, v00, w00, fx00, fy00, fz00
        common /stats/ U_cl, &                           ! U_cl
                       u00(my), v00(my), w00(my), &      ! kx = kz = 0 modes of velocity
                       fx00(my), fy00(my), fz00(my)      ! kx = kz = 0 modes of forcing
        save /stats/

        ! velocities in physical space
        real(kind=cp), dimension(mgalx,mgalz) :: u, v, w
        ! derivatives in physical space
        real(kind=cp), dimension(mgalx,mgalz) :: ddx, ddy, ddz
        ! f in physical space
        real(kind=cp), dimension(mgalx,mgalz) :: f


        ! ifft2 of uf, vf, wf into physical spac
        call ifft2( uf, u )
        call ifft2( vf, v )
        call ifft2( wf, w )

        ! ----------------------------- Compute fx -----------------------------
        ! Compute the 3 derivatives of u
        call ifft2( kx_derivative(uf), ddx )
        call ifft2( dudyf            , ddy )
        call ifft2( kz_derivative(uf), ddz )
        ! Compute fx in physical space
        f = - u * ddx - v * ddy - w * ddz
        ! fft2 of fx back into the Fourier space
        call fft2( f, fxf )
        ! Mean subtract
        if ( meansubtract1D_forcing .eqv. .true. ) then
            fxf(1,1) = fxf(1,1) - fx00(yplane)
        endif
        ! kx = 0 subtract
        if ( subtract_kx0_forcing .eqv. .true. ) then
            fxf(1,:) = (0.0_cp, 0.0_cp)
        endif
        ! kx truncation
        call truncate_kx_2D( fxf, truncate_kx_forcing )

        ! ----------------------------- Compute fy -----------------------------
        ! Compute the 3 derivatives of v
        call ifft2( kx_derivative(vf), ddx )
        call ifft2( dvdyf            , ddy )
        call ifft2( kz_derivative(vf), ddz )
        ! Compute fy in physical space
        f = - u * ddx - v * ddy - w * ddz
        ! fft2 of fy back into the Fourier space
        call fft2( f, fyf )
        ! Mean subtract
        if ( meansubtract1D_forcing .eqv. .true. ) then
            fyf(1,1) = fyf(1,1) - fy00(yplane)
        endif
        ! kx = 0 subtract
        if ( subtract_kx0_forcing .eqv. .true. ) then
            fyf(1,:) = (0.0_cp, 0.0_cp)
        endif
        ! kx truncation
        call truncate_kx_2D( fyf, truncate_kx_forcing )

        ! ----------------------------- Compute fz -----------------------------
        ! Compute the 3 derivatives of w
        call ifft2( kx_derivative(wf), ddx )
        call ifft2( dwdyf            , ddy )
        call ifft2( kz_derivative(wf), ddz )
        ! Compute fz in physical space
        f = - u * ddx - v * ddy - w * ddz
        ! fft2 of fz back into the Fourier space
        call fft2( f, fzf )
        ! Mean subtract
        if ( meansubtract1D_forcing .eqv. .true. ) then
            fzf(1,1) = fzf(1,1) - fz00(yplane)
        endif
        ! kx = 0 subtract
        if ( subtract_kx0_forcing .eqv. .true. ) then
            fzf(1,:) = (0.0_cp, 0.0_cp)
        endif
        ! kx truncation
        call truncate_kx_2D( fzf, truncate_kx_forcing )

    end subroutine compute_forcing


end module ComputeForcing