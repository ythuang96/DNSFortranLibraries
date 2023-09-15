module resolvent
    use types, only: dp
    use wallnormalgrid, only: d_dy_complexvector
    implicit none
    private
#   include "parameters"

    ! 1D Mean and its first and second derivatives
    real(kind=dp), dimension(myf) :: U, dUdy, dUdy2
    ! Reynolds number
    real(kind=dp) :: ReyNum
    ! Identity matrix
    real(kind=dp), dimension(myf,myf) :: IdentityMat
    ! complex number
    complex(kind=dp) :: complex_i
    ! First and second derivative matrices
    real(kind=dp), dimension(myf,myf) :: D1, D2

    ! Temporary vectors
    complex(kind=dp),dimension(myf) :: temp1, temp2, temp3, temp4

    public :: initialize_resolvent_module
    public :: resolvent_A, resolvent_B
    public :: resolvent_LAP, resolvent_LOS, resolvent_LSQ
    public :: resolvent_H, resolvent_H00
    public :: vec_X_A1HB, vec_X_A2HB, vec_X_A3HB, AHB_X_vec

contains
! *************************** Initialize This Moudle ***************************
    ! initialize this module
    ! Arguments
    !   U_input, dUdy_input, dUdy2_input [real, size myf, Input] 1D streamwise mean and its 1st and 2nd derivative
    !   ReyNum_input                     [real, Input] Reynolds number (based on Ucl and h, as everything else)
    !   D1D2File                         [string] h5 file name for reading the 1st and 2nd derivative matrices
    !   myid                             [integer] MPI processor ID
    subroutine initialize_resolvent_module( U_input, dUdy_input, dUdy2_input, ReyNum_input, D1D2File, myid)
        use mpi
        use h5load, only: h5load_R2
        real(kind=dp), dimension(myf), intent(in) :: U_input, dUdy_input, dUdy2_input
        real(kind=dp), intent(in) :: ReyNum_input
        character(len=*), intent(in) :: D1D2File
        integer, intent(in) :: myid

        integer :: ii, ierr

        ! set as global variables for this module
        U      = U_input
        dUdy   = dUdy_input
        dUdy2  = dUdy2_input
        ReyNum = ReyNum_input

        ! complex number
        complex_i = ( 0.0_dp, 1.0_dp )
        ! Identity matrix
        IdentityMat = 0.0_dp
        DO ii = 1, myf
            IdentityMat(ii,ii) = 1.0_dp
        END DO

        ! Master read D1 and D2 matrices from file
        if (myid .eq. 0) then
            D1 = h5load_R2( D1D2File, "D1" )
            D2 = h5load_R2( D1D2File, "D2" )
        endif
        ! Broadcast D1 and D2 to all slaves
        call MPI_BCAST(D1, myf*myf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(D2, myf*myf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    end subroutine initialize_resolvent_module

! **************************** Compute Sub Matrices ****************************
    ! Compute the A matrix (v, o2) -> (u, v, w)
    ! Arguments
    !   kx, kz [real, Input] kx and kz wavenumbers
    ! Return
    !   A      [complex, 3Ny X 2Ny, output] the A matrix, with the special case of kx = kz = 0 taken care of
    function resolvent_A(kx, kz) result(A)
        real(kind=dp), intent(in) :: kx, kz
        complex(kind=dp), dimension(3*myf, 2*myf) :: A

        real(kind=dp) :: k2
        k2 = kx**2 + kz**2

        ! Initialize all elements to 0
        A = ( 0.0_dp, 0.0_dp )

        if ((kx == 0) .and. (kz == 0)) then ! special case for kx = kz = 0
            A(      1:  myf,     1:  myf) = IdentityMat
            A(2*myf+1:3*myf, myf+1:2*myf) = IdentityMat
        else
            A(      1:  myf,     1:  myf) = + (complex_i*kx/k2)*D1
            A(      1:  myf, myf+1:2*myf) = - (complex_i*kz/k2)*IdentityMat
            A(  myf+1:2*myf,     1:  myf) = IdentityMat
            A(2*myf+1:3*myf,     1:  myf) = + (complex_i*kz/k2)*D1
            A(2*myf+1:3*myf, myf+1:2*myf) = + (complex_i*kx/k2)*IdentityMat
        endif
    end function resolvent_A

    ! Compute the B matrix (fx, fy, fz) -> (fv, feta)
    ! Arguments
    !   kx, kz [real, Input] kx and kz wavenumbers
    ! Return
    !   B      [complex, 2Ny X 3Ny, output] the B matrix, with the proper boundary conditions applied
    !                                       and the special case of kx = kz = 0 taken care of
    function resolvent_B(kx,kz) result(B)
        real(kind=dp), intent(in) :: kx, kz
        complex(kind=dp), dimension(2*myf, 3*myf) :: B

        real(kind=dp) :: k2
        k2 = kx**2 + kz**2

        ! Initialize all elements to 0
        B = ( 0.0_dp, 0.0_dp )

        if ((kx == 0) .and. (kz == 0)) then ! special case for kx = kz = 0
            B(    1:  myf,       1:  myf) = IdentityMat
            B(myf+1:2*myf, 2*myf+1:3*myf) = IdentityMat

            ! Apply boundary conditions
            B(    1,:) = (0.0_dp, 0.0_dp) !   u( 1 ) = 0
            B(myf  ,:) = (0.0_dp, 0.0_dp) !   u(myf) = 0
            B(myf+1,:) = (0.0_dp, 0.0_dp) !   w( 1 ) = 0
            B(myf*2,:) = (0.0_dp, 0.0_dp) !   w(myf) = 0
        else
            B(    1:  myf,       1:  myf) = - (complex_i*kx)*D1
            B(    1:  myf,   myf+1:2*myf) = - (k2)*IdentityMat
            B(    1:  myf, 2*myf+1:3*myf) = - (complex_i*kz)*D1
            B(myf+1:2*myf,       1:  myf) = + (complex_i*kz)*IdentityMat
            B(myf+1:2*myf, 2*myf+1:3*myf) = - (complex_i*kx)*IdentityMat

            ! Apply boundary conditions
            B(    1,:) = (0.0_dp, 0.0_dp) !    v( 1 ) = 0
            B(myf  ,:) = (0.0_dp, 0.0_dp) !    v(myf) = 0
            B(    2,:) = (0.0_dp, 0.0_dp) ! dvdy( 1 ) = 0
            B(myf-1,:) = (0.0_dp, 0.0_dp) ! dvdy(myf) = 0
            B(myf+1,:) = (0.0_dp, 0.0_dp) !   o2( 1 ) = 0
            B(myf*2,:) = (0.0_dp, 0.0_dp) !   o2(myf) = 0
        endif
    end function resolvent_B

    ! Compute the Laplacian ( D2 - k^2 )
    ! Arguments
    !   kx, kz [real, Input] kx and kz wavenumbers
    ! Return
    !   LAP    [real, Ny X Ny, output] the Laplacian matrix, not used for the speical case of kx = kz = 0
    function resolvent_LAP(kx,kz) result(LAP)
        real(kind=dp), intent(in) :: kx, kz
        real(kind=dp), dimension(myf, myf) :: LAP
        integer :: ii

        real(kind=dp) :: k2
        k2 = kx**2 + kz**2

        LAP = D2
        DO ii = 1, myf
            LAP(ii,ii) = LAP(ii,ii) - k2
        END DO
    end function resolvent_LAP

    ! Compute the OS operator ( i kx U'' - i kx U LAP + LAP^2/Re )
    ! Arguments
    !   kx  [real, Input] kx wavenumber
    !   LAP [real, Ny X Ny, Input] Laplacian matrix
    ! Return
    !   LOS [complex, Ny X Ny, output] the Orr-Sommerfeld matrix, not used for the speical case of kx = kz = 0
    function resolvent_LOS(kx,LAP) result(LOS)
        real(kind=dp), intent(in) :: kx
        real(kind=dp), intent(in), dimension(myf,myf) :: LAP
        complex(kind=dp), dimension(myf, myf) :: LOS
        integer :: ii

        ! LAP^2 / ReyNum
        LOS = MATMUL(LAP, LAP)/ReyNum
        ! Subtract i kx U LAP (note that U is a diagonal matrix)
        DO ii = 1, myf
            LOS(ii,:) = LOS(ii,:) - (complex_i*kx*U(ii))*LAP(ii,:)
        END DO
        ! Add i kx U'' (only changes the diagonal terms)
        DO ii = 1, myf
            LOS(ii,ii) = LOS(ii,ii) + complex_i*kx*dUdy2(ii)
        END DO

    end function resolvent_LOS

    ! Compute the SQ operator ( i kx U- LAP/Re )
    ! Arguments
    !   kx  [real, Input] kx wavenumber
    !   LAP [real, Ny X Ny, Input] Laplacian matrix
    ! Return
    !   LSQ [complex, Ny X Ny, output] the Squire matrix, not used for the speical case of kx = kz = 0
    function resolvent_LSQ(kx,LAP) result(LSQ)
        real(kind=dp), intent(in) :: kx
        real(kind=dp), intent(in), dimension(myf,myf) :: LAP
        complex(kind=dp), dimension(myf, myf) :: LSQ
        integer :: ii

        ! - LAP / ReyNum
        LSQ = - LAP/ReyNum
        ! Add i kx U (only changes the diagonal terms)
        DO ii = 1, myf
            LSQ(ii,ii) = LSQ(ii,ii) + complex_i*kx*U(ii)
        END DO
    end function resolvent_LSQ

! ***************************** Compute Resolvent *****************************
    ! Compute the Resolvent operator
    ! Arguments
    !   kz            [real, Input] kz wavenumber
    !   omega         [real, Input] temporal frequency
    !   LAP           [real, Ny X Ny, Input] Laplacian matrix
    !   LOS, LSQ      [complex Ny X Ny, Input] Orr-Sommerfeld and Squirre operators
    !   Hvv, Hev, Hee [complex Ny X Ny, Output] three blocks of resolvent matrix with appropriate boundary conditions
    subroutine resolvent_H( kz, omega, LAP, LOS, LSQ, Hvv, Hev, Hee)
        use matrixinv, only: MatrixInverse
        real(kind=dp), intent(in) :: kz, omega
        real(kind=dp), intent(in), dimension(myf,myf) :: LAP
        complex(kind=dp), intent( in), dimension(myf,myf) :: LOS, LSQ
        complex(kind=dp), intent(out), dimension(myf,myf) :: Hvv, Hev, Hee
        complex(kind=dp), dimension(myf,myf) :: temp

        integer :: ii

        ! ----------------- Build Hvv and Hee before inversion -----------------
        Hvv = (-complex_i*omega)*LAP - LOS

        Hee = LSQ
        DO ii = 1,myf
            Hee(ii,ii) = Hee(ii,ii) - complex_i*omega ! - i omega + LSQ
        END DO

        ! --------------------- Apply boundary conditions ---------------------
        ! v(1) = 0
        Hvv(  1,2:myf  ) = (0.0_dp, 0.0_dp)
        Hvv(  1,1      ) = (1.0_dp, 0.0_dp)
        ! v(myf) = 0
        Hvv(myf,1:myf-1) = (0.0_dp, 0.0_dp)
        Hvv(myf,  myf  ) = (1.0_dp, 0.0_dp)

        ! dvdy(1) = 0
        Hvv(    2,:) = D1(  1,:)
        ! dvdy(myf) = 0
        Hvv(myf-1,:) = D1(myf,:)

        ! o2(1) = 0
        Hee(  1,2:myf  ) = (0.0_dp, 0.0_dp)
        Hee(  1,1      ) = (1.0_dp, 0.0_dp)
        ! o2(myf) = 0
        Hee(myf,1:myf-1) = (0.0_dp, 0.0_dp)
        Hee(myf,  myf  ) = (1.0_dp, 0.0_dp)
        ! Note, since v(1,myf) = 0, the off-diagonal term does not effect the o2 boundary conditions

        ! -------------------------- Compute inverse --------------------------
        call MatrixInverse( Hvv )
        call MatrixInverse( Hee )

        ! --------------------- Compute off-diagonal term ---------------------
        DO ii = 1,myf
            temp(ii,:) = ( - complex_i*kz*dUdy(ii) )*Hvv(ii,:) ! -i kz U' Hvv (where U' is a diagonal matrix)
        ENDDO
        Hev = matmul(Hee, temp) ! Hee X (-i kz U' Hvv), matmul is matrix multiplication

    end subroutine resolvent_H


    ! Compute the Resolvent operator for the special case of kx = kz = 0
    ! Arguments
    !   omega    [real, Input] temporal frequency
    !   Hvv, Hee [complex Ny X Ny, Output] two diagonal blocks of resolvent matrix with appropriate boundary conditions
    !   Hev      [complex Ny X Ny, Output] zero matrix to keep the code consistent with resolvent_H
    subroutine resolvent_H00( omega, Hvv, Hev, Hee)
        use matrixinv, only: MatrixInverse
        real(kind=dp), intent(in) :: omega
        complex(kind=dp), intent(out), dimension(myf,myf) :: Hvv, Hev, Hee

        ! --------------------- Build Hvv before inversion ---------------------
        Hvv = (-complex_i*omega*IdentityMat) - D2/ReyNum

        ! --------------------- Apply boundary conditions ---------------------
        ! u(1) = 0
        Hvv(  1,2:myf  ) = (0.0_dp, 0.0_dp)
        Hvv(  1,1      ) = (1.0_dp, 0.0_dp)
        ! u(myf) = 0
        Hvv(myf,1:myf-1) = (0.0_dp, 0.0_dp)
        Hvv(myf,  myf  ) = (1.0_dp, 0.0_dp)

        ! -------------------------- Compute inverse --------------------------
        call MatrixInverse( Hvv )

        Hee = Hvv

        Hev = (0.0_dp, 0.0_dp)

    end subroutine resolvent_H00


! *********************** Multiplication with Resolvent ***********************
! ******** (OLD Unused code, replaced by more effecient combined code) ********
    ! This function computes fve = ve X H (left multiplication of H by vector ve)
    ! Arguments
    !   ve            [complex, size 2*Ny, Input] velocity vorticity vector
    !   Hvv, Hev, Hee [complex, Ny X Ny, Input] Resolvent sub blocks
    ! Return
    !   fve           [complex, size 2*Ny, Output] fve = ve X H (left multiplication of H by vector ve)
    ! function ve_X_H( ve, Hvv, Hev, Hee ) result( fve )
    !     complex(kind=dp), dimension(2*myf), intent(in) :: ve
    !     complex(kind=dp), dimension(myf, myf), intent(in) :: Hvv, Hev, Hee
    !     complex(kind=dp), dimension(2*myf) :: fve

    !     ! Note: matmul does not need to transpose the vector
    !     fve(    1:  myf) = matmul( ve(    1:  myf), Hvv ) + matmul( ve(myf+1:2*myf), Hev )
    !     fve(myf+1:2*myf) = matmul( ve(myf+1:2*myf), Hee )
    ! end function


    ! This function computes ve = H X fve (right multiplication of H by vector fve)
    ! Arguments
    !   fve           [complex, size 2*Ny, Input] velocity vorticity forcing vector
    !   Hvv, Hev, Hee [complex, Ny X Ny, Input] Resolvent sub blocks
    ! Return
    !   ve            [complex, size 2*Ny, Output] ve = H X fve (right multiplication of H by vector fve)
    ! function H_X_fve( Hvv, Hev, Hee, fve ) result( ve )
    !     complex(kind=dp), dimension(2*myf), intent(in) :: fve
    !     complex(kind=dp), dimension(myf, myf), intent(in) :: Hvv, Hev, Hee
    !     complex(kind=dp), dimension(2*myf) :: ve

    !     ve(    1:  myf) = matmul( Hvv, fve(1:myf) ) ! v
    !     ve(myf+1:2*myf) = matmul( Hev, fve(1:myf) ) + matmul( Hee, fve(myf+1:2*myf) ) ! o2
    ! end function

! ************************ Multiplication with A and B ************************
! ******** (OLD Unused code, replaced by more effecient combined code) ********
    ! LEFT multiplication with the first Ny rows of A
    ! the special case of kx = kz = 0 is taken care of
    ! Arguments
    !   vec    [complex, size(myf), input] left multiplication vector
    !   kx, kz [real, input] kx and kz wavenumbers
    ! Return
    !   vecout [complex, size(2*myf)] result vector
    ! function vec_X_A1(vec, kx, kz) result(vecout)
    !     real(kind=dp), intent(in) :: kx, kz
    !     complex(kind=dp), intent(in), dimension(myf) :: vec
    !     complex(kind=dp), dimension(2*myf) :: vecout

    !     real(kind=dp) :: k2
    !     k2 = kx**2 + kz**2

    !     if ((kx == 0) .and. (kz == 0)) then ! special case for kx = kz = 0
    !         vecout(    1:  myf) = vec
    !         vecout(myf+1:2*myf) = (0.0_dp, 0.0_dp)
    !     else
    !         vecout(    1:  myf) = + (complex_i*kx/k2)*matmul( vec, D1 ) ! note this is left multiplication
    !         vecout(myf+1:2*myf) = - (complex_i*kz/k2)*vec
    !     endif
    ! end function vec_X_A1

    ! LEFT multiplication with the second Ny rows of A
    ! the special case of kx = kz = 0 is taken care of
    ! Arguments
    !   vec    [complex, size(myf), input] left multiplication vector
    !   kx, kz [real, input] kx and kz wavenumbers
    ! Return
    !   vecout [complex, size(2*myf)] result vector
    ! function vec_X_A2(vec, kx, kz) result(vecout)
    !     real(kind=dp), intent(in) :: kx, kz
    !     complex(kind=dp), intent(in), dimension(myf) :: vec
    !     complex(kind=dp), dimension(2*myf) :: vecout

    !     real(kind=dp) :: k2
    !     k2 = kx**2 + kz**2

    !     if ((kx == 0) .and. (kz == 0)) then ! special case for kx = kz = 0
    !         vecout = (0.0_dp, 0.0_dp)
    !     else
    !         vecout(    1:  myf) = vec
    !         vecout(myf+1:2*myf) = (0.0_dp, 0.0_dp)
    !     endif
    ! end function vec_X_A2

    ! LEFT multiplication with the third Ny rows of A
    ! the special case of kx = kz = 0 is taken care of
    ! Arguments
    !   vec    [complex, size(myf), input] left multiplication vector
    !   kx, kz [real, input] kx and kz wavenumbers
    ! Return
    !   vecout [complex, size(2*myf)] result vector
    ! function vec_X_A3(vec, kx, kz) result(vecout)
    !     real(kind=dp), intent(in) :: kx, kz
    !     complex(kind=dp), intent(in), dimension(myf) :: vec
    !     complex(kind=dp), dimension(2*myf) :: vecout

    !     real(kind=dp) :: k2
    !     k2 = kx**2 + kz**2

    !     if ((kx == 0) .and. (kz == 0)) then ! special case for kx = kz = 0
    !         vecout(    1:  myf) = (0.0_dp, 0.0_dp)
    !         vecout(myf+1:2*myf) = vec
    !     else
    !         vecout(    1:  myf) = + (complex_i*kz/k2)*matmul( vec, D1 ) ! note this is left multiplication
    !         vecout(myf+1:2*myf) = + (complex_i*kx/k2)*vec
    !     endif
    ! end function vec_X_A3


    ! RIGHT multiplication with the entire of A
    ! the special case of kx = kz = 0 is taken care of
    ! Arguments
    !   vec    [complex, size(2*myf), input] right multiplication vector
    !   kx, kz [real, input] kx and kz wavenumbers
    ! Return
    !   vecout [complex, size(3*myf)] result vector
    ! function A_X_vec(kx, kz, vec) result(vecout)
    !     real(kind=dp), intent(in) :: kx, kz
    !     complex(kind=dp), intent(in), dimension(2*myf) :: vec
    !     complex(kind=dp), dimension(3*myf) :: vecout
    !     complex(kind=dp), dimension(myf) :: dvecdy

    !     real(kind=dp) :: k2
    !     k2 = kx**2 + kz**2

    !     if ((kx == 0) .and. (kz == 0)) then ! special case for kx = kz = 0
    !         vecout(      1:  myf) = vec(    1:  myf)
    !         vecout(  myf+1:2*myf) = (0.0_dp, 0.0_dp)
    !         vecout(2*myf+1:3*myf) = vec(myf+1:2*myf)
    !     else
    !         ! y derivative of vec
    !         dvecdy = d_dy_complexvector( vec(1:myf) )
    !         vecout(      1:  myf) = + (complex_i*kx/k2)*dvecdy - (complex_i*kz/k2)*vec(myf+1:2*myf)
    !         vecout(  myf+1:2*myf) = vec(    1:  myf)
    !         vecout(2*myf+1:3*myf) = + (complex_i*kz/k2)*dvecdy + (complex_i*kx/k2)*vec(myf+1:2*myf)
    !     endif
    ! end function A_X_vec


    ! RIGHT multiplication with the entire of B
    ! the special case of kx = kz = 0 is taken care of
    ! the boundary conditions are also taken care of
    ! Arguments
    !   vec    [complex, size(3*myf), input] right multiplication vector
    !   kx, kz [real, input] kx and kz wavenumbers
    ! Return
    !   vecout [complex, size(2*myf)] result vector
    ! function B_X_vec(kx, kz, vec) result(vecout)
    !     real(kind=dp), intent(in) :: kx, kz
    !     complex(kind=dp), intent(in), dimension(3*myf) :: vec
    !     complex(kind=dp), dimension(2*myf) :: vecout
    !     complex(kind=dp), dimension(myf) :: temp

    !     real(kind=dp) :: k2
    !     k2 = kx**2 + kz**2

    !     if ((kx == 0) .and. (kz == 0)) then ! special case for kx = kz = 0
    !         vecout(      1:  myf) = vec(      1:  myf)
    !         vecout(  myf+1:2*myf) = vec(2*myf+1:3*myf)
    !         ! Apply Boundary conditions
    !         vecout(    1) = (0.0_dp, 0.0_dp) !   u( 1 ) = 0
    !         vecout(myf  ) = (0.0_dp, 0.0_dp) !   u(myf) = 0
    !         vecout(myf+1) = (0.0_dp, 0.0_dp) !   w( 1 ) = 0
    !         vecout(myf*2) = (0.0_dp, 0.0_dp) !   w(myf) = 0
    !     else
    !         ! y derivative of vec
    !         temp                = - (complex_i*kx)*vec(      1:  myf) - (complex_i*kz)*vec(2*myf+1:3*myf)
    !         vecout(    1:  myf) =   d_dy_complexvector(temp)          -            k2 *vec(  myf+1:2*myf)
    !         vecout(myf+1:2*myf) = + (complex_i*kz)*vec(      1:  myf) - (complex_i*kx)*vec(2*myf+1:3*myf)
    !         ! Apply Boundary conditions
    !         vecout(    1) = (0.0_dp, 0.0_dp) !    v( 1 ) = 0
    !         vecout(myf  ) = (0.0_dp, 0.0_dp) !    v(myf) = 0
    !         vecout(    2) = (0.0_dp, 0.0_dp) ! dvdy( 1 ) = 0
    !         vecout(myf-1) = (0.0_dp, 0.0_dp) ! dvdy(myf) = 0
    !         vecout(myf+1) = (0.0_dp, 0.0_dp) !   o2( 1 ) = 0
    !         vecout(myf*2) = (0.0_dp, 0.0_dp) !   o2(myf) = 0
    !     endif
    ! end function B_X_vec

    ! LEFT multiplication with the entire of B
    ! the special case of kx = kz = 0 is taken care of
    ! the boundary conditions are also taken care of
    ! Arguments
    !   vec    [complex, size(2*myf), input] left multiplication vector
    !   kx, kz [real, input] kx and kz wavenumbers
    ! Return
    !   vecout [complex, size(3*myf)] result vector
    ! function vec_X_B(vec, kx, kz) result(vecout)
    !     real(kind=dp), intent(in) :: kx, kz
    !     complex(kind=dp), intent(in), dimension(2*myf) :: vec
    !     complex(kind=dp), dimension(3*myf) :: vecout
    !     complex(kind=dp), dimension(myf) :: temp
    !     complex(kind=dp), dimension(2*myf) :: tempvec

    !     real(kind=dp) :: k2
    !     k2 = kx**2 + kz**2

    !     if ((kx == 0) .and. (kz == 0)) then ! special case for kx = kz = 0
    !         vecout(      1:  myf) = vec(    1:  myf)
    !         vecout(  myf+1:2*myf) = (0.0_dp, 0.0_dp)
    !         vecout(2*myf+1:3*myf) = vec(myf+1:2*myf)
    !         ! Apply Boundary conditions
    !         ! Since this is left multiplication, the boundary conditions should be applied to the input vector
    !         ! but for this special case, B are made up of identity matrices, so it does not matter
    !         vecout(      1) = (0.0_dp, 0.0_dp) !   u( 1 ) = 0
    !         vecout(  myf  ) = (0.0_dp, 0.0_dp) !   u(myf) = 0
    !         vecout(2*myf+1) = (0.0_dp, 0.0_dp) !   w( 1 ) = 0
    !         vecout(3*myf  ) = (0.0_dp, 0.0_dp) !   w(myf) = 0
    !     else
    !         tempvec = vec
    !         ! Apply Boundary conditions
    !         ! Since this is left multiplication, the boundary conditions should be applied to the input vector
    !         tempvec(    1) = (0.0_dp, 0.0_dp) !    v( 1 ) = 0
    !         tempvec(myf  ) = (0.0_dp, 0.0_dp) !    v(myf) = 0
    !         tempvec(    2) = (0.0_dp, 0.0_dp) ! dvdy( 1 ) = 0
    !         tempvec(myf-1) = (0.0_dp, 0.0_dp) ! dvdy(myf) = 0
    !         tempvec(myf+1) = (0.0_dp, 0.0_dp) !   o2( 1 ) = 0
    !         tempvec(myf*2) = (0.0_dp, 0.0_dp) !   o2(myf) = 0

    !         ! left multiplication with the derivative matrix
    !         temp =  matmul( tempvec(      1:  myf), D1 ) ! note this is a left multiplication

    !         vecout(      1:  myf) = - (complex_i*kx)*temp + (complex_i*kz)*tempvec(myf+1:2*myf)
    !         vecout(  myf+1:2*myf) =                       -            k2 *tempvec(    1:  myf)
    !         vecout(2*myf+1:3*myf) = - (complex_i*kz)*temp - (complex_i*kx)*tempvec(myf+1:2*myf)
    !     endif
    ! end function vec_X_B


! ************* Combination of all Left and Right Multiplications *************
    function vec_X_A1HB(vec, kx, kz, Hvv, Hev, Hee) result(vecout)
        real(kind=dp), intent(in) :: kx, kz
        complex(kind=dp), intent(in), dimension(myf) :: vec
        complex(kind=dp), dimension(myf, myf), intent(in) :: Hvv, Hev, Hee
        complex(kind=dp), dimension(3*myf) :: vecout

        real(kind=dp) :: k2
        k2 = kx**2 + kz**2

        if ((kx == 0) .and. (kz == 0)) then ! special case for kx = kz = 0
            ! vec_X_A1HB
            vecout(      1:  myf) = matmul( vec, Hvv )
            vecout(  myf+1:3*myf) = (0.0_dp, 0.0_dp)
            ! Apply Boundary conditions
            ! Since this is left multiplication, the boundary conditions should be applied to the input vector
            ! but for this special case, B are made up of identity matrices, so it does not matter
            vecout(      1) = (0.0_dp, 0.0_dp) !   u( 1 ) = 0
            vecout(  myf  ) = (0.0_dp, 0.0_dp) !   u(myf) = 0
        else
            ! vec_X_A1
            temp1 = + (complex_i*kx/k2)*matmul( vec, D1 ) ! note this is left multiplication
            temp2 = - (complex_i*kz/k2)*vec
            ! ve_X_H
            temp3 = matmul( temp1, Hvv ) + matmul( temp2, Hev )
            temp4 = matmul( temp2, Hee )

            ! vec_X_B
            ! Apply Boundary conditions
            ! Since this is left multiplication, the boundary conditions should be applied to the input vector
            temp3(    1) = (0.0_dp, 0.0_dp) !    v( 1 ) = 0
            temp3(myf  ) = (0.0_dp, 0.0_dp) !    v(myf) = 0
            temp3(    2) = (0.0_dp, 0.0_dp) ! dvdy( 1 ) = 0
            temp3(myf-1) = (0.0_dp, 0.0_dp) ! dvdy(myf) = 0
            temp4(    1) = (0.0_dp, 0.0_dp) !   o2( 1 ) = 0
            temp4(myf  ) = (0.0_dp, 0.0_dp) !   o2(myf) = 0

            ! left multiplication with the derivative matrix
            temp1 =  matmul( temp3, D1 ) ! note this is a left multiplication

            vecout(      1:  myf) = - (complex_i*kx)*temp1 + (complex_i*kz)*temp4
            vecout(  myf+1:2*myf) =                        -            k2 *temp3
            vecout(2*myf+1:3*myf) = - (complex_i*kz)*temp1 - (complex_i*kx)*temp4
        endif
    end function vec_X_A1HB

    function vec_X_A2HB(vec, kx, kz, Hvv, Hev, Hee) result(vecout)
        real(kind=dp), intent(in) :: kx, kz
        complex(kind=dp), intent(in), dimension(myf) :: vec
        complex(kind=dp), dimension(myf, myf), intent(in) :: Hvv, Hev, Hee
        complex(kind=dp), dimension(3*myf) :: vecout

        real(kind=dp) :: k2
        k2 = kx**2 + kz**2

        if ((kx == 0) .and. (kz == 0)) then ! special case for kx = kz = 0
            vecout = (0.0_dp, 0.0_dp)
        else
            ! ve_X_H( vec_X_A2( vec ))
            temp1 = matmul( vec, Hvv )

            ! Apply Boundary conditions
            ! Since this is left multiplication, the boundary conditions should be applied to the input vector
            temp1(    1) = (0.0_dp, 0.0_dp) !    v( 1 ) = 0
            temp1(myf  ) = (0.0_dp, 0.0_dp) !    v(myf) = 0
            temp1(    2) = (0.0_dp, 0.0_dp) ! dvdy( 1 ) = 0
            temp1(myf-1) = (0.0_dp, 0.0_dp) ! dvdy(myf) = 0

            ! left multiplication with the derivative matrix
            temp2 =  matmul( temp1, D1 ) ! note this is a left multiplication

            vecout(      1:  myf) = - (complex_i*kx)*temp2
            vecout(  myf+1:2*myf) = -            k2 *temp1
            vecout(2*myf+1:3*myf) = - (complex_i*kz)*temp2
        endif
    end function vec_X_A2HB

    function vec_X_A3HB(vec, kx, kz, Hvv, Hev, Hee) result(vecout)
        real(kind=dp), intent(in) :: kx, kz
        complex(kind=dp), intent(in), dimension(myf) :: vec
        complex(kind=dp), dimension(myf, myf), intent(in) :: Hvv, Hev, Hee
        complex(kind=dp), dimension(3*myf) :: vecout

        real(kind=dp) :: k2
        k2 = kx**2 + kz**2

        if ((kx == 0) .and. (kz == 0)) then ! special case for kx = kz = 0
            vecout(      1:  myf) = matmul( vec, Hev )
            vecout(  myf+1:2*myf) = (0.0_dp, 0.0_dp)
            vecout(2*myf+1:3*myf) = matmul( vec, Hee )
            ! Apply Boundary conditions
            ! Since this is left multiplication, the boundary conditions should be applied to the input vector
            ! but for this special case, B are made up of identity matrices, so it does not matter
            vecout(      1) = (0.0_dp, 0.0_dp) !   u( 1 ) = 0
            vecout(  myf  ) = (0.0_dp, 0.0_dp) !   u(myf) = 0
            vecout(2*myf+1) = (0.0_dp, 0.0_dp) !   w( 1 ) = 0
            vecout(3*myf  ) = (0.0_dp, 0.0_dp) !   w(myf) = 0
        else
            ! vec_X_A3
            temp1 = + (complex_i*kz/k2)*matmul( vec, D1 ) ! note this is left multiplication
            temp2 = + (complex_i*kx/k2)*vec
            ! vec_X_H
            temp3 = matmul( temp1, Hvv ) + matmul( temp2, Hev )
            temp4 = matmul( temp2, Hee )

            ! Apply Boundary conditions
            ! Since this is left multiplication, the boundary conditions should be applied to the input vector
            temp3(    1) = (0.0_dp, 0.0_dp) !    v( 1 ) = 0
            temp3(myf  ) = (0.0_dp, 0.0_dp) !    v(myf) = 0
            temp3(    2) = (0.0_dp, 0.0_dp) ! dvdy( 1 ) = 0
            temp3(myf-1) = (0.0_dp, 0.0_dp) ! dvdy(myf) = 0
            temp4(    1) = (0.0_dp, 0.0_dp) !   o2( 1 ) = 0
            temp4(myf  ) = (0.0_dp, 0.0_dp) !   o2(myf) = 0

            ! left multiplication with the derivative matrix
            temp1 =  matmul( temp3, D1 ) ! note this is a left multiplication

            vecout(      1:  myf) = - (complex_i*kx)*temp1 + (complex_i*kz)*temp4
            vecout(  myf+1:2*myf) =                        -            k2 *temp3
            vecout(2*myf+1:3*myf) = - (complex_i*kz)*temp1 - (complex_i*kx)*temp4
        endif
    end function vec_X_A3HB

    function AHB_X_vec(kx, kz, Hvv, Hev, Hee, vec) result(vecout)
        real(kind=dp), intent(in) :: kx, kz
        complex(kind=dp), intent(in), dimension(3*myf) :: vec
        complex(kind=dp), dimension(myf, myf), intent(in) :: Hvv, Hev, Hee
        complex(kind=dp), dimension(3*myf) :: vecout

        real(kind=dp) :: k2
        k2 = kx**2 + kz**2

        if ((kx == 0) .and. (kz == 0)) then ! special case for kx = kz = 0
            temp1 = vec(      1:  myf)
            temp2 = vec(2*myf+1:3*myf)
            ! Apply Boundary conditions
            temp1(    1) = (0.0_dp, 0.0_dp) !   u( 1 ) = 0
            temp1(myf  ) = (0.0_dp, 0.0_dp) !   u(myf) = 0
            temp2(    1) = (0.0_dp, 0.0_dp) !   w( 1 ) = 0
            temp2(myf  ) = (0.0_dp, 0.0_dp) !   w(myf) = 0

            vecout(      1:  myf) = matmul( Hvv, temp1 ) ! v
            vecout(  myf+1:2*myf) = (0.0_dp, 0.0_dp)
            vecout(2*myf+1:3*myf) = matmul( Hev, temp1 ) + matmul( Hee, temp2 ) ! o2
        else
            ! B_X_vec
            temp3 = - (complex_i*kx)*vec(      1:  myf) - (complex_i*kz)*vec(2*myf+1:3*myf)
            temp1 =   d_dy_complexvector(temp3)         -            k2 *vec(  myf+1:2*myf)
            temp2 = + (complex_i*kz)*vec(      1:  myf) - (complex_i*kx)*vec(2*myf+1:3*myf)
            ! Apply Boundary conditions
            temp1(    1) = (0.0_dp, 0.0_dp) !    v( 1 ) = 0
            temp1(myf  ) = (0.0_dp, 0.0_dp) !    v(myf) = 0
            temp1(    2) = (0.0_dp, 0.0_dp) ! dvdy( 1 ) = 0
            temp1(myf-1) = (0.0_dp, 0.0_dp) ! dvdy(myf) = 0
            temp2(    1) = (0.0_dp, 0.0_dp) !   o2( 1 ) = 0
            temp2(myf  ) = (0.0_dp, 0.0_dp) !   o2(myf) = 0

            ! H_X_vec
            temp3 = matmul( Hvv, temp1 ) ! v
            temp4 = matmul( Hev, temp1 ) + matmul( Hee, temp2 ) ! o2

            ! A_X_vec
            ! y derivative of vec
            temp1 = d_dy_complexvector( temp3 )
            vecout(      1:  myf) = + (complex_i*kx/k2)*temp1 - (complex_i*kz/k2)*temp4
            vecout(  myf+1:2*myf) = temp3
            vecout(2*myf+1:3*myf) = + (complex_i*kz/k2)*temp1 + (complex_i*kx/k2)*temp4
        endif
    end function AHB_X_vec


end module resolvent