module resolvent
    use types, only: dp
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

    public :: initialize_resolvent_module
    public :: resolvent_A, resolvent_B
    public :: resolvent_LAP, resolvent_LOS, resolvent_LSQ
    public :: resolvent_H, resolvent_H00
    public :: ve_X_H, H_X_fve

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
    !   A      [complex, 3Ny X 2Ny, output] the A matrix
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
    !   LAP    [real, Ny X Ny, output] the Laplacian matrix
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
    !   LOS [complex, Ny X Ny, output] the Orr-Sommerfeld matrix
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
    !   LSQ [complex, Ny X Ny, output] the Squire matrix
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
        complex(kind=dp), dimension(myf,myf) :: temp

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
    ! This function computes fve = ve X H (left multiplication of H by vector ve)
    ! Arguments
    !   ve            [complex, size 2*Ny, Input] velocity vorticity vector
    !   Hvv, Hev, Hee [complex, Ny X Ny, Input] Resolvent sub blocks
    ! Return
    !   fve           [complex, size 2*Ny, Output] fve = ve X H (left multiplication of H by vector ve)
    function ve_X_H( ve, Hvv, Hev, Hee ) result( fve )
        complex(kind=dp), dimension(2*myf), intent(in) :: ve
        complex(kind=dp), dimension(myf, myf), intent(in) :: Hvv, Hev, Hee
        complex(kind=dp), dimension(2*myf) :: fve

        ! Note: matmul does not need to transpose the vector
        fve(    1:  myf) = matmul( ve(    1:  myf), Hvv ) + matmul( ve(myf+1:2*myf), Hev )
        fve(myf+1:2*myf) = matmul( ve(myf+1:2*myf), Hee )
    end function


    ! This function computes ve = H X fve (right multiplication of H by vector fve)
    ! Arguments
    !   fve           [complex, size 2*Ny, Input] velocity vorticity forcing vector
    !   Hvv, Hev, Hee [complex, Ny X Ny, Input] Resolvent sub blocks
    ! Return
    !   ve            [complex, size 2*Ny, Output] ve = H X fve (right multiplication of H by vector fve)
    function H_X_fve( Hvv, Hev, Hee, fve ) result( ve )
        complex(kind=dp), dimension(2*myf), intent(in) :: fve
        complex(kind=dp), dimension(myf, myf), intent(in) :: Hvv, Hev, Hee
        complex(kind=dp), dimension(2*myf) :: ve

        ve(    1:  myf) = matmul( Hvv, fve(1:myf) ) ! v
        ve(myf+1:2*myf) = matmul( Hev, fve(1:myf) ) + matmul( Hee, fve(myf+1:2*myf) ) ! o2
    end function

end module resolvent