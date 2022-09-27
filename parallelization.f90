module parallelization
    use types, only: sp, dp
    use mpi
    implicit none
    private
#   include "parameters"

    public :: pointers, distribute_velocity, change_dim2slice_to_dim3slice


contains
    ! subroutine pointers( begin_vec, end_vec, totalgridpoints )
    ! This function divides a grid for all the processors
    ! If the number of grid points are not divisible by the number of processors,
    ! some processors will have one extra grid point
    ! Arguments:
    !   begin_vec      : [Integer, size numerop, Ouput] the start index for each processor
    !   end_vec        : [Integer, size numerop, Ouput] the end   index for each processor
    !   totalgridpoints: [Integer, Input] the total number of grid points for this coordinate
    subroutine pointers( begin_vec, end_vec, totalgridpoints )
        integer, intent(out), dimension(numerop) :: begin_vec, end_vec
        integer, intent(in) :: totalgridpoints

        integer n,n1,n2


        ! divide all the points
        n1 = totalgridpoints/numerop
        ! integer division, truncated towards 0, this is the minimum amount of grid points each processor will get
        n2 = totalgridpoints-numerop*n1 ! the first n2 processors will get one extra grid point

        ! first processor starts at 1
        begin_vec(1)=1

        ! first n2 processors have n1+1 pointss
        do n = 1, n2
            end_vec(n)     = begin_vec(n) + n1
            begin_vec(n+1) = end_vec(n) + 1
        enddo

        ! remaining processors get n1 points
        do n = n2+1, numerop-1
            end_vec(n)     = begin_vec(n) + n1 - 1
            begin_vec(n+1) = end_vec(n) + 1
        enddo

        ! last processor end point
        end_vec(numerop) = begin_vec(numerop) + n1 - 1

    end subroutine pointers


    ! subroutine distribute_velocity( myid, sourceid, vel_full, vel_slice )
    ! This function distribute a velocity field from sourceid to each processor
    ! Each processor will get full xy planes (each get a few kz wavenumbers)
    ! Arguments:
    !   myid:      [integer, Input]
    !              processor ID
    !   sourceid:  [integer, Input]
    !              the processor ID where the full data field is saved
    !   vel_full:  [double/single complex, size (mxf,  mzf,myf), Input]
    !              the full velocity field, read from h5, interpolated and normalized
    !   vel_slice: [double/single complex, size (mxf,kb:ke,myf), Output]
    !              velocity data for each processor, full xy planes with a few kz
    subroutine distribute_velocity( myid, sourceid, vel_full, vel_slice )
        ! Global Pointers
        integer :: jbeg,jend,kbeg,kend,jb,je,kb,ke
        common /point/ jbeg(0:numerop-1),jend(0:numerop-1), &
                       kbeg(0:numerop-1),kend(0:numerop-1), &
                       jb,je,kb,ke
        save /point/
        ! Input/Output
        integer, intent(in) :: myid, sourceid
        complex(kind=cp), intent( in), dimension(mxf,  mzf,myf) :: vel_full
        complex(kind=cp), intent(out), dimension(mxf,kb:ke,myf) :: vel_slice

        ! send/receive data size
        integer :: nsend, nreceive
        ! send/receive buffers
        real(kind=cp), dimension(:,:,:,:), allocatable :: send_buffer, receive_buffer
        ! kb, ke for the target processor
        integer :: kb_target, ke_target
        ! Loop index
        integer :: iproc

        ! MPI
        integer :: istat(MPI_STATUS_SIZE), ierr


        if ( myid .eq. sourceid ) then
            ! if I am the source, then copy the data I need
            vel_slice = vel_full(:,kb:ke,:)

            ! sourceid send data to each processor except for itself
            ! each processor gets (mxf,kb_target:ke_target,myf) (full xy planes)
            do iproc = numerop-1,0, -1
                if ( iproc .ne. sourceid ) then ! no need to send data to itself
                    ! kb, ke for the target processor
                    kb_target = kbeg(iproc)
                    ke_target = kend(iproc)
                    ! number of data sent (*2 for real and imaginary parts)
                    nsend = mxf*(ke_target-kb_target+1)*myf*2

                    ! allocate a buffer to contain both the real and imaginary parts
                    ALLOCATE(send_buffer(mxf,kb_target:ke_target,myf,2))
                    send_buffer(:,:,:,1) =  real( vel_full(:,kb_target:ke_target,:), cp)
                    send_buffer(:,:,:,2) = aimag( vel_full(:,kb_target:ke_target,:) )
                    ! sourceid send data to iproc
                    if ( cp .eq. dp ) then
                        call MPI_SEND(send_buffer, nsend, MPI_DOUBLE_PRECISION, iproc, 0, MPI_COMM_WORLD, ierr )
                    else if ( cp .eq. sp ) then
                        call MPI_SEND(send_buffer, nsend, MPI_REAL            , iproc, 0, MPI_COMM_WORLD, ierr )
                    endif
                    ! deallocate buffer
                    DEALLOCATE(send_buffer)
                endif
            enddo
        else
            ! All other processors receive data from sourceid
            ! number of data received (*2 for real and imaginary parts)
            nreceive = mxf*(ke-kb+1)*myf*2
            ! allocate a buffer to contain both the real and imaginary parts
            ALLOCATE(receive_buffer(mxf,kb:ke,myf,2))
            ! receive data from sourceid
            if ( cp .eq. dp ) then
                call MPI_RECV(receive_buffer, nreceive, MPI_DOUBLE_PRECISION, sourceid, 0, MPI_COMM_WORLD, istat, ierr )
            else if ( cp .eq. sp ) then
                call MPI_RECV(receive_buffer, nreceive, MPI_REAL            , sourceid, 0, MPI_COMM_WORLD, istat, ierr )
            endif
            ! build output from both real and imaginary parts
            vel_slice = CMPLX( receive_buffer(:,:,:,1), receive_buffer(:,:,:,2), cp)
            DEALLOCATE(receive_buffer)
        endif
    end subroutine distribute_velocity


    ! subroutine change_dim2slice_to_dim3slice( dim2slice, myid, dim3slice )
    ! This function exchange the data between all processors
    ! each processor originally has all the data for dimension 1 and 3 ( but only a few dim 2 grid points )
    ! after the exchange, each processor has all the data for dimension 1 and 2 ( but only a few dim 3 grid points )
    ! Arguments:
    !   dim2slice: [double/single, size (size_dim1, size_dim2slice, size_dim3full ), Input]
    !              all the data for dimension 1 and 3 ( but only a few dim 2 grid points )
    !   myid:      [integer, Input]
    !              processor ID
    !   dim3slice: [double/single, size (size_dim1, size_dim2full , size_dim3slice), Output]
    !              all the data for dimension 1 and 2 ( but only a few dim 3 grid points )
    subroutine change_dim2slice_to_dim3slice( dim2slice, myid, dim3slice )
        ! Input/outpus
        complex(kind=cp), intent( in), dimension(:,:,:) :: dim2slice ! size (size_dim1, size_dim2slice, size_dim3full )
        complex(kind=cp), intent(out), dimension(:,:,:) :: dim3slice ! size (size_dim1, size_dim2full , size_dim3slice)
        integer, intent(in) :: myid

        ! Grid sizes
        integer, dimension(3) :: size_matrix
        integer :: size_dim1
        integer :: size_dim2slice, size_dim2full
        integer :: size_dim3slice, size_dim3full

        ! Pointers
        integer :: b2v(0:numerop-1), e2v(0:numerop-1)
        integer :: b3v(0:numerop-1), e3v(0:numerop-1)
        integer :: b2, e2, b3, e3

        ! temp b2 e2 b3 e3 pointers
        integer :: b2_target, e2_target
        integer :: b3_target, e3_target

        ! buffer to package real with imaginary
        real(kind=cp), dimension(:,:,:,:), allocatable :: send_buffer
        real(kind=cp), dimension(:,:,:,:), allocatable :: recv_buffer

        ! send and recieve data size
        integer :: nsend, nrecv

        ! MPI variables
        integer istat(MPI_STATUS_SIZE), ierr

        ! Loop variable
        integer iproc


        ! --------------------- Get Full Matrix Dimensions ---------------------
        size_matrix    = shape(dim2slice) ! Input matrix size
        size_dim1      = size_matrix(1)   ! size of dimension 1
        size_dim2slice = size_matrix(2)   ! size of dimension 2 slice
        size_dim3full  = size_matrix(3)   ! full size of dimension 3
        ! compute the full dimension 2 size by adding all the dimension 2 slice sizes
        call MPI_ALLREDUCE(size_dim2slice, size_dim2full, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)


        ! ------------------ Compute Dimension 2 and 3 Slices ------------------
        ! compute the vector pointers for dim 2 and 3
        ! these vectors contain the begin and end indices for all processors
        call pointers( b2v, e2v, size_dim2full ) ! dimension 2 begin and end indices
        call pointers( b3v, e3v, size_dim3full ) ! dimension 3 begin and end indices
        ! begin and end indices for the current processor
        b2 = b2v(myid)
        e2 = e2v(myid)
        b3 = b3v(myid)
        e3 = e3v(myid)


        ! ----------------- Allocate send and receive buffers -----------------
        ALLOCATE( send_buffer( size_dim1, b2:e2, size_dim3full, 2 ) )
        ALLOCATE( recv_buffer( size_dim1, size_dim2full, b3:e3, 2 ) )


        ! ------------- Package Real and Imaginary Part of My Data -------------
        send_buffer(:,b2:e2,:,1) =  real( dim2slice(:,:,:), cp)
        send_buffer(:,b2:e2,:,2) = aimag( dim2slice(:,:,:) )


        ! ----------------- Get Part of Data from Myself First -----------------
        recv_buffer(:, b2:e2, b3:e3, :) = send_buffer(:, b2:e2, b3:e3, :)


        ! ------------- Loop Over All Processors Except For Myself -------------
        ! each processor will exchange data with me
        do iproc = 0, numerop - 1
            if ( iproc .ne. myid ) then
                ! I send data to iproc
                ! and recieve data from iproc

                ! b2, e2, b3, e3 for iproc
                b2_target = b2v(iproc)
                e2_target = e2v(iproc)
                b3_target = b3v(iproc)
                e3_target = e3v(iproc)

                ! data sent from myid to iproc is size (size_dim1, b2:e2, b3_target:e3_target)
                ! I only have data for b2:e2,
                ! Although I have data for all dim 3, iproc only need data for b3_target:e3_target

                ! data recieved from iproc is size (size_dim1, b2_target:e2_target, b3:e3)
                ! I only need data for b3:e3
                ! and iproc only has data for b2_target:e2_target

                ! send and recieve data size, *2 for real and imaginary
                nsend = size_dim1 * (e2-b2+1) * (e3_target-b3_target+1) * 2
                nrecv = size_dim1 * (e2_target-b2_target+1) * (e3-b3+1) * 2

                ! exchange data
                if ( cp .eq. dp ) then
                    call MPI_SENDRECV( &
                        send_buffer(:, b2:e2, b3_target:e3_target,:) ,nsend,MPI_DOUBLE_PRECISION,iproc,0, &
                        recv_buffer(:, b2_target:e2_target, b3:e3,:) ,nrecv,MPI_DOUBLE_PRECISION,iproc,0, &
                        MPI_COMM_WORLD,istat,ierr)
                else if ( cp .eq. sp ) then
                    call MPI_SENDRECV( &
                        send_buffer(:, b2:e2, b3_target:e3_target,:) ,nsend,MPI_REAL            ,iproc,0, &
                        recv_buffer(:, b2_target:e2_target, b3:e3,:) ,nrecv,MPI_REAL            ,iproc,0, &
                        MPI_COMM_WORLD,istat,ierr)
                endif
            endif
        enddo


        ! ------------------------ Build Complex Output ------------------------
        dim3slice = CMPLX( recv_buffer(:,:,b3:e3,1), recv_buffer(:,:,b3:e3,2), cp)


        ! -------------------------- Deallocate buffer --------------------------
        DEALLOCATE( send_buffer )
        DEALLOCATE( recv_buffer )

    end subroutine


end module parallelization