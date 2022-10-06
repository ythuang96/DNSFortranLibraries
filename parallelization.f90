module parallelization
    use types, only: sp, dp
    use mpi
    implicit none
    private
#   include "parameters"

    public :: pointers, distribute_velocity_slicedim2, distribute_velocity_slicedim3
    public :: change_dim2slice_to_dim3slice, change_dim3slice_to_dim2slice
    public :: changeR_dim3slice_to_dim2slice


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


    ! subroutine distribute_velocity_slicedim2( myid, sourceid, vel_full, vel_slice )
    ! This function distribute a velocity field from sourceid to each processor
    ! Each processor will get full data for dimension 1 and 3,
    ! but only get a few grid points for dimension 2
    ! Arguments:
    !   myid:      [integer, Input]
    !              processor ID
    !   sourceid:  [integer, Input]
    !              the processor ID where the full data field is saved
    !   vel_full:  [double/single complex, size (size_dim1,size_dim2,size_dim3), Input]
    !              the full velocity field, read from h5, interpolated and normalized
    !   vel_slice: [double/single complex, size (size_dim1,b2:e2,size_dim3), Output]
    !              velocity data for each processor, full data for dimension 1 and 3,
    !              but only a few grid points for dimension 2
    subroutine distribute_velocity_slicedim2( myid, sourceid, vel_full, vel_slice )
        ! Input/Output
        integer, intent(in) :: myid, sourceid
        complex(kind=cp), intent( in), dimension(:,:,:) :: vel_full
        complex(kind=cp), intent(out), dimension(:,:,:) :: vel_slice

        ! Matrix sizes
        integer, dimension(3) :: size_matrix
        integer :: size_dim1, size_dim2, size_dim3

        ! Pointers for dimension 2 distribution
        integer :: b2v(0:numerop-1), e2v(0:numerop-1)
        integer :: b2, e2

        ! temp b2 e2 pointers for the target processor
        integer :: b2_target, e2_target

        ! send/receive data size
        integer :: nsend, nreceive
        ! send/receive buffers
        real(kind=cp), dimension(:,:,:,:), allocatable :: send_buffer, receive_buffer
        ! Loop index
        integer :: iproc

        ! MPI
        integer :: istat(MPI_STATUS_SIZE), ierr


        ! --------------------- Get Full Matrix Dimensions ---------------------
        size_matrix = shape(vel_full) ! Input matrix size
        size_dim1   = size_matrix(1)  ! size of dimension 1
        size_dim2   = size_matrix(2)  ! size of dimension 2
        size_dim3   = size_matrix(3)  ! size of dimension 3


        ! --------------------- Compute Dimension 2 Slices ---------------------
        ! compute the vector pointers for dim 2
        ! these vectors contain the begin and end indices for all processors
        call pointers( b2v, e2v, size_dim2 ) ! dimension 2 begin and end indices
        ! begin and end indices for the current processor
        b2 = b2v(myid)
        e2 = e2v(myid)


        ! ----------------------- Send and Receive Data -----------------------
        if ( myid .eq. sourceid ) then
            ! if I am the source, then copy the data I need
            vel_slice = vel_full(:,b2:e2,:)

            ! sourceid send data to each processor except for itself
            ! each processor gets (size_dim1,b2_target:e2_target,size_dim3)
            ! (full dim 1 and 3 data, but only at a few dim 2 grid points)
            do iproc = numerop-1, 0, -1
                if ( iproc .ne. sourceid ) then ! no need to send data to itself
                    ! b2, e2 for the target processor
                    b2_target = b2v(iproc)
                    e2_target = e2v(iproc)
                    ! number of data sent (*2 for real and imaginary parts)
                    nsend = size_dim1 * (e2_target-b2_target+1) * size_dim3 * 2

                    ! allocate a buffer to contain both the real and imaginary parts
                    ALLOCATE(send_buffer(size_dim1, b2_target:e2_target, size_dim3, 2))
                    send_buffer(:,:,:,1) =  real( vel_full(:,b2_target:e2_target,:), cp)
                    send_buffer(:,:,:,2) = aimag( vel_full(:,b2_target:e2_target,:) )
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
            nreceive = size_dim1 * (e2-b2+1) * size_dim3 * 2
            ! allocate a buffer to contain both the real and imaginary parts
            ALLOCATE(receive_buffer(size_dim1, b2:e2, size_dim3, 2))
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
    end subroutine distribute_velocity_slicedim2


    ! subroutine distribute_velocity_slicedim3( myid, sourceid, vel_full, vel_slice )
    ! This function distribute a velocity field from sourceid to each processor
    ! Each processor will get full data for dimension 1 and 2,
    ! but only get a few grid points for dimension 3
    ! Arguments:
    !   myid:      [integer, Input]
    !              processor ID
    !   sourceid:  [integer, Input]
    !              the processor ID where the full data field is saved
    !   vel_full:  [double/single complex, size (size_dim1,size_dim2,size_dim3), Input]
    !              the full velocity field, read from h5, interpolated and normalized
    !   vel_slice: [double/single complex, size (size_dim1,size_dim2,b3:e3), Output]
    !              velocity data for each processor, full data for dimension 1 and 2,
    !              but only a few grid points for dimension 3
    subroutine distribute_velocity_slicedim3( myid, sourceid, vel_full, vel_slice )
        ! Input/Output
        integer, intent(in) :: myid, sourceid
        complex(kind=cp), intent( in), dimension(:,:,:) :: vel_full
        complex(kind=cp), intent(out), dimension(:,:,:) :: vel_slice

        ! Matrix sizes
        integer, dimension(3) :: size_matrix
        integer :: size_dim1, size_dim2, size_dim3

        ! Pointers for dimension 2 distribution
        integer :: b3v(0:numerop-1), e3v(0:numerop-1)
        integer :: b3, e3

        ! temp b2 e2 pointers for the target processor
        integer :: b3_target, e3_target

        ! send/receive data size
        integer :: nsend, nreceive
        ! send/receive buffers
        real(kind=cp), dimension(:,:,:,:), allocatable :: send_buffer, receive_buffer
        ! Loop index
        integer :: iproc

        ! MPI
        integer :: istat(MPI_STATUS_SIZE), ierr


        ! --------------------- Get Full Matrix Dimensions ---------------------
        size_matrix = shape(vel_full) ! Input matrix size
        size_dim1   = size_matrix(1)  ! size of dimension 1
        size_dim2   = size_matrix(2)  ! size of dimension 2
        size_dim3   = size_matrix(3)  ! size of dimension 3


        ! --------------------- Compute Dimension 3 Slices ---------------------
        ! compute the vector pointers for dim 3
        ! these vectors contain the begin and end indices for all processors
        call pointers( b3v, e3v, size_dim3 ) ! dimension 3 begin and end indices
        ! begin and end indices for the current processor
        b3 = b3v(myid)
        e3 = e3v(myid)


        ! ----------------------- Send and Receive Data -----------------------
        if ( myid .eq. sourceid ) then
            ! if I am the source, then copy the data I need
            vel_slice = vel_full(:,:,b3:e3)

            ! sourceid send data to each processor except for itself
            ! each processor gets (size_dim1,size_dim2,b3_target:e3_target)
            ! (full dim 1 and 2 data, but only at a few dim 3 grid points)
            do iproc = numerop-1, 0, -1
                if ( iproc .ne. sourceid ) then ! no need to send data to itself
                    ! b3, e3 for the target processor
                    b3_target = b3v(iproc)
                    e3_target = e3v(iproc)
                    ! number of data sent (*2 for real and imaginary parts)
                    nsend = size_dim1 * size_dim2 * (e3_target-b3_target+1) * 2

                    ! allocate a buffer to contain both the real and imaginary parts
                    ALLOCATE(send_buffer(size_dim1, size_dim2, b3_target:e3_target, 2))
                    send_buffer(:,:,:,1) =  real( vel_full(:,:,b3_target:e3_target), cp)
                    send_buffer(:,:,:,2) = aimag( vel_full(:,:,b3_target:e3_target) )
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
            nreceive = size_dim1 * size_dim2 * (e3-b3+1) * 2
            ! allocate a buffer to contain both the real and imaginary parts
            ALLOCATE(receive_buffer(size_dim1, size_dim2, b3:e3, 2))
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
    end subroutine distribute_velocity_slicedim3


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


    ! subroutine change_dim3slice_to_dim2slice( dim3slice, myid, dim2slice )
    ! This function exchange the data between all processors
    ! each processor originally has all the data for dimension 1 and 2 ( but only a few dim 3 grid points )
    ! after the exchange, each processor has all the data for dimension 1 and 3 ( but only a few dim 2 grid points )
    ! Arguments:
    !   dim3slice: [double/single, size (size_dim1, size_dim2full , size_dim3slice), Input]
    !              all the data for dimension 1 and 2 ( but only a few dim 3 grid points )
    !   myid:      [integer, Input]
    !              processor ID
    !   dim2slice: [double/single, size (size_dim1, size_dim2slice, size_dim3full ), Output]
    !              all the data for dimension 1 and 3 ( but only a few dim 2 grid points )
    subroutine change_dim3slice_to_dim2slice( dim3slice, myid, dim2slice )
        ! Input/outpus
        complex(kind=cp), intent( in), dimension(:,:,:) :: dim3slice ! size (size_dim1, size_dim2full , size_dim3slice)
        complex(kind=cp), intent(out), dimension(:,:,:) :: dim2slice ! size (size_dim1, size_dim2slice, size_dim3full )
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
        size_matrix    = shape(dim3slice) ! Input matrix size
        size_dim1      = size_matrix(1)   ! size of dimension 1
        size_dim2full  = size_matrix(2)   ! full size of dimension 2 slice
        size_dim3slice = size_matrix(3)   ! size of dimension 3 slice
        ! compute the full dimension 3 size by adding all the dimension 3 slice sizes
        call MPI_ALLREDUCE(size_dim3slice, size_dim3full, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)


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
        ALLOCATE( send_buffer( size_dim1, size_dim2full, b3:e3, 2 ) )
        ALLOCATE( recv_buffer( size_dim1, b2:e2, size_dim3full, 2 ) )


        ! ------------- Package Real and Imaginary Part of My Data -------------
        send_buffer(:,:,b3:e3,1) =  real( dim3slice(:,:,:), cp)
        send_buffer(:,:,b3:e3,2) = aimag( dim3slice(:,:,:) )


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

                ! data sent from myid to iproc is size (size_dim1, b2_target:e2_targer, b3:e3)
                ! I only have data for b3:e3,
                ! Although I have data for all dim 2, iproc only need data for b2_target:e2_target

                ! data recieved from iproc is size (size_dim1, b2:e2, b3_target:e3_target)
                ! I only need data for b2:e2
                ! and iproc only has data for b3_target:e3_target

                ! send and recieve data size, *2 for real and imaginary
                nsend = size_dim1 * (e2_target-b2_target+1) * (e3-b3+1) * 2
                nrecv = size_dim1 * (e2-b2+1) * (e3_target-b3_target+1) * 2

                ! exchange data
                if ( cp .eq. dp ) then
                    call MPI_SENDRECV( &
                        send_buffer(:, b2_target:e2_target, b3:e3,:) ,nsend,MPI_DOUBLE_PRECISION,iproc,0, &
                        recv_buffer(:, b2:e2, b3_target:e3_target,:) ,nrecv,MPI_DOUBLE_PRECISION,iproc,0, &
                        MPI_COMM_WORLD,istat,ierr)
                else if ( cp .eq. sp ) then
                    call MPI_SENDRECV( &
                        send_buffer(:, b2_target:e2_target, b3:e3,:) ,nsend,MPI_REAL            ,iproc,0, &
                        recv_buffer(:, b2:e2, b3_target:e3_target,:) ,nrecv,MPI_REAL            ,iproc,0, &
                        MPI_COMM_WORLD,istat,ierr)
                endif
            endif
        enddo


        ! ------------------------ Build Complex Output ------------------------
        dim2slice = CMPLX( recv_buffer(:,b2:e2,:,1), recv_buffer(:,b2:e2,:,2), cp)


        ! -------------------------- Deallocate buffer --------------------------
        DEALLOCATE( send_buffer )
        DEALLOCATE( recv_buffer )

    end subroutine


    ! subroutine changeR_dim3slice_to_dim2slice( dim3slice, myid, dim2slice )
    ! This function exchange the data between all processors
    ! each processor originally has all the data for dimension 1 and 2 ( but only a few dim 3 grid points )
    ! after the exchange, each processor has all the data for dimension 1 and 3 ( but only a few dim 2 grid points )
    !
    ! This funciton is for real data only, complex data should use change_dim3slice_to_dim2slice
    !
    ! Arguments:
    !   dim3slice: [double/single, size (size_dim1, size_dim2full, size_dim3slice), Input]
    !              all the data for dimension 1 and 2 ( but only a few dim 3 grid points )
    !   myid:      [integer, Input]
    !              processor ID
    !   dim2slice: [double/single, size (size_dim1, size_dim2slice, size_dim3full), Output]
    !              all the data for dimension 1 and 3 ( but only a few dim 2 grid points )
    subroutine changeR_dim3slice_to_dim2slice( dim3slice, myid, dim2slice )
        ! Input/outpus
        real(kind=cp), intent( in), dimension(:,:,:) :: dim3slice ! size (size_dim1, size_dim2full , size_dim3slice)
        real(kind=cp), intent(out), dimension(:,:,:) :: dim2slice ! size (size_dim1, size_dim2slice, size_dim3full )
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
        real(kind=cp), dimension(:,:,:), allocatable :: send_buffer
        real(kind=cp), dimension(:,:,:), allocatable :: recv_buffer

        ! send and recieve data size
        integer :: nsend, nrecv

        ! MPI variables
        integer istat(MPI_STATUS_SIZE), ierr

        ! Loop variable
        integer iproc


        ! --------------------- Get Full Matrix Dimensions ---------------------
        size_matrix    = shape(dim3slice) ! Input matrix size
        size_dim1      = size_matrix(1)   ! size of dimension 1
        size_dim2full  = size_matrix(2)   ! full size of dimension 2 slice
        size_dim3slice = size_matrix(3)   ! size of dimension 3 slice
        ! compute the full dimension 3 size by adding all the dimension 3 slice sizes
        call MPI_ALLREDUCE(size_dim3slice, size_dim3full, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)


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
        ALLOCATE( send_buffer( size_dim1, size_dim2full, b3:e3) )
        ALLOCATE( recv_buffer( size_dim1, b2:e2, size_dim3full) )


        ! -------------------- Package Real Part of My Data --------------------
        send_buffer(:,:,b3:e3) =  dim3slice(:,:,:)


        ! ----------------- Get Part of Data from Myself First -----------------
        recv_buffer(:, b2:e2, b3:e3) = send_buffer(:, b2:e2, b3:e3)


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

                ! data sent from myid to iproc is size (size_dim1, b2_target:e2_targer, b3:e3)
                ! I only have data for b3:e3,
                ! Although I have data for all dim 2, iproc only need data for b2_target:e2_target

                ! data recieved from iproc is size (size_dim1, b2:e2, b3_target:e3_target)
                ! I only need data for b2:e2
                ! and iproc only has data for b3_target:e3_target

                ! send and recieve data size
                nsend = size_dim1 * (e2_target-b2_target+1) * (e3-b3+1)
                nrecv = size_dim1 * (e2-b2+1) * (e3_target-b3_target+1)

                ! exchange data
                if ( cp .eq. dp ) then
                    call MPI_SENDRECV( &
                        send_buffer(:, b2_target:e2_target, b3:e3) ,nsend,MPI_DOUBLE_PRECISION,iproc,0, &
                        recv_buffer(:, b2:e2, b3_target:e3_target) ,nrecv,MPI_DOUBLE_PRECISION,iproc,0, &
                        MPI_COMM_WORLD,istat,ierr)
                else if ( cp .eq. sp ) then
                    call MPI_SENDRECV( &
                        send_buffer(:, b2_target:e2_target, b3:e3) ,nsend,MPI_REAL            ,iproc,0, &
                        recv_buffer(:, b2:e2, b3_target:e3_target) ,nrecv,MPI_REAL            ,iproc,0, &
                        MPI_COMM_WORLD,istat,ierr)
                endif
            endif
        enddo


        ! ------------------------ Build Complex Output ------------------------
        dim2slice = recv_buffer(:,b2:e2,:)


        ! -------------------------- Deallocate buffer --------------------------
        DEALLOCATE( send_buffer )
        DEALLOCATE( recv_buffer )

    end subroutine

end module parallelization