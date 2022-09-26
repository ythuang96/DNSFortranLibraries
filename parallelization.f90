module parallelization
    use types, only: sp, dp
    use mpi
    implicit none
    private
#   include "parameters"

    public :: pointers, distribute_velocity, change_xyplane2xzplane


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


    ! subroutine change_xyplane2xzplane( xy, myid, xz )
    ! This function exchange the data between all processors
    ! each processor originally has all the data in a few xy planes (a few kz wavenumbers)
    ! after the exchange, each processor has all the data in a few xz planes (a few y grid points)
    ! Arguments:
    !   xy:   [double/single, size (mxf,kb:ke,myf), Input]
    !         all the data in a few xy planes for this processor (a few kz wavenumbers)
    !   myid: [integer, Input]
    !         processor ID
    !   xz:   [double/single, size (mxf,mzf,jb:je), Output]
    !         all the data in a few xz planes for this processor (a few y grid points)
    subroutine change_xyplane2xzplane( xy, myid, xz )
        ! global pointer variables
        integer :: jbeg,jend,kbeg,kend,jb,je,kb,ke
        common /point/ jbeg(0:numerop-1),jend(0:numerop-1), &
                       kbeg(0:numerop-1),kend(0:numerop-1), &
                       jb,je,kb,ke
        save /point/
        ! Input/outpus
        complex(kind=cp), intent( in), dimension(mxf,kb:ke,myf) :: xy
        complex(kind=cp), intent(out), dimension(mxf,mzf,jb:je) :: xz
        integer, intent(in) :: myid

        ! buffer to package real with imaginary
        real(kind=cp), dimension(mxf,kb:ke,myf,2) :: send_buffer
        real(kind=cp), dimension(mxf,mzf,jb:je,2) :: recv_buffer

        ! temp jb je kb ke pointers
        integer :: jb_target, je_target, kb_target, ke_target

        ! send and recieve data size
        integer :: nsend, nrecv

        ! MPI variables
        integer istat(MPI_STATUS_SIZE), ierr

        ! Loop variable
        integer iproc


        ! package real and imaginary part of my data
        send_buffer(:,kb:ke,:,1) =  real( xy(:,kb:ke,:), cp)
        send_buffer(:,kb:ke,:,2) = aimag( xy(:,kb:ke,:) )

        ! get part of data from myself first
        recv_buffer(:, kb:ke, jb:je, :) = send_buffer(:, kb:ke, jb:je, :)

        ! Loop over all processors except for myself, each will exchange data with me
        do iproc=0,numerop-1
            if ( iproc .ne. myid ) then
                ! I send data to iproc
                ! and recieve data from iproc

                ! jb, je, kb, ke for iproc
                jb_target = jbeg(iproc)
                je_target = jend(iproc)
                kb_target = kbeg(iproc)
                ke_target = kend(iproc)

                ! data sent from myid to iproc is (mxf, kb:ke, jb_target:je_target)
                ! I only have data for kb:ke,
                ! Although I have data for all y, iproc only need data for jb_target:je_target

                ! data recieved from iproc is (mxf, kb_target:ke_target, jb:je)
                ! I only need data for jb:je
                ! and iproc only has data for kb_target:ke_target

                ! send and recieve data size, *2 for real and imaginary
                nsend = mxf*(ke-kb+1)*(je_target-jb_target+1)*2
                nrecv = mxf*(ke_target-kb_target+1)*(je-jb+1)*2

                ! exchange data
                if ( cp .eq. dp ) then
                    call MPI_SENDRECV( &
                        send_buffer(:, kb:ke, jb_target:je_target,:) ,nsend,MPI_DOUBLE_PRECISION,iproc,0, &
                        recv_buffer(:, kb_target:ke_target, jb:je,:) ,nrecv,MPI_DOUBLE_PRECISION,iproc,0, &
                        MPI_COMM_WORLD,istat,ierr)
                else if ( cp .eq. sp ) then
                    call MPI_SENDRECV( &
                        send_buffer(:, kb:ke, jb_target:je_target,:) ,nsend,MPI_REAL            ,iproc,0, &
                        recv_buffer(:, kb_target:ke_target, jb:je,:) ,nrecv,MPI_REAL            ,iproc,0, &
                        MPI_COMM_WORLD,istat,ierr)
                endif
            endif
        enddo

        ! Build complex output
        xz = CMPLX( recv_buffer(:,:,jb:je,1), recv_buffer(:,:,jb:je,2), cp)
    end subroutine


end module parallelization