    ! subroutine h5save_C3Serial_dim3( filename, varname, myid, matrix )
    ! save a complex rank 3 matrix to h5 file in serial.
    ! with parallization in dimension 3 (Each processor has data for a few
    ! dimension 3 grid points, while having all data for dimension 1 and 2.)
    ! Data will be saved on dim 3 grid point at a time.
    !
    ! Funciton wise is the same as h5save_C3Parallel_dim3
    !
    ! All processors should call this function.
    !
    ! Arguments:
    !   filename: [string, Input] h5 filename with path
    !   varname : [string, Input] variable name
    !   myid    : [integer, Input] processor ID
    !   matrix  : [double/single complex 3d matrix, Input] data to be saved
    subroutine h5save_C3Serial_dim3( filename, varname, myid, matrix)
        use parallelization, only: pointers
        ! Inputs
        character(len=*), intent(in) :: filename, varname
        integer, intent(in) :: myid
        complex(kind=cp), intent(in), dimension(:,:,:) :: matrix

        ! data dimensions for full data and slice data
        integer, dimension(3) :: full_data_dim, slice_data_dim
        integer :: dim3slice, dim3full

        ! work assignment pointers for dimension 3
        integer :: b3v(0:numerop-1), e3v(0:numerop-1)
        integer :: b3, e3

        ! MPI variables
        integer :: ierr
        ! Loop index
        integer :: ii
        ! Data consistency check
        logical :: dataagree


        ! --------------------- Get Full Matrix Dimensions ---------------------
        ! get matrix dimension for the slices that this current processor has
        slice_data_dim = shape(matrix)
        dim3slice = slice_data_dim(3)
        ! compute the full dimension 3 size by adding all the dimension 3 slice sizes
        call MPI_ALLREDUCE(dim3slice, dim3full, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
        full_data_dim(1) = slice_data_dim(1)
        full_data_dim(2) = slice_data_dim(2)
        full_data_dim(3) = dim3full

        ! -------------- Compute the distribution in dimension 3 --------------
        call pointers( b3v, e3v, dim3full)
        b3 = b3v(myid)
        e3 = e3v(myid)

        ! --------------- Initialize variable for partial saving ---------------
        if (myid .eq. 0) then
            call h5save_CPartial_Init( filename, varname, full_data_dim)
        endif

        ! ------------- Save data to h5 and check data consistency -------------
        ! This ensures the data saving is executed at least once
        dataagree = .false.

        ! Loop until h5 file passes data consistency check
        DO WHILE ( dataagree .eqv. .false. )

            ! Loop over dim3 and save each dim 3
            DO ii = 1, dim3full
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                if ( (ii .ge. b3) .and. (ii .le. e3) ) then
                    ! if this plane is part of the data
                    call h5save_C3Partial_SingleDim3( filename, varname, matrix(:,:,ii-b3+1), ii )
                endif
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            ENDDO

            ! verify data consistency
            dataagree = verifysave_C3_dim3(filename, varname, myid, matrix, full_data_dim)
            if ( myid .eq. 0 ) then
                if ( dataagree ) then
                    write(*,*) "   Data consistency check passed."
                else
                    write(*,*) "   Data consistency check failed, re-saving data."
                endif
            endif
        ENDDO
    end subroutine h5save_C3Serial_dim3


    ! function verifysave_C3_dim3(filename, varname, myid, matrix, full_data_dim) result( dataagree_all )
    ! Read data from the h5 file and compare to the data in memory to ensure
    ! no data corruption in the saving process
    !
    ! Arguments:
    !   filename     : [string, Input] h5 filename with path
    !   varname      : [string, Input] variable name
    !   myid         : [integer, Input] processor ID
    !   matrix       : [double/single complex 3d matrix, Input] data saved
    !   full_data_dim: [integer, size (3)] size of the full rank 3 matrix
    ! Return:
    !   dataagree_all: [logical] true if the data in the h5 file matches with what was saved
    !                            false if data corruption occured
    function verifysave_C3_dim3(filename, varname, myid, matrix, full_data_dim) result( dataagree_all )
        use parallelization, only: pointers, distribute_C3_slicedim3
        use h5load, only: h5load_C3
        ! Inputs
        character(len=*), intent(in) :: filename, varname
        integer, intent(in) :: myid
        complex(kind=cp), intent(in), dimension(:,:,:) :: matrix
        integer, intent(in), dimension(3) :: full_data_dim
        ! Output
        logical :: dataagree_local, dataagree_all

        ! work assignment pointers for dimension 3
        integer :: b3v(0:numerop-1), e3v(0:numerop-1)
        integer :: b3, e3
        ! full data read from the file
        complex(kind=cp), dimension(:,:,:), allocatable :: full_data, slice_data
        ! MPI variables
        integer :: ierr


        ! -------------- Compute the distribution in dimension 3 --------------
        call pointers( b3v, e3v, full_data_dim(3))
        b3 = b3v(myid)
        e3 = e3v(myid)

        ! ---------------- Read from h5 and distribute in dim 3 ----------------
        allocate(  full_data( full_data_dim(1),full_data_dim(2),full_data_dim(3) ) )
        allocate( slice_data( full_data_dim(1),full_data_dim(2),b3:e3            ) )
        ! Master reads data from the h5 file
        if (myid .eq. 0) then
            full_data = h5load_C3(filename, varname )
        endif
        ! distribute data in dimension 3
        call distribute_C3_slicedim3( myid, 0, full_data, slice_data )
        deallocate( full_data )

        ! ---------------------------- Compare data ----------------------------
        if ( all(slice_data .eq. matrix) ) then
            dataagree_local = .true.
        else
            write(*,"(A,I2.2)") "   Data corruption at processor: ", myid
            dataagree_local = .false.
        endif
        deallocate( slice_data )

        ! ------------------ Collect data from all processors ------------------
        call MPI_ALLREDUCE(dataagree_local, dataagree_all, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    end function verifysave_C3_dim3


    ! subroutine h5save_C3Serial_dim2( filename, varname, myid, matrix )
    ! save a complex rank 3 matrix to h5 file in serial.
    ! with parallization in dimension 2 (Each processor has data for a few
    ! dimension 2 grid points, while having all data for dimension 1 and 3.)
    ! Data will be saved one dim 2 grid point at a time.
    !
    ! Function wise is the same as h5save_C3Parallel_dim2
    !
    ! All processors should call this function.
    !
    ! Arguments:
    !   filename: [string, Input] h5 filename with path
    !   varname : [string, Input] variable name
    !   myid    : [integer, Input] processor ID
    !   matrix  : [double/single complex 3d matrix, Input] data to be saved
    subroutine h5save_C3Serial_dim2( filename, varname, myid, matrix)
        use parallelization, only: pointers
        ! Inputs
        character(len=*), intent(in) :: filename, varname
        integer, intent(in) :: myid
        complex(kind=cp), intent(in), dimension(:,:,:) :: matrix

        ! data dimensions for full data and slice data
        integer, dimension(3) :: full_data_dim, slice_data_dim
        integer :: dim2slice, dim2full

        ! work assignment pointers for dimension 2
        integer :: b2v(0:numerop-1), e2v(0:numerop-1)
        integer :: b2, e2

        ! MPI variables
        integer :: ierr
        ! Loop index
        integer :: ii
        ! Data consistency check
        logical :: dataagree


        ! --------------------- Get Full Matrix Dimensions ---------------------
        ! get matrix dimension for the slices that this current processor has
        slice_data_dim = shape(matrix)
        dim2slice = slice_data_dim(2)
        ! compute the full dimension 2 size by adding all the dimension 2 slice sizes
        call MPI_ALLREDUCE(dim2slice, dim2full, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
        full_data_dim(1) = slice_data_dim(1)
        full_data_dim(2) = dim2full
        full_data_dim(3) = slice_data_dim(3)

        ! -------------- Compute the distribution in dimension 2 --------------
        call pointers( b2v, e2v, dim2full)
        b2 = b2v(myid)
        e2 = e2v(myid)

        ! --------------- Initialize variable for partial saving ---------------
        if (myid .eq. 0) then
            call h5save_CPartial_Init( filename, varname, full_data_dim)
        endif

        ! ------------- Save data to h5 and check data consistency -------------
        ! This ensures the data saving is executed at least once
        dataagree = .false.

        ! Loop until h5 file passes data consistency check
        DO WHILE ( dataagree .eqv. .false. )

            ! Loop over dim3 and save each dim 3
            DO ii = 1, dim2full
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                if ( (ii .ge. b2) .and. (ii .le. e2) ) then
                    ! if this plane is part of the data the current processor has
                    call h5save_C3Partial_SingleDim2( filename, varname, matrix(:,ii-b2+1,:), ii )
                endif
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            ENDDO

            ! verify data consistency
            dataagree = verifysave_C3_dim2(filename, varname, myid, matrix, full_data_dim)
            if ( myid .eq. 0 ) then
                if ( dataagree ) then
                    write(*,*) "   Data consistency check passed."
                else
                    write(*,*) "   Data consistency check failed, re-saving data."
                endif
            endif
        ENDDO
    end subroutine h5save_C3Serial_dim2


    ! function verifysave_C3_dim2(filename, varname, myid, matrix, full_data_dim) result( dataagree_all )
    ! Read data from the h5 file and compare to the data in memory to ensure
    ! no data corruption in the saving process
    !
    ! Arguments:
    !   filename     : [string, Input] h5 filename with path
    !   varname      : [string, Input] variable name
    !   myid         : [integer, Input] processor ID
    !   matrix       : [double/single complex 3d matrix, Input] data saved
    !   full_data_dim: [integer, size (3)] size of the full rank 3 matrix
    ! Return:
    !   dataagree_all: [logical] true if the data in the h5 file matches with what was saved
    !                            false if data corruption occured
    function verifysave_C3_dim2(filename, varname, myid, matrix, full_data_dim) result( dataagree_all )
        use parallelization, only: pointers, distribute_C3_slicedim2
        use h5load, only: h5load_C3
        ! Inputs
        character(len=*), intent(in) :: filename, varname
        integer, intent(in) :: myid
        complex(kind=cp), intent(in), dimension(:,:,:) :: matrix
        integer, intent(in), dimension(3) :: full_data_dim
        ! Output
        logical :: dataagree_local, dataagree_all

        ! work assignment pointers for dimension 3
        integer :: b2v(0:numerop-1), e2v(0:numerop-1)
        integer :: b2, e2
        ! full data read from the file
        complex(kind=cp), dimension(:,:,:), allocatable :: full_data, slice_data
        ! MPI variables
        integer :: ierr


        ! -------------- Compute the distribution in dimension 2 --------------
        call pointers( b2v, e2v, full_data_dim(2))
        b2 = b2v(myid)
        e2 = e2v(myid)

        ! ---------------- Read from h5 and distribute in dim 2 ----------------
        allocate(  full_data( full_data_dim(1),full_data_dim(2),full_data_dim(3) ) )
        allocate( slice_data( full_data_dim(1),b2:e2           ,full_data_dim(3) ) )
        ! Master reads data from the h5 file
        if (myid .eq. 0) then
            full_data = h5load_C3(filename, varname )
        endif
        ! distribute data in dimension 2
        call distribute_C3_slicedim2( myid, 0, full_data, slice_data )
        deallocate( full_data )

        ! ---------------------------- Compare data ----------------------------
        if ( all(slice_data .eq. matrix) ) then
            dataagree_local = .true.
        else
            write(*,"(A,I2.2)") "   Data corruption at processor: ", myid
            dataagree_local = .false.
        endif
        deallocate( slice_data )

        ! ------------------ Collect data from all processors ------------------
        call MPI_ALLREDUCE(dataagree_local, dataagree_all, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    end function verifysave_C3_dim2