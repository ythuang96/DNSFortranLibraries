    ! subroutine h5save_C3Parallel_dim3( filename, varname, myid, matrix )
    ! save a complex rank 3 matrix to h5 file in parallel.
    ! with parallization in dimension 3 (Each processor has data for a few
    ! dimension 3 grid points, while having all data for dimension 1 and 2.)
    !
    ! All processors should call this function.
    !
    ! Arguments:
    !   filename: [string, Input] h5 filename with path
    !   varname : [string, Input] variable name
    !   myid    : [integer, Input] processor ID
    !   matrix  : [double/single complex 3d matrix, Input] data to be saved
    subroutine h5save_C3Parallel_dim3( filename, varname, myid, matrix)
        use parallelization, only: pointers
        ! Inputs
        character(len=*), intent(in) :: filename, varname
        integer, intent(in) :: myid
        complex(kind=cp), intent(in), dimension(:,:,:) :: matrix

        real(kind=cp), dimension(:,:,:), allocatable :: temp

        character(len=100) :: dset_name ! dataset name

        ! data dimensions for full data and slice data
        integer(HSIZE_T), dimension(3) :: full_data_dim, slice_data_dim
        integer(HSIZE_T), dimension(3) :: slabOffset ! hyperslab offset
        integer :: dim3slice, dim3full

        integer :: error ! error flag
        INTEGER(HID_T) :: file_id   ! file id
        INTEGER(HID_T) :: group_id  ! group id
        INTEGER(HID_T) :: dspace_id ! dataspace id
        INTEGER(HID_T) :: mspace_id ! memeory space id
        INTEGER(HID_T) :: dset_id   ! dataset id
        INTEGER(HID_T) :: PropertyList_id ! property list for parallel access

        logical :: file_exists

        ! work assignment pointers for dimension 3
        integer :: b3v(0:numerop-1), e3v(0:numerop-1)
        integer :: b3, e3

        ! MPI variables
        integer :: ierr


        ! --------------------- Get Full Matrix Dimensions ---------------------
        ! get matrix dimension for the slices that this current processor has
        slice_data_dim = shape(matrix)
        dim3slice = slice_data_dim(3)
        ! compute the full dimension 3 size by adding all the dimension 3 slice sizes
        call MPI_ALLREDUCE(dim3slice, dim3full, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
        full_data_dim(1) = slice_data_dim(1)
        full_data_dim(2) = slice_data_dim(2)
        full_data_dim(3) = INT( dim3full, HID_T)

        call pointers( b3v, e3v, dim3full)
        b3 = b3v(myid)
        e3 = e3v(myid)

        ! Allocate temp buffer
        allocate(temp( slice_data_dim(1), slice_data_dim(2), slice_data_dim(3)))

        ! ------------------------ Setup File and Group ------------------------
        ! Initialize hdf5 interface
        call h5open_f(error)
        ! Create h5 property list for parallel access
        call h5pcreate_f(H5P_FILE_ACCESS_F, PropertyList_id, error)
        call h5pset_fapl_mpio_f(PropertyList_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
        ! Check if file exist
        INQUIRE(FILE=filename, EXIST=file_exists)
        if ( file_exists ) then
            ! If file already exist, then open file with read and write access
            ! With parallel access property list
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error, access_prp=PropertyList_id)
        else
            ! Create new file
            ! With parallel access property list
            CALL h5fcreate_f(filename, H5F_ACC_EXCL_F, file_id, error, access_prp=PropertyList_id)
        endif
        ! Close Property list
        call h5pclose_f(PropertyList_id, error)

        ! Create a group in the HDF5 file with the variable name
        CALL h5gcreate_f(file_id, varname, group_id, error)
        ! Close the group
        CALL h5gclose_f(group_id, error)

        ! ----------------------------- Real part -----------------------------
        dset_name = varname // "/" // varname // "_REAL"

        ! Create disk dataspace with rank 3 and size full_data_dim
        call h5screate_simple_f(3, full_data_dim, dspace_id, error)
        ! select that hyperlsab of the disk dataspace, which is the xz planes of this processor (offset b3-1)
        slabOffset = (/ 0, 0, b3 - 1 /)
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)
        ! Create memory dataspace with rank 3 and size slice_data_dim
        call h5screate_simple_f(3, slice_data_dim, mspace_id, error)
        ! select that hyperlsab of the memory dataspace, (no offset)
        slabOffset = (/ 0, 0, 0 /)
        call h5sselect_hyperslab_f(mspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)

        ! Create h5 property list for parallel data transfer
        call h5pcreate_f(H5P_DATASET_XFER_F, PropertyList_id, error)
        call h5pset_dxpl_mpio_f(PropertyList_id, H5FD_MPIO_COLLECTIVE_F, error)
        ! Create double/single precision dataset with path '/var/var_REAL' and write data
        temp = real(matrix, cp)
        if ( cp .eq. dp ) then
            call h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
            call h5dwrite_f(dset_id, H5T_IEEE_F64LE, temp, slice_data_dim, error, &
                    mem_space_id=mspace_id, file_space_id=dspace_id, xfer_prp=PropertyList_id)
        else if ( cp .eq. sp ) then
            call h5dcreate_f(file_id, dset_name, H5T_IEEE_F32LE, dspace_id, dset_id, error)
            call h5dwrite_f(dset_id, H5T_IEEE_F32LE, temp, slice_data_dim, error, &
                    mem_space_id=mspace_id, file_space_id=dspace_id, xfer_prp=PropertyList_id)
        endif
        ! Close property list
        call h5pclose_f(PropertyList_id, error)
        ! Close dataset
        call h5dclose_f(dset_id, error)

        ! Close disk dataspace
        call h5sclose_f(dspace_id, error)
        ! Close memory dataspace
        call h5sclose_f(mspace_id, error)
        ! Flush buffers to disk. Not flushing buffers may lead to corrupted hdf5 files
        call h5fflush_f(file_id, H5F_SCOPE_GLOBAL_F, error)

        ! --------------------------- Imaginary part ---------------------------
        dset_name = varname // "/" // varname // "_IMAG"

        ! Create disk dataspace with rank 3 and size full_data_dim
        call h5screate_simple_f(3, full_data_dim, dspace_id, error)
        ! select that hyperlsab of the disk dataspace, which is the xz planes of this processor (offset b3-1)
        slabOffset = (/ 0, 0, b3 - 1 /)
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)
        ! Create memory dataspace with rank 3 and size slice_data_dim
        call h5screate_simple_f(3, slice_data_dim, mspace_id, error)
        ! select that hyperlsab of the memory dataspace, (no offset)
        slabOffset = (/ 0, 0, 0 /)
        call h5sselect_hyperslab_f(mspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)

        ! Create h5 property list for parallel data transfer
        call h5pcreate_f(H5P_DATASET_XFER_F, PropertyList_id, error)
        call h5pset_dxpl_mpio_f(PropertyList_id, H5FD_MPIO_COLLECTIVE_F, error)
        ! Create double/single precision dataset with path '/var/var_REAL' and write data
        temp = aimag(matrix)
        if ( cp .eq. dp ) then
            call h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
            call h5dwrite_f(dset_id, H5T_IEEE_F64LE, temp, slice_data_dim, error, &
                    mem_space_id=mspace_id, file_space_id=dspace_id, xfer_prp=PropertyList_id)
        else if ( cp .eq. sp ) then
            call h5dcreate_f(file_id, dset_name, H5T_IEEE_F32LE, dspace_id, dset_id, error)
            call h5dwrite_f(dset_id, H5T_IEEE_F32LE, temp, slice_data_dim, error, &
                    mem_space_id=mspace_id, file_space_id=dspace_id, xfer_prp=PropertyList_id)
        endif
        ! Close property list
        call h5pclose_f(PropertyList_id, error)
        ! Close dataset
        call h5dclose_f(dset_id, error)

        ! Close disk dataspace
        call h5sclose_f(dspace_id, error)
        ! Close memory dataspace
        call h5sclose_f(mspace_id, error)
        ! Flush buffers to disk. Not flushing buffers may lead to corrupted hdf5 files
        call h5fflush_f(file_id, H5F_SCOPE_GLOBAL_F, error)

        ! ------------------------------ Clean up ------------------------------
        ! Close the file
        CALL h5fclose_f(file_id, error)
        ! Close FORTRAN interface
        CALL h5close_f(error)
        ! Deallocate temp buffer
        DEALLOCATE(temp)
    end subroutine h5save_C3Parallel_dim3


    ! subroutine h5save_C3Parallel_dim2( filename, varname, myid, matrix )
    ! save a complex rank 3 matrix to h5 file in parallel.
    ! with parallization in dimension 2 (Each processor has data for a few
    ! dimension 2 grid points, while having all data for dimension 1 and 3.)
    !
    ! All processors should call this function.
    !
    ! Arguments:
    !   filename: [string, Input] h5 filename with path
    !   varname : [string, Input] variable name
    !   myid    : [integer, Input] processor ID
    !   matrix  : [double/single complex 3d matrix, Input] data to be saved
    subroutine h5save_C3Parallel_dim2( filename, varname, myid, matrix)
        use parallelization, only: pointers
        ! Inputs
        character(len=*), intent(in) :: filename, varname
        integer, intent(in) :: myid
        complex(kind=cp), intent(in), dimension(:,:,:) :: matrix

        real(kind=cp), dimension(:,:,:), allocatable :: temp

        character(len=100) :: dset_name ! dataset name

        ! data dimensions for full data and slice data
        integer(HSIZE_T), dimension(3) :: full_data_dim, slice_data_dim
        integer(HSIZE_T), dimension(3) :: slabOffset ! hyperslab offset
        integer :: dim2slice, dim2full

        integer :: error ! error flag
        INTEGER(HID_T) :: file_id   ! file id
        INTEGER(HID_T) :: group_id  ! group id
        INTEGER(HID_T) :: dspace_id ! dataspace id
        INTEGER(HID_T) :: mspace_id ! memeory space id
        INTEGER(HID_T) :: dset_id   ! dataset id
        INTEGER(HID_T) :: PropertyList_id ! property list for parallel access

        logical :: file_exists

        ! work assignment pointers for dimension 2
        integer :: b2v(0:numerop-1), e2v(0:numerop-1)
        integer :: b2, e2

        ! MPI variables
        integer :: ierr


        ! --------------------- Get Full Matrix Dimensions ---------------------
        ! get matrix dimension for the slices that this current processor has
        slice_data_dim = shape(matrix)
        dim2slice = slice_data_dim(2)
        ! compute the full dimension 2 size by adding all the dimension 2 slice sizes
        call MPI_ALLREDUCE(dim2slice, dim2full, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
        full_data_dim(1) = slice_data_dim(1)
        full_data_dim(2) = INT( dim2full, HID_T)
        full_data_dim(3) = slice_data_dim(3)

        call pointers( b2v, e2v, dim2full)
        b2 = b2v(myid)
        e2 = e2v(myid)

        ! Allocate temp buffer
        allocate(temp( slice_data_dim(1), slice_data_dim(2), slice_data_dim(3)))

        ! ------------------------ Setup File and Group ------------------------
        ! Initialize hdf5 interface
        call h5open_f(error)
        ! Create h5 property list for parallel access
        call h5pcreate_f(H5P_FILE_ACCESS_F, PropertyList_id, error)
        call h5pset_fapl_mpio_f(PropertyList_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
        ! Check if file exist
        INQUIRE(FILE=filename, EXIST=file_exists)
        if ( file_exists ) then
            ! If file already exist, then open file with read and write access
            ! With parallel access property list
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error, access_prp=PropertyList_id)
        else
            ! Create new file
            ! With parallel access property list
            CALL h5fcreate_f(filename, H5F_ACC_EXCL_F, file_id, error, access_prp=PropertyList_id)
        endif
        ! Close Property list
        call h5pclose_f(PropertyList_id, error)

        ! Create a group in the HDF5 file with the variable name
        CALL h5gcreate_f(file_id, varname, group_id, error)
        ! Close the group
        CALL h5gclose_f(group_id, error)

        ! ----------------------------- Real part -----------------------------
        dset_name = varname // "/" // varname // "_REAL"

        ! Create disk dataspace with rank 3 and size full_data_dim
        call h5screate_simple_f(3, full_data_dim, dspace_id, error)
        ! select that hyperlsab of the disk dataspace, which is the xz planes of this processor (offset b2-1)
        slabOffset = (/ 0, b2 - 1, 0 /)
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)
        ! Create memory dataspace with rank 3 and size slice_data_dim
        call h5screate_simple_f(3, slice_data_dim, mspace_id, error)
        ! select that hyperlsab of the memory dataspace, (no offset)
        slabOffset = (/ 0, 0, 0 /)
        call h5sselect_hyperslab_f(mspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)

        ! Create h5 property list for parallel data transfer
        call h5pcreate_f(H5P_DATASET_XFER_F, PropertyList_id, error)
        call h5pset_dxpl_mpio_f(PropertyList_id, H5FD_MPIO_COLLECTIVE_F, error)
        ! Create double/single precision dataset with path '/var/var_REAL' and write data
        temp = real(matrix, cp)
        if ( cp .eq. dp ) then
            call h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
            call h5dwrite_f(dset_id, H5T_IEEE_F64LE, temp, slice_data_dim, error, &
                    mem_space_id=mspace_id, file_space_id=dspace_id, xfer_prp=PropertyList_id)
        else if ( cp .eq. sp ) then
            call h5dcreate_f(file_id, dset_name, H5T_IEEE_F32LE, dspace_id, dset_id, error)
            call h5dwrite_f(dset_id, H5T_IEEE_F32LE, temp, slice_data_dim, error, &
                    mem_space_id=mspace_id, file_space_id=dspace_id, xfer_prp=PropertyList_id)
        endif
        ! Close property list
        call h5pclose_f(PropertyList_id, error)
        ! Close dataset
        call h5dclose_f(dset_id, error)

        ! Close disk dataspace
        call h5sclose_f(dspace_id, error)
        ! Close memory dataspace
        call h5sclose_f(mspace_id, error)
        ! Flush buffers to disk. Not flushing buffers may lead to corrupted hdf5 files
        call h5fflush_f(file_id, H5F_SCOPE_GLOBAL_F, error)

        ! --------------------------- Imaginary part ---------------------------
        dset_name = varname // "/" // varname // "_IMAG"

        ! Create disk dataspace with rank 3 and size full_data_dim
        call h5screate_simple_f(3, full_data_dim, dspace_id, error)
        ! select that hyperlsab of the disk dataspace, which is the xz planes of this processor (offset b2-1)
        slabOffset = (/ 0, b2 - 1, 0 /)
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)
        ! Create memory dataspace with rank 3 and size slice_data_dim
        call h5screate_simple_f(3, slice_data_dim, mspace_id, error)
        ! select that hyperlsab of the memory dataspace, (no offset)
        slabOffset = (/ 0, 0, 0 /)
        call h5sselect_hyperslab_f(mspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)

        ! Create h5 property list for parallel data transfer
        call h5pcreate_f(H5P_DATASET_XFER_F, PropertyList_id, error)
        call h5pset_dxpl_mpio_f(PropertyList_id, H5FD_MPIO_COLLECTIVE_F, error)
        ! Create double/single precision dataset with path '/var/var_REAL' and write data
        temp = aimag(matrix)
        if ( cp .eq. dp ) then
            call h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
            call h5dwrite_f(dset_id, H5T_IEEE_F64LE, temp, slice_data_dim, error, &
                    mem_space_id=mspace_id, file_space_id=dspace_id, xfer_prp=PropertyList_id)
        else if ( cp .eq. sp ) then
            call h5dcreate_f(file_id, dset_name, H5T_IEEE_F32LE, dspace_id, dset_id, error)
            call h5dwrite_f(dset_id, H5T_IEEE_F32LE, temp, slice_data_dim, error, &
                    mem_space_id=mspace_id, file_space_id=dspace_id, xfer_prp=PropertyList_id)
        endif
        ! Close property list
        call h5pclose_f(PropertyList_id, error)
        ! Close dataset
        call h5dclose_f(dset_id, error)

        ! Close disk dataspace
        call h5sclose_f(dspace_id, error)
        ! Close memory dataspace
        call h5sclose_f(mspace_id, error)
        ! Flush buffers to disk. Not flushing buffers may lead to corrupted hdf5 files
        call h5fflush_f(file_id, H5F_SCOPE_GLOBAL_F, error)

        ! ------------------------------ Clean up ------------------------------
        ! Close the file
        CALL h5fclose_f(file_id, error)
        ! Close FORTRAN interface
        CALL h5close_f(error)
        ! Deallocate temp buffer
        DEALLOCATE(temp)
    end subroutine h5save_C3Parallel_dim2