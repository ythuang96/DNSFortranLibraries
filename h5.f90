module h5
    use types, only: sp, dp
    use mpi
    use hdf5
    implicit none
    private
    include 'parameters'

    public :: h5load_real_vector, h5readVelocities, h5save_complex_2d_matrix
    public :: h5save_real_2d_matrix, h5save_complex_3d_matrix_parallel


contains
    ! function ouput = h5load_real_vector
    ! Arguments:
    !   filename: [string, Input]
    !             h5 filename with path
    !   varname : [string, Input]
    !             variable name in h5 file, must be real numerical vector
    ! Output:
    !   output:   [double vector]
    function h5load_real_vector(filename, varname ) result(vector)
        character(len=*), intent(in) :: filename, varname
        real(kind=dp), dimension(:,:), allocatable :: matrix
        real(kind=dp), dimension(:), allocatable :: vector

        INTEGER(HID_T) :: file_id        ! File identifier
        INTEGER(HID_T) :: dset_id        ! Dataset identifier
        INTEGER(HID_T) :: space_id       ! Dataspace identifier
        INTEGER(HID_T) :: dtype_id       ! Dataspace identifier

        INTEGER :: error ! Error flag
        INTEGER :: dim1, dim2 ! matrix dimensions

        INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
        INTEGER(HSIZE_T), DIMENSION(2) :: max_dims


        ! Initialize FORTRAN interface.
        CALL h5open_f(error)
        ! Open an existing file with read only
        CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)
        ! Open an existing dataset.
        CALL h5dopen_f(file_id, varname, dset_id, error)
        !Get dataspace ID
        CALL h5dget_space_f(dset_id, space_id,error)
        !Get dataspace dims
        CALL h5sget_simple_extent_dims_f(space_id, data_dims, max_dims, error)
        dim1 = data_dims(1)
        dim2 = data_dims(2)

        ! Allocate dimensions to dset_data for reading
        ALLOCATE(matrix(dim1,dim2))
        ! Get data
        ! H5T_IEEE_F64LE (double) or H5T_IEEE_F32LE (single) has to be
        ! consistent with the variable type of matrix
        CALL h5dread_f(dset_id, H5T_IEEE_F64LE, matrix, data_dims, error)

        ! Reshape matrix into a vector
        vector = reshape( matrix, (/ dim1 /) )

        ! close dataset
        CALL h5dclose_f(dset_id, error)
        ! close file
        CALL h5fclose_f(file_id, error)
        ! close h5 interface
        CALL h5close_f(error)
    end function h5load_real_vector

    ! subroutine h5readVelocities( filename, varname, Buffer )
    ! Arguments:
    !   filename: [string, Input]
    !             h5 filename with path
    !   varname:  [string, Input]
    !             "u" or "v" or "w" indicating which velocity component to read
    !   Buffer:   [3d double matrix, size (mxf,mzf,myf), Output]
    !             velocity read from the h5 file
    subroutine h5readVelocities( filename, varname, Buffer )
        character(len=*), intent(in) :: filename
        character(len=1), intent(in) :: varname
        complex(kind=dp), intent(out), dimension(mxf,mzf,myf) :: Buffer
        ! temp arrays for real and imaginary parts
        real(kind=sp), dimension(mxf,mzf,myf) :: rtemp, itemp

        INTEGER(HID_T) :: file_id        ! File identifier
        INTEGER(HID_T) :: dset_id        ! Dataset identifier

        INTEGER :: error ! Error flag

        INTEGER(HSIZE_T), DIMENSION(3) :: data_dims ! data dimensions


        ! data dimensions
        data_dims = (/ mxf,mzf,myf /)

        ! Initialize FORTRAN interface.
        CALL h5open_f(error)
        ! Open an existing file with read only
        CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)
        if ( error .eq. -1) then
            write(*,*) "Failed to open file: ", filename
            write(*,*)
            stop
        endif

        ! Open an existing dataset and get data
        ! H5T_IEEE_F64LE (double) or H5T_IEEE_F32LE (single) has to be
        ! consistent with the variable type of output
        ! ----------------------------- Real part -----------------------------
        ! Open dataset
        CALL h5dopen_f(file_id, ("/velocityFieldsFourier/" // varname // "RealPart"), dset_id, error)
        ! Read dataset
        CALL h5dread_f(dset_id, H5T_IEEE_F32LE, rtemp, data_dims, error)
        ! Close dataset
        CALL h5dclose_f(dset_id, error)
        ! --------------------------- Imaginary part ---------------------------
        ! Open dataset
        CALL h5dopen_f(file_id, ("/velocityFieldsFourier/" // varname // "ImaginaryPart"), dset_id, error)
        ! Read dataset
        CALL h5dread_f(dset_id, H5T_IEEE_F32LE, itemp, data_dims, error)
        ! Close dataset
        CALL h5dclose_f(dset_id, error)


        ! Build data from real and imaginary part
        Buffer = CMPLX( DBLE(rtemp), DBLE(itemp) , dp )

        ! close file
        CALL h5fclose_f(file_id, error)
        ! close h5 interface
        CALL h5close_f(error)
    end subroutine h5readVelocities


    ! subroutine h5save_complex_2d_matrix( filename, varname, matrix )
    ! Arguments:
    !   filename: [string, Input] h5 filename with path
    !   varname : [string, Input] variable name (complex 2d matrix)
    !   matrix  : [double complex 2d matrix, Input] data to be saved
    ! The following standard for complex variables are used:
    ! if varibale 'var' is complex, save as '/var/var_REAL' and '/var/var_IMAG'
    ! (one group with two datasets)
    ! Same as my matlab h5 libaries
    subroutine h5save_complex_2d_matrix( filename, varname, matrix )
        character(len=*), intent(in) :: filename, varname
        complex(kind=dp), intent(in), dimension(:,:) :: matrix

        character(len=100) :: dset_name ! dataset name
        integer(HSIZE_T), dimension(2) :: data_dim ! data dimensions

        integer :: error ! error flag
        INTEGER(HID_T) :: file_id  ! file id
        INTEGER(HID_T) :: group_id ! group id
        INTEGER(HID_T) :: dspace_id ! dataspace id
        INTEGER(HID_T) :: dset_id ! dataset id

        logical :: file_exists

        ! get matrix dimensions
        data_dim = shape(matrix)


        ! ------------------------ Setup File and Group ------------------------
        ! Initialize hdf5 interface
        call h5open_f(error)
        ! Check if file exist
        INQUIRE(FILE=filename, EXIST=file_exists)
        if ( file_exists ) then
            ! If file already exist, then open file with read and write access
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
        else
            ! Create new file
            CALL h5fcreate_f(filename, H5F_ACC_EXCL_F, file_id, error)
        endif

        ! Create a group in the HDF5 file with the variable name
        CALL h5gcreate_f(file_id, varname, group_id, error)
        ! Close the group
        CALL h5gclose_f(group_id, error)

        ! ----------------------------- Real part -----------------------------
        dset_name = varname // "/" // varname // "_REAL"
        ! Create dataspace with rank 2 and size data_dim
        CALL h5screate_simple_f(2, data_dim, dspace_id, error)
        ! Create double precision dataset with path '/var/var_REAL'
        CALL h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
        ! Write dataset
        CALL h5dwrite_f(dset_id, H5T_IEEE_F64LE, real(matrix, dp), data_dim, error)
        ! Close dataset
        CALL h5dclose_f(dset_id, error)
        ! Close dataspace
        CALL h5sclose_f(dspace_id, error)

        ! --------------------------- Imaginary part ---------------------------
        dset_name = varname // "/" // varname // "_IMAG"
        ! Create dataspace with rank 2 and size data_dim
        CALL h5screate_simple_f(2, data_dim, dspace_id, error)
        ! Create double precision dataset with path '/var/var_IMAG'
        CALL h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
        ! Write dataset
        CALL h5dwrite_f(dset_id, H5T_IEEE_F64LE, aimag(matrix), data_dim, error)
        ! Close dataset
        CALL h5dclose_f(dset_id, error)
        ! Close dataspace
        CALL h5sclose_f(dspace_id, error)

        ! ------------------------------ Clean up ------------------------------
        ! Close the file
        CALL h5fclose_f(file_id, error)
        ! Close FORTRAN interface
        CALL h5close_f(error)
    end subroutine h5save_complex_2d_matrix


    ! subroutine h5save_real_2d_matrix( filename, varname, matrix )
    ! Arguments:
    !   filename: [string, Input] h5 filename with path
    !   varname : [string, Input] variable name (real 2d matrix)
    !   matrix  : [double 2d matrix, Input] data to be saved
    ! The following standard for complex variables are used:
    ! if varibale 'var' is complex, save as '/var/var_REAL' and '/var/var_IMAG'
    ! (one group with two datasets)
    ! Same as my matlab h5 libaries
    subroutine h5save_real_2d_matrix( filename, varname, matrix )
        character(len=*), intent(in) :: filename, varname
        real(kind=dp), intent(in), dimension(:,:) :: matrix

        character(len=100) :: dset_name ! dataset name
        integer(HSIZE_T), dimension(2) :: data_dim ! data dimensions

        integer :: error ! error flag
        INTEGER(HID_T) :: file_id  ! file id
        INTEGER(HID_T) :: group_id ! group id
        INTEGER(HID_T) :: dspace_id ! dataspace id
        INTEGER(HID_T) :: dset_id ! dataset id

        logical :: file_exists

        ! get matrix dimensions
        data_dim = shape(matrix)


        ! Initialize hdf5 interface
        call h5open_f(error)
        ! Check if file exist
        INQUIRE(FILE=filename, EXIST=file_exists)
        if ( file_exists ) then
            ! If file already exist, then open file with read and write access
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
        else
            ! Create new file
            CALL h5fcreate_f(filename, H5F_ACC_EXCL_F, file_id, error)
        endif


        dset_name = varname
        ! Create dataspace with rank 2 and size data_dim
        CALL h5screate_simple_f(2, data_dim, dspace_id, error)
        ! Create double precision dataset with path '/var'
        CALL h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
        ! Write dataset
        CALL h5dwrite_f(dset_id, H5T_IEEE_F64LE, matrix, data_dim, error)
        ! Close dataset
        CALL h5dclose_f(dset_id, error)
        ! Close dataspace
        CALL h5sclose_f(dspace_id, error)


        ! Close the file
        CALL h5fclose_f(file_id, error)
        ! Close FORTRAN interface
        CALL h5close_f(error)
    end subroutine h5save_real_2d_matrix


    ! subroutine h5save_complex_3d_matrix_parallel( filename, varname, matrix )
    ! This function saves a complex 3d matrix in parallel. Each processor has
    ! data for a few wall parallel planes. And the full data has myf as the
    ! full third dimension. All processors should call this function.
    !
    ! Arguments:
    !   filename: [string, Input] h5 filename with path
    !   varname : [string, Input] variable name
    !   matrix  : [double complex 3d matrix, Input] data to be saved
    ! The following standard for complex variables are used:
    ! if varibale 'var' is complex, save as '/var/var_REAL' and '/var/var_IMAG'
    ! (one group with two datasets)
    ! Same as my matlab h5 libaries
    subroutine h5save_complex_3d_matrix_parallel( filename, varname, matrix)
        ! Global variables
        integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
        common /point/ jbeg(0:numerop-1),jend(0:numerop-1), &
                       kbeg(0:numerop-1),kend(0:numerop-1), &
                       jb,je,kb,ke,mmy,mmz
        save /point/
        ! Inputs
        character(len=*), intent(in) :: filename, varname
        complex(kind=dp), intent(in), dimension(:,:,:) :: matrix

        character(len=100) :: dset_name ! dataset name

        ! data dimensions for full data and slice data
        integer(HSIZE_T), dimension(3) :: full_data_dim, slice_data_dim
        integer(HSIZE_T), dimension(3) :: slabOffset ! hyperslab offset

        integer :: error ! error flag
        INTEGER(HID_T) :: file_id   ! file id
        INTEGER(HID_T) :: group_id  ! group id
        INTEGER(HID_T) :: dspace_id ! dataspace id
        INTEGER(HID_T) :: mspace_id ! memeory space id
        INTEGER(HID_T) :: dset_id   ! dataset id
        INTEGER(HID_T) :: PropertyList_id ! property list for parallel access

        logical :: file_exists

        ! get full matrix dimensions
        full_data_dim = shape(matrix)
        full_data_dim(3) = myf ! the third dimension should always be the y grid size
        ! get matrix dimension for the slices that this current processor has
        slice_data_dim = shape(matrix)
        slice_data_dim(3) = je - jb + 1


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
        ! select that hyperlsab of the disk dataspace, which is the xz planes of this processor (offset jb-1)
        slabOffset = (/ 0, 0, jb - 1 /)
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)
        ! Create memory dataspace with rank 3 and size slice_data_dim
        call h5screate_simple_f(3, slice_data_dim, mspace_id, error)
        ! select that hyperlsab of the memory dataspace, (no offset)
        slabOffset = (/ 0, 0, 0 /)
        call h5sselect_hyperslab_f(mspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)

        ! Create double precision dataset with path '/var/var_REAL'
        call h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
        ! Create h5 property list for parallel data transfer
        call h5pcreate_f(H5P_DATASET_XFER_F, PropertyList_id, error)
        call h5pset_dxpl_mpio_f(PropertyList_id, H5FD_MPIO_COLLECTIVE_F, error)
        ! Write dataset
        call h5dwrite_f(dset_id, H5T_IEEE_F64LE, real(matrix, dp), slice_data_dim, error, &
                mem_space_id=mspace_id, file_space_id=dspace_id, xfer_prp=PropertyList_id)
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
        ! select that hyperlsab of the disk dataspace, which is the xz planes of this processor (offset jb-1)
        slabOffset = (/ 0, 0, jb - 1 /)
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)
        ! Create memory dataspace with rank 3 and size slice_data_dim
        call h5screate_simple_f(3, slice_data_dim, mspace_id, error)
        ! select that hyperlsab of the memory dataspace, (no offset)
        slabOffset = (/ 0, 0, 0 /)
        call h5sselect_hyperslab_f(mspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)

        ! Create double precision dataset with path '/var/var_REAL'
        call h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
        ! Create h5 property list for parallel data transfer
        call h5pcreate_f(H5P_DATASET_XFER_F, PropertyList_id, error)
        call h5pset_dxpl_mpio_f(PropertyList_id, H5FD_MPIO_COLLECTIVE_F, error)
        ! Write dataset
        call h5dwrite_f(dset_id, H5T_IEEE_F64LE, aimag(matrix), slice_data_dim, error, &
                mem_space_id=mspace_id, file_space_id=dspace_id, xfer_prp=PropertyList_id)
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

    end subroutine h5save_complex_3d_matrix_parallel

end module h5

