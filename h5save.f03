module h5save
    use types, only: sp, dp
    use mpi
    use hdf5
    implicit none
    private
#   include "parameters"

    public :: check_filename, h5save_logical, h5save_string, h5save_R, h5save_R1, h5save_C2, h5save_R2
    public :: h5save_C3Partial_Init, h5save_C3Partial_SingleDim3, h5save_C3Partial_SingleDim2
    public :: h5save_C3Serial_dim3, h5save_C3Serial_dim2

    ! The following standard for complex variables are used:
    ! if varibale 'var' is complex, save as '/var/var_REAL' and '/var/var_IMAG'
    ! (one group with two datasets)
    ! Same as my matlab h5 libaries
contains
    ! subroutine check_filename( filename, myid )
    ! This function check if the file with filename already exist
    ! if it exist, will append the current date and time to the filename to
    ! to create a new filename and send to all processors
    ! Arguments:
    !   filename: [string, Input/Output] h5 filename with path, it is appended
    !                                    with date and time if file already exist
    !   myid    : [integer, Input] processor ID
    subroutine check_filename( filename, myid )
        character(len=*), intent(inout) :: filename
        integer, intent(in) :: myid

        logical :: file_exists

        integer :: values(8)
        character(len=10) :: current_data_and_time

        integer :: loc, count

        integer :: ierr


        if ( myid .eq. 0 ) then
            ! master checks if file exist
            INQUIRE(FILE=filename, EXIST=file_exists)

            if ( file_exists ) then
                ! get current date and time
                call date_and_time(VALUES=values)
                write(current_data_and_time,"(5I2.2)") values(2), values(3), values(5), values(6), values(7)

                ! append date and time to current filename if file already exist
                loc = INDEX( filename, ".h5" )
                filename = filename(1:loc-1) // '_' // current_data_and_time // '.h5'
            endif
        endif

        ! send filename to all processors regardless of whether it is appended
        ! with the date and time or not
        count = LEN( filename )
        call MPI_BCAST(filename, count, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

        ! write to output file
        if ( myid .eq. 0 ) then
            write(*,*) "  Saving data to file: ", filename
        endif

    end subroutine check_filename


    ! subroutine h5save_logical( filename, varname, logicvar )
    ! save a logical to h5 file by saving as a "true" or "false" string
    ! Arguments:
    !   filename: [string, Input] h5 filename with path
    !   varname : [string, Input] variable name saved in h5
    !   logicvar: [logical, Input] logical variable to be saved
    subroutine h5save_logical( filename, varname, logicvar )
        character(len=*), intent(in) :: filename, varname
        logical, intent(in) :: logicvar

        if ( logicvar .eqv. .true. ) then
            call h5save_string( filename, varname, "true"  )
        else if ( logicvar .eqv. .false. ) then
            call h5save_string( filename, varname, "false" )
        endif

    end subroutine h5save_logical


    ! subroutine h5save_string( filename, varname, string )
    ! save a string to h5 file
    ! Arguments:
    !   filename: [string, Input] h5 filename with path
    !   varname : [string, Input] variable name saved in h5
    !   string  : [string, Input] single string to be saved
    subroutine h5save_string( filename, varname, string )
        character(len=*), intent(in) :: filename, varname
        character(len=*), intent(in) :: string

        character(len=100) :: dset_name ! dataset name

        INTEGER(SIZE_T)  :: slength ! string length
        integer(HSIZE_T), dimension(1) :: data_dim ! data dimensions

        integer :: error ! error flag
        INTEGER(HID_T) :: file_id   ! file id
        INTEGER(HID_T) :: dspace_id ! dataspace id
        INTEGER(HID_T) :: dset_id   ! dataset id
        INTEGER(HID_T) :: filetype  ! custom type for the fixed length string

        logical :: file_exists


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

        ! Get string length
        slength = LEN( string )
        ! Create string type with length slenght from the single length tyupe H5T_FORTRAN_S1
        CALL H5TCOPY_F(H5T_FORTRAN_S1,filetype, error)
        CALL H5TSET_SIZE_F(filetype, slength, error)

        ! String array of size 1 (single string)
        data_dim = 1

        dset_name = varname
        ! Create dataspace with rank 0 and size 1 (single string)
        CALL h5screate_simple_f(0, data_dim, dspace_id, error)
        ! Create string dataset with path '/var' and write data
        CALL h5dcreate_f(file_id, dset_name, filetype, dspace_id, dset_id, error)
        CALL h5dwrite_f(dset_id,  filetype, string, data_dim, error)
        ! Close dataset
        CALL h5dclose_f(dset_id, error)
        ! Close dataspace
        CALL h5sclose_f(dspace_id, error)

        ! Close the file
        CALL h5fclose_f(file_id, error)
        ! Close FORTRAN interface
        CALL h5close_f(error)
    end subroutine h5save_string


    ! subroutine h5save_R( filename, varname, scalar )
    ! save a real number to h5 file
    ! Arguments:
    !   filename: [string, Input] h5 filename with path
    !   varname : [string, Input] variable name saved in h5
    !   scalar  : [double/single scalar, Input] data to be saved
    subroutine h5save_R( filename, varname, scalar )
        character(len=*), intent(in) :: filename, varname
        real(kind=cp), intent(in) :: scalar

        character(len=100) :: dset_name ! dataset name
        integer(HSIZE_T), dimension(1) :: data_dim ! data dimensions

        integer :: error ! error flag
        INTEGER(HID_T) :: file_id  ! file id
        INTEGER(HID_T) :: dspace_id ! dataspace id
        INTEGER(HID_T) :: dset_id ! dataset id

        logical :: file_exists


        data_dim = 1

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
        ! Create dataspace with rank 0 and size 1
        CALL h5screate_simple_f(0, data_dim, dspace_id, error)
        ! Create double/single precision dataset with path '/var' and write data
        if ( cp .eq. dp ) then
            CALL h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
            CALL h5dwrite_f(dset_id, H5T_IEEE_F64LE, scalar, data_dim, error)
        else if ( cp .eq. sp ) then
            CALL h5dcreate_f(file_id, dset_name, H5T_IEEE_F32LE, dspace_id, dset_id, error)
            CALL h5dwrite_f(dset_id, H5T_IEEE_F32LE, scalar, data_dim, error)
        endif
        ! Close dataset
        CALL h5dclose_f(dset_id, error)
        ! Close dataspace
        CALL h5sclose_f(dspace_id, error)

        ! Close the file
        CALL h5fclose_f(file_id, error)
        ! Close FORTRAN interface
        CALL h5close_f(error)
    end subroutine h5save_R


    ! subroutine h5save_C2( filename, varname, matrix )
    ! save a complex rank 2 matrix to h5 file
    ! Arguments:
    !   filename: [string, Input] h5 filename with path
    !   varname : [string, Input] variable name (complex 2d matrix)
    !   matrix  : [double/single complex 2d matrix, Input] data to be saved
    subroutine h5save_C2( filename, varname, matrix )
        character(len=*), intent(in) :: filename, varname
        complex(kind=cp), intent(in), dimension(:,:) :: matrix

        real(kind=cp), dimension(:,:), allocatable :: temp

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
        ! Allocate temp buffer
        ALLOCATE(temp(data_dim(1), data_dim(2)))

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
        ! Create double/single precision dataset with path '/var/var_REAL' and write data
        temp = real(matrix, cp)
        if ( cp .eq. dp ) then
            CALL h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
            CALL h5dwrite_f(dset_id, H5T_IEEE_F64LE, temp, data_dim, error)
        else if ( cp .eq. sp ) then
            CALL h5dcreate_f(file_id, dset_name, H5T_IEEE_F32LE, dspace_id, dset_id, error)
            CALL h5dwrite_f(dset_id, H5T_IEEE_F32LE, temp, data_dim, error)
        endif
        ! Close dataset
        CALL h5dclose_f(dset_id, error)
        ! Close dataspace
        CALL h5sclose_f(dspace_id, error)

        ! --------------------------- Imaginary part ---------------------------
        dset_name = varname // "/" // varname // "_IMAG"
        ! Create dataspace with rank 2 and size data_dim
        CALL h5screate_simple_f(2, data_dim, dspace_id, error)
        ! Create double/single precision dataset with path '/var/var_IMAG' and write data
        temp = aimag(matrix)
        if ( cp .eq. dp ) then
            CALL h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
            CALL h5dwrite_f(dset_id, H5T_IEEE_F64LE, temp, data_dim, error)
        else if ( cp .eq. sp ) then
            CALL h5dcreate_f(file_id, dset_name, H5T_IEEE_F32LE, dspace_id, dset_id, error)
            CALL h5dwrite_f(dset_id, H5T_IEEE_F32LE, temp, data_dim, error)
        endif
        ! Close dataset
        CALL h5dclose_f(dset_id, error)
        ! Close dataspace
        CALL h5sclose_f(dspace_id, error)

        ! ------------------------------ Clean up ------------------------------
        ! Close the file
        CALL h5fclose_f(file_id, error)
        ! Close FORTRAN interface
        CALL h5close_f(error)
        ! Deallocate temp buffer
        DEALLOCATE(temp)
    end subroutine h5save_C2


    ! subroutine h5save_R1( filename, varname, vector )
    ! save a real rank 1 vector to h5 file
    ! Arguments:
    !   filename: [string, Input] h5 filename with path
    !   varname : [string, Input] variable name (real 1d vector)
    !   vector  : [double/single 1d vector, Input] data to be saved
    subroutine h5save_R1( filename, varname, vector )
        character(len=*), intent(in) :: filename, varname
        real(kind=cp), intent(in), dimension(:) :: vector

        character(len=100) :: dset_name ! dataset name
        integer(HSIZE_T), dimension(1) :: data_dim ! data dimensions

        integer :: error ! error flag
        INTEGER(HID_T) :: file_id  ! file id
        INTEGER(HID_T) :: dspace_id ! dataspace id
        INTEGER(HID_T) :: dset_id ! dataset id

        logical :: file_exists


        ! get vector dimensions
        data_dim = shape(vector)

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
        ! Create dataspace with rank 1 and size data_dim
        CALL h5screate_simple_f(1, data_dim, dspace_id, error)
        ! Create double/single precision dataset with path '/var' and write data
        if ( cp .eq. dp ) then
            CALL h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
            CALL h5dwrite_f(dset_id, H5T_IEEE_F64LE, vector, data_dim, error)
        else if ( cp .eq. sp ) then
            CALL h5dcreate_f(file_id, dset_name, H5T_IEEE_F32LE, dspace_id, dset_id, error)
            CALL h5dwrite_f(dset_id, H5T_IEEE_F32LE, vector, data_dim, error)
        endif
        ! Close dataset
        CALL h5dclose_f(dset_id, error)
        ! Close dataspace
        CALL h5sclose_f(dspace_id, error)

        ! Close the file
        CALL h5fclose_f(file_id, error)
        ! Close FORTRAN interface
        CALL h5close_f(error)
    end subroutine h5save_R1


    ! subroutine h5save_R2( filename, varname, matrix )
    ! save a real rank 2 matrix to h5 file
    ! Arguments:
    !   filename: [string, Input] h5 filename with path
    !   varname : [string, Input] variable name (real 2d matrix)
    !   matrix  : [double/single 2d matrix, Input] data to be saved
    subroutine h5save_R2( filename, varname, matrix )
        character(len=*), intent(in) :: filename, varname
        real(kind=cp), intent(in), dimension(:,:) :: matrix

        character(len=100) :: dset_name ! dataset name
        integer(HSIZE_T), dimension(2) :: data_dim ! data dimensions

        integer :: error ! error flag
        INTEGER(HID_T) :: file_id  ! file id
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
        ! Create double/single precision dataset with path '/var' and write data
        if ( cp .eq. dp ) then
            CALL h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
            CALL h5dwrite_f(dset_id, H5T_IEEE_F64LE, matrix, data_dim, error)
        else if ( cp .eq. sp ) then
            CALL h5dcreate_f(file_id, dset_name, H5T_IEEE_F32LE, dspace_id, dset_id, error)
            CALL h5dwrite_f(dset_id, H5T_IEEE_F32LE, matrix, data_dim, error)
        endif
        ! Close dataset
        CALL h5dclose_f(dset_id, error)
        ! Close dataspace
        CALL h5sclose_f(dspace_id, error)

        ! Close the file
        CALL h5fclose_f(file_id, error)
        ! Close FORTRAN interface
        CALL h5close_f(error)
    end subroutine h5save_R2


    ! subroutine h5save_C3Partial_Init( filename, varname, full_data_dim )
    ! initialize a complex rank 3 matrix to h5 file in preparation for partial data saving
    !
    ! Arguments:
    !   filename     : [string, Input] h5 filename with path
    !   varname      : [string, Input] variable name
    !   full_data_dim: [integer, size 3, Input] the dimension of the full data matrix
    !
    ! Note:
    !   This function needs to be called once and only once before any partial saving is performed
    subroutine h5save_C3Partial_Init( filename, varname, full_data_dim)
        ! Inputs
        character(len=*), intent(in) :: filename, varname
        integer, intent(in), dimension(3) :: full_data_dim

        integer(HSIZE_T), dimension(3) :: full_data_dim2
        character(len=100) :: dset_name ! dataset name

        integer :: error ! error flag
        INTEGER(HID_T) :: file_id   ! file id
        INTEGER(HID_T) :: group_id  ! group id
        INTEGER(HID_T) :: dspace_id ! dataspace id
        INTEGER(HID_T) :: dset_id   ! dataset id

        logical :: file_exists


        ! Integer type conversion for the full data dimension
        full_data_dim2 = INT( full_data_dim, HSIZE_T)

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

        ! Create disk dataspace with rank 3 and size full_data_dim
        call h5screate_simple_f(3, full_data_dim2, dspace_id, error)

        ! Create double/single precision dataset with path '/var/var_REAL'
        if ( cp .eq. dp ) then
            call h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
        else if ( cp .eq. sp ) then
            call h5dcreate_f(file_id, dset_name, H5T_IEEE_F32LE, dspace_id, dset_id, error)
        endif
        ! Close dataset
        call h5dclose_f(dset_id, error)

        ! Close disk dataspace
        call h5sclose_f(dspace_id, error)
        ! Flush buffers to disk. Not flushing buffers may lead to corrupted hdf5 files
        call h5fflush_f(file_id, H5F_SCOPE_GLOBAL_F, error)

        ! --------------------------- Imaginary part ---------------------------
        dset_name = varname // "/" // varname // "_IMAG"

        ! Create disk dataspace with rank 3 and size full_data_dim
        call h5screate_simple_f(3, full_data_dim2, dspace_id, error)

        ! Create double/single precision dataset with path '/var/var_REAL' and write data
        if ( cp .eq. dp ) then
            call h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
        else if ( cp .eq. sp ) then
            call h5dcreate_f(file_id, dset_name, H5T_IEEE_F32LE, dspace_id, dset_id, error)
        endif
        ! Close dataset
        call h5dclose_f(dset_id, error)

        ! Close disk dataspace
        call h5sclose_f(dspace_id, error)
        ! Flush buffers to disk. Not flushing buffers may lead to corrupted hdf5 files
        call h5fflush_f(file_id, H5F_SCOPE_GLOBAL_F, error)

        ! ------------------------------ Clean up ------------------------------
        ! Close the file
        CALL h5fclose_f(file_id, error)
        ! Close FORTRAN interface
        CALL h5close_f(error)

    end subroutine h5save_C3Partial_Init


    ! subroutine h5save_C3Partial_SingleDim3( filename, varname, matrix, dim3index )
    ! save a complex rank 2 matrix to h5 file as a single plane in the complex rank 3 matrix
    !
    ! Arguments:
    !   filename : [string, Input] h5 filename with path
    !   varname  : [string, Input] variable name
    !   matrix   : [double/single complex 2d matrix, Input] data to be saved
    !   dim3index: [integer, Input] the index for the third dimension
    !
    ! Note:
    !   The variable in this h5 file must already be initilaized using the function h5save_C3Partial_Init
    subroutine h5save_C3Partial_SingleDim3( filename, varname, matrix, dim3index)
        ! Inputs
        character(len=*), intent(in) :: filename, varname
        complex(kind=cp), intent(in), dimension(:,:) :: matrix
        integer, intent(in) :: dim3index

        real(kind=cp), dimension(:,:), allocatable :: temp

        character(len=100) :: dset_name ! dataset name

        ! data dimensions for full data and slice data
        integer(HSIZE_T), dimension(2) :: matrix_dim
        integer(HSIZE_T), dimension(3) :: slice_data_dim
        integer(HSIZE_T), dimension(3) :: slabOffset ! hyperslab offset

        integer :: error ! error flag
        INTEGER(HID_T) :: file_id   ! file id
        INTEGER(HID_T) :: dspace_id ! dataspace id
        INTEGER(HID_T) :: mspace_id ! memeory space id
        INTEGER(HID_T) :: dset_id   ! dataset id


        ! get matrix dimension
        matrix_dim = shape(matrix)
        slice_data_dim(1) = matrix_dim(1)
        slice_data_dim(2) = matrix_dim(2)
        slice_data_dim(3) = 1
        ! Allocate temp buffer
        allocate(temp( slice_data_dim(1), slice_data_dim(2)))

        ! ------------------------ Setup File and Group ------------------------
        ! Initialize hdf5 interface
        call h5open_f(error)
        ! Open file with read and write access
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
        ! Group is already created by h5save_C3Partial_Init

        ! ----------------------------- Real part -----------------------------
        dset_name = varname // "/" // varname // "_REAL"

        ! open dataset
        call h5dopen_f(file_id, dset_name, dset_id, error)

        ! Get disk dataspace ID from the dataset
        call h5dget_space_f(dset_id, dspace_id, error)
        ! select that hyperlsab of the disk dataspace, which is the xz planes at this dim3index (offset dim3index-1)
        slabOffset = (/ 0, 0, dim3index - 1 /)
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)
        ! Create memory dataspace with rank 3 and size slice_data_dim
        call h5screate_simple_f(3, slice_data_dim, mspace_id, error)
        ! select that hyperlsab of the memory dataspace, (no offset)
        slabOffset = (/ 0, 0, 0 /)
        call h5sselect_hyperslab_f(mspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)

        ! Create double/single precision dataset with path '/var/var_REAL' and write data
        temp = real(matrix, cp)
        if ( cp .eq. dp ) then
            call h5dwrite_f(dset_id, H5T_IEEE_F64LE, temp, slice_data_dim, error, &
                    mem_space_id=mspace_id, file_space_id=dspace_id)
        else if ( cp .eq. sp ) then
            call h5dwrite_f(dset_id, H5T_IEEE_F32LE, temp, slice_data_dim, error, &
                    mem_space_id=mspace_id, file_space_id=dspace_id)
        endif
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

        ! open dataset
        call h5dopen_f(file_id, dset_name, dset_id, error)

        ! Get disk dataspace ID from the dataset
        call h5dget_space_f(dset_id, dspace_id, error)
        ! select that hyperlsab of the disk dataspace, which is the xz planes at this dim3index (offset dim3index-1)
        slabOffset = (/ 0, 0, dim3index - 1 /)
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)
        ! Create memory dataspace with rank 3 and size slice_data_dim
        call h5screate_simple_f(3, slice_data_dim, mspace_id, error)
        ! select that hyperlsab of the memory dataspace, (no offset)
        slabOffset = (/ 0, 0, 0 /)
        call h5sselect_hyperslab_f(mspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)

        ! Create double/single precision dataset with path '/var/var_IMAG' and write data
        temp = aimag(matrix)
        if ( cp .eq. dp ) then
            call h5dwrite_f(dset_id, H5T_IEEE_F64LE, temp, slice_data_dim, error, &
                    mem_space_id=mspace_id, file_space_id=dspace_id)
        else if ( cp .eq. sp ) then
            call h5dwrite_f(dset_id, H5T_IEEE_F32LE, temp, slice_data_dim, error, &
                    mem_space_id=mspace_id, file_space_id=dspace_id)
        endif
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
    end subroutine h5save_C3Partial_SingleDim3


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
            call h5save_C3Partial_Init( filename, varname, full_data_dim)
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


    ! subroutine h5save_C3Partial_SingleDim2( filename, varname, matrix, dim2index )
    ! save a complex rank 2 matrix to h5 file as a single plane of dimension 2
    ! of the full complex rank 3 matrix
    !
    ! Arguments:
    !   filename : [string, Input] h5 filename with path
    !   varname  : [string, Input] variable name
    !   matrix   : [double/single complex 3d matrix, size (dim1size, dim3size) Input] data to be saved
    !   dim2index: [integer, Input] the index for the second dimension
    !
    ! Note:
    !   The variable in this h5 file must already be initilaized using the function h5save_C3Partial_Init
    subroutine h5save_C3Partial_SingleDim2( filename, varname, matrix, dim2index)
        ! Inputs
        character(len=*), intent(in) :: filename, varname
        complex(kind=cp), intent(in), dimension(:,:) :: matrix
        integer, intent(in) :: dim2index

        real(kind=cp), dimension(:,:,:), allocatable :: temp

        character(len=100) :: dset_name ! dataset name

        ! data dimensions for full data and slice data
        integer(HSIZE_T), dimension(2) :: matrix_dim
        integer(HSIZE_T), dimension(3) :: slice_data_dim
        integer(HSIZE_T), dimension(3) :: slabOffset ! hyperslab offset

        integer :: error ! error flag
        INTEGER(HID_T) :: file_id   ! file id
        INTEGER(HID_T) :: dspace_id ! dataspace id
        INTEGER(HID_T) :: mspace_id ! memeory space id
        INTEGER(HID_T) :: dset_id   ! dataset id


        ! get matrix dimension
        matrix_dim = shape(matrix)
        slice_data_dim(1) = matrix_dim(1)
        slice_data_dim(2) = 1
        slice_data_dim(3) = matrix_dim(2)
        ! Allocate temp buffer
        allocate(temp( slice_data_dim(1), slice_data_dim(2), slice_data_dim(3)))

        ! ------------------------ Setup File and Group ------------------------
        ! Initialize hdf5 interface
        call h5open_f(error)
        ! Open file with read and write access
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
        ! Group is already created by h5save_C3Partial_Init

        ! ----------------------------- Real part -----------------------------
        dset_name = varname // "/" // varname // "_REAL"

        ! open dataset
        call h5dopen_f(file_id, dset_name, dset_id, error)

        ! Get disk dataspace ID from the dataset
        call h5dget_space_f(dset_id, dspace_id, error)
        ! select that hyperlsab of the disk dataspace, which is the xz planes at this dim2index (offset dim2index-1)
        slabOffset = (/ 0, dim2index - 1, 0 /)
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)
        ! Create memory dataspace with rank 3 and size slice_data_dim
        call h5screate_simple_f(3, slice_data_dim, mspace_id, error)
        ! select that hyperlsab of the memory dataspace, (no offset)
        slabOffset = (/ 0, 0, 0 /)
        call h5sselect_hyperslab_f(mspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)

        ! Create double/single precision dataset with path '/var/var_REAL' and write data
        temp(:,1,:) = real(matrix, cp)
        if ( cp .eq. dp ) then
            call h5dwrite_f(dset_id, H5T_IEEE_F64LE, temp, slice_data_dim, error, &
                    mem_space_id=mspace_id, file_space_id=dspace_id)
        else if ( cp .eq. sp ) then
            call h5dwrite_f(dset_id, H5T_IEEE_F32LE, temp, slice_data_dim, error, &
                    mem_space_id=mspace_id, file_space_id=dspace_id)
        endif
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

        ! open dataset
        call h5dopen_f(file_id, dset_name, dset_id, error)

        ! Get disk dataspace ID from the dataset
        call h5dget_space_f(dset_id, dspace_id, error)
        ! select that hyperlsab of the disk dataspace, which is the xz planes at this dim2index (offset dim2index-1)
        slabOffset = (/ 0, dim2index - 1, 0 /)
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)
        ! Create memory dataspace with rank 3 and size slice_data_dim
        call h5screate_simple_f(3, slice_data_dim, mspace_id, error)
        ! select that hyperlsab of the memory dataspace, (no offset)
        slabOffset = (/ 0, 0, 0 /)
        call h5sselect_hyperslab_f(mspace_id, H5S_SELECT_SET_F, slabOffset, slice_data_dim, error)

        ! Create double/single precision dataset with path '/var/var_IMAG' and write data
        temp(:,1,:) = aimag(matrix)
        if ( cp .eq. dp ) then
            call h5dwrite_f(dset_id, H5T_IEEE_F64LE, temp, slice_data_dim, error, &
                    mem_space_id=mspace_id, file_space_id=dspace_id)
        else if ( cp .eq. sp ) then
            call h5dwrite_f(dset_id, H5T_IEEE_F32LE, temp, slice_data_dim, error, &
                    mem_space_id=mspace_id, file_space_id=dspace_id)
        endif
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
    end subroutine h5save_C3Partial_SingleDim2


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
            call h5save_C3Partial_Init( filename, varname, full_data_dim)
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


end module h5save

