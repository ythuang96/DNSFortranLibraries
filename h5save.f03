module h5save
    use types, only: sp, dp
    use mpi
    use hdf5
    implicit none
    private
#   include "parameters"

    public :: check_filename, h5save_logical, h5save_string
    public :: h5save_R, h5save_R1, h5save_R2
    public :: h5save_C2, h5save_C3_sp, h5save_C3_dp

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


    ! subroutine h5save_C3_dp( filename, varname, matrix )
    ! save a complex rank 3 matrix to h5 file
    ! Arguments:
    !   filename: [string, Input] h5 filename with path
    !   varname : [string, Input] variable name (complex 3d matrix)
    !   matrix  : [double complex 3d matrix, Input] data to be saved
    subroutine h5save_C3_dp( filename, varname, matrix )
        character(len=*), intent(in) :: filename, varname
        complex(kind=dp), intent(in), dimension(:,:,:) :: matrix

        real(kind=dp), dimension(:,:,:), allocatable :: temp

        character(len=100) :: dset_name ! dataset name
        integer(HSIZE_T), dimension(3) :: data_dim ! data dimensions

        integer :: error ! error flag
        INTEGER(HID_T) :: file_id  ! file id
        INTEGER(HID_T) :: group_id ! group id
        INTEGER(HID_T) :: dspace_id ! dataspace id
        INTEGER(HID_T) :: dset_id ! dataset id

        logical :: file_exists


        ! get matrix dimensions
        data_dim = shape(matrix)
        ! Allocate temp buffer
        ALLOCATE(temp(data_dim(1), data_dim(2), data_dim(3)))

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
        ! Create dataspace with rank 3 and size data_dim
        CALL h5screate_simple_f(3, data_dim, dspace_id, error)
        ! Create double/single precision dataset with path '/var/var_REAL' and write data
        temp = real(matrix, dp)
        CALL h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
        CALL h5dwrite_f(dset_id, H5T_IEEE_F64LE, temp, data_dim, error)
        ! Close dataset
        CALL h5dclose_f(dset_id, error)
        ! Close dataspace
        CALL h5sclose_f(dspace_id, error)

        ! --------------------------- Imaginary part ---------------------------
        dset_name = varname // "/" // varname // "_IMAG"
        ! Create dataspace with rank 3 and size data_dim
        CALL h5screate_simple_f(3, data_dim, dspace_id, error)
        ! Create double/single precision dataset with path '/var/var_IMAG' and write data
        temp = aimag(matrix)
        CALL h5dcreate_f(file_id, dset_name, H5T_IEEE_F64LE, dspace_id, dset_id, error)
        CALL h5dwrite_f(dset_id, H5T_IEEE_F64LE, temp, data_dim, error)
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
    end subroutine h5save_C3_dp


    ! subroutine h5save_C3_sp( filename, varname, matrix )
    ! save a complex rank 3 matrix to h5 file
    ! Arguments:
    !   filename: [string, Input] h5 filename with path
    !   varname : [string, Input] variable name (complex 3d matrix)
    !   matrix  : [single complex 3d matrix, Input] data to be saved
    subroutine h5save_C3_sp( filename, varname, matrix )
        character(len=*), intent(in) :: filename, varname
        complex(kind=sp), intent(in), dimension(:,:,:) :: matrix

        real(kind=sp), dimension(:,:,:), allocatable :: temp

        character(len=100) :: dset_name ! dataset name
        integer(HSIZE_T), dimension(3) :: data_dim ! data dimensions

        integer :: error ! error flag
        INTEGER(HID_T) :: file_id  ! file id
        INTEGER(HID_T) :: group_id ! group id
        INTEGER(HID_T) :: dspace_id ! dataspace id
        INTEGER(HID_T) :: dset_id ! dataset id

        logical :: file_exists


        ! get matrix dimensions
        data_dim = shape(matrix)
        ! Allocate temp buffer
        ALLOCATE(temp(data_dim(1), data_dim(2), data_dim(3)))

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
        ! Create dataspace with rank 3 and size data_dim
        CALL h5screate_simple_f(3, data_dim, dspace_id, error)
        ! Create double/single precision dataset with path '/var/var_REAL' and write data
        temp = real(matrix, sp)
        CALL h5dcreate_f(file_id, dset_name, H5T_IEEE_F32LE, dspace_id, dset_id, error)
        CALL h5dwrite_f(dset_id, H5T_IEEE_F32LE, temp, data_dim, error)
        ! Close dataset
        CALL h5dclose_f(dset_id, error)
        ! Close dataspace
        CALL h5sclose_f(dspace_id, error)

        ! --------------------------- Imaginary part ---------------------------
        dset_name = varname // "/" // varname // "_IMAG"
        ! Create dataspace with rank 3 and size data_dim
        CALL h5screate_simple_f(3, data_dim, dspace_id, error)
        ! Create double/single precision dataset with path '/var/var_IMAG' and write data
        temp = aimag(matrix)
        CALL h5dcreate_f(file_id, dset_name, H5T_IEEE_F32LE, dspace_id, dset_id, error)
        CALL h5dwrite_f(dset_id, H5T_IEEE_F32LE, temp, data_dim, error)
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
    end subroutine h5save_C3_sp


end module h5save

