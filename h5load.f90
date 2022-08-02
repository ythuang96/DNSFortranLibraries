module h5load
    use types, only: sp, dp
    use mpi
    use hdf5
    implicit none
    private
    include 'parameters'

    public :: h5load_R1, h5loadVelocities


contains
    ! function ouput = h5load_R1(filename, varname )
    ! Read real rank 1 matrix (vector) from h5 file
    ! Arguments:
    !   filename: [string, Input]
    !             h5 filename with path
    !   varname : [string, Input]
    !             variable name in h5 file, must be real numerical vector
    ! Output:
    !   output:   [double vector]
    function h5load_R1(filename, varname ) result(vector)
        character(len=*), intent(in) :: filename, varname
        real(kind=dp), dimension(:,:), allocatable :: matrix
        real(kind=cp), dimension(:), allocatable :: vector

        INTEGER(HID_T) :: file_id        ! File identifier
        INTEGER(HID_T) :: dset_id        ! Dataset identifier
        INTEGER(HID_T) :: space_id       ! Dataspace identifier

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
        CALL h5sclose_f(space_id, error)

        ! Allocate dimensions to dset_data for reading
        ALLOCATE(matrix(dim1,dim2))
        ! Get data
        ! H5T_IEEE_F64LE (double) or H5T_IEEE_F32LE (single) has to be
        ! consistent with the variable type of matrix
        ! this will default to a double precision read, and later converted to cp
        CALL h5dread_f(dset_id, H5T_IEEE_F64LE, matrix, data_dims, error)

        ! Convert to cp and reshape matrix into a vector
        vector = reshape( real(matrix, cp), (/ dim1 /) )

        ! close dataset
        CALL h5dclose_f(dset_id, error)
        ! close file
        CALL h5fclose_f(file_id, error)
        ! close h5 interface
        CALL h5close_f(error)
    end function h5load_R1

    ! subroutine h5loadVelocities( filename, varname, Buffer )
    ! Arguments:
    !   filename: [string, Input]
    !             h5 filename with path
    !   varname:  [string, Input]
    !             "u" or "v" or "w" indicating which velocity component to read
    !   Buffer:   [3d double matrix, size (mxf,mzf,myf), Output]
    !             velocity read from the h5 file
    subroutine h5loadVelocities( filename, varname, Buffer )
        character(len=*), intent(in) :: filename
        character(len=1), intent(in) :: varname
        complex(kind=cp), intent(out), dimension(mxf,mzf,myf) :: Buffer
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


        ! Build data from real and imaginary part with conversion to cp
        Buffer = CMPLX( REAL(rtemp, cp), REAL(itemp, cp) , cp )

        ! close file
        CALL h5fclose_f(file_id, error)
        ! close h5 interface
        CALL h5close_f(error)
    end subroutine h5loadVelocities


end module h5load

