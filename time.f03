module time
    use types, only: sp, dp
    use mpi
    implicit none
    private

    real(kind=dp) :: elapsed_time
    integer :: progress_counter, total_loops

    public :: timer_continue, timer_stop, progress_init, progress_startloop, progress_endloop

contains
    ! subroutine timer_continue( timer )
    ! Continue a timer
    ! Arguments:
    !   timer: [double, IN/OUT] a timer variable
    subroutine timer_continue( timer )
        real(kind=dp), intent(inout) :: timer
        timer = timer - MPI_WTIME()
    end subroutine timer_continue


    ! subroutine timer_stop( timer )
    ! Stop a timer
    ! Arguments:
    !   timer: [double, IN/OUT] a timer variable
    subroutine timer_stop( timer )
        real(kind=dp), intent(inout) :: timer
        timer = timer + MPI_WTIME()
    end subroutine timer_stop


    ! subroutine sec2hms( sec, h, m, s )
    ! convert seconds to hour minutes and seconds
    ! Arguments
    !   sec: [double, input] seconds
    !   h  : [integer, output] hours
    !   m  : [integer, output] minutes
    !   s  : [integer, output] seconds
    subroutine sec2hms( sec, h, m, s )
        real(kind=dp), intent(in) :: sec
        integer, intent(out) :: h, m, s

        h = floor( sec/3600 )
        m = floor( (sec - h*3600)/60 )
        s = floor(  sec - h*3600 - m*60 )
    end subroutine sec2hms


    ! subroutine progress_init( nloops )
    ! Initialize the progress indicator, call before loop
    ! Argument:
    !   nloops: [integer, input] number of do loops
    subroutine progress_init( nloops )
        integer, intent(in) :: nloops

        progress_counter = 1
        total_loops = nloops
        elapsed_time = 0.0d0
    end subroutine progress_init


    ! subroutine progress_startloop
    ! Call at the very beginning of the do loop, it will start the timer
    subroutine progress_startloop
        call timer_continue( elapsed_time )
    end subroutine progress_startloop


    ! subroutine progress_endloop(myid)
    ! Call at the very end of the do loop, it will print the progress at each loop
    ! Arguments
    !   myid: [integer, input] processor ID
    subroutine progress_endloop(myid)
        use iso_fortran_env, only: output_unit
        integer, intent(in) :: myid

        integer :: hhe, mme, sse, hhr, mmr, ssr
        integer :: month, day, hh, mm, ss
        real(kind=dp) :: remaining_time
        integer :: values(8)

        ! stop timer
        call timer_stop( elapsed_time )

        if (myid .eq. 0) then
            ! compute remaining time
            remaining_time = elapsed_time / progress_counter * total_loops - elapsed_time
            ! convert eslapsed time and remaining time to hms
            call sec2hms(   elapsed_time, hhe, mme, sse )
            call sec2hms( remaining_time, hhr, mmr, ssr )

            ! get current date and time
            call date_and_time(VALUES=values)
            ! compute finishing date and time
            month = values(2)
            call sec2hms( values(5)*3600 + values(6)*60 + values(7) + remaining_time, hh, mm, ss )
            day   = values(3) + hh/24
            hh    = hh - (hh/24)*24

            ! write progress
            write(*,"(A,I4.4,A,I4.4,2(A,I3.3,A,I2.2,A,I2.2),5(A,I2.2))") &
            "Loop ", progress_counter, "/", total_loops, &
            "  Elapsed: "  , hhe, ":", mme, ":", sse, &
            "  Remaining: ", hhr, ":", mmr, ":", ssr, &
            "  ETC: ", month, "/", day, "  ", hh, ":", mm, ":", ss

            ! Increment progress counter, which counts the finished loops
            progress_counter = progress_counter + 1

            flush(output_unit)
        endif

    end subroutine progress_endloop

end module time