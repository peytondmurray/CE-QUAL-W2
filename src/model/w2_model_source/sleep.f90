! Call the C sleep function. Modified from
! https://cyber.dabamos.de/programming/modernfortran/sleep.html
module posix
    interface
        ! int sleep(seconds_t seconds)
        subroutine c_sleep(seconds) bind(c, name='sleep')
            use, intrinsic :: iso_c_binding, only: c_int
            implicit none
            integer(kind=c_int), value :: seconds
        end subroutine c_sleep
    end interface
end module posix
