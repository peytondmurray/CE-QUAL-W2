subroutine GREGORIAN_DATE()
    use GDAYC
    implicit none
    integer :: INCR

! Determine if new year (regular or leap) and increment year

    do while (JDAYG >= 366)
        if (.not. LEAP_YEAR .and. JDAYG >= 366) then
            JDAYG = JDAYG - 365
            YEAR = YEAR + 1
            LEAP_YEAR = MOD(YEAR, 4) == 0
        else
            if (JDAYG >= 367) then
                JDAYG = JDAYG - 366
                YEAR = YEAR + 1
                LEAP_YEAR = MOD(YEAR, 4) == 0
            else
                exit
            end if
        end if
    end do
    INCR = 0
    if (LEAP_YEAR) then
        INCR = 1
    end if

! Determine month and day of year

    if (JDAYG >= 1 .and. JDAYG < 32) then
        GDAY = JDAYG
        DAYM = 31.0
        MONTH = "  January"
        IMON = 1
    else
        if (JDAYG >= 32 .and. JDAYG < 60 + INCR) then
            GDAY = JDAYG - 31
            DAYM = 29.0
            MONTH = " February"
            IMON = 2
        else
            if (JDAYG >= 60 .and. JDAYG < 91 + INCR) then
                GDAY = JDAYG - 59 - INCR
                DAYM = 31.0
                MONTH = "    March"
                IMON = 3
            else
                if (JDAYG >= 91 .and. JDAYG < 121 + INCR) then
                    GDAY = JDAYG - 90 - INCR
                    DAYM = 30.0
                    MONTH = "    April"
                    IMON = 4
                else
                    if (JDAYG >= 121 .and. JDAYG < 152 + INCR) then
                        GDAY = JDAYG - 120 - INCR
                        DAYM = 31.0
                        MONTH = "      May"
                        IMON = 5
                    else
                        if (JDAYG >= 152 .and. JDAYG < 182 + INCR) then
                            GDAY = JDAYG - 151 - INCR
                            DAYM = 30.0
                            MONTH = "     June"
                            IMON = 6
                        else
                            if (JDAYG >= 182 .and. JDAYG < 213 + INCR) then
                                GDAY = JDAYG - 181 - INCR
                                DAYM = 31.0
                                MONTH = "     July"
                                IMON = 7
                            else
                                if (JDAYG >= 213 .and. JDAYG < 244 + INCR) then
                                    GDAY = JDAYG - 212 - INCR
                                    DAYM = 31.0
                                    MONTH = "   August"
                                    IMON = 8
                                else
                                    if (JDAYG >= 244 .and. JDAYG < 274 + INCR) then
                                        GDAY = JDAYG - 243 - INCR
                                        DAYM = 30.0
                                        MONTH = "September"
                                        IMON = 9
                                    else
                                        if (JDAYG >= 274 .and. JDAYG < 305 + INCR) then
                                            GDAY = JDAYG - 273 - INCR
                                            DAYM = 31.0
                                            MONTH = "  October"
                                            IMON = 10
                                        else
                                            if (JDAYG >= 305 .and. JDAYG < 335 + INCR) then
                                                GDAY = JDAYG - 304 - INCR
                                                DAYM = 30.0
                                                MONTH = " November"
                                                IMON = 11
                                            else
                                                if (JDAYG >= 335 .and. JDAYG < 366 + INCR) then
                                                    GDAY = JDAYG - 334 - INCR
                                                    DAYM = 31.0
                                                    MONTH = " December"
                                                    IMON = 12
                                                end if
                                            end if
                                        end if
                                    end if
                                end if
                            end if
                        end if
                    end if
                end if
            end if
        end if
    end if
end subroutine GREGORIAN_DATE
