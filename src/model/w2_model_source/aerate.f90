subroutine AERATE()
    use GLOBAL;     use MAIN;     use KINETIC;     use TRANS;     use SCREENC
    implicit none
    real, allocatable, dimension(:) :: DZMULTA, SMASS, ATIMON, ATIMOFF, CUMDOMASS
    real, allocatable, dimension(:,:) :: DZMULT
    integer, allocatable, dimension(:) :: IASEG, KTOPA, KBOTA, IPRB, KPRB
    real, allocatable, dimension(:) :: DOOFF, DOON
    real, allocatable, dimension(:) :: ACTUAL_MASS
    logical, allocatable, dimension(:) :: AERATEO2
    integer :: NAER, NLAYERS, KTOP, KBOT
    logical :: CSVFORMAT
    character(len=16) :: CONAER
    character(len=30) :: HEADER
    save

    CSVFORMAT = .false.
    allocate(DZMULT(KMX, IMX))
    DZMULT = 1.0
    open(AERATEFN, FILE="W2_AERATE.NPT", STATUS="OLD")
    read(AERATEFN, "(//A)") HEADER
    do J = 1, 30
        if (HEADER(J:J) == ",") then
            CSVFORMAT = .true.
            exit
        end if
    end do
    rewind(AERATEFN)

    if (CSVFORMAT) then
        read(AERATEFN, *)
        read(AERATEFN, *)
        read(AERATEFN, *) NAER, CONAER
    else
        read(AERATEFN, "(//I8,A16)") NAER, CONAER
    end if

    if (NAER == 0) then
        NAER = 1
    end if
    allocate(DZMULTA(NAER), IASEG(NAER), KTOPA(NAER), KBOTA(NAER), SMASS(NAER), ATIMON(NAER), ATIMOFF(NAER), DOOFF(NAER), DOON(NAER), IPRB(NAER), KPRB(NAER), ACTUAL_MASS(NAER), AERATEO2(NAER), CUMDOMASS(NAER))
    DZMULTA = 0.0;     ACTUAL_MASS = 0.0;     CUMDOMASS = 0.0
    read(AERATEFN, 1013)
    1013 format(/)
    do I = 1, NAER
        if (CSVFORMAT) then
            read(AERATEFN, *) IASEG(I), KTOPA(I), KBOTA(I), SMASS(I), ATIMON(I), ATIMOFF(I), DZMULTA(I), DOOFF(I), DOON(I), IPRB(I), KPRB(I)
        else
            read(AERATEFN, "(I8,I8,I8,F8.0,F8.0,F8.0,3F8.0,2I8)") IASEG(I), KTOPA(I), KBOTA(I), SMASS(I), ATIMON(I), ATIMOFF(I), DZMULTA(I), DOOFF(I), DOON(I), IPRB(I), KPRB(I)
        end if
    end do

    close(AERATEFN)

    if (RESTART_IN) then
        open(AERATEFN, FILE=CONAER, POSITION="APPEND")
        JDAY1 = 0.0
        rewind(AERATEFN)
        read(AERATEFN, "(//)")
        do while (JDAY1 < JDAY)
            read(AERATEFN, "(F9.0)", END=106) JDAY1
        end do
        backspace(AERATEFN)
        106 JDAY1 = 0.0
    else
        open(AERATEFN, FILE=CONAER, STATUS="UNKNOWN")
        write(AERATEFN, "(A,I3,A)") "OUTPUT FILE FOR AERATION INPUT WITH", NAER, " INPUT(S)."
        write(AERATEFN, "(A114)") "JDAY,INSTMASSRATE#1(KGO2/D),CUMMASS#1(KGO2),DOPROBE#1(MG/L),INSTMASSRATE#2(KGO2/D),CUMMASS#2(KGO2),DOPROBE#2(MG/L)"
    end if
    DZMULT = 1.0 ! ALWAYS RESET MIXING COEFFICIENT FOR AERATION SW IPC 2/01/01
    return

    entry DZAERATE() ! FROM W2_ MAIN CODE
    do I = 1, NAER
        DZ(KTOPA(I):KBOTA(I), IASEG(I)) = DZ(KTOPA(I):KBOTA(I), IASEG(I))*DZMULT(KTOPA(I):KBOTA(I), IASEG(I))
    end do
    return

    entry AERATEMASS() ! FROM WQ_CONSTITUENTS

! SECTION FOR HYPOLIMNETIC AERATION 
    DZMULT = 1.0 ! ALWAYS RESET TO 1.0 IN CASE NO AERATION
    do II = 1, NAER
        if (JDAY >= ATIMON(II) .and. JDAY <= ATIMOFF(II)) then

! FIND BRANCH AND WATERBODY FOR ISEG
            do JJB = 1, NBR
                if (BR_INACTIVE(JJB)) then
                    cycle
                end if ! SW 6/12/2017
                if (IASEG(II) >= US(JJB) .and. IASEG(II) <= DS(JJB)) then
                    exit
                end if
            end do

            if (JJB == JB) then ! IF THIS BRANCH DOESN'T HAVE AERATION SKIP IT

                KTOP = MAX(KTWB(JW), KTOPA(II))
                KBOT = MIN(KB(IASEG(II)), KBOTA(II))

                NLAYERS = KBOT - KTOP + 1

                AERATEO2(II) = .true.

                if (O2(KPRB(II), IPRB(II)) > DOON(II) .and. O2(KPRB(II), IPRB(II)) > DOOFF(II)) then
                    AERATEO2(II) = .false.
                    ACTUAL_MASS(II) = 0.0
                end if

                if (AERATEO2(II)) then
                    ACTUAL_MASS(II) = SMASS(II) ! SW 12/28/01
                    CUMDOMASS(II) = CUMDOMASS(II) + SMASS(II)*DLT/86400.
                    do K = KTOP, KBOT

                        DZMULT(K, IASEG(II)) = DZMULTA(II) ! THIS MEANS THE INCREASE IN DZ IS LAGGED ONE TIME STEP
                        CSSB(K, IASEG(II), NDO) = CSSB(K, IASEG(II), NDO) + SMASS(II)/(86.4*REAL(NLAYERS))
! UNITS OF SMASS ARE IN KG/DAY
! TYPICAL UNITS OF CSSB: (MG/L)*(M3/S) TO OONVERT MULTIPLY BY (1KG/10^6MG)*(1000L/1M3)*(86400S/DAY)==86.4

                    end do
                end if
            end if
        else
            AERATEO2(II) = .false.
            ACTUAL_MASS(II) = 0.0
        end if

    end do

    return

    entry AERATEOUTPUT() ! FROM W2_MAIN CODE

    write(AERATEFN, "(F9.3,<NAER>(1X,E12.3,1X,E12.3,1X,F8.3))") JDAY, (ACTUAL_MASS(II), CUMDOMASS(II), O2(KPRB(II), IPRB(II)), II = 1, NAER)

    return
    entry DEALLOCATE_AERATE()
    deallocate(DZMULTA, IASEG, KTOPA, KBOTA, SMASS, ATIMON, ATIMOFF, DOOFF, DOON, IPRB, KPRB, ACTUAL_MASS, AERATEO2, CUMDOMASS)
    deallocate(DZMULT)
    close(AERATEFN)
    return
end subroutine AERATE
