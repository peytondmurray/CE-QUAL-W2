subroutine OUTPUTINIT()

    use MAIN
    use GLOBAL;     use NAMESC;     use GEOMC;     use LOGICC;     use PREC;     use SURFHE;     use KINETIC;     use SHADEC;     use EDDY
    use STRUCTURES;     use TRANS;     use TVDC;     use SELWC;     use GDAYC;     use SCREENC;     use TDGAS;     use RSTART
    use MACROPHYTEC;     use POROSITYC;     use ZOOPLANKTONC;     use BIOENERGETICS;     use CEMAVars, only: SEDIMENT_DIAGENESIS
    use ALGAE_TOXINS
    implicit none
    external :: RESTART_OUTPUT

    real :: DIST
    real(4) :: ELC, DLX_OLD
    integer(4), allocatable, dimension(:,:) :: ICOMP

    integer :: JN, IW, JO, IC, JJ
    integer(4) :: IFLAG
    character(len=60) :: TITLEWITH2
    character(len=100) :: TITLEWITH
    character(len=3) :: ICHAR3
    integer, allocatable, dimension(:) :: IBR

    allocate(ICOMP(KMX, IMX))

!***********************************************************************************************************************************
!*                                                           Task 1.5: Outputs                                                    **
!***********************************************************************************************************************************

    allocate(IBR(NBR))
    inquire(FILE="w2_tecplotbr.csv", EXIST=BR_NOTECPLOT(1))
    if (BR_NOTECPLOT(1)) then
        open(CON, FILE="w2_tecplotbr.csv", STATUS="OLD")
        read(CON, *)
        read(CON, *) IC ! NUMBER OF BRANCHES TO PLOT, MUST BE LESS THAN NBR
        read(CON, *)
        read(CON, *) (IBR(J), J = 1, IC)
        do J = 1, IC
            BR_NOTECPLOT(IBR(J)) = .false.
        end do
        close(CON)
    else
        BR_NOTECPLOT = .false.
    end if
    deallocate(IBR)


! Open output files

    if (RESTART_IN) then

        do JW = 1, NWB ! INITIALIZE OUTPUT FOR TECPLOT FOR RESTART  9/4/2019
            if (JW == 1) then
                DIST = 0.0
            end if
            do JB = BS(JW), BE(JW)
                X1(US(JB)) = DIST + DLX(US(JB))/2.
                do I = US(JB) + 1, DS(JB)
                    DIST = DIST + (DLX(I) + DLX(I - 1))/2.0
                    X1(I) = DIST
                end do
                DIST = DIST + DLX(DS(JB))
                X1(DS(JB) + 1) = DIST
            end do
        end do

        do JW = 1, NWB
            if (SNAPSHOT(JW)) then
                open(SNP(JW), FILE=SNPFN(JW), POSITION="APPEND")
            end if
!IF (VECTOR(JW))      OPEN (VPL(JW),FILE=VPLFN(JW),POSITION='APPEND')   *** DSI
            if (SPREADSHEET(JW)) then
                open(SPR(JW), FILE=SPRFN(JW), POSITION="APPEND")
            end if
            if (SPRC(JW) == "     ONV") then
                open(SPRV(JW), FILE=SPRVFN(JW), POSITION="APPEND")
            end if ! SW 9/28/2018
            if (CONTOUR(JW)) then
                open(CPL(JW), FILE=CPLFN(JW), POSITION="APPEND")
            end if
            if (PROFILE(JW)) then
                open(PRF(JW), FILE=PRFFN(JW), POSITION="APPEND")
            end if
            if (FLUX(JW)) then
                open(FLX(JW), FILE=FLXFN(JW), POSITION="APPEND")
                JDAY1 = 0.0
                rewind(FLX(JW))
                do while (JDAY1 < JDAY)
                    read(FLX(JW), "(A72)", END=12) LINE
                    if (LINE(1:8) == "New date") then
                        backspace(FLX(JW))
                        read(FLX(JW), "(8X,F10.0)", END=12) JDAY1
                    end if
                end do
                do J = 1, 14
                    backspace(FLX(JW))
                end do

                12 continue
                write(SEGNUM, "(I0)") JW
                SEGNUM = ADJUSTL(SEGNUM)
                L = LEN_TRIM(SEGNUM)
                open(FLX2(JW), FILE="kflux_wb" // SEGNUM(1:L) // ".csv", POSITION="APPEND")
                JDAY1 = 0.0
                rewind(FLX2(JW))
                read(FLX2(JW), "(/)", END=13)
                do while (JDAY1 < JDAY)
                    read(FLX2(JW), "(F10.0)", END=13) JDAY1
                end do
                backspace(FLX2(JW))
                13 JDAY1 = 0.0
            end if
            if (SNAPSHOT(JW)) then
                rewind(SNP(JW))
                do while (.true.)
                    read(SNP(JW), "(A72)", END=100) LINE
                    if (LINE(26:28) == "NIT") then
                        backspace(SNP(JW))
                        read(SNP(JW), "(31X,I10)", END=100) NIT1
                        if (NIT1 > NIT) then
                            do J = 1, 24
                                backspace(SNP(JW))
                            end do
                            exit
                        end if
                    end if
                end do
            end if
            100 continue
            if (SPREADSHEET(JW)) then
                rewind(SPR(JW))
                read(SPR(JW), *)
                do while (JDAY1 < JDAY)
                    read(SPR(JW), *, END=101) LINE(1:38), JDAY1 !'(A,F10.0)'
                end do
                backspace(SPR(JW))
                101 continue
                JDAY1 = 0.0
                if (SPRC(JW) == "     ONV") then
                    rewind(SPRV(JW))
                    read(SPRV(JW), *)
                    do while (JDAY1 < JDAY)
                        read(SPRV(JW), *, END=104) LINE(1:38), JDAY1 !'(A,F10.0)'
                    end do
                    backspace(SPRV(JW))
                    104 continue
                    JDAY1 = 0.0
                end if
            end if
            if (PROFILE(JW) .and. iprf(1, 1) /= -1) then ! SW 4/1/2016
                rewind(PRF(JW))
                read(PRF(JW), "(A)") (LINE, J = 1, 11)
                read(PRF(JW), "(8I8)") I
                read(PRF(JW), "(10I8)") (I, J = 1, NIPRF(JW))
                read(PRF(JW), "(20(1X,A))") LINE(1:8), (LINE(1:3), JC = 1, NCT), (LINE(1:3), JD = 1, NDC)
                read(PRF(JW), "(2A)") LINE(1:26), (LINE(1:26), JC = 1, NCT), (LINE(1:43), JD = 1, NDC)
                do while (JDAY1 < JDAY)
                    read(PRF(JW), "(A72)", END=102) LINE
                    L1 = 0
                    L1 = SCAN(LINE, ",")
                    if (L1 /= 0) then
                        backspace(PRF(JW))
                        read(PRF(JW), "(F8.0)", END=102) JDAY1
                    end if
                end do
                backspace(PRF(JW))
                JDAY1 = 0.0
            end if
            102 continue
            JDAY1 = 0.0
            if (CONTOUR(JW)) then
                if (TECPLOT(JW) /= "      ON") then
                    rewind(CPL(JW))
                    do while (JDAY1 < JDAY)
                        read(CPL(JW), "(A72)", END=103) LINE
                        if (LINE(1:8) == "New date") then
                            backspace(CPL(JW))
                            read(CPL(JW), "(A,F12.4)", END=103) LINE(1:9), JDAY1
                        end if
                    end do
                    backspace(CPL(JW))
                    JDAY1 = 0.0
                else
                    rewind(CPL(JW))
                    do while (JDAY1 < JDAY)
                        read(CPL(JW), "(A72)", END=103) LINE
                        if (LINE(1:8) == 'ZONE T="') then
                            backspace(CPL(JW))
                            read(CPL(JW), "(A,F9.0)", END=103) LINE(1:8), JDAY1
                        end if
                    end do
                    backspace(CPL(JW))
                    JDAY1 = 0.0
                end if
            end if
            103 continue
        end do

        if (VECTOR(1)) then
            open(VPL(1), FILE=VPLFN(1), STATUS="UNKNOWN", ACCESS="SEQUENTIAL", FORM="BINARY", POSITION="APPEND")
        end if

        if (DOWNSTREAM_OUTFLOW) then
            JFILE = 0
            L1 = SCAN(WDOFN, ".", BACK=.true.) ! SW 8/22/14 CHECKS FROM RIGHTHAND SIDE NOT LEFTHANDSIDE
            do JWD = 1, NIWDO
                write(SEGNUM, "(I0)") IWDO(JWD)
                SEGNUM = ADJUSTL(SEGNUM)
                L = LEN_TRIM(SEGNUM)
                open(WDO(JWD, 1), FILE="qwo_" // SEGNUM(1:L) // WDOFN(L1:L1 + 4), POSITION="APPEND") ! '.opt' SW 4/14/2017
                rewind(WDO(JWD, 1))
                read(WDO(JWD, 1), "(//)", END=106)
                do while (JDAY1 < JDAY)
                    read(WDO(JWD, 1), *, END=106) JDAY1 !'(F8.0)'
                end do
                backspace(WDO(JWD, 1))
                106 JDAY1 = 0.0
                open(WDO(JWD, 2), FILE="two_" // SEGNUM(1:L) // WDOFN(L1:L1 + 4), POSITION="APPEND") ! '.opt' SW 4/14/2017 repeated for all .opt below for Downstream outflow
                rewind(WDO(JWD, 2))
                read(WDO(JWD, 2), "(//)", END=107)
                do while (JDAY1 < JDAY)
                    read(WDO(JWD, 2), *, END=107) JDAY1 !'(F8.0)'
                end do
                backspace(WDO(JWD, 2))
                107 JDAY1 = 0.0
                if (CONSTITUENTS) then
                    open(WDO(JWD, 3), FILE="cwo_" // SEGNUM(1:L) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                    rewind(WDO(JWD, 3))
                    read(WDO(JWD, 3), "(//)", END=108)
                    do while (JDAY1 < JDAY)
                        read(WDO(JWD, 3), *, END=108) JDAY1 ! '(F8.0)'
                    end do
                    backspace(WDO(JWD, 3))
                    108 continue
                    JDAY1 = 0.0
                end if
                if (DERIVED_CALC) then
                    open(WDO(JWD, 4), FILE="dwo_" // SEGNUM(1:L) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                    rewind(WDO(JWD, 4))
                    read(WDO(JWD, 4), "(//)")
                    do while (JDAY1 < JDAY)
                        read(WDO(JWD, 4), *, END=109) JDAY1 ! '(F8.0)'
                    end do
                    backspace(WDO(JWD, 4))
                    109 continue
                    JDAY1 = 0.0
                end if

! Determine the # of withdrawals at the WITH SEG
                do JB = 1, NBR ! structures
                    if (IWDO(JWD) == DS(JB) .and. NSTR(JB) /= 0) then
                        do JS = 1, NSTR(JB)
                            JFILE = JFILE + 1
                            write(SEGNUM2, "(I0)") IWDO(JWD)
                            SEGNUM2 = ADJUSTL(SEGNUM2)
                            L2 = LEN_TRIM(SEGNUM2)
                            TITLEWITH2 = " STR WITHDRAWAL AT SEG" // SEGNUM2(1:L2)
                            write(SEGNUM, "(F10.0)") ESTR(JS, JB)
                            SEGNUM = ADJUSTL(SEGNUM)
                            L = LEN_TRIM(SEGNUM)
                            TITLEWITH = ADJUSTL(TITLEWITH2) // " ELEV OF WITHDRAWAL CENTERLINE:" // SEGNUM(1:L)
                            write(SEGNUM, "(I0)") JS
                            SEGNUM = ADJUSTL(SEGNUM)
                            L = LEN_TRIM(SEGNUM)
                            WDO2(JFILE, 1) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 1), FILE="qwo_str" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                            rewind(WDO2(JFILE, 1))
                            read(WDO2(JFILE, 1), "(//)")
                            do while (JDAY1 < JDAY)
                                read(WDO2(JFILE, 1), *, END=110) JDAY1 ! '(F8.0)'
                            end do
                            backspace(WDO2(JFILE, 1))
                            110 continue
                            JDAY1 = 0.0

                            WDO2(JFILE, 2) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 2), FILE="two_str" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                            rewind(WDO2(JFILE, 2))
                            read(WDO2(JFILE, 2), "(//)")
                            do while (JDAY1 < JDAY)
                                read(WDO2(JFILE, 2), *, END=111) JDAY1 ! '(F8.0)'
                            end do
                            backspace(WDO2(JFILE, 2))
                            111 continue
                            JDAY1 = 0.0
                            if (CONSTITUENTS) then
                                WDO2(JFILE, 3) = NUNIT;                                 NUNIT = NUNIT + 1
                                open(WDO2(JFILE, 3), FILE="cwo_str" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                                rewind(WDO2(JFILE, 3))
                                read(WDO2(JFILE, 3), "(//)")
                                do while (JDAY1 < JDAY)
                                    read(WDO2(JFILE, 3), *, END=112) JDAY1 ! '(F8.0)'
                                end do
                                backspace(WDO2(JFILE, 3))
                                112 continue
                                JDAY1 = 0.0
                            end if

                            if (DERIVED_CALC) then
                                WDO2(JFILE, 4) = NUNIT;                                 NUNIT = NUNIT + 1
                                open(WDO2(JFILE, 4), FILE="dwo_str" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                                rewind(WDO2(JFILE, 4))
                                read(WDO2(JFILE, 4), "(//)")
                                do while (JDAY1 < JDAY)
                                    read(WDO2(JFILE, 4), *, END=113) JDAY1 !'(F8.0)'
                                end do
                                backspace(WDO2(JFILE, 4))
                                113 continue
                                JDAY1 = 0.0
                            end if
                        end do
                    end if
                end do

                do JS = 1, NWD ! withdrawals
                    if (IWDO(JWD) == IWD(JS)) then
                        JFILE = JFILE + 1
                        write(SEGNUM2, "(I0)") IWDO(JWD)
                        SEGNUM2 = ADJUSTL(SEGNUM2)
                        L2 = LEN_TRIM(SEGNUM2)
                        TITLEWITH2 = " WITHDRAWAL AT SEG" // SEGNUM2(1:L2)
                        write(SEGNUM, "(F10.0)") EWD(JS)
                        SEGNUM = ADJUSTL(SEGNUM)
                        L = LEN_TRIM(SEGNUM)
                        TITLEWITH = ADJUSTL(TITLEWITH2) // " ELEV OF WITHDRAWAL CENTERLINE:" // SEGNUM(1:L)
                        write(SEGNUM, "(I0)") JS
                        SEGNUM = ADJUSTL(SEGNUM)
                        L = LEN_TRIM(SEGNUM)
                        WDO2(JFILE, 1) = NUNIT;                         NUNIT = NUNIT + 1
                        open(WDO2(JFILE, 1), FILE="qwo_wd" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                        rewind(WDO2(JFILE, 1))
                        read(WDO2(JFILE, 1), "(//)")
                        do while (JDAY1 < JDAY)
                            read(WDO2(JFILE, 1), *, END=114) JDAY1 !'(F8.0)'
                        end do
                        backspace(WDO2(JFILE, 1))
                        114 continue
                        JDAY1 = 0.0
                        WDO2(JFILE, 2) = NUNIT;                         NUNIT = NUNIT + 1
                        open(WDO2(JFILE, 2), FILE="two_wd" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                        rewind(WDO2(JFILE, 2))
                        read(WDO2(JFILE, 2), "(//)")
                        do while (JDAY1 < JDAY)
                            read(WDO2(JFILE, 2), *, END=115) JDAY1 !'(F8.0)'
                        end do
                        backspace(WDO2(JFILE, 2))
                        115 continue
                        JDAY1 = 0.0
                        if (CONSTITUENTS) then
                            WDO2(JFILE, 3) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 3), FILE="cwo_wd" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                            rewind(WDO2(JFILE, 3))
                            read(WDO2(JFILE, 3), "(//)")
                            do while (JDAY1 < JDAY)
                                read(WDO2(JFILE, 3), *, END=116) JDAY1 !'(F8.0)'
                            end do
                            backspace(WDO2(JFILE, 3))
                            116 continue
                            JDAY1 = 0.0
                        end if
                        if (DERIVED_CALC) then
                            WDO2(JFILE, 4) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 4), FILE="dwo_wd" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                            rewind(WDO2(JFILE, 4))
                            read(WDO2(JFILE, 4), "(//)")
                            do while (JDAY1 < JDAY)
                                read(WDO2(JFILE, 4), *, END=117) JDAY1 !'(F8.0)'
                            end do
                            backspace(WDO2(JFILE, 4))
                            117 continue
                            JDAY1 = 0.0
                        end if
                    end if
                end do

                do JS = 1, NSP ! spillways
                    if (IWDO(JWD) == IUSP(JS)) then
                        JFILE = JFILE + 1
                        write(SEGNUM2, "(I0)") IWDO(JWD)
                        SEGNUM2 = ADJUSTL(SEGNUM2)
                        L2 = LEN_TRIM(SEGNUM2)
                        TITLEWITH2 = " SPILLWAY WITHDRAWAL AT SEG" // SEGNUM2(1:L2)
                        write(SEGNUM, "(F10.0)") ESP(JS)
                        SEGNUM = ADJUSTL(SEGNUM)
                        L = LEN_TRIM(SEGNUM)
                        TITLEWITH = ADJUSTL(TITLEWITH2) // " ELEV OF WITHDRAWAL CENTERLINE:" // SEGNUM(1:L)
                        write(SEGNUM, "(I0)") JS
                        SEGNUM = ADJUSTL(SEGNUM)
                        L = LEN_TRIM(SEGNUM)
                        WDO2(JFILE, 1) = NUNIT;                         NUNIT = NUNIT + 1
                        open(WDO2(JFILE, 1), FILE="qwo_sp" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                        rewind(WDO2(JFILE, 1))
                        read(WDO2(JFILE, 1), "(//)")
                        do while (JDAY1 < JDAY)
                            read(WDO2(JFILE, 1), *, END=118) JDAY1 !'(F8.0)'
                        end do
                        backspace(WDO2(JFILE, 1))
                        118 continue
                        JDAY1 = 0.0
                        WDO2(JFILE, 2) = NUNIT;                         NUNIT = NUNIT + 1
                        open(WDO2(JFILE, 2), FILE="two_sp" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                        rewind(WDO2(JFILE, 2))
                        read(WDO2(JFILE, 2), "(//)")
                        do while (JDAY1 < JDAY)
                            read(WDO2(JFILE, 2), *, END=119) JDAY1 !'(F8.0)'
                        end do
                        backspace(WDO2(JFILE, 2))
                        119 continue
                        JDAY1 = 0.0
                        if (CONSTITUENTS) then
                            WDO2(JFILE, 3) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 3), FILE="cwo_sp" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                            rewind(WDO2(JFILE, 3))
                            read(WDO2(JFILE, 3), "(//)")
                            do while (JDAY1 < JDAY)
                                read(WDO2(JFILE, 3), *, END=120) JDAY1 !'(F8.0)'
                            end do
                            backspace(WDO2(JFILE, 3))
                            120 continue
                            JDAY1 = 0.0
                        end if
                        if (DERIVED_CALC) then
                            WDO2(JFILE, 4) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 4), FILE="dwo_sp" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                            rewind(WDO2(JFILE, 4))
                            read(WDO2(JFILE, 4), "(//)")
                            do while (JDAY1 < JDAY)
                                read(WDO2(JFILE, 4), *, END=121) JDAY1 !'(F8.0)'
                            end do
                            backspace(WDO2(JFILE, 4))
                            121 continue
                            JDAY1 = 0.0
                        end if
                    end if
                end do


                do JS = 1, NPU ! pumps
                    if (IWDO(JWD) == IUPU(JS)) then
                        JFILE = JFILE + 1
                        write(SEGNUM2, "(I0)") IWDO(JWD)
                        SEGNUM2 = ADJUSTL(SEGNUM2)
                        L2 = LEN_TRIM(SEGNUM2)
                        TITLEWITH2 = "PUMP WITHDRAWAL AT SEG" // SEGNUM2(1:L2)
                        write(SEGNUM, "(F10.0)") EPU(JS)
                        SEGNUM = ADJUSTL(SEGNUM)
                        L = LEN_TRIM(SEGNUM)
                        TITLEWITH = ADJUSTL(TITLEWITH2) // " ELEV OF WITHDRAWAL CENTERLINE:" // SEGNUM(1:L)
                        write(SEGNUM, "(I0)") JS
                        SEGNUM = ADJUSTL(SEGNUM)
                        L = LEN_TRIM(SEGNUM)
                        WDO2(JFILE, 1) = NUNIT;                         NUNIT = NUNIT + 1
                        open(WDO2(JFILE, 1), FILE="qwo_pmp" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                        rewind(WDO2(JFILE, 1))
                        read(WDO2(JFILE, 1), "(//)")
                        do while (JDAY1 < JDAY)
                            read(WDO2(JFILE, 1), *, END=122) JDAY1 !'(F8.0)'
                        end do
                        backspace(WDO2(JFILE, 1))
                        122 continue
                        JDAY1 = 0.0
                        WDO2(JFILE, 2) = NUNIT;                         NUNIT = NUNIT + 1
                        open(WDO2(JFILE, 2), FILE="two_pmp" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                        rewind(WDO2(JFILE, 2))
                        read(WDO2(JFILE, 2), "(//)")
                        do while (JDAY1 < JDAY)
                            read(WDO2(JFILE, 2), *, END=123) JDAY1 !'(F8.0)'
                        end do
                        backspace(WDO2(JFILE, 2))
                        123 continue
                        JDAY1 = 0.0
                        if (CONSTITUENTS) then
                            WDO2(JFILE, 3) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 3), FILE="cwo_pmp" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                            rewind(WDO2(JFILE, 3))
                            read(WDO2(JFILE, 3), "(//)")
                            do while (JDAY1 < JDAY)
                                read(WDO2(JFILE, 3), *, END=124) JDAY1 !'(F8.0)'
                            end do
                            backspace(WDO2(JFILE, 3))
                            124 continue
                            JDAY1 = 0.0
                        end if
                        if (DERIVED_CALC) then
                            WDO2(JFILE, 4) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 4), FILE="dwo_pmp" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                            rewind(WDO2(JFILE, 4))
                            read(WDO2(JFILE, 4), "(//)")
                            do while (JDAY1 < JDAY)
                                read(WDO2(JFILE, 4), *, END=125) JDAY1
                            end do
                            backspace(WDO2(JFILE, 4))
                            125 continue
                            JDAY1 = 0.0
                        end if
                    end if
                end do


                do JS = 1, NPI ! pipes
                    if (IWDO(JWD) == IUPI(JS)) then
                        JFILE = JFILE + 1
                        write(SEGNUM2, "(I0)") IWDO(JWD)
                        SEGNUM2 = ADJUSTL(SEGNUM2)
                        L2 = LEN_TRIM(SEGNUM2)
                        TITLEWITH2 = "PIPE WITHDRAWAL AT SEG" // SEGNUM2(1:L2)
                        write(SEGNUM, "(F10.0)") EUPI(JS)
                        SEGNUM = ADJUSTL(SEGNUM)
                        L = LEN_TRIM(SEGNUM)
                        TITLEWITH = ADJUSTL(TITLEWITH2) // " ELEV OF WITHDRAWAL CENTERLINE:" // SEGNUM(1:L)
                        write(SEGNUM, "(I0)") JS
                        SEGNUM = ADJUSTL(SEGNUM)
                        L = LEN_TRIM(SEGNUM)
                        WDO2(JFILE, 1) = NUNIT;                         NUNIT = NUNIT + 1
                        open(WDO2(JFILE, 1), FILE="qwo_pipe" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                        rewind(WDO2(JFILE, 1))
                        read(WDO2(JFILE, 1), "(//)")
                        do while (JDAY1 < JDAY)
                            read(WDO2(JFILE, 1), *, END=126) JDAY1
                        end do
                        backspace(WDO2(JFILE, 1))
                        126 continue
                        JDAY1 = 0.0
                        WDO2(JFILE, 2) = NUNIT;                         NUNIT = NUNIT + 1
                        open(WDO2(JFILE, 2), FILE="two_pipe" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                        rewind(WDO2(JFILE, 2))
                        read(WDO2(JFILE, 2), "(//)")
                        do while (JDAY1 < JDAY)
                            read(WDO2(JFILE, 2), *, END=127) JDAY1
                        end do
                        backspace(WDO2(JFILE, 2))
                        127 continue
                        JDAY1 = 0.0
                        if (CONSTITUENTS) then
                            WDO2(JFILE, 3) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 3), FILE="cwo_pipe" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                            rewind(WDO2(JFILE, 3))
                            read(WDO2(JFILE, 3), "(//)")
                            do while (JDAY1 < JDAY)
                                read(WDO2(JFILE, 3), *, END=128) JDAY1
                            end do
                            backspace(WDO2(JFILE, 3))
                            128 continue
                            JDAY1 = 0.0
                        end if
                        if (DERIVED_CALC) then
                            WDO2(JFILE, 4) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 4), FILE="dwo_pipe" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                            rewind(WDO2(JFILE, 4))
                            read(WDO2(JFILE, 4), "(//)")
                            do while (JDAY1 < JDAY)
                                read(WDO2(JFILE, 4), *, END=129) JDAY1
                            end do
                            backspace(WDO2(JFILE, 4))
                            129 continue
                            JDAY1 = 0.0
                        end if
                    end if
                end do

                do JS = 1, NGT ! gates
                    if (IWDO(JWD) == IUGT(JS)) then
                        JFILE = JFILE + 1
                        write(SEGNUM2, "(I0)") IWDO(JWD)
                        SEGNUM2 = ADJUSTL(SEGNUM2)
                        L2 = LEN_TRIM(SEGNUM2)
                        TITLEWITH2 = "GATE WITHDRAWAL AT SEG" // SEGNUM2(1:L2)
                        write(SEGNUM, "(F10.0)") EGT(JS)
                        SEGNUM = ADJUSTL(SEGNUM)
                        L = LEN_TRIM(SEGNUM)
                        TITLEWITH = ADJUSTL(TITLEWITH2) // " ELEV OF WITHDRAWAL CENTERLINE:" // SEGNUM(1:L)
                        write(SEGNUM, "(I0)") JS
                        SEGNUM = ADJUSTL(SEGNUM)
                        L = LEN_TRIM(SEGNUM)
                        WDO2(JFILE, 1) = NUNIT;                         NUNIT = NUNIT + 1
                        open(WDO2(JFILE, 1), FILE="qwo_gate" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                        rewind(WDO2(JFILE, 1))
                        read(WDO2(JFILE, 1), "(//)")
                        do while (JDAY1 < JDAY)
                            read(WDO2(JFILE, 1), *, END=130) JDAY1
                        end do
                        backspace(WDO2(JFILE, 1))
                        130 continue
                        JDAY1 = 0.0
                        WDO2(JFILE, 2) = NUNIT;                         NUNIT = NUNIT + 1
                        open(WDO2(JFILE, 2), FILE="two_gate" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                        rewind(WDO2(JFILE, 2))
                        read(WDO2(JFILE, 2), "(//)")
                        do while (JDAY1 < JDAY)
                            read(WDO2(JFILE, 2), *, END=131) JDAY1
                        end do
                        backspace(WDO2(JFILE, 2))
                        131 continue
                        JDAY1 = 0.0
                        if (CONSTITUENTS) then
                            WDO2(JFILE, 3) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 3), FILE="cwo_gate" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                            rewind(WDO2(JFILE, 3))
                            read(WDO2(JFILE, 3), "(//)")
                            do while (JDAY1 < JDAY)
                                read(WDO2(JFILE, 3), *, END=132) JDAY1
                            end do
                            backspace(WDO2(JFILE, 3))
                            132 continue
                            JDAY1 = 0.0
                        end if
                        if (DERIVED_CALC) then
                            WDO2(JFILE, 4) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 4), FILE="dwo_gate" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), POSITION="APPEND")
                            rewind(WDO2(JFILE, 4))
                            read(WDO2(JFILE, 4), "(//)")
                            do while (JDAY1 < JDAY)
                                read(WDO2(JFILE, 4), *, END=133) JDAY1
                            end do
                            backspace(WDO2(JFILE, 4))
                            133 continue
                            JDAY1 = 0.0
                        end if
                    end if
                end do

            end do
        end if
! BIOENERGETICS mlm
        if (BIOEXP) then
            JDAY1 = 0.0
            L1 = SCAN(BIOFN, ".")
            do J = 1, NIBIO
                write(SEGNUM, "(I0)") J
                SEGNUM = ADJUSTL(SEGNUM)
                L = LEN_TRIM(SEGNUM)
                BIOFN = BIOFN(1:L1 - 1) // "_" // SEGNUM(1:L) // ".OPT"
                open(BIOEXPFN(J), FILE=BIOFN, POSITION="APPEND")
                rewind(BIOEXPFN(J))
                read(bioexpfn(j), "(/)", end=1111)
                do while (JDAY1 < JDAY)
                    read(bioexpfn(j), "(F8.0)", end=1111) JDAY1
                end do
                backspace(BIOEXPFN(J))
                1111 continue
                JDAY1 = 0.0
            end do

            L1 = SCAN(WEIGHTFN, ".")
            do J = 1, NIBIO
                write(SEGNUM, "(I0)") J
                SEGNUM = ADJUSTL(SEGNUM)
                L = LEN_TRIM(SEGNUM)
                WEIGHTFN = WEIGHTFN(1:L1 - 1) // "_" // SEGNUM(1:L) // ".OPT"
                open(WEIGHTNUM(J), FILE=WEIGHTFN, POSITION="APPEND")
                rewind(WEIGHTNUM(J))
                read(weightnum(j), "(/)", end=1112)
                do while (JDAY1 < JDAY)
                    read(weightnum(j), "(F10.0)", end=1112) JDAY1
                end do
                backspace(WEIGHTNUM(J))
                1112 continue
                JDAY1 = 0.0
            end do
        end if

        if (TIME_SERIES) then
            L1 = SCAN(TSRFN1, ".", BACK=.true.) ! SW 8/22/14 CHECKS FROM RIGHTHAND SIDE NOT LEFTHANDSIDE
            do J = 1, NIKTSR
                write(SEGNUM, "(I0)") ITSR(J)
                SEGNUM = ADJUSTL(SEGNUM)
                L = LEN_TRIM(SEGNUM)
                write(SEGNUM2, "(I0)") J
                SEGNUM2 = ADJUSTL(SEGNUM2)
                L2 = LEN_TRIM(SEGNUM2)
                TSRFN = TSRFN1(1:L1 - 1) // "_" // SEGNUM2(1:L2) // "_seg" // SEGNUM(1:L) // "." // TSRFN1(L1 + 1:L1 + 4) !'.opt'
                open(TSR(J), FILE=TSRFN, POSITION="APPEND")
                rewind(TSR(J))
! READ   (TSR(J),'(A72)',END=140)   (LINE,I=1,11)
                read(TSR(J), "(/F10.3)", END=140) JDAYTS
                do while (JDAYTS < JDAY)
                    read(TSR(J), "(F10.0)", END=140) JDAYTS
                end do
                backspace(TSR(J))
                140 continue
            end do
        end if
        do JW = 1, NWB
            if (FLOWBALC == "      ON") then
                open(FLOWBFN, FILE="flowbal.csv", POSITION="APPEND")
                JDAY1 = 0.0
                rewind(FLOWBFN)
                read(FLOWBFN, "(/)", END=141)
                do while (JDAY1 < JDAY)
                    read(FLOWBFN, *, END=141) JDAY1
                end do
                backspace(FLOWBFN)
                141 JDAY1 = 0.0
                exit
            end if
        end do
        do JW = 1, NWB
            if (NPBALC == "      ON") then
                open(MASSBFN, FILE="massbal.csv", POSITION="APPEND")
                JDAY1 = 0.0
                rewind(MASSBFN)
                read(MASSBFN, "(/)", END=143)
                do while (JDAY1 < JDAY)
                    read(MASSBFN, *, END=143) JDAY1
                end do
                backspace(MASSBFN)
                143 JDAY1 = 0.0
                exit
            end if
        end do
        if (WLC == "      ON") then
            open(WLFN, FILE="wl.csv", POSITION="APPEND")
            JDAY1 = 0.0
            rewind(WLFN)
            read(WLFN, "(//)", END=142)
            do while (JDAY1 < JDAY)
                read(WLFN, "(F10.0)", END=142) JDAY1
            end do
            backspace(WLFN)
            142 JDAY1 = 0.0
        end if
    else
! *** Cold Start Output File Initialization

        do JW = 1, NWB ! moved so that x1 get initialized even if tecplot (and contour) OFF for first water body but ON for others   ! cb 10/10/10
!c calculating longitudinal distance of segments
            if (JW == 1) then
                DIST = 0.0
            end if
            do JB = BS(JW), BE(JW)
                X1(US(JB)) = DIST + DLX(US(JB))/2.
                do I = US(JB) + 1, DS(JB)
                    DIST = DIST + (DLX(I) + DLX(I - 1))/2.0
                    X1(I) = DIST
                end do
                DIST = DIST + DLX(DS(JB))
                X1(DS(JB) + 1) = DIST
            end do
        end do

        inquire(FILE="w2_lake_river_contour.csv", EXIST=LAKE_RIVER_CONTOURC)

        if (LAKE_RIVER_CONTOURC) then
            open(LAKE_RIVER_CONTOUR, FILE="w2_lake_river_contour.csv", STATUS="OLD", IOSTAT=I)
            read(LAKE_RIVER_CONTOUR, *)
            read(LAKE_RIVER_CONTOUR, *) LAKE_RIVER_CONTOUR_ON
            read(LAKE_RIVER_CONTOUR, *)
            read(LAKE_RIVER_CONTOUR, *) NUM_LAKE_CONTOUR, LAKE_CONTOUR_FORMAT
            read(LAKE_RIVER_CONTOUR, *)
            do JJ = 1, NUM_LAKE_CONTOUR
                read(LAKE_RIVER_CONTOUR, *) LAKE_CONTOUR_SEG(JJ), LAKE_CONTOUR_START(JJ), LAKE_CONTOUR_FREQ(JJ)
                if (JDAY >= LAKE_CONTOUR_START(JJ)) then
                    NXT_LAKE_CONTOUR(JJ) = JDAY
                else
                    NXT_LAKE_CONTOUR(JJ) = LAKE_CONTOUR_START(JJ)
                end if

                do JW = 1, NWB
                    if (US(BS(JW)) <= LAKE_CONTOUR_SEG(JJ) .and. DS(BE(JW)) >= LAKE_CONTOUR_SEG(JJ)) then
                        JW_LAKE_CONTOUR(JJ) = JW
                        exit
                    end if
                end do
            end do
            read(LAKE_RIVER_CONTOUR, *)
            read(LAKE_RIVER_CONTOUR, *) NUM_RIVER_CONTOUR, RIVER_CONTOUR_FORMAT
            read(LAKE_RIVER_CONTOUR, *)
            do JJ = 1, NUM_RIVER_CONTOUR
                read(LAKE_RIVER_CONTOUR, *) RIVER_CONTOUR_BR1(JJ), RIVER_CONTOUR_BR2(JJ), RIVER_CONTOUR_START(JJ), RIVER_CONTOUR_FREQ(JJ)
                if (JDAY >= RIVER_CONTOUR_START(JJ)) then
                    NXT_RIVER_CONTOUR(JJ) = JDAY
                else
                    NXT_RIVER_CONTOUR(JJ) = RIVER_CONTOUR_START(JJ)
                end if
                do JW = 1, NWB
                    if (BS(JW) <= RIVER_CONTOUR_BR1(JJ) .and. BE(JW) >= RIVER_CONTOUR_BR1(JJ)) then
                        JW_RIVER_CONTOUR(JJ) = JW
                        exit
                    end if
                end do
            end do
            close(LAKE_RIVER_CONTOUR)
            if (LAKE_RIVER_CONTOUR_ON == "ON") then
                do JJ = 1, NUM_LAKE_CONTOUR
                    write(ICHAR3, "(I3)") LAKE_CONTOUR_SEG(JJ)
                    FILE_LAKE_CONTOUR_T(JJ) = "LakeContour_T_Seg" // ADJUSTL(trim(ICHAR3)) // ".csv"
                    open(LAKE_RIVER_CONTOUR + JJ - 1, FILE=FILE_LAKE_CONTOUR_T(JJ), status="unknown")
                    if (LAKE_CONTOUR_FORMAT == 1) then
                        write(LAKE_RIVER_CONTOUR + JJ - 1, *) "JDAY,ELEVATION(m),TEMPERATURE(C)"
                    else
                        write(LAKE_RIVER_CONTOUR + JJ - 1, '("TIME,",*(F8.2,","))') ((EL(K, LAKE_CONTOUR_SEG(JJ)) + EL(K, LAKE_CONTOUR_SEG(JJ)))*0.5, K = 2, KB(LAKE_CONTOUR_SEG(JJ)))
                    end if
                    if (OXYGEN_DEMAND) then
                        write(ICHAR3, "(I3)") LAKE_CONTOUR_SEG(JJ)
                        FILE_LAKE_CONTOUR_DO(JJ) = "LakeContour_DO_Seg" // ADJUSTL(trim(ICHAR3)) // ".csv"
                        open(LAKE_RIVER_CONTOUR + 10 + JJ - 1, FILE=FILE_LAKE_CONTOUR_DO(JJ), status="unknown")
                        if (LAKE_CONTOUR_FORMAT == 1) then
                            write(LAKE_RIVER_CONTOUR + 10 + JJ - 1, *) "JDAY,ELEVATION(m),DisslvedOxygen(mg/l)"
                        else
                            write(LAKE_RIVER_CONTOUR + 10 + JJ - 1, '("TIME,",*(F8.2,","))') ((EL(K, LAKE_CONTOUR_SEG(JJ)) + EL(K, LAKE_CONTOUR_SEG(JJ)))*0.5, K = 2, KB(LAKE_CONTOUR_SEG(JJ)))
                        end if
                    end if
                end do

                do JJ = 1, NUM_RIVER_CONTOUR
                    write(ICHAR3, "(I3)") RIVER_CONTOUR_BR1(JJ)
                    FILE_RIVER_CONTOUR_T(JJ) = "RiverContour_T_Br" // ADJUSTL(trim(ICHAR3)) // ".csv"
                    open(LAKE_RIVER_CONTOUR + 20 + JJ - 1, FILE=FILE_RIVER_CONTOUR_T(JJ), status="unknown")
                    if (RIVER_CONTOUR_FORMAT == 1) then
                        write(LAKE_RIVER_CONTOUR + 20 + JJ - 1, *) "JDAY,DISTANCE(M),TEMPERATURE(C)"
                    else
                        write(LAKE_RIVER_CONTOUR + 20 + JJ - 1, '("TIME,",*(F8.2,","))') ((X1(I), I = US(JB), DS(JB)), JB = RIVER_CONTOUR_BR1(JJ), RIVER_CONTOUR_BR2(JJ))
                    end if
                    if (OXYGEN_DEMAND) then
                        write(ICHAR3, "(I3)") LAKE_CONTOUR_SEG(JJ)
                        FILE_RIVER_CONTOUR_DO(JJ) = "RiverContour_DO_Br" // ADJUSTL(trim(ICHAR3)) // ".csv"
                        open(LAKE_RIVER_CONTOUR + 30 + JJ - 1, FILE=FILE_RIVER_CONTOUR_DO(JJ), status="unknown")
                        if (RIVER_CONTOUR_FORMAT == 1) then
                            write(LAKE_RIVER_CONTOUR + 30 + JJ - 1, *) "JDAY,ELEVATION(M),DisslvedOxygen(mg/l)"
                        else
                            write(LAKE_RIVER_CONTOUR + 30 + JJ - 1, '("TIME,",*(F8.2,","))') ((EL(K, LAKE_CONTOUR_SEG(JJ)) + EL(K, LAKE_CONTOUR_SEG(JJ)))*0.5, K = 2, KB(LAKE_CONTOUR_SEG(JJ)))
                        end if
                    end if
                end do
            end if
        end if

        do JW = 1, NWB
            if (SNAPSHOT(JW)) then
                open(SNP(JW), FILE=SNPFN(JW), STATUS="UNKNOWN")
            end if
!IF (VECTOR(JW))      OPEN (VPL(JW),FILE=VPLFN(JW),STATUS='UNKNOWN')
            if (PROFILE(JW)) then
                open(PRF(JW), FILE=PRFFN(JW), STATUS="UNKNOWN")
            end if
            if (SPREADSHEET(JW)) then
                open(SPR(JW), FILE=SPRFN(JW), STATUS="UNKNOWN")
            end if
            if (SPRC(JW) == "     ONV") then
                open(SPRV(JW), FILE=SPRVFN(JW), STATUS="UNKNOWN")
            end if ! SW 9/28/2018
            if (CONTOUR(JW)) then
                open(CPL(JW), FILE=CPLFN(JW), STATUS="UNKNOWN")
            end if
            if (FLUX(JW)) then
                open(FLX(JW), FILE=FLXFN(JW), STATUS="UNKNOWN")
                write(SEGNUM, "(I0)") JW
                SEGNUM = ADJUSTL(SEGNUM)
                L = LEN_TRIM(SEGNUM)
                open(FLX2(JW), FILE="kflux_wb" // SEGNUM(1:L) // ".csv", STATUS="UNKNOWN") ! SW 7/1/2019
                write(FLX2(JW), '("JDAY,  ELTM,",*(A,","))') (KFNAME2(KFCN(JF, JW)), JF = 1, NAF(JW))
            end if

!**** Output files

            if (PROFILE(JW) .and. iprf(1, 1) /= -1) then ! SW 4/1/2016
                TTIME = TMSTRT
                do while (TTIME <= TMEND)
                    NDSP = NDSP + 1
                    TTIME = TTIME + PRFF(PRFDP(JW), JW)
                    if (TTIME >= PRFD(PRFDP(JW) + 1, JW)) then
                        PRFDP(JW) = PRFDP(JW) + 1
                    end if
                end do
                PRFDP(JW) = 1
                write(PRF(JW), "(A)") TITLE
                write(PRF(JW), "(8I8,L2)") KMX, NIPRF(JW), NDSP, NCT, NDC, NAC + NACD(JW) + 1, PRFDP(JW), KTWB(JW), CONSTITUENTS
                write(PRF(JW), "(10I8)") IPRF(1:NIPRF(JW), JW)
                write(PRF(JW), "(20(1X,A))") " ON", CPRWBC(:, JW)(6:8), CDWBC(:, JW)(6:8)
                write(PRF(JW), "(2A)") "Temperature,oC                            ", ADJUSTL(CNAME), ADJUSTL(CDNAME)
                write(PRF(JW), "(20I4)") 1, CN(1:NAC) + 1, CDN(1:NACD(JW), JW) + NCT + 1
                write(PRF(JW), "(10F8.0)") 1.0, CMULT, CDMULT
                write(PRF(JW), "(20I4)") (KB(IPRF(I, JW)), I = 1, NIPRF(JW)) ! KB(IPRF(1:NIPRF(JW),JW))
                write(PRF(JW), "(10F8.2)") H
                do JP = 1, NIPRF(JW)
                    NRS = KB(IPRF(JP, JW)) - KTWB(JW) + 1
                    write(PRF(JW), "(A8,I4/(8(F10.2)))") "TEMP    ", NRS, (T2(K, IPRF(JP, JW)), K = KTWB(JW), KB(IPRF(JP, JW)))
                end do
                do JC = 1, NAC
                    if (PRINT_CONST(CN(JC), JW)) then
                        do JP = 1, NIPRF(JW)
                            NRS = KB(IPRF(JP, JW)) - KTWB(JW) + 1
                            write(PRF(JW), "(A,I4/(8(E13.6,1x)))") ADJUSTL(CNAME2(CN(JC))), NRS, (C2(K, IPRF(JP, JW), CN(JC))*CMULT(CN(JC)), K = KTWB(JW), KB(IPRF(JP, JW)))
                        end do
                    end if
                end do
                do JD = 1, NACD(JW)
                    do JP = 1, NIPRF(JW)
                        NRS = KB(IPRF(JP, JW)) - KTWB(JW) + 1
                        write(PRF(JW), "(A,I4/(8(E13.6,1x)))") ADJUSTL(CDNAME2(CDN(JD, JW))), NRS, (CD(K, IPRF(JP, JW), CDN(JD, JW))*CDMULT(CDN(JD, JW)), K = KTWB(JW), KB(IPRF(JP, JW)))
                    end do
                end do
            end if
            if (SPREADSHEET(JW)) then
                do J = 1, NISPR(JW)
                    write(SEGNUM, "(I0)") ISPR(J, JW)
                    SEGNUM = ADJUSTL(SEGNUM)
                    L = LEN_TRIM(SEGNUM)
                    SEG(J) = "Seg_" // SEGNUM(1:L)
                end do
                if (SPRC(JW) == "     ONV") then
                    write(SPRV(JW), '(A,A,*(A7,","))') "Constituent,", "Julian_day,", (SEG(J), J = 1, NISPR(JW)) ! SW 9/28/2018
                end if
                write(SPR(JW), '(A,A,A,*("Elevation,",A7,","))') "Constituent,", "Julian_day,", "Depth,", (SEG(J), J = 1, NISPR(JW))

            end if
            if (CONTOUR(JW)) then
                if (TECPLOT(JW) /= "      ON") then
                    write(CPL(JW), "(A)") TITLE
                    write(CPL(JW), "(8(I8,2X))") NBR
                    write(CPL(JW), "(8(I8,2X))") IMX, KMX
                    do JB = BS(JW), BE(JW)
                        write(CPL(JW), "(9(I8,2X))") US(JB), DS(JB)
                        write(CPL(JW), "(9(I8,2X))") KB(US(JB):DS(JB))
                    end do
                    write(CPL(JW), "(8(E13.6,2X))") DLX
                    write(CPL(JW), "(8(E13.6,2X))") H
                    write(CPL(JW), "(8(I8,2X))") NAC
                    write(CPL(JW), "(A)") (CNAME1(CN(JN)), JN = 1, NAC) ! SW 3/1/2017
                else
!c calculating longitudinal distance of segments
!           IF(JW == 1)DIST=0.0
!           do jb=BS(JW),BE(JW)
!               x1(us(jb))=dist+dlx(us(jb))/2.
!               DO I=US(JB)+1,DS(JB)
!                   DIST=DIST+(DLX(I)+dlx(i-1))/2.0
!                   X1(I)=DIST
!               END DO
!               DIST=DIST+DLX(DS(JB))
!               X1(DS(JB)+1)=DIST
!           ENDDO
                    write(CPL(JW), *) 'TITLE="CE-QUAL-W2"'
                    if (HABTATC == "      ON") then
                        if (NAC == 0) then
                            write(CPL(JW), 19231)
                        else
                            if (NACD(JW) == 0) then
                                write(CPL(JW), 19232) (CNAME2(CN(JN)), JN = 1, NAC) !WRITE (CPL(JW),19233)(CNAME2(CN(JN)),JN=1,NAC)    SW 1/17/17
                            else
                                if (NACD(JW) /= 0) then
                                    write(CPL(JW), 19233) (CNAME2(CN(JN)), JN = 1, NAC), (CDNAME2(CDN(JD, JW)), JD = 1, NACD(JW))
                                end if
                            end if
                        end if
                    else
                        if (NAC == 0) then
                            write(CPL(JW), 19230)
                        else
                            if (NACD(JW) == 0) then
                                write(CPL(JW), 19234) (CNAME2(CN(JN)), JN = 1, NAC) !    1/17/17          !WRITE (CPL(JW),19233)(CNAME2(CN(JN)),JN=1,NAC)    SW 1/17/17
                            else
                                if (NACD(JW) /= 0) then
                                    write(CPL(JW), 19235) (CNAME2(CN(JN)), JN = 1, NAC), (CDNAME2(CDN(JD, JW)), JD = 1, NACD(JW)) !    1/17/17      !WRITE (CPL(JW),19234)(CNAME2(CN(JN)),JN=1,NAC)  SW 9/28/13
                                end if
                            end if
                        end if
                    end if
                    19230 format('VARIABLES="Distance, m","Elevation, m","U(m/s)","W(m/s)","T(C)","RHO" ')
                    19231 format('VARIABLES="Distance, m","Elevation, m","U(m/s)","W(m/s)","T(C)","RHO", "HABITAT" ')
                    19232 format('VARIABLES="Distance, m","Elevation, m","U(m/s)","W(m/s)","T(C)","RHO", "HABITAT" ',*(',"',A8,'"'))
                    19233 format('VARIABLES="Distance, m","Elevation, m","U(m/s)","W(m/s)","T(C)","RHO", "HABITAT" ',*(',"',A8,'"'),*(',"',A8,'"'))
                    19234 format('VARIABLES="Distance, m","Elevation, m","U(m/s)","W(m/s)","T(C)","RHO" ',*(',"',A8,'"')) ! SW 9/28/13
                    19235 format('VARIABLES="Distance, m","Elevation, m","U(m/s)","W(m/s)","T(C)","RHO" ',*(',"',A8,'"'),*(',"',A8,'"')) ! SW 9/28/13
                end if
            end if

!IF (VECTOR(JW)) THEN     *** DSI
!  WRITE (VPL(JW),*)  TITLE
!  WRITE (VPL(JW),*)  H,KB,US,DS,DLX
!END IF

        end do

! ***
        if (TIME_SERIES) then
            L1 = SCAN(TSRFN1, ".", BACK=.true.) ! SW 8/22/14
            do J = 1, NIKTSR
                write(SEGNUM, "(I0)") ITSR(J)
                SEGNUM = ADJUSTL(SEGNUM)
                L = LEN_TRIM(SEGNUM)
                write(SEGNUM2, "(I0)") J
                SEGNUM2 = ADJUSTL(SEGNUM2)
                L2 = LEN_TRIM(SEGNUM2)
                TSRFN = TSRFN1(1:L1 - 1) // "_" // SEGNUM2(1:L2) // "_seg" // SEGNUM(1:L) // "." // TSRFN1(L1 + 1:L1 + 4) !'.opt'
                open(TSR(J), FILE=TSRFN, STATUS="UNKNOWN")
!WRITE (TSR(J),'(A)') (TITLE(I),I=1,11)
                I = ITSR(J) ! SR 5/10/05
                do JW = 1, NWB
                    if (I >= US(BS(JW)) .and. I <= DS(BE(JW))) then
                        exit
                    end if
                end do

                if (ICE_COMPUTATION) then
                    if (SEDIMENT_CALC(JW)) then
                        write(TSR(J), '(*(A,","))') "JDAY", "DLT(s)", "ELWS(m)", "T2(C)", "U(ms-1)", "Q(m3s-1)", "SRON(Wm-2)", "EXT(m-1)", "DEPTH(m)", "WIDTH(m)", "SHADE", "ICETH(m)", "Tvolavg(C)", "NetRad(Wm-2)", "SWSolar(Wm-2)", "LWRad(Wm-2)", "BackRad(Wm-2)", "EvapF(Wm-2)", "ConducF(Wm-2)", "ReaerationCoeff(day-1)", (CNAME2(CN(JC)), JC = 1, NAC), ("     EPI", JE = 1, NEP), ("     MAC", JM = 1, NMC), "     SED(Organic matter in sediments g/m3)", "    SEDP(OrgP in sediments gP/m3)", "    SEDN(OrgN in sediments gN/m3)", "    SEDC(OrgCSediments gC/m3)", (CDNAME2(CDN(JD, JW)), JD = 1, NACD(JW)), (KFNAME2(KFCN(JF, JW)), JF = 1, NAF(JW)), ("PLIM_" // ADJUSTL(CNAME2(NAS + JA - 1)), JA = 1, NAL), ("NLIM_" // ADJUSTL(CNAME2(NAS + JA - 1)), JA = 1, NAL), ("LLIM_" // ADJUSTL(CNAME2(NAS + JA - 1)), JA = 1, NAL)
                    else
                        write(TSR(J), '(*(A,","))') "JDAY", "DLT(s)", "ELWS(m)", "T2(C)", "U(ms-1)", "Q(m3s-1)", "SRON(Wm-2)", "EXT(m-1)", "DEPTH(m)", "WIDTH(m)", "SHADE", "ICETH(m)", "Tvolavg(C)", "NetRad(Wm-2)", "SWSolar(Wm-2)", "LWRad(Wm-2)", "BackRad(Wm-2)", "EvapF(Wm-2)", "ConducF(Wm-2)", "ReaerationCoeff(day-1)", (CNAME2(CN(JC)), JC = 1, NAC), ("     EPI", JE = 1, NEP), ("     MAC", JM = 1, NMC), (CDNAME2(CDN(JD, JW)), JD = 1, NACD(JW)), (KFNAME2(KFCN(JF, JW)), JF = 1, NAF(JW)), ("PLIM_" // ADJUSTL(CNAME2(NAS + JA - 1)), JA = 1, NAL), ("NLIM_" // ADJUSTL(CNAME2(NAS + JA - 1)), JA = 1, NAL), ("LLIM_" // ADJUSTL(CNAME2(NAS + JA - 1)), JA = 1, NAL) ! SW 10/20/15
                    end if
                else
                    if (SEDIMENT_CALC(JW)) then !mlm 7/25/06
                        write(TSR(J), '(*(A,","))') "JDAY", "DLT(s)", "ELWS(m)", "T2(C)", "U(ms-1)", "Q(m3s-1)", "SRON(Wm-2)", "EXT(m-1)", "DEPTH(m)", "WIDTH(m)", "SHADE", "Tvolavg(C)", "NetRad(Wm-2)", "SWSolar(Wm-2)", "LWRad(Wm-2)", "BackRad(Wm-2)", "EvapF(Wm-2)", "ConducF(Wm-2)", "ReaerationCoeff(day-1)", (CNAME2(CN(JC)), JC = 1, NAC), ("     EPI", JE = 1, NEP), ("     MAC", JM = 1, NMC), "     SED(Organic matter in sediments g/m3)", "    SEDP(OrgP in sediments gP/m3)", "    SEDN(OrgN in sediments gN/m3)", "    SEDC(OrgCSediments gC/m3)", (CDNAME2(CDN(JD, JW)), JD = 1, NACD(JW)), (KFNAME2(KFCN(JF, JW)), JF = 1, NAF(JW)), ("PLIM_" // ADJUSTL(CNAME2(NAS + JA - 1)), JA = 1, NAL), ("NLIM_" // ADJUSTL(CNAME2(NAS + JA - 1)), JA = 1, NAL), ("LLIM_" // ADJUSTL(CNAME2(NAS + JA - 1)), JA = 1, NAL)
                    else
                        write(TSR(J), '(*(A,","))') "JDAY", "DLT(s)", "ELWS(m)", "T2(C)", "U(ms-1)", "Q(m3s-1)", "SRON(Wm-2)", "EXT(m-1)", "DEPTH(m)", "WIDTH(m)", "SHADE", "Tvolavg(C)", "NetRad(Wm-2)", "SWSolar(Wm-2)", "LWRad(Wm-2)", "BackRad(Wm-2)", "EvapF(Wm-2)", "ConducF(Wm-2)", "ReaerationCoeff(day-1)", (CNAME2(CN(JC)), JC = 1, NAC), ("     EPI", JE = 1, NEP), ("     MAC", JM = 1, NMC), (CDNAME2(CDN(JD, JW)), JD = 1, NACD(JW)), (KFNAME2(KFCN(JF, JW)), JF = 1, NAF(JW)), ("PLIM_" // ADJUSTL(CNAME2(NAS + JA - 1)), JA = 1, NAL), ("NLIM_" // ADJUSTL(CNAME2(NAS + JA - 1)), JA = 1, NAL), ("LLIM_" // ADJUSTL(CNAME2(NAS + JA - 1)), JA = 1, NAL)
                    end if
                end if

            end do
        end if
        if (ALGAE_TOXIN) then
            if (ATOX_DEBUG == "ON") then
                open(ATOXIN_DEBUG_FN, FILE="algae_toxin_debug.csv", STATUS="UNKNOWN")
                write(ATOXIN_DEBUG_FN, '(A,<NUMATOXINS>("EX_TOXIN_",I1,","),<NUMATOXINS>("IN_TOXIN_",I1,","),<NUMATOXINS>("CTESS_TOXIN_",I1,","),<NAL>("ALG_",I1,","))') "JDAY,K,I,", (J, J = 1, numatoxins), (J, J = 1, numatoxins), (J, J = 1, numatoxins), (J, J = 1, NAL)
            end if
        end if
! BIOENERGETICS mlm
        if (BIOEXP) then
            L1 = SCAN(BIOFN, ".")
            do J = 1, NIBIO
                write(SEGNUM, "(I0)") J
                SEGNUM = ADJUSTL(SEGNUM)
                L = LEN_TRIM(SEGNUM)
                BIOFN = BIOFN(1:L1 - 1) // "_" // SEGNUM(1:L) // ".opt"
                open(BIOEXPFN(J), FILE=BIOFN, STATUS="UNKNOWN")
                write(bioexpfn(J), "(12A8)") (adjustr(bhead(ii)), ii = 1, 11)
            end do
            L1 = SCAN(WEIGHTFN, ".")
            do J = 1, NIBIO
                write(SEGNUM, "(I0)") J
                SEGNUM = ADJUSTL(SEGNUM)
                L = LEN_TRIM(SEGNUM)
                WEIGHTFN = WEIGHTFN(1:L1 - 1) // "_" // SEGNUM(1:L) // ".opt"
                open(WEIGHTNUM(J), FILE=WEIGHTFN, STATUS="UNKNOWN")
                write(weightnum(J), '(50(A10,","))') adjustr("jday"), (adjustr(cname2(cn(jc))), jc = 1, nac), adjustr("pH"), adjustr("TP")
            end do

        end if
        if (DOWNSTREAM_OUTFLOW) then
            L1 = SCAN(WDOFN, ".", BACK=.true.) ! SW 8/22/14 CHECKS FROM RIGHTHAND SIDE NOT LEFTHANDSIDE  ! SW 8/3/2018
            do JWD = 1, NIWDO
                write(SEGNUM, "(I0)") IWDO(JWD)
                SEGNUM = ADJUSTL(SEGNUM)
                L = LEN_TRIM(SEGNUM)
                open(WDO(JWD, 1), FILE="qwo_" // SEGNUM(1:L) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !      OPEN  (WDO(JWD,1),FILE='qwo_'//SEGNUM(1:L)//'.opt',STATUS='UNKNOWN')
                open(WDO(JWD, 2), FILE="two_" // SEGNUM(1:L) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                write(WDO(JWD, 1), "(A,I0/A/A)") "$Flow file for segment ", IWDO(JWD), "To the right of the sum of flows are individual flows starting with QWD then QSTR", "JDAY,QWD(m3s-1),"
                write(WDO(JWD, 2), "(A,I0/A/A)") "$Temperature file for segment ", IWDO(JWD), "To the right of the sum of temperatures are individual temperatures starting with QWD then QSTR", "JDAY,T(C),"
                do JW = 1, NWB
                    if (IWDO(JWD) >= US(BS(JW)) .and. IWDO(JWD) <= DS(BE(JW))) then
                        exit
                    end if
                end do
                if (CONSTITUENTS) then
                    open(WDO(JWD, 3), FILE="cwo_" // SEGNUM(1:L) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                    write(WDO(JWD, 3), '(A,I0//(*(A,",")))') "$Concentration file for segment ", IWDO(JWD), "JDAY", (CNAME2(CN(J)), J = 1, NAC) !CNAME2(CN(1:NAC))
                end if
                if (DERIVED_CALC) then
                    open(WDO(JWD, 4), FILE="dwo_" // SEGNUM(1:L) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                    write(WDO(JWD, 4), '(A,I0//(*(A,",")))') "Derived constituent file for segment ", IWDO(JWD), "JDAY", (CDNAME2(CDN(J, JW)), J = 1, NACD(JW)) !CDNAME2(CDN(1:NACD(JW),JW))
                end if
            end do
        end if

        if (DOWNSTREAM_OUTFLOW) then
            JFILE = 0
            do JWD = 1, NIWDO

! Determine the # of withdrawals at the WITH SEG
                do JB = 1, NBR ! structures
                    if (IWDO(JWD) == DS(JB) .and. NSTR(JB) /= 0) then
                        do JS = 1, NSTR(JB)
                            JFILE = JFILE + 1
                            write(SEGNUM2, "(I0)") IWDO(JWD)
                            SEGNUM2 = ADJUSTL(SEGNUM2)
                            L2 = LEN_TRIM(SEGNUM2)
                            TITLEWITH2 = "$STR WITHDRAWAL AT SEG" // SEGNUM2(1:L2)
                            write(SEGNUM, "(F10.0)") ESTR(JS, JB)
                            SEGNUM = ADJUSTL(SEGNUM)
                            L = LEN_TRIM(SEGNUM)
                            TITLEWITH = ADJUSTL(TITLEWITH2) // " ELEV OF WITHDRAWAL CENTERLINE:" // SEGNUM(1:L)
                            write(SEGNUM, "(I0)") JS
                            SEGNUM = ADJUSTL(SEGNUM)
                            L = LEN_TRIM(SEGNUM)
                            WDO2(JFILE, 1) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 1), FILE="qwo_str" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                            write(WDO2(JFILE, 1), "(A//A)") TITLEWITH, "JDAY,QWD(m3s-1)"
                            WDO2(JFILE, 2) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 2), FILE="two_str" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                            write(WDO2(JFILE, 2), "(A//A)") TITLEWITH, "JDAY,T(C)"
                            if (CONSTITUENTS) then
                                WDO2(JFILE, 3) = NUNIT;                                 NUNIT = NUNIT + 1
                                open(WDO2(JFILE, 3), FILE="cwo_str" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                                write(WDO2(JFILE, 3), '(A//(*(A,",")))') TITLEWITH, "JDAY", (CNAME2(CN(J)), J = 1, NAC) !CNAME2(CN(1:NAC))
                            end if
                            if (DERIVED_CALC) then
                                WDO2(JFILE, 4) = NUNIT;                                 NUNIT = NUNIT + 1
                                open(WDO2(JFILE, 4), FILE="dwo_str" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                                do JW = 1, NWB
                                    if (IWDO(JWD) >= US(BS(JW)) .and. IWDO(JWD) <= DS(BE(JW))) then
                                        exit
                                    end if
                                end do
                                write(WDO2(JFILE, 4), '(A//(*(A,",")))') TITLEWITH, "JDAY", (CDNAME2(CDN(J, JW)), J = 1, NACD(JW)) !CDNAME2(CDN(1:NACD(JW),JW))
                            end if
                        end do
                    end if
                end do

                do JS = 1, NWD ! withdrawals
                    if (IWDO(JWD) == IWD(JS)) then
                        JFILE = JFILE + 1
                        write(SEGNUM2, "(I0)") IWDO(JWD)
                        SEGNUM2 = ADJUSTL(SEGNUM2)
                        L2 = LEN_TRIM(SEGNUM2)
                        TITLEWITH2 = "$WITHDRAWAL AT SEG" // SEGNUM2(1:L2)
                        write(SEGNUM, "(F10.0)") EWD(JS)
                        SEGNUM = ADJUSTL(SEGNUM)
                        L = LEN_TRIM(SEGNUM)
                        TITLEWITH = ADJUSTL(TITLEWITH2) // " ELEV OF WITHDRAWAL CENTERLINE:" // SEGNUM(1:L)
                        write(SEGNUM, "(I0)") JS
                        SEGNUM = ADJUSTL(SEGNUM)
                        L = LEN_TRIM(SEGNUM)
                        WDO2(JFILE, 1) = NUNIT;                         NUNIT = NUNIT + 1
                        open(WDO2(JFILE, 1), FILE="qwo_wd" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                        write(WDO2(JFILE, 1), "(A//A)") TITLEWITH, "JDAY,QWD(m3s-1)"
                        WDO2(JFILE, 2) = NUNIT;                         NUNIT = NUNIT + 1
                        open(WDO2(JFILE, 2), FILE="two_wd" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                        write(WDO2(JFILE, 2), "(A//A)") TITLEWITH, "JDAY,T(C)"
                        if (CONSTITUENTS) then
                            WDO2(JFILE, 3) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 3), FILE="cwo_wd" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                            write(WDO2(JFILE, 3), '(A//(*(A,",")))') TITLEWITH, "JDAY", (CNAME2(CN(J)), J = 1, NAC) !CNAME2(CN(1:NAC))
                        end if
                        if (DERIVED_CALC) then
                            WDO2(JFILE, 4) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 4), FILE="dwo_wd" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                            do JW = 1, NWB
                                if (IWDO(JWD) >= US(BS(JW)) .and. IWDO(JWD) <= DS(BE(JW))) then
                                    exit
                                end if
                            end do
                            write(WDO2(JFILE, 4), '(A//(*(A,",")))') TITLEWITH, "JDAY", (CDNAME2(CDN(J, JW)), J = 1, NACD(JW)) !CDNAME2(CDN(1:NACD(JW),JW))
                        end if
                    end if
                end do

                do JS = 1, NSP ! spillways
                    if (IWDO(JWD) == IUSP(JS)) then
                        JFILE = JFILE + 1
                        write(SEGNUM2, "(I0)") IWDO(JWD)
                        SEGNUM2 = ADJUSTL(SEGNUM2)
                        L2 = LEN_TRIM(SEGNUM2)
                        TITLEWITH2 = "$SPILLWAY WITHDRAWAL AT SEG" // SEGNUM2(1:L2)
                        write(SEGNUM, "(F10.0)") ESP(JS)
                        SEGNUM = ADJUSTL(SEGNUM)
                        L = LEN_TRIM(SEGNUM)
                        TITLEWITH = ADJUSTL(TITLEWITH2) // " ELEV OF WITHDRAWAL CENTERLINE:" // SEGNUM(1:L)
                        write(SEGNUM, "(I0)") JS
                        SEGNUM = ADJUSTL(SEGNUM)
                        L = LEN_TRIM(SEGNUM)
                        WDO2(JFILE, 1) = NUNIT;                         NUNIT = NUNIT + 1
                        open(WDO2(JFILE, 1), FILE="qwo_sp" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                        write(WDO2(JFILE, 1), "(A//A)") TITLEWITH, "JDAY,QWD(m3s-1)"
                        WDO2(JFILE, 2) = NUNIT;                         NUNIT = NUNIT + 1
                        open(WDO2(JFILE, 2), FILE="two_sp" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                        write(WDO2(JFILE, 2), "(A//A)") TITLEWITH, "JDAY,T(C)"
                        if (CONSTITUENTS) then
                            WDO2(JFILE, 3) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 3), FILE="cwo_sp" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                            write(WDO2(JFILE, 3), '(A//(*(A,",")))') TITLEWITH, "JDAY", (CNAME2(CN(J)), J = 1, NAC) !CNAME2(CN(1:NAC))
                        end if
                        if (DERIVED_CALC) then
                            WDO2(JFILE, 4) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 4), FILE="dwo_sp" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                            do JW = 1, NWB
                                if (IWDO(JWD) >= US(BS(JW)) .and. IWDO(JWD) <= DS(BE(JW))) then
                                    exit
                                end if
                            end do
                            write(WDO2(JFILE, 4), '(A//(*(A,",")))') TITLEWITH, "JDAY", (CDNAME2(CDN(J, JW)), J = 1, NACD(JW)) !CDNAME2(CDN(1:NACD(JW),JW))
                        end if
                    end if
                end do


                do JS = 1, NPU ! pumps
                    if (IWDO(JWD) == IUPU(JS)) then
                        JFILE = JFILE + 1
                        write(SEGNUM2, "(I0)") IWDO(JWD)
                        SEGNUM2 = ADJUSTL(SEGNUM2)
                        L2 = LEN_TRIM(SEGNUM2)
                        TITLEWITH2 = "$PUMP WITHDRAWAL AT SEG" // SEGNUM2(1:L2)
                        write(SEGNUM, "(F10.0)") EPU(JS)
                        SEGNUM = ADJUSTL(SEGNUM)
                        L = LEN_TRIM(SEGNUM)
                        TITLEWITH = ADJUSTL(TITLEWITH2) // " ELEV OF WITHDRAWAL CENTERLINE:" // SEGNUM(1:L)
                        write(SEGNUM, "(I0)") JS
                        SEGNUM = ADJUSTL(SEGNUM)
                        L = LEN_TRIM(SEGNUM)
                        WDO2(JFILE, 1) = NUNIT;                         NUNIT = NUNIT + 1
                        open(WDO2(JFILE, 1), FILE="qwo_pmp" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                        write(WDO2(JFILE, 1), "(A//A)") TITLEWITH, "JDAY,QWD(m3s-1)"
                        WDO2(JFILE, 2) = NUNIT;                         NUNIT = NUNIT + 1
                        open(WDO2(JFILE, 2), FILE="two_pmp" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                        write(WDO2(JFILE, 2), "(A//A)") TITLEWITH, "JDAY,T(C)"
                        if (CONSTITUENTS) then
                            WDO2(JFILE, 3) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 3), FILE="cwo_pmp" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                            write(WDO2(JFILE, 3), '(A//(*(A,",")))') TITLEWITH, "JDAY", (CNAME2(CN(J)), J = 1, NAC) !CNAME2(CN(1:NAC))
                        end if
                        if (DERIVED_CALC) then
                            WDO2(JFILE, 4) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 4), FILE="dwo_pmp" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                            do JW = 1, NWB
                                if (IWDO(JWD) >= US(BS(JW)) .and. IWDO(JWD) <= DS(BE(JW))) then
                                    exit
                                end if
                            end do
                            write(WDO2(JFILE, 4), '(A//(*(A,",")))') TITLEWITH, "JDAY", (CDNAME2(CDN(J, JW)), J = 1, NACD(JW)) !CDNAME2(CDN(1:NACD(JW),JW))
                        end if
                    end if
                end do


                do JS = 1, NPI ! pipes
                    if (IWDO(JWD) == IUPI(JS)) then
                        JFILE = JFILE + 1
                        write(SEGNUM2, "(I0)") IWDO(JWD)
                        SEGNUM2 = ADJUSTL(SEGNUM2)
                        L2 = LEN_TRIM(SEGNUM2)
                        TITLEWITH2 = "$PIPE WITHDRAWAL AT SEG" // SEGNUM2(1:L2)
                        write(SEGNUM, "(F10.0)") EUPI(JS)
                        SEGNUM = ADJUSTL(SEGNUM)
                        L = LEN_TRIM(SEGNUM)
                        TITLEWITH = ADJUSTL(TITLEWITH2) // " ELEV OF WITHDRAWAL CENTERLINE:" // SEGNUM(1:L)
                        write(SEGNUM, "(I0)") JS
                        SEGNUM = ADJUSTL(SEGNUM)
                        L = LEN_TRIM(SEGNUM)
                        WDO2(JFILE, 1) = NUNIT;                         NUNIT = NUNIT + 1
                        open(WDO2(JFILE, 1), FILE="qwo_pipe" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                        write(WDO2(JFILE, 1), "(A//A)") TITLEWITH, "JDAY,QWD(m3s-1)"
                        WDO2(JFILE, 2) = NUNIT;                         NUNIT = NUNIT + 1
                        open(WDO2(JFILE, 2), FILE="two_pipe" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                        write(WDO2(JFILE, 2), "(A//A)") TITLEWITH, "JDAY,T(C)"
                        if (CONSTITUENTS) then
                            WDO2(JFILE, 3) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 3), FILE="cwo_pipe" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                            write(WDO2(JFILE, 3), '(A//(*(A,",")))') TITLEWITH, "JDAY", (CNAME2(CN(J)), J = 1, NAC) !CNAME2(CN(1:NAC))
                        end if
                        if (DERIVED_CALC) then
                            WDO2(JFILE, 4) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 4), FILE="dwo_pipe" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                            do JW = 1, NWB
                                if (IWDO(JWD) >= US(BS(JW)) .and. IWDO(JWD) <= DS(BE(JW))) then
                                    exit
                                end if
                            end do
                            write(WDO2(JFILE, 4), '(A//(*(A,",")))') TITLEWITH, "JDAY", (CDNAME2(CDN(J, JW)), J = 1, NACD(JW)) !CDNAME2(CDN(1:NACD(JW),JW))
                        end if
                    end if
                end do

                do JS = 1, NGT ! gates
                    if (IWDO(JWD) == IUGT(JS)) then
                        JFILE = JFILE + 1
                        write(SEGNUM2, "(I0)") IWDO(JWD)
                        SEGNUM2 = ADJUSTL(SEGNUM2)
                        L2 = LEN_TRIM(SEGNUM2)
                        TITLEWITH2 = "$GATE WITHDRAWAL AT SEG" // SEGNUM2(1:L2)
                        write(SEGNUM, "(F10.0)") EGT(JS)
                        SEGNUM = ADJUSTL(SEGNUM)
                        L = LEN_TRIM(SEGNUM)
                        TITLEWITH = ADJUSTL(TITLEWITH2) // " ELEV OF WITHDRAWAL CENTERLINE:" // SEGNUM(1:L)
                        write(SEGNUM, "(I0)") JS
                        SEGNUM = ADJUSTL(SEGNUM)
                        L = LEN_TRIM(SEGNUM)
                        WDO2(JFILE, 1) = NUNIT;                         NUNIT = NUNIT + 1
                        open(WDO2(JFILE, 1), FILE="qwo_gate" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                        write(WDO2(JFILE, 1), "(A//A)") TITLEWITH, "JDAY,QWD(m3s-1)"
                        WDO2(JFILE, 2) = NUNIT;                         NUNIT = NUNIT + 1
                        open(WDO2(JFILE, 2), FILE="two_gate" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                        write(WDO2(JFILE, 2), "(A//A)") TITLEWITH, "JDAY,T(C)"
                        if (CONSTITUENTS) then
                            WDO2(JFILE, 3) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 3), FILE="cwo_gate" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                            write(WDO2(JFILE, 3), '(A//(*(A,",")))') TITLEWITH, "JDAY", (CNAME2(CN(J)), J = 1, NAC) !CNAME2(CN(1:NAC))
                        end if
                        if (DERIVED_CALC) then
                            WDO2(JFILE, 4) = NUNIT;                             NUNIT = NUNIT + 1
                            open(WDO2(JFILE, 4), FILE="dwo_gate" // SEGNUM(1:L) // "_seg" // SEGNUM2(1:L2) // WDOFN(L1:L1 + 4), STATUS="UNKNOWN") !'.opt',STATUS='UNKNOWN')
                            do JW = 1, NWB
                                if (IWDO(JWD) >= US(BS(JW)) .and. IWDO(JWD) <= DS(BE(JW))) then
                                    exit
                                end if
                            end do
                            write(WDO2(JFILE, 4), '(A//(*(A,",")))') TITLEWITH, "JDAY", (CDNAME2(CDN(J, JW)), J = 1, NACD(JW)) !CDNAME2(CDN(1:NACD(JW),JW))
                        end if
                    end if
                end do

            end do
        end if

! OUTPUT FLOW LOADING AND POLLUTANT LOADING DEBUGGING INFORMATION OUTPUT FREQUENCY AT OUTPUT OF CPL
        do JW = 1, NWB
            if (FLOWBALC == "      ON") then
                open(FLOWBFN, FILE="flowbal.csv", STATUS="UNKNOWN")
                if (VOLUME_BALANCE(JW)) then
                    write(FLOWBFN, "(a102)") "JDAY,WB,VOLIN(m3),VOLPR(m3),VOLOUT(m3),VOLWD(m3),VOLEV(m3),VOLDT(m3),VOLTRB(m3),VOLICE(m3),%VOLerror,"
                else
                    write(FLOWBFN, "(a93)") "JDAY,WB,VOLIN(m3),VOLPR(m3),VOLOUT(m3),VOLWD(m3),VOLEV(m3),VOLDT(m3),VOLTRB(m3),VOLICE(m3),"
                end if
                exit
            end if
        end do
        do JW = 1, NWB
            if (NPBALC == "      ON") then
                open(MASSBFN, FILE="massbal.csv", STATUS="UNKNOWN")
                if (SEDIMENT_DIAGENESIS) then
                    write(MASSBFN, "(A574)") "JDAY,WB,TP-Waterbody(kg),TP-Sediment(kg),TP-Plants(kg),OutflowTP(kg),TributaryTP(kg),DistributedTributaryTP(kg),WithdrawalTP(kg),PrecipitationTP(kg),InflowTP(kg),AtmosphericDepositionP(kg),SED+SOD_PRelease(kg),PFluxtoSediments(kg),SedimentDiagenesisPFlux(kg),TN-Waterbody(kg),TN-Sediment(kg),TN-Plants(kg),OutflowTN(kg),TributaryTN(kg),DistributedTributaryTN(kg),WithdrawalTN(kg),PrecipitationTN(kg),InflowTN(kg),AtmosphericDepositionN(kg),NH3GasLoss(kg),SED+SOD_NRelease(kg),NFluxtoSediments(kg),SedimentDiagenesisNH4Flux(kg),SedimentDiagenesisNO3Flux(kg)" ! 'JDAY,WB,TP-Waterbody(kg),TP-Sediment(kg),TP-Plants(kg),OutflowTP(kg),TributaryTP(kg),DistributedTributaryTP(kg),WithdrawalTP(kg),PrecipitationTP(kg),InflowTP(kg),BurialP_SED(kg),SED+SOD_PRelease(kg),TN-Waterbody(kg),TN-Sediment(kg),TN-Plants(kg),OutflowTN(kg),TributaryTN(kg),DistributedTributaryTN(kg),WithdrawalTN(kg),PrecipitationTN(kg),InflowTN(kg),BurialN_SED(kg),SED+SOD_NRelease(kg)'
                else
                    write(MASSBFN, "(A467)") "JDAY,WB,TP-Waterbody(kg),TP-Sediment(kg),TP-Plants(kg),OutflowTP(kg),TributaryTP(kg),DistributedTributaryTP(kg),WithdrawalTP(kg),PrecipitationTP(kg),InflowTP(kg),AtmosphericDepositionP(kg),SED+SOD_PRelease(kg),PFluxtoSediments(kg),TN-Waterbody(kg),TN-Sediment(kg),TN-Plants(kg),OutflowTN(kg),TributaryTN(kg),DistributedTributaryTN(kg),WithdrawalTN(kg),PrecipitationTN(kg),InflowTN(kg),AtmosphericDepositionN(kg),NH3GasLoss(kg),SED+SOD_NRelease(kg),NFluxtoSediments(kg)" ! 'JDAY,WB,TP-Waterbody(kg),TP-Sediment(kg),TP-Plants(kg),OutflowTP(kg),TributaryTP(kg),DistributedTributaryTP(kg),WithdrawalTP(kg),PrecipitationTP(kg),InflowTP(kg),BurialP_SED(kg),SED+SOD_PRelease(kg),TN-Waterbody(kg),TN-Sediment(kg),TN-Plants(kg),OutflowTN(kg),TributaryTN(kg),DistributedTributaryTN(kg),WithdrawalTN(kg),PrecipitationTN(kg),InflowTN(kg),BurialN_SED(kg),SED+SOD_NRelease(kg)'
                end if

                exit
            end if
        end do

        if (WLC == "      ON") then
            open(WLFN, FILE="wl.csv", STATUS="UNKNOWN")
            write(WLFN, '("JDAY,",*("SEG",i4,","))') ((i, i = us(jb), ds(jb)), jb = 1, nbr)
        end if

!**** DSI W2 Linkage File (W2L) (Supercedes Old Velocity vectors)
        if (VECTOR(1)) then
! *** Apply the same linkage settings for all waterbodies
            open(VPL(1), FILE=VPLFN(1), STATUS="UNKNOWN", ACCESS="SEQUENTIAL", FORM="BINARY")

! *** W2 Version
            write(VPL(1)) W2VER

            write(VPL(1)) TITLE

            write(VPL(1)) INT4(NWB), INT4(NBR), INT4(IMX), INT4(KMX), INT4(NCT), INT4(NAC)

! *** Flag the output file if using Outlet time series
            if (TIME_SERIES) then
                IFLAG = -1
            else
                IFLAG = 0
            end if
            write(VPL(1)) IFLAG

! *** MODEL CONFIGURATION
            do JW = 1, NWB
                write(VPL(1)) INT4(BS(JW)), INT4(BE(JW))
            end do

            write(VPL(1)) (REAL(DLX(I), 4), I = 1, IMX)
            write(VPL(1)) (REAL(PHI0(I), 4), I = 1, IMX)
            write(VPL(1)) (INT4(US(K)), K = 1, NBR)
            write(VPL(1)) (INT4(DS(K)), K = 1, NBR)
            write(VPL(1)) (INT4(UHS(K)), K = 1, NBR)
            write(VPL(1)) (INT4(DHS(K)), K = 1, NBR)
            write(VPL(1)) ((REAL(H(K, JW), 4), K = 1, KMX), JW = 1, NWB)

! *** INITIALIZE ALL CELLS AS INACTIVE (SET FLAG = 0)
            ICOMP = 0

! *** Set up the array to handle cell types : Outflows
            do JB = 1, NBR
                if (NSTR(JB) > 0) then
                    do JS = 1, NSTR(JB)
                        do K = KMX - 1, 2, -1
                            if (ESTR(JS, JB) <= EL(K, DS(JB))) then
                                ICOMP(K, DS(JB)) = 4
                            end if
                        end do
                    end do
!    ELSE      ! SW NEEDS FIXED SINCE WE DO NOT HAVE NOUT AND KOUT DEFINED 9/28/13
!     DO JO = 1, NOUT(JB)
!        ICOMP(KOUT(JO,JB),DS(JB)) = 3
!      END DO
                end if
            end do

! *** Set up the array to handle cell types : Withdrawals
            do IW = 1, NWD
!DO K=KTW(IW),KBW(IW),-1
                do K = KMX, 2, -1
                    if (EWD(IW) <= EL(K, IWD(IW))) then
                        ICOMP(K - 1, IWD(IW)) = 2
                        exit
                    end if
                end do
            end do

! *** Set up the array to handle cell types : Tributaries
            if (NTR > 0) then
                do JT = 1, NTR
                    do K = 2, KB(ITR(JT))
                        if (ICOMP(K, ITR(JT)) == 0) then
                            ICOMP(K, ITR(JT)) = 5
                        end if
                    end do
                end do
            end if

! *** FLAG ACTIVE CELLS
! *** (IF BATHYMETRY(B)>0, FLAG CELL WITH VALUE OF 1, IF IT HAS
! *** NOT BEEN SET ABOVE TO A VALUE OF 2,3, OR 4 )
            do JB = 1, NBR
                IU = US(JB)
                ID = DS(JB)
                do I = IU, ID
                    do K = 2, KB(I)
                        if (B(K, I) > 0. .and. ICOMP(K, I) == 0) then
                            ICOMP(K, I) = 1
                        end if
                    end do
                end do
            end do

! *** NOW WRITE SEGMENT WIDTHS (B)
            do I = 1, IMX
                do K = 1, KMX
                    if (ICOMP(K, I) > 0) then
                        write(VPL(1)) REAL(B(K, I), 4)
                    else
                        write(VPL(1)) 0.0_4
                    end if
                end do
            end do

! *** Output X -
! *** THIS SEQUENCE SWEEPS ACROSS COLUMNS WITH BLOCK CENTERED
! *** X COORDS WORKING UPSTREAM FOR EACH BRANCH
            do K = 1, KMX
                do JB = 1, NBR
                    IU = US(JB) - 1
                    ID = DS(JB) + 1
                    do I = ID, IU, -1
                        if (I == ID) then
                            DLX_OLD = 0.5*DLX(I)
                        else
                            DLX_OLD = DLX_OLD + 0.5*DLX(I) + 0.5*DLX(I + 1)
                        end if
                        write(VPL(1)) REAL(DLX_OLD, 4)
                    end do
                end do
            end do

! *** Output Y
! *** THIS SEQUENCE SWEEPS ACROSS COLUMNS WITH BLOCK CENTERED
! *** ELEVATIONS WORKING FROM THE BOTTOM LAYER TO THE TOP
            do K = 1, KMX
                do JB = 1, NBR
! *** BOTTOM INACTIVE CELL
!ELC = ELBOT(JW)-0.5*H(KMX,JW)
!WRITE(VPL(1)) ELC

                    IU = US(JB) - 1
                    ID = DS(JB) + 1
                    do I = ID, IU, -1
                        write(VPL(1)) REAL(EL(K, I), 4)
                    end do
                end do
            end do

            write(VPL(1)) ((ICOMP(K, I), K = 1, KMX), I = 1, IMX)

! *** CONSTITUENT NUMBER LIST
            write(VPL(1)) (INT4(CN(JC)), JC = 1, NAC)

! *** CONSTITUENT NAME
            write(VPL(1)) (CNAME1(CN(JC)), JC = 1, NAC)

! *** CONSTITUENT UNITS
            write(VPL(1)) (CUNIT2(CN(JC)), JC = 1, NAC)

            if (TIME_SERIES) then
!OPEN(ITSRNAG,FILE='ITSRNAG.TMP',STATUS='SCRATCH', ACCESS='SEQUENTIAL',FORM='BINARY')    PMC - FIX LATER
            end if

        end if

    end if

    return
end subroutine OUTPUTINIT
