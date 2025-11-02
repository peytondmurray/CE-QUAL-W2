subroutine OUTPUT(JDAY, IUPR, IDPR, KBR, ISNP, BL, NBL)
    use GLOBAL;     use GDAYC;     use GEOMC;     use KINETIC;     use TVDC;     use NAMESC;     use LOGICC
    use MACROPHYTEC
    use CEMAVars
    implicit none

! Type declaration

    real :: JDAY, LIMIT
    integer :: NBL, IUPR, IDPR, JH, L, K, J, NLINES, KBR, JAC, JD, JE, JT, JJ, JA
    integer, dimension(IMX) :: BL
    integer, dimension(IMX,NWB) :: ISNP
    logical :: NEW_PAGE
    character(len=8) :: LFAC
    character(len=10) :: BLANK = "          "

! Variable initialization

    NEW_PAGE = .true.

! Blank inactive cells

    NBL = 1
    JB = 1
    IUPR = 1
    do I = 1, IDPR - 1
        if (CUS(JB) > ISNP(I, JW)) then
            BL(NBL) = I
            NBL = NBL + 1
            if (JB == 1) then
                IUPR = I + 1
            end if
        end if
        if (ISNP(I + 1, JW) > DS(JB)) then
            JB = JB + 1
        end if
    end do
    NBL = NBL - 1

! Water surface elevation, water surface deviation, ice cover, and sediment oxygen demand

    CONV(1, :) = BLANK
    do I = IUPR, IDPR
        do JJB = 1, NBR
            if (ISNP(I, JW) >= US(JJB) - 1 .and. ISNP(I, JW) <= DS(JJB) + 1) then
                exit
            end if
        end do
        write(CONV(1, I), "(F10.3)") EL(KTWB(JW), ISNP(I, JW)) - Z(ISNP(I, JW))*COSA(JJB)
    end do
    write(SNP(JW), "(/A//2X,1000I10)") "          Water Surface, m", (ISNP(I, JW), I = IUPR, IDPR)
    write(SNP(JW), "(2X,1000A10/)") (CONV(1, I), I = IUPR, IDPR)
    do I = IUPR, IDPR
        write(CONV(1, I), "(F10.4)") SNGL(Z(ISNP(I, JW)))
    end do
    write(SNP(JW), "(/A//2X,1000I10)") "          Water Surface Deviation (positive downwards), m", (ISNP(I, JW), I = IUPR, IDPR)
    write(SNP(JW), "(2X,1000A10/)") (CONV(1, I), I = IUPR, IDPR)
    if (ICE_CALC(JW)) then
        do I = IUPR, IDPR
            write(CONV(1, I), "(F10.3)") ICETH(ISNP(I, JW))
        end do
        write(SNP(JW), "(/A//3X,1000A10)") "          Ice Thickness, m", (CONV(1, I), I = IUPR, IDPR)
    end if
    if (CONSTITUENTS) then
        do I = IUPR, IDPR
            write(CONV(1, I), "(F10.3)") SOD(ISNP(I, JW))*DAY
        end do
        if (OXYGEN_DEMAND) then
            write(SNP(JW), "(/A//3X,1000A10/)") "          Sediment Oxygen Demand, g/m^2/day", (CONV(1, I), I = IUPR, IDPR)
        end if
    end if

! Hydrodynamic variables and temperature

    CONV = BLANK
    do JH = 1, NHY
        L = LEN_TRIM(FMTH(JH))
        if (PRINT_HYDRO(JH, JW)) then
            do I = IUPR, IDPR
                if (JH == 1) then
                    do K = KTWB(JW), KB(ISNP(I, JW))
                        write(CONV(K, I), FMTH(JH)) INT(HYD(K, ISNP(I, JW), JH))
                    end do
                else
                    if (JH > 6) then
                        do K = KTWB(JW), KB(ISNP(I, JW))
                            write(CONV(K, I), FMTH(JH)(1:L)) HYD(K, ISNP(I, JW), JH)*DLT
                        end do
                    else
                        do K = KTWB(JW), KB(ISNP(I, JW))
                            write(CONV(K, I), FMTH(JH)(1:L)) HYD(K, ISNP(I, JW), JH)*HMULT(JH)
                        end do
                    end if
                end if
            end do
            if (NEW_PAGE) then
                write(SNP(JW), '("1",11(A/1X))') (TITLE(J), J = 1, 11)
                NLINES = KMX - KTWB(JW) + 14
            end if
            NLINES = NLINES + KMX - KTWB(JW) + 11
            NEW_PAGE = NLINES > 72
            write(SNP(JW), "(/1X,3(A,1X,I0),A,F0.2,A/)") MONTH, GDAY, ",", YEAR, "    Julian day = ", INT(JDAY), " days ", (JDAY - INT(JDAY))*24.0, " hours   " // HNAME(JH)
            write(SNP(JW), "(1X,A,1000I10)") "Layer  Depth", (ISNP(I, JW), I = IUPR, IDPR)
            do K = KTWB(JW), KBR
                write(SNP(JW), "(1X,I4,F8.2,1000A10)") K, DEPTHM(K, DS(BS(JW))), (CONV(K, I), I = IUPR, IDPR)
            end do
        end if
    end do

! Constituent concentrations

    if (CONSTITUENTS) then
        do JAC = 1, NAC
            JC = CN(JAC)
            L = LEN_TRIM(FMTC(JC))
            if (PRINT_CONST(JC, JW)) then
                do I = IUPR, IDPR
                    do K = KTWB(JW), KB(ISNP(I, JW))
                        write(CONV(K, I), FMTC(JC)(1:L)) C2(K, ISNP(I, JW), JC)*CMULT(JC)
                    end do
                end do
                if (NEW_PAGE) then
                    write(SNP(JW), '("1",11(A/1X))') (TITLE(J), J = 1, 11)
                    NLINES = KMX - KTWB(JW) + 14
                end if
                NLINES = NLINES + KMX - KTWB(JW) + 11
                NEW_PAGE = NLINES > 72
                write(SNP(JW), "(/1X,3(A,1X,I0),A,F0.2,A/)") MONTH, GDAY, ",", YEAR, "    Julian Date ", INT(JDAY), " days ", (JDAY - INT(JDAY))*24.0, " hours   " // CNAME(JC)
                write(SNP(JW), "(1X,A,1000I10)") "Layer  Depth", (ISNP(I, JW), I = IUPR, IDPR)
                do K = KTWB(JW), KBR
                    write(SNP(JW), "(1X,I4,F8.2,1000A10)") K, DEPTHM(K, DS(BS(JW))), (CONV(K, I), I = IUPR, IDPR)
                end do
            end if
        end do

!** Derived constituent concentrations

        do JD = 1, NDC
            L = LEN_TRIM(FMTCD(JD))
            if (PRINT_DERIVED(JD, JW)) then
                do I = IUPR, IDPR
                    do K = KTWB(JW), KB(ISNP(I, JW))
                        write(CONV(K, I), FMTCD(JD)(1:L)) CD(K, ISNP(I, JW), JD)*CDMULT(JD)
                    end do
                end do
                if (NEW_PAGE) then
                    write(SNP(JW), '("1",11(A/1X))') (TITLE(J), J = 1, 11)
                    NLINES = KMX - KTWB(JW) + 14
                end if
                NLINES = NLINES + KMX - KTWB(JW) + 11
                NEW_PAGE = NLINES > 72
                write(SNP(JW), "(/1X,3(A,1X,I0),A,F0.2,A/)") MONTH, GDAY, ",", YEAR, "    Julian Date ", INT(JDAY), " days ", (JDAY - INT(JDAY))*24.0, " hours    " // CDNAME(JD)
                write(SNP(JW), "(1X,A,1000I10)") "Layer  Depth", (ISNP(I, JW), I = IUPR, IDPR)
                do K = KTWB(JW), KBR
                    write(SNP(JW), "(1X,I4,F8.2,1000A10)") K, DEPTHM(K, DS(BS(JW))), (CONV(K, I), I = IUPR, IDPR)
                end do
            end if
        end do

!** Sediment

        if (PRINT_SEDIMENT(JW)) then
            do I = IUPR, IDPR
                do K = KTWB(JW), KB(ISNP(I, JW))
                    write(CONV(K, I), "(F10.2)") SED(K, ISNP(I, JW))
                end do
            end do
            if (NEW_PAGE) then
                write(SNP(JW), '("1",11(A/1X))') (TITLE(J), J = 1, 11)
                NLINES = KMX - KTWB(JW) + 14
            end if
            NLINES = NLINES + KMX - KTWB(JW) + 11
            NEW_PAGE = NLINES > 72
            write(SNP(JW), "(/1X,3(A,1X,I0),A,F0.2,A/)") MONTH, GDAY, ",", YEAR, "    Julian Date ", INT(JDAY), " days ", (JDAY - INT(JDAY))*24.0, " hours     Organic sediments, g/m^3"
            write(SNP(JW), "(1X,A,1000I10)") "Layer  Depth", (ISNP(I, JW), I = IUPR, IDPR)
            do K = KTWB(JW), KBR
                write(SNP(JW), "(1X,I4,F8.2,1000A10)") K, DEPTHM(K, DS(BS(JW))), (CONV(K, I), I = IUPR, IDPR)
            end do
        end if

        if (PRINT_SEDIMENT(JW)) then
            do I = IUPR, IDPR
                do K = KTWB(JW), KB(ISNP(I, JW))
                    write(CONV(K, I), "(F10.2)") SEDp(K, ISNP(I, JW))
                end do
            end do
            if (NEW_PAGE) then
                write(SNP(JW), '("1",11(A/1X))') (TITLE(J), J = 1, 11)
                NLINES = KMX - KTWB(JW) + 14
            end if
            NLINES = NLINES + KMX - KTWB(JW) + 11
            NEW_PAGE = NLINES > 72
            write(SNP(JW), "(/1X,3(A,1X,I0),A,F0.2,A/)") MONTH, GDAY, ",", YEAR, "    Julian Date ", INT(JDAY), " days ", (JDAY - INT(JDAY))*24.0, " hours     Organic phosphorus sediments, g/m^3"
            write(SNP(JW), "(1X,A,1000I10)") "Layer  Depth", (ISNP(I, JW), I = IUPR, IDPR)
            do K = KTWB(JW), KBR
                write(SNP(JW), "(1X,I4,F8.2,1000A10)") K, DEPTHM(K, DS(BS(JW))), (CONV(K, I), I = IUPR, IDPR)
            end do
        end if

        if (PRINT_SEDIMENT(JW)) then
            do I = IUPR, IDPR
                do K = KTWB(JW), KB(ISNP(I, JW))
                    write(CONV(K, I), "(F10.2)") SEDn(K, ISNP(I, JW))
                end do
            end do
            if (NEW_PAGE) then
                write(SNP(JW), '("1",11(A/1X))') (TITLE(J), J = 1, 11)
                NLINES = KMX - KTWB(JW) + 14
            end if
            NLINES = NLINES + KMX - KTWB(JW) + 11
            NEW_PAGE = NLINES > 72
            write(SNP(JW), "(/1X,3(A,1X,I0),A,F0.2,A/)") MONTH, GDAY, ",", YEAR, "    Julian Date ", INT(JDAY), " days ", (JDAY - INT(JDAY))*24.0, " hours     Organic nitrogen sediments, g/m^3"
            write(SNP(JW), "(1X,A,1000I10)") "Layer  Depth", (ISNP(I, JW), I = IUPR, IDPR)
            do K = KTWB(JW), KBR
                write(SNP(JW), "(1X,I4,F8.2,1000A10)") K, DEPTHM(K, DS(BS(JW))), (CONV(K, I), I = IUPR, IDPR)
            end do
        end if

        if (PRINT_SEDIMENT(JW)) then
            do I = IUPR, IDPR
                do K = KTWB(JW), KB(ISNP(I, JW))
                    write(CONV(K, I), "(F10.2)") SEDc(K, ISNP(I, JW))
                end do
            end do
            if (NEW_PAGE) then
                write(SNP(JW), '("1",11(A/1X))') (TITLE(J), J = 1, 11)
                NLINES = KMX - KTWB(JW) + 14
            end if
            NLINES = NLINES + KMX - KTWB(JW) + 11
            NEW_PAGE = NLINES > 72
            write(SNP(JW), "(/1X,3(A,1X,I0),A,F0.2,A/)") MONTH, GDAY, ",", YEAR, "    Julian Date ", INT(JDAY), " days ", (JDAY - INT(JDAY))*24.0, " hours     Organic carbon sediments, g/m^3"
            write(SNP(JW), "(1X,A,1000I10)") "Layer  Depth", (ISNP(I, JW), I = IUPR, IDPR)
            do K = KTWB(JW), KBR
                write(SNP(JW), "(1X,I4,F8.2,1000A10)") K, DEPTHM(K, DS(BS(JW))), (CONV(K, I), I = IUPR, IDPR)
            end do
        end if

! Amaila Start
        if (PRINT_SEDIMENT1(JW)) then
            do I = IUPR, IDPR
                do K = KTWB(JW), KB(ISNP(I, JW))
                    write(CONV(K, I), "(F10.2)") SED1(K, ISNP(I, JW))
                end do
            end do
            if (NEW_PAGE) then
                write(SNP(JW), '("1",11(A/1X))') (TITLE(J), J = 1, 11)
                NLINES = KMX - KTWB(JW) + 14
            end if
            NLINES = NLINES + KMX - KTWB(JW) + 11
            NEW_PAGE = NLINES > 72
            write(SNP(JW), "(/1X,3(A,1X,I0),A,F0.2,A/)") MONTH, GDAY, ",", YEAR, "    Julian Date ", INT(JDAY), " days ", (JDAY - INT(JDAY))*24.0, " hours     Tree Organic Matter-Labile, g/m^3"
            write(SNP(JW), "(1X,A,1000I10)") "Layer  Depth", (ISNP(I, JW), I = IUPR, IDPR)
            do K = KTWB(JW), KBR
                write(SNP(JW), "(1X,I4,F8.2,1000A10)") K, DEPTHM(K, DS(BS(JW))), (CONV(K, I), I = IUPR, IDPR)
            end do
        end if

        if (PRINT_SEDIMENT2(JW)) then
            do I = IUPR, IDPR
                do K = KTWB(JW), KB(ISNP(I, JW))
                    write(CONV(K, I), "(F10.2)") SED2(K, ISNP(I, JW))
                end do
            end do
            if (NEW_PAGE) then
                write(SNP(JW), '("1",11(A/1X))') (TITLE(J), J = 1, 11)
                NLINES = KMX - KTWB(JW) + 14
            end if
            NLINES = NLINES + KMX - KTWB(JW) + 11
            NEW_PAGE = NLINES > 72
            write(SNP(JW), "(/1X,3(A,1X,I0),A,F0.2,A/)") MONTH, GDAY, ",", YEAR, "    Julian Date ", INT(JDAY), " days ", (JDAY - INT(JDAY))*24.0, " hours     Tree Organic Matter-Refactory, g/m^3"
            write(SNP(JW), "(1X,A,1000I10)") "Layer  Depth", (ISNP(I, JW), I = IUPR, IDPR)
            do K = KTWB(JW), KBR
                write(SNP(JW), "(1X,I4,F8.2,1000A10)") K, DEPTHM(K, DS(BS(JW))), (CONV(K, I), I = IUPR, IDPR)
            end do
        end if

! Amaila End

!** Epiphyton

        do JE = 1, NEP
            if (PRINT_EPIPHYTON(JW, JE)) then
                do I = IUPR, IDPR
                    do K = KTWB(JW), KB(ISNP(I, JW))
                        write(CONV(K, I), "(F10.2)") EPD(K, ISNP(I, JW), JE)
                    end do
                end do
                if (NEW_PAGE) then
                    write(SNP(JW), '("1",11(A/1X))') (TITLE(J), J = 1, 11)
                    NLINES = KMX - KTWB(JW) + 14
                end if
                NLINES = NLINES + KMX - KTWB(JW) + 11
                NEW_PAGE = NLINES > 72
                write(SNP(JW), "(/1X,3(A,1X,I0),A,F0.2,A/)") MONTH, GDAY, ",", YEAR, "    Julian Date ", INT(JDAY), " days ", (JDAY - INT(JDAY))*24.0, " hours     Epiphyton, g/m^2"
                write(SNP(JW), "(1X,A,1000I10)") "Layer  Depth", (ISNP(I, JW), I = IUPR, IDPR)
                do K = KTWB(JW), KBR
                    write(SNP(JW), "(1X,I4,F8.2,1000A10)") K, DEPTHM(K, DS(BS(JW))), (ADJUSTR(CONV(K, I)), I = IUPR, IDPR)
                end do
            end if
        end do

!********* macrophytes
        do L = 1, nmc
            CONV = BLANK
            if (PRINT_macrophyte(jw, L)) then
                do I = IUPR, IDPR
                    do K = ktwb(jw), KB(ISNP(I, jw))
                        write(CONV(K, I), "(F10.2)") mac(K, ISNP(I, jw), L)
                    end do
                end do
                if (NEW_PAGE) then
                    write(SNP(JW), '("1",11(A/1X))') (TITLE(J), J = 1, 11)
                    NLINES = KMX - KTWB(JW) + 14
                end if
                NLINES = NLINES + KMX - KTWB(JW) + 11
                NEW_PAGE = NLINES > 72
                write(SNP(jw), "(/1X,3(A,1X,I0),A,F0.2,2A,i0,A/)") MONTH, GDAY, ",", YEAR, "    Julian Date ", INT(JDAY), " days ", (JDAY - INT(JDAY))*24.0, " hours   ", "  Macrophyte Group #", L, " g/m^3"
                write(SNP(JW), "(1X,A,1000I10)") "Layer  Depth", (ISNP(I, JW), I = IUPR, IDPR)
                do K = ktwb(jw), KBR
                    write(SNP(JW), "(1X,I4,F8.2,1000A10)") K, DEPTHM(K, DS(BS(JW))), (CONV(K, I), I = IUPR, IDPR)
                end do
                do i = iupr, idpr
                    conv2 = blank
                    do K = ktwb(jw), KB(ISNP(I, jw))
                        if (k == ktwb(jw)) then
                            jt = kti(isnp(i, jw))
                        else
                            jt = k
                        end if
                        je = kb(isnp(i, jw))
                        do jj = jt, je
                            write(CONV2(K, jj), "(F10.2)") macrc(jj, K, ISNP(I, jw), L)
                        end do
                    end do
                    if (NEW_PAGE) then
                        write(SNP(JW), '("1",11(A/1X))') (TITLE(J), J = 1, 11)
                        NLINES = KMX - KTWB(JW) + 14
                    end if
                    NLINES = NLINES + KMX - KTWB(JW) + 11
                    NEW_PAGE = NLINES > 72
!          WRITE (SNP(jw),'(/1X,3(A,1X,I0),A,F0.2,A/)')  MONTH,GDAY,',',YEAR,'    Julian Date ',INT(JDAY),' days ',(JDAY-INT(JDAY))    &
!                                                     *24.0,' hours   ','  Macrophyte Columns g/m^3'
                    write(SNP(JW), "(/1X,3(A,1X,I0),A,F0.2,2A,i0,A/)") MONTH, GDAY, ",", YEAR, "    Julian Date ", INT(JDAY), " days ", (JDAY - INT(JDAY))*24.0, " hours   ", "  Macrophyte Group #", L, " Columns g/m^3"
                    write(SNP(JW), 3052) ISNP(I, JW)
                    3052 format(7X,"SEGMENT   ",I8)
                    JT = KTI(ISNP(I, JW))
                    JE = KB(ISNP(I, JW))
                    write(SNP(JW), "(A,200I10)") " LAYER  DEPTH", (JJ, JJ = JT, JE)
                    do K = KTWB(JW), KB(ISNP(I, JW))
                        write(SNP(JW), "(1X,I4,F8.2,1000A10)") K, DEPTHM(K, DS(BS(JW))), (CONV2(K, JJ), JJ = JT, JE)
                    end do
                    write(SNP(JW), 3053) L, ISNP(I, JW)
                    3053 format(7X,"MACROPHYTE GROUP #",I0," LIMITATION SEGMENT   ",I8) ! CB 12/20/11
                    write(SNP(JW), "(A,200I10)") " LAYER  DEPTH", (JJ, JJ = JT, JE)
                    MLFPR = BLANK
                    do K = KTWB(JW), KB(ISNP(I, JW))
                        if (K == KTWB(JW)) then
                            JT = KTI(ISNP(I, JW))
                        else
                            JT = K
                        end if
                        do JJ = JT, JE
                            LIMIT = MIN(MPLIM(K, ISNP(I, JW), L), MNLIM(K, ISNP(I, JW), L), MCLIM(K, ISNP(I, JW), L), MLLIM(JJ, K, ISNP(I, JW), L))
                            if (LIMIT == MPLIM(K, ISNP(I, JW), L)) then
                                write(LFAC, "(F8.4)") MPLIM(K, ISNP(I, JW), L)
                                MLFPR(JJ, K, ISNP(I, JW), L) = " P" // LFAC
                            else
                                if (LIMIT == MNLIM(K, I, L)) then
                                    write(LFAC, "(F8.4)") MNLIM(K, ISNP(I, JW), L)
                                    MLFPR(JJ, K, ISNP(I, JW), L) = " N" // LFAC
                                else
                                    if (LIMIT == MCLIM(K, I, L)) then
                                        write(LFAC, "(F8.4)") MCLIM(K, ISNP(I, JW), L)
                                        MLFPR(JJ, K, ISNP(I, JW), L) = " C" // LFAC
                                    else
                                        if (LIMIT == MLLIM(JJ, K, ISNP(I, JW), L)) then
                                            write(LFAC, "(F8.4)") MLLIM(JJ, K, ISNP(I, JW), L)
                                            MLFPR(JJ, K, ISNP(I, JW), L) = " L" // LFAC
                                        end if
                                    end if
                                end if
                            end if
                        end do
                    end do
                    JT = KTI(ISNP(I, JW))
                    JE = KB(ISNP(I, JW))
                    do K = KTWB(JW), KB(ISNP(I, JW))
                        write(SNP(JW), "(1X,I4,F8.2,1000A10)") K, DEPTHM(K, DS(BS(JW))), (MLFPR(JJ, K, ISNP(I, JW), L), JJ = JT, JE)
                    end do
                end do
            end if
        end do

!** Algal nutrient limitations

        do JA = 1, NAL
            if (LIMITING_FACTOR(JA)) then
                if (NEW_PAGE) then
                    write(SNP(JW), '("1",11(A/1X))') (TITLE(J), J = 1, 11)
                    NLINES = KMX - KTWB(JW) + 14
                end if
                NLINES = NLINES + KMX - KTWB(JW) + 11
                NEW_PAGE = NLINES > 72
                do I = IUPR, IDPR !mlm 6/30/06
                    do K = KTWB(JW), KB(ISNP(I, JW)) !mlm  6/30/06
                        LIMIT = MIN(APLIM(K, ISNP(I, JW), JA), ANLIM(K, ISNP(I, JW), JA), ASLIM(K, ISNP(I, JW), JA), ALLIM(K, ISNP(I, JW), JA))
                        if (LIMIT == APLIM(K, ISNP(I, JW), JA)) then
                            write(LFAC, "(F8.4)") APLIM(K, ISNP(I, JW), JA)
                            LFPR(K, ISNP(I, JW)) = " P" // LFAC
                        else
                            if (LIMIT == ANLIM(K, ISNP(I, JW), JA)) then
                                write(LFAC, "(F8.4)") ANLIM(K, ISNP(I, JW), JA)
                                LFPR(K, ISNP(I, JW)) = " N" // LFAC
                            else
                                if (LIMIT == ASLIM(K, ISNP(I, JW), JA)) then
                                    write(LFAC, "(F8.4)") ASLIM(K, ISNP(I, JW), JA)
                                    LFPR(K, ISNP(I, JW)) = " S" // LFAC
                                else
                                    if (LIMIT == ALLIM(K, ISNP(I, JW), JA)) then
                                        write(LFAC, "(F8.4)") ALLIM(K, ISNP(I, JW), JA)
                                        LFPR(K, ISNP(I, JW)) = " L" // LFAC
                                    end if
                                end if
                            end if
                        end if
                    end do
                end do
                write(SNP(JW), "(/1X,3(A,1X,I0),A,F0.2,A,I0,A/)") MONTH, GDAY, ",", YEAR, "    Julian Date", INT(JDAY), " days ", (JDAY - INT(JDAY))*24.0, " hours    Algal group ", JA, " limiting factor"
                write(SNP(JW), "(1X,A,1000I10)") "Layer  Depth", (ISNP(I, JW), I = IUPR, IDPR)
                do K = KTWB(JW), KBR
                    write(SNP(JW), "(1X,I4,F8.2,1000A)") K, DEPTHM(K, DS(BS(JW))), (LFPR(K, ISNP(I, JW)), I = IUPR, IDPR)
                end do
            end if
        end do

!** Epiphyton nutrient limitations

        do JE = 1, NEP
            if (PRINT_EPIPHYTON(JW, JE)) then
                if (NEW_PAGE) then
                    write(SNP(JW), '("1",11(A/1X))') (TITLE(J), J = 1, 11)
                    NLINES = KMX - KTWB(JW) + 14
                end if
                NLINES = NLINES + KMX - KTWB(JW) + 11
                NEW_PAGE = NLINES > 72
                do I = IUPR, IDPR !mlm   6/30/2006
                    do K = KTWB(JW), KB(ISNP(I, JW)) !mlm  6/30/2006
                        LIMIT = MIN(EPLIM(K, ISNP(I, JW), JE), ENLIM(K, ISNP(I, JW), JE), ESLIM(K, ISNP(I, JW), JE), ELLIM(K, ISNP(I, JW), JE))
                        if (LIMIT == EPLIM(K, ISNP(I, JW), JE)) then
                            write(LFAC, "(F8.4)") EPLIM(K, ISNP(I, JW), JE)
                            LFPR(K, ISNP(I, JW)) = " P" // LFAC
                        else
                            if (LIMIT == ENLIM(K, ISNP(I, JW), JE)) then
                                write(LFAC, "(F8.4)") ENLIM(K, ISNP(I, JW), JE)
                                LFPR(K, ISNP(I, JW)) = " N" // LFAC
                            else
                                if (LIMIT == ESLIM(K, ISNP(I, JW), JE)) then
                                    write(LFAC, "(F8.4)") ESLIM(K, ISNP(I, JW), JE)
                                    LFPR(K, ISNP(I, JW)) = " S" // LFAC
                                else
                                    if (LIMIT == ELLIM(K, ISNP(I, JW), JE)) then
                                        write(LFAC, "(F8.4)") ELLIM(K, ISNP(I, JW), JE)
                                        LFPR(K, ISNP(I, JW)) = " L" // LFAC
                                    end if
                                end if
                            end if
                        end if
                    end do
                end do
                write(SNP(JW), "(/1X,3(A,1X,I0),A,F0.2,A,I0,A/)") MONTH, GDAY, ",", YEAR, "    Julian Date", INT(JDAY), " days ", (JDAY - INT(JDAY))*24.0, " hours    Epiphyton group ", JE, " limiting factor"
                write(SNP(JW), "(1X,A,1000I10)") "Layer  Depth", (ISNP(I, JW), I = IUPR, IDPR)
                do K = KTWB(JW), KBR
                    write(SNP(JW), "(1X,I4,F8.2,1000A10)") K, DEPTHM(K, DS(BS(JW))), (LFPR(K, ISNP(I, JW)), I = IUPR, IDPR)
                end do
            end if
        end do
    end if
end subroutine OUTPUT
