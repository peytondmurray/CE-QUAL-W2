subroutine OUTPUTA()

    use MAIN
    use GLOBAL;     use NAMESC;     use GEOMC;     use LOGICC;     use PREC;     use SURFHE;     use KINETIC;     use SHADEC;     use EDDY
    use STRUCTURES;     use TRANS;     use TVDC;     use SELWC;     use GDAYC;     use SCREENC;     use TDGAS;     use RSTART
    use MACROPHYTEC;     use POROSITYC;     use ZOOPLANKTONC;     use BIOENERGETICS
    use CEMAVars
    use CEMAOutputRoutines;     use ALGAE_TOXINS
    implicit none

    external :: RESTART_OUTPUT

    integer :: JAD, JAF, IFLAG, JWWD, NLINES, ITOT, JJ, NUMOUTLETS, JSSS(100), JJC, IC, KK
    real :: TVOLAVG, QSUMM, CGASD, XDUM, QOUTLET(100), TOUTLET(100), VOLTOT, TBLANK ! SW 2/28/2020
    real(R8) :: DLVBR, DLE, CGAS, TGATE, TSPILL

! *** DSI W2_TOOL LINKAGE
    real(4), save, allocatable, dimension(:) :: WSEL
    real(4), save, allocatable, dimension(:,:) :: WDSI

    if (VECTOR(1)) then
        if (.not. ALLOCATED(WSEL)) then
            allocate(WSEL(IMX))
            allocate(WDSI(KMX, IMX))
        end if
    end if


!***********************************************************************************************************************************
!*                                                    Task 2.8: Output Results                                                    **
!***********************************************************************************************************************************

    if (LAKE_RIVER_CONTOUR_ON == "ON") then ! SW 2/28/2020
        do JJ = 1, NUM_LAKE_CONTOUR
            if (JDAY >= NXT_LAKE_CONTOUR(JJ)) then
                NXT_LAKE_CONTOUR(JJ) = NXT_LAKE_CONTOUR(JJ) + LAKE_CONTOUR_FREQ(JJ)
                TBLANK = -99.0
                if (LAKE_CONTOUR_FORMAT == 1) then
                    write(LAKE_RIVER_CONTOUR + JJ - 1, '(F10.4,",",*(F8.2,","))') JDAY, ELWS(LAKE_CONTOUR_SEG(JJ)) + 0.1, TBLANK
                    write(LAKE_RIVER_CONTOUR + JJ - 1, '(F10.4,",",*(F8.2,","))') JDAY, ELWS(LAKE_CONTOUR_SEG(JJ)), T2(KTWB(JW_LAKE_CONTOUR(JJ)), LAKE_CONTOUR_SEG(JJ))
                    do K = KTWB(JW_LAKE_CONTOUR(JJ)) + 1, KB(LAKE_CONTOUR_SEG(JJ))
                        write(LAKE_RIVER_CONTOUR + JJ - 1, '(F10.4,",",*(F8.2,","))') JDAY, (EL(K, LAKE_CONTOUR_SEG(JJ)) + EL(K, LAKE_CONTOUR_SEG(JJ)))*0.5, T2(K, LAKE_CONTOUR_SEG(JJ))
                    end do
                else
                    write(LAKE_RIVER_CONTOUR + JJ - 1, '(F10.4,",",*(F8.2,","))') JDAY, (TBLANK, K = 2, KTWB(JW_LAKE_CONTOUR(JJ)) - 1), (T2(K, LAKE_CONTOUR_SEG(JJ)), K = KTWB(JW_LAKE_CONTOUR(JJ)), KB(LAKE_CONTOUR_SEG(JJ)))
                end if
                if (OXYGEN_DEMAND) then
                    if (LAKE_CONTOUR_FORMAT == 1) then
                        write(LAKE_RIVER_CONTOUR + 10 + JJ - 1, '(F10.4,",",*(F8.2,","))') JDAY, ELWS(LAKE_CONTOUR_SEG(JJ)) + 0.1, TBLANK
                        write(LAKE_RIVER_CONTOUR + 10 + JJ - 1, '(F10.4,",",*(F8.2,","))') JDAY, ELWS(LAKE_CONTOUR_SEG(JJ)), C2(KTWB(JW_LAKE_CONTOUR(JJ)), LAKE_CONTOUR_SEG(JJ), NDO)
                        do K = KTWB(JW_LAKE_CONTOUR(JJ)) + 1, KB(LAKE_CONTOUR_SEG(JJ))
                            write(LAKE_RIVER_CONTOUR + 10 + JJ - 1, '(F10.4,",",*(F8.2,","))') JDAY, (EL(K, LAKE_CONTOUR_SEG(JJ)) + EL(K, LAKE_CONTOUR_SEG(JJ)))*0.5, C2(K, LAKE_CONTOUR_SEG(JJ), NDO)
                        end do
                    else
                        write(LAKE_RIVER_CONTOUR + 10 + JJ - 1, '(F10.4,",",*(F8.2,","))') JDAY, (TBLANK, K = 2, KTWB(JW_LAKE_CONTOUR(JJ)) - 1), (C2(K, LAKE_CONTOUR_SEG(JJ), NDO), K = KTWB(JW_LAKE_CONTOUR(JJ)), KB(LAKE_CONTOUR_SEG(JJ)))
                    end if
                end if
            end if
        end do
        do JJ = 1, NUM_RIVER_CONTOUR
            if (JDAY >= NXT_RIVER_CONTOUR(JJ)) then
                NXT_RIVER_CONTOUR(JJ) = NXT_RIVER_CONTOUR(JJ) + RIVER_CONTOUR_FREQ(JJ)
                if (RIVER_CONTOUR_FORMAT == 1) then
                    do JB = RIVER_CONTOUR_BR1(JJ), RIVER_CONTOUR_BR2(JJ)
                        do I = US(JB), DS(JB)
                            write(LAKE_RIVER_CONTOUR + 20 + JJ - 1, '(F10.4,",",*(F8.2,","))') JDAY, X1(I), T2(KTWB(JW_RIVER_CONTOUR(JJ)), I)
                        end do
                    end do
                else
                    write(LAKE_RIVER_CONTOUR + 20 + JJ - 1, '(F10.4,",",*(F8.2,","))') JDAY, ((T2(KTWB(JW_RIVER_CONTOUR(JJ)), I), I = US(JB), DS(JB)), JB = RIVER_CONTOUR_BR1(JJ), RIVER_CONTOUR_BR2(JJ))
                end if
                if (OXYGEN_DEMAND) then
                    if (RIVER_CONTOUR_FORMAT == 1) then
                        do JB = RIVER_CONTOUR_BR1(JJ), RIVER_CONTOUR_BR2(JJ)
                            do I = US(JB), DS(JB)
                                write(LAKE_RIVER_CONTOUR + 30 + JJ - 1, '(F10.4,",",*(F8.2,","))') JDAY, X1(I), C2(KTWB(JW_RIVER_CONTOUR(JJ)), I, NDO)
                            end do
                        end do
                    else
                        write(LAKE_RIVER_CONTOUR + 30 + JJ - 1, '(F10.4,",",*(F8.2,","))') JDAY, ((C2(KTWB(JW_RIVER_CONTOUR(JJ)), I, NDO), I = US(JB), DS(JB)), JB = RIVER_CONTOUR_BR1(JJ), RIVER_CONTOUR_BR2(JJ))
                    end if
                end if

            end if
        end do
    end if

! BIOEXP MLM
    if (BIOEXP) then
        if (JDAY >= NXTBIO) then
            NXTBIO = NXTBIO + BIOF(1) ! mlm3
            do J = 1, NIBIO

                I = IBIO(J)
                do JW = 1, NWB
                    if (I >= US(BS(JW)) - 1 .and. I <= DS(BE(JW)) + 1) then
                        exit
                    end if
                end do
! bioexp
                do K = KTWB(JW), KBI(I)
                    do JC = NZOOS, NZOOE
                        C2ZOO(K, I, JC) = C2(K, I, JC)*CMULT(JC)
                    end do
                    write(BIOEXPFN(J), '(F8.2,",",I8,",",3(F8.2,","),<NZOOE-NZOOS+1>(F8.3,","),I8,",",2(F8.2,","),A,",",I0,",",I0)') JDAY, IBIO(J), DEPTHM(K, I), T1(K, I), GAMMA(K, I), (C2ZOO(K, I, JC), JC = NZOOS, NZOOE), K, BH(K, I), EL(K, I), MONTH, GDAY, YEAR

                end do
! VOLUME WEIGHTING OF ACTIVE CONSTITUENTS     !MLM 18.07.06
                KLIM = MIN(KB(I), KTWB(JW) + 5)
                do K = KTWB(JW), KLIM
                    do JJC = 1, NAC
                        C2W(I, JJC) = C2W(I, JJC) + C2(K, I, CN(JJC))*CMULT(CN(JJC))*BH(K, I)
                    end do
                    C2W(I, NAC + 1) = C2W(I, NAC + 1) + CD(K, I, PH_DER)*BH(K, I) !PH
                    C2W(I, NAC + 2) = C2W(I, NAC + 2) + CD(K, I, TP_DER)*BH(K, I) ! TOTAL PHOSPHOROUS
                    VOLROOS(I) = VOLROOS(I) + BH(K, I) ! VOLUME WEIGHTED
                end do
                do JJC = 1, NAC + 2
                    C2W(I, JJC) = C2W(I, JJC)/VOLROOS(I)
                end do

                write(WEIGHTNUM(J), '(F10.3,",",<NAC>(F9.3,","),2(F9.3,","))') JDAY, (C2W(I, JJC), JJC = 1, NAC + 2)
                VOLROOS = 0.0

            end do
        end if
    end if

    if (WLC == "      ON") then
        if (JDAY >= NXWL) then
            NXWL = NXWL + WLF
! write out water level File
            write(WLFN, '(f10.3,",",*(f8.3,","))') jday, ((elws(i), i = us(jb), ds(jb)), jb = 1, nbr)
        end if
    end if

!** Time series

    if (TIME_SERIES) then
        if (JDAY >= NXTMTS .or. JDAY >= TSRD(TSRDP + 1)) then
            if (JDAY >= TSRD(TSRDP + 1)) then
                TSRDP = TSRDP + 1
                NXTMTS = TSRD(TSRDP)
            end if
            NXTMTS = NXTMTS + TSRF(TSRDP)

            do J = 1, NIKTSR
                I = ITSR(J)
! find out if segment is inactive OR cell is inactive for fixed layer - do not write out tsr file  ! SW 7/24/2018
                if (BR_INACTIVE(JBTSR(J))) then
                    cycle
                end if
                if (CUS(JBTSR(J)) > I) then
                    cycle
                end if
                do JW = 1, NWB
                    if (I >= US(BS(JW)) - 1 .and. I <= DS(BE(JW)) + 1) then
                        exit
                    end if
                end do
                if (ETSR(J) < 0) then ! SW 7/24/2018
                    if (INT(ABS(ETSR(J))) < KTWB(JW)) then
                        cycle
                    end if
                end if

! TEMP VOL WEIGHTED AVERAGE  SW 6/1/2015
                TVOLAVG = 0.0 ! SW 6/1/2015
                VOLTOT = 0.0
                do K = KTWB(JW), KB(I)
                    VOLTOT = VOLTOT + VOL(K, I)
                    TVOLAVG = TVOLAVG + T1(K, I)*VOL(K, I)
                end do
                if (KB(I) >= KTWB(JW)) then
                    tvolavg = tvolavg/voltot
                end if

                if (ETSR(J) < 0) then
                    K = INT(ABS(ETSR(J)))
                else
                    do K = KTWB(JW), KB(I)
                        if (DEPTHB(K, I) > ETSR(J)) then
                            exit
                        end if
                    end do
                    if (K > KB(I)) then
                        cycle
                    end if
                end if

                if (K >= KTWB(JW)) then ! SW 4/4/2018 ADDED TO ELIMINATE THE LAST VALUE OF A VARIABLE BEING USED FOR CASE WHEN ESTR IS NEGATIVE

                    do JAC = 1, NAC
                        L = LEN_TRIM(FMTC(CN(JAC)))
                        write(C2CH(JAC), FMTC(CN(JAC))(1:L)) C2(K, I, CN(JAC))*CMULT(CN(JAC))
                    end do
                    do JF = 1, NAF(JW)
!WRITE(KFCH(JF),'(E10.3)')KF(K,I,(KFCN(JF,JW)))*VOL(K,I)/1000./DAY
                        write(KFCH(JF), "(E10.3)") KF(K, I, KFCN(JF, JW))*VOL(K, I)/1000.*DAY ! BG 12/09/18
                    end do
                    do JAD = 1, NACD(JW)
                        L = LEN_TRIM(FMTCD(CDN(JAD, JW)))
                        write(CDCH(JAD), FMTCD(CDN(JAD, JW))(1:L)) CD(K, I, CDN(JAD, JW))*CDMULT(CDN(JAD, JW))
                    end do
                    do JE = 1, NEP
                        write(EPCH(JE), "(F10.3)") EPD(K, I, JE) ! SW 8/13/06
                    end do
                    do JA = 1, NAL
                        write(APCH(JA), "(F10.3)") APLIM(K, I, JA) ! SW 8/13/06
                        write(ANCH(JA), "(F10.3)") ANLIM(K, I, JA) ! SW 8/13/06
                        write(ALCH(JA), "(F10.3)") ALLIM(K, I, JA) ! SW 8/13/06
                    end do
                    do JM = 1, NMC
                        write(macCH(Jm), "(F10.3)") mac(K, I, Jm) ! SW 8/13/06
                    end do
                    if (SEDIMENT_CALC(JW)) then
                        write(sedch, "(F10.3)") sed(K, I) ! SW 8/13/06
                        write(sedpch, "(F10.3)") sedp(K, I)
                        write(sednch, "(F10.3)") sedn(K, I)
                        write(sedcch, "(F10.3)") sedc(K, I)
                    end if
                    if (ICE_COMPUTATION) then
                        if (SEDIMENT_CALC(JW)) then
                            write(TSR(J), '(f10.3,",",19(F10.3,","),*(A,","))') JDAY, DLT, ELWS(I), T1(K, I), U(K, I), QC(I), SRON(JW)*1.06, GAMMA(K, I), DEPTHB(KB(I), I), BI(KTWB(JW), I), SHADE(I), ICETH(I), TVOLAVG, rn(i), rs(i), ranlw(jw), rb(i), re(i), rc(i), REAER(I)*86400., (ADJUSTR(C2CH(JAC)), JAC = 1, NAC), (ADJUSTR(EPCH(JE)), JE = 1, NEP), (ADJUSTR(MACCH(JM)), JM = 1, NMC), SEDCH, SEDPCH, SEDNCH, SEDCCH, (ADJUSTR(CDCH(JAD)), JAD = 1, NACD(JW)), (ADJUSTR(KFCH(JF)), JF = 1, NAF(JW)), (ADJUSTR(APCH(JA)), JA = 1, NAL), (ADJUSTR(ANCH(JA)), JA = 1, NAL), (ADJUSTR(ALCH(JA)), JA = 1, NAL) ! SW 10/20/15
                        else
                            write(TSR(J), '(f10.3,",",19(F10.3,","),*(A,","))') JDAY, DLT, ELWS(I), T1(K, I), U(K, I), QC(I), SRON(JW)*1.06, GAMMA(K, I), DEPTHB(KB(I), I), BI(KTWB(JW), I), SHADE(I), ICETH(I), TVOLAVG, rn(i), rs(i), ranlw(jw), rb(i), re(i), rc(i), REAER(I)*86400., (ADJUSTR(C2CH(JAC)), JAC = 1, NAC), (ADJUSTR(EPCH(JE)), JE = 1, NEP), (ADJUSTR(MACCH(JM)), JM = 1, NMC), (ADJUSTR(CDCH(JAD)), JAD = 1, NACD(JW)), (ADJUSTR(KFCH(JF)), JF = 1, NAF(JW)), (ADJUSTR(APCH(JA)), JA = 1, NAL), (ADJUSTR(ANCH(JA)), JA = 1, NAL), (ADJUSTR(ALCH(JA)), JA = 1, NAL) ! SW 10/20/15
                        end if
                    else
                        if (SEDIMENT_CALC(JW)) then
                            write(TSR(J), '(f10.3,",",18(F10.3,","),*(A,","))') JDAY, DLT, ELWS(I), T1(K, I), U(K, I), QC(I), SRON(JW)*1.06, GAMMA(K, I), DEPTHB(KB(I), I), BI(KTWB(JW), I), SHADE(I), TVOLAVG, rn(i), rs(i), ranlw(jw), rb(i), re(i), rc(i), REAER(I)*86400., (ADJUSTR(C2CH(JAC)), JAC = 1, NAC), (ADJUSTR(EPCH(JE)), JE = 1, NEP), (ADJUSTR(MACCH(JM)), JM = 1, NMC), SEDCH, SEDPCH, SEDNCH, SEDCCH, (ADJUSTR(CDCH(JAD)), JAD = 1, NACD(JW)), (ADJUSTR(KFCH(JF)), JF = 1, NAF(JW)), (ADJUSTR(APCH(JA)), JA = 1, NAL), (ADJUSTR(ANCH(JA)), JA = 1, NAL), (ADJUSTR(ALCH(JA)), JA = 1, NAL) ! SW 10/20/15
                        else
                            write(TSR(J), '(f10.3,",",18(F10.3,","),*(A,","))') JDAY, DLT, ELWS(I), T1(K, I), U(K, I), QC(I), SRON(JW)*1.06, GAMMA(K, I), DEPTHB(KB(I), I), BI(KTWB(JW), I), SHADE(I), TVOLAVG, rn(i), rs(i), ranlw(jw), rb(i), re(i), rc(i), REAER(I)*86400., (ADJUSTR(C2CH(JAC)), JAC = 1, NAC), (ADJUSTR(EPCH(JE)), JE = 1, NEP), (ADJUSTR(MACCH(JM)), JM = 1, NMC), (ADJUSTR(CDCH(JAD)), JAD = 1, NACD(JW)), (ADJUSTR(KFCH(JF)), JF = 1, NAF(JW)), (ADJUSTR(APCH(JA)), JA = 1, NAL), (ADJUSTR(ANCH(JA)), JA = 1, NAL), (ADJUSTR(ALCH(JA)), JA = 1, NAL) ! SW 10/20/15
                        end if
                    end if
                else ! SW 4/4/2018
                    if (ICE_COMPUTATION) then
                        if (SEDIMENT_CALC(JW)) then
                            write(TSR(J), '(f10.3,",",19(F10.3,","),*(F10.1,","))') JDAY, DLT, ELWS(I), -99., -99., QC(I), SRON(JW)*1.06, -99., DEPTHB(KB(I), I), BI(KTWB(JW), I), SHADE(I), ICETH(I), TVOLAVG, rn(i), rs(i), ranlw(jw), rb(i), re(i), rc(i), REAER(I)*86400., (-99., JAC = 1, NAC), (-99., JE = 1, NEP), (-99., JM = 1, NMC), SED(K, I), SEDP(K, I), SEDN(K, I), SEDC(K, I), (-99., JAD = 1, NACD(JW)), (-99., JF = 1, NAF(JW)), (-99., JA = 1, NAL), (-99., JA = 1, NAL), (-99., JA = 1, NAL) ! SW 10/20/15
                        else
                            write(TSR(J), '(f10.3,",",19(F10.3,","),*(F10.1,","))') JDAY, DLT, ELWS(I), -99., -99., QC(I), SRON(JW)*1.06, -99., DEPTHB(KB(I), I), BI(KTWB(JW), I), SHADE(I), ICETH(I), TVOLAVG, rn(i), rs(i), ranlw(jw), rb(i), re(i), rc(i), REAER(I)*86400., (-99., JAC = 1, NAC), (-99., JE = 1, NEP), (-99., JM = 1, NMC), (-99., JAD = 1, NACD(JW)), (-99., JF = 1, NAF(JW)), (-99., JA = 1, NAL), (-99., JA = 1, NAL), (-99., JA = 1, NAL) ! SW 10/20/15
                        end if
                    else
                        if (SEDIMENT_CALC(JW)) then
                            write(TSR(J), '(f10.3,",",18(F10.3,","),*(F10.1,","))') JDAY, DLT, ELWS(I), -99., -99., QC(I), SRON(JW)*1.06, -99., DEPTHB(KB(I), I), BI(KTWB(JW), I), SHADE(I), TVOLAVG, rn(i), rs(i), ranlw(jw), rb(i), re(i), rc(i), REAER(I)*86400., (-99., JAC = 1, NAC), (-99., JE = 1, NEP), (-99., JM = 1, NMC), SED(K, I), SEDP(K, I), SEDN(K, I), SEDC(K, I), (-99., JAD = 1, NACD(JW)), (-99., JF = 1, NAF(JW)), (-99., JA = 1, NAL), (-99., JA = 1, NAL), (-99., JA = 1, NAL) ! SW 10/20/15
                        else
                            write(TSR(J), '(f10.3,",",18(F10.3,","),*(F10.1,","))') JDAY, DLT, ELWS(I), -99., -99., QC(I), SRON(JW)*1.06, -99., DEPTHB(KB(I), I), BI(KTWB(JW), I), SHADE(I), TVOLAVG, rn(i), rs(i), ranlw(jw), rb(i), re(i), rc(i), REAER(I)*86400., (-99., JAC = 1, NAC), (-99., JE = 1, NEP), (-99., JM = 1, NMC), (-99., JAD = 1, NACD(JW)), (-99., JF = 1, NAF(JW)), (-99., JA = 1, NAL), (-99., JA = 1, NAL), (-99., JA = 1, NAL) ! SW 10/20/15
                        end if
                    end if
                end if ! SW 4/4/2018

                if (ALGAE_TOXIN) then
                    if (ATOX_DEBUG == "ON") then
                        write(ATOXIN_DEBUG_FN, '(F10.3,",",I4,",",I4,",",<NUMATOXINS>(E12.4,","),<NUMATOXINS>(E12.4,","),<NUMATOXINS>(E12.4,","),<NAL>(E12.4,","))') JDAY, K, I, (EX_TOXIN(K, I, JA), JA = 1, numatoxins), (IN_TOXIN(K, I, JA), JA = 1, numatoxins), (CTESS(K, I, JA), JA = 1, numatoxins), (ALG(K, I, JA), JA = 1, NAL)
                    end if
                end if

            end do
        end if
    end if

! SEDIMENT DIAGENESIS FREQUENCY OUTPUT
    if (constituents .and. sediment_diagenesis .and. JDAY >= NXTSEDIAG) then
        NXTSEDIAG = NXTSEDIAG + SEDIAGFREQ
        call WriteCEMASedimentModelOutput()
        call WriteCEMASedimentFluxOutput()
    end if

    do JW = 1, NWB

!**** Inactive segments

        JB = BS(JW)
        NBL(JW) = 1
        IBPR(JW) = 1
        do I = 1, NISNP(JW) - 1
            if (CUS(JB) > ISNP(I, JW)) then
                BL(NBL(JW), JW) = I
                NBL(JW) = NBL(JW) + 1
                IBPR(JW) = I + 1
            end if
            if (ISNP(I + 1, JW) > DS(JB)) then
                JB = JB + 1
            end if
        end do
        NBL(JW) = NBL(JW) - 1

!**** Snapshots

        if (SNAPSHOT(JW)) then
            if (JDAY >= NXTMSN(JW) .or. JDAY >= SNPD(SNPDP(JW) + 1, JW)) then
                if (JDAY >= SNPD(SNPDP(JW) + 1, JW)) then
                    SNPDP(JW) = SNPDP(JW) + 1
                    NXTMSN(JW) = SNPD(SNPDP(JW), JW)
                end if
                NXTMSN(JW) = NXTMSN(JW) + SNPF(SNPDP(JW), JW)
                write(SNP(JW), 10490) W2VER, (TITLE(J), J = 1, 10)
                write(SNP(JW), 10500) "Time Parameters", MONTH, GDAY, YEAR, INT(JDAY), (JDAY - INT(JDAY))*24.0, INT(ELTMJD), (ELTMJD - INT(ELTMJD))*24.0, INT(DLTS1), KLOC, ILOC, INT(MINDLT), INT(JDMIN), (JDMIN - INT(JDMIN))*24.0, KMIN, IMIN
                if (LIMITING_DLT(JW)) then
                    write(SNP(JW), 10510) KMIN, IMIN
                end if
                write(SNP(JW), 10520) INT(DLTAV), NIT, NV
                write(SNP(JW), 10530) "Meteorological Parameters"
                write(SNP(JW), 10540) TAIR(JW), DEG, TDEW(JW), DEG, PHI(JW), CLOUD(JW), ET(DS(1)), DEG, CSHE(DS(1)), SRON(JW), DEG
                write(SNP(JW), 10550) "Inflows", "Upstream inflows"
                do JB = BS(JW), BE(JW)
                    if (UP_FLOW(JB)) then
                        write(SNP(JW), 10560) JB, KTQIN(JB), KBQIN(JB), QIN(JB), TIN(JB), DEG
                    end if
                end do
                do JB = BS(JW), BE(JW)
                    if (DIST_TRIBS(JB)) then
                        write(SNP(JW), 10570)
                        write(SNP(JW), 10580) JB, QDTR(JB), TDTR(JB), DEG
                    end if
                end do
                if (TRIBUTARIES) then
                    write(SNP(JW), 10590) (ITR(JT), JT = 1, JTT)
                    write(SNP(JW), 10600) (KTTR(JT), KBTR(JT), JT = 1, JTT)
                    write(SNP(JW), 10610) (QTR(JT), JT = 1, JTT)
                    write(SNP(JW), 10620) (TTR(JT), JT = 1, JTT)
                end if
                write(SNP(JW), 10630)
                do JB = BS(JW), BE(JW)
                    if (DN_FLOW(JB)) then
                        write(SNP(JW), 10640) JB, (QSTR(JS, JB), JS = 1, JSS(JB))
                        write(SNP(JW), 10650) QSUM(JB), (K, K = KTWB(JW), KB(DS(JB)))
                        write(SNP(JW), 10660) (QOUT(K, JB), K = KTWB(JW), KB(DS(JB)))
                        write(SNP(JW), *)
                        write(SNP(JW), "(A)") "  LAYER   DEPTH(m)   T(C) DENSITY(kg/m3)  U(m/s)     Q(m3/s)"
                        do K = KTWB(JW), KB(DS(JB))
                            write(SNP(JW), "(I5,1X,F10.2,1X,F8.2,1X,F10.3,1X,F10.5,1X,F10.3)") K, DEPTHM(K, DS(JB)), T2(K, DS(JB)), RHO(K, DS(JB)), U(K, DS(JB)), QOUT(K, JB)
                        end do
                    end if
                end do
                if (WITHDRAWALS) then
                    do JWD = 1, JWW
                        write(SNP(JW), 10670) MAX(CUS(JBWD(JWD)), IWD(JWD)), QWD(JWD)
                        if (QWD(JWD) /= 0.0) then
                            write(SNP(JW), 10680) (K, K = KTW(JWD), KBW(JWD))
                            write(SNP(JW), 10690) (QSW(K, JWD), K = KTW(JWD), KBW(JWD))

                            write(SNP(JW), *)
                            write(SNP(JW), "(A)") "  LAYER   DEPTH(m)   T(C) DENSITY(kg/m3)  Q(m3/s)"
                            do K = KTWB(JW), KB(IWD(JWD))
                                write(SNP(JW), "(I5,1X,F10.2,1X,F8.2,1X,F10.3,1X,F10.3)") K, DEPTHM(K, IWD(JWD)), T2(K, IWD(JWD)), RHO(K, IWD(JWD)), QSW(K, JWD)
                            end do

                        else
                            write(SNP(JW), 10680)
                            write(SNP(JW), 10690) QWD(JWD)
                        end if
                    end do
                end if
                if (CONSTITUENTS) then
                    write(SNP(JW), 10700) "Constituent Inflow Concentrations"
                    do JB = BS(JW), BE(JW)
                        if (UP_FLOW(JB) .and. NACIN(JB) > 0) then
                            write(SNP(JW), 10710) JB, (CNAME1(INCN(JC, JB))(1:18), CIN(INCN(JC, JB), JB), CUNIT2(INCN(JC, JB)), JC = 1, NACIN(JB))
                        end if
                        if (DIST_TRIBS(JB) .and. NACDT(JB) > 0) then
                            write(SNP(JW), 10730) JB, (CNAME1(DTCN(JC, JB))(1:18), CDTR(DTCN(JC, JB), JB), CUNIT2(DTCN(JC, JB)), JC = 1, NACDT(JB))
                        end if
                    end do
                    do JT = 1, NTR
                        if (NACTR(JT) > 0) then
                            write(SNP(JW), 10720) JT, (CNAME1(TRCN(JC, JT))(1:18), CTR(TRCN(JC, JT), JT), CUNIT2(TRCN(JC, JT)), JC = 1, NACTR(JT))
                        end if
                    end do
                end if
                if (EVAPORATION(JW) .or. PRECIPITATION(JW)) then
                    write(SNP(JW), 10740)
                end if
                if (EVAPORATION(JW)) then
                    write(SNP(JW), 10750) (JB, EVBR(JB), JB = BS(JW), BE(JW)) ! SW 9/15/05
                    write(SNP(JW), 10755) (JB, -VOLEV(JB), JB = BS(JW), BE(JW))
                end if
                if (PRECIPITATION(JW)) then
                    write(SNP(JW), 10760) (JB, PR(JB), JB = BS(JW), BE(JW))
                end if
                if (HEAD_BOUNDARY(JW)) then
                    write(SNP(JW), 10770)
                    do JB = BS(JW), BE(JW)
                        if (UH_EXTERNAL(JB)) then
                            write(SNP(JW), 10780) JB, ELUH(JB)
                        end if
                        if (DH_EXTERNAL(JB)) then
                            write(SNP(JW), 10790) JB, ELDH(JB)
                        end if
                    end do
                end if
                if (VOLUME_BALANCE(JW)) then
                    write(SNP(JW), 10800)
                    write(SNP(JW), 10810) JW, VOLSR(JW), VOLTR(JW), VOLTR(JW) - VOLSR(JW), DLVR(JW)
                    do JB = BS(JW), BE(JW)
                        if (VOLSBR(JB) /= 0.0) then
                            DLVBR = (VOLTBR(JB) - VOLSBR(JB))/VOLSBR(JB)
                        end if
                        write(SNP(JW), 10820) JB, VOLSBR(JB), VOLTBR(JB), VOLTBR(JB) - VOLSBR(JB), DLVBR*100.0
                    end do
                end if
                if (ENERGY_BALANCE(JW)) then
                    write(SNP(JW), 10830)
                    if (ESR(JW) /= 0.0) then
                        DLE = (ESR(JW) - ETR(JW))/ESR(JW)
                    end if
                    write(SNP(JW), 10840) JW, ESR(JW)*4.184E3, ETR(JW)*4.184E3, (ESR(JW) - ETR(JW))*4.184E3, DLE*100.0
                    do JB = BS(JW), BE(JW)
                        write(SNP(JW), 10870) JB
                        if (ESBR(JB) /= 0.0) then
                            DLE = (ESBR(JB) - ETBR(JB))/ESBR(JB)
                        end if
                        write(SNP(JW), 10850) ESBR(JB)*4.184E3, ETBR(JB)*4.1843E3, (ESBR(JB) - ETBR(JB))*4.1843E3, DLE*100.0
                    end do
                end if
                if (MASS_BALANCE(JW)) then
                    write(SNP(JW), 10860)
                    do JB = BS(JW), BE(JW)
                        write(SNP(JW), 10870) JB
                        do JC = 1, NAC
                            if (CMBRS(CN(JC), JB) /= 0.0) then
                                DLMR = (CMBRT(CN(JC), JB) - CMBRS(CN(JC), JB))/(CMBRS(CN(JC), JB) + NONZERO)*100.0
                            end if
                            write(SNP(JW), 10880) CNAME1(CN(JC)), CMBRS(CN(JC), JB), CUNIT1(CN(JC)), CMBRT(CN(JC), JB), CUNIT1(CN(JC)), CMBRT(CN(JC), JB) - CMBRS(CN(JC), JB), CUNIT1(CN(JC)), DLMR
                        end do

                        do M = 1, NMC
                            if (MACROPHYTE_CALC(JW, M)) then
                                if (MACMBRS(JB, M) /= 0.0) then
                                    DLMR = (MACMBRT(JB, M) - MACMBRS(JB, M))/(MACMBRS(JB, M) + NONZERO)
                                end if
                                write(SNP(JW), 3312) M, MACMBRS(JB, M), MACMBRT(JB, M), MACMBRT(JB, M) - MACMBRS(JB, M), DLMR*100.0
! 3312                            format(5X,'Macrophyte spec ',i2,/7X,'Spatially integrated mass [MACMBRS] = ',1PE15.8E2,1X,'g ',/7X, 'Temporally integrated mass [MACMBRT] = ',1PE15.8E2,1X,'g ',/7X,'Mass error = ', 1PE15.8E2,1X,'g ',/7X,'Percent error                      = ',1PE15.8E2,' %')
                            end if
                        end do

                    end do
                end if
                write(SNP(JW), 10890) "Geometry", KTWB(JW), ELKT(JW)
                write(SNP(JW), 10900) (JB, CUS(JB), JB = BS(JW), BE(JW))
                call OUTPUT(JDAY, IBPR(JW), NISNP(JW), KBR(JW), ISNP, BL(1, JW), NBL(JW))

            end if
        end if

!**** Vertical profiles

        if (PROFILE(JW)) then
            if (JDAY >= NXTMPR(JW) .or. JDAY >= PRFD(PRFDP(JW) + 1, JW)) then
                if (JDAY >= PRFD(PRFDP(JW) + 1, JW)) then
                    PRFDP(JW) = PRFDP(JW) + 1
                    NXTMPR(JW) = PRFD(PRFDP(JW), JW)
                end if
                NXTMPR(JW) = NXTMPR(JW) + PRFF(PRFDP(JW), JW)
                NSPRF(JW) = NSPRF(JW) + 1
                if (iprf(1, 1) /= -1) then ! SW 4/1/2016
                    write(PRF(JW), "(F8.3,1X,A3,I3,A,2I4,F8.4,I8)") JDAY, ADJUSTL(MONTH), GDAY, ", ", YEAR, KTWB(JW), SNGL(Z(DS(BS(JW)))), NSPRF(JW)
                    do JP = 1, NIPRF(JW)
                        NRS = KB(IPRF(JP, JW)) - KTWB(JW) + 1
                        write(PRF(JW), "(A8,I4/(8F10.2))") "TEMP    ", NRS, (T2(K, IPRF(JP, JW)), K = KTWB(JW), KB(IPRF(JP, JW)))
                    end do
                    do JC = 1, NAC
                        if (PRINT_CONST(CN(JC), JW)) then
                            do JP = 1, NIPRF(JW)
                                NRS = KB(IPRF(JP, JW)) - KTWB(JW) + 1
                                write(PRF(JW), "(A8,I4/(8(E13.6,2X)))") ADJUSTL(CNAME2(CN(JC))), NRS, (C2(K, IPRF(JP, JW), CN(JC))*CMULT(CN(JC)), K = KTWB(JW), KB(IPRF(JP, JW)))
                            end do
                        end if
                    end do
                    if (CONSTITUENTS) then
                        do JD = 1, NACD(JW)
                            do JP = 1, NIPRF(JW)
                                NRS = KB(IPRF(JP, JW)) - KTWB(JW) + 1
                                write(PRF(JW), "(A8,I4/(8(E13.6,2X)))") ADJUSTL(CDNAME2(CDN(JD, JW))), NRS, (CD(K, IPRF(JP, JW), CDN(JD, JW))*CDMULT(CDN(JD, JW)), K = KTWB(JW), KB(IPRF(JP, JW)))
                            end do
                        end do
                    end if

                else
                    if (jw == 1) then ! write out individual files on these days
                        write(SEGNUM, "(F8.2)") JDAY !   '(I0)'    INT(JDAY)
                        SEGNUM = ADJUSTL(SEGNUM)
                        L = LEN_TRIM(SEGNUM)
                        open(NUNIT, FILE="ProfLongJD" // SEGNUM(1:L) // ".csv", STATUS="UNKNOWN")
                        if (CONSTITUENTS) then
                            write(NUNIT, '(*(A,","))') "Seg#", "ElevWaterSurf(m)", "Q(m3/s)", "SurfaceTemp(oC)", "Depth(m)", "Width(m)", "VolWeighTemp(oC)", (CNAME2(CN(JC)), JC = 1, NAC), (CDNAME2(CDN(JD, JW)), JD = 1, NACD(JW))
                        else
                            write(NUNIT, '(*(A,","))') "Seg#", "ElevWaterSurf(m)", "Q(m3/s)", "SurfaceTemp(oC)", "Depth(m)", "Width(m)", "VolWeighTemp(oC)"
                        end if

                        do JJ = 1, NWB
                            K = KTWB(JJ)
                            do JB = BS(JJ), BE(JJ)
                                do I = CUS(JB), DS(JB)

! TEMP VOL WEIGHTED AVERAGE  SW 8/30/2018
                                    TVOLAVG = 0.0 ! SW 8/30/2018
                                    VOLTOT = 0.0
                                    do KK = KTWB(JW), KB(I)
                                        VOLTOT = VOLTOT + VOL(KK, I)
                                        TVOLAVG = TVOLAVG + T1(KK, I)*VOL(KK, I)
                                    end do
                                    if (KB(I) >= KTWB(JW)) then
                                        tvolavg = tvolavg/voltot
                                    end if

                                    if (CONSTITUENTS) then
                                        write(NUNIT, '(I5,",",*(F12.4,","))') I, ELWS(I), QC(I), T2(K, I), DEPTHB(KB(I), I), B(KTI(I), I), TVOLAVG, (C2(K, I, CN(JAC))*CMULT(CN(JAC)), JAC = 1, NAC), (CD(K, I, CDN(JD, JW))*CDMULT(CDN(JD, JW)), JD = 1, NACD(JW))
                                    else
                                        write(NUNIT, '(I5,",",*(F12.3,","))') I, ELWS(I), QC(I), T2(K, I), DEPTHB(KB(I), I), B(KTI(I), I), TVOLAVG
                                    end if
                                end do
                            end do
                        end do
                        close(NUNIT)
                    end if
                end if ! end of iprf /= -1

            end if
        end if

!**** Spreadsheet

        if (SPREADSHEET(JW)) then
            if (JDAY >= NXTMSP(JW) .or. JDAY >= SPRD(SPRDP(JW) + 1, JW)) then
                if (JDAY >= SPRD(SPRDP(JW) + 1, JW)) then
                    SPRDP(JW) = SPRDP(JW) + 1
                    NXTMSP(JW) = SPRD(SPRDP(JW), JW)
                end if
                CONV1 = BLANK1
                NXTMSP(JW) = NXTMSP(JW) + SPRF(SPRDP(JW), JW)
!DO J=1,NISPR(JW)
!  KBMAX(JW) = MAX(KB(ISPR(J,JW)),KBMAX(JW))
                if (SPRC(JW) == "     ONV") then ! SW 8/31/2018
                    do J = 1, NISPR(JW)
!KBMAX(JW) = MAX(KB(ISPR(J,JW)),KBMAX(JW))
                        K = KTWB(JW)
! TEMP VOL WEIGHTED AVERAGE  SW 8/30/2018
                        TVOLAVG = 0.0 ! SW 8/30/2018
                        VOLTOT = 0.0
                        do KK = KTWB(JW), KB(ISPR(J, JW)) ! cb 11/13/18
                            VOLTOT = VOLTOT + VOL(KK, ISPR(J, JW))
                            TVOLAVG = TVOLAVG + T1(KK, ISPR(J, JW))*VOL(KK, ISPR(J, JW))
                        end do
                        if (KB(ISPR(J, JW)) >= KTWB(JW)) then
                            tvolavg = tvolavg/voltot
                        else
                            TVOLAVG = -99.0
                        end if
                        write(CONV1(K, ISPR(J, JW)), "(F12.3)") TVOLAVG
                    end do
                    write(SPRV(JW), '(A,",",F10.3,",",*(A,","))') "Temperature", JDAY, (CONV1(K, ISPR(J, JW)), J = 1, NISPR(JW))

                    do JC = 1, NAC
                        if (PRINT_CONST(CN(JC), JW)) then
                            do J = 1, NISPR(JW)
                                TVOLAVG = 0.0 ! SW 8/30/2018
                                VOLTOT = 0.0
                                do KK = KTWB(JW), KB(ISPR(J, JW)) ! cb 11/13/18
                                    VOLTOT = VOLTOT + VOL(KK, ISPR(J, JW))
                                    TVOLAVG = TVOLAVG + C2(KK, ISPR(J, JW), CN(JC))*CMULT(CN(JC))*VOL(KK, ISPR(J, JW))
                                end do
                                if (KB(ISPR(J, JW)) >= KTWB(JW)) then
                                    tvolavg = tvolavg/voltot
                                else
                                    TVOLAVG = -99.0
                                end if
                                write(CONV1(K, ISPR(J, JW)), "(F12.4)") TVOLAVG ! SW 8/31/18
                            end do
                            write(SPRV(JW), '(A38,",",F10.3,",",*(A,","))') CNAME3(CN(JC)), JDAY, (CONV1(K, ISPR(J, JW)), J = 1, NISPR(JW))
                        end if
                    end do
                    if (CONSTITUENTS) then
                        do JD = 1, NACD(JW)
                            if (PRINT_DERIVED(CDN(JD, JW), JW)) then
                                do J = 1, NISPR(JW)
                                    TVOLAVG = 0.0 ! SW 8/30/2018
                                    VOLTOT = 0.0
                                    do KK = KTWB(JW), KB(ISPR(J, JW)) ! cb 11/13/18
                                        VOLTOT = VOLTOT + VOL(KK, ISPR(J, JW))
                                        TVOLAVG = TVOLAVG + CD(KK, ISPR(J, JW), CDN(JD, JW))*CDMULT(CDN(JD, JW))*VOL(KK, ISPR(J, JW))
                                    end do
                                    if (KB(ISPR(J, JW)) >= KTWB(JW)) then
                                        tvolavg = tvolavg/voltot
                                    else
                                        TVOLAVG = -99.0
                                    end if
                                    write(CONV1(K, ISPR(J, JW)), "(F12.4)") TVOLAVG ! SW 8/31/18
                                end do
                                write(SPRV(JW), '(A38,",",F10.3,",",*(A,","))') CDNAME3(CDN(JD, JW)), JDAY, (CONV1(K, ISPR(J, JW)), J = 1, NISPR(JW))
                            end if
                        end do
                    end if
                end if ! END SECTION ONV
                do J = 1, NISPR(JW)
                    KBMAX(JW) = MAX(KB(ISPR(J, JW)), KBMAX(JW))
                    do K = KTWB(JW), KB(ISPR(J, JW))
                        write(CONV1(K, J), "(F12.3)") T2(K, ISPR(J, JW))
                    end do

                end do
                do K = KTWB(JW), KBMAX(JW)
                    write(SPR(JW), '(A,",",2(F10.3,","),*(F10.3,",",A,","))') "Temperature", JDAY, DEPTHM(K, DS(BS(JW))), (ELWS(ISPR(J, JW)) - DEPTHM(K, ISPR(J, JW)), CONV1(K, J), J = 1, NISPR(JW))
                end do
                do JC = 1, NAC
                    if (PRINT_CONST(CN(JC), JW)) then
                        do J = 1, NISPR(JW)
                            do K = KTWB(JW), KB(ISPR(J, JW))
                                write(CONV1(K, J), "(F12.4)") C2(K, ISPR(J, JW), CN(JC))*CMULT(CN(JC)) ! SW 8/13/06
                            end do
                        end do
                        do K = KTWB(JW), KBMAX(JW)
                            write(SPR(JW), '(A38,",",2(F10.3,","),*(F10.3,",",A,","))') CNAME3(CN(JC)), JDAY, DEPTHM(K, DS(BS(JW))), (ELWS(ISPR(J, JW)) - DEPTHM(K, ISPR(J, JW)), CONV1(K, J), J = 1, NISPR(JW))
                        end do
                    end if
                end do
                if (CONSTITUENTS) then
                    do JD = 1, NACD(JW)
                        if (PRINT_DERIVED(CDN(JD, JW), JW)) then
                            do J = 1, NISPR(JW)
                                do K = KTWB(JW), KB(ISPR(J, JW))
                                    write(CONV1(K, J), "(F12.4)") CD(K, ISPR(J, JW), CDN(JD, JW))*CDMULT(CDN(JD, JW)) ! SW 8/13/06
                                end do
                            end do
                            do K = KTWB(JW), KBMAX(JW)
                                write(SPR(JW), '(A38,",",2(F10.3,","),*(F10.3,",",A,","))') CDNAME3(CDN(JD, JW)), JDAY, DEPTHM(K, DS(BS(JW))), (ELWS(ISPR(J, JW)) - DEPTHM(K, ISPR(J, JW)), CONV1(K, J), J = 1, NISPR(JW))
                            end do
                        end if
                    end do
                end if
            end if
        end if

!**** Contours

        if (CONTOUR(JW)) then
            if (JDAY >= NXTMCP(JW) .or. JDAY >= CPLD(CPLDP(JW) + 1, JW)) then
                if (JDAY >= CPLD(CPLDP(JW) + 1, JW)) then
                    CPLDP(JW) = CPLDP(JW) + 1
                    NXTMCP(JW) = CPLD(CPLDP(JW), JW)
                end if
                NXTMCP(JW) = NXTMCP(JW) + CPLF(CPLDP(JW), JW)

                if (TECPLOT(JW) /= "      ON") then
                    write(CPL(JW), "(A,F12.4,5X,A9,5X,I2,5X,I4)") "New date ", JDAY, MONTH, GDAY, YEAR
                    write(CPL(JW), "(9(I8,2X))") KTWB(JW)
                    write(CPL(JW), "(9(E13.6,2X))") (QTR(JT), JT = 1, NTR)
                    write(CPL(JW), "(9(E13.6,2X))") (TTR(JT), JT = 1, NTR)
                    do JT = 1, NTR
                        do JAC = 1, NACTR(JT)
                            if (PRINT_CONST(TRCN(JAC, JT), JW)) then
                                write(CPL(JW), "(9(E13.6,2X))") CTR(TRCN(JAC, JT), JT)
                            end if
                        end do
                    end do
                    do JB = BS(JW), BE(JW)
                        if (BR_INACTIVE(JB)) then
                            cycle
                        end if ! SW 7/1/2019
                        write(CPL(JW), "(9(I8,2X))") CUS(JB)
                        write(CPL(JW), "(9(E13.6,2X))") QIN(JB), QSUM(JB)
                        do I = CUS(JB), DS(JB)
                            write(CPL(JW), "(A38/(9(E13.6,2X)))") "BHR", (BHR1(K, I), K = KTWB(JW) + 1, KB(I))
                        end do
                        do I = CUS(JB), DS(JB)
                            write(CPL(JW), "(A38/(9(E13.6,2X)))") "U", (U(K, I), K = KTWB(JW), KB(I))
                        end do
                        write(CPL(JW), "(A38/(9(E13.6,2X)))") "QC", (QC(I), I = CUS(JB), DS(JB))
                        write(CPL(JW), "(A38/(9(E13.6,2X)))") "Z", (Z(I), I = CUS(JB), DS(JB))
                        write(CPL(JW), "(A38/(9(I8,2X)))") "KTI", (kti(I), I = CUS(JB), DS(JB)) ! v3.5
                        do I = CUS(JB), DS(JB)
                            write(CPL(JW), "(A38/(9(E13.6,2X)))") "Temperature", (T2(K, I), K = KTWB(JW), KB(I))
                        end do
                        do JC = 1, NAC
                            if (PRINT_CONST(CN(JC), JW)) then
                                do I = CUS(JB), DS(JB)
                                    write(CPL(JW), "(A38/(9(E13.6,2X)))") CNAME(CN(JC)), (C2(K, I, CN(JC))*CMULT(CN(JC)), K = KTWB(JW), KB(I))
                                end do
                            end if
                        end do
                        do JE = 1, NEP
                            do I = CUS(JB), DS(JB)
                                if (PRINT_EPIPHYTON(JW, JE)) then
                                    write(CPL(JW), "(A38/(9(E13.6,2X)))") "Epiphyton", (EPD(K, I, JE), K = KTWB(JW), KB(I))
                                end if
                            end do
                        end do
                        if (PRINT_SEDIMENT(JW)) then
                            do I = CUS(JB), DS(JB)
                                write(CPL(Jw), "(A38/(9(E13.6,2X)))") "Sediment", (seD(K, I), K = KTWB(JW), KB(I))
                            end do
                            do I = CUS(JB), DS(JB)
                                write(CPL(Jw), "(A38/(9(E13.6,2X)))") "Sediment P", (seDp(K, I), K = KTWB(JW), KB(I))
                            end do
                            do I = CUS(JB), DS(JB)
                                write(CPL(Jw), "(A38/(9(E13.6,2X)))") "Sediment N", (seDn(K, I), K = KTWB(JW), KB(I))
                            end do
                            do I = CUS(JB), DS(JB)
                                write(CPL(Jw), "(A38/(9(E13.6,2X)))") "Sediment C", (seDc(K, I), K = KTWB(JW), KB(I))
                            end do
                        end if
                        do M = 1, NMC
                            if (PRINT_MACROPHYTE(JW, M)) then
                                do I = CUS(JB), DS(JB)
                                    write(CPL(Jw), "(A38/(9(E13.6,2X)))") "Macrophytes", ((macrc(j, K, I, m), j = kti(i), kb(i)), K = KTwb(Jw), KB(I))
                                end do
                            end if
                        end do
                        if (CONSTITUENTS) then
                            do JD = 1, NACD(JW)
                                if (PRINT_DERIVED(CDN(JD, JW), JW)) then
                                    do I = CUS(JB), DS(JB)
                                        write(CPL(JW), "(A38/(9(F10.3,2X)))") CDNAME(CDN(JD, JW)), (CD(K, I, CDN(JD, JW))*CDMULT(CDN(JD, JW)), K = KTWB(JW), KB(I)) ! cb 6/28/13
                                    end do
!WRITE (CPL(JW),'(A38/(9(F10.3,2X)))') CDNAME(CDN(JD,JW)),((CD(K,I,CDN(JD,JW))*CDMULT(CDN(JD,JW)),             &        ! SW 8/12/06
!K=KTWB(JW),KB(I)),I=CUS(JB),DS(JB))  ! CB 1/03/05
                                end if
                            end do
                        end if
                    end do
                else
!         ICPL=ICPL+1
                    ICPL(JW) = ICPL(JW) + 1 ! cb 1/26/09
                    ITOT = 0
!         do jb=1,nbr
                    do JB = BS(JW), BE(JW) ! cb 1/26/09
                        if (BR_INACTIVE(JB)) then
                            cycle
                        end if ! SW 7/1/2019
                        if (BR_NOTECPLOT(JB)) then
                            cycle
                        end if ! SW 8/27/2019
                        ITOT = ITOT + DS(JB) - CUS(JB) + 2
                    end do
                    write(CPL(JW), 9864) JDAY, KMX - KTWB(JW) + 2, ITOT
                    9864 format('ZONE T="',f9.3,'"',' I=',I3,' J=',I3,' F=POINT')
                    do JB = BS(JW), BE(JW)
                        if (BR_INACTIVE(JB)) then
                            cycle
                        end if ! SW 7/1/2019
                        if (BR_NOTECPLOT(JB)) then
                            cycle
                        end if ! SW 8/27/2019
                        do I = CUS(JB), DS(JB) + 1
                            K = KTWB(JW) ! PRINT AN EXTRA LINE FOR THE SURFACE
                            if (I /= DS(JB) + 1) then
                                if (HABTATC == "      ON") then
                                    write(CPL(JW), 9999) X1(I), ELWS(I), U(K, I), -W(K, I), T1(K, I), RHO(K, I), HAB(K, I), (C2(K, I, CN(JC)), JC = 1, NAC), (CD(K, I, CDN(JD, JW)), JD = 1, NACD(JW)) ! SW 1/17/17
                                else
                                    write(CPL(JW), 9999) X1(I), ELWS(I), U(K, I), -W(K, I), T1(K, I), RHO(K, I), (C2(K, I, CN(JC)), JC = 1, NAC), (CD(K, I, CDN(JD, JW)), JD = 1, NACD(JW)) ! SW 1/17/17
                                end if

                            else
                                XDUM = -99.0
                                if (HABTATC == "      ON") then ! SW 7/15/14
                                    write(CPL(JW), 9999) X1(I), ELWS(I), XDUM, XDUM, XDUM, XDUM, XDUM, (XDUM, JJ = 1, NAC), (XDUM, JJ = 1, NACD(JW))
                                else ! SW 7/15/14
                                    write(CPL(JW), 9999) X1(I), ELWS(I), XDUM, XDUM, XDUM, XDUM, (XDUM, JJ = 1, NAC), (XDUM, JJ = 1, NACD(JW)) ! SW 1/17/17 7/15/14
                                end if ! SW 7/15/14
                            end if
                            do K = KTWB(JW), KMX - 1
                                if (I /= DS(JB) + 1 .and. K <= KB(I)) then
                                    if (HABTATC == "      ON") then
                                        write(CPL(JW), 9999) X1(I), ELWS(I) - DEPTHM(K, I), U(K, I), -W(K, I), T1(K, I), RHO(K, I), HAB(K, I), (C2(K, I, CN(JC)), JC = 1, NAC), (CD(K, I, CDN(JD, JW)), JD = 1, NACD(JW)) ! SW 1/17/17
                                    else
                                        write(CPL(JW), 9999) X1(I), ELWS(I) - DEPTHM(K, I), U(K, I), -W(K, I), T1(K, I), RHO(K, I), (C2(K, I, CN(JC)), JC = 1, NAC), (CD(K, I, CDN(JD, JW)), JD = 1, NACD(JW)) ! SW 1/17/17
                                    end if

                                    if (K == KB(I)) then
                                        if (HABTATC == "      ON") then
                                            write(CPL(JW), 9999) X1(I), ELWS(I) - DEPTHB(K, I), U(K, I), -W(K, I), T1(K, I), RHO(K, I), HAB(K, I), (C2(K, I, CN(JC)), JC = 1, NAC), (CD(K, I, CDN(JD, JW)), JD = 1, NACD(JW)) ! SW 1/17/17
                                        else
                                            write(CPL(JW), 9999) X1(I), ELWS(I) - DEPTHB(K, I), U(K, I), -W(K, I), T1(K, I), RHO(K, I), (C2(K, I, CN(JC)), JC = 1, NAC), (CD(K, I, CDN(JD, JW)), JD = 1, NACD(JW)) ! SW 1/17/17
                                        end if

                                    end if
                                else
                                    XDUM = -99.0
                                    if (HABTATC == "      ON") then
                                        write(CPL(JW), 9999) X1(I), ELWS(I - 1) - DEPTHM(K, I - 1), XDUM, XDUM, XDUM, XDUM, XDUM, (XDUM, JJ = 1, NAC), (XDUM, JJ = 1, NACD(JW))
                                    else
                                        write(CPL(JW), 9999) X1(I), ELWS(I - 1) - DEPTHM(K, I - 1), XDUM, XDUM, XDUM, XDUM, (XDUM, JJ = 1, NAC), (XDUM, JJ = 1, NACD(JW)) ! SW 7/15/14
                                    end if

                                    if (K == KB(I)) then
                                        if (HABTATC == "      ON") then
                                            write(CPL(JW), 9999) X1(I), ELWS(I - 1) - DEPTHB(K, I - 1), XDUM, XDUM, XDUM, XDUM, XDUM, (XDUM, JJ = 1, NAC), (XDUM, JJ = 1, NACD(JW))
                                        else
                                            write(CPL(JW), 9999) X1(I), ELWS(I - 1) - DEPTHB(K, I - 1), XDUM, XDUM, XDUM, XDUM, (XDUM, JJ = 1, NAC), (XDUM, JJ = 1, NACD(JW)) ! SW 7/15/14
                                        end if

                                    end if
                                end if
                            end do
                        end do
                    end do
                    write(CPL(JW), 9899) ICPL(JW), IMON, GDAY, YEAR ! cb 1/26/09
                    9899 format('TEXT X=0.75, y=0.85, H=2.8,ZN=',i4,',',' C=BLACK,','T= "',i2,'/',i2,'/',i4,'"')
                    write(CPL(JW), 9863) ICPL(JW), JDAY ! cb 1/26/09
                    9863 format('TEXT X=0.75, y=0.90, H=2.8,ZN=',i4,',',' C=BLACK,','T= "Julian Day ',f9.3,'"')
                    9999 format(200(e13.6,1x))
                end if
            end if
        end if

!**** Fluxes   KF is the instantaneous flux in g/m3/s, KFS is the summed flux in g eventually divided by elapsed time between calls to FLUX output and converted to kg below


        if (FLUX(JW)) then
            if (JDAY >= NXTMFL(JW) .or. JDAY >= FLXD(FLXDP(JW) + 1, JW)) then
                if (JDAY >= FLXD(FLXDP(JW) + 1, JW)) then
                    FLXDP(JW) = FLXDP(JW) + 1
                    NXTMFL(JW) = FLXD(FLXDP(JW), JW)
                end if

                NLINES = 0 ! SW 3/8/16
                NXTMFL(JW) = NXTMFL(JW) + FLXF(FLXDP(JW), JW)
                CONV = BLANK
                do JAF = 1, NAF(JW)
                    if (ELTMF(JW) > 0.0) then
                        do JB = BS(JW), BE(JW)
                            do I = CUS(JB), DS(JB)
                                do K = KTWB(JW), KB(I)
                                    KFS(K, I, KFCN(JAF, JW)) = DAY*KFS(K, I, KFCN(JAF, JW))/(1000.*ELTMF(JW)) ! KFS IN G, 86400 S/D * G /ELAPSED TIME IN S/1000 G/KG == KG/D
                                    KFJW(JW, KFCN(JAF, JW)) = KFJW(JW, KFCN(JAF, JW)) + KFS(K, I, KFCN(JAF, JW)) ! SUM UP FOR ENTIRE WATERBODY
                                end do
                            end do
                        end do
                    end if
                    do I = 1, NISNP(JW)
                        do K = KTWB(JW), KB(ISNP(I, JW))
                            write(CONV(K, I), "(E10.3)") KFS(K, ISNP(I, JW), KFCN(JAF, JW)) ! KG/D
                        end do
                    end do
                    if (NEW_PAGE) then
                        write(FLX(JW), "(/(A72))") (TITLE(J), J = 1, 11)
                        NLINES = KMX - KTWB(JW) + 14
                        NEW_PAGE = .false.
                    end if
                    NLINES = NLINES + KMX - KTWB(JW) + 11
                    NEW_PAGE = NLINES > 72
                    write(FLX(JW), "(/A,F10.3,X,3(A,I0),A,F0.2,A/)") "New date ", JDAY, MONTH // " ", GDAY, ", ", YEAR, "   Julian Date = ", INT(JDAY), " days ", (JDAY - INT(JDAY))*24.0, " hours           " // KFNAME(KFCN(JAF, JW))
                    write(FLX(JW), "(3X,2000I10)") (ISNP(I, JW), I = 1, NISNP(JW))
                    do K = KTWB(JW), KBR(JW)
                        write(FLX(JW), "(1X,I2,200A)") K, (CONV(K, I), I = 1, NISNP(JW))
                    end do
                end do
                write(FLX2(JW), '(F10.3,",",f8.3,",",*(E12.4,","))') JDAY, ELTMF(JW)/DAY, (KFJW(JW, KFCN(K, JW)), K = 1, NAF(JW))
                ELTMF(JW) = 0.0
                KF(:, CUS(BS(JW)):DS(BE(JW)), KFCN(1:NAF(JW), JW)) = 0.0
                KFS(:, CUS(BS(JW)):DS(BE(JW)), KFCN(1:NAF(JW), JW)) = 0.0
                KFJW(JW, KFCN(1:NAF(JW), JW)) = 0.0
            end if
        end if

    end do

!** Downstream flow, temperature, and constituent files

    if (DOWNSTREAM_OUTFLOW) then
        if (JDAY >= NXTMWD .or. JDAY >= WDOD(WDODP + 1)) then
            if (JDAY >= WDOD(WDODP + 1)) then
                WDODP = WDODP + 1
                NXTMWD = WDOD(WDODP)
            end if
            if (WDOC == "      ON") then
                NXTMWD = NXTMWD + WDOF(WDODP)
            else
                if (WDOC == "     ONS") then
                    NXTMWD_SEC = NXTMWD_SEC + WDOF(WDODP) ! cb 4/6/18 frequency test seconds
                    NXTMWD = NXTMWD_SEC/86400.0
                else
                    if (WDOC == "     ONH") then
                        NXTMWD_SEC = NXTMWD_SEC + WDOF(WDODP)*3600. ! cb 4/6/18 frequency test hourly
                        NXTMWD = NXTMWD_SEC/86400.0
                    end if
                end if
            end if

            JFILE = 0

            do J = 1, NIWDO
                QWDO(J) = 0.0
                TWDO(J) = 0.0
                CWDO(:, J) = 0.0
                CDWDO(:, J) = 0.0
                CDTOT = 0.0
                NUMOUTLETS = 0
                JWD = 0 ! SW 5/17/13

                do JW = 1, NWB
                    do JB = BS(JW), BE(JW)
                        if (DS(JB) == IWDO(J)) then

                            do JS = 1, JSS(JB) !NSTR(JB)
                                NUMOUTLETS = NUMOUTLETS + 1
                                if (QSTR(JS, JB) == 0.0) then
                                    TAVG(JS, JB) = -99.0
                                    CAVG(JS, JB, :) = -99.0
                                    CDAVG(JS, JB, :) = -99.0
                                end if
                                QOUTLET(NUMOUTLETS) = QSTR(JS, JB)
                                TOUTLET(NUMOUTLETS) = TAVG(JS, JB)
                            end do

! cb 1/16/13 removed old code

! OUTPUT INDIVIDUAL FILES
                            do JS = 1, NSTR(JB)
                                JFILE = JFILE + 1
                                qwdo(j) = qwdo(j) + qstr(js, jb) ! cb 1/16/13
                                twdo(j) = twdo(j) + qstr(js, jb)*tavg(js, jb)
                                write(WDO2(JFILE, 1), '(F10.3,",",F10.4)') JDAY, QSTR(JS, JB)
                                write(WDO2(JFILE, 2), '(F10.3,",",F8.2)') JDAY, TAVG(JS, JB)
                                if (CONSTITUENTS) then
                                    cwdo(cn(1:nac), j) = cwdo(cn(1:nac), j) + qstr(js, jb)*cavg(js, jb, cn(1:nac)) ! cb 1/16/13
                                    write(WDO2(JFILE, 3), '(F10.3,",",*(F10.4,","))') JDAY, (CAVG(JS, JB, CN(JC)), JC = 1, NAC)
                                end if
                                if (DERIVED_CALC) then
                                    cdwdo(cdn(1:nacd(jw), jw), j) = cdwdo(cdn(1:nacd(jw), jw), j) + qstr(js, jb)*cdavg(js, jb, cdn(1:nacd(jw), jw)) ! cb 1/16/13
                                    write(WDO2(JFILE, 4), '(F10.3,",",*(F10.4,","))') JDAY, (CDAVG(JS, JB, CDN(JD, JW)), JD = 1, NACD(JW))
                                end if
                            end do
                            JSSS(JB) = NSTR(JB)
                        end if
                    end do
                end do

! OUTPUT INDIVIDUAL FILES
! Order: spillways NSP, pumps NPU, gates NGT, pipes NPI
!     JWD=0                        ! sw 5/17/13
                do jw = 1, nwb ! cb 1/16/13
                    if (iwdo(j) >= us(bs(jw)) .and. iwdo(j) <= ds(be(jw))) then
                        exit
                    end if
                end do
                do JJ = 1, NWD
                    JWD = JWD + 1
                    if (QWD(JWD) == 0.0) then
                        TAVGW(JWD) = -99.0
                        CAVGW(JWD, :) = -99.0
                        CDAVGW(JWD, :) = -99.0
                    end if
                    if (IWD(JWD) == IWDO(J)) then
                        JFILE = JFILE + 1
                        qwdo(j) = qwdo(j) + qwd(jwd) ! cb 1/16/13
                        twdo(j) = twdo(j) + qwd(jwd)*tavgw(jwd) ! cb 1/16/13
                        write(WDO2(JFILE, 1), '(F10.3,",",F10.4)') JDAY, QWD(JWD)
                        write(WDO2(JFILE, 2), '(F10.3,",",F8.2)') JDAY, TAVGW(JWD)
                        if (CONSTITUENTS) then
                            cwdo(cn(1:nac), j) = cwdo(cn(1:nac), j) + qwd(jwd)*cavgw(jwd, cn(1:nac)) ! cb 1/16/13
                            write(WDO2(JFILE, 3), '(F10.3,",",*(F10.4,","))') JDAY, (CAVGW(JWD, CN(JC)), JC = 1, NAC)
                        end if
                        if (DERIVED_CALC) then
!DO JW=1,NWB
!  IF (IWDO(J) >= US(BS(JW)) .AND. IWDO(J) <= DS(BE(JW))) EXIT
!END DO
                            cdwdo(cdn(1:nacd(jw), jw), j) = cdwdo(cdn(1:nacd(jw), jw), j) + qwd(jwd)*cdavgw(jwd, cdn(1:nacd(jw), jw)) ! cb 1/16/13
                            write(WDO2(JFILE, 4), '(F10.3,",",*(F10.4,","))') JDAY, (CDAVGW(JWD, CDN(JD, JW)), JD = 1, NACD(JW))
                        end if
                    end if
                end do
                do JS = 1, NSP ! spillways
                    if (LATERAL_SPILLWAY(JS)) then
                        JWD = JWD + 1
                    else
                        JSSS(JBUSP(JS)) = JSSS(JBUSP(JS)) + 1
                    end if

                    if (IWDO(J) == IUSP(JS)) then
                        JFILE = JFILE + 1
                        write(WDO2(JFILE, 1), '(F10.3,",",F10.4)') JDAY, QSP(JS)
                        if (LATERAL_SPILLWAY(JS)) then
!  JWD=JWD+1
                            qwdo(j) = qwdo(j) + qsp(js) ! cb 1/16/13
                            twdo(j) = twdo(j) + qsp(js)*tavgw(jwd) ! cb 1/16/13
                            if (QSP(JS) == 0.0) then
                                TAVGW(JWD) = -99.0
                                CAVGW(JWD, :) = -99.0
                                CDAVGW(JWD, :) = -99.0
                            end if
                            write(WDO2(JFILE, 2), '(F10.3,",",F8.2)') JDAY, TAVGW(JWD)
                            if (CONSTITUENTS) then
                                cwdo(cn(1:nac), j) = cwdo(cn(1:nac), j) + qsp(js)*cavgw(jwd, cn(1:nac)) ! cb 1/16/13
                                write(WDO2(JFILE, 3), '(F10.3,",",*(F10.4,","))') JDAY, (CAVGW(JWD, CN(JC)), JC = 1, NAC)
                            end if
                            if (DERIVED_CALC) then
                                cdwdo(cdn(1:nacd(jw), jw), j) = cdwdo(cdn(1:nacd(jw), jw), j) + qsp(js)*cdavgw(jwd, cdn(1:nacd(jwusp(js)), jwusp(js))) ! cb 1/16/13
                                write(WDO2(JFILE, 4), '(F10.3,",",*(F10.4,","))') JDAY, (CDAVGW(JWD, CDN(JD, JWUSP(JS))), JD = 1, NACD(JWUSP(JS)))
                            end if
                        else
!  JSSS(JBUSP(JS))=JSSS(JBUSP(JS))+1
                            qwdo(j) = qwdo(j) + qsp(js) ! cb 1/16/13
                            twdo(j) = twdo(j) + qsp(js)*tavg(jsss(jbusp(js)), jbusp(js)) ! cb 1/16/13
                            if (QSP(JS) == 0.0) then
                                TAVG(JSSS(JBUSP(JS)), JBUSP(JS)) = -99.0
                                CAVG(JSSS(JBUSP(JS)), JBUSP(JS), :) = -99.0
                                CDAVG(JSSS(JBUSP(JS)), JBUSP(JS), :) = -99.0
                            end if
                            write(WDO2(JFILE, 2), '(F10.3,",",F8.2)') JDAY, TAVG(JSSS(JBUSP(JS)), JBUSP(JS))
                            if (CONSTITUENTS) then
                                cwdo(cn(1:nac), j) = cwdo(cn(1:nac), j) + qsp(js)*cavg(jsss(jbusp(js)), jbusp(js), cn(1:nac)) ! cb 1/16/13
                                write(WDO2(JFILE, 3), '(F10.3,",",*(F10.4,","))') JDAY, (CAVG(JSSS(JBUSP(JS)), JBUSP(JS), CN(JC)), JC = 1, NAC)
                            end if
                            if (DERIVED_CALC) then
                                cdwdo(cdn(1:nacd(jw), jw), j) = cdwdo(cdn(1:nacd(jw), jw), j) + qsp(js)*cdavg(jsss(jbusp(js)), jbusp(js), cdn(1:nacd(jwusp(js)), jwusp(js))) ! cb 1/16/13
                                write(WDO2(JFILE, 4), '(F10.3,",",*(F10.4,","))') JDAY, (CDAVG(JSSS(JBUSP(JS)), JBUSP(JS), CDN(JD, JWUSP(JS))), JD = 1, NACD(JWUSP(JS)))
                            end if
                        end if
                    end if
                end do
                do JS = 1, NPU ! PUMP
                    if (LATERAL_PUMP(JS)) then
                        JWD = JWD + 1
                    else
                        JSSS(JBUPU(JS)) = JSSS(JBUPU(JS)) + 1
                    end if
                    if (IWDO(J) == IUPU(JS)) then
                        JFILE = JFILE + 1
                        if (PUMPON(JS)) then
                            write(WDO2(JFILE, 1), '(F10.3,",",F10.4)') JDAY, QPU(JS)
                        else
                            write(WDO2(JFILE, 1), '(F10.3,",",F8.3)') JDAY, 0.0
                        end if
                        if (LATERAL_PUMP(JS)) then
!  JWD=JWD+1
                            if (QPU(JS) == 0.0) then
                                TAVGW(JWD) = -99.0
                                CAVGW(JWD, :) = -99.0
                                CDAVGW(JWD, :) = -99.0
                            end if
                            if (pumpon(js)) then ! cb 1/16/13
                                qwdo(j) = qwdo(j) + qpu(js)
                                twdo(j) = twdo(j) + qpu(js)*tavgw(jwd)
                            end if
                            write(WDO2(JFILE, 2), '(F10.3,",",F8.2)') JDAY, TAVGW(JWD)
! Debug
!WRITE(WDO2(JFILE,2),'(F10.3,",",F8.2,",",f8.3,",",i5,",",i5)')JDAY,TAVGW(JWD),qpu(js),js,jwd      ! Debug
                            if (CONSTITUENTS) then
                                if (pumpon(js)) then ! cb 1/16/13
                                    cwdo(cn(1:nac), j) = cwdo(cn(1:nac), j) + qpu(js)*cavgw(jwd, cn(1:nac))
                                end if
                                write(WDO2(JFILE, 3), '(F10.3,",",*(F10.4,","))') JDAY, (CAVGW(JWD, CN(JC)), JC = 1, NAC)
                            end if
                            if (DERIVED_CALC) then
                                if (pumpon(js)) then ! cb 1/16/13
                                    cdwdo(cdn(1:nacd(jw), jw), j) = cdwdo(cdn(1:nacd(jw), jw), j) + qpu(js)*cdavgw(jwd, cdn(1:nacd(jwupu(js)), jwupu(js)))
                                end if
                                write(WDO2(JFILE, 4), '(F10.3,",",*(F10.4,","))') JDAY, (CDAVGW(JWD, CDN(JD, JWUPU(JS))), JD = 1, NACD(JWUPU(JS)))
                            end if
                        else
!  JSSS(JBUPU(JS))=JSSS(JBUPU(JS))+1
                            if (QPU(JS) == 0.0) then
                                TAVG(JSSS(JBUPU(JS)), JBUPU(JS)) = -99.0
                                CAVG(JSSS(JBUPU(JS)), JBUPU(JS), :) = -99.0
                                CDAVG(JSSS(JBUPU(JS)), JBUPU(JS), :) = -99.0
                            end if
                            if (pumpon(js)) then ! cb 1/16/13
                                qwdo(j) = qwdo(j) + qpu(js)
                                twdo(j) = twdo(j) + qpu(js)*tavg(jsss(jbupu(js)), jbupu(js))
                            end if
                            write(WDO2(JFILE, 2), '(F10.3,",",F8.2)') JDAY, TAVG(JSSS(JBUPU(JS)), JBUPU(JS))
                            if (CONSTITUENTS) then
                                if (pumpon(js)) then ! cb 1/16/13
                                    cwdo(cn(1:nac), j) = cwdo(cn(1:nac), j) + qpu(js)*cavg(jsss(jbupu(js)), jbupu(js), cn(1:nac))
                                end if
                                write(WDO2(JFILE, 3), '(F10.3,",",*(F10.4,","))') JDAY, (CAVG(JSSS(JBUPU(JS)), JBUPU(JS), CN(JC)), JC = 1, NAC)
                            end if
                            if (DERIVED_CALC) then
                                if (pumpon(js)) then ! cb 1/16/13
                                    cdwdo(cdn(1:nacd(jw), jw), j) = cdwdo(cdn(1:nacd(jw), jw), j) + qpu(js)*cdavg(jsss(jbupu(js)), jbupu(js), cdn(1:nacd(jwupu(js)), jwupu(js)))
                                end if
                                write(WDO2(JFILE, 4), '(F10.3,",",*(F10.4,","))') JDAY, (CDAVG(JSSS(JBUPU(JS)), JBUPU(JS), CDN(JD, JWUPU(JS))), JD = 1, NACD(JWUPU(JS)))
                            end if
                        end if
                    end if
                end do
                do JS = 1, NPI ! pipes
                    if (LATERAL_PIPE(JS)) then
                        JWD = JWD + 1
                    else
                        JSSS(JBUPI(JS)) = JSSS(JBUPI(JS)) + 1
                    end if
                    if (IWDO(J) == IUPI(JS)) then
                        JFILE = JFILE + 1
                        write(WDO2(JFILE, 1), '(F10.3,",",F10.4)') JDAY, QPI(JS)
                        if (LATERAL_PIPE(JS)) then
!   JWD=JWD+1
                            qwdo(j) = qwdo(j) + qpi(js) ! cb 1/16/13
                            twdo(j) = twdo(j) + qpi(js)*tavgw(jwd)
                            if (QPI(JS) == 0.0) then
                                TAVGW(JWD) = -99.0
                                CAVGW(JWD, :) = -99.0
                                CDAVGW(JWD, :) = -99.0
                            end if
                            write(WDO2(JFILE, 2), '(F10.3,",",F8.2)') JDAY, TAVGW(JWD)
                            if (CONSTITUENTS) then
                                cwdo(cn(1:nac), j) = cwdo(cn(1:nac), j) + qpi(js)*cavgw(jwd, cn(1:nac)) ! cb 1/16/13
                                write(WDO2(JFILE, 3), '(F10.3,",",*(F10.4,","))') JDAY, (CAVGW(JWD, CN(JC)), JC = 1, NAC)
                            end if
                            if (DERIVED_CALC) then
                                cdwdo(cdn(1:nacd(jw), jw), j) = cdwdo(cdn(1:nacd(jw), jw), j) + qpi(js)*cdavgw(jwd, cdn(1:nacd(jwupi(js)), jwupi(js))) ! cb 1/16/13
                                write(WDO2(JFILE, 4), '(F10.3,*(F10.4,","))') JDAY, (CDAVGW(JWD, CDN(JD, JWUPI(JS))), JD = 1, NACD(JWUPI(JS)))
                            end if
                        else
!  JSSS(JBUPI(JS))=JSSS(JBUPI(JS))+1
                            qwdo(j) = qwdo(j) + qpi(js) ! cb 1/16/13
                            twdo(j) = twdo(j) + qpi(js)*tavg(jsss(jbupi(js)), jbupi(js))
                            if (QPI(JS) == 0.0) then
                                TAVG(JSSS(JBUPI(JS)), JBUPI(JS)) = -99.0
                                CAVG(JSSS(JBUPI(JS)), JBUPI(JS), :) = -99.0
                                CDAVG(JSSS(JBUPI(JS)), JBUPI(JS), :) = -99.0
                            end if
                            write(WDO2(JFILE, 2), '(F10.3,",",F8.2)') JDAY, TAVG(JSSS(JBUPI(JS)), JBUPI(JS))
                            if (CONSTITUENTS) then
                                cwdo(cn(1:nac), j) = cwdo(cn(1:nac), j) + qpi(js)*cavg(jsss(jbupi(js)), jbupi(js), cn(1:nac)) ! cb 1/16/13
                                write(WDO2(JFILE, 3), '(F10.3,",",*(F10.4,","))') JDAY, (CAVG(JSSS(JBUPI(JS)), JBUPI(JS), CN(JC)), JC = 1, NAC)
                            end if
                            if (DERIVED_CALC) then
                                cdwdo(cdn(1:nacd(jw), jw), j) = cdwdo(cdn(1:nacd(jw), jw), j) + qpi(js)*cdavg(jsss(jbupi(js)), jbupi(js), cdn(1:nacd(jwupi(js)), jwupi(js))) ! cb 1/16/13
                                write(WDO2(JFILE, 4), '(F10.3,",",*(F10.4,","))') JDAY, (CDAVG(JSSS(JBUPI(JS)), JBUPI(JS), CDN(JD, JWUPI(JS))), JD = 1, NACD(JWUPI(JS)))
                            end if
                        end if
                    end if
                end do
                do JS = 1, NGT ! gates
                    if (LATERAL_GATE(JS)) then
                        JWD = JWD + 1
                    else
                        JSSS(JBUGT(JS)) = JSSS(JBUGT(JS)) + 1
                    end if

                    if (IWDO(J) == IUGT(JS)) then
                        JFILE = JFILE + 1
                        write(WDO2(JFILE, 1), '(F10.3,",",F10.4)') JDAY, QGT(JS)
                        if (LATERAL_GATE(JS)) then
!    JWD=JWD+1
                            qwdo(j) = qwdo(j) + qgt(js) ! cb 1/16/13
                            twdo(j) = twdo(j) + qgt(js)*tavgw(jwd)
                            if (QGT(JS) == 0.0) then
                                TAVGW(JWD) = -99.0
                                CAVGW(JWD, :) = -99.0
                                CDAVGW(JWD, :) = -99.0
                            end if
                            write(WDO2(JFILE, 2), '(F10.3,",",F8.2)') JDAY, TAVGW(JWD)
                            if (CONSTITUENTS) then
                                cwdo(cn(1:nac), j) = cwdo(cn(1:nac), j) + qgt(js)*cavgw(jwd, cn(1:nac)) ! cb 1/16/13
                                write(WDO2(JFILE, 3), '(F10.3,",",*(F10.4,","))') JDAY, (CAVGW(JWD, CN(JC)), JC = 1, NAC)
                            end if
                            if (DERIVED_CALC) then
                                cdwdo(cdn(1:nacd(jw), jw), j) = cdwdo(cdn(1:nacd(jw), jw), j) + qgt(js)*cdavgw(jwd, cdn(1:nacd(jwugt(js)), jwugt(js))) ! cb 1/16/13
                                write(WDO2(JFILE, 4), '(F10.3,",",*(F10.4,","))') JDAY, (CDAVGW(JWD, CDN(JD, JWUGT(JS))), JD = 1, NACD(JWUGT(JS)))
                            end if
                        else
!   JSSS(JBUGT(JS))=JSSS(JBUGT(JS))+1
                            qwdo(j) = qwdo(j) + qgt(js) ! cb 1/16/13
                            twdo(j) = twdo(j) + qgt(js)*tavg(jsss(jbugt(js)), jbugt(js))
                            if (QGT(JS) == 0.0) then
                                TAVG(JSSS(JBUGT(JS)), JBUGT(JS)) = -99.0
                                CAVG(JSSS(JBUGT(JS)), JBUGT(JS), :) = -99.0
                                CDAVG(JSSS(JBUGT(JS)), JBUGT(JS), :) = -99.0
                            end if
                            write(WDO2(JFILE, 2), '(F10.3,",",F8.2)') JDAY, TAVG(JSSS(JBUGT(JS)), JBUGT(JS))
                            if (CONSTITUENTS) then
                                cwdo(cn(1:nac), j) = cwdo(cn(1:nac), j) + qgt(js)*cavg(jsss(jbugt(js)), jbugt(js), cn(1:nac)) ! cb 1/16/13
                                write(WDO2(JFILE, 3), '(F10.3,",",*(F10.4,","))') JDAY, (CAVG(JSSS(JBUGT(JS)), JBUGT(JS), CN(JC)), JC = 1, NAC)
                            end if
                            if (DERIVED_CALC) then
                                cdwdo(cdn(1:nacd(jw), jw), j) = cdwdo(cdn(1:nacd(jw), jw), j) + qgt(js)*cdavg(jsss(jbugt(js)), jbugt(js), cdn(1:nacd(jwugt(js)), jwugt(js))) ! cb 1/16/13
                                write(WDO2(JFILE, 4), '(F10.3,",",*(F10.4,","))') JDAY, (CDAVG(JSSS(JBUGT(JS)), JBUGT(JS), CDN(JD, JWUGT(JS))), JD = 1, NACD(JWUGT(JS)))
                            end if
                        end if
                    end if
                end do

! cb 1/16/13 deleted old withdrawal output code

                if (QWDO(J) /= 0.0) then
                    TWDO(J) = TWDO(J)/QWDO(J)
                end if
                do JC = 1, NAC
                    if (QWDO(J) /= 0.0) then
                        CWDO(CN(JC), J) = CWDO(CN(JC), J)/QWDO(J)
                    end if
                    write(CWDOC(CN(JC)), "(F10.4)") CWDO(CN(JC), J) ! SW 9/23/13 Changed format from G8.3 to F8.3 to avoid format overflow
                    CWDOC(CN(JC)) = ADJUSTR(CWDOC(CN(JC)))
                end do
                do JW = 1, NWB
                    if (IWDO(J) >= US(BS(JW)) .and. IWDO(J) <= DS(BE(JW))) then
                        exit
                    end if
                end do
                do JD = 1, NACD(JW)
                    if (QWDO(J) /= 0.0) then
                        CDWDO(CDN(JD, JW), J) = CDWDO(CDN(JD, JW), J)/QWDO(J)
                    end if
                    write(CDWDOC(CDN(JD, JW)), "(F10.4)") CDWDO(CDN(JD, JW), J) ! SW 9/23/13 Changed format from G8.3 to F8.3 to avoid format overflow
                    CDWDOC(CDN(JD, JW)) = ADJUSTR(CDWDOC(CDN(JD, JW)))
                end do
                2499 write(WDO(J, 1), '(F10.3,",",F9.3,",",8X,*(F9.3,","))', ERR=2500) JDAY, QWDO(J), (QOUTLET(I), I = 1, NUMOUTLETS) ! sw 3/2019 This code was necessary during multiple WB read/write possible error with one file copying and another writing at the same time
                go to 2501
                2500 call SLEEPQQ(100)
                write(9911, "(A,f12.3,A,f12.3,A,f12.3)") "ERROR writing qwo file output on JDAY", JDAY, " retrying read after pausing 0.1 s"
                go to 2499
                2501 continue
                write(WDO(J, 2), '(F10.3,",",F8.2,",",8X,*(F8.2,","))', ERR=2502) JDAY, TWDO(J), (TOUTLET(I), I = 1, NUMOUTLETS)
                go to 2503
                2502 call SLEEPQQ(100)
                write(9911, "(A,f12.3,A,f12.3,A,f12.3)") "ERROR writing two file output on JDAY", JDAY, " retrying read after pausing 0.1 s"
                go to 2501
                2503 continue
                if (CONSTITUENTS) then
                    write(WDO(J, 3), '(F10.3,",",*(A10,","))', ERR=2504) JDAY, (CWDOC(CN(JC)), JC = 1, NAC)
                end if
                go to 2505
                2504 call SLEEPQQ(100)
                write(9911, "(A,f12.3,A,f12.3,A,f12.3)") "ERROR writing cwo file output on JDAY", JDAY, " retrying read after pausing 0.1 s"
                go to 2503
                2505 continue
                if (DERIVED_CALC) then
                    write(WDO(J, 4), '(F10.3,",",*(A10,","))') JDAY, (CDWDOC(CDN(JD, JW)), JD = 1, NACD(JW))
                end if
            end do
        end if
    end if

!**** DSI W2 Linkage File (W2L) (Supercedes Old Velocity vectors)
    if (VECTOR(1)) then
! *** Apply the same linkage settings for all waterbodies
        if (JDAY >= NXTMVP(1) .or. JDAY >= VPLD(VPLDP(1) + 1, 1)) then
            if (JDAY >= VPLD(VPLDP(1) + 1, 1)) then
                VPLDP(1) = VPLDP(1) + 1
                NXTMVP(1) = VPLD(VPLDP(1), 1)
            end if
            NXTMVP(1) = NXTMVP(1) + VPLF(VPLDP(1), 1)

! *** Write the W2L Snapshot (DSI)
            write(VPL(1)) REAL(JDAY, 4)

! *** Compute the elevation, with ELWS zeroed for segments < CUS
            do JW = 1, NWB
                do JB = BS(JW), BE(JW)
                    do I = US(JB) - 1, DS(JB) + 1
                        if (I < CUS(JB)) then
                            WSEL(I) = -9999
                        else
                            WSEL(I) = ELWS(I)
                        end if
                    end do
                end do
            end do

            write(VPL(1)) (WSEL(I), I = 1, IMX)
            write(VPL(1)) ((REAL(U(K, I), 4), K = 1, KMX), I = 1, IMX)
            write(VPL(1)) ((REAL(W(K, I), 4), K = 1, KMX), I = 1, IMX)

            do I = 1, IMX
! *** ABOVE ACTIVE LAYER
                do K = 1, KTI(I) - 1
                    WDSI(K, I) = -9999.
                end do
! *** ACTIVE LAYERS
                do K = KTI(I), KB(I)
                    WDSI(K, I) = T2(K, I)
                end do
! *** BELOW ACTIVE LAYER
                do K = KB(I) + 1, KMX
                    WDSI(K, I) = -9999.
                end do
            end do
            write(VPL(1)) ((WDSI(K, I), K = 1, KMX), I = 1, IMX)

! *** Constituent data
            do JC = 1, NAC
                IC = CN(JC)
                do I = 1, IMX
! *** ABOVE ACTIVE LAYER
                    do K = 1, KTI(I) - 1
                        WDSI(K, I) = -9999.
                    end do
! *** ACTIVE LAYERS
                    do K = KTI(I), KB(I)
                        WDSI(K, I) = C2(K, I, IC)
                    end do
! *** BELOW ACTIVE LAYER
                    do K = KB(I) + 1, KMX
                        WDSI(K, I) = -9999.
                    end do
                end do
                write(VPL(1)) ((WDSI(K, I), K = 1, KMX), I = 1, IMX)
            end do

        end if
    end if

!** Restart

    if (RESTART_OUT) then
        if (JDAY >= NXTMRS .or. JDAY >= RSOD(RSODP + 1)) then
            if (JDAY >= RSOD(RSODP + 1)) then
                RSODP = RSODP + 1
                NXTMRS = RSOD(RSODP)
            end if
            NXTMRS = NXTMRS + RSOF(RSODP)
            if (RSOF(RSODP) >= 1.0) then
                write(EXT, "(I0)") INT(JDAY)
            else
                write(EXT, '(I0,"_",I2)') INT(JDAY), INT(100.*(JDAY - INT(JDAY))) ! Allows for writing out file names with JDAY FREQ less than 1 day
            end if
            EXT = ADJUSTL(EXT)
            L = LEN_TRIM(EXT)
            RSOFN = "rso" // EXT(1:L) // ".opt"
            call RESTART_OUTPUT(RSOFN)
        end if
    end if

! Snapshot formats

    10490 format('CE-QUAL-W2 VERSION',F4.2/                                                                                          (1X,A72))
    10500 format(/1X,A/                                                                                                               3X,'Gregorian date      [GDAY] =',A19,1X,I0,', ',I0/                                                                 3X,'Julian date         [JDAY] =',I10,' days',F6.2,' hours'/                                                         3X,'Elapsed time      [ELTMJD] =',I10,' days',F6.2,' hours'/                                                         3X,'Timestep             [DLT] =',I10,' sec'/                                                                        3X,'  at location  [KLOC,ILOC] = (',I0,',',I0,')'/                                                                   3X,'Minimum timestep  [MINDLT] =',I10,' sec '/                                                                       3X,'  at Julian day    [JDMIN] =',I10,' days',F6.2,' hours'/                                                         3X,'  at location  [KMIN,IMIN] = (',I0,',',I0,')')
    10510 format(3X,'Limiting timestep'/                                                                                              3X,'  at location  [KMIN,IMIN] = (',I0,',',I0,')')
    10520 format(3X,'Average timestep   [DLTAV] =',I10,' sec'/                                                                        3X,'Number of iterations [NIT] =',I10/                                                                               3X,'Number of violations  [NV] =',I10/)
    10530 format(1X,A)
    10540 format(3X,'Input'/                                                                                                          3X,'  Air temperature          [TAIR] =',F9.2,1X,A/                                                                  3X,'  Dewpoint temperature     [TDEW] =',F9.2,1X,A/                                                                  3X,'  Wind direction            [PHI] =',F9.2,' rad'/                                                                3X,'  Cloud cover             [CLOUD] =',F9.2/                                                                       3X,'  Calculated'/                                                                                                   5X,'  Equilibrium temperature    [ET] =',F9.2,1X,A/                                                                  5X,'  Surface heat exchange    [CSHE] =',E9.2,' m/sec'/                                                              5X,'  Net short wave radiation [SRON] =',E9.2,1X,A,' W/m^2'/)
    10550 format(1X,A/                                                                                                                3X,A)
    10560 format(5X,'Branch ',I0/                                                                                                     5X,'  Layer       [KQIN] = ',I0,'-',I0/                                                                              5X,'  Inflow       [QIN] =',F8.2,' m^3/sec'/                                                                         5X,'  Temperature  [TIN] =',F8.2,1X,A)
    10570 format(/3X,'Distributed Tributaries')
    10580 format(5X,'Branch ',I0/                                                                                                     5X,'  Inflow      [QDTR] =',F8.2,' m^3/sec'/                                                                         5X,'  Temperature [TDTR] =',F8.2,1X,A)
    10590 format(:/3X,'Tributaries'/                                                                                                  5X,'Segment     [ITR] =',11I8:/                                                                                      (T25,11I8))
    10600 format(:5X,'Layer      [KTWB] = ',11(I0,'-',I0,2X):/                                                                        (T25,11(I0,'-',I0)))
    10610 format(:5X,'Inflow      [QTR] =',11F8.2:/                                                                                   (T25,11F8.1))
    10620 format(:5X,'Temperature [TTR] =',11F8.2:/                                                                                   (T25,11F8.1))
    10630 format(/1X,'Outflows')
    10640 format(3X,'Structure outflows [QSTR]'/                                                                                      3X,'  Branch ',I0,' = ',11F8.2:/                                                                                     (T16,11F8.2))
    10650 format(:/3X,'Total outflow [QOUT] =',F8.2,' m^3/s'/                                                                         5X,'Outlets'/                                                                                                        5X,'  Layer             [KOUT] =',12I7:/                                                                             (33X,12I7))
    10660 format(:7X,'Outflow (m^3/sec) [QOUT] =',12F7.2:/                                                                            (33X,12F7.2))
    10670 format(:5X,'Withdrawals'/                                                                                                   5X,'  Segment            [IWD] =',I7/                                                                                5X,'  Outflow (m^3/sec)  [QWD] =',F7.2)
    10680 format(5X,'  Layer              [KWD] =',12I7/                                                                              (33X,12I7))
    10690 format(:5X,'  Outflow (m^3/sec)  [QSW] =',12F7.2/                                                                           (33X,12F7.2))
    10700 format(/'1',A)
    10710 format(3X,'Branch ',I0,' [CIN]'/                                                                                            (5X,A,T25,'=',F9.3,1X,A))
    10720 format(3X,'Tributary ',I0,' [CTR]'/                                                                                         (5X,A,T25,'=',F9.3,1X,A))
    10730 format(3X,'Distributed tributary ',I0,' [CDT]'/                                                                             (5X,A,T25,'=',F9.3,1X,A))
    10740 format(/'Surface calculations')
    10750 format(3X,'Evaporation rate [EV]'/                                                                                          (:3X,'  Branch ',I0,' = ',E10.3,' m^3/s')) ! SW 9/15/05 4/21/10
    10755 format(3x,'Cumulative evaporation [VOLEV]'/                                                                                 (:3X,'  Branch ',I0,' = ',F0.1,' m^3'))
    10760 format(3X,'Precipitation [PR]'/                                                                                             (3X,'  Branch ',I0,' = ',F8.6),' m/s')
    10770 format(/1X,'External head boundary elevations'/)
    10780 format(3X,'Branch ',I0/5X,'Upstream elevation   [ELUH] =',F8.3,' m')
    10790 format(3X,'Branch ',I0/5X,'Downstream elevation [ELDH] =',F8.3,' m')
    10800 format(/'Water Balance')
    10810 format(3X,'Waterbody ',I0/                                                                                                  3X,'  Spatial change  [VOLSR]  = ',E15.8,' m^3'/                                                                     3X,'  Temporal change [VOLTR]  = ',E15.8,' m^3'/                                                                     3X,'  Volume error             = ',E15.8,' m^3'/                                                                     3X,'  Percent error            = ',E15.8,' %')
    10820 format(3X,'Branch ',I0/                                                                                                     3X,'  Spatial change  [VOLSBR] = ',E15.8,' m^3'/                                                                     3X,'  Temporal change [VOLTBR] = ',E15.8,' m^3'/                                                                     3X,'  Volume error             = ',E15.8,' m^3'/                                                                     3X,'  Percent error            = ',E15.8,' %')
    10830 format(/1X,'Energy Balance')
    10840 format(3X,'Waterbody ',I0/                                                                                                  3X,'  Spatially integrated energy   [ESR] = ',E15.8,' kJ'/                                                           3X,'  Temporally integrated energy  [ETR] = ',E15.8,' kJ'/                                                           3X,'  Energy error                        = ',E15.8,' kJ'/                                                           3X,'  Percent error                       = ',E15.8,' %')
    10850 format(3X,'  Spatially integrated energy  [ESBR] = ',E15.8,' kJ'/                                                           3X,'  Temporally integrated energy [ETBR] = ',E15.8,' kJ'/                                                           3X,'  Energy error                        = ',E15.8,' kJ'/                                                           3X,'  Percent error                       = ',E15.8,' %')
    10860 format(/1X,'Mass Balance')
    10870 format(3X,'Branch ',I0)
    10880 format(5X,A/                                                                                                                5X,'  Spatially integrated mass  [CMBRS] = ',E15.8,1X,A/                                                             5X,'  Temporally integrated mass [CMBRT] = ',E15.8,1X,A/                                                             5X,'  Mass error                         = ',E15.8,1X,A/                                                             5X,'  Percent error                      = ',E15.8,' %')
    10890 format(/1X,A/                                                                                                               3X,'Surface layer [KT] = ',I0/                                                                                       3X,'Elevation   [ELKT] =',F10.3,' m')
    10900 format(/3X,'Current upstream segment [CUS]'/                                                                                (3X,'  Branch ',I0,' = ',I0))

    return

end subroutine OUTPUTA
