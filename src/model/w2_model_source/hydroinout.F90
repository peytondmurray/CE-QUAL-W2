subroutine HYDROINOUT()

    use MAIN
    use GLOBAL;     use NAMESC;     use GEOMC;     use LOGICC;     use PREC;     use SURFHE;     use KINETIC;     use SHADEC;     use EDDY
    use STRUCTURES;     use TRANS;     use TVDC;     use SELWC;     use GDAYC;     use SCREENC;     use TDGAS;     use RSTART
    use MACROPHYTEC;     use POROSITYC;     use ZOOPLANKTONC
    use modSYSTDG, only: GTNAME, SYSTDG_TDG, UPDATE_TDGC, TDG_TDG, TDG_ROSP, POWNO, POWGTNO, FLNO, FLGTNO, TDGLOC, ip, ifl ! systdg
    implicit none
    external :: RESTART_OUTPUT
    integer :: JBU, JBD, JLAT, JWU
    real(R8) :: ELW, CGAS, TM, VPTG, DTVL, RHOTR, VQTR, VQTRI, QTRFR, AKBR, FW
    real(R8) :: TSUM, QSUMM
    real(R8) :: dosat, n2sat ! systdg - update powerhouse release tdg
!***********************************************************************************************************************************
!**                                            Task 2.1: Hydrodynamic sources/sinks                                               **
!***********************************************************************************************************************************

    QINSUM = 0.0;     TINSUM = 0.0;     CINSUM = 0.0;     UXBR = 0.0;     UYBR = 0.0;     tdgon = .false. ! cb 1/16/13
    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            if (BR_INACTIVE(JB)) then
                cycle
            end if
            IU = CUS(JB)
            ID = DS(JB)
            TSUM = 0.0;             CSUM = 0.0;             QSUM(JB) = 0.0;             QOUT(:, JB) = 0.0;             TOUT(JB) = 0.0;             COUT(:, JB) = 0.0

!****** Densities

            do I = IU - 1, ID + 1
                do K = KT, KB(I)
                    TISS(K, I) = 0.0
                    do JS = 1, NSS
                        TISS(K, I) = TISS(K, I) + SS(K, I, JS)
                    end do
                    RHO(K, I) = DENSITY(T2(K, I), DMAX1(TDS(K, I), 0.0D0), DMAX1(TISS(K, I), 0.0D0))
                end do
            end do
! v3.5 deleted pumpback code from v3.2
            do JS = 1, NSTR(JB)
                if (QSTR(JS, JB) /= 0.0) then
                    call DOWNSTREAM_WITHDRAWAL(JS)
                end if
            end do
            do K = KT, KB(ID)
                QSUM(JB) = QSUM(JB) + QOUT(K, JB)
                TSUM = TSUM + QOUT(K, JB)*T2(K, ID)
                CSUM(CN(1:NAC)) = CSUM(CN(1:NAC)) + QOUT(K, JB)*C2(K, ID, CN(1:NAC))
            end do
            if (QSUM(JB) /= 0.0) then
                TOUT(JB) = TSUM/QSUM(JB)
                COUT(CN(1:NAC), JB) = CSUM(CN(1:NAC))/QSUM(JB)
            end if
            if (QSUM(JB) /= 0.0 .and. DAM_OUTFLOW(JB)) then !TC 08/03/04
                TINSUM(JBDAM(JB)) = (TSUM + QINSUM(JBDAM(JB))*TINSUM(JBDAM(JB)))/(QSUM(JB) + QINSUM(JBDAM(JB)))
                CINSUM(CN(1:NAC), JBDAM(JB)) = (CSUM(CN(1:NAC)) + QINSUM(JBDAM(JB))*CINSUM(CN(1:NAC), JBDAM(JB)))/(QSUM(JB) + QINSUM(JBDAM(JB)))
                QINSUM(JBDAM(JB)) = QINSUM(JBDAM(JB)) + QSUM(JB)
            end if
        end do
    end do
    ILAT = 0
    JWW = NWD
    withdrawals = jww > 0
    if (nwdt > nwd) then
        qwd(nwd + 1:nwdt) = 0.0
    end if ! SW 10/30/2017
    JTT = NTR
    tributaries = jtt > 0
    if (ntrt > ntr) then
        qtr(ntr + 1:ntrt) = 0.0
    end if ! SW 10/30/2017
    JSS = NSTR
    if (SPILLWAY) then
        call SPILLWAY_FLOW()
        do JS = 1, NSP

!****** Positive flows

            JLAT = 0
            JBU = JBUSP(JS)
            JBD = JBDSP(JS)
            tdgon = .false. ! cb 1/16/13
            jsg = js
            nnsg = 0
            if (cac(ndo) == "      ON" .and. gasspc(js) == "      ON") then
                tdgon = .true.
            end if
            if (QSP(JS) >= 0.0) then
                if (LATERAL_SPILLWAY(JS)) then
                    JWW = JWW + 1
                    IWD(JWW) = IUSP(JS)
                    QWD(JWW) = QSP(JS)
                    KTWD(JWW) = KTUSP(JS)
                    KBWD(JWW) = KBUSP(JS)
                    EWD(JWW) = ESP(JS)
                    JBWD(JWW) = JBU
                    I = MAX(CUS(JBWD(JWW)), IWD(JWW))
                    JB = JBWD(JWW)
                    JW = JWUSP(JS)
                    KT = KTWB(JW)
                    jwd = jww
                    call LATERAL_WITHDRAWAL() !(JWW)
                    do K = KTW(JWW), KBW(JWW)
                        QSS(K, I) = QSS(K, I) - QSW(K, JWW)
                    end do
                    if (IDSP(JS) /= 0) then ! cb 9/11/13
                        JTT = JTT + 1
                        QTR(JTT) = QSP(JS)
                        ITR(JTT) = IDSP(JS)
                        PLACE_QTR(JTT) = PDSPC(JS) == " DENSITY"
                        SPECIFY_QTR(JTT) = PDSPC(JS) == " SPECIFY"
                        if (SPECIFY_QTR(JTT)) then
                            ELTRT(JTT) = ETDSP(JS)
                            ELTRB(JTT) = EBDSP(JS)
                        end if
                        JBTR(JTT) = JBD
                    end if ! cb 9/11/13
                    if (IDSP(JS) /= 0 .and. QSP(JS) > 0.0) then
                        TSUM = 0.0;                         QSUMM = 0.0;                         CSUM = 0.0
                        do K = KTW(JWW), KBW(JWW)
                            QSUMM = QSUMM + QSW(K, JWW)
                            TSUM = TSUM + QSW(K, JWW)*T2(K, IWD(JWW))
                            CSUM(CN(1:NAC)) = CSUM(CN(1:NAC)) + QSW(K, JWW)*C2(K, IWD(JWW), CN(1:NAC))
                        end do
                        TTR(JTT) = TSUM/QSUMM
                        do JC = 1, NAC
                            CTR(CN(JC), JTT) = CSUM(CN(JC))/QSUMM
                            if (CN(JC) == NDO .and. GASSPC(JS) == "      ON" .and. QSP(JS) > 0.0) then
                                TDG_SPILLWAY(JWW, JS) = .true.
                                call TOTAL_DISSOLVED_GAS(0, PALT(I), 0, JS, TTR(JTT), CTR(CN(JC), JTT)) ! DO
                            end if
                            if (CN(JC) == NN2 .and. GASSPC(JS) == "      ON" .and. QSP(JS) > 0.0) then
                                TDG_SPILLWAY(JWW, JS) = .true.
                                call TOTAL_DISSOLVED_GAS(1, PALT(I), 0, JS, TTR(JTT), CTR(CN(JC), JTT)) ! N2
                            end if
                            if (CN(JC) == NDGP .and. GASSPC(JS) == "      ON" .and. QSP(JS) > 0.0) then !8/2020 TDGP
                                TDG_SPILLWAY(JWW, JS) = .true.
                                call TOTAL_DISSOLVED_GAS(2, PALT(I), 0, JS, TTR(JTT), CTR(CN(JC), JTT))
                            end if
                        end do
                    else
                        if (CAC(NDO) == "      ON" .and. GASSPC(JS) == "      ON" .and. QSP(JS) > 0.0) then
                            TDG_SPILLWAY(JWW, JS) = .true.
                        end if
                    end if

                else
                    JSS(JBU) = JSS(JBU) + 1
                    KTSW(JSS(JBU), JBU) = KTUSP(JS)
                    KBSW(JSS(JBU), JBU) = KBUSP(JS)
                    JB = JBU
                    POINT_SINK(JSS(JBU), JBU) = .true.
                    ID = IUSP(JS)
                    QSTR(JSS(JBU), JBU) = QSP(JS)
                    ESTR(JSS(JBU), JBU) = ESP(JS)
                    KT = KTWB(JWUSP(JS))
                    JW = JWUSP(JS)
                    call DOWNSTREAM_WITHDRAWAL(JSS(JBU))
                    QSUM(JB) = 0.0;                     TSUM = 0.0;                     CSUM = 0.0
                    do K = KT, KB(ID)
                        QSUM(JB) = QSUM(JB) + QOUT(K, JB)
                        TSUM = TSUM + QOUT(K, JB)*T2(K, ID)
                        do JC = 1, NAC
                            if (CN(JC) == NDO .and. CAC(NDO) == "      ON" .and. GASSPC(JS) == "      ON" .and. QSP(JS) > 0.0) then ! MM 5/21/2009
                                T2R4 = T2(K, ID)
                                CGAS = C2(K, ID, CN(JC)) ! MM 5/21/2009
                                call TOTAL_DISSOLVED_GAS(0, PALT(ID), 0, JS, T2R4, CGAS) ! O2
                                CSUM(CN(JC)) = CSUM(CN(JC)) + QOUT(K, JB)*CGAS
!ELSEIF (CN(JC)==NGN2 .AND. CAC(NGN2) == '      ON' .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN     ! SW 10/27/15
                            else
                                if (CN(JC) == NN2 .and. GASSPC(JS) == "      ON" .and. QSP(JS) > 0.0) then ! SW 10/27/15
                                    if (CAC(NN2) == "      ON") then ! cb 1/13/16
                                        T2R4 = T2(K, ID)
                                        CGAS = C2(K, ID, CN(JC)) ! 
                                        call TOTAL_DISSOLVED_GAS(1, PALT(ID), 0, JS, T2R4, CGAS) ! N2
                                        CSUM(CN(JC)) = CSUM(CN(JC)) + QOUT(K, JB)*CGAS
                                    end if
                                else
                                    if (CN(JC) == NDGP .and. GASSPC(JS) == "      ON" .and. QSP(JS) > 0.0) then !8/2020 TDGP
                                        if (CAC(NDGP) == "      ON") then
                                            T2R4 = T2(K, ID)
                                            CGAS = C2(K, ID, CN(JC))
                                            call TOTAL_DISSOLVED_GAS(2, PALT(ID), 0, JS, T2R4, CGAS)
                                            CSUM(CN(JC)) = CSUM(CN(JC)) + QOUT(K, JB)*CGAS
                                        end if
                                    else
                                        CSUM(CN(JC)) = CSUM(CN(JC)) + QOUT(K, JB)*C2(K, ID, CN(JC))
                                    end if
                                end if
                            end if
                        end do
                    end do
                    if (QSUM(JB) /= 0.0) then
                        TOUT(JB) = TSUM/QSUM(JB)
                        COUT(CN(1:NAC), JB) = CSUM(CN(1:NAC))/QSUM(JB)
                    end if
                    if (IDSP(JS) /= 0 .and. US(JBD) == IDSP(JS)) then
                        QSUMM = 0.0;                         TSUM = 0.0;                         CSUM = 0.0
                        do K = KT, KB(ID)
                            QSUMM = QSUMM + QNEW(K)
                            TSUM = TSUM + QNEW(K)*T2(K, ID)
                            do JC = 1, NAC
                                if (CN(JC) == NDO .and. CAC(NDO) == "      ON" .and. GASSPC(JS) == "      ON" .and. QSP(JS) > 0.0) then ! MM 5/21/2009
                                    T2R4 = T2(K, ID)
                                    CGAS = C2(K, ID, CN(JC)) ! MM 5/21/2009
                                    call TOTAL_DISSOLVED_GAS(0, PALT(ID), 0, JS, T2R4, CGAS)
                                    CSUM(CN(JC)) = CSUM(CN(JC)) + QNEW(K)*CGAS
!ELSEIF (CN(JC)==NGN2 .AND. CAC(NGN2) == '      ON' .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN             ! SW 10/27/15
                                else
                                    if (CN(JC) == NN2 .and. GASSPC(JS) == "      ON" .and. QSP(JS) > 0.0) then ! SW 10/27/15
                                        if (CAC(NN2) == "      ON") then ! cb 1/13/16
                                            T2R4 = T2(K, ID)
                                            CGAS = C2(K, ID, CN(JC)) ! 
                                            call TOTAL_DISSOLVED_GAS(1, PALT(ID), 0, JS, T2R4, CGAS)
                                            CSUM(CN(JC)) = CSUM(CN(JC)) + QNEW(K)*CGAS
                                        end if
                                    else
                                        if (CN(JC) == NDGP .and. GASSPC(JS) == "      ON" .and. QSP(JS) > 0.0) then !8/2020 TDGP        
                                            if (CAC(NDGP) == "      ON") then
                                                T2R4 = T2(K, ID)
                                                CGAS = C2(K, ID, CN(JC))
                                                call TOTAL_DISSOLVED_GAS(2, PALT(ID), 0, JS, T2R4, CGAS)
                                                CSUM(CN(JC)) = CSUM(CN(JC)) + QNEW(K)*CGAS
                                            end if
                                        else
                                            CSUM(CN(JC)) = CSUM(CN(JC)) + QNEW(K)*C2(K, ID, CN(JC))
                                        end if
                                    end if
                                end if
                            end do
                        end do
                        if (QSUMM /= 0.0) then
                            TINSUM(JBD) = (TSUM + QINSUM(JBD)*TINSUM(JBD))/(QSUMM + QINSUM(JBD))
                            CINSUM(CN(1:NAC), JBD) = (CSUM(CN(1:NAC)) + QINSUM(JBD)*CINSUM(CN(1:NAC), JBD))/(QSUMM + QINSUM(JBD))
                            QINSUM(JBD) = QINSUM(JBD) + QSUMM
                        end if
                    else
                        if (IDSP(JS) /= 0) then
                            JTT = JTT + 1
                            QTR(JTT) = QSP(JS)
                            ITR(JTT) = IDSP(JS)
                            PLACE_QTR(JTT) = PDSPC(JS) == " DENSITY"
                            SPECIFY_QTR(JTT) = PDSPC(JS) == " SPECIFY"
                            if (SPECIFY_QTR(JTT)) then
                                ELTRT(JTT) = ETDSP(JS)
                                ELTRB(JTT) = EBDSP(JS)
                            end if
                            JBTR(JTT) = JBD
                            TTR(JTT) = TOUT(JB)
                            do JC = 1, NAC
                                CTR(CN(JC), JTT) = COUT(CN(JC), JB)
                                if (CN(JC) == NDO .and. GASSPC(JS) == "      ON" .and. QSP(JS) > 0.0) then
                                    call TOTAL_DISSOLVED_GAS(0, PALT(ITR(JTT)), 0, JS, TTR(JTT), CTR(CN(JC), JTT))
                                end if
                                if (CN(JC) == NN2 .and. GASSPC(JS) == "      ON" .and. QSP(JS) > 0.0) then ! SW 10/27/15
                                    call TOTAL_DISSOLVED_GAS(1, PALT(ITR(JTT)), 0, JS, TTR(JTT), CTR(CN(JC), JTT))
                                end if
                                if (CN(JC) == NDGP .and. GASSPC(JS) == "      ON" .and. QSP(JS) > 0.0) then
                                    call TOTAL_DISSOLVED_GAS(2, PALT(ITR(JTT)), 0, JS, TTR(JTT), CTR(CN(JC), JTT))
                                end if
                            end do
                        end if
                    end if
                end if
            else
                if (QSP(JS) < 0.0) then
                    JTT = JTT + 1
                    JWW = JWW + 1
                    IWD(JWW) = IDSP(JS)
                    ITR(JTT) = IUSP(JS)
                    QTR(JTT) = -QSP(JS)
                    QWD(JWW) = -QSP(JS)
                    KTWD(JWW) = KTDSP(JS)
                    KBWD(JWW) = KBDSP(JS)
                    EWD(JWW) = ESP(JS)
                    PLACE_QTR(JTT) = PUSPC(JS) == " DENSITY"
                    SPECIFY_QTR(JTT) = PUSPC(JS) == " SPECIFY"
                    if (SPECIFY_QTR(JTT)) then
                        ELTRT(JTT) = ETUSP(JS)
                        ELTRB(JTT) = EBUSP(JS)
                    end if
                    JBTR(JTT) = JBU
                    JBWD(JWW) = JBD
                    I = MAX(CUS(JBWD(JWW)), IWD(JWW))
                    JB = JBWD(JWW)
                    JW = JWDSP(JS)
                    KT = KTWB(JW)
                    jwd = jww
                    call LATERAL_WITHDRAWAL() !(JWW)
                    do K = KTW(JWW), KBW(JWW)
                        QSS(K, I) = QSS(K, I) - QSW(K, JWW)
                    end do
                    if (IDSP(JS) /= 0) then
                        TSUM = 0.0;                         QSUMM = 0.0;                         CSUM = 0.0
                        do K = KTW(JWW), KBW(JWW)
                            QSUMM = QSUMM + QSW(K, JWW)
                            TSUM = TSUM + QSW(K, JWW)*T2(K, IWD(JWW))
                            CSUM(CN(1:NAC)) = CSUM(CN(1:NAC)) + QSW(K, JWW)*C2(K, IWD(JWW), CN(1:NAC))
                        end do
                        TTR(JTT) = TSUM/QSUMM
                        do JC = 1, NAC
                            CTR(CN(JC), JTT) = CSUM(CN(JC))/QSUMM
                            if (CN(JC) == NDO .and. GASSPC(JS) == "      ON" .and. QSP(JS) > 0.0) then
                                TDG_SPILLWAY(JWW, JS) = .true.
                                call TOTAL_DISSOLVED_GAS(0, PALT(I), 0, JS, TTR(JTT), CTR(CN(JC), JTT)) ! O2
                            end if
                            if (CN(JC) == NN2 .and. GASSPC(JS) == "      ON" .and. QSP(JS) > 0.0) then
                                TDG_SPILLWAY(JWW, JS) = .true.
                                call TOTAL_DISSOLVED_GAS(1, PALT(I), 0, JS, TTR(JTT), CTR(CN(JC), JTT)) ! N2
                            end if
                            if (CN(JC) == NDGP .and. GASSPC(JS) == "      ON" .and. QSP(JS) > 0.0) then !8/2020 TDGP
                                TDG_SPILLWAY(JWW, JS) = .true.
                                call TOTAL_DISSOLVED_GAS(2, PALT(I), 0, JS, TTR(JTT), CTR(CN(JC), JTT))
                            end if
                        end do
                    else
                        if (CAC(NDO) == "      ON" .and. GASSPC(JS) == "      ON" .and. QSP(JS) > 0.0) then
                            TDG_SPILLWAY(JWW, JS) = .true.
                        end if
                    end if
                end if
            end if
        end do
    end if
    if (PUMPS) then
        do JP = 1, NPU
            JLAT = 0
            JWU = JWUPU(JP)
            JBU = JBUPU(JP)
            JBD = JBDPU(JP)
            tdgon = .false. ! cb 1/16/13
            if (LATERAL_PUMP(JP)) then
                if (PUMP_DOWNSTREAM(JP)) then
                    ELW = EL(KTWB(JWDPU(JP)), IDPU(JP)) - Z(IDPU(JP))*COSA(JBD)
                    JWW = JWW + 1 ! SW 10/30/2017
                    JBWD(JWW) = JBU
                    IWD(JWW) = IUPU(JP)
                else
                    ELW = EL(KTWB(JWU), IUPU(JP)) - Z(IUPU(JP))*COSA(JBU)
                    JWW = JWW + 1 ! SW 10/30/2017
                    JBWD(JWW) = JBU
                    IWD(JWW) = IUPU(JP)
                end if
            else
                if (PUMP_DOWNSTREAM(JP)) then
                    ELW = EL(KTWB(JWDPU(JP)), IDPU(JP)) - Z(IDPU(JP))*COSA(JBD) - SINA(JBD)*DLX(IDPU(JP))*0.5
                    JSS(JBU) = JSS(JBU) + 1 ! SW 10/30/2017
                else
                    ELW = EL(KTWB(JWU), IUPU(JP)) - Z(IUPU(JP))*COSA(JBU) - SINA(JBU)*DLX(IUPU(JP))*0.5
                    JSS(JBU) = JSS(JBU) + 1 ! SW 10/30/2017
                end if
            end if

            if (JDAY >= ENDPU(JP)) then
                PUMPON(JP) = .false.
            end if !  CB 1/13/06
            if (JDAY >= STRTPU(JP) .and. JDAY < ENDPU(JP)) then
                if (PUMP_DOWNSTREAM(JP)) then ! IF BASED ON DOWNSTREAM WATER LEVEL AND NOT UPSTREAM
                    if (ELW >= EOFFPU(JP)) then
                        PUMPON(JP) = .false.
                    end if ! CB 1/13/06
                    if (ELW < EOFFPU(JP) .and. QPU(JP) > 0.0) then
                        if (ELW <= EONPU(JP)) then
                            PUMPON(JP) = .true.
                        end if
                        if (PUMPON(JP)) then
                            if (LATERAL_PUMP(JP)) then
                                JLAT = 1
                                QWD(JWW) = QPU(JP)
                                KTWD(JWW) = KTPU(JP)
                                KBWD(JWW) = KBPU(JP)
                                EWD(JWW) = EPU(JP)
                                I = MAX(CUS(JBWD(JWW)), IWD(JWW))
                                JB = JBWD(JWW)
                                JW = JWU
                                KT = KTWB(JW)
                                jwd = jww
                                call LATERAL_WITHDRAWAL()
                                do K = KTW(JWW), KBW(JWW)
                                    QSS(K, I) = QSS(K, I) - QSW(K, JWW)
                                end do
                                if (IDPU(JP) /= 0) then
                                    JTT = JTT + 1
                                    QTR(JTT) = QPU(JP)
                                    ITR(JTT) = IDPU(JP)
                                    PLACE_QTR(JTT) = PPUC(JP) == " DENSITY"
                                    SPECIFY_QTR(JTT) = PPUC(JP) == " SPECIFY"
                                    if (SPECIFY_QTR(JTT)) then
                                        ELTRT(JTT) = ETPU(JP)
                                        ELTRB(JTT) = EBPU(JP)
                                    end if
                                    JBTR(JTT) = JBD
                                    TSUM = 0.0;                                     QSUMM = 0.0;                                     CSUM(CN(1:NAC)) = 0.0
                                    do K = KTW(JWW), KBW(JWW)
                                        QSUMM = QSUMM + QSW(K, JWW)
                                        TSUM = TSUM + QSW(K, JWW)*T2(K, IWD(JWW))
                                        CSUM(CN(1:NAC)) = CSUM(CN(1:NAC)) + QSW(K, JWW)*C2(K, IWD(JWW), CN(1:NAC))
                                    end do
                                    if (QSUMM > 0.0) then
                                        TTR(JTT) = TSUM/QSUMM
                                        CTR(CN(1:NAC), JTT) = CSUM(CN(1:NAC))/QSUMM
                                    end if
                                end if
                            else
!JSS(JBU)                 =  JSS(JBU)+1     ! SW 9/25/13
                                KTSW(JSS(JBU), JBU) = KTPU(JP)
                                KBSW(JSS(JBU), JBU) = KBPU(JP)
                                JB = JBU
                                POINT_SINK(JSS(JBU), JBU) = .true.
                                ID = IUPU(JP)
                                QSTR(JSS(JBU), JBU) = QPU(JP)
                                ESTR(JSS(JBU), JBU) = EPU(JP)
                                KT = KTWB(JWU)
                                JW = JWU
                                call DOWNSTREAM_WITHDRAWAL(JSS(JBU))
                                if (IDPU(JP) /= 0 .and. US(JBD) == IDPU(JP)) then
                                    QSUMM = 0.0;                                     TSUM = 0.0;                                     CSUM = 0.0
                                    do K = KT, KB(ID)
                                        QSUMM = QSUMM + QNEW(K)
                                        TSUM = TSUM + QNEW(K)*T2(K, ID)
                                        CSUM(CN(1:NAC)) = CSUM(CN(1:NAC)) + QNEW(K)*C2(K, ID, CN(1:NAC))
                                    end do
                                    if (QSUMM /= 0.0) then
                                        TINSUM(JBD) = (TSUM + TINSUM(JBD)*QINSUM(JBD))/(QSUMM + QINSUM(JBD))
                                        CINSUM(CN(1:NAC), JBD) = (CSUM(CN(1:NAC)) + CINSUM(CN(1:NAC), JBD)*QINSUM(JBD))/(QSUMM + QINSUM(JBD))
                                        QINSUM(JBD) = QINSUM(JBD) + QSUMM
                                    end if
                                end if
                                QSUM(JB) = 0.0;                                 TSUM = 0.0;                                 CSUM = 0.0
                                do K = KT, KB(ID)
                                    QSUM(JB) = QSUM(JB) + QOUT(K, JB)
                                    TSUM = TSUM + QOUT(K, JB)*T2(K, ID)
                                    CSUM(CN(1:NAC)) = CSUM(CN(1:NAC)) + QOUT(K, JB)*C2(K, ID, CN(1:NAC))
                                end do
                                if (QSUM(JB) /= 0.0) then
                                    TOUT(JB) = TSUM/QSUM(JB)
                                    COUT(CN(1:NAC), JB) = CSUM(CN(1:NAC))/QSUM(JB)
                                end if
                                if (IDPU(JP) /= 0) then ! SW 9/25/13 Moved code start
                                    if (US(JBD) /= IDPU(JP) .or. HEAD_FLOW(JBD) .or. UP_HEAD(JBD)) then
                                        JTT = JTT + 1
                                        QTR(JTT) = QPU(JP)
                                        ITR(JTT) = IDPU(JP)
                                        PLACE_QTR(JTT) = PPUC(JP) == " DENSITY"
                                        SPECIFY_QTR(JTT) = PPUC(JP) == " SPECIFY"
                                        if (SPECIFY_QTR(JTT)) then
                                            ELTRT(JTT) = ETPU(JP)
                                            ELTRB(JTT) = EBPU(JP)
                                        end if
                                        JBTR(JTT) = JBD
                                        TTR(JTT) = TOUT(JB)
                                        CTR(CN(1:NAC), JTT) = COUT(CN(1:NAC), JB)
                                    end if
                                end if
                            end if
                        end if
                    end if



                else
                    if (ELW <= EOFFPU(JP)) then
                        PUMPON(JP) = .false.
                    end if ! CB 1/13/06
                    if (ELW > EOFFPU(JP) .and. QPU(JP) > 0.0) then
                        if (ELW >= EONPU(JP)) then
                            PUMPON(JP) = .true.
                        end if
                        if (PUMPON(JP)) then
                            if (LATERAL_PUMP(JP)) then
                                JLAT = 1
!JWW       = JWW+1               ! SW 9/25/13
!JBWD(JWW) = JBU  
!IWD(JWW)  = IUPU(JP)
                                QWD(JWW) = QPU(JP)
                                KTWD(JWW) = KTPU(JP)
                                KBWD(JWW) = KBPU(JP)
                                EWD(JWW) = EPU(JP)
                                I = MAX(CUS(JBWD(JWW)), IWD(JWW))
                                JB = JBWD(JWW)
                                JW = JWU
                                KT = KTWB(JW)
                                jwd = jww
                                call LATERAL_WITHDRAWAL() ! (JWW)
                                do K = KTW(JWW), KBW(JWW)
                                    QSS(K, I) = QSS(K, I) - QSW(K, JWW)
                                end do
                                if (IDPU(JP) /= 0) then ! MOVED CODE SW 9/25/13
                                    JTT = JTT + 1
                                    QTR(JTT) = QPU(JP)
                                    ITR(JTT) = IDPU(JP)
                                    PLACE_QTR(JTT) = PPUC(JP) == " DENSITY"
                                    SPECIFY_QTR(JTT) = PPUC(JP) == " SPECIFY"
                                    if (SPECIFY_QTR(JTT)) then
                                        ELTRT(JTT) = ETPU(JP)
                                        ELTRB(JTT) = EBPU(JP)
                                    end if
                                    JBTR(JTT) = JBD
                                    TSUM = 0.0;                                     QSUMM = 0.0;                                     CSUM(CN(1:NAC)) = 0.0
                                    do K = KTW(JWW), KBW(JWW)
                                        QSUMM = QSUMM + QSW(K, JWW)
                                        TSUM = TSUM + QSW(K, JWW)*T2(K, IWD(JWW))
                                        CSUM(CN(1:NAC)) = CSUM(CN(1:NAC)) + QSW(K, JWW)*C2(K, IWD(JWW), CN(1:NAC))
                                    end do
                                    if (QSUMM > 0.0) then
                                        TTR(JTT) = TSUM/QSUMM
                                        CTR(CN(1:NAC), JTT) = CSUM(CN(1:NAC))/QSUMM
                                    end if
                                end if ! SW 9/25/13 END MOVED CODE
                            else
!JSS(JBU)                 =  JSS(JBU)+1     ! SW 9/25/13
                                KTSW(JSS(JBU), JBU) = KTPU(JP)
                                KBSW(JSS(JBU), JBU) = KBPU(JP)
                                JB = JBU
                                POINT_SINK(JSS(JBU), JBU) = .true.
                                ID = IUPU(JP)
                                QSTR(JSS(JBU), JBU) = QPU(JP)
                                ESTR(JSS(JBU), JBU) = EPU(JP)
                                KT = KTWB(JWU)
                                JW = JWU
                                call DOWNSTREAM_WITHDRAWAL(JSS(JBU))
                                if (IDPU(JP) /= 0 .and. US(JBD) == IDPU(JP)) then
                                    QSUMM = 0.0;                                     TSUM = 0.0;                                     CSUM = 0.0
                                    do K = KT, KB(ID)
                                        QSUMM = QSUMM + QNEW(K)
                                        TSUM = TSUM + QNEW(K)*T2(K, ID)
                                        CSUM(CN(1:NAC)) = CSUM(CN(1:NAC)) + QNEW(K)*C2(K, ID, CN(1:NAC))
                                    end do
                                    if (QSUMM /= 0.0) then
                                        TINSUM(JBD) = (TSUM + TINSUM(JBD)*QINSUM(JBD))/(QSUMM + QINSUM(JBD))
                                        CINSUM(CN(1:NAC), JBD) = (CSUM(CN(1:NAC)) + CINSUM(CN(1:NAC), JBD)*QINSUM(JBD))/(QSUMM + QINSUM(JBD))
                                        QINSUM(JBD) = QINSUM(JBD) + QSUMM
                                    end if
                                end if
                                QSUM(JB) = 0.0;                                 TSUM = 0.0;                                 CSUM = 0.0
                                do K = KT, KB(ID)
                                    QSUM(JB) = QSUM(JB) + QOUT(K, JB)
                                    TSUM = TSUM + QOUT(K, JB)*T2(K, ID)
                                    CSUM(CN(1:NAC)) = CSUM(CN(1:NAC)) + QOUT(K, JB)*C2(K, ID, CN(1:NAC))
                                end do
                                if (QSUM(JB) /= 0.0) then
                                    TOUT(JB) = TSUM/QSUM(JB)
                                    COUT(CN(1:NAC), JB) = CSUM(CN(1:NAC))/QSUM(JB)
                                end if
                                if (IDPU(JP) /= 0) then ! SW 9/25/13 Moved code start
                                    if (US(JBD) /= IDPU(JP) .or. HEAD_FLOW(JBD) .or. UP_HEAD(JBD)) then
                                        JTT = JTT + 1
                                        QTR(JTT) = QPU(JP)
                                        ITR(JTT) = IDPU(JP)
                                        PLACE_QTR(JTT) = PPUC(JP) == " DENSITY"
                                        SPECIFY_QTR(JTT) = PPUC(JP) == " SPECIFY"
                                        if (SPECIFY_QTR(JTT)) then
                                            ELTRT(JTT) = ETPU(JP)
                                            ELTRB(JTT) = EBPU(JP)
                                        end if
                                        JBTR(JTT) = JBD
                                        TTR(JTT) = TOUT(JB)
                                        CTR(CN(1:NAC), JTT) = COUT(CN(1:NAC), JB)
                                    end if
                                end if ! Moved code end SW 9/25/13
                            end if
                        end if
                    end if
                end if
            end if
        end do
    end if
    if (PIPES) then
        YSS = YS
        VSS = VS
        VSTS = VST
        YSTS = YST
        DTPS = DTP
        QOLDS = QOLD
        call PIPE_FLOW() ! (NIT)
        do JP = 1, NPI

            if (dynpipe(jp) == "      ON") then ! SW 5/10/10
                qpi(jp) = qpi(jp)*bp(jp)
            end if


!****** Positive flows

            JLAT = 0
            JBU = JBUPI(JP)
            JBD = JBDPI(JP)
            tdgon = .false. ! cb 1/16/13
            if (QPI(JP) >= 0.0) then
                if (LATERAL_PIPE(JP)) then
                    JLAT = 1
                    JWW = JWW + 1
                    IWD(JWW) = IUPI(JP)
                    QWD(JWW) = QPI(JP)
                    KTWD(JWW) = KTUPI(JP)
                    KBWD(JWW) = KBUPI(JP)
                    EWD(JWW) = EUPI(JP)
                    JBWD(JWW) = JBU
                    I = MAX(CUS(JBWD(JWW)), IWD(JWW))
                    JB = JBWD(JWW)
                    JW = JWUPI(JP)
                    KT = KTWB(JW)
                    jwd = jww
                    call LATERAL_WITHDRAWAL() !(JWW)
                    do K = KTW(JWW), KBW(JWW)
                        QSS(K, I) = QSS(K, I) - QSW(K, JWW)
                    end do
                else
                    JSS(JBU) = JSS(JBU) + 1
                    KTSW(JSS(JBU), JBU) = KTDPI(JP)
                    KBSW(JSS(JBU), JBU) = KBDPI(JP)
                    JB = JBU
                    POINT_SINK(JSS(JBU), JBU) = .true.
                    ID = IUPI(JP)
                    QSTR(JSS(JBU), JBU) = QPI(JP)
                    ESTR(JSS(JBU), JBU) = EUPI(JP)
                    KT = KTWB(JWUPI(JP))
                    JW = JWUPI(JP)
                    call DOWNSTREAM_WITHDRAWAL(JSS(JBU))
                    if (IDPI(JP) /= 0 .and. US(JBD) == IDPI(JP)) then
                        QSUMM = 0.0;                         TSUM = 0.0;                         CSUM = 0.0
                        do K = KT, KB(ID)
                            QSUMM = QSUMM + QNEW(K)
                            TSUM = TSUM + QNEW(K)*T2(K, ID)
                            CSUM(CN(1:NAC)) = CSUM(CN(1:NAC)) + QNEW(K)*C2(K, ID, CN(1:NAC))
                        end do
                        if (QSUMM /= 0.0) then
                            TINSUM(JBD) = (TSUM + QINSUM(JBD)*TINSUM(JBD))/(QSUMM + QINSUM(JBD))
                            CINSUM(CN(1:NAC), JBD) = (CSUM(CN(1:NAC)) + QINSUM(JBD)*CINSUM(CN(1:NAC), JBD))/(QSUMM + QINSUM(JBD))
                            QINSUM(JBD) = QINSUM(JBD) + QSUMM
                        end if
                    end if
                    QSUM(JB) = 0.0;                     TSUM = 0.0;                     CSUM = 0.0
                    do K = KT, KB(ID)
                        QSUM(JB) = QSUM(JB) + QOUT(K, JB)
                        TSUM = TSUM + QOUT(K, JB)*T2(K, ID)
                        CSUM(CN(1:NAC)) = CSUM(CN(1:NAC)) + QOUT(K, JB)*C2(K, ID, CN(1:NAC))
                    end do
                    if (QSUM(JB) /= 0.0) then
                        TOUT(JB) = TSUM/QSUM(JB)
                        COUT(CN(1:NAC), JB) = CSUM(CN(1:NAC))/QSUM(JB)
                    end if
                end if
                if (IDPI(JP) /= 0) then
                    if (US(JBD) /= IDPI(JP) .or. HEAD_FLOW(JBD) .or. UP_HEAD(JBD)) then
                        JTT = JTT + 1
                        QTR(JTT) = QPI(JP)
                        ITR(JTT) = IDPI(JP)
                        PLACE_QTR(JTT) = PDPIC(JP) == " DENSITY"
                        SPECIFY_QTR(JTT) = PDPIC(JP) == " SPECIFY"
                        if (SPECIFY_QTR(JTT)) then
                            ELTRT(JTT) = ETDPI(JP)
                            ELTRB(JTT) = EBDPI(JP)
                        end if
                        JBTR(JTT) = JBD
                        if (JLAT == 1) then
                            TSUM = 0.0;                             QSUMM = 0.0;                             CSUM(CN(1:NAC)) = 0.0
                            do K = KTW(JWW), KBW(JWW)
                                QSUMM = QSUMM + QSW(K, JWW)
                                TSUM = TSUM + QSW(K, JWW)*T2(K, IWD(JWW))
                                CSUM(CN(1:NAC)) = CSUM(CN(1:NAC)) + QSW(K, JWW)*C2(K, IWD(JWW), CN(1:NAC))
                            end do
                            TTR(JTT) = TSUM/QSUMM
                            CTR(CN(1:NAC), JTT) = CSUM(CN(1:NAC))/QSUMM
                        else
                            TTR(JTT) = TOUT(JB)
                            CTR(CN(1:NAC), JTT) = COUT(CN(1:NAC), JB)
                        end if
                    else
                        if (LATERAL_PIPE(JP)) then
                            TSUM = 0.0;                             QSUMM = 0.0;                             CSUM = 0.0
                            ILAT(JWW) = 1
                            JB = JBD
                            do K = KTW(JWW), KBW(JWW)
                                QSUMM = QSUMM + QSW(K, JWW)
                                TSUM = TSUM + QSW(K, JWW)*T2(K, IWD(JWW))
                                CSUM(CN(1:NAC)) = CSUM(CN(1:NAC)) + QSW(K, JWW)*C2(K, IWD(JWW), CN(1:NAC))
                            end do
                            TINSUM(JB) = (TINSUM(JB)*QINSUM(JB) + TSUM)/(QSUMM + QINSUM(JB))
                            CINSUM(CN(1:NAC), JB) = (CINSUM(CN(1:NAC), JB)*QINSUM(JB) + CSUM(CN(1:NAC)))/(QSUMM + QINSUM(JB))
                            QINSUM(JB) = QSUMM + QINSUM(JB)
                        end if
                    end if
                end if
            else
                JTT = JTT + 1
                JWW = JWW + 1
                IWD(JWW) = IDPI(JP)
                ITR(JTT) = IUPI(JP)
                QTR(JTT) = -QPI(JP)
                QWD(JWW) = -QPI(JP)
                KTWD(JWW) = KTDPI(JP)
                KBWD(JWW) = KBDPI(JP)
                EWD(JWW) = EDPI(JP)
                PLACE_QTR(JTT) = PUPIC(JP) == " DENSITY"
                SPECIFY_QTR(JTT) = PUPIC(JP) == " SPECIFY"
                if (SPECIFY_QTR(JTT)) then
                    ELTRT(JTT) = ETUPI(JP)
                    ELTRB(JTT) = EBUPI(JP)
                end if
                JBTR(JTT) = JBU
                JBWD(JWW) = JBD
                I = MAX(CUS(JBWD(JWW)), IWD(JWW))
                JB = JBWD(JWW)
                JW = JWDPI(JP)
                KT = KTWB(JW)
                jwd = jww
                call LATERAL_WITHDRAWAL() !(JWW)
                do K = KTW(JWW), KBW(JWW)
                    QSS(K, I) = QSS(K, I) - QSW(K, JWW)
                end do
                if (IDPI(JP) /= 0) then
                    TSUM = 0.0;                     QSUMM = 0.0;                     CSUM = 0.0
                    do K = KTW(JWW), KBW(JWW)
                        QSUMM = QSUMM + QSW(K, JWW)
                        TSUM = TSUM + QSW(K, JWW)*T2(K, IWD(JWW))
                        CSUM(CN(1:NAC)) = CSUM(CN(1:NAC)) + QSW(K, JWW)*C2(K, IWD(JWW), CN(1:NAC))
                    end do
                    TTR(JTT) = TSUM/QSUMM
                    CTR(CN(1:NAC), JTT) = CSUM(CN(1:NAC))/QSUMM
                end if
            end if
        end do
    end if
    if (GATES) then
        call GATE_FLOW()
        if (TDGTA) then
            call TDGtarget()
        end if ! tdgtarget 
        if (SYSTDG) then
            call SYSTDG_TDG()
        end if ! SYSTDG - CALCULATE TDG_TDG (%)
        do JG = 1, NGT

!****** Positive flows

            JLAT = 0
            JBU = JBUGT(JG)
            JBD = JBDGT(JG)
            tdgon = .false. ! cb 1/16/13
            jsg = jg
            nnsg = 1
            if (cac(ndo) == "      ON" .and. gasgtc(jg) == "      ON") then
                tdgon = .true.
            end if
            if (QGT(JG) >= 0.0) then
                if (LATERAL_GATE(JG)) then
!           JLAT      = 1
                    JWW = JWW + 1
                    IWD(JWW) = IUGT(JG)
                    QWD(JWW) = QGT(JG)
                    KTWD(JWW) = KTUGT(JG)
                    KBWD(JWW) = KBUGT(JG)
                    EWD(JWW) = EGT(JG)
                    if (DYNGTC(JG) == "     ZGT" .and. GT2CHAR == "EGT2ELEV") then ! SW 2/25/11
                        if (EGT2(JG) /= 0.0) then
                            EWD(JWW) = EGT2(JG)
                        end if
                    end if
                    JBWD(JWW) = JBU
                    I = MAX(CUS(JBWD(JWW)), IWD(JWW))
                    JW = JWUGT(JG)
                    JB = JBWD(JWW)
                    KT = KTWB(JW)
                    jwd = jww
                    call LATERAL_WITHDRAWAL() !(JWW)
                    do K = KTW(JWW), KBW(JWW)
                        QSS(K, I) = QSS(K, I) - QSW(K, JWW)
                    end do
                    if (IDGT(JG) /= 0) then
                        CSUM(CN(1:NAC)) = 0.0;                         TSUM = 0.0;                         QSUMM = 0.0
                        JTT = JTT + 1 !  SW 4/1/09
                        QTR(JTT) = QGT(JG) !  SW 4/1/09
                        ITR(JTT) = IDGT(JG) !  SW 4/1/09
                        PLACE_QTR(JTT) = PDGTC(JG) == " DENSITY" !  SW 4/1/09
                        SPECIFY_QTR(JTT) = PDGTC(JG) == " SPECIFY" !  SW 4/1/09
                        if (SPECIFY_QTR(JTT)) then !  SW 4/1/09
                            ELTRT(JTT) = ETDGT(JG) !  SW 4/1/09
                            ELTRB(JTT) = EBDGT(JG) !  SW 4/1/09
                        end if
                        JBTR(JTT) = JBD !  SW 4/1/09
                        do K = KTW(JWW), KBW(JWW)
                            QSUMM = QSUMM + QSW(K, JWW)
                            TSUM = TSUM + QSW(K, JWW)*T2(K, IWD(JWW))
                            CSUM(CN(1:NAC)) = CSUM(CN(1:NAC)) + QSW(K, JWW)*C2(K, IWD(JWW), CN(1:NAC))
                        end do
                        if (QSUMM == 0.0) then
                            TTR(JTT) = 0.0
                            CTR(:, JTT) = 0.0
                        else
                            TTR(JTT) = TSUM/QSUMM
                            do JC = 1, NAC
                                CTR(CN(JC), JTT) = CSUM(CN(JC))/QSUMM
                                if (CN(JC) == NDO .and. GASGTC(JG) == "      ON" .and. QGT(JG) > 0.0) then
                                    TDG_GATE(JWW, JG) = .true.
!
! systdg 
                                    if (SYSTDG) then
                                        if (GTNAME(JG)) then
                                            call UPDATE_TDGC(0, PALT(ID), JG, TTR(JTT), CTR(CN(JC), JTT)) ! O2        
                                        else
                                            call TOTAL_DISSOLVED_GAS(0, PALT(ID), 1, JG, TTR(JTT), CTR(CN(JC), JTT)) ! O2
                                        end if
                                    else
                                        call TOTAL_DISSOLVED_GAS(0, PALT(ID), 1, JG, TTR(JTT), CTR(CN(JC), JTT)) ! O2
                                    end if
!
                                end if
                                if (CN(JC) == NN2 .and. GASGTC(JG) == "      ON" .and. QGT(JG) > 0.0) then
                                    TDG_GATE(JWW, JG) = .true.
!
! systdg 
                                    if (SYSTDG) then
                                        if (GTNAME(JG)) then
                                            call UPDATE_TDGC(1, PALT(ID), JG, TTR(JTT), CTR(CN(JC), JTT)) ! N2         
                                        else
                                            call TOTAL_DISSOLVED_GAS(1, PALT(ID), 1, JG, TTR(JTT), CTR(CN(JC), JTT)) ! N2
                                        end if
                                    else
                                        call TOTAL_DISSOLVED_GAS(1, PALT(ID), 1, JG, TTR(JTT), CTR(CN(JC), JTT)) ! N2
                                    end if
!
                                end if
                                if (CN(JC) == NDGP .and. GASGTC(JG) == "      ON" .and. QGT(JG) > 0.0) then !8/2020 TDGP
                                    TDG_GATE(JWW, JG) = .true.
!
                                    if (SYSTDG) then
                                        if (GTNAME(JG)) then
                                            call UPDATE_TDGC(2, PALT(ID), JG, TTR(JTT), CTR(CN(JC), JTT))
                                        else
                                            call TOTAL_DISSOLVED_GAS(2, PALT(ID), 1, JG, TTR(JTT), CTR(CN(JC), JTT))
                                        end if
                                    else
                                        call TOTAL_DISSOLVED_GAS(2, PALT(ID), 1, JG, TTR(JTT), CTR(CN(JC), JTT))
                                    end if
!
                                end if
                            end do
                        end if
                    else
                        if (CAC(NDO) == "      ON" .and. GASGTC(JG) == "      ON" .and. QGT(JG) > 0.0) then
                            TDG_GATE(JWW, JG) = .true.
                        end if
                    end if
                else
                    JSS(JBU) = JSS(JBU) + 1
                    KTSW(JSS(JBU), JBU) = KTUGT(JG)
                    KBSW(JSS(JBU), JBU) = KBUGT(JG)
                    JB = JBU
                    POINT_SINK(JSS(JBU), JBU) = .true.
                    ID = IUGT(JG)
                    ESTR(JSS(JBU), JBU) = EGT(JG)
                    if (DYNGTC(JG) == "     ZGT" .and. GT2CHAR == "EGT2ELEV") then ! SW 2/25/11
                        if (EGT2(JG) /= 0.0) then
                            ESTR(JSS(JBU), JBU) = EGT2(JG)
                        end if
                    end if
                    QSTR(JSS(JBU), JBU) = QGT(JG)
                    KT = KTWB(JWUGT(JG))
                    JW = JWUGT(JG)
                    call DOWNSTREAM_WITHDRAWAL(JSS(JBU))
                    QSUM(JB) = 0.0;                     TSUM = 0.0;                     CSUM = 0.0
                    do K = KT, KB(ID)
                        QSUM(JB) = QSUM(JB) + QOUT(K, JB)
                        TSUM = TSUM + QOUT(K, JB)*T2(K, ID)
                        do JC = 1, NAC
                            if (CN(JC) == NDO .and. CAC(NDO) == "      ON" .and. GASGTC(JG) == "      ON" .and. QGT(JG) > 0.0) then ! MM 5/21/2009
                                T2R4 = T2(K, ID)
                                CGAS = C2(K, ID, CN(JC)) ! MM 5/21/2009                  
!
! systdg 
                                if (SYSTDG) then
                                    if (GTNAME(JG)) then
                                        call UPDATE_TDGC(0, PALT(ID), JG, T2R4, CGAS)
                                    else
                                        call TOTAL_DISSOLVED_GAS(0, PALT(ID), 1, JG, T2R4, CGAS)
                                    end if
                                else ! MM 5/21/2009                  
                                    call TOTAL_DISSOLVED_GAS(0, PALT(ID), 1, JG, T2R4, CGAS)
                                end if
!
                                CSUM(CN(JC)) = CSUM(CN(JC)) + QOUT(K, JB)*CGAS
!ELSEIF (CN(JC) == NGN2 .AND. CAC(NGN2) == '      ON' .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN   ! SW 10/27/15
                            else
                                if (CN(JC) == NN2 .and. GASGTC(JG) == "      ON" .and. QGT(JG) > 0.0) then ! SW 10/27/15
                                    if (CAC(NN2) == "      ON") then ! cb 1/13/16
                                        T2R4 = T2(K, ID)
                                        CGAS = C2(K, ID, CN(JC))
!
! systdg 
                                        if (SYSTDG) then
                                            if (GTNAME(JG)) then
                                                call UPDATE_TDGC(1, PALT(ID), JG, T2R4, CGAS)
                                            else
                                                call TOTAL_DISSOLVED_GAS(1, PALT(ID), 1, JG, T2R4, CGAS)
                                            end if
                                        else
                                            call TOTAL_DISSOLVED_GAS(1, PALT(ID), 1, JG, T2R4, CGAS)
                                        end if
!
                                        CSUM(CN(JC)) = CSUM(CN(JC)) + QOUT(K, JB)*CGAS
                                    end if
                                else
                                    if (CN(JC) == NDGP .and. GASGTC(JG) == "      ON" .and. QGT(JG) > 0.0) then !8/2020 TDGP
                                        if (CAC(NDGP) == "      ON") then
                                            T2R4 = T2(K, ID)
                                            CGAS = C2(K, ID, CN(JC))
!               
                                            if (SYSTDG) then
                                                if (GTNAME(JG)) then
                                                    call UPDATE_TDGC(2, PALT(ID), JG, T2R4, CGAS)
                                                else
                                                    call TOTAL_DISSOLVED_GAS(2, PALT(ID), 1, JG, T2R4, CGAS)
                                                end if
                                            else
                                                call TOTAL_DISSOLVED_GAS(2, PALT(ID), 1, JG, T2R4, CGAS)
                                            end if
!
                                            CSUM(CN(JC)) = CSUM(CN(JC)) + QOUT(K, JB)*CGAS
                                        end if
                                    else
                                        CSUM(CN(JC)) = CSUM(CN(JC)) + QOUT(K, JB)*C2(K, ID, CN(JC))
                                    end if
                                end if
                            end if
                        end do
                    end do
                    if (QSUM(JB) /= 0.0) then
                        TOUT(JB) = TSUM/QSUM(JB)
                        COUT(CN(1:NAC), JB) = CSUM(CN(1:NAC))/QSUM(JB)
                    end if
                    if (IDGT(JG) /= 0 .and. US(JBD) == IDGT(JG)) then
                        QSUMM = 0.0
                        TSUM = 0.0
                        CSUM = 0.0
                        do K = KT, KB(ID)
                            QSUMM = QSUMM + QNEW(K)
                            TSUM = TSUM + QNEW(K)*T2(K, ID)
                            do JC = 1, NAC
                                if (CN(JC) == NDO .and. CAC(NDO) == "      ON" .and. GASGTC(JG) == "      ON" .and. QGT(JG) > 0.0) then ! MM 5/21/2009
                                    T2R4 = T2(K, ID)
                                    CGAS = C2(K, ID, CN(JC)) ! MM 5/21/2009
!
! systdg 
                                    if (SYSTDG) then
                                        if (GTNAME(JG)) then
                                            call UPDATE_TDGC(0, PALT(ID), JG, T2R4, CGAS) ! O2       
                                        else
                                            call TOTAL_DISSOLVED_GAS(0, PALT(ID), 1, JG, T2R4, CGAS) ! O2
                                        end if
                                    else ! MM 5/21/2009
                                        call TOTAL_DISSOLVED_GAS(0, PALT(ID), 1, JG, T2R4, CGAS) ! O2
                                    end if
!
                                    CSUM(CN(JC)) = CSUM(CN(JC)) + QNEW(K)*CGAS
!ELSEIF (CN(JC) == NGN2 .AND. CAC(NGN2) == '      ON' .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN   ! SW 10/27/15
                                else
                                    if (CN(JC) == NN2 .and. GASGTC(JG) == "      ON" .and. QGT(JG) > 0.0) then ! SW 10/27/15
                                        if (CAC(NN2) == "      ON") then
                                            T2R4 = T2(K, ID)
                                            CGAS = C2(K, ID, CN(JC))
!
! systdg 
                                            if (SYSTDG) then
                                                if (GTNAME(JG)) then
                                                    call UPDATE_TDGC(1, PALT(ID), JG, T2R4, CGAS) ! N2           
                                                else
                                                    call TOTAL_DISSOLVED_GAS(1, PALT(ID), 1, JG, T2R4, CGAS) ! N2
                                                end if
                                            else
                                                call TOTAL_DISSOLVED_GAS(1, PALT(ID), 1, JG, T2R4, CGAS) ! N2
                                            end if
!
                                            CSUM(CN(JC)) = CSUM(CN(JC)) + QNEW(K)*CGAS
                                        end if
                                    else
                                        if (CN(JC) == NDGP .and. GASGTC(JG) == "      ON" .and. QGT(JG) > 0.0) then !8/2020 TDGP
                                            if (CAC(NDGP) == "      ON") then
                                                T2R4 = T2(K, ID)
                                                CGAS = C2(K, ID, CN(JC))
!
                                                if (SYSTDG) then
                                                    if (GTNAME(JG)) then
                                                        call UPDATE_TDGC(2, PALT(ID), JG, T2R4, CGAS)
                                                    else
                                                        call TOTAL_DISSOLVED_GAS(2, PALT(ID), 1, JG, T2R4, CGAS)
                                                    end if
                                                else
                                                    call TOTAL_DISSOLVED_GAS(2, PALT(ID), 1, JG, T2R4, CGAS)
                                                end if
!
                                                CSUM(CN(JC)) = CSUM(CN(JC)) + QNEW(K)*CGAS
                                            end if
                                        else
                                            CSUM(CN(JC)) = CSUM(CN(JC)) + QNEW(K)*C2(K, ID, CN(JC))
                                        end if
                                    end if
                                end if
                            end do
                        end do
                        if (QSUMM /= 0.0) then
                            TINSUM(JBD) = (TSUM + QINSUM(JBD)*TINSUM(JBD))/(QSUMM + QINSUM(JBD))
                            CINSUM(CN(1:NAC), JBD) = (CSUM(CN(1:NAC)) + QINSUM(JBD)*CINSUM(CN(1:NAC), JBD))/(QSUMM + QINSUM(JBD))
                            QINSUM(JBD) = QINSUM(JBD) + QSUMM
                        end if
                    else
                        if (IDGT(JG) /= 0) then
                            JTT = JTT + 1
                            QTR(JTT) = QGT(JG)
                            ITR(JTT) = IDGT(JG)
                            PLACE_QTR(JTT) = PDGTC(JG) == " DENSITY"
                            SPECIFY_QTR(JTT) = PDGTC(JG) == " SPECIFY"
                            if (SPECIFY_QTR(JTT)) then
                                ELTRT(JTT) = ETDGT(JG)
                                ELTRB(JTT) = EBDGT(JG)
                            end if
                            JBTR(JTT) = JBD
                            TTR(JTT) = TOUT(JB)
                            do JC = 1, NAC
                                CTR(CN(JC), JTT) = COUT(CN(JC), JB)
                                if (CN(JC) == NDO .and. GASGTC(JG) == "      ON" .and. QGT(JG) > 0.0) then
!
! systdg 
                                    if (SYSTDG) then
                                        if (GTNAME(JG)) then
                                            call UPDATE_TDGC(0, PALT(ID), JG, TTR(JTT), CTR(CN(JC), JTT)) ! O2          
                                        else
                                            call TOTAL_DISSOLVED_GAS(0, PALT(ID), 0, JS, TTR(JTT), CTR(CN(JC), JTT)) ! O2
                                        end if
                                    else
                                        call TOTAL_DISSOLVED_GAS(0, PALT(ID), 0, JS, TTR(JTT), CTR(CN(JC), JTT)) ! O2
                                    end if
!
                                end if
                                if (CN(JC) == NN2 .and. GASGTC(JG) == "      ON" .and. QGT(JG) > 0.0) then
! 
! systdg 
                                    if (SYSTDG) then
                                        if (GTNAME(JG)) then
                                            call UPDATE_TDGC(1, PALT(ID), JG, TTR(JTT), CTR(CN(JC), JTT)) ! N2       
                                        else
                                            call TOTAL_DISSOLVED_GAS(1, PALT(ID), 0, JS, TTR(JTT), CTR(CN(JC), JTT)) ! N2
                                        end if
                                    else
                                        call TOTAL_DISSOLVED_GAS(1, PALT(ID), 0, JS, TTR(JTT), CTR(CN(JC), JTT)) ! N2
                                    end if
!
                                end if
                                if (CN(JC) == NDGP .and. GASGTC(JG) == "      ON" .and. QGT(JG) > 0.0) then !8/2020 TDGP
! 
                                    if (SYSTDG) then
                                        if (GTNAME(JG)) then
                                            call UPDATE_TDGC(2, PALT(ID), JG, TTR(JTT), CTR(CN(JC), JTT))
                                        else
                                            call TOTAL_DISSOLVED_GAS(2, PALT(ID), 0, JS, TTR(JTT), CTR(CN(JC), JTT))
                                        end if
                                    else
                                        call TOTAL_DISSOLVED_GAS(2, PALT(ID), 0, JS, TTR(JTT), CTR(CN(JC), JTT))
                                    end if
!
                                end if
                            end do
                        end if
                    end if
                end if
            else
                if (QGT(JG) < 0.0) then
                    JTT = JTT + 1
                    JWW = JWW + 1
                    IWD(JWW) = IDGT(JG)
                    ITR(JTT) = IUGT(JG)
                    QTR(JTT) = -QGT(JG)
                    QWD(JWW) = -QGT(JG)
                    KTWD(JWW) = KTDGT(JG)
                    KBWD(JWW) = KBDGT(JG)
                    EWD(JWW) = EGT(JG)
                    PLACE_QTR(JTT) = PUGTC(JG) == " DENSITY"
                    SPECIFY_QTR(JTT) = PUGTC(JG) == " SPECIFY"
                    if (SPECIFY_QTR(JTT)) then
                        ELTRT(JTT) = ETUGT(JG)
                        ELTRB(JTT) = EBUGT(JG)
                    end if
                    JBTR(JTT) = JBU
                    JBWD(JWW) = JBD
                    I = MAX(CUS(JBWD(JWW)), IWD(JWW))
                    JW = JWDGT(JG)
                    JB = JBWD(JWW)
                    KT = KTWB(JW)
                    jwd = jww
                    call LATERAL_WITHDRAWAL() !(JWW)
                    do K = KTW(JWW), KBW(JWW)
                        QSS(K, I) = QSS(K, I) - QSW(K, JWW)
                    end do
                    if (IDGT(JG) /= 0) then
                        CSUM(CN(1:NAC)) = 0.0;                         TSUM = 0.0;                         QSUMM = 0.0
                        do K = KTW(JWW), KBW(JWW)
                            QSUMM = QSUMM + QSW(K, JWW)
                            TSUM = TSUM + QSW(K, JWW)*T2(K, IWD(JWW))
                            CSUM(CN(1:NAC)) = CSUM(CN(1:NAC)) + QSW(K, JWW)*C2(K, IWD(JWW), CN(1:NAC))
                        end do
                        TTR(JTT) = TSUM/QSUMM
                        do JC = 1, NAC
                            CTR(CN(JC), JTT) = CSUM(CN(JC))/QSUMM
                            if (CN(JC) == NDO .and. GASGTC(JG) == "      ON" .and. QGT(JG) > 0.0) then
                                TDG_GATE(JWW, JG) = .true.
!
! systdg 
                                if (SYSTDG) then
                                    if (GTNAME(JG)) then
                                        call UPDATE_TDGC(0, PALT(ID), JG, TTR(JTT), CTR(CN(JC), JTT))
                                    else
                                        call TOTAL_DISSOLVED_GAS(0, PALT(ID), 1, JG, TTR(JTT), CTR(CN(JC), JTT))
                                    end if
                                else
                                    call TOTAL_DISSOLVED_GAS(0, PALT(ID), 1, JG, TTR(JTT), CTR(CN(JC), JTT))
                                end if
!
                            end if
                            if (CN(JC) == NN2 .and. GASGTC(JG) == "      ON" .and. QGT(JG) > 0.0) then
                                TDG_GATE(JWW, JG) = .true.
!
! systdg 
                                if (SYSTDG) then
                                    if (GTNAME(JG)) then
                                        call UPDATE_TDGC(1, PALT(ID), JG, TTR(JTT), CTR(CN(JC), JTT))
                                    else
                                        call TOTAL_DISSOLVED_GAS(1, PALT(ID), 1, JG, TTR(JTT), CTR(CN(JC), JTT))
                                    end if
                                else
                                    call TOTAL_DISSOLVED_GAS(1, PALT(ID), 1, JG, TTR(JTT), CTR(CN(JC), JTT))
                                end if
!
                            end if
                            if (CN(JC) == NDGP .and. GASGTC(JG) == "      ON" .and. QGT(JG) > 0.0) then !8/2020 TDGP
                                TDG_GATE(JWW, JG) = .true.
!
                                if (SYSTDG) then
                                    if (GTNAME(JG)) then
                                        call UPDATE_TDGC(2, PALT(ID), JG, TTR(JTT), CTR(CN(JC), JTT))
                                    else
                                        call TOTAL_DISSOLVED_GAS(2, PALT(ID), 1, JG, TTR(JTT), CTR(CN(JC), JTT))
                                    end if
                                else
                                    call TOTAL_DISSOLVED_GAS(2, PALT(ID), 1, JG, TTR(JTT), CTR(CN(JC), JTT))
                                end if
!
                            end if
                        end do
                    else
                        if (CAC(NDO) == "      ON" .and. GASGTC(JG) == "      ON" .and. QGT(JG) > 0.0) then
                            TDG_GATE(JWW, JG) = .true.
                        end if
                    end if
                end if
            end if
        end do
! systdg - Add power house release tdg update
        if (SYSTDG) then
            if (POWNO > 0 .and. TDGLOC == "     REL") then
                do ip = 1, POWNO
                    if (QGT(POWGTNO(ip)) > 0.0 .and. TDG_ROSP > 0.0) then
                        if (LATERAL_GATE(POWGTNO(ip))) then
                            if (CONSTITUENTS) then
                                call UPDATE_TDGC(0, palt(IUGT(POWGTNO(ip))), POWGTNO(ip), tavgw(POWGTNO(ip)), CAVGW(POWGTNO(ip), NDO))
                                call UPDATE_TDGC(1, palt(IUGT(POWGTNO(ip))), POWGTNO(ip), tavgw(POWGTNO(ip)), CAVGW(POWGTNO(ip), NN2))
                                call UPDATE_TDGC(2, palt(IUGT(POWGTNO(ip))), POWGTNO(ip), tavgw(POWGTNO(ip)), CAVGW(POWGTNO(ip), NDGP))
                            end if
                            if (DERIVED_CALC) then
                                CDAVGW(POWGTNO(ip), O2DG_DER) = tdg_tdg
                                CDAVGW(POWGTNO(ip), TDG_DER) = tdg_tdg
                            end if
                        else
                            if (CONSTITUENTS) then
                                call UPDATE_TDGC(0, PALT(IUGT(POWGTNO(ip))), POWGTNO(ip), Tavg(POWGTNO(ip), JBUGT(POWGTNO(ip))), CAVG(POWGTNO(ip), JBUGT(POWGTNO(ip)), NDO))
                                call UPDATE_TDGC(1, PALT(IUGT(POWGTNO(ip))), POWGTNO(ip), Tavg(POWGTNO(ip), JBUGT(POWGTNO(ip))), CAVG(POWGTNO(ip), JBUGT(POWGTNO(ip)), NN2))
                                call UPDATE_TDGC(2, PALT(IUGT(POWGTNO(ip))), POWGTNO(ip), Tavg(POWGTNO(ip), JBUGT(POWGTNO(ip))), CAVG(POWGTNO(ip), JBUGT(POWGTNO(ip)), NDGP))
                            end if
                            if (DERIVED_CALC) then
                                CDAVG(POWGTNO(ip), JBUGT(POWGTNO(ip)), O2DG_DER) = tdg_tdg
                                CDAVG(POWGTNO(ip), JBUGT(POWGTNO(ip)), TDG_DER) = tdg_tdg
                            end if
                        end if
                    else
                        if (QGT(POWGTNO(ip)) < 0.0 .and. TDG_ROSP > 0.0) then ! QGT<0.0
                            if (CONSTITUENTS) then
                                call UPDATE_TDGC(0, palt(IUGT(POWGTNO(ip))), POWGTNO(ip), tavgw(POWGTNO(ip)), CAVGW(POWGTNO(ip), NDO))
                                call UPDATE_TDGC(1, palt(IUGT(POWGTNO(ip))), POWGTNO(ip), tavgw(POWGTNO(ip)), CAVGW(POWGTNO(ip), NN2))
                                call UPDATE_TDGC(2, palt(IUGT(POWGTNO(ip))), POWGTNO(ip), tavgw(POWGTNO(ip)), CAVGW(POWGTNO(ip), NDGP))
                            end if
                            if (DERIVED_CALC) then
                                CDAVGW(POWGTNO(ip), O2DG_DER) = tdg_tdg
                                CDAVGW(POWGTNO(ip), TDG_DER) = tdg_tdg
                            end if
                        end if
                    end if
                end do
            end if
            if (POWNO > 0 .and. TDGLOC == "     REL") then
                do ifl = 1, FLNO
                    if (QGT(FLGTNO(ifl)) > 0.0 .and. TDG_ROSP > 0.0) then
                        if (CONSTITUENTS) then
                            call UPDATE_TDGC(0, PALT(IUGT(FLGTNO(ifl))), FLGTNO(ifl), Tavg(FLGTNO(ifl), JBUGT(FLGTNO(ifl))), CAVG(FLGTNO(ifl), JBUGT(FLGTNO(ifl)), NDO))
                            call UPDATE_TDGC(1, PALT(IUGT(FLGTNO(ifl))), FLGTNO(ifl), Tavg(FLGTNO(ifl), JBUGT(FLGTNO(ifl))), CAVG(FLGTNO(ifl), JBUGT(FLGTNO(ifl)), NN2))
                            call UPDATE_TDGC(2, PALT(IUGT(FLGTNO(ifl))), FLGTNO(ifl), Tavg(FLGTNO(ifl), JBUGT(FLGTNO(ifl))), CAVG(FLGTNO(ifl), JBUGT(FLGTNO(ifl)), NDGP))
                        end if
                        if (DERIVED_CALC) then
                            CDAVG(FLGTNO(ifl), JBUGT(FLGTNO(ifl)), O2DG_DER) = tdg_tdg
                            CDAVG(FLGTNO(ifl), JBUGT(FLGTNO(ifl)), TDG_DER) = tdg_tdg
                        end if
                    else
                        if (QGT(FLGTNO(ifl)) < 0.0 .and. TDG_ROSP > 0.0) then
                            if (CONSTITUENTS) then
                                call UPDATE_TDGC(0, palt(IUGT(FLGTNO(ifl))), FLGTNO(ifl), tavgw(FLGTNO(ifl)), CAVGW(FLGTNO(ifl), NDO))
                                call UPDATE_TDGC(1, palt(IUGT(FLGTNO(ifl))), FLGTNO(ifl), tavgw(FLGTNO(ifl)), CAVGW(FLGTNO(ifl), NN2))
                                call UPDATE_TDGC(2, palt(IUGT(FLGTNO(ifl))), FLGTNO(ifl), tavgw(FLGTNO(ifl)), CAVGW(FLGTNO(ifl), NDGP))
                            end if
                            if (DERIVED_CALC) then
                                CDAVGW(FLGTNO(ifl), O2DG_DER) = tdg_tdg
                                CDAVGW(FLGTNO(ifl), TDG_DER) = tdg_tdg
                            end if
                        end if
                    end if
                end do
            end if
        end if
! systdg - Add power house release tdg update
    end if

    tdgon = .false. ! cb 1/17/13
    tributaries = jtt > 0
    withdrawals = jww > 0


    do JW = 1, NWB
        do JB = BS(JW), BE(JW)
            if (BR_INACTIVE(JB)) then
! CONVERT INFLOWS TO TRIBS SET TO THE CUS(1) LOCATION
                JTT = JTT + 1
!ITR(JTT)         =  CUS(1)   ! HARDWIRED TO FIRST BRANCH
                ITR(JTT) = CUS(jbdn(jw)) ! changed HARDWIRE TO JBDN BRANCH  ! cb 11/20/19
                QTR(JTT) = QIN(JB)
                TTR(JTT) = TIN(JB)
                do JC = 1, NAC
                    CTR(CN(JC), JTT) = CIN(CN(JC), JB)
                end do
!  PLACE_QTR(JTT)   =  '   DISTR'
                PLACE_QTR(JTT) = .false. !SR 01/22/2018
                JBTR(JTT) = 1
            end if
        end do
    end do



    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            if (BR_INACTIVE(JB)) then
                cycle
            end if ! SW 6/12/2017
            IU = CUS(JB)
            ID = DS(JB)
            if (EVAPORATION(JW)) then
                EVBR(JB) = 0.0
                do I = IU, ID
                    FW = AFW(JW) + BFW(JW)*WIND2(I)**CFW(JW)
                    if (RH_EVAP(JW)) then
                        EA = EXP(2.3026*(7.5*TDEW(JW)/(TDEW(JW) + 237.3) + 0.6609))
                        ES = EXP(2.3026*(7.5*T2(KT, I)/(T2(KT, I) + 237.3) + 0.6609))
                        if (TDEW(JW) < 0.0) then
                            EA = EXP(2.3026*(9.5*TDEW(JW)/(TDEW(JW) + 265.5) + 0.6609))
                        end if
                        if (T2(KT, I) < 0.0) then
                            ES = EXP(2.3026*(9.5*T2(KT, I)/(T2(KT, I) + 265.5) + 0.6609))
                        end if
                        TAIRV = (TAIR(JW) + 273.0)/(1.0 - 0.378*EA/760.0)
                        DTV = (T2(KT, I) + 273.0)/(1.0 - 0.378*ES/760.0) - TAIRV
                        DTVL = 0.0084*WIND2(I)**3
                        if (DTV < DTVL) then
                            DTV = DTVL
                        end if
                        FW = 3.59*DTV**0.3333333 + 4.26*WIND2(I)
                    end if
                    TM = (T2(KT, I) + TDEW(JW))*0.5
                    VPTG = 0.35 + 0.015*TM + 0.0012*TM*TM
                    EV(I) = VPTG*(T2(KT, I) - TDEW(JW))*FW*BI(KT, I)*DLX(I)/2.45E9
                    if (EV(I) < 0.0 .or. ICE(I)) then
                        EV(I) = 0.0
                    end if
                    QSS(KT, I) = QSS(KT, I) - EV(I)
                    EVBR(JB) = EVBR(JB) + EV(I)
                end do
            end if
            if (PRECIPITATION(JW)) then
                QPRBR(JB) = 0.0
                do I = IU, ID
                    QPR(I) = PR(JB)*BI(KT, I)*DLX(I)
                    QPRBR(JB) = QPRBR(JB) + QPR(I)
                    QSS(KT, I) = QSS(KT, I) + QPR(I)
                end do
            end if
            if (TRIBUTARIES) then
                do JT = 1, JTT

!********** Inflow fractions

                    if (JB == JBTR(JT)) then
                        I = MAX(ITR(JT), IU)
                        QTRF(KT:KB(I), JT) = 0.0
                        if (PLACE_QTR(JT)) then

!************** Inflow layer

                            SSTOT = 0.0
                            do J = NSSS, NSSE
                                SSTOT = SSTOT + CTR(J, JT)
                            end do
                            RHOTR = DENSITY(TTR(JT), CTR(NTDS, JT), SSTOT)
                            K = KT
                            do while (RHOTR > RHO(K, I) .and. K < KB(I))
                                K = K + 1
                            end do
                            KTTR(JT) = K
                            KBTR(JT) = K

!************** Layer inflows

                            VQTR = QTR(JT)*DLT
                            VQTRI = VQTR
                            QTRFR = 1.0
                            INCR = -1
                            do while (QTRFR > 0.0)
                                if (K <= KB(I)) then
                                    V1 = VOL(K, I)
                                    if (VQTR > 0.5*V1) then
                                        QTRF(K, JT) = 0.5*V1/VQTRI
                                        QTRFR = QTRFR - QTRF(K, JT)
                                        VQTR = VQTR - QTRF(K, JT)*VQTRI
                                        if (K == KT) then
                                            K = KBTR(JT)
                                            INCR = 1
                                        end if
                                    else
                                        QTRF(K, JT) = QTRFR
                                        QTRFR = 0.0
                                    end if
                                    if (INCR < 0) then
                                        KTTR(JT) = K
                                    end if
                                    if (INCR > 0) then
                                        KBTR(JT) = MIN(KB(I), K)
                                    end if
                                    K = K + INCR
                                else
                                    QTRF(KT, JT) = QTRF(KT, JT) + QTRFR
                                    QTRFR = 0.0
                                end if
                            end do
                        else
                            if (SPECIFY_QTR(JT)) then
                                KTTR(JT) = 2
!             DO WHILE (EL(KTTR(JT),I) > ELTRT(JT))
                                do while (EL(KTTR(JT), I) > ELTRT(JT) .and. EL(KTTR(JT) + 1, I) > ELTRT(JT)) ! SW 10/3/13
                                    KTTR(JT) = KTTR(JT) + 1
                                end do
                                KBTR(JT) = KMX - 1
                                do while (EL(KBTR(JT), I) < ELTRB(JT))
                                    KBTR(JT) = KBTR(JT) - 1
                                end do
                            else
                                KTTR(JT) = KT
                                KBTR(JT) = KB(I)
                            end if
                            KTTR(JT) = MAX(KT, KTTR(JT))
                            KBTR(JT) = MIN(KB(I), KBTR(JT))
                            if (KBTR(JT) < KTTR(JT)) then
                                KBTR(JT) = KTTR(JT)
                            end if
                            BHSUM = 0.0
                            do K = KTTR(JT), KBTR(JT)
                                BHSUM = BHSUM + BH2(K, I)
                            end do
                            do K = KTTR(JT), KBTR(JT)
                                QTRF(K, JT) = BH2(K, I)/BHSUM
                            end do
                        end if
                        do K = KTTR(JT), KBTR(JT)
                            QSS(K, I) = QSS(K, I) + QTR(JT)*QTRF(K, JT)
                        end do
                    end if
                end do
            end if
            if (DIST_TRIBS(JB)) then
                AKBR = 0.0
                do I = IU, ID
                    AKBR = AKBR + BI(KT, I)*DLX(I)
                end do
                do I = IU, ID
                    QDT(I) = QDTR(JB)*BI(KT, I)*DLX(I)/AKBR
                    QSS(KT, I) = QSS(KT, I) + QDT(I)
                end do
            end if
            if (WITHDRAWALS) then
                do JWD = 1, NWD
                    if (JB == JBWD(JWD)) then
                        I = MAX(CUS(JBWD(JWD)), IWD(JWD))
                        call LATERAL_WITHDRAWAL() !(JWD)
                        do K = KTW(JWD), KBW(JWD)
                            QSS(K, I) = QSS(K, I) - QSW(K, JWD)
                        end do
                    end if
                end do
            end if
            if (UH_INTERNAL(JB)) then
                if (UHS(JB) /= DS(JBUH(JB)) .or. DHS(JBUH(JB)) /= US(JB)) then
                    if (JBUH(JB) >= BS(JW) .and. JBUH(JB) <= BE(JW)) then
                        do K = KT, KB(IU - 1)
                            QSS(K, UHS(JB)) = QSS(K, UHS(JB)) - VOLUH2(K, JB)/DLT
                        end do
                    else
                        call UPSTREAM_FLOW()
                    end if
                end if
            end if
            if (DH_INTERNAL(JB)) then
                if (DHS(JB) /= US(JBDH(JB)) .or. UHS(JBDH(JB)) /= DS(JB)) then
                    if (JBDH(JB) >= BS(JW) .and. JBDH(JB) <= BE(JW)) then
                        do K = KT, KB(ID + 1)
                            QSS(K, CDHS(JB)) = QSS(K, CDHS(JB)) + VOLDH2(K, JB)/DLT
                        end do
                    else
                        call DOWNSTREAM_FLOW()
                    end if
                end if
            end if
        end do
    end do

! including tributary flows for inactive branches
    do JW = 1, NWB ! cb 11/20/19
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            if (BR_INACTIVE(JB)) then
                IU = CUS(JB)
                ID = DS(JB)
                if (TRIBUTARIES) then
                    do JT = 1, JTT

!********** Inflow fractions

                        if (JB == JBTR(JT)) then
                            I = cus(jbdn(jw)) ! placing tributary flows in upstream end of main branch
                            QTRF(KT:KB(I), JT) = 0.0
                            KTTR(JT) = KT
                            KBTR(JT) = KB(I)
                            KTTR(JT) = MAX(KT, KTTR(JT))
                            KBTR(JT) = MIN(KB(I), KBTR(JT))
                            if (KBTR(JT) < KTTR(JT)) then
                                KBTR(JT) = KTTR(JT)
                            end if
                            BHSUM = 0.0
                            do K = KTTR(JT), KBTR(JT)
                                BHSUM = BHSUM + BH2(K, I)
                            end do
                            do K = KTTR(JT), KBTR(JT)
                                QTRF(K, JT) = BH2(K, I)/BHSUM
                            end do
                            do K = KTTR(JT), KBTR(JT)
                                QSS(K, I) = QSS(K, I) + QTR(JT)*QTRF(K, JT)
                            end do
                        end if
                    end do
                end if
            end if
        end do
    end do

!** Compute tributary contribution to cross-shear

    if (TRIBUTARIES) then
        do JW = 1, NWB
            do JB = BS(JW), BE(JW)
                do JT = 1, JTT
                    if (JB == JBTR(JT)) then
                        I = MAX(CUS(JB), ITR(JT))
                        do K = KTWB(JW), KBMIN(I)
                            UYBR(K, I) = UYBR(K, I) + ABS(QTR(JT))*QTRF(K, JT)
                        end do
                    end if
                end do
            end do
        end do
    end if

    return
end subroutine HYDROINOUT
