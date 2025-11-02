subroutine BALANCES()

    use MAIN
    use GLOBAL;     use NAMESC;     use GEOMC;     use LOGICC;     use PREC;     use SURFHE;     use KINETIC;     use SHADEC;     use EDDY
    use STRUCTURES;     use TRANS;     use TVDC;     use SELWC;     use GDAYC;     use SCREENC;     use TDGAS;     use RSTART
    use MACROPHYTEC;     use POROSITYC;     use ZOOPLANKTONC
    use CEMAVars
    implicit none
    external :: RESTART_OUTPUT
    real :: VOLINJW, VOLPRJW, VOLOUTJW, VOLWDJW, VOLEVJW, VOLDTJW, VOLTRBJW, VOLICEJW, TPWB, TPSED, TNWB, TNSED, TPPLANT, TNPLANT

!***********************************************************************************************************************************
!*                                                    TASK 2.6: BALANCES                                                          **
!***********************************************************************************************************************************

!QINT  = 0.0
!QOUTT = 0.0
    VOLSR = 0.0
    VOLTR = 0.0

    do JW = 1, NWB
        VOLINJW = 0.0
        VOLPRJW = 0.0
        VOLOUTJW = 0.0
        VOLWDJW = 0.0
        VOLEVJW = 0.0
        VOLDTJW = 0.0
        VOLTRBJW = 0.0
        VOLICEJW = 0.0

        TPWB = 0.0
        TPSED = 0.0
        TNWB = 0.0
        TNSED = 0.0
        TNPLANT = 0.0
        TPPLANT = 0.0

        KT = KTWB(JW)
        if (VOLUME_BALANCE(JW)) then
            do JB = BS(JW), BE(JW)
                if (.not. BR_INACTIVE(JB)) then ! SW 8/8/2018 
                    VOLSBR(JB) = VOLSBR(JB) + DLVOL(JB)
                    VOLTBR(JB) = VOLEV(JB) + VOLPR(JB) + VOLTRB(JB) + VOLDT(JB) + VOLWD(JB) + VOLUH(JB) + VOLDH(JB) + VOLIN(JB) + VOLOUT(JB) + VOLICE(JB)
                    if (sediment_diagenesis) then
                        if (CEMARelatedCode .and. IncludeBedConsolidation) then
                            VOLTBR(JB) = VOLTBR(JB) + VOLCEMA(JB)
                        end if
                    end if
                    VOLSR(JW) = VOLSR(JW) + VOLSBR(JB)
                    VOLTR(JW) = VOLTR(JW) + VOLTBR(JB)
                end if
!QINT(JW)   = QINT(JW) +VOLIN(JB)+VOLTRB(JB)+VOLDT(JB)+VOLPR(JB)
!QOUTT(JW)  = QOUTT(JW)-VOLEV(JB)-VOLWD(JB) -VOLOUT(JB)
                if (ABS(VOLSBR(JB) - VOLTBR(JB)) > VTOL .and. VOLTBR(JB) > 100.0*VTOL) then
                    if (VOLUME_WARNING) then
                        write(WRN, "(A,F0.4,/A,I0,A,I0,A,I0,3(:/A,E15.8,A))") "COMPUTATIONAL WARNING AT JULIAN DAY = ", JDAY, "WATERBODY=", JW, ", BRANCH=", JB, ", KT=", KT, "SPATIAL CHANGE  =", VOLSBR(JB), " M^3", "TEMPORAL CHANGE =", VOLTBR(JB), " M^3", "VOLUME ERROR    =", VOLSBR(JB) - VOLTBR(JB), " M^3" !SR 11/16/19
                        write(WRN, *) "LAYER CHANGE:", LAYERCHANGE(JW)
                        write(WRN, *) "SZ", SZ(CUS(JB):DS(JB)), "Z", Z(CUS(JB):DS(JB)), "H2KT", H2(KT, CUS(JB):DS(JB)), "H1KT", H1(KT, CUS(JB):DS(JB)), "WSE", ELWS(CUS(JB):DS(JB)), "Q", Q(CUS(JB):DS(JB)), "QC", QC(CUS(JB):DS(JB)), "T1", T1(KT, CUS(JB):DS(JB)), "T2", T2(KT, CUS(JB):DS(JB)), "SUKT", SU(KT, CUS(JB):DS(JB)), "UKT", U(KT, CUS(JB):DS(JB)), "QIN", QINSUM(JB), "QTR", QTR, "QWD", QWD !SR 11/16/19
                        WARNING_OPEN = .true.
                        VOLUME_WARNING = .false.
                    end if
                end if
                if (VOLSR(JW) /= 0.0) then
                    DLVR(JW) = (VOLTR(JW) - VOLSR(JW))/VOLSR(JW)*100.0
                end if
                VOLINJW = VOLINJW + VOLIN(JB)
                VOLPRJW = VOLPRJW + VOLPR(JB)
                VOLOUTJW = VOLOUTJW + VOLOUT(JB)
                VOLWDJW = VOLWDJW + VOLWD(JB)
                VOLEVJW = VOLEVJW + VOLEV(JB)
                VOLDTJW = VOLDTJW + VOLDT(JB)
                VOLTRBJW = VOLTRBJW + VOLTRB(JB)
                VOLICEJW = VOLICEJW + VOLICE(JB)
            end do

            if (FLOWBALC == "      ON") then
                if (JDAY >= NXFLOWBAL) then
!NXFLOWBAL = NXFLOWBAL+FLOWBALF  

                    if (VOLUME_BALANCE(JW)) then
                        write(FLOWBFN, '(F10.3,",",1X,I3,",",11(E16.8,",",1X))') JDAY, JW, VOLINJW, VOLPRJW, VOLOUTJW, VOLWDJW, VOLEVJW, VOLDTJW, VOLTRBJW, VOLICEJW, DLVR(JW)
                    else
                        write(FLOWBFN, '(F10.3,",",1X,I3,",",10(E16.8,",",1X))') JDAY, JW, VOLINJW, VOLPRJW, VOLOUTJW, VOLWDJW, VOLEVJW, VOLDTJW, VOLTRBJW, VOLICEJW
                    end if
                end if

            end if ! CONTOUR INTERVAL FOR WRITING OUT FLOW BALANCE 
        end if ! VOLUME BALANCE
        if (ENERGY_BALANCE(JW)) then
            ESR(JW) = 0.0
            ETR(JW) = 0.0
            do JB = BS(JW), BE(JW)
                if (BR_INACTIVE(JB)) then
                    cycle
                end if ! SW 8/8/2018 
                ETBR(JB) = EBRI(JB) + TSSEV(JB) + TSSPR(JB) + TSSTR(JB) + TSSDT(JB) + TSSWD(JB) + TSSUH(JB) + TSSDH(JB) + TSSIN(JB) + TSSOUT(JB) + TSSS(JB) + TSSB(JB) + TSSICE(JB)
                ESBR(JB) = 0.0
                do I = CUS(JB), DS(JB)
                    do K = KT, KB(I)
                        ESBR(JB) = ESBR(JB) + T1(K, I)*DLX(I)*BH1(K, I)
                    end do
                end do
                ETR(JW) = ETR(JW) + ETBR(JB)
                ESR(JW) = ESR(JW) + ESBR(JB)
            end do
        end if
        if (MASS_BALANCE(JW)) then
            do JB = BS(JW), BE(JW)
                if (BR_INACTIVE(JB)) then
                    cycle
                end if ! SW 8/8/2018 
                do JC = 1, NAC
                    CMBRS(CN(JC), JB) = 0.0
                    do I = CUS(JB), DS(JB)
                        do K = KT, KB(I)
                            CMBRS(CN(JC), JB) = CMBRS(CN(JC), JB) + C1(K, I, CN(JC))*DLX(I)*BH1(K, I)
                            CMBRT(CN(JC), JB) = CMBRT(CN(JC), JB) + (CSSB(K, I, CN(JC)) + CSSK(K, I, CN(JC))*BH1(K, I)*DLX(I))*DLT
                        end do
                    end do
                end do
                if (DERIVED_CALC) then
                    do I = CUS(JB), DS(JB)
                        do K = KT, KB(I)
                            TPWB = TPWB + TP(K, I)*VOL(K, I)*0.001 !/1000.   ! kg
                            TPSED = TPSED + SEDP(K, I)*VOL(K, I)*0.001 !/1000.   ! kg        
                            TNWB = TNWB + TN(K, I)*VOL(K, I)*0.001 !/1000.  ! kg
                            TNSED = TNSED + SEDN(K, I)*VOL(K, I)*0.001 !/1000.  ! kg     
                            PFLUXIN(JW) = PFLUXIN(JW) + SEDPINFLUX(K, I)*VOL(K, I)*0.001 !/1000.   ! kg
                            NFLUXIN(JW) = NFLUXIN(JW) + SEDNINFLUX(K, I)*VOL(K, I)*0.001 !/1000.   ! kg
                            do M = 1, NMC
                                TNPLANT = TNPLANT + MAC(K, I, M)*VOL(K, I)*MN(M)*0.001 !/1000.
                                TPPLANT = TPPLANT + MAC(K, I, M)*VOL(K, I)*MP(M)*0.001 !/1000.
                            end do
                            do M = 1, NEP
                                TNPLANT = TNPLANT + EPM(K, I, M)*EN(M)*0.001 !/1000.
                                TPPLANT = TPPLANT + EPM(K, I, M)*EP(M)*0.001 !/1000.
                            end do
                        end do
                    end do
                end if

! MACROPHYTES
                do M = 1, NMC
                    if (MACROPHYTE_CALC(JW, M)) then
                        MACMBRS(JB, M) = 0.0
                        do I = CUS(JB), DS(JB)
                            if (KTICOL(I)) then
                                JT = KTI(I)
                            else
                                JT = KTI(I) + 1
                            end if
                            JE = KB(I)
                            do J = JT, JE
                                if (J < KT) then
                                    COLB = EL(J + 1, I)
                                else
                                    COLB = EL(KT + 1, I)
                                end if
!COLDEP=ELWS(I)-COLB
                                coldep = EL(KT, i) - Z(i)*COSA(JB) - colb ! cb 3/7/16
                                MACMBRS(JB, M) = MACMBRS(JB, M) + MACRM(J, KT, I, M)
                                MACMBRT(JB, M) = MACMBRT(JB, M) + MACSS(J, KT, I, M)*COLDEP*CW(J, I)*DLX(I)*DLT
                            end do
                            do K = KT + 1, KB(I)
                                JT = K
                                JE = KB(I)
                                do J = JT, JE
                                    MACMBRS(JB, M) = MACMBRS(JB, M) + MACRM(J, K, I, M)
!                    MACMBRT(JB,M) = MACMBRT(JB,M)+(MACSS(J,K,I,M)*H2(K,I)*CW(J,I)*DLX(I))*DLT
                                    MACMBRT(JB, M) = MACMBRT(JB, M) + MACSS(J, K, I, M)*CW(J, I)/B(K, I)*BH1(K, I)*DLX(I)*DLT
                                end do
                            end do
                        end do
                    end if
                end do
! END MACROPHYTES
            end do

            if (NPBALC == "      ON") then
                if (JDAY >= NXNPBAL) then
                    if (SEDIMENT_DIAGENESIS) then
                        write(MASSBFN, '(F10.3,",",1X,I3,",",41(E16.8,",",1X))') JDAY, JW, TPWB, TPSED, TPPLANT, TPOUT(JW), TPTRIB(JW), TPDTRIB(JW), TPWD(JW), TPPR(JW), TPIN(JW), ATMDEP_P(JW), TP_SEDSOD_PO4(JW), PFLUXIN(JW), SDPFLUX(JW), TNWB, TNSED, TNPLANT, TNOUT(JW), TNTRIB(JW), TNDTRIB(JW), TNWD(JW), TNPR(JW), TNIN(JW), ATMDEP_N(JW), NH3GASLOSS(JW), TN_SEDSOD_NH4(JW), NFLUXIN(JW), SDNH4FLUX(JW), SDNO3FLUX(JW) ! TP_SEDBURIAL(JW),TN_SEDBURIAL(JW),
                    else
                        write(MASSBFN, '(F10.3,",",1X,I3,",",31(E16.8,",",1X))') JDAY, JW, TPWB, TPSED, TPPLANT, TPOUT(JW), TPTRIB(JW), TPDTRIB(JW), TPWD(JW), TPPR(JW), TPIN(JW), ATMDEP_P(JW), TP_SEDSOD_PO4(JW), PFLUXIN(JW), TNWB, TNSED, TNPLANT, TNOUT(JW), TNTRIB(JW), TNDTRIB(JW), TNWD(JW), TNPR(JW), TNIN(JW), ATMDEP_N(JW), NH3GASLOSS(JW), TN_SEDSOD_NH4(JW), NFLUXIN(JW) ! TP_SEDBURIAL(JW),TN_SEDBURIAL(JW),
                    end if
                end if

            end if ! CONTOUR INTERVAL FOR WRITING OUT FLOW BALANCE 


        end if ! MASS BALANCE
    end do

    if (JDAY >= NXNPBAL) then
        NXNPBAL = NXNPBAL + NPBALF
    end if
    if (FLOWBALC == "      ON") then ! cb 8/22/21
        if (JDAY >= NXFLOWBAL) then
            NXFLOWBAL = NXFLOWBAL + FLOWBALF
        end if
    end if

    return
end subroutine BALANCES
