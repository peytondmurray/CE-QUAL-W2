subroutine WQCONSTITUENTS()

    use MAIN
    use GLOBAL;     use NAMESC;     use GEOMC;     use LOGICC;     use PREC;     use SURFHE;     use KINETIC;     use SHADEC;     use EDDY
    use STRUCTURES;     use TRANS;     use TVDC;     use SELWC;     use GDAYC;     use SCREENC;     use TDGAS;     use RSTART
    use MACROPHYTEC;     use POROSITYC;     use ZOOPLANKTONC;     use TRIDIAG_V
    use CEMAVars
    use CEMASedimentDiagenesis, only: SedimentFlux;     use ALGAE_TOXINS

    implicit none
    external :: RESTART_OUTPUT
    real :: TPALG, TNALG, TPZ, TNZ, TPBOD, TNBOD

    if (MACROPHYTE_ON .and. UPDATE_KINETICS) then
        call POROSITY()
    end if
    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            if (BR_INACTIVE(JB)) then
                cycle
            end if ! SW 6/12/2017
            IU = CUS(JB)
            ID = DS(JB)

!******** Kinetic sources/sinks

            if (SEDIMENT_CALC(JW)) then
                call SEDIMENT()
                call SEDIMENTP()
                call SEDIMENTN()
                call SEDIMENTC()
                if (DYNSEDK(JW) == "      ON") then
                    call SEDIMENT_DECAY_RATE()
                end if
            end if
! Standing biomass decay
            if (STANDING_BIOMASS_DECAY) then ! SW 5/26/15
                if (SEDIMENT_CALC1(JW)) then
                    call SEDIMENT1()
                end if
                if (SEDIMENT_CALC2(JW)) then
                    call SEDIMENT2()
                end if
            end if

            do M = 1, NMC
                if (MACROPHYTE_CALC(JW, M)) then
                    call MACROPHYTE(M)
                end if
            end do

            if (UPDATE_KINETICS) then
                if (UPDATE_RATES) then
                    call TEMPERATURE_RATES()
                    call KINETIC_RATES()
                end if
                do JAC = 1, NAC
                    JC = CN(JAC)
                    if (JC == NPO4) then
                        call PHOSPHORUS()
                    end if
                    if (JC == NWAGE) then
                        call WATER_AGE()
                    end if
                    if (JC == NBACT) then
                        call BACTERIA()
                    end if
                    if (JC == NDGP) then
                        call DISSOLVED_GAS()
                    end if
                    if (JC == NN2) then
                        call DISSOLVED_N2()
                    end if
                    if (JC == NH2S) then
                        call SULFIDE()
                    end if
                    if (JC == NCH4) then
                        call METHANE()
                    end if
                    if (JC == NSO4) then
                        call SULFATE()
                    end if
                    if (JC == NFEII) then
                        call FERROUS()
                    end if
                    if (JC == NFEOOH) then
                        call OXIDIZEDFE()
                    end if
                    if (JC == NMNII) then
                        call BIVALENTMN()
                    end if
                    if (JC == NMNO2) then
                        call OXIDIZEDMN()
                    end if

                    if (JC == NNH4) then
                        call AMMONIUM()
                    end if
                    if (JC == NNO3) then
                        call NITRATE()
                    end if
                    if (JC == NDSI) then
                        call DISSOLVED_SILICA()
                    end if
                    if (JC == NPSI) then
                        call PARTICULATE_SILICA()
                    end if
!
                    if (ORGC_CALC) then
                        if (JC == NLDOMC) then
                            call LABILE_DOM_C()
                        end if
                        if (JC == NRDOMC) then
                            call REFRACTORY_DOM_C()
                        end if
                        if (JC == NLPOMC) then
                            call LABILE_POM_C()
                        end if
                        if (JC == NRPOMC) then
                            call REFRACTORY_POM_C()
                        end if
                    else
                        if (JC == NLDOM) then
                            call LABILE_DOM()
                        end if
                        if (JC == NRDOM) then
                            call REFRACTORY_DOM()
                        end if
                        if (JC == NLPOM) then
                            call LABILE_POM()
                        end if
                        if (JC == NRPOM) then
                            call REFRACTORY_POM()
                        end if
                    end if
                    if (JC == NDO) then
                        call DISSOLVED_OXYGEN()
                    end if
                    if (JC >= NGCS .and. JC <= NGCE) then
                        call GENERIC_CONST(JC - NGCS + 1)
                    end if
                    if (JC >= NSSS .and. JC <= NSSE) then
                        call SUSPENDED_SOLIDS(JC - NSSS + 1)
                    end if
                    if (JC >= NAS .and. JC <= NAE) then
                        if (ALG_CALC(JC - NAS + 1)) then
                            call ALGAE(JC - NAS + 1)
                        end if
                    end if
                    if (JC >= NBODS .and. JC <= NBODE) then
                        do JCB = 1, NBOD ! VARIABLE STOICHIOMETRY FOR CBOD, CB 6/6/10
                            if (BOD_CALC(JCB)) then
                                if (JC == NBODC(JCB)) then
                                    call BIOCHEMICAL_O2_DEMAND(JCB)
                                end if
                                if (JC == NBODP(JCB) .and. BOD_CALCP(JCB)) then
                                    call BIOCHEMICAL_O2_DEMAND_P(JCB)
                                end if ! CB 5/19/2011
                                if (JC == NBODN(JCB) .and. BOD_CALCN(JCB)) then
                                    call BIOCHEMICAL_O2_DEMAND_N(JCB)
                                end if ! CB 5/19/2011
                            end if
                        end do
                    end if
                    if (JC >= NZOOS .and. JC <= NZOOE .and. ZOOPLANKTON_CALC) then
                        call ZOOPLANKTON()
                    end if
                    if (JC == NLDOMP) then
                        call LABILE_DOM_P()
                    end if
                    if (JC == NRDOMP) then
                        call REFRACTORY_DOM_P()
                    end if
                    if (JC == NLPOMP) then
                        call LABILE_POM_P()
                    end if
                    if (JC == NRPOMP) then
                        call REFRACTORY_POM_P()
                    end if
                    if (JC == NLDOMN) then
                        call LABILE_DOM_N()
                    end if
                    if (JC == NRDOMN) then
                        call REFRACTORY_DOM_N()
                    end if
                    if (JC == NLPOMN) then
                        call LABILE_POM_N()
                    end if
                    if (JC == NRPOMN) then
                        call REFRACTORY_POM_N()
                    end if
!IF (JC == NALK .and. NONCON_ALKALINITY)                  CALL alkalinity
                    if (JC == NALK) then
                        call alkalinity()
                    end if ! NW 2/11/16
                    if (JC >= NATS .and. JC <= NATE .and. ALGAE_TOXIN) then
                        call INTRACELLULAR_TOXIN(JC - NATS + 1)
                        call EXTRACELLULAR_TOXIN(JC - NATS + 1)
                    end if

                end do
                if (PH_CALC(JW)) then
                    call INORGANIC_CARBON()
                end if
                if (PH_CALC(JW)) then
                    if (ph_buffering) then ! enhanced pH buffering                
                        call pH_CO2_new()
                    else
                        call PH_CO2()
                    end if
                end if
                if (CEMARelatedCode .and. IncludeCEMASedDiagenesis) then
                    call SedimentFlux()
                end if

            end if
            do JE = 1, NEP ! sw 5/16/06
                if (EPIPHYTON_CALC(JW, JE)) then
                    call EPIPHYTON(JE)
                end if
            end do

!******** External sources/sinks

            if (AERATEC == "      ON") then
                call AERATEMASS()
            end if
            if (EVAPORATION(JW) .and. WATER_AGE_ACTIVE) then ! CORRECT WATER AGE FOR EVAPORATION SR 7/27/2017
                do I = IU, ID
!JC=NGCS+JG_AGE-1
                    CSSB(KT, I, NWAGE) = CSSB(KT, I, NWAGE) - EV(I)*WAGE(KT, I) ! SW 10/17/2019 CG(KT,I,JC)           
                end do
            end if

            do JAC = 1, NAC
                JC = CN(JAC)
                if (TRIBUTARIES) then
                    do JT = 1, JTT
                        if (JB == JBTR(JT)) then
                            I = ITR(JT)
                            if (I < CUS(JB)) then
                                I = CUS(JB)
                            end if
                            do K = KTTR(JT), KBTR(JT)
                                if (QTR(JT) < 0.0) then
                                    CSSB(K, I, JC) = CSSB(K, I, JC) + C1(K, I, JC)*QTR(JT)*QTRF(K, JT)
                                else
                                    CSSB(K, I, JC) = CSSB(K, I, JC) + CTR(JC, JT)*QTR(JT)*QTRF(K, JT)
                                end if
                            end do
                        end if
                    end do
                end if
                if (DIST_TRIBS(JB)) then
                    do I = IU, ID
                        if (QDT(I) < 0.0) then
                            CSSB(KT, I, JC) = CSSB(KT, I, JC) + C1(KT, I, JC)*QDT(I)
                        else
                            CSSB(KT, I, JC) = CSSB(KT, I, JC) + CDTR(JC, JB)*QDT(I)
                        end if
                    end do
                end if
                if (WITHDRAWALS) then
                    do JWD = 1, JWW
                        if (QWD(JWD) /= 0.0) then
                            if (JB == JBWD(JWD)) then
                                I = MAX(CUS(JBWD(JWD)), IWD(JWD))
                                do K = KTW(JWD), KBW(JWD) !CONCURRENT(K=KTW(JWD):KBW(JWD))                 ! FORALL
                                    CSSB(K, I, JC) = CSSB(K, I, JC) - C1S(K, I, JC)*QSW(K, JWD)
                                end do
                            end if
                        end if
                    end do
                end if
                if (PRECIPITATION(JW)) then
                    do I = IU, ID !CONCURRENT (I=IU:ID)                                  !FORALL
                        CSSB(KT, I, JC) = CSSB(KT, I, JC) + CPR(JC, JB)*QPR(I)
                    end do
                end if
                if (UP_FLOW(JB)) then
                    do K = KT, KB(IU)
                        if (.not. HEAD_FLOW(JB)) then
                            CSSB(K, IU, JC) = CSSB(K, IU, JC) + QINF(K, JB)*QIN(JB)*CIN(JC, JB)
                        else
                            if (U(K, IU - 1) >= 0.0) then
                                CSSB(K, IU, JC) = CSSB(K, IU, JC) + U(K, IU - 1)*BHR1(K, IU - 1)*C1S(K, IU - 1, JC)
                            else
                                CSSB(K, IU, JC) = CSSB(K, IU, JC) + U(K, IU - 1)*BHR1(K, IU - 1)*C1S(K, IU, JC)
                            end if
                        end if
                    end do
                end if
                if (DN_FLOW(JB)) then
                    CSSB(KT:KB(ID), ID, JC) = CSSB(KT:KB(ID), ID, JC) - QOUT(KT:KB(ID), JB)*C1S(KT:KB(ID), ID, JC)
                end if
                if (UP_HEAD(JB)) then
                    do K = KT, KB(IU)
                        IUT = IU
                        if (QUH1(K, JB) >= 0.0) then
                            IUT = IU - 1
                        end if
                        CSSUH1(K, JC, JB) = C1S(K, IUT, JC)*QUH1(K, JB)
                        CSSB(K, IU, JC) = CSSB(K, IU, JC) + CSSUH1(K, JC, JB)
                    end do
                    if (UH_INTERNAL(JB)) then
                        if (UHS(JB) /= DS(JBUH(JB)) .or. DHS(JBUH(JB)) /= US(JB)) then
                            if (JBUH(JB) >= BS(JW) .and. JBUH(JB) <= BE(JW)) then
                                I = UHS(JB)
!dir$ ivdep
                                do K = KT, KB(IU)
                                    CSSB(K, I, JC) = CSSB(K, I, JC) - CSSUH2(K, JC, JB)/DLT
                                end do
                            else
                                call UPSTREAM_CONSTITUENT(C2(:, :, JC), CSSB(:, :, JC))
                            end if
                        end if
                    end if
                end if
                if (DN_HEAD(JB)) then
                    do K = KT, KB(ID + 1)
                        IDT = ID + 1
                        if (QDH1(K, JB) >= 0.0) then
                            IDT = ID
                        end if
                        CSSDH1(K, JC, JB) = C1S(K, IDT, JC)*QDH1(K, JB)
                        CSSB(K, ID, JC) = CSSB(K, ID, JC) - CSSDH1(K, JC, JB)
                    end do
                    if (DH_INTERNAL(JB)) then
                        if (DHS(JB) /= US(JBDH(JB)) .or. UHS(JBDH(JB)) /= DS(JB)) then
                            if (JBDH(JB) >= BS(JW) .and. JBDH(JB) <= BE(JW)) then
                                I = DHS(JB)
                                do K = KT, KB(ID + 1)
                                    CSSB(K, I, JC) = CSSB(K, I, JC) + CSSDH2(K, JC, JB)/DLT
                                end do
                            else
                                call DOWNSTREAM_CONSTITUENT(C2(:, :, JC), CSSB(:, :, JC))
                            end if
                        end if
                    end if
                end if
            end do

            if (NPBALC == "      ON") then !IF(MASS_BALANCE(JW).AND.CONTOUR(JW).AND.DERIVED_CALC)THEN     ! TO COMPUTE TP AND TN INFLOWS AND OUTFLOWS FOR MASSBAL.OPT FILE
                if (TRIBUTARIES) then
                    do JT = 1, JTT
                        if (JB == JBTR(JT)) then
                            I = ITR(JT)
                            if (I < CUS(JB)) then
                                I = CUS(JB)
                            end if
                            do K = KTTR(JT), KBTR(JT)
                                TPALG = 0.0;                                 TNALG = 0.0;                                 TPZ = 0.0;                                 TNZ = 0.0;                                 TPBOD = 0.0;                                 TNBOD = 0.0
                                if (QTR(JT) < 0.0) then
                                    do JA = 1, NAL
                                        TPALG = TPALG + ALG(K, I, JA)*AP(JA)
                                        TNALG = TNALG + ALG(K, I, JA)*AN(JA)
                                    end do
                                    do JCB = 1, NBOD
                                        TPBOD = TPBOD + CBODP(K, I, JCB)
                                        TNBOD = TNBOD + CBODN(K, I, JCB)
                                    end do
                                    do JZ = 1, NZP
                                        TPZ = TPZ + ZOO(K, I, JZ)*ZP(JZ)
                                        TNZ = TNZ + ZOO(K, I, JZ)*ZN(JZ)
                                    end do
                                    TPOUT = TPOUT - (TPALG + TPBOD + TPZ + PO4(K, I) + LDOMP(K, I) + LPOMP(K, I) + RDOMP(K, I) + RPOMP(K, I))*QTR(JT)*QTRF(K, JT)*DLT/1000.
                                    TNOUT = TNOUT - (TNALG + TNBOD + TNZ + NO3(K, I) + NH4(K, I) + LDOMN(K, I) + LPOMN(K, I) + RDOMN(K, I) + RPOMN(K, I))*QTR(JT)*QTRF(K, JT)*DLT/1000.
                                else
                                    do JC = NAS, NAE
                                        TPALG = TPALG + CTR(JC, JT)*AP(JC - NAS + 1)
                                        TNALG = TNALG + CTR(JC, JT)*AN(JC - NAS + 1)
                                    end do
                                    do JC = NBODS, NBODE, 3
                                        TPBOD = TPBOD + CTR(JC + 1, JT)
                                        TNBOD = TNBOD + CTR(JC + 2, JT)
                                    end do
                                    do JC = NZOOS, NZOOE
                                        TPZ = TPZ + CTR(JC, JT)*ZP(JC - NZOOS + 1)
                                        TNZ = TNZ + CTR(JC, JT)*ZN(JC - NZOOS + 1)
                                    end do
                                    TPTRIB(JW) = TPTRIB(JW) + (TPALG + TPBOD + TPZ + CTR(NPO4, JT) + CTR(NLDOMP, JT) + CTR(NRDOMP, JT) + CTR(NLPOMP, JT) + CTR(NRPOMP, JT))*QTR(JT)*QTRF(K, JT)*DLT/1000.
                                    TNTRIB(JW) = TNTRIB(JW) + (TNALG + TNBOD + TNZ + CTR(NNH4, JT) + CTR(NNO3, JT) + CTR(NLDOMN, JT) + CTR(NRDOMN, JT) + CTR(NLPOMN, JT) + CTR(NRPOMN, JT))*QTR(JT)*QTRF(K, JT)*DLT/1000.
                                end if
                            end do
                        end if
                    end do
                end if
                if (DIST_TRIBS(JB)) then
                    do I = IU, ID
                        TPALG = 0.0;                         TNALG = 0.0;                         TPZ = 0.0;                         TNZ = 0.0;                         TPBOD = 0.0;                         TNBOD = 0.0
                        if (QDT(I) < 0.0) then
                            do JA = 1, NAL
                                TPALG = TPALG + ALG(KT, I, JA)*AP(JA)
                                TNALG = TNALG + ALG(KT, I, JA)*AN(JA)
                            end do
                            do JCB = 1, NBOD
                                TPBOD = TPBOD + CBODP(KT, I, JCB)
                                TNBOD = TNBOD + CBODN(KT, I, JCB)
                            end do
                            do JZ = 1, NZP
                                TPZ = TPZ + ZOO(KT, I, JZ)*ZP(JZ)
                                TNZ = TNZ + ZOO(KT, I, JZ)*ZN(JZ)
                            end do
                            TPOUT = TPOUT - (TPALG + TPBOD + TPZ + PO4(KT, I) + LDOMP(KT, I) + LPOMP(KT, I) + RDOMP(KT, I) + RPOMP(KT, I))*QDT(I)*DLT/1000.
                            TNOUT = TNOUT - (TNALG + TNBOD + TNZ + NO3(KT, I) + NH4(KT, I) + LDOMN(KT, I) + LPOMN(KT, I) + RDOMN(KT, I) + RPOMN(KT, I))*QDT(I)*DLT/1000.
                        else
                            do JC = NAS, NAE
                                TPALG = TPALG + CDTR(JC, JB)*AP(JC - NAS + 1)
                                TNALG = TNALG + CDTR(JC, JB)*AN(JC - NAS + 1)
                            end do
                            do JC = NBODS, NBODE, 3
                                TPBOD = TPBOD + CDTR(JC + 1, JB)
                                TNBOD = TNBOD + CDTR(JC + 2, JB)
                            end do
                            do JC = NZOOS, NZOOE
                                TPZ = TPZ + CDTR(JC, JB)*ZP(JC - NZOOS + 1)
                                TNZ = TNZ + CDTR(JC, JB)*ZN(JC - NZOOS + 1)
                            end do
                            TPDTRIB(JW) = TPDTRIB(JW) + (TPALG + TPBOD + TPZ + CDTR(NPO4, JB) + CDTR(NLDOMP, JB) + CDTR(NRDOMP, JB) + CDTR(NLPOMP, JB) + CDTR(NRPOMP, JB))*QDT(I)*DLT/1000.
                            TNDTRIB(JW) = TNDTRIB(JW) + (TNALG + TNBOD + TNZ + CDTR(NNH4, JB) + CDTR(NNO3, JB) + CDTR(NLDOMN, JB) + CDTR(NRDOMN, JB) + CDTR(NLPOMN, JB) + CDTR(NRPOMN, JB))*QDT(I)*DLT/1000.
                        end if
                    end do
                end if
                if (WITHDRAWALS) then
                    do JWD = 1, JWW
                        if (QWD(JWD) /= 0.0) then
                            if (JB == JBWD(JWD)) then
                                I = MAX(CUS(JBWD(JWD)), IWD(JWD))
                                do K = KTW(JWD), KBW(JWD) !CONCURRENT(K=KTW(JWD):KBW(JWD))                 ! FORALL
                                    TPALG = 0.0;                                     TNALG = 0.0;                                     TPZ = 0.0;                                     TNZ = 0.0;                                     TPBOD = 0.0;                                     TNBOD = 0.0
                                    do JC = NAS, NAE
                                        TPALG = TPALG + C1S(K, I, JC)*AP(JC - NAS + 1)
                                        TNALG = TNALG + C1S(K, I, JC)*AN(JC - NAS + 1)
                                    end do
                                    do JC = NBODS, NBODE, 3
                                        TPBOD = TPBOD + C1S(K, I, JC)
                                        TNBOD = TNBOD + C1S(K, I, JC)
                                    end do
                                    do JC = NZOOS, NZOOE
                                        TPZ = TPZ + C1S(K, I, JC)*ZP(JC - NZOOS + 1)
                                        TNZ = TNZ + C1S(K, I, JC)*ZN(JC - NZOOS + 1)
                                    end do
                                    TPWD(JW) = TPWD(JW) + (TPALG + TPBOD + TPZ + C1S(K, I, NPO4) + C1S(K, I, NLDOMP) + C1S(K, I, NRDOMP) + C1S(K, I, NLPOMP) + C1S(K, I, NRPOMP))*QSW(K, JWD)*DLT/1000.
                                    TNWD(JW) = TNWD(JW) + (TNALG + TNBOD + TNZ + C1S(K, I, NNH4) + C1S(K, I, NNO3) + C1S(K, I, NLDOMN) + C1S(K, I, NRDOMN) + C1S(K, I, NLPOMN) + C1S(K, I, NRPOMN))*QSW(K, JWD)*DLT/1000.
                                end do
                            end if
                        end if
                    end do
                end if
                if (PRECIPITATION(JW)) then
                    do I = IU, ID !CONCURRENT (I=IU:ID)                                  !FORALL
                        TPALG = 0.0;                         TNALG = 0.0;                         TPZ = 0.0;                         TNZ = 0.0;                         TPBOD = 0.0;                         TNBOD = 0.0
                        do JC = NAS, NAE
                            TPALG = TPALG + CPR(JC, JB)*AP(JC - NAS + 1)
                            TNALG = TNALG + CPR(JC, JB)*AN(JC - NAS + 1)
                        end do
                        do JC = NBODS, NBODE, 3
                            TPBOD = TPBOD + CPR(JC, JB)
                            TNBOD = TNBOD + CPR(JC, JB)
                        end do
                        do JC = NZOOS, NZOOE
                            TPZ = TPZ + CPR(JC, JB)*ZP(JC - NZOOS + 1)
                            TNZ = TNZ + CPR(JC, JB)*ZN(JC - NZOOS + 1)
                        end do
                        TPPR(JW) = TPPR(JW) + (TPALG + TPBOD + TPZ + CPR(NPO4, JB) + CPR(NLDOMP, JB) + CPR(NRDOMP, JB) + CPR(NLPOMP, JB) + CPR(NRPOMP, JB))*QPR(I)*DLT/1000.
                        TNPR(JW) = TNPR(JW) + (TNALG + TNBOD + TNZ + CPR(NNH4, JB) + CPR(NNO3, JB) + CPR(NLDOMN, JB) + CPR(NRDOMN, JB) + CPR(NLPOMN, JB) + CPR(NRPOMN, JB))*QPR(I)*DLT/1000.
                    end do
                end if
                if (UP_FLOW(JB)) then
                    do K = KT, KB(IU)
                        TPALG = 0.0;                         TNALG = 0.0;                         TPZ = 0.0;                         TNZ = 0.0;                         TPBOD = 0.0;                         TNBOD = 0.0
                        if (.not. HEAD_FLOW(JB)) then
                            do JC = NAS, NAE
                                TPALG = TPALG + CIN(JC, JB)*AP(JC - NAS + 1)
                                TNALG = TNALG + CIN(JC, JB)*AN(JC - NAS + 1)
                            end do
                            do JC = NBODS, NBODE, 3
                                TPBOD = TPBOD + CIN(JC, JB)
                                TNBOD = TNBOD + CIN(JC, JB)
                            end do
                            do JC = NZOOS, NZOOE
                                TPZ = TPZ + CIN(JC, JB)*ZP(JC - NZOOS + 1)
                                TNZ = TNZ + CIN(JC, JB)*ZN(JC - NZOOS + 1)
                            end do
                            TPIN(JW) = TPIN(JW) + (TPALG + TPBOD + TPZ + CIN(NPO4, JB) + CIN(NLDOMP, JB) + CIN(NRDOMP, JB) + CIN(NLPOMP, JB) + CIN(NRPOMP, JB))*QINF(K, JB)*QIN(JB)*DLT/1000.
                            TNIN(JW) = TNIN(JW) + (TNALG + TNBOD + TNZ + CIN(NNH4, JB) + CIN(NNO3, JB) + CIN(NLDOMN, JB) + CIN(NRDOMN, JB) + CIN(NLPOMN, JB) + CIN(NRPOMN, JB))*QINF(K, JB)*QIN(JB)*DLT/1000.
                        end if
                    end do
                end if
                if (DN_FLOW(JB)) then
                    do K = KT, KB(ID)
                        TPALG = 0.0;                         TNALG = 0.0;                         TPZ = 0.0;                         TNZ = 0.0;                         TPBOD = 0.0;                         TNBOD = 0.0
                        do JC = NAS, NAE
                            TPALG = TPALG + C1S(K, ID, JC)*AP(JC - NAS + 1)
                            TNALG = TNALG + C1S(K, ID, JC)*AN(JC - NAS + 1)
                        end do
                        do JC = NBODS, NBODE, 3
                            TPBOD = TPBOD + C1S(K, ID, JC)
                            TNBOD = TNBOD + C1S(K, ID, JC)
                        end do
                        do JC = NZOOS, NZOOE
                            TPZ = TPZ + C1S(K, ID, JC)*ZP(JC - NZOOS + 1)
                            TNZ = TNZ + C1S(K, ID, JC)*ZN(JC - NZOOS + 1)
                        end do
                        TPOUT(JW) = TPOUT(JW) + (TPALG + TPBOD + TPZ + C1S(K, ID, NPO4) + C1S(K, ID, NLDOMP) + C1S(K, ID, NRDOMP) + C1S(K, ID, NLPOMP) + C1S(K, ID, NRPOMP))*QOUT(K, JB)*DLT/1000. ! C1S(KT:KB(ID),ID,JC)
                        TNOUT(JW) = TNOUT(JW) + (TNALG + TNBOD + TNZ + C1S(K, ID, NNH4) + C1S(K, ID, NNO3) + C1S(K, ID, NLDOMN) + C1S(K, ID, NRDOMN) + C1S(K, ID, NLPOMN) + C1S(K, ID, NRPOMN))*QOUT(K, JB)*DLT/1000.
                    end do
                end if

!IF (UP_HEAD(JB)) THEN
!    DO K=KT,KB(IU)
!    IUT = IU
!    IF (QUH1(K,JB) >= 0.0) IUT = IU-1
!    !CSSUH1(K,JC,JB) = C1S(K,IUT,JC)*QUH1(K,JB)
!
!  END DO
!  IF (UH_INTERNAL(JB)) THEN
!    IF (UHS(JB) /= DS(JBUH(JB)) .OR. DHS(JBUH(JB)) /= US(JB)) THEN
!      IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
!        I = UHS(JB)
!        DO K=KT,KB(IU)
!          !CSSB(K,I,JC) = CSSB(K,I,JC)-CSSUH2(K,JC,JB)/DLT
!        END DO
!      ELSE
!        !CALL UPSTREAM_CONSTITUENT(C2(:,:,JC),CSSB(:,:,JC))
!      END IF
!    END IF
!  END IF
!END IF
!IF (DN_HEAD(JB)) THEN
!  DO K=KT,KB(ID+1)
!    IDT = ID+1
!    IF (QDH1(K,JB) >= 0.0) IDT = ID
!    !CSSDH1(K,JC,JB) = C1S(K,IDT,JC)*QDH1(K,JB)
!
!  END DO
!  IF (DH_INTERNAL(JB)) THEN
!    IF (DHS(JB) /= US(JBDH(JB)) .OR. UHS(JBDH(JB)) /= DS(JB)) THEN
!      IF (JBDH(JB) >= BS(JW) .AND. JBDH(JB) <= BE(JW)) THEN
!        I = DHS(JB)
!        DO K=KT,KB(ID+1)
!          !CSSB(K,I,JC) = CSSB(K,I,JC)+CSSDH2(K,JC,JB)/DLT
!        END DO
!      ELSE
!        !CALL DOWNSTREAM_CONSTITUENT(C2(:,:,JC),CSSB(:,:,JC))
!      END IF
!    END IF
!  END IF
!END IF
            end if ! END OF TP AND TN MASS BALANCES   


        end do ! JB loop

! Atmospheric Depsition original unit kg/km2/year
        if (ATM_DEPOSITION(JW)) then
            do JB = BS(JW), BE(JW)
                do I = CUS(JB), DS(JB)
                    do JAC = 1, NACATD(JW)
                        CSSB(KT, I, ATMDCN(JAC, JW)) = CSSB(KT, I, ATMDCN(JAC, JW)) + ATM_DEP_LOADING(ATMDCN(JAC, JW), JW)*BI(KT, I)*DLX(I)*3.17098E-11 ! Conversion: kg/km2/year to g/m2/s 1000/(365*86400*1000*1000)=3.17098E-11
                        if (ATMDCN(JAC, JW) == NPO4 .or. ATMDCN(JAC, JW) == NLPOMP .or. ATMDCN(JAC, JW) == NRPOMP) then
                            ATMDEP_P(JW) = ATMDEP_P(JW) + ATM_DEP_LOADING(ATMDCN(JAC, JW), JW)*BI(KT, I)*DLX(I)*3.17098E-11*DLT/1000. ! P MASS BALANCE IN KG - CUMULATIVE
                        else
                            if (ATMDCN(JAC, JW) == NNO3 .or. ATMDCN(JAC, JW) == NLPOMN .or. ATMDCN(JAC, JW) == NRPOMN .or. ATMDCN(JAC, JW) == NNH4) then
                                ATMDEP_N(JW) = ATMDEP_N(JW) + ATM_DEP_LOADING(ATMDCN(JAC, JW), JW)*BI(KT, I)*DLX(I)*3.17098E-11*DLT/1000. ! N MASS BALANCE IN KG - CUMULATIVE
                            end if
                        end if
                    end do
                end do
            end do
        end if


    end do ! JW Loop

!**** Kinetic fluxes

    do JW = 1, NWB
        KT = KTWB(JW) ! SW 10/25/2017
        if (FLUX(JW)) then
            call KINETIC_FLUXES()
        end if
    end do

!SP CEMA
    if (UPDATE_KINETICS) then
        if (CEMARelatedCode .and. IncludeBedConsolidation) then
            call CEMASedimentModel()
        end if
!If(CEMARelatedCode .and. IncludeCEMASedDiagenesis)Call SedimentFlux
        if (CEMARelatedCode .and. IncludeCEMASedDiagenesis .and. BUBBLES_CALCULATION) then
            call CEMACalculateRiseVelocity()
        end if
        if (CEMARelatedCode .and. IncludeCEMASedDiagenesis .and. ApplyBubbTurb .and. BUBBLES_CALCULATION) then
            call CEMABubblesReleaseTurbulence()
        end if
        if (CEMARelatedCode .and. IncludeCEMASedDiagenesis .and. BUBBLES_CALCULATION) then
            call CEMABubblesTransport()
        end if
        if (CEMARelatedCode .and. IncludeCEMASedDiagenesis .and. BUBBLES_CALCULATION) then
            call CEMABubbWatTransfer()
        end if
        if (CEMARelatedCode .and. IncludeCEMASedDiagenesis .and. BUBBLES_CALCULATION) then
            call CEMABubblesRelease()
        end if
        if (IncludeFFTLayer) then
            call CEMAFFTLayerCode()
        end if
    end if
!End SP CEMA

!**** Constituent transport

!!$OMP PARALLEL DO PRIVATE(I,JC,KT,JB,JW,DT,K,IU,ID,BTA1,GMA1)    !I,JC,KT,JW,JB,CNEW,SSB,SSK,COLD,AT,VT,CT,DT) 

    do JAC = 1, NAC !CONCURRENT(JAC=1:NAC)              !JAC=1,NAC
        JC = CN(JAC)
        COLD => C1S(:, :, JC)
        do JW = 1, NWB
            KT = KTWB(JW)
            do JB = BS(JW), BE(JW)
                if (BR_INACTIVE(JB)) then
                    cycle
                end if ! SW 6/12/2017
                IU = CUS(JB)
                ID = DS(JB)
!    DO JAC=1,NAC
!      JC   =  CN(JAC)
!COLD => C1S(:,:,JC)
                call HORIZONTAL_MULTIPLIERS()
                call VERTICAL_MULTIPLIERS()

!       CNEW => C1(:,:,JC)
!       SSB  => CSSB(:,:,JC)
!       SSK  => CSSK(:,:,JC)
!       CALL HORIZONTAL_TRANSPORT
                do I = IU, ID
                    do K = KT, KB(I)
                        DT(K, I) = (C1S(K, I, JC)*BH2(K, I)/DLT + (ADX(K, I)*BHR1(K, I) - ADX(K, I - 1)*BHR1(K, I - 1))/DLX(I) + (1.0D0 - THETA(JW))*(ADZ(K, I)*BB(K, I) - ADZ(K - 1, I)*BB(K - 1, I)) + CSSB(K, I, JC)/DLX(I))*DLT/BH1(K, I) + CSSK(K, I, JC)*DLT
                    end do
                end do
                do I = IU, ID
!      CALL TRIDIAG(AT(:,I),VT(:,I),CT(:,I),DT(:,I),KT,KB(I),KMX,CNEW(:,I))
                    BTA1(KT) = VT(KT, I)
                    GMA1(KT) = DT(KT, I)
                    do K = KT + 1, KB(I)
                        BTA1(K) = VT(K, I) - AT(K, I)/BTA1(K - 1)*CT(K - 1, I)
                        GMA1(K) = DT(K, I) - AT(K, I)/BTA1(K - 1)*GMA1(K - 1)
                    end do
                    C1(KB(I), I, JC) = GMA1(KB(I))/BTA1(KB(I))
                    do K = KB(I) - 1, KT, -1
                        C1(K, I, JC) = (GMA1(K) - CT(K, I)*C1(K + 1, I, JC))/BTA1(K)
                    end do
                end do
            end do
        end do
    end do
!!$OMP END PARALLEL DO
    if (DERIVED_CALC) then
        call DERIVED_CONSTITUENTS()
    end if

end subroutine WQCONSTITUENTS
