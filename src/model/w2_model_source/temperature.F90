subroutine temperature()

    use MAIN
    use GLOBAL;     use NAMESC;     use GEOMC;     use LOGICC;     use PREC;     use SURFHE;     use KINETIC;     use SHADEC;     use EDDY
    use STRUCTURES;     use TRANS;     use TVDC;     use SELWC;     use GDAYC;     use SCREENC;     use TDGAS;     use RSTART
    use MACROPHYTEC;     use POROSITYC;     use ZOOPLANKTONC
    use CEMAVars
    implicit none
    external :: RESTART_OUTPUT

    real(R8) :: BTA1(1000), GMA1(1000) ! places a limit of 1000 vertical layers
    real :: RN1

    do JW = 1, NWB
        if (READ_EXTINCTION(JW)) then
            GAMMA(:, US(BS(JW)):DS(BE(JW))) = EXH2O(JW)
        end if ! SW 1/28/13
        KT = KTWB(JW)
        if (.not. NO_HEAT(JW)) then
            if (.not. READ_RADIATION(JW)) then
                call SHORT_WAVE_RADIATION(JDAY)
            end if
            if (TERM_BY_TERM(JW)) then ! SW 1/25/05
                if (TAIR(JW) >= 5.0) then
!RANLW(JW) = 5.31D-13*(273.15D0+TAIR(JW))**6*(1.0D0+0.0017D0*CLOUD(JW)**2)*0.97D0
                    RANLW(JW) = 5.31D-13*(273.15D0 + TAIR(JW))**6*(1.0D0 + 0.0017D0*CLOUD(JW)*CLOUD(JW))*0.97D0 ! SW 4/20/16 SPEED
                else
!RANLW(JW) = 5.62D-8*(273.15D0+TAIR(JW))**4*(1.D0-0.261D0*DEXP(-7.77D-4*TAIR(JW)**2))*(1.0D0+0.0017D0*CLOUD(JW)**2)*0.97D0
                    RANLW(JW) = 5.62D-8*(273.15D0 + TAIR(JW))**4*(1.D0 - 0.261D0*DEXP((-7.77D-4)*TAIR(JW)*TAIR(JW)))*(1.0D0 + 0.0017D0*CLOUD(JW)*CLOUD(JW))*0.97D0 ! SW 4/20/16 SPEED
                end if
            end if
        end if
        do JB = BS(JW), BE(JW)
            if (BR_INACTIVE(JB)) then
                cycle
            end if
            IU = CUS(JB)
            ID = DS(JB)

!****** Heat exchange

            if (.not. NO_HEAT(JW)) then
                do I = IU, ID
                    if (DYNAMIC_SHADE(I)) then
                        call SHADING()
                    end if

!********** Surface

                    if (.not. ICE(I)) then
                        if (TERM_BY_TERM(JW)) then
                            call SURFACE_TERMS(T2(KT, I))
                            RS(I) = SRON(JW)*SHADE(I)
                            RN(I) = RS(I) + RANLW(JW) - RB(I) - RE(I) - RC(I)
                            HEATEX = RN(I)/RHOWCP*BI(KT, I)*DLX(I)
                        else
                            call EQUILIBRIUM_TEMPERATURE()
                            HEATEX = (ET(I) - T2(KT, I))*CSHE(I)*BI(KT, I)*DLX(I)
                        end if
                        TSS(KT, I) = TSS(KT, I) + HEATEX
                        TSSS(JB) = TSSS(JB) + HEATEX*DLT
                        SROOUT = (1.0D0 - BETA(JW))*SRON(JW)*SHADE(I)/RHOWCP*BI(KT, I)*DLX(I)*DEXP((-GAMMA(KT, I))*DEPTHB(KT, I))
                        TSS(KT, I) = TSS(KT, I) - SROOUT
                        TSSS(JB) = TSSS(JB) - SROOUT*DLT
                        if (KT == KB(I)) then ! SW 4/18/07
                            SROSED = SROOUT*TSEDF(JW)
                        else
                            SROSED = SROOUT*(1.0D0 - BI(KT + 1, I)/BI(KT, I))*TSEDF(JW)
                        end if
                        TSS(KT, I) = TSS(KT, I) + SROSED
                        TSSS(JB) = TSSS(JB) + SROSED*DLT
                        SROIN = SROOUT*B(KT + 1, I)/BI(KT, I)
                        do K = KT + 1, KB(I)
                            SROOUT = SROIN*DEXP((-GAMMA(K, I))*H1(K, I))
                            SRONET = SROIN - SROOUT
                            if (K /= KB(I)) then ! SW 1/18/08
                                SROSED = SROOUT*(1.0D0 - BI(K + 1, I)/BI(K, I))*TSEDF(JW)
                            else
                                SROSED = SROOUT*TSEDF(JW)
                            end if
                            TSS(K, I) = TSS(K, I) + SRONET + SROSED
                            TSSS(JB) = TSSS(JB) + (SRONET + SROSED)*DLT
                            SROIN = SROOUT*B(K + 1, I)/B(K, I)
                        end do
                    end if

!********** Sediment/water

                    do K = KT, KB(I)
                        if (K == KB(I)) then ! SW 4/18/07
                            TFLUX = CBHE(JW)/RHOWCP*(TSED(JW) - T2(K, I))*BI(K, I)*DLX(I)
                        else
                            TFLUX = CBHE(JW)/RHOWCP*(TSED(JW) - T2(K, I))*(BI(K, I) - BI(K + 1, I))*DLX(I)
                        end if
                        TSS(K, I) = TSS(K, I) + TFLUX
                        TSSB(JB) = TSSB(JB) + TFLUX*DLT
                    end do
                end do

!******** Ice cover

                if (ICE_CALC(JW)) then
!           HIA = 0.2367*CSHE(I)/5.65E-8    ! SW 10/20/09 Duplicate line of code
                    do I = IU, ID
                        ALLOW_ICE(I) = .true.
                        if (T2(KT, I) > ICET2(JW)) then
                            ALLOW_ICE(I) = .false.
                        end if ! RC/SW 4/28/11
!             DO K=KT,KB(I)                                                          ! RC 4/28/11: no reason to loop over all layers
!               IF (T2(K,I) > ICET2(JW)) ALLOW_ICE(I) = .FALSE.
!             END DO
                    end do
!          ICE_IN(JB) = .TRUE.                                                      ! RC/SW 4/28/11 eliminate ICE_IN
!          DO I=IU,ID
!            IF (ICETH(I) < ICEMIN(JW)) ICE_IN(JB) = .FALSE.
!          END DO
                    do I = IU, ID
                        if (SALT_WATER(JW)) then ! SW/RC 4/28/11
                            if (TDS(KT, I) < 35.) then
                                RIMT = (-0.0545)*TDS(KT, I) ! REGRESSION FOR TDS BETWEEN 0 AND 35 PPT
                            else
                                RIMT = (-0.31462) - 0.04177*TDS(KT, I) - 0.000166*TDS(KT, I)*TDS(KT, I) ! REGRESSION EQN FOR TDS>35 PPT
                            end if
                        else
                            RIMT = 0.0
                        end if
                        if (DETAILED_ICE(JW)) then
                            if (T2(KT, I) < 0.0) then
                                if (.not. ICE(I)) then
                                    ICETH2 = (-T2(KT, I))*RHO(KT, I)*CP*H2(KT, I)/RHOIRL1
                                    if (ICETH2 < ICE_TOL) then
                                        ICETH2 = 0.0D0
                                    else
                                        TFLUX = T2(KT, I)*RHO(KT, I)*CP*H2(KT, I)*BI(KT, I)/(RHOWCP*DLT)*DLX(I)
                                        TSS(KT, I) = TSS(KT, I) - TFLUX
                                        TSSICE(JB) = TSSICE(JB) - TFLUX*DLT
                                    end if
                                end if
                            end if

!************** Ice balance

                            if (ICE(I)) then
                                TICE = TAIR(JW)
                                DEL = 2.0D0
                                J = 1
                                if (TAIR(JW) >= 5.0) then
                                    RANLW(JW) = 5.31D-13*(273.15D0 + TAIR(JW))**6*(1.0D0 + 0.0017D0*CLOUD(JW)**2)*0.97D0
                                else
                                    RANLW(JW) = 5.62D-8*(273.15D0 + TAIR(JW))**4*(1.D0 - 0.261D0*DEXP((-7.77D-4)*TAIR(JW)**2))*(1.0D0 + 0.0017D0*CLOUD(JW)**2)*0.97D0
                                end if
                                RN1 = SRON(JW)/REFL*SHADE(I)*(1.0D0 - ALBEDO(JW))*BETAI(JW) + RANLW(JW) ! SW 4/19/10 eliminate spurious divsion of SRO by RHOCP
                                do while (ABS(DEL) > 1.0 .and. J < 500) ! SW 4/21/10 Should have been ABS of DEL
                                    call SURFACE_TERMS(TICE)
                                    RN(I) = RN1 - RB(I) - RE(I) - RC(I) ! 4/19/10 
!                    RN(I) = SRON(JW)/(REFL*RHOWCP)*SHADE(I)*(1.0-ALBEDO(JW))*BETAI(JW)+RANLW(JW)-RB(I)-RE(JW)-RC(I)
                                    DEL = RN(I) + RK1*(RIMT - TICE)/ICETH(I) ! RK1 is ice conductivity 2.12 W/m/oC
                                    if (ABS(DEL) > 1.0) then
                                        TICE = TICE + DEL/500.0D0
                                    end if
                                    J = J + 1
                                end do

!**************** Solar radiation attenuation

                                TFLUX = DLX(I)*SRON(JW)/(RHOWCP*REFL)*SHADE(I)*(1.0D0 - ALBEDO(JW))*(1.0D0 - BETAI(JW))*DEXP((-GAMMAI(JW))*ICETH(I))*BI(KT, I)
                                TSS(KT, I) = TSS(KT, I) + TFLUX
                                TSSICE(JB) = TSSICE(JB) + TFLUX*DLT
                                if (TICE > 0.0) then
                                    HICE = RHOICP*0.5D0*TICE*0.5D0*ICETH(I)*BI(KT, I)/(RHOWCP*DLT)
                                    ICETHU = (-DLT)*HICE/B(KTI(I), I)*RHOWCP/RHOIRL1
                                    TICE = 0.0D0
                                end if

!**************** Ice growth

                                if (TICE < 0.0) then
                                    ICETH1 = DLT*RK1*(RIMT - TICE)/ICETH(I)/RHOIRL1
                                end if

!**************** Ice melt from water-ice interface

                                if (T2(KT, I) > 0.0) then
                                    ICETH2 = (-DLT)*HWI(JW)*(T2(KT, I) - RIMT)/RHOIRL1
                                    TFLUX = 2.392D-7*HWI(JW)*(RIMT - T2(KT, I))*BI(KT, I)*DLX(I)
                                    TSS(KT, I) = TSS(KT, I) + TFLUX
                                    TSSICE(JB) = TSSICE(JB) + TFLUX*DLT
                                end if
                            end if

!************** Ice thickness

                            ICETH(I) = ICETH(I) + ICETHU + ICETH1 + ICETH2
! SW 9/4/15
                            if (ICEC(JW) == "    ONWB") then
                                IceThicknessChange = ICETHU + ICETH1 + ICETH2
                                if (IceThicknessChange >= 0.0 .and. ICETH(I) > ICE_TOL) then
                                    IceQSS(I) = (-IceThicknessChange)*BI(KT, I)*DLX(I)*0.917/DLT ! CEMA: removal of 0.917 m3 of water for every 1 m3 of ice formed
                                    IceBank(I) = IceBank(I) + IceThicknessChange*BI(KT, I)*DLX(I) ! Icebank - always + - volume of ice in m3
                                else
                                    if (ICETH(I) < ICE_TOL) then
                                        ICETH(I) = 0.0D0
                                        ICEQSS(I) = ICEBANK(I)*0.917/DLT
                                        ICEBANK(I) = 0.0
                                        if (I == IU .and. US(JB) < IU) then ! CHECK SUBTRACTED SEGMENTS THAT MAY HAVE ICE, MELT THEM AND PUT WATER IN SEGMENT IU
                                            do II = US(JB), IU - 1
                                                ICETH(II) = 0.0D0
                                                ICEQSS(IU) = ICEQSS(IU) + ICEBANK(II)*0.917/DLT
                                                ICEBANK(II) = 0.0
                                                ICE(II) = .false.
                                            end do
                                        end if

                                    else
                                        ICEQSS(I) = (-ICETHICKNESSCHANGE)/ICETH(I)*ICEBANK(I)*0.917/DLT ! SW 9/29/15
                                        ICEBANK(I) = (1. + ICETHICKNESSCHANGE/ICETH(I))*ICEBANK(I) ! Note: ICETHICKNESSCHANGE is negative
                                        if (I == IU .and. US(JB) < IU) then ! CHECK SUBTRACTED SEGMENTS THAT MAY HAVE ICE, MELT THEM AND PUT WATER IN SEGMENT IU
                                            do II = US(JB), IU - 1
                                                ICEQSS(IU) = ICEQSS(IU) - ICETHICKNESSCHANGE/ICETH(IU)*ICEBANK(II)*0.917/DLT
                                                ICEBANK(II) = (1. + ICETHICKNESSCHANGE/ICETH(IU))*ICEBANK(II)
                                            end do
                                        end if
                                    end if
                                end if
!VolIce(jb)=VolIce(jb)+iceqss(i)*dlt   ! since this flow is not exercised until the next time step - moved code to main program calculation of qss
                            end if

                            if (ICETH(I) < ICE_TOL) then
                                ICETH(I) = 0.0D0
                            end if
!           IF (WINTER .AND. (.NOT. ICE_IN(JB))) THEN            ! RC 4/28/11 No reason for this
!             IF (.NOT. ALLOW_ICE(I)) ICETH(I) = 0.0
!           END IF
                            ICE(I) = ICETH(I) > 0.0
                            if (ICE(I)) then ! 3/27/08 SW
                                ICESW(I) = 0.0
                            else
                                ICESW(I) = 1.0
                            end if
                            ICETHU = 0.0
                            ICETH1 = 0.0
                            ICETH2 = 0.0
                            if (ICETH(I) < ICE_TOL .and. ICETH(I) > 0.0) then
                                ICETH(I) = ICE_TOL
                            end if
                        else ! IF no ice the preceding time step
                            if (TERM_BY_TERM(JW)) then
                                call EQUILIBRIUM_TEMPERATURE()
                            end if ! SW 10/20/09 Must call this first otherwise ET and CSHE are 0
                            HIA = 0.2367D0*CSHE(I)/5.65D-8 ! JM 11/08 convert SI units of m/s to English (btu/ft2/d/F) and then back to SI W/m2/C
!                ICETH(I) = MAX(0.0,ICETH(I)+DLT*((RIMT-ET(I))/(ICETH(I)/RK1+1.0/HIA)-(T2(KT,I)-RIMT))/RHOIRL1)
                            ICETH(I) = MAX(0.0, ICETH(I) + DLT*((RIMT - ET(I))/(ICETH(I)/RK1 + 1.0D0/HIA) - HWI(JW)*(T2(KT, I) - RIMT))/RHOIRL1) ! SW 10/20/09 Revised missing HWI(JW)
                            ICE(I) = ICETH(I) > 0.0
                            ICESW(I) = 1.0
                            if (ICE(I)) then
!                  TFLUX      = 2.392E-7*(RIMT-T2(KT,I))*BI(KT,I)*DLX(I)
                                TFLUX = 2.392D-7*HWI(JW)*(RIMT - T2(KT, I))*BI(KT, I)*DLX(I) ! SW 10/20/09 Revised missing HWI(JW)
                                TSS(KT, I) = TSS(KT, I) + TFLUX
                                TSSICE(JB) = TSSICE(JB) + TFLUX*DLT
                                ICESW(I) = 0.0
                            end if
                        end if
                    end do
                end if
            end if
!DO I=IU,ID
!  icebank_all=icebank(i)+icebank_all
!END DO
!write(25600,*)jday,icebank_all
!****** Heat sources/sinks and total inflow/outflow

            if (EVAPORATION(JW)) then
                do I = IU, ID
                    TSS(KT, I) = TSS(KT, I) - EV(I)*T2(KT, I)
                    TSSEV(JB) = TSSEV(JB) - EV(I)*T2(KT, I)*DLT
                    VOLEV(JB) = VOLEV(JB) - EV(I)*DLT
                end do
            end if
            if (PRECIPITATION(JW)) then
                do I = IU, ID
                    TSS(KT, I) = TSS(KT, I) + QPR(I)*TPR(JB)
                    TSSPR(JB) = TSSPR(JB) + QPR(I)*TPR(JB)*DLT
                    VOLPR(JB) = VOLPR(JB) + QPR(I)*DLT
                end do
            end if
            if (TRIBUTARIES) then
                do JT = 1, JTT
                    if (JB == JBTR(JT)) then
                        I = ITR(JT)
                        if (I < CUS(JB)) then
                            I = CUS(JB)
                        end if
                        do K = KTTR(JT), KBTR(JT)
                            if (QTR(JT) < 0) then
                                TSS(K, I) = TSS(K, I) + T2(K, I)*QTR(JT)*QTRF(K, JT)
                                TSSTR(JB) = TSSTR(JB) + T2(K, I)*QTR(JT)*QTRF(K, JT)*DLT
                            else
                                TSS(K, I) = TSS(K, I) + TTR(JT)*QTR(JT)*QTRF(K, JT)
                                TSSTR(JB) = TSSTR(JB) + TTR(JT)*QTR(JT)*QTRF(K, JT)*DLT
                            end if
                        end do
                        VOLTRB(JB) = VOLTRB(JB) + QTR(JT)*DLT
                    end if
                end do
            end if
            if (DIST_TRIBS(JB)) then
                do I = IU, ID
                    if (QDT(I) < 0) then
                        TSS(KT, I) = TSS(KT, I) + T2(KT, I)*QDT(I)
                        TSSDT(JB) = TSSDT(JB) + T2(KT, I)*QDT(I)*DLT
                    else
                        TSS(KT, I) = TSS(KT, I) + TDTR(JB)*QDT(I)
                        TSSDT(JB) = TSSDT(JB) + TDTR(JB)*QDT(I)*DLT
                    end if
                    VOLDT(JB) = VOLDT(JB) + QDT(I)*DLT
                end do
            end if
            if (WITHDRAWALS) then
                do JWD = 1, JWW
                    if (QWD(JWD) /= 0.0) then
                        if (JB == JBWD(JWD)) then
                            I = MAX(CUS(JBWD(JWD)), IWD(JWD))
                            do K = KTW(JWD), KBW(JWD)
                                TSS(K, I) = TSS(K, I) - T2(K, I)*QSW(K, JWD)
                                TSSWD(JB) = TSSWD(JB) - T2(K, I)*QSW(K, JWD)*DLT
                            end do
                            VOLWD(JB) = VOLWD(JB) - QWD(JWD)*DLT
                        end if
                    end if
                end do
            end if
            if (UP_FLOW(JB)) then
                do K = KT, KB(IU)
                    if (.not. HEAD_FLOW(JB)) then
                        TSS(K, IU) = TSS(K, IU) + QINF(K, JB)*QIN(JB)*TIN(JB)
                        TSSIN(JB) = TSSIN(JB) + QINF(K, JB)*QIN(JB)*TIN(JB)*DLT
                    else
                        if (U(K, IU - 1) >= 0.0) then
                            TSS(K, IU) = TSS(K, IU) + U(K, IU - 1)*BHR1(K, IU - 1)*T1(K, IU - 1)
                            TSSIN(JB) = TSSIN(JB) + U(K, IU - 1)*BHR1(K, IU - 1)*T1(K, IU - 1)*DLT
                        else
                            TSS(K, IU) = TSS(K, IU) + U(K, IU - 1)*BHR1(K, IU - 1)*T1(K, IU)
                            TSSIN(JB) = TSSIN(JB) + U(K, IU - 1)*BHR1(K, IU - 1)*T1(K, IU)*DLT
                        end if
                    end if
                end do
                VOLIN(JB) = VOLIN(JB) + QIN(JB)*DLT
            end if
            if (DN_FLOW(JB)) then
                do K = KT, KB(ID)
                    TSS(K, ID) = TSS(K, ID) - QOUT(K, JB)*T2(K, ID + 1)
                    TSSOUT(JB) = TSSOUT(JB) - QOUT(K, JB)*T2(K, ID + 1)*DLT
                    VOLOUT(JB) = VOLOUT(JB) - QOUT(K, JB)*DLT
                end do
            end if
            if (UP_HEAD(JB)) then
                do K = KT, KB(IU)
                    IUT = IU
                    if (QUH1(K, JB) >= 0.0) then
                        IUT = IU - 1
                    end if
                    TSSUH1(K, JB) = T2(K, IUT)*QUH1(K, JB)
                    TSS(K, IU) = TSS(K, IU) + TSSUH1(K, JB)
                    TSSUH(JB) = TSSUH(JB) + TSSUH1(K, JB)*DLT
                    VOLUH(JB) = VOLUH(JB) + QUH1(K, JB)*DLT
                end do
            end if
            if (UH_INTERNAL(JB)) then
                if (UHS(JB) /= DS(JBUH(JB)) .or. DHS(JBUH(JB)) /= US(JB)) then
                    if (JBUH(JB) >= BS(JW) .and. JBUH(JB) <= BE(JW)) then
                        do K = KT, KB(IU - 1)
                            TSS(K, UHS(JB)) = TSS(K, UHS(JB)) - TSSUH2(K, JB)/DLT
                            TSSUH(JBUH(JB)) = TSSUH(JBUH(JB)) - TSSUH2(K, JB)
                            VOLUH(JBUH(JB)) = VOLUH(JBUH(JB)) - VOLUH2(K, JB)
                        end do
                    else
                        call UPSTREAM_CONSTITUENT(T2, TSS)
                        do K = KT, KB(IU - 1)
                            TSSUH(JBUH(JB)) = TSSUH(JBUH(JB)) - TSSUH2(K, JB)
                            VOLUH(JBUH(JB)) = VOLUH(JBUH(JB)) - VOLUH2(K, JB)
                        end do
                    end if
                end if
            end if
            if (DN_HEAD(JB)) then
                do K = KT, KB(ID + 1)
                    IDT = ID + 1
                    if (QDH1(K, JB) >= 0.0) then
                        IDT = ID
                    end if
                    TSSDH1(K, JB) = T2(K, IDT)*QDH1(K, JB)
                    TSS(K, ID) = TSS(K, ID) - TSSDH1(K, JB)
                    TSSDH(JB) = TSSDH(JB) - TSSDH1(K, JB)*DLT
                    VOLDH(JB) = VOLDH(JB) - QDH1(K, JB)*DLT
                end do
            end if
            if (DH_INTERNAL(JB)) then
                if (DHS(JB) /= US(JBDH(JB)) .or. UHS(JBDH(JB)) /= DS(JB)) then
                    if (JBDH(JB) >= BS(JW) .and. JBDH(JB) <= BE(JW)) then
                        do K = KT, KB(ID + 1)
                            TSS(K, CDHS(JB)) = TSS(K, CDHS(JB)) + TSSDH2(K, JB)/DLT
                            TSSDH(JBDH(JB)) = TSSDH(JBDH(JB)) + TSSDH2(K, JB)
                            VOLDH(JBDH(JB)) = VOLDH(JBDH(JB)) + VOLDH2(K, JB)
                        end do
                    else
                        call DOWNSTREAM_CONSTITUENT(T2, TSS)
                        do K = KT, KB(ID + 1)
                            TSSDH(JBDH(JB)) = TSSDH(JBDH(JB)) + TSSDH2(K, JB)
                            VOLDH(JBDH(JB)) = VOLDH(JBDH(JB)) + VOLDH2(K, JB)
                        end do
                    end if
                end if
            end if
        end do
    end do

    if (NIT == 0) then
        if (DERIVED_CALC) then
            do JW = 1, NWB
                KT = KTWB(JW)
                do JB = BS(JW), BE(JW)
                    IU = CUS(JB)
                    ID = DS(JB)
                    call TEMPERATURE_RATES()
                    call KINETIC_RATES()
                    if (CDWBC(PH_DER, JW) == "      ON") then
                        call PH_CO2()
                    end if
                end do
                call DERIVED_CONSTITUENTS()
            end do
        end if
        call OUTPUTA() ! INITIALIZE OUTPUT FOR FIRST TIME STEP IF OUTPUT IS SET FOR FIRST TIME STEP ONLY ISSUE IS VELOCITY AND FLOW IS AT T+DT RATHER THAN AT T=0
    end if

!** Temperature transport
    COLD => HYD(:, :, 4)
    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            if (BR_INACTIVE(JB)) then
                cycle
            end if
            IU = CUS(JB)
            ID = DS(JB)
!COLD => HYD(:,:,4)
            call HORIZONTAL_MULTIPLIERS1()
            call VERTICAL_MULTIPLIERS1()
            call HORIZONTAL_MULTIPLIERS()
            call VERTICAL_MULTIPLIERS()
!   CNEW => T1(:,:)
!   SSB  => TSS(:,:)
!   SSK  => CSSB(:,:,1)
!   CALL HORIZONTAL_TRANSPORT

            do I = IU, ID !CONCURRENT(I=IU:ID)   !
                do K = KT, KB(I) !CONCURRENT(K=KT:KB(I))      !FORALL                                         !DO K=KT,KB(I)
                    DT(K, I) = (COLD(K, I)*BH2(K, I)/DLT + (ADX(K, I)*BHR1(K, I) - ADX(K, I - 1)*BHR1(K, I - 1))/DLX(I) + (1.0D0 - THETA(JW))*(ADZ(K, I)*BB(K, I) - ADZ(K - 1, I)*BB(K - 1, I)) + TSS(K, I)/DLX(I))*DLT/BH1(K, I)
                end do
            end do


            do I = IU, ID !CONCURRENT(I=IU:ID)   !
                do K = KT, KB(I)
                    AT(K, I) = 0.0D0;                     CT(K, I) = 0.0D0;                     VT(K, I) = 0.0D0 !; DT(:,I) = 0.0D0    SW CODE SPEEDUP 6/15/13
                end do
                do K = KT, KB(I) !CONCURRENT(K=KT:KB(I))          !FORALL(K=KT:KB(I))                                                 !DO K=KT,KB(I)
                    AT(K, I) = (-DLT)/BH1(K, I)*BB(K - 1, I)*(DZ(K - 1, I)/AVH1(K - 1, I) + THETA(JW)*0.5D0*W(K - 1, I))
                    CT(K, I) = DLT/BH1(K, I)*BB(K, I)*(THETA(JW)*0.5D0*W(K, I) - DZ(K, I)/AVH1(K, I))
                    VT(K, I) = 1.0D0 + DLT/BH1(K, I)*(BB(K, I)*(DZ(K, I)/AVH1(K, I) + THETA(JW)*0.5D0*W(K, I)) + BB(K - 1, I)*(DZ(K - 1, I)/AVH1(K - 1, I) - THETA(JW)*0.5D0*W(K - 1, I)))
!     DT(K,I) =  CNEW(K,I)
                end do
! CALL TRIDIAG(AT(:,I),VT(:,I),CT(:,I),DT(:,I),KT,KB(I),KMX,CNEW(:,I))
                BTA1(KT) = VT(KT, I)
                GMA1(KT) = DT(KT, I)
                do K = KT + 1, KB(I)
                    BTA1(K) = VT(K, I) - AT(K, I)/BTA1(K - 1)*CT(K - 1, I)
                    GMA1(K) = DT(K, I) - AT(K, I)/BTA1(K - 1)*GMA1(K - 1)
                end do
                T1(KB(I), I) = GMA1(KB(I))/BTA1(KB(I))
                do K = KB(I) - 1, KT, -1
                    T1(K, I) = (GMA1(K) - CT(K, I)*T1(K + 1, I))/BTA1(K)
                end do
            end do
        end do
    end do


end subroutine temperature
