subroutine INITCOND()
    use MAIN
    use GLOBAL;     use NAMESC;     use GEOMC;     use LOGICC;     use PREC;     use SURFHE;     use KINETIC;     use SHADEC;     use EDDY
    use STRUCTURES;     use TRANS;     use TVDC;     use SELWC;     use GDAYC;     use SCREENC;     use TDGAS;     use RSTART
    use MACROPHYTEC;     use POROSITYC;     use ZOOPLANKTONC
    use CEMAVars, only: CEMARelatedCode, IncludeCEMASedDiagenesis ! cb 07/23/18
    use CEMASedimentDiagenesis
    implicit none
    external :: RESTART_OUTPUT

    real :: TMAC, XSAR
    character(len=1) :: ICHAR
    character(len=8) :: IBLANK

    BIC = B

    do JW = 1, NWB
        KT = KTWB(JW)
        if (VERT_PROFILE(JW)) then

!**** Temperature and water quality

            open(VPR(JW), FILE=VPRFN(JW), STATUS="OLD")
            read(VPR(JW), "(A1)") ICHAR

            if (ICHAR == "$") then
                read(VPR(JW), "(/)")
                if (VERT_TEMP(JW)) then
                    read(VPR(JW), *) IBLANK, (TVP(K, JW), K = KT, KBMAX(JW))
                end if
                if (CONSTITUENTS) then
                    do JC = 1, NCT
                        if (VERT_CONC(JC, JW)) then
                            read(VPR(JW), *) IBLANK, (CVP(K, JC, JW), K = KT, KBMAX(JW))
                        end if
                    end do
                    do JE = 1, NEP
                        if (VERT_EPIPHYTON(JW, JE)) then
                            read(VPR(JW), *) IBLANK, (EPIVP(K, JW, JE), K = KT, KBMAX(JW))
                        end if
                    end do
                    do m = 1, nmc ! cb 8/21/15
                        if (VERT_macrophyte(JW, m)) then
                            read(VPR(JW), *) IBLANK, (macrcvp(K, JW, m), K = KT, KBMAX(JW))
                        end if
                    end do
                    if (VERT_SEDIMENT(JW)) then
                        read(VPR(JW), *) IBLANK, (SEDVP(K, JW), K = KT, KBMAX(JW))
                    end if
                end if
            else
                if (VERT_TEMP(JW)) then
                    read(VPR(JW), "(//(8X,9F8.0))") (TVP(K, JW), K = KT, KBMAX(JW))
                end if
                if (CONSTITUENTS) then
                    do JC = 1, NCT
                        if (VERT_CONC(JC, JW)) then
                            read(VPR(JW), "(//(8X,9F8.0))") (CVP(K, JC, JW), K = KT, KBMAX(JW))
                        end if
                    end do
                    do JE = 1, NEP
                        if (VERT_EPIPHYTON(JW, JE)) then
                            read(VPR(JW), "(//(8X,9F8.0))") (EPIVP(K, JW, JE), K = KT, KBMAX(JW))
                        end if
                    end do
                    do m = 1, nmc ! cb 8/21/15
                        if (VERT_macrophyte(JW, m)) then
                            read(VPR(JW), "(//(8X,9F8.0))") (macrcvp(K, JW, m), K = KT, KBMAX(JW))
                        end if
                    end do
                    if (VERT_SEDIMENT(JW)) then
                        read(VPR(JW), "(//(8X,9F8.0))") (SEDVP(K, JW), K = KT, KBMAX(JW))
                    end if
                end if
            end if
        end if

!** Longitudinal/vertical initial profiles

        if (LONG_PROFILE(JW)) then
            open(LPR(JW), FILE=LPRFN(JW), STATUS="OLD")
            read(LPR(JW), "(A1)") ICHAR
            if (ICHAR == "$") then
                read(LPR(JW), *)
            end if
        end if

!** Branch related variables

        if (.not. RESTART_IN) then
            if (LONG_TEMP(JW) .and. ICHAR == "$") then
                read(LPR(JW), *)
            end if
            do JB = BS(JW), BE(JW)

!****** Temperature

                do I = CUS(JB), DS(JB)
                    if (LONG_TEMP(JW)) then
                        if (ICHAR == "$") then
                            read(LPR(JW), *) IBLANK, (T1(K, I), K = KT, KB(I))
                        else
                            read(LPR(JW), "(//(8X,9F8.0))") (T1(K, I), K = KT, KB(I))
                        end if
                    end if
                    do K = KT, KB(I)
                        if (ISO_TEMP(JW)) then
                            T1(K, I) = T2I(JW)
                        end if
                        if (VERT_TEMP(JW)) then
                            T1(K, I) = TVP(K, JW)
                        end if
                        T2(K, I) = T1(K, I)
                    end do
                end do
            end do

!**** Constituents

            do JC = 1, NAC
                if (LONG_CONC(CN(JC), JW) .and. ICHAR == "$") then
                    read(LPR(JW), *)
                end if
                do JB = BS(JW), BE(JW)
                    do I = CUS(JB), DS(JB)
                        JAC = CN(JC)
                        if (LONG_CONC(JAC, JW)) then
                            if (ICHAR == "$") then
                                read(LPR(JW), *) IBLANK, (C2(K, I, JAC), K = KT, KB(I))
                            else
                                read(LPR(JW), "(//(8X,9F8.0))") (C2(K, I, JAC), K = KT, KB(I))
                            end if
                        end if
                        do K = KT, KB(I)
                            if (ISO_CONC(JAC, JW)) then
                                C2(K, I, JAC) = C2I(JAC, JW)
                            end if
                            if (VERT_CONC(JAC, JW)) then
                                C2(K, I, JAC) = CVP(K, JAC, JW)
                            end if
                            C1(K, I, JAC) = C2(K, I, JAC)
                            C1S(K, I, JAC) = C1(K, I, JAC)
                        end do
                    end do
                end do
            end do

!**** Epiphyton


            do JE = 1, NEP
                if (EPIPHYTON_CALC(JW, JE)) then
                    if (LONG_EPIPHYTON(JW, JE) .and. ICHAR == "$") then
                        read(LPR(JW), *)
                    end if
                    do JB = BS(JW), BE(JW)
                        do I = CUS(JB), DS(JB)
                            if (LONG_EPIPHYTON(JW, JE)) then
                                if (ICHAR == "$") then
                                    read(LPR(JW), *) IBLANK, (EPD(K, I, JE), K = KT, KB(I))
                                else
                                    read(LPR(JW), "(//(8X,9F8.0))") (EPD(K, I, JE), K = KT, KB(I))
                                end if
                            end if
                            if (ISO_EPIPHYTON(JW, JE)) then
                                EPD(:, I, JE) = EPICI(JW, JE)
                            end if
                            if (VERT_EPIPHYTON(JW, JE)) then
                                EPD(:, I, JE) = EPIVP(:, JW, JE)
                            end if ! CB 5/16/2009
                        end do
                    end do
                end if
            end do

!**** macrophytes - added 8/21/15


            do m = 1, nmc
                if (macrophyte_CALC(JW, m)) then
                    if (LONG_macrophyte(JW, m) .and. ICHAR == "$") then
                        read(LPR(JW), *)
                    end if
                    do JB = BS(JW), BE(JW)
                        do I = CUS(JB), DS(JB)
                            if (LONG_macrophyte(JW, m)) then
                                if (ICHAR == "$") then
                                    read(LPR(JW), *) IBLANK, (macrclp(K, I, m), K = KT, KB(I))
                                else
                                    read(LPR(JW), "(//(8X,9F8.0))") (macrclp(K, I, m), K = KT, KB(I))
                                end if
                            end if
!IF (ISO_macrophyte(JW,m))  macrc(:,I,m) = macwbci(JW,m)
!IF (VERT_macrophyte(JW,m)) macrc(:,I,m) = macrcvp(:,JW,m)    
                        end do
                    end do
                end if
            end do

!**** Sediments

            do JB = BS(JW), BE(JW)
!     SDKV(:,US(JB):DS(JB))=SDK(JW)
                SDKV(:, US(JB) - 1:DS(JB) + 1) = SDK(JW) ! SW 9/28/13
                if (SEDIMENT_CALC(JW)) then
                    if (LONG_SEDIMENT(JW) .and. JB == BS(JW)) then
                        read(LPR(JW), *)
                    end if
                    do I = CUS(JB), DS(JB)
                        if (LONG_SEDIMENT(JW)) then
                            if (ICHAR == "$") then
                                read(LPR(JW), *) IBLANK, (SED(K, I), K = KT, KB(I))
                            else
                                read(LPR(JW), "(//(8X,9F8.0))") (SED(K, I), K = KT, KB(I))
                            end if
                        end if
                        do K = KT, KB(I)
                            if (ISO_SEDIMENT(JW)) then
                                SED(K, I) = SEDCI(JW)
                            end if
                            if (VERT_SEDIMENT(JW)) then
                                SED(K, I) = SEDVP(K, JW)
                            end if
                        end do
                        SED(KT, I) = SED(KT, I)/H2(KT, I)
                        if (CEMARelatedCode .and. IncludeCEMASedDiagenesis) then ! cb 07/23/18
                            SED(KT + 1:KB(I) - 1, I) = SED(KT + 1:KB(I) - 1, I)/H2(KT + 1:KB(I) - 1, I)
                            sed(kb(i), i) = 0.0
                        else
                            SED(KT + 1:KB(I), I) = SED(KT + 1:KB(I), I)/H2(KT + 1:KB(I), I)
                        end if
                    end do
                end if
            end do
            do JB = BS(JW), BE(JW)
                if (SEDIMENT_CALC(JW)) then
                    do I = CUS(JB), DS(JB)
                        do K = KT, KB(I)
                            if (ISO_SEDIMENT(JW)) then
                                SEDP(K, I) = ORGP(JW)*SEDCI(JW)
                            end if
                            if (VERT_SEDIMENT(JW)) then
                                SEDP(K, I) = SEDVP(K, JW)*ORGP(JW)
                            end if
                            if (LONG_SEDIMENT(JW)) then
                                SEDP(K, I) = ORGP(JW)*SED(K, I)
                            end if
                        end do
                        SEDP(KT, I) = SEDP(KT, I)/H2(KT, I)
                        if (CEMARelatedCode .and. IncludeCEMASedDiagenesis) then ! cb 07/23/18
                            SEDp(KT + 1:KB(I) - 1, I) = SEDp(KT + 1:KB(I) - 1, I)/H2(KT + 1:KB(I) - 1, I)
                            sedp(kb(i), i) = 0.0
                        else
                            SEDP(KT + 1:KB(I), I) = SEDP(KT + 1:KB(I), I)/H2(KT + 1:KB(I), I)
                        end if
                    end do
                end if
            end do
            do JB = BS(JW), BE(JW)
                if (SEDIMENT_CALC(JW)) then
                    do I = CUS(JB), DS(JB)
                        do K = KT, KB(I)
                            if (ISO_SEDIMENT(JW)) then
                                SEDn(K, I) = orgn(JW)*sedci(jw)
                            end if
                            if (VERT_SEDIMENT(JW)) then
                                SEDn(K, I) = SEDVP(K, JW)*orgn(jw)
                            end if
                            if (LONG_SEDIMENT(JW)) then
                                sedn(k, i) = orgn(jw)*sed(k, i)
                            end if
                        end do
                        SEDn(KT, I) = SEDn(KT, I)/H2(KT, I)
                        if (CEMARelatedCode .and. IncludeCEMASedDiagenesis) then ! cb 07/23/18                
                            SEDn(KT + 1:KB(I) - 1, I) = SEDn(KT + 1:KB(I) - 1, I)/H2(KT + 1:KB(I) - 1, I)
                            sedn(kb(i), i) = 0.0
                        else
                            SEDn(KT + 1:KB(I), I) = SEDn(KT + 1:KB(I), I)/H2(KT + 1:KB(I), I)
                        end if
                    end do
                end if
            end do
            do JB = BS(JW), BE(JW)
                if (SEDIMENT_CALC(JW)) then
                    do I = CUS(JB), DS(JB)
                        do K = KT, KB(I)
                            if (ISO_SEDIMENT(JW)) then
                                SEDc(K, I) = SEDCI(JW)*orgc(jw)
                            end if
                            if (VERT_SEDIMENT(JW)) then
                                SEDc(K, I) = SEDVP(K, JW)*orgc(jw)
                            end if
                            if (LONG_SEDIMENT(JW)) then
                                sedc(k, i) = orgc(jw)*sed(k, i)
                            end if
                        end do
                        SEDc(KT, I) = SEDc(KT, I)/H2(KT, I)
                        if (CEMARelatedCode .and. IncludeCEMASedDiagenesis) then ! cb 07/23/18                
                            SEDc(KT + 1:KB(I) - 1, I) = SEDc(KT + 1:KB(I) - 1, I)/H2(KT + 1:KB(I) - 1, I)
                            sedc(kb(i), i) = 0.0
                        else
                            SEDc(KT + 1:KB(I), I) = SEDc(KT + 1:KB(I), I)/H2(KT + 1:KB(I), I)
                        end if
                    end do
                end if
            end do

            SED(:, US(BS(JW)):DS(BE(JW))) = SED(:, US(BS(JW)):DS(BE(JW)))*FSED(JW)
            SEDp(:, US(BS(JW)):DS(BE(JW))) = SEDp(:, US(BS(JW)):DS(BE(JW)))*FSED(JW)
            SEDn(:, US(BS(JW)):DS(BE(JW))) = SEDn(:, US(BS(JW)):DS(BE(JW)))*FSED(JW)
            SEDc(:, US(BS(JW)):DS(BE(JW))) = SEDc(:, US(BS(JW)):DS(BE(JW)))*FSED(JW)

!  Amaila start Additional sediment compartments
            do JB = BS(JW), BE(JW)
                if (SEDIMENT_CALC1(JW)) then
!IF(LONG_SEDIMENT(JW).AND.JB==BS(JW))READ (LPR(JW),*)
                    if (LONG_SEDIMENT(JW) .and. ICHAR == "$") then
                        read(LPR(JW), *)
                    end if ! cb 6/10/13
!DO I=CUS(JB),DS(JB)
                    do I = US(JB), DS(JB) ! cb 6/17/17
                        if (LONG_SEDIMENT1(JW)) then
                            if (ICHAR == "$") then
!READ (LPR(JW),*)IBLANK, (SED1(K,I),K=KT,KB(I)) 
                                read(LPR(JW), *) IBLANK, (SED1(K, I), K = 2, KB(I)) ! cb 6/17/17
                            else
!READ (LPR(JW),'(//(8X,9F8.0))') (SED1(K,I),K=KT,KB(I))
                                read(LPR(JW), "(//(8X,9F8.0))") (SED1(K, I), K = 2, KB(I)) ! cb 6/17/17
                            end if
                        end if
                        do K = KT, KB(I)
                            if (ISO_SEDIMENT1(JW)) then
                                SED1(K, I) = SEDCI1(JW)
                            end if
                            if (VERT_SEDIMENT1(JW)) then
                                SED1(K, I) = SEDVP1(K, JW)
                            end if
                        end do
!SED1(KT,I)         = SED1(KT,I)/H2(KT,I)   ! intial conditions for "tree" sediment compartments are given in g/m^3
!SED1(KT+1:KB(I),I) = SED1(KT+1:KB(I),I)/H2(KT+1:KB(I),I)             
                    end do
                end if
            end do

            do JB = BS(JW), BE(JW)
                if (SEDIMENT_CALC2(JW)) then
!IF(LONG_SEDIMENT(JW).AND.JB==BS(JW))READ (LPR(JW),*)
                    if (LONG_SEDIMENT(JW) .and. ICHAR == "$") then
                        read(LPR(JW), *)
                    end if ! cb 6/10/13
!DO I=CUS(JB),DS(JB)
                    do I = US(JB), DS(JB) ! cb 6/17/17
                        if (LONG_SEDIMENT2(JW)) then
                            if (ICHAR == "$") then
!READ (LPR(JW),*)IBLANK, (SED2(K,I),K=KT,KB(I)) 
                                read(LPR(JW), *) IBLANK, (SED2(K, I), K = 2, KB(I)) ! cb 6/17/17
                            else
!READ (LPR(JW),'(//(8X,9F8.0))') (SED2(K,I),K=KT,KB(I))
                                read(LPR(JW), "(//(8X,9F8.0))") (SED2(K, I), K = 2, KB(I)) ! cb 6/17/17
                            end if
                        end if
                        do K = KT, KB(I)
                            if (ISO_SEDIMENT2(JW)) then
                                SED2(K, I) = SEDCI2(JW)
                            end if
                            if (VERT_SEDIMENT2(JW)) then
                                SED2(K, I) = SEDVP2(K, JW)
                            end if
                        end do
!SED2(KT,I)         = SED2(KT,I)/H2(KT,I)      ! intial conditions for "tree" sediment compartments are given in g/m^3
!SED2(KT+1:KB(I),I) = SED2(KT+1:KB(I),I)/H2(KT+1:KB(I),I)             
                    end do
                end if
            end do

            do JB = BS(JW), BE(JW) ! 9/3/17      
                do I = US(JB), DS(JB)
                    do k = kt, kb(i)
                        sdfirstadd(k, i) = .false.
                    end do
                end do
            end do

            SED1(:, US(BS(JW)):DS(BE(JW))) = SED1(:, US(BS(JW)):DS(BE(JW)))*FSEDc1(JW) ! cb 6/7/17
            SED2(:, US(BS(JW)):DS(BE(JW))) = SED2(:, US(BS(JW)):DS(BE(JW)))*FSEDc2(JW)
            sed1ic = sed1 ! cb 6/17/17
            sed2ic = sed2 ! cb 6/17/17


! Amaila end

            do JB = BS(JW), BE(JW)
                do M = 1, NMC
                    if (MACROPHYTE_CALC(JW, M)) then

!C DISTRIBUTING INITIAL MACROPHYTE CONC TO BOTTOM COLUMN CELLS; MACWBCI = G/M^3
                        do I = CUS(JB), DS(JB)

                            DEPKTI = ELWS(I) - EL(KTI(I) + 1, I)

                            if (DEPKTI >= THRKTI) then
                                KTICOL(I) = .true.
                                JT = KTI(I)
                            else
                                KTICOL(I) = .false.
                                JT = KTI(I) + 1
                            end if

                            JE = KB(I)
                            do J = JT, JE
                                if (J <= KT) then
                                    K = KT
                                else
                                    K = J
                                end if
!MACRC(J,K,I,M) = MACWBCI(JW,M)
!SMACRC(J,K,I,M) = MACWBCI(JW,M)                
                                if (ISO_macrophyte(JW, m)) then
                                    macrc(j, k, I, m) = macwbci(JW, m)
                                end if ! cb 8/24/15
                                if (VERT_macrophyte(JW, m)) then
                                    macrc(j, k, I, m) = macrcvp(k, JW, m)
                                end if
                                if (long_macrophyte(JW, m)) then
                                    macrc(j, k, I, m) = macrclp(K, I, m)
                                end if
                                SMACRC(J, K, I, M) = macrc(j, k, I, m)
                            end do
                        end do

                        do I = CUS(JB), DS(JB)
                            TMAC = 0.0
                            XSAR = 0.0
                            do K = KTI(I), KT
                                JT = K
                                JE = KB(I)
                                COLB = EL(K + 1, I)
                                COLDEP = ELWS(I) - COLB
                                do J = JT, JE
                                    TMAC = TMAC + MACRC(J, KT, I, M)*CW(J, I)*COLDEP
                                    XSAR = XSAR + CW(J, I)*COLDEP
                                end do
                            end do
                            MAC(KT, I, M) = TMAC/XSAR
                            SMAC(KT, I, M) = MAC(KT, I, M)

                            do K = KT + 1, KB(I)
                                JT = K
                                JE = KB(I)
                                TMAC = 0.0
                                do J = JT, JE
                                    TMAC = TMAC + MACRC(J, K, I, M)*CW(J, I)
                                end do
                                MAC(K, I, M) = TMAC/B(K, I)
                                SMAC(K, I, M) = MAC(K, I, M)
                            end do
                        end do

                        do I = CUS(JB), DS(JB)
                            JT = KTI(I)
                            JE = KB(I)
                            do J = JT, JE
                                if (J < KT) then
                                    COLB = EL(J + 1, I)
                                else
                                    COLB = EL(KT + 1, I)
                                end if
                                COLDEP = ELWS(I) - COLB
                                MACRM(J, KT, I, M) = MACRC(J, KT, I, M)*COLDEP*CW(J, I)*DLX(I)
                                SMACRM(J, KT, I, M) = MACRM(J, KT, I, M)
                            end do

                            do K = KT + 1, KB(I)

                                JT = K
                                JE = KB(I)

                                do J = JT, JE

                                    MACRM(J, K, I, M) = MACRC(J, K, I, M)*H2(K, I)*CW(J, I)*DLX(I)
                                    SMACRM(J, K, I, M) = MACRM(J, K, I, M)
                                end do

                            end do
                        end do

                    end if
                end do
            end do
! V3.5 END

!**** ENERGY

            do JB = BS(JW), BE(JW)
                if (BR_INACTIVE(JB)) then
                    cycle
                end if ! SW 12/18/2018
                do I = CUS(JB), DS(JB)
                    if (ENERGY_BALANCE(JW)) then
                        do K = KT, KB(I)
                            EBRI(JB) = EBRI(JB) + T2(K, I)*DLX(I)*BH2(K, I)
                        end do
                    end if
                    do K = KT, KB(I)
                        CMBRT(CN(1:NAC), JB) = CMBRT(CN(1:NAC), JB) + C2(K, I, CN(1:NAC))*DLX(I)*BH2(K, I)
                    end do
                end do

! V3.5 START
!C   INITIALIZING MACROPHYTE TEMPORAL MASS BALANCE TERM....
                do M = 1, NMC
                    if (MACROPHYTE_CALC(JW, M)) then
                        do I = CUS(JB), DS(JB)
                            if (KTICOL(I)) then
                                JT = KTI(I)
                            else
                                JT = KTI(I) + 1
                            end if
                            JE = KB(I)
                            do J = JT, JE
                                MACMBRT(JB, M) = MACMBRT(JB, M) + MACRM(J, KT, I, M)
                            end do
                            do K = KT + 1, KB(I)
                                JT = K
                                JE = KB(I)
                                do J = JT, JE
                                    MACMBRT(JB, M) = MACMBRT(JB, M) + MACRM(J, K, I, M)
                                end do
                            end do
                        end do
                    end if
                end do

!****** Ice cover

                if (ICE_CALC(JW)) then
                    ICETH(CUS(JB):DS(JB)) = ICETHI(JW) ! SW 9/29/15 only initialize CUS to DS
                    ICE(US(JB):DS(JB)) = ICETH(US(JB):DS(JB)) > 0.0
                    do I = CUS(JB), DS(JB) ! SW 9/29/15
                        ICEBANK(I) = ICETH(I)*BI(KT, I)*DLX(I) ! Initial volume of ice in m3
                    end do
                end if

!****** Vertical eddy viscosity

                IUT = CUS(JB)
                IDT = DS(JB) - 1
                if (UP_HEAD(JB)) then
                    IUT = IU - 1
                end if
                if (DN_HEAD(JB)) then
                    IDT = ID
                end if
                do I = IUT, IDT
                    do K = KT, KB(I) - 1
                        AZ(K, I) = AZMIN
                        TKE(K, I, 1) = 1.25E-7
                        TKE(K, I, 2) = 1.0E-9
                    end do
                end do
                do JWR = 1, NIW
                    if (WEIR_CALC) then
                        AZ(MAX(KT, KTWR(JWR) - 1):KBWR(JWR), IWR(JWR)) = 0.0
                    end if
                end do
            end do
        end if

!** Horizontal diffusivities

        do JB = BS(JW), BE(JW)
            do I = CUS(JB), DS(JB) - 1
                do K = KT, KBMIN(I)
                    DX(K, I) = ABS(DXI(JW)) ! SW 8/2/2017 FIRST TIME STEP EVEN IF NEGATIVE USE AS ABS OF DX SINCE IT WILL ALWAYS BE LESS THAN 1     
                    if (INTERNAL_WEIR(K, I)) then
                        DX(K, I) = 0.0
                    end if
                end do
            end do
        end do
        if (VERT_PROFILE(JW)) then
            close(VPR(JW))
        end if
        if (LONG_PROFILE(JW)) then
            close(LPR(JW))
        end if
    end do

! Atmospheric pressure
    if (CONSTITUENTS) then
        PALT(:) = (1.0 - ELWS(:)/1000.0/44.3)**5.25 ! SW 2/3/08
        YEAROLD = YEAR
        if (CO2YEARLYPPM == "      ON") then
            if (YEAR < 1980) then
                PCO2 = (0.000041392*REAL(YEAR*YEAR*YEAR) - 0.231409975*REAL(YEAR*YEAR) + 430.804190829*REAL(YEAR) - 266735.857433224)*PALT(DS(BE(1)))*1.0E-6 ! PPM CO2 AND ALTITUDE CORRECTION 
            else
                PCO2 = (0.015903*YEAR**2 - 61.799598*YEAR + 60357.055057)*PALT(DS(BE(1)))*1.0E-6
            end if
        else
            PCO2 = PCO2ATMPPM*PALT(DS(BE(1)))*1.0E-6 ! IN ATM
        end if
    end if

    ELWS_INI(:) = ELWS(:) ! systdg - Add ELWS_INI
    return

end subroutine INITCOND
