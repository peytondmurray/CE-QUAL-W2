subroutine LAYERADDSUB()
    use MAIN
    use GLOBAL;     use NAMESC;     use GEOMC;     use LOGICC;     use PREC;     use SURFHE;     use KINETIC;     use SHADEC;     use EDDY
    use STRUCTURES;     use TRANS;     use TVDC;     use SELWC;     use GDAYC;     use SCREENC;     use TDGAS;     use RSTART
    use MACROPHYTEC;     use POROSITYC;     use ZOOPLANKTONC
    implicit none

    integer :: KTMAX, JBIZ, KKB, I_BR_NUM
    real(R8) :: TMAC
    real(R8) :: W1, W2, W3, DUMMY

!***********************************************************************************************************************************
!**                                       Task 2.5: Layer - Segment Additions and Subtractions                                    **
!***********************************************************************************************************************************

!** Water surface minimum thickness

    do JW = 1, NWB
        KT = KTWB(JW)
        ZMIN(JW) = -1000.0
        KTMAX = 2 ! SR 10/17/05
        do JB = BS(JW), BE(JW)
            if (BR_INACTIVE(JB)) then
                if (DS(JB) - US(JB) + 1 >= 3) then
                    I_BR_NUM = DS(JB) - 2
                else
                    I_BR_NUM = DS(JB) - 1 ! FOR BRANCHES WITH LESS THAN 3 SEGMENTS
                end if

                if (CUS(JBDH(JB)) <= DHS(JB) .and. ELWS(DHS(JB)) > EL(KB(I_BR_NUM), I_BR_NUM)) then ! ***
                    BR_INACTIVE(JB) = .false.
                    if (SNAPSHOT(JW)) then
                        write(SNP(JW), '(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,13("*"))') "   Branch Active: ", jb, " at Julian day = ", JDAY, "   NIT = ", NIT
                    end if
                    write(WRN, '(1X,13("*"),1X,A,I0,A,F0.3,A,I0)') "   Branch Active: ", jb, " at Julian day = ", JDAY, "   NIT = ", NIT ! SW 7/24/2018
                    if (SNAPSHOT(JW)) then
                        write(SNP(JW), "(/17X,2(A,I0))") " Add segments ", DS(JB) - 1, " through ", DS(JB)
                    end if
                    write(WRN, "(/17X,2(A,I0))") " Add segments ", DS(JB) - 1, " through ", DS(JB)
                    CUS(JB) = DS(JB) - 1
                    do I = DS(JB) - 1, DS(JB)
                        Z(I) = Z(DHS(JB))
                        KTI(I) = KTI(DHS(JB))
                        H1(KT + 1, I) = H(KT + 1, JW)
                        AVH1(KT + 1, I) = (H1(KT + 1, I) + H1(KT + 2, I))*0.5
                        AVHR(KT + 1, I) = H1(KT + 1, I) + (H1(KT + 1, I + 1) - H1(KT + 1, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04
                        if (.not. TRAPEZOIDAL(JW)) then
                            BH1(KT + 1, I) = B(KT + 1, I)*H(KT + 1, JW)
                            H1(KT, I) = H(KT, JW) - Z(I)
                            BI(KT:KB(I), I) = B(KT:KB(I), I) ! SW 4/18/07
                            BI(KT, I) = B(KTI(I), I)
                            AVH1(KT, I) = (H1(KT, I) + H1(KT + 1, I))*0.5
                            AVHR(KT, I) = H1(KT, I) + (H1(KT, I + 1) - H1(KT, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04
                            BH1(KT, I) = BI(KT, I)*(EL(KT, I) - Z(I)*COSA(JB) - EL(KTI(I) + 1, I))/COSA(JB)
                            if (KTI(I) >= KB(I)) then
                                BH1(KT, I) = B(KT, I)*H1(KT, I)
                            end if
                            do K = KTI(I) + 1, KT
                                BH1(KT, I) = BH1(KT, I) + BH(K, I)
                            end do
                        else
                            call GRID_AREA1(EL(KT, I) - Z(I), EL(KT + 1, I), BH1(KT, I), BI(KT, I)) !SW 08/03/04
                            BH1(KT + 1, I) = 0.25*H(KT + 1, JW)*(BB(KT, I) + 2.*B(KT + 1, I) + BB(KT + 1, I))
                            H1(KT, I) = H(KT, JW) - Z(I)
                            AVH1(KT, I) = (H1(KT, I) + H1(KT + 1, I))*0.5
                            AVHR(KT, I) = H1(KT, I) + (H1(KT, I + 1) - H1(KT, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04
                        end if
                        BKT(I) = BH1(KT, I)/H1(KT, I)
                        DEPTHB(KT, I) = H1(KT, I)
                        DEPTHM(KT, I) = H1(KT, I)*0.5
                        do K = KT + 1, KB(I)
                            DEPTHB(K, I) = DEPTHB(K - 1, I) + H1(K, I)
                            DEPTHM(K, I) = DEPTHM(K - 1, I) + (H1(K - 1, I) + H1(K, I))*0.5
                        end do
                    end do
                    do I = DS(JB) - 1, DS(JB)
                        BHR1(KT + 1, I) = BH1(KT + 1, I) + (BH1(KT + 1, I + 1) - BH1(KT + 1, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04
                        BHR1(KT, I) = BH1(KT, I) + (BH1(KT, I + 1) - BH1(KT, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04
                        if (CONSTRICTION(KT, I)) then ! SW 6/26/2018
                            if (BHR1(KT, I) > BCONSTRICTION(I)*H1(KT, I)) then
                                BHR1(KT, I) = BCONSTRICTION(I)*H1(KT, I)
                            end if
                            if (BHR1(KT + 1, I) > BCONSTRICTION(I)*H1(KT + 1, I)) then
                                BHR1(KT + 1, I) = BCONSTRICTION(I)*H1(KT + 1, I)
                            end if
                        end if
                    end do
                    do I = DS(JB) - 1, DS(JB)
                        WIND2(I) = WIND2(DHS(JB))
                        if (DYNAMIC_SHADE(I)) then
                            call SHADING()
                        end if
                        do K = KT, KB(I)
                            U(K, I) = 0.0
                            SDKV(K, I) = SDK(JW)
                            if (DXI(JW) >= 0.0) then
                                DX(K, I) = DXI(JW)
                            else
                                DX(K, I) = ABS(U(K, I))*ABS(DXI(JW))*H(K, JW) ! SW 8/2/2017
                            end if
                            if (INTERNAL_WEIR(K, I)) then
                                DX(K, I) = 0.0
                                U(K, I) = 0.0
                            end if
                            T1(K, I) = T1(K, DHS(JB))
                            T2(K, I) = T1(K, DHS(JB))
                            SU(K, I) = 0.0
                            C1(K, I, CN(1:NAC)) = C1(K, DHS(JB), CN(1:NAC))
                            C2(K, I, CN(1:NAC)) = C1(K, DHS(JB), CN(1:NAC))
                            do JE = 1, NEP
                                EPD(K, I, JE) = 0.01
                                EPC(K, I, JE) = 0.01/H1(K, I)
                            end do
                            CMBRT(CN(1:NAC), JB) = CMBRT(CN(1:NAC), JB) + C1(K, DHS(JB), CN(1:NAC))*DLX(I)*BH1(K, I)
                            EBRI(JB) = EBRI(JB) + T1(K, DHS(JB))*DLX(I)*BH1(K, I)
                        end do
                        do K = KT, KB(I) - 1
                            AZ(K, I) = AZ(K, DHS(JB))
                            TKE(K, I, 1) = TKE(K, DHS(JB), 1) !sg 10/4/07
                            TKE(K, I, 2) = TKE(K, DHS(JB), 2) !sg 10/4/07
                            SAZ(K, I) = AZ(K, DHS(JB))
                            if (INTERNAL_WEIR(K, I)) then
                                AZ(K, I) = 0.0
                                TKE(K, I, 1) = 0.0 !sg 10/4/07
                                TKE(K, I, 2) = 0.0 !sg 10/4/07
                                SAZ(K, I) = 0.0
                            end if
                        end do
                    end do
                end if
            end if
            if (.not. BR_INACTIVE(JB)) then
                do I = CUS(JB), DS(JB)
                    if (KB(I) > KTMAX) then
                        KTMAX = KB(I)
                    end if ! SR 10/17/05
                    if (Z(I) > ZMIN(JW)) then
                        IZMIN(JW) = I
                        JBIZ = JB
                    end if
                    ZMIN(JW) = MAX(ZMIN(JW), Z(I))
                end do
            end if
        end do
        ADD_LAYER = ZMIN(JW) < (-0.85)*H(KT - 1, JW) .and. KT /= 2
        SUB_LAYER = ZMIN(JW) > 0.60*H(KT, JW) .and. KT < KTMAX ! SR 10/17/05
        if (KTWB(JW) == KMX - 1 .and. SLOPE(JBIZ) > 0.0 .and. SUB_LAYER .and. ONE_LAYER(IZMIN(JW))) then
            if (ZMIN(JW) > 0.99*H(KT, JW)) then
                write(WRN, "(A,I0,2(A,F0.3))") "Low water in segment ", IZMIN(JW), " water surface deviation" // " = ", ZMIN(JW), " at day ", JDAY
                WARNING_OPEN = .true.
            end if
            SUB_LAYER = .false.
        end if

        if (ADD_LAYER == .true. .or. SUB_LAYER == .true.) then
            LAYERCHANGE(JW) = .true.
        else
            LAYERCHANGE(JW) = .false.
        end if

!**** Add layers

        do while (ADD_LAYER)
            if (SNAPSHOT(JW)) then
                write(SNP(JW), '(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,13("*"))') "   Add layer ", KT - 1, " at Julian day = ", JDAY, "   NIT = ", NIT, " IZMIN =", IZMIN(JW)
            end if ! SW 1/23/06
            WARNING_OPEN = .true.
            write(WRN, '(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,13("*"))') "   Add layer ", KT - 1, " at Julian day = ", JDAY, "   NIT = ", NIT, " IZMIN =", IZMIN(JW) ! SW 1/23/06

!****** Variable initialization

            KTWB(JW) = KTWB(JW) - 1
            KT = KTWB(JW)
            ilayer = 0

! RECOMPUTE INTERNAL WEIR FOR FLOATING WEIR      
            if (WEIR_CALC) then !  SW 3/16/18
                do JWR = 1, NIW
                    if (IWR(JWR) >= US(BS(JW)) .and. IWR(JWR) <= DS(BE(JW))) then
                        if (EKTWR(JWR) == 0.0) then
                            KTWR(JWR) = KTWB(JW)
                        else
                            KTWR(JWR) = INT(EKTWR(JWR))
                        end if
                        if (EKBWR(JWR) <= 0.0) then
                            do K = KTWR(JWR), KB(IWR(JWR))
                                if (DEPTHB(K, IWR(JWR)) > ABS(EKBWR(JWR))) then
                                    KBWR(JWR) = K
                                    exit
                                end if
                            end do
                        else
                            KBWR(JWR) = INT(EKBWR(JWR))
                        end if

                        do K = 2, KMX - 1
                            if (K >= KTWR(JWR) .and. K <= KBWR(JWR)) then
                                INTERNAL_WEIR(K, IWR(JWR)) = .true.
                            end if
                        end do
                    end if
                end do
            end if


            do JB = BS(JW), BE(JW)
                if (BR_INACTIVE(JB)) then
                    cycle
                end if ! SW 6/12/2017
                IU = CUS(JB)
                ID = DS(JB)
                do I = IU - 1, ID + 1
                    Z(I) = H(KT, JW) + Z(I)
                    H1(KT, I) = H(KT, JW) - Z(I)
                    H1(KT + 1, I) = H(KT + 1, JW)
                    H2(KT + 1, I) = H(KT + 1, JW)
                    AVH1(KT, I) = (H1(KT, I) + H1(KT + 1, I))*0.5D0
                    AVH1(KT + 1, I) = (H1(KT + 1, I) + H1(KT + 2, I))*0.5D0
                    if (.not. TRAPEZOIDAL(JW)) then
                        BH1(KT, I) = BH1(KT + 1, I) - Bnew(KT + 1, I)*H1(KT + 1, I) ! SW 1/23/06
                        BH1(KT + 1, I) = BNEW(KT + 1, I)*H1(KT + 1, I) ! SW 1/23/06
                    else
                        call GRID_AREA1(EL(KT, I) - Z(I), EL(KT + 1, I), BH1(KT, I), DUMMY) !SW 08/03/04
                        BH1(KT + 1, I) = 0.25D0*H1(KT + 1, JW)*(BB(KT, I) + 2.D0*B(KT + 1, I) + BB(KT + 1, I))
                    end if
                    VOL(KT, I) = BH1(KT, I)*DLX(I)
                    VOL(KT + 1, I) = BH1(KT + 1, I)*DLX(I)
                    BKT(I) = BH1(KT, I)/H1(KT, I)
                    DEPTHB(KT, I) = H1(KT, I)
                    DEPTHM(KT, I) = H1(KT, I)*0.5
                    BI(KT:KB(I), I) = B(KT:KB(I), I) ! SW 8/26/05
                    BI(KT, I) = B(KTI(I), I)
                    T1(KT, I) = T1(KT + 1, I)
                    SDKV(KT, I) = SDKV(KT + 1, I) ! SW 1/18/08
                    SED(KT, I) = SED(KT + 1, I) ! SW 1/18/08
                    SEDN(KT, I) = SEDN(KT + 1, I) ! SW 1/18/08
                    SEDP(KT, I) = SEDP(KT + 1, I) ! SW 1/18/08
                    SEDC(KT, I) = SEDC(KT + 1, I) ! SW 1/18/08
!            RHO(KT,I)     = DENSITY(T1(KT,I),MAX(TDS(KT,I),0.0),MAX(TISS(KT,I),0.0))    ! SR 5/15/06
                    if (sdfirstadd(kt, i)) then
                        sed1(kt, i) = sed1ic(kt, i) ! cb 6/17/17
                        sed2(kt, i) = sed2ic(kt, i) ! cb 6/17/17
                        sdfirstadd(kt, i) = .false.
                    end if
                    do K = KT + 1, KMX
                        DEPTHB(K, I) = DEPTHB(K - 1, I) + H1(K, I)
                        DEPTHM(K, I) = DEPTHM(K - 1, I) + (H1(K - 1, I) + H1(K, I))*0.5D0
                    end do
                    C1(KT, I, CN(1:NAC)) = C1(KT + 1, I, CN(1:NAC))
                    CSSK(KT, I, CN(1:NAC)) = CSSK(KT + 1, I, CN(1:NAC))
                    KF(KT, I, KFCN(1:NAF(JW), JW)) = KF(KT + 1, I, KFCN(1:NAF(JW), JW))
                    KFS(KT, I, KFCN(1:NAF(JW), JW)) = KFS(KT + 1, I, KFCN(1:NAF(JW), JW)) !KF(KT+1,I,KFCN(1:NAF(JW),JW))   CODE ERROR FIX SW 10/24/2017
                    KF(KT + 1, I, KFCN(1:NAF(JW), JW)) = 0.0
                    KFS(KT + 1, I, KFCN(1:NAF(JW), JW)) = 0.0
                    if (KT >= KBI(I)) then ! CB 5/24/06
                        ADX(KT + 1, I) = 0.0 ! CB 5/15/06
                        C1(KT + 1, I, CN(1:NAC)) = 0.0 ! CB 5/15/06
                        CSSK(KT + 1, I, CN(1:NAC)) = 0.0 ! CB 5/15/06
                    end if ! CB 5/15/06
                    do JE = 1, NEP
                        if (kt < kbi(i)) then ! CB 4/28/06
                            EPD(KT, I, JE) = EPD(KT + 1, I, JE)
                            EPM(KT, I, JE) = EPD(KT, I, JE)*(Bi(KT, I) - B(KT + 1, I) + 2.0*H1(KT, I))*DLX(I) ! SR 5/15/06
                            EPM(KT + 1, I, JE) = EPM(KT + 1, I, JE) - EPM(KT, I, JE)
                            EPC(KT, I, JE) = EPM(KT, I, JE)/VOL(KT, I)
                            EPC(KT + 1, I, JE) = EPM(KT + 1, I, JE)/VOL(KT + 1, I)
                        else
                            EPD(KT, I, JE) = EPD(KT + 1, I, JE) ! SW 5/15/06
                            EPM(KT, I, JE) = EPM(KT + 1, I, JE)
                            EPC(KT, I, JE) = EPC(KT + 1, I, JE)
                            EPD(KT + 1, I, JE) = 0.0
                            EPM(KT + 1, I, JE) = 0.0
                            EPC(KT + 1, I, JE) = 0.0
                        end if ! CB 4/28/06
                    end do
                end do

!********macrophytes...
                do I = IU, ID
                    JT = KTI(I)
                    JE = KB(I)
                    do J = JT, JE
                        if (J < KT) then
                            COLB = EL(J + 1, I)
                        else
                            COLB = EL(KT + 1, I)
                        end if
!     COLDEP=ELWS(I)-COLB
                        coldep = EL(KT, i) - Z(i)*COSA(JB) - colb ! cb 3/7/16
!     MACT(J,KT,I)=MACT(J,KT+1,I)
                        if (MACROPHYTE_ON) then
                            MACT(J, KT, I) = MACT(J, KT + 1, I)
                        end if ! SW 9/28/13
                        do M = 1, NMC
                            if (MACROPHYTE_CALC(JW, M)) then
                                MACRC(J, KT, I, M) = MACRC(J, KT + 1, I, M)
                                MACRM(J, KT, I, M) = MACRC(J, KT, I, M)*CW(J, I)*COLDEP*DLX(I)
                            end if
                        end do
                    end do

                    JT = KT + 1
                    JE = KB(I)
                    do J = JT, JE
                        do M = 1, NMC
                            if (MACROPHYTE_CALC(JW, M)) then
                                MACRM(J, KT + 1, I, M) = MACRC(J, KT + 1, I, M)*CW(J, I)*H(KT + 1, JW)*DLX(I)
                            end if
                        end do
                    end do
                end do
                do I = IU - 1, ID
                    AVHR(KT + 1, I) = H1(KT + 1, I) + (H1(KT + 1, I + 1) - H1(KT + 1, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04 (H1(KT+1,I+1) +H1(KT+1,I))*0.5
                    AVHR(KT, I) = H1(KT, I) + (H1(KT, I + 1) - H1(KT, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04 (H1(KT,I+1)   +H1(KT,I))*0.5
                    BHR1(KT, I) = BH1(KT, I) + (BH1(KT, I + 1) - BH1(KT, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04  (BH1(KT,I+1)  +BH1(KT,I))*0.5
                    BHR1(KT + 1, I) = BH1(KT + 1, I) + (BH1(KT + 1, I + 1) - BH1(KT + 1, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04  (BH1(KT+1,I+1)+BH1(KT+1,I))*0.5  
                    if (CONSTRICTION(KT, I)) then ! SW 6/26/2018 Valid for all K
                        if (BHR1(KT, I) > BCONSTRICTION(I)*H1(KT, I)) then
                            BHR1(KT, I) = BCONSTRICTION(I)*H1(KT, I)
                        end if
                        if (BHR1(KT + 1, I) > BCONSTRICTION(I)*H1(KT + 1, I)) then
                            BHR1(KT + 1, I) = BCONSTRICTION(I)*H1(KT + 1, I)
                        end if
                    end if
                    U(KT, I) = U(KT + 1, I)
                end do
                do I = IU, ID
                    if (ONE_LAYER(I)) then
                        W(KT, I) = 0.0
                    else
                        W1 = W(KT + 1, I)*BB(KT + 1, I)
                        W2 = (BHR1(KT + 1, I)*U(KT + 1, I) - BHR1(KT + 1, I - 1)*U(KT + 1, I - 1))/DLX(I)
                        W3 = (-QSS(KT + 1, I))*BH1(KT + 1, I)/(BH1(KT + 1, I) + BH1(KT, I))/DLX(I)
                        W(KT, I) = (W1 + W2 + W3)/BB(KT, I)
                    end if
                end do
                if (UP_HEAD(JB)) then
                    BHSUM = BHR1(KT, IU - 1) + BHR1(KT + 1, IU - 1)
                    QUH1(KT, JB) = QUH1(KT + 1, JB)*BHR1(KT, IU - 1)/BHSUM
                    QUH1(KT + 1, JB) = QUH1(KT + 1, JB)*BHR1(KT + 1, IU - 1)/BHSUM
                    TSSUH1(KT, JB) = TSSUH1(KT + 1, JB)*BHR1(KT, IU - 1)/BHSUM
                    TSSUH1(KT + 1, JB) = TSSUH1(KT + 1, JB)*BHR1(KT + 1, IU - 1)/BHSUM
                    CSSUH1(KT, CN(1:NAC), JB) = CSSUH1(KT + 1, CN(1:NAC), JB)*BHR1(KT, IU - 1)/BHSUM
                    CSSUH1(KT + 1, CN(1:NAC), JB) = CSSUH1(KT + 1, CN(1:NAC), JB)*BHR1(KT + 1, IU - 1)/BHSUM
                end if
                if (DN_HEAD(JB)) then
                    BHSUM = BHR1(KT, ID) + BHR1(KT + 1, ID)
                    QDH1(KT, JB) = QDH1(KT + 1, JB)*BHR1(KT, ID)/BHSUM
                    QDH1(KT + 1, JB) = QDH1(KT + 1, JB)*BHR1(KT + 1, ID)/BHSUM
                    TSSDH1(KT, JB) = TSSDH1(KT + 1, JB)*BHR1(KT, ID)/BHSUM
                    TSSDH1(KT + 1, JB) = TSSDH1(KT + 1, JB)*BHR1(KT + 1, ID)/BHSUM
                    CSSDH1(KT, CN(1:NAC), JB) = CSSDH1(KT + 1, CN(1:NAC), JB)*BHR1(KT, ID)/BHSUM
                    CSSDH1(KT + 1, CN(1:NAC), JB) = CSSDH1(KT + 1, CN(1:NAC), JB)*BHR1(KT + 1, ID)/BHSUM
                end if
                do I = IU, ID - 1
                    if (DXI(JW) >= 0.0) then
                        DX(KT, I) = DXI(JW)
                    else
                        DX(KT, I) = ABS(U(KT, I))*ABS(DXI(JW))*H(K, JW) ! SW 8/2/2017
                    end if
                    if (INTERNAL_WEIR(KT, I)) then
                        DX(KT, I) = 0.0
                    end if
                end do
                IUT = IU
                IDT = ID - 1
                if (UP_HEAD(JB)) then
                    IUT = IU - 1
                end if
                if (DN_HEAD(JB)) then
                    IDT = ID
                end if
                do I = IUT, IDT
                    AZ(KT, I) = AZMIN
                    TKE(KT, I, 1) = 1.25E-7 !sg 10/4/07
                    TKE(KT, I, 2) = 1.0E-9 !sg 10/4/07
                    SAZ(KT, I) = AZMIN
                    if (INTERNAL_WEIR(KT, I)) then
                        AZ(KT, I) = 0.0
                        TKE(KT, I, 1) = 0.0 !sg  10/4/07
                        TKE(KT, I, 2) = 0.0 !sg  10/4/07
                        SAZ(KT, I) = 0.0
                    end if
                end do
!IF (CONSTITUENTS) THEN     ! SW 6/26/2019 No need to call twice - just do at end of BRANCH loop
!  CALL TEMPERATURE_RATES
!  CALL KINETIC_RATES
!END IF

!******** Upstream active segment

                IUT = US(JB)
                if (SLOPE(JB) == 0.0) then
                    do I = US(JB), DS(JB)
                        if (KB(I) - KT < NL(JB) - 1) then
                            IUT = I + 1
                        end if
                    end do
                else
                    do I = US(JB) - 1, DS(JB) + 1
                        if (KB(I) > KBI(I)) then
                            BNEW(KB(I), I) = B(KB(I), I) ! SW 1/23/06   ! SW 3/2/05
                            DX(KB(I), I) = 0.0
                            KB(I) = KB(I) - 1
                            ILAYER(I) = 1
                            U(KB(I) + 1, I) = 0.0 ! SW 1/23/06
                            write(WRN, "(2(A,I8),A,F0.3)") "Raising bottom layer at segment ", I, " at iteration ", NIT, " at Julian day ", JDAY
                            WARNING_OPEN = .true.
                        end if
                    end do ! SW 1/23/06
                    do I = US(JB) - 1, DS(JB) + 1 ! SW 1/23/06
!              IF (KB(I)-KT < NL(JB)-1) IUT = I+1    ! SW 1/23/06
!                IF (I /= DS(JB)+1) KBMIN(I)   = MIN(KB(I),KB(I+1))                    ! SW 1/23/06                                 ! SW 3/2/05
                        if (I /= US(JB) - 1) then
                            KBMIN(I - 1) = MIN(KB(I - 1), KB(I))
                        end if ! SW 1/23/06
                        if (KBI(I) < KB(I)) then
                            BKT(I) = BH1(KT, I)/(H1(KT, I) - (EL(KBI(I) + 1, I) - EL(KB(I) + 1, I))) ! SW 1/23/06
                            DEPTHB(KTWB(JW), I) = H1(KTWB(JW), I) - (EL(KBI(I) + 1, I) - EL(KB(I) + 1, I)) ! SW 1/23/06
                            DEPTHM(KTWB(JW), I) = (H1(KTWB(JW), I) - (EL(KBI(I) + 1, I) - EL(KB(I) + 1, I)))*0.5 ! SW 1/23/06
                            if (I <= DS(JB)) then ! SW 8/6/2018
                                AVHR(KT, I) = H1(KT, I) - (EL(KBI(I) + 1, I) - EL(KB(I) + 1, I)) + (H1(KT, I + 1) - (EL(KBI(I) + 1, I + 1) - EL(KB(I) + 1, I + 1)) - H1(KT, I) + EL(KBI(I) + 1, I) - EL(KB(I) + 1, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I)
                            else
                                AVHR(KT, I) = AVHR(KT, I - 1)
                            end if
! SW 1/23/06
                        end if
                    end do
                    do I = US(JB) - 1, DS(JB) + 1 ! SW 1/23/06
                        do K = KBMIN(I) + 1, KB(I)
                            U(K, I) = 0.0
                        end do
                    end do ! SW 11/9/07

                    do I = US(JB), DS(JB) ! SW 11/9/07
                        if (ILAYER(I) == 1 .and. ILAYER(I + 1) == 0) then
                            BHRSUM = 0.0
                            Q(I) = 0.0
                            do K = KT, KBMIN(I)
                                if (.not. INTERNAL_WEIR(K, I)) then
                                    BHRSUM = BHRSUM + BHR1(K, I)
                                    Q(I) = Q(I) + U(K, I)*BHR1(K, I)
                                end if
                            end do
                            do K = KT, KBMIN(I)
                                if (INTERNAL_WEIR(K, I)) then
                                    U(K, I) = 0.0
                                else
                                    U(K, I) = U(K, I) + (QC(I) - Q(I))/BHRSUM
                                end if
                            end do
                        else
                            if (ILAYER(I) == 1 .and. ILAYER(I - 1) == 0) then
                                BHRSUM = 0.0
                                Q(I - 1) = 0.0
                                do K = KT, KBMIN(I - 1)
                                    if (.not. INTERNAL_WEIR(K, I - 1)) then
                                        BHRSUM = BHRSUM + BHR1(K, I - 1)
                                        Q(I - 1) = Q(I - 1) + U(K, I - 1)*BHR1(K, I - 1)
                                    end if
                                end do
                                do K = KT, KBMIN(I - 1)
                                    if (INTERNAL_WEIR(K, I - 1)) then
                                        U(K, I - 1) = 0.0
                                    else
                                        U(K, I - 1) = U(K, I - 1) + (QC(I - 1) - Q(I - 1))/BHRSUM
                                    end if
                                end do
                            end if
                        end if

                    end do
                end if

!******** Segment addition

                if (IUT /= IU) then
                    if (SNAPSHOT(JW)) then
                        write(SNP(JW), "(/17X,2(A,I0))") " Add segments ", IUT, " through ", IU - 1
                    end if
                    WARNING_OPEN = .true.
                    write(WRN, "(/17X,2(A,I0))") " Add segments ", IUT, " through ", IU - 1
                    do I = IUT - 1, IU - 1
                        Z(I) = Z(IU)
                        KTI(I) = KTI(IU)
                        H1(KT + 1, I) = H(KT + 1, JW)
                        AVH1(KT + 1, I) = (H1(KT + 1, I) + H1(KT + 2, I))*0.5
                        AVHR(KT + 1, I) = H1(KT + 1, I) + (H1(KT + 1, I + 1) - H1(KT + 1, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04
                        if (.not. TRAPEZOIDAL(JW)) then
                            BH1(KT + 1, I) = B(KT + 1, I)*H(KT + 1, JW)
                            H1(KT, I) = H(KT, JW) - Z(I)
                            BI(KT:KB(I), I) = B(KT:KB(I), I) ! SW 4/18/07
                            BI(KT, I) = B(KTI(I), I)
                            AVH1(KT, I) = (H1(KT, I) + H1(KT + 1, I))*0.5
                            AVHR(KT, I) = H1(KT, I) + (H1(KT, I + 1) - H1(KT, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04
                            BH1(KT, I) = BI(KT, I)*(EL(KT, I) - Z(I)*COSA(JB) - EL(KTI(I) + 1, I))/COSA(JB)
                            if (KTI(I) >= KB(I)) then
                                BH1(KT, I) = B(KT, I)*H1(KT, I)
                            end if
                            do K = KTI(I) + 1, KT
                                BH1(KT, I) = BH1(KT, I) + BH(K, I)
                            end do
                        else
                            call GRID_AREA1(EL(KT, I) - Z(I), EL(KT + 1, I), BH1(KT, I), BI(KT, I)) !SW 08/03/04
                            BH1(KT + 1, I) = 0.25*H(KT + 1, JW)*(BB(KT, I) + 2.*B(KT + 1, I) + BB(KT + 1, I))
                            H1(KT, I) = H(KT, JW) - Z(I)
                            AVH1(KT, I) = (H1(KT, I) + H1(KT + 1, I))*0.5
                            AVHR(KT, I) = H1(KT, I) + (H1(KT, I + 1) - H1(KT, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04
                        end if
                        BKT(I) = BH1(KT, I)/H1(KT, I)
                        DEPTHB(KT, I) = H1(KT, I)
                        DEPTHM(KT, I) = H1(KT, I)*0.5
                        do K = KT + 1, KB(I)
                            DEPTHB(K, I) = DEPTHB(K - 1, I) + H1(K, I)
                            DEPTHM(K, I) = DEPTHM(K - 1, I) + (H1(K - 1, I) + H1(K, I))*0.5
                        end do
                    end do
                    do I = IUT - 1, IU - 1
                        BHR1(KT + 1, I) = BH1(KT + 1, I) + (BH1(KT + 1, I + 1) - BH1(KT + 1, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04
                        BHR1(KT, I) = BH1(KT, I) + (BH1(KT, I + 1) - BH1(KT, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04
                        if (CONSTRICTION(KT, I)) then ! SW 6/26/2018 Valid for all K
                            if (BHR1(KT, I) > BCONSTRICTION(I)*H1(KT, I)) then
                                BHR1(KT, I) = BCONSTRICTION(I)*H1(KT, I)
                            end if
                            if (BHR1(KT + 1, I) > BCONSTRICTION(I)*H1(KT + 1, I)) then
                                BHR1(KT + 1, I) = BCONSTRICTION(I)*H1(KT + 1, I)
                            end if
                        end if
                    end do
                    do I = IUT, IU - 1
!ICE(I)   = ICE(IU)   ! SW 9/29/15
!ICETH(I) = ICETH(IU)
                        WIND2(I) = WIND2(IU)
                        if (DYNAMIC_SHADE(I)) then
                            call SHADING()
                        end if
                        do K = KT, KBMIN(I) !KB(I)   SW 12/18/2018
                            U(K, I) = U(K, IU)
                            if (DXI(JW) >= 0.0) then
                                DX(K, I) = DXI(JW)
                            else
                                DX(K, I) = ABS(U(K, I))*ABS(DXI(JW))*H(K, JW) ! SW 8/2/2017
                            end if
                            SDKV(K, I) = SDK(JW) ! SW 1/18/08
                            if (INTERNAL_WEIR(K, I)) then
                                DX(K, I) = 0.0
                                U(K, I) = 0.0
                            end if
                        end do

                        do K = KT, KB(I) ! SW 12/18/2018
                            T1(K, I) = T1(K, IU)
                            T2(K, I) = T1(K, IU)
                            SU(K, I) = U(K, IU)
                            C1(K, I, CN(1:NAC)) = C1(K, IU, CN(1:NAC))
                            C2(K, I, CN(1:NAC)) = C1(K, IU, CN(1:NAC))
                            do JE = 1, NEP
                                EPD(K, I, JE) = 0.01
                                EPC(K, I, JE) = 0.01/H1(K, I)
                            end do
                            CMBRT(CN(1:NAC), JB) = CMBRT(CN(1:NAC), JB) + C1(K, IU, CN(1:NAC))*DLX(I)*BH1(K, I)
                            EBRI(JB) = EBRI(JB) + T1(K, IU)*DLX(I)*BH1(K, I)
                        end do
                        do K = KT, KB(I) - 1
                            AZ(K, I) = AZ(K, IU)
                            TKE(K, I, 1) = TKE(K, IU, 1) !sg 10/4/07
                            TKE(K, I, 2) = TKE(K, IU, 2) !sg 10/4/07
                            SAZ(K, I) = AZ(K, IU)
                            if (INTERNAL_WEIR(K, I)) then
                                AZ(K, I) = 0.0
                                TKE(K, I, 1) = 0.0 !sg 10/4/07
                                TKE(K, I, 2) = 0.0 !sg 10/4/07
                                SAZ(K, I) = 0.0
                            end if
                        end do
!END DO
!*********macrophytes
                        do M = 1, NMC
                            if (MACROPHYTE_CALC(JW, M)) then
                                JT = KTI(I)
                                JE = KB(I)
                                do J = JT, JE
                                    if (J < KT) then
                                        COLB = EL(J + 1, I)
                                    else
                                        COLB = EL(KT + 1, I)
                                    end if
!COLDEP=ELWS(I)-COLB
                                    coldep = EL(KT, i) - Z(i)*COSA(JB) - colb ! cb 3/7/16
!MACRC(J,KT,I,M)=MACWBCI(JW,M)
                                    if (ISO_macrophyte(JW, m)) then
                                        macrc(j, kt, I, m) = macwbci(JW, m)
                                    end if ! cb 3/7/16
                                    if (VERT_macrophyte(JW, m)) then
                                        macrc(j, kt, I, m) = 0.1
                                    end if
                                    if (long_macrophyte(JW, m)) then
                                        macrc(j, kt, I, m) = 0.1
                                    end if
                                    MACRM(J, KT, I, M) = MACRC(J, KT, I, M)*CW(J, I)*COLDEP*DLX(I)
                                    MACMBRT(JB, M) = MACMBRT(JB, M) + MACRM(J, KT, I, M)
                                end do
                                do K = KT + 1, KB(I)
                                    JT = K
                                    JE = KB(I)
                                    do J = JT, JE
!MACRC(J,K,I,M)=MACWBCI(JW,M)                      
                                        if (ISO_macrophyte(JW, m)) then
                                            macrc(j, k, I, m) = macwbci(JW, m)
                                        end if ! cb 3/7/16
                                        if (VERT_macrophyte(JW, m)) then
                                            macrc(j, k, I, m) = 0.1
                                        end if
                                        if (long_macrophyte(JW, m)) then
                                            macrc(j, k, I, m) = 0.1
                                        end if
                                        MACRM(J, K, I, M) = MACRC(J, K, I, M)*CW(J, I)*H2(K, I)*DLX(I)
                                        MACMBRT(JB, M) = MACMBRT(JB, M) + MACRM(J, K, I, M)
                                    end do
                                end do
                            end if
                        end do
                    end do ! cb 3/7/16  moved enddo to include macrophytes
                    U(KB(IUT):KB(IU), IU - 1) = 0.0
                    SU(KB(IUT):KB(IU), IU - 1) = 0.0
                    ADX(KB(IUT):KB(IU), IU) = 0.0
                    IU = IUT
                    CUS(JB) = IU
                    if (UH_EXTERNAL(JB)) then
                        KB(IU - 1) = KB(IU)
                    end if
                    if (UH_INTERNAL(JB)) then
                        if (JBUH(JB) >= BS(JW) .and. JBUH(JB) <= BE(JW)) then
                            KB(IU - 1) = MIN(KB(UHS(JB)), KB(IU))
                        else
                            do KKB = KT, KMX
                                if (EL(KKB, IU) <= EL(KB(UHS(JB)), UHS(JB))) then
                                    exit
                                end if
                            end do
                            KB(IU - 1) = MIN(KKB, KB(IU))
                        end if
                    end if
                    if (UP_HEAD(JB)) then
                        AZ(KT:KB(IU - 1) - 1, IU - 1) = AZMIN
                        TKE(KT:KB(IU - 1) - 1, IU - 1, 1) = 1.25E-7 !SG 10/4/07
                        TKE(KT:KB(IU - 1) - 1, IU - 1, 2) = 1.0E-9 !SG 10/4/07
                        SAZ(KT:KB(IU - 1) - 1, IU - 1) = AZMIN
                    end if
                end if
                if (CONSTITUENTS) then
                    call TEMPERATURE_RATES()
                    call KINETIC_RATES()
                end if

!******** Total active cells and single layers

                do I = IU, ID
                    NTAC = NTAC + 1
                    ONE_LAYER(I) = KTWB(JW) == KB(I)
                end do
                NTACMX = MAX(NTAC, NTACMX)
            end do
            call INTERPOLATION_MULTIPLIERS()

!****** Additional layers

            ZMIN(JW) = -1000.0
            do JB = BS(JW), BE(JW)
                if (BR_INACTIVE(JB)) then
                    cycle
                end if ! SW 6/12/2017
                do I = CUS(JB), DS(JB)
                    ZMIN(JW) = MAX(ZMIN(JW), Z(I))
                end do
            end do
            ADD_LAYER = ZMIN(JW) < (-0.80)*H(KT - 1, JW) .and. KT /= 2
        end do

!**** Subtract layers

        do while (SUB_LAYER)
            if (SNAPSHOT(JW)) then
                write(SNP(JW), '(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,1x,13("*"))') "Subtract layer ", KT, " at Julian day = ", JDAY, " NIT = ", NIT, " IZMIN =", IZMIN(JW)
            end if ! SW 1/23/06   
            write(WRN, '(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,1x,13("*"))') "Subtract layer ", KT, " at Julian day = ", JDAY, " NIT = ", NIT, " IZMIN =", IZMIN(JW)

!****** Variable initialization

            KTWB(JW) = KTWB(JW) + 1
            KT = KTWB(JW)
            ILAYER = 0 ! SW 1/23/06  11/7/07

! RECOMPUTE INTERNAL WEIR FOR FLOATING WEIR      
            if (WEIR_CALC) then !  SW 3/16/18
                do JWR = 1, NIW
                    if (IWR(JWR) >= US(BS(JW)) .and. IWR(JWR) <= DS(BE(JW))) then
                        if (EKTWR(JWR) == 0.0) then
                            KTWR(JWR) = KTWB(JW)
                        else
                            KTWR(JWR) = INT(EKTWR(JWR))
                        end if
                        if (EKBWR(JWR) <= 0.0) then
                            do K = KTWR(JWR), KB(IWR(JWR))
                                if (DEPTHB(K, IWR(JWR)) > ABS(EKBWR(JWR))) then
                                    KBWR(JWR) = K
                                    exit
                                end if
                            end do
                        else
                            KBWR(JWR) = INT(EKBWR(JWR))
                        end if

                        do K = 2, KMX - 1
                            if (K >= KTWR(JWR) .and. K <= KBWR(JWR)) then
                                INTERNAL_WEIR(K, IWR(JWR)) = .true.
                            end if
                        end do
                    end if
                end do
            end if



            do JB = BS(JW), BE(JW)
                if (BR_INACTIVE(JB)) then
                    cycle
                end if ! SW 6/12/2017
                IU = CUS(JB)
                ID = DS(JB)
                if (CONSTITUENTS) then
                    DO1(KT - 1, IU - 1:ID + 1) = 0.0
                end if
                do I = IU - 1, ID + 1
                    Z(I) = Z(I) - H(KT - 1, JW)
                    H1(KT - 1, I) = H(KT - 1, JW)
                    H1(KT, I) = H(KT, JW) - Z(I)
                    BI(KT, I) = B(KTI(I), I)
                    BI(KT - 1, I) = B(KT - 1, I)
                    AVH1(KT - 1, I) = (H(KT - 1, JW) + H(KT, JW))*0.5
                    AVH1(KT, I) = (H1(KT, I) + H1(KT + 1, I))*0.5
                    if (.not. TRAPEZOIDAL(JW)) then
                        if (KB(I) >= KT) then ! SW 1/23/06
                            BH1(KT, I) = BH1(KT - 1, I) + BH1(KT, I) ! SW 1/23/06
                            BH1(KT - 1, I) = B(KT - 1, I)*H(KT - 1, JW)
                        else ! SW 1/23/06
                            BH1(KT, I) = BH1(KT - 1, I) ! SW 1/23/06
                        end if ! SW 1/23/06
                    else
                        call GRID_AREA1(EL(KT, I) - Z(I), EL(KT + 1, I), BH1(KT, I), DUMMY) !SW 08/03/04
                        BH1(KT - 1, I) = 0.25*H1(KT - 1, JW)*(BB(KT - 2, I) + 2.*B(KT - 1, I) + BB(KT - 1, I))
                    end if
                    VOL(KT, I) = BH1(KT, I)*DLX(I)
                    VOL(KT - 1, I) = BH1(KT - 1, I)*DLX(I)
                    BKT(I) = BH1(KT, I)/H1(KT, I)
                    if (kb(i) >= kt) then ! SW 1/23/06
                        U(KT, I) = (U(KT - 1, I)*BHR1(KT - 1, I) + U(KT, I)*BHR(KT, I))/(BHR1(KT - 1, I) + BHR(KT, I))
                        T1(KT, I) = (T1(KT - 1, I)*(BH1(KT, I) - BH(KT, I)) + T1(KT, I)*BH(KT, I))/BH1(KT, I)
                    else
!              EBRI(JB) = EBRI(JB)-T1(KT,I)*VOL(KT,I)   1/23/06
                        u(kt, i) = u(kt - 1, i) ! SW 1/23/06
                        t1(kt, i) = t1(kt - 1, i) ! SW 1/23/06
                    end if
                    if (KB(I) >= KT) then ! SW 1/23/06
                        C1(KT, I, CN(1:NAC)) = (C1(KT - 1, I, CN(1:NAC))*(BH1(KT, I) - BH(KT, I)) + C1(KT, I, CN(1:NAC))*BH(KT, I))/BH1(KT, I)
                        CSSK(KT, I, CN(1:NAC)) = (CSSK(KT - 1, I, CN(1:NAC))*(BH1(KT, I) - BH(KT, I)) + CSSK(KT, I, CN(1:NAC))*BH(KT, I))/BH1(KT, I)
                    else
                        C1(KT, I, CN(1:NAC)) = C1(KT - 1, I, CN(1:NAC)) ! SW 1/23/06
                        CSSK(KT, I, CN(1:NAC)) = CSSK(KT - 1, I, CN(1:NAC)) ! SW 1/23/06
                    end if ! SW 1/23/06
                    CSSB(KT, I, CN(1:NAC)) = CSSB(KT - 1, I, CN(1:NAC)) + CSSB(KT, I, CN(1:NAC))
!KF(KT,I,KFCN(1:NAF(JW),JW))    =  KF(KT-1,I,KFCN(1:NAF(JW),JW))
!KFS(KT,I,KFCN(1:NAF(JW),JW))   =  KFS(KT-1,I,KFCN(1:NAF(JW),JW))
                    KF(KT, I, KFCN(1:NAF(JW), JW)) = (KF(KT - 1, I, KFCN(1:NAF(JW), JW))*VOL(KT - 1, I) + KF(KT, I, KFCN(1:NAF(JW), JW))*VOL(KT, I))/(VOL(KT - 1, I) + VOL(KT, I)) ! SW Fix suggested by Taylor Adams Hydros 25Feb2021  ! KF is in units of g/m3/s
                    KFS(KT, I, KFCN(1:NAF(JW), JW)) = KFS(KT - 1, I, KFCN(1:NAF(JW), JW)) + KFS(KT, I, KFCN(1:NAF(JW), JW)) ! SW Fix suggested by Taylor Adams Hydros 25Feb2021 ! KFS is in units of g KFS=KF*VOL*DT          C1(KT-1,I,CN(1:NAC))           =  0.0
                    C1(KT - 1, I, CN(1:NAC)) = 0.0
                    C2(KT - 1, I, CN(1:NAC)) = 0.0
                    CSSB(KT - 1, I, CN(1:NAC)) = 0.0
                    CSSK(KT - 1, I, CN(1:NAC)) = 0.0
                    KF(KT - 1, I, KFCN(1:NAF(JW), JW)) = 0.0
                    KFS(KT - 1, I, KFCN(1:NAF(JW), JW)) = 0.0
                    do JE = 1, NEP
                        if (KT <= KBI(I)) then ! CB 4/28/06
                            EPM(KT, I, JE) = EPM(KT - 1, I, JE) + EPM(KT, I, JE)
                            EPD(KT, I, JE) = EPM(KT, I, JE)/((BI(KT, I) - BI(KT + 1, I) + 2.0*H1(KT, I))*DLX(I))
                            EPC(KT, I, JE) = EPM(KT, I, JE)/VOL(KT, I)
                            EPM(KT - 1, I, JE) = 0.0
                            EPD(KT - 1, I, JE) = 0.0
                            EPC(KT - 1, I, JE) = 0.0
                        else ! SW 5/15/06
                            EPM(KT, I, JE) = EPM(KT - 1, I, JE)
                            EPD(KT, I, JE) = EPD(KT - 1, I, JE)
                            EPC(KT, I, JE) = EPC(KT - 1, I, JE)
                            EPM(KT - 1, I, JE) = 0.0
                            EPD(KT - 1, I, JE) = 0.0
                            EPC(KT - 1, I, JE) = 0.0
                        end if ! CB 4/28/06
                    end do
                end do
                do I = IU, ID
                    do M = 1, NMC
                        if (MACROPHYTE_CALC(JW, M)) then
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
                                if (J < KT) then
                                    MACRM(J, KT, I, M) = MACRM(J, KT - 1, I, M)
                                else
                                    MACRM(J, KT, I, M) = MACRM(J, KT - 1, I, M) + MACRM(J, KT, I, M)
                                end if
                                if (CW(J, I) > 0.0) then
                                    MACRC(J, KT, I, M) = MACRM(J, KT, I, M)/(CW(J, I)*COLDEP*DLX(I))
                                else
                                    MACRC(J, KT, I, M) = 0.0
                                end if
                                MACRM(J, KT - 1, I, M) = 0.0
                                MACRC(J, KT - 1, I, M) = 0.0

                            end do
                        end if
                    end do
                    JT = KTI(I)
                    JE = KB(I)
                    do J = JT, JE
                        MACT(J, KT, I) = 0.0
                        MACT(J, KT - 1, I) = 0.0
                    end do
                    do M = 1, NMC
                        if (MACROPHYTE_CALC(JW, M)) then
                            do J = JT, JE
                                MACT(J, KT, I) = MACRC(J, KT, I, M) + MACT(J, KT, I)
                            end do
                        end if
                    end do
                    do M = 1, NMC
                        TMAC = 0.0
                        if (MACROPHYTE_CALC(JW, M)) then
                            JT = KTI(I)
                            JE = KB(I)
                            do J = JT, JE
                                TMAC = TMAC + MACRM(J, KT, I, M)
                            end do
                        end if
                        MAC(KT, I, M) = TMAC/(BH1(KT, I)*DLX(I))
                    end do
                    do M = 1, NMC
                        if (MACROPHYTE_CALC(JW, M)) then
                            MAC(KT - 1, I, M) = 0.0
                        end if
                    end do
                end do

                do I = IU - 1, ID
                    AVHR(KT - 1, I) = H1(KT - 1, I) + (H1(KT - 1, I + 1) - H1(KT - 1, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04
                    AVHR(KT, I) = H1(KT, I) + (H1(KT, I + 1) - H1(KT, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04
                    BHR1(KT, I) = BH1(KT, I) + (BH1(KT, I + 1) - BH1(KT, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04
                    BHR1(KT - 1, I) = BH1(KT - 1, I) + (BH1(KT - 1, I + 1) - BH1(KT - 1, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04
                    if (CONSTRICTION(KT, I)) then ! SW 6/26/2018 Valid for all K
                        if (BHR1(KT, I) > BCONSTRICTION(I)*H1(KT, I)) then
                            BHR1(KT, I) = BCONSTRICTION(I)*H1(KT, I)
                        end if
                        if (BHR1(KT - 1, I) > BCONSTRICTION(I)*H1(KT - 1, I)) then
                            BHR1(KT - 1, I) = BCONSTRICTION(I)*H1(KT - 1, I)
                        end if
                    end if
                end do
                U(KT - 1, IU - 1:ID + 1) = 0.0
                W(KT - 1, IU - 1:ID + 1) = 0.0
                SW(KT - 1, IU - 1:ID + 1) = 0.0 !TC 3/9/05
                P(KT - 1, IU - 1:ID + 1) = 0.0
                AZ(KT - 1, IU - 1:ID + 1) = 0.0
                TKE(KT - 1, IU - 1:ID + 1, 1) = 0.0 !sg 10/4/07
                TKE(KT - 1, IU - 1:ID + 1, 2) = 0.0 !sg 10/4/07
                DZ(KT - 1, IU - 1:ID + 1) = 0.0
                ADMZ(KT - 1, IU - 1:ID + 1) = 0.0
                ADZ(KT - 1, IU - 1:ID + 1) = 0.0
                DECAY(KT - 1, IU - 1:ID + 1) = 0.0
                if (UP_HEAD(JB)) then
                    QUH1(KT, JB) = QUH1(KT, JB) + QUH1(KT - 1, JB)
                    TSSUH1(KT, JB) = TSSUH1(KT - 1, JB) + TSSUH1(KT, JB)
                    CSSUH1(KT, CN(1:NAC), JB) = CSSUH1(KT - 1, CN(1:NAC), JB) + CSSUH1(KT, CN(1:NAC), JB)
                end if
                if (DN_HEAD(JB)) then
                    QDH1(KT, JB) = QDH1(KT, JB) + QDH1(KT - 1, JB)
                    TSSDH1(KT, JB) = TSSDH1(KT - 1, JB) + TSSDH1(KT, JB)
                    CSSDH1(KT, CN(1:NAC), JB) = CSSDH1(KT - 1, CN(1:NAC), JB) + CSSDH1(KT, CN(1:NAC), JB)
                end if

!******** Upstream active segment

                IUT = US(JB)
                if (SLOPE(JB) /= 0.0) then
                    do I = US(JB) - 1, DS(JB) + 1
                        if (KB(I) < KT) then ! SR 10/17/05
                            KB(I) = KT
                            Bnew(KB(I), I) = 0.000001 ! sw 1/23/06
                            if (DXI(JW) >= 0.0) then
                                DX(KB(I), I) = DXI(JW)
                            else
                                DX(KB(I), I) = ABS(U(KB(I), I))*ABS(DXI(JW))*H(K, JW) ! SW 8/2/2017
                            end if
                            ilayer(i) = 1
                            T1(KB(I), I) = T1(KT, I) !    SW 5/15/06    T1(KB(I)-1,I)
                            C1(KB(I), I, CN(1:NAC)) = C1(KT, I, CN(1:NAC)) !    SW 5/15/06    C1(KB(I)-1,I,CN(1:NAC))

                            if (SEDIMENT_CALC(JW)) then ! SW 5/26/2022
                                SED(KB(I), I) = SED(KT - 1, I);                                 SEDC(KB(I), I) = SEDC(KT - 1, I);                                 SEDN(KB(I), I) = SEDN(KT - 1, I);                                 SEDP(KB(I), I) = SEDP(KT - 1, I)
                            end if

                            write(WRN, "(2(A,I8),A,F0.3,A,F0.3)") "Lowering bottom segment ", I, " at iteration ", NIT, " at Julian day ", JDAY, " Z(I)=", Z(I)
                            WARNING_OPEN = .true.
                        end if
                    end do ! SW 1/23/06
                    do I = US(JB) - 1, DS(JB) + 1 ! SW 1/23/06
!               IF (I /= DS(JB)+1) KBMIN(I)   = MIN(KB(I),KB(I+1))  ! SW 1/23/06
                        if (I /= US(JB) - 1) then
                            KBMIN(I - 1) = MIN(KB(I - 1), KB(I))
                        end if ! SW 1/23/06
                        if (KBI(I) < KB(I)) then
                            BKT(I) = BH1(KT, I)/(H1(KT, I) - (EL(KBI(I) + 1, I) - EL(KB(I) + 1, I))) ! SW 1/23/06
                            DEPTHB(KTWB(JW), I) = H1(KTWB(JW), I) - (EL(KBI(I) + 1, I) - EL(KB(I) + 1, I)) ! SW 1/23/06
                            DEPTHM(KTWB(JW), I) = (H1(KTWB(JW), I) - (EL(KBI(I) + 1, I) - EL(KB(I) + 1, I)))*0.5 ! SW 1/23/06
                            AVHR(KT, I) = H1(KT, I) - (EL(KBI(I) + 1, I) - EL(KB(I) + 1, I)) + (H1(KT, I + 1) - (EL(KBI(I) + 1, I + 1) - EL(KB(I) + 1, I + 1)) - H1(KT, I) + EL(KBI(I) + 1, I) - EL(KB(I) + 1, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) ! SW 1/23/06
                        end if
                    end do

                    do I = US(JB), DS(JB) ! SW 1/23/06   11/13/07   US(JB)-1,DS(JB)+1

                        if (ILAYER(I) == 1 .and. ILAYER(I + 1) == 0) then ! SW 1/23/06
                            BHRSUM = 0.0
                            Q(I) = 0.0
                            do K = KT, KBMIN(I)
                                if (.not. INTERNAL_WEIR(K, I)) then
                                    BHRSUM = BHRSUM + BHR1(K, I)
                                    Q(I) = Q(I) + U(K, I)*BHR1(K, I)
                                end if
                            end do
                            do K = KT, KBMIN(I)
                                if (INTERNAL_WEIR(K, I)) then
                                    U(K, I) = 0.0
                                else
                                    U(K, I) = U(K, I) + (QC(I) - Q(I))/BHRSUM
                                end if
                            end do
                        else
                            if (ILAYER(I) == 1 .and. ILAYER(I - 1) == 0) then
                                BHRSUM = 0.0
                                Q(I - 1) = 0.0
                                do K = KT, KBMIN(I - 1)
                                    if (.not. INTERNAL_WEIR(K, I - 1)) then
                                        BHRSUM = BHRSUM + BHR1(K, I - 1)
                                        Q(I - 1) = Q(I - 1) + U(K, I - 1)*BHR1(K, I - 1)
                                    end if
                                end do
                                do K = KT, KBMIN(I - 1)
                                    if (INTERNAL_WEIR(K, I - 1)) then
                                        U(K, I - 1) = 0.0
                                    else
                                        U(K, I - 1) = U(K, I - 1) + (QC(I - 1) - Q(I - 1))/BHRSUM
                                    end if
                                end do
                            end if
                        end if ! SW 1/23/06

                    end do
                end if
                do I = US(JB), DS(JB)
                    if (KB(I) - KT < NL(JB) - 1) then
                        IUT = I + 1
                    end if
                    ONE_LAYER(I) = KTWB(JW) == KB(I)
                end do
                if (IUT > DS(JB)) then
                    if (JB == 1) then
                        write(W2ERR, "(A,I0/A,F0.2,2(A,I0))") "Fatal error - insufficient segments in branch ", JB, "Julian day = ", JDAY, " at iteration ", NIT, " with water surface layer = ", KT
                        write(W2ERR, "(2(A,I0))") "Minimum water surface located at segment ", IZMIN(JW), " with bottom layer at ", KB(IZMIN(JW))
                        TEXT = "Runtime error - see w2.err"
                        ERROR_OPEN = .true.
                        return
                    else
                        BR_INACTIVE(JB) = .true. ! SW 6/12/2017
                        if (SNAPSHOT(JW)) then
                            write(SNP(JW), '(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,13("*"))') "   Branch Inactive: ", jb, " at Julian day = ", JDAY, "   NIT = ", NIT
                        end if
                        WARNING_OPEN = .true.
                        write(WRN, '(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,13("*"))') "   Branch Inactive: ", jb, " at Julian day = ", JDAY, "   NIT = ", NIT ! SW 11/16/2018
                        do I = IU, DS(JB) ! SW 12/17/2018
                            do K = KT - 1, KB(I) ! SW 12/18/2018  KT,KB(I)
                                EBRI(JB) = EBRI(JB) - T1(K, I)*BH1(K, I)*DLX(I) ! VOL(K,I)   SW 12/18/2018
                                CMBRT(CN(1:NAC), JB) = CMBRT(CN(1:NAC), JB) - C1(K, I, CN(1:NAC))*BH1(K, I)*DLX(I) !VOL(K,I)              !+(CSSB(K,I,CN(1:NAC))+CSSK(K,I,CN(1:NAC))*VOL(K,I))*DLT
                            end do
                        end do
                        cycle
                    end if

! Go to 230
                end if

!******** Segment subtraction

                if (IUT /= IU) then
                    if (SNAPSHOT(JW)) then
                        write(SNP(JW), "(/17X,A,I0,A,I0)") " Subtract segments ", IU, " through ", IUT - 1
                    end if
                    WARNING_OPEN = .true.
                    write(WRN, "(/17X,A,I0,A,I0)") " Subtract segments ", IU, " through ", IUT - 1
                    do I = IU, IUT - 1
                        do K = KT - 1, KB(I) ! SW 12/18/2018 KT,KB(I)
                            EBRI(JB) = EBRI(JB) - T1(K, I)*BH1(K, I)*DLX(I) !*VOL(K,I)    SW 12/18/2018
                            CMBRT(CN(1:NAC), JB) = CMBRT(CN(1:NAC), JB) - C1(K, I, CN(1:NAC))*BH1(K, I)*DLX(I) !*VOL(K,I)        !+(CSSB(K,I,CN(1:NAC))+CSSK(K,I,CN(1:NAC))*VOL(K,I))*DLT
                        end do
                    end do

                    do I = IU, IUT - 1
                        do M = 1, NMC
                            if (MACROPHYTE_CALC(JW, M)) then
                                JT = KTI(I)
                                JE = KB(I)
                                do J = JT, JE
                                    if (J < KT) then
                                        COLB = EL(J + 1, I)
                                    else
                                        COLB = EL(KT + 1, I)
                                    end if
!COLDEP=ELWS(I)-COLB
                                    coldep = EL(KT, i) - Z(i)*COSA(JB) - colb ! cb 3/7/16
!                    MACMBRT(JB,M) = MACMBRT(JB,M)-MACRM(J,KT,I,M)+(MACSS(J,KT,I,M)*COLDEP*CW(J,I)*DLX(I))*DLT
                                    MACMBRT(JB, M) = MACMBRT(JB, M) - MACRM(J, KT, I, M)
                                end do
                                do K = KT + 1, KB(I)
                                    JT = K
                                    JE = KB(I)
                                    do J = JT, JE
!                      MACMBRT(JB,M) = MACMBRT(JB,M)-MACRM(J,K,I,M)+(MACSS(J,K,I,M)*H2(K,I)*CW(J,I)*DMX(I))*DLT
                                        MACMBRT(JB, M) = MACMBRT(JB, M) - MACRM(J, K, I, M)
                                    end do
                                end do
                            end if
                        end do
                    end do

                    F(IU - 1:IUT - 1) = 0.0
                    Z(IU - 1:IUT - 1) = 0.0
!ICETH(IU-1:IUT-1) =  0.0    ! SW 9/29/15
                    BHRHO(IU - 1:IUT - 1) = 0.0
!ICE(IU-1:IUT-1)   = .FALSE.  ! SW 9/29/15
                    do K = KT, KB(IUT)
                        ADX(K, IU - 1:IUT - 1) = 0.0
                        DX(K, IU - 1:IUT - 1) = 0.0
                        AZ(K, IU - 1:IUT - 1) = 0.0
                        TKE(K, IU - 1:IUT - 1, 1) = 0.0 !SG  10/4/07
                        TKE(K, IU - 1:IUT - 1, 2) = 0.0 !SG  10/4/07
                        SAZ(K, IU - 1:IUT - 1) = 0.0
                        U(K, IU - 1:IUT - 1) = 0.0
                        SU(K, IU - 1:IUT - 1) = 0.0
                        T1(K, IU - 1:IUT - 1) = 0.0
                        TSS(K, IU - 1:IUT - 1) = 0.0
                        QSS(K, IU - 1:IUT - 1) = 0.0
                        C1(K, IU - 1:IUT - 1, CN(1:NAC)) = 0.0
                        C2(K, IU - 1:IUT - 1, CN(1:NAC)) = 0.0
                        C1S(K, IU - 1:IUT - 1, CN(1:NAC)) = 0.0
                        CSSB(K, IU - 1:IUT - 1, CN(1:NAC)) = 0.0
                        CSSK(K, IU - 1:IUT - 1, CN(1:NAC)) = 0.0
                    end do

                    do M = 1, NMC
                        if (MACROPHYTE_CALC(JW, M)) then
                            MAC(K, I, M) = 0.0
                            MACT(J, K, I) = 0.0
                        end if
                    end do
                    JT = KTI(I)
                    JE = KB(I)
                    do J = JT, JE
                        do M = 1, NMC
                            if (MACROPHYTE_CALC(JW, M)) then
                                MACRC(J, K, I, M) = 0.0
                            end if
                        end do
                    end do

                    IU = IUT
                    CUS(JB) = IU
                    Z(IU - 1) = (EL(KT, IU - 1) - (EL(KT, IU) - Z(IU)*COSA(JB)))/COSA(JB)
                    SZ(IU - 1) = Z(IU)
                    KTI(IU - 1) = KTI(IU)
                    if (.not. TRAPEZOIDAL(JW)) then
                        BI(KT, IU - 1) = B(KTI(IU - 1), I)
                        H1(KT, IU - 1) = H(KT, JW) - Z(IU - 1)
                        BH1(KT, IU - 1) = Bnew(KTI(IU - 1), IU - 1)*(EL(KT, IU - 1) - EL(KTI(IU - 1) + 1, IU - 1) - Z(IU - 1)*COSA(JB))/COSA(JB) ! sw 1/23/06  Bnew(KTI(IU-1),IU-1)*(EL(KT,IU-1)-EL(KTI(IU-1)+1,IU-1)-Z(IU-1)*COSA(JB))/COSA(JB)     ! SR 10/17/05
                        if (KT >= KB(IU - 1)) then
                            BH1(KT, IU - 1) = Bnew(KT, IU - 1)*H1(KT, IU - 1)
                        end if ! sw 1/23/06
                        do K = KTI(IU - 1) + 1, KT
                            BH1(KT, IU - 1) = BH1(KT, IU - 1) + BH1(K, IU - 1)
                        end do
                    else
                        call GRID_AREA1(EL(KT, I) - Z(I), EL(KT + 1, IU - 1), BH1(KT, IU - 1), BI(KT, IU - 1)) !SW 08/03/04
                        BH1(KT, I) = 0.25*H(KT, JW)*(BB(KT - 1, I) + 2.*B(KT, I) + BB(KT, I))
                        H1(KT, I) = H(KT, JW) - Z(I)
                    end if
                    BKT(IU - 1) = BH1(KT, IU - 1)/H1(KT, IU - 1)
                    BHR1(KT, IU - 1) = BH1(KT, IU - 1) + (BH1(KT, IU) - BH1(KT, IU - 1))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04
                    if (CONSTRICTION(KT, IU - 1)) then ! SW 6/26/2018 Valid for all K
                        if (BHR1(KT, IU - 1) > BCONSTRICTION(IU - 1)*H1(KT, IU - 1)) then
                            BHR1(KT, IU - 1) = BCONSTRICTION(IU - 1)*H1(KT, IU - 1)
                        end if
                    end if
                    if (UH_EXTERNAL(JB)) then
                        KB(IU - 1) = KB(IU)
                    end if
                    if (UH_INTERNAL(JB)) then
                        if (JBUH(JB) >= BS(JW) .and. JBUH(JB) <= BE(JW)) then
                            KB(IU - 1) = MIN(KB(UHS(JB)), KB(IU))
                        else
                            do KKB = KT, KMX
                                if (EL(KKB, IU) <= EL(KB(UHS(JB)), UHS(JB))) then
                                    exit
                                end if
                            end do
                            KB(IU - 1) = MIN(KKB, KB(IU))
                        end if
                    end if
                end if
                if (CONSTITUENTS) then ! SW 5/15/06
                    call TEMPERATURE_RATES()
                    call KINETIC_RATES()
                end if


!******** Total active cells

                do I = IU, ID
                    NTAC = NTAC - 1
                end do
                NTACMN = MIN(NTAC, NTACMN)
            end do
            call INTERPOLATION_MULTIPLIERS()

!****** Additional layer subtractions

            ZMIN(JW) = -1000.0
            do JB = BS(JW), BE(JW)
                if (BR_INACTIVE(JB)) then
                    cycle
                end if ! SW 6/12/2017
                do I = CUS(JB), DS(JB)
                    ZMIN(JW) = MAX(ZMIN(JW), Z(I))
                end do
            end do
            SUB_LAYER = ZMIN(JW) > 0.60*H(KT, JW) .and. KT < KTMAX ! SR 10/17/05
        end do
    end do

!** Temporary downstream head segment

    do JB = 1, NBR
        if (BR_INACTIVE(JB)) then
            cycle
        end if ! SW 6/12/2017
        if (DHS(JB) > 0) then
            do JJB = 1, NBR
                if (DHS(JB) >= US(JJB) .and. DHS(JB) <= DS(JJB)) then
                    exit
                end if
            end do
            if (CUS(JJB) > DHS(JB)) then
                CDHS(JB) = CUS(JJB)
            end if
        end if
    end do
    return
end subroutine LAYERADDSUB
