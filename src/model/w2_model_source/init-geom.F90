subroutine INITGEOM()
    use MAIN
    use GLOBAL;     use NAMESC;     use GEOMC;     use LOGICC;     use PREC;     use SURFHE;     use KINETIC;     use SHADEC;     use EDDY
    use STRUCTURES;     use TRANS;     use TVDC;     use SELWC;     use GDAYC;     use SCREENC;     use POROSITYC;     use MACROPHYTEC
    use RSTART
    implicit none

    integer :: NNBP, NCBP, NINTERNAL, NUP, KTMAX, JJW, IEXIT, KUP, KDN, K1, ICON, ISEG
    real :: ELL1, ELL2, EL1, EL2, B11, ERR1, ERR2, ELR, ELL, ELR2
    character(len=2) :: ICHAR2

    NPOINT = 0
    do JW = 1, NWB
        if (ZERO_SLOPE(JW)) then
            do I = US(BS(JW)) - 1, DS(BE(JW)) + 1
                EL(KMX, I) = ELBOT(JW)
                do K = KMX - 1, 1, -1
                    EL(K, I) = EL(K + 1, I) + H(K, JW)
                end do
            end do
        else
            EL(KMX, DS(JBDN(JW)) + 1) = ELBOT(JW)
            JB = JBDN(JW)
            NPOINT(JB) = 1
            NNBP = 1
            NCBP = 0
            NINTERNAL = 0
            NUP = 0
            do while (NNBP <= BE(JW) - BS(JW) + 1)
                NCBP = NCBP + 1
                if (NINTERNAL == 0) then
                    if (NUP == 0) then
                        do I = DS(JB), US(JB), -1
                            if (I /= DS(JB)) then
                                EL(KMX, I) = EL(KMX, I + 1) + SINA(JB)*(DLX(I) + DLX(I + 1))*0.5
                            else
                                EL(KMX, I) = EL(KMX, I + 1)
                            end if
                            do K = KMX - 1, 1, -1
                                EL(K, I) = EL(K + 1, I) + H(K, JW)*COSA(JB)
                            end do
                        end do
                    else
                        do I = US(JB), DS(JB)
                            if (I /= US(JB)) then
                                EL(KMX, I) = EL(KMX, I - 1) - SINA(JB)*(DLX(I) + DLX(I - 1))*0.5
                            else
                                EL(KMX, I) = EL(KMX, I - 1)
                            end if
                            do K = KMX - 1, 1, -1
                                EL(K, I) = EL(K + 1, I) + H(K, JW)*COSA(JB)
                            end do
                        end do
                        NUP = 0
                    end if
                    do K = KMX, 1, -1
                        if (UP_HEAD(JB)) then
                            EL(K, US(JB) - 1) = EL(K, US(JB)) + SINA(JB)*DLX(US(JB))
                        else
                            EL(K, US(JB) - 1) = EL(K, US(JB))
                        end if
                        if (DN_HEAD(JB)) then
                            EL(K, DS(JB) + 1) = EL(K, DS(JB)) - SINA(JB)*DLX(DS(JB))
                        else
                            EL(K, DS(JB) + 1) = EL(K, DS(JB))
                        end if
                    end do
                else
                    do K = KMX - 1, 1, -1
                        EL(K, UHS(JJB)) = EL(K + 1, UHS(JJB)) + H(K, JW)*COSA(JB)
                    end do
                    do I = UHS(JJB) + 1, DS(JB)
                        EL(KMX, I) = EL(KMX, I - 1) - SINA(JB)*(DLX(I) + DLX(I - 1))*0.5
                        do K = KMX - 1, 1, -1
                            EL(K, I) = EL(K + 1, I) + H(K, JW)*COSA(JB)
                        end do
                    end do
                    do I = UHS(JJB) - 1, US(JB), -1
                        EL(KMX, I) = EL(KMX, I + 1) + SINA(JB)*(DLX(I) + DLX(I + 1))*0.5
                        do K = KMX - 1, 1, -1
                            EL(K, I) = EL(K + 1, I) + H(K, JW)*COSA(JB)
                        end do
                    end do
                    NINTERNAL = 0
                end if
                if (NNBP == BE(JW) - BS(JW) + 1) then
                    exit
                end if

!****** Find next branch connected to furthest downstream branch

                do JB = BS(JW), BE(JW)
                    if (NPOINT(JB) /= 1) then
                        do JJB = BS(JW), BE(JW)
                            if (DHS(JB) >= US(JJB) .and. DHS(JB) <= DS(JJB) .and. NPOINT(JJB) == 1) then
                                NPOINT(JB) = 1
                                EL(KMX, DS(JB) + 1) = EL(KMX, DHS(JB)) + SINA(JB)*(DLX(DS(JB)) + DLX(DHS(JB)))*0.5
                                NNBP = NNBP + 1;                                 exit
                            end if
                            if (UHS(JJB) == DS(JB) .and. NPOINT(JJB) == 1) then
                                NPOINT(JB) = 1
                                EL(KMX, DS(JB) + 1) = EL(KMX, US(JJB)) + (SINA(JJB)*DLX(US(JJB)) + SINA(JB)*DLX(DS(JB)))*0.5
                                NNBP = NNBP + 1;                                 exit
                            end if
                            if (UHS(JJB) >= US(JB) .and. UHS(JJB) <= DS(JB) .and. NPOINT(JJB) == 1) then
                                NPOINT(JB) = 1
                                EL(KMX, UHS(JJB)) = EL(KMX, US(JJB)) + SINA(JJB)*DLX(US(JJB))*0.5
                                NNBP = NNBP + 1
                                NINTERNAL = 1;                                 exit
                            end if
                            if (UHS(JB) >= US(JJB) .and. UHS(JB) <= DS(JJB) .and. NPOINT(JJB) == 1) then
                                NPOINT(JB) = 1
                                EL(KMX, US(JB) - 1) = EL(KMX, UHS(JB)) - SINA(JB)*DLX(US(JB))*0.5
                                NNBP = NNBP + 1
                                NUP = 1;                                 exit
                            end if
                        end do
                        if (NPOINT(JB) == 1) then
                            exit
                        end if
                    end if
                end do
            end do
        end if
    end do

! Minimum/maximum layer heights

    do JW = 1, NWB
        do K = KMX - 1, 1, -1
            HMIN = DMIN1(H(K, JW), HMIN)
            HMAX = DMAX1(H(K, JW), HMAX)
        end do
    end do
    HMAX2 = HMAX**2

! Water surface and bottom layers

    do JW = 1, NWB
        do JB = BS(JW), BE(JW)
            do I = US(JB) - 1, DS(JB) + 1
                if (.not. RESTART_IN) then
                    KTI(I) = 2
                    do while (EL(KTI(I), I) > ELWS(I))
                        KTI(I) = KTI(I) + 1
                        if (kti(i) == kmx) then ! cb 7/7/2010 if elws below grid, setting to elws to just within grid ! SW 5/27/17 DELETED
                            kti(i) = kmx - 1 ! 2
                            elws(i) = el(kti(i), i)
                            exit
                        end if
                    end do

                    Z(I) = (EL(KTI(I), I) - ELWS(I))/COSA(JB)


                    ZMIN(JW) = DMAX1(ZMIN(JW), Z(I))
!KTI(I)   =  MAX(KTI(I)-1,2)   ! MOVED SW 5/27/17 
                    KTMAX = MAX(2, KTI(I))
                    KTWB(JW) = MAX(KTMAX, KTWB(JW))
                    KTI(I) = MAX(KTI(I) - 1, 2) ! original    IF(KTI(I) /= KMX)
                    if (Z(I) > ZMIN(JW)) then
                        IZMIN(JW) = I
                    end if
                end if
                K = 2
                do while (B(K, I) > 0.0)
                    KB(I) = K
                    K = K + 1
                end do
                KBMAX(JW) = MAX(KBMAX(JW), KB(I))
            end do
            KB(US(JB) - 1) = KB(US(JB))
            KB(DS(JB) + 1) = KB(DS(JB))
        end do


!** Correct for water surface going over several layers

        if (.not. RESTART_IN) then
            KT = KTWB(JW)
            do JB = BS(JW), BE(JW)
                do I = US(JB) - 1, DS(JB) + 1
                    H2(KT, I) = H(KT, JW) - Z(I)
                    K = KTI(I) + 1
                    do while (KT > K)
                        Z(I) = Z(I) - H(K, JW)
                        H2(KT, I) = H(KT, JW) - Z(I)
                        K = K + 1
                    end do
                end do
            end do
        end if
        ELKT(JW) = EL(KTWB(JW), DS(BE(JW))) - Z(DS(BE(JW)))*COSA(BE(JW))
    end do
    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            IU = US(JB)
            ID = DS(JB)

!**** Boundary bottom layers

            if (UH_EXTERNAL(JB)) then
                KB(IU - 1) = KB(IU)
            end if ! CB 6/12/07
            if (DH_EXTERNAL(JB)) then
                KB(ID + 1) = KB(ID)
            end if

!**** Branch numbers corresponding to tributaries, withdrawals, and head

            if (TRIBUTARIES) then
                do JT = 1, NTR
                    if (ITR(JT) >= US(JB) .and. ITR(JT) <= DS(JB)) then
                        JBTR(JT) = JB
                    end if
                end do
            end if
            if (WITHDRAWALS) then
                do JWD = 1, NWD
                    if (IWD(JWD) >= US(JB) .and. IWD(JWD) <= DS(JB)) then
                        JBWD(JWD) = JB
                    end if
                end do
            end if
            if (UH_INTERNAL(JB)) then
                do JJB = 1, NBR
                    JBUH(JB) = JJB
                    if (UHS(JB) >= US(JJB) .and. UHS(JB) <= DS(JJB)) then
                        exit
                    end if
                end do
                do JJW = 1, NWB
                    JWUH(JB) = JJW
                    if (JBUH(JB) >= BS(JJW) .and. JBUH(JB) <= BE(JJW)) then
                        exit
                    end if
                end do
            end if
            if (INTERNAL_FLOW(JB)) then
                do JJB = 1, NBR
                    JBUH(JB) = JJB
                    if (UHS(JB) >= US(JJB) .and. UHS(JB) <= DS(JJB)) then
                        exit
                    end if
                end do
                do JJW = 1, NWB
                    JWUH(JB) = JJW
                    if (JBUH(JB) >= BS(JJW) .and. JBUH(JB) <= BE(JJW)) then
                        exit
                    end if
                end do
            end if
            if (DH_INTERNAL(JB)) then
                do JJB = 1, NBR
                    JBDH(JB) = JJB
                    if (DHS(JB) >= US(JJB) .and. DHS(JB) <= DS(JJB)) then
                        exit
                    end if
                end do
                do JJW = 1, NWB
                    JWDH(JB) = JJW
                    if (JBDH(JB) >= BS(JJW) .and. JBDH(JB) <= BE(JJW)) then
                        exit
                    end if
                end do
            end if

!**** Bottom boundary cells

            if (UH_INTERNAL(JB)) then
                if (JBUH(JB) >= BS(JW) .and. JBUH(JB) <= BE(JW)) then
                    KB(IU - 1) = MIN(KB(UHS(JB)), KB(IU)) ! CB 6/12/07
                else
                    if (EL(KB(IU), IU) >= EL(KB(UHS(JB)), UHS(JB))) then ! CB 6/12/07
                        KB(IU - 1) = KB(IU) ! CB 6/12/07
                    else
                        do K = KT, KB(IU) ! CB 6/12/07
                            if (EL(KB(UHS(JB)), UHS(JB)) >= EL(K, IU)) then ! CB 6/12/07
                                KB(IU - 1) = K;                                 exit ! CB 6/12/07
                            end if
                        end do
                    end if
                end if
            end if
            if (DH_INTERNAL(JB)) then
                if (JBDH(JB) >= BS(JW) .and. JBDH(JB) <= BE(JW)) then
                    KB(ID + 1) = MIN(KB(DHS(JB)), KB(ID))
                else
                    if (EL(KB(ID), ID) >= EL(KB(DHS(JB)), DHS(JB))) then
                        KB(ID + 1) = KB(ID)
                    else
                        do K = KT, KB(ID)
                            if (EL(KB(DHS(JB)), DHS(JB)) >= EL(K, ID)) then
                                KB(ID + 1) = K;                                 exit
                            end if
                        end do
                    end if
                end if
            end if

!**** Boundary segment lengths

            DLX(IU - 1) = DLX(IU)
            DLX(ID + 1) = DLX(ID)

!**** Minimum bottom layers and average segment lengths

            do I = IU - 1, ID
                KBMIN(I) = MIN(KB(I), KB(I + 1))
                DLXR(I) = (DLX(I) + DLX(I + 1))*0.5
            end do
            KBMIN(ID + 1) = KBMIN(ID)
            DLXR(ID + 1) = DLX(ID)

!**** Minimum/maximum segment lengths

            do I = IU, ID
                DLXMIN = DMIN1(DLXMIN, DLX(I))
                DLXMAX = DMAX1(DLXMAX, DLX(I))
            end do
        end do
    end do


! Constrictions                       ! SW 6/26/2018
    CONSTRICTION = .false. ! INIITALIZE VARIABLES
    BCONSTRICTION = 0.0
    inquire(FILE="w2_constriction.csv", EXIST=Constriction(1, 1)) !  Just use the first element for file existance
    if (CONSTRICTION(1, 1)) then
        open(CON, FILE="w2_constriction.csv", STATUS="OLD")
        read(CON, *)
        read(CON, *) ICHAR2, ICON ! NUMBER OF CONSTRICTIONS
        if (ICHAR2 /= "ON") then
            CONSTRICTION = .false.
        else
            read(CON, *)
            do J = 1, ICON
                read(CON, *) ISEG, BCONSTRICTION(ISEG)
                CONSTRICTION(:, ISEG) = .true.
            end do
        end if
    end if


! Boundary widths

    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            IU = US(JB)
            ID = DS(JB)
            do I = IU - 1, ID + 1
                B(1, I) = B(2, I)
                do K = KB(I) + 1, KMX
                    B(K, I) = B(KB(I), I)
                end do
            end do
        end do
    end do

    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            IU = US(JB)
            ID = DS(JB)
            IEXIT = 0
            do K = 1, KMX - 1
                B(K, IU - 1) = B(K, IU)
                if (UH_INTERNAL(JB) .or. HEAD_FLOW(JB)) then
                    if (JBUH(JB) >= BS(JW) .and. JBUH(JB) <= BE(JW)) then
                        B(K, IU - 1) = B(K, UHS(JB))
                    else
                        ELR = EL(K, IU) + SINA(JB)*DLX(IU)*0.5
                        ELL = EL(2, UHS(JB)) - SINA(JBUH(JB))*DLX(UHS(JB))*0.5
                        if (ELR >= ELL) then
                            B(K, IU - 1) = B(2, UHS(JB))
                        else
                            do KUP = 2, KMX - 1
                                ELL1 = EL(KUP, UHS(JB)) - SINA(JBUH(JB))*DLX(UHS(JB))*0.5
                                ELL2 = EL(KUP + 1, UHS(JB)) - SINA(JBUH(JB))*DLX(UHS(JB))*0.5
                                if (ELL1 > ELR .and. ELL2 <= ELR) then
                                    if (KUP > KB(UHS(JB))) then
                                        KB(IU - 1) = K - 1
                                        KBMIN(IU - 1) = MIN(KB(IU), KB(IU - 1))
                                        IEXIT = 1
                                        exit
                                    end if
                                    ELR2 = EL(K + 1, IU) + SINA(JB)*DLX(IU)*0.5
                                    if (ELR2 >= ELL2) then
                                        B(K, IU - 1) = B(KUP, UHS(JB));                                         exit
                                    else
                                        K1 = KUP + 1
                                        if (K1 > KMX) then
                                            exit
                                        end if
                                        B11 = 0.0
                                        EL1 = ELR
                                        EL2 = EL(K1, UHS(JB)) - SINA(JBUH(JB))*DLX(IU)*0.5
                                        do while (ELR2 <= EL2)
                                            B11 = B11 + (EL1 - EL2)*B(K1 - 1, UHS(JB))
                                            EL1 = EL2
                                            K1 = K1 + 1
                                            if (K1 >= KMX + 1 .or. EL2 == ELR2) then
                                                exit
                                            end if
                                            EL2 = EL(K1, UHS(JB)) - SINA(JBUH(JB))*DLX(UHS(JB))*0.5
                                            if (EL2 <= ELR2) then
                                                EL2 = ELR2
                                            end if
                                        end do
                                        B(K, IU - 1) = B11/H(K, JW);                                         exit
                                    end if
                                end if
                            end do
                            if (EL(KMX, UHS(JB)) > EL(K, IU)) then
                                B(K, IU - 1) = B(K - 1, IU - 1)
                            end if
                            if (B(K, IU - 1) == 0.0) then
                                B(K, IU - 1) = B(K - 1, IU - 1)
                            end if
                            if (IEXIT == 1) then
                                exit
                            end if
                        end if
                    end if
                end if
            end do
            IEXIT = 0
            do K = 1, KMX - 1
                B(K, ID + 1) = B(K, ID)
                if (DH_INTERNAL(JB)) then
                    if (JBDH(JB) >= BS(JW) .and. JBDH(JB) <= BE(JW)) then
                        B(K, ID + 1) = B(K, DHS(JB))
                    else
                        ELL = EL(K, ID) - SINA(JB)*DLX(ID)*0.5
                        ELR = EL(2, DHS(JB)) + SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                        if (ELL >= ELR) then
                            B(K, ID + 1) = B(2, DHS(JB))
                        else
                            do KDN = 2, KMX - 1
                                ERR1 = EL(KDN, DHS(JB)) + SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                                ERR2 = EL(KDN + 1, DHS(JB)) + SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                                if (ERR1 >= ELL .and. ERR2 < ELL) then
                                    if (KDN > KB(DHS(JB))) then
                                        KB(ID + 1) = K - 1
                                        KBMIN(ID) = MIN(KB(ID), KB(ID + 1))
                                        IEXIT = 1
                                        exit
                                    end if
                                    ELL2 = EL(K + 1, ID) - SINA(JB)*DLX(ID)*0.5
                                    if (ELL2 >= ERR2) then
                                        B(K, ID + 1) = B(KDN, DHS(JB));                                         exit
                                    else
                                        K1 = KDN + 1
                                        if (K1 > KMX) then
                                            exit
                                        end if
                                        B11 = 0.0
                                        EL2 = ELL
                                        EL1 = EL(K1, DHS(JB)) + SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                                        do while (ELL2 <= EL1)
                                            B11 = B11 + (EL2 - EL1)*B(K1 - 1, DHS(JB))
                                            EL2 = EL1
                                            K1 = K1 + 1
                                            if (K1 >= KMX + 1 .or. EL1 == ELL2) then
                                                exit
                                            end if
                                            EL1 = EL(K1, DHS(JB)) + SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                                            if (EL1 <= ELL2) then
                                                EL1 = ELL2
                                            end if
                                        end do
                                        B(K, ID + 1) = B11/H(K, JW);                                         exit
                                    end if
                                end if
                            end do
                            if (EL(KMX, DHS(JB)) > EL(K, ID)) then
                                B(K, ID + 1) = B(K - 1, ID + 1)
                            end if
                            if (B(K, ID + 1) == 0.0) then
                                B(K, ID + 1) = B(K - 1, ID + 1)
                            end if
                            if (IEXIT == 1) then
                                exit
                            end if
                        end if
                    end if
                end if
            end do
        end do ! SW 1/23/06
    end do ! SW 1/23/06
    BNEW = B ! SW 1/23/06
    KBI = KB ! SW 10/29/2010

!**** Upstream active segment and single layer  ! 1/23/06 entire section moved SW
    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            IU = US(JB)
            ID = DS(JB)
            IEXIT = 0
            if (SLOPE(JB) /= 0.0) then
                do I = US(JB) - 1, DS(JB) + 1
                    if (KBi(I) < KT) then ! SW 1/23/06
                        do k = kbi(i) + 1, kt
                            Bnew(K, I) = 0.000001 ! SW 1/23/06
                        end do
                        KB(I) = KT
                    end if
                end do
            end if
            IUT = IU
            do I = IU, ID
                if (KB(I) - KT < NL(JB) - 1) then
                    IUT = I + 1
                end if
                ONE_LAYER(I) = KT == KB(I)
            end do
            do I = IU - 1, ID ! recompute kbmin after one_layer computation   11/12/07
                KBMIN(I) = MIN(KB(I), KB(I + 1))
            end do
            KBMIN(ID + 1) = KBMIN(ID)

            CUS(JB) = IUT
            if (IUT >= DS(JB)) then
                BR_INACTIVE(JB) = .true.
            end if ! SW 6/12/2017


!**** Areas and bottom widths

            if (.not. TRAPEZOIDAL(JW)) then !SW 07/16/04
                do I = IU - 1, ID + 1
                    do K = 1, KMX - 1
                        BH2(K, I) = B(K, I)*H(K, JW)
                        BH(K, I) = B(K, I)*H(K, JW)
                        BB(K, I) = B(K, I) - (B(K, I) - B(K + 1, I))/(0.5*(H(K, JW) + H(K + 1, JW)))*H(K, JW)*0.5 !SW 08/02/04
                    end do
                    BH(KMX, I) = BH(KMX - 1, I)
                end do

! column widths
                do I = IU - 1, ID + 1
                    CW(KB(I), I) = B(KB(I), I)
                    do K = 1, KB(I) - 1
                        CW(K, I) = B(K, I) - B(K + 1, I)
                    end do
                end do

!****** Derived geometry

                do I = IU - 1, ID + 1
                    BH2(KT, I) = B(KTI(I), I)*(EL(KT, I) - EL(KTI(I) + 1, I) - Z(I)*COSA(JB))/COSA(JB)
                    if (KT == KTI(I)) then
                        BH2(KT, I) = H2(KT, I)*B(KT, I)
                    end if
                    do K = KTI(I) + 1, KT
                        BH2(KT, I) = BH2(KT, I) + Bnew(K, I)*H(K, JW) ! sw 1/23/06    BH(K,I)
                    end do
                    BKT(I) = BH2(KT, I)/H2(KT, I)
                    BI(KT, I) = B(KTI(I), I)
                end do
            else !SW 07/16/04
                do I = IU - 1, ID + 1
                    do K = 1, KMX - 1
                        BB(K, I) = B(K, I) - (B(K, I) - B(K + 1, I))/(0.5*(H(K, JW) + H(K + 1, JW)))*H(K, JW)*0.5
                    end do
                    BB(KB(I), I) = B(KB(I), I)*0.5
                    BH2(1, I) = B(1, I)*H(1, JW)
                    BH(1, I) = BH2(1, I)
                    do K = 2, KMX - 1
                        BH2(K, I) = 0.25*H(K, JW)*(BB(K - 1, I) + 2.*B(K, I) + BB(K, I))
                        BH(K, I) = BH2(K, I)
                    end do
                    BH(KMX, I) = BH(KMX - 1, I)
                end do
                do I = IU - 1, ID + 1
                    call GRID_AREA1(EL(KT, I) - Z(I), EL(KT + 1, I), BH2(KT, I), BI(KT, I))
                    BKT(I) = BH2(KT, I)/H2(KT, I)
                end do
            end if
            do I = IU - 1, ID
                do K = 1, KMX - 1
                    AVH2(K, I) = (H2(K, I) + H2(K + 1, I))*0.5
                    AVHR(K, I) = H2(K, I) + (H2(K, I + 1) - H2(K, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04
                end do
                AVH2(KMX, I) = H2(KMX, I)
                do K = 1, KMX
                    BR(K, I) = B(K, I) + (B(K, I + 1) - B(K, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04
                    BHR(K, I) = BH(K, I) + (BH(K, I + 1) - BH(K, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04
                    BHR2(K, I) = BH2(K, I) + (BH2(K, I + 1) - BH2(K, I))/(0.5*(DLX(I) + DLX(I + 1)))*0.5*DLX(I) !SW 07/29/04
                    if (CONSTRICTION(K, I)) then ! SW 6/26/2018
                        if (BR(K, I) > BCONSTRICTION(I)) then
                            BR(K, I) = BCONSTRICTION(I)
                        end if
                        if (BHR(K, I) > BCONSTRICTION(I)*H(K, JW)) then
                            BHR(K, I) = BCONSTRICTION(I)*H(K, JW)
                        end if
                        if (BHR2(K, I) > BCONSTRICTION(I)*H(K, JW)) then
                            BHR2(K, I) = BCONSTRICTION(I)*H(K, JW)
                        end if
                    end if

                end do
            end do
            do K = 1, KMX - 1
                AVH2(K, ID + 1) = (H2(K, ID + 1) + H2(K + 1, ID + 1))*0.5
                BR(K, ID + 1) = B(K, ID + 1)
                BHR(K, ID + 1) = BH(K, ID + 1)
            end do
            AVH2(KMX, ID + 1) = H2(KMX, ID + 1)
            AVHR(KT, ID + 1) = H2(KT, ID + 1)
            BHR2(KT, ID + 1) = BH2(KT, ID + 1)
            IUT = IU
            if (UP_HEAD(JB)) then
                IUT = IU - 1
            end if
            do I = IUT, ID
                do K = 1, KMX - 1
                    VOL(K, I) = B(K, I)*H2(K, I)*DLX(I)
                end do
                VOL(KT, I) = BH2(KT, I)*DLX(I)
                DEPTHB(KT, I) = H2(KT, I)
                DEPTHM(KT, I) = H2(KT, I)*0.5
                do K = KT + 1, KMX
                    DEPTHB(K, I) = DEPTHB(K - 1, I) + H2(K, I)
                    DEPTHM(K, I) = DEPTHM(K - 1, I) + (H2(K - 1, I) + H2(K, I))*0.5
                end do
            end do
        end do
    end do
    H1 = H2
    BH1 = BH2
    BHR1 = BHR2
    AVH1 = AVH2

! Temporary downstream head segment

    do JB = 1, NBR
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

! Total active cells

    do JW = 1, NWB
        do JB = BS(JW), BE(JW)
            if (BR_INACTIVE(JB)) then
                cycle
            end if ! SW 6/12/2017
            do I = CUS(JB), DS(JB)
                do K = KTWB(JW), KB(I)
                    NTAC = NTAC + 1
                end do
            end do
            NTACMX = NTAC
            NTACMN = NTAC

!**** Wind fetch lengths

            do I = US(JB), DS(JB)
                FETCHD(I, JB) = FETCHD(I - 1, JB) + DLX(I)
            end do
            do I = DS(JB), US(JB), -1
                FETCHU(I, JB) = FETCHU(I + 1, JB) + DLX(I)
            end do
        end do
    end do

! Segment heights

    do JW = 1, NWB
        do JB = BS(JW), BE(JW)
            do I = US(JB) - 1, DS(JB) + 1
                do K = MIN(KMX - 1, KB(I)), 2, -1
                    HSEG(K, I) = HSEG(K + 1, I) + H2(K, I)
                end do
            end do
        end do
    end do

! Beginning and ending segments/layers for snapshots

    do JW = 1, NWB
        do I = 1, NISNP(JW)
            KBR(JW) = MAX(KB(ISNP(I, JW)), KBR(JW))
        end do
    end do

    return
end subroutine INITGEOM
