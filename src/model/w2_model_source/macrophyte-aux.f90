subroutine POROSITY()

    use GEOMC;     use GLOBAL;     use MACROPHYTEC;     use POROSITYC;     use LOGICC
    use SCREENC
    implicit none

    integer :: K, M, KUP, IEXIT, K1, KDN, IUT
    real :: B11, VSTOT, ELR, ELL, ELL1, ELL2, ELR2, EL1, EL2, ERR1, ERR2

    if (NIT == 0) then

        1040 format((8X,I8,3F8.0))

        do JW = 1, NWB
            KT = KTWB(JW)
            do JB = BS(JW), BE(JW)
                IU = CUS(JB)
                ID = DS(JB)
                do I = IU, ID
                    do K = 2, KB(I)
                        VOLI(K, I) = BH(K, I)*DLX(I)
                    end do
                    VOLI(KT, I) = BH2(KT, I)*DLX(I)
                end do
            end do
        end do

    end if

    do JB = 1, NBR
        COSA(JB) = COS(ALPHA(JB))
    end do

!C  CALCULATING # OF MACROPHYTE STEMS IN EACH CELL

    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            do I = IU, ID
!          HKTI  = H(KT,JW)-Z(I)    REPLACED BY H1(KT,I)
                if (KT == KTI(I)) then
                    VOLKTI(I) = H1(KT, I)*BIC(KT, I)*DLX(I)
                else
                    VOLKTI(I) = BIC(KTI(I), I)*(EL(KT, I) - EL(KTI(I) + 1, I) - Z(I)*COSA(JB))/COSA(JB)*DLX(I)
                end if
                do K = KTI(I) + 1, KT
                    VOLKTI(I) = VOLKTI(I) + VOLI(K, I)
                end do

                do M = 1, NMC
                    VSTEMKT(I, M) = MAC(KT, I, M)*VOLKTI(I)/DWV(M) !CB 6/29/06
                end do

                do K = KT + 1, KB(I)
                    do M = 1, NMC
                        VSTEM(K, I, M) = MAC(K, I, M)*VOLI(K, I)/DWV(M) !CB 6/29/06
                    end do
                end do
            end do
        end do
    end do

    POR = 1.0
    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            IU = CUS(JB)
            ID = DS(JB)

            do I = IU, ID
                do K = KT, KB(I)
                    if (K == KT) then
                        VSTOT = 0.0
                        do M = 1, NMC
                            VSTOT = VSTOT + VSTEMKT(I, M)
                        end do
                        POR(KT, I) = (VOLKTI(I) - VSTOT)/VOLKTI(I)
                    else
                        VSTOT = 0.0
                        do M = 1, NMC
                            VSTOT = VSTOT + VSTEM(K, I, M)
                        end do
                        POR(K, I) = (VOLI(K, I) - VSTOT)/VOLI(K, I)
                    end if
                end do
            end do

            do I = IU, ID
                do K = KTI(I), KB(I)
                    if (K <= KT) then
                        B(K, I) = POR(KT, I)*BIC(K, I)
                    else
                        B(K, I) = POR(K, I)*BIC(K, I)
                    end if

                end do
            end do

        end do
    end do



! BOUNDARY WIDTHS

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

!**** AREAS AND BOTTOM WIDTHS

            if (.not. TRAPEZOIDAL(JW)) then !SW 07/16/04
                do I = IU - 1, ID + 1
                    do K = 1, KMX - 1
                        BH2(K, I) = B(K, I)*H(K, JW)
                        BH(K, I) = B(K, I)*H(K, JW)
                        BB(K, I) = B(K, I) - (B(K, I) - B(K + 1, I))/(0.5*(H(K, JW) + H(K + 1, JW)))*H(K, JW)*0.5 !SW 08/02/04
                    end do
                    BH(KMX, I) = BH(KMX - 1, I)
                end do
!****** DERIVED GEOMETRY

                do I = IU - 1, ID + 1
                    BH2(KT, I) = B(KTI(I), I)*(EL(KT, I) - EL(KTI(I) + 1, I) - Z(I)*COSA(JB))/COSA(JB)
                    if (KT == KTI(I)) then
                        BH2(KT, I) = H2(KT, I)*B(KT, I)
                    end if
                    do K = KTI(I) + 1, KT
                        BH2(KT, I) = BH2(KT, I) + BH(K, I)
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

    10 continue

    return
end subroutine POROSITY

!************************************************************************
!**               S U B R O U T I N E    MACROPHYTE_FRICTION           **
!************************************************************************


subroutine MACROPHYTE_FRICTION(HRAD, BEDFR, EFFRIC, K, II)

    use GEOMC;     use GLOBAL;     use MACROPHYTEC;     use POROSITYC
    implicit none
    integer :: K, M, II
    real(R8) :: BEDFR, EFFRIC, HRAD
    real :: SAVOLRAT, XSAREA, TSAREA, ARTOT, SCTOT, CDAVG, FRIN

    do M = 1, NMC
        SAVOLRAT = DWV(M)/DWSA(M) !CB 6/29/2006
        if (K == KT) then
!      SAREA(M)=VSTEMKT(II,M)*SAVOLRAT/PI
            SAREA(M) = VSTEMKT(II, M)*SAVOLRAT*ANORM(M) !CB 6/29/2006
        else
!      SAREA(M)=VSTEM(K,II,M)*SAVOLRAT/PI
            SAREA(M) = VSTEM(K, II, M)*SAVOLRAT*ANORM(M) !CB 6/29/2006
        end if
    end do
    XSAREA = BH2(K, II)

    TSAREA = 0.0
    ARTOT = 0.0
    SCTOT = 0.0
    do M = 1, NMC
        ARTOT = ARTOT + SAREA(M)
        SCTOT = SCTOT + CDDRAG(M)*SAREA(M)
        TSAREA = TSAREA + SAREA(M)
    end do

    if (ARTOT > 0.0) then
        CDAVG = SCTOT/ARTOT
        FRIN = CDAVG*TSAREA*HRAD**(4./3.)/(2.0*G*XSAREA*DLX(II)*BEDFR**2)
        EFFRIC = BEDFR*SQRT(1.0 + FRIN)
    else
        EFFRIC = BEDFR
    end if

    return
end subroutine MACROPHYTE_FRICTION
