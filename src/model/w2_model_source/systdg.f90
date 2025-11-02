module modSYSTDG
    use PREC, only: R8
    use TDGAS;     use STRUCTURES;     use GLOBAL;     use MAIN, only: EA, GTTYP, GTPC, Q, WBSEG, TEXT, ERROR_OPEN, SYSTDGC, N2BNDC, DOBNDC, TDGTAC, TMSTRT, CONTDG, TDG2BNDC
    use TVDC, only: TAIR, TDEW;     use KINETIC, only: TDG
    implicit none
    real(R8), allocatable, dimension(:) :: TDG_PHS, TDG_FLS, TDG_TDP
    real(R8), allocatable, dimension(:) :: BAYC, QBAY
    real(R8) :: qs, TDG_ROSP, TDG_TDG, TDG_REL
    real(R8) :: TDGP1, TDGP2, TDGP3, TDGP4, TDGP12, TDGP22, TDGP32, TDGP42, TDGE1, TDGE2, TDGE12, TDGE22, ROP1, ROP2, ROP3, ROP4
    real(R8) :: TWE, TWCE, TWE_TS, FBE, QSPILL, TDGSPMN
    integer, allocatable, dimension(:) :: POWGTNO, FLGTNO
    integer :: POWNO, FLNO, BEGNO, ENDNO, NRO, NBAY
    integer :: TWEMOD, TDGEQ, TDGROEQ, TDGENTEQ
    integer :: ig, ip, ifl, ib
    logical, allocatable, dimension(:) :: GTNAME
    character(len=72) :: TITLESYSTDG(10)
    character(len=8) :: TWETSC, TDGLOC
    character(len=72) :: TWEFN
    real :: NXTSPLIT3
!

contains

    subroutine INPUT_SYSTDG()
!USE MAIN, ONLY: TMSTRT, CONTDG       !, GTTYP, GTPC, SYSTDGC, N2BNDC, DOBNDC, TDGTAC
        implicit none
        character(len=8) :: AID1
        logical :: CSVFORMAT
        integer :: I, IG, N_POW, N_FLD, N_SPB
        NXTSPLIT3 = TMSTRT
        open(88888, FILE="TDG_output.csv", status="unknown")
        write(88888, '(A, <NGT>("QGT-",I2,","))') "JDAY,TDG_TDG,SUM_QGT2,", (IG, IG = 1, NGT)
        NRO = 0
        POWNO = 0
        FLNO = 0
        NBAY = 0

        CSVFORMAT = .false.
        read(CONTDG, "(A)") TITLESYSTDG(1)
        if (TITLESYSTDG(1)(1:1) == "$") then
            CSVFORMAT = .true.
        end if

        if (.not. CSVFORMAT) then
            read(CONTDG, "(A)") TITLESYSTDG(1) ! READ NEXT LINE - IF COMMAS IN FIRST FEW FIELDS IT IS IN CSV FORMAT
            do I = 1, 7
                if (TITLESYSTDG(1)(I:I) == ",") then
                    CSVFORMAT = .true.
                    exit
                end if
            end do
            rewind(CONTDG)
        end if

        if (CSVFORMAT) then
            read(CONTDG, *)
            read(CONTDG, *)
            read(CONTDG, *)
            do I = 1, 10
                read(CONTDG, *) AID1, TITLESYSTDG(I)
            end do
            read(CONTDG, *)
            read(CONTDG, *)
            read(CONTDG, *) AID1, SYSTDGC, N2BNDC, DOBNDC, TDG2BNDC, TDGTAC;             SYSTDGC = ADJUSTR(SYSTDGC);             N2BNDC = ADJUSTR(N2BNDC);             DOBNDC = ADJUSTR(DOBNDC);             TDGTAC = ADJUSTR(TDGTAC)
            TDG2BNDC = ADJUSTR(TDG2BNDC)
            read(CONTDG, *)
            read(CONTDG, *)
            do IG = 1, NGT
                read(CONTDG, *) AID1, GTTYP(IG), GTPC(IG)
            end do
            GTTYP = ADJUSTR(GTTYP)

            do IG = 1, NGT
                if (GTTYP(IG) == "     POW") then
                    POWNO = POWNO + 1
                end if
                if (GTTYP(IG) == "     FLD") then
                    FLNO = FLNO + 1
                end if
                if (GTTYP(IG) == "      RO") then
                    NRO = NRO + 1
                end if
                if (GTTYP(IG) == "     SPB") then
                    NBAY = NBAY + 1
                end if
            end do
            read(CONTDG, *)
            read(CONTDG, *)

            read(CONTDG, *) AID1, FBE, TWCE, TWEMOD, TWE, TWETSC, TDGLOC, QSPILL, TDGSPMN;             TWETSC = ADJUSTR(TWETSC);             TDGLOC = ADJUSTR(TDGLOC)
            read(CONTDG, *)
            read(CONTDG, *)
            read(CONTDG, *) AID1, TDGEQ, TDGP1, TDGP2, TDGP3, TDGP4, TDGP12, TDGP22, TDGP32, TDGP42
            if (NRO > 0) then
                read(CONTDG, *) AID1, TDGROEQ, ROP1, ROP2, ROP3, ROP4
            end if
            read(CONTDG, *)
            read(CONTDG, *)

            read(CONTDG, *) AID1, TDGENTEQ, TDGE1, TDGE2, TDGE12, TDGE22
            read(CONTDG, *)
            read(CONTDG, *)

            read(CONTDG, *) AID1, TWEFN
            close(CONTDG)
        else
            read(CONTDG, "(///(8X,A72))") (TITLESYSTDG(i), i = 1, 10)
            read(CONTDG, "(//8x,5A8)") SYSTDGC, N2BNDC, DOBNDC, TDG2BNDC, TDGTAC
            read(CONTDG, "(//(:8X,A8,F8.2))") (GTTYP(IG), GTPC(IG), IG = 1, NGT)

            do IG = 1, NGT
                if (GTTYP(IG) == "     POW") then
                    POWNO = POWNO + 1
                end if
                if (GTTYP(IG) == "     FLD") then
                    FLNO = FLNO + 1
                end if
                if (GTTYP(IG) == "      RO") then
                    NRO = NRO + 1
                end if
                if (GTTYP(IG) == "     SPB") then
                    NBAY = NBAY + 1
                end if
            end do

            read(CONTDG, "(//8X,2F8.3,I8,F8.3,2A8,2F8.3)") FBE, TWCE, TWEMOD, TWE, TWETSC, TDGLOC, QSPILL, TDGSPMN
            read(CONTDG, "(//8X,I8,8F8.3)") TDGEQ, TDGP1, TDGP2, TDGP3, TDGP4, TDGP12, TDGP22, TDGP32, TDGP42
            if (NRO > 0) then
                read(CONTDG, "(8X,I8,4F8.5)") TDGROEQ, ROP1, ROP2, ROP3, ROP4
            end if
            read(CONTDG, "(//8X,I8,4F8.3)") TDGENTEQ, TDGE1, TDGE2, TDGE12, TDGE22
            read(CONTDG, "(//(8X,A72))") TWEFN
            close(CONTDG)
        end if
        if (SYSTDGC == "     OFF") then
            go to 100
        end if ! DO NOT ALLOCATE ARRAYS IF WE ARE NOT USING SYSTDG
        if (POWNO > 0) then
            allocate(POWGTNO(POWNO), TDG_PHS(POWNO))
        else
            allocate(POWGTNO(1), TDG_PHS(1))
            POWGTNO(1) = 0
        end if
        if (FLNO > 0) then
            allocate(FLGTNO(FLNO), TDG_FLS(FLNO))
        else
            allocate(FLGTNO(1), TDG_FLS(1))
            FLGTNO(1) = 0
        end if
        N_POW = 0
        N_FLD = 0
        do ig = 1, NGT
            if (GTTYP(ig) == "     POW") then
                N_POW = N_POW + 1
                POWGTNO(N_POW) = ig
            end if
            if (GTTYP(ig) == "     FLD") then
                N_FLD = N_FLD + 1
                FLGTNO(N_FLD) = ig
            end if
        end do
        do ig = 1, NGT
            if (GTTYP(ig) == "     SPB") then
                BEGNO = ig
                exit
            end if
        end do
        ENDNO = BEGNO + NBAY - 1
        if (NBAY /= 0) then
            allocate(BAYC(NBAY), QBAY(NBAY))
        else
            allocate(BAYC(1), QBAY(1))
        end if
        N_SPB = 0
        do ig = 1, NGT
            if (GTTYP(ig) == "     SPB") then
                N_SPB = N_SPB + 1
                BAYC(N_SPB) = GTPC(ig)
            end if
        end do
!
        allocate(GTNAME(NGT), TDG_TDP(NGT))
        GTNAME = .false.
        do ig = 1, NGT
            if (GTTYP(ig) == "      RO") then
                GTNAME(ig) = .true.
            end if ! C IN RO IS UPDATE FROM TDG
            if (ig >= BEGNO .and. ig <= ENDNO) then
                GTNAME(ig) = .true.
            end if ! GTNAME = GATE IS A SPILLBAY== .TRUE.
        end do
        100 return
    end subroutine INPUT_SYSTDG

!===========================================================================================================================
! allocate and initialize all input parameter

    subroutine SYSTDG_qs()
        implicit none
        real(R8) :: SUMQ, SUMQ1
        SUMQ = 0.0
        SUMQ1 = 0.0
        do ib = 1, NBAY
            SUMQ = SUMQ + QBAY(ib)**BAYC(ib)
            SUMQ1 = SUMQ1 + QBAY(ib)**(BAYC(ib) - 1.0)
        end do
        if (SUMQ1 /= 0.0) then
            qs = SUMQ/SUMQ1*35.3147/1000.0
        end if ! CMS TO KCFS 
    end subroutine SYSTDG_qs

!===========================================================================================================================
! TDG production calculation in SYSTDG

    subroutine UPDATE_TDGC(NSAT, P, N, T, TDGC)
        implicit none
        integer :: NSAT, N
        real(R8) :: P, T, TDGC
        real(R8) :: SAT
! CALCULATE SAT
        if (NSAT == 0) then ! O2 saturation
            SAT = EXP(7.7117 - 1.31403*LOG(T + 45.93))*P
        else
            if (NSAT == 1) then ! N2 saturation
                EA = DEXP(2.3026D0*(7.5D0*TDEW(WBSEG(IUGT(N)))/(TDEW(WBSEG(IUGT(N))) + 237.3D0) + 0.6609D0))*0.001316 ! mmHg     
                SAT = 1.5568D06*0.79*(P - EA)*(1.8816D-5 - 4.116D-7*T + 4.6D-9*T*T)
            else
                SAT = P
            end if
        end if
        TDGC = TDG_TDP(N)*SAT/100.0
        if (POWNO > 0 .and. TDGLOC == "     REL") then
            TDGC = TDG_TDG*SAT/100.0
        end if
    end subroutine UPDATE_TDGC


    subroutine SYSTDG_TDG()
        use SCREENC, only: JDAY; 
        implicit none
        real(R8) :: P1, P2, P3, P4, E1, E2
        real(R8) :: qs_RO, W2FBE, Q_SUM, TEMP_TW, SUM_TDGPHK
        real(R8) :: Q_ROSP, QRO, QSP, QPH, QTOT, TDG_QROSP, TDG_QPH, TDG_QTOT, TDG_QENT
        real(R8) :: SUM_TDG_ROS, SUM_TDG_SPS, SUM_TDG_PHS
        real(R8) :: TDG_RO, TDG_SP, TDG_PH
        integer(8) :: SUM_K, IK
        real :: SUM_QGT2
        Q_SUM = 0.0
! ADD QRO AND QSP
        do ig = 1, NGT
            if (GTTYP(ig) == "      RO") then
                Q_SUM = Q_SUM + QGT(ig)
            else
                if (ig >= BEGNO .and. ig <= ENDNO) then
                    Q_SUM = Q_SUM + QGT(ig)
                end if
            end if
        end do
        TDG_ROSP = 0.0
        if (Q_SUM /= 0.0) then
            TDG_TDP(:) = 0.0 ! INITIAL TDG FOR GATES
!
            SUM_TDG_ROS = 0.0 ! INITIAL SUM TDG FOR RO
            QRO = 0.0 ! INITIAL SUM Q FOR RO
            do ig = 1, NGT
                if (GTTYP(ig) == "      RO") then
                    QRO = QRO + QGT(ig)
                end if
            end do
            qs_RO = QRO*35.3147/1000.0 ! CMS TO KCFS 
            do ig = 1, NGT
                if (GTTYP(ig) == "      RO" .and. QGT(ig) > 0.0) then
                    TDG_TDP(ig) = (ROP1*(1 - EXP(ROP3*qs_RO)) + PALT(IUGT(ig))*760.0)/(PALT(IUGT(ig))*760.0)*100.0 ! RO USE EQ 1     ! (mmgh) TO TDG (%)
                    if (TDG_TDP(ig) > 145.0) then
                        TDG_TDP(ig) = 145.0
                    end if ! TDG <= 145.0
                    SUM_TDG_ROS = SUM_TDG_ROS + TDG_TDP(ig)*QGT(ig) ! SUM TDG FOR RO
                end if
            end do
            if (QRO /= 0.0) then
                TDG_RO = SUM_TDG_ROS/QRO ! FLOW AVERAGED TDG_RO (%)
            else
                TDG_RO = 0.0 ! ZERO FLOW --> RO TDG IS 0.0
            end if
! SET E1 E2 INCASE ONLY RO FLOW TO CALCULATE TDG_REL
            if (POWNO > 0 .and. TDGLOC == "     REL") then
                do ip = 1, POWNO
                    if (QGT(POWGTNO(ip)) /= 0.0) then
                        W2FBE = Q(IUGT(ip) - US(JBUGT(ip)) + 1) ! FLOW DISCHARGE BEFORE DAM
                        exit
                    end if
                end do
                if (FBE < 0.0) then
                    E1 = TDGE1
                    E2 = TDGE2
                else
                    if (W2FBE < FBE) then
                        E1 = TDGE1
                        E2 = TDGE2
                    else
                        if (W2FBE >= FBE) then
                            E1 = TDGE12
                            E2 = TDGE22
                        end if
                    end if
                end if
            end if
!
            if (NBAY /= 0) then
                QBAY(1:NBAY) = QGT(BEGNO:ENDNO)
            end if ! BAY FLOW QBAY FROM GATE FLOW QGT
            call SYSTDG_qs() ! CALCULATE qs                              
! TWE Recalculation
            if (TWETSC == "      ON") then
                TWE = TWE_TS
            end if ! UPDATE TWE TO TWE_TS ACCORDING TO CONTROL VARIABLE TWETSC 
            if (TWEMOD == 1) then
                TWE = TWE*0.934 + 4.94
            end if ! UPDATE TWE ACCORDING TO TWEMOD
            SUM_TDG_SPS = 0.0 ! INITIAL TDG*QSP FOR SPILL
            QSP = 0.0 ! INITIAL SUM Q FOR SPILL
            TDG_SP = 0.0 ! INITIAL TDG FOR SPILL
! CALCULATE TDG_TDP FOR EACH BAY FROM BEGNO TO ENDNO
            do ig = BEGNO, ENDNO
                if (QGT(ig) /= 0.0) then
! READ W2FBE
                    if (IUGT(ig) >= US(JBUGT(ig)) + 1) then ! GATE IS NOT LOCATED IN THE 1ST SEGMENT                 
                        W2FBE = Q(IUGT(ig) - US(JBUGT(ig)) + 1) ! FLOW DISCHARGE BEFORE DAM
                    else
                        W2FBE = 0.0
                    end if
                    if (FBE < 0.0) then
                        P1 = TDGP1
                        P2 = TDGP2
                        P3 = TDGP3
                        P4 = TDGP4
                        E1 = TDGE1
                        E2 = TDGE2
                    else
                        if (W2FBE < FBE) then
                            P1 = TDGP1
                            P2 = TDGP2
                            P3 = TDGP3
                            P4 = TDGP4
                            E1 = TDGE1
                            E2 = TDGE2
                        else
                            if (W2FBE >= FBE) then
                                P1 = TDGP12
                                P2 = TDGP22
                                P3 = TDGP32
                                P4 = TDGP42
                                E1 = TDGE12
                                E2 = TDGE22
                            end if
                        end if
                    end if
! TDG EQUATIONS 1/2/3/4/5
                    if (TDGEQ == 2) then
                        TDG_TDP(ig) = P1*(TWE - TWCE)**P2*(1.0 - EXP(P3*qs)) + P4 + PALT(IUGT(ig))*760.0 ! mmHg
!
                    else
                        if (TDGEQ == 3) then
                            TDG_TDP(ig) = P1*(TWE - TWCE)**P2*qs**P3 + P4 + PALT(IUGT(ig))*760.0 ! mmHg
!
                        else
                            if (TDGEQ == 4) then
                                TDG_TDP(ig) = P1*(TWE - TWCE) + P2*qs**P3 + P4 + PALT(IUGT(ig))*760.0 ! mmHg
!
                            else
                                if (TDGEQ == 5) then
                                    SUM_K = 0
                                    TEMP_TW = 0.0
                                    do IK = KT, KB(IUGT(ig) + 1)
                                        if (T2(IK, IUGT(ig) + 1) /= -99.0 .and. T2(IK, IUGT(ig) + 1) > 0.0) then ! WATER TEMP IS >0.0
                                            SUM_K = SUM_K + 1
                                            TEMP_TW = TEMP_TW + T2(IK, IUGT(ig) + 1) ! SUM EVERY LAYER OF TAIL WATER
                                        end if
                                    end do
                                    if (SUM_K /= 0) then
                                        TEMP_TW = TEMP_TW/SUM_K
                                    end if ! AVERAGED TEMP OF TAIL WATER
                                    TDG_TDP(ig) = P1*(1.0 - EXP(P2*qs)) + P3*(TEMP_TW - P4) + PALT(IUGT(ig))*760.0 ! mmHg
!
                                else
                                    if (TDGEQ == 1) then
                                        TDG_TDP(ig) = P1*(1.0 - EXP(P3*qs)) + PALT(IUGT(ig))*760.0 ! mmHg
!
                                    else
                                        TEXT = "TDGEQ INPUT ERROR"
                                        ERROR_OPEN = .true.
!
                                    end if
                                end if
                            end if
                        end if
                    end if
                    TDG_TDP(ig) = TDG_TDP(ig)/(PALT(IUGT(ig))*760.0)*100.0 ! (mmgh) TO TDG (%)
                    if (TDG_TDP(ig) > 145.0) then
                        TDG_TDP(ig) = 145.0
                    end if ! TDG <= 145.0
                    SUM_TDG_SPS = SUM_TDG_SPS + TDG_TDP(ig)*QGT(ig) ! SUM TDG_SP*QGT 
                    QSP = QSP + QGT(ig) ! SUM Q OF SPILL
! If total spill <= 50 kcfs, TDG % saturation = 110 % in the spillway outlets.  
                    if (QSP*35.3147/1000.0 > 0.0 .and. QSP*35.3147/1000.0 <= QSPILL) then
                        TDG_TDP(ig) = TDGSPMN
                    end if
                end if
            end do
            if (QSP /= 0.0) then
                TDG_SP = SUM_TDG_SPS/QSP
            end if ! FLOW AVERAGED TDG FOR SPILL
            Q_ROSP = QRO + QSP ! ADD QSP INTO Q_ROSP
            TDG_ROSP = (SUM_TDG_ROS + SUM_TDG_SPS)/Q_ROSP ! FLOW AVERAGED TDG FOR SPILL WITH RO
            TDG_TDG = TDG_ROSP ! OUTPUT TDG = TDG_SP
! If total spill <= 50 kcfs, TDG % saturation = 110 % in the spillway outlets.  
            if (QSP*35.3147/1000.0 > 0.0 .and. QSP*35.3147/1000.0 <= QSPILL) then
                TDG_TDG = TDGSPMN
            end if
!
! QENT CALCULATIONS UPDATE TDG_TDG TO TDG_REL
            if (POWNO > 0 .and. TDGLOC == "     REL") then
! TDG POWER HOUSE(S) AND Q POWER HOUSE(S)
                QPH = 0.0 ! INITIAL Q OF POWER HOUSE(S)
                TDG_PHS(:) = 0.0 ! INITIAL TDG OF POWER HOUSE(S)                               
                SUM_TDG_PHS = 0.0 ! INITIAL SUM TDG*QPH
                do ip = 1, POWNO
                    if (QGT(POWGTNO(ip)) /= 0.0) then
                        QPH = QPH + QGT(POWGTNO(ip)) ! Q POWER HOUSE(S)
                        SUM_K = 0 ! INITIAL SUM K
                        SUM_TDGPHK = 0.0 ! INITIAL SUM TDG OF POWER HOUSE(S) SEGMENT
                        do IK = 1, KMX
                            if (TDG(IK, IUGT(POWGTNO(ip))) > 0.0) then ! TDG > 0.0 LAYER 
                                SUM_K = SUM_K + 1 ! SUM LAYER COUNT
                                SUM_TDGPHK = SUM_TDGPHK + TDG(IK, IUGT(POWGTNO(ip))) ! SUM TDG OF POWER HOUSE(S)
                            end if
                        end do
                        if (SUM_K /= 0) then
                            TDG_PHS(ip) = SUM_TDGPHK/SUM_K ! LAYER AVERAGED TDG FOR POWER HOUSE(S)
                        else
                            TDG_PHS(ip) = 0.0
                        end if
                        SUM_TDG_PHS = SUM_TDG_PHS + TDG_PHS(ip)*QGT(POWGTNO(ip)) ! SUM TDG*QPH
                    end if
                end do
                if (QPH /= 0.0) then
                    TDG_PH = SUM_TDG_PHS/QPH
                end if ! FLOW AVERAGED TDG FOR POWER HOUSE(S)
                QTOT = 0.0 ! INITIAL QTOT
                do ig = 1, NGT
                    QTOT = QTOT + QGT(ig) ! SUM QTOT
                end do
                if (FLNO > 0) then
                    do ifl = 1, FLNO
                        QTOT = QTOT - QGT(FLGTNO(ifl)) ! QTOT = QTOT - Q FISH LADDER
                    end do
                end if
                TDG_QROSP = Q_ROSP*35.3147/1000.0 ! CMS TO KCFS
                TDG_QPH = QPH*35.3147/1000.0 ! CMS TO KCFS
                TDG_QTOT = QTOT*35.3147/1000.0 ! CMS TO KCFS
! CALCULATE QENT
                if (TDGENTEQ == 1) then
                    TDG_QENT = E1*TDG_QROSP + E2
!
                else
                    if (TDGENTEQ == 2) then
                        TDG_QENT = MIN(TDG_QTOT/60, 1.0)*E1*TDG_QROSP + E2
!
                    else
                        if (TDGENTEQ == 3) then
                            TDG_QENT = MIN(TDG_QROSP/20, 1.0)*E1*TDG_QROSP + E2
!
                        else
                            TEXT = "TDGENTEQ INPUT ERROR"
                            ERROR_OPEN = .true.
                        end if
                    end if
                end if
                TDG_QENT = MIN(TDG_QENT, TDG_QTOT - TDG_QROSP)
                TDG_REL = (TDG_ROSP*(TDG_QROSP + TDG_QENT) + TDG_PH*(TDG_QPH - TDG_QENT))/(TDG_QPH + TDG_QROSP) ! TDG RELEASE CALCULATION
                if (TDG_REL > 145.0) then
                    TDG_REL = 145.0
                end if ! TDG RELEASE <=145.0
                TDG_TDG = TDG_REL ! OUTPUT TDG = TDG_REL
            end if ! END IF POWNO >0 
            if (TDG_TDG > 145.0) then
                TDG_TDG = 145.0
            end if ! TDG <=145.0
            if (JDAY >= NXTSPLIT3) then
                do ig = 1, NGT
                    SUM_QGT2 = SUM_QGT2 + QGT(ig)
                end do
                write(88888, '(A, F10.3, 2A, F10.3, A, F9.3, A, <NGT>(F9.3,","))') " ", JDAY, ",  ", ", ", TDG_TDG, ",  ", SUM_QGT2, ",  ", (QGT(ig), ig = 1, NGT)
                NXTSPLIT3 = NXTSPLIT3 + 1.0
            end if
        end if ! END IF Q_SUM/=0.0
    end subroutine SYSTDG_TDG
!===========================================================================================================================

    subroutine DEALLOCATE_SYSTDG()
        implicit none
        deallocate(TDG_PHS, TDG_FLS, TDG_TDP)
        deallocate(BAYC, QBAY)
        deallocate(POWGTNO, FLGTNO)
        deallocate(GTNAME)
    end subroutine DEALLOCATE_SYSTDG
!

end module modSYSTDG
