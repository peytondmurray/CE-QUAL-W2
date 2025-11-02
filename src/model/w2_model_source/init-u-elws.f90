subroutine initial_water_level()

    use MAIN
    use GLOBAL;     use NAMESC;     use GEOMC;     use LOGICC;     use PREC;     use SURFHE;     use KINETIC;     use SHADEC;     use EDDY
    use STRUCTURES;     use TRANS;     use TVDC;     use SELWC;     use GDAYC;     use SCREENC;     use TDGAS;     use RSTART
    use MACROPHYTEC;     use POROSITYC;     use ZOOPLANKTONC
    use INITIALVELOCITY
    implicit none
    external :: RESTART_OUTPUT
    integer :: JBU, JBD, JJW
    real :: QGATE, WSUP, BRLEN, WLSLOPE, DIST, WLDIFF

    LOOP_BRANCH = .false.

! estimating initial flows in each segment
    QSSI = 0.0

! first considering specified flows: upstream inflows, tributaries, distributed tribs, structural withdrawals and withdrawals    
    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            do I = IU, ID
                if (I == IU .and. UP_FLOW(JB)) then
                    QSSI(I) = QIN(JB) + QSSI(I)
                end if
! if downstream of structure          
!IF (i == iu .and. UHS(JB) < 0) then
                if (I == IU .and. DAM_INFLOW(JB)) then ! CB 4/27/2011
                    do JJW = 1, NWB
                        do JJB = BS(JJW), BE(JJW)
                            if (DS(JJB) == ABS(UHS(JB))) then
                                do JS = 1, NSTR(JJB)
                                    QSSI(I) = QSSI(I) + QSTR(JS, JJB)
                                end do
                            end if
                        end do
                    end do
                end if

                if (TRIBUTARIES) then
                    do JT = 1, NTR
                        if (ITR(JT) == I) then
                            QSSI(I) = QSSI(I) + QTR(JT)
                        end if
                    end do
                end if
                if (DIST_TRIBS(JB)) then
                    QSSI(I) = QSSI(I) + QDTR(JB)/REAL(ID - IU + 1) ! SINCE INITIAL WL UNKNOWN, DISTRIBUTING FLOW EVENLY BTW. SEGS.              
                end if
                if (WITHDRAWALS) then
                    do JWD = 1, NWD
                        if (IWD(JWD) == I) then
                            QSSI(I) = QSSI(I) - QWD(JWD)
                        end if
                    end do
                end if
                if (I == ID) then
                    do JS = 1, NSTR(JB)
                        QSSI(I) = QSSI(I) - QSTR(JS, JB)
                    end do
                end if
            end do
        end do
    end do

! INCLUDING FLOWS WITHIN BRANCH UPSTEAM OF SEGMENT
    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            do I = IU + 1, ID
                QSSI(I) = QSSI(I) + QSSI(I - 1)
            end do
        end do
    end do

    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            do I = IU, ID
! DETERMINING IF SEGMENT IS DOWNSTREAM INTERNAL HEAD BOUNDARY OF ANOTHER *UPSTREAM* BRANCH, AND ADDING FLOW TO SEGMENT AND SEGMENTS DOWNSTREAM
                do JJW = 1, NWB
                    do JJB = BS(JJW), BE(JJW)
                        if (DHS(JJB) == I) then
                            do II = I, ID
                                QSSI(II) = QSSI(II) + QSSI(DS(JJB))
                            end do
                        end if
                    end do
                end do
! DETERMINING IF SEGMENT IS DOWNSTREAM OF SPILLWAY BELOW ANOTHER BRANCH
                do JS = 1, NSP
                    if (ESP(JS) < EL(2, I)) then ! DISREGARDING IF CREST ABOVE GRID
                        if (I == IDSP(JS)) then
                            do II = I, ID
                                QSSI(II) = QSSI(II) + QSSI(IUSP(JS))
                            end do
                        end if
                    end if
                end do
! DETERMINING IF SEGMENT IS DOWNSTREAM OF GATE BELOW ANOTHER BRANCH
                do JG = 1, NGT
                    if (EGT(JG) < EL(2, I)) then ! DISREGARDING IF CREST ABOVE GRID
                        if (I == IDGT(JG)) then
                            if (DYNGTC(JG) == "    FLOW") then
                                QGATE = BGT(JG)
                                do II = I, ID
                                    QSSI(II) = QSSI(II) + QGATE
                                end do
                            else
                                do II = I, ID
                                    QSSI(II) = QSSI(II) + QSSI(IUGT(JG))
                                end do
                            end if
                        end if
                    end if
                end do
            end do
        end do
    end do

! DETERMINING IF BRANCH IS A SECONDARY BRANCH THAT LOOPS AROUND AN ISLAND WITH INTERNAL HEAD BOUNDARIES AT UPSTREAM AND DOWNSTREAM BC
! ATTACHED TO A SINGLE BRANCH
    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            if (UH_INTERNAL(JB) .and. DH_INTERNAL(JB)) then
                do JJW = 1, NWB
                    do JJB = BS(JJW), BE(JJW)
                        if (UHS(JB) > US(JJB) .and. UHS(JB) < DS(JJB)) then
                            JBU = JJB
                        end if
                    end do
                end do
                do JJW = 1, NWB
                    do JJB = BS(JJW), BE(JJW)
                        if (DHS(JB) > US(JJB) .and. DHS(JB) < DS(JJB)) then ! WW 8/19/2013
                            JBD = JJB
                        end if
                    end do
                end do
                LOOP_BRANCH(JB) = JBU == JBD
            end if
        end do
    end do


! GIVEN ESTIMATED FLOWS FOR EACH SEGMENT, ESTIMATING WL WITH NORMAL DEPTH EQUATION
    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
! ONLY CONSIDERING BRANCHES WITH SLOPES > 0      
            if (SLOPE(JB) > 0.0 .and. .not. LOOP_BRANCH(JB)) then
                IU = CUS(JB)
                ID = DS(JB)
                do I = IU, ID
                    call NORMAL_DEPTH(QSSI(I))
                end do
            end if
        end do
    end do

! SMOOTHING WATER SURFACE ELEVATIONS, FIRST WITHIN BRANCHES
    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            if (SLOPE(JB) > 0.0 .and. .not. LOOP_BRANCH(JB)) then
                IU = CUS(JB)
                ID = DS(JB)
                do I = ID - 1, IU, -1
                    if (ELWS(I + 1) > ELWS(I)) then
                        ELWS(I) = ELWS(I + 1)
                    end if
                end do
            end if
        end do
    end do

! SMOOTHING WATER SURFACE ELEVATIONS AT INTERNAL HEAD BOUNDARIES
    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            if (SLOPE(JB) > 0.0 .and. .not. LOOP_BRANCH(JB) .and. DH_INTERNAL(JB)) then
                IU = CUS(JB)
                ID = DS(JB)
                if (ELWS(DHS(JB)) > ELWS(ID)) then
                    ELWS(ID) = ELWS(DHS(JB))
                    do I = ID - 1, IU, -1
                        if (ELWS(I + 1) > ELWS(I)) then
                            ELWS(I) = ELWS(I + 1)
                        end if
                    end do
                end if
            end if
        end do
    end do

! IF SPILLWAY OR GATE AT DOWNSTREAM END OF BRANCH, MAKING SURE WATER LEVEL IS ABOVE CREST ELEVATION
    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            if (SLOPE(JB) > 0.0 .and. .not. LOOP_BRANCH(JB) .and. NSTR(JB) == 0) then
                IU = CUS(JB)
                ID = DS(JB)
! SPILLWAYS          
                do JS = 1, NSP
                    if (ESP(JS) < EL(2, I)) then ! DISREGARDING IF CREST ABOVE GRID
                        if (ID == IUSP(JS)) then
                            WSUP = ESP(JS) + (QSSI(ID)/A1SP(JS))**(1.0/B1SP(JS)) ! ESTIMATING UPSTREAM WS ELEV
                            if (ELWS(ID) < WSUP) then
                                if (IDSP(JS) /= 0) then ! CB 8/10/10
                                    if (ELWS(IDSP(JS)) > WSUP) then
                                        WSUP = ELWS(IDSP(JS))
                                    end if ! CHECKING TO SEE IF DOWNSTREAM WS ELEVATION ISN'T ALREADY 'HIGH'
                                    ELWS(ID) = WSUP
                                end if ! CB 8/10/10
                                do I = ID - 1, IU, -1
                                    if (ELWS(I + 1) > ELWS(I)) then
                                        ELWS(I) = ELWS(I + 1)
                                    end if
                                end do
                            end if
                        end if
                    end if
                end do
! GATES
                do JG = 1, NGT
                    if (EGT(JG) < EL(2, I) .and. DYNGTC(JG) /= "    FLOW") then ! DISREGARDING IF CREST ABOVE GRID
                        if (ID == IUGT(JG)) then
                            if (DYNGTC(JG) == "       B") then
                                WSUP = EGT(JG) + (QSSI(ID)/(A1GT(JG)*BGT(JG)**G1GT(JG)))**(1.0/B1GT(JG))
                            end if
                            if (DYNGTC(JG) == "     ZGT") then
                                WSUP = EGT(JG) + (QSSI(ID)/A1GT(JG))**(1.0/B1GT(JG))
                            end if
                            if (ELWS(ID) < WSUP) then
                                if (ELWS(IDGT(JG)) > WSUP) then
                                    WSUP = ELWS(IDGT(JG))
                                end if ! CHECKING TO SEE IF DOWNSTREAM WS ELEVATION ISN'T ALREADY 'HIGH'  WX 8/21/13
                                ELWS(ID) = WSUP
                                do I = ID - 1, IU, -1
                                    if (ELWS(I + 1) > ELWS(I)) then
                                        ELWS(I) = ELWS(I + 1)
                                    end if
                                end do
                            end if
                        end if
                    end if
                end do
            end if
        end do
    end do

! SMOOTHING WATER LEVEL AROUND LOOP BRANCHES
    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            if (LOOP_BRANCH(JB)) then
                WLDIFF = ELWS(UHS(JB)) - ELWS(UHS(JB))
                BRLEN = 0.0
                do I = IU, ID
                    BRLEN = BRLEN + DLX(I)
                end do
                WLSLOPE = WLDIFF/BRLEN
                DIST = DLX(IU)/2.0
                ELWS(IU) = ELWS(UHS(JB)) - WLSLOPE*DIST
                do I = IU + 1, ID
                    DIST = DIST + (DLX(I - 1) + DLX(I))/2.0
                    ELWS(IU) = ELWS(UHS(JB)) - WLSLOPE*DIST
                end do
            end if
        end do
    end do

    return
end subroutine initial_water_level

!***********************************************************************************************************************************
!**        S U B R O U T I N E    N O R M A L    D E P T H                                                                        **
!***********************************************************************************************************************************


subroutine NORMAL_DEPTH(FLOW)

    use GLOBAL;     use GEOMC

    integer, parameter :: JMAX = 40
    real(R8) :: FLOW
    real :: X1, X2, FUNCVAL1, FUNCVAL2, XACC, FMID, FUNC1, RTBIS, DX, XMID
    integer :: JJ, J

! FIRST, BRACKETING ROOT
    X1 = 0.001
    X2 = 1.0
    call MANNINGS_EQN(FLOW, X1, FUNCVAL1)
    call MANNINGS_EQN(FLOW, X2, FUNCVAL2)

    do JJ = 1, JMAX
        if (FUNCVAL1*FUNCVAL2 > 0.0) then
            if (ABS(FUNCVAL1) < ABS(FUNCVAL2)) then
                X1 = X1/2.0
                call MANNINGS_EQN(FLOW, X1, FUNCVAL1)
            else
                X2 = X2 + 1.5*(X2 - X1)
                call MANNINGS_EQN(FLOW, X2, FUNCVAL2)
            end if
        else
            exit
        end if
    end do

! FINDING ROOT BY BISECTION      
    XACC = 0.01
    call MANNINGS_EQN(FLOW, X2, FMID)
    call MANNINGS_EQN(FLOW, X1, FUNC1)
!    IF(FUNC1*FMID.GE.0.) PAUSE 'ROOT MUST BE BRACKETED IN RTBIS'
    if (FUNC1 < 0.) then
        RTBIS = X1
        DX = X2 - X1
    else
        RTBIS = X2
        DX = X1 - X2
    end if
    do J = 1, JMAX
        DX = DX*.5
        XMID = RTBIS + DX
        call MANNINGS_EQN(FLOW, XMID, FMID)
        if (FMID <= 0.) then
            RTBIS = XMID
        end if
        if (ABS(DX) < XACC .or. FMID == 0.) then
            ELWS(I) = RTBIS + EL(KB(I) + 1, I) ! SW 4/5/13
            return
        end if
    end do
!    PAUSE 'TOO MANY BISECTIONS IN RTBIS'
end subroutine NORMAL_DEPTH


!***********************************************************************************************************************************
!**        S U B R O U T I N E    M A N N I N G S    E Q U A T I O N                                                              **
!***********************************************************************************************************************************


subroutine MANNINGS_EQN(FLOW, DEPTH, FUNCVALUE)

    use GLOBAL;     use GEOMC;     use EDDY;     use LOGICC

    real(R8) :: FLOW
    real :: WSURF, DEPTH, XAREA, WPER, HRAD, FMANN, FUNCVALUE

!     WSURF=EL(KB(I)-1,I)+DEPTH
    WSURF = EL(KB(I) + 1, I) + DEPTH ! CB 7/7/10
    call XSECTIONAL_AREA(WSURF, XAREA)
    WPER = B(KTI(I), I) + 2.0*DEPTH
    HRAD = XAREA/WPER
    if (MANNINGS_N(JW)) then
        FMANN = FRIC(I)
    else
        FMANN = HRAD**0.166666667/FRIC(I)
    end if
    FUNCVALUE = FLOW - XAREA*HRAD**0.6667*SLOPEC(JB)**0.5/FMANN ! SW 4/5/2013

    return
end subroutine MANNINGS_EQN

!***********************************************************************************************************************************
!**        S U B R O U T I N E    C R O S S    S E C T I O N A L    A R E A                                                       **
!***********************************************************************************************************************************


subroutine XSECTIONAL_AREA(WSURF, XAREA)

    use GLOBAL;     use GEOMC;     use MAIN;     use INITIALVELOCITY

    real :: XAREA, WSURF ! 4/5/13 SW
    integer :: KTTOP ! 4/5/13 SW

!     KTTOP = 2
    do K = 2, KMX - 1
        if (EL(K, I) < WSURF) then
            KTTOP = K - 1
            exit
        end if
        KTTOP = K ! CB 8/10/10
    end do
!     DO WHILE (EL(KTTOP,I) > WSURF)
!        KTTOP = KTTOP+1
!     END DO            
    XAREA = (WSURF - EL(KTTOP + 1, I))*BSAVE(KTTOP, I)
    do K = KTTOP + 1, KBI(I)
        XAREA = XAREA + BSAVE(K, I)*H(K, JW)
    end do

    return
end subroutine XSECTIONAL_AREA


!***********************************************************************************************************************************
!**        S U B R O U T I N E    I N I T I A L    H O R I Z O N T A L    V E L O C I T Y                                         **
!***********************************************************************************************************************************


subroutine INITIAL_U_VELOCITY()

    use GLOBAL;     use GEOMC
    use INITIALVELOCITY

    real :: XAREA, WSURF
    integer :: K

    do JW = 1, NWB
        KT = KTWB(JW)

        do JB = BS(JW), BE(JW)
            if (SLOPE(JB) > 0.0 .and. .not. LOOP_BRANCH(JB)) then
                IU = CUS(JB)
                ID = DS(JB)
                do I = IU, ID
                    WSURF = ELWS(I)
                    call XSECTIONAL_AREA(WSURF, XAREA)
                    UAVG(I) = QSSI(I)/XAREA
                    do K = KT, KB(I)
                        U(K, I) = UAVG(I)
                    end do
                end do
            end if
        end do
    end do

    return
end subroutine INITIAL_U_VELOCITY
