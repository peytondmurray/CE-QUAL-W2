subroutine InitTDGtarget()
    use Selective1TDGtarget;     use MAIN;     use modSYSTDG, only: POWNO, FLNO, NBAY, BEGNO, ENDNO, POWGTNO, FLGTNO, NRO, TDGLOC, GTNAME
    use GLOBAL;     use NAMESC;     use GEOMC;     use LOGICC;     use PREC;     use SURFHE;     use KINETIC;     use SHADEC;     use EDDY
    use STRUCTURES;     use TRANS;     use TVDC;     use SELWC;     use GDAYC;     use SCREENC;     use TDGAS;     use RSTART
!
    implicit none
    integer :: it, n, ig
    character(len=8) :: AID1
    character(len=72) :: TGAFN

    targetfnno = 8888
    open(targetfnno + 1, file="w2_TDGtarget.csv", status="old")
!read (targetfnno+1, '(///(8X,A72))') (TITLETDGTARGET(it), it=1,10)
!read (targetfnno+1,'(//8x,2f8.3)') tsfreq, tsconv
!read (targetfnno+1,'(//8x,a8,3f8.3,a8,i8,a8,i8)')tsyearly, tstsrt, tstend, tstarget, tsdynsel, tsiteration, dyupdate, dygroup   
    read(targetfnno + 1, *)
    read(targetfnno + 1, *)
    read(targetfnno + 1, *)
    do it = 1, 10
        read(targetfnno + 1, *) TITLETDGTARGET(it)
    end do
    read(targetfnno + 1, *)
    read(targetfnno + 1, *)
    read(targetfnno + 1, *) AID1, tsfreq, tsconv
    read(targetfnno + 1, *)
    read(targetfnno + 1, *)
    read(targetfnno + 1, *) AID1, tsyearly, tstsrt, tstend, tstarget, tsdynsel, tsiteration, dyupdate, dygroup

    tsyearly = ADJUSTR(tsyearly);     tsdynsel = ADJUSTR(tsdynsel);     dyupdate = ADJUSTR(dyupdate)
    if (tsyearly == "     OFF" .and. tstsrt < TMSTRT) then
        tstsrt = TMSTRT
    end if
    NXTSPLIT = TMSTRT
    NXTSPLIT2 = TMSTRT
    if (tsyearly == "     OFF") then
        DAYTEST = jday
    else
        DAYTEST = real(jdayg) + jday - int(jday)
    end if
    if (NXTSPLIT > TMSTRT) then
        NXTSPLIT = TMSTRT
    end if
    if (DAYTEST <= tstsrt) then
        DAYTEST = tstsrt
    end if
    if (tstsrt > TMSTRT) then
        NXTSPLIT = tstsrt
    end if
    if (NXTSPLIT2 > TMSTRT) then
        NXTSPLIT2 = TMSTRT
    end if
    if (tstsrt > TMSTRT) then
        NXTSPLIT2 = tstsrt
    end if
    NGSP = NBAY + NRO
    NGPH = POWNO
    NGFL = FLNO
    NOUTS = NGT
    NGSPPH = NGSP + NGPH
    allocate(SPGTNO(NGSP), SPPRIOR(NGSP), SPMINFRAC(NGSP))
    if (NGPH > 0) then
        allocate(PHGTNO(NGPH), PHMAXFLOW(NGPH))
    end if
    it = 0
    do ig = 1, NGT
        if (GTNAME(ig)) then
            it = it + 1
            SPGTNO(it) = ig
        end if
    end do

    if (NGPH > 0) then
        it = 0
        do ig = 1, NGT
            if (GTTYP(ig) == "     POW") then
                it = it + 1
                PHGTNO(it) = ig
            end if
        end do
    end if

!read (targetfnno+1,'(//8x,<NGSP>i8)') (SPPRIOR(n),n=1,NGSP)  
!read (targetfnno+1,'(//8x,<NGSP>f8.3)') (SPMINFRAC(n),n=1,NGSP)
!if (NGPH>0) read (targetfnno+1,'(//8x,<NGPH>f8.3)') (PHMAXFLOW(n),n=1,NGPH) 
    read(targetfnno + 1, *)
    read(targetfnno + 1, *)
    read(targetfnno + 1, *) AID1, (SPPRIOR(n), n = 1, NGSP)
    read(targetfnno + 1, *)
    read(targetfnno + 1, *)
    read(targetfnno + 1, *) AID1, (SPMINFRAC(n), n = 1, NGSP)
    if (NGPH > 0) then
        read(targetfnno + 1, *)
        read(targetfnno + 1, *)
        read(targetfnno + 1, *) AID1, (PHMAXFLOW(n), n = 1, NGPH)
    end if

    do n = 1, NGSP
        if (SPMINFRAC(n) > 1.0) then
            SPMINFRAC(n) = 1.0
        end if ! remove unrealistic input value
    end do
    do n = 1, NGPH
        if (PHMAXFLOW(n) < 0.0) then
            PHMAXFLOW(n) = 0.0
        end if ! remove unrealistic input value
    end do
    if (tsconv < 0.001) then
        tsconv = 0.001
    end if ! constrain the convergence criterion to be >= 0.001 and <= 5.0
    if (tsconv > 5.0) then
        tsconv = 5.0
    end if

! OPEN DYNAMIC TDG TARGET FILES    
    if (tsdynsel == "      ON") then
!read (targetfnno+1,'(//(8X,A72))') TGAFN
        read(targetfnno + 1, *)
        read(targetfnno + 1, *)
        read(targetfnno + 1, *) AID1, TGAFN
        open(targetfnno + 3, file=TGAFN, status="old")
        read(targetfnno + 3, *)
        read(targetfnno + 3, *)
        read(targetfnno + 3, *)
        read(targetfnno + 3, *) nxtjday, tstarget2
        tstarget = tstarget2
        read(targetfnno + 3, *) nxtjday, tstarget2
    end if
    close(targetfnno + 1)

!!  Initial output file  
    open(targetfnno, FILE="TDGTarget_output.csv", status="unknown")
!if (NGT>0) write (targetfnno,'(3A,<NGT>(A,i2))')'     JDAY','         TDG' ,'      SUM Q   ',('      Q',n, n=1,NGT)
    if (NGT > 0) then
        write(targetfnno, '("     JDAY,", " C,", "       TDG,", "    SUM Q,", <NGT>(A,i2,","))') ("      Q", n, n = 1, NGT)
    end if
    open(targetfnno + 2, FILE="TDGTarget_warning.opt", action="READWRITE", status="unknown")
    return
end subroutine InitTDGtarget

!***********************************************************************************************************************************
!**                                                 T D G    T A R G E T   S U B R O U T I N E                                    **
!***********************************************************************************************************************************

subroutine TDGtarget()
    use Selective1TDGtarget;     use modSYSTDG, only: TDG_TDG, TDGLOC, SYSTDG_TDG;     use MAIN, only: targetfnno, WARNING_OPEN
    use GLOBAL;     use NAMESC;     use GEOMC;     use LOGICC;     use PREC;     use SURFHE;     use KINETIC;     use SHADEC;     use EDDY
    use STRUCTURES;     use TRANS;     use TVDC;     use SELWC;     use GDAYC;     use SCREENC;     use TDGAS;     use RSTART; 
! 
    implicit none
!
    real, allocatable, dimension(:) :: QGTSAVE
    logical, allocatable, dimension(:) :: SP_ACTIVE, PH_ACTIVE
    integer :: ig, ii, ITERATION, prior_top, priortop_n
    real :: Q_ALL, Q_SP, Q_PH, QSP_AVL, QPH_AVL, Q_TEMP
    real :: SUM_SP_FRAC, SUM_PH_MAXFLOW
    real :: Q_MAX, Q_MIN, Q_CUT, Q_CUTTED, Q_LEFT, QPH_ADDED
    real :: SUM_TOP_FLOW, SUM_TOP_MINFRAC
    real :: SUM_QGT, SUM_QGT2
    real :: MINV = 0.0000000001
    integer, allocatable, dimension(:) :: priortop_spno
    character(len=1) :: CO

    allocate(QGTSAVE(NGT), SP_ACTIVE(NGSP))
    if (NGPH > 0) then
        allocate(PH_ACTIVE(NGPH))
    end if
    QGTSAVE = QGT
    SUM_QGT = 0.0
    do ig = 1, NGT
        SUM_QGT = SUM_QGT + QGT(ig)
    end do
    SP_ACTIVE = .false.
    PH_ACTIVE = .false.
!  
    SUM_SP_FRAC = 0.0
    do ig = 1, NGSP
        if (QGT(SPGTNO(ig)) > 0.0) then
            SUM_SP_FRAC = SUM_SP_FRAC + SPMINFRAC(ig)
        end if
    end do
    if (SUM_SP_FRAC > 1.0) then ! IF SUM SP MIN FRAC IS BIGGER THAN 1.0 THEN ADJUST THE SPMINFRAC BY PERCENTAGE
        do ig = 1, NGSP
            if (QGT(SPGTNO(ig)) > 0.0) then
                SPMINFRAC(ig) = SPMINFRAC(ig)*1.0/SUM_SP_FRAC
            end if
        end do
        SUM_SP_FRAC = 1.0
    end if
!  
    SUM_PH_MAXFLOW = 0.0
    do ig = 1, NGPH
        if (QGT(PHGTNO(ig)) > 0.0) then
            SUM_PH_MAXFLOW = SUM_PH_MAXFLOW + PHMAXFLOW(ig)
        end if
    end do
!  
    Q_SP = 0.0
    do ig = 1, NGSP
        Q_SP = Q_SP + QGT(SPGTNO(ig))
    end do

    if (Q_SP > 0.0) then
        if (tsyearly == "     OFF") then
            DAYTEST = JDAY
        else
            DAYTEST = real(JDAYG) + JDAY - int(JDAY)
        end if

        call SYSTDG_TDG() ! CALCULATE CURRENT TDG
        if (tsdynsel == "      ON") then
            if (DAYTEST > tstsrt) then
                tstarget = tstarget2
            end if
            do while (JDAY >= nxtjday)
                read(targetfnno + 3, *) nxtjday, tstarget2
                tstarget = tstarget2
            end do
        end if

        CO = " "
        if (DAYTEST >= tstsrt .and. DAYTEST <= tstend .and. TDG_TDG > tstarget + tsconv) then
! INITIAL VARIABLES FOR ITERATIONS
            CO = "R"
            Q_SP = 0.0
            Q_PH = 0.0
            Q_ALL = 0.0
            do ig = 1, NGT
                Q_ALL = Q_ALL + QGT(ig)
            end do
            do ig = 1, NGSP
                Q_SP = Q_SP + QGT(SPGTNO(ig))
                if (QGT(SPGTNO(ig)) > 0.0 .and. QGT(SPGTNO(ig)) > Q_ALL*SPMINFRAC(SPGTNO(ig)) + MINV) then
                    SP_ACTIVE(ig) = .true.
                end if
            end do
            do ig = 1, NGPH
                Q_PH = Q_PH + QGT(PHGTNO(ig))
                if (QGT(PHGTNO(ig)) > 0.0 .and. QGT(SPGTNO(ig)) + MINV < PHMAXFLOW(ig)) then
                    PH_ACTIVE(ig) = .true.
                end if
            end do

            QSP_AVL = Q_SP - Q_ALL*SUM_SP_FRAC
            if (SUM_PH_MAXFLOW > 0.0) then
                QPH_AVL = SUM_PH_MAXFLOW - Q_PH
            end if
            if (QPH_AVL < 0.0) then
                QPH_AVL = 0.0
            end if
            Q_MAX = Q_SP
            if (TDGLOC == "     REL" .and. SUM_PH_MAXFLOW /= 0.0 .and. QPH_AVL > 0.0) then
                Q_MAX = MIN(Q_SP, QPH_AVL)
            end if ! IF RELEASE TDG, FLOW CUT TO POWERHOUSE MUST LESS THAN MAX FLOW
            Q_MIN = MAX(0.0, Q_ALL*SUM_SP_FRAC)
            if (TDGLOC == "     REL" .and. SUM_PH_MAXFLOW /= 0.0) then
                Q_MIN = MAX(Q_ALL*SUM_SP_FRAC, Q_SP - QPH_AVL, 0.0)
            end if
            Q_CUT = 0.5*(Q_MAX + Q_MIN)
            ITERATION = 1
            if (JDAY >= NXTSPLIT2 .and. dyupdate == "      ON") then
                call Dy_Priority()
            end if ! DYNAMIC PRIORITY UPDATE

            do while (ABS(TDG_TDG - tstarget) > tsconv .and. ITERATION <= tsiteration .and. Q_CUT > 0.0)
                Q_CUTTED = 0.0
                Q_LEFT = Q_SP - Q_CUTTED
                do ig = 1, NGSP
                    if (QGT(SPGTNO(ig)) > 0.0 .and. QGT(SPGTNO(ig)) > Q_ALL*SPMINFRAC(SPGTNO(ig)) + MINV) then
                        SP_ACTIVE(ig) = .true.
                    end if
                end do
                do ig = 1, NGPH
                    if (QGT(PHGTNO(ig)) > 0.0 .and. QGT(SPGTNO(ig)) + MINV < PHMAXFLOW(ig)) then
                        PH_ACTIVE(ig) = .true.
                    end if
                end do
                Q_TEMP = 0.0

                do while (Q_CUTTED + MINV < Q_CUT .and. Q_LEFT > Q_ALL*SUM_SP_FRAC + MINV) ! CUT FLOW TO SP UNTILL Q_CUTTED = Q_CUT
! UPDATE PRIOR TOP
                    prior_top = -999 ! highest prior
                    do ig = 1, NGSP
                        if (prior_top == -999 .or. SPPRIOR(ig) < prior_top) then
                            if (SP_ACTIVE(ig) .and. SPPRIOR(ig) > 0) then
                                prior_top = SPPRIOR(ig)
                            end if
                        end if
                    end do
                    priortop_n = 0
                    do ig = 1, ngsp ! initial highest prior out NO
                        if (SPPRIOR(ig) == prior_top .and. SP_ACTIVE(ig)) then
                            priortop_n = priortop_n + 1
                        end if
                    end do
                    if (priortop_n > 0) then
                        allocate(priortop_spno(priortop_n))
                    end if
                    ii = 0
                    SUM_TOP_FLOW = 0.0
                    SUM_TOP_MINFRAC = 0.0
                    do ig = 1, NGSP ! initial highest prior GTNO, SUM FLOW AND SUM MINFRAC
                        if (SPPRIOR(ig) == prior_top .and. SP_ACTIVE(ig)) then
                            ii = ii + 1
                            priortop_spno(ii) = ig
                            SUM_TOP_FLOW = SUM_TOP_FLOW + QGT(SPGTNO(ig))
                            SUM_TOP_MINFRAC = SUM_TOP_MINFRAC + SPMINFRAC(ig)
                        end if
                    end do

                    if (SUM_TOP_FLOW - Q_ALL*SUM_TOP_MINFRAC <= Q_CUT - Q_CUTTED) then ! top prior flow is not enough to cut all Q_CUT
                        do ig = 1, priortop_n
                            Q_CUTTED = Q_CUTTED + QGT(SPGTNO(priortop_spno(ig))) - Q_ALL*SPMINFRAC(priortop_spno(ig))
                            QGT(SPGTNO(priortop_spno(ig))) = Q_ALL*SPMINFRAC(priortop_spno(ig))
                            SP_ACTIVE(priortop_spno(ig)) = .false.
                        end do
                        Q_LEFT = Q_SP - Q_CUTTED
                    else
                        Q_TEMP = Q_CUT - Q_CUTTED ! top prior available flow if more than Q_CUT-Q_CUTTED, so cut it by flow percentage
                        do ig = 1, priortop_n
                            Q_CUTTED = Q_CUTTED + Q_TEMP*(QGT(SPGTNO(priortop_spno(ig))) - Q_ALL*SPMINFRAC(priortop_spno(ig)))/(SUM_TOP_FLOW - Q_ALL*SUM_TOP_MINFRAC)
                            QGT(SPGTNO(priortop_spno(ig))) = QGT(SPGTNO(priortop_spno(ig))) - Q_TEMP*(QGT(SPGTNO(priortop_spno(ig))) - Q_ALL*SPMINFRAC(priortop_spno(ig)))/(SUM_TOP_FLOW - Q_ALL*SUM_TOP_MINFRAC)
                        end do
                        Q_LEFT = Q_SP - Q_CUTTED
                    end if

                    do ig = 1, priortop_n
                        if (QGT(SPGTNO(priortop_spno(ig))) > Q_ALL*SPMINFRAC(priortop_spno(ig)) + MINV) then
                            SP_ACTIVE(priortop_spno(ig)) = .true.
                        end if
                    end do
                    if (priortop_n > 0) then
                        deallocate(priortop_spno)
                    end if
! THIS PRIOR TOP ENDED
                end do ! END DO WHILE LOOP FOR Q_CUT

                if (Q_LEFT > Q_SP - Q_CUT) then ! KEEP FLOW BALANCE
                    do ig = 1, NGSP
                        if (QGT(SPGTNO(ig)) > Q_ALL*SPMINFRAC(ig) + Q_LEFT + Q_CUT - Q_SP .and. SP_ACTIVE(ig)) then
                            QGT(SPGTNO(ig)) = QGT(SPGTNO(ig)) - (Q_LEFT + Q_CUT - Q_SP)
                            Q_LEFT = Q_SP - Q_CUT
                            Q_CUTTED = Q_CUT
                            exit
                        end if
                    end do
                end if

! CUT FLOW TO POWERHOUSE
                QPH_ADDED = 0.0
                if (SUM_PH_MAXFLOW == 0.0 .and. Q_CUTTED >= Q_CUT) then ! NO MAX FLOW LIMIT FOR ALL POWERHOUSE, FLOW ADDED BY PERCENTAGE        
                    do ig = 1, NGPH
                        if (QGT(PHGTNO(ig)) > 0.0) then
                            QPH_ADDED = QPH_ADDED + Q_CUTTED*QGT(PHGTNO(ig))/Q_PH
                            QGT(PHGTNO(ig)) = QGT(PHGTNO(ig)) + Q_CUTTED*QGT(PHGTNO(ig))/Q_PH
                        end if
                    end do
                else
                    if (SUM_PH_MAXFLOW /= 0.0 .and. Q_CUTTED >= Q_CUT) then
                        do ig = 1, NGPH
                            if (PHMAXFLOW(ig) /= 0.0 .and. QGT(PHGTNO(ig)) > 0.0) then
                                QPH_ADDED = QPH_ADDED + Q_CUTTED*(PHMAXFLOW(ig) - QGT(PHGTNO(ig)))/QPH_AVL
                                QGT(PHGTNO(ig)) = QGT(PHGTNO(ig)) + Q_CUTTED*(PHMAXFLOW(ig) - QGT(PHGTNO(ig)))/QPH_AVL
                            else
                                if (QGT(PHGTNO(ig)) > 0.0) then
                                    QGT(PHGTNO(ig)) = QGT(PHGTNO(ig)) + Q_CUTTED
                                    QPH_ADDED = Q_CUTTED
                                end if
                            end if
                        end do
                    end if
                end if

                if (QPH_ADDED < Q_CUTTED) then ! FLOW BALANCE
                    do ig = 1, NGPH
                        if (QGT(PHGTNO(ig)) > 0.0 .and. QGT(PHGTNO(ig)) + Q_CUTTED - QPH_ADDED < PHMAXFLOW(ig)) then
                            QGT(PHGTNO(ig)) = QGT(PHGTNO(ig)) + Q_CUTTED - QPH_ADDED
                            QPH_ADDED = Q_CUTTED
                            exit
                        end if
                    end do
                end if

                call SYSTDG_TDG()
                if (ABS(TDG_TDG - tstarget) <= tsconv) then
                    exit
                end if
                if (TDG_TDG - tstarget > tsconv) then
                    Q_MIN = Q_CUT
                    QGT = QGTSAVE ! NOT CONVERGENT
                    call SYSTDG_TDG()
                else
                    if (tstarget - TDG_TDG > tsconv) then
                        Q_MAX = Q_CUT
                        QGT = QGTSAVE ! NOT CONVERGENT
                        call SYSTDG_TDG()
                    end if
                end if
                Q_CUT = 0.5*(Q_MAX + Q_MIN)
                ITERATION = ITERATION + 1
!
            end do ! END DO WHILE FOR DICHONOMY FLOW CUT

            if (ITERATION == 1 .and. Q_SP == 0.0) then
                write(targetfnno + 2, "(A,F12.3)") "SPILL FLOW IS ZERO ON JDAY", JDAY
            end if
            if (ITERATION == 1 .and. TDG_TDG < tstarget + tsconv) then
                write(targetfnno + 2, "(A,F12.3)") "TDG IS LOWER THAN TARGET, NO ITERATION CALCULATION NEEDED ON JDAY", JDAY
            end if
            if (TDG_TDG - tstarget > tsconv .and. ITERATION > tsiteration .and. Q_SP > 0.0) then
                write(targetfnno + 2, "(A,F12.3)") "THE ITERATION LIMIT IS EXCEEDED, AND TDG HAS STILL NOT CONVERGED TO TARGET ON JDAY", JDAY
                WARNING_OPEN = .true.
                QGT = QGTSAVE ! UNDO THE FLOW CUT
                call SYSTDG_TDG()
            end if
!
        end if

        if (JDAY >= NXTSPLIT) then
            if (TDG_TDG - tstarget > tsconv) then
                CO = "U"
            else
                if (CO /= "R") then
                    CO = " "
                end if
            end if
            SUM_QGT2 = 0.0
            do ig = 1, NGT
                SUM_QGT2 = SUM_QGT2 + QGT(ig)
            end do
!WRITE (targetfnno, '(A, F10.3, 2A, F10.3, A, F9.3, A, <NGT>(F9.3))')' ',JDAY,'  ', CO,TDG_TDG,'  ',SUM_QGT2,'  ',(QGT(ig), ig = 1, NGT)
            write(targetfnno, '(F10.3, ",", A, ",", F10.3, ",", F9.3, ",", <NGT>(F9.3,","))') JDAY, CO, TDG_TDG, SUM_QGT2, (QGT(ig), ig = 1, NGT)
            NXTSPLIT = NXTSPLIT + tsfreq
        end if
        if (JDAY >= NXTSPLIT2) then
            NXTSPLIT2 = NXTSPLIT2 + tsfreq
        end if
!
    end if
    deallocate(QGTSAVE, SP_ACTIVE)
    if (NGPH > 0) then
        deallocate(PH_ACTIVE)
    end if
end subroutine TDGtarget


subroutine DEALLOCATE_TDGtarget()
    use Selective1TDGtarget;     use MAIN, only: targetfnno
    implicit none
!
    close(targetfnno)
    close(targetfnno + 2)
    close(targetfnno + 3)
    if (tsdynsel == "      ON") then
        close(targetfnno + 3)
    end if
    deallocate(SPGTNO, SPPRIOR, SPMINFRAC)
    if (NGPH > 0) then
        deallocate(PHGTNO, PHMAXFLOW)
    end if
end subroutine DEALLOCATE_TDGtarget


subroutine Dy_Priority()
    use Selective1TDGtarget;     use STRUCTURES, only: QGT
    implicit none
!  
    integer :: ib, IBB, ig, TOP_SPNO, TOP_PRIOR, COUNT, CONTU, NGROUP, NLEFT
    integer, dimension(NGSP) :: PRIOR_SP
    logical, dimension(NGSP) :: PRIOR_UPDATED
    real, allocatable, dimension(:) :: QSP_TEMP
    integer, allocatable, dimension(:) :: GTNO
    TOP_PRIOR = 1
    PRIOR_SP = NGSP + 1
    PRIOR_UPDATED = .false.

    COUNT = 0
    do ib = 1, NGSP
        if (QGT(SPGTNO(ib)) > 0.0) then
            COUNT = COUNT + 1
        else
            PRIOR_UPDATED(ib) = .true.
        end if
    end do
    if (COUNT > 0) then
        allocate(GTNO(COUNT))
    end if
    if (COUNT > 0) then
        allocate(QSP_TEMP(COUNT))
    end if
    GTNO = 0
    COUNT = 0
    do ib = 1, NGSP
        if (QGT(SPGTNO(ib)) > 0.0) then
            COUNT = COUNT + 1
            GTNO(COUNT) = ib
            QSP_TEMP(COUNT) = QGT(SPGTNO(ib))
        end if
    end do
    if (COUNT == 1) then
        PRIOR_SP(GTNO(1)) = 1
    end if
    if (COUNT > 1) then
        call BUBBLE_SORT(QSP_TEMP, GTNO, COUNT)
        PRIOR_SP(GTNO(1)) = 1
        do ib = 1, COUNT - 1
            if (QSP_TEMP(ib) /= QSP_TEMP(ib + 1)) then
                TOP_PRIOR = TOP_PRIOR + 1
                PRIOR_SP(GTNO(ib + 1)) = TOP_PRIOR
            else
                PRIOR_SP(GTNO(ib + 1)) = TOP_PRIOR
            end if
        end do
    end if

    if (dygroup > 1) then ! GROUP THE PRIORITY BY DYGROUP
        NGROUP = INT(COUNT/dygroup)
        NLEFT = MOD(COUNT, dygroup)
        do ib = 1, NGROUP
            do ig = 1, dygroup
                PRIOR_SP(GTNO((ib - 1)*dygroup + ig)) = ib
            end do
        end do
        if (NLEFT > 0) then
            do ig = 1, NLEFT
                PRIOR_SP(GTNO(NGROUP*dygroup + ig)) = NGROUP + 1
            end do
        end if
    end if
    if (COUNT > 0) then
        deallocate(GTNO, QSP_TEMP)
    end if
    SPPRIOR = PRIOR_SP
    return
end subroutine Dy_Priority


subroutine BUBBLE_SORT(QSP_TEMP, GTNO, COUNT)
    implicit none
    integer :: COUNT
    real, dimension(COUNT) :: QSP_TEMP
    integer, dimension(COUNT) :: GTNO
    integer :: I, J
    real :: TEMP, TEMP_NO

    do I = COUNT - 1, 1, -1
        do J = 1, I
            if (QSP_TEMP(J) < QSP_TEMP(J + 1)) then
                TEMP = QSP_TEMP(J)
                QSP_TEMP(J) = QSP_TEMP(J + 1)
                QSP_TEMP(J + 1) = TEMP

                TEMP_NO = GTNO(J)
                GTNO(J) = GTNO(J + 1)
                GTNO(J + 1) = TEMP_NO
            end if
        end do
    end do
    return
end subroutine BUBBLE_SORT
