subroutine ENDSIMULATION()

    use MAIN
    use GLOBAL;     use NAMESC;     use GEOMC;     use LOGICC;     use PREC;     use SURFHE;     use KINETIC;     use SHADEC;     use EDDY
    use STRUCTURES;     use TRANS;     use TVDC;     use SELWC;     use GDAYC;     use SCREENC;     use TDGAS;     use RSTART
    use MACROPHYTEC;     use POROSITYC;     use ZOOPLANKTONC;     use INITIALVELOCITY;     use BIOENERGETICS;     use TRIDIAG_V
    use modSYSTDG, only: DEALLOCATE_SYSTDG !  systdg
! CEMA testing start
    use CEMAVars
! CEMA testing end
    implicit none
    external :: RESTART_OUTPUT
    integer :: IFILE
!***********************************************************************************************************************************
!*                                                    Task 3: End Simulation                                                      **
!***********************************************************************************************************************************

! CEMA testing start
!real jcs
!!jcs=jcinzz*2.67/1000.0  ! converting C flux to DO flux, assuming 2.67 gO/gC
!jcs=SD_jctest*2.67  ! converting C flux to DO flux, assuming 2.67 gO/gC
! !write(1081,'(5g12.5)')xjnh4,Jcinzz/1000.0,Jcs,MFTSedFlxVars(2,26)
! write(1081,'(5g12.5)')xjnh4,SD_Jctest,Jcs,MFTSedFlxVars(2,26)
! CEMA testing end

    call DATE_AND_TIME(CDATE, CCTIME)
    if (.not. ERROR_OPEN) then
        TEXT = "Normal termination at " // CCTIME(1:2) // ":" // CCTIME(3:4) // ":" // CCTIME(5:6) // " on " // CDATE(5:6) // "/" // CDATE(7:8) // "/" // CDATE(3:4)
    end if
    TEXT = ADJUSTL(TRIM(TEXT))
    call CPU_TIME(CURRENT)
    do JW = 1, NWB
        if (SNAPSHOT(JW)) then
            write(SNP(JW), "(/A/)") ADJUSTL(TRIM(TEXT))
            write(SNP(JW), "(A)") "Runtime statistics"
            write(SNP(JW), "(2(A,I0))") "  Grid                 = ", IMX, " x ", KMX
            write(SNP(JW), "(A,I0)") "  Maximum active cells = ", NTACMX, "  Minimum active cells = ", NTACMN
            write(SNP(JW), "(3(A,F0.1))") "  Segment lengths, m   = ", DLXMIN, "-", DLXMAX
            write(SNP(JW), "(3(A,F0.1))") "  Layer heights, m     = ", HMIN, "-", HMAX
            write(SNP(JW), "(A)") "  Timestep"
            write(SNP(JW), "(A,I0)") "    Total iterations   = ", NIT
            write(SNP(JW), "(A,I0)") "    # of violations    = ", NV
            write(SNP(JW), "(A,F0.2)") "    % violations       = ", FLOAT(NV)/FLOAT(NIT)*100.0
            write(SNP(JW), "(A,I0,A)") "    Average timestep   = ", INT(DLTAV), " sec"
            write(SNP(JW), "(A,I0,A,F0.2,A)") "  Simulation time      = ", INT(ELTMJD), " days ", (ELTMJD - INT(ELTMJD))*24.0, " hours"
            write(SNP(JW), "(A,F0.2,A)") "  Total CPU runtime    = ", (CURRENT - START)/60.0, " min"
            close(SNP(JW))
        end if
!IF (VECTOR(JW))      CLOSE (VPL(JW))
        if (PROFILE(JW)) then
            close(PRF(JW))
        end if
        if (SPREADSHEET(JW)) then
            close(SPR(JW))
        end if
        if (CONTOUR(JW)) then
            close(CPL(JW))
        end if
    end do

! *** DSI W2_TOOL LINKAGE
    if (VECTOR(1)) then
        close(VPL(1))
    end if

    if (TIME_SERIES) then
        do J = 1, NIKTSR
            close(TSR(J))
        end do
        close(WLFN) ! WL output file  ! SW 9/25/13
    end if
    if (WARNING_OPEN) then
        close(WRN)
    else
        close(WRN, STATUS="DELETE")
    end if
    if (ERROR_OPEN) then
        close(W2ERR)
    else
        close(W2ERR, STATUS="DELETE")
    end if
    do J = 40, NOPEN
        close(J)
    end do

    if (FLOWBALC == "      ON") then
        close(FLOWBFN) ! flowbal file
    end if

    if (NPBALC == "      ON") then
        close(MASSBFN) ! MASS BALANCE file
    end if

    if (SELECTC == "      ON") then ! SW 9/25/13 New Section on closing files
        ifile = 1949
        do jb = 1, nbr
            if (nstr(jb) > 0) then
                ifile = ifile + 1
                close(ifile)
            end if
        end do
        if (nwd > 0) then
            ifile = ifile + 1
            close(ifile)
        end if
        do jw = 1, nwb ! sw 4/20/15
            ifile = ifile + 1
            close(ifile)
        end do

    end if

    if (DOWNSTREAM_OUTFLOW) then
        JFILE = 0
        do JWD = 1, NIWDO
            close(WDO(JWD, 1))
            close(WDO(JWD, 2))
            if (CONSTITUENTS) then
                close(WDO(JWD, 3))
            end if
            if (DERIVED_CALC) then
                close(WDO(JWD, 4))
            end if

! Determine the # of withdrawals at the WITH SEG  
            do JB = 1, NBR ! structures
                if (IWDO(JWD) == DS(JB) .and. NSTR(JB) /= 0) then
                    do JS = 1, NSTR(JB)
                        JFILE = JFILE + 1
                        close(WDO2(JFILE, 1))
                        close(WDO2(JFILE, 2))
                        if (CONSTITUENTS) then
                            close(WDO2(JFILE, 3))
                        end if
                        if (DERIVED_CALC) then
                            close(WDO2(JFILE, 4))
                        end if
                    end do
                end if
            end do

            do JS = 1, NWD ! withdrawals
                if (IWDO(JWD) == IWD(JS)) then
                    JFILE = JFILE + 1
                    close(WDO2(JFILE, 1))
                    close(WDO2(JFILE, 2))
                    if (CONSTITUENTS) then
                        close(WDO2(JFILE, 3))
                    end if
                    if (DERIVED_CALC) then
                        close(WDO2(JFILE, 4))
                    end if
                end if
            end do

            do JS = 1, NSP ! spillways
                if (IWDO(JWD) == IUSP(JS)) then
                    JFILE = JFILE + 1
                    close(WDO2(JFILE, 1))
                    close(WDO2(JFILE, 2))
                    if (CONSTITUENTS) then
                        close(WDO2(JFILE, 3))
                    end if
                    if (DERIVED_CALC) then
                        close(WDO2(JFILE, 4))
                    end if
                end if
            end do

            do JS = 1, NPU ! pumps
                if (IWDO(JWD) == IUPU(JS)) then
                    JFILE = JFILE + 1
                    close(WDO2(JFILE, 1))
                    close(WDO2(JFILE, 2))
                    if (CONSTITUENTS) then
                        close(WDO2(JFILE, 3))
                    end if
                    if (DERIVED_CALC) then
                        close(WDO2(JFILE, 4))
                    end if
                end if
            end do

            do JS = 1, NPI ! pipes
                if (IWDO(JWD) == IUPI(JS)) then
                    JFILE = JFILE + 1
                    close(WDO2(JFILE, 1))
                    close(WDO2(JFILE, 2))
                    if (CONSTITUENTS) then
                        close(WDO2(JFILE, 3))
                    end if
                    if (DERIVED_CALC) then
                        close(WDO2(JFILE, 4))
                    end if
                end if
            end do

            do JS = 1, NGT ! gates
                if (IWDO(JWD) == IUGT(JS)) then
                    JFILE = JFILE + 1
                    close(WDO2(JFILE, 1))
                    close(WDO2(JFILE, 2))
                    if (CONSTITUENTS) then
                        close(WDO2(JFILE, 3))
                    end if
                    if (DERIVED_CALC) then
                        close(WDO2(JFILE, 4))
                    end if
                end if
            end do

        end do
    end if



!OPEN(W2ERR,FILE='W2Errordump.opt',status='unknown')
!WRITE(w2err,*)'JDAY',jday,'SZ',sz,'Z',z,'H2KT',h2(kt,1:imx),'H1KT',h1(kt,1:imx),'BHR1',bhr1(kt,1:imx),'BHR2',bhr2(kt,1:imx),'WSE',elws,'Q',q,'QC',qc,'QERR',qerr,'T1',t1(kt,1:imx),'T2',t2(kt,1:imx),'SUKT',su(kt,1:imx),&
!                        'UKT',u(kt,1:imx),'QIN',qin,'QTR',qtr,'QWD',qwd
    if (ERROR_OPEN) then ! modified to be more organized and comma-delimited  !SR 12/26/2019
        open(W2ERR, FILE="W2Errordump.csv", status="unknown") ! changed to csv                                     !SR 12/26/2019
        write(W2ERR, *) "JDAY = ", JDAY !SR 12/26/2019
        write(W2ERR, '(A,1000(",",F0.6))') "QIN:", (QIN(J), J = 1, NBR) !SR 12/26/2019
        write(W2ERR, '(A,1000(",",F0.6))') "QTR:", (QTR(J), J = 1, NTRT) !SR 12/26/2019
        write(W2ERR, '(A,1000(",",F0.6))') "QDT:", (QDTR(J), J = 1, NBR) !SR 12/26/2019
        write(W2ERR, '(A,1000(",",F0.6))') "QWD:", (QWD(J), J = 1, NWDT) !SR 12/26/2019
        write(W2ERR, "(/A)") "SEG,BRANCH,KT,WSE,SZ,Z,Q,QC,QERR,H2KT,H1KT,BHR1,BHR2,T1,T2,SUKT,UKT" !SR 12/26/2019
        do JW = 1, NWB !SR 12/26/2019
            KT = KTWB(JW) !SR 12/26/2019
            do JB = BS(JW), BE(JW) !SR 12/26/2019
                do I = US(JB) - 1, DS(JB) + 1 !SR 12/26/2019
                    write(W2ERR, '(I0,",",I0,",",I0,14(",",F0.6))') I, JB, KT, ELWS(I), SZ(I), Z(I), Q(I), QC(I), QERR(I), H2(KT, I), H1(KT, I), BHR1(KT, I), BHR2(KT, I), T1(KT, I), T2(KT, I), SU(KT, I), U(KT, I) !SR 12/26/2019
                end do !SR 12/26/2019
            end do !SR 12/26/2019
        end do !SR 12/26/2019
        close(W2ERR)
    end if

    deallocate(LAYERCHANGE, TECPLOT, X1, HAB)
    deallocate(TSR, WDO, WDO2, ETSR, IWDO, ITSR, TITLE, CDAC, WSC, ESTR, WSTR, QSTR, KTSW, KBSW, SINKC, JBTSR, CONSTRICTION, BCONSTRICTION) ! SW 7/24/2018  8/5/2018
    deallocate(EBC, MBC, PQC, EVC, PRC, WINDC, QINC, QOUTC, HEATC, SLHTC, QINIC, DTRIC, TRIC, WDIC)
    deallocate(EXC, EXIC, VBC, METIC, SLTRC, THETA, FRICC, NAF, ELTMF, ZMIN, IZMIN, C2CH, CDCH, EPCH, KFCH, APCH, ANCH, ALCH)
    deallocate(CPLTC, HPLTC, CMIN, CMAX, HYMIN, HYMAX, CDMIN, CDMAX, JBDAM, ILAT, CDPLTC, QINSUM, TINSUM, TIND)
    deallocate(QOLD, DTP, DTPS, QOLDS, QIND, JSS, HDIC, QNEW, YSS, VSS, YS, VS, VSTS, NSPRF)
    deallocate(LATGTC, LATSPC, LATPIC, DYNPIPE, LATPUC, DYNGTC, OPT, CIND, CINSUM, CDWBC, KFWBC, CPRWBC, CINBRC, CTRTRC, CDTBRC, DYNPUMP) ! SW 5/10/10
    deallocate(YSTS, YST, VST, ALLIM, APLIM, ANLIM, ASLIM, ELLIM, EPLIM, ENLIM, ESLIM, CSSK, C1, C2, Z0)
    deallocate(KFS, AF, EF, HYD, KF, AZSLC, STRIC, CPRBRC, CD, KBMAX, ELKT, WIND2, VISC, CELC, DLTADD)
    deallocate(QOAVR, QIMXR, QOMXR, REAERC, LAT, LONGIT, ELBOT, BTH, VPR, LPR, NISNP, NIPRF, NISPR, DECL)
    deallocate(A00, HH, T2I, KTWB, KBR, IBPR, DLVR, ESR, ETR, NBL, LPRFN, EXTFN, BTHFN, METFN)
    deallocate(SNPFN, PRFFN, SPRFN, CPLFN, VPLFN, FLXFN, FLXFN2, VPRFN, AFW, BFW, CFW, WINDH, RHEVC, FETCHC, JBDN, SPRVFN) ! SW 9/28/2018
    deallocate(KBI, MACCH, GRIDC, GMA, BTA, QTOT, SEDCIP, SEDCIN, SEDCIC, SEDCIS) ! SW 9/27/2007
    deallocate(SDK, FSOD, FSED, SEDCI, SEDCC, SEDPRC, ICEC, SLICEC, ICETHI, ALBEDO, HWI, BETAI, GAMMAI, ICEMIN)
    deallocate(SEDS, SEDB) !CB 11/28/06
    deallocate(EXH2O, BETA, EXOM, EXSS, DXI, CBHE, TSED, TSEDF, FI, ICET2, AZC, AZMAX) ! QINT,   QOUTT
    deallocate(AX, WTYPEC, TAIR, TDEW, WIND, PHI, CLOUD, CSHE, SRON, RANLW, RB, RC, RE, SHADE)
    deallocate(ET, RS, RN, SNPC, SCRC, PRFC, SPRC, CPLC, VPLC, FLXC, NXTMCP, NXTMVP, NXTMFL, GAMMA, F_NH3)
    deallocate(NXTMSN, NXTMSC, NXTMPR, NXTMSP, SNPDP, SCRDP, PRFDP, SPRDP, CPLDP, VPLDP, FLXDP, NCPL, NVPL, NFLX)
    deallocate(NSNP, NSCR, NPRF, NSPR, NEQN, PO4R, PARTP, NH4DK, NH4R, NO3DK, NO3S, CDSUM)
    deallocate(SSCS, CH4R, H2SR, FEIIR, MNIIR, SO4R)
    deallocate(BACTQ10, BACT1DK, BACTLDK, BACTS)
    deallocate(CoeffA_Turb, CoeffB_Turb, SECC_PAR)
    deallocate(H2SQ10, H2S1DK, CH4Q10, CH41DK)
    deallocate(KFE_OXID, KFE_RED, KFEOOH_HalfSat, FeSetVel)
    deallocate(KMN_OXID, KMN_RED, KMNO2_HalfSat, MnSetVel, DGPO2, MINKL)
    deallocate(LPOMHK, RPOMHK)
    deallocate(LDOMPDK, LRDOMPDK, RDOMPDK, LDOMNDK, LRDOMNDK, RDOMNDK, LDOMCDK, LRDOMCDK, RDOMCDK)
    deallocate(LPOMPDK, LRPOMPDK, RPOMPDK, LPOMNDK, LRPOMNDK, RPOMNDK, LPOMCDK, LRPOMCDK, RPOMCDK)
    deallocate(LDOMCMP, LPOMCMP, RPOMCMP)
    deallocate(LPZOOINC, LPZOOOUTC)
    deallocate(LDOP, RDOP, LPOP, RPOP, LDON, RDON, LPON, RPON)
    deallocate(LDOC, RDOC, LPOC, RPOC)
    deallocate(PSIEM, SEDEB)
    deallocate(LPOMPEP, LPOMNEP, LPOMCEP)
    deallocate(LPOMHD, RPOMHD)
    deallocate(LDOMCAP, LDOMCEP, LPOMCAP, LPOMCNS, RPOMCNS)
    deallocate(LDOMPD, LRDOMPD, RDOMPD, LPOMPD, LRPOMPD, RPOMPD, LPOMPHD, RPOMPHD)
    deallocate(LDOMND, LRDOMND, RDOMND, LPOMND, LRPOMND, RPOMND, LPOMNHD, RPOMNHD)
    deallocate(LDOMCD, LRDOMCD, RDOMCD, LPOMCD, LRPOMCD, RPOMCD, LPOMCHD, RPOMCHD)

    deallocate(CO2R, SROC, O2ER, O2EG, CAQ10, CADK, CAS, BODP, BODN, BODC, KBOD, TBOD, RBOD, DTRC)
    deallocate(LDOMDK, RDOMDK, LRDDK, OMT1, OMT2, OMK1, OMK2, LPOMDK, RPOMDK, LRPDK, POMS, ORGP, ORGN, ORGC)
    deallocate(RCOEF1, RCOEF2, RCOEF3, RCOEF4, ORGSI, NH4T1, NH4T2, NH4K1, NH4K2, NO3T1, NO3T2, NO3K1, NO3K2, NSTR)
    deallocate(DSIR, PSIS, PSIDK, PARTSI, SODT1, SODT2, SODK1, SODK2, O2NH4, O2OM, O2AR, O2AG, CG1DK, CGS)
    deallocate(CGQ10, CG0DK, CGLDK, CGKLF, CGCS, CGR, CUNIT, CUNIT1, CUNIT2, CAC, INCAC, TRCAC, DTCAC, PRCAC, CNAME, CNAME1, CNAME2, CMULT) !LCJ 2/26/15
    deallocate(CN, INCN, DTCN, PRCN, CSUM, DLTMAX, QWDO, TWDO, SSS, SEDRC, TAUCR, XBR, FNO3SED, DYNSTRUC)
!  DEALLOCATE (SSFLOC, FLOCEQN)                                                 
!DEALLOCATE (SEDCC1,SEDCC2, ICEQSS,SDK1,sdk2,SEDCI1,SEDCI2,SEDPRC1,SEDPRC2,SEDVP1,SEDVP2,SED1,SED2) 
!DEALLOCATE (SEDCC1,SEDCC2, ICEQSS,SDK1,sdk2,SEDCI1,SEDCI2,SEDPRC1,SEDPRC2,SEDVP1,SEDVP2,SED1,SED2,fsedc1,fsedc2,pbiom,nbiom,cbiom)   ! Amaila, cb 6/7/17
    deallocate(SEDCC1, SEDCC2, ICEQSS, SDK1, sdk2, SEDCI1, SEDCI2, SEDPRC1, SEDPRC2, SEDVP1, SEDVP2, SED1, SED2, fsedc1, fsedc2, pbiom, nbiom, cbiom, sed1ic, sed2ic, sdfirstadd) ! cb 9/3/17
    deallocate(ISO_SEDIMENT1, VERT_SEDIMENT1, LONG_SEDIMENT1, ISO_SEDIMENT2, VERT_SEDIMENT2, LONG_SEDIMENT2, PRINT_SEDIMENT1, PRINT_SEDIMENT2)
    deallocate(SEDIMENT_CALC1, SEDIMENT_CALC2)
    deallocate(QTAVB, QTMXB, BS, BE, JBUH, JBDH, TSSS, TSSB, TSSICE, ESBR, ETBR, EBRI, QDTR, EVBR)
    deallocate(QIN, PR, QPRBR, TIN, TOUT, TPR, TDTR, TPB, NACPR, NACIN, NACDT, NACTR, NACD, ELDH)
    deallocate(QSUM, NOUT, KTQIN, KBQIN, ELUH, NL, NPOINT, SLOPE, SLOPEC, ALPHA, COSA, SINA, SINAC, TDHFN, QOTFN, PREFN)
    deallocate(CPRFN, EUHFN, TUHFN, CUHFN, EDHFN, QINFN, TINFN, CINFN, CDHFN, QDTFN, TDTFN, CDTFN, TPRFN, VOLEV)
    deallocate(VOLWD, VOLSBR, VOLTBR, DLVOL, VOLG, VOLSR, VOLTR, VOLB, VOLPR, VOLTRB, VOLDT, VOLUH, VOLDH, VOLIN, VOLICE, ICEBANK)
    deallocate(US, DS, CUS, UHS, DHS, UQB, DQB, CDHS, VOLOUT, TSSWD, TSSUH, TSSDH, TSSIN, TSSOUT)
    deallocate(TSSEV, TSSPR, TSSTR, TSSDT, SOD, ELWS, BKT, REAER, ICETH, ICE, ICESW, Q, QC, QERR)
    deallocate(KTI, SROSH, SEG, DLXRHO, QSSUM, DLX, DLXR, QUH1, QDH1, BI, JWUH, JWDH)
    deallocate(A, C, D, F, V, SKTI, KBMIN, EV, QDT, QPR, SBKT, BHRHO)
    deallocate(SZ, WSHX, WSHY, WIND10, CZ, FETCH, PHI0, FRIC, ADZ, HMULT, FMTC, FMTCD, CNAME3, CDNAME3)
    deallocate(Z, KB, VNORM, ANPR, ANEQN, APOM, ACHLA, AHSP, AHSN, AHSSI)
    deallocate(AC, ASI, AT1, AT2, AT3, AT4, AK1, AK2, AK3, AK4, EXA, ASAT, AP, AN, AVERTM)
    deallocate(AG, AR, AE, AM, AS, ENPR, ENEQN, EG, ER, EE, EM, EB, ESAT, EP)
    deallocate(EC, ESI, ECHLA, EHSP, EHSN, EHSSI, EPOM, EHS, EN, ET4, EK1, EK2, EK3, EK4)
    deallocate(ET1, ET2, ET3, HNAME, FMTH, KFAC, KFNAME, KFNAME2, KFCN, C2I, TRCN, CDN, CDNAME, CDNAME2, CDMULT)
    deallocate(CMBRS, CMBRT, FETCHU, FETCHD, IPRF, ISNP, ISPR, BL, LFPR, DO3, SED, TKE, PALT)
    deallocate(ADX, DO1, DO2, B, CONV, CONV1, EL, DZ, DZQ, DX, SAZ, T1, TSS, QSS, BNEW, ILAYER) ! SW 1/23/06
    deallocate(P, SU, SW, BB, BR, BH, BHR, VOL, HSEG, DECAY, FPFE, FRICBR, UXBR, UYBR)
    deallocate(DEPTHB, DEPTHM, FPSS, TUH, TDH, TSSUH1, TSSUH2, TSSDH1, TSSDH2, SEDVP, H, EPC)
    deallocate(TVP, QINF, QOUT, KOUT, VOLUH2, VOLDH2, CWDO, CDWDO, CWDOC, CDWDOC, CDTOT, CPR, CPB, COUT)
    deallocate(CIN, CDTR, RSOD, RSOF, DLTD, DLTF, TSRD, TSRF, WDOD, WDOF, SNPD, SNPF, SPRD, SPRF)
    deallocate(SCRD, SCRF, PRFD, PRFF, CPLD, CPLF, VPLD, VPLF, FLXD, FLXF, EPIC, EPICI, EPIPRC, EPIVP)
    deallocate(CUH, CDH, EPM, EPD, C1S, CSSB, CVP, CSSUH1, CSSUH2, CSSDH2, CSSDH1, LNAME, IWR, KTWR, EKTWR, EKBWR)
    deallocate(JWUSP, JWDSP, QSP, KBWR, KTWD, KBWD, JBWD, GTA1, GTB1, GTA2, GTB2, BGT, IUGT, IDGT)
    deallocate(QTR, TTR, KTTR, KBTR, EGT, EGT2, AGASGT, BGASGT, CGASGT, GASGTC, PUGTC, ETUGT, EBUGT, KTUGT, KBUGT)
    deallocate(PDGTC, ETDGT, EBDGT, KTDGT, KBDGT, A1GT, B1GT, G1GT, A2GT, B2GT, G2GT, JWUGT, JWDGT, QGT)
    deallocate(EQGT, JBUGT, JBDGT, JBUPI, JBDPI, JWUPI, JWDPI, QPI, IUPI, IDPI, EUPI, EDPI, WPI, DLXPI, BP) ! SW 5/5/10
    deallocate(ETUPI, EBUPI, KTUPI, KBUPI, PDPIC, ETDPI, EBDPI, KTDPI, KBDPI, FPI, FMINPI, PUPIC, ETDSP, EBDSP)
    deallocate(PUSPC, ETUSP, EBUSP, KTUSP, KBUSP, PDSPC, KTDSP, KBDSP, IUSP, IDSP, ESP, A1SP, B1SP, A2SP)
    deallocate(B2SP, AGASSP, BGASSP, CGASSP, EQSP, GASSPC, JBUSP, JBDSP, STRTPU, ENDPU, EONPU, EOFFPU, QPU, PPUC)
    deallocate(IUPU, IDPU, EPU, ETPU, EBPU, KTPU, KBPU, JWUPU, JWDPU, JBUPU, JBDPU, PUMPON, KTW, KBW, PUMP_DOWNSTREAM)
    deallocate(IWD, KWD, QWD, EWD, ITR, QTRFN, TTRFN, CTRFN, ELTRT, ELTRB, TRC, JBTR, QTRF, CLRB)
    deallocate(TTLB, TTRB, CLLB, SRLB1, SRRB1, SRLB2, SRRB2, SRFJD1, SHADEI, SRFJD2, TOPO, QSW, CTR) ! SW 10/17/05
    deallocate(H1, H2, BH1, BH2, BHR1, BHR2, AVH1, AVH2, SAVH2, AVHR, SAVHR, CBODD)
    deallocate(POINT_SINK, HPRWBC, READ_EXTINCTION, READ_RADIATION)
    deallocate(DIST_TRIBS, UPWIND, ULTIMATE, FRESH_WATER, SALT_WATER, LIMITING_FACTOR)
    deallocate(UH_EXTERNAL, DH_EXTERNAL, UH_INTERNAL, DH_INTERNAL, UQ_INTERNAL, DQ_INTERNAL)
    deallocate(UQ_EXTERNAL, DQ_EXTERNAL, UP_FLOW, DN_FLOW, UP_HEAD, DN_HEAD)
    deallocate(INTERNAL_FLOW, DAM_INFLOW, DAM_OUTFLOW, HEAD_FLOW, HEAD_BOUNDARY) !TC 08/03/04
    deallocate(ISO_CONC, VERT_CONC, LONG_CONC, VERT_SEDIMENT, LONG_SEDIMENT)
    deallocate(ISO_SEDIMENT, VISCOSITY_LIMIT, CELERITY_LIMIT, IMPLICIT_AZ, ONE_LAYER, IMPLICIT_VISC)
    deallocate(FETCH_CALC, LIMITING_DLT, TERM_BY_TERM, MANNINGS_N, PLACE_QTR, SPECIFY_QTR)
    deallocate(PLACE_QIN, PRINT_CONST, PRINT_HYDRO, PRINT_SEDIMENT, ENERGY_BALANCE, MASS_BALANCE)
    deallocate(VOLUME_BALANCE, DETAILED_ICE, ICE_CALC, ALLOW_ICE, PH_CALC, BR_INACTIVE) ! ICE_IN,       RC/SW 4/28/11
    deallocate(BOD_CALCP, BOD_CALCN)
    deallocate(EVAPORATION, PRECIPITATION, RH_EVAP, NO_INFLOW, NO_OUTFLOW, NO_HEAT, BR_NOTECPLOT) ! SW 8/27/2019
    deallocate(ISO_TEMP, VERT_TEMP, LONG_TEMP, VERT_PROFILE, LONG_PROFILE, NO_WIND)
    deallocate(SNAPSHOT, PROFILE, VECTOR, CONTOUR, SPREADSHEET, INTERNAL_WEIR)
    deallocate(SCREEN_OUTPUT, FLUX, DYNAMIC_SHADE, TRAPEZOIDAL, BOD_CALC, ALG_CALC)
    deallocate(SEDIMENT_CALC, EPIPHYTON_CALC, PRINT_DERIVED, PRINT_EPIPHYTON, TDG_SPILLWAY, TDG_GATE, DYNSEDK)
    deallocate(ISO_EPIPHYTON, VERT_EPIPHYTON, LONG_EPIPHYTON, LATERAL_SPILLWAY, LATERAL_GATE, LATERAL_PUMP)
    deallocate(iso_macrophyte, vert_macrophyte, long_macrophyte, macrcvp, macrclp) ! cb 8/21/15
    deallocate(INTERP_HEAD, INTERP_WITHDRAWAL, INTERP_EXTINCTION, INTERP_DTRIBS, LATERAL_PIPE, INTERP_TRIBS)
    deallocate(INTERP_OUTFLOW, INTERP_INFLOW, INTERP_METEOROLOGY, ZERO_SLOPE)
    deallocate(SEDIMENT_RESUSPENSION, ACTIVE_RULE_W2SELECTIVE) !HYDRO_PLOT, CONSTITUENT_PLOT, DERIVED_PLOT,        
    deallocate(ORGPLD, ORGPRD, ORGPLP, ORGPRP, ORGNLD, ORGNRD, ORGNLP)
    deallocate(ICPL, TAVG, TAVGW, CAVG, CAVGW, CDAVG, CDAVGW)
    deallocate(ORGNRP, KG_H2O_CONSTANT)
    deallocate(PRINT_MACROPHYTE, MACROPHYTE_CALC, MACWBC, CONV2)
    deallocate(MAC, MACRC, MACT, MACRM, MACSS)
    deallocate(MGR, MMR, MRR)
    deallocate(SMACRC, SMACRM)
    deallocate(SMACT, SMAC)
    deallocate(MT1, MT2, MT3, MT4, MK1, MK2, MK3, MK4, MG, MR, MM)
    deallocate(MP, MN, MC, PSED, NSED, MHSP, MHSN, MHSC, MSAT)
    deallocate(CDDRAG, KTICOL, ARMAC, MACWBCI, ANORM, DWV, DWSA)
    deallocate(MBMP, MMAX, MPOM, LRPMAC, O2MR, O2MG)
    deallocate(MACMBRS, MACMBRT, SSMACMB)
    deallocate(CW, BIC)
    deallocate(MACTRMR, MACTRMF, MACTRM)
    deallocate(MLFPR)
    deallocate(MLLIM, MPLIM, MCLIM, MNLIM)
    deallocate(GAMMAJ)
    deallocate(POR, VOLKTI, VOLI, VSTEM, VSTEMKT, SAREA)
    deallocate(IWIND) ! MLM 08/12/05
    deallocate(ZG, ZM, ZEFF, PREFP, ZR, ZOOMIN, ZS2P, EXZ, ZT1, ZT2, ZT3, ZT4, ZK1, ZK2)
    deallocate(LDOMPMP, LDOMNMP, LPOMPMP, LPOMNMP, RPOMPMP, RPOMNMP, O2ZR) ! MLM 06/10/06
    deallocate(MPRWBC) ! MLM 06/10/06
    deallocate(EXM) ! MLM 06/10/06
    deallocate(USTARBTKE, E, EROUGH, ARODI, STRICK, TKELATPRDCONST, AZT, DZT)
    deallocate(FIRSTI, LASTI, TKELATPRD, STRICKON, WALLPNT, IMPTKE, TKEBC)
    deallocate(ZK3, ZK4, ZP, ZN, ZC, PREFA, ZMU, TGRAZE, ZRT, ZMT, ZOORM, ZOORMR, ZOORMF, ZSR, ZS) ! POINTERS ,ZOO,ZOOSS,   SW 1/28/2019
    deallocate(LPZOOOUT, LPZOOIN, PO4ZR, NH4ZR, DOZR, TICZR, AGZ, AGZT)
    deallocate(GTIC, BGTO, EGTO) ! CB 8/13/2010
    deallocate(INTERP_GATE) ! CB 8/13/2010  
    deallocate(ZGZ, PREFZ) !OMNIVOROUS ZOOPLANKTON
    deallocate(LPZOOINP, LPZOOINN, LPZOOOUTP, LPZOOOUTN)
    deallocate(SEDC, SEDN, SEDP, SEDNINFLUX, SEDPINFLUX, PFLUXIN, NFLUXIN)
    deallocate(SEDVPC, SEDVPP, SEDVPN)
    deallocate(SDKV, SEDDKTOT)
    deallocate(CBODS, KFJW)
    deallocate(BSAVE, GMA1, BTA1, ATMDEP_P, ATMDEP_N)
    deallocate(C_ATM_DEPOSITION, ATM_DEPOSITION, ATM_DEPOSITIONC, ATM_DEP_LOADING, ATM_DEPOSITION_INTERPOLATION, ATMDEPFN, ATMDCN, NACATD)
    deallocate(TN_SEDSOD_NH4, TP_SEDSOD_PO4, TPOUT, TPTRIB, TPDTRIB, TPWD, TPPR, TPIN, TNOUT, TNTRIB, TNDTRIB, TNWD, TNPR, TNIN, NH3GASLOSS) ! TP_SEDBURIAL,TN_SEDBURIAL,
    if (NBOD > 0) then
        deallocate(NBODC, NBODN, NBODP)
    end if
!IF(MWB_EXIST)THEN
    deallocate(WAIT_FOR_TRIB_INPUT, WAIT_FOR_BRANCH_INPUT, TR_FILEDIR, BR_FILEDIR) !SR 11/26/19
!ENDIF
    if (WAIT_FOR_INFLOW_RESULTS) then !SR 11/26/19
        deallocate(WAIT_TYPE, WAIT_INDEX, FILEDIR) !SR 11/26/19
        close(9911) !SR 11/26/19
    end if

    if (FISHBIO) then
        deallocate(BIOEXPFN, WEIGHTNUM, C2ZOO, VOLROOS, C2W, IBIO)
        deallocate(BIOD, BIOF, BIODP)
    end if

    call DEALLOCATE_TIME_VARYING_DATA()
    call DEALLOCATE_TRANSPORT()
    if (CONSTITUENTS) then
        call DEALLOCATE_KINETICS()
    end if
    call DEALLOCATE_WATERBODY()
    call DEALLOCATE_PIPE_FLOW()
    call DEALLOCATE_OPEN_CHANNEL()
    if (CONSTITUENTS .and. AERATEC == "      ON") then
        call DEALLOCATE_AERATE()
    end if
    if (SELECTC == "      ON") then
        call DEALLOCATE_SELECTIVE()
    end if
    if (SELECTC == "    USGS") then
        call DEALLOCATE_SELECTIVEUSGS()
    end if
    deallocate(WBSEG, PALT_JW, ELWS_INI, GTTYP, GTPC) ! systdg 
    if (SYSTDG) then
        call DEALLOCATE_SYSTDG()
    end if ! systdg
    if (TDGTA) then
        close(targetfnno)
    end if ! systdg TDGtarget
    if (TDGTA) then
        call DEALLOCATE_TDGtarget()
    end if ! systdg TDGtarget
    if (CONSTITUENTS .and. PHBUFF_EXIST) then
        deallocate(SDENI, PKI, PKSD)
        if (OM_BUFFERING) then
            if (OMTYPE == "    DIST") then
                deallocate(SDEN, PK, FRACT)
            else
                deallocate(SDEN, PK)
            end if
        end if
    end if
    if (CEMARelatedCode) then
        call DEALLOCATE_CEMA()
    end if

    return

end subroutine ENDSIMULATION
