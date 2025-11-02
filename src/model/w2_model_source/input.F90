subroutine INPUT()

    use MAIN
    use GLOBAL;     use NAMESC;     use GEOMC;     use LOGICC;     use PREC;     use SURFHE;     use KINETIC;     use SHADEC;     use EDDY
    use STRUCTURES;     use TRANS;     use TVDC;     use SELWC;     use GDAYC;     use SCREENC;     use TDGAS;     use RSTART
    use MACROPHYTEC;     use POROSITYC;     use ZOOPLANKTONC;     use INITIALVELOCITY;     use BIOENERGETICS;     use TRIDIAG_V;     use MSCLIB, only: RESTART_PUSHED
    use modSYSTDG ! SYSTDG
    use ALGAE_TOXINS
    implicit none
    external :: RESTART_OUTPUT

    real :: sum ! enhanced pH buffering
    integer :: NPROC, NNDC, N, NSTT, NIDUM, NDUM, JJ, NEPTT, NZPTT, NALT, NMCTT ! SW 7/13/09   9/28/2018
    character(len=1) :: CHAR1
    character(len=8) :: AID
    character(len=8) :: CDUM
    character(len=8) :: ORGCC
    integer, allocatable, dimension(:) :: IDUM
    real, allocatable, dimension(:) :: DDUM, FDUM

! Title and array dimensions

    allocate(TITLE(11))

    if (CONFN == "w2_con.npt") then
        read(CON, "(///(8X,A72))") (TITLE(J), J = 1, 10)
        read(CON, "(//8X,5I8,2A8)") NWB, NBR, IMX, KMX, NPROC, CLOSEC ! SW 7/31/09
        read(CON, "(//8X,8I8)") NTR, NST, NIW, NWD, NGT, NSP, NPI, NPU
        read(CON, "(//8X,7I8,a8)") NGC, NSS, NAL, NEP, NBOD, nmc, nzp
        read(CON, "(//8X,I8,7A8)") NOD, SELECTC, HABTATC, ENVIRPC, AERATEC, inituwl, ORGCC, SED_DIAG !, SYSTDGC, N2BNDC, DOBNDC, TDGTAC       ! systdg - Add control variables
    else
        read(CON, *)
        read(CON, *)
        read(CON, *)
        do J = 1, 10
            read(CON, *) TITLE(J)
        end do
        read(CON, *)
        read(CON, *)
        read(CON, *) NWB, NBR, IMX, KMX, NPROC, CLOSEC;         CLOSEC = ADJUSTR(CLOSEC) !'(A,5I0,A)'  
        read(CON, *)
        read(CON, *)
        read(CON, *) NTR, NST, NIW, NWD, NGT, NSP, NPI, NPU
        read(CON, *)
        read(CON, *)
        read(CON, *) NGC, NSS, NAL, NEP, NBOD, NMC, NZP
        read(CON, *)
        read(CON, *)
        read(CON, *) NOD, SELECTC, HABTATC, ENVIRPC, AERATEC, INITUWL, ORGCC, SED_DIAG !'(I0,5(A))'   
        SELECTC = ADJUSTR(SELECTC);         HABTATC = ADJUSTR(HABTATC);         ENVIRPC = ADJUSTR(ENVIRPC);         AERATEC = ADJUSTR(AERATEC);         INITUWL = ADJUSTR(INITUWL)
        ORGCC = ADJUSTR(ORGCC);         SED_DIAG = ADJUSTR(SED_DIAG)
    end if

    if (NPROC == 0) then
        NPROC = 1
    end if ! SW 7/31/09
!call omp_set_num_threads(NPROC)   ! set # of processors to NPROC  Moved to INPUT subroutine  TOGGLE FOR DEBUG
    if (SELECTC == "        ") then
        SELECTC = "     OFF"
    end if
!
    ORGC_CALC = ORGCC == "      ON"
! Fix the values here
    O2CH4 = 5.33
    O2H2S = 1.88
    O2FE2 = 0.143
    O2MN2 = 0.291

! Constituent numbers

    NTDS = 1
    NGCS = 2
    NGCE = NGCS + NGC - 1
    NSSS = NGCE + 1
    NSSE = NSSS + NSS - 1
    NWAGE = NSSE + 1
    NBACT = NWAGE + 1
    NDGP = NBACT + 1
    NN2 = NDGP + 1
    NH2S = NN2 + 1
    NCH4 = NH2S + 1
    NSO4 = NCH4 + 1
    NFEII = NSO4 + 1
    NFEOOH = NFEII + 1
    NMNII = NFEOOH + 1
    NMNO2 = NMNII + 1
    NPO4 = NMNO2 + 1
    NNH4 = NPO4 + 1
    NNO3 = NNH4 + 1
    NDSI = NNO3 + 1
    NPSI = NDSI + 1
    if (ORGC_CALC) then
        NLDOMC = NPSI + 1
        NRDOMC = NLDOMC + 1
        NLPOMC = NRDOMC + 1
        NRPOMC = NLPOMC + 1
        NBODS = NRPOMC + 1
    else
        NLDOM = NPSI + 1
        NRDOM = NLDOM + 1
        NLPOM = NRDOM + 1
        NRPOM = NLPOM + 1
        NBODS = NRPOM + 1
    end if
    if (NBOD > 0) then ! VARIABLE STOICHIOMETRY FOR CBOD    ! CB 6/6/10
        allocate(NBODC(NBOD), NBODP(NBOD), NBODN(NBOD))
        IBOD = NBODS
        NBODCS = IBOD
        do JCB = 1, NBOD
            NBODC(JCB) = IBOD
            IBOD = IBOD + 1
        end do
        NBODCE = IBOD - 1
        NBODPS = IBOD
        do JCB = 1, NBOD
            NBODP(JCB) = IBOD
            IBOD = IBOD + 1
        end do
        NBODPE = IBOD - 1
        NBODNS = IBOD
        do JCB = 1, NBOD
            NBODN(JCB) = IBOD
            IBOD = IBOD + 1
        end do
        NBODNE = IBOD - 1
    else
        NBODNS = 1;         NBODNE = 1;         NBODPS = 1;         NBODPE = 1;         NBODCS = 1;         NBODCE = 1

    end if
    NBODE = NBODS + NBOD*3 - 1 ! each BOD group has C, N and P groups
    NAS = NBODE + 1
    NAE = NAS + NAL - 1
    NDO = NAE + 1
    NTIC = NDO + 1
    NALK = NTIC + 1
    NZOOS = NALK + 1
    NZOOE = NZOOS + NZP - 1
    NLDOMP = NZOOE + 1
    NRDOMP = NLDOMP + 1
    NLPOMP = NRDOMP + 1
    NRPOMP = NLPOMP + 1
    NLDOMN = NRPOMP + 1
    NRDOMN = NLDOMN + 1
    NLPOMN = NRDOMN + 1
    NRPOMN = NLPOMN + 1
    NATS = NRPOMN + 1
    NATE = NATS + 3

    NCT = NATE !NRPOMN

! Constituent, tributary, and widthdrawal totals
    if (NTR == 0) then
        NTR1 = 1
    else
        NTR1 = NTR
    end if

    NTRT = NTR + NGT + NSP + NPI + NPU + NBR - 1 ! ADDING NBR FOR RESERVOIR FILLING    SW 6/12/2017
    NWDT = NWD + NGT + NSP + NPI + NPU
    NEPT = MAX(NEP, 1)
    NMCT = MAX(NMC, 1)
    NZPT = MAX(NZP, 1)

    if (NST > 5) then
        NSTT = NST ! FIXED # OF ROWS IN W2_CON.CSV MIN IS 5 INCREASE IF NECESSARY
    else
        NSTT = 5
    end if
    if (NEPT > 5) then
        NEPTT = NEPT
    else
        NEPTT = 5
    end if

    if (NMCT > 5) then
        NMCTT = NMCT
    else
        NMCTT = 5
    end if
    if (NAL > 5) then
        NALT = NAL
    else
        NALT = 5
    end if
    if (NZPT > 5) then
        NZPTT = NZPT
    else
        NZPTT = 5
    end if


    allocate(CDAC(NDC), X1(IMX), TECPLOT(NWB))
    allocate(BTA1(KMX), GMA1(KMX))
    allocate(WSC(IMX), KBI(IMX))
    allocate(VBC(NWB), EBC(NWB), MBC(NWB), PQC(NWB), EVC(NWB), PRC(NWB))
    allocate(WINDC(NWB), QINC(NWB), QOUTC(NWB), HEATC(NWB), SLHTC(NWB))
    allocate(QINIC(NBR), DTRIC(NBR), TRIC(NTR), WDIC(NWD), HDIC(NBR), METIC(NWB))
    allocate(EXC(NWB), EXIC(NWB))
    allocate(SLTRC(NWB), THETA(NWB), FRICC(NWB), NAF(NWB), ELTMF(NWB), Z0(NWB))
    allocate(ZMIN(NWB), IZMIN(NWB))
    allocate(C2CH(NCT), CDCH(NDC), EPCH(NEPT), macch(nmct), KFCH(NFL), APCH(NAL), ANCH(NAL), ALCH(NAL))
    allocate(CPLTC(NCT), HPLTC(NHY), CDPLTC(NDC))
    allocate(CMIN(NCT), CMAX(NCT), HYMIN(NHY), HYMAX(NHY), CDMIN(NDC), CDMAX(NDC))
    allocate(JBDAM(NBR), ILAT(NWDT))
    allocate(QINSUM(NBR), TINSUM(NBR), TIND(NBR), JSS(NBR), QIND(NBR))
    allocate(QOLD(NPI), DTP(NPI), DTPS(NPI), QOLDS(NPI))
    allocate(LATGTC(NGT), LATSPC(NSP), LATPIC(NPI), DYNPIPE(NPI), DYNPUMP(NPU), LATPUC(NPU), DYNGTC(NGT))
    allocate(GTIC(NGT), BGTO(NGT), EGTO(NGT)) ! cb 8/13/2010
    allocate(INTERP_GATE(NGT)) ! cb 8/13/2010  
    allocate(OPT(NWB, 8), CIND(NCT, NBR), CINSUM(NCT, NBR))
    allocate(CDWBC(NDC, NWB), KFWBC(NFL, NWB), CPRWBC(NCT, NWB), CINBRC(NCT, NBR), CTRTRC(NCT, NTR1))
    allocate(CDTBRC(NCT, NBR), CPRBRC(NCT, NBR))
    allocate(YSS(NNPIPE, NPI), VSS(NNPIPE, NPI), YS(NNPIPE, NPI), VS(NNPIPE, NPI), VSTS(NNPIPE, NPI))
    allocate(YSTS(NNPIPE, NPI), YST(NNPIPE, NPI), VST(NNPIPE, NPI))
    allocate(CBODD(KMX, IMX, NBOD))
    allocate(ALLIM(KMX, IMX, NAL), APLIM(KMX, IMX, NAL), ANLIM(KMX, IMX, NAL), ASLIM(KMX, IMX, NAL))
    allocate(ELLIM(KMX, IMX, NEP), EPLIM(KMX, IMX, NEP), ENLIM(KMX, IMX, NEP), ESLIM(KMX, IMX, NEP))
    allocate(CSSK(KMX, IMX, NCT), C1(KMX, IMX, NCT), C2(KMX, IMX, NCT), CD(KMX, IMX, NDC), KF(KMX, IMX, NFL))
    allocate(KFS(KMX, IMX, NFL), AF(KMX, IMX, NAL, 5), EF(KMX, IMX, NEP, 5), HYD(KMX, IMX, NHY), KFJW(NWB, NFL))
    allocate(TKE(KMX, IMX, 3), AZT(KMX, IMX), DZT(KMX, IMX))
    allocate(USTARBTKE(IMX), E(IMX), EROUGH(NWB), ARODI(NWB), STRICK(NWB), TKELATPRDCONST(NWB))
    allocate(FIRSTI(NWB), LASTI(NWB), TKELATPRD(NWB), STRICKON(NWB), WALLPNT(NWB), IMPTKE(NWB), TKEBC(NWB))
!  ALLOCATE (HYDRO_PLOT(NHY),    CONSTITUENT_PLOT(NCT), DERIVED_PLOT(NDC))
    allocate(ZERO_SLOPE(NWB), DYNAMIC_SHADE(IMX))
    allocate(AZSLC(NWB))
    allocate(NSPRF(NWB))
    allocate(KBMAX(NWB), ELKT(NWB), WIND2(IMX))
    allocate(VISC(NWB), CELC(NWB), DLTADD(NWB), REAERC(NWB))
    allocate(QOAVR(NWB), QIMXR(NWB), QOMXR(NWB))
    allocate(LAT(NWB), LONGIT(NWB), ELBOT(NWB))
    allocate(BTH(NWB), VPR(NWB), LPR(NWB))
    allocate(NISNP(NWB), NIPRF(NWB), NISPR(NWB))
    allocate(DDUM(NOD), FDUM(NOD), IDUM(IMX))
    allocate(ICPL(NWB))
    allocate(TN_SEDSOD_NH4(NWB), NH3GASLOSS(NWB), TP_SEDSOD_PO4(NWB), TPOUT(NWB), TPTRIB(NWB), TPDTRIB(NWB), TPWD(NWB), TPPR(NWB), TPIN(NWB), TNOUT(NWB), TNTRIB(NWB), TNDTRIB(NWB), TNWD(NWB), TNPR(NWB), TNIN(NWB)) ! TP_SEDBURIAL(NWB),TN_SEDBURIAL(NWB),
    allocate(A00(NWB), HH(NWB), DECL(NWB))
    allocate(T2I(NWB), KTWB(NWB), KBR(NWB), IBPR(NWB))
    allocate(DLVR(NWB), ESR(NWB), ETR(NWB), NBL(NWB))
    allocate(LPRFN(NWB), EXTFN(NWB), BTHFN(NWB), METFN(NWB), VPRFN(NWB))
    allocate(SNPFN(NWB), PRFFN(NWB), SPRFN(NWB), CPLFN(NWB), VPLFN(NWB), FLXFN(NWB), FLXFN2(NWB), SPRVFN(NWB)) ! SW 9/28/2018
    allocate(AFW(NWB), BFW(NWB), CFW(NWB), WINDH(NWB), RHEVC(NWB), FETCHC(NWB))
    allocate(SDK(NWB), FSOD(NWB), FSED(NWB), SEDCI(NWB), SEDCC(NWB), SEDPRC(NWB), SEDS(NWB), SEDB(NWB), DYNSEDK(NWB)) !cb 11/28/06
    allocate(SDK1(NWB), sdk2(nwb), SEDCI1(NWB), SEDCI2(NWB), SEDPRC1(NWB), SEDPRC2(NWB), SEDCC1(NWB), SEDCC2(NWB), fsedc1(nwb), fsedc2(nwb)) ! cb 6/17/17
    allocate(ICEC(NWB), SLICEC(NWB), ICETHI(NWB), ALBEDO(NWB), HWI(NWB), BETAI(NWB), GAMMAI(NWB), ICEMIN(NWB), ICET2(NWB))
    allocate(EXH2O(NWB), BETA(NWB), EXOM(NWB), EXSS(NWB), DXI(NWB), CBHE(NWB), TSED(NWB), TSEDF(NWB), FI(NWB))
    allocate(AX(NWB), WTYPEC(NWB), JBDN(NWB), AZC(NWB), AZMAX(NWB), GRIDC(NWB)) !SW 07/14/04    !  QINT(NWB),   QOUTT(NWB),  
    allocate(TAIR(NWB), TDEW(NWB), WIND(NWB), PHI(NWB), CLOUD(NWB), CSHE(IMX), SRON(NWB), RANLW(NWB))
    allocate(SNPC(NWB), SCRC(NWB), PRFC(NWB), SPRC(NWB), CPLC(NWB), VPLC(NWB), FLXC(NWB))
    allocate(NXTMSN(NWB), NXTMSC(NWB), NXTMPR(NWB), NXTMSP(NWB), NXTMCP(NWB), NXTMVP(NWB), NXTMFL(NWB))
    allocate(SNPDP(NWB), SCRDP(NWB), PRFDP(NWB), SPRDP(NWB), CPLDP(NWB), VPLDP(NWB), FLXDP(NWB))
    allocate(NSNP(NWB), NSCR(NWB), NPRF(NWB), NSPR(NWB), NCPL(NWB), NVPL(NWB), NFLX(NWB))
    allocate(NEQN(NWB), PO4R(NWB), PARTP(NWB))
    allocate(NH4DK(NWB), NH4R(NWB))
    allocate(CH4R(NWB), H2SR(NWB), FEIIR(NWB), MNIIR(NWB), SO4R(NWB))
    allocate(BACTQ10(NWB), BACT1DK(NWB), BACTLDK(NWB), BACTS(NWB))
!ALLOCATE (A_DISG(NWB),B_DISG(NWB),C_DISG(NWB))
    allocate(CoeffA_Turb(NWB), CoeffB_Turb(NWB), SECC_PAR(NWB))
    allocate(H2SQ10(NWB), H2S1DK(NWB), CH4Q10(NWB), CH41DK(NWB))
    allocate(KFE_OXID(NWB), KFE_RED(NWB), KFEOOH_HalfSat(NWB), FeSetVel(NWB))
    allocate(KMN_OXID(NWB), KMN_RED(NWB), KMNO2_HalfSat(NWB), MnSetVel(NWB))
    allocate(NO3DK(NWB), NO3S(NWB), FNO3SED(NWB))
    allocate(CO2R(NWB), SROC(NWB))
    allocate(O2ER(NEPT), O2EG(NEPT))
    allocate(CAQ10(NWB), CADK(NWB), CAS(NWB))
    allocate(BODP(NBOD), BODN(NBOD), BODC(NBOD))
    allocate(KBOD(NBOD), TBOD(NBOD), RBOD(NBOD))
    allocate(LDOMDK(NWB), RDOMDK(NWB), LRDDK(NWB))
    allocate(OMT1(NWB), OMT2(NWB), OMK1(NWB), OMK2(NWB))
    allocate(LPOMDK(NWB), RPOMDK(NWB), LRPDK(NWB), POMS(NWB))
    allocate(ORGP(NWB), ORGN(NWB), ORGC(NWB), ORGSI(NWB))
    allocate(Pbiom(NWB), Nbiom(NWB), Cbiom(NWB)) ! Amaila, cb 6/8/17
    allocate(RCOEF1(NWB), RCOEF2(NWB), RCOEF3(NWB), RCOEF4(NWB), DGPO2(NWB), MINKL(NWB))
    allocate(NH4T1(NWB), NH4T2(NWB), NH4K1(NWB), NH4K2(NWB), KG_H2O_CONSTANT(NWB))
    allocate(NO3T1(NWB), NO3T2(NWB), NO3K1(NWB), NO3K2(NWB))
    allocate(DSIR(NWB), PSIS(NWB), PSIDK(NWB), PARTSI(NWB))
    allocate(SODT1(NWB), SODT2(NWB), SODK1(NWB), SODK2(NWB))
    allocate(O2NH4(NWB), O2OM(NWB))
    allocate(O2AR(NAL), O2AG(NAL))
    allocate(CGQ10(NGC), CG0DK(NGC), CG1DK(NGC), CGS(NGC), CGLDK(NGC), CGKLF(NGC), CGCS(NGC), CGR(NGC)) !LCJ 2/26/15
    allocate(CUNIT(NCT), CUNIT1(NCT), CUNIT2(NCT))
    allocate(CAC(NCT), INCAC(NCT), TRCAC(NCT), DTCAC(NCT), PRCAC(NCT))
    allocate(CNAME(NCT), CNAME1(NCT), CNAME2(NCT), CNAME3(NCT), CMULT(NCT), CSUM(NCT))
    allocate(CN(NCT))
    allocate(SSS(NSS), TAUCR(NSS), SEDRC(NSS), SSCS(NSS)) !,  SSFLOC(NSS), FLOCEQN(NSS))                                           !SR 04/21/13
    allocate(CDSUM(NDC))
    allocate(DTRC(NBR))
    allocate(NSTR(NBR), XBR(NBR), DYNSTRUC(NBR))
    allocate(QTAVB(NBR), QTMXB(NBR))
    allocate(BS(NWB), BE(NWB), JBUH(NBR), JBDH(NBR), JWUH(NBR), JWDH(NBR))
    allocate(TSSS(NBR), TSSB(NBR), TSSICE(NBR))
    allocate(ESBR(NBR), ETBR(NBR), EBRI(NBR))
    allocate(QIN(NBR), PR(NBR), QPRBR(NBR), QDTR(NBR), EVBR(NBR))
    allocate(TIN(NBR), TOUT(NBR), TPR(NBR), TDTR(NBR), TPB(NBR))
    allocate(NACPR(NBR), NACIN(NBR), NACDT(NBR), NACTR(NTR), NACD(NWB))
    allocate(QSUM(NBR), NOUT(NBR), KTQIN(NBR), KBQIN(NBR), ELUH(NBR), ELDH(NBR))
    allocate(NL(NBR), NPOINT(NBR), SLOPE(NBR), SLOPEC(NBR), ALPHA(NBR), COSA(NBR), SINA(NBR), SINAC(NBR), ilayer(imx))
    allocate(CPRFN(NBR), EUHFN(NBR), TUHFN(NBR), CUHFN(NBR), EDHFN(NBR), TDHFN(NBR), QOTFN(NBR), PREFN(NBR))
    allocate(QINFN(NBR), TINFN(NBR), CINFN(NBR), CDHFN(NBR), QDTFN(NBR), TDTFN(NBR), CDTFN(NBR), TPRFN(NBR))
    allocate(VOLWD(NBR), VOLSBR(NBR), VOLTBR(NBR), DLVOL(NBR), VOLG(NWB), VOLSR(NWB), VOLTR(NWB), VOLEV(NBR), VOLICE(NBR), ICEBANK(IMX))
    allocate(VOLB(NBR), VOLPR(NBR), VOLTRB(NBR), VOLDT(NBR), VOLUH(NBR), VOLDH(NBR), VOLIN(NBR), VOLOUT(NBR))
    allocate(US(NBR), DS(NBR), CUS(NBR), UHS(NBR), DHS(NBR), UQB(NBR), DQB(NBR), CDHS(NBR))
    allocate(TSSEV(NBR), TSSPR(NBR), TSSTR(NBR), TSSDT(NBR), TSSWD(NBR), TSSUH(NBR), TSSDH(NBR), TSSIN(NBR), TSSOUT(NBR))
    allocate(ET(IMX), RS(IMX), RN(IMX), RB(IMX), RC(IMX), RE(IMX), SHADE(IMX))
    allocate(DLTMAX(NOD), QWDO(IMX), TWDO(IMX)) ! SW 1/24/05
    allocate(SOD(IMX), ELWS(IMX), BKT(IMX), REAER(IMX))
    allocate(ICETH(IMX), ICE(IMX), ICESW(IMX))
    allocate(Q(IMX), QC(IMX), QERR(IMX), QSSUM(IMX))
    allocate(KTI(IMX), SKTI(IMX), SROSH(IMX), SEG(IMX), DLXRHO(IMX))
    allocate(DLX(IMX), DLXR(IMX))
    allocate(A(IMX), C(IMX), D(IMX), F(IMX), V(IMX), BTA(IMX), GMA(IMX))
    allocate(KBMIN(IMX), EV(IMX), QDT(IMX), QPR(IMX), SBKT(IMX), BHRHO(IMX))
    allocate(SZ(IMX), WSHX(IMX), WSHY(IMX), WIND10(IMX), CZ(IMX), FETCH(IMX), PHI0(IMX), FRIC(IMX))
    allocate(Z(IMX), KB(IMX), PALT(IMX))
    allocate(VNORM(KMX))
    allocate(ANPR(NAL), ANEQN(NAL), APOM(NAL))
    allocate(AC(NAL), ASI(NAL), ACHLA(NAL), AHSP(NAL), AHSN(NAL), AHSSI(NAL))
    allocate(AT1(NAL), AT2(NAL), AT3(NAL), AT4(NAL), AK1(NAL), AK2(NAL), AK3(NAL), AK4(NAL), AVERTM(NAL))
    allocate(AG(NAL), AR(NAL), AE(NAL), AM(NAL), AS(NAL), EXA(NAL), ASAT(NAL), AP(NAL), AN(NAL))
    allocate(ENPR(NEPT), ENEQN(NEPT))
    allocate(EG(NEPT), ER(NEPT), EE(NEPT), EM(NEPT), EB(NEPT), ESAT(NEPT), EP(NEPT), EN(NEPT))
    allocate(EC(NEPT), ESI(NEPT), ECHLA(NEPT), EHSP(NEPT), EHSN(NEPT), EHSSI(NEPT), EPOM(NEPT), EHS(NEPT))
    allocate(ET1(NEPT), ET2(NEPT), ET3(NEPT), ET4(NEPT), EK1(NEPT), EK2(NEPT), EK3(NEPT), EK4(NEPT))
    allocate(HNAME(NHY), FMTH(NHY), HMULT(NHY), FMTC(NCT), FMTCD(NDC))
    allocate(KFAC(NFL), KFNAME(NFL), KFNAME2(NFL), KFCN(NFL, NWB))
    allocate(ATMDCN(NCT, NWB), NACATD(NWB))
    allocate(C2I(NCT, NWB), TRCN(NCT, NTR), C_ATM_DEPOSITION(NCT, NWB), ATM_DEPOSITIONC(NWB), ATM_DEPOSITION(NWB), ATM_DEP_LOADING(NCT, NWB), ATM_DEPOSITION_INTERPOLATION(NWB), ATMDEPFN(NWB), ATMDEP_P(NWB), ATMDEP_N(NWB))
    allocate(CDN(NDC, NWB), CDNAME(NDC), CDNAME2(NDC), CDNAME3(NDC), CDMULT(NDC))
    allocate(CMBRS(NCT, NBR), CMBRT(NCT, NBR), INCN(NCT, NBR), DTCN(NCT, NBR), PRCN(NCT, NBR))
    allocate(FETCHU(IMX, NBR), FETCHD(IMX, NBR))
    allocate(IPRF(IMX, NWB), ISNP(IMX, NWB), ISPR(IMX, NWB), BL(IMX, NWB))
    allocate(H1(KMX, IMX), H2(KMX, IMX), BH1(KMX, IMX), BH2(KMX, IMX), BHR1(KMX, IMX), BHR2(KMX, IMX), QTOT(KMX, IMX))
    allocate(SAVH2(KMX, IMX), AVH1(KMX, IMX), AVH2(KMX, IMX), AVHR(KMX, IMX), SAVHR(KMX, IMX))
    allocate(LFPR(KMX, IMX), BI(KMX, IMX), BNEW(KMX, IMX)) ! SW 1/23/06
    allocate(ADX(KMX, IMX), ADZ(KMX, IMX), DO1(KMX, IMX), DO2(KMX, IMX), DO3(KMX, IMX), SED(KMX, IMX))
    allocate(B(KMX, IMX), CONV(KMX, IMX), CONV1(KMX, IMX), EL(KMX, IMX), DZ(KMX, IMX), DZQ(KMX, IMX), DX(KMX, IMX))
    allocate(P(KMX, IMX), SU(KMX, IMX), SW(KMX, IMX), SAZ(KMX, IMX), T1(KMX, IMX), TSS(KMX, IMX), QSS(KMX, IMX))
    allocate(BB(KMX, IMX), BR(KMX, IMX), BH(KMX, IMX), BHR(KMX, IMX), VOL(KMX, IMX), HSEG(KMX, IMX), DECAY(KMX, IMX), CONSTRICTION(KMX, IMX), BCONSTRICTION(IMX)) ! SW 6/26/2018
    allocate(DEPTHB(KMX, IMX), DEPTHM(KMX, IMX), FPSS(KMX, IMX), FPFE(KMX, IMX), FRICBR(KMX, IMX), UXBR(KMX, IMX), UYBR(KMX, IMX))
    allocate(QUH1(KMX, NBR), QDH1(KMX, NBR), VOLUH2(KMX, NBR), VOLDH2(KMX, NBR), TUH(KMX, NBR), TDH(KMX, NBR))
    allocate(TSSUH1(KMX, NBR), TSSUH2(KMX, NBR), TSSDH1(KMX, NBR), TSSDH2(KMX, NBR))
    allocate(TVP(KMX, NWB), SEDVP(KMX, NWB), H(KMX, NWB))
    allocate(SEDVP1(KMX, NWB), SEDVP2(KMX, NWB), SED1(KMX, IMX), SED2(KMX, IMX), sed1ic(kmx, imx), sed2ic(kmx, imx), sdfirstadd(kmx, imx)) !  cb 9/3/17
    allocate(QINF(KMX, NBR), QOUT(KMX, NBR), KOUT(KMX, NBR))
    allocate(CT(KMX, IMX), AT(KMX, IMX), VT(KMX, IMX), DT(KMX, IMX), GAMMA(KMX, IMX), F_NH3(KMX, IMX))
    allocate(CWDO(NCT, NOD), CDWDO(NDC, NOD), CWDOC(NCT), CDWDOC(NDC), CDTOT(NDC))
    allocate(CIN(NCT, NBR), CDTR(NCT, NBR), CPR(NCT, NBR), CPB(NCT, NBR), COUT(NCT, NBR))
    allocate(RSOD(NOD), RSOF(NOD), DLTD(NOD), DLTF(NOD))
    allocate(TSRD(NOD), TSRF(NOD), WDOD(NOD), WDOF(NOD))
    allocate(SNPD(NOD, NWB), SNPF(NOD, NWB), SPRD(NOD, NWB), SPRF(NOD, NWB))
    allocate(SCRD(NOD, NWB), SCRF(NOD, NWB), PRFD(NOD, NWB), PRFF(NOD, NWB))
    allocate(CPLD(NOD, NWB), CPLF(NOD, NWB), VPLD(NOD, NWB), VPLF(NOD, NWB), FLXD(NOD, NWB), FLXF(NOD, NWB))
    allocate(EPIC(NWB, NEPTT), EPICI(NWB, NEPTT), EPIPRC(NWB, NEPTT))
    allocate(EPIVP(KMX, NWB, NEP), macrcvp(KMX, NWB, nmc), macrclp(KMX, imx, nmc)) ! cb 8/21/15
    allocate(CUH(KMX, NCT, NBR), CDH(KMX, NCT, NBR))
    allocate(EPM(KMX, IMX, NEPT), EPD(KMX, IMX, NEPT), EPC(KMX, IMX, NEPT))
    allocate(C1S(KMX, IMX, NCT), CSSB(KMX, IMX, NCT), CVP(KMX, NCT, NWB))
    allocate(CSSUH1(KMX, NCT, NBR), CSSUH2(KMX, NCT, NBR), CSSDH2(KMX, NCT, NBR), CSSDH1(KMX, NCT, NBR))
    allocate(READ_EXTINCTION(NWB), READ_RADIATION(NWB))
    allocate(DIST_TRIBS(NBR), LIMITING_FACTOR(NAL))
    allocate(UPWIND(NWB), ULTIMATE(NWB))
    allocate(STRIC(NSTT, NBR), ESTRT(NSTT, NBR), WSTRT(NSTT, NBR), KTSWT(NSTT, NBR), KBSWT(NSTT, NBR), SINKCT(NSTT, NBR))
    allocate(FRESH_WATER(NWB), SALT_WATER(NWB), TRAPEZOIDAL(NWB)) !SW 07/16/04
    allocate(UH_EXTERNAL(NBR), DH_EXTERNAL(NBR), UH_INTERNAL(NBR), DH_INTERNAL(NBR))
    allocate(UQ_EXTERNAL(NBR), DQ_EXTERNAL(NBR), UQ_INTERNAL(NBR), DQ_INTERNAL(NBR))
    allocate(UP_FLOW(NBR), DN_FLOW(NBR), UP_HEAD(NBR), DN_HEAD(NBR))
    allocate(INTERNAL_FLOW(NBR), DAM_INFLOW(NBR), DAM_OUTFLOW(NBR), HEAD_FLOW(NBR), HEAD_BOUNDARY(NWB)) !TC 08/03/04
    allocate(ISO_CONC(NCT, NWB), VERT_CONC(NCT, NWB), LONG_CONC(NCT, NWB))
    allocate(ISO_SEDIMENT(NWB), VERT_SEDIMENT(NWB), LONG_SEDIMENT(NWB))
    allocate(ISO_SEDIMENT1(NWB), VERT_SEDIMENT1(NWB), LONG_SEDIMENT1(NWB)) !Amaila
    allocate(ISO_SEDIMENT2(NWB), VERT_SEDIMENT2(NWB), LONG_SEDIMENT2(NWB)) !Amaila
    allocate(VISCOSITY_LIMIT(NWB), CELERITY_LIMIT(NWB), IMPLICIT_AZ(NWB))
    allocate(FETCH_CALC(NWB), ONE_LAYER(IMX), IMPLICIT_VISC(NWB))
    allocate(LIMITING_DLT(NWB), TERM_BY_TERM(NWB), MANNINGS_N(NWB))
    allocate(PLACE_QIN(NWB), PLACE_QTR(NTRT), SPECIFY_QTR(NTRT))
    allocate(PRINT_CONST(NCT, NWB), PRINT_HYDRO(NHY, NWB), PRINT_SEDIMENT(NWB))
    allocate(PRINT_SEDIMENT1(NWB), PRINT_SEDIMENT2(NWB)) ! Amaila
    allocate(VOLUME_BALANCE(NWB), ENERGY_BALANCE(NWB), MASS_BALANCE(NWB))
    allocate(DETAILED_ICE(NWB), ICE_CALC(NWB), ALLOW_ICE(IMX), BR_INACTIVE(NBR), BR_NOTECPLOT(NBR)) !   ICE_IN(NBR),    RC/SW 4/28/11 SW 8/27/2019
    allocate(EVAPORATION(NWB), PRECIPITATION(NWB), RH_EVAP(NWB), PH_CALC(NWB))
    allocate(NO_INFLOW(NWB), NO_OUTFLOW(NWB), NO_HEAT(NWB), NO_WIND(NWB))
    allocate(ISO_TEMP(NWB), VERT_TEMP(NWB), LONG_TEMP(NWB), VERT_PROFILE(NWB), LONG_PROFILE(NWB))
    allocate(SNAPSHOT(NWB), PROFILE(NWB), VECTOR(NWB), CONTOUR(NWB), SPREADSHEET(NWB))
    allocate(SCREEN_OUTPUT(NWB), FLUX(NWB))
    allocate(PRINT_DERIVED(NDC, NWB), PRINT_EPIPHYTON(NWB, NEPT))
    allocate(SEDIMENT_CALC(NWB), EPIPHYTON_CALC(NWB, NEPT), SEDIMENT_RESUSPENSION(NSS), BOD_CALC(NBOD), ALG_CALC(NAL))
    allocate(SEDIMENT_CALC1(NWB), SEDIMENT_CALC2(NWB)) ! Amaila
    allocate(BOD_CALCP(NBOD), BOD_CALCN(NBOD)) ! cb 5/19/2011
    allocate(TDG_SPILLWAY(NWDT, NSP), TDG_GATE(NWDT, NGT), INTERNAL_WEIR(KMX, IMX))
    allocate(ISO_EPIPHYTON(NWB, NEPT), VERT_EPIPHYTON(NWB, NEPT), LONG_EPIPHYTON(NWB, NEPT))
    allocate(iso_macrophyte(NWB, nmc), vert_macrophyte(NWB, nmc), long_macrophyte(NWB, nmc)) ! cb 8/21/15
    allocate(LATERAL_SPILLWAY(NSP), LATERAL_GATE(NGT), LATERAL_PUMP(NPU), LATERAL_PIPE(NPI))
    allocate(INTERP_HEAD(NBR), INTERP_WITHDRAWAL(NWD), INTERP_EXTINCTION(NWB), INTERP_DTRIBS(NBR))
    allocate(INTERP_OUTFLOW(NST, NBR), INTERP_INFLOW(NBR), INTERP_METEOROLOGY(NWB), INTERP_TRIBS(NTR))
    allocate(LNAME(NCT + NHY + NDC))
    allocate(IWR(NIW), KTWR(NIW), KBWR(NIW), EKTWR(NIW), EKBWR(NIW)) ! SW 3/18/16
    allocate(JWUSP(NSP), JWDSP(NSP), QSP(NSP))
    allocate(KTWD(NWDT), KBWD(NWDT), JBWD(NWDT))
    allocate(GTA1(NGT), GTB1(NGT), GTA2(NGT), GTB2(NGT))
    allocate(BGT(NGT), IUGT(NGT), IDGT(NGT), EGT(NGT), EGT2(NGT))
    allocate(QTR(NTRT), TTR(NTRT), KTTR(NTRT), KBTR(NTRT))
    allocate(AGASGT(NGT), BGASGT(NGT), CGASGT(NGT), GASGTC(NGT))
    allocate(PUGTC(NGT), ETUGT(NGT), EBUGT(NGT), KTUGT(NGT), KBUGT(NGT))
    allocate(PDGTC(NGT), ETDGT(NGT), EBDGT(NGT), KTDGT(NGT), KBDGT(NGT))
    allocate(A1GT(NGT), B1GT(NGT), G1GT(NGT), A2GT(NGT), B2GT(NGT), G2GT(NGT))
    allocate(EQGT(NGT), JBUGT(NGT), JBDGT(NGT), JWUGT(NGT), JWDGT(NGT), QGT(NGT))
    allocate(JBUPI(NPI), JBDPI(NPI), JWUPI(NPI), JWDPI(NPI), QPI(NPI), BP(NPI)) ! SW 5/10/10
    allocate(IUPI(NPI), IDPI(NPI), EUPI(NPI), EDPI(NPI), WPI(NPI), DLXPI(NPI), FPI(NPI), FMINPI(NPI), PUPIC(NPI))
    allocate(ETUPI(NPI), EBUPI(NPI), KTUPI(NPI), KBUPI(NPI), PDPIC(NPI), ETDPI(NPI), EBDPI(NPI), KTDPI(NPI), KBDPI(NPI))
    allocate(PUSPC(NSP), ETUSP(NSP), EBUSP(NSP), KTUSP(NSP), KBUSP(NSP), PDSPC(NSP), ETDSP(NSP), EBDSP(NSP))
    allocate(KTDSP(NSP), KBDSP(NSP), IUSP(NSP), IDSP(NSP), ESP(NSP), A1SP(NSP), B1SP(NSP), A2SP(NSP))
    allocate(B2SP(NSP), AGASSP(NSP), BGASSP(NSP), CGASSP(NSP), EQSP(NSP), GASSPC(NSP), JBUSP(NSP), JBDSP(NSP))
    allocate(IUPU(NPU), IDPU(NPU), EPU(NPU), STRTPU(NPU), ENDPU(NPU), EONPU(NPU), EOFFPU(NPU), QPU(NPU), PPUC(NPU))
    allocate(ETPU(NPU), EBPU(NPU), KTPU(NPU), KBPU(NPU), JWUPU(NPU), JWDPU(NPU), JBUPU(NPU), JBDPU(NPU), PUMPON(NPU), PUMP_DOWNSTREAM(NPU))
    allocate(IWD(NWDT), KWD(NWDT), QWD(NWDT), EWD(NWDT), KTW(NWDT), KBW(NWDT))
    allocate(ITR(NTRT), QTRFN(NTR), TTRFN(NTR), CTRFN(NTR), ELTRT(NTRT), ELTRB(NTRT), TRC(NTRT), JBTR(NTRT), QTRF(KMX, NTRT))
    allocate(TTLB(IMX), TTRB(IMX), CLLB(IMX), CLRB(IMX))
    allocate(SRLB1(IMX), SRRB1(IMX), SRLB2(IMX), SRRB2(IMX), SRFJD1(IMX), SHADEI(IMX), SRFJD2(IMX))
    allocate(TOPO(IMX, IANG)) ! SW 10/17/05
    allocate(QSW(KMX, NWDT), CTR(NCT, NTRT), HPRWBC(NHY, NWB))
    allocate(RATZ(KMX, NWB), CURZ1(KMX, NWB), CURZ2(KMX, NWB), CURZ3(KMX, NWB)) ! SW 5/15/06
    allocate(ZG(NZPT), ZM(NZPT), ZEFF(NZPT), PREFP(NZPT), ZR(NZPT), ZOOMIN(NZPT), ZS2P(NZPT), EXZ(NZPT), PREFZ(NZPTT, NZPT), ZS(NZPT)) ! SW 1/29/2018
    allocate(ZT1(NZPT), ZT2(NZPT), ZT3(NZPT), ZT4(NZPT), ZK1(NZPT), ZK2(NZPT), ZK3(NZPT), ZK4(NZPT), O2ZR(NZPT))
    allocate(ZP(NZPT), ZN(NZPT), ZC(NZPT))
    allocate(PREFA(NALT, NZPT))
    allocate(PO4ZR(KMX, IMX), NH4ZR(KMX, IMX))
    allocate(ZMU(KMX, IMX, NZP), TGRAZE(KMX, IMX, NZP), ZRT(KMX, IMX, NZP), ZMT(KMX, IMX, NZP), ZSR(KMX, IMX, NZP)) ! MLM POINTERS:,ZOO(KMX,IMX,NZP),ZOOSS(KMX,IMX,NZP))   SW 1/29/2019
    allocate(ZOORM(KMX, IMX, NZP), ZOORMR(KMX, IMX, NZP), ZOORMF(KMX, IMX, NZP))
    allocate(LPZOOOUT(KMX, IMX), LPZOOIN(KMX, IMX), DOZR(KMX, IMX), TICZR(KMX, IMX))
    allocate(AGZ(KMX, IMX, NAL, NZP), ZGZ(KMX, IMX, NZP, NZP), AGZT(KMX, IMX, NAL)) !OMNIVOROUS ZOOPLANKTON
    allocate(ORGPLD(KMX, IMX), ORGPRD(KMX, IMX), ORGPLP(KMX, IMX), ORGPRP(KMX, IMX), ORGNLD(KMX, IMX), ORGNRD(KMX, IMX), ORGNLP(KMX, IMX))
    allocate(ORGNRP(KMX, IMX))
    allocate(LDOMPMP(KMX, IMX), LDOMNMP(KMX, IMX), LPOMPMP(KMX, IMX), LPOMNMP(KMX, IMX), RPOMPMP(KMX, IMX), RPOMNMP(KMX, IMX))
    allocate(LPZOOINP(KMX, IMX), LPZOOINN(KMX, IMX), LPZOOOUTP(KMX, IMX), LPZOOOUTN(KMX, IMX))
    allocate(SEDVPP(KMX, NWB), SEDVPC(KMX, NWB), SEDVPN(KMX, NWB))
    allocate(SEDP(KMX, IMX), SEDN(KMX, IMX), SEDC(KMX, IMX), SEDNINFLUX(KMX, IMX), SEDPINFLUX(KMX, IMX), PFLUXIN(NWB), NFLUXIN(NWB))
    allocate(SDKV(KMX, IMX), SEDDKTOT(KMX, IMX))
    allocate(SEDCIP(NWB), SEDCIN(NWB), SEDCIC(NWB), SEDCIS(NWB))
    allocate(CBODS(NBOD)) !, CBODNS(KMX,IMX))     POINTERS SW 3/2019 , SEDCB(KMX,IMX), SEDCBP(KMX,IMX), SEDCBN(KMX,IMX), SEDCBC(KMX,IMX))
    allocate(PRINT_MACROPHYTE(NWB, NMCT), MACROPHYTE_CALC(NWB, NMCT), MACWBC(NWB, NMCTT), CONV2(KMX, KMX), MPRWBC(NWB, NMCTT))
    allocate(MAC(KMX, IMX, NMCT), MACRC(KMX, KMX, IMX, NMCT), MACT(KMX, KMX, IMX), MACRM(KMX, KMX, IMX, NMCT), MACSS(KMX, KMX, IMX, NMCT))
    allocate(MGR(KMX, KMX, IMX, NMCT), MMR(KMX, IMX, NMCT), MRR(KMX, IMX, NMCT))
    allocate(SMACRC(KMX, KMX, IMX, NMCT), SMACRM(KMX, KMX, IMX, NMCT))
    allocate(SMACT(KMX, KMX, IMX), SMAC(KMX, IMX, NMCT))
    allocate(MT1(NMCT), MT2(NMCT), MT3(NMCT), MT4(NMCT), MK1(NMCT), MK2(NMCT), MK3(NMCT), MK4(NMCT), MG(NMCT), MR(NMCT), MM(NMCT))
    allocate(MBMP(NMCT), MMAX(NMCT), CDDRAG(NMCT), DWV(NMCT), DWSA(NMCT), ANORM(NMCT))
    allocate(MP(NMCT), MN(NMCT), MC(NMCT), PSED(NMCT), NSED(NMCT), MHSP(NMCT), MHSN(NMCT), MHSC(NMCT), MSAT(NMCT), EXM(NMCT))
    allocate(O2MG(NMCT), O2MR(NMCT), LRPMAC(NMCT), MPOM(NMCT))
    allocate(KTICOL(IMX), ARMAC(IMX), MACWBCI(NWB, NMCTT))
    allocate(MACMBRS(NBR, NMCT), MACMBRT(NBR, NMCT), SSMACMB(NBR, NMCT))
    allocate(CW(KMX, IMX), BIC(KMX, IMX))
    allocate(MACTRMR(KMX, IMX, NMCT), MACTRMF(KMX, IMX, NMCT), MACTRM(KMX, IMX, NMCT))
    allocate(MLFPR(KMX, KMX, IMX, NMCT))
    allocate(MLLIM(KMX, KMX, IMX, NMCT), MPLIM(KMX, IMX, NMCT), MCLIM(KMX, IMX, NMCT), MNLIM(KMX, IMX, NMCT))
    allocate(GAMMAJ(KMX, KMX, IMX))
    allocate(POR(KMX, IMX), VOLKTI(IMX), VOLI(KMX, IMX), VSTEM(KMX, IMX, NMCT), VSTEMKT(IMX, NMCT), SAREA(NMCT))
    allocate(IWIND(NWB))
    allocate(LAYERCHANGE(NWB))
!ALLOCATE (CBODP(KMX,IMX,NBOD), CBODN(KMX,IMX,NBOD))      ! CB 6/6/10  ! SW 3/2019 CODE ERROR MENORY LEAK SINCE POINTER NO NEED TO ALLOCATE
    allocate(HAB(KMX, IMX))
!
!
    allocate(LPOMHK(NWB), RPOMHK(NWB))
    allocate(LDOMPDK(NWB), LRDOMPDK(NWB), RDOMPDK(NWB), LDOMNDK(NWB), LRDOMNDK(NWB), RDOMNDK(NWB), LDOMCDK(NWB), LRDOMCDK(NWB), RDOMCDK(NWB))
    allocate(LPOMPDK(NWB), LRPOMPDK(NWB), RPOMPDK(NWB), LPOMNDK(NWB), LRPOMNDK(NWB), RPOMNDK(NWB), LPOMCDK(NWB), LRPOMCDK(NWB), RPOMCDK(NWB))
    allocate(LDOMCMP(KMX, IMX), LPOMCMP(KMX, IMX), RPOMCMP(KMX, IMX))
    allocate(LPZOOINC(KMX, IMX), LPZOOOUTC(KMX, IMX))
    allocate(LDOP(KMX, IMX), RDOP(KMX, IMX), LPOP(KMX, IMX), RPOP(KMX, IMX), LDON(KMX, IMX), RDON(KMX, IMX), LPON(KMX, IMX), RPON(KMX, IMX))
    allocate(LDOC(KMX, IMX), RDOC(KMX, IMX), LPOC(KMX, IMX), RPOC(KMX, IMX))
    allocate(PSIEM(KMX, IMX), SEDEB(KMX, IMX))
    allocate(LPOMPEP(KMX, IMX), LPOMNEP(KMX, IMX), LPOMCEP(KMX, IMX))
    allocate(LPOMHD(KMX, IMX), RPOMHD(KMX, IMX))
    allocate(LDOMCAP(KMX, IMX), LDOMCEP(KMX, IMX), LPOMCAP(KMX, IMX), LPOMCNS(KMX, IMX), RPOMCNS(KMX, IMX))
    allocate(LDOMPD(KMX, IMX), LRDOMPD(KMX, IMX), RDOMPD(KMX, IMX), LPOMPD(KMX, IMX), LRPOMPD(KMX, IMX), RPOMPD(KMX, IMX), LPOMPHD(KMX, IMX), RPOMPHD(KMX, IMX))
    allocate(LDOMND(KMX, IMX), LRDOMND(KMX, IMX), RDOMND(KMX, IMX), LPOMND(KMX, IMX), LRPOMND(KMX, IMX), RPOMND(KMX, IMX), LPOMNHD(KMX, IMX), RPOMNHD(KMX, IMX))
    allocate(LDOMCD(KMX, IMX), LRDOMCD(KMX, IMX), RDOMCD(KMX, IMX), LPOMCD(KMX, IMX), LRPOMCD(KMX, IMX), RPOMCD(KMX, IMX), LPOMCHD(KMX, IMX), RPOMCHD(KMX, IMX))
!
!
    allocate(IceQSS(IMX)) ! CEMA
    allocate(WBSEG(IMX), PALT_JW(NWB), ELWS_INI(IMX), GTTYP(NGT), GTPC(NGT)) ! systdg 
! systdg 

!WAIT_FOR_TRIB_INPUT   -- logical array (trib index) used to determine which tributary has awaited input
!WAIT_FOR_BRANCH_INPUT -- logical array (branch index) used to determine which branch has awaited input
!TR_FILEDIR            -- character array (tributary index) to hold the directory names of any awaited tributary input files
!BR_FILEDIR            -- character array (branch index) to hold the directory names of any awaited branch input files
    allocate(WAIT_FOR_TRIB_INPUT(NTR), WAIT_FOR_BRANCH_INPUT(NBR), TR_FILEDIR(NTR), BR_FILEDIR(NBR)) !SR 11/26/19

! BIOENERGETICS !mlm
    if (FISHBIO) then
        allocate(BIOD(NOD), BIOF(NOD), biodp(NOD))
! BIOENERGETICS OUTPUT CARDS
        read(FISHBIOFN, "(//(8X,A8,2I8))") BIOC, NBIO, NIBIO
        allocate(BIOEXPFN(NIBIO), WEIGHTNUM(NIBIO), C2ZOO(KMX, IMX, NCT), VOLROOS(IMX), C2W(IMX, NCT + 2), IBIO(NIBIO))
        C2W = 0.0
        read(FISHBIOFN, "(//(:8X,9F8.0))") (BIOD(II), II = 1, NBIO)
        read(FISHBIOFN, "(//(:8X,9F8.0))") (BIOF(II), II = 1, NBIO)
        read(FISHBIOFN, "(//(:8X,9I8))") (IBIO(II), II = 1, NIBIO)
        read(FISHBIOFN, "(//(8x,a72))") biofn
        close(FISHBIOFN)
        weightfn = "weight.opt" ! for output filename generation
        bioexp = BIOC == "      ON"
        if (bioexp) then
            bhead(1) = "Jday99"
            bhead(2) = "Segmnt"
            bhead(3) = "Dpth_m"
            bhead(4) = "T_C"
            bhead(5) = "gamma"

            do j = 1, nzooe - nzoos + 1
                write(segnum, "(i0)") j
                SEGNUM = ADJUSTL(SEGNUM)
                L = LEN_TRIM(SEGNUM)
                bhead(5 + j) = "Zoo" // SEGNUM(1:L)
            end do
!bhead(6) = 'Zoo1'
!bhead(7) = 'Zoo2'	
            bhead(5 + j) = "K"
            bhead(6 + j) = "BH"
            bhead(7 + j) = "EL"
            bhead(8 + j) = "date"
        end if
! NVIOL CARD
! read(1222,'(//(8x,a8))') nviolc
! NVIOL_PRINT = NVIOLC        == '      ON'
! if(nviol_print) then
!   open(1333,file='nviol.dat',status='unknown')
!	allocate(nviol_loc(kmx,imx))
!	nviol_loc = 0
! end if
    end if
! Allocate subroutine variables

    call TRANSPORT()
    call WATERBODY()
    call OPEN_CHANNEL_INITIALIZE()
    call PIPE_FLOW_INITIALIZE()

! State variables

    TDS => C2(:, :, 1);     PO4 => C2(:, :, NPO4);     NH4 => C2(:, :, NNH4);     NO3 => C2(:, :, NNO3);     DSI => C2(:, :, NDSI)
    N2 => C2(:, :, NN2);     H2S => C2(:, :, NH2S);     CH4 => C2(:, :, NCH4);     SO4 => C2(:, :, NSO4)
    FEII => C2(:, :, NFEII);     FEOOH => C2(:, :, NFEOOH);     MNII => C2(:, :, NMNII);     MNO2 => C2(:, :, NMNO2)
    WAGE => C2(:, :, NWAGE);     BACT => C2(:, :, NBACT);     DGP => C2(:, :, NDGP)
    PSI => C2(:, :, NPSI)
    if (ORGC_CALC) then
        LDOMC => C2(:, :, NLDOMC);         RDOMC => C2(:, :, NRDOMC);         LPOMC => C2(:, :, NLPOMC);         RPOMC => C2(:, :, NRPOMC)
    else
        LDOM => C2(:, :, NLDOM);         RDOM => C2(:, :, NRDOM);         LPOM => C2(:, :, NLPOM);         RPOM => C2(:, :, NRPOM)
    end if
    O2 => C2(:, :, NDO);     TIC => C2(:, :, NTIC);     ALK => C2(:, :, NALK)
    CG => C2(:, :, NGCS:NGCE);     SS => C2(:, :, NSSS:NSSE);     ALG => C2(:, :, NAS:NAE)
    CBOD => C2(:, :, NBODCS:NBODCE);     CBODP => C2(:, :, NBODPS:NBODPE);     CBODN => C2(:, :, NBODNS:NBODNE) ! CB 6/6/10
    ZOO => C2(:, :, NZOOS:NZOOE)
    LDOMP => C2(:, :, NLDOMP);     RDOMP => C2(:, :, NRDOMP);     LPOMP => C2(:, :, NLPOMP);     RPOMP => C2(:, :, NRPOMP)
    LDOMN => C2(:, :, NLDOMN);     RDOMN => C2(:, :, NRDOMN);     LPOMN => C2(:, :, NLPOMN);     RPOMN => C2(:, :, NRPOMN)
    EX_TOXIN => C2(:, :, NATS:NATE) ! EXTRACELLULAR ALAGE TOXIN

! State variable source/sinks

    CGSS => CSSK(:, :, NGCS:NGCE);     SSSS => CSSK(:, :, NSSS:NSSE);     PO4SS => CSSK(:, :, NPO4);     NH4SS => CSSK(:, :, NNH4)
    N2SS => CSSK(:, :, NN2);     H2SSS => CSSK(:, :, NH2S);     CH4SS => CSSK(:, :, NCH4);     SO4SS => CSSK(:, :, NSO4)
    FEIISS => CSSK(:, :, NFEII);     FEOOHSS => CSSK(:, :, NFEOOH);     MNIISS => CSSK(:, :, NMNII);     MNO2SS => CSSK(:, :, NMNO2)
    AGESS => CSSK(:, :, NWAGE);     BACTSS => CSSK(:, :, NBACT);     DISGSS => CSSK(:, :, NDGP)
    NO3SS => CSSK(:, :, NNO3);     DSISS => CSSK(:, :, NDSI);     PSISS => CSSK(:, :, NPSI)
    if (ORGC_CALC) then
        LDOMCSS => CSSK(:, :, NLDOMC);         RDOMCSS => CSSK(:, :, NRDOMC);         LPOMCSS => CSSK(:, :, NLPOMC);         RPOMCSS => CSSK(:, :, NRPOMC)
    else
        LDOMSS => CSSK(:, :, NLDOM);         RDOMSS => CSSK(:, :, NRDOM);         LPOMSS => CSSK(:, :, NLPOM);         RPOMSS => CSSK(:, :, NRPOM)
    end if
    ASS => CSSK(:, :, NAS:NAE);     DOSS => CSSK(:, :, NDO);     TICSS => CSSK(:, :, NTIC)
    CBODSS => CSSK(:, :, NBODCS:NBODCE);     CBODPSS => CSSK(:, :, NBODPS:NBODPE);     CBODNSS => CSSK(:, :, NBODNS:NBODNE) ! CB 6/6/10
    ZOOSS => CSSK(:, :, NZOOS:NZOOE)
    LDOMPSS => CSSK(:, :, NLDOMP);     RDOMPSS => CSSK(:, :, NRDOMP);     LPOMPSS => CSSK(:, :, NLPOMP);     RPOMPSS => CSSK(:, :, NRPOMP)
    LDOMNSS => CSSK(:, :, NLDOMN);     RDOMNSS => CSSK(:, :, NRDOMN);     LPOMNSS => CSSK(:, :, NLPOMN);     RPOMNSS => CSSK(:, :, NRPOMN)
    alkss => CSSK(:, :, nalk) ! enhanced pH buffering
    CTESS => CSSK(:, :, NATS:NATE) ! ALGAE TOXINS

! Derived variables
    DOC_DER = 1;     POC_DER = 2;     TOC_DER = 3;     DON_DER = 4;     PON_DER = 5;     TON_DER = 6
    TKN_DER = 7;     TN_DER = 8;     NH3_DER = 9;     DOP_DER = 10;     POP_DER = 11;     TOP_DER = 12
    TP_DER = 13;     APR_DER = 14;     CHLA_DER = 15;     ATOT_DER = 16;     O2DG_DER = 17
    TDG_DER = 18;     TURB_DER = 19;     TOTSS_DER = 20
    TISS_DER = 21;     CBODU_DER = 22;     PH_DER = 23;     CO2_DER = 24;     HCO3_DER = 25;     CO3_DER = 26;     SECCHI_DER = 27

    DOC => CD(:, :, DOC_DER);     POC => CD(:, :, POC_DER);     TOC => CD(:, :, TOC_DER);     DON => CD(:, :, DON_DER)
    PON => CD(:, :, PON_DER);     TON => CD(:, :, TON_DER)
    TKN => CD(:, :, TKN_DER);     TN => CD(:, :, TN_DER);     NH3 => CD(:, :, NH3_DER)
    DOP => CD(:, :, DOP_DER);     POP => CD(:, :, POP_DER);     TOP => CD(:, :, TOP_DER);     TP => CD(:, :, TP_DER)
    APR => CD(:, :, APR_DER);     CHLA => CD(:, :, CHLA_DER);     ATOT => CD(:, :, ATOT_DER);     O2DG => CD(:, :, O2DG_DER); 
    TDG => CD(:, :, TDG_DER);     TURB => CD(:, :, TURB_DER)
    TOTSS => CD(:, :, TOTSS_DER);     TISS => CD(:, :, TISS_DER)
    CBODU => CD(:, :, CBODU_DER);     PH => CD(:, :, PH_DER);     CO2 => CD(:, :, CO2_DER)
    HCO3 => CD(:, :, HCO3_DER);     CO3 => CD(:, :, CO3_DER); ;     SECCHID => CD(:, :, SECCHI_DER)

! Kinetic fluxes

    KF_PO4_SD = 12;     KF_PO4_SR = 13
    KF_NH4_SD = 25;     KF_NH4_SR = 26;     KF_NH3GAS = 27;     KF_NO3D = 28;     KF_NO3AG = 29;     KF_NO3EG = 30;     KF_NO3SED = 31
    KF_DO_SED = 64;     KF_DO_SOD = 65
    KF_SED_PBURIAL = 115;     KF_SED_NBURIAL = 116;     KF_CO2X = 120;     KF_DOH2S = 121;     KF_DOCH4 = 125;     KF_FE2D = 129;     KF_MN2D = 133;     KF_SDINC = 137
    KF_SEDD = 141

    SSSI => KF(:, :, 1);     SSSO => KF(:, :, 2);     PO4AR => KF(:, :, 3);     PO4AG => KF(:, :, 4);     PO4AP => KF(:, :, 5)
    PO4ER => KF(:, :, 6);     PO4EG => KF(:, :, 7);     PO4EP => KF(:, :, 8);     PO4POM => KF(:, :, 9);     PO4DOM => KF(:, :, 10)
    PO4OM => KF(:, :, 11);     PO4SD => KF(:, :, KF_PO4_SD);     PO4SR => KF(:, :, KF_PO4_SR);     PO4NS => KF(:, :, 14);     NH4D => KF(:, :, 15)
    NH4AR => KF(:, :, 16);     NH4AG => KF(:, :, 17);     NH4AP => KF(:, :, 18);     NH4ER => KF(:, :, 19);     NH4EG => KF(:, :, 20)
    NH4EP => KF(:, :, 21);     NH4POM => KF(:, :, 22);     NH4DOM => KF(:, :, 23);     NH4OM => KF(:, :, 24);     NH4SD => KF(:, :, KF_NH4_SD)
    NH4SR => KF(:, :, KF_NH4_SR);     NH3GAS => KF(:, :, KF_NH3GAS)

    NO3D => KF(:, :, KF_NO3D);     NO3AG => KF(:, :, KF_NO3AG);     NO3EG => KF(:, :, KF_NO3EG);     NO3SED => KF(:, :, KF_NO3SED)
    DSIAG => KF(:, :, 32);     DSIEG => KF(:, :, 33);     DSID => KF(:, :, 34);     DSISD => KF(:, :, 35);     DSISR => KF(:, :, 36)
    DSIS => KF(:, :, 37);     PSIAM => KF(:, :, 38);     PSINS => KF(:, :, 39);     PSID => KF(:, :, 40)
    LDOMD => KF(:, :, 41);     LRDOMD => KF(:, :, 42);     RDOMD => KF(:, :, 43);     LDOMAP => KF(:, :, 44)
    LDOMEP => KF(:, :, 45);     LPOMD => KF(:, :, 46);     LRPOMD => KF(:, :, 47);     RPOMD => KF(:, :, 48);     LPOMAP => KF(:, :, 49)
    LPOMEP => KF(:, :, 50);     LPOMNS => KF(:, :, 51);     RPOMNS => KF(:, :, 52);     CBODDK => KF(:, :, 53);     DOAP => KF(:, :, 54)
    DOEP => KF(:, :, 55);     DOAR => KF(:, :, 56);     DOER => KF(:, :, 57);     DOPOM => KF(:, :, 58);     DODOM => KF(:, :, 59)
    DOOM => KF(:, :, 60);     DONIT => KF(:, :, 61);     DOBOD => KF(:, :, 62);     DOAE => KF(:, :, 63);     DOSED => KF(:, :, KF_DO_SED)
    DOSOD => KF(:, :, KF_DO_SOD);     TICAP => KF(:, :, 66);     TICEP => KF(:, :, 67);     SEDD => KF(:, :, 68);     SEDAS => KF(:, :, 69)
    SEDOMS => KF(:, :, 70);     SEDNS => KF(:, :, 71);     SODD => KF(:, :, 72)

    LDOMPAP => KF(:, :, 73);     LDOMPeP => KF(:, :, 74);     LPOMpAP => KF(:, :, 75);     LPOMPNS => KF(:, :, 76);     RPOMPNS => KF(:, :, 77)
    LDOMnAP => KF(:, :, 78);     LDOMneP => KF(:, :, 79);     LPOMnAP => KF(:, :, 80);     LPOMnNS => KF(:, :, 81);     RPOMnNS => KF(:, :, 82)
    SEDDp => KF(:, :, 83);     SEDASp => KF(:, :, 84);     SEDOMSp => KF(:, :, 85);     SEDNSp => KF(:, :, 86);     lpomepp => KF(:, :, 87)
    SEDDn => KF(:, :, 88);     SEDASn => KF(:, :, 89);     SEDOMSn => KF(:, :, 90);     SEDNSn => KF(:, :, 91);     lpomepn => KF(:, :, 92)
    SEDDc => KF(:, :, 93);     SEDASc => KF(:, :, 94);     SEDOMSc => KF(:, :, 95);     SEDNSc => KF(:, :, 96);     lpomepc => KF(:, :, 97)
    SEDNO3 => KF(:, :, 98)

    PO4MR => KF(:, :, 99);     PO4MG => KF(:, :, 100);     NH4MR => KF(:, :, 101);     NH4MG => KF(:, :, 102);     LDOMMAC => KF(:, :, 103)
    RPOMMAC => KF(:, :, 104);     LPOMMAC => KF(:, :, 105);     DOMP => KF(:, :, 106);     DOMR => KF(:, :, 107);     TICMC => KF(:, :, 108)
    CBODNS => KF(:, :, 109);     SEDCB => KF(:, :, 110);     SEDCBP => KF(:, :, 111);     SEDCBN => KF(:, :, 112);     SEDCBC => KF(:, :, 113)
    SEDBR => KF(:, :, 114);     SEDBRP => KF(:, :, KF_SED_PBURIAL);     SEDBRN => KF(:, :, KF_SED_NBURIAL);     SEDBRC => KF(:, :, 117); 
    CBODNSP => KF(:, :, 118);     CBODNSN => KF(:, :, 119); ;     CO2REAER => KF(:, :, KF_CO2X)

    DOH2S => KF(:, :, KF_DOH2S);     H2SREAER => KF(:, :, 122);     H2SD => KF(:, :, 123);     H2SSR => KF(:, :, 124)
    DOCH4 => KF(:, :, KF_DOCH4);     CH4REAER => KF(:, :, 126);     CH4D => KF(:, :, 127);     CH4SR => KF(:, :, 128)
    FE2D => KF(:, :, KF_FE2D);     DOFE2 => KF(:, :, 130);     FEIISR => KF(:, :, 131);     SDINFEOOH => KF(:, :, 132)
    MN2D => KF(:, :, KF_MN2D);     DOMN2 => KF(:, :, 134);     MNIISR => KF(:, :, 135);     SDINMNO2 => KF(:, :, 136)

    SDINC => KF(:, :, KF_SDINC);     SDINN => KF(:, :, 138);     SDINP => KF(:, :, 139);     DOSEDIA => KF(:, :, 140); 
    SEDD1 => KF(:, :, KF_SEDD);     SEDD2 => KF(:, :, 142)


! Algal rate variables

    AGR => AF(:, :, :, 1);     ARR => AF(:, :, :, 2);     AER => AF(:, :, :, 3);     AMR => AF(:, :, :, 4);     ASR => AF(:, :, :, 5)
    EGR => EF(:, :, :, 1);     ERR => EF(:, :, :, 2);     EER => EF(:, :, :, 3);     EMR => EF(:, :, :, 4);     EBR => EF(:, :, :, 5)

! Hydrodynamic variables

    DLTLIM => HYD(:, :, 1);     U => HYD(:, :, 2);     W => HYD(:, :, 3);     T2 => HYD(:, :, 4);     RHO => HYD(:, :, 5);     AZ => HYD(:, :, 6)
    VSH => HYD(:, :, 7);     ST => HYD(:, :, 8);     SB => HYD(:, :, 9);     ADMX => HYD(:, :, 10);     DM => HYD(:, :, 11);     HDG => HYD(:, :, 12)
    ADMZ => HYD(:, :, 13);     HPG => HYD(:, :, 14);     GRAV => HYD(:, :, 15)

! I/O units

    SNP => OPT(:, 1);     PRF => OPT(:, 2);     VPL => OPT(:, 3);     CPL => OPT(:, 4);     SPR => OPT(:, 5);     FLX => OPT(:, 6);     FLX2 => OPT(:, 7);     SPRV => OPT(:, 8) ! SW 9/27/2018

! Zero variables

    ITR = 0;     JBTR = 0;     KTTR = 0;     KBTR = 0;     QTR = 0.0;     TTR = 0.0;     CTR = 0.0;     QTRF = 0.0;     SNPD = 0.0;     TSRD = 0.0
    PRFD = 0.0;     SPRD = 0.0;     CPLD = 0.0;     VPLD = 0.0;     SCRD = 0.0;     FLXD = 0.0;     WDOD = 0.0;     RSOD = 0.0;     ELTRB = 0.0;     ELTRT = 0.0

    KFNAME2 = "     " ! SW 9/27/13 INITIALIZE ENTIRE ARRAY
    KFWBC = "     " ! SW 9/27/13 INITIALIZE ENTIRE ARRAY

! Input file unit numbers

    NUNIT = 40
    do JW = 1, NWB
        BTH(JW) = NUNIT
        VPR(JW) = NUNIT + 1
        LPR(JW) = NUNIT + 2
        NUNIT = NUNIT + 3
    end do
    GRF = NUNIT;     NUNIT = NUNIT + 1

! Time control cards
    if (CONFN == "w2_con.npt") then
        read(CON, "(//8X,2F8.0,I8)") TMSTRT, TMEND, YEAR
        read(CON, "(//8X,I8,F8.0,a8)") NDLT, DLTMIN, DLTINTER;         DLTD = 0.0 ! SW 9/28/13 INITIALIZE ARRAY TO NOD SINCE ONLY NDLT ASSIGNED
        read(CON, "(//(:8X,9F8.0))") (DLTD(J), J = 1, NDLT)
        read(CON, "(//(:8X,9F8.0))") (DLTMAX(J), J = 1, NDLT)
        read(CON, "(//(:8X,9F8.0))") (DLTF(J), J = 1, NDLT)
        read(CON, "(//(8X,3A8))") (VISC(JW), CELC(JW), DLTADD(JW), JW = 1, NWB)

! Grid definition cards

        read(CON, "(//(8X,5I8,F8.0,F8.0))") (US(JB), DS(JB), UHS(JB), DHS(JB), NL(JB), SLOPE(JB), SLOPEC(JB), JB = 1, NBR)
        read(CON, "(//(8X,3F8.0,3I8))") (LAT(JW), LONGIT(JW), ELBOT(JW), BS(JW), BE(JW), JBDN(JW), JW = 1, NWB)

! Initial condition cards

        read(CON, "(//(8X,2F8.0,2A8))") (T2I(JW), ICETHI(JW), WTYPEC(JW), GRIDC(JW), JW = 1, NWB)
        read(CON, "(//(8X,6A8))") (VBC(JW), EBC(JW), MBC(JW), PQC(JW), EVC(JW), PRC(JW), JW = 1, NWB)
        read(CON, "(//(8X,4A8))") (WINDC(JW), QINC(JW), QOUTC(JW), HEATC(JW), JW = 1, NWB)
        read(CON, "(//(8X,3A8))") (QINIC(JB), DTRIC(JB), HDIC(JB), JB = 1, NBR)
        read(CON, "(//(8X,5A8,4F8.0))") (SLHTC(JW), SROC(JW), RHEVC(JW), METIC(JW), FETCHC(JW), AFW(JW), BFW(JW), CFW(JW), WINDH(JW), JW = 1, NWB)
        read(CON, "(//(8X,2A8,6F8.0))") (ICEC(JW), SLICEC(JW), ALBEDO(JW), HWI(JW), BETAI(JW), GAMMAI(JW), ICEMIN(JW), ICET2(JW), JW = 1, NWB)
        read(CON, "(//(8X,A8,F8.0))") (SLTRC(JW), THETA(JW), JW = 1, NWB)
        read(CON, "(//(8X,6F8.0,A8,F8.0))") (AX(JW), DXI(JW), CBHE(JW), TSED(JW), FI(JW), TSEDF(JW), FRICC(JW), Z0(JW), JW = 1, NWB)
        read(CON, "(//(8X,2A8,F8.0,I8,F8.0,F8.0,F8.0,F8.0,A8))") (AZC(JW), AZSLC(JW), AZMAX(JW), TKEBC(JW), EROUGH(JW), ARODI(JW), STRICK(JW), TKELATPRDCONST(JW), IMPTKE(JW), JW = 1, NWB) !,PHISET(JW
    else ! CSV INPUT FILE 
        read(CON, *)
        read(CON, *)
        read(CON, *) TMSTRT, TMEND, YEAR
        read(CON, *)
        read(CON, *)
        read(CON, *) NDLT, DLTMIN, DLTINTER;         DLTD = 0.0;         DLTINTER = ADJUSTR(DLTINTER)
        read(CON, *)
        read(CON, *)
        read(CON, *) (DLTD(J), J = 1, NDLT)
        read(CON, *)
        read(CON, *)
        read(CON, *) (DLTMAX(J), J = 1, NDLT)
        read(CON, *)
        read(CON, *)
        read(CON, *) (DLTF(J), J = 1, NDLT)
        read(CON, *)
        read(CON, *)
        read(CON, *) (VISC(JW), JW = 1, NWB)
        read(CON, *) (CELC(JW), JW = 1, NWB)
        read(CON, *) (DLTADD(JW), JW = 1, NWB)
        VISC = ADJUSTR(VISC);         CELC = ADJUSTR(CELC);         DLTADD = ADJUSTR(DLTADD)

! Grid definition cards
        read(CON, *)
        read(CON, *)
        read(CON, *) (US(JB), JB = 1, NBR)
        read(CON, *) (DS(JB), JB = 1, NBR)
        read(CON, *) (UHS(JB), JB = 1, NBR)
        read(CON, *) (DHS(JB), JB = 1, NBR)
        read(CON, *) (NL(JB), JB = 1, NBR)
        read(CON, *) (SLOPE(JB), JB = 1, NBR)
        read(CON, *) (SLOPEC(JB), JB = 1, NBR)
        read(CON, *)
        read(CON, *)


        read(CON, *) (LAT(JW), JW = 1, NWB)
        read(CON, *) (LONGIT(JW), JW = 1, NWB)
        read(CON, *) (ELBOT(JW), JW = 1, NWB)
        read(CON, *) (BS(JW), JW = 1, NWB)
        read(CON, *) (BE(JW), JW = 1, NWB)
        read(CON, *) (JBDN(JW), JW = 1, NWB)
        read(CON, *)
        read(CON, *)

! Initial condition cards

        read(CON, *) (T2I(JW), JW = 1, NWB)
        read(CON, *) (ICETHI(JW), JW = 1, NWB)
        read(CON, *) (WTYPEC(JW), JW = 1, NWB);         WTYPEC = ADJUSTR(WTYPEC)
        read(CON, *) (GRIDC(JW), JW = 1, NWB);         GRIDC = ADJUSTR(GRIDC)
        read(CON, *)
        read(CON, *)

        read(CON, *) (VBC(JW), JW = 1, NWB);         VBC = ADJUSTR(VBC)
        read(CON, *) (EBC(JW), JW = 1, NWB);         EBC = ADJUSTR(EBC)
        read(CON, *) (MBC(JW), JW = 1, NWB);         MBC = ADJUSTR(MBC)
        read(CON, *) (PQC(JW), JW = 1, NWB);         PQC = ADJUSTR(PQC)
        read(CON, *) (EVC(JW), JW = 1, NWB);         EVC = ADJUSTR(EVC)
        read(CON, *) (PRC(JW), JW = 1, NWB);         PRC = ADJUSTR(PRC)
        read(CON, *)
        read(CON, *)

        read(CON, *) (WINDC(JW), JW = 1, NWB);         WINDC = ADJUSTR(WINDC)
        read(CON, *) (QINC(JW), JW = 1, NWB);         QINC = ADJUSTR(QINC)
        read(CON, *) (QOUTC(JW), JW = 1, NWB);         QOUTC = ADJUSTR(QOUTC)
        read(CON, *) (HEATC(JW), JW = 1, NWB);         HEATC = ADJUSTR(HEATC)
        read(CON, *)
        read(CON, *)

        read(CON, *) (QINIC(JB), JB = 1, NBR);         QINIC = ADJUSTR(QINIC)
        read(CON, *) (DTRIC(JB), JB = 1, NBR);         DTRIC = ADJUSTR(DTRIC)
        read(CON, *) (HDIC(JB), JB = 1, NBR);         HDIC = ADJUSTR(HDIC)
        read(CON, *)
        read(CON, *)

        read(CON, *) (SLHTC(JW), JW = 1, NWB);         SLHTC = ADJUSTR(SLHTC)
        read(CON, *) (SROC(JW), JW = 1, NWB);         SROC = ADJUSTR(SROC)
        read(CON, *) (RHEVC(JW), JW = 1, NWB);         RHEVC = ADJUSTR(RHEVC)
        read(CON, *) (METIC(JW), JW = 1, NWB);         METIC = ADJUSTR(METIC)
        read(CON, *) (FETCHC(JW), JW = 1, NWB);         FETCHC = ADJUSTR(FETCHC)
        read(CON, *) (AFW(JW), JW = 1, NWB)
        read(CON, *) (BFW(JW), JW = 1, NWB)
        read(CON, *) (CFW(JW), JW = 1, NWB)
        read(CON, *) (WINDH(JW), JW = 1, NWB)
        read(CON, *)
        read(CON, *)

        read(CON, *) (ICEC(JW), JW = 1, NWB);         ICEC = ADJUSTR(ICEC)
        read(CON, *) (SLICEC(JW), JW = 1, NWB);         SLICEC = ADJUSTR(SLICEC)
        read(CON, *) (ALBEDO(JW), JW = 1, NWB)
        read(CON, *) (HWI(JW), JW = 1, NWB)
        read(CON, *) (BETAI(JW), JW = 1, NWB)
        read(CON, *) (GAMMAI(JW), JW = 1, NWB)
        read(CON, *) (ICEMIN(JW), JW = 1, NWB)
        read(CON, *) (ICET2(JW), JW = 1, NWB)
        read(CON, *)
        read(CON, *)

        read(CON, *) (SLTRC(JW), JW = 1, NWB);         SLTRC = ADJUSTR(SLTRC)
        read(CON, *) (THETA(JW), JW = 1, NWB)
        read(CON, *)
        read(CON, *)

        read(CON, *) (AX(JW), JW = 1, NWB)
        read(CON, *) (DXI(JW), JW = 1, NWB)
        read(CON, *) (CBHE(JW), JW = 1, NWB)
        read(CON, *) (TSED(JW), JW = 1, NWB)
        read(CON, *) (FI(JW), JW = 1, NWB)
        read(CON, *) (TSEDF(JW), JW = 1, NWB)
        read(CON, *) (FRICC(JW), JW = 1, NWB);         FRICC = ADJUSTR(FRICC)
        read(CON, *) (Z0(JW), JW = 1, NWB)
        read(CON, *)
        read(CON, *)

        read(CON, *) (AZC(JW), JW = 1, NWB);         AZC = ADJUSTR(AZC)
        read(CON, *) (AZSLC(JW), JW = 1, NWB);         AZSLC = ADJUSTR(AZSLC)
        read(CON, *) (AZMAX(JW), JW = 1, NWB)
        read(CON, *) (TKEBC(JW), JW = 1, NWB)
        read(CON, *) (EROUGH(JW), JW = 1, NWB)
        read(CON, *) (ARODI(JW), JW = 1, NWB)
        read(CON, *) (STRICK(JW), JW = 1, NWB)
        read(CON, *) (TKELATPRDCONST(JW), JW = 1, NWB)
        read(CON, *) (IMPTKE(JW), JW = 1, NWB);         IMPTKE = ADJUSTR(IMPTKE)

    end if

    do JW = 1, NWB
        if (Z0(JW) <= 0.0) then
            Z0(JW) = 0.001
        end if ! SW 11/28/07
        do JB = BS(JW), BE(JW)
            do I = US(JB), DS(JB)
                E(I) = EROUGH(JW)
            end do
        end do
    end do

    if (CONFN == "w2_con.npt") then
! Inflow-outflow cards

        read(CON, "(//(8X,I8,A8))") (NSTR(JB), DYNSTRUC(JB), JB = 1, NBR)
        read(CON, "(/)")
        do JB = 1, NBR
            read(CON, "(:8X,9A8)") (STRIC(JS, JB), JS = 1, NSTR(JB))
        end do
        read(CON, "(/)")
        do JB = 1, NBR
            read(CON, "(:8X,9I8)") (KTSWT(JS, JB), JS = 1, NSTR(JB))
        end do
        read(CON, "(/)")
        do JB = 1, NBR
            read(CON, "(:8X,9I8)") (KBSWT(JS, JB), JS = 1, NSTR(JB))
        end do
        read(CON, "(/)")
        do JB = 1, NBR
            read(CON, "(:8X,9A8)") (SINKCT(JS, JB), JS = 1, NSTR(JB))
        end do
        read(CON, "(/)")
        do JB = 1, NBR
            read(CON, "(:8X,9F8.0)") (ESTRT(JS, JB), JS = 1, NSTR(JB))
        end do
        read(CON, "(/)")
        do JB = 1, NBR
            read(CON, "(:8X,9F8.0)") (WSTRT(JS, JB), JS = 1, NSTR(JB))
        end do
        read(CON, "(//(:8X,2I8,6F8.0,A8,A8))") (IUPI(JP), IDPI(JP), EUPI(JP), EDPI(JP), WPI(JP), DLXPI(JP), FPI(JP), FMINPI(JP), LATPIC(JP), DYNPIPE(JP), JP = 1, NPI)
        read(CON, "(//(:8X,A8,2F8.0,2I8))") (PUPIC(JP), ETUPI(JP), EBUPI(JP), KTUPI(JP), KBUPI(JP), JP = 1, NPI)
        read(CON, "(//(:8X,A8,2F8.0,2I8))") (PDPIC(JP), ETDPI(JP), EBDPI(JP), KTDPI(JP), KBDPI(JP), JP = 1, NPI)
        read(CON, "(//(:8X,2I8,5F8.0,A8))") (IUSP(JS), IDSP(JS), ESP(JS), A1SP(JS), B1SP(JS), A2SP(JS), B2SP(JS), LATSPC(JS), JS = 1, NSP)
        read(CON, "(//(:8X,A8,2F8.0,2I8))") (PUSPC(JS), ETUSP(JS), EBUSP(JS), KTUSP(JS), KBUSP(JS), JS = 1, NSP)
        read(CON, "(//(:8X,A8,2F8.0,2I8))") (PDSPC(JS), ETDSP(JS), EBDSP(JS), KTDSP(JS), KBDSP(JS), JS = 1, NSP)
        read(CON, "(//(:8X,A8,I8,3F8.0))") (GASSPC(JS), EQSP(JS), AGASSP(JS), BGASSP(JS), CGASSP(JS), JS = 1, NSP)
        read(CON, "(//(:8X,2I8,7F8.0,A8))") (IUGT(JG), IDGT(JG), EGT(JG), A1GT(JG), B1GT(JG), G1GT(JG), A2GT(JG), B2GT(JG), G2GT(JG), LATGTC(JG), JG = 1, NGT)
        read(CON, "(//(:8X,4F8.0,2A8))") (GTA1(JG), GTB1(JG), GTA2(JG), GTB2(JG), DYNGTC(JG), GTIC(JG), JG = 1, NGT) ! cb 8/13/2010
        read(CON, "(//(:8X,A8,2F8.0,2I8))") (PUGTC(JG), ETUGT(JG), EBUGT(JG), KTUGT(JG), KBUGT(JG), JG = 1, NGT)
        read(CON, "(//(:8X,A8,2F8.0,2I8))") (PDGTC(JG), ETDGT(JG), EBDGT(JG), KTDGT(JG), KBDGT(JG), JG = 1, NGT)
        read(CON, "(//(:8X,A8,I8,3F8.0))") (GASGTC(JG), EQGT(JG), AGASGT(JG), BGASGT(JG), CGASGT(JG), JG = 1, NGT)

        read(CON, "(//(:8X,2I8,6F8.0,2A8))") (IUPU(JP), IDPU(JP), EPU(JP), STRTPU(JP), ENDPU(JP), EONPU(JP), EOFFPU(JP), QPU(JP), LATPUC(JP), DYNPUMP(JP), JP = 1, NPU)
    else ! w2_con.csv file format
! Inflow-outflow cards

        read(CON, *)
        read(CON, *)
        read(CON, *) (NSTR(JB), JB = 1, NBR)
        read(CON, *) (DYNSTRUC(JB), JB = 1, NBR);         DYNSTRUC = adjustr(DYNSTRUC)

        do JS = 1, NSTT
            read(CON, *) (STRIC(JS, JB), JB = 1, NBR)
        end do
        STRIC = adjustr(STRIC)
        do JS = 1, NSTT
            read(CON, *) (KTSWT(JS, JB), JB = 1, NBR)
        end do
        do JS = 1, NSTT
            read(CON, *) (KBSWT(JS, JB), JB = 1, NBR)
        end do
        do JS = 1, NSTT
            read(CON, *) (SINKCT(JS, JB), JB = 1, NBR)
        end do
        SINKCT = adjustr(SINKCT)
        do JS = 1, NSTT
            read(CON, *) (ESTRT(JS, JB), JB = 1, NBR)
        end do
        do JS = 1, NSTT
            read(CON, *) (WSTRT(JS, JB), JB = 1, NBR)
        end do
        read(CON, *)
        read(CON, *)
        read(CON, *) (IUPI(JP), JP = 1, NPI)
        read(CON, *) (IDPI(JP), JP = 1, NPI)
        read(CON, *) (EUPI(JP), JP = 1, NPI)
        read(CON, *) (EDPI(JP), JP = 1, NPI)
        read(CON, *) (WPI(JP), JP = 1, NPI)
        read(CON, *) (DLXPI(JP), JP = 1, NPI)
        read(CON, *) (FPI(JP), JP = 1, NPI)
        read(CON, *) (FMINPI(JP), JP = 1, NPI)
        read(CON, *) (LATPIC(JP), JP = 1, NPI);         LATPIC = adjustr(LATPIC)
        read(CON, *) (DYNPIPE(JP), JP = 1, NPI);         DYNPIPE = adjustr(DYNPIPE)

        read(CON, *) (PUPIC(JP), JP = 1, NPI);         PUPIC = adjustr(PUPIC)
        read(CON, *) (ETUPI(JP), JP = 1, NPI)
        read(CON, *) (EBUPI(JP), JP = 1, NPI)
        read(CON, *) (KTUPI(JP), JP = 1, NPI)
        read(CON, *) (KBUPI(JP), JP = 1, NPI)

        read(CON, *) (PDPIC(JP), JP = 1, NPI);         PDPIC = adjustr(PDPIC)
        read(CON, *) (ETDPI(JP), JP = 1, NPI)
        read(CON, *) (EBDPI(JP), JP = 1, NPI)
        read(CON, *) (KTDPI(JP), JP = 1, NPI)
        read(CON, *) (KBDPI(JP), JP = 1, NPI)

        read(CON, *)
        read(CON, *)
        read(CON, *) (IUSP(JS), JS = 1, NSP)
        read(CON, *) (IDSP(JS), JS = 1, NSP)
        read(CON, *) (ESP(JS), JS = 1, NSP)
        read(CON, *) (A1SP(JS), JS = 1, NSP)
        read(CON, *) (B1SP(JS), JS = 1, NSP)
        read(CON, *) (A2SP(JS), JS = 1, NSP)
        read(CON, *) (B2SP(JS), JS = 1, NSP)
        read(CON, *) (LATSPC(JS), JS = 1, NSP);         LATSPC = ADJUSTR(LATSPC)

        read(CON, *) (PUSPC(JS), JS = 1, NSP);         PUSPC = adjustr(PUSPC)
        read(CON, *) (ETUSP(JS), JS = 1, NSP)
        read(CON, *) (EBUSP(JS), JS = 1, NSP)
        read(CON, *) (KTUSP(JS), JS = 1, NSP)
        read(CON, *) (KBUSP(JS), JS = 1, NSP)

        read(CON, *) (PDSPC(JS), JS = 1, NSP);         PDSPC = ADJUSTR(PDSPC)
        read(CON, *) (ETDSP(JS), JS = 1, NSP)
        read(CON, *) (EBDSP(JS), JS = 1, NSP)
        read(CON, *) (KTDSP(JS), JS = 1, NSP)
        read(CON, *) (KBDSP(JS), JS = 1, NSP)

        read(CON, *) (GASSPC(JS), JS = 1, NSP);         GASSPC = ADJUSTR(GASSPC)
        read(CON, *) (EQSP(JS), JS = 1, NSP)
        read(CON, *) (AGASSP(JS), JS = 1, NSP)
        read(CON, *) (BGASSP(JS), JS = 1, NSP)
        read(CON, *) (CGASSP(JS), JS = 1, NSP)

        read(CON, *)
        read(CON, *)
        read(CON, *) (IUGT(JG), JG = 1, NGT)
        read(CON, *) (IDGT(JG), JG = 1, NGT)
        read(CON, *) (EGT(JG), JG = 1, NGT)
        read(CON, *) (A1GT(JG), JG = 1, NGT)
        read(CON, *) (B1GT(JG), JG = 1, NGT)
        read(CON, *) (G1GT(JG), JG = 1, NGT)
        read(CON, *) (A2GT(JG), JG = 1, NGT)
        read(CON, *) (B2GT(JG), JG = 1, NGT)
        read(CON, *) (G2GT(JG), JG = 1, NGT)
        read(CON, *) (LATGTC(JG), JG = 1, NGT);         LATGTC = ADJUSTR(LATGTC)

        read(CON, *) (GTA1(JG), JG = 1, NGT)
        read(CON, *) (GTB1(JG), JG = 1, NGT)
        read(CON, *) (GTA2(JG), JG = 1, NGT)
        read(CON, *) (GTB2(JG), JG = 1, NGT)
        read(CON, *) (DYNGTC(JG), JG = 1, NGT);         DYNGTC = ADJUSTR(DYNGTC)
        read(CON, *) (GTIC(JG), JG = 1, NGT);         GTIC = ADJUSTR(GTIC)

        read(CON, *) (PUGTC(JG), JG = 1, NGT);         PUGTC = ADJUSTR(PUGTC)
        read(CON, *) (ETUGT(JG), JG = 1, NGT)
        read(CON, *) (EBUGT(JG), JG = 1, NGT)
        read(CON, *) (KTUGT(JG), JG = 1, NGT)
        read(CON, *) (KBUGT(JG), JG = 1, NGT)
        read(CON, *) (PDGTC(JG), JG = 1, NGT);         PDGTC = ADJUSTR(PDGTC)
        read(CON, *) (ETDGT(JG), JG = 1, NGT)
        read(CON, *) (EBDGT(JG), JG = 1, NGT)
        read(CON, *) (KTDGT(JG), JG = 1, NGT)
        read(CON, *) (KBDGT(JG), JG = 1, NGT)

        read(CON, *) (GASGTC(JG), JG = 1, NGT);         GASGTC = ADJUSTR(GASGTC)
        read(CON, *) (EQGT(JG), JG = 1, NGT)
        read(CON, *) (AGASGT(JG), JG = 1, NGT)
        read(CON, *) (BGASGT(JG), JG = 1, NGT)
        read(CON, *) (CGASGT(JG), JG = 1, NGT)

        read(CON, *)
        read(CON, *)

        read(CON, *) (IUPU(JP), JP = 1, NPU)
        read(CON, *) (IDPU(JP), JP = 1, NPU)
        read(CON, *) (EPU(JP), JP = 1, NPU)
        read(CON, *) (STRTPU(JP), JP = 1, NPU)
        read(CON, *) (ENDPU(JP), JP = 1, NPU)
        read(CON, *) (EONPU(JP), JP = 1, NPU)
        read(CON, *) (EOFFPU(JP), JP = 1, NPU)
        read(CON, *) (QPU(JP), JP = 1, NPU)
        read(CON, *) (LATPUC(JP), JP = 1, NPU);         LATPUC = adjustr(LATPUC)
        read(CON, *) (DYNPUMP(JP), JP = 1, NPU);         DYNPUMP = ADJUSTR(DYNPUMP)

    end if

! Pump level based on downstream location SW 2/19/2020 rather than upstream one
    PUMP_DOWNSTREAM = .false.
    do JP = 1, NPU
        if (IDPU(JP) < 0) then
            PUMP_DOWNSTREAM(JP) = .true.
            IDPU(JP) = ABS(IDPU(JP))
        end if
    end do

    if (CONFN == "w2_con.npt") then

        read(CON, "(//(:8X,A8,2F8.0,2I8))") (PPUC(JP), ETPU(JP), EBPU(JP), KTPU(JP), KBPU(JP), JP = 1, NPU)
        read(CON, "(//(:8X,9I8))") (IWR(JW), JW = 1, NIW)
        read(CON, "(//(:8X,9F8.0))") (EKTWR(JW), JW = 1, NIW) ! SW 3/18/16
        read(CON, "(//(:8X,9F8.0))") (EKBWR(JW), JW = 1, NIW) ! SW 3/18/16
        read(CON, "(//(:8X,9A8))") (WDIC(JW), JW = 1, NWD)
        read(CON, "(//(:8X,9I8))") (IWD(JW), JW = 1, NWD)
        read(CON, "(//(:8X,9F8.0))") (EWD(JW), JW = 1, NWD)
        read(CON, "(//(:8X,9I8))") (KTWD(JW), JW = 1, NWD)
        read(CON, "(//(:8X,9I8))") (KBWD(JW), JW = 1, NWD);         TRC = "      " ! SW 9/27/13 INITIALIZATION SINCE ALLOCATION IS TO NTRT
        read(CON, "(//(:8X,9A8))") (TRC(JT), JT = 1, NTR)
        read(CON, "(//(:8X,9A8))") (TRIC(JT), JT = 1, NTR)
        read(CON, "(//(:8X,9I8))") (ITR(JT), JT = 1, NTR)
        read(CON, "(//(:8X,9F8.0))") (ELTRT(JT), JT = 1, NTR)
        read(CON, "(//(:8X,9F8.0))") (ELTRB(JT), JT = 1, NTR)
        read(CON, "(//(8X,A8))") (DTRC(JB), JB = 1, NBR)

! Output control cards (excluding constituents)

        read(CON, "(/)")
        do JH = 1, NHY
            read(CON, "(:8X,9A8)") (HPRWBC(JH, JW), JW = 1, NWB)
        end do
        read(CON, "(//(8X,A8,2I8))") (SNPC(JW), NSNP(JW), NISNP(JW), JW = 1, NWB)
        read(CON, "(/)")
        do JW = 1, NWB
            read(CON, "(:8X,9F8.0)") (SNPD(J, JW), J = 1, NSNP(JW))
        end do
        read(CON, "(/)")
        do JW = 1, NWB
            read(CON, "(:8X,9F8.0)") (SNPF(J, JW), J = 1, NSNP(JW))
        end do
        read(CON, "(/)")
        do JW = 1, NWB
            read(CON, "(:8X,9I8)") (ISNP(I, JW), I = 1, NISNP(JW))
        end do
        read(CON, "(//(8X,A8,I8))") (SCRC(JW), NSCR(JW), JW = 1, NWB)
        read(CON, "(/)")
        do JW = 1, NWB
            read(CON, "(:8X,9F8.0)") (SCRD(J, JW), J = 1, NSCR(JW))
        end do
        read(CON, "(/)")
        do JW = 1, NWB
            read(CON, "(:8X,9F8.0)") (SCRF(J, JW), J = 1, NSCR(JW))
        end do
        read(CON, "(//(8X,A8,2I8))") (PRFC(JW), NPRF(JW), NIPRF(JW), JW = 1, NWB)
        read(CON, "(/)")
        do JW = 1, NWB
            read(CON, "(:8X,9F8.0)") (PRFD(J, JW), J = 1, NPRF(JW))
        end do
        read(CON, "(/)")
        do JW = 1, NWB
            read(CON, "(:8X,9F8.0)") (PRFF(J, JW), J = 1, NPRF(JW))
        end do
        read(CON, "(/)")
        do JW = 1, NWB
            read(CON, "(:8X,9I8)") (IPRF(J, JW), J = 1, NIPRF(JW))
        end do
        read(CON, "(//(8X,A8,2I8))") (SPRC(JW), NSPR(JW), NISPR(JW), JW = 1, NWB)
        read(CON, "(/)")
        do JW = 1, NWB
            read(CON, "(:8X,9F8.0)") (SPRD(J, JW), J = 1, NSPR(JW))
        end do
        read(CON, "(/)")
        do JW = 1, NWB
            read(CON, "(:8X,9F8.0)") (SPRF(J, JW), J = 1, NSPR(JW))
        end do
        read(CON, "(/)")
        do JW = 1, NWB
            read(CON, "(:8X,9I8)") (ISPR(J, JW), J = 1, NISPR(JW))
        end do
        read(CON, "(//(8X,A8,I8))") (VPLC(JW), NVPL(JW), JW = 1, NWB)
        read(CON, "(/)")
        do JW = 1, NWB
            read(CON, "(:8X,9F8.0)") (VPLD(J, JW), J = 1, NVPL(JW))
        end do
        read(CON, "(/)")
        do JW = 1, NWB
            read(CON, "(:8X,9F8.0)") (VPLF(J, JW), J = 1, NVPL(JW))
        end do
        read(CON, "(//(8X,A8,I8,A8))") (CPLC(JW), NCPL(JW), TECPLOT(JW), JW = 1, NWB)
        read(CON, "(/)")
        do JW = 1, NWB
            read(CON, "(:8X,9F8.0)") (CPLD(J, JW), J = 1, NCPL(JW))
        end do
        read(CON, "(/)")
        do JW = 1, NWB
            read(CON, "(:8X,9F8.0)") (CPLF(J, JW), J = 1, NCPL(JW))
        end do
        read(CON, "(//(8X,A8,I8))") (FLXC(JW), NFLX(JW), JW = 1, NWB)
        read(CON, "(/)")
        do JW = 1, NWB
            read(CON, "(:8X,9F8.0)") (FLXD(J, JW), J = 1, NFLX(JW))
        end do
        read(CON, "(/)")
        do JW = 1, NWB
            read(CON, "(:8X,9F8.0)") (FLXF(J, JW), J = 1, NFLX(JW))
        end do
        read(CON, "(//8X,A8,2I8)") TSRC, NTSR, NIKTSR;         allocate(ITSR(MAX(1, NIKTSR)), ETSR(MAX(1, NIKTSR)), JBTSR(MAX(1, NIKTSR)))
        read(CON, "(//(:8X,9F8.0))") (TSRD(J), J = 1, NTSR)
        read(CON, "(//(:8X,9F8.0))") (TSRF(J), J = 1, NTSR)
        read(CON, "(//(:8X,9I8))") (ITSR(J), J = 1, NIKTSR)
        read(CON, "(//(:8X,9F8.0))") (ETSR(J), J = 1, NIKTSR)

        read(CON, "(//(8X,A,F8.0))") WLC, WLF
        read(CON, "(//(8X,A,F8.0))") FLOWBALC, FLOWBALF
        read(CON, "(//(8X,A,F8.0))") NPBALC, NPBALF

        read(CON, "(//8X,A8,2I8)") WDOC, NWDO, NIWDO;         allocate(IWDO(MAX(1, NIWDO)))
        read(CON, "(//(:8X,9F8.0))") (WDOD(J), J = 1, NWDO)
        read(CON, "(//(:8X,9F8.0))") (WDOF(J), J = 1, NWDO)
        read(CON, "(//(8X,9I8))") (IWDO(J), J = 1, NIWDO)
        read(CON, "(//8X,A8,I8,A8)") RSOC, NRSO, RSIC;         RSOD = 0.0 ! SW 9/27/13 INITIALIZE SINCE ALLOCATED AS NOD BUT ONLY NRSO USED
        read(CON, "(//(:8X,9F8.0))") (RSOD(J), J = 1, NRSO)
        read(CON, "(//(:8X,9F8.0))") (RSOF(J), J = 1, NRSO)
    else ! w2_con.csv file format

        read(CON, *) (PPUC(JP), JP = 1, NPU);         PPUC = ADJUSTR(PPUC)
        read(CON, *) (ETPU(JP), JP = 1, NPU)
        read(CON, *) (EBPU(JP), JP = 1, NPU)
        read(CON, *) (KTPU(JP), JP = 1, NPU)
        read(CON, *) (KBPU(JP), JP = 1, NPU)
        read(CON, *)
        read(CON, *)

        read(CON, *) (IWR(JW), JW = 1, NIW)
        read(CON, *) (EKTWR(JW), JW = 1, NIW)
        read(CON, *) (EKBWR(JW), JW = 1, NIW)
        read(CON, *)
        read(CON, *)

        read(CON, *) (WDIC(JW), JW = 1, NWD);         WDIC = adjustr(WDIC)
        read(CON, *) (IWD(JW), JW = 1, NWD)
        read(CON, *) (EWD(JW), JW = 1, NWD)
        read(CON, *) (KTWD(JW), JW = 1, NWD)
        read(CON, *) (KBWD(JW), JW = 1, NWD);         TRC = "      " ! SW 9/27/13 INITIALIZATION SINCE ALLOCATION IS TO NTRT
        read(CON, *)
        read(CON, *)

        read(CON, *) (TRC(JT), JT = 1, NTR)
        read(CON, *) (TRIC(JT), JT = 1, NTR);         TRC = ADJUSTR(TRC)
        read(CON, *) (ITR(JT), JT = 1, NTR);         TRIC = ADJUSTR(TRIC)
        read(CON, *) (ELTRT(JT), JT = 1, NTR)
        read(CON, *) (ELTRB(JT), JT = 1, NTR)
        read(CON, *) (QTRFN(JT), JT = 1, NTR)
        read(CON, *) (TTRFN(JT), JT = 1, NTR)
        read(CON, *) (CTRFN(JT), JT = 1, NTR)
        read(CON, *)
        read(CON, *)

        read(CON, *) (DTRC(JB), JB = 1, NBR);         DTRC = adjustr(DTRC)
        read(CON, *)
        read(CON, *)

! Output control cards (excluding constituents)

        do JH = 1, NHY
            read(CON, *) HNAME(JH), FMTH(JH), HMULT(JH), (HPRWBC(JH, JW), JW = 1, NWB);             HPRWBC = ADJUSTR(HPRWBC)
        end do
        read(CON, *)
        read(CON, *)

        read(CON, *) SNPC(1)
        read(CON, *) NSNP(1)
        SNPC(2:NWB) = SNPC(1);         SNPC = ADJUSTR(SNPC)
        NSNP(2:NWB) = NSNP(1)

!  READ (CON,*)        (NISNP(JW), JW=1,NWB)   In contrast to w2_con.npt, all segments are used for SNP output, no user input required.

        NISNP = 0 ! SW 3/31/2020
        do JW = 1, NWB
            do JB = BS(JW), BE(JW)
                do I = US(JB), DS(JB)
                    NISNP(JW) = NISNP(JW) + 1
                    ISNP(NISNP(JW), JW) = I
                end do
            end do
        end do
        read(CON, *) (SNPD(J, 1), J = 1, NSNP(1))
        read(CON, *) (SNPF(J, 1), J = 1, NSNP(1))
        do J = 1, NSNP(1)
            SNPD(J, 2:NWB) = SNPD(J, 1)
            SNPF(J, 2:NWB) = SNPF(J, 1)
        end do

        read(CON, *)
        read(CON, *)
        read(CON, *) SCRC(1);         SCRC(1) = ADJUSTR(SCRC(1))
        read(CON, *) NSCR(1)
        read(CON, *) (SCRD(J, 1), J = 1, NSCR(1))
        read(CON, *) (SCRF(J, 1), J = 1, NSCR(1))
        read(CON, *)
        read(CON, *)

        if (NWB > 1) then
            SCRC(2:NWB) = SCRC(1)
            NSCR(2:NWB) = NSCR(1)
            do J = 1, NSCR(1)
                SCRD(J, 2:NWB) = SCRD(J, 1)
                SCRF(J, 2:NWB) = SCRD(J, 1)
            end do
        end if

        CDUM = "     OFF"
        read(CON, *) CDUM;         CDUM = ADJUSTR(CDUM)
        read(CON, *) NDUM ! PRFC(1), NPRF(1), NIPRF(1)    
        read(CON, *) NIDUM ! PRFC(1), NPRF(1), NIPRF(1)    

        read(CON, *) (DDUM(J), J = 1, NDUM) !(PRFD(J,1),J=1,NPRF(1))
        read(CON, *) (FDUM(J), J = 1, NDUM) !(PRFF(J,1),J=1,NPRF(1))
        read(CON, *) (IDUM(J), J = 1, NIDUM) !(IPRF(J,1),J=1,NIPRF(1))
        read(CON, *)
        read(CON, *)

        PRFC = "     OFF"
        NPRF = 0
        NIPRF = 0
        if (CDUM == "      ON") then !NWB>1 .AND. 
            do J = 1, NDUM
                PRFD(J, 1:NWB) = DDUM(J)
                PRFF(J, 1:NWB) = FDUM(J)
            end do
            NPRF(1:NWB) = NDUM
            do JW = 1, NWB
                JJ = 0
                do J = 1, NIDUM
                    if (IDUM(J) >= US(BS(JW)) .and. IDUM(J) <= DS(BE(JW))) then
                        JJ = JJ + 1
                        IPRF(JJ, JW) = IDUM(J)
                        NIPRF(JW) = JJ
                        PRFC(JW) = "      ON"
                    else
                        if (IDUM(J) == -1) then
                            JJ = JJ + 1
                            IPRF(JJ, JW) = IDUM(J)
                            NIPRF(JW) = JJ
                            PRFC(JW) = "      ON"
                            exit
                        end if
                    end if
                end do
            end do
        end if

        CDUM = "     OFF"
        read(CON, *) CDUM;         CDUM = ADJUSTR(CDUM)
        read(CON, *) NDUM
        read(CON, *) NIDUM
        read(CON, *) (DDUM(J), J = 1, NDUM)
        read(CON, *) (FDUM(J), J = 1, NDUM)
        read(CON, *) (IDUM(J), J = 1, NIDUM)
        read(CON, *)
        read(CON, *)

        SPRC = "     OFF"
        NSPR = 0
        NISPR = 0
        if (CDUM == "      ON" .or. CDUM == "     ONV") then
            do J = 1, NDUM
                SPRD(J, 1:NWB) = DDUM(J)
                SPRF(J, 1:NWB) = FDUM(J)
            end do
            NSPR(1:NWB) = NDUM
            do JW = 1, NWB
                JJ = 0
                do J = 1, NIDUM
                    if (IDUM(J) >= US(BS(JW)) .and. IDUM(J) <= DS(BE(JW))) then
                        JJ = JJ + 1
                        ISPR(JJ, JW) = IDUM(J)
                        NISPR(JW) = JJ
                        SPRC(JW) = CDUM
                    end if
                end do
            end do
        end if

        VPLC = "     OFF"
        NVPL = 0
        read(CON, *) VPLC(1);         VPLC(1) = ADJUSTR(VPLC(1))
        read(CON, *) NVPL(1)
        read(CON, *) (VPLD(J, 1), J = 1, NVPL(1))
        read(CON, *) (VPLF(J, 1), J = 1, NVPL(1))
        read(CON, *)
        read(CON, *)

        read(CON, *) CPLC(1);         CPLC(1) = ADJUSTR(CPLC(1))
        read(CON, *) NCPL(1)
        read(CON, *) TECPLOT(1);         TECPLOT(1) = ADJUSTR(TECPLOT(1))

        CPLC(2:NWB) = CPLC(1)
        NCPL(2:NWB) = NCPL(1)
        TECPLOT(2:NWB) = TECPLOT(1)
        read(CON, *) (CPLD(J, 1), J = 1, NCPL(1))
        read(CON, *) (CPLF(J, 1), J = 1, NCPL(1))
        read(CON, *)
        read(CON, *)
        do J = 1, NCPL(1)
            CPLD(J, 2:NWB) = CPLD(J, 1)
            CPLF(J, 2:NWB) = CPLF(J, 1)
        end do

        read(CON, *) FLXC(1);         FLXC(1) = ADJUSTR(FLXC(1))
        read(CON, *) NFLX(1)
        FLXC(2:NWB) = FLXC(1)
        NFLX(2:NWB) = NFLX(1)

        read(CON, *) (FLXD(J, 1), J = 1, NFLX(1))
        read(CON, *) (FLXF(J, 1), J = 1, NFLX(1))
        read(CON, *)
        read(CON, *)
        do J = 1, NFLX(1)
            FLXD(J, 2:NWB) = FLXD(J, 1)
            FLXF(J, 2:NWB) = FLXF(J, 1)
        end do

        read(CON, *) TSRC;         TSRC = ADJUSTR(TSRC)
        read(CON, *) NTSR
        read(CON, *) NIKTSR
        read(CON, *) TSRFN1

        allocate(ITSR(MAX(1, NIKTSR)), ETSR(MAX(1, NIKTSR)), JBTSR(MAX(1, NIKTSR)))
        read(CON, *) (TSRD(J), J = 1, NTSR)
        read(CON, *) (TSRF(J), J = 1, NTSR)
        read(CON, *) (ITSR(J), J = 1, NIKTSR)
        read(CON, *) (ETSR(J), J = 1, NIKTSR)
        read(CON, *)
        read(CON, *)
        read(CON, *) WLC;         WLC = ADJUSTR(WLC)
        read(CON, *) WLF

        read(CON, *)
        read(CON, *)

        read(CON, *) FLOWBALC;         FLOWBALC = ADJUSTR(FLOWBALC)
        read(CON, *) FLOWBALF

        read(CON, *)
        read(CON, *)

        read(CON, *) NPBALC;         NPBALC = ADJUSTR(NPBALC)
        read(CON, *) NPBALF
        read(CON, *)
        read(CON, *)


        read(CON, *) WDOC;         WDOC = ADJUSTR(WDOC)
        read(CON, *) NWDO
        read(CON, *) NIWDO
        read(CON, *) WDOFN

        allocate(IWDO(MAX(1, NIWDO)))
        read(CON, *) (WDOD(J), J = 1, NWDO)
        read(CON, *) (WDOF(J), J = 1, NWDO)
        read(CON, *) (IWDO(J), J = 1, NIWDO)
        read(CON, *)
        read(CON, *)

        read(CON, *) RSOC;         RSOC = ADJUSTR(RSOC)
        read(CON, *) NRSO
        read(CON, *) RSIC;         RSIC = ADJUSTR(RSIC)
        read(CON, *) RSIFN

        RSOD = 0.0 ! SW 9/27/13 INITIALIZE SINCE ALLOCATED AS NOD BUT ONLY NRSO USED
        read(CON, *) (RSOD(J), J = 1, NRSO)
        read(CON, *) (RSOF(J), J = 1, NRSO)
        read(CON, *)
        read(CON, *)

    end if


! DETERMINE BRANCH FOR EACH TSR FILE   ! SW 7/24/2018
    do J = 1, NIKTSR
        do JB = 1, NBR
            if (ITSR(J) >= US(JB) .and. ITSR(J) <= DS(JB)) then
                JBTSR(J) = JB
                exit
            end if
        end do

    end do

    if (CONFN == "w2_con.npt") then

! Constituent control cards

        read(CON, "(//8X,2A8,I8,F8.0,A8)") CCC, LIMC, CUF, PCO2ATMPPM, CO2YEARLYPPM

        read(CON, "(//8x,A8,A8)") ATM_DEPOSITIONC(1), ATM_DEPOSITION_INTERPOLATION(1)
        do JW = 2, NWB
            read(CON, "(8x,A8,A8)") ATM_DEPOSITIONC(JW), ATM_DEPOSITION_INTERPOLATION(JW)
        end do

        read(CON, "(//(2A8))") (CNAME2(JC), CAC(JC), JC = 1, NCT)
        read(CON, "(/)")

        do JD = 1, NDC
            if (nwb < 10) then
                read(CON, "(A8,(:9A8))") CDNAME2(JD), (CDWBC(JD, JW), JW = 1, NWB)
            end if
            if (nwb >= 10) then
                read(CON, "(A8,9A8,/(:8X,9A8))") CDNAME2(JD), (CDWBC(JD, JW), JW = 1, NWB)
            end if !cb 9/13/12  sw 2/18/13  Foramt 6/16/13 8/13/13
        end do

        read(CON, "(/)")
        do jf = 1, 72
            if (nwb < 10) then
                read(CON, "(A8,(:9A8))") KFNAME2(JF), (KFWBC(JF, JW), JW = 1, NWB)
            end if
            if (nwb >= 10) then
                read(CON, "(A8,9A8,/(:8X,9A8))") KFNAME2(JF), (KFWBC(JF, JW), JW = 1, NWB)
            end if !cb 9/13/12  sw2/18/13  Foramt 6/16/13 8/13/13
            KFNAME2(JF) = KFNAME2(JF)(1:8) // "(kg/d)"
        end do

        read(CON, "(/)")
        do JC = 1, NCT
            read(CON, "(:8X,9F8.0)") (C2I(JC, JW), JW = 1, NWB)
        end do
        read(CON, "(/)")
        do JC = 1, NCT
            read(CON, "(:8X,9A8)") (CPRWBC(JC, JW), JW = 1, NWB)
        end do
        read(CON, "(/)")
        do JC = 1, NCT
            read(CON, "(:8X,9A8)") (CINBRC(JC, JB), JB = 1, NBR)
        end do
        read(CON, "(/)")
        do JC = 1, NCT
            read(CON, "(:8X,9A8)") (CTRTRC(JC, JT), JT = 1, NTR)
        end do
        read(CON, "(/)")
        do JC = 1, NCT
            read(CON, "(:8X,9A8)") (CDTBRC(JC, JB), JB = 1, NBR)
        end do
        read(CON, "(/)")
        do JC = 1, NCT
            read(CON, "(:8X,9A8)") (CPRBRC(JC, JB), JB = 1, NBR)
        end do

! Kinetics coefficients

        read(CON, "(//(8X,4F8.0,2A8))") (EXH2O(JW), EXSS(JW), EXOM(JW), BETA(JW), EXC(JW), EXIC(JW), JW = 1, NWB)
        read(CON, "(//(8X,9F8.0))") (EXA(JA), JA = 1, NAL)
        read(CON, "(//(8X,9F8.0))") (EXZ(JZ), JZ = 1, NZPT)
        read(CON, "(//(8X,9F8.0))") (EXM(JM), JM = 1, NMCT)
        read(CON, "(//(8X,8F8.0))") (CGQ10(JG), CG0DK(JG), CG1DK(JG), CGS(JG), CGLDK(JG), CGKLF(JG), CGCS(JG), CGR(JG), JG = 1, NGC) !LCJ 2/26/15
        read(CON, "(//(8X,F8.0,A,2F8.0))") (SSS(JS), SEDRC(JS), TAUCR(JS), SSCS(JS), JS = 1, NSS) ! READ (CON,'(//(8X,F8.0,A8,2F8.0,I8))') (SSS(JS), SEDRC(JS),  TAUCR(JS),  SSFLOC(JS), FLOCEQN(JS),            JS=1,NSS) !SR 04/21/13

        read(CON, "(//(8X,4F8.0))") (BACTQ10(JW), BACT1DK(JW), BACTS(JW), BACTLDK(JW), JW = 1, NWB)
!READ (CON,'(//(8X,3F8.0))')         (A_DISG(JW), B_DISG(JW), C_DISG(JW),         JW=1,NWB)
        read(CON, "(//(8X,4F8.0))") (H2SR(JW), H2SQ10(JW), H2S1DK(JW), SO4R(JW), JW = 1, NWB)
        read(CON, "(//(8X,3F8.0))") (CH4R(JW), CH4Q10(JW), CH41DK(JW), JW = 1, NWB)
        read(CON, "(//(8X,5F8.0))") (FEIIR(JW), KFE_OXID(JW), KFE_RED(JW), KFEOOH_HalfSat(JW), FeSetVel(JW), JW = 1, NWB)
        read(CON, "(//(8X,5F8.0))") (MNIIR(JW), KMN_OXID(JW), KMN_RED(JW), KMNO2_HalfSat(JW), MnSetVel(JW), JW = 1, NWB)

        read(CON, "(//(8X,9F8.0))") (AG(JA), AR(JA), AE(JA), AM(JA), AS(JA), AHSP(JA), AHSN(JA), AHSSI(JA), ASAT(JA), JA = 1, NAL)
        read(CON, "(//(8X,8F8.0))") (AT1(JA), AT2(JA), AT3(JA), AT4(JA), AK1(JA), AK2(JA), AK3(JA), AK4(JA), JA = 1, NAL)
        read(CON, "(//(8X,6F8.0,I8,F8.0,A8))") (AP(JA), AN(JA), AC(JA), ASI(JA), ACHLA(JA), APOM(JA), ANEQN(JA), ANPR(JA), AVERTM(JA), JA = 1, NAL)
        read(CON, "(//(8X,9A8))") (EPIC(JW, 1), JW = 1, NWB)
        do JE = 2, NEPT
            read(CON, "(8X,9A8)") (EPIC(JW, JE), JW = 1, NWB)
        end do
        read(CON, "(//(8X,9A8))") (EPIPRC(JW, 1), JW = 1, NWB)
        do JE = 2, NEPT
            read(CON, "(8X,9A8)") (EPIPRC(JW, JE), JW = 1, NWB)
        end do
        read(CON, "(//(8X,9F8.0))") (EPICI(JW, 1), JW = 1, NWB)
        do JE = 2, NEPT
            read(CON, "(8X,9F8.0)") (EPICI(JW, JE), JW = 1, NWB)
        end do
        read(CON, "(//(8X,8F8.0))") (EG(JE), ER(JE), EE(JE), EM(JE), EB(JE), EHSP(JE), EHSN(JE), EHSSI(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, "(//(8X,2F8.0,I8,F8.0))") (ESAT(JE), EHS(JE), ENEQN(JE), ENPR(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, "(//(8X,8F8.0))") (ET1(JE), ET2(JE), ET3(JE), ET4(JE), EK1(JE), EK2(JE), EK3(JE), EK4(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, "(//(8X,6F8.0))") (EP(JE), EN(JE), EC(JE), ESI(JE), ECHLA(JE), EPOM(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, "(//(8X,8F8.0))") (ZG(JZ), ZR(JZ), ZM(JZ), ZEFF(JZ), PREFP(JZ), ZOOMIN(JZ), ZS2P(JZ), ZS(JZ), JZ = 1, NZPT)

        read(CON, "(//(8X,8F8.0))") (PREFA(JA, 1), JA = 1, NAL) ! MM 7/13/06
        do JZ = 2, NZPT
            read(CON, "((8X,8F8.0))") (PREFA(JA, JZ), JA = 1, NAL)
        end do
        read(CON, "(//(8X,8F8.0))") (PREFZ(JJZ, 1), JJZ = 1, NZPT)
        do JZ = 2, NZPT
            read(CON, "((8X,8F8.0))") (PREFZ(JJZ, JZ), JJZ = 1, NZPT) ! MM 7/13/06
        end do
        read(CON, "(//(8X,8F8.0))") (ZT1(JZ), ZT2(JZ), ZT3(JZ), ZT4(JZ), ZK1(JZ), ZK2(JZ), ZK3(JZ), ZK4(JZ), JZ = 1, NZPT)
        read(CON, "(//(8X,3F8.0))") (ZP(JZ), ZN(JZ), ZC(JZ), JZ = 1, NZPT)
        read(CON, "(//(8X,9A8))") (MACWBC(JW, 1), JW = 1, NWB)
        do JM = 2, NMCT
            read(CON, "(8X,9A8)") (MACWBC(JW, JM), JW = 1, NWB)
        end do
        read(CON, "(//(8X,9A8))") (MPRWBC(JW, 1), JW = 1, NWB)
        do JM = 2, NMCT
            read(CON, "(8X,9A8)") (MPRWBC(JW, JM), JW = 1, NWB)
        end do
        read(CON, "(//(8X,9F8.0))") (MACWBCI(JW, 1), JW = 1, NWB)
        do JM = 2, NMCT
            read(CON, "(8X,9F8.0)") (MACWBCI(JW, JM), JW = 1, NWB)
        end do
        read(CON, "(//(8X,9F8.0))") (MG(JM), MR(JM), MM(JM), MSAT(JM), MHSP(JM), MHSN(JM), MHSC(JM), MPOM(JM), LRPMAC(JM), JM = 1, NMCT)
        read(CON, "(//(8X,2F8.0))") (PSED(JM), NSED(JM), JM = 1, NMCT)
        read(CON, "(//(8X,2F8.0))") (MBMP(JM), MMAX(JM), JM = 1, NMCT)
        read(CON, "(//(8X,4F8.0))") (CDDRAG(JM), DWV(JM), DWSA(JM), ANORM(JM), JM = 1, NMCT) !CB 6/29/06
        read(CON, "(//(8X,8F8.0))") (MT1(JM), MT2(JM), MT3(JM), MT4(JM), MK1(JM), MK2(JM), MK3(JM), MK4(JM), JM = 1, NMCT)
        read(CON, "(//(8X,3F8.0))") (MP(JM), MN(JM), MC(JM), JM = 1, NMCT)
        if (ORGC_CALC) then
            read(CON, "(//(8X,12F8.0))") (LDOMDK(JW), RDOMDK(JW), LRDDK(JW), LDOMPDK(JW), RDOMPDK(JW), LRDOMPDK(JW), LDOMNDK(JW), RDOMNDK(JW), LRDOMNDK(JW), LDOMCDK(JW), RDOMCDK(JW), LRDOMCDK(JW), JW = 1, NWB)
        else
            read(CON, "(//(8X,3F8.0))") (LDOMDK(JW), RDOMDK(JW), LRDDK(JW), JW = 1, NWB)
            LDOMPDK = LDOMDK;             LDOMNDK = LDOMDK;             LDOMCDK = LDOMDK
            RDOMPDK = RDOMDK;             RDOMNDK = RDOMDK;             RDOMCDK = RDOMDK
            LRDOMPDK = LRDDK;             LRDOMNDK = LRDDK;             LRDOMCDK = LRDDK
        end if
        if (ORGC_CALC) then
            read(CON, "(//(8X,15F8.0))") (LPOMDK(JW), RPOMDK(JW), LRPDK(JW), LPOMHK(JW), RPOMHK(JW), POMS(JW), LPOMPDK(JW), RPOMPDK(JW), LRPOMPDK(JW), LPOMNDK(JW), RPOMNDK(JW), LRPOMNDK(JW), LPOMCDK(JW), RPOMCDK(JW), LRPOMCDK(JW), JW = 1, NWB)
        else
            read(CON, "(//(8X,4F8.0))") (LPOMDK(JW), RPOMDK(JW), LRPDK(JW), POMS(JW), JW = 1, NWB)
            LPOMHK = 0.0;             RPOMHK = 0.0
            LPOMPDK = LPOMDK;             LPOMNDK = LPOMDK;             LPOMCDK = LPOMDK
            RPOMPDK = RPOMDK;             RPOMNDK = RPOMDK;             RPOMCDK = RPOMDK
            LRPOMPDK = LRPDK;             LRPOMNDK = LRPDK;             LRPOMCDK = LRPDK
        end if
        read(CON, "(//(8X,4F8.0))") (ORGP(JW), ORGN(JW), ORGC(JW), ORGSI(JW), JW = 1, NWB)
        read(CON, "(//(8X,4F8.0))") (OMT1(JW), OMT2(JW), OMK1(JW), OMK2(JW), JW = 1, NWB)
        read(CON, "(//(8X,3F8.0))") (CoeffA_Turb(JW), CoeffB_Turb(JW), SECC_PAR(JW), JW = 1, NWB)
        read(CON, "(//(8X,4F8.0))") (KBOD(JB), TBOD(JB), RBOD(JB), CBODS(JB), JB = 1, NBOD)
        read(CON, "(//(8X,3F8.0))") (BODP(JB), BODN(JB), BODC(JB), JB = 1, NBOD)
        read(CON, "(//(8X,2F8.0))") (PO4R(JW), PARTP(JW), JW = 1, NWB)
        read(CON, "(//(8X,3F8.0))") (NH4R(JW), NH4DK(JW), KG_H2O_CONSTANT(JW), JW = 1, NWB)
        read(CON, "(//(8X,4F8.0))") (NH4T1(JW), NH4T2(JW), NH4K1(JW), NH4K2(JW), JW = 1, NWB)
        read(CON, "(//(8X,3F8.0))") (NO3DK(JW), NO3S(JW), FNO3SED(JW), JW = 1, NWB)
        read(CON, "(//(8X,4F8.0))") (NO3T1(JW), NO3T2(JW), NO3K1(JW), NO3K2(JW), JW = 1, NWB)
        read(CON, "(//(8X,4F8.0))") (DSIR(JW), PSIS(JW), PSIDK(JW), PARTSI(JW), JW = 1, NWB)
        read(CON, "(//(8X,F8.0))") (CO2R(JW), JW = 1, NWB)
        read(CON, "(//(8X,2F8.0))") (O2NH4(JW), O2OM(JW), JW = 1, NWB)
        read(CON, "(//(8X,2F8.0))") (O2AR(JA), O2AG(JA), JA = 1, NAL)
        read(CON, "(//(8X,2F8.0))") (O2ER(JE), O2EG(JE), JE = 1, NEPT)
        read(CON, "(//(8X,F8.0))") (O2ZR(JZ), JZ = 1, NZPT)
        read(CON, "(//(8X,2F8.0))") (O2MR(JM), O2MG(JM), JM = 1, NMCT)
        read(CON, "(//(8X,F8.0))") KDO
        if (KDO == 0.0) then
            KDO = 0.01
        end if ! SW 10/24/15 ERROR TRAPPING
        read(CON, "(//(8X,2A8,6F8.0,A8))") (SEDCC(JW), SEDPRC(JW), SEDCI(JW), SDK(JW), SEDS(JW), FSOD(JW), FSED(JW), SEDB(JW), DYNSEDK(JW), JW = 1, NWB) ! cb 11/28/06
        read(CON, "(//(8X,4F8.0))") (SODT1(JW), SODT2(JW), SODK1(JW), SODK2(JW), JW = 1, NWB)
        read(CON, "(//(8X,9F8.0))") (SOD(I), I = 1, IMX)
        read(CON, "(//(8X,A8,I8,4F8.2))") (REAERC(JW), NEQN(JW), RCOEF1(JW), RCOEF2(JW), RCOEF3(JW), RCOEF4(JW), DGPO2(JW), MINKL(JW), JW = 1, NWB)

! Input filenames

        read(CON, "(//(8X,A72))") RSIFN
        read(CON, "(//(8X,A72))") QWDFN
        read(CON, "(//(8X,A72))") QGTFN
        read(CON, "(//(8X,A72))") WSCFN
        read(CON, "(//(8X,A72))") SHDFN
        read(CON, "(//(8X,A72))") (BTHFN(JW), JW = 1, NWB)
        read(CON, "(//(8X,A72))") (METFN(JW), JW = 1, NWB)
        read(CON, "(//(8X,A72))") (EXTFN(JW), JW = 1, NWB)
        read(CON, "(//(8X,A72))") (ATMDEPFN(JW), JW = 1, NWB)
        read(CON, "(//(8X,A72))") (VPRFN(JW), JW = 1, NWB)
        read(CON, "(//(8X,A72))") (LPRFN(JW), JW = 1, NWB)
        read(CON, "(//(8X,A72))") (QINFN(JB), JB = 1, NBR)
        read(CON, "(//(8X,A72))") (TINFN(JB), JB = 1, NBR)
        read(CON, "(//(8X,A72))") (CINFN(JB), JB = 1, NBR)
        read(CON, "(//(8X,A72))") (QOTFN(JB), JB = 1, NBR)
        read(CON, "(//(8X,A72))") (QTRFN(JT), JT = 1, NTR)
        read(CON, "(//(8X,A72))") (TTRFN(JT), JT = 1, NTR)
        read(CON, "(//(8X,A72))") (CTRFN(JT), JT = 1, NTR)
        read(CON, "(//(8X,A72))") (QDTFN(JB), JB = 1, NBR)
        read(CON, "(//(8X,A72))") (TDTFN(JB), JB = 1, NBR)
        read(CON, "(//(8X,A72))") (CDTFN(JB), JB = 1, NBR)
        read(CON, "(//(8X,A72))") (PREFN(JB), JB = 1, NBR)
        read(CON, "(//(8X,A72))") (TPRFN(JB), JB = 1, NBR)
        read(CON, "(//(8X,A72))") (CPRFN(JB), JB = 1, NBR)
        read(CON, "(//(8X,A72))") (EUHFN(JB), JB = 1, NBR)
        read(CON, "(//(8X,A72))") (TUHFN(JB), JB = 1, NBR)
        read(CON, "(//(8X,A72))") (CUHFN(JB), JB = 1, NBR)
        read(CON, "(//(8X,A72))") (EDHFN(JB), JB = 1, NBR)
        read(CON, "(//(8X,A72))") (TDHFN(JB), JB = 1, NBR)
        read(CON, "(//(8X,A72))") (CDHFN(JB), JB = 1, NBR)

! Output filenames

        read(CON, "(//(8X,A72))") (SNPFN(JW), JW = 1, NWB)
        read(CON, "(//(8X,A72))") (PRFFN(JW), JW = 1, NWB)
        read(CON, "(//(8X,A72))") (VPLFN(JW), JW = 1, NWB)
        read(CON, "(//(8X,A72))") (CPLFN(JW), JW = 1, NWB)
        read(CON, "(//(8X,A72))") (SPRFN(JW), JW = 1, NWB)

        read(CON, "(//(8X,A72))") (FLXFN(JW), JW = 1, NWB)
        read(CON, "(//(8X,A72))") TSRFN1
        read(CON, "(//(8X,A72))") WDOFN
        close(CON)
    else

! Constituent control cards

        read(CON, *) CCC, LIMC, CUF, PCO2ATMPPM, CO2YEARLYPPM;         CCC = ADJUSTR(CCC);         LIMC = ADJUSTR(LIMC);         CO2YEARLYPPM = ADJUSTR(CO2YEARLYPPM)
        read(CON, *)
        read(CON, *)
        read(CON, *) (ATM_DEPOSITIONC(JW), JW = 1, NWB);         ATM_DEPOSITIONC = ADJUSTR(ATM_DEPOSITIONC)
        read(CON, *) (ATM_DEPOSITION_INTERPOLATION(JW), JW = 1, NWB);         ATM_DEPOSITION_INTERPOLATION = ADJUSTR(ATM_DEPOSITION_INTERPOLATION)
        read(CON, *)
        read(CON, *)

        do JC = 1, NCT
            read(CON, *) CNAME2(JC), CNAME(JC), CAC(JC), FMTC(JC), CMULT(JC), (C2I(JC, JW), JW = 1, NWB), (CPRWBC(JC, JW), JW = 1, NWB), (C_ATM_DEPOSITION(JC, JW), JW = 1, NWB), (CINBRC(JC, JB), JB = 1, NBR), (CTRTRC(JC, JT), JT = 1, NTR1), (CDTBRC(JC, JB), JB = 1, NBR), (CPRBRC(JC, JB), JB = 1, NBR)
        end do
        CAC = ADJUSTR(CAC);         CPRWBC = ADJUSTR(CPRWBC);         CINBRC = ADJUSTR(CINBRC);         CTRTRC = ADJUSTR(CTRTRC);         CDTBRC = ADJUSTR(CDTBRC);         CPRBRC = ADJUSTR(CPRBRC);         C_ATM_DEPOSITION = ADJUSTR(C_ATM_DEPOSITION)
        read(CON, *)
        read(CON, *)
        do JD = 1, NDC
            read(CON, *) CDNAME2(JD), CDNAME(JD), FMTCD(JD), CDMULT(JD), (CDWBC(JD, JW), JW = 1, NWB)
        end do
        CDWBC = ADJUSTR(CDWBC)
        read(CON, *)
        read(CON, *)

!  DO JF=1,NFL
        do JF = 1, 72 ! THESE ARE THE NUMBER IN THE CONTROL FILE FOR READING
            read(CON, *) KFNAME2(JF), (KFWBC(JF, JW), JW = 1, NWB)
            KFNAME2(JF) = KFNAME2(JF)(1:8) // "(kg/d)"
        end do
        KFWBC = ADJUSTR(KFWBC)
! Kinetics coefficients
        read(CON, *)
        read(CON, *)

        read(CON, *) (EXH2O(JW), JW = 1, NWB)
        read(CON, *) (EXSS(JW), JW = 1, NWB)
        read(CON, *) (EXOM(JW), JW = 1, NWB)
        read(CON, *) (BETA(JW), JW = 1, NWB)
        read(CON, *) (EXC(JW), JW = 1, NWB);         EXC = ADJUSTR(EXC)
        read(CON, *) (EXIC(JW), JW = 1, NWB);         EXIC = ADJUSTR(EXIC)

        read(CON, *)
        read(CON, *)

        read(CON, *) (EXA(JA), JA = 1, NAL)
        read(CON, *)
        read(CON, *)

        read(CON, *) (EXZ(JZ), JZ = 1, NZPT)
        read(CON, *)
        read(CON, *)

        read(CON, *) (EXM(JM), JM = 1, NMCT)
        read(CON, *)
        read(CON, *)

        read(CON, *) (CGQ10(JG), JG = 1, NGC)
        read(CON, *) (CG0DK(JG), JG = 1, NGC)
        read(CON, *) (CG1DK(JG), JG = 1, NGC)
        read(CON, *) (CGS(JG), JG = 1, NGC)
        read(CON, *) (CGLDK(JG), JG = 1, NGC)
        read(CON, *) (CGKLF(JG), JG = 1, NGC)
        read(CON, *) (CGCS(JG), JG = 1, NGC)
        read(CON, *) (CGR(JG), JG = 1, NGC)

        read(CON, *)
        read(CON, *)

        read(CON, *) (SSS(JS), JS = 1, NSS) ! READ (CON,'(//(8X,F8.0,A8,2F8.0,I8))') (SSS(JS), SEDRC(JS),  TAUCR(JS),  SSFLOC(JS), FLOCEQN(JS),            JS=1,NSS) !SR 04/21/13
        read(CON, *) (SEDRC(JS), JS = 1, NSS) ! READ (CON,'(//(8X,F8.0,A8,2F8.0,I8))') (SSS(JS), SEDRC(JS),  TAUCR(JS),  SSFLOC(JS), FLOCEQN(JS),            JS=1,NSS) !SR 04/21/13
        read(CON, *) (TAUCR(JS), JS = 1, NSS) ! READ (CON,'(//(8X,F8.0,A8,2F8.0,I8))') (SSS(JS), SEDRC(JS),  TAUCR(JS),  SSFLOC(JS), FLOCEQN(JS),            JS=1,NSS) !SR 04/21/13
        read(CON, *) (SSCS(JS), JS = 1, NSS)
        SEDRC = ADJUSTR(SEDRC)
        read(CON, *)
        read(CON, *)

        read(CON, *) (BACTQ10(JW), JW = 1, NWB)
        read(CON, *) (BACT1DK(JW), JW = 1, NWB)
        read(CON, *) (BACTLDK(JW), JW = 1, NWB)
        read(CON, *) (BACTS(JW), JW = 1, NWB)
        read(CON, *)
        read(CON, *)

!READ (CON,*)     (A_DISG(JW),  JW=1,NWB)
!READ (CON,*)     (B_DISG(JW),  JW=1,NWB)
!READ (CON,*)     (C_DISG(JW),  JW=1,NWB)
!READ (CON,*)
!READ (CON,*)

        read(CON, *) (H2SR(JW), JW = 1, NWB)
        read(CON, *) (H2SQ10(JW), JW = 1, NWB)
        read(CON, *) (H2S1DK(JW), JW = 1, NWB)
        read(CON, *) (SO4R(JW), JW = 1, NWB)
        read(CON, *)
        read(CON, *)

        read(CON, *) (CH4R(JW), JW = 1, NWB)
        read(CON, *) (CH4Q10(JW), JW = 1, NWB)
        read(CON, *) (CH41DK(JW), JW = 1, NWB)
        read(CON, *)
        read(CON, *)

        read(CON, *) (FEIIR(JW), JW = 1, NWB)
        read(CON, *) (KFE_OXID(JW), JW = 1, NWB)
        read(CON, *) (KFE_RED(JW), JW = 1, NWB)
        read(CON, *) (KFEOOH_HalfSat(JW), JW = 1, NWB)
        read(CON, *) (FeSetVel(JW), JW = 1, NWB)
        read(CON, *)
        read(CON, *)

        read(CON, *) (MNIIR(JW), JW = 1, NWB)
        read(CON, *) (KMN_OXID(JW), JW = 1, NWB)
        read(CON, *) (KMN_RED(JW), JW = 1, NWB)
        read(CON, *) (KMNO2_HalfSat(JW), JW = 1, NWB)
        read(CON, *) (MNSetVel(JW), JW = 1, NWB)
        read(CON, *)
        read(CON, *)

        read(CON, *) (AG(JA), JA = 1, NAL)
        read(CON, *) (AR(JA), JA = 1, NAL)
        read(CON, *) (AE(JA), JA = 1, NAL)
        read(CON, *) (AM(JA), JA = 1, NAL)
        read(CON, *) (AS(JA), JA = 1, NAL)
        read(CON, *) (AHSP(JA), JA = 1, NAL)
        read(CON, *) (AHSN(JA), JA = 1, NAL)
        read(CON, *) (AHSSI(JA), JA = 1, NAL)
        read(CON, *) (ASAT(JA), JA = 1, NAL)

        read(CON, *) (AT1(JA), JA = 1, NAL)
        read(CON, *) (AT2(JA), JA = 1, NAL)
        read(CON, *) (AT3(JA), JA = 1, NAL)
        read(CON, *) (AT4(JA), JA = 1, NAL)
        read(CON, *) (AK1(JA), JA = 1, NAL)
        read(CON, *) (AK2(JA), JA = 1, NAL)
        read(CON, *) (AK3(JA), JA = 1, NAL)
        read(CON, *) (AK4(JA), JA = 1, NAL)

        read(CON, *) (AP(JA), JA = 1, NAL)
        read(CON, *) (AN(JA), JA = 1, NAL)
        read(CON, *) (AC(JA), JA = 1, NAL)
        read(CON, *) (ASI(JA), JA = 1, NAL)
        read(CON, *) (ACHLA(JA), JA = 1, NAL)
        read(CON, *) (APOM(JA), JA = 1, NAL)
        read(CON, *) (ANEQN(JA), JA = 1, NAL)
        read(CON, *) (ANPR(JA), JA = 1, NAL)

        read(CON, *) (O2AR(JA), JA = 1, NAL)
        read(CON, *) (O2AG(JA), JA = 1, NAL)
        read(CON, *) (AVERTM(JA), JA = 1, NAL);         AVERTM = ADJUSTR(AVERTM)
        read(CON, *)
        read(CON, *)

        if (NEPT < 6) then
            NEPTT = 5
        else
            NEPTT = NEPT
        end if

        do JE = 1, NEPTT
            read(CON, *) (EPIC(JW, JE), JW = 1, NWB)
            read(CON, *) (EPIPRC(JW, JE), JW = 1, NWB)
            read(CON, *) (EPICI(JW, JE), JW = 1, NWB)
        end do
        EPIC = ADJUSTR(EPIC);         EPIPRC = ADJUSTR(EPIPRC)

        read(CON, *)
        read(CON, *)

        read(CON, *) (EG(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (ER(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (EE(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (EM(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (EB(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (EHSP(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (EHSN(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (EHSSI(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13  

        read(CON, *) (ESAT(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (EHS(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (ENEQN(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (ENPR(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13  

        read(CON, *) (ET1(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (ET2(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (ET3(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (ET4(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (EK1(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (EK2(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (EK3(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (EK4(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13

        read(CON, *) (EP(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (EN(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (EC(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (ESI(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (ECHLA(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (EPOM(JE), JE = 1, NEPT) !JE=1,NEP)  SW 9/27/13
        read(CON, *) (O2ER(JE), JE = 1, NEPT)
        read(CON, *) (O2EG(JE), JE = 1, NEPT)

        read(CON, *)
        read(CON, *)

        read(CON, *) (ZG(JZ), JZ = 1, NZPT)
        read(CON, *) (ZR(JZ), JZ = 1, NZPT)
        read(CON, *) (ZM(JZ), JZ = 1, NZPT)
        read(CON, *) (ZEFF(JZ), JZ = 1, NZPT)
        read(CON, *) (PREFP(JZ), JZ = 1, NZPT)
        read(CON, *) (ZOOMIN(JZ), JZ = 1, NZPT)
        read(CON, *) (ZS2P(JZ), JZ = 1, NZPT)
        read(CON, *) (ZS(JZ), JZ = 1, NZPT)
        read(CON, *) (ZT1(JZ), JZ = 1, NZPT)
        read(CON, *) (ZT2(JZ), JZ = 1, NZPT)
        read(CON, *) (ZT3(JZ), JZ = 1, NZPT)
        read(CON, *) (ZT4(JZ), JZ = 1, NZPT)
        read(CON, *) (ZK1(JZ), JZ = 1, NZPT)
        read(CON, *) (ZK2(JZ), JZ = 1, NZPT)
        read(CON, *) (ZK3(JZ), JZ = 1, NZPT)
        read(CON, *) (ZK4(JZ), JZ = 1, NZPT)

        read(CON, *) (ZP(JZ), JZ = 1, NZPT)
        read(CON, *) (ZN(JZ), JZ = 1, NZPT)
        read(CON, *) (ZC(JZ), JZ = 1, NZPT)

        read(CON, *) (O2ZR(JZ), JZ = 1, NZPT)

        do JA = 1, NALT
            read(CON, *) (PREFA(JA, JZ), JZ = 1, NZPT)
        end do

        do JZ = 1, NZPTT
            read(CON, *) (PREFZ(JZ, JJZ), JJZ = 1, NZPT)
        end do

        read(CON, *)
        read(CON, *)

        do JM = 1, NMCTT
            read(CON, *) (MACWBC(JW, JM), JW = 1, NWB)
        end do
        do JM = 1, NMCTT
            read(CON, *) (MPRWBC(JW, JM), JW = 1, NWB)
        end do
        do JM = 1, NMCTT
            read(CON, *) (MACWBCI(JW, JM), JW = 1, NWB)
        end do
        MACWBC = ADJUSTR(MACWBC);         MPRWBC = ADJUSTR(MPRWBC)

        read(CON, *)
        read(CON, *)

        read(CON, *) (MG(JM), JM = 1, NMCT)
        read(CON, *) (MR(JM), JM = 1, NMCT)
        read(CON, *) (MM(JM), JM = 1, NMCT)
        read(CON, *) (MSAT(JM), JM = 1, NMCT)
        read(CON, *) (MHSP(JM), JM = 1, NMCT)
        read(CON, *) (MHSN(JM), JM = 1, NMCT)
        read(CON, *) (MHSC(JM), JM = 1, NMCT)
        read(CON, *) (MPOM(JM), JM = 1, NMCT)
        read(CON, *) (LRPMAC(JM), JM = 1, NMCT)

        read(CON, *) (PSED(JM), JM = 1, NMCT)
        read(CON, *) (NSED(JM), JM = 1, NMCT)

        read(CON, *) (MBMP(JM), JM = 1, NMCT)
        read(CON, *) (MMAX(JM), JM = 1, NMCT)
        read(CON, *) (CDDRAG(JM), JM = 1, NMCT) !CB 6/29/06
        read(CON, *) (DWV(JM), JM = 1, NMCT) !CB 6/29/06
        read(CON, *) (DWSA(JM), JM = 1, NMCT) !CB 6/29/06
        read(CON, *) (ANORM(JM), JM = 1, NMCT) !CB 6/29/06  

        read(CON, *) (MT1(JM), JM = 1, NMCT)
        read(CON, *) (MT2(JM), JM = 1, NMCT)
        read(CON, *) (MT3(JM), JM = 1, NMCT)
        read(CON, *) (MT4(JM), JM = 1, NMCT)
        read(CON, *) (MK1(JM), JM = 1, NMCT)
        read(CON, *) (MK2(JM), JM = 1, NMCT)
        read(CON, *) (MK3(JM), JM = 1, NMCT)
        read(CON, *) (MK4(JM), JM = 1, NMCT)

        read(CON, *) (MP(JM), JM = 1, NMCT)
        read(CON, *) (MN(JM), JM = 1, NMCT)
        read(CON, *) (MC(JM), JM = 1, NMCT)

        read(CON, *) (O2MR(JM), JM = 1, NMCT)
        read(CON, *) (O2MG(JM), JM = 1, NMCT)

        read(CON, *)
        read(CON, *)

        if (ORGC_CALC) then
            read(CON, *) (LDOMDK(JW), JW = 1, NWB)
            read(CON, *) (RDOMDK(JW), JW = 1, NWB)
            read(CON, *) (LRDDK(JW), JW = 1, NWB)
!READ (CON,*)         (LDOMPDK(JW), JW=1,NWB)   !************Temporary FIX****net version add these lines  ! SW 12/18/2021
!READ (CON,*)         (RDOMPDK(JW), JW=1,NWB)
!READ (CON,*)         (LRDOMPDK(JW),JW=1,NWB)
!READ (CON,*)         (LDOMNDK(JW), JW=1,NWB)
!READ (CON,*)         (RDOMNDK(JW), JW=1,NWB)
!READ (CON,*)         (LRDOMNDK(JW),JW=1,NWB)
!READ (CON,*)         (LDOMCDK(JW), JW=1,NWB)
!READ (CON,*)         (RDOMCDK(JW), JW=1,NWB)
!READ (CON,*)         (LRDOMCDK(JW),JW=1,NWB)
            LDOMPDK = LDOMDK;             LDOMNDK = LDOMDK;             LDOMCDK = LDOMDK ! remove these once the above lines are read in
            RDOMPDK = RDOMDK;             RDOMNDK = RDOMDK;             RDOMCDK = RDOMDK
            LRDOMPDK = LRDDK;             LRDOMNDK = LRDDK;             LRDOMCDK = LRDDK
        else
            read(CON, *) (LDOMDK(JW), JW = 1, NWB)
            read(CON, *) (RDOMDK(JW), JW = 1, NWB)
            read(CON, *) (LRDDK(JW), JW = 1, NWB)
            LDOMPDK = LDOMDK;             LDOMNDK = LDOMDK;             LDOMCDK = LDOMDK
            RDOMPDK = RDOMDK;             RDOMNDK = RDOMDK;             RDOMCDK = RDOMDK
            LRDOMPDK = LRDDK;             LRDOMNDK = LRDDK;             LRDOMCDK = LRDDK
        end if
        read(CON, *)
        read(CON, *)


        if (ORGC_CALC) then
            read(CON, *) (LPOMDK(JW), JW = 1, NWB)
            read(CON, *) (RPOMDK(JW), JW = 1, NWB)
            read(CON, *) (LRPDK(JW), JW = 1, NWB)
!READ (CON,*)         (LPOMHK(JW),   JW=1,NWB)
!READ (CON,*)         (RPOMHK(JW),   JW=1,NWB)
            read(CON, *) (POMS(JW), JW = 1, NWB)
!READ (CON,*)         (LPOMPDK(JW),  JW=1,NWB)     !************TEMPORARY fix SW 12/18/2021
!READ (CON,*)         (RPOMPDK(JW),  JW=1,NWB)
!READ (CON,*)         (LRPOMPDK(JW), JW=1,NWB)
!READ (CON,*)         (LPOMNDK(JW),  JW=1,NWB)
!READ (CON,*)         (RPOMNDK(JW),  JW=1,NWB)
!READ (CON,*)         (LRPOMNDK(JW), JW=1,NWB)
!READ (CON,*)         (LPOMCDK(JW),  JW=1,NWB)
!READ (CON,*)         (RPOMCDK(JW),  JW=1,NWB)
!READ (CON,*)         (LRPOMCDK(JW), JW=1,NWB) 
            LPOMHK = 0.0;             RPOMHK = 0.0 ! remove once the above lines are added back
            LPOMPDK = LPOMDK;             LPOMNDK = LPOMDK;             LPOMCDK = LPOMDK
            RPOMPDK = RPOMDK;             RPOMNDK = RPOMDK;             RPOMCDK = RPOMDK
            LRPOMPDK = LRPDK;             LRPOMNDK = LRPDK;             LRPOMCDK = LRPDK
        else
            read(CON, *) (LPOMDK(JW), JW = 1, NWB)
            read(CON, *) (RPOMDK(JW), JW = 1, NWB)
            read(CON, *) (LRPDK(JW), JW = 1, NWB)
            read(CON, *) (POMS(JW), JW = 1, NWB)
            LPOMHK = 0.0;             RPOMHK = 0.0
            LPOMPDK = LPOMDK;             LPOMNDK = LPOMDK;             LPOMCDK = LPOMDK
            RPOMPDK = RPOMDK;             RPOMNDK = RPOMDK;             RPOMCDK = RPOMDK
            LRPOMPDK = LRPDK;             LRPOMNDK = LRPDK;             LRPOMCDK = LRPDK
        end if

        read(CON, *)
        read(CON, *)

        read(CON, *) (ORGP(JW), JW = 1, NWB)
        read(CON, *) (ORGN(JW), JW = 1, NWB)
        read(CON, *) (ORGC(JW), JW = 1, NWB)
        read(CON, *) (ORGSI(JW), JW = 1, NWB)
        read(CON, *) (O2OM(JW), JW = 1, NWB)

        read(CON, *) (OMT1(JW), JW = 1, NWB)
        read(CON, *) (OMT2(JW), JW = 1, NWB)
        read(CON, *) (OMK1(JW), JW = 1, NWB)
        read(CON, *) (OMK2(JW), JW = 1, NWB)

        read(CON, *)
        read(CON, *)
        read(CON, *) (CoeffA_Turb(JW), JW = 1, NWB)
        read(CON, *) (CoeffB_Turb(JW), JW = 1, NWB)
        read(CON, *) (SECC_PAR(JW), JW = 1, NWB)
        read(CON, *)
        read(CON, *)

        read(CON, *) (KBOD(JB), JB = 1, NBOD)
        read(CON, *) (TBOD(JB), JB = 1, NBOD)
        read(CON, *) (RBOD(JB), JB = 1, NBOD)
        read(CON, *) (CBODS(JB), JB = 1, NBOD)
        read(CON, *) (BODP(JB), JB = 1, NBOD)
        read(CON, *) (BODN(JB), JB = 1, NBOD)
        read(CON, *) (BODC(JB), JB = 1, NBOD)

        read(CON, *)
        read(CON, *)

        read(CON, *) (PO4R(JW), JW = 1, NWB)
        read(CON, *) (PARTP(JW), JW = 1, NWB)
        read(CON, *) (NH4R(JW), JW = 1, NWB)
        read(CON, *) (NH4DK(JW), JW = 1, NWB)
        read(CON, *) (NH4T1(JW), JW = 1, NWB)
        read(CON, *) (NH4T2(JW), JW = 1, NWB)
        read(CON, *) (NH4K1(JW), JW = 1, NWB)
        read(CON, *) (NH4K2(JW), JW = 1, NWB)
        read(CON, *) (KG_H2O_CONSTANT(JW), JW = 1, NWB)
        read(CON, *) (O2NH4(JW), JW = 1, NWB)
        read(CON, *) (NO3DK(JW), JW = 1, NWB)
        read(CON, *) (NO3S(JW), JW = 1, NWB)
        read(CON, *) (FNO3SED(JW), JW = 1, NWB)
        read(CON, *) (NO3T1(JW), JW = 1, NWB)
        read(CON, *) (NO3T2(JW), JW = 1, NWB)
        read(CON, *) (NO3K1(JW), JW = 1, NWB)
        read(CON, *) (NO3K2(JW), JW = 1, NWB)
        read(CON, *) (DSIR(JW), JW = 1, NWB)
        read(CON, *) (PSIS(JW), JW = 1, NWB)
        read(CON, *) (PSIDK(JW), JW = 1, NWB)
        read(CON, *) (PARTSI(JW), JW = 1, NWB)

        read(CON, *)
        read(CON, *)

        read(CON, *) (CO2R(JW), JW = 1, NWB)
        read(CON, *)
        read(CON, *)

        read(CON, *) KDO
        if (KDO == 0.0) then
            KDO = 0.01
        end if ! SW 10/24/15 ERROR TRAPPING

        read(CON, *)
        read(CON, *)

        read(CON, *) (SEDCC(JW), JW = 1, NWB);         SEDCC = adjustr(SEDCC)
        read(CON, *) (SEDPRC(JW), JW = 1, NWB);         SEDPRC = adjustr(SEDPRC)
        read(CON, *) (SEDCI(JW), JW = 1, NWB)
        read(CON, *) (SDK(JW), JW = 1, NWB)
        read(CON, *) (SEDS(JW), JW = 1, NWB)
        read(CON, *) (FSOD(JW), JW = 1, NWB)
        read(CON, *) (FSED(JW), JW = 1, NWB)
        read(CON, *) (SEDB(JW), JW = 1, NWB)
        read(CON, *) (DYNSEDK(JW), JW = 1, NWB);         DYNSEDK = adjustr(DYNSEDK)
        read(CON, *) (SODT1(JW), JW = 1, NWB)
        read(CON, *) (SODT2(JW), JW = 1, NWB)
        read(CON, *) (SODK1(JW), JW = 1, NWB)
        read(CON, *) (SODK2(JW), JW = 1, NWB)
        read(CON, *)
        read(CON, *)

        read(CON, *) (SOD(I), I = 1, IMX)

        read(CON, *)
        read(CON, *)

        read(CON, *) (REAERC(JW), JW = 1, NWB);         REAERC = adjustr(REAERC)
        read(CON, *) (NEQN(JW), JW = 1, NWB)
        read(CON, *) (RCOEF1(JW), JW = 1, NWB)
        read(CON, *) (RCOEF2(JW), JW = 1, NWB)
        read(CON, *) (RCOEF3(JW), JW = 1, NWB)
        read(CON, *) (RCOEF4(JW), JW = 1, NWB)
        read(CON, *) (DGPO2(JW), JW = 1, NWB)
        read(CON, *) (MINKL(JW), JW = 1, NWB)
        read(CON, *)
        read(CON, *)

! Input filenames

        read(CON, *) QWDFN
        read(CON, *) QGTFN
        read(CON, *) WSCFN
        read(CON, *) SHDFN
        read(CON, *) VPLFN(1)
        VPLFN(2:NWB) = VPLFN(1)

        read(CON, *)
        read(CON, *)

        read(CON, *) (BTHFN(JW), JW = 1, NWB)
        read(CON, *) (METFN(JW), JW = 1, NWB)
        read(CON, *) (EXTFN(JW), JW = 1, NWB)
        read(CON, *) (ATMDEPFN(JW), JW = 1, NWB)
        read(CON, *) (VPRFN(JW), JW = 1, NWB)
        read(CON, *) (LPRFN(JW), JW = 1, NWB)

! Output filenames

        read(CON, *) (SNPFN(JW), JW = 1, NWB)
        read(CON, *) (PRFFN(JW), JW = 1, NWB)
        read(CON, *) (CPLFN(JW), JW = 1, NWB)
        read(CON, *) (SPRFN(JW), JW = 1, NWB)
        read(CON, *) (FLXFN(JW), JW = 1, NWB)

        read(CON, *)
        read(CON, *)

        read(CON, *) (QINFN(JB), JB = 1, NBR)
        read(CON, *) (TINFN(JB), JB = 1, NBR)
        read(CON, *) (CINFN(JB), JB = 1, NBR)
        read(CON, *) (QOTFN(JB), JB = 1, NBR)
        read(CON, *) (QDTFN(JB), JB = 1, NBR)
        read(CON, *) (TDTFN(JB), JB = 1, NBR)
        read(CON, *) (CDTFN(JB), JB = 1, NBR)
        read(CON, *) (PREFN(JB), JB = 1, NBR)
        read(CON, *) (TPRFN(JB), JB = 1, NBR)
        read(CON, *) (CPRFN(JB), JB = 1, NBR)
        read(CON, *) (EUHFN(JB), JB = 1, NBR)
        read(CON, *) (TUHFN(JB), JB = 1, NBR)
        read(CON, *) (CUHFN(JB), JB = 1, NBR)
        read(CON, *) (EDHFN(JB), JB = 1, NBR)
        read(CON, *) (TDHFN(JB), JB = 1, NBR)
        read(CON, *) (CDHFN(JB), JB = 1, NBR)


        close(CON)

    end if

    do JW = 1, NWB ! SW 9/28/2018
        if (SPRC(JW) == "     ONV") then
            do N = 1, 70
                if (SPRFN(JW)(N:N) == ".") then
                    SPRVFN(JW) = SPRFN(JW)(1:N - 1) // "_volw.csv"
                    exit
                end if
            end do
        end if
    end do

    KFNAME2(73) = "ALDOMPM(kg/d)" !'LDOM P algal mortality - source, kg/day      '           ! 1-72 are defined in the control file
    KFNAME2(74) = "ELDOMPM(kg/d)" !'LDOM P epiphyton mortality - source, kg/day  '
    KFNAME2(75) = "ALPOMP(kg/d)" !'LPOM P algal production- source, kg/day      '
    KFNAME2(76) = "LPOMPSet(kg/d)" !'LPOM P net settling - source/sink, kg/day    '
    KFNAME2(77) = "RPOMPSet(kg/d)" ! 'RPOM P net settling - source/sink, kg/day    '
    KFNAME2(78) = "LDOMPAM(kg/d)" ! 'LDOM P algal mortality - source, kg/day      '
    KFNAME2(79) = "LDOMPEM(kg/d)" ! 'LDOM P epiphyton mortality - source, kg/day  '
    KFNAME2(80) = "LPOMPAP(kg/d)" ! 'LPOM P algal production- source, kg/day      '
    KFNAME2(81) = "LPOMPNetSet(kg/d)" ! 'LPOM P net settling - source/sink, kg/day    '
    KFNAME2(82) = "RPOMPNetSet(kg/d)" ! 'RPOM P net settling - source/sink, kg/day    '
    KFNAME2(83) = "SedPD(kg/d)" !'Sediment P decay - sink, kg/day              '
    KFNAME2(84) = "SedASet(kg/d)" !'Sediment algal P settling - source, kg/day   '
    KFNAME2(85) = "SedLPOMPSet(kg/d)" !'Sediment P LPOM settling - source,kg/day     '
    KFNAME2(86) = "SedNetPSet(kg/d)" !'Sediment net P settling - source/sink, kg/day'
    KFNAME2(87) = "SedEpiPSett(kg/d)" !'Sediment epiphyton P settling - source,kg/day'
    KFNAME2(88) = "SedND(kg/d)" !'Sediment N decay - sink, kg/day              '
    KFNAME2(89) = "SedANSett(kg/d)" !'Sediment algal N settling - source, kg/day   '
    KFNAME2(90) = "SedLPOMNSett(kg/d)" !'Sediment N LPOM settling - source,kg/day     '
    KFNAME2(91) = "SedNetNSett(kg/d)" !'Sediment net N settling - source/sink, kg/day'
    KFNAME2(92) = "SedEpiNSett(kg/d)" !'Sediment epiphyton N settling - source,kg/day'
    KFNAME2(93) = "SedCD(kg/d)" !'Sediment C decay - sink, kg/day              '
    KFNAME2(94) = "SedACSett(kg/d)" !'Sediment algal C settling - source, kg/day   '
    KFNAME2(95) = "SedLPOMCSett(kg/d)" !'Sediment C LPOM settling - source,kg/day     '
    KFNAME2(96) = "SedNetCSett(kg/d)" !'Sediment net C settling - source/sink, kg/day'
    KFNAME2(97) = "SedEpiCSett(kg/d)" !'Sediment epiphyton C settling - source,kg/day'
    KFNAME2(98) = "SedDeN(kg/d)" !'Sediment N denitrification - source, kg/day  '
    KFNAME2(99) = "PO4MacR(kg/d)" !'PO4 macrophyte resp - source, kg/day         '
    KFNAME2(100) = "PO4MacG(kg/d)" !'PO4 macrophyte growth - sink, kg/day         '
    KFNAME2(101) = "NH4MacR(kg/d)" ! 'NH4 macrophyte resp - source, kg/day         '
    KFNAME2(102) = "NH4MacG(kg/d)" !'NH4 macrophyte growth - sink, kg/day         '
    KFNAME2(103) = "LDOMMacM(kg/d)" !'LDOM macrophyte mort  - source, kg/day       '
    KFNAME2(104) = "LDOMMacM(kg/d)" !'LPOM macrophyte mort  - source, kg/day       '
    KFNAME2(105) = "RPOMMacM(kg/d)" !'RPOM macrophyte mort  - source, kg/day       '
    KFNAME2(106) = "DOMacP(kg/d)" !'DO  macrophyte production  - source, kg/day  '
    KFNAME2(107) = "DOMacR(kg/d)" !'DO  macrophyte respiration - sink, kg/day    '
    KFNAME2(108) = "TICMacG(kg/d)" !'TIC macrophyte growth/resp  - S/S, kg/day    '
    KFNAME2(109) = "CBODS(kg/d)" !'CBOD settling - sink, kg/day                 '
    KFNAME2(110) = "SEDCBOD(kg/d)" !'Sediment CBOD settling - source, kg/day      '
    KFNAME2(111) = "SEDCBODP(kg/d)" !'Sediment CBOD P settling - source, kg/day    '
    KFNAME2(112) = "SEDCBODN(kg/d)" !'Sediment CBOD N settling - source, kg/day    '
    KFNAME2(113) = "SEDCBODC(kg/d)" !'Sediment CBOD C settling - source, kg/day    '
    KFNAME2(114) = "SEDB(kg/d)" !'Sediment Burial - sink, kg/day               '
    KFNAME2(KF_SED_PBURIAL) = "SEDPB(kg/d)" !'Sediment P Burial - sink, kg/day             '
    KFNAME2(KF_SED_NBURIAL) = "SEDNB(kg/d)" !'Sediment N Burial - sink, kg/day             '
    KFNAME2(117) = "SEDCB(kg/d)" !'Sediment C Burial - sink, kg/day             '
    KFNAME2(118) = "CBODPS(kg/d)" !'CBOD P settling - sink, kg/day               '
    KFNAME2(119) = "CBODNS(kg/d)" !'CBOD N settling - sink, kg/day               '
    KFNAME2(KF_CO2X) = "CO2GASX(kg/d)"
    KFNAME2(KF_DOH2S) = "DOH2S(kg/d)"
    KFNAME2(122) = "H2SGASX(kg/d)"
    KFNAME2(123) = "H2SDK(kg/d)"
    KFNAME2(124) = "H2SSOD(kg/d)"
    KFNAME2(KF_DOCH4) = "DOCH4(kg/d)"
    KFNAME2(126) = "CH4GASX(kg/d)"
    KFNAME2(127) = "CH4DK(kg/d)"
    KFNAME2(128) = "CH4SOD((kg/d)"
    KFNAME2(KF_FE2D) = "Fe2D(kg/d)"
    KFNAME2(130) = "DOFe2(kg/d)"
    KFNAME2(131) = "FEIISOD(kg/d)" !'FeOOH settling from water col. - sink, kg/day'
    KFNAME2(132) = "SDINFeOOH(kg/d)" !'Settling of FeOOH into water layer, kg/d     '
    KFNAME2(KF_MN2D) = "Mn2d(kg/d)"
    KFNAME2(134) = "DOMn2(kg/d)"
    KFNAME2(135) = "MNIISOD(kg/d)" !'MnO2 settling from water col. - sink, kg/day '
    KFNAME2(136) = "SDINMnO2(kg/d)" !  'Settling of MnO2  into water layer, kg/d     '
    KFNAME2(KF_SDINC) = "SD_C_IN(kg/d)"
    KFNAME2(138) = "SD_N_IN(kg/d)"
    KFNAME2(139) = "SD_P_IN(kg/d)"
    KFNAME2(140) = "DOSEDIA(kg/d)"

    KFNAME2(KF_SEDD) = "SEDD1(kg/d)" ! 'Labile standing biomass decay- sink, kg/day  '    ! OPTIONAL VARIABLE FOR STANDING ORGANIC MATTER LIKE TREES IN A WATER COLUMN
    KFNAME2(142) = "SEDD2(kg/d)" !'Refract. stand. biomass decay- sink, kg/day  '

! INITIALIZE WATER QUALITY

    ALGAE_TOXIN = .false.
    ALGAE_SETTLING_EXIST = .false.
    CONSTITUENTS = CCC == "      ON"
    if (CONSTITUENTS) then
        do J = NATS, NATE
            if (CAC(J) == "      ON") then
                ALGAE_TOXIN = .true.
                exit
            end if
        end do
        do JA = 1, NAL
            if (AVERTM(JA) == "      ON") then
                ALGAE_SETTLING_EXIST = .true.
                exit
            end if
        end do

        call KINETICS()
    end if

! Bathymetry file

    do JW = 1, NWB
        open(BTH(JW), FILE=BTHFN(JW), STATUS="OLD")
        read(BTH(JW), "(a1)") char1 ! New Bathymetry format option SW 6/22/09
        if (CHAR1 == "$") then
            read(BTH(JW), *)
            read(BTH(JW), *) AID, (DLX(I), I = US(BS(JW)) - 1, DS(BE(JW)) + 1)
            read(BTH(JW), *) AID, (ELWS(I), I = US(BS(JW)) - 1, DS(BE(JW)) + 1)
            read(BTH(JW), *) AID, (PHI0(I), I = US(BS(JW)) - 1, DS(BE(JW)) + 1)
            read(BTH(JW), *) AID, (FRIC(I), I = US(BS(JW)) - 1, DS(BE(JW)) + 1)
            read(BTH(JW), *)
            do K = 1, KMX
                read(BTH(JW), *) H(K, JW), (B(K, I), I = US(BS(JW)) - 1, DS(BE(JW)) + 1)
            end do
            do I = US(BS(JW)) - 1, DS(BE(JW)) + 1
                H2(:, I) = H(:, JW)
            end do
        else
            read(BTH(JW), "(//(10F8.0))") (DLX(I), I = US(BS(JW)) - 1, DS(BE(JW)) + 1)
            read(BTH(JW), "(//(10F8.0))") (ELWS(I), I = US(BS(JW)) - 1, DS(BE(JW)) + 1)
            read(BTH(JW), "(//(10F8.0))") (PHI0(I), I = US(BS(JW)) - 1, DS(BE(JW)) + 1)
            read(BTH(JW), "(//(10F8.0))") (FRIC(I), I = US(BS(JW)) - 1, DS(BE(JW)) + 1)
            read(BTH(JW), "(//(10F8.0))") (H(K, JW), K = 1, KMX)
            do I = US(BS(JW)) - 1, DS(BE(JW)) + 1
                read(BTH(JW), "(//(10F8.0))") (B(K, I), K = 1, KMX)
                H2(:, I) = H(:, JW)
            end do
        end if
! Set water surface of inactive segments to those active cells next to them   SW 8/6/2018
        do JB = BS(JW), BE(JW)
            ELWS(US(JB) - 1) = ELWS(US(JB))
            ELWS(DS(JB) + 1) = ELWS(DS(JB))
            DLX(US(JB) - 1) = DLX(US(JB))
            DLX(DS(JB) + 1) = DLX(DS(JB))
            PHI0(US(JB) - 1) = PHI0(US(JB))
            PHI0(DS(JB) + 1) = PHI0(DS(JB))
            FRIC(US(JB) - 1) = FRIC(US(JB))
            FRIC(DS(JB) + 1) = FRIC(DS(JB))
        end do
!      
        close(BTH(JW))
    end do
    H1 = H2
    BI = B

    allocate(BSAVE(KMX, IMX))
    BSAVE = 0.0
    BSAVE = B

!  Amaila start - reading additional sediment compartments coefficients
    STANDING_BIOMASS_DECAY = .false.
    SEDPRC1 = "     OFF"
    SEDPRC2 = "     OFF"
    inquire(FILE="w2_amaila.npt", EXIST=STANDING_BIOMASS_DECAY) ! SW 4/30/15
    if (STANDING_BIOMASS_DECAY) then
        open(NUNIT, file="w2_amaila.npt", status="old")
        read(NUNIT, "(//(8X,2A8,3F8.0))") (SEDCC1(JW), SEDPRC1(JW), SEDCI1(JW), SDK1(JW), FSEDC1(JW), JW = 1, NWB) ! cb 6/7/17
        read(NUNIT, "(//(8X,2A8,3F8.0))") (SEDCC2(JW), SEDPRC2(JW), SEDCI2(JW), SDK2(JW), FSEDC2(JW), JW = 1, NWB)
        read(NUNIT, "(//(8X,3F8.0))") (pbiom(JW), nbiom(JW), cbiom(JW), JW = 1, NWB) ! cb 6/7/17
        close(NUNIT)
    end if

! SYSTDG INPUT FILE

    SYSTDG = .false.
    N2BND = .false.
    DOBND = .false.
    DGPBND = .false.
    TDGTA = .false.
    inquire(FILE="w2_systdg.npt", EXIST=SYSTDG)
    if (SYSTDG) then
        CONTDG = NUNIT
        open(CONTDG, FILE="w2_systdg.npt", STATUS="OLD")
        call INPUT_SYSTDG()
        SYSTDG = SYSTDGC == "      ON"
        N2BND = N2BNDC == "      ON"
        DOBND = DOBNDC == "      ON"
        DGPBND = TDG2BNDC == "      ON"
        TDGTA = TDGTAC == "      ON"
    end if

! End SYSTDG

! Output file unit numbers

    allocate(TSR(NIKTSR))
    allocate(WDO(NIWDO, 4), WDO2(NWD + NST + NGT + NSP + NPU + NPI, 4))
    do J = 1, 8 ! SW 9/28/2018
        do JW = 1, NWB
            OPT(JW, J) = NUNIT;             NUNIT = NUNIT + 1
        end do
    end do
    do J = 1, NIKTSR
        TSR(J) = NUNIT;         NUNIT = NUNIT + 1
    end do
    do JW = 1, NIWDO
        WDO(JW, 1) = NUNIT;         NUNIT = NUNIT + 1
        WDO(JW, 2) = NUNIT;         NUNIT = NUNIT + 1
        WDO(JW, 3) = NUNIT;         NUNIT = NUNIT + 1
        WDO(JW, 4) = NUNIT;         NUNIT = NUNIT + 1
    end do

! BIOENERGETICS bioexp mlm output filenumber assigment
    if (FISHBIO) then
        do J = 1, NIBIO
            BIOEXPFN(J) = NUNIT;             NUNIT = NUNIT + 1
        end do
        do J = 1, NIBIO
            WEIGHTNUM(J) = NUNIT;             NUNIT = NUNIT + 1
        end do
    end if
! Variable names, formats, multipliers, and Compaq Visual FORTRAN array viewer controls
    if (CONFN == "w2_con.npt") then
        open(GRF, FILE="graph.npt", STATUS="OLD")
        read(GRF, "(///(A43,1X,A9,3F8.0,A8))") (HNAME(J), FMTH(J), HMULT(J), HYMIN(J), HYMAX(J), HPLTC(J), J = 1, NHY)
        read(GRF, "(// (A43,1X,A9,3F8.0,A8))") (CNAME(J), FMTC(J), CMULT(J), CMIN(J), CMAX(J), CPLTC(J), J = 1, NCT)
        read(GRF, "(// (A43,1X,A9,3F8.0,A8))") (CDNAME(J), FMTCD(J), CDMULT(J), CDMIN(J), CDMAX(J), CDPLTC(J), J = 1, NDC) ! SW 10/20/15 INTERNAL TDG
        close(GRF)
    end if
!CDNAME(NDC)='TDG(%)'      ! SW 10/17/15
!FMTCD(NDC)=' (F10.3)'
!CDMULT(NDC)=1.0


    do JC = 1, NCT
        L3 = 1
        L1 = SCAN(CNAME(JC), ",") + 2
        if (L1 == 2) then
            L1 = 43
        end if ! SW 12/3/2012   Implies no comma found
        L2 = SCAN(CNAME(JC)(L1:43), "  ") + L1
        if (L2 > 43) then
            L2 = 43
        end if ! SW 12/3/2012
        CUNIT(JC) = CNAME(JC)(L1:L2)
        CNAME1(JC) = CNAME(JC)(1:L1 - 3)
        CNAME3(JC) = CNAME1(JC)
        do while (L3 < L1 - 3)
            if (CNAME(JC)(L3:L3) == " ") then
                CNAME3(JC)(L3:L3) = "_"
            end if
            L3 = L3 + 1
        end do
        CUNIT1(JC) = CUNIT(JC)(1:1)
        CUNIT2(JC) = CUNIT(JC)
        if (CUNIT(JC)(1:2) == "mg") then
            CUNIT1(JC) = "g"
            CUNIT2(JC) = "g/m^3"
        end if
        if (CUNIT(JC)(1:2) /= "g/" .and. CUNIT(JC)(1:2) /= "mg") then
            CUNIT1(JC) = "  "
        end if
    end do
    do JC = 1, NDC
        L1 = 1
        L2 = MAX(4, SCAN(CDNAME(JC), ",") - 1)
        CDNAME3(JC) = CDNAME(JC)(1:L2)
        do while (L1 < L2)
            if (CDNAME(JC)(L1:L1) == " ") then
                CDNAME3(JC)(L1:L1) = "_"
            end if
            L1 = L1 + 1
        end do
    end do
    FMTH(1:NHY) = ADJUSTL(FMTH(1:NHY))

! Initialize logical control variables

    VERT_PROFILE = .false.
    LONG_PROFILE = .false.
    do JW = 1, NWB
        ISO_TEMP(JW) = T2I(JW) >= 0
        VERT_TEMP(JW) = T2I(JW) == -1
        LONG_TEMP(JW) = T2I(JW) < -1
        if (CONSTITUENTS) then ! CB 12/04/08
            ISO_SEDIMENT(JW) = SEDCI(JW) >= 0 .and. SEDCC(JW) == "      ON"
            VERT_SEDIMENT(JW) = SEDCI(JW) == -1.0 .and. SEDCC(JW) == "      ON"
            LONG_SEDIMENT(JW) = SEDCI(JW) < -1.0 .and. SEDCC(JW) == "      ON"
! Amaila Start
            if (STANDING_BIOMASS_DECAY) then
                ISO_SEDIMENT1(JW) = SEDCI1(JW) >= 0 .and. SEDCC1(JW) == "      ON"
                VERT_SEDIMENT1(JW) = SEDCI1(JW) == -1.0 .and. SEDCC1(JW) == "      ON"
                LONG_SEDIMENT1(JW) = SEDCI1(JW) < -1.0 .and. SEDCC1(JW) == "      ON"
                ISO_SEDIMENT2(JW) = SEDCI2(JW) >= 0 .and. SEDCC2(JW) == "      ON"
                VERT_SEDIMENT2(JW) = SEDCI2(JW) == -1.0 .and. SEDCC2(JW) == "      ON"
                LONG_SEDIMENT2(JW) = SEDCI2(JW) < -1.0 .and. SEDCC2(JW) == "      ON"
            end if
! Amaila End
            ISO_EPIPHYTON(JW, :) = EPICI(JW, :) >= 0 .and. EPIC(JW, :) == "      ON"
            VERT_EPIPHYTON(JW, :) = EPICI(JW, :) == -1.0 .and. EPIC(JW, :) == "      ON"
            LONG_EPIPHYTON(JW, :) = EPICI(JW, :) < -1.0 .and. EPIC(JW, :) == "      ON"
            iso_macrophyte(JW, :) = MACWBCI(JW, :) >= 0 .and. MACWBC(JW, :) == "      ON" ! cb 8/21/15
            vert_macrophyte(JW, :) = MACWBCI(JW, :) == -1.0 .and. MACWBC(JW, :) == "      ON" ! cb 8/21/15
            long_macrophyte(JW, :) = MACWBCI(JW, :) < -1.0 .and. MACWBC(JW, :) == "      ON" ! cb 8/21/15
            do JC = 1, NCT
                ISO_CONC(JC, JW) = C2I(JC, JW) >= 0.0
                VERT_CONC(JC, JW) = C2I(JC, JW) == -1.0 .and. CAC(JC) == "      ON"
                LONG_CONC(JC, JW) = C2I(JC, JW) < -1.0 .and. CAC(JC) == "      ON"
                if (VERT_CONC(JC, JW)) then
                    VERT_PROFILE(JW) = .true.
                end if
                if (LONG_CONC(JC, JW)) then
                    LONG_PROFILE(JW) = .true.
                end if
            end do
            if (VERT_SEDIMENT(JW)) then
                VERT_PROFILE(JW) = .true.
            end if
            if (VERT_SEDIMENT1(JW)) then
                VERT_PROFILE(JW) = .true.
            end if ! amaila
            if (VERT_SEDIMENT2(JW)) then
                VERT_PROFILE(JW) = .true.
            end if ! amaila
            if (LONG_SEDIMENT(JW)) then
                LONG_PROFILE(JW) = .true.
            end if
            if (LONG_SEDIMENT1(JW)) then
                LONG_PROFILE(JW) = .true.
            end if ! amaila
            if (LONG_SEDIMENT2(JW)) then
                LONG_PROFILE(JW) = .true.
            end if ! amaila
            if (ANY(VERT_EPIPHYTON(JW, :))) then
                VERT_PROFILE(JW) = .true.
            end if
            if (ANY(LONG_EPIPHYTON(JW, :))) then
                LONG_PROFILE(JW) = .true.
            end if
            if (ANY(VERT_macrophyte(JW, :))) then
                VERT_PROFILE(JW) = .true.
            end if ! cb 8/21/15
            if (ANY(LONG_macrophyte(JW, :))) then
                LONG_PROFILE(JW) = .true.
            end if ! cb 8/21/15
        end if ! cb 12/04/08
        if (VERT_TEMP(JW)) then
            VERT_PROFILE(JW) = .true.
        end if
        if (LONG_TEMP(JW)) then
            LONG_PROFILE(JW) = .true.
        end if
        do M = 1, NMC
!MACROPHYTE_CALC(JW,M) = CONSTITUENTS.AND.MACWBC(JW,M).EQ.' ON'
!PRINT_MACROPHYTE(JW,M) = MACROPHYTE_CALC(JW,M).AND.MPRWBC(JW,M).EQ.' ON'
            MACROPHYTE_CALC(JW, M) = CONSTITUENTS .and. MACWBC(JW, M) == "      ON" ! cb 8/24/15
            PRINT_MACROPHYTE(JW, M) = MACROPHYTE_CALC(JW, M) .and. MPRWBC(JW, M) == "      ON" ! cb 8/24/15
            if (MACROPHYTE_CALC(JW, M)) then
                MACROPHYTE_ON = .true.
            end if
        end do
    end do


! Initialize variables for enhanced pH buffering ! entire section ! SR 01/01/12
    PHBUFF_EXIST = .false.
    inquire(FILE="pH_buffering.npt", EXIST=PHBUFF_EXIST)
    if (CONSTITUENTS .and. PHBUFF_EXIST) then
        open(NUNIT, FILE="ph_buffering.npt", STATUS="OLD")
        read(NUNIT, "(///8X,2A8)") PHBUFC, NCALKC
        read(NUNIT, "(//8X,3A8)") NH4BUFC, PO4BUFC, OMBUFC
        read(NUNIT, "(//8X,A8,I8,A8)") OMTYPE, NAGI, POMBUFC
        allocate(SDENI(NAGI), PKI(NAGI), PKSD(NAGI))
        read(NUNIT, "(//(:8X,9F8.0))") (SDENI(J), J = 1, NAGI)
        read(NUNIT, "(//(:8X,9F8.0))") (PKI(J), J = 1, NAGI)
        read(NUNIT, "(//(:8X,9F8.0))") (PKSD(J), J = 1, NAGI)
        close(NUNIT)
        pH_BUFFERING = PHBUFC == "      ON"
        NONCON_ALKALINITY = NCALKC == "      ON"
        AMMONIA_BUFFERING = NH4BUFC == "      ON"
        PHOSPHATE_BUFFERING = PO4BUFC == "      ON"
        OM_BUFFERING = OMBUFC == "      ON"
        POM_BUFFERING = POMBUFC == "      ON" .and. OM_BUFFERING
        if (OM_BUFFERING) then
            SDENI = ABS(SDENI)
            if (OMTYPE == "    DIST") then
                if (ANY(PKSD <= 0)) then
                    WARNING_OPEN = .true.
                    write(WRN, "(A)") "WARNING -- PKSD inputs in the ph_buffering.npt file must be greater than zero."
                    write(WRN, "(A/)") "Please fix your inputs. For now, PKSD values of zero will be set to 1."
                end if
                do JA = 1, NAGI
                    if (PKSD(JA) <= 0) then
                        PKSD(JA) = 1.0
                    end if
                end do
                NAG = 27
                allocate(SDEN(NAG), PK(NAG), FRACT(NAG))
                SDEN = 0.0
                do J = 1, NAG
                    PK(J) = 0.5*J
                end do
                do JA = 1, NAGI
                    SUM = 0.0
                    do J = 1, NAG
                        FRACT(J) = EXP((-0.5)*((PK(J) - PKI(JA))/PKSD(JA))**2)
                        SUM = SUM + FRACT(J)
                    end do
                    do J = 1, NAG
                        SDEN(J) = SDEN(J) + SDENI(JA)*FRACT(J)/SUM
                    end do
                end do
            else
                allocate(SDEN(NAGI), PK(NAGI))
                NAG = NAGI
                SDEN = SDENI
                PK = PKI
                OMTYPE = " MONO"
            end if
        end if
        open(NUNIT, FILE="ph_buffering.opt", STATUS="UNKNOWN")
        write(NUNIT, "(A/)") "Enhanced pH buffering output file"
        write(NUNIT, "(2A)") "Ammonia buffering: ", ADJUSTL(TRIM(NH4BUFC))
        write(NUNIT, "(2A)") "Phosphate buffering: ", ADJUSTL(TRIM(PO4BUFC))
        write(NUNIT, "(2A)") "OM buffering: ", ADJUSTL(TRIM(OMBUFC))
        if (OM_BUFFERING) then
            write(NUNIT, "(2A)") "POM buffering: ", ADJUSTL(TRIM(POMBUFC))
            write(NUNIT, "(2A)") "OM buffer type: ", ADJUSTL(TRIM(OMTYPE))
            write(NUNIT, "(/A/A)") "Inputs:", "Group Site density pKa std.dev."
            do JA = 1, NAGI
                if (OMTYPE == " DIST") then
                    write(NUNIT, "(1X,I3,5X,F8.4,4X,F6.3,3X,F6.3)") JA, SDENI(JA), PKI(JA), PKSD(JA)
                else
                    write(NUNIT, "(1X,I3,5X,F8.4,4X,F6.3,3X,A)") JA, SDENI(JA), PKI(JA), " N/A"
                end if
            end do
            write(NUNIT, "(/A/A)") "Modeled:", "Group Site density pKa"
            do JA = 1, NAG
                write(NUNIT, "(1X,I3,5X,F8.4,4X,F6.3)") JA, SDEN(JA), PK(JA)
            end do
        end if
        close(NUNIT)
    else
        AMMONIA_BUFFERING = .false.
        PHOSPHATE_BUFFERING = .false.
        OM_BUFFERING = .false.
        POM_BUFFERING = .false.
    end if
!NGCTDG = 0
!NGN2=0
!NDC=NDC-1     ! SW 10/20/15
!DO JG=1,NGC   ! SW 10/16/2015
!     IF(CGCS(JG)==-1.0)THEN
!         NGCTDG=JG    ! JG COUNTER
!         NDC=NDC+1    ! THERE CAN BE ONLY 1 TDG AMONG THE CONSITUENTS
!         NGN2=NGCS+JG-1        ! GLOBAL COUNTER
!         EXIT
!     ENDIF
! ENDDO
! IF (NGCTDG == 0) CDWBC(NDC+1,:)= '     OFF'     ! SR 7/15/17
    return
end subroutine INPUT
