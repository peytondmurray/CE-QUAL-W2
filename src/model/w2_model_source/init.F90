subroutine INIT()

    use MAIN
    use GLOBAL;     use NAMESC;     use GEOMC;     use LOGICC;     use PREC;     use SURFHE;     use KINETIC;     use SHADEC;     use EDDY
    use STRUCTURES;     use TRANS;     use TVDC;     use SELWC;     use GDAYC;     use SCREENC;     use TDGAS;     use RSTART
    use MACROPHYTEC;     use POROSITYC;     use ZOOPLANKTONC
    use INITIALVELOCITY;     use BIOENERGETICS;     use CEMAVARS;     use ALGAE_TOXINS
    use CEMASedimentDiagenesis, only: InitCond_SedFlux
    implicit none
    external :: RESTART_OUTPUT
    integer :: NB

!***********************************************************************************************************************************
!**                                             Task 1.1: Variable Initialization                                                 **
!***********************************************************************************************************************************

!***********************************************************************************************************************************
!**                                                 Task 1.1.1: Zero Variables                                                    **
!***********************************************************************************************************************************
    IceQSS = 0.0d00;     WATER_AGE_ACTIVE = .false. ! SR 7/27/2017
    ATM_DEP_LOADING = 0.0;     IN_TOXIN = 0.0
    KB = 0;     KBR = 0;     NAC = 0;     NTAC = 0;     NACD = 0;     NACIN = 0;     NACTR = 0;     NACDT = 0;     NACPR = 0
    NDSP = 0;     HMAX = 0;     KBMAX = 0;     DLXMAX = 0;     KBQIN = 0;     KTQIN = 0;     QGT = 0.0;     QSP = 0.0 ! SW 8/26/15 Initialize Qgt and Qsp for screen output on restart
    NAF = 0;     TISS = 0.0;     CSHE = 0.0;     CIN = 0.0;     TIN = 0.0;     EV = 0.0;     NACATD = 0
    DZ = 0.0D0;     ET = 0.0;     CSHE = 0.0;     A = 0.0D0;     F = 0.0D0;     D = 0.0D0;     C = 0.0D0;     ELTMF = 0.0
    EL = 0.0;     DX = 0.0D0;     ST = 0.0D0;     SB = 0.0D0;     DZQ = 0.0D0;     TSS = 0.0
    HSEG = 0.0;     QSS = 0.0D0;     HPG = 0.0D0;     HDG = 0.0D0;     VSH = 0.0D0;     QDH1 = 0.0D0;     ADMX = 0.0D0;     DECAY = 0.0D0
    ADMZ = 0.0D0;     UYBR = 0.0D0;     GRAV = 0.0D0;     FETCH = 0.0D0;     FETCHU = 0.0D0;     FETCHD = 0.0D0;     DLTTVD = 0.0;     ICETHU = 0.0;     ICETH1 = 0.0
    ICETH2 = 0.0;     P = 0.0D0;     CELRTY = 0.0;     TAU1 = 0.0D0;     TAU2 = 0.0D0;     VOLSR = 0.0;     VOLTR = 0.0;     AF = 0.0;     EF = 0.0
    ELTMS = 0.0;     DM = 0.0D0;     QIN = 0.0D0;     REAER = 0.0;     ST = 0.0D0;     SB = 0.0D0;     ADMX = 0.0D0
    ADMZ = 0.0D0;     HPG = 0.0D0;     HDG = 0.0D0;     RHO = 0.0D0;     JDAYTS = 0.0;     JDAY1 = 0.0;     DEPTHB = 0.0D0;     DEPTHM = 0.0D0;     UXBR = 0.0D0
    BHRHO = 0.0D0;     DLMR = 0.0D0;     SRON = 0.0;     CSSB = 0.0;     Q = 0.0D0;     BH1 = 0.0D0;     BH2 = 0.0D0;     BHR1 = 0.0D0;     BHR2 = 0.0D0
    AVHR = 0.0D0;     GRAV = 0.0D0;     KBP = 0;     DZT = 0.0D0;     AZT = 0.0D0;     KFJW = 0.0;     QC = 0.0D0;     YSS = 0.0;     YSTS = 0.0
    QWD = 0.0D0;     QDTR = 0.0D0;     TTR = 0.0;     CTR = 0.0;     TDTR = 0.0;     QOLDS = 0.0D0;     DTPS = 0.0;     VSTS = 0.0;     VSS = 0.0
    EGT2 = 0.0;     HAB = 100.0;     sedpinflux = 0.0;     sedninflux = 0.0;     FPSS = 0.0;     FPFE = 0.0 ! SR 3/2019
    RS = 0.0;     RN = 0.0;     RB = 0.0;     RE = 0.0;     RC = 0.0;     RANLW = 0.0;     TICAP = 0.0;     TICZR = 0.0;     TICEP = 0.0;     TICMC = 0.0;     VOL = 0.0
    sdfirstadd = .true. ! cb 9/3/17
    BR_NOTECPLOT = .true. ! SW 8/27/2019
    if (.not. RESTART_IN) then
        BR_INACTIVE = .false.;         WARNING_OPEN = .false.;         JDMIN = 0;         EPC = 0.0
        NSPRF = 0;         IZMIN = 0;         KTWB = 2;         KMIN = 1;         IMIN = 1;         NH3GASLOSS = 0.0
        T1 = 0.0D0;         T2 = 0.0D0;         C1 = 0.0D0;         C2 = 0.0D0;         CD = 0.0;         CIN = 0.0;         C1S = 0.0;         KF = 0.0;         CMBRT = 0.0
        KFS = 0.0;         U = 0.0D0;         W = 0.0D0;         SU = 0.0D0;         SW = 0.0D0;         SAZ = 0.0D0;         AZ = 0.0D0;         ESBR = 0.0;         EPD = 0.0
        ETBR = 0.0;         EBRI = 0.0;         DLTLIM = 0.0;         VOLEV = 0.0;         VOLPR = 0.0;         VOLDT = 0.0;         VOLWD = 0.0;         CURRENT = 0.0;         VOLICE = 0.0;         ICEBANK = 0.0
        VOLUH = 0.0;         VOLDH = 0.0;         VOLIN = 0.0;         VOLOUT = 0.0;         VOLSBR = 0.0;         VOLTRB = 0.0;         TSSS = 0.0;         TSSB = 0.0;         EF = 0.0
        TSSEV = 0.0;         TSSPR = 0.0;         TSSTR = 0.0;         TSSDT = 0.0;         TSSWD = 0.0;         TSSUH = 0.0;         TSSDH = 0.0;         TSSIN = 0.0;         CSSK = 0.0
        TSSOUT = 0.0;         TSSICE = 0.0;         TSSUH1 = 0.0;         TSSUH2 = 0.0;         CSSUH1 = 0.0;         CSSUH2 = 0.0;         TSSDH1 = 0.0;         TSSDH2 = 0.0
        CSSDH1 = 0.0;         CSSDH2 = 0.0;         QIND = 0.0;         TIND = 0.0;         CIND = 0.0;         SAVH2 = 0.0;         SAVHR = 0.0;         VOLUH2 = 0.0
        AVH1 = 0.0;         AVH2 = 0.0;         VOLDH2 = 0.0;         Z = 0.0D0;         QUH1 = 0.0D0;         SED = 0.0;         SEDC = 0.0;         SEDN = 0.0
        VS = 0.0;         YS = 0.0;         YST = 0.0;         VST = 0.0;         DTP = 0.0;         QOLD = 0.0;         QSUM = 0.0;         VOLTBR = 0.0;         DLVOL = 0.0;         EVBR = 0.0 ! SW 7/24/2017
        TPOUT = 0.0;         TPTRIB = 0.0;         TPDTRIB = 0.0;         TPWD = 0.0;         TPPR = 0.0;         TPIN = 0.0;         TNOUT = 0.0;         TNTRIB = 0.0;         TNDTRIB = 0.0;         TNWD = 0.0;         TNPR = 0.0;         TNIN = 0.0;         TN_SEDSOD_NH4 = 0.0;         TP_SEDSOD_PO4 = 0.0 ! SW 2/19/16  TP_SEDBURIAL=0.0;TN_SEDBURIAL=0.0;
        ATMDEP_P = 0.0;         ATMDEP_N = 0.0
        SEDP = 0.0;         ICETH = 0.0;         PFLUXIN = 0.0;         NFLUXIN = 0.0 ! SW 4/19/10
        ZMIN = -1000.0
        TKE = 0.0 ! SG 10/4/07
        SEDP = 0.0;         SEDC = 0.0;         SEDN = 0.0;         DLTAV = 0.0;         ELTMJD = 0.0
        MACMBRT = 0.0;         MACRC = 0.0;         SMACRC = 0.0;         MAC = 0.0;         SMAC = 0.0;         MACRM = 0.0;         EPM = 0.0;         macss = 0.0 ! cb 3/8/16
        KTICOL = .false.
    end if
    ANLIM = 1.0;     APLIM = 1;     ASLIM = 1.0;     ALLIM = 1.0;     ENLIM = 1.0;     EPLIM = 1;     ESLIM = 1.0;     ELLIM = 1.0;     KLOC = 1;     ILOC = 1
    MNLIM = 1.0;     MPLIM = 1;     MCLIM = 1.0;     MLLIM = 1.0
    ICESW = 1.0
    HMIN = 1.0E10
    DLXMIN = 1.0E10
    LFPR = BLANK
    CONV = BLANK
    CONV1 = BLANK1
    CNAME2 = ADJUSTR(CNAME2)
    CDNAME2 = ADJUSTR(CDNAME2)
    KFNAME2 = ADJUSTR(KFNAME2)
    TITLE(11) = " "
    TEXT = " "
    ICPL = 0
    if (.not. CONSTITUENTS) then
        if (NBOD > 0) then
            deallocate(NBODC, NBODN, NBODP)
        end if
        NAL = 0;         NEP = 0;         NSS = 0;         NBOD = 0; 
    end if
    do JW = 1, NWB
        GAMMA(:, US(BS(JW)):DS(BE(JW))) = EXH2O(JW)
    end do
!***********************************************************************************************************************************
!**                                            Task 1.1.2: Miscellaneous Variables                                                **
!***********************************************************************************************************************************

! Logical controls
    NEW_PAGE = .true.;     VOLUME_WARNING = .true.;     INITIALIZE_GRAPH = .true.;     UPDATE_GRAPH = .true.
    ICE = .false.;     FLUX = .false.;     PUMPON = .false.
    TDG_GATE = .false.;     TDG_SPILLWAY = .false.;     INTERNAL_WEIR = .false.;     SURFACE_WARNING = .false.
    PRINT_CONST = .false.;     PRINT_DERIVED = .false.;     ERROR_OPEN = .false.
    LIMITING_FACTOR = .false.
    HEAD_BOUNDARY = .false.;     PRINT_HYDRO = .false.;     ONE_LAYER = .false.;     ZERO_SLOPE = .true.
    INTERNAL_FLOW = .false.;     DAM_INFLOW = .false.;     DAM_OUTFLOW = .false.;     HEAD_FLOW = .false. !TC 08/03/04
    UPDATE_RATES = .false. !TC 08/03/04
    WEIR_CALC = NIW > 0;     GATES = NGT > 0;     PIPES = NPI > 0
    PUMPS = NPU > 0;     SPILLWAY = NSP > 0;     TRIBUTARIES = NTR > 0
    WITHDRAWALS = NWD > 0
    VOLUME_BALANCE = VBC == "      ON"
    PLACE_QIN = PQC == "      ON";     EVAPORATION = EVC == "      ON"
    ENERGY_BALANCE = EBC == "      ON";     RH_EVAP = RHEVC == "      ON"
    PRECIPITATION = PRC == "      ON";     RESTART_OUT = RSOC == "      ON"
    INTERP_TRIBS = TRIC == "      ON";     INTERP_DTRIBS = DTRIC == "      ON"
    INTERP_HEAD = HDIC == "      ON";     INTERP_INFLOW = QINIC == "      ON"
    INTERP_OUTFLOW = STRIC == "      ON";     INTERP_WITHDRAWAL = WDIC == "      ON"
    INTERP_GATE = GTIC == "      ON" ! cb 8/13/2010
!INTERP_METEOROLOGY    = METIC       == '      ON'; DOWNSTREAM_OUTFLOW = WDOC   == '      ON'
    INTERP_METEOROLOGY = METIC == "      ON"
    if (WDOC == "      ON" .or. WDOC == "     ONH" .or. WDOC == "     ONS") then
        DOWNSTREAM_OUTFLOW = .true.
    end if ! cb 4/11/18
    CELERITY_LIMIT = CELC == "      ON";     VISCOSITY_LIMIT = VISC == "      ON"
    PRINT_HYDRO = HPRWBC == "      ON" ! HYDRO_PLOT            = HPLTC       == '      ON'; 
    LIMITING_DLT = HPRWBC(1, :) == "      ON";     FETCH_CALC = FETCHC == "      ON"
    SCREEN_OUTPUT = SCRC == "      ON";     SNAPSHOT = SNPC == "      ON"
    CONTOUR = CPLC == "      ON";     VECTOR = VPLC == "      ON"
    PROFILE = PRFC == "      ON" !; SPREADSHEET        = SPRC   == '      ON'; 
    SPREADSHEET = .false. ! INITIALIZE SW 2/10/2019
    ATM_DEPOSITION = ATM_DEPOSITIONC == "      ON"

    GAS_TRANSFER_UPDATE = .false.
    if (CAC(NDO) == "      ON" .or. CAC(NCH4) == "      ON" .or. CAC(NH2S) == "      ON" .or. CAC(NN2) == "      ON" .or. CAC(NTIC) == "      ON" .or. CAC(NDGP) == "      ON") then
        GAS_TRANSFER_UPDATE = .true.
    else
        do JG = 1, NGC
            if (CGKLF(JG) > 0.0) then
                GAS_TRANSFER_UPDATE = .true.
                exit
            end if
        end do
    end if




    do JW = 1, NWB
        if (SPRC(JW) == "      ON") then ! SW 9/28/2018
            SPREADSHEET(JW) = .true.
        else
            if (SPRC(JW) == "     ONV") then
                SPREADSHEET(JW) = .true.
            end if
        end if
    end do

    TIME_SERIES = TSRC == "      ON";     READ_RADIATION = SROC == "      ON"
    ICE_CALC = ICEC == "      ON" .or. ICEC == "    ONWB"
    INTERP_EXTINCTION = EXIC == "      ON";     READ_EXTINCTION = EXC == "      ON"
    NO_INFLOW = QINC == "     OFF";     NO_OUTFLOW = QOUTC == "     OFF"
    NO_HEAT = HEATC == "     OFF";     NO_WIND = WINDC == "     OFF"
    SPECIFY_QTR = TRC == " SPECIFY";     DIST_TRIBS = DTRC == "      ON"
    IMPLICIT_VISC = AZSLC == "     IMP";     UPWIND = SLTRC == "  UPWIND"
    ULTIMATE = SLTRC == "ULTIMATE";     TERM_BY_TERM = SLHTC == "    TERM"
    MANNINGS_N = FRICC == "    MANN";     PLACE_QTR = TRC == " DENSITY"
    LATERAL_SPILLWAY = LATSPC /= "    DOWN";     LATERAL_PUMP = LATPUC /= "    DOWN"
    LATERAL_GATE = LATGTC /= "    DOWN";     LATERAL_PIPE = LATPIC /= "    DOWN"
    TRAPEZOIDAL = GRIDC == "    TRAP" !SW 07/16/04
    EPIPHYTON_CALC = CONSTITUENTS .and. EPIC == "      ON"
    MASS_BALANCE = CONSTITUENTS .and. MBC == "      ON"
    SUSP_SOLIDS = CONSTITUENTS .and. CAC(NSSS) == "      ON"
    OXYGEN_DEMAND = CONSTITUENTS .and. CAC(NDO) == "      ON"
    WATER_AGE_ACTIVE = CONSTITUENTS .and. CAC(NWAGE) == "      ON"
    SEDIMENT_CALC = CONSTITUENTS .and. SEDCc == "      ON"
    zooplankton_CALC = CONSTITUENTS .and. cac(nzooS) == "      ON"
    SEDIMENT_RESUSPENSION = CONSTITUENTS .and. SEDRC == "      ON"
!  DERIVED_PLOT          = CONSTITUENTS .AND. CDPLTC      == '      ON'
    DERIVED_CALC = CONSTITUENTS .and. ANY(CDWBC == "      ON")
    PH_CALC = CONSTITUENTS .and. CDWBC(PH_DER, :) == "      ON"
    if (.not. PH_CALC(1)) then
        PH_CALC = CAC(NTIC) == "      ON"
    end if ! CALL THIS ROUTINE EVEN IF Ph IS OFF FOR CO2 GAS CALCULATION
    PRINT_EPIPHYTON = CONSTITUENTS .and. EPIPRC == "      ON" .and. EPIPHYTON_CALC
    PRINT_SEDIMENT = CONSTITUENTS .and. SEDPRC == "      ON" .and. SEDIMENT_CALC
    SEDIMENT_CALC1 = CONSTITUENTS .and. SEDCc1 == "      ON"
    SEDIMENT_CALC2 = CONSTITUENTS .and. SEDCc2 == "      ON"
    PRINT_SEDIMENT1 = CONSTITUENTS .and. SEDPRC1 == "      ON" .and. SEDIMENT_CALC1
    PRINT_SEDIMENT2 = CONSTITUENTS .and. SEDPRC2 == "      ON" .and. SEDIMENT_CALC2
    FRESH_WATER = CONSTITUENTS .and. WTYPEC == "   FRESH" .and. CAC(NTDS) == "      ON"
    SALT_WATER = CONSTITUENTS .and. WTYPEC == "    SALT" .and. CAC(NTDS) == "      ON"
!  CONSTITUENT_PLOT      = CONSTITUENTS .AND. CPLTC       == '      ON' .AND. CAC       == '      ON'
    DETAILED_ICE = ICE_CALC .and. SLICEC == "  DETAIL"
    LEAP_YEAR = MOD(YEAR, 4) == 0
    ICE_COMPUTATION = ANY(ICE_CALC)
    END_RUN = JDAY > TMEND
    UPDATE_KINETICS = CONSTITUENTS
    where (READ_EXTINCTION)
        EXOM = 0.0
        EXSS = 0.0
    end where
    if (CONSTITUENTS) then
        SUSP_SOLIDS = .false.
        FLUX = FLXC == "      ON"
        PRINT_CONST = CPRWBC == "      ON"
        PRINT_DERIVED = CDWBC == "      ON"
        if (ANY(CAC(NSSS:NSSE) == "      ON")) then
            SUSP_SOLIDS = .true.
        end if
        if (ANY(CAC(NSSS:NCT) == "      ON")) then
            UPDATE_RATES = .true.
        end if
        do JA = 1, NAL
            LIMITING_FACTOR(JA) = CONSTITUENTS .and. CAC(NAS - 1 + JA) == "      ON" .and. LIMC == "      ON"
            ALG_CALC(JA) = CAC(NAS - 1 + JA) == "      ON"
        end do
        do NB = 1, NBOD
            BOD_CALC(NB) = CAC(NBODS - 1 + NB) == "      ON"
            BOD_CALCP(NB) = CAC(NBODS - 1 + NBOD + NB) == "      ON" ! cb 5/19/2011
            BOD_CALCN(NB) = CAC(NBODS - 1 + 2*NBOD + NB) == "      ON" ! cb 5/19/2011
        end do
        DSI_CALC = CAC(NDSI) == "      ON" ! cb 10/12/11
        PO4_CALC = CAC(NPO4) == "      ON" ! cb 10/12/11
        N_CALC = CAC(NNH4) == "      ON" .or. CAC(NNO3) == "      ON" ! cb 10/12/11

    end if
    JBDAM = 0
    CDHS = DHS
    do JB = 1, NBR
        UP_FLOW(JB) = UHS(JB) == 0
        DN_FLOW(JB) = DHS(JB) == 0
        UP_HEAD(JB) = UHS(JB) /= 0
        UH_INTERNAL(JB) = UHS(JB) > 0
        if (UP_HEAD(JB)) then
            do JJB = 1, NBR
                if (ABS(UHS(JB)) >= US(JJB) .and. ABS(UHS(JB)) <= DS(JJB)) then
                    if (ABS(UHS(JB)) == DS(JJB)) then
                        if (DHS(JJB) == US(JB)) then
                            UP_FLOW(JB) = .true.
                            HEAD_FLOW(JB) = .true.
                            INTERNAL_FLOW(JB) = .true.
                            UP_HEAD(JB) = .false.
                            UH_INTERNAL(JB) = .false.
                        end if
                        if (UHS(JB) < 0) then
                            do JJJB = 1, NBR
                                if (ABS(UHS(JB)) == DS(JJJB)) then
                                    exit
                                end if ! CB 1/2/05
                            end do
                            UP_FLOW(JB) = .true.
                            DAM_INFLOW(JB) = .true. !TC 08/03/04
                            DAM_OUTFLOW(JJJB) = .true. !TC 08/03/04
                            INTERNAL_FLOW(JB) = .true.
                            UP_HEAD(JB) = .false.
                            UHS(JB) = ABS(UHS(JB))
                            JBDAM(JJJB) = JB
                        end if
                    end if
                    exit
                end if
            end do
        end if
        DH_INTERNAL(JB) = DHS(JB) > 0;         DN_HEAD(JB) = DHS(JB) /= 0;         UH_EXTERNAL(JB) = UHS(JB) == -1
        DH_EXTERNAL(JB) = DHS(JB) == -1;         UQ_EXTERNAL(JB) = UHS(JB) == 0;         DQ_EXTERNAL(JB) = DHS(JB) == 0
        DQ_INTERNAL(JB) = DQB(JB) > 0;         UQ_INTERNAL(JB) = UQB(JB) > 0 .and. .not. DAM_INFLOW(JB) !TC 08/03/04
    end do
    do JW = 1, NWB
        if (TKELATPRDCONST(JW) > 0.0) then
            TKELATPRD(JW) = .true.
        end if
        if (STRICK(JW) > 0.0) then
            STRICKON(JW) = .true.
        end if
        do JB = BS(JW), BE(JW)
            if (UH_EXTERNAL(JB) .or. DH_EXTERNAL(JB)) then
                HEAD_BOUNDARY(JW) = .true.
            end if
            if (SLOPE(JB) /= 0.0) then
                ZERO_SLOPE(JW) = .false.
            end if
        end do
    end do
    where (CAC == "     OFF")
        CPLTC = "     OFF"
    end where
! systdg - Add WBSEG
    WBSEG(1:IMX) = 1
    do JW = 1, NWB
        do JB = BS(JW), BE(JW)
            do I = 1, IMX
                if (I >= US(JB) .and. I <= DS(JB)) then
                    WBSEG(I) = JW
                end if
            end do
        end do
    end do
! systdg - Add WBSEG END
! Kinetic flux variables

    KFNAME(1) = "TISS settling in - source, kg/day            ";     KFNAME(2) = "TISS settling out - sink, kg/day             "
    KFNAME(3) = "PO4 algal respiration - source, kg/day       ";     KFNAME(4) = "PO4 algal growth - sink, kg/day              "
    KFNAME(5) = "PO4 algal net- source/sink, kg/day           ";     KFNAME(6) = "PO4 epiphyton respiration - source, kg/day   "
    KFNAME(7) = "PO4 epiphyton growth - sink, kg/day          ";     KFNAME(8) = "PO4 epiphyton net- source/sink, kg/day       "
    KFNAME(9) = "PO4 POM decay - source, kg/day               ";     KFNAME(10) = "PO4 DOM decay - source, kg/day               "
    KFNAME(11) = "PO4 OM decay - source, kg/day                ";     KFNAME(KF_PO4_SD) = "PO4 sediment decay - source, kg/day          "
    KFNAME(KF_PO4_SR) = "PO4 SOD release - source, kg/day             ";     KFNAME(14) = "PO4 net settling  - source/sink, kg/day      "
    KFNAME(15) = "NH4 nitrification - sink, kg/day             ";     KFNAME(16) = "NH4 algal respiration - source, kg/day       "
    KFNAME(17) = "NH4 algal growth - sink, kg/day              ";     KFNAME(18) = "NH4 algal net - source/sink, kg/day          "
    KFNAME(19) = "NH4 epiphyton respiration - source, kg/day   ";     KFNAME(20) = "NH4 epiphyton growth - sink, kg/day          "
    KFNAME(21) = "NH4 epiphyton net - source/sink, kg N/day    ";     KFNAME(22) = "NH4 POM decay - source, kg N/day             "
    KFNAME(23) = "NH4 DOM decay  - source, kg N/day            ";     KFNAME(24) = "NH4 OM decay - source, kg N/day              "
    KFNAME(KF_NH4_SD) = "NH4 sediment decay - source, kg N/day        ";     KFNAME(KF_NH4_SR) = "NH4 SOD release - source, kg N/day           "
    KFNAME(KF_NH3GAS) = "NH3 gas loss - sink, kg N/day                "

    KFNAME(KF_NO3D) = "NO3 denitrification - sink, kg/day           ";     KFNAME(KF_NO3AG) = "NO3 algal growth - sink, kg/day              "
    KFNAME(KF_NO3EG) = "NO3 epiphyton growth - sink, kg/day          ";     KFNAME(KF_NO3SED) = "NO3 sediment uptake - sink, kg/day           "
    KFNAME(32) = "DSi algal growth - sink, kg/day              ";     KFNAME(33) = "DSi epiphyton growth - sink, kg/day          "
    KFNAME(34) = "DSi PBSi decay - source, kg/day              ";     KFNAME(35) = "DSi sediment decay - source, kg/day          "
    KFNAME(36) = "DSi SOD release  - source, kg/day            ";     KFNAME(37) = "DSi net settling - source/sink, kg/day       "
    KFNAME(38) = "PBSi algal mortality  - source, kg/day       ";     KFNAME(39) = "PBSi net settling - source/sink, kg/day      "
    KFNAME(40) = "PBSi decay - sink, kg/day                    "

    KFNAME(41) = "LDOM decay - sink, kg/day                    "
    KFNAME(42) = "LDOM decay to RDOM - sink, kg/day            ";     KFNAME(43) = "RDOM decay - sink, kg/day                    "
    KFNAME(44) = "LDOM algal mortality - source, kg/day        ";     KFNAME(45) = "LDOM epiphyton mortality - source, kg/day    "
    KFNAME(46) = "LPOM decay - sink, kg/day                    ";     KFNAME(47) = "LPOM decay to RPOM - sink, kg/day            "
    KFNAME(48) = "RPOM decay - sink, kg/day                    ";     KFNAME(49) = "LPOM algal production - source, kg/day       "
    KFNAME(50) = "LPOM epiphyton production - source, kg/day   ";     KFNAME(51) = "LPOM net settling - source/sink, kg/day      "
    KFNAME(52) = "RPOM net settling - source/sink, kg/day      ";     KFNAME(53) = "CBOD decay - sink, kg/day                    "
    KFNAME(54) = "DO algal production  - source, kg/day        ";     KFNAME(56) = "DO algal respiration - sink, kg/day          " ! cb 6/2/2009
    KFNAME(55) = "DO epiphyton production  - source, kg/day    ";     KFNAME(57) = "DO epiphyton respiration - sink, kg/day      " ! cb 6/2/2009
    KFNAME(58) = "DO POM decay - sink, kg/day                  ";     KFNAME(59) = "DO DOM decay - sink, kg/day                  "
    KFNAME(60) = "DO OM decay - sink, kg/day                   ";     KFNAME(61) = "DO nitrification - sink, kg/day              "
    KFNAME(62) = "DO CBOD uptake - sink, kg/day                ";     KFNAME(63) = "DO reaeration - source/sink, kg/day          "
    KFNAME(KF_DO_SED) = "DO sediment uptake - sink, kg/day            ";     KFNAME(KF_DO_SOD) = "DO SOD uptake - sink, kg/day                 "
    KFNAME(66) = "TIC algal uptake - sink, kg/day              ";     KFNAME(67) = "TIC epiphyton uptake - sink, kg/day          "
    KFNAME(68) = "Sediment decay - sink, kg/day                ";     KFNAME(69) = "Sediment algal settling - sink, kg/day       "
    KFNAME(70) = "Sediment LPOM settling - source,kg/day       ";     KFNAME(71) = "Sediment net settling - source/sink, kg/day  "
    KFNAME(72) = "SOD decay - sink, kg/day                     "

    KFNAME(73) = "LDOM P algal mortality - source, kg/day      ";     KFNAME(74) = "LDOM P epiphyton mortality - source, kg/day  "
    KFNAME(75) = "LPOM P algal production- source, kg/day      ";     KFNAME(76) = "LPOM P net settling - source/sink, kg/day    "
    KFNAME(77) = "RPOM P net settling - source/sink, kg/day    "
    KFNAME(78) = "LDOM P algal mortality - source, kg/day      ";     KFNAME(79) = "LDOM P epiphyton mortality - source, kg/day  "
    KFNAME(80) = "LPOM P algal production- source, kg/day      ";     KFNAME(81) = "LPOM P net settling - source/sink, kg/day    "
    KFNAME(82) = "RPOM P net settling - source/sink, kg/day    "
    KFNAME(83) = "Sediment P decay - sink, kg/day              ";     KFNAME(84) = "Sediment algal P settling - source, kg/day   "
    KFNAME(85) = "Sediment P LPOM settling - source,kg/day     ";     KFNAME(86) = "Sediment net P settling - source/sink, kg/day"
    KFNAME(87) = "Sediment epiphyton P settling - source,kg/day"
    KFNAME(88) = "Sediment N decay - sink, kg/day              ";     KFNAME(89) = "Sediment algal N settling - source, kg/day   "
    KFNAME(90) = "Sediment N LPOM settling - source,kg/day     ";     KFNAME(91) = "Sediment net N settling - source/sink, kg/day"
    KFNAME(92) = "Sediment epiphyton N settling - source,kg/day"
    KFNAME(93) = "Sediment C decay - sink, kg/day              ";     KFNAME(94) = "Sediment algal C settling - source, kg/day   "
    KFNAME(95) = "Sediment C LPOM settling - source,kg/day     ";     KFNAME(96) = "Sediment net C settling - source/sink, kg/day"
    KFNAME(97) = "Sediment epiphyton C settling - source,kg/day"
    KFNAME(98) = "Sediment N denitrification - source, kg/day  "
    KFNAME(99) = "PO4 macrophyte resp - source, kg/day         "
    KFNAME(100) = "PO4 macrophyte growth - sink, kg/day         "
    KFNAME(101) = "NH4 macrophyte resp - source, kg/day         "
    KFNAME(102) = "NH4 macrophyte growth - sink, kg/day         "
    KFNAME(103) = "LDOM macrophyte mort  - source, kg/day       "
    KFNAME(104) = "LPOM macrophyte mort  - source, kg/day       "
    KFNAME(105) = "RPOM macrophyte mort  - source, kg/day       "
    KFNAME(106) = "DO  macrophyte production  - source, kg/day  "
    KFNAME(107) = "DO  macrophyte respiration - sink, kg/day    "
    KFNAME(108) = "TIC macrophyte growth/resp  - S/S, kg/day    "
    KFNAME(109) = "CBOD settling - sink, kg/day                 "
    KFNAME(110) = "Sediment CBOD settling - source, kg/day      "
    KFNAME(111) = "Sediment CBOD P settling - source, kg/day    "
    KFNAME(112) = "Sediment CBOD N settling - source, kg/day    "
    KFNAME(113) = "Sediment CBOD C settling - source, kg/day    "
    KFNAME(114) = "Sediment Burial - sink, kg/day               "
    KFNAME(KF_SED_PBURIAL) = "Sediment P Burial - sink, kg/day             "
    KFNAME(KF_SED_NBURIAL) = "Sediment N Burial - sink, kg/day             "
    KFNAME(117) = "Sediment C Burial - sink, kg/day             "
    KFNAME(118) = "CBOD P settling - sink, kg/day               "
    KFNAME(119) = "CBOD N settling - sink, kg/day               "
    KFNAME(KF_CO2X) = "CO2 gas exchange air/water interface, kg/day "
    KFNAME(KF_DOH2S) = "DO H2S decay - sink, kg/day                  "
    KFNAME(122) = "H2S gas exchange air/water interface, kg/day "
    KFNAME(123) = "H2S decay - sink, kg/day                     "
    KFNAME(124) = "H2S release 0 order model   - source, kg/day "
    KFNAME(KF_DOCH4) = "DO CH4 decay - sink, kg/day                  "
    KFNAME(126) = "CH4 gas exchange air/water interface, kg/day "
    KFNAME(127) = "CH4 decay - sink, kg/day                     "
    KFNAME(128) = "CH4 release 0 order model   - source, kg/day "
    KFNAME(KF_FE2D) = "Fe(II) oxidation water column - sink, kg/day "
    KFNAME(130) = "DO Fe(II) oxidation water col.- sink, kg/day "
    KFNAME(131) = "FeOOH settling from water col. - sink, kg/day"
    KFNAME(132) = "Settling of FeOOH into water layer, kg/d     "
    KFNAME(KF_MN2D) = "Mn(II) oxidation water column - sink, kg/day "
    KFNAME(134) = "DO Mn(II) oxidation water col.- sink, kg/day "
    KFNAME(135) = "MnO2 settling from water col. - sink, kg/day "
    KFNAME(136) = "Settling of MnO2  into water layer, kg/d     "
    KFNAME(KF_SDINC) = "C to Sed. Diagenesis module - source, kg/day "
    KFNAME(138) = "N to Sed. Diagenesis module - source, kg/day "
    KFNAME(139) = "P to Sed. Diagenesis module - source, kg/day "
    KFNAME(140) = "DO sediment diagenesis uptake - sink, kg/day "

    KFNAME(KF_SEDD) = "Labile standing biomass decay- sink, kg/day  " ! OPTIONAL VARIABLE FOR STANDING ORGANIC MATTER LIKE TREES IN A WATER COLUMN
    KFNAME(142) = "Refract. stand. biomass decay- sink, kg/day  "

! Convert rates from per-day to per-second

    if (CONSTITUENTS) then
        AE = AE/DAY;         AM = AM/DAY;         AR = AR/DAY;         AG = AG/DAY;         AS = AS/DAY
        EE = EE/DAY;         EM = EM/DAY;         ER = ER/DAY;         EG = EG/DAY;         EB = EB/DAY
        CGS = CGS/DAY;         CG0DK = CG0DK/DAY;         CG1DK = CG1DK/DAY;         SSS = SSS/DAY
        H2S1DK = H2S1DK/DAY;         CH41DK = CH41DK/DAY
        KFE_OXID = KFE_OXID/DAY;         KFE_RED = KFE_RED/DAY;         FeSetVel = FeSetVel/DAY
        KMN_OXID = KMN_OXID/DAY;         KMN_RED = KMN_RED/DAY;         MnSetVel = MnSetVel/DAY
        BACT1DK = BACT1DK/DAY;         BACTLDK = BACTLDK/DAY;         BACTS = BACTS/DAY
        PSIS = PSIS/DAY;         POMS = POMS/DAY;         SDK = SDK/DAY;         NH4DK = NH4DK/DAY;         NO3DK = NO3DK/DAY
        NO3S = NO3S/DAY;         PSIDK = PSIDK/DAY;         LRDDK = LRDDK/DAY;         LRPDK = LRPDK/DAY;         LDOMDK = LDOMDK/DAY
        LPOMDK = LPOMDK/DAY;         RDOMDK = RDOMDK/DAY;         RPOMDK = RPOMDK/DAY;         KBOD = KBOD/DAY;         seds = seds/day !v3.5
!   
        LPOMHK = LPOMHK/DAY;         RPOMHK = RPOMHK/DAY
        LDOMPDK = LDOMPDK/DAY;         RDOMPDK = RDOMPDK/DAY;         LRDOMPDK = LRDOMPDK/DAY
        LPOMPDK = LPOMPDK/DAY;         RPOMPDK = RPOMPDK/DAY;         LRPOMPDK = LRPOMPDK/DAY
        LDOMNDK = LDOMNDK/DAY;         RDOMNDK = RDOMNDK/DAY;         LRDOMNDK = LRDOMNDK/DAY
        LPOMNDK = LPOMNDK/DAY;         RPOMNDK = RPOMNDK/DAY;         LRPOMNDK = LRPOMNDK/DAY
        LDOMCDK = LDOMCDK/DAY;         RDOMCDK = RDOMCDK/DAY;         LRDOMCDK = LRDOMCDK/DAY
        LPOMCDK = LPOMCDK/DAY;         RPOMCDK = RPOMCDK/DAY;         LRPOMCDK = LRPOMCDK/DAY
!
        SDK1 = SDK1/DAY;         SDK2 = SDK2/DAY ! Amaila
        SEDB = SEDB/DAY !CB 11/27/06
        CBODS = CBODS/DAY !CB 7/23/07
!    SSFLOC = SSFLOC/DAY                                                                                                 !SR 04/21/13
        do JW = 1, NWB
            SOD(US(BS(JW)) - 1:DS(BE(JW)) + 1) = SOD(US(BS(JW)) - 1:DS(BE(JW)) + 1)/DAY*FSOD(JW)
        end do
        do J = 1, NEP
            EBR(:, :, J) = EB(J)
        end do

        MG = MG/DAY
        MR = MR/DAY
        MM = MM/DAY
        ZG = ZG/DAY
        ZR = ZR/DAY
        ZM = ZM/DAY


    end if

! Convert slope to angle alpha in radians

    ALPHA = ATAN(SLOPE)
    SINA = SIN(ALPHA)
    SINAC = SIN(ATAN(SLOPEC))
    COSA = COS(ALPHA)

! Time and printout control variables

    if (.not. RESTART_IN) then
        JDAY = TMSTRT
        ELTM = TMSTRT*DAY
        DLT = DLTMAX(1)
        DLTS = DLT
        MINDLT = DLT
        NIT = 0
        NV = 0
        DLTDP = 1;         RSODP = 1;         TSRDP = 1;         SNPDP = 1;         VPLDP = 1;         PRFDP = 1
        SPRDP = 1;         CPLDP = 1;         SCRDP = 1;         FLXDP = 1;         WDODP = 1
        NXTSEDIAG = TMSTRT
        NXWL = TMSTRT;         NXFLOWBAL = TMSTRT;         NXNPBAL = TMSTRT
        do JW = 1, NWB
            do J = 1, NOD
                if (TMSTRT > SNPD(J, JW)) then
                    SNPD(J, JW) = TMSTRT
                end if;                 if (TMSTRT > PRFD(J, JW)) then
                    PRFD(J, JW) = TMSTRT
                end if
                if (TMSTRT > SPRD(J, JW)) then
                    SPRD(J, JW) = TMSTRT
                end if;                 if (TMSTRT > CPLD(J, JW)) then
                    CPLD(J, JW) = TMSTRT
                end if
                if (TMSTRT > VPLD(J, JW)) then
                    VPLD(J, JW) = TMSTRT
                end if;                 if (TMSTRT > SCRD(J, JW)) then
                    SCRD(J, JW) = TMSTRT
                end if
                if (TMSTRT > FLXD(J, JW)) then
                    FLXD(J, JW) = TMSTRT
                end if
            end do
            NXTMSN(JW) = SNPD(SNPDP(JW), JW);             NXTMPR(JW) = PRFD(PRFDP(JW), JW);             NXTMSP(JW) = SPRD(SPRDP(JW), JW)
            NXTMCP(JW) = CPLD(CPLDP(JW), JW);             NXTMVP(JW) = VPLD(VPLDP(JW), JW);             NXTMSC(JW) = SCRD(SCRDP(JW), JW)
            NXTMFL(JW) = FLXD(FLXDP(JW), JW)
        end do
        do J = 1, NOD
            if (TMSTRT > TSRD(J)) then
                TSRD(J) = TMSTRT
            end if;             if (TMSTRT > WDOD(J)) then
                WDOD(J) = TMSTRT
            end if
            if (TMSTRT > RSOD(J)) then
                RSOD(J) = TMSTRT
            end if;             if (TMSTRT > DLTD(J)) then
                DLTD(J) = TMSTRT
            end if
        end do
        NXTMTS = TSRD(TSRDP);         NXTMRS = RSOD(RSODP);         NXTMWD = WDOD(WDODP)

        NXTMWD_SEC = WDOD(WDODP)*86400. ! cb 4/6/18 frequency test seconds 

        if (bioexp) then
            BIODP = 1 ! MLM BIOEXP -BIOENERGETICS
            NXBIO = BIOD(1) !BIOD(BIODP) ! MLM BIOEXP  BIOENERGETICS
            NXTBIO = BIOD(1) ! MLM INITIALIZE THE OUTPUT VARIABLES
            BIOD(NBIO + 1:NOD) = TMEND + 1.0 ! MLM BIOEXP BIOENERGETICS
        end if
    end if
    TSRD(NTSR + 1:NOD) = TMEND + 1.0;     WDOD(NWDO + 1:NOD) = TMEND + 1.0;     RSOD(NRSO + 1:NOD) = TMEND + 1.0;     DLTD(NDLT + 1:NOD) = TMEND + 1.0

    do JW = 1, NWB
        SNPD(NSNP(JW) + 1:NOD, JW) = TMEND + 1.0;         PRFD(NPRF(JW) + 1:NOD, JW) = TMEND + 1.0;         SPRD(NSPR(JW) + 1:NOD, JW) = TMEND + 1.0
        VPLD(NVPL(JW) + 1:NOD, JW) = TMEND + 1.0;         CPLD(NCPL(JW) + 1:NOD, JW) = TMEND + 1.0;         SCRD(NSCR(JW) + 1:NOD, JW) = TMEND + 1.0
        FLXD(NFLX(JW) + 1:NOD, JW) = TMEND + 1.0
    end do
    JDAYG = JDAY
    JDAYNX = JDAYG + 1
    NXTVD = JDAY
    do J = 1, NOD ! SW 12/9/2016
        if (DLTD(J) == DLTD(J + 1)) then
            DLTDP = DLTDP + 1
        else
            exit
        end if
    end do

    DLTMAXX = DLTMAX(DLTDP)
    DLTFF = DLTF(DLTDP)
    CURMAX = DLTMAX(DLTDP)/DLTF(DLTDP)

! Hydraulic structures

    if (SPILLWAY) then
        do JS = 1, NSP
            if (LATERAL_SPILLWAY(JS)) then
                if (IDSP(JS) /= 0) then
                    TRIBUTARIES = .true.
                    WITHDRAWALS = .true.
                else
                    WITHDRAWALS = .true.
                end if
            end if
            do JB = 1, NBR
                if (IUSP(JS) >= US(JB) .and. IUSP(JS) <= DS(JB)) then
                    exit
                end if
            end do
            JBUSP(JS) = JB
            if (IUSP(JS) == DS(JBUSP(JS)) .and. .not. LATERAL_SPILLWAY(JS)) then
                NST = NST + 1
            end if
            do JW = 1, NWB
                if (JB >= BS(JW) .and. JB <= BE(JW)) then
                    exit
                end if
            end do
            JWUSP(JS) = JW
            if (IDSP(JS) > 0) then
                do JB = 1, NBR
                    if (IDSP(JS) >= US(JB) .and. IDSP(JS) <= DS(JB)) then
                        exit
                    end if
                end do
                JBDSP(JS) = JB
                do JW = 1, NWB
                    if (JB >= BS(JW) .and. JB <= BE(JW)) then
                        exit
                    end if
                end do
                JWDSP(JS) = JW
            else
                JBDSP(JS) = 1
                JWDSP(JS) = 1
            end if
        end do
    end if
    if (PIPES) then
        do JP = 1, NPI
            if (LATERAL_PIPE(JP)) then
                if (IDPI(JP) /= 0) then
                    TRIBUTARIES = .true.
                    WITHDRAWALS = .true.
                else
                    WITHDRAWALS = .true.
                end if
            end if
            do JB = 1, NBR
                if (IUPI(JP) >= US(JB) .and. IUPI(JP) <= DS(JB)) then
                    exit
                end if
            end do
            JBUPI(JP) = JB
            if (IUPI(JP) == DS(JBUPI(JP)) .and. .not. LATERAL_PIPE(JP)) then
                NST = NST + 1
            end if
            do JW = 1, NWB
                if (JB >= BS(JW) .and. JB <= BE(JW)) then
                    exit
                end if
            end do
            JWUPI(JP) = JW
            if (IDPI(JP) > 0) then
                do JB = 1, NBR
                    if (IDPI(JP) >= US(JB) .and. IDPI(JP) <= DS(JB)) then
                        exit
                    end if
                end do
                JBDPI(JP) = JB
                do JW = 1, NWB
                    if (JB >= BS(JW) .and. JB <= BE(JW)) then
                        exit
                    end if
                end do
                JWDPI(JP) = JW
            else
                JBDPI(JP) = 1
                JWDPI(JP) = 1
            end if
        end do
    end if
    if (GATES) then
        do JG = 1, NGT
            if (LATERAL_GATE(JG)) then
                if (IDGT(JG) /= 0) then
                    TRIBUTARIES = .true.
                    WITHDRAWALS = .true.
                else
                    WITHDRAWALS = .true.
                end if
            end if
            do JB = 1, NBR
                if (IUGT(JG) >= US(JB) .and. IUGT(JG) <= DS(JB)) then
                    exit
                end if
            end do
            JBUGT(JG) = JB
            if (IUGT(JG) == DS(JBUGT(JG)) .and. .not. LATERAL_GATE(JG)) then
                NST = NST + 1
            end if
            do JW = 1, NWB
                if (JB >= BS(JW) .and. JB <= BE(JW)) then
                    exit
                end if
            end do
            JWUGT(JG) = JW
            if (IDGT(JG) > 0) then
                do JB = 1, NBR
                    if (IDGT(JG) >= US(JB) .and. IDGT(JG) <= DS(JB)) then
                        exit
                    end if
                end do
                JBDGT(JG) = JB
                do JW = 1, NWB
                    if (JB >= BS(JW) .and. JB <= BE(JW)) then
                        exit
                    end if
                end do
                JWDGT(JG) = JW
            else
                JBDGT(JG) = 1 ! SW 3/24/10
                JWDGT(JG) = 1 ! SW 3/24/10
            end if
        end do
    end if
    if (PUMPS) then
        do JP = 1, NPU
            if (LATERAL_PUMP(JP)) then
                if (IDPU(JP) /= 0) then
                    TRIBUTARIES = .true.
                    WITHDRAWALS = .true.
                else
                    WITHDRAWALS = .true.
                end if
            end if
            do JB = 1, NBR
                if (IUPU(JP) >= US(JB) .and. IUPU(JP) <= DS(JB)) then
                    exit
                end if
            end do
            JBUPU(JP) = JB
            if (IUPU(JP) == DS(JBUPU(JP)) .and. .not. LATERAL_PUMP(JP)) then
                NST = NST + 1
            end if
            do JW = 1, NWB
                if (JB >= BS(JW) .and. JB <= BE(JW)) then
                    exit
                end if
            end do
            JWUPU(JP) = JW
            if (IDPU(JP) > 0) then
                do JB = 1, NBR
                    if (IDPU(JP) >= US(JB) .and. IDPU(JP) <= DS(JB)) then
                        exit
                    end if
                end do
                JBDPU(JP) = JB
                do JW = 1, NWB
                    if (JB >= BS(JW) .and. JB <= BE(JW)) then
                        exit
                    end if
                end do
                JWDPU(JP) = JW
            else
                JBDPU(JP) = 1
                JWDPU(JP) = 1
            end if
        end do
    end if

    allocate(ESTR(NST, NBR), WSTR(NST, NBR), QSTR(NST, NBR), KTSW(NST, NBR), KBSW(NST, NBR), SINKC(NST, NBR), POINT_SINK(NST, NBR), QNEW(KMX), tavg(nst, nbr), tavgw(NWD + NSP + NGT + NPI + NPU), CAVG(NST, NBR, NCT), CDAVG(NST, NBR, NDC), CAVGW(NWD + NSP + NGT + NPI + NPU, NCT), CDAVGW(NWD + NSP + NGT + NPI + NPU, NDC))
    allocate(ACTIVE_RULE_W2SELECTIVE(NST, NBR));     ACTIVE_RULE_W2SELECTIVE = .false.
    TAVGW = 0.0
    TAVG = 0.0
    CAVG = 0.0
    CDAVG = 0.0
    CAVGW = 0.0
    CDAVGW = 0.0

    QSTR = 0.0
    do JB = 1, NBR
        ESTR(1:NSTR(JB), JB) = ESTRT(1:NSTR(JB), JB)
        KTSW(1:NSTR(JB), JB) = KTSWT(1:NSTR(JB), JB)
        KBSW(1:NSTR(JB), JB) = KBSWT(1:NSTR(JB), JB)
        WSTR(1:NSTR(JB), JB) = WSTRT(1:NSTR(JB), JB)
        SINKC(1:NSTR(JB), JB) = SINKCT(1:NSTR(JB), JB)
        POINT_SINK(1:NSTR(JB), JB) = SINKC(1:NSTR(JB), JB) == "   POINT" ! SW 9/27/13
    end do
    deallocate(ESTRT, KBSWT, KTSWT, WSTRT, SINKCT)

! Active constituents, derived constituents, and fluxes

    if (CONSTITUENTS) then
        do JC = 1, NCT
            if (CAC(JC) == "      ON") then
                NAC = NAC + 1
                CN(NAC) = JC
            end if
            do JB = 1, NBR
                if (CINBRC(JC, JB) == "      ON") then
                    NACIN(JB) = NACIN(JB) + 1
                    INCN(NACIN(JB), JB) = JC
                end if
                if (CDTBRC(JC, JB) == "      ON") then
                    NACDT(JB) = NACDT(JB) + 1
                    DTCN(NACDT(JB), JB) = JC
                end if
                if (CPRBRC(JC, JB) == "      ON") then
                    NACPR(JB) = NACPR(JB) + 1
                    PRCN(NACPR(JB), JB) = JC
                end if
            end do
            do JT = 1, NTR
                if (CTRTRC(JC, JT) == "      ON") then
                    NACTR(JT) = NACTR(JT) + 1
                    TRCN(NACTR(JT), JT) = JC
                end if
            end do
        end do
        do JW = 1, NWB
            do JD = 1, NDC
                if (CDWBC(JD, JW) == "      ON") then
                    NACD(JW) = NACD(JW) + 1
                    CDN(NACD(JW), JW) = JD
                end if
            end do
            do JF = 1, NFL
                if (KFWBC(JF, JW) == "      ON") then
                    NAF(JW) = NAF(JW) + 1
                    KFCN(NAF(JW), JW) = JF
                else
                    if (PH_CALC(jw) .and. JF == KF_CO2X) then
                        NAF(JW) = NAF(JW) + 1
                        KFCN(NAF(JW), JW) = JF
                    else
                        if (CAC(NH2S) == "      ON" .and. JF >= KF_DOH2S .and. JF < KF_DOCH4) then
                            NAF(JW) = NAF(JW) + 1
                            KFCN(NAF(JW), JW) = JF
                        else
                            if (CAC(NCH4) == "      ON" .and. JF >= KF_DOCH4 .and. JF < KF_FE2D) then
                                NAF(JW) = NAF(JW) + 1
                                KFCN(NAF(JW), JW) = JF
                            else
                                if (CAC(NFEII) == "      ON" .and. JF >= KF_FE2D .and. JF < KF_MN2D) then
                                    NAF(JW) = NAF(JW) + 1
                                    KFCN(NAF(JW), JW) = JF
                                else
                                    if (CAC(NMNII) == "      ON" .and. JF >= KF_MN2D .and. JF < KF_SDINC) then
                                        NAF(JW) = NAF(JW) + 1
                                        KFCN(NAF(JW), JW) = JF
                                    else
                                        if (JF >= KF_SDINC .and. JF < KF_SEDD .and. CEMARelatedCode) then ! CEMA turning on flux output FOR SEDIMENT DIAGENESIS
                                            NAF(JW) = NAF(JW) + 1
                                            KFCN(NAF(JW), JW) = JF
                                        else
                                            if (JF >= KF_SEDD .and. STANDING_BIOMASS_DECAY) then ! STANDING ORGANIC MATTER
                                                NAF(JW) = NAF(JW) + 1
                                                KFCN(NAF(JW), JW) = JF
                                            end if
                                        end if
                                    end if
                                end if
                            end if
                        end if
                    end if
                end if
            end do
            if (ATM_DEPOSITION(JW)) then
                do JC = 1, NCT
                    if (C_ATM_DEPOSITION(JC, JW) == "      ON") then
                        NACATD(JW) = NACATD(JW) + 1
                        ATMDCN(NACATD(JW), JW) = JC
                    end if
                end do
            end if
        end do ! JW LOOP
    end if

! Starting time

    DEG = CHAR(248) // "C"
    ESC = CHAR(27)
    call DATE_AND_TIME(CDATE, CCTIME)
    do JW = 1, NWB
        TITLE(11) = "Model run at " // CCTIME(1:2) // ":" // CCTIME(3:4) // ":" // CCTIME(5:6) // " on " // CDATE(5:6) // "/" // CDATE(7:8) // "/" // CDATE(3:4)
        if (RESTART_IN) then
            TITLE(11) = "Model restarted at " // CCTIME(1:2) // ":" // CCTIME(3:4) // ":" // CCTIME(5:6) // " on " // CDATE(5:6) // "/" // CDATE(7:8) // "/" // CDATE(3:4)
        end if
    end do

    call INITGEOM() ! Call initial geometry

! Density related derived constants

    RHOWCP = RHOW*CP
    RHOICP = RHOI*CP
    RHOIRL1 = RHOI*RL1
    DLXRHO = 0.0D0
    do JW = 1, NWB
        do JB = BS(JW), BE(JW)
            DLXRHO(US(JB):DS(JB)) = 0.5D0/(DLXR(US(JB):DS(JB))*RHOW)
            if (UP_HEAD(JB)) then
                DLXRHO(US(JB) - 1) = 0.5D0/(DLXR(US(JB))*RHOW)
            end if
        end do
    end do

! Transport interpolation multipliers

    do JW = 1, NWB
        call INTERPOLATION_MULTIPLIERS()
    end do

    call INITCOND() ! CALL ROUTINE TO SET UP IC

    if (SD_GLOBAL) then
        call INIT_CEMA()
    end if ! SW 2/18/2019
    if (IncludeCEMASedDiagenesis) then
        call InitCond_SedFlux()
    end if

    if (WEIR_CALC) then ! MOVED FROM ABOVE AFTER GEOMETRY SETUP  SW 3/16/18
        do JWR = 1, NIW
            if (EKTWR(JWR) == 0.0) then
                do JW = 1, NWB
                    if (IWR(JWR) >= US(BS(JW)) .and. IWR(JWR) <= DS(BE(JW))) then
                        KTWR(JWR) = KTWB(JW)
                        exit
                    end if
                end do
            else
                KTWR(JWR) = INT(EKTWR(JWR))
            end if
            if (EKBWR(JWR) <= 0.0) then
                do K = KTWR(JWR), KB(IWR(JWR))
                    if (DEPTHB(K, IWR(JWR)) >= ABS(EKBWR(JWR))) then
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
        end do
    end if


! Saved variables for autostepping

    if (.not. RESTART_IN) then
        SZ = Z
        SU = U
        SW = W
        SAZ = AZ
        SKTI = KTI
        SBKT = BKT
        SAVH2 = AVH2
        SAVHR = AVHR
    end if
    call GREGORIAN_DATE()
!  CALL TIME_VARYING_DATA
!  if(.not. once_through)
    call TIME_VARYING_DATA() ! w2-ressim
!  CALL READ_INPUT_DATA (NXTVD)
!  if(.not. once_through)
    call READ_INPUT_DATA(NXTVD) ! w2-ressim
!IF (CONSTITUENTS) THEN   !SW 6/26/2019 No need to compute since done at the beginning of WQ for NIT=0
!  DO JW=1,NWB
!    KT = KTWB(JW)
!    DO JB=BS(JW),BE(JW)
!      IU = US(JB)
!      ID = DS(JB)
!      CALL TEMPERATURE_RATES
!      CALL KINETIC_RATES
!    END DO
!  END DO
!END IF

    return
end subroutine INIT
