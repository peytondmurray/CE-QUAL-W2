subroutine TIME_VARYING_DATA()
    use GLOBAL;     use SURFHE;     use SCREENC;     use TVDC;     use LOGICC;     use SELWC;     use STRUCTURES;     use NAMESC
    use KINETIC, only: EXH2O;     use SHADEC;     use MAIN, only: PUMPS, SEGNUM, WBSEG, N2BND, DOBND, DGPBND, NDO, NDGP, ITR, EA, SYSTDG, NN2, TMEND, ATM_DEPOSITION, ATMDEPFN, ATM_DEP_LOADING, ATMDCN, NACATD, ATM_DEPOSITION_INTERPOLATION, GASGTC, GASSPC ! systdg ADD WBSEG, N2BND, DOBND, NDO, NGN2, ITR, EA, SYSTDG
    use modSYSTDG, only: TWETSC, TWEFN, TWE_TS;     use TDGAS ! systdg
    use IFPORT ! for SLEEPQQ and SYSTEMQQ commands                                       !SR 11/28/19
    implicit none

! Type declaration

    character(len=1) :: INFORMAT
    character(len=2) :: INFORMAT2
    real :: NXQGT2
    real :: NXQWD1, NXQWD2, NXQGT, NXTVD, NXQPT
    real :: NXWSC
    real(R8) :: RATIO, QRATIO, TRATIO, CRATIO, HRATIO
    real(R8), allocatable, dimension(:) :: QDTRO, TDTRO, ELUHO, ELDHO, QWDO, QTRO, TTRO, QINO, TINO
    real, allocatable, dimension(:) :: TAIRNX, TDEWNX, PHINX, WINDNX, SRONX, CLOUDNX, BGTNX, XX, PALT_JWNX ! systdg - PALT_JWNX
    real, allocatable, dimension(:) :: NXEXT1, NXEXT2, EXTNX, EXTO, NXATMD, NXATMD2
    real, allocatable, dimension(:) :: TAIRO, TDEWO, PHIO, WINDO, SROO, CLOUDO, PALT_JWO ! systdg - PALT_JWO
    real(R8), allocatable, dimension(:) :: QDTRNX, TDTRNX, PRNX, TPRNX, ELUHNX, ELDHNX, QWDNX, QTRNX, TTRNX, QINNX, TINNX
    real, allocatable, dimension(:) :: NXQTR1, NXTTR1, NXCTR1, NXQIN1, NXTIN1, NXCIN1, NXQDT1, NXTDT1, NXCDT1, NXDYNS
    real, allocatable, dimension(:) :: NXPR1, NXTPR1, NXCPR1, NXEUH1, NXTUH1, NXCUH1, NXEDH1, NXTDH1, NXCDH1, NXQOT1, NXMET1
    real, allocatable, dimension(:) :: NXQTR2, NXTTR2, NXCTR2, NXQIN2, NXTIN2, NXCIN2, NXQDT2, NXTDT2, NXCDT2, NXSPDO, NXGTDO
    real, allocatable, dimension(:) :: NXPR2, NXTPR2, NXCPR2, NXEUH2, NXTUH2, NXCUH2, NXEDH2, NXTDH2, NXCDH2, NXQOT2, NXMET2
    real, allocatable, dimension(:) :: WSCNX, BPNX, EPU2, EONPU2, EOFFPU2, QPU2, NXPUMP, NXATMDEP, AGASSPNX, AGASGTNX
    real, allocatable, dimension(:,:) :: NXESTRT, ATM_DEP_LOADINGNX, ATM_DEP_LOADING0
    real(R8), allocatable, dimension(:,:) :: CTRO, CINO, QOUTO, CDTRO, TUHO, TDHO, QSTRO
    real(R8), allocatable, dimension(:,:) :: CTRNX, CINNX, QOUTNX, CDTRNX, CPRNX, TUHNX, TDHNX, QSTRNX
    real, allocatable, dimension(:,:,:) :: CUHO, CDHO, CUHNX, CDHNX
    integer :: WDQ, GTQ, WSH, SHD, PIPED, IOPENPIPE, L, NJS
    real(R8) :: TWERATIO, NXTWE1, NXTWE2, TWE_TSNX, TWE_TSO ! systdg - TWE
    integer :: TWEFNNO ! systdg - TWE
    logical :: TWEF ! systdg - TWE
    real(R8), allocatable, dimension(:) :: DO_SAT, N2_SAT, DO_SATJ, N2_SATJ, DO_SATD, N2_SATD, DO_SATP, N2_SATP ! systdg - DO_SATJ, N2_SATJ, DO_SATD, N2_SATD, DO_SATP, N2_SATP 
    integer :: NPT, J, JT, JAC, JS, K, JG, JWD
    integer, allocatable, dimension(:) :: TRQ, TRT, TRC, INQ, DTQ, PRE, UHE, DHE, INFT, DTT, PUMPD, JJS, ATMDEP, FGASSP, FGASGT
    integer, allocatable, dimension(:) :: PRT, UHT, DHT, INC, DTC, PRC, UHC, DHC, OTQ, MET, EXT, ODYNS, DYNPUMPF
    logical, allocatable, dimension(:) :: INFLOW_CONST, TRIB_CONST, DTRIB_CONST, PRECIP_CONST, OTQF, TRCF, DTCF, INCF, PRCF, METF, DYNEF, TRQF, TRTF, DTTF, DTQF, INQF, INTF, PRQF, PRTF, EXTF
    logical :: WDQF, WSHF, GATEF, ATMDEPCSV ! SW 9/26/2017
    integer, allocatable, dimension(:) :: EUHF, TUHF, CUHF, EDHF, TDHF, CDHF ! =0 Old format for head BCs, =1 Time series format no vertical variation, =2 csv format vertical variation             SW 2/28/17
    character(len=240), allocatable, dimension(:) :: FILE_GAS_GT, FILE_GAS_SP

!RESULT1               -- integer variable holding success or failure result from SYSTEMQQ command
!ITER                  -- integer variable holding the number of iterations for an awaited file read (see code)
!LAST_JDAY             -- real variable typically used to hold the last available date in awaited input file
!FULL_FILE_NAME        -- character string to hold the full path and file name of interest
!INPUT_FILE_EXISTS     -- logical variable holding the results of an inquiry into a file's existence
!
!WAIT_TIME             -- integer specifying how many days of data to await in the awaited input files
!TIME_BUFFER           -- integer number of seconds to wait each time the model waits for new input
!WAIT_FOR_TRIB_INPUT   -- logical array (trib index) used to determine which tributary has awaited input
!WAIT_FOR_BRANCH_INPUT -- logical array (branch index) used to determine which branch has awaited input
!TR_FILEDIR            -- character array (tributary index) to hold the directory names of any awaited tributary input files
!BR_FILEDIR            -- character array (branch index) to hold the directory names of any awaited branch input files

    integer :: RESULT1, ITER, N !SR 11/28/19
    real :: LAST_JDAY, GET_LAST_JDAY !SR 11/28/19
    character(len=240) :: FULL_FILE_NAME !SR 11/28/19
    logical :: INPUT_FILE_EXISTS !SR 11/28/19
    external :: GET_LAST_JDAY !SR 11/28/19
    external :: PRINT_ERROR_AND_STOP !SR 11/28/19

    save
! Allocation declarations

    if (NAC < 1) then
        allocate(XX(1))
    else
        allocate(XX(NCT))
    end if
!
! systdg - time series input
    allocate(DO_SAT(NTR), N2_SAT(NTR))
    allocate(DO_SATJ(NBR), N2_SATJ(NBR))
    allocate(DO_SATD(NBR), N2_SATD(NBR))
    allocate(DO_SATP(NBR), N2_SATP(NBR))
!
    allocate(NXQTR1(NTR), NXTTR1(NTR), NXCTR1(NTR), NXQIN1(NBR), NXTIN1(NBR), NXCIN1(NBR), NXQDT1(NBR), NXTDT1(NBR), NXCDT1(NBR))
    allocate(NXPR1(NBR), NXTPR1(NBR), NXCPR1(NBR), NXEUH1(NBR), NXTUH1(NBR), NXCUH1(NBR), NXEDH1(NBR), NXTDH1(NBR), NXCDH1(NBR))
    allocate(NXQOT1(NBR), NXMET1(NWB), NXQTR2(NTR), NXTTR2(NTR), NXCTR2(NTR), NXQIN2(NBR), NXTIN2(NBR), NXCIN2(NBR), NXQDT2(NBR))
    allocate(NXTDT2(NBR), NXCDT2(NBR), NXPR2(NBR), NXTPR2(NBR), NXCPR2(NBR), NXEUH2(NBR), NXTUH2(NBR), NXCUH2(NBR), NXEDH2(NBR))
    allocate(NXTDH2(NBR), NXCDH2(NBR), NXQOT2(NBR), NXMET2(NWB), DYNPUMPF(NPU))
    allocate(WSCNX(IMX), PRCF(NBR), METF(NWB), ODYNS(NBR), NXDYNS(NBR), NXESTRT(NST, NBR), PRTF(NBR), PRQF(NBR), EXTF(NWB))
    allocate(QDTRO(NBR), TDTRO(NBR), ELUHO(NBR), ELDHO(NBR), QWDO(NWD), QTRO(NTR), TTRO(NTR), QINO(NBR), TINO(NBR))
    allocate(QDTRNX(NBR), TDTRNX(NBR), PRNX(NBR), TPRNX(NBR), ELUHNX(NBR), ELDHNX(NBR), QWDNX(NWD), QTRNX(NTR), TTRNX(NTR))
    allocate(QINNX(NBR), TINNX(NBR), SROO(NWB), TAIRO(NWB), TDEWO(NWB), CLOUDO(NWB), PHIO(NWB), WINDO(NWB), TAIRNX(NWB), PALT_JWO(NWB), PALT_JWNX(NWB)) ! systdg - time series PALT
    allocate(TDEWNX(NWB), CLOUDNX(NWB), PHINX(NWB), WINDNX(NWB), SRONX(NWB), BGTNX(NGT), BPNX(NPI))
    allocate(TRQ(NTR), TRT(NTR), TRC(NTR), INQ(NBR), DTQ(NBR), PRE(NBR), UHE(NBR), DHE(NBR), INFT(NBR))
    allocate(DTT(NBR), PRT(NBR), UHT(NBR), DHT(NBR), INC(NBR), DTC(NBR), PRC(NBR), UHC(NBR), DHC(NBR))
    allocate(OTQ(NBR), MET(NWB), EXT(NWB), NXATMDEP(NWB))
    allocate(NXEXT1(NWB), NXEXT2(NWB), EXTNX(NWB), EXTO(NWB))
    allocate(CTRO(NCT, NTR), CINO(NCT, NBR), QOUTO(KMX, NBR), CDTRO(NCT, NBR), TUHO(KMX, NBR), TDHO(KMX, NBR), QSTRO(NST, NBR))
    allocate(CTRNX(NCT, NTR), CINNX(NCT, NBR), QOUTNX(KMX, NBR), CDTRNX(NCT, NBR), CPRNX(NCT, NBR), TUHNX(KMX, NBR), TDHNX(KMX, NBR))
    allocate(QSTRNX(NST, NBR), PUMPD(NPU), NXPUMP(NPU), EPU2(NPU), EONPU2(NPU), EOFFPU2(NPU), QPU2(NPU), DYNEF(NBR))
    allocate(CUHO(KMX, NCT, NBR), CDHO(KMX, NCT, NBR), CUHNX(KMX, NCT, NBR), CDHNX(KMX, NCT, NBR), JJS(NST))
    allocate(INFLOW_CONST(NBR), TRIB_CONST(NTR), DTRIB_CONST(NBR), PRECIP_CONST(NBR), OTQF(NBR), TRCF(NTR), DTCF(NBR), INCF(NBR), TRQF(NTR), TRTF(NTR), DTTF(NBR), DTQF(NBR), INQF(NBR), INTF(NBR))
    allocate(EUHF(NBR), TUHF(NBR), CUHF(NBR), EDHF(NBR), TDHF(NBR), CDHF(NBR))
    allocate(ATM_DEP_LOADINGNX(NCT, NWB), ATMDEP(NWB), NXATMD(NWB), NXATMD2(NWB), ATM_DEP_LOADING0(NCT, NWB))

    if (NGT > 0) then
        allocate(NXGTDO(NGT), FGASGT(NGT), FILE_GAS_GT(NGT), AGASGTNX(NGT))
    end if
    if (NSP > 0) then
        allocate(NXSPDO(NSP), FGASSP(NSP), FILE_GAS_SP(NSP), AGASSPNX(NSP))
    end if

    NXPR1 = 0.0;     NXQTR1 = 0.0;     NXTTR1 = 0.0;     NXCTR1 = 0.0;     NXQIN1 = 0.0;     NXTIN1 = 0.0;     NXCIN1 = 0.0;     NXQDT1 = 0.0;     NXTDT1 = 0.0
    NXCDT1 = 0.0;     NXTPR1 = 0.0;     NXCPR1 = 0.0;     NXEUH1 = 0.0;     NXTUH1 = 0.0;     NXCUH1 = 0.0;     NXEDH1 = 0.0;     NXTDH1 = 0.0;     NXCDH1 = 0.0
    NXQOT1 = 0.0;     NXMET1 = 0.0;     QSTRNX = 0.0;     CDTRNX = 0.0;     CTRNX = 0.0;     CINNX = 0.0;     CPRNX = 0.0;     CUHNX = 0.0;     CDHNX = 0.0
    QINNX = 0.0;     TINNX = 0.0;     CINNX = 0.0;     NXWSC = 0.0;     NXATMDEP = 0.0;     ATM_DEP_LOADINGNX = 0.0
    NXTWE1 = 0.0;     TWE_TSNX = 0.0;  ! systdg - time series TWE
    TWEF = .false. ! systdg - time series TWE
! Set logical variables

    INFLOW_CONST = CONSTITUENTS .and. NACIN > 0;     TRIB_CONST = CONSTITUENTS .and. NACTR > 0
    DTRIB_CONST = CONSTITUENTS .and. NACDT > 0;     PRECIP_CONST = CONSTITUENTS .and. NACPR > 0
    OTQF = .false.
    WSHF = .false.
    WDQF = .false.
    TRCF = .false.
    TRQF = .false.
    TRTF = .false.
    DTCF = .false.
    INCF = .false.
    PRCF = .false.
    PRTF = .false.
    PRQF = .false.
    METF = .false.
    DYNEF = .false.
    DTTF = .false.
    DTQF = .false.
    INQF = .false.
    INTF = .false.
    EXTF = .false.
    DYNPUMPF = .false.;     GATEF = .false. ! SW 9/26/2017
    EUHF = 0;     TUHF = 0;     CUHF = 0;     EDHF = 0;     TDHF = 0;     CDHF = 0 ! SW 2/28/17

! Open input files

    NPT = NUNIT
    SHD = NPT;     NPT = NPT + 1
    open(SHD, FILE=SHDFN, STATUS="OLD")
    read(SHD, "(A1)") INFORMAT
    if (INFORMAT == "$") then
        read(SHD, "(/)")
        do I = 1, IMX
            read(SHD, *) J, SHADEI(I) ! SW 3/14/2018 ADDED TO BE COMPATIBLE WITH PREPROCESSOR
            if (SHADEI(I) < 0.0) then
                backspace(SHD)
                read(SHD, *) J, SHADEI(I), TTLB(I), TTRB(I), CLLB(I), CLRB(I), SRLB1(I), SRLB2(I), SRRB1(I), SRRB2(I), (TOPO(I, J), J = 1, IANG), SRFJD1(I), SRFJD2(I)
            end if
        end do
    else
        read(SHD, "(//(8X,29F8.0))") (SHADEI(I), TTLB(I), TTRB(I), CLLB(I), CLRB(I), SRLB1(I), SRLB2(I), SRRB1(I), SRRB2(I), (TOPO(I, J), J = 1, IANG), SRFJD1(I), SRFJD2(I), I = 1, IMX)
    end if
    SHADE = SHADEI
    WSH = NPT;     NPT = NPT + 1
    open(WSH, FILE=WSCFN, STATUS="OLD")
    read(WSH, "(A1)") INFORMAT
    if (INFORMAT == "$") then
        WSHF = .true.
    end if
    if (WSHF) then
        read(WSH, "(/)")
        read(WSH, *) NXWSC, (WSCNX(I), I = 1, IMX)
        WSC = WSCNX
        read(WSH, *) NXWSC, (WSCNX(I), I = 1, IMX)
    else
        read(WSH, "(//10F8.0:/(8X,9F8.0))") NXWSC, (WSCNX(I), I = 1, IMX)
        WSC = WSCNX
        read(WSH, "(10F8.0:/(8X,9F8.0))") NXWSC, (WSCNX(I), I = 1, IMX)
    end if
    do JW = 1, NWB
        MET(JW) = NPT;         NPT = NPT + 1
        open(MET(JW), FILE=METFN(JW), STATUS="OLD")
        read(MET(JW), "(A1)") INFORMAT
        if (INFORMAT == "$") then
            METF(JW) = .true.
        end if
        if (SYSTDG) then
            if (READ_RADIATION(JW)) then
                if (METF(JW)) then
                    read(MET(JW), "(/)")
                    read(MET(JW), *) NXMET2(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW), SRONX(JW), PALT_JWNX(JW) ! systdg - time series input PALT
                else
                    read(MET(JW), "(//8F8.0)") NXMET2(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW), SRONX(JW), PALT_JWNX(JW) ! systdg - time series input PALT
                end if
                SRONX(JW) = SRONX(JW)*REFL
                SRON(JW) = SRONX(JW)
                SROO(JW) = SRON(JW)
            else
                if (METF(JW)) then
                    read(MET(JW), "(/)")
                    read(MET(JW), *) NXMET2(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW), PALT_JWNX(JW) ! systdg - time series input PALT      
                else
                    read(MET(JW), "(//7F8.0)") NXMET2(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW), PALT_JWNX(JW) ! systdg - time series input PALT     
                end if
            end if
        else
            if (READ_RADIATION(JW)) then
                if (METF(JW)) then
                    read(MET(JW), "(/)")
                    read(MET(JW), *) NXMET2(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW), SRONX(JW)
                else
                    read(MET(JW), "(//8F8.0)") NXMET2(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW), SRONX(JW)
                end if
                SRONX(JW) = SRONX(JW)*REFL
                SRON(JW) = SRONX(JW)
                SROO(JW) = SRON(JW)
            else
                if (METF(JW)) then
                    read(MET(JW), "(/)")
                    read(MET(JW), *) NXMET2(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW)
                else
                    read(MET(JW), "(//7F8.0)") NXMET2(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW)
                end if
            end if

        end if

        TAIR(JW) = TAIRNX(JW)
        TDEW(JW) = TDEWNX(JW)
        WIND(JW) = WINDNX(JW)
        PHI(JW) = PHINX(JW)
        CLOUD(JW) = CLOUDNX(JW)
        if (SYSTDG) then
            PALT_JW(JW) = PALT_JWNX(JW) ! systdg - time series input PALT
            if (PALT_JW(JW) <= 0.0) then
                PALT_JW(JW) = 760.0
            end if ! systdg - time series input PALT
            do I = US(BS(JW)) - 1, DS(BE(JW)) + 1
                PALT(I) = PALT_JW(JW)/760.0*(1.0 - ELWS_INI(I)/1000.0/44.3)**5.25 ! systdg - time series input PALT
            end do
            PALT_JWO(JW) = PALT_JWNX(JW) ! systdg - time series input PALT
        end if
        TAIRO(JW) = TAIRNX(JW)
        TDEWO(JW) = TDEWNX(JW)
        WINDO(JW) = WINDNX(JW)
        PHIO(JW) = PHINX(JW)
        CLOUDO(JW) = CLOUDNX(JW)
!IF (PHISET > 0) PHI(JW)  = PHISET
!IF (PHISET > 0) PHIO(JW) = PHISET
        do N = 1, NSP ! SW 1/18/2022
            if (GASSPC(N) == "      ON" .and. EQSP(N) == 4 .and. BGASSP(N) == 1) then
                FGASSP(N) = NPT;                 NPT = NPT + 1
                write(SEGNUM, "(I0)") N
                SEGNUM = ADJUSTL(SEGNUM)
                L = LEN_TRIM(SEGNUM)
                FILE_GAS_SP(N) = "w2_sp" // SEGNUM(1:L) // "DO.csv"
                open(FGASSP(N), FILE=FILE_GAS_SP(N), STATUS="OLD")
                read(FGASSP(N), *);                 read(FGASSP(N), *);                 read(FGASSP(N), *) ! SKIP 3 LINES
                read(FGASSP(N), *) NXSPDO(N), AGASSPNX(N)
                AGASSP(N) = AGASSPNX(N)
            end if
        end do
        do N = 1, NGT
            if (GASGTC(N) == "      ON" .and. EQGT(N) == 4 .and. BGASGT(N) == 1.0) then
                FGASGT(N) = NPT;                 NPT = NPT + 1
                write(SEGNUM, "(I0)") N
                SEGNUM = ADJUSTL(SEGNUM)
                L = LEN_TRIM(SEGNUM)
                FILE_GAS_GT(N) = "w2_gt" // SEGNUM(1:L) // "DO.csv"
                open(FGASGT(N), FILE=FILE_GAS_GT(N), STATUS="OLD")
                read(FGASGT(N), *);                 read(FGASGT(N), *);                 read(FGASGT(N), *) ! SKIP 3 LINES
                read(FGASGT(N), *) NXGTDO(N), AGASGTNX(N)
                AGASGT(N) = AGASGTNX(N)
            end if
        end do

        if (SYSTDG) then
            if (READ_RADIATION(JW)) then
                if (METF(JW)) then
                    read(MET(JW), *) NXMET1(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW), SRONX(JW), PALT_JWNX(JW) ! systdg - time series input PALT
                else
                    read(MET(JW), "(8F8.0)") NXMET1(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW), SRONX(JW), PALT_JWNX(JW) ! systdg - time series input PALT
                end if
                SRONX(JW) = SRONX(JW)*REFL
            else
                if (METF(JW)) then
                    read(MET(JW), *) NXMET1(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW), PALT_JWNX(JW) ! systdg - time series input PALT
                else
                    read(MET(JW), "(7F8.0)") NXMET1(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW), PALT_JWNX(JW) ! systdg - time series input PALT
                end if
            end if
        else
            if (READ_RADIATION(JW)) then
                if (METF(JW)) then
                    read(MET(JW), *) NXMET1(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW), SRONX(JW)
                else
                    read(MET(JW), "(8F8.0)") NXMET1(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW), SRONX(JW)
                end if
                SRONX(JW) = SRONX(JW)*REFL
            else
                if (METF(JW)) then
                    read(MET(JW), *) NXMET1(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW)
                else
                    read(MET(JW), "(7F8.0)") NXMET1(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW)
                end if
            end if
        end if

        if (READ_EXTINCTION(JW)) then
            EXT(JW) = NPT;             NPT = NPT + 1
            open(EXT(JW), FILE=EXTFN(JW), STATUS="OLD")
            read(EXT(JW), "(A1)") INFORMAT
            if (INFORMAT == "$") then
                EXTF(JW) = .true.
            end if
            if (EXTF(JW)) then
                read(EXT(JW), "(/)")
                read(EXT(JW), *) NXEXT2(JW), EXTNX(JW)
            else
                read(EXT(JW), "(///2F8.0)") NXEXT2(JW), EXTNX(JW)
            end if

            EXH2O(JW) = EXTNX(JW)
            EXTO(JW) = EXTNX(JW)

            if (EXTF(JW)) then
                read(EXT(JW), *) NXEXT1(JW), EXTNX(JW)
            else
                read(EXT(JW), "(2F8.0)") NXEXT1(JW), EXTNX(JW)
            end if
        end if
!DO I=CUS(BS(JW)),DS(BE(JW))   ! SW CODE FIX 5-21-15
!  WIND2(I) = WIND(JW)*WSC(I)*DLOG(2.0D0/Z0(JW))/DLOG(WINDH(JW)/Z0(JW))
!END DO
    end do
    if (NWD > 0) then
        WDQ = NPT;         NPT = NPT + 1
        open(WDQ, FILE=QWDFN, STATUS="OLD")
        read(WDQ, "(A1)") INFORMAT
        if (INFORMAT == "$") then
            WDQF = .true.
        end if
        if (WDQF) then
            read(WDQ, "(/)")
            read(WDQ, *) NXQWD2, (QWDNX(JW), JW = 1, NWD)
            do JW = 1, NWD
                QWD(JW) = QWDNX(JW)
                QWDO(JW) = QWDNX(JW)
            end do
            read(WDQ, *) NXQWD1, (QWDNX(JW), JW = 1, NWD)
        else
            read(WDQ, "(//10F8.0:/(8X,9F8.0))") NXQWD2, (QWDNX(JW), JW = 1, NWD)
            do JW = 1, NWD
                QWD(JW) = QWDNX(JW)
                QWDO(JW) = QWDNX(JW)
            end do
            read(WDQ, "(10F8.0:/(8X,9F8.0))") NXQWD1, (QWDNX(JW), JW = 1, NWD)
        end if
    end if
    if (TRIBUTARIES) then
        do JT = 1, NTR
            TRQ(JT) = NPT;             NPT = NPT + 1
            TRT(JT) = NPT;             NPT = NPT + 1

            if (WAIT_FOR_TRIB_INPUT(JT)) then ! Wait for input from tributary flow file !SR 11/28/19
                FULL_FILE_NAME = TRIM(ADJUSTL(TR_FILEDIR(JT))) // "\" // TRIM(ADJUSTL(QTRFN(JT))) !SR 11/28/19
                inquire(FILE=TRIM(FULL_FILE_NAME), EXIST=INPUT_FILE_EXISTS) ! Check that file exists                  !SR 11/28/19
                do while (.not. INPUT_FILE_EXISTS) !SR 11/28/19
                    write(*, "(3A)") "Input file ", TRIM(FULL_FILE_NAME), " does not exist (yet).  Waiting..." !SR 11/28/19
                    write(9911, "(3A)") "Input file ", TRIM(FULL_FILE_NAME), " does not exist (yet).  Waiting..." !SR 11/28/19
                    call SLEEPQQ(WAIT_TIME*1000) !SR 11/28/19
                    inquire(FILE=TRIM(FULL_FILE_NAME), EXIST=INPUT_FILE_EXISTS) !SR 11/28/19
                end do !SR 11/28/19
                RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory        !SR 11/28/19
                write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                if (.not. RESULT1) then
                    call PRINT_ERROR_AND_STOP("QTR", JT)
                end if ! Problem with copy; write msg and stop   !SR 11/28/19

                LAST_JDAY = GET_LAST_JDAY(QTRFN(JT)) ! Find last JDAY in input file            !SR 11/28/19
                do while (LAST_JDAY <= TMEND - 0.5 .and. JDAY > LAST_JDAY - TIME_BUFFER) ! Not enough data in input file           !SR 11/28/19
                    write(9911, "(A,I0,3(A,F0.4))") "WAIT: Input QTR", JT, " DAY= ", LAST_JDAY, " JDAY= ", JDAY, " TMEND= ", TMEND !SR 11/28/19
                    call SLEEPQQ(WAIT_TIME*1000) !SR 11/28/19
                    RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory        !SR 11/28/19
                    write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                    if (.not. RESULT1) then
                        call PRINT_ERROR_AND_STOP("QTR", JT)
                    end if ! Problem with copy; write msg and stop   !SR 11/28/19
                    LAST_JDAY = GET_LAST_JDAY(QTRFN(JT)) !SR 11/28/19
                end do !SR 11/28/19
            end if !SR 11/28/19
            open(TRQ(JT), FILE=QTRFN(JT), STATUS="OLD")
            read(TRQ(JT), "(A1)") INFORMAT
            if (INFORMAT == "$") then
                TRQF(JT) = .true.
            end if
            if (TRQF(JT)) then
                read(TRQ(JT), "(/)")
                read(TRQ(JT), *) NXQTR2(JT), QTRNX(JT)
            else
                read(TRQ(JT), "(//2F8.0)") NXQTR2(JT), QTRNX(JT)
            end if

            if (WAIT_FOR_TRIB_INPUT(JT)) then ! Wait for input from tributary temp file !SR 11/28/19
                FULL_FILE_NAME = TRIM(ADJUSTL(TR_FILEDIR(JT))) // "\" // TRIM(ADJUSTL(TTRFN(JT))) !SR 11/28/19
                inquire(FILE=TRIM(FULL_FILE_NAME), EXIST=INPUT_FILE_EXISTS) ! Check that file exists                  !SR 11/28/19
                do while (.not. INPUT_FILE_EXISTS) !SR 11/28/19
                    write(*, "(3A)") "Input file ", TRIM(FULL_FILE_NAME), " does not exist (yet).  Waiting..." !SR 11/28/19
                    write(9911, "(3A)") "Input file ", TRIM(FULL_FILE_NAME), " does not exist (yet).  Waiting..." !SR 11/28/19
                    call SLEEPQQ(WAIT_TIME*1000) !SR 11/28/19
                    inquire(FILE=TRIM(FULL_FILE_NAME), EXIST=INPUT_FILE_EXISTS) !SR 11/28/19
                end do !SR 11/28/19
                RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory        !SR 11/28/19
                write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                if (.not. RESULT1) then
                    call PRINT_ERROR_AND_STOP("TTR", JT)
                end if ! Problem with copy; write msg and stop   !SR 11/28/19

                LAST_JDAY = GET_LAST_JDAY(TTRFN(JT)) ! Find last JDAY in input file            !SR 11/28/19
                do while (LAST_JDAY <= TMEND - 0.5 .and. JDAY > LAST_JDAY - TIME_BUFFER) ! Not enough data in input file           !SR 11/28/19
                    write(9911, "(A,I0,3(A,F0.4))") "WAIT: Input TTR", JT, " DAY= ", LAST_JDAY, " JDAY= ", JDAY, " TMEND= ", TMEND !SR 11/28/19
                    call SLEEPQQ(WAIT_TIME*1000) !SR 11/28/19
                    RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory        !SR 11/28/19
                    write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                    if (.not. RESULT1) then
                        call PRINT_ERROR_AND_STOP("TTR", JT)
                    end if ! Problem with copy; write msg and stop   !SR 11/28/19
                    LAST_JDAY = GET_LAST_JDAY(TTRFN(JT)) !SR 11/28/19
                end do !SR 11/28/19
            end if !SR 11/28/19
            open(TRT(JT), FILE=TTRFN(JT), STATUS="OLD")
            read(TRT(JT), "(A1)") INFORMAT
            if (INFORMAT == "$") then
                TRTF(JT) = .true.
            end if
            if (TRTF(JT)) then
                read(TRT(JT), "(/)")
                read(TRT(JT), *) NXTTR2(JT), TTRNX(JT)
            else
                read(TRT(JT), "(//2F8.0)") NXTTR2(JT), TTRNX(JT)
            end if

            if (TRIB_CONST(JT)) then
                TRC(JT) = NPT;                 NPT = NPT + 1

                if (WAIT_FOR_TRIB_INPUT(JT)) then ! Wait for input from tributary WQ file  !SR 11/28/19
                    FULL_FILE_NAME = TRIM(ADJUSTL(TR_FILEDIR(JT))) // "\" // TRIM(ADJUSTL(CTRFN(JT))) !SR 11/28/19
                    inquire(FILE=TRIM(FULL_FILE_NAME), EXIST=INPUT_FILE_EXISTS) ! Check that file exists                 !SR 11/28/19
                    do while (.not. INPUT_FILE_EXISTS) !SR 11/28/19
                        write(*, "(3A)") "Input file ", TRIM(FULL_FILE_NAME), " does not exist (yet).  Waiting..." !SR 11/28/19
                        write(9911, "(3A)") "Input file ", TRIM(FULL_FILE_NAME), " does not exist (yet).  Waiting..." !SR 11/28/19
                        call SLEEPQQ(WAIT_TIME*1000) !SR 11/28/19
                        inquire(FILE=TRIM(FULL_FILE_NAME), EXIST=INPUT_FILE_EXISTS) !SR 11/28/19
                    end do !SR 11/28/19
                    RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory       !SR 11/28/19
                    write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                    if (.not. RESULT1) then
                        call PRINT_ERROR_AND_STOP("CTR", JT)
                    end if ! Problem with copy; write msg and stop  !SR 11/28/19

                    LAST_JDAY = GET_LAST_JDAY(CTRFN(JT)) ! Find last JDAY in input file           !SR 11/28/19
                    do while (LAST_JDAY <= TMEND - 0.5 .and. JDAY > LAST_JDAY - TIME_BUFFER) ! Not enough data in input file          !SR 11/28/19
                        write(9911, "(A,I0,3(A,F0.4))") "WAIT: Input CTR", JT, " DAY= ", LAST_JDAY, " JDAY= ", JDAY, " TMEND= ", TMEND !SR 11/28/19
                        call SLEEPQQ(WAIT_TIME*1000) !SR 11/28/19
                        RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory       !SR 11/28/19
                        write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                        if (.not. RESULT1) then
                            call PRINT_ERROR_AND_STOP("CTR", JT)
                        end if ! Problem with copy; write msg and stop  !SR 11/28/19
                        LAST_JDAY = GET_LAST_JDAY(CTRFN(JT)) !SR 11/28/19
                    end do !SR 11/28/19
                end if !SR 11/28/19

                open(TRC(JT), FILE=CTRFN(JT), STATUS="OLD")
                read(TRC(JT), "(A1)") INFORMAT
                if (INFORMAT == "$") then
                    TRCF(JT) = .true.
                end if
                if (TRCF(JT)) then
                    read(TRC(JT), "(/)")
                    read(TRC(JT), *) NXCTR2(JT), (CTRNX(TRCN(JAC, JT), JT), JAC = 1, NACTR(JT))
                else
                    read(TRC(JT), "(//1000F8.0)") NXCTR2(JT), (CTRNX(TRCN(JAC, JT), JT), JAC = 1, NACTR(JT))
                end if

!     OPEN (TRQ(JT),FILE=QTRFN(JT),STATUS='OLD')
!     OPEN (TRT(JT),FILE=TTRFN(JT),STATUS='OLD')
!     
!      READ( TRQ(JT),'(A1)')INFORMAT
!      IF(INFORMAT=='$')TRQF(JT)=.TRUE.
!           IF(TRQF(JT))THEN
!           READ (TRQ(JT),'(/)')
!           READ (TRQ(JT),*) NXQTR2(JT),QTRNX(JT)
!           ELSE
!           READ (TRQ(JT),'(//2F8.0)') NXQTR2(JT),QTRNX(JT)
!           ENDIF
!     READ( TRT(JT),'(A1)')INFORMAT
!      IF(INFORMAT=='$')TRTF(JT)=.TRUE.
!           IF(TRTF(JT))THEN
!           READ (TRT(JT),'(/)')
!           READ (TRT(JT),*) NXTTR2(JT),TTRNX(JT)
!           ELSE
!           READ (TRT(JT),'(//2F8.0)') NXTTR2(JT),TTRNX(JT)
!           ENDIF
!     
!!     READ (TRQ(JT),'(///2F8.0)') NXQTR2(JT),QTRNX(JT)
!!     READ (TRT(JT),'(///2F8.0)') NXTTR2(JT),TTRNX(JT)
!     IF (TRIB_CONST(JT)) THEN
!       TRC(JT) = NPT; NPT = NPT+1
!       OPEN (TRC(JT),FILE=CTRFN(JT),STATUS='OLD')
!       READ( TRC(JT),'(A1)')INFORMAT
!       IF(INFORMAT=='$')TRCF(JT)=.TRUE.
!           IF(TRCF(JT))THEN
!           READ (TRC(JT),'(/)')
!           READ (TRC(JT),*) NXCTR2(JT),(CTRNX(TRCN(JAC,JT),JT),JAC=1,NACTR(JT))
!           ELSE
!           READ (TRC(JT),'(//1000F8.0)') NXCTR2(JT),(CTRNX(TRCN(JAC,JT),JT),JAC=1,NACTR(JT))
!           ENDIF
            end if
        end do
        QTR(1:NTR) = QTRNX(1:NTR)
        QTRO(1:NTR) = QTRNX(1:NTR)
        TTR(1:NTR) = TTRNX(1:NTR)
        TTRO(1:NTR) = TTRNX(1:NTR)
        CTR(:, 1:NTR) = CTRNX(:, 1:NTR)
        CTRO(:, 1:NTR) = CTRNX(:, 1:NTR)
!
! systdg - time series input 
        if (SYSTDG) then
            if (DOBND .or. N2BND .or. DGPBND) then
                do JT = 1, NTR
                    if (DOBND) then
                        DO_SAT(JT) = EXP(7.7117 - 1.31403*LOG(TTR(JT) + 45.93))*PALT(ITR(JT))
                        CTR(NDO, JT) = DO_SAT(JT)*CTR(NDO, JT)
                    end if
                    if (N2BND) then
                        EA = DEXP(2.3026D0*(7.5D0*TDEW(WBSEG(ITR(JT)))/(TDEW(WBSEG(ITR(JT))) + 237.3D0) + 0.6609D0))*0.001316
                        N2_SAT(JT) = 1.5568D06*0.79*(PALT(ITR(JT)) - EA)*(1.8816D-5 - 4.116D-7*TTR(JT) + 4.6D-9*TTR(JT)*TTR(JT))
                        CTR(NN2, JT) = N2_SAT(JT)*CTR(NN2, JT)
                    end if
!
                    if (DGPBND) then
                        CTR(NDGP, JT) = PALT(ITR(JT))*CTR(NDGP, JT)
                    end if
                end do
            end if
        end if

!
! systdg - add for time series input 
        do JT = 1, NTR
            if (TRQF(JT)) then
                read(TRQ(JT), *) NXQTR1(JT), QTRNX(JT) ! cb 5/22/14
            else
                read(TRQ(JT), "(2F8.0)") NXQTR1(JT), QTRNX(JT)
            end if
            if (TRTF(JT)) then
                read(TRT(JT), *) NXTTR1(JT), TTRNX(JT) ! cb 5/22/14
            else
                read(TRT(JT), "(2F8.0)") NXTTR1(JT), TTRNX(JT)
            end if
            if (TRIB_CONST(JT)) then
                if (TRCF(JT)) then
                    read(TRC(JT), *) NXCTR1(JT), (CTRNX(TRCN(JAC, JT), JT), JAC = 1, NACTR(JT))
                else
                    read(TRC(JT), "(1000F8.0)") NXCTR1(JT), (CTRNX(TRCN(JAC, JT), JT), JAC = 1, NACTR(JT))
                end if
            end if
        end do
    end if
!
! systdg - time series input TWE
    if (SYSTDG) then
        if (TWETSC == "      ON") then
            TWEFNNO = NPT;             NPT = NPT + 1
            open(TWEFNNO, FILE=TWEFN, STATUS="OLD")
            read(TWEFNNO, "(A1)") INFORMAT
            if (INFORMAT == "$") then
                TWEF = .true.
            end if
            if (TWEF) then
                read(TWEFNNO, "(/)")
                read(TWEFNNO, *) NXTWE2, TWE_TSNX
            else
                read(TWEFNNO, "(//2F8.3)") NXTWE2, TWE_TSNX
            end if
            TWE_TS = TWE_TSNX
            TWE_TSO = TWE_TSNX
            if (TWEF) then
                read(TWEFNNO, *) NXTWE1, TWE_TSNX
            else
                read(TWEFNNO, "(2F8.3)") NXTWE1, TWE_TSNX
            end if
        end if
    end if
!
! systdg - time series input TWE 
    do JW = 1, NWB
        do JB = BS(JW), BE(JW)
            if (UP_FLOW(JB)) then
                if (.not. INTERNAL_FLOW(JB) .and. .not. DAM_INFLOW(JB)) then !TC 08/03/04 RA 1/13/06
                    INQ(JB) = NPT;                     NPT = NPT + 1
                    INFT(JB) = NPT;                     NPT = NPT + 1

                    if (WAIT_FOR_BRANCH_INPUT(JB)) then ! Wait for input, branch flow file    !SR 11/28/19
                        FULL_FILE_NAME = TRIM(ADJUSTL(BR_FILEDIR(JB))) // "\" // TRIM(ADJUSTL(QINFN(JB))) !SR 11/28/19
                        inquire(FILE=TRIM(FULL_FILE_NAME), EXIST=INPUT_FILE_EXISTS) ! Check that file exists              !SR 11/28/19
                        do while (.not. INPUT_FILE_EXISTS) !SR 11/28/19
                            write(*, "(3A)") "Input file ", TRIM(FULL_FILE_NAME), " does not exist (yet).  Waiting..." !SR 11/28/19
                            write(9911, "(3A)") "Input file ", TRIM(FULL_FILE_NAME), " does not exist (yet).  Waiting..." !SR 11/28/19
                            call SLEEPQQ(WAIT_TIME*1000) !SR 11/28/19
                            inquire(FILE=TRIM(FULL_FILE_NAME), EXIST=INPUT_FILE_EXISTS) !SR 11/28/19
                        end do !SR 11/28/19
                        RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory    !SR 11/28/19
                        write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                        if (.not. RESULT1) then
                            call PRINT_ERROR_AND_STOP("QIN", JB)
                        end if ! Problem with copy; write msg, stop  !SR 11/28/19

                        LAST_JDAY = GET_LAST_JDAY(QINFN(JB)) ! Find last JDAY in input file        !SR 11/28/19
                        do while (LAST_JDAY <= TMEND - 0.5 .and. JDAY > LAST_JDAY - TIME_BUFFER) ! Not enough data in input file       !SR 11/28/19
                            write(9911, "(A,I0,3(A,F0.4))") "WAIT: Input QIN", JB, " DAY= ", LAST_JDAY, " JDAY= ", JDAY, " TMEND= ", TMEND !SR 11/28/19
                            call SLEEPQQ(WAIT_TIME*1000) !SR 11/28/19
                            RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory    !SR 11/28/19
                            write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                            if (.not. RESULT1) then
                                call PRINT_ERROR_AND_STOP("QIN", JB)
                            end if ! Problem with copy; write msg, stop  !SR 11/28/19
                            LAST_JDAY = GET_LAST_JDAY(QINFN(JB)) !SR 11/28/19
                        end do !SR 11/28/19
                    end if !SR 11/28/19
                    open(INQ(JB), FILE=QINFN(JB), STATUS="OLD")
                    read(INQ(JB), "(A1)") INFORMAT
                    if (INFORMAT == "$") then
                        INQF(JB) = .true.
                    end if
                    if (INQF(JB)) then
                        read(INQ(JB), "(/)")
                        read(INQ(JB), *) NXQIN2(JB), QINNX(JB)
                    else
                        read(INQ(JB), "(//2F8.0)") NXQIN2(JB), QINNX(JB)
                    end if

                    if (WAIT_FOR_BRANCH_INPUT(JB)) then ! Wait for input, branch temp file    !SR 11/28/19
                        FULL_FILE_NAME = TRIM(ADJUSTL(BR_FILEDIR(JB))) // "\" // TRIM(ADJUSTL(TINFN(JB))) !SR 11/28/19
                        inquire(FILE=TRIM(FULL_FILE_NAME), EXIST=INPUT_FILE_EXISTS) ! Check that file exists              !SR 11/28/19
                        do while (.not. INPUT_FILE_EXISTS) !SR 11/28/19
                            write(*, "(3A)") "Input file ", TRIM(FULL_FILE_NAME), " does not exist (yet).  Waiting..." !SR 11/28/19
                            write(9911, "(3A)") "Input file ", TRIM(FULL_FILE_NAME), " does not exist (yet).  Waiting..." !SR 11/28/19
                            call SLEEPQQ(WAIT_TIME*1000) !SR 11/28/19
                            inquire(FILE=TRIM(FULL_FILE_NAME), EXIST=INPUT_FILE_EXISTS) !SR 11/28/19
                        end do !SR 11/28/19
                        RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory    !SR 11/28/19
                        write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                        if (.not. RESULT1) then
                            call PRINT_ERROR_AND_STOP("TIN", JB)
                        end if ! Problem with copy; write msg, stop  !SR 11/28/19

                        LAST_JDAY = GET_LAST_JDAY(TINFN(JB)) ! Find last JDAY in input file        !SR 11/28/19
                        do while (LAST_JDAY <= TMEND - 0.5 .and. JDAY > LAST_JDAY - TIME_BUFFER) ! Not enough data in input file       !SR 11/28/19
                            write(9911, "(A,I0,3(A,F0.4))") "WAIT: Input TIN", JB, " DAY= ", LAST_JDAY, " JDAY= ", JDAY, " TMEND= ", TMEND !SR 11/28/19
                            call SLEEPQQ(WAIT_TIME*1000) !SR 11/28/19
                            RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory    !SR 11/28/19
                            write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                            if (.not. RESULT1) then
                                call PRINT_ERROR_AND_STOP("TIN", JB)
                            end if ! Problem with copy; write msg, stop  !SR 11/28/19
                            LAST_JDAY = GET_LAST_JDAY(TINFN(JB)) !SR 11/28/19
                        end do !SR 11/28/19
                    end if !SR 11/28/19
                    open(INFT(JB), FILE=TINFN(JB), STATUS="OLD")
                    read(INFT(JB), "(A1)") INFORMAT
                    if (INFORMAT == "$") then
                        INTF(JB) = .true.
                    end if
                    if (INTF(JB)) then
                        read(INFT(JB), "(/)")
                        read(INFT(JB), *) NXTIN2(JB), TINNX(JB)
                    else
                        read(INFT(JB), "(//2F8.0)") NXTIN2(JB), TINNX(JB)
                    end if

                    if (INFLOW_CONST(JB)) then
                        INC(JB) = NPT;                         NPT = NPT + 1

                        if (WAIT_FOR_BRANCH_INPUT(JB)) then ! Wait for input, branch WQ file     !SR 11/28/19
                            FULL_FILE_NAME = TRIM(ADJUSTL(BR_FILEDIR(JB))) // "\" // TRIM(ADJUSTL(CINFN(JB))) !SR 11/28/19
                            inquire(FILE=TRIM(FULL_FILE_NAME), EXIST=INPUT_FILE_EXISTS) ! Check that file exists             !SR 11/28/19
                            do while (.not. INPUT_FILE_EXISTS) !SR 11/28/19
                                write(*, "(3A)") "Input file ", TRIM(FULL_FILE_NAME), " does not exist (yet).  Waiting..." !SR 11/28/19
                                write(9911, "(3A)") "Input file ", TRIM(FULL_FILE_NAME), " does not exist (yet).  Waiting..." !SR 11/28/19
                                call SLEEPQQ(WAIT_TIME*1000) !SR 11/28/19
                                inquire(FILE=TRIM(FULL_FILE_NAME), EXIST=INPUT_FILE_EXISTS) !SR 11/28/19
                            end do !SR 11/28/19
                            RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory   !SR 11/28/19
                            write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                            if (.not. RESULT1) then
                                call PRINT_ERROR_AND_STOP("CIN", JB)
                            end if ! Problem with copy; write msg, stop !SR 11/28/19

                            LAST_JDAY = GET_LAST_JDAY(CINFN(JB)) ! Find last JDAY in input file       !SR 11/28/19
                            do while (LAST_JDAY <= TMEND - 0.5 .and. JDAY > LAST_JDAY - TIME_BUFFER) ! Not enough data in input file      !SR 11/28/19
                                write(9911, "(A,I0,3(A,F0.4))") "WAIT: Input CIN", JB, " DAY= ", LAST_JDAY, " JDAY= ", JDAY, " TMEND= ", TMEND !SR 11/28/19
                                call SLEEPQQ(WAIT_TIME*1000) !SR 11/28/19
                                RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory   !SR 11/28/19
                                write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                                if (.not. RESULT1) then
                                    call PRINT_ERROR_AND_STOP("CIN", JB)
                                end if ! Problem with copy; write msg, stop !SR 11/28/19
                                LAST_JDAY = GET_LAST_JDAY(CINFN(JB)) !SR 11/28/19
                            end do !SR 11/28/19
                        end if !SR 11/28/19
                        open(INC(JB), FILE=CINFN(JB), STATUS="OLD")
!DO JW=1,NWB
!  DO JB=BS(JW),BE(JW)
!    IF (UP_FLOW(JB)) THEN
!      IF (.NOT. INTERNAL_FLOW(JB) .AND. .NOT. DAM_INFLOW(JB)) THEN                                                  !TC 08/03/04 RA 1/13/06
!        INQ(JB)  = NPT; NPT = NPT+1
!        INFT(JB) = NPT; NPT = NPT+1
!        OPEN (INQ(JB) ,FILE=QINFN(JB),STATUS='OLD')
!        OPEN (INFT(JB),FILE=TINFN(JB),STATUS='OLD')
!     READ( INQ(JB),'(A1)')INFORMAT
!     IF(INFORMAT=='$')INQF(JB)=.TRUE.
!          IF(INQF(JB))THEN
!          READ (INQ(JB),'(/)')
!          READ (INQ(JB),*) NXQIN2(JB),QINNX(JB)
!          ELSE
!          READ (INQ(JB), '(//2F8.0)') NXQIN2(JB),QINNX(JB)
!          ENDIF
!     READ( INFT(JB),'(A1)')INFORMAT
!     IF(INFORMAT=='$')INTF(JB)=.TRUE.
!          IF(INTF(JB))THEN
!          READ (INFT(JB),'(/)')
!          READ (INFT(JB),*) NXTIN2(JB),TINNX(JB)
!          ELSE
!          READ (INFT(JB),'(//2F8.0)') NXTIN2(JB),TINNX(JB)
!          ENDIF      
!        
!   !     READ (INQ(JB), '(///2F8.0)') NXQIN2(JB),QINNX(JB)
!   !     READ (INFT(JB),'(///2F8.0)') NXTIN2(JB),TINNX(JB)
!        IF (INFLOW_CONST(JB)) THEN
!          INC(JB) = NPT; NPT = NPT+1
!          OPEN (INC(JB),FILE=CINFN(JB),STATUS='OLD')
                        read(INC(JB), "(A1)") INFORMAT
                        if (INFORMAT == "$") then
                            INCF(JB) = .true.
                        end if
                        if (INCF(JB)) then
                            read(INC(JB), "(/)")
                            read(INC(JB), *) NXCIN2(JB), (CINNX(INCN(JC, JB), JB), JC = 1, NACIN(JB))
                        else
                            read(INC(JB), "(//1000F8.0)") NXCIN2(JB), (CINNX(INCN(JC, JB), JB), JC = 1, NACIN(JB))
                        end if
                    end if
                end if
                QIN(JB) = QINNX(JB)
                QIND(JB) = QINNX(JB)
                QINO(JB) = QINNX(JB)
                TIN(JB) = TINNX(JB)
                TIND(JB) = TINNX(JB)
                TINO(JB) = TINNX(JB)
                CIN(:, JB) = CINNX(:, JB)
                CIND(:, JB) = CINNX(:, JB)
                CINO(:, JB) = CINNX(:, JB)
!
! systdg - time series input 
                if (SYSTDG) then
                    if (DOBND) then
                        DO_SATJ(JB) = EXP(7.7117 - 1.31403*LOG(TIN(JB) + 45.93))*PALT(CUS(JB))
                        CIN(NDO, JB) = DO_SATJ(JB)*CIN(NDO, JB)
                        CIND(NDO, JB) = DO_SATJ(JB)*CIND(NDO, JB)
                    end if
                    if (N2BND) then
                        EA = DEXP(2.3026D0*(7.5D0*TDEW(WBSEG(CUS(JB)))/(TDEW(WBSEG(CUS(JB))) + 237.3D0) + 0.6609D0))*0.001316
                        N2_SATJ(JB) = 1.5568D06*0.79*(PALT(CUS(JB)) - EA)*(1.8816D-5 - 4.116D-7*TIN(JB) + 4.6D-9*TIN(JB)*TIN(JB))
                        CIN(NN2, JB) = N2_SATJ(JB)*CIN(NN2, JB)
                        CIND(NN2, JB) = N2_SATJ(JB)*CIND(NN2, JB)
                    end if
!
                    if (DGPBND) then
                        CIND(NDGP, JB) = PALT(CUS(JB))*CIND(NDGP, JB)
                    end if
                end if
!
! systdg - time series input                                              
                if (.not. INTERNAL_FLOW(JB) .and. .not. DAM_INFLOW(JB)) then !TC 08/03/04  RA 1/13/06
                    if (INQF(JB)) then
                        read(INQ(JB), *) NXQIN1(JB), QINNX(JB)
                    else
                        read(INQ(JB), "(2F8.0)") NXQIN1(JB), QINNX(JB)
                    end if
                    if (INTF(JB)) then
                        read(INFT(JB), *) NXTIN1(JB), TINNX(JB)
                    else
                        read(INFT(JB), "(2F8.0)") NXTIN1(JB), TINNX(JB)
                    end if

                    if (INFLOW_CONST(JB)) then
                        if (INCF(JB)) then
                            read(INC(JB), *) NXCIN1(JB), (CINNX(INCN(JC, JB), JB), JC = 1, NACIN(JB))
                        else
                            read(INC(JB), "(1000F8.0)") NXCIN1(JB), (CINNX(INCN(JC, JB), JB), JC = 1, NACIN(JB))
                        end if
                    end if
                end if
            end if
            if (DN_FLOW(JB)) then
                if (NSTR(JB) > 0) then
                    OTQ(JB) = NPT;                     NPT = NPT + 1
                    open(OTQ(JB), FILE=QOTFN(JB), STATUS="OLD")
                    read(OTQ(JB), "(A1)") INFORMAT
                    if (INFORMAT == "$") then
                        OTQF(JB) = .true.
                    end if
                    if (OTQF(JB)) then
                        read(OTQ(JB), "(/)")
                        read(OTQ(JB), *) NXQOT2(JB), (QSTRNX(JS, JB), JS = 1, NSTR(JB))
                        QSTR(:, JB) = QSTRNX(:, JB)
                        QSTRO(:, JB) = QSTRNX(:, JB)
                        read(OTQ(JB), *) NXQOT1(JB), (QSTRNX(JS, JB), JS = 1, NSTR(JB))
                    else
                        read(OTQ(JB), "(//10F8.0:/(8X,9F8.0))") NXQOT2(JB), (QSTRNX(JS, JB), JS = 1, NSTR(JB))
                        QSTR(:, JB) = QSTRNX(:, JB)
                        QSTRO(:, JB) = QSTRNX(:, JB)
                        read(OTQ(JB), "(10F8.0:/(8X,9F8.0))") NXQOT1(JB), (QSTRNX(JS, JB), JS = 1, NSTR(JB))
                    end if

                    if (DYNSTRUC(JB) == "      ON") then
                        ODYNS(JB) = NPT;                         NPT = NPT + 1
                        write(SEGNUM, "(I0)") JB
                        SEGNUM = ADJUSTL(SEGNUM)
                        L = LEN_TRIM(SEGNUM)
                        open(ODYNS(JB), FILE="dynselev" // SEGNUM(1:L) // ".npt", STATUS="OLD")
                        read(ODYNS(JB), "(A1)") INFORMAT
                        if (INFORMAT == "$") then
                            DYNEF(JB) = .true.
                        end if
                        read(ODYNS(JB), *) NJS
                        do J = 1, NJS
                            read(ODYNS(JB), *) JJS(J)
                        end do
                        if (DYNEF(JB)) then
                            read(ODYNS(JB), *)
                            read(ODYNS(JB), *) NXDYNS(JB), (ESTR(JJS(J), JB), J = 1, NJS)
                            read(ODYNS(JB), *) NXDYNS(JB), (NXESTRT(JJS(J), JB), J = 1, NJS)
                        else
                            read(ODYNS(JB), "(/10F8.0:/(8X,9F8.0))") NXDYNS(JB), (ESTR(JJS(J), JB), J = 1, NJS)
                            read(ODYNS(JB), "(10F8.0:/(8X,9F8.0))") NXDYNS(JB), (NXESTRT(JJS(J), JB), J = 1, NJS)
                        end if
                    end if
                end if
            end if
            if (PRECIPITATION(JW)) then
                PRE(JB) = NPT;                 NPT = NPT + 1
                PRT(JB) = NPT;                 NPT = NPT + 1
                open(PRE(JB), FILE=PREFN(JB), STATUS="OLD")
                open(PRT(JB), FILE=TPRFN(JB), STATUS="OLD")
                read(PRE(JB), "(A1)") INFORMAT
                if (INFORMAT == "$") then
                    PRQF(JB) = .true.
                end if
                if (PRQF(JB)) then
                    read(PRE(JB), "(/)")
                    read(PRE(JB), *) NXPR2(JB), PRNX(JB)
                else
                    read(PRE(JB), "(//2F8.0)") NXPR2(JB), PRNX(JB)
                end if
                read(PRT(JB), "(A1)") INFORMAT
                if (INFORMAT == "$") then
                    PRTF(JB) = .true.
                end if
                if (PRTF(JB)) then
                    read(PRT(JB), "(/)")
                    read(PRT(JB), *) NXTPR2(JB), TPRNX(JB)
                else
                    read(PRT(JB), "(//2F8.0)") NXTPR2(JB), TPRNX(JB)
                end if

!READ (PRE(JB),'(///2F8.0)') NXPR2(JB), PRNX(JB)
!READ (PRT(JB),'(///2F8.0)') NXTPR2(JB),TPRNX(JB)
                if (PRECIP_CONST(JB)) then
                    PRC(JB) = NPT;                     NPT = NPT + 1
                    open(PRC(JB), FILE=CPRFN(JB), STATUS="OLD")
                    read(PRC(JB), "(A1)") INFORMAT
                    if (INFORMAT == "$") then
                        PRCF(JB) = .true.
                    end if
                    if (PRCF(JB)) then
                        read(PRC(JB), "(/)")
                        read(PRC(JB), *) NXCPR2(JB), (CPRNX(PRCN(JAC, JB), JB), JAC = 1, NACPR(JB))
                    else
                        read(PRC(JB), "(//1000F8.0)") NXCPR2(JB), (CPRNX(PRCN(JAC, JB), JB), JAC = 1, NACPR(JB))
                    end if
                end if
                PR(JB) = PRNX(JB)
                TPR(JB) = TPRNX(JB)
                CPR(:, JB) = CPRNX(:, JB)
!
! systdg - time series input 
                if (SYSTDG) then
                    if (DOBND) then
                        DO_SATP(JB) = EXP(7.7117 - 1.31403*LOG(TPR(JB) + 45.93))*PALT(CUS(JB))
                        CPR(NDO, JB) = DO_SATP(JB)*CPR(NDO, JB)
                    end if
                    if (N2BND) then
                        EA = DEXP(2.3026D0*(7.5D0*TDEW(WBSEG(CUS(JB)))/(TDEW(WBSEG(CUS(JB))) + 237.3D0) + 0.6609D0))*0.001316
                        N2_SATP(JB) = 1.5568D06*0.79*(PALT(CUS(JB)) - EA)*(1.8816D-5 - 4.116D-7*TPR(JB) + 4.6D-9*TPR(JB)*TPR(JB))
                        CPR(NN2, JB) = N2_SATP(JB)*CPR(NN2, JB)
                    end if
!
                    if (DGPBND) then
                        CPR(NDGP, JB) = PALT(CUS(JB))*CPR(NDGP, JB)
                    end if
                end if
!
! systdg - time series input 
                if (PRQF(JB)) then
                    read(PRE(JB), *) NXPR1(JB), PRNX(JB)
                else
                    read(PRE(JB), "(2F8.0)") NXPR1(JB), PRNX(JB)
                end if

                if (PRTF(JB)) then
                    read(PRT(JB), *) NXTPR1(JB), TPRNX(JB)
                else
                    read(PRT(JB), "(2F8.0)") NXTPR1(JB), TPRNX(JB)
                end if
!READ (PRE(JB),'(2F8.0)') NXPR1(JB), PRNX(JB)
!READ (PRT(JB),'(2F8.0)') NXTPR1(JB),TPRNX(JB)
                if (PRECIP_CONST(JB)) then
                    if (PRCF(JB)) then
                        read(PRC(JB), *) NXCPR1(JB), (CPRNX(PRCN(JAC, JB), JB), JAC = 1, NACPR(JB))
                    else
                        read(PRC(JB), "(1000F8.0)") NXCPR1(JB), (CPRNX(PRCN(JAC, JB), JB), JAC = 1, NACPR(JB))
                    end if
                end if
            end if
            if (DIST_TRIBS(JB)) then
                DTQ(JB) = NPT;                 NPT = NPT + 1
                DTT(JB) = NPT;                 NPT = NPT + 1
                open(DTQ(JB), FILE=QDTFN(JB), STATUS="OLD")
                open(DTT(JB), FILE=TDTFN(JB), STATUS="OLD")

                read(DTQ(JB), "(A1)") INFORMAT
                if (INFORMAT == "$") then
                    DTQF(JB) = .true.
                end if
                if (DTQF(JB)) then
                    read(DTQ(JB), "(/)")
                    read(DTQ(JB), *) NXQDT2(JB), QDTRNX(JB)
                else
                    read(DTQ(JB), "(//2F8.0)") NXQDT2(JB), QDTRNX(JB)
                end if
                read(DTT(JB), "(A1)") INFORMAT
                if (INFORMAT == "$") then
                    DTTF(JB) = .true.
                end if
                if (DTTF(JB)) then
                    read(DTT(JB), "(/)")
                    read(DTT(JB), *) NXTDT2(JB), TDTRNX(JB)
                else
                    read(DTT(JB), "(//2F8.0)") NXTDT2(JB), TDTRNX(JB)
                end if


!      READ (DTQ(JB),'(///2F8.0)') NXQDT2(JB),QDTRNX(JB)
!      READ (DTT(JB),'(///2F8.0)') NXTDT2(JB),TDTRNX(JB)
                if (DTRIB_CONST(JB)) then
                    DTC(JB) = NPT;                     NPT = NPT + 1
                    open(DTC(JB), FILE=CDTFN(JB), STATUS="OLD")
                    read(DTC(JB), "(A1)") INFORMAT
                    if (INFORMAT == "$") then
                        DTCF(JB) = .true.
                    end if
                    if (DTCF(JB)) then
                        read(DTC(JB), "(/)")
                        read(DTC(JB), *) NXCDT2(JB), (CDTRNX(DTCN(JAC, JB), JB), JAC = 1, NACDT(JB))
                    else
                        read(DTC(JB), "(//1000F8.0)") NXCDT2(JB), (CDTRNX(DTCN(JAC, JB), JB), JAC = 1, NACDT(JB))
                    end if
                end if
                QDTR(JB) = QDTRNX(JB)
                QDTRO(JB) = QDTRNX(JB)
                TDTR(JB) = TDTRNX(JB)
                TDTRO(JB) = TDTRNX(JB)
                CDTR(:, JB) = CDTRNX(:, JB)
                CDTRO(:, JB) = CDTRNX(:, JB)
!
! systdg - time series input 
                if (SYSTDG) then
                    if (DOBND) then
                        DO_SATD(JB) = EXP(7.7117 - 1.31403*LOG(TDTR(JB) + 45.93))*PALT(CUS(JB))
                        CDTR(NDO, JB) = DO_SATD(JB)*CDTR(NDO, JB)
                    end if
                    if (N2BND) then
                        EA = DEXP(2.3026D0*(7.5D0*TDEW(WBSEG(CUS(JB)))/(TDEW(WBSEG(CUS(JB))) + 237.3D0) + 0.6609D0))*0.001316
                        N2_SATD(JB) = 1.5568D06*0.79*(PALT(CUS(JB)) - EA)*(1.8816D-5 - 4.116D-7*TDTR(JB) + 4.6D-9*TDTR(JB)*TDTR(JB))
                        CDTR(NN2, JB) = N2_SATD(JB)*CDTR(NN2, JB)
                    end if
!
                    if (DGPBND) then
                        CDTR(NDGP, JB) = PALT(CUS(JB))*CDTR(NDGP, JB)
                    end if
                end if
!
! systdg - time series input                                         
                if (DTQF(JB)) then
                    read(DTQ(JB), *) NXQDT1(JB), QDTRNX(JB)
                else
                    read(DTQ(JB), "(2F8.0)") NXQDT1(JB), QDTRNX(JB)
                end if

                if (DTTF(JB)) then
                    read(DTT(JB), *) NXTDT1(JB), TDTRNX(JB)
                else
                    read(DTT(JB), "(2F8.0)") NXTDT1(JB), TDTRNX(JB)
                end if

                if (DTRIB_CONST(JB)) then
                    if (DTCF(JB)) then
                        read(DTC(JB), *) NXCDT1(JB), (CDTRNX(DTCN(JAC, JB), JB), JAC = 1, NACDT(JB))
                    else
                        read(DTC(JB), "(1000F8.0)") NXCDT1(JB), (CDTRNX(DTCN(JAC, JB), JB), JAC = 1, NACDT(JB))
                    end if
                end if
            end if
            if (UH_EXTERNAL(JB)) then
                UHE(JB) = NPT;                 NPT = NPT + 1
                UHT(JB) = NPT;                 NPT = NPT + 1
                open(UHE(JB), FILE=EUHFN(JB), STATUS="OLD")
                open(UHT(JB), FILE=TUHFN(JB), STATUS="OLD")

                read(UHE(JB), "(A1)") INFORMAT
                if (INFORMAT == "$") then
                    EUHF(JB) = 1
                end if

                read(UHT(JB), "(A2)") INFORMAT2
                if (INFORMAT2 == "$T") then
                    TUHF(JB) = 1
                else
                    if (INFORMAT2(1:1) == "$") then
                        TUHF(JB) = 2
                    end if
                end if

                if (EUHF(JB) > 0) then
                    read(UHE(JB), "(/)")
                    read(UHE(JB), *) NXEUH2(JB), ELUHNX(JB)
                else
                    read(UHE(JB), "(//2F8.0)") NXEUH2(JB), ELUHNX(JB)
                end if

                if (TUHF(JB) == 1) then
                    read(UHT(JB), "(/)")
                    read(UHT(JB), *) NXTUH2(JB), XX(1)
                    TUHNX(2:KB(US(JB)), JB) = XX(1)
                else
                    if (TUHF(JB) == 2) then
                        read(UHT(JB), "(/)")
                        read(UHT(JB), *) NXTUH2(JB), (TUHNX(K, JB), K = 2, KB(US(JB)))
                    else
                        read(UHT(JB), "(//10F8.0:/(8X,9F8.0))") NXTUH2(JB), (TUHNX(K, JB), K = 2, KB(US(JB)))
                    end if
                end if

! READ (UHE(JB),'(///2F8.0)')              NXEUH2(JB), ELUHNX(JB)
! READ (UHT(JB),'(///10F8.0:/(8X,9F8.0))') NXTUH2(JB),(TUHNX(K,JB),K=2,KB(US(JB)))
                if (CONSTITUENTS) then
                    UHC(JB) = NPT;                     NPT = NPT + 1
                    open(UHC(JB), FILE=CUHFN(JB), STATUS="OLD")

                    read(UHC(JB), "(A2)") INFORMAT2
                    if (INFORMAT2 == "$T") then
                        CUHF(JB) = 1
                    else
                        if (INFORMAT2(1:1) == "$") then
                            CUHF(JB) = 2
                        end if
                    end if

                    read(UHC(JB), "(/)")

!  READ (UHC(JB),'(//)')
                    if (CUHF(JB) == 1) then
                        read(UHC(JB), *) NXCUH2(JB), (XX(CN(JAC)), JAC = 1, NAC)
                        do JAC = 1, NAC
                            CUHNX(2:KB(US(JB)), CN(JAC), JB) = XX(CN(JAC))
                        end do

                    else


                        do JAC = 1, NAC
!IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (UHC(JB),'(10F8.0:/(8X,9F8.0))') NXCUH2(JB),(CUHNX(K,CN(JAC),JB),     &
!                                                  K=2,KB(US(JB)))

                            if (CUHF(JB) == 2) then
                                read(UHC(JB), *) NXCUH2(JB), (CUHNX(K, CN(JAC), JB), K = 2, KB(US(JB)))
                            else
                                read(UHC(JB), "(10F8.0:/(8X,9F8.0))") NXCUH2(JB), (CUHNX(K, CN(JAC), JB), K = 2, KB(US(JB)))
                            end if


                        end do
                    end if

                end if
                ELUH(JB) = ELUHNX(JB)
                ELUHO(JB) = ELUHNX(JB)
                TUH(:, JB) = TUHNX(:, JB)
                TUHO(:, JB) = TUHNX(:, JB)
                CUH(:, :, JB) = CUHNX(:, :, JB)
                CUHO(:, :, JB) = CUHNX(:, :, JB)

                if (EUHF(JB) > 0) then
                    read(UHE(JB), *) NXEUH1(JB), ELUHNX(JB)
                else
                    read(UHE(JB), "(2F8.0)") NXEUH1(JB), ELUHNX(JB)
                end if

!READ (UHE(JB),'(2F8.0)')              NXEUH1(JB), ELUHNX(JB)

                if (TUHF(JB) == 1) then
                    read(UHT(JB), *) NXTUH1(JB), XX(1)
                    TUHNX(2:KB(US(JB)), JB) = XX(1)
                else
                    if (TUHF(JB) == 2) then
                        read(UHT(JB), *) NXTUH1(JB), (TUHNX(K, JB), K = 2, KB(US(JB)))
                    else
                        read(UHT(JB), "(10F8.0:/(8X,9F8.0))") NXTUH1(JB), (TUHNX(K, JB), K = 2, KB(US(JB)))
                    end if
                end if

! READ (UHT(JB),'(10F8.0:/(8X,9F8.0))') NXTUH1(JB),(TUHNX(K,JB),K=2,KB(US(JB)))
                if (CONSTITUENTS) then

                    if (CUHF(JB) == 1) then
                        read(UHC(JB), *) NXCUH1(JB), (XX(CN(JAC)), JAC = 1, NAC)
                        do JAC = 1, NAC
                            CUHNX(2:KB(US(JB)), CN(JAC), JB) = XX(CN(JAC))
                        end do

                    else


                        do JAC = 1, NAC

                            if (CUHF(JB) == 2) then
                                read(UHC(JB), *) NXCUH1(JB), (CUHNX(K, CN(JAC), JB), K = 2, KB(US(JB)))
                            else
                                read(UHC(JB), "(10F8.0:/(8X,9F8.0))") NXCUH1(JB), (CUHNX(K, CN(JAC), JB), K = 2, KB(US(JB)))
                            end if

                        end do
                    end if


!DO JAC=1,NAC
!  IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (UHC(JB),'(10F8.0:/(8X,9F8.0))') NXCUH1(JB),(CUHNX(K,CN(JAC),JB),     &
!                                                    K=2,KB(US(JB)))
!END DO
                end if
            end if
            if (DH_EXTERNAL(JB)) then
                DHE(JB) = NPT;                 NPT = NPT + 1
                DHT(JB) = NPT;                 NPT = NPT + 1
                open(DHE(JB), FILE=EDHFN(JB), STATUS="OLD")
                open(DHT(JB), FILE=TDHFN(JB), STATUS="OLD")

                read(DHE(JB), "(A1)") INFORMAT
                if (INFORMAT == "$") then
                    EDHF(JB) = 1
                end if

                read(DHT(JB), "(A2)") INFORMAT2
                if (INFORMAT2 == "$T") then
                    TDHF(JB) = 1
                else
                    if (INFORMAT2(1:1) == "$") then
                        TDHF(JB) = 2
                    end if
                end if

                if (EDHF(JB) > 0) then
                    read(DHE(JB), "(/)")
                    read(DHE(JB), *) NXEDH2(JB), ELDHNX(JB)
                else
                    read(DHE(JB), "(//2F8.0)") NXEDH2(JB), ELDHNX(JB)
                end if

                if (TDHF(JB) == 1) then
                    read(DHT(JB), "(/)")
                    read(DHT(JB), *) NXTDH2(JB), XX(1)
                    TDHNX(2:KB(DS(JB)), JB) = XX(1)
                else
                    if (TDHF(JB) == 2) then
                        read(DHT(JB), "(/)")
                        read(DHT(JB), *) NXTDH2(JB), (TDHNX(K, JB), K = 2, KB(DS(JB)))
                    else
                        read(DHT(JB), "(//10F8.0:/(8X,9F8.0))") NXTDH2(JB), (TDHNX(K, JB), K = 2, KB(DS(JB)))
                    end if
                end if

!READ (DHE(JB),'(///10F8.0)')             NXEDH2(JB),ELDHNX(JB)
!READ (DHT(JB),'(///10F8.0:/(8X,9F8.0))') NXTDH2(JB),(TDHNX(K,JB),K=2,KB(DS(JB)))
                if (CONSTITUENTS) then
                    DHC(JB) = NPT;                     NPT = NPT + 1
                    open(DHC(JB), FILE=CDHFN(JB), STATUS="OLD")

                    read(DHC(JB), "(A2)") INFORMAT2
                    if (INFORMAT2 == "$T") then
                        CDHF(JB) = 1
                    else
                        if (INFORMAT2(1:1) == "$") then
                            CDHF(JB) = 2
                        end if
                    end if

                    read(DHC(JB), "(/)")

                    if (CDHF(JB) == 1) then
                        read(DHC(JB), *) NXCDH2(JB), (XX(CN(JAC)), JAC = 1, NAC)
                        do JAC = 1, NAC
                            CDHNX(2:KB(DS(JB)), CN(JAC), JB) = XX(CN(JAC))
                        end do
                    else
                        do JAC = 1, NAC
                            if (CDHF(JB) == 2) then
                                read(DHC(JB), *) NXCDH2(JB), (CDHNX(K, CN(JAC), JB), K = 2, KB(DS(JB)))
                            else
                                read(DHC(JB), "(10F8.0:/(8X,9F8.0))") NXCDH2(JB), (CDHNX(K, CN(JAC), JB), K = 2, KB(DS(JB)))
                            end if
                        end do
                    end if

!READ (DHC(JB),'(//)')
!DO JAC=1,NAC
!  IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (DHC(JB),'(10F8.0:/(8X,9F8.0))') NXCDH2(JB),(CDHNX(K,CN(JAC),JB),     &
!                                                    K=2,KB(DS(JB)))
!END DO
                end if
                ELDH(JB) = ELDHNX(JB)
                ELDHO(JB) = ELDHNX(JB)
                TDH(:, JB) = TDHNX(:, JB)
                TDHO(:, JB) = TDHNX(:, JB)
                CDH(:, :, JB) = CDHNX(:, :, JB)
                CDHO(:, :, JB) = CDHNX(:, :, JB)

                if (EDHF(JB) > 0) then
                    read(DHE(JB), *) NXEDH1(JB), ELDHNX(JB)
                else
                    read(DHE(JB), "(2F8.0)") NXEDH1(JB), ELDHNX(JB)
                end if

                if (TDHF(JB) == 1) then
                    read(DHT(JB), *) NXTDH1(JB), XX(1)
                    TDHNX(2:KB(DS(JB)), JB) = XX(1)
                else
                    if (TDHF(JB) == 2) then
                        read(DHT(JB), *) NXTDH1(JB), (TDHNX(K, JB), K = 2, KB(DS(JB)))
                    else
                        read(DHT(JB), "(10F8.0:/(8X,9F8.0))") NXTDH1(JB), (TDHNX(K, JB), K = 2, KB(DS(JB)))
                    end if
                end if

!READ (DHE(JB),'(10F8.0)')             NXEDH1(JB),ELDHNX(JB)
!READ (DHT(JB),'(10F8.0:/(8X,9F8.0))') NXTDH1(JB),(TDHNX(K,JB),K=2,KB(DS(JB)))
                if (CONSTITUENTS) then

                    if (CDHF(JB) == 1) then
                        read(DHC(JB), *) NXCDH1(JB), (XX(CN(JAC)), JAC = 1, NAC)
                        do JAC = 1, NAC
                            CDHNX(2:KB(DS(JB)), CN(JAC), JB) = XX(CN(JAC))
                        end do
                    else
                        do JAC = 1, NAC
                            if (CDHF(JB) == 2) then
                                read(DHC(JB), *) NXCDH1(JB), (CDHNX(K, CN(JAC), JB), K = 2, KB(DS(JB)))
                            else
                                read(DHC(JB), "(10F8.0:/(8X,9F8.0))") NXCDH1(JB), (CDHNX(K, CN(JAC), JB), K = 2, KB(DS(JB)))
                            end if
                        end do
                    end if

!DO JAC=1,NAC
!  IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (DHC(JB),'(10F8.0:/(8X,9F8.0))') NXCDH1(JB),(CDHNX(K,CN(JAC),JB),     &
!                                                    K=2,KB(DS(JB)))
!END DO
                end if
            end if
        end do ! JB LOOP
        if (CONSTITUENTS) then
            if (ATM_DEPOSITION(JW)) then

                ATMDEP(JW) = NPT;                 NPT = NPT + 1
                open(ATMDEP(JW), FILE=ATMDEPFN(JW), STATUS="OLD")

                read(ATMDEP(JW), "(A1)") INFORMAT
                if (INFORMAT == "$") then
                    ATMDEPCSV = .true.
                end if
                read(ATMDEP(JW), "(/)")

                if (ATMDEPCSV) then
                    read(ATMDEP(JW), *) NXATMD2(JW), (ATM_DEP_LOADINGNX(ATMDCN(JAC, JW), JW), JAC = 1, NACATD(JW))
                else
                    read(ATMDEP(JW), "(100F8.0)") NXATMD2(JW), (ATM_DEP_LOADINGNX(ATMDCN(JAC, JW), JW), JAC = 1, NACATD(JW))
                end if
                ATM_DEP_LOADING = ATM_DEP_LOADINGNX
                ATM_DEP_LOADING0 = ATM_DEP_LOADINGNX
                if (ATMDEPCSV) then
                    read(ATMDEP(JW), *) NXATMD(JW), (ATM_DEP_LOADINGNX(ATMDCN(JAC, JW), JW), JAC = 1, NACATD(JW))
                else
                    read(ATMDEP(JW), "(100F8.0)") NXATMD(JW), (ATM_DEP_LOADINGNX(ATMDCN(JAC, JW), JW), JAC = 1, NACATD(JW))
                end if
            end if
        end if
    end do ! JW LOOP
    if (GATES) then
        GTQ = NPT;         NPT = NPT + 1
        open(GTQ, FILE=QGTFN, STATUS="OLD")
        read(GTQ, "(A1)") INFORMAT
        if (INFORMAT == "$") then
            GATEF = .true.
        end if
! Added code is to allow for computing flow based on BGT elevation - this can be a target water level in the reservoir - but allowing that flow to be removed from a different elevation
!READ(GTQ,*)                           ! SW 2/25/11
        if (GATEF) then
            read(GTQ, "(A8)") GT2CHAR

            if (GT2CHAR == "EGT2ELEV") then
                rewind(GTQ)
                read(GTQ, *)
                read(GTQ, *) GT2CHAR, (EGT2(JG), JG = 1, NGT)
            end if

            read(GTQ, *)
            read(GTQ, *) NXQGT2, (BGTNX(JG), JG = 1, NGT)
            where (DYNGTC == "     ZGT")
                EGT = BGTNX
                egto = bgtnx
                BGT = 1.0
                G1GT = 1.0
                G2GT = 1.0
            elsewhere
                BGT = BGTNX
                bgto = bgtnx
            end where
            read(GTQ, *) NXQGT, (BGTNX(JG), JG = 1, NGT)

        else

            read(GTQ, "(A8)") GT2CHAR

            if (GT2CHAR == "EGT2ELEV") then
                rewind(GTQ)
                read(GTQ, *)
                read(GTQ, "(8X,1000F8.0)") (EGT2(JG), JG = 1, NGT)
            end if

            read(GTQ, *)
            read(GTQ, "(1000F8.0)") NXQGT2, (BGTNX(JG), JG = 1, NGT)
!    READ (GTQ,'(///1000F8.0)') NXQGT2,(BGTNX(JG),JG=1,NGT)                                       
            where (DYNGTC == "     ZGT")
                EGT = BGTNX
                egto = bgtnx
                BGT = 1.0
                G1GT = 1.0
                G2GT = 1.0
            elsewhere
                BGT = BGTNX
                bgto = bgtnx
            end where
            read(GTQ, "(1000F8.0)") NXQGT, (BGTNX(JG), JG = 1, NGT)
        end if

    end if
    if (PIPES) then ! SW 5/5/10
        IOPENPIPE = 0
        do j = 1, npi
            if (DYNPIPE(j) == "      ON") then
                iopenpipe = 1
                exit
            end if
        end do
        if (iopenpipe == 1) then
            PIPED = NPT;             NPT = NPT + 1
            open(PIPED, FILE="dynpipe.npt", STATUS="OLD")
            read(PIPED, "(///1000F8.0)") NXQPT, (BPNX(J), J = 1, NPI)
            BP = BPNX
            read(PIPED, "(1000F8.0)") NXQPT, (BPNX(J), J = 1, NPI)
        end if
    end if
    if (PUMPS) then
        do J = 1, NPU
            if (DYNPUMP(j) == "      ON") then
                write(SEGNUM, "(I0)") J
                SEGNUM = ADJUSTL(SEGNUM)
                L = LEN_TRIM(SEGNUM)
                PUMPD(J) = NPT;                 NPT = NPT + 1
                open(PUMPD(J), FILE="dynpump" // SEGNUM(1:L) // ".npt", STATUS="OLD")
                read(PUMPD(J), "(A1)") INFORMAT
                if (INFORMAT == "$") then
                    DYNPUMPF(J) = .true.
                end if
                if (DYNPUMPF(J)) then
                    read(PUMPD(J), "(/)")
                    read(PUMPD(J), *) NXPUMP(J), EPU2(J), EONPU2(J), EOFFPU2(J), QPU2(J)
                    EPU(J) = EPU2(J)
                    EONPU(J) = EONPU2(J)
                    EOFFPU(J) = EOFFPU2(J)
                    QPU(J) = QPU2(J)
                    read(PUMPD(J), *) NXPUMP(J), EPU2(J), EONPU2(J), EOFFPU2(J), QPU2(J)
                else
                    read(PUMPD(J), "(//1000F8.0)") NXPUMP(J), EPU2(J), EONPU2(J), EOFFPU2(J), QPU2(J)
                    EPU(J) = EPU2(J)
                    EONPU(J) = EONPU2(J)
                    EOFFPU(J) = EOFFPU2(J)
                    QPU(J) = QPU2(J)
                    read(PUMPD(J), "(1000F8.0)") NXPUMP(J), EPU2(J), EONPU2(J), EOFFPU2(J), QPU2(J)
                end if
!READ (PUMPD(J),'(///1000F8.0)') NXPUMP(J),EPU2(J),EONPU2(J),EOFFPU2(J),QPU2(J)
!  EPU(J)=EPU2(J)
!  EONPU(J)=EONPU2(J)
!  EOFFPU(J)=EOFFPU2(J)
!  QPU(J)=QPU2(J)
!READ (PUMPD(J),'(1000F8.0)') NXPUMP(J),EPU2(J),EONPU2(J),EOFFPU2(J),QPU2(J)
            end if
        end do
    end if

    NOPEN = NPT - 1
    DYNAMIC_SHADE = SHADEI < 0
    NUNIT = NPT
    return

!***********************************************************************************************************************************
!**                                                  R E A D  I N P U T  D A T A                                                  **
!***********************************************************************************************************************************

    entry READ_INPUT_DATA(NXTVD)
    NXTVD = 1.0E10

! Meteorological data

    do while (JDAY >= NXWSC)
        WSC = WSCNX
        if (WSHF) then
            read(WSH, *) NXWSC, (WSCNX(I), I = 1, IMX)
        else
            read(WSH, "(10F8.0:/(8X,9F8.0))") NXWSC, (WSCNX(I), I = 1, IMX)
        end if
    end do
    do JW = 1, NWB
        do while (JDAY >= NXMET1(JW))
            TDEW(JW) = TDEWNX(JW)
            TDEWO(JW) = TDEWNX(JW)
            WIND(JW) = WINDNX(JW)
            WINDO(JW) = WINDNX(JW)
            PHI(JW) = PHINX(JW)
            PHIO(JW) = PHINX(JW)
!IF (PHISET > 0) PHI(JW)  = PHISET
!IF (PHISET > 0) PHIO(JW) = PHISET
            TAIR(JW) = TAIRNX(JW)
            TAIRO(JW) = TAIRNX(JW)
            CLOUD(JW) = CLOUDNX(JW)
            CLOUDO(JW) = CLOUDNX(JW)
            NXMET2(JW) = NXMET1(JW)
!
            if (SYSTDG) then
! systdg - time series input PALT
                PALT_JW(JW) = PALT_JWNX(JW)
                PALT_JWO(JW) = PALT_JWNX(JW)
                if (PALT_JW(JW) <= 0.0) then
                    PALT_JW(JW) = 760.0
                end if
                do I = US(BS(JW)) - 1, DS(BE(JW)) + 1
                    PALT(I) = PALT_JW(JW)/760.0*(1.0 - ELWS_INI(I)/1000.0/44.3)**5.25 ! systdg - time series input PALT
                end do
!PALT(:) = PALT_JW(JW)/760.0*(1.0-ELWS_INI(:)/1000.0/44.3)**5.25
! systdg - add time series input PALT

!
                if (READ_RADIATION(JW)) then
                    SRON(JW) = SRONX(JW)
                    SROO(JW) = SRON(JW)
                    if (METF(JW)) then
                        read(MET(JW), *) NXMET1(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW), SRONX(JW), PALT_JWNX(JW) ! systdg - time series input PALT_JW
                    else
                        read(MET(JW), "(8F8.0)") NXMET1(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW), SRONX(JW), PALT_JWNX(JW) ! systdg - time series input PALT_JW
                    end if
                    SRONX(JW) = SRONX(JW)*REFL
                else
                    if (METF(JW)) then
                        read(MET(JW), *) NXMET1(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW), PALT_JWNX(JW) ! systdg - time series input PALT_JW
                    else
                        read(MET(JW), "(7F8.0)") NXMET1(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW), PALT_JWNX(JW) ! systdg - time series input PALT_JW
                    end if
                end if
            else
                if (READ_RADIATION(JW)) then
                    SRON(JW) = SRONX(JW)
                    SROO(JW) = SRON(JW)
                    if (METF(JW)) then
                        read(MET(JW), *) NXMET1(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW), SRONX(JW)
                    else
                        read(MET(JW), "(8F8.0)") NXMET1(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW), SRONX(JW)
                    end if
                    SRONX(JW) = SRONX(JW)*REFL
                else
                    if (METF(JW)) then
                        read(MET(JW), *) NXMET1(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW)
                    else
                        read(MET(JW), "(7F8.0)") NXMET1(JW), TAIRNX(JW), TDEWNX(JW), WINDNX(JW), PHINX(JW), CLOUDNX(JW)
                    end if
                end if

            end if

        end do
        NXTVD = MIN(NXTVD, NXMET1(JW))
        if (READ_EXTINCTION(JW)) then
            do while (JDAY >= NXEXT1(JW))
                EXH2O(JW) = EXTNX(JW)
                EXTO(JW) = EXTNX(JW)
                NXEXT2(JW) = NXEXT1(JW)
                if (EXTF(JW)) then
                    read(EXT(JW), *) NXEXT1(JW), EXTNX(JW)
                else
                    read(EXT(JW), "(2F8.0)") NXEXT1(JW), EXTNX(JW)
                end if
            end do
        end if
!DO I=CUS(BS(JW)),DS(BE(JW))   ! SW CODE FIX 5-21-15
!  WIND2(I) = WIND(JW)*WSC(I)*DLOG(2.0D0/Z0(JW))/DLOG(WINDH(JW)/Z0(JW))    ! old value  z0 == 0.003
!END DO
    end do

! Withdrawals

    if (NWD > 0) then
        do while (JDAY >= NXQWD1)
            NXQWD2 = NXQWD1
            do JWD = 1, NWD
                QWD(JWD) = QWDNX(JWD)
                QWDO(JWD) = QWDNX(JWD)
            end do
            if (WDQF) then
                read(WDQ, *) NXQWD1, (QWDNX(JWD), JWD = 1, NWD)
            else
                read(WDQ, "(10F8.0:/(8X,9F8.0))") NXQWD1, (QWDNX(JWD), JWD = 1, NWD)
            end if
        end do
        NXTVD = MIN(NXTVD, NXQWD1)
    end if

! Spillways DO gas
    do N = 1, NSP ! SW 1/18/2022
        if (GASSPC(N) == "      ON" .and. EQSP(N) == 4 .and. BGASSP(N) == 1) then
            do while (JDAY >= NXSPDO(N))
                AGASSP(N) = AGASSPNX(N)
                read(FGASSP(N), *) NXSPDO(N), AGASSPNX(N)
            end do
        end if
    end do

! gates DO gas

    do N = 1, NGT ! SW 1/18/2022
        if (GASGTC(N) == "      ON" .and. EQGT(N) == 4 .and. BGASGT(N) == 1.0) then
            do while (JDAY >= NXGTDO(N))
                AGASGT(N) = AGASGTNX(N)
                read(FGASGT(N), *) NXGTDO(N), AGASGTNX(N)
            end do
        end if
    end do

! Tributaries

    if (TRIBUTARIES) then
        do JT = 1, NTR

!**** Inflow
            do while (JDAY >= NXQTR1(JT))
                QTR(JT) = QTRNX(JT)
                QTRO(JT) = QTRNX(JT)
                NXQTR2(JT) = NXQTR1(JT)

                if (TRQF(JT)) then
                    read(TRQ(JT), *, END=8710) NXQTR1(JT), QTRNX(JT) !SR 11/28/19
                else
                    read(TRQ(JT), "(2F8.0)", END=8710) NXQTR1(JT), QTRNX(JT) !SR 11/28/19
                end if
                go to 8712 ! Isolate error instructions            !SR 11/28/19
                8710 if (EOF(TRQ(JT))) then ! End of file, but more data needed     !SR 11/28/19
                    if (WAIT_FOR_TRIB_INPUT(JT)) then ! Additional data might be available    !SR 11/28/19
                        close(TRQ(JT)) ! Must close the file to get new copy   !SR 11/28/19
                        FULL_FILE_NAME = TRIM(ADJUSTL(TR_FILEDIR(JT))) // "\" // TRIM(ADJUSTL(QTRFN(JT))) !SR 11/28/19
                        RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory      !SR 11/28/19
                        write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                        if (.not. RESULT1) then
                            call PRINT_ERROR_AND_STOP("QTR", JT)
                        end if ! Problem with copy; write msg, stop    !SR 11/28/19

                        LAST_JDAY = GET_LAST_JDAY(QTRFN(JT)) ! Find last JDAY in input file          !SR 11/28/19
                        ITER = 0 !SR 11/28/19
                        do while (LAST_JDAY <= TMEND - 0.5 .and. JDAY > LAST_JDAY - TIME_BUFFER .or. LAST_JDAY > TMEND - 0.5 .and. LAST_JDAY < TMEND .and. ITER < 3) ! Not enough data in input file  !SR 11/28/19
                            ITER = ITER + 1 !SR 11/28/19
                            write(9911, "(A,I0,3(A,F0.4))") "WAIT: Input QTR", JT, " DAY= ", LAST_JDAY, " JDAY= ", JDAY, " TMEND= ", TMEND !SR 11/28/19
                            call SLEEPQQ(WAIT_TIME*1000) !SR 11/28/19
                            RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory      !SR 11/28/19
                            write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                            if (.not. RESULT1) then
                                call PRINT_ERROR_AND_STOP("QTR", JT)
                            end if ! Problem with copy; write msg, stop    !SR 11/28/19
                            LAST_JDAY = GET_LAST_JDAY(QTRFN(JT)) !SR 11/28/19
                        end do !SR 11/28/19
                        if (LAST_JDAY > TMEND - 0.5 .and. ITER >= 3 .and. LAST_JDAY <= NXQTR2(JT)) then ! near TMEND, no new data   !SR 11/28/19
                            open(TRQ(JT), FILE=QTRFN(JT), STATUS="OLD", POSITION="APPEND") ! Open file and push pointer to end     !SR 11/28/19
                            NXQTR1(JT) = TMEND + 1.0 ! Push input date past TMEND            !SR 11/28/19
                            write(9911, "(A,I0,2(A,F0.4))") "MOVING ON: Input QTR", JT, " DAY= ", LAST_JDAY, " close to TMEND= ", TMEND !SR 11/28/19
                        else !SR 11/28/19
                            LAST_JDAY = NXQTR2(JT) ! Save the last date read               !SR 11/28/19
                            open(TRQ(JT), FILE=QTRFN(JT), STATUS="OLD") ! Open the newly copied file and read   !SR 11/28/19
                            read(TRQ(JT), "(A1)") INFORMAT !SR 11/28/19
                            TRQF(JT) = .false. !SR 11/28/19
                            if (INFORMAT == "$") then
                                TRQF(JT) = .true.
                            end if !SR 11/28/19
                            if (TRQF(JT)) then !SR 11/28/19
                                read(TRQ(JT), "(/)") !SR 11/28/19
                                read(TRQ(JT), *) NXQTR1(JT) ! Just read the date                    !SR 11/28/19
                                do while (LAST_JDAY > NXQTR1(JT) .and. .not. EOF(TRQ(JT))) ! Get file ptr to previous position     !SR 11/28/19
                                    read(TRQ(JT), *) NXQTR1(JT) ! Just read the date                    !SR 11/28/19
                                end do !SR 11/28/19
                                if (EOF(TRQ(JT))) then
                                    backspace(TRQ(JT))
                                end if ! A bit of insurance                    !SR 11/28/19
                                read(TRQ(JT), *) NXQTR1(JT), QTRNX(JT) ! Read new data point                   !SR 11/28/19
                            else !SR 11/28/19
                                read(TRQ(JT), "(//F8.0)") NXQTR1(JT) ! Just read the date                    !SR 11/28/19
                                do while (LAST_JDAY > NXQTR1(JT) .and. .not. EOF(TRQ(JT))) ! Get file ptr to previous position     !SR 11/28/19
                                    read(TRQ(JT), "(F8.0)") NXQTR1(JT) ! Just read the date                    !SR 11/28/19
                                end do !SR 11/28/19
                                if (EOF(TRQ(JT))) then
                                    backspace(TRQ(JT))
                                end if ! A bit of insurance                    !SR 11/28/19
                                read(TRQ(JT), "(2F8.0)") NXQTR1(JT), QTRNX(JT) ! Read new data point                   !SR 11/28/19
                            end if !SR 11/28/19
                        end if !SR 11/28/19
                    else ! Not waiting for input from this file. Stop run. File has no more data. !SR 11/28/19
                        write(W2ERR, "(2A/A,F0.4)") "ERROR-- End of tributary input file ", QTRFN(JT), "Simulation terminated at JDAY ", JDAY !SR 11/28/19
                        write(*, "(2A/A,F0.4)") "ERROR-- End of tributary input file ", QTRFN(JT), "Simulation terminated at JDAY ", JDAY !SR 11/28/19
                        stop !SR 11/28/19
                    end if !SR 11/28/19
                end if !SR 11/28/19
                8712 continue !SR 11/28/19
            end do
!    DO WHILE (JDAY >= NXQTR1(JT))
!      QTR(JT)    = QTRNX(JT)
!      QTRO(JT)   = QTRNX(JT)
!      NXQTR2(JT) = NXQTR1(JT)
!      
!        IF(TRQF(JT))THEN
!        READ (TRQ(JT),*) NXQTR1(JT),QTRNX(JT)
!        ELSE
!      READ (TRQ(JT),'(2F8.0)') NXQTR1(JT),QTRNX(JT)
!        ENDIF
!    END DO
            NXTVD = MIN(NXTVD, NXQTR1(JT))

!**** Inflow temperatures
            do while (JDAY >= NXTTR1(JT))
                TTR(JT) = TTRNX(JT)
                TTRO(JT) = TTRNX(JT)
                NXTTR2(JT) = NXTTR1(JT)
                if (TRTF(JT)) then
                    read(TRT(JT), *, END=8720) NXTTR1(JT), TTRNX(JT) !SR 11/28/19
                else
                    read(TRT(JT), "(2F8.0)", END=8720) NXTTR1(JT), TTRNX(JT) !SR 11/28/19
                end if
                go to 8722 ! Isolate error instructions            !SR 11/28/19
                8720 if (EOF(TRT(JT))) then ! End of file, but more data needed     !SR 11/28/19
                    if (WAIT_FOR_TRIB_INPUT(JT)) then ! Additional data might be available    !SR 11/28/19
                        close(TRT(JT)) ! Must close the file to get new copy   !SR 11/28/19
                        FULL_FILE_NAME = TRIM(ADJUSTL(TR_FILEDIR(JT))) // "\" // TRIM(ADJUSTL(TTRFN(JT))) !SR 11/28/19
                        RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory      !SR 11/28/19
                        write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                        if (.not. RESULT1) then
                            call PRINT_ERROR_AND_STOP("TTR", JT)
                        end if ! Problem with copy; write msg, stop    !SR 11/28/19

                        LAST_JDAY = GET_LAST_JDAY(TTRFN(JT)) ! Find last JDAY in input file          !SR 11/28/19
! ITER = 0                                                            ! Carry over iterations from QTR        !SR 11/28/19
                        do while (LAST_JDAY <= TMEND - 0.5 .and. JDAY > LAST_JDAY - TIME_BUFFER .or. LAST_JDAY > TMEND - 0.5 .and. LAST_JDAY < TMEND .and. ITER < 3) ! Not enough data in input file  !SR 11/28/19
                            ITER = ITER + 1 !SR 11/28/19
                            write(9911, "(A,I0,3(A,F0.4))") "WAIT: Input TTR", JT, " DAY= ", LAST_JDAY, " JDAY= ", JDAY, " TMEND= ", TMEND !SR 11/28/19
                            call SLEEPQQ(WAIT_TIME*1000) !SR 11/28/19
                            RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory      !SR 11/28/19
                            write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                            if (.not. RESULT1) then
                                call PRINT_ERROR_AND_STOP("TTR", JT)
                            end if ! Problem with copy; write msg, stop    !SR 11/28/19
                            LAST_JDAY = GET_LAST_JDAY(TTRFN(JT)) !SR 11/28/19
                        end do !SR 11/28/19
                        if (LAST_JDAY > TMEND - 0.5 .and. ITER >= 3 .and. LAST_JDAY <= NXTTR2(JT)) then ! near TMEND, no new data   !SR 11/28/19
                            open(TRT(JT), FILE=TTRFN(JT), STATUS="OLD", POSITION="APPEND") ! Open file and push pointer to end     !SR 11/28/19
                            NXTTR1(JT) = TMEND + 1.0 ! Push input date past TMEND            !SR 11/28/19
                            write(9911, "(A,I0,2(A,F0.4))") "MOVING ON: Input TTR", JT, " DAY= ", LAST_JDAY, " close to TMEND= ", TMEND !SR 11/28/19
                        else !SR 11/28/19
                            LAST_JDAY = NXTTR2(JT) ! Save the last date read               !SR 11/28/19
                            open(TRT(JT), FILE=TTRFN(JT), STATUS="OLD") ! Open the newly copied file and read   !SR 11/28/19
                            read(TRT(JT), "(A1)") INFORMAT !SR 11/28/19
                            TRTF(JT) = .false. !SR 11/28/19
                            if (INFORMAT == "$") then
                                TRTF(JT) = .true.
                            end if !SR 11/28/19
                            if (TRTF(JT)) then !SR 11/28/19
                                read(TRT(JT), "(/)") !SR 11/28/19
                                read(TRT(JT), *) NXTTR1(JT) ! Just read the date                    !SR 11/28/19
                                do while (LAST_JDAY > NXTTR1(JT) .and. .not. EOF(TRT(JT))) ! Get file ptr to previous position     !SR 11/28/19
                                    read(TRT(JT), *) NXTTR1(JT) ! Just read the date                    !SR 11/28/19
                                end do !SR 11/28/19
                                if (EOF(TRT(JT))) then
                                    backspace(TRT(JT))
                                end if ! A bit of insurance                    !SR 11/28/19
                                read(TRT(JT), *) NXTTR1(JT), TTRNX(JT) ! Read new data point                   !SR 11/28/19
                            else !SR 11/28/19
                                read(TRT(JT), "(//F8.0)") NXTTR1(JT) ! Just read the date                    !SR 11/28/19
                                do while (LAST_JDAY > NXTTR1(JT) .and. .not. EOF(TRT(JT))) ! Get file ptr to previous position     !SR 11/28/19
                                    read(TRT(JT), "(F8.0)") NXTTR1(JT) ! Just read the date                    !SR 11/28/19
                                end do !SR 11/28/19
                                if (EOF(TRT(JT))) then
                                    backspace(TRT(JT))
                                end if ! A bit of insurance                    !SR 11/28/19
                                read(TRT(JT), "(2F8.0)") NXTTR1(JT), TTRNX(JT) ! Read new data point                   !SR 11/28/19
                            end if !SR 11/28/19
                        end if !SR 11/28/19
                    else ! Not waiting for input from this file. Stop run. File has no more data. !SR 11/28/19
                        write(W2ERR, "(2A/A,F0.4)") "ERROR-- End of tributary input file ", TTRFN(JT), "Simulation terminated at JDAY ", JDAY !SR 11/28/19
                        write(*, "(2A/A,F0.4)") "ERROR-- End of tributary input file ", TTRFN(JT), "Simulation terminated at JDAY ", JDAY !SR 11/28/19
                        stop !SR 11/28/19
                    end if !SR 11/28/19
                end if !SR 11/28/19
                8722 continue !SR 11/28/19
            end do

!   IF (JDAY >= NXTTR1(JT)) THEN
!     DO WHILE (JDAY >= NXTTR1(JT))
!       TTR(JT)    = TTRNX(JT)
!       TTRO(JT)   = TTRNX(JT)
!       NXTTR2(JT) = NXTTR1(JT)
!       
!       IF(TRTF(JT))THEN
!       READ (TRT(JT),*) NXTTR1(JT),TTRNX(JT)
!       ELSE
!       READ (TRT(JT),'(2F8.0)') NXTTR1(JT),TTRNX(JT)
!       ENDIF
!       
!       
!!       READ (TRT(JT),'(2F8.0)') NXTTR1(JT),TTRNX(JT)
!     END DO
!   END IF
            NXTVD = MIN(NXTVD, NXTTR1(JT))

!**** Inflow constituent concentrations

            if (TRIB_CONST(JT)) then
                do while (JDAY >= NXCTR1(JT))
                    CTR(TRCN(1:NACTR(JT), JT), JT) = CTRNX(TRCN(1:NACTR(JT), JT), JT)
                    CTRO(TRCN(1:NACTR(JT), JT), JT) = CTRNX(TRCN(1:NACTR(JT), JT), JT)
!
! systdg - time series input
                    if (SYSTDG) then
                        if (DOBND) then
                            DO_SAT(JT) = EXP(7.7117 - 1.31403*LOG(TTR(JT) + 45.93))*PALT(ITR(JT))
                            CTR(NDO, JT) = CTRNX(NDO, JT)*DO_SAT(JT)
                        end if
                        if (N2BND) then
                            EA = DEXP(2.3026D0*(7.5D0*TDEW(WBSEG(ITR(JT)))/(TDEW(WBSEG(ITR(JT))) + 237.3D0) + 0.6609D0))*0.001316

                            N2_SAT(JT) = 1.5568D06*0.79*(PALT(ITR(JT)) - EA)*(1.8816D-5 - 4.116D-7*TTR(JT) + 4.6D-9*TTR(JT)*TTR(JT))
                            CTR(NN2, JT) = CTRNX(NN2, JT)*N2_SAT(JT)
                        end if
                        if (DGPBND) then
                            CTR(NDGP, JT) = PALT(ITR(JT))*CTRNX(NDGP, JT)
                        end if
                    end if

! systdg - time series input
!
                    NXCTR2(JT) = NXCTR1(JT)
                    if (TRCF(JT)) then
                        read(TRC(JT), *, END=8730) NXCTR1(JT), (CTRNX(TRCN(JAC, JT), JT), JAC = 1, NACTR(JT)) !SR 11/28/19
                    else
                        read(TRC(JT), "(1000F8.0)", END=8730) NXCTR1(JT), (CTRNX(TRCN(JAC, JT), JT), JAC = 1, NACTR(JT)) !SR 11/28/19
                    end if
                    go to 8732 ! Isolate error instructions          !SR 11/28/19
                    8730 if (EOF(TRC(JT))) then ! End of file, but more data needed   !SR 11/28/19
                        if (WAIT_FOR_TRIB_INPUT(JT)) then ! Additional data might be available  !SR 11/28/19
                            close(TRC(JT)) ! Must close the file to get new copy !SR 11/28/19
                            FULL_FILE_NAME = TRIM(ADJUSTL(TR_FILEDIR(JT))) // "\" // TRIM(ADJUSTL(CTRFN(JT))) !SR 11/28/19
                            RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory    !SR 11/28/19
                            write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                            if (.not. RESULT1) then
                                call PRINT_ERROR_AND_STOP("CTR", JT)
                            end if ! Problem with copy; write msg, stop  !SR 11/28/19

                            LAST_JDAY = GET_LAST_JDAY(CTRFN(JT)) ! Find last JDAY in input file        !SR 11/28/19
! ITER = 0                                                            ! Carry over iterations from QTR,TTR  !SR 11/28/19
                            do while (LAST_JDAY <= TMEND - 0.5 .and. JDAY > LAST_JDAY - TIME_BUFFER .or. LAST_JDAY > TMEND - 0.5 .and. LAST_JDAY < TMEND .and. ITER < 3) ! Not enough data in input file!SR 11/28/19
                                ITER = ITER + 1 !SR 11/28/19
                                write(9911, "(A,I0,3(A,F0.4))") "WAIT: Input CTR", JT, " DAY= ", LAST_JDAY, " JDAY= ", JDAY, " TMEND= ", TMEND !SR 11/28/19
                                call SLEEPQQ(WAIT_TIME*1000) !SR 11/28/19
                                RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory    !SR 11/28/19
                                write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                                if (.not. RESULT1) then
                                    call PRINT_ERROR_AND_STOP("CTR", JT)
                                end if ! Problem with copy; write msg, stop  !SR 11/28/19
                                LAST_JDAY = GET_LAST_JDAY(CTRFN(JT)) !SR 11/28/19
                            end do !SR 11/28/19
                            if (LAST_JDAY > TMEND - 0.5 .and. ITER >= 3 .and. LAST_JDAY <= NXCTR2(JT)) then ! near TMEND, no new data !SR 11/28/19
                                open(TRC(JT), FILE=CTRFN(JT), STATUS="OLD", POSITION="APPEND") ! Open file and push pointer to end   !SR 11/28/19
                                NXCTR1(JT) = TMEND + 1.0 ! Push input date past TMEND          !SR 11/28/19
                                write(9911, "(A,I0,2(A,F0.4))") "MOVING ON: Input CTR", JT, " DAY= ", LAST_JDAY, " close to TMEND= ", TMEND !SR 11/28/19
                            else !SR 11/28/19
                                LAST_JDAY = NXCTR2(JT) ! Save the last date read             !SR 11/28/19
                                open(TRC(JT), FILE=CTRFN(JT), STATUS="OLD") ! Open the newly copied file and read !SR 11/28/19
                                read(TRC(JT), "(A1)") INFORMAT !SR 11/28/19
                                TRCF(JT) = .false. !SR 11/28/19
                                if (INFORMAT == "$") then
                                    TRCF(JT) = .true.
                                end if !SR 11/28/19
                                if (TRCF(JT)) then !SR 11/28/19
                                    read(TRC(JT), "(/)") !SR 11/28/19
                                    read(TRC(JT), *) NXCTR1(JT) ! Just read the date                  !SR 11/28/19
                                    do while (LAST_JDAY > NXCTR1(JT) .and. .not. EOF(TRC(JT))) ! Get file ptr to previous position   !SR 11/28/19
                                        read(TRC(JT), *) NXCTR1(JT) ! Just read the date                  !SR 11/28/19
                                    end do !SR 11/28/19
                                    if (EOF(TRC(JT))) then
                                        backspace(TRC(JT))
                                    end if ! A bit of insurance                  !SR 11/28/19
                                    read(TRC(JT), *) NXCTR1(JT), (CTRNX(TRCN(JAC, JT), JT), JAC = 1, NACTR(JT)) ! Read new data point      !SR 11/28/19
                                else !SR 11/28/19
                                    read(TRC(JT), "(//F8.0)") NXCTR1(JT) ! Just read the date                  !SR 11/28/19
                                    do while (LAST_JDAY > NXCTR1(JT) .and. .not. EOF(TRC(JT))) ! Get file ptr to previous position   !SR 11/28/19
                                        read(TRC(JT), "(F8.0)") NXCTR1(JT) ! Just read the date                  !SR 11/28/19
                                    end do !SR 11/28/19
                                    if (EOF(TRC(JT))) then
                                        backspace(TRC(JT))
                                    end if ! A bit of insurance                  !SR 11/28/19
                                    read(TRC(JT), "(1000F8.0)") NXCTR1(JT), (CTRNX(TRCN(JAC, JT), JT), JAC = 1, NACTR(JT)) ! Read new data      !SR 11/28/19
                                end if !SR 11/28/19
                            end if !SR 11/28/19
                        else ! Not waiting for input from this file. Stop run. File has no more data. !SR 11/28/19
                            write(W2ERR, "(2A/A,F0.4)") "ERROR-- End of tributary input file ", CTRFN(JT), "Simulation terminated at JDAY ", JDAY !SR 11/28/19
                            write(*, "(2A/A,F0.4)") "ERROR-- End of tributary input file ", CTRFN(JT), "Simulation terminated at JDAY ", JDAY !SR 11/28/19
                            stop !SR 11/28/19
                        end if !SR 11/28/19
                    end if !SR 11/28/19
                    8732 continue !SR 11/28/19

!  READ (TRC(JT),*) NXCTR1(JT),(CTRNX(TRCN(JAC,JT),JT),JAC=1,NACTR(JT))
!  ELSE
!  READ (TRC(JT),'(1000F8.0)') NXCTR1(JT),(CTRNX(TRCN(JAC,JT),JT),JAC=1,NACTR(JT))
!  ENDIF
                end do
                NXTVD = MIN(NXTVD, NXCTR1(JT))
            end if
        end do
    end if
!
! systdg - time series input TWE
    if (SYSTDG) then
        if (TWETSC == "      ON") then
            do while (JDAY >= NXTWE1)
                TWE_TS = TWE_TSNX
                TWE_TSO = TWE_TSNX
                NXTWE2 = NXTWE1
                if (TWEF) then
                    read(TWEFNNO, *) NXTWE1, TWE_TSNX
                else
                    read(TWEFNNO, "(2F8.3)") NXTWE1, TWE_TSNX
                end if
            end do
            NXTVD = MIN(NXTVD, NXTWE1)
        end if
    end if
! systdg - time series input TWE

! Branch related inputs

    do JW = 1, NWB
        do JB = BS(JW), BE(JW)

!**** Inflow

            if (UP_FLOW(JB)) then
                if (.not. INTERNAL_FLOW(JB) .and. .not. DAM_INFLOW(JB)) then !TC 08/03/04 RA 1/13/06
                    do while (JDAY >= NXQIN1(JB))
                        QIND(JB) = QINNX(JB)
                        QINO(JB) = QINNX(JB)
                        NXQIN2(JB) = NXQIN1(JB)
                        if (INQF(JB)) then
                            read(INQ(JB), *, END=8810) NXQIN1(JB), QINNX(JB) !SR 11/28/19
                        else
                            read(INQ(JB), "(2F8.0)", END=8810) NXQIN1(JB), QINNX(JB) !SR 11/28/19
                        end if
                        go to 8812 ! Isolate error instructions         !SR 11/28/19
                        8810 if (EOF(INQ(JB))) then ! End of file, but more data needed  !SR 11/28/19
                            if (WAIT_FOR_BRANCH_INPUT(JB)) then ! Additional data might be available !SR 11/28/19
                                close(INQ(JB)) ! Must close file to get new copy    !SR 11/28/19
                                FULL_FILE_NAME = TRIM(ADJUSTL(BR_FILEDIR(JB))) // "\" // TRIM(ADJUSTL(QINFN(JB))) !SR 11/28/19
                                RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory   !SR 11/28/19
                                write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                                if (.not. RESULT1) then
                                    call PRINT_ERROR_AND_STOP("QIN", JB)
                                end if ! Problem with copy; write msg, stop !SR 11/28/19

                                LAST_JDAY = GET_LAST_JDAY(QINFN(JB)) ! Find last JDAY in input file       !SR 11/28/19
                                ITER = 0 !SR 11/28/19
                                do while (LAST_JDAY <= TMEND - 0.5 .and. JDAY > LAST_JDAY - TIME_BUFFER .or. LAST_JDAY > TMEND - 0.5 .and. LAST_JDAY < TMEND .and. ITER < 3) ! Not enough data in file    !SR 11/28/19
                                    ITER = ITER + 1 !SR 11/28/19
                                    write(9911, "(A,I0,3(A,F0.4))") "WAIT: Input QIN", JB, " DAY= ", LAST_JDAY, " JDAY= ", JDAY, " TMEND= ", TMEND !SR 11/28/19
                                    call SLEEPQQ(WAIT_TIME*1000) !SR 11/28/19
                                    RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory   !SR 11/28/19
                                    write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                                    if (.not. RESULT1) then
                                        call PRINT_ERROR_AND_STOP("QIN", JB)
                                    end if ! Problem with copy; write msg, stop !SR 11/28/19
                                    LAST_JDAY = GET_LAST_JDAY(QINFN(JB)) !SR 11/28/19
                                end do !SR 11/28/19
                                if (LAST_JDAY > TMEND - 0.5 .and. ITER >= 3 .and. LAST_JDAY <= NXQIN2(JB)) then ! near TMEND, no new data !SR 11/28/19
                                    open(INQ(JB), FILE=QINFN(JB), STATUS="OLD", POSITION="APPEND") ! Open file and push pointer to end  !SR 11/28/19
                                    NXQIN1(JB) = TMEND + 1.0 ! Push input date past TMEND         !SR 11/28/19
                                    write(9911, "(A,I0,2(A,F0.4))") "MOVING ON: Input QIN", JB, " DAY= ", LAST_JDAY, " close to TMEND= ", TMEND !SR 11/28/19
                                else !SR 11/28/19
                                    LAST_JDAY = NXQIN2(JB) ! Save the last date read            !SR 11/28/19
                                    open(INQ(JB), FILE=QINFN(JB), STATUS="OLD") ! Open newly copied file and read    !SR 11/28/19
                                    read(INQ(JB), "(A1)") INFORMAT !SR 11/28/19
                                    INQF(JB) = .false. !SR 11/28/19
                                    if (INFORMAT == "$") then
                                        INQF(JB) = .true.
                                    end if !SR 11/28/19
                                    if (INQF(JB)) then !SR 11/28/19
                                        read(INQ(JB), "(/)") !SR 11/28/19
                                        read(INQ(JB), *) NXQIN1(JB) ! Just read the date                 !SR 11/28/19
                                        do while (LAST_JDAY > NXQIN1(JB) .and. .not. EOF(INQ(JB))) ! Get file ptr to previous position  !SR 11/28/19
                                            read(INQ(JB), *) NXQIN1(JB) ! Just read the date                 !SR 11/28/19
                                        end do !SR 11/28/19
                                        if (EOF(INQ(JB))) then
                                            backspace(INQ(JB))
                                        end if ! A bit of insurance                 !SR 11/28/19
                                        read(INQ(JB), *) NXQIN1(JB), QINNX(JB) ! Read new data point                !SR 11/28/19
                                    else !SR 11/28/19
                                        read(INQ(JB), "(//F8.0)") NXQIN1(JB) ! Just read the date                 !SR 11/28/19
                                        do while (LAST_JDAY > NXQIN1(JB) .and. .not. EOF(INQ(JB))) ! Get file ptr to previous position  !SR 11/28/19
                                            read(INQ(JB), "(F8.0)") NXQIN1(JB) ! Just read the date                 !SR 11/28/19
                                        end do !SR 11/28/19
                                        if (EOF(INQ(JB))) then
                                            backspace(INQ(JB))
                                        end if ! A bit of insurance                 !SR 11/28/19
                                        read(INQ(JB), "(2F8.0)") NXQIN1(JB), QINNX(JB) ! Read new data point                !SR 11/28/19
                                    end if !SR 11/28/19
                                end if !SR 11/28/19
                            else ! Not waiting for input from this file. Stop run. File has no more data.!SR 11/28/19
                                write(W2ERR, "(2A/A,F0.4)") "ERROR-- End of branch input file ", QINFN(JB), "Simulation terminated at JDAY ", JDAY !SR 11/28/19
                                write(*, "(2A/A,F0.4)") "ERROR-- End of branch input file ", QINFN(JB), "Simulation terminated at JDAY ", JDAY !SR 11/28/19
                                stop !SR 11/28/19
                            end if !SR 11/28/19
                        end if !SR 11/28/19
                        8812 continue !SR 11/28/19
                    end do
!  IF(INQF(JB))THEN
!  READ (INQ(JB),*) NXQIN1(JB),QINNX(JB)
!  ELSE
!    READ (INQ(JB),'(2F8.0)') NXQIN1(JB),QINNX(JB)
!  ENDIF 
!  END DO
                    NXTVD = MIN(NXTVD, NXQIN1(JB))

!******** Inflow temperature

                    do while (JDAY >= NXTIN1(JB))
                        TIND(JB) = TINNX(JB)
                        TINO(JB) = TINNX(JB)
                        NXTIN2(JB) = NXTIN1(JB)
                        if (INTF(JB)) then
                            read(INFT(JB), *, END=8820) NXTIN1(JB), TINNX(JB) !SR 11/28/19
                        else
                            read(INFT(JB), "(2F8.0)", END=8820) NXTIN1(JB), TINNX(JB) !SR 11/28/19
                        end if
                        go to 8822 ! Isolate error instructions         !SR 11/28/19
                        8820 if (EOF(INFT(JB))) then ! End of file, but more data needed  !SR 11/28/19
                            if (WAIT_FOR_BRANCH_INPUT(JB)) then ! Additional data might be available !SR 11/28/19
                                close(INFT(JB)) ! Must close file to get new copy    !SR 11/28/19
                                FULL_FILE_NAME = TRIM(ADJUSTL(BR_FILEDIR(JB))) // "\" // TRIM(ADJUSTL(TINFN(JB))) !SR 11/28/19
                                RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory   !SR 11/28/19
                                write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                                if (.not. RESULT1) then
                                    call PRINT_ERROR_AND_STOP("TIN", JB)
                                end if ! Problem with copy; write msg, stop !SR 11/28/19

                                LAST_JDAY = GET_LAST_JDAY(TINFN(JB)) ! Find last JDAY in input file       !SR 11/28/19
! ITER = 0                                                           ! Carry over iterations from QIN     !SR 11/28/19
                                do while (LAST_JDAY <= TMEND - 0.5 .and. JDAY > LAST_JDAY - TIME_BUFFER .or. LAST_JDAY > TMEND - 0.5 .and. LAST_JDAY < TMEND .and. ITER < 3) ! Not enough data in file    !SR 11/28/19
                                    ITER = ITER + 1 !SR 11/28/19
                                    write(9911, "(A,I0,3(A,F0.4))") "WAIT: Input TIN", JB, " DAY= ", LAST_JDAY, " JDAY= ", JDAY, " TMEND= ", TMEND !SR 11/28/19
                                    call SLEEPQQ(WAIT_TIME*1000) !SR 11/28/19
                                    RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory   !SR 11/28/19
                                    write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                                    if (.not. RESULT1) then
                                        call PRINT_ERROR_AND_STOP("TIN", JB)
                                    end if ! Problem with copy; write msg, stop !SR 11/28/19
                                    LAST_JDAY = GET_LAST_JDAY(TINFN(JB)) !SR 11/28/19
                                end do !SR 11/28/19
                                if (LAST_JDAY > TMEND - 0.5 .and. ITER >= 3 .and. LAST_JDAY <= NXTIN2(JB)) then ! near TMEND, no new data !SR 11/28/19
                                    open(INFT(JB), FILE=TINFN(JB), STATUS="OLD", POSITION="APPEND") ! Open file and push pointer to end  !SR 11/28/19
                                    NXTIN1(JB) = TMEND + 1.0 ! Push input date past TMEND         !SR 11/28/19
                                    write(9911, "(A,I0,2(A,F0.4))") "MOVING ON: Input TIN", JB, " DAY= ", LAST_JDAY, " close to TMEND= ", TMEND !SR 11/28/19
                                else !SR 11/28/19
                                    LAST_JDAY = NXTIN2(JB) ! Save the last date read            !SR 11/28/19
                                    open(INFT(JB), FILE=TINFN(JB), STATUS="OLD") ! Open newly copied file and read    !SR 11/28/19
                                    read(INFT(JB), "(A1)") INFORMAT !SR 11/28/19
                                    INTF(JB) = .false. !SR 11/28/19
                                    if (INFORMAT == "$") then
                                        INTF(JB) = .true.
                                    end if !SR 11/28/19
                                    if (INTF(JB)) then !SR 11/28/19
                                        read(INFT(JB), "(/)") !SR 11/28/19
                                        read(INFT(JB), *) NXTIN1(JB) ! Just read the date                 !SR 11/28/19
                                        do while (LAST_JDAY > NXTIN1(JB) .and. .not. EOF(INFT(JB))) ! Get file ptr to previous position  !SR 11/28/19
                                            read(INFT(JB), *) NXTIN1(JB) ! Just read the date                 !SR 11/28/19
                                        end do !SR 11/28/19
                                        if (EOF(INFT(JB))) then
                                            backspace(INFT(JB))
                                        end if ! A bit of insurance                 !SR 11/28/19
                                        read(INFT(JB), *) NXTIN1(JB), TINNX(JB) ! Read new data point                !SR 11/28/19
                                    else !SR 11/28/19
                                        read(INFT(JB), "(//F8.0)") NXTIN1(JB) ! Just read the date                 !SR 11/28/19
                                        do while (LAST_JDAY > NXTIN1(JB) .and. .not. EOF(INFT(JB))) ! Get file ptr to previous position  !SR 11/28/19
                                            read(INFT(JB), "(F8.0)") NXTIN1(JB) ! Just read the date                 !SR 11/28/19
                                        end do !SR 11/28/19
                                        if (EOF(INFT(JB))) then
                                            backspace(INFT(JB))
                                        end if ! A bit of insurance                 !SR 11/28/19
                                        read(INFT(JB), "(2F8.0)") NXTIN1(JB), TINNX(JB) ! Read new data point                !SR 11/28/19
                                    end if !SR 11/28/19
                                end if !SR 11/28/19
                            else ! Not waiting for input from this file. Stop run. File has no more data.!SR 11/28/19
                                write(W2ERR, "(2A/A,F0.4)") "ERROR-- End of branch input file ", TINFN(JB), "Simulation terminated at JDAY ", JDAY !SR 11/28/19
                                write(*, "(2A/A,F0.4)") "ERROR-- End of branch input file ", TINFN(JB), "Simulation terminated at JDAY ", JDAY !SR 11/28/19
                                stop !SR 11/28/19
                            end if !SR 11/28/19
                        end if !SR 11/28/19
                        8822 continue !SR 11/28/19
                    end do

!     IF(INTF(JB))THEN
!     READ (INFT(JB),*) NXTIN1(JB),TINNX(JB)
!     ELSE
!       READ (INFT(JB),'(2F8.0)') NXTIN1(JB),TINNX(JB)
!     ENDIF 
!       
!!       READ (INFT(JB),'(2F8.0)') NXTIN1(JB),TINNX(JB)
!     END DO
                    NXTVD = MIN(NXTVD, NXTIN1(JB))

!******** Inflow constituent concentrations

                    if (INFLOW_CONST(JB)) then
                        do while (JDAY >= NXCIN1(JB))
                            CIND(INCN(1:NACIN(JB), JB), JB) = CINNX(INCN(1:NACIN(JB), JB), JB)
                            CINO(INCN(1:NACIN(JB), JB), JB) = CINNX(INCN(1:NACIN(JB), JB), JB)
!
! systdg - time series input
                            if (SYSTDG) then
                                if (DOBND) then
                                    DO_SATJ(JB) = EXP(7.7117 - 1.31403*LOG(TIND(JB) + 45.93))*PALT(CUS(JB))
                                    CIND(NDO, JB) = CIND(NDO, JB)*DO_SATJ(JB)
                                end if
                                if (N2BND) then
                                    EA = DEXP(2.3026D0*(7.5D0*TDEW(WBSEG(CUS(JB)))/(TDEW(WBSEG(CUS(JB))) + 237.3D0) + 0.6609D0))*0.001316

                                    N2_SATJ(JB) = 1.5568D06*0.79*(PALT(CUS(JB)) - EA)*(1.8816D-5 - 4.116D-7*TIND(JB) + 4.6D-9*TIND(JB)*TIND(JB))
                                    CIND(NN2, JB) = CIND(NN2, JB)*N2_SATJ(JB)
                                end if
                                if (DGPBND) then
                                    CIND(NDGP, JB) = CIND(NDGP, JB)*PALT(CUS(JB))
                                end if
                            end if
! systdg - time series input

                            NXCIN2(JB) = NXCIN1(JB)
                            if (INCF(JB)) then
                                read(INC(JB), *, END=8830) NXCIN1(JB), (CINNX(INCN(JAC, JB), JB), JAC = 1, NACIN(JB)) !SR 11/28/19
                            else
                                read(INC(JB), "(1000F8.0)", END=8830) NXCIN1(JB), (CINNX(INCN(JAC, JB), JB), JAC = 1, NACIN(JB)) !SR 11/28/19
                            end if
                            go to 8832 ! Isolate error instructions         !SR 11/28/19
                            8830 if (EOF(INC(JB))) then ! End of file, but more data needed  !SR 11/28/19
                                if (WAIT_FOR_BRANCH_INPUT(JB)) then ! Additional data might be available !SR 11/28/19
                                    close(INC(JB)) ! Must close file to get new copy    !SR 11/28/19
                                    FULL_FILE_NAME = TRIM(ADJUSTL(BR_FILEDIR(JB))) // "\" // TRIM(ADJUSTL(CINFN(JB))) !SR 11/28/19
                                    RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory   !SR 11/28/19
                                    write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                                    if (.not. RESULT1) then
                                        call PRINT_ERROR_AND_STOP("CIN", JB)
                                    end if ! Problem with copy; write msg, stop !SR 11/28/19

                                    LAST_JDAY = GET_LAST_JDAY(CINFN(JB)) ! Find last JDAY in input file       !SR 11/28/19
! ITER = 0                                                         ! Carry over iterations from QIN,TIN !SR 11/28/19
                                    do while (LAST_JDAY <= TMEND - 0.5 .and. JDAY > LAST_JDAY - TIME_BUFFER .or. LAST_JDAY > TMEND - 0.5 .and. LAST_JDAY < TMEND .and. ITER < 3) ! Not enough data in file  !SR 11/28/19
                                        ITER = ITER + 1 !SR 11/28/19
                                        write(9911, "(A,I0,3(A,F0.4))") "WAIT: Input CIN", JB, " DAY= ", LAST_JDAY, " JDAY= ", JDAY, " TMEND= ", TMEND !SR 11/28/19
                                        call SLEEPQQ(WAIT_TIME*1000) !SR 11/28/19
                                        RESULT1 = SYSTEMQQ("COPY " // TRIM(FULL_FILE_NAME)) ! Copy file into current directory   !SR 11/28/19
                                        write(9911, "(F12.4,2X,A,A)") JDAY, "COPY ", TRIM(FULL_FILE_NAME) !SR 11/28/19
                                        if (.not. RESULT1) then
                                            call PRINT_ERROR_AND_STOP("CIN", JB)
                                        end if ! Problem with copy; write msg, stop !SR 11/28/19
                                        LAST_JDAY = GET_LAST_JDAY(CINFN(JB)) !SR 11/28/19
                                    end do !SR 11/28/19
                                    if (LAST_JDAY > TMEND - 0.5 .and. ITER >= 3 .and. LAST_JDAY <= NXCIN2(JB)) then !near TMEND, no new data!SR 11/28/19
                                        open(INC(JB), FILE=CINFN(JB), STATUS="OLD", POSITION="APPEND") ! Open file and push pointer to end  !SR 11/28/19
                                        NXCIN1(JB) = TMEND + 1.0 ! Push input date past TMEND         !SR 11/28/19
                                        write(9911, "(A,I0,2(A,F0.4))") "MOVING ON: Input CIN", JB, " DAY= ", LAST_JDAY, " close to TMEND= ", TMEND !SR 11/28/19
                                    else !SR 11/28/19
                                        LAST_JDAY = NXCIN2(JB) ! Save the last date read            !SR 11/28/19
                                        open(INC(JB), FILE=CINFN(JB), STATUS="OLD") ! Open newly copied file and read    !SR 11/28/19
                                        read(INC(JB), "(A1)") INFORMAT !SR 11/28/19
                                        INCF(JB) = .false. !SR 11/28/19
                                        if (INFORMAT == "$") then
                                            INCF(JB) = .true.
                                        end if !SR 11/28/19
                                        if (INCF(JB)) then !SR 11/28/19
                                            read(INC(JB), "(/)") !SR 11/28/19
                                            read(INC(JB), *) NXCIN1(JB) ! Just read the date                 !SR 11/28/19
                                            do while (LAST_JDAY > NXCIN1(JB) .and. .not. EOF(INC(JB))) ! Get file ptr to previous position  !SR 11/28/19
                                                read(INC(JB), *) NXCIN1(JB) ! Just read the date                 !SR 11/28/19
                                            end do !SR 11/28/19
                                            if (EOF(INC(JB))) then
                                                backspace(INC(JB))
                                            end if ! A bit of insurance                 !SR 11/28/19
                                            read(INC(JB), *) NXCIN1(JB), (CINNX(INCN(JAC, JB), JB), JAC = 1, NACIN(JB)) ! Read new data point      !SR 11/28/19
                                        else !SR 11/28/19
                                            read(INC(JB), "(//F8.0)") NXCIN1(JB) ! Just read the date                 !SR 11/28/19
                                            do while (LAST_JDAY > NXCIN1(JB) .and. .not. EOF(INC(JB))) ! Get file ptr to previous position  !SR 11/28/19
                                                read(INC(JB), "(F8.0)") NXCIN1(JB) ! Just read the date                 !SR 11/28/19
                                            end do !SR 11/28/19
                                            if (EOF(INC(JB))) then
                                                backspace(INC(JB))
                                            end if ! A bit of insurance                 !SR 11/28/19
                                            read(INC(JB), "(1000F8.0)") NXCIN1(JB), (CINNX(INCN(JAC, JB), JB), JAC = 1, NACIN(JB)) ! Read new data  !SR 11/28/19
                                        end if !SR 11/28/19
                                    end if !SR 11/28/19
                                else ! Not waiting on this file. Stop run. File has no more data.     !SR 11/28/19
                                    write(W2ERR, "(2A/A,F0.4)") "ERROR-- End of branch input file ", CINFN(JB), "Simulation terminated at JDAY ", JDAY !SR 11/28/19
                                    write(*, "(2A/A,F0.4)") "ERROR-- End of branch input file ", CINFN(JB), "Simulation terminated at JDAY ", JDAY !SR 11/28/19
                                    stop !SR 11/28/19
                                end if !SR 11/28/19
                            end if !SR 11/28/19
                            8832 continue !SR 11/28/19
                        end do
!  IF(INCF(JB))THEN
!  READ (INC(JB),*) NXCIN1(JB),(CINNX(INCN(JAC,JB),JB),JAC=1,NACIN(JB))    
!  ELSE
!  READ (INC(JB),'(1000F8.0)') NXCIN1(JB),(CINNX(INCN(JAC,JB),JB),JAC=1,NACIN(JB))
!  ENDIF
!END DO
                        NXTVD = MIN(NXTVD, NXCIN1(JB))
                    end if
                end if
            end if

!**** Outflow

            if (DN_FLOW(JB) .and. NSTR(JB) > 0) then
                do while (JDAY >= NXQOT1(JB))
                    QSTR(1:NSTR(JB), JB) = QSTRNX(1:NSTR(JB), JB)
                    QSTRO(1:NSTR(JB), JB) = QSTRNX(1:NSTR(JB), JB)
                    NXQOT2(JB) = NXQOT1(JB)
                    if (OTQF(JB)) then
                        read(OTQ(JB), *) NXQOT1(JB), (QSTRNX(JS, JB), JS = 1, NSTR(JB))
                    else
                        read(OTQ(JB), "(10F8.0:/(8X,9F8.0))") NXQOT1(JB), (QSTRNX(JS, JB), JS = 1, NSTR(JB))
                    end if
                end do
                if (DYNSTRUC(JB) == "      ON") then
                    do while (JDAY >= NXDYNS(JB))
                        do J = 1, NJS
                            if (.not. ACTIVE_RULE_W2SELECTIVE(J, JB)) then
                                ESTR(JJS(J), JB) = NXESTRT(JJS(J), JB)
                            end if
                        end do
                        if (DYNEF(JB)) then
                            read(ODYNS(JB), *) NXDYNS(JB), (NXESTRT(JJS(J), JB), J = 1, NJS)
                        else
                            read(ODYNS(JB), "(10F8.0:/(8X,9F8.0))") NXDYNS(JB), (NXESTRT(JJS(J), JB), J = 1, NJS)
                        end if
                    end do
                    NXTVD = MIN(NXTVD, NXDYNS(JB))
                end if
                NXTVD = MIN(NXTVD, NXQOT1(JB))
            end if

!**** Distributed tributaries

            if (DIST_TRIBS(JB)) then

!****** Inflow

                do while (JDAY >= NXQDT1(JB))
                    QDTR(JB) = QDTRNX(JB)
                    QDTRO(JB) = QDTRNX(JB)
                    NXQDT2(JB) = NXQDT1(JB)

                    if (DTQF(JB)) then
                        read(DTQ(JB), *) NXQDT1(JB), QDTRNX(JB)
                    else
                        read(DTQ(JB), "(2F8.0)") NXQDT1(JB), QDTRNX(JB)
                    end if

!        READ (DTQ(JB),'(2F8.0)') NXQDT1(JB),QDTRNX(JB)
                end do
                NXTVD = MIN(NXTVD, NXQDT1(JB))

!****** Temperature

                do while (JDAY >= NXTDT1(JB))
                    TDTR(JB) = TDTRNX(JB)
                    TDTRO(JB) = TDTRNX(JB)
                    NXTDT2(JB) = NXTDT1(JB)

                    if (DTTF(JB)) then
                        read(DTT(JB), *) NXTDT1(JB), TDTRNX(JB)
                    else
                        read(DTT(JB), "(2F8.0)") NXTDT1(JB), TDTRNX(JB)
                    end if

!         READ (DTT(JB),'(2F8.0)') NXTDT1(JB),TDTRNX(JB)
                end do
                NXTVD = MIN(NXTVD, NXTDT1(JB))

!****** Constituent concentrations

                if (DTRIB_CONST(JB)) then
                    do while (JDAY >= NXCDT1(JB))
                        CDTR(DTCN(1:NACDT(JB), JB), JB) = CDTRNX(DTCN(1:NACDT(JB), JB), JB)
                        CDTRO(DTCN(1:NACDT(JB), JB), JB) = CDTRNX(DTCN(1:NACDT(JB), JB), JB)
!
! systdg - time series input
                        if (SYSTDG) then
                            if (DOBND) then
                                DO_SATD(JB) = EXP(7.7117 - 1.31403*LOG(TDTR(JB) + 45.93))*PALT(CUS(JB))
                                CDTR(NDO, JB) = CDTR(NDO, JB)*DO_SATD(JB)
                            end if
                            if (N2BND) then
                                EA = DEXP(2.3026D0*(7.5D0*TDEW(WBSEG(CUS(JB)))/(TDEW(WBSEG(CUS(JB))) + 237.3D0) + 0.6609D0))*0.001316
                                N2_SATD(JB) = 1.5568D06*0.79*(PALT(CUS(JB)) - EA)*(1.8816D-5 - 4.116D-7*TDTR(JB) + 4.6D-9*TDTR(JB)*TDTR(JB))
                                CDTR(NN2, JB) = CDTR(NN2, JB)*N2_SATD(JB)
                            end if
                            if (DGPBND) then
                                CDTR(NDGP, JB) = CDTR(NDGP, JB)*PALT(CUS(JB))
                            end if
                        end if

! systdg - time series input

                        NXCDT2(JB) = NXCDT1(JB)
                        if (DTCF(JB)) then
                            read(DTC(JB), *) NXCDT1(JB), (CDTRNX(DTCN(JAC, JB), JB), JAC = 1, NACDT(JB))
                        else
                            read(DTC(JB), "(1000F8.0)") NXCDT1(JB), (CDTRNX(DTCN(JAC, JB), JB), JAC = 1, NACDT(JB))
                        end if
                    end do
                    NXTVD = MIN(NXTVD, NXCDT1(JB))
                end if
            end if

!**** Precipitation

            if (PRECIPITATION(JW)) then
                do while (JDAY >= NXPR1(JB))
                    PR(JB) = PRNX(JB)
                    NXPR2(JB) = NXPR1(JB)

                    if (PRQF(JB)) then
                        read(PRE(JB), *) NXPR1(JB), PRNX(JB)
                    else
                        read(PRE(JB), "(2F8.0)") NXPR1(JB), PRNX(JB)
                    end if

!READ (PRE(JB),'(2F8.0)') NXPR1(JB),PRNX(JB)
                end do
                NXTVD = MIN(NXTVD, NXPR1(JB))

!****** Temperature

                do while (JDAY >= NXTPR1(JB))
                    TPR(JB) = TPRNX(JB)
                    NXTPR2(JB) = NXTPR1(JB)

                    if (PRTF(JB)) then
                        read(PRT(JB), *) NXTPR1(JB), TPRNX(JB)
                    else
                        read(PRT(JB), "(2F8.0)") NXTPR1(JB), TPRNX(JB)
                    end if

!READ (PRT(JB),'(2F8.0)') NXTPR1(JB),TPRNX(JB)
                end do
                NXTVD = MIN(NXTVD, NXTPR1(JB))

!****** Constituent concentrations

                if (PRECIP_CONST(JB)) then
                    do while (JDAY >= NXCPR1(JB))
                        CPR(PRCN(1:NACPR(JB), JB), JB) = CPRNX(PRCN(1:NACPR(JB), JB), JB)
!
! systdg - time series input
                        if (SYSTDG) then
                            if (DOBND) then
                                DO_SATP(JB) = EXP(7.7117 - 1.31403*LOG(TPR(JB) + 45.93))*PALT(CUS(JB))
                                CPR(NDO, JB) = CPR(NDO, JB)*DO_SATP(JB)
                            end if
                            if (N2BND) then
                                EA = DEXP(2.3026D0*(7.5D0*TDEW(WBSEG(CUS(JB)))/(TDEW(WBSEG(CUS(JB))) + 237.3D0) + 0.6609D0))*0.001316
                                N2_SATP(JB) = 1.5568D06*0.79*(PALT(CUS(JB)) - EA)*(1.8816D-5 - 4.116D-7*TPR(JB) + 4.6D-9*TPR(JB)*TPR(JB))
                                CPR(NN2, JB) = CPR(NN2, JB)*N2_SATP(JB)
                            end if
                            if (DGPBND) then
                                CPR(NDGP, JB) = CPR(NDGP, JB)*PALT(CUS(JB))
                            end if
                        end if

! systdg - time series input

                        NXCPR2(JB) = NXCPR1(JB)
                        if (PRCF(JB)) then
                            read(PRC(JB), *) NXCPR1(JB), (CPRNX(PRCN(JAC, JB), JB), JAC = 1, NACPR(JB))
                        else
                            read(PRC(JB), "(1000F8.0)") NXCPR1(JB), (CPRNX(PRCN(JAC, JB), JB), JAC = 1, NACPR(JB))
                        end if
                    end do
                    NXTVD = MIN(NXTVD, NXCPR1(JB))
                end if
            end if

!**** Upstream head conditions

            if (UH_EXTERNAL(JB)) then

!****** Elevations

                do while (JDAY >= NXEUH1(JB))
                    ELUH(JB) = ELUHNX(JB)
                    ELUHO(JB) = ELUHNX(JB)
                    NXEUH2(JB) = NXEUH1(JB)

                    if (EUHF(JB) > 0) then
                        read(UHE(JB), *) NXEUH1(JB), ELUHNX(JB)
                    else
                        read(UHE(JB), "(2F8.0)") NXEUH1(JB), ELUHNX(JB)
                    end if
!READ (UHE(JB),'(2F8.0)') NXEUH1(JB),ELUHNX(JB)
                end do
                NXTVD = MIN(NXTVD, NXEUH1(JB))

!****** Temperatures

                do while (JDAY >= NXTUH1(JB))
                    do K = 2, KMX - 1
                        TUH(K, JB) = TUHNX(K, JB)
                        TUHO(K, JB) = TUHNX(K, JB)
                    end do
                    NXTUH2(JB) = NXTUH1(JB)

                    if (TUHF(JB) == 1) then
                        read(UHT(JB), *) NXTUH1(JB), XX(1)
                        TUHNX(2:KB(US(JB)), JB) = XX(1)
                    else
                        if (TUHF(JB) == 2) then
                            read(UHT(JB), *) NXTUH1(JB), (TUHNX(K, JB), K = 2, KB(US(JB)))
                        else
                            read(UHT(JB), "(10F8.0:/(8X,9F8.0))") NXTUH1(JB), (TUHNX(K, JB), K = 2, KB(US(JB)))
                        end if
                    end if


!READ (UHT(JB),'(10F8.0:/(8X,9F8.0))') NXTUH1(JB),(TUHNX(K,JB),K=2,KB(US(JB)))
                end do
                NXTVD = MIN(NXTVD, NXTUH1(JB))

!****** Constituent concentrations

                if (CONSTITUENTS) then
                    do while (JDAY >= NXCUH1(JB))
                        do K = 2, KMX - 1
                            CUH(K, CN(1:NAC), JB) = CUHNX(K, CN(1:NAC), JB)
                            CUHO(K, CN(1:NAC), JB) = CUHNX(K, CN(1:NAC), JB)
                        end do
                        NXCUH2(JB) = NXCUH1(JB)

                        if (CUHF(JB) == 1) then
                            read(UHC(JB), *) NXCUH1(JB), (XX(CN(JAC)), JAC = 1, NAC)
                            do JAC = 1, NAC
                                CUHNX(2:KB(US(JB)), CN(JAC), JB) = XX(CN(JAC))
                            end do
                        else

                            do JAC = 1, NAC
                                if (CUHF(JB) == 2) then
                                    read(UHC(JB), *) NXCUH1(JB), (CUHNX(K, CN(JAC), JB), K = 2, KB(US(JB)))
                                else
                                    read(UHC(JB), "(10F8.0:/(8X,9F8.0))") NXCUH1(JB), (CUHNX(K, CN(JAC), JB), K = 2, KB(US(JB)))
                                end if
                            end do
                        end if
!DO JAC=1,NAC
!  IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (UHC(JB),'(10F8.0:/(8X,9F8.0))') NXCUH1(JB),(CUHNX(K,CN(JAC),JB),   &
!                                                    K=2,KB(US(JB)))
!END DO
                    end do
                    NXTVD = MIN(NXTVD, NXCUH1(JB))
                end if
            end if

!**** Downstream head

            if (DH_EXTERNAL(JB)) then

!****** Elevation

                do while (JDAY >= NXEDH1(JB))
                    ELDH(JB) = ELDHNX(JB)
                    ELDHO(JB) = ELDHNX(JB)
                    NXEDH2(JB) = NXEDH1(JB)

                    if (EDHF(JB) > 0) then
                        read(DHE(JB), *) NXEDH1(JB), ELDHNX(JB)
                    else
                        read(DHE(JB), "(2F8.0)") NXEDH1(JB), ELDHNX(JB)
                    end if
!READ (DHE(JB),'(2F8.0)') NXEDH1(JB),ELDHNX(JB)
                end do
                NXTVD = MIN(NXTVD, NXEDH1(JB))

!****** Temperature

                do while (JDAY >= NXTDH1(JB))
                    do K = 2, KMX - 1
                        TDH(K, JB) = TDHNX(K, JB)
                        TDHO(K, JB) = TDHNX(K, JB)
                    end do
                    NXTDH2(JB) = NXTDH1(JB)
                    if (TDHF(JB) == 1) then
                        read(DHT(JB), *) NXTDH1(JB), XX(1)
                        TDHNX(2:KB(DS(JB)), JB) = XX(1)
                    else
                        if (TDHF(JB) == 2) then
                            read(DHT(JB), *) NXTDH1(JB), (TDHNX(K, JB), K = 2, KB(DS(JB)))
                        else
                            read(DHT(JB), "(10F8.0:/(8X,9F8.0))") NXTDH1(JB), (TDHNX(K, JB), K = 2, KB(DS(JB)))
                        end if
                    end if
!READ (DHT(JB),'(10F8.0:/(8X,9F8.0))') NXTDH1(JB),(TDHNX(K,JB),K=2,KB(DS(JB)))
                end do
                NXTVD = MIN(NXTVD, NXTDH1(JB))

!****** Constituents

                if (CONSTITUENTS) then
                    do while (JDAY >= NXCDH1(JB))
                        do K = 2, KMX - 1
                            CDH(K, CN(1:NAC), JB) = CDHNX(K, CN(1:NAC), JB)
                            CDHO(K, CN(1:NAC), JB) = CDHNX(K, CN(1:NAC), JB)
                        end do
                        NXCDH2(JB) = NXCDH1(JB)
                        if (CDHF(JB) == 1) then
                            read(DHC(JB), *) NXCDH1(JB), (XX(CN(JAC)), JAC = 1, NAC)
                            do JAC = 1, NAC
                                CDHNX(2:KB(DS(JB)), CN(JAC), JB) = XX(CN(JAC))
                            end do
                        else
                            do JAC = 1, NAC
                                if (CDHF(JB) == 2) then
                                    read(DHC(JB), *) NXCDH1(JB), (CDHNX(K, CN(JAC), JB), K = 2, KB(DS(JB)))
                                else
                                    read(DHC(JB), "(10F8.0:/(8X,9F8.0))") NXCDH1(JB), (CDHNX(K, CN(JAC), JB), K = 2, KB(DS(JB)))
                                end if
                            end do
                        end if
!DO JAC=1,NAC
!  IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (DHC(JB),'(10F8.0:/(8X,9F8.0))') NXCDH1(JB),(CDHNX(K,CN(JAC),JB),   &
!                                                    K=2,KB(DS(JB)))
!END DO
                    end do
                    NXTVD = MIN(NXTVD, NXCDH1(JB))
                end if
            end if
        end do
! ATMOSPHERIC LOADING
        if (CONSTITUENTS) then
            if (ATM_DEPOSITION(JW)) then
                do while (JDAY >= NXATMD(JW))
                    NXATMD2(JW) = NXATMD(JW)
                    ATM_DEP_LOADING = ATM_DEP_LOADINGNX
                    ATM_DEP_LOADING0 = ATM_DEP_LOADINGNX
                    if (ATMDEPCSV) then
                        read(ATMDEP(JW), *) NXATMD(JW), (ATM_DEP_LOADINGNX(ATMDCN(JAC, JW), JW), JAC = 1, NACATD(JW))
                    else
                        read(ATMDEP(JW), "(100F8.0)") NXATMD(JW), (ATM_DEP_LOADINGNX(ATMDCN(JAC, JW), JW), JAC = 1, NACATD(JW))
                    end if
                end do
                NXTVD = MIN(NXTVD, NXATMD(JW))
            end if
        end if
    end do

! Gate height opening

    if (GATES) then
        do while (JDAY >= NXQGT)
            nxqgt2 = nxqgt
            where (DYNGTC == "     ZGT")
                EGT = BGTNX
                egto = bgtnx
                BGT = 1.0
                G1GT = 1.0
                G2GT = 1.0
            elsewhere
                BGT = BGTNX
                bgto = bgtnx
            end where
            if (GATEF) then
                read(GTQ, *) NXQGT, (BGTNX(JG), JG = 1, NGT)
            else
                read(GTQ, "(1000F8.0)") NXQGT, (BGTNX(JG), JG = 1, NGT)
            end if
        end do
        NXTVD = MIN(NXTVD, NXQGT)
    end if

! Pipe reduction factor

    if (PIPES .and. iopenpipe == 1) then
        do while (JDAY >= NXQPT)
            BP = BPNX
            read(PIPED, "(1000F8.0)") NXQPT, (BPNX(J), J = 1, NPI)
        end do
        NXTVD = MIN(NXTVD, NXQPT)
    end if

! DYNAMIC PUMPS

    if (PUMPS) then

        do J = 1, NPU
            if (DYNPUMP(J) == "      ON") then
                do while (JDAY >= NXPUMP(J))
                    EPU(J) = EPU2(J)
                    EONPU(J) = EONPU2(J)
                    EOFFPU(J) = EOFFPU2(J)
                    QPU(J) = QPU2(J)
                    if (DYNPUMPF(J)) then
                        read(PUMPD(J), *) NXPUMP(J), EPU2(J), EONPU2(J), EOFFPU2(J), QPU2(J)
                    else
                        read(PUMPD(J), "(1000F8.0)") NXPUMP(J), EPU2(J), EONPU2(J), EOFFPU2(J), QPU2(J)
                    end if
!     READ (PUMPD(J),'(1000F8.0)') NXPUMP(J),EPU2(J),EONPU2(J),EOFFPU2(J),QPU2(J)
                end do
                NXTVD = MIN(NXTVD, NXPUMP(J))
            end if
        end do
    end if

! Dead sea case

    do JW = 1, NWB
        if (NO_INFLOW(JW)) then
            QIN(BS(JW):BE(JW)) = 0.0
            QINO(BS(JW):BE(JW)) = 0.0
            QIND(BS(JW):BE(JW)) = 0.0
            QINNX(BS(JW):BE(JW)) = 0.0
            QDTR(BS(JW):BE(JW)) = 0.0
            QDTRO(BS(JW):BE(JW)) = 0.0
            QDTRNX(BS(JW):BE(JW)) = 0.0
            PR(BS(JW):BE(JW)) = 0.0
            PRNX(BS(JW):BE(JW)) = 0.0
        end if
        if (NO_OUTFLOW(JW)) then
            QSTR(:, BS(JW):BE(JW)) = 0.0
            QSTRO(:, BS(JW):BE(JW)) = 0.0
            QSTRNX(:, BS(JW):BE(JW)) = 0.0
        end if
    end do
    where (NO_WIND)
        WIND = 0.0
        WINDO = 0.0
        WINDNX = 0.0
    end where
    where (READ_RADIATION .and. NO_HEAT)
        SRON = 0.0
        SROO = 0.0
        SRONX = 0.0
    end where
    if (ANY(NO_INFLOW)) then
        QTR = 0.0
        QTRO = 0.0
        QTRNX = 0.0
        QWD = 0.0
        QWDO = 0.0
        QWDNX = 0.0
    end if
    return

!***********************************************************************************************************************************
!**                                              I N T E R P O L A T E  I N P U T S                                               **
!***********************************************************************************************************************************

    entry INTERPOLATE_INPUTS()

! Meteorological/light extinction data

    do JW = 1, NWB
        if (INTERP_METEOROLOGY(JW)) then
            RATIO = (NXMET1(JW) - JDAY)/(NXMET1(JW) - NXMET2(JW))
            TDEW(JW) = (1.0 - RATIO)*TDEWNX(JW) + RATIO*TDEWO(JW)
            WIND(JW) = (1.0 - RATIO)*WINDNX(JW) + RATIO*WINDO(JW)
! CONVERT PHIO AND PHINX TO LESS THAN 2*PI     SW 2/13/15
            do while (PHIO(JW) > 2.*PI)
                PHIO(JW) = PHIO(JW) - 2.*PI
            end do
            do while (PHINX(JW) > 2.*PI)
                PHINX(JW) = PHINX(JW) - 2.*PI
            end do
            if (PHIO(JW) - PHINX(JW) > PI) then
                PHI(JW) = (1.0 - RATIO)*(PHINX(JW) + 2.0*PI) + RATIO*PHIO(JW)
            else
                if (PHIO(JW) - PHINX(JW) < -PI) then ! WX 2/13/15
                    PHI(JW) = (1.0 - RATIO)*PHINX(JW) + RATIO*(PHIO(JW) + 2.0*PI) ! WX 2/13/15
                else
                    PHI(JW) = (1.0 - RATIO)*PHINX(JW) + RATIO*PHIO(JW)
                end if
            end if

!IF (ABS(PHIO(JW)-PHINX(JW)) > PI) THEN
!  PHI(JW) = (1.0-RATIO)*(PHINX(JW)+2.0*PI)+RATIO*PHIO(JW)
!ELSE
!  PHI(JW) = (1.0-RATIO)*PHINX(JW)+RATIO*PHIO(JW)
!END IF
            TAIR(JW) = (1.0 - RATIO)*TAIRNX(JW) + RATIO*TAIRO(JW)
            CLOUD(JW) = (1.0 - RATIO)*CLOUDNX(JW) + RATIO*CLOUDO(JW)
!
            if (SYSTDG) then
! systdg - time series input PALT_JW
                PALT_JW(JW) = (1.0 - RATIO)*PALT_JWNX(JW) + RATIO*PALT_JWO(JW)
                if (PALT_JW(JW) <= 0.0) then
                    PALT_JW(JW) = 760.0
                end if
                do I = US(BS(JW)) - 1, DS(BE(JW)) + 1
                    PALT(I) = PALT_JW(JW)/760.0*(1.0 - ELWS_INI(I)/1000.0/44.3)**5.25 ! systdg - time series input PALT
                end do
!PALT(:) = PALT_JW(JW)/760.0*(1.0-ELWS_INI(:)/1000.0/44.3)**5.25
! systdg - time series input PALT_JW
            end if

            if (READ_RADIATION(JW)) then
                SRON(JW) = (1.0 - RATIO)*SRONX(JW) + RATIO*SROO(JW)
            end if
        end if
        if (READ_EXTINCTION(JW) .and. INTERP_EXTINCTION(JW)) then ! 6/30/15 SW
            RATIO = (NXEXT1(JW) - JDAY)/(NXEXT1(JW) - NXEXT2(JW))
            EXH2O(JW) = (1.0 - RATIO)*EXTNX(JW) + RATIO*EXTO(JW)
        end if
! ATMOSPHERIC DEPOSITION
        if (CONSTITUENTS) then
            if (ATM_DEPOSITION(JW) .and. ATM_DEPOSITION_INTERPOLATION(JW) == "      ON") then
                RATIO = (NXATMD(JW) - JDAY)/(NXATMD(JW) - NXATMD2(JW))
                ATM_DEP_LOADING(:, JW) = (1.0 - RATIO)*ATM_DEP_LOADINGNX(:, JW) + RATIO*ATM_DEP_LOADING0(:, JW)
            end if
        end if

    end do

! Withdrawals

    if (NWD > 0) then
        QRATIO = (NXQWD1 - JDAY)/(NXQWD1 - NXQWD2)
        do JWD = 1, NWD
            if (INTERP_WITHDRAWAL(JWD)) then
                QWD(JWD) = (1.0 - QRATIO)*QWDNX(JWD) + QRATIO*QWDO(JWD)
            end if
        end do
    end if

! Gates  adding interpolation cb 8/13/2010  
    if (gates) then
        QRATIO = (NXQgt - JDAY)/(NXQgt - NXQgt2)
        do Jg = 1, ngt
            if (INTERP_gate(Jg)) then
                if (DYNGTC(Jg) == "     ZGT") then
                    egt(jg) = (1.0 - QRATIO)*bgtNX(jg) + QRATIO*egtO(jg)
                else
                    bgt(jg) = (1.0 - QRATIO)*bgtNX(Jg) + QRATIO*bgtO(Jg)
                end if
            end if
        end do
    end if

! Tributaries

    if (NTR > 0) then
        do JT = 1, NTR
            if (INTERP_TRIBS(JT)) then
                QRATIO = (NXQTR1(JT) - JDAY)/(NXQTR1(JT) - NXQTR2(JT))
                TRATIO = (NXTTR1(JT) - JDAY)/(NXTTR1(JT) - NXTTR2(JT))
                if (TRIB_CONST(JT)) then
                    CRATIO = (NXCTR1(JT) - JDAY)/(NXCTR1(JT) - NXCTR2(JT))
                end if
                QTR(JT) = (1.0 - QRATIO)*QTRNX(JT) + QRATIO*QTRO(JT)
                TTR(JT) = (1.0 - TRATIO)*TTRNX(JT) + TRATIO*TTRO(JT)
                CTR(TRCN(1:NACTR(JT), JT), JT) = (1.0 - CRATIO)*CTRNX(TRCN(1:NACTR(JT), JT), JT) + CRATIO*CTRO(TRCN(1:NACTR(JT), JT), JT)
!
! systdg - time series input
                if (SYSTDG) then
                    if (DOBND) then
                        DO_SAT(JT) = EXP(7.7117 - 1.31403*LOG(TTR(JT) + 45.93))*PALT(ITR(JT))
                        CTR(NDO, JT) = CTR(NDO, JT)*DO_SAT(JT)
                    end if
                    if (N2BND) then
                        EA = DEXP(2.3026D0*(7.5D0*TDEW(WBSEG(ITR(JT)))/(TDEW(WBSEG(ITR(JT))) + 237.3D0) + 0.6609D0))*0.001316
                        N2_SAT(JT) = 1.5568D06*0.79*(PALT(ITR(JT)) - EA)*(1.8816D-5 - 4.116D-7*TTR(JT) + 4.6D-9*TTR(JT)*TTR(JT))
                        CTR(NN2, JT) = CTR(NN2, JT)*N2_SAT(JT)
                    end if
!
                    if (DGPBND) then
                        CTR(NDGP, JT) = CTR(NDGP, JT)*PALT(ITR(JT))
                    end if
                end if

! systdg - time series input
            end if
        end do
    end if
!
! systdg - time series input TWE
    if (SYSTDG) then
        if (TWETSC == "      ON") then
            TWERATIO = (NXTWE1 - JDAY)/(NXTWE1 - NXTWE2)
            TWE_TS = (1.0 - TWERATIO)*TWE_TSNX + TWERATIO*TWE_TSO
        end if
    end if
! systdg - time series input TWE

! Branch related inputs

    do JB = 1, NBR

!** Inflow

        if (UP_FLOW(JB)) then
            if (.not. INTERNAL_FLOW(JB) .and. .not. DAM_INFLOW(JB)) then !TC 08/03/04 RA 1/13/06
                if (INTERP_INFLOW(JB)) then
                    QRATIO = (NXQIN1(JB) - JDAY)/(NXQIN1(JB) - NXQIN2(JB))
                    TRATIO = (NXTIN1(JB) - JDAY)/(NXTIN1(JB) - NXTIN2(JB))
                    if (INFLOW_CONST(JB)) then
                        CRATIO = (NXCIN1(JB) - JDAY)/(NXCIN1(JB) - NXCIN2(JB))
                    end if
                    QIND(JB) = (1.0 - QRATIO)*QINNX(JB) + QRATIO*QINO(JB)
                    TIND(JB) = (1.0 - TRATIO)*TINNX(JB) + TRATIO*TINO(JB)
                    CIND(INCN(1:NACIN(JB), JB), JB) = (1.0 - CRATIO)*CINNX(INCN(1:NACIN(JB), JB), JB) + CRATIO*CINO(INCN(1:NACIN(JB), JB), JB)
!
! systdg - time series input
                    if (SYSTDG) then
                        if (DOBND) then
                            DO_SATJ(JB) = EXP(7.7117 - 1.31403*LOG(TIND(JB) + 45.93))*PALT(CUS(JB))
                            CIND(NDO, JB) = CIND(NDO, JB)*DO_SATJ(JB)
                        end if
                        if (N2BND) then
                            EA = DEXP(2.3026D0*(7.5D0*TDEW(WBSEG(CUS(JB)))/(TDEW(WBSEG(CUS(JB))) + 237.3D0) + 0.6609D0))*0.001316

                            N2_SATJ(JB) = 1.5568D06*0.79*(PALT(CUS(JB)) - EA)*(1.8816D-5 - 4.116D-7*TIND(JB) + 4.6D-9*TIND(JB)*TIND(JB))
                            CIND(NN2, JB) = CIND(NN2, JB)*N2_SATJ(JB)
                        end if
!
                        if (DGPBND) then
                            CIND(NDGP, JB) = CIND(NDGP, JB)*PALT(CUS(JB))
                        end if
                    end if
! systdg - time series input
                end if
            end if
        end if

!** Outflow

        if (DN_FLOW(JB) .and. NSTR(JB) > 0) then
            QRATIO = (NXQOT1(JB) - JDAY)/(NXQOT1(JB) - NXQOT2(JB))
            do JS = 1, NSTR(JB)
                if (INTERP_OUTFLOW(JS, JB)) then
                    QSTR(JS, JB) = (1.0 - QRATIO)*QSTRNX(JS, JB) + QRATIO*QSTRO(JS, JB)
                end if
            end do
        end if

!** Distributed tributaries

        if (DIST_TRIBS(JB)) then
            if (INTERP_DTRIBS(JB)) then
                QRATIO = (NXQDT1(JB) - JDAY)/(NXQDT1(JB) - NXQDT2(JB))
                TRATIO = (NXTDT1(JB) - JDAY)/(NXTDT1(JB) - NXTDT2(JB))
                if (DTRIB_CONST(JB)) then
                    CRATIO = (NXCDT1(JB) - JDAY)/(NXCDT1(JB) - NXCDT2(JB))
                end if
                QDTR(JB) = (1.0 - QRATIO)*QDTRNX(JB) + QRATIO*QDTRO(JB)
                TDTR(JB) = (1.0 - TRATIO)*TDTRNX(JB) + TRATIO*TDTRO(JB)
                CDTR(DTCN(1:NACDT(JB), JB), JB) = (1.0 - CRATIO)*CDTRNX(DTCN(1:NACDT(JB), JB), JB) + CRATIO*CDTRO(DTCN(1:NACDT(JB), JB), JB)
!
! systdg - time series input
                if (SYSTDG) then
                    if (DOBND) then
                        DO_SATD(JB) = EXP(7.7117 - 1.31403*LOG(TDTR(JB) + 45.93))*PALT(CUS(JB))
                        CDTR(NDO, JB) = CDTR(NDO, JB)*DO_SATD(JB)
                    end if
                    if (N2BND) then
                        EA = DEXP(2.3026D0*(7.5D0*TDEW(WBSEG(CUS(JB)))/(TDEW(WBSEG(CUS(JB))) + 237.3D0) + 0.6609D0))*0.001316
                        N2_SATD(JB) = 1.5568D06*0.79*(PALT(CUS(JB)) - EA)*(1.8816D-5 - 4.116D-7*TDTR(JB) + 4.6D-9*TDTR(JB)*TDTR(JB))
                        CDTR(NN2, JB) = CDTR(NN2, JB)*N2_SATD(JB)
                    end if
!
                    if (DGPBND) then
                        CDTR(NDGP, JB) = CDTR(NDGP, JB)*PALT(CUS(JB))
                    end if
                end if
! systdg - time series input
            end if
        end if

!** Upstream head

        if (UH_EXTERNAL(JB)) then
            if (INTERP_HEAD(JB)) then
                HRATIO = (NXEUH1(JB) - JDAY)/(NXEUH1(JB) - NXEUH2(JB))
                TRATIO = (NXTUH1(JB) - JDAY)/(NXTUH1(JB) - NXTUH2(JB))
                if (CONSTITUENTS) then
                    CRATIO = (NXCUH1(JB) - JDAY)/(NXCUH1(JB) - NXCUH2(JB))
                end if
                ELUH(JB) = (1.0 - HRATIO)*ELUHNX(JB) + HRATIO*ELUHO(JB)
                do K = 2, KMX - 1
                    TUH(K, JB) = (1.0 - TRATIO)*TUHNX(K, JB) + TRATIO*TUHO(K, JB)
                    CUH(K, CN(1:NAC), JB) = (1.0 - CRATIO)*CUHNX(K, CN(1:NAC), JB) + CRATIO*CUHO(K, CN(1:NAC), JB)
                end do
            end if
        end if

!** Downstream head

        if (DH_EXTERNAL(JB)) then
            if (INTERP_HEAD(JB)) then
                HRATIO = (NXEDH1(JB) - JDAY)/(NXEDH1(JB) - NXEDH2(JB))
                TRATIO = (NXTDH1(JB) - JDAY)/(NXTDH1(JB) - NXTDH2(JB))
                if (CONSTITUENTS) then
                    CRATIO = (NXCDH1(JB) - JDAY)/(NXCDH1(JB) - NXCDH2(JB))
                end if
                ELDH(JB) = (1.0 - HRATIO)*ELDHNX(JB) + HRATIO*ELDHO(JB)
                do K = 2, KMX - 1
                    TDH(K, JB) = (1.0 - TRATIO)*TDHNX(K, JB) + TRATIO*TDHO(K, JB)
                    CDH(K, CN(1:NAC), JB) = (1.0 - CRATIO)*CDHNX(K, CN(1:NAC), JB) + CRATIO*CDHO(K, CN(1:NAC), JB)
                end do
            end if
        end if
    end do
    return
    entry DEALLOCATE_TIME_VARYING_DATA()
!
! systdg
    deallocate(DO_SAT, N2_SAT)
    deallocate(DO_SATJ, N2_SATJ)
    deallocate(DO_SATD, N2_SATD)
    deallocate(DO_SATP, N2_SATP)
!
    deallocate(NXQTR1, NXTTR1, NXCTR1, NXQIN1, NXTIN1, NXCIN1, NXQDT1, NXTDT1, NXCDT1, NXPR1, NXTPR1, NXCPR1, NXEUH1, NXTUH1)
    deallocate(NXCUH1, NXEDH1, NXTDH1, NXCDH1, NXQOT1, NXMET1, NXQTR2, NXTTR2, NXCTR2, NXQIN2, NXTIN2, NXCIN2, NXQDT2, NXTDT2)
    deallocate(NXCDT2, NXPR2, NXTPR2, NXCPR2, NXEUH2, NXTUH2, NXCUH2, NXEDH2, NXTDH2, NXCDH2, NXQOT2, NXMET2, WSCNX, DYNPUMPF)
    deallocate(QDTRO, TDTRO, ELUHO, ELDHO, QWDO, QTRO, TTRO, QINO, TINO, QDTRNX, TDTRNX, PRNX, TPRNX, ELUHNX)
    deallocate(ELDHNX, QWDNX, QTRNX, TTRNX, QINNX, TINNX, SROO, TAIRO, TDEWO, CLOUDO, PHIO, WINDO, TAIRNX, BGTNX, BPNX, PALT_JWO, PALT_JWNX) ! systdg  PALT_JWO, PALT_JWNX
    deallocate(TDEWNX, CLOUDNX, PHINX, WINDNX, SRONX, TRQ, TRT, TRC, INQ, DTQ, PRE, UHE, DHE, INFT)
    deallocate(DTT, PRT, UHT, DHT, INC, DTC, PRC, UHC, DHC, OTQ, MET, EXT, EXTNX, EXTO, EXTF)
    deallocate(NXEXT1, NXEXT2, CTRO, CINO, QOUTO, CDTRO, TUHO, TDHO, QSTRO, CTRNX, CINNX, QOUTNX, CDTRNX, CPRNX)
    deallocate(TUHNX, TDHNX, QSTRNX, CUHO, CDHO, CUHNX, CDHNX, INCF, DTCF, PRCF, METF, ODYNS, NXDYNS, NXESTRT, DYNEF, JJS, PRTF, PRQF)
    deallocate(EUHF, TUHF, CUHF, EDHF, TDHF, CDHF, XX)
    deallocate(ATMDEP, ATM_DEP_LOADINGNX, NXATMD, NXATMD2, ATM_DEP_LOADING0, NXATMDEP)
    deallocate(INFLOW_CONST, TRIB_CONST, DTRIB_CONST, PRECIP_CONST, PUMPD, NXPUMP, EPU2, EONPU2, EOFFPU2, QPU2, OTQF, TRCF, TRQF, TRTF, DTTF, DTQF, INQF, INTF)
    return
end subroutine TIME_VARYING_DATA

!***********************************************************************************************************************************
!**                                       F U N C T I O N   G E T _ L A S T _ J D A Y                                             **
!***********************************************************************************************************************************

! Function to read a standard input file and return the last available JDAY                           Entire function:  !SR 11/28/19
! File name is passed to the function, and file is guaranteed to exist.


real function GET_LAST_JDAY(FNAME)
    use MAIN, only: CON ! Control file already read and closed
    real :: LAST_JDAY
    character(len=1) :: INFORMAT
    character(len=72) :: FNAME

    open(CON, FILE=FNAME, STATUS="OLD") ! File was previously tested for existence
    read(CON, "(A1)") INFORMAT ! Assumption that 3-line header is in place, with at least one line of data
    read(CON, "(/)")
    if (INFORMAT == "$") then
        do while (.not. EOF(CON))
            read(CON, *) LAST_JDAY ! Free format, probably comma-delimited
        end do
    else
        do while (.not. EOF(CON))
            read(CON, "(F8.0)") LAST_JDAY ! Fixed format, needed if fixed-format file did not have space after JDAY
        end do
    end if
    close(CON)
    GET_LAST_JDAY = LAST_JDAY
end function GET_LAST_JDAY


!***********************************************************************************************************************************
!**                                   S U B R O U T I N E   P R I N T _ E R R O R _ A N D _ S T O P                               **
!***********************************************************************************************************************************

! Function to print a customized error message from a recent system command                           Entire function:  !SR 11/28/19
! Executing this function will also stop the program


subroutine PRINT_ERROR_AND_STOP(TYPE_TEXT, TYPE_INDEX)
    use IFPORT ! For GETLASTERRORXX()

    integer :: TYPE_INDEX, IRESULT
    character(len=3) :: TYPE_TEXT
    character(len=8) :: ERR_TEXT

    IRESULT = GETLASTERRORQQ() ! USE GETLASTERRORQQ TO GET THE ERROR CODE
    write(ERR_TEXT, "(A,I0)") TYPE_TEXT, TYPE_INDEX
    write(*, *) TRIM(ERR_TEXT), ": ERROR EXECUTING BATCH PROGRAM: ERROR CODE:", IRESULT
    write(9911, *) TRIM(ERR_TEXT), ": ERROR EXECUTING BATCH PROGRAM: ERROR CODE:", IRESULT
    close(9911)
    stop
end subroutine PRINT_ERROR_AND_STOP
