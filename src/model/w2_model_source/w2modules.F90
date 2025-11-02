module MSCLIB
    integer :: HTHREAD
    logical :: STOP_PUSHED, STOPPED, RESTART_PUSHED, RESTART_EXISTS
! INCLUDE "RESOURCE.FD"
! INTERFACE
!   FUNCTION $BEGINTHREADEX (SECURITY,STACK_SIZE,START_ADDRESS,ARGLIST,INITFLAG,THRDADDR)
!     USE DFWINTY, RENAMED => DLT
!     !DEC$ ATTRIBUTES C,ALIAS : "__BEGINTHREADEX" :: $BEGINTHREADEX
!     !DEC$ ATTRIBUTES REFERENCE,ALLOW_NULL        :: SECURITY
!     !DEC$ ATTRIBUTES REFERENCE,IGNORE_LOC        :: THRDADDR
!     INTEGER(UINT)                                :: $BEGINTHREADEX
!     INTEGER(UINT),               INTENT(IN)      :: STACK_SIZE, INITFLAG
!     INTEGER(PVOID),              INTENT(IN)      :: START_ADDRESS, ARGLIST
!     INTEGER(UINT),               INTENT(OUT)     :: THRDADDR
!     TYPE(T_SECURITY_ATTRIBUTES), INTENT(IN)      :: SECURITY
!   END FUNCTION $BEGINTHREADEX
! END INTERFACE
! INTERFACE
!   SUBROUTINE $ENDTHREADEX (RETVAL)
!     USE DFWINTY, RENAMED => DLT
!     !DEC$ ATTRIBUTES C, ALIAS : "__ENDTHREADEX" :: $ENDTHREADEX
!     INTEGER(UINT), INTENT(IN) :: RETVAL
!   END SUBROUTINE $ENDTHREADEX
! END INTERFACE
end module MSCLIB

module PREC
    integer, parameter :: I2 = SELECTED_INT_KIND(3)
    integer, parameter :: R8 = SELECTED_REAL_KIND(15)
end module PREC

module RSTART
    use PREC
    real(R8) :: DLTS, CURMAX, DLTFF, DLTMAXX ! SW 7/13/2010
    integer :: RSODP, DLTDP, TSRDP, WDODP, CUF, RSO = 31
    integer, allocatable, dimension(:) :: SNPDP, VPLDP, CPLDP, PRFDP, SCRDP, SPRDP, FLXDP, NSPRF
!REAL                                               :: NXTMRS, NXTMWD, NXTMTS
    real(R8) :: NXTMRS, NXTMWD, NXTMTS, NXTMWD_sec ! cb 4/6/17
    real, allocatable, dimension(:) :: NXTMSN, NXTMPR, NXTMSP, NXTMCP, NXTMVP, NXTMSC, NXTMFL
    real(R8), allocatable, dimension(:) :: SBKT, ELTMF
    real(R8), allocatable, dimension(:,:) :: TSSUH2, TSSDH2, SAVH2, SAVHR, SU, SW, SAZ
    real(R8), allocatable, dimension(:,:,:) :: CSSUH2, CSSDH2
    real(R8) :: ELTM
    real(R8), allocatable, dimension(:) :: VOLIN, VOLOUT, VOLUH, VOLDH, VOLPR, VOLTRB, VOLDT, VOLWD, VOLEV, VOLICE, ICEBANK
    real(R8), allocatable, dimension(:) :: VOLSBR, VOLTBR, VOLSR, VOLTR
    real(R8), allocatable, dimension(:) :: TSSEV, TSSPR, TSSTR, TSSDT, TSSWD, TSSUH, TSSDH, TSSIN, TSSOUT
    real(R8), allocatable, dimension(:) :: TSSS, TSSB, TSSICE
    real(R8), allocatable, dimension(:) :: ESBR, ETBR, EBRI, SZ
    real(R8), allocatable, dimension(:,:) :: CMBRT
end module RSTART

module GLOBAL
    use PREC
    real(4) :: W2VER = 4.5
    real(R8), parameter :: DAY = 86400.0D0, NONZERO = 1.0D-20, REFL = 0.94D0, FRAZDZ = 0.14D0, DZMIN = 1.4D-7
    real(R8), parameter :: AZMIN = 1.4D-6, DZMAX = 1.0D3, RHOW = 1000.0D0
    real(R8) :: DLT, DLTMIN, DLTTVD
    real(R8) :: BETABR, START, HMAX2, CURRENT
    real(R8), pointer, dimension(:,:) :: U, W, T2, AZ, RHO, ST, SB
    real(R8), pointer, dimension(:,:) :: DLTLIM, VSH, ADMX, DM, ADMZ, HDG, HPG, GRAV
    real(R8), target, allocatable, dimension(:,:) :: T1, TSS
    real(R8), target, allocatable, dimension(:,:,:) :: C1, C2, C1S, CSSB, CSSK
    real(R8), target, allocatable, dimension(:,:,:) :: KF, CD
    real(R8), target, allocatable, dimension(:,:,:) :: HYD
    real, target, allocatable, dimension(:,:,:,:) :: AF, EF
    real(R8), allocatable, dimension(:) :: ICETH, ELKT, HMULT, CMULT, CDMULT, WIND2, AZMAX, PALT, Z0
    real, allocatable, dimension(:) :: TN_SEDSOD_NH4, NH3GASLOSS, TP_SEDSOD_PO4, TPOUT, TPTRIB, TPDTRIB, TPWD, TPPR, TPIN, TNOUT, TNTRIB, TNDTRIB, TNWD, TNPR, TNIN, ATMDEP_P, ATMDEP_N !TP_SEDBURIAL,TN_SEDBURIAL,
    real(R8), allocatable, dimension(:,:) :: QSS, VOLUH2, VOLDH2, QUH1, QDH1, UXBR, UYBR, VOL
    real, allocatable, dimension(:,:,:) :: ALLIM, APLIM, ANLIM, ASLIM, KFS
    real, allocatable, dimension(:,:,:) :: ELLIM, EPLIM, ENLIM, ESLIM
    integer :: W2ERR = 33, WRN = 32
    integer :: IMX, KMX, NBR, NTR, NWD, NWB, NCT, NBOD, NTR1 ! SW 2/17/2021
    integer :: NST, NSP, NGT, NPI, NPU, NWDO, NIKTSR, NUNIT
    integer :: JW, JB, JC, IU, ID, KT, I, JJB
    integer :: NOD, NDC = 27, NAL, NSS, NHY = 15, NFL = 142, NEP, NEPT
    integer :: NZP, NZPT, JZ, NZOOS, NZOOE, NMC, NMCT ! number of zooplankton groups, CONSTIUENT NUMBER FOR ZOOPLANKTON, START AND END
    integer :: NGCS, NGCE
    integer, pointer, dimension(:) :: SNP, PRF, VPL, CPL, SPR, FLX, FLX2, SPRV ! SE 9/28/2018
    integer, allocatable, dimension(:) :: BS, BE, US, CUS, DS, JBDN
    integer, allocatable, dimension(:) :: KB, KTI, SKTI, KTWB, KBMIN, CDHS
    integer, allocatable, dimension(:) :: UHS, DHS, UQB, DQB
    integer, target, allocatable, dimension(:,:) :: OPT
    integer, allocatable, dimension(:) :: NBODC, NBODN, NBODP ! CB 6/6/10
    logical, allocatable, dimension(:) :: ICE, ICE_CALC, LAYERCHANGE, BR_INACTIVE, BR_NOTECPLOT
    character(len=10) :: CCTIME
    character(len=12) :: CDATE
    character(len=72) :: RSIFN
    character(len=180) :: MODDIR ! CURRENT WORKING DIRECTORY
    real(R8), save, allocatable, dimension(:,:) :: RATZ, CURZ1, CURZ2, CURZ3 ! SW 5/15/06
    real(R8) :: G = 9.81D0, PI = 3.14159265359D0
    real(R8) :: DENSITY
    external :: DENSITY
end module GLOBAL

module GEOMC
    use PREC
    integer, allocatable, dimension(:) :: JBUH, JBDH, JWUH, JWDH
    real(R8), allocatable, dimension(:) :: ALPHA, SINA, COSA, SLOPE, BKT, DLX, DLXR, SLOPEC, SINAC
    real(R8), allocatable, dimension(:,:) :: H, H1, H2, BH1, BH2, BHR1, BHR2, AVHR
    real(R8), allocatable, dimension(:,:) :: B, BI, BB, BH, BHR, BR, EL, AVH1, AVH2, BNEW ! SW 1/23/06
    real(R8), allocatable, dimension(:,:) :: DEPTHB, DEPTHM, FETCHU, FETCHD
    real(R8), allocatable, dimension(:) :: Z, ELWS
    real(R8), allocatable, dimension(:) :: BCONSTRICTION
    logical, allocatable, dimension(:,:) :: CONSTRICTION
end module GEOMC

module NAMESC
    integer, allocatable, dimension(:) :: LNAME
    character(len=6), allocatable, dimension(:) :: CUNIT, CUNIT2
    character(len=8), allocatable, dimension(:) :: CNAME2, CDNAME2
    character(len=9), allocatable, dimension(:) :: FMTH, FMTC, FMTCD
    character(len=19), allocatable, dimension(:) :: CNAME1
    character(len=43), allocatable, dimension(:) :: CNAME, CNAME3, CDNAME, CDNAME3, HNAME
    character(len=72), allocatable, dimension(:) :: TITLE
    character(len=10), allocatable, dimension(:,:) :: CONV
end module NAMESC

module STRUCTURES
    use PREC
    real(R8) :: DIA, FMAN, CLEN, CLOSS, UPIE, DNIE
    real(R8), allocatable, dimension(:) :: QOLD, QOLDS, VMAX, DTP, DTPS
    real(R8), allocatable, dimension(:) :: EGT, A1GT, B1GT, G1GT, A2GT, B2GT, G2GT, EGT2
    real(R8), allocatable, dimension(:) :: QGT, GTA1, GTB1, GTA2, GTB2, BGT
    real(R8), allocatable, dimension(:) :: QSP, A1SP, B1SP, A2SP, B2SP, ESP
    real(R8), allocatable, dimension(:) :: EUPI, EDPI, WPI, DLXPI, FPI, FMINPI, QPI, BP
    real(R8), allocatable, dimension(:,:) :: YS, VS, YSS, VSS, YST, VST, YSTS, VSTS
    integer, allocatable, dimension(:) :: IUPI, IDPI, JWUPI, JWDPI, JBDPI, JBUPI
    integer, allocatable, dimension(:) :: IUSP, IDSP, JWUSP, JWDSP, JBUSP, JBDSP
    integer, allocatable, dimension(:) :: IUGT, IDGT, JWUGT, JWDGT, JBUGT, JBDGT
    integer, allocatable, dimension(:) :: IWR, KTWR, KBWR
    real, allocatable, dimension(:) :: EKTWR, EKBWR ! SW 3/18/16
    logical, allocatable, dimension(:) :: LATERAL_SPILLWAY, LATERAL_PIPE, LATERAL_GATE, LATERAL_PUMP, BEGIN, WLFLAG, PUMP_DOWNSTREAM
    character(len=8), allocatable, dimension(:) :: LATGTC, LATSPC, LATPIC, LATPUC, DYNGTC, DYNPIPE, DYNPUMP ! SW 5/10/10
    character(len=8) :: GT2CHAR
    real(R8), allocatable, dimension(:) :: EPU, STRTPU, ENDPU, EONPU, EOFFPU, QPU
    integer, allocatable, dimension(:) :: IUPU, IDPU, KTPU, KBPU, JWUPU, JWDPU, JBUPU, JBDPU
    real(R8) :: THR = 0.01D0, OMEGA = 0.8D0, EPS2 = 0.0001D0
    integer :: NN = 19, NNPIPE = 19, NC = 7
    real, allocatable, dimension(:) :: EGTO, BGTO
    character(len=8), allocatable, dimension(:) :: GTIC
    logical, allocatable, dimension(:,:) :: ACTIVE_RULE_W2SELECTIVE
end module STRUCTURES

module TRANS
    use PREC
    real(R8), allocatable, dimension(:) :: THETA
    real(R8), pointer, dimension(:,:) :: COLD, CNEW, SSB, SSK
    real(R8), allocatable, dimension(:,:) :: DX, DZ, DZQ
    real(R8), allocatable, dimension(:,:) :: ADX, ADZ, AT, VT, CT, DT
end module TRANS

module SURFHE
    use PREC
    real(R8) :: RHOWCP !, PHISET
    real(R8), allocatable, dimension(:) :: ET, CSHE, LAT, LONGIT, SHADE, RB, RE, RC
    real(R8), allocatable, dimension(:) :: WIND, WINDH, WSC, AFW, BFW, CFW, PHI0
    logical, allocatable, dimension(:) :: RH_EVAP
    integer, allocatable, dimension(:) :: IWIND !MLM 08/12/05
end module SURFHE

module TVDC
    use PREC
    real(R8), allocatable, dimension(:) :: QIN, QTR, QDTR, PR, ELUH, ELDH, QWD, QSUM
    real(R8), allocatable, dimension(:) :: TIN, TTR, TDTR, TPR, TOUT, TWDO, TIND, QIND
    real(R8), allocatable, dimension(:) :: TAIR, TDEW, CLOUD, PHI, SRON, PALT_JW, ELWS_INI ! systdg PALT_JW, ELWS_INI
    real(R8), allocatable, dimension(:,:) :: CIN, CTR, CDTR, CPR, CIND, TUH, TDH, QOUT
    real(R8), allocatable, dimension(:,:,:) :: CUH, CDH
    integer :: NAC, NOPEN
    integer, allocatable, dimension(:) :: NACPR, NACIN, NACDT, NACTR, NACD, CN
    integer, allocatable, dimension(:,:) :: TRCN, INCN, DTCN, PRCN
    logical :: CONSTITUENTS
    character(len=72) :: QGTFN, QWDFN, WSCFN, SHDFN
    character(len=72), allocatable, dimension(:) :: METFN, QOTFN, QINFN, TINFN, CINFN, QTRFN, TTRFN, CTRFN, QDTFN
    character(len=72), allocatable, dimension(:) :: TDTFN, CDTFN, PREFN, TPRFN, CPRFN, EUHFN, TUHFN, CUHFN, EDHFN
    character(len=72), allocatable, dimension(:) :: EXTFN, CDHFN, TDHFN
    integer :: WAIT_TIME, TIME_BUFFER !SR 11/26/19
    logical, allocatable, dimension(:) :: WAIT_FOR_TRIB_INPUT, WAIT_FOR_BRANCH_INPUT !SR 11/26/19
    character(len=240), allocatable, dimension(:) :: TR_FILEDIR, BR_FILEDIR
end module TVDC

module KINETIC
    use PREC
    real :: KDO, PCO2, PCO2ATMPPM ! SW 8/16/2020
    real(R8) :: O2CH4, O2H2S, O2FE2, O2MN2
    real(R8), allocatable, dimension(:) :: CoeffA_Turb, CoeffB_Turb, SECC_PAR
    real(R8), allocatable, dimension(:) :: H2SQ10, H2S1DK, CH4Q10, CH41DK
    real(R8), allocatable, dimension(:) :: CH4R, H2SR, FEIIR, MNIIR, SO4R
    real(R8), allocatable, dimension(:) :: KFE_OXID, KFE_RED, KFEOOH_HalfSat, FeSetVel, KMN_OXID, KMN_RED, KMNO2_HalfSat, MnSetVel
    real(R8), allocatable, dimension(:) :: BACTQ10, BACT1DK, BACTLDK, BACTS !A_DISG, B_DISG, C_DISG
    real(R8), pointer, dimension(:,:) :: TDS, COL, NH4, NO3, PO4, DSI, PSI, LDOM
    real(R8), pointer, dimension(:,:) :: N2, H2S, CH4, SO4, FEII, FEOOH, MNII, MNO2
    real(R8), pointer, dimension(:,:) :: WAGE, BACT, DGP
    real(R8), pointer, dimension(:,:) :: RDOM, LPOM, RPOM, O2, TIC, ALK
    real(R8), pointer, dimension(:,:) :: COLSS, NH4SS, NO3SS, PO4SS, FESS, DSISS, PSISS, LDOMSS
    real(R8), pointer, dimension(:,:) :: N2SS, H2SSS, CH4SS, SO4SS, FEIISS, FEOOHSS, MNIISS, MNO2SS
    real(R8), pointer, dimension(:,:) :: AGESS, BACTSS, DISGSS
    real(R8), pointer, dimension(:,:) :: RDOMSS, LPOMSS, RPOMSS, DOSS, TICSS, CASS
    real(R8), pointer, dimension(:,:) :: ALKSS ! enhanced pH buffering
    real, pointer, dimension(:,:) :: PH, CO2, HCO3, CO3
    real, pointer, dimension(:,:) :: TN, TP, TKN
    real, pointer, dimension(:,:) :: DON, DOP, DOC, NH3
    real, pointer, dimension(:,:) :: PON, POP, POC
    real, pointer, dimension(:,:) :: TON, TOP, TOC
    real, pointer, dimension(:,:) :: APR, CHLA, ATOT
    real, pointer, dimension(:,:) :: O2DG, TDG, TURB, SECCHID
    real, pointer, dimension(:,:) :: SSSI, SSSO, TISS, TOTSS
    real, pointer, dimension(:,:) :: PO4AR, PO4AG, PO4AP, PO4SD, PO4SR, PO4NS, PO4POM, PO4DOM, PO4OM
    real, pointer, dimension(:,:) :: CH4SR, H2SSR, FEIISR, MNIISR
    real, pointer, dimension(:,:) :: PO4ER, PO4EG, PO4EP, TICEP, DOEP, DOER
    real, pointer, dimension(:,:) :: NH4ER, NH4EG, NH4EP, NO3EG, DSIEG, LDOMEP, LPOMEP
    real, pointer, dimension(:,:) :: NH4AR, NH4AG, NH4AP, NH4SD, NH4SR, NH4D, NH4POM, NH4DOM, NH4OM, NH3GAS
    real, pointer, dimension(:,:) :: NO3AG, NO3D, NO3SED
    real, pointer, dimension(:,:) :: DSIAG, DSID, DSISD, DSISR, DSIS
    real, pointer, dimension(:,:) :: PSIAM, PSID, PSINS
    real, pointer, dimension(:,:) :: FENS, FESR
    real, pointer, dimension(:,:) :: LDOMAP, LDOMD, LRDOMD, RDOMD
    real, pointer, dimension(:,:) :: LPOMAP, LPOMD, LRPOMD, RPOMD, LPOMNS, RPOMNS
    real, pointer, dimension(:,:) :: DOAP, DOAR, DODOM, DOPOM, DOOM, DONIT
    real, pointer, dimension(:,:) :: DOSED, DOSOD, DOBOD, DOAE
    real, pointer, dimension(:,:) :: CBODU, CBODDK, TICAP
    real, pointer, dimension(:,:) :: SEDD, SODD, SEDAS, SEDOMS, SEDNS
    real, pointer, dimension(:,:) :: SEDD1, SEDD2
    real(R8), pointer, dimension(:,:,:) :: SS, ALG, CBOD, CG
    real(R8), pointer, dimension(:,:,:) :: SSSS, ASS, CBODSS, CGSS
    real, pointer, dimension(:,:,:) :: AGR, ARR, AER, AMR, ASR
    real, pointer, dimension(:,:,:) :: EGR, ERR, EER, EMR, EBR
    real(R8), pointer, dimension(:,:) :: LDOMP, RDOMP, LPOMP, RPOMP, LDOMN, RDOMN, LPOMN, RPOMN
    real(R8), pointer, dimension(:,:) :: LDOMPSS, RDOMPSS, LPOMPSS, RPOMPSS, LDOMNSS, RDOMNSS
    real(R8), pointer, dimension(:,:) :: LPOMNSS, RPOMNSS
    real, pointer, dimension(:,:) :: LDOMPAP, LDOMPEP, LPOMPAP, LPOMPNS, RPOMPNS
    real, pointer, dimension(:,:) :: LDOMNAP, LDOMNEP, LPOMNAP, LPOMNNS, RPOMNNS
    real, pointer, dimension(:,:) :: SEDDP, SEDASP, SEDOMSP, SEDNSP, LPOMEPP
    real, pointer, dimension(:,:) :: SEDDN, SEDASN, SEDOMSN, SEDNSN, LPOMEPN, SEDNO3
    real, pointer, dimension(:,:) :: SEDDC, SEDASC, SEDOMSC, SEDNSC, LPOMEPC
    real, pointer, dimension(:,:) :: CBODNS, SEDCB, SEDCBP, SEDCBN, SEDCBC
    real, pointer, dimension(:,:) :: SEDBR, SEDBRP, SEDBRC, SEDBRN, CO2REAER !CB 11/30/06
    real(R8), pointer, dimension(:,:,:) :: CBODP, CBODN ! CB 6/6/10
    real(R8), pointer, dimension(:,:,:) :: CBODPSS, CBODNSS ! CB 6/6/10
    real, pointer, dimension(:,:) :: CBODNSP, CBODNSN ! cb 6/6/10
    real, pointer, dimension(:,:) :: DOH2S, DOCH4, DOSEDIA, DOFE2
    real, pointer, dimension(:,:) :: H2SREAER, CH4REAER, CH4D, H2SD, SDINC, SDINP, SDINN, H2SSEDD
    real, pointer, dimension(:,:) :: FE2D, SDINFEOOH, SDINMNO2, MN2D, DOMN2
    real, allocatable, dimension(:,:,:) :: EPM, EPD, EPC
    real, allocatable, dimension(:) :: CGQ10, CG0DK, CG1DK, CGS, CGLDK, CGKLF, CGCS, CGR !LCJ 2/26/15 SW 10/16/15
    real, allocatable, dimension(:) :: SOD, SDK, LPOMDK, RPOMDK, LDOMDK, RDOMDK, LRDDK, LRPDK
    real, allocatable, dimension(:) :: SDK1, sdk2 ! amaila
    real, allocatable, dimension(:) :: SSS, TAUCR, POMS, seds, sedb, SSCS !  , SSFLOC   cb 11/27/06   SR 04/21/13
    real, allocatable, dimension(:) :: AG, AR, AE, AM, AS, AHSN, AHSP, AHSSI, ASAT
    real, allocatable, dimension(:) :: AP, AN, AC, ASI, ACHLA, APOM, ANPR
    real, allocatable, dimension(:) :: EG, ER, EE, EM, EB
    real, allocatable, dimension(:) :: EHSN, EHSP, EHSSI, ESAT, EHS, ENPR
    real, allocatable, dimension(:) :: EP, EN, EC, ESI, ECHLA, EPOM
    real(R8), allocatable, dimension(:) :: BETA, EXH2O, EXSS, EXOM, EXA
    real, allocatable, dimension(:) :: DSIR, PSIS, PSIDK, PARTSI
    real, allocatable, dimension(:) :: ORGP, ORGN, ORGC, ORGSI
    real, allocatable, dimension(:) :: pbiom, nbiom, cbiom ! Amaila, cb 6/8/17
    real, allocatable, dimension(:) :: BODP, BODN, BODC
    real, allocatable, dimension(:) :: PO4R, PARTP
    real, allocatable, dimension(:) :: NH4DK, NH4R, NO3DK, NO3S, FNO3SED
    real, allocatable, dimension(:) :: O2AG, O2AR, O2OM, O2NH4
    real, allocatable, dimension(:) :: O2EG, O2ER
    real, allocatable, dimension(:) :: CO2R
    real, allocatable, dimension(:) :: KBOD, TBOD, RBOD
    real, allocatable, dimension(:) :: CAQ10, CADK, CAS
    real, allocatable, dimension(:) :: OMT1, OMT2, SODT1, SODT2, NH4T1, NH4T2, NO3T1, NO3T2
    real, allocatable, dimension(:) :: OMK1, OMK2, SODK1, SODK2, NH4K1, NH4K2, NO3K1, NO3K2, KG_H2O_CONSTANT
    real, allocatable, dimension(:) :: AT1, AT2, AT3, AT4
    real, allocatable, dimension(:) :: AK1, AK2, AK3, AK4
    real, allocatable, dimension(:) :: ET1, ET2, ET3, ET4
    real, allocatable, dimension(:) :: EK1, EK2, EK3, EK4
    real(R8), allocatable, dimension(:) :: WIND10, CZ, QC, QERR
    real, allocatable, dimension(:) :: REAER, RCOEF1, RCOEF2, RCOEF3, RCOEF4, DGPO2, MINKL
    real, allocatable, dimension(:,:) :: DO1, DO2, DO3, GAMMA, F_NH3
    real, allocatable, dimension(:,:) :: SED, FPSS, FPFE, FE
    real, allocatable, dimension(:,:) :: SED1, sed2, SED1ic, sed2ic ! cb 6/17/17
    real, allocatable, dimension(:,:,:) :: CBODD
    real, allocatable, dimension(:) :: CBODS, PFLUXIN, NFLUXIN
    real, allocatable, dimension(:,:) :: ORGPLD, ORGPRD, ORGPLP, ORGPRP, ORGNLD, ORGNRD, ORGNLP, ORGNRP
    real, allocatable, dimension(:,:) :: LDOMPMP, LDOMNMP, LPOMPMP, LPOMNMP, RPOMPMP, RPOMNMP
    real, allocatable, dimension(:,:) :: LPZOOINP, LPZOOINN, LPZOOOUTP, LPZOOOUTN
    real, allocatable, dimension(:,:) :: SEDC, SEDN, SEDP, SEDNINFLUX, SEDPINFLUX ! SW 4/2016
    real, allocatable, dimension(:,:) :: SEDVPC, SEDVPP, SEDVPN
    real, allocatable, dimension(:,:) :: SDKV, SEDDKTOT
    integer :: NLDOMP, NRDOMP, NLPOMP, NRPOMP, NLDOMN, NRDOMN, NLPOMN, NRPOMN
    integer, allocatable, dimension(:) :: NAF, NEQN, ANEQN, ENEQN !, FLOCEQN   ! SR 04/21/13
    integer, allocatable, dimension(:,:) :: KFCN
    logical, allocatable, dimension(:) :: SEDIMENT_RESUSPENSION
    character(len=8), allocatable, dimension(:) :: CAC, REAERC, AVERTM
    character(len=10), allocatable, dimension(:,:) :: LFPR
! enhanced pH buffering start
    character(len=8) :: nh4bufc, po4bufc, ombufc, omtype, pombufc, phbufc, ncalkc, CO2YEARLYPPM
    integer :: nag, nagi
    logical :: ammonia_buffering, phosphate_buffering, om_buffering, pom_buffering, pH_buffering, NONCON_ALKALINITY, ALGAE_SETTLING_EXIST
    logical, allocatable, dimension(:,:) :: sdfirstadd
    real, allocatable, dimension(:) :: sdeni, pki, pksd, sden, pk, fract
! enhanced pH buffering end
!
    real(R8), pointer, dimension(:,:) :: LDOMC, RDOMC, LPOMC, RPOMC ! W2V3.8 NEW STATE VARIABLES
    real(R8), pointer, dimension(:,:) :: LDOMCSS, RDOMCSS, LPOMCSS, RPOMCSS ! W2V3.8 NEW STATE VARIABLES SOURCE AND SINK
    real, allocatable, dimension(:,:) :: PSIEM, SEDEB ! W2V3.8 NEW FLUX
    real, allocatable, dimension(:,:) :: LPOMPEP, LPOMNEP, LPOMCEP
    real, allocatable, dimension(:,:) :: LPOMHD, RPOMHD
    real, allocatable, dimension(:,:) :: LDOMCAP, LDOMCEP, LPOMCAP, LPOMCNS, RPOMCNS
    real, allocatable, dimension(:,:) :: LDOMPD, LRDOMPD, RDOMPD, LPOMPD, LRPOMPD, RPOMPD, LPOMPHD, RPOMPHD
    real, allocatable, dimension(:,:) :: LDOMND, LRDOMND, RDOMND, LPOMND, LRPOMND, RPOMND, LPOMNHD, RPOMNHD
    real, allocatable, dimension(:,:) :: LDOMCD, LRDOMCD, RDOMCD, LPOMCD, LRPOMCD, RPOMCD, LPOMCHD, RPOMCHD
    real, allocatable, dimension(:) :: LPOMHK, RPOMHK
    real, allocatable, dimension(:) :: LDOMPDK, LRDOMPDK, RDOMPDK, LDOMNDK, LRDOMNDK, RDOMNDK, LDOMCDK, LRDOMCDK, RDOMCDK
    real, allocatable, dimension(:) :: LPOMPDK, LRPOMPDK, RPOMPDK, LPOMNDK, LRPOMNDK, RPOMNDK, LPOMCDK, LRPOMCDK, RPOMCDK
    real, allocatable, dimension(:,:) :: LDOMCMP, LPOMCMP, RPOMCMP
    real, allocatable, dimension(:,:) :: LPZOOINC, LPZOOOUTC
    real, allocatable, dimension(:,:) :: LDOP, RDOP, LPOP, RPOP, LDON, RDON, LPON, RPON, LDOC, RDOC, LPOC, RPOC
    integer :: NLDOMC, NRDOMC, NLPOMC, NRPOMC
!

contains

    real function SATO(T, SAL, P, SALT_WATER)
        real(R8) :: T, SAL
        real(R8) :: P
        logical :: SALT_WATER
        SATO = EXP(7.7117 - 1.31403*LOG(T + 45.93))*P
        if (SALT_WATER) then
            SATO = EXP(LOG(SATO) - SAL*(1.7674E-2 - 1.0754E1/(T + 273.15) + 2.1407E3/(T + 273.15)**2)) ! SAL is in ppt
        else
            if (SAL > 100.) then
                SATO = EXP(LOG(SATO) - SAL/1000.*(1.7674E-2 - 1.0754E1/(T + 273.15) + 2.1407E3/(T + 273.15)**2)) ! SAL is in mg/l
            end if
        end if
    end function SATO

    real function FR(TT, TT1, TT2, SK1, SK2)
        real(R8) :: TT
        real :: TT1, TT2, SK1, SK2
        FR = SK1*EXP(LOG(SK2*(1.0 - SK1)/(SK1*(1.0 - SK2)))/(TT2 - TT1)*(TT - TT1))
    end function FR

    real function FF(TT, TT3, TT4, SK3, SK4)
        real :: TT3, TT4, SK3, SK4
        real(R8) :: TT
        FF = SK4*EXP(LOG(SK3*(1.0 - SK4)/(SK4*(1.0 - SK3)))/(TT4 - TT3)*(TT4 - TT))
    end function FF

end module KINETIC

module SELWC
    use PREC
    real(R8), allocatable, dimension(:) :: VNORM, QNEW
!REAL,                  ALLOCATABLE, DIMENSION(:)    :: EWD, TAVGW                   ! cb 1/16/13
    real, allocatable, dimension(:) :: ewd
    real(R8), allocatable, dimension(:,:) :: QSTR, QSW
!REAL,                  ALLOCATABLE, DIMENSION(:,:)  :: ESTR,   WSTR, TAVG            ! SW Selective 7/30/09
    real, allocatable, dimension(:,:) :: estr, wstr ! cb 1/16/13
    real(r8), allocatable, dimension(:,:) :: tavg ! cb 1/16/13
    real(r8), allocatable, dimension(:) :: tavgw ! cb 1/16/13
    real(R8), allocatable, dimension(:,:) :: CAVGW, CDAVGW
    real(R8), allocatable, dimension(:,:,:) :: CAVG, CDAVG
    integer, allocatable, dimension(:) :: NSTR, NOUT, KTWD, KBWD, KTW, KBW
    integer, allocatable, dimension(:,:) :: KTSW, KBSW, KOUT
    character(len=8), allocatable, dimension(:) :: DYNSTRUC
end module SELWC

module GDAYC
    real :: DAYM, EQTNEW
    integer :: JDAYG, IMON, YEAR, GDAY, YEAROLD
    logical :: LEAP_YEAR
    character(len=9) :: MONTH
end module GDAYC

module SCREENC
    use PREC
    real :: JDAY, DLTS1, JDMIN, MINDLT, DLTAV, ELTMJD
    real, allocatable, dimension(:) :: ZMIN, CMIN, CMAX, HYMIN, HYMAX, CDMIN, CDMAX
    integer :: ILOC, KLOC, IMIN, KMIN, NIT, NV, JTT, JWW
    integer, allocatable, dimension(:) :: IZMIN
    character(len=8), allocatable, dimension(:) :: ACPRC, AHPRC, ACDPRC
end module SCREENC

module TRIDIAG_V
    use PREC
!INTEGER,                             INTENT(IN)  :: S, E, N
!  REAL(R8),              DIMENSION(:), INTENT(IN)  :: A(E),V(E),C(E),D(E)
!  REAL(R8),              DIMENSION(:), INTENT(OUT) :: U(N)
    real(R8), allocatable, dimension(:) :: BTA1, GMA1
!REAL(R8), DIMENSION(1000)              :: BTA, GMA
!  INTEGER                                          :: I
end module TRIDIAG_V

module TDGAS
    real, allocatable, dimension(:) :: AGASSP, BGASSP, CGASSP, AGASGT, BGASGT, CGASGT
    integer, allocatable, dimension(:) :: EQSP, EQGT
end module TDGAS

module LOGICC
    logical :: SUSP_SOLIDS, OXYGEN_DEMAND, UPDATE_GRAPH, INITIALIZE_GRAPH
    logical :: WITHDRAWALS, TRIBUTARIES, GATES, PIPES
    logical, allocatable, dimension(:) :: NO_WIND, NO_INFLOW, NO_OUTFLOW, NO_HEAT
    logical, allocatable, dimension(:) :: UPWIND, ULTIMATE, FRESH_WATER, SALT_WATER
    logical, allocatable, dimension(:) :: LIMITING_DLT, TERM_BY_TERM, MANNINGS_N, PH_CALC
    logical, allocatable, dimension(:) :: ONE_LAYER, DIST_TRIBS, PRECIPITATION
    logical, allocatable, dimension(:) :: PRINT_SEDIMENT, LIMITING_FACTOR, READ_EXTINCTION, READ_RADIATION
    logical, allocatable, dimension(:) :: PRINT_SEDIMENT1, PRINT_SEDIMENT2 ! Amaila
    logical, allocatable, dimension(:) :: UH_INTERNAL, DH_INTERNAL, UH_EXTERNAL, DH_EXTERNAL
    logical, allocatable, dimension(:) :: UQ_INTERNAL, DQ_INTERNAL, UQ_EXTERNAL, DQ_EXTERNAL
    logical, allocatable, dimension(:) :: UP_FLOW, DN_FLOW, INTERNAL_FLOW
    logical, allocatable, dimension(:) :: DAM_INFLOW, DAM_OUTFLOW !TC 08/03/04
    logical, allocatable, dimension(:) :: INTERP_METEOROLOGY, INTERP_INFLOW, INTERP_DTRIBS, INTERP_TRIBS
    logical, allocatable, dimension(:) :: INTERP_WITHDRAWAL, INTERP_HEAD, INTERP_EXTINCTION
    logical, allocatable, dimension(:) :: VISCOSITY_LIMIT, CELERITY_LIMIT, IMPLICIT_AZ, TRAPEZOIDAL !SW 07/16/04
!  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: HYDRO_PLOT,         CONSTITUENT_PLOT, DERIVED_PLOT
    logical, allocatable, dimension(:) :: INTERP_GATE ! cb 8/13/2010
    logical, allocatable, dimension(:,:) :: PRINT_DERIVED, PRINT_HYDRO, PRINT_CONST, PRINT_EPIPHYTON
    logical, allocatable, dimension(:,:) :: POINT_SINK, INTERNAL_WEIR, INTERP_OUTFLOW
end module LOGICC

module SHADEC
    integer, parameter :: IANG = 18
    real, parameter :: GAMA = 3.1415926*2./REAL(IANG) ! SW 10/17/05
    real, dimension(IANG) :: ANG ! SW 10/17/05
    real, allocatable, dimension(:) :: A00, DECL, HH, TTLB, TTRB, CLLB, CLRB ! SW 10/17/05
    real, allocatable, dimension(:) :: SRLB1, SRRB1, SRLB2, SRRB2, SRFJD1, SRFJD2, SHADEI
    real, allocatable, dimension(:,:) :: TOPO
    logical, allocatable, dimension(:) :: DYNAMIC_SHADE
end module SHADEC

module EDDY
    use PREC
    character(len=8), allocatable, dimension(:) :: AZC, IMPTKE
    real(R8), allocatable, dimension(:) :: WSHY, FRIC
    real(R8), allocatable, dimension(:,:) :: FRICBR, DECAY
    real(R8), allocatable, dimension(:,:,:) :: TKE
    real(R8), allocatable, dimension(:,:) :: AZT, DZT
    real(R8), allocatable, dimension(:) :: USTARBTKE, E
    real(R8), allocatable, dimension(:) :: EROUGH, ARODI, TKELATPRDCONST, STRICK
    integer, allocatable, dimension(:) :: FIRSTI, LASTI, WALLPNT, TKEBC
    logical, allocatable, dimension(:) :: STRICKON, TKELATPRD
end module EDDY

module MACROPHYTEC
    real, pointer, dimension(:,:) :: NH4MR, NH4MG, LDOMMAC, RPOMMAC, LPOMMAC, DOMP, DOMR, TICMC
    real, pointer, dimension(:,:) :: PO4MR, PO4MG
    real, allocatable, dimension(:) :: MG, MR, MM, MMAX, MBMP
    real, allocatable, dimension(:) :: MT1, MT2, MT3, MT4, MK1, MK2, MK3, MK4
    real, allocatable, dimension(:) :: MP, MN, MC
    real, allocatable, dimension(:) :: PSED, NSED, MHSP, MHSN, MHSC, MSAT, EXM
    real, allocatable, dimension(:) :: CDDRAG, DWV, DWSA, ANORM
    real, allocatable, dimension(:) :: ARMAC
    real, allocatable, dimension(:) :: O2MG, O2MR, LRPMAC, MPOM
    real, allocatable, dimension(:,:) :: MACMBRS, MACMBRT, SSMACMB
    real, allocatable, dimension(:,:) :: CW, BIC, MACWBCI
    real, allocatable, dimension(:,:,:) :: MACTRMR, MACTRMF, MACTRM
    real, allocatable, dimension(:,:,:) :: MMR, MRR
    real, allocatable, dimension(:,:,:) :: MAC, MACT
    real, allocatable, dimension(:,:,:) :: MPLIM, MNLIM, MCLIM
    real, allocatable, dimension(:,:,:) :: SMAC, SMACT
    real, allocatable, dimension(:,:,:) :: GAMMAJ
    real, allocatable, dimension(:,:,:,:) :: MGR
    real, allocatable, dimension(:,:,:,:) :: MACRC, MACRM
    real, allocatable, dimension(:,:,:,:) :: MLLIM
    real, allocatable, dimension(:,:,:,:) :: MACSS
    real, allocatable, dimension(:,:,:,:) :: SMACRC, SMACRM
    logical, allocatable, dimension(:) :: KTICOL
    logical, allocatable, dimension(:,:) :: PRINT_MACROPHYTE, MACROPHYTE_CALC
    logical :: MACROPHYTE_ON
    character(len=8), allocatable, dimension(:,:) :: MPRWBC, MACWBC ! cb 8/24/15
    character(len=10), allocatable, dimension(:,:) :: CONV2
    character(len=10), allocatable, dimension(:,:,:,:) :: MLFPR
end module MACROPHYTEC

module POROSITYC
    real, allocatable, dimension(:) :: SAREA, VOLKTI
    real, allocatable, dimension(:,:) :: POR, VOLI, VSTEMKT
    real, allocatable, dimension(:,:,:) :: VSTEM
    logical, allocatable, dimension(:) :: HEAD_FLOW
    logical, allocatable, dimension(:) :: UP_HEAD
end module POROSITYC

module ZOOPLANKTONC
    use PREC
    logical :: ZOOPLANKTON_CALC
    real, allocatable, dimension(:) :: ZG, ZM, ZEFF, PREFP, ZR, ZOOMIN, ZS2P, EXZ, ZS ! SW 1/29/2019
    real, allocatable, dimension(:) :: ZT1, ZT2, ZT3, ZT4, ZK1, ZK2, ZK3, ZK4
    real, allocatable, dimension(:) :: ZP, ZN, ZC, O2ZR
    real, allocatable, dimension(:,:) :: PREFA, PREFZ ! OMNIVOROUS ZOOPLANKTON
    real, allocatable, dimension(:,:) :: PO4ZR, NH4ZR, DOZR, TICZR, LPZOOOUT, LPZOOIN
    real(R8), pointer, dimension(:,:,:) :: ZOO, ZOOSS
    real, allocatable, dimension(:,:,:) :: ZMU, TGRAZE, ZRT, ZMT, ZSR
    real, allocatable, dimension(:,:,:) :: ZOORM, ZOORMR, ZOORMF
    real, allocatable, dimension(:,:,:) :: AGZT
    real, allocatable, dimension(:,:,:,:) :: AGZ, ZGZ ! OMNIVOROUS ZOOPLANKTON
end module ZOOPLANKTONC


module INITIALVELOCITY
    use PREC
    logical :: INIT_VEL, ONCE_THROUGH
    real(R8), allocatable, dimension(:) :: QSSI, ELWSS, UAVG
    real(R8), allocatable, dimension(:,:) :: BSAVE
    logical, allocatable, dimension(:) :: LOOP_BRANCH
end module INITIALVELOCITY


module ENVIRPMOD
    character(len=3), save, allocatable, dimension(:) :: CC_E, CD_E
    character(len=3), save :: VEL_VPR, TEMP_VPR, SELECTIVEC, DEPTH_VPR
    real, save :: TEMP_INT, TEMP_TOP, VEL_INT, VEL_TOP, TIMLAST, DLTT, TEMP_C
    real, save :: D_INT, D_TOP, SJDAY1, SJDAY2
    real, save :: T_SUM, V_SUM, D_SUM, D_AVG, D_C, V_AVG, T_AVG, T_CRIT, V_CRIT, C_CRIT, CD_CRIT, D_CRIT
    real, save, allocatable, dimension(:) :: T_CNT, T_TOT, SUMVOLT, V_CNT, V_TOT, VOLGL, VEL_C, D_CNT, D_TOT, C_SUM, CD_SUM, C_AVG, CD_AVG
    real, allocatable, save, dimension(:,:) :: C_CNT, CD_CNT, C_TOT, CD_TOT, T_CLASS, V_CLASS, D_CLASS
    real, allocatable, save, dimension(:,:,:) :: C_CLASS, CD_CLASS
    real, allocatable, save, dimension(:,:) :: CONC_C, CONC_CD
    real, allocatable, save, dimension(:) :: C_INT, C_TOP, CD_INT, CD_TOP, CN_E, CDN_E
    integer :: CONE = 1500, NUMCLASS, IOPENFISH, NAC_E, NACD_E, JJ, JACD, I_SEGINT
    integer :: ISTART(9), IEND(9)
end module ENVIRPMOD


module ALGAE_TOXINS
    use PREC
    integer, parameter :: NUMATOXINS = 4
    integer :: NATS, NATE, ATOXIN_DEBUG_FN = 2501
    logical :: ALGAE_TOXIN !,ALGAE_TOXIN_FILE
    real, allocatable, dimension(:,:) :: CTP, CTB
    real, dimension(:) :: CTREL(NUMATOXINS), CTD(NUMATOXINS)
    real(R8), pointer, dimension(:,:,:) :: EX_TOXIN, CTESS
    real(R8), allocatable, dimension(:,:,:) :: IN_TOXIN
    character(len=2) :: ATOX, ATOX_DEBUG
end module ALGAE_TOXINS


module MAIN
    use PREC
! Variable declaration
    real(8), allocatable, dimension(:) :: IceQSS
!INTEGER       :: J,NIW,NGC,NGCS,NTDS,NCCS,NGCE,NSSS,NSSE,NPO4,NNH4
    integer :: J, NIW, NGC, NTDS, NCCS, NSSS, NSSE, NPO4, NNH4 ! CEMA -placed NGCS and NGCE in global module  SW 10/16/2015
    integer :: NN2, NH2S, NCH4, NSO4, NFEII, NFEOOH, NMNII, NMNO2, NMFT
    integer :: NWAGE, NBACT, NDGP
    integer :: NNO3, NDSI, NPSI, NFE, NLDOM, NRDOM, NLPOM, NRPOM, NBODS
    integer :: KF_DO_SED, KF_DO_SOD, KF_SED_PBURIAL, KF_SED_NBURIAL, KF_NH4_SD, KF_NH4_SR, KF_PO4_SD, KF_PO4_SR, KF_CO2X
    integer :: KF_DOH2S, KF_SDINC, KF_DOCH4, KF_FE2D, KF_MN2D, KF_SEDD
    integer :: KF_NO3D, KF_NO3AG, KF_NO3EG, KF_NO3SED, KF_NH3GAS
    integer :: DOC_DER, POC_DER, TOC_DER, DON_DER, PON_DER, TON_DER, TKN_DER, TN_DER, NH3_DER, DOP_DER, POP_DER, TOP_DER, TP_DER, APR_DER
    integer :: CHLA_DER, ATOT_DER, O2DG_DER, TDG_DER, TURB_DER, TOTSS_DER, TISS_DER, CBODU_DER, PH_DER, CO2_DER, HCO3_DER, CO3_DER, SECCHI_DER
    integer :: NBODE, NAS, NAE, NDO, NTIC, NALK, NTRT, NWDT
    integer :: NDT, JS, JP, JG, JT, JH, NTSR, NIWDO, NRSO, JD
    integer :: NBODCS, NBODCE, NBODPS, NBODPE, NBODNS, NBODNE, IBOD, JCB ! cb 6/6/10
    integer :: JF, JA, JM, JE, JJZ, K, L3, L1, L2, NTAC, NDSP, NTACMX, NTACMN, JFILE, M
    integer :: KBP, JWR, JJJB, JDAYNX, NIT1, JWD, L, IUT, IDT, KTIP
    integer :: INCRIS, IE, II, NDLT, NRS, INCR, IS, JAC
    real :: JDAYTS, JDAY1, TMSTRT, TMEND, HMAX, DLXMAX, CELRTY, NXTVD, TTIME
    real(R8) :: DLMR, TICE ! SW 4/19/10
    real(R8) :: TAU1, TAU2, ELTMS, HMIN, DLXMIN, RHOICP
    real(R8) :: WLF, FLOWBALF, NPBALF, NXWL, NXFLOWBAL, NXNPBAL
    real(R8) :: HRAD
    real(R8) :: ZB, WWT, DFC, GC2, UDR, UDL, AB, EFFRIC
    real(R8) :: DEPKTI, COLB, COLDEP, SSTOT, RHOIN, VQIN, VQINI
    real(R8) :: QINFR, ELT, RHOIRL1, V1, BHSUM, BHRSUM, WT1, WT2
    real(R8) :: ICETHU, ICETH1, ICETH2, ICE_TOL = 0.0050D0, DEL, HICE ! SW 4/19/10
    real(R8) :: DLTCAL, HEATEX, SROOUT, SROSED, SROIN, SRONET, TFLUX, HIA
    real(R8) :: TAIRV, EA, ES, DTV
    real(R8) :: T2R4
    integer :: CON = 10, RSI, GRF, NDG = 16, FLOWBFN = 9500, WLFN = 9510, AERATEFN = 9520, FISHHABFN = 9530, MASSBFN = 9501, LAKE_RIVER_CONTOUR = 9540
    integer :: VSF, SIF, JG_AGE ! SR 7/27/2017
    integer :: JSG, NNSG ! 1/16/13
    logical :: ADD_LAYER, SUB_LAYER, WARNING_OPEN, ERROR_OPEN, VOLUME_WARNING, SURFACE_WARNING
    logical :: END_RUN, BRANCH_FOUND, NEW_PAGE, UPDATE_KINETICS, UPDATE_RATES
    logical :: WEIR_CALC, DERIVED_CALC, RESTART_IN, RESTART_OUT
    logical :: SPILLWAY, PUMPS, MWB_EXIST ! SW 12/2/2019
    logical :: TIME_SERIES, DOWNSTREAM_OUTFLOW, ICE_COMPUTATION
    logical :: DSI_CALC, PO4_CALC, N_CALC ! cb 10/12/11
    logical :: TDGON, GAS_TRANSFER_UPDATE
    logical :: FISH_PARTICLE_EXIST, WAIT_FOR_INFLOW_RESULTS ! SW 4/30/15, 2/9/2019
    character(len=1) :: ESC
    character(len=2) :: DEG
    character(len=3) :: GDCH
    character(len=8) :: WLC, FLOWBALC, NPBALC, SED_DIAG
    character(len=8) :: RSOC, RSIC, CCC, LIMC, WDOC, TSRC, EXT, SELECTC, CLOSEC, HABTATC, ENVIRPC, AERATEC, INITUWL, DLTINTER ! SW 7/31/09; 8/24/09
    character(len=10) :: BLANK = "          ", BLANK1 = "    -99.00", SEDCH, SEDPCH, SEDNCH, SEDCCH
    character(len=72) :: WDOFN, RSOFN, TSRFN, SEGNUM, LINE, SEGNUM2, TSRFN1
    logical :: RETLOG, STANDING_BIOMASS_DECAY, PHBUFF_EXIST, WATER_AGE_ACTIVE ! SW 5/26/15  SR 7/27/2017
    logical :: DYNPIPEADJUST ! SW 2/18/2020
    character(len=2) :: DYNPAD
    integer :: DYNPAD_SEG, DYNPAD_PIPE, DYNPIPELOG = 9505

    real :: DYNPAD_WL, DYNPAD_MAXRATE, DYNPAD_PERCENTCHANGE
    integer :: N_WAITS, NWAIT !SR 11/26/19
    integer, allocatable, dimension(:) :: WAIT_INDEX !SR 11/26/19
    character(len=2), allocatable, dimension(:) :: WAIT_TYPE !SR 11/26/19
    character(len=240), allocatable, dimension(:) :: FILEDIR

    logical :: ORGC_CALC

! Allocatable array declarations

    real, allocatable, dimension(:) :: ETUGT, EBUGT, ETDGT, EBDGT
    real, allocatable, dimension(:) :: ETUSP, EBUSP, ETDSP, EBDSP
    real, allocatable, dimension(:) :: ETUPI, EBUPI, ETDPI, EBDPI, ETPU, EBPU, TSEDF
    real, allocatable, dimension(:) :: CSUM, CDSUM, X1
    real, allocatable, dimension(:) :: RSOD, RSOF, DLTD, DLTF, DLTMAX
    real(R8), allocatable, dimension(:) :: QWDO
    real(R8), allocatable, dimension(:) :: ICETHI, ALBEDO, HWI, BETAI, GAMMAI, ICEMIN, ICET2, CBHE, TSED
    real(R8), allocatable, dimension(:) :: FI, SEDCI, FSOD, FSED, AX, RANLW, T2I, ELBOT, DXI
    real(R8), allocatable, dimension(:) :: SEDCI1, SEDCI2, fsedc1, fsedc2 ! cb 6/7/17, Amaila
    real(R8), allocatable, dimension(:) :: WSHX ! QINT,   QOUTT,
    real(R8), allocatable, dimension(:) :: SROSH, EV, RS, RN
    real(R8), allocatable, dimension(:) :: QDT, QPR, ICESW
    real(R8), allocatable, dimension(:) :: XBR, QPRBR, EVBR, TPB
    real(R8), allocatable, dimension(:) :: DLXRHO, Q, QSSUM
    real, allocatable, dimension(:) :: ELTRT, ELTRB
    real, allocatable, dimension(:) :: TSRD, TSRF, WDOD, WDOF
    real, allocatable, dimension(:) :: QOAVR, QIMXR, QOMXR, QTAVB, QTMXB
    real(R8), allocatable, dimension(:) :: FETCH, ETSR
    real(R8), allocatable, dimension(:) :: QINSUM, TINSUM
    real, allocatable, dimension(:) :: CDTOT
    real, allocatable, dimension(:) :: SEDCIP, SEDCIN, SEDCIC, SEDCIS
    real, allocatable, dimension(:,:) :: ESTRT, WSTRT, CINSUM, HAB
    real(R8), allocatable, dimension(:,:) :: P, HSEG, QTOT
    real, allocatable, dimension(:,:) :: CPB, COUT, CWDO, CDWDO, KFJW
    real(R8), allocatable, dimension(:,:) :: C2I, EPICI
    real(R8), allocatable, dimension(:,:) :: QTRF
    real, allocatable, dimension(:,:) :: SNPD, SCRD, PRFD, SPRD, CPLD, VPLD, FLXD
    real, allocatable, dimension(:,:) :: SNPF, SCRF, PRFF, SPRF, CPLF, VPLF, FLXF
    real, allocatable, dimension(:,:) :: TVP, SEDVP, QINF
    real, allocatable, dimension(:,:) :: SEDVP1, SEDVP2 ! Amaila
    real(R8), allocatable, dimension(:,:) :: TSSUH1, TSSDH1, ATM_DEP_LOADING
    real(R8), allocatable, dimension(:,:,:) :: CSSUH1, CSSDH1
    real, allocatable, dimension(:,:,:) :: EPIVP, CVP, macrcvp, macrclp ! cb 8/21/15
    real(R8), allocatable, dimension(:) :: VOLB
    real(R8), allocatable, dimension(:) :: DLVOL, VOLG
    real(R8), allocatable, dimension(:) :: A, C, D, F, V, BTA, GMA, BHRHO
    real(R8), allocatable, dimension(:) :: DLVR, ESR, ETR
    real, allocatable, dimension(:,:) :: CMBRS
    integer, allocatable, dimension(:) :: KTUGT, KBUGT, KTDGT, KBDGT
    integer, allocatable, dimension(:) :: KTUSP, KBUSP, KTDSP, KBDSP
    integer, allocatable, dimension(:) :: KTUPI, KBUPI, KTDPI, KBDPI
    integer, allocatable, dimension(:) :: NSNP, NSCR, NSPR, NVPL, NFLX, NCPL, BTH
    integer, allocatable, dimension(:) :: VPR, LPR, NIPRF, NISPR, NPRF
    integer, allocatable, dimension(:) :: NISNP
    integer, allocatable, dimension(:) :: NBL, KBMAX, KBI
    integer, allocatable, dimension(:) :: KBR, IBPR
    integer, allocatable, dimension(:) :: TSR
    integer, allocatable, dimension(:) :: NPOINT, NL, KTQIN, KBQIN, ilayer ! SW 1/23/06
    integer, allocatable, dimension(:) :: ITR, KTTR, KBTR, JBTR
    integer, allocatable, dimension(:) :: IWD, KWD, JBWD
    integer, allocatable, dimension(:) :: IWDO, ITSR, JBTSR
    integer, allocatable, dimension(:) :: ILAT, JBDAM, JSS
    integer, allocatable, dimension(:) :: ICPL, NACATD
    integer, allocatable, dimension(:,:) :: KTSWT, KBSWT, ATMDCN
    integer, allocatable, dimension(:,:) :: IPRF, ISPR, ISNP, BL, WDO, CDN, WDO2
    logical, allocatable, dimension(:) :: ALLOW_ICE, PUMPON, FETCH_CALC, ATM_DEPOSITION ! ICE_IN,     RC/SW 4/28/11
    logical, allocatable, dimension(:) :: DN_HEAD, HEAD_BOUNDARY
    logical, allocatable, dimension(:) :: PLACE_QIN, PLACE_QTR, SPECIFY_QTR
    logical, allocatable, dimension(:) :: ISO_TEMP, VERT_TEMP, LONG_TEMP, VERT_PROFILE, LONG_PROFILE
    logical, allocatable, dimension(:) :: SEDIMENT_CALC, DETAILED_ICE, IMPLICIT_VISC, SNAPSHOT, PROFILE
    logical, allocatable, dimension(:) :: SEDIMENT_CALC1, SEDIMENT_CALC2 ! Amaila
    logical, allocatable, dimension(:) :: VECTOR, CONTOUR, SPREADSHEET, SCREEN_OUTPUT
    logical, allocatable, dimension(:) :: FLUX, EVAPORATION, ZERO_SLOPE
    logical, allocatable, dimension(:) :: ISO_SEDIMENT, VERT_SEDIMENT, LONG_SEDIMENT
    logical, allocatable, dimension(:) :: ISO_SEDIMENT1, VERT_SEDIMENT1, LONG_SEDIMENT1, ISO_SEDIMENT2, VERT_SEDIMENT2, LONG_SEDIMENT2 ! Amaila
    logical, allocatable, dimension(:) :: VOLUME_BALANCE, ENERGY_BALANCE, MASS_BALANCE, BOD_CALC, ALG_CALC
    logical, allocatable, dimension(:) :: BOD_CALCP, BOD_CALCN ! cb 5/19/2011
    logical, allocatable, dimension(:,:) :: ISO_EPIPHYTON, VERT_EPIPHYTON, LONG_EPIPHYTON, EPIPHYTON_CALC
    logical, allocatable, dimension(:,:) :: ISO_CONC, VERT_CONC, LONG_CONC, TDG_SPILLWAY, TDG_GATE
    logical, allocatable, dimension(:,:) :: iso_macrophyte, vert_macrophyte, long_macrophyte ! cb 8/21/15
    character(len=4), allocatable, dimension(:) :: CUNIT1
    character(len=8), allocatable, dimension(:) :: SEG, SEDRC, TECPLOT
    character(len=8), allocatable, dimension(:) :: HPLTC, CPLTC, CDPLTC
    character(len=8), allocatable, dimension(:) :: EXC, EXIC
    character(len=8), allocatable, dimension(:) :: GASGTC, GASSPC, ATM_DEPOSITIONC, ATM_DEPOSITION_INTERPOLATION
    character(len=10), allocatable, dimension(:) :: CWDOC, CDWDOC
    character(len=8), allocatable, dimension(:) :: ICEC, SEDCc, SEDPRC, SNPC, SCRC, SPRC, PRFC, DYNSEDK
    character(len=8), allocatable, dimension(:) :: SEDCc1, SEDPRC1, SEDCc2, SEDPRC2 ! Amaila
    character(len=8), allocatable, dimension(:) :: RHEVC, VPLC, CPLC, AZSLC, FETCHC
    character(len=8), allocatable, dimension(:) :: DTRC, SROC, KFAC, CDAC
    character(len=8), allocatable, dimension(:) :: INCAC, TRCAC, DTCAC, PRCAC
    character(len=8), allocatable, dimension(:) :: WTYPEC, GRIDC !SW 07/16/04
    character(len=8), allocatable, dimension(:) :: PUSPC, PDSPC, PUGTC, PDGTC, PDPIC, PUPIC, PPUC, TRC
    character(len=8), allocatable, dimension(:) :: SLICEC, FLXC
    character(len=8), allocatable, dimension(:) :: VBC, MBC, EBC
    character(len=8), allocatable, dimension(:) :: PQC, EVC, PRC
    character(len=8), allocatable, dimension(:) :: QINC, QOUTC, WINDC, HEATC
    character(len=8), allocatable, dimension(:) :: VISC, CELC, DLTADD
    character(len=8), allocatable, dimension(:) :: SLTRC, SLHTC, FRICC
    character(len=8), allocatable, dimension(:) :: QINIC, TRIC, DTRIC, WDIC, HDIC, METIC !, KFNAME2
    character(len=14), allocatable, dimension(:) :: KFNAME2
    character(len=10), allocatable, dimension(:) :: C2CH, CDCH, EPCH, macch, KFCH, APCH, ANCH, ALCH ! SW 10/20/15
    character(len=45), allocatable, dimension(:) :: KFNAME
    character(len=72), allocatable, dimension(:) :: SNPFN, PRFFN, VPLFN, CPLFN, SPRFN, FLXFN, FLXFN2, BTHFN, VPRFN, LPRFN, SPRVFN, ATMDEPFN ! SW 9/28/2018
    character(len=8), allocatable, dimension(:,:) :: SINKC, SINKCT
    character(len=8), allocatable, dimension(:,:) :: CPRBRC, CDTBRC, CPRWBC, CINBRC, CTRTRC, HPRWBC, STRIC, CDWBC, KFWBC
    character(len=8), allocatable, dimension(:,:) :: EPIC, EPIPRC, C_ATM_DEPOSITION
    character(len=12), allocatable, dimension(:,:) :: CONV1
    character(len=72) :: CONFN = "w2_con.npt"
    character(len=72) :: TEXT
    integer, allocatable, dimension(:) :: WBSEG !  systdg
    integer :: CONTDG ! systdg
    integer :: TARGETFNNO ! systdg TDGtarget
    real(R8), allocatable, dimension(:) :: GTPC ! systdg
    logical :: SYSTDG, N2BND, DOBND, DGPBND ! systdg
    logical :: TDGTA ! systdg TDGtarget
    character(len=8) :: SYSTDGC, N2BNDC, DOBNDC, TDG2BNDC ! systdg
    character(len=8) :: TDGTAC ! systdg TDGtarget
    character(len=8), allocatable, dimension(:) :: GTTYP ! systdg
!  CHARACTER(72), PARAMETER              :: CONTDGFN='w2_systdg.npt'   ! systdg
! Data declarations
    real(R8) :: RK1 = 2.12D0, RL1 = 333507.0D0, RIMT = 0.0D0, RHOA = 1.25D0, RHOI = 916.0D0, VTOL = 1.0D3, CP = 4186.0D0, THRKTI = 0.10D0
    logical :: LAKE_RIVER_CONTOURC ! SW 2/27/2020
    character(len=2) :: LAKE_RIVER_CONTOUR_ON
    character(len=80) :: FILE_LAKE_CONTOUR_T(20), FILE_LAKE_CONTOUR_DO(20), FILE_RIVER_CONTOUR_T(20), FILE_RIVER_CONTOUR_DO(20)
    integer :: NUM_LAKE_CONTOUR, NUM_RIVER_CONTOUR, LAKE_CONTOUR_SEG(20), RIVER_CONTOUR_BR1(20), RIVER_CONTOUR_BR2(20)
    integer :: JW_RIVER_CONTOUR(20), JW_LAKE_CONTOUR(20), LAKE_CONTOUR_FORMAT, RIVER_CONTOUR_FORMAT
    real :: LAKE_CONTOUR_START(20), LAKE_CONTOUR_FREQ(20), RIVER_CONTOUR_START(20), RIVER_CONTOUR_FREQ(20), NXT_RIVER_CONTOUR(20), NXT_LAKE_CONTOUR(20)

!DATA CON   /10/                      !,  RSI /11/
!DATA FLOWBFN /9500/, WLFN /9510/, AERATEFN /9520/, FISHHABFN /9530/, MASSBFN /9501/        ! SW 5/25/15 NOTE THAT FISHHABFN INCREMENTS 3 TIMES SO 9531,9532,9533 ARE RESERVED
!DATA DYNPIPELOG /9505/, LAKE_RIVER_CONTOUR /9540/  ! ALSO 9540 TO 9550 to 9560 to 9570 to 9580 ARE TAKEN
! for fish habitat output filenames
    integer :: JBFILE1 = 9581, JWFILE1 = 9620 !ALSO JWFILE1=9620+NWB  JBFILE1=9581+NBR


end module MAIN

module BIOENERGETICS
    integer, allocatable, dimension(:) :: IBIO, BIODP, BIOEXPFN, WEIGHTNUM
    real, allocatable, dimension(:) :: BIOD, BIOF, VOLROOS
    real, allocatable, dimension(:,:) :: C2W
    real, allocatable, dimension(:,:,:) :: C2ZOO
    character(len=8) :: BIOC
    logical :: BIOEXP, FISHBIO
    integer :: NBIO, NIBIO, KLIM, FISHBIOFN = 9502
    real(8) :: NXBIO, NXTBIO, GAMMAB
    character(len=72) :: BIOFN, WEIGHTFN
    character(len=8) :: BHEAD(20)

end module BIOENERGETICS


module CEMAVars
    use PREC
    integer(4), allocatable, dimension(:) :: ConsolidationType, ConstPoreWtrRate, NumCEMAPWInst
    integer(4), allocatable, dimension(:) :: ConsRegSegSt, ConsRegSegEn, ConsolidRegnNum
    character(len=256) :: ConsolidRateRegnFil
    real(R8), allocatable, dimension(:) :: ConsolidRateTemp
    real(R8), allocatable, dimension(:) :: BedElevation, BedElevationLayer, BedPorosity
    real(R8), allocatable, dimension(:,:) :: CellArea
    real(R8), allocatable, dimension(:) :: BedConsolidRate, PorewaterRelRate, ConstConsolidRate
    real(R8), allocatable, dimension(:) :: CEMACumPWRelease, CEMACumPWReleaseRate, CEMACumPWToRelease, CEMACumPWReleased
    real(R8), allocatable, dimension(:,:) :: CEMASedConc, BubbleRelWB ! SW 7/1/2017
    real(R8), allocatable, dimension(:) :: VOLCEMA
    logical, allocatable, dimension(:) :: CEMALayerAdded, CEMASSApplied
    logical, allocatable, dimension(:) :: EndBedConsolidation, BedConsolidationSeg ! cb 6/28/17
    logical, allocatable, dimension(:) :: ApplyCEMAPWRelease
    logical :: DYNAMIC_SD

    real(R8), allocatable, dimension(:) :: SDRegnPOC_T, SDRegnPON_T, SDRegnPOP_T, SDRegnSul_T
    real(R8), allocatable, dimension(:) :: SDRegnPOC_L_Fr, SDRegnPOC_R_Fr, SDRegnPON_L_Fr
    real(R8), allocatable, dimension(:) :: SDRegnPON_R_Fr, SDRegnPW_DiffCoeff, SDRegnOx_Threshold
    real(R8), allocatable, dimension(:) :: SDRegnPOP_L_Fr, SDRegnPOP_R_Fr
    real(R8), allocatable, dimension(:) :: SDRegnAe_NH3_NO3_L, SDRegnAe_NH3_NO3_H, SDRegnAe_NO3_N2_L
    real(R8), allocatable, dimension(:) :: SDRegnAe_NO3_N2_H, SDRegnAn_NO3_N2, SDRegnAe_CH4_CO2
    real(R8), allocatable, dimension(:) :: SDRegnAe_HS_NH4_Nit, SDRegnAe_HS_O2_Nit, SDRegn_Theta_PW, SDRegn_Theta_PM
    real(R8), allocatable, dimension(:) :: SDRegn_Theta_NH3_NO3, SDRegn_Theta_NO3_N2, SDRegn_Theta_CH4_CO2
    real(R8), allocatable, dimension(:) :: SDRegn_Sulfate_CH4_H2S, SDRegnAe_H2S_SO4, SDRegn_Theta_H2S_SO4
    real(R8), allocatable, dimension(:) :: SDRegn_NormConst_H2S_SO4, SDRegn_MinRate_PON_Lab, SDRegn_MinRate_PON_Ref
    real(R8), allocatable, dimension(:) :: SDRegn_MinRate_PON_Ine, SDRegn_MinRate_POC_Lab, SDRegn_MinRate_POC_Ref
    real(R8), allocatable, dimension(:) :: SDRegn_MinRate_POC_Ine, SDRegn_Theta_PON_Lab, SDRegn_Theta_PON_Ref
    real(R8), allocatable, dimension(:) :: SDRegn_Theta_PON_Ine, SDRegn_Theta_POC_Lab, SDRegn_Theta_POC_Ref
    real(R8), allocatable, dimension(:) :: SDRegn_Theta_POC_Ine ! cb 10/8/13    !Real(8), Allocatable, Dimension(:) :: SDRegn_Theta_POC_Ine, SDRegn_CH4CompMethod
    real(R8), allocatable, dimension(:) :: SDRegn_Theta_POP_Lab, SDRegn_Theta_POP_Ref, SDRegn_Theta_POP_Ine
    real(R8), allocatable, dimension(:) :: Kdp1, Kdp2, KdFe1, KdFe2, KdMn1, KdMn2, KdNH31, KdNH32, KdH2S1, KdH2S2
    real(R8), allocatable, dimension(:) :: PartMixVel, BurialVel, POCr, delta_kpo41, DOcr, KsOxch
    real(R8), allocatable, dimension(:) :: SDRegn_MinRate_POP_Lab, SDRegn_MinRate_POP_Ref, SDRegn_MinRate_POP_Ine
    real(R8), allocatable, dimension(:) :: SD_NO3p2, SD_NH3p2, SD_NH3Tp2, SD_CH4p2, SD_PO4p2, SD_PO4Tp2
    real(R8), allocatable, dimension(:) :: SD_HSp2, SD_HSTp2
    real(R8), allocatable, dimension(:) :: SD_poc2, SD_pon2, SD_pop2, SD_NH3Tp, SD_NO3p, SD_PO4Tp, SD_HSTp
    real(R8), allocatable, dimension(:) :: SD_fpon, SD_fpoc, SD_kdiaPON, SD_ThtaPON, SD_kdiaPOC, SD_ThtaPOC
    real(R8), allocatable, dimension(:) :: SD_kdiaPOP, SD_ThtaPOP, SD_NH3T, SD_PO4, SD_FPOP
    real(R8), allocatable, dimension(:) :: SD_JPOC, SD_JPON, SD_JPOP, SD_TDS

    real(R8), allocatable, dimension(:) :: SD_Denit, SD_JDenit, SD_JO2NO3, SD_HS ! cb 7/26/18
    real(R8), allocatable, dimension(:) :: SD_Fe2
    real(R8), allocatable, dimension(:) :: SD_Mn2
    real(R8), allocatable, dimension(:) :: SD_pHValue ! cb 7/26/18    !Real(8), Allocatable, Dimension(:) :: SD_SO4Conc, SD_pHValue
    real(R8), allocatable, dimension(:) :: SD_EPOC, SD_EPON, SD_EPOP
    real(R8), allocatable, dimension(:) :: SD_AerLayerThick

    real(R8), allocatable, dimension(:,:,:) :: MFTSedFlxVars, CEMA_SD_Vars
    real(R8), allocatable, dimension(:,:) :: CEMATSSCopy

    integer(4), allocatable, dimension(:) :: CEMAMFT_RandC_RegN, CEMAMFT_InCond_RegN
    integer(4), allocatable, dimension(:) :: SedBedInitRegSegSt, SedBedInitRegSegEn
    integer(4), allocatable, dimension(:) :: SedBedDiaRCRegSegSt, SedBedDiaRCRegSegEn

    integer(4), allocatable, dimension(:) :: FFTActPrdSt, FFTActPrdEn
    integer, allocatable, dimension(:) :: SDRegn_CH4CompMethod, SDRegn_POMResuspMethod
    real(R8), allocatable, dimension(:) :: FFTLayConc

    real(R8), allocatable, dimension(:) :: H2SDis, H2SGas, CH4Dis, CH4Gas, NH4Dis, NH4Gas, CO2Dis, CO2Gas
    real(R8), allocatable, dimension(:) :: BubbleRadiusSed, PresBubbSed, PresCritSed
    real(R8), allocatable, dimension(:) :: CgSed, C0Sed, CtSed
    real(R8), allocatable, dimension(:,:,:) :: TConc, TConcP, SConc
    real(R8), allocatable, dimension(:,:,:) :: DissolvedGasSediments
    integer(8), allocatable, dimension(:) :: MFTBubbReleased, LastDiffVolume
    integer(4), allocatable, dimension(:,:) :: BubblesLNumber, BubblesStatus
    real(R8), allocatable, dimension(:,:) :: BubblesRadius, BubblesRiseV, BubblesCarried
    real(R8), allocatable, dimension(:,:,:) :: BubblesGasConc, BRVoluAGas, BRRateAGas
    real(R8), allocatable, dimension(:,:) :: BubblesReleaseAllValue, BRRateAGasNet
    real(R8), allocatable, dimension(:) :: BottomTurbulence
    logical, allocatable, dimension(:) :: CrackOpen
    logical, allocatable, dimension(:,:) :: FirstBubblesRelease, BubblesAtSurface


    integer(4) :: CEMAFilN, NumConsolidRegns, CEMASedimentType
    integer(4) :: CEMASNPOutFilN = 2411, CEMATSR1OutFilN = 2412, SegNumI, LayerNum
    integer(4) :: CEMABtmLayFilN = 2414, TempCntr1
    integer(4) :: CEMASedFlxFilN1 = 2415, CEMASedFlxFilN2 = 2416, CEMASedFlxFilN3 = 2417, CEMALogFilN = 2418, CEMASedFlxFilN4 = 2419
    integer(4) :: CEMASedFlxFilN5 = 2420, CEMASedFlxFilN6 = 2421, CEMASedFlxFilN7 = 2422, CEMASedFlxFilN8 = 2423, CEMASedFlxFilN9 = 2424, CEMASedFlxFilN10 = 2441, CEMASedFlxFilN11 = 2442
    integer(4) :: CEMASedFlxFilN12 = 3483, CEMASedFlxFilN13 = 3484, CEMASedFlxFilN14 = 3485, CEMASedFlxFilN15 = 3486, CEMASedFlxFilN16 = 3487, CEMASedFlxFilN17 = 3488, CEMASedFlxFilN18 = 3489
    integer(4) :: CEMASedFlxFilN19 = 3490, CEMASedFlxFilN20 = 3491, CEMASedFlxFilN21 = 3492, CEMASedFlxFilN22 = 3493
    integer(4) :: CEMASedFlxFilN23 = 3494, CEMASedFlxFilN24 = 3495, CEMASedFlxFilN25 = 3496, CEMASedFlxFilN26 = 3497, CEMASedFlxFilN27 = 3498, CEMASedFlxFilN28 = 3499
    integer(4) :: CEMASedFlxFilN29 = 3500, CEMASedFlxFilN30 = 3501, CEMASedFlxFilN31 = 3502
    integer(4) :: CEMASedFlxFilN32 = 3503, CEMASedFlxFilN33 = 3504, CEMASedFlxFilN34 = 3505
    integer(4) :: CEMASedFlxFilN35 = 3506, CEMASedFlxFilN36 = 3507, CEMASedFlxFilN37 = 3508
    integer(4) :: CEMAOutFilN1 = 2426, CEMAOutFilN2 = 2429, CEMAOutFilN3 = 2435, CEMAOutFilN4 = 2437, CEMAOutFilBub = 2440
    integer(4) :: CEMAOutFilN5 = 2438, CEMAOutFilN6 = 2439
    integer(4) :: NumRegnsSedimentDiagenesis, NumRegnsSedimentBedComposition
    integer(4) :: NumFFTActivePrds, FFTActPrd, NumBubRelArr, NumGas = 4
    real(R8) :: LayerAddThkFrac, BedElevationInit, BedPorosityInit, CEMASedimentDensity
    real(R8) :: CEMAParticleSize, TotalPoreWatVolume, TotalSedimentsInBed, TotalPoreWatRemoved
    real(R8) :: CEMASedimentSVelocity, CEMAPWpH
    real(R8) :: NH4_NH3_Eqb_Const, HS_H2S_Eqb_Const, VolumeIncreasedConsolid
    real(R8) :: GasConst_R = 0.0821 !L.atm/mol/K
    real(R8) :: HenryConst_NH3, HenryConst_CH4, HenryConst_H2S, HenryConst_CO2
    real(R8) :: InitFFTLayerConc, FFTLayerSettVel
    real(R8) :: GasDiff_Sed, CalibParam_R1, YoungModulus, CritStressIF
    real(R8) :: BubbRelScale, CrackCloseFraction, MaxBubbRad, BubbRelFraction, BubbAccFraction
    real(R8) :: BubbRelFractionAtm, BubbWatGasExchRate
    real(R8) :: CEMATurbulenceScaling
    real(R8) :: IceThicknessChange ! cb 2/5/13
    real(R8) :: TAUCRPOM, crshields, spgrav_POM, dia_POM, GasReleaseCH4 !,GasReleaseCO2   ! SW 10/10/2017   ! SW 10/19/2017
    real(R8) :: NXTSEDIAG, SEDIAGFREQ ! SW 5/25/2017

    logical :: CEMARelatedCode, IncludeBedConsolidation, IncludeCEMASedDiagenesis, IncludeFFTLayer, FFTActive, FirstTimeInFFTCode
    logical :: IncludeIron, IncludeManganese, IncludeDynamicpH, IncludeAlkalinity, SD_global ! cb 5/22/15
    logical :: CEMASedimentProcessesInc, WriteBESnp, WritePWSnp, WriteCEMAMFTSedFlx, CEMA_POM_Resuspension
    logical :: FirstTimeinCEMAMFTSedDiag, MoveFFTLayerDown
    logical :: LimBubbSize, UseReleaseFraction, FirstTimeInBubbles, ApplyBubbTurb
    logical :: sediment_diagenesis, cao_method, Bubbles_CalculatioN

    real(R8), allocatable, dimension(:) :: SDRegnH2S_T, SDRegnNH3_T, SDRegnCH4_T, SDRegnNO3_T
    real(R8), allocatable, dimension(:) :: SDRegnTIC_T, SDRegnALK_T, SDRegnPO4_T
    real(R8), allocatable, dimension(:) :: SDRegnFe2_T, SDRegnFeOOH_T, SDRegnMn2_T, SDRegnMnO2_T
    real(R8), allocatable, dimension(:) :: SDRegnT_T, SDRegnpH
    real(R8), allocatable, dimension(:) :: SDPFLUX, SDNH4FLUX, SDNO3FLUX

end module CEMAVars


module Selective1TDGtarget
    character(len=72) :: TITLETDGTARGET(10)
    integer :: NGFL, NGPH, NGSP, NOUTS, NGSPPH
    integer :: tsiteration, dygroup
    integer, allocatable, dimension(:) :: SPGTNO, PHGTNO, SPPRIOR
    real :: tsfreq, tsconv, NXTSPLIT, DAYTEST, NXTSPLIT2
    real :: tstsrt, tstend, tstarget, tstarget2, nxtjday
    real, allocatable, dimension(:) :: SPMINFRAC, PHMAXFLOW
    character(len=8) :: tsyearly, tsdynsel, dyupdate
end module Selective1TDGtarget
