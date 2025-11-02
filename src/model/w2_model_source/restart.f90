subroutine RESTART_OUTPUT(RSOFN)
    use GLOBAL;     use SCREENC;     use RSTART;     use GDAYC;     use GEOMC;     use KINETIC, only: EPM, EPD, SEDC, SEDN, SEDP, PH, SDKV;     use TVDC, only: QSUM
    use KINETIC, only: SED, PFLUXIN, NFLUXIN;     use ZOOPLANKTONC, only: ZOO;     use EDDY, only: TKE;     use MAIN, only: ENVIRPC, WARNING_OPEN, NXWL, NXFLOWBAL, NXNPBAL;     use ENVIRPMOD;     use LOGICC;     use STRUCTURES
    use MACROPHYTEC;     use CEMAVars;     use CEMASedimentDiagenesis, only: C2SF
    implicit none

    character(len=*) :: RSOFN
    open(RSO, FILE=RSOFN, FORM="UNFORMATTED", STATUS="UNKNOWN")
    write(RSO) NIT, NV, KMIN, IMIN, NSPRF, CMBRT, ZMIN, IZMIN, START, CURRENT
    write(RSO) DLTDP, SNPDP, TSRDP, VPLDP, PRFDP, CPLDP, SPRDP, RSODP, SCRDP, FLXDP, WDODP
    write(RSO) JDAY, ELTM, ELTMF, DLT, DLTAV, DLTS, MINDLT, JDMIN, CURMAX
    write(RSO) NXTMSN, NXTMTS, NXTMPR, NXTMCP, NXTMVP, NXTMRS, NXTMSC, NXTMSP, NXTMFL, NXTMWD, NXWL, NXFLOWBAL, NXNPBAL, NXTMWD_SEC
    write(RSO) VOLIN, VOLOUT, VOLUH, VOLDH, VOLPR, VOLTRB, VOLDT, VOLWD, VOLEV, VOLSBR, VOLTR, VOLSR, VOLICE, ICEBANK
    write(RSO) TSSEV, TSSPR, TSSTR, TSSDT, TSSWD, TSSIN, TSSOUT, TSSS, TSSB, TSSICE
    write(RSO) TSSUH, TSSDH, TSSUH2, TSSDH2, CSSUH2, CSSDH2, VOLUH2, VOLDH2, QUH1
    write(RSO) ESBR, ETBR, EBRI
    write(RSO) Z, SZ, ELWS, SAVH2, SAVHR, H2
    write(RSO) KTWB, KTI, SKTI, SBKT
    write(RSO) ICE, ICETH, CUF, QSUM
    write(RSO) U, W, SU, SW, AZ, SAZ, DLTLIM
    write(RSO) T1, T2, C1, C2, C1S, SED, KFS, CSSK
    write(RSO) EPD, EPM
    write(RSO) MACMBRT, MACRC, SMACRC, MAC, SMAC, MACRM, MACSS
    write(RSO) SEDC, SEDN, SEDP, ZOO, CD ! mlm 10/06
    write(RSO) SDKV ! MLM 6/10/07
    write(RSO) TKE ! SW 10/4/07
    write(RSO) BR_INACTIVE, WARNING_OPEN ! SW 8/1/2018
    if (ENVIRPC == "      ON") then
        write(RSO) T_CLASS, V_CLASS, C_CLASS, CD_CLASS, T_TOT, T_CNT, SUMVOLT, V_CNT, V_TOT, C_TOT, C_CNT, CD_TOT, CD_CNT
    end if
    if (PIPES) then
        write(RSO) YS, VS, VST, YST, DTP, QOLD
    end if
    write(RSO) TPOUT, TPTRIB, TPDTRIB, TPWD, TPPR, TPIN, TP_SEDSOD_PO4, PFLUXIN, TNOUT, TNTRIB, TNDTRIB, TNWD, TNPR, TNIN, TN_SEDSOD_NH4, NFLUXIN, ATMDEP_P, ATMDEP_N, NH3GASLOSS !TP_SEDBURIAL,TN_SEDBURIAL,
!IF(SEDIMENT_DIAGENESIS)THEN
!WRITE(RSO)MFTSedFlxVars,CellArea,TConc,SConc,CrackOpen,BubbleRelWB,BedPorosity,MFTBubbReleased,GasReleaseCH4,SedGenBODConc   
!ENDIF
    if (IncludeCEMASedDiagenesis) then
        write(RSO) C2SF, CellArea, BedPorosity
        if (Bubbles_Calculation) then
            write(RSO) TConc, SConc, CrackOpen, BubbleRelWB, MFTBubbReleased, GasReleaseCH4
        end if
    end if

    close(RSO)
end subroutine RESTART_OUTPUT
