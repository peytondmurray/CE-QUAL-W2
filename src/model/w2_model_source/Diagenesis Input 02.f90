subroutine CEMA_W2_Input()
    use MAIN
    use GLOBAL
    use KINETIC
    use GEOMC
    use CEMAVars
    use SCREENC, only: JDAY

! Type declarations
    implicit none

    logical :: SkipLoop !, file_exists
    character(len=256) :: MessageTemp
    character(len=20) :: ADUMMY ! SW 2/2019

    integer :: monzz, dayzz, yearzz, ninp, JSKIP

    SD_global = .false.
    IncludeIron = .false.
    IncludeManganese = .false.
    IncludeDynamicpH = .false.
    IncludeAlkalinity = .false.
    Bubbles_Calculation = .false.

!INQUIRE(FILE="W2_diagenesis.npt", EXIST=file_exists)   ! file_exists will be TRUE if the file
    if (SED_DIAG /= "      ON") then
        CEMARelatedCode = .false.
        IncludeBedConsolidation = .false.
        return
    end if
    CEMARelatedCode = .true.

    CEMAFilN = NUNIT;     NUNIT = NUNIT + 1 ! SW 7/8/2019
    open(CEMAFilN, File="W2_diagenesis.npt", STATUS="OLD")
    open(CEMALogFilN, File="DiagenesisLogFile.opt", STATUS="UNKNOWN")

!Read Header
    SkipLoop = .false.
    do while (.not. SkipLoop)
        read(CEMAFilN, "(a)") MessageTemp
        if (index(MessageTemp, "$") == 0) then
            SkipLoop = .true.
        end if
    end do
    backspace(CEMAFilN)
!
! GROUP 1: Global Control
    read(CEMAFilN, *) MessageTemp, SD_global
    if (.not. SD_global) then
        CEMARelatedCode = .false.
        IncludeFFTLayer = .false.
        IncludeBedConsolidation = .false.
        IncludeCEMASedDiagenesis = .false.
        return
    end if

! GROUP 2: FFT Layer
    read(CEMAFilN, *) MessageTemp, IncludeFFTLayer
    if (IncludeFFTLayer) then
        FirstTimeInFFTCode = .true.
        read(CEMAFilN, *) MessageTemp, NumFFTActivePrds
        allocate(FFTActPrdSt(NumFFTActivePrds), FFTActPrdEn(NumFFTActivePrds))
        allocate(FFTLayConc(IMX))
        read(CEMAFilN, *) MessageTemp, (FFTActPrdSt(i), i = 1, NumFFTActivePrds)
        read(CEMAFilN, *) MessageTemp, (FFTActPrdEn(i), i = 1, NumFFTActivePrds)
        read(CEMAFilN, *) MessageTemp, InitFFTLayerConc
        read(CEMAFilN, *) MessageTemp, FFTLayerSettVel
        FFTLayerSettVel = FFTLayerSettVel/DAY !m/d --> m/s
        FFTLayConc = 0.d00
        FFTActPrd = 1
        MoveFFTLayerDown = .false.
        read(CEMAFilN, *) MessageTemp, MoveFFTLayerDown
    else
        FirstTimeInFFTCode = .false.
        MoveFFTLayerDown = .false.
        do JSKIP = 1, 6
            read(CEMAFilN, *)
        end do
    end if
!
! GROUP 3: Bed Consolidation
    read(CEMAFilN, *) MessageTemp, IncludeBedConsolidation
    if (IncludeBedConsolidation) then
        read(CEMAFilN, *) MessageTemp, LayerAddThkFrac
        read(CEMAFilN, *) MessageTemp, NumConsolidRegns
        allocate(ConsolidationType(NumConsolidRegns), ConstConsolidRate(NumConsolidRegns))
        allocate(ConstPoreWtrRate(NumConsolidRegns), ConsolidRateTemp(NumConsolidRegns))
        allocate(ConsRegSegSt(NumConsolidRegns), ConsRegSegEn(NumConsolidRegns))
        read(CEMAFilN, *) MessageTemp, (ConsRegSegSt(i), i = 1, NumConsolidRegns)
        read(CEMAFilN, *) MessageTemp, (ConsRegSegEn(i), i = 1, NumConsolidRegns)
        read(CEMAFilN, *) MessageTemp, (ConsolidationType(i), i = 1, NumConsolidRegns)
        read(CEMAFilN, *) MessageTemp, (ConstConsolidRate(i), i = 1, NumConsolidRegns)
        read(CEMAFilN, *) MessageTemp, ConsolidRateRegnFil
        read(CEMAFilN, *) MessageTemp, WriteBESnp
        read(CEMAFilN, *) MessageTemp, WritePWSnp
    else
        CEMASedimentProcessesInc = .false.
        do JSKIP = 1, 9
            read(CEMAFilN, *)
        end do
    end if
! GROUP 4: Sediment Diagenesis
    read(CEMAFilN, *) MessageTemp, IncludeCEMASedDiagenesis
    if (IncludeCEMASedDiagenesis) then
        SEDIMENT_CALC = .false.
        SOD = 0.0
    end if
    read(CEMAFilN, *) MessageTemp, BedElevationInit
    read(CEMAFilN, *) MessageTemp, BedPorosityInit
    read(CEMAFilN, *) MessageTemp, CEMAParticleSize
    CEMAParticleSize = 1.d-6*CEMAParticleSize !Microns to m
    read(CEMAFilN, *) MessageTemp, CEMASedimentType
    read(CEMAFilN, *) MessageTemp, CEMASedimentDensity
    read(CEMAFilN, *) MessageTemp, CEMASedimentSVelocity
    CEMASedimentSVelocity = CEMASedimentSVelocity/DAY !m/d to m/s
    read(CEMAFilN, *) MessageTemp, CEMASedimentProcessesInc
!
    allocate(BedElevation(IMX), BedElevationLayer(IMX), BedPorosity(IMX))
    allocate(ConsolidRegnNum(IMX), BedConsolidRate(IMX), PorewaterRelRate(IMX))
    allocate(CEMASedConc(IMX, KMX))
    allocate(CEMACumPWRelease(IMX), CEMALayerAdded(IMX), CEMASSApplied(IMX))
    allocate(CEMACumPWToRelease(IMX), CEMACumPWReleased(IMX))
    allocate(NumCEMAPWInst(IMX))
    allocate(ApplyCEMAPWRelease(IMX))
    allocate(CEMACumPWReleaseRate(IMX))
    allocate(EndBedConsolidation(IMX), BedConsolidationSeg(IMX))
    allocate(CEMATSSCopy(KMX, IMX))
    allocate(VOLCEMA(NBR))
!
!
    if (IncludeCEMASedDiagenesis) then
        sediment_diagenesis = .true.
        FirstTimeinCEMAMFTSedDiag = .true.
        read(CEMAFilN, *) MessageTemp, Bubbles_Calculation

! GROUP 5: Bubbles
!IF (.NOT. IncludeCEMASedDiagenesis) Bubbles_Calculation = .FALSE.
        if (Bubbles_Calculation) then
            read(CEMAFilN, *) MessageTemp, GasDiff_Sed ! in m^2/s
            read(CEMAFilN, *) MessageTemp, CalibParam_R1
            read(CEMAFilN, *) MessageTemp, YoungModulus
            read(CEMAFilN, *) MessageTemp, CritStressIF
            read(CEMAFilN, *) MessageTemp, BubbRelScale
            read(CEMAFilN, *) MessageTemp, CrackCloseFraction
            read(CEMAFilN, *) MessageTemp, LimBubbSize
            read(CEMAFilN, *) MessageTemp, MaxBubbRad
            read(CEMAFilN, *) MessageTemp, UseReleaseFraction
            read(CEMAFilN, *) MessageTemp, BubbRelFraction
            read(CEMAFilN, *) MessageTemp, BubbAccFraction
            read(CEMAFilN, *) MessageTemp, NumBubRelArr
            read(CEMAFilN, *) MessageTemp, BubbRelFractionAtm
            read(CEMAFilN, *) MessageTemp, BubbWatGasExchRate
            read(CEMAFilN, *) MessageTemp, ApplyBubbTurb
            read(CEMAFilN, *) MessageTemp, CEMATurbulenceScaling
            allocate(BubblesCarried(IMX, NumBubRelArr), BubblesRadius(IMX, NumBubRelArr))
            allocate(BubblesLNumber(IMX, NumBubRelArr), BubblesStatus(IMX, NumBubRelArr))
            allocate(BubblesRiseV(IMX, NumBubRelArr))
            allocate(BubblesGasConc(IMX, NumBubRelArr, NumGas))
            allocate(BRVoluAGas(IMX, NumBubRelArr, NumGas), BRRateAGas(IMX, NumBubRelArr, NumGas))
!Allocate(FirstBubblesRelease(IMX,NumBubRelArr), BubblesReleaseAllValue(IMX,NumBubRelArr))    ! SW 7/1/2017
            allocate(BubblesAtSurface(IMX, NumBubRelArr))
        else
            do JSKIP = 1, 16
                read(CEMAFilN, *)
            end do
            LimBubbSize = .false.
            UseReleaseFraction = .false.
            ApplyBubbTurb = .false.
        end if

        read(CEMAFilN, *) MessageTemp, CEMA_POM_Resuspension

        if (CEMA_POM_Resuspension) then
            read(CEMAFilN, *) MessageTemp, TAUCRPOM
            read(CEMAFilN, *) MessageTemp, crshields
            read(CEMAFilN, *) MessageTemp, cao_method
            read(CEMAFilN, *) MessageTemp, spgrav_POM
            read(CEMAFilN, *) MessageTemp, dia_POM
        else
            do JSKIP = 1, 5
                read(CEMAFilN, *)
            end do
        end if

        read(CEMAFilN, *) MessageTemp, IncludeAlkalinity
        read(CEMAFilN, *) MessageTemp, IncludeIron
        read(CEMAFilN, *) MessageTemp, IncludeManganese
!
        if (IncludeAlkalinity) then
            IncludeDynamicpH = .true.
        end if
!
!IF(.NOT. IncludeBedConsolidation) THEN
!    Read(CEMAFilN,*)MessageTemp, BedElevationInit
!    Read(CEMAFilN,*)MessageTemp, BedPorosityInit
!    Read(CEMAFilN,*)MessageTemp, CEMASedimentDensity
!    Allocate(BedElevation(IMX), BedPorosity(IMX),EndBedConsolidation(IMX), PorewaterRelRate(IMX))
!END IF
!
        read(CEMAFilN, *) MessageTemp, NumRegnsSedimentBedComposition
        allocate(SDRegnPOC_T(NumRegnsSedimentBedComposition), SDRegnPON_T(NumRegnsSedimentBedComposition), SDRegnSul_T(NumRegnsSedimentBedComposition))
        allocate(SDRegnPOP_T(NumRegnsSedimentBedComposition))
        allocate(SDRegnH2S_T(NumRegnsSedimentBedComposition), SDRegnNH3_T(NumRegnsSedimentBedComposition), SDRegnCH4_T(NumRegnsSedimentBedComposition))
        allocate(SDRegnTIC_T(NumRegnsSedimentBedComposition), SDRegnPO4_T(NumRegnsSedimentBedComposition), SDRegnNO3_T(NumRegnsSedimentBedComposition))
        if (IncludeAlkalinity) then
            allocate(SDRegnALK_T(NumRegnsSedimentBedComposition))
        end if
        if (IncludeIron) then
            allocate(SDRegnFe2_T(NumRegnsSedimentBedComposition), SDRegnFeOOH_T(NumRegnsSedimentBedComposition))
        end if
        if (IncludeManganese) then
            allocate(SDRegnMn2_T(NumRegnsSedimentBedComposition), SDRegnMnO2_T(NumRegnsSedimentBedComposition))
        end if
        allocate(SDRegnT_T(NumRegnsSedimentBedComposition))
        if (.not. IncludeDynamicpH) then
            allocate(SDRegnpH(NumRegnsSedimentBedComposition))
        end if
        allocate(SedBedInitRegSegSt(NumRegnsSedimentBedComposition), SedBedInitRegSegEn(NumRegnsSedimentBedComposition))
!
        read(CEMAFilN, *) ! skip line for header
        read(CEMAFilN, *) MessageTemp, (SedBedInitRegSegSt(i), i = 1, NumRegnsSedimentBedComposition)
        read(CEMAFilN, *) MessageTemp, (SedBedInitRegSegEn(i), i = 1, NumRegnsSedimentBedComposition)
        read(CEMAFilN, *) MessageTemp, (SDRegnT_T(i), i = 1, NumRegnsSedimentBedComposition)
        if (.not. IncludeDynamicpH) then
            read(CEMAFilN, *) MessageTemp, (SDRegnpH(i), i = 1, NumRegnsSedimentBedComposition)
        else
            read(CEMAFilN, *)
        end if
        read(CEMAFilN, *) MessageTemp, (SDRegnPOC_T(i), i = 1, NumRegnsSedimentBedComposition)
        read(CEMAFilN, *) MessageTemp, (SDRegnPON_T(i), i = 1, NumRegnsSedimentBedComposition)
        read(CEMAFilN, *) MessageTemp, (SDRegnPOP_T(i), i = 1, NumRegnsSedimentBedComposition)
        read(CEMAFilN, *) MessageTemp, (SDRegnSul_T(i), i = 1, NumRegnsSedimentBedComposition)
        read(CEMAFilN, *) MessageTemp, (SDRegnNH3_T(i), i = 1, NumRegnsSedimentBedComposition)
        read(CEMAFilN, *) MessageTemp, (SDRegnNO3_T(i), i = 1, NumRegnsSedimentBedComposition)
        read(CEMAFilN, *) MessageTemp, (SDRegnPO4_T(i), i = 1, NumRegnsSedimentBedComposition)
        read(CEMAFilN, *) MessageTemp, (SDRegnH2S_T(i), i = 1, NumRegnsSedimentBedComposition)
        read(CEMAFilN, *) MessageTemp, (SDRegnCH4_T(i), i = 1, NumRegnsSedimentBedComposition)
        read(CEMAFilN, *) MessageTemp, (SDRegnTIC_T(i), i = 1, NumRegnsSedimentBedComposition)
        if (IncludeAlkalinity) then
            read(CEMAFilN, *) MessageTemp, (SDRegnALK_T(i), i = 1, NumRegnsSedimentBedComposition)
        else
            read(CEMAFilN, *)
        end if
!
        if (IncludeIron) then
            read(CEMAFilN, *) MessageTemp, (SDRegnFe2_T(i), i = 1, NumRegnsSedimentBedComposition)
            read(CEMAFilN, *) MessageTemp, (SDRegnFeOOH_T(i), i = 1, NumRegnsSedimentBedComposition)
        else
            do JSKIP = 1, 2
                read(CEMAFilN, *)
            end do
        end if
        if (IncludeManganese) then
            read(CEMAFilN, *) MessageTemp, (SDRegnMn2_T(i), i = 1, NumRegnsSedimentBedComposition)
            read(CEMAFilN, *) MessageTemp, (SDRegnMnO2_T(i), i = 1, NumRegnsSedimentBedComposition)
        else
            do JSKIP = 1, 2
                read(CEMAFilN, *)
            end do
        end if
!
        read(CEMAFilN, *) MessageTemp, NumRegnsSedimentDiagenesis
        read(CEMAFilN, *) ! SKIP LINE FOR HEADER
        allocate(SDRegnPOC_L_Fr(NumRegnsSedimentDiagenesis), SDRegnPOC_R_Fr(NumRegnsSedimentDiagenesis), SDRegnPON_L_Fr(NumRegnsSedimentDiagenesis))
        allocate(SDRegnPON_R_Fr(NumRegnsSedimentDiagenesis), SDRegnPW_DiffCoeff(NumRegnsSedimentDiagenesis), SDRegnOx_Threshold(NumRegnsSedimentDiagenesis))
        allocate(SDRegnPOP_L_Fr(NumRegnsSedimentDiagenesis), SDRegnPOP_R_Fr(NumRegnsSedimentDiagenesis))
        allocate(SDRegnAe_NH3_NO3_L(NumRegnsSedimentDiagenesis), SDRegnAe_NH3_NO3_H(NumRegnsSedimentDiagenesis), SDRegnAe_NO3_N2_L(NumRegnsSedimentDiagenesis))
        allocate(SDRegnAe_NO3_N2_H(NumRegnsSedimentDiagenesis), SDRegnAn_NO3_N2(NumRegnsSedimentDiagenesis), SDRegnAe_CH4_CO2(NumRegnsSedimentDiagenesis))
        allocate(SDRegnAe_HS_NH4_Nit(NumRegnsSedimentDiagenesis), SDRegnAe_HS_O2_Nit(NumRegnsSedimentDiagenesis), SDRegn_Theta_PW(NumRegnsSedimentDiagenesis), SDRegn_Theta_PM(NumRegnsSedimentDiagenesis))
        allocate(SDRegn_Theta_NH3_NO3(NumRegnsSedimentDiagenesis), SDRegn_Theta_NO3_N2(NumRegnsSedimentDiagenesis), SDRegn_Theta_CH4_CO2(NumRegnsSedimentDiagenesis))
        allocate(SDRegn_Sulfate_CH4_H2S(NumRegnsSedimentDiagenesis), SDRegnAe_H2S_SO4(NumRegnsSedimentDiagenesis), SDRegn_Theta_H2S_SO4(NumRegnsSedimentDiagenesis))
        allocate(SDRegn_NormConst_H2S_SO4(NumRegnsSedimentDiagenesis), SDRegn_MinRate_PON_Lab(NumRegnsSedimentDiagenesis), SDRegn_MinRate_PON_Ref(NumRegnsSedimentDiagenesis))
        allocate(SDRegn_MinRate_PON_Ine(NumRegnsSedimentDiagenesis), SDRegn_MinRate_POC_Lab(NumRegnsSedimentDiagenesis), SDRegn_MinRate_POC_Ref(NumRegnsSedimentDiagenesis))
        allocate(SDRegn_MinRate_POC_Ine(NumRegnsSedimentDiagenesis), SDRegn_Theta_PON_Lab(NumRegnsSedimentDiagenesis), SDRegn_Theta_PON_Ref(NumRegnsSedimentDiagenesis))
        allocate(SDRegn_Theta_PON_Ine(NumRegnsSedimentDiagenesis), SDRegn_Theta_POC_Lab(NumRegnsSedimentDiagenesis), SDRegn_Theta_POC_Ref(NumRegnsSedimentDiagenesis))
        allocate(SDRegn_Theta_POC_Ine(NumRegnsSedimentDiagenesis), SDRegn_CH4CompMethod(NumRegnsSedimentDiagenesis), SDRegn_POMResuspMethod(NumRegnsSedimentDiagenesis))
        allocate(SDRegn_Theta_POP_Lab(NumRegnsSedimentDiagenesis), SDRegn_Theta_POP_Ref(NumRegnsSedimentDiagenesis), SDRegn_Theta_POP_Ine(NumRegnsSedimentDiagenesis))
        allocate(SDRegn_MinRate_POP_Lab(NumRegnsSedimentDiagenesis), SDRegn_MinRate_POP_Ref(NumRegnsSedimentDiagenesis), SDRegn_MinRate_POP_Ine(NumRegnsSedimentDiagenesis))
        allocate(SedBedDiaRCRegSegSt(NumRegnsSedimentDiagenesis), SedBedDiaRCRegSegEn(NumRegnsSedimentDiagenesis))
        allocate(Kdp2(NumRegnsSedimentDiagenesis), KdNH31(NumRegnsSedimentDiagenesis), KdNH32(NumRegnsSedimentDiagenesis))
        allocate(delta_kpo41(NumRegnsSedimentDiagenesis), DOcr(NumRegnsSedimentDiagenesis))
        allocate(KsOxch(NumRegnsSedimentDiagenesis))
        allocate(KdH2S1(NumRegnsSedimentDiagenesis), KdH2S2(NumRegnsSedimentDiagenesis))
        allocate(KdFe1(NumRegnsSedimentDiagenesis), KdFe2(NumRegnsSedimentDiagenesis), KdMn1(NumRegnsSedimentDiagenesis), KdMn2(NumRegnsSedimentDiagenesis))
        allocate(PartMixVel(NumRegnsSedimentDiagenesis), BurialVel(NumRegnsSedimentDiagenesis), POCr(NumRegnsSedimentDiagenesis))
!
        DYNAMIC_SD = .false.
        read(CEMAFilN, *) MessageTemp, (SedBedDiaRCRegSegSt(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SedBedDiaRCRegSegEn(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegnPOC_L_Fr(i), i = 1, NumRegnsSedimentDiagenesis)
        if (SDRegnPOC_L_Fr(1) < 0.0) then
            DYNAMIC_SD = .true.
            SDRegnPOC_L_Fr(1) = ABS(SDRegnPOC_L_Fr(1))
        end if
        read(CEMAFilN, *) MessageTemp, (SDRegnPOC_R_Fr(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegnPON_L_Fr(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegnPON_R_Fr(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegnPOP_L_Fr(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegnPOP_R_Fr(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegnPW_DiffCoeff(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (PartMixVel(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (BurialVel(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (POCr(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_CH4CompMethod(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegnOx_Threshold(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegnAe_NH3_NO3_L(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegnAe_NH3_NO3_H(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegnAe_NO3_N2_L(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegnAe_NO3_N2_H(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegnAn_NO3_N2(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegnAe_CH4_CO2(i), i = 1, NumRegnsSedimentDiagenesis) !Eq. 10.35
        read(CEMAFilN, *) MessageTemp, (KsOxch(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegnAe_HS_NH4_Nit(i), i = 1, NumRegnsSedimentDiagenesis) !Eq. 3.3
        read(CEMAFilN, *) MessageTemp, (SDRegnAe_HS_O2_Nit(i), i = 1, NumRegnsSedimentDiagenesis) !Eq. 3.3
        read(CEMAFilN, *) MessageTemp, (SDRegn_Theta_PW(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_Theta_PM(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_Theta_NH3_NO3(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_Theta_NO3_N2(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_Theta_CH4_CO2(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_Sulfate_CH4_H2S(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegnAe_H2S_SO4(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_Theta_H2S_SO4(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_NormConst_H2S_SO4(i), i = 1, NumRegnsSedimentDiagenesis) !Eq. 9.6
        read(CEMAFilN, *) MessageTemp, (SDRegn_MinRate_POC_Lab(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_MinRate_POC_Ref(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_MinRate_POC_Ine(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_MinRate_PON_Lab(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_MinRate_PON_Ref(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_MinRate_PON_Ine(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_MinRate_POP_Lab(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_MinRate_POP_Ref(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_MinRate_POP_Ine(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_Theta_POC_Lab(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_Theta_POC_Ref(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_Theta_POC_Ine(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_Theta_PON_Lab(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_Theta_PON_Ref(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_Theta_PON_Ine(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_Theta_POP_Lab(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_Theta_POP_Ref(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_Theta_POP_Ine(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (Kdp2(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (delta_kpo41(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (DOcr(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (KdNH31(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (KdNH32(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (KdH2S1(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (KdH2S2(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (SDRegn_POMResuspMethod(i), i = 1, NumRegnsSedimentDiagenesis)

        read(CEMAFilN, *) MessageTemp, (KdFe1(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (KdFe2(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (KdMn1(i), i = 1, NumRegnsSedimentDiagenesis)
        read(CEMAFilN, *) MessageTemp, (KdMn2(i), i = 1, NumRegnsSedimentDiagenesis)
! GROUP 6: Output
        read(CEMAFilN, *) MessageTemp, WriteCEMAMFTSedFlx
        read(CEMAFilN, *) MessageTemp, SEDIAGFREQ ! FREQUENCY OF OUTPUT SW 5/25/2017
    else
        IncludeDynamicpH = .false.
        IncludeAlkalinity = .false.
        CEMA_POM_Resuspension = .false.
        IncludeIron = .false.
        IncludeManganese = .false.
        cao_method = .false.
    end if

    close(CEMAFilN)
!
    if (IncludeCEMASedDiagenesis .and. WriteCEMAMFTSedFlx) then
        if (RESTART_IN) then
            open(CEMASedFlxFilN4, File="DiagenesisSOD.csv", POSITION="APPEND")
            JDAY1 = 0.0
            rewind(CEMASedFlxFilN4)
            read(CEMASedFlxFilN4, "(/)", END=101)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN4, "(A,F12.0)", END=101) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN4)
            101 JDAY1 = 0.0
            open(CEMASedFlxFilN5, File="Diagenesis_POCG1.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN5)
            read(CEMASedFlxFilN5, "(/)", END=102)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN5, "(A,F12.0)", END=102) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN5)
            102 JDAY1 = 0.0
            open(CEMASedFlxFilN6, File="Diagenesis_POCG2.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN6)
            read(CEMASedFlxFilN6, "(/)", END=103)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN6, "(A,F12.0)", END=103) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN6)
            103 JDAY1 = 0.0
            open(CEMASedFlxFilN7, File="Diagenesis_JC.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN7)
            read(CEMASedFlxFilN7, "(/)", END=104)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN7, "(A,F12.0)", END=104) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN7)
            104 JDAY1 = 0.0
            open(CEMASedFlxFilN8, File="Diagenesis_JN.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN8)
            read(CEMASedFlxFilN8, "(/)", END=105)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN8, "(A,F12.0)", END=105) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN8)
            105 JDAY1 = 0.0
            open(CEMASedFlxFilN9, File="Diagenesis_PONG1.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN9)
            read(CEMASedFlxFilN9, "(/)", END=106)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN9, "(A,F12.0)", END=106) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN9)
            106 JDAY1 = 0.0
            open(CEMASedFlxFilN10, File="Diagenesis_PONG2.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN10)
            read(CEMASedFlxFilN10, "(/)", END=107)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN10, "(A,F12.0)", END=107) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN10)
            107 JDAY1 = 0.0
            open(CEMASedFlxFilN11, File="Diagenesis_SD_JCH4.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN11)
            read(CEMASedFlxFilN11, "(/)", END=108)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN11, "(A,F12.0)", END=108) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN11)
            108 JDAY1 = 0.0
            open(CEMASedFlxFilN12, File="Diagenesis_SD_JNH4.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN12)
            read(CEMASedFlxFilN12, "(/)", END=109)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN12, "(A,F12.0)", END=109) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN12)
            109 JDAY1 = 0.0
            open(CEMASedFlxFilN13, File="Diagenesis_SD_JNO3.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN13)
            read(CEMASedFlxFilN13, "(/)", END=110)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN13, "(A,F12.0)", END=110) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN13)
            110 JDAY1 = 0.0
            open(CEMASedFlxFilN14, File="Diagenesis_SD_JPO4.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN14)
            read(CEMASedFlxFilN14, "(/)", END=111)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN14, "(A,F12.0)", END=111) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN14)
            111 JDAY1 = 0.0
            open(CEMASedFlxFilN15, File="Diagenesis_POPG1.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN15)
            read(CEMASedFlxFilN15, "(/)", END=112)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN15, "(A,F12.0)", END=112) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN15)
            112 JDAY1 = 0.0
            open(CEMASedFlxFilN16, File="Diagenesis_POPG2.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN16)
            read(CEMASedFlxFilN16, "(/)", END=113)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN16, "(A,F12.0)", END=113) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN16)
            113 JDAY1 = 0.0
            open(CEMASedFlxFilN17, File="DiagenesisCSOD.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN17)
            read(CEMASedFlxFilN17, "(/)", END=114)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN17, "(A,F12.0)", END=114) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN17)
            114 JDAY1 = 0.0
            open(CEMASedFlxFilN18, File="DiagenesisNSOD.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN18)
            read(CEMASedFlxFilN18, "(/)", END=115)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN18, "(A,F12.0)", END=115) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN18)
            115 JDAY1 = 0.0
            open(CEMASedFlxFilN19, File="Diagenesis_JP.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN19)
            read(CEMASedFlxFilN19, "(/)", END=116)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN19, "(A,F12.0)", END=116) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN19)
            116 JDAY1 = 0.0
            open(CEMASedFlxFilN20, File="DiagenesisAerobicLayer.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN20)
            read(CEMASedFlxFilN20, "(/)", END=117)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN20, "(A,F12.0)", END=117) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN20)
            117 JDAY1 = 0.0
            open(CEMASedFlxFilN21, File="Diagenesis_TemperatureAerobicLayer.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN21)
            read(CEMASedFlxFilN21, "(/)", END=118)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN21, "(A,F12.0)", END=118) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN21)
            118 JDAY1 = 0.0
            open(CEMASedFlxFilN22, File="Diagenesis_TemperatureAnaerobicLayer.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN22)
            read(CEMASedFlxFilN22, "(/)", END=119)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN22, "(A,F12.0)", END=119) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN22)
            119 JDAY1 = 0.0

            open(CEMASedFlxFilN23, File="Diagenesis_NO3AerobicLayer.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN23)
            read(CEMASedFlxFilN23, "(/)", END=120)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN23, "(A,F12.0)", END=120) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN23)
            120 JDAY1 = 0.0
            open(CEMASedFlxFilN24, File="Diagenesis_NO3AnaerobicLayer.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN24)
            read(CEMASedFlxFilN24, "(/)", END=121)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN24, "(A,F12.0)", END=121) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN24)
            121 JDAY1 = 0.0
            open(CEMASedFlxFilN25, File="Diagenesis_NH3AerobicLayer.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN25)
            read(CEMASedFlxFilN25, "(/)", END=122)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN25, "(A,F12.0)", END=122) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN25)
            122 JDAY1 = 0.0
            open(CEMASedFlxFilN26, File="Diagenesis_NH3AnaerobicLayer.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN26)
            read(CEMASedFlxFilN26, "(/)", END=123)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN26, "(A,F12.0)", END=123) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN26)
            123 JDAY1 = 0.0
            open(CEMASedFlxFilN27, File="Diagenesis_PO4AerobicLayer.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN27)
            read(CEMASedFlxFilN27, "(/)", END=124)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN27, "(A,F12.0)", END=124) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN27)
            124 JDAY1 = 0.0
            open(CEMASedFlxFilN28, File="Diagenesis_PO4AnaerobicLayer.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN28)
            read(CEMASedFlxFilN28, "(/)", END=125)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN28, "(A,F12.0)", END=125) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN28)
            125 JDAY1 = 0.0
            open(CEMASedFlxFilN29, File="Diagenesis_SO4AerobicLayer.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN29)
            read(CEMASedFlxFilN29, "(/)", END=126)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN29, "(A,F12.0)", END=126) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN29)
            126 JDAY1 = 0.0
            open(CEMASedFlxFilN30, File="Diagenesis_SO4AnaerobicLayer.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN30)
            read(CEMASedFlxFilN30, "(/)", END=127)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN30, "(A,F12.0)", END=127) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN30)
            127 JDAY1 = 0.0

            open(CEMASedFlxFilN31, File="Diagenesis_FeIIAerobicLayer.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN31)
            read(CEMASedFlxFilN31, "(/)", END=128)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN31, "(A,F12.0)", END=128) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN31)
            128 JDAY1 = 0.0
            open(CEMASedFlxFilN32, File="Diagenesis_FeIIAnaerobicLayer.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN32)
            read(CEMASedFlxFilN32, "(/)", END=129)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN32, "(A,F12.0)", END=129) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN32)
            129 JDAY1 = 0.0
            open(CEMASedFlxFilN33, File="Diagenesis_MnIIAerobicLayer.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN33)
            read(CEMASedFlxFilN33, "(/)", END=130)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN33, "(A,F12.0)", END=130) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN33)
            130 JDAY1 = 0.0
            open(CEMASedFlxFilN34, File="Diagenesis_MnIIAnaerobicLayer.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN34)
            read(CEMASedFlxFilN34, "(/)", END=131)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN34, "(A,F12.0)", END=131) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN34)
            131 JDAY1 = 0.0
            open(CEMASedFlxFilN35, File="Diagenesis_CH4AerobicLayer.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN35)
            read(CEMASedFlxFilN35, "(/)", END=132)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN35, "(A,F12.0)", END=132) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN35)
            132 JDAY1 = 0.0
            open(CEMASedFlxFilN36, File="Diagenesis_CH4AnaerobicLayer.csv", POSITION="APPEND")
            rewind(CEMASedFlxFilN36)
            read(CEMASedFlxFilN36, "(/)", END=133)
            do while (JDAY1 < JDAY)
                read(CEMASedFlxFilN36, "(A,F12.0)", END=133) ADUMMY, JDAY1
            end do
            backspace(CEMASedFlxFilN36)
            133 JDAY1 = 0.0

        else

            open(CEMASedFlxFilN4, File="DiagenesisSOD.csv", STATUS="unknown")
            write(CEMASedFlxFilN4, '("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))') (SegNumI, SegNumI = 1, IMX), (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN5, File="Diagenesis_POCG1.csv", STATUS="unknown")
            write(CEMASedFlxFilN5, '("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))') (SegNumI, SegNumI = 1, IMX), (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN6, File="Diagenesis_POCG2.csv", STATUS="unknown")
            write(CEMASedFlxFilN6, '("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))') (SegNumI, SegNumI = 1, IMX), (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN7, File="Diagenesis_JC.csv", STATUS="unknown")
            write(CEMASedFlxFilN7, '("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))') (SegNumI, SegNumI = 1, IMX), (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN8, File="Diagenesis_JN.csv", STATUS="unknown")
            write(CEMASedFlxFilN8, '("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))') (SegNumI, SegNumI = 1, IMX), (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN9, File="Diagenesis_PONG1.csv", STATUS="unknown")
            write(CEMASedFlxFilN9, '("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))') (SegNumI, SegNumI = 1, IMX), (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN10, File="Diagenesis_PONG2.csv", STATUS="unknown")
            write(CEMASedFlxFilN10, '("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))') (SegNumI, SegNumI = 1, IMX), (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN11, File="Diagenesis_SD_JCH4.csv", STATUS="unknown")
            write(CEMASedFlxFilN11, '("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))') (SegNumI, SegNumI = 1, IMX), (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN12, File="Diagenesis_SD_JNH4.csv", STATUS="unknown")
            write(CEMASedFlxFilN12, '("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))') (SegNumI, SegNumI = 1, IMX), (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN13, File="Diagenesis_SD_JNO3.csv", STATUS="unknown")
            write(CEMASedFlxFilN13, '("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))') (SegNumI, SegNumI = 1, IMX), (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN14, File="Diagenesis_SD_JPO4.csv", STATUS="unknown")
            write(CEMASedFlxFilN14, '("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))') (SegNumI, SegNumI = 1, IMX), (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN15, File="Diagenesis_POPG1.csv", STATUS="unknown")
            write(CEMASedFlxFilN15, '("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))') (SegNumI, SegNumI = 1, IMX), (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN16, File="Diagenesis_POPG2.csv", STATUS="unknown")
            write(CEMASedFlxFilN16, '("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))') (SegNumI, SegNumI = 1, IMX), (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN17, File="DiagenesisCSOD.csv", STATUS="unknown")
            write(CEMASedFlxFilN17, '("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))') (SegNumI, SegNumI = 1, IMX), (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN18, File="DiagenesisNSOD.csv", STATUS="unknown")
            write(CEMASedFlxFilN18, '("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))') (SegNumI, SegNumI = 1, IMX), (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN19, File="Diagenesis_JP.csv", STATUS="unknown")
            write(CEMASedFlxFilN19, '("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))') (SegNumI, SegNumI = 1, IMX), (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN20, File="DiagenesisAerobicLayer.csv", STATUS="unknown")
            write(CEMASedFlxFilN20, '("Variable,JDAY,",<IMX>(i5,","))') (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN21, File="Diagenesis_TemperatureAerobicLayer.csv", STATUS="unknown")
            write(CEMASedFlxFilN21, '("Variable,JDAY,",<IMX>(i5,","))') (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN22, File="Diagenesis_TemperatureAnaerobicLayer.csv", STATUS="unknown")
            write(CEMASedFlxFilN22, '("Variable,JDAY,",<IMX>(i5,","))') (SegNumI, SegNumI = 1, IMX)

            open(CEMASedFlxFilN23, File="Diagenesis_NO3AerobicLayer.csv", STATUS="unknown")
            write(CEMASedFlxFilN23, '("Variable,JDAY,",<IMX>(i5,","))') (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN24, File="Diagenesis_NO3AnaerobicLayer.csv", STATUS="unknown")
            write(CEMASedFlxFilN24, '("Variable,JDAY,",<IMX>(i5,","))') (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN25, File="Diagenesis_NH3AerobicLayer.csv", STATUS="unknown")
            write(CEMASedFlxFilN25, '("Variable,JDAY,",<IMX>(i5,","))') (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN26, File="Diagenesis_NH3AnaerobicLayer.csv", STATUS="unknown")
            write(CEMASedFlxFilN26, '("Variable,JDAY,",<IMX>(i5,","))') (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN27, File="Diagenesis_PO4AerobicLayer.csv", STATUS="unknown")
            write(CEMASedFlxFilN27, '("Variable,JDAY,",<IMX>(i5,","))') (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN28, File="Diagenesis_PO4AnaerobicLayer.csv", STATUS="unknown")
            write(CEMASedFlxFilN28, '("Variable,JDAY,",<IMX>(i5,","))') (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN29, File="Diagenesis_SO4AerobicLayer.csv", STATUS="unknown")
            write(CEMASedFlxFilN29, '("Variable,JDAY,",<IMX>(i5,","))') (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN30, File="Diagenesis_SO4AnaerobicLayer.csv", STATUS="unknown")
            write(CEMASedFlxFilN30, '("Variable,JDAY,",<IMX>(i5,","))') (SegNumI, SegNumI = 1, IMX)

            open(CEMASedFlxFilN31, File="Diagenesis_FeIIAerobicLayer.csv", STATUS="unknown")
            write(CEMASedFlxFilN31, '("Variable,JDAY,",<IMX>(i5,","))') (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN32, File="Diagenesis_FeIIAnaerobicLayer.csv", STATUS="unknown")
            write(CEMASedFlxFilN32, '("Variable,JDAY,",<IMX>(i5,","))') (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN33, File="Diagenesis_MnIIAerobicLayer.csv", STATUS="unknown")
            write(CEMASedFlxFilN33, '("Variable,JDAY,",<IMX>(i5,","))') (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN34, File="Diagenesis_MnIIAnaerobicLayer.csv", STATUS="unknown")
            write(CEMASedFlxFilN34, '("Variable,JDAY,",<IMX>(i5,","))') (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN35, File="Diagenesis_CH4AerobicLayer.csv", STATUS="unknown")
            write(CEMASedFlxFilN35, '("Variable,JDAY,",<IMX>(i5,","))') (SegNumI, SegNumI = 1, IMX)
            open(CEMASedFlxFilN36, File="Diagenesis_CH4AnaerobicLayer.csv", STATUS="unknown")
            write(CEMASedFlxFilN36, '("Variable,JDAY,",<IMX>(i5,","))') (SegNumI, SegNumI = 1, IMX)

        end if

    end if
!
! Fix the values here
    NH4_NH3_Eqb_Const = 9.1
    HS_H2S_Eqb_Const = 9.0
    HenryConst_NH3 = 0.0179
    HenryConst_CH4 = 469.0
    HenryConst_H2S = 10.0
    HenryConst_CO2 = 29.0
!
!Allocate other variablesnd 
    allocate(CellArea(KMX, IMX))
    if (IncludeCEMASedDiagenesis) then
        allocate(CEMAMFT_RandC_RegN(IMX), CEMAMFT_InCond_RegN(IMX), MFTSedFlxVars(KMX, IMX, 59), CEMA_SD_Vars(KMX, IMX, 22))
        allocate(SD_NO3p2(2), SD_NH3p2(2), SD_NH3Tp2(2), SD_CH4p2(2), SD_PO4p2(2), SD_PO4Tp2(2), SD_PO4(2))
        allocate(SD_HSp2(2), SD_HSTp2(2))
        allocate(SD_poc2(3), SD_pon2(3), SD_pop2(3), SD_NH3Tp(2), SD_NO3p(2), SD_PO4Tp(2), SD_HSTp(2))
        allocate(SD_fpon(3), SD_fpoc(3), SD_kdiaPON(3), SD_ThtaPON(3), SD_kdiaPOC(3), SD_ThtaPOC(3))
        allocate(SD_JPOC(3), SD_JPON(3), SD_JPOP(3), SD_TDS(2))
        allocate(SD_EPOC(3), SD_EPON(3), SD_EPOP(3))
        allocate(SD_Denit(2), SD_JDenit(2), SD_JO2NO3(2), SD_HS(2))
        if (IncludeIron) then
            allocate(SD_Fe2(2))
        end if
        if (IncludeManganese) then
            allocate(SD_Mn2(2))
        end if
        allocate(SD_kdiaPOP(3), SD_ThtaPOP(3), SD_NH3T(2), SD_FPOP(3))
        allocate(SD_pHValue(IMX))
        allocate(SD_AerLayerThick(IMX))
        if (Bubbles_Calculation) then
            allocate(H2SDis(IMX), H2SGas(IMX), CH4Dis(IMX), CH4Gas(IMX))
            allocate(NH4Dis(IMX), NH4Gas(IMX), CO2Dis(IMX), CO2Gas(IMX))
            allocate(BubbleRadiusSed(IMX), PresBubbSed(IMX), PresCritSed(IMX))
            allocate(CgSed(IMX), C0Sed(IMX), CtSed(IMX))
            allocate(TConc(NumGas, KMX, IMX), TConcP(NumGas, KMX, IMX), SConc(NumGas, KMX, IMX))
            allocate(DissolvedGasSediments(NumGas, KMX, IMX))
            allocate(CrackOpen(IMX), MFTBubbReleased(IMX), LastDiffVolume(IMX))
            allocate(BottomTurbulence(IMX))
            allocate(FirstBubblesRelease(IMX, NumBubRelArr), BubblesReleaseAllValue(IMX, NumBubRelArr), BubbleRelWB(NWB, NumGas)) ! SW 7/1/2017
            allocate(BRRateAGasNet(IMX, NumGas))
        end if
        allocate(SDPFLUX(NWB), SDNH4FLUX(NWB), SDNO3FLUX(NWB))
    end if

    return
end subroutine CEMA_W2_Input



subroutine INIT_CEMA()
    use CEMAVars;     use MAIN
    implicit none

    if (IncludeCEMASedDiagenesis) then
        SD_NO3p2 = 0.d00;         SD_NH3p2 = 0.d00;         SD_NH3Tp2 = 0.d00;         SD_CH4p2 = 0.d00
        SD_PO4p2 = 0.d00;         SD_PO4Tp2 = 0.d00;         SD_HSp2 = 0.d00;         SD_HSTp2 = 0.d00
        SD_POC2 = 0.d00;         SD_PON2 = 0.d00;         SD_POP2 = 0.d00;         SD_NH3Tp = 0.d00
        SD_NO3p = 0.d00;         SD_PO4Tp = 0.d00;         SD_HSTp = 0.d00
        SD_FPON = 0.d00;         SD_FPOC = 0.d00;         SD_kdiaPON = 0.d00;         SD_ThtaPON = 0.d00
        SD_kdiaPOC = 0.d00;         SD_ThtaPOC = 0.d00
        SD_JPOC = 0.d00;         SD_JPON = 0.d00;         SD_JPOP = 0.d00
        SD_EPOC = 0.d00;         SD_EPON = 0.d00;         SD_EPOP = 0.d00
        SD_Denit = 0.d00;         SD_JDenit = 0.d00;         SD_JO2NO3 = 0.d00
        SD_PO4 = 0.d00;         SD_FPOP = 0.d00;         SD_HS = 0.d00
!
        if (IncludeIron) then
            SD_Fe2 = 0.d00
        end if
        if (IncludeManganese) then
            SD_Mn2 = 0.d00
        end if
        SD_kdiaPOP = 0.d00;         SD_ThtaPOP = 0.d00;         SD_NH3T = 0.d00
        SD_AerLayerThick = 0.d00
!
        if (Bubbles_Calculation) then
            H2SDis = 0.d00;             H2SGas = 0.d00;             CH4Dis = 0.d00;             CH4Gas = 0.d00
            NH4Dis = 0.d00;             NH4Gas = 0.d00;             CO2Dis = 0.d00;             CO2Gas = 0.d00
            BubbleRadiusSed = 0.d00;             PresBubbSed = 0.d00;             PresCritSed = 0.d00
            CgSed = 0.d00;             C0Sed = 0.d00;             CtSed = 0.d00;             TConcP = 0.d00
            LastDiffVolume = 0.d00
            BubblesCarried = 0;             BubblesLNumber = 0;             BubblesStatus = 0
            BubblesRadius = 0.d00;             BubblesRiseV = 0.d00;             BubblesGasConc = 0.d00
            BubblesReleaseAllValue = 0.d00
            BRVoluAGas = 0.d00;             BRRateAGas = 0.d00;             BRRateAGasNet = 0.d00
            BottomTurbulence = 0.d00
            DissolvedGasSediments = 0.d00
            FirstTimeInBubbles = .true.
            FirstBubblesRelease = .true.
            BubblesAtSurface = .false.
        end if
        CEMAMFT_RandC_RegN = 0
        CEMA_SD_Vars = 0.d00
        SDPFLUX = 0.0
        SDNH4FLUX = 0.0
        SDNO3FLUX = 0.0
    end if
!
    if (IncludeBedConsolidation) then
        BedElevationLayer = 0.d00
        BedConsolidRate = 0.d00
        PorewaterRelRate = 0.d00
        CEMASedConc = 0.d00
        CEMACumPWRelease = 0.d00
        CEMACumPWReleaseRate = 0.d00
        CEMACumPWToRelease = 0.d00
        CEMACumPWReleased = 0.d00
        EndBedConsolidation = .false.
        BedConsolidationSeg = .false. ! cb 6/28/18
        VOLCEMA = 0.d00
        NumCEMAPWInst = 0
        ApplyCEMAPWRelease = .false.
    end if
!
    if (IncludeCEMASedDiagenesis .and. .not. IncludeBedConsolidation) then
        EndBedConsolidation = .false.
        PorewaterRelRate = 0.d00
    end if
    if (IncludeCEMASedDiagenesis .or. IncludeBedConsolidation) then
        BedElevation = BedElevationInit
    end if
!
    if (.not. RESTART_IN) then
        CellArea = 0.0
        if (IncludeCEMASedDiagenesis) then
            MFTSedFlxVars = 0.d00
            BedPorosity = BedPorosityInit
            if (Bubbles_Calculation) then
                MFTBubbReleased = 0
                TConc = 0.d00;                 SConc = 0.d00
                CrackOpen = .false.;                 BubbleRelWB = 0.0
                GasReleaseCH4 = 0.0
            end if
        end if
    end if
end subroutine INIT_CEMA


subroutine Deallocate_CEMA()
    use CEMAVars
    implicit none
!
    deallocate(CellArea)
    if (IncludeBedConsolidation) then
        deallocate(ConsolidationType, ConstConsolidRate)
        deallocate(ConstPoreWtrRate, ConsolidRateTemp)
        deallocate(ConsRegSegSt, ConsRegSegEn)
    end if
    deallocate(ConsolidRegnNum, BedConsolidRate, PorewaterRelRate)
    deallocate(CEMASedConc)
    deallocate(CEMACumPWRelease, CEMALayerAdded, CEMASSApplied)
    deallocate(CEMACumPWToRelease, CEMACumPWReleased)
    deallocate(NumCEMAPWInst)
    deallocate(ApplyCEMAPWRelease)
    deallocate(CEMACumPWReleaseRate)
    deallocate(EndBedConsolidation)
    deallocate(BedConsolidationSeg) ! cb 6/28/18
    deallocate(CEMATSSCopy)
    deallocate(VOLCEMA)
    deallocate(BedElevationLayer)
!END IF
!
    if (IncludeCEMASedDiagenesis .or. IncludeBedConsolidation) then
        deallocate(BedElevation, BedPorosity)
    end if
    if (IncludeFFTLayer) then
        deallocate(FFTActPrdSt, FFTActPrdEn)
        deallocate(FFTLayConc)
    end if
    if (IncludeCEMASedDiagenesis) then
        deallocate(SDRegnPOC_T, SDRegnPON_T, SDRegnSul_T)
        deallocate(SDRegnPOP_T)
        deallocate(SDRegnH2S_T, SDRegnNH3_T, SDRegnCH4_T, SDRegnNO3_T)
        if (IncludeAlkalinity) then
            deallocate(SDRegnALK_T)
        end if
        if (.not. IncludeDynamicpH) then
            deallocate(SDRegnpH)
        end if
        deallocate(SDRegnTIC_T, SDRegnPO4_T)
        if (IncludeIron) then
            deallocate(SDRegnFe2_T, SDRegnFeOOH_T, SD_Fe2, KdFe1, KdFe2)
        end if
        if (IncludeManganese) then
            deallocate(SDRegnMn2_T, SDRegnMnO2_T, SD_Mn2, KdMn1, KdMn2)
        end if
        deallocate(SDRegnT_T)
        deallocate(SedBedInitRegSegSt, SedBedInitRegSegEn)
        deallocate(SDRegnPOC_L_Fr, SDRegnPOC_R_Fr, SDRegnPON_L_Fr)
        deallocate(SDRegnPON_R_Fr, SDRegnPW_DiffCoeff, SDRegnOx_Threshold)
        deallocate(SDRegnPOP_L_Fr, SDRegnPOP_R_Fr)
        deallocate(SDRegnAe_NH3_NO3_L, SDRegnAe_NH3_NO3_H, SDRegnAe_NO3_N2_L)
        deallocate(SDRegnAe_NO3_N2_H, SDRegnAn_NO3_N2, SDRegnAe_CH4_CO2)
        deallocate(SDRegnAe_HS_NH4_Nit, SDRegnAe_HS_O2_Nit, SDRegn_Theta_PW, SDRegn_Theta_PM)
        deallocate(SDRegn_Theta_NH3_NO3, SDRegn_Theta_NO3_N2, SDRegn_Theta_CH4_CO2)
        deallocate(SDRegn_Sulfate_CH4_H2S, SDRegnAe_H2S_SO4, SDRegn_Theta_H2S_SO4)
        deallocate(SDRegn_NormConst_H2S_SO4, SDRegn_MinRate_PON_Lab, SDRegn_MinRate_PON_Ref)
        deallocate(SDRegn_MinRate_PON_Ine, SDRegn_MinRate_POC_Lab, SDRegn_MinRate_POC_Ref)
        deallocate(SDRegn_MinRate_POC_Ine, SDRegn_Theta_PON_Lab, SDRegn_Theta_PON_Ref)
        deallocate(SDRegn_Theta_PON_Ine, SDRegn_Theta_POC_Lab, SDRegn_Theta_POC_Ref)
        deallocate(SDRegn_Theta_POC_Ine, SDRegn_CH4CompMethod, SDRegn_POMResuspMethod)
        deallocate(SDRegn_Theta_POP_Lab, SDRegn_Theta_POP_Ref, SDRegn_Theta_POP_Ine)
        deallocate(Kdp2, KdNH31, KdNH32, KdH2S1, KdH2S2, delta_kpo41, DOcr)
        deallocate(PartMixVel, BurialVel, POCr, KsOxch)
        deallocate(SDRegn_MinRate_POP_Lab, SDRegn_MinRate_POP_Ref, SDRegn_MinRate_POP_Ine)
        deallocate(SedBedDiaRCRegSegSt, SedBedDiaRCRegSegEn)
        deallocate(CEMAMFT_RandC_RegN, CEMAMFT_InCond_RegN, MFTSedFlxVars, CEMA_SD_Vars)
        deallocate(SD_NO3p2, SD_NH3p2, SD_NH3Tp2, SD_CH4p2, SD_PO4p2, SD_PO4Tp2, SD_PO4)
        deallocate(SD_HSp2, SD_HSTp2)
        deallocate(SD_poc2, SD_pon2, SD_pop2, SD_NH3Tp, SD_NO3p, SD_PO4Tp, SD_HSTp)
        deallocate(SD_fpon, SD_fpoc, SD_kdiaPON, SD_ThtaPON, SD_kdiaPOC, SD_ThtaPOC)
        deallocate(SD_JPOC, SD_JPON, SD_JPOP, SD_TDS)
        deallocate(SD_EPOC, SD_EPON, SD_EPOP)
        deallocate(SD_Denit, SD_JDenit, SD_JO2NO3, SD_HS) ! cb 7/26/18
        deallocate(SD_kdiaPOP, SD_ThtaPOP, SD_NH3T, SD_FPOP)
        deallocate(SD_pHValue)
        deallocate(SD_AerLayerThick)
        if (Bubbles_Calculation) then
            deallocate(H2SDis, H2SGas, CH4Dis, CH4Gas)
            deallocate(NH4Dis, NH4Gas, CO2Dis, CO2Gas)
            deallocate(BubbleRadiusSed, PresBubbSed, PresCritSed)
            deallocate(CgSed, C0Sed, CtSed)
            deallocate(TConc, TConcP, SConc)
            deallocate(DissolvedGasSediments)
            deallocate(CrackOpen, MFTBubbReleased, LastDiffVolume)
            deallocate(BubblesCarried, BubblesRadius)
            deallocate(BubblesLNumber, BubblesStatus)
            deallocate(BubblesRiseV)
            deallocate(BubbleRelWB)
            deallocate(BubblesGasConc)
            deallocate(BRVoluAGas, BRRateAGas)
            deallocate(FirstBubblesRelease, BubblesReleaseAllValue)
            deallocate(BRRateAGasNet)
            deallocate(BubblesAtSurface)
            deallocate(BottomTurbulence)
        end if
    end if

    return
end subroutine Deallocate_CEMA
