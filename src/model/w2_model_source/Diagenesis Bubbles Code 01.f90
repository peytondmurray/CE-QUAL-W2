subroutine GasBubblesFormation(Radius, DeltaT, Volume)

    use GLOBAL
    use SCREENC
    use CEMAVars
    use CEMASedimentDiagenesis, only: SD_T1

! Type declarations
    implicit none

    integer :: nGas, nRelArr, temp ! nTry
    real(8), allocatable, dimension(:) :: Ctot, CgB, C0B, C1B, Source0
    real(8), allocatable, dimension(:) :: K, Mw, Henry
    real(8) :: Volume, Porosity
    real(8) :: DeltaT, Ro, CgT
    real(8) :: RSI, P0, Nbubbles, NbubblesP, NbubbLost
    real(8) :: Pcrit, Pbubb, PbubbT
    real(8) :: NetMass, DisMass, GasMass
    real(8) :: C1T, C0T, CtT, Radius, Vbub, SourceT
    real(8) :: Vbubbles, DiffVolume, BubSedT
    real(8) :: Source
    logical :: FoundOpenArray

    allocate(Ctot(NumGas), CgB(NumGas), C0B(NumGas), C1B(NumGas))
    allocate(Source0(NumGas), Henry(NumGas), K(NumGas), Mw(NumGas))

    Porosity = BedPorosity(SegNumI)
    Henry(1) = HenryConst_H2S !L atm/M  H2S
    Henry(2) = HenryConst_CH4 !L atm/M  CH4
    Henry(3) = HenryConst_NH3 !L atm/M  NH3
    Henry(4) = HenryConst_CO2 !L atm/M  CO2
    BubSedT = 273.15 + SD_T1 !K           ! cb 5/22/15
    K(1) = Henry(1)/GasConst_R/BubSedT
    K(2) = Henry(2)/GasConst_R/BubSedT
    K(3) = Henry(3)/GasConst_R/BubSedT
    K(4) = Henry(4)/GasConst_R/BubSedT
    Ro = 0.d0 !m
    RSI = 8318.78 !l-N/m�/mol/K
    P0 = 9800. !N/m�
    Mw(1) = 36. !H2S gm/mol
    Mw(2) = 16. !CH4 gm/mol
    Mw(3) = 17. !NH3 gm/mol
    Mw(4) = 44. !CO2 gm/mol
    if (CrackOpen(SegNumI)) then
        NbubbLost = MFTBubbReleased(SegNumI)
    end if

    do nGas = 1, NumGas
        Ctot(nGas) = TConc(nGas, LayerNum, SegNumI)
        Source0(nGas) = SConc(nGas, LayerNum, SegNumI)
    end do !nGas

    if (FirstTimeInBubbles) then

        CgT = 0.d0
        C1T = 0.d0
        C0T = 0.d0
        CtT = 0.d0
        do nGas = 1, NumGas
            C0B(nGas) = Ctot(nGas)/(1 + K(nGas))
            CgB(nGas) = C0B(nGas)*K(nGas)
            C1B(nGas) = Ctot(nGas)
            CgT = CgT + CgB(nGas)
            C1T = C1T + C1B(nGas)
            C0T = C0T + C0B(nGas)
            CtT = CtT + Ctot(nGas)
        end do !nGas
        Radius = sqrt(2.0*Porosity*GasDiff_Sed*DeltaT*(C1T - C0T)/CgT + Ro**2)
        Vbub = 4.0/3.0*3.1415927*Radius**3
        NetMass = CtT*Volume*Porosity
        DisMass = C0T*Volume*Porosity
        GasMass = NetMass - DisMass
        Nbubbles = GasMass/(Vbub*CgT)
        Pcrit = 1.32*(CritStressIF**6/(YoungModulus*Nbubbles*Vbub))**0.2 + P0
        PbubbT = 0.d0
        do nGas = 1, NumGas
            Pbubb = CgB(nGas)*RSI*0.001*BubSedT/Mw(nGas)
            PbubbT = PbubbT + Pbubb
        end do !nGas

    else

        SourceT = 0.d0
        CgT = 0.d0
        C1T = 0.d0
        C0T = 0.d0
        CtT = 0.d0
        do nGas = 1, NumGas

            Source = Source0(nGas)
            SourceT = SourceT + Source
            Ctot(nGas) = Ctot(nGas) + Source*DeltaT
            C0B(nGas) = Ctot(nGas)/(1 + K(nGas))
            CgB(nGas) = C0B(nGas)*K(nGas)
            C1B(nGas) = C0B(nGas)
            CgT = CgT + CgB(nGas)
            C1T = C1T + C1B(nGas)
            C0T = C0T + C0B(nGas)
            CtT = CtT + Ctot(nGas)

        end do !nGas

        Radius = Radius + Porosity*GasDiff_Sed/(Radius*CgT)*(SourceT*CalibParam_R1**2/(6*GasDiff_Sed) + C1T - C0T)*DeltaT

        if (LimBubbSize) then
            if (Radius > MaxBubbRad/1000.0) then
                Radius = MaxBubbRad/1000.0
            end if
        end if

    end if

    Vbub = 4.0/3.0*3.1415927*Radius**3
    NetMass = 0.d0
    DisMass = 0.d0
    CgT = 0.d0
    do nGas = 1, NumGas
        NetMass = NetMass + Ctot(nGas)*Volume*Porosity
        DisMass = DisMass + C0B(nGas)*Volume*Porosity
        CgT = CgT + CgB(nGas)
    end do !nGas
    GasMass = NetMass - DisMass
    Nbubbles = GasMass/(Vbub*CgT)

    Pcrit = 1.324*(CritStressIF**6/(YoungModulus*Nbubbles*Vbub))**0.2 + P0
    PbubbT = 0.d0
    do nGas = 1, NumGas
        Pbubb = CgB(nGas)*RSI*0.001*BubSedT/Mw(nGas)
        PbubbT = PbubbT + Pbubb
    end do !nGas

    PresBubbSed(SegNumI) = PbubbT
    PresCritSed(SegNumI) = Pcrit

    if (PbubbT < Pcrit*CrackCloseFraction) then
        CrackOpen(SegNumI) = .false.
        NbubbLost = 0
    end if

    if (LastDiffVolume(SegNumI) < 0) then
        LastDiffVolume(SegNumI) = 0.d00
    end if

    if (PbubbT > Pcrit .and. .not. CrackOpen(SegNumI)) then

        CrackOpen(SegNumI) = .true.
        Vbubbles = CritStressIF**6/(YoungModulus*((PbubbT - P0)/1.32)**5)
        DiffVolume = Vbub*Nbubbles - Vbubbles
        DiffVolume = DiffVolume*BubbRelScale
        if (UseReleaseFraction) then
            Vbubbles = CritStressIF**6/(YoungModulus*(Pcrit*CrackCloseFraction/1.32)**5)
            DiffVolume = BubbRelFraction*(Vbub*Nbubbles - Vbubbles)
        end if
        LastDiffVolume(SegNumI) = DiffVolume
        NbubbLost = (DiffVolume)/Vbub
        NbubblesP = Nbubbles
        Nbubbles = Nbubbles - NbubbLost
        do nGas = 1, NumGas
            CgB(nGas) = CgB(nGas)*(Vbub*NbubblesP - DiffVolume)/(Vbub*NbubblesP)
            C0B(nGas) = CgB(nGas)/K(nGas)
            Ctot(nGas) = C0B(nGas)*(1 + K(nGas))
        end do !nGas

    end if

!!!!!!!!!!!!!!!!! debug    
!    CrackOpen(SegNumI)=.false.
!!!!!!!!!!!!!!!!!!!!!!! debug
    if (CrackOpen(SegNumI)) then

        NbubblesP = Nbubbles
!Nbubbles = Nbubbles - NbubbLost  ! cb 2/21/13
        do nGas = 1, NumGas
            CgB(nGas) = CgB(nGas)*(Vbub*NbubblesP - LastDiffVolume(SegNumI))/(Vbub*NbubblesP)
            C0B(nGas) = CgB(nGas)/K(nGas)
            Ctot(nGas) = C0B(nGas)*(1 + K(nGas))

            TConcP(nGas, LayerNum, SegNumI) = Ctot(nGas)
            TConc(nGas, LayerNum, SegNumI) = Ctot(nGas)

        end do !nGas

        temp = INT4(NbubbLost)

        FoundOpenArray = .false.
        MFTBubbReleased(SegNumI) = KIDINT(NbubbLost)
!nTry = 0
        do nRelArr = 1, NumBubRelArr
!nTry = nTry + 1
            if (BubblesStatus(SegNumI, nRelArr) == 0) then
                FoundOpenArray = .true.
                BubblesStatus(SegNumI, nRelArr) = 1
                BubblesCarried(SegNumI, nRelArr) = MFTBubbReleased(SegNumI)
                BubblesRadius(SegNumI, nRelArr) = Radius
                BubblesLNumber(SegNumI, nRelArr) = KB(SegNumI)
                do nGas = 1, NumGas
                    BubblesGasConc(SegNumI, nRelArr, nGas) = CgB(nGas)
                end do
                exit
            end if
        end do
        if (.not. FoundOpenArray) then
            write(CEMALogFilN, *) "Insufficient array size for bubbles release at JDAY = ", JDAY
            write(w2err, *) "Insufficient array size for bubbles release at JDAY = ", JDAY
            stop
        end if
    end if

    CgSed(SegNumI) = CgT
    C0Sed(SegNumI) = C0T
    CtSed(SegNumI) = CtT

    return
end subroutine GasBubblesFormation


subroutine CEMACalculateRiseVelocity()

    use MAIN
    use GLOBAL
    use GEOMC
    use CEMAVars
    implicit none

    integer :: nGas, nRelArr
    real(8) :: Rhog

    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            do SegNumI = IU, ID
                do nRelArr = 1, NumBubRelArr
                    if (BubblesStatus(SegNumI, nRelArr) == 0) then

                        do nGas = 1, NumGas
                            BRVoluAGas(SegNumI, nRelArr, nGas) = 0.d00
                            BRRateAGas(SegNumI, nRelArr, nGas) = 0.d00
                            BRRateAGasNet(SegNumI, nGas) = 0.d00
                        end do !nGas  

                    end if
                end do
            end do
        end do
    end do


    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            do SegNumI = IU, ID
                do nRelArr = 1, NumBubRelArr
                    if (BubblesStatus(SegNumI, nRelArr) == 1) then

                        Rhog = 0.d0
                        do nGas = 1, NumGas
                            Rhog = Rhog + BubblesGasConc(SegNumI, nRelArr, nGas)/1000.0 !kg/m�
                        end do !nGas    

                        call CEMABubblesRiseVelocity(BubblesRadius(SegNumI, nRelArr), Rhog, BubblesRiseV(SegNumI, nRelArr))

                    end if
                end do
            end do
        end do
    end do

    return
end subroutine CEMACalculateRiseVelocity


subroutine CEMABubblesRiseVelocity(Radius, Rhog, RiseVelocity)

    implicit none

    real(8) :: Radius, Rhog, RiseVelocity
    real(8) :: Rhow, DynVisc
    real(8) :: Sigma, Nd, W, Reynolds, M, Eo
    real(8) :: J, H
    real(8) :: Radius1, RiseVelocity1, RiseVelocity2

!Calculate Rise Velocity
    Rhow = 1000 !kg/m3
    DynVisc = 0.001002 !kg/m/s
    Sigma = 0.0725 !N/m

    if (Radius*1000 <= 1) then !<= 1 mm
        Nd = 4.*Rhow*(Rhow - Rhog)*9.8*Radius**3/(3.*DynVisc**2)
        W = dlog10(Nd)
        if (Nd <= 73.) then
            Reynolds = Nd/24. - 1.7569d-4*Nd**2 + 6.9252d-7*Nd**3 - 2.3027d-10*Nd**4
        end if
        if (Nd > 73.0 .and. Nd <= 580.0) then
            Reynolds = 10**((-1.7095) + 1.33438*W - 0.11591*W**2)
        end if
        if (Nd > 580.) then
            Reynolds = 10**((-1.81391) + 1.34671*W - 0.12427*W**2 + 0.006344*W**3)
        end if
        RiseVelocity = Reynolds*DynVisc/(Rhow*Radius)
    end if

    if (Radius*1000.0 <= 15.0 .and. Radius*1000.0 > 1.0) then !<= 15 mm
        M = 9.8*DynVisc**4*(Rhow - Rhog)/(Rhow**2*Sigma**3)
        Eo = 9.8*(Rhow - Rhog)*Radius**2/Sigma
        H = 4./3.*Eo*M**(-0.149)*(DynVisc/DynVisc)**(-0.14)
        if (H < 59.3) then
            J = 0.94*H**0.757
        else
            J = 3.42*H**0.441
        end if
        RiseVelocity = DynVisc/(Rhow*Radius)*M**(-0.149)*(J - 0.857)
    end if

    if (Radius*1000.0 > 15.0 .and. Radius*1000.0 <= 18.0) then !<= 15 mm to 18 mm

        Radius1 = Radius
        Radius = 0.015
        M = 9.8*DynVisc**4*(Rhow - Rhog)/(Rhow**2*Sigma**3)
        Eo = 9.8*(Rhow - Rhog)*Radius**2/Sigma
        H = 4./3.*Eo*M**(-0.149)*(DynVisc/DynVisc)**(-0.14)
        if (H < 59.3) then
            J = 0.94*H**0.757
        else
            J = 3.42*H**0.441
        end if
        RiseVelocity1 = DynVisc/(Rhow*Radius)*M**(-0.149)*(J - 0.857)

        Radius = 0.018
        RiseVelocity2 = 0.711*sqrt(9.8*Radius*(Rhow - Rhog)/Rhow)

        Radius = Radius1

        RiseVelocity = ((Radius - 0.015)*RiseVelocity2 + RiseVelocity1*(0.018 - Radius))/(0.018 - 0.015)

    end if

    if (Radius*1000.0 > 18.0) then !> 18 mm
        RiseVelocity = 0.711*sqrt(9.8*Radius*(Rhow - Rhog)/Rhow)
    end if

    return
end subroutine CEMABubblesRiseVelocity



subroutine CEMABubblesTransport()

    use MAIN
    use GLOBAL
    use GEOMC
    use CEMAVars
    implicit none

    integer :: BubbLayer, nRelArr
    real(8) :: VLocationBubble, VDistTravBubble

    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            do SegNumI = IU, ID
                do nRelArr = 1, NumBubRelArr
                    if (BubblesStatus(SegNumI, nRelArr) == 1 .and. .not. BubblesAtSurface(SegNumI, nRelArr)) then

                        BubbLayer = BubblesLNumber(SegNumI, nRelArr)
                        VLocationBubble = 0.5*(el(BubbLayer + 1, SegNumI) + el(BubbLayer, SegNumI))
                        VDistTravBubble = BubblesRiseV(SegNumI, nRelArr)*dlt
                        VLocationBubble = VLocationBubble + VDistTravBubble

!Locate vertical location
                        BubblesLNumber(SegNumI, nRelArr) = KT
                        do K = KT, KB(SegNumI)
                            if (VLocationBubble < el(K, SegNumI)) then
                                BubblesLNumber(SegNumI, nRelArr) = K
                            end if
                        end do !K

                        if (BubblesLNumber(SegNumI, nRelArr) == KT) then
                            BubblesAtSurface(SegNumI, nRelArr) = .true.
                            FirstBubblesRelease(SegNumI, nRelArr) = .true.
                            BubblesReleaseAllValue(SegNumI, nRelArr) = BubbRelFractionAtm*BubblesCarried(SegNumI, nRelArr)
                        end if

                    end if
                end do
            end do
        end do
    end do

    return
end subroutine CEMABubblesTransport


subroutine CEMABubblesRelease()

    use MAIN
    use GLOBAL
    use GEOMC
    use CEMAVars
    use Screenc
!IMPLICIT NONE    
    integer :: nRelArr
    real(8) :: TempBubblesRelVolume


    BRRateAGasNet = 0.d00

    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            do SegNumI = IU, ID
                do nRelArr = 1, NumBubRelArr
                    if (BubblesStatus(SegNumI, nRelArr) == 1 .and. BubblesAtSurface(SegNumI, nRelArr) .and. .not. ICE(SegNumI)) then

                        do nGas = 1, NumGas
!TempBubblesRelVolume = 4/3*3.14*BubblesRadius(SegNumI, nRelArr)**3
                            TempBubblesRelVolume = 4./3.*3.14*BubblesRadius(SegNumI, nRelArr)**3 ! SW 10/10/2017
                            BRVoluAGas(SegNumI, nRelArr, nGas) = BubblesReleaseAllValue(SegNumI, nRelArr)*TempBubblesRelVolume*BubblesGasConc(SegNumI, nRelArr, nGas) !gm
                            BRRateAGas(SegNumI, nRelArr, nGas) = BRVoluAGas(SegNumI, nRelArr, nGas)/dlt !gm/s
                            BRRateAGasNet(SegNumI, nGas) = BRRateAGasNet(SegNumI, nGas) + BRRateAGas(SegNumI, nRelArr, nGas) !gm/s
                            BubbleRelWB(JW, nGas) = BubbleRelWB(JW, nGas) + DLT*BRRateAGasNet(SegNumI, nGas)/1000. ! SW 7/1/2017 Convert from gm/s to kg
                        end do !nGas  

                        BubblesCarried(SegNumI, nRelArr) = BubblesCarried(SegNumI, nRelArr) - BubblesReleaseAllValue(SegNumI, nRelArr)
                        if (BubblesCarried(SegNumI, nRelArr) <= 0.d00) then
                            BubblesReleaseAllValue(SegNumI, nRelArr) = 0.d00
                            BubblesCarried(SegNumI, nRelArr) = 0.d00
                            BubblesGasConc(SegNumI, nRelArr, :) = 0.d00
                            BubblesAtSurface(SegNumI, nRelArr) = .false.
                            BubblesStatus(SegNumI, nRelArr) = 0
                        end if

                    end if
                end do
            end do
        end do
    end do

    return
end subroutine CEMABubblesRelease



subroutine CEMABubbWatTransfer()

    use MAIN
    use GLOBAL
    use GEOMC
    use SCREENC
    use KINETIC
    use CEMAVars
    use CEMASedimentDiagenesis, only: SD_T1
    implicit none
    real(8), allocatable, dimension(:) :: KValue, Mw, Henry
    real(8) :: BubbDissSrcSnk, BubSedT, EqbDissConcentration
    integer :: BubbLNumber, nRelArr, ngasconst, ngas

    allocate(Henry(NumGas), KValue(NumGas), Mw(NumGas))

    Henry(1) = HenryConst_H2S !L atm/M  H2S
    Henry(2) = HenryConst_CH4 !L atm/M  CH4
    Henry(3) = HenryConst_NH3 !L atm/M  NH3
    Henry(4) = HenryConst_CO2 !L atm/M  CO2
    BubSedT = 273.15 + SD_T1 !K
    Mw(1) = 36. !H2S gm/mol
    Mw(2) = 16. !CH4 gm/mol
    Mw(3) = 17. !NH3 gm/mol
    Mw(4) = 44. !CO2 gm/mol

    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            do SegNumI = IU, ID
                do nRelArr = 1, NumBubRelArr
                    if (BubblesStatus(SegNumI, nRelArr) == 1) then

                        do nGas = 1, NumGas
!                        Do nGas = 1, NumGas-1 ! debug

                            if (ngas == 1) then
                                ngasconst = NH2S
                            end if ! cb 2/18/13
                            if (ngas == 2) then
                                ngasconst = NCH4
                            end if
                            if (ngas == 3) then
                                ngasconst = NSO4
                            end if
                            if (ngas == 4) then
                                ngasconst = NTIC
                            end if

                            BubbLNumber = BubblesLNumber(SegNumI, nRelArr)
!KValue(1)        = Henry(1)/GasConst_R/(T1(BubbLNumber,SegNumI) + 273.15)
!KValue(2)        = Henry(2)/GasConst_R/(T1(BubbLNumber,SegNumI) + 273.15)
!KValue(3)        = Henry(3)/GasConst_R/(T1(BubbLNumber,SegNumI) + 273.15)
!KValue(4)        = Henry(4)/GasConst_R/(T1(BubbLNumber,SegNumI) + 273.15)
                            KValue(ngas) = Henry(ngas)/GasConst_R/(T1(BubbLNumber, SegNumI) + 273.15)

                            EqbDissConcentration = BubblesGasConc(SegNumI, nRelArr, nGas)/KValue(nGas)
!BubbDissSrcSnk = BubbWatGasExchRate*(EqbDissConcentration - C1(BubbLNumber,SegNumI,1+nGas))    !g/m�/s
                            BubbDissSrcSnk = BubbWatGasExchRate*(EqbDissConcentration - C1(BubbLNumber, SegNumI, ngasconst)) !g/m�/s  cb 2/18/13
!CGSS(BubbLNumber,SegNumI,nGasconst) = CGSS(BubbLNumber,SegNumI,nGasconst) + BubbDissSrcSnk    !BubbDissSrcSnk > 0 Bubbles --> Water
                            C1(BubbLNumber, SegNumI, nGasconst) = C1(BubbLNumber, SegNumI, nGasconst) + BubbDissSrcSnk !BubbDissSrcSnk > 0 Bubbles --> Water
                            BubblesGasConc(SegNumI, nRelArr, nGas) = BubblesGasConc(SegNumI, nRelArr, nGas) - BubbDissSrcSnk*dlt !BubbDissSrcSnk > 0 Bubbles --> Water
                            if (BubblesGasConc(SegNumI, nRelArr, nGas) < 0.d0) then
                                BubblesGasConc(SegNumI, nRelArr, nGas) = 0.d0
                            end if

                        end do !nGas  

                    end if
                end do
            end do
        end do
    end do

    return
end subroutine CEMABubbWatTransfer


subroutine CEMABubblesTurbulence()
    use MAIN
    use GLOBAL
    use GEOMC
    use SCREENC
    use KINETIC
    use CEMAVars
    implicit none
    real(8) :: TempBubbDiam, TempRelVelocity
    integer :: BubbCntr, nRelArr, BubbLNumber

    SegNumI = I
    do k = KT, KBMIN(SegNumI) - 1
        BubbCntr = 0
        TempBubbDiam = 0
        TempRelVelocity = 0
        do nRelArr = 1, NumBubRelArr
            if (BubblesStatus(SegNumI, nRelArr) == 1) then
                BubbLNumber = BubblesLNumber(SegNumI, nRelArr)

                if (BubbLNumber == k) then
                    BubbCntr = BubbCntr + 1
                    TempBubbDiam = TempBubbDiam + 2.0*BubblesRadius(SegNumI, nRelArr)
                    TempRelVelocity = TempRelVelocity + BubblesRiseV(SegNumI, nRelArr) + W(K - 1, SegNumI) !Rise velocity is +ve upwards and W is +ve downwards
                end if

            end if
        end do !nRelArr

        if (BubbCntr > 0) then
            TempBubbDiam = TempBubbDiam/BubbCntr
            TempRelVelocity = TempRelVelocity/BubbCntr
            if (k == KB(SegNumI)) then
                AZ(K, SegNumI) = AZ(K, SegNumI) + BottomTurbulence(SegNumI)
            else
!AZ(K,SegNumI) = AZ(K,SegNumI)  + TempBubbDiam/TempRelVelocity
                AZ(K, SegNumI) = AZ(K, SegNumI) + TempBubbDiam*TempRelVelocity ! cb 2/7/13
            end if
        end if

        AZ(KB(SegNumI) - 1, SegNumI) = AZ(KB(SegNumI) - 1, SegNumI) + BottomTurbulence(SegNumI)

    end do

    return
end subroutine CEMABubblesTurbulence


subroutine CEMABubblesReleaseTurbulence()
    use MAIN
    use GLOBAL
    use GEOMC
    use SCREENC
    use KINETIC
    use CEMAVars
    implicit none
    real(8) :: TempBubbDiam, TempRelVelocity
    integer :: BubbCntr, nRelArr, BubbLNumber

    BottomTurbulence = 0.d00
    do JW = 1, NWB
        k = KB(SegNumI)
        do JB = BS(JW), BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            do SegNumI = IU, ID

                BubbCntr = 0
                TempBubbDiam = 0
                TempRelVelocity = 0
                do nRelArr = 1, NumBubRelArr
                    if (BubblesStatus(SegNumI, nRelArr) == 1) then
                        BubbLNumber = BubblesLNumber(SegNumI, nRelArr)

                        if (BubbLNumber == k) then
                            BubbCntr = BubbCntr + 1
                            TempBubbDiam = TempBubbDiam + 2.0*BubblesRadius(SegNumI, nRelArr)
                            TempRelVelocity = TempRelVelocity + BubblesRiseV(SegNumI, nRelArr) + W(K - 1, SegNumI) !Rise velocity is +ve upwards and W is +ve downwards
                        end if

                    end if
                end do !nRelArr

                if (BubbCntr > 0) then
                    TempBubbDiam = TempBubbDiam/BubbCntr
                    TempRelVelocity = TempRelVelocity/BubbCntr
                    BottomTurbulence(SegNumI) = CEMATurbulenceScaling*TempBubbDiam/TempRelVelocity
                end if

            end do !SegNumI
        end do !JB       
    end do !JW

    do JW = 1, NWB
        k = KB(SegNumI)
        do JB = BS(JW), BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            do SegNumI = IU, ID

                if (BottomTurbulence(SegNumI) > 0) then
                    continue
                end if

            end do
        end do
    end do

    return
end subroutine CEMABubblesReleaseTurbulence
