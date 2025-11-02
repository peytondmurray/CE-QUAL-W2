subroutine CEMAFFTLayerCode()

    use MAIN
    use GLOBAL
    use GEOMC
    use SCREENC
    use RSTART
    use PREC
    use EDDY
    use LOGICC
    use TVDC
    use KINETIC
    use CEMAVars

    implicit none

    integer :: iTemp

    NMFT = 0
    do JG = 1, NSS
        if (SSCS(JG) == -1.0) then
            NMFT = NSSS + JG - 1
            exit
        end if
    end do

    if (FirstTimeInFFTCode) then
        FirstTimeInFFTCode = .false.
        FFTActive = .true.
        do JW = 1, NWB
            do JB = BS(JW), BE(JW)
                do I = CUS(JB), DS(JB)
!JAC = NSSS
                    JAC = NMFT ! cb 2/18/13
                    K = KB(I)
                    C1(K, I, JAC) = InitFFTLayerConc
                    C1S(K, I, JAC) = InitFFTLayerConc
                    C2(K, I, JAC) = InitFFTLayerConc
                end do
            end do
        end do
        return
    end if


    do iTemp = 1, NumFFTActivePrds

        if (.not. FFTActive) then
            if (JDAY > FFTActPrdSt(iTemp) .and. JDAY < FFTActPrdEn(iTemp)) then

                FFTActive = .true.
                do JW = 1, NWB
                    do JB = BS(JW), BE(JW)
                        do I = CUS(JB), DS(JB)
!JAC = NSSS
                            JAC = NMFT ! cb 2/18/13
                            K = KB(I)
                            C1(K, I, JAC) = FFTLayConc(I)
                            C1S(K, I, JAC) = FFTLayConc(I)
                            C2(K, I, JAC) = FFTLayConc(I)
                        end do
                    end do
                end do
                FFTActPrd = iTemp
                exit
            end if
        end if

    end do

    if (FFTActive) then
        iTemp = FFTActPrd
        if (JDAY > FFTActPrdEn(iTemp)) then
            FFTActive = .false.
            do JW = 1, NWB
                do JB = BS(JW), BE(JW)
                    do I = CUS(JB), DS(JB)
!JAC = NSSS
                        JAC = NMFT ! cb 2/18/13
                        K = KB(I)
                        FFTLayConc(I) = C1(K, I, JAC)
                        C1(K, I, JAC) = 0.d00
                        C1S(K, I, JAC) = 0.d00
                        C2(K, I, JAC) = 0.d00
                    end do
                end do
            end do
        end if
    end if

    return

    entry MoveFFTLayerConsolid()

    do JW = 1, NWB
        KT = KTWB(JW)
        do JB = BS(JW), BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            do SegNumI = IU, ID
                if (CEMALayerAdded(SegNumI)) then
!JAC = NSSS
                    jac = NMFT ! cb 2/18/13
                    do K = KB(SegNumI) - 1, KT, -1
                        C1(K + 1, SegNumI, JAC) = C1(K, SegNumI, JAC)
                        C1S(K + 1, SegNumI, JAC) = C1S(K, SegNumI, JAC)
                        C2(K + 1, SegNumI, JAC) = C2(K, SegNumI, JAC)
                    end do !K
                    C1(KT, SegNumI, JAC) = 0.d00
                    C1S(KT, SegNumI, JAC) = 0.d00
                    C2(KT, SegNumI, JAC) = 0.d00
                end if
            end do !SegNumI
        end do !JB
    end do !JW

    return
end subroutine CEMAFFTLayerCode
