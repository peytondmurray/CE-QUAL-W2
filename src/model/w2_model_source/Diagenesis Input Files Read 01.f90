subroutine InitializeBedConsolidationFiles(TempFilNum, TempFilName)
    use MAIN
    use GLOBAL
    use SCREENC
    use CEMAVars

    implicit none

    logical :: SkipLoop
    integer(4) :: TempFilNum
    character(len=256) :: MessageTemp, TempFilName

!Open File	
    open(TempFilNum, File=TempFilName(1:len_trim(TempFilName) - 1))

!Read Header
    SkipLoop = .false.
    do while (.not. SkipLoop)
        read(TempFilNum, "(a)") MessageTemp
        if (index(MessageTemp, "$") == 0) then
            SkipLoop = .true.
        end if
    end do

end subroutine InitializeBedConsolidationFiles


subroutine ReadBedConsolidationFiles(TempFilNum)
    use MAIN
    use GLOBAL
    use SCREENC
    use CEMAVars

    logical :: SkipLoop
    integer(4) :: TempFilNum
    real(8) :: TimeJD1
    real(8) :: TimeJD2
    real(8) :: FactorInterp
    real(8) :: ConsolidRateTemp11(NumConsolidRegns), ConsolidRateTemp1(NumConsolidRegns), ConsolidRateTemp2(NumConsolidRegns)

!Read Data
    SkipLoop = .false.
    do while (.not. SkipLoop .or. EOF(TempFilNum))
        read(TempFilNum, "(F8.0,<NumConsolidRegns>F8.0)") TimeJD1, (ConsolidRateTemp1(i), i = 1, NumConsolidRegns)
        read(TempFilNum, "(F8.0,<NumConsolidRegns>F8.0)") TimeJD2, (ConsolidRateTemp2(i), i = 1, NumConsolidRegns)

        if (JDay >= TimeJD1 .and. JDay <= TimeJD2) then
            SkipLoop = .true.
        end if
        backspace(TempFilNum)
    end do

    backspace(TempFilNum)

    FactorInterp = (JDay - TimeJD1)/(TimeJD2 - TimeJD1)

    do i = 1, NumConsolidRegns
        ConsolidRateTemp11(i) = ConsolidRateTemp1(i)*(1 - FactorInterp) + ConsolidRateTemp2(i)*FactorInterp
    end do

end subroutine ReadBedConsolidationFiles
