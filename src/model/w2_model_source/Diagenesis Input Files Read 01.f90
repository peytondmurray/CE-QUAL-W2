Subroutine InitializeBedConsolidationFiles(TempFilNum, TempFilName)
    Use MAIN
    Use GLOBAL
    Use SCREENC
    Use CEMAVars

    Implicit None

    Logical SkipLoop
    Integer(4) TempFilNum
    Character(256) MessageTemp, TempFilName

    !Open File
	Open(TempFilNum, File = TempFilName(1:len_trim(TempFilName)-1))

	!Read Header
	SkipLoop = .FALSE.
	Do While(.NOT. SkipLoop)
		Read(TempFilNum,'(a)')MessageTemp
		If(index(MessageTemp, "$") == 0)SkipLoop = .TRUE.
	End Do

End Subroutine

Subroutine ReadBedConsolidationFiles(TempFilNum)
    Use MAIN
    Use GLOBAL
    Use SCREENC
    Use CEMAVars

    Logical    :: SkipLoop
    Integer :: ios
    Integer(4) :: TempFilNum
    Real(8) :: TimeJD1
    Real(8) :: TimeJD2
    Real(8) :: FactorInterp
	Real(8) :: ConsolidRateTemp11(NumConsolidRegns),ConsolidRateTemp1(NumConsolidRegns),ConsolidRateTemp2(NumConsolidRegns)
    Character(:), allocatable :: fmt

    !Read Data
	SkipLoop = .FALSE.
    fmt = '(F8.0'//Repeat(',F8.0', NumConsolidRegns)//')'

	Do While(.NOT. SkipLoop .or. EOF(TempFilNum))
        Read(TempFilNum, fmt, ios=ios) TimeJD1, (ConsolidRateTemp1(i), i=1, NumConsolidRegns)
        if (is_iostat_end(ios)) exit
        if (ios /= 0) stop "Error reading TimeJD1 from TempFilNum"

        Read(TempFilNum, fmt, ios=ios) TimeJD2, (ConsolidRateTemp2(i), i=1, NumConsolidRegns)
        if (is_iostat_end(ios)) exit
        if (ios /= 0) stop "Error reading TimeJD2 from TempFilNum"

	    If(JDay >= TimeJD1 .and. JDay <= TimeJD2)then
		    SkipLoop = .TRUE.
	    End If
	    BackSpace(TempFilNum)
	End Do

	BackSpace(TempFilNum)

	FactorInterp = (JDay - TimeJD1)/(TimeJD2 - TimeJD1)

  DO i=1,NumConsolidRegns
    ConsolidRateTemp11(i) = ConsolidRateTemp1(i)*(1-FactorInterp) + ConsolidRateTemp2(i)*FactorInterp
  END DO

End Subroutine
