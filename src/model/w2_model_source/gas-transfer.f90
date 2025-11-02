subroutine GAS_TRANSFER()
    use GLOBAL;     use GEOMC;     use KINETIC
    implicit none
    real, parameter :: THETA_REAERATION = 1.024, M_TO_FT = 3.2808
    real :: AREA, ADEPTH, UAVG, HDEPTH, S, USTAR, A, BCOEF, DMO2
    integer :: K

    if (REAERC(JW) == "   RIVER") then

!** Average depth in ft

        AREA = 0.0
        do K = KT, KBMIN(I)
            AREA = AREA + BHR1(K, I)
        end do
        ADEPTH = AREA/BR(KTI(I), I)*M_TO_FT

!** Average velocity in feet/second

        UAVG = ABS(QC(I))/AREA*M_TO_FT

!** Reaeration factor

        if (NEQN(JW) == 0) then
            if (ADEPTH <= 2.0) then
                REAER(I) = 21.64*UAVG**0.67/ADEPTH**1.85
            else
                if (UAVG <= 1.8) then
                    REAER(I) = 12.96*SQRT(UAVG)/ADEPTH**1.5
                else
                    HDEPTH = (-11.875)*UAVG + 23.375
                    if (HDEPTH >= ADEPTH) then
                        REAER(I) = 12.96*SQRT(UAVG)/ADEPTH**1.5
                    else
                        REAER(I) = 11.57*UAVG**0.969/ADEPTH**1.673
                    end if
                end if
            end if
        else
            if (NEQN(JW) == 1) then !O'connor-Dobbins
                REAER(I) = 12.96*SQRT(UAVG)/ADEPTH**1.5 ! units: day-1
            else
                if (NEQN(JW) == 2) then !Churchill
                    REAER(I) = 11.57*UAVG**0.969/ADEPTH**1.673 ! units: day-1
                else
                    if (NEQN(JW) == 3) then !Tsivoglou
                        S = SLOPEC(JB)*5280.0
                        if (ABS(QC(I))*35.5 >= 10.0) then
                            REAER(I) = 0.88*S*UAVG
                        else
                            REAER(I) = 1.8*S*UAVG
                        end if
                    else
                        if (NEQN(JW) == 4) then !Owens
                            REAER(I) = 21.64*UAVG**0.67/ADEPTH**1.85
                        else
                            if (NEQN(JW) == 5) then !Thackston and Krenkel
                                USTAR = SQRT(ADEPTH*SLOPEC(JB)*32.2) ! SR 5/10/05
                                REAER(I) = 24.88*(1.0 + SQRT(0.176*UAVG/SQRT(ADEPTH)))*USTAR/ADEPTH ! SR 5/10/05
                            else
                                if (NEQN(JW) == 6) then !Langbien and Durum
                                    REAER(I) = 7.60*UAVG/ADEPTH**1.33
                                else
                                    if (NEQN(JW) == 7) then !Melching and Flores
                                        UAVG = UAVG/M_TO_FT
                                        if (QC(I) == 0.0) then
                                            REAER(I) = 0.0
                                        else
                                            if (ABS(QC(I)) < 0.556) then
                                                REAER(I) = 517.0*(UAVG*SLOPEC(JB))**0.524*ABS(QC(I))**(-0.242)
                                            else
                                                REAER(I) = 596.0*(UAVG*SLOPEC(JB))**0.528*ABS(QC(I))**(-0.136)
                                            end if
                                        end if
                                    else
                                        if (NEQN(JW) == 8) then !Melching and Flores
                                            UAVG = UAVG/M_TO_FT
                                            ADEPTH = ADEPTH/M_TO_FT
                                            if (ABS(QC(I)) < 0.556) then
                                                REAER(I) = 88.0*(UAVG*SLOPEC(JB))**0.313*ADEPTH**(-0.353)
                                            else
                                                REAER(I) = 142.0*(UAVG*SLOPEC(JB))**0.333*ADEPTH**(-0.66)*BI(KT, I)**(-0.243)
                                            end if
                                        else
                                            if (NEQN(JW) == 9) then !User defined SI units
                                                UAVG = UAVG/M_TO_FT
                                                ADEPTH = ADEPTH/M_TO_FT
                                                REAER(I) = RCOEF1(JW)*UAVG**RCOEF2(JW)*ADEPTH**RCOEF3(JW)*SLOPEC(JB)**RCOEF4(JW)
                                            else
                                                if (NEQN(JW) == 10) then ! Thackston and Krenkel - updated
                                                    USTAR = SQRT(ADEPTH*SLOPEC(JB)*32.2) ! SR 5/10/05
                                                    REAER(I) = 4.99*(1.0 + 9.0*(0.176*UAVG/SQRT(ADEPTH))**0.25)*USTAR/ADEPTH ! SR 5/10/05
                                                end if
                                            end if
                                        end if
                                    end if
                                end if
                            end if
                        end if
                    end if
                end if
            end if
        end if
!REAER(I) = REAER(I)*ADEPTH/M_TO_FT - now eliminated since changed computation in DO water quality section
        if (MINKL(JW) > REAER(I)) then
            REAER(I) = MINKL(JW)
        end if ! minimum value units day-1
    else
        if (REAERC(JW) == "    LAKE") then
            if (NEQN(JW) == 1) then !Broecker
                REAER(I) = 0.864*WIND10(I)
            else
                if (NEQN(JW) == 2) then
                    if (WIND10(I) <= 3.5) then !Gelda
                        A = 0.2
                        BCOEF = 1.0
                    else
                        A = 0.057
                        BCOEF = 2.0
                    end if
                    REAER(I) = A*WIND10(I)**BCOEF
                else
                    if (NEQN(JW) == 3) then !Banks & Herrera
                        REAER(I) = 0.728*SQRT(WIND10(I)) - 0.317*WIND10(I) + 0.0372*WIND10(I)**2 !units of m/day
                    else
                        if (NEQN(JW) == 4) then !Wanninkhof
                            REAER(I) = 0.0986*WIND10(I)**1.64 !units of m/day
                        else
                            if (NEQN(JW) == 5) then !Chen & Kanwisher
                                DMO2 = 2.04E-9
                                REAER(I) = DAY*DMO2/((200.0 - 60.0*SQRT(MIN(WIND10(I), 11.0)))*1.E-6)
                            else
                                if (NEQN(JW) == 6) then !Cole & Buchak
                                    REAER(I) = 0.5 + 0.05*WIND10(I)*WIND10(I)
                                else
                                    if (NEQN(JW) == 7) then !Banks
                                        if (WIND10(I) <= 5.5) then
                                            REAER(I) = 0.362*SQRT(WIND10(I))
                                        else
                                            REAER(I) = 0.0277*WIND10(I)**2
                                        end if
                                    else
                                        if (NEQN(JW) == 8) then !Smith
                                            REAER(I) = 0.64 + 0.128*WIND10(I)**2
                                        else
                                            if (NEQN(JW) == 9) then !Liss
                                                if (WIND10(I) <= 4.1) then
                                                    REAER(I) = 0.156*WIND10(I)**0.63
                                                else
                                                    REAER(I) = 0.0269*WIND10(I)**1.9
                                                end if
                                            else
                                                if (NEQN(JW) == 10) then !Downing and Truesdale
                                                    REAER(I) = 0.0276*WIND10(I)**2
                                                else
                                                    if (NEQN(JW) == 11) then !Kanwisher
                                                        REAER(I) = 0.0432*WIND10(I)**2
                                                    else
                                                        if (NEQN(JW) == 12) then !Yu, et al
                                                            REAER(I) = 0.319*WIND10(I)
                                                        else
                                                            if (NEQN(JW) == 13) then !Weiler
                                                                if (WIND10(I) <= 1.6) then
                                                                    REAER(I) = 0.398
                                                                else
                                                                    REAER(I) = 0.155*WIND10(I)**2
                                                                end if
                                                            else
                                                                if (NEQN(JW) == 14) then !User defined
                                                                    REAER(I) = RCOEF1(JW) + RCOEF2(JW)*WIND10(I)**RCOEF3(JW)
                                                                end if
                                                            end if
                                                        end if
                                                    end if
                                                end if
                                            end if
                                        end if
                                    end if
                                end if
                            end if
                        end if
                    end if
                end if
            end if
            if (MINKL(JW) > REAER(I)) then
                REAER(I) = MINKL(JW)
            end if ! minimum value units m day-1
            REAER(I) = REAER(I)*BI(KT, I)/BH2(KT, I) ! conversion from m/d to 1/d for all wind based equations
        else
            if (REAERC(JW) == " ESTUARY") then
                AREA = 0.0
                do K = KT, KBMIN(I)
                    AREA = AREA + BHR1(K, I)
                end do

                ADEPTH = AREA/BR(KTI(I), I)
                UAVG = ABS(QC(I))/AREA

!** Reaeration factor

                if (NEQN(JW) == 0) then
                    ADEPTH = ADEPTH*M_TO_FT
                    UAVG = UAVG*M_TO_FT
                    if (ADEPTH <= 2.0) then
                        REAER(I) = 21.64*UAVG**0.67/ADEPTH**1.85 ! units of day-1
                    else
                        if (UAVG <= 1.8) then
                            REAER(I) = 12.96*SQRT(UAVG)/ADEPTH**1.5
                        else
                            HDEPTH = (-11.875)*UAVG + 23.375
                            if (HDEPTH >= ADEPTH) then
                                REAER(I) = 12.96*SQRT(UAVG)/ADEPTH**1.5
                            else
                                REAER(I) = 11.57*UAVG**0.969/ADEPTH**1.673
                            end if
                        end if
                    end if
                else
                    if (NEQN(JW) == 1) then !Thomann and Fitzpatrick
                        REAER(I) = 0.728*SQRT(WIND10(I)) - 0.317*WIND10(I) + 0.0372*WIND10(I)**2 + 3.93*SQRT(UAVG)/(ADEPTH)**0.5 ! units of m/day
                        REAER(I) = REAER(I)/ADEPTH ! units of 1/day
                    else
                        if (NEQN(JW) == 2) then !User defined
                            REAER(I) = RCOEF1(JW)*UAVG**RCOEF2(JW)*ADEPTH**RCOEF3(JW) + 0.5 + RCOEF4(JW)*WIND10(I)*WIND10(I) ! units of 1/day
                        end if
                    end if
                end if
                if (MINKL(JW) > REAER(I)) then
                    REAER(I) = MINKL(JW)
                end if ! minimum value units day-1

            end if
        end if
    end if
!IF(RCOEF1(JW)>REAER(I))REAER(I)=RCOEF1(JW)    ! minimum value units day-1
    REAER(I) = REAER(I)*THETA_REAERATION**(T1(KT, I) - 20.0) ! units of day-1
    REAER(I) = REAER(I)/DAY
end subroutine GAS_TRANSFER
