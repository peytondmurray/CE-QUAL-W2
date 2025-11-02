subroutine SHADING()
    use SHADEC;     use GLOBAL;     use GDAYC;     use SURFHE;     use GEOMC;     use SCREENC;     use LOGICC
    implicit none
    character(len=1) :: BANK
    real :: LOCAL, STANDARD, HOUR, TAUD, SINAL, A02, AZ00, A0, AX, ANG1, ANG2, TOPOANG, SFACT, HT, CLINE, SRED, STLEN, EDGE, EDAZ, SN, AZT
    integer :: IDAY, J

! Calculate solar altitude, declination, and local hour angle when short-wave solar radiation is provided as input

    if (READ_RADIATION(JW)) then
        LOCAL = LONGIT(JW)
        STANDARD = 15.0*INT(LONGIT(JW)/15.0)
        HOUR = (JDAY - INT(JDAY))*24.0
        IDAY = JDAY - INT(JDAY/365)*365
!IDAY     =  IDAY+INT(INT(JDAY/365)/4)    SR 12/2018 IDAY FIX
        IDAY = IDAY - INT(INT(JDAY/365)/4)
        TAUD = 2*PI*(IDAY - 1)/365
        EQTNEW = 0.170*SIN(4*PI*(IDAY - 80)/373) - 0.129*SIN(2*PI*(IDAY - 8)/355)
        HH(JW) = 0.261799*(HOUR - (LOCAL - STANDARD)*0.0666667 + EQTNEW - 12.0)
        DECL(JW) = 0.006918 - 0.399912*COS(TAUD) + 0.070257*SIN(TAUD) - 0.006758*COS(2*TAUD) + 0.000907*SIN(2*TAUD) - 0.002697*COS(3*TAUD) + 0.001480*SIN(3*TAUD)
        SINAL = SIN(LAT(JW)*.0174533)*SIN(DECL(JW)) + COS(LAT(JW)*.0174533)*COS(DECL(JW))*COS(HH(JW))
        A00(JW) = 57.2957795*ASIN(SINAL)
    end if

! If the sun is below the horizon, set SHADE(I) to 0

    if (A00(JW) < 0.0) then
        SHADE(I) = 0.0
    else

!** Calculate solar azimuth angle

        A02 = A00(JW)/57.2957795
        AX = (SIN(DECL(JW))*COS(LAT(JW)*0.017453) - COS(DECL(JW))*COS(HH(JW))*SIN(LAT(JW)*0.017453))/COS(A02)
        if (AX > 1.0) then
            AX = 1.0
        end if
        if (AX < -1.0) then
            AX = -1.0
        end if
        AZT = ACOS(AX)
        if (HH(JW) < 0.0) then
            AZ00 = AZT
        else
            AZ00 = 2.0*PI - AZT
        end if
        A0 = A02

!** Interpolate the topographic shade angle

        do J = 1, IANG - 1
            if (AZ00 > ANG(J) .and. AZ00 <= ANG(J + 1)) then
                ANG1 = AZ00 - ANG(J)
                ANG2 = (TOPO(I, J + 1) - TOPO(I, J))/GAMA ! SW 10/17/05
                TOPOANG = TOPO(I, J) + ANG2*ANG1
            end if
        end do
        if (AZ00 > ANG(IANG) .and. AZ00 <= 2*PI) then
            ANG1 = AZ00 - ANG(IANG)
            ANG2 = (TOPO(I, 1) - TOPO(I, IANG))/GAMA ! SW 10/17/05
            TOPOANG = TOPO(I, IANG) + ANG2*ANG1
        end if

!** Complete topographic shading if solar altitude less than topo angle

        if (A0 <= TOPOANG) then
            SFACT = 0.90
            go to 100
        end if

!** No vegetative shading if azimuth angle is oriented parallel to stream

        if (AZ00 == PHI0(I) .or. AZ00 == PHI0(I) + PI .or. AZ00 + PI == PHI0(I)) then
            SFACT = 0.0
            go to 100
        end if

!** Bank with the controlling vegetation

        if (PHI0(I) > 0.0 .and. PHI0(I) <= PI) then
            if (AZ00 > PHI0(I) .and. AZ00 <= PHI0(I) + PI) then
                BANK = "L"
            end if
            if (AZ00 > 0.0 .and. AZ00 <= PHI0(I)) then
                BANK = "R"
            end if
            if (AZ00 > PHI0(I) + PI .and. AZ00 < 2.0*PI) then
                BANK = "R"
            end if
        else
            if (PHI0(I) > PI .and. PHI0(I) <= 2.0*PI) then
                if (AZ00 >= PHI0(I) .and. AZ00 < 2.0*PI) then
                    BANK = "L"
                end if
                if (AZ00 >= 0.0 .and. AZ00 < PHI0(I) - PI) then
                    BANK = "L"
                end if
                if (AZ00 >= PHI0(I) - PI .and. AZ00 < PHI0(I)) then
                    BANK = "R"
                end if
            end if
        end if

!** No topographic shading

        if (BANK == "L") then
            if (TTLB(I) < ELWS(I)) then
                SFACT = 0.0
                go to 100
            else
                HT = TTLB(I) - ELWS(I)
                CLINE = CLLB(I)
                SRED = SRLB2(I)
                if (JDAYG > SRFJD1(I) .and. JDAYG <= SRFJD2(I)) then
                    SRED = SRLB1(I)
                end if
            end if
        else
            if (TTRB(I) < ELWS(I)) then
                SFACT = 0.0
                go to 100
            else
                HT = TTRB(I) - ELWS(I)
                CLINE = CLRB(I)
                SRED = SRRB2(I)
                if (JDAYG > SRFJD1(I) .and. JDAYG <= SRFJD2(I)) then
                    SRED = SRRB1(I)
                end if
            end if
        end if
        STLEN = HT/TAN(A0)
        EDGE = MAX(0.0, CLINE - BI(KT, I)/2.0)

!** Distance from vegetation to water edge on line parallel to azimuth

        EDAZ = EDGE/ABS(SIN(PHI0(I) - AZ00))
        if (STLEN <= EDAZ) then
            SFACT = 0.0
            go to 100
        end if

!** Distance shadow extends over water (perpendicular to segment orientation)

        SN = MIN(HT*ABS(SIN(ABS(PHI0(I) - AZ00)))/TAN(A0) - EDGE, BI(KT, I))
        SFACT = SRED*SN/BI(KT, I)
        100 continue
        SHADE(I) = MAX(0.0, 1.0 - SFACT)
        SHADE(I) = MIN(ABS(SHADEI(I)), SHADE(I)) ! SW 10/2/2017 Allows for fixed canopy cover over top of channel - only used if shade is less than shadei only valid for -0.99 and 0.0
    end if
    return
end subroutine SHADING
