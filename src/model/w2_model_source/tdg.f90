subroutine TOTAL_DISSOLVED_GAS(NSAT, P, NSG, N, T, C)
    use TDGAS;     use STRUCTURES;     use GLOBAL;     use MAIN, only: EA;     use TVDC, only: TDEW
    implicit none
    integer :: N, NSG, NSAT
    real(R8) :: T, P, C
    real :: SAT, DB, DA, TDG

    if (NSAT == 0) then ! DISSOLVED OXYGEN
        SAT = EXP(7.7117 - 1.31403*LOG(T + 45.93))*P
    else
        if (NSAT == 1) then ! N2 GAS
            EA = DEXP(2.3026D0*(7.5D0*TDEW(JW)/(TDEW(JW) + 237.3D0) + 0.6609D0))*0.001316 ! in mm Hg   0.0098692atm=7.5006151mmHg     
            SAT = 1.5568D06*0.79*(P - EA)*(1.8816D-5 - 4.116D-7*T + 4.6D-9*T*T) ! SW 10/27/15    4/20/16 SPEED 
        else
            SAT = p
        end if
    end if

    if (NSG == 0) then
        if (EQSP(N) == 1) then
            TDG = AGASSP(N)*.035313*QSP(N) + BGASSP(N)
            if (TDG > 145.0) then
                TDG = 145.0
            end if
            C = SAT
            if (TDG >= 100.0) then
                C = TDG*SAT/100.0
            end if
        else
            if (EQSP(N) == 2) then
                TDG = AGASSP(N) + BGASSP(N)*EXP(0.035313*QSP(N)*CGASSP(N))
                if (TDG > 145.0) then
                    TDG = 145.0
                end if
                C = SAT
                if (TDG >= 100.0) then
                    C = TDG*SAT/100.0
                end if
            else
                if (EQSP(N) == 3) then
                    DA = SAT - C ! MM 5/21/2009 DA: Deficit upstream
                    DB = DA/(1.0 + 0.38*AGASSP(N)*BGASSP(N)*CGASSP(N)*(1.0 - 0.11*CGASSP(N))*(1.0 + 0.046*T)) ! DB: deficit downstream
                    C = SAT - DB
                else
                    if (EQSP(N) == 4) then
                        C = AGASSP(N)*C
                        if (CGASSP(N) == 1 .and. C > SAT) then
                            C = SAT
                        end if
                    end if
                end if
            end if
        end if
    else
        if (EQGT(N) == 1) then
            TDG = AGASGT(N)*0.035313*QGT(N) + BGASGT(N)
            if (TDG > 145.0) then
                TDG = 145.0
            end if
            C = SAT
            if (TDG >= 100.0) then
                C = TDG*SAT/100.0
            end if
        else
            if (EQGT(N) == 2) then
                TDG = AGASGT(N) + BGASGT(N)*EXP(.035313*QGT(N)*CGASGT(N))
                if (TDG > 145.0) then
                    TDG = 145.0
                end if
                C = SAT
                if (TDG >= 100.0) then
                    C = TDG*SAT/100.0
                end if
            else
                if (EQGT(N) == 3) then
                    DA = SAT - C ! MM 5/21/2009 DA: Deficit upstream
                    DB = DA/(1.0 + 0.38*AGASGT(N)*BGASGT(N)*CGASGT(N)*(1.0 - 0.11*CGASGT(N))*(1.0 + 0.046*T)) ! DB: deficit downstream
                    C = SAT - DB
                else
                    if (EQGT(N) == 4) then
                        C = AGASGT(N)*C
                        if (CGASGT(N) == 1 .and. C > SAT) then
                            C = SAT
                        end if
                    end if
                end if
            end if
        end if
    end if
end subroutine TOTAL_DISSOLVED_GAS
