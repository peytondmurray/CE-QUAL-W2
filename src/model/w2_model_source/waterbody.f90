subroutine WATERBODY()
    use GLOBAL;     use GEOMC;     use TVDC;     use LOGICC;     use PREC

! Type declarations
    implicit none

    real, save, allocatable, dimension(:) :: ELL, ELR, CL
    real, save, allocatable, dimension(:,:) :: QU, QD
    real(R8) :: C(KMX,IMX), SS(KMX,IMX)
    real :: ELW, EL1
    real :: Q1, HT, FRAC, BRTOT, T1L, T2L, U1, B1, B2
    integer :: IUT, IDT, JJW, K, KL, KR

! Allocation declarations

    allocate(ELL(KMX), ELR(KMX), CL(NCT), QU(KMX, IMX), QD(KMX, IMX))

! Variable initialization

    ELL = 0.0;     ELR = 0.0;     CL = 0.0;     QU = 0.0;     QD = 0.0

! Debug variable
!  ncount=0
! End debug

    return

!***********************************************************************************************************************************
!**                                               U P S T R E A M   V E L O C I T Y                                               **
!***********************************************************************************************************************************

    entry UPSTREAM_VELOCITY()
    do JJB = 1, NBR
        if (UHS(JB) >= US(JJB) .and. UHS(JB) <= DS(JJB)) then
            exit
        end if
    end do
    do JJW = 1, NWB
        if (JJB >= BS(JJW) .and. JJB <= BE(JJW)) then
            exit
        end if
    end do
    do K = KTWB(JJW), KB(UHS(JB)) + 1
        ELL(K) = EL(K, UHS(JB)) - SINA(JJB)*DLX(UHS(JB))*0.5
    end do
    do K = KT, KB(IU) + 1
        ELR(K) = EL(K, IU) + SINA(JB)*DLX(IU)*0.5
    end do
    ELW = ELWS(UHS(JB)) !EL(KTWB(JJW),UHS(JB))-Z(UHS(JB))*COSA(JJB)
    EL1 = ELW - SINA(JJB)*DLX(UHS(JB))*0.5
!   IF(SLOPE(JB) /= 0.0)THEN
!   EL1  = ELW-(ELWS(UHS(JB)+1)-ELW)/(0.5*(DLX(UHS(JB)+DLX(JB)+1))*DLX(UHS(JB))*0.5   ! SW 7/17/09
!   ELSE
!   EL1=ELW
!   ENDIF
    ELW = EL1
    KL = KTWB(JJW) + 1
    do K = KT + 1, KB(IU) + 1
        if (ELR(K) >= ELL(KL)) then
            if (KL == KTWB(JJW) + 1) then
                Q1 = U(KL - 1, UHS(JB))*BHR1(KTWB(JJW), UHS(JB))
                if (KL == KB(UHS(JB)) + 1 .and. ELL(KL) < ELR(KB(IU) + 1)) then
                    HT = ELW - ELR(KB(IU) + 1)
                else
                    HT = H1(KTWB(JJW), UHS(JB))
                end if
            else
                Q1 = U(KL - 1, UHS(JB))*BHR1(KL - 1, UHS(JB))
                HT = H1(KL - 1, UHS(JB))
            end if
            if (K == KT + 1) then
                EL1 = EL(KTWB(JJW), UHS(JB)) - Z(UHS(JB))*COSA(JJB) - SINA(JJB)*DLX(UHS(JB))*0.5
                U(K - 1, IU - 1) = Q1*(EL1 - ELR(K))/HT/BHR1(KT, IU - 1)
            else
                U(K - 1, IU - 1) = Q1*(EL1 - ELR(K))/HT/BHR1(K - 1, IU - 1)
            end if
            EL1 = ELR(K)
            if (ELR(K) == ELL(KL)) then
                KL = KL + 1
            end if
        else
            Q1 = 0.0
            do while (ELR(K) <= ELL(KL))
                if (KL == KTWB(JJW) + 1 .and. K == KT + 1) then
                    Q1 = Q1 + U(KL - 1, UHS(JB))*BHR1(KTWB(JJW), UHS(JB))
                else
                    if (KL == KTWB(JJW) + 1) then
                        Q1 = Q1 + U(KL - 1, UHS(JB))*BHR1(KTWB(JJW), UHS(JB))*(EL1 - ELL(KL))/H1(KTWB(JJW), UHS(JB))
                    else
                        FRAC = (EL1 - ELL(KL))/H1(KL - 1, UHS(JB))
                        Q1 = Q1 + U(KL - 1, UHS(JB))*BHR1(KL - 1, UHS(JB))*FRAC
                    end if
                end if
                EL1 = ELL(KL)
                KL = KL + 1
                if (KL > KB(UHS(JB))) then
                    exit
                end if
            end do
            if (K == KT + 1) then
                BRTOT = BHR1(KT, IU - 1)
            else
                BRTOT = BHR1(K - 1, IU - 1)
            end if
            FRAC = 0.0
!      IF (KL < KMX) THEN                             ! SW 6/29/06
            HT = H1(KL - 1, UHS(JB))
            FRAC = (EL1 - ELR(K))/HT
            if (KB(UHS(JB)) >= KB(IU - 1) .and. K > KB(IU - 1)) then
                FRAC = (EL1 - ELL(KL))/HT
            end if ! SW 6/29/06
            if (KL >= KMX .and. ELR(K) < ELL(KL)) then
                FRAC = 1.0
            end if ! SW 6/29/06
            Q1 = Q1 + U(KL - 1, UHS(JB))*BHR1(KL - 1, UHS(JB))*FRAC
!      ELSE                                           ! SW 6/29/06
!        Q1 = Q1+U(KL-1,UHS(JB))*BHR1(KL-1,UHS(JB))   ! SW 6/29/06
!      END IF                                         ! SW 6/29/06
            U(K - 1, IU - 1) = Q1/BRTOT
            if (KL > KB(UHS(JB))) then
                if (FRAC < 1.0 .and. FRAC /= 0.0) then
                    if (K == KB(IU) + 1) then
                        U(K - 1, IU - 1) = (Q1 + U(KL - 1, UHS(JB))*BHR1(KL - 1, UHS(JB))*(1.0 - FRAC))/BRTOT
                    else
                        U(K, IU - 1) = U(KL - 1, UHS(JB))*BHR1(KL - 1, UHS(JB))*(1.0 - FRAC)/BHR(K, IU - 1)
                    end if
                end if
                go to 100
            end if
            EL1 = ELR(K)
        end if
    end do
    100 continue

! Debug
!    QL=0.0; QR=0.0
!  do k=ktwb(jjw),kb(uhs(jb))
!    QL=QL+u(k,uhs(jb))*BHR1(k,uhs(jb))
!  enddo
!  do k=kt,kb(iu-1)
!    QR=QR+u(k,iu-1)*BHR1(k,iu-1)
!  enddo
!  if((QR-QL)/QL > 0.02)then
!     if(ncount.eq.0)open(1299,file='debug_out.txt',status='unknown')
!     ncount=ncount+1
!     write(1299,*)'***QL=',QL,' QR=',QR
!     write(1299,*)'UHS(JB)=',uhs(jb),'  iu-1=',iu-1,' ELWS=', elw
!     write(1299,*)'LEFT: K      U     BHR1    ELEV      H1     BHR    ELL'
!     do k=ktwb(jjw),kb(uhs(jb))
!       write(1299,'(i8,f8.3,f8.3,5f8.3)')k,u(k,uhs(jb)),BHR1(k,uhs(jb)),el(k,uhs(jb)),h1(k,uhs(jb)),bhr(k,uhs(jb)),ell(k)
!     enddo
!     write(1299,*)'RIGHT: K     U     BHR1    ELEV      H1     BHR    ELR'
!     do k=kt,kb(iu-1)
!       write(1299,'(i8,f8.3,f8.3,5f8.3)')k,u(k,iu-1),BHR1(k,iu-1),el(k,iu-1),h1(k,iu-1),bhr(k,iu-1),elr(k)
!     enddo
!  end if
! End debug

    return

!***********************************************************************************************************************************
!**                                              U P S T R E A M   W A T E R B O D Y                                              **
!***********************************************************************************************************************************

    entry UPSTREAM_WATERBODY()
    do JJB = 1, NBR
        if (UHS(JB) >= US(JJB) .and. UHS(JB) <= DS(JJB)) then
            exit
        end if
    end do
    do JJW = 1, NWB
        if (JJB >= BS(JJW) .and. JJB <= BE(JJW)) then
            exit
        end if
    end do
    do K = KTWB(JJW), KB(UHS(JB)) + 1
        ELL(K) = EL(K, UHS(JB)) - SINA(JJB)*DLX(UHS(JB))*0.5
    end do
    do K = KT, KB(IU) + 1
        ELR(K) = EL(K, IU) + SINA(JB)*DLX(IU)*0.5
    end do
    ELW = ELWS(UHS(JB)) !EL(KTWB(JJW),UHS(JB))-Z(UHS(JB))*COSA(JJB)
    EL1 = ELW - SINA(JJB)*DLX(UHS(JB))*0.5
! EL1  = ELW-(ELWS(UHS(JB)+1)-ELW)/(0.5*(DLX(UHS(JB)+DLX(JB)+1))*DLX(UHS(JB))*0.5   ! SW 7/17/09
    KL = KTWB(JJW) + 1
    do K = KT + 1, KB(IU) + 1
        if (ELR(K) >= ELL(KL)) then
            T1(K - 1, IU - 1) = T1(KL - 1, UHS(JB))
            T2(K - 1, IU - 1) = T2(KL - 1, UHS(JB))
            C1S(K - 1, IU - 1, CN(1:NAC)) = C1S(KL - 1, UHS(JB), CN(1:NAC))
            C1(K - 1, IU - 1, CN(1:NAC)) = C1S(KL - 1, UHS(JB), CN(1:NAC))
            C2(K - 1, IU - 1, CN(1:NAC)) = C1S(KL - 1, UHS(JB), CN(1:NAC))
            EL1 = ELR(K)
            if (ELR(K) == ELL(KL)) then
                KL = KL + 1
            end if
        else
            BRTOT = 0.0
            CL = 0.0
            T1L = 0.0
            T2L = 0.0
            do while (ELR(K) <= ELL(KL))
                if (KL == KTWB(JJW) + 1 .and. K == KT + 1) then
                    B1 = BH2(KTWB(JJW), UHS(JB))
                else
                    B1 = B(KL - 1, UHS(JB))*(EL1 - ELL(KL))
                end if
                BRTOT = BRTOT + B1
                T1L = T1L + B1*T1(KL - 1, UHS(JB))
                T2L = T2L + B1*T2(KL - 1, UHS(JB))
                CL(CN(1:NAC)) = CL(CN(1:NAC)) + B1*C1S(KL - 1, UHS(JB), CN(1:NAC))
                EL1 = ELL(KL)
                KL = KL + 1
                if (KL > KB(UHS(JB)) + 1) then
                    exit
                end if
            end do
            if (KL <= KB(UHS(JB)) + 1) then
                B1 = B(KL - 1, UHS(JB))*(EL1 - ELR(K))
                BRTOT = BRTOT + B1
                if (BRTOT > 0.0) then
                    T1(K - 1, IU - 1) = (T1L + B1*T1(KL - 1, UHS(JB)))/BRTOT
                    T2(K - 1, IU - 1) = (T2L + B1*T2(KL - 1, UHS(JB)))/BRTOT
                    C1S(K - 1, IU - 1, CN(1:NAC)) = (CL(CN(1:NAC)) + B1*C1S(KL - 1, UHS(JB), CN(1:NAC)))/BRTOT
                    C1(K - 1, IU - 1, CN(1:NAC)) = C1S(K - 1, IU - 1, CN(1:NAC))
                    C2(K - 1, IU - 1, CN(1:NAC)) = C1S(K - 1, IU - 1, CN(1:NAC))
                else
                    T1(K - 1, IU - 1) = T1(KL - 1, UHS(JB))
                    T2(K - 1, IU - 1) = T2(KL - 1, UHS(JB))
                    C1S(K - 1, IU - 1, CN(1:NAC)) = C1S(KL - 1, UHS(JB), CN(1:NAC))
                    C1(K - 1, IU - 1, CN(1:NAC)) = C1S(K - 1, IU - 1, CN(1:NAC))
                    C2(K - 1, IU - 1, CN(1:NAC)) = C1S(K - 1, IU - 1, CN(1:NAC))
                end if
            else
                if (BRTOT > 0.0) then
                    T1(K - 1, IU - 1) = T1L/BRTOT
                    T2(K - 1, IU - 1) = T2L/BRTOT
                    C1S(K - 1, IU - 1, CN(1:NAC)) = CL(CN(1:NAC))/BRTOT
                    C1(K - 1, IU - 1, CN(1:NAC)) = C1S(K - 1, IU - 1, CN(1:NAC))
                    C2(K - 1, IU - 1, CN(1:NAC)) = C1S(K - 1, IU - 1, CN(1:NAC))
                else
                    T1(K - 1, IU - 1) = T1(KL - 1, UHS(JB))
                    T2(K - 1, IU - 1) = T2(KL - 1, UHS(JB))
                    C1S(K - 1, IU - 1, CN(1:NAC)) = C1S(KL - 1, UHS(JB), CN(1:NAC))
                    C1(K - 1, IU - 1, CN(1:NAC)) = C1S(K - 1, IU - 1, CN(1:NAC))
                    C2(K - 1, IU - 1, CN(1:NAC)) = C1S(K - 1, IU - 1, CN(1:NAC))
                end if
                exit
            end if
            EL1 = ELR(K)
        end if
    end do
    return

!***********************************************************************************************************************************
!**                                            D O W N S T R E A M   W A T E R B O D Y                                            **
!***********************************************************************************************************************************

    entry DOWNSTREAM_WATERBODY()
    do JJB = 1, NBR
        if (CDHS(JB) >= CUS(JJB) .and. CDHS(JB) <= DS(JJB)) then
            exit
        end if
    end do
    do JJW = 1, NWB
        if (JJB >= BS(JJW) .and. JJB <= BE(JJW)) then
            exit
        end if
    end do
    KR = KTWB(JJW) + 1
    do K = KT, KB(ID) + 1
        ELL(K) = EL(K, ID) - SINA(JB)*DLX(ID)*0.5
    end do
    do K = KTWB(JJW), KB(CDHS(JB)) + 1
        ELR(K) = EL(K, CDHS(JB)) + SINA(JJB)*DLX(CDHS(JB))*0.5
    end do
    ELW = ELWS(ID) !EL(KTWB(JW),ID)-Z(ID)*COSA(JB)
    EL1 = ELW - SINA(JB)*DLX(ID)*0.5
!   EL1  = ELW+(ELW-ELWS(ID-1))/(0.5*(DLX(ID)+DLX(ID-1))*DLX(ID)*0.5   ! SW 7/17/09
    IDT = ID + 1
    do K = KT + 1, KB(ID) + 1
        if (ELL(K) >= ELR(KR)) then
            T1(K - 1, IDT) = T1(KR - 1, CDHS(JB))
            T2(K - 1, IDT) = T2(KR - 1, CDHS(JB))
            C1S(K - 1, IDT, CN(1:NAC)) = C1S(KR - 1, CDHS(JB), CN(1:NAC))
            C1(K - 1, IDT, CN(1:NAC)) = C1S(KR - 1, CDHS(JB), CN(1:NAC))
            C2(K - 1, IDT, CN(1:NAC)) = C1S(KR - 1, CDHS(JB), CN(1:NAC))
            EL1 = ELL(K)
            if (ELL(K) == ELR(KR)) then
                KR = KR + 1
            end if
            if (KR > KB(CDHS(JB)) + 1) then
                exit
            end if
        else
            BRTOT = 0.0
            CL = 0.0
            T1L = 0.0
            T2L = 0.0
            do while (ELL(K) <= ELR(KR))
                if (KR == KTWB(JJW) + 1 .and. K == KT + 1) then
                    B1 = BH2(KTWB(JJW), CDHS(JB))
                    BRTOT = BRTOT + B1
                else
                    B1 = B(KR - 1, CDHS(JB))*(EL1 - ELR(KR))
                    BRTOT = BRTOT + B1
                end if
                T1L = T1L + B1*T1(KR - 1, CDHS(JB))
                T2L = T2L + B1*T2(KR - 1, CDHS(JB))
                CL(CN(1:NAC)) = CL(CN(1:NAC)) + B1*C1S(KR - 1, CDHS(JB), CN(1:NAC))
                EL1 = ELR(KR)
                KR = KR + 1
                if (KR > KB(CDHS(JB)) + 1) then
                    exit
                end if
            end do
            if (KR <= KB(CDHS(JB) + 1)) then
                B1 = B(KR - 1, CDHS(JB))*(EL1 - ELL(K))
            else
                B1 = 0.0
            end if
            BRTOT = BRTOT + B1
            if (BRTOT == 0.0) then
                exit
            end if
            T1(K - 1, IDT) = (T1L + B1*T1(KR - 1, CDHS(JB)))/BRTOT
            T2(K - 1, IDT) = (T2L + B1*T2(KR - 1, CDHS(JB)))/BRTOT
            C1S(K - 1, IDT, CN(1:NAC)) = (CL(CN(1:NAC)) + B1*C1S(KR - 1, CDHS(JB), CN(1:NAC)))/BRTOT
            C1(K - 1, IDT, CN(1:NAC)) = C1S(K - 1, IDT, CN(1:NAC))
            C2(K - 1, IDT, CN(1:NAC)) = C1S(K - 1, IDT, CN(1:NAC))
            EL1 = ELL(K)
        end if
    end do
    return

!***********************************************************************************************************************************
!**                                                 U P S T R E A M   B R A N C H                                                 **
!***********************************************************************************************************************************

    entry UPSTREAM_BRANCH()
    do JJW = 1, NWB
        if (JJB >= BS(JJW) .and. JJB <= BE(JJW)) then
            exit
        end if
    end do
    do K = KT, KB(I) + 1
        ELL(K) = EL(K, I)
    end do
    do K = KTWB(JWUH(JB)), KB(CUS(JJB)) + 1
        ELR(K) = EL(K, CUS(JJB)) + SINA(JJB)*DLX(CUS(JJB))*0.5
    end do
    ELW = ELWS(CUS(JJB)) !EL(KTWB(JWUH(JB)),CUS(JJB))-Z(CUS(JJB))*COSA(JJB)
    EL1 = ELW + SINA(JJB)*DLX(CUS(JJB))*0.5
!   EL1  = ELW-(ELWS(CUS(JJB)+1)-ELW)/(0.5*(DLX(CUS(JJB)+DLX(CUS(JJB)+1))*DLX(CUS(JJB))*0.5   ! SW 7/17/09
    KR = KTWB(JWUH(JB)) + 1
    do K = KT + 1, KB(I) + 1
        if (ELL(K) >= ELR(KR)) then
            Q1 = VOLUH2(KR - 1, JJB)/DLT
            B1 = BHR(KR - 1, CUS(JJB) - 1)
            if (KR == KTWB(JWUH(JB)) + 1) then
                B1 = BHR2(KT, CUS(JJB) - 1)
            end if
            U1 = U(KR - 1, CUS(JJB) - 1)*B1
            HT = H2(KR - 1, JWUH(JB))
            if (KR == KTWB(JWUH(JB)) + 1) then
                HT = H2(KT, CUS(JJB))
            end if
            Q1 = Q1*(EL1 - ELL(K))/HT
            B2 = BHR(K - 1, I)
            if (K == KTWB(JW) + 1) then
                B2 = BHR2(KT, I)
            end if
            UXBR(K - 1, I) = UXBR(K - 1, I) + ABS(U1/B2*COS(BETABR)*Q1)/DLX(I)
            UYBR(K - 1, I) = UYBR(K - 1, I) + ABS(Q1*SIN(BETABR))
            EL1 = ELL(K)
            if (ELL(K) == ELR(KR)) then
                KR = KR + 1
            end if
        else
            U1 = 0.0
            Q1 = 0.0
            BRTOT = 0.0
            do while (ELL(K) <= ELR(KR))
                if (KR /= KTWB(JWUH(JB)) + 1) then
                    FRAC = (EL1 - ELR(KR))/H(KR - 1, JWUH(JB))
                    B1 = BHR(KR - 1, CUS(JJB) - 1)
                else
                    FRAC = (EL1 - ELR(KR))/H2(KT, CUS(JJB))
                    B1 = BHR2(KT, CUS(JJB) - 1)
                end if
                U1 = U1 + U(KR - 1, CUS(JJB) - 1)*B1*FRAC
                Q1 = Q1 + VOLUH2(KR - 1, JJB)/DLT*FRAC
                EL1 = ELR(KR)
                KR = KR + 1
                if (KR > KB(CUS(JJB) + 1)) then
                    exit
                end if
            end do
            if (K == KTWB(JW) + 1) then
                B2 = BHR2(KT, I)
            else
                B2 = BHR(K - 1, I)
            end if
            if (H(KR - 1, JWUH(JB)) /= 0.0) then
                if (KR - 1 == KTWB(JWUH(JB))) then
                    HT = H2(KT, CUS(JJB))
                else
                    HT = H2(KR - 1, JWUH(JB))
                end if
                FRAC = (EL1 - ELL(KR - 1))/HT
                Q1 = Q1 + FRAC*VOLUH2(KR - 1, JJB)/DLT
                UXBR(K - 1, I) = UXBR(K - 1, I) + ABS(U1/B2*COS(BETABR)*Q1)/DLX(I)
                UYBR(K - 1, I) = UYBR(K - 1, I) + ABS(Q1*SIN(BETABR))
            end if
            if (KR > KB(CUS(JJB) + 1)) then
                exit
            end if
            EL1 = ELL(K)
        end if
    end do
    return

!***********************************************************************************************************************************
!**                                               D O W N S T R E A M   B R A N C H                                               **
!***********************************************************************************************************************************

    entry DOWNSTREAM_BRANCH()
    do JJW = 1, NWB
        if (JJB >= BS(JJW) .and. JJB <= BE(JJW)) then
            exit
        end if
    end do
    do K = KT, KB(I) + 1
        ELR(K) = EL(K, I)
    end do
    do K = KTWB(JJW), KB(DS(JJB)) + 1
        ELL(K) = EL(K, DS(JJB)) + SINA(JJB)*DLX(DS(JJB))*0.5
    end do
    ELW = ELWS(DS(JJB)) !EL(KTWB(JJW),DS(JJB))-Z(DS(JJB))*COSA(JJB)
    EL1 = ELW - SINA(JJB)*DLX(DS(JJB))*0.5
!  IF(SLOPE(JJB) /= 0.0)THEN
!  EL1  = ELW+(ELW-ELWS(DS(JJB)-1))/(0.5*(DLX(DS(JJB)+DLX(DS(JJB)+1))*DLX(DS(JJB))*0.5   ! SW 7/17/09
!  ELSE
!  EL1=ELW
!  ENDIF
    KL = KTWB(JJW) + 1
    do K = KT + 1, KB(I) + 1
        if (ELR(K) >= ELL(KL)) then
            Q1 = VOLDH2(KL - 1, JJB)/DLT
            B1 = BHR(KL - 1, DS(JJB) - 1)
            if (KL == KTWB(JJW) + 1) then
                B1 = BHR2(KT, DS(JJB))
            end if
            U1 = U(KL - 1, DS(JJB) - 1)*B1
            HT = H2(KL - 1, JJW)
            if (KL == KTWB(JJW) + 1) then
                HT = H2(KT, DS(JJB))
            end if
            Q1 = Q1*(EL1 - ELR(K))/HT
            B2 = BHR(K - 1, I)
            if (K == KTWB(JW) + 1) then
                B2 = BHR2(KT, I)
            end if
            UXBR(K - 1, I) = UXBR(K - 1, I) + ABS(U1/B2*COS(BETABR)*Q1)/DLX(I)
            UYBR(K - 1, I) = UYBR(K - 1, I) + ABS(Q1*SIN(BETABR))
            EL1 = ELR(K)
            if (ELR(K) == ELL(KL)) then
                KL = KL + 1
            end if
        else
            U1 = 0.0
            Q1 = 0.0
            BRTOT = 0.0
            do while (ELR(K) <= ELL(KL))
                if (KL /= KTWB(JJW) + 1) then
                    FRAC = (EL1 - ELL(KL))/H(KL - 1, JJW)
                    B1 = BHR(KL - 1, DS(JJB) - 1)
                else
                    FRAC = (EL1 - ELL(KL))/H2(KT, DS(JJB))
                    B1 = BHR2(KT, DS(JJB))
                end if
                U1 = U1 + U(KL - 1, DS(JJB))*B1*FRAC
                Q1 = Q1 + VOLDH2(KL - 1, JJB)/DLT*FRAC
                EL1 = ELL(KL)
                KL = KL + 1
                if (KL > KB(DS(JJB) + 1)) then
                    exit
                end if
            end do
            if (K == KTWB(JW) + 1) then
                B2 = BHR2(KT, I)
            else
                B2 = BHR(K - 1, I)
            end if
            if (H(KL - 1, JJW) /= 0.0) then
                if (KL - 1 == KTWB(JJW)) then
                    HT = H2(KT, DS(JJB))
                else
                    HT = H2(KL - 1, JJW)
                end if
                FRAC = (EL1 - ELL(KL - 1))/HT
                Q1 = Q1 + FRAC*VOLDH2(KL - 1, JJB)/DLT
                UXBR(K - 1, I) = UXBR(K - 1, I) + ABS(U1/B2*COS(BETABR)*Q1)/DLX(I)
                UYBR(K - 1, I) = UYBR(K - 1, I) + ABS(Q1*SIN(BETABR))
            end if
            if (KL > KB(DS(JJB) + 1)) then
                exit
            end if
            EL1 = ELR(K)
        end if
    end do
    return

!***********************************************************************************************************************************
!**                                                   U P S T R E A M   F L O W                                                   **
!***********************************************************************************************************************************

    entry UPSTREAM_FLOW()
    do JJB = 1, NBR
        if (UHS(JB) >= US(JJB) .and. UHS(JB) <= DS(JJB)) then
            exit
        end if
    end do
    do JJW = 1, NWB
        if (JJB >= BS(JJW) .and. JJB <= BE(JJW)) then
            exit
        end if
    end do
    do K = KTWB(JWUH(JB)), KB(UHS(JB)) + 1
        ELL(K) = EL(K, UHS(JB))
    end do
    do K = KT, KB(IU) + 1
        ELR(K) = EL(K, IU) + SINA(JB)*DLX(IU)*0.5
    end do
    ELW = ELWS(IU) !EL(KTWB(JW),IU)-Z(IU)*COSA(JB)
    EL1 = ELW + SINA(JB)*DLX(IU)*0.5
!  IF(SLOPE(JB) /= 0.0)THEN
!  EL1  = ELW-(ELWS(IU+1)-ELW)/(0.5*(DLX(IU)+DLX(IU+1))*DLX(IU)*0.5   ! SW 7/17/09
!  ELSE
!  EL1=ELW
!  ENDIF
    KR = KT + 1
    do K = KTWB(JWUH(JB)) + 1, KB(UHS(JB)) + 1
        if (ELL(K) >= ELR(KR)) then
            Q1 = VOLUH2(KR - 1, JB)/DLT
            HT = H2(KR - 1, JW)
            if (KR == KTWB(JW) + 1) then
                QU(K - 1, UHS(JB)) = Q1
                QSS(K - 1, UHS(JB)) = QSS(K - 1, UHS(JB)) - Q1
            else
                QU(K - 1, UHS(JB)) = Q1*(EL1 - ELL(K))/HT
                QSS(K - 1, UHS(JB)) = QSS(K - 1, UHS(JB)) - QU(K - 1, UHS(JB))
            end if
            EL1 = ELL(K)
            if (ELL(K) == ELR(KR)) then
                KR = KR + 1
            end if
        else
            Q1 = 0.0
            do while (ELL(K) <= ELR(KR))
                if (KR /= KTWB(JW) + 1) then
                    FRAC = (EL1 - ELR(KR))/H(KR - 1, JW)
                else
                    FRAC = (EL1 - ELR(KR))/H2(KT, IU)
                end if
                Q1 = Q1 + VOLUH2(KR - 1, JB)/DLT*FRAC
                EL1 = ELR(KR)
                KR = KR + 1
                if (KR > KB(IU) + 1) then
                    exit
                end if
            end do
            if (H(KR - 1, JW) /= 0.0) then
                FRAC = (EL1 - ELL(K))/H(KR - 1, JW)
                QU(K - 1, UHS(JB)) = Q1 + FRAC*VOLUH2(KR - 1, JB)/DLT
            else
                QU(K - 1, UHS(JB)) = Q1
            end if
            QSS(K - 1, UHS(JB)) = QSS(K - 1, UHS(JB)) - QU(K - 1, UHS(JB))
            if (KR > KB(IU) + 1) then
                exit
            end if
            EL1 = ELL(K)
        end if
    end do
    return

!***********************************************************************************************************************************
!**                                                 D O W N S T R E A M   F L O W                                                 **
!***********************************************************************************************************************************

    entry DOWNSTREAM_FLOW()
    do JJB = 1, NBR
        if (CDHS(JB) >= US(JJB) .and. CDHS(JB) <= DS(JJB)) then
            exit
        end if
    end do
    do JJW = 1, NWB
        if (JJB >= BS(JJW) .and. JJB <= BE(JJW)) then
            exit
        end if
    end do
    do K = KTWB(JJW), KB(CDHS(JB)) + 1
        ELR(K) = EL(K, CDHS(JB))
    end do
    do K = KT, KB(ID) + 1
        ELL(K) = EL(K, ID) - SINA(JB)*DLX(ID)*0.5
    end do
    ELW = ELWS(ID) !EL(KTWB(JW),ID)-Z(ID)*COSA(JB)
    EL1 = ELW - SINA(JB)*DLX(ID)*0.5
!  IF(SLOPE(JB) /= 0.0)THEN
!  EL1  = ELW+(ELW-ELWS(ID-1))/(0.5*(DLX(ID)+DLX(ID-1))*DLX(ID)*0.5   ! SW 7/17/09
!  ELSE
!  EL1=ELW
!  ENDIF
    KL = KT + 1
    do K = KTWB(JJW) + 1, KB(CDHS(JB)) + 1
        if (ELR(K) >= ELL(KL)) then
            Q1 = VOLDH2(KL - 1, JB)/DLT
            HT = H2(KL - 1, JW)
            if (KL == KTWB(JW) + 1) then
                HT = H2(KT, ID)
            end if
            QD(K - 1, CDHS(JB)) = Q1*(EL1 - ELR(K))/HT
            QSS(K - 1, CDHS(JB)) = QSS(K - 1, CDHS(JB)) + QD(K - 1, DHS(JB))
            EL1 = ELR(K)
            if (ELR(K) == ELL(KL)) then
                KL = KL + 1
            end if
        else
            Q1 = 0.0
            do while (ELR(K) <= ELL(KL))
                if (KL /= KTWB(JW) + 1) then
                    FRAC = (EL1 - ELL(KL))/H(KL - 1, JW)
                else
                    FRAC = (EL1 - ELL(KL))/H2(KT, ID)
                end if
                Q1 = Q1 + VOLDH2(KL - 1, JB)/DLT*FRAC
                EL1 = ELL(KL)
                KL = KL + 1
                if (KL > KB(ID) + 1) then
                    exit
                end if
            end do
            if (H(KL - 1, JW) /= 0.0) then
                FRAC = (EL1 - ELR(K))/H(KL - 1, JW)
                QD(K - 1, CDHS(JB)) = Q1 + FRAC*VOLDH2(KL - 1, JB)/DLT
            else
                QD(K - 1, CDHS(JB)) = Q1
            end if
            QSS(K - 1, CDHS(JB)) = QSS(K - 1, CDHS(JB)) + QD(K - 1, CDHS(JB))
            if (KL > KB(ID) + 1) then
                exit
            end if
            EL1 = ELR(K)
        end if
    end do
    return

!***********************************************************************************************************************************
!**                                            U P S T R E A M   C O N S T I T U E N T                                            **
!***********************************************************************************************************************************

    entry UPSTREAM_CONSTITUENT(C, SS)
    do JJB = 1, NBR
        if (UHS(JB) >= US(JJB) .and. UHS(JB) <= DS(JJB)) then
            exit
        end if
    end do
    do JJW = 1, NWB
        if (JJB >= BS(JJW) .and. JJB <= BE(JJW)) then
            exit
        end if
    end do
    do K = KTWB(JWUH(JB)), KB(UHS(JB)) + 1
        ELL(K) = EL(K, UHS(JB))
    end do
    do K = KT, KB(IU) + 1
        ELR(K) = EL(K, IU) + SINA(JB)*DLX(IU)*0.5
    end do
    ELW = ELWS(IU) !EL(KTWB(JW),IU)-Z(IU)*COSA(JB)
    EL1 = ELW + SINA(JB)*DLX(IU)*0.5
!  IF(SLOPE(JB) /= 0.0)THEN
!  EL1  = ELW-(ELWS(IU+1)-ELW)/(0.5*(DLX(IU)+DLX(IU+1))*DLX(IU)*0.5   ! SW 7/17/09
!  ELSE
!  EL1=ELW
!  ENDIF
    KR = KT + 1
    do K = KTWB(JWUH(JB)) + 1, KB(UHS(JB)) + 1
        IUT = IU
        if (QU(K - 1, UHS(JB)) >= 0.0) then
            IUT = IU - 1
        end if
        if (ELL(K) >= ELR(KR)) then
            T1L = C(KR - 1, IUT)
            SS(K - 1, UHS(JB)) = SS(K - 1, UHS(JB)) - T1L*QU(K - 1, UHS(JB))
            EL1 = ELL(K)
            if (ELL(K) == ELR(KR)) then
                KR = KR + 1
            end if
        else
            T1L = 0.0
            BRTOT = 0.0
            do while (ELL(K) <= ELR(KR))
                if (KR == KT + 1 .and. K == KTWB(JWUH(JB)) + 1) then
                    B1 = BH2(KT, IU)
                    BRTOT = BRTOT + B1
                else
                    B1 = B(KR - 1, IU)*(EL1 - ELR(KR))
                    BRTOT = BRTOT + B1
                end if
                IUT = IU
                if (QU(K - 1, UHS(JB)) >= 0.0) then
                    IUT = IU - 1
                end if
                T1L = T1L + B1*C(KR - 1, IUT)
                EL1 = ELR(KR)
                KR = KR + 1
                if (KR > KB(IU) + 1) then
                    exit
                end if
            end do
            IUT = IU
            if (QU(K - 1, UHS(JB)) >= 0.0) then
                IUT = IU - 1
            end if
            B1 = B(KR - 1, IU)*(EL1 - ELL(K))
            BRTOT = BRTOT + B1
            T1L = (T1L + B1*C(KR - 1, IUT))/BRTOT
            SS(K - 1, UHS(JB)) = TSS(K - 1, UHS(JB)) - T1L*QU(K - 1, UHS(JB))
            if (KR > KB(IU) + 1) then
                exit
            end if
            EL1 = ELL(K)
        end if
    end do
    return

!***********************************************************************************************************************************
!**                                          D O W N S T R E A M   C O N S T I T U E N T                                          **
!***********************************************************************************************************************************

    entry DOWNSTREAM_CONSTITUENT(C, SS)
    do JJB = 1, NBR
        if (CDHS(JB) >= US(JJB) .and. CDHS(JB) <= DS(JJB)) then
            exit
        end if
    end do
    do JJW = 1, NWB
        if (JJB >= BS(JJW) .and. JJB <= BE(JJW)) then
            exit
        end if
    end do
    do K = KTWB(JJW), KB(CDHS(JB)) + 1
        ELR(K) = EL(K, CDHS(JB))
    end do
    do K = KT, KB(ID) + 1
        ELL(K) = EL(K, ID) + SINA(JB)*DLX(ID)*0.5
    end do
    ELW = ELWS(ID) !EL(KTWB(JW),ID)-Z(ID)*COSA(JB)
    EL1 = ELW + SINA(JB)*DLX(ID)*0.5
!  IF(SLOPE(JB) /= 0.0)THEN
!  EL1  = ELW+(ELW-ELWS(ID-1))/(0.5*(DLX(ID)+DLX(ID-1))*DLX(ID)*0.5   ! SW 7/17/09
!  ELSE
!  EL1=ELW
!  ENDIF
    KL = KT + 1
    do K = KTWB(JJW) + 1, KB(CDHS(JB)) + 1
        IDT = ID + 1
        if (QD(K - 1, CDHS(JB)) >= 0.0) then
            IDT = ID
        end if
        if (ELR(K) >= ELL(KL)) then
            T1L = C(KL - 1, IDT)
            SS(K - 1, CDHS(JB)) = SS(K - 1, CDHS(JB)) + T1L*QD(K - 1, CDHS(JB))
            EL1 = ELR(K)
            if (ELR(K) == ELL(KL)) then
                KL = KL + 1
            end if
        else
            T1L = 0.0
            BRTOT = 0.0
            do while (ELR(K) <= ELL(KL))
                if (KL == KTWB(JW) + 1 .and. K == KTWB(JJW) + 1) then
                    B1 = BH2(KT, ID)
                    BRTOT = BRTOT + B1
                else
                    B1 = B(KL - 1, ID)*(EL1 - ELL(KL))
                    BRTOT = BRTOT + B1
                end if
                T1L = T1L + B1*C(KL - 1, IDT)
                EL1 = ELL(KL)
                KL = KL + 1
                if (KL > KB(ID)) then
                    B1 = B(KL - 1, ID)*(EL1 - ELR(K))
                    BRTOT = BRTOT + B1
                    T1L = (T1L + B1*C(KL - 1, IDT))/BRTOT
                    SS(K - 1, CDHS(JB)) = SS(K - 1, CDHS(JB)) + T1L*QD(K - 1, CDHS(JB))
                    go to 200
                end if
            end do
            B1 = B(KL - 1, ID)*(EL1 - ELR(K))
            BRTOT = BRTOT + B1
            T1L = (T1L + B1*C(KL - 1, IDT))/BRTOT
            SS(K - 1, CDHS(JB)) = SS(K - 1, CDHS(JB)) + T1L*QD(K - 1, CDHS(JB))
            EL1 = ELL(KL)
        end if
    end do
    200 continue
    return
    entry DEALLOCATE_WATERBODY()
    deallocate(ELL, ELR, CL, QU, QD)
    return
end subroutine WATERBODY
