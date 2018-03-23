program dp1
implicit real*8 (a-h,o-z)
external force
common/mu/amu/tru/lgcl/kk/Kcounter(5)
dimension x(14)
logical start,lgcl,test
character(1) :: choice
    !считаем из файла начальные услови§
    open(100,file = "1.txt")
    read(100,*) xnach, xkonza, delx, cinit, cfin, delc, nv, ss, ni, t0, tp0, tf, Gm1, Gm2, P1P2, RCh
    close(100)
    
    !высчитаем ню 
    amu = Gm2/(Gm1 + Gm2)
    !зададим значение угла наклона
    AnglI = 0.d0
    RChinConvUnits = RCh/P1P2
    C_cr = 0.3620423879030583d0
    
    !нек-е координаты точек равновеси§(Оагранджа)
    x_L1 = 0.5929868613955974d0
    x_L2 = 1.262525159065689d0
    x_L3 = -1.045151336842185d0
    y_L4 = 0.8660254037844386d0

    !посчитаем машинный эпсилон
    eps = 1.d0
    do while (eps/2.d0 + 1.d0 > 1.d0)
        eps = eps/2.d0
    enddo
    write(*,*) 'computer epsilon:', eps 
    
    !блок выбора типа решени§
    write(*,*) 'select count with or without variation Y/N'
    read(*,*) choice
    if (choice == 'Y' .or. choice == 'y') then 
        lgcl = .true.
        nv = 14
    endif
    
    Kcounter = 0
    open(111, file = 'The orbit is ejected into the Pluto sphere through L1.dat')
    open(222, file = 'the orbit is ejected into the outer sphere by L2.dat')
    open(333, file = 'The orbit is ejected into the outer sphere by L1.dat')
    open(444, file = 'clash with Charon.dat')
    open(555, file = 'Orbit complete integration.dat')
   
    !основной цикл по посто§нной §коби
    do while (cinit <= cfin)
        write(*,*) 'c = ',cinit,' : '
        xnach1 = xnach
        !основной цикл по оси х
        do while (xnach1 <= xkonza)
            !проверим на допустимость движени§
            if (2.d0*Sigm(xnach1,0.d0,0.d0) - cinit >= 0.d0) then
                !инициализаци§ начального вектора х
                x = 0.d0
                call xinit1(xnach1,cinit,AnglI,x)
                tm = t0
                tp = tp0
                start = .true.
                !вложенный цикл по времени
                test = .true.
                do while(test)
                    if (tm + tp > tf) then
                        tp = tf - tm
                        call rada15(tm,x,nv,tp,ss,ni,ns,nf,force,start)
                    else
                        call rada15(tm,x,nv,tp,ss,ni,ns,nf,force,start)
                    endif
                    !проверка всех четырех условий
                    call VerificationOfConditions(x,tm,tf,eps,RChinConvUnits,x_L1,x_L2,x_L3,y_L4,test,KeyValue)
                enddo
                
                !расчет значени§ функции устойчивости по Тиллу, индикатора хаоса
                call ValuesOfIndicators(x,tm,C_cr,Omegno,FuncHill)
                !вывод данных в файлы
                call PrintInFile(x,KeyValue,tm,Omegno,FuncHill,cinit,xnach1)
            endif
            xnach1 = xnach1 + delx
        enddo
        write(*,*) Kcounter(1),';',Kcounter(2),';',Kcounter(3),';',Kcounter(4),';',Kcounter(5)
        cinit = cinit + delc
    enddo
    
    close(111)
    close(222)
    close(333)
    close(444)
    close(555)
    end program dp1
    
    subroutine force(tm,x,f)
    implicit real*8(a-h,o-z)
    dimension x(14),f(14)
    common/mu/amu/tru/lgcl
        !подсчет правых частей основных уравнений
        r1 = r1_res(x(1),x(2),x(3))
        r2 = r2_res(x(1),x(2),x(3))
        f(1) = x(4)
        f(2) = x(5)
        f(3) = x(6)
        f(4) = 2.d0*x(5) + x(1) - (1.d0 - amu)*(x(1) + amu)/r1**3 - amu*((x(1) - 1.d0 + amu)/r2**3)
        f(5) = -2.d0*x(4) + x(2) - (1.d0 - amu)*x(2)/r1**3 - amu*(x(2)/r2**3)
        f(6) = -(1.d0 - amu)*x(3)/r1**3 - amu*x(3)/r2**3
        
        !условие на подсчет вариации
        if (lgcl) then
            !вычисление нек-х общих членов выражений вариации
            call CoefOfVar_tion(x(1),r1,r2,A,B,D,E,ff,G)
            !подсчет правых частей ур-й вариаций
            f(7) = x(10)
            f(8) = x(11)
            f(9) = x(12)
            f(10) = ff*x(7) + G*(x(2)*x(8) + x(3)*x(9)) + 2.d0*x(11)
            f(11) = (D + E*x(2)**2)*x(8) + E*x(2)*x(3)*x(9) + G*x(2)*x(7) - 2.d0*x(10)
            f(12) = (D - 1.d0 + E*x(3)**2)*x(9) + E*x(2)*x(3)*x(8) + G*x(3)*x(7)
            
            !подсчет компанент дл§ индикатора
            SqVar_ionRate = 0.d0
            Prod_D_F = 0.d0
            SqFRate = 0.d0
            do i = 1,6
                SqVar_ionRate = SqVar_ionRate +  x(i + 6)**2 
                Prod_D_F = Prod_D_F + x(i + 6)*f(i)
                SqFRate = SqFRate + f(i)**2
            enddo
            
            !диф ур-§ дл§ индикатора хаоса
            f(13) = dlog(dsqrt(SqVar_ionRate - Prod_D_F**2/SqFRate))
            if (tm == 0.d0) then
                f(14) = 0.d0
            else
                f(14) = f(13)/tm 
            endif
        endif
    return
    end
    
    subroutine VerificationOfConditions(x,tm,tf,eps,RCh,x_L1,x_L2,x_L3,y_L4,test,KeyValue)
    implicit real*8 (a-h,o-z)
    common/mu/amu/kk/Kcounter(5)
    dimension x(14)
    logical test
        !проверим остаетс§ ли орбита в сфере второго тела(Тарона)
        if (dabs(tm - tf) >= eps) then
            !подсчет радиусов
            r1 = R1_res(x(1), x(2), x(3))
            r2 = R2_res(x(1), x(2), x(3))
            !зададим нек-е пределы точности
            delt1 = 0.12d0
            delt2 = 0.10d0
            !посчитаем значение координаты по х дл§ Њлутона
            x_P1 = -amu
            !орбита выбрасываетс§ в сферу Њлутона через L1
            if ((x(1) < x_L1 - delt1).and.(r1 <= dabs(x_L3 - x_P1))) then
                KeyValue = 1
                Kcounter(1) = Kcounter(1) + 1
                test = .false.
            !орбита выбрасываетс§ во внешнию сферу через L2
            else if ((x(1) > x_L2 + delt2).or.(dsqrt(x(2)**2 + x(3)**2) > y_L4)) then
                KeyValue = 2
                Kcounter(2) = Kcounter(2) + 1
                test = .false.
            !орбита выбрасываетс§ во внешнию сферу через L1
            else if ((x(1) < x_L1 - delt1).and.(r1 > dabs(x_L3 - x_P1))) then
                KeyValue = 3
                Kcounter(3) = Kcounter(3) + 1
                test = .false.
            !столкновение с Тароном
            else if (r2 <= RCh) then
                KeyValue = 4
                Kcounter(4) = Kcounter(4) + 1
                test = .false.
            endif
        else
            KeyValue = 5
            Kcounter(5) = Kcounter(5) + 1
            test = .false.
        endif
    return  
    end
    
    subroutine ValuesOfIndicators(x,tm,C_cr,Omegno,FuncHill)
    implicit real*8 (a-h,o-z)
    common/mu/amu
    dimension x(14)
        !значение индикатора хаоса
        Omegno = 2.d0*(x(13) - x(14))/tm
        !значени§ функции устойчивости по Тиллу
            !обобщенный потенциал равен:
            GenPotent = Sigm(x(1),x(2),x(3)) + amu*(1.d0 - amu)/2.d0
            !обобщенное значение интеграа §коби
            Gen_C = 2.d0*GenPotent - x(4)**2 - x(5)**2 - x(6)**2
            !значение функции устойчивости по Тиллу
            FuncHill = (Gen_C - C_cr)/C_cr
    return
    end
    
    subroutine PrintInFile(x,KeyValue,tm,Omegno,FuncHill,cinit,xnach1)
    implicit real*8 (a-h,o-z)
    CHARACTER(LEN=63) :: FMT = '(F24.16,3X,F24.16,3X,F24.16,3X,F24.16,3X,F24.16,3X,F24.16,3X)'
    !CHARACTER(LEN=30) :: FMT = "7(F20.16,3X)"
    dimension x(14)
        CJ = 2.d0*Sigm(x(1),x(2),x(3)) - x(4)**2 - x(5)**2 - x(6)**2
        deltC = (CJ - cinit)/CJ   
        
        select case (KeyValue)
        case (1)
            write(111,FMT)tm,deltc,cinit,xnach1,omegno,funchill
        case (2)
            write(222,FMT)tm,deltc,cinit,xnach1,omegno,funchill
        case (3)
            write(333,FMT)tm,deltc,cinit,xnach1,omegno,funchill
        case (4)
            write(444,FMT)tm,deltc,cinit,xnach1,omegno,funchill
        case (5)
            write(555,FMT)tm,deltc,cinit,xnach1,omegno,funchill
        end select
    return
    end
    
    function R1_res(x, y, z)
    implicit real*8 (a-h,o-z)    
    common/mu/amu
        R1_res = dsqrt((x + amu)**2 + y**2 + z**2)
    return
    end
    
    function R2_res (x, y, z)
    implicit real*8 (a-h,o-z)        
    common/mu/amu
        R2_res = dsqrt((x - 1.d0 + amu)**2 + y**2 + z**2)
    return
    end
    
    subroutine CalcVar_tion(x,r1,r2)
    implicit real*8 (a-h,o-z)
    dimension x(14)
    common/mu/amu
        
        DifC_DifX = 2.d0*(x(1) - (1.d0 - amu)*(x(1) + amu)/r1**3 - amu*(x(1) - 1.d0 + amu)/r2**3)
        DifC_DifY = 2.d0*(x(2) - (1.d0 - amu)*x(2)/r1**3 - amu*x(2)/r2**3)
        DifC_DifZ = - 2.d0*((1.d0 - amu)*x(3)/r1*3 + amu*x(3)/r2**3)
        DifC_DifXp = - 2.d0*x(4)
        DifC_DifYp = - 2.d0*x(5)
        DifC_DifZp = - 2.d0*x(6)
        
        GradNorm = dsqrt(DifC_DifX**2 + DifC_DifY**2 + DifC_DifZ**2 + DifC_DifXp**2 + DifC_DifYp**2 + DifC_DifZp**2)
    
        x(7) = DifC_DifX/GradNorm
        x(8) = DifC_DifY/GradNorm
        x(9) = DifC_DifZ/GradNorm
        x(10) = DifC_DifXp/GradNorm
        x(11) = DifC_DifYp/GradNorm
        x(12) = DifC_DifZp/GradNorm
    return
    end
    
    subroutine CoefOfVar_tion(x,r1,r2,A,B,D,E,F,G)
    implicit real*8 (a-h,o-z)
    common/mu/amu
        A = 3.d0*(1.d0 - amu)/r1**5
        B = 3.d0*amu/r2**5
        D= 1.d0 - (1.d0 - amu)/(r1*r1*r1) - amu/(r2*r2*r2)
        D = 1.d0 - (1.d0 - amu)/r1**3 - amu/r2**3
        E = A + B
        F = D + A*(x + amu)**2 + B*(x - 1.d0 + amu)**2
        G = A*(x + amu) + B*(x - 1.d0 + amu)
    return
    end
    
    subroutine xinit1(xnach,c0,AnglI,xinit)
    implicit real*8 (a-h,o-z)
    common/mu/amu/tru/lgcl
    dimension xinit(14)
        xinit(1) = xnach
        xinit(2) = 0.d0
        xinit(3) = 0.d0
        xinit(4) = 0.d0
        Vstr = dsqrt(2.d0*Sigm(xinit(1),xinit(2),xinit(3)) - c0)
        xinit(5) = Vstr*dcos(AnglI)
        xinit(6) = Vstr*dsin(AnglI)
        if (lgcl) then
            r1 = r1_res(xinit(1),xinit(2),xinit(3))
            r2 = r2_res(xinit(1),xinit(2),xinit(3))
            call CalcVar_tion(xinit,r1,r2)
        endif
    return
    end
    
    function Sigm (x,y,z)
    implicit real*8 (a-h,o-z)
    common/mu/amu
        r1 = r1_res(x,y,z)
        r2 = r2_res(x,y,z)
        Sigm = x**2/2.d0 + y**2/2.d0 + (1.d0 - amu)/r1 + amu/r2       
    return
    end
    
    SUBROUTINE RADA15(TM,X,NV,TP,SS,NI,NS,NF,FORCE,START)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NV),BD(7),H(8),W(8),C(21),D(21),R(21)
      COMMON/F1/F1(14)/FJ/FJ(14)/Y/Y(14)/ST/T/pw/pw!по размеру x(14)
      COMMON/B/B(7,14)/G/G(7,14)/E/E(7,14)!по размеру x(14)
      LOGICAL START,NSF
      DATA ZERO,ONE,SR/0.D0,1.D0,1.4D0/
      DATA H/0.D0,.05626256053692215D0,.18024069173689236D0,&
      .35262471711316964D0, .54715362633055538D0, .73421017721541053D0,&
      .88532094683909577D0, .97752061356128750D0/
      DATA W/.5D0,.3333333333333333D0,.25D0,.2D0,.1666666666666667D0,&
            .1428571428571429D0,.125D0,.1111111111111111D0/
      DATA C/-.5626256053692215D-1,.1014080283006363D-1,&
     -.2365032522738145D+0,-.3575897729251617D-2, .9353769525946207D-1,&
     -.5891279693869842D+0, .1956565409947221D-2,-.5475538688906869D-1,&
      .4158812000823069D+0,-.1136281595717540D+1,-.1436530236370892D-2,&
      .4215852772126871D-1,-.3600995965020568D+0, .1250150711840691D+1,&
     -.1870491772932950D+1, .1271790309026868D-2,-.3876035791590677D-1,&
      .3609622434528460D+0,-.1466884208400427D+1, .2906136259308429D+1,&
     -.2755812719772046D1/
      DATA D/.5626256053692215D-1,.3165475718170830D-2,&
      .2365032522738145D+0, .1780977692217434D-3, .4579298550602792D-1,&
      .5891279693869842D+0, .1002023652232913D-4, .8431857153525702D-2,&
      .2535340690545693D+0, .1136281595717540D+1, .5637641639318209D-6,&
      .1529784002500466D-2, .9783423653244401D-1, .8752546646840911D+0,&
      .1870491772932950D+1, .3171881540176138D-7, .2762930909826477D-3,&
      .3602855398373646D-1, .5767330002770787D+0, .2248588760769160D+1,&
      .2755812719772046D1/
      DATA R/.8065938648381888D1,.3374249976962636D1,&
      .5801001559264062D1, .2037111835358584D1, .2725442211808226D1,&
      .5140624105810932D1, .1475040217560412D1, .1805153580140251D1,&
      .2620644926387035D1, .5345976899871110D1, .1206187666058446D1,&
      .1418278263734739D1, .1877242496186810D1, .2957116017290456D1,&
      .6617662013702422D1, .1085472193938643D1, .1254264622281878D1,&
      .1600266549490816D1, .2323598300219695D1, .4109975778344558D1,&
      .1084602619023685D2/
      IF(.NOT.START) GO TO 2
      NS=0
      NF=0
      NL=NI+4
      PW=W(8)
      DO 1 N=1,7
      BD(N)=ZERO
      DO 1 K=1,NV
      B(N,K)=ZERO
    1 G(N,K)=ZERO
      NSF=.FALSE.
      START=.FALSE.
      GO TO 6
    2 Q=DABS(TP/T)
      DO 5 K=1,NV
      IF(NS.EQ.1) GO TO 4
      DO 3 J=1,7
    3 BD(J)=B(J,K)-E(J,K)
    4 E(1,K)=    Q*(B(1,K)+2.D0*B(2,K) +3.D0*B(3,K)+4.D0*B(4,K)+5.D0*B(5,K) +6.D0*B(6,K) +7.D0*B(7,K))
      E(2,K)=             Q**2*(B(2,K) +3.D0*B(3,K)+6.D0*B(4,K)+10.D0*B(5,K)+15.D0*B(6,K)+21.D0*B(7,K))
      E(3,K)=                          Q**3*(B(3,K)+4.D0*B(4,K)+10.D0*B(5,K)+20.D0*B(6,K)+35.D0*B(7,K))
      E(4,K)=Q**4*(B(4,K)+ 5.D0*B(5,K)+15.D0*B(6,K)+35.D0*B(7,K))
      E(5,K)=             Q**5*(B(5,K)+ 6.D0*B(6,K)+21.D0*B(7,K))
      E(6,K)=                          Q**6*(B(6,K)+ 7.D0*B(7,K))
      E(7,K)=                                        Q**7*B(7,K)
      DO 5 L=1,7
    5 B(L,K)=E(L,K)+BD(L)
      NL=NI
    6 DIR=TP/DABS(TP)
      CALL FORCE(TM,X,F1)
      NF=NF+1
      IF(NS.EQ.0) GO TO 8
      DO 7 K=1,NV
      G(1,K)=B(1,K)+D(1)*B(2,K)+D(2)*B(3,K)+D(4)*B(4,K)+D( 7)*B(5,K)+D(11)*B(6,K)+D(16)*B(7,K)
      G(2,K)=            B(2,K)+D(3)*B(3,K)+D(5)*B(4,K)+D( 8)*B(5,K)+D(12)*B(6,K)+D(17)*B(7,K)
      G(3,K)=B(3,K)+D(6)*B(4,K)+D( 9)*B(5,K)+D(13)*B(6,K)+D(18)*B(7,K)
      G(4,K)=            B(4,K)+D(10)*B(5,K)+D(14)*B(6,K)+D(19)*B(7,K)
      G(5,K)=                         B(5,K)+D(15)*B(6,K)+D(20)*B(7,K)
      G(6,K)=                                      B(6,K)+D(21)*B(7,K)
    7 G(7,K)=                                                   B(7,K)
    8 T=TP
      TVAL=DABS(T)
      DO 26 M=1,NL
      DO 25 J=2,8
      JD=J-1
      S=H(J)
      DO 9 K=1,NV
      A=W(3)*B(3,K)+S*(W(4)*B(4,K)+S*(W(5)*B(5,K)+S*(W(6)*B(6,K)+S*W(7)*B(7,K))))
      Y(K)=X(K)+T*S*(F1(K)+S*(W(1)*B(1,K)+S*(W(2)*B(2,K)+S*A)))
    9 CONTINUE
      CALL FORCE(TM+S*T,Y,FJ)
      NF=NF+1
      DO 24 K=1,NV
      TEMP=G(JD,K)
      GK=(FJ(K)-F1(K))/S
      GO TO (10,11,12,13,14,15,16),JD
   10 G(1,K)=      GK
      GO TO 17
   11 G(2,K)=     (GK-G(1,K))*R(1)
      GO TO 17
   12 G(3,K)=    ((GK-G(1,K))*R(2)-G(2,K))*R(3)
      GO TO 17
   13 G(4,K)=   (((GK-G(1,K))*R(4)-G(2,K))*R(5)-G(3,K))*R(6)
      GO TO 17
   14 G(5,K)=  ((((GK-G(1,K))*R(7)-G(2,K))*R(8)-G(3,K))*R(9)-G(4,K))*R(10)
      GO TO 17
   15 G(6,K)= (((((GK-G(1,K))*R(11)-G(2,K))*R(12)-G(3,K))*R(13)-G(4,K))*R(14)-G(5,K))*R(15)
      GO TO 17
   16 G(7,K)=((((((GK-G(1,K))*R(16)-G(2,K))*R(17)-G(3,K))*R(18)-G(4,K))*R(19)-G(5,K))*R(20)-G(6,K))*R(21)
   17 TEMP=G(JD,K)-TEMP
      B(JD,K)=B(JD,K)+TEMP
      GO TO (24,24,18,19,20,21,22,23),J
   18 B(1,K)=B(1,K)+C(1)*TEMP
      GO TO 24
   19 B(1,K)=B(1,K)+C(2)*TEMP
      B(2,K)=B(2,K)+C(3)*TEMP
      GO TO 24
   20 B(1,K)=B(1,K)+C(4)*TEMP
      B(2,K)=B(2,K)+C(5)*TEMP
      B(3,K)=B(3,K)+C(6)*TEMP
      GO TO 24
   21 B(1,K)=B(1,K)+C(7)*TEMP
      B(2,K)=B(2,K)+C(8)*TEMP
      B(3,K)=B(3,K)+C(9)*TEMP
      B(4,K)=B(4,K)+C(10)*TEMP
      GO TO 24
   22 B(1,K)=B(1,K)+C(11)*TEMP
      B(2,K)=B(2,K)+C(12)*TEMP
      B(3,K)=B(3,K)+C(13)*TEMP
      B(4,K)=B(4,K)+C(14)*TEMP
      B(5,K)=B(5,K)+C(15)*TEMP
      GO TO 24
   23 B(1,K)=B(1,K)+C(16)*TEMP
      B(2,K)=B(2,K)+C(17)*TEMP
      B(3,K)=B(3,K)+C(18)*TEMP
      B(4,K)=B(4,K)+C(19)*TEMP
      B(5,K)=B(5,K)+C(20)*TEMP
      B(6,K)=B(6,K)+C(21)*TEMP
   24 CONTINUE
   25 CONTINUE
   26 CONTINUE
      HV=ZERO
      !DO 27 K=1,NV
      DO 27 K=1,6
   27 HV=DMAX1(HV,DABS(B(7,K)))
      TVAL1=TVAL**3
      HV=HV*W(7)/TVAL1/TVAL1/TVAL
      IF(NSF) GO TO 29
      TP=SS**PW/HV**PW*DIR
      IF(TP/T.GT.ONE) GO TO 28
      TP=0.8D0*TP
      GO TO 8
   28 NSF=.TRUE.
   29 DO 30 K=1,NV
      X(K)=X(K)+T*(F1(K)+W(1)*B(1,K)+W(2)*B(2,K)+W(3)*B(3,K)+W(4)*B(4,K)+W(5)*B(5,K)+W(6)*B(6,K)+W(7)*B(7,K))
   30 CONTINUE
      TM=TM+T
      NS=NS+1
      TP=DIR*SS**PW/HV**PW
      IF(TP/T.GT.SR) TP=T*SR

    RETURN
    END

