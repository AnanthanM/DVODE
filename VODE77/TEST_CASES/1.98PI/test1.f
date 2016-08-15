c234567
      PROGRAM MAIN
      IMPLICIT NONE

      EXTERNAL  F,JEX !Subroutine which supplies ydot(i) 
      DOUBLE PRECISION ATOL,RPAR,RTOL,RWORK,T,TOUT,Y,PI,L,DELX,XAXIS,
     +                 TIMESCALE,J  
      INTEGER IWORK,ITOL,NEQ,ITASK,ISTATE,IOPT,LRW,LIW,MF,IOUT,IPAR,
     +   I,ML,MU
      PARAMETER   ( PI = 3.141592654,NEQ = 60,
c     + L = 2*SQRT(2.0)*PI,
     + L = 1.98*PI,
     +  ML = 2, MU = 2  )
      
      DIMENSION Y(NEQ),ATOL(NEQ), RWORK(22+(11*NEQ)+((3*ML+2*MU)*NEQ)),
     + IWORK(30 + NEQ),XAXIS(NEQ),TIMESCALE(20)

      DELX = L/(NEQ-1)

      IWORK(1) = 2
      IWORK(2) = 2

      DO 100 I = 1,NEQ
        Y(I) = 1.0 + 0.01*RAND(0) 
  100 CONTINUE
      
     
c L = 2sqrt(2)PI      
c      TIMESCALE(1) = 0.0 
c      TIMESCALE(2) = 0.01
c      TIMESCALE(3) = 5.0
c      TIMESCALE(4) = 10.0
c      TIMESCALE(5) = 15.0
c      TIMESCALE(6) = 25.0
c      TIMESCALE(7) = 28.0
c      TIMESCALE(8) = 28.1
c      TIMESCALE(9) = 28.2
c      TIMESCALE(10) = 28.5
c
c      TIMESCALE(11) = 28.1
c      TIMESCALE(12) = 28.15
c      TIMESCALE(13) = 28.2
c      TIMESCALE(14) = 28.3
c      TIMESCALE(15) = 28.5
c      TIMESCALE(16) = 28.5
c      TIMESCALE(17) = 28.5
c      TIMESCALE(18) = 28.5
c      TIMESCALE(19) = 28.5
c      TIMESCALE(20) = 28.5
 

c L = 2.01PI

c      TIMESCALE(1) = 0.0 
c      TIMESCALE(2) = 0.01
c      TIMESCALE(3) = 0.02
c      TIMESCALE(4) = .4
c      TIMESCALE(5) = 20.0
c      TIMESCALE(6) = 100.0
c      TIMESCALE(7) = 110
c      TIMESCALE(8) = 120
c      TIMESCALE(9) = 125
c      TIMESCALE(10) = 130
c
c      TIMESCALE(11) = 150 
c      TIMESCALE(12) = 160
c      TIMESCALE(13) = 170
c      TIMESCALE(14) = 180
c      TIMESCALE(15) = 190
c      TIMESCALE(16) = 191
c      TIMESCALE(17) = 192
c      TIMESCALE(18) = 195
c      TIMESCALE(19) = 197
c      TIMESCALE(20) = 198.17

c L = 1.98PI
c
      TIMESCALE(1) = 0.0 
      TIMESCALE(2) = 0.01
      TIMESCALE(3) = 0.03
      TIMESCALE(4) = .2
      TIMESCALE(5) = 1
      TIMESCALE(6) = 10.0
      TIMESCALE(7) = 50.0
      TIMESCALE(8) = 150.0
      TIMESCALE(9) = 200.0
      TIMESCALE(10) = 500.0

      T = 0.0D0
      TOUT = TIMESCALE(1)
      ITOL = 2

      RTOL = 1.0D-4
   
      DO 101  I = 1,NEQ
        ATOL(I) = 1.0D-4
  101 CONTINUE  

      ITASK = 1
      ISTATE = 1
      IOPT = 1
      LRW = (22+(11*NEQ)+((3*ML+2*MU)*NEQ)) 
      LIW = 30 + NEQ
      MF = 24

      DO 102 I = 5,10
        RWORK(I) = 0.0
        IWORK(I) = 0    
  102 CONTINUE

      IWORK(6) = 4000000000000

      DO 103 I = 1,NEQ
        XAXIS(I) = (I - 1)*DELX
  103 CONTINUE
       
      DO 104 I = 1,NEQ
         WRITE(*,'(F14.6)',ADVANCE = 'NO') XAXIS(I)
  104 CONTINUE
      
      WRITE(*,*) ' ' 
      
      DO 105 IOUT = 2,11
        CALL DVODE(F,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,
     +             IOPT,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR) 
        DO 106 I = 1,NEQ
          WRITE(*,'(F14.6)',ADVANCE = 'NO') Y(I)       
  106   CONTINUE
        WRITE(*,*) ' '  
c        WRITE(*,*) IOUT 
c        IF (ISTATE .LT. 0 ) GO TO 30
        TOUT = TIMESCALE(IOUT)
  105 CONTINUE
c   30 WRITE(*,40) ISTATE 
c   40 FORMAT(' ERROR ISTATE = ' , I3 )
      STOP
      END

      SUBROUTINE F(NEQ,T,Y,YDOT,RPAR,IPAR)
      DOUBLE PRECISION RPAR,T,Y,YDOT,X,L,PI,pY
      INTEGER NEQ,IPAR,I
      DIMENSION Y(NEQ),YDOT(NEQ) ,pY(NEQ + 4)
      PARAMETER  (PI = 3.141592654,
c     + L = 2*SQRT(2.0)*PI)
     + L = 1.98*PI)
      X = L/(NEQ-1)
      
      DO 1001 I = 1,NEQ
        pY(I + 2) = Y(I)
 1001 CONTINUE       
 
      pY(NEQ + 3) = pY(3) 
      pY(NEQ + 4) = pY(4)
      pY(1) = pY(NEQ + 1)
      pY(2) = pY(NEQ + 2)

      DO 1002 I = 1,NEQ

        YDOT(I) = -((((((pY(I+2+1)+pY(I+2))/2)**3)*                   
     +                      ((P(I+1,NEQ,Y)-P(I+0,NEQ,Y))/X))         
     +           -    ((((pY(I+2)+pY(I+2-1))/2)**3)*                  
     +                       ((P(I+0,NEQ,Y)-P(I-1,NEQ,Y))/X))        
     +           )/X) 

 1002 CONTINUE    

      RETURN
      END  

      REAL*4 FUNCTION P(I,NEQ,Y)
      REAL*8 Y(NEQ),L,X,PI,sL,pY(NEQ+4)
      INTEGER I,NEQ,C
      PI = 3.141592654
      sL = .02
c      L = 2*SQRT(2.0)*PI
      L = 1.98*PI
      X = L/(NEQ-1)
      
      DO 1003 C = 1,NEQ
        pY(C + 2) = Y(C)
 1003 CONTINUE       
       
      pY(NEQ + 3) = pY(3) 
      pY(NEQ + 4) = pY(4)
      pY(1) = pY(NEQ + 1)
      pY(2) = pY(NEQ + 2)
      
      P = ((pY(I+1+2)-2*pY(I+2)+pY(I-1+2))/X**2) -                  
     +       ((1-(sL/pY(I+2))**6)/(3*(pY(I+2))**3))
        
      RETURN 
      END
     

      SUBROUTINE JEX (NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
      DOUBLE PRECISION PD, RPAR, T, Y ,PY, X , L , PI , BL
      DIMENSION Y(NEQ), PD(NRPD,NEQ) , PY(NEQ+4)
      PARAMETER ( PI = 3.141592654,
c     + L = 2*SQRT(2.0)*PI )  
     + L = 1.98*PI )  
      X = L/(NEQ-1)
      BL = 0.02

      DO 10001  I = 1,NEQ
        PY(I+2) = Y(I)
10001 CONTINUE

      PY(NEQ+3) = PY(3)
      PY(NEQ+4) = PY(4)
      PY(1)     = PY(NEQ+1)
      PY(2)     = PY(NEQ+2)
      
      DO 10002 I = 1,NEQ 
      PD(I,I) = (
     +           (-1/(8*(X**2)))*
     +           (
     +           (
     +           ((PY(I+2)+PY(I+3))**3)*
     +     ((3/(X**2))-(1/((PY(I+2)**4)))+((3*(BL**6))/(PY(I+2)**10))) 
     +           ) 
     +           +
     +           (
     +           (3*((PY(I+2)+PY(I+3))**2))*
     +           (
     +           (
     +           ((PY(I+4)-2*PY(I+3)+PY(I+2))/(X**2)) -
     +           ((1/(3*(PY(I+3)**3)))*(1-((BL/PY(I+3))**6))) 
     +           ) -
     +           (
     +           ((PY(I+3)-2*PY(I+2)+PY(I+1))/(X**2)) -
     +           ((1/(3*(PY(I+2)**3)))*(1-((BL/PY(I+2))**6))) 
     +           ) 
     +           )
     +           )
     +           -
     +           (
     +           ((PY(I+2)+PY(I+1))**3)*
     +     ((-3/(X**2))+(1/((PY(I+2)**4)))-((3*(BL**6))/(PY(I+2)**10))) 
     +           ) 
     +           -
     +           (
     +           (3*((PY(I+2)+PY(I+1))**2))*
     +           (
     +           (
     +           ((PY(I+3)-2*PY(I+2)+PY(I+1))/(X**2)) -
     +           ((1/(3*(PY(I+2)**3)))*(1-((BL/PY(I+2))**6))) 
     +           ) -
     +           (
     +           ((PY(I+2)-2*PY(I+1)+PY(I+0))/(X**2)) -
     +           ((1/(3*(PY(I+1)**3)))*(1-((BL/PY(I+1))**6))) 
     +           ) 
     +           )
     +           )
     +           )
     +           )

10002 CONTINUE
      
      DO 10003 I =1,(NEQ-1)
      
        PD(I+1,I) = (      
     +          (-1/(8*(X**2)))*
     +          (
     +          ((-((PY(I+4)+PY(I+3))**3))/(X**2))
     +          -
     +          (
     +          (3*((PY(I+3)+PY(I+2))**2))*         
     +          (
     +          (
     +          ((PY(I+4)-2*PY(I+3)+PY(I+2))/(X**2))-
     +          ((1/(3*(PY(I+3)**3)))*(1-((BL/PY(I+3))**6))) 
     +          ) - 
     +          (
     +          ((PY(I+3)-2*PY(I+2)+PY(I+1))/(X**2))-
     +          ((1/(3*(PY(I+2)**3)))*(1-((BL/PY(I+2))**6))) 
     +          )  
     +          )
     +          ) 
     +          -
     +          (
     +           ((PY(I+3)+PY(I+2))**3)*
     +     ((3/(X**2))-(1/((PY(I+2)**4)))+((3*(BL**6))/(PY(I+2)**10))) 
     +          ) 
     +          )   
     +          )
 
10003 CONTINUE   

      DO 10004 I = 1,(NEQ-2)

        PD(I+2,I) = (
     +          (-1/(8*(X**2)))*
     +          (((PY(I+4)+PY(I+3))**3)/(X**2))
     +              )
       
10004 CONTINUE      

      DO 10005 I = 2,NEQ

        PD(I-1,I) = (
     +          (-1/(8*(X**2)))*
     +          (   
     +          (
     +           ((PY(I+2)+PY(I+1))**3)*
     +    ((-3/(X**2))+(1/((PY(I+2)**4)))-((3*(BL**6))/(PY(I+2)**10))) 
     +          ) 
     +          +
     +          (
     +          (3*((PY(I+2)+PY(I+1))**2))*
     +          (
     +          (
     +          ((PY(I+3)-2*PY(I+2)+PY(I+1))/(X**2))-
     +          ((1/(3*(PY(I+2)**3)))*(1-((BL/PY(I+2))**6))) 
     +          ) - 
     +          (
     +          ((PY(I+2)-2*PY(I+1)+PY(I+0))/(X**2))-
     +          ((1/(3*(PY(I+1)**3)))*(1-((BL/PY(I+1))**6))) 
     +          )  
     +          )
     +          )
     +          -
     +          (((PY(I+1)+PY(I+0))**3)/(X**2))
     +          )
     +          )   
                  
10005 CONTINUE      
      
      DO 10006 I = 3,NEQ

        PD(I-2,I) = ( 
     +          (-1/(8*(X**2)))*
     +          (((PY(I+1)+PY(I+0))**3)/(X**2))
     +              )

10006 CONTINUE      

      RETURN
      END

C-----------------------------------------------------------------------

