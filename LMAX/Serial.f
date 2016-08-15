c234567
c gfortran  -cpp -dM  -fno-range-check  Parallel.f Solver.f
#define comments     0     
#define comments2     0     
#define NEQ  60 
#define L    1.0*8.885765876
#define PI   3.141592653589793238462643383276   
c#define NO_OF_ITERATIONS  500000 
#define NO_OF_ITERATIONS  5
      PROGRAM MAIN
      
      implicit none 

      EXTERNAL  F,JEX !Subroutine which supplies ydot(i) 

      DOUBLE PRECISION , DIMENSION(:),ALLOCATABLE::ATOL,RWORK,Y,XAXIS
      INTEGER , DIMENSION(:) , ALLOCATABLE :: IWORK 
      
      DOUBLE PRECISION RPAR,RTOL,T,TOUT,DELX,pYT(NEQ + 4)                 
      INTEGER ITOL,ITASK,ISTATE,IOPT,LRW,LIW,MF,IOUT,IPAR,
     +   I,ML,MU,NEQL
 
      common / PERIODIC_Y / pYT

      PARAMETER   (   ML = 2, MU = 2  )
      
      ALLOCATE ( Y(NEQ) ,ATOL(NEQ), XAXIS(NEQ),
     + RWORK(22+(11*NEQ)+((3*ML+2*MU)*NEQ)),
     + IWORK(30 + NEQ))

      DELX = L/(NEQ-1)
      NEQL = NEQ 
      IWORK(1) = 2
      IWORK(2) = 2

      DO 100 I = 1,NEQ
        Y(I) = 1.0 + 0.01*RAND(0) 
  100 CONTINUE
      
     

      T = 0.0D0
      TOUT = 1.0
      ITOL = 2

      RTOL = 1.0D-4
   
      DO 101  I = 1,NEQ
        ATOL(I) = 1.0D-4
  101 CONTINUE  

      ITASK = 2
      ISTATE = 1
      IOPT = 1
      LRW = (22+(11*NEQ)+((3*ML+2*MU)*NEQ)) 
      LIW = 30 + NEQ
      MF = 24

      DO 102 I = 5,10
        RWORK(I) = 0.0
        IWORK(I) = 0    
  102 CONTINUE

      RWORK(5) = 1.0D-7
      IWORK(6) = 40000000000

      DO 103 I = 1,NEQ
        XAXIS(I) = (I - 1)*DELX
  103 CONTINUE
       
      DO 104 I = 1,NEQ
         WRITE(*,'(F14.6)',ADVANCE = 'NO') XAXIS(I)
  104 CONTINUE
      
      WRITE(*,*) ' ' 
      
      DO 105 IOUT = 1,NO_OF_ITERATIONS

        CALL DVODE(F,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,
     +             IOPT,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR) 
        print *,'At time ', T , 'at ' , IOUT 
        DO 106 I = 1,NEQ
          WRITE(*,'(F14.6)',ADVANCE = 'NO') Y(I)       
  106   CONTINUE
        WRITE(*,*) ' '  
c        print *, "The HNEXT IS", RWORK(12)

  105 CONTINUE
  
      STOP 
      END

      SUBROUTINE F(NEQL,T,Y,YDOT,RPAR,IPAR)
      implicit none 
      double precision, external :: P
      DOUBLE PRECISION RPAR,T,Y,YDOT,X,pY,pYT
      INTEGER IPAR,I,NEQL 
      DIMENSION Y(NEQ),YDOT(NEQ),pY(NEQ + 4),pYT(NEQ + 4)
      
      common / PERIODIC_Y / pYT

      X = L/(NEQ-1)

      DO 1001 I = 1,NEQ
        pY(I + 2) = Y(I)
 1001 CONTINUE       
 
      pY(NEQ + 3) = pY(3)
      pY(NEQ + 4) = pY(4)
      pY(1) = pY(NEQ + 1) 
      pY(2) = pY(NEQ + 2) 

      DO 1004 I = 1,(NEQ+4)
        pYT(I) = pY(I) 
 1004 CONTINUE      


      print*," Y(10) is ", Y(10)

      DO 1002 I = 1,NEQ

        YDOT(I) = -((((((pY(I+2+1)+pY(I+2))/2)**3)*                   
     +                      ((P(I+1,NEQ,Y)-P(I+0,NEQ,Y))/X))         
     +           -    ((((pY(I+2)+pY(I+2-1))/2)**3)*                  
     +                       ((P(I+0,NEQ,Y)-P(I-1,NEQ,Y))/X))        
     +           )/X) 

 1002 CONTINUE    

      RETURN
      END  

      DOUBLE PRECISION  FUNCTION P(I,NEQL,Y)
     
      implicit none 
      DOUBLE PRECISION  Y(NEQ),X,sL,pY(NEQ+4)
     + ,pYT(NEQ+4)
      INTEGER I,C,NEQL
      common / PERIODIC_Y / pYT
      sL = .02
      X = L/(NEQ-1)
      
      DO 1003 C = 1,(NEQ+4)
        pY(C) = pYT(C)
 1003 CONTINUE       
       
      
      P = ((pY(I+1+2)-2*pY(I+2)+pY(I-1+2))/X**2) -                  
     +       ((1-(sL/pY(I+2))**6)/(3*(pY(I+2))**3))
        
      RETURN 
      END
     

      SUBROUTINE JEX (NEQL, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
      
      DOUBLE PRECISION PD, RPAR, T, Y ,PY, X ,  BL
     + ,pYT(NEQ+4) 
      INTEGER I, NEQL 
     + ,IPAR,MU,ML,NRPD 
      
      DIMENSION Y(NEQ), PD(NRPD,NEQ) , PY(NEQ+4)
      common / PERIODIC_Y / pYT
      
      X = L/(NEQ-1)
      BL = 0.02

      DO 10001  I = 1,(NEQ+4)
        PY(I) = pYT(I)
10001 CONTINUE

      print*," from Jac Y(10) is ", Y(10)
      
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

