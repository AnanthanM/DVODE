c***********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c Maximum Syntax Length is 72
c***********************************************************************

c***********************************************************************
c***********************************************************************
C     THIN FILM EQUATION SOLVER
c
C This code solves the non dimensionalised Thin film Equation in one
C dimension ( Eqn 6 in "Pattern Formation in Spontaneous Dewetting of
C Thin Apolar Films" - R.Khanna and A.Sharma(IITK) ) of length L(which
C from Dispersion relation has a critical length of 2*PI).
c Successive central differencing in space on N nodes in L.Numerical
c experiments have shown 30 nodes are enough for grid independency so 
c we have to solve 30 O.D.E's with time as the independent variable.
c O.D.E's were integrated in time using VODE Solver in double precision.
c (Brown,Byrne and Hindmarsh - SIAM J.Sci.Stat.Comput(1989) )
c This particular solver was preferred because of stiffness introduced
c by widely seperated length and time scales,especially close to breakup
c time.
c 
c*********************************************************************** 
c
c NOTE ON VODE SOLVER
c
c 
c
c
c







c234567
      PROGRAM MAIN
      implicit none

      EXTERNAL  F !Subroutine which supplies ydot(i) 
      DOUBLE PRECISION ATOL,RPAR,RTOL,RWORK,T,TOUT,Y,PI,L,DELX,XAXIS,
     +                 TIMESCALE,J  
      INTEGER IWORK,ITOL,NEQ,ITASK,ISTATE,IOPT,LRW,LIW,MF,IOUT,IPAR,
     +   A,D,E,G
      PARAMETER   ( PI = 3.14159,NEQ = 30,L = 2*SQRT(2.0)*PI )
      
      DIMENSION Y(NEQ), ATOL(NEQ), RWORK(22 + (9*NEQ) + (2*NEQ**2)),
     + IWORK(30 + NEQ),XAXIS(NEQ),TIMESCALE(10)

      DELX = L/(NEQ-1)

      DO 80 E = 1,NEQ
      Y(E) = 1.0 + 0.01*RAND(0) 

c      Y(E) = 1 +  0.01*SIN(DELX*(E-1)/SQRT(2.0))
   80 CONTINUE
      
c      TIMESCALE(1) = 0.0
c      DO 2000 G = 2,200
c        TIMESCALE(G) = TIMESCALE(G-1) + 0.1 
c 2000 CONTINUE
      
      
      TIMESCALE(1) = 0.0 
      TIMESCALE(2) = 0.01
      TIMESCALE(3) = 5.0
      TIMESCALE(4) = 10.0
      TIMESCALE(5) = 15.0
      TIMESCALE(6) = 20.0
      TIMESCALE(7) = 25.0
      TIMESCALE(8) = 27.0
      TIMESCALE(9) = 28.0
      TIMESCALE(10) = 28.16033077

      T = 0.0D0
      TOUT = TIMESCALE(1)
      ITOL = 2
      RTOL = 1.0D-4
   
      DO 70  D = 1,NEQ
      ATOL(D) = 1.0D-4
   70 CONTINUE  


      ITASK = 1
      ISTATE = 1
      IOPT = 1
      LRW = (22 + (9*NEQ) + (2*NEQ**2))
      LIW = 30 + NEQ
      MF = 22

      DO 50 A = 5,10
      RWORK(A) = 0.0
   50 IWORK(A) = 0    
      
      IWORK(6) = 4000

      DO 91 G = 1,NEQ
        XAXIS(G) = (G - 1)*DELX
   91 CONTINUE
       
      DO 92 G = 1,NEQ
         WRITE(*,'(F14.6)',ADVANCE = 'NO') XAXIS(G)
   92 CONTINUE
      
      WRITE(*,*) ' ' 
      
      DO 10 IOUT = 2,11
        CALL DVODE(F,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,
     +             IOPT,RWORK,LRW,IWORK,LIW,J,MF,RPAR,IPAR) 
c        WRITE(*,'(D14.6)') T
        DO 90 G = 1,NEQ
            
           WRITE(*,'(F14.6)',ADVANCE = 'NO') Y(G)
           
           
   90   CONTINUE
        WRITE(*,*) ' '  
c        WRITE(*,*) IOUT 
c        IF (ISTATE .LT. 0 ) GO TO 30
   10   TOUT = TIMESCALE(IOUT)
c   30 WRITE(*,40) ISTATE 
c   40 FORMAT(' ERROR ISTATE = ' , I3 )
      STOP
      END

      SUBROUTINE F(NEQ,T,Y,YDOT,RPAR,IPAR)
      DOUBLE PRECISION RPAR,T,Y,YDOT,DELX,L,PI
      INTEGER NEQ,IPAR,C
      DIMENSION Y(NEQ),YDOT(NEQ)
      PARAMETER  (PI = 3.14159,L = 2*SQRT(2.0)*PI)
      
      DELX = L/(NEQ-1)



      YDOT(1) = -(
     +  ( ( ( ( (Y(1+1)+Y(1))/2 )**3) )*( (P(1+1,Y) - P(1,Y))/DELX ) )
     +                   -
     +  ( ( ( ( (Y(1)+Y(NEQ))/2 )**3) )*( (P(1,Y) - P(NEQ,Y))/DELX ) )
     +          )/DELX
      
      DO 60 C = 2,( NEQ-1)
      YDOT(C) = -(
     +  ( ( ( ( (Y(C+1)+Y(C))/2 )**3) )*( (P(C+1,Y) - P(C,Y))/DELX ) )
     +                   -
     +  ( ( ( ( (Y(C)+Y(C-1))/2 )**3) )*( (P(C,Y) - P(C-1,Y))/DELX ) )
     +          )/DELX
   60 CONTINUE  
      YDOT(NEQ) = -(
     +  ( ( ( ( (Y(1)+Y(NEQ))/2 )**3) )*( (P(1,Y) - P(NEQ,Y))/DELX ) )
     +                   -
     +  ( ( ( ( (Y(NEQ)+Y(NEQ-1))/2 )**3) )*
     +  ( (P(NEQ,Y) - P(NEQ-1,Y))/DELX ) )
     +          )/DELX
      RETURN
      END  

      REAL*4 FUNCTION P(I,Y)
      REAL*8 Y(*),L,DELX,PI
      INTEGER I,B,NEQ
      NEQ = 30
      PI = 3.14159
      L = 2*SQRT(2.0)*PI
      
      DELX = L/(NEQ-1)

      B = 0

        
      IF (I.EQ.(1)) THEN
       P = ( ( Y(2)-2*Y(1)+Y(NEQ)  ) / (DELX)**2  )  +  B*Y(I)  -
     +    (  (1/( 3*( (Y(I))**3 ) ))*(1-(0.137/Y(I))**6) )

      ELSEIF (I.EQ.(NEQ)) THEN
        P = ( ( Y(1)-2*Y(I)+Y(I-1)  ) / (DELX)**2  )  +  B*Y(I)  -
     +    (  (1/( 3*( (Y(I))**3 ) ))*(1-(0.137/Y(I))**6) )
          
      ELSE
       P = ( ( Y(I+1)-2*Y(I)+Y(I-1)  ) / (DELX)**2  )  +  B*Y(I)  -
     +    (  (1/( 3*( (Y(I))**3 ) ))*(1-(0.137/Y(I))**6) )
      
      ENDIF
      RETURN
 
      END
          

      
        
      


