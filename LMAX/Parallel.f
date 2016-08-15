c mpif90  -cpp -dM  -fno-range-check  Parallel.f Solver.f
c234567
#define commentsF     1     
#define test     1     
#define NEQ  1*60 
#define L    1*8.885765876
#define PI   3.141592653589793238462643383276  
#define petest  0      
#define testF  1      
      PROGRAM MAIN
       
      implicit none 

      include 'mpif.h'
      EXTERNAL  F,JEX !Subroutine which supplies ydot(i) 
      DOUBLE PRECISION , DIMENSION(:),ALLOCATABLE::ATOL,RWORK,Y,XAXIS

      INTEGER , DIMENSION(:) , ALLOCATABLE :: IWORK 
      
      DOUBLE PRECISION RPAR,RTOL,T,TOUT,DELX
     +  ,HCUR,TNEXT  
     + ,pYT(NEQ+4)
      INTEGER ITOL,ITASK,ISTATE,IOPT,LRW,LIW,MF,IOUT,IPAR,
     +   I,ML,MU,NEQL,mype,npes,ierror,neqperpe,nrem,mybase 
      common / PERIODIC_Y / pYT
      common / ParallelInfo / mype , npes
      PARAMETER   (   ML = 2, MU = 2  )
      
      call MPI_INIT(ierror)

      call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierror)

      call MPI_COMM_SIZE(MPI_COMM_WORLD,npes,ierror)
      
      neqperpe = NEQ/npes

      nrem = NEQ  - npes*neqperpe

      if(mype.LT.nrem) then
        NEQL = neqperpe + 1
      else 
        NEQL = neqperpe
      endif


      if(mype.LT.nrem) then
        mybase = mype*NEQL + 1
      else 
        mybase = mype*neqperpe + nrem + 1
      endif
      ALLOCATE ( Y(NEQL) ,ATOL(NEQL), XAXIS(NEQL),
     + RWORK(22+(11*NEQL)+((3*ML+2*MU)*NEQL)),
     + IWORK(30 + NEQL))

      DELX = L/(NEQ-1)
      IWORK(1) = 2
      IWORK(2) = 2
c CHANGE THE RANDOM PERTUBATION FOR EACH PROCESSOR
      DO 100 I = 1,NEQL
        Y(I) = 1.0 + 0.01*RAND(0) 
  100 CONTINUE
      
     

      T = 0.0D0
      TOUT = 10
      ITOL = 2

      RTOL = 1.0D-4
   
      DO 101  I = 1,NEQL
        ATOL(I) = 1.0D-4
  101 CONTINUE  

      ITASK = 2
      ISTATE = 1
      IOPT = 1
      LRW = (22+(11*NEQL)+((3*ML+2*MU)*NEQL)) 
      LIW = 30 + NEQL
      MF = 24

      DO 102 I = 5,10
        RWORK(I) = 0.0
        IWORK(I) = 0    
  102 CONTINUE

      IWORK(6) = 4.0D10
      RWORK(5) = 1.0D-7

      DO 10300 I = 1,NEQL 
        XAXIS(I) = (mybase - 1)*DELX + (I - 1)*DELX 
10300 CONTINUE 
      

      if(mype.eq.petest ) then 
      print*,"From Procees", mype  
      DO 104 I = 1,NEQL
         WRITE(*,'(F14.6)',ADVANCE = 'NO') XAXIS(I)
  104 CONTINUE
      
      WRITE(*,*) ' ' 
      WRITE(*,*) ' ' 
      WRITE(*,*) ' ' 
      
      print*,"From Procees", mype  
      DO 10610 I = 1,NEQL
        WRITE(*,'(F14.6)',ADVANCE = 'NO') Y(I)       
10610 CONTINUE
      WRITE(*,*) ' '  
      WRITE(*,*) ' '  
      WRITE(*,*) ' '  
      endif 

      DO 105 IOUT = 1,4
      

        CALL DVODE(F,NEQL,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,
     +             IOPT,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR) 
        
      if(mype.eq.petest) then 
        print *,'At time ', T , 'at ' , IOUT ,"from ", mype 
        DO 106 I = 1,NEQL
          WRITE(*,'(F14.10)',ADVANCE = 'NO') Y(I)       
  106   CONTINUE
        WRITE(*,*) ' '  
        WRITE(*,*) ' '  
        WRITE(*,*) ' '  

        print* , 'This is from', mype ,'process with
     +   RWORK',RWORK(12)      
      endif 
  105 CONTINUE

      call MPI_FINALIZE(ierror) 
      STOP 
      END

      SUBROUTINE F(NEQL ,T,Y,YDOT,RPAR,IPAR)
      implicit none 
      include 'mpif.h'
      double precision, external :: P
      DOUBLE PRECISION RPAR,T,Y,YDOT,X,pY,
     + pYT(NEQ+4) 
      INTEGER IPAR,I,NEQL,mype,npes,lastpe,ierr,
     +  ierror,mypeleft,myperight  
      integer status(MPI_STATUS_SIZE)
      DIMENSION Y(NEQL),YDOT(NEQL) ,pY(NEQL + 4)
      
      common / PERIODIC_Y / pYT 

      common / ParallelInfo / mype , npes

      X = L/(NEQ-1)

#if commentsF      

      if(mype.eq.petest ) then 
      print* , 'This is from inside F beginning in ', mype , 'process'  
      endif 
      CALL FLUSH(6)
#endif 
      lastpe = npes - 1
#if commentsF      

      if(testF.eq.1 ) then 
      print* , 'This is from inside F beginning in ', mype , 'process'  
      endif 
      CALL FLUSH(6)
#endif 
      
      DO 10010 I = 1,NEQL
        pY(I + 2) = Y(I)
10010 CONTINUE       
      
        
        if(mype.NE.lastpe) then   
          myperight  = mype + 1 
        else
          myperight = 0 
        endif 
  
        if(mype.NE.0) then 
          mypeleft  = mype - 1 
        else
          mypeleft = lastpe 
        endif 
  
#if commentsF      

      if(mype.eq.petest ) then 
      print* , 'This is from inside F line 183 in ', mype , 'process'  
      endif 
      CALL FLUSH(6)
#endif 
  
#if commentsF      

      if(testF.eq.1 ) then 
      print* , 'This is from inside F line 183 in ', mype , 'process'  
      endif 
      CALL FLUSH(6)
#endif 
  
  
        call MPI_SEND( pY(3), 2, MPI_DOUBLE_PRECISION, mypeleft, 
     +                   3024, MPI_COMM_WORLD, ierr ) 
  
        call MPI_RECV(pY(NEQL+3), 2, MPI_DOUBLE_PRECISION, myperight, 
     +                  3024, MPI_COMM_WORLD, status, ierr ) 
  
        call MPI_SEND( pY(NEQL+1), 2, MPI_DOUBLE_PRECISION, myperight, 
     +                   3512, MPI_COMM_WORLD, ierr ) 
  
  
        call MPI_RECV(pY(1), 2, MPI_DOUBLE_PRECISION, mypeleft, 
     +                  3512, MPI_COMM_WORLD, status, ierr ) 
  
  
#if commentsF      

      if(testF.eq.1 ) then 
      print* , 'This is from inside F line 221 in ', mype , 'process'  
      endif 
      CALL FLUSH(6)
#endif 
       
      do 629 I = 1,(NEQL +4)
         pYT(I) = pY(I)
  629 continue 





      DO 1002 I = 1,NEQL

        YDOT(I) = -((((((pY(I+2+1)+pY(I+2))/2)**3)*                   
     +                      ((P(I+1,NEQ,Y)-P(I+0,NEQ,Y))/X))         
     +           -    ((((pY(I+2)+pY(I+2-1))/2)**3)*                  
     +                       ((P(I+0,NEQ,Y)-P(I-1,NEQ,Y))/X))        
     +           )/X) 

 1002 CONTINUE    
#if commentsF

      if(mype.eq.petest ) then 
        print*, 'From F end from ', mype ,'process' 
      endif 
#endif

#if commentsF

      if(testF.eq.1 ) then 
        print*, 'From F end from ', mype ,'process' 
      endif 
#endif        
      RETURN
      END  

      DOUBLE PRECISION  FUNCTION P(I,NEQL,Y)
     
      implicit none 
      include 'mpif.h'
      DOUBLE PRECISION  Y(NEQL),X,sL,pY(NEQL+4)
     + ,pY1,pY2,pY3,pY4
     + ,pYT(NEQ+4) 
      INTEGER I,C,NEQL,mype,ierror,npes 
      common / PERIODIC_Y / pYT 
      common / ParallelInfo / mype , npes
      sL = .02
      X = L/(NEQ-1)
      
      call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierror)
       
      DO 1003 C = 1,NEQL+4
        pY(C) = pYT(C)
 1003 CONTINUE       
      
      P = ((pY(I+1+2)-2*pY(I+2)+pY(I-1+2))/X**2) -                  
     +       ((1-(sL/pY(I+2))**6)/(3*(pY(I+2))**3))
        
      RETURN 
      END
     

      SUBROUTINE JEX (NEQL, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
      
      include 'mpif.h'
      DOUBLE PRECISION PD, RPAR, T, Y ,PY, X ,  BL
     + ,pYT(NEQ+4)
      INTEGER I, NEQL,mype,ierror,npes 
     + ,IPAR,MU,ML,NRPD 
      
      DIMENSION Y(NEQL), PD(NRPD,NEQL) , PY(NEQL+4)
      common / PERIODIC_Y / pYT
      common / ParallelInfo / mype , npes
      
      X = L/(NEQ-1)
      BL = 0.02

     
      DO 10001  I = 1,NEQL+4
        PY(I) = pYT(I)
10001 CONTINUE
      
      DO 10002 I = 1,NEQL 
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
      
      DO 10003 I =1,(NEQL-1)
      
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

      DO 10004 I = 1,(NEQL-2)

        PD(I+2,I) = (
     +          (-1/(8*(X**2)))*
     +          (((PY(I+4)+PY(I+3))**3)/(X**2))
     +              )
       
10004 CONTINUE      

      DO 10005 I = 2,NEQL

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
      
      DO 10006 I = 3,NEQL

        PD(I-2,I) = ( 
     +          (-1/(8*(X**2)))*
     +          (((PY(I+1)+PY(I+0))**3)/(X**2))
     +              )

10006 CONTINUE      

      RETURN
      END
      
C-----------------------------------------------------------------------

