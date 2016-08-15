c234567
#define comments     1     
#define comments2     0     
#define NEQ  10*60 
#define L    10*8.885765876
#define PI   3.141592653589793238462643383276   
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

      IWORK(6) = 40000000000
      RWORK(5) = 1.0D-7

      DO 10300 I = 1,NEQL 
        XAXIS(I) = (mybase - 1)*DELX + (I - 1)*DELX 
10300 CONTINUE 
      print*,"From Procees", mype  
      DO 104 I = 1,NEQL
         WRITE(*,'(F14.6)',ADVANCE = 'NO') XAXIS(I)
  104 CONTINUE
      
      WRITE(*,*) ' ' 
      WRITE(*,*) ' ' 
      WRITE(*,*) ' ' 
      
        DO 10610 I = 1,NEQL
          WRITE(*,'(F14.6)',ADVANCE = 'NO') Y(I)       
10610   CONTINUE
        WRITE(*,*) ' '  

      DO 105 IOUT = 2,3
        print *,'pYT(3) is',pYT(3)
        CALL DVODE(F,NEQL,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,
     +             IOPT,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR) 
        print *,'At time ', T , 'at ' , IOUT ,"from ", mype 
        print *,'pYT(1) is',pYT(1)
        DO 106 I = 1,NEQL
          WRITE(*,'(F14.6)',ADVANCE = 'NO') Y(I)       
  106   CONTINUE
        WRITE(*,*) ' '  
        print *, "The HNEXT IS", RWORK(12)

        print* , 'This is from', mype ,'process with
     +   RWORK',RWORK(12)      
        call MPI_ALLREDUCE(HCUR,TNEXT,1,MPI_DOUBLE_PRECISION,MPI_MIN,
     +   MPI_COMM_WORLD)
  105 CONTINUE

      call MPI_FINALIZE(ierror) 
      STOP 
      END

      SUBROUTINE F(NEQL ,T,Y,YDOT,RPAR,IPAR)
      implicit none 
      include 'mpif.h'
      double precision, external :: P
      DOUBLE PRECISION RPAR,T,Y,YDOT,X,pY,
c     + pY1,pY2,pY3,pY4
     + pYT(NEQ+4) 
      INTEGER IPAR,I,NEQL,mype,npes,lastpe,ierr,
     +  ierror,mypeleft,myperight  
      integer status(MPI_STATUS_SIZE)
      DIMENSION Y(NEQL),YDOT(NEQL) ,pY(NEQL + 4)
      
c      common / PERIODIC_Y / pY1,pY2,pY3,pY4
      common / PERIODIC_Y / pYT 

      X = L/(NEQ-1)
c      print*, 'From F beginning'

      call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierror)

      call MPI_COMM_SIZE(MPI_COMM_WORLD,npes,ierror)
#if comments      

      print* , 'This is from inside F beginning ', mype , 'process'  
      CALL FLUSH(6)
#endif 
      lastpe = npes - 1
      
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
  
  
  
  
        call MPI_SEND( pY(3), 2, MPI_DOUBLE_PRECISION, mypeleft, 
     +                   3024, MPI_COMM_WORLD, ierr ) 
  
        call MPI_RECV(pY(NEQL+3), 2, MPI_DOUBLE_PRECISION, myperight, 
     +                  3024, MPI_COMM_WORLD, status, ierr ) 
c        print* , 'This is from inside F beginning3 ', mype , 'process'    
  
        call MPI_SEND( pY(NEQL+1), 2, MPI_DOUBLE_PRECISION, myperight, 
     +                   3512, MPI_COMM_WORLD, ierr ) 
  
c        print* , 'This is from inside F beginning4 ', mype , 'process'    
  
        call MPI_RECV(pY(1), 2, MPI_DOUBLE_PRECISION, mypeleft, 
     +                  3512, MPI_COMM_WORLD, status, ierr ) 
  
c        print* , 'This is from inside F beginning5 ', mype , 'process'    
  
       
      do 629 I = 1,(NEQL +4)
         pYT(I) = pY(I)
  629 continue 



      print*, "Y2 is ", Y(2), "from F in  ",mype
      print*, "pY4 is ", pY(4), "from F in  ",mype
      print*, "pYT4 is ", pYT(4), "from F in  ",mype


      DO 1002 I = 1,NEQL

        YDOT(I) = -((((((pY(I+2+1)+pY(I+2))/2)**3)*                   
     +                      ((P(I+1,NEQ,Y)-P(I+0,NEQ,Y))/X))         
     +           -    ((((pY(I+2)+pY(I+2-1))/2)**3)*                  
     +                       ((P(I+0,NEQ,Y)-P(I-1,NEQ,Y))/X))        
     +           )/X) 

 1002 CONTINUE    
#if comments
        print*, 'From F end from ', mype 
#endif        
      RETURN
      END  

       DOUBLE PRECISION  FUNCTION P(I,NEQL,Y)
     
      implicit none 
      include 'mpif.h'
      DOUBLE PRECISION  Y(NEQL),X,sL,pY(NEQL+4)
     + ,pY1,pY2,pY3,pY4
     + ,pYT(NEQ+4) 
      INTEGER I,C,NEQL,mype,ierror 
c      common / PERIODIC_Y / pY1,pY2,pY3,pY4 
      common / PERIODIC_Y / pYT 
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
      INTEGER I, NEQL,mype,ierror 
     + ,IPAR,MU,ML,NRPD 
      
      DIMENSION Y(NEQL), PD(NRPD,NEQL) , PY(NEQL+4)
      common / PERIODIC_Y / pYT
      
      X = L/(NEQ-1)
      BL = 0.02

      call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierror)
     
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

c        print*, 'From J end '
      RETURN
      END
      
C-----------------------------------------------------------------------

