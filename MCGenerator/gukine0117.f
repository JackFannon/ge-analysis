C************************************************************************
C
      SUBROUTINE GUKINE
C     (initialization of incident particle)
C************************************************************************     
      INCLUDE 'ugeant.h'
      DIMENSION PLAB(3),VERTEX(3),VERRN(4)
      REAL UPLAB
      COMMON UPLAB
C
C---------> set vertex
C

  100 call grndm(VERRN,3)
      
C Lluis     XTRG = 1.2*sqrt(VERRN(1))*cos(2.*3.1415927*VERRN(2))
C Lluis     YTRG = 1.2*sqrt(VERRN(1))*sin(2.*3.1415927*VERRN(2))
      XTRG = 0.0     ! Lluis needle-like beam
      YTRG = 0.0    ! Lluis needle-like beam
c      if (YTRG.ge.0.5) go to 100
      VERTEX(1) = XTRG
      VERTEX(2) = YTRG
      VERTEX(3) = 15.
C

      CALL GSVERT(VERTEX,0,0,0,0,NVERT)
C
C---------> set particle type & momentum (0.15% fluctuation)
C

      PLAB(1) = 0.
      PLAB(2) = 0.
      PLAB(3) = -1*PKINE(1)
      UPLAB = PLAB(3)
C

      CALL GSKINE(PLAB,IKINE,NVERT,0,0,NT)
C
C     This part adds a second electron going into the Ge detector (by Lluis)
C---------> set vertex
C

c      XTRG = 0.1     ! Lluis needle-like beam
c      YTRG = 0.1     ! Lluis needle-like beam
c      VERTEX(1) = XTRG
c      VERTEX(2) = YTRG
c      VERTEX(3) = 15.
C

c      CALL GSVERT(VERTEX,0,0,0,0,NVERT)

C
C---------> set particle type & momentum (0.15% fluctuation)
C

      PLAB(1) = 0.
      PLAB(2) = 0.
      PLAB(3) = -1*PKINE(1)
      UPLAB = PLAB(3)
C

      CALL GSKINE(PLAB,IKINE,NVERT,0,0,NT)
C

      RETURN
      END

