C**************************************************************************
C
      SUBROUTINE GSUM2
C
C**************************************************************************
      INCLUDE 'ugeant.h'
      PARAMETER(IMAX2=60)
      DIMENSION M2(IMAX2)
      INTEGER MAXM2,IP
C
      DO 10 I=1,IMAX2
       M2(I)=HI(300,I)
 10   CONTINUE
C
        MAXM2 = M2(1)
      DO 20 I=2,IMAX2
       IF(M2(I).GT.MAXM2) MAXM2 = M2(I)
 20   CONTINUE
C
      DO 30 I=1,IMAX2
       IF(M2(I).EQ.MAXM2) THEN
           IP = I 
         GOTO 40
       ENDIF
 30   CONTINUE 
C 
 40   WRITE(6,*) 'peak ',IP
C
      RETURN
      END
