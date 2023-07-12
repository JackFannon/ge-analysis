
C***************************************************************************
C
      SUBROUTINE GSUM
C
C***************************************************************************
      INCLUDE 'ugeant.h'
      PARAMETER (IMAX=180)
      DIMENSION M(IMAX)
      INTEGER MAXM,K,L,IP1,IP2,IP3,NFDN1,NFDN2,NEXR,ITOT
      REAL IEFF
C
      DO 20 I=1,IMAX
         M(I) = HI(200,I)
 20   CONTINUE
C
       MAXM = M(1)
      DO 30 I=2,IMAX
        IF(M(I).GT.MAXM) MAXM = M(I)
 30   CONTINUE
C
      K = MOD(MAXM,10)
      L = (MAXM-K)/10
        IF(K.GE.5) THEN
           L = L + 1
        ENDIF 
C    
      DO 40 I=1,IMAX
        IF(M(I).EQ.MAXM) THEN
            IP1 = I
            IP2 = I
            IP3 = I
          GOTO 50
        ENDIF
 40   CONTINUE
C
       NFDN1=0
 50   DO 60 I=1,IMAX
         IP1 = IP1 + 1
           IF(M(IP1).EQ.0) THEN
               GOTO 70
           ENDIF
         NFDN1 = NFDN1 + M(IP1)
 60   CONTINUE   
C
       NFDN2=0
 70   DO 80 I=1,IMAX
         IP2 = IP2 - 1
           IF(M(IP2).LT.L) THEN
               GOTO 90
           ENDIF
         NFDN2 = NFDN2 + M(IP2)
 80   CONTINUE   
C 
 90   NFDN = NFDN1 + NFDN2 + MAXM
C
       NEXR = 0
      DO 100 I=1,IMAX 
        IF(IP2.EQ.0) THEN
          GOTO 110
        ENDIF
         NEXR = NEXR + M(IP2)
          IP2 = IP2 - 1
 100  CONTINUE    
C
 110  ITOT = NFDN + NEXR
       RITOT = REAL(ITOT)
       RNFDN = REAL(NFDN)
      IEFF = RNFDN/RITOT
C
      WRITE(6,*) IP3,ITOT,'   efficiency=',IEFF
C
      RETURN
      END
