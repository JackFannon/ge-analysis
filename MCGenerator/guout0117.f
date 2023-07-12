C*********************************************************************
C
      SUBROUTINE GUOUT
C
C       GEANT user routine
C       (called at the end of each event)
C       1997 07 30 Change LUN 100 -> 11, 110 -> 12
C       1998 11 10 Add Histogram( # fo bin = 4096,1024)
C*********************************************************************
      INCLUDE 'ugeant.h'
      INTEGER BIN
      REAL rBIN , DEOBS, DELOSS
      REAL DTSIGMA,DTMEAN,GESIGMA,PRESIGMA,AMPSIGMA
      COMMON/ENERGY/EELE
      COMMON/TOTE/DETOT,DETOT0,DETOT1,DETOT2,DETOT3,DETOT4
      COMMON/EXIT/EXTOT,EXMOM,BETOT,BEMOM
      COMMON/ESCEN/ESCA
      COMMON UPLAB
C
C--------> fill histogram
C

      IF(DETOT.GT.0.) THEN
         DETOT = DETOT * 1000
C--------> consider the detector resolution
c         GESIGMA = (( 0.15 * (DETOT*1e6) * 2.96)**0.5) * (1e-6)
c         PRESIGMA = 0.2894 * (1e-3)
c         AMPSIGMA = 2.119 * (1e-3)
c         AMPSIGMA = 0.
c         DTSIGMA = ((GESIGMA**2)+(PRESIGMA**2)+(AMPSIGMA**2))**0.5
c     The detector resolution obtained from Ni real data
c         pni = 8.999
c         a = (4.2011e-3) / sqrt(pni)
         DTSIGMA = 4.098739e-3
         DTMEAN = DETOT
         DEOBS = rngaus(DTSIGMA,DTMEAN)
C--------> Use the same binning as MCA
         BIN = ( DEOBS*1000 - 34.01 ) / 5.224
         rBIN = BIN
C
C--------> Convert GeV to MeV
         DETOT0 = DETOT0 * 1000
         DETOT1 = DETOT1 * 1000
         DETOT2 = DETOT2 * 1000
         DETOT3 = DETOT3 * 1000
         DETOT4 = DETOT4 * 1000
         EXTOT = EXTOT * 1000
         EXMOM = EXMOM * 1000
         BETOT = BETOT * 1000
         BEMOM = BEMOM * 1000

         DELOSS = DETOT1 +  DETOT3 + DETOT4
         SKOBS = DETOT0 + DETOT2 + DETOT
 
         ESCA = ESCA * 1000
         write(40) DETOT,EXTOT,BETOT,DETOT0,DETOT1,DETOT2,
     &        DETOT3,DETOT4,ESCA
         CALL HF1(11,DETOT,1.)   !Total Energy Deposit in Ge (2048)
         CALL HF1(21,DETOT,1.)   !Total Energy Deposit in Ge (4096)
         CALL HF1(31,DETOT,1.)   !Total Energy Deposit in Ge (1024)
         CALL HF1(41,DETOT,1.)   !Total Energy Deposit in Ge(1keV/bin)
         CALL HF1(12,DEOBS,1.)   !Total Energy Deposit in Ge(resol.)
         CALL HF1(120,EXTOT,1.)  !Total Energy after Ti
         CALL HF1(130,EXMOM,1.)  !Momentum after Ti
         CALL HF1(140,BETOT,1.)  !Total Energy after Be
         CALL HF1(150,BEMOM,1.)  !Momentum after Be
         CALL HF1(200,SKOBS,1.)  !E-loss in Be + Inactive. + Ge
         CALL HF1(300,DETOT0,1.) !E-loss in Be
         CALL HF1(400,DETOT1,1.) !E-loss in Ti
         CALL HF1(500,DETOT2,1.) !E-loss in Inactive region
         CALL HF1(600,DETOT3,1.) !E-loss in Scinti
         CALL HF1(700,DETOT4,1.) !E-loss in Alx2
         CALL HF1(800,DELOSS,1.) !E-loss in Trg.counter
         IF(ESCA.GT.0) THEN
            CALL HF1(900,ESCA,1.) !Escape energy from Ge
         ENDIF
      END IF
C
      RETURN
      END
