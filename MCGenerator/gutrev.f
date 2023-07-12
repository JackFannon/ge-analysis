C**********************************************************************
      SUBROUTINE GUTREV
C
C       GEANT3 user routine
C       (control the tracking of each event)
C**********************************************************************
      INCLUDE 'ugeant.h'
      INTEGER Geflag
      COMMON/ENERGY/EELE
      COMMON/TOTE/DETOT,DETOT0,DETOT1,DETOT2,DETOT3,DETOT4
      COMMON/EXIT/EXTOT,EXMOM,BETOT,BEMOM
      COMMON/ESCEN/ESCA
      COMMON Geflag
      REAL DETOT,ESCA,EXTOT,EXMOM,BETOT,BEMOM
C
C--------> initialization before each event
C
      EELE = 0.
      DETOT = 0.
      DETOT0 = 0.
      DETOT1 = 0.
      DETOT2 = 0.
      DETOT3 = 0.
      DETOT4 = 0.

      ESCA = 0.

      Geflag = 0
      EXTOT = 0.
      EXMOM = 0.
      BETOT = 0.
      BEMOM = 0.
C
C--------> start tracking
C
      CALL GTREVE
C
      RETURN
      END
