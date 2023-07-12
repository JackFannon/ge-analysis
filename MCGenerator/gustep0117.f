C**********************************************************************
C
      SUBROUTINE GUSTEP
C
C       GEANT user routine
C       (called at the end of each step)
C
C     1997 06.06 Modify 'NUMED' according to ugeom0117.f
C     1997 06.09 Output the electron momentum and total energy after Ti window
C     1998 06.01 Output the electron momentum and total energy after Be window
C**********************************************************************
      INCLUDE 'ugeant.h'
      REAL DETOT
      COMMON/ENERGY/EELE
      COMMON/TOTE/DETOT,DETOT0,DETOT1,DETOT2,DETOT3,DETOT4
      COMMON/EXIT/EXTOT,EXMOM,BETOT,BEMOM
      COMMON/ESCEN/ESCA
      COMMON Geflag
C
C----------> Store secondary info in GEANT memory structure
C
      IF(NGKINE.GT.0) THEN 
         CALL GSKING(0)
      ENDIF

C----------> Energy that escapes from Ge crystal
      IF(NUMED.EQ.21) THEN
         IF(Geflag.EQ.1) THEN
            IF(VECT(7).GT.0) THEN
               ESCA = ESCA + GETOT
            ENDIF
            ISTOP = 1.
c            Geflag = 0.
            GOTO 1000
         ENDIF
      ENDIF
C
C----------> Calculate the energy loss inside Ge
C
      IF(NUMED.EQ.24) THEN
         IF (DESTEP.GT.0.) THEN
            DETOT = DETOT + DESTEP
C----------> Geflag: Flag that a particle enters Ge
            Geflag = 1
         ENDIF
      ENDIF
C
C----------> Calculate the energy loss inside Be
C
      IF(NUMED.EQ.25) THEN
         IF(DESTEP.GT.0) THEN
              DETOT0 = DETOT0 + DESTEP
         ENDIF
      ENDIF
C----------> Calculate the energy loss inside Ti
C
      IF(NUMED.EQ.27) THEN
         IF(DESTEP.GT.0) THEN
              DETOT1 = DETOT1 + DESTEP
         ENDIF
      ENDIF
C----------> Calculate the energy loss inside Ge inactive region
C
      IF(NUMED.EQ.28) THEN
         IF(DESTEP.GE.0) THEN
            DETOT2 = DETOT2 + DESTEP
         ENDIF
      ENDIF
C----------> Calculate the energy loss inside scintillator
C
      IF(NUMED.EQ.29) THEN
         IF(DESTEP.GE.0) THEN
            DETOT3 = DETOT3 + DESTEP
         ENDIF
      ENDIF

C----------> Calculate the energy loss inside Al foil
C
      IF(NUMED.EQ.26) THEN
         IF(DESTEP.GE.0) THEN
            DETOT4 = DETOT4 + DESTEP
         ENDIF
      ENDIF

C----------> electron momentum and total energy after Ti window
C
      if(IPART.EQ.3.and.INWVOL.EQ.2.and.NUMED.EQ.27.
     +   and.VECT(6).le.0.)then
         EXTOT  = GETOT
         EXMOM  = VECT(7)
c         write(*,*) 'exit mom=',VECT(7)*1000,'  total=',GETOT
      endif
C----------> electron momentum and total energy after Be window
C
      if(IPART.EQ.3.and.INWVOL.EQ.2.and.NUMED.EQ.25.
     +   and.VECT(6).le.0.)then
         BETOT  = GETOT
         BEMOM  = VECT(7)
      endif
C-----------> draw track
C
      IF(LPLOT(1).NE.0)CALL GDCXYZ
C-----------> debug event
C
 1000 CONTINUE
      CALL GDEBUG
C
      END
