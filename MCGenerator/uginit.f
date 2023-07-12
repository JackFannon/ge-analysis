C************************************************************************
C
      SUBROUTINE UGINIT(fname)
C
C       Initialization
C************************************************************************

      integer IEnergy
      character*100  fname

      INCLUDE 'ugeant.h'
      PARAMETER(ILIN=10,ILOUT=6)
      PARAMETER(IWKTYP=1)
C--------> GEANT initialization
C
      CALL GINIT
C
C--------> Represent the GEANT version
C 
C      WRITE(ILOUT,100)
C 100  FORMAT(/,3X,'GEANT VERSION 3_21'/)
C
C--------> Initialize event number
C
      NEVENTS = 0
      MAXNST=100000000
C
C--------> set key from data record
C
      CALL UFFKEY
C
C--------> read data records
C
      OPEN(UNIT=ILIN,FILE=fname,STATUS='OLD')
      CALL FFSET('LINP',ILIN)
C
      CALL GFFGO
C
      CLOSE(ILIN)
      CALL FFSET('LOUT',ILOUT)
C
C--------> initialize data structure
C
      CALL GZINIT
C
C-----> open graphic device
C
      IF(LPLOT(1).NE.0)THEN
         CALL HPLINT(IWKTYP)
         OPEN(3,FILE='out.ps',FORM='formatted',STATUS='unknown')
         CALL HPLCAP(3)
      ENDIF
C
C--------> initialize graphic device

      IF(LPLOT(1).NE.0) CALL GDINIT

C--------> readout particle & standard material data
C   
      CALL GPART
c      PAUSE
      CALL GMATE
c      PAUSE
C
C--------> set material & geometry
C
      CALL UGEOM
      IF(LPLOT(1).NE.0) CALL VIEWYZ(1)
C
C--------> compute cross-section table
C
      CALL GPHYSI
C
C--------> banks
C   
      CALL GPRINT('MATE',0)
c      PAUSE
      CALL GPRINT('TMED',0)
c      PAUSE
      CALL GPRINT('VOLU',0)
c      PAUSE
C
C--------> set HBOOK
C
      CALL UHINIT
C
C
C--------> Initialises the random number generator
c      CALL GRNDMQ(0.,0.,1,' ')
C
C--------> Open output file
      open(40,file='geout.dat',form='unformatted',status='unknown')
      RETURN
      END