C**********************************************************************
C
      SUBROUTINE UGLAST(fname)
C
C       final process
C**********************************************************************
 
      integer IEnergy
      character*100  fname

      INCLUDE 'ugeant.h'
C-------->
C
      CALL GSUM
C
C-------->
C
      CALL GSUM2
C
C--------> save histogram
C
      CALL HRPUT(0,fname,'N')
C
C--------> final process of GEANT
C
      CALL GLAST
C
C---------> close graphic device
C
      CALL IGEND
C
      RETURN
      END
