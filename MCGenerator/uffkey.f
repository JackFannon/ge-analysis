C*************************************************************************
C
      SUBROUTINE UFFKEY
C       set the key of data record
C*************************************************************************

      INCLUDE 'ugeant.h'
C
      CALL FFKEY('SMAX',MAXNST,1,'INTEGER')
      CALL FFKEY('ZOOM',ZOOM,2,'REAL')
      RETURN
      END
