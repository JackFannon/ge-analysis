C*************************************************************************
C
      SUBROUTINE VIEWYZ(IVIEW)
C
C       graphic display
C*************************************************************************
      INCLUDE 'ugeant.h'
C
C-----> tree display of tracking medium configuration
C
      CALL GDTREE('VACB',0,000001)
c      PAUSE
      CALL ICLRWK(0,0)
C
C-----> display of overview
C
      CALL GDSPEC('VACB')
c      PAUSE
      CALL ICLRWK(0,0)
C
C-----> display of cross section
C
C      CALL GDRAWC('VACB',1,0.,-30.,10.,.01,.01)
      CALL GDRAWC('VACB',1,0.,10.,10.,ZOOM(1),ZOOM(2))
c      PAUSE
C
      RETURN
      END
