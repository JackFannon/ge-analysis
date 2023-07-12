C************************************************************************
C
      SUBROUTINE UGINIT(fname)
C
C       �������
C************************************************************************

      integer IEnergy
      character*100  fname

      INCLUDE 'ugeant.h'
      PARAMETER(ILIN=10,ILOUT=6)
      PARAMETER(IWKTYP=1)
C--------> GEANT �ν����
C
      CALL GINIT
C
C--------> GEANT �ΥС�������ɽ��
C 
C      WRITE(ILOUT,100)
C 100  FORMAT(/,3X,'GEANT VERSION 3_21'/)
C
C--------> ���٥���ֹ�ν����
C
      NEVENTS = 0
      MAXNST=100000000
C
C--------> data record �����ɤ߹��� key ������
C
      CALL UFFKEY
C
C--------> data records ���ɤ�
C
      OPEN(UNIT=ILIN,FILE=fname,STATUS='OLD')
      CALL FFSET('LINP',ILIN)
C
      CALL GFFGO
C
      CLOSE(ILIN)
      CALL FFSET('LOUT',ILOUT)
C
C--------> data structure �ν����
C
      CALL GZINIT
C
C-----> ����ե��å��ǥХ����򳫤�
C
      IF(LPLOT(1).NE.0)THEN
         CALL HPLINT(IWKTYP)
         OPEN(3,FILE='out.ps',FORM='formatted',STATUS='unknown')
         CALL HPLCAP(3)
      ENDIF
C
C--------> ����ե��å��ǥХ����ν����

      IF(LPLOT(1).NE.0) CALL GDINIT

C--------> particle & standard material data ���ɤ߹���
C   
      CALL GPART
c      PAUSE
      CALL GMATE
c      PAUSE
C
C--------> material �� geometry ������
C
      CALL UGEOM
      IF(LPLOT(1).NE.0) CALL VIEWYZ(1)
C
C--------> cross-section table �� compute ���롣
C
      CALL GPHYSI
C
C--------> banks ��ɽ��
C   
      CALL GPRINT('MATE',0)
c      PAUSE
      CALL GPRINT('TMED',0)
c      PAUSE
      CALL GPRINT('VOLU',0)
c      PAUSE
C
C--------> HBOOK ������
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
