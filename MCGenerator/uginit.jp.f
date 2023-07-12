C************************************************************************
C
      SUBROUTINE UGINIT(fname)
C
C       初期設定
C************************************************************************

      integer IEnergy
      character*100  fname

      INCLUDE 'ugeant.h'
      PARAMETER(ILIN=10,ILOUT=6)
      PARAMETER(IWKTYP=1)
C--------> GEANT の初期化
C
      CALL GINIT
C
C--------> GEANT のバージョンを表示
C 
C      WRITE(ILOUT,100)
C 100  FORMAT(/,3X,'GEANT VERSION 3_21'/)
C
C--------> イベント番号の初期化
C
      NEVENTS = 0
      MAXNST=100000000
C
C--------> data record から読み込む key の設定
C
      CALL UFFKEY
C
C--------> data records を読む
C
      OPEN(UNIT=ILIN,FILE=fname,STATUS='OLD')
      CALL FFSET('LINP',ILIN)
C
      CALL GFFGO
C
      CLOSE(ILIN)
      CALL FFSET('LOUT',ILOUT)
C
C--------> data structure の初期化
C
      CALL GZINIT
C
C-----> グラフィックデバイスを開く
C
      IF(LPLOT(1).NE.0)THEN
         CALL HPLINT(IWKTYP)
         OPEN(3,FILE='out.ps',FORM='formatted',STATUS='unknown')
         CALL HPLCAP(3)
      ENDIF
C
C--------> グラフィックデバイスの初期化

      IF(LPLOT(1).NE.0) CALL GDINIT

C--------> particle & standard material data の読み込み
C   
      CALL GPART
c      PAUSE
      CALL GMATE
c      PAUSE
C
C--------> material と geometry の設定
C
      CALL UGEOM
      IF(LPLOT(1).NE.0) CALL VIEWYZ(1)
C
C--------> cross-section table を compute する。
C
      CALL GPHYSI
C
C--------> banks の表示
C   
      CALL GPRINT('MATE',0)
c      PAUSE
      CALL GPRINT('TMED',0)
c      PAUSE
      CALL GPRINT('VOLU',0)
c      PAUSE
C
C--------> HBOOK の設定
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
