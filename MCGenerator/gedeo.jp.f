C********************************************************************
C
      PROGRAM NAI
C
C       Ge detecor simulator
C
C********************************************************************
      PARAMETER(NG=500000,NH=1000000)
      COMMON/PAWC/H(NH)
      COMMON/GCBANK/Q(NG)
      PARAMETER(TIME=1000000.0)

      character*100 argv, fname_card, fname_hbk

c*** check arguments
      narg = iargc()
      if (narg .lt. 2) then
         print *, 'Usage : gedeo [card name] f_hbk'
         call exit(1)
      endif

      call getarg(1, fname_card)
      call getarg(2, fname_hbk)
C
C-----> GEANT と HBOOK のメモリの初期化
C
      CALL TIMEST(TIME)
      CALL GZEBRA(NG)
      CALL HLIMIT(-NH)
C
C-----> GEANT の初期化
C
      CALL UGINIT(fname_card)
C
C-----> イベントを生成
C
      CALL GRUN
C
C-----> 終了処理
C
      CALL UGLAST(fname_hbk)
C
      STOP 
      END
