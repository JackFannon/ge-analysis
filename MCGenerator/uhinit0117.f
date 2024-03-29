C**********************************************************************
      SUBROUTINE UHINIT
C       generate HBOOK
C
C     ! Change number of bin of lun 500 ( 30 -> 300 )
C     1997 07 30 change LUN of 'Total Energy Deposit w/o resolution'
C     1998 11 10 Add Histogram( # fo bin = 4096,1024)
C**********************************************************************
      INCLUDE 'ugeant.h'
      real offset,slope,hmin,hmax
      parameter (offset=1.12369)
      parameter (slope=6.74495)

      hmin = offset/1000.
      hmax = (slope*4096.+offset)/1000.

      CALL HBOOK1(11,' Total Energy Deposit(MeV) $'
     &     ,2048,hmin,hmax,0.)
      CALL HBOOK1(21,' Total Energy Deposit(MeV) $'
     &     ,4096,hmin,hmax,0.)
      CALL HBOOK1(31,' Total Energy Deposit(MeV) $'
     &     ,1024,hmin,hmax,0.)
      CALL HBOOK1(41,' Total Energy Deposit(MeV) $'
     &     ,20000,0.0,20.,0.)
      CALL HBOOK1(12,' Total Energy Deposit(MeV) with resolution$'
     &     ,2048,hmin,hmax,0.)
      CALL HBOOK1(120,' Total Energy after Ti window$'
     &     ,900,1.,19.,0.)
      CALL HBOOK1(130,'Momentum after Ti window (MeV/c)$'
     &     ,900,1.,19.,0.)
      CALL HBOOK1(140,' Total Energy after Be window$'
     &     ,900,1.,19.,0.)
      CALL HBOOK1(150,'Momentum after Be window (MeV/c)$'
     &     ,900,1.,19.,0.)
      CALL HBOOK1(200,'Be window + Inactive region + Ge crystal $'
     &     ,2250,0.,18.,0.)
      CALL HBOOK1(300,' Total Energy Deposit in Be (MeV) $'
     &     ,300,0.,3.E-1,0.)
      CALL HBOOK1(400,' Total Energy Deposit in Ti (MeV) $'
     &     ,200,0.,2.E-1,0.)
      CALL HBOOK1(500,' Total Energy Deposit in Inactive region (MeV) $'
     &     ,300,0.,1.E-1,0.)
      CALL HBOOK1(600,' Total Energy Deposit in Scintilator (MeV) $'
     &     ,200,0.,5.e-1,0.)
      CALL HBOOK1(700,' Total Energy Deposit in Al x 2 (MeV) $'
     &     ,1000,0.,1.E-1,0.)
      CALL HBOOK1(800,'Energy Deposit in Trg.counter (MeV) $'
     &     ,500,0.,5.e-1,0.)
      CALL HBOOK1(900,' Energy that escapes from Ge(MeV) $'
     &     ,3600,0.,18.,0.)

      RETURN
      END
