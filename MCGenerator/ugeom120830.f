C**********************************************************************
C
      SUBROUTINE UGEOM
C
C       define geometry
C
C     ! 23.Sep 1998
C     (for the endcap used then)
C
C     1998 11.10 based on ugeom0117.f
C                Simulation of the measurement of the energy loss
C                at various materials which was done on 23 Sep 1998
C
C**********************************************************************
      INCLUDE 'ugeant.h'
      DIMENSION PAR(3),ASC(2),ZSC(2),WSC(2)
      DIMENSION AMYLAR(3),ZMYLAR(3),WMYLAR(3)
      DIMENSION AELECT(3),ZELECT(3),WELECT(3)
      real length,diam,inact,hole_d,hole_l
      real bethick,althick,tithick,scthick
      real gepos1,gepos2,igepos,bepos,ecpos
      real alpos1,alpos2,tipos,scpos,pippos
      double precision gelen1,gelen2
      real skcutg,skcute

      real rpip1i,rpip1o,pip1l,rpip2i,rpip2o,pip2l
      real rfra,twin,tfra,tmidfra,scinz
      real pip1pos,pip2pos,mfrapos
      real brass_d
C
C--------> define material
C  A:  mass number
C  Z:  atomic number
C  D:  density [g cm^-3]
C  RL: radiation length [cm]
C  AL: absorption length (not used) [cm]
C  W:  the atomic number ratio of each element in a compound
C
C--------->  Ge
      DATA AGE/72.59/
      DATA ZGE/32./
      DATA DGE/5.323/
      DATA RLGE/2.30/
C--------->  Ti
      DATA ATI/47.88/
      DATA ZTI/22./
      DATA DTI/4.54/
      DATA RLTI/3.56/
C---------> Scintilator (Plastic)
C
      DATA DSC/1.032/
      DATA RLSC/42.4/      
C--------->  Electrode
C is made from brass
      DATA DELECT/8.41/
      DATA brass_d/0.635/
C--------> Mylar 
      DATA DMYLAR/1.38/
C--------> geometory of Ge crystal
      DATA length/5.96/
      DATA diam/5.73/
C inactive thickness written in the spec. sheet
      DATA inact/0.00003/
C inactive thickness tuned to reproduce Bi207 measurement
c      DATA inact/0.0041/
      DATA hole_d/1.08/
      DATA hole_l/5.19/
C--------> geometory of Be window
      DATA bethick/0.05/
C--------> geometory of Al hoil
      DATA althick/0.08/
C--------> geometory of scintilator
      DATA scthick/0.1/
C--------> geometry of detector capsul
      DATA rcapsul/7.0/
      DATA tcapsul/0.12/
      DATA lcapsul/13.4/
C--------> geometry of Beam pipe
      DATA rpipe_o/3.3/
      DATA rpipe_i/3./
      DATA lpipe/20./
C--------> geometry of Ti window
      DATA tithick/0.003/
      DATA tithick_test/0.01/
C--------> geometry of plastic part of the aluminized mylar 
      DATA plmythick/0.0030/
C--------> geometry of Aluminum part of the aluminized mylar 
      DATA almythick/0.0030/

C--------> position of materials
      gelen1 = length - hole_l - inact
      gelen2 = hole_l
      gepos1 = gelen1/2. + gelen2/2.
      gepos2 = 0.
      igepos = inact/2. + 0.001 + gelen1 + gelen2/2.
      bepos =  bethick/2. + 0.3 + inact + gelen1 + gelen2/2.
      ecpos =  -1.*lcapsul/2.+ bethick + 0.3 + inact + gelen1 + gelen2/2.
      tipos =  bepos + bethick/2. + 0.7 + tithick/2.
      bppos =  tipos + tithick/2. + lpipe/2.

      tipos_test = bepos + bethick/2. + 0.2 + tithick_test/2.
      scpos = tipos_test + tithick_test/2. + scthick/2.
      alpos1 = scpos + scthick/2. + althick/2.
      alpos2 = alpos1 + althick/2. + althick/2.
c  Mylar is 1mm above the Ge crystal
      plmypos = plmythick/2.+ 0.1 + inact + 0.001+gelen1+gelen2/2.
c put 10 micro m gap between Al and Mylar for safety
      almypos = plmypos + plmythick/2. + 0.001
C--------> Definition of scintilator
      ASC(1)= 12.01
      ASC(2)= 1.01
      ZSC(1)= 6.
      ZSC(2)= 1.
      WSC(1)= 0.9224
      WSC(2)= 0.0776
C--------> Definition of Mylar
C  C:H:O =  10:8:4 = weight ratio = 10x12:8x1:4x16 = 0.63:0.04:0.33 


      AMYLAR(1)= 12.01
      AMYLAR(2)= 1.01
      AMYLAR(3)= 16.00
      ZMYLAR(1)= 6.
      ZMYLAR(2)= 1.
      ZMYLAR(3)= 8.
      WMYLAR(1)= 0.63
      WMYLAR(2)= 0.04
      WMYLAR(3)= 0.33

C--------> Definition of brass for electrode
C  60% Cu, 39.25% Zn, 0.75% Sn

      AELECT(1)= 63.546
      AELECT(2)= 65.38
      AELECT(3)= 118.71
      ZELECT(1)= 29.
      ZELECT(2)= 30.
      ZELECT(3)= 50.
      WELECT(1)= 0.60
      WELECT(2)= 0.3925
      WELECT(3)= 0.0075

C
C---------> define material
C
      CALL GSMATE(21,'GE$',AGE,ZGE,DGE,RLGE,0.,0.,0.)
      CALL GSMATE(22,'TI$',ATI,ZTI,DTI,RLTI,0.,0.,0.)
      CALL GSMATE(23,'GI$',AGE,ZGE,DGE,RLGE,0.,0.,0.)

      CALL GSMIXT(24,'SC$',ASC,ZSC,DSC,2,WSC)
      CALL GSMIXT(25,'EL$',AELECT,ZELECT,DELECT,3,WELECT)
      CALL GSMIXT(26,'MY$',AMYLAR,ZMYLAR,DMYLAR,3,WMYLAR)

C
C---------> tracking medium parameter ¤ÎÀßÄê
C
      ISVOL = 0
      IFIELD = 0
      FIELDM = 0.
      EPSIL = 0.0001

      TMAXMD = 0.
      STEMAX = 1.0 
      DEEMAX = 0.01
      STMIN = 0.1

      CALL GPTMED(0)
C     
      CALL GSTMED(21,'AIR$',15,ISVOL,IFIELD,FIELDM,TMAXMD,
     +           STEMAX,DEEMAX,EPSIL,STMIN,0,0)
      CALL GSTMED(22,'VAC$',16,ISVOL,IFIELD,FIELDM,TMAXMD,
     +           STEMAX,DEEMAX,EPSIL,STMIN,0,0)
      CALL GSTMED(23,'EL$',25,ISVOL,IFIELD,FIELDM,TMAXMD,
     +           STEMAX,DEEMAX,EPSIL,STMIN,0,0)

      ISVOL = 1
      CALL GSTMED(24,'GE$',21,1,IFIELD,FIELDM,TMAXMD,
     +           STEMAX,DEEMAX,EPSIL,STMIN,0,0)
      ISVOL = 0
      TMAXMD = 0.
      STEMAX = 0.001 
      DEEMAX = 0.01
      STMIN = 0.000001

      CALL GSTMED(25,'BE$',5,ISVOL,IFIELD,FIELDM,TMAXMD,
     +           STEMAX,DEEMAX,EPSIL,STMIN,0,0)
      CALL GSTMED(26,'Al$',9,ISVOL,IFIELD,FIELDM,TMAXMD,
     +           STEMAX,DEEMAX,EPSIL,STMIN,0,0)
      CALL GSTMED(27,'TI$',22,ISVOL,IFIELD,FIELDM,TMAXMD,
     +           STEMAX,DEEMAX,EPSIL,STMIN,0,0)
      CALL GSTMED(28,'GEIN$',23,ISVOL,IFIELD,FIELDM,TMAXMD,
     +           STEMAX,DEEMAX,EPSIL,STMIN,0,0)
      CALL GSTMED(29,'SCIN$',24,ISVOL,IFIELD,FIELDM,TMAXMD,
     +           STEMAX,DEEMAX,EPSIL,STMIN,0,0)

      DO 200 III=30,33
      CALL GSTMED(III,'IRON$',10,ISVOL,IFIELD,FIELDM,TMAXMD,
     +           STEMAX,DEEMAX,EPSIL,STMIN,0,0)
 200     continue

      CALL GSTMED(34,'LGDE$',24,ISVOL,IFIELD,FIELDM,TMAXMD,
     +           STEMAX,DEEMAX,EPSIL,STMIN,0,0)

      CALL GSTMED(35,'PLMY$',26,ISVOL,IFIELD,FIELDM,TMAXMD,
     +           STEMAX,DEEMAX,EPSIL,STMIN,0,0)
      CALL GSTMED(36,'ALMY$',9,ISVOL,IFIELD,FIELDM,TMAXMD,
     +           STEMAX,DEEMAX,EPSIL,STMIN,0,0)
Cc Ge inactive region (inner surface of the hole) Ge/Li layer
      CALL GSTMED(37,'GELI$',23,ISVOL,IFIELD,FIELDM,TMAXMD,
     +           STEMAX,DEEMAX,EPSIL,STMIN,0,0)

C
C--------> generate volumes
C
C--------> laboratory
      PAR(1) = 30./2.
      PAR(2) = 30./2.
      PAR(3) = 100./2.
      CALL GSVOLU('VACB','BOX ',21,PAR,3,ICALO)

C**************** ENDCAP ********************

C--------> Beampipe (vacuum)
      PAR(1) = 0.
      PAR(2) = rpipe_i/2.
      PAR(3) = lpipe/2.
      CALL GSVOLU('PVAC','TUBE',22,PAR,3,IVOLU)
C--------> Beampipe (iron)
      PAR(1) = rpipe_i/2.
      PAR(2) = rpipe_o/2.
      PAR(3) = lpipe/2.
      CALL GSVOLU('PIPE','TUBE',30,PAR,3,IVOLU)
C--------> Ti-Window
      PAR(1) = 0.
      PAR(2) = rpipe_i/2.
      PAR(3) = tithick/2.
      CALL GSVOLU('TIWN','TUBE',27,PAR,3,IABS)

C--------> Ti-Window(E-loss test sample)
      PAR(1) = 0.
      PAR(2) = 5.0/2.
      PAR(3) = tithick_test/2.
      CALL GSVOLU('TIWT','TUBE',27,PAR,3,IABS)

C--------> Scintilator(E-loss test sample)
      PAR(1) = 0.
      PAR(2) = 5.0/2.
      PAR(3) = scthick/2.
      CALL GSVOLU('SCTI','TUBE',29,PAR,3,IABS)
C--------> Al(E-loss test sample)
      PAR(1) = 0.
      PAR(2) = 5.0/2.
      PAR(3) = althick/2.
      CALL GSVOLU('ALC1','TUBE',26,PAR,3,IABS)
      PAR(1) = 0.
      PAR(2) = 5.0/2.
      PAR(3) = althick/2.
      CALL GSVOLU('ALC2','TUBE',26,PAR,3,IABS)

C****************** Germanium Detector ***********************

C--------> Detecor capsul (vacuum)
      PAR(1) = 0.
      PAR(2) = rcapsul/2.-tcapsul
      PAR(3) = lcapsul/2.
      CALL GSVOLU('ECVA','TUBE',22,PAR,3,ICALO)
C--------> Detecor capsul (Al)
      PAR(1) = rcapsul/2.-tcapsul
      PAR(2) = rcapsul/2.
      PAR(3) = lcapsul/2.
      CALL GSVOLU('ECAL','TUBE',26,PAR,3,ICALO)
C--------> Be-Window in front of Ge crystal
      PAR(1) = 0.
      PAR(2) = rcapsul/2.-tcapsul
      PAR(3) = bethick/2.
      CALL GSVOLU('BEWN','TUBE',25,PAR,3,IABS)
C--------> Ge Crystal
      PAR(1) = 0.
      PAR(2) = diam/2.
      PAR(3) = gelen1/2.
      CALL GSVOLU('GEC1','TUBE',24,PAR,3,IABS)
      PAR(1) = hole_d/2.
      PAR(2) = diam/2.
      PAR(3) = gelen2/2.
      CALL GSVOLU('GEC2','TUBE',24,PAR,3,IABS)
C--------> Inactive region in front of the Ge crystal
      PAR(1) = 0.
      PAR(2) = diam/2.
      PAR(3) = inact/2.
      CALL GSVOLU('GEIN','TUBE',28,PAR,3,IABS)
C--------> Inactive region (inner surface of the hole)
      PAR(1) = hole_d/2. - 0.07
      PAR(2) = hole_d/2.
      PAR(3) = hole_l/2.
      CALL GSVOLU('GELI','TUBE',37,PAR,3,IABS)
C--------> Electorode hole in the Ge crystal 
      PAR(1) =0.
      PAR(2) = brass_d/2.
      PAR(3) = hole_l/2.
      CALL GSVOLU('ELND','TUBE',23,PAR,3,IABS)
C--------> Plastic part of the Aluminized mylar
      PAR(1) = 0.
      PAR(2) = diam/2.
      PAR(3) = plmythick/2.
      CALL GSVOLU('PLMY','TUBE',35,PAR,3,IABS)
C--------> Plastic part of the Aluminized mylar
      PAR(1) = 0.
      PAR(2) = diam/2.
      PAR(3) = almythick/2.
      CALL GSVOLU('ALMY','TUBE',36,PAR,3,IABS)

C---------> configuration of volume
C
C---------> Endcap
      CALL GSPOS('PIPE',1,'VACB',0.,0., bppos,  0,'MANY')
      CALL GSPOS('PVAC',1,'VACB',0.,0., bppos,  0,'MANY')
      CALL GSPOS('TIWN',1,'VACB',0.,0., tipos,  0,'ONLY')

C---------> Test sample
c      CALL GSPOS('SCTI',1,'VACB',0.,0., scpos,  0,'ONLY')
c      CALL GSPOS('TIWT',1,'VACB',0.,0., tipos_test,  0,'ONLY')
c      CALL GSPOS('ALC1',1,'VACB',0.,0., alpos1,  0,'ONLY')
c      CALL GSPOS('ALC2',1,'VACB',0.,0., alpos2,  0,'ONLY')
C---------> Germanium detector
      CALL GSPOS('ECVA',1,'VACB',0.,0., ecpos,  0,'MANY')
      CALL GSPOS('ECAL',1,'VACB',0.,0., ecpos,  0,'MANY')
      CALL GSPOS('BEWN',1,'VACB',0.,0., bepos , 0,'ONLY')
      CALL GSPOS('GEIN',1,'VACB',0.,0., igepos, 0,'ONLY')
      CALL GSPOS('GEC1',1,'VACB',0.,0., gepos1, 0,'ONLY')
      CALL GSPOS('GEC2',1,'VACB',0.,0., gepos2, 0,'ONLY')
      CALL GSPOS('GELI',1,'VACB',0.,0., gepos2, 0,'ONLY')
      CALL GSPOS('ELND',1,'VACB',0.,0., gepos2, 0,'ONLY')
      CALL GSPOS('PLMY',1,'VACB',0.,0., plmypos, 0,'ONLY')
      CALL GSPOS('ALMY',1,'VACB',0.,0., almypos, 0,'ONLY')

C---------> redifine tracking parameters individually.
C
      skcutg = 1.e-4
      skcute = 2.64064e-4
      CALL GSTPAR(26,'CUTGAM',skcutg)
      CALL GSTPAR(26,'CUTELE',skcute)
      CALL GSTPAR(27,'CUTGAM',skcutg)
      CALL GSTPAR(27,'CUTELE',skcute)
      CALL GSTPAR(29,'CUTGAM',skcutg)
      CALL GSTPAR(29,'CUTELE',skcute)
      CALL GSTPAR(34,'CUTGAM',skcutg)
      CALL GSTPAR(34,'CUTELE',skcute)
      CALL GSTPAR(35,'CUTGAM',skcutg)
      CALL GSTPAR(35,'CUTELE',skcute)

c      CALL GSTPAR(25,'STRA',1.)
c      CALL GSTPAR(26,'STRA',1.)
c      CALL GSTPAR(27,'STRA',1.)
c      CALL GSTPAR(28,'STRA',1.)
C
C---------> close geometry banks
C
      CALL GGCLOS
C
      RETURN
      END
