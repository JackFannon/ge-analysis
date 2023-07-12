************************************************************************
*     ---------------------
      FUNCTION RNGAUS(DTSIGMA,DTMEAN)
*     ---------------------
*
*     (Purpose)
*       MAKE GAUSSIAN DISTRIBUTED RANDUM NUMBER
*
*     (Input)
*       DTMEAN : MEAN
*       DTSIGMA : SIGMA      
*
*     (Output)
*       RNGAUS : GAUSSIAN 
*      
*     (Creation Date and Author)
*       1997.01.24 ; N.Sakurai    Programed for SCINTI
*     (Comment)
*     
************************************************************************
      INCLUDE 'ugeant.h'
      REAL rtest,factor,RANDOM1,RANDOM2,R,SEL
      DIMENSION RN(3)
c      DATA PI/3.14159265318979/
      CALL GRNDM(RN,3)
      factor = ( -2 * log(RN(1)))**0.5
      RANDOM1 = factor * cos(2 * PI * RN(2))
      RANDOM2 = factor * sin(2 * PI * RN(2))
      IF(SEL.EQ.0) THEN
         R = RANDOM1
         SEL = 1
      ELSEIF(SEL.EQ.1) THEN
         R = RANDOM2
         SEL = 0
      ENDIF
      RNGAUS=R*DTSIGMA+DTMEAN
      RETURN
      END

