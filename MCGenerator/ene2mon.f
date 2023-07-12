      implicit none

c     energy -> momentum

      integer ene, narg
      integer enes, enee 
      real*8 mom
      character*200 fname, cenes, cenee

c*** check arguments
      narg = iargc()
      if (narg .lt. 2) then
         print *, 'Usage : ene2mon [energy start] [energy end]'
         call exit(1)
      endif

c*** prepare input files
      call getarg(1, cenes)
      call getarg(2, cenee)
           
      read(cenes,*) enes
      read(cenee,*) enee

*     13640 13679
*     8889 8850
*     5000 5999

c         ene = 4670
c         do while(ene.le.4700)
c         ene = 5080
c         do while(ene.le.5119)
c         ene = 8850
c         do while(ene.le.8899)
c         ene = 13640
c         do while(ene.le.13699)
c         ene = 18900
c         do while(ene.le.18920)
      
      ene = enes
      do while(ene.le.enee)
         write(fname,100) ene
 100     format('./ene/',i5.5)
         
         open(10,file=fname)
         mom = sqrt((real(ene)*1.E-3)**2 - 0.511**2)
         write(10,"(A,F15.10)") 'KINE	3', mom*1.E-3
         close(10)

         ene = ene + 1
         enddo

      stop
      end
