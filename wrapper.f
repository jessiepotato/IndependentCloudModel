      program wrapper
C      use interp
      use fluffycloud
      use clouddata
      implicit none
!
      teplota = (/ 113.5, 112.9, 111.4, 110.9, 111.7, 112.9, 114.1,
     & 115.7, 117.4, 119.5, 121.2, 123.8, 127.2, 130.1, 132.5, 136.9,
     & 142.4, 146.9, 151.2, 157.5, 165. /)
      telak = (/ 0.1, 0.1122, 0.12589, 0.14125, 0.15849, 0.17783,
     & 0.19953, 0.22387, 0.25119, 0.28184, 0.31623, 0.35481, 0.39811,
     & 0.44668, 0.50119, 0.56234, 0.63096, 0.70795, 0.79433, 0.89125,
     & 1. /)
!
      open (unit=20,file='tuesclouds.txt',action="read",status="old")
      do i=1,128
      read (20, '(F5.1, F10.3)' ) INTEMP(i), INPRES(i)
      end do
      close(UNIT=20)
C
C      OPEN (1, file = 'tuesclouds.txt')
C      DO
C          READ (1,*, END=10)
C          nlines = nlines + 1
C      END DO
C   10 CLOSE (1)
C      print *, 'lines in file: ', nlines
C
C      teplota = INTEMP
C      telak = INPRES*1.e-6
C
      call FLUFFY(stredni,cislo)
      open (unit=20,file="outfluff.txt",action="write",status="replace")
C      do i=1, size(telak)! layers-1
C        write (20,*) i, telak(i), teplota(i)
C      enddo
!
C      write (20,*)'i','Pz(i)','Tz(i)','Qc(i)','Qv(i)','Qt(i)','Qs(i)'
      do i=1, size(zed)
C        write (20,*) i, Pz(i), Tz(i), stredni(i)
        write (20,*) i, Pz(i),Tz(i),Qc(i),Qv(i),Qt(i),Qs(i)
C        write (20,*) i, zed(i)
      enddo
      close (20)
      print *, 'input size: ', size(teplota)
      print *, 'calculation size: ', size(Tz)
!
      end program wrapper
