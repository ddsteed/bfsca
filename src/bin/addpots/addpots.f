C * AUG 20 2005 - RDS - GENERALIZE TO HETERONUCLEAR DIATOMIC MOLECULE AND INCLUDE DSG POTENTIAL
C *
C * JAN 01 1984 - WKT - ADD STATIC POTENTIAL AND POLARIZATION POTENTIAL
C *-
C * PROGRAM ADDPOTS
C *
C *-/
      double precision stat,r,pol
      character title(80)
      integer i
      integer npt, ilam, lambda
      integer stlam, pollam, lammin, lamdeg

      PARAMETER (NREG=10)
      real*8 rstart(NREG), rstep(NREG), rstop(NREG)
      real*8 za, ra, zb, rb
      integer ncent

 499  format(80a1)
 500  format(i5)
 501  format(3f10.5)
 502  format(2e15.6)
 503  format(f10.5,d24.16)
C503  format(f10.5,d23.16)
 504  format(2i5)

      read (10,499) title
      read (11,*)
      write (12,499) title

      read (10,500) npt
      read (11,500) npt
      write (12,500)npt
      do i = 1, npt
         read (10,501) rstart(i), rstep(i), rstop(i)
         read (11,501) rstart(i), rstep(i), rstop(i)
         write (12,501)rstart(i), rstep(i), rstop(i)
      end do
      read (10,500) ncent
      read (11,500) ncent
      read (10,502) za, ra
      read (10,502) zb, rb
      read (11,502) za, ra
      read (11,502) zb, rb

      write (12,500) ncent
      write (12,502) za, ra
      write (12,502) zb, rb

      read (5,*) stlam, pollam, lamdeg
      if (stlam .le. pollam) then
         lammin = stlam 
      else
         lammin = pollam 
      endif

      write (6,*) " Lambda of static potential = ", stlam
      write (6,*) " Lambda of polar. potential = ", pollam

      do ilam = 0, lammin, lamdeg
         read (10,504) lambda, npt
         read (11,504) lambda, npt
         write (12,504)lambda, npt
         do i = 1, npt
           read (10,503) r, stat
           read (11,503) r, pol
           write (12,503)r, stat+pol
        end do
      end do

      if (lammin .lt. stlam) then
         do ilam = lammin+lamdeg, stlam, lamdeg
            read (10,504,end=999)lambda, npt
            write (12,504)lambda, npt
            do i=1,npt
               read (10,503) r, stat
               write (12,503)r, stat
            end do
         end do
      endif
 999  continue
      stop
      end
