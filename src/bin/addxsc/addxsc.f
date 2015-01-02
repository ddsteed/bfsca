C * APR 06 2008 - RDS - ADD INTEGRAL/MOMENTUM CROSS SECTIONS FROM 0-> 0 TO VMAX
C *-
C * PROGRAM ADDXSC
C * 
C *-/
      double precision energy, xsc, msc, txsc, mxsc
      integer i, j, numv, nume, ifile, tfile, mfile
      character titlex(80), titlem(80)

      read(5,*) numv, nume, ifile

      do i = 1, numv 
         tfile = ifile + i - 1
         mfile = tfile + numv + 1

         read(tfile,200) titlex
         read(tfile,*)
         read(mfile,200) titlem
         read(mfile,*)
      enddo

      write(6,200)    titlex
      write(6,200)    
      write(7,200)    titlem
      write(7,200)    
      
      do j = 1, nume
         txsc = 0.0
         mxsc = 0.0
         do i = 1, numv 
            tfile = ifile + i - 1
            mfile = tfile + numv + 1
            read(tfile,*) energy, xsc
            read(mfile,*) energy, msc
            txsc = txsc + xsc
            mxsc = mxsc + msc
         enddo
         write(6,100) energy, txsc
         write(7,100) energy, mxsc
      enddo

 200  format(80a1)
 100  format(1x, f10.5, e20.10)

      end
