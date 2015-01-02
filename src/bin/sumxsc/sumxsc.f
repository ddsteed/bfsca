C * APR 05 2008 - RDS - SUM UP THE CROSS SECTIONS OF DIFFERENT SYMMETRIES
C * 
C *-/
      parameter (maxsym = 10, maxexi = 30)
      double precision xsc(maxexi,maxsym), energy
      double precision totxsc(maxexi), tot
      integer numsym, num
      integer i,j
      
      read(5,*) num, numsym

      do j = 1, numsym
         do i = 1, num
            read(5,*) 
            read(5,*) energy, xsc(i,j)
         enddo
      enddo

      tot = 0.0
      do i = 1, num
         totxsc(i) = 0.0
         do j = 1, numsym
            totxsc(i) = totxsc(i) + xsc(i,j)
            tot = tot + xsc(i,j)
         enddo
         write(i+10, 100) energy, totxsc(i)
      enddo
      write(10,100) energy,tot

 100  format(1x,f10.5,d20.10)

      end
      
