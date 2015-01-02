C * SEP 01 2006 - RDS - CALCULATE THE COORDINATES OF SCATTERED ELECTRON. (XIE,YIE,ZIE)
C *                   - CALCULATE THE COEFFICIENTS OF LEGENDRE POLYNOMIAL OF POLARIZATION
C *                   - POTENTIAL USING GAUSSIAN QUADRATURE.
C *
C * v_lam = \int_0^{\pi} V(rsina,0,rcosa)*P_lam(cosa)*sina*da
C * 

      program main
      implicit none

      integer nmax
      parameter (nmax=100)

      integer M
      integer NN, num, lambda, i
      real*8  Pi, re, xie, yie, zie, theta
      real*8  xg(nmax), wg(nmax)
      real*8  V1(nmax), V2(nmax)
      real*8  P1, P2
      real*8  fx, ss, vL, aL, factor
      real*8  PML
      external PML

C=========================================================================
      read(5,*) M
      read(5,*) NN,num,lambda,re

C-------------------------------------------------
      PI = 3.1415926535897932384626d0
      call gauleg(0,PI/2.0,xg,wg,NN)

      if (M .eq. 1) then
         theta = xg(num)
C * remember that we convert theta = alpha - pi/2
         xie = re*dcos(theta)
         yie = 0.0d0
         zie = re*dsin(theta)
         write(6,*)   xie,yie,zie,-zie
         write(7,500) xie,yie, zie
         write(8,500) xie,yie,-zie

      else if(M .eq. 2) then
         do i =  1, NN
            read(56,*) re, V1(i)
            read(56,*) re, V2(i)
         enddo
         ss = 0.0d0
         do i = 1, NN
            theta = xg(i)
            fx = wg(i)*dcos(theta)
            P1 = PML(lambda,0,-theta+PI/2.0)
            P2 = PML(lambda,0, theta+PI/2.0)
            fx = fx*(V1(i)*P1+V2(i)*P2)
            ss = ss+fx
         enddo

C * remember that PML is NOT normalized function
         factor = 0.5d0*(lambda*2.0+1.0d0)
         vL = ss*factor
         aL = -2.0d0*vL*re**4
         write(6,510) re,vL
         write(7,510) re,aL
      endif

 500  format(1x,3(f20.15))
 510  format(1x,f10.6,e34.26)
      END


