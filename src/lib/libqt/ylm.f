      double precision function ylm(l,m,th)
c - spherical harmonics without the Exp(i phi) factor
      implicit none
      integer l,m
      double precision th
      double precision plgndr,factrl

      integer mm, phase
      double precision fourpi
      data fourpi / 12.566370614359172954d0/

      if(m.lt.0) then
	mm = -m
	phase = (-1)**mm
      else
	mm = m
	phase = 1
      endif

      ylm = phase * 
     &        sqrt( ((2*l+1)/fourpi)*(factrl(l-mm)/factrl(l+mm))) *
     &        plgndr(l,mm,cos(th))
      
      return
      end

      double precision function factrl(n)
      implicit none
      integer n
CU    USES gammln
      integer j,ntop
      double precision a(33),gammln
      save ntop,a
      data ntop,a(1)/0,1./
      if (n.lt.0) then
        stop 'negative factorial in factrl'
      else if (n.le.ntop) then
        factrl=a(n+1)
      else if (n.le.32) then
        do 11 j=ntop+1,n
          a(j+1)=j*a(j)
11      continue
        ntop=n
        factrl=a(n+1)
      else
        factrl=exp(gammln(n+1.d0))
      endif
      return
      end

      double precision function gammln(xx)
      implicit none
      double precision xx
      integer j
      double precision ser,stp,tmp,x,y,cof(6)
      save cof,stp
      data cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      end

      double precision function plgndr(l,m,x)
      implicit none
      integer l,m
      double precision x
      integer i,ll
      double precision fact,pll,pmm,pmmp1,somx2
      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.)then
        write(6,*)'bad arguments in plgndr'
        write(6,*)'l,m,x = '
        write(6,*)l,m,x
	stop
      endif
      pmm=1.d0
      if(m.gt.0) then
        somx2=sqrt((1.d0-x)*(1.d0+x))
        fact=1.d0
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.d0
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
          do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
      return
      end
