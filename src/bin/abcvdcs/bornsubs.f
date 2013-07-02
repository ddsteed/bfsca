c
c  This is a set of subroutines taken from the r-closure code of
c  the lavib package.  Almost reusable code!  Some modification
c  was done for specific usage here.  The routines here are:
c
c  bornkel     - calculates each Born K matrix element
c  rlam        - returns the (analytic) radial integral
c  facx        - factorial function
c  gam         - gamma function for half integer values only
c  hyper       - hypergeometric function 2F1
c
c  we have also written the function called glam and appended it
c
      double precision function bornkel(l,lp,k,kp,alpha0,alpha2,q
     $     ,lamda)
c-----------------------------------------------------------------
c  we return the Born K matrix element analytic result
c-----------------------------------------------------------------
      implicit none
c  define the arguments
      integer l,lp,lamda
      double precision k,kp,alpha0,alpha2,q
c  internal variables
      double precision bkt,glamda,sqrkkp,a,pi
      integer s
c  we call rlam and glam
      double precision rlam,glam

c  We should ignore l=0,lp=0 since the radial integral diverges
c  ... but we won't for our testing...
c      if (l.eq.0 .and. lp.eq.0) return
c      write(6,*) 'Debug: l,lp,alpha0,k,kp=',l,lp,alpha0,k,kp
      if (l.eq.0 .and. lp.eq.0) then
         pi = 4.d0*atan(1.d0)
         bornkel = -pi*alpha0*sqrt(kp/k)/12.d0*(kp*kp+3.d0*k*k)
         return
      endif
      bkt = 0.d0
      sqrkkp = sqrt(k*kp)
c  lambda(potential) = 0
      glamda = glam(0,l,lp,lamda)
c      write(6,*) 'Debug: glambda for lambda = 0: l,lp,lam,glam =',l,lp,
c     1           lamda,glamda
      if (glamda .ne. 0.d0) then
         s = 4
         a = -alpha0/2.d0
         bkt = bkt - 2*sqrkkp*glamda*rlam(l,lp,k,kp,a,s)
      endif
c  lambda(potential) = 2
      glamda = glam(2,l,lp,lamda)
c      write(6,*) 'Debug: glambda for lambda = 2: l,lp,lam,glam =',l,lp,
c     1           lamda,glamda
      if (glamda .ne. 0.d0) then
         s = 4
         a = -alpha2/2.d0
         bkt = bkt - 2*sqrkkp*glamda*rlam(l,lp,k,kp,a,s)
         s = 3
         a = -q
         bkt = bkt - 2*sqrkkp*glamda*rlam(l,lp,k,kp,a,s)
      endif
c      write(6,*) 'Debug: returning from bornkel, bkt=',bkt
      bornkel = bkt
      return
      end

      double precision function rlam(lin,lpin,kin,kpin,al,s)
c-----------------------------------------------------------------
c  This calculates the analytic radial integral which appears in the
c  definition of the Born K matrix.  We only use the quadrupole moment
c  and polarizabilities
c-----------------------------------------------------------------
      implicit none
c  define the arguments
      integer lin,lpin,s
      double precision kin,kpin,al
c  internal variables
      logical first 
      DATA first /.true./

      integer maxel,l,lp
      parameter (maxel=36)
      integer it,l0,l1,l2,num,ii
      double precision gamma(maxel+2)
      double precision k,kp,kpow,pi,ratk,gam0,gam1,gam2,a,b,c
      double precision ratk2,ratklp,fnum,gamrat,f,tol
      save first,gamma
c  we call facx
      double precision facx
      k = kin
      kp = kpin
      l = lin
      lp = lpin
c  use the symmetry of the integral as needed
      if (kin .lt. kpin) then
         k = kpin
         kp = kin
         l = lpin
         lp = lin
      endif
      rlam = 0.d0
      if ((l+lp-s+3).le.0) return
      pi = 4.d0*atan(1.d0)
c  tol is used to terminate the series when calculating the
c  hypergeometric function
      tol = 1.d-8
c  here we fill the matrix gamma with the needed gamma function values
      if (first) call gam(maxel+1,gamma)
      first=.false.
      kpow=k**(s-3)
      ratk=kp/k
      it=l-lp+s
c      write(6,*) 'Debug: l,lp,s,ratk=',l,lp,s,ratk
      if (ratk.lt.0.99999d+00) then
         if (it/2*2.ne.it) then
            l0=(l+lp-s+1)/2
            l1=(it-1)/2
            gam0=facx(l0)
            gam1=gamma(l1+1)
         else
            l1=it/2-1
            l0=(l+lp-s)/2+1
            if (l0.gt.maxel) stop 'error in function rlam: l0 > maxel'
            gam1=facx(l1)
            gam0=gamma(l0+1)
         endif
         l2=lp+1
         if (l2.gt.maxel) stop 'error in function rlam: l2 > maxel'
         gam2=gamma(l2+1)
         a=float(l+lp-s+3)/2.d0
         b=float(lp-l-s+2)/2.d0
         c=float(lp+lp+3)/2.d0
         ratk2=ratk**2
         ratklp=ratk**lp
         call hyper(a,b,c,ratk2,f,tol)
         rlam=al*pi/(2**s)*k**(s-3)*ratklp*gam0/gam1/gam2*f
      else
         if (it/2*2.ne.it) then
            l1=(it-1)/2
            l2=(lp-l+s-1)/2
            gam1=gamma(l1+1)
            gam2=gamma(l2+1)
         else
            l1=it/2-1
            l2=(lp-l+s)/2-1
            gam1=facx(l1)
            gam2=facx(l2)
         endif
         num=l+lp+s-1
         fnum=1.d0
         do ii=1,s-1
            fnum=fnum*num
            num=num-2
         enddo
c	 write(6,*) 'Debug: fnum,gam1,gam2=',fnum,gam1,gam2
         gamrat=1.d0/fnum
         rlam=al*pi/2.d0*k**(s-3)*facx(s-2)*gamrat/gam1/gam2
      endif
      return
      end

      double precision function facx(i)
      implicit none
c-----------------------------------------------------------------
c     This is a factorial function  i!
c     It is done in floating point so large i can be used
c-----------------------------------------------------------------
      integer i
      double precision fj,table(15)
      integer j
      data table/1.0d+00,2.0d+00,6.0d+00,24.0d+00,120.0d+00,720.0d+00,
     1      5040.0d+00,40320.0d+00,362880.0d+00,
     2      36288.0d+02,399168.0d+02,4790016.0d+02,62270208.0d+02,
     3      871782912.0d+02,1307674368.0d+03/
      if (i.lt.0) then
         facx=0.d0
      else if (i.eq.0) then
         facx=1.0d0
      else if (i.le.15) then
         facx=table(i)
      else
         fj=16.0d+00
         facx=table(15)
         do j=16,i
            facx=facx*fj
            fj=fj+1.0d+00
         enddo
      endif
      return
      end

      subroutine gam(n,g)
      implicit none
c-------------------------------------------------------
c     This function computes the gamma function gamma(n+1/2)
c     and stores it in vector gamma(n+1).  one call to gam
c     with n set to the maximum is sufficient.
c-------------------------------------------------------
      integer n
      double precision g(n+1)
      integer nn,i,ii
      g(1)=sqrt(4.d0*atan(1.d0))
      if (n.lt.1) return
      nn=n+1
      do i=2,nn
         ii=i-1
         g(i)=(float(ii)-0.5d0)*g(ii)
      enddo
      return
      end

      subroutine hyper(a,b,c,x,f12,tol)
      implicit none
c-----------------------------------------------------------------
c        this routine calculates hypergeometric functions by
c        series sum until some term is less than tol times the
c        current summed value.
c-----------------------------------------------------------------
      double precision a,b,c,x,f12,tol
c  internal variables
      double precision f,pocha,pochb,pochc,pabc,xpow,term,plus
      integer i
      f=1.0d0+a/c*b*x
      pocha=a
      pochb=b
      pochc=c
      pabc=pocha/pochc*pochb
      i=2
      xpow=x
      term=100.0d0*tol*f
      do while (abs(term/f).gt.tol)
         plus=float(i-1)
         pabc=pabc*(a+plus)/(c+plus)*(b+plus)/float(i)
         xpow=xpow*x
         term=pabc*xpow
         f=f+term
         i=i+1
      enddo
      f12=f
      return
      end

      double precision function glam(lam,l,lp,lambda)
      implicit none
c-----------------------------------------------------------------
c     This returns the result of:
c        (-1)^lambda*sqrt((2l+1)*(2l'+1))*
c        f3j(l' l lam / 0 0 0) * f3j(l' l lam / -lambda lambda 0)
c     which is used in producing Born K matrix elements
c     We use the f3j routine and assume that facset has been
c     called by our main routine.
c-----------------------------------------------------------------
      integer lam,l,lp,lambda
      double precision f3j

      glam = (-1)**lambda*sqrt((2.d0*l+1.d0)*(2.d0*lp+1.d0))*
     $                    f3j(lp,l,lam,0,0,0)*
     $                    f3j(lp,l,lam,-lambda,lambda,0)
      return
      end
