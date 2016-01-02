C *-
C * PROGRAM R-CLOSURE
C *
C * This is a set of subroutines needed to implement r-closure
C * they are in a separate file because this aided the original debugging
C * 
C * V1.0: A. Feldt (OU), August, 2004
C *
C *-/

C * ====================================================================

C *-
C * This subroutine takes the kmatrix input at radial distance r
C * and adds the Born contribution from r to infinity to it and
C * returns the modified k matrix
C * as written, this considers homonuclear systems only with only
C * the quadrupole and polarizability terms used for the potential

      subroutine addbornk(kmatrix,k,r,l0,chin,nch,alpha0,alpha2,q
     $     ,vibdim,fmat,npw,nlamda)

      implicit none

C * argument declarations
      integer nch,l0,vibdim,npw,nlamda
CARR
      integer chanmax,vibmax,pwavemax,lammax
      parameter(chanmax=220,vibmax=20,pwavemax=11,lammax=21)
C      double precision kmatrix(nch,nch),k(nch),r
C      double precision alpha0(vibdim,vibdim),alpha2(vibdim,vibdim)
C      double precision q(vibdim,vibdim)
C      double precision fmat(npw,npw,nlamda)
C      integer chin(nch,2)
      double precision kmatrix(chanmax,chanmax),k(chanmax),r
      double precision alpha0(vibmax,vibmax),alpha2(vibmax,vibmax)
      double precision q(vibmax,vibmax)
      double precision fmat(pwavemax,pwavemax,lammax)
      integer chin(chanmax,2)

C * internal variables
      integer i,j,l,lp,il,ilp,indexv,indexvp,v,vp
      double precision temp,a0vvp,a2vvp,qvvp,fullborn
      double precision bornkadd(chanmax,chanmax)

C * we call this function
      double precision bornkel

      write(6,*)'Born r-closure from r = ',r,' to infinity K matrices:'
      write(6,*)' v vp  l lp Born: 0 -> Inf Born: r -> Inf final kmat'

C * fill the part of the matrix where kp .le. k (as required by
C * the radial integral)

      do i = 1,nch
         indexv = chin(i,1)
         v  =  indexv-1
         il  =  chin(i,2)
         l = 2*(il-1)+l0
         do j = 1,nch
            if (k(j).le.k(i)) then
               indexvp = chin(j,1)
               vp  =  indexvp-1
               a0vvp = alpha0(indexv,indexvp)
               a2vvp = alpha2(indexv,indexvp)
               qvvp = q(indexv,indexvp)
               ilp = chin(j,2)
               lp = 2*(ilp-1)+l0
               bornkadd(i,j) = bornkel(il,ilp,l,lp,k(i),k(j),r,a0vvp
     $              ,a2vvp,qvvp,fmat,npw,nlamda,fullborn)
               kmatrix(i,j)  =  kmatrix(i,j) + bornkadd(i,j)
               write(6,500)v,vp,l,lp,fullborn,bornkadd(i,j),kmatrix(i,j)
            endif
         end do
      end do

C * now supplement the rest of the matrix by symmetry
      do i = 1,nch
         do j = 1,nch
            if (k(i).lt.k(j)) then
               kmatrix(i,j)  =  kmatrix(i,j) + bornkadd(j,i)
            endif
         enddo
      enddo
 500  format(4i3,3(1pe12.4))
      return
      end

C *-
C * We return the Born K matrix element contribution from r to Inf
C * We also put the full (r=0 to Inf) analytic result in fullborn

      double precision function bornkel(il,ilp,l,lp,k,kp,r,alpha0
     $     ,alpha2,q,fmat,npw,nlamda, fullborn)

      implicit none

C * define the arguments
      integer il,ilp,l,lp,npw,nlamda
      double precision k,kp,r,alpha0,alpha2,q,fullborn
CARR
      integer pwavemax,lammax
      parameter(pwavemax=11,lammax=21)
C      double precision fmat(npw,npw,nlamda)
      double precision fmat(pwavemax,pwavemax,lammax)

C *  internal variables
      double precision glam,sqrkkp,rint,abserr
      double precision rmin,epsabs,epsrel,a,kc,kpc,rlanalytic
      integer lambda,s,limit,lenw,ier,neval,last,lc,lpc
      parameter (limit = 1000,lenw=4000)
      double precision work(lenw),iwork(limit)

C * We have to have a common to get info into rlami which must
C * have a single argument (r) to be passed into the integrator
      common /rlparams/ kc,kpc,a,s,lc,lpc

C * We call rlam and rlami
      double precision rlam
      external rlami

C * fmat used here is equal to (-1)^m*sqrt((2l+1)*(2l'+1))*
C * f3j(l' l lam / 0 0 0) * f3j(l' l lam / -m m 0)
C * (It can be shown that the form in fmat can be transformed
C * to the above which is the g_lamda in my derivations.)

      rmin = 0.d0
      epsabs = 0.d0
      epsrel = 1.d-6
      bornkel = 0.d0
      fullborn = 0.d0

C * Here we set part of the common rlparams values needed in the
C * function rlami.  The other two (a & s) are set as needed for
C * each lamba and potential part later.
      kc = k
      kpc = kp
      lc = l
      lpc = lp

C * We must ignore l=0,lp=0 since the radial integral diverges
      if (l.eq.0 .and. lp.eq.0) return
      sqrkkp = sqrt(k*kp)

C * lambda = 0
      glam = fmat(il,ilp,1)
      if (glam .ne. 0.d0) then
         s = 4
         a = -alpha0/2.d0
C * dqags is an adaptive quadrature routine.  It integrates the
C * function in the first argument between the min and max values
C * given in the second and third arguments.  See the comments
C * in dqags for further information.
         call dqags(rlami,rmin,r,epsabs,epsrel,rint,abserr,neval,ier
     $        ,limit,lenw,last,iwork,work)
         if (ier.ne.0) stop 'alpha0 - ier from dqags is non-zero'
         rlanalytic = rlam(l,lp,k,kp,a,s)
         fullborn = fullborn - 2*sqrkkp*glam*rlanalytic
         bornkel = bornkel - 2*sqrkkp*glam*(rlanalytic-rint)
      endif

C * lambda = 2
      glam = fmat(il,ilp,2)
      if (glam .ne. 0.d0) then
         s = 4
         a = -alpha2/2.d0
         call dqags(rlami,rmin,r,epsabs,epsrel,rint,abserr,neval,ier
     $        ,limit,lenw,last,iwork,work)
         if (ier.ne.0) stop 'alpha2 - ier from dqags is non-zero'
         rlanalytic = rlam(l,lp,k,kp,a,s)
         fullborn = fullborn - 2*sqrkkp*glam*rlanalytic
         bornkel = bornkel - 2*sqrkkp*glam*(rlanalytic-rint)
         s = 3
         a = -q
         call dqags(rlami,rmin,r,epsabs,epsrel,rint,abserr,neval,ier
     $        ,limit,lenw,last,iwork,work)
         if (ier.ne.0) stop 'q - ier from dqags is non-zero'
         rlanalytic = rlam(l,lp,k,kp,a,s)
         fullborn = fullborn - 2*sqrkkp*glam*rlanalytic
         bornkel = bornkel - 2*sqrkkp*glam*(rlanalytic-rint)
      endif
      return
      end

C *-
C * This calculates the analytic radial integral which appears in the
C * definition of the Born K matrix.  We assume that kp <= k here
C * and, we only use the quadrupole moment and polarizabilities

      double precision function rlam(l,lp,k,kp,al,s)
      implicit none

C * define the arguments
      integer l,lp,s
      double precision k,kp,al

C * internal variables
      logical first 
      DATA first /.true./

      integer maxel
      parameter (maxel=96)
      integer it,l0,l1,l2,num,ii
      double precision gamma(maxel+3)
      double precision kpow,pi,ratk,gam0,gam1,gam2,a,b,c
      double precision ratk2,ratklp,fnum,gamrat,f,tol
      save first,gamma

C * we call facx
      double precision facx
      if (k .lt. kp) then
         stop 'error in function rlam:  k < kp'
      endif
      pi = 4.d0*atan(1.d0)

C * tol is used to terminate the series when calculating the
C * hypergeometric function
      tol = 1.d-8

C * here we fill the matrix gamma with the needed gamma function values
      if (first) call gam(39,gamma)
      first = .false.
      kpow = k**(s-3)
      ratk = kp/k
      it = l-lp+s
      if (ratk.lt.0.99999d+00) then
         if (it/2*2.ne.it) then
            l0 = (l+lp-s+1)/2
            l1 = (it-1)/2
            gam0 = facx(l0)
            gam1 = gamma(l1+1)
         else
            l1 = it/2-1
            l0 = (l+lp-s)/2+1
            if (l0.gt.maxel) stop 'error in function rlam: l0 > maxel'
            gam1 = facx(l1)
            gam0 = gamma(l0+1)
         endif
         l2 = lp+1
         if (l2.gt.maxel) stop 'error in function rlam: l2 > maxel'
         gam2 = gamma(l2+1)
         a = float(l+lp-s+3)/2.d0
         b = float(lp-l-s+2)/2.d0
         c = float(lp+lp+3)/2.d0
         ratk2 = ratk**2
         ratklp = ratk**lp
         call hyper(a,b,c,ratk2,f,tol)
         rlam = al*pi/(2**s)*k**(s-3)*ratklp*gam0/gam1/gam2*f
      else
         if (it/2*2.ne.it) then
            l1 = (it-1)/2
            l2 = (lp-l+s-1)/2
            gam1 = gamma(l1+1)
            gam2 = gamma(l2+1)
         else
            l1 = it/2-1
            l2 = (lp-l+s)/2-1
            gam1 = facx(l1)
            gam2 = facx(l2)
         endif
         num = l+lp+s-1
         fnum = 1.d0
         do ii = 1,s-1
            fnum = fnum*num
            num = num-2
         enddo
         gamrat = 1.d0/fnum
         rlam = al*pi/2.d0*k**(s-3)*facx(s-2)*gamrat/gam1/gam2
      endif
      return
      end

C *-
C * This is a factorial function  i!
C * It is done in floating point so large i can be used

      double precision function facx(i)
      implicit none

      integer i
      double precision fj,table(15)
      integer j
      data table/1.0d+00,2.0d+00,6.0d+00,24.0d+00,120.0d+00,720.0d+00,
     $     5040.0d+00,40320.0d+00,362880.0d+00,36288.0d+02,399168.0d+02
     $     ,4790016.0d+02,62270208.0d+02,871782912.0d+02
     $     ,1307674368.0d+03/

      if (i.lt.0) then
         facx = 0.d0
      else if (i.eq.0) then
         facx = 1.0d0
      else if (i.le.15) then
         facx = table(i)
      else
         fj = 16.0d+00
         facx = table(15)
         do j = 16,i
            facx = facx*fj
            fj = fj+1.0d+00
         enddo
      endif
      return
      end

C *-
C * This function computes the gamma function gamma(n+1/2)
C * and stores it in vector gamma(n+1).  one call to gam
C * with n set to the maximum is sufficient.

      subroutine gam(n,g)
      implicit none

      integer n
      double precision g(n+1)
      integer nn,i,ii

      g(1) = sqrt(4.d0*atan(1.d0))
      if (n.lt.1) return
      nn = n+1
      do i = 2,nn
         ii = i-1
         g(i) = (float(ii)-0.5d0)*g(ii)
      enddo
      return
      end

C *-
C * this routine calculates hypergeometric functions by
C * series sum until some term is less than tol times the
C * current summed value.

      subroutine hyper(a,b,c,x,f12,tol)
      implicit none

      double precision a,b,c,x,f12,tol

C * internal variables
      double precision f,pocha,pochb,pochc,pabc,xpow,term,plus
      integer i

      f = 1.0d0+a/c*b*x
      pocha = a
      pochb = b
      pochc = c
      pabc = pocha/pochc*pochb
      i = 2
      xpow = x
      term = 100.0d0*tol*f
      do while (abs(term/f).gt.tol)
         plus = float(i-1)
         pabc = pabc*(a+plus)/(c+plus)*(b+plus)/float(i)
         xpow = xpow*x
         term = pabc*xpow
         f = f+term
         i = i+1
      enddo
      f12 = f
      return
      end

C *-
C * This is the integrand of the radial integral used to determine
C * the Born K matrix

      double precision function rlami(r)

      implicit none
      double precision k,kp,r,a
      integer s,l,lp
      double precision besselj1,besselj2,tmp1,tmp2,tmp3
      double precision pi,arg1,arg2
      common /rlparams/ k,kp,a,s,l,lp

      rlami = 0.d0
      if (r .le. 0.d0) return
      pi = 4*atan(1.0d0)
      arg1 = kp*r
      arg2 = float(lp)+0.5d0
      call bessjy(arg1,arg2,besselj1,tmp1,tmp2,tmp3)
      arg1 = k*r
      arg2 = float(l)+0.5d0
      call bessjy(arg1,arg2,besselj2,tmp1,tmp2,tmp3)
      rlami = a*pi/(2.d0*sqrt(k*kp))*besselj1*besselj2*r**(1-s)
      return
      end
