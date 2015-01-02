C * SEP 01 2005 - AF  - CAL. DCS USING ABCVIB FORMALISM
c
c  This is a program to calculate the differential cross-section for
c  a homonuclear system utilizing the abcvib formalism of Feldt
c  and Morrison (extended from the elastic abc formalism of Isaacs
c  and Morrison).
c
      implicit none

      integer lmax
c
c  ldim is the largest value of lmax we will allow by any dimensioning
c  this will limit l,l0,Lambda, etc. as well
c
      integer ldim
      parameter (ldim=20)
c
c  these are the T matrices - we expect them as T(Lam,l,l0)
c
      double precision tre(0:ldim,0:ldim,0:ldim)
      double precision tim(0:ldim,0:ldim,0:ldim)

      integer lsymmax(0:ldim)
      double precision worst(2,0:ldim)
      double precision k0,k,alpha0,alpha2,q
      double precision d(2)
      double precision Einp
      integer tform,debug,noabc,lbext,lend,i,i2,i3,nchs
      character*75 title
c
c  we always want the dcs at degree intervals, but will allow
c  parameterization just in case
c
      integer ntheta
      parameter (ntheta=181)
      double precision dcs(0:ntheta-1)
      double precision mominteg(0:ntheta-1)
c
c  note that setup must be called before we ever use the thrj routine
c  see the note in clebg.f and we call thrj from our f3j utility
c
      call setup
c
c  read the title and then
c  read the six basic inputs - the rest depend on how we input
c  the T-matrices
c      
c  lmax is the maximum l that is *not* Born - use -1 to obtain
c       pure Born results. -- but see noabc = 2 below
c
c  tform = 0 for using internal Born t-matrices (for testing)
c          1 for reading standard lavib K-matrices.
c          2 for reading bvfcc K-matrices
c        Other forms can be added later if needed.
c
c  noabc is a flag, which has the following meanings:
c  noabc = 0 -> perform the abcv cross section calculation
c          1 -> turn off abcv but still sum to a maximum
c               partial wave
c          2 -> turn off abcv and sum to a maximum symmetry
c               for all available partial waves
c               *** For noabc = 2, we change the meaning
c               *** of lmax and lbext so that lmax is the
c               *** maximum Lambda (symmetry) for which
c               *** we use c.c. T-matrices and lbext is
c               *** the maximum Lambda (symmetry) to which
c               *** we extend these with Born values
c               See also nchs
c
c  lbext is the maximum value of l to extend our T matrices
c        with Born values.  This should normally be set to lmax
c        so that we automatically fill in any missing c.c.
c        elements.  Be aware that you must verify that lmax
c        is large enough by checking the 'goodness' values.
c        It can also be usefull when noabc is non-zero.
c        Also, see noabc = 2 for special meaning for lbext then
c  
c  debug can be set non-zero for debugging purposes (lots of
c        extra output)
c
c  nchs is equal to the number of partial waves to be used for
c       each symmetry, intended for use with noabc = 2, but
c       will restrict the c.c. input in any case if non-zero
c       When noabc = 2 and nchs is set to 0, each symmetry is
c       filled to all partial waves implied by the maximum
c       symmetry.
c
      read(5,'(a75)') title
      read(5,*) lmax,tform,noabc,lbext,debug,nchs
      write(6,*) 'Abcvib differential cross section code'
      write(6,*)
      write(6,'(a75)') title
      write(6,*)
      write(6,*) 'Input values:'
      write(6,*) '   lmax   = ',lmax
      write(6,*) '   tform  = ',tform
      write(6,*) '   noabc  = ',noabc
      if (noabc .ne. 0) then
         write(6,*) '   No abc formalism will be used'
         if (noabc .eq. 1) then
            write(6,*) ' Summations will be limited by partial wave'
         else if (noabc .eq. 2) then
            write(6,*) ' Summations will be limited by maximum symmetry'
         endif
      endif
      write(6,*) '   lbext  = ',lbext
      write(6,*) '   debug  = ',debug
      write(6,*) '   nchs   = ',nchs
      if (nchs .eq. 0 .and. noabc .eq. 2) then
         write(6,*) '   All l values to the maximum Lambda will be used'
      else if (nchs .gt. 0) then
         write(6,*) '   Will limit each symmetry to nchs channels'
      endif

c
c  zero our matrices - not all compilers do it for us
c
      lend = max(lmax,lbext)
      do i=0,lend
         do i2=0,lend
            do i3=0,lend
               tre(i,i2,i3) = 0.d0
               tim(i,i2,i3) = 0.d0
            enddo
         enddo
      enddo

      call readtmats(Einp,d,tre,tim,ldim,lmax,tform,k0,k,alpha0,alpha2,q
     $     ,lbext,lsymmax, worst,noabc,debug,nchs)

      if (noabc .ne. 0) then
         alpha0 = 0.d0
         alpha2 = 0.d0
         q = 0.d0
      endif

      if (noabc .ne. 2) then
         call abcvib(lend,k0,k,alpha0,alpha2,q,tre,tim,ldim, dcs,ntheta
     $        ,noabc,debug)
      else
         if (debug .ne. 0) then
            do i=0,lend
               write(6,*) 'Debug: Lambda,lsymmax(Lambda)=',i,lsymmax(i)
            enddo
         endif
         call tradxsec(lend,k0,k,tre,tim,ldim, dcs,ntheta,lsymmax,debug)
      endif
      call results(Einp,lend,k0,k,alpha0,alpha2,q,d,worst, dcs,mominteg
     $     ,ntheta,tform)

      end

      subroutine readtmats(Einp,d,tre,tim,ldim,lmax,tform,k0,k,alpha0
     $     ,alpha2,q,lbext,lsymmax,worst,noabc,debug,nchs)
c
c  read and store the T-matrix elements
c
      implicit none
      integer ldim,lmax,tform,lbext,noabc,debug,nchs
      double precision tre(0:ldim,0:ldim,0:ldim)
      double precision tim(0:ldim,0:ldim,0:ldim)
      integer i,lsymmax(0:ldim)
      double precision worst(2,0:ldim)
      double precision k0,k,alpha0,alpha2,q
      double precision d(2)
      double precision Einp

c
c  initialize lsymmax to zero
c
      do i=0,ldim
         lsymmax(i) = 0
      enddo
c
c  d is the constant defined in Mike's 'Gun and Camera' document
c       S = 1 - 2 Pi i d T
c  and will automatically be set to   -i/(2 Pi) for the internal
c  Born T matrices and for our lavib code T matrices
c  d(1) is the real part of d
c  d(2) is the imaginary part of d
c  we will actually return d*T in the T-matrix used by the
c  abcvdcs routine.  This allows us to read in a T-matrix
c  with one normalization and still augment it with Born
c  values if we are turning off the abcvib formalism and
c  just using the traditional form.
c
      if (tform.eq.0) then
c
c  k0, k are the channel wavenumbers
c  alpha0, alpha2 and q are the moments (in a.u.)
c
         read(5,*) k0,k
         read(5,*) alpha0,alpha2,q
         write(6,*) '   Born T-matrices used. (tform=0)'
         call borntmats(d,tre,tim,ldim,0,lmax,lsymmax,worst,k0,k,alpha0
     $        ,alpha2,q,noabc,debug,nchs)
      else if (tform .eq. 1 .or. tform .eq. 2) then
c
c  this is for the K matrices output by the lavib/bfvcc codes.
c  they are almost identical, so we use the same basic routine.
c
         if (tform .eq. 1) then
            write(6,*) '   lavib T-matrices used. (tform=1)'
         else
            write(6,*) '   bfvcc T-matrices used. (tform=2)'
         endif
         write(6,*)
         write(6,*) 'Further input for this form:'
         write(6,*)
         call lavibtmats(Einp,d,tre,tim,ldim,lmax,lsymmax,worst,k0,k
     $        ,alpha0,alpha2,q,lbext,noabc, tform,debug,nchs)
      else
        stop 'Unknown T-matrix form.  Please write a new module!'
      endif
c
c  now verify that we have all the T matrices we need
c
      call verifyt(tre,tim,ldim,lmax)

      return
      end

      subroutine abcvib(lmax,k0,k,alpha0,alpha2,q,tre,tim,ldim,dcs
     $     ,ntheta,noabc,debug)
c
c  This calculates the differential cross section based on
c  the abcvib formalism - the only difference from the
c  notes is that we assume that the constant 'd' has been
c  absorbed into the T-matrix prior to this routine.  Hence,
c  we won't refer to 'd' here.
c  When noabc is no-zero, we will only use the T-matrices
c  given, but will still limit our sums by the maximum partial
c  wave.  See tradxsec for a routine limiting the sums based
c  on the maximum symmetry.
c
      implicit none
      integer lmax,ldim,ntheta,noabc,debug
      double precision k0,k,alpha0,alpha2,q
      double precision tre(0:ldim,0:ldim,0:ldim)
      double precision tim(0:ldim,0:ldim,0:ldim)
      double precision dcs(0:ntheta-1)

      double precision term3,term4,term5,term6

      integer i,m
      double precision theta,pi,f0t,f2t(-2:2)
      double precision f0,f2
      double precision dt5

      pi = 4.d0*atan(1.d0)
      do i = 0,ntheta-1
         theta = i*pi/(ntheta-1)
         dcs(i) = 0.d0
         if (noabc .eq. 0) then
            f0t = f0(pi,alpha0,k,k0,lmax,theta)
            do m = -2,2
               f2t(m) = f2(m,pi,alpha2,q,k,k0,lmax,theta)
            enddo
            if (debug.gt.0) then
               write(6,*) 'Debug: f0=',f0t
               write(6,*) 'Debug: f2=',f2t
               write(6,*) 'Debug: term4=',term4(f2t)
               dt5 = term5(f2t,pi,k,k0,tre,tim,ldim,lmax,theta)
               write(6,*) 'Debug: term5=',dt5
               write(6,*) 'Debug: term6=', term6(pi,k,k0,tre,tim,ldim
     $              ,lmax,theta)
            endif
            dcs(i) = f0t**2 + term3(f0t,pi,k,k0,tre,tim,ldim,lmax,theta)
     $                      + term4(f2t) 
     $                      + term5(f2t,pi,k,k0,tre,tim,ldim,lmax,theta)
         endif
         dcs(i) = dcs(i) + term6(pi,k,k0,tre,tim,ldim,lmax,theta)
         dcs(i) = k/k0 * dcs(i)
         if (dcs(i) .lt. 0.d0) dcs(i) = 0.d0
      enddo
      return
      end

      double precision function f0(pi,alpha0,k,k0,lmax,theta)
c
c  From the derivation notes - f0 depends only on the lambda=0
c  component of the long range potential
c
      implicit none
      integer lmax,l
      double precision pi,alpha0,k,k0,theta
      double precision plgndr,smallq
      double precision f0hold

      f0hold = smallq(k,k0,theta)/4
      if (lmax.ge.0) then
         do l=0,lmax
            f0hold = f0hold + k0/4*(k/k0)**l*(1.d0/(2.d0*l-1.d0)-(k/k0)
     $           **2/(2.d0*l+3.d0))* plgndr(l,0,cos(theta))
         enddo
      endif
      f0 = pi*alpha0*f0hold
      return
      end

      double precision function f2(m,pi,alpha2,q,k,k0,lmax,theta)
c
c  From the derivation notes - f2 depends only on the lambda=2
c  component(s) of the long range potential
c
      implicit none
      integer m,lmax
      double precision pi,alpha2,q,k,k0,theta
      double precision f2hold
      double precision ylm,smallq,thetaq,f3j,bigm
      double precision fjs
      integer l,l0

      f2hold = sqrt(4.d0*pi/5.d0)*(-1)**m*ylm(2,m,thetaq(pi,k,k0
     $     ,theta))* (alpha2*pi*smallq(k,k0,theta)/16.d0 + 2.d0*q/3.d0)
      if (lmax.ge.0) then
         do l = abs(m),lmax
            do l0 = 0,lmax
               if (mod(l-l0,2).eq.0) then
                  fjs = f3j(l,l0,2,0,0,0)*f3j(l,l0,2,-m,0,m)
                  if (fjs .ne. 0.d0) then
                     f2hold = f2hold + sqrt(4.d0*pi)*(-1)**((l-l0)/2)*(
     $                    -1)**m*sqrt((2.d0*l+1.d0))*(2.d0*l0+1.d0)* fjs
     $                    * (alpha2*bigm(4,l,l0,k,k0)+ 2.d0*q*bigm(3,l
     $                    ,l0,k,k0))* ylm(l,m,theta)
                  endif
               endif
            enddo
         enddo
      endif
      f2 = f2hold
      return
      end

      double precision function bigm(eta,l,l0,k,k0)
c
c  Note that I am using 'bigm' to match the derivation notes, but
c  we already have the rlam function used to produce the test
c  Born matrices and they are essentially the same (requiring
c  only a 'moment' of 1 in rlam to make them match)
c  
      implicit none
      double precision k,k0
      integer eta,l,l0
      double precision rlam
      bigm = rlam(l0,l,k0,k,1.d0,eta)
      return
      end

      double precision function littlec(l,l0,m,lam)
c
c  From the derivation notes - a "little" c
c
      implicit none
      integer l,l0,m,lam
      double precision f3j
      littlec = (-1)**(m-lam+(l0-l)/2)*sqrt(2.d0*l0+1)*f3j(l,l0,2,m,0,
     $     -m)*f3j(l,l0,2,lam,-lam,0)
      return
      end

      double precision function smallq(k,k0,theta)
c
c  The wavevector equivalent of the momentum transfer
c  This is the magnitude
c
      implicit none
      double precision k,k0,theta
      double precision radical
      radical = k0**2-2.d0*k*k0*cos(theta)+k**2
      if (radical.lt.0.d0) radical = 0.d0
      smallq = sqrt(radical)
      return
      end

      double precision function thetaq(pi,k,k0,theta)
c
c  and this is the angle associated with q
c
      implicit none
      double precision pi,k,k0,theta
      double precision term
      double precision smallq
      if (k.ne.k0) then
         term = (k0-k*cos(theta))/smallq(k,k0,theta)
         if (term .lt. -1.d0) term=-1.d0
         if (term .gt. 1.d0) term=1.d0
         thetaq = acos(term)
      else
         thetaq = 0.5d0*(pi-theta)
      endif
      return
      end

      double precision function term3(f0t,pi,k,k0,tre,tim,ldim,lmax
     $     ,theta)
c
c  The third term from the derivation - it has cross terms
c  between the lamda=0 part of the Born potential and
c  the full T-matrix elements
c
      implicit none
      double precision f0t,pi,k,k0,theta
      integer ldim,lmax
      double precision tre(0:ldim,0:ldim,0:ldim)
      double precision tim(0:ldim,0:ldim,0:ldim)
      double precision plgndr
      double precision t3,temp,tmatr,tmati
      integer l,lam
      t3 = 0.d0
      if (lmax.ge.0) then
         do l = 0,lmax
            temp = 2.d0*pi*f0t/sqrt(k*k0)*plgndr(l,0,cos(theta))
            do lam = -l,l
               if (lam .lt. 0) then
                  tmatr = tre(-lam,l,l)
                  tmati = tim(-lam,l,l)
               else
                  tmatr = tre(lam,l,l)
                  tmati = tim(lam,l,l)
               endif
               t3 = t3 + temp*tmatr
            enddo
         enddo
      endif
      term3 = t3
      return
      end

      double precision function term4(f2t)
c
c  The fourth term from the derivation - it depends only on
c  the lambda=2 part of the Born potential.
c
      implicit none
      double precision f2t(-2:2)
      double precision t4
      integer m
      t4 = 0.d0
      do m = -2,2
         t4 = t4 + f2t(m)**2/5.d0
      enddo
      term4 = t4
      return
      end

      double precision function term5(f2t,pi,k,k0,tre,tim,ldim,lmax
     $     ,theta)
c
c  The fifth term from the derivation - it has cross terms
c  between the lamda=2 part of the Born potential and
c  the full T-matrix elements
c
      implicit none
      double precision f2t(-2:2),pi,k,k0,theta
      integer ldim,lmax
      double precision tre(0:ldim,0:ldim,0:ldim)
      double precision tim(0:ldim,0:ldim,0:ldim)
      double precision t5,tmatr,tmati,pre
      integer m,l,l0,lam,lammax,mmax
      double precision ylm,littlec
      t5 = 0.d0
      if (lmax.ge.0) then
         pre = 4.d0*pi*sqrt(pi/k/k0)
         do l = 0,lmax
            mmax = min(l,2)
            do l0 = 0,lmax
               if (mod(l-l0,2).eq.0) then
                  lammax = min(l,l0)
                  do lam = -lammax,lammax
                     if (lam .lt. 0) then
                        tmatr = tre(-lam,l,l0)
                        tmati = tim(-lam,l,l0)
                     else
                        tmatr = tre(lam,l,l0)
                        tmati = tim(lam,l,l0)
                     endif
                     do m = -mmax,mmax
                        t5 = t5 + pre*tmatr*f2t(m)*ylm(l,m,theta)*
     $                       littlec(l,l0,m,lam)
                     enddo
                  enddo
               endif
            enddo
         enddo
      endif
      term5 = t5
      return
      end

      double precision function term6(pi,k,k0,tre,tim,ldim,lmax,theta)
c
c  This is the term which corresponds to the traditional
c  differential cross section.  Note, however, that because
c  of the abcv formalism, we organize our sums so that they
c  are dictated by the maximum partial wave, not the maximum
c  symmetry.  To get a symmetry-truncated cross section, we
c  will use a different subroutine - tradxsec.
c
      implicit none
      double precision pi,k,k0,theta
      integer ldim,lmax
      double precision tre(0:ldim,0:ldim,0:ldim)
      double precision tim(0:ldim,0:ldim,0:ldim)
      double precision t6,pre,lfac
      integer l,l0,lp,l0p,lam,lamp,bigl
      integer lammax,lampmax,biglmax,biglmin
      double precision c1,c2,c3,c4,tstart,plgndr
      double precision f3j
      t6 = 0.d0
      if (lmax.ge.0) then
         pre = pi**2 / k / k0
         do l = 0,lmax
            do l0 = 0,lmax
               if (mod(l+l0,2).eq.0) then
                  lammax = min(l,l0)
                  do lp = 0,lmax
                     do l0p = 0,lmax
                        if (mod(lp+l0p,2).eq.0) then
                           biglmax = min(l0+l0p,l+lp)
                           biglmin = max(abs(l0-l0p),abs(l-lp))
                           lampmax = min(lp,l0p)
                           lfac = (-1)**((l0p-lp-l0+l)/2)* sqrt((2.d0*l0
     $                          +1.d0)*(2.d0*l0p+1.d0)* (2.d0*l+1.d0)
     $                          *(2.d0*lp+1.d0))
                           do lam = -lammax,lammax
                              do lamp = -lampmax,lampmax
                                 tstart = tre(abs(lam),l,l0)
     $                                *tre(abs(lamp),lp,l0p)
     $                                +tim(abs(lam),l,l0)*tim(abs(lamp)
     $                                ,lp,l0p)
                                 do bigl = biglmin,biglmax
                                    c1 = f3j(l0,l0p,bigl,lam,-lamp,lamp
     $                                   -lam)
                                    c2 = f3j(l,lp,bigl,lam,-lamp,lamp
     $                                   -lam)
                                    c3 = f3j(l0,l0p,bigl,0,0,0)
                                    c4 = f3j(l,lp,bigl,0,0,0)
                                    if (c1*c2*c3*c4 .ne. 0.d0) then
                                       t6 = t6+pre*lfac*(2*bigl+1)*c1*c2
     $                                      *c3*c4*tstart*plgndr(bigl,0
     $                                      ,cos(theta))
                                    endif
                                 enddo
                              enddo
                           enddo
                        endif
                     enddo
                  enddo
               endif
            enddo
         enddo
      endif
      term6 = t6
      return
      end

      subroutine tradxsec(lammax,k,k0,tre,tim,ldim,dcs,ntheta,lsymmax
     $     ,debug)
c
c  This is the traditional differential cross section.
c  The summations here are symmetry-limited rather than
c  partial wave-limited.
c
      implicit none
      double precision k,k0,theta
      integer lammax,ldim,ntheta,debug
      double precision tre(0:ldim,0:ldim,0:ldim)
      double precision tim(0:ldim,0:ldim,0:ldim)
      double precision dcs(0:ntheta-1)
      integer lsymmax(0:ldim)
      double precision pre,lfac,costheta
      integer i,l,l0,lp,l0p,lam,lamp,bigl
      integer biglmax,biglmin,lmin,lpmin
      double precision pi,c1,c2,c3,c4,cs,tstart,plgndr
      double precision f3j
      pi = 4.d0*atan(1.d0)
      pre = pi**2 / k / k0
      do i = 0,ntheta-1
         theta = i*pi/(ntheta-1)
         costheta = cos(theta)
         dcs(i) = 0.d0
         do lam = -lammax,lammax
            lmin = abs(lam)
            do lamp = -lammax,lammax
               lpmin = abs(lamp)
               do l = lmin,lsymmax(lmin)
                  do l0 = lmin,lsymmax(lmin)
                     if (mod(l+l0,2).eq.0) then
                        do lp = lpmin,lsymmax(lpmin)
                           do l0p = lpmin,lsymmax(lpmin)
                              if (mod(lp+l0p,2).eq.0) then
                                 biglmax = min(l0+l0p,l+lp)
                                 biglmin = max(abs(l0-l0p),abs(l-lp))
                                 lfac = (-1)**((l0p-lp-l0+l)/2)*
     $                                sqrt((2.d0*l0+1.d0)*(2.d0*l0p
     $                                +1.d0)* (2.d0*l+1.d0)*(2.d0*lp
     $                                +1.d0))
                                 tstart = tre(lmin,l,l0)*tre(lpmin,lp
     $                                ,l0p)+tim(lmin,l,l0)*tim(lpmin,lp
     $                                ,l0p)
                                 do bigl = biglmin,biglmax
                                    c1 = f3j(l0,l0p,bigl,lam,-lamp,lamp
     $                                   -lam)
                                    c2 = f3j(l,lp,bigl,lam,-lamp,lamp
     $                                   -lam)
                                    c3 = f3j(l0,l0p,bigl,0,0,0)
                                    c4 = f3j(l,lp,bigl,0,0,0)
                                    cs = c1*c2*c3*c4
                                    if (cs .ne. 0.d0) then
                                       dcs(i) = dcs(i)+pre*lfac*(2*bigl
     $                                      +1)*cs*tstart*plgndr(bigl,0
     $                                      ,costheta)
                                    endif
                                 enddo
                              endif
                           enddo
                        enddo
                     endif
                  enddo
               enddo
            enddo
         enddo
         dcs(i) = k/k0 * dcs(i)
         if (dcs(i) .lt. 0.d0) dcs(i) = 0.d0
      enddo
      return
      end

      subroutine results(Einp,lmax,k0,k,alpha0,alpha2,q,d,worst,dcs
     $     ,mominteg,ntheta,tform)
c
c  the fruit of our labors...
c
      implicit none
      integer lmax,ntheta,tform
      double precision k0,k,alpha0,alpha2,q
      double precision d(2)
      double precision Einp
      double precision dcs(0:ntheta-1)
      double precision mominteg(0:ntheta-1)
      double precision worst(2,0:lmax)

      integer i
      double precision pi,theta,xsec,momxsec,integrate

      open(10,file='integral_xsc',status='new',form='formatted')
      open(11,file='momentum_msc',status='new',form='formatted')

      write(6,*)
      write(6,*) 'Summary of the key abcv parameters:'
      write(6,*)
      write(6,*) '   k0     = ',k0
      write(6,*) '   k      = ',k 
      write(6,*) '   q      = ',q 
      write(6,*) '   alpha0 = ',alpha0 
      write(6,*) '   alpha2 = ',alpha2 
      write(6,*) '   Re(d)  = ',d(1)
      write(6,*) '   Im(d)  = ',d(2)
      write(6,*)
      write(6,*) '   Maximum l (or Lambda) value used in sums = ',lmax
      write(6,*)
      write(6,*) '   The worst "goodness" for each l value:'
      write(6,*) ' l lam=0  other lambdas'
      do i=1,lmax
         write(6,'(i2,1x,f6.3,1x,f6.3)') i,worst(1,i),worst(2,i)
      enddo
      write(6,*)

      pi = 4.d0*atan(1.d0)
      do i=0,ntheta-1
         theta = i*pi/(ntheta-1)
         mominteg(i) = dcs(i)*cos(theta)
      enddo
      xsec = integrate(dcs,ntheta)
      momxsec = xsec - integrate(mominteg,ntheta)

      write(10,100) Einp, xsec
      write(11,100) Einp, momxsec

      write(6,*) 'Integrated cross section (a0**2): ',xsec
      write(6,*) 'Momentum transfer cross section (a0**2): ',momxsec
      write(6,*)
      write(6,*) 'Differential cross section:'
      write(6,*)
      write(6,*) 'Theta   DCS(a.u.)'
      do i = 0,ntheta-1
         theta = i*180.d0/(ntheta-1)
         write(6,'(f5.1,1x,g13.6)') theta,dcs(i)
         write(7,'(f5.1,1x,g13.6)') theta,dcs(i)
      enddo

 100  format(1x, f10.5, e20.10)
      return
      end

      double precision function integrate(dcs,ntheta)
c
c  We integrate the differential cross section * sin(theta)
c  using the trapezoidal rule.  Because we assume we are
c  always going from 0 to Pi, where sin(theta) is zero,
c  we ignore the endpoints...
c
      integer ntheta
      double precision dcs(0:ntheta-1)

      double precision theta,step,pi,xs
      pi = 4.d0*atan(1.d0)
      step = pi/(ntheta-1)
      xs = 0.d0
      do i = 1,ntheta-2
         theta = i*step
         xs = xs + dcs(i)*sin(theta)
      enddo
      integrate = xs*2.d0*pi*step
      return
      end

      double precision function f3j(j1,j2,j3,m1,m2,m3)
c
c  this utility routine allows us to write our calls
c  for 3j symbols more naturally since the thrj routine
c  require the angular momenta be twice their actual
c  value so that half-integer momenta can be represented
c  by integers
c
      implicit none
      integer j1,j2,j3,m1,m2,m3
      double precision thrj
      f3j = thrj(2*j1,2*j2,2*j3,2*m1,2*m2,2*m3)
      return
      end

      subroutine verifyt(tre,tim,ldim,lmax)
c
c  this routine attempts to verify that we have all the
c  necessary t matrix elements.  It will warn if not,
c  but will not stop the calculation - so watch the output!
c  we assume that if the necessary diagonal (in l) elements
c  have been found, the rest have, too
c
      implicit none
      integer ldim,lmax
      double precision tre(0:ldim,0:ldim,0:ldim)
      double precision tim(0:ldim,0:ldim,0:ldim)

      integer lam,symcount
      symcount = 0
      do lam=0,lmax
         if (tre(lam,lmax,lmax).ne.0.d0) symcount = symcount + 1
         if (lmax.gt.0) then
            if (tre(lam,lmax-1,lmax-1).ne.0.d0) symcount = symcount + 1
         endif
      enddo
      if (symcount .ne. (2*lmax+1)) then
         write(6,*) 'WARNING: only ',symcount,' symmetries found',
     $        ' but expected ',2*lmax+1
         if (symcount .le. 0) then
            stop 'No T matrix elements found to match the input'
         endif
      endif

      return
      end

      subroutine borntmats(d,tre,tim,ldim,lmin,lmax,lsymmax,worst,k0,k
     $     ,alpha0,alpha2,q,noabc,debug,nchs)
c
c  This calculates the Born T matrix elements.  It also assumes
c  that it could be called after the close coupling elements
c  are loaded.  If so, it then calculates how "close" to
c  the Born values they are by the "goodness" measure.
c  See the function goodness (gracious!)
c
      implicit none
      integer ldim,lmin,lmax,noabc,debug,nchs
      double precision tre(0:ldim,0:ldim,0:ldim)
      double precision tim(0:ldim,0:ldim,0:ldim)
      integer lsymmax(0:ldim)
      double precision worst(2,0:ldim)
      double precision k0,k,alpha0,alpha2,q
      double precision d(2)

      double precision bornkel,goodness,g,tr
      integer l0,l,lam,llmax,lg,lamg

      d(1) = 0.d0
      d(2) = -1.d0/(2.d0 * 4.d0*atan(1.d0))

      llmax = lmax
      do lam = 0,lmax
         if (noabc .eq. 2) then
            if (nchs .eq. 0) then
               llmax = min(lmax+1,ldim)
            else
               llmax = lam + 2*nchs - 1
            endif
         endif
         if (llmax .gt. lsymmax(lam)) then
            lsymmax(lam) = llmax
         endif
         if (debug.gt.0) then
            write(6,*) 'Debug: borntmats: lam,lsymmax(lam)=', lam
     $           ,lsymmax(lam)
         endif
         do l = lam,llmax
            do l0 = lam,llmax
               if (l.ge.lmin .or. l0.ge.lmin) then
                  if (tim(lam,l,l0) .eq. 0.d0 .and. tre(lam,l,l0) .eq.
     $                 0.d0) then
                     tre(lam,l,l0) = 2.d0*d(2)* bornkel(l0,l,k0,k,alpha0
     $                    ,alpha2,q,lam)
                     if (tre(lam,l,l0).ne.0.d0 .and. debug.gt.0) then
                        write(6,*) 'Debug: lam=',lam,' l=',l,' l0=',l0,
     $                       ' Re T=',tre(lam,l,l0)
                     endif
                  else
c  now we have a non-zero c.c. T-matrix element, how good is it?
c  we treat Lambda=0 separately because we don't really expect
c  Born for high energies for that symmetry
                     tr = 2.d0*d(2)*bornkel(l0,l,k0,k,alpha0,alpha2,q
     $                    ,lam)
                     if (tr .ne. 0.d0) then
                        lg = max(l,l0)
                        lamg = min(lam,1)+1
                        g = goodness(tre(lam,l,l0),tim(lam,l,l0),tr
     $                       ,0.d0)
                        if (g .gt. worst(lamg,lg)) worst(lamg,lg) = g
                        if (debug.gt.0) then
                           write(6,*) 'lam,l,l0,lg,g=',lam,l,l0,lg,g
                           write(6,*) 'tre,tim,tr=',tre(lam,l,l0),
     $                          tim(lam,l,l0),tr
                        endif
                     endif
                  endif
               endif
            enddo
         enddo
      enddo
      return
      end

      double precision function goodness(tccr,tcci,tbr,tbi)
c
c  This returns the ratio of the distance between the two T matrix
c  elements and the magnitude of the second (treating them in
c  the complex plane) as a measure of how "close" the first is
c  to the second.  The smaller this number, the closer they are.
c  
      implicit none
      double precision tccr,tcci,tbr,tbi
      goodness = sqrt(((tccr-tbr)**2+(tcci-tbi)**2)/(tbr**2+tbi**2))
      return
      end
