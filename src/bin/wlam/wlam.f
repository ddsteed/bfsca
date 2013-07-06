C * JUL 06 2013 - RDS - AUGMENT COUPLING NUMBER
C * MAY 29 2006 - RDS - DO E-E AND E-N COUPLING POTENTIALS SEPERATELY
C * DEC 31 2005 - RDS - AUGMENT THE ORDER OF GAUSS-HERMITE QUADRATURE FROM 20 -> 50
C *                     GENERALIZE NUMFUNCEX TO COPE WITH E-N2
C * JUN 18 1982 - AF  - VIB WFN CALCULATION SCAVANGED FROM PARKER'S VIVAS CODE
C *-
C * Code to calculate vibrational matrix elementts of:
C *   1)  the vlams (i.e. produce the wlams)
C *   2)  any numerical function
C *   3)  body frame t-matrices
C *   4)  any analytic function

      program nwlam
      implicit double precision (a-h,o-z)

C * logical variables
      logical iprint,ibasis,idbug,wrtvib,redvib,sho
      logical outvib

C * problem-dependent dimension statements

C * arrays dimensioned as vectors

C * the following should have dimension of nhermt
C     dimension energs(20),xsq(20),xk(20)
      dimension energs(30),xsq(30),xk(30)

C * the following should have length of at least nstats where nstats
C * is the total number of quantum states for the system regardless of
C * symmetry considerations (rotational times vibrational states.)

C * note that nstats must be specified in a data statement
      dimension mvib(60),jrot(60)
      dimension energy(60),msrot(60),msvib(60)

C * the following should have length of at least maxvib where
C * maxvib is the maximum number of oscillator basis
C * functions considered.

C     dimension jmin(20),jmax(20)
      dimension jmin(30),jmax(30)

C * the following should have length of at least nrot where nrot is
C * the maximum number of rotational j values.
      dimension numax(10),numin(10)

C * the following should have dimension of noscil*noscil
C     dimension tstore(400)
      dimension tstore(2500)

C * the following should have length of roughly nrot*nvib*nhermt,
C * where nvib is the number of vibrational states considered.
C     dimension chi(400)
      dimension chi(2500)

C * these are limited by the maximum nhermt allowed by
C * subroutine hermit
      dimension vquad(50),wquad(50),rquad(50)

      common/swchs/ibasis,idbug,wrtvib,redvib,sho
      common/readvib/outvib,lvib

      common/splin/ bpar(11),csplin(5000)
      common/prntop/iprint(23)
      common/silly/ vecold(20,20),noscil

C *   call link('unit5=(i5,open,text),unit6=(o6,create,text)//')
C *   call timeused(ocpu,osys,io,mem)

      level = 1
      nstats = 60
      zero = 0.0e+00

      call READIT(chi,iwlam,noff,nanl,inucl,nquad,rquad,wquad,minvib
     $     ,maxvib,tstore,vecold,numax,numin,jmax,jmin,xsq,xk,energs
     $     ,energy,msrot,jrot,msvib,mvib,nstats,noscil,ivibtr,ihoex)

      ivm = maxvib - minvib + 1
      if(ivm.ne.ivibtr) WRITE (6,*) ' truncating vibrational matrix by '
     $     ,ivibtr
      if (ihoex.ne.0) WRITE (6,*) ' expand kernel by har. osc. ',ihoex

      WRITE (6,9129)
      WRITE (6,9131) (rquad(i),i=1,nquad)

      if(iwlam.ne.0) call wlam(tstore,energs,xsq,xk,jmin,inucl,energy,
     $     chi,nquad,rquad,vquad,wquad,ivibtr)

      if(noff.gt.0) call numfun(noff,ivm,nquad,energy,vecold,chi,rquad
     $     ,vquad,wquad,minvib,ivibtr,ihoex)

      if(noff.lt.0) call tmatav(ivibtr,nquad,tstore,energy,vecold,chi
     $     ,rquad,vquad,wquad)

      if(nanl.ne.0) call anlfun(chi,rquad,wquad,nquad,ivm,minvib,tstore)

      stop

 9129 format(///' quadrature r-points (au)'/)
 9131 format(7x,f10.5)
      end

C *-
      SUBROUTINE READIT(chi,iwlam,noff,nanl,inucl,nhermt,xpth,wpth,
     $     minvib,maxvib,tstore,vecold,numax,numin,jmax,jmin,xsq,xk,
     $     energs,energy,msrot,jrot,msvib,mvib,nstats,noscil,ivibtr
     $     ,ihoex)

      implicit double precision (a-h,o-z)
      logical iprint,ibasis,idbug,wrtvib,redvib,sho
      logical outvib

      character title(72)

C * the following need only be given a single dimension

      dimension tstore(1),xsq(1),xk(1)

C * the following need to be dimensioned as they are in the main
C * program so that if redvib or wrtvib are true, the whole
C * matrix for each will be read or written

C     dimension chi(400),xpth(50),wpth(50),numax(10),numin(10)
      dimension chi(2500),xpth(50),wpth(50),numax(10),numin(10)

      dimension vecold(20,20)

C     dimension jmax(20),jmin(20),energs(20)
      dimension jmax(30),jmin(30),energs(30)

      dimension energy(60),msrot(60),jrot(60),msvib(60),mvib(60)

      common/swchs/ibasis,idbug,wrtvib,redvib,sho
      common/readvib/outvib,lvib

      common/masses/bmass,cmass,ubc,uml,umr

      common/need/rx,alpha,re

      common/diatom/we,wexe,be,alfe,de
      common/prntop/iprint(23)

C *   title -- of course
C *   iwlam  = 0  do not calculate wlam (i.e. <v|vlam|vp>)
C *          = 1  do calc. wlam
C *   noff   = -1 calculate <v|t|vp> from body frame t matrices
C *          =  0 calculate no matrix elements of numerical functions
C *          = +n calculate v,vp matrix elements for n numerical funcs
C *   nanl   =  0 calculate no matrix elements of analytic funcs
C *          =  1 calculate matrix elements of function given by
C *               subr. fanl
C *   inucl  =  0 homonuclear molecule
C *          =  1 heteronuclear molecule
C *   lvib   =    logical unit number of file for storing or reading
C *  vibfun  =  name of the file to be read or written on channel lvib
C *             the vibrational wavefunctions
C *  note: there is currently a problem with reading this
C *        data properly after it has been written so until
C *        it is fixed, always recalculate the wfns as needed
C *   idbug  =    logical switch - if true print lots of intermediate stuff
C *   ibasis =      '      '     - if true print the sho basis used
C *   wrtvib =      '      '     - if true write the wfn to unit lvib
C *   redvib =      '      '     - if true read the wfn from unit lvib; unformatted data
C *   sho    =      '      '     - if true use simple harm. osc. wfn.
C *   outvib =      '      '     - if true read vibrational from outer file; formatted data
C *            The FORMATS are:
C *         Nvib        Nrw      - number of vib. states, number of pair data each Nvib
C *         r(1,1)      vw(1,1)
C *         ...         ...
C *         r(Nrw,1)    vw(Nrw,1)
C *         ...         ...
C *         r(Nrw,Nvib) vw(Nrw,Nvib)
C *         ENER(1)              - Eigen Energy of each vib. state
C *         ...
C *         ENER(Nvib)

      READ (5,9000) TITLE

C     READ (5,9010) IWLAM,NOFF,NANL,INUCL,LVIB
      READ (5,*)    iwlam,noff,nanl,inucl,lvib

C     READ (5,2001) IDBUG,IBASIS,WRTVIB,REDVIB,SHO
      READ (5,2001) idbug,ibasis,wrtvib,redvib,sho,outvib

      IPRINT(1) = SHO
      WRITE (6,9031) TITLE
      IF(SHO) WRITE (6,9041)
      WRITE (6,9051) IWLAM,NOFF,NANL
      IF(INUCL.EQ.0) WRITE (6,9053)
      IF(INUCL.EQ.1) WRITE (6,9055)
      IF(REDVIB) GO TO 100

C * note: the following are not read if redvib is true
C *       bmass  - the mass in a.u. of the 'left' atom
C *       cmass  - the mass in a.u. of the 'right' atom
C *       we     - the vibrational constant (a.u.)
C *       re     - the equilibrium separation (a0)
C *       minvib - the minimum v + 1
C *       maxvib - the maximum v + 1
C *       nhermt - the order of the gauss-hermite quadrature
C *                used in evaluating integrals in the determ.
C *                of the vib wfn and in any matrix element
C *                determinations
C *       ivibtr - number of vib. channels used for homogeneous terms.
C *                if zero ivibtr=maxvib-minvib+1. only needed when more
C *                inhomogeneous terms than homo terms
C *       ihoex  - no. osc. states to expand exch kernel.
C *                use zero for normal vib state run
C *       ralfa  - a parameter near unity which must be adjusted
C *                upon the first calculations of the vib. wfns.
C *                to maximize the agreement with the experimental
C *                energy levels
C *                For H2, try 0.9
C *       rx     - another parameter near unity which must
C *                similarly be adjusted
C *                For H2, try 1.1
C *       noscil - no. of harmonic oscillators in the basis set
C *                for calc. of the vib. wfn.
C *       jmin(i),jmax(i) - are the minimum and maximum values of
C *                rotor state for each vib manifold.  the code
C *                currently only allows a single value for j so
C *                that jmin must equal jmax.  this could be
C *                fairly easily adapted.
C *       wexe,be,alfe,de - rovibrational experimental constants
C *                (along with we) for calculating exp. energies
C *                to compare with the theoretical values

C     READ (5,2031) bmass,cmass,we,re
      READ (5,*)    bmass,cmass,We,Re

C     READ (5,2021) minvib,maxvib,nhermt,ivibtr,ihoex
      READ (5,*)    minvib,maxvib,nhermt,ivibtr,ihoex

      if (ivibtr .eq. 0) ivibtr = maxvib - minvib + 1

C     READ (5,2011) ralfa,rx,noscil
      READ (5,*)    ralfa,rx,noscil

C     READ (5,2021) (jmin(i),i=1,maxvib)
C     READ (5,2021) (jmax(i),i=1,maxvib)
      READ (5,*)    (jmin(i),i=1,maxvib)
      READ (5,*)    (jmax(i),i=1,maxvib)

C     READ (5,2031) wexe,be,alfe,de
      READ (5,*)    WeXe,Be,Alfe,De

      if(maxvib-minvib+1 .gt. noscil) noscil = maxvib - minvib + 1

      if(.not.SHO) go to 100
      rx = 1.e+00
      ralfa = 1.e+00
      noscil = maxvib - minvib + 1

 100  continue

      do i = 11, 23
         iprint(i) = ibasis
      end do
      iprint(15) = idbug .and. ibasis

C * read vibrational wavefunction and quantum numbers.

      if(redvib) open(lvib, file='vibfun', status='old', form
     $     ='unformatted')

      if(redvib) READ (lvib) mvib,jrot,n,alpha,rx,re,nhermt,noscil,chi
     $     ,xpth,wpth,vecold,jmin,jmax,minvib,maxvib,energs,energy
     $     ,numax,numin,msrot,msvib,nstats,bmass,cmass,ubc,ralfa,we
     $     ,wexe,be,alfe,de

      WRITE (6,1049)

C     WRITE (6,1048) redvib,wrtvib,ibasis
      WRITE (6,1048) redvib,wrtvib,ibasis,outvib

      WRITE (6,1031) noscil,nhermt
      if(redvib) go to 35
      ubc = bmass*cmass/(bmass+cmass)
      uml = ubc/bmass
      umr = ubc/cmass

C * 1822.88625 is the "proton"/electron mass ratio, when you input
C * proton mass =/= electron mass = 1 (a.u.), you MUST use next line !
C * 1822.828 is A. Feldt's choice!

      ubc = (1822.828d0)*ubc
C     ubc = (1822.88625d0)*ubc

C *  calculate alpha

      alpha = ralfa*sqrt(we*ubc)*re

C * potint reads in data for the potential

 35   continue
      WRITE (6,1010) bmass,cmass,ubc
      WRITE (6,1041) ralfa,alpha,rx,re,we,wexe,be,alfe,de
      if (redvib) go to 37
      if ((.not. sho) .and. (.not. outvib)) call potint

      call input(mvib,jrot,n,vecold,tstore,alpha,rx,re,nhermt,noscil
     $     ,chi,xpth,wpth,xsq,xk,jmin,jmax,minvib,maxvib,energs,energy
     $     ,numax,numin,msrot,msvib,nstats)

C * write vibrational wavefunction and quantum numbers.

 37   continue
      if(wrtvib) open (lvib, file='vibfun', status='unknown', form
     $     ='unformatted')


      if(wrtvib) WRITE (lvib) mvib,jrot,n,alpha,rx,re,nhermt,noscil,chi
     $     ,xpth,wpth,vecold,jmin,jmax,minvib,maxvib,energs,energy
     $     ,numax,numin,msrot,msvib,nstats,bmass,cmass,ubc, ralfa,we
     $     ,wexe,be,alfe,de

      do i = 1, nhermt
         xpth(i) = (xpth(i)/alpha+rx)*re
      end do

 1010 format(1x,'bmass =',1pe14.7,/1x,'cmass =',1pe14.7,/3x,'ubc ='
     $     ,1pe14.7/)
 1031 format(1x,'noscil =',i5,',   nhermt =',i5/)
 1041 format(1x,'ralfa =',1pe14.7,/1x,'alpha =',1pe14.7,/4x,'rx =',1pe14
     $     .7,/4x,'re =',1pe14.7,/4x,'we =',1pe14.7,/2x,'wexe =',1pe14.7
     $     ,/4x,'be =',1pe14.7, /2x ,'alfe =',1pe14.7,/4x,'de =',1pe14.7
     $     /)
 1049 format(1x,'redvib wrtvib ibasis outvib')
 1048 format(1x,l6,l7,l7,l7/)
 1030 format(1x,i7,i7/)
 1040 format(1x,9e13.6/)
C2001 format(5l2)
 2001 format(6l2)
 2011 format(2f10.5,i5)
 2021 format(14i5)
 2031 format(4e18.8)
 9000 format(72a1)
 9010 format(5i5,a8)
 9031 format(1x,72a1)
 9041 format(//' SHO wavefunctions chosen'//)
 9051 format(10x,'options: iwlam =',i5/20x,'noff =',i5/20x,'nanl =',i5/)
 9053 format(/' homonuclear molecule'/)
 9055 format(/' heteronuclear molecule'/)
      return
      end

C *-
      double precision function fanl(r)
      implicit double precision (a-h,o-z)
      common/anlyt/parms(10)

C *  this one will evaluate q(r) by eq. 2 of truhlar or any similar
C *   polynomial
C *    i.e. f= parms(1)+parms(2)*rt+parms(3)*rt**2+parms(4)*rt**3
C *            where rt = r - req
C *  parms(5) is his req

      req=parms(5)
      rtran=r-req
      j=3
      f=parms(4)
      do i=1,3
        f=f*rtran+parms(j)
        j=j-1
        end do
      fanl=f
      return
      end

C *-
      subroutine wlam(rs,rstrt,rstop,delp,nstep,inucl,vrsep,chis,nquad
     $     ,rquad,vquad,wquad,ivm)

      implicit double precision (a-h,o-z)
      logical twice,second,test,tstchi,tstpot,tstint,tstpi
      logical outvib

      character*8 iflnm(20)
      integer ivtemp

C     dimension wlamd(10,6,6,60),rall(600),lam(4),rtst(20),ltst(4)
      dimension wlamd(10,50,50,60),rall(600),lam(4),rtst(20),ltst(4)
C     dimension wlamda(10,10),nmunit(16)
      dimension wlamda(50,50),nmunit(20)

      dimension rs(1),rstrt(1),rstop(1),delp(1),nstep(1),vrsep(1),
     $     chis(nquad,ivm),rquad(1),vquad(1),wquad(1),rnl(15),rnr(15)

      integer vibmax
      integer geomax
      parameter(geomax=15)
      double precision venuc(geomax),venq(50)

      common /splin/ bpar(11),csplin(5000)
      common /vnpar/znuc,umol,relc,lamb
      common/masses/bmass,cmass,ubc,uml,umr

      common/need/ rx,alpha,re
      common/readvib/outvib,lvib
 
      common/store/ h(20),xx(21),num(21)
      common/vvp/iv,ivp,rss(15),vrs(15),ndo

      aerr = 1.0e-06

C *  this subroutine performs the vibrational averaging of
C *  of the radial components of the electron-molecule pot.
C *  the total pot. is read in and
C *  the electron-nuclear part is subtracted and averaged 
C *      by adaptive-romberg quadrature for SHO/SPFD vib. wave
C *  and by gauss-hermite quadrature for outer reading vib.
C *  the electron-electron part is averaged by
C *      gauss-hermite quadrature.  
C *  for lamda values for which only the electr.-nucl. part 
C *  contributes, no pot. is read in and only the electr.-nucl. 
C *  pot. averaging is done.

C *  input variables:
c
C *  test     if true, then only perform certain tests of this
C *           subroutine  - see tstchi,tstpot and tstint
C *  nofr     first dimensioned size of wlamd.
C *  lelctr   maximum value of lamda for the electronic part of pot.
C *  lammax   maximum value of lambda for the vlams and wlams.
C *  nsepar   number of internuclear separations for which the
C *           vlams exist.
C *  rerr     relative error parameter for the imsl function dcadre
C *           must be small enough or a smooth wlam will not appear.
C *  rs       the values of internuclear sep. at which the vlams
C *           exist.
C *  nreg     no. of integration regions for the collision coord.
C *  rstrt    beginning value of r in each region. (in bohr)
C *  rstop    ending value of r in each region.
C *  delp     incremental value of r in each region.
C *  zl,zr    left and right nuclear charges.
C *  rnl,rnr  left and right nuclear locations.
C *  iprint=1 to print wlamda on unit 6 for each lamda.
C *  irad     frequency of such printing w.r.t. radial distance.
C *  iflnm    names of files containing the fixed nuclear r vlamdas
C *           and the name of the averaged potential.

      ixil = 0
      ixir = 1

C     READ (5,9070) test,nofr,lelctr,lammax,nsepar,iprint,irad,rerr
      READ (5,9070) test
      READ (5,*)    nofr,lelctr,lammax,nsepar,iprint,irad,rerr

C     READ (5,9080) (rnl(i),i=1,nsepar)
C     READ (5,9080) (rnr(i),i=1,nsepar)
      READ (5,*) (rnl(i),rnr(i), i = 1, nsepar)

C     READ (5,9100) nreg
C     READ (5,9110) (rstrt(i),rstop(i),delp(i),i=1,nreg)
C     READ (5,9110) zl,zr
      READ (5,*)    nreg
      READ (5,*)    (rstrt(i),rstop(i),delp(i),i=1,nreg)
      READ (5,*)    zl,zr

      write (6,*)
      WRITE (6,*) '    nofr  lelctr lammax   rerr'
      WRITE (6,9071) nofr,lelctr,lammax,rerr
      WRITE (6,9091) (i,rnl(i),rnr(i),i=1,nsepar)
      WRITE (6,9119)
      WRITE (6,9121) (rstrt(i),delp(i),rstop(i),i=1,nreg)

      nmax = 0
      do i = 1, nreg
        nstep(i) = (rstop(i)-rstrt(i))/delp(i)+5.0e-02
        nmax = nmax + nstep(i)
      end do

      nsep1 = nsepar + 1
C     READ (5,9111) (iflnm(i),i=1,nsep1)
      do i = 1, nsep1
         READ (5,9111) iflnm(i)
      enddo

      WRITE (6,9111) (iflnm(i),i=1,nsep1)
      do i = 1, nsepar
        iunit = 10 + i
        nmunit(i) = iunit
        open(iunit,file=iflnm(i),status='old', form='formatted')
      end do
      iunit1 = 10 + nsep1
      open(iunit1,file=iflnm(nsep1),status='unknown',form='unformatted')

C * the fixed r vlams as prepared by wlut are assumed to be in
C * order of increasing rs assigned to files 11,12,...

      do i = 1, nsepar
        rs(i) = rnr(i) + rnl(i)
        idev = nmunit(i)
        READ (idev,9140) rsep
        if(abs(rsep-rs(i)).lt.1.e-08) go to 30
        WRITE (6,9151) idev,rsep,rs(i)
        stop
 30     continue
      end do

      WRITE (6,9153) nsepar
      lamstp = 1
      if(inucl.eq.0) lamstp = 2
      ls = lammax/lamstp + 1

C *       tstchi - if true, calculate vib. wfns on a dense mesh and
C *                write them to file 9
C *       tstpot - if true, calculate the spline fit to the potential
C *                as a function of internuclear separation for a set
C *                of r interaction (see nofrt,noflmt) and
C *                write to file 9
C *       tstint - if true, write the integrand for <v|vlam|vp> to
C *                file 9 on a mesh ( in internuclear sep) determined
C *                roughly by stepr

      tstpot=.false.
      tstint=.false.
      tstchi=.false.
      tstpi=.false.
      if(.not.test) go to 34
      READ (5,9504) tstchi,tstpot,tstint
      WRITE (6,9508) tstchi,tstpot,tstint
      if(.not.tstpot.and..not.tstint) go to 34
      tstpi=.true.

C *    read only if tstpot or tstint are true
C *       nofrt  - no. of r interaction (rint) values to test
C *       noflmt - no. of lambdas to test
C *       stepr  - see tstint
C *       rtst   - list of rint to test
C *       ltst   - list of lambdas to test

      READ (5,9505) nofrt,noflmt,stepr
      READ (5,9506) (rtst(i),i=1,nofrt)
      READ (5,9507) (ltst(i),i=1,noflmt)
      WRITE (6,9509) (rtst(i),i=1,nofrt)
      WRITE (6,9510) (ltst(i),i=1,noflmt)
 34   continue

      if(tstchi) open(9,file='file9',status='unknown',form='formatted')
      if(tstchi) call chitst(rs(1),rs(nsepar),ivm)
      if(test.and..not.tstpi) go to 900

C *  some parameter initializations.

      lamend = lammax + 1
      del = delp(1)
      nstp = nstep(1)
      r = del
      xnrm = alpha/re
      idex = 1
      istep = 1
      nas = 0
      if(iprint.eq.1) WRITE (6,503)
      rlow = (0.8d0)*(rnr(1)+rnl(1))
      rup = (1.1d0)*(rnr(nsepar)+rnl(nsepar))
      r0 = re*rx

C *  loop for the electronic co-ordinate.
      do i = 1, nmax
         relc = r
         rmidl = r/uml
         rmidr = r/umr
         jmpl = 1
         jmpr = 1
         if(rmidl .lt. rlow) jmpl = 2
         if(rmidl .gt. rup) jmpl = 3
         if(rmidr .lt. rlow) jmpr = 2
         if(rmidr .gt. rup) jmpr = 3
         nas = nas+1
         k = mod(i,nofr)
         if(k .eq. 0) k = nofr
         rall(i) = r
         rx2 = 2.e+00*r

C *  loop for lamda.

         do l = 1, lamend, lamstp
            lind = (l-1)/lamstp + 1
            lamb = l - 1
            fctr = (-1.0d0)**lamb

C *  read the potential and subtract the nuclear part.
            if(lamb .gt. lelctr) go to 81
            do id = 1, nsepar
               idev = nmunit(id)
C              READ (idev,9160) lpot, rpot, vrsep(id)
               READ (idev,9161) lpot, rpot, vrsep(id)
               ratl = r/rnl(id)
               if(ratl .le. 1.0d0) go to 36
               vnuc = -fctr*zl*((1.0d0/ratl)**lamb)/r
               go to 37
 36            vnuc = -fctr*zl*(ratl**lamb)/rnl(id)
 37            ratr = r/rnr(id)
               if(ratr.le.1.0d0) go to 38
               vnuc = vnuc - zr*((1.0d0/ratr)**lamb)/r
               go to 39
 38            vnuc = vnuc-zr*(ratr**lamb)/rnr(id)
 39            continue
               vrs(id) = vrsep(id)-vnuc
               rss(id) = rs(id)
               if(lpot.eq.lamb.and.abs(r-rpot).lt.1.e-05) go to 40
               WRITE (6,9171) idev,lpot,rpot,lamb,r
               stop
 40            continue
            end do

C *  test section.
            if(.not.tstpi) go to 903
            do id = 1, nofrt
               if(abs(r-rtst(id)).lt.1.e-05) go to 901
            end do
            go to 80
 901        continue
            do id = 1, noflmt
               if(lamb.eq.ltst(id)) go to 903
            end do
            go to 80

C *  perform electronic potential interpolation (and some tests).
 903        continue
            nsep = nsepar
            ndo = nsepar
            second = .false.
            ier = 0

C * imsl v9
C *     call icsicu(rss,vrs,nsep,bpar,csplin,29,ier)
C *     call icsevu(rss,vrs,nsep,csplin,29,rquad,vquad,nquad,ier)
C * imsl v10
C *      call dcsdec(nsep,rss,vrs,0,0.d0,0,0.d0,bpar,csplin)
C *      do 51 ix=1,nquad
C *      vquad(ix)=dcsval(rquad(ix),nsep-1,bpar,csplin)
C *   51 continue
C * numerical recipes version:

            call spline(rss,vrs,nsep,1.d30,1.d30,csplin)
            do ix = 1, nquad
               call splint(rss,vrs,csplin,nsep,rquad(ix),vquad(ix))
            end do

            if(tstint) call inttst(r,lamb,rss(1),rss(ndo),stepr,ivm
     $           ,second)
            if(tstpot) call pottst(r,lamb,rss,vrs,ndo,second,twice)
            if(test) go to 80

C *  loops for the vibrational quantum numbers.

            if (outvib) then
               pref = re/alpha
            else
               pref = 1.0d0
            endif

 81         do iv = 1, ivm
               do ivp = iv, ivm
                  sum1 = 0.e+00
                  sum2 = 0.e+00

C *  guass-hermite integration of electronic part of pot.

                  if(lamb.gt.lelctr) go to 82
                  do j = 1, nquad
                     if (outvib) then
                        zj = (rquad(j)/re - rx)*alpha
                        temp = chis(j,iv)*chis(j,ivp)*vquad(j)*wquad(j)
                        sum1 = sum1 + temp*dexp(zj*zj)
                      else
                         sum1 = sum1 + chis(j,iv)*chis(j,ivp)*vquad(j)
     $                        *wquad(j)
                      endif
                  end do
                  sum1 = pref*sum1

C *  adaptive-romberg integration of nuclear part of pot.
C *  if outvib = .true., still do gauss-hermite integration

 82               continue
C * for outer reading vib. NOT do adaptive-romberg integration
                  if(outvib) goto 79
                  if(inucl.eq.0) go to 74
                  znuc = fctr*zl
                  umol = uml
                  if(jmpl.eq.2) go to 72
                  if(jmpl.eq.3) go to 73
                  sum2 = sum2 + quanc8(ixil,rlow,rmidl,aerr,rerr,error
     $                 ,nofun,flag)
                  sum2 = sum2 + quanc8(ixir,rmidl,rup,aerr,rerr,error
     $                 ,nofun,flag)
                  go to 74
 72               sum2 = sum2 + quanc8(ixir,rlow,r0,aerr,rerr,error
     $                 ,nofun,flag)
                  sum2 = sum2 + quanc8(ixir,r0,rup,aerr,rerr,error,nofun
     $                 ,flag)
                  go to 74
 73               sum2 = sum2 + quanc8(ixil,rlow,r0,aerr,rerr,error
     $                 ,nofun,flag)
                  sum2 = sum2 + quanc8(ixil,r0,rup,aerr,rerr,error,nofun
     $                 ,flag)

 74               znuc = zr
                  umol = umr
                  if(jmpr.eq.2) go to 76
                  if(jmpr.eq.3) go to 77
                  sum2 = sum2 + quanc8(ixil,rlow,rmidr,aerr,rerr,error
     $                 ,nofun,flag)
                  sum2 = sum2 + quanc8(ixir,rmidr,rup,aerr,rerr,error
     $                 ,nofun,flag)
                  go to 78
 76               sum2 = sum2 + quanc8(ixir,rlow,r0,aerr,rerr,error
     $                 ,nofun,flag)
                  sum2 = sum2 + quanc8(ixir,r0,rup,aerr,rerr,error,nofun
     $                 ,flag)
                  go to 78
 77               sum2 = sum2 + quanc8(ixil,rlow,r0,aerr,rerr,error
     $                 ,nofun,flag)
                  sum2 = sum2 + quanc8(ixil,r0,rup,aerr,rerr,error,nofun
     $                 ,flag)

C *  final matrix element assignment.

 78               if(inucl.eq.0) sum2 = 2.0d0*sum2
C                 wlamd(k,iv,ivp,lind) = sum2*xnrm+sum1
C                 wlamd(k,ivp,iv,lind) = sum2*xnrm+sum1

 79               continue
                  if (outvib) then
C * do ele-nucl part coupling by gauss-hermite integration
                     nsep = nsepar
                     do id = 1, nsepar
                        venuc(id) = potenu(r,rnl(id),rnr(id),zl,zr,lamb)
                        rss(id) = rs(id)
                     enddo
                     call spline(rss,venuc,nsep,1.d30,1.d30,csplin)
                     do ix = 1, nquad
                        call splint(rss,venuc,csplin,nsep,rquad(ix)
     $                       ,venq(ix))
                     enddo
                     pref = re/alpha
                     sum2 = 0.0d0
                     do j = 1, nquad
                        zj = (rquad(j)/re - rx)*alpha
                        temp = chis(j,iv)*chis(j,ivp)*venq(j)*wquad(j)
                        sum2 = sum2 + temp*dexp(zj*zj)
                     enddo
                     sum2 = sum2*pref
C * sum up ele-ele and ele-nuc
                     temp = sum1 + sum2
                  else
                     temp = sum2*xnrm+sum1
                  endif

                  wlamd(k,iv,ivp,lind) = temp
                  wlamd(k,ivp,iv,lind) = temp

                  if(iprint.ne.1) go to 70
C                 wlamda(iv,ivp) = sum2*xnrm+sum1
C                 wlamda(ivp,iv) = sum2*xnrm+sum1
                  wlamda(iv,ivp) = temp
                  wlamda(ivp,iv) = temp
 70               continue
C * end ivp loop
               enddo
C * end iv loop
            enddo

C *  to print wlams on unit 6 if iprint=1

            if(iprint.eq.0) go to 80
C           if(nas.ne.irad) go to 80
            nas = 0
            ivtemp = ivm
            if (ivm .gt. 5) ivtemp = 5
            do iv = 1,ivm
               ivmin1 = iv-1
c$$$               WRITE (6,502) rall(i),lamb,ivmin1,(wlamda(iv,ivp),ivp = 1
c$$$     $              ,ivm)
C *               WRITE (6,502) rall(i),lamb,iv,(wlamda(iv,ivp),ivp =1,ivm)
               WRITE (6,502) rall(i),lamb,iv,(wlamda(iv,ivp),ivp =1
     $              ,ivtemp)

            end do
 80         continue
C * end l loop
         enddo

         if(k.eq.nofr.and..not.test) WRITE (iunit1) ((((wlamd(jr,iv,jv,l
     $        ),jr=1,nofr),iv=1,ivm),jv=1,ivm),l=1,ls)
         do jr = 1, nofr
            WRITE (50,5001)(wlamd(jr,jv,jv,1),jv=1,ivm)
         end do
 5001    format(4(1pe12.4))

c     if(k.ne.nofr) WRITE (iunit1) wlamd
         if(idex.ne.nstp) go to 90
         istep = istep + 1
         del = delp(istep)
         nstp = nstep(istep)
         idex = 0
 90      r = r + del
         idex = idex + 1
C * end i loop
      enddo

      do i = 1, nsepar
         iunit = nmunit(i)
         close(iunit)
      end do
      close(iunit1)

      if(test) go to 900
      WRITE (6,3001)
      WRITE (6,3003) lamend, lamstp
      do l = 1, lamend, lamstp
         lind = (l-1)/lamstp + 1
         lam(lind) = l - 1
      end do
      
      WRITE (6,9173)
 900  return

 502  format(f10.5,2i4,15e18.6)
 503  format(/6x,'r',3x,'lamda',3x,'vi',6x,'wlamda(vi,vf),vf=1,ivm',/)
 3001 format(' main loop finished')
 3003 format(' all of wlamd written -- lamend=',i5,' lamstp=',i5)
 9070 format(l5,6i5,f15.10)
 9071 format(3i7,f15.10)
 9080 format(7f10.5)
 9091 format(//10x,'rnucl(',i2,') = ',2f10.5)
 9100 format(5i5)
 9110 format(3f10.5)
 9111 format(10a8)
 9119 format(//' rmesh:   start     step     stop'/)
 9121 format(7x,3f10.5)
 9140 format(10x,f10.5)
 9151 format(//' tsk-tsk-tsk ... file ',i2,' read rsep=',f10.5,' while t
     1he program had rsep=',f10.5,' from file 5'/1x,90('*'))
 9153 format(/' rseps in ',i2,' files 11 on are okay'/)
 9160 format(i5,f10.5,e23.16)
 9161 format(i5,f10.5,e24.16)
 9171 format(///' come now!!  error in file ',i2,' which had lambda=',i5
     1,' and r=',f10.5/15x,'whereas file  5     had lambda=',i5,' and r=
     2',f10.5/1x,90('*'))
 9173 format(/' wlamd calculation successful'/)
 9504 format(3l5)
 9505 format(2i5,f10.5)
 9506 format(7f10.5)
 9507 format(10i5)
 9508 format(/2x,'tstchi tstpot tstint',/2x,l5,l7,l7)
 9509 format(/2x,' rtst=',7f10.5,(7x,7f10.5))
 9510 format(/2x,' ltst=',10i5)
      end

C *-
C * this is a modification of the numfun routine in the wlam codes.
C *  this routine give the vib matrix elements for the exchange
C *  kernel as used by the lavib code.

      subroutine numfunex(ivm,nquad,fnum,chis,rquad,vquad,wquad,ivibtr
     $     ,ihoex)

      implicit double precision (a-h,o-z)
      integer pwav,vibs,points,seps
      parameter (maxqd=20,maxho=20,pwav=11,vibs=25,points=350,seps=14)
C * maxqd is max quadrature, maxho=max harmonic osc.
      dimension fin(pwav,pwav,points,seps), fout(pwav*vibs,pwav*vibs
     $     ,points)
      double precision finone(pwav,pwav,points)
      dimension chis(nquad,ivm),fnum(1),rquad(1),vquad(1),wquad(1)
     $     ,nmunit(16),rss(15)
      dimension ho(maxqd,maxho),hox(maxho)
      common/need/rx,alpha,re
 
C ***** build harmonic oscillators if needed
      if (ihoex.gt.0) then
         do ix=1,nquad
            zold=alpha*(rquad(ix)/re-rx)
            call dhep(hox,zold,ihoex)
            do iy=1,ihoex
               ho(ix,iy)=hox(iy)
            end do
         end do
      endif
*****
c
C *  this subroutine vibrationally averages functions given
C *  at a set of internuclear separations.  the values can
C *  all be in one file or in the case of radial orbitals
C *  can be in nsepar files in unformatted form.
c
C *  input variables:
c
C *       nsepar     number of internuclear seps at which the numerical
C *                  double precision function exists.
C *  the next two are read in by exinout now, but documented here
C *  nfile      if 1, file has values to be averaged which are
C *             read into fnum, or if .gt. 1, equals the number
C *             (nsepar) of files containing values to averaged.
C *  *** iprint,ivibopt,norbi and norbo are deadwood when used for
C *  *** the exchange kernel, but are kept for input consistency
C *  *** with the case of moment matrix elements - see subr. numfun
C *  iprint     if nonzero, results printed on unit 6.
C *  ivibopt    if nonzero, do overlap of function and vib state
C *             if zero, do vib. avg. <v|f(r)|vp>
C *  norbi      no. of orbitals input.
C *  norbo      no. of orbitals output.
C *  nlexi      no. of l values input. - must match nexdim in vibker input
C *  nlexo      no. of l values output. - wlam4 must match this and ivm
c

C     READ (5,9180) nsepar,nfile,iprint,ivibopt
      READ (5,*)    nsepar,nfile,iprint,ivibopt
C     READ (5,9180) norbi,norbo,nlexi,nlexo
      READ (5,*)    norbi,norbo,nlexi,nlexo
      WRITE (6,*) ' nsepar nfile iprint'
      WRITE (6,9400) nsepar,nfile,iprint
      WRITE (6,*) '  norbi norbo nlexi nlexo'
      WRITE (6,9400) norbi,norbo,nlexi,nlexo

      READ (5,*)     npts
      WRITE (6,*) '  npts'
      WRITE (6,9400) npts

C *  the remainder of this has been modified by A. Feldt to
C *  allow the unformatted I/O size to be dynamic (up to the
C *  compiled in maxima above)

C *  do some sanity checks here...

      if (nlexi.gt.pwav .or. nlexo.gt.pwav) then
         WRITE (6,*) 'You must recompile with pwav > nlexi and'
         WRITE (6,*) 'pwav > nlexo'
         stop
      endif      
      if (nlexi.lt.nlexo) then
         WRITE (6,*) 'nlexi cannot be less than nlexo'
         stop
      endif
      if (ivm.gt.vibs) then
         WRITE (6,*) 'You must recompile with vibs > ivm'
         stop
      endif
      if (nsepar.ne.nfile) then
        WRITE (6,*) 'nsepar must equal nfile - try again'
        stop
      endif

C     npts=points

      call exinout(ivm,nquad,fnum,chis,rquad,vquad,wquad,ivibtr,ihoex
     $     ,nsepar,nlexi,nlexo,npts,fin,fout,finone)

      return
 9180 format(5i5)
 9400 format(6i7)
      end

c  ***
C *  This routine added by A. Feldt to modularize the unformatted I/O
C *  The code is essentially unchanged from the original, warts and all
C *
C *  Hao Feng modified to cope with other vibrational wavefunctions
c  ***
      subroutine exinout(ivm,nquad,fnum,chis,rquad,vquad,wquad,ivibtr
     $     ,ihoex,nsepar,nlexi,nlexo,npts,fin,fout,finone)

      implicit double precision (a-h,o-z)

      dimension fin(nlexi,nlexi,npts,nsepar)
      dimension finone(nlexi,nlexi,npts)
      dimension fout(nlexo*ivm,nlexo*ivm,npts)
      dimension chis(nquad,ivm),fnum(1),rquad(1),vquad(1),wquad(1)
     $     ,nmunit(16),rss(15)
      character*8 iflnm(20)

      logical outvib

      common /splin/ bpar(11),csplin(5000)

      common/need/rx,alpha,re
      common/readvib/outvib,lvib

      zero = 0.0e+00
      if (outvib) then
         pref = re/alpha
      else
         pref = 1.0d0
      endif

C *       rss,fnum   the vectors holding the internuc. sep. values and
C *                  the values of the numerical function (e.g. the
C *                  polarizabilities given by poll and wolniewicz).

C     READ (5,9080) (rss(j),j=1,nsepar)
      READ (5,*)    (rss(j),j=1,nsepar)
      WRITE (6,*) 'the nuclear separations are:'
      WRITE (6,9080) (rss(i),i=1,nsepar)

C *  iflnm      names of the input and output files (.le. 6 char.).
C *             there should be nsepar + 1 names.

      nfil1 = nsepar + 1
C     READ (5,9111) (iflnm(i),i=1,nfil1)
      do i = 1, nfil1
         READ (5,9111) iflnm(i)
      enddo

      WRITE (6,9111) (iflnm(i),i=1,nfil1)
      do i = 1, nsepar
         iunit = 10 + i
         nmunit(i) = iunit
         open(iunit,file=iflnm(i),status='old', form='unformatted')
      end do

      iunit1 = 10 + nfil1
      open(iunit1,file=iflnm(nfil1),status='unknown',form='unformatted')

      do ir = 1, npts
         do  i = 1, nsepar
            idev = nmunit(i)
            READ (idev) finone
            do j = 1, nlexi
               do l = 1, nlexi
                  do m = 1, npts
                     fin(j,l,m,i) = finone(j,l,m)
                  end do
               end do
            end do
         end do
         do i = 1, nlexo
            do k = 1, nlexo
               do l = 1, npts
                  do j = 1, nsepar
                     fnum(j)=fin(i,k,l,j)
                  end do
                  ier = 0
C * numerical recipes version:
                  call spline(rss,fnum,nsepar,1.d30,1.d30,csplin)
                  do ix = 1, nquad
                     call splint(rss,fnum,csplin,nsepar,rquad(ix)
     $                    ,vquad(ix))
                  end do

                  if (ihoex .eq. 0) then
                     do iv = 1, ivm
                        ich = (iv-1)*nlexo + i
                        do ivp = 1, ivm
                           ichp = (ivp-1)*nlexo + k
                           if (iv.gt.ivibtr.and.ivp.gt.ivibtr) goto 220
                           sum = zero
                           do j = 1, nquad
                              if (outvib) then
                                 zj = (rquad(j)/re - rx)*alpha
                                 temp = chis(j,iv)*chis(j,ivp)*vquad(j)
     $                                *wquad(j)
                                 sum = sum + temp*dexp(zj*zj)
                              else
                                 sum=sum+chis(j,iv)*chis(j,ivp)*vquad(j)
     $                                *wquad(j)
                              endif
                           end do
                           sum = pref*sum
                           fout(ich,ichp,l)=sum
 220                       continue
                        enddo
                     enddo
                     iout = ivm
                  endif
C * end l loop
               enddo
C * end k loop
            enddo
C * end i loop
         enddo

         WRITE (iunit1) fout

C * end ir loop
      enddo

c
      do i = 1, nsepar
         iunit = nmunit(i)
         close(iunit)
      end do
      close(iunit1)
c
      return
 9111 format(10a8)
 9080 format(7f10.5)
      end
********************************************************************
      subroutine tmatav(ivm,nquad,rs,fnum,flin,chis,rquad,vquad,
     1  wquad)
      implicit double precision (a-h,o-z)
        logical kmats
      character*8  iflnm(20)
      dimension chis(nquad,ivm),rs(1),fnum(1),flin(20,11),rquad(1),
     1  vquad(1),wquad(1)
      dimension tmat(20,20,3,3,2)
      dimension nmunit(16)
      common /splin/ bpar(11),csplin(5000)
       zero=0.0e+00
c
C *       nsepar - the number of internuclear seps. at which the body
C *                frame t matrices exist
C *       rs     - the values of these internuc. seps.
C *       kmats  - if true, use k matrices instead... a stupid idea!
c
      READ (5,9180) nsepar
      WRITE (6,*) nsepar,' separations'
      READ (5,9080) (rs(j),j=1,nsepar)
      WRITE (6,9080) (rs(j),j=1,nsepar)
      READ (5,9200) kmats
      WRITE (6,*) ' kmats ? ',kmats
      READ (5,9111) (iflnm(i),i=1,nsepar+1)
      WRITE (6,*) 'file names -'
      WRITE (6,9111) (iflnm(i),i=1,nsepar+1)
      ind=(ivm*ivm-ivm)/2+ivm

      do i=1,nsepar
        iunit=10+i
        nmunit(i)=iunit
        open(iunit,file=iflnm(i),status='old',form='formatted')
        end do
      open(nsepar+11,file=iflnm(nsepar+1),status='unknown',
     +         form='formatted')

c
C *       the t matrices are assumed to be assigned in order of increasing
C *  orders to files 11,12,13,...
c
   10 continue
      do i=1,nsepar
        iunit=10+i
        READ (iunit,*,end=9999) m,li,lf,nchan,e
        end do
      if(nchan.gt.20) stop 'nchan>20'
      WRITE (nsepar+11,9300) m,li,lf,ivm*nchan,e

      nri=2
      if(kmats) nri=1

      do iri=1,nri
      do kk=1,nchan
      ii=0
      do i=1,nsepar
        iunit=10+i
        ii=ii+1
        READ (iunit,9310) (flin(j,ii),j=1,nchan)
        end do
      do jj=1,nchan
      do i=1,nsepar
        fnum(i)=flin(jj,i)
        end do
      ier=0
C * imsl v9
C *     call icsicu(rs,fnum,nsepar,bpar,csplin,29,ier)
C *     call icsevu(rs,fnum,nsepar,csplin,29,rquad,vquad,nquad,ier)
C * imsl v10
C *      call dcsdec(nsepar,rs,fnum,0,0.d0,0,0.d0,bpar,csplin)
C *      do 51 ix=1,nquad
C *      vquad(ix)=dcsval(rquad(ix),nsepar-1,bpar,csplin)
C *   51 continue
C * numerical recipes version:
      call spline(rs,fnum,nsepar,1.d30,1.d30,csplin)
      do ix=1,nquad
         call splint(rs,fnum,csplin,nsepar,rquad(ix),vquad(ix))
      end do
      ind=0
      do iv=1,ivm
        do ivp=iv,ivm
          ind=ind+1
          sum=zero
          do j=1,nquad
            sum=sum+chis(j,iv)*chis(j,ivp)*vquad(j)*wquad(j)
            end do
          tmat(jj,kk,ivp,iv,iri)=sum
          tmat(jj,kk,iv,ivp,iri)=sum
          enddo
        enddo
C * end jj loop
      enddo
C * end kk loop
      enddo
C * end iri loop
      enddo

*********** print avg. t-matrix and calc cross-section
      do iri=1,nri
        do iv=1,ivm
          do i=1,nchan
            WRITE (nsepar+11,9310)
     x     ((tmat(j,i,iv,ivp,iri),j=1,nchan),ivp=1,ivm)
            end do
          end do
        end do

      do iv=1,ivm
        do ivp=1,ivm
          cfm=0.d0
          do i=1,nchan
            eksq=e
            coef=3.14159265e+00/eksq
            do k=1,nchan
              treal=tmat(i,k,iv,ivp,1)
              timag=tmat(i,k,iv,ivp,2)
              tabsq=treal*treal+timag*timag
              temp=coef*tabsq
              cfm=cfm+temp
              end do
            end do
          WRITE (6,475) iv-1,ivp-1,cfm,e
  475 format(' xsec',2i4,2x,e14.7,f10.6)
          end do
        end do
      go to 10
c
 9999 WRITE (6,*) ' end of run'
      return

 9180 format(5i5)
 9080 format(7f10.5)
 9200 format(l5)
 9300 format(4i5,f15.8)
 9111 format(10a8)
 9310 format(4(1pe23.16))
      end
********************************************************************
      subroutine anlfun(chis,rquad,wquad,nquad,ivm,minvib,fmat)
      implicit double precision (a-h,o-z)
      dimension chis(nquad,ivm),fmat(ivm,ivm),rquad(1),wquad(1)
      common /splin/ bpar(11),csplin(5000)
      common/anlyt/parms(10)
      zero=0.0e+00
c
C *       npar   - number of parameters to pass to the function fanl
C *       parms  - the values of these parameters
c
      READ (5,9210) npar
      if(npar.le.0) go to 310
      READ (5,9220) (parms(i),i=1,npar)
      WRITE (6,9221)
      WRITE (6,9223) (parms(i),i=1,npar)
  310 do iv=1,ivm
        do ivp=iv,ivm
          sum=zero
          do j=1,nquad
            sum=sum+chis(j,iv)*chis(j,ivp)*fanl(rquad(j))*wquad(j)
            end do
          fmat(iv,ivp)=sum
          end do
        end do
      WRITE (6,9231)
      do iv=1,ivm
        ivr=iv-1+minvib-1
        do ivp=iv,ivm
          ivpr=ivp-1+minvib-1
          WRITE (6,9201) ivr,ivpr,fmat(iv,ivp)
          end do
        end do
      return
 9201 format(1x,2i5,2x,1pd14.7)
 9210 format(5i5)
 9220 format(4e20.12)
 9221 format(//' values of parms for fanl'/)
 9223 format(5x,e20.12)
 9231 format(///' vibrational matrix elements for analytic function'///
     15x,'v',3x,'vp',9x,'<f>')
      end

C *-
      subroutine input(mvib,mrot,negvl,vect,haml,alpha,rx,re,ngaush
     $     ,noscl,chi,xpth,wpth,e,h,jmin,jmax,minvib,maxvib,energs
     $     ,energy,numax,numin,msrot,msvib,nstats)

      implicit double precision (a-h,o-z)

C * this routine determine the quantum numbers for a particular
C * symmetry and the unweighted vibrational wavefunctions
C * for that symmetry.
C * written by g. a. parker. modified by jim lill.  and probably destroyed
C * by a. feldt

      logical iprint

C *  iprint(11): prints quantum numbers for external region

      dimension vect(20,20),haml(1),mvib(1),mrot(1),chi(1),xpth(1)
     $     ,wpth(1),e(1),h(1),jmin(1),jmax(1)

      dimension energs(1),energy(1),numax(1),numin(1),msrot(1),msvib(1)

      common /prntop/ iprint(23)

C * branch is made.

      mprnt = minvib - 1
      nprnt = maxvib - 1
      write (6,380) mprnt,nprnt
      ntot = 0
      write (6,400) (jmax(i),i=minvib,maxvib)
      write (6,410) (jmin(i),i=minvib,maxvib)
      jp1 = 0

C * determine the maximum rotational state.

      do i = minvib, maxvib
         if (jmax(i) .gt. jp1) jp1 = jmax(i)
      end do
      jp1 = jp1 + 1
      negvl = 1

C * loop over the vibrations.

      do k = minvib, maxvib
         jhigh = jmax(k)
         jlow = jmin(k)

C * determine quantum numbers for each vibrational state

         call  figure(nrotk,jlow,jhigh,mrot,ntot)
         do l = 1, nrotk
            mvib(l+ntot) = k - 1
         end do
         ntot = ntot + nrotk
      end do
      negvl = ntot
      if(negvl .gt. nstats) go to 370

C * check to see if negvl is bigger than the dimensions
C * statements allow.

      do i = 1, negvl
         msrot(i) = mrot(i)
         msvib(i) = mvib(i)
      end do
      if (.not.iprint(11)) go to 50
      write (6,420)
      write (6,430) negvl

C * all quantum numbers are stored in msrot,msvib

      write (6,450)
      write (6,440) (msvib(i),msrot(i),i=1,negvl)
 50   continue
      nsgvl = negvl

C * determine the maximum little j used.
C * determine the minimum little j used

      jlmax = 0
      jlmin = 10000
      do jl = minvib, maxvib
        if (jmax(jl).gt.jlmax) jlmax = jmax(jl)
        if (jmin(jl).lt.jlmin) jlmin = jmin(jl)
      end do

C * store the maximum vibration for each little j in numax
C * store the minimum vibration for each little j in numin

      jlmax = jlmax + 1
      jlmin = jlmin + 1
      do jlp1 = jlmin, jlmax
         jl = jlp1 - 1
         numax(jlp1) = 0
         numin(jlp1) = 10000
         do i = 1, negvl
            if (jl.ne.msrot(i)) go to 70
            if (msvib(i).gt.numax(jlp1)) numax(jlp1) = msvib(i)
            if (msvib(i).lt.numin(jlp1)) numin(jlp1) = msvib(i)
 70         continue
         end do
      end do

C *-
C * Calculate wavefunctions and store them in chi
C *         xpth & wpth: store the nodes & weights of hermite polynomials.
C *                      when jmin(i), jmax(i) are inputted as 0 for all
C *                      rotational states, jlmin=jlmax=1 .
C * nulow, numin(jlp1) : quantum number of the lowest state;
C * nuhigh,numax(jlp1) : quantum number of the highest state.
C *                ivm : total number of vibrational states.  ivm =< nudim.
C *              nudim : total number of states.   nudim=ivm when jlmin=jlmax=1.
C * when vibfun returns, chi contains the unweighted eigenfunctions. 
C *-/
      index = 1
      do jlp1 = jlmin, jlmax
         jl = jlp1 - 1
         nuhigh = numax(jlp1) + 1
         nulow = numin(jlp1) + 1
         nudim = nuhigh-nulow + 1
         mgh = ngaush*nudim

         call  vibfun(noscl,jl,ngaush,haml,vect,chi(index),nulow,nuhigh
     $        ,nudim,energs,alpha,rx,xpth,wpth,re,e,h)

C * store internal vibrational-rotational energies in energy.

         do nvp1 = 1, nudim
            nu = nvp1 + nulow - 2
            do i = 1, negvl
               if (jl.ne.msrot(i)) go to 90
               if (nu.ne.msvib(i)) go to 90
               energy(i) = energs(nvp1)
 90            continue
            end do
         end do
         index = index + mgh
      end do
      return
 370  continue
      WRITE (6,510) negvl
      stop

 380  format(1x,'vibrational states are from nu =',i3,1x,'to nu =',i3)
  400 format(/1x,'jmax(i) =',20i5)
  410 format(1x,'jmin(i) =',20i5)
  420 format(/1x,'quantum numbers.')
  430 format(3x,'n = ',i5)
  440 format(1x,2i10)
  450 format(/6x,'msvib',5x,'msrot')
  510 format(1x,'***** error nstats is too small   negvl=',i10)
      end

C *-
      subroutine dhep(h,x,n)
      implicit double precision (a-h,o-z)

C * written by g. a. parker
C * recursion to calculate the hermite polynomials h(i+1) for
C * i=0 to i=n at x.

      dimension h(n+1)
      data sqr2/1.4142135623730950488016887d0/
      data pi4i/7.511255444649424828587030d-1/
      dx=x
      sqr2dx=sqr2*dx
      hnm1=pi4i
      hn=sqr2dx*pi4i
      h(1)=hnm1
      h(2)=hn
      if(n.ge.2)go to 10
      return

C * start recursion

 10   do i=2,n
         hnp1=(sqr2dx*hn-sqrt(float(i-1))*hnm1)/sqrt(float(i))
         hnm1=hn
         h(i+1)=hnp1
         hn=hnp1
      end do
      return
      end
C *-/

C *-
      subroutine figure(nrotk,jlow,jhigh,mrot,ntot)
      implicit double precision (a-h,o-z)

C * determines the body-fixed quantum numbers.

      dimension  mrot(1)
      jhigh1 = jhigh + 1
      jlow1 = jlow + 1
      icount = 1 + ntot
      do i = jlow1, jhigh1
        jlitl = i - 1
        i21 = 1
        do k = 1, i21
          mrot(icount) = jlitl
          icount = icount + 1
          end do
        end do
      nrotk = icount - 1 - ntot
      return
      end

C *-
      subroutine hermit(npt,xpt,wht)
      implicit double precision (a-h,o-z)

C * purpose
C *   supplies the gauss hermite abscissas and weights

C *   description of parameters
C *      npt-number of points desired
C *      xpt-resultant array containing the gauss hermite abscissas.
C *      wht-resultant array containing the gauss hermite weights
C *   someone removed use of kase from this routine, it now acts
C *      only as if kase=0 is true - found by A. Feldt, 2-17-04
C *      kase=0 if integrating exp(-x**2)*f(x)*dx
C *      kase=1 if integrating f(x)*dx

C *   subroutine and functions used
C *      none

C *   written by g. a. parker , all zeros and weights were calculated
C *   on the ibm 7030 stretch using double precision (128 bits).

      dimension xpt(1),wht(1)
c     dimension x(649),w(649)
      dimension x(700),w(700)
      dimension x2(1),w2(1)
      dimension x3(2),w3(2)
      dimension x4(2),w4(2)
      dimension x5(3),w5(3)
      dimension x6(3),w6(3)
      dimension x7(4),w7(4)
      dimension x8(4),w8(4)
      dimension x9(5),w9(5)
      dimension x10(5),w10(5)
      dimension x11(6),w11(6)
      dimension x12(6),w12(6)
      dimension x13(7),w13(7)
      dimension x14(7),w14(7)
      dimension x15(8),w15(8)
      dimension x16(8),w16(8)
      dimension x17(9),w17(9)
      dimension x18(9),w18(9)
      dimension x19(10),w19(10)
      dimension x20(10),w20(10)
      DIMENSION X21(11),W21(11)
      DIMENSION X22(11),W22(11)
      DIMENSION X23(12),W23(12)
      DIMENSION X24(12),W24(12)
      DIMENSION X25(13),W25(13)
      DIMENSION X26(13),W26(13)
      DIMENSION X27(14),W27(14)
      DIMENSION X28(14),W28(14)
      DIMENSION X29(15),W29(15)
      DIMENSION X30(15),W30(15)
      DIMENSION X31(16),W31(16)
      DIMENSION X32(16),W32(16)
      DIMENSION X33(17),W33(17)
      DIMENSION X34(17),W34(17)
      DIMENSION X35(18),W35(18)
      DIMENSION X36(18),W36(18)
      DIMENSION X37(19),W37(19)
      DIMENSION X38(19),W38(19)
      DIMENSION X39(20),W39(20)
      DIMENSION X40(20),W40(20)
      DIMENSION X41(21),W41(21)
      DIMENSION X42(21),W42(21)
      DIMENSION X43(22),W43(22)
      DIMENSION X44(22),W44(22)
      DIMENSION X45(23),W45(23)
      DIMENSION X46(23),W46(23)
      DIMENSION X47(24),W47(24)
      DIMENSION X48(24),W48(24)
      DIMENSION X49(25),W49(25)
      DIMENSION X50(25),W50(25)

      equivalence(x2(1),x(  1)),(w2(1),w(  1))
      equivalence(x3(1),x(  2)),(w3(1),w(  2))
      equivalence(x4(1),x(  4)),(w4(1),w(  4))
      equivalence(x5(1),x(  6)),(w5(1),w(  6))
      equivalence(x6(1),x(  9)),(w6(1),w(  9))
      equivalence(x7(1),x( 12)),(w7(1),w( 12))
      equivalence(x8(1),x( 16)),(w8(1),w( 16))
      equivalence(x9(1),x( 20)),(w9(1),w( 20))
      equivalence(x10(1),x( 25)),(w10(1),w( 25))
      equivalence(x11(1),x( 30)),(w11(1),w( 30))
      equivalence(x12(1),x( 36)),(w12(1),w( 36))
      equivalence(x13(1),x( 42)),(w13(1),w( 42))
      equivalence(x14(1),x( 49)),(w14(1),w( 49))
      equivalence(x15(1),x( 56)),(w15(1),w( 56))
      equivalence(x16(1),x( 64)),(w16(1),w( 64))
      equivalence(x17(1),x( 72)),(w17(1),w( 72))
      equivalence(x18(1),x( 81)),(w18(1),w( 81))
      equivalence(x19(1),x( 90)),(w19(1),w( 90))
      equivalence(x20(1),x(100)),(w20(1),w(100))
      EQUIVALENCE(X21(1),X(110)),(W21(1),W(110))
      EQUIVALENCE(X22(1),X(121)),(W22(1),W(121))
      EQUIVALENCE(X23(1),X(132)),(W23(1),W(132))
      EQUIVALENCE(X24(1),X(144)),(W24(1),W(144))
      EQUIVALENCE(X25(1),X(156)),(W25(1),W(156))
      EQUIVALENCE(X26(1),X(169)),(W26(1),W(169))
      EQUIVALENCE(X27(1),X(182)),(W27(1),W(182))
      EQUIVALENCE(X28(1),X(196)),(W28(1),W(196))
      EQUIVALENCE(X29(1),X(210)),(W29(1),W(210))
      EQUIVALENCE(X30(1),X(225)),(W30(1),W(225))
      EQUIVALENCE(X31(1),X(240)),(W31(1),W(240))
      EQUIVALENCE(X32(1),X(256)),(W32(1),W(256))
      EQUIVALENCE(X33(1),X(272)),(W33(1),W(272))
      EQUIVALENCE(X34(1),X(289)),(W34(1),W(289))
      EQUIVALENCE(X35(1),X(306)),(W35(1),W(306))
      EQUIVALENCE(X36(1),X(324)),(W36(1),W(324))
      EQUIVALENCE(X37(1),X(342)),(W37(1),W(342))
      EQUIVALENCE(X38(1),X(361)),(W38(1),W(361))
      EQUIVALENCE(X39(1),X(380)),(W39(1),W(380))
      EQUIVALENCE(X40(1),X(400)),(W40(1),W(400))
      EQUIVALENCE(X41(1),X(420)),(W41(1),W(420))
      EQUIVALENCE(X42(1),X(441)),(W42(1),W(441))
      EQUIVALENCE(X43(1),X(462)),(W43(1),W(462))
      EQUIVALENCE(X44(1),X(484)),(W44(1),W(484))
      EQUIVALENCE(X45(1),X(506)),(W45(1),W(506))
      EQUIVALENCE(X46(1),X(529)),(W46(1),W(529))
      EQUIVALENCE(X47(1),X(552)),(W47(1),W(552))
      EQUIVALENCE(X48(1),X(576)),(W48(1),W(576))
      EQUIVALENCE(X49(1),X(600)),(W49(1),W(600))
      EQUIVALENCE(X50(1),X(625)),(W50(1),W(625))


C * abscissas and weights for the gauss hermite integration order=2

       data x2(1)/.70710678118654d+00/,w2(1)/.88622692545275d+00/

C * abscissas and weights for the gauss hermite integration order=3

       data x3(1)/.12247448713916d+01/,w3(1)/.29540897515092d-00/
       data x3(2)/.00000000000000d 00/,w3(2)/.11816359006037d+01/

C * abscissas and weights for the gauss hermite integration order=4

       data x4(1)/.16506801238858d+01/,w4(1)/.81312835447243d-01/
       data x4(2)/.52464762327529d+00/,w4(2)/.80491409000550d+00/

C * abscissas and weights for the gauss hermite integration order=5

       data x5(1)/.20201828704561d+01/,w5(1)/.19953242059046d-01/
       data x5(2)/.95857246461381d+00/,w5(2)/.39361932315223d-00/
       data x5(3)/.00000000000000d 00/,w5(3)/.94530872048294d+00/

C * abscissas and weights for the gauss hermite integration order=6

       data x6(1)/.23506049736745d+01/,w6(1)/.45300099055088d-02/
       data x6(2)/.13358490740137d+01/,w6(2)/.15706732032285d-00/
       data x6(3)/.43607741192761d-00/,w6(3)/.72462959522438d+00/

C * abscissas and weights for the gauss hermite integration order=7

       data x7(1)/.26519613568352d+01/,w7(1)/.97178124509949d-03/
       data x7(2)/.16735516287675d+01/,w7(2)/.54515582819126d-01/
       data x7(3)/.81628788285896d+00/,w7(3)/.42560725261012d-00/
       data x7(4)/.00000000000000d 00/,w7(4)/.81026461755680d+00/

C * abscissas and weights for the gauss hermite integration order=8

       data x8(1)/.29306374202572d+01/,w8(1)/.19960407221136d-03/
       data x8(2)/.19816567566958d+01/,w8(2)/.17077983007413d-01/
       data x8(3)/.11571937124468d+01/,w8(3)/.20780232581489d-00/
       data x8(4)/.38118699020732d-00/,w8(4)/.66114701255822d+00/

C * abscissas and weights for the gauss hermite integration order=9

       data x9(1)/.31909932017815d+01/,w9(1)/.39606977263262d-04/
       data x9(2)/.22665805845318d+01/,w9(2)/.49436242755368d-02/
       data x9(3)/.14685532892167d+01/,w9(3)/.88474527394374d-01/
       data x9(4)/.72355101875284d+00/,w9(4)/.43265155900254d-00/
       data x9(5)/.00000000000000d 00/,w9(5)/.72023521560605d+00/

C * abscissas and weights for the gauss hermite integration order=10

       data x10( 1)/.34361591188377d+01/,w10( 1)/.76404328552324d-05/
       data x10( 2)/.25327316742328d+01/,w10( 2)/.13436457467812d-02/
       data x10( 3)/.17566836492999d+01/,w10( 3)/.33874394455480d-01/
       data x10( 4)/.10366108297895d+01/,w10( 4)/.24013861108231d-00/
       data x10( 5)/.34290132722370d-00/,w10( 5)/.61086263373530d+00/

C * abscissas and weights for the gauss hermite integration order=11

       data x11( 1)/.36684708465596d+01/,w11( 1)/.14395603937142d-05/
       data x11( 2)/.27832900997816d+01/,w11( 2)/.34681946632334d-03/
       data x11( 3)/.20259480158257d+01/,w11( 3)/.11911395444911d-01/
       data x11( 4)/.13265570844949d+01/,w11( 4)/.11722787516770d-00/
       data x11( 5)/.65680956688210d+00/,w11( 5)/.42935975235611d-00/
       data x11( 6)/.00000000000000d 00/,w11( 6)/.65475928691459d+00/

C * abscissas and weights for the gauss hermite integration order=12

       data x12( 1)/.38897248978698d+01/,w12( 1)/.26585516843562d-06/
       data x12( 2)/.30206370251209d+01/,w12( 2)/.85736870435875d-04/
       data x12( 3)/.22795070805010d+01/,w12( 3)/.39053905846290d-02/
       data x12( 4)/.15976826351526d+01/,w12( 4)/.51607985615883d-01/
       data x12( 5)/.94778839124016d+00/,w12( 5)/.26049231026415d-00/
       data x12( 6)/.31424037625436d-00/,w12( 6)/.57013523626246d+00/

C * abscissas and weights for the gauss hermite integration order=13

       data x13( 1)/.41013375961786d+01/,w13( 1)/.48257318500728d-07/
       data x13( 2)/.32466089783724d+01/,w13( 2)/.20430360402706d-04/
       data x13( 3)/.25197356856782d+01/,w13( 3)/.12074599927193d-02/
       data x13( 4)/.18531076516015d+01/,w13( 4)/.20862775296169d-01/
       data x13( 5)/.12200550365907d+01/,w13( 5)/.14032332068702d-00/
       data x13( 6)/.60576387917106d+00/,w13( 6)/.42161629689853d-00/
       data x13( 7)/.00000000000000d 00/,w13( 7)/.60439318792116d+00/

C * abscissas and weights for the gauss hermite integration order=14

       data x14( 1)/.43044485704736d+01/,w14( 1)/.86285911681250d-08/
       data x14( 2)/.34626569336023d+01/,w14( 2)/.47164843550187d-05/
       data x14( 3)/.27484707249854d+01/,w14( 3)/.35509261355191d-03/
       data x14( 4)/.20951832585077d+01/,w14( 4)/.78500547264575d-02/
       data x14( 5)/.14766827311411d+01/,w14( 5)/.68505534223462d-01/
       data x14( 6)/.87871378732939d+00/,w14( 6)/.27310560906423d-00/
       data x14( 7)/.29174551067256d-00/,w14( 7)/.53640590971207d+00/

C * abscissas and weights for the gauss hermite integration order=15

       data x15( 1)/.44999907073094d+01/,w15( 1)/.15224758042535d-08/
       data x15( 2)/.36699503734044d+01/,w15( 2)/.10591155477110d-05/
       data x15( 3)/.29671669279056d+01/,w15( 3)/.10000444123250d-03/
       data x15( 4)/.23257324861738d+01/,w15( 4)/.27780688429126d-02/
       data x15( 5)/.17199925751865d+01/,w15( 5)/.30780033872545d-01/
       data x15( 6)/.11361155852109d+01/,w15( 6)/.15848891579593d-00/
       data x15( 7)/.56506958325557d+00/,w15( 7)/.41202868749887d-00/
       data x15( 8)/.00000000000000d 00/,w15( 8)/.56410030872642d+00/

C * abscissas and weights for the gauss hermite integration order=16

       data x16( 1)/.46887389393058d+01/,w16( 1)/.26548074740110d-09/
       data x16( 2)/.38694479048601d+01/,w16( 2)/.23209808448651d-06/
       data x16( 3)/.31769991619799d+01/,w16( 3)/.27118600925378d-04/
       data x16( 4)/.25462021578475d+01/,w16( 4)/.93228400862413d-03/
       data x16( 5)/.19517879909163d+01/,w16( 5)/.12880311535509d-01/
       data x16( 6)/.13802585391989d+01/,w16( 6)/.83810041398983d-01/
       data x16( 7)/.82295144914465d+00/,w16( 7)/.28064745852852d-00/
       data x16( 8)/.27348104613815d-00/,w16( 8)/.50792947901659d+00/

C * abscissas and weights for the gauss hermite integration order=17

       data x17( 1)/.48713451936744d+01/,w17( 1)/.45805789307984d-10/
       data x17( 2)/.40619466758755d+01/,w17( 2)/.49770789816306d-07/
       data x17( 3)/.33789320911415d+01/,w17( 3)/.71122891400209d-05/
       data x17( 4)/.27577629157039d+01/,w17( 4)/.29864328669774d-03/
       data x17( 5)/.21735028266666d+01/,w17( 5)/.50673499576273d-02/
       data x17( 6)/.16129243142212d+01/,w17( 6)/.40920034149754d-01/
       data x17( 7)/.10676487257434d+01/,w17( 7)/.17264829767009d-00/
       data x17( 8)/.53163300134265d+00/,w17( 8)/.40182646947039d-00/
       data x17( 9)/.00000000000000d 00/,w17( 9)/.53091793762486d+00/

C * abscissas and weights for the gauss hermite integration order=18

       data x18( 1)/.50483640088745d+01/,w18( 1)/.78281997721153d-11/
       data x18( 2)/.42481178735681d+01/,w18( 2)/.10467205795791d-07/
       data x18( 3)/.35737690684862d+01/,w18( 3)/.18106544810933d-05/
       data x18( 4)/.29613775055316d+01/,w18( 4)/.91811268679288d-04/
       data x18( 5)/.23862990891667d+01/,w18( 5)/.18885226302683d-02/
       data x18( 6)/.18355316042616d+01/,w18( 6)/.18640042387543d-01/
       data x18( 7)/.13009208583896d+01/,w18( 7)/.97301747641309d-01/
       data x18( 8)/.77668291926741d+00/,w18( 8)/.28480728566996d-00/
       data x18( 9)/.25826775051909d-00/,w18( 9)/.48349569472542d-00/

C * abscissas and weights for the gauss hermite integration order=19

       data x19( 1)/.52202716905375d+01/,w19( 1)/.13262970944984d-11/
       data x19( 2)/.44285328066038d+01/,w19( 2)/.21630510098635d-08/
       data x19( 3)/.37621873519640d+01/,w19( 3)/.44882431472230d-06/
       data x19( 4)/.31578488183476d+01/,w19( 4)/.27209197763160d-04/
       data x19( 5)/.25911337897945d+01/,w19( 5)/.67087752140714d-03/
       data x19( 6)/.20492317098506d+01/,w19( 6)/.79888667777225d-02/
       data x19( 7)/.15241706193935d+01/,w19( 7)/.50810386909049d-01/
       data x19( 8)/.10103683871343d+01/,w19( 8)/.18363270130699d-00/
       data x19( 9)/.50352016342389d+00/,w19( 9)/.39160898861301d-00/
       data x19(10)/.00000000000000d 00/,w19(10)/.50297488827618d+00/

C * abscissas and weights for the gauss hermite integration order=20

       data x20( 1)/.53874808900112d+01/,w20( 1)/.22293936455340d-12/
       data x20( 2)/.46036824495507d+01/,w20( 2)/.43993409922731d-09/
       data x20( 3)/.39447640401156d+01/,w20( 3)/.10860693707692d-06/
       data x20( 4)/.33478545673832d+01/,w20( 4)/.78025564785314d-05/
       data x20( 5)/.27888060584281d+01/,w20( 5)/.22833863601634d-03/
       data x20( 6)/.22549740020893d+01/,w20( 6)/.32437733422377d-02/
       data x20( 7)/.17385377121166d+01/,w20( 7)/.24810520887462d-01/
       data x20( 8)/.12340762153953d+01/,w20( 8)/.10901720602002d-00/
       data x20( 9)/.73747372854539d+00/,w20( 9)/.28667550536281d-00/
       data x20(10)/.24534070830090d-00/,w20(10)/.46224366960059d-00/

C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=21
C-----------------------------------------------------------------------
       DATA X21( 1)/.55503518732646D+01/,W21( 1)/.37203650701359D-13/
       DATA X21( 2)/.47739923434112D+01/,W21( 2)/.88186112420496D-10/
       DATA X21( 3)/.41219955474918D+01/,W21( 3)/.25712301800593D-07/
       DATA X21( 4)/.35319728771377D+01/,W21( 4)/.21718848980565D-05/
       DATA X21( 5)/.29799912077046D+01/,W21( 5)/.74783988673097D-04/
       DATA X21( 6)/.24535521245128D+01/,W21( 6)/.12549820417263D-02/
       DATA X21( 7)/.19449629491863D+01/,W21( 7)/.11414065837434D-01/
       DATA X21( 8)/.14489342506507D+01/,W21( 8)/.60179646658908D-01/
       DATA X21( 9)/.96149963441837D+00/,W21( 9)/.19212032406698D-00/
       DATA X21(10)/.47945070707910D-00/,W21(10)/.38166907361348D-00/
       DATA X21(11)/.00000000000000D 00/,W21(11)/.47902370312017D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=22
C-----------------------------------------------------------------------
       DATA X22( 1)/.57092013532053D+01/,W22( 1)/.61671834244035D-14/
       DATA X22( 2)/.49398341310601D+01/,W22( 2)/.17443390075478D-10/
       DATA X22( 3)/.42943124805932D+01/,W22( 3)/.59669909860593D-08/
       DATA X22( 4)/.37107015328778D+01/,W22( 4)/.58842875633007D-06/
       DATA X22( 5)/.31652659092021D+01/,W22( 5)/.23655128552509D-04/
       DATA X22( 6)/.26456374410582D+01/,W22( 6)/.46488505088422D-03/
       DATA X22( 7)/.21442335927985D+01/,W22( 7)/.49783993350513D-02/
       DATA X22( 8)/.16558743732864D+01/,W22( 8)/.31140370884423D-01/
       DATA X22( 9)/.11767139584812D+01/,W22( 9)/.11910236095878D-00/
       DATA X22(10)/.70368609717001D+00/,W22(10)/.28697143324689D-00/
       DATA X22(11)/.23417913993099D-00/,W22(11)/.44354522643494D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=23
C-----------------------------------------------------------------------
       DATA X23( 1)/.58643094989845D+01/,W23( 1)/.10160384620636D-14/
       DATA X23( 2)/.51015346104767D+01/,W23( 2)/.34083140980304D-11/
       DATA X23( 3)/.44620911737400D+01/,W23( 3)/.13596296504028D-08/
       DATA X23( 4)/.38844727081061D+01/,W23( 4)/.15553393291457D-06/
       DATA X23( 5)/.33451271599412D+01/,W23( 5)/.72492959180019D-05/
       DATA X23( 6)/.28318037871261D+01/,W23( 6)/.16556169914186D-03/
       DATA X23( 7)/.23370162114744D+01/,W23( 7)/.20695678749605D-02/
       DATA X23( 8)/.18556770376714D+01/,W23( 8)/.15207084004483D-01/
       DATA X23( 9)/.13840395856825D+01/,W23( 9)/.68890289429083D-01/
       DATA X23(10)/.91915146544256D+00/,W23(10)/.19864489857801D-00/
       DATA X23(11)/.45853835006810D-00/,W23(11)/.37214382487754D-00/
       DATA X23(12)/.00000000000000D 00/,W23(12)/.45819658559321D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=24
C-----------------------------------------------------------------------
       DATA X24( 1)/.60159255614257D+01/,W24( 1)/.16643684964890D-15/
       DATA X24( 2)/.52593829276680D+01/,W24( 2)/.65846202430777D-12/
       DATA X24( 3)/.46256627564238D+01/,W24( 3)/.30462542699874D-09/
       DATA X24( 4)/.40536644024481D+01/,W24( 4)/.40189711749411D-07/
       DATA X24( 5)/.35200068130345D+01/,W24( 5)/.21582457049022D-05/
       DATA X24( 6)/.30125461375656D+01/,W24( 6)/.56886916364041D-04/
       DATA X24( 7)/.25238810170114D+01/,W24( 7)/.82369248268837D-03/
       DATA X24( 8)/.20490035736617D+01/,W24( 8)/.70483558100722D-02/
       DATA X24( 9)/.15842500109617D+01/,W24( 9)/.37445470503229D-01/
       DATA X24(10)/.11267608176112D+01/,W24(10)/.12773962178455D-00/
       DATA X24(11)/.67417110703721D+00/,W24(11)/.28617953534642D-00/
       DATA X24(12)/.22441454747251D-00/,W24(12)/.42693116386868D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=25
C-----------------------------------------------------------------------
       DATA X25( 1)/.61642724340524D+01/,W25( 1)/.27119235140383D-16/
       DATA X25( 2)/.54136363552800D+01/,W25( 2)/.12588149877465D-12/
       DATA X25( 3)/.47853203673522D+01/,W25( 3)/.67196384177059D-10/
       DATA X25( 4)/.42186094443865D+01/,W25( 4)/.10170382503018D-07/
       DATA X25( 5)/.36902828769983D+01/,W25( 5)/.62570324996908D-06/
       DATA X25( 6)/.31882949244251D+01/,W25( 6)/.18915972957339D-04/
       DATA X25( 7)/.27053202371730D+01/,W25( 7)/.31508363874546D-03/
       DATA X25( 8)/.22364201302673D+01/,W25( 8)/.31157087201254D-02/
       DATA X25( 9)/.17780011243371D+01/,W25( 9)/.19243098965407D-01/
       DATA X25(10)/.13272807020731D+01/,W25(10)/.76888995175803D-01/
       DATA X25(11)/.88198275621381D+00/,W25(11)/.20362113667811D-00/
       DATA X25(12)/.44014729864531D-00/,W25(12)/.36308898927587D-00/
       DATA X25(13)/.00000000000000D 00/,W25(13)/.43986872216948D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=26
C-----------------------------------------------------------------------
       DATA X26( 1)/.63095503856257D+01/,W26( 1)/.43969160947533D-17/
       DATA X26( 2)/.55645249819501D+01/,W26( 2)/.23831486593720D-13/
       DATA X26( 3)/.49413249572414D+01/,W26( 3)/.14609999339815D-10/
       DATA X26( 4)/.43796026629833D+01/,W26( 4)/.25244940344904D-08/
       DATA X26( 5)/.38562884199091D+01/,W26( 5)/.17701063373973D-06/
       DATA X26( 6)/.33594271823508D+01/,W26( 6)/.61032917173956D-05/
       DATA X26( 7)/.28817622195431D+01/,W26( 7)/.11622970160310D-03/
       DATA X26( 8)/.24184157647738D+01/,W26( 8)/.13190647223238D-02/
       DATA X26( 9)/.19658547856411D+01/,W26( 9)/.93979012911589D-02/
       DATA X26(10)/.15213615166519D+01/,W26(10)/.43598227217247D-01/
       DATA X26(11)/.10827330110779D+01/,W26(11)/.13511332791178D-00/
       DATA X26(12)/.64809521399344D+00/,W26(12)/.28463224117677D-00/
       DATA X26(13)/.21577785624346D-00/,W26(13)/.41204365059034D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=27
C-----------------------------------------------------------------------
       DATA X27( 1)/.64519401407535D+01/,W27( 1)/.70957792970506D-18/
       DATA X27( 2)/.57122555528165D+01/,W27( 2)/.44707724573928D-14/
       DATA X27( 3)/.50939100031132D+01/,W27( 3)/.31341176136228D-11/
       DATA X27( 4)/.45369066633724D+01/,W27( 4)/.61550315782314D-09/
       DATA X27( 5)/.40183186704087D+01/,W27( 5)/.48954004096992D-07/
       DATA X27( 6)/.35262753401343D+01/,W27( 6)/.19152809005951D-05/
       DATA X27( 7)/.30535824198222D+01/,W27( 7)/.41467580043839D-04/
       DATA X27( 8)/.25954163389108D+01/,W27( 8)/.53676961568807D-03/
       DATA X27( 9)/.21482966453616D+01/,W27( 9)/.43812798357922D-02/
       DATA X27(10)/.17095607392603D+01/,W27(10)/.23415933625340D-01/
       DATA X27(11)/.12770668173399D+01/,W27(11)/.84173081084045D-01/
       DATA X27(12)/.84901134206010D+00/,W27(12)/.20737048075099D-00/
       DATA X27(13)/.42380790054385D-00/,W27(13)/.35451730409972D-00/
       DATA X27(14)/.00000000000000D 00/,W27(14)/.42357728801506D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=28
C-----------------------------------------------------------------------
       DATA X28( 1)/.65916054423677D+01/,W28( 1)/.11401393479036D-18/
       DATA X28( 2)/.58570146413828D+01/,W28( 2)/.83159379512061D-15/
       DATA X28( 3)/.52432853732029D+01/,W28( 3)/.66394367149093D-12/
       DATA X28( 4)/.46907565239431D+01/,W28( 4)/.14758531682776D-09/
       DATA X28( 5)/.41766367421293D+01/,W28( 5)/.13256825015416D-07/
       DATA X28( 6)/.36891342384617D+01/,W28( 6)/.58577197209925D-06/
       DATA X28( 7)/.32211120765614D+01/,W28( 7)/.14345504229713D-04/
       DATA X28( 8)/.27677953529136D+01/,W28( 8)/.21061810002402D-03/
       DATA X28( 9)/.23257498426564D+01/,W28( 9)/.19573312944088D-02/
       DATA X28(10)/.18923604968377D+01/,W28(10)/.11968423214354D-01/
       DATA X28(11)/.14655372634574D+01/,W28(11)/.49514889289895D-01/
       DATA X28(12)/.10435352737542D+01/,W28(12)/.14139460978694D-00/
       DATA X28(13)/.62483671950521D+00/,W28(13)/.28256139125936D-00/
       DATA X28(14)/.20806738269074D-00/,W28(14)/.39860471782642D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=29
C-----------------------------------------------------------------------
       DATA X29( 1)/.67286951986088D+01/,W29( 1)/.18244608527673D-19/
       DATA X29( 2)/.59989712894638D+01/,W29( 2)/.15345004446052D-15/
       DATA X29( 3)/.53896405219667D+01/,W29( 3)/.13901072714494D-12/
       DATA X29( 4)/.48413636510592D+01/,W29( 4)/.34841301613081D-10/
       DATA X29( 5)/.43314782938191D+01/,W29( 5)/.35203123276004D-08/
       DATA X29( 6)/.38482667922136D+01/,W29( 6)/.17492291299498D-06/
       DATA X29( 7)/.33846451410922D+01/,W29( 7)/.48230734976473D-05/
       DATA X29( 8)/.29358825042901D+01/,W29( 8)/.79909203545214D-04/
       DATA X29( 9)/.24985856910194D+01/,W29( 9)/.84079250614020D-03/
       DATA X29(10)/.20701810760534D+01/,W29(10)/.58455035452709D-02/
       DATA X29(11)/.16486229138923D+01/,W29(11)/.27639655592021D-01/
       DATA X29(12)/.12322157550847D+01/,W29(12)/.90768842215570D-01/
       DATA X29(13)/.81949868127091D+00/,W29(13)/.21014269444919D-00/
       DATA X29(14)/.40916463639493D-00/,W29(14)/.34641893907164D-00/
       DATA X29(15)/.00000000000000D 00/,W29(15)/.40897117463523D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=30
C-----------------------------------------------------------------------
       DATA X30( 1)/.68633452935299D+01/,W30( 1)/.29082547001310D-20/
       DATA X30( 2)/.61382792201239D+01/,W30( 2)/.28103336027506D-16/
       DATA X30( 3)/.55331471515674D+01/,W30( 3)/.28786070805485D-13/
       DATA X30( 4)/.49889189685899D+01/,W30( 4)/.81061862974624D-11/
       DATA X30( 5)/.44830553570925D+01/,W30( 5)/.91785804243775D-09/
       DATA X30( 6)/.40039086038612D+01/,W30( 6)/.51085224507756D-07/
       DATA X30( 7)/.35444438731553D+01/,W30( 7)/.15790948873246D-05/
       DATA X30( 8)/.30999705295864D+01/,W30( 8)/.29387252289227D-04/
       DATA X30( 9)/.26671321245356D+01/,W30( 9)/.34831012431866D-03/
       DATA X30(10)/.22433914677615D+01/,W30(10)/.27379224730675D-02/
       DATA X30(11)/.18267411436037D+01/,W30(11)/.14703829704825D-01/
       DATA X30(12)/.14155278001982D+01/,W30(12)/.55144176870230D-01/
       DATA X30(13)/.10083382710467D+01/,W30(13)/.14673584754088D-00/
       DATA X30(14)/.60392105862555D+00/,W30(14)/.28013093083919D-00/
       DATA X30(15)/.20112857654887D-00/,W30(15)/.38639488954178D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=31
C-----------------------------------------------------------------------
       DATA X31( 1)/.69956801237185D+01/,W31( 1)/.46189683944640D-21/
       DATA X31( 2)/.62750787049428D+01/,W31( 2)/.51106090079269D-17/
       DATA X31( 3)/.56739614446186D+01/,W31( 3)/.58995564987534D-14/
       DATA X31( 4)/.51335955771124D+01/,W31( 4)/.18603735214519D-11/
       DATA X31( 5)/.46315595063128D+01/,W31( 5)/.23524920032085D-09/
       DATA X31( 6)/.41562717558181D+01/,W31( 6)/.14611988344910D-07/
       DATA X31( 7)/.37007434032315D+01/,W31( 7)/.50437125589393D-06/
       DATA X31( 8)/.32603207323135D+01/,W31( 8)/.10498602757675D-04/
       DATA X31( 9)/.28316804533902D+01/,W31( 9)/.13952090395046D-03/
       DATA X31(10)/.24123177054804D+01/,W31(10)/.12336833073068D-02/
       DATA X31(11)/.20002585489356D+01/,W31(11)/.74827999140346D-02/
       DATA X31(12)/.15938858604721D+01/,W31(12)/.31847230731298D-01/
       DATA X31(13)/.11918269983500D+01/,W31(13)/.96717948160863D-01/
       DATA X31(14)/.79287697691531D+00/,W31(14)/.21213278866874D-00/
       DATA X31(15)/.39594273647142D-00/,W31(15)/.33877265789409D-00/
       DATA X31(16)/.00000000000000D 00/,W31(16)/.39577855609861D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=32
C-----------------------------------------------------------------------
       DATA X32( 1)/.71258139098307D+01/,W32( 1)/.73106764273835D-22/
       DATA X32( 2)/.64094981492696D+01/,W32( 2)/.92317365365176D-18/
       DATA X32( 3)/.58122259495159D+01/,W32( 3)/.11973440170928D-14/
       DATA X32( 4)/.52755509865159D+01/,W32( 4)/.42150102113259D-12/
       DATA X32( 5)/.47771645035026D+01/,W32( 5)/.59332914633963D-10/
       DATA X32( 6)/.43055479533512D+01/,W32( 6)/.40988321647706D-08/
       DATA X32( 7)/.38537554854714D+01/,W32( 7)/.15741677925455D-06/
       DATA X32( 8)/.34171674928186D+01/,W32( 8)/.36505851295621D-05/
       DATA X32( 9)/.29924908250024D+01/,W32( 9)/.54165840618195D-04/
       DATA X32(10)/.25772495377323D+01/,W32(10)/.53626836552793D-03/
       DATA X32(11)/.21694991836061D+01/,W32(11)/.36548903266541D-02/
       DATA X32(12)/.17676541094632D+01/,W32(12)/.17553428831572D-01/
       DATA X32(13)/.13703764109529D+01/,W32(13)/.60458130955907D-01/
       DATA X32(14)/.97650046358968D+00/,W32(14)/.15126973407663D-00/
       DATA X32(15)/.58497876543593D+00/,W32(15)/.27745814230251D-00/
       DATA X32(16)/.19484074156940D-00/,W32(16)/.37523835259277D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=33
C-----------------------------------------------------------------------
       DATA X33( 1)/.72538518220152D+01/,W33( 1)/.11533162185459D-22/
       DATA X33( 2)/.65416554457381D+01/,W33( 2)/.16570947415335D-18/
       DATA X33( 3)/.59480711820871D+01/,W33( 3)/.24077856795577D-15/
       DATA X33( 4)/.54149290026142D+01/,W33( 4)/.94348141590141D-13/
       DATA X33( 5)/.49200285205950D+01/,W33( 5)/.14739809370924D-10/
       DATA X33( 6)/.44519111488328D+01/,W33( 6)/.11289222471082D-08/
       DATA X33( 7)/.40036716099569D+01/,W33( 7)/.48077456763227D-07/
       DATA X33( 8)/.35707219802327D+01/,W33( 8)/.12376933672011D-05/
       DATA X33( 9)/.31497966817038D+01/,W33( 9)/.20423684051422D-04/
       DATA X33(10)/.27384458243513D+01/,W33(10)/.22544277059631D-03/
       DATA X33(11)/.23347511515295D+01/,W33(11)/.17184546377608D-02/
       DATA X33(12)/.19371545818222D+01/,W33(12)/.92656899706843D-02/
       DATA X33(13)/.15443482612431D+01/,W33(13)/.35987982318574D-01/
       DATA X33(14)/.11552002041268D+01/,W33(14)/.10206907999553D-00/
       DATA X33(15)/.76870137975886D+00/,W33(15)/.21349393115027D-00/
       DATA X33(16)/.38392601450841D-00/,W33(16)/.33155200075071D-00/
       DATA X33(17)/.00000000000000D 00/,W33(17)/.38378526651986D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=34
C-----------------------------------------------------------------------
       DATA X34( 1)/.73798909504812D+01/,W34( 1)/.18138001119595D-23/
       DATA X34( 2)/.66716591360701D+01/,W34( 2)/.29567089223603D-19/
       DATA X34( 3)/.60816169939363D+01/,W34( 3)/.47999019978943D-16/
       DATA X34( 4)/.55518613309888D+01/,W34( 4)/.20878403731158D-13/
       DATA X34( 5)/.50602960186057D+01/,W34( 5)/.36098817474857D-11/
       DATA X34( 6)/.45955197481081D+01/,W34( 6)/.30562520419152D-09/
       DATA X34( 7)/.41506656029708D+01/,W34( 7)/.14387732912538D-07/
       DATA X34( 8)/.37211752324761D+01/,W34( 8)/.40979740352242D-06/
       DATA X34( 9)/.33038084315644D+01/,W34( 9)/.74934487833017D-05/
       DATA X34(10)/.28961389431744D+01/,W34(10)/.91861189828716D-04/
       DATA X34(11)/.24962719408165D+01/,W34(11)/.77981759962315D-03/
       DATA X34(12)/.21026736904673D+01/,W34(12)/.46984636292662D-02/
       DATA X34(13)/.17140625533873D+01/,W34(13)/.20473151727018D-01/
       DATA X34(14)/.13293355518848D+01/,W34(14)/.65445134108747D-01/
       DATA X34(15)/.94751645803344D+00/,W34(15)/.15511041662329D-00/
       DATA X34(16)/.56771726855487D+00/,W34(16)/.27462771563511D-00/
       DATA X34(17)/.18910806052714D-00/,W34(17)/.36499244699662D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=35
C-----------------------------------------------------------------------
       DATA X35( 1)/.75040211464489D+01/,W35( 1)/.28441134657250D-24/
       DATA X35( 2)/.67996094132841D+01/,W35( 2)/.52456527291741D-20/
       DATA X35( 3)/.62129737476337D+01/,W35( 3)/.94898848794720D-17/
       DATA X35( 4)/.56864689480904D+01/,W35( 4)/.45704441904179D-14/
       DATA X35( 5)/.51980993461977D+01/,W35( 5)/.87225260125362D-12/
       DATA X35( 6)/.47365184774132D+01/,W35( 6)/.81406956866768D-10/
       DATA X35( 7)/.42948958144927D+01/,W35( 7)/.42239761300090D-08/
       DATA X35( 8)/.38687007309691D+01/,W35( 8)/.13269405173065D-06/
       DATA X35( 9)/.34547164957520D+01/,W35( 9)/.26798156546133D-05/
       DATA X35(10)/.30505384204304D+01/,W35(10)/.36352769174154D-04/
       DATA X35(11)/.26542927811972D+01/,W35(11)/.34234009319958D-03/
       DATA X35(12)/.22644675010426D+01/,W35(12)/.22950283295839D-02/
       DATA X35(13)/.18798039887309D+01/,W35(13)/.11166806590265D-01/
       DATA X35(14)/.14992244886117D+01/,W35(14)/.40024775133098D-01/
       DATA X35(15)/.11217809907203D+01/,W35(15)/.10687290695545D-00/
       DATA X35(16)/.74661763987986D+00/,W35(16)/.21434719059605D-00/
       DATA X35(17)/.37294171704962D-00/,W35(17)/.32472872157454D-00/
       DATA X35(18)/.00000000000000D 00/,W35(18)/.37281997319072D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=36
C-----------------------------------------------------------------------
       DATA X36( 1)/.76263257540039D+01/,W36( 1)/.44471534175756D-25/
       DATA X36( 2)/.69255989902599D+01/,W36( 2)/.92564030835575D-21/
       DATA X36( 3)/.63422433309944D+01/,W36( 3)/.18615979878079D-17/
       DATA X36( 4)/.58188632795055D+01/,W36( 4)/.99029016869057D-15/
       DATA X36( 5)/.53335601071130D+01/,W36( 5)/.20809126742867D-12/
       DATA X36( 6)/.48750399724671D+01/,W36( 6)/.21353541817843D-10/
       DATA X36( 7)/.44365069701928D+01/,W36( 7)/.12178678209556D-08/
       DATA X36( 8)/.40134565677495D+01/,W36( 8)/.42074913231087D-07/
       DATA X36( 9)/.36026938571485D+01/,W36( 9)/.93555700897725D-06/
       DATA X36(10)/.32018339457882D+01/,W36(10)/.13997029220328D-04/
       DATA X36(11)/.28090222351311D+01/,W36(11)/.14569615962635D-03/
       DATA X36(12)/.24227660420535D+01/,W36(12)/.10825347566867D-02/
       DATA X36(13)/.20418271835544D+01/,W36(13)/.58564259865590D-02/
       DATA X36(14)/.16651500018434D+01/,W36(14)/.23425767509729D-01/
       DATA X36(15)/.12918109588209D+01/,W36(15)/.70104750105171D-01/
       DATA X36(16)/.92098180157075D+00/,W36(16)/.15835545375115D-00/
       DATA X36(17)/.55190143329042D+00/,W36(17)/.27170124700950D-00/
       DATA X36(18)/.18385336710581D-00/,W36(18)/.35554007427368D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=37
C-----------------------------------------------------------------------
       DATA X37( 1)/.77468822496494D+01/,W37( 1)/.69350835505832D-26/
       DATA X37( 2)/.70497138557782D+01/,W37( 2)/.16249853271917D-21/
       DATA X37( 3)/.64695200365240D+01/,W37( 3)/.36247269512440D-18/
       DATA X37( 4)/.59491472174620D+01/,W37( 4)/.21248851941313D-15/
       DATA X37( 5)/.54667903359685D+01/,W37( 5)/.49047393641849D-13/
       DATA X37( 6)/.50112061385731D+01/,W37( 6)/.55203972646368D-11/
       DATA X37( 7)/.45756317486673D+01/,W37( 7)/.34519062272117D-09/
       DATA X37( 8)/.41555872811264D+01/,W37( 8)/.13079698821801D-07/
       DATA X37( 9)/.37478982064755D+01/,W37( 9)/.31928854810029D-06/
       DATA X37(10)/.33501978949725D+01/,W37(10)/.52521435771195D-05/
       DATA X37(11)/.29606491813033D+01/,W37(11)/.60226951652668D-04/
       DATA X37(12)/.25777768581133D+01/,W37(12)/.49417150121465D-03/
       DATA X37(13)/.22003609340093D+01/,W37(13)/.29608001846026D-02/
       DATA X37(14)/.18273652487636D+01/,W37(14)/.13160179977210D-01/
       DATA X37(15)/.14578876468742D+01/,W37(15)/.43931362299535D-01/
       DATA X37(16)/.10911237649759D+01/,W37(16)/.11117916510719D-00/
       DATA X37(17)/.72633961660512D+00/,W37(17)/.21478887596477D-00/
       DATA X37(18)/.36284990505066D-00/,W37(18)/.31827467975440D-00/
       DATA X37(19)/.00000000000000D 00/,W37(19)/.36274375769908D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=38
C-----------------------------------------------------------------------
       DATA X38( 1)/.78657628033800D+01/,W38( 1)/.10787188820740D-26/
       DATA X38( 2)/.71720339353200D+01/,W38( 2)/.28387387532092D-22/
       DATA X38( 3)/.65948913272655D+01/,W38( 3)/.70078433305578D-19/
       DATA X38( 4)/.60774160035375D+01/,W38( 4)/.45173712236042D-16/
       DATA X38( 5)/.55978935141847D+01/,W38( 5)/.11428616750256D-13/
       DATA X38( 6)/.51451293207408D+01/,W38( 6)/.14076333438451D-11/
       DATA X38( 7)/.47123921320849D+01/,W38( 7)/.96270292926046D-10/
       DATA X38( 8)/.42952254197496D+01/,W38( 8)/.39906528856885D-08/
       DATA X38( 9)/.38904737609633D+01/,W38( 9)/.10665865704408D-06/
       DATA X38(10)/.34957874548356D+01/,W38(10)/.19234656867322D-05/
       DATA X38(11)/.31093453117179D+01/,W38(11)/.24223415120933D-04/
       DATA X38(12)/.27296879628883D+01/,W38(12)/.21875759533937D-03/
       DATA X38(13)/.23556117330355D+01/,W38(13)/.14463071215625D-02/
       DATA X38(14)/.19860977780391D+01/,W38(14)/.71152368261332D-02/
       DATA X38(15)/.16202627556330D+01/,W38(15)/.26380533520484D-01/
       DATA X38(16)/.12573231317007D+01/,W38(16)/.74443983368231D-01/
       DATA X38(17)/.89656834619313D+00/,W38(17)/.16108797149546D-00/
       DATA X38(18)/.53733981087098D+00/,W38(18)/.26872376002890D-00/
       DATA X38(19)/.17901372329588D-00/,W38(19)/.34678411786873D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=39
C-----------------------------------------------------------------------
       DATA X39( 1)/.79830347727198D+01/,W39( 1)/.16737892457495D-27/
       DATA X39( 2)/.72926336708657D+01/,W39( 2)/.49359073375899D-23/
       DATA X39( 3)/.67184385064441D+01/,W39( 3)/.13457268201343D-19/
       DATA X39( 4)/.62037579977281D+01/,W39( 4)/.95193204821151D-17/
       DATA X39( 5)/.57269654517821D+01/,W39( 5)/.26340961107221D-14/
       DATA X39( 6)/.52769133152304D+01/,W39( 6)/.35426508019788D-12/
       DATA X39( 7)/.48469005687435D+01/,W39( 7)/.26440048769928D-10/
       DATA X39( 8)/.44324928825930D+01/,W39( 8)/.11961664237666D-08/
       DATA X39( 9)/.40305528146024D+01/,W39( 9)/.34915088439104D-07/
       DATA X39(10)/.36387464248745D+01/,W39(10)/.68844204093396D-06/
       DATA X39(11)/.32552672359922D+01/,W39(11)/.94942995361406D-05/
       DATA X39(12)/.28786703113749D+01/,W39(12)/.94077279467663D-04/
       DATA X39(13)/.25077666938913D+01/,W39(13)/.68405759521850D-03/
       DATA X39(14)/.21415530119869D+01/,W39(14)/.37112302317138D-02/
       DATA X39(15)/.17791625828543D+01/,W39(15)/.15221980837288D-01/
       DATA X39(16)/.14198301576857D+01/,W39(16)/.47690071331630D-01/
       DATA X39(17)/.10628655672812D+01/,W39(17)/.11503504579976D-00/
       DATA X39(18)/.70763327334857D+00/,W39(18)/.21489607677744D-00/
       DATA X39(19)/.35353584699633D-00/,W39(19)/.31216284886714D-00/
       DATA X39(20)/.00000000000000D 00/,W39(20)/.35344263570680D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=40
C-----------------------------------------------------------------------
       DATA X40( 1)/.80987611392508D+01/,W40( 1)/.25910437138468D-28/
       DATA X40( 2)/.74115825314854D+01/,W40( 2)/.85440569637745D-24/
       DATA X40( 3)/.68402373052493D+01/,W40( 3)/.25675933654114D-20/
       DATA X40( 4)/.63282553512201D+01/,W40( 4)/.19891810121162D-17/
       DATA X40( 5)/.58540950560304D+01/,W40( 5)/.60083587894902D-15/
       DATA X40( 6)/.54066542479701D+01/,W40( 6)/.88057076452151D-13/
       DATA X40( 7)/.49792609785452D+01/,W40( 7)/.71565280526894D-11/
       DATA X40( 8)/.45675020728444D+01/,W40( 8)/.35256207913652D-09/
       DATA X40( 9)/.41682570668325D+01/,W40( 9)/.11212360832275D-07/
       DATA X40(10)/.37792067534352D+01/,W40(10)/.24111441636703D-06/
       DATA X40(11)/.33985582658596D+01/,W40(11)/.36315761506926D-05/
       DATA X40(12)/.30248798839013D+01/,W40(12)/.39369339810921D-04/
       DATA X40(13)/.26569959984429D+01/,W40(13)/.31385359454130D-03/
       DATA X40(14)/.22939171418751D+01/,W40(14)/.18714968295977D-02/
       DATA X40(15)/.19347914722823D+01/,W40(15)/.84608880082571D-02/
       DATA X40(16)/.15788698949316D+01/,W40(16)/.29312565536169D-01/
       DATA X40(17)/.12254801090463D+01/,W40(17)/.78474605865395D-01/
       DATA X40(18)/.87400661235708D+00/,W40(18)/.16337873271326D-00/
       DATA X40(19)/.52387471383228D+00/,W40(19)/.26572825187735D-00/
       DATA X40(20)/.17453721459758D-00/,W40(20)/.33864327742555D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=41
C-----------------------------------------------------------------------
       DATA X41( 1)/.82130008955983D+01/,W41( 1)/.40019596646658D-29/
       DATA X41( 2)/.75289454645396D+01/,W41( 2)/.14726537286519D-24/
       DATA X41( 3)/.69603584006367D+01/,W41( 3)/.48687379360127D-21/
       DATA X41( 4)/.64509845971747D+01/,W41( 4)/.41234085375311D-18/
       DATA X41( 5)/.59793650041651D+01/,W41( 5)/.13569875295860D-15/
       DATA X41( 6)/.55344413406134D+01/,W41( 6)/.21629872471748D-13/
       DATA X41( 7)/.51095696265331D+01/,W41( 7)/.19103833646807D-11/
       DATA X41( 8)/.47003568963041D+01/,W41( 8)/.10226798927781D-09/
       DATA X41( 9)/.43036987671546D+01/,W41( 9)/.35356816483804D-08/
       DATA X41(10)/.39172898548378D+01/,W41(10)/.82726584187489D-07/
       DATA X41(11)/.35393499373637D+01/,W41(11)/.13573781404869D-05/
       DATA X41(12)/.31684594539420D+01/,W41(12)/.16055965736676D-04/
       DATA X41(13)/.28034549614843D+01/,W41(13)/.13993025659672D-03/
       DATA X41(14)/.24433595531234D+01/,W41(14)/.91423425632859D-03/
       DATA X41(15)/.20873346819187D+01/,W41(15)/.45403929873028D-02/
       DATA X41(16)/.17346456088220D+01/,W41(16)/.17330881362115D-01/
       DATA X41(17)/.13846357891600D+01/,W41(17)/.51289872397089D-01/
       DATA X41(18)/.10367072529242D+01/,W41(18)/.11848439124683D-00/
       DATA X41(19)/.69030505233021D+00/,W41(19)/.21473086449006D-00/
       DATA X41(20)/.34490446301543D-00/,W41(20)/.30636781693782D-00/
       DATA X41(21)/.00000000000000D 00/,W41(21)/.34482208361639D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=42
C-----------------------------------------------------------------------
       DATA X42( 1)/.83258093895669D+01/,W42( 1)/.61678589258105D-30/
       DATA X42( 2)/.76447832957047D+01/,W42( 2)/.25278698640532D-25/
       DATA X42( 3)/.70788678730491D+01/,W42( 3)/.91778906956913D-22/
       DATA X42( 4)/.65720171713874D+01/,W42( 4)/.84821520800855D-19/
       DATA X42( 5)/.61028523343815D+01/,W42( 5)/.30358903478104D-16/
       DATA X42( 6)/.56603575812830D+01/,W42( 6)/.52533377155679D-14/
       DATA X42( 7)/.52379158850176D+01/,W42( 7)/.50327055821835D-12/
       DATA X42( 8)/.48311536291282D+01/,W42( 8)/.29217288372330D-10/
       DATA X42( 9)/.44369817058810D+01/,W42( 9)/.10958052288077D-08/
       DATA X42(10)/.40531077444247D+01/,W42(10)/.27834715265488D-07/
       DATA X42(11)/.36777633163885D+01/,W42(11)/.49636593935792D-06/
       DATA X42(12)/.33095400965109D+01/,W42(12)/.63902459677345D-05/
       DATA X42(13)/.29472857823055D+01/,W42(13)/.60719621077877D-04/
       DATA X42(14)/.25900348706171D+01/,W42(14)/.43341227172120D-03/
       DATA X42(15)/.22369607870543D+01/,W42(15)/.23571613945960D-02/
       DATA X42(16)/.18873416205435D+01/,W42(16)/.98795240531873D-02/
       DATA X42(17)/.15405348009155D+01/,W42(17)/.32202101288904D-01/
       DATA X42(18)/.11959577943778D+01/,W42(18)/.82211269303286D-01/
       DATA X42(19)/.85307290916055D+00/,W42(19)/.16528800127465D-00/
       DATA X42(20)/.51137491831547D+00/,W42(20)/.26273890678226D-00/
       DATA X42(21)/.17038058556182D-00/,W42(21)/.33104891389082D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=43
C-----------------------------------------------------------------------
       DATA X43( 1)/.84372386310833D+01/,W43( 1)/.94863063412557D-31/
       DATA X43( 2)/.77591530847325D+01/,W43( 2)/.43221460599334D-26/
       DATA X43( 3)/.71958276123464D+01/,W43( 3)/.17203368260280D-22/
       DATA X43( 4)/.66914198727121D+01/,W43( 4)/.17320753790912D-19/
       DATA X43( 5)/.62246289668942D+01/,W43( 5)/.67308117032423D-17/
       DATA X43( 6)/.57844803140775D+01/,W43( 6)/.12622032209756D-14/
       DATA X43( 7)/.53643829011515D+01/,W43( 7)/.13092104166583D-12/
       DATA X43( 8)/.49599816751949D+01/,W43( 8)/.82271061208154D-11/
       DATA X43( 9)/.45682020755441D+01/,W43( 9)/.33407179153172D-09/
       DATA X43(10)/.41867640213660D+01/,W43(10)/.91931960675394D-08/
       DATA X43(11)/.38139101240654D+01/,W43(11)/.17777830445478D-06/
       DATA X43(12)/.34482424822200D+01/,W43(12)/.24851377424112D-05/
       DATA X43(13)/.30886190396004D+01/,W43(13)/.25680573077054D-04/
       DATA X43(14)/.27340846945379D+01/,W43(14)/.19972435671321D-03/
       DATA X43(15)/.23838237072843D+01/,W43(15)/.11860700927673D-02/
       DATA X43(16)/.20371256888642D+01/,W43(16)/.54414552932410D-02/
       DATA X43(17)/.16933605309399D+01/,W43(17)/.19468160151099D-01/
       DATA X43(18)/.13519593686709D+01/,W43(18)/.54724809081869D-01/
       DATA X43(19)/.10123996845633D+01/,W43(19)/.12156737988637D-00/
       DATA X43(20)/.67419327674231D+00/,W43(20)/.21434349700186D-00/
       DATA X43(21)/.33687619662553D-00/,W43(21)/.30086599386762D-00/
       DATA X43(22)/.00000000000000D 00/,W43(22)/.33680296539275D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=44
C-----------------------------------------------------------------------
       DATA X44( 1)/.85473375667355D+01/,W44( 1)/.14561153081760D-31/
       DATA X44( 2)/.78721084427748D+01/,W44( 2)/.73621261042957D-27/
       DATA X44( 3)/.73112956789163D+01/,W44( 3)/.32072216919047D-23/
       DATA X44( 4)/.68092552715358D+01/,W44( 4)/.35121449558363D-20/
       DATA X44( 5)/.63447621648968D+01/,W44( 5)/.14794099758170D-17/
       DATA X44( 6)/.59068817590275D+01/,W44( 6)/.30015214029805D-15/
       DATA X44( 7)/.54890481838807D+01/,W44( 7)/.33650356127919D-13/
       DATA X44( 8)/.50869242297582D+01/,W44( 8)/.22848078971217D-11/
       DATA X44( 9)/.46974492230014D+01/,W44( 9)/.10026000778908D-09/
       DATA X44(10)/.43183547234421D+01/,W44(10)/.29831210190393D-08/
       DATA X44(11)/.39478937115247D+01/,W44(11)/.62427245062806D-07/
       DATA X44(12)/.35846779931052D+01/,W44(12)/.94545385422844D-06/
       DATA X44(13)/.32275749742774D+01/,W44(13)/.10600171772516D-04/
       DATA X44(14)/.28756390821704D+01/,W44(14)/.89597973640606D-04/
       DATA X44(15)/.25280644271477D+01/,W44(15)/.57942043271559D-03/
       DATA X44(16)/.21841510195665D+01/,W44(16)/.29012937752962D-02/
       DATA X44(17)/.18432799206176D+01/,W44(17)/.11357866939068D-01/
       DATA X44(18)/.15048944344829D+01/,W44(18)/.35033757109683D-01/
       DATA X44(19)/.11684854631902D+01/,W44(19)/.85670187300236D-01/
       DATA X44(20)/.83357976100207D+00/,W44(20)/.16686716582562D-00/
       DATA X44(21)/.49973020830067D-00/,W44(21)/.25977339611585D-00/
       DATA X44(22)/.16650746707363D-00/,W44(22)/.32394262884195D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=45
C-----------------------------------------------------------------------
       DATA X45( 1)/.86561523259903D+01/,W45( 1)/.22308126620262D-32/
       DATA X45( 2)/.79836998162220D+01/,W45( 2)/.12494829120557D-27/
       DATA X45( 3)/.74253266258562D+01/,W45( 3)/.59481747159219D-24/
       DATA X45( 4)/.69255820732776D+01/,W45( 4)/.70737246502691D-21/
       DATA X45( 5)/.64633149436642D+01/,W45( 5)/.32248298707323D-18/
       DATA X45( 6)/.60276294722537D+01/,W45( 6)/.70674619157535D-16/
       DATA X45( 7)/.56119841218397D+01/,W45( 7)/.85500874663691D-14/
       DATA X45( 8)/.52120588631620D+01/,W45( 8)/.62620383158347D-12/
       DATA X45( 9)/.48248063087032D+01/,W45( 9)/.29642049893154D-10/
       DATA X45(10)/.44479690736019D+01/,W45(10)/.95182669412326D-09/
       DATA X45(11)/.40798099079254D+01/,W45(11)/.21512839074192D-07/
       DATA X45(12)/.37189496893415D+01/,W45(12)/.35225192415061D-06/
       DATA X45(13)/.33642645958565D+01/,W45(13)/.42754558863968D-05/
       DATA X45(14)/.30148178197492D+01/,W45(14)/.39183295448815D-04/
       DATA X45(15)/.26698124652513D+01/,W45(15)/.27524242709949D-03/
       DATA X45(16)/.23285579714532D+01/,W45(16)/.15001278829296D-02/
       DATA X45(17)/.19904454586466D+01/,W45(17)/.64070706562032D-02/
       DATA X45(18)/.16549291191548D+01/,W45(18)/.21617644433763D-01/
       DATA X45(19)/.13215117856392D+01/,W45(19)/.57992746631622D-01/
       DATA X45(20)/.98973342486571D+00/,W45(20)/.12432050192892D-00/
       DATA X45(21)/.65916168887414D+00/,W45(21)/.21377487850199D-00/
       DATA X45(22)/.32938389969667D-00/,W45(22)/.29563565196659D-00/
       DATA X45(23)/.00000000000000D 00/,W45(23)/.32931845505069D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=46
C-----------------------------------------------------------------------
       DATA X46( 1)/.87637264425763D+01/,W46( 1)/.34113792299091D-33/
       DATA X46( 2)/.80939747412671D+01/,W46( 2)/.21132014962693D-28/
       DATA X46( 3)/.75379717872258D+01/,W46( 3)/.10976558019966D-24/
       DATA X46( 4)/.70404554428052D+01/,W46( 4)/.14155042016753D-21/
       DATA X46( 5)/.65803464347374D+01/,W46( 5)/.69737995966938D-19/
       DATA X46( 6)/.61467867545054D+01/,W46( 6)/.16484447389505D-16/
       DATA X46( 7)/.57332584418499D+01/,W46( 7)/.21486472823676D-14/
       DATA X46( 8)/.53354580360572D+01/,W46( 8)/.16947076167748D-12/
       DATA X46( 9)/.49503508867398D+01/,W46( 9)/.86391402914471D-11/
       DATA X46(10)/.45756901345704D+01/,W46(10)/.29885455879920D-09/
       DATA X46(11)/.42097477613270D+01/,W46(11)/.72816347565630D-08/
       DATA X46(12)/.38511531505025D+01/,W46(12)/.12865324130562D-06/
       DATA X46(13)/.34987906412463D+01/,W46(13)/.16869468397353D-05/
       DATA X46(14)/.31517315196385D+01/,W46(14)/.16725935920652D-04/
       DATA X46(15)/.28091871354773D+01/,W46(15)/.12731983755105D-03/
       DATA X46(16)/.24704755132281D+01/,W46(16)/.75339673454306D-03/
       DATA X46(17)/.21349968729707D+01/,W46(17)/.35010578505224D-02/
       DATA X46(18)/.18022151950575D+01/,W46(18)/.12883487559776D-01/
       DATA X46(19)/.14716439717321D+01/,W46(19)/.37795830992242D-01/
       DATA X46(20)/.11428349031264D+01/,W46(20)/.88868227565993D-01/
       DATA X46(21)/.81536867934468D+00/,W46(21)/.16816012736210D-00/
       DATA X46(22)/.48884723585117D-00/,W46(22)/.25684454350649D-00/
       DATA X46(23)/.16288702795526D-00/,W46(23)/.31727438491813D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=47
C-----------------------------------------------------------------------
       DATA X47( 1)/.88701010540232D+01/,W47( 1)/.52074550906397D-34/
       DATA X47( 2)/.82029780727978D+01/,W47( 2)/.35619817807810D-29/
       DATA X47( 3)/.76492795366290D+01/,W47( 3)/.20158558166545D-25/
       DATA X47( 4)/.71539272947757D+01/,W47( 4)/.28149547423917D-22/
       DATA X47( 5)/.66959122110632D+01/,W47( 5)/.14966299166255D-19/
       DATA X47( 6)/.62644130149705D+01/,W47( 6)/.38101463479225D-17/
       DATA X47( 7)/.58529346161254D+01/,W47( 7)/.53428423961479D-15/
       DATA X47( 8)/.54571895556816D+01/,W47( 8)/.45312676467473D-13/
       DATA X47( 9)/.50741554169850D+01/,W47( 9)/.24836135431245D-11/
       DATA X47(10)/.47015953849679D+01/,W47(10)/.92402978883634D-10/
       DATA X47(11)/.43377901884924D+01/,W47(11)/.24228134937241D-08/
       DATA X47(12)/.39813772107065D+01/,W47(12)/.46104364856629D-07/
       DATA X47(13)/.36312484092926D+01/,W47(13)/.65181055230394D-06/
       DATA X47(14)/.32864825713247D+01/,W47(14)/.69771079617897D-05/
       DATA X47(15)/.29462986351922D+01/,W47(15)/.57425938238314D-04/
       DATA X47(16)/.26100224744061D+01/,W47(16)/.36806169587481D-03/
       DATA X47(17)/.22770625758595D+01/,W47(17)/.18562527788730D-02/
       DATA X47(18)/.19468917579992D+01/,W47(18)/.74296915378891D-02/
       DATA X47(19)/.16190330813590D+01/,W47(19)/.23765555762403D-01/
       DATA X47(20)/.12930487170109D+01/,W47(20)/.61094386221793D-01/
       DATA X47(21)/.96853092273781D+00/,W47(21)/.12677670015499D-00/
       DATA X47(22)/.64509452302768D+00/,W47(22)/.21305846026512D-00/
       DATA X47(23)/.32237044915615D-00/,W47(23)/.29065687385121D-00/
       DATA X47(24)/.00000000000000D 00/,W47(24)/.32231167941131D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=48
C-----------------------------------------------------------------------
       DATA X48( 1)/.89753150819316D+01/,W48( 1)/.79355514607730D-35/
       DATA X48( 2)/.83107521907047D+01/,W48( 2)/.59846126933129D-30/
       DATA X48( 3)/.77592955197658D+01/,W48( 3)/.36850360801502D-26/
       DATA X48( 4)/.72660465541643D+01/,W48( 4)/.55645774689014D-23/
       DATA X48( 5)/.68100645780741D+01/,W48( 5)/.31883873235044D-20/
       DATA X48( 6)/.63805640961864D+01/,W48( 6)/.87301596011857D-18/
       DATA X48( 7)/.59710722250135D+01/,W48( 7)/.13151596226583D-15/
       DATA X48( 8)/.55773169812237D+01/,W48( 8)/.11975898654790D-13/
       DATA X48( 9)/.51962877187923D+01/,W48( 9)/.70469325815448D-12/
       DATA X48(10)/.48257572281332D+01/,W48(10)/.28152965378377D-10/
       DATA X48(11)/.44640145469345D+01/,W48(11)/.79304674951642D-09/
       DATA X48(12)/.41097046035606D+01/,W48(12)/.16225141358956D-07/
       DATA X48(13)/.37617264902283D+01/,W48(13)/.24686589936694D-06/
       DATA X48(14)/.34191659693639D+01/,W48(14)/.28472586917345D-05/
       DATA X48(15)/.30812489886451D+01/,W48(15)/.25285990277481D-04/
       DATA X48(16)/.27473086248224D+01/,W48(16)/.17515043180115D-03/
       DATA X48(17)/.24167609048732D+01/,W48(17)/.95639231981928D-03/
       DATA X48(18)/.20890866609443D+01/,W48(18)/.41530049119770D-02/
       DATA X48(19)/.17638175798953D+01/,W48(19)/.14444961574979D-01/
       DATA X48(20)/.14405252201376D+01/,W48(20)/.40479676984599D-01/
       DATA X48(21)/.11188121524022D+01/,W48(21)/.91822297079269D-01/
       DATA X48(22)/.79830462777856D+00/,W48(22)/.16920447194561D-00/
       DATA X48(23)/.47864633759449D-00/,W48(23)/.25396154266473D-00/
       DATA X48(24)/.15949293584886D-00/,W48(24)/.31100103037792D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=49
C-----------------------------------------------------------------------
       DATA X49( 1)/.90794053951994D+01/,W49( 1)/.12072872881282D-35/
       DATA X49( 2)/.84173371862679D+01/,W49( 2)/.10023613427435D-30/
       DATA X49( 3)/.78680628640816D+01/,W49( 3)/.67063563353452D-27/
       DATA X49( 4)/.73768593906319D+01/,W49( 4)/.10936809436126D-23/
       DATA X49( 5)/.69228528349597D+01/,W49( 5)/.67446743456580D-21/
       DATA X49( 6)/.64952925650076D+01/,W49( 6)/.19836416602791D-18/
       DATA X49( 7)/.60877272810547D+01/,W49( 7)/.32059589524000D-16/
       DATA X49( 8)/.56958999850352D+01/,W49( 8)/.31301391441239D-14/
       DATA X49( 9)/.53168113741835D+01/,W49( 9)/.19744849360629D-12/
       DATA X49(10)/.49482434430049D+01/,W49(10)/.84575241714300D-11/
       DATA X49(11)/.45884931402745D+01/,W49(11)/.25554702439985D-09/
       DATA X49(12)/.42362125300688D+01/,W49(12)/.56118428733868D-08/
       DATA X49(13)/.38903074057662D+01/,W49(13)/.91729613211690D-07/
       DATA X49(14)/.35498700372121D+01/,W49(14)/.11378411996772D-05/
       DATA X49(15)/.32141328684568D+01/,W49(15)/.10881805485673D-04/
       DATA X49(16)/.28824356106456D+01/,W49(16)/.81291046453144D-04/
       DATA X49(17)/.25542011939519D+01/,W49(17)/.47952088600471D-03/
       DATA X49(18)/.22289177453859D+01/,W49(18)/.22536737753385D-02/
       DATA X49(19)/.19061247609884D+01/,W49(19)/.85018075629536D-02/
       DATA X49(20)/.15854022514550D+01/,W49(20)/.25900304438339D-01/
       DATA X49(21)/.12663620236874D+01/,W49(21)/.64032496212200D-01/
       DATA X49(22)/.94864050761025D+00/,W49(22)/.12896559830847D-00/
       DATA X49(23)/.63189269951096D+00/,W49(23)/.21222171925814D-00/
       DATA X49(24)/.31578690052375D-00/,W49(24)/.28591145169463D-00/
       DATA X49(25)/.00000000000000D 00/,W49(25)/.31573389003557D-00/
C-----------------------------------------------------------------------
C-----ABSCISSAS AND WEIGHTS FOR THE GAUSS HERMITE INTEGRATION ORDER=50
C-----------------------------------------------------------------------
       DATA X50( 1)/.91824069581293D+01/,W50( 1)/.18337940485732D-36/
       DATA X50( 2)/.85227710309178D+01/,W50( 2)/.16738016679076D-31/
       DATA X50( 3)/.79756223682056D+01/,W50( 3)/.12152441234043D-27/
       DATA X50( 4)/.74864094298642D+01/,W50( 4)/.21376583083598D-24/
       DATA X50( 5)/.70343235097706D+01/,W50( 5)/.14170935995732D-21/
       DATA X50( 6)/.66086479738553D+01/,W50( 6)/.44709843654072D-19/
       DATA X50( 7)/.62029525192747D+01/,W50( 7)/.77423829570423D-17/
       DATA X50( 8)/.58129946754204D+01/,W50( 8)/.80942618934640D-15/
       DATA X50( 9)/.54357860872249D+01/,W50( 9)/.54659440318148D-13/
       DATA X50(10)/.50691175849172D+01/,W50(10)/.25066555238993D-11/
       DATA X50(11)/.47112936661690D+01/,W50(11)/.81118773649292D-10/
       DATA X50(12)/.43609731604545D+01/,W50(12)/.19090405438116D-08/
       DATA X50(13)/.40170681728581D+01/,W50(13)/.33467934040209D-07/
       DATA X50(14)/.36786770625153D+01/,W50(14)/.44570299668173D-06/
       DATA X50(15)/.33450383139379D+01/,W50(15)/.45816827079549D-05/
       DATA X50(16)/.30154977695745D+01/,W50(16)/.36840190537803D-04/
       DATA X50(17)/.26894847022677D+01/,W50(17)/.23426989210923D-03/
       DATA X50(18)/.23664939042987D+01/,W50(18)/.11890117817495D-02/
       DATA X50(19)/.20460719686864D+01/,W50(19)/.48532638261714D-02/
       DATA X50(20)/.17278065475159D+01/,W50(20)/.16031941068410D-01/
       DATA X50(21)/.14113177548983D+01/,W50(21)/.43079159156760D-01/
       DATA X50(22)/.10962511289577D+01/,W50(22)/.94548935477074D-01/
       DATA X50(23)/.78227172955460D+00/,W50(23)/.17003245567714D-00/
       DATA X50(24)/.46905905667824D-00/,W50(24)/.25113085633197D-00/
       DATA X50(25)/.15630254688947D-00/,W50(25)/.30508512920436D-00/
C *-/

C     if(npt.lt.2.or.npt.gt.20)go to 20
      if(npt.lt.2.or.npt.gt.50)go to 20
      ic = npt*npt/4
      n2 = (npt+1)/2
      i1 = 1
      i2 = npt
      do i = 1,n2
        xpt(i1) = -x(ic)
        xpt(i2) = x(ic)
        wht(i1) = w(ic)
        wht(i2) = wht(i1)
        i1 = i1 + 1
        i2 = i2 - 1
        ic = ic + 1
        end do
      return
   20 WRITE (6,30) npt
      stop
   30 format(1x,'***** error, gauss-hermite requested npt = ',i6)
      end

C * *-
      subroutine mxoutd(a,n,m,ms)
      implicit double precision (a-h,o-z)
      dimension a(1),b(11)

C * matrix printing routine.

      data isp/1/
      data ipos/132/
      data lins/59/
      j=1

C * write heading

      nend=ipos/12-1
      lend=(lins/isp)-2
      ipage=1
   10 lstrt=1
   20 continue
      jnt=j+nend-1
      ipage=ipage+1
      if(jnt.le.m) go to 50
      jnt=m
   50 continue
      WRITE (6,270)(jcur,jcur=j,jnt)
      ltend=lstrt+lend-1
      do l=lstrt,ltend
c
C * form output row line
c
      do k=1,nend
      kk=k
c
C * determine index
c
      jt=j+k-1
      if(ms.gt.1) go to 120
      if(ms.eq.1) go to 90
      ijnt=n*(jt-1)+l
      go to 140
   90 if(l.ge.jt) go to 110
      ijnt=l+(jt*jt-jt)/2
      go to 140
  110 ijnt=jt+(l*l-l)/2
      go to 140
  120 ijnt=0
      if(l.ne.jt) go to 140
      ijnt=l
  140 b(k)=0.0
      if(ijnt.le.0) go to 160
      b(k)=a(ijnt)
  160 continue
c
C * check if last column. if yes branch
c
      if(jt.ge.m) go to 180
C * end k loop
      enddo
c
C * end of line, now write
c
  180 if(isp.gt.1) go to 200
      WRITE (6,290)l,(b(jw),jw=1,kk)
      go to 210
  200 WRITE (6,300)l,(b(jw),jw=1,kk)
c
C * if end of rows, go check columns
  210 if(n.le.l) go to 230
c
C * end l loop
      enddo
c
C * end of page, now check for more output
c
      lstrt=lstrt+lend
      go to 20
c
C * end of column then return
c
  230 if(jt.ge.m) go to 250
      j=jt+1
      go to 10
  250 return
  270 format(1x,3x,'column ',10(3x,i3,6x))
  290 format(1x,1x,'row ',i3,2x,10(1pe12.4))
  300 format(1x,1x,'row ',i3,2x,10(1pe12.4))
      end
********************************************************************
      double precision function potbc(rvib)
      implicit double precision (a-h,o-z)
c
C * calculates the b-c potential.
c
      common/bcpar/parbc(20),idbc
      call potcf(rvib,parbc,idbc,v)
c
C * shift the zero-point energy of the diatom.
      potbc = v + parbc(2)
c
      return
      end
********************************************************************
      subroutine potcf(r,coef,ident,pot)
      implicit double precision (a-h,o-z)
c
C * written by g. a. parker   modified by a. feldt
C * evaluates a central field potential
C * on entering
C * r is the radial coordinate
C * coef contains the potential parameters
C * ident determines the potential type as follows:
C * ident=1 lennard-jones
C *       2 12-6-8
C *       3 morse
C *       4 morse-spline-van der waals
C *       5 buckingham-corner
C *       6 morse-morse-spline-van der waals
C *       7 eps*exp(-alpha*r-beta*r*r)
C *       8 simons-parr-finlan-dunham
C *       9 dunham
C *      10 n-6
C *      11 h-f-d
C *      12 morse-morse
C * on return
C * pot is the value of the potential a distance r.
C * in comments, x=r-rm
c
      dimension coef(20)
      if(ident.eq.1) go to 10
      if(ident.eq.2) go to 20
      if(ident.eq.3) go to 30
      if(ident.eq.4) go to 50
      if(ident.eq.5) go to 80
      if(ident.eq.6) go to 100
      if(ident.eq.7) go to 110
      if(ident.eq.8) go to 120
      if(ident.eq.9) go to 150
      if(ident.eq.10) go to 170
      if(ident.eq.11) go to 180
      if(ident.eq.12) go to 205
c
C * lennard-jones potential
C *   v= c/x**12 - 2c/x**6
c
   10 eps=coef(2)
      rm=coef(1)
      redr=r/rm
      r6i=(1.0/redr)**6
      pot=eps*(r6i-2.0)*r6i
      return
c
C * 12-6-8 potential
C *   v = c1/x**12 - 2c1/x**6 - c2/x**8
c
   20 eps=coef(2)
      rm=coef(1)
      c8=coef(3)
      redr=r/rm
      r8i=(1.0/r)**8.0
      red6i=(1.0/redr)**6.0
      pot=eps*red6i*(red6i-2.0)-c8*r8i
      return
c
C * morse potential
C *   v = d( exp(-2gx) -2exp(-gx) )
c
   30 eps=coef(2)
      rm=coef(1)
      beta=coef(3)
   40 expx=exp(-beta*(r-rm)/rm)
      pot=eps*expx*(expx-2.0)
      return
c
C * morse-spline-van der waals
C *   v = d( exp(-2gx) - 2exp(-gx) )
C *     = c3(x**3) + c2(x**2) + c1(x) + c
C *     = -a/x**6 - b/x**8
c
   50 eps=coef(2)
      rm=coef(1)
      beta=coef(3)
      c6=coef(4)
      c8=coef(5)
      rho=r/rm
      rho2=1.6
      rho1=1.d0+0.693147180559945d0/beta
      if(rho.ge.rho2)go to 70
      if(rho.ge.rho1)go to 60
c
C * morse branch of the msv or m2sv potential
c
      go to 40
c
C * spline branch of the msv or m2sv potential
c
   60 c=(c6+c8/(1.6*rm)**2)/(1.6*rm)**6.0
      cp=-(c6*6.0+8.0*c8/(1.6*rm)**2)/(1.6*rm)**7.0*rm
      del1=rho-rho1
      del2=rho2-rho
      exp1=exp(-beta*(rho1-1.0))
      cf1=eps*exp1*(exp1-2.0)
      df1=-2.0*beta*eps*exp1*(exp1-1.0)
      cf2=-c
      df2=-cp
      rho21=rho2-rho1
      a=(3.0*cf1+df1*rho21)/rho21
      b=(3.0*cf2-df2*rho21)/rho21
      s3=(2.0*a-b)/3.0
      s4=(2.0*b-a)/3.0
      s1=(cf1-s3*rho21)/(rho21**3.0)
      s2=(cf2-s4*rho21)/(rho21**3.0)
      pot=del2*(s1*del2*del2+s3)+del1*(s2*del1*del1+s4)
      return
c
C * van der waals branch of the msv or m2sv potential
c
   70 r2i=(1.0/r)**2
      r6i=r2i**3
      pot=-r6i*(c6+c8*r2i)
      return
c
C * buckingham-corner with c8 as a variable
C *   see mark kiel's thesis page 18
c
   80 eps=coef(2)
      rm=coef(1)
      beta=4.0
      c6=coef(4)/(eps*rm**6.0)
      c8=coef(3)
      b=(6.0*c6+8.0*c8)/(c6+c8-1.0)
      a=(8.0-2.0*c6)*exp(b)/(b-8.0)
      redr=r/rm
      rsqi=1.0/(redr*redr)
      r6i=rsqi**3
      if(redr.gt.1.0)go to 90
      rho3i=(abs(1.0/redr-1.0))**3
      if(redr.gt.1.0)rho3i=-rho3i
      pot=eps*(a*exp(-b*redr)-exp(-beta*(rho3i))*(c6+c8*rsqi)*r6i)
      return
c
C * r.gt.rm branch for the buckingham-corner potential
c
   90 dexpbr=exp(-b*redr)
      pot=eps*(a*dexpbr-(c6+c8*rsqi)*r6i)
      return
c
C * morse-morse-spline-van der waals.
C *   v = d1( exp(-g1x) - 2exp(-g1x) )
C *     = d2( exp(-g2x) - 2exp(-g2x) )
C *     = c3(x**3) + c2(x**2) + c1(x) + c
C *     = - a/x**6 - b/x**8
c
  100 eps=coef(2)
      rm=coef(1)
      beta=coef(3)
      c6=coef(4)
      c8=coef(5)
      rho=r/rm
      rho2=1.6
      rho1=1.d0+0.693147180559945d0/beta
      if(rho.ge.rho2)go to 70
      if(rho.ge.rho1)go to 60
      sigma=1.0+alog(0.5)/beta
      if(rho.gt.sigma)go to 40
      betap=coef(6)
      rmp=sigma*rm/(1.0+alog(0.5)/betap)
      epsp=eps*beta*rmp/(betap*rm)
      beta=betap
      rm=rmp
      eps=epsp
      if(rho.ne.0.0)go to 40
      return
c
C * pot=eps*exp(-alpha*r-beta*r*r)
c
  110 eps=coef(2)
      alpha=coef(1)
      beta=coef(3)
      pot=exp(eps-r*(alpha+beta*r))
      return
c
C * simons-parr-finlan-dunham
C *   v = sum( cn(x/r)**n )
c
  120 rm=coef(1)
      eps=coef(2)
      b0=coef(3)
      pot=1.d0
      rexp=(r-rm)/r
      n=coef(4)+0.005
      rn=rexp
      do i=5,n
        pot=pot+coef(i)*rn
        rn=rn*rexp
        end do
      bsq=rexp*rexp
      pot=pot*b0*bsq
      pot=pot-eps
      return
c
C * dunham
C *   v = sum( cn(x/rm)**n )
c
  150 rm=coef(1)
      eps=coef(2)
      b0=coef(3)
      n=coef(4)+3.005
      rexp=(r-rm)/rm
      rmi=1.d0/rm
      s1=1.d0
      if(n.lt.5)go to 161
      r1=1.d0
      do i=5,n
        r1=r1*rexp
        s1=s1+coef(i)*r1
        end do
  161 c1=b0*rexp*rexp
      pot=c1*s1-eps
      return
c
C * n-6
C *   v = a( (rm/r)**n - b(rm/r)**6 )
c
  170 rm=coef(1)
      eps=coef(2)
      alpha=coef(3)
      alpha6=alpha/6.d0
      rho=r/rm
      rhoalf=rho**alpha
      rho6=rho**6.d0
      pot=eps*(rhoalf-rho6*alpha6)
      return
c
C * h-f-d
C *   see jcp 71(6) 2637
c
  180 rm=coef(1)
      a=coef(2)
      alfa=coef(3)
      c6=coef(4)
      c8=coef(5)
      c10=coef(6)
      x=r/rm
      if(x.ge.1.28d0)go to 190
      b=1.28*rm
      c=b/r-1.d0
      d=exp(-c*c)
      rsqi=1.d0/(r*r)
      fac=2.d0*b*c*rsqi
      dd=d*fac
      ddd=d*(fac*fac-(6.d0*b/r-4.d0)*b*rsqi/r)
      go to 200
  190 d=1.d0
      dd=0.d0
      ddd=0.d0
  200 sr=a*exp(-alfa*r)
      dsr=-alfa*sr
      ddsr=-alfa*dsr
      rsqi=1.d0/(r*r)
      r6i=rsqi*rsqi*rsqi
      vdw=-(c6+(c8+c10*rsqi)*rsqi)*r6i
      pot=sr+vdw*d
      return
c
C * morse-morse.
C *   v = d1( exp(-g1x) - 2exp(-g1x) )
C *     = d2( exp(-g2x) - 2exp(-g2x) )
c
  205 eps=coef(2)
      rm=coef(1)
      beta=coef(3)
      betap=coef(4)
      rho=r/rm
      sigma=1.0+alog(0.5)/beta
      if(rho.gt.sigma)go to 40
      rmp=sigma*rm/(1.0+alog(0.5)/betap)
      epsp=eps*beta*rmp/(betap*rm)
      beta=betap
      rm=rmp
      eps=epsp
      if(rho.ne.0.0)go to 40
      return
      end

C *-
      subroutine potint
      implicit double precision (a-h,o-z)

C * reads in data to intialize diatom potentials.

      dimension ivlim(20)
      common/bcpar/parbc(20),idbc
C *   namelist/potbc/parbc,idbc

C * set default data

      do i=1,20
         parbc(i)=0.d0
      end do
      ivlim(1)=2
      ivlim(2)=3
      ivlim(3)=3
      ivlim(4)=5
      ivlim(5)=4
      ivlim(6)=6
      ivlim(7)=3
      ivlim(8)=4
      ivlim(9)=10
      ivlim(10)=3
      ivlim(11)=6
      ivlim(12)=4

C * read in data for the bc diatom
C *  idbc: determines which potential is used 
C *        see the corresponding format statement for the computed goto below
C * parbc: holds the parameters which determine these potentials

C     READ (5,3001) idbc
      READ (5,*)    idbc
      WRITE (6,330)
      lim = ivlim(idbc)

C     READ (5,3011) (parbc(i),i=1,lim)
      READ (5,*)    (parbc(i),i=1,lim)

      if (idbc .eq. 8) then
         nterm = parbc(4) + .0005
C        READ (5,3011) (parbc(i),i=5,nterm)
         READ (5,*)    (parbc(i),i=5,nterm)
      endif
      
 3001 format(5i5)
 3011 format(4e18.8)
      if(idbc.eq.1) go to 40
      if(idbc.eq.2) go to 50
      if(idbc.eq.3) go to 60
      if(idbc.eq.4) go to 70
      if(idbc.eq.5) go to 80
      if(idbc.eq.6) go to 90
      if(idbc.eq.7) go to 100
      if(idbc.eq.8) go to 110
      if(idbc.eq.9) go to 120
      if(idbc.eq.10) go to 130
      if(idbc.eq.11) go to 140
      if(idbc.eq.12) go to 145
   40 continue
      WRITE (6,350)
      WRITE (6,2460)
      WRITE (6,2620) (parbc(i),i=1,lim)
      go to 150
   50 continue
      WRITE (6,360)
      WRITE (6,2470)
      WRITE (6,2620) (parbc(i),i=1,lim)
      go to 150
   60 continue
      WRITE (6,370)
      WRITE (6,2480)
      WRITE (6,2620) (parbc(i),i=1,lim)
      go to 150
   70 continue
      WRITE (6,380)
      WRITE (6,2490)
      WRITE (6,2620) (parbc(i),i=1,lim)
      go to 150
   80 continue
      WRITE (6,390)
      WRITE (6,2500)
      WRITE (6,2620) (parbc(i),i=1,lim)
      go to 150
   90 continue
      WRITE (6,400)
      WRITE (6,2510)
      WRITE (6,2620) (parbc(i),i=1,lim)
      go to 150
  100 continue
      WRITE (6,410)
      WRITE (6,2520)
      WRITE (6,2620) (parbc(i),i=1,lim)
      go to 150
  110 continue
      WRITE (6,420)
      WRITE (6,2530) (parbc(i),i=1,nterm)
C     WRITE (6,2640) (parbc(i),i=1,nterm)
      go to 150
  120 continue
      WRITE (6,430)
      WRITE (6,2540)
      WRITE (6,2640) (parbc(i),i=1,lim)
      go to 150
  130 continue
      WRITE (6,440)
      WRITE (6,2550)
      WRITE (6,2620) (parbc(i),i=1,lim)
      go to 150
  140 continue
      WRITE (6,450)
      WRITE (6,2560)
      WRITE (6,2620) (parbc(i),i=1,lim)
      go to 150
  145 continue
      WRITE (6,455)
      WRITE (6,2565)
      WRITE (6,2620) (parbc(i),i=1,lim)
  150 continue
      return

C * format statements.

  330 format(1x,'parameters for the bc diatom')
  350 format(1x,'lennard-jones parameters')
  360 format(1x,'12-6-8 parameters')
  370 format(1x,'morse parameters')
  380 format(1x,'morse-spline-van der waals parameters')
  390 format(1x,'buckingham-corner parameters')
  400 format(1x,'morse-morse-spline-van der waals parameters')
  410 format(1x,'repulsive exponential parameters')
  420 format(1x,'simons-parr-finlan-dunham parameters')
  430 format(1x,'dunham parameters')
  440 format(1x,'n-6 parameters')
  450 format(1x,'h-f-d parameters')
  455 format(1x,'morse-morse parameters')
 2460 format(1x,5x,'rmin',8x,'epsilon')
 2470 format(1x,5x,'rmin',8x,'epsilon',10x,'c8')
 2480 format(1x,5x,'rmin',8x,'epsilon',9x,'beta')
 2490 format(1x,5x,'rmin',8x,'epsilon',9x,
     >'beta',11x,'c6',12x,'c8')
 2500 format(1x,5x,'rmin',8x,'epsilon',10x,'c6',
     >12x,'c8')
 2510 format(1x,5x,'rmin',8x,'epsilon',9x,'beta',11x,
     1'c6',12x,'c8',11x,'betap')
 2520 format(1x,4x,'alpha',8x,'epsilon',9x,'beta')
 2530 format(5x,'rmin = ',1pe15.8,/2x,'epsilon = ',1pe15.8,/1x,
     $     'expansion coefficients = ',1pe15.8,7(/26x,1pe15.8),/)
 2540 format(1x,5x,'rmin',8x,'epsilon',11x,
     >'expansion coefficients')
 2550 format(1x,5x,'rmin',8x,'epsilon',8x,'alpha')
 2560 format(1x,5x,'rmin',11x,'a',12x,'alpha',10x,
     1'c6',12x,'c8',11x,'c10')
 2565 format(1x,5x,'rmin',8x,'epsilon',9x,'beta',8x,'betap')
 2620 format(1x,6(1pe14.7))
 2640 format(1x,10(1pe13.6))
      end

C *-
      subroutine tdiag(v, e, t, subd, n, ntdim)
      implicit double precision (a-h,o-z)

C *   <tdiag>--calls tred3, tql2, and trbak3 to diagonalize a matrix
C *   v(n*(n+1)/2) -- matrix to be diagonalized, must be symmetric
C *   e(n) -- output, the eigenvalues of <v>
C *   t(n,n) -- output, the eigenvectors of <v>
C *   n  -- appropriate dimension variables   (integer*4)
C *   subd(n) -- scratch array used by tred2 and tql2
c
   10 format ( ' diagonalization error #', i3, ' returned by tql2' )
      dimension e(n), t(ntdim,ntdim), subd(n), v(1)
      double precision v, e, t, subd
      nv=(n*(n+1))/2
      call tred3(n,nv,v,e,subd,subd)
      do i=1,n
        do j=1,n
          t(j,i)=0.d0
          end do
        t(i,i)=1.d0
        end do
      call tql2(ntdim,n,e,subd,t,ierr)
      if(ierr.ne.0) go to 40
      call trbak3(ntdim,n,nv,v,n,t)
      return
   40 continue
      if ( ierr .ne. 0 ) print 10, ierr
      return
      end
      
C *-
      subroutine tql2(nm,n,d,e,z,ierr)
      implicit double precision (a-h,o-z)
      integer i,j,k,l,m,n,ii,l1,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision b,c,f,g,h,p,r,s,machep

C * machep is a machine dependent parameter specifying
C * the relative precision of floating point arithmetic.
C * machep = 16.0d0**(-13) for long form arithmetic
C * on s360

      machep=1.0e-14
      ierr = 0
      if (n .eq. 1) go to 160
      do i = 2, n
         e(i-1) = e(i)
      end do
      f = 0.0d0
      b = 0.0d0
      e(n) = 0.0d0
      do l = 1, n
         j = 0
         h = machep * (abs(d(l)) + abs(e(l)))
         if (b .lt. h) b = h

C * look for small sub-diagonal element

         do m = l, n
            if (abs(e(m)) .le. b) go to 30

C * e(n) is always zero, so there is no exit
C * through the bottom of the loop

         end do
 30      if (m .eq. l) go to 100
 40      if (j .eq. 30) go to 150
         j = j + 1

C * form shift

         l1 = l + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = sqrt(p*p+1.0d0)
         d(l) = e(l) / (p + sign(r,p))
         h = g - d(l)
         do i = l1, n
            d(i) = d(i) - h
         end do
         f = f + h

C * ql transformation

         p = d(m)
         c = 1.0d0
         s = 0.0d0
         mml = m - l

C * for i=m-1 step -1 until l do --

         do ii = 1, mml
            i = m - ii
            g = c * e(i)
            h = c * p
            if (abs(p) .lt. abs(e(i))) go to 60
            c = e(i) / p
            r = sqrt(c*c+1.0d0)
            e(i+1) = s * p * r
            s = c / r
            c = 1.0d0 / r
            go to 70
 60         c = p / e(i)
            r = sqrt(c*c+1.0d0)
            e(i+1) = s * e(i) * r
            s = 1.0d0 / r
            c = c * s
 70         p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))

C * form vector
            do k = 1, n

               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
            end do
C * end ii loop
         enddo
         e(l) = s * p
         d(l) = c * p
         if (abs(e(l)) .gt. b) go to 40
 100     d(l) = d(l) + f
C * end l loop
      enddo

C * order eigenvalues and eigenvectors

      do ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
         do j = ii, n
            if (d(j) .ge. p) go to 120
            k = j
            p = d(j)
 120        continue
         end do
         if (k .eq. i) go to 140
         d(k) = d(i)
         d(i) = p
         do j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
         end do
 140     continue
      end do
      go to 160

C * set error -- no convergence to an
C * eigenvalue after 30 iterations

 150  ierr = l
 160  return
      end

C *-
      subroutine trbak3(nm,n,nv,a,m,z)
      implicit double precision (a-h,o-z)
      integer i,j,k,l,m,n,ik,iz,nm,nv
      double precision a(nv),z(nm,m)
      double precision h,s

C *     this routine is a translation of the algol procedure trbak3,
C *     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
C *     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).

C *     this routine forms the eigenvectors of a double precision symmetric
C *     matrix by back transforming those of the corresponding
C *     symmetric tridiagonal matrix determined bytred3.
C *     on input:
C *     nm must be set to the row dimension of two-dimensional
C *     array parameters as declared in the calling program
C *     dimension statement;
C *     n is the order of the matrix;
C *     nv must be set to the dimension of the array parameter a
C *     as declared in the calling program dimension statement;
C *     a contains information about the orthogonal transformations
C *     used in the reduction bytred3in its first
C *     n*(n+1)/2 positions;
C *     m is the number of eigenvectors to be back transformed;
C *     z contains the eigenvectors to be back transformed
C *     in its first m columns.
C *     on output:
C *     z contains the transformed eigenvectors
C *     in its first m columns.
C *     note that trbak3 preserves vector euclidean norms.
C *     questions and comments should be directed to b. s. garbow,
C *     applied mathematics division, argonne national laboratory

      if (m .eq. 0) go to 50
      if (n .eq. 1) go to 50
      do i = 2, n
         l = i - 1
         iz = (i * l) / 2
         ik = iz + i
         h = a(ik)
         if (h .eq. 0.0d0) go to 40
         do j = 1, m
            s = 0.0d0
            ik = iz
            do k = 1, l
               ik = ik + 1
               s = s + a(ik) * z(k,j)
            end do
            s = (s / h) / h
            ik = iz
            do k = 1, l
               ik = ik + 1
               z(k,j) = z(k,j) - s * a(ik)
            end do
         end do
 40      continue
      end do
 50   return
      end

C *-
      subroutine tred3(n,nv,a,d,e,e2)
      implicit double precision (a-h,o-z)
      integer i,j,k,l,n,ii,iz,jk,nv
      double precision a(nv),d(n),e(n),e2(n)
      double precision f,g,h,hh,scale

C * this routine is a translation of the algol procedure tred3,
C * num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
C * handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
C * this routine reduces a double precision symmetric matrix, stored as
C * a one-dimensional array, to a symmetric tridiagonal matrix
C * using orthogonal similarity transformations.
C * on input:
C * n is the order of the matrix;
C * nv must be set to the dimension of the array parameter a
C * as declared in the calling program dimension statement;
C * a contains the lower triangle of the double precision symmetric
C * input matrix, stored row-wise as a one-dimensional
C * array, in its first n*(n+1)/2 positions.
C * on output:
C * a contains information about the orthogonal
C * transformations used in the reduction;
C * d contains the diagonal elements of the tridiagonal matrix;
C * e contains the subdiagonal elements of the tridiagonal
C * matrix in its last n-1 positions.e(1) is set to zero;
C * e2 contains the squares of the corresponding elements of e.
C * e2 may coincide with e if the squares are not needed.
C * questions and comments should be directed to b. s. garbow,
C * applied mathematics division, argonne national laboratory

      do ii = 1, n
         i = n + 1 - ii
         l = i - 1
         iz = (i * l) / 2
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 1) go to 20
         do k = 1, l
            iz = iz + 1
            d(k) = a(iz)
            scale = scale + abs(d(k))
         end do
         if (scale .ne. 0.0d0) go to 30
 20      e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 80
 30      do k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
         end do
         e2(i) = scale * scale * h
         f = d(l)
         g = -sign(sqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
         a(iz) = scale * d(l)
         if (l .eq. 1) go to 80
         f = 0.0d0
         do j = 1, l
            g = 0.0d0
            jk = (j * (j-1)) / 2
            do k = 1, l
               jk = jk + 1
               if (k .gt. j) jk = jk + k - 2
               g = g + a(jk) * d(k)
            end do
            e(j) = g / h
            f = f + e(j) * d(j)
         end do
         hh = f / (h + h)
         jk = 0
         do j = 1, l
            f = d(j)
            g = e(j) - hh * f
            e(j) = g
            do k = 1, j
               jk = jk + 1
               a(jk) = a(jk) - f * e(k) - g * d(k)
            end do
         end do
 80      d(i) = a(iz+1)
         a(iz+1) = scale * sqrt(h)
C * end ii loop
      enddo
      return
      end

C *-
      subroutine vibfun(noscil,jang,ngaush,haml,vect,chi,nulow,nuhigh
     $     ,nudim,eint,alpha,rx,z,w,re,e,h)

      implicit double precision (a-h,o-z)

C * this routine determines the vibrational eigenfunctions and
C * eigenvalues of a diatomic molecule in a fixed rotational
C * state using a harmonic oscillator basis set to diagonalize
C * the hamiltonian.

C * ON ENTERING
C * noscil: the number of harmonic oscilators to be used
C *   jang: the orbital angular momentum quantum of the diatomic molecule.
C * ngaush: the number of gauss hermite integration points.
C *  nulow: label the minimun vibrational levels
C * nuhigh: label the maximun vibrational levels
C *     rx: the ratio of the harmonic oscillator position to the equilibrium position of the diatomic molecule.
C *  alpha:

C * ON RETURN
C * vect: contains the eigenvectors.
C *  chi: contains the unweighted eigenfunctions.
C *    z: gauss-hermite integration points.
C *    w: gauss-hermite weights.

C * this program was written by g. a. parker.  and adapted by a. feldt

      logical iprint,sho,outvib

      dimension vect(20,20), chi(ngaush,nudim), dumvec(400)

      integer Nvib, Nrw
      real*8  rr(5000,50), vw(5000,50), rsca(50)
      real*8  rtemp(5000), vtemp(5000), vsca(50), csplin(5000)

      dimension eint(1), e(1), h(1)
      dimension haml(1), w(1), z(1)
      common/diatom/we,wexe,be,alfe,de
      common/prntop/ iprint(23)
      common/masses/bmass,cmass,ubc,uml,umr

      common/readvib/outvib,lvib

C *  iprint(21): prints alpha,rx,noscil,ngaush, the lowest
C *              and highest nu's, hs, and hsinv
C *  iprint(22): prints jang, and the vibrational energies
C *  iprint(23): prints the eigenvectors, the non-weighted
C *              wavefunctions, and the weighted wavefunctions

      itime = 0
      nvec = 20
      nvec = 20
      cminv = 2.1947462e5
      sho = iprint(1)

C * determine zeros and weights

      call  hermit(ngaush,z,w)
      sqrj = jang*(jang+1)
      covms = 2.d0*ubc
      hs = 1.0d0/(covms*re*re)

C * calculate experimental energies for a comparison.

      do np1 = nulow, nuhigh
         vph = float(np1) - 0.5d0
         e(np1) = (we*vph - wexe*vph*vph + be*sqrj - de*sqrj*sqrj - alfe
     $        *vph*sqrj)
      end do
      alphas = alpha*alpha
      alphai = 1.0d0/alpha
      hsinv = 1.0d0/hs

C * if desired print the input parameters.

      if (.not.iprint(21)) go to 60
      if (itime.ne.0) go to 50
      itime = 1
      write (6,280) alpha,rx
      write (6,290) noscil
      write (6,300) ngaush
      mprnt = nulow - 1
      nprnt = nuhigh - 1
      write (6,310) mprnt,nprnt
      write (6,200) hs,hsinv

 50   if (iprint(22).or.iprint(23)) write (6,270) jang
      if (iprint(22)) write (6,320)
 60   continue
      if(sho) go to 81

C * initialize the hamiltonian matrix to zero

      npos = noscil*(noscil+1)/2
      do i = 1, npos
         haml(i) = 0.0d0
      end do

C * Preform the integrations
C * call dhep to calculate the hermite polynomials h(i+1) for node z.
C * defines hamiltonian for spfd potential potbc which supplies the
C * simons-parr-finlan-dunham (or morse) potential for ab molecule.

      do i = 1,ngaush
         call  dhep(h,z(i),noscil)
         r = alphai*z(i) + rx
         x = r - 1.0d0
         rre = r*re
         cterm = -alphas*z(i)**2 + sqrj/r**2 + potbc(rre)*hsinv
         wi = w(i)
         nnp = 0
         do n = 1,noscil
            term = wi*(float(2*n-1)*alphas + cterm)*h(n)
            do np = 1,n
               nnp = nnp + 1
               haml(nnp) = haml(nnp) + h(np)*term
            end do
         end do
      end do

C * diagonalize the hamiltonian to determine the eigenvectors.

      call  tdiag(haml,eint,vect,h,noscil,nvec)
      do j = nulow, nuhigh
         if(vect(j,j).gt.0.e+00) go to 84
         do k = 1,noscil
            vect(k,j) = -vect(k,j)
         end do
 84      continue
      end do
 81   continue

C * Determine the unweighted eigenfunctions at the gaussian
C * integration points.
C * Store vib. functions in array CHI.
C *   If SHO is true, use SHO (-> H(z) ) as vib. functions. These Hermite
C * functions H(z) are calculated at the nodes z of Hermite polynomials,
C * so they are NOT scaled and unweighted.

      nos = noscil
      if(.not.sho) go to 89
      nos = nuhigh
      do i = 1,noscil
         vph = float(i-1) + nulow - 0.5e+00
         vect(i,i) = 1.e+00
         eint(i) = we*vph/hs
      end do
 89   continue

C * find negative energies and skip them
      do iskip = nulow, nuhigh
         if (eint(iskip).gt.0d0) goto 102
      end do
 102  iskip = iskip - nulow
      if(iskip.gt.0) WRITE (6,*) ' ** warning - ',iskip
     $     ,' energies were negative and skipped'

C     if(iskip.gt.noscil-(nuhigh-nulow+1)) then
      if(iskip.gt.0) then
         WRITE (6,*) '** error - too many energies were negative to',
     $        ' recover run....stop...increase noscil'
         stop
      endif

      do i = 1, ngaush
         call  dhep(h,z(i),nos)
         do j = nulow+iskip, nuhigh+iskip
            jnew = j - nulow - iskip + 1
            if(sho) go to 95
            sum = 0.0d0
            do k = 1, noscil
               sum = sum + vect(k,j)*h(k)
            end do
            go to 100
 95         sum = h(j)
 100        chi(i,jnew) = sum
         end do
      end do

C *-
C * For NUMERICAL vib. wavefunctions only. READ them from outer files.
C *-/
      if (outvib) then
         write(6,*)
         write(6,*) " Read vibrational wavefunctions from outer files:"

C * Read raw vibrational wavefunctions data and eigenenergies

         open(lvib, file='vibfun', status='old', form ='formatted')
         read (lvib,*) Nvib, Nrw
         do i = 1, Nvib
            do j = 1, Nrw
               read (lvib,*) rr(j,i), vw(j,i)
            enddo
         enddo

C        write (6,355)
         do i = 1, Nvib
            read (lvib,*) eint(i)
C           write (6,365) i, eint(i)
         enddo

C * Interpolate to get the unweighted eigenfunctions

C        write (6,370)
         do i = 1, NGAUSH
            rsca(i) = (z(i)/alpha+rx)*re
C           write (6,380) i, rsca(i)
         enddo

         do i = 1, Nvib
            do j = 1, Nrw
               rtemp(j) = rr(j,i)
               vtemp(j) = vw(j,i)
            enddo
            call spline(rtemp,vtemp,Nrw,1.d30,1.d30,csplin)
            do j = 1, ngaush
               call splint(rtemp,vtemp,csplin,Nrw,rsca(j),vsca(j))
               chi(j,i) = vsca(j)
            enddo
         enddo

         write(6,*) " Check the following r's !"
         write(6,*)
         write(6,*) "   The first and last r from reading are: ",rr(1,1)
     $        ,rr(Nrw,1)
         write(6,*) " The first and last r needed in wlam are: ",rsca(1)
     $        ,rsca(NGAUSH)
         write(6,*) 

      endif

C * if desired compare with the experimental energies.

      if (.not.iprint(22)) go to 130
      write (6,330)
 130  continue
      do i = nulow, nuhigh
         im1 = i - 1
         if (outvib) then
            einti = eint(i)
            eintwn = einti
         else
            einti = eint(i+iskip)*hs
         endif
         err = e(i) - einti
         if(e(i).eq.0.d0) perr = 0.d0
         if(e(i).ne.0.d0) perr = 100.0d0*err/e(i)

         eintwn = einti*cminv
         if (.not.iprint(22)) go to 140
         write (6,210) im1,eintwn,einti,e(i),err,perr
 140     continue
         inew = i - nulow + 1
         eint(inew) = einti
      end do

C * if desired print eigenvectors and eigenfunctions.

      if (.not.iprint(23)) go to 180
      ms = 0
      write (6,260)
      index = 0
      do i = 1,noscil
         do j = 1,noscil
            index = index + 1
            dumvec(index) = vect(j,i)
         end do
      end do

C     call  mxoutd(dumvec,noscil,noscil,ms)
      if (.not. outvib) call  mxoutd(dumvec,noscil,noscil,ms)

      write (6,220)
      do k = 1, ngaush
        x = alphai*z(k) + rx - 1.0d0
        r = (1.d0+x)*re
        write (6,230) z(k),x,r,w(k)
        end do
      write (6,250)
      call  mxoutd(chi,ngaush,nudim,ms)
      index = 0
      do i = 1, nudim
        do j = 1, ngaush
          index = index + 1
          haml(index) = chi(j,i)*sqrt(w(j))
          end do
        end do
      write (6,240)
      call  mxoutd(haml,ngaush,nudim,ms)
  180 continue
      return
c
C * formats.
c
  200 format(/4x,'hs =',e14.7,/1x,'hsinv =',e14.7)
  210 format(1x,i2,3x,5f15.8)
  220 format(1x,'      z                  x                  r'
     1,'                weight')
  230 format(1x,e14.7,5x,e14.7,5x,e14.7,5x,e14.7)
  240 format(/1x,'weighted wavefunctions')
  250 format(/1x,'non-weighted wavefunctions')
  260 format(/1x,'vibrational part of vib-rot eigenvectors')
  270 format(/1x,'j =',i3)
  280 format(/1x,'alpha =',f10.5,/4x,'rx =',f10.5)
  290 format(/1x,'number of oscillator basis functions =',i5)
  300 format(1x,'number of integration points =',i5)
  310 format(1x,'vibrational states are from nu =',i3,3x,'to nu =',i3)
  320 format(/1x,'the vibrational energies:')
  330 format(/1x,'nu',6x,'E-calc(w.n.)',3x,'  E-calc(a.u.)',2x,
     1'E-Herzberg(a.u.)',2x,'error',4x,'    percent error')
  355 format(/2x
     $     ,'Eigenvalues of vibrational functions (read from file):',/5x
     $     ,'  i    E(i)  ')
  365 format(5x,i3,1x,1pe16.8)
  370 format(/5x,'The r quadrature points used are ',
     #/5x,'(Interpolated from outer file): ',
     #/5x,'   k         r(k)  ')
  380 format(5x,i4,5x,f9.5)

      end

C *-
      double precision function xil(r)
      implicit double precision (a-h,o-z)

C *  function for electron-nuclear pot.: relc .gt. rn.

      common /silly/ vect(20,20),noscil
      common/vvp/iv,ivp,rss(15),vrs(15),ndo

      common /need/ rx,alpha,re
 
      common /vnpar/znuc,umol,relc,lamb
      common/store/ h(20),xx(21),num(21)
      ier=0
      rn=umol*r
      v=-znuc*((rn/relc)**lamb)/relc
      z=alpha*(r/re-rx)
      call dhep(h,z,noscil)
      chi1=0.e+00
      chi2=0.e+00
      do i=1,noscil
        chi1=chi1+h(i)*vect(i,iv)
        chi2=chi2+h(i)*vect(i,ivp)
        end do
      xil=chi1*chi2*v*exp(-z*z)
      return
      end

C *-
      double precision function xir(r)
      implicit double precision (a-h,o-z)

C *  function for electron-nuclear pot.: relc .lt. rn.

      common /silly/ vect(20,20),noscil
      common/vvp/iv,ivp,rss(15),vrs(15),ndo

      common /need/ rx,alpha,re
 
      common /vnpar/znuc,umol,relc,lamb
      common/store/ h(20),xx(21),num(21)
      ier=0
      rn=umol*r
      v=-znuc*((relc/rn)**lamb)/rn
      z=alpha*(r/re-rx)
      call dhep(h,z,noscil)
      chi1=0.e+00
      chi2=0.e+00
      do i=1,noscil
        chi1=chi1+h(i)*vect(i,iv)
        chi2=chi2+h(i)*vect(i,ivp)
        end do
      xir=chi1*chi2*v*exp(-z*z)
      return
      end

C *-
      subroutine chitst(rbeg,rend,ivmax)
      implicit double precision (a-h,o-z)
      common/silly/ vect(20,20),noscil
      common/store/ h(20),xx(21),num(21)

      common/need/ rx,alpha,re

      rdel=(rend-rbeg)/250.e+00
      r=rbeg
      do i=1,ivmax
        num(i)=i-1
        end do
      ist=ivmax+1
      do i=1,ivmax
        do j=i,ivmax
          num(ist)=10*(j-1)+i-1
          ist=ist+1
          end do
        end do
      iend=ist-1
      WRITE (9,1000) (num(i),i=1,iend)
   10 continue
      z=(r/re-rx)*alpha
      zex=exp(-z*z/2.0e+00)*sqrt(alpha/re)
      do i=1,ivmax
        xx(i)=0.0e+00
        end do
      call dhep(h,z,noscil)
      do j=1,ivmax
        do i=1,noscil
          xx(j)=xx(j)+h(i)*vect(i,j)
          end do
        xx(j)=xx(j)*zex
        end do
      ist=ivmax+1
      do i=1,ivmax
        do j=i,ivmax
          xx(ist)=xx(i)*xx(j)
          ist=ist+1
          end do
        end do
      iend=ist-1
      WRITE (9,1010) r,(xx(i),i=1,iend)
      r=r+rdel
      if(r.le.rend) go to 10
      return
 1000 format(/2x,'wavefunctions ',/6x,'r',7(5x,i5)/(7x,7i10))
 1010 format(f10.5,7(1pe10.3),(10x,7e10.3))
      end

C *-
      subroutine pottst(r,lam,rss,vrs,ndo,second,twice)
      implicit double precision (a-h,o-z)
      logical second,twice
      dimension rss(1),vrs(1),rs(102),vr(102)
      common/splin/ bpar(11),csplin(5000)
      if(.not.second) ir=1
      irs=ir
      rdel=0.025e+00
      rt=rss(1)
      do i=1,101
        rs(ir)=rt
        ir=ir+1
        ii=i
        rt=rt+rdel
        if(rt.gt.rss(ndo)+1.e-05) go to 20
        end do
   20 id=1
      if(second) id=irs
C * imsl v9
C *     call icsevu(rss,vrs,ndo,csplin,29,rs(id),vr(id),ii,ier)
C * imsl v10
C *      do 51 ix=id,ii
C *      vr(ix)=dcsval(rs(ix),ndo-1,bpar,csplin)
C *   51 continue
      do ix=id,ii
        call splint(rss,vrs,csplin,ndo,rs(ix),vr(ix))
        end do
      if(twice.and..not.second) return
      WRITE (9,1000) r, lam
      ir=ir-1
      WRITE (9,1010) (rs(i),vr(i),i=1,ir)
      return
 1000 format(/2x,'potential for rint =',f10.5,5x,'lamda =',i3)
 1010 format(f10.5,e23.16)
      end

C *-
      subroutine inttst(r,lam,rbeg,rend,stepr,ivm,second)
      implicit double precision (a-h,o-z)
      logical second
      common/store/ h(20),xx(21),num(21)
      common/vvp/iv,ivp,rss(15),vrs(15),ndo

      common/need/rx,alpha,re

      if(second) go to 5
      WRITE (9,1010) r, lam
      ivi=1
      do i=1,ivm
        do j=i,ivm
          num(ivi)=(j-1)*10+i-1
          ivi=ivi+1
          end do
        end do
      ivi=ivi-1
      WRITE (9,1020) (num(i),i=1,ivi)
    5 continue
      nofr=(rend-rbeg)/stepr
      step=(rend-rbeg)/float(nofr)
      if(.not.second) nofr=nofr+1
      WRITE (9,1040) nofr,step
      rt=rbeg
      if(second) rt=rt+step
      do i=1,nofr
        ivi=1
        do iv=1,ivm
          do ivp=iv,ivm
            xx(ivi)=xil(rt)*alpha/re
            ivi=ivi+1
            end do
          end do
        ive=ivi-1
        WRITE (9,1000) rt,(xx(j),j=1,ive)
        rt=rt+step
        end do
      return
 1000 format(f10.5,107e10.3/(10x,7e10.3))
 1010 format(/2x,'integrand for rint=',f10.5,5x,'lamda=',i3)
 1020 format(5x,'r',7i10,(6x,7i10))
 1040 format(i5,e23.16)
      end

C *-
      double precision function quanc8(ifun,a,b,relerr,abserr,errest
     $     ,nofun,flag)

      implicit double precision (a-h,o-z)

      double precision a,b,abserr,relerr,result,errest,flag
      integer nofun
      double precision w0,w1,w2,w3,w4,area,x0,f0,stone,step,cor11,temp
      double precision qprev,qnow,qdiff,qleft,esterr,tolerr
      double precision qright(31),f(16),x(16),fsave(8,30),xsave(8,30)
      integer levmin,levmax,levout,nomax,nofin,lev,nim,i,j

c---->stage 1 *** general initialization

      levmin=1
      levmax=30
      levout=6
      nomax=5000
      nofin=nomax-8*(levmax-levout+2**(levout+1))
      w0=3.956d03/1.4175d04
      w1=2.3552d04/1.4175d04
      w2=-3.712d03/1.4175d04
      w3=4.1984d04/1.4175d04
      w4=-1.816d04/1.4175d04
      flag=0.d00
      result=0.d00
      cor11=0.d00
      errest=0.d00
      area=0.d00
      nofun=0
      quanc8=result
      if(a.eq.b) return

c---->stage 2 *** initialization for first interval

      lev=0
      nim=1
      x0=a
      x(16)=b
      qprev=0.d00
      if (ifun.eq.0) then
        f0=xil(x0)
      else
        f0=xir(x0)
      endif
      stone=(b-a)/1.6d01
      x(8)=(x0+x(16))/2.d00
      x(4)=(x0+x(8))/2.d00
      x(12)=(x(8)+x(16))/2.d00
      x(2)=(x0+x(4))/2.d00
      x(6)=(x(4)+x(8))/2.d00
      x(10)=(x(8)+x(12))/2.d00
      x(14)=(x(12)+x(16))/2.d00
      do j=2,16,2
        if (ifun.eq.0) then
          f(j)=xil(x(j))
        else
          f(j)=xir(x(j))
        endif
C       f(j)=fun(x(j))
        end do
      nofun=9

c---->stage 3 *** central calculation

   30 x(1)=(x0+x(2))/2.d00
C     f(1)=fun(x(1))
      if (ifun.eq.0) then
        f(1)=xil(x(1))
      else
        f(1)=xir(x(1))
        endif
      do j=3,15,2
        x(j)=(x(j-1)+x(j+1))/2.d00
        if (ifun.eq.0) then
          f(j)=xil(x(j))
        else
          f(j)=xir(x(j))
        endif
C       f(j)=fun(x(j))
        end do
      nofun=nofun+8
      step=(x(16)-x0)/1.6d01
      qleft=(w0*(f0+f(8))+w1*(f(1)+f(7))+w2*(f(2)+f(6))
     1 +w3*(f(3)+f(5))+w4*f(4))*step
      qright(lev+1)=(w0*(f(8)+f(16))+w1*(f(9)+f(15))+w2*(f(10)+f(14))
     1 +w3*(f(11)+f(13))+w4*f(12))*step
      qnow=qleft+qright(lev+1)
      qdiff=qnow-qprev
      area=area+qdiff

c---->stage 4 *** interval convergence test

      esterr=abs(qdiff)/1.023d03
      tolerr=max(abserr,relerr*abs(area))*(step/stone)
      if(lev.lt.levmin) goto 50
      if(lev.ge.levmax) goto 62
      if(nofun.gt.nofin) goto 60
      if(esterr.le.tolerr) goto 70

c---->stage 5 *** no convergence

   50 nim=2*nim
      lev=lev+1

      do i=1,8
        fsave(i,lev)=f(i+8)
        xsave(i,lev)=x(i+8)
        end do
      qprev=qleft
      do i=1,8
        j=-i
        f(2*j+18)=f(j+9)
        x(2*j+18)=x(j+9)
        end do
      goto 30

c---->stage 6 *** trouble section

   60 nofin=2*nofin
      levmax=levout
      flag=flag+(b-x0)/(b-a)
      goto 70
   62 flag=flag+1.d00

c---->stage 7 *** interval converged

   70 result=result+qnow
      errest=errest+esterr
      cor11=cor11+qdiff/1.023d03
   72 if(nim.eq.2*(nim/2)) goto 75
      nim=nim/2
      lev=lev-1
      goto 72
   75 nim=nim+1
      if(lev.le.0) goto 80
      qprev=qright(lev)
      x0=x(16)
      f0=f(16)
      do i=1,8
        f(2*i)=fsave(i,lev)
        x(2*i)=xsave(i,lev)
        end do
      goto 30

c---->stage 8 *** finalize and return

   80 result=result+cor11
      quanc8=result
      if(errest.eq.0.d00) return
   82 temp=abs(result)+errest
      if(temp.ne.abs(result)) return
      errest=2.d00*errest
      go to 82
      end

C *-
      subroutine numfun(noff,ivm,nquad,fnum,fmat,chis,rquad,vquad,wquad
     $     ,minvib,ivibtr,ihoex)

      implicit double precision (A-H,O-Z)
      parameter (maxqd=20,maxho=20)
C * maxqd is max quadrature   maxho=max harmonic osc.
      dimension fin(35,21,11),fout(20,18,5,12)
      character*8 iflnm(20)
      integer doexch

      logical outvib

      dimension chis(nquad,ivm),fnum(1),rquad(1),fmat(ivm,ivm), vquad(1)
     $     ,wquad(1),nmunit(16),rss(15)
      dimension ho(maxqd,maxho),hox(maxho)
      common /splin/ bpar(11),csplin(5000)

      common/need/rx,alpha,re
      common/readvib/outvib,lvib
 
      zero = 0.0e+00
      if (outvib) then
         pref = re/alpha
      else
         pref = 1.0d0
      endif

      READ (5,*) doexch
      if (doexch .gt. 0) then
         WRITE (6,*) ' '
         WRITE (6,*) 'Calculating vibrational matrix elements of ',
     $        'the exchange potential.'
         WRITE (6,*) ' '
         call numfunex(ivm,nquad,fnum,chis,rquad,vquad,wquad,ivibtr
     $        ,ihoex)
         return
      endif

C * build harmonic oscillators if needed
      if (ihoex.gt.0) then
        do ix=1,nquad
          zold=alpha*(rquad(ix)/re-rx)
          call dhep(hox,zold,ihoex)
          do iy=1,ihoex
            ho(ix,iy)=hox(iy)
          enddo
        enddo
      endif
C *-/

C *  this subroutine vibrationally averages functions given
C *  at a set of internuclear separations.  the values can
C *  all be in one file or in the case of radial orbitals
C *  can be in nsepar files in unformatted form.

C *  input variables:

C *       nsepar     number of internuclear seps at which the numerical
C *                  double precision function exists.
C *       rss,fnum   the vectors holding the internuc. sep. values and
C *                  the values of the numerical function (e.g. the
C *                  polarizabilities given by poll and wolniewicz).
C *  nfile      if 1, file has values to be averaged which are
C *             read into fnum, or if .gt. 1, equals the number
C *             (nsepar) of files containing values to averaged.
C *  iprint     if nonzero, results printed on unit 6.
C *  ivibopt    if nonzero, do overlap of function and vib state
C *             if zero, do vib. avg. <v|f(R)|vp>
C *  norbi      no. of orbitals input.
C *  norbo      no. of orbitals output.
C *  nlexi      no. of l values input.
C *  nlexo      no. of l values output.
C *  iflnm      names of the input and output files (.le. 8 char.).
C *             there should be nfile + 1 names.

C     READ (5,9180) nsepar,nfile,iprint,ivibopt
      READ (5,*)    nsepar,nfile,iprint,ivibopt
C     READ (5,9180) norbi,norbo,nlexi,nlexo
      READ (5,*)    norbi,norbo,nlexi,nlexo
      write (6,*)
      WRITE (6,*) ' nsepar = ', nsepar
      write (6,*) '  nfile = ', nfile 
      write (6,*) ' iprint = ', iprint
C     WRITE (6,9400) nsepar,nfile,iprint

      WRITE (6,*) '  norbi = ', norbi 
      WRITE (6,*) '  norbo = ', norbo
      write (6,*) '  nlexi = ', nlexi
      write (6,*) '  nlexo = ', nlexo
C     WRITE (6,9400) norbi,norbo,nlexi,nlexo

C     READ (5,9080) (rss(j),j=1,nsepar)
      READ (5,*)    (rss(j),j=1,nsepar)
      write (6,*)
      WRITE (6,*) 'The nuclear separations are:'
      WRITE (6,9080) (rss(i),i=1,nsepar)
c
      nfil1 = nfile + 1
      do i = 1, nfil1
         READ (5,9111) iflnm(i)
      enddo

      write (6,*)
      write (6,*) ' Names of input/output files are :'
      WRITE (6,9111) (iflnm(i),i=1,nfil1)
      do i = 1, nfile
         iunit = 10 + i
         nmunit(i) = iunit
         if(nfile.eq.1) then
            open(iunit,file=iflnm(i),status='old', form='formatted')
         else
            open(iunit,file=iflnm(i),status='old', form='unformatted')
         endif
      enddo
      iunit1 = 10 + nfil1
      if(nfile.eq.1) then
         open(iunit1,file=iflnm(nfil1),status='unknown', form
     $        ='formatted')
      else
         open(iunit1,file=iflnm(nfil1),status='unknown', form
     $        ='unformatted')
      endif

      do ir = 1, noff
         if(nfile.ne.1) then
            do i=1,nfile
               idev=nmunit(i)
               do k=1,nlexi
                  READ (idev) (fin(j,k,i),j=1,norbi)
               enddo
            enddo
         else
C           READ (11,9300,end=950) (fnum(i),i=1,nsepar)
            do i = 1, nsepar
               READ (11,9301,end=950) RT, fnum(i)
            enddo
         endif
         do i = 1, nlexo
            do k = 1, norbo
               if(nfile.eq.1) go to 41
               do j = 1, nsepar
                  fnum(j) = fin(k,i,j)
               enddo
 41            ier = 0
C * imsl v9
C *     call icsicu(rss,fnum,nsepar,bpar,csplin,29,ier)
C *     call icsevu(rss,fnum,nsepar,csplin,29,rquad,vquad,nquad,ier)
C * imsl v10
C *      call dcsdec(nsepar,rss,fnum,0,0.d0,0,0.d0,bpar,csplin)
C *      do 51 ix=1,nquad
C *      vquad(ix)=dcsval(rquad(ix),nsepar-1,bpar,csplin)
C *   51 continue
C * numerical recipes version:
               call spline(rss,fnum,nsepar,1.d30,1.d30,csplin)
               do ix = 1, nquad
                  call splint(rss,fnum,csplin,nsepar,rquad(ix),vquad(ix)
     $                 )
               enddo
               
               if (ihoex .eq. 0) then
                  if(ivibopt .eq. 1) then
                     do iv = 1, ivibtr
                        sum = zero
                        do j = 1, nquad
                           sum = sum + chis(j,iv)*vquad(j)*wquad(j)
                        enddo
                        fout(k,1,iv,i) = sum
                     enddo
                     iout=1
                  else
                     do iv = 1, ivm
                        do ivp = iv, ivm
                           if (iv.gt.ivibtr.and.ivp.gt.ivibtr) goto 220
                           sum = zero
                           do j = 1, nquad
                              if (outvib) then
                                 zj = (rquad(j)/re - rx)*alpha
                                 temp = chis(j,iv)*chis(j,ivp)*vquad(j)
     $                                *wquad(j)
                                 sum = sum + temp*dexp(zj*zj)
                              else
                                 sum=sum+chis(j,iv)*chis(j,ivp)*vquad(j)
     $                                *wquad(j)
                              endif
                           enddo
                           sum = pref*sum
                           if(iprint.eq.1) fmat(iv,ivp) = sum
                           if(iv.le.ivibtr)  fout(k,ivp,iv,i) = sum
                           if(ivp.le.ivibtr) fout(k,iv,ivp,i) = sum
 220                       continue
                        enddo
                     enddo
                     iout=ivm
                  endif
               else
                  do iv=1,ihoex
                     do ivp=1,ivm
                        sum=zero
                        do j=1,nquad
                          sum=sum+ho(j,iv)*chis(j,ivp)*vquad(j)*wquad(j)
                        enddo
                        if(iprint.eq.1) fmat(iv,ivp)=sum
                        fout(k,iv,ivp,i)=sum
                     enddo
                  enddo
                  iout=ihoex
               endif

               if(iprint.ne.1) go to 240
               WRITE (6,9191) i
               do iv=1,ivm
                  ivr=iv-1+minvib-1
                  do ivp=iv,ivm
                     ivpr=ivp-1+minvib-1
                     WRITE (6,9201) ivr,ivpr,fmat(iv,ivp)
                  enddo
               enddo
 240           continue
C * end k loop
            enddo
C * end i loop
         enddo

         if(nfile.ne.1) then
            do l=1,nlexo
               do k=1,ivibtr
                  do j=1,iout
                     WRITE (iunit1) (fout(i,j,k,l),i=1,norbo)
                  enddo
               enddo
            enddo
         else
            do i=1,ivm
               WRITE (iunit1,9500) (fout(1,i,j,1),j=1,ivibtr)
            enddo
         endif
C * end ir loop
      enddo

      do i=1,nfile
         iunit=nmunit(i)
         close(iunit)
      enddo
      close(iunit1)
c
 950  return
 9111 format(10a8)
 9180 format(5i5)
 9080 format(7f10.5)
 9191 format(///' vibrational matrix elements for function',i3//
     15x,'v',3x,'vp',9x,'<f>')
 9201 format(1x,2i5,2x,1pd14.7)
 9300 format(e16.8)
 9301 format(f10.5,f10.5)
 9400 format(6i7)
 9500 format(5(1pe15.6))
      end
C *-/

C *- 
C * electron-nuclear potential
      double precision function potenu(r,rnl,rnr,zl,zr,lamb)
      implicit none
      
      double precision r,rnl,rnr,zl,zr,factor
      integer lamb
      
      factor = -1.0d0
      rnl = abs(rnl)
      rnr = abs(rnr)
      potenu = -zl*(min(r,rnl)**lamb)/(max(r,rnl)**(lamb+1))
      potenu = potenu*(factor**lamb)
      potenu = potenu - zr*(min(r,rnr)**lamb)/(max(r,rnr)**(lamb+1))

      end
