C * JAN 01 2015 - RDS - REMOVE SYMMETRY NUMBER DEPENDENCY
C * JUN 12 2006 - RDS - ADD GAUSS-LEGENDRE NODES IN EXCHANGE KERNEL INTEGRAL
C * JUN 06 2006 - RDS - UPDATE SPHRJ
C * APR 25 2006 - RDS - FIX BUG OF KERNCALC
C * DEC 28 2005 - RDS - REMOVE THE NAMES OF I/O FILES
C *                     REMOVE # OF GEOMETRIES SO THAT ONLY ONE GIVEN GEOMETRY ONE TIME
C *                     UPDATE CLEBEXCH/KERNCALC TO HANDLE E-N2
C * JAN 01 1990 - WT  - FIRST VERSION
C *-
C * PROGRAM EXKER
C *
C * program to calculate the exchange kernels at a given internuclear
C * geometries so they may be used in the construction of a
C * vibrational exchange kernel

C * Now, it is ONLY for homonuclear diatomic molecule! (Hao Feng)
C *                     ***********
      implicit none
      integer chanmax,ptsmax,lammax
      integer boundmax,nrgmax

      integer nlmomax
      parameter (boundmax=6, chanmax=21, ptsmax=1200, lammax=60,
     $           nlmomax=16, nrgmax=100)

      integer isymind
      integer nreg
      double precision rgs(ptsmax),wtt(ptsmax)
      double precision kernel(chanmax,chanmax,ptsmax)

C *   NBOUND: the number of bound target orbitals
C *     NPTS: the number of radial mesh points in the MO coefficients
C * NPTSEXCH: the number of radial mesh points in the exchange region
      integer nbound,npts,nptsexch

      integer l0, symlam

C * the exchange stuff

C * KERNLFIL: an array of the names of the files containing the
C *           kernel (the '.ker' will be appended internally here) for 
C *           all symmetries
      character*8 kernlfil
      character*8 sphnam

C * EXUNIT: the unit associated with the exchang kernel. I set it to 55
C * NEXDIM: the number of partial waves used in the CREATION of the
C *         Exchange kernel.  It is important to get this number right 
C *         so that the unformatted reads of the exchange kernel work
C *         out ok.
      integer nexdim

C * SPHRJ(): the spherical projections of the molecular orbitals
C *          of the target.
      double precision sphrj(boundmax,chanmax,ptsmax)

C * NLPROJ: the number of projections for the orbital
C * MLPROJ: the Ml for that orbital
C * L0PROJ: the L0 for the orbital
      integer nlproj(boundmax),mlproj(boundmax),l0proj(boundmax)

C     double precision clebx1(boundmax,chanmax,lammax,chanmax)
C     double precision clebx2(boundmax,chanmax,lammax,chanmax)
C     double precision clebx3(boundmax,chanmax,lammax,chanmax)
C     double precision clebx4(boundmax,chanmax,lammax,chanmax)
C     double precision clebx5(boundmax,chanmax,lammax,chanmax)
      double precision clebx1(boundmax,chanmax,lammax,nlmomax)
      double precision clebx2(boundmax,chanmax,lammax,nlmomax)
      double precision clebx3(boundmax,chanmax,lammax,nlmomax)
      double precision clebx4(boundmax,chanmax,lammax,nlmomax)
      double precision clebx5(boundmax,chanmax,lammax,nlmomax)

      double precision faclog(500)
      double precision rac,u9
      integer ia,ib,ic,id,ie,if,l9(9),kselc
      common/faclg/faclog,rac,u9,ia,ib,ic,id,ie,if,l9,kselc
      character*4 dat,ker
      integer i,iunit,ounit
      integer iwhich

C * set up the factorial table
      dat='.dat'
      ker='.ker'
      kselc=1

      call facset

      iunit=10
      ounit=11

      call indata(l0,symlam,nbound,sphnam,kernlfil,nexdim,nlproj,mlproj
     $     ,l0proj,npts,nptsexch,rgs,wtt)

C * judge if the number of partial waves is enough
      call judgenpwav(nbound,nlproj,nexdim)

      call sphjread(sphrj,nlproj,nbound,npts,nptsexch,rgs,nexdim,iunit)

      if(l0.eq.0 .and. symlam.eq.0) then
         isymind = 1
      else if(l0.eq.1 .and. symlam.eq.0) then
         isymind = 2
      else if(l0.eq.1 .and. symlam.eq.1) then
         isymind = 3
      else if(l0.eq.2 .and. symlam.eq.1) then
         isymind = 4
      else
         WRITE (6,*) 'l0=',l0, 'symlam=', symlam
         stop 'error---symmetry outside range'
      endif

C * fkerncalc is just a choice to test clebexch/kernel. It is NOT optimized! 
C * Apr. 25, 2006 --- Hao Feng
      iwhich = 0
      if (iwhich .eq. 0) then
         call clebexch(clebx1,clebx2,clebx3,clebx4,clebx5
     $        ,symlam,l0,nlproj,mlproj,l0proj,nbound
     $        ,nexdim)


C * Calculate the exchange Kernel for this symmetry 
C * of course, you should only do this once for all energies

         call kerncalc(sphrj,kernel,clebx1,clebx2,clebx3,clebx4
     $        ,clebx5,nlproj,l0proj,ounit,l0,nexdim
     $        ,nbound,rgs,nptsexch)
      else if (iwhich .eq. 1) then
         call fkerncalc(sphrj,kernel,symlam,nlproj,mlproj,l0proj
     $        ,ounit,l0,nexdim,nbound,rgs,nptsexch)
      endif
      
      close(unit=11)

      close(unit=10)
      stop
      end

C *- 
C * This subroutine reads in the input

      subroutine indata(l0,symlam,nbound,sphnam,kernfil,nexdim
     $     ,nlproj,mlproj,l0proj,npts,nptsexch,rgs,wtt)

      implicit none
      integer boundmax, nrgmax, ptsmax
      parameter(boundmax=6, nrgmax=100, ptsmax=1200)
      integer l0,symlam
      integer nbound,nexdim,nreg
      character*8 sphnam,kernfil
      integer nlproj(boundmax),mlproj(boundmax),l0proj(boundmax)
      integer i, l
      character*4 dat

      integer ity,nps,jz,ic,npts,nptsexch
      double precision rs,rfn,rfs,del(nrgmax),dell,r1,wx,wttp,rttp
      double precision wtt(ptsmax),rgs(ptsmax)
      double precision wnode(ptsmax),rnode(ptsmax)

      dat='.dat'

      READ (5,*) nbound
      WRITE (6,*)' There are',nbound,' target molecular orbitals'
      WRITE (6,*)

      READ (5,*) npts
      WRITE(6,*) ' There are',npts,' points in every molecular orbitals'
      WRITE (6,*) 

C * nreg = no. of regions
C *
C * ity = 0 - Gauss-Legendre pts and wts
C *       1 - Trapezodial pts and wts
C * nps = No. of pts. in a region
C *  rs = starting radius of region
C * rfn = ending radius of region

      READ (5,*) nreg
      WRITE (6,*) ' There are',nreg,' mesh regions'
      WRITE (6,*)
      WRITE (6,*)'ity nps    rs     rfn '
      ic = 0
      do i = 1, nreg
         READ (5,*) ity, nps, rs, rfn
         WRITE (6,'(i2,i5,f8.3,f8.3)')ity, nps, rs, rfn
         rfs = rfn - rs
         if (ity .eq. 0) then
            del(i) = rfs/nps
            r1 = rs - del(i)
            if (i .eq. 1) then
               r1 = rs
               nps = nps - 1
            endif
            do jz = 1, nps
               dell = del(i)
               ic = ic + 1
               r1 = r1 + dell
               wx = 1.0d0
               if (ic .eq. 1) wx = 0.5d0
               if (jz .eq. 1 .and. i .ne. 1) then
                  wx = 0.5d0
                  dell = del(i) + del(i-1)
               endif
               wtt(ic) = wx*dell
               rgs(ic) = r1
            enddo
         else if (ity .eq. 1) then
            do jz = 1, nps
               ic = ic + 1
               if (nps .gt. 3) then
                  call lgndrx(nps,jz,wttp,rttp)
               else if (nps .le. 3) then
                  call lgndx2(nps,jz,wttp,rttp)
               endif
               wtt(ic) = wttp*rfs
               rgs(ic) = rttp*rfs + rs
            enddo
         else if (ity .eq. 2) then
            call gauleg(rs,rfn,rnode,wnode,nps)
            do jz = 1, nps
               ic = ic + 1
               rgs(ic) = rnode(jz)
               wtt(ic) = wnode(jz)
            enddo
         endif
      end do
      nptsexch = ic + 1
      rgs(nptsexch) = rfn
      wtt(nptsexch) = 0.0d0
      write(6,*) 
      write(6,*) 'Now the nodes and weights are: (',nptsexch,' points)'
      do jz = 1, nptsexch
         write(6,'(i5,f9.5,f14.10)') jz,rgs(jz),wtt(jz)
      enddo

      WRITE (6,*)
      WRITE (6,*)'   l0   Symlam'
      READ (5,*) l0, symlam
      WRITE (6,*)l0, symlam

      WRITE (6,*)
      WRITE (6,*)'nlproj  mlproj  l0proj'

      do i = 1, nbound
         READ (5,*) nlproj(i),mlproj(i),l0proj(i)
         WRITE (6,*)nlproj(i),mlproj(i),l0proj(i)
      end do

      WRITE (6,*)
      WRITE (6,*) 'ExChans'

      READ (5,*) nexdim
      WRITE (6,*) nexdim
      write (6,*)

      return
      end
C *-/

C *-
C * This routine reads in the spherical projections of the Target MO's
C * as from the output of the ALAM code.

      Subroutine sphjread(sphrj,nlproj,nbound,npts,nptsexch,rgs,nexdim
     $     ,iunit)
      implicit none
      integer nexdim,npts,nptsexch,nbound,iunit
      integer i,j,k,nlproj(nbound)
      double precision sphrj(nbound,nexdim,nptsexch),r(npts),sj(npts)
     $     ,sy(nptsexch)
      double precision rgs(nptsexch),csplin(npts)

C * read in the r mesh points from the sphj file
C * go to npts+1 because sphj file has point at r=0.0

      WRITE (6,*) 'in sphrj'
      WRITE (6,*) 'nbound    nlproj(1)         npts'
      WRITE (6,*) nbound, nlproj(1), npts

      read (iunit) r(1)
      write (60,*) r(1)
      do i = 1, npts 
         READ (iunit) r(i)
         WRITE (60,*) r(i)
      end do

C * read in the spherical projs from the ALAM output file

      do i = 1, nbound
         do j = 1, nlproj(i)
C * as above, we read in the points at r=0.0 but we don't use it since
C * we don't have it on the potential mesh
            READ (iunit) sj(1)
            do k = 1, npts
               READ (iunit) sj(k)
               write(61,*) r(k),sj(k)
            end do
C * interpolate sphrj(nptsexch) through original MO coefficients sj(npts)
            call spline(r,sj,npts,1.d30,1.d30,csplin)
            do k = 1, nptsexch
               call splint(r,sj,csplin,npts,rgs(k),sy(k))
               sphrj(i,j,k) = sy(k)
               write(62,*) rgs(k),sphrj(i,j,k)
            enddo
         end do
      end do

      return
      end
C *-/


C *-
C * calculate the Clebsch Gordon elements to be used
C * when calculating the exchange potential

      Subroutine clebexch(clebx1,clebx2,clebx3,clebx4,clebx5,symlam,l0
     $     ,nlproj,mlproj,l0proj,nbound,nexdim)
      implicit none
      integer nexdim,nbound
      integer symlam,l0
      integer nlproj(nbound),l0proj(nbound),mlproj(nbound)

      integer boundmax,chanmax,lammax,nlmomax
      parameter (boundmax=6,chanmax=21,lammax=60,nlmomax=16)

      double precision clebx1(boundmax,chanmax,lammax,nlmomax)
      double precision clebx2(boundmax,chanmax,lammax,nlmomax)
      double precision clebx3(boundmax,chanmax,lammax,nlmomax)
      double precision clebx4(boundmax,chanmax,lammax,nlmomax)
      double precision clebx5(boundmax,chanmax,lammax,nlmomax)

      integer m,mi,ibnd,i,j,k,l,lp,lam,lamlb,lamub

C * connect to the Cleb Gordon routines:

      double precision faclog(500)
      double precision rac,u9
      integer ia,ib,ic,id,ie,if,l9(9),kselc
      common/faclg/faclog,rac,u9,ia,ib,ic,id,ie,if,l9,kselc

C *- Now we only deal with sigma and pi orbitals (Hao Feng, Dec. 28, 2005)

      m = symlam
      do ibnd = 1, nbound
         mi = mlproj(ibnd)
         call iferr(mi)
         do i = 1, nexdim
            l = 2*(i-1) + l0
            do k = 1, nlproj(ibnd)
               lp = 2*(k-1) + l0proj(ibnd)
               lamlb = iabs(l - lp)
               lamub = l + lp
               if (lamub .gt. lammax) then
                  write(0,*) "lammax = ",lammax ," and should be >/="
     $                 ,lamub
                  stop ' INCREASE lammax in source!'
               endif
               do lam = lamlb, lamub
                  clebx1(ibnd,i,lam,k) = 0.0d0
                  clebx2(ibnd,i,lam,k) = 0.0d0
                  clebx3(ibnd,i,lam,k) = 0.0d0
                  clebx4(ibnd,i,lam,k) = 0.0d0
                  clebx5(ibnd,i,lam,k) = 0.0d0

                  if (mod(l+lp+lam,2) .eq. 0) then
C *-
C * Hao Feng added clebx4 and clebx5 to calculate pi-orbital. Dec. 28, 2005
C *-/
                     call cleb(l,lam,lp,-m,m-mi,-mi)
                     clebx1(ibnd,i,lam,k) = rac
                     if (mi .eq. 1) then
                        call cleb(l,lam,lp,-m,m+mi,mi)
                        clebx4(ibnd,i,lam,k) = rac
                     else if (mi .eq. 0) then
                        clebx4(ibnd,i,lam,k) = 0.0
                     endif
                     call cleb(l,lam,lp,0,0,0)
                     clebx2(ibnd,i,lam,k) = rac
                     call cleb(l,lam,lp,m,mi-m,mi)
                     clebx3(ibnd,i,lam,k) = rac
                     if (mi .eq. 1) then
                        call cleb(l,lam,lp,m,-mi-m,-mi)
                        clebx5(ibnd,i,lam,k) = rac
                     else if (mi .eq. 0) then
                        clebx5(ibnd,i,lam,k) = 0.0
                     endif
                  endif
               end do
            end do
         end do
      end do
      return
      end
C *-/

C *-
C * this routine calculates the exchange kernel an writes it to
C * a file.  It will be used later to calculate the
C * exchange potential

      Subroutine Kerncalc (sphrj,kernel,clebx1,clebx2,clebx3,clebx4
     $     ,clebx5,nlproj,l0proj,exunit,l0,nexdim,nbound,rgs,nptsexch)

      implicit none
      integer exunit
      integer nexdim,nptsexch,nbound
      integer nlproj(nbound),l0proj(nbound)
      double precision rgs(nptsexch)
      double precision kernel(nexdim,nexdim,nptsexch)
      double precision sphrj(nbound,nexdim,nptsexch)
      integer l0
      integer i,j,k,il,ilp
      integer l,lp,l2p,l3p
      integer ipts,jpts,ibnd
      double precision ralfa,rbeta,rless,rgrtr,rfact

      integer boundmax,chanmax,lammax,nlmomax
      parameter (boundmax=6,chanmax=21,lammax=60,nlmomax=16)

      double precision clebx1(boundmax,chanmax,lammax,nlmomax)
      double precision clebx2(boundmax,chanmax,lammax,nlmomax)
      double precision clebx3(boundmax,chanmax,lammax,nlmomax)
      double precision clebx4(boundmax,chanmax,lammax,nlmomax)
      double precision clebx5(boundmax,chanmax,lammax,nlmomax)

      double precision numer
      double precision lamterm,l2term,l3term,ibndterm
      double precision clebtemp
      integer lam,lamlb,lamub

      do ipts = 1, nptsexch
         rbeta = rgs(ipts)
         do jpts = 1, nptsexch
            ralfa = rgs(jpts)
            if(ralfa.gt.rbeta) then
               rgrtr = ralfa
               rless = rbeta
            else
               rgrtr = rbeta
               rless = ralfa
            endif
            do il = 1, nexdim
               l = 2*(il-1) + l0
               do ilp = 1, nexdim
                  lp = 2*(ilp-1) + l0
                  numer = dsqrt((2.d0*l+1.d0)*(2.d0*lp+1.d0))
                  
                  ibndterm = 0.d0
                  do ibnd = 1, nbound

                     l2term = 0.d0
                     do i = 1, nlproj(ibnd)
                        l2p = 2*(i-1) + l0proj(ibnd)

                        l3term = 0.d0
                        do j = 1, nlproj(ibnd)
                           l3p = 2*(j-1) + l0proj(ibnd)

                           lamterm = 0.d0
                           lamlb = max(iabs(l-l2p), iabs(lp-l3p))
                           lamub = min(l+l2p, lp+l3p)
                           do lam = lamlb, lamub
                              if (mod(l+l2p+lam,2) .eq. 0 .and.
     $                             mod(lp+l3p+lam,2) .eq. 0) then 
                                 rfact = rless**(lam)/rgrtr**(lam+1)
                                 clebtemp = (clebx1(ibnd,il,lam,i)
     $                                *clebx3(ibnd,ilp,lam,j)
     $                                +clebx4(ibnd,il,lam,i)
     $                                *clebx5(ibnd,ilp,lam,j))
     $                                *clebx2(ibnd,il,lam,i)
     $                                *clebx2(ibnd,ilp,lam,j)

                                 lamterm = lamterm + clebtemp *rfact
                              endif
                           end do
                           l3term = l3term + lamterm*sphrj(ibnd,j,jpts)
     $                          /dsqrt(2.d0*l3p+1.d0)
                        end do
                        l2term = l2term + l3term*sphrj(ibnd,i,ipts)
     $                       /dsqrt(2.d0*l2p+1.d0)
                     end do
                     ibndterm = ibndterm + l2term
                  end do
                  kernel(il,ilp,jpts) = ibndterm*numer
               end do
            end do
         end do
         WRITE (exunit) kernel
      end do

      return
      end
C *-/

C  may not want to use the array below as the technique
C  for storing the arrays
C
C  Don't forget to call FACSET as the first line f the scattering
C  code so that the factorials are set up before we get to the
C  Clebsch Gordon Coeffs
C
C       double precision faclog(500)
C       common /clebG/ icleb
C       common/clebgor/clebgor
C       COMMON/FACLG/ FACLOG,RAC,U9,IA,IB,IC,ID,IE,IF,L9(9),KSELC
C
C********************************************************************

C *-
      SUBROUTINE CLEB(j1,j2,j,m1,m2,m)
C-----------------------------------------------------------------
C      REFERENCE   T. TAMURA   COMPUTER PHYS. COMM. VOL. 1, P. 337
C      ERRATUM  COMPUTER PHYS. COMM. VOL. 2, P. 174
C
C      CALLING SEQUENCE...
C      CALL FACSET ONCE BEFORE ANY CALLS TO CLEB, RACAH, OR NINEJ
C      SET IA=2*J1, IB=2*J2, IC=2*J
C      SET ID=2*M1, IE=2*M2, IF=2*M
C      CALL CLEB
C      THE CLEBSCH-GORDAN COEFFICIENT (J1 M1 J2 M2/ J M) WILL BE
C      RETURNED IN THE VARIABLE RAC
C      CALLING PARAMETERS ARE PRESERVED
C-----------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer j1,j2,j,m1,m2,m
      common /clebG/ icleb
      COMMON/FACLG/ FACLOG(500),RAC,U9,IA,IB,IC,ID,IE,IF,L9(9),KSELC
      icleb=icleb+1
      ia=2*j1
      ib=2*j2
      ic=2*j
      id=2*m1
      ie=2*m2
      if=2*m
      RAC=0.0D+00
      IG4=IA+IB+IC
      IF(KSELC.EQ.0) GO TO 10
      IF(ID+IE.NE.IF) GO TO 70
      IF(MOD(IG4,2).NE.0) GO TO 70
      IF(IA+IB-IC.LT.0.OR.IC-IABS(IA-IB).LT.0) GO TO 70
      IF(MIN0(IA-IABS(ID),IB-IABS(IE),IC-IABS(IF)).LT.0) GO TO 70
      IF(MOD(IB+IE,2).NE.0.OR.MOD(IC+IF,2).NE.0) GO TO 70
 10   IF(IA.EQ.0.OR.IB.EQ.0) GO TO 20
      IF(IC.GT.0) GO TO 40
      IF(IC.EQ.0) GO TO 30
      IF(IC.LT.0) GO TO 70
 20   RAC=1.0D+00
      GO TO 70
 30   FB=IB+1
      RAC=((-1.0D+00)**((IA-ID)/2))/DSQRT(FB)
      GO TO 70
 40   IF(ID.NE.0.OR.IE.NE.0) GO TO 50
      IG2=IG4/2
      IF(MOD(IG2,2).NE.0) GO TO 70
      IG=IG2/2
      I1=(IA+IB-IC)/2+1
      I2=(IC+IA-IB)/2+1
      I3=(IC+IB-IA)/2+1
      I4=IG2+2
      I5=IG+1
      I6=(IG2-IA)/2+1
      I7=(IG2-IB)/2+1
      I8=(IG2-IC)/2+1
      F1=DEXP(0.5D+00*(FACLOG(I1)+FACLOG(I2)+FACLOG(I3)-FACLOG(I4))
     $     +(FACLOG(I5)-FACLOG(I6)-FACLOG(I7)-FACLOG(I8)))
      F2=IC+1
      F2=DSQRT(F2)
      S1=1-2*MOD((IG2+IC)/2,2)
      RAC=S1*F1*F2
      GO TO 70
 50   FC2=IC+1
      IABCP=IG4/2+1
      IABC=IABCP-IC
      ICAB=IABCP-IB
      IBCA=IABCP-IA
      IAPD=(IA+ID)/2+1
      IAMD=IAPD-ID
      IBPE=(IB+IE)/2+1
      IBME=IBPE-IE
      ICPF=(IC+IF)/2+1
      ICMF=ICPF-IF
      SQFCLG=0.5D+00*(DLOG(FC2)-FACLOG(IABCP+1)+FACLOG(IABC)
     $     +FACLOG(ICAB)+FACLOG(IBCA)+FACLOG(IAPD)+FACLOG(IAMD)
     $     +FACLOG(IBPE)+FACLOG(IBME)+FACLOG(ICPF)+FACLOG(ICMF))
      NZMIC2=(IB-IC-ID)/2
      NZMIC3=(IA-IC+IE)/2
      NZMI=MAX0(0,NZMIC2,NZMIC3)+1
      NZMX=MIN0(IABC,IAMD,IBPE)
      S1=(-1.0D+00)**(NZMI-1)
      DO NZ=NZMI,NZMX
         NZM1=NZ-1
         NZT1=IABC-NZM1
         NZT2=IAMD-NZM1
         NZT3=IBPE-NZM1
         NZT4=NZ-NZMIC2
         NZT5=NZ-NZMIC3
         TERMLG=SQFCLG-FACLOG(NZ)-FACLOG(NZT1)-FACLOG(NZT2)-FACLOG(NZT3)
     $        -FACLOG(NZT4)-FACLOG(NZT5)
         SSTERM=S1*DEXP(TERMLG)
         RAC=RAC+SSTERM
         S1=-S1
      ENDDO
      IF(DABS(RAC).LT.1.0D-12) RAC=0.0D+00
 70   RETURN
      END
C *-/

C *-
      SUBROUTINE FACSET
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/FACLG/ FACLOG(500),RAC,U9,IA,IB,IC,ID,IE,IF,L9(9),KSELC

      FACLOG(1)=0.0D+00
      FACLOG(2)=0.0D+00
      F1=1.0D+00
      DO 10 N=3,500
         F1=F1+1.0D+00
         FACLOG(N)=FACLOG(N-1)+DLOG(F1)
 10   CONTINUE
      RETURN
      END
C *-/

C *-
C * another choice to cal. exchange kernel
C * this routine calculates the exchange kernel and writes it to
C * a file.  It will be used later to calculate the
C * exchange potential

      Subroutine fkerncalc(sphrj,kernel,symlam,nlproj,mlproj,l0proj
     $     ,exunit,l0,nexdim,nbound,rgs,nptsexch)

      implicit none
      integer exunit
      integer nexdim,nptsexch,nbound
      integer nlproj(nbound),mlproj(nbound),l0proj(nbound)
      integer symlam
      double precision rgs(nptsexch)
      double precision kernel(nexdim,nexdim,nptsexch)
      double precision sphrj(nbound,nexdim,nptsexch)
      integer l0
      integer i,j,k,il,ilp
      integer l,lp,l2p,l3p
      integer ipts,jpts,ibnd
      double precision ralfa,rbeta,rless,rgrtr,rfact
      double precision numer
      double precision lamterm,l2term,l3term,ibndterm

      integer nlam,lmda,lam0, lamlb, lamub
      integer mi,itrig1,itrig2,lmdale,lmdagr
      double precision kerntemp
      double precision excoef

      do ipts = 1, nptsexch
         rbeta = rgs(ipts)
         do jpts = 1, nptsexch 
            ralfa = rgs(jpts)
            if(ralfa.gt.rbeta) then
               rgrtr = ralfa
               rless = rbeta
            else
               rgrtr = rbeta
               rless = ralfa
            endif
            do il = 1, nexdim
               l = 2*(il-1) + l0
               do ilp = 1, nexdim
                  lp = 2*(ilp-1) + l0
                  
                  ibndterm = 0.0d0
                  numer = dsqrt((2.d0*l+1.d0)*(2.d0*lp+1.d0))
                  do ibnd = 1, nbound
                     mi = mlproj(ibnd)
                     call iferr(mi)
                     
                     l2term = 0.0d0
                     do i = 1, nlproj(ibnd)
                        l2p = 2*(i-1) + l0proj(ibnd)
                        
                        l3term = 0.0d0
                        do j = 1, nlproj(ibnd)
                           l3p = 2*(j-1) + l0proj(ibnd)
                           
                           lamterm = 0.0d0
                           lamlb = iabs(l-l2p)
                           lam0 = iabs(lp-l3p)
                           if (lamlb .lt. lam0) lamlb = lam0
                           lamub = l + l2p
                           lam0 = lp + l3p
                           if (lamub .gt. lam0) lamub = lam0
                           do lmda = lamlb, lamub
                              if (mod(l+l2p+lmda,2) .eq. 0 .and. mod(lp
     $                             +l3p+lmda,2) .eq. 0) then
                                 rfact = rless**(lmda)/rgrtr **(lmda+1)
                                 lamterm = lamterm + rfact *excoef(l,lp
     $                                ,l2p,l3p,lmda ,symlam,mi)
                              endif
                           enddo
                           l3term = l3term + lamterm*sphrj(ibnd,j,jpts)
     $                          /dsqrt(2.0d0*l3p+1.0d0)
                        enddo
                        l2term = l2term + l3term*sphrj(ibnd,i,ipts)
     $                       /dsqrt(2.0d0*l2p+1.0d0)
                     enddo
                     ibndterm = ibndterm + l2term
                  enddo
                  kernel(il,ilp,jpts) = ibndterm*numer
               enddo
            enddo
         enddo
         write(exunit) kernel
      enddo
      
      return
      end
C *-/

C *-
C * coefficients of exchange kernel 
C * A(l,lp,l2p,l3p;0,LAM) ---- sigma orbital
C * B(l,lp,l2p,l3p;0,LAM) ---- pi orbital

      double precision function excoef(l,lp,l2p,l3p,lmda,symlam,mi)
      implicit none
      integer l,lp,l2p,l3p,lmda,symlam,mi
      double precision glam

      if (mi .eq. 0) then
         excoef = glam(l,lp,l2p,l3p,lmda,symlam,0)
      else if (mi .eq. 1) then
         excoef = glam(l,lp,l2p,l3p,lmda,symlam,1) 
     $          + glam(l,lp,l2p,l3p,lmda,symlam,-1)
      else
         write (0,*) " Now mi = 0, 1 could be workd, and your mi = ",mi
         stop
      endif

      return 
      end
C *-/

C *-
      double precision function glam(l,lp,l2p,l3p,lmda,symlam,mi)
      implicit none
      integer l,lp,l2p,l3p,lmda,symlam,mi
      double precision numer,CG

      double precision faclog(500)
      double precision rac,u9
      integer ia,ib,ic,id,ie,if,l9(9),kselc
      common/faclg/faclog,rac,u9,ia,ib,ic,id,ie,if,l9,kselc

C     numer = (2.0d0*l+1.0d0)*(2.0d0*lp+1.0d0)/((2.0d0*l2p+1.0d0)*(2.0d0
C    $     *l3p+1.0d0))
C     numer = dsqrt(numer)

      glam = CG(l,lmda,l2p,-symlam,symlam-mi,-mi)
     $      *CG(l,lmda,l2p,0,0,0)
     $      *CG(lp,lmda,l3p,symlam,mi-symlam,mi)
     $      *CG(lp,lmda,l3p,0,0,0)

C     call cleb(l,lmda,l2p,-symlam,symlam-mi,-mi)
C     glam = numer*rac
C     call cleb(l,lmda,l2p,0,0,0)
C     glam = glam*rac
C     call cleb(lp,lmda,l3p,symlam,mi-symlam,mi)
C     glam = glam*rac
C     call cleb(lp,lmda,l3p,0,0,0)
C     glam = glam*rac

      return
      end
C *-/

C *-
C * Clebsh-Gordon Coefficients; Use 3-j subroutine
C * 
      DOUBLE PRECISION FUNCTION CG (L1,L2,L3,M1,M2,M3)
      IMPLICIT REAL*8 (A-H,O-Z)

      ll1 = 2*l1
      ll2 = 2*l2
      ll3 = 2*l3
      mm1 = 2*m1
      mm2 = 2*m2
      mm3 = 2*m3
      v1 = f3j(ll1,ll2,ll3,mm1,mm2,-mm3)
      cmult = 1.0
      nexp = l2 - l1 - m3
      nexp = iabs(nexp)
      if (mod(nexp,2).eq.0) go to 10
      cmult = -1.0
 10   factor = sqrt(2.*float(l3)+1.)*cmult
      cg = factor*v1
      return
      end
C *-/

C *-
C * 3-j subroutnie

      REAL*8 FUNCTION F3J (JD1,JD2,JD3,MD1,MD2,MD3)
      IMPLICIT REAL*8 (A-H,O-Z)
C *
C * NBS ROUTINE FOR CALCULATING THE WIGNER 3-J SYMBOLS
C * F3J VERSION II FORTRAN IV
C *
      dimension mtri(9)
      common /fact/ fl(322), ncall

      dfloat(i) = i
      j1 = jd1
      j2 = jd2
      j3 = jd3
      m1 = md1
      m2 = md2
      m3 = md3

C * Determine whether to calculate fl(n) s
      if (ncall+1867) 10,30,10
 10   ncall = -1867
      fl(1) = 0.d0
      fl(2) = 0.d0
      do 20 n = 3, 322
         fn = n-1
 20      fl(n) = fl(n-1) + dlog(fn)

 30   i = j1 + j2 - j3
      i1 = i/2
      if (i-2*i1) 340,40,340
 40   mtri(1) = i1
      i = j1 - j2 + j3
      i1 = i/2
      if (i-2*i1) 340,50,340
 50   mtri(2) = i1
      i = -j1 + j2 + j3
      i1 = i/2
      if (i-2*i1) 340,60,340
 60   mtri(3) = i1
      if (m1+m2+m3) 340,70,340
 70   i = j1+m1
      i1 = i/2
      if (i-2*i1) 340,80,340
 80   mtri(4) = i1
      mtri(5) = (j1-m1)/2
      i = j2 + m2
      i1 = i/2
      if (i-2*i1) 340,90,340
 90   mtri(6) = i1
      mtri(7) = (j2-m2)/2
      i = j3 + m3
      i1 = i/2
      if (i-2*i1) 340,100,340
 100  mtri(8) = i1
      mtri(9) = (j3-m3)/2
      do 110 n = 1,9
         if (mtri(n)) 340,110,110
 110  continue
      if (j3-j2+m1) 120,130,130
 120  kmin = -j3+j2-m1
      go to 140
 130  kmin = 0
 140  if (-j3+j1+m2-kmin) 160,160,150
 150  kmin = -j3 + j1 + m2
 160  kmin = kmin/2
      if (j2-j3+m1) 170,180,180
 170  kmax = j1 + j2 - j3
      go to 190
 180  kmax = j1 - m1
 190  if (j2+m2-kmax) 200,210,210
 200  kmax = j2 + m2
 210  kmax = kmax/2
      mini = mtri(1) - kmin + 1
      min2 = mtri(5) - kmin + 1
      min3 = mtri(6) - kmin + 1
      min4 = (j3-j2+m1)/2 + kmin
      min5 = (j3-j1-m2)/2 + kmin

C * sum series in double precision

      uk = 1.d-10
      s = 1.d-10
      ncut = 0
      kmax = kmax - kmin
      if (kmax) 260,260,220
 220  do 250 k = 1, kmax
         uk = uk*dfloat(mini-k)*dfloat(min2-k)*dfloat(min3-k)
     $        /(dfloat(kmin+k)*dfloat(min4+k)*dfloat(min5+k))
         uk = -uk
         if (dabs(uk)-1.d30) 240,230,230
 230     uk = 1.d-10*uk
         s = 1.d-10*s
         ncut = ncut+1
 240     if (dabs(uk)-1.d-20) 260,250,250
 250     s = s + uk

C * Calculate Delta Functions

 260  delog = 0.d0
      do 270 n = 1,9
         num = mtri(n)
 270     delog = delog+fl(num+1)
      num = (j1+j2+j3)/2+2
      delog = 0.5d0*(delog-fl(num))
      ulog = -fl(kmin+1) - fl(mini) - fl(min2) - fl(min3) - fl(min4+1) -
     $     fl(min5+1)
      plog = delog + ulog
      if (sngl(plog)+80.0) 290,280,280
 280  if (ncut) 300,300,290
 290  sig = dsign(1.d0,s)
      s = dabs(s)
      slog = dlog(s) + dfloat(ncut+1)*dlog(1.d+10)
      f3j = sig*dexp(slog+plog)
      go to 310
 300  s = s*1.d+10
      p = dexp(plog)
      f3j = p*s
 310  num = kmin + (j1-j2-m3)/2
      if (mod(num,2)) 320,330,320
 320  f3j = -f3j
 330  continue
      go to 350
 340  f3j = 0.0
 350  return
      end
C *-/

C *-
      subroutine iferr(mi)
      integer mi

      if (mi .ne. 0 .and. mi .ne. 1) then
         write(0,*) "Only sigma and pi orbitals (mi = 0, 1) "
         write(0,*) "could be worked. Now your orbital is: mi = ", mi
         stop
      endif

      return
      end
C *-/

C *-
C * judge if the number of partial waves is enough

      subroutine judgenpwav(nbound,nlproj,nexdim)
      implicit none
      
      integer boundmax
      parameter(boundmax=6)
      integer nexdim,nlproj(boundmax)

      integer nbound,i,nlmax

      nlmax = 0
      do i = 1, nbound
         if(nlproj(i) .gt. nlmax) nlmax = nlproj(i)
      enddo

      if(nexdim .lt. nlmax) then
         write (0,*)
         write(0,*) 'nexdim is TOO SMALL and should be greater than',
     $        nlmax
         write(0,*) 'or nlproj is TOO LARGE and should be less than'
     $        ,nexdim
         stop
      endif

      end
C *-/
