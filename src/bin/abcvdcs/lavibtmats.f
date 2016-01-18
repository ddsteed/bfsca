C * JUN 19 2006 - RDS - FIX A BUG OF READING K MATRIX IF THE K-MAT DOES NOT INCLUDE THE 
C *                   - INPUT ENERGY
C *

      subroutine lavibtmats(Einp,d,tre,tim,ldim,lmax,lsymmax,worst,k0,k
     $     ,alpha0,alpha2,q,lbext,noabc,tform,debug,nchs)
c
c  this module reads in the data and K matrices from the
c  lavib code and returns the necessary stuff for abcvdcs
c
      implicit none
      integer ldim,lmax,lbext,noabc,tform,debug,nchs
      double precision tre(0:ldim,0:ldim,0:ldim)
      double precision tim(0:ldim,0:ldim,0:ldim)
      integer lsymmax(0:ldim)
      double precision worst(2,0:ldim)
      double precision k0,k,alpha0,alpha2,q
      double precision d(2)
      double precision Einp

      double precision scatten,we,wexe,sce,cmtoau,ev,ev0
      integer iswit
      integer v0,v,numv,numsym,maxsym,maxv,i,j,jj,maxchan
      integer lam,l0,lm,nch,nv,lend
      parameter (maxsym=100)
      parameter (maxv=25)
      parameter (maxchan=500)
      integer mflen,kflen(maxsym)
      character*512 momfile,kmatfile(maxsym)
      integer iread

      double precision a0(0:maxv-1,0:maxv-1)
      double precision a2(0:maxv-1,0:maxv-1)
      double precision qmat(0:maxv-1,0:maxv-1)
      double precision kmat(maxchan,maxchan)
      double precision tmtre(maxchan,maxchan)
      double precision tmtim(maxchan,maxchan)
      double precision temp(maxchan,maxchan)

      read(5,*) scatten,v0,v,numv,numsym
      Einp = scatten

      if (numsym.gt.maxsym) then
         write(6,*) 'numsym of',numsym,' > maxsym of',maxsym,
     $        ' in lavibtmats - check or edit and recompile'
         stop 'parameter check error'
      endif

      cmtoau = 1.239841d-4/(2.d0*13.6058d0)
      read(5,*) iswit,we,wexe
      if (iswit .eq. 0) then
         we = we/cmtoau
         wexe = wexe/cmtoau
      endif
c
c  we will require that our filenames and path have no spaces
c
      read(5,*) momfile
      mflen = index(momfile,' ')
      do i=1,numsym
        read(5,*) kmatfile(i)
        kflen(i) = index(kmatfile(i),' ')
      enddo
      write(6,*) '  Scattering energy (eV) is ',scatten
      write(6,*) '  Choosing v0,v = ',v0,v
      write(6,*) '  Assuming ',numv,' vibrational levels'
      write(6,*) '     (this affects input of moments)'
      write(6,*) '  Will read in ',numsym,' symmetries'
      write(6,*) '  we,wexe = ',we,wexe
      write(6,*) '  Moment matrix file is: ',momfile(1:mflen)
      write(6,*) '  K-matrix files:'
      do i=1,numsym
        write(6,*) '     ',kmatfile(i)(1:kflen(i))
      enddo
c
c  convert input scattering energy to Rydbergs for comparison
c  with the K matrix files
c
      scatten = scatten*2.d0/27.2116d0
      write(6,*) '  Scattering energy input (Ry) =',scatten
c
c  read the matrix elements of the moments as used by lavib
c
      open(7,file=momfile(1:mflen),form='formatted',status='old')
      if (tform .eq. 1) then
         do i=0,numv-1
            read(7,'(5e15.5)') (a0(i,j),j=0,numv-1)
         enddo
         do i=0,numv-1
            read(7,'(5e15.5)') (a2(i,j),j=0,numv-1)
         enddo
         do i=0,numv-1
            read(7,'(5e15.5)') (qmat(i,j),j=0,numv-1)
         enddo
      else
         do i=0,numv-1
C           read(7,'(5e14.6)') (a0(i,j),j=0,numv-1)
            read(7,*) (a0(i,j),j=0,numv-1)
         enddo
         do i=0,numv-1
C           read(7,'(5e14.6)') (a2(i,j),j=0,numv-1)
            read(7,*) (a2(i,j),j=0,numv-1)
         enddo
         do i=0,numv-1
C           read(7,'(5e14.6)') (qmat(i,j),j=0,numv-1)
            read(7,*) (qmat(i,j),j=0,numv-1)
         enddo
      endif
      close(7)
      alpha0 = a0(v0,v)
      alpha2 = a2(v0,v)
      q = qmat(v0,v)
      write(6,*) '  Moment matrix elements chosen are:'
      write(6,*) '  alpha0 = ',alpha0
      write(6,*) '  alpha2 = ',alpha2
      write(6,*) '       q = ',q
c
c  the lavib code has the same d value as our Born T matrices
c  but let's specify them again here so this subroutine
c  can be more easily used as a template for other T matrices
c
      d(1) = 0.d0
      d(2) = -1.d0/(2.d0 * 4.d0*atan(1.d0))
c
c  read the K matrices - search for the correct scattering
c  energy, then transform to T matrix and store it
c  we need to do this in a subroutine so that we
c  can get the dimensioning to work right
c
      do i=1,numsym
         iread = 0
         open(7,file=kmatfile(i)(1:kflen(i)), form='formatted',status
     $        ='old')
 100     continue
         if (tform .eq. 1) then
            read(7,'(4i5,d15.7,i5)',end=200,err=200) lam,l0,lm,nch,sce
     $           ,nv
         else
            read(7,'(4i5,d15.7,6x,i4)',end=200,err=200) lam,l0,lm,nch
     $           ,sce,nv
         endif

         call lavibreadk(kmat,nch)
         if (abs((sce-scatten)/scatten) .lt. 1.d-3) then
            if ((lam.gt.lmax) .or. (noabc.lt.2 .and. l0.gt.lmax)) then
               write(6,*) 'Skipping K matrix due to lmax for:'
               write(6,*) 'lam =',lam,', l0 =',l0,', lf =',lm
               go to 100
            endif
c
c  calculate k0 and k - use the value from the Kmat file
c  we'll divide by 2. in using the lavib formulae for
c  obtaining k from the vibrational levels since this
c  wants the scattering energy in Hartrees
c
            ev0 = sce/2.d0 - cmtoau*(we-2.d0*wexe*v0)*v0
            ev = sce/2.d0 - cmtoau*(we-2.d0*wexe*v)*v
            k0 = sqrt(2.d0*ev0)
            k = sqrt(2.d0*ev)
            if (lm .gt. lsymmax(lam)) then
               if (nchs.gt.0 .and. lm.le.(lam + 2*nchs - 1)) then
                  lsymmax(lam) = lm
               endif
            endif
            if (debug.ne.0) then
               write(6,*) 'Debug: Found K matrix:'
               write(6,'(4i5,d15.7,i5)') lam,l0,lm,nch,sce,nv
               write(6,*) 'Debug: lam,lsymmax(lam)=',lam,lsymmax(lam)
               call lavibmatout(kmat,nch)
            endif
            call tmatrix(kmat,tmtim,tmtre,temp,nch)
            if (debug.ne.0) then
               write(6,*) 'Converted to T matrix:'
               write(6,*) 'Re T:'
               call lavibmatout(tmtre,nch)
               write(6,*) 'Im T:'
               call lavibmatout(tmtim,nch)
            endif
            call lavibtselect(tmtre,tmtim,nch,lam,l0,lm,nv, tre,tim,ldim
     $           ,lmax,v0,v,d)

            iread = 1
         endif

         if (iread .eq. 0) go to 100

 200     close(7)

         if (iread .eq. 0) then
            write(0,*) 'the input K matrix ---',kmatfile(i)(1:kflen(i)),
     $           'does NOT include your input energy'
            stop
         endif

      enddo
c
c  now supplement Born matrix elements
c  from l=0 up to lbext - we do 0 here instead
c  of lmax because we may only have higher sym c.c. matrices
c  available - we will check for
c  non-zero T matrices in the borntmats routine and not
c  over-write any we have already filled.
c
      if (lbext .ge. 0) then
         call  borntmats(d,tre,tim,ldim,0,lbext,lsymmax,worst, k0,k
     $        ,alpha0,alpha2,q,noabc,debug,nchs)
      endif
      if (debug.ne.0) then
         write(6,*) '  Debug: k0 from Kmat file =',k0
         write(6,*) '  Debug: k from Kmat file =',k
         write(6,*) 'Debug: T after select is:'
         lend = max(lmax,lbext)
         do i=0,lend
            do jj=0,lend
               do j=0,lend
                  if (tre(i,j,jj).ne.0.d0) write(6,*) 'lam,l,l0,tre=', i
     $                 ,j,jj,tre(i,j,jj)
                  if (tim(i,j,jj).ne.0.d0) write(6,*) 'lam,l,l0,tim=', i
     $                 ,j,jj,tim(i,j,jj)
               enddo
            enddo
         enddo
      endif

      return
      end

      subroutine lavibreadk(kmat,nch)
      implicit none
      integer nch,i,j
      double precision kmat(nch,nch)
      do i=1,nch
         read(7,'(4d23.16)') (kmat(i,j),j=1,nch)
C        read(7,'(4d24.16)') (kmat(i,j),j=1,nch)
      enddo
      return
      end

      subroutine lavibmatout(kmat,nch)
      implicit none
      integer nch,i,j
      double precision kmat(nch,nch)
      do i=1,nch
         write(6,'(4d23.16)') (kmat(i,j),j=1,nch)
C        write(6,'(4d24.16)') (kmat(i,j),j=1,nch)
      enddo
      return
      end

      subroutine lavibtselect(tmtre,tmtim,nch,lam,ls,lm,nv, tre,tim,ldim
     $     ,lmax,v0,v,d)
      implicit none
      integer nch,lam,ls,lm,nv,ldim,lmax,v0,v
      double precision tre(0:ldim,0:ldim,0:ldim)
      double precision tim(0:ldim,0:ldim,0:ldim)
      double precision tmtre(nch,nch)
      double precision tmtim(nch,nch)
      double precision d(2)
      integer i,j,iv,iv0,l,l0,lind1,lind2
      if ((v .gt. nv - 1)  .or. (v0 .gt. nv - 1)) then
         write(6,*) 'Lambda,v0,v,nv = ',lam,v0,v,nv
         write(6,*) 'We found not enough open channels.  Will ',
     $        'use Born for this symmetry.'
         return
      endif
      do i=1,nch
         iv0 = (i*nv-1)/nch
         lind1 = mod(i,nch/nv)
         if (lind1.eq.0) lind1 = nch/nv
         l0 = ls + (lind1-1)*2
         do j=1,nch
            iv = (j*nv-1)/nch
            lind2 = mod(j,nch/nv)
            if (lind2.eq.0) lind2 = nch/nv
            l = ls + (lind2-1)*2
            if (iv0.eq.v0 .and. iv.eq.v) then
               tre(lam,l,l0) = d(1)*tmtre(i,j)-d(2)*tmtim(i,j)
               tim(lam,l,l0) = d(1)*tmtim(i,j)+d(2)*tmtre(i,j)
            endif
         enddo
      enddo
      return
      end
