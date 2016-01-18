c
c  This set of routines is taken verbatim out of the lavib
c  code for use in the abcvdcs code.  Well, ok, we bumped
c  up chanmax so we shouldn't have to worry about it...
c
       subroutine tmatrix(kmat,tmati,tmatr,temp,nch)
         implicit none
         integer chanmax
         parameter (chanmax=500) 
         integer nch
C numerical recipes inverse routines
         integer indx(chanmax)
         double precision temp(nch,nch)
         double precision kmat(nch,nch)
         double precision tmatr(nch,nch),tmati(nch,nch)
         integer i,j
         call mymatmul(kmat,kmat,tmati,nch,nch)
         do i=1,nch
           tmati(i,i)=tmati(i,i)+1.0d0
           end do
         call nrinverse(tmati,tmatr,temp,nch,nch,indx)
C kmat * (1+k*k)^-1  ---> tmati
         call mymatmul(kmat,tmatr,tmati,nch,nch)
C  tmati = -2 K (1 + K * K)^-1
         do i=1,nch
           do j=1,nch
             tmati(i,j)=-2.d0*tmati(i,j)
             end do
           end do
         call mymatmul(kmat,tmati,tmatr,nch,nch)
         do i=1,nch
           do j=1,nch
             tmatr(i,j)=-tmatr(i,j)
             end do
           end do
         return
         end

       subroutine mymatmul(m1,m2,m3,nsize,ndim)
         implicit none
         integer nsize,ndim
         double precision m1(ndim,ndim),m2(ndim,ndim),m3(ndim,ndim)
         integer i,j,k
         do i=1,nsize
           do j=1,nsize
             m3(i,j)=0.d0
             do k=1,nsize
               m3(i,j)=m3(i,j)+m1(i,k)*m2(k,j)
               end do
             end do
           end do
         return
         end

C
C this routine combines the numerical recipes routines
C  LUDCMP, and LUBKSB to get the inverse of a matrix.
C  I have modified the routines to do double precision
C
       subroutine nrinverse(a,b,temp,n,np,indx)
       implicit none
       integer n,np,indx(np),i,j
       double precision a(np,np),temp(np,np),b(np,np)
       double precision d
       do i=1,n
         do j=1,n
           temp(i,j)=a(i,j)
           b(i,j)=0.d0
           end do
         b(i,i)=1.d0
         end do
       call ludcmp(temp,n,np,indx,d)
       do j=1,n
         call lubksb(temp,n,np,indx,b(1,j))
         end do
       return
       end
C
C this routine carries out a LU decomposition of
C a matrix 'a' and overwrites 'a' in the process.
C  See the numerical recipes book p. 38
      SUBROUTINE ludcmp(a,n,np,indx,d)
      implicit none
      INTEGER n,np,indx(n),NMAX
      DOUBLE PRECISION d,a(np,np),TINY
      PARAMETER (NMAX=800,TINY=1.0e-20)
      INTEGER i,imax,j,k
      DOUBLE PRECISION aamax,dum,sum,vv(NMAX)
      d=1.d0
      do i=1,n
        aamax=0.d0
        do j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
          end do
        if (aamax.eq.0.d0) stop 'singular matrix in ludcmp'
        vv(i)=1.d0/aamax
        end do
      do j=1,n
        do i=1,j-1
          sum=a(i,j)
          do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
            end do
          a(i,j)=sum
          end do
        aamax=0.d0
        do i=j,n
          sum=a(i,j)
          do k=1,j-1
            sum=sum-a(i,k)*a(k,j)
            end do
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
          end do
        if (j.ne.imax)then
          do k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
            end do
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.d0)a(j,j)=TINY
        if(j.ne.n)then
          dum=1.d0/a(j,j)
          do i=j+1,n
            a(i,j)=a(i,j)*dum
            end do
        endif
        end do
      return
      END
C
C numerical recipes routine to calculate the
C inverse of a matrix.  This routine
C  requires that ludcmp be run first
C  
      SUBROUTINE lubksb(a,n,np,indx,b)
      implicit none
      INTEGER n,np,indx(n)
      DOUBLE PRECISION a(np,np),b(n)
      INTEGER i,ii,j,ll
      DOUBLE PRECISION sum
      ii=0
      do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do j=ii,i-1
            sum=sum-a(i,j)*b(j)
            end do
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
        end do
      do i=n,1,-1
        sum=b(i)
        do j=i+1,n
          sum=sum-a(i,j)*b(j)
          end do
        b(i)=sum/a(i,i)
        end do
      return
      END
