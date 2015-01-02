C * MAY 25 2006 - RDS - AUGMENT THE DIMENSION OF LAM AND V TO COPE WITH E-N2
C * JAN 01 1984 - AF  - READS IN THE POTENTIAL (VLAM) PRODUCED BY VLAM
C *
C * file 10 is input  : file 11 output

      program wlut
      implicit none

      integer lammax,ptsmax
      parameter(lammax=40,ptsmax=600)

      DOUBLE PRECISION V(ptsmax,lammax)
      double precision R(ptsmax),RNUC(3),rsep
      integer LAM(lammax)
      CHARACTER DUM(4)
      integer nlam,l,i,nreg,nnuc,nofr

      READ(10,1010) DUM
      READ(10,1060) NREG
      DO I=1,NREG
         READ(10,1010) DUM
      ENDDO
      READ(10,1060) NNUC
      DO I=1,NNUC
         READ(10,1070) RNUC(I)
      ENDDO
      RSEP=ABS(RNUC(1))+ABS(RNUC(NNUC))
      L=1
 25   CONTINUE
      READ(10,1020,END=30) LAM(L),NOFR
      DO I=1,NOFR
         READ(10,1030) R(I),V(I,L)
      ENDDO
      L=L+1
      GO TO 25
 30   CONTINUE
      NLAM=L-1
      WRITE(11,1041) RSEP
      DO I=1,NOFR
         DO L=1,NLAM
            WRITE(11,1051) LAM(L),R(I),V(I,L)
         ENDDO
      ENDDO
      STOP
 1010 FORMAT(4A1)
 1020 FORMAT(2I5)
 1030 FORMAT(F10.5,D24.16)
C1030 FORMAT(F10.5,D23.16)
 1041 FORMAT(5X,'RSEP=',F10.5)
 1051 FORMAT(I5,F10.5,D24.16)
C1051 FORMAT(I5,F10.5,D23.16)
 1060 FORMAT(I5)
 1070 FORMAT(15X,D15.6)
      END
