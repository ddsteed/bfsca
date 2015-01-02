C * JAN 01 1984 - AF  - FIT SIMONS-PARR-FINLAN-DUNHAM PARAMETERS TO INPUT POTENTIAL
C *
C * The position of the minimum is assumed to be known or
C * can be guessed and an accurate minimum found by minimizing the
C * chisquare
C *

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C *
C * TITLE: 72 CHARACTER TITLE
C * FORMT: FORMAT FOR READING THE POTENTIAL MAX. LENGTH 72 CHAR.
C *   RMS: INPUT VALUE FOR R AT THE MINIMUM OF THE POTENTIAL
C *        A GOOD GUESS IS HELPFUL
C *  NOFR: NUMBER OF POTENTIAL VALUES INPUT
C *  NFIT: ORDER OF S-P-F-D FIT (MAX. N OF B(N))
C * PERCNT: INITIAL PERCENTAGE CHANGE IN RMS TO MAKE FIRST PARABOLA
C *   TOL: REQUIRED ACCURACY IN RM IF OPT(4) IS TRUE
C * OPT(1): .TRUE. PRINTS THE INPUT POTENTIAL
C * OPT(2): .TRUE. PRINTS THE INPUT TO THE LINEAR EQ. SOLVER
C * OPT(3): .TRUE. PRINTS CHISQUARE FOR SOLUTION(S)
C * OPT(4): .TRUE. FINDS BEST RMS BY MINIMIZING CHISQUARE
C *          ALSO SETS OPT(3)=.TRUE.
C *
      LOGICAL OPT(4),IEND
      CHARACTER TITLE(72),FORMT(72)
      DIMENSION AVEC(36),BVEC(8),POT(110),R(110),X(16),CHI(3)
C     DIMENSION AVEC(36),BVEC(8),WKAREA(52),POT(110),R(110),X(16),CHI(3)

      DIMENSION RSAV(3)
      DATA ZERO,ONE,TWO,M,IB/0.D+00,1.D+00,2.D+00,1,8/

C* ELIMINATE UNDERFLOW MESSAGES

      READ (5,1000) (TITLE(I),I=1,72)
C     READ (5,1010) (FORMT(I),I=1,72)

C     READ (5,1020) RMS,NOFR,NFIT,PERCNT,TOL
      READ (5,*) RMS,NOFR,NFIT,PERCNT,TOL
      READ (5,1030) (OPT(I),I=1,4)
      WRITE (6,2000) (TITLE(I),I=1,72)
      WRITE (6,2010) RMS,NOFR,NFIT,PERCNT,TOL
      WRITE (6,2020) (OPT(I),I=1,4)

C * READ POTENTIAL

      IF(OPT(1)) WRITE (6,2030)
      DO I=1, NOFR
         READ (5,*) R(I),POT(I)
         IF(OPT(1)) WRITE (6,2040) R(I),POT(I)
      ENDDO
      WGT=FLOAT(NOFR-1)
      XN=FLOAT(NFIT+4)
      NFITP=NFIT+2
      NX=NFIT*2+4
      NOFA=NFITP*(NFITP+1)/2
      IF(NOFR.GE.NFITP) GO TO 20
      WRITE (6,2050)
      STOP
 20   CONTINUE
      RM=RMS
      IEND=.TRUE.
      IS=3
      IRM=3
      IF(.NOT.OPT(4)) GO TO 30
      OPT(3)=.TRUE.
      IEND=.FALSE.
      IS=1
      RMI=PERCNT*0.01D+00*RMS
      ARMI=TWO*DABS(RMI)
      RM=RMS-RMI
 30   CONTINUE

C * LOOP OVER RM

      DO IR=IS,IRM
         DO I=1,NFITP
            BVEC(I)=ZERO
         ENDDO
         DO I=1,NOFA
            AVEC(I)=ZERO
         ENDDO
         DO I=1,NOFR
            VNOW=POT(I)
            X(1)=(R(I)-RM)/R(I)
            DO J=2,NX
               X(J)=X(J-1)*X(1)
            ENDDO
            IAV=2
            IK=2
            BVEC(1)=-VNOW+BVEC(1)
            DO J=2,NFITP
               BVEC(J)=X(J)*VNOW+BVEC(J)
               KK=0
               DO K=1,IK
                  IF(K.NE.1) KK=K
                  TEMP=X(J+KK)
                  IF(K.EQ.1) TEMP=-TEMP
                  AVEC(IAV)=TEMP+AVEC(IAV)
                  IAV=IAV+1
               ENDDO
               IK=IK+1
            ENDDO
         ENDDO
         AVEC(1)=FLOAT(NOFR)
         IF(.NOT.OPT(2)) GO TO 100
         WRITE (6,2060) (AVEC(I),I=1,NOFA)
         WRITE (6,2070) (BVEC(I),I=1,NFITP)
 100     CONTINUE
         IER=0

C *
C * SOLVE LEAST SQUARES FIT BY IMSL LINEAR EQ. SOLN.
C * SOLUTION OF AM*X=BVEC -- X IS RETURNED IN BVEC
C * AVEC IS VECTOR FORM OF AM FOR POS.DEFINITE SYMMETRIC MATRIX
C * SIZE OF AVEC >= NFITP*(NFITP+1)/2
C * M= # OF INDEP. SOLNS (COLUMN DIMENSION OF BVEC) =1 HERE ALWAYS
C * IB= DIMENSIONED ROW SIZE OF BVEC (= DIM. SIZE OF VECTOR BVEC)
C * D1,D2= NOT USED BY US
C * WKAREA= MUST BE DIM. >= NFITP*(NFITP+1)/2+2*NFITP
C * IER= ERROR INDICATOR
C *
        CALL LEQT1P(AVEC,M,NFITP,BVEC,IB,D1,D2,IER)
        IF(IER.GE.64) STOP
C *
C * TEST FOR : (1+SUM OF B(I))*B(0)=EPS   (REQUIRED BY R = INFINITY)
C *
        EPS=BVEC(1)
        TEST=ONE
        DO I=3,NFITP
           BVEC(I)=BVEC(I)/BVEC(2)
           TEST=TEST+BVEC(I)
        ENDDO
        TEST=TEST*BVEC(2)-EPS
        IF(IEND) THEN 
           WRITE (6,2110)
           WRITE (11,2091) RM,EPS,BVEC(2),XN,(BVEC(I),I=3,NFITP  )
        ENDIF
        WRITE (6,2080)
        WRITE (6,2090) RM,EPS,BVEC(2),XN,(BVEC(I),I=3,NFITP)
        WRITE (6,2095) TEST
        IF(.NOT.OPT(3)) GO TO 140

C *
C * IF DESIRED, CALC. CHISQUARE
C *
        CHISQ=ZERO
        B0=BVEC(2)
        DO I=1,NOFR
           VCALC=ONE
           REXP=(R(I)-RM)/R(I)
           RN=REXP
           DO J=3,NFITP
              VCALC=VCALC+BVEC(J)*RN
              RN=RN*REXP
           ENDDO
           VCALC=VCALC*B0*REXP*REXP-EPS
           CHISQ=CHISQ+(POT(I)-VCALC)**2/DABS(VCALC)
        ENDDO
        CHISQ=CHISQ/WGT
        CHI(IR)=CHISQ
        WRITE (6,2100) CHISQ
 140    CONTINUE
        RSAV(IR)=RM
        RM=RM+RMI
      ENDDO
      IF(IEND) GO TO 160

C *
C * FIT CHISQARE VS RM TO PARABOLA AND CALC FINAL RM
C *
      RAT=(RSAV(1)-RSAV(3))/(RSAV(2)-RSAV(3))
      R1S=RSAV(1)**2
      R2S=RSAV(2)**2
      R3S=RSAV(3)**2
      A=(CHI(1)-CHI(3)-RAT*(CHI(2)-CHI(3)))/(R1S-R3S-RAT*(R2S-R3S))
      B=(CHI(2)-CHI(3)-A*(R2S-R3S))/(RSAV(2)-RSAV(3))
      RM=-B/(TWO*A)
      RI=RM-RSAV(3)
      ARI=DABS(RI)
      IF(ARI.GT.ARMI) RM=RSAV(3)+RI/ARI*ARMI
      CHI(1)=CHI(2)
      CHI(2)=CHI(3)
      RSAV(1)=RSAV(2)
      RSAV(2)=RSAV(3)
      IS=3
      IF(ARI.LT.TOL) IEND=.TRUE.
      GO TO 30
 160  CONTINUE
      STOP
 1000 FORMAT(72A1)
 1010 FORMAT(72A1)
 1020 FORMAT(F10.5,2I5,F5.1,D10.3)
 1030 FORMAT(4L2)
 2000 FORMAT(1X,72A1//)
 2010 FORMAT(1X,'INPUT VALUES:     RM',6X,'NOFR NFIT    PERCNT  TOL'/
     $     14X,F10.5,I6,I5,F10.2,2X,1PD10.3//)
 2020 FORMAT(1X,'WRTPOT WRTVEC DOCHI DOLOOP'/1X,L4,L7,L6,L7//)
 2030 FORMAT(6X,'R',10X,'V'/)
 2040 FORMAT(1X,F10.5,2X,1PD15.7)
 2050 FORMAT(//1X,25('*')/1X,'NOT ENOUGH POINTS FOR FIT'/1X,25('*'))
 2060 FORMAT(1X,'AVEC =',10(1X,D10.3)/3(8X,D10.3,9D11.3))
 2070 FORMAT(1X,'BVEC =',7(1X,1PD10.3))
 2080 FORMAT(/1X,'SIMONS-PARR-FINLAN-DUNHAM EXPANSION')
 2090 FORMAT(5X,'RM',10X,'EPS',8X,'COEFFICIENTS'/1X,10(1PD13.6))
 2091 format(4(1pe18.8))
 2095 FORMAT(1X,'TEST =',1PD11.3)
 2100 FORMAT(1X,'CHISQUARE =',1PD10.3)
 2110 FORMAT(/' FINAL RESULTS'/)
      END

C *-
C * IMSL ROUTINE NAME   - LEQT1P
C *
C *
C * COMPUTER            - VAX/DOUBLE
C *
C * LATEST REVISION     - JANUARY 1, 1978
C *
C * PURPOSE             - LINEAR EQUATION SOLUTION - POSITIVE DEFINITE
C *                       MATRIX - SYMMETRIC STORAGE MODE - SPACE
C *                       ECONOMIZER SOLUTION
C *
C * USAGE               - CALL LEQT1P (A,M,N,B,IB,D1,D2,IER)
C *
C * ARGUMENTS    A      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING THE
C *                       N BY N COEFFICIENT MATRIX OF THE EQUATION
C *                       AX = B. A IS A POSITIVE DEFINITE SYMMETRIC
C *                       MATRIX STORED IN SYMMETRIC STORAGE MODE.
C *                       ON OUTPUT, A IS REPLACED BY THE LOWER
C *                       TRIANGULAR MATRIX L WHERE A = L*L-TRANSPOSE.
C *                       L IS STORED IN SYMMETRIC STORAGE MODE WITH
C *                       THE DIAGONAL ELEMENTS OF L IN RECIPROCAL
C *                       FORM.
C *              M      - NUMBER OF RIGHT HAND SIDES (COLUMNS IN B).
C *                       (INPUT)
C *              N      - ORDER OF A AND NUMBER OF ROWS IN B. (INPUT)
C *              B      - INPUT MATRIX OF DIMENSION N BY M CONTAINING
C *                       THE RIGHT-HAND SIDES OF THE EQUATION
C *                       AX = B.
C *                       ON OUTPUT, THE N BY M SOLUTION MATRIX X
C *                       REPLACES B.
C *              IB     - ROW DIMENSION OF B EXACTLY AS SPECIFIED IN
C *                       THE DIMENSION STATEMENT IN THE CALLING
C *                       PROGRAM. (INPUT)
C *              D1,D2  - COMPONENTS OF THE DETERMINANT OF A.
C *                       DETERMINANT(A) = D1*2.**D2. (OUTPUT)
C *              IER    - ERROR PARAMETER. (OUTPUT)
C *                       TERMINAL ERROR
C *                       IER = 129 INDICATES THAT THE INPUT MATRIX
C *                       A IS ALGORITHMICALLY NOT POSITIVE
C *                       DEFINITE. (SEE THE CHAPTER L PRELUDE).
C *
C * PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C *                     - SINGLE/H36,H48,H60
C *
C * REQD. IMSL ROUTINES - LUDECP,LUELMP,UERTST,UGETIO
C *
C * NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C *                       CONVENTIONS IS AVAILABLE IN THE MANUAL
C *                       INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C *
C * COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C *
C * WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C *                       APPLIED TO THIS CODE. NO OTHER WARRANTY,
C *                       EXPRESSED OR IMPLIED, IS APPLICABLE.
C *
C *

      SUBROUTINE LEQT1P (A,M,N,B,IB,D1,D2,IER)

      DIMENSION          A(1),B(IB,1)
      DOUBLE PRECISION   A,B,D1,D2
C * FIRST EXECUTABLE STATEMENT
C * INITIALIZE IER
      IER = 0
C * DECOMPOSE A
      CALL LUDECP (A,A,N,D1,D2,IER)
      IF (IER.NE.0) GO TO 9000
C * PERFORM ELIMINATION
      DO I = 1,M
         CALL LUELMP (A,B(1,I),N,B(1,I))
      ENDDO
      GO TO 9005
 9000 CONTINUE
      stop 'error in leqt1p'
 9005 RETURN
      END

C *-
C *  IMSL ROUTINE NAME   - LUDECP
C *
C *
C * COMPUTER            - VAX/DOUBLE
C *
C * LATEST REVISION     - JANUARY 1, 1978
C *
C * PURPOSE             - DECOMPOSITION OF A POSITIVE DEFINITE MATRIX -
C *                       SYMMETRIC STORAGE MODE
C *
C * USAGE               - CALL LUDECP (A,UL,N,D1,D2,IER)
C *
C * ARGUMENTS    A      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING
C *                       THE N BY N POSITIVE DEFINITE SYMMETRIC
C *                       MATRIX STORED IN SYMMETRIC STORAGE MODE.
C *              UL     - OUTPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING
C *                       THE DECOMPOSED MATRIX L SUCH THAT A = L*
C *                       L-TRANSPOSE. L IS STORED IN SYMMETRIC
C *                       STORAGE MODE. THE DIAGONAL OF L CONTAINS THE
C *                       RECIPROCALS OF THE ACTUAL DIAGONAL ELEMENTS.
C *              N      - ORDER OF A. (INPUT)
C *              D1,D2  - COMPONENTS OF THE DETERMINANT OF A.
C *                       DETERMINANT(A) = D1*2.**D2. (OUTPUT)
C *              IER    - ERROR PARAMETER. (OUTPUT)
C *                       TERMINAL ERROR
C *                       IER = 129 INDICATES THAT MATRIX A IS
C *                       ALGORITHMICALLY NOT POSITIVE DEFINITE.
C *                       (SEE THE CHAPTER L PRELUDE).
C *
C * PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C *                     - SINGLE/H36,H48,H60
C *
C * REQD. IMSL ROUTINES - UERTST,UGETIO
C *
C * NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C *                       CONVENTIONS IS AVAILABLE IN THE MANUAL
C *                       INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C *
C * COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C *
C * WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C *                       APPLIED TO THIS CODE. NO OTHER WARRANTY,
C *                       EXPRESSED OR IMPLIED, IS APPLICABLE.
C *

      SUBROUTINE LUDECP (A,UL,N,D1,D2,IER)

      DIMENSION          A(1),UL(1)
      DOUBLE PRECISION   A,UL,D1,D2,ZERO,ONE,FOUR,SIXTN,SIXTH,X,RN
      DATA               ZERO,ONE,FOUR,SIXTN,SIXTH/
     *                   0.0D0,1.D0,4.D0,16.D0,.0625D0/
C * FIRST EXECUTABLE STATEMENT
      D1=ONE
      D2=ZERO
      RN = ONE/(N*SIXTN)
      IP = 1
      IER=0
      DO I = 1,N
         IQ = IP
         IR = 1
         DO J = 1,I
            X = A(IP)
            IF (J .EQ. 1) GO TO 10
            DO K=IQ,IP1
               X = X - UL(K) * UL(IR)
               IR = IR+1
            ENDDO
 10         IF (I.NE.J) GO TO 30
            D1 = D1*X
            IF (A(IP) + X*RN .LE. A(IP)) GO TO 50
 15         IF (DABS(D1).LE.ONE) GO TO 20
            D1 = D1 * SIXTH
            D2 = D2 + FOUR
            GO TO 15
 20         IF (DABS(D1) .GE. SIXTH) GO TO 25
            D1 = D1 * SIXTN
            D2 = D2 - FOUR
            GO TO 20
 25         UL(IP) = ONE/DSQRT(X)
            GO TO 35
 30         UL(IP) = X * UL(IR)
 35         IP1 = IP
            IP = IP+1
            IR = IR+1
         ENDDO
      ENDDO
      GO TO 9005
 50   IER = 129
      stop 'error in ludecp'
 9005 RETURN
      END

C *-
C * IMSL ROUTINE NAME   - LUELMP
C *
C * COMPUTER            - VAX/DOUBLE
C *
C * LATEST REVISION     - JANUARY 1, 1978
C *
C * PURPOSE             - ELIMINATION PART OF THE SOLUTION OF AX=B -
C *                       POSITIVE DEFINITE MATRIX - SYMMETRIC
C *                       STORAGE MODE
C *
C * USAGE               - CALL LUELMP (A,B,N,X)
C *
C * ARGUMENTS    A      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING
C *                       THE N BY N MATRIX L WHERE A = L*L-TRANSPOSE.
C *                       L IS A LOWER TRIANGULAR MATRIX STORED IN
C *                       SYMMETRIC STORAGE MODE. THE MAIN DIAGONAL
C *                       ELEMENTS OF L ARE STORED IN RECIPROCAL
C *                       FORM. MATRIX L MAY BE OBTAINED FROM IMSL
C *                       ROUTINE LUDECP.
C *              B      - VECTOR OF LENGTH N CONTAINING THE RIGHT HAND
C *                       SIDE OF THE EQUATION AX = B. (INPUT)
C *              N      - ORDER OF A AND THE LENGTH OF B AND X. (INPUT)
C *              X      - VECTOR OF LENGTH N CONTAINING THE SOLUTION TO
C *                       THE EQUATION AX = B. (OUTPUT)
C *
C * PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C *                     - SINGLE/H36,H48,H60
C *
C * REQD. IMSL ROUTINES - NONE REQUIRED
C *
C * NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C *                       CONVENTIONS IS AVAILABLE IN THE MANUAL
C *                       INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C *
C * COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C *
C * WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C *                       APPLIED TO THIS CODE. NO OTHER WARRANTY,
C *                       EXPRESSED OR IMPLIED, IS APPLICABLE.
C *

      SUBROUTINE LUELMP (A,B,N,X)

      DIMENSION          A(1),B(1),X(1)
      DOUBLE PRECISION   A,B,X,T,ZERO
      DATA               ZERO/0.0D0/

C * FIRST EXECUTABLE STATEMENT
C * SOLUTION OF LY = B

      IP=1
      IW = 0
      DO I=1,N
         T=B(I)
         IM1 = I-1
         IF (IW .EQ. 0) GO TO 9
         IP=IP+IW-1
         DO K=IW,IM1
            T = T-A(IP)*X(K)
            IP=IP+1
         ENDDO
         GO TO 10
    9    IF (T .NE. ZERO) IW = I
         IP = IP+IM1
   10    X(I)=T*A(IP)
         IP=IP+1
      ENDDO

C * SOLUTION OF UX = Y

      N1 = N+1
      DO I = 1,N
         II = N1-I
         IP=IP-1
         IS=IP
         IQ=II+1
         T=X(II)
         IF (N.LT.IQ) GO TO 25
         KK = N
         DO K=IQ,N
            T = T - A(IS) * X(KK)
            KK = KK-1
            IS = IS-KK
         ENDDO
   25    X(II)=T*A(IS)
      ENDDO
      RETURN
      END
