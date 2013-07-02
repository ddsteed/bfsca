C=======================================================================
C     PML = P[l,m,cos(x)]
C
      DOUBLE PRECISION FUNCTION PML(LX,MX,XF)
C
C     THIS SUBFUNCTION CALCULATES ASSOCIATED LEGENDRE POLYNOMIALS OF
C     ANY ORDER P SUB LX SUPER MX OF XF
C     THE ARGUMENT XF MUST BE THE ANGLE IN RADIANS
C     THE CODE USES EXPLICIT FORMULAS FOR L THROUGH 8
C     NOTE THAT FOR L GREATER THAN 4 ONLY M = 0,1,AND 2 ARE INCLUDED
C     THE RECURRENCE FORMULAS APPLY ONLY TO M = 0,1, AND 2
C     THE CONVENTIONS OF MCLEAN AND YOSHIMINE ARE USED FOR SIGNS
C
      IMPLICIT REAL*8 (A-H,O-Z)
      X= DCOS(XF)
      XP = DSIN(XF)
      X2 = X * X
      X3 = X2 * X
      X4 = X3 * X
      XP2 = XP * XP
      XP3 = XP2 * XP
      XP4 = XP3 * XP
      K = LX+1
      M = MX+1
C     SPECIAL CASES FOR SMALL ARGUMENTS
      IF(DABS(X-1.0D+00) .LE. 1.0D-06)  GO TO 215
C
C     K = ORDER OF POLYNOMIAL + 1
C     IF ORDER GREATER THAN 8 USE STANDARD RECURRENCE RELATION
      IF (K .GT. 9) GO TO 110
C
      GO TO (10,20,30,40,50,70,80,90,100),K
   10 PML = 1.0D+00
      GO TO 60
 20   GO TO (21,22),M
 21   PML = X
      GO TO 60
 22   PML = XP
      GO TO 60
   30 GO TO (31,32,33),M
   31 PML = 1.5D0*X2-0.5D0
      GO TO 60
 32   PML = 3.0D+00 * X * XP
      GO TO 60
 33   PML = 3.0D+00 * XP2
      GO TO 60
 40   GO TO (41,42,43,44),M
 41   PML = 2.5D0 * X3-1.5D0 * X
      GO TO 60
 42   PML = XP * (7.5D+00 * X2-1.5D+00)
      GO TO 60
 43   PML = 15.0D+00 * XP2 * X
      GO TO 60
 44   PML = 15.0D+00 * XP3
      GO TO 60
 50   GO TO (51,52,53,54,55),M
 51   PML = (35.0D0*X4-30.D0*X2+3.0D0)/8.0D0
      GO TO 60
 52   PML = XP * (35.0D+00 * X3-15.0D+00 * X)/2.0D+00
      GO TO 60
  53  PML = XP2 * (105.0D+00 * X2-15.0D+00)/2.0D+00
      GO TO 60
 54   PML = 105.0D+00 * XP3 * X
      GO TO 60
 55   PML=105.0D+00 * XP4
 70   GO TO (71,72,73),M
 71   PML=(63.D0*X4*X-70.D0*X3+15.D0*X)/8.D0
      GO TO 60
   72 PML=XP*(315.D0*X4-210.D0*X2+15.D0)/8.D0
      GO TO 60
   73 PML = (105.0D+00 * XP2 * (3.0D+00 * X2-1.0D+00) * X)/2.0D+00
      GO TO 60
 80   GO TO (81,82,83),M
 81   PML=(231.D0*X2*X4-315.D0*X4+105.D0*X2-5.0D0)/16.D0
      GO TO 60
   82 PML=X*XP*(693.D0*X4-630.D0*X2+105.D0)/8.D0
      GO TO 60
 83   PML=(105.0D+00*XP2*(33.0D+00*X4-18.0D+00*X2+1.0D+00))/8.0D+00
      GO TO 60
 90   GO TO (91,92,93),M
 91   PML=(429.D0*X4*X3-693.0D0*X4*X+315.0D0*X3-35.0D0*X)/16.D0
      GO TO 60
   92 PML=XP*((3003.D0*X4-3465.D0*X2+945.D0)*X2-35.D0)/16.D0
      GO TO 60
 93   PML=63.D0*(143.D0*X4-110.D0*X2+15.D0)*X*(XP2/8.D0)
      GO TO 60
 100  GO TO (101,102,103),M
 101  PML=(6435.D0*X4*X4-12012.D0*X4*X2+6930.D0*X4-1260.D0*X2+35.D0)
     1/128.D0
      GO TO 60
  102 PML=X*XP*((6435.D0*X4-9009.D0*X2+3465.D0)*X2-315.D0)/16.D0
      GO TO 60
 103  PML=630.D0*(143.D0*X4*(X2-1.D0)+33.D0*X2-1.0D0)*(XP2/32.D0)
      GO TO 60
  110 LTEMP=8
C
C     SET PML EQUAL TO PM8 AND PMO EQUAL TO PM7 INITIALLY
C
      GO TO (200,201,202), M
  200 PML=(6435.D0*X4*X4-12012.D0*X4*X2+6930.D0*X4-1260.D0*X2+35.D0)
     1/128.D0
      PMO=(429.D0*X4*X3-693.0D0*X4*X+315.0D0*X3-35.0D0*X)/16.D0
      GO TO 210
  201 PML=X*XP*((6435.D0*X4-9009.D0*X2+3465.D0)*X2-315.D0)/16.D0
      PMO=XP*((3003.D0*X4-3465.D0*X2+945.D0)*X2-35.D0)/16.D0
      GO TO 210
  202 PML=XP2*630.D0*(143.D0*X4*(X2-1.D0)+33.D0*X2-1.D0)/32.D0
      PMO=XP2*63.D0*(143.D0*X4-110.D0*X2+15.D0)*X/8.D0
  210 FACTR1=DFLOAT(2*LTEMP+1)/DFLOAT(LTEMP-MX+1)
      FACTR2=DFLOAT(LTEMP+MX)/DFLOAT(LTEMP-MX+1)
C     PMLNEW EQUALS P OF ORDER LTEMP+1
      PMLNEW=FACTR1*X *PML-FACTR2*PMO
      PMO=PML
      PML=PMLNEW
      LTEMP=LTEMP+1
      IF (LTEMP.EQ.LX) GO TO 60
      GO TO 210
  215 PML=1.0D+00
      IF(MX .NE. 0) PML=0.0D+00
   60 CONTINUE
      RETURN
      END