C * SEP 17 2005 - RDS - CALCULATE MULTIPOLE MOMENTS

C *         D_dipole(r;R) = - (r**2)*V_static(Lambda=1;r,R)
C *     Q_quatrupole(r;R) = - (r**3)*V_static(Lambda=2;r,R)
C *         A_alpha0(r;R) = - (2*r**4)*V_pola(Lambda=0;r,R)
C *         A_alpha2(r;R) = - (2*r**4)*V_pola(Lambda=2;r,R)

C *  For Lambda > 2 :
C *       O_octopole(r;R) = - (r**4)*V_static(Lambda=3;r,R)
C *       H_hexapole(r;R) = - (r**5)*V_static(Lambda=4;r,R)
C *               ....................
C *       G(r;R) = - (r**L1)*V_static(Lambda;r,R),  L1=Lambda+1 .

      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER INP, OUT, LAB, LAMBDA
      REAL*8 RM, MPMOM, R, V
      DIMENSION X(1000),Y(1000),X1(40),X2(40),X3(40),C(40),D(40)
      DIMENSION TITLE(72)

C *
C * Read in data :
C *     M = 1, Calc. Alpha_0 -- isotropic polarizability
C *       = 2, Calc. Alpha_2 -- anisotropic polarizability

C *       = 3, Calc. D -- dipole moment
C *       = 4, Calc. Q -- quadrupole moment

C * V_static potential are readed in from fort.2;
C * V_polar. potential are readed in from fort.3.

 12   FORMAT (I5)
 13   FORMAT (2I5)
 21   FORMAT (2D15.6)
 27   FORMAT (F10.5,D24.16)
 30   FORMAT (F10.5,F10.5)
 40   FORMAT (3F10.5)
 50   FORMAT (72a1)

      READ (5,*) M

      OUT = 6

      IF (M .EQ. 3 .OR. M .EQ. 4) THEN
         INP = 2
      ELSE IF (M .EQ. 1 .OR. M .EQ. 2) THEN
         INP = 3
      ELSE
         WRITE (OUT,*) " READ FILE DOES NOT EXIST!"
         STOP
      ENDIF

C *- READ auxiliary data
 10   READ (INP,50) TITLE
      READ (INP,12) NRG
      DO I=1,NRG
         READ (INP,40) X1(I),X2(I),X3(I)
      ENDDO 

      READ (INP,12) NATOM
      RM = 0.0D0
      DO I = 1, NATOM
         READ (INP,21) C(I),D(I)
         RM = RM + DABS(D(I))
      ENDDO 
      
      IF (M .EQ. 1) THEN
         LAMBDA = 0
      ELSE IF (M .EQ. 2) THEN
         LAMBDA = 2
      ELSE IF (M .EQ. 3) THEN
         LAMBDA = 1
      ELSE IF (M .EQ. 4) THEN
         LAMBDA = 2
      ENDIF

C *- READ R AND V UNTIL THE DATA TO CALCULATE MOMENT ARE FOUND
      IREAD = 0

 100  READ(INP,13,END=800) LAB, NX
      DO I = 1, NX
         READ(INP,27,END=900) R, V
      ENDDO

      IF (LAMBDA .EQ. LAB) THEN
         IREAD = 1
         IF (M .EQ. 1 .OR. M .EQ. 2) THEN
            MPMOM = -2.0D0*R**4*V
         ELSE IF (M .EQ. 3 .OR. M .EQ. 4) THEN
            MPMOM = -R**(LAB + 1)*V
         ENDIF
      ENDIF
      GO TO 850

 800  IREAD = 1
 850  IF (IREAD .EQ. 0) GOTO 100

      WRITE (OUT,30) RM, MPMOM
      GOTO 990

 900  WRITE (OUT,*) "DATA ARE NOT SUFFICIENT IN READ FILE!"
 990  stop
      end
