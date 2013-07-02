C * AUG 20 2005 - RDS - INTERPOLATE POLARIZATION POTENTIAL BY CUBIC FUNCTION

      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION X(1000), XP(1000), Y(1000)
      DIMENSION Y2(1000)
      DIMENSION RI(10), RF(10), RS(10), C1(10), C2(10), LAB(30)
      DIMENSION TITLE(72)

C *
C * Read in data 
C *
C * Card 1:  x1, x2 --- starting and ending cubic spline points of x dimension
C *            nin  --- number of data pairs in the input data set
C *
C * Card 2:  Specific xp values for which you want to have y=f(x). These y's
C *          are calculated from your fitted results y(i)'s. xp(i)
C *
C * Card 3:  The (input) data to be fitted, x(i), y(i), i=1, nin.
C *
C *
      READ (5,*) X1, X2, NIN
C *
C * Read potential mesh
C * 
      READ (3,600) TITLE
      READ (5,*)   LABDA
      WRITE (8,600) TITLE
      READ (3,'(I5)')  NRG
      WRITE (8,'(I5)') NRG
      DO I = 1, NRG
         READ (3,'(3F10.5)')   RI(I),  RF(I), RS(I)
         WRITE (8,'(3F10.5)')  RI(I),  RF(I), RS(I)
      ENDDO
      READ (3,'(I5)') NATOM
      WRITE (8,'(I5)') NATOM
      DO I = 1, NATOM
         READ (3,'(2D15.6)')  C1(I), C2(I)
         WRITE (8,'(2D15.6)') C1(I), C2(I)
      ENDDO

      READ (3,'(2I5)') LAB(K), NX
      DO I = 1,NX
         READ (3,'(F10.5,D24.16)') XP(I), AX
      ENDDO

C * grid ready, read data x,y (which will be fitted)

      NN = NIN - 1
      WRITE (6,625) X1,X2,NIN

      DO I = 1, NIN
         READ (5,*)       X(I), Y(I)
         WRITE (6,630) I, X(I), Y(I)
      ENDDO

C *
C * Using cubic spline to INTERPOLATE data
C *
      WRITE (6,640) NIN
      CALL SPLINE(X,Y,NIN,YP1,YPN,Y2)

      DO I = 1, NIN
         WRITE (6,630) I, X(I), Y(I), Y2(I)
      ENDDO

C * Calculate y=f(x) for a specific xp(i) 

      WRITE (6,650) NX
      WRITE (8,'(2I5)')  LABDA, NX
      WRITE (9,'(2I5)')  LABDA, NX

      DO I = 1, NX
         X3 = XP(I)
         CALL SPLINT(X,Y,Y2,NIN,X3,Y3)

C * DON'T do extrapolation
         IF (X3 .LT. X1 .OR. X3 .GT. X2) Y3  =  0.0

         WRITE (6,630) I, X3, Y3
         WRITE (8,'(F10.5,D24.16)')  X3, Y3
         WRITE (9,'(F10.5,D24.16)')  X3, Y3
      ENDDO

 600  format(72a1)
 625  format(/,' Using cubic spline fitter to fit a set of', ' dat
     $     a',//' Interpolating from x1 =',f8.4,'  to x2 =', f8.4,/
     $     ,'         Input n =',i4,' data ',//, 3x
     $     ,'  i      x(i)          y(i) '/)
 630  format(3x,i3,1x,f10.5,1pe18.8,1pe18.8)
 640  format(//,'  *****  Output from spline fitter  ***** ',//,
     $ '         Output n  = ',i4,' data ',//, 3x,
     $ '  i      x(i)          y(i)         y2(i)(second derivative) '/)
 650  format(//,'  *****  Output for specific (x, y)s  ***** ',//,
     $     '         Output n  = ',i4,' data ',//, 3x
     $     ,'  i      x(i)          y(i) '/)

      stop
      end
