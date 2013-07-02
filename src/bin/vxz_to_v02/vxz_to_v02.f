C * vxz_to_v02
C *
C * v0 = (vz + 2vx)/3
C * v2 = (vz -  vx)/3 * 2
C*
      PROGRAM MAIN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      READ(5,*) NUM
      DO I = 1, NUM
         READ(7,*) RE, VX
         READ(8,*) RE, VZ

         READ(9,*) RE, AX
         READ(10,*)RE, AZ

         V0 = (VZ + 2*VX)/3.0
         V2 = (VZ -   VX)/3.0 * 2.0

         A0 = (AZ + 2*AX)/3.0
         A2 = (AZ -   AX)/3.0 * 2.0

         WRITE(11,100) RE, V0
         WRITE(12,100) RE, V2
         WRITE(13,100) RE, A0
         WRITE(14,100) RE, A2
      ENDDO
 100  FORMAT(F10.6,3X,E34.26)

      END

