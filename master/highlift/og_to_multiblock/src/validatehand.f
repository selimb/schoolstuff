      SUBROUTINE VALIDATEHAND(NB, IJK)
      use globals, only: X3D, ni, nj, nk
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NB
      INTEGER, INTENT(OUT) :: IJK
      INTEGER, PARAMETER :: RIGHT = 1
      INTEGER, PARAMETER :: LEFT = -1
      INTEGER :: IJK2, HAND
      INTEGER :: I, J, K
      REAL*8, DIMENSION(3) :: VI, VJ, VK, CROSS
      IJK2 = 0
      DO 11 K = 1, nk(NB) - 1
      DO 11 J = 1, nj(NB) - 1
      DO 11 I = 1, ni(NB) - 1
!         VI(:) = X(NB,I+1,J,K,:) - X(NB,I,J,K,:)
!         VJ(:) = X(NB,I,J+1,K,:) - X(NB,I,J,K,:)
!         VK(:) = X(NB,I,J,K+1,:) - X(NB,I,J,K,:)
          VI(:) = X3D(:,I+1,J,K) - X3D(:,I,J,K)
          VJ(:) = X3D(:,I,J+1,K) - X3D(:,I,J,K)
          VK(:) = X3D(:,I,J,K+1) - X3D(:,I,J,K)
          CROSS(1) = VI(2)*VJ(3) - VI(3)*VJ(2)
          CROSS(2) = VI(3)*VJ(1) - VI(1)*VJ(3)
          CROSS(3) = VI(1)*VJ(2) - VI(2)*VJ(1)
          IF (DOT_PRODUCT(CROSS, VK) .GT. 0) THEN
              HAND = RIGHT
          ELSE
              HAND = LEFT
          END IF
          IF (IJK2 .EQ. 0) THEN
              IJK2 = HAND
          ELSE IF (IJK2 .NE. HAND) THEN
              write(*,*) "INCONSTANT CALCULATED HAND FOR BLOCK", NB
              IJK = -666
              RETURN
          END IF
   11 CONTINUE
      IJK = IJK2
!     IF (IJK2.NE.IJK) THEN
!         write(*,*) "HANDEDNESS FOR BLOCK", NB, "DOES NOT MATCH"
!         write(*,*) "INPUT, CALCULATED", IJK, IJK2
!         VALID = .FALSE.
!     END IF
      END SUBROUTINE
