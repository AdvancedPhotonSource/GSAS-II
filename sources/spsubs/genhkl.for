      SUBROUTINE GENHKL(XH,NSYM,RT,ICEN,NCV,CEN,JHK,HKL,IABSNT,MULP)

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 
!
! The HKL index generation S.R. - Generate equivs. even for spgp absent refl.

! Input data list

      REAL*4        XH(4)               ! Input Miller indices
      INTEGER*4     NSYM                ! Number of symmetry matrices
      INTEGER*4     RT(3,5,24)          ! The symmetry matrices
      INTEGER*4     ICEN                ! Flag indicating 1bar
      INTEGER*4     NCV                 ! The number of lattice centering vectors
      REAL*4        CEN(3,NCV)            ! The lattice centering vectors

!Output data list

      INTEGER*4     JHK                 ! Number of equivalent indices generated
      REAL*4        HKL(4,24)           ! The generated Miller indices
      INTEGER*4     IHKL(4,24)          ! The generated Miller indices
      INTEGER*4     IABSNT              ! Space group absence flag
      INTEGER*4     MULP                ! Multiplicity for powder line intensities

!CODE

      MULP = 0
      JHK = 1                                              ! Set generated reflection count to 1
      IABSNT = 0                                          ! Assume NOT Space Group Extinct

      IF ( ABS(NINT(XH(1))-XH(1))+
     1  ABS(NINT(XH(2))-XH(2))+
     1  ABS(NINT(XH(3))-XH(3)).GT.0.05 ) THEN                  !Check for non-integral indices, leave a bit of slop
        IABSNT = 1                                        !Non-integral indices are not allowed
      ELSE
        DO I=2,NCV                                          ! First check for lattice type extinctions
          K = 0
          DO J=1,3
            K = K+NINT(XH(J)*CEN(J,I)*12.0)
          END DO
          IF ( MOD(K,12).NE.0 ) IABSNT=1
        END DO
      END IF

      I = 1
      DO WHILE ( I.LE.NSYM )                                    ! Generate the equivalent index set
        DO J=1,4
          IHKL(J,JHK) = 0.0
          DO K=1,3
            IHKL(J,JHK) = IHKL(J,JHK)+IFIX(XH(K))*RT(K,J,I)
          END DO
        END DO

        NEW = 1
        NEWX = 1
        IF ( JHK.GT.1 ) THEN                                    ! Check for previous generation of this index
          J = 1
          DO WHILE ( J.LT.JHK .AND. NEW.EQ.1 )
            IF ( IHKL(1,J).EQ.IHKL(1,JHK) ) THEN
              IF ( IHKL(2,J).EQ.IHKL(2,JHK) ) THEN
                IF ( IHKL(3,J).EQ.IHKL(3,JHK) ) THEN
                  NEW = 0
                  NEWX = 0
                  IF ( MOD(IHKL(4,JHK)-IHKL(4,J)+960,12).NE.0 ) THEN
                    IABSNT = 1
                  END IF
                END IF
              END IF
            END IF
            IF ( NEW.EQ.1 ) THEN
              IF ( IHKL(1,J).EQ.-IHKL(1,JHK) ) THEN                  ! Check -h,k,l)
                IF ( IHKL(2,J).EQ.-IHKL(2,JHK) ) THEN
                  IF ( IHKL(3,J).EQ.-IHKL(3,JHK) ) THEN
                    NEWX = 0
                    IF ( ICEN.GT.0 ) THEN
                      NEW = 0
                      IF ( MOD(IHKL(4,JHK)-IHKL(4,J)+960,12).NE.0 ) THEN
                        IABSNT = 1
                      END IF
                    END IF
                  END IF
                END IF
              END IF
            END IF
            J = J+1
          END DO
        END IF
        MULP = MULP+NEWX
        IF ( NEW.EQ.1 ) THEN
!         IMAT(JHK) = I
          JHK = JHK+NEW
        END IF
        I = I+1
      END DO
      JHK = JHK-1

      DO I=1,JHK                                          ! All is OK. Convert to F.P.
        HKL(1,I) = IHKL(1,I)
        HKL(2,I) = IHKL(2,I)
        HKL(3,I) = IHKL(3,I)
        HKL(4,I) = MOD(FLOAT(IHKL(4,I))/12.0,1.0)
      END DO
      RETURN
      END
