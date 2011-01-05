      SUBROUTINE SGTRCF(M,RT,N,M2,LCENT,LAUENO,IER)

!Purpose:      Search for translation conflicts

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 

C       This program was developed for
C                    The Division of Chemistry
C                               of
C               The National Research Council of Canada
C                               by
C       Allen C. Larson, 14 Cerrado Loop, Santa Fe, NM 87505-8832, USA

      INTEGER*4     M                   !
      REAL*4        RT(5,4,25)          !Matrices being generated
      INTEGER*4     N                   !Sequence no. of matrix 1
      INTEGER*4     M2                  !Sequence no. of matrix 2
      INTEGER*4     LCENT               !Number of Lattice centering vectors
      INTEGER*4     LAUENO              !Laue group flag
      INTEGER*4     IER                 !Error flag

      DIMENSION     ICENV(3,5),NCVT(7),JCVT(7)
      DATA ICENV/0,0,0, 0,6,6, 6,0,6, 6,6,0, 6,6,6/
      DATA NCVT/1,2,3,4,5,4,1/
      DATA JCVT/1,1,2,3,4,1,1/

      IER = 0
      IRX = 12.0*MOD((RT(1,4,N)-RT(1,4,M2)),1.0)
      IRY = 12.0*MOD((RT(2,4,N)-RT(2,4,M2)),1.0)
      IRZ = 12.0*MOD((RT(3,4,N)-RT(3,4,M2)),1.0)
      NCV = NCVT(LCENT)
      JCV = JCVT(LCENT)

      ICV = 1-JCV
      TOTTR = 1
      DO WHILE ( TOTTR.NE.0 .AND. ICV.LT.NCV )            !Loop over the lattice centering vectors
        ICV = ICV+JCV
        IRX1 = MOD(IRX+ICENV(1,ICV),12)
        IRY1 = MOD(IRY+ICENV(2,ICV),12)
        IRZ1 = MOD(IRZ+ICENV(3,ICV),12)

        IF ( RT(5,1,N)+RT(5,1,M2).EQ.0 ) THEN                  !Does this pair generate 1bar?
          M2Z = 1                                  ! No
        ELSE
          M2Z = M2                                    !Yes, they do generate 1Bar
        END IF
        IF ( RT(3,3,N)+RT(3,3,M2Z).LE.0 ) THEN            !Is Z constrained in the unit cell
          IRZ1 = 0                                    ! Yes
        END IF
        IF ( LAUENO.LE.3 .OR. M.NE.4 ) THEN                  ! Does this operator operate along the face diagonal?
          IF ( RT(1,1,N)+RT(1,1,M2Z).LE.0 ) IRX1=0            ! No
          IF ( RT(2,2,N)+RT(2,2,M2Z).LE.0 ) IRY1=0
        ELSE
          IRX1 = MOD(IRX1+IRY1,12)                        ! Yes
          IRY1 = 0
        END IF
        TOTTR = 144*IRX1+12*IRY1+IRZ1
      END DO
      IF ( TOTTR.NE.0 ) THEN
        IER = 18
!        IF ( LPT.GT.0 ) THEN
!          WRITE (LPT,2991) M,N,M2
!2991      FORMAT (' Operator ',I2,' generates Matrix',I3,
!     1      ' which has a translation conflict with',
!     1      ' matrix ',I2)
!          WRITE (LPT,'(A,I3,A,3(I4,2I3),3F5.2,F8.1)') '  Matrix',N,
!     1      ' is',((NINT(RT(I,J,N)),J=1,3),I=1,3),(RT(I,4,N),I=1,3)
!     1      ,RT(5,2,N)
!          WRITE (LPT,'(A,I3,A,3(I4,2I3),3F5.2,F8.1)') '  Matrix',M2,
!     1      ' is',((NINT(RT(I,J,M2)),J=1,3),I=1,3),(RT(I,4,M2),I=1,3)
!     1      ,RT(5,2,M2)
!        END IF
      END IF
      RETURN
      END
