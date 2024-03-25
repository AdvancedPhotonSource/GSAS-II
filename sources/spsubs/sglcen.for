      SUBROUTINE SGLCEN(LCENT,CEN,NCV)

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 

C       This program was developed for
C                    The Division of Chemistry
C                               of
C               The National Research Council of Canada
C                               by
C       Allen C. Larson, 14 Cerrado Loop, Santa Fe, NM 87505, USA

      INTEGER*4     LCENT               !Lattice centering type flag
      REAL*4        CEN(3,4)            !List of lattice centering vectors
      INTEGER*4     NCV                 !Number of lattcie centering vectors

      REAL*4        CENV(3,6)           
      INTEGER*4     NCVT(7)             

      DATA NCVT/1,2,2,2,2,4,3/
      DATA CENV/  0,0.5,0.5,  0.5,0,0.5,  0.5,0.5,0,  0.5,0.5,0.5,
     1  0.3333333,0.6666667,0.6666667,  0.6666667,0.3333333,0.3333333/

      NCV = NCVT(LCENT)
      CEN(1,1) = 0.0
      CEN(2,1) = 0.0
      CEN(3,1) = 0.0
      IF ( NCV.GT.1 ) THEN
        J = LCENT-1
        IF ( LCENT.EQ.6 ) J=1
        IF ( LCENT.EQ.7 ) J=5
        DO I=2,NCV                                          !Copy the lattice centering vectors
          CEN(1,I) = CENV(1,J)
          CEN(2,I) = CENV(2,J)
          CEN(3,I) = CENV(3,J)
          J = J+1
        END DO
      END IF
      RETURN
      END
