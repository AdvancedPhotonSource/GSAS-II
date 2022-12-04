      SUBROUTINE SPCON2(ICOND,EXP1,EXP2,EXP3,ISPFIX,ELE1,ELE2,ELE3,ELE4,
     $                  ELE5,ELE6)
C**
      COMMON /CK02/ ICK021,ICK022,ICK023,ICK024,ICK025,ICK026,ICK027,
     $              ICK028,ICK029
C*
      DATA PC2/0.0005/
C
C
C     SUBROUTINE 'SPCON2' ...
C        CHECK AND, IF NECESSARY, SATISFY SPECIAL CONDITIONS
C        FOR REDUCTION
C
C
C
C     --- WHEN A SPECIAL RELATIONSHIP IN THE CELL MATRIX OCCURS,
C         CHECK, AND IF NECESSARY SATISFY, THE SPECIAL CONDITION
      IF(EXP1.GT.PC2) GO TO 100
C**
C        --- FOR CHECKING, WRITE EXECUTION POINT, SPECIAL CONDITION
C            NUMBER, INTERMEDIATE VARIABLES
         IF(ICK027.EQ.1) CALL CKPT02(30+ICOND)
C*
C
         IF(EXP2.GE.ABS(EXP3)) GO TO 100
C
C           --- FAILED ... NOW SATISFY SPECIAL CONDITION
            ELE1 = -1.0
            ELE2 = -1.0
            ELE3 = -1.0
            ELE4 =  1.0
            ELE5 =  1.0
            ELE6 =  1.0
C
C           --- UPDATE THE TOTAL TRANSFORMATION MATRIX AND
C               APPLY THE RESULTING MATRIX TO THE INPUT CELL
            CALL MULTIP
            CALL TRANS(0)
C**
C           --- FOR CHECKING, WRITE EXECUTION POINT AND INTERMEDIATE
C               VARIABLES
            IF(ICK028.EQ.1) CALL CKPT02(11)
C*
C
            CALL NORMAL
            ISPFIX = 1
  100 CONTINUE
      RETURN
      END