      SUBROUTINE MNCOND
      COMMON /DOTP/ S11,S22,S33,S23,S13,S12
      COMMON /MATR2/ U(9),T(9)
      COMMON /TYPE/ ITYPE
      COMMON /UNIT2/ IUNITB
C**
      COMMON /CK02/ ICK021,ICK022,ICK023,ICK024,ICK025,ICK026,ICK027,
     $              ICK028,ICK029
C*
      DATA PC1/0.00006/
C
C
C     SUBROUTINE 'MNCOND' ...
C        SATISFY MAIN CONDITIONS FOR REDUCTION
C        (THE THREE SHORTEST LATTICE TRANSLATIONS)
C
C
C
C     --- INITIALIZE VARIABLES
      ICT1 = 0
      ICT2 = 0
      ICT3 = 0
C
C
C      --- START CALCULATIONS TO SATISFY MAIN CONDITIONS
   50 CONTINUE
C**
C     --- FOR CHECKING, WRITE EXECUTION POINT AND INTERMEDIATE VARIABLES
      IF(ICK021.EQ.1) CALL CKPT02(1)
C*
C
      A23 = ABS(S23)
      A13 = ABS(S13)
      A12 = ABS(S12)
C
C     --- PROTECTIVE CODE ... THIS CODE SHOULD BE ACTIVE ONLY WHEN
C         THE PROGRAM CONSTANT (PC1) IS INCORRECTLY SET FOR THE
C         SPECIFIC COMPUTER
      ICT1 = ICT1 + 1
      IF(ICT1.GT.21) WRITE(IUNITB,6000)
      IF(ICT1.GT.21) GO TO 300
C
C     --- CONTINUE WHEN HAVE THE SHORTEST VECTORS IN THE BC PLANE;
C         OTHERWISE, CALCULATE
      IF(.NOT.((2.0*A23-S22-PC1).GT.0.0.OR.(2.0*A23-S33-PC1).GT.0.0))
     $   GO TO 100
C**
C        --- FOR CHECKING, WRITE EXECUTION POINT
         IF(ICK022.EQ.1) CALL CKPT02(2)
C*
C
C        --- CALCULATE SHORTEST VECTORS IN BC PLANE,
C            UPDATE THE TOTAL TRANSFORMATION MATRIX AND APPLY
C            THE RESULTING MATRIX TO THE INPUT CELL
         CALL SHORTV(S22,S23,U(1),U(5),U(9),U(8))
         CALL MULTIP
         CALL TRANS(0)
         CALL SHORTV(S33,S23,U(1),U(5),U(9),U(6))
         CALL MULTIP
         CALL TRANS(0)
         GO TO 50
  100 CONTINUE
C
C     --- CONTINUE WHEN HAVE THE SHORTEST VECTORS IN THE AC PLANE;
C         OTHERWISE, CALCULATE
      IF(.NOT.((2.0*A13-S11-PC1).GT.0.0.OR.(2.0*A13-S33-PC1).GT.0.0))
     $   GO TO 200
C**
C        --- FOR CHECKING, WRITE EXECUTION POINT
         IF(ICK022.EQ.1) CALL CKPT02(3)
C*
C
C        --- CALCULATE SHORTEST VECTORS IN AC PLANE,
C            UPDATE THE TOTAL TRANSFORMATION MATRIX AND APPLY
C            THE RESULTING MATRIX TO THE INPUT CELL
         CALL SHORTV(S11,S13,U(1),U(5),U(9),U(7))
         CALL MULTIP
         CALL TRANS(0)
         CALL SHORTV(S33,S13,U(1),U(5),U(9),U(3))
         CALL MULTIP
         CALL TRANS(0)
         GO TO 50
  200 CONTINUE
C
C     --- CONTINUE WHEN HAVE THE SHORTEST VECTORS IN THE AB PLANE;
C         OTHERWISE, CALCULATE
      IF(.NOT.((2.0*A12-S11-PC1).GT.0.0.OR.(2.0*A12-S22-PC1).GT.0.0))
     $   GO TO 300
C**
C        --- FOR CHECKING, WRITE EXECUTION POINT
         IF(ICK022.EQ.1) CALL CKPT02(4)
C*
C
C        --- CALCULATE SHORTEST VECTORS IN AB PLANE,
C            UPDATE THE TOTAL TRANSFORMATION MATRIX AND APPLY
C            THE RESULTING MATRIX TO THE INPUT CELL
         CALL SHORTV(S11,S12,U(1),U(5),U(9),U(4))
         CALL MULTIP
         CALL TRANS(0)
         CALL SHORTV(S22,S12,U(1),U(5),U(9),U(2))
         CALL MULTIP
         CALL TRANS(0)
         GO TO 50
  300 CONTINUE
C
  400 CONTINUE
C
C     --- ARRANGE SO A IS <= TO B AND B IS <= TO C
C
C     --- PROTECTIVE CODE ... THIS CODE SHOULD BE ACTIVE ONLY WHEN
C         THE PROGRAM CONSTANT (PC1) IS INCORRECTLY SET FOR THE
C         SPECIFIC COMPUTER
      ICT2 = ICT2 + 1
      IF(ICT2.GT.11) WRITE(IUNITB,6100)
      IF(ICT2.GT.11) GO TO 500
C
C     --- CHECK THAT A IS <= TO B
      IF(S11.LE.(S22+PC1)) GO TO 420
C
C        --- INTERCHANGE VECTORS A AND B
         CALL SET
         U(2) =  1.0
         U(4) =  1.0
         U(9) = -1.0
         CALL MULTIP
         CALL TRANS(0)
C**
C        --- FOR CHECKING, WRITE EXECUTION POINT AND INTERMEDIATE
C            VARIABLES
         IF(ICK023.EQ.1) CALL CKPT02(5)
C*
C
  420 CONTINUE
C
C     --- CHECK THAT B IS <= TO C
      IF(S22.LE.(S33+PC1)) GO TO 440
C
C        --- INTERCHANGE VECTORS B AND C
         CALL SET
         U(1) =  1.0
         U(6) =  1.0
         U(8) = -1.0
         CALL MULTIP
         CALL TRANS(0)
C**
C        --- FOR CHECKING, WRITE EXECUTION POINT AND INTERMEDIATE
C            VARIABLES
         IF(ICK023.EQ.1) CALL CKPT02(6)
C*
  440 CONTINUE
C
C     --- CHECK THAT A IS <= TO B AND B IS <= TO C
      IF(.NOT.(S11.LE.(S22+PC1).AND.S22.LE.(S33+PC1))) GO TO 400
  500 CONTINUE
C
C     --- PUT CELL MATRIX IN NORMAL REPRESENTATION
      CALL NORMAL
      IF(ITYPE.EQ.1) GO TO 700
C
C        --- PROTECTIVE CODE ... THIS CODE SHOULD BE ACTIVE ONLY
C            WHEN THE PROGRAM CONSTANT (PC1) IS INCORRECTLY SET
C            FOR THE SPECIFIC COMPUTER
         ICT3 = ICT3 + 1
         IF(ICT3.GT.5) WRITE(IUNITB,6200)
         IF(ICT3.GT.5) GO TO  700
C
C           --- TYPE 2 CELL MATRIX ... CHECK PART 2 OF MAIN CONDITIONS
            TX1 = S11 + S22 + PC1
            TX2 = 2.0*(ABS(S23) + ABS(S13) + ABS(S12))
            IF(TX1.GE.TX2) GO TO 700
               CALL SET
               U(1) = 1.0
               U(5) = 1.0
               U(7) = 1.0
               U(8) = 1.0
               U(9) = 1.0
               CALL MULTIP
               CALL TRANS(0)
C**
C              --- FOR CHECKING, WRITE EXECUTION POINT AND
C                  INTERMEDIATE VARIABLES
               IF(ICK024.EQ.1) CALL CKPT02(7)
C*
C
               ICT2 = 0
               GO TO 400
  700 CONTINUE
C
C     --- MAIN CONDITIONS SHOULD BE SATISFIED ...
C         RE-CHECK, WRITE WARNING IF NOT SATISFIED AND CONTINUE
      Q1 = 2.0*ABS(S23) - 2.0*PC1
      Q2 = 2.0*ABS(S13) - 2.0*PC1
      Q3 = 2.0*ABS(S12) - 2.0*PC1
      IF(.NOT.(S11.LE.(S22+2.0*PC1).AND.S22.LE.(S33+2.0*PC1).AND.Q1.LE.
     $   S22.AND.Q2.LE.S11.AND.Q3.LE.S11)) WRITE(IUNITB,6300)
      IF(ITYPE.EQ.1) GO TO 800
C
C        --- TYPE 2 CELL MATRIX (---), TEST EXTRA MAIN CONDITION
         TX1 = S11 + S22 + 2.0*PC1
         TX2 = 2.0*(ABS(S23) + ABS(S13) + ABS(S12))
         IF(TX1.LT.TX2) WRITE(IUNITB,6400)
  800 CONTINUE
      RETURN
 6000 FORMAT(/1X,'*MNCOND* WARNING ... Cell may not be reduced - Check a
     $ll conditions for reduction.'/1X,21X,'Program error when determini
     $ing the shortest vectors in the BC, AC, AB planes.'
     $/1X,21X,'The program constant (PC1) may be incorrectly set for thi
     $s computer.'/)
 6100 FORMAT(/1X,'*MNCOND* WARNING ... Cell may not be reduced - Check a
     $ll conditions for reduction.'/1X,21X,'Program error when ordering
     $the cell edges so that a<=b<=c.'
     $/1X,21X,'The program constant (PC1) may be incorrectly set for thi
     $s computer.'/)
 6200 FORMAT(/1X,'*MNCOND* WARNING ... Cell may not be reduced - Check a
     $ll conditions for reduction.'/1X,21X,'Program error when checking
     $Part II of Main Conditions for a type II cell.'
     $/1X,21X,'The program constant (PC1) may be incorrectly set for thi
     $s computer.'/)
 6300 FORMAT(/1X,'*MNCOND* WARNING ... Cell may not be reduced - Check a
     $ll conditions for reduction.'/1X,21X,'Program error in the final c
     $heck of Part I of the Main Conditions.'
     $/1X,21X,'The program constant (PC1) may be incorrectly set for thi
     $s computer.'/)
 6400 FORMAT(/1X,'*MNCOND* WARNING ... Cell may not be reduced - Check a
     $ll conditions for reduction.'/1X,21X,'Program error in the final c
     $heck of Part II of the Main Conditions.'
     $/1X,21X,'The program constant (PC1) may be incorrectly set for thi
     $s computer.'/)
      END
