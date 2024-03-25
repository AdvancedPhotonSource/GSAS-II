      SUBROUTINE SPCOND
      COMMON /COSANG/ COSA,COSB,COSG
      COMMON /DOTP/ S11,S22,S33,S23,S13,S12
      COMMON /MATR2/ U(9),T(9)
      COMMON /TYPE/ ITYPE
      COMMON /UNIT2/ IUNITB
      COMMON /VAR1/ VAR90
      DATA PC1/0.00006/
C
C
C     SUBROUTINE 'SPCOND' ...
C        CHECK AND, IF NECESSARY, SATISFY SPECIAL CONDITIONS
C        FOR REDUCTION
C
C
C
C     --- INITIALIZE VARIABLES
      ICOUNT = 0
      DUM1 = 0.0
      DUM2 = 0.0
      DUM3 = 0.0
C
  100 CONTINUE
C
C     --- STARTING POINT TO CHECK, AND IF NECESSARY SATISFY, EACH
C         SPECIAL CONDITION
C
C     --- INITIALIZE VARIABLES
      CALL SET
      ISPFIX = 0
C
C     --- PROTECTIVE CODE ... THIS CODE SHOULD BE ACTIVE ONLY
C         WHEN THE PROGRAM CONSTANT (PC1) IS INCORRECTLY SET
C         FOR THE SPECIFIC COMPUTER
      ICOUNT = ICOUNT + 1
      IF(ICOUNT.LE.20) GO TO 200
         WRITE(IUNITB,6000)
         GO TO  600
  200 CONTINUE
C
C     --- CALCULATE TEMPORARY VARIABLES USED TO CHECK SPECIAL CONDITIONS
      TS1  = ABS(S11 - S22)/S11
      TS2  = ABS(S22 - S33)/S22
      TS3  = ABS(ABS(2.0*S23) - S22)/S22
      TS4  = ABS(ABS(2.0*S13) - S11)/S11
      TS5  = ABS(ABS(2.0*S12) - S11)/S11
      TEM1 = S11 + S22
      TEM2 = 2.0*(ABS(S23) + ABS(S13) + ABS(S12))
      TS6  = ABS(TEM1 - TEM2)/TEM1
      TS7  = 2.0*ABS(S13) + ABS(S12) + PC1
C
C
C     --- SPECIAL CONDITION (A) ...
C         TYPE 1 CELL - IF A.A = B.B THEN  B.C  <=  A.C
C         TYPE 2 CELL - IF A.A = B.B THEN /B.C/ <= /A.C/
      ICOND = 1
      CALL SPCON2(ICOND,TS1,ABS(S13)+PC1,S23,ISPFIX,U(2),U(4),U(9),
     $            DUM1,DUM2,DUM3)
      IF(ISPFIX.EQ.1) GO TO 100
C
C     --- SPECIAL CONDITION (B) ...
C         TYPE 1 CELL - IF B.B = C.C THEN  A.C  <=  A.B
C         TYPE 2 CELL - IF B.B = C.C THEN /A.C/ <= /A.B/
      ICOND = 2
      CALL SPCON2(ICOND,TS2,ABS(S12)+PC1,S13,ISPFIX,U(1),U(6),U(8),
     $            DUM1,DUM2,DUM3)
      IF(ISPFIX.EQ.1) GO TO 100
C
C     --- SPECIAL CONDITION (C) ...
C         TYPE 1 CELL - IF  B.C  = 1/2B.B THEN A.B <= 2A.C
C         TYPE 2 CELL - IF /B.C/ = 1/2B.B THEN A.B = 0
      ICOND = 3
      IF(ITYPE.EQ.1) CALL SPCON2(ICOND,TS3,2.0*S13+PC1,S12,ISPFIX,
     $                           U(1),U(5),U(8),U(9),DUM1,DUM2)
      IF(ITYPE.EQ.2) CALL SPCON2(ICOND,TS3,VAR90,COSG,ISPFIX,
     $                           U(5),U(8),U(9),U(1),DUM1,DUM2)
      IF(ISPFIX.EQ.1) GO TO 100
C
C     --- SPECIAL CONDITION (D) ...
C         TYPE 1 CELL - IF  A.C  = 1/2A.A THEN A.B <= 2B.C
C         TYPE 2 CELL - IF /A.C/ = 1/2A.A THEN A.B = 0
      ICOND = 4
      IF(ITYPE.EQ.1) CALL SPCON2(ICOND,TS4,2.0*S23+PC1,S12,ISPFIX,
     $                           U(1),U(5),U(7),U(9),DUM1,DUM2)
      IF(ITYPE.EQ.2) CALL SPCON2(ICOND,TS4,VAR90,COSG,ISPFIX,
     $                           U(1),U(7),U(9),U(5),DUM1,DUM2)
      IF(ISPFIX.EQ.1) GO TO 100
C
C     --- SPECIAL CONDITION (E) ...
C         TYPE 1 CELL - IF  A.B  = 1/2A.A THEN A.C <= 2B.C
C         TYPE 2 CELL - IF /A.B/ = 1/2A.A THEN A.C = 0
      ICOND = 5
      IF(ITYPE.EQ.1) CALL SPCON2(ICOND,TS5,2.0*S23+PC1,S13,ISPFIX,
     $                           U(1),U(4),U(9),U(5),DUM1,DUM2)
      IF(ITYPE.EQ.2) CALL SPCON2(ICOND,TS5,VAR90,COSB,ISPFIX,
     $                           U(1),U(4),U(5),U(9),DUM1,DUM2)
      IF(ISPFIX.EQ.1) GO TO 100
C
C     --- SPECIAL CONDITIONS (A)-(E) HAVE BEEN SATISFIED.
C         REDUCTION IS COMPLETE FOR A TYPE 1 CELL.
      IF(ITYPE.EQ.1) GO TO 600
C
C        --- SPECIAL CONDITION (F) ... APPLIES ONLY FOR TYPE 2 CELL
C            IF (/B.C/ + /A.C/ + /A.B/) = 1/2(A.A + B.B) ,
C            THEN A.A <= (2/A.C/ + A.B)
         ICOND = 6
         CALL SPCON2(ICOND,TS6,TS7,S11,ISPFIX,
     $              U(1),U(5),DUM1,U(7),U(8),U(9))
         IF(ISPFIX.EQ.1) GO TO 100
C
  600 CONTINUE
      RETURN
 6000 FORMAT(/1X,'*SPCOND* WARNING ... Cell may not be reduced - Check a
     $ll conditions for reduction.'/1X,21X,'Program error when satisfyin
     $g Special Conditions.'
     $/1X,21X,'The program constant (PC1) may be incorrectly set for thi
     $s computer.'/)
      END
