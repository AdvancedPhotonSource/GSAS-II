      SUBROUTINE QMATRI
      COMMON /DELTA2/ IDEL1,IDEL2
      COMMON /DELTA4/ IDERIV,IFINIS,IQMATF
      COMMON /QMAT/ ISQ11(519),ISQ12(519),ISQ13(519),
     $              ISQ22(519),ISQ23(519),ISQ33(519)
C
C
C     SUBROUTINE 'QMATRI' ...
C        UPPER TRIANGULAR MATRIX GENERATION FOR THE RANGE OF
C        DELTAS SUPPLIED
C
C
C
C     --- IDEL IS THE VALUE OF THE DETERMINANT AND IMATT IS THE TOTAL
C         NON-EQUIVALENT MATRICES CONSISTENT WITH A GIVEN IDEL
C     --- CALCULATE ALL THE UNIQUE MATRICES (IQMATF) CONSISTENT WITH
C         WITH A GIVEN DELTA (IDEL) AND STORE THEM IN THE ARRAYS.
C         EACH UNIQUE MATRIX IS IN THE FORM :
C               Q11     Q12     Q13
C               0       Q22     Q23
C               0       0       Q33
C
C     --- A TOTAL OF 519 UNIQUE MATRICES WILL BE GENERATED FOR
C         DELTAS 2 THROUGH 9 (I.E. 7+13+35+31+91+57+155+130 = 519).
C
      DO 300  IDEL = IDEL1,IDEL2
         DO 200 I = 1,IDEL
         DO 200 J = 1,IDEL
         DO 200 K = 1,IDEL
            IPROD = I*J*K
            IF(IPROD.NE.IDEL) GO TO 200
C
C              --- SET UP DIAGONAL ELEMENTS
               IQ11 = I
               IQ22 = J
               IQ33 = K
C
C
C
C
C
C
C              --- SET TOP TRIANGULAR ELEMENTS EQUAL TO ZERO
               IQ12 = 0
               IQ13 = 0
               IQ23 = 0
C
C              --- GENERATE ALL MATRICES CONSISTENT WITH GIVEN DIAGONAL
               DO 100 II = 1,IDEL
               DO 100 JJ = 1,IDEL
               DO 100 KK = 1,IDEL
                  IQ12 = II - 1
                  IQ13 = JJ - 1
                  IQ23 = KK - 1
C
C                 --- THE VALUE OF EACH NON-DIAGONAL ELEMENT IS RESTRICTED
C                     TO BE LESS THAN THE DIAGONAL ELEMENT IN THE SAME
C                     COLUMN.  IN DOING SO, A UNIQUE UPPER TRIANGULAR
C                     MATRIX IS GENERATED FOR A GIVEN DELTA.
                  IF(.NOT.(IQ12.LT.IQ22.AND.IQ13.LT.IQ33.AND.IQ23.LT.
     $               IQ33)) GO TO 100
C
C                    --- MATRIX IS ACCEPTABLE
                     IQMATF = IQMATF + 1
C
C                    --- STOP PROGRAM EXECUTION IF TOO MANY UPPER
C                        TRIANGULAR MATRICES HAVE BEEN GENERATED
C                        (SHOULD NOT OCCUR)
                     IF(IQMATF.GT.519) STOP  '*QMATRI* ... Error in the
     $generation of upper triangular matrices.'
C
                     ISQ11(IQMATF) = IQ11
                     ISQ12(IQMATF) = IQ12
                     ISQ13(IQMATF) = IQ13
                     ISQ22(IQMATF) = IQ22
                     ISQ23(IQMATF) = IQ23
                     ISQ33(IQMATF) = IQ33
  100          CONTINUE
  200    CONTINUE
  300 CONTINUE
      RETURN
      END
