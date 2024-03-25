      SUBROUTINE HMATRX
C
C     * PATENT PENDING *
C     ALL RIGHTS RESERVED.  NO PART OF THIS CODE, OR THE ALGORITHMS
C     UPON WHICH IT IS BASED, MAY BE REPRODUCED, COPIED, TRANSFORMED
C     OR TRANSLATED IN ANY FORM WITHOUT WRITTEN PERMISSION OF THE
C     AUTHORS.
C
      DIMENSION AU(2000),AV(2000),AW(2000)
      DIMENSION BU(2000),BV(2000),BW(2000)
      DIMENSION CU(2000),CV(2000),CW(2000)
      DIMENSION DA(2000),DB(2000),DC(2000)
      DIMENSION TZA(2000),TZB(2000),TZC(2000)
C
      COMMON /HCELL1/ YA,YB,YC,YAL,YBE,YGA,YV
      COMMON /HCELL2/ ZA,ZB,ZC,ZAL,ZBE,ZGA,ZV
      COMMON /HDET1/ HDET
      COMMON /HELEM1/ NHEL
      COMMON /HELEM2/ HEL(99)
      COMMON /HTOL1/ TOLI1,TOLI2,TOLI3,TOLI4,TOLI5,TOLI6
      COMMON /HTOL2/ TOLA,TOLB,TOLC,TOLAL,TOLBE,TOLGA
C
      COMMON /CONST1/ RADIAN
      COMMON /INVER1/ UI1,VI1,WI1,UI2,VI2,WI2,UI3,VI3,WI3
      COMMON /MATR1/ U1,V1,W1,U2,V2,W2,U3,V3,W3
      COMMON /UNIT2/ IUNITB
C**
      COMMON /CK05/ ICK051
C*
C
C
C     SUBROUTINE 'HMATRX' ...
C
C        IMPLEMENTS THE CONVERSE-TRANSFORMATION OPERATOR, DEFINED AS
C
C                 CT ]  (Y,Z) = { (H,T)   }
C                    H,T               I,J  I=0,N
C                                           J=0,M
C
C               ]
C         WHERE H,T REPRESENT THE DOMAINS OF H,T, RESPECTIVELY, AND
C         Y,Z,H,T REPRESENT VECTOR TRIPLES IN R3 (KAREN AND MIGHELL,
C         U.S. PATENT PENDING).
C
C
C
C     COMPUTING NOTE (1987) - THIS ROUTINE COULD BE SHORTENED BY
C        REPLACING CODE WITH CALLS TO SUBROUTINES.  HOWEVER, THESE
C        CALLS *SIGNIFICANTLY* INCREASE THE CPU TIME.
C
      Z11  = ZA*ZA
      Z22  = ZB*ZB
      Z33  = ZC*ZC
      Z23  = ZB*ZC*COS(ZAL/RADIAN)
      Z13  = ZA*ZC*COS(ZBE/RADIAN)
      Z12  = ZA*ZB*COS(ZGA/RADIAN)
      ICTMA = 0
      ICTMB = 0
      ICTMC = 0
      DO 270 JJ = 1,NHEL
         DO 260 KK = 1,NHEL
            DO 250 LL = 1,NHEL
               TZEE =  HEL(JJ)*HEL(JJ)*Z11 + HEL(KK)*HEL(KK)*Z22 +
     $                 HEL(LL)*HEL(LL)*Z33 +
     $            2.0*(HEL(KK)*HEL(LL)*Z23 + HEL(JJ)*HEL(LL)*Z13 +
     $                 HEL(JJ)*HEL(KK)*Z12)
               IF(TZEE.LE.0.0) GO TO 240
                  TZEDG = SQRT(TZEE)
                  IF((TOLI1-ABS(TZEDG-YA)).LT.0.0) GO TO 200
                     ICTMA = ICTMA + 1
                     AU(ICTMA) = HEL(JJ)
                     AV(ICTMA) = HEL(KK)
                     AW(ICTMA) = HEL(LL)
                     DA(ICTMA) = TZEDG - YA
                     TZA(ICTMA) = TZEDG
  200             CONTINUE
                  IF((TOLI2-ABS(TZEDG-YB)).LT.0.0) GO TO 210
                     ICTMB = ICTMB + 1
                     BU(ICTMB) = HEL(JJ)
                     BV(ICTMB) = HEL(KK)
                     BW(ICTMB) = HEL(LL)
                     DB(ICTMB) = TZEDG - YB
                     TZB(ICTMB) = TZEDG
  210             CONTINUE
                  IF((TOLI3-ABS(TZEDG-YC)).LT.0.0) GO TO 220
                     ICTMC = ICTMC + 1
                     CU(ICTMC) = HEL(JJ)
                     CV(ICTMC) = HEL(KK)
                     CW(ICTMC) = HEL(LL)
                     DC(ICTMC) = TZEDG - YC
                     TZC(ICTMC) = TZEDG
  220             CONTINUE
  240          CONTINUE
  250       CONTINUE
  260    CONTINUE
  270 CONTINUE
C
C**
C     --- FOR CHECKING ... PRINT NUMBER OF MATRIX ROWS SATISFYING
C         TRANSFORMATION OF CELL EDGES A,B,C (ONLY), RESPECTIVELY
      IF(ICK051.EQ.1) WRITE(IUNITB,9500) ICTMA,ICTMB,ICTMC
C*
      IF(ICTMA.LE.0.OR.ICTMB.LE.0.OR.ICTMC.LE.0) GO TO 900
      IF(ICTMA.LE.2000.AND.ICTMB.LE.2000.AND.ICTMC.LE.2000) GO TO 300
          WRITE(IUNITB,6000)
          STOP
  300 CONTINUE
      DO 565 J = 1,ICTMA
         DO 555 K = 1,ICTMB
            TZEF = AU(J)*BU(K)*Z11 + AV(J)*BV(K)*Z22 +
     $             AW(J)*BW(K)*Z33 +
     $            (AV(J)*BW(K) + AW(J)*BV(K))*Z23 +
     $            (AU(J)*BW(K) + AW(J)*BU(K))*Z13 +
     $            (AU(J)*BV(K) + AV(J)*BU(K))*Z12
            COS6 = TZEF/(TZA(J)*TZB(K))
            IF(COS6.GE.1.0.OR.COS6.LE.-1.0) GO TO 550
               TZGA = (ACOS(COS6))*RADIAN
               IF((TOLI6-ABS(TZGA-YGA)).LT.0.0) GO TO 550
               DO 545 L = 1,ICTMC
                  HDET = AU(J)*BV(K)*CW(L) + AV(J)*BW(K)*CU(L) +
     $                   AW(J)*BU(K)*CV(L) - CU(L)*BV(K)*AW(J) -
     $                   CV(L)*BW(K)*AU(J) - CW(L)*BU(K)*AV(J)
                  IF(HDET.LE.0.0) GO TO 540
                     TZEF = BU(K)*CU(L)*Z11 + BV(K)*CV(L)*Z22 +
     $                      BW(K)*CW(L)*Z33 +
     $                     (BV(K)*CW(L) + BW(K)*CV(L))*Z23 +
     $                     (BU(K)*CW(L) + BW(K)*CU(L))*Z13 +
     $                     (BU(K)*CV(L) + BV(K)*CU(L))*Z12
                     COS4 = TZEF/(TZB(K)*TZC(L))
                     IF(COS4.GE.1.0.OR.COS4.LE.-1.0) GO TO 540
                        TZAL = (ACOS(COS4))*RADIAN
                        IF((TOLI4-ABS(TZAL-YAL)).LT.0.0) GO TO 540
                        TZEF = AU(J)*CU(L)*Z11 + AV(J)*CV(L)*Z22 +
     $                         AW(J)*CW(L)*Z33 +
     $                        (AV(J)*CW(L) + AW(J)*CV(L))*Z23 +
     $                        (AU(J)*CW(L) + AW(J)*CU(L))*Z13 +
     $                        (AU(J)*CV(L) + AV(J)*CU(L))*Z12
                        COS5 = TZEF/(TZA(J)*TZC(L))
                        IF(COS5.GE.1.0.OR.COS5.LE.-1.0) GO TO 540
                           TZBE = (ACOS(COS5))*RADIAN
                           IF((TOLI5-ABS(TZBE-YBE)).LT.0.0) GO TO 540
                           U1 = AU(J)
                           V1 = AV(J)
                           W1 = AW(J)
                           U2 = BU(K)
                           V2 = BV(K)
                           W2 = BW(K)
                           U3 = CU(L)
                           V3 = CV(L)
                           W3 = CW(L)
                           TOLA = DA(J)
                           TOLB = DB(K)
                           TOLC = DC(L)
                           TOLAL = TZAL - YAL
                           TOLBE = TZBE - YBE
                           TOLGA = TZGA - YGA
                           UI1 =  (V2*W3 - V3*W2)/HDET
                           VI1 = -(V1*W3 - V3*W1)/HDET
                           WI1 =  (V1*W2 - V2*W1)/HDET
                           UI2 = -(U2*W3 - U3*W2)/HDET
                           VI2 =  (U1*W3 - U3*W1)/HDET
                           WI2 = -(U1*W2 - U2*W1)/HDET
                           UI3 =  (U2*V3 - U3*V2)/HDET
                           VI3 = -(U1*V3 - U3*V1)/HDET
                           WI3 =  (U1*V2 - U2*V1)/HDET
                           CALL OUTPT2(7)
  540                      CONTINUE
  545          CONTINUE
  550       CONTINUE
  555    CONTINUE
  565 CONTINUE
C
  900 CONTINUE
      RETURN
 6000 FORMAT(////1X,'*HMATRX* ERROR ... More than 2000 possible matrix r
     $ows have been generated.'/1X,19X,'Input parameters and/or program
     $array sizes must be changed.')
 9500 FORMAT(1X,'ICTMA =',I5,5X,'ICTMB =',I5,5X,'ICTMC =',I5//)
      END
