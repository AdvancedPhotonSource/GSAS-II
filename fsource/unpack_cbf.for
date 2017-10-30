      SUBROUTINE UNPACK_CBF(N,CMPR,MXY,IMG)

Cf2py intent(in) N
Cf2py intent(in) CMPR
Cf2py depend(N) CMPR
Cf2py intent(in) MXY
Cf2py intent(in,out) IMG
Cf2py depend(MXY) IMG

      IMPLICIT NONE
      INTEGER*4 N,MXY
      CHARACTER*1 CMPR(0:N-1)
      INTEGER*4 IMG(0:MXY-1),BASEPIXEL
      INTEGER*4 I,J,ISIZE
      CHARACTER*1 C1,E1
      CHARACTER*2 C2,E2
      CHARACTER*4 C4,E4
      INTEGER*1 IONEBYTE
      INTEGER*2 ITWOBYTES
      INTEGER*4 IFOURBYTES

      E1 = CHAR(128)
      E2 = CHAR(0)//CHAR(128)
      E4 = CHAR(0)//CHAR(0)//CHAR(0)//CHAR(128)

      I = 0
      J = 0
      BASEPIXEL = 0
      DO WHILE ( I.LT.N )
        C1 = CMPR(I)
        ISIZE = 1
        IF ( C1.EQ.E1 ) THEN
           ISIZE = 2
           I = I+1
           C2 = CMPR(I)//CMPR(I+1)
           IF ( C2.EQ.E2 ) THEN
              ISIZE = 4
              I = I+2
              C4 = CMPR(I)//CMPR(I+1)//CMPR(I+2)//CMPR(I+3)
              IF ( C4.EQ.E4 ) THEN
                 ISIZE = 8
                 I = I+4
              END IF
           END IF
        END IF
        IF ( ISIZE .EQ. 1 ) THEN
           IONEBYTE = ICHAR(CMPR(I))
           I = I+1
           BASEPIXEL = BASEPIXEL+IONEBYTE
        ELSE IF ( ISIZE .EQ. 2 ) THEN
           ITWOBYTES = ICHAR(CMPR(I))
           ITWOBYTES = ITWOBYTES+ISHFT(ICHAR(CMPR(I+1)),8)
           I = I+2
           BASEPIXEL = BASEPIXEL+ITWOBYTES
        ELSE IF ( ISIZE.EQ.4 ) THEN
           IFOURBYTES = ICHAR(CMPR(I))
           IFOURBYTES = IFOURBYTES+ISHFT(ICHAR(CMPR(I+1)),8)
           IFOURBYTES = IFOURBYTES+ISHFT(ICHAR(CMPR(I+2)),16)
           IFOURBYTES = IFOURBYTES+ISHFT(ICHAR(CMPR(I+3)),24)
           I = I+4
           BASEPIXEL = BASEPIXEL+IFOURBYTES
        END IF
c        IF ( MOD(J,100000).EQ.0 ) PRINT *,I,J,BASEPIXEL
        IMG(J) = BASEPIXEL
        J = J+1
      END DO
      RETURN
      END

      SUBROUTINE UNPACK_CBF3(N,CMPR,MXY,IMG)

Cf2py intent(in) N
Cf2py intent(in) CMPR
Cf2py depend(N) CMPR
Cf2py intent(in) MXY
Cf2py intent(in,out) IMG
Cf2py depend(MXY) IMG

      IMPLICIT NONE
      INTEGER*4 N,MXY
      INTEGER*1 CMPR(0:N-1)
      INTEGER*4 IMG(0:MXY-1),BASEPIXEL
      INTEGER*4 I,J,ISIZE
      CHARACTER*1 C1,E1
      CHARACTER*2 C2,E2
      CHARACTER*4 C4,E4
      INTEGER*1 IONEBYTE
      INTEGER*2 ITWOBYTES
      INTEGER*4 IFOURBYTES

      E1 = CHAR(128)
      E2 = CHAR(0)//CHAR(128)
      E4 = CHAR(0)//CHAR(0)//CHAR(0)//CHAR(128)

      I = 0
      J = 0
      BASEPIXEL = 0
      DO WHILE ( I.LT.N )
        C1 = CHAR(CMPR(I))
        ISIZE = 1
        IF ( C1.EQ.E1 ) THEN
           ISIZE = 2
           I = I+1
           C2 = CHAR(CMPR(I))//CHAR(CMPR(I+1))
           IF ( C2.EQ.E2 ) THEN
              ISIZE = 4
              I = I+2
              C4 = CHAR(CMPR(I))//CHAR(CMPR(I+1))//
     1            CHAR(CMPR(I+2))//CHAR(CMPR(I+3))
              IF ( C4.EQ.E4 ) THEN
                 ISIZE = 8
                 I = I+4
              END IF
           END IF
        END IF
        IF ( ISIZE .EQ. 1 ) THEN
           IONEBYTE = ICHAR(CHAR(CMPR(I)))
           I = I+1
           BASEPIXEL = BASEPIXEL+IONEBYTE
        ELSE IF ( ISIZE .EQ. 2 ) THEN
           ITWOBYTES = ICHAR(CHAR(CMPR(I)))
           ITWOBYTES = ITWOBYTES+ISHFT(ICHAR(CHAR(CMPR(I+1))),8)
           I = I+2
           BASEPIXEL = BASEPIXEL+ITWOBYTES
        ELSE IF ( ISIZE.EQ.4 ) THEN
           IFOURBYTES = ICHAR(CHAR(CMPR(I)))
           IFOURBYTES = IFOURBYTES+ISHFT(ICHAR(CHAR(CMPR(I+1))),8)
           IFOURBYTES = IFOURBYTES+ISHFT(ICHAR(CHAR(CMPR(I+2))),16)
           IFOURBYTES = IFOURBYTES+ISHFT(ICHAR(CHAR(CMPR(I+3))),24)
           I = I+4
           BASEPIXEL = BASEPIXEL+IFOURBYTES
        END IF
c        IF ( MOD(J,100000).EQ.0 ) PRINT *,I,J,BASEPIXEL
        IMG(J) = BASEPIXEL
        J = J+1
      END DO
      RETURN
      END
