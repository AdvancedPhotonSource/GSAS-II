      SUBROUTINE PACK_F(N,CMPR,MX,MY,IMG)

Cf2py intent(in) N
Cf2py intent(in) CMPR
Cf2py depend(N) CMPR
Cf2py intent(in) MX
Cf2py intent(in) MY
Cf2py intent(in,out) IMG
Cf2py depend(MX,MY) IMG

      IMPLICIT NONE
      INTEGER*4 BITDECODE(0:7),SETBITS(0:16),IN,N,MX,MY,BITNUM
      INTEGER*4 PIXEL,SPILLBITS,USEDBITS,VALIDS,WINDOW,TOTAL
      INTEGER*4 IMG(0:MX-1,0:MY-1),NEXTINT
      INTEGER*4 SPILL,ROW,COL,PIXNUM,MM1
      INTEGER*2 TMP
      CHARACTER*1 CMPR(0:N-1)
      DATA BITDECODE /0,4,5,6,7,8,16,32/
      DATA SETBITS /Z'0000',Z'0001',Z'0003',Z'0007',
     1  Z'000F',Z'001F',Z'003F',Z'007F',Z'00FF',
     1  Z'01FF',Z'03FF',Z'07FF',Z'0FFF',Z'1FFF',
     1  Z'3FFF',Z'7FFF',Z'FFFF'/

      PIXEL = 0
      SPILLBITS = 0
      SPILL = 0
      USEDBITS = 0
      VALIDS = 0
      WINDOW = 0
      ROW = 0
      COL = 0
      TOTAL = MX*MY
      MM1 = MX-1
      IN = 0
      DO WHILE (PIXEL .LT. TOTAL)
        IF (VALIDS .LT. 6) THEN
          IF (SPILLBITS .GT. 0) THEN
            WINDOW = IOR(WINDOW,ISHFT(SPILL,VALIDS))
            VALIDS = VALIDS + SPILLBITS
            SPILLBITS = 0
          ELSE
            SPILL = ICHAR(CMPR(IN))
            IN = IN+1
            SPILLBITS = 8
          END IF
        ELSE
          PIXNUM = ISHFT(1,IAND(WINDOW,SETBITS(3)))
          WINDOW = ISHFT(WINDOW,-3)
          BITNUM = BITDECODE(IAND(WINDOW,SETBITS(3)))
          WINDOW = ISHFT(WINDOW,-3)
          VALIDS = VALIDS-6
          DO WHILE ( (PIXNUM .GT. 0) .AND. (PIXEL .LT. TOTAL) )
            IF ( VALIDS .LT. BITNUM ) THEN
              IF ( SPILLBITS .GT. 0 ) THEN
                WINDOW = IOR(WINDOW,ISHFT(SPILL,VALIDS))
                IF ( (32-VALIDS) .GT. SPILLBITS ) THEN
                  VALIDS = VALIDS + SPILLBITS
                  SPILLBITS = 0
                ELSE
                  USEDBITS = 32-VALIDS
                  SPILL = ISHFT(SPILL,-USEDBITS)
                  SPILLBITS = SPILLBITS-USEDBITS
                  VALIDS = 32
                END IF
              ELSE
                SPILL = ICHAR(CMPR(IN))
                IN = IN+1
                SPILLBITS = 8
              END IF                
            ELSE
              PIXNUM = PIXNUM-1
              IF ( BITNUM .EQ. 0 ) THEN
                NEXTINT = 0
              ELSE
                NEXTINT = IAND(WINDOW,SETBITS(BITNUM))
                VALIDS = VALIDS-BITNUM
                WINDOW = ISHFT(WINDOW,-BITNUM)
                IF ( BTEST(NEXTINT,BITNUM-1) ) 
     1            NEXTINT = IOR(NEXTINT,NOT(SETBITS(BITNUM)))
              END IF

              ROW = PIXEL/MX
              COL = MOD(PIXEL,MX)
              IF ( PIXEL .GT. MX ) THEN
                IF ( COL .EQ. 0 ) THEN
                  TMP = NEXTINT +
     1              (IMG(MM1,ROW-1)+IMG(COL+1,ROW-1)+
     1              IMG(COL,ROW-1)+IMG(MM1,ROW-2) +2)/4
                ELSE IF ( COL.EQ.MM1 ) THEN
                  TMP = NEXTINT +
     1              (IMG(COL-1,ROW)+IMG(0,ROW)+
     1              IMG(MM1,ROW-1)+IMG(MM1-1,ROW-1) +2)/4
                ELSE
                  TMP = NEXTINT + 
     1              (IMG(COL-1,ROW)+IMG(COL+1,ROW-1)+
     1              IMG(COL,ROW-1)+IMG(COL-1,ROW-1) +2)/4
                END IF
              ELSE IF (PIXEL .NE. 0) THEN
                TMP = IMG(COL-1,ROW)+NEXTINT
              ELSE
                TMP = NEXTINT
              END IF
              IMG(COL,ROW) = TMP
              PIXEL = PIXEL+1
            END IF
          END DO
        END IF      
      END DO
      DO ROW=0,MM1
        DO COL=0,MM1
            IF ( IMG(COL,ROW).LT.0 ) IMG(COL,ROW) = IMG(COL,ROW)+65536
        END DO
      END DO
      
      RETURN
      END

      
      SUBROUTINE PACK_F3(N,CMPR,MX,MY,IMG)

Cf2py intent(in) N
Cf2py intent(in) CMPR
Cf2py depend(N) CMPR
Cf2py intent(in) MX
Cf2py intent(in) MY
Cf2py intent(in,out) IMG
Cf2py depend(MX,MY) IMG

      IMPLICIT NONE
      INTEGER*4 BITDECODE(0:7),SETBITS(0:16),IN,N,MX,MY,BITNUM
      INTEGER*4 PIXEL,SPILLBITS,USEDBITS,VALIDS,WINDOW,TOTAL
      INTEGER*4 IMG(0:MX-1,0:MY-1),NEXTINT
      INTEGER*4 SPILL,ROW,COL,PIXNUM,MM1
      INTEGER*2 TMP
      INTEGER*1 CMPR(0:N-1)
      DATA BITDECODE /0,4,5,6,7,8,16,32/
      DATA SETBITS /Z'0000',Z'0001',Z'0003',Z'0007',
     1  Z'000F',Z'001F',Z'003F',Z'007F',Z'00FF',
     1  Z'01FF',Z'03FF',Z'07FF',Z'0FFF',Z'1FFF',
     1  Z'3FFF',Z'7FFF',Z'FFFF'/

      PIXEL = 0
      SPILLBITS = 0
      SPILL = 0
      USEDBITS = 0
      VALIDS = 0
      WINDOW = 0
      ROW = 0
      COL = 0
      TOTAL = MX*MY
      MM1 = MX-1
      IN = 0
      DO WHILE (PIXEL .LT. TOTAL)
        IF (VALIDS .LT. 6) THEN
          IF (SPILLBITS .GT. 0) THEN
            WINDOW = IOR(WINDOW,ISHFT(SPILL,VALIDS))
            VALIDS = VALIDS + SPILLBITS
            SPILLBITS = 0
          ELSE
            SPILL = ICHAR(CHAR(CMPR(IN)))
            IN = IN+1
            SPILLBITS = 8
          END IF
        ELSE
          PIXNUM = ISHFT(1,IAND(WINDOW,SETBITS(3)))
          WINDOW = ISHFT(WINDOW,-3)
          BITNUM = BITDECODE(IAND(WINDOW,SETBITS(3)))
          WINDOW = ISHFT(WINDOW,-3)
          VALIDS = VALIDS-6
          DO WHILE ( (PIXNUM .GT. 0) .AND. (PIXEL .LT. TOTAL) )
            IF ( VALIDS .LT. BITNUM ) THEN
              IF ( SPILLBITS .GT. 0 ) THEN
                WINDOW = IOR(WINDOW,ISHFT(SPILL,VALIDS))
                IF ( (32-VALIDS) .GT. SPILLBITS ) THEN
                  VALIDS = VALIDS + SPILLBITS
                  SPILLBITS = 0
                ELSE
                  USEDBITS = 32-VALIDS
                  SPILL = ISHFT(SPILL,-USEDBITS)
                  SPILLBITS = SPILLBITS-USEDBITS
                  VALIDS = 32
                END IF
              ELSE
                SPILL = ICHAR(CHAR(CMPR(IN)))
                IN = IN+1
                SPILLBITS = 8
              END IF                
            ELSE
              PIXNUM = PIXNUM-1
              IF ( BITNUM .EQ. 0 ) THEN
                NEXTINT = 0
              ELSE
                NEXTINT = IAND(WINDOW,SETBITS(BITNUM))
                VALIDS = VALIDS-BITNUM
                WINDOW = ISHFT(WINDOW,-BITNUM)
                IF ( BTEST(NEXTINT,BITNUM-1) ) 
     1            NEXTINT = IOR(NEXTINT,NOT(SETBITS(BITNUM)))
              END IF

              ROW = PIXEL/MX
              COL = MOD(PIXEL,MX)
              IF ( PIXEL .GT. MX ) THEN
                IF ( COL .EQ. 0 ) THEN
                  TMP = NEXTINT +
     1              (IMG(MM1,ROW-1)+IMG(COL+1,ROW-1)+
     1              IMG(COL,ROW-1)+IMG(MM1,ROW-2) +2)/4
                ELSE IF ( COL.EQ.MM1 ) THEN
                  TMP = NEXTINT +
     1              (IMG(COL-1,ROW)+IMG(0,ROW)+
     1              IMG(MM1,ROW-1)+IMG(MM1-1,ROW-1) +2)/4
                ELSE
                  TMP = NEXTINT + 
     1              (IMG(COL-1,ROW)+IMG(COL+1,ROW-1)+
     1              IMG(COL,ROW-1)+IMG(COL-1,ROW-1) +2)/4
                END IF
              ELSE IF (PIXEL .NE. 0) THEN
                TMP = IMG(COL-1,ROW)+NEXTINT
              ELSE
                TMP = NEXTINT
              END IF
              IMG(COL,ROW) = TMP
              PIXEL = PIXEL+1
            END IF
          END DO
        END IF      
      END DO
      DO ROW=0,MM1
        DO COL=0,MM1
            IF ( IMG(COL,ROW).LT.0 ) IMG(COL,ROW) = IMG(COL,ROW)+65536
        END DO
      END DO
      
      RETURN
      END

