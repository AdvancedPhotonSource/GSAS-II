      FUNCTION SGOPRN(MVAL)

!Purpose: Determine the symmetry operation flags needed for control of
!            spin flip operations

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 

!CALL SEQUENCE PARAMETERS:

      REAL*4        MVAL                !Packed matrix value

!Local variables:

      LOGICAL*4     NOTFOUND            !Loop control flag
      REAL*4        MATVALS(64)         !Packed matrix values
      REAL*4        OPRFLGS(64)         !Operation flags

!Data statements:

      DATA      MATVALS( 1), OPRFLGS( 1) / 193,  0 /      !1              0+ 0+0+0+0+0
      DATA      MATVALS( 2), OPRFLGS( 2) /-193,  7 /      !1bar                0+ 0+0+4+2+1
      DATA      MATVALS( 3), OPRFLGS( 3) / 185,  1 /      !m[001]          0+ 0+0+0+0+1
      DATA      MATVALS( 4), OPRFLGS( 4) /-185,  6 /      !2[001]              0+ 0+0+4+2+0
      DATA      MATVALS( 5), OPRFLGS( 5) / 139,  2 /      !m[010]               0+ 0+0+0+2+0
      DATA      MATVALS( 6), OPRFLGS( 6) /-139,  5 /      !2[010]              0+ 0+0+4+0+1
      DATA      MATVALS( 7), OPRFLGS( 7) /-131,  4 /      !m[100]              0+ 0+0+4+0+0
      DATA      MATVALS( 8), OPRFLGS( 8) / 131,  3 /      !2[100]              0+ 0+0+0+2+1
      DATA      MATVALS( 9), OPRFLGS( 9) /-257, 98 /      !m[110]              0+16+8+0+0+0
      DATA      MATVALS(10), OPRFLGS(10) / 257, 99 /      !2[110]              0+16+0+0+0+1
      DATA      MATVALS(11), OPRFLGS(11) / 265,  8 /      !m[1-10]        0+16+8+0+0+0
      DATA      MATVALS(12), OPRFLGS(12) /-265, 12 /      !2[1-10]        0+16+0+0+0+1
      DATA      MATVALS(13), OPRFLGS(13) /-299,  8 /      !m[101]              0+16+8+0+0+0
      DATA      MATVALS(14), OPRFLGS(14) / 299, 12 /      !2[101]              0+16+0+0+0+1
      DATA      MATVALS(15), OPRFLGS(15) / 353,  8 /      !m[10-1]        0+16+8+0+0+0
      DATA      MATVALS(16), OPRFLGS(16) /-353, 12 /      !2[10-1]        0+16+0+0+0+1
      DATA      MATVALS(17), OPRFLGS(17) / 123,  8 /      !m[011]              0+16+8+0+0+0
      DATA      MATVALS(18), OPRFLGS(18) /-123, 12 /      !2[011]              0+16+0+0+0+1
      DATA      MATVALS(19), OPRFLGS(19) / 201,  8 /      !m[01-1]        0+16+8+0+0+0
      DATA      MATVALS(20), OPRFLGS(20) /-201, 12 /      !2[01-1]        0+16+0+0+0+1
      DATA      MATVALS(21), OPRFLGS(21) /-221, 16 /      !4+[001]       32+16+0+0+0+0
      DATA      MATVALS(22), OPRFLGS(22) / 221, 23 /      !-4+[001]       32+16+0+0+0+0
      DATA      MATVALS(23), OPRFLGS(23) / 229, 16 /      !4-[001]       32+16+0+0+0+0
      DATA      MATVALS(24), OPRFLGS(24) /-229, 23 /      !-4-[001]       32+16+0+0+0+0
      DATA      MATVALS(25), OPRFLGS(25) /-295, 16 /      !4+[010]       32+16+0+0+0+0
      DATA      MATVALS(26), OPRFLGS(26) / 295, 23 /      !-4+[010]       32+16+0+0+0+0
      DATA      MATVALS(27), OPRFLGS(27) / 349, 16 /      !4-[010]       32+16+0+0+0+0
      DATA      MATVALS(28), OPRFLGS(28) /-349, 23 /      !-4-[010]       32+16+0+0+0+0
      DATA      MATVALS(29), OPRFLGS(29) / 129, 16 /      !4+[100]       32+16+0+0+0+0
      DATA      MATVALS(30), OPRFLGS(30) /-129, 23 /      !-4+[100]       32+16+0+0+0+0
      DATA      MATVALS(31), OPRFLGS(31) / 195, 16 /      !4-[100]       32+16+0+0+0+0
      DATA      MATVALS(32), OPRFLGS(32) /-195, 23 /      !-4-[100]       32+16+0+0+0+0
      DATA      MATVALS(33), OPRFLGS(33) / 345,  0 /      !3+[111]        0+ 0+0+0+0+0
      DATA      MATVALS(34), OPRFLGS(34) /-345,  7 /      !-3+[111]        0+ 0+0+4+2+1
      DATA      MATVALS(35), OPRFLGS(35) / 281,  0 /      !3-[111]        0+ 0+0+0+0+0
      DATA      MATVALS(36), OPRFLGS(36) /-281,  7 /      !-3-[111]        0+ 0+0+4+2+1
      DATA      MATVALS(37), OPRFLGS(37) /-309,  0 /      !3+[11-1]        0+ 0+0+0+0+0
      DATA      MATVALS(38), OPRFLGS(38) / 309,  7 /      !-3+[11-1]        0+ 0+0+4+2+1
      DATA      MATVALS(39), OPRFLGS(39) / 205,  0 /      !3-[11-1]        0+ 0+0+0+0+0
      DATA      MATVALS(40), OPRFLGS(40) /-205,  7 /      !-3-[11-1]        0+ 0+0+4+2+1
      DATA      MATVALS(41), OPRFLGS(41) / 303,  0 /      !3+[1-1]        0+ 0+0+0+0+0
      DATA      MATVALS(42), OPRFLGS(42) /-303,  7 /      !-3+[1-11]        0+ 0+0+4+2+1
      DATA      MATVALS(43), OPRFLGS(43) /-277,  0 /      !3-[1-11]        0+ 0+0+0+0+0
      DATA      MATVALS(44), OPRFLGS(44) / 277,  7 /      !-3-[1-11]        0+ 0+0+4+2+1
      DATA      MATVALS(45), OPRFLGS(45) /-339,  0 /      !3+[-111]        0+ 0+0+0+0+0
      DATA      MATVALS(46), OPRFLGS(46) / 339,  7 /      !-3+[-111]        0+ 0+0+4+2+1
      DATA      MATVALS(47), OPRFLGS(47) /-209,  0 /      !3-[-111]        0+ 0+0+0+0+0
      DATA      MATVALS(48), OPRFLGS(48) / 209,  7 /      !-3-[-111]        0+ 0+0+4+2+1
      DATA      MATVALS(49), OPRFLGS(49) /-248,  0 /      !3+[001]        0+ 0+0+0+0+0
      DATA      MATVALS(50), OPRFLGS(50) / 248,  7 /      !-3+[001]        0+ 0+0+4+2+1
      DATA      MATVALS(51), OPRFLGS(51) /  67,  0 /      !3-[001]        0+ 0+0+0+0+0
      DATA      MATVALS(52), OPRFLGS(52) / -67,  7 /      !-3-[001]        0+ 0+0+4+2+1
      DATA      MATVALS(53), OPRFLGS(53) / -59, 16 /      !6+[001]       32+16+0+0+0+0
      DATA      MATVALS(54), OPRFLGS(54) /  59, 23 /      !-6+[001]       32+16+0+4+2+1
      DATA      MATVALS(55), OPRFLGS(55) / 256, 16 /      !6-[001]       32+16+0+0+0+0
      DATA      MATVALS(56), OPRFLGS(56) /-256, 23 /      !-6-[001]       32+16+0+4+2+1
      DATA      MATVALS(57), OPRFLGS(57) / 112,  4 /      !m[100] hex       32+ 0+0+0+0+0
      DATA      MATVALS(58), OPRFLGS(58) /-112,  3 /      !2[100] hex       32+ 0+0+0+0+1
      DATA      MATVALS(59), OPRFLGS(59) / 157,  4 /      !m[010] hex       32+ 0+0+0+0+0
      DATA      MATVALS(60), OPRFLGS(60) /-157,  3 /      !2[010] hex       32+ 0+0+0+0+1
      DATA      MATVALS(61), OPRFLGS(61) /-149,  8 /      !m[210]               0+16+0+0+0+0
      DATA      MATVALS(62), OPRFLGS(62) / 149, 12 /      !2[210]               0+16+0+0+0+1
      DATA      MATVALS(63), OPRFLGS(63) /-104,  8 /      !m[120]               0+16+0+0+0+0
      DATA      MATVALS(64), OPRFLGS(64) / 104, 12 /      !2[120]               0+16+0+0+0+1

!Code:

      SGOPRN = -1                                        !Set negative as an error flag
      NOTFOUND = .TRUE.
      I = 0
      DO WHILE ( I.LT.64 .AND. NOTFOUND )
        I = I+1
        IF ( MVAL.EQ.MATVALS(I) ) THEN                        !Search for the operation in the list above
          NOTFOUND = .FALSE.
          SGOPRN = OPRFLGS(I)
        END IF
      END DO
      RETURN
      END
