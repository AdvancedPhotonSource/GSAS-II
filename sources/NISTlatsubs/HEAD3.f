      SUBROUTINE HEAD3
      COMMON /UNIT2/ IUNITB
C
C
C     SUBROUTINE 'HEAD3' ...
C        WRITE SHORT DESCRIPTION OF PROGRAM OUTPUT
C        TRANS FUNCTION (CELL TRANSFORMATION)
C
C
C     --- WRITE HEADING FOR OUTPUT
      WRITE(IUNITB,6100)
      WRITE(IUNITB,6110)
      WRITE(IUNITB,6120)
      WRITE(IUNITB,6130)
      WRITE(IUNITB,6140)
      WRITE(IUNITB,6150)
      WRITE(IUNITB,6160)
      WRITE(IUNITB,6170)
      WRITE(IUNITB,6180)
      WRITE(IUNITB,6190)
      WRITE(IUNITB,9020)
      RETURN
 6100 FORMAT(19X,'** CELL TRANSFORMATION **'/)
 6110 FORMAT(5X,'    CELL 1   =   Input cell.')
 6120 FORMAT(5X,'    CELL 2   =   Transformed cell.')
 6130 FORMAT(5X,'    T 2      =   Input transformation matrix.')
 6140 FORMAT(5X,'    T 2 INV  =   Inverse matrix for T 2.'/)
 6150 FORMAT(5X,'For CELL 1 and CELL 2, the output parameters given are:
     $')
 6160 FORMAT(5X,'    a, b, c, alpha, beta, gamma and volume.  Cell edges
     $')
 6170 FORMAT(5X,'    are in angstroms and angles in degrees.'/)
 6180 FORMAT(5X,'The cell matrix is of the form:      a.a   b.b   c.c')
 6190 FORMAT(34X,'        b.c   a.c   a.b     .')
 9020 FORMAT(//)
      END
