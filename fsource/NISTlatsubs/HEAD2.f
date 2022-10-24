      SUBROUTINE HEAD2
      COMMON /UNIT2/ IUNITB
C
C
C     SUBROUTINE 'HEAD2' ...
C        WRITE SHORT DESCRIPTION OF PROGRAM OUTPUT
C        RSS FUNCTION (REDUCTION AND DERIVATIVE LATTICE)
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
      WRITE(IUNITB,6200)
      WRITE(IUNITB,6210)
      WRITE(IUNITB,6220)
      WRITE(IUNITB,6230)
      WRITE(IUNITB,6240)
      WRITE(IUNITB,6250)
      WRITE(IUNITB,6260)
      WRITE(IUNITB,6270)
      WRITE(IUNITB,9000)
      RETURN
 6100 FORMAT(2X,25X,'** REDUCTION AND DERIVATIVE LATTICE **'/)
 6110 FORMAT(2X,'These calculations fall into two categories:'/)
 6120 FORMAT(2X,' I. Reduction of an input cell.'/)
 6130 FORMAT(2X,'    CELL 1 = Input cell.  This cell may be primitive or
     $ centered (A,B,C,I,F,RR,RH).')
 6140 FORMAT(2X,'    CELL 2 = Reduced primitive cell of the lattice.')
 6150 FORMAT(2X,'    T 1    = A matrix that transforms CELL 1 to a primi
     $tive cell of the lattice.')
 6160 FORMAT(2X,'    T 2    = A matrix that transforms CELL 1 to CELL 2.
     $'//)
 6170 FORMAT(2X,'II. Calculation and reduction of a series of derivative
     $ supercells and/or subcells.')
 6180 FORMAT(2X,'    These derivative cells are calculated from the redu
     $ced cell of the lattice')
 6190 FORMAT(2X,'    (i.e. to carry out the Type II calculation, the pro
     $gram first carries out'/2X,'    the Type I calculation).'/)
 6200 FORMAT(2X,'    CELL 1 = Reduced primitive cell (i.e. CELL 2 from P
     $art I).')
 6210 FORMAT(2X,'    CELL 2 = Reduced supercell or subcell.')
 6220 FORMAT(2X,'    T 1    = A matrix that transforms CELL 1 to a super
     $cell or subcell of the lattice.')
 6230 FORMAT(2X,'    T 2    = A matrix that transforms CELL 1 to CELL 2.
     $'/)
 6240 FORMAT(2X,'For CELL 1 and CELL 2, the output parameters given are:
     $ a, b, c, alpha, beta, gamma')
 6250 FORMAT(2X,'    and volume.  Cell edges are in angstroms and angles
     $ in degrees.'/)
 6260 FORMAT(2X,'The reduced cell matrix is of the form:           a.a
     $ b.b   c.c')
 6270 FORMAT(1X,51X,'b.c   a.c   a.b     .')
 9000 FORMAT(1H1)
      END
