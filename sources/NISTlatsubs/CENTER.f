      SUBROUTINE CENTER
      COMMON /CNTR1/ ICNTR1
      COMMON /MATR1/ U1,V1,W1,U2,V2,W2,U3,V3,W3
C
C
C     SUBROUTINE 'CENTER' ...
C        FIND MATRIX TO TRANSFORM THE INPUT CELL TO A PRIMITIVE CELL
C        (I.E. MATRIX TO REMOVE CELL CENTERING)
C
C
C
C     --- (1/2/3/4/5/6/7)    (P OR RR/A/B/C/F/I/RH)
      GO TO (100,200,300,400,500,600,700) ICNTR1
C
  100 CONTINUE
C
C     --- PRIMITIVE LATTICE
      U1 = 1.0
      V1 = 0.0
      W1 = 0.0
      U2 = 0.0
      V2 = 1.0
      W2 = 0.0
      U3 = 0.0
      V3 = 0.0
      W3 = 1.0
      GO TO 800
  200 CONTINUE
C
C     --- A-CENTERED LATTICE
      U1 =  0.0
      V1 =  0.5
      W1 = -0.5
      U2 =  0.0
      V2 =  0.5
      W2 =  0.5
      U3 =  1.0
      V3 =  0.0
      W3 =  0.0
      GO TO 800
  300 CONTINUE
C
C     --- B-CENTERED LATTICE
      U1 = -0.5
      V1 =  0.0
      W1 =  0.5
      U2 =  0.5
      V2 =  0.0
      W2 =  0.5
      U3 =  0.0
      V3 =  1.0
      W3 =  0.0
      GO TO 800
  400 CONTINUE
C
C     --- C-CENTERED LATTICE
      U1 =  0.5
      V1 = -0.5
      W1 =  0.0
      U2 =  0.5
      V2 =  0.5
      W2 =  0.0
      U3 =  0.0
      V3 =  0.0
      W3 =  1.0
      GO TO 800
  500 CONTINUE
C
C     --- F-CENTERED LATTICE
      U1 = 0.5
      V1 = 0.5
      W1 = 0.0
      U2 = 0.0
      V2 = 0.5
      W2 = 0.5
      U3 = 0.5
      V3 = 0.0
      W3 = 0.5
      GO TO 800
  600 CONTINUE
C
C     --- I-CENTERED LATTICE
      U1 =  0.5
      V1 =  0.5
      W1 = -0.5
      U2 = -0.5
      V2 =  0.5
      W2 =  0.5
      U3 =  0.5
      V3 = -0.5
      W3 =  0.5
      GO TO 800
  700 CONTINUE
C
C     --- RHOMBOHEDRAL CELL SET ON HEXAGONAL AXES
      U1 =  1.0/3.0
      V1 = -1.0/3.0
      W1 = -1.0/3.0
      U2 = -2.0/3.0
      V2 = -1.0/3.0
      W2 = -1.0/3.0
      U3 =  1.0/3.0
      V3 =  2.0/3.0
      W3 = -1.0/3.0
  800 CONTINUE
C
C     --- WRITE MATRIX TO TRANSFORM THE INPUT CELL TO
C         A PRIMITIVE CELL OF THE LATTICE
      CALL OUTPT1(10)
      RETURN
      END
