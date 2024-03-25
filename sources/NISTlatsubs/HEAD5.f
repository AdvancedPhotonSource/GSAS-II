      SUBROUTINE HEAD5
      COMMON /UNIT2/ IUNITB
C
C
C     SUBROUTINE 'HEAD5' ...
C        WRITE SHORT DESCRIPTION OF PROGRAM OUTPUT
C        REL FUNCTION (RELATE TWO CELLS VIA THE CONVERSE-TRANSFORM)
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
      WRITE(IUNITB,6280)
      WRITE(IUNITB,9000)
      RETURN
 6100 FORMAT(5X,10X,'** RELATE TWO LATTICES via the Converse-Transform *
     $*'//)
 6110 FORMAT(5X, 5X,'This program function makes use of the Converse-Tra
     $nsformation')
 6120 FORMAT(5X,'operator, defined as'/)
 6130 FORMAT(5X,18X,'CT ]  (Y,Z)  =  { (H,T)   }'      )
 6140 FORMAT(5X,18X,'   h,t                 i,j  i=0,n')
 6150 FORMAT(5X,18X,'                            j=0,m'/)
 6160 FORMAT(5X,'where h,t represent the domains of H,T, respectively, a
     $nd Y,Z,H,T represent')
 6170 FORMAT(5X,'vector triples in R3 (Karen and Mighell, U.S. Patent Pe
     $nding).'/)
 6180 FORMAT(5X, 5X,'In this application of the CT operator, the REL pro
     $gram function allows')
 6190 FORMAT(5X,'the user to relate ANY two cells using the specified ma
     $trix elements, h, and')
 6200 FORMAT(5X,'input tolerances, t.  Specifically, all transformation
     $matrices, H, relating')
 6210 FORMAT(5X,'input CELL 2 to input CELL 1 are generated.  The output
     $ tolerance matrix,')
 6220 FORMAT(5X,'associated with a matrix H, represents how closely the
     $transformed CELL 2')
 6230 FORMAT(5X,'agrees with CELL 1, and is of the form'/)
 6240 FORMAT(5X,15X,'tol a      tol b      tol c')
 6250 FORMAT(5X,15X,'tol alpha  tol beta   tol gamma  .'/)
 6260 FORMAT(5X,'Thus, if CELL 2 is defined by A2, B2, C2, ALPHA2, BETA2
     $, GAMMA2, then the')
 6270 FORMAT(5X,'transformation of CELL 2 by a matrix H will give a tran
     $sformed cell having')
 6280 FORMAT(5X,'lattice parameters A1 + tol a, ... GAMMA1 + tol gamma.'
     $  )
 9000 FORMAT(1H1)
      END
