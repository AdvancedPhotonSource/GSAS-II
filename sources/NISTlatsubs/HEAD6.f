      SUBROUTINE HEAD6
      COMMON /UNIT2/ IUNITB
C
C
C     SUBROUTINE 'HEAD6' ...
C        WRITE SHORT DESCRIPTION OF PROGRAM OUTPUT
C        SYM FUNCTION (SYMMETRY DETERMINATION THROUGH
C                      CONVERSE TRANSFORMATION ANALYSIS)
C
C
C     --- WRITE HEADING FOR OUTPUT
      WRITE(IUNITB,9020)
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
      WRITE(IUNITB,6290)
      WRITE(IUNITB,6300)
      WRITE(IUNITB,6310)
      WRITE(IUNITB,6320)
      WRITE(IUNITB,6330)
      WRITE(IUNITB,6340)
      WRITE(IUNITB,6350)
      WRITE(IUNITB,6360)
      WRITE(IUNITB,6370)
      WRITE(IUNITB,6380)
      WRITE(IUNITB,6390)
      WRITE(IUNITB,6400)
      WRITE(IUNITB,6410)
      WRITE(IUNITB,6420)
      WRITE(IUNITB,6430)
      WRITE(IUNITB,6440)
      WRITE(IUNITB,6450)
      WRITE(IUNITB,6460)
      WRITE(IUNITB,6470)
      WRITE(IUNITB,6480)
      WRITE(IUNITB,6490)
      WRITE(IUNITB,6500)
      WRITE(IUNITB,6510)
      WRITE(IUNITB,6520)
      WRITE(IUNITB,6530)
      WRITE(IUNITB,6540)
      WRITE(IUNITB,6550)
      WRITE(IUNITB,6560)
      WRITE(IUNITB,9000)
      RETURN
 6100 FORMAT(5X, 4X,'** SYMMETRY DETERMINATION through Converse-Transfor
     $mation Analysis **'//)
 6110 FORMAT(5X, 5X,'In sharp contrast to other methods which focus on t
     $he consequences of')
 6120 FORMAT(5X,'symmetry (such as dot products, d-spacings, etc.), this
     $ approach deals with')
 6130 FORMAT(5X,'symmetry in its most abstract form - represented as mat
     $rices.  The basis of')
 6140 FORMAT(5X,'the SYM program function is to generate a group of matr
     $ices reflecting the')
 6150 FORMAT(5X,'holohedry of the lattice.  This is accomplished through
     $ a specialized')
 6160 FORMAT(5X,'application of the Converse-Transformation operator, de
     $fined as'/)
 6170 FORMAT(5X,18X,'CT ]  (Y,Z)  =  { (H,T)   }'      )
 6180 FORMAT(5X,18X,'   h,t                 i,j  i=0,n')
 6190 FORMAT(5X,18X,'                            j=0,m'/)
 6200 FORMAT(5X,'where h,t represent the domains of H,T, respectively, a
     $nd Y,Z,H,T represent')
 6210 FORMAT(5X,'vector triples in R3 (Karen and Mighell, U.S. Patent Pe
     $nding).'/)
 6220 FORMAT(5X, 5X,'The SYM program function generates the matrices, H,
     $ that relate ANY')
 6230 FORMAT(5X,'primitive cell of the lattice to itself.  In theory, it
     $ is the nature of the')
 6240 FORMAT(5X,'matrices themselves that defines the set to be analyzed
     $ (i.e., those')
 6250 FORMAT(5X,'defining a symmetry group).  In practice, however, the
     $usual result is that')
 6260 FORMAT(5X,'the tolerance matrices alone clearly define the groups
     $and all that is')
 6270 FORMAT(5X,'required to determine metric lattice symmetry and pseud
     $osymmetry is to')
 6280 FORMAT(5X,'count.  The numbers of matrices for the seven lattice m
     $etric symmetries')
 6290 FORMAT(5X,'are: triclinic, 1; monoclinic, 2; orthorhombic, 4; rhom
     $bohedral, 6;')
 6300 FORMAT(5X,'tetragonal, 8; hexagonal, 12; and cubic, 24.'/)
 6310 FORMAT(5X, 5X,'Generated with each symmetry matrix is a tolerance
     $matrix of the form')
 6320 FORMAT(5X,15X,'tol a      tol b      tol c')
 6330 FORMAT(5X,15X,'tol alpha  tol beta   tol gamma  .')
 6340 FORMAT(5X,'Simply by averaging the group of tolerance matrices, an
     $ error matrix is')
 6350 FORMAT(5X,'calculated that may be compared directly to the e.s.d.'
     $'s for the primitive')
 6360 FORMAT(5X,'cell.  More importantly, this error matrix (= averaged
     $tolerance matrix) may')
 6370 FORMAT(5X,'be used to calculate an idealized cell reflecting exact
     $ metric symmetry')
 6380 FORMAT(5X,'(i.e., idealized cell = A1 + avg.tol.a, ... GAMMA1 + av
     $g.tol.gamma).  By')
 6390 FORMAT(5X,'incorporating information from both sets of matrices, H
     $ and T, the idealized')
 6400 FORMAT(5X,'cell provides a unique means of evaluating experimental
     $ error based on ALL')
 6410 FORMAT(5X,'the symmetry operations of the lattice.  This is crucia
     $l when determining')
 6420 FORMAT(5X,'standard or conventional cells.'/)
 6430 FORMAT(5X, 5X,'When evaluating symmetry using CT analysis, the exp
     $erimentalist need')
 6440 FORMAT(5X,'not rely solely on metric information.  The group of H
     $matrices may be')
 6450 FORMAT(5X,'viewed as sets of equivalent (h,k,l)''s represented in
     $matrix form.  Write')
 6460 FORMAT(5X,'the indices of a known reflection as a column matrix an
     $d premultiply by H to')
 6470 FORMAT(5X,'generate sets of reflections that should have equivalen
     $t intensities if the')
 6480 FORMAT(5X,'metric and crystal symmetry agree.'/)
 6490 FORMAT(5X, 5X,'Once the symmetry and error have been evaluated, ob
     $taining a')
 6500 FORMAT(5X,'transformation matrix to a standard or conventional cel
     $l may be accomplished')
 6510 FORMAT(5X,'either through analysis of the H matrices themselves, o
     $r by inputing the')
 6520 FORMAT(5X,'idealized cell reflecting exact metric symmetry into th
     $e RSS program')
 6530 FORMAT(5X,'function with subsequent Table look-up (*Note* This is
     $possible ONLY after')
 6540 FORMAT(5X,'the symmetry and experimental error have been evaluated
     $ through CT analysis).'/)
 6550 FORMAT(5X, 5X,'For additional details as well as discussions of ot
     $her theoretical and')
 6560 FORMAT(5X,'practical applications in crystallography, see referenc
     $e cited below.')
 9000 FORMAT(1H1)
 9020 FORMAT(//)
      END
