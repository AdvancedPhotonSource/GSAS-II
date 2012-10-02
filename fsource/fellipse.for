      SUBROUTINE FELLIPSE(NMAX,X,Y,BB,TAU,RNORM)

Cf2py intent(in) NMAX
Cf2py intent(in) X,Y
Cf2py intent(in) Y
Cf2py depend(NMAX) X,Y
Cf2py intent(in)  TAU
Cf2py intent(out) BB
Cf2py intent(out) RNORM

C
C     NMAX is the number of data points
C     A, B : work arrays
C     X, Y : data points coordinates 
C
C     Daniel Pfenniger, Geneva Observatory, 6/1991
C
C
      REAL*8    X(NMAX),Y(NMAX)
      REAL*4    A(5,NMAX),B(NMAX),BB(0:4),TAU
      INTEGER*4 KRANK,NMAX
      REAL*4    RNORM
C
C
C     
      CALL FITQDR(NMAX,NMAX,A,B,X,Y,TAU,KRANK,RNORM)

C     Coefficients of the quadratic form 
C
      BB(0) = B(1)
      BB(1) = B(2)
      BB(2) = B(3)
      BB(3) = B(4)
      BB(4) = B(5)
C
      RETURN
      END
C-----------------------------------------------------------------------
C
      SUBROUTINE FITQDR(N,NMAX,A,B,X,Y,TAU,KRANK,RNORM)
C
C     Fortran subroutine using the linear least squares routine HFTI 
C     finding the parameters  a, c, d, e, f  defining the best 
C     quadratic form 
C
C            (1+a) X^2 + (1-a) Y^2 + c X Y + d X + e Y + f = 0
C
C     passing through a set of data points {Xi,Yi} i=1,..N 
C     The result is invariant by rotation and translation.
C     If N<6, this program finds an exact solution passing through the
C             points
C     If N>5, it finds the quadratic form which minimizes
C
C            SUM(i=1,N)  [ (1+a)X^2 + (1-a)Y^2 + c X Y + d X + e Y + f ]^2
C
C     N     : Actual number of data points
C     NMAX  : Maximum number of data points
C     NPAR  : Number of parameters to find (5)
C     A, B  : Work arrays
C     X, Y  : Data point coordinates
C     KRANK : Rank of the data matrix A (should be 5)
C     RNORM : Norm of the residual
C
C     TAU   : tolerance for a zero in HFTI
C
C     On outpout  a,c,d,e,f  are in B(1),B(2), ... B(5) respectively
C
C     Daniel Pfenniger, Geneva Observatory, 6/1991
C
      PARAMETER (NPAR=5)
C
      INTEGER*4 NMAX,KRANK
      REAL*4 A(NMAX,NPAR), B(NMAX),TAU, RNORM
      REAL*8 X(NMAX), Y(NMAX)
C
C     H, G, IP : Work arrays for HFTI
C
      INTEGER*4 H(NPAR), IP(NPAR),I,N
      REAL*4    G(NPAR)
C
C     Building the matrices A and B 
C
      DO 10 I = 1, N
        A(I,1) = X(I)**2-Y(I)**2
        A(I,2) = X(I)*Y(I)
        A(I,3) = X(I)
        A(I,4) = Y(I)
        A(I,5) = 1. 
   10   B(I) = -(X(I)**2+Y(I)**2)
C
C     Solving the LS problem: ||A Z - B|| minimum
C     Call HFTI, the result Z is in the first NPAR rows of B
C     TAU is the tolerance for zero
C     KRANK is the rank of A (should be NPAR)
C     RNORM is the norm of the residual
C     See Lawson & Hanson (1974) for detail
      CALL HFTI(A,NMAX,N,NPAR,B,NMAX,1,TAU,KRANK,RNORM,H,G,IP)
C
      RETURN
      END

C     Some LSQ routines from Lawson & Hanson 
C
      SUBROUTINE HFTI (A,MDA,M,N,B,MDB,NB,TAU,KRANK,RNORM,H,G,IP)
C     C.L.LAWSON & R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUN 12
C     TO APPEAR IN 'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL 1974
C     SOLVE LEAST SQUARES PROBLEM USING ALGORITHM HFTI.
C
      INTEGER*4 MDA,M,N,MDB,NB,KRANK,H
      REAL*4 A,B,TAU,RNORM,G
      DIMENSION A(MDA,N),B(MDB,NB),H(N),G(N),RNORM(NB)
      INTEGER IP(N)
      DOUBLE PRECISION SM,DZERO
      PARAMETER (SZERO=0., DZERO=0.D0, FACTOR = 0.001 )
C
      K = 0
      LDIAG = MIN0(M,N)
      IF (LDIAG.LE.0) GOTO 270
        DO 80 J=1,LDIAG
          IF (J.EQ.1) GOTO 20
C
C     UPDATE SQUARED COLUMN LENGTHS AND FIND LMAX
C
          LMAX = J
            DO 10 L = J,N
              H(L) = H(L) - A(J-1,L)**2
              IF (H(L).GT.H(LMAX)) LMAX = L
   10       CONTINUE
          IF (DIFF(HMAX+FACTOR*H(LMAX),HMAX)) 20,20,50
C
C     COMPUTE SQUARED COLUMN LENGTHS AND FIND LMAX
C
   20     LMAX = J
            DO 40 L = J,N
              H(L) = SZERO
              DO 30 I = J,M
   30           H(L) = H(L) + A(I,L)**2
              IF (H(L).GT.H(LMAX)) LMAX = L
   40       CONTINUE
          HMAX = H(LMAX)
C
C     LMAX HAS BEEN DETERMINED
C
C     DO COLUMN INTERCHANGES IF NEEDED.
C
   50     CONTINUE
          IP(J) = LMAX
          IF (IP(J).EQ.J) GOTO 70
          DO 60 I = 1,M
            TMP = A(I,J)
            A(I,J) = A(I,LMAX)
   60       A(I,LMAX) = TMP
          H(LMAX) = H(J)
C
C     COMPUTE THE J-TH TRANSFORMATION AND APPLY IT TO A AND B.
C
   70     CALL H12 (1,J,J+1,M,A(1,J),1,H(J),A(1,J+1),1,MDA,N-J)
   80     CALL H12 (2,J,J+1,M,A(1,J),1,H(J),B,1,MDB,NB)
C
C     DETERMINE THE PSEUDORANK, K, USING THE TOLERANCE, TAU.
C
          DO 90 J = 1,LDIAG
            IF (ABS(A(J,J)).LE.TAU) GOTO 100
   90     CONTINUE
      K = LDIAG
      GOTO 110
  100 K = J - 1
  110 KP1 = K + 1
C
C     COMPUTE THE NORMS OF THE RESIDUAL VECTORS.
C
      IF (NB.LE.0) GOTO 140
        DO 130 JB = 1,NB
          TMP = SZERO
          IF (KP1.GT.M) GOTO 130
          DO 120 I = KP1,M
  120       TMP = TMP + B(I,JB)**2
  130     RNORM(JB) = SQRT(TMP)
  140 CONTINUE
C     SPECIAL FOR PSEUDORANK = 0
      IF (K.GT.0) GOTO 160
      IF (NB.LE.0) GOTO 270
        DO 150 JB = 1,NB
          DO 150 I = 1,N
  150       B(I,JB) = SZERO
      GOTO 270
C
C     IF THE PSEUDORANK IS LESS THAN N COMPUTE HOUSHOLDER
C     DECOMPOSITION OF FIRST K ROWS.
C
  160 IF (K.EQ.N) GOTO 180
        DO 170 II = 1,K
          I = KP1 - II
  170     CALL H12 (1,I,KP1,N,A(I,1),MDA,G(I),A,MDA,1,I-1)
  180 CONTINUE
C
C
      IF (NB.LE.0) GOTO 270
        DO 260 JB = 1,NB
C
C     SOLVE THE K BY K TRIANGULAR SYSTEM
C
          DO 210 L = 1,K
            SM = DZERO
            I = KP1 - L
            IF (I.EQ.K) GOTO 200
            IP1 = I + 1
            DO 190 J = IP1,K
  190         SM = SM + A(I,J)*DBLE(B(J,JB))
  200       SM1 = SM
  210       B(I,JB) = (B(I,JB)-SM1)/A(I,I)
C
          IF (K.EQ.N) GOTO 240
            DO 220 J = KP1,N
  220         B(J,JB) = SZERO
            DO 230 I = 1,K
  230         CALL H12 (2,I,KP1,N,A(I,1),MDA,G(I),B(1,JB),1,MDB,1)
C
C     RE-ORDER THE SOLUTION VECTOR TO COMPENSATE FOR THE
C     COLUMN INTERCHANGES.
C
C     COMPLETE COMPUTATION OF SOLUTION VECTOR
C
  240     DO 250 JJ = 1,LDIAG
            J = LDIAG + 1 - JJ
            IF (IP(J).EQ.J) GOTO 250
            L = IP(J)
            TMP = B(L,JB)
            B(L,JB) = B(J,JB)
            B(J,JB) = TMP
  250     CONTINUE
  260   CONTINUE
C
C     THE SOLUTION VECTORS, X, ARE NOW
C     IN THE FIRST N ROWS OF THE ARRAY B(.).
C
  270 KRANK = K
      RETURN
      END
C-----------------------------------------------------------------------
C
C     SUBROUTINE H12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)
C     C.L.LAWSON & R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUN 12
C     TO APPEAR IN 'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974
C
C     CONSTRUCTION AND/OR APPLICATION OF A SINGLE
C     HOUSHOLDER TRANSFORMATION:  Q = I + U*(U**T)/B
C
C     MODE   = 1 OR 2  TO SELECT ALGORITHM H1 OR H2
C     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT.
C     L1,M   IF L1 .LE. M  THE TRANSFORMATION WILL BE CONSTRUCTED TO
C            ZERO ELEMENTS INDEXED FROM L1 THROUGH M.  IF L1 .GT. M
C            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
C     U(),IUE,UP    ON ENTRY TO H1 U() CONTAINS THE PIVOT VECTOR.
C                   IUE IS THE STORAGE INCREMENT BETWEEN ELEMENTS.
C                   ON EXIT FROM H1 U() AND UP
C                   CONTAIN QUANTITES DEFINING THE VECTOR U OF THE
C                   HOUSHOLDER TRANSFORMATION.  ON ENTRY TO H2 U()
C                   AND UP SHOULD CONTAIN QUATITIES PREVIOUSLY COMPUTED
C                   BY H1. THESE WILL NOT BE MODIFIED BY H2.
C     C()    ON ENTRY TO H1 OR H2 C() CONTAINS A MATRIX WHICH WILL BE
C            REGARDED AS A SET OF VECTORS TO WHICH THE HOUSHOLDER 
C            TRANSFORMATION IS TO BE APPLIED. ON EXIT C() CONTAINS THE
C            SET OF TRANSFORMED VECTORS.
C     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C().
C     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C().
C     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED. IF NCV .LE. 0
C            NO OPERATIONS WILL BE DONE ON C().
C
      SUBROUTINE H12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)
      INTEGER*4 MODE,LPIVOT,L1,M,IUE,ICE,ICV,NCV
      REAL*4 U,UP,C,X,Y
      DIMENSION U(IUE,M),C(*)
      DOUBLE PRECISION SM,B
      PARAMETER (ONE = 1.)
C
      IF (0.GE.LPIVOT.OR.LPIVOT.GE.L1.OR.L1.GT.M) RETURN
      CL = ABS(U(1,LPIVOT))
      IF (MODE.EQ.2) GOTO 60
C     *** CONSTRUCT THE TRANSFORMATION ***
      DO 10 J = L1,M
   10   CL = AMAX1(ABS(U(1,J)),CL)
      IF (CL) 130,130,20
   20 CLINV = ONE/CL
      SM = (DBLE(U(1,LPIVOT))*CLINV)**2
      DO 30 J = L1,M
   30   SM = SM + (DBLE(U(1,J))*CLINV)**2
C       CONVERT DBLE PREC SM TO SNGL. PREC. SM1
      SM1 = SM
      CL = CL*SQRT(SM1)
      IF (U(1,LPIVOT)) 50,50,40
   40 CL = -CL
   50 UP = U(1,LPIVOT) - CL
      U(1,LPIVOT) = CL
      GOTO 70
C     *** APPLY THE TRANSFORMATION I + U*(U**T)*B TO C ***
C
   60 IF (CL) 130,130,70
   70 IF (NCV.LE.0) RETURN
      B = DBLE(UP)*U(1,LPIVOT)
C       B MUST BE NONPOSITIVE HERE. IF B=0.,RETURN
C
      IF (B) 80,130,130
   80 B = ONE/B
      I2 = 1 - ICV + ICE*(LPIVOT-1)
      INCR = ICE*(L1-LPIVOT)
      DO 120 J = 1,NCV
        I2 = I2 +ICV
        I3 = I2 + INCR
        I4 = I3
        SM = C(I2)*DBLE(UP)
        DO 90 I = L1,M
          SM = SM + C(I3)*DBLE(U(1,I))
   90     I3 = I3 + ICE
          IF (SM) 100,120,100
  100     SM = SM*B
          C(I2) = C(I2) + SM*DBLE(UP)
          DO 110 I = L1,M
            C(I4) = C(I4) + SM*DBLE(U(1,I))
  110       I4 = I4 + ICE
  120 CONTINUE
  130 RETURN
      END
C--------------------------- this routine is not a joke!
C
      FUNCTION DIFF(X,Y)
C     C.L.LAWSON & R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUNE 7
C     TO APPEAR IN 'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974
      REAL X,Y
      DIFF = X-Y
      RETURN
      END
