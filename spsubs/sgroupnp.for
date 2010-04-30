      SUBROUTINE SGROUPNP(SPG,LAUENO,NAXIS,NCENT,LCENT,NSYM,NPOL,JRT,
     1  CEN,NCV,RT,IER)

!Purpose:      S.R. which generates a space group from the symbol  - no printing

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 

!       This program was developed for
!                    The Division of Chemistry
!                               of
!               The National Research Council of Canada
!                               by
!       Allen C. Larson, 14 Cerrado Loop, Santa Fe, NM 87505, USA
!
!
!                         This SR interprets the space group symbol
!               Data in the calling sequence are
!       SPG    Input  20 Characters containing the space group symbol
!       LAUENO Output The Laue Group no. where
!                1=1BAR, 2=2/M, 3=MMM, 4=4/M, 5=4/MM, 6=R3R, 7=R3MR,
!                8=3, 9=3M1, 10=31M, 11=6/M, 12=6/MMM, 13=M3 AND 14=M3M
!       NAXIS  Output Unique axis in monoclinic space groups
!               = 4 on error exits; = -1 for rhombahedral in hexagonal setting
!       NCENT  Output 1Bar flag  (0/1) for (acentric/centric)
!       LCENT  Output Lattice centering no.
!                1=P, 2=A, 3=B, 4=C, 5=I, 6=F AND 7=R
!       NSYM   Output The no. of matrices generated
!       NPOL   Output The polar axis flag 
!                1=x, 2=y, 3=x y, 4=z, 5=x z, 6=y z, 7=xyz, 8=111
!       JRT    Output The NSYM (3,5,NSYM) matrices
!       CEN    Output The lattice centering vectors
!       NCV    Output The no. of lattice centering vectors
!       RT     Scratch array of 500 words needed by sgroup
!       IER    Error flag no.

      INTEGER*4     JRT(3,5,24)         !Output matrices, with flags
      CHARACTER*20  SPG                 !Input stribg to be parced
      REAL*4        CEN(3,4)            !Lattice centering vectors
      REAL*4        RT(5,4,25)          !Raw trial matrices with some flags
      REAL*4        D(3,3)              !Origin definition data
      CHARACTER*33  CHR                 !List of characters which will be recognized
      INTEGER*4     LCEN(7)             !Latice centering flags
      INTEGER*4     L(4,4)              !First parcing output, Characters converted to numbers

!               C B A P F I R
      DATA LCEN/4,3,2,1,6,5,7/

!                        111111111122222222223333
!               123456789012345678901234567890123
      DATA CHR/' CBAPFIRMND123456-/H.cbapfirmndh '/
      CHR(33:33) = CHAR(9)                                !Set to "tab"

      IM = 0
      DO I=1,20
        IF ( SPG(I:I).EQ.CHAR(9) ) SPG(I:I) = ' '         !Remove tabs; set to ' '
        IF ( SPG(I:I).NE.' ' ) IM = I
      END DO
      I = 1
      DO WHILE ( I.LE.IM )                                !Squeeze out extra spaces
        IF ( I.LT.20 .AND. SPG(I:I+1).EQ.'  ' ) THEN
          DO J=I+1,IM
            SPG(J:J) = SPG(J+1:J+1)
          END DO
          IM = IM-1
        ELSE
          I = I+1
        END IF
      END DO
      DO I=1,4                                          !Clear the L-array
        DO J=1,4
          L(I,J) = 0
        END DO
      END DO
      K = 1                                              !The number of operator fields
      M = 0                                                !The number of elements in a single field
      IER = 0                                                !General error flag
      NCENT = 0                                          !Set the centric/acentric flag to acentric
      LAUENO = 0                                          !Laue Group number
      NAXIS = 0                                          !Unique axis
      IERX = 0                                          !Error flag of type 2
      N = 0                                                !Matrix count
      J = 1
      DO WHILE ( IER.EQ.0 .AND. J.LE.20 .AND. K.LE.4 )            !Break the space group symbol into the 4 fields
        I = 1                                              !Code as numerical values for manipulation
        DO WHILE ( I.LE.33 .AND. SPG(J:J).NE.CHR(I:I) )            !Search for this character among the legal chars
          I = I+1
        END DO
        IF ( I.LE.33 ) THEN                                    !If character was a legal one
          IF ( I.EQ.32 ) THEN
            I = 20                                          !Convert h to H
          ELSE IF ( I.GT.21 .AND. I.LT.33 ) THEN
            I = I-20                                          !Lower case letters are to be treated as u.c.
          END IF
          IF ( I.GT.1 .AND. I.LT.33 ) THEN                        !We ignore extra spaces
            M = M+1
            L(M,K) = I
            IF ( I.LT.12 .OR. M.GE.4 ) M = 0
            IF ( M.EQ.0 ) K=K+1
          ELSE
            IF ( M.GT.0 ) THEN
              M = 0
              K = K+1
            END IF
          END IF
        ELSE
          IER = 29
        END IF
        J = J+1                                        !Count the input characters
      END DO
      IF ( IER.EQ.0 ) THEN
        K = K-1
        IF ( K.LE.1 ) THEN
          IER = 1                                        !If only 1 field was found.  There is an error.
        END IF

        IF ( IER.EQ.0 ) THEN
          IF ( L(1,1).GT.8 ) THEN
            IER = 2                                          !If the first character was not a P, A, B, C,
          END IF

          IF ( IER.EQ.0 ) THEN
            J = 1
            DO WHILE ( J.LT.4 .AND. IER.EQ.0 )
              J = J+1
              IF ( L(1,J).EQ.18 ) CALL SGLPAK(L(1,J),IER)            !Convert the -N notation to the Nb(ar) notation
            END DO
          END IF
        END IF
      END IF
      IF ( IER.GT.0 ) GO TO 500

      DO I=1,3
        DO J=1,3
          D(I,J) = 0.0                                    !Clear the origin definition translation flags
        END DO
      END DO

      N = 2                                                !Set the matrix count N to 2
      I209 = 0                                          !Clear the body diagonal 3-axis flag
      LCENT = L(1,1)-1                                  !Set the lattice centering flag.
      LCENT = LCEN(LCENT)                                    !   1=P, 2=A, 3=B, 4=C, 5=I, 6=F, 7=R
      IF ( LCENT.NE.7 ) THEN
        CALL SGLATC(K,L,D,LCENT,LAUENO,NAXIS,IER,I209,ID)            !Call a S.R. to determine LAUENO and some
        IF ( IER.GT.0 ) GO TO 500                              !  preliminary data
      ELSE
        IF ( L(1,2).NE.14 ) THEN                              !Rhombohedral lattice.
          IER = 3                                          !Make sure that there is a 3-axis.
          GO TO 500
        ELSE
          IF ( L(1,K).NE.8 ) THEN
            IF ( L(1,K).EQ.20 ) K=K-1                            !Hexagonal axes. R centering. Set LAUENO to 8 or 9
            LAUENO = K+6
          ELSE                                          !Rhombohedral axes.
            LCENT = 1                                        !Delete R centering. Set LAUENO to 6 or 7
            K = K-1
            LAUENO = K+4
            I209 = 1
          END IF
        END IF
      END IF
      CALL SGLCEN(LCENT,CEN,NCV)                              !Establish the list of lattice centering vectors

      IOP = 0                                                !Set the matrix generator flag to 0
      CALL SGRMAT(IOP,RT,1,1.,0.,0.,0.,1.,0.,0.,0.,1.)            !Generate the Idenity operator

      IF ( I209.GT.0 ) THEN
        CALL SGRMAT(IOP,RT,2,0.,0.,1.,1.,0.,0.,0.,1.,0.)            !Cubic or rhombohedral cell. Generate z,x,y
        CALL SGRMAT(IOP,RT,3,0.,1.,0.,0.,0.,1.,1.,0.,0.)            !   and y,z,x
        N = 4
      END IF

      DO MF=2,K                  !Old 3000 loop                  !Decode the last 3 fields of the symbol
        IF ( L(1,MF).EQ.0 ) THEN
          IER = 6
          GO TO 500
        END IF
        IFLD = 1
        DO WHILE ( IFLD.LT.4 .AND. L(IFLD,MF).GT.0 )
          IF ( IFLD.GT.1 ) THEN
            DO WHILE ( IFLD.LE.3 .AND. L(IFLD,MF).NE.19 )
              IF ( L(IFLD,MF).EQ.0 ) THEN
                IFLD = 4
              ELSE
                IF ( L(IFLD,MF).LT.12 ) IER=16
                IF ( IER.GT.0 ) GO TO 500
                IFLD = IFLD+1
              END IF
            END DO
            IFLD = IFLD+1
            IF ( IFLD.LT.5 .AND. L(IFLD,MF).LE.1 ) IER=17
            IF ( IER.GT.0 ) GO TO 500
          END IF
          IF ( IFLD.LT.5 ) THEN
            I = ABS(L(IFLD,MF)-5)
            IF ( I.LE.0 .OR. I.GT.15 ) THEN
              IER = 7
              GO TO 500
            END IF
            NDELT = 1
            NXI = N                                          !Set first matrix pointer
            IF ( I.LE.5 ) THEN                              !Character was A, B, C, M or N
              IF ( MF.EQ.2 .AND. LAUENO.LE.3 ) THEN
                IF ( K.EQ.2 ) THEN                              !Monoclinic B-axis unique
                  IF ( I.EQ.2 ) IER=9
                  IF ( IER.GT.0 ) GO TO 500
                  IOP = 32+2
                  CALL SGRMAT(IOP,RT,N,1.,0.,0.,0.,-1.,0.,0.,0.,1.)      !A B-axis mirror
                  RT(2,4,N) = D(2,2)
                  IF ( I.EQ.1 .OR. I.EQ.5 ) RT(1,4,N) = 0.5
                  IF ( I.EQ.3 .OR. I.EQ.5 ) RT(3,4,N) = 0.5
                ELSE
                  IF ( I.EQ.1 ) IER=8
                  IF ( IER.GT.0 ) GO TO 500
                  IOP = 32+4
                  CALL SGRMAT(IOP,RT,N,-1.,0.,0.,0.,1.,0.,0.,0.,1.)      !An A-axis mirror
                  RT(1,4,N) = D(1,1)
                  IF ( I.EQ.2 .OR. I.EQ.5 ) RT(2,4,N)=0.5
                  IF ( I.EQ.3 .OR. I.EQ.5 ) RT(3,4,N)=0.5
                END IF
              ELSE IF ( MF.EQ.3 .AND. LAUENO.NE.7 ) THEN            !Third field and not a Rombohedral lattice
                IF ( L(1,2).EQ.14 .OR. L(1,2).EQ.17 ) THEN
                  IOP = 32+4
                  CALL SGRMAT(IOP,RT,N,-1.,1.,0.,0.,1.,0.,0.,0.,1.)      !Mirror normal to [100] in hex cell
                  IF ( I.EQ.3 ) RT(3,4,N)=0.5
                ELSE
                  IF ( L(1,2).EQ.15 ) THEN                        !It is not trigonal or hexagonal
                    IF ( I.EQ.1 ) IER=8
                    IF ( IER.GT.0 ) GO TO 500
                    IOP = 32+4
                    CALL SGRMAT(IOP,RT,N,-1.,0.,0.,0.,1.,0.,0.,0.,1.)      !An A-axis mirror
                    RT(1,4,N) = D(1,1)
                    IF ( I.EQ.2 .OR. I.EQ.5 ) RT(2,4,N)=0.5
                    IF ( I.EQ.3 .OR. I.EQ.5 ) RT(3,4,N)=0.5
                  ELSE
                    IF ( I.EQ.2 ) IER=9
                    IF ( IER.GT.0 ) GO TO 500
                    IOP = 32+2
                    CALL SGRMAT(IOP,RT,N,1.,0.,0.,0.,-1.,0.,0.,0.,1.)      !A B-axis mirror
                    RT(2,4,N) = D(2,2)
                    IF ( I.EQ.1 .OR. I.EQ.5 ) RT(1,4,N) = 0.5
                    IF ( I.EQ.3 .OR. I.EQ.5 ) RT(3,4,N) = 0.5
                  END IF
                END IF
              ELSE IF ( MF.EQ.4 .OR. LAUENO.GT.3 ) THEN
                IF ( (MF.EQ.4 .OR. LAUENO.EQ.7) .AND.
     1            (L(1,3).EQ.14 .OR. L(1,2).EQ.15 .OR.
     1            L(1,2).EQ.14 .OR. L(1,2).EQ.17) ) THEN            !It is not cubic or tetragonal
                  IOP = 16+8                                    !Set the op flag to 24
                  CALL SGRMAT(IOP,RT,N,0.,1.,0.,1.,0.,0.,0.,0.,1.)      !A diagonal mirrror normal to [-110]
                  RT(1,4,N) = D(2,2)
                  RT(2,4,N) = -D(2,2)
                  IF ( I.EQ.3 .OR. I.EQ.5 ) RT(3,4,N) = 0.5
                  IF ( (LAUENO.EQ.7 .AND. I.EQ.3) .OR.
     1             (I.LT.3 .OR. I.GT.4) ) THEN
                    IF ( LCENT.EQ.6 .OR. LCENT.EQ.4 ) THEN
                      RT(1,4,N) = 0.25+RT(1,4,N)                  !Either F or C-centered tetragonal.
                      RT(2,4,N) = 0.25+RT(2,4,N)                  !   Glides are 1/4,1/4
                    ELSE
                      RT(1,4,N) = 0.5+RT(1,4,N)
                      RT(2,4,N) = 0.5+RT(2,4,N)
                    END IF
                  END IF
                ELSE
                  IF ( I.EQ.3 ) IER=10
                  IF ( IER.GT.0 ) GO TO 500
                  IF ( LAUENO.GT.12 ) THEN
                    IOP = 32+4
                  ELSE
                    IOP = 1
                  END IF
                  CALL SGRMAT(IOP,RT,N,1.,0.,0.,0.,1.,0.,0.,0.,-1.)      !A C-axis mirror
                  RT(3,4,N) = D(3,3)
                  IF ( I.EQ.1 .OR. I.EQ.5 ) RT(1,4,N) = 0.5
                  IF ( I.EQ.2 .OR. I.EQ.5 ) RT(2,4,N) = 0.5
                  IF ( MF.EQ.2 .AND. L(1,2).EQ.17 .AND. L(2,2).EQ.14 )
     1              RT(3,4,N)=0.5                              !If this a 63-axis the mirror is at 1/4
                END IF
              END IF
            ELSE IF ( I.EQ.6 ) THEN                              !d glide type mirror
              IF ( LCENT.LE.1 ) IER=11
              IF ( IER.GT.0 ) GO TO 500
              ICV = 2
              IF ( MF.EQ.2 .AND. LAUENO.LE.3 ) THEN
                IF ( K.EQ.2 ) THEN
                  IF ( NCV.EQ.4 ) ICV=3
                  IOP = 32+2
                  CALL SGRMAT(IOP,RT,N,1.,0.,0.,0.,-1.,0.,0.,0.,1.)
                  RT(1,4,N) = CEN(1,ICV)/2.0
!                  IF ( LAUENO.EQ.5 ) RT(2,4,N) = D(2,1)
                  RT(3,4,N) = CEN(3,ICV)/2.0
                ELSE
                  IOP = 32+4
                  CALL SGRMAT(IOP,RT,N,-1.,0.,0.,0.,1.,0.,0.,0.,1.)
                  IF ( ID.EQ.2 ) RT(1,4,N)=0.25
                  RT(2,4,N) = CEN(2,ICV)/2.0
                  RT(3,4,N) = CEN(3,ICV)/2.0
                END IF
              ELSE IF ( MF.EQ.3 ) THEN
                IF ( NCV.EQ.4 ) ICV=3
                IOP = 32+2
                CALL SGRMAT(IOP,RT,N,1.,0.,0.,0.,-1.,0.,0.,0.,1.)
                RT(1,4,N) = CEN(1,ICV)/2.0
                IF ( ID.EQ.2 ) RT(2,4,N)=0.25
                IF ( LAUENO.EQ.5 ) RT(2,4,N) = D(2,1)
                RT(3,4,N) = CEN(3,ICV)/2.0
              ELSE IF ( MF.EQ.4 .OR. LAUENO.GT.3 ) THEN
                IF ( MF.EQ.4 .AND. (L(1,2).EQ.15 .OR. L(1,3).EQ.14) )
     1            THEN
                  IOP = 16+8                                    !Set the op flag to 24
                  CALL SGRMAT(IOP,RT,N,0.,1.,0.,1.,0.,0.,0.,0.,1.)      !Cubic or tetragonal. D-glide along diagonal
                  IF ( L(1,3).EQ.13 ) THEN
                    RT(1,4,N) = 0.0
                    RT(2,4,N) = 0.5
                  ELSE
                    RT(1,4,N) = 0.25
                    RT(2,4,N) = 0.25
                  END IF
                  RT(3,4,N) = 0.25
                ELSE
                  IF ( NCV.EQ.4 ) ICV=4
                  IF ( LAUENO.GT.12 ) THEN
                    IOP = 32+4
                  ELSE
                    IOP = 1
                  END IF
                  CALL SGRMAT(IOP,RT,N,1.,0.,0.,0.,1.,0.,0.,0.,-1.)
                  RT(1,4,N) = CEN(1,ICV)/2.0
                  RT(2,4,N) = CEN(2,ICV)/2.0
                  IF ( ID.EQ.2 ) RT(3,4,N)=0.25
                END IF
              END IF
            ELSE IF ( I.EQ.7 ) THEN                              ! 1-fold axis
              NDELT = 0
              IF ( L(2,MF).EQ.18 ) THEN
                NCENT = 1                                  !We have a center of symmetry
                IFLD = IFLD+1
              END IF
            ELSE IF ( I.EQ.8 ) THEN                              !2 fold rotation axis
              IF ( L(2,MF).EQ.18 ) IER=19                        !We will not allow a -2 axis.
              IF ( IER.GT.0 ) GO TO 500
              IF ( MF.EQ.2 ) THEN                              !First rotation operator
                IF ( K.EQ.2 ) THEN
                  IOP = 6
                  CALL SGRMAT(IOP,RT,N,-1.,0.,0.,0.,1.,0.,0.,0.,-1.)      !Rotation about the B-axis
                  RT(1,4,N) = D(1,2)
                  RT(3,4,N) = D(3,2)
                  IF ( L(2,MF).EQ.12 ) RT(2,4,N)=0.5
                ELSE
                  IOP = 32+3
                  CALL SGRMAT(IOP,RT,N,1.,0.,0.,0.,-1.,0.,0.,0.,-1.)      !Rotation about the A-axis.
                  RT(2,4,N) = D(2,1)
                  RT(3,4,N) = D(3,1)
                  IF ( IABS(L(2,MF)-13).EQ.1 ) RT(1,4,N) = 0.5
                END IF
              ELSE IF ( MF.EQ.3 ) THEN                        !Second rotation operator
                IF ( LAUENO.EQ.7 ) THEN
                  IOP = 16+1
                  CALL SGRMAT(IOP,RT,N,0.,-1.,0.,-1.,0.,0.,0.,0.,-1.)      !2-axis along [1-10]
                ELSE IF ( L(1,2).EQ.17 .AND. L(1,4).NE.12 ) THEN
                  IOP = 32+3
                  CALL SGRMAT(IOP,RT,N,1.,-1.,0.,0.,-1.,0.,0.,0.,-1.)      !2-axis along [100] used for the P 6n22 groups
                ELSE IF ( L(1,2).EQ.14 ) THEN
                  IOP = 16+1                                  !op flag will be 9
                  CALL SGRMAT(IOP,RT,N,0.,1.,0.,1.,0.,0.,0.,0.,-1.)      !2-axis along [110] trig
                  RT(1,4,N) = D(2,1)                              ! Also used for the P 3n21 groups
                  IF ( L(2,MF).EQ.12 ) RT(1,4,N)=RT(1,4,N)+0.5
                  RT(2,4,N) = -D(2,1)
                  RT(3,4,N) = D(3,1)
                ELSE                                          !It is not a hexagonal or trigonal space group
                  IOP = 32+5
                  CALL SGRMAT(IOP,RT,N,-1.,0.,0.,0.,1.,0.,0.,0.,-1.)      !Rotation about the B-axis
                  IF ( L(1,2).EQ.9 .AND. L(1,4).EQ.10 ) THEN
                    RT(1,4,N) = 0.5
                  ELSE
                    RT(1,4,N) = D(1,2)
                  END IF
                  RT(3,4,N) = D(3,2)
                  IF ( L(2,MF).EQ.12 ) RT(2,4,N)=0.5
                END IF
              ELSE IF ( MF.EQ.4 ) THEN
                IF ( L(1,2).GE.14 .OR. L(1,3).EQ.14 ) THEN
                  IF ( L(1,2).EQ.15 ) THEN
                    IOP = 32+5                              !op flag should be 37
                    CALL SGRMAT(IOP,RT,N,0.,1.,0.,1.,0.,0.,0.,0.,-1.)      !2-axis along [110] tetrag
                    RT(1,4,N) = D(2,1)
                    IF ( L(2,MF).EQ.12 ) RT(1,4,N)=RT(1,4,N)+0.5
                    RT(2,4,N) = -D(2,1)
                    RT(3,4,N) = D(3,1)
                  ELSE
                    IOP = 16+1
                    CALL SGRMAT(IOP,RT,N,1.,0.,0.,1.,-1.,0.,0.,0.,-1.)!2-axis along [210]
                  END IF
                ELSE
                  IOP = 6
                  CALL SGRMAT(IOP,RT,N,-1.,0.,0.,0.,-1.,0.,0.,0.,1.)      !2-Fold rotation about the C-axis
                  RT(1,4,N) = D(1,3)
                  RT(2,4,N) = D(2,3)
                  IF ( IABS(L(2,MF)-13).EQ.1 ) RT(3,4,N) = 0.5
                  IF ( L(2,MF).EQ.16 ) RT(3,4,N) = 0.5
                END IF
              END IF
            ELSE IF ( I.EQ.9 ) THEN                              !3-fold axis
              IF ( MF.EQ.2 .AND. LAUENO.GT.7 ) THEN
                IOP = 0
                CALL SGRMAT(IOP,RT,N,0.,-1.,0.,1.,-1.,0.,0.,0.,1.)
                IF ( L(2,MF).EQ.12 ) RT(3,4,N)=0.33333333
                IF ( L(2,MF).EQ.13 ) RT(3,4,N)=0.66666667
                IF ( L(2,MF).EQ.18 ) THEN
                  NCENT = 1
                  IFLD = IFLD+1
                 END IF
              ELSE IF ( MF.EQ.3 .OR. LAUENO.LE.7 ) THEN
                NDELT = 0
                IF ( L(2,MF).EQ.18 ) THEN
                  NCENT=1
                  IFLD = IFLD+1
                END IF
              ELSE
                IER = 25
                GO TO 500
              END IF
            ELSE IF ( I.EQ.10 ) THEN
              IF ( MF.NE.2 ) IER=12                              !Four fold axis
              IF ( IER.GT.0 ) GO TO 500
              IF ( L(2,MF).EQ.18 ) THEN
                IOP = 32+16+1
                CALL SGRMAT(IOP,RT,N,0.,1.,0.,-1.,0.,0.,0.,0.,-1.)      !4-bar axis
                RT(1,4,N) = D(1,3)
                RT(2,4,N) = D(2,3)
                RT(3,4,N) = D(3,3)
                IFLD = IFLD+1
              ELSE
                IOP = 32+16
                CALL SGRMAT(IOP,RT,N,0.,-1.,0.,1.,0.,0.,0.,0.,1.)      !4-axis
                RT(1,4,N) = D(1,3)
                RT(2,4,N) = D(2,3)
                IF ( L(2,2).EQ.12 ) RT(3,4,N) = 0.25                  !41 axis
                IF ( L(2,2).EQ.13 ) RT(3,4,N) = 0.5                  !42 axis
                IF ( L(2,2).EQ.14 ) RT(3,4,N) = 0.75                  !43 axis
              END IF
            ELSE IF ( I.EQ.12 ) THEN
              IF ( MF.NE.2 ) IER=13                              !6-axis
              IF ( IER.GT.0 ) GO TO 500
              IF ( L(2,MF).EQ.18 ) THEN
                IOP = 32+16+1
                CALL SGRMAT(IOP,RT,N,-1.,1.,0.,-1.,0.,0.,0.,0.,-1.)      !6-bar operation
                IF ( L(1,3).EQ.2 .OR. L(1,4).EQ.2 ) RT(3,4,N)=0.5
                IFLD = IFLD+1
              ELSE
                IOP = 32+16
                CALL SGRMAT(IOP,RT,N,1.,-1.,0.,1.,0.,0.,0.,0.,1.)      !6 operation
                IF ( L(2,2).GT.11 .AND. L(2,2).LT.17 )
     1            RT(3,4,N)=(L(2,2)-11)/6.0
              END IF
            END IF
            IF ( NDELT.EQ.1 ) THEN
              RT(1,4,N) = MOD(RT(1,4,N)+7.0,1.0)
              RT(2,4,N) = MOD(RT(2,4,N)+7.0,1.0)
              RT(3,4,N) = MOD(RT(3,4,N)+7.0,1.0)
              RT(5,2,N) = 1728*RT(1,4,N)+144*RT(2,4,N)+12*RT(3,4,N)
              RT(5,2,N) = NINT(RT(5,2,N))
              M2 = 1
              IERZ = 0
              DO WHILE ( M2.LT.N .AND. IERZ.EQ.0 )
                IF ( RT(5,1,M2).EQ.RT(5,1,N) ) THEN
                  IERZ = 1                                  !Duplicate rotation matrices
                  IF ( RT(5,2,N).NE.RT(5,2,M2) ) THEN
                    CALL SGTRCF(MF,RT,N,M2,LCENT,LAUENO,IER)      !Different translations
                    IF ( IER.GT.0 ) IERX = IER
                    IER = 0
                  END IF
                ELSE IF ( RT(5,1,M2).EQ.-RT(5,1,N) ) THEN            !New matrix defines a center of symmetry
                  IF ( RT(5,2,N).NE.RT(5,2,M2) ) THEN
                    CALL SGTRCF(MF,RT,N,M2,LCENT,LAUENO,IER)      !Different translations
                    IF ( IER.GT.0 ) IERX = IER
                    IER = 0
                  END IF
                  IERZ = 1
                  NCENT = 1
                END IF
                M2 = M2+1
              END DO
              IF ( IERZ.EQ.0 ) THEN                              !Now if no error has been detected
                N = N+1                                  !Increment the matrix count
                IF ( N.GT.25 ) IER=14
                IF ( IER.GT.0 ) GO TO 500                        !Should never be more than 24
                NXL = N-1                                  !Set NXL to the last currently defined matrix
                DO WHILE ( NXI.LE.NXL )                        !We will repeat this loop until no new matrices
                  DO NX=NXI,NXL
                    DO M1=2,NX
                      CALL SGMTML(RT,NX,M1,N)                        !Apply NX to M1 to generate matrix N
                      IERZ = 0
                      M2 = 1
                      DO WHILE ( M2.LT.N .AND. IERZ.EQ.0 )            !Check for duplication of previous matrix
                        IF ( RT(5,1,N).EQ.RT(5,1,M2) ) THEN
                          IERZ = 1                            !A duplicate
                          IF ( RT(5,2,N).NE.RT(5,2,M2) ) THEN            !Check the translation vectors
                            CALL SGTRCF(MF,RT,N,M2,LCENT,LAUENO,IER)                  !Different translations
                            IF ( IER.GT.0 ) IERX = IER
                            IER = 0
                          END IF
!     PRINT '(a,4i3,a,2i3)','  Duplicate matrix.',NX,M1,N,M2,
!    1    ' Flags are',nint(RT(5,3,N)),nint(RT(5,3,M2))
                        ELSE IF ( RT(5,1,N).EQ.-RT(5,1,M2) ) THEN      !Matrix N is related to M2 by 1bar
                          IERZ = 1
                          NCENT = 1
                        END IF
                        M2 = M2+1
                      END DO
                      IF ( IERZ.EQ.0 ) THEN                        !A new matrix
!     PRINT '(3(a,i3))',' Matrix ',N,' is ',NX,' times ',M1
                        N = N+1                            !Increment the NEW matrix pinter
                        IF ( N.GT.25 ) IER=14
                        IF ( IER.GT.0 ) GO TO 500                  !This pointer should never be larger than 25
                      END IF
                    END DO
                  END DO
                  NXI = NXL+1                                  !Set first matrix to first new matrix
                  NXL = N-1                                  !Set last matrix
                END DO
              END IF
            END IF
          END IF
          IFLD = IFLD+1
        END DO
      END DO                        !end of the old 3000 loop
      NSYM = N-1
      DO K=1,NSYM
        DO I=1,3
          DO J=1,3
            JRT(I,J,K) = RT(I,J,K)
          END DO
          JRT(I,4,K) = 12*RT(I,4,K)+144.1
          JRT(I,4,K) = JRT(I,4,K)-12*(JRT(I,4,K)/12)
          JRT(I,5,K) = RT(5,I,K)
        END DO
        JRT(3,5,K) = SGOPRN(RT(5,1,K))
!        IF ( JRT(3,5,K).LT.0 ) THEN
!          PRINT '(A,I3)',' ***** ERROR in defining operation flags'
!     1      ,K
!        END IF
      END DO
      NPX = 1                                              !Assume X is indeterminate
      NPY = 2                                                !Assume Y is indeterminate
      NPZ = 4                                                !Assume Z is indeterminate
      NPXYZ = 0                                          !Assume no 3-axis along [1,1,1]
      NPYXZ = 1                                        !Assume origin undefined along [1,1,1]
      DO I=1,NSYM                                          !Determine presence of indeterminate origin
        IF ( JRT(1,1,I).LE.0 ) NPX=0                              !Origin is defined along X
        IF ( JRT(2,2,I).LE.0 ) NPY=0                              !Origin is defined along Y
        IF ( JRT(3,3,I).LE.0 ) NPZ=0                              !Origin is defined along Z
        IF ( JRT(1,3,I).GT.0 ) NPXYZ=8                        !There is a 3-axis along [1,1,1]
        IF ( JRT(1,3,I).LT.0 ) NPYXZ=0                        !Origin is defined along [1,1,1]
      END DO
      NPOL = (NPX+NPY+NPZ+NPXYZ*NPYXZ)*(1-NCENT)                  !Set the indeterminate origin flag
!      CALL SGPRNT(SPG,JRT,LAUENO,NAXIS,NCENT,LCENT,NSYM,NPOL,CEN,
!     1  NCV,LPT)
      IF ( LCENT.EQ.7 ) NAXIS = -1
      IF ( IERX.EQ.0 ) RETURN
      IER = IERX
500   CONTINUE
!      IF ( LPTX.GT.0 ) CALL SGERRS(SPG,IER,LPTX)
      NAXIS = 4
      RETURN
      END
