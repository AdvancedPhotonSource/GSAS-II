      SUBROUTINE SGLATC(K,L,D,LCENT,LAUENO,NAXIS,IER,I209,ID)

!Purpose:      Determine Laue group and some other preliminary data

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 

!       This program was developed for
!                    the Division of Chemistry
!                               of
!               The National Research Council of Canada
!                               by
!       Allen C. Larson, 14 Cerrado Loop, SANTA FE, NM 87505, USA

!Calling sequence variables

      INTEGER*4     K                   !Number of fields found in the space group symbol
      INTEGER*4     L(4,4)              !Integer values for the characters in the symbol
      REAL*4        D(3,3)              !Location of some key elements
      INTEGER*4     LCENT               !Lattice centering flag
      INTEGER*4     LAUENO              !Laue Group number
      INTEGER*4     NAXIS               !Unique axis flag for monoclinic cells
      INTEGER*4     IER                 !Error flag
      INTEGER*4     I209                !Diagonal 3-axis flag
      INTEGER*4     ID                  !Number of D-glides

!Local variables:

!Code:

      ID = 0                                !Set no. of d-glides to zero
                                          !Now let us determine the Laue group and unique axis if monoclinic
      IF ( K.EQ.2 ) THEN                                    ! Only 2 fields were read
        IF ( L(1,2).EQ.17 ) THEN
          LAUENO = 11                                        !  6/M
        ELSE IF ( L(1,2).EQ.14 ) THEN
          LAUENO = 8                                          !  3Bar
        ELSE IF ( L(1,2).EQ.15 ) THEN
          LAUENO = 4                                          !  4/M
          IF ( LCENT.GE.5 ) THEN                              !Is it I-centered or F-centered?
            IF ( L(4,2).NE.4 .AND. L(4,2).NE.11 ) GO TO 1553      !Is there either an A-glide or a D-glide normal to C?
            D(1,3) = 0.75                                    !Yes.
            IF ( LCENT.EQ.5 ) D(2,3) = 0.25
          ELSE IF ( LCENT.EQ.4 ) THEN                              !Is it C-centered?   C-centered 4/m tetragonal
            IF ( L(3,2).NE.4 .AND. L(4,2).NE.4 ) GO TO 210            !If there is no A-glide normal to C we are through
            D(1,3) = 0.25
            D(2,3) = 0.25
            IF ( L(4,2).EQ.4 ) D(2,3)=0.75
          ELSE
            IF ( L(3,2).EQ.10 ) THEN                        !No.  Is there a N-glide normal to C?
              D(1,3) = 0.5                                    !  P 4n/n * *
              GO TO 210
            END IF
            IF ( L(4,2).EQ.10 ) D(2,3)=0.5
          END IF
        ELSE IF ( L(1,2).EQ.12 ) THEN
          LAUENO = 1                                        !1Bar
        ELSE IF ( L(1,2).EQ.16 ) THEN
          IER = 5                                            !bad 5-fold
          GO TO 500
        ELSE
          IM = 2                                          !2/M, B-axis unique
          GO TO 1419
        END IF
        GO TO 210
      ELSE IF ( K.EQ.3 ) THEN                  !Only 3 Fields were read.  Must be M3 cubic. (R3r has been taken care of)
        IF ( L(1,3).NE.14 ) THEN
          IER = 20
          GO TO 500
        END IF
        LAUENO = 13
        IF ( L(2,2).EQ.12 ) D(2,1)=0.5                        !Set the B-axis flag if a 21 along A
        IF ( L(1,2).EQ.3 .OR. L(1,2).EQ.4 ) D(3,3)=0.5            !Set the C-axis flag if an A-glide normal to C
        GO TO 209
      ELSE                                                !Four fields were read
        IF ( L(1,3).EQ.14 ) THEN                              !It is m3m cubic
          LAUENO = 14
          IF ( L(1,2).EQ.3 .OR. L(1,2).EQ.4 ) D(3,3)=0.5            !Set the C-axis translation if an A or B normal to C
          IF ( L(1,2).EQ.15 ) THEN                              !a 4n-axis specified
            IF ( L(2,2).EQ.18 ) THEN                              !It is 4bar 3 *
              IF ( L(1,4).NE.9 ) THEN                              !It is not 4bar 3 m
                IF ( L(1,4).EQ.11 ) THEN                        !It is 4bar 3 d
                  IF ( LCENT.NE.5 ) THEN                        ! I 4bar 3 d, we hope
                    IER = 21
                    GO TO 500
                  END IF
                  D(1,3) = 0.75
                  D(2,3) = 0.25
                  D(3,3) = 0.75
                ELSE
                  D(1,3) = 0.5
                  D(2,3) = 0.5
                  D(3,3) = 0.5
                END IF
              END IF
            ELSE IF ( L(2,2).EQ.12 ) THEN                        !41-axis.
              IF ( LCENT.EQ.6 ) THEN                              !  F 41 3 2
                D(1,3) = 0.75
                D(2,3) = 0.75
                D(3,3) = 0.25
              ELSE                                          !IT IS EITHER P 41 3 2 OR I 41 3 2
                D(1,3) = 0.25
                D(2,3) = 0.75
                D(3,3) = 0.25
              END IF
            ELSE IF ( L(2,2).EQ.13 ) THEN                        !  P 42 3 2
              D(1,3) = 0.5
              D(2,3) = 0.5
              D(3,3) = 0.5
            ELSE IF ( L(2,2).EQ.14 ) THEN                        !It is 43 3 2
              D(1,3) = 0.75
              D(2,3) = 0.25
              D(3,3) = 0.75
            END IF
            GO TO 209
          END IF
          GO TO 209
        ELSE IF ( L(1,2).EQ.17 ) THEN                              !It is hexagonal
          IF ( L(1,3).EQ.12 .AND. L(1,4).EQ.12 ) THEN                  !We have something like P 6n 1 *
            LAUENO = 11                                  ! 6/M
          ELSE
            LAUENO = 12                                    ! 6/MMM
          END IF
          GO TO 210
        ELSE IF ( L(1,2).EQ.14 ) THEN                              !It is trigonal
          IF ( L(1,3).EQ.12 ) THEN                              ! P3**
            IF ( L(1,4).EQ.12 ) THEN                              ! 31*
              LAUENO = 8                                    ! 3BAR
            ELSE
              LAUENO = 10                                    ! 31m
            END IF
          ELSE IF ( L(1,4).NE.12 ) THEN
            LAUENO = 12                                    ! 6/MMM
          ELSE
            LAUENO = 9                                    ! 3M1
          END IF
          GO TO 210
        ELSE IF ( L(1,2).EQ.15 ) THEN                              !It is tetragonal 4/MMM
          LAUENO = 5
                                                      !If there is an N-glide normal to C place any
          IF ( L(3,2).EQ.10 .OR. L(4,2).EQ.10 ) D(1,1)=0.5            ! mirror normal to A at 1/4
                                                      !If there is an A-glide normal to C place any
          IF ( L(3,2).EQ.4 .OR. L(4,2).EQ.4 ) D(2,2)=0.25            ! mirror normal to (110) at 1/4
          IF ( L(1,3).EQ.13 .AND. L(2,3).EQ.12 ) D(1,2)=0.5            !If there is a 21 along B move place it at x=1/4
                                                      !If there is a B- or N-glide normal to the A-axis
          IF ( L(1,3).EQ.3 .OR. L(1,3).EQ.10 ) D(1,1)=D(1,1)+0.5      ! shift the mirror by 1/4 along the A-axis
                                                      !If there is either a B- or N-glide normal to (110)
          IF ( L(1,4).EQ.3 .OR. L(1,4).EQ.10 ) D(2,2)=D(2,2)+0.25      ! shift the mirror by 1/4 along the A-axis
          IF ( LCENT.EQ.1 .AND.                              !If Primative
     1      L(2,2).GT.11 .AND. L(2,2).LT.15 .AND.                  ! and this is a 41, 42 or 43
     1      L(2,3).NE.12 )                                    ! and not 4n 21 2
     1      D(3,1)=-(L(2,2)-11)/4.0                              ! Set the Z-location for 2-axes along (110)
          IF ( L(1,4).EQ.13 .AND. L(2,4).EQ.12 .AND.                  !If fourth field is 21
     1      L(2,2).GT.11 .AND. L(2,2).LT.15 )                  ! and this is a 41, 42 or 43
     1      D(3,1)=(L(2,2)-11)/4.0
          IF ( L(1,3).EQ.13 .AND. L(2,3).EQ.12 .AND.                  !Set the Z-translation for 21-axes along B
     1      L(2,2).GT.11 .AND. L(2,2).LT.15 )
     1      D(3,2)=(L(2,2)-11)/4.0
          IF ( L(1,3)+L(3,2).EQ.11 .AND. LCENT.EQ.6 ) D(2,1)=0.75      !Place the D in F 4* D * at Y=7/8
          IF ( L(1,4).EQ.2 .AND. LCENT.EQ.6 ) D(1,1)=0.5            !Set M in F 4** * * at X=1/8 If a C along (110)
          IF ( L(2,2).EQ.18 ) GO TO 1556                        !Is this a 4bar?
          IF ( LCENT.GT.1 ) GO TO 1553                        !Is the lattice primative?
          IF ( L(3,2).EQ.10 .OR. L(4,2).EQ.10 ) GO TO 1552            !Yes.  Do we have a N-glide normal to C?
          IF ( L(1,3).EQ.13 .AND. L(2,3).EQ.12 ) GO TO 1551          !No.  Do we have a 21 along B?
          IF ( L(1,3).NE.10 ) GO TO 210                        !No. Do we have a N-glide normal to A?
          IF ( L(2,2).LE.0 ) GO TO 210
          IF ( L(2,2).GT.15 ) GO TO 210
1551      CONTINUE
          D(1,3) = 0.5
          D(2,3) = 0.5
          GO TO 210
1552      CONTINUE
          D(1,3) = 0.5                                    !  P 4n/n * *
          GO TO 210
1553      CONTINUE
          IF ( LCENT.LT.5 ) GO TO 1555                        !Is the lattice I or F-centered?
                                                !  YES.
          IF ( L(1,4).EQ.2 ) D(2,1)=D(2,1)+0.5                  !If there is a C along (110) place the D at Y=1/4
          IF ( L(4,2).NE.4 .AND. L(4,2).NE.11 ) GO TO 1554            !IS THIS I 41/A * * OR F 41/D * * ?
                                                !  YES.
          D(1,3) = 0.25
          IF ( LCENT.EQ.5 ) D(2,3) = 0.75
          GO TO 210
1554      CONTINUE
          IF ( L(2,2).NE.12 ) GO TO 210                        !Is there a 41 present?
                                                !  YES.
          IF ( LCENT.EQ.6 ) GO TO 1558                        !If F-centered go to 1558
          D(2,3) = 0.5                                        !  SET THE B-AXIS TRANSLATION FLAGS FOR I 41 2 2
          GO TO 1557
1555      CONTINUE
          IF ( LCENT.NE.4 ) IER=23                              !Is the lattice C-centered?
          IF ( IER.GT.0 ) GO TO 500
          IF ( L(3,2).EQ.4 .OR. L(4,2).EQ.4 ) GO TO 1559            !C-Centered.  an A normal to C
!         IF ( L(3,2).EQ.0 ) D(1,1)=2.0*D(2,2)+D(1,1)
          IF ( D(1,1).EQ.0.0 ) D(1,1)=2.0*D(2,2)
          IF ( L(1,4).EQ.13 .AND. L(2,4).EQ.12 ) GO TO 1552            !Is there a 21 on the diagonal?
          IF ( L(2,2).LE.0 ) GO TO 210
          IF ( L(1,4).NE.10 ) GO TO 210                        !Is there a N-glide normal to (110)?
          IF ( L(2,2).GT.15 ) GO TO 210
          D(1,1) = D(1,1)-2.0*D(2,2)
          GO TO 1552
1556      CONTINUE
                                          !  ACCOUNT FOR TRANSLATIONS DUE TO DIAGONAL SYMMETRY OPERATION
          IF ( L(1,3).EQ.11 .AND. LCENT.EQ.6 ) D(3,1)=0.25            !  IF F 4B D 2 WE WANT THE 2 AT Z=1/8
          IF ( L(1,4).EQ.13 .AND. L(2,4).EQ.12 ) D(1,1)=0.5            !  IF * 4B * 21 WE WANT THE MIRROR AT X=1/4
          IF ( L(1,4).EQ.2 .OR. L(1,4).EQ.10 ) D(3,2)=0.5            !If a C- or a N-glide (110) set 2-axis at Z=1/4
          IF ( L(1,4).EQ.3 .OR. L(1,4).EQ.10 ) D(1,2)=0.5            !If a B- or a N-glide (110) SET 2 AT X=1/4
          IF ( L(1,4).NE.11 ) GO TO 210
1557      CONTINUE
          IF ( LCENT.EQ.5 ) D(1,2) = 0.5
          D(3,2) = 0.75
          GO TO 210
1558      CONTINUE
          D(1,3) = 0.25                                    !  F 41 * *
          D(2,3) = 0.75
          GO TO 210
1559      CONTINUE
          D(1,3) = 0.25                                    !  C 4*/A * *
          D(2,3) = 0.25
          IF ( L(1,4).EQ.3 .OR. L(1,4).EQ.10 ) D(1,1)=0.5
          GO TO 210
        ELSE IF ( L(1,2).EQ.12 ) THEN
141       IF ( L(1,3).EQ.12 ) GO TO 143                     !  IT IS NOT C-AXIS UNIQUE MONOCLINIC
          IF ( L(1,4).NE.12 ) GO TO 1399
          IM = 3
1419      CONTINUE
          LAUENO = 2                                          !IT IS B-AXIS UNIQUE MONOCLINIC. (FULL SYMBOL USED)
          NAXIS = 2
          IA = 4
          IC = 2
          NA = 1
          NB = 2
          NC = 3
          GO TO 1430
        ELSE IF ( L(1,3).EQ.12 ) THEN
142       IF ( L(1,4).NE.12 ) GO TO 1399                !  IT IS A-AXIS UNIQUE MONOCLINIC
          LAUENO = 2
          NAXIS = 1
          IA = 3
          IC = 2
          NA = 2
          NB = 1
          NC = 3
          IM = 2
          GO TO 1430
143       IF ( L(1,4).EQ.12 ) THEN
            LAUENO = 1                                  !  1BAR
          ELSE
            LAUENO = 2                                    !  IT IS C-AXIS UNIQUE MONOCLINIC
            NAXIS = 3
            IA = 4
            IC = 3
            NA = 1
            NB = 3
            NC = 2
            IM = 4
1430        CONTINUE
            IF ( L(2,IM).EQ.12 ) D(NB,NAXIS)=0.5
            IF ( L(3,IM).EQ.IA .OR. L(3,IM).EQ.10 ) D(NA,NAXIS)=0.5
            IF ( L(3,IM).EQ.IC .OR. L(3,IM).EQ.10 ) D(NC,NAXIS)=0.5
            IF ( L(4,IM).EQ.IA .OR. L(4,IM).EQ.10 ) D(NA,NAXIS)=0.5
            IF ( L(4,IM).EQ.IC .OR. L(4,IM).EQ.10 ) D(NC,NAXIS)=0.5
          END IF
          GO TO 210
        ELSE
                                                      !It may be orthorhombic
1399      CONTINUE
                                                      !It is orthorhombic
          LAUENO = 3
                                                      !Set up counts of the various types of mirrors.
          IM = 0
          IR = 0
          IA = 0
          IB = 0
          IC = 0
          ID = 0
          I21 = 0
          IF ( L(1,2).NE.13 ) GO TO 1400                        !Do we have a 2-axis along A
          IF ( L(2,2).NE.12 ) GO TO 1401                      !Yes, is it a 21?
          D(1,2) = 0.5
          D(1,3) = 0.5
          I21 = 4
          GO TO 1401
1400      CONTINUE
          IR = 1
          IF ( L(1,2).EQ.9 ) IM=4
          IF ( L(1,2).EQ.3 ) IB=1
          IF ( L(1,2).EQ.2 ) IC=1
          IF ( L(1,2).EQ.11 ) ID=1
          IF ( L(1,3).EQ.4 .OR. L(1,3).EQ.10 ) D(1,1)=0.5
          IF ( L(1,4).EQ.4 .OR. L(1,4).EQ.10 ) D(1,1)=D(1,1)+0.5

1401      CONTINUE
          IF ( L(1,3).NE.13 ) GO TO 1402                        !Do we have a 2-axis along B
          IF ( L(2,3).NE.12 ) GO TO 1403                        !Yes, is it a 21?
          D(2,1) = 0.5                                    !Yes, it is a 21
          D(2,3) = 0.5
          I21 = I21+2
          GO TO 1403
1402      CONTINUE
          IR = IR+1
          IF ( L(1,3).EQ.9 ) IM=IM+2
          IF ( L(1,3).EQ.4 ) IA=1
          IF ( L(1,3).EQ.2 ) IC=IC+1
          IF ( L(1,3).EQ.11 ) ID=ID+1
          IF ( L(1,2).EQ.3 .OR. L(1,2).EQ.10 ) D(2,2)=0.5
          IF ( L(1,4).EQ.3 .OR. L(1,4).EQ.10 ) D(2,2)=D(2,2)+0.5

1403      CONTINUE
          IF ( L(1,4).NE.13 ) GO TO 1404                        !Do we have a 2-axis along C
          IF ( L(2,4).NE.12 ) GO TO 1405                        !Yes, is it a 21?
          D(3,1) = 0.5
          D(3,2) = 0.5
          I21 = I21+1
          GO TO 1405

1404      CONTINUE
          IR = IR+1
          IF ( L(1,4).EQ.9 ) IM=IM+1
          IF ( L(1,4).EQ.4 ) IA=IA+1
          IF ( L(1,4).EQ.3 ) IB=IB+1
!         IF ( L(1,4).EQ.2 ) GO TO 500
          IF ( L(1,4).EQ.11 ) ID=ID+1
          IF ( L(1,2).EQ.2 .OR. L(1,2).EQ.10 ) D(3,3)=0.5
          IF ( L(1,3).EQ.2 .OR. L(1,3).EQ.10 ) D(3,3)=D(3,3)+0.5
1405      CONTINUE
                                                      !If there are 3 mirrors check for centering
                                                      !     which may alter the origin location
          IF ( IR.EQ.3 ) THEN                                    !  3 mirrors present.  Is the lattice centered?
            IF ( LCENT.EQ.1 ) THEN                              !No
                                                      !Yes.  Is it A-centered?
            ELSE IF ( LCENT.EQ.2 ) THEN                        !An A-centered lattice.
              IF ( IB+IC.EQ.1 .AND. IA.NE.2 ) THEN            !If only one B or C glide present relocate the mirrors by A
                D(2,2) = D(2,2)+0.5
                D(3,3) = D(3,3)+0.5
              END IF
            ELSE IF ( LCENT.EQ.3 ) THEN                        !A B-centered lattice
              IF ( IA+IC.EQ.1 .AND. IB.NE.2 ) THEN
                D(1,1) = D(1,1)+0.5
                D(3,3) = D(3,3)+0.5
              END IF
            ELSE IF ( LCENT.EQ.4 ) THEN                        !A C-centered lattice
              IF ( IA+IB.EQ.1 .AND. IC.NE.2 ) THEN
                D(1,1) = D(1,1)+0.5
                D(2,2) = D(2,2)+0.5
            END IF
            ELSE IF ( LCENT.EQ.5 ) THEN                        !It is I-centered
              IF ( IA+IB+IC.EQ.1 ) THEN                        !Yes.  if only 1 glide plane shift the mirrors by I
                D(1,1) = D(1,1)+0.5
                D(2,2) = D(2,2)+0.5
                D(3,3) = D(3,3)+0.5
              END IF
            END IF
          ELSE                                          !Less than 3 mirrors. set up the 2-axes locations
            IF ( I21.EQ.4 .OR. I21.EQ.5 .OR. I21.EQ.7 ) D(1,2)=0.0
            IF ( I21.EQ.6 .OR. I21.EQ.7 ) D(1,3)=0.0
            IF ( I21.EQ.3 ) D(2,1)=0.0
            IF ( I21.EQ.2 .OR. I21.EQ.6 .OR. I21.EQ.7 ) D(2,3)=0.0
            IF ( I21.EQ.1 .OR. I21.EQ.3 .OR. I21.EQ.7 ) D(3,1)=0.0
            IF ( I21.EQ.5 ) D(3,2)=0.0
            IF ( IM.LE.0 ) THEN
            ELSE IF ( IM.EQ.1 .AND. (I21.EQ.4 .OR. I21.EQ.2)
     1        .AND. D(3,3).NE.0.0 ) THEN
              D(3,3) = 0.0
              D(3,2) = D(3,2)+0.5
            ELSE IF ( IM.EQ.2 .AND. (I21.EQ.4 .OR. I21.EQ.1)
     1        .AND. D(2,2).NE.0.0 ) THEN
              D(2,2) = 0.0
              D(2,1) = D(2,1)+0.5
            ELSE IF ( IM.EQ.4 .AND. (I21.EQ.2 .OR. I21.EQ.1)
     1        .AND. D(1,1).NE.0.0 ) THEN
              D(1,1) = 0.0
              D(1,3) = D(1,3)+0.5
            END IF
          END IF
          GO TO 210
        END IF
      END IF
209   CONTINUE
      I209 = 1
210   CONTINUE
      RETURN
500   CONTINUE
      IF ( IER.EQ.0 ) IER=5
      RETURN
      END
