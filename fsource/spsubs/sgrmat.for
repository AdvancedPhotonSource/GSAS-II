      SUBROUTINE SGRMAT(IOP,RT,N,A,B,C,D,E,F,G,H,O)

!Purpose:      S.R. to create a 3,3 matrix from 9 scalers

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 

C       THIS PROGRAM WAS DEVELOPED FOR
C                    THE DIVISION OF CHEMISTRY
C                               OF
C               THE NATIONAL RESEARCH COUNCIL OF CANADA
C                               BY
C       ALLEN C. LARSON, 14 CERRADO LOOP, SANTA FE, NM 87505, USA

!Calling sequence parameters:

      INTEGER*4     IOP                 !Matrix generator count
      REAL*4        RT(5,4,25)          !Matrix to be generated
      INTEGER*4     N                   !Number of the matrix to be generated
      REAL*4        A,B,C,D,E,F,G,H,O   !Matrix terms

!Local varaibles:

!Code:

      RT(1,1,N) = A
      RT(1,2,N) = B
      RT(1,3,N) = C
      RT(1,4,N) = 0.0
      RT(2,1,N) = D
      RT(2,2,N) = E
      RT(2,3,N) = F
      RT(2,4,N) = 0.0
      RT(3,1,N) = G
      RT(3,2,N) = H
      RT(3,3,N) = O
      RT(3,4,N) = 0.0
      RT(4,1,N) = 0.0
      RT(4,2,N) = 0.0
      RT(4,3,N) = 0.0
      RT(4,4,N) = 1.0
      RT(5,1,N) = 81*(2*RT(1,1,N)+3*RT(1,2,N)+4*RT(1,3,N))
     1  +9*(2*RT(2,1,N)+3*RT(2,2,N)+4*RT(2,3,N))
     1  +2*RT(3,1,N)+3*RT(3,2,N)+4*RT(3,3,N)
      RT(5,2,N) = 0.0                                          !Clear the translation info
      RT(5,3,N) = IOP
      RT(5,4,N) = 20.

      RETURN
      END
