      SUBROUTINE HEAD0
      COMMON /UNIT2/ IUNITB
C
C
C     SUBROUTINE 'HEAD0' ...
C        WRITE SHORT DESCRIPTION OF PROGRAM OUTPUT
C
C
C     --- WRITE HEADING FOR OUTPUT
      WRITE(IUNITB,6100)
      WRITE(IUNITB,6110)
      WRITE(IUNITB,6120)
      WRITE(IUNITB,6122)
      WRITE(IUNITB,6130)
      WRITE(IUNITB,6140)
      WRITE(IUNITB,6150)
      WRITE(IUNITB,6160)
      WRITE(IUNITB,9020)
      WRITE(IUNITB,6170)
      WRITE(IUNITB,6180)
      WRITE(IUNITB,6190)
      WRITE(IUNITB,9040)
      RETURN
 6100 FORMAT(1H1,////T45,'*** NIST*LATTICE ***'/)
 6110 FORMAT(1X,T34,'A PROGRAM TO ANALYZE LATTICE RELATIONSHIPS')
 6120 FORMAT(1X,T44,'Version of Spring, 1991'/)
 6122 FORMAT(1X,/T54,'by'//)
 6130 FORMAT(1X,T36,'Vicky Lynn Karen  and  Alan D. Mighell'/)
 6140 FORMAT(1X,T32,'National Institute of Standards and Technology')
 6150 FORMAT(1X,T33,'Materials Science and Engineering Laboratory')
 6160 FORMAT(1X,T41,'Gaithersburg, Maryland 20899')
 6170 FORMAT(1X,T22,'*** Patent Pending on certain algorithms used in th
     $is program ***')
 6180 FORMAT(1X,T17,'Copyright by the U.S. Secretary of Commerce on beha
     $lf of the United States')
 6190 FORMAT(1X,T45,'All Rights Reserved.')
 9020 FORMAT(//)
 9040 FORMAT(////)
      END
