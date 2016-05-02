      SUBROUTINE PYLOADSCF(NATP,ATYPES,SFDAT,DEBG)
        
Cf2py intent(in) NATP
Cf2py intent(in) ATYPES
Cf2py intent(in) SFDAT
cf2py depend(NATP) ATYPES,SFDAT
cf2py intent(in) DEBG
            
      INCLUDE 'DIFFaXsubs/DIFFaX.par'
      INCLUDE 'DIFFaXsubs/DIFFaX.inc'
                
      INTEGER*4 NATP,I,J
      CHARACTER*4 ATYPES(NATP)
      REAL*4  SFDAT(9,NATP)
      LOGICAL DEBG
                
C fill common x-ray scattering factors
      debug = DEBG
      DO J=1,NATP
        WRITE(atom_l(J),'(A4)') ATYPES(J)
        DO I=1,9
          x_sf(I,J) = SFDAT(I,J)
        END DO
        if (debug) print '(1x,a4,9f10.6)',ATYPES(J),(x_sf(I,J),I=1,9)
      END DO
      intp_F = .TRUE.
      n_atoms = NATP
      RETURN
      END
        
      SUBROUTINE PYGETINST(LAMB,TTHMIN,TTHMAX,DELTTH)
        
Cf2py intent(in) LAMB
Cf2py intent(in) TTHMIN
Cf2py intent(in) TTHMAX
Cf2py intent(in) DELTTH
    
      INCLUDE 'DIFFaXsubs/DIFFaX.par'
      INCLUDE 'DIFFaXsubs/DIFFaX.inc'

      REAL*8 LAMB,TTHMIN,TTHMAX,DELTTH
        
      lambda = lamb
      th2_min = TTHMIN*DEG2RAD
      th2_max = TTHMAX*DEG2RAD
      d_theta = DELTTH*HALF*DEG2RAD
        
      RETURN
      END
        
      SUBROUTINE PYGETCLAY(CNTRLS,LAUESYM,WDTH,NST,STSEQ)
        
Cf2py intent(in) CNTRLS
Cf2py intent(in) LAUESYM
Cf2py intent(in) WDTH
Cf2py intent(in) NST
Cf2py intent(in) STSEQ
cf2py depend(NST) STSEQ
      
      INCLUDE 'DIFFaXsubs/DIFFaX.par'
      INCLUDE 'DIFFaXsubs/DIFFaX.inc'

      CHARACTER*12 LAUESYM
      INTEGER*4 CNTRLS(7),NST,STSEQ(NST),I
      REAL*8 WDTH(2)                  
      LOGICAL*4 ok,GETLAY
      EXTERNAL GETLAY
                                      
      PI = FOUR * atan(ONE)
      PI2 = TWO * PI
      DEG2RAD = PI / ONE_EIGHTY
      RAD2DEG = ONE_EIGHTY / PI
      rad_type = X_RAY
      lambda = 0.1
      trim_origin = .TRUE.
      blurring = NONE
      loglin = 1
      tolerance = 0.01
      finite_width = .TRUE.
      Wa = WDTH(1)*10000.
      Wb = WDTH(2)*10000.
      IF (Wa.GE.10000.) finite_width = .FALSE.
      WRITE(pnt_grp,'(A12)') LAUESYM
      SymGrpNo = CNTRLS(1)
      check_sym = .TRUE.
      full_shrp = 1
      full_brd = 1
C CNTRLS = [laueId,planeId,lmax,mult,StkType,StkParm,ranSeed]
      bitdepth = 16
      ok = .TRUE.
      scaleint = FLOAT(CNTRLS(4))
C fill in stacking seq stuff                  
      IF (CNTRLS(5).NE.0) THEN
        xplcit = .TRUE.
        recrsv = .FALSE.
        IF (CNTRLS(6).EQ.1) THEN
            rndm = .TRUE.
        ELSE
            rndm = .FALSE.
            l_cnt = NST
            DO I=1,NST
              l_seq(I) = STSEQ(I)
            END DO        
        END IF    
      ELSE 
        recrsv = .TRUE.
        xplcit = .FALSE.
        IF (CNTRLS(6).NE.0) THEN
            l_cnt = CNTRLS(6)
            inf_thick = .FALSE.
        ELSE
            inf_thick = .TRUE.
        END IF
      END IF
      IF (rndm) ok = GETLAY()
      RETURN
      END
            
      SUBROUTINE PYCELLAYER(CELL,NATM,ATMTP,ATMXOU,NU,LSYM,NL,LNUM)
                    
Cf2py intent(in) CELL
Cf2py intent(in) NATM
Cf2py intent(in) ATMTP
Cf2py intent(in) ATMXOU
cf2py depend(NATM) ATMTP,ATMXOU
Cf2py intent(in) NU
Cf2py intent(in) LSYM
Cf2py depend(NU) LSYM
Cf2py intent(in) NL
Cf2py intent(in) LNUM 
Cf2py depend(NL) LNUM
                       
      INCLUDE 'DIFFaXsubs/DIFFaX.par'
      INCLUDE 'DIFFaXsubs/DIFFaX.inc'

      INTEGER*4 NATM,NL,LNUM(NL),NU,LSYM(NU)
      CHARACTER*4 ATMTP(NATM)
      REAL*8  CELL(4),ATMXOU(8,NATM)
      INTEGER*4 I,J,K,IL,IA

C fill Common - cell stuff & finish symmetry stuff
      cell_a = CELL(1)
      cell_b = CELL(2)
      cell_c = CELL(3)
      cell_gamma = CELL(4)*DEG2RAD
C fill common layer stuff - atoms & symm
      DO I=1,NATM
        IL = NINT(ATMXOU(1,I))
        IA = NINT(ATMXOU(2,I))        
        a_type(IA,IL) = NINT(ATMXOU(3,I))
        a_number(IA,IL) = IA
        WRITE(a_name(IA,IL),'(A4)') ATMTP(I)
        DO K=1,3
            a_pos(K,IA,IL) = ATMXOU(K+3,I)
        END DO
        high_atom(IL) = max(high_atom(IL),a_pos(3,IA,IL))
        low_atom(IL) = min(low_atom(IL),a_pos(3,IA,IL))
        IF (LSYM(IL).EQ.CENTRO) THEN
            high_atom(IL) = MAX(high_atom(IL),-a_pos(3,IA,IL))
            low_atom(IL) = MIN(low_atom(IL),-a_pos(3,IA,IL))
        END IF
        a_occup(IA,IL) = ATMXOU(7,I)
        a_B(IA,IL) = ATMXOU(8,I)
        l_n_atoms(IL) = IA
        l_symmetry(IL) = LSYM(IL)
      END DO
      n_actual = IL
      n_layers = NL
      DO I=1,NL
        l_actual(I) = LNUM(I)
        DO J=1,NL
            Bs_zero(J,I) = .TRUE.
        END DO
      END DO
      all_Bs_zero = .TRUE.
      RETURN
      END

      SUBROUTINE PYGETTRANS(NL,TRP,TRX)
      
Cf2py intent(in) NL
Cf2py intent(in) TRP
Cf2py intent(in) TRX
Cf2py depend(NL) TRP,TRX
      
      
      INCLUDE 'DIFFaXsubs/DIFFaX.par'
      INCLUDE 'DIFFaXsubs/DIFFaX.inc'
        
      INTEGER*4 I,J,K
      INTEGER*4 NL
      REAL*4  TRP(NL,NL),TRX(NL,NL,3)
                               
C fill common transitions stuff
      DO J=1,NL
        DO I=1,NL
          l_alpha(J,I) = TRP(I,J)
          DO K=1,3
            l_r(K,J,I) = TRX(J,I,K)
          END DO
        END DO
      END DO
      RETURN
      END
        
      SUBROUTINE PYGETSPC(CNTRLS,NSADP,SADP)
        
Cf2py intent(in) CNTRLS
Cf2py intent(in) NSADP
Cf2py intent(in/out) SADP
Cf2py depend(NSADP) SADP
            
      INCLUDE 'DIFFaXsubs/DIFFaX.par'
      INCLUDE 'DIFFaXsubs/DIFFaX.inc'

      INTEGER*4 CNTRLS(7),NSADP,I,j,k
      REAL*8 SADP(NSADP),AGLQ16
      LOGICAL GETSPC,ok
        
      EXTERNAL AGLQ16,GETSPC

      DoSymDump = .FALSE.
    
      ok = .TRUE.
      CALL SPHCST()
      CALL DETUN()
      CALL OPTIMZ('GSAS-II',ok)
      If (debug) then
        print *,cell_a,cell_b,cell_c,cell_gamma,pnt_grp,SymGrpNo
        DoSymDump = .TRUE.
        print *,n_actual,(l_n_atoms(i),i=1,n_actual)
        do j=1,n_actual
          do i=1,l_n_atoms(j)
            print *,a_name(i,j),(a_pos(k,i,j),k=1,3)
          end do
        end do
        do i=1,n_layers
        print *,' layer',i
           do j=1,n_layers
              print *,'layer',j,l_alpha(i,j),(l_r(k,i,j),k=1,3)
           end do
        end do
        print *, recrsv,inf_thick,xplcit,rndm,l_cnt,has_l_mirror
        
C      print *,lambda,max_angle,h_bnd,k_bnd,l_bnd,no_trials,
C     1  rad_type,X_RAY,n_atoms
C      print *,(l_g(j),j=1,n_layers)
C      do j=1,n_layers
C        print *,(hx_ky(i,j),i=1,l_n_atoms(j))
C        print *,(mat(i,j),i=1,n_layers)
C        print *,(mat1(i,j),i=1,n_layers)
C        print *,(l_phi(i,j),i=1,n_layers)
C      end do
      end if
        
      ok = GETSPC(AGLQ16,'GSAS-II')              
      DO I=1,NSADP
        SADP(I) = spec(I)
      END DO
      RETURN
      END
        
      SUBROUTINE PYPROFILE(U,V,W,HW,BLUR,NBRD,BRDSPC)
        
Cf2py intent(in) U
Cf2py intent(in) V
Cf2py intent(in) W
Cf2py intent(in) HW
Cf2py intent(in) NBRD
Cf2py intent(in/out) BRDSPC
Cf2py depend(NBRD) BRDSPC
                
      INCLUDE 'DIFFaXsubs/DIFFaX.par'
      INCLUDE 'DIFFaXsubs/DIFFaX.inc'
        
      INTEGER*4 BLUR,i,NBRD    
      REAL*8 U,V,W,HW,BRDSPC(NBRD),tth_min
        
      tth_min = ZERO
            
      if (blur.eq.GAUSS) then
        FWHM = HW
        call GAUSSN(tth_min)
      else if (blur.eq.PV_GSS) then
        pv_u = U
        pv_v = V
        pv_w = W
        pv_gamma = ZERO
        call PV(tth_min)
      end if
      do i=1,NBRD
        BRDSPC(i) = brd_spc(i)
      end do
        
      RETURN
      END
        
      SUBROUTINE PYGETSADP(CNTRLS,NSADP,SADP,HKLIM,INCR,NBLK)
        
Cf2py intent(in) CNTRLS
Cf2py intent(in) NSADP
Cf2py intent(in/out) SADP
Cf2py depend(NSADP) SADP
Cf2py intent(out) HKLIM
Cf2py intent(out) INCR
Cf2py intent(out) NBLK
    
      INCLUDE 'DIFFaXsubs/DIFFaX.par'
      INCLUDE 'DIFFaXsubs/DIFFaX.inc'

      INTEGER*4 CNTRLS(7),NSADP,i_plane,hk_lim,i,j,k
      INTEGER*4 HKLIM,NBLK
      REAL*8 SADP(NSADP),AGLQ16,l_upper,INCR
      LOGICAL ok
        
      EXTERNAL AGLQ16                 
                    
      i_plane = CNTRLS(2)
      l_upper = CNTRLS(3)
      DoSymDump = .FALSE.
      if (debug) then
          print *,cell_a,cell_b,cell_c,cell_gamma
          print *,pnt_grp,SymGrpNo
          DoSymDump = .TRUE.
          print *,n_actual,(l_n_atoms(i),i=1,n_actual)
          do j=1,n_actual
            do i=1,l_n_atoms(j)
              print *,a_name(i,j),(a_pos(k,i,j),k=1,3)
            end do
          end do
          do i=1,n_layers
          print *,' layer',i
             do j=1,n_layers
                print *,'layer',j,l_alpha(i,j),(l_r(k,i,j),k=1,3)
             end do
          end do
          print *, recrsv,inf_thick,xplcit,rndm,l_cnt,has_l_mirror
      end if
      ok = .TRUE.
        
      CALL SPHCST()
      CALL DETUN()
      CALL OPTIMZ('GSAS-II',ok)
C      print *,lambda,max_angle,h_bnd,k_bnd,l_bnd,no_trials,
C     1  rad_type,X_RAY,n_atoms
C      print *,(l_g(j),j=1,n_layers)
C      do j=1,n_layers
C        print *,(hx_ky(i,j),i=1,l_n_atoms(j))
C        print *,(mat(i,j),i=1,n_layers)
C        print *,(mat1(i,j),i=1,n_layers)
C        print *,(l_phi(i,j),i=1,n_layers)
C      end do
      CALL GETSAD(AGLQ16,i_plane,l_upper,hk_lim,'GSAS-II',ok)
      NBLK = sadblock
      HKLIM = hk_lim+1
      INCR = dble(SADSIZE/2)/l_upper
      if (i_plane.eq.1) then
        INCR = INCR*sqrt(a0/c0)
      else if (i_plane.eq.2) then
        INCR = INCR*sqrt(b0/c0)
      else if (i_plane.eq.3) then
        INCR = INCR*sqrt((a0+b0+d0)/c0)
      else if (i_plane.eq.4) then
        INCR = INCR*sqrt((a0+b0-d0)/c0)
      end if
      do I=1,NSADP
        SADP(i) = spec(i)
      end do
      RETURN
      END

            
